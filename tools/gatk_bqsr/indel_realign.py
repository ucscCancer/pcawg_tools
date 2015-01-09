#!/usr/bin/env python

import sys
import re
import os
import shutil
import subprocess
import tempfile
import argparse
import logging
from string import Template
from multiprocessing import Pool

def which(cmd):
    cmd = ["which",cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: return None
    return res

def fai_chunk(path, blocksize):
    seq_map = {}
    with open( path ) as handle:
        for line in handle:
            tmp = line.split("\t")
            seq_map[tmp[0]] = long(tmp[1])

    for seq in seq_map:
        l = seq_map[seq]
        for i in xrange(1, l, blocksize):
            yield (seq, i, min(i+blocksize-1, l))

def cmd_caller(cmd):
    logging.info("RUNNING: %s" % (cmd))
    print "running", cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode

def cmds_runner(cmds, cpus):
    p = Pool(cpus)
    values = p.map(cmd_caller, cmds, 1)
    return values

def call_scan(java, gatk, ncpus, ref_seq, input_bam, output_base, known_vcfs):

    known_str = " ".join(list("-known %s" % (a) for a in known_vcfs))
    out = "%s.intervals" % (output_base)

    template = Template("""
${JAVA}
-Xmx8g -XX:ParallelGCThreads=2 -jar ${GATK}
-T RealignerTargetCreator
-nt ${NCPUS}
-R ${REF_SEQ}
-I ${INPUT_BAM}
${KNOWN_STR}
--downsampling_type NONE
-o ${OUT}
""".replace("\n", " "))
    cmd = template.substitute(
        dict(
            JAVA=java,
            REF_SEQ=ref_seq,
            GATK=gatk,
            INPUT_BAM=input_bam,
            OUTPUT_BASE=output_base,
            KNOWN_STR=known_str,
            NCPUS=ncpus,
            OUT=out
    ))
    return cmd, out


def call_realign_iter(java, gatk, ref_seq, block_size, input_bam, output_base, target_intervals, known_vcfs):
    known_str = " ".join(list("-known %s" % (a) for a in known_vcfs))

    template = Template("""
${JAVA}
-Xmx8g -XX:ParallelGCThreads=2 -jar ${GATK}
-T IndelRealigner
-R ${REF_SEQ}
-I ${INPUT_BAM}
-L ${INTERVAL}
-targetIntervals ${TARGET_INTERVALS}
--downsampling_type NONE
${KNOWN_STR}
-maxReads 720000 -maxInMemory 5400000 \
-o  ${OUTPUT_BASE}.${BLOCK_NUM}.bam
""".replace("\n", " "))

    for i, block in enumerate(fai_chunk( ref_seq + ".fai", block_size ) ):
        cmd = template.substitute(
            dict(
                JAVA=java,
                REF_SEQ=ref_seq,
                GATK=gatk,
                BLOCK_NUM=i,
                INTERVAL="%s:%s-%s" % (block[0], block[1], block[2]) ),
                INPUT_BAM=input_bam,
                OUTPUT_BASE=output_base,
                KNOWN_STR=known_str,
                TARGET_INTERVALS=target_intervals
        )
        yield cmd, "%s.%s.bam" % (output_base, i)



def run_indel_realign(args):

    workdir = tempfile.mkdtemp(dir=args['workdir'], prefix="indel_work_")

    input_bam = os.path.join(workdir, "input.bam")
    os.symlink(os.path.abspath(args["input_file"]), input_bam)

    samtools = which("samtools")
    if args['input_file:index'] is not None:
        os.symlink(os.path.abspath(args["input_file:index"]), input_bam + ".bai")
    elif os.path.exists(os.path.abspath(args["input_file"]) + ".bai"):
        os.symlink(os.path.abspath(args["input_file"]) + ".bai", input_bam + ".bai")
    else:
        logging.info("Indexing Input %s" % (input_bam))
        subprocess.check_call( [samtools, "index", input_bam] )

    ref_seq = os.path.join(workdir, "ref_genome.fasta")
    ref_dict = os.path.join(workdir, "ref_genome.dict")
    os.symlink(os.path.abspath(args['reference_sequence']), ref_seq)
    subprocess.check_call( [samtools, "faidx", ref_seq] )
    subprocess.check_call( [args['java'], "-jar",
        args['dict_jar'],
        "R=%s" % (ref_seq),
        "O=%s" % (ref_dict)
    ])

    cmd, intervals = call_scan(java=args['java'],
        gatk=args['gatk_jar'],
        ncpus=args['ncpus'],
        ref_seq=ref_seq,
        input_bam=input_bam, output_base=os.path.join(workdir, "output.file"),
        known_vcfs=args['known'])
    if cmd_caller(cmd) != 0:
        raise Exception("Program fail")

    cmds = list(call_realign_iter(ref_seq=ref_seq,
        java=args['java'],
        gatk=args['gatk_jar'],
        block_size=args['b'],
        input_bam=input_bam,
        output_base=os.path.join(workdir, "output.file"),
        known_vcfs=args['known'],
        target_intervals=intervals
        )
    )

    rvals = cmds_runner(list(a[0] for a in cmds), args['ncpus'])
    if any( rvals ):
        raise Exception("Program fail")

    cmd = [samtools, "merge", args['out'] ] + list( a[1] for a in cmds )
    logging.info("Running: %s" % " ".join(cmd))
    subprocess.check_call(cmd)

    if not args['no_clean']:
        shutil.rmtree(workdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-Ii", "--input_file:index")
    parser.add_argument("-I", "--input_file", required=True)
    parser.add_argument("-R", "--reference-sequence", required=True)
    parser.add_argument("--ncpus", type=int, default=8)
    parser.add_argument("--workdir", default="/tmp")
    parser.add_argument("--known", action="append", default=[])
    parser.add_argument("-o", "--out")
    parser.add_argument("--no-clean", action="store_true", default=False)
    parser.add_argument("--java", default="/usr/bin/java")

    parser.add_argument("-b", type=long, help="Parallel Block Size", default=250000000)

    parser.add_argument("--gatk-jar", default="/opt/GenomeAnalysisTK.jar")
    parser.add_argument("--dict-jar", default="/opt/picard/CreateSequenceDictionary.jar")

    args = parser.parse_args()
    run_indel_realign(vars(args))
