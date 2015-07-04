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



def get_bam_seq(inputBamFile):
	samtools = which("samtools")
	cmd = [samtools, "idxstats", inputBamFile]
	process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE)
	stdout, stderr = process.communicate()
	seqs = []
	for line in stdout.split("\n"):
		tmp = line.split("\t")
		if len(tmp) == 4 and tmp[2] != "0":
			seqs.append(tmp[0])
	return seqs


def cmd_caller(cmd):
    logging.info("RUNNING: %s" % (cmd))
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print "Failed job: %s" % (cmd)
        print "--stdout--"
        print stdout
        print "--stderr--"
        print stderr
    return p.returncode

def cmds_runner(cmds, cpus):
    p = Pool(cpus)
    values = p.map(cmd_caller, cmds, 1)
    return values

def call_scan(java, gatk, ncpus, ref_seq, input_list, known_vcfs, intervals, mem="8g"):

    known_str = " ".join(list("-known %s" % (a) for a in known_vcfs))
    #out = "%s.intervals" % (output_base)
    out = intervals
    template = Template("""
${JAVA}
-Xmx${MEM} -XX:ParallelGCThreads=2 -jar ${GATK}
-T RealignerTargetCreator
-disable_auto_index_creation_and_locking_when_reading_rods
-nt ${NCPUS}
-R ${REF_SEQ}
-I ${INPUT_LIST}
${KNOWN_STR}
-o ${OUT}
""".replace("\n", " "))
#output_base
    cmd = template.substitute(
        dict(
            JAVA=java,
            REF_SEQ=ref_seq,
            GATK=gatk,
            INPUT_LIST=input_list,
        #    OUTPUT_BASE=output_base,
            KNOWN_STR=known_str,
            NCPUS=ncpus,
            OUT=out,
            MEM=mem
    ))
    return cmd #, out


def call_realign_iter(java, gatk, ref_seq, blocks, input_list, work_dir, target_intervals, known_vcfs, mem="8g"):
    known_str = " ".join(list("-known %s" % (a) for a in known_vcfs))

    input_files = []
    with open(input_list) as handle:
        for line in handle:
            input_files.append(line.rstrip())

    template = Template("""
${JAVA}
-Xmx${MEM} -XX:ParallelGCThreads=2 -jar ${GATK}
-T IndelRealigner
-R ${REF_SEQ}
-I ${INPUT_LIST}
-L '${INTERVAL}'
-disable_auto_index_creation_and_locking_when_reading_rods
-targetIntervals ${TARGET_INTERVALS}
${KNOWN_STR}
-nWayOut ${OUTPUT_MAP}
""".replace("\n", " "))

    #for i, block in enumerate(fai_chunk( ref_seq + ".fai", block_size ) ):
    for i, block in enumerate(blocks):
        outdir = os.path.join(work_dir, "block.%s" % (i))
        os.mkdir( outdir )
        output_map = os.path.join(work_dir, "output.%s.map" % (i))
        outputs = []
        with open(output_map, "w") as handle:
            for j, o in enumerate(input_files):
                out = os.path.join(work_dir, "output.%s.%d.bam" % (i, j))
                handle.write("input.%s.bam\t%s\n" % (j, out))
                outputs.append(out)

        cmd = template.substitute(
            dict(
                JAVA=java,
                REF_SEQ=ref_seq,
                GATK=gatk,
                BLOCK_NUM=i,
                #INTERVAL="%s:%s-%s" % (block[0], block[1], block[2]),
                INTERVAL="%s" % (block),
                INPUT_LIST=input_list,
                KNOWN_STR=known_str,
                OUTPUT_MAP=output_map,
                TARGET_INTERVALS=target_intervals,
                MEM=mem
            )
        )
        yield cmd, outputs

def call_realign_single(java, gatk, ref_seq, blocks, input_list, output_map, target_intervals, known_vcfs, mem="8g"):
    known_str = " ".join(list("-known %s" % (a) for a in known_vcfs))

    template = Template("""
${JAVA}
-Xmx${MEM} -XX:ParallelGCThreads=2 -jar ${GATK}
-T IndelRealigner
-R ${REF_SEQ}
-disable_auto_index_creation_and_locking_when_reading_rods
-I ${INPUT_LIST}
-targetIntervals ${TARGET_INTERVALS}
${KNOWN_STR}
-nWayOut ${OUTPUT_MAP}
""".replace("\n", " "))

    #for i, block in enumerate(fai_chunk( ref_seq + ".fai", block_size ) ):
    #for i, block in enumerate(blocks):
    cmd = template.substitute(
        dict(
            JAVA=java,
            REF_SEQ=ref_seq,
            GATK=gatk,
            #BLOCK_NUM=i,
            #INTERVAL="%s:%s-%s" % (block[0], block[1], block[2]),
            #INTERVAL="%s" % (block),
            INPUT_LIST=input_list,
            #OUTPUT_BASE=output_base,
            OUTPUT_MAP=output_map,
            KNOWN_STR=known_str,
            TARGET_INTERVALS=target_intervals,
            MEM=mem
        )
    )
    return cmd



def run_indel_realign(args):

    samtools = which("samtools")
    workdir = tempfile.mkdtemp(dir=args['workdir'], prefix="indel_work_")

    input_files = []
    seqs = set()
    for i, p in enumerate(args["input_files"]):
        input_bam = os.path.join(workdir, "input.%s.bam" % (i))
        os.symlink(os.path.abspath(p), input_bam)

        if len(args['input_files_index']) > i:
            logging.info("Index for %s from %s" % (p, args["input_files_index"][i]))
            os.symlink(os.path.abspath(args["input_files_index"][i]), input_bam + ".bai")
        elif os.path.exists(os.path.abspath(p) + ".bai"):
            os.symlink(os.path.abspath(p) + ".bai", input_bam + ".bai")
        else:
            logging.info("Indexing Input %s" % (input_bam))
            subprocess.check_call( [samtools, "index", input_bam] )
        input_files.append(input_bam)
        for a in get_bam_seq(input_bam):
            seqs.add(a)

    logging.info("Creating input list")
    input_list = os.path.join(workdir, "input.list")
    with open(input_list, "w") as handle:
        for p in input_files:
            handle.write(p + "\n")

    ref_seq = os.path.join(workdir, "ref_genome.fasta")
    ref_dict = os.path.join(workdir, "ref_genome.dict")
    os.symlink(os.path.abspath(args['reference_sequence']), ref_seq)
    logging.info("Indexing reference genome")
    subprocess.check_call( [samtools, "faidx", ref_seq] )
    subprocess.check_call( [args['java'], "-jar",
        args['dict_jar'],
        "R=%s" % (ref_seq),
        "O=%s" % (ref_dict)
    ])

    known_vcfs = []
    for i, v in enumerate(args['known']):
        lpath = os.path.join(workdir, "known.%s.vcf" % (i))
        os.symlink(os.path.abspath(v), lpath)
        known_vcfs.append(lpath)

    intervals = os.path.join(workdir, "scan.intervals")
    cmd = call_scan(java=args['java'],
        gatk=args['gatk_jar'],
        ncpus=args['ncpus'],
        ref_seq=ref_seq,
        input_list=input_list,
        #output_base=os.path.join(workdir, "output.file"),
        known_vcfs=known_vcfs,
        intervals=intervals,
        mem="%sg" % (args['mem']))
    logging.info("Calling RealignerTargetCreator")
    if cmd_caller(cmd) != 0:
        raise Exception("RealignerTargetCreator failed")

    if args['parallel_realign']:
        cmds = list(call_realign_iter(ref_seq=ref_seq,
            java=args['java'],
            gatk=args['gatk_jar'],
            blocks=seqs,
            input_list=input_list,
            work_dir=workdir,
            known_vcfs=known_vcfs,
            target_intervals=intervals,
            mem="%sg" % (args['mem'])
            )
        )
        logging.info("Calling IndelRealigner")
        rvals = cmds_runner(list(a[0] for a in cmds), args['ncpus'])
        if any( rvals ):
            raise Exception("IndelRealigner failed")

        merge_cmds = []
        for i, o in enumerate(args['out']):
            merge_cpus = min(4, args['ncpus'])
            cmd = [samtools, "merge", "-@%s" % (merge_cpus), o ] + list( a[1][i] for a in cmds )
            merge_cmds.append( " ".join(cmd) )
        rvals = cmds_runner(merge_cmds, args['ncpus'])
        if any( rvals ):
            raise Exception("samtools merge failed")

    else:
        output_map = os.path.join(workdir, "output.map")
        with open(output_map, "w") as handle:
            for i, o in enumerate(args['out']):
                handle.write("input.%s.bam\t%s\n" % (i, o))
        cmd = call_realign_single(ref_seq=ref_seq,
            java=args['java'],
            gatk=args['gatk_jar'],
            blocks=seqs,
            input_list=input_list,
            output_map=output_map,
            #output_base=os.path.join(workdir, "output.file"),
            known_vcfs=known_vcfs,
            target_intervals=intervals,
            mem="%sg" % (args['mem'])
            )
        logging.info("Calling IndelRealigner")
        if cmd_caller(cmd) != 0:
            raise Exception("IndelRealigner failed")

    if not args['no_clean']:
        shutil.rmtree(workdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-Ii", "--input_file:index", dest="input_files_index", default=[], action="append")
    parser.add_argument("-I", "--input_file", dest="input_files", action="append", required=True)
    parser.add_argument("-R", "--reference-sequence", required=True)
    parser.add_argument("--ncpus", type=int, default=8)
    parser.add_argument("--mem", type=int, default=8)
    parser.add_argument("--workdir", default="/tmp")
    parser.add_argument("--known", action="append", default=[])
    parser.add_argument("-o", "--out", action="append", required=True)
    parser.add_argument("--no-clean", action="store_true", default=False)
    parser.add_argument("--java", default="/usr/bin/java")
    parser.add_argument("--parallel-realign", action="store_true", default=False)

    #parser.add_argument("-b", type=long, help="Parallel Block Size", default=250000000)

    parser.add_argument("--gatk-jar", default="/opt/GenomeAnalysisTK.jar")
    parser.add_argument("--dict-jar", default="/opt/picard/CreateSequenceDictionary.jar")

    logging.basicConfig(level=logging.INFO)
    args = parser.parse_args()
    run_indel_realign(vars(args))
