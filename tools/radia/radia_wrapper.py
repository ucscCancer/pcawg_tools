#!/usr/bin/env python

import sys
import re
import os
import shutil
import subprocess
import tempfile
import vcf
import argparse
import logging
from string import Template
from multiprocessing import Pool


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


def radia_cmd():
    RADIA_CMD = """python /opt/radia-1.1.0/scripts/radia.pyc
    ${SAMPLE_ID} ${CHR}
    -n ${DNA_NORMAL_BAM}
    -t ${DNA_TUMOR_BAM}
    -f ${REF_SEQ}
    -i ${REF_NAME}
    -o ${OUT_PATH}
    """
    template = Template(RADIA_CMD.replace("\n", " "))

    cmd = template.substitute(
        dict(
            REF_SEQ=ref_seq,
            DNA_TUMOR_BAM=tumor_bam,
            DNA_NORMAL_BAM=normal_bam
        )
    yield cmd, "%s.%s" % (output_base, i)



def run_radia(args):

    workdir = tempfile.mkdtemp(dir=args['workdir'], prefix="radia_work_")

    tumor_bam = os.path.join(workdir, "tumor.bam")
    normal_bam = os.path.join(workdir, "normal.bam")
    os.symlink(os.path.abspath(args["input_file:normal"]), normal_bam)
    os.symlink(os.path.abspath(args['input_file:tumor']),  tumor_bam)

    if args['input_file:index:normal'] is not None:
        os.symlink(os.path.abspath(args["input_file:index:normal"]), normal_bam + ".bai")
    elif os.path.exists(os.path.abspath(args["input_file:normal"]) + ".bai"):
        os.symlink(os.path.abspath(args["input_file:normal"]) + ".bai", normal_bam + ".bai")
    else:
        subprocess.check_call( ["/usr/bin/samtools", "index", normal_bam] )

    if args['input_file:index:tumor'] is not None:
        os.symlink(os.path.abspath(args["input_file:index:tumor"]), tumor_bam + ".bai")
    elif os.path.exists(os.path.abspath(args["input_file:tumor"]) + ".bai"):
        os.symlink(os.path.abspath(args["input_file:tumor"]) + ".bai", tumor_bam + ".bai")
    else:
        subprocess.check_call( ["/usr/bin/samtools", "index", tumor_bam] )

    ref_seq = os.path.join(workdir, "ref_genome.fasta")
    ref_dict = os.path.join(workdir, "ref_genome.dict")
    os.symlink(os.path.abspath(args['reference_sequence']), ref_seq)
    subprocess.check_call( ["/usr/bin/samtools", "faidx", ref_seq] )
    subprocess.check_call( [args['java'], "-jar",
        args['dict_jar'],
        "R=%s" % (ref_seq),
        "O=%s" % (ref_dict)
    ])

    cmds = list(call_cmd_iter(
        ref_seq=ref_seq,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        )
    )

    rvals = cmds_runner(list(a[0] for a in cmds), args['ncpus'])

    vcf_writer = None
    for cmd, file in cmds:
        vcf_reader = vcf.Reader(filename=file + ".vcf")
        if vcf_writer is None:
            vcf_writer = vcf.Writer(open(os.path.join(args['vcf']), "w"), vcf_reader)
        for record in vcf_reader:
            vcf_writer.write_record(record)
    vcf_writer.close()

    if args['out'] is not None:
        with open(args['out'], "w") as handle:
            for cmd, file in cmds:
                with open(file + ".out") as ihandle:
                    for line in ihandle:
                        handle.write(line)

    first_file = True
    if args['coverage_file'] is not None:
        with open(args['coverage_file'], "w") as handle:
            for cmd, file in cmds:
                with open(file + ".coverage") as ihandle:
                    first_line = True
                    for line in ihandle:
                        if first_line:
                            if first_file:
                                handle.write(line)
                                first_line = False
                                first_file = False
                        else:
                            handle.write(line)


    if not args['no_clean']:
        shutil.rmtree(workdir)
