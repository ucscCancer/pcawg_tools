#!/usr/bin/env python

import os
import sys
import re
import argparse
from glob import glob
import subprocess
import tempfile
import shutil

BROAD_README = """
This archive contains intermediate results from the Broad PCAWG variant calling
pipeline. They will be further processed by the Broad to generate the final
pipeline outputs that PCAWG groups will use for downstream analysis.
They are being shared with all PCAWG members to foster current and future
collaboration. 

Some of the files contain results from algorithms that are still experimental.
Other files might look like they are good for downstream analysis when really
they will give bad results.  For instance, our intermediate germline calls
require further steps such as joint-calling and filtering in order to become
good quality.

Please do not use these files without first consulting with the Broad.
"""

raw_temp_paths = [
    "collapse_het_sites_intervals/.*.het_sites_2_forcecall.intervals.gz$",
    "genotype_gvcf/.*.single_sample.vcf.gz$",
    "mutect_sg_gather/.*.call_stats.txt.gz$",
    "BreakPointer_Normal_sg_gather/.*.breakpoints.txt.gz$",
    "BreakPointer_Normal_sg_gather/.*.matched.sam.gz$",
    "snowman/.*.contigs.bam$",
    "snowman/.*.cigarmap.txt.gz$",
    "snowman/.*.r2c.bam$",
    "snowman/.*.contigs.bam.bai$",
    "snowman/.*.discordant.txt.gz$",
    "snowman/.*.contigs_all.sam.gz$",
    "snowman/.*.bps.txt.gz$",
    "snowman/.*.alignments.txt.gz$",
    "dRangerPreProcess_Normal_sg_gather/.*.all.isz.gz$",
    "fragcounter_tumor/cov.rds.gz$",
    "haplotypecaller_sg_sg_gather/.*.gvcf.gz$",
    "haplotypecaller_sg_sg_gather/.*.gvcf.gz.tbi$",
    "re_capseg_coverage_tumor_merged/.*.tumor.uncorrected_target_order.coverage.gz$",
    "dRangerPreProcess_Tumor_sg_gather/.*.all.isz.gz$",
    "getdRangerSupportingReads/.*.dRanger.supporting_reads.txt.gz$",
    "BreakPointer_Tumor_sg_gather/.*.breakpoints.txt.gz$",
    "BreakPointer_Tumor_sg_gather/.*.matched.sam.gz$",
    "dRanger_Finalize/.*.dRanger_results.detail.all.mat.gz$",
    "dRanger_Finalize/.*.dRanger_results.detail.all.txt.gz$",
    "dRanger_Finalize/.*.dRanger_results.somatic.txt.gz$",
    "dRanger_Finalize/.*.dRanger_results.detail.somatic.txt.gz$",
    "fragcounter_normal/cov.rds.gz$",
    "tokens/.*.tok.gz$",
    "mutect_het_sites_sg_gather/.*.call_stats.txt.gz$",
    "mutect_het_sites_sg_gather/.*.coverage.wig.txt.gz$",
    "re_capseg_coverage_normal_merged/.*.normal.uncorrected_target_order.coverage.gz$",
    "dRangerRun/.*.dRanger_results.mat.gz$",
    "dRangerRun/.*.dRanger_results.forBP.txt.gz$",
    "dRanger2VCF/.*.dRanger_results.detail.txt.gz"
]

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("workdir")
parser.add_argument("outdir")
parser.add_argument("--rename", nargs=2, default=None)

args = parser.parse_args()
input = os.path.abspath(args.input)
workdir = tempfile.mkdtemp( prefix="broad_repack_", dir=os.path.abspath(args.workdir) )
outdir = os.path.abspath(args.outdir)
if not os.path.exists(workdir):
    os.mkdir(workdir)

subprocess.check_call("tar xvf %s" % (input), shell=True, cwd=workdir)

for root, dirs, files in os.walk(workdir):
    if os.path.basename(root) == 'links_for_broad':
        sample_name = os.path.basename(os.path.dirname(root))
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for fdir in dirs:
            d = os.path.join(root, fdir)
            for file in glob(os.path.join(d, "*")):
                dstdir = os.path.join(outdir, fdir)
                if not os.path.exists(dstdir):
                    os.mkdir(dstdir)
                mark_temp = False
                for t in raw_temp_paths:
                    if re.search(t, file):
                        mark_temp = True
                if mark_temp:
                    dst_path = os.path.join(dstdir, os.path.basename(file) + ".raw.tmp")
                else:
                    dst_path = os.path.join(dstdir, os.path.basename(file))
                
                if args.rename is not None:
                    dst_path = re.sub(r"/%s\." % (args.rename[0]), "/%s." % (args.rename[1]), dst_path)

                os.symlink(file, dst_path)

with open(os.path.join(outdir, "README"), "w") as handle:
    handle.write(BROAD_README)

subprocess.check_call("tar -C %s -cvhf - %s | tee %s.tar | md5sum - | awk '{print$1}' > %s.tar.md5" % (os.path.dirname(outdir),os.path.basename(outdir),outdir,outdir), shell=True, cwd=workdir)
shutil.rmtree( workdir )
shutil.rmtree( outdir )
