#!/usr/bin/env python

import sys, os, re, argparse, shutil
import traceback
import os
import hashlib
import os
import logging
import multiprocessing
import subprocess
import urllib2
from glob import glob
from xml.etree import ElementTree

default_logger = logging.getLogger()
default_logger.setLevel(logging.DEBUG)
default_logger.addHandler(logging.StreamHandler())

FORCE_RUN=True

class HeaderException(Exception):
    pass

class BAMException(Exception):
    pass

def get_bam_file(dirpath, logger=default_logger):
    """Get the BAM file in the directory - assumes one BAM file"""

    if not os.path.isdir(dirpath):
        msg = "Directory for BAM file does not exist: %s" % dirpath
        logger.error(msg)
        raise BAMException(msg)

    for filename in os.listdir(dirpath):
        if filename.endswith(".bam"):
            return os.path.join(dirpath, filename)


def lane_level_bam_from_fastq(RG, workdir, logger=default_logger):
    """Create lane level bam, assuming bamtofastq has already been run"""
    default_logger.info("Starting fastq2bam")
    first_reads = os.path.join(workdir, "%s_1.fq.gz" % RG["ID"])
    second_reads = os.path.join(workdir, "%s_2.fq.gz" % RG["ID"])
    unmatched_first_reads = os.path.join(workdir, "%s_o1.fq.gz" % RG["ID"])
    unmatched_second_reads = os.path.join(workdir, "%s_o2.fq.gz" % RG["ID"])
    single_reads = os.path.join(workdir, "%s_s.fq" % RG["ID"])

    bam_filename = "%s/%s.paired.bam" % (workdir, RG["ID"])

    if not os.path.exists(first_reads):
        logger.warning("skipping RG that has no paired reads: %s" % first_reads)
        return False, bam_filename

    if not os.path.exists(second_reads):
        logger.warning("skipping RG that has no paired reads %s" % second_reads)
        return False, bam_filename

    #don't md5sum at this point because we go ahead and reheader
    cmd = "fastqtobam I=%s I=%s gz=1 level=1 threads=3 RGID=%s:%s RGCN=%s RGPL=%s RGLB=%s RGPI=%s RGSM=%s RGPU=%s RGDT=%s > %s" % (first_reads, second_reads,
                   RG["CN"], RG["ID"], RG["CN"], RG["PL"], RG["LB"], RG["PI"], RG["SM"], RG["PU"], RG["DT"],
                   bam_filename)
    default_logger.info("Running Command: %s" % (cmd))
    exit_code = os.system(cmd)

    if exit_code != 0:
        #remove the bam file if something goes wrong
        if os.path.exists(bam_filename):
            os.remove(bam_filename)
        logger.warning("removing bam file %s after fastqtobam returned error: %d" % (bam_filename, exit_code))
        return False, bam_filename

    return True, bam_filename

def process_rg(rg_dict, work_dir, output_dir, logger=default_logger):
    valid_bam, lane_level_bam = lane_level_bam_from_fastq(rg_dict, work_dir)
    if not valid_bam:
        return False, None
    return True, lane_level_bam


def gen_unaligned_data(bam_filename, metadata, work_dir, output_dir, header_output, num_processes=4, to_bam=False, logger=default_logger ):
    """
    The bulk of the work, calls splitting, generates new headers, generates initial
    unaligned BAM, reheaders with new headers
    """

    #get the read groups from the original sample level BAM
    exit_code = os.system("samtools view -H %s | grep \"@RG\" > %s" %(bam_filename, header_output))
    if exit_code != 0:
        print "Failure in bam splitting during read group extraction from %s" % bam_filename
        return 1

    if to_bam:
        #create the read group fastqs
        try:
            cmd = "bamtofastq outputperreadgroup=1 gz=1 level=1 inputbuffersize=2097152000 tryoq=1 outputdir=%s T=`mktemp -p %s bamtofastq_XXXXXXXXX` < %s" %(work_dir, work_dir, bam_filename)
            logger.info("Running %s" % cmd)
            subprocess.check_call(cmd, shell=True)
        except:
            print "Failure in bam splitting"
            return 1
    else:
        #create the read group fastqs
        try:
            cmd = "bamtofastq outputperreadgroup=1 level=1 inputbuffersize=2097152000 tryoq=1 outputdir=%s T=`mktemp -p %s bamtofastq_XXXXXXXXX` < %s" %(output_dir, work_dir, bam_filename)
            logger.info("Running %s" % cmd)
            subprocess.check_call(cmd, shell=True)
            for a in glob(os.path.join(output_dir, "*.fq")):

                shutil.move( a, a[:-3] + ".fastq" )
            return 0
        except:
            print "Failure in bam splitting"
            return 1


    rg_file = open(header_output, "r")
    pool = multiprocessing.Pool(processes=num_processes)
    results = []
    for line in rg_file:
        rg_dict = get_read_group_info(line)
        r = pool.apply_async(process_rg, (rg_dict, work_dir, output_dir))
        results.append(r)
    rg_file.close()

    out = []
    for r in results:
        out.append(r.get())

    clean_up_dir(output_dir)
    if not all( a[0] for a in out ):
        #one of the read group bamtofastq failed
        return 1
    with open(os.path.join(output_dir, "results.list"), "w") as out_handle:
        for ok, file_name in out:
            out_handle.write("%s\n" % (file_name))

    return 0


def clean_up_dir(dirname):
    """Remove all the fastq files and intermediate bams/headers"""

    for filename in os.listdir(dirname):
        if(filename.endswith(".fq")
            or filename.endswith(".paired.bam")
            or filename.endswith(".sam")):
            filepath = os.path.join(dirname, filename)
            os.remove(filepath)
            assert(not(os.path.exists(filepath)))

def make_new_dir(path):
    """Create a new directory if it doesn't exist"""

    if(not(os.path.isdir(path))):
        os.makedirs(path)

    return path

def calc_md5sum(bam_filename):
    """What it says"""

    md5worker = hashlib.md5()
    bam_file = open(bam_filename)
    block_size = 128
    while True:
        data = bam_file.read(block_size)
        if not data:
            break
        md5worker.update(data)

    return md5worker.hexdigest()

if __name__ == '__main__':
    basedir = os.path.abspath(os.path.dirname( __file__ ))

    parser = argparse.ArgumentParser(prog='biobambam_split.py', description='Create unaligned FASTQ/BAM files from original BAM file')
    parser.add_argument('--bam_path', type=str, help='path/to/tcga/data/labeled_by_analysis_id', required=True)
    parser.add_argument('--bam', action="store_true", default=False)
    #parser.add_argument('--normal_id', type=str, help='UUID for normal analysis (analysis_id)', default=None)
    #parser.add_argument('--tumor_id', type=str, help='Comma separated list of tumor analysis UUIDs (analysid_id(s))', default=None)
    parser.add_argument('--work_dir', type=str, help='path/to/output/directory', default=None)
    parser.add_argument('--output_dir', type=str, help='path/to/output/directory', required=True)
    parser.add_argument('--header', type=str, required=True)

    args = parser.parse_args()

    if args.work_dir is None:
        args.work_dir = args.output_dir
    exit_code = 0
    output_dir = make_new_dir(args.output_dir)
    work_dir = make_new_dir(args.work_dir)
    try:
        exit_code = gen_unaligned_data(bam_filename=args.bam_path, metadata={}, work_dir=work_dir, output_dir=output_dir, header_output=args.header, to_bam=args.bam)

    	"""
        if args.tumor_id is None and args.normal_id is not None:
            metadata = parse_cghub_metadata(args.normal_id)
            metadata['use_cntl'] = 'N/A'
            exit_code = bam_utils.gen_unaligned_bam(args.bam_path, args.normal_id, metadata, specimen_dict, work_dir, output_dir)
        elif args.tumor_id is not None and args.normal_id is not None:
            metadata = header_utils.parse_cghub_metadata(args.tumor_id)
            metadata['use_cntl'] = args.normal_id
            exit_code = bam_utils.gen_unaligned_bam(args.bam_path, args.tumor_id, metadata, specimen_dict, work_dir, output_dir)
        else:
            print "Please define --normal_id or (--normal_id and --tumor_id)"
            sys.exit(1)
        """
    except:
        print "PCAP SPLIT Failure!!!"
        traceback.print_exc(file=sys.stderr)
        exit_code = 1
    if output_dir != work_dir:
        shutil.rmtree(work_dir)

    sys.exit(exit_code)
