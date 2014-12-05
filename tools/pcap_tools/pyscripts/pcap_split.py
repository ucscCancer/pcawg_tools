#!/usr/bin/env python

import sys, os, re, argparse, shutil
import dateutil.parser
import header_utils, bam_utils, utils
import traceback

basedir = os.path.abspath(os.path.dirname( __file__ ))

def parse_specimen_dict(spec_filename, header_exists=True):
    """Parse the mapping between TCGA sample codes and ICGC specimen vocab"""

    spec_file = open(spec_filename, "r")
    specimen_dict = dict()

    if header_exists:
        header = spec_file.readline()

    for line in spec_file:
        line = line.rstrip()
        line = line.split("\t")

        if(not(line[0] == "")):
            specimen_dict[line[0]] = (line[1], line[2])

    spec_file.close()
    return specimen_dict

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='pcap_split.py', description='Create unaligned BAM files')
    parser.add_argument('--bam_path', type=str, help='path/to/tcga/data/labeled_by_analysis_id', required=True)
    parser.add_argument('--normal_id', type=str, help='UUID for normal analysis (analysis_id)', default=None)
    parser.add_argument('--tumor_id', type=str, help='Comma separated list of tumor analysis UUIDs (analysid_id(s))', default=None)
    parser.add_argument('--work_dir', type=str, help='path/to/output/directory', default=None)    
    parser.add_argument('--output_dir', type=str, help='path/to/output/directory', required=True)
    parser.add_argument('--specimen_map', type=str, default=os.path.join(basedir,'tcga_dcc_specimen_type.txt'), help='path/to/tcga/icgc/sample_code_specimen_mapping')

    args = parser.parse_args()
    
    specimen_dict = parse_specimen_dict(args.specimen_map)
    
    if args.work_dir is None:
        args.work_dir = args.output_dir
    exit_code = 0
    output_dir = utils.make_new_dir(args.output_dir)
    work_dir = utils.make_new_dir(args.work_dir)
    try:
        if args.tumor_id is None and args.normal_id is not None:
            metadata = header_utils.parse_cghub_metadata(args.normal_id)
            metadata['use_cntl'] = 'N/A'
            exit_code = bam_utils.gen_unaligned_bam(args.bam_path, args.normal_id, metadata, specimen_dict, work_dir, output_dir)
        elif args.tumor_id is not None and args.normal_id is not None:
            metadata = header_utils.parse_cghub_metadata(args.tumor_id)
            metadata['use_cntl'] = args.normal_id
            exit_code = bam_utils.gen_unaligned_bam(args.bam_path, args.tumor_id, metadata, specimen_dict, work_dir, output_dir)
        else:
            print "Please define --normal_id or (--normal_id and --tumor_id)"
            sys.exit(1)
    except:
        print "PCAP SPLIT Failure!!!"
        traceback.print_exc(file=sys.stderr)
        exit_code = 1
    if output_dir != work_dir:
        shutil.rmtree(work_dir)

    sys.exit(exit_code)

