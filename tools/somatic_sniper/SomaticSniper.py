#!/usr/bin/env python

import argparse
import pysam
import os
import logging
import subprocess
import os.path
from datetime import date
from urlparse import urlparse, urlunparse

"""This script runs SomaticSniper
Required Option:
        -f FILE   REQUIRED reference sequence in the FASTA format

Options:
        -v        Display version information

        -q INT    filtering reads with mapping quality less than INT [0]
        -Q INT    filtering somatic snv output with somatic quality less than  INT [15]
        -L FLAG   do not report LOH variants as determined by genotypes
        -G FLAG   do not report Gain of Reference variants as determined by genotypes
        -p FLAG   disable priors in the somatic calculation. Increases sensitivity for solid tumors
        -J FLAG   Use prior probabilities accounting for the somatic mutation rate
        -s FLOAT  prior probability of a somatic mutation (implies -J) [0.010000]
        -T FLOAT  theta in maq consensus calling model (for -c/-g) [0.850000]
        -N INT    number of haplotypes in the sample (for -c/-g) [2]
        -r FLOAT  prior of a difference between two haplotypes (for -c/-g) [0.001000]
        -n STRING normal sample id (for VCF header) [NORMAL]
        -t STRING tumor sample id (for VCF header) [TUMOR]
        -F STRING select output format [classic]
           Available formats:
             classic
             vcf
             bed
"""

def sniper_argparser():
    parser = argparse.ArgumentParser(description='Run SomaticSniper')
    parser.add_argument('-f', metavar='REFERENCE', required=True, help='reference sequence in the FASTA format')
    parser.add_argument('-q', required=False, metavar='MAPQ', type=int, help='filtering reads with mapping quality less than [%(default)s]', default=0)
    parser.add_argument('-Q', required=False, metavar='SOMATICQ', type=int, help='filtering somatic snv output with somatic quality less than [%(default)s]', default=15)
    parser.add_argument('-L', required=False, action='store_true', help='do not report LOH variants as determined by genotypes')
    parser.add_argument('-G', required=False, action='store_true', help='do not report Gain of Reference variants as determined by genotypes')
    parser.add_argument('-p', required=False, action='store_true', help='disable priors in the somatic calculation. Increases sensitivity for solid tumors')
    parser.add_argument('-J', required=False, action='store_true', help='Use prior probabilities accounting for the somatic mutation rate')
    parser.add_argument('-s', required=False, metavar='SOMATIC_PRIOR',type=float, help='prior probability of a somatic mutation (implies -J) [%(default)s]', default=0.01)
    #parser.add_argument('-T', required=False, metavar='THETA', type=float, help='theta in maq consensus calling model [%(default)s]', default=0.85)
    #parser.add_argument('-N', required=False, metavar='NUM_HAPLOTYPES', type=int, help='number of haplotypes in the sample [%(default)s]', default=2)
    parser.add_argument('-r', required=False, metavar='GERMLINE_PRIOR', type=float, help='prior of a difference between two haplotypes [%(default)s]', default=0.001)
    parser.add_argument('-n', required=True, metavar='NORMAL_NAME', help='normal sample id (for VCF header) [%(default)s]', default='NORMAL')
    parser.add_argument('-t', required=True, metavar='TUMOR_NAME', help='tumor sample id (for VCF header) [%(default)s]', default='TUMOR')
    parser.add_argument('-F', required=True, metavar='FORMAT', choices=('classic', 'vcf', 'bed'), help='select output format [%(default)s]', default='vcf')
    parser.add_argument('tumor_bam', help='tumor BAM file')
    parser.add_argument('normal_bam', help='normal BAM file')
    parser.add_argument('output', help='Output file')
 
    group = parser.add_argument_group('wrapper specific options')
    group.add_argument('--workdir', default='/tmp/', help='Working directory of the wrapper')
    group.add_argument('--reference-id', dest='reference_id', help='TCGA name of the reference', choices=['hg18', 'hg19', 'GRCh37', 'GRCh37-lite', '36', '36.1', '37'])
    group.add_argument('--tumor-uuid', dest='tumor_uuid', help='Tumor Sample uuid')
    group.add_argument('--tumor-barcode', dest='tumor_barcode', help='TCGA Tumor Sample barcode')
    group.add_argument('--tumor-accession', dest='tumor_accession', help='Tumor CGHub analysis id')
    group.add_argument('--tumor-platform', dest='tumor_platform', help='Tumor sequencing platform')
    group.add_argument('--normal-uuid', dest='normal_uuid', help='Normal Sample uuid')
    group.add_argument('--normal-barcode', dest='normal_barcode', help='TCGA Normal Sample barcode')
    group.add_argument('--normal-accession', dest='normal_accession', help='Normal CGHub analysis id')
    group.add_argument('--normal-platform', dest='normal_platform', help='Normal sequencing platform')
    group.add_argument('--individual', dest='individual', help='Individual barcode being analyzed')
    group.add_argument('--center', dest='center', help='Center name')
    group.add_argument('--sniper-exe', dest='sniper_exe', default='bam-somaticsniper', help='Normal CGHub analysis id')
    return parser

def wrapper_specific_arguments():
    return set(('f', 'tumor_bam', 'normal_bam', 'output', 'workdir', 'sniper_exe', 'reference_id', 'center',
        'tumor_uuid', 'tumor_barcode', 'tumor_accession', 'tumor_platform',
        'normal_uuid', 'normal_barcode', 'normal_accession', 'normal_platform',
        'individual'))

def create_sniper_opts(namespace_dict):
    args = []
    flag_opts = set(('L','G','p','J'))
    wrapper_arguments = wrapper_specific_arguments()
    
    for option, value in namespace_dict.items():
        if option in wrapper_arguments:
            continue
        if (option in flag_opts):
            if value:
                args.append("-" + option)
        else:
            args.append("-" + option + " " + str(value)) 
    return " ".join(args)

def create_sniper_cmdline(namespace_dict, reference, tumor_bam, normal_bam, temp_output_file):
    return " ".join([
        namespace_dict['sniper_exe'],
        create_sniper_opts(namespace_dict),
        "-f " + reference,
        tumor_bam, 
        normal_bam,
        temp_output_file])

def create_workspace(workdir, reference, tumor_bam, normal_bam):
    #need to symlink in stuff
    symlink_workspace_file(workdir, reference + ".fai" , "ref_genome.fasta.fai"),
    return ( 
            symlink_workspace_file(workdir, tumor_bam, "tumor.bam"),
            symlink_workspace_file(workdir, normal_bam, "normal.bam"),
            symlink_workspace_file(workdir, reference, "ref_genome.fasta"),
            )
    
def symlink_workspace_file(workdir, original_file, new_file):
    symlink_name = os.path.join(workdir, new_file)
    os.symlink(os.path.abspath(original_file), symlink_name)
    return symlink_name

def sample_from_bam(bam_file):
    file = pysam.Samfile(bam_file) #this may need to be AlignmentFile depending on pysam version
    readgroups = file.header['RG']
    samples = set()
    for readgroup in readgroups:
        samples.add(readgroup['SM'])
    assert(len(samples) == 1)
    return samples.pop()

def reheader_vcf(options, original_vcf_file, new_vcf_file):
    #need to add fields for the samples to the VCF header? as well as other info we add in Genome::Model::Tools::Vcf::Convert::Base.pm
    outfile = open(new_vcf_file, 'w')
    outfile.write("##fileformat=VCFv4.1\n")
    outfile.write("##fileDate=" + str(date.today().strftime("%Y%m%d")) + "\n")
    outfile.write("##tcgaversion=1.1.1\n")
    outfile.write(reference_header_line(options['reference_id'], options['f']))
    outfile.write("##phasing=none\n")
    outfile.write('##center="' + options['center'] + "\"\n")
    outfile.write(sample_header_line(
        options['t'],
        options['tumor_uuid'],
        options['tumor_barcode'],
        options['individual'],
        options['tumor_bam'],
        options['tumor_platform'],
        options['tumor_accession'],))
    outfile.write(sample_header_line(
        options['n'],
        options['normal_uuid'],
        options['normal_barcode'],
        options['individual'],
        options['normal_bam'],
        options['normal_platform'],
        options['normal_accession'],))
    outfile.write(vcfprocesslog_header_line(options))

    original = open(original_vcf_file, 'r')
    for line in original:
        if line.startswith("##"):
            if line.startswith("##FORMAT"):
                outfile.write(line)
        else:
            outfile.write(line)
    outfile.close()

def sample_header_line(sample_id, uuid, barcode, individual, file_name, platform, accession):
    return "##SAMPLE=<ID={0},SampleUUID={1},SampleTCGABarcode={2},Individual={3},File={4},Platform={5},Source=CGHub,Accession={6}>\n".format(sample_id, uuid, barcode, individual, os.path.basename(file_name), platform, accession)

def vcfprocesslog_header_line(options):
    input_vcf = "InputVCF=<.>"  #no VCF is put into this program so empty
    source = "InputVCFSource=<bam-somaticsniper>"
    version = "InputVCFVer=<1.0.4>" #Sniper version, could be different if exe specified differently
    param = 'InputVCFParam=<"' + create_sniper_opts(options) + '">'
    anno = "InputVCFgeneAnno=<.>"
    return "##vcfProcessLog=<" + ",".join([input_vcf, source, version, param, anno]) + ">\n"

def reference_header_line(reference_name, reference_path):
    reference_url = None
    reference_parse_result = urlparse(reference_path)
    if(reference_parse_result.scheme == ''):
        #assuming we've been given a file path
        reference_url = urlunparse(tuple(["file"]) + reference_parse_result[1:])
    else:
        reference_url = reference_path #could just unparse this, but that might reconstruct the string oddly

    header_line = "##reference=<ID=" + reference_name + ',Source=' + reference_url + ">\n"
    return header_line

def execute(options, ref, tumor, normal, temp_output_file):
    cmd = create_sniper_cmdline(options, ref, tumor, normal, temp_output_file)

    logging.info("RUNNING: %s" % (cmd))
    print "running", cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode

if __name__ == "__main__":
    parser = sniper_argparser()
    args = parser.parse_args()
    #opt_string = create_sniper_opts(vars(args))

    #tumor_sample_name = sample_from_bam(args.tumor_bam)
    #normal_sample_name = sample_from_bam(args.normal_bam)

    (workspace_tumor, workspace_normal, workspace_ref) = create_workspace(args.workdir, args.f, args.tumor_bam, args.normal_bam)
    temp_output_file = os.path.join(args.workdir, "raw_output.vcf")
    
    execute(vars(args), workspace_ref, workspace_tumor, workspace_normal, temp_output_file)
    reheader_vcf(vars(args), temp_output_file, args.output)

