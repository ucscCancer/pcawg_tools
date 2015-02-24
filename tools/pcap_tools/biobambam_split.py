#!/usr/bin/env python

import sys, os, re, argparse, shutil
import dateutil.parser
import header_utils
import traceback
import os
import hashlib
import os
import header_utils
import utils
import logging
import multiprocessing
import subprocess

import os
import urllib2
import logging
import dateutil.parser
from xml.etree import ElementTree

default_logger = logging.getLogger()
default_logger.setLevel(logging.DEBUG)
default_logger.addHandler(logging.StreamHandler())

FORCE_RUN=True

FORCE_RUN=True

class HeaderException(Exception):
    pass

def get_value(result, field, logger=default_logger):
    """Get the value of the field from the XML tree root"""

    if result.find(field) == None:
        msg = "missing field in XML: %s" % field
        logger.error(msg)
        raise HeaderException(msg)

    val = result.find(field).text
    if(val == None):
        #raise HeaderException("no value for field in XML: %s" % field)
        return "None"
    else:
        return val


def create_header(analysis_outdir, metadata, rg_dict, specimen_dict, logger=default_logger):
    """Generate headers in output dir for the new unaligned BAM named by RG ID"""

    rgid = rg_dict["ID"].replace(".", "_")
    header = "%s/header-%s.sam" %(analysis_outdir, rg_dict["ID"])
    header_file = open(header, "w")
    header_file.write("@HD\tVN:1.4\n")
    PI_STR = ""
    if len(rg_dict["PI"]):
        PI_STR="PI:%s\t" % (rg_dict["PI"])
    header_file.write("@RG\tID:%s:%s\tCN:%s\tPL:%s\tPM:%s\tLB:%s:%s:%s\t%sSM:%s\tPU:%s:%s\tDT:%s\n"
                %(metadata["center_name"], rgid,metadata["center_name"], metadata["platform"],metadata["platform_model"], metadata["seqtype"],
                metadata["center_name"], rg_dict["LB"], PI_STR, metadata["aliquot_id"], rg_dict["CN"], rg_dict["PU"], getUTCDate(rg_dict["DT"])))
    header_file.write("@CO\tdcc_project_code:%s-US\n" %metadata["disease"])
    header_file.write("@CO\tsubmitter_donor_id:%s\n" %metadata["participant_id"])
    header_file.write("@CO\tsubmitter_specimen_id:%s\n" %metadata["sample_id"])
    header_file.write("@CO\tsubmitter_sample_id:%s\n" %metadata["aliquot_id"])

    if metadata["sample_type"] not in specimen_dict:
        msg = "sample_type %s not found in specimen mapping" % metadata["sample_type"]
        logger.error(msg)
        if not FORCE_RUN:
            raise HeaderException(msg)

    if "sample_type" in metadata and metadata["sample_type"] in specimen_dict:
        (icgc_type, sample_class) = specimen_dict[metadata["sample_type"]]
    else:
        icgc_type = "unknown"
        sample_class = "unknown"

    #Sanity check about use_cntl
    if "use_cntl" in metadata:
        if metadata["use_cntl"] == "N/A" and sample_class == "tumour":
            msg = "Tumour sample requires use_cntl, set to %s. Are your IDs in the wrong order?" % metadata["use_cntl"]
            logger.error(msg)
            raise HeaderException(msg)
        if sample_class == "normal" and metadata["use_cntl"] != "N/A":
            msg = "Normal sample requires N/A use_cntl, set to %s. Are your IDs in the wrong order?" % metadata["use_cntl"]
            logger.error(msg)
            raise HeaderException(msg)

    header_file.write("@CO\tdcc_specimen_type:%s\n" % icgc_type)
    header_file.write("@CO\tuse_cntl:%s\n" %(metadata.get("use_cntl", "NA")))
    header_file.close()
    return header


def parse_cghub_metadata(analysis_id, logger=default_logger):
    """Get the metadata from the XML returned by cgquery"""

    xml_url = "https://cghub.ucsc.edu/cghub/metadata/analysisFull/%s" % analysis_id
    xml_resp = urllib2.urlopen(xml_url).read()

    root = ElementTree.fromstring(xml_resp)

    #error if no or > 1 results
    results = root.findall("Result")
    if len(results) == 0:
        msg = "no CGHub metadata found at %s" % xml_url
        logger.error(msg)
        raise HeaderException(msg)

    if len(results) > 1:
        msg = "more than one CGHub result found at %s" % xml_url
        logger.error(msg)
        raise HeaderException(msg)

    metadata = dict()
    for result in root.iter("Result"):
        metadata["participant_id"] = get_value(result, "participant_id")
        metadata["sample_id"] = get_value(result, "sample_id")
        metadata["disease"] = get_value(result, "disease_abbr")
        metadata["tss_id"] = get_value(result, "tss_id")
        metadata["seqtype"] = get_value(result, "library_strategy")
        metadata["analyte_code"] = get_value(result, "analyte_code")
        metadata["sample_type"] = get_value(result, "sample_type")
        metadata["platform"] = get_value(result, "platform")
        metadata["aliquot_id"] = get_value(result, "aliquot_id")
        metadata["center_name"] = get_value(result, "center_name")
        metadata["legacy_sample_id"] = get_value(result, "legacy_sample_id")
    metadata["platform_model"] = get_platform_model(root, metadata["platform"])

    return metadata

def get_platform_model(root, platform, logger=default_logger):
    """For the PM: code in the RG"""

    for result in root.iter("Result"):
        if(not(result == None)):
            for elem in result.iter("experiment_xml"):
                if(not(elem == None)):
                    for expr_set in elem.iter("EXPERIMENT_SET"):
                        if(not(expr_set == None)):
                            for plat in expr_set.iter(platform):
                                platform_model = get_value(plat, "INSTRUMENT_MODEL")
    return platform_model

def get_read_group_info(line, logger=default_logger):
    """Retrieve required read group information"""

    rg_dict = dict()
    #Initialize the dictionary, so we know if any fields are missing
    rg_dict["PI"] = ""
    rg_dict["CN"] = "UNKNOWN"
    rg_dict["ID"] = ""
    rg_dict["PL"] = "UNKNOWN"
    rg_dict["LB"] = ""
    rg_dict["SM"] = ""
    rg_dict["PU"] = ""
    rg_dict["DT"] = ""
    sline = line.split('\t')

    for item in sline:
        item = item.strip()

        if(item.startswith("ID:")):
            rg_dict["ID"] = item[3:]
        elif(item.startswith("PL:")):
            rg_dict["PL"] = item[3:]
        elif(item.startswith("PU:")):
            item = item.replace(".", "_") #to agree with ICGC SOP
            rg_dict["PU"] = item[3:]
        elif(item.startswith("LB:")):
            rg_dict["LB"] = item[3:]
        elif(item.startswith("DT:")):
            rg_dict["DT"] = item[3:]
        elif(item.startswith("SM:")):
            rg_dict["SM"] = item[3:]
        elif(item.startswith("CN:")):
            rg_dict["CN"] = item[3:]
        elif(item.startswith("PI:")):
            rg_dict["PI"] = item[3:]
        else:
            pass

    for key,value in rg_dict.items():
        if value == "":
            logger.warning("missing RG field %s" % key)

    return rg_dict

def get_cghub_xml(dirname, analysis_id, logger=default_logger):
    """Get the XML file using cgquery"""

    os.system("cgquery analysis_id=%s -a -o %s/metadata.xml" %(analysis_id, dirname))
    return "%s/metadata.xml" %(dirname)

def is_valid_analysis(info, logger=default_logger):
    """Check fields with controlled vocabulary"""

    if(info["center_name"] not in ["BCM", "BCCAGSC", "BI", "HMS-RK", "UNC-LCCC", "WUGSC", "USC-JHU"]):
        logger.error("The center %s is not in the defined center vocabulary" %(info.get("CN", "None")))
        return False
    if(info["platform"] not in ["CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "PACBIO"]):
        logger.error("The platform %s is not in the defined platform vocabulary" %(info.get("PL", "None")))
        return False
    if(info["platform_model"] not in ["Illumina Genome Analyzer II", "Illumina HiSeq", "Illumina HiSeq 2000", "Illumina HiSeq 2500"]):
        logger.error("The platform unit %s is not in the defined platform vocabulary" %(info.get("PU", "None")))
        return False
    return True

def rehead(input_dir, laneLevelBam, header, rgid, analysis_id, logger=default_logger):
    """
    Use samtools to reheader the unaligned BAMs,
    feel it's better to keep the metadata with the BAM itself as well
    """

    filename = "%s_%s.cleaned.bam" %(analysis_id, rgid)
    output_file = os.path.join(input_dir, filename)
    os.system("samtools reheader %s %s > %s" %(header, laneLevelBam,  output_file))
    if not os.path.exists(output_file):
        msg = "after reheader output_file does not exist: %s" % output_file
        logger.error(msg)
        raise HeaderException(msg)

    return output_file

def getUTCDate(dateString, logger=default_logger):
    """Actual UTC with timezone does not validate, change the time to Z"""

    if(dateString):
        date = dateutil.parser.parse(dateString)
        date = date.astimezone(dateutil.tz.tzutc())
        return date.strftime('%Y-%m-%dT%H:%M:%SZ')
    return dateString



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

def process_rg(analysis_id, rg_dict, header, work_dir, output_dir, logger=default_logger):

    valid_bam, lane_level_bam = lane_level_bam_from_fastq(rg_dict, work_dir)
    if not valid_bam:
        return False, None

    final_file = header_utils.rehead(output_dir, lane_level_bam, header, rg_dict["ID"], analysis_id)
    md5_file = open(final_file + '.md5', "w")
    md5_file.write(utils.calc_md5sum(final_file))
    md5_file.close()
    return True, final_file


def gen_unaligned_bam(bam_filename, analysis_id, metadata, specimen_dict, work_dir, output_dir, num_processes=4, logger=default_logger ):
    """
    The bulk of the work, calls splitting, generates new headers, generates initial
    unaligned BAM, reheaders with new headers
    """

    read_group_sam = os.path.join(output_dir, 'rg_header.sam')

    #get the read groups from the original sample level BAM
    exit_code = os.system("samtools view -H %s | grep \"@RG\" > %s" %(bam_filename, read_group_sam))
    if exit_code != 0:
        print "Failure in bam splitting during read group extraction from %s" % bam_filename
        return 1


    rg_file = open(read_group_sam, "r")

    #create the read group fastqs
    try:
        cmd = "bamtofastq outputperreadgroup=1 gz=1 level=1 inputbuffersize=2097152000 tryoq=1 outputdir=%s T=`mktemp -p %s bamtofastq_XXXXXXXXX` < %s" %(work_dir, work_dir, bam_filename)
        logger.info("Running %s" % cmd)
        subprocess.check_call(cmd, shell=True)
    except:
        print "Failure in bam splitting"
        return 1


    if header_utils.is_valid_analysis(metadata) or FORCE_RUN:
        pool = multiprocessing.Pool(processes=num_processes)
        results = []
        for line in rg_file:
            rg_dict = header_utils.get_read_group_info(line)
            header = header_utils.create_header(output_dir, metadata, rg_dict, specimen_dict)
            r = pool.apply_async(process_rg, (analysis_id, rg_dict, header, work_dir, output_dir))
            results.append(r)

        rg_file.close()

        out = []
        for r in results:
            out.append(r.get())

        utils.clean_up_dir(output_dir)
        if not all( a[0] for a in out ):
            #one of the read group bamtofastq failed
            return 1
        with open(os.path.join(output_dir, "results.list"), "w") as out_handle:
            for ok, file_name in out:
                out_handle.write("%s\n" % (file_name))

    else:
        print "Invalid header/metadata for BAM" % bam_filename
        return 1
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

    parser = argparse.ArgumentParser(prog='biobambam_split.py', description='Create unaligned BAM files')
    parser.add_argument('--bam_path', type=str, help='path/to/tcga/data/labeled_by_analysis_id', required=True)
    parser.add_argument('--normal_id', type=str, help='UUID for normal analysis (analysis_id)', default=None)
    parser.add_argument('--tumor_id', type=str, help='Comma separated list of tumor analysis UUIDs (analysid_id(s))', default=None)
    parser.add_argument('--work_dir', type=str, help='path/to/output/directory', default=None)
    parser.add_argument('--output_dir', type=str, help='path/to/output/directory', required=True)

    args = parser.parse_args()

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
