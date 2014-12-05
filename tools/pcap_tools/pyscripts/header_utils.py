
import os
import urllib2
import logging
import dateutil.parser
from xml.etree import ElementTree

default_logger = logging.getLogger(name='pcap_split')

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
