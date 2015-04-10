#!/usr/bin/env python

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import logging
import glob
import gzip
import os.path
from itertools import izip

'''
'   Amie Radenbaugh - 02/09/2015
'   UCSC - RADIA
'   Program name: "mageTab.py"
'''

def get_read_fileHandler(aFilename):
    '''
    ' Open aFilename for reading and return
    ' the file handler.  The file can be
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'rb')
    else:
        return open(aFilename,'r')


def get_write_fileHandler(aFilename):
    '''
    ' Open aFilename for writing and return
    ' the file handler.  The file can be
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'wb')
    else:
        return open(aFilename,'w')


def createSDRFFilename(aVCFDir, aProtectedIDFFilename, aProtectedSDRFFilename, aProtectedArchiveName, anOpenAccessIDFFilename, anOpenAccessSDRFFilename, anOpenAccessArchiveName, aCmdLineOptionsDict, anIsDebug):

    protectedSDRFFileHandler = get_write_fileHandler(aProtectedSDRFFilename)
    openAccessSDRFFileHandler = get_write_fileHandler(anOpenAccessSDRFFilename)

    # output the header lines
    protectedHeaderList = ["Extract Name", "Comment [TCGA Barcode]", "Comment [is tumor]", "Material Type", "Annotation REF", "Comment [TCGA Genome Reference]"]
    # some default ones that will be empty for us
    protectedHeaderList += ["Protocol REF", "Parameter Value [Vendor]", "Parameter Value [Catalog Name]", "Parameter Value [Catalog Number]"]
    # protocol for bam file name
    protectedHeaderList += ["Protocol REF", "Comment [Derived Data File REF]", "Comment [TCGA CGHub ID]", "Comment [TCGA Include for Analysis]"]

    # so far the protected and open access headers are the same
    openAccessHeaderList = list(protectedHeaderList)

    # protocol somatic_variant_detection
    protectedHeaderList += ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]", "Comment [TCGA Include for Analysis]"]
    protectedHeaderList += ["Comment [TCGA Data Type]", "Comment [TCGA Data Level]", "Comment [TCGA Archive Name]"]
    # protocol vcf2maf for protected access MAF
    protectedHeaderList += ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]", "Comment [TCGA Include for Analysis]"]
    protectedHeaderList += ["Comment [TCGA Data Type]", "Comment [TCGA Data Level]", "Comment [TCGA Archive Name]"]

    # protocol vcf2maf for open access MAF
    openAccessHeaderList += ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]", "Comment [TCGA Include for Analysis]"]
    openAccessHeaderList += ["Comment [TCGA Data Type]", "Comment [TCGA Data Level]", "Comment [TCGA Archive Name]"]

    # add the header lines
    protectedSDRFFileHandler.write("\t".join(protectedHeaderList) + "\n")
    openAccessSDRFFileHandler.write("\t".join(openAccessHeaderList) + "\n")

    # hard-coded values
    emptyValue = "->"
    mafVersion = "2.4.1"

    # for each vcf file
    for vcfFile in (glob.glob(aVCFDir + "TCGA*_TCGA*.vcf*")):

        # open the file
        vcfFileHandler = get_read_fileHandler(vcfFile)

        rnaFasta = "GRCh37-lite"
        dnaFasta = None
        dnaNormalBarcode = None
        dnaNormalUUID = None
        dnaNormalBamFilename = None
        dnaNormalCgHubId = None
        dnaTumorBarcode = None
        dnaTumorUUID = None
        dnaTumorBamFilename = None
        dnaTumorCgHubId = None
        rnaNormalBarcode = None
        rnaNormalUUID = None
        rnaNormalBamFilename = None
        rnaNormalCgHubId = None
        rnaTumorBarcode = None
        rnaTumorUUID = None
        rnaTumorBamFilename = None
        rnaTumorCgHubId = None

        for line in vcfFileHandler:

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            if (anIsDebug):
                logging.debug("vcfLine: %s", line)

            # if it is an empty line, then just continue
            if (line.isspace()):
                continue;

            # we need to extract the tcga spec that was used
            elif (line.startswith("##tcgaversion")):
                ##tcgaversion=1.0

                (key, value) = line.split("=")
                tcgaSpecVersion = value

            # we need to extract info from the reference tag in the header
            elif (line.startswith("##reference")):
                ##reference=<ID=GRCh37,Source=file:/inside/depot/fa/Homo_sapiens_assembly19.fasta>
                ##assembly=file:/inside/depot/fa/Homo_sapiens_assembly19.fasta

                if (anIsDebug):
                    logging.debug(line)

                if ("GRCh37-lite" in line):
                    dnaFasta = "GRCh37-lite"
                elif ("GRCh37" in line):
                    dnaFasta = "GRCh37"
                elif ("Homo_sapiens_assembly19" in line):
                    dnaFasta = "Homo_sapiens_assembly19.fasta"
                elif ("hg19" in line):
                    dnaFasta = "hg19"
                elif ("hg18" in line):
                    dnaFasta = "hg18"
                elif ("NCBI36" in line):
                    dnaFasta = "36.1"
                elif ("NCBI37" in line):
                    dnaFasta = "37"

            # we need to extract info from the SAMPLE tag in the header
            elif (line.startswith("##SAMPLE")):

                ##SAMPLE=<ID=DNA_NORMAL,SampleUUID=5218f2a6-5b5d-4e1a-aaf1-a0cd5e536571,SampleTCGABarcode=TCGA-IB-7646-10A-01D-2154-08,Individual=TCGA-IB-7646,
                #Description="Blood Derived Normal DNA",File="/inside/depot4/users/aradenba/data/hg19/paad/bams/wxs/TCGA-IB-7646-10A-01D-2154-08.bam",
                #Platform="Illumina HiSeq 2000",Source=CGHub,Accession=a6fd3469-b442-416d-aae9-991092fcb579,SequenceSource=WXS>
                ##SAMPLE=<ID=DNA_TUMOR,SampleUUID=e26b2473-9399-4eb8-9574-5b934b560740,SampleTCGABarcode=TCGA-IB-7646-01A-11D-2154-08,Individual=TCGA-IB-7646,
                #Description="Primary solid Tumor DNA",File="/inside/depot4/users/aradenba/data/hg19/paad/bams/wxs/TCGA-IB-7646-01A-11D-2154-08.bam",
                #Platform="Illumina HiSeq 2000",Source=CGHub,Accession=58287578-15b7-41d2-a5ac-70aaeac1fa40,SequenceSource=WXS>
                ##SAMPLE=<ID=RNA_TUMOR,SampleUUID=c3d463e2-8625-403f-af7b-63029bcb6eab,SampleTCGABarcode=TCGA-IB-7646-01A-11R-2156-07,Individual=TCGA-IB-7646,
                #Description="Primary solid Tumor RNA",File="/inside/depot4/users/aradenba/data/hg19/paad/bams/rna/TCGA-IB-7646-01A-11R-2156-07.bam",
                #Platform="Illumina HiSeq 2000",Source=CGHub,Accession=f4b705d6-29e2-4b83-9a99-3ef36d23d39b,SequenceSource=RNA-Seq>

                sampleLine = line[0:(len(line)-1)]
                sampleLine = sampleLine[len("##SAMPLE=<"):len(sampleLine)]
                sampleParamsList = sampleLine.split(",")
                sampleType = -1
                sampleBarcode = ""
                sampleUUID = ""
                sampleBamFilename = ""
                sampleCgHubId = ""
                isRNA = False

                # create a dictionary of existing params
                for param in sampleParamsList:
                    if ("software" in param):
                        continue;

                    (key, value) = param.split("=")

                    # remove the quotes
                    if (value.startswith("\"")):
                        value = value[1:(len(value)-1)]

                    # keep track of the barcode, uuid, and sequence source
                    if (key == "SampleTCGABarcode"):
                        sampleBarcode = value
                        sampleBarcodeList = sampleBarcode.split("-")
                        sampleType = int(sampleBarcodeList[3][:2])
                        analyteType = sampleBarcodeList[4][2:]
                        if (analyteType == "R"):
                            isRNA = True
                    elif (key == "SampleUUID"):
                        sampleUUID = value
                    elif (key == "File"):
                        sampleBamFilename = value
                    elif (key == "Accession"):
                        sampleCgHubId = value

                if (anIsDebug):
                    logging.debug("SampleTCGABarcode=%s, SampleUUID=%s, SampleFile=%s, CGHubId=%s", sampleBarcode, sampleUUID, sampleBamFilename, sampleCgHubId)

                # the sample types from 0-9 are tumor and 10-19 are normal, 20 and above are control samples
                if (sampleType != -1):
                    if (sampleType < 10):
                        if (isRNA):
                            rnaTumorBarcode = sampleBarcode
                            rnaTumorUUID = sampleUUID
                            rnaTumorBamFilename = sampleBamFilename
                            rnaTumorCgHubId = sampleCgHubId
                        else:
                            dnaTumorBarcode = sampleBarcode
                            dnaTumorUUID = sampleUUID
                            dnaTumorBamFilename = sampleBamFilename
                            dnaTumorCgHubId = sampleCgHubId
                    elif (sampleType >= 10 and sampleType < 20):
                        if (isRNA):
                            rnaNormalBarcode = sampleBarcode
                            rnaNormalUUID = sampleUUID
                            rnaNormalBamFilename = sampleBamFilename
                            rnaNormalCgHubId = sampleCgHubId
                        else:
                            dnaNormalBarcode = sampleBarcode
                            dnaNormalUUID = sampleUUID
                            dnaNormalBamFilename = sampleBamFilename
                            dnaNormalCgHubId = sampleCgHubId
                    else:
                        logging.critical("Traceback:  Unexpected sample type %s", sampleType)
                        sys.exit(1)

                    if (anIsDebug):
                        logging.debug("normalBarcode=%s, normalUUID=%s, tumorBarcode=%s, tumorUUID=%s, rnaNormalBarcode=%s, rnaNormalUUID=%s, rnaTumorBarcode=%s, rnaTumorUUID=%s", dnaNormalBarcode, dnaNormalUUID, dnaTumorBarcode, dnaTumorUUID, rnaNormalBarcode, rnaNormalUUID, rnaTumorBarcode, rnaTumorUUID)

            # these are other header lines that we don't need, so just continue
            elif (line.startswith("#")):
                continue;

            # now we are to the data
            else:

                # get some variables from the command line params or defaults
                mafFilenameProtectedAccess = aProtectedArchiveName + ".protected.maf"
                mafFilenameOpenAccess = anOpenAccessArchiveName + ".somatic.maf"

                protocolSomaticVariants = None
                protocolVCFToMAF = None
                protocolList = aCmdLineOptionsDict["protocolNames"].split(",")
                for protocol in protocolList:
                    if "variant_calling" in protocol:
                        protocolSomaticVariants = protocol
                    if "vcf2maf" in protocol:
                        protocolVCFToMAF = protocol

                vcfFilename = os.path.basename(vcfFile)

                # output one line per patient sample
                # normal and tumor DNA are required
                # check if RNA was available

                #################
                # Normal DNA
                #################
                # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
                dnaNormalProtectedOutputList = [dnaNormalUUID, dnaNormalBarcode, "no", "DNA", emptyValue, dnaFasta]
                # protocolRef:default, vendor, catalogName, CatalogNumber
                dnaNormalProtectedOutputList += [emptyValue, emptyValue, emptyValue, emptyValue]
                # protocolRef:bamfile, bamfile, cgHubId, used?
                dnaNormalProtectedOutputList += [emptyValue, dnaNormalBamFilename, dnaNormalCgHubId, "yes"]

                # so far the protected and open access results are the same
                dnaNormalOpenAccessOutputList = list(dnaNormalProtectedOutputList)

                # back to the protected
                # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
                dnaNormalProtectedOutputList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]
                # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                dnaNormalProtectedOutputList += [protocolVCFToMAF, mafFilenameProtectedAccess, mafVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]

                # back to the open access
                # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                dnaNormalOpenAccessOutputList += [protocolVCFToMAF, mafFilenameOpenAccess, mafVersion, "yes", "Mutations", "Level 2", anOpenAccessArchiveName]

                # write to file
                protectedSDRFFileHandler.write("\t".join(dnaNormalProtectedOutputList) + "\n")
                openAccessSDRFFileHandler.write("\t".join(dnaNormalOpenAccessOutputList) + "\n")

                #################
                # Tumor DNA
                #################
                # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
                dnaTumorProtectedOutputList = [dnaTumorUUID, dnaTumorBarcode, "yes", "DNA", emptyValue, dnaFasta]
                # protocolRef:default, vendor, catalogName, CatalogNumber
                dnaTumorProtectedOutputList += [emptyValue, emptyValue, emptyValue, emptyValue]
                # protocolRef:bamfile, bamfile, cgHubId, used?
                dnaTumorProtectedOutputList += [emptyValue, dnaTumorBamFilename, dnaTumorCgHubId, "yes"]

                # so far the protected and open access results are the same
                dnaTumorOpenAccessOutputList = list(dnaTumorProtectedOutputList)

                # back to the protected
                # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
                dnaTumorProtectedOutputList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]
                # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                dnaTumorProtectedOutputList += [protocolVCFToMAF, mafFilenameProtectedAccess, mafVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]

                # back to the open access
                # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                dnaTumorOpenAccessOutputList += [protocolVCFToMAF, mafFilenameOpenAccess, mafVersion, "yes", "Mutations", "Level 2", anOpenAccessArchiveName]
                protectedSDRFFileHandler.write("\t".join(dnaTumorProtectedOutputList) + "\n")
                openAccessSDRFFileHandler.write("\t".join(dnaTumorOpenAccessOutputList) + "\n")

                #################
                # Normal RNA
                #################
                if (rnaNormalBarcode != None):
                    # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
                    rnaNormalProtectedOutputList = [rnaNormalUUID, rnaNormalBarcode, "no", "RNA", emptyValue, rnaFasta]
                    # protocolRef:default, vendor, catalogName, CatalogNumber
                    rnaNormalProtectedOutputList += [emptyValue, emptyValue, emptyValue, emptyValue]
                    # protocolRef:bamfile, bamfile, cgHubId, used?
                    rnaNormalProtectedOutputList += [emptyValue, rnaNormalBamFilename, rnaNormalCgHubId, "yes"]

                    # so far the protected and open access results are the same
                    rnaNormalOpenAccessOutputList = list(rnaNormalProtectedOutputList)

                    # back to the protected
                    # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
                    rnaNormalProtectedOutputList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]
                    # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                    rnaNormalProtectedOutputList += [protocolVCFToMAF, mafFilenameProtectedAccess, mafVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]

                    # back to the open access
                    # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                    rnaNormalOpenAccessOutputList += [protocolVCFToMAF, mafFilenameOpenAccess, mafVersion, "yes", "Mutations", "Level 2", anOpenAccessArchiveName]
                    protectedSDRFFileHandler.write("\t".join(rnaNormalProtectedOutputList) + "\n")
                    openAccessSDRFFileHandler.write("\t".join(rnaNormalOpenAccessOutputList) + "\n")

                #################
                # Tumor RNA
                #################
                if (rnaTumorBarcode != None):
                    # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
                    rnaTumorProtectedOutputList = [rnaTumorUUID, rnaTumorBarcode, "yes", "RNA", emptyValue, rnaFasta]
                    # protocolRef:default, vendor, catalogName, CatalogNumber
                    rnaTumorProtectedOutputList += [emptyValue, emptyValue, emptyValue, emptyValue]
                    # protocolRef:bamfile, bamfile, cgHubId, used?
                    rnaTumorProtectedOutputList += [emptyValue, rnaTumorBamFilename, rnaTumorCgHubId, "yes"]

                    # so far the protected and open access results are the same
                    rnaTumorOpenAccessOutputList = list(rnaTumorProtectedOutputList)

                    # back to the protected
                    # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
                    rnaTumorProtectedOutputList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]
                    # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                    rnaTumorProtectedOutputList += [protocolVCFToMAF, mafFilenameProtectedAccess, mafVersion, "yes", "Mutations", "Level 2", aProtectedArchiveName]

                    # back to the open access
                    # protocol:vcf2maf, mafFilename, tcgaSpecVerion, used?, dataType, dataLevel, archiveName
                    rnaTumorOpenAccessOutputList += [protocolVCFToMAF, mafFilenameOpenAccess, mafVersion, "yes", "Mutations", "Level 2", anOpenAccessArchiveName]
                    protectedSDRFFileHandler.write("\t".join(rnaTumorProtectedOutputList) + "\n")
                    openAccessSDRFFileHandler.write("\t".join(rnaTumorOpenAccessOutputList) + "\n")

                break;

        # close the file
        vcfFileHandler.close()

    # close the file
    protectedSDRFFileHandler.close()
    openAccessSDRFFileHandler.close()

    return


def createIDFFilename(aProtectedIDFFilename, aProtectedSDRFFilename, anOpenAccessIDFFilename, anOpenAccessSDRFFilename, aCmdLineOptionsDict, anIsDebug):

    # create two files:  one for the protected archive and one for the open-access archive
    idfFilenames = [aProtectedIDFFilename, anOpenAccessIDFFilename]
    sdrfFilenames = [aProtectedSDRFFilename, anOpenAccessSDRFFilename]

    # for each idf and sdrf file
    for (idfFilename, sdrfFilename) in izip(idfFilenames, sdrfFilenames):
        idfFileHandler = get_write_fileHandler(idfFilename)

        # output the experimental design lines
        idfFileHandler.write("\t".join(["Investigation Title", "Analysis of TCGA " + aCmdLineOptionsDict["disease"] + " Whole-Exome Sequencing (WES) and RNA-Seq data"]) + "\n")
        idfFileHandler.write("\t".join(["Experimental Design", aCmdLineOptionsDict["expDesign"]]) + "\n")
        idfFileHandler.write("\t".join(["Experimental Design Term Source REF", aCmdLineOptionsDict["expDesignOntology"]]) + "\n")
        idfFileHandler.write("\t".join(["Experimental Factor Name", aCmdLineOptionsDict["expDesignFactorName"]]) + "\n")
        idfFileHandler.write("\t".join(["Experimental Factor Type", aCmdLineOptionsDict["expDesignFactorType"]]) + "\n")
        idfFileHandler.write("\n")

        # output the person lines
        idfFileHandler.write("\t".join(["Person Last Name", aCmdLineOptionsDict["personLastName"]]) + "\n")
        idfFileHandler.write("\t".join(["Person First Name", aCmdLineOptionsDict["personFirstName"]]) + "\n")
        idfFileHandler.write("\t".join(["Person Mid Initials", aCmdLineOptionsDict["personMidInitial"]]) + "\n")
        idfFileHandler.write("\t".join(["Person Email", aCmdLineOptionsDict["personEmail"]]) + "\n")
        idfFileHandler.write("\t".join(["Person Address", aCmdLineOptionsDict["personAddress"]]) + "\n")
        idfFileHandler.write("\t".join(["Person Affiliation", aCmdLineOptionsDict["personAffiliation"]]) + "\n")
        idfFileHandler.write("\t".join(["Person Roles", aCmdLineOptionsDict["personRole"]]) + "\n")
        idfFileHandler.write("\n")

        # output the publication lines
        idfFileHandler.write("\t".join(["PubMed ID", aCmdLineOptionsDict["pubMedId"]]) + "\n")
        idfFileHandler.write("\t".join(["Publication Author List", aCmdLineOptionsDict["pubAuthors"]]) + "\n")
        idfFileHandler.write("\t".join(["Publication Title", aCmdLineOptionsDict["pubTitle"]]) + "\n")
        idfFileHandler.write("\t".join(["Publication Status", aCmdLineOptionsDict["pubStatus"]]) + "\n")
        idfFileHandler.write("\t".join(["Experiment Description", aCmdLineOptionsDict["expDescription"]]) + "\n")
        idfFileHandler.write("\n")

        # output the protocol lines
        idfFileHandler.write("\t".join(["Protocol Name", "\t".join(aCmdLineOptionsDict["protocolNames"].split(","))]) + "\n")
        idfFileHandler.write("\t".join(["Protocol Type", "\t".join(aCmdLineOptionsDict["protocolTypes"].split(","))]) + "\n")
        idfFileHandler.write("\t".join(["Protocol Description", "\t".join(aCmdLineOptionsDict["protocolDescriptions"].split(","))]) + "\n")
        idfFileHandler.write("\t".join(["Protocol Term Source REF", "\t".join(aCmdLineOptionsDict["protocolOntologies"].split(","))]) + "\n")
        idfFileHandler.write("\t".join(["Protocol Parameters", "\t".join(aCmdLineOptionsDict["protocolParameters"].split(","))]) + "\n")
        idfFileHandler.write("\n")

        # output the sdrf line
        sdrfBasename = os.path.basename(sdrfFilename)
        idfFileHandler.write("\t".join(["SDRF Files", sdrfBasename]) + "\n")
        idfFileHandler.write("\n")

        # output the ontology lines
        idfFileHandler.write("\t".join(["Term Source Name", aCmdLineOptionsDict["ontologyName"]]) + "\n")
        idfFileHandler.write("\t".join(["Term Source File", aCmdLineOptionsDict["ontologyFile"]]) + "\n")
        idfFileHandler.write("\t".join(["Term Source Version", aCmdLineOptionsDict["ontologyVersion"]]) + "\n")

        # close the file
        idfFileHandler.close()

    return


def main():

    # create the usage statement
    usage = "usage: python %prog vcfDir protectedIDFFilename protectedSDRFFilename protectedArchiveName openAccessIDFFilename openAccessSDRFFilename openAccessArchiveName [Options]"
    i_cmdLineParser = OptionParser(usage=usage)

    i_cmdLineParser.add_option("", "--expDesign", dest="expDesign", default="individual_genetic_characteristics_design", metavar="EXP_DESIGN", help="the experimental design tag")
    i_cmdLineParser.add_option("", "--expDesignOntology", dest="expDesignOntology", default="MGED Ontology", metavar="EXP_DESIGN_ONTOLOGY", help="the experimental design term source reference (e.g. MGED, HUGO, etc.")
    i_cmdLineParser.add_option("", "--expDesignFactorName", dest="expDesignFactorName", default="somatic_variant", metavar="EXP_DESIGN_FACTOR_NAME", help="the experimental design factor name")
    i_cmdLineParser.add_option("", "--expDesignFactorType", dest="expDesignFactorType", default="disease_state", metavar="EXP_DESIGN_FACTOR_TYPE", help="the experimental design factor type")

    i_cmdLineParser.add_option("", "--personLastName", dest="personLastName", default="Radenbaugh", metavar="LAST_NAME", help="the submitter's last name")
    i_cmdLineParser.add_option("", "--personFirstName", dest="personFirstName", default="Amie", metavar="FIRST_NAME", help="the submitter's first name")
    i_cmdLineParser.add_option("", "--personMidInitial", dest="personMidInitial", default="J", metavar="MID_INITIAL", help="the submitter's middle initial")
    i_cmdLineParser.add_option("", "--personEmail", dest="personEmail", default="aradenba@soe.ucsc.edu", metavar="EMAIL", help="the submitter's email")
    i_cmdLineParser.add_option("", "--personAddress", dest="personAddress", default="University of California Santa Cruz, 1156 High St, Mail Stop CBSE, Santa Cruz, CA 95064, USA", metavar="ADDRESS", help="the submitter's institutional address")
    i_cmdLineParser.add_option("", "--personAffiliation", dest="personAffiliation", default="University of California Santa Cruz Genome Institute", metavar="AFFILIATION", help="the submitter's affiliation")
    i_cmdLineParser.add_option("", "--personRole", dest="personRole", default="submitter", metavar="ROLE", help="the submitter's role")

    i_cmdLineParser.add_option("", "--pubMedId", dest="pubMedId", default="25405470", metavar="PUB_MED_ID", help="the PubMed ID")
    i_cmdLineParser.add_option("", "--pubAuthors", dest="pubAuthors", default="Radenbaugh AJ, Ma S, Ewing A, Stuart JM, Collisson EA, Zhu J, Haussler D", metavar="PUB_AUTHORS", help="the publication author list")
    i_cmdLineParser.add_option("", "--pubTitle", dest="pubTitle", default="RADIA: RNA and DNA Integrated Analysis for Somatic Mutation Detection", metavar="PUB_TITLE", help="the publication title")
    i_cmdLineParser.add_option("", "--pubStatus", dest="pubStatus", default="published", metavar="PUB_STATUS", help="the publication status")
    i_cmdLineParser.add_option("", "--expDescription", dest="expDescription", default="Detection of somatic variants from the TCGA Whole-Exome Sequencing (WES) and RNA-Seq data using RADIA", metavar="EXP_DESCRIPTION", help="the description of the experiment")

    i_cmdLineParser.add_option("", "--protocolNames", dest="protocolNames", default="ucsc.edu:variant_calling:Illumina_DNASeq:01,ucsc.edu:vcf2maf:Illumina_DNASeq:01", metavar="PROTOCOL_NAMES", help="the protocol names")
    i_cmdLineParser.add_option("", "--protocolTypes", dest="protocolTypes", default="Variant Calling,MAF Generation", metavar="PROTOCOL_TYPES", help="the protocol types")
    i_cmdLineParser.add_option("", "--protocolDescriptions", dest="protocolDescriptions", default="Somatic Variant Calling Pipeline: RADIA v1.1.1 (https://github.com/aradenbaugh/radia/),Annotation: SnpEff v3.3_GRCh37.69 (http://snpeff.sourceforge.net/)", metavar="PROTOCOL_DESCRIPTIONS", help="the protocol descriptions")
    i_cmdLineParser.add_option("", "--protocolOntologies", dest="protocolOntologies", default="MGED Ontology,MGED Ontology", metavar="PROTOCOL_ONTOLOGIES", help="the protocol ontologies")
    i_cmdLineParser.add_option("", "--protocolParameters", dest="protocolParameters", default="", metavar="PROTOCOL_PARAMS", help="the protocol parameters")

    i_cmdLineParser.add_option("", "--ontologyName", dest="ontologyName", default="MGED Ontology", metavar="ONTOLOGY_NAME", help="the ontology name")
    i_cmdLineParser.add_option("", "--ontologyFile", dest="ontologyFile", default="http://mged.sourceforge.net/ontologies/MGEDontology.php", metavar="ONTOLOGY_FILE", help="the ontology file")
    i_cmdLineParser.add_option("", "--ontologyVersion", dest="ontologyVersion", default="1.3.1.1", metavar="ONTOLOGY_VERSION", help="the ontology version")

    i_cmdLineParser.add_option("-d", "--disease", dest="disease", default="", metavar="DISEASE", help="a disease abbreviation (e.g. BRCA), will be taken from VCF header if available")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(6,37,1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_vcfDir = str(i_cmdLineArgs[0])
    i_protectedIDFFilename = str(i_cmdLineArgs[1])
    i_protectedSDRFFilename = str(i_cmdLineArgs[2])
    i_protectedArchiveName = str(i_cmdLineArgs[3])
    i_openAccessIDFFilename = str(i_cmdLineArgs[4])
    i_openAccessSDRFFilename = str(i_cmdLineArgs[5])
    i_openAccessArchiveName = str(i_cmdLineArgs[6])
    writeFilenameList = [i_protectedIDFFilename, i_protectedSDRFFilename, i_openAccessIDFFilename, i_openAccessSDRFFilename]

    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_logFilename = None
    # try to get any optional parameters with no defaults
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        writeFilenameList += [i_logFilename]

    # assuming loglevel is bound to the string value obtained from the
    # command line argument. Convert to upper case to allow the user to
    # specify --log=DEBUG or --log=debug
    i_numericLogLevel = getattr(logging, i_logLevel.upper(), None)
    if not isinstance(i_numericLogLevel, int):
        raise ValueError("Invalid log level: '%s' must be one of the following:  DEBUG, INFO, WARNING, ERROR, CRITICAL", i_logLevel)

    # set up the logging
    if (i_logFilename != None):
        logging.basicConfig(level=i_numericLogLevel, filename=i_logFilename, filemode='w', format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=i_numericLogLevel, format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # set the debug
    i_debug = (i_numericLogLevel < logging.WARNING)

    # do some debugging
    if (i_debug):
        logging.debug("i_vcfDir=%s", i_vcfDir)
        logging.debug("i_protectedIDFFilename=%s", i_protectedIDFFilename)
        logging.debug("i_protectedSDRFFilename=%s", i_protectedSDRFFilename)
        logging.debug("i_openAccessIDFFilename=%s", i_openAccessIDFFilename)
        logging.debug("i_openAccessSDRFFilename=%s", i_openAccessSDRFFilename)
        logging.debug("i_logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)

    # check for any errors
    readFilenameList = []

    i_cmdLineOptionsDict = vars(i_cmdLineOptions)
    createIDFFilename(i_protectedIDFFilename, i_protectedSDRFFilename, i_openAccessIDFFilename, i_openAccessSDRFFilename, i_cmdLineOptionsDict, i_debug)
    createSDRFFilename(i_vcfDir, i_protectedIDFFilename, i_protectedSDRFFilename, i_protectedArchiveName, i_openAccessIDFFilename, i_openAccessSDRFFilename, i_openAccessArchiveName, i_cmdLineOptionsDict, i_debug)

    return

main()
sys.exit(0)
