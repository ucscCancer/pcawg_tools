#!/usr/bin/env python

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import rnaEditingUtil               # utility functions for rna editing
import logging
import glob
import re
import collections
from itertools import izip
import gzip

'''
'   Amie Radenbaugh - 02/10/2014
'   UCSC - RADIA  
'   Program name: "vcfToMaf.py"
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
    

def get_entrez_id(anEntrezDbDict, anEnsemblDbDict, aHugoSymbol, anEnsemblId, anIsDebug):
    
    if (aHugoSymbol in anEntrezDbDict):
        return anEntrezDbDict[aHugoSymbol]
    elif (anEnsemblId in anEnsemblDbDict):
        return anEnsemblDbDict[anEnsemblId]
    return "0"


def load_entrez_db(anEntrezDbFilename, anIsDebug):

    entrezDbDict = dict()
    ensemblDbDict = dict()
    entrezFileHandler = get_read_fileHandler(anEntrezDbFilename)
    
    for line in entrezFileHandler:
        # here is the header and example line:
        # Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description type_of_gene Symbol_from_nomenclature_authority 
        # Full_name_from_nomenclature_authority Nomenclature_status Other_designations Modification_date
        # 9606    1    A1BG    -    A1B|ABG|GAB|HYST2477    HGNC:5|MIM:138670|Ensembl:ENSG00000121410|HPRD:00726|Vega:OTTHUMG00000183507    19    19q13.4    alpha-1-B glycoprotein
        # protein-coding    A1BG    alpha-1-B glycoprotein    O    HEL-S-163pA|alpha-1B-glycoprotein|epididymis secretory sperm binding protein Li 163pA    20140209
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        #if (anIsDebug):
        #    logging.debug("entrezLine: %s", line)
            
        # if it is an empty line or comment, then just continue
        if (line.isspace() or line.startswith("#")):
            continue;
        else:
            # split the line on the tab
            splitLine = line.split("\t")
            
            # parse the fields
            #taxId = splitLine[0]
            entrezId = splitLine[1] 
            hugoSymbol = splitLine[2]
            entrezDbDict[hugoSymbol] = entrezId
            #locusTag = splitLine[3]
            #synonymsList = splitLine[4].split("|")
            
            # if there are no ids, then it's just "-"
            dbIdsList = splitLine[5].split("|")
            if (len(dbIdsList) == 1 and dbIdsList[0] == "-"):
                continue;
            
            for dbAndId in dbIdsList:
                (dbType, dbId) = dbAndId.split(":")
                if (dbType == "Ensembl"):
                    ensemblDbDict[dbId] = entrezId
             
    entrezFileHandler.close()
    return entrezDbDict, ensemblDbDict


def convert_events(aVCFDir, anEntrezDbFilename, aProtectedFilename, anOpenFilename, anIsDebug):

    if (aProtectedFilename != None):
        protectedFileHandler = get_write_fileHandler(aProtectedFilename)
    if (anOpenFilename != None):
        openFileHandler = get_write_fileHandler(anOpenFilename)
        
    # output the header lines
    headerVersion = "#version 2.4"
    headerList = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand"]
    headerList += ["Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status"]
    headerList += ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"]
    headerList += ["Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2"]
    headerList += ["Verification_Status", "Validation_Status", "Mutation_Status", "Sequencing_Phase", "Sequence_Source", "Validation_Method"]
    headerList += ["Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID"]
    headerList += ["RNA_Tumor_Sample_Barcode", "RNA_Tumor_Sample_UUID", "RNA_Tumor_Seq_Allele1", "RNA_Tumor_Seq_Allele2"]
    headerList += ["RNA_Normal_Sample_Barcode", "RNA_Normal_Sample_UUID", "RNA_Normal_Seq_Allele1", "RNA_Normal_Seq_Allele2"]
    headerList += ["Match_Norm_Ref_Count", "Match_Norm_Alt_Count", "Tumor_Ref_Count", "Tumor_Alt_Count"]
    headerList += ["RNA_Norm_Ref_Count", "RNA_Norm_Alt_Count", "RNA_Tumor_Ref_Count", "RNA_Tumor_Alt_Count"]    
    
    if (anOpenFilename != None):
        openFileHandler.write(headerVersion + "\n")
        openFileHandler.write("\t".join(headerList) + "\n")
        
    if (aProtectedFilename != None):
        protectedFileHandler.write(headerVersion + "\n")
        protectedFileHandler.write("\t".join(headerList) + "\n")
    else:
        print >> sys.stdout, headerVersion
        print >> sys.stdout, "\t".join(headerList)

    # load the entrez DB to convert HUGO names to entrez Ids
    (entrezDbDict, ensemblDbDict) = load_entrez_db(anEntrezDbFilename, anIsDebug)

    effectRegEx = re.compile("(\\w).*\\({1}")
    
    # hard-coded values for all lines
    center = "ucsc.edu"
    
    # for each vcf file
    for vcfFile in (glob.glob(aVCFDir + "*.vcf*")):
        # open the file
        vcfFileHandler = get_read_fileHandler(vcfFile)
        
        reference = None
        normalBarcode = ""
        normalUUID = ""
        tumorBarcode = ""
        tumorUUID = ""
        rnaNormalBarcode = ""
        rnaNormalUUID = ""
        rnaTumorBarcode = ""
        rnaTumorUUID = ""
        dnaTumorSequenceSource = ""
        rnaTumorSequenceSource = ""
                
        for line in vcfFileHandler:
            
            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            #if (anIsDebug):
            #    logging.debug("vcfLine: %s", line)
                
            # if it is an empty line, then just continue
            if (line.isspace()):
                continue;
            
            # we need to extract info from the reference tag in the header
            elif (line.startswith("##reference")):
                ##reference=<ID=GRCh37,Source=file:/inside/depot/fa/Homo_sapiens_assembly19.fasta>
                ##assembly=file:/inside/depot/fa/Homo_sapiens_assembly19.fasta
                
                if (anIsDebug):
                    logging.debug(line)
                
                if ("GRCh37-lite" in line):
                    reference = "GRCh37-lite"
                elif ("GRCh37" in line):
                    reference = "GRCh37"
                elif ("hg19" in line):
                    reference = "hg19"
                elif ("hg18" in line):
                    reference = "hg18"
                elif ("NCBI36" in line):
                    reference = "36.1"    
                elif ("NCBI37" in line):
                    reference = "37"    
                
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
                if (anIsDebug):
                    logging.debug(line)
                
                sampleLine = line[0:(len(line)-1)]
                sampleLine = sampleLine[len("##SAMPLE=<"):len(sampleLine)]
                sampleParamsList = sampleLine.split(",")
                sampleType = -1
                sampleBarcode = ""
                sampleUUID = ""
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
                    elif (key == "SequenceSource"):
                        sampleSource = value
                    elif (key == "Platform"):
                        samplePlatform = value
                    elif (key == "SampleName"):
                        tumorBarcode = value
                        
                if (isRNA):
                    rnaTumorSequenceSource = sampleSource
                else:
                    dnaTumorSequenceSource = sampleSource
                    
                if (anIsDebug):
                    logging.debug("SampleTCGABarcode=%s, SampleUUID=%s, SequenceSource=%s, Platform=%s", sampleBarcode, sampleUUID, sampleSource, samplePlatform)
                    
                # translate the sequencer into MAF acceptable format
                if (samplePlatform == "Illumina HiSeq 2000"):
                    platform = "Illumina HiSeq"
                elif (samplePlatform == "Illumina HiSeq 2500"):
                    platform = "Illumina HiSeq 2500"
                elif (samplePlatform == "Illumina Genome Analyzer IIx"):
                    platform = "Illumina GAIIx"
                elif (samplePlatform == "Illumina MiSeq"):
                    platform = "Illumina MiSeq"
                elif (samplePlatform == "Illumina"):
                    platform = "Illumina GAIIx"
                else:
                    logging.critical("Traceback:  Unexpected platform type %s", samplePlatform)
                    sys.exit(1)
                
                # the sample types from 0-9 are tumor and 10-19 are normal, 20 and above are control samples
                if (sampleType != -1):
                    if (sampleType < 10):
                        if (isRNA):
                            rnaTumorBarcode = sampleBarcode
                            rnaTumorUUID = sampleUUID
                        else:
                            tumorBarcode = sampleBarcode
                            tumorUUID = sampleUUID
                    elif (sampleType >= 10 and sampleType < 20):
                        if (isRNA):
                            rnaNormalBarcode = sampleBarcode
                            rnaNormalUUID = sampleUUID
                        else:
                            normalBarcode = sampleBarcode
                            normalUUID = sampleUUID
                    else:
                        logging.critical("Traceback:  Unexpected sample type %s", sampleType)
                        sys.exit(1)
            
                    if (anIsDebug):
                        logging.debug("normalBarcode=%s, normalUUID=%s, tumorBarcode=%s, tumorUUID=%s, rnaNormalBarcode=%s, rnaNormalUUID=%s, rnaTumorBarcode=%s, rnaTumorUUID=%s, SampleType=%s, dnaTumorSequenceSource=%s, rnaTumorSequenceSource=%s, Platform=%s", normalBarcode, normalUUID, tumorBarcode, tumorUUID, rnaNormalBarcode, rnaNormalUUID, rnaTumorBarcode, rnaTumorUUID, sampleType, dnaTumorSequenceSource, rnaTumorSequenceSource, platform)
            
            # if we find the vcfGenerator line, then create the dict of params
            elif ("vcfGenerator" in line):
                generatorLine = line[0:(len(line)-1)]
                generatorLine = generatorLine[16:len(generatorLine)]
                generatorParamsList = generatorLine.split(",")
                generatorParamsDict = {}
                
                # create a dictionary of existing params
                for param in generatorParamsList:
                    (key, value) = param.split("=")
                    value = value.rstrip(">")
                    value = value.lstrip("<")
                    generatorParamsDict[key] = value
                
                if (reference == None):
                    normalFasta = generatorParamsDict["dnaNormalFastaFilename"]
                    if ("GRCh37-lite" in normalFasta):
                        reference = "GRCh37-lite"
                    elif ("GRCh37" in normalFasta):
                        reference = "GRCh37"
                    elif ("hg19" in normalFasta):
                        reference = "hg19"
                    elif ("19" in normalFasta):
                        reference = "hg19"
                    elif ("hg18" in normalFasta):
                        reference = "hg18"
                    elif ("NCBI36" in normalFasta):
                        reference = "36.1"    
                    elif ("NCBI37" in normalFasta):
                        reference = "37"
                        
                    if (anIsDebug):
                        logging.debug("normalFasta=%s, reference=%s", normalFasta, reference)
                
            # if we find the column headers
            elif ("#CHROM" in line):
                columnsLine = line.lstrip("#")
                columnsList = columnsLine.split("\t")
                columnsList = columnsList[9:len(columnsList)]
                
            # these are other header lines that we don't need, so just continue
            elif (line.startswith("#")):
                continue;
            
            # now we are to the data
            else:
                
                # split the line on the tab
                splitLine = line.split("\t")
                
                # all calls start as open access and get eliminated below
                openAccess = True
                
                # get the fields to yield
                #columnHeaders = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                chrom = splitLine[0]
                stopCoordinate = int(splitLine[1]) 
                idList = splitLine[2].split(";")
                # they only allow one, so let's take the first one
                idList = [idList[0]]
            
                # if the list is empty, clear it
                if (len(idList) == 1 and "." in idList):
                    idList = ["novel"]
                    
                refList = splitLine[3].split(",")
                altList = splitLine[4].split(",")
                #score = float(splitLine[5])
                #filterSet = set(splitLine[6].split(";"))
                
                # parse the info column and create a dict
                infoList = splitLine[7].split(";")
                infoDict = collections.defaultdict(list)
                for info in infoList:
                    keyValueList = info.split("=")
                    # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
                    if (len(keyValueList) == 1):
                        infoDict[keyValueList[0]] = ["True"]
                    else:
                        # the value can be a comma separated list
                        infoDict[keyValueList[0]] = keyValueList[1].split(",")  
                        
                # Anything in dbSNP that is not in COSMIC is considered protected access
                if ("DB" in infoDict and "COSMIC" not in infoDict):
                    openAccess = False
                        
                # find out if the mutation was called from the DNA or RNA
                originList = infoDict["ORIGIN"]
                sequenceSource = ""
                if "DNA" in originList:
                    sequenceSource = dnaTumorSequenceSource
                elif "RNA" in originList:
                    sequenceSource = rnaTumorSequenceSource
                
                effectList = infoDict["EFF"]
                ignoreEffectsList = ["UPSTREAM", "DOWNSTREAM"]
                
                geneDict = collections.defaultdict(int)
                geneFamilySet = set()
                variantEffectSet = set()
                ensembleId = ""
                for rawEffect in effectList:
                    rawEffect = rawEffect.rstrip(")")
                    iterator = effectRegEx.finditer(rawEffect)
                        
                    # for each match object in the iterator
                    for match in iterator:
                        effect = match.group()
                        rawEffect = rawEffect.replace(effect, "")
                        effect = effect.rstrip("(")
                                
                    if (effect in ignoreEffectsList):
                        continue
                    
                    effectParts = rawEffect.split("|")
                    #effectImpact = effectParts[0]
                    functionalClass = effectParts[1]
                    #codonChange = effectParts[2]
                    #aaChange = effectParts[3]
                    #aaLength = effectParts[4]
                    geneName = effectParts[5]
                    transcriptBiotype = effectParts[6]
                    #geneCoding = effectParts[7]
                    ensembleId = effectParts[8]
                    #exonNumber = effectParts[9]
                    #genotypeNumber = effectParts[10]
                    
                    # These have MISSENSE as their functional impact, so they will be processed below
                    # START_GAINED, START_LOST, NON_SYNONYMOUS_START = MISSENSE 
                    # These have SILENT as their functional impact, so they will be processed below
                    # SYNONYMOUS_START, SYNONYMOUS_STOP = SILENT
                    # No place to really put these, but they are generally caught by the other info
                    # CDS, GENE, TRANSCRIPT, EXON
                    
                    if (effect == "INTERGENIC" or effect == "INTERGENIC_CONSERVED"):
                        variantEffectSet.add("IGR")
                    if (effect == "INTRON" or effect == "INTRON_CONSERVED"):
                        variantEffectSet.add("Intron")
                    elif (effect == "UTR_5_PRIME" or effect == "UTR_5_DELETED"):
                        variantEffectSet.add("5'UTR")
                    elif (effect == "UTR_3_PRIME" or effect == "UTR_3_DELETED"):
                        variantEffectSet.add("3'UTR")
                    elif (effect == "SPLICE_SITE_ACCEPTOR" or effect == "SPLICE_SITE_DONOR"):
                        variantEffectSet.add("Splice_Site")
                    elif (effect == "SPLICE_SITE_ACCEPTOR" or effect == "SPLICE_SITE_DONOR" or effect == "SPLICE_SITE_BRANCH"):
                        variantEffectSet.add("Splice_Site")
                    elif (effect == "FRAME_SHIFT"):
                        variantEffectSet.add("Frame_Shift_Del")
                        variantEffectSet.add("Frame_Shift_Ins")
                    elif (effect == "STOP_LOST"):
                        variantEffectSet.add("Nonstop_Mutation")    
                    elif (functionalClass == "MISSENSE"):
                        variantEffectSet.add("Missense_Mutation")
                    elif (functionalClass == "SILENT"):
                        variantEffectSet.add("Silent")
                    elif (functionalClass == "NONSENSE"):
                        variantEffectSet.add("Nonsense_Mutation")
                    elif (transcriptBiotype == "retained_intron"):
                        variantEffectSet.add("Intron")
                    
                    geneDict[geneName] += 1
                    geneFamilySet.add(transcriptBiotype)
                
                # can only specify one variant effect, so pick the order
                finalVariantEffect = None
                if ("Missense_Mutation" in variantEffectSet):
                    finalVariantEffect = "Missense_Mutation"
                elif ("Nonsense_Mutation" in variantEffectSet):
                    finalVariantEffect = "Nonsense_Mutation"
                elif ("Nonstop_Mutation" in variantEffectSet):
                    finalVariantEffect = "Nonstop_Mutation"
                elif ("Splice_Site" in variantEffectSet):
                    finalVariantEffect = "Splice_Site"
                elif ("Frame_Shift_Del" in variantEffectSet):
                    finalVariantEffect = "Frame_Shift_Del"
                elif ("Silent" in variantEffectSet):
                    finalVariantEffect = "Silent"
                elif ("5'UTR" in variantEffectSet):
                    finalVariantEffect = "5'UTR"
                    openAccess = False
                elif ("3'UTR" in variantEffectSet):
                    finalVariantEffect = "3'UTR"
                    openAccess = False
                elif ("Intron" in variantEffectSet):
                    finalVariantEffect = "Intron"
                    openAccess = False
                elif ("IGR" in variantEffectSet):
                    finalVariantEffect = "IGR"
                    openAccess = False
                else:
                    # if there is no snpEff output, then it couldn't find any annotation for this position, so let's make it IGR
                    finalVariantEffect = "IGR"
                    openAccess = False
                
                geneMaxCount = 0
                geneMaxName = "Unknown"
                for (gene, count) in geneDict.iteritems():
                    if (count > geneMaxCount):
                        geneMaxCount = count
                        geneMaxName = gene
                        
                entrezId = get_entrez_id(entrezDbDict, ensemblDbDict, geneMaxName, ensembleId, anIsDebug)
                
                '''
                # get the first non-zero entrez id for the hugo gene symbols
                entrezId = "0"
                for gene in geneSet:
                    tmpId = get_entrez_id(entrezDbDict, ensemblDbDict, gene, ensembleId, anIsDebug)
                    if (tmpId != "0"):
                        entrezId = tmpId
                        break 
                   
                # add the hugo name and entrez id to the output
                if (len(geneSet) > 0):
                    mafOutputList = [",".join(geneSet), entrezId, center, reference, chrom, str(stopCoordinate), str(stopCoordinate), "+", finalVariantEffect, "SNP"]
                else:
                    mafOutputList = ["Unknown", "0", center, reference, chrom, str(stopCoordinate), str(stopCoordinate), "+", finalVariantEffect, "SNP"]
                '''
                
                mafOutputList = [geneMaxName, entrezId, center, reference, chrom, str(stopCoordinate), str(stopCoordinate), "+", finalVariantEffect, "SNP"]
                
                # get the event format list
                formatList = splitLine[8].split(":")
                
                # initialize the optional columns to none
                dnaNormalList = None
                rnaNormalList = None
                dnaTumorList = None
                rnaTumorList = None
    
                # if we have a 9th column, figure out which dataset it is
                if (len(splitLine) > 9):
                    if (columnsList[0] == "DNA_NORMAL"):
                        dnaNormalList = splitLine[9].split(":")
                    elif (columnsList[0] == "RNA_NORMAL"):
                        rnaNormalList = splitLine[9].split(":")
                    elif (columnsList[0] == "DNA_TUMOR"):
                        dnaTumorList = splitLine[9].split(":")
                    elif (columnsList[0] == "RNA_TUMOR"):
                        rnaTumorList = splitLine[9].split(":")
                # if we have a 10th column, figure out which dataset it is
                if (len(splitLine) > 10):
                    if (columnsList[1] == "RNA_NORMAL"):
                        rnaNormalList = splitLine[10].split(":")
                    elif (columnsList[1] == "DNA_TUMOR"):
                        dnaTumorList = splitLine[10].split(":")
                    elif (columnsList[1] == "RNA_TUMOR"):
                        rnaTumorList = splitLine[10].split(":") 
                # if we have a 11th column, figure out which dataset it is
                if (len(splitLine) > 11):
                    if (columnsList[2] == "DNA_TUMOR"):
                        dnaTumorList = splitLine[11].split(":")
                    elif (columnsList[2] == "RNA_TUMOR"):
                        rnaTumorList = splitLine[11].split(":")
                # if we have a 12th column, figure out which dataset it is
                if (len(splitLine) > 12):
                    if (columnsList[3] == "RNA_TUMOR"):
                        rnaTumorList = splitLine[12].split(":")
                
                haveDnaNormData = True
                haveRnaNormData = True
                haveDnaTumData = True
                haveRnaTumData = True
                
                # if there is no data, then set the flag
                if (dnaNormalList == None or dnaNormalList[0] == "." or dnaNormalList[0] == ".:.:.:.:.:.:.:.:.:.:."):
                    haveDnaNormData = False
                # if there is no data, then set the flag
                if (rnaNormalList == None or rnaNormalList[0] == "." or rnaNormalList[0] == ".:.:.:.:.:.:.:.:.:.:."):
                    haveRnaNormData = False
                # if there is no data, then set the flag
                if (dnaTumorList == None or dnaTumorList[0] == "." or dnaTumorList[0] == ".:.:.:.:.:.:.:.:.:.:."):
                    haveDnaTumData = False
                # if there is no data, then set the flag
                if (rnaTumorList == None or rnaTumorList[0] == "." or rnaTumorList[0] == ".:.:.:.:.:.:.:.:.:.:."):
                    haveRnaTumData = False
                             
                # parse the dna and rna columns and create dicts for each
                dnaNormalDict = collections.defaultdict(list)
                rnaNormalDict = collections.defaultdict(list)
                dnaTumorDict = collections.defaultdict(list)
                rnaTumorDict = collections.defaultdict(list)
                
                try:
                    index = 0
                    for formatItem in formatList:
                        if (formatItem == "GT"):
                            sep = "/"
                        else:
                            sep = ","
                            
                        if (haveDnaNormData):
                            dnaNormalItem = dnaNormalList[index]
                            dnaNormalDict[formatItem] = dnaNormalItem.split(sep)
                        if (haveRnaNormData):
                            rnaNormalItem = rnaNormalList[index]
                            rnaNormalDict[formatItem] = rnaNormalItem.split(sep)
                        if (haveDnaTumData):
                            dnaTumorItem = dnaTumorList[index]
                            dnaTumorDict[formatItem] = dnaTumorItem.split(sep)
                        if (haveRnaTumData):
                            rnaTumorItem = rnaTumorList[index]
                            rnaTumorDict[formatItem] = rnaTumorItem.split(sep)
                        index += 1
                except:
                    print "Error", line
                    
                refPlusAltList = refList + altList
                
                dnaNormalRefDepth = ""
                dnaNormalAltDepth= ""
                dnaTumorRefDepth = ""
                dnaTumorAltDepth= ""
                rnaNormalRefDepth = ""
                rnaNormalAltDepth= ""
                rnaTumorRefDepth = ""
                rnaTumorAltDepth= ""
                normalAllele1 = ""
                normalAllele2 = ""
                tumorAllele1 = ""
                tumorAllele2 = ""
                rnaNormalAllele1 = ""
                rnaNormalAllele2 = ""
                rnaTumorAllele1 = ""
                rnaTumorAllele2 = ""
                for (modType, modChange) in izip(infoDict["MT"], infoDict["MC"]):    
                    
                    if (modType == "SOM" or modType == "TUM_EDIT"):
                        # get the source and target alleles
                        (source, target) = modChange.split(">")
                        sourceIndex = refPlusAltList.index(source)     
                        targetIndex = refPlusAltList.index(target)
                        
                        referenceAllele = source
                        
                        if (haveDnaNormData):
                            dnaNormalRefDepth = int(dnaNormalDict["AD"][sourceIndex])
                            dnaNormalAltDepth = int(dnaNormalDict["AD"][targetIndex])
                            normalAllele1 = source
                            normalAllele2 = source
                            '''
                            # chrom Y only has one genotype
                            if (len(dnaNormalDict["GT"]) == 2):
                                normalAllele1 = refPlusAltList[int(dnaNormalDict["GT"][0])]
                                normalAllele2 = refPlusAltList[int(dnaNormalDict["GT"][1])]
                            else:
                                #normalAllele1 = refPlusAltList[int(dnaNormalDict["GT"][0])]
                                #normalAllele2 = refPlusAltList[int(dnaNormalDict["GT"][0])]
                                normalAllele1 = source
                                normalAllele2 = source
                            '''    
                        if (haveDnaTumData):
                            dnaTumorRefDepth = int(dnaTumorDict["AD"][sourceIndex])
                            dnaTumorAltDepth = int(dnaTumorDict["AD"][targetIndex])
                            tumorAllele1 = source
                            tumorAllele2 = target    
                            '''
                            if (len(dnaTumorDict["GT"]) == 2):
                                tumorAllele1 = refPlusAltList[int(dnaTumorDict["GT"][0])]
                                tumorAllele2 = refPlusAltList[int(dnaTumorDict["GT"][1])]
                            else:
                                #tumorAllele1 = refPlusAltList[int(dnaTumorDict["GT"][0])]
                                #tumorAllele2 = refPlusAltList[int(dnaTumorDict["GT"][0])]
                                tumorAllele1 = source
                                tumorAllele2 = target
                            '''
                        if (haveRnaNormData):
                            rnaNormalRefDepth = int(rnaNormalDict["AD"][sourceIndex])
                            rnaNormalAltDepth = int(rnaNormalDict["AD"][targetIndex])
                            rnaNormalAllele1 = source
                            rnaNormalAllele2 = source
                            '''
                            if (len(rnaNormalDict["GT"]) == 2):
                                rnaNormalAllele1 = refPlusAltList[int(rnaNormalDict["GT"][0])]
                                rnaNormalAllele2 = refPlusAltList[int(rnaNormalDict["GT"][1])]
                            else:
                                #rnaNormalAllele1 = refPlusAltList[int(rnaNormalDict["GT"][0])]
                                #rnaNormalAllele2 = refPlusAltList[int(rnaNormalDict["GT"][0])]
                                rnaNormalAllele1 = source
                                rnaNormalAllele2 = source
                            '''
                        if (haveRnaTumData):
                            rnaTumorRefDepth = int(rnaTumorDict["AD"][sourceIndex])
                            rnaTumorAltDepth = int(rnaTumorDict["AD"][targetIndex])
                            rnaTumorAllele1 = source
                            if (rnaTumorAltDepth > 3):
                                rnaTumorAllele2 = target
                            else:
                                rnaTumorAllele2 = source
                            '''
                            if (len(dnaNormalDict["GT"]) == 2):
                                rnaTumorAllele1 = refPlusAltList[int(rnaTumorDict["GT"][0])]
                                rnaTumorAllele2 = refPlusAltList[int(rnaTumorDict["GT"][1])]
                            else:
                                #rnaTumorAllele1 = refPlusAltList[int(rnaTumorDict["GT"][0])]
                                #rnaTumorAllele2 = refPlusAltList[int(rnaTumorDict["GT"][0])]
                                rnaTumorAllele1 = source
                                if (rnaTumorAltDepth > 3):
                                    rnaTumorAllele2 = target
                                else:
                                    rnaTumorAllele2 = source
                            '''
                        #print line
                        #print normalAllele1, normalAllele2, tumorAllele1, tumorAllele2, rnaNormalAllele1, rnaNormalAllele2, rnaTumorAllele1, rnaTumorAllele2
               
                #mafOutputList += [referenceAllele, tumorAllele1, tumorAllele2, ";".join(idList), ";".join(["bySubmitter"] * len(idList))]
                mafOutputList += [referenceAllele, tumorAllele1, tumorAllele2, ";".join(idList), ""]     
                # tumorValAllele1, tumorValAllele2, normalValAllele1, normalValAllele2, VerificationStatus, ValidationStatus, MutationStatus
                mafOutputList += [tumorBarcode, normalBarcode, normalAllele1, normalAllele2, "", "", "", "", "Unknown", "Untested", "Somatic"]
                # sequencingPhase, sequenceSource, validationMethod, score, bamFile, sequencer, tumorUUID, normalUUID
                mafOutputList += ["", sequenceSource, "none", "", "", platform, tumorUUID, normalUUID]
                #["RNA_Tumor_Sample_Barcode", "RNA_Tumor_Sample_UUID", "RNA_Tumor_Seq_Allele1", "RNA_Tumor_Seq_Allele2"]
                mafOutputList += [rnaTumorBarcode, rnaTumorUUID, rnaTumorAllele1, rnaTumorAllele2]
                #["RNA_Norm_Sample_Barcode", "RNA_Norm_Sample_UUID", "RNA_Norm_Seq_Allele1", "RNA_Norm_Seq_Allele2"]
                mafOutputList += [rnaNormalBarcode, rnaNormalUUID, rnaNormalAllele1, rnaNormalAllele2]
                #["Match_Norm_Ref_Count", "Match_Norm_Alt_Count", "Tumor_Ref_Count", "Tumor_Alt_Count"]
                mafOutputList += [str(dnaNormalRefDepth), str(dnaNormalAltDepth), str(dnaTumorRefDepth), str(dnaTumorAltDepth)]
                #["RNA_Norm_Ref_Count", "RNA_Norm_Alt_Count", "RNA_Tumor_Ref_Count", "RNA_Tumor_Alt_Count"
                mafOutputList += [str(rnaNormalRefDepth), str(rnaNormalAltDepth), str(rnaTumorRefDepth), str(rnaTumorAltDepth)]
                
                if (aProtectedFilename != None):
                    protectedFileHandler.write("\t".join(mafOutputList) + "\n")
                else:
                    print >> sys.stdout, "\t".join(mafOutputList)
                    
                if (openAccess and anOpenFilename != None):
                    openFileHandler.write("\t".join(mafOutputList) + "\n")
            
        # close the file
        vcfFileHandler.close()
        
    if (anOpenFilename != None):
        openFileHandler.close()
    if (aProtectedFilename != None):
        protectedFileHandler.close()
        
    return


def main():
    
    # create the usage statement
    usage = "usage: python %prog vcfDir entrezDbFilename [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-p", "--protectedFilename", dest="protectedFilename", metavar="PROTECTED_FILE", help="the name of the protected access MAF output file, STDOUT by default")
    i_cmdLineParser.add_option("-o", "--openAccessFilename", dest="openAccessFilename", metavar="OPEN_ACCESS_FILE", help="the name of the open access MAF output file")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(2,11,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_vcfDir = str(i_cmdLineArgs[0])
    i_entrezDbFilename = str(i_cmdLineArgs[1])
    
    # get the optional params with default values   
    i_logLevel = i_cmdLineOptions.logLevel
    
    # try to get any optional parameters with no defaults
    i_protectedFilename = None
    i_openFilename = None    
    i_logFilename = None
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
    if (i_cmdLineOptions.protectedFilename != None):
        i_protectedFilename = str(i_cmdLineOptions.protectedFilename)
    if (i_cmdLineOptions.openAccessFilename != None):
        i_openFilename = str(i_cmdLineOptions.openAccessFilename)
         
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
        logging.debug("input=%s", i_vcfDir)
        logging.debug("protected=%s", i_protectedFilename)
        logging.debug("open=%s", i_openFilename)
        logging.debug("i_logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        
    # check for any errors
    writeFilenameList = []
    if (i_protectedFilename != None):
        writeFilenameList += [i_protectedFilename]
    if (i_openFilename != None):
        writeFilenameList += [i_openFilename]
    if (i_logFilename != None):
        writeFilenameList += [i_logFilename]
        
    readFilenameList = [i_entrezDbFilename]        
    if (not rnaEditingUtil.check_for_argv_errors([i_vcfDir], readFilenameList, writeFilenameList)):
        sys.exit(1)           
    
    convert_events(i_vcfDir, i_entrezDbFilename, i_protectedFilename, i_openFilename, i_debug)
       
    return

main()
sys.exit(0)
