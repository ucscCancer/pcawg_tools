#!/usr/bin/env python

import argparse, os, shutil, subprocess, tempfile, time
from multiprocessing import Pool

def execute(cmd, output=None):
    import sys, shlex
    # function to execute a cmd and report if an error occur
    print(cmd)
    try:
        process = subprocess.Popen(args=shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = process.communicate()
    except Exception, e: # une erreur de ma commande : stderr
        sys.stderr.write("problem doing : %s\n%s\n" %(cmd, e))
        return 1
    if output:
        output = open(output, 'w')
        output.write(stdout)
        output.close()
    if stderr != '': # une erreur interne au programme : stdout (sinon, souvent des warning arrete les programmes)
        sys.stdout.write("warning or error while doing : %s\n-----\n%s-----\n\n" %(cmd, stderr))
        return 1
    return 0


def indexBam(workdir, prefix, inputBamFile, inputBamFileIndex=None):
    inputBamLink = os.path.join(os.path.abspath(workdir), prefix + ".bam" )
    os.symlink(inputBamFile, inputBamLink)
    if inputBamFileIndex is None:
        cmd = "samtools index %s" %(inputBamLink)
        execute(cmd)
    else:
        os.symlink(inputBamFileIndex, inputBamLink + ".bai")
    return inputBamLink


def indexFasta(workdir, inputFastaFile, inputFastaFileIndex=None, prefix="dna"):
    inputFastaLink = os.path.join(os.path.abspath(workdir), prefix + "_reference.fa" )
    os.symlink(inputFastaFile, inputFastaLink)
    if inputFastaFileIndex is None:
        cmd = "samtools faidx %s" %(inputFastaLink)
        execute(cmd)
    else:
        os.symlink(inputFastaFileIndex, inputFastaLink + ".fai")
    return inputFastaLink

def bamChrScan(bamfile):
    samtools = which("samtools")
    cmd = [samtools, "idxstats", bamfile]
    process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    found = False
    for line in stdout.split("\n"):
        tmp = line.split("\t")
        if len(tmp) == 4 and tmp[0].startswith("chr"):
            found = True
    return found

def radia(chrom, args,
          dnaNormalFilename=None, rnaNormalFilename=None, dnaTumorFilename=None, rnaTumorFilename=None,
          dnaNormalFastaFilename=None, rnaNormalFastaFilename=None, dnaTumorFastaFilename=None, rnaTumorFastaFilename=None):

    # python radia.pyc id chrom [Options]

    # python radia.pyc TCGA-02-0047 1
    # -n TCGA-02-0047-10A-01D-1490-08.bam
    # -t TCGA-02-0047-01A-01D-1490-08.bam
    # -r TCGA-02-0047-01A*.bam
    # --rnaTumorUseChr
    # --rnaTumorFasta=GRCh37-lite_w_chr_prefix.fa
    # -f Homo_sapiens_assembly19.fasta
    # -o TCGA-02-0047-10A-01D-1490-08_TCGA-02-0047-01A-01D-1490-08_chr1.vcf.gz
    # -i GRCh37
    # -m Homo_sapiens_assembly19.fasta
    # -d CGHub
    # -q Illumina
    # --disease GBM
    # --log=INFO

    # quadruplets
    if (rnaNormalFilename != None and rnaTumorFilename != None):
        cmd = "python %s/radia.pyc %s %s -n %s -x % -t %s -r %s --dnaNormalFasta %s --rnaNormalFasta %s --dnaTumorFasta %s --rnaTumorFasta %s " %(
                args.scriptsDir,
                args.patientId, chrom,
                dnaNormalFilename, rnaNormalFilename, dnaTumorFilename, rnaTumorFilename,
                dnaNormalFastaFilename, rnaNormalFastaFilename, dnaTumorFastaFilename, rnaTumorFastaFilename)
    # triplets
    elif (rnaTumorFilename != None):
        cmd = "python %s/radia.pyc %s %s -n %s -t %s -r %s --dnaNormalFasta %s --dnaTumorFasta %s --rnaTumorFasta %s " %(
                args.scriptsDir,
                args.patientId, chrom,
                dnaNormalFilename, dnaTumorFilename, rnaTumorFilename,
                dnaNormalFastaFilename, dnaTumorFastaFilename, rnaTumorFastaFilename)
    # pairs
    else:
        cmd = "python %s/radia.pyc %s %s -n %s -t %s --dnaNormalFasta %s --dnaTumorFasta %s " % (
                args.scriptsDir,
                args.patientId, chrom,
                dnaNormalFilename, dnaTumorFilename,
                dnaNormalFastaFilename, dnaTumorFastaFilename
        )

    if dnaNormalFilename is not None and bamChrScan(dnaNormalFilename):
        cmd += ' --dnaNormalUseChr '
    if rnaNormalFilename is not None and bamChrScan(rnaNormalFilename):
        cmd += ' --rnaNormalUseChr '
    if dnaTumorFilename is not None and bamChrScan(dnaTumorFilename):
        cmd += ' --dnaTumorUseChr '
    if rnaTumorFilename is not None and bamChrScan(rnaTumorFilename):
        cmd += ' --rnaTumorUseChr '
    if args.gzip:
        cmd += ' --gzip '
    return cmd, args.outputFilename


def radiaFilter(args, chrom, rawFile, outputDir, scriptsDir, blatFastaFilename):

    # python filterRadia.pyc id chrom inputFile outputDir scriptsDir [Options]

    # python filterRadia.pyc TCGA-02-0047-10A-01D-1490-08_TCGA-02-0047-01A-01D-1490-08 1
    # TCGA-02-0047-10A-01D-1490-08_TCGA-02-0047-01A-01D-1490-08_chr1.vcf.gz
    # /radia/finalChromVCFs/
    # /rnaEditing/scripts/
    # --blacklistDir /rnaEditing/data/hg19/blacklists/1000Genomes/phase1/
    # --dbSnpDir /rnaEditing/data/hg19/snp135/
    # --retroGenesDir /rnaEditing/data/hg19/retroGenes/
    # --pseudoGenesDir /rnaEditing/data/hg19/pseudoGenes/
    # --cosmicDir /rnaEditing/data/hg19/cosmic/
    # --targetDir /rnaEditing/data/hg19/broadTargets/
    # --snpEffDir /snpEff/
    # --rnaGeneBlckFile /rnaEditing/data/rnaGeneBlacklist.tab
    # --rnaGeneFamilyBlckFile /rnaEditing/data/rnaGeneFamilyBlacklist.tab
    # --blatFastaFilename hg19.fasta
    # --canonical
    # --log=INFO
    # --gzip

    cmd = "python %s/filterRadia.pyc %s %s %s %s %s --blatFastaFilename %s " % (
        args.scriptsDir,
        args.patientId, chrom, rawFile,
        outputDir, scriptsDir, blatFastaFilename)

    if args.canonical:
        cmd += ' --canonical '
    if args.noBlacklist:
        cmd += ' --noBlacklist '
    if args.noTargets:
        cmd += ' --noTargets '
    if args.noDbSnp:
        cmd += ' --noDbSnp '
    if args.noRetroGenes:
        cmd += ' --noRetroGenes '
    if args.noPseudoGenes:
        cmd += ' --noPseudoGenes '
    if args.noCosmic:
        cmd += ' --noCosmic '
    if args.noBlat:
        cmd += ' --noBlat '
    if args.noPositionalBias:
        cmd += ' --noPositionalBias '
    if args.noRnaBlacklist:
        cmd += ' --noRnaBlacklist '
    if args.noSnpEff:
        cmd += ' --noSnpEff '
    if args.dnaOnly:
        cmd += ' --dnaOnly '
    if args.rnaOnly:
        cmd += ' --rnaOnly '
    if args.gzip:
        cmd += ' --gzip '

    filteredOut = outputDir + args.patientId + "_" + chrom + ".vcf"
    return cmd, filteredOut


def move(avant, apres):
    if os.path.exists(avant):
        execute("mv %s %s" %(avant, apres))

def which(cmd):
    cmd = ["which",cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: return None
    return res

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


def __main__():
    time.sleep(1) #small hack, sometimes it seems like docker file systems are avalible instantly
    parser = argparse.ArgumentParser(description="RNA and DNA Integrated Analysis (RADIA)")


    #############################
    #    RADIA params    #
    #############################
    parser.add_argument("-o", "--outputFilename", dest="outputFilename", required=True, metavar="OUTPUT_FILE", help="the name of the output file")
    parser.add_argument("--patientId", dest="patientId", required=True, metavar="PATIENT_ID", help="a unique patient Id that will be used to name the output file")

    parser.add_argument("-b", "--batchSize", type=int, dest="batchSize", default=int(250000000), metavar="BATCH_SIZE", help="the size of the samtool selections that are loaded into memory at one time, %default by default")
    parser.add_argument("-c", "--chromSizesFilename", dest="chromSizesFilename", metavar="CHROM_SIZES_FILE", help="the name of the file with the chromosome sizes")
    parser.add_argument("-f", "--fastaFilename", dest="fastaFilename", metavar="FASTA_FILE", help="the name of the fasta file that can be used on all .bams, see below for specifying individual fasta files for each .bam file")
    parser.add_argument("-p", "--useChrPrefix", action="store_true", default=False, dest="useChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for all .bams, see below for specifying the prefix for individual .bam files")
    parser.add_argument("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    parser.add_argument("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    parser.add_argument("-i", "--refId", dest="refId", metavar="REF_ID", help="the reference Id - used in the reference VCF meta tag")
    parser.add_argument("-u", "--refUrl", dest="refUrl", metavar="REF_URL", help="the URL for the reference - used in the reference VCF meta tag")
    parser.add_argument("-m", "--refFilename", dest="refFilename", metavar="REF_FILE", help="the location of the reference - used in the reference VCF meta tag")
    parser.add_argument("-a", "--startCoordinate", type=int, default=int(1), dest="startCoordinate", metavar="START_COORDINATE", help="the start coordinate for testing small regions, %default by default")
    parser.add_argument("-z", "--stopCoordinate", type=int, default=int(0), dest="stopCoordinate", metavar="STOP_COORDINATE", help="the stop coordinate for testing small regions, %default by default")
    parser.add_argument("-d", "--dataSource", dest="dataSource", metavar="DATA_SOURCE", help="the source of the data - used in the sample VCF meta tag")
    parser.add_argument("-q", "--sequencingPlatform", dest="sequencingPlatform", metavar="SEQ_PLATFORM", help="the sequencing platform - used in the sample VCF meta tag")
    parser.add_argument("-s", "--statsDir", dest="statsDir", metavar="STATS_DIR", help="a stats directory where some basic stats can be output")
    parser.add_argument("--disease", dest="disease", metavar="DISEASE", help="a disease abbreviation (i.e. BRCA) for the header")

    parser.add_argument("--genotypeMinDepth", type=int, default=int(2), dest="genotypeMinDepth", metavar="GT_MIN_DP", help="the minimum number of bases required for the genotype, %default by default")
    parser.add_argument("--genotypeMinPct", type=float, default=float(.10), dest="genotypeMinPct", metavar="GT_MIN_PCT", help="the minimum percentage of reads required for the genotype, %default by default")

    # params for normal DNA
    parser.add_argument("-n", "--dnaNormalFilename", dest="dnaNormalFilename", metavar="DNA_NORMAL_FILE", help="the name of the normal DNA .bam file")
    parser.add_argument("--dnaNormalBaiFilename", dest="dnaNormalBaiFilename", metavar="DNA_NORMAL_BAI_FILE", help="the name of the normal DNA .bai file")
    parser.add_argument("--dnaNormalMinTotalBases", type=int, default=int(4), dest="dnaNormalMinTotalNumBases", metavar="DNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal DNA reads covering a position, %default by default")
    parser.add_argument("--dnaNormalMinAltBases", type=int, default=int(2), dest="dnaNormalMinAltNumBases", metavar="DNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal DNA reads supporting a variant at a position, %default by default")
    parser.add_argument("--dnaNormalBaseQual", type=int, default=int(10), dest="dnaNormalMinBaseQuality", metavar="DNA_NOR_BASE_QUAL", help="the minimum normal DNA base quality, %default by default")
    parser.add_argument("--dnaNormalMapQual", type=int, default=int(10), dest="dnaNormalMinMappingQuality", metavar="DNA_NOR_MAP_QUAL", help="the minimum normal DNA mapping quality, %default by default")
    parser.add_argument("--dnaNormalUseChr", action="store_true", default=False, dest="dnaNormalUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the normal DNA .bam file")
    parser.add_argument("--dnaNormalFasta", dest="dnaNormalFastaFilename", metavar="DNA_NOR_FASTA_FILE", help="the name of the fasta file for the normal DNA .bam file")
    parser.add_argument("--dnaNormalMitochon", default = "M", dest="dnaNormalMitochon", metavar="DNA_NOR_MITOCHON", help="the short name for the mitochondrial DNA (e.g 'M' or 'MT'), %default by default")
    parser.add_argument("--dnaNormalDescription", default = "Normal DNA Sample", dest="dnaNormalDesc", metavar="DNA_NOR_DESC", help="the description for the sample in the VCF header, %default by default")

    # params for normal RNA
    parser.add_argument("-x", "--rnaNormalFilename", dest="rnaNormalFilename", metavar="RNA_NORMAL_FILE", help="the name of the normal RNA-Seq .bam file")
    parser.add_argument("--rnaNormalBaiFilename", dest="rnaNormalBaiFilename", metavar="RNA_NORMAL_BAI_FILE", help="the name of the normal RNA .bai file")
    parser.add_argument("--rnaNormalMinTotalBases", type=int, default=int(4), dest="rnaNormalMinTotalNumBases", metavar="RNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal RNA-Seq reads covering a position, %default by default")
    parser.add_argument("--rnaNormalMinAltBases", type=int, default=int(2), dest="rnaNormalMinAltNumBases", metavar="RNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    parser.add_argument("--rnaNormalBaseQual", type=int, default=int(10), dest="rnaNormalMinBaseQuality", metavar="RNA_NOR_BASE_QUAL", help="the minimum normal RNA-Seq base quality, %default by default")
    parser.add_argument("--rnaNormalMapQual", type=int, default=int(10), dest="rnaNormalMinMappingQuality", metavar="RNA_NOR_MAP_QUAL", help="the minimum normal RNA-Seq mapping quality, %default by default")
    parser.add_argument("--rnaNormalUseChr", action="store_true", default=False, dest="rnaNormalUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the normal RNA .bam file")
    parser.add_argument("--rnaNormalFasta", dest="rnaNormalFastaFilename", metavar="RNA_NOR_FASTA_FILE", help="the name of the fasta file for the normal RNA .bam file")
    parser.add_argument("--rnaNormalMitochon", default = "M", dest="rnaNormalMitochon", metavar="RNA_NOR_MITOCHON", help="the short name for the mitochondrial RNA (e.g 'M' or 'MT'), %default by default")
    parser.add_argument("--rnaNormalDescription", default = "Normal RNA Sample", dest="rnaNormalDesc", metavar="RNA_NOR_DESC", help="the description for the sample in the VCF header, %default by default")

    # params for tumor DNA
    parser.add_argument("-t", "--dnaTumorFilename", dest="dnaTumorFilename", metavar="DNA_TUMOR_FILE", help="the name of the tumor DNA .bam file")
    parser.add_argument("--dnaTumorBaiFilename", dest="dnaTumorBaiFilename", metavar="DNA_TUMOR_BAI_FILE", help="the name of the tumor DNA .bai file")
    parser.add_argument("--dnaTumorMinTotalBases", type=int, default=int(4), dest="dnaTumorMinTotalNumBases", metavar="DNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor DNA reads covering a position, %default by default")
    parser.add_argument("--dnaTumorMinAltBases", type=int, default=int(2), dest="dnaTumorMinAltNumBases", metavar="DNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor DNA reads supporting a variant at a position, %default by default")
    parser.add_argument("--dnaTumorBaseQual", type=int, default=int(10), dest="dnaTumorMinBaseQuality", metavar="DNA_TUM_BASE_QUAL", help="the minimum tumor DNA base quality, %default by default")
    parser.add_argument("--dnaTumorMapQual", type=int, default=int(10), dest="dnaTumorMinMappingQuality", metavar="DNA_TUM_MAP_QUAL", help="the minimum tumor DNA mapping quality, %default by default")
    parser.add_argument("--dnaTumorUseChr", action="store_true", default=False, dest="dnaTumorUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the tumor DNA .bam file")
    parser.add_argument("--dnaTumorFasta", dest="dnaTumorFastaFilename", metavar="DNA_TUM_FASTA_FILE", help="the name of the fasta file for the tumor DNA .bam file")
    parser.add_argument("--dnaTumorMitochon", default = "M", dest="dnaTumorMitochon", metavar="DNA_TUM_MITOCHON", help="the short name for the mitochondrial DNA (e.g 'M' or 'MT'), %default by default")
    parser.add_argument("--dnaTumorDescription", default = "Tumor DNA Sample", dest="dnaTumorDesc", metavar="DNA_TUM_DESC", help="the description for the sample in the VCF header, %default by default")

    # params for tumor RNA
    parser.add_argument("-r", "--rnaTumorFilename", dest="rnaTumorFilename", metavar="RNA_TUMOR_FILE", help="the name of the tumor RNA-Seq .bam file")
    parser.add_argument("--rnaTumorBaiFilename", dest="rnaTumorBaiFilename", metavar="RNA_TUMOR_BAI_FILE", help="the name of the tumor RNA .bai file")
    parser.add_argument("--rnaTumorMinTotalBases", type=int, default=int(4), dest="rnaTumorMinTotalNumBases", metavar="RNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor RNA-Seq reads covering a position, %default by default")
    parser.add_argument("--rnaTumorMinAltBases", type=int, default=int(2), dest="rnaTumorMinAltNumBases", metavar="RNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    parser.add_argument("--rnaTumorBaseQual", type=int, default=int(10), dest="rnaTumorMinBaseQuality", metavar="RNA_TUM_BASE_QUAL", help="the minimum tumor RNA-Seq base quality, %default by default")
    parser.add_argument("--rnaTumorMapQual", type=int, default=int(10), dest="rnaTumorMinMappingQuality", metavar="RNA_TUM_MAP_QUAL", help="the minimum tumor RNA-Seq mapping quality, %default by default")
    parser.add_argument("--rnaTumorUseChr", action="store_true", default=False, dest="rnaTumorUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the tumor RNA .bam file")
    parser.add_argument("--rnaTumorFasta", dest="rnaTumorFastaFilename", metavar="RNA_TUM_FASTA_FILE", help="the name of the fasta file for the tumor RNA .bam file")
    parser.add_argument("--rnaTumorMitochon", default = "M", dest="rnaTumorMitochon", metavar="RNA_TUM_MITOCHON", help="the short name for the mitochondrial RNA (e.g 'M' or 'MT'), %default by default")
    parser.add_argument("--rnaTumorDescription", default = "Tumor RNA Sample", dest="rnaTumorDesc", metavar="RNA_TUM_DESC", help="the description for the sample in the VCF header, %default by default")


    #############################
    #    RADIA filter params    #
    #############################
    parser.add_argument("--outputDir", dest="outputDir", required=True, metavar="FILTER_OUT_DIR", help="the directory where temporary and final filtered output should be stored")
    parser.add_argument("--scriptsDir", dest="scriptsDir", required=True, metavar="SCRIPTS_DIR", help="the directory that contains the RADIA filter scripts")

    parser.add_argument("--blacklistDir", dest="blacklistDir", metavar="BLACKLIST_DIR", help="the path to the blacklist directory")
    parser.add_argument("--targetDir", dest="targetDir", metavar="TARGET_DIR", help="the path to the exon capture targets directory")
    parser.add_argument("--dbSnpDir", dest="dbSnpDir", metavar="SNP_DIR", help="the path to the dbSNP directory")
    parser.add_argument("--retroGenesDir", dest="retroGenesDir", metavar="RETRO_DIR", help="the path to the retrogenes directory")
    parser.add_argument("--pseudoGenesDir", dest="pseudoGenesDir", metavar="PSEUDO_DIR", help="the path to the pseudogenes directory")
    parser.add_argument("--cosmicDir", dest="cosmicDir", metavar="COSMIC_DIR", help="the path to the cosmic directory")
    parser.add_argument("--snpEffDir", dest="snpEffDir", metavar="SNP_EFF_DIR", help="the path to the snpEff directory")
    parser.add_argument("--snpEffGenome", dest="snpEffGenome", default="GRCh37.69", metavar="SNP_EFF_GENOME", help="the snpEff Genome, %default by default")
    parser.add_argument("--blatFastaFilename", dest="blatFastaFilename", metavar="FASTA_FILE", help="the fasta file that can be used during the BLAT filtering, default is the one specified in the VCF header")
    parser.add_argument("--canonical", action="store_true", default=False, dest="canonical", help="include this argument if only the canonical transcripts from snpEff should be used, %default by default")

    parser.add_argument("--rnaGeneBlckFile", dest="rnaGeneBlckFile", metavar="RNA_GENE_FILE", help="the RNA gene blacklist file")
    parser.add_argument("--rnaGeneFamilyBlckFile", dest="rnaGeneFamilyBlckFile", metavar="RNA_GENE_FAMILY_FILE", help="the RNA gene family blacklist file")

    # we do all filtering by default, so it's better for the user to specify --no flags to disable some filters
    # but internally, the code is nicer if we can avoid the double negatives, so store true by default and drop the "no" in the flag name
    parser.add_argument("--noBlacklist", action="store_false", default=True, dest="noBlacklist", help="include this argument if the blacklist filter should not be applied")
    parser.add_argument("--noTargets", action="store_false", default=True, dest="noTargets", help="include this argument if the target filter should not be applied")
    parser.add_argument("--noDbSnp", action="store_false", default=True, dest="noDbSnp", help="include this argument if the dbSNP info/filter should not be applied")
    parser.add_argument("--noRetroGenes", action="store_false", default=True, dest="noRetroGenes", help="include this argument if the info/retrogenes filter should not be applied")
    parser.add_argument("--noPseudoGenes", action="store_false", default=True, dest="noPseudoGenes", help="include this argument if the info/pseudogenes filter should not be applied")
    parser.add_argument("--noCosmic", action="store_false", default=True, dest="noCosmic", help="include this argument if the cosmic annotation should not be applied")
    parser.add_argument("--noBlat", action="store_false", default=True, dest="noBlat", help="include this argument if the blat filter should not be applied")
    parser.add_argument("--noPositionalBias", action="store_false", default=True, dest="noPositionalBias", help="include this argument if the positional bias filter should not be applied")
    parser.add_argument("--noRnaBlacklist", action="store_false", default=True, dest="noRnaBlacklist", help="include this argument if the RNA blacklist filter should not be applied")
    parser.add_argument("--noSnpEff", action="store_false", default=True, dest="noSnpEff", help="include this argument if the snpEff annotation should not be applied (without the snpEff annotation, filtering of RNA blacklisted genes will also not be applied")
    parser.add_argument("--dnaOnly", action="store_true", default=False, dest="dnaOnly", help="include this argument if you only have DNA or filtering should only be done on the DNA")
    parser.add_argument("--rnaOnly", action="store_true", default=False, dest="rnaOnly", help="include this argument if the filtering should only be done on the RNA")
    parser.add_argument("--gzip", action="store_true", default=False, dest="gzip", help="include this argument if the final VCF should be compressed with gzip")


    # some extra stuff
    parser.add_argument('--number_of_threads', dest='number_of_threads', type=int, default='1')
    parser.add_argument('--number_of_procs', dest='procs', type=int, default=1)
    parser.add_argument('--workdir', default="./")
    parser.add_argument('--no_clean', action="store_true", default=False)

    args = parser.parse_args()
    tempDir = tempfile.mkdtemp(dir="./", prefix="radia_work_")

    try:
        # if a universal fasta file is specified, then use it
        if (args.fastaFilename != None):
            universalFastaFile = indexFasta(args.workdir, args.fastaFilename, args.fastaFilename + ".fai", "universal")

        # if individual fasta files are specified, they over-ride the universal one
        if (args.dnaNormalFastaFilename != None):
            i_dnaNormalFastaFilename = indexFasta(args.workdir, args.dnaNormalFastaFilename, args.dnaNormalFastaFilename + ".fai", "dna")
        else:
            i_dnaNormalFastaFilename = universalFastaFile
        if (args.rnaNormalFastaFilename != None):
            i_rnaNormalFastaFilename = indexFasta(args.workdir, args.rnaNormalFastaFilename, args.rnaNormalFastaFilename + ".fai", "rna")
        else:
            i_rnaNormalFastaFilename = universalFastaFile
        if (args.dnaTumorFastaFilename != None):
            i_dnaTumorFastaFilename = indexFasta(args.workdir, args.dnaTumorFastaFilename, args.dnaTumorFastaFilename + ".fai", "dna")
        else:
            i_dnaTumorFastaFilename = universalFastaFile
        if (args.rnaTumorFastaFilename != None):
            i_rnaTumorFastaFilename = indexFasta(args.workdir, args.rnaTumorFastaFilename, args.rnaTumorFastaFilename + ".fai", "rna")
        else:
            i_rnaTumorFastaFilename = universalFastaFile


        if (args.dnaNormalFilename != None):
            i_dnaNormalFilename = indexBam(workdir=args.workdir, inputBamFile=args.dnaNormalFilename, inputBamFileIndex=args.dnaNormalBaiFilename, prefix="dnaNormal")
        else:
            i_dnaNormalFilename = None

        if (args.rnaNormalFilename != None):
            i_rnaNormalFilename = indexBam(workdir=args.workdir, inputBamFile=args.rnaNormalFilename, inputBamFileIndex=args.rnaNormalBaiFilename, prefix="rnaNormal")
        else:
            i_rnaNormalFilename = None

        if (args.dnaTumorFilename != None):
            i_dnaTumorFilename = indexBam(workdir=args.workdir, inputBamFile=args.dnaTumorFilename, inputBamFileIndex=args.dnaTumorBaiFilename, prefix="dnaTumor")
        else:
            i_dnaTumorFilename = None

        if (args.rnaTumorFilename != None):
            i_rnaTumorFilename = indexBam(workdir=args.workdir, inputBamFile=args.rnaTumorFilename, inputBamFileIndex=args.rnaTumorBaiFilename, prefix="rnaTumor")
        else:
            i_rnaTumorFilename = None

        # the BLAT fasta file is different from the other ones, it should contain the extra contigs and hla genes
        if (args.blatFastaFilename != None):
            blatFastaFile = indexFasta(args.workdir, args.blatFastaFilename, args.blatFastaFilename + ".fai", "blat")
        else:
            blatFastaFile = None


        if args.procs == 1:
            chroms = get_bam_seq(i_dnaNormalFilename)
            for chrom in chroms:
                # first run the RADIA command
                cmd, rawOutput = radia(chrom, args,
                            i_dnaNormalFilename, i_rnaNormalFilename, i_dnaTumorFilename, i_rnaTumorFilename,
                            i_dnaNormalFastaFilename, i_rnaNormalFastaFilename, i_dnaTumorFastaFilename, i_rnaTumorFastaFilename)
                if execute(cmd):
                    raise Exception("Radia Call failed")
                # then filter it
                cmd, out = radiaFilter(args=args, chrom=chrom, rawFile=rawOutput, outputDir=args.outputDir, scriptsDir=args.scriptsDir, blatFastaFilename=blatFastaFile)
                if execute(cmd):
                    raise Exception("RadiaFilter Call failed")
        else:
            chroms = get_bam_seq(i_dnaNormalFilename)
            cmds = []
            rawOuts = []
            for chrom in chroms:
                # create the RADIA commands
                cmd, rawOutput = radia(chrom, args,
                            i_dnaNormalFilename, i_rnaNormalFilename, i_dnaTumorFilename, i_rnaTumorFilename,
                            i_dnaNormalFastaFilename, i_rnaNormalFastaFilename, i_dnaTumorFastaFilename, i_rnaTumorFastaFilename)
                cmds.append(cmd)
                rawOuts.append(rawOutput)

            p = Pool(args.procs)
            values = p.map(execute, cmds, 1)

            cmds = []
            filteredOuts = []
            for rawOutput in rawOuts:
                # create the filter commands
                cmd, filteredOut = radiaFilter(args=args, chrom=chrom, rawFile=rawOutput, outputDir=args.outputDir, scriptsDir=args.scriptsDir, blatFastaFilename=blatFastaFile)
                cmds.append(cmd)
                filteredOuts.append(filteredOut)
            values = p.map(execute, cmds, 1)

            vcf_writer = None
            for file in filteredOuts:
                print file
                vcf_reader = vcf.Reader(filename=file)
                if vcf_writer is None:
                    vcf_writer = vcf.Writer(open(args.outputVcfFile, "w"), vcf_reader)
                for record in vcf_reader:
                    vcf_writer.write_record(record)
            vcf_writer.close()

    finally:
        if not args.no_clean and os.path.exists(tempDir):
            shutil.rmtree(tempDir)

if __name__=="__main__":
    __main__()
