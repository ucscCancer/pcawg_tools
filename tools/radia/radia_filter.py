#!/usr/bin/env python

import argparse, os, shutil, subprocess, tempfile, time, gzip, zipfile, sys
from multiprocessing import Pool

def execute(cmd, output=None):
    import shlex
    # function to execute a cmd and report if an error occurs
    print(cmd)
    try:
        process = subprocess.Popen(args=shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = process.communicate()
    except Exception, e: # error from my command : stderr
        sys.stderr.write("problem doing : %s\n%s\n" %(cmd, e))
        return 1
    if output:
        output = open(output, 'w')
        output.write(stdout)
        output.close()
    if stderr != '': # internal program error : stdout 
        sys.stdout.write("warning or error while doing : %s\n-----\n%s-----\n\n" %(cmd, stderr))
        return 1
    return 0

def move(avant, apres):
    if os.path.exists(avant):
        execute("mv %s %s" %(avant, apres))

def correctLineCount(number, vcfFile):
    """Checks number of non header lines in input VCF."""
    f = get_read_fileHandler(vcfFile)
    count=0
    for line in f:
        if line.startswith('#'):
            continue
        count += 1
    f.close()
    if count != number:
        sys.stdout.write("ERROR expected %d non header lines in %s, got %d\n" % (vcfFile, count, number))
        return False
    return True

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
    """Checks if fasta index exists. If so, creates link. If not, creates index"""
    inputFastaLink = os.path.join(os.path.abspath(workdir), prefix + "_reference.fa" )
    os.symlink(inputFastaFile, inputFastaLink)
    inputFastaFileIndex = inputFastaFile + ".fai"
    if os.path.exists(inputFastaFileIndex):
        os.symlink(inputFastaFileIndex, inputFastaLink + ".fai")
    else:
        cmd = "samtools faidx %s" %(inputFastaLink)
        execute([cmd])
    return inputFastaLink

def get_read_fileHandler(aFilename):
    """ Open aFilename for reading and return the file handler.  The file can be gzipped or not."""
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'rb')
    else:
        return open(aFilename,'r')

def get_write_fileHandler(aFilename):
    """ Open aFilename for writing and return the file handler.  The file can be gzipped or not."""
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'wb')
    else:
        return open(aFilename,'w')

def rewriteVcfGenerator(vcfline, files):
    """Replace filenames in vcfGenerator field to match local files."""
    fields = vcfline.split(',')
    for i in xrange(len(fields)):
        try:
            key, value = fields[i].split('=')
        except:
            continue
        if key == 'dnaNormalFilename':
            if files.dnaNormalFilename == None:
                sys.stderr.write("VCF header contains DNA normal bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.dnaNormalFilename, ">"])
        elif key == 'dnaTumorFilename':
            if files.dnaTumorFilename == None:
                sys.stderr.write("VCF header contains DNA tumor bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.dnaTumorFilename, ">"])
        elif key == 'rnaNormalFilename':
            if files.rnaNormalFilename == None:
                sys.stderr.write("VCF header contains RNA normal bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.rnaNormalFilename, ">"])
        elif key == 'rnaTumorFilename':
            if files.rnaTumorFilename == None:
                sys.stderr.write("VCF header contains RNA tumor bam, please input corresponding bam file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.rnaTumorFilename, ">"])
        elif key == 'dnaNormalFastaFilename':
            if files.dnaNormalFastaFilename == None:
                sys.stderr.write("VCF header contains DNA normal fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.dnaNormalFastaFilename, ">"])
        elif key == 'dnaTumorFastaFilename':
            if files.dnaTumorFastaFilename == None:
                sys.stderr.write("VCF header contains DNA tumor fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.dnaTumorFastaFilename, ">"])
        elif key == 'rnaNormalFastaFilename':
            if files.rnaNormalFastaFilename == None:
                sys.stderr.write("VCF header contains RNA normal fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.rnaNormalFastaFilename, ">"])
        elif key == 'rnaTumorFastaFilename':
            if files.rnaTumorFastaFilename == None:
                sys.stderr.write("VCF header contains RNA tumor fasta, please input corresponding fasta file")
                sys.exit(1)
            fields[i] = ('').join([key, "=<", files.rnaTumorFastaFilename, ">"])
    newline = (',').join(fields)
    return newline
    
def splitVcf(infile, outdir, files=None, expected=None):
    """Splits up VCF file in chromosome files, a dict of chromosome names with non header linecounts and a dict of chromosome names (without chr) and corresponding vcf files (with chr)"""
    chrNames=dict()
    chrLines=dict()
    header=""
    f = get_read_fileHandler(infile)
    headFlag = True
    for line in f:
        if headFlag:
	    if files != None and line.startswith('##vcfGenerator'):
                line = rewriteVcfGenerator(line, files)
            header += line
            if line.startswith('#CHROM'):
                headFlag = False
        else:
            fields = line.split("\t")
            chrom = fields[0].replace('chr', '')	# remove chr if present
            if not chrom in chrNames:
                try:
                    o.close()
                except:
                    pass
                outfile = os.path.join(outdir, 'chr' + chrom + ".vcf.gz" )
                o = get_write_fileHandler(outfile)      # append not necessary
                o.write(header)
                chrLines[chrom] = 0
                chrNames[chrom] = outfile
            o.write(line)
            chrLines[chrom] += 1
    o.close
    f.close()
    if expected != None:
        wanted = expected.keys()
        wanted = [s.replace('chr','') for s in wanted]	# remove chr
        created = set(chrNames.keys())
        for chrom in set(wanted).difference(created):
            sys.stderr.write("WARNING, missing chromosome %s in filter %s, ignoring...\n" % (chrom, outdir))
            # creating empty file
            outfile = os.path.join(outdir, "chr" + chrom + ".vcf" )
            o = get_write_fileHandler(outfile)
            o.close
    return chrNames, chrLines

def splitBed(bedfile, outdir, chromDict):
    """Splits bed file in chromosome files, issued warnings on missing files and creates empties"""
    wanted = chromDict.keys()
    wanted = [s.replace('chr','') for s in wanted]	# remove chr
    wanted = ['chr'+s for s in wanted]	# add chr
    chrNames = set()
    header=""
    f = get_read_fileHandler(bedfile)
    for line in f:
            fields = line.split("\t")
            chrom = fields[0]
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            if not chrom in chrNames:
                try:
                    o.close()
                except:
                    pass
                outfile = os.path.join(outdir, chrom + ".bed.gz" )
                o = get_write_fileHandler(outfile)      # append not necessary
                chrNames.add(chrom)
            o.write(line)
    o.close
    f.close()
    for chrom in set(wanted).difference(chrNames):
        sys.stderr.write("WARNING, missing chromosome %s in filter %s, ignoring...\n" % (chrom, outdir))
        # creating empty file
        outfile = os.path.join(outdir, chrom + ".bed") 
        o = get_write_fileHandler(outfile)
        o.close

def identicalName(inputList):
    """returns duplicate name if two inputs have the same name and are not None"""
    dup=set(x for x in inputList if inputList.count(x) >= 2)
    dup.discard(None)   # this doesn't complain if None is not in the set
    if dup:
        print "ERROR: found duplicate input %s" % dup.pop()
        return True
    return False

def makeSnpEffConfig(workdir, genome, datadir):
    """Creates a short config file for snpEff. Assumes a human genome."""
    configFile = os.path.join(workdir, "snpEff.config")
    f = open(configFile, 'w')
    f.write("data.dir = %s\n" % datadir)
    f.write("lof.ignoreProteinCodingAfter : 0.95\n")
    f.write("lof.ignoreProteinCodingBefore : 0.05\n")
    f.write("lof.deleteProteinCodingBases : 0.50\n")
    f.write("codon.Standard : TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\n")
    f.write("codon.Vertebrate_Mitochondrial : TTT/F, TTC/F, TTA/L, TTG/L, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/W, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I+, ATC/I+, ATA/M+, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/*, AGG/*, GTT/V, GTC/V, GTA/V, GTG/V+, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\n")
    f.write("%s.genome : Homo_sapiens\n" % genome)
    f.close()
    return configFile

def radiaFilter(filterDirs, snpEffGenome, snpEffConfig, args, chrom, inputVcf, outputDir, logFile ):

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
    cmd = "python %s/filterRadia.py %s %s %s %s %s --gzip --log=WARNING -g %s" % (
        args.scriptsDir,
        args.patientId, chrom, inputVcf,
        outputDir, args.scriptsDir, 
        logFile)

    if args.blatFastaFilename != None:
        cmd += " --blatFastaFilename %s" % args.blatFastaFilename
    else:
        cmd += ' --noBlat'

    if "blacklist" in filterDirs:
        cmd += " --blacklistDir %s" % filterDirs["blacklist"]
    else:
        cmd += ' --noBlacklist'

    if "target" in filterDirs:
        cmd += " --targetDir %s" % filterDirs["target"]
    else:
        cmd += ' --noTargets'

    if "snp" in filterDirs:
        cmd += " --dbSnpDir %s" % filterDirs["snp"]
    else:
        cmd += ' --noDbSnp'

    if "pseudoGenes" in filterDirs:
        cmd += " --pseudoGenesDir %s" % filterDirs["pseudoGenes"]
    else:
        cmd += ' --noPseudoGenes'

    if "retroGenes" in filterDirs:
        cmd += " --retroGenesDir %s" % filterDirs["retroGenes"]
    else:
        cmd += ' --noRetroGenes'

    if "cosmic" in filterDirs:
        cmd += " --cosmicDir %s" % filterDirs["cosmic"]
    else:
        cmd += ' --noCosmic'

    if snpEffGenome != None:
        cmd += " --snpEffDir %s --snpEffGenome %s --snpEffConfig %s" % (args.snpEffDir, snpEffGenome, snpEffConfig)
        if args.canonical:
            cmd += ' --canonical '
        if args.rnaGeneBlckFile:
            cmd += ' --rnaGeneBlckFile %s --rnaGeneFamilyBlckFile %s' % (args.rnaGeneBlckFile, args.rnaGeneFamilyBlckFile)
        else:
            cmd += ' --noRnaBlacklist'
    else:
        cmd += ' --noSnpEff --noRnaBlacklist'

    if args.noPositionalBias:
        cmd += ' --noPositionalBias'
    if args.dnaOnly:
        cmd += ' --dnaOnly'
#    if args.rnaOnly:
#        cmd += ' --rnaOnly '
#    if args.gzip:
#        cmd += ' --gzip '
    outfile = os.path.join(outputDir, args.patientId + "_chr" + chrom + ".vcf.gz")
    return cmd, outfile

def radiaMerge(args, inputDir):
    """Merges vcf files if they follow the pattern patientID_chr<N>.vcf(.gz)"""

    # python mergeChroms.py patientId /radia/filteredChroms/ /radia/filteredPatients/ --gzip
    #  -h, --help            show this help message and exit
    #  -o OUTPUT_FILE, --outputFilename=OUTPUT_FILE
    #                   the name of the output file, <id>.vcf(.gz) by default

    #  -l LOG, --log=LOG     the logging level (DEBUG, INFO, WARNING, ERROR,
    #                        CRITICAL), WARNING by default
    #  -g LOG_FILE, --logFilename=LOG_FILE
    #                        the name of the log file, STDOUT by default
    #  --gzip                include this argument if the final VCF should be
    #                        compressed with gzip
    # radia works in the workdir
    cmd = "python %s/mergeChroms.py %s %s %s -o %s" % (
        args.scriptsDir,
        args.patientId, inputDir, args.workdir,
        args.outputFilename)
#    if args.gzip:
#        cmd += ' --gzip'
    return cmd


class localFiles(object):
    """Formats input fasta and bam files and creates local filenames"""
    def __init__(self):
        self.dnaNormalFilename, self.dnaTumorFilename, self.rnaNormalFilename, self.rnaTumorFilename, self.dnaNormalFastaFilename, self.dnaTumorFastaFilename, self.rnaNormalFastaFilename, self.rnaTumorFastaFilename = (None for i in xrange(8))


    def universalFasta(self, args):
        if (args.fastaFilename != None):
            universalFastaFile = indexFasta(args.workdir, args.fastaFilename, prefix="universal")
            self.dnaNormalFastaFilename, self.dnaTumorFastaFilename, self.rnaNormalFastaFilename, self.rnaTumorFastaFilename = (universalFastaFile for i in xrange(4))


    def doFasta(self, args):
        # if individual fasta files are specified, they over-ride the universal one
        if (args.dnaNormalFastaFilename != None):
            self.dnaNormalFastaFilename = indexFasta(args.workdir, args.dnaNormalFastaFilename, prefix="dnaN")

        if (args.rnaNormalFastaFilename != None):
            self.rnaNormalFastaFilename = indexFasta(args.workdir, args.rnaNormalFastaFilename, prefix="rnaN")
        if (args.dnaTumorFastaFilename != None):
            self.dnaTumorFastaFilename = indexFasta(args.workdir, args.dnaTumorFastaFilename, prefix="dnaT")
        if (args.rnaTumorFastaFilename != None):
            self.rnaTumorFastaFilename = indexFasta(args.workdir, args.rnaTumorFastaFilename, prefix="rnaT")


    def doBam(self, args):
        # index bam files and return True if there's only DNA files
        if (args.dnaNormalFilename != None):
            self.dnaNormalFilename = indexBam(workdir=args.workdir, inputBamFile=args.dnaNormalFilename, inputBamFileIndex=args.dnaNormalBaiFilename, prefix="dnaNormal")

        if (args.rnaNormalFilename != None):
            self.rnaNormalFilename = indexBam(workdir=args.workdir, inputBamFile=args.rnaNormalFilename, inputBamFileIndex=args.rnaNormalBaiFilename, prefix="rnaNormal")

        if (args.dnaTumorFilename != None):
            self.dnaTumorFilename = indexBam(workdir=args.workdir, inputBamFile=args.dnaTumorFilename, inputBamFileIndex=args.dnaTumorBaiFilename, prefix="dnaTumor")

        if (args.rnaTumorFilename != None):
            self.rnaTumorFilename = indexBam(workdir=args.workdir, inputBamFile=args.rnaTumorFilename, inputBamFileIndex=args.rnaTumorBaiFilename, prefix="rnaTumor")
        if (args.rnaTumorFilename == None and args.rnaNormalFilename == None):
            return True
        return False

def __main__():
    time.sleep(1) #small hack, sometimes it seems like docker file systems are avalible instantly
    parser = argparse.ArgumentParser(description="RADIA filter")


    #############################
    #    RADIA filter params    #
    #############################
    parser.add_argument("--inputVCF", dest="inputVCF", required=True, metavar="INPUT_VCF", help="The input Radia vcf file")
    parser.add_argument("--patientId", dest="patientId", required=True, metavar="PATIENT_ID", help="a unique patient Id that will be used to name the output file")
    parser.add_argument("-o", "--outputFilename", dest="outputFilename", required=True, metavar="OUTPUT_FILE", default='out.vcf', help="the name of the output file")
    parser.add_argument("--outputDir", dest="outputDir", required=True, metavar="FILTER_OUT_DIR", help="the directory where temporary and final filtered output should be stored")
    parser.add_argument("--scriptsDir", dest="scriptsDir", required=True, metavar="SCRIPTS_DIR", help="the directory that contains the RADIA filter scripts")
    parser.add_argument("-f", "--fastaFilename", dest="fastaFilename", metavar="FASTA_FILE", help="the name of the fasta file that can be used on all .bams, see below for specifying individual fasta files for each .bam file")

    # normal DNA
    parser.add_argument("-n", "--dnaNormalFilename", dest="dnaNormalFilename", metavar="DNA_NORMAL_FILE", help="the name of the normal DNA .bam file")
    parser.add_argument("--dnaNormalBaiFilename", dest="dnaNormalBaiFilename", metavar="DNA_NORMAL_BAI_FILE", help="the name of the normal DNA .bai file")
    parser.add_argument ("--dnaNormalFastaFilename", dest="dnaNormalFastaFilename", metavar="DNA_NORMAL_FASTA_FILE", help="the name of the fasta file that was used to create the BAM alignments")

    # tumor DNA
    parser.add_argument("-t", "--dnaTumorFilename", dest="dnaTumorFilename", metavar="DNA_TUMOR_FILE", help="the name of the tumor DNA .bam file")
    parser.add_argument("--dnaTumorBaiFilename", dest="dnaTumorBaiFilename", metavar="DNA_TUMOR_BAI_FILE", help="the name of the tumor DNA .bai file")
    parser.add_argument("--dnaTumorFastaFilename", dest="dnaTumorFastaFilename", metavar="DNA_TUMOR_FASTA_FILE", help="the name of the fasta file that was used to create the BAM alignments")

    # normal RNA
    parser.add_argument("-x", "--rnaNormalFilename", dest="rnaNormalFilename", metavar="RNA_NORMAL_FILE", help="the name of the normal RNA .bam file")
    parser.add_argument("--rnaNormalBaiFilename", dest="rnaNormalBaiFilename", metavar="RNA_NORMAL_BAI_FILE", help="the name of the normal RNA .bai file")
    parser.add_argument ("--rnaNormalFastaFilename", dest="rnaNormalFastaFilename", metavar="RNA_NORMAL_FASTA_FILE", help="the name of the fasta file that was used to create the BAM alignments")

    # tumor RNA
    parser.add_argument("-r", "--rnaTumorFilename", dest="rnaTumorFilename", metavar="RNA_TUMOR_FILE", help="the name of the tumor RNA .bam file")
    parser.add_argument("--rnaTumorBaiFilename", dest="rnaTumorBaiFilename", metavar="RNA_TUMOR_BAI_FILE", help="the name of the tumor RNA .bai file")
    parser.add_argument("--rnaTumorFastaFilename", dest="rnaTumorFastaFilename", metavar="RNA_TUMOR_FASTA_FILE", help="the name of the fasta file that was used to create the BAM alignments")

    parser.add_argument("--blacklistFilename", dest="blacklistFilename", metavar="BLACKLIST_FILE", help="the name of the blacklist bed file")

    parser.add_argument("--targetFilename", dest="targetFilename", metavar="TARGET_FILE", help="the name of the exon capture targets file")
    parser.add_argument("--snpFilename", dest="snpFilename", metavar="SNP_FILE", help="dbSNP vcf file")
    parser.add_argument("--retroGenesFilename", dest="retroGenesFilename", metavar="RETRO_FILE", help="the name of the retrogenes bed file")
    parser.add_argument("--pseudoGenesFilename", dest="pseudoGenesFilename", metavar="PSEUDO_FILE", help="the name of the pseudogenes bed file")
    parser.add_argument("--cosmicFilename", dest="cosmicFilename", metavar="COSMIC_FILE", help="the name of the Catalogue Of Somatic Mutations In Cancer (COSMIC) annotations file")
    parser.add_argument("--snpEffDir", dest="snpEffDir", metavar="SNP_EFF_DIR", help="the path to the snpEff directory")
    parser.add_argument("--snpEffFilename", dest="snpEffFilename", metavar="SNP_EFF_FILE", help="the snpEff input database zip file")
    parser.add_argument("--canonical", action="store_true", default=False, dest="canonical", help="include this argument if only the canonical transcripts from snpEff should be used, %default by default")
    parser.add_argument("--rnaGeneBlckFile", dest="rnaGeneBlckFile", metavar="RNA_GENE_FILE", help="the RNA gene blacklist file")
    parser.add_argument("--rnaGeneFamilyBlckFile", dest="rnaGeneFamilyBlckFile", metavar="RNA_GENE_FAMILY_FILE", help="the RNA gene family blacklist file")
    parser.add_argument("--blatFastaFilename", dest="blatFastaFilename", metavar="FASTA_FILE", help="the fasta file that can be used during the BLAT filtering")
    parser.add_argument("--noPositionalBias", action="store_false", default=True, dest="noPositionalBias", help="include this argument if the positional bias filter should not be applied")

    parser.add_argument("--dnaOnly", action="store_true", default=False, dest="dnaOnly", help="include this argument if you only have DNA or filtering should only be done on the DNA")
#    parser.add_argument("--rnaOnly", action="store_true", default=False, dest="rnaOnly", help="include this argument if the filtering should only be done on the RNA")
#    parser.add_argument("--gzip", action="store_true", default=False, dest="gzip", help="include this argument if the final VCF should be compressed with gzip")


    # some extra stuff
    parser.add_argument('--number_of_procs', dest='procs', type=int, default=1)
    parser.add_argument('--workdir', default="./")
    parser.add_argument('--no_clean', action="store_true", default=False)

    args = parser.parse_args()
    tempDir = tempfile.mkdtemp(dir="./", prefix="radia_work_")

    # sanity checks
    if identicalName([args.dnaNormalFilename, args.dnaTumorFilename, args.rnaNormalFilename, args.rnaTumorFilename]):
            raise Exception("ERROR: Found duplicate input bam file")
    if identicalName([args.blacklistFilename, args.targetFilename, args.retroGenesFilename, 
        args.pseudoGenesFilename, args.cosmicFilename]):
            raise Exception("ERROR: Found duplicate input bed file")
    if (args.rnaGeneFamilyBlckFile and not args.rnaGeneBlckFile) or (args.rnaGeneBlckFile and not args.rnaGeneFamilyBlckFile):
            raise Exception("ERROR: Must input two RNA blacklist files")
    if identicalName([args.rnaGeneFamilyBlckFile, args.rnaGeneBlckFile]):
            raise Exception("ERROR: Found duplicate input RNA blacklist file")

    files = localFiles()	# prepares and holds fasta and bam files
    filterDirs = dict()		# holds filters and file locations
    try:
        files.universalFasta(args)
        files.doFasta(args)
        args.dnaOnly = files.doBam(args)
	# split vcf in chromosomes
	chromDict, chromLines=splitVcf(args.inputVCF, args.workdir, files=files)

        # All files come in as complete genome files, so first split them
	# split blacklist
        if (args.blacklistFilename != None):
            blacklistDir = os.path.join(args.workdir, "blacklistDir")
            os.mkdir(blacklistDir)
            splitBed(args.blacklistFilename, blacklistDir, chromDict)
            filterDirs["blacklist"] = blacklistDir
	# split target
        if (args.targetFilename != None):
            targetDir = os.path.join(args.workdir, "targetDir")
            os.mkdir(targetDir)
            splitBed(args.targetFilename, targetDir, chromDict)
            filterDirs["target"] = targetDir
	# split snp
        if (args.snpFilename != None):
            snpDir = os.path.join(args.workdir, "snpDir")
            os.mkdir(snpDir)
            splitVcf(args.snpFilename, snpDir, expected=chromDict)
            filterDirs["snp"] = snpDir
	# split retrogenes
        if (args.retroGenesFilename != None):
            retroGenesDir = os.path.join(args.workdir, "retroGenesDir")
            os.mkdir(retroGenesDir)
            splitBed(args.retroGenesFilename, retroGenesDir, chromDict)
            filterDirs["retroGenes"] = retroGenesDir
	# split pseudogenes
        if (args.pseudoGenesFilename != None):
            pseudoGenesDir = os.path.join(args.workdir, "pseudoGenesDir")
            os.mkdir(pseudoGenesDir)
            splitBed(args.pseudoGenesFilename, pseudoGenesDir, chromDict)
            filterDirs["pseudoGenes"] = pseudoGenesDir
	# split cosmic
        if (args.cosmicFilename != None):
            cosmicDir = os.path.join(args.workdir, "cosmicDir")
            os.mkdir(cosmicDir)
            splitBed(args.cosmicFilename, cosmicDir, chromDict)
            filterDirs["cosmic"] = cosmicDir
        # setup snpEff database
	if (args.snpEffFilename):
            with zipfile.ZipFile(args.snpEffFilename, "r") as z:
                z.extractall(args.workdir)
            # this creates a directory named data, which holds the genome directory
            # the name of that directory is used by snpEff
            datadir = os.path.join(args.workdir, "data")
            snpEffGenome = os.listdir(datadir)[0]
            snpEffConfig = makeSnpEffConfig(args.workdir, snpEffGenome, datadir)
        else:
            snpEffGenome = None
            snpEffConfig = None

        rfOuts = []
        if args.procs == 1:
            for chrom in chromDict:
                logFile = os.path.join(args.workdir, "log." + chrom)
                cmd,outfile = radiaFilter(filterDirs, snpEffGenome, snpEffConfig, args, chrom, chromDict[chrom], tempDir, logFile)
		# the output is generated by the filter, not on stdout
                if execute(cmd):
                    raise Exception("RadiaFilter Call failed")
                if not correctLineCount(chromLines[chrom], outfile):
                    raise Exception("RadiaFilter sanity check failed")
                with open(logFile, 'r') as f:
                    print >>sys.stderr, f.read()
        else:
            cmds = []
            rawOuts = dict()
            for chrom in chromDict:
                logFile = os.path.join(args.workdir, "log." + chrom)
                cmd, outfile = radiaFilter(filterDirs, snpEffGenome, snpEffConfig, args, chrom, chromDict[chrom], tempDir, logFile)
                cmds.append(cmd)
                rawOuts[chrom] = outfile
            p = Pool(args.procs)
            values = p.map(execute, cmds, 1)
            # check if all output files are the same size as inputs
            # and print logging info to stderr
            for chrom in chromDict:
                logFile = os.path.join(args.workdir, "log." + chrom)
                with open(logFile, 'r') as f:
                    print >>sys.stderr, f.read()
                if not correctLineCount(chromLines[chrom], rawOuts[chrom]):
                    raise Exception("RadiaFilter sanity check failed")




        # the radiaMerge command only uses the output directory and patient name
        cmd = radiaMerge(args, tempDir)
        if execute(cmd):
            raise Exception("RadiaMerge Call failed")


    finally:
        args.no_clean = True
        if not args.no_clean and os.path.exists(tempDir):
            shutil.rmtree(tempDir)

if __name__=="__main__":
    __main__()
