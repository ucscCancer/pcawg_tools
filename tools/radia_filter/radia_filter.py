#!/usr/bin/env python

import argparse, os, shutil, subprocess, tempfile, time, gzip
from multiprocessing import Pool

def execute(cmd, output=None):
    import sys, shlex
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

def rewriteVcfGenerator(vcfline, dnaNormalFilename, dnaTumorFilename, universalFastaFile):
    """Replace filenames in vcfGenerator field to match local files."""
    # NOTE: This is NOT complete, needs rna files too
    fields = vcfline.split(',')
    for i in xrange(len(fields)):
        try:
            key, value = fields[i].split('=')
            if key == 'dnaNormalFilename':
                fields[i] = ('').join([key, "=<", dnaNormalFilename, ">"])
            if key == 'dnaNormalFastaFilename':
                fields[i] = ('').join([key, "=<", universalFastaFile, ">"])
            if key == 'dnaTumorFilename':
                fields[i] = ('').join([key, "=<", dnaTumorFilename, ">"])
            if key == 'dnaTumorFastaFilename':
                fields[i] = ('').join([key, "=<", universalFastaFile, ">"])
        except:
            continue
    newline = (',').join(fields)
    return newline
    
def splitVcf(args, outdir, dnaNormalFilename, dnaTumorFilename, universalFastaFile):
    """Splits up VCF file in chromosome files, returns dict of chromosome names and corresponding vcf files"""
    chrNames=dict()
    header=""
    f = get_read_fileHandler(args.inputVCF)
    headFlag = True
    for line in f:
        if headFlag:
	    if line.startswith('##vcfGenerator'):
                line = rewriteVcfGenerator(line, dnaNormalFilename, dnaTumorFilename, universalFastaFile)
            header += line
            if line.startswith('#CHROM'):
                headFlag = False
        else:
            fields = line.split("\t")
            if not fields[0] in chrNames:
                try:
                    o.close()
                except:
                    pass
                outfile = os.path.join(outdir, fields[0] + ".vcf" )
                o = open(outfile, 'w')      # append not necessary
                o.write(header)
                chrNames[fields[0]] = outfile
            o.write(line)
    o.close
    f.close()
    return chrNames

def radiaFilter(args, chrom, inputVcf, outputDir):

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
    justSayNo='--noBlacklist --noTargets --noDbSnp --noRetroGenes --noPseudoGenes ' + \
              '--noCosmic --noBlat --noPositionalBias --noRnaBlacklist --noSnpEff'
#    cmd = "python %s/filterRadia.py %s %s %s %s %s --blatFastaFilename %s " % (
#        args.scriptsDir,
#        args.patientId, chrom, rawFile,
#        outputDir, scriptsDir, blatFastaFilename)
    cmd = "python %s/filterRadia.py %s %s %s %s %s %s" % (
        args.scriptsDir,
        args.patientId, chrom, inputVcf,
        outputDir, args.scriptsDir, 
        justSayNo)

#    if args.canonical:
#        cmd += ' --canonical '
#    if args.noBlacklist:
#        cmd += ' --noBlacklist '
#    if args.noTargets:
#        cmd += ' --noTargets '
#    if args.noDbSnp:
#        cmd += ' --noDbSnp '
#    if args.noRetroGenes:
#        cmd += ' --noRetroGenes '
#    if args.noPseudoGenes:
#        cmd += ' --noPseudoGenes '
#    if args.noCosmic:
#        cmd += ' --noCosmic '
#    if args.noBlat:
#        cmd += ' --noBlat '
#    if args.noPositionalBias:
#        cmd += ' --noPositionalBias '
#    if args.noRnaBlacklist:
#        cmd += ' --noRnaBlacklist '
#    if args.noSnpEff:
#        cmd += ' --noSnpEff '
#    if args.dnaOnly:
#        cmd += ' --dnaOnly '
#    if args.rnaOnly:
#        cmd += ' --rnaOnly '
#    if args.gzip:
#        cmd += ' --gzip '
    outfile = os.path.join(outputDir, args.patientId + "_chr" + chrom + ".vcf")
    return cmd, outfile


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

    # tumor DNA
    parser.add_argument("-t", "--dnaTumorFilename", dest="dnaTumorFilename", metavar="DNA_TUMOR_FILE", help="the name of the tumor DNA .bam file")
    parser.add_argument("--dnaTumorBaiFilename", dest="dnaTumorBaiFilename", metavar="DNA_TUMOR_BAI_FILE", help="the name of the tumor DNA .bai file")

#    parser.add_argument("--blacklistFile", dest="blacklistFile", metavar="BLACKLIST_FILE", help="the path to the blacklist directory")

#    parser.add_argument("--targetFile", dest="targetFile", metavar="TARGET_FILE", help="the path to the exon capture targets directory")
#    parser.add_argument("--dbSnpFile", dest="dbSnpFile", metavar="SNP_FILE", help="the path to the dbSNP directory")
#    parser.add_argument("--retroGenesFile", dest="retroGenesFile", metavar="RETRO_FILE", help="the path to the retrogenes directory")
#    parser.add_argument("--pseudoGenesFile", dest="pseudoGenesFile", metavar="PSEUDO_FILE", help="the path to the pseudogenes directory")
#    parser.add_argument("--cosmicFile", dest="cosmicFile", metavar="COSMIC_FILE", help="the path to the cosmic directory")
##    parser.add_argument("--snpEffDir", dest="snpEffDir", metavar="SNP_EFF_DIR", help="the path to the snpEff directory")
#    parser.add_argument("--snpEffGenome", dest="snpEffGenome", default="GRCh37.69", metavar="SNP_EFF_GENOME", help="the snpEff Genome, %default by default")
#    parser.add_argument("--blatFastaFilename", dest="blatFastaFilename", metavar="FASTA_FILE", help="the fasta file that can be used during the BLAT filtering, default is the one specified in the VCF header")
#    parser.add_argument("--canonical", action="store_true", default=False, dest="canonical", help="include this argument if only the canonical transcripts from snpEff should be used, %default by default")
#
#    parser.add_argument("--rnaGeneBlckFile", dest="rnaGeneBlckFile", metavar="RNA_GENE_FILE", help="the RNA gene blacklist file")
#    parser.add_argument("--rnaGeneFamilyBlckFile", dest="rnaGeneFamilyBlckFile", metavar="RNA_GENE_FAMILY_FILE", help="the RNA gene family blacklist file")
#
#    # we do all filtering by default, so it's better for the user to specify --no flags to disable some filters
#    # but internally, the code is nicer if we can avoid the double negatives, so store true by default and drop the "no" in the flag name
#    parser.add_argument("--noBlacklist", action="store_false", default=True, dest="noBlacklist", help="include this argument if the blacklist filter should not be applied")
#    parser.add_argument("--noTargets", action="store_false", default=True, dest="noTargets", help="include this argument if the target filter should not be applied")
#    parser.add_argument("--noDbSnp", action="store_false", default=True, dest="noDbSnp", help="include this argument if the dbSNP info/filter should not be applied")
#    parser.add_argument("--noRetroGenes", action="store_false", default=True, dest="noRetroGenes", help="include this argument if the info/retrogenes filter should not be applied")
#    parser.add_argument("--noPseudoGenes", action="store_false", default=True, dest="noPseudoGenes", help="include this argument if the info/pseudogenes filter should not be applied")
#    parser.add_argument("--noCosmic", action="store_false", default=True, dest="noCosmic", help="include this argument if the cosmic annotation should not be applied")
#    parser.add_argument("--noBlat", action="store_false", default=True, dest="noBlat", help="include this argument if the blat filter should not be applied")
#    parser.add_argument("--noPositionalBias", action="store_false", default=True, dest="noPositionalBias", help="include this argument if the positional bias filter should not be applied")
#    parser.add_argument("--noRnaBlacklist", action="store_false", default=True, dest="noRnaBlacklist", help="include this argument if the RNA blacklist filter should not be applied")
#    parser.add_argument("--noSnpEff", action="store_false", default=True, dest="noSnpEff", help="include this argument if the snpEff annotation should not be applied (without the snpEff annotation, filtering of RNA blacklisted genes will also not be applied")
#    parser.add_argument("--dnaOnly", action="store_true", default=False, dest="dnaOnly", help="include this argument if you only have DNA or filtering should only be done on the DNA")
#    parser.add_argument("--rnaOnly", action="store_true", default=False, dest="rnaOnly", help="include this argument if the filtering should only be done on the RNA")
#    parser.add_argument("--gzip", action="store_true", default=False, dest="gzip", help="include this argument if the final VCF should be compressed with gzip")


    # some extra stuff
    parser.add_argument('--number_of_threads', dest='number_of_threads', type=int, default='1')
    parser.add_argument('--number_of_procs', dest='procs', type=int, default=1)
    parser.add_argument('--workdir', default="./")
    parser.add_argument('--no_clean', action="store_true", default=False)

    args = parser.parse_args()
    tempDir = tempfile.mkdtemp(dir="./", prefix="radia_work_")

    # All files come in as complete genome files, so first split them
        
    try:
        # if a universal fasta file is specified, then use it
        if (args.fastaFilename != None):
            universalFastaFile = indexFasta(args.workdir, args.fastaFilename, args.fastaFilename + ".fai", "universal")

        # if individual fasta files are specified, they over-ride the universal one
#        if (args.dnaNormalFastaFilename != None):
#            i_dnaNormalFastaFilename = indexFasta(args.workdir, args.dnaNormalFastaFilename, args.dnaNormalFastaFilename + ".fai", "dna")
#        else:
#            i_dnaNormalFastaFilename = universalFastaFile
#        if (args.rnaNormalFastaFilename != None):
#            i_rnaNormalFastaFilename = indexFasta(args.workdir, args.rnaNormalFastaFilename, args.rnaNormalFastaFilename + ".fai", "rna")
#        else:
#            i_rnaNormalFastaFilename = universalFastaFile
#        if (args.dnaTumorFastaFilename != None):
#            i_dnaTumorFastaFilename = indexFasta(args.workdir, args.dnaTumorFastaFilename, args.dnaTumorFastaFilename + ".fai", "dna")
#        else:
#            i_dnaTumorFastaFilename = universalFastaFile
#        if (args.rnaTumorFastaFilename != None):
#            i_rnaTumorFastaFilename = indexFasta(args.workdir, args.rnaTumorFastaFilename, args.rnaTumorFastaFilename + ".fai", "rna")
#        else:
#            i_rnaTumorFastaFilename = universalFastaFile


        if (args.dnaNormalFilename != None):
            i_dnaNormalFilename = indexBam(workdir=args.workdir, inputBamFile=args.dnaNormalFilename, inputBamFileIndex=args.dnaNormalBaiFilename, prefix="dnaNormal")
        else:
            i_dnaNormalFilename = None

#        if (args.rnaNormalFilename != None):
#            i_rnaNormalFilename = indexBam(workdir=args.workdir, inputBamFile=args.rnaNormalFilename, inputBamFileIndex=args.rnaNormalBaiFilename, prefix="rnaNormal")
#        else:
#            i_rnaNormalFilename = None

        if (args.dnaTumorFilename != None):
            i_dnaTumorFilename = indexBam(workdir=args.workdir, inputBamFile=args.dnaTumorFilename, inputBamFileIndex=args.dnaTumorBaiFilename, prefix="dnaTumor")
        else:
            i_dnaTumorFilename = None

#        if (args.rnaTumorFilename != None):
#            i_rnaTumorFilename = indexBam(workdir=args.workdir, inputBamFile=args.rnaTumorFilename, inputBamFileIndex=args.rnaTumorBaiFilename, prefix="rnaTumor")
#        else:
#            i_rnaTumorFilename = None

	# split vcf in chromosomes
	chromDict=splitVcf(args, args.workdir, i_dnaNormalFilename, i_dnaTumorFilename, universalFastaFile)

#        time.sleep(3600)
        rfOuts = []
        if args.procs == 1:
            for chrom in chromDict:
                chrom = '21'	# DEBUG setting
                cmd, rfOutput = radiaFilter(args, chrom, chromDict[chrom], args.workdir)
		# the output is generated by the filter, not on stdout
                if execute(cmd):
                    raise Exception("RadiaFilter Call failed")
		rfOuts.append(rfOutput)
                move(rfOutput, args.outputFilename)
                break
        else:
            cmds = []
            rawOuts = []
#                cmds.append(cmd)
#                rawOuts.append(rawOutput)

#            p = Pool(args.procs)
#            values = p.map(execute, cmds, 1)

        # even though we have a list of radia output files, we don't really need it:
        # the radiaMerge command only uses the working directory and patient name

#        cmd = radiaMerge(args)
#        if execute(cmd):
#            raise Exception("RadiaMerge Call failed")

    finally:
        if not args.no_clean and os.path.exists(tempDir):
            shutil.rmtree(tempDir)

if __name__=="__main__":
    __main__()
