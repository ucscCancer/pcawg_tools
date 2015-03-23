#!/usr/bin/env python

import argparse, os, shutil, subprocess, sys, tempfile, time, shlex
from multiprocessing import Pool
import vcf

def execute(cmd, output=None):
	import subprocess, sys, shlex
	# function to execute a cmd and report if an error occur
	print(cmd)
	try:
		process = subprocess.Popen(args=shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout,stderr = process.communicate()
	except Exception, e: # une erreur de ma commande : stderr
		sys.stderr.write("problem doing : %s\n%s\n" %(cmd, e))
		return
	if output:
		output = open(output, 'w')
		output.write(stdout)
		output.close()
	if stderr != '': # une erreur interne au programme : stdout (sinon, souvent des warning arrete les programmes)
		sys.stdout.write("warning or error while doing : %s\n-----\n%s-----\n\n" %(cmd, stderr))


def indexBam(workdir, inputFastaFile, inputBamFile, inputBamFileIndex=None):
	inputFastaLink = os.path.join(os.path.abspath(workdir), "reference.fa" )
	os.symlink(inputFastaFile, inputFastaLink)
	inputBamLink = os.path.join(os.path.abspath(workdir), "sample.bam" )
	os.symlink(inputBamFile, inputBamLink)
	if inputBamFileIndex is None:
		cmd = "samtools index %s" %(inputBamLink)
		execute(cmd)
	else:
		os.symlink(inputBamFileIndex, inputBamLink + ".bai")
	cmd = "samtools faidx %s" %(inputFastaLink)
	execute(cmd)
	return inputFastaLink, inputBamLink


def config(inputBamFile, meanInsertSize, tag, tempDir):
	configFile = tempDir+"/pindel_configFile"
	fil = open(configFile, 'w')
	fil.write("%s\t%s\t%s\n" %(inputBamFile, meanInsertSize, tag))
	fil.close()
	return configFile


def pindel(reference, configFile, args, tempDir, chrome=None):
	if chrome is None:
		pindelTempDir=tempDir + "/pindel"
	else:
		pindelTempDir=tempDir + "/pindel_" + chrome

	cmd = "pindel -f %s -i %s -o %s " %(reference, configFile, pindelTempDir)
	cmd += " --number_of_threads %d --max_range_index %d --window_size %d --sequencing_error_rate %f --sensitivity %f" %(args.number_of_threads, args.max_range_index, args.window_size, args.sequencing_error_rate, args.sensitivity)
	cmd += " -u %f -n %d -a %d -m %d -v %d -d %d -B %d -A %d -M %d " %(args.maximum_allowed_mismatch_rate, args.NM, args.additional_mismatch, args.min_perfect_match_around_BP, args.min_inversion_size, args.min_num_matched_bases, args.balance_cutoff, args.anchor_quality, args.minimum_support_for_event)

	if chrome is not None:
		cmd += "-c %s " % (chrome)

	if args.report_long_insertions:
		cmd += ' --report_long_insertions '
	if args.report_duplications:
		cmd += ' --report_duplications '
	if args.report_inversions:
		cmd += ' --report_inversions '
	if args.report_breakpoints:
		cmd += ' --report_breakpoints '

	if args.report_close_mapped_reads:
		cmd += ' --report_close_mapped_reads '
	if args.report_only_close_mapped_reads:
		cmd += ' --report_only_close_mapped_reads '
	if args.report_interchromosomal_events:
		cmd += ' --report_interchromosomal_events '
	if args.IndelCorrection:
		cmd += ' --IndelCorrection '
	if args.NormalSamples:
		cmd += ' --NormalSamples '
	if args.input_SV_Calls_for_assembly:
		cmd += ' --input_SV_Calls_for_assembly %s ' %(args.input_SV_Calls_for_assembly)

	if args.exclude is not None:
		cmd += '--exclude %s' % (args.exclude)

	if args.detect_DD:
		cmd += ' -q '
		cmd += ' --MAX_DD_BREAKPOINT_DISTANCE '+str(args.MAX_DD_BREAKPOINT_DISTANCE)
		cmd += ' --MAX_DISTANCE_CLUSTER_READS '+str(args.MAX_DISTANCE_CLUSTER_READS)
		cmd += ' --MIN_DD_CLUSTER_SIZE '+str(args.MIN_DD_CLUSTER_SIZE)
		cmd += ' --MIN_DD_BREAKPOINT_SUPPORT '+str(args.MIN_DD_BREAKPOINT_SUPPORT)
		cmd += ' --MIN_DD_MAP_DISTANCE '+str(args.MIN_DD_MAP_DISTANCE)
		if args.DD_REPORT_DUPLICATION_READS:
			cmd += ' --DD_REPORT_DUPLICATION_READS '

	return (cmd, pindelTempDir)


def move(avant, apres):
	if os.path.exists(avant):
		execute("mv %s %s" %(avant, apres))


def pindel2vcf(inputFastaFile, name, pindelTempDir, chrome=None):
	date = str(time.strftime('%d/%m/%y',time.localtime()))
	if chrome is None:
		cmd = "pindel2vcf -P %s -r %s -R %s -d %s" %(pindelTempDir, inputFastaFile, name, date)
		return (cmd, pindelTempDir+".vcf")
	else:
		output = "%s.%s.vcf" % (pindelTempDir, chrome)
		cmd = "pindel2vcf -P %s -r %s -R %s -d %s -v %s" %(pindelTempDir, inputFastaFile, name, date, output)
		return (cmd, output)



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


# def getMeanInsertSize(bamFile, tempDir):
# 	samFile = tempDir+"/sam"
# 	histogramFile = tempDir+"/histogram"
# 	outputFile = tempDir+"/output"
# 	cmd = "samtools view -h %s -o %s" % (bamFile, samFile)
# 	execute(cmd, )
# 	cmd = "picard-tools CollectInsertSizeMetrics H=%s O=%s I=%s" % (histogramFile, outputFile, samFile)
# 	execute(cmd)
# 	try:
# 		f = open(outputFile, 'r')
# 		lignes = f.readlines()
# 		f.close()
# 		meanInsertSize = lignes[7].split('\t')[0]
# 	except:
# 		sys.stderr.write("problem while opening file %s " %(outputFile))
# 		return 0
# 	return meanInsertSize


def __main__():
	time.sleep(1) #small hack, sometimes it seems like docker file systems aren't avalible instantly
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-r', dest='inputFastaFile', required=True, help='the reference file')
	parser.add_argument('-b', dest='inputBamFile', required=True, help='the bam file')
	parser.add_argument('-bi', dest='inputBamFileIndex', default=None, help='the bam file')
	parser.add_argument('-s', dest='insert_size', type=int, default=None, required=False, help='the insert size')
	parser.add_argument('-t', dest='sampleTag', default='default', required=False, help='the sample tag')
	# parser.add_argument('-o1', dest='outputRaw', help='the output raw', required=True)
	parser.add_argument('-o2', dest='outputVcfFile', help='the output vcf', required=True)
	parser.add_argument('--number_of_threads', dest='number_of_threads', type=int, default='1')
	parser.add_argument('--number_of_procs', dest='procs', type=int, default=1)

	parser.add_argument('-x', '--max_range_index', dest='max_range_index', type=int, default='4')
	parser.add_argument('--window_size', dest='window_size', type=int, default='5')
	parser.add_argument('--sequencing_error_rate', dest='sequencing_error_rate', type=float, default='0.01')
	parser.add_argument('--sensitivity', dest='sensitivity', default='0.95', type=float)
	parser.add_argument('--report_long_insertions', dest='report_long_insertions', action='store_true', default=False)
	parser.add_argument('--report_duplications', dest='report_duplications', action='store_true', default=False)
	parser.add_argument('--report_inversions', dest='report_inversions', action='store_true', default=False)
	parser.add_argument('--report_breakpoints', dest='report_breakpoints', action='store_true', default=False)
	parser.add_argument('-u', '--maximum_allowed_mismatch_rate', dest='maximum_allowed_mismatch_rate', type=float, default='0.02')
	parser.add_argument('--report_close_mapped_reads', dest='report_close_mapped_reads', action='store_true', default=False)
	parser.add_argument('--report_only_close_mapped_reads', dest='report_only_close_mapped_reads', action='store_true', default=False)
	parser.add_argument('--report_interchromosomal_events', dest='report_interchromosomal_events', action='store_true', default=False)
	parser.add_argument('--IndelCorrection', dest='IndelCorrection', action='store_true', default=False)
	parser.add_argument('--NormalSamples', dest='NormalSamples', action='store_true', default=False)
	parser.add_argument('-a', '--additional_mismatch', dest='additional_mismatch', type=int, default='1')
	parser.add_argument('-m', '--min_perfect_match_around_BP', dest='min_perfect_match_around_BP', type=int, default='3')
	parser.add_argument('-v', '--min_inversion_size', dest='min_inversion_size', type=int, default='50')
	parser.add_argument('-d', '--min_num_matched_bases', dest='min_num_matched_bases', type=int, default='30')
	parser.add_argument('-B', '--balance_cutoff', dest='balance_cutoff', type=int, default='0')
	parser.add_argument('-A', '--anchor_quality', dest='anchor_quality', type=int, default='0')
	parser.add_argument('-M', '--minimum_support_for_event', dest='minimum_support_for_event', type=int, default='3')
	parser.add_argument('-n', '--NM', dest='NM', type=int, default='2')
	parser.add_argument('--detect_DD', dest='detect_DD', action='store_true', default=False)
	parser.add_argument('--MAX_DD_BREAKPOINT_DISTANCE', dest='MAX_DD_BREAKPOINT_DISTANCE', type=int, default='350')
	parser.add_argument('--MAX_DISTANCE_CLUSTER_READS', dest='MAX_DISTANCE_CLUSTER_READS', type=int, default='100')
	parser.add_argument('--MIN_DD_CLUSTER_SIZE', dest='MIN_DD_CLUSTER_SIZE', type=int, default='3')
	parser.add_argument('--MIN_DD_BREAKPOINT_SUPPORT', dest='MIN_DD_BREAKPOINT_SUPPORT', type=int, default='3')
	parser.add_argument('--MIN_DD_MAP_DISTANCE', dest='MIN_DD_MAP_DISTANCE', type=int, default='8000')
	parser.add_argument('--DD_REPORT_DUPLICATION_READS', dest='DD_REPORT_DUPLICATION_READS', action='store_true', default=False)

	parser.add_argument("-J", "--exclude", dest="exclude", default=None)

	parser.add_argument('-z', '--input_SV_Calls_for_assembly', dest='input_SV_Calls_for_assembly', action='store_true', default=False)

	parser.add_argument('--workdir', default="./")
	parser.add_argument('--no_clean', action="store_true", default=False)

	args = parser.parse_args()
	tempDir = tempfile.mkdtemp(dir="./", prefix="pindel_work_")

	try:
		inputFastaFile, inputBamFile = indexBam(args.workdir, args.inputFastaFile, args.inputBamFile)
		# if args.insert_size==None:
		# 	meanInsertSize = int(getMeanInsertSize(args.inputBamFile, tempDir))
		# else:
		meanInsertSize=args.insert_size
		configFile = config(inputBamFile, meanInsertSize, args.sampleTag, tempDir)
		if args.procs == 1:
			cmd, pindelTempDir = pindel(inputFastaFile, configFile, args, tempDir)
			execute(cmd)
			cmd, pindelTmpVCF = pindel2vcf(inputFastaFile, args.sampleTag, pindelTempDir)
			execute(cmd)
			shutil.copy(pindelTmpVCF, args.outputVcfFile)
		else:
			seqs = get_bam_seq(inputBamFile)
			cmds = []
			runs = []
			for a in seqs:
				cmd, pindelTempDir = pindel(inputFastaFile, configFile, args, tempDir, a)
				cmds.append(cmd)
				runs.append(pindelTempDir)
			p = Pool(args.procs)
			values = p.map(execute, cmds, 1)
			cmds = []
			outs = []
			for a, b in zip(runs, seqs):
				cmd, out = pindel2vcf(inputFastaFile, args.sampleTag, a, chrome=b)
				cmds.append(cmd)
				outs.append(out)
			values = p.map(execute, cmds, 1)

			vcf_writer = None
			for file in outs:
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
