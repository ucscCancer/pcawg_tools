#!/usr/bin/python

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '-b', dest='bam', help='the bam file', required=True )
parser.add_argument ( '-t', dest='types', help='SV analysis type (DEL, DUP, INV, TRA)', nargs='+', required=True )
parser.add_argument ( '-q', dest='map_qual', help='min. paired-end mapping quality', default='0' )
parser.add_argument ( '-s', dest='mad_cutoff', help='insert size cutoff, median+s*MAD (deletions only)', default=5 )
parser.add_argument ( '-g', dest='genome', help='genome fasta file' )
parser.add_argument ( '-m', dest='min_flank', help='minimum flanking sequence size', default=13 )
parser.add_argument ( '-e', dest='epsilon', help='allowed epsilon deviation of PE vs. SR deletion', default=0.1 )
parser.add_argument ( '-o', dest='output', help='output file' )
parser.add_argument ( '--procs', dest='procs', type=int, default=1, help='Number of Parallel Processes' )


dellyPath="/delly/src/delly"

def execute( cmd, output="" ):
	tmp_dir = tempfile.mkdtemp()
	try:
		err = open(tmp_dir+"/errorLog", 'a')
		if output != "":
			out = open(output, 'w')
		else:
			out = subprocess.PIPE
		sys_env = dict(os.environ)
		sys_env['OMP_NUM_THREADS'] = "3"
		print "Running", cmd
		process = subprocess.Popen( args=shlex.split(cmd), stdout=out, stderr=err, env=sys_env )
		process.wait()
		err.close()
		if out != subprocess.PIPE:
			out.close()
	except Exception, e:
		sys.stderr.write("problem doing : %s\n" %(cmd))
		sys.stderr.write( '%s\n\n' % str(e) )


def delly(args, bamfile, tempDir):
	tempOutputs=[]
	cmds = []
	for typ in args.types:
		output=str(tempDir)+"/"+str(typ)
		tempOutputs.append(output)
		cmd = "%s %s -t %s -o %s -q %s -s %s -m %s -e %s " % (dellyPath, bamfile, typ, output, args.map_qual, args.mad_cutoff, args.min_flank, args.epsilon)
		if args.genome!=None and args.genome!="" and os.path.exists(args.genome):
			cmd += "-g %s" % (args.genome)
		cmds.append(cmd)
	p = Pool(args.procs)
	p.map(execute, cmds, 1)
	return tempOutputs


def merge(outputs, outputFile):
	template = vcf.Reader(filename=outputs[0])
	print "Merging to %s " % (outputFile)
	vcfWriter = vcf.Writer(open(outputFile, 'w'), template)
	for output in outputs:
		vcfReader = vcf.Reader(filename=output)
		for record in vcfReader:
			vcfWriter.write_record(record)
	return 0


def getVersion(program):
	try:
		tmp = tempfile.NamedTemporaryFile().name
		tmp_stdout = open( tmp, 'wb' )
		proc = subprocess.Popen( args=program, shell=True, stdout=tmp_stdout )
		tmp_stdout.close()
		returncode = proc.wait()
		stdout = None
		for line in open( tmp_stdout.name, 'rb' ):
			if line.lower().find( 'version' ) >= 0:
				stdout = line.strip()
				break
		if stdout:
			sys.stdout.write( '%s\n' % stdout )
	except:
		sys.stdout.write( 'Could not determine %s version\n' % (program) )


def __main__():
	print(os.path.dirname(os.path.realpath(__file__)))
	args = parser.parse_args()

	tempDir = tempfile.mkdtemp();
	getVersion(dellyPath)

	bamFile = os.path.join(tempDir, "input.bam")
	try:
		print "Indexing BAM"
		os.symlink(args.bam, bamFile)
		execute("samtools index " + bamFile)

	except Exception, e:
		print("problem while indexing bam file " + str(e))

	try:
		tempOutputs = delly(args, bamFile, tempDir)
	except Exception, e:
		sys.stdout.write("problem while runing delly " + str(e))

	try:
		if args.output:
			outputs=[]
			for output in tempOutputs:
				if os.path.exists(output):
					if os.path.getsize(output)>0:
						outputs.append(output)
			if len(outputs)>0:
				merge(outputs, args.output)
	except Exception, e:
		sys.stdout.write("problem while merging files " + str(e))

	finally:
		try:
			if os.path.exists(tempDir):
				shutil.rmtree(tempDir)
		except:
			pass

if __name__=="__main__":
	__main__()
