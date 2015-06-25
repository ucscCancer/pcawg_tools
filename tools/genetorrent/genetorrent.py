#!/usr/bin/env python

import sys
import subprocess
import argparse
import os
import shutil
from glob import glob
import json

def which(cmd):
	cmd = ["which",cmd]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	res = p.stdout.readline().rstrip()
	if len(res) == 0: return None
	return res

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("uuid")
	parser.add_argument("output")
	parser.add_argument("output_id")
	parser.add_argument("--cred-file", default="https://cghub.ucsc.edu/software/downloads/cghub_public.key")
	parser.add_argument("--server", default="cghub.ucsc.edu")
	args = parser.parse_args()
	
	gtdownload = which("gtdownload")
	if gtdownload is None:
		sys.stderr.write("gtdownload not found\n")
		sys.exit(1)
	
	cred_dir = os.path.join(os.path.dirname(os.path.dirname(gtdownload)), "share", "GeneTorrent")
	serr = open("std.error", "w")
	
	#proc = subprocess.Popen( [gtdownload, "-c", cred_file, "-C", cred_dir, "-d", uuid, "-vv"], stderr=serr )
	proc = subprocess.Popen( [gtdownload, "-c", args.cred_file, "-d", "https://%s/cghub/data/analysis/download/%s" % (args.server, args.uuid), "-vv"], stderr=serr )
	proc.communicate()
	serr.close()

	copied = False
	for f in glob(os.path.join(args.uuid, "*.bam")):
		copied = True
		shutil.move(f, args.output)
	
	if not copied:
		handle = open("std.error")
		sys.stderr.write("Not found in: %s" % (",".join( glob(os.path.join(args.uuid, "*" ))) ))
		sys.stderr.write(handle.read())
		handle.close()
		sys.exit(1)

	json_file = open( 'galaxy.json', 'w' )
	info = dict( type = 'dataset',
				 ext = "bam",
				 name = args.uuid + " CGHub BAM",
				 uuid = args.uuid,
				dataset_id = args.output_id 
			)
	json_file.write( json.dumps( info ) )
	json_file.close()
