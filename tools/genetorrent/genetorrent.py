#!/usr/bin/env python

import sys
import subprocess
import argparse
import os
import shutil
from glob import glob
import json
import random

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
	parser.add_argument("--retry", type=int, default=3)
	args = parser.parse_args()

	gtdownload = which("gtdownload")
	if gtdownload is None:
		sys.stderr.write("gtdownload not found\n")
		sys.exit(1)

	cred_dir = os.path.join(os.path.dirname(os.path.dirname(gtdownload)), "share", "GeneTorrent")
	serr = open("std.error", "w")
	
	server = random.choice(args.server.split(","))

	for i in range(args.retry):
		serr.write("Download Attempt %d\n" % (i))
		timeout = "10"
		if i != 0:
			#gtdownload doesn't factor in the file check time when it sets up timeout check
			#of it is possible to get a network timeout before ever opening a connection
			#the 120 minute timeout on retries is an attempt to get around that
			timeout = "120"
		#proc = subprocess.Popen( [gtdownload, "-c", cred_file, "-C", cred_dir, "-d", uuid, "-vv"], stderr=serr )
		cmd = [gtdownload, "-c", args.cred_file, "-k", timeout, "-d", "https://%s/cghub/data/analysis/download/%s" % (server, args.uuid), "-vv"]
		print "Running:", " ".join(cmd)
		proc = subprocess.Popen( cmd, stderr=serr )
		proc.communicate()
		#if the output directory is there and gtdownload finishied correctly
		if os.path.exists(args.uuid) and proc.returncode == 0:
			break
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
