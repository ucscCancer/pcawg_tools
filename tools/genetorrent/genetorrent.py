#!/usr/bin/env python

import sys
import subprocess
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
	uuid = sys.argv[1]
	output = sys.argv[2]
	output_id = int(sys.argv[3])
	if len(sys.argv) > 4:
		cred_file = sys.argv[4]
	else:
		cred_file = "https://cghub.ucsc.edu/software/downloads/cghub_public.key"

	gtdownload = which("gtdownload")
	if gtdownload is None:
		sys.stderr.write("gtdownload not found\n")
		sys.exit(1)
	
	cred_dir = os.path.join(os.path.dirname(os.path.dirname(gtdownload)), "share", "GeneTorrent")
	serr = open("std.error", "w")
	
	#proc = subprocess.Popen( [gtdownload, "-c", cred_file, "-C", cred_dir, "-d", uuid, "-vv"], stderr=serr )
	proc = subprocess.Popen( [gtdownload, "-c", cred_file, "-d", uuid, "-vv"], stderr=serr )
	proc.communicate()
	serr.close()	

	copied = False
	for f in glob(os.path.join(uuid, "*.bam")):
		copied = True
		shutil.move(f, output)
	
	if not copied:
		handle = open("std.error")
		sys.stderr.write("Not found in: %s" % (",".join( glob(os.path.join(uuid, "*" ))) ))
		sys.stderr.write(handle.read())
		handle.close()
		sys.exit(1)

	json_file = open( 'galaxy.json', 'w' )
	info = dict( type = 'dataset',
				 ext = "bam",
				 name = uuid + " CGHub BAM",
				 uuid = uuid,
				dataset_id = output_id 
			)
	json_file.write( json.dumps( info ) )
	json_file.close()
