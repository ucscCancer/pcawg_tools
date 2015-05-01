#!/usr/bin/env python

import os
import re
import sys
import yaml
import subprocess
import tempfile
from StringIO import StringIO

def which(cmd):
    cmd = ["which",cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: return None
    return res

def parse_rg(line):
    res = re.search("^@RG\t(.*)$", line)
    text = res.group(1)
    fields = {}
    for f in text.split("\t"):
        res = re.search("^(\w+):(.*)$", f)
        fields[res.group(1)] = res.group(2)
    return fields

def write_rg(fields):
    rgm = []
    rgm.append( "ID:" + fields['ID'])
    for k,v in fields.items():
        if k != 'ID':
            rgm.append( "%s:%s" % (k,v) )
    return "@RG\t%s" % ( "\t".join( rgm ))


if __name__ == '__main__':
    config_file = sys.argv[1]
    in_bam      = sys.argv[2]
    out_bam     = sys.argv[3]

    #get original header
    sam_tools = which("samtools")
    if sam_tools is None:
        raise Exception("samtools not found")
    proc = subprocess.Popen( [sam_tools, "view", "-H", in_bam], stdout=subprocess.PIPE )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise Exception("Failed to read input BAM")

    with open(config_file) as handle:
        config = yaml.load(handle.read())

    header = StringIO()
    for line in stdout.split("\n"):
        if len(line):
            if 'RG' in config:
                if line.startswith("@RG"):
                    rg = parse_rg(line)
                    for k,v in config['RG'].items():
                        rg[k] = v
                    line = write_rg(rg)
            header.write(line + "\n")
    (t, tfile) = tempfile.mkstemp()
    os.close(t)
    with open(tfile, "w") as handle:
        handle.write(header.getvalue())

    proc = subprocess.Popen( "samtools reheader %s %s > %s" % (tfile, in_bam, out_bam), shell=True )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise Exception("Failed to read input BAM")
