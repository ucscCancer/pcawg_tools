#!/usr/bin/env python

import os
import re
import json
import argparse
import subprocess
import shutil
import synapseclient
from nebula.dag import Target, TargetFile
from nebula.docstore import FileDocStore

def run_synapse(docstore, parent, workdir):
    doc = FileDocStore(file_path=docstore)

    syn = synapseclient.Synapse()
    syn.login()

    for id, entry in doc.filter(visible=True, data_type='galaxy.datatypes.tabular.Vcf'):
        if 'tags' in entry:
            sample = None
            for s in entry['tags']:
                tmp = s.split(":")
                if tmp[0] == 'sample':
                    sample = tmp[1]
            name = entry['name']
            name = re.sub(r'.vcf$', '', name)
            file_name = sample + "." + name + ".snv_mnv.vcf"
            target = Target(uuid=entry['uuid'])
            if doc.size(target) > 0:
                src_file = doc.get_filename(target)
                dst_file = os.path.join(workdir, file_name)
                query = "select * from entity where parentId=='%s' and name=='%s'" % (parent, file_name + ".gz")
                r = syn.query(query)['results']
                if len(r) == 0:
                    #print r
                    print dst_file
                    shutil.copy(src_file, dst_file)
                    subprocess.check_call("bgzip %s" % (dst_file), shell=True)
                    f = synapseclient.File(dst_file + ".gz", parentId = parent, name=file_name + ".gz" )
                    f.fileType = 'vcf'
                    f.pipeline = 'UCSC'
                    f.variant_type = "snv"
                    f = syn.store(f,
                        executed="https://github.com/ucsccancer/pcawg_tools"
                    )
                else:
                    print "Skipping", file_name

def run_audit(docstore, sample_list):
    doc = FileDocStore(file_path=docstore)

    master_list = []
    with open(sample_list) as handle:
        for line in handle:
            master_list.append(line.rstrip())

    results = {}
    pending = {}
    for id, entry in doc.filter(visible=True, data_type='galaxy.datatypes.tabular.Vcf'):
        if 'tags' in entry:
            sample = None
            for s in entry['tags']:
                tmp = s.split(":")
                if tmp[0] == 'sample':
                    sample = tmp[1]
            if doc.size(entry) > 0:
                results[sample] = results.get(sample, []) + [entry['name']]
    for sample, files in results.items():
        print "%s (%s) %s" % (sample, len(files), "\t".join(files))

    for sample in master_list:
        if sample not in results or len(results[sample]) < 3:
            print "missing (%s)" % (len(results.get(sample, []))), sample



def run_get(docstore, uuid, outpath):
    doc = FileDocStore(file_path=docstore)
    print doc.get_filename(Target(uuid=uuid))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--docstore", help="DocStore", required=True)
    subparsers = parser.add_subparsers(title="subcommand")

    parser_query = subparsers.add_parser('synapse')
    parser_query.set_defaults(func=run_synapse)
    parser_query.add_argument("--parent", required=True)
    parser_query.add_argument("--workdir", default="work")

    parser_audit = subparsers.add_parser('audit')
    parser_audit.add_argument("sample_list")
    parser_audit.set_defaults(func=run_audit)


    args = parser.parse_args()
    func = args.func

    vargs = vars(args)
    del vargs['func']

    func(**vargs)
