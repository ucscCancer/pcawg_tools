#!/usr/bin/env python

import os
import re
import json
import csv
import argparse
import subprocess
import shutil
import synapseclient
import datetime
from nebula.target import Target, TargetFile
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

def run_scan(docstore, workdir, keyfile, upload_url, manifest):
    doc = FileDocStore(file_path=docstore)

    file_map = {
        'broad' : {},
        'muse' : {}
    }

    wl_map = {}
    with open(manifest) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            wl_map[row['Donor_ID']] = row

    for id, entry in doc.filter(visible=True):
        if entry.get('extension', None) in ["vcf", "vcf_bgzip"]:
            if 'tags' in entry:
                sample = None
                for s in entry['tags']:
                    tmp = s.split(":")
                    if tmp[0] == 'sample':
                        sample = tmp[1]

                pipeline = None
                method = None
                call_type = None
                variant_type = None
                if entry['name'] in ['MUSE_1.0rc', 'MUSE_0.9.9.5']:
                    pipeline = "muse"
                    method = entry['name'].replace(".", "-")
                    variant_type = 'somatic'
                    call_type = 'snv_mnv'
                elif entry['name'].split(".")[0] in ['broad-dRanger', 'broad-dRanger_snowman', 'broad-snowman', 'broad-mutect' ]:
                    pipeline = "broad"
                    method = entry['name'].split(".")[0]
                    if 'somatic' in entry['name']:
                        variant_type = 'somatic'
                    elif 'germline' in entry['name']:
                        variant_type = 'germline'
                    else:
                        raise Exception("Unknown variant type")
                    if 'snv_mnv.vcf' in entry['name']:
                        call_type = 'snv_mnv'
                    elif 'sv.vcf' in entry['name']:
                        call_type = 'sv'
                    elif 'indel.vcf' in entry['name']:
                        call_type = 'indel'
                    else:
                        raise Exception("Unknown call type: %s" % (entry['name']))
                else:
                    raise Exeception("Unknown pipeline %s" % (entry['name']))

                datestr = datetime.datetime.now().strftime("%Y%m%d")
                name = "%s.%s.%s.%s.%s" % (sample, method, datestr, variant_type, call_type )

                name = re.sub(r'.vcf$', '', name)
                if entry['extension'] == 'vcf':
                    file_name = name + ".vcf"
                elif entry['extension'] == 'vcf_bgzip':
                    file_name = name + ".vcf.gz"
                print file_name
                target = Target(uuid=entry['uuid'])
                if doc.size(target) > 0:
                    src_file = doc.get_filename(target)
                    dst_file = os.path.join(workdir, file_name)
                    shutil.copy(src_file, dst_file)

                    if entry['extension'] == 'vcf':
                        subprocess.check_call( "bgzip %s" % dst_file, shell=True )
                        dst_file = dst_file + ".gz"

                    subprocess.check_call("tabix -p vcf %s" % (dst_file), shell=True)
                    shutil.move("%s.tbi" % (dst_file), "%s.idx" % (dst_file))
                    subprocess.check_call("md5sum %s | awk '{print$1}' > %s.md5" % (dst_file, dst_file), shell=True)
                    subprocess.check_call("md5sum %s.idx | awk '{print$1}' > %s.idx.md5" % (dst_file, dst_file), shell=True)

                    if sample not in file_map[pipeline]:
                        file_map[pipeline][sample] = []

                    input_file = os.path.basename(dst_file)
                    file_map[pipeline][sample].append(input_file)

    for pipeline, samples in file_map.items():
        for sample, files in samples.items():
            with open( os.path.join(workdir, "%s.%s.sh" %(pipeline, sample)), "w" ) as handle:
                input_file = os.path.basename(dst_file)
                urls = [
                    "%scghub/metadata/analysisFull/%s" % (wl_map[sample]['Normal_GNOS_endpoint'], wl_map[sample]['Normal_Analysis_ID']),
                    "%scghub/metadata/analysisFull/%s" % (wl_map[sample]['Tumour_GNOS_endpoint'], wl_map[sample]['Tumour_Analysis_ID'])
                ]
                cmd_str = "perl /opt/vcf-uploader/gnos_upload_vcf.pl"
                cmd_str += " --metadata-urls %s" % (",".join(urls))
                cmd_str += " --vcfs %s " % (",".join(files))
                cmd_str += " --vcf-md5sum-files %s " % ((",".join( ("%s.md5" % i for i in files) )))
                cmd_str += " --vcf-idxs %s" % ((",".join( ("%s.idx" % i for i in files) )))
                cmd_str += " --vcf-idx-md5sum-files %s" % ((",".join( ("%s.idx.md5" % i for i in files) )))
                cmd_str += " --outdir %s.%s.dir" % (pipeline, sample)
                cmd_str += " --key %s " % (keyfile)
                cmd_str += " --upload-url %s" % (upload_url)
                cmd_str += " --study-refname-override tcga_pancancer_vcf_test"

                handle.write("""#!/bin/bash

%s

""" % (cmd_str) )


def run_errors(docstore):
    doc = FileDocStore(file_path=docstore)

    results = {}
    pending = {}
    for id, entry in doc.filter(visible=True):
        if entry.get('state', 'ok') in ['error']:
            print entry


def run_get(docstore, uuid, outpath):
    doc = FileDocStore(file_path=docstore)
    print doc.get_filename(Target(uuid=uuid))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("docstore", help="DocStore")
    subparsers = parser.add_subparsers(title="subcommand")

    parser_query = subparsers.add_parser('synapse')
    parser_query.set_defaults(func=run_synapse)
    parser_query.add_argument("--parent", required=True)
    parser_query.add_argument("--workdir", default="work")

    parser_audit = subparsers.add_parser('audit')
    parser_audit.add_argument("sample_list")
    parser_audit.set_defaults(func=run_audit)

    parser_scan = subparsers.add_parser('scan')
    parser_scan.set_defaults(func=run_scan)
    parser_scan.add_argument("--workdir", default="work")
    parser_scan.add_argument("--keyfile", default="cghub.pem")
    parser_scan.add_argument("--upload_url", default="https://gtrepo-osdc-tcga.annailabs.com")
    parser_scan.add_argument("--manifest")


    parser_audit = subparsers.add_parser('errors')
    parser_audit.set_defaults(func=run_errors)

    args = parser.parse_args()
    func = args.func

    vargs = vars(args)
    del vargs['func']

    func(**vargs)
