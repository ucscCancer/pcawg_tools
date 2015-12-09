#!/usr/bin/env python

import re
import os
import json
import uuid
import argparse
import synqueue
import shutil
import synapseclient
import subprocess
import string
from math import isnan
from nebula.service import GalaxyService
from nebula.galaxy import GalaxyWorkflow
from nebula.docstore import from_url
from nebula.target import Target
from nebula.tasks import TaskGroup, GalaxyWorkflowTask
import tempfile
import datetime
from urlparse import urlparse
from urllib2 import urlopen

REFDATA_PROJECT="syn3241088"

config = {
  #"table_id" : "syn3498886",
  #"state_col" : "Processing State",
  #"primary_col" : "Donor_ID",
  #"table_id" : "syn4556289",
  "table_id" : "syn4951939",
  "primary_col" : "Submitter_donor_ID",
  "assignee_col" : "Assignee",
  "state_col" : "State"
}

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

key_map = {
    "cghub.ucsc.edu" : "/tool_data/files/cghub.key",
    "gtrepo-ebi.annailabs.com" : "/tool_data/files/icgc.key",
    "gtrepo-bsc.annailabs.com" : "/tool_data/files/icgc.key",
    "gtrepo-osdc-icgc.annailabs.com" : "/tool_data/files/icgc.key",
    "gtrepo-riken.annailabs.com" : "/tool_data/files/icgc.key",
    "gtrepo-dkfz.annailabs.com" : "/tool_data/files/icgc.key",
    "gtrepo-etri.annailabs.com" : "/tool_data/files/icgc.key",
}

def run_gen(args):
    args = parser.parse_args()

    syn = synapseclient.Synapse()
    syn.login()

    docstore = from_url(args.out_base)

    data_mapping = {
        "reference_genome" : "genome.fa",
        "dbsnp" : "dbsnp_132_b37.leftAligned.vcf",
        "cosmic" : "b37_cosmic_v54_120711.vcf",
        "gold_indels" : "Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf",
        "phase_one_indels" : "1000G_phase1.indels.hg19.sites.fixed.vcf",
        "centromere" : "centromere_hg19.bed"
    }

    if args.ref_download:
        #download reference files from Synapse and populate the document store
        for a in syn.chunkedQuery('select * from entity where parentId=="%s"' % (REFDATA_PROJECT)):
            print "found",  a['entity.name']
            if a['entity.name'] in data_mapping.values() or a['entity.name'].replace(".gz", "") in data_mapping.values():
                print "loading"
                ent = syn.get(a['entity.id'])
                id = ent.annotations['uuid'][0]
                t = Target(uuid=id)
                docstore.create(t)
                path = docstore.get_filename(t)
                name = ent.name
                if 'dataPrep' in ent.annotations:
                    if ent.annotations['dataPrep'][0] == 'gunzip':
                        subprocess.check_call("gunzip -c %s > %s" % (ent.path, path), shell=True)
                        name = name.replace(".gz", "")
                    else:
                        print "Unknown DataPrep"
                else:
                    shutil.copy(ent.path, path)
                docstore.update_from_file(t)
                meta = {}
                meta['name'] = name
                meta['uuid'] = id
                if 'dataPrep' in meta:
                    del meta['dataPrep']
                docstore.put(id, meta)

    dm = {}
    for k,v in data_mapping.items():
        hit = None
        for a in docstore.filter(name=v):
            hit = a[0]
        if hit is None:
            raise Exception("%s not found" % (v))
        dm[k] = { "uuid" : hit }

    workflow = GalaxyWorkflow(ga_file="workflows/Galaxy-Workflow-PCAWG_CGHUB.ga")
    tasks = TaskGroup()
    for ent in synqueue.listAssignments(syn, **config):
        #print "'%s'" % (ent['state']), ent['state'] == 'nan', type(ent['state']), type('nan')
        if not isinstance(ent['state'], basestring) and isnan(ent['state']):
            if "," not in ent['meta']['Tumour_WGS_alignment_GNOS_analysis_IDs']:
                normal_gnos_endpoints = ",".join(list(
                    urlparse(a).netloc for a in ent['meta']['Normal_WGS_alignment_GNOS_repos'].split('|')
                ))
                tumour_gnos_endpoints = ",".join(list(
                    urlparse(a).netloc for a in ent['meta']['Tumour_WGS_alignment_GNOS_repos'].split('|')
                ))

                gnos_endpoint = urlparse(ent['meta']['Normal_WGS_alignment_GNOS_repos']).netloc
                task = GalaxyWorkflowTask("workflow_%s" % (ent['id']),
                    workflow,
                    inputs=dm,
                    parameters={
                        'normal_bam_download' : {
                            "uuid" : ent['meta']['Normal_WGS_alignment_GNOS_analysis_ID'],
                            "gnos_endpoint" : normal_gnos_endpoints,
                            "cred_file" : key_map[normal_gnos_endpoints.split(",")[0]]
                        },
                        'tumor_bam_download' : {
                            "uuid" : ent['meta']['Tumour_WGS_alignment_GNOS_analysis_IDs'],
                            "gnos_endpoint" : tumour_gnos_endpoints,
                            "cred_file" : key_map[tumour_gnos_endpoints.split(",")[0]]
                        },
                        'broad_variant_pipeline' : {
                            "broad_ref_dir" : "/tool_data/files/refdata",
                            "sample_id" : ent['meta']['Submitter_donor_ID']
                        }
                    },
                    tags=[ "donor:%s" % (ent['meta']['Submitter_donor_ID']) ]
                )
                tasks.append(task)

    if not os.path.exists("%s.tasks" % (args.out_base)):
        os.mkdir("%s.tasks" % (args.out_base))

    for data in tasks:
        with open("%s.tasks/%s" % (args.out_base, data.task_id), "w") as handle:
            handle.write(json.dumps(data.to_dict()))
            state_file = "%s.tasks/%s.state" % (args.out_base, data.task_id)
            if os.path.exists( state_file ):
                os.unlink( state_file )

    print "Tasks Created: %s" % (len(tasks))

    if args.create_service:
        service = GalaxyService(
            docstore=docstore,
            galaxy=args.galaxy,
            sudo=args.sudo,
            timeout=args.timeout,
            tool_data=os.path.abspath("tool_data"),
            tool_dir=os.path.abspath("tools"),
            work_dir=args.work_dir,
            smp=[
                ["MuSE", 8],
                ["pindel", 8],
                ["muTect", 8],
                ["delly", 4],
                ["gatk_bqsr", 12],
                ["gatk_indel", 24],
                ["bwa_mem", 12],
                ["broad_variant_pipline", 24]
            ]
        )
        with open("%s.service" % (args.out_base), "w") as handle:
            s = service.get_config()
            if args.scratch:
                print "Using scratch", args.scratch
                s.set_docstore_config(cache_path=args.scratch, open_perms=True)
            s.store(handle)

def run_audit(args):
    args = parser.parse_args()

    syn = synapseclient.Synapse()
    syn.login()

    docstore = from_url(args.out_base)

    donor_map = {}
    for id, ent in docstore.filter( state="ok" ):
        if ent['visible']:
            if docstore.size(Target(id)) > 0:
                donor = None
                for i in ent['tags']:
                    t = i.split(":")
                    if t[0] == "donor":
                        donor = t[1]
                if donor not in donor_map:
                    donor_map[donor] = {}
                
                if 'update_time' in ent:
                    dt = datetime.datetime.strptime(ent['update_time'], "%Y-%m-%dT%H:%M:%S.%f")
                else:
                    dt = None
                donor_map[donor][ent['name']] = { 'id' : id, 'update_time' : dt }
                
                

    for ent in synqueue.listAssignments(syn, list_all=True, **config):
        if ent['meta']['Submitter_donor_ID'] in donor_map:
            print "%s\t%s\t%s" % (
                ent['meta']['Submitter_donor_ID'], 
                len(donor_map[ent['meta']['Submitter_donor_ID']]), 
                max(v['update_time'] for k,v in donor_map[ent['meta']['Submitter_donor_ID']].items() if v['update_time'] is not None )
            )

def run_gnos_audit(args):
    args = parser.parse_args()
    syn = synapseclient.Synapse()
    syn.login()

    """
    r_map = synqueue.getValues(syn, "Normal_WGS_alignment_GNOS_repos", **config)
    repos = {}
    for i in r_map.values():
        repos[i] = True
    print repos.keys()
    """
    server_list = ["https://gtrepo-osdc-tcga.annailabs.com", "https://gtrepo-osdc-icgc.annailabs.com"]

    uuid_map = {}
    uuid_map['broad'] = dict( (a[1], a[0]) for a in synqueue.getValues(syn, "Broad_VCF_UUID",  **config).items() )
    uuid_map['muse']  = dict( (a[1], a[0]) for a in synqueue.getValues(syn, "Muse_VCF_UUID",  **config).items() )
    uuid_map['broad_tar'] = dict( (a[1], a[0]) for a in synqueue.getValues(syn, "Broad_TAR_UUID", **config).items() )

    analysis_re = re.compile(r'<analysis_id>(.*)</analysis_id>')

    found = {}
    for server in server_list:
        for study in ['tcga_pancancer_vcf', 'icgc_pancancer_vcf']:
            handle = urlopen(server + "/cghub/metadata/analysisId?study=%s" % (study))
            for line in handle:
                res = analysis_re.search(line)
                if res:
                    gid = res.group(1)
                    for p in ['broad', 'muse','broad_tar']:
                        if gid in uuid_map[p]:
                            pid = uuid_map[p][gid]
                            if pid not in found:
                                found[pid] = {}
                            found[pid][p] = [server, gid]

    print "\n\n"
    for p in found:
        print p, "\t".join( "%s\t%s" % (found[p].get(a, ["NA","NA"])[0], found[p].get(a, ["NA","NA"])[1]) for a in ['broad', 'muse','broad_tar'])


def run_register(args):

    syn = synapseclient.Synapse()
    syn.login()
    synqueue.registerAssignments(syn, args.count, display=True, force=args.force, **config)

upload_remap = {
    "cghub.ucsc.edu" : "gtrepo-osdc-tcga.annailabs.com"
}

upload_study_map = {
    "cghub.ucsc.edu" : "tcga_pancancer_vcf",
    "gtrepo-osdc-tcga.annailabs.com" : "tcga_pancancer_vcf",
    "gtrepo-ebi.annailabs.com" : "icgc_pancancer_vcf",
    "gtrepo-bsc.annailabs.com" : "icgc_pancancer_vcf",
    "gtrepo-osdc-icgc.annailabs.com" : "icgc_pancancer_vcf",
    "gtrepo-riken.annailabs.com" : "icgc_pancancer_vcf",
    "gtrepo-dkfz.annailabs.com" : "icgc_pancancer_vcf",
    "gtrepo-etri.annailabs.com" : "icgc_pancancer_vcf"
}

upload_key_map = {
    "gtrepo-osdc-tcga.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/bionimbus.key"),
    "gtrepo-ebi.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/icgc.key"),
    "gtrepo-bsc.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/icgc.key"),
    "gtrepo-osdc-icgc.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/icgc.key"),
    "gtrepo-riken.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/icgc.key"),
    "gtrepo-dkfz.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/icgc.key"),
    "gtrepo-etri.annailabs.com" : os.path.join(BASEDIR, "tool_data/files/icgc.key")
}

def run_uploadprep(args):

    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)
    doc = from_url(args.out_base)
    file_map = {
        'broad' : {},
        'muse' : {},
        'broad_tar' : {}
    }

    syn = synapseclient.Synapse()
    syn.login()

    wl_map = {}
    job_map = {}
    for ent in synqueue.listAssignments(syn, list_all=True, **config):
        if len(args.ids) == 0 or ent['id'] in args.ids:
            wl_map[ent['id']] = ent['meta']

    uuid_map = {}
    uuid_map['broad'] = synqueue.getValues(syn, "Broad_VCF_UUID", orSet=lambda x: str(uuid.uuid4()), **config)
    uuid_map['muse']  = synqueue.getValues(syn, "Muse_VCF_UUID", orSet=lambda x: str(uuid.uuid4()), **config)
    uuid_map['broad_tar'] = synqueue.getValues(syn, "Broad_TAR_UUID", orSet=lambda x: str(uuid.uuid4()), **config)

    #scan through all of the docs
    print "Scanning"
    for id, entry in doc.filter():
        donor = None
        #look for docs with donor tags
        if 'tags' in entry and 'state' in entry and entry['state'] == 'ok':
            for s in entry['tags']:
                tmp = s.split(":")
                if tmp[0] == 'donor':
                    donor = tmp[1]
        if donor is not None and donor in wl_map:
            if donor not in job_map:
                job_map[donor] = {}
                
            upload_host = urlparse(wl_map[donor]['Normal_WGS_alignment_GNOS_repos']).netloc
            for k,v in upload_remap.items():
                upload_host = re.sub(k,v,upload_host)
            if args.server_filter is not None and upload_host != args.server_filter:
                print "skipping", upload_host,  args.server_filter
                continue
            
            #scan out the job metrics for this job
            if 'job' in entry and 'job_metrics' in entry['job']:
                job_id = entry['job']['id']
                tool_id = entry['job']['tool_id']
                job_info = { tool_id : {} }
                for met in entry['job']['job_metrics']:
                    job_info[tool_id][met['name']] = met['raw_value']
                job_map[donor][job_id] = job_info
            donor_tumor = wl_map[donor]['Tumour_WGS_aliquot_IDs']

            #look for the vcf output files
            if entry.get('visible', False) and entry.get('extension', None) in ["vcf", "vcf_bgzip"]:
                pipeline = None
                method = None
                call_type = None
                variant_type = None
                #fill out the info depending on which caller created the file
                if entry['name'].split('.')[0] in ['MUSE_1']:
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
                    raise Exception("Unknown pipeline %s" % (entry['name']))

                datestr = datetime.datetime.now().strftime("%Y%m%d")
                name = "%s.%s.%s.%s.%s" % (donor_tumor, method, datestr, variant_type, call_type )

                name = re.sub(r'.vcf$', '', name)
                if entry['extension'] == 'vcf':
                    file_name = name + ".vcf"
                elif entry['extension'] == 'vcf_bgzip':
                    file_name = name + ".vcf.gz"
                target = Target(uuid=entry['uuid'])
                if doc.size(target) > 0:
                    src_file = doc.get_filename(target)
                    dst_dir = os.path.join(args.workdir, upload_host, donor, uuid_map[pipeline][donor])
                    if not os.path.exists(dst_dir):
                        os.makedirs(dst_dir)
                    dst_file = os.path.join(dst_dir, file_name)

                    shutil.copy(src_file, dst_file)
                    #if the files wasn't compressed already, go ahead and do that
                    if entry['extension'] == 'vcf':
                        subprocess.check_call( "bgzip -c %s > %s.gz" % (src_file, dst_file), shell=True )
                        dst_file = dst_file + ".gz"

                    #add file to output map
                    if donor not in file_map[pipeline]:
                        file_map[pipeline][donor] = []
                    input_file = os.path.basename(dst_file)
                    file_map[pipeline][donor].append(input_file)
            else:
                if entry['name'] == "broad.tar.gz":
                    pipeline = "broad_tar"
                    target = Target(uuid=entry['uuid'])
                    src_file = doc.get_filename(target)
                    file_map[pipeline][donor] = [ src_file ]
                    dst_dir = os.path.join(args.workdir, upload_host, donor, uuid_map[pipeline][donor])
                    if not os.path.exists(dst_dir):
                        os.makedirs(dst_dir)


    timing_map = {}
    for donor in job_map:
        timing_map[donor] = {}
        for job_id in job_map[donor]:
            for tool_id in job_map[donor][job_id]:
                if tool_id not in timing_map[donor]:
                    timing_map[donor][tool_id] = []
                timing_map[donor][tool_id].append( job_map[donor][job_id][tool_id] )

    result_counts = {}
    for pipeline, donors in file_map.items():
        for donor in donors:
            result_counts[donor] = result_counts.get(donor, 0) + 1

    #go through every pipeline
    print "Building Wrappers"
    for pipeline, donors in file_map.items():
        print " ...", pipeline
        #for that pipeline go through every donor
        for donor, files in donors.items():
            #we're only outputing data for donors on the work list
            if donor in wl_map and result_counts[donor] == 3:
                upload_host = urlparse(wl_map[donor]['Normal_WGS_alignment_GNOS_repos']).netloc
                for k,v in upload_remap.items():
                    upload_host = re.sub(k,v,upload_host)
                if args.server_filter is not None and upload_host != args.server_filter:
                    print "skipping", upload_host,  args.server_filter
                    continue

                #output the timing json
                timing_json = os.path.abspath(os.path.join(args.workdir, upload_host, donor, uuid_map[pipeline][donor], "timing.json" ))
                with open( timing_json, "w" ) as handle:
                    handle.write(json.dumps( timing_map[donor] ) )

                #output meta-data file
                with open( os.path.join(args.workdir, upload_host, donor, uuid_map[pipeline][donor], "meta.json"), "w" ) as handle:
                    urls = [
                        "%scghub/metadata/analysisFull/%s" % (wl_map[donor]['Normal_WGS_alignment_GNOS_repos'].split("|")[0], wl_map[donor]['Normal_WGS_alignment_GNOS_analysis_ID']),
                        "%scghub/metadata/analysisFull/%s" % (wl_map[donor]['Tumour_WGS_alignment_GNOS_repos'].split("|")[0], wl_map[donor]['Tumour_WGS_alignment_GNOS_analysis_IDs'])
                    ]
                    related_uuids = []
                    for p in uuid_map:
                        if p != pipeline:
                            related_uuids.append(uuid_map[p][donor])
                    meta = {
                        "donor" : donor,
                        "upload_host" : upload_host,
                        "uuid" : uuid_map[pipeline][donor],
                        "pipeline" : pipeline,
                        "related_uuids" : related_uuids,
                        "metadata-urls" : urls,
                        "timing" : timing_map[donor]
                    }
                    handle.write(json.dumps(meta))

                #output the prep script
                with open( os.path.join(args.workdir, upload_host, donor, uuid_map[pipeline][donor], "prep.sh"), "w" ) as handle:
                    input_file = os.path.basename(dst_file)
                    donor_tumor = wl_map[donor]['Tumour_WGS_aliquot_IDs']
                    upload_key = upload_key_map[upload_host]

                    if pipeline in ['broad', 'muse']:
                        prep_cmd_str = ""
                        for vcf in files:
                            prep_cmd_str += "tabix -p vcf %s\n" % (vcf)
                            prep_cmd_str += "mv %s.tbi %s.idx\n" % (vcf,vcf)
                            prep_cmd_str += "md5sum %s | awk '{print$1}' > %s.md5\n" % (vcf, vcf)
                            prep_cmd_str += "md5sum %s.idx | awk '{print$1}' > %s.idx.md5\n\n" % (vcf, vcf)

                    if pipeline in ['broad_tar']:
                        prep_cmd_str = ""
                        new_files = []
                        for tar in files:
                            basename = donor_tumor + ".broad.intermediate"
                            prep_cmd_str = "%s/remap_broad_tar.py %s %s %s --rename %s %s" % (
                                os.path.dirname(os.path.abspath(__file__)),
                                tar,
                                "./",
                                basename,
                                donor, donor_tumor
                            )
                            new_files.append( basename + ".tar" )
                    handle.write("""#!/bin/bash
cd `dirname $0`
set -ex
%s
echo $? > `basename $0`.ready
""" % (prep_cmd_str))

                #output the upload script
                with open( os.path.join(args.workdir, upload_host, donor, uuid_map[pipeline][donor], "upload.sh"), "w" ) as handle:
                    handle.write("""#!/bin/bash
cd `dirname $0`
set -ex
""")
                    file_rename_cmd = ''
                    file_rename_map = {}
                    update_analysis_xml = ''
                    broad_docker_version_string = "-10"
                    # Handle the Broad output files
                    if pipeline in ['broad']:
                        # First, we need to add "SB-10" to all file names, right before that date-part (after dRanger., mutect., snowman.)
                        file_rename_cmd = ''
                        file_rename_map = {}
                        for f in files:
                            parts = re.split('(broad-dRanger_snowman|broad-dRanger|broad-snowman|broad-mutect)',f)
                            new_file_name = parts[0]+parts[1]+broad_docker_version_string+parts[2]
                            file_rename_map[f]=new_file_name
                            
                        for k in file_rename_map:
                            file_rename_cmd += 'mv '+k+'\t\t'+file_rename_map[k]+'\n'
                            # We also have to rename the broad files vcf.gz.idx, vcf.gz.md5, vcf.gz.idx.md5
                            # But since these files aren't in the list `files`, we can't rely on that list to properly
                            # generate the renaming commands, so we'll do it here.
                            parts = re.split('(broad-dRanger_snowman|broad-dRanger|broad-snowman|broad-mutect)',k)
                            end_without_suffix = re.split('(\.vcf)',parts[2])
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf.gz.md5' + '\t\t' + parts[0] + parts[1] + broad_docker_version_string + end_without_suffix[0] + '.vcf.gz.md5\n'
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf.gz.idx' + '\t\t' + parts[0] + parts[1] + broad_docker_version_string + end_without_suffix[0] + '.vcf.gz.idx\n'
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf.gz.idx.md5' + '\t\t' + parts[0] + parts[1] + broad_docker_version_string + end_without_suffix[0] + '.vcf.gz.idx.md5\n'
                        
                        submit_cmd_str = "perl -I /opt/gt-download-upload-wrapper/gt-download-upload-wrapper-2.0.12/lib"
                        submit_cmd_str += " /opt/vcf-uploader/vcf-uploader-2.0.6/gnos_upload_vcf.pl"
                        submit_cmd_str += " --metadata-urls %s" % (",".join(urls))
                        submit_cmd_str += " --vcfs %s " % (",".join(file_rename_map.values()))
                        submit_cmd_str += " --vcf-md5sum-files %s " % ((",".join( ("%s.md5" % i for i in file_rename_map.values()) )))
                        submit_cmd_str += " --vcf-idxs %s" % ((",".join( ("%s.idx" % i for i in file_rename_map.values()) )))
                        submit_cmd_str += " --vcf-idx-md5sum-files %s" % ((",".join( ("%s.idx.md5" % i for i in file_rename_map.values()) )))
                        submit_cmd_str += " --outdir %s.%s.dir" % (pipeline, donor_tumor)
                        submit_cmd_str += " --key %s " % (upload_key)
                        submit_cmd_str += " --k-timeout-min 10"
                        submit_cmd_str += " --upload-url https://%s" % (upload_host)
                        submit_cmd_str += " --study-refname-override %s" % (upload_study_map[upload_host])
                        submit_cmd_str += " --workflow-url '%s'" % args.pipeline_src
                        submit_cmd_str += " --workflow-src-url '%s'" % args.pipeline_src
                        submit_cmd_str += " --workflow-name '%s'" % args.pipeline_name
                        submit_cmd_str += " --workflow-version '%s'" % args.pipeline_version
                        submit_cmd_str += " --vm-instance-type '%s'" % args.vm_instance_type
                        submit_cmd_str += " --vm-instance-cores %s" % args.vm_instance_cores
                        submit_cmd_str += " --vm-instance-mem-gb %s" % args.vm_instance_mem_gb
                        submit_cmd_str += " --vm-location-code %s" % args.vm_location_code
                        submit_cmd_str += " --timing-metrics-json %s" % (timing_json)
                        submit_cmd_str += " --workflow-file-subset %s" % (pipeline)
                        submit_cmd_str += " --related-file-subset-uuids %s" % (",".join(related_uuids))
                        submit_cmd_str += " --uuid %s" % (uuid_map[pipeline][donor])
                        if args.rsync:
                            submit_cmd_str += " --skip-upload --skip-validate"
                    
                    # Handle the MuSE output files
                    if pipeline in ['muse']:
                        # First, we need to add the MuSE v1.0rc_submission_b391201 SHA1 hash sum to the file names.
                        # ATTENTION: If you change the MuSE version, you will need to put the correct hashsum in here. 
                        muse_1_0rc_identifier='b391201'
                        file_rename_cmd = ''
                        file_rename_map = {}
                        for f in files:
                            parts = re.split('(MUSE_1-0rc)',f)
                            new_file_name = parts[0] + parts[1] + '-' + muse_1_0rc_identifier + parts[2]
                            file_rename_map[f] = new_file_name
                            
                        for k in file_rename_map:
                            file_rename_cmd += 'mv ' + k + '\t\t' + file_rename_map[k] + '\n'
                            # We also have to rename the MUSE vcf, vcf.gz.idx, vcf.gz.md5, vcf.gz.idx.md5
                            # But since these files aren't in the list `files`, we can't rely on that list to properly
                            # generate the renaming commands, so we'll do it here.
                            parts = re.split('(MUSE_1-0rc)',k)
                            end_without_suffix = re.split('(\.vcf)',parts[2])
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf' + '\t\t' + parts[0] + parts[1] + muse_1_0rc_identifier + end_without_suffix[0] + '.vcf\n'                            
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf.gz.md5' + '\t\t' + parts[0] + parts[1] + muse_1_0rc_identifier + end_without_suffix[0] + '.vcf.gz.md5\n'
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf.gz.idx' + '\t\t' + parts[0] + parts[1] + muse_1_0rc_identifier + end_without_suffix[0] + '.vcf.gz.idx\n'
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.vcf.gz.idx.md5' + '\t\t' + parts[0] + parts[1] + muse_1_0rc_identifier + end_without_suffix[0] + '.vcf.gz.idx.md5\n'
                            
                        submit_cmd_str = "perl -I /opt/gt-download-upload-wrapper/gt-download-upload-wrapper-2.0.12/lib"
                        submit_cmd_str += " /opt/vcf-uploader/vcf-uploader-2.0.6/gnos_upload_vcf.pl"
                        submit_cmd_str += " --metadata-urls %s" % (",".join(urls))
                        submit_cmd_str += " --vcfs %s " % (",".join(file_rename_map.values()))
                        submit_cmd_str += " --vcf-md5sum-files %s " % ((",".join( ("%s.md5" % i for i in  file_rename_map.values()) )))
                        submit_cmd_str += " --vcf-idxs %s" % ((",".join( ("%s.idx" % i for i in  file_rename_map.values()) )))
                        submit_cmd_str += " --vcf-idx-md5sum-files %s" % ((",".join( ("%s.idx.md5" % i for i in  file_rename_map.values()) )))
                        submit_cmd_str += " --outdir %s.%s.dir" % (pipeline, donor_tumor)
                        submit_cmd_str += " --key %s " % (upload_key)
                        submit_cmd_str += " --k-timeout-min 10"
                        submit_cmd_str += " --upload-url https://%s" % (upload_host)
                        submit_cmd_str += " --study-refname-override %s" % (upload_study_map[upload_host])
                        submit_cmd_str += " --workflow-url '%s'" % args.pipeline_src
                        submit_cmd_str += " --workflow-src-url '%s'" % args.pipeline_src
                        submit_cmd_str += " --workflow-name '%s'" % args.pipeline_name
                        submit_cmd_str += " --workflow-version '%s'" % args.pipeline_version
                        submit_cmd_str += " --vm-instance-type '%s'" % args.vm_instance_type
                        submit_cmd_str += " --vm-instance-cores %s" % args.vm_instance_cores
                        submit_cmd_str += " --vm-instance-mem-gb %s" % args.vm_instance_mem_gb
                        submit_cmd_str += " --vm-location-code %s" % args.vm_location_code
                        submit_cmd_str += " --timing-metrics-json %s" % (timing_json)
                        submit_cmd_str += " --workflow-file-subset %s" % (pipeline)
                        submit_cmd_str += " --related-file-subset-uuids %s" % (",".join(related_uuids))
                        submit_cmd_str += " --uuid %s" % (uuid_map[pipeline][donor])
                        if args.rsync:
                            submit_cmd_str += " --skip-upload --skip-validate"

                    # Handle the Broad intermediate tar file
                    if pipeline in ['broad_tar']:
                        # First, we need to add "SB-10" to all file names, right before that date-part (after dRanger., mutect., snowman.)
                        file_rename_cmd = ''
                        file_rename_map = {}
                        for f in new_files:
                            parts = re.split('(broad)',f)
                            new_file_name = parts[0] + parts[1] + broad_docker_version_string + parts[2]
                            file_rename_map[f] = new_file_name
                            
                        for k in file_rename_map:
                            file_rename_cmd += 'mv ' + k + '\t\t' + file_rename_map[k] + '\n'
                            # We also have to rename the broad intermediate tar md5 file.
                            # But since these files aren't in the list `new_files`, we can't rely on that list to properly
                            # generate the renaming commands, so we'll do it here.
                            parts = re.split('(broad)',k)
                            end_without_suffix =re.split('(\.tar)', parts[2])
                            file_rename_cmd += 'mv ' + parts[0] + parts[1] + end_without_suffix[0] + '.tar' + '\t\t' + parts[0] + parts[1] + broad_docker_version_string + end_without_suffix[0] + '.tar.md5\n'
                            
                        submit_cmd_str = "perl -I /opt/gt-download-upload-wrapper/gt-download-upload-wrapper-2.0.12/lib"
                        submit_cmd_str += " /opt/vcf-uploader/vcf-uploader-2.0.6/gnos_upload_vcf.pl"
                        submit_cmd_str += " --metadata-urls %s" % (",".join(urls))
                        submit_cmd_str += " --tarballs %s " % (",".join(file_rename_map.values()))
                        submit_cmd_str += " --tarball-md5sum-files %s " % ((",".join( ("%s.md5" % i for i in file_rename_map.values()) )))
                        submit_cmd_str += " --outdir %s.%s.dir" % (pipeline, donor_tumor)
                        submit_cmd_str += " --key %s " % (upload_key)
                        submit_cmd_str += " --k-timeout-min 10"
                        submit_cmd_str += " --upload-url https://%s" % (upload_host)
                        submit_cmd_str += " --study-refname-override %s" % (upload_study_map[upload_host])
                        submit_cmd_str += " --workflow-url '%s'" % args.pipeline_src
                        submit_cmd_str += " --workflow-src-url '%s'" % args.pipeline_src
                        submit_cmd_str += " --workflow-name '%s'" % args.pipeline_name
                        submit_cmd_str += " --workflow-version '%s'" % args.pipeline_version
                        submit_cmd_str += " --vm-instance-type '%s'" % args.vm_instance_type
                        submit_cmd_str += " --vm-instance-cores %s" % args.vm_instance_cores
                        submit_cmd_str += " --vm-instance-mem-gb %s" % args.vm_instance_mem_gb
                        submit_cmd_str += " --workflow-file-subset %s" % (pipeline)
                        submit_cmd_str += " --timing-metrics-json %s" % (timing_json)
                        submit_cmd_str += " --related-file-subset-uuids %s" % (",".join(related_uuids))
                        submit_cmd_str += " --uuid %s" % (uuid_map[pipeline][donor])
                        if args.rsync:
                            submit_cmd_str += " --skip-upload --skip-validate"

                    if args.rsync:
                        handle.write("""
{rename_cmd}

{update_analysis_xml}

{prep_cmd} --out ./
mv vcf/*/*.xml ./
rm -rf vcf xml2
ssh -i {key} -o StrictHostKeyChecking=no {remote_host} mkdir -p {remote_dir}/{upload_host}/{donor}/
RSYNC_RSH="ssh -i {key} -o StrictHostKeyChecking=no" {script_dir}/rsync_wrap -av {local_dir} {rsync}/{upload_host}/{donor}/
echo $? > `basename $0`.submitted
""".format(
    rename_cmd=file_rename_cmd,
    update_analysis_xml=update_analysis_xml,
    prep_cmd=submit_cmd_str, 
    key=os.path.abspath(args.rsync_key),
    rsync=args.rsync,
    donor=donor,
    pipeline=pipeline,
    upload_host=upload_host,
    remote_host=args.rsync.split(":")[0],
    remote_dir=args.rsync.split(":")[1],
    script_dir=os.path.dirname(os.path.abspath(__file__)),
    local_dir= os.path.join(os.path.abspath(args.workdir), upload_host, donor, uuid_map[pipeline][donor])
    ))
                    else:
                        handle.write("""
%s
echo $? > `basename $0`.submitted
""" % (submit_cmd_str))

def run_list(args):
    syn = synapseclient.Synapse()
    syn.login()
    synqueue.listAssignments(syn, display=True, **config)

def run_set(args):
    syn = synapseclient.Synapse()
    syn.login()

    synqueue.setStates(syn, args.state, args.ids, **config)

def check_within(datestr, max_hours):
    delta = datetime.datetime.now() - datetime.datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S.%f")
    age =  delta.days * 24 + delta.seconds / 3600
    if age < max_hours:
        return True
    return False

def run_errors(args):

    doc = from_url(args.out_base)

    for id, entry in doc.filter():
        if entry.get('state', '') == 'error':
            if args.within is None or 'update_time' not in entry or check_within(entry['update_time'], args.within):
                print "Dataset", id, entry.get("job", {}).get("tool_id", ""), entry.get('update_time', ''), entry.get("tags", "")
                if args.full:
                    if 'provenance' in entry:
                        print "tool:", entry['provenance']['tool_id']
                        print "-=-=-=-=-=-=-"
                    print entry['job']['stdout']
                    print "-------------"
                    print entry['job']['stderr']
                    print "-=-=-=-=-=-=-"

def run_timing(args):
    doc = from_url(args.out_base)

    job_map = {}
    for id, entry in doc.filter(state='ok'):
        donor = None
        #look for docs with donor tags
        if 'tags' in entry and 'state' in entry and entry['state'] == 'ok':
            for s in entry['tags']:
                tmp = s.split(":")
                if tmp[0] == 'donor':
                    donor = tmp[1]
        if donor is not None:   
            if 'job' in entry and 'job_metrics' in entry['job']:
                job_id = entry['job']['id']
                tool_id = entry['job']['tool_id']
                job_info = { }
                for met in entry['job']['job_metrics']:
                    job_info[met['name']] = met['raw_value']
                if tool_id not in job_map:
                    job_map[tool_id] = {}
                job_map[tool_id][donor] = job_info
    
    for tool in job_map:
        for donor in job_map[tool]:
            if 'runtime_seconds' in job_map[tool][donor]:
                print "%s\t%s\t%s" % (tool, donor, job_map[tool][donor]['runtime_seconds'])
            
def run_clean(args):
    doc = from_url(args.out_base)
    error_count = 0
    for id, entry in doc.filter(state='error'):
        if entry.get('state', '') == 'error':
            doc.delete(Target(uuid=id))
            error_count += 1
    print "Error count", error_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers( title='commands')

    parser_gen = subparsers.add_parser('gen')
    parser_gen.add_argument("--out-base", default="pcawg_data")
    parser_gen.add_argument("--ref-download", action="store_true", default=False)
    parser_gen.add_argument("--create-service", action="store_true", default=False)
    parser_gen.add_argument("--scratch", default=None)
    parser_gen.add_argument("--work-dir", default=None)
    parser_gen.add_argument("--sudo", action="store_true", default=False)
    parser_gen.add_argument("--galaxy", default="bgruening/galaxy-stable")
    parser_gen.add_argument("--timeout", type=int, default=60)
    parser_gen.add_argument("--tool-data", default=os.path.abspath("tool_data"))
    parser_gen.add_argument("--tool-dir", default=os.path.abspath("tools"))

    parser_gen.set_defaults(func=run_gen)

    parser_audit = subparsers.add_parser('audit')
    parser_audit.add_argument("--out-base", default="pcawg_data")
    parser_audit.set_defaults(func=run_audit)

    parser_gnos = subparsers.add_parser('gnos-audit')
    parser_gnos.add_argument("--study", default="tcga_pancancer_vcf")
    parser_gnos.set_defaults(func=run_gnos_audit)

    parser_register = subparsers.add_parser('register',
                                       help='Returns set of new assignments')
    parser_register.add_argument("-c", "--count", help="Number of assignments to get", type=int, default=1)
    parser_register.add_argument("-f", "--force", help="Force Register specific ID", default=None)
    parser_register.set_defaults(func=run_register)

    parser_upload = subparsers.add_parser('upload-prep')
    parser_upload.add_argument("--workdir", default="upload")
    #parser_upload.add_argument("--keyfile", default="/keys/cghub.key")
    parser_upload.add_argument("--out-base", default="pcawg_data")
    parser_upload.add_argument("--pipeline-src", default="https://github.com/ucscCancer/pcawg_tools")
    parser_upload.add_argument("--pipeline-version", default="1.1.0")
    parser_upload.add_argument("--pipeline-name", default="BROAD_MUSE_PIPELINE")
    parser_upload.add_argument("--server-filter", default=None)
    parser_upload.add_argument("--study", default="tcga_pancancer_vcf")
    parser_upload.add_argument("--vm-instance-type", default="d2.8xlarge")
    parser_upload.add_argument("--vm-instance-cores", default="32")
    parser_upload.add_argument("--vm-instance-mem-gb", default="256")
    parser_upload.add_argument("--vm-location-code", default="27")
    parser_upload.add_argument("--rsync", default=None)
    parser_upload.add_argument("--rsync-key", default=None)

    parser_upload.add_argument("ids", nargs="*")
    parser_upload.set_defaults(func=run_uploadprep)

    parser_list = subparsers.add_parser('list')
    parser_list.set_defaults(func=run_list)

    parser_errors = subparsers.add_parser('errors')
    parser_errors.add_argument("--within", type=int, default=None)
    parser_errors.add_argument("--full", action="store_true", default=False)
    parser_errors.add_argument("--out-base", default="pcawg_data")
    parser_errors.set_defaults(func=run_errors)

    parser_clean = subparsers.add_parser('clean')
    parser_clean.add_argument("--out-base", default="pcawg_data")
    parser_clean.set_defaults(func=run_clean)

    parser_set = subparsers.add_parser('set')
    parser_set.add_argument("state")
    parser_set.add_argument("ids", nargs="+")
    parser_set.set_defaults(func=run_set)

    parser_timing = subparsers.add_parser('timing')
    parser_timing.add_argument("--out-base", default="pcawg_data")
    parser_timing.set_defaults(func=run_timing)

    args = parser.parse_args()
    args.func(args)
