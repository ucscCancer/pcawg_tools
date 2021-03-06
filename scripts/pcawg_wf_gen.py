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
  "table_id" : "syn4556289",
  "primary_col" : "Submitter_donor_ID",
  "assignee_col" : "Assignee",
  "state_col" : "State"
}

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
            gnos_endpoint = urlparse(ent['meta']['Normal_WGS_alignment_GNOS_repos']).netloc
            task = GalaxyWorkflowTask("workflow_%s" % (ent['id']),
                workflow,
                inputs=dm,
                parameters={
                    'normal_bam_download' : {
                        "uuid" : ent['meta']['Normal_WGS_alignment_GNOS_analysis_ID'],
                        "gnos_endpoint" : gnos_endpoint,
                        "cred_file" : key_map[gnos_endpoint]
                    },
                    'tumor_bam_download' : {
                        "uuid" : ent['meta']['Tumour_WGS_alignment_GNOS_analysis_IDs'],
                        "gnos_endpoint" : gnos_endpoint,
                        "cred_file" : key_map[gnos_endpoint]
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
            galaxy="bgruening/galaxy-stable",
            sudo=args.sudo,
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
                donor_map[donor][ent['name']] = id

    for ent in synqueue.listAssignments(syn, list_all=True, **config):
        if ent['meta']['Submitter_donor_ID'] in donor_map:
            print ent['meta']['Submitter_donor_ID'], len(donor_map[ent['meta']['Submitter_donor_ID']])

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
    server_list = ["https://gtrepo-osdc-tcga.annailabs.com"]

    uuid_map = {}
    uuid_map['broad'] = dict( (a[1], a[0]) for a in synqueue.getValues(syn, "Broad_VCF_UUID",  **config).items() )
    uuid_map['muse']  = dict( (a[1], a[0]) for a in synqueue.getValues(syn, "Muse_VCF_UUID",  **config).items() )
    uuid_map['broad_tar'] = dict( (a[1], a[0]) for a in synqueue.getValues(syn, "Broad_TAR_UUID", **config).items() )

    analysis_re = re.compile(r'<analysis_id>(.*)</analysis_id>')

    found = {}
    for server in server_list:
        handle = urlopen(server + "/cghub/metadata/analysisId?study=tcga_pancancer_vcf")
        for line in handle:
            res = analysis_re.search(line)
            if res:
                gid = res.group(1)
                for p in uuid_map:
                    if gid in uuid_map[p]:
                        pid = uuid_map[p][gid]
                        if pid not in found:
                            found[pid] = {}
                        found[pid][p] = [server, gid]
    print "\n\n"
    for p in found:
        print p, "\t".join( "%s (%s : %s)" % (a[0], a[1][0], a[1][1]) for a in sorted(found[p].items(), key=lambda x:x[0]) )


def run_register(args):

    syn = synapseclient.Synapse()
    syn.login()
    synqueue.registerAssignments(syn, args.count, display=True, force=args.force, **config)

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
        wl_map[ent['id']] = ent['meta']

    uuid_map = {}
    uuid_map['broad'] = synqueue.getValues(syn, "Broad_VCF_UUID", orSet=lambda x: str(uuid.uuid4()), **config)
    uuid_map['muse']  = synqueue.getValues(syn, "Muse_VCF_UUID", orSet=lambda x: str(uuid.uuid4()), **config)
    uuid_map['broad_tar'] = synqueue.getValues(syn, "Broad_TAR_UUID", orSet=lambda x: str(uuid.uuid4()), **config)

    #scan through all of the docs
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
            #scan out the job metrics for this job
            if 'job' in entry and 'job_metrics' in entry['job']:
                job_id = entry['job']['id']
                tool_id = entry['job']['tool_id']
                job_info = { tool_id : {} }
                for met in entry['job']['job_metrics']:
                    job_info[tool_id][met['name']] = met['raw_value']
                job_map[donor][job_id] = job_info
            donor_tumor = wl_map[donor]['Tumour_WGS_alignment_GNOS_analysis_IDs']
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
                    dst_file = os.path.join(args.workdir, file_name)

                    shutil.copy(src_file, dst_file)
                    #if the files wasn't compressed already, go ahead and do that
                    if entry['extension'] == 'vcf':
                        subprocess.check_call( "bgzip -c %s > %s.gz" % (dst_file, dst_file), shell=True )
                        dst_file = dst_file + ".gz"

                    #add file to output map
                    if donor not in file_map[pipeline]:
                        file_map[pipeline][donor] = []
                    input_file = os.path.basename(dst_file)
                    file_map[pipeline][donor].append(input_file)
            else:
                if entry['name'] == "broad.tar.gz":
                    target = Target(uuid=entry['uuid'])
                    src_file = doc.get_filename(target)
                    file_map['broad_tar'][donor] = [ src_file ]


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
    for pipeline, donors in file_map.items():
        #for that pipeline go through every donor
        for donor, files in donors.items():
            #we're only outputing data for donors on the work list
            if donor in wl_map and result_counts[donor] == 3:
                #output the timing json
                timing_json = os.path.abspath(os.path.join(args.workdir, "%s.%s.timing.json" %(pipeline, donor)))
                with open( timing_json, "w" ) as handle:
                    handle.write(json.dumps( timing_map[donor] ) )

                #output the uploader script
                with open( os.path.join(args.workdir, "%s.%s.sh" %(pipeline, donor)), "w" ) as handle:
                    input_file = os.path.basename(dst_file)
                    urls = [
                        "%scghub/metadata/analysisFull/%s" % (wl_map[donor]['Normal_WGS_alignment_GNOS_repos'], wl_map[donor]['Normal_WGS_alignment_GNOS_analysis_ID']),
                        "%scghub/metadata/analysisFull/%s" % (wl_map[donor]['Tumour_WGS_alignment_GNOS_repos'], wl_map[donor]['Tumour_WGS_alignment_GNOS_analysis_IDs'])
                    ]
                    donor_tumor = wl_map[donor]['Tumour_WGS_alignment_GNOS_analysis_IDs']

                    if pipeline in ['broad', 'muse']:
                        prep_cmd_str = ""
                        for vcf in files:
                            prep_cmd_str += "tabix -p vcf %s\n" % (vcf)
                            prep_cmd_str += "mv %s.tbi %s.idx\n" % (vcf,vcf)
                            prep_cmd_str += "md5sum %s | awk '{print$1}' > %s.md5\n" % (vcf, vcf)
                            prep_cmd_str += "md5sum %s.idx | awk '{print$1}' > %s.idx.md5\n\n" % (vcf, vcf)

                        related_uuids = []
                        for p in uuid_map:
                            if p != pipeline:
                                related_uuids.append(uuid_map[p][donor])

                        submit_cmd_str = "perl -I /opt/gt-download-upload-wrapper/gt-download-upload-wrapper-2.0.11/lib"
                        submit_cmd_str += " /opt/vcf-uploader/vcf-uploader-2.0.5/gnos_upload_vcf.pl"
                        submit_cmd_str += " --metadata-urls %s" % (",".join(urls))
                        submit_cmd_str += " --vcfs %s " % (",".join(files))
                        submit_cmd_str += " --vcf-md5sum-files %s " % ((",".join( ("%s.md5" % i for i in files) )))
                        submit_cmd_str += " --vcf-idxs %s" % ((",".join( ("%s.idx" % i for i in files) )))
                        submit_cmd_str += " --vcf-idx-md5sum-files %s" % ((",".join( ("%s.idx.md5" % i for i in files) )))
                        submit_cmd_str += " --outdir %s.%s.dir" % (pipeline, donor_tumor)
                        submit_cmd_str += " --key %s " % (args.keyfile)
                        submit_cmd_str += " --upload-url %s" % (args.upload_url)
                        submit_cmd_str += " --study-refname-override %s" % (args.study)
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
                        #submit_cmd_str += " --skip-upload"

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

                        related_uuids = []
                        for p in uuid_map:
                            if p != pipeline:
                                related_uuids.append(uuid_map[p][donor])

                        submit_cmd_str = "perl -I /opt/gt-download-upload-wrapper/gt-download-upload-wrapper-2.0.11/lib"
                        submit_cmd_str += " /opt/vcf-uploader/vcf-uploader-2.0.5/gnos_upload_vcf.pl"
                        submit_cmd_str += " --metadata-urls %s" % (",".join(urls))
                        submit_cmd_str += " --tarballs %s " % (",".join(new_files))
                        submit_cmd_str += " --tarball-md5sum-files %s " % ((",".join( ("%s.md5" % i for i in new_files) )))
                        submit_cmd_str += " --outdir %s.%s.dir" % (pipeline, donor_tumor)
                        submit_cmd_str += " --key %s " % (args.keyfile)
                        submit_cmd_str += " --upload-url %s" % (args.upload_url)
                        submit_cmd_str += " --study-refname-override %s" % (args.study)
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
                        #submit_cmd_str += " --skip-upload"

                    handle.write(string.Template("""#!/bin/bash
set -ex
${PREP}
${SUBMIT}
echo $$? > $$0.submitted
#pushd ${SUBMIT_DIR}
#gtupload -v -c ${KEY} -u ./manifest.xml
#ECODE=$$?
#popd
#echo $$ECODE > $$0.uploaded
""").substitute(PREP=prep_cmd_str, SUBMIT=submit_cmd_str,
                            SUBMIT_DIR=os.path.join(os.path.abspath(args.workdir), "vcf", pipeline + "." + donor_tumor + ".dir", uuid_map[pipeline][donor] ),
                            KEY=args.keyfile
                    ) )


def run_list(args):
    syn = synapseclient.Synapse()
    syn.login()
    synqueue.listAssignments(syn, display=True, **config)

def run_set(args):
    syn = synapseclient.Synapse()
    syn.login()

    synqueue.setStates(syn, args.state, args.ids, **config)



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
    parser_gen.add_argument("--tool-data", default=os.path.abspath("tool_data"))
    parser_gen.add_argument("--tool-dir", default=os.path.abspath("tools"))

    parser_gen.set_defaults(func=run_gen)

    parser_submit = subparsers.add_parser('audit')
    parser_submit.add_argument("--out-base", default="pcawg_data")
    parser_submit.set_defaults(func=run_audit)

    parser_submit = subparsers.add_parser('gnos-audit')
    parser_submit.set_defaults(func=run_gnos_audit)

    parser_register = subparsers.add_parser('register',
                                       help='Returns set of new assignments')
    parser_register.add_argument("-c", "--count", help="Number of assignments to get", type=int, default=1)
    parser_register.add_argument("-f", "--force", help="Force Register specific ID", default=None)
    parser_register.set_defaults(func=run_register)

    parser_upload = subparsers.add_parser('upload-prep')
    parser_upload.add_argument("--workdir", default="work")
    parser_upload.add_argument("--keyfile", default="/keys/cghub.pem")
    parser_upload.add_argument("--upload_url", default="https://gtrepo-osdc-tcga.annailabs.com")
    parser_upload.add_argument("--out-base", default="pcawg_data")
    parser_upload.add_argument("--pipeline-src", default="https://github.com/ucscCancer/pcawg_tools")
    parser_upload.add_argument("--pipeline-version", default="1.0.0")
    parser_upload.add_argument("--pipeline-name", default="BROAD_MUSE_PIPELINE")
    parser_upload.add_argument("--study", default="tcga_pancancer_vcf_test")
    parser_upload.add_argument("--vm-instance-type", default="d2.8xlarge")
    parser_upload.add_argument("--vm-instance-cores", default="32")
    parser_upload.add_argument("--vm-instance-mem-gb", default="256")
    parser_upload.add_argument("--vm-location-code", default="27")
    parser_upload.set_defaults(func=run_uploadprep)

    parser_list = subparsers.add_parser('list')
    parser_list.set_defaults(func=run_list)

    parser_set = subparsers.add_parser('set')
    parser_set.add_argument("state")
    parser_set.add_argument("ids", nargs="+")
    parser_set.set_defaults(func=run_set)

    args = parser.parse_args()
    args.func(args)
