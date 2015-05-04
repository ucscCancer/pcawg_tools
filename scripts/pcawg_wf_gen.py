#!/usr/bin/env python

import os
import json
import argparse
import synqueue
import shutil
import synapseclient
import subprocess
from nebula.service import GalaxyService
from nebula.galaxy import GalaxyWorkflow
from nebula.docstore import from_url
from nebula.target import Target
from nebula.tasks import TaskGroup, GalaxyWorkflowTask
import tempfile

REFDATA_PROJECT="syn3241088"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("docstore", help="DocStore")
    parser.add_argument("--ref-download", action="store_true", default=False)
    parser.add_argument("--create-service", action="store_true", default=False)

    args = parser.parse_args()

    syn = synapseclient.Synapse()
    syn.login()

    docstore = from_url(args.docstore)

    if args.ref_download:
        #download reference files from Synapse and populate the document store
        for a in syn.chunkedQuery('select * from entity where parentId=="%s"' % (REFDATA_PROJECT)):
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

    data_mapping = {
        "reference_genome" : "genome.fa",
        "dbsnp" : "dbsnp_132_b37.leftAligned.vcf",
        "cosmic" : "b37_cosmic_v54_120711.vcf",
        "gold_indels" : "Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf",
        "phase_one_indels" : "1000G_phase1.indels.hg19.sites.fixed.vcf"
    }

    dm = {}
    for k,v in data_mapping.items():
        hit = None
        for a in docstore.filter(name=v):
            hit = a[0]
        if hit is None:
            raise Exception("%s not found" %s (v))
        dm[k] = { "uuid" : hit }

    workflow = GalaxyWorkflow(ga_file="workflows/Galaxy-Workflow-PCAWG_BROAD_MUSE.ga")

    config = synqueue.find_config()

    tasks = TaskGroup()

    for ent in synqueue.listAssignments(syn, **config):
        task = GalaxyWorkflowTask("workflow_%s" % (ent['id']),
            workflow,
            inputs=dm,
            parameters={
                'normal_bam_download' : {
                    "uuid" : ent['meta']['Normal_Analysis_ID'],
                    "gnos_endpoint" : "https://cghub.ucsc.edu",
                    "cred_file" : "/tool_data/files/cghub.key"
                },
                'tumor_bam_download' : {
                    "uuid" : ent['meta']['Tumour_Analysis_ID'],
                    "gnos_endpoint" : "https://cghub.ucsc.edu",
                    "cred_file" : "/tool_data/files/cghub.key"
                },
                'broad_variant_pipeline' : {
                    "broad_ref_dir" : "/tool_data/files/refdata",
                    "sample_id" : ent['meta']['Donor_ID']
                }
            },
            tags=[ "sample:%s" % (ent['meta']['Donor_ID']) ]
        )
        tasks.append(task)

    if not os.path.exists("pcawg.tasks"):
        os.mkdir("pcawg.tasks")
    for data in tasks:
        with open("pcawg.tasks/%s" % (data.task_id), "w") as handle:
            handle.write(json.dumps(data.to_dict()))

    if args.create_service:
        service = GalaxyService(
            docstore=docstore,
            galaxy="bgruening/galaxy-stable:dev",
            sudo=True,
            tool_data=os.path.abspath("tool_data"),
            tool_dir=os.path.abspath("tools"),
            smp=[
                ["MuSE", 8],
                ["pindel", 8],
                ["muTect", 8],
                ["delly", 4],
                ["gatk_bqsr", 12],
                ["gatk_indel", 12],
                ["bwa_mem", 12],
                ["broad_variant_pipline", 28]
            ]
        )
        with open("pcawg.service", "w") as handle:
            service.get_config().store(handle)
