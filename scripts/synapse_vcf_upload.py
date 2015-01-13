#!/usr/bin/env python

import sys
import json
import synapseclient
from synapseclient import File, Activity, Wiki
syn = synapseclient.login()

input_path = sys.argv[1]

with open(input_path + ".json") as handle:
    meta_data = json.loads(handle.read())

DST_FOLDER = 'syn3079564' #test upload folder

#Create Provenance log
provenance = Activity(name=meta_data['activity'],
                      desciption=meta_data['description'],
                      used = meta_data['used']
                      exectuted = meta_data['used']
                )
#prov = syn.store(prov)

name  = of.path.basename(input_path)
#Add metadata to files to be uploaded
f = File(input_path, name = name, parentId=DST_FOLDER)
f.dataType = meta_data['dataType']
f.fileType = meta_data['dataType']
f.variant_workflow = meta_data['workflow']
f.variant_workflow_version = meta_data['workflowVersion']
f.call_type = call_type
f.reference_build = meta_data['referenceBuild']
f.center_name = meta_data['center_name']
f.file_md5 = synapseclient.utils.md5_for_file(input_path)
f.study = 'PCAWG 2.0'
f.submitter_donor_id = meta_data['donor_id']
f.alignment_workflow_name='Workflow_Bundle_BWA (UCSC Implementation)'
f.alignment_workflow_source_url='https://github.com/kellrott/tcga_realign'
f.alignment_workflow_version='2.6.0'


#Store metadata and file to Synapse
#f = syn.store(f, activity = provenance)

#Add Description
#wiki = synapseclient.Wiki(TITLE, f, DESCRIPTION)
#wiki = syn.store(wiki)