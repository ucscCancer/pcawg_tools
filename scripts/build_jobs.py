#!/usr/bin/env python

import os
import sys
import csv
import json

ds_map = {
    "reference_genome" :
        {"uuid" : "6a9cf2d3-03e0-48ea-881f-e98ce685ddca"},
    "dbsnp" :
        {"uuid" : "c1abf048-7116-11e4-a62f-3417ebb385f7"},
    "cosmic" :
        {"uuid" : "5194d0a8-7117-11e4-871f-3417ebb385f7"},
    "pop_hap_vcf" :
        {"uuid" : "9e41d576-066f-43c8-9f2e-5cafa5bc30fd"},
    "platform_intervals" :
        {"uuid" : "48afa8c6-d74a-4c60-9002-fb82dfe0e449"},
    "gene_bed" :
        {"uuid" : "a2b3abb6-97a7-40fc-846c-9d7e8148ac46"},
    "gold_indels" :
        {"uuid" : "95a54557-1893-4410-879c-41c10bea8991"},
    "phase_one_indels" :
        {"uuid" : "039f2e11-e2a6-4270-8947-2b6c0f0752c4"}
}

with open(sys.argv[1]) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row['Normal GNOS endpoint'] == "https://cghub.ucsc.edu":
            pipeline = "cghub"
            out = {
                "parameters" : {
                    'normal_bam_download' : {
                        "uuid" : row['Normal Analysis ID'],
                        "gnos_endpoint" : "https://cghub.ucsc.edu",
                        "cred_file" : "/tool_data/files/cghub.key"
                    },
                    'tumor_bam_download' : {
                        "uuid" : row['Tumour Analysis ID'],
                        "gnos_endpoint" : "https://cghub.ucsc.edu",
                        "cred_file" : "/tool_data/files/cghub.key"
                    }
                },
                "ds_map" : ds_map,
                "tags" : [
                    "sample:%s" % (row['Donor ID'])
                ]
            }
        else:
            pipeline = "irods"
            out = {
                "parameters": {
                    "tumor_bam_download": {
                        "cred_file": "/tool_data/files/bazaar.yaml", 
                        "query|method" : "filter",
                        "query|filters_0|method|method_type" : "meta",
                        "query|filters_0|method|name" : "ANALYSIS_ID",
                        "query|filters_0|method|value" : row['Tumour Analysis ID'],
                        "query|filters_1|method|method_type" : "meta",
                        "query|filters_1|method|name" : "FILE_TYPE",
                        "query|filters_1|method|value" : "bam"
                    },
                    "normal_bam_download": {
                        "cred_file": "/tool_data/files/bazaar.yaml", 
                        "cred_file": "/tool_data/files/bazaar.yaml", 
                        "query|method" : "filter",
                        "query|filters_0|method|method_type" : "meta",
                        "query|filters_0|method|name" : "ANALYSIS_ID",
                        "query|filters_0|method|value" : row['Normal Analysis ID'],
                        "query|filters_1|method|method_type" : "meta",
                        "query|filters_1|method|name" : "FILE_TYPE",
                        "query|filters_1|method|value" : "bam"
                    }
                }, 
                "ds_map" : ds_map,
                "tags": [
                    "sample:%s" % (row['Donor ID'])
                ]
            }
            
        if out is not None:
            with open(os.path.join(sys.argv[2], pipeline + "." + row['Donor ID'] + ".json" ), "w") as ohandle:
                ohandle.write(json.dumps(out, indent=4))
