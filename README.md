pcawg_tools
===========

Docker builds, wrapper scripts and Galaxy tool config files related to the ICGC/TCGA PCAWG work ( http://pancancer.info/ )


Deployment
==========


Get Nebula for deployment:
```
git clone https://github.com/kellrott/nebula.git
```

Add to PYTHONPATH:
```
export PYTHONPATH=`pwd`/nebula
```

Obtain GATK then:
```
cp GenomeAnalysisTK.jar tools/gatk_bqsr/
```

Build Images:
```
python nebula/nebula/warpdrive.py build -o images/ tools/
```

Fetch Prebuilt Galaxy Image:
```
docker pull bgruening/galaxy-stable:dev
```


Cache Docker Galaxy image:
```
docker save bgruening/galaxy-stable:dev > images/galaxy.tar
```

Generate jobs:
```
./scripts/pcawg_wf_gen.py --ref-download /path/to/pcawg/docstore --create-service
```


Submit job:
```
qsub sge_qsub_runworkflow.sh pcawg.service pcawg.tasks/workflow_id
```


Debugging and Development
=========================
To manually instance galaxy (with PCAWG Tools loaded) for interactive analysis and testing
```
python ./nebula/nebula/warpdrive.py up -t tools -l data -a -p 8080 --sudo -c run_config.yaml -td tool_data -f -d

```


Tools
=====


pcap_tools
----------
A docker build with all of the tools needed to run the PCAP-Core scripts ( https://github.com/ICGC-TCGA-PanCancer/PCAP-core/ ).

contest
-------

genetorrent
-----------

sam_pileup
----------

delly
-----

muse
----

pindel
------

varscan
-------

mutect
------

synapse_interface
-----------------
