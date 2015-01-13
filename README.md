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
docker export bgruening/galaxy-stable:dev > images/galaxy.tar
```

Manually instance galaxy (with PCAWG Tools loaded)
```
python nebula/nebula/warpdrive.py up -x tools/ -l data/ -a -f -w . -c
```

Generate jobs:
```
./scripts/build_jobs.py data/PCAWG_Data_Freeze_Train_2.0_Pilot_64.tsv jobs/
```


Submit job:
```
qsub sge_qsub_runworkflow.sh jobs/00ad0ffe-2105-4829-a495-1c2aceb5bb31.json workflows/Galaxy-Workflow-PCAWG_Pilot.ga
```

Tools
====


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
