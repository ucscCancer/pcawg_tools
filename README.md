pcawg_tools
===========

Docker builds, wrapper scripts and Galaxy tool config files related to the ICGC/TCGA PCAWG work ( http://pancancer.info/ )


Main code Base
https://github.com/ucscCancer/pcawg_tools

Runner Project
https://github.com/kellrott/nebula

Workflow Engine Docker build
https://github.com/bgruening/docker-galaxy-stable

Reference Data
https://www.synapse.org/#!Synapse:syn3241088

WhiteList/Assignment Tracking
https://www.synapse.org/#!Synapse:syn4556289


All the tools are written as Galaxy tool wrappers (in the tool directory), with Dockerfile builds (except for the broad tool). The workflows are the .ga files in the workflow directory (in Galaxy format).

These can be used manually in a properly configured Galaxy environment. The trick is to automate deployment and instancing of the workflows. Nebula ( https://github.com/kellrott/nebula ) manages spinning up docker based galaxy instances, which in turn runs the workflow and manages the running of the tool containers. Once the work is done, the result data is copied out and the instance terminated.

There is also a workflow generator script included in each project that will check out tasks from the central whitelist (using the synapse table API), download reference files (from https://www.synapse.org/#!Synapse:syn3241088) and build the workflow request files.
I also use a 'doc store' to store large files and associated json meta-data documents, each with a uuid. Right now, its just a directory with predictable file names and associated json files.

The generator script will also scan the doc-store after the workflows are done, to build upload submissions (I haven't connected the data submission to the workflow, because I needed to review results).


Installation

-1 - System Dependencies (Ubuntu 14.04)

sudo apt-get install python-virtualenv python-dev

0 - Python dependencies
 If there are any questions about the python environment, create a 'virtualenv' one and install packages

virtualenv venv
. venv/bin/activate
pip install synapseclient
pip install pandas


1 - Download PCAWG tools
	git clone https://github.com/ucscCancer/pcawg_tools.git
This directory will be referred to at $PCAWG_DIR
You probably want to put $PCAWG_DIR in your .bashrc file
2- Checkout deployment version:
cd $PCAWG_DIR && git checkout 1.0.0
3 - Get Nebula for deployment (master branch should work for deployment)j:

git clone https://github.com/kellrott/nebula.git


This directory will be referred to as $NEBULA

If you are running any programs that use this code, you'll need to set PYTHONPATH:
export PYTHONPATH=$NEBULA

3.5 - Set Environmental Variables

Change for your system:

echo "export PCAWG_DIR=/home/ubuntu/gitroot/pcawg_tools">> ~/.bashrc
echo "export NEBULA=/home/ubuntu/gitroot/nebula">> ~/.bashrc
echo "export PYTHONPATH=\$NEBULA">> ~/.bashrc
source ~/.bashrc

4 - Build tool containers
  4.1 - Obtain GATK manually from https://www.broadinstitute.org/gatk/ then:
cp GenomeAnalysisTK.jar $PCAWG_DIR/tools/gatk_bqsr/


  4.2 - Run builder. This calls 'docker build' with the correct tag and dumps a copy of the image onto the file system.
In $PCAWG_DIR run the command
python $NEBULA/nebula/warpdrive.py build -o images/ tools/


4.3 - Obtain Broad Docker image
The Broad docker image needs to be built separately. Gordon will provide the download and decrypt instructions (its an encrypted file on the jamboree site).
Jamboree, the URL is  https://tcga-data-secure.nci.nih.gov/tcgafiles/pancan/forucsc/
And files to download:
/pancan/forucsc/refdata9
/pancan/forucsc/docker9
You will need to decrypt:
openssl aes-128-cbc -d -salt -in docker9 -out docker9.tar.gz
openssl aes-128-cbc -d -salt -in refdata9 -out refdata9.tar.gz
And Gordon will send you decryption passphrase.
You then build the docker image:
tar zxf refdata9.tar.gz; tar zxf docker9.tar.gz; cd docker && docker build -t broad_variant_pipeline .
Save the image to image directory
docker save broad_variant_pipeline > $PCAWG_DIR/images/docker_broad_variant_pipeline.tar
Need to check on the name and version used here...
4.4 Final list of containers needed
Broad
GATK co-cleaning
MuSE
???

5 - Build workflow engine
You can obtain a copy of the workflow engine from the Docker registry (docker pull bgruening/galaxy-stable), but its recommended you build it yourself, so that you can change the deployed user and group IDs to match your own account, which removes file access issues.
git clone --recursive https://github.com/bgruening/docker-galaxy-stable.git
Edit galaxy/Dockerfile to match your User and Group IDs (GALAXY_UID and GALAXY_GID on line 37 and 38). (You can find out your account number with 'id -u' and 'id -g')  
Then run
docker build -t bgruening/galaxy-stable galaxy
There’s also the DockerHub version of the above but it doesn’t allow you to change the running user (since it’s pre-built).  https://hub.docker.com/r/bgruening/galaxy-stable/

docker pull bgruening/galaxy-stable:dev

Cache Docker Galaxy image:
docker save bgruening/galaxy-stable > $PCAWG_DIR/images/galaxy.tar


6 - Set up configuration Files
GNOS keys and the BROAD reference data directory go under the $PCAWG_DIR/tool_data/files directory
They should be
$PCAWG_DIR/tool_data/files/cghub.key
$PCAWG_DIR/tool_data/files/icgc.key
$PCAWG_DIR/tool_data/files/refdata

7 - Generate Jobs
Register for jobs (note, you'll have to have your Synapse login with a ~/.synapseConfig file, see https://www.synapse.org/#!Synapse:syn2375225/wiki/ , you'll also have to request write access to the tracking table)
./scripts/pcawg_wf_gen.py register --count 15


You can confirm this worked by querying the Synapse table (https://www.synapse.org/#!Synapse:syn4556289 ). In this case I’m querying for my username “briandoconnor” and it finds 15 entries.

Notes about possible failure:
I tried the above on Ubuntu 14.04 and had to do the following to get it to work:
change the following in $PCAWG_DIR/scripts/synqueue.py

import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
to
import urllib3

sudo apt-get install python-dev
sudo pip install pandas

Now run the job creation script

Run the job creation script:
./scripts/pcawg_wf_gen.py gen --ref-download --create-service


Note: to use a separate directory for the workspace, add the parameter --work-dir </path/to/local/dir>, otherwise, the workflow will execute in a directory under /var/lib/docker, and may fill up the '/' partition
Submit job:
qsub sge_qsub_runworkflow.sh pcawg.service pcawg.tasks/workflow_id



Result files Uploads

Debugging and Development
To manually instance galaxy (with PCAWG Tools loaded) for interactive analysis and testing
In the $PCAWG_DIR run the command
python $NEBULA/nebula/warpdrive.py up -t tools -l data -a -p 8080 --sudo -c run_config.yaml -td tool_data -f -d



Runtime Notes
We're running on 32 core machines, scheduling a single donor per machine. For the PCAWG work, its been an average of 9 hours for GATK Indel realignment, 23 hours for BQSR, and then Broad pipeline averages 32 hours and the MuSE averages 11 hours (but they run in parallel). So a wall time average of 64 hours.
