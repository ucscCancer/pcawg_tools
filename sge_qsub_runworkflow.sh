#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 32
#$ -l scratch=4000g


#Where Nebula is installed
NEBULA=`pwd`/nebula/

#The docker tag for the Galaxy docker image
GALAXY_TAG=galaxy

#input JSON block
INPUT=$1

#workflow to run
WORKFLOW=$2

SUDO="" # set to "--sudo" if docker must be called with sudo

hostname > ${INPUT}.host
echo "LOADING" > ${INPUT}.state
for a in images/*.tar; do
	./scripts/docker_check_load.sh $a
done

export PYTHONPATH=$NEBULA

echo "RUNNING" > ${INPUT}.state
$NEBULA/scripts/run_galaxy_workflow.py \
-g $GALAXY_TAG \
-w ${WORKFLOW} \
-d tool_data/ -t tools/ $SUDO \
$INPUT \
-td tool_data/ -d data/ --config run_config.yaml 2> ${INPUT}.err > ${INPUT}.out

echo "DONE" > ${INPUT}.state
