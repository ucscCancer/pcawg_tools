#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 32
#$ -l scratch=4000g


#Where Nebula is installed
NEBULA=`pwd`/nebula/

#The docker tag for the Galaxy docker image
GALAXY_TAG=bgruening/galaxy-stable:dev

SERVICE=$1
TASK=$2

SUDO="--sudo" # set to "--sudo" if docker must be called with sudo

hostname > ${TASK}.host
echo "LOADING" > ${TASK}.state
for a in images/*.tar; do
	./scripts/docker_check_load.sh $a
done

export PYTHONPATH=$NEBULA

echo "RUNNING" > ${TASK}.state

$NEBULA/bin/nebula run --hold \
$SERVICE \
$TASK 2> ${TASK}.err > ${TASK}.out

echo "DONE" > ${TASK}.state
