#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N pcawg
#$ -pe smp 32
#$ -l scratch=4000g


#Where Nebula is installed
NEBULA=`pwd`/nebula/

#The docker tag for the Galaxy docker image
GALAXY_TAG=bgruening/galaxy-stable:dev

SERVICE=$1
TASK=$2

hostname > ${TASK}.host
echo "LOADING" > ${TASK}.state
for a in images/*.tar; do
	./scripts/docker_check_load.sh $a
done

export PYTHONPATH=$NEBULA

echo "RUNNING" > ${TASK}.state

$NEBULA/bin/nebula run \
$SERVICE \
$TASK 2> ${TASK}.err > ${TASK}.out

rm -rf /scratch/nebula_work /scratch/nebula_cache /scratch/nebula_scratch
echo "DONE" > ${TASK}.state
