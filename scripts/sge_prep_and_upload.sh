#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N pcawg_upload
#$ -pe smp 32
#$ -l scratch=4000g

./scripts/docker_check_load.sh images/pancancer_upload_download.tar

for a in $1/*; do
  if [ -d $a ]; then 
    ./scripts/run_vcfupload.sh ./scripts/prep_and_upload.sh $a
  fi
done
