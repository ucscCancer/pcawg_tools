#!/bin/bash
export INDIVIDUAL_ID=$1
export BAM_TUMOR=$2
export BAM_TUMOR_BAI=$3
export BAM_NORMAL=$4
export BAM_NORMAL_BAI=$5
export REF_DIR=$6
export FINALRESULTSDIR=$7

ln -s $BAM_TUMOR tumor.bam
ln -s $BAM_TUMOR_BAI tumor.bam.bai
ln -s $BAM_NORMAL normal.bam
ln -s $BAM_NORMAL_BAI normal.bam.bai

ln -s $REF_DIR /cga/fh/pcawg_pipeline/refdata

PIPELINE=/cga/fh/pcawg_pipeline/pipelines/pcawg_pipeline_v5.py

export PIPETTE_SERVER_DIR=/cga/fh/pcawg_pipeline/utils/pipette_server

export COMMDIR=/cga/fh/pcawg_pipeline/jobResults_pipette/status
#OUTDIR contains the intermediate files
export OUTDIR=/cga/fh/pcawg_pipeline/jobResults_pipette/jobs/$INDIVIDUAL_ID
#FINALRESULTSDIR contains all the files that should be kept after the pipeline completes
#export FINALRESULTSDIR=/cga/fh/pcawg_pipeline/jobResults_pipette/results

rm -rf $COMMDIR
mkdir -p $COMMDIR

set -e

python3 $PIPETTE_SERVER_DIR/pipetteSynchronousRunner.py $COMMDIR $OUTDIR $PIPELINE $COMMDIR $OUTDIR $INDIVIDUAL_ID `pwd`/tumor.bam `pwd`/normal.bam # > run.out 2> run.err

find $OUTDIR -name pipette.module.usage.txt  | xargs  sh -c 'for f; do cat "$f" ; done' true |sort | uniq > $FINALRESULTSDIR/summary.usage.txt

mkdir $FINALRESULTSDIR/gnos_vcfs
cp $OUTDIR/links_for_gnos/*/*.vcf.gz $FINALRESULTSDIR/gnos_vcfs/

#collect the file outputs for return back to Broad
tar -cvhf $FINALRESULTSDIR/broad.tar.gz $OUTDIR/links_for_broad

#display any failing modules
cat summary.usage.txt
