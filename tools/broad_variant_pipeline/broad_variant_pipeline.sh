#!/bin/bash


export NAMESPACE=gsaksena
export PIPELINE=/cga/fh/pcawg_pipeline/pipelines/pcawg_pipeline_v4.py
export INDIVIDUAL_ID=$1
export BAM_TUMOR=$2
export BAM_NORMAL=$3

export PIPETTE_SERVER_DIR=/cga/fh/pcawg_pipeline/utils/pipette/pipette/server
export TIMESTAMP=`date +'%Y-%m-%d__%H-%M-%S'`
export COMMDIR=/cga/fh/pcawg_pipeline/tmp/pipette_commdirs/$NAMESPACE/$INDIVIDUAL_ID/run
export OUTDIR=/cga/fh/pcawg_pipeline/jobResults_pipette/$NAMESPACE/$INDIVIDUAL_ID/run

mkdir -p $COMMDIR

python3 $PIPETTE_SERVER_DIR/pipetteSynchronousRunner.py $COMMDIR $OUTDIR $PIPELINE $COMMDIR $OUTDIR $INDIVIDUAL_ID $BAM_TUMOR $BAM_NORMAL

find $OUTDIR -name pipette.module.usage.txt  | xargs  sh -c 'for f; do cat "$f" ; done' true |sort | uniq > total.summary.usage.txt

gunzip -c ${OUTDIR}/links_for_gnos/tabix_snowman_germline_indel/PCAWG-test-01.broad-snowman.DATECODE.germline.indel.vcf.gz > snowman.germline.indel.vcf
gunzip -c ${OUTDIR}/links_for_gnos/tabix_snowman_germline_sv/PCAWG-test-01.broad-snowman.DATECODE.germline.sv.vcf.gz > snowman.germline.sv.vcf
gunzip -c ${OUTDIR}/links_for_gnos/tabix_snowman_somatic_indel/PCAWG-test-01.broad-snowman.DATECODE.somatic.indel.vcf.gz > snowman.somatic.indel.vcf
gunzip -c ${OUTDIR}/links_for_gnos/tabix_snowman_somatic_sv/PCAWG-test-01.broad-snowman.DATECODE.somatic.sv.vcf.gz > snowman.somatic.sv.vcf
