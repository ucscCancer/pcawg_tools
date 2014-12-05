#!/bin/bash


#genomic data
wget -r ftp://ftp.sanger.ac.uk/pub/project/PanCancer/

#data for ContEst
wget wget http://www.broadinstitute.org/~kcibul/contest/hg19_population_stratified_af_hapmap_3.3.vcf.gz
gunzip hg19_population_stratified_af_hapmap_3.3.vcf.gz


wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/dbsnp_132_b37.leftAligned.vcf.gz
gunzip dbsnp_132_b37.leftAligned.vcf.gz