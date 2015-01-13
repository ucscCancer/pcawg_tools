#!/bin/bash

#genomic data
wget -r ftp://ftp.sanger.ac.uk/pub/project/PanCancer/
mv ftp.sanger.ac.uk/pub/project/PanCancer/* ./

#data for ContEst
#wget http://www.broadinstitute.org/~kcibul/contest/hg19_population_stratified_af_hapmap_3.3.vcf.gz
#gunzip hg19_population_stratified_af_hapmap_3.3.vcf.gz
wget http://www.broadinstitute.org/~gsaksena/arrayfree_ContEst/arrayfree_ContEst/SNP6.hg19.interval_list
wget http://www.broadinstitute.org/~gsaksena/arrayfree_ContEst/arrayfree_ContEst/hg19_population_stratified_af_hapmap_3.3.fixed.vcf
wget http://www.broadinstitute.org/~gsaksena/arrayfree_ContEst/arrayfree_ContEst/gaf_20111020+broad_wex_1.1_hg19.bed

#data for MuTect
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/dbsnp_132_b37.leftAligned.vcf.gz
gunzip dbsnp_132_b37.leftAligned.vcf.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf.gz

gzcat Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz | perl -pe 's/^chr//' > Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf
gzcat 1000G_phase1.indels.hg19.sites.vcf.gz | perl -pe 's/^chr//' > 1000G_phase1.indels.hg19.sites.vcf.fixed.vcf

