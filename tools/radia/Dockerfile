FROM java

RUN   apt-get update && apt-get install -y wget python python-pip samtools python-dev zlib1g-dev unzip
RUN		pip install pysam
RUN   pip install PyVCF

WORKDIR /opt

#Download and unpack radia, then remove reference data files
RUN   wget https://github.com/aradenbaugh/radia/archive/v1.1.0.tar.gz && tar xvzf v1.1.0.tar.gz && rm -rf radia-1.1.0/data/hg19
RUN		wget -O /usr/bin/blat http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat && chmod +x /usr/bin/blat

RUN   wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip && unzip snpEff_latest_core.zip
