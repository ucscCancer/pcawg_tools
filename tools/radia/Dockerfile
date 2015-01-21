FROM ubuntu


RUN     apt-get update
RUN     apt-get install -y wget python python-pip samtools python-dev zlib1g-dev
RUN		pip install pysam

WORKDIR /opt

RUN     wget https://github.com/aradenbaugh/radia/archive/v1.1.0.tar.gz
RUN		tar tvzf v1.1.0.tar.gz 
RUN		wget -O /usr/bin/blat http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat && chmod +x /usr/bin/blat 