FROM java

RUN   apt-get update && apt-get install -y wget python python-pip samtools python-dev zlib1g-dev unzip
RUN   apt-get install -y vim
RUN   pip install pysam
RUN   pip install PyVCF

WORKDIR /opt

#Download and unpack radia
RUN   wget https://github.com/ucscCancer/radia/archive/v1.1.5.tar.gz && tar xvzf v1.1.5.tar.gz 

# blat
RUN   wget -O /usr/bin/blat http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat && chmod +x /usr/bin/blat

# Download and unpack snpEff
RUN   wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip && unzip snpEff_latest_core.zip

# cleanup
RUN   rm -rf v1.1.5.tar.gz snpEff_latest_core.zip snpEff/galaxy snpEff/examples 

