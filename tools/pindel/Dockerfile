FROM ubuntu


RUN     apt-get update
RUN     apt-get install -y git wget make gcc zlib1g-dev ncurses-dev g++ python python-pip
RUN     pip install PyVCF

WORKDIR /opt

RUN		mkdir /opt/bin
RUN     wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
RUN     tar xvjf samtools-0.1.19.tar.bz2
RUN     cd /opt/samtools-0.1.19 && make
RUN     cp /opt/samtools-0.1.19/*.a /usr/local/lib/
RUN     mkdir /usr/local/include/bam
RUN     cp /opt/samtools-0.1.19/*.h /usr/local/include/bam/
RUN		cp /opt/samtools-0.1.19/samtools /usr/local/bin/

RUN     git clone https://github.com/genome/pindel.git
RUN     cd pindel && git checkout v0.2.5a7
RUN     cd pindel && ./INSTALL /opt/samtools-0.1.19
RUN     cp /opt/pindel/pindel* /usr/local/bin/