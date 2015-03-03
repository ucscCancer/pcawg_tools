#!/bin/bash

SAM_VER=0.1.16

wget http://sourceforge.net/projects/samtools/files/samtools/$SAM_VER/samtools-$SAM_VER.tar.bz2
tar xvjf samtools-$SAM_VER.tar.bz2
rm -rf samtools-$SAM_VER.tar.bz2
cd samtools-$SAM_VER && make
#cp /opt/samtools-$SAM_VER/*.a /usr/local/lib/
#mkdir /usr/local/include/bam
#cp /opt/samtools-$SAM_VER/*.h /usr/local/include/bam/
#cp /opt/samtools-$SAM_VER/samtools /usr/local/bin/
