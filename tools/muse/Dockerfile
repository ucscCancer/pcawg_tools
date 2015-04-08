FROM    ubuntu

RUN     apt-get update
RUN     apt-get install -y python python-pip wget samtools tabix
RUN     pip install PyVCF
RUN     mkdir /opt/bin
RUN     wget -O /opt/bin/MuSEv0.9.9.5 http://bioinformatics.mdanderson.org/Software/MuSE/MuSEv0.9.9.5

RUN     chmod +x /opt/bin/MuSE*
ADD     ./muse.py /opt/bin/

ENV     PATH /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/bin
