FROM    ubuntu

RUN     apt-get update
RUN     apt-get install -y python wget samtools tabix
RUN     mkdir /opt/bin
RUN     wget -O /opt/bin/MuSE http://odin.mdacc.tmc.edu/~wwang7/MuSE
RUN     chmod +x /opt/bin/MuSE
ADD     ./muse.py /opt/bin/

ENV     PATH /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/bin
