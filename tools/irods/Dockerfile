FROM ubuntu

RUN    apt-get update
RUN    apt-get install -y git libbz2-dev python-pip wget
RUN    pip install git+git://github.com/iPlantCollaborativeOpenSource/pycommands.git
RUN    pip install pyyaml
RUN    wget -P /opt/ ftp://ftp.renci.org/pub/irods/releases/4.0.3/irods-icommands-4.0.3-64bit.deb

RUN    apt-get install -y `dpkg -I /opt/irods-icommands-4.0.3-64bit.deb | sed -n 's/^ Depends: //p' | sed 's/,//g'`
RUN    dpkg -i /opt/irods-icommands-4.0.3-64bit.deb 
