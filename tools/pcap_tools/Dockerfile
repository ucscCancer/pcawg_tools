
FROM	ubuntu:12.04

RUN		apt-get update
RUN		apt-get install -y zlib1g-dev g++ dh-autoreconf libncurses-dev pkg-config git libgd2-noxpm-dev unzip python python-dateutil libipc-system-simple-perl python-lxml libjson-perl wget libxml-simple-perl curl python-yaml

#Get PCAP Tools
RUN		mkdir /opt/ICGC
RUN		cd /opt && git clone https://github.com/kellrott/PCAP-core.git
RUN		cd /opt/PCAP-core/ && ./setup.sh /opt/ICGC/

ADD		cgquery /opt/ICGC/bin/
ADD		pyscripts /opt/pyscripts

ENV		PERL5LIB /opt/ICGC/lib/perl5/
ENV		PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/ICGC/bin:/opt/pyscripts
