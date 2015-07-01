
FROM	ubuntu:12.04

RUN		apt-get update && apt-get install -y wget curl
RUN		cd /opt && wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.7/genetorrent-download_3.8.7-ubuntu2.207-12.04_amd64.deb
RUN		cd /opt && wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.7/genetorrent-common_3.8.7-ubuntu2.207-12.04_amd64.deb
RUN		cd /opt && wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.7/genetorrent-upload_3.8.7-ubuntu2.207-12.04_amd64.deb
RUN		apt-get install -y libcurl3 libxqilla6 python
RUN		apt-get install -y libboost-program-options1.48.0 libboost-system1.48.0  libboost-filesystem1.48.0 libboost-regex1.48.0

RUN		cd /opt && dpkg --install genetorrent-download_3.8.7-ubuntu2.207-12.04_amd64.deb genetorrent-common_3.8.7-ubuntu2.207-12.04_amd64.deb genetorrent-upload_3.8.7-ubuntu2.207-12.04_amd64.deb