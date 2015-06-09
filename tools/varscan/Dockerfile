FROM    java


RUN     apt-get update
RUN     apt-get install -y wget samtools

WORKDIR /opt
ENV     JAVA_JAR_PATH   /opt
RUN     wget http://downloads.sourceforge.net/project/varscan/VarScan.v2.3.8.jar
