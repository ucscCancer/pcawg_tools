
FROM java

RUN apt-get update
RUN apt-get install -y samtools wget zip


# We'll be working in /opt from now on
WORKDIR /opt

# Link the picard tools to /opt/picard
RUN wget https://github.com/broadinstitute/picard/releases/download/1.122/picard-tools-1.122.zip
RUN unzip picard-tools-1.122.zip
RUN ln -s picard-tools-1.122 picard

ADD GenomeAnalysisTK.jar /opt/