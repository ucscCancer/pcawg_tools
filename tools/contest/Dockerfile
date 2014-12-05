
FROM java

RUN apt-get update
RUN apt-get install -y zip wget

# We'll be working in /opt from now on
WORKDIR /opt

#
# Download and unpack ConTest
#

RUN wget http://www.broadinstitute.org/~kcibul/contest/contest-1.0.24530-bin.zip
RUN unzip contest-1.0.24530-bin.zip && mv contest-1.0.24530-bin/ContEst.jar /opt/

