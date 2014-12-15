
FROM java

RUN apt-get update
RUN apt-get install -y zip wget

# We'll be working in /opt from now on
WORKDIR /opt

#
# Download and unpack ConTest
#

RUN wget http://www.broadinstitute.org/~gsaksena/arrayfree_ContEst/arrayfree_ContEst/Queue-1.4-437-g6b8a9e1-svn-35362.jar
RUN wget http://www.broadinstitute.org/~gsaksena/arrayfree_ContEst/arrayfree_ContEst/ContaminationPipeline.scala


