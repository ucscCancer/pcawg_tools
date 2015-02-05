FROM ubuntu


RUN apt-get update && apt-get install -y gcc git wget make \
zlib1g-dev ncurses-dev cmake g++ unzip

WORKDIR /opt

RUN git clone --recursive https://github.com/genome/breakdancer.git && cd breakdancer && rm -rf .git
ENV PERL_MM_USE_DEFAULT 1
RUN apt-get install -y libstatistics-descriptive-perl libgd-dev libgd-graph-perl
RUN  perl -MCPAN -e 'install GD::Graph::histogram'

RUN cd breakdancer && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local && make install
