FROM ubuntu


RUN apt-get update && apt-get install -y gcc git wget make \
zlib1g-dev ncurses-dev cmake g++


WORKDIR /opt

RUN git clone --recursive https://github.com/genome/tigra-sv.git && cd tigra-sv && rm -rf .git
ADD build_samtools.sh /opt/
RUN bash /opt/build_samtools.sh
ENV SAMTOOLS_ROOT /opt/samtools-0.1.16
# RUN cd /usr/src/gtest && cmake . && make && mv libg* /usr/lib/

RUN cd /opt/tigra-sv && sed -i s'/\/gsc\/pkg\/bio\/samtools\/samtools-0.1.6/\/opt\/samtools-0.1.16/' CMakeLists.txt
RUN cd /opt/tigra-sv && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install
