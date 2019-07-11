# Build on top of ubuntu:xenial, isn't this a little old?
FROM ubuntu:xenial

# Install dependencies for extend if building from sources
#RUN apt-get install git-core build-essential binutils-dev cmake flex \
#        zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
#        libxt-dev rpm mercurial graphviz python python-dev

# Else, install binaries and continue
RUN wget https://sourceforge.net/projects/foam-extend/files/foam-extend-4.0/Ubuntu_16.04/foam-extend-4.0_amd64_Ubuntu1604_f500917.deb/download -O extend.deb
RUN apt ./extend.deb

# Download & install openrsr sourses
RUN git clone -r https://github.com/FoamScience/OpenRSR.git
RUN cd OpenRSR && ./Allwmake

RUN cd OpenRSR.tutorials/BL-Case1A && pSwCoupledFoam
