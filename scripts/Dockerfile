# Build on top of ubuntu:xenial, isn't this a little old?
FROM ubuntu:xenial

# Install dependencies for extend if building from sources
RUN apt-get update && apt install -qq -y wget vim git-core build-essential binutils-dev cmake flex zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev libxt-dev rpm mercurial graphviz python3 python3-dev

# Else, install binaries and continue
RUN wget https://sourceforge.net/projects/foam-extend/files/foam-extend-4.0/Ubuntu_16.04/foam-extend-4.0_amd64_Ubuntu1604_f500917.deb/download -O extend.deb
RUN apt-get install -y ./extend.deb

# Create a user and act as him
RUN useradd -ms /bin/bash foam
USER foam
WORKDIR /home/foam

# Clone Repo
RUN git clone -b dev --single-branch https://github.com/FoamScience/OpenRSR
WORKDIR  OpenRSR

# Source FE4 installation and Build
RUN /bin/bash -c "source /opt/foam/foam-extend-4.0/etc/bashrc; ./Allwmake"
