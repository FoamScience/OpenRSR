FROM ubuntu:xenial

RUN apk add git-core build-essential binutils-dev cmake flex \
        zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
        libxt-dev rpm mercurial graphviz python python-dev
