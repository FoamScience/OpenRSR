FROM alpine:latest

RUN apk add git build-base binutils-dev cmake flex \
        zlib-dev libncurses5-dev libiberty-dev \
        libxt-dev rpm mercurial graphviz python python-dev
