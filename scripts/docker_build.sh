#!/bin/bash
# Build and Test the toolkit in a container based off of foamscience/openrsr  docker image

# Source FE4
source /opt/foam/foam-extend-4.0/etc/bashrc
set -ev

# Get to the directory, compile libraries (Opt mode)
cd /home/foam/OpenRSR; ./Allwmake

# A full test of all library classes
tests=`find /home/foam/OpenRSR -iname "test" -type d`
for t in $tests; do
    echo "Testing $t:"
    pushd . > /dev/null
    cd $t
    wmake && ./*Test
    popd > /dev/null
done
