#!/bin/bash
# Build and Test the toolkit in a container based off of foamscience/openrsr  docker image

set -ev

# Foam-Extend-4 is installed in a user (foam) home directory
. $HOME/foam/foam-4.0/etc/bashrc

# build libraries and solvers
cd $HOME/foam/OpenRSR; ./Allwmake

# For testing
cd solvers/pSwCoupledFoam; wmake
cd ../../tutorials; ls

# Run tests
