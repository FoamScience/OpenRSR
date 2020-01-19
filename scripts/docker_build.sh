#!/bin/bash
# Build and Test the toolkit in a container based off of foamscience/openrsr  docker image

set -ev

# System-wide access to FE4
source /opt/foam/foam-extend-4.0/etc/bashrc

# build libraries and solvers
cd $HOME/foam/OpenRSR; ./Allwmake

# For testing
cd solvers/pSwCoupledFoam; wmake
cd ../../tutorials; ls

# Run tests
