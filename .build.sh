#!/bin/bash
set -ev

source /home/foam/foam/foam-4.0/etc/bashrc
cd /home/foam/OpenRSR; ./Allwmake
cd solvers/pSwCoupledFoam; wmake
cd ../../tutorials; ls
