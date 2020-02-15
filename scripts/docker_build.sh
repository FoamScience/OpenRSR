#!/bin/bash
# Build and Test the toolkit in a container based off of foamscience/openrsr  docker image

check_errs()
{
  # Function. Parameter 1 is the return code
  # Para. 2 is text to display on failure.
  if [ "${1}" -ne "0" ]; then
    echo "ERROR # ${1} : ${2}"
    # as a bonus, make our script exit with the right error code.
    exit ${1}
  fi
}

# Source FE4
source /opt/foam/foam-extend-4.0/etc/bashrc
pip install --user cpp-coveralls
set -ev

# Get to the directory, compile libraries (Opt mode)
cd /home/foam/OpenRSR; ./Allwmake

# A full test of all library classes
tests=`find /home/foam/OpenRSR -iname "test" -type d`
for t in $tests; do
    echo "Testing $t:"
    echo "------------------------------------------"
    pushd . > /dev/null
    cd $t
    wmake
    check_errs $? "Test didn't compile ..."
    ./*Test
    popd > /dev/null
done
coveralls --gcov-options '\-lp'
