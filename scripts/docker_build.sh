#!/bin/bash
# Build and Test the toolkit in a container based off of foamscience/openrsr  docker image

check_errs()
{
  # Function. Parameter 1 is the return code
  # Para. 2 is text to display on failure.
  if [ "${1}" -ne "0" ]; then
    echo "ERROR # ${1} : ${2}"
    # as a bonus, make our script exit with the right error code.
    exit "${1}"
  fi
}

# Source FE4 in debug mode
source /opt/foam/foam-extend-4.0/etc/bashrc
set -ev

# Install Lcov for coverage reports
apt install -qq -y lcov curl jq

# Get to the directory, compile libraries (Opt mode)
cd /home/foam/OpenRSR; ./Allwmake

# A full test of all library classes
# While generating partial coverage reports
echo "####################################################################"
ls /home/foam/OpenRSR/libs/wellModels/lnInclude
cat /home/foam/OpenRSR/libs/wellModels/Make/files
echo "####################################################################"
tests=$(find /home/foam/OpenRSR -iname "test" -type d)
for t in $tests; do
    echo "Testing $t:"
    echo "------------------------------------------"
    pushd . > /dev/null
    cd "$t"
    wmake
    check_errs $? "Test didn't compile ..."
    ./*Test
    lcov -c --directory Make/linux64GccDPOpt --output-file coverage-test.info
    popd > /dev/null
done

# Combine all generated reports
find . -iname coverage-test.info -exec cat {} >> main-coverage.info +

# Filter coverage results
# We don't care for testing Foam-Extend and Catch code
lcov --extract  main-coverage.info "/home/foam/OpenRSR/*" -o lcov.info
lcov --remove  lcov.info "*catch2*" -o lcov.info
#genhtml lcov.info --output-directory coverage-output

# Send final report to Codacy
