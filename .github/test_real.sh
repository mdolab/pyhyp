#!/bin/bash
set -e
./reg_tests/ref/get-ref-files.sh
cd tests
testflo . -v -n 1

cd ../examples
ls -l */*.cgns
