#!/bin/bash
set -e
./reg_tests/ref/get-ref-files.sh
testflo . -v -n 1
cd pyhyp_examples
ls -l */*.cgns
