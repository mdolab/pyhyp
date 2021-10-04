#!/bin/bash
set -e
cd tests
./ref/get-ref-files.sh
testflo . -v -n 1  --coverage --coverpkg pyhyp

cd ../examples
ls -l */*.cgns
