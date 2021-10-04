#!/bin/bash
set -e
cd tests
./ref/get-ref-files.sh
testflo . -v -n 1

cd ../examples
ls -l */*.cgns
