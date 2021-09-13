#!/bin/bash
set -e

cd reg_tests
testflo . -v -n 1

cd ../examples
ls -l */*.cgns
