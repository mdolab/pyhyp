#!/bin/bash
set -e

cd tests
testflo . -v -n 1

cd ../examples
ls -l */*.cgns
