#!/bin/bash
set -e
testflo . -v -n 1  --coverage --coverpkg pyhyp
cd examples
ls -l */*.cgns
