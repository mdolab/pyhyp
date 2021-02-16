#!/bin/bash
set -e
testflo . -v -n 1
cd examples
ls -l */*.cgns
