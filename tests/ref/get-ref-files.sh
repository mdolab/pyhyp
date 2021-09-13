#!/bin/bash
# This file will download the ref files for the pyHyp regression tests and extract them to the right place.

DIR=$(dirname $0)
TAR="ref_files.tar.gz"
wget -O $DIR/$TAR http://umich.edu/~mdolaboratory/repo_files/pyHyp/pyhyp_ref_files.tar.gz
tar -xzf $DIR/$TAR -C $DIR/../
rm $DIR/$TAR
