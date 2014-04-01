#!/bin/bash

SRCDIR="`pwd`/src"
DATADIR="`pwd`/data"
RESDIR="`pwd`/results"

echo "SRCDIR = '$SRCDIR';" > src/load_paths.m
echo "DATADIR = '$DATADIR';" >> src/load_paths.m
echo "RESDIR = '$RESDIR';" >> src/load_paths.m

