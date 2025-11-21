#!/bin/bash
# USAGE: ./tree_simplifier.sh pathToTree
# remove distances from tree

SAMPLE=$1

cat $SAMPLE | perl -ne '$_=~s/:[0-9\.]+//g; $_=~s/\)[0-9\.]+/\)/g; print "$_"' > "$SAMPLE-simple"
