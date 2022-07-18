#!/bin/bash
#USAGE: ./4.fasttree.sh pathToFasttreemp fastaMultipleAlignment

#FASTTREEMP v.2.1.10 SSE3, OpenMP
FASTTREEMPPATH=$1
SAMPLE=$2

$FASTTREEMPPATH -nt -fastest $SAMPLE > "$SAMPLE.tre"
cat "$SAMPLE.tre" | perl -ne '$_=~s/:[0-9\.]+//g; $_=~s/\)[0-9\.]+/\)/g; print "$_"' > "$SAMPLE.tre-simple"
