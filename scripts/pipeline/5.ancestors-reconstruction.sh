#!/bin/bash
#USAGE: ./6.ancestors-reconstruction.sh pathToPrank fastaMultipleAlignment binarySimplifiedTree
##USAGE: ./6.ancestors-reconstruction.sh pathToIqtree fastaMultipleAlignment binarySimplifiedTree

#PRANK v.170427
PRANKPATH=$1

#IQTREE v.1.6.11
#CAUTION! Long computation time.
#IQTREEPATH=$1

ALIGN=$2
BINSIMTREE=$3

$PRANKPATH -d=$ALIGN -t=$BINSIMTREE -o="$ALIGN.prank" -DNA -keep -showanc -uselogs -once -nomissing

#IQTREEPATH -s $ALIGN -te $BINSIMTREE -asr -m GTR+FO+R6+I -pre "$ALIGN.iqtree"
