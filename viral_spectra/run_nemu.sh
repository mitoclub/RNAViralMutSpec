#!/bin/bash

# NEMU=/home/kpotoh/mutational_spectra_of_SSU/scripts/nemu-core.nf
# CONFIG=/home/kpotoh/mutational_spectra_of_SSU/scripts/ssu.config

PATH_TO_INPUT=./nemu_input
PATH_TO_OUTPUT=./nemu_output

MAX_NJOBS=32
SLEEP_TIME=300  # secs
COUNTER=0

done_file=$(mktemp)
inp_lst_file=$(mktemp)

ls $PATH_TO_OUTPUT > $done_file
ls $PATH_TO_INPUT | cut -d '.' -f 1 > $inp_lst_file

total_samples=$(grep -vf $done_file $inp_lst_file)

for sample in $total_samples;
do 
    input_fasta_abs=$(realpath $PATH_TO_INPUT/$sample.fasta)
    workdir=$PATH_TO_OUTPUT/$sample

    if [ -e $workdir ]; then 
        echo "'$workdir' already exist: pass"
        continue
    fi
    let COUNTER++

    echo -e "$COUNTER) workdir: $workdir"
    mkdir -p $workdir
    cd $workdir
    cp $input_fasta_abs .

    # -resume
    # -q -bg
    nextflow -bg -q -c $CONFIG run $NEMU -with-trace --sequence $sample.fasta --outdir .
    cd - >/dev/null
    sleep 1

    if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
        echo sleep $SLEEP_TIME secs
        sleep $SLEEP_TIME
    fi

done

# rm -r nemu_output/*/work
