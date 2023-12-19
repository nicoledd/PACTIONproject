#!/bin/bash

for noise in 0;
    do
    for nsamples in 1;
        do
        for seed in 1;
            do
            prefix=${seed}'_'${nsamples}'_'${noise};
            echo $prefix;
            VAR='/scratch/data/nsdong2/projectPACTION/'
            VAR0=${VAR}'newpaction/src/test_suite/example_input/props_mode0.out'
            VAR1=${VAR}'newpaction/src/test_suite/example_input/props_mode1.out'
            VAR2=${VAR}'newpaction/src/test_suite/example_input/props_mode2.out'
            VAR3=${VAR}'newpaction/src/test_suite/example_input/tree_mode0.tsv'
            echo $VAR1

            python3 /scratch/data/nsdong2/projectPACTION/newpaction/src/paction.py \
                -p ${VAR0} ${VAR1} ${VAR2} \
                -t ${VAR3} None None \
                -m simultaneous \
                -o /scratch/data/nsdong2/projectPACTION/newpaction/src/test_suite/example_output/prog_
            done
        done
    done
