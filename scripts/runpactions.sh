for numModes in 3; do #{2..6}; do
    for nsamples in 1; do #{1..5}; do
        for seed in {4..10}; do #{1..3}; do
            prefix=${seed}'_'${nsamples}'_'${numModes};
            echo $prefix;
            PACTIONPREFIX='/scratch/data/nsdong2/projectPACTION/'

            TREES=$PACTIONPREFIX"simulations/multipaction/"$prefix"_tree_mode0.tsv"
            PROPS=$PACTIONPREFIX"simulations/multipaction/"$prefix"_props_mode0.out"
            for ((i=1; i<numModes; i++)); do
                TREES+=" None"
                PROPS+=" "$PACTIONPREFIX"simulations/multipaction/"$prefix"_props_mode"$i".out"
            done

            echo "python3 /scratch/data/nsdong2/projectPACTION/newpaction/src/paction.py -p "$PROPS" -t "$TREES" -m both -o /scratch/data/nsdong2/projectPACTION/results/multipaction/"${prefix}"_" >> runPactions.txt
        done
    done
done
