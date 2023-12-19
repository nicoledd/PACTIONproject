mkdir -p inputs
for seed in 1 2 3;
do
	for snv in 10 50 100;
	do
		for segments in 1 10 100;
		do
			echo $seed $snv $segments;
			mkdir -p ../data/k${segments}\_n${snv}\_s${seed};
			../clonesim/build/simulate -r -n ${snv} -k ${segments} -s ${seed} -S cnatrees.txt -dot ../data/k${segments}_n${snv}_s${seed}/T.dot > ../data/k${segments}_n${snv}_s${seed}/T.txt;
        done
	done
done