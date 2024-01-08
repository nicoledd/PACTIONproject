for k in 10; do #{10 100}; do
    for n in 10; do #{10 50 100}; do
        for s in 1; do #{1..3}; do
            python3 ../src/paction.py -o ../results/ -g ../data/k${k}_n${n}_s${s}/T.dot;
        done
    done
done
