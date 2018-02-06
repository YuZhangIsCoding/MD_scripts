#!/bin/bash
for i in `seq 1 12`; do
    if [ "$i" = "10" ]; then
        name="0.1""_"
        charge="0.101"
    elif [ "$i" -lt "10" ]; then
        name="0.0$i"
        charge="0.0$i""0$i"
    else
        name="0.$i"
        charge="0.$i$i"
    fi
    for j in TEA*; do
        if [[ $j == *"$name"* ]]; then
            cd $j
            python /Users/yuzhang/simulation/work/code/calc_potential.py -i Density.txt -s $charge
            cd ../
        fi
    done
done
