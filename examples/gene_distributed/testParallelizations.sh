#!/bin/bash

#declare -a p=("p = 1 1 1 1 2 1" "p = 1 1 1 2 1 1" "p = 1 1 2 1 1 1" "p = 1 2 1 1 1 1" "p = 2 1 1 1 1 1")
declare -a p=("p = 1 1 1 1 2 1" "p = 1 1 1 2 1 1" "p = 1 1 2 1 1 1")
n=${#p[@]} 
mkdir $1
cp manager $1/
cp -r template $1/
cp ctparam $1/
cp run.sh $1/
cp preproc.py $1/
cp errorCalc $1/
cp -r lib $1
cd $1/
for (( i=0; i<$n; i++ ))
do
    mkdir $i
    cp manager $i/
    cp -r template $i/
    cp -r lib $i/
    cp ctparam $i/
    cp run.sh $i/
    cp preproc.py $i/
    cd $i
    sed -i.bak "0,/p =.*/s/p =.*/${p[$i]}/1" ctparam
    ./run.sh
    cd ..
done

for (( i=0; i<$((n - 1)); i++ ))
do
    for (( j=$((i+1)); j<$n; j++ ))
    do
        ./errorCalc ff $i/plot.dat0 $j/plot.dat0 error.out "$i , $j"
    done
done
