#!/bin/bash

paramfile='ctparam'
log='log'
#create ctparam

for ((lmax=2;lmax<=6;lmax++)); do
    echo -e "[ct] \ndim=3 \nlmin=$lmax $lmax $lmax \nlmax=$lmax $lmax $lmax \nleval=$lmax $lmax $lmax 
    p= 1 1 1\nncombi=1 \n DO_COMBINE=true \nFE=FE_DGQ \n \n[application]\ndt=1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=true" > ctparam
#this generates the reference, fine discretized solution
sh  ./run2.sh


    for ((lmin=1;lmin<=lmax;lmin ++)); do
        echo "lmin=$lmin" >> log
        #create ctparam
        echo -e "[ct] \ndim=3 \nlmin=$lmin $lmin $lmin \nlmax=$lmax $lmax $lmax \nleval=$lmax $lmax $lmax 
        p= 1 1 1\nncombi=10 \nDO_COMBINE=false \nFE=FE_DGQ \n \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact]\n makeExact=false" > ctparam
        sh  ./run2.sh
        for i in 0 1 2 3 4 5 6 7; do
            fileOne="out/FE_DGQ/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$lmax,$lmax,${lmax}_combi_false.dat$i"
            referencesolution="out/FE_DGQ/exact_solution_level$lmax,$lmax,$lmax.dat$i"
            #now calculate the error
            ../../tools/errorCalc abs  $referencesolution $fileOne l2_error_DG.dat ${lmin}_${lmax}_$i
        done
        #copy the output from the error file to errors.
    done
done