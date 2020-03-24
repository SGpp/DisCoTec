#!/bin/bash
leval=7
paramfile='ctparam'
log='log'
#create ctparam
echo -e "[ct] \ndim=3 \nlmin=$leval $leval $leval \nlmax=$leval $leval $leval \nleval=$leval $leval $leval 
p= 1 1 1\nncombi=10 \nFE=FE_DGQ \n \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1" > ctparam
#this generates the reference, fine discretized solution
sh  ./run2.sh


for ((lmax=2;lmax<=leval;lmax++)); do
    echo "lmax=$lmax"
    for ((lmin=2;lmin<=lmax;lmin ++)); do
        echo "lmin=$lmin" >> log
        #create ctparam
        echo -e "[ct] \ndim=3 \nlmin=$lmin $lmin $lmin \nlmax=$lmax $lmax $lmax \nleval=$leval $leval $leval 
        p= 1 1 1\nncombi=10 \nFE=FE_DGQ \n \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1" > ctparam
        sh  ./run2.sh
        for i in 0 1 2 3 4 5 6 7; do
            fileOne="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$leval,$leval,$leval.dat$i"
            referencesolution="out/FE_DGQ/csmi_$leval,$leval,${leval}_ma_$leval,$leval,${leval}_ev_$leval,$leval,$leval.dat$i"
            #now calculate the error
            ../../tools/errorCalc abs  $referencesolution $fileOne l2_error.dat ${lmin}_${lmax}_$i
        done
        #copy the output from the error file to errors.
    done
done