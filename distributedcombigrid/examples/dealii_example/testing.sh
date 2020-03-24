#!/bin/bash
leval=7
paramfile='ctparam'
errors='errors.txt'
log='log'
#create ctparam
echo -e "[ct] \ndim=3 \nlmin=$leval $leval $leval \nlmax=$leval $leval $leval \nleval=$leval $leval $leval 
p= 1 1 1\nncombi=10 \nFE=FE_Q \n \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1" > ctparam
#this generates the reference, fine discretized solution
sh . ./run2.sh

referencesolution="out/FE_Q/csmi_$leval,$leval,${leval}_ma_$leval,$leval,${leval}_ev_$leval,$leval,$leval.dat0"
echo $referencesolution

for ((lmax=2;lmax<=10;lmax++)); do
    echo "lmax=$lmax"
    for ((lmin=2;lmin<=lmax;lmin ++)); do
        if (($lmax==7))
        then 
        if (($lmin <6))
        then
        continue
        fi
        fi
        echo "lmin=$lmin" >> log
        #create ctparam
        echo -e "[ct] \ndim=3 \nlmin=$lmin $lmin $lmin \nlmax=$lmax $lmax $lmax \nleval=$leval $leval $leval 
        p= 1 1 1\nncombi=10 \nFE=FE_Q \n \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1" > ctparam
        sh  ./run2.sh
        fileOne="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$leval,$leval,$leval.dat0"
        #now calculate the error
        ../../tools/errorCalc abs   $fileOne $referencesolution l2_error.dat ${lmin}_${lmax}_$leval

        #copy the output from the error file to errors.
    done
done