#!/bin/bash

paramfile='ctparam'
errors='errors.txt'
log='log'
#create ctparam
echo "\n DG with combination now ncombi=10 und dt=0.1 and 2d">> l2_error_combi.dat


for ((lmax=2;lmax<=7;lmax++)); do
    echo "lmax=$lmax"
    echo -e "[ct] \ndim=2 \nlmin= $lmax $lmax \nlmax= $lmax $lmax  \nleval= $lmax $lmax 
    p= 1 1\nncombi=10 \nFE=FE_DGQ \n DO_COMBINE=False \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=true" > ctparam
    sh  ./run2.sh
    referencesolution="out/FE_DGQ/exact_solution_level$lmax,$lmax.dat0"
    for ((lmin=1;lmin<=lmax;lmin ++)); do
        if ((lmin>3 ))
        then
        if ((lmin<lmax))
        then
        continue
        fi
        fi
        echo "lmin=$lmin" >> log
        #create ctparam
        echo -e "[ct] \ndim=2 \nlmin= $lmin $lmin \nlmax=$lmax $lmax \nleval=$lmax $lmax 
        p= 1 1\nncombi=10 \nFE=FE_DGQ \n DO_COMBINE=true \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=false" > ctparam
        sh  ./run2.sh
        fileOne="out/FE_DGQ/csmi_$lmin,${lmin}_ma_$lmax,${lmax}_ev_$lmax,${lmax}_combi_true.dat0"
        
        #now calculate the error
         
        ../../tools/errorCalc abs $fileOne $referencesolution l2_error_combi.dat ${lmin}_${lmax}

        #copy the output from the error file to errors.
    done
done

