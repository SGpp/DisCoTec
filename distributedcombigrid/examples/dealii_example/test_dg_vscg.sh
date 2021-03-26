#!/bin/bash

paramfile='ctparam'
errors='errors.txt'
log='log'
#create ctparam
echo "\n DG with combination now ncombi=10 und dt=0.1, e=9">> l2_error_cgdg.dat

for ((lmax=2;lmax<=4;lmax++)); do
    echo "lmax=$lmax"
    echo -e "[ct] \ndim=3 \nlmin= $lmax $lmax $lmax\nlmax= $lmax $lmax $lmax \nleval= $lmax $lmax $lmax
    p= 1 1 1\nncombi=100 \nFE=FE_DGQ \n DO_COMBINE=True \n[application]\ndt=0.01\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=true" > ctparam
    sh  ./run2.sh
    referencesolution="out/FE_DGQ/exact_solution_level$lmax,$lmax,$lmax.dat5"
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
        echo -e "[ct] \ndim=3 \nlmin= $lmin $lmin $lmin\nlmax=$lmax $lmax $lmax\nleval=$lmax $lmax $lmax
        p= 1 1 1\nncombi=10 \nFE=FE_DGQ \n DO_COMBINE=true \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=false" > ctparam
        sh  ./run2.sh
        fileOne="out/FE_DGQ/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$lmax,$lmax,${lmax}_combi_true.dat5"
        
        #now calculate the error
         
        ../../tools/errorCalc abs $fileOne $referencesolution l2_error_cgdg.dat ${lmin}_${lmax}

        #copy the output from the error file to errors.
    done
done

echo "\n CG with combination now ncombi=10 und dt=0.1, e=9">> l2_error_cgdg.dat

for ((lmax=2;lmax<=4;lmax++)); do
    echo "lmax=$lmax"
    echo -e "[ct] \ndim=3 \nlmin= $lmax $lmax $lmax\nlmax= $lmax $lmax $lmax \nleval= $lmax $lmax $lmax
    p= 1 1 1\nncombi=100 \nFE=FE_Q \n DO_COMBINE=True \n[application]\ndt=0.01\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=true" > ctparam
    sh  ./run2.sh
    referencesolution="out/FE_Q/exact_solution_level$lmax,$lmax,$lmax.dat0"
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
        echo -e "[ct] \ndim=3 \nlmin= $lmin $lmin $lmin\nlmax=$lmax $lmax $lmax\nleval=$lmax $lmax $lmax
        p= 1 1 1\nncombi=10 \nFE=FE_Q \n DO_COMBINE=true \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=false" > ctparam
        sh  ./run2.sh
        fileOne="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$lmax,$lmax,${lmax}_combi_true.dat0"
        
        #now calculate the error
         
        ../../tools/errorCalc abs $fileOne $referencesolution l2_error_cgdg.dat ${lmin}_${lmax}

        #copy the output from the error file to errors.
    done
done
