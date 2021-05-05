#!/bin/bash

paramfile='ctparam'
errors='errors.txt'
log='log'
#create ctparam

echo "\n CG without combination now ncombi=100 und dt=0.01">> l2_error_combi.dat

for ((lmax=2;lmax<=7;lmax++)); do
    echo "lmax=$lmax"
    echo -e "[ct] \ndim=3 \nlmin= $lmax $lmax $lmax\nlmax= $lmax $lmax $lmax \nleval= $lmax $lmax $lmax
    p= 1 1 1\nncombi=100 \nFE=FE_Q \n DO_COMBINE=False \n[application]\ndt=0.01\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=true" > ctparam
    #sh  ./run2.sh
    referencesolution="out/FE_Q/exact_solution_level$lmax,$lmax,${lmax}_0.dat"
    for ((lmin=1;lmin<=lmax;lmin ++)); do
        if ((lmin>1 ))
        then
        if ((lmin<lmax))
        then
        continue
        fi
        fi
        echo "lmin=$lmin" >> log
        #create ctparam
        echo -e "[ct] \ndim=3 \nlmin= $lmin $lmin $lmin\nlmax=$lmax $lmax $lmax\nleval=$lmax $lmax $lmax
        p= 1 1 1\nncombi=10 \nFE=FE_Q \n DO_COMBINE=false \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=false" > ctparam
        sh  ./run2.sh
        fileOne="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$lmax,$lmax,${lmax}_combi_false10_0.dat"
        
        #now calculate the error
         
        ../../tools/errorCalc abs $fileOne $referencesolution L2_errors.dat Q_false_10_${lmin}_${lmax} CG

        #copy the output from the error file to errors.
    done
done
echo "\n CG without combination now ncombi=10 und dt=0.1">> l2_error_combi.dat

for ((lmax=2;lmax<=7;lmax++)); do
    echo "lmax=$lmax"
    echo -e "[ct] \ndim=3 \nlmin= $lmax $lmax $lmax\nlmax= $lmax $lmax $lmax \nleval= $lmax $lmax $lmax
    p= 1 1 1\nncombi=10 \nFE=FE_Q \n DO_COMBINE=True \n[application]\ndt=0.1\n \n[manager]\nngroup = 1\nnprocs = 1 \n [exact] \n makeExact=true" > ctparam
    #sh  ./run2.sh
    referencesolution="out/FE_Q/exact_solution_level$lmax,$lmax,${lmax}_0.dat"
    for ((lmin=1;lmin<=lmax;lmin ++)); do
        if ((lmin>1 ))
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
        fileOne="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$lmax,$lmax,${lmax}_combi_true10_0.dat"
        
        #now calculate the error
         
        ../../tools/errorCalc abs $fileOne $referencesolution L2_errors.dat Q_true_10_${lmin}_${lmax}= CG

        #copy the output from the error file to errors.
    done
done

