#!/bin/bash

paramfile='ctparam'
errors='errors.txt'
log='log'

dt=1

max_proc=16
p="1 1 1"
for comp_type_num in 0 1; do
    if (($comp_type_num == 0))
    then 
    comp_type="Q"
    elif (($comp_type_num == 1))
    then 
    comp_type="DGQ"
    fi

    for processes in 1 2 8 16; do
        if ((processes>1))
        then 
        if (($comp_type_num == 1))
        then
        continue
        fi
        fi
        if ((processes==1))
        then 
        p="1 1 1"
        elif ((processes==2))
        then 
        p="1 2 1"
        elif ((processes==8))
        then
        p="2 2 2"
        elif ((processes==16))
        then
        p="4 2 2"
        fi
        (( proc_1=$processes+1))
        echo -e "rm  solution/*
        rm deal_config/*
        mpirun -np $proc_1 ./dealii_example" > run2.sh
                
        for ncombi in 10 100; do
            if ((ncombi == 10)) 
            then
                dt=0.1
            else 
                dt=0.01
            fi
            for DO_COMBINE in "true" "false"; do
                echo "$comp_type DO_COMBINE=$DO_COMBINE now ncombi=$ncombi, p=$p">> L2_errors.dat
                for ((lmax=2;lmax<=9;lmax++)); do
                    echo "lmax=$lmax"
                    echo -e "[ct] \ndim=3 \nlmin= $lmax $lmax $lmax\nlmax= $lmax $lmax $lmax \nleval= $lmax $lmax $lmax
                    p= $p\nncombi=$ncombi \nFE=FE_$comp_type \n DO_COMBINE=$DO_COMBINE \n[application]\ndt=$dt\n \n[manager]\nngroup = 1\nnprocs = $processes \n [exact] \n makeExact=true" > ctparam
                    sh  ./run2.sh
                    referencesolution="out/FE_$comp_type/exact_solution_level$lmax,$lmax,$lmax.dat0"
                    for ((lmin=1;lmin<=lmax;lmin ++)); do
                        if ((lmin>3 ))
                        then
                        if ((lmin<lmax))
                        then
                        continue
                        fi
                        if ((lmax>7))
                        then 
                        continue
                        fi
                        fi
                        
                        echo "lmin=$lmin" >> log
                        #create ctparam
                        echo -e "[ct] \ndim=3 \nlmin= $lmin $lmin $lmin\nlmax=$lmax $lmax $lmax\nleval=$lmax $lmax $lmax
                        p= $p\nncombi=$ncombi \nFE=FE_$comp_type \n DO_COMBINE=$DO_COMBINE \n[application]\ndt=$dt\n \n[manager]\nngroup = 1\nnprocs = $processes \n [exact] \n makeExact=false" > ctparam
                        sh  ./run2.sh
                        fileOne="out/FE_$comp_type/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$lmax,$lmax,${lmax}_combi_$DO_COMBINE.dat0"
                        
                        ../../tools/errorCalc abs $fileOne $referencesolution L2_errors.dat ${comp_type}_${DO_COMBINE}_${ncombi}_${lmin}_${lmax}=
                    done
                done
            done
        done
    done
done