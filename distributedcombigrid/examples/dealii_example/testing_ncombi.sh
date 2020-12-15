leval=5
paramfile='ctparam'
errors='errors.txt'
log='log'
#create ctparam
ncombi=100
    for ((lmax=4;lmax<=5;lmax++)); do
        for ((lmin=3;lmin<=lmax;lmin ++)); do
            if (($lmax==$lmin || $lmax==5))
            then
            continue
            
            fi
            #create ctparam
            dt=$((1 /ncombi))
            echo -e  "[ct] \ndim=3 \nlmin=$lmin $lmin $lmin \nlmax=$lmax $lmax $lmax \nleval=$leval $leval $leval 
            p= 1 1 1\nncombi=$ncombi \nFE=FE_Q \n DO_COMBINE=true \n[application]\ndt=0.01\n \n[manager]\nngroup = 1\nnprocs = 1" > ctparam
            sh  ./run2.sh
            fileOne="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$leval,$leval,${leval}_combi_true.dat0"
            echo -e "[ct] \ndim=3 \nlmin=$lmin $lmin $lmin \nlmax=$lmax $lmax $lmax \nleval=$leval $leval $leval 
            p= 1 1 1\nncombi=$ncombi \nFE=FE_Q \n DO_COMBINE=false \n[application]\ndt=0.01\n \n[manager]\nngroup = 1\nnprocs = 1" > ctparam
            sh  ./run2.sh
            fileTwo="out/FE_Q/csmi_$lmin,$lmin,${lmin}_ma_$lmax,$lmax,${lmax}_ev_$leval,$leval,${leval}_combi_false.dat0"
            #now calculate the error
            ../../tools/errorCalc abs   $fileOne  l2_error_combi.dat ${lmin}_${lmax}_$leval

            #copy the output from the error file to errors.
        done
    done
