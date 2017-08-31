for i in `seq 300 332`; do
	aprun ./errorCalc fg /lustre/cray/ws8/ws/ipvober-awesome/test_31555_31777_512procs/test${i}/plot.dat /lustre/cray/ws8/ws/ipvmario-geneIV/distributedcombigrid/examples/gene_distributed/ref31888_6000/checkpoint error31555_31777_512procs_toRef_newWeibull_many.out ${i}

done

