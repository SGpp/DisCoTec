for i in `seq 301 301`; do
	aprun ./errorCalc fg /lustre/cray/ws8/ws/ipvober-awesome/test_31444_31888_512procs/test${i}/plot.dat /lustre/cray/ws8/ws/ipvmario-geneIV/distributedcombigrid/examples/gene_distributed/ref31888_6000/checkpoint error31444_31888_512procs_toRef_newWeibull_veryFew_new.out ${i}

done

