for i in `seq 75 79`; do
	aprun ./errorCalc ff /lustre/cray/ws8/ws/ipvober-awesome/test_31666_31888_512procs/test${i}/plot.dat /lustre/cray/ws8/ws/ipvober-awesome/test_31444_31888_512procs_noFaults/test50/plot.dat error31666_31888_512procs.out ${i}

done
