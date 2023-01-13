#!/bin/bash

TEMPLATE_DIR=./template
NNODESSYSTEM=128

paramfile="ctparam"
# allows to read the parameter file from the arguments.
# Useful for testing the third level combination on a single system
if [ $# -ge 1 ] ; then
   paramfile=$1
fi

runfile="submit.sh"

basisfunctions=(hat biorthogonal_periodic fullweighting_periodic)

for dim in {2..6}; do
   for test_maxlevel in {2..11}; do # assume that cfl condition will keep us from running any absurdly fine simulations
	for b in "${basisfunctions[@]}"; do
	FOLDER=advection_${b}_${dim}D_2-${test_maxlevel}
	# copy parameter file to new directory
	cp -rT $TEMPLATE_DIR $FOLDER
	# replace resolutions and parallelizations
	cd $FOLDER
	
	lmin=()
	lmax=()
	p=()
	for (( j=0; j<$dim; j++ )) do
		lmin+=(2)
		lmax+=($test_maxlevel)
		p+=(1)
	done
	
	p[0]=4
	p[1]=2
	p=${p[@]}
	lmin=${lmin[@]}
	leval=${lmax[@]}
	lmax=${lmax[@]}

        sed -i "s/dim.*/dim = $dim/g" $paramfile
	sed -i "s/lmin.*/lmin = $lmin/g" $paramfile
        sed -i "s/lmax.*/lmax = $lmax/g" $paramfile
        sed -i "s/leval.*/leval = $leval/g" $paramfile
	sed -i "s/basis.*/basis = $b/g" $paramfile
	sed -i "s/p =.*/p = $p/g" $paramfile
	sed -i "s/ngroup =.*/ngroup = 1/g" $paramfile
	#submit
	echo "submit $FOLDER"
	sbatch $runfile
	
	cd -
	done
	FOLDER=advection_fg_${dim}D_${test_maxlevel}-${test_maxlevel}
        # copy parameter file to new directory
        cp -rT $TEMPLATE_DIR $FOLDER
        # replace resolutions and parallelizations
        cd $FOLDER

        lmin=()
        lmax=()
	p=()
        for (( j=0; j<$dim; j++ )) do
                lmin+=($test_maxlevel)
                lmax+=($test_maxlevel)
		p+=(1)
        done

	p[0]=4
        p[1]=2
        p=${p[@]}
        lmin=${lmin[@]}
	leval=${lmax[@]}
        lmax=${lmax[@]}

        sed -i "s/dim.*/dim = $dim/g" $paramfile
        sed -i "s/lmin.*/lmin = $lmin/g" $paramfile
        sed -i "s/lmax.*/lmax = $lmax/g" $paramfile
        sed -i "s/leval.*/leval = $leval/g" $paramfile
        sed -i "s/basis.*/basis = $b/g" $paramfile
	sed -i "s/p =.*/p = $p/g" $paramfile
	sed -i "s/ngroup =.*/ngroup = 1/g" $paramfile


        #submit
        echo "submit $FOLDER"
        sbatch $runfile

        cd -
    done
done
