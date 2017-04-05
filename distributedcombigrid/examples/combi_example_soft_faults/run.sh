#!/bin/bash
read_properties ()
{
	file="$1"
	while IFS=" = " read -r key value; 
	do
		case "$key" in
			"ngroup") 
				ngroup=$value;;
			"nprocs") 
				nprocs=$value;;
		esac
	done < "$file"
}

ngroup=0
nprocs=0
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
read_properties "$SCRIPT_DIR/settings.ini"
mpiprocs=$((ngroup*nprocs+1))

mpirun.mpich -n "$mpiprocs" ./combi_example_soft_faults
