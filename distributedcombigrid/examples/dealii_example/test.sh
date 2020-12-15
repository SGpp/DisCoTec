echo -e "rm  solution/*
rm deal_config/*
mpirun -np 2 ./dealii_example" > run2.sh

#. ./test_cg_nocombi.sh
#. ./test_cg_combi.sh
#. ./test_dg_nocombi.sh
#. ./test_dg_combi.sh

. ./test_cg_parallel.sh
