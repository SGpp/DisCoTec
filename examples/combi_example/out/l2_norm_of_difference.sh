#!/bin/bash

make -C ../../../tools
#TODO: replace these by different outputs of the same size
# e.g., to get a full grid solution, set all levels (lmin, lmax, leval) in ctparam to be the same
# and compare with a combi solution that is written at the same resolution (=same leval)
fileOne = solution_10.dat0
fileTwo = solution_10.dat0
../../../tools/errorCalc abs $fileOne $fileTwo l2_error.dat same_field_errors_should_be_zero
