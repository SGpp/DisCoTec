#!/usr/bin/env python3

"""Script for visualizing the sparse grid min/max output (introduced in e9e287b9ea3)
"""

import sys
import os
import pandas as pd
from itertools import chain, combinations
from icecream import ic
# use pysgpp visualization, can clone SGPP (no installation required) and set
# export PYTHONPATH=$SGPP_PATH/combigrid/python:$PYTHONPATH
# I am currently using SGpp branch `extend_adaptive_combigrid_convenience`
from plotDeltas3d import plotDeltas3D

if __name__ == "__main__":
    """pass in the sg output file name to get the """
    if len(sys.argv) == 1:
        raise RuntimeError("no input file specified")
    if len(sys.argv) > 2:
        raise RuntimeError("too many command line arguments")
    filename = sys.argv[1]
    sgMinMax = pd.read_csv(
        filename, sep='[\s\[\]:,]+', engine='python', header=None, index_col=False)
    sgMinMax = sgMinMax.dropna(axis=1)
    numColumns = sgMinMax.shape[1]
    column_names = ["l_" + str(i) for i in range(numColumns - 2)]
    s = column_names.copy()
    column_names += ["min", "max"]
    sgMinMax.columns = column_names
    sgMinMax["min"] = sgMinMax["min"].abs()
    sgMinMax["max"] = sgMinMax["max"].abs()
    sgMinMax["delta"] = sgMinMax[["min", "max"]].max(axis=1)
    ic(sgMinMax)
    combinations = list(chain.from_iterable(combinations(s, r)
                        for r in range(len(s)+1)))
    combinations = [list(c) for c in list(combinations) if len(c) == 3]
    for c in combinations:
        ic(c)
        plotDeltas3D(sgMinMax, c, os.path.basename(filename)+"_"+str(c)+".pdf")
