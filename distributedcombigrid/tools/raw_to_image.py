#!/usr/bin/env python3

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import functools
from icecream import ic


def read_ndarray_from_file(filename):
    with open(filename+"_header", "r") as f:
        for line in f.readlines():
            if line.lower().startswith("dimensionality"):
                dimensionality = int(line.split()[-1])
            elif line.lower().startswith("extents"):
                extents = [int(x) for x in line.split()[1:]]
            elif line.lower().startswith("data type size"):
                data_type_size = int(line.split()[-1])
    ic(dimensionality, extents, data_type_size)
    assert(len(extents) == dimensionality)

    dtype = np.dtype(np.float64)
    assert(data_type_size == 8)
    filelength = functools.reduce(lambda a, b: a*b, extents)

    with open(filename, "rb") as f:
        numpy_data = np.fromfile(f, dtype)
        assert(len(numpy_data) == filelength)
    numpy_data = np.reshape(numpy_data, extents)
    return numpy_data

""" Pass the name of the .raw file that you want to visualize"""
if __name__ == "__main__":
    assert(len(sys.argv) == 2)
    path = sys.argv[1]
    ic(path)
    fg_data = read_ndarray_from_file(path)
    shape_ref = np.shape(fg_data)
    dim = len(shape_ref)
    assert(dim == 2)  # for plotting, assume 2D data (or adapt to plot slices of higher dimensional domain)
    extent = [0., 1., 0., 1.]

    fig = plt.figure()
    m = plt.imshow(fg_data, extent=extent)
    # plt.title("analytical")
    plt.colorbar(orientation='vertical', format='%.0e')
    # plt.show()
    plt.savefig(path[:-4] + ".pdf", format="pdf",bbox_inches='tight' , pad_inches=0.)

