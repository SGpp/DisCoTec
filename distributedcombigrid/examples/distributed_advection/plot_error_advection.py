#!/usr/bin/env python3

from argparse import ArgumentError
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
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


if __name__ == "__main__":
    if (len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        raise ArgumentError
    if (len(sys.argv) > 2):
        t = sys.argv[2]
    else:
        # assume time 1.
        t = 1.

    fg_data = read_ndarray_from_file(filename)
    ic(np.max(fg_data))

    shape_ref = np.shape(fg_data)
    dim = len(shape_ref)
    assert(dim == 2)  # for plotting, assume 2D data
    extent = [0., 1., 0., 1.]

    sigma_squared_inv = 1./(1./3.*1./3.)

    def coordinate_shift_velocity1(coord):
        return math.fmod(1.+math.fmod(coord-t, 1.), 1.)-0.5

    def analytic_sln(coords):
        return math.exp(sigma_squared_inv*(
            -math.pow(coordinate_shift_velocity1(coords[0]), 2)
            - math.pow(coordinate_shift_velocity1(coords[1]), 2)))

    analytic_data = np.array([[analytic_sln([i/shape_ref[0], j/shape_ref[1]])
                               for j in range(shape_ref[1])] for i in range(shape_ref[0])])

    relative = False
    diff_a_fg = (analytic_data-fg_data)/(analytic_data if relative else 1.)
    ic(np.max(diff_a_fg))

    m = plt.imshow(analytic_data, extent=extent)
    plt.title("analytical")
    plt.colorbar(m, orientation='vertical')
    plt.show()
    m = plt.imshow(diff_a_fg, extent=extent)
    plt.title("error analytical-fg")
    plt.colorbar(m, orientation='vertical')
    plt.show()
