#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import functools
from icecream import ic
from os.path import splitext


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
    assert (len(extents) == dimensionality)

    dtype = np.dtype(np.float64)
    assert (data_type_size == 8)
    filelength = functools.reduce(lambda a, b: a*b, extents)

    with open(filename, "rb") as f:
        numpy_data = np.fromfile(f, dtype)
        assert (len(numpy_data) == filelength)
    numpy_data = np.reshape(numpy_data, extents)
    return numpy_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file_name",
        type=str
    )
    args = parser.parse_args()
    filename = args.file_name
    ic(filename)

    nd_data = read_ndarray_from_file(filename)
    ic(np.min(nd_data), np.max(nd_data))

    nd_shape = np.shape(nd_data)
    dim_full = len(nd_shape)
    assert dim_full > 1

    longest_dim = np.argmax(nd_shape)
    nd_shape_shortened = list(nd_shape)
    nd_shape_shortened[longest_dim] = 0
    second_longest_dim = np.argmax(nd_shape_shortened)
    ic(longest_dim, second_longest_dim, nd_shape)

    # for the other dimensions, we assume the last index (for now)
    # TODO should be parameterizable
    nd_slice = [slice(-1, None)]*dim_full
    nd_slice[longest_dim] = slice(0, nd_shape[longest_dim])
    nd_slice[second_longest_dim] = slice(0, nd_shape[second_longest_dim])
    nd_slice = tuple(nd_slice)
    # ic(nd_slice)

    two_d_data = nd_data[nd_slice]
    ic(two_d_data.shape)
    dim = sum([1 if s > 1 else 0 for s in two_d_data.shape])
    assert (dim == 2)
    two_d_data = two_d_data.reshape(
        two_d_data.shape[longest_dim], two_d_data.shape[second_longest_dim])
    dim = len(two_d_data.shape)
    assert (dim == 2)  # for plotting, require 2D data
    extent = [0., 1., 0., 1.]

    fig = plt.figure()
    m = plt.imshow(two_d_data, extent=extent)
    plt.colorbar(orientation='vertical', format='%.2e')
    plt.savefig(splitext(filename)[0]+".pdf", format="pdf",
                bbox_inches='tight', pad_inches=0.)
