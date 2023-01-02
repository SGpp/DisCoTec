#!/usr/bin/env python3

"""Script for evaluating errors and norms from hdf5 output (like generated from writeInterpolatedValues, writeInterpolationCoordinates) """

import argparse
import h5py
from icecream import ic
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "values_file",
        nargs=1,
        type=str,
    )
    parser.add_argument(
        "--solution",
        choices=('advection', 'paraboloid', 'hyperplane', None),
        nargs=1,
        type=str,
        default=None,
    )
    parser.add_argument(
        "--coordinates",
        nargs=1,
        type=str,
        default=None,
    )

    args = parser.parse_args()
    reference_solution = args.solution
    values_file = args.values_file[0]
    calculate_error_norms = reference_solution is not None
    ic(values_file, reference_solution)

    values_h5 = h5py.File(values_file, 'r')

    # for now, assume a single group and a single solution for a single point in time
    assert (len(values_h5.keys()) == 1)
    group_name = list(values_h5.keys())[0]
    assert(isinstance(values_h5[group_name], h5py.Group))
    assert(len(values_h5[group_name].keys()) == 1)

    # get the first and only dataset in the group
    values_dataset_name = list(values_h5[group_name].keys())[0]
    values_dataset = values_h5[group_name][values_dataset_name]
    ic(values_dataset)
    assert (isinstance(values_dataset, h5py.Dataset))
    values_np = np.array(values_dataset)

    # for this dataset, also get the time attribute
    assert("simulation_time" in values_dataset.attrs)
    simulation_time = values_dataset.attrs["simulation_time"]
    ic(simulation_time)

    # the interpolation coordinates are only required if we compute error norms
    if calculate_error_norms:
        assert args.coordinates is not None
        coordinates_h5 = h5py.File(args.coordinates[0], 'r')
        assert (len(coordinates_h5.keys()) == 1)
        group_name = list(coordinates_h5.keys())[0]
        assert(isinstance(coordinates_h5[group_name], h5py.Group))
        assert(len(coordinates_h5[group_name].keys()) == 1)
        coordinates_dataset_name = list(coordinates_h5[group_name].keys())[0]
        coordinates_dataset = coordinates_h5[group_name][coordinates_dataset_name]
        ic(coordinates_dataset)
        assert (isinstance(coordinates_dataset, h5py.Dataset))
        assert values_dataset.shape[0] == coordinates_dataset.shape[0]
        coordinates_np = np.array(coordinates_dataset)

        raise NotImplementedError("Error norms not implemented yet")

    else:
        # otherwise, we just compute the Monte Carlo norms of the interpolated values
        # (assuming domain to be an arbitrary hypercube of volume 1)
        for p in [np.inf, 1, 2]:
            p_norm = np.linalg.norm(values_np, ord=p)
            print("Monte Carlo norm of order "+ str(p) +" : " + str(p_norm))