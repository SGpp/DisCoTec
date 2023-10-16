#!/usr/bin/env python3

"""Script for evaluating errors and norms from hdf5 output (like generated from writeInterpolatedValues, writeInterpolationCoordinates) """

import argparse
import math
import h5py
from icecream import ic
import numpy as np


def paraboloid_function(coords):
    dim = len(coords)
    sign = 1. if dim % 2 == 1 else -1.
    result = sign
    for coord in coords:
        result *= coord * (coord - 1.)
    return result


def hyperplane_function(coords, simulation_time):
    result = 0.
    for d in range(len(coords)):
        result += coords[d] * 10**d
    # assume simulation_time combinations were performed with TaskCount
    result *= simulation_time
    return result


def gaussian_function(coords, time):
    exponent = 0.
    sigma = 1./3.
    sigma_squared = sigma * sigma
    sigma_squared_inv = 1. / sigma_squared
    for dim in range(len(coords)):
        # transform coordinates acccording to periodic boundary condition
        coords[dim] = math.fmod(1.0 + math.fmod(coords[dim] - time, 1.0), 1.0)
        assert(coords[dim] >= 0.)
        assert(coords[dim] <= 1.)
        exponent -= pow(coords[dim] - 0.5, 2)

    exponent *= sigma_squared_inv
    # leave out normalization, such that maximum is always 1
    return np.exp(exponent)  # / sqrt( pow(2*pi*sigma_squared, len(coords)))


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
    parser.add_argument(
        "--relative",
        default=False,
        action='store_true',
    )

    args = parser.parse_args()
    if 'relative' in vars(args) and 'solution' not in vars(args):
        raise ValueError("Relative error norms require a reference solution")

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
    num_samples = values_np.shape[0]

    # for this dataset, also get the time attribute
    assert("simulation_time" in values_dataset.attrs)
    simulation_time = values_dataset.attrs["simulation_time"]
    ic(simulation_time)

    if calculate_error_norms:
        # the interpolation coordinates are only required if we compute error norms
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
        number_dimensions = coordinates_dataset.shape[1]
        coordinates_np = np.array(coordinates_dataset)

        if reference_solution[0] == 'advection':
            reference_values_np = np.array(
                [gaussian_function(coords, simulation_time) for coords in coordinates_np])
        elif reference_solution[0] == 'paraboloid':
            reference_values_np = np.array(
                [paraboloid_function(coords) for coords in coordinates_np])
        elif reference_solution[0] == 'hyperplane':
            reference_values_np = np.array(
                [hyperplane_function(coords, simulation_time) for coords in coordinates_np])

        error_np = values_np - reference_values_np
        for p in [np.inf, 1, 2]:
            p_error_norm = np.linalg.norm(error_np, ord=p)
            if args.relative:
                exact_norm = np.linalg.norm(reference_values_np, ord=p)
                # ic(exact_norm, p_error_norm)
                p_error_norm = p_error_norm / exact_norm
            elif p in [1, 2]:
                p_error_norm = p_error_norm / num_samples
            print("Monte Carlo error integral of order " +
                  str(p) + " : " + str(p_error_norm))

    else:
        # otherwise, we just compute the Monte Carlo norms of the interpolated values
        # (assuming domain to be an arbitrary hypercube of volume 1)
        for p in [np.inf, 1, 2]:
            p_norm = np.linalg.norm(values_np, ord=p)
            if p in [1, 2]:
                p_norm = p_norm / num_samples
            print("Monte Carlo norm of order " + str(p) + " : " + str(p_norm))
