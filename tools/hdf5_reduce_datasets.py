#!/usr/bin/env python3

"""Script for evaluating errors and norms from hdf5 output (like generated from writeInterpolatedValues, writeInterpolationCoordinates) """

import argparse
import math
import h5py
from icecream import ic
import numpy as np



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "values_files",
        nargs='*',
        type=str,
    )

    args = parser.parse_args()
    values_files = args.values_files
    ic(values_files)
    first_simulation_time = None
    first_dataset_shape = None
    values_np_accumulated = None
    file_name_accumulated = ""
    for values_file in values_files:
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
        if first_dataset_shape is None:
            first_dataset_shape = values_dataset.shape
        assert(values_dataset.shape == first_dataset_shape)
        values_np = np.array(values_dataset)
        num_samples = values_np.shape[0]

        # for this dataset, also get the time attribute
        assert("simulation_time" in values_dataset.attrs)
        simulation_time = values_dataset.attrs["simulation_time"]
        ic(simulation_time)
        if first_simulation_time is None:
            first_simulation_time = simulation_time
        assert(simulation_time == first_simulation_time)

        if values_np_accumulated is None:
            values_np_accumulated = values_np
        else:
            values_np_accumulated += values_np
        ic(values_np_accumulated[:5])

        if file_name_accumulated == "":
            file_name_accumulated = values_file.split("/")[-1].split(".")[0]
        else:
            file_name_accumulated += "_and_" + values_file.split("/")[-1].split(".")[0]

    # if everything is summed up, write out as h5 file
    # one hdf5 dataset in a group
    file_name_accumulated += ".h5"
    ic(file_name_accumulated)
    group_name = "summed_values"
    values_h5 = h5py.File(file_name_accumulated, 'w')
    values_h5.create_group(group_name)
    values_h5[group_name].create_dataset("values", data=values_np_accumulated)
    # set attribute simulation time
    values_h5[group_name]["values"].attrs["simulation_time"] = first_simulation_time
