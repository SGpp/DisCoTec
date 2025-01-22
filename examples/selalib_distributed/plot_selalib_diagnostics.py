#!/usr/bin/env python3

"""Script for plotting the selalib electrical energy and other QoIs / diagnostics
"""

import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess
import argparse

column_names = [
    "time",
    "mass $m$",
    "L2",
    "norm of el. charge density rho",
    "norm of el. potential $||\phi||$",
    "norm of el. field in x $||E_1||$",
    "norm of el. field in y $||E_2||$",
    "norm of el. field in z $||E_3||$",
    "norm of momentum in x $||v_1 \cdot f||$",
    "norm of momentum in y $||v_2 \cdot f||$",
    "norm of momentum in z $||v_3 \cdot f||$",
    "$||v_1^2||$",
    "$||v_2^2||$",
    "$||v_3^2||$",
]

# QoIs=["norm of el. potential $||\phi||$",  "mass $m$", "L2", "norm of el. charge density $|| \rho ||$", "norm(v1^2)" ]#, "norm(v2^2)", "norm(v3^2)"] #"norm(rho)" #"mass" #"norm(phi)"

QoI_indices = [4, 1, 2, 3, 11, 12, 13]
QoIs = [column_names[i] for i in QoI_indices]


def plot_selalib_qois(existing_dat_file_filenames, to_generate_ctparam_file):
    dataframes_existing = [
        pd.read_csv(name, sep=r"\s{1,}", names=column_names)
        for name in existing_dat_file_filenames
    ]
    if to_generate_ctparam_file is not None:
        # try to generate the postprocessing-combi-solution
        subprocess.run(
            ["./combine_selalib_diagnostics", to_generate_ctparam_file], check=True
        )
        # there are two lines of output for every time step in the .dat file;
        # sort them into different data frames
        dataframes_post = [
            pd.read_csv(
                "combined_vp_B2_3d3v.dat",
                sep=r"\s{1,}",
                names=column_names,
                skiprows=lambda x: x % 2 != m,
            )
            for m in range(3)
        ]
        m = 1
    else:
        dataframes_post = None

    for i in QoI_indices:
        QoI = column_names[i]
        if i == 4 or i == 5 or i == 6 or i == 7:
            for f, dataframe in enumerate(dataframes_existing):
                plt.semilogy(
                    dataframe["time"], dataframe[QoI], label="existing " + str(f)
                )
            if dataframes_post is not None:
                plt.semilogy(
                    dataframes_post[m]["time"],
                    dataframes_post[m][QoI],
                    label="combined" + str(m),
                )

        else:
            for f, dataframe in enumerate(dataframes_existing):
                plt.plot(dataframe["time"], dataframe[QoI], label="existing " + str(f))
            if dataframes_post is not None:
                plt.plot(
                    dataframes_post[m]["time"],
                    dataframes_post[m][QoI],
                    label="combined " + str(m),
                )

        plt.ylabel(QoI)
        plt.xlabel("simulation time $t$")

        # plt.xlim([0, 0.1])

        plt.legend()
        plt.savefig(
            "selalib_diagnostics_"
            + str(i)
            + "_"
            + os.path.split(os.getcwd())[-1]
            + ".svg",
            format="svg",
        )
        plt.show()

        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file_names",
        type=str,
        nargs="+",
        help="existing selalib .dat files",
        default=[],
        required=False,
    )
    parser.add_argument(
        "--ctparam_to_use",
        type=str,
        nargs=1,
        help="path to the ctparam file to generate .dat file for",
        default=None,
        required=False,
    )
    args = parser.parse_args()
    plot_selalib_qois(args.file_names, args.ctparam_to_use[0])
