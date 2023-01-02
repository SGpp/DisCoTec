#!/usr/bin/env python3

"""Script for plotting the selalib electrical energy for Landau damping
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import subprocess

column_names = ["time", "mass $m$", "L2", "norm of el. charge density rho", "norm of el. potential $||\phi||$", "norm of el. field in x $||E_1||$", "norm of el. field in y $||E_2||$",
                "norm of el. field in z $||E_3||$", "norm of momentum in x $||v_1 \cdot f||$", "norm of momentum in y $||v_2 \cdot f||$", "norm of momentum in z $||v_3 \cdot f||$", "$||v_1^2||$", "$||v_2^2||$", "$||v_3^2||$"]

# QoIs=["norm of el. potential $||\phi||$",  "mass $m$", "L2", "norm of el. charge density $|| \rho ||$", "norm(v1^2)" ]#, "norm(v2^2)", "norm(v3^2)"] #"norm(rho)" #"mass" #"norm(phi)"

QoI_indices = [4, 1, 2, 3, 11, 12, 13]
QoIs = [column_names[i] for i in QoI_indices]

# dataframe_full = pd.read_csv("/import/epyc1.scratch/sse/pollinta/selalib_lower_bin/vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names)
dataframe_full = pd.read_csv(
    "/import/epyc1.scratch/sse/pollinta/selalib_bin/vp_B2_3d3v.dat_10", sep=r'\s{1,}', names=column_names)
dataframe_full_2 = pd.read_csv(
    "/import/epyc1.scratch/sse/pollinta/selalib_bin/vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names)
dataframe_full_2["time"] = dataframe_full_2["time"] + 10.
dataframe_full = dataframe_full.append(dataframe_full_2, ignore_index=True)

dataframe_combi = pd.read_csv(
    "./seladisco_leval_0/vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names)
    # "/import/epyc1.scratch/sse/pollinta/DisCoTec-selalib-example/seladisco_leval_0/vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names)
# try to generate the postprocessing-combi-solution
subprocess.run(["/home/pollinta/epyc/DisCoTec-selalib/distributedcombigrid/examples/selalib_distributed/combine_selalib_diagnostics", "ctparam"], check=True)
dataframes_post = [pd.read_csv("combined_vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names, skiprows=lambda x: x % 3 != m) for m in range(3)]
# dataframe_combi = pd.read_csv("./seladisco_leval_0/vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names)
# dataframe_combi_2 = pd.read_csv("/import/epyc1.scratch/sse/pollinta/DisCoTec-selalib-example-2/seladisco_leval_0/vp_B2_3d3v.dat", sep=r'\s{1,}', names=column_names)
# dataframe_combi_2["time"] = dataframe_combi_2["time"] + 10.
# dataframe_combi = dataframe_combi.append(dataframe_combi_2, ignore_index=True)

# print(dataframe_full)

for i in QoI_indices:
    QoI = column_names[i]
    if i == 4 or i == 5 or i == 6 or i == 7:
        plt.semilogy(dataframe_full["time"],
                     dataframe_full[QoI], label="full grid")
        plt.semilogy(dataframe_combi["time"],
                     dataframe_combi[QoI], label="combi")
        m = 2
        plt.semilogy(dataframes_post[m]["time"],
                     dataframes_post[m][QoI], label="post_" + str(m))
        m = 0
        #for dataframe_post in dataframes_post:
        #    plt.semilogy(dataframe_post["time"],
        #             dataframe_post[QoI], label="post_" + str(m))
        #    m += 1

    else:
        plt.plot(dataframe_full["time"],
                 dataframe_full[QoI], label="full grid")
        plt.plot(dataframe_combi["time"], dataframe_combi[QoI], label="combi")
        m = 2
        plt.plot(dataframes_post[m]["time"],
                     dataframes_post[m][QoI], label="post_" + str(m))
        #m = 0
        #for dataframe_post in dataframes_post:
        #    plt.plot(dataframe_post["time"], dataframe_post[QoI], label="post_" + str(m))
        #    m += 1
        # plt.plot(dataframe_combi["time"][:10], dataframe_combi[QoI][:10], label="combi")

    plt.ylabel(QoI)
    plt.xlabel("simulation time $t$")

    # plt.xlim([0, 0.1])

    plt.legend()
    plt.savefig("landau_"+str(i)+"_"+os.path.split(os.getcwd())[-1]+".svg", format="svg")
    plt.show()

    plt.clf()
