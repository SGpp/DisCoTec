#!/usr/bin/env python3

"""
Script to extract best, worst, mean, and std values from csv files
"""

import pandas as pd
import os
try:
    from icecream import ic
    ic("test ic")
except:
    def ic(*args):
        pass

import re

columns = ["t", "n", "nsamples", "best", "worst", "mean", "std"]
df_hawk = pd.DataFrame(columns=columns)
df_NG = pd.DataFrame(columns=columns)

for root, dirs, files in os.walk('.'):
    for filename in files:
        if (filename.find("broker_") == 0 and filename.find(".csv") > 17):
            ic(os.path.join(root, filename))
            temp = re.findall(r'\d+', filename)
            t, n = list(map(int, temp))
            small_series = pd.read_csv(filename, header=None)[0]
            small_series = small_series[small_series != 0]
            new_row = [t, n, small_series.size, small_series.min(), small_series.max(),
                       small_series.mean(), small_series.std()]
            ic(new_row)
            if (filename.find("_ng_") > 5):
                df_NG.loc[len(df_NG.index)] = new_row
            else:
                df_hawk.loc[len(df_hawk.index)] = new_row

df_NG = df_NG.astype({'t': 'int', 'n': 'int', 'nsamples': 'int'})
df_hawk = df_hawk.astype({'t': 'int', 'n': 'int', 'nsamples': 'int'})
df_NG.sort_values(by=["n","t"],inplace=True)
df_hawk.sort_values(by=["n","t"],inplace=True)
df_NG.to_csv("copy_uftp_NG_to_hawk.csv", index=False)
df_hawk.to_csv("copy_uftp_hawk_to_NG.csv", index=False)

