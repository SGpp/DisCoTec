#!/usr/bin/env python3

import numpy as np
import functools
import math
import pandas as pd
from icecream import ic

import os
import re
# cf https://stackoverflow.com/questions/39293968/how-do-i-search-directories-and-find-files-that-match-regex
rootdir = "./"
regex = re.compile('(advect*.*out$)')


if __name__ == "__main__":
    df = pd.DataFrame(columns=['basis','dim', 'lmin', 'lmax', 't', 'L0error', 'L1error', 'L2error', 'dof'])
    for root, dirs, files in os.walk(rootdir):
        for file_found in files:
            if regex.match(file_found):
                ic(root, file_found)
                dof = 0
                with open(os.path.join(rootdir, root, file_found), "r") as f:
                    for line in f:
                        if line.startswith("["):
                            level =  re.split('[\[\W\]]+', line)
                            level = list(filter(None, level))
                            #ic(level)
                            level_sum = functools.reduce(lambda a,b : int(a)+int(b), level)
                            dof = dof + 2**level_sum
                        pass
                    last_line = line
                    #ic(last_line)
                    try:
                        numbers = [float(x) for x in last_line.split(', ')]
                        if numbers[0] != 1.0:
                            raise ValueError
                    except ValueError:
                        print("Not finished: ", root+"/"+file_found)
                        continue
                    #ic(numbers)

                parameters = re.split('[_\-D]', root)
                if (parameters[2] == "periodic"):
                    parameters[1] = parameters[1] + "_periodic"
                    parameters.pop(2)
                #ic(parameters)
                ls = [parameters[1], int(parameters[2]), int(parameters[4]), int(parameters[5])]
                ls = ls + numbers + [dof]
                row = pd.Series(ls, index=df.columns)
                #ic(row)
                df = df.append(row, ignore_index=True)

    df.sort_values(by=['basis', 'dim', 'dof'], ignore_index=True, inplace=True)
    df.to_csv("errors_advection.csv")

    df.sort_values(by=['dim', 'lmax'], ignore_index=True, inplace=True)
    for row in df.iterrows():
        #ic(row)

        pass
    df.to_csv("errors_advection.csv")
