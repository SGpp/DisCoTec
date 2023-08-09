#!/usr/bin/python3

import argparse
import json
import numpy as np
import pandas as pd
from tqdm import tqdm

io_only_event = "write SG"

# parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("--input_files", nargs='+', help="the input file", required=True)
parser.add_argument("--no_compute_per_rank_statistics", help="compute the per rank statistics", action="store_true", default=False)
parser.add_argument("--no_compute_per_step_statistics", help="compute the per step statistics", action="store_true", default=False)
parser.add_argument("--no_compute_total_statistics", help="compute the total statistics", action="store_true", default=False)
parser.add_argument("--io_only", help="only output IO ranks; a rank is an IO rank if it has the event \"{}\"".format(io_only_event), action="store_true", default=False)
parser.add_argument("--no_tqdm",help="display no progress bars", action="store_true", default=False)
#parser.add_argument("--groups_from_path",help="infer groups from json file name instead of file content", action="store_true", default=False)

args = parser.parse_args()

if args.no_compute_per_rank_statistics and args.no_compute_per_step_statistics and args.no_compute_total_statistics:
    raise ValueError("Not all of \"--no_compute_per_rank_statistics\", \"--no_compute_per_step_statistics\", and \"--no_compute_total_statistics\" must be set!")

data_per_rank = {}

for file in args.input_files:
    print("Reading data from file: \"{}\"".format(file))
    input_data = json.load(open(file))

    # iterate JSON object
    for rank in range(len(input_data)):
        input_rank_name = "rank" + str(rank)
        group = int(input_data[input_rank_name]["attributes"]["group"])
        output_rank_name = "rank" + str(len(input_data) * group + rank)

        # collect statistics per rank
        if output_rank_name not in data_per_rank:
            data_per_rank[output_rank_name] = {}

        for event in input_data[input_rank_name]["events"]:
            event_data = np.asarray(input_data[input_rank_name]["events"][event])
            event_data = np.apply_along_axis(lambda x: (x[1] - x[0]) * 1e-6, 1, event_data)

            # collect data for statistic across all ranks
            if event not in data_per_rank[output_rank_name]:
                data_per_rank[output_rank_name][event] = []
            data_per_rank[output_rank_name][event].append(event_data)

solver_event = None
if not args.no_compute_per_rank_statistics:
    print("Computing statistics per rank:")

    processed_data_per_rank = pd.DataFrame(columns=["rank", "event", "min", "max", "mean", "median", "std", "sum"])

    for rank in tqdm(data_per_rank, disable=args.no_tqdm):
        if args.io_only and io_only_event not in data_per_rank[rank]:
            continue

        for event in data_per_rank[rank]:
            data = np.vstack(data_per_rank[rank][event])
            processed_data_per_rank.loc[len(processed_data_per_rank.index)] = \
                [rank, event, np.min(data), np.max(data), np.mean(data), np.median(data), np.std(data), np.sum(data)]
            if "run" in event:
                if solver_event is None:
                    solver_event = event
                    solver_measurements = []
                assert(solver_event == event)
                solver_measurements.append(np.mean(data))

    if solver_event is not None:
        solver_load_imbalance = np.max(solver_measurements) / np.mean(solver_measurements)

    processed_data_per_rank.to_csv("processed_data_per_rank{}.csv".format("_io_only" if args.io_only else ""), index=False)


if not args.no_compute_per_step_statistics:
    print("Computing statistics per step:")

    processed_data_per_step = pd.DataFrame(columns=["step", "event", "min", "max", "mean", "median", "std", "sum"])

    max_number_of_timesteps = 0
    for rank in data_per_rank:
        for event in data_per_rank[rank]:
            max_number_of_timesteps = max(max_number_of_timesteps, len(data_per_rank[rank][event]))
    print("max. number of time steps: {}".format(max_number_of_timesteps))

    for step in tqdm(range(max_number_of_timesteps), disable=args.no_tqdm):
        data = {}
        for rank in data_per_rank:
            if args.io_only and io_only_event not in data_per_rank[rank]:
                continue

            for event in data_per_rank[rank]:
                if event not in data:
                    data[event] = np.array([])
                if step < len(data_per_rank[rank][event]):
                    data[event] = np.append(data[event], data_per_rank[rank][event][step])

        for event in data:
            if len(data[event]) > 0:
                processed_data_per_step.loc[len(processed_data_per_step.index)] = \
                    [step, event, np.min(data[event]), np.max(data[event]), np.mean(data[event]),
                     np.median(data[event]), np.std(data[event]), np.sum(data[event])]

    processed_data_per_step.to_csv("processed_data_per_step{}.csv".format("_io_only" if args.io_only else ""), index=False)

if not args.no_compute_total_statistics:
    print("Computing total statistics:")

    processed_data_total = pd.DataFrame(columns=["event", "min", "max", "mean", "median", "std", "sum"])

    data_per_event = {}
    for rank in tqdm(data_per_rank, disable=args.no_tqdm):
        if args.io_only and io_only_event not in data_per_rank[rank]:
            continue

        for event in data_per_rank[rank]:
            if event not in data_per_event:
                data_per_event[event] = np.array([])
            data_per_event[event] = np.append(data_per_event[event], np.vstack(data_per_rank[rank][event]))

    for event in tqdm(data_per_event, disable=args.no_tqdm):
        processed_data_total.loc[len(processed_data_total.index)] = \
            [event, np.min(data_per_event[event]), np.max(data_per_event[event]), np.mean(data_per_event[event]),
             np.median(data_per_event[event]), np.std(data_per_event[event]), np.sum(data_per_event[event])]
#        if "run" in event:
#            solver_event = event
#            solver_load_imbalance = np.max(data_per_event[event]) / np.mean(data_per_event[event])

    processed_data_total.to_csv("processed_data_per_total{}.csv".format("_io_only" if args.io_only else ""), index=False)
    # output data in PGF plots format
    pgf_plots_header_str = "numRanks,"
    pgf_plots_data_str = "{},".format(len(data_per_rank))
    for row_index, row in processed_data_total.iterrows():
        event_name = ""
        for col_name, col in row.items():
            if col_name == "event":
                event_name = col
            else:
                pgf_plots_header_str += "{}-{},".format(event_name, col_name)
                pgf_plots_data_str += str(col) + ","
    if solver_event is not None:
        pgf_plots_header_str += "{}-loadImbalance,".format(solver_event)
        pgf_plots_data_str += str(solver_load_imbalance) + ","
    print(pgf_plots_header_str[:-1].replace(" ", "").replace("/", ""))
    print(pgf_plots_data_str[:-1])
