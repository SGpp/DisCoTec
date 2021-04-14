#!/usr/bin/env python3

"""Script for plotting the data gathered by the Stats class
"""

import sys
import os
import json
import matplotlib.patches as mpatch
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

# list of distinguishable colors, see:
# https://graphicdesign.stackexchange.com/revisions/3815/8
color_list = ["#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6",
              "#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762",
              "#004D43","#8FB0FF","#997D87","#5A0007","#809693","#FEFFE6",
              "#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A",
              "#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA",
              "#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018",
              "#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED",
              "#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C",
              "#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1",
              "#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459",
              "#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA",
              "#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329",
              "#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1",
              "#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C",
              "#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625",
              "#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534",
              "#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72",
              "#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6",
              "#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D",
              "#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5",
              "#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379",
              "#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6",
              "#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD",
              "#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD",
              "#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546",
              "#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87",
              "#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160",
              "#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA",
              "#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001",
              "#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749",
              "#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891",
              "#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190",
              "#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400",
              "#797868","#C6005A","#3B000A","#C86240","#29607C","#402334",
              "#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469",
              "#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9",
              "#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5",
              "#76912F","#003109","#0060CD","#D20096","#895563","#29201D",
              "#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85",
              "#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2",
              "#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D",
              "#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D",
              "#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88",
              "#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99",
              "#00CCFF","#674E60","#FC009C","#92896B"]

def color_pool(proc):
    """assigns a specific color to each event
    """
    color_cycler = cycle(color_list)
    color_map = {}
    for i in range(len(proc)):
        data = "rank" + str(i)
        for event in proc[data]["events"]:
            if event not in color_map:
                color_map[event] = next(color_cycler)
    return color_map

def timeline_plot_all(proc):
    """plots the timeline of all processes
    """
    colors = color_pool(proc)

    labels = []

    fig, ax = plt.subplots()

    ax.set_xlabel("time (s)")
    ax.set_ylabel('process rank')
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_yticklabels([str(i) for i in range(len(proc))])

    yticks = []
    ylim = 0
    offset = 1
    rank = 0
    try:
        for i in range(len(proc)):
            data = "rank" + str(i)
            group_offset = 3 * int(proc[data]["attributes"]["group"]) + offset
            ylim = group_offset + 3
            yticks.append(group_offset + 1)
            for j in proc[data]["events"]:
                bars = []
                for k in proc[data]["events"][j]:
                    bars.append((k[0], k[1]-k[0]))
                ax.broken_barh(bars, (group_offset, 2), facecolor=colors[j],
                               edgecolor="black", linewidth=1)
                if j not in labels:
                    labels.append(j)
            offset += 3
            rank = i
    except Exception as err:
        raise RuntimeError("rank " + str(rank) +
                           " is missing attribute " + str(err))
    ax.set_yticks(yticks)
    ax.set_ylim(0, ylim)

    # proxy artist workaround for broken_barh labels with older matplotlib versions
    ax.legend([mpatch.Rectangle((0,0),1,1,fc=colors[i]) for i in labels],
              labels, loc=1).get_frame().set_alpha(0.75)

    # use seconds as unit
    scale_x = 1e6
    ticks_x = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x/scale_x))
    ax.xaxis.set_major_formatter(ticks_x)

def timeline_plot_group_managers(proc):
    """plots the timeline of the group managers
    """
    colors = color_pool(proc)
    labels = []

    fig, ax = plt.subplots()
    ax.set_xlabel('time (s)')
    ax.set_ylabel('process rank')
    ax.grid(True)
    ax.set_axisbelow(True)

    ylables = []
    yticks = []
    ylim = 0
    offset = 1
    rank = 0
    try:
        for i in range(len(proc)):
            data = "rank" + str(i)
            group = int(proc[data]["attributes"]["group"])
            if bool(int(proc[data]["attributes"]["group_manager"])):
                ylim = offset + 3
                yticks.append(offset + 1)
                ylables.append(rank)
                for i in proc[data]["events"]:
                    bars = []
                    for j in proc[data]["events"][i]:
                        bars.append((j[0], j[1]-j[0]))
                    ax.broken_barh(bars, (offset, 2), facecolor=colors[i],
                                   edgecolor="black", linewidth=1)
                    if i not in labels:
                        labels.append(i)
                offset += 3
            rank = i
    except Exception as err:
        raise RuntimeError("rank " + str(rank) +
                           " is missing attribute " + str(err))
    ax.set_yticks(yticks)
    ax.set_ylim(0, ylim)
    ax.set_yticklabels(ylables)

    # proxy artist workaround for broken_barh labels with older matplotlib versions
    ax.legend([mpatch.Rectangle((0,0),1,1,fc=colors[i]) for i in labels],
              labels).get_frame().set_alpha(0.75)

    # use seconds as unit
    scale_x = 1e6
    ticks_x = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x/scale_x))
    ax.xaxis.set_major_formatter(ticks_x)

def bar_plot_all(proc, maxAverage):
    """plots the average times of all processes
    """
    colors = color_pool(proc)

    fig, ax = plt.subplots()
    ax.set_xlabel('global')
    ax.set_ylabel('time (s)')
    ax.grid(True)
    ax.set_axisbelow(True)

    ax.set_xticks([])

    times = {}
    timesSum = {}
    currentSum = {}
    for data in proc:
        for i in proc[data]["events"]:
            if i not in times:
                times[i] = []
                timesSum[i] = []
                currentSum[i] = 0
            for j in proc[data]["events"][i]:
                times[i].append(j[1] - j[0])
    
            timesSum[i].append(np.sum(times[i]) - currentSum[i])
            currentSum[i] = np.sum(times[i])
    offset = 1
    for i in times:
        value = 0
        if(maxAverage):
            value = np.max(timesSum[i])
        else:
            value = np.mean(times[i])
        ax.bar(offset, value , 2, color=colors[i],
               edgecolor="black", linewidth=1,
               yerr=np.std(times[i]), label=i,
               error_kw=dict(elinewidth=1,ecolor='black',
                             capsize=2,capthick=1))
        offset += 2
    ax.legend(loc=2).get_frame().set_alpha(0.75)

    ax.set_ylim(ymin=0)

    # use seconds as unit
    scale_y = 1e6
    ticks_y = ticker.FuncFormatter(lambda y, pos: "{0:g}".format(y/scale_y))
    ax.yaxis.set_major_formatter(ticks_y)

def bar_plot_group_managers(proc):
    """plots the average times of the process groups
    """
    colors = color_pool(proc)
    labels = set()

    fig, ax = plt.subplots()
    ax.set_xlabel('process group')
    ax.set_ylabel('time (s)')
    ax.grid(True)
    ax.set_axisbelow(True)

    times = {}
    rank = 0
    try:
        for i in range(len(proc)):
            data = "rank" + str(i)
            group = int(proc[data]["attributes"]["group"])
            if bool(int(proc[data]["attributes"]["group_manager"])):
                times[group] = {}
            rank += 1
        for data in proc:
            group = int(proc[data]["attributes"]["group"])
            for i in proc[data]["events"]:
                if i not in times[group]:
                    times[group][i] = []
                for j in proc[data]["events"][i]:
                    times[group][i].append(j[1] - j[0])
    except Exception as err:
        raise RuntimeError("rank " + str(rank) +
                           " is missing attribute " + str(err))

    xticks = []
    xlables = []
    offset = 0
    group = 0
    for t in times:
        xticks.append(offset + len(times[t]) - 1)
        xlables.append(group)
        for i in times[t]:
            ax.bar(offset, np.mean(times[t][i]), 2, color=colors[i],
                   edgecolor="black", linewidth=1,
                   yerr=np.std(times[t][i]),
                   label=i if i not in labels else "_nolegend_",
                   error_kw=dict(elinewidth=1,ecolor='black',
                                 capsize=2,capthick=1))
            labels.add(i)
            offset += 2
        offset += 8
        group += 1
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlables)
    ax.legend(loc=2).get_frame().set_alpha(0.75)

    ax.set_ylim(ymin=0)

    # use seconds as unit
    scale_y = 1e6
    ticks_y = ticker.FuncFormatter(lambda y, pos: "{0:g}".format(y/scale_y))
    ax.yaxis.set_major_formatter(ticks_y)

try:
    if len(sys.argv) == 1:
        raise RuntimeError("no input file specified")
    if len(sys.argv) > 2:
        raise RuntimeError("too many command line arguments")
    proc = json.load(open(sys.argv[1]))

    print("Choose type of plot:")
    print("1 (timeline all processes),")
    print("2 (timeline group managers),")
    print("3 (average time all processes),")
    print("4 (max total-times of all processes),")
    print("5 (average time process groups)")
    plot_type = int(input("\n"))

    if plot_type == 1:
        timeline_plot_all(proc)
    elif plot_type == 2:
        timeline_plot_group_managers(proc)
    elif plot_type == 3:
        bar_plot_all(proc,False)
    elif plot_type == 4:
        bar_plot_all(proc,True)
    elif plot_type == 5:
        bar_plot_group_managers(proc)
    else:
        raise ValueError("invalid type number")

    plt.tight_layout()
    plt.show()
except Exception as err:
    print("Error: " + str(err))
