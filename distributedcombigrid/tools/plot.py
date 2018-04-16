#!/usr/bin/env python2

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

def color_pool(data):
    """assigns a specific color to each event
    """
    color_cycler = cycle(color_list)
    color_map = {}
    for i in range(len(data)):
        rank = "rank" + str(i)
        for event in data[rank]["events"]:
            if event not in color_map:
                color_map[event] = next(color_cycler)
    return color_map

def overlap(i, j):
    return not(i[0] >= j[1] or i[1] <= j[0])

class Process:
    def __init__(self, rank, data):
        self.rank = rank
        data = data["rank" + str(rank)]
        self.intervals = []
        self.events = []
        for e in data["events"]:
            for i in data["events"][e]:
                # start, end, bottom, height, name
                self.intervals.append([i[0], i[1], 0, e])
                self.events.append(e)

        # build conflict graph
        # quadratic time complexity, can be improved with interval tree
        conflict = {}
        for i in range(len(self.intervals)):
            if not i in conflict:
                conflict[i] = []
            for j in range(i+1,len(self.intervals)):
                if overlap(self.intervals[i], self.intervals[j]):
                    if not j in conflict:
                        conflict[j] = []
                    conflict[i].append(j)
                    conflict[j].append(i)

        # graph coloring heuristic
        colors = {}
        maximum = 0
        for i in range(len(self.intervals)):
            c = 0
            jc = [colors[j] for j in conflict[i] if j in colors]
            while True:
                if len(jc) != 0:
                    minimum = min(jc)
                    if c < minimum:
                        break
                    elif c == minimum:
                        while len(jc) != 0 and min(jc) == minimum:
                            jc.remove(minimum)
                        c += 1
                else:
                    break
            colors[i] = c
            maximum = max(maximum, c)
        self.lanes = maximum + 1

        # assign lane to each interval
        for i in range(len(self.intervals)):
            self.intervals[i][2] =  colors[i]

    def timeline_plot(self, ax, colors, labels, offset):
        h = 1.0 / self.lanes
        for e in self.events:
            for l in range(self.lanes):
                bars = [(i[0], i[1]-i[0]) for i in self.intervals
                                          if i[2] == l and i[3] == e]
                if len(bars) > 0:
                    ax.broken_barh(bars, (offset + l*h, h), facecolor=colors[e],
                                   edgecolor="black", linewidth=1)
                    if e not in labels:
                        labels.append(e)

class ProcessGroup:
    def __init__(self, group_id, data):
        self.group_id = group_id
        self.group_worker = []
        for i in range(len(data)):
            r = "rank" + str(i)
            if int(data[r]["attributes"]["group"]) == group_id:
                t = Process(i, data)
                self.group_worker.append(t)
                if int(data[r]["attributes"]["group_manager"]) == 1:
                    self.group_manager = t

    def timeline_plot_all(self, ax, colors, labels, yticks, yticklables, offset):
        for i in range(len(self.group_worker)):
            self.group_worker[i].timeline_plot(ax, colors, labels, offset+i*1.1)
            yticks.append(offset+i*1.1)
            yticklables.append(self.group_worker[i].rank)
        offset += len(self.group_worker)*1.1 + 1.1
        return offset

    def timeline_plot_manager(self, ax, colors, labels, yticks, yticklables, offset):
        if self.group_manager != None:
            self.group_manager.timeline_plot(ax, colors, labels, offset)
            yticks.append(offset)
            yticklables.append(self.group_worker[i].rank)
            offset += 1.1
        return offset

def timeline_plot(data, groups, only_manager):
    """plots the timeline of the group or only the group manager
    """
    colors = color_pool(data)
    fig, ax = plt.subplots()
    ax.set_xlabel("time (s)")
    ax.set_ylabel('process rank')
    ax.grid(True)
    ax.set_axisbelow(True)
    yticklables = []
    yticks = []
    labels = []
    offset = 0
    for i in range(len(groups)):
        if only_manager == True:
            offset = groups[i].timeline_plot_manager(ax, colors, labels, yticks,
                                                     yticklables, offset)
        else:
            offset = groups[i].timeline_plot_all(ax, colors, labels, yticks,
                                                 yticklables, offset)
    ax.set_yticklabels(yticklables)
    ax.set_yticks(yticks)
    # proxy artist workaround for broken_barh labels with older matplotlib versions
    ax.legend([mpatch.Rectangle((0,0),1,1,fc=colors[i]) for i in labels],
              labels, loc=1).get_frame().set_alpha(0.75)
    # use seconds as unit
    scale_x = 1e6
    ticks_x = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x/scale_x))
    ax.xaxis.set_major_formatter(ticks_x)

def bar_plot_all(data, groups):
    """plots the average times of all processes or only the group managers
    """
    colors = color_pool(data)

    fig, ax = plt.subplots()
    ax.set_xlabel('global')
    ax.set_ylabel('time (s)')
    ax.grid(True)
    ax.set_axisbelow(True)

    ax.set_xticks([])

    times = {}
    for g in groups:
        for w in g.group_worker:
            for i in w.intervals:
                if i[3] not in times:
                    times[i[3]] = []
                else:
                    times[i[3]].append(i[1] - i[0])
    offset = 1
    for i in times:
        ax.bar(offset, np.mean(times[i]), 2, color=colors[i],
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

def bar_plot_groups(data, groups):
    """plots the average times of the process groups
    """
    colors = color_pool(data)
    labels = set()

    fig, ax = plt.subplots()
    ax.set_xlabel('process group')
    ax.set_ylabel('time (s)')
    ax.grid(True)
    ax.set_axisbelow(True)

    times = {}
    for g in groups:
        times[g.group_id] = {}
        for w in g.group_worker:
            for i in w.intervals:
                if i[3] not in times[g.group_id]:
                    times[g.group_id][i[3]] = []
                else:
                    times[g.group_id][i[3]].append(i[1] - i[0])

    xticks = []
    xlables = []
    offset = 0
    for t in times:
        xticks.append(offset - 1)
        xlables.append(t)
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
    data = json.load(open(sys.argv[1]))

    # get the number of groups
    group_count = 0
    for i in range(len(data)):
        r = "rank" + str(i)
        group_count = max(group_count, int(data[r]["attributes"]["group"]))
    group_count += 1

    groups = []
    for i in range(group_count):
        groups.append(ProcessGroup(i, data))

    print("Choose type of plot:")
    print("1 (timeline all processes),")
    print("2 (timeline group managers),")
    print("3 (average time all processes),")
    print("4 (average time process groups)")
    plot_type = int(input("\n"))

    if plot_type == 1:
        timeline_plot(data, groups, False)
    elif plot_type == 2:
        timeline_plot(data, groups, True)
    elif plot_type == 3:
        bar_plot_all(data, groups)
    elif plot_type == 4:
        bar_plot_groups(data, groups)
    else:
        raise ValueError("invalid type number")

    plt.tight_layout()
    plt.show()
except Exception as err:
    print("Error: " + str(err))
