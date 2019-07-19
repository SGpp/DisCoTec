import json
import sys
import os
import statistics


if len(sys.argv) == 1:
	print("using default file")
	filename="distributedcombigrid/examples/shared_example/timers.json"
else:
	filename=sys.argv[1]

data = json.load(open(filename))
htimes=[]
durations=[]
dtimes=[]
ddurations=[]
o_htimes=[]
o_durations=[]
o_dtimes=[]
o_ddurations=[]
for i in range(len(data)-1):
	r = "rank" + str(i)
	htimes.append(data[r]["events"]["combine hierarchize"])
	durations.append([])
	dtimes.append(data[r]["events"]["combine dehierarchize"])
	ddurations.append([])
	o_htimes.append(data[r]["events"]["only hierarchize"])
	o_durations.append([])
	o_dtimes.append(data[r]["events"]["only dehierarchize"])
	o_ddurations.append([])

for i in range(len(htimes)):
	for j in range(len(htimes[i])):
		durations[i].append(htimes[i][j][1]-htimes[i][j][0])
	for j in range(len(dtimes[i])):
		ddurations[i].append(dtimes[i][j][1]-dtimes[i][j][0])
	for j in range(len(o_htimes[i])):
		o_durations[i].append(o_htimes[i][j][1]-o_htimes[i][j][0])
	for j in range(len(o_dtimes[i])):
		o_ddurations[i].append(o_dtimes[i][j][1]-o_dtimes[i][j][0])

print("total time",(data["rank" +str(len(data)-1)]["events"]["total time"][0][1]-data["rank"+str(len(data)-1)]["events"]["total time"][0][0])/1000000,"s")

print("Hierachize")
for i in range(len(htimes)):
	print(i,": ", statistics.mean(durations[i])/1000,"ms\t",max(durations[i])/1000,"ms")

print("Dehierachize")
for i in range(len(dtimes)):
	print(i,": ", statistics.mean(ddurations[i])/1000,"ms\t",max(ddurations[i])/1000,"ms")

print("Only Hierachize")
for i in range(len(o_htimes)):
	print(i,": ", sum(o_durations[i])/(1000*len(htimes[i])),"ms\t")

print("Only Dehierachize")
for i in range(len(o_dtimes)):
	print(i,": ", sum(o_ddurations[i])/(1000*len(dtimes[i])),"ms\t")
