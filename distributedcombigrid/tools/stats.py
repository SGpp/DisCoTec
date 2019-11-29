import json
import sys
import os
import statistics
import getopt
import numpy as np
import matplotlib.pyplot as plt

argv=sys.argv
argc=len(argv)

filenames=[]
#arg parsing
argi=1
while argi<argc and argv[argi][0]!="-":
	filenames.append(argv[argi])
	argi+=1

#print ("filenames",filenames)

useMaxMax=False
combine =False
plotstuff=False
specialevents=[]
sortingId=1
options,events= getopt.getopt(argv[argi:],"ACDEGHOmps")
for i,j in options:
	if i=='-H':
		events.append("combine hierarchize")
	elif i=='-D':
		events.append("combine dehierarchize")
	elif i=='-O':
		specialevents.append("only hierarchize")
		specialevents.append("only dehierarchize")
		onlyhier=True
	elif i=='-A':
		events.append("combine_allreduce")
	elif i=='-G':
		events.append("combine global reduce")
	elif i=='-C':
		#events.append("combine")
		combine=True
	elif i=='-m':
		useMaxMax=True
	elif i=='-s':
		sortingId=0#sorts after file name
	elif i=="-p":
		sortingId=0
		plotstuff=True
	elif i=='-E':
		specialevents.append("extSG")
		specialevents.append("addSG")
steps=10

# define functions
		 




def gather_stats(tagname,outerlist,durationlist,data):
	for i in range(len(data)-1):
		r = "rank" + str(i)
		outerlist.append(data[r]["events"][tagname])
		durationlist.append([])

	for i in range(len(outerlist)):
		for j in range(len(outerlist[i])):
			durationlist[i].append(outerlist[i][j][1]-outerlist[i][j][0])




		



def printstat(tagname,outlist,durlist):
	print(tagname)
	for i in range(len(outlist)):
		print(i,": ", statistics.mean(durlist[i])/1000,"ms\t",max(durlist[i])/1000,"ms")

def collectmax(outlist,durlist):
	tmp=[]
	for i in range(len(outlist)):
		tmp.append(max(durlist[i]))
	return max(tmp)

def collectavgmax(outlist,durlist):
	tmp=[]
	for j in range(len(durlist[0])):
		imax=0
		for i in range(len(durlist)):
			imax =max((imax,durlist[i][j]))
		tmp.append(imax)
	return statistics.mean(tmp)

def collectavgmaxspecial(outlist,durlist,steps_):
	tmp=[]
	dlist=[]
	for i in range(len(durlist)):
		tasksperp=len(durlist[i])//steps_
		chunks = [durlist[i][x:x+tasksperp] for x in range(0, len(durlist[i]), tasksperp)]
		if len(chunks) != steps_:
			print("Error:",len(chunks))
		dlist.append([sum(chunks[x]) for x in range(len(chunks))])
	

	for j in range(len(dlist[0])):
		imax=0
		for i in range(len(dlist)):
			imax =max((imax,dlist[i][j]))
		tmp.append(imax)
	return statistics.mean(tmp)






dataset=[]
for f in filenames:
	dataset.append((f,json.load(open(f))))

totaltimes=[]
combitime=[]
for f,data in dataset:
	totaltimes.append((f,(data["rank" +str(len(data)-1)]["events"]["total time"][0][1]-data["rank"+str(len(data)-1)]["events"]["total time"][0][0])/1000000))
	ctemp= [0]*steps
	for i in range (steps):
		ctemp[i]=(data["rank" +str(len(data)-1)]["events"]["combine"][i][1]-data["rank"+str(len(data)-1)]["events"]["combine"][i][0])/1000
	combitime.append((f,statistics.mean(ctemp)))




print("Total time")
totaltimes.sort(key=lambda x:x[sortingId])
for f,time in totaltimes:
	print(f,"  ",time,"s")


if(combine):
	print("\n combine:")
	combitime.sort(key=lambda x:x[sortingId])
	for f,time in combitime:
		print(f,"  ",time,"ms")
	

slist=[]
for e in events:
	flist=[]
	for f,data in dataset:
		ext1=[]
		ext2=[]
		gather_stats(e,ext1,ext2,data)
		if useMaxMax:
			flist.append((f,collectmax(ext1,ext2)))
		else:
			flist.append((f,collectavgmax(ext1,ext2)))
	slist.append((e,flist))

for e in specialevents:
	flist=[]
	for f,data in dataset:
		ext1=[]
		ext2=[]
		gather_stats(e,ext1,ext2,data)
		if useMaxMax:
			flist.append((f,collectmax(ext1,ext2)))
			print("maxmax not supported for",e)
		else:
			flist.append((f,collectavgmaxspecial(ext1,ext2,steps)))
	slist.append((e,flist))






for e,fl in slist:
	print("\n",e)
	res=sorted(fl,key=lambda x:x[sortingId])
	for f,time in res:
		print(f,"  ",time/1000,"ms")

if plotstuff:
	N=len(filenames)
	ind = np.arange(N)    # the x locations for the groups
	width = 0.35


	elist=list(zip(slist))
	
	plots=[]
	i=0
	rsum=np.zeros(N)
	for ex in elist:
		for (e,reslist) in ex:
			rs= np.array([r[1]/1000 for r in reslist])
			if i==0:
				plots.append(plt.bar(ind,rs))
			else:
				plots.append(plt.bar(ind,rs,bottom=rsum))
			rsum=rsum+rs
			i+=1
		
	plt.ylabel('Time in ms')
	#plt.title('Scores by group and gender')
	#plt.xticks(ind, filenames)
	#plt.legend([p[0] for p in plots],[e[0] for e in slist])
	plt.xticks(ind, ["3G 8N 14T","6G 4N 14T","12G 2N 14T","25G 1N 14T","11G 32N 1T","5G 64N 1T","22G 16N 1T"])
	abc=["Global Reduction","Hierarchization","Dehierarchization","Extraction from SG","Addition to SG"]
	plt.legend([p[0] for p in reversed(plots)],list(reversed(abc)))
	plt.show()

