import json
import sys
import os
import statistics
import getopt

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
options,events= getopt.getopt(argv[argi:],"ACDGHOm")
for i,j in options:
	if i=='-H':
		events.append("combine hierarchize")
	elif i=='-D':
		events.append("combine dehierarchize")
	elif i=='-O':
		events.append("only hierarchize")
		events.append("only dehierarchize")
	elif i=='-A':
		events.append("combine_allreduce")
	elif i=='-G':
		events.append("combine global reduce")
	elif i=='-C':
		events.append("combine")
	elif i=='-m':
		useMaxMax=True


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



dataset=[]
for f in filenames:
	dataset.append((f,json.load(open(f))))

totaltimes=[]
for f,data in dataset:
	totaltimes.append((f,(data["rank" +str(len(data)-1)]["events"]["total time"][0][1]-data["rank"+str(len(data)-1)]["events"]["total time"][0][0])/1000000))

print("Total times")
totaltimes.sort(key=lambda x:x[1])
for f,time in totaltimes:
	print(f,"  ",time,"s")

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

for e,fl in slist:
	print("\n",e)
	res=sorted(fl,key=lambda x:x[1])
	for f,time in res:
		print(f,"  ",time/1000,"ms")
