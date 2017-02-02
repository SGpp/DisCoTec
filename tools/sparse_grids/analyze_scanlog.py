import sys, os
class analyze_scanlog:
    def __init__(self,path):
        if os.path.exists(path+'/scan.log'):
            sfile=open(path+'/scan.log')
        else:
            sfile=open(path+'/neoclassic.log')

        firstline=sfile.readline()
        firstline=firstline.split('|')
        self.scanpars=[]
        for entry in firstline[1:]:
            if 'val1' in entry.strip(' '):
                break
            else:
                self.scanpars.append(entry.strip(' '))
        self.scandim=len(self.scanpars)
        self.coords=[]
        for ind in range(self.scandim):
            self.coords.append([])
        for line in sfile:
            vals=line.split('|')
            for ind in range(self.scandim):
                self.coords[ind].append(float(vals[ind+1]))
        sfile.close()
