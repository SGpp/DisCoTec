
import re
import numpy as np

class AutoparFile(object):
    def __init__(self):
        self.data=[]
        self.paral_label=[]
        self.paral_index = 0
        self.perf_vec_index = 0
        self.perf_vec_label = {}

    def read(self,filename):
        afile=open(filename,"rt")
        skip_next_line=True

        for line in afile:
            if skip_next_line:
                # skip first line
                skip_next_line=False
                continue
            mo=re.search(r"parallelization:\s*([0-9 ]+?)\s+nblocks:\s*\d+",line)
            if mo:
                if self.paral_index>=1 and len(self.data[-1])==0:
                    # last parallelization had no valid perf_vecs (due to memory)
                    # remove the whole parallelization
                    self.paral_label.pop()
                    self.paral_index -= 1
                    self.data.pop()

                # now add the new parallelization
                self.paral_label.append(mo.group(1))
                self.paral_index += 1
                # append a new empty block to the data list
                self.data.append([])
                continue
    
            # this is a perf_vec line
            mo=re.search(r"\s*([0-9 ]+?)\s*:\s*([0-9.]+)\s+MB, t =\s*([0-9.E+-]+)\s*sec",line)
            if mo:
                this_perf_vec_label=mo.group(1)
                if this_perf_vec_label not in self.perf_vec_label:
                    self.perf_vec_label.setdefault(this_perf_vec_label,self.perf_vec_index)
                    self.perf_vec_index += 1

                # append to the last block (the actual one) a new performance run
                self.data[-1].append({'perf_vec_index':self.perf_vec_label[this_perf_vec_label],'WCT': float(mo.group(3)),'memory':float(mo.group(2))})
                continue

            mo=re.search(r"\s*best perf_vec",line)
            if mo:
                skip_next_line=True
                continue
    
            mo=re.search(r"\s*([0-9 ]+?)\s*:\s*skipped due to estimated memory requirement of\s*([0-9.]+)\s+MB per core",line)
            if mo:
                continue

            # nothing matches, then exit
            #print line,
            continue

        afile.close()
        if (len(self.data[-1])==0):
            print "Last element is empty, it will be removed."
            del self.data[-1]
            self.paral_index -= 1

        #templist = [x for x in self.data if len(x)!=0]
        #self.data=templist
        

        print "We have ",len(self.data)," parallelizations."
        print "We have ",self.perf_vec_index," different perf_vecs."

        self.data_array=np.zeros((self.paral_index,self.perf_vec_index),order='F')
        # initialize the data_array with None
        for icol in range(self.perf_vec_index):
            for irow in range(self.paral_index):
                self.data_array[irow,icol]=None

        self.parall_array  = np.zeros_like(self.data_array)
        self.perfvec_array = np.zeros_like(self.data_array)
        self.mem_array     = np.zeros_like(self.data_array)

        # the next two loops are used to build a numpy array
        # which contains the grid coordinates (parallelization index and
        # perf_vec index
        # loop over all parallelizations
        for iparall in range(self.paral_index):
            # for one parallelization loop over the perf_vecs
            for pv in self.data[iparall]:
                self.data_array[iparall,pv['perf_vec_index']]=pv['WCT']
                self.mem_array[iparall,pv['perf_vec_index']]=pv['memory']

        # loop over all perf_vec indices
        for pv in range(self.perf_vec_index):
            self.parall_array[iparall,pv]=iparall
            self.perfvec_array[iparall,pv]=pv


    def getBestParallelization(self):
        min_time = 100000.0
        for (par,plabel) in zip(self.data,self.paral_label):
            for pv in par:
                if pv['WCT']<min_time:
                    min_time=pv['WCT']
                    min_pv = pv['perf_vec_index']
                    min_par = plabel

        for pv in self.perf_vec_label:
            print pv

        for (key,val) in self.perf_vec_label.iteritems():
            if min_pv==val:
                perf_vec = key

        return (min_par,perf_vec)
