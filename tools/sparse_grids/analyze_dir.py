import os, linecache, sys
from logfile import *

class analyze_dir:
    def __init__(self,in_spf,continuation,par=[],nsp=0,reanalyze=False):
        global nspec, evlog, timelog, nclog

        self.sims_per_folder=in_spf
        self.neoclassic = ('NC' in par.pardict['comp_type'])
        
        if par!=[]:
            nspec=int(par.pardict['n_spec'])
            #initialize the log files
            if self.neoclassic:
                nclog=logfile(par,'neoclassic.log',continuation)
            else:
                evlog=logfile(par,'scan.log',continuation)
                self.qllog=logfile(par,'quasilinear.log',continuation)
            timelog=logfile(par,'runtimes.log',continuation)
            self.log=True
            self.extfile=par.external_file
        else:
            nspec=nsp
            self.log=False
        self.reanalyze=reanalyze

    def get(self,ddir,parpoints,noldpoints,target=1):

        self.noldpoints=noldpoints

        if self.neoclassic:
            res=self.nc(ddir,parpoints)
        else:
            res=self.ql(ddir,parpoints,target)

        return res

    def ql(self,ddir,parpoints,target=1):

        if self.extfile:
            ga_q,ev=self.ga_q_ev(ddir,parpoints)
        else:
            ev=self.ev(ddir,parpoints)
            ga_q=self.ga_q(ddir)
            self.log_runtimes(ddir,parpoints)

	weight_exp=6

        qlmodel=[]
        for parind in range(len(ev)):
            if ev[parind]!=[]:
                #combine the transport values
                numerator=0.
                denominator=0.
                weight=0.
                for evind in range(len(ev[parind])):
                    if ev[parind][evind].real>0:
                        numerator+=ga_q[parind][evind]*ev[parind][evind].real**weight_exp
                        denominator+=ev[parind][evind].real**(weight_exp)
                        weight+=ev[parind][evind].real
                if target==0:
                    if denominator>0.:
                        qlmodel.append(weight*numerator/denominator)
                    else:
                        qlmodel.append(0.)
                elif target==1:
                    qlmodel.append(numerator)
                elif target==2:
                    qlmodel.append(weight)
                elif target==3:
                    qlmodel.append(denominator)
                elif target==4:
                    #find non-Hermitian degeneracies
                    qlmodel.append(abs(ev[parind][0]-ev[parind][1]))
                elif target==5:
                    #test
                    qlmodel.append(abs(ev[parind][0].real+ev[parind][1].real))
                elif target==6:
                    #test
                    qlmodel.append(ev[parind][1].imag)
            else:
                qlmodel.append([])

        if self.log:
            self.qllog.write(parpoints,qlmodel,self.noldpoints)
        return qlmodel

        #this is the growth rate of the zeroth mode
        #gr=[val[0].real for val in ev]
        #self.qllog.write(parpoints,gr,self.noldpoints)
        #return gr       

    def ev(self,ddir,parpoints):
        ev=[]

        if os.path.exists('%s/0k' % ddir):
            bind=0
            subdir='%dk'%bind
        else:
            bind=-1
            subdir=''

        pnum=1
        while os.path.exists('%s/%s' % (ddir,subdir)):
            filelist=os.listdir('%s/%s'%(ddir,subdir))
            if not self.reanalyze:
                pnum=1
            while 'parameters_%d' % pnum in filelist:
            #read eigenvalues                
                if 'eigenvalues_%d' % pnum in filelist:
                    evfile = open('%s/%s/eigenvalues_%d' % (ddir,subdir,pnum), 'r')
                #discard header
                    evfile.readline()
                #read eigenvalues
                    ev.append([])
                    for line in evfile:
                        strev=line.split()
                        ev[-1].append(float(strev[0])+float(strev[1])*1j)
                    evfile.close()
                else:
                    ev.append([])
                pnum+=1
            bind+=1
            subdir='%dk'%bind

        if ev==[]:
            print 'ERROR: Problem with GENE run, no results computed'
            sys.exit()

        if self.log:
            evlog.write(parpoints,ev,self.noldpoints)
        return ev
            
    def ga_q(self,ddir):
        ga_q=[]
        if os.path.exists('%s/0k' % ddir):
            bind=0
            subdir='%dk'%bind
        else:
            bind=-1
            subdir=''

        pnum=1
        while os.path.exists('%s/%s' % (ddir,subdir)):
            filelist=os.listdir('%s/%s'%(ddir,subdir))
            if not self.reanalyze:
                pnum=1
            while 'parameters_%d' % pnum in filelist:
                if 'nrg_%d' % pnum in filelist:
                    nrgfile = open('%s/%s/nrg_%d' % (ddir,subdir,pnum), 'r')
                    ga_q.append([])
                    linenum=0
                    for line in nrgfile:
                        linenum+=1

                    #use only the last species
                        if linenum%(nspec+1)==0:
                            strnrg=line.split()
                            ga_q[-1].append(float(strnrg[4])/float(strnrg[6]))

                    nrgfile.close()
                else:
                    ga_q.append([])
                pnum+=1
            bind+=1
            subdir='%dk'%bind

        return ga_q

    def log_runtimes(self,ddir,parpoints):
        i_t=[]

        if os.path.exists('%s/0k' % ddir):
            bind=0
            subdir='%dk'%bind
        else:
            bind=-1
            subdir=''

        pnum=1
        while os.path.exists('%s/%s' % (ddir,subdir)):
            filelist=os.listdir('%s/%s'%(ddir,subdir))
            if not self.reanalyze:
                pnum=1
            while 'parameters_%d' % pnum in filelist:
                parfile = open('%s/%s/parameters_%d' % (ddir,subdir,pnum), 'r')
                iterations=0
                time=0.
                for line in parfile:
                    if line.find('time for') >=0:
                        strti=line.split('=')
                        time=float(strti[1])
                    if line.find('iterations') >=0:
                        strti=line.split('=')
                        iterations=float(strti[1])
                i_t.append([iterations, time])
                parfile.close()
                pnum+=1
            bind+=1
            subdir='%dk'%bind

        timelog.write(parpoints,i_t,self.noldpoints)
        
    def ga_q_ev(self,ddir,parpoints):
        parsims=2
        ga_q=[]
        ev=[]
        i_t=[]
        for ind in range(len(parpoints)):
            ga_q.append([])
            ev.append([])
            i_t.append([])

        if os.path.exists('%s/0k' % ddir):
            bind=0
            subdir='%dk'%bind
        else:
            bind=0
            subdir=''
        while os.path.exists('%s/%s' % (ddir,subdir)):
            fileext=0
            while os.path.exists('%s/%s/scan_out_%d' % (ddir,subdir,fileext)):
                scanfile= open('%s/%s/scan_out_%d' % (ddir,subdir,fileext), 'r')
                for line in scanfile:
                    if 'parameters_' in line:
                        index=int(line.split('parameters_')[1])+bind*self.sims_per_folder
                        restype=1
                        nrgcount=0
                        continue
                    elif 'eigenvalues' in line:
                        restype=2
                        continue
                    elif '&parallelization' in line:
                        restype=3
                        continue
                    else:
                        if restype==1:
                            nrgcount=nrgcount+1
                            if nrgcount%(nspec+1)==0:
                                strnrg=line.split()
                                ga_q[index-1].append(float(strnrg[4])/float(strnrg[6]))
                        elif restype==2:
                            strev=line.split()
                            ev[index-1].append(float(strev[0])+float(strev[1])*1j)
                        elif restype==3:
                            if 'iterations' in line:
                                val=int(line.split('=')[1])
                                i_t[index-1].append(val)
                            elif 'time for' in line:
                                val=float(line.split('=')[1])
                                i_t[index-1].append(val)
                scanfile.close()
                fileext+=1
            bind+=1
            subdir='%dk'%bind
        if self.log:
            evlog.write(parpoints,ev,self.noldpoints)
            timelog.write(parpoints,i_t,self.noldpoints)            

        return ga_q,ev

    def nc(self,ddir,parpoints):
        nc=[]

        if os.path.exists('%s/0k' % ddir):
            bind=0
            subdir='%dk'%bind
        else:
            bind=-1
            subdir=''

        pnum=1
        while os.path.exists('%s/%s' % (ddir,subdir)):
            filelist=os.listdir('%s/%s'%(ddir,subdir))
            if not self.reanalyze:
                pnum=1
            while 'parameters_%d' % pnum in filelist:
            #read nc transport                
                if 'neoclass_%d' % pnum in filelist:
                    ncfile = open('%s/%s/neoclass_%d' % (ddir,subdir,pnum), 'r')
                #discard header
                    ncfile.readline()
                #read nc transport
                    line=ncfile.readline()
                    nc.append(float(line.split()[0]))
                    ncfile.close()
                else:
                    nc.append([])
                pnum+=1
            bind+=1
            subdir='%dk'%bind

        if nc==[]:
            print 'ERROR: Problem with GENE run, no results computed'
            sys.exit()

        if self.log:
            nclog.write(parpoints,nc,self.noldpoints)
        return nc
