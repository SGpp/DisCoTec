import os, sys, shutil

from analyze_dir import *
from checkpoints import *

class rungene:
    
    def __init__(self,par_in,syscall_in,continuation):
        global par, results, evlog, qllog
        global n_spec, in_par, syscall, genelogs
        par=par_in
        syscall=syscall_in

        n_spec=int(par.pardict['n_spec'])
        in_par = par.ddir+'/in_par'
        self.genefiles = par.ddir+'/genefiles'
        genelogs=self.genefiles+'/logs'
        self.done = par.ddir+'/done'
        self.chptfiles= par.ddir+'/checkpoints'
  
        if not continuation:
            dirlist=[in_par,self.genefiles,self.done,self.chptfiles,genelogs]
            for dir in dirlist:
                if os.path.exists(dir):
                    shutil.rmtree(dir)
                os.mkdir(dir)

        self.noldpars=0
        self.refnum=0
        self.kmax=-1

        #initialize some parameters file entries, create them if necessary       
        if par.use_checkpoints:
            par.nmldict['read_checkpoint']='in_out'
      
        if 'scan' not in par.namelists:
            par.namelists.append('scan')

        par.pardict['chpt_in']= "''"
        par.nmldict['chpt_in']='scan'

        par.pardict['par_in_dir']="'"+in_par+"'"
        par.nmldict['par_in_dir']='scan'

        par.Write_Pars('./parameters')
        shutil.copy('./parameters', '%s/parameters' % par.ddir)

        par.pardict['diagdir']="'"+self.done+"'"
        if 'chptdir' in par.pardict:
            del(par.pardict['chptdir'])

        self.sims_per_folder=1000

        #initialize the checkpoint and analysis routines
        if par.use_checkpoints:
            self.chpt=checkpoints(self.chptfiles,par)
        results=analyze_dir(self.sims_per_folder,continuation,par=par)

    def comp(self,parpoints):

        if os.path.exists(in_par):
            shutil.rmtree(in_par)
        os.mkdir(in_par)
        if os.path.exists(self.done):
            shutil.rmtree(self.done)
        os.mkdir(self.done)

        #create parameter files
        #split if parpoints list is too long
        nsets=(len(parpoints)-1)/self.sims_per_folder+1

        for setnum in range(nsets):
            lower_ind=setnum*self.sims_per_folder
            upper_ind=min(self.sims_per_folder,len(parpoints)-lower_ind)
            if nsets==1:
                subdir=''
            else:
                subdir='/%dk' % setnum
                os.mkdir(in_par+subdir)
                os.mkdir(self.done+subdir)
                
            par.pardict['diagdir']="'%s%s'" % (self.done,subdir)
            par.pardict['par_in_dir']="'%s%s'" % (in_par,subdir)

            for parind in range(0,upper_ind):
            #look for closest checkpoint
                if par.use_checkpoints:
                    clchpt=self.chpt.get_closest_checkpoint(parpoints[parind+lower_ind])
                else:
                    clchpt=0
                self.update_dict(parpoints[parind+lower_ind],clchpt)   
                par.Write_Pars('%s%s/parameters_%d' % (in_par,subdir,parind+1))

        #run scangene
        for setnum in range(nsets):
            if nsets==1:
                subdir=''
            else:
                subdir='%dk' % setnum
            par.set_par_dirs('/%s' % subdir)
            os.system('%s > %s/%s/gene_out'%(syscall, self.done,subdir))
            if nsets>1:
                open('%s/finished_%s' % (in_par,subdir), 'w').close() 

        #extract results
        quasilinear=results.get(self.done,parpoints,self.noldpars)

        #move files
        for setnum in range(nsets):
            if nsets==1:
                subdir=''
                os.rename('%s/gene_out'%self.done,'%s/gene_%d' % (genelogs,self.refnum))
            else:
                subdir='%dk' % setnum
                os.rename('%s/%s/gene_out'%(self.done,subdir),'%s/gene_%d_%s' % (genelogs,self.refnum,subdir))

            filelist=os.listdir('%s/%s'%(self.done,subdir))
            for file in filelist:
                fnstr=file.split('_')
                #number in the present refinement stage
                num=int(fnstr[-1])
                #global run number
                gnum=(setnum*self.sims_per_folder+num+self.noldpars)
                knum=(gnum-1)/self.sims_per_folder
                if knum>self.kmax:
                    #create new subdirectory if necessary
                    self.kmax=knum
                    if not par.external_file:
                        os.mkdir('%s/%dk' % (self.genefiles,self.kmax))
                    if par.use_checkpoints:
                        os.mkdir('%s/%dk' % (self.chptfiles,self.kmax))
                newfname=file.replace(fnstr[-1],str(gnum))
                if fnstr[0]=='checkpoint':
                    if par.use_checkpoints:
                        os.rename('%s/%s/%s' % (self.done,subdir,file), '%s/%dk/%s' % (self.chptfiles,knum,newfname))
                elif not par.external_file:
                    os.rename('%s/%s/%s' % (self.done,subdir,file), '%s/%dk/%s' % (self.genefiles,knum,newfname))

        if par.use_checkpoints:        
            self.chpt.update(parpoints)

        if par.external_file:
            os.rename(self.done,'%s/ref_%d' % (self.genefiles,self.refnum))
        
        self.noldpars+=len(quasilinear)
        self.refnum+=1
        return quasilinear


    def update_dict(self,parpoint,chpt):
        #updates the dictionary to write the parameters file
        ind=0
        for var in par.scanpars:
            if var in par.spec_nl:
                for spec in range(1, n_spec+1):
                    par.pardict[var+str(spec)]=str(parpoint[ind])
            else:
                par.pardict[var]=str(parpoint[ind])
            ind+=1

        if (chpt!=0):
            par.pardict['read_checkpoint']='.t.'
            par.pardict['chpt_in'] = "'%s/%dk/checkpoint_%s'" % (self.chptfiles,(chpt-1)/self.sims_per_folder,str(chpt)) 
        else:
            par.pardict['read_checkpoint']='.f.'


    #only for debugging    
    def print_scanpars(self):
        for var in par.scanpars:
            if var in par.spec_nl:
                for spec in range(1, n_spec+1):
                    print var+str(spec)+': '+par.pardict[var+str(spec)]
            else:
                print var+': '+par.pardict[var]

                
