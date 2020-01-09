import sys,os

sys.path.append('../tools/python')
sys.path.append('../python')
from ParIO import Parameters

class parameters_scan(Parameters):
    
    def __init__(self): pass

    def read_scanpars(self,path):
        self.path=path
        self.Read_Pars(path)
        
        self.scanpars=[]
        self.scanrange=[]
        self.stepwidth=[]
        #look for scan keyword
        for key,val in self.pardict.iteritems():
            if 'scan:' in val:
                if 'omn' in key:
                    self.scanpars.append('omn')
                else:
                    self.scanpars.append(key)
                scanbound=val.split('scan:')[1]
                sb=scanbound.split(',')
                if len(sb)==2:
                    self.scanrange.append([float(sb[0]),float(sb[1]),float(sb[1])-float(sb[0])])
                elif len(sb)==3:
                    self.scanrange.append([float(sb[0]),float(sb[2]),float(sb[2])-float(sb[0])])
                    self.stepwidth.append(float(sb[1]))

        self.scandim=len(self.scanpars)
        self.ddir=self.pardict['diagdir'].strip("'").split('/done/')[0]
        if self.pardict['read_checkpoint'] in ['.f.','.false.','.F.','.FALSE.','f','F']:
            self.use_checkpoints=False
        else:
            self.use_checkpoints=True

        if 'cat_output' in self.pardict:
            if self.pardict['cat_output'] in ['.t.','.true.','.T.','.TRUE.']:
                self.external_file=True
            else:
                self.external_file=False
        else:
            self.external_file=False

    def set_par_dirs(self,subdir):
        par=open(self.path)
        tpar='%s.tmp' % self.path
        temp=open(tpar,'w')
        for line in par:
            if 'par_in_dir' in line:
                line="par_in_dir='%s/in_par%s'\n" % (self.ddir,subdir)
            elif 'diagdir' in line:
                line="diagdir='%s/done%s'\n" % (self.ddir,subdir)
            temp.write(line)
        par.close()
        temp.close()
        os.rename(tpar,self.path)
