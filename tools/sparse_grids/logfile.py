import os, sys

class logfile:
    
    def __init__(self,par_in,filename,continuation):

        if filename=='scan.log':
            #logfile for the eigenvalues: complex
            self.n_res=int(par_in.pardict['n_ev'])
            self.float_per_res=2
        elif filename=='runtimes.log':
            self.n_res=2
            self.float_per_res=1
        else:
            self.n_res=1
            self.float_per_res=1

        if continuation:
            self.lfile=open(par_in.ddir+'/'+filename,'a')
        else:
            self.lfile=open(par_in.ddir+'/'+filename,'w')
            self.lfile.write('#   nr ')
            for param in par_in.scanpars:    
                self.lfile.write('|%15s ' % param)

            if filename=='scan.log':
                for ind in range(self.n_res):
                    self.lfile.write('|%30s%d ' % ('val',ind+1))
            elif filename=='runtimes.log':
                self.lfile.write('|%15s ' % ('iterations'))
                self.lfile.write('|%15s ' % ('runtime'))
            else:
                for ind in range(self.n_res):
                    self.lfile.write('|%14s%d ' % ('val',ind+1))

            self.lfile.write('\n')

    def write(self,parpoints,vals,noldpoints):

        for parset in range(len(parpoints)):
            self.lfile.write('%06d ' % (noldpoints+1+parset))
            for param in parpoints[parset]:
                self.lfile.write('|%15.8g ' % param)
            if type(vals[parset])==list:
                for ind in range(self.n_res):
                    if ind<len(vals[parset]):
                        if self.float_per_res==2:
                            self.lfile.write('|%15.8g %15.8g ' % (vals[parset][ind].real,vals[parset][ind].imag))
                        else:
                            self.lfile.write('|%15.8g ' % vals[parset][ind])
                    else:
                        if self.float_per_res==2:
                            self.lfile.write('|%15s %15s ' % ('NaN','NaN'))
                        else:
#                            self.lfile.write('|%15s ' % 'NaN')
                            self.lfile.write('|%15.8g ' % 0.0)
                self.lfile.write('\n')
            else:
                self.lfile.write('|%15.8g \n' % vals[parset])
        self.lfile.flush()

    def __del__(self):
        self.lfile.close()
