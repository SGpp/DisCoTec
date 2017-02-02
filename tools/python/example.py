#!/usr/bin/env python
from ParIO import Parameters

#initialize parameters class instance
par=Parameters()

#at the moment: two ways to get parameters dictionary (1. as return value of Read_Pars 2. as dictionary within par: par.pardict)
par.Read_Pars('./parameters')

#dictionary which allocates parameters to namelists
#nml=par.nmldict

#for item in par.pardict.keys(): print item,par.pardict[item],par.nmldict[item]
#print par.namelists

par.Write_Pars('./testpars')
