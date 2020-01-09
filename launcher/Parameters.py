from Tkinter import *
from Tooltips import *
from tkFileDialog   import askopenfilename,askdirectory,asksaveasfilename
from tkMessageBox import askokcancel
import os
import re

class Parameters:
    def __init__(self,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod

    def clearcomments(self,variable):
        regex=re.compile(r'\s*([-+\'\".,/a-zA-Z0-9_\s*]*)\s*!?\s*(.*)')
        result=regex.search(variable)
        if result and result.group(2)[:4]!='scan': 
            return regex.search(variable).group(1)
        else:
            vars.scan.set(1)
            return variable

    def Read_Pars(self):
        try: parin=askopenfilename(initialdir=launcher.startdir)
        except:
            launcher.Msg('Error opening parameters file.')
            return
        if not parin: return
        vars.clear()
        (path,parname)=os.path.split(parin)
        dum,jobpath=os.path.split(path)
        vars.jobpath.set(jobpath)
        parfile=open(parin,"r")
        #dictionary for parameters: {parameter: value}
        self.pardict={}
        #dictionary for recording the namelists the parameters belong to: {parameter: namelist}
        self.nmldict={}
        #counts species namelists
        countspec=0
        #Search file for parameters using regular expressions 
        for line in parfile:
            # Exclude commented lines
            if re.search(r'\s*!\w*\s*=.*',line)==None:
                # Check for and count species namelists
                if re.search(r'^\s*&(.*)',line):
                    #if namelist belongs to a species, append its number to the namelist
                    if re.search(r'^\s*&(.*)',line).group(1)=='species':
                        countspec+=1
                        self.nml=re.search(r'^\s*&(.*)',line).group(1)+str(countspec)
                        vars.listofnamelists.append('species'+str(countspec))
                    else: self.nml=re.search(r'^\s*&(.*)',line).group(1)
                # Search lines for <parameter> = <value> patterns
                p=re.compile(r'^\s*(\w*)\s*=\s*(.*)')
                m=p.search(line)
                # Pick matching lines and build parameter dictionary
                if m:
                    if m.group(1) in ('omn','omt','mass','charge','dens','temp','name','passive',
                                      'kappa_n','kappa_T','LT_center','Ln_center','LT_width','Ln_width',
                                      'prof_type','src_prof_type','src_amp','src_width','src_x0',
                                      'delta_x_n','delta_x_T'):
                        self.pardict[m.group(1)+str(countspec)]=m.group(2)
                        self.nmldict[m.group(1)+str(countspec)]=self.nml
                    else:
                        self.pardict[m.group(1)]=m.group(2)
                        self.nmldict[m.group(1)]=self.nml
        #clear the comments from all variables
        for item in self.pardict:
            self.pardict[item]=self.clearcomments(self.pardict[item])
        #allocate the parameters to the Tkinter variables
        self.Allocate_Parvals(self.pardict)
        #form is not empty anymore
        vars.form_cleared=0
        #after reading a parameters file, it has to be saved again before submit is possible
        vars.job_saved[vars.currentjob]=0
            
    def Allocate_Parvals(self,pars):
        # -------------------------------------------------------------------------------------------
        # this list will contain every parameter that is hard-wired and will appear in an entry field;
        # all other parameters will appear on the development frame
        listofpars=[]
                
        # this dictionary contains all parameters that are just numbers or strings and tells the launcher
        # to which entry field they should be allocated
        simpleparams={#inout
                    'diagdir': vars.diagdirV, 'chptdir': vars.chptdirV,
                    'istep_nrg': vars.istepnrgV,'istep_field': vars.istepfldV,'istep_mom': vars.istepmomV,
                    'istep_prof': vars.istepprofV,'istep_vsp': vars.istepvspV,'istep_schpt': vars.istepschptV,
                    'istep_omega': vars.istepomegaV,'istep_energy': vars.istepenergyV,
                    'istep_energy3d': vars.istepenergy3dV,
                    'istep_neoclass': vars.istepncV, 'iterdb_file': vars.iterdbfileV,
                    'iterdb_time': vars.iterdbtimeV,
                    
                    #box
                    'nx0': vars.nxV, 'nky0': vars.nyV, 'nw0': vars.nwV,'n_spec': vars.numspecV,
                    'nv0': vars.nvV, 'nz0': vars.nzV, 'kymin': vars.kyminV,
                    'lv': vars.lvV, 'lw': vars.lwV, 'lx': vars.lxV,
                    'x0': vars.x0V, 'ky0_ind': vars.ky0indV, 'kx_center': vars.kxcentV,

                    #general
                    'ntimesteps': vars.ntimeV, 'timelim': vars.timelimV, 'simtimelim': vars.simtimelimV,
                    'init_cond': vars.initcondV, 'init_aux_x': vars.auxxV, 'init_aux_y': vars.auxyV,
                    'init_aux_z': vars.auxzV, 'dt_max': vars.dtmaxV,
                    'timescheme': vars.timeschV, 'parscheme': vars.parschV, 'hyp_z': vars.hypzV,
                    'hyp_v': vars.hypvV, 'hyp_x': vars.hypxV, 'hyp_z_order': vars.hypzordV,
                    'hyp_x_order': vars.hypxordV, 'hyp_y': vars.hypyV, 'hyp_y_order': vars.hypyordV,
                    'ev_prec': vars.omprV,
                    'omega_prec': vars.omprV, 'overflow_limit': vars.oflV, 'underflow_limit': vars.uflV,
                    'beta': vars.betaV, 'courant': vars.courantV, 'coll': vars.collV,
                    'coll_cons_model':vars.collcmV, 'Zeff': vars.ZeffV,
                    'debye2': vars.debye2V, 'n_ev': vars.nevV, 'ev_max_it': vars.maxitV,
                    'ev_shift': vars.shiftV,
                    
                    #geometry
                    'magn_geometry': vars.geomtype, 'geomdir': vars.geomdirV, 'geomfile': vars.geomfileV,
                    'shat': vars.shatV, 'minor_r': vars.minorrV, 'major_R': vars.lengthV,
                    'parscale': vars.lengthV, 'q0': vars.q0V, 'trpeps': vars.trpepsV,
                    'amhd': vars.amhdV, 'kappa': vars.kappaV,'s_kappa': vars.s_kappaV,'delta': vars.deltaV,
                    's_delta': vars.s_deltaV,'zeta': vars.zetaV,'s_zeta': vars.s_zetaV,'drR': vars.drRV,'x_def': vars.fluxlabelV,
                    'flux_pos': vars.fluxposV, 'rref': vars.rrefV,
                    'q_coeffs': vars.q0V, 'rhostar': vars.rhostarV,

                    #nonlocal_x
                    'rad_bc_type': vars.radbctypeV, 'ck_heat': vars.ckheatV, 'l_buffer_size': vars.lbufsizeV,
                    'u_buffer_size': vars.ubufsizeV, 'lcoef_krook': vars.lcoefkrookV, 'ucoef_krook': vars.ucoefkrookV,
                    'lpow_krook': vars.lpowkrookV, 'upow_krook': vars.upowkrookV, 'reset_limit': vars.resetlimitV,

                    #external_contr
                    'ExBrate': vars.ExBrateV, 'pfsrate': vars.pfsrateV, 'ExB_stime': vars.ExBtimeV,
                    'kxind_phi_ext': vars.kxExBV, 'phi0_ext': vars.PhikxV,

                    #units
                    'Tref': vars.TrefV, 'nref': vars.nrefV, 'Bref': vars.BrefV, 'Lref': vars.LrefV, 'mref': vars.mrefV
                    }

        for item in simpleparams:
            if item in pars:
                simpleparams[item].set(pars[item].strip())
                listofpars.append(item)

        # ----------------------------------------------------------------------------
        # This dictionary contains boolean parameters. The 'value' part (after the colon) consists of a list, of which
        # the first item gives the variable in which the parameter will be saved, the second gives the default value,
        # and the third allows the treatment of inverted variables (e.g. pressure_off=.f. means vars.pressV=True in Launcher).
        # In this case, set the third item in the list to 1. 

        boolparams={'write_checkpoint': [vars.wrtchptV,1], 'nonlinear': [vars.nl,0], 'x_local': [vars.nonlocal,0,1],
                    'delzonal': [vars.delzonalV,0],# 'pressure_off': [vars.pressV,0,1],#'quasilinear_tm': [vars.quasilintmV,0],
                    'arakawa': [vars.arakawaV,1],'mag_prof': [vars.magprofV,1],#'spacediff_off': [vars.spacediffV,0,1],
                    'arakawa_zv': [vars.arakawazvV,1], 'shifted_metric': [vars.shiftedV,0],
                    'read_checkpoint': [vars.rdchptV,0],'chpt_h5': [vars.chpth5V,0], 'write_h5': [vars.wrth5V,0],
                    'write_std': [vars.wrtstdV,1], 'include_f0_contr':[vars.incf0V,0], 'del_phi':[vars.delphiV,0],
                    'coll_f_fm_on':[vars.collffmV,0], 'spacediff': [vars.spacediffV,0], 'pressure_term': [vars.pressV,0]}

        #recognize the settings of boolean parameters and transfer these to the appropriate Tkinter variables
        #specified above
        for item in boolparams:
            if item in pars:
                if pars[item].strip() in ('.f.','.F.','F','f','.false.'):
                    try:
                        if boolparams[item][2]:
                            boolparams[item][0].set(1)
                    except:
                        boolparams[item][0].set(0)
                if pars[item].strip() in ('.t.','.T.','T','t','.true.'):
                    try:
                        if boolparams[item][2]:
                            boolparams[item][0].set(0)
                    except:
                        boolparams[item][0].set(1)
            else:
                boolparams[item][0].set(boolparams[item][1])
            listofpars.append(item)
                       

        parallelizationpars={'n_procs_sim': [vars.nprocsV], 'n_procs_s': [vars.sV,vars.sautoV],
                             'n_procs_v': [vars.vV,vars.vautoV],'n_procs_w': [vars.wV,vars.wautoV],
                             'n_procs_x': [vars.xV,vars.xautoV], 'n_procs_y': [vars.yV,vars.yautoV],
                             'n_procs_z': [vars.zV,vars.zautoV], 'n_parallel_sims': [vars.nparsimV]}


        # if parameter is present, set entry field accordingly, if not => auto-parallelization
        for item in parallelizationpars:
            if item in pars:
                parallelizationpars[item][0].set(pars[item].strip())
                if int(pars[item].strip())>0:
                    try:
                        parallelizationpars[item][1].set(0)
                    except: pass
                else:
                    try:
                        parallelizationpars[item][1].set(1)
                    except: pass
            else:
                try: parallelizationpars[item][1].set(1)                                
                except: pass
            listofpars.append(item)
        launcher.autosw
        #if int(vars.nparsimV.get())!=1: vars.scan.set(1)

        if 'comp_type' in pars:
            if pars['comp_type'].strip() in ("'IV'","'iv'"):
                vars.ivev.set(1)
            elif pars['comp_type'].strip() in ("'EV'","'ev'"):
                vars.ivev.set(0)
            elif pars['comp_type'].strip() in ("'NC'","'nc'"):
                vars.ivev.set(2)
        else:
            vars.ivev.set(1)
        listofpars.append('comp_type')

        if 'nexc' in pars: 
            vars.nexcV.set(pars['nexc'].strip())
            listofpars.append('nexc')
        else: vars.nexcV.set('')
                
        if 'adapt_lx' in pars:
            if pars['adapt_lx'].strip() in ('.t.','.T.','T'):
                vars.adaptlx.set(1)
            if pars['adapt_lx'].strip() in ('.f.','.F.','F'):
                vars.adaptlx.set(0)
            listofpars.append('adapt_lx')
        else:
            if vars.nl.get()==0: vars.adaptlx.set(1)
            if vars.nl.get()==1: vars.adaptlx.set(0)

        if 'adapt_ly' in pars:
            if pars['adapt_ly'].strip() in ('.t.','.T.','T'):
                vars.adaptly.set(1)
            if pars['adapt_ly'].strip() in ('.f.','.F.','F'):
                vars.adaptly.set(0)
            listofpars.append('adapt_ly')
        else:
            if vars.nl.get()==0: vars.adaptly.set(1)
            if vars.nl.get()==1: vars.adaptly.set(0)

        if 'calc_dt' in pars:
            if pars['calc_dt'].strip() in ('.t.','.T.','T'):
                vars.calcdtV.set(1) 
            if pars['calc_dt'].strip() in ('.f.','.F.','F'):
                vars.calcdtV.set(0) 
            listofpars.append('calc_dt')
        else:
            if 'dt_max' in pars: vars.calcdtV.set(0)
            else: vars.calcdtV.set(1)

        if 'collision_op' in pars:
            listofpars.append('collision_op')
            vars.collisV.set(pars['collision_op'].strip())
            if pars['collision_op'].strip()=="'none'":
                vars.collonoffV.set(0)
            elif pars['collision_op'].strip() in ("'pitch-angle'","'landau'"):
                vars.collonoffV.set(1)
            else:
                vars.collisV.set("'landau'")
                vars.collonoffV.set(1)
            if not 'coll_cons_model' in pars:
                vars.collcmV.set("'xu_rosenbluth'")
            else: 
                vars.collcmV.set(pars['coll_cons_model'])

        #clear species arrays
        vars.omn=[]
        vars.omt=[]
        vars.mass=[]
        vars.charge=[]
        vars.temp=[]
        vars.dens=[]
        vars.specname=[]
        vars.passive=[]
        vars.fromfile=[]
        #the following ones are for nonlocal runs
        vars.proftype=[]
        vars.kappa_T=[]
        vars.kappa_n=[]
        vars.LT_center=[]
        vars.Ln_center=[]
        vars.LT_width=[]
        vars.Ln_width=[]
        vars.delta_x_T=[]
        vars.delta_x_n=[]
        vars.src_amp=[]
        vars.src_width=[]
        vars.src_x0=[]
        vars.src_proftype=[]
        nonlocalspecpars={'kappa_n': [vars.kappa_n,0.],'kappa_T': [vars.kappa_T,0.],'LT_center': [vars.LT_center,0.5],
                          'Ln_center': [vars.Ln_center,0.5],'LT_width': [vars.LT_width,0.1],'Ln_width': [vars.Ln_width,0.1],
                          'prof_type': [vars.proftype,2],'src_prof_type': [vars.src_proftype,0],
                          'src_amp': [vars.src_amp,0.],'src_width': [vars.src_width,0.1],'src_x0': [vars.src_x0,0.5],
                          'delta_x_n': [vars.delta_x_n,0.4],'delta_x_T': [vars.delta_x_T,0.4]}
        localspecpars={'omt': vars.omt, 'omn': vars.omn, 'mass': vars.mass, 'charge': vars.charge,
                       'temp': vars.temp, 'dens': vars.dens, 'name': vars.specname}

        #Build lists containing a single parameter for each species
        for i in xrange(1,10):
            if ('name'+str(i)) in pars and i<=int(vars.numspecV.get()):
                if vars.nonlocal.get():
                    for item in nonlocalspecpars:
                        if (item+str(i)) in pars:
                            nonlocalspecpars[item][0].append(pars[(item+str(i))].strip())
                        else:
                            nonlocalspecpars[item][0].append(nonlocalspecpars[item][1])
                        listofpars.append(item+str(i))
                    if ('omt'+str(i)) not in pars: vars.omt.append(0.)
                    if ('omn'+str(i)) not in pars: vars.omn.append(0.)
                else:
                    #in the local case, prof_type=-1,-2,-3 corresponds to file input, so we catch this here
                    if ('prof_type'+str(i)) in pars:
                        var=int(pars['prof_type'+str(i)])
                        if var in [-1,-2,-3]:
                            vars.proftype.append(var)
                        else:
                            vars.proftype.append(0)
                        for item in nonlocalspecpars:
                            if item!='prof_type'+str(i):
                                nonlocalspecpars[item][0].append(-100)
                        listofpars.append('prof_type'+str(i))
                    else:
                        vars.proftype.append(0)
                        for item in nonlocalspecpars:
                            nonlocalspecpars[item][0].append(-100)
                #take proftype from first species as we don't allow mixed proftypes
                vars.proftypeV.set(vars.proftype[0])
  
                for item in localspecpars:
                    if (item+str(i)) in pars:
                        localspecpars[item].append(pars[item+str(i)].strip())
                        listofpars.append(item+str(i))
                    else:
                        localspecpars[item].append(-1)
                #passive currently has its own treatment, since its boolean
                if ('passive'+str(i)) in pars:
                    if pars['passive'+str(i)].strip() in ['T','.t.','.T.']:
                        vars.passive.append(1)
                    elif pars['passive'+str(i)].strip() in ['F','.f.','.F.']:
                        vars.passive.append(0)
                        listofpars.append('passive'+str(i))
                else:
                    if ('name'+str(i)) in pars: vars.passive.append(0)
        #set number of species
        if 'n_spec' in pars: 
            vars.numspecV.set(pars['n_spec'].strip())
            listofpars.append('n_spec')
            launcher.ch_numspec()
                    
        #display and save the first species
        vars.ch_spec1(1)
        vars.save_spec1()
                
        #same for second species if present
        if int(vars.numspecV.get())>=2: 
            vars.ch_spec2(2)
            vars.save_spec2()


        #put all additional parameters into the development frame
        addparcount=0
        for par in pars:
            if par not in listofpars: 
                vars.addedpars[addparcount].set(par)
                vars.addedvalues[addparcount].set(pars[par])
                vars.namelist[addparcount].set(self.nmldict[par])
                addparcount+=1
        launcher.switch_nonlocal(weak=1)
        launcher.switch_lin_nonlin(weak=1)
        launcher.sel_geom()
        launcher.switch_ev_iv(weak=1)
        launcher.switch_collop()
        launcher.adaptlx_sel()
        launcher.calcdtclick()
        launcher.autosw()

    def set_output_parfile(self,mode):
        try: 
            if mode==0:
                if vars.jobpath.get()=='':
                    vars.parout=askdirectory(initialdir=launcher.startdir,title="Please select a 'prob' directory")
                    dum,jobpath=os.path.split(vars.parout)
                    vars.jobpath.set(jobpath)
                else: vars.parout=os.path.join(launcher.startdir,vars.jobpath.get())
                vars.parout=os.path.join(vars.parout,'parameters')
            elif mode==1: 
                vars.parout=asksaveasfilename(initialdir=launcher.startdir,title="Save parameters as...")
                dum1,dum2=os.path.split(vars.parout)
                dum3,jobpath=os.path.split(dum1)
                vars.jobpath.set(jobpath)
            elif mode==2:
                #vars.jobpath is already defined in this case
                vars.parout=os.path.join(launcher.startdir,vars.jobpath.get()+'/temp_pars')
            elif mode==3:
                #vars.jobpath is already defined in this case
                vars.parout=os.path.join(launcher.startdir,vars.jobpath.get()+'/parameters')
        except:
            return 

    def propagate_profile_data(self):
        if vars.proftypeV.get() in [-1,-2,-3]:
            strlist=[]
            if vars.proftoTrefV.get(): vars.TrefV.set(-1); strlist+=['Tref']
            if vars.proftonrefV.get(): vars.nrefV.set(-1); strlist+=['nref']
            if vars.proftobetaV.get(): vars.betaV.set(-1); strlist+=['beta']
            if vars.proftocollV.get(): vars.collV.set(-1); strlist+=['coll']
            if vars.proftodebyeV.get(): vars.debye2V.set(-1); strlist+=['debye2']
            if vars.proftorhostarV.get(): vars.rhostarV.set(-1); strlist+=['rhostar']
            if vars.proftypeV.get() in [-2,-3]:
                if vars.proftoExBV.get():
                    vars.ExBrateV.set(-1111); strlist+=['ExBrate']
                #set ExBrate explicitly to 0 if it is still set on auto to prevent accidental usage
                elif float(vars.ExBrateV.get())==-1111.: vars.ExBrateV.set(0)
                if vars.proftopfsV.get():
                    vars.pfsrateV.set(-1111); strlist+=['pfsrate']
                #ditto for pfs
                elif float(vars.pfsrateV.get())==-1111.: vars.pfsrateV.set(0)
                    
            if strlist:
                string=strlist[0]
                for item in range(1,len(strlist)):
                    string+=', '+strlist[item]
                launcher.Msg("Values of %s will be obtained from profile file input." %(string))
        



    #save a parameters file
    #mode=0: fixed filename 'parameters'
    #mode=1: user can choose filename
    #mode=2: write temporary parameters file
    #mode=3: write to 'parameters' when jobpath is already defined
    def Write_Pars(self,mode):
        if mode in [0,1]: self.Check_Pars()
        if vars.pars_okay>0: launcher.Msg('Warning: Found possible errors in parameters.')
        self.set_output_parfile(mode)
        if not vars.parout: return
        parfileout=open(vars.parout,"w")
        #routine that writes parameters only if they are specified
        def writeiftrue(name,var):
            if var!='':
                parfileout.write(name+var+'\n')

        #routine for writing parameters from the development tab
        def writeaddedpars(nml):
            for par in vars.addedpars:
                ind=vars.addedpars.index(par)
                x=par.get()
                y=vars.namelist[ind].get()
                vars.nmldict[str(x)]=y
                if x and vars.nmldict[str(x)]==nml:
                    parfileout.write(x.strip()+' = '+vars.addedvalues[ind].get().strip()+'\n')

        #routine that writes floating point parameters only if nonzero
        def writeifnonzero(name,var):
            try:
                if float(var)!=0.:
                    parfileout.write(name+var+'\n')
            except:
                #check if we have a scan if float() failed
                if 'scan' in var:
                    parfileout.write(name+var+'\n')
                
        
        #write parallelization namelist
        parfileout.write('&parallelization\n')
        parfileout.write('n_procs_s = '+vars.sV.get()+'\n')
        parfileout.write('n_procs_v = '+vars.vV.get()+'\n')
        parfileout.write('n_procs_w = '+vars.wV.get()+'\n')
        parfileout.write('n_procs_x = '+vars.xV.get()+'\n')
        parfileout.write('n_procs_y = '+vars.yV.get()+'\n')
        parfileout.write('n_procs_z = '+vars.zV.get()+'\n')
        if int(vars.nparsimV.get()) not in [0,1]: parfileout.write('n_parallel_sims = '+str(vars.nparsimV.get())+'\n')
        if vars.npsautoV.get()==0: parfileout.write('n_procs_sim = '+str(int(vars.nprocsV.get())/int(vars.nparsimV.get()))+'\n')
        writeaddedpars('parallelization')
        parfileout.write('/'+'\n\n')

        #write box namelist
        parfileout.write('&box\n')
        parfileout.write('n_spec = '+vars.numspecV.get()+'\n')
        parfileout.write('nx0 = '+vars.nxV.get()+'\n')
        parfileout.write('nky0 = '+vars.nyV.get()+'\n')
        parfileout.write('nz0 = '+vars.nzV.get()+'\n')
        parfileout.write('nv0 = '+vars.nvV.get()+'\n')
        parfileout.write('nw0 = '+vars.nwV.get()+'\n')
        parfileout.write('\n')
        parfileout.write('kymin = '+vars.kyminV.get()+'\n')
        if vars.ky0indV.get()!='1' and vars.nl.get()==0: writeiftrue('ky0_ind = ',vars.ky0indV.get())
        if vars.kxcentV.get()!='0.0': writeiftrue('kx_center = ',vars.kxcentV.get())
        if vars.write_lx: writeiftrue('lx = ',vars.lxV.get())
        if vars.adaptlx.get()==0 and vars.write_nexc:
            writeiftrue('nexc = ',vars.nexcV.get())
        parfileout.write('lv = '+vars.lvV.get()+'\n')
        parfileout.write('lw = '+vars.lwV.get()+'\n')
        if vars.nonlocal.get() or vars.geomtype.get() in ["'tracer_efit'","'chease'"]:
            parfileout.write('x0 = '+vars.x0V.get()+'\n')
        if vars.nl.get()!=1 and vars.adaptlx.get()==0: parfileout.write('adapt_lx = .F.\n')
        if vars.adaptly.get()==1: parfileout.write('adapt_ly = .T.\n')
        writeaddedpars('box')
        parfileout.write('/'+'\n\n')

        #write in_out namelist
        parfileout.write('&in_out\n')
        parfileout.write('diagdir = '+vars.diagdirV.get()+'\n')
        if vars.chptdirV.get().strip() not in ("''",""): parfileout.write('chptdir = '+vars.chptdirV.get()+'\n')
        
        if vars.rdchptV.get()==1: parfileout.write('read_checkpoint = .T.\n')
        else: parfileout.write('read_checkpoint = .F.\n')
        if vars.wrtchptV.get()==0: parfileout.write('write_checkpoint = .F.\n')
        if vars.chpth5V.get()==1: parfileout.write('chpt_h5 = .T.\n')
        if vars.wrth5V.get()==1: parfileout.write('write_h5 = .T.\n')
        if vars.wrtstdV.get()==0: parfileout.write('write_std = .F.\n')
#       else: parfileout.write('write_checkpoint = .F.\n')
        parfileout.write('istep_nrg = '+vars.istepnrgV.get()+'\n')
        if vars.istepomegaV.get()!='20': writeiftrue('istep_omega = ',vars.istepomegaV.get())
        parfileout.write('istep_field = '+vars.istepfldV.get()+'\n')
        parfileout.write('istep_mom = '+vars.istepmomV.get()+'\n')
        writeiftrue('istep_energy = ',vars.istepenergyV.get())
        writeiftrue('istep_energy3d = ',vars.istepenergy3dV.get())
        if vars.nonlocal.get(): writeiftrue('istep_prof = ',vars.istepprofV.get())
        parfileout.write('istep_vsp = '+vars.istepvspV.get()+'\n')
        writeiftrue('istep_neoclass = ',vars.istepncV.get())
        parfileout.write('istep_schpt = '+vars.istepschptV.get()+'\n')
        if vars.proftypeV.get() in [-2,-3]:
            parfileout.write('iterdb_file = '+vars.iterdbfileV.get()+'\n')
            writeiftrue('iterdb_time = ',vars.iterdbtimeV.get())
        writeaddedpars('in_out')
        parfileout.write('/'+'\n\n')


        #write general namelist
        parfileout.write('&general\n')
        if vars.nl.get()==1: parfileout.write('nonlinear = .T.\n')
        else: parfileout.write('nonlinear = .F.\n')

        if vars.nonlocal.get(): parfileout.write('x_local = .F.\n')
        if not vars.arakawazvV.get(): parfileout.write('arakawa_zv = .F.\n')
        
        if vars.ivev.get()==1: 
            parfileout.write("comp_type = 'IV'\n")
            if vars.calcdtV.get()==0: writeiftrue('dt_max = ',vars.dtmaxV.get())
            try:
                if float(vars.courantV.get())!=1.25:
                     writeiftrue('courant = ',vars.courantV.get())
            except:
                pass
            if vars.timeschV.get()!="'RK4'": parfileout.write('timescheme = '+vars.timeschV.get()+'\n')
        elif vars.ivev.get()==2:
            parfileout.write("comp_type = 'NC'\n")
        else: parfileout.write("comp_type = 'EV'\n")

        if vars.calcdtV.get()==1: parfileout.write('calc_dt = .T.\n')
        else: parfileout.write('calc_dt = .F.\n')

        if vars.pressV.get()==1: parfileout.write('pressure_term = .T.\n')
        if vars.delzonalV.get()==1: parfileout.write('delzonal = .T.\n')
#        if vars.quasilintmV.get()==1: parfileout.write('quasilinear_tm = .T.\n')
        
        if vars.ivev.get()==1: 
            writeiftrue('ntimesteps = ',vars.ntimeV.get())
        else: parfileout.write('ntimesteps = 0\n')

        parfileout.write('timelim = '+vars.timelimV.get()+'\n')
        writeiftrue('simtimelim = ',vars.simtimelimV.get())
        if vars.parschV.get()!="'c4th'": parfileout.write('parscheme = '+vars.parschV.get()+'\n')
        writeiftrue('n_ev = ',vars.nevV.get())
        if vars.ivev.get()==0 and vars.nl.get()==0:
            if vars.shiftV.get()!='(20,0)': writeiftrue('ev_shift = ',vars.shiftV.get())
            if vars.maxitV.get()!='0': writeiftrue('ev_max_it = ',vars.maxitV.get())
        if vars.omprV.get()!='0.001':
            if vars.ivev.get()==1:
                writeiftrue('omega_prec = ',vars.omprV.get())
            else:
                writeiftrue('ev_prec = ',vars.omprV.get())              
        if vars.oflV.get()!='1e-09': writeiftrue('overflow_limit = ',vars.oflV.get())
        if vars.uflV.get()!='1e-09': writeiftrue('underflow_limit = ',vars.uflV.get())
        
        if vars.incf0V.get()==1: parfileout.write('include_f0_contr = .T.\n')
        if vars.delphiV.get()==1: parfileout.write('del_phi = .T.\n')

        if vars.collonoffV.get()==1:
            parfileout.write('collision_op = '+vars.collisV.get()+'\n')
            parfileout.write('coll_cons_model = '+vars.collcmV.get()+'\n')
            if vars.collffmV.get()==1: parfileout.write('coll_f_fm_on = .T.\n')
            else: parfileout.write('coll_f_fm_on = .F.\n')
            #parameter is inverted
            if vars.spacediffV.get()==1: parfileout.write('spacediff = .T.\n')
            writeiftrue('coll = ',vars.collV.get())
            writeiftrue('Zeff = ',vars.ZeffV.get())
        writeiftrue('beta = ',vars.betaV.get())
        if vars.debye2V.get()!='0.0': writeiftrue('debye2 = ',vars.debye2V.get())
        if vars.hypxV.get()!='0.0': writeiftrue('hyp_x = ',vars.hypxV.get())
        if vars.hypxordV.get()!='4': writeiftrue('hyp_x_order = ',vars.hypxordV.get())
        if vars.hypyV.get()!='0.0': writeiftrue('hyp_y = ',vars.hypyV.get())
        if vars.hypyordV.get()!='4': writeiftrue('hyp_y_order = ',vars.hypyordV.get())
        writeiftrue('hyp_z = ',vars.hypzV.get())
        if vars.hypzordV.get()!='4': writeiftrue('hyp_z_order = ',vars.hypzordV.get())
        writeiftrue('hyp_v = ',vars.hypvV.get())
        parfileout.write('init_cond = '+vars.initcondV.get()+'\n')
        if vars.auxxV.get()!='-100': writeiftrue('init_aux_x = ',vars.auxxV.get())
        if vars.auxyV.get()!='-100': writeiftrue('init_aux_y = ',vars.auxyV.get())
        if vars.auxzV.get()!='-100': writeiftrue('init_aux_z = ',vars.auxzV.get())
        writeaddedpars('general')
        parfileout.write('/'+'\n\n')
        

        #write nonlocal_x namelist
        parfileout.write('&nonlocal_x\n')
        if vars.nonlocal.get():
            parfileout.write('rad_bc_type = '+vars.radbctypeV.get()+'\n')
            if not vars.arakawaV.get():
                parfileout.write('arakawa = .F.\n')
            if vars.shiftedV.get():
                parfileout.write('shifted_metric = .T.\n')
            if float(vars.ckheatV.get())>0.0: parfileout.write('ck_heat = '+vars.ckheatV.get()+'\n')
            parfileout.write('l_buffer_size = '+vars.lbufsizeV.get()+'\n')
            parfileout.write('lcoef_krook = '+vars.lcoefkrookV.get()+'\n')
            parfileout.write('lpow_krook = '+vars.lpowkrookV.get()+'\n')
            parfileout.write('u_buffer_size = '+vars.ubufsizeV.get()+'\n')
            parfileout.write('ucoef_krook = '+vars.ucoefkrookV.get()+'\n')
            parfileout.write('upow_krook = '+vars.upowkrookV.get()+'\n')
            if float(vars.resetlimitV.get())!=1000.:
                parfileout.write('reset_limit = '+vars.resetlimitV.get()+'\n')
        writeaddedpars('nonlocal_x')
        parfileout.write('/'+'\n\n')


        #write external_contr namelist
        parfileout.write('&external_contr\n')
        writeifnonzero('ExBrate = ',vars.ExBrateV.get())
        writeifnonzero('pfsrate = ',vars.pfsrateV.get())
        writeifnonzero('ExB_stime = ',vars.ExBtimeV.get())
        writeifnonzero('kxind_phi_ext = ',vars.kxExBV.get())
        writeifnonzero('phi0_ext = ',vars.PhikxV.get())
        writeaddedpars('external_contr')
        parfileout.write('/'+'\n\n')


        #write geometry namelist
        parfileout.write('&geometry\n')

        parfileout.write('magn_geometry = '+vars.geomtype.get()+'\n')
        if vars.geomtype.get()=="'tracer'":
            parfileout.write('geomdir = '+vars.geomdirV.get()+'\n')
            parfileout.write('geomfile = '+vars.geomfileV.get()+'\n')
        if vars.geomtype.get()=="'tracer_efit'":
            parfileout.write('geomdir = '+vars.geomdirV.get()+'\n')
            parfileout.write('geomfile = '+vars.geomfileV.get()+'\n')
        if vars.geomtype.get()=="'slab'":
            writeiftrue('parscale = ',vars.lengthV.get())
            writeiftrue('trpeps = ',vars.trpepsV.get())
            writeiftrue('shat = ',vars.shatV.get())
        if vars.geomtype.get() in ("'s_alpha'","'circular'"):
            if vars.nonlocal.get():
                parfileout.write('q_coeffs = '+vars.q0V.get()+'\n')
                writeiftrue('minor_r = ',vars.minorrV.get())
            else:
                writeiftrue('trpeps = ',vars.trpepsV.get())
                parfileout.write('q0 = '+vars.q0V.get()+'\n')
                writeiftrue('shat = ',vars.shatV.get())
            if vars.lengthV.get()!='1.0': writeiftrue('major_R = ',vars.lengthV.get())
        if vars.geomtype.get() in ["'s_alpha'","'miller'"]: writeiftrue('amhd = ',vars.amhdV.get())
        if vars.geomtype.get()=="'miller'":
            writeiftrue('trpeps = ',vars.trpepsV.get())
            parfileout.write('q0 = '+vars.q0V.get()+'\n')
            writeiftrue('shat = ',vars.shatV.get())
            writeiftrue('kappa = ',vars.kappaV.get())
            writeiftrue('s_kappa = ',vars.s_kappaV.get())
            writeiftrue('delta = ',vars.deltaV.get())
            writeiftrue('s_delta = ',vars.s_deltaV.get())
            writeiftrue('zeta = ',vars.zetaV.get())
            writeiftrue('s_zeta = ',vars.s_zetaV.get())
            writeiftrue('drR = ',vars.drRV.get())
            writeiftrue('minor_r = ',vars.minorrV.get())
            writeiftrue('major_R = ',vars.lengthV.get())
        if vars.geomtype.get()=="'chease'":
            parfileout.write('x_def = '+vars.fluxlabelV.get()+'\n')
            parfileout.write('flux_pos = '+vars.x0V.get()+'\n')
            if vars.geomdirV.get()!="''": parfileout.write('geomdir = '+vars.geomdirV.get()+'\n')
            if vars.geomfileV.get()!="''": parfileout.write('geomfile = '+vars.geomfileV.get()+'\n')
            if vars.nonlocal.get():
                parfileout.write('minor_r = '+vars.minorrV.get()+'\n')
            if vars.lengthV.get()!='1.0': writeiftrue('major_R = ',vars.lengthV.get())
        writeiftrue('rhostar = ', vars.rhostarV.get())
        if vars.nonlocal.get():
            if vars.magprofV.get(): parfileout.write('mag_prof = .T.\n')
        writeaddedpars('geometry')
        parfileout.write('/'+'\n\n')
            
        #write species namelists
        vars.save_spec1()
        if int(vars.numspecV.get())>1: vars.save_spec2()
        for i in range(int(vars.numspecV.get())):
            parfileout.write('&species\n')
            parfileout.write('name = '+vars.specname[i]+'\n')
            parfileout.write('mass = '+str(vars.mass[i])+'\n')
            parfileout.write('charge = '+str(vars.charge[i])+'\n')
            if vars.proftype[i]==0:
                parfileout.write('\n')
                parfileout.write('temp = '+str(vars.temp[i])+'\n')
                parfileout.write('omt = '+str(vars.omt[i])+'\n')
                parfileout.write('dens = '+str(vars.dens[i])+'\n')
                parfileout.write('omn = '+str(vars.omn[i])+'\n')
            elif vars.proftype[i] in range(1,6):
                parfileout.write('\n')
                parfileout.write('temp = '+str(vars.temp[i])+'\n')
                parfileout.write('kappa_T = '+str(vars.kappa_T[i])+'\n')
                parfileout.write('LT_center = '+str(vars.LT_center[i])+'\n')
                parfileout.write('LT_width = '+str(vars.LT_width[i])+'\n')
                if int(vars.proftypeV.get()) in (3,4):
                    parfileout.write('delta_x_T = '+str(vars.delta_x_T[i])+'\n')
                parfileout.write('\n')
                parfileout.write('dens = '+str(vars.dens[i])+'\n')
                parfileout.write('kappa_n = '+str(vars.kappa_n[i])+'\n')
                parfileout.write('Ln_center = '+str(vars.Ln_center[i])+'\n')
                parfileout.write('Ln_width = '+str(vars.Ln_width[i])+'\n')
                if int(vars.proftypeV.get()) in (3,4):
                     parfileout.write('delta_x_n = '+str(vars.delta_x_n[i])+'\n')
                parfileout.write('\n')
                parfileout.write('prof_type = '+str(vars.proftype[i])+'\n')
                parfileout.write('src_prof_type = '+str(vars.src_proftype[i])+'\n')
                if int(vars.src_proftype[i])!=0:
                    parfileout.write('src_amp = '+str(vars.src_amp[i])+'\n')
                    parfileout.write('src_width = '+str(vars.src_width[i])+'\n')
                    parfileout.write('src_x0 = '+str(vars.src_x0[i])+'\n')
            else:
                parfileout.write('prof_type = '+str(vars.proftype[i])+'\n')                    
            if vars.passive[i]: parfileout.write('passive = .T.\n')
            writeaddedpars('species'+str(i+1))
            parfileout.write('/'+'\n\n')
        vars.job_saved[vars.currentjob]=1

        #write units namelist
        for item in [vars.TrefV,vars.nrefV,vars.BrefV,vars.LrefV,vars.mrefV]:
            #if at least one reference value is given, write the namelist
            if item.get():
                parfileout.write('&units\n')
                if vars.TrefV.get():
                    parfileout.write('Tref = '+vars.TrefV.get()+'\n')
                if vars.nrefV.get():
                        parfileout.write('nref = '+vars.nrefV.get()+'\n')
                if vars.BrefV.get() and float(vars.BrefV.get())>0.:
                    parfileout.write('Bref = '+vars.BrefV.get()+'\n')
                if vars.LrefV.get() and float(vars.LrefV.get())>0.:
                    parfileout.write('Lref = '+vars.LrefV.get()+'\n')
                if vars.mrefV.get() and float(vars.mrefV.get())>0.:
                    parfileout.write('mref = '+vars.mrefV.get()+'\n')
                writeaddedpars('units')
                parfileout.write('/'+'\n\n')
                #break the for loop after values have been written
                break

        parfileout.close()
        self.Check_profile_input()
        launcher.Msg('Wrote parameters file to %s.' %(os.path.join(launcher.startdir,vars.jobpath.get())))

    def append_scan_nml(self):
        parfile=open(os.path.join(launcher.startdir,vars.jobpath.get()+
                                  '/parameters'),'a')
        def writeaddedpars(nml):
            for par in vars.addedpars:
                ind=vars.addedpars.index(par)
                x=par.get()
                y=vars.namelist[ind].get()
                vars.nmldict[str(x)]=y
                if x and vars.nmldict[str(x)]==nml:
                    parfile.write(x.strip()+' = '+
                    vars.addedvalues[ind].get().strip()+'\n')

        parfile.write('&scan\n')
        writeaddedpars('scan')
        parfile.write('/'+'\n\n')
        parfile.close()


    def Check_Pars(self):
        vars.pars_okay=0
        vars.save_spec1()
        if int(vars.numspecV.get())>1:
            vars.save_spec2()
        vars.Tkparlist=[vars.diagdirV,vars.nxV,vars.nyV,vars.nzV,vars.nvV,vars.nwV,vars.kyminV,vars.lvV,
                        vars.lwV,vars.dtmaxV,vars.q0V,vars.shatV,vars.trpepsV,vars.geomdirV,
                        vars.geomfileV,vars.x0V,vars.fluxlabelV,vars.lengthV]
        vars.parlist=['Output directory','Number of x points','Number of y modes','Number of z points',
                      'Number of v points','Number of mu points','Minimum ky','Extent of v space',
                      'Extent of mu space','Maximum timestep','Safety factor','Shear',
                      'Inverse aspect ratio','Geometry directory','Geometry file','Radial position',
                      'Flux label','Major radius']
        if vars.form_cleared==1:
            return
        launcher.MsgClear()
        checklist=[vars.diagdirV,vars.nxV,vars.nyV,vars.nzV,vars.nvV,vars.nwV,vars.kyminV,vars.lvV,vars.lwV]
        if vars.calcdtV.get()==0 and vars.ivev.get()==1:
            checklist.extend([vars.dtmaxV])
        if vars.geomtype.get() in ("'s_alpha'","'circular'"):
            checklist.extend([vars.q0V,vars.shatV])
        if vars.geomtype.get()=="'circular'":
            checklist.extend([vars.trpepsV])
        if vars.geomtype.get()=="'tracer'":
            checklist.extend([vars.geomdirV,vars.geomfileV,vars.q0V,vars.shatV])
        if vars.geomtype.get()=="'tracer_efit'":
            checklist.extend([vars.geomdirV,vars.geomfileV,vars.x0V])
        if vars.geomtype.get()=="'chease'":
            checklist.extend([vars.geomdirV,vars.geomfileV,vars.fluxlabelV,vars.lengthV,vars.x0V])
        if not vars.nonlocal.get():
            for item in checklist:
                if item.get() in ['',"''",'""']: 
                    launcher.Msg(vars.parlist[vars.Tkparlist.index(item)]+' is missing!')
                    vars.pars_okay=vars.pars_okay+1

        # Filter comments in n_procs lines, since these raise exceptions when doing the parallelization check
        regex=re.compile(r'\s*(-?\d*)\s*!?.*')
        if regex.search(vars.sV.get()):
            nps=int(regex.search(vars.sV.get()).group(1))
        else:
            nps=int(vars.sV.get())
        if regex.search(vars.vV.get()):
            npv=int(regex.search(vars.vV.get()).group(1))
        else:
            npv=int(vars.vV.get())
        if regex.search(vars.wV.get()):
            npw=int(regex.search(vars.wV.get()).group(1))
        else:
            npw=int(vars.wV.get())
        if regex.search(vars.xV.get()):
            npx=int(regex.search(vars.xV.get()).group(1))
        else:
            npx=int(vars.xV.get())
        if regex.search(vars.yV.get()):
            npy=int(regex.search(vars.yV.get()).group(1))
        else:
            npy=int(vars.yV.get())
        if regex.search(vars.zV.get()):
            npz=int(regex.search(vars.zV.get()).group(1))
        else:
            npz=int(vars.zV.get())
        #remove comments or !scan from grid node entries for parallelization checks
        if regex.search(vars.numspecV.get()):
            ns=int(regex.search(vars.numspecV.get()).group(1))
        else:
            ns=int(vars.numspecV.get())
        if regex.search(vars.nvV.get()):
            nv=int(regex.search(vars.nvV.get()).group(1))
        else:
            nv=int(vars.nvV.get())
        if regex.search(vars.nwV.get()):
            nw=int(regex.search(vars.nwV.get()).group(1))
        else:
            nw=int(vars.nwV.get())
        if regex.search(vars.nxV.get()):
            nx=int(regex.search(vars.nxV.get()).group(1))
        else:
            nx=int(vars.nxV.get())
        if regex.search(vars.nyV.get()):
            ny=int(regex.search(vars.nyV.get()).group(1))
        else:
            ny=int(vars.nyV.get())
        if regex.search(vars.nzV.get()):
            nz=int(regex.search(vars.nzV.get()).group(1))
        else:
            nz=int(vars.nzV.get())
        nparsims=int(vars.nparsimV.get())
        prod=nps*npv*npw*npx*npy*npz*nparsims
        # Check parallelization 
        if vars.sautoV.get()==0:
            if nps==0:
                launcher.Msg('Check species parallelization!')
                launcher.Opframe.sC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif ns%nps!=0:
                launcher.Msg('Check species parallelization!')
                launcher.Opframe.sC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif ns%nps==0 and ns/nps>1: launcher.Msg('Advice: Increase species parallelization.')
            else: launcher.Opframe.sC.config(bg=vars.stdbg)
        else: 
            launcher.Opframe.sC.config(bg=vars.stdbg)
        if vars.vautoV.get()==0:
            if npv==0:
                launcher.Msg('Check v parallelization!')
                launcher.Opframe.vC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif nv%npv!=0:
                launcher.MsgV('Check v parallelization!')
                launcher.Opframe.vC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif abs(nv/npv)<2 and nv%npv==0 and npv!=1:
                launcher.Msg('Warning: Such a high v parallelization is not recommended.')      
            else: launcher.Opframe.vC.config(bg=vars.stdbg)
        else: 
            launcher.Opframe.vC.config(bg=vars.stdbg)
        if vars.wautoV.get()==0:
            if npw==0:
                launcher.Msg('Check mu parallelization!')
                launcher.Opframe.wC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif nw%npw!=0:
                launcher.Msg('Check mu parallelization!')
                launcher.Opframe.wC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
#            elif vars.collonoffV.get()==1 and abs(nw/npw)<2 and nw%npw==0 and npw!=1:
#                launcher.Msg('Warning: With collisions, such a high mu parallelization is not recommended.') 
            else: launcher.Opframe.wC.config(bg=vars.stdbg)
        else: 
            launcher.Opframe.wC.config(bg=vars.stdbg)
        if vars.xautoV.get()==0:
            if npx==0:
                launcher.Msg('Check x parallelization!')
                launcher.Opframe.xC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif nx%npx!=0:
                launcher.Msg('Check x parallelization!')
                launcher.Opframe.xC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif npx>1:
                if not vars.nonlocal.get():
                    launcher.Msg('x parallelization not possible in local GENE version!')
                    launcher.Opframe.xC.config(bg="red")
                    vars.pars_okay=vars.pars_okay+1
            else: launcher.Opframe.xC.config(bg=vars.stdbg)
        else: 
            launcher.Opframe.xC.config(bg=vars.stdbg)
        if vars.yautoV.get()==0:
            if npy==0:
                launcher.Msg('Check y parallelization!')
                launcher.Opframe.yC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif ny%npy!=0:
                launcher.Msg('Check y parallelization!')
                launcher.Opframe.yC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            else: launcher.Opframe.yC.config(bg=vars.stdbg)
        else: 
            launcher.Opframe.yC.config(bg=vars.stdbg)
        if vars.zautoV.get()==0:
            if npz==0:
                launcher.Msg('Check z parallelization!')
                launcher.Opframe.zC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif nz%npz!=0:
                launcher.Msg('Check z parallelization!')
                launcher.Opframe.zC.config(bg="red")
                vars.pars_okay=vars.pars_okay+1
            elif abs(nz/npz)<2 and nz%npz==0 and npz!=1:
                launcher.Msg('Warning: Such a high z parallelization is not recommended.')      
            else: launcher.Opframe.zC.config(bg=vars.stdbg)
        else: 
            launcher.Opframe.zC.config(bg=vars.stdbg)

        if vars.nprocsV.get()==0 and prod==0:
            launcher.Msg('Total number of processes must be given for autoparallelization.')
            vars.pars_okay=vars.pars_okay+1
        if prod!=0:
            vars.nprocsV.set(prod)
            launcher.Msg('Note: Changed total number of processes to match parallelization settings.')


#       #check whether all integer numbers are really integer
#       listofints=[vars.nexcV,self,vars.istepnrgV,vars.istepfldV,vars.istepmomV,vars.istepschptV,vars.istepomegaV,vars.istepvspV,vars.hypxordV,vars.hypzordV,vars.ntimeV,vars.timelimV,vars.nevV,vars.maxitV]
#       for par in listofints:
#           try: int(par.get())
#           except:
#               launcher.Msg('Parameter '+str(par)+' must be given as an integer number.')
#               vars.pars_okay+=1
        
        if vars.nl.get()==0:
            if not vars.nonlocal.get():
                if ny>1: launcher.Msg('Warning: Using more than one ky mode for linear simulation.')
                if vars.adaptlx.get()==0: launcher.Msg('Warning: Radial box length adaptation is switched off for linear simulation.')
            if (not vars.istepomegaV.get() or int(vars.istepomegaV.get())==0):
                if not vars.incf0V.get() or int(vars.incf0V.get())==0:
                    vars.istepomegaV.set(20)
                    launcher.Msg('Switching on frequency diagnostics for linear simulation.')
        if vars.nl.get()==1:
            if ny==1: 
                launcher.Msg('Nonlinear simulation with only one ky mode not possible!')
                vars.pars_okay=vars.pars_okay+1

        if vars.geomtype.get()=="'tracer'":
            if not os.path.exists(os.path.join(vars.geomdirV.get().strip("'"),vars.geomfileV.get().strip("'"))):
                launcher.Msg('Geometry file does not exist!')
                vars.pars_okay+=1
            if vars.geomread:
                if int(vars.gridpoints)!=nz: 
                    launcher.Msg('Number of z points does not match geometry file!')
                    vars.pars_okay=vars.pars_okay+1
        if vars.geomtype.get()=="'circular'" and vars.trpepsV.get()=='0.0':
            launcher.Msg('Circular model does not work with zero inverse aspect ratio!')
            vars.pars_okay=vars.pars_okay+1
        #Check whether filename variables are enclosed in quotes and add those if not
        parlist=[vars.diagdirV,vars.chptdirV,vars.iterdbfileV,vars.geomdirV,vars.geomfileV]
        for item in parlist:
            if item.get():
                if item.get()[0] not in ("'",'"'): item.set("'"+item.get())
                if item.get()[-1] not in ("'",'"'): item.set(item.get()+"'")
        # Check whether diagdir exists
        self.set_output_parfile(0)
        if not os.path.isdir(vars.diagdirV.get().strip("'")):
            os.chdir(os.path.join(launcher.startdir,vars.jobpath.get()))
            if self.AskYesNo('Output directory', 'Output directory does not exist. Create?'):
                os.system("mkdir -p "+vars.diagdirV.get())
                if os.path.isdir(vars.diagdirV.get().strip("'")): 
                    launcher.Msg("Output directory created.")
                    isdir=1
                else:
                    launcher.Msg("Could not create directory; please check path and/or permissions.")
                    isdir=0
            else: 
                launcher.Msg("Output directory not created.")
                isdir=0
        else: 
            isdir=1
        if isdir:
            #check whether there is a previous run in the output dir
            filelist=[]
            if vars.wrth5V.get()==1:
                filelist+=['all_params.dat.h5','nrg.dat.h5','field.dat.h5']
            if vars.wrtstdV.get()==1:
                filelist+=['parameters.dat','nrg.dat','field.dat']
            if vars.nl.get()==1:
                count=0
                for file in filelist:
                    if os.path.exists(os.path.join(vars.diagdirV.get().strip("'"),file)):
                        count+=1
                if count>=3:
                    launcher.Msg('Warning: Old run data present in output directory.')
            if vars.ivev.get()==1:
                if vars.rdchptV.get()==1:
                    if vars.chpth5V.get()==1 and not os.path.exists(os.path.join(vars.diagdirV.get().strip("'"),'checkpoint.h5')):
                        launcher.Msg("Warning: No 'checkpoint.h5' file in output directory.")
                    elif not vars.chpth5V.get() and not os.path.exists(os.path.join(vars.diagdirV.get().strip("'"),'checkpoint')):
                        launcher.Msg("Warning: No 'checkpoint' file in output directory.")
        os.chdir(launcher.startdir)
        if int(vars.numspecV.get())>len(vars.specname):
            launcher.Msg(vars.numspecV.get()+' species requested, but only '+str(len(vars.specname))+' species saved.')
            vars.pars_okay+=1
        # Check quasineutrality, but only for the first numspec entries
        chargesum=0.
        adiabatic_electrons=1
        if ns>=2:
            failed=0
            for i in range(int(vars.numspecV.get())):
                try:
                    if not vars.passive[i]:
                        if float(vars.charge[i])==-1:
                            adiabatic_electrons=0
                        #check only for local profiles
                        if vars.proftypeV.get()==0:
                            chargesum=chargesum+float(vars.charge[i])*float(vars.dens[i])
                except: 
                    launcher.Msg('Warning: Error calculating quasineutrality.')
                    failed=1
                    break
            if not failed:
                if adiabatic_electrons and abs(chargesum-1.)<=1e-6:
                    pass
                elif not adiabatic_electrons and abs(chargesum+1.)<=1e-6:
                    pass
                elif not adiabatic_electrons and abs(chargesum)<=1e-6:
                    pass
                else:
                    launcher.Msg('Quasineutrality violated, please check densities and charges.')
                    vars.pars_okay=vars.pars_okay+1
        if ns==1 and float(vars.charge[0])==-1: adiabatic_electrons=0
        if vars.ivev.get()==1 and int(vars.ntimeV.get())==0:
            launcher.Msg('Warning: Number of timesteps for initial value calculation is zero!')
        self.propagate_profile_data()
        if adiabatic_electrons and float(vars.betaV.get())!=0.:
            launcher.Msg('Note: Setting beta=0 due to adiabatic electrons.')
            vars.betaV.set(0.)
        
        if vars.ivev.get()==2 and vars.incf0V.get()==0:
            launcher.Msg('Note: Switching on neoclassical term for NC equilibrium computation.')
            vars.incf0V.set(1)


    def Check_profile_input(self):
        error_1=0; error_2=0
        for i in range(int(vars.numspecV.get())):
            file='profiles_'+vars.specname[i].strip("'")
            if vars.proftypeV.get()==-1 and not os.path.exists(os.path.join(vars.jobpath.get(),file)):
                error_1+=1
            if vars.proftypeV.get() in[-2,-3] and not os.path.exists(vars.iterdbfileV.get().strip("'")):
                error_2+=1
        if error_1: launcher.Msg('Error: Profile file(s) missing in specified prob directory!')
        if error_2: launcher.Msg('Error: IterDB file not found!')
        return 
            
    def AskYesNo(self, title, message):
        return askokcancel( title, message )

