from Tkinter import *

class Multijob:

#initialize variables for multiple-job management
    def __init__(self,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        self.init_multi()

    def init_multi(self):
        self.mt_geomread=[]
        self.mt_currentjob=[]
        self.mt_jobpath=[]
        self.mt_numspec=[]
        self.mt_omn=[]
        self.mt_omt=[]
        self.mt_mass=[]
        self.mt_charge=[]
        self.mt_temp=[]
        self.mt_dens=[]
        self.mt_specname=[]
        self.mt_passive=[]
        self.mt_kappa_T=[]
        self.mt_kappa_n=[]
        self.mt_LT_center=[]
        self.mt_Ln_center=[]
        self.mt_LT_width=[]
        self.mt_Ln_width=[]
        self.mt_delta_x_T=[]
        self.mt_delta_x_n=[]
        self.mt_proftype=[]
        self.mt_src_proftype=[]
        self.mt_src_amp=[]
        self.mt_src_width=[]
        self.mt_src_x0=[]
        
        self.mt_nl=[] 
        self.mt_nonlocal=[] 
        self.mt_scan=[]
        self.mt_adaptlx=[]
        self.mt_adaptly=[]
        self.mt_timesch=[]
        self.mt_parsch=[]
        self.mt_s=[]
        self.mt_v=[]
        self.mt_w=[]
        self.mt_x=[]
        self.mt_y=[]
        self.mt_z=[]
        self.mt_nparsim=[]
        self.mt_sauto=[]
        self.mt_vauto=[]
        self.mt_wauto=[]
        self.mt_xauto=[]
        self.mt_yauto=[]
        self.mt_zauto=[]
        self.mt_npsauto=[]
        self.mt_nprocs=[]
        self.mt_hypz=[]
        self.mt_hypzord=[]
        self.mt_hypx=[]
        self.mt_hypxord=[]
        self.mt_hypy=[]
        self.mt_hypyord=[]
        self.mt_hypv=[]
        self.mt_arakawazv=[]
        self.mt_rdchpt=[]
        self.mt_wrtchpt=[]
        self.mt_chpth5=[]
        self.mt_wrth5=[]
        self.mt_wrtstd=[]
        self.mt_diagdir=[]
        self.mt_chptdir=[]
        self.mt_istepnrg=[]
        self.mt_istepomega=[]
        self.mt_istepfld=[]
        self.mt_istepmom=[]
        self.mt_istepenergy=[]
        self.mt_istepenergy3d=[]
        self.mt_istepprof=[]
        self.mt_istepneoclass=[]
        self.mt_istepvsp=[]
        self.mt_istepschpt=[]
        self.mt_iterdbfile=[]
        self.mt_iterdbtime=[]
        self.mt_nx=[]
        self.mt_ny=[]
        self.mt_nz=[]
        self.mt_nv=[]
        self.mt_nw=[]
        self.mt_kymin=[]
        self.mt_lx=[]
        self.mt_lv=[]
        self.mt_lw=[]
        self.mt_x0=[]
        self.mt_ky0ind=[]
        self.mt_nexc=[]
        self.mt_kxcent=[]
        self.mt_initcond=[]
        self.mt_auxx=[]
        self.mt_auxy=[]
        self.mt_auxz=[] 
        self.mt_ntime=[]
        self.mt_timelim=[]
        self.mt_simtimelim=[]
        self.mt_dtmax=[]
        self.mt_calcdt=[]
        self.mt_ompr=[]
        self.mt_ofl=[]
        self.mt_ufl=[]
        self.mt_beta=[]
        self.mt_courant=[]
        self.mt_incf0=[]
        self.mt_delphi=[]
        self.mt_coll=[]
        self.mt_Zeff=[]
        self.mt_collis=[]
        self.mt_collonoff=[]
        self.mt_collffm=[]
        self.mt_collcm=[]
        self.mt_spacediff=[]
        self.mt_debye2=[]
        self.mt_ivev=[]
        self.mt_whichev=[]
        self.mt_nev=[]
        self.mt_maxit=[]
        self.mt_shift=[]
        self.mt_press=[]
        self.mt_delzonal=[]
#        self.mt_quasilintm=[]
        self.mt_ExBrate=[]
        self.mt_pfsrate=[]
        self.mt_ExBtime=[]
        self.mt_kxExB=[]
        self.mt_Phikx=[]        
        self.mt_geomtype=[]
        self.mt_minorr=[]
        self.mt_length=[]
        self.mt_shat=[]
        self.mt_trpeps=[]
        self.mt_q0=[]
        self.mt_amhd=[]
        self.mt_kappa=[]
        self.mt_s_kappa=[]
        self.mt_delta=[]
        self.mt_s_delta=[]
        self.mt_zeta=[]
        self.mt_s_zeta=[]
        self.mt_drR=[]
        self.mt_geomdir=[]
        self.mt_geomfile=[]
        self.mt_fluxlabel=[]
        self.mt_fluxpos=[]
        self.mt_rref=[] 
        self.mt_rhostar=[]
        self.mt_magprof=[]
        self.mt_x0=[]
        self.mt_radbctype=[]
        self.mt_lcoefkrook=[]
        self.mt_lbufsize=[]
        self.mt_lpow=[]
        self.mt_ucoefkrook=[]
        self.mt_ubufsize=[]
        self.mt_upow=[]
        self.mt_arakawa=[]
        self.mt_ckheat=[]
        self.mt_resetlimit=[]
        self.mt_shifted=[]
        self.mt_addedpars=[]
        self.mt_addedvalues=[]
        self.mt_namelist=[]
        self.mt_Tref=[]
        self.mt_nref=[]
        self.mt_Lref=[]
        self.mt_Bref=[]
        self.mt_mref=[]
        self.mt_proftoTref=[]
        self.mt_proftonref=[]
        self.mt_proftobeta=[]
        self.mt_proftocoll=[]
        self.mt_proftodebye=[]
        self.mt_proftoExB=[]
        self.mt_proftopfs=[]

        
    #check whether the current job has been saved before and set self.index accordingly
    def mt_find(self,a,b):
        try: ind=a.index(b)
        except: ind=-1
        if ind>=0:
            self.index=ind
            self.append=0
        else:
            self.append=1

    #saves variables for multiple-job management
    def mt_save(self,a,b): 
        if not self.append:
            try: a[self.index]=b
            except: self.Msg("Error saving into multi-job array.")
        else: a.append(b)
             
    #this routine saves an array of Tkinter variables as is needed for the entries in the development frame
    def mt_saveTkarray(self,a,b): 
        if not self.append:
            try: 
                for m in range(len(b)):
                    a[self.index][m]=b[m].get()
            except: self.Msg("Error saving into multi-job array.")
        else: 
            c=[]
            for m in range(len(b)):
                c.append(b[m].get())
            a.append(c)

    #loads variables from another job
    def mt_load(self,a):
        return a[self.index]

    #loads Tk variables from another job
    def mt_loadTk(self,a,b):
        try:    b.set(a[self.index])
        except: pass
    
    def multiswitch(self,event,i):
        #Button with number i lights up
        vars.save_spec1()
        if int(vars.numspecV.get())>=2: vars.save_spec2()
        for j in range(len(launcher.multibutt)):
            if j==i: launcher.multibutt[j].config(bg='#eee')
            else: launcher.multibutt[j].config(bg=vars.stdbg)
        #upon changing to different job number:
        if i!=vars.currentjob:
            #check length of multitasking lists (mt_currentjob as indicator)
            #print len(self.mt_currentjob)
            #is the current job present twice? (debugging)
            #print self.mt_currentjob.count(self.currentjob)

            #save Tk variables into normal python lists
            self.append=0
            self.index=-1
            self.mt_find(self.mt_currentjob,vars.currentjob)
            self.mt_save(self.mt_currentjob,vars.currentjob)
            self.mt_save(self.mt_jobpath,vars.jobpath.get())
            self.mt_save(self.mt_geomread,vars.geomread)
            self.mt_save(self.mt_numspec,vars.numspecV.get())
            self.mt_save(self.mt_omn,vars.omn)
            self.mt_save(self.mt_omt,vars.omt)
            self.mt_save(self.mt_mass,vars.mass)
            self.mt_save(self.mt_charge,vars.charge)
            self.mt_save(self.mt_temp,vars.temp)
            self.mt_save(self.mt_dens,vars.dens)
            self.mt_save(self.mt_passive,vars.passive)
            self.mt_save(self.mt_specname,vars.specname)
            self.mt_save(self.mt_kappa_T,vars.kappa_T)
            self.mt_save(self.mt_kappa_n,vars.kappa_n)
            self.mt_save(self.mt_LT_center,vars.LT_center)
            self.mt_save(self.mt_Ln_center,vars.Ln_center)
            self.mt_save(self.mt_LT_width,vars.LT_width)
            self.mt_save(self.mt_Ln_width,vars.Ln_width)
            self.mt_save(self.mt_delta_x_T,vars.delta_x_T)
            self.mt_save(self.mt_delta_x_n,vars.delta_x_n)
            self.mt_save(self.mt_proftype,vars.proftype)
            self.mt_save(self.mt_src_proftype,vars.src_proftype)
            self.mt_save(self.mt_src_amp,vars.src_amp)
            self.mt_save(self.mt_src_width,vars.src_width)
            self.mt_save(self.mt_src_x0,vars.src_x0)
            
            self.mt_save(self.mt_nl,vars.nl.get())
            self.mt_save(self.mt_nonlocal,vars.nonlocal.get())
            self.mt_save(self.mt_scan,vars.scan.get())
            self.mt_save(self.mt_adaptlx,vars.adaptlx.get())
            self.mt_save(self.mt_adaptly,vars.adaptly.get())
            self.mt_save(self.mt_timesch,vars.timeschV.get())
            self.mt_save(self.mt_parsch,vars.parschV.get())
            self.mt_save(self.mt_s,vars.sV.get())
            self.mt_save(self.mt_v,vars.vV.get())
            self.mt_save(self.mt_w,vars.wV.get())
            self.mt_save(self.mt_x,vars.xV.get())
            self.mt_save(self.mt_y,vars.yV.get())
            self.mt_save(self.mt_z,vars.zV.get())
            self.mt_save(self.mt_nparsim,vars.nparsimV.get())
            self.mt_save(self.mt_sauto,vars.sautoV.get())
            self.mt_save(self.mt_vauto,vars.vautoV.get())
            self.mt_save(self.mt_wauto,vars.wautoV.get())
            self.mt_save(self.mt_xauto,vars.xautoV.get())
            self.mt_save(self.mt_yauto,vars.yautoV.get())
            self.mt_save(self.mt_zauto,vars.zautoV.get())
            self.mt_save(self.mt_npsauto,vars.npsautoV.get())
            self.mt_save(self.mt_nprocs,vars.nprocsV.get())
            self.mt_save(self.mt_hypz,vars.hypzV.get())
            self.mt_save(self.mt_hypzord,vars.hypzordV.get())
            self.mt_save(self.mt_hypx,vars.hypxV.get())
            self.mt_save(self.mt_hypxord,vars.hypxordV.get())
            self.mt_save(self.mt_hypy,vars.hypyV.get())
            self.mt_save(self.mt_hypyord,vars.hypyordV.get())
            self.mt_save(self.mt_hypv,vars.hypvV.get())
            self.mt_save(self.mt_arakawazv,vars.arakawazvV.get())
            self.mt_save(self.mt_rdchpt,vars.rdchptV.get())
            self.mt_save(self.mt_wrtchpt,vars.wrtchptV.get())
            self.mt_save(self.mt_chpth5,vars.chpth5V.get())
            self.mt_save(self.mt_wrth5,vars.wrth5V.get())
            self.mt_save(self.mt_wrtstd,vars.wrtstdV.get())
            self.mt_save(self.mt_diagdir,vars.diagdirV.get())
            self.mt_save(self.mt_chptdir,vars.chptdirV.get())
            self.mt_save(self.mt_istepnrg,vars.istepnrgV.get())
            self.mt_save(self.mt_istepomega,vars.istepomegaV.get())
            self.mt_save(self.mt_istepfld,vars.istepfldV.get())
            self.mt_save(self.mt_istepmom,vars.istepmomV.get())
            self.mt_save(self.mt_istepenergy,vars.istepenergyV.get())
            self.mt_save(self.mt_istepenergy3d,vars.istepenergy3dV.get())
            self.mt_save(self.mt_istepprof,vars.istepprofV.get())
            self.mt_save(self.mt_istepvsp,vars.istepvspV.get())
            self.mt_save(self.mt_istepneoclass,vars.istepncV.get())
            self.mt_save(self.mt_istepschpt,vars.istepschptV.get())
            self.mt_save(self.mt_iterdbfile,vars.iterdbfileV.get())
            self.mt_save(self.mt_iterdbtime,vars.iterdbtimeV.get())
            self.mt_save(self.mt_nx,vars.nxV.get())
            self.mt_save(self.mt_ny,vars.nyV.get())
            self.mt_save(self.mt_nz,vars.nzV.get())
            self.mt_save(self.mt_nv,vars.nvV.get())
            self.mt_save(self.mt_nw,vars.nwV.get())
            self.mt_save(self.mt_kymin,vars.kyminV.get())
            self.mt_save(self.mt_lx,vars.lxV.get())
            self.mt_save(self.mt_lv,vars.lvV.get())
            self.mt_save(self.mt_lw,vars.lwV.get())
            self.mt_save(self.mt_x0,vars.x0V.get())
            self.mt_save(self.mt_ky0ind,vars.ky0indV.get())
            self.mt_save(self.mt_nexc,vars.nexcV.get())
            self.mt_save(self.mt_kxcent,vars.kxcentV.get())
            self.mt_save(self.mt_initcond,vars.initcondV.get())
            self.mt_save(self.mt_auxx,vars.auxxV.get())
            self.mt_save(self.mt_auxy,vars.auxyV.get())
            self.mt_save(self.mt_auxz,vars.auxzV.get()) 
            self.mt_save(self.mt_ntime,vars.ntimeV.get())
            self.mt_save(self.mt_timelim,vars.timelimV.get())
            self.mt_save(self.mt_simtimelim,vars.simtimelimV.get())
            self.mt_save(self.mt_dtmax,vars.dtmaxV.get())
            self.mt_save(self.mt_calcdt,vars.calcdtV.get())
            self.mt_save(self.mt_ompr,vars.omprV.get())
            self.mt_save(self.mt_ofl,vars.oflV.get())
            self.mt_save(self.mt_ufl,vars.uflV.get())
            self.mt_save(self.mt_beta,vars.betaV.get())
            self.mt_save(self.mt_courant,vars.courantV.get())
            self.mt_save(self.mt_incf0,vars.incf0V.get())
            self.mt_save(self.mt_delphi,vars.delphiV.get())
            self.mt_save(self.mt_coll,vars.collV.get())
            self.mt_save(self.mt_Zeff,vars.ZeffV.get())
            self.mt_save(self.mt_collis,vars.collisV.get())
            self.mt_save(self.mt_collonoff,vars.collonoffV.get())
            self.mt_save(self.mt_collffm,vars.collffmV.get())
            self.mt_save(self.mt_collcm,vars.collcmV.get())
            self.mt_save(self.mt_spacediff,vars.spacediffV.get())
            self.mt_save(self.mt_debye2,vars.debye2V.get())
            self.mt_save(self.mt_ivev,vars.ivev.get())
            self.mt_save(self.mt_whichev,vars.whichevV.get())
            self.mt_save(self.mt_nev,vars.nevV.get())
            self.mt_save(self.mt_maxit,vars.maxitV.get())
            self.mt_save(self.mt_shift,vars.shiftV.get())
            self.mt_save(self.mt_press,vars.pressV.get())
            self.mt_save(self.mt_ExBrate,vars.ExBrateV.get())
            self.mt_save(self.mt_pfsrate,vars.pfsrateV.get())
            self.mt_save(self.mt_ExBtime,vars.ExBtimeV.get())
            self.mt_save(self.mt_kxExB,vars.kxExBV.get())
            self.mt_save(self.mt_Phikx,vars.PhikxV.get())
            self.mt_save(self.mt_delzonal,vars.delzonalV.get())
#            self.mt_save(self.mt_quasilintm,vars.quasilintmV.get())
            self.mt_save(self.mt_geomtype,vars.geomtype.get())
            self.mt_save(self.mt_minorr,vars.minorrV.get())
            self.mt_save(self.mt_length,vars.lengthV.get())
            self.mt_save(self.mt_shat,vars.shatV.get())
            self.mt_save(self.mt_trpeps,vars.trpepsV.get())
            self.mt_save(self.mt_q0,vars.q0V.get())
            self.mt_save(self.mt_amhd,vars.amhdV.get())
            self.mt_save(self.mt_kappa,vars.kappaV.get())
            self.mt_save(self.mt_s_kappa,vars.s_kappaV.get())
            self.mt_save(self.mt_delta,vars.deltaV.get())
            self.mt_save(self.mt_s_delta,vars.s_deltaV.get())
            self.mt_save(self.mt_zeta,vars.zetaV.get())
            self.mt_save(self.mt_s_zeta,vars.s_zetaV.get())
            self.mt_save(self.mt_drR,vars.drRV.get())
            self.mt_save(self.mt_geomdir,vars.geomdirV.get())
            self.mt_save(self.mt_geomfile,vars.geomfileV.get())
            self.mt_save(self.mt_fluxlabel,vars.fluxlabelV.get())
            self.mt_save(self.mt_fluxpos,vars.fluxposV.get())
            self.mt_save(self.mt_rref,vars.rrefV.get())
            self.mt_save(self.mt_rhostar,vars.rhostarV.get())
            self.mt_save(self.mt_magprof,vars.magprofV.get())
            self.mt_save(self.mt_radbctype,vars.radbctypeV.get())
            self.mt_save(self.mt_lcoefkrook,vars.lcoefkrookV.get())
            self.mt_save(self.mt_lbufsize,vars.lbufsizeV.get())
            self.mt_save(self.mt_lpow,vars.lpowkrookV.get())
            self.mt_save(self.mt_ucoefkrook,vars.ucoefkrookV.get())
            self.mt_save(self.mt_ubufsize,vars.ubufsizeV.get())
            self.mt_save(self.mt_upow,vars.upowkrookV.get())
            self.mt_save(self.mt_arakawa,vars.arakawaV.get())
            self.mt_save(self.mt_ckheat,vars.ckheatV.get())
            self.mt_save(self.mt_resetlimit,vars.resetlimitV.get())
            self.mt_save(self.mt_shifted,vars.shiftedV.get())
            self.mt_saveTkarray(self.mt_addedpars,vars.addedpars)           
            self.mt_saveTkarray(self.mt_addedvalues,vars.addedvalues)
            self.mt_saveTkarray(self.mt_namelist,vars.namelist)
            self.mt_save(self.mt_Tref,vars.TrefV.get())
            self.mt_save(self.mt_nref,vars.nrefV.get())
            self.mt_save(self.mt_Lref,vars.LrefV.get())
            self.mt_save(self.mt_Bref,vars.BrefV.get())
            self.mt_save(self.mt_mref,vars.mrefV.get())
            self.mt_save(self.mt_proftoTref,vars.proftoTrefV.get())
            self.mt_save(self.mt_proftonref,vars.proftonrefV.get())
            self.mt_save(self.mt_proftobeta,vars.proftobetaV.get())
            self.mt_save(self.mt_proftocoll,vars.proftocollV.get())
            self.mt_save(self.mt_proftodebye,vars.proftodebyeV.get())
            self.mt_save(self.mt_proftoExB,vars.proftoExBV.get())
            self.mt_save(self.mt_proftopfs,vars.proftopfsV.get())
            self.append=0
            vars.currentjob=i
            vars.job_saved[i]=0
            self.index=-1
            self.mt_find(self.mt_currentjob,vars.currentjob)
        if self.index!=-1 and not self.append:
            vars.geomread=self.mt_load(self.mt_geomread)
            self.mt_loadTk(self.mt_numspec,vars.numspecV)
            vars.kappa_T=self.mt_load(self.mt_kappa_T)
            vars.kappa_n=self.mt_load(self.mt_kappa_n)
            vars.LT_center=self.mt_load(self.mt_LT_center)
            vars.Ln_center=self.mt_load(self.mt_Ln_center)
            vars.LT_width=self.mt_load(self.mt_LT_width)
            vars.Ln_width=self.mt_load(self.mt_Ln_width)
            vars.delta_x_T=self.mt_load(self.mt_delta_x_T)
            vars.delta_x_n=self.mt_load(self.mt_delta_x_n)
            vars.proftype=self.mt_load(self.mt_proftype)
            vars.src_proftype=self.mt_load(self.mt_src_proftype)
            vars.src_amp=self.mt_load(self.mt_src_amp)
            vars.src_width=self.mt_load(self.mt_src_width)
            vars.src_x0=self.mt_load(self.mt_src_x0)
            vars.omn=self.mt_load(self.mt_omn)
            vars.omt=self.mt_load(self.mt_omt)
            vars.mass=self.mt_load(self.mt_mass)
            vars.charge=self.mt_load(self.mt_charge)
            vars.temp=self.mt_load(self.mt_temp)
            vars.dens=self.mt_load(self.mt_dens)
            vars.passive=self.mt_load(self.mt_passive)
            vars.specname=self.mt_load(self.mt_specname)
            for k in range(len(vars.EntryDevPar)):
                vars.addedpars[k].set(self.mt_addedpars[self.index][k])
                vars.addedvalues[k].set(self.mt_addedvalues[self.index][k])
                vars.namelist[k].set(self.mt_namelist[self.index][k])
            self.mt_loadTk(self.mt_nl,vars.nl)
            self.mt_loadTk(self.mt_nonlocal,vars.nonlocal)
            self.mt_loadTk(self.mt_scan,vars.scan)
            self.mt_loadTk(self.mt_adaptlx,vars.adaptlx)
            self.mt_loadTk(self.mt_adaptly,vars.adaptly)
            self.mt_loadTk(self.mt_jobpath,vars.jobpath)
            self.mt_loadTk(self.mt_timesch,vars.timeschV)
            self.mt_loadTk(self.mt_parsch,vars.parschV)
            self.mt_loadTk(self.mt_s,vars.sV)
            self.mt_loadTk(self.mt_v,vars.vV)
            self.mt_loadTk(self.mt_w,vars.wV)
            self.mt_loadTk(self.mt_x,vars.xV)
            self.mt_loadTk(self.mt_y,vars.yV)
            self.mt_loadTk(self.mt_z,vars.zV)
            self.mt_loadTk(self.mt_nparsim,vars.nparsimV)
            self.mt_loadTk(self.mt_sauto,vars.sautoV)
            self.mt_loadTk(self.mt_vauto,vars.vautoV)
            self.mt_loadTk(self.mt_wauto,vars.wautoV)
            self.mt_loadTk(self.mt_xauto,vars.xautoV)
            self.mt_loadTk(self.mt_yauto,vars.yautoV)
            self.mt_loadTk(self.mt_zauto,vars.zautoV)
            self.mt_loadTk(self.mt_nprocs,vars.nprocsV)
            self.mt_loadTk(self.mt_npsauto,vars.npsautoV)
            self.mt_loadTk(self.mt_hypz,vars.hypzV)
            self.mt_loadTk(self.mt_hypzord,vars.hypzordV)
            self.mt_loadTk(self.mt_hypx,vars.hypxV)
            self.mt_loadTk(self.mt_hypxord,vars.hypxordV)
            self.mt_loadTk(self.mt_hypy,vars.hypyV)
            self.mt_loadTk(self.mt_hypyord,vars.hypyordV)
            self.mt_loadTk(self.mt_hypv,vars.hypvV)
            self.mt_loadTk(self.mt_arakawazv,vars.arakawazvV)
            self.mt_loadTk(self.mt_rdchpt,vars.rdchptV)
            self.mt_loadTk(self.mt_wrtchpt,vars.wrtchptV)
            self.mt_loadTk(self.mt_chpth5,vars.chpth5V)
            self.mt_loadTk(self.mt_wrth5,vars.wrth5V)
            self.mt_loadTk(self.mt_wrtstd,vars.wrtstdV)
            self.mt_loadTk(self.mt_diagdir,vars.diagdirV)
            self.mt_loadTk(self.mt_chptdir,vars.chptdirV)
            self.mt_loadTk(self.mt_istepnrg,vars.istepnrgV)
            self.mt_loadTk(self.mt_istepomega,vars.istepomegaV)
            self.mt_loadTk(self.mt_istepfld,vars.istepfldV)
            self.mt_loadTk(self.mt_istepmom,vars.istepmomV)
            self.mt_loadTk(self.mt_istepenergy,vars.istepenergyV)
            self.mt_loadTk(self.mt_istepenergy3d,vars.istepenergy3dV)
            self.mt_loadTk(self.mt_istepprof,vars.istepprofV)
            self.mt_loadTk(self.mt_istepvsp,vars.istepvspV)
            self.mt_loadTk(self.mt_istepneoclass,vars.istepncV)
            self.mt_loadTk(self.mt_istepschpt,vars.istepschptV)
            self.mt_loadTk(self.mt_iterdbfile,vars.iterdbfileV)
            self.mt_loadTk(self.mt_iterdbtime,vars.iterdbtimeV)
            self.mt_loadTk(self.mt_nx,vars.nxV)
            self.mt_loadTk(self.mt_ny,vars.nyV)
            self.mt_loadTk(self.mt_nz,vars.nzV)
            self.mt_loadTk(self.mt_nv,vars.nvV)
            self.mt_loadTk(self.mt_nw,vars.nwV)
            self.mt_loadTk(self.mt_kymin,vars.kyminV)
            self.mt_loadTk(self.mt_lx,vars.lxV)
            self.mt_loadTk(self.mt_lv,vars.lvV)
            self.mt_loadTk(self.mt_lw,vars.lwV)
            self.mt_loadTk(self.mt_x0,vars.x0V)
            self.mt_loadTk(self.mt_ky0ind,vars.ky0indV)
            self.mt_loadTk(self.mt_nexc,vars.nexcV)
            self.mt_loadTk(self.mt_kxcent,vars.kxcentV)
            self.mt_loadTk(self.mt_initcond,vars.initcondV)
            self.mt_loadTk(self.mt_auxx,vars.auxxV)
            self.mt_loadTk(self.mt_auxy,vars.auxyV)
            self.mt_loadTk(self.mt_auxz,vars.auxzV) 
            self.mt_loadTk(self.mt_ntime,vars.ntimeV)
            self.mt_loadTk(self.mt_timelim,vars.timelimV)
            self.mt_loadTk(self.mt_simtimelim,vars.simtimelimV)
            self.mt_loadTk(self.mt_dtmax,vars.dtmaxV)
            self.mt_loadTk(self.mt_calcdt,vars.calcdtV)
            self.mt_loadTk(self.mt_ompr,vars.omprV)
            self.mt_loadTk(self.mt_ofl,vars.oflV)
            self.mt_loadTk(self.mt_ufl,vars.uflV)
            self.mt_loadTk(self.mt_beta,vars.betaV)
            self.mt_loadTk(self.mt_courant,vars.courantV)
            self.mt_loadTk(self.mt_incf0,vars.incf0V)
            self.mt_loadTk(self.mt_delphi,vars.delphiV)
            self.mt_loadTk(self.mt_coll,vars.collV)
            self.mt_loadTk(self.mt_Zeff,vars.ZeffV)
            self.mt_loadTk(self.mt_collis,vars.collisV)
            self.mt_loadTk(self.mt_collonoff,vars.collonoffV)
            self.mt_loadTk(self.mt_collcm,vars.collcmV)
            self.mt_loadTk(self.mt_collffm,vars.collffmV)
            self.mt_loadTk(self.mt_spacediff,vars.spacediffV)
            self.mt_loadTk(self.mt_debye2,vars.debye2V)
            self.mt_loadTk(self.mt_ivev,vars.ivev)
            self.mt_loadTk(self.mt_whichev,vars.whichevV)
            self.mt_loadTk(self.mt_nev,vars.nevV)
            self.mt_loadTk(self.mt_maxit,vars.maxitV)
            self.mt_loadTk(self.mt_shift,vars.shiftV)
            self.mt_loadTk(self.mt_ExBrate,vars.ExBrateV)
            self.mt_loadTk(self.mt_pfsrate,vars.pfsrateV)
            self.mt_loadTk(self.mt_ExBtime,vars.ExBtimeV)
            self.mt_loadTk(self.mt_kxExB,vars.kxExBV)
            self.mt_loadTk(self.mt_Phikx,vars.PhikxV)
            self.mt_loadTk(self.mt_press,vars.pressV)
            self.mt_loadTk(self.mt_delzonal,vars.delzonalV)
#            self.mt_loadTk(self.mt_quasilintm,vars.quasilintmV)
            self.mt_loadTk(self.mt_geomtype,vars.geomtype)
            self.mt_loadTk(self.mt_minorr,vars.minorrV)
            self.mt_loadTk(self.mt_length,vars.lengthV)
            self.mt_loadTk(self.mt_shat,vars.shatV)
            self.mt_loadTk(self.mt_trpeps,vars.trpepsV)
            self.mt_loadTk(self.mt_q0,vars.q0V)
            self.mt_loadTk(self.mt_amhd,vars.amhdV)
            self.mt_loadTk(self.mt_kappa,vars.kappaV)
            self.mt_loadTk(self.mt_s_kappa,vars.s_kappaV)
            self.mt_loadTk(self.mt_delta,vars.deltaV)
            self.mt_loadTk(self.mt_s_delta,vars.s_deltaV)
            self.mt_loadTk(self.mt_zeta,vars.zetaV)
            self.mt_loadTk(self.mt_s_zeta,vars.s_zetaV)
            self.mt_loadTk(self.mt_drR,vars.drRV)
            self.mt_loadTk(self.mt_geomdir,vars.geomdirV)
            self.mt_loadTk(self.mt_geomfile,vars.geomfileV)
            self.mt_loadTk(self.mt_fluxlabel,vars.fluxlabelV)
            self.mt_loadTk(self.mt_fluxpos,vars.fluxposV)
            self.mt_loadTk(self.mt_rref,vars.rrefV)
            self.mt_loadTk(self.mt_rhostar,vars.rhostarV)
            self.mt_loadTk(self.mt_magprof,vars.magprofV)
            self.mt_loadTk(self.mt_radbctype,vars.radbctypeV)
            self.mt_loadTk(self.mt_lcoefkrook,vars.lcoefkrookV)
            self.mt_loadTk(self.mt_lbufsize,vars.lbufsizeV)
            self.mt_loadTk(self.mt_lpow,vars.lpowkrookV)
            self.mt_loadTk(self.mt_ucoefkrook,vars.ucoefkrookV)
            self.mt_loadTk(self.mt_ubufsize,vars.ubufsizeV)
            self.mt_loadTk(self.mt_upow,vars.upowkrookV)
            self.mt_loadTk(self.mt_arakawa,vars.arakawaV)
            self.mt_loadTk(self.mt_ckheat,vars.ckheatV)
            self.mt_loadTk(self.mt_resetlimit,vars.resetlimitV)
            self.mt_loadTk(self.mt_shifted,vars.shiftedV)
            self.mt_loadTk(self.mt_Tref,vars.TrefV)
            self.mt_loadTk(self.mt_nref,vars.nrefV)
            self.mt_loadTk(self.mt_Bref,vars.BrefV)
            self.mt_loadTk(self.mt_Lref,vars.LrefV)
            self.mt_loadTk(self.mt_mref,vars.mrefV)
            self.mt_loadTk(self.mt_proftoTref,vars.proftoTrefV)
            self.mt_loadTk(self.mt_proftonref,vars.proftonrefV)
            self.mt_loadTk(self.mt_proftobeta,vars.proftobetaV)
            self.mt_loadTk(self.mt_proftocoll,vars.proftocollV)
            self.mt_loadTk(self.mt_proftodebye,vars.proftodebyeV)
            self.mt_loadTk(self.mt_proftoExB,vars.proftoExBV)
            self.mt_loadTk(self.mt_proftopfs,vars.proftopfsV)
            launcher.switch_lin_nonlin(weak=1)
            launcher.switch_ev_iv(weak=1)
            launcher.switch_nonlocal(weak=1)
            launcher.switch_collop()
            launcher.sel_geom()
            launcher.ch_numspec()
            vars.ch_spec1(1)
            if int(vars.numspecV.get())>=2:
                if len(vars.specname)>=2:
                    curr2=2
                else:
                    curr2=1
                    vars.ch_spec2(curr2)
        else:
            self.append=0
            self.index=-1
            vars.nonlocal.set(0)
            vars.clear()
            vars.newspec()
