PRO global_vars

  COMMON global_vars, cfg, gui, series, par, spec, nrg, mom, mom_time, $
    cdiag0, file_par, vm_sav, scanlog, vsp, vsp_time, chpt, chpt_time, $
    energy, energy_time, eterms, eterms_time, add_norm, show_global, $
    srcmom, srcmom_time,nltterms_time,nltterms,nlt2d,nlt_time,$
    nlp_gdt,nlp_kxind,nlp_kyind

  ; additional renormalization (mass and temperature) currently restricted
  ; since they still need to be tested
  add_norm = TOTAL(GETENV('USER') EQ ['dtold','flm','mirjam','mjpuesch',$
     'tbg']) GT 0 ? 1 : 0
  ; additional features for global version are currently restricted
  ; since they still need to be tested
  show_global = TOTAL(GETENV('USER') EQ ['global','dtold','flm','fsj','lapillon',$
    'mirjam','mjpuesch','sbrunner','tbg','tisd','xal','patma','hkd','tbird']) GT 0 ? 1 : 0

  out_format_str = ['data: ASCII','data: H5','postscript','png'] ;,'mpg']

  cfg = {$ ; diag configuration file settings (for customization)
    ps_viewer      : 'gv',$ ; program to view postscript files
    vm_lowcase     : 0,$ ; set to 'yes' (=1) for lowercase VM output file names
    short_gui      : 0,$ ; set to 'yes' (=1) for reduced height of GUI window
    speedup_gui    : 0,$ ; use an alternate GUI process that may be faster
    always_refresh : 0,$ ; refresh data not only for dat/act but for all runs
    scanlog_repair : 1,$ ; attempt to repair/create scanlogs automatically
    info_str       : 2,$ ; add info string (0: none, 1: time, 2: full) to bottom of page
    charsize       : 1,$ ; set font size
    charthick      : 1,$ ; set font thickness
    linethick      : 1,$ ; set line thickness
    bg_color       : 'gray',$ ; widget background color, e.g. lightgray, purple, #D0A0FF
    fg_color       : 'black',$ ; widget foreground color
    kpar_filter    : -1} ; zero or positive value: filter k_parallel 

  gui = {$ ; these variables pertain to the graphical user interface
    window         : {$
      main         : 0L,$
      var_list     : 0L,$
      sortmeth	   : 0L,$
      cdiags       : 0L,$
      help         : 0L},$
    tab     	   : {$
      base  	   : 0L,$
      arr 	   : PTR_NEW()},$
    text           : {$
      data_path    : 0L,$
      out_path     : 0L,$
      runs         : 0L,$
      start_t      : 0L,$
      end_t        : 0L,$
      sortmeth	   : 0L,$
      sparsefac    : 0L,$
      message      : 0L},$
    button         : {$
      data_path    : 0L,$
      path_bbnav   : 0L,$
      path_backnav : 0L,$
      path_forwnav : 0L,$
      path_ffnav   : 0L,$
      tmp1_path    : 0L,$
      tmp2_path    : 0L,$
      out_path     : 0L,$
      runs         : 0L,$
      nrg_data     : 0L,$
      series_info  : 0L,$
      refvalues_info : 0L,$
      first_step   : 0L,$
      last_step    : 0L,$
      geometry     : 0L,$
      profiles	   : 0L,$
      start        : 0L,$
      save         : 0L,$
      save_as      : 0L,$
      clear        : 0L,$
      load         : 0L,$
      var_list     : 0L,$
      cdiags       : 0L,$
      recent_h5    : 0L,$
      h5_files     : 0L,$
      recent_ps    : 0L,$
      ps_files     : 0L,$
      exit         : 0L},$
    bgroup         : {$
      species      : 0L,$
      out_format   : 0L,$
      resolve_steps: 0L},$
    droplist       : {$
      Lref_base    : 0L,$
      Lref         : 0L,$
      mref_base    : 0L,$
      mref  	   : 0L,$
      Tref_base    : 0L,$
      Tref  	   : 0L},$
    table          : {$
      var_list     : 0L,$
;      diags        : 0L,$
      cdiags       : 0L},$
    hist           : {$
      pos          : 0,$
      nhist        : 0,$
      ; change max number of stored previous settings here
      store        : REPLICATE({$
        dpath      : '',$
        runs       : '',$
        stime      : 0.0D,$
        etime      : 0.0D},20)},$
    misc           : {$
      n_diags      : 0,$
      n_cdiags     : 0,$
      n_sel_diags  : 0,$
      tmp1_path    : '',$
      tmp2_path    : '',$
      recent_h5    : '',$
      recent_ps    : '',$
      wcolors_set  : 0},$
    out            : {$
      data_path    : '',$
      out_path     : '',$
      start_t      : 0.0D,$
      current_gui  : '',$
      end_t        : 0.0D,$
      spec_select  : PTR_NEW(),$
      n_spec_sel   : 0,$
      sortmeth     : 0,$
      sparsefac    : 1L,$
      Lref_sel	   : 0,$
      res_steps    : 0,$
      out_format_str : out_format_str,$
      out_format   : INTARR(N_ELEMENTS(out_format_str))}}

  par = {$ ; these variables contain values from the parameter file
    n_spec         : 0,$
    n_procs_y      : 1,$
    n_procs_z      : 1,$
    n_procs_v      : 1,$
    n_procs_w      : 1,$
    nonlinear      : 0,$
    x_local 	   : 1,$
    y_local 	   : 1,$
    lilo           : 0,$
    comp_type      : 'IV',$ ; switch between EigenValue and Initial Value
    ntimesteps	   : 0L,$
    write_h5       : 0,$
    chpt_h5        : 0,$
    istep_nrg      : PTR_NEW(),$
    istep_field    : PTR_NEW(),$
    istep_mom      : PTR_NEW(),$
    istep_energy   : PTR_NEW(),$
    istep_energy3d : PTR_NEW(),$
;nonlinear transfer diagnostics
    istep_nlt      : PTR_NEW(),$
    num_nlt_modes  : 0,$
    kx_nlt_ind     : PTR_NEW(),$
    ky_nlt_ind     : PTR_NEW(),$
    nlt_old        : 0,$
    num_nlt_pod_modes  : 0,$
    nlt_pod        : 0,$
    nlp_gdt        : 0,$
    nlp_kxind        : 0,$
    nlp_kyind        : 0,$
;nonlinear transfer diagnostics
    n_ev           : 0,$
    which_ev       : '',$
    nkx0           : 0,$
    nx0            : 0,$
    nky0           : 0,$
    n0_global      : 0,$
    ny0            : 0,$
    nz0            : 0,$
    nv0            : 0,$
    nw0            : 0,$
    kx0_ind        : 0,$
    ky0_ind        : 0,$
    dt_max         : 0.0,$
    kx_center	   : 0.0,$
    kx	    	   : PTR_NEW(),$
    ky	    	   : PTR_NEW(),$
    adapt_lx	   : 0,$
    arakawa_cons_bc: 1,$
    arakawa_zv     : 1,$
    lxarr   	   : PTR_NEW(),$
    z	    	   : PTR_NEW(),$
    lx             : 0.0,$
    kymin          : 0.0,$
    ly             : 0.0,$
    lv             : 0.0,$
    lw             : 0.0,$
    mu_grid_type   : 'gau_leg',$
    n_conn    	   : 0,$
    n_pol   	   : 1,$
    dx             : 0.0,$
    dy             : 0.0,$
    dz             : 0.0,$
    beta           : 0.0,$
    debye          : 0.0,$
    collision_op   : '',$
    coll           : 0.0,$
    shifted_metric : 0,$
    ExBrate        : 0.0,$
    Omega_tor      : 0.0,$
    ExB_stime      : 100000,$
    pfsrate        : 0.0,$
    ehat           : 0.0,$
    q0	    	   : 0.0,$
    shat           : 0.0,$
    curv           : 0.0,$
    major_R 	   : 0.0,$
    major_Z        : 0.0,$
    minor_r        : 0.0,$
    parscale	   : 0.0,$
    amhd           : 0.0,$
    trpeps         : 0.0,$
    kappa          : 1.0,$
    s_kappa        : 0.0,$
    delta          : 0.0,$
    s_delta        : 0.0,$
    zeta           : 0.0,$
    s_zeta         : 0.0,$
    drR            : 0.0,$
    drZ            : 0.0,$
    mag_prof	   : 0,$
    hyp_x          : 0.0,$
    hyp_y          : 0.0,$
    hyp_z          : 0.0,$
    hyp_v          : 0.0,$
    hyp_perp       : 0.0,$
    delzonal       : 0,$
    n_fields       : 0,$
    n_moms         : 0,$
    n_energies     : 6,$
    n_srcmoms      : 0,$
    nrgcols        : 8,$
    magn_geometry  : '',$
    geomfile	   : '',$
    diag_trap_levels : 0,$
; --- experimental reference values
    Tref    	   : 1.0D,$
    mref    	   : 1.0D,$
    nref    	   : 1.0D,$
    Lref    	   : 1.0D,$
    Qref    	   : 1.0D,$
    Bref    	   : 1.0D,$
    omegatorref    : 0.0D,$
; --- derived expt. reference values
    cref           : 1.0D,$
    rhoref         : 1.0D,$
; --- for x nonlocal version
    rhostar        : 0.0D,$
    x0             : 0.5,$
    rad_bc_type    : 0,$
    reset_limit    : 0.0,$
; ---
    prec_single    : 0,$
    gene_version   : 11,$
    in_data 	   : 'sxky',$
    norm08  	   : 0,$ ; new normalization introduced in 2008
    norm_flux_projection : 1,$
    svn_rev        : ''}

  spec = {$
    name           : '',$
    passive        : 0,$
    omn            : 0.0,$
    omt            : 0.0,$
    mass           : 0.0,$
    charge         : 0,$
    temp           : 0.0,$
    dens           : 0.0,$
    prof           : PTR_NEW()}

  series = {$
    n_runs         : 0,$
    run_labels     : '',$
    par_updated    : 0,$
    filename_fmt   : PTR_NEW(),$ ; file name format
    request_mom    : INTARR(200,4),$ ;(3+48,4),$
    doppler_corr_mom : 0.0,$
    request_nrg    : 0,$
    request_energy : INTARR(6,4),$
    request_srcmom : INTARR(3,2),$
    request_eterms : 0,$
; --- reference values which are used for normalization
    Tref    	   : 1.0D,$
    Tref_str	   : 'T!Dref!N',$
    mref    	   : 1.0D,$
    mref_str	   : 'm!Dref!N',$
    nref    	   : 1.0D,$
    Lref    	   : 1.0D,$
    Lref_str	   : 'L!Dref!N',$
    Qref    	   : 1.0D,$
    Bref    	   : 1.0D,$
; --- normalized lengths
    kx	    	   : PTR_NEW(),$
    ky	    	   : PTR_NEW(),$
    lxarr   	   : PTR_NEW(),$
    lx             : 0.0,$
    kymin          : 0.0,$
    ly             : 0.0,$
    dy             : 0.0,$
    dx             : 0.0,$
; ---
    geom    	   : PTR_NEW(),$
    swap_endian    : 0,$
    step_count     : 0L,$
    h5_nst_field   : 0L,$
    h5_nst_mom     : 0L,$
    h5_nst_vsp     : 0L,$
    h5_count_field : 0L,$
    h5_count_mom   : 0L,$
    h5_count_vsp   : 0L}

  scanlog = PTR_NEW() ; see data_struct.pro: read_scanlog

  file_par = {exist:0,path:''}

  nrg = PTR_NEW()
  vsp = PTR_NEW()
  vsp_time = 0.0D

  chpt = PTR_NEW()
  chpt_time = 0.0D

  mom = {kxky : PTR_NEW(),$
         kxsy : PTR_NEW(),$
         sxky : PTR_NEW(),$
         sxsy : PTR_NEW()}

  energy = {kxky : PTR_NEW(),$
            kxsy : PTR_NEW(),$
            sxky : PTR_NEW(),$
            sxsy : PTR_NEW()}

  srcmom = PTR_NEW()

  mom_time = 0.0D
  energy_time = 0.0D
  srcmom_time = 0.0D

  eterms = PTR_NEW()
  eterms_time = PTR_NEW()

  nltterms = PTR_NEW()
  nltterms_time = PTR_NEW()
  nlt2d=PTR_NEW()
  nlt_time= 0.0D

  cdiag0 = PTR_NEW()

  vm_sav = 0 ; by default, no virtual machine sav support

END
