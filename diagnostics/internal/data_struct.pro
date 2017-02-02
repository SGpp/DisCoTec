;##########################################################################
;# this library contains almost all data structures which had too few     #
;# related internal functions to justify their own library                #
;##########################################################################

PRO read_par, archive

  COMMON global_vars

;  file_par.exist = 0

  run_labels = STRSPLIT(series.run_labels,',',/EXTRACT)
  n_runs = N_ELEMENTS(run_labels)

  ; reset all optional or later introduced parameters
  par.arakawa_cons_bc=0
  par.arakawa_zv=0
  par.delzonal = 0
  par.n_spec = 0
  par.magn_geometry = ''
  par.n_fields = 0
  par.n_moms = 0
  par.n_energies = 6
  par.n_srcmoms = 0
  par.nrgcols = 8
  par.svn_rev = ''
  par.prec_single = 0
  par.kx0_ind = 0
  par.ky0_ind = 0
  par.n0_global = 0
  par.write_h5 = 0
  par.chpt_h5 = 0
  istep_field = INTARR(n_runs) - 1 ; distinguish old runs from runs with no field file
  istep_mom = INTARR(n_runs)
  istep_nrg = INTARR(n_runs)
  istep_energy = INTARR(n_runs)
  istep_energy3d = INTARR(n_runs)
  istep_nlt = INTARR(n_runs)
  num_nlt_modes = 0
  nlt_old = 0
  kx_nlt_ind = INTARR(60)
  ky_nlt_ind = INTARR(60)
  par.num_nlt_pod_modes=0
  par.nlt_pod = 0
  par.nlp_gdt = 0
  par.nlp_kxind = 0
  par.nlp_kyind = 0
  par.n_ev = 0
  which_ev = ''
  par.lx = 0
  par.kx_center = 0.0
  par.mu_grid_type = 'gau_leg'
  nexc = 0
  par.n_conn = 0
  par.n_pol = 1
  par.ehat = 0
  par.curv = 0
  par.amhd = 0.0
  par.hyp_x = 0
  par.hyp_y = 0
  par.hyp_z = 0
  par.hyp_v = 0
  par.hyp_perp = 0
  par.debye = 0
  par.coll = 0
  par.ExBrate = 0
  par.ExB_stime = 100000
  par.pfsrate = 0
  par.shifted_metric = 0
  endianness = -1
  old_spec_standard = 0
  gene10_1_standard = 0
  dsmooth_x = ''
  dsmooth_z = ''
  par.gene_version = 11
  par.in_data = 'kxky'
  par.x_local = 1
  par.y_local = 1
  par.lilo = 0
  par.comp_type = 'IV'
  par.norm08 = 0
  par.norm_flux_projection = -1
  par.q0 = 0
  par.major_R = 0
  par.major_Z = 0
  par.minor_r = 0
  par.kappa = 1.0
  par.s_kappa = 0.0
  par.delta = 0.0
  par.s_delta = 0.0
  par.zeta = 0.0
  par.s_zeta = 0.0
  par.drR = 0.0
  par.drZ = 0.0
  par.parscale = 0
  par.mag_prof = 0
  par.diag_trap_levels = 0
  ; some estimated experimental values if none given in parameter file
  par.Bref = 1.0 ; Tesla
  par.Lref = 1.65 ; m
  par.mref = 3.3445e-27 ; m_D
  par.Qref = 1.6022e-19 ; e
  par.nref = 3.5  ; 10e19/m^3
  par.Tref = 0.35 ; keV
  par.omegatorref = 0.0 ;1/s = toroidal angular velocity
  ; variable for nonlocal x version
  par.rhostar = 0.0
  par.x0 = -1.0
  par.rad_bc_type = 0
  par.reset_limit = 0.0

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

  IF KEYWORD_SET(archive) THEN BEGIN
    old_file_par = {exist:0,path:''} ; store the original data path
    old_file_par = file_par
    file_par.exist = arx_file_par.exist
    file_par.path = arx_file_par.path[archive-1] ; do 'read_par' only for the current run, as in the archive, the runs are not necessarily follow-up runs
    dummy_labels = STRSPLIT(arx.run_labels,',',/EXTRACT)
    run_labels = dummy_labels[archive-1]
  ENDIF ELSE run_labels = STRSPLIT(series.run_labels,',',/EXTRACT)
  n_runs = N_ELEMENTS(run_labels)

  FOR run = 0, n_runs - 1 DO BEGIN
    OPENR, par_lun, file_par.path[run], /GET_LUN, ERROR=err
    IF err NE 0 THEN BEGIN
      PRINT, file_par.path[run] + ' not found'
      message = DIALOG_MESSAGE('par file(s) not found',/ERROR,$
        DIALOG_PARENT=gui.window.main)
      RETURN
    ENDIF ELSE IF !QUIET NE 1 THEN PRINT, 'reading ', file_par.path[run] ELSE $
      PRINT, 'reading run ', run_labels[run], ': par'

    old_spec_standard = -1
    n_lines = FILE_LINES(file_par.path[run])
    lines = STRARR(n_lines)
    READF, par_lun, lines

    FREE_LUN, par_lun

    n_spec_pos = STRPOS(lines,'n_spec')
    n_spec_ind = (WHERE(n_spec_pos EQ 0))[0]
    IF (n_spec_ind GE 0) THEN BEGIN
      IF n_spec_pos[n_spec_ind] EQ 0 THEN BEGIN
        temp = STRSPLIT(lines[n_spec_ind],'=',/EXTRACT)
        IF (STRTRIM(temp[0],2) NE 'n_spec') OR (N_ELEMENTS(temp) LT 2) THEN $
          printerror, 'n_spec in parameter file is not key in key_val pair' $
          ELSE BEGIN

          par.n_spec = FIX(STRTRIM(temp[1],2))
          old_spec_standard = 0
        ENDELSE
      ENDIF
    ENDIF ELSE BEGIN
      nspec_pos = STRPOS(lines,'nspec')
      nspec_ind = (WHERE(nspec_pos EQ 0))[0]
      IF nspec_pos[nspec_ind] EQ 0 THEN BEGIN
        temp = STRSPLIT(lines[nspec_ind],'=',/EXTRACT)
        IF (STRTRIM(temp[0],2) NE 'nspec') OR (N_ELEMENTS(temp) LT 2) THEN $
          printerror, 'nspec in parameter file is not key in key_val pair' $
          ELSE BEGIN

          par.n_spec = FIX(STRTRIM(temp[1],2))
          old_spec_standard = 1
        ENDELSE
      ENDIF
    ENDELSE

    IF run EQ 0 THEN spec = REPLICATE(spec[0],par.n_spec > 1) ELSE BEGIN
      oldpar = par
      oldspec = spec
    ENDELSE

    l = 0
    in_spec_block = 0
    spec_nr = -1
    intvar = 0L
    fltvar = 0.0
    strvar = ''
    WHILE l LT n_lines DO BEGIN
      line = lines[l]

      IF gene10_1_standard THEN BEGIN
        temp1 = STRSPLIT(line,':',/EXTRACT)
	IF N_ELEMENTS(temp1) EQ 2 THEN BEGIN
	  in_spec_block = 1
	  spec_nr = temp1[0]
    	  line = temp1[1]
	  spec[spec_nr].name = rm0es(spec_nr)
	ENDIF
      ENDIF

      temp = STRSPLIT(line,'=',/EXTRACT)
      key = STRTRIM(temp[0],2)
      IF N_ELEMENTS(temp) GT 1 THEN value = STRTRIM(temp[1],2) $
        ELSE value = ''
      IF value EQ 'T' THEN value = '1'
      IF value EQ 'F' THEN value = '0'

      IF (key EQ 'itime') THEN itime = FIX(value)

;      IF run EQ 0 THEN BEGIN
        IF in_spec_block EQ 1 THEN BEGIN
          CASE key OF
            'name'      : spec[spec_nr].name = STRMID(value,1,STRLEN(value)-2)
            'passive'   : spec[spec_nr].passive = FIX(value)
            'omn'       : spec[spec_nr].omn = FLOAT(value)
            'omt'       : spec[spec_nr].omt = FLOAT(value)
            'mass'      : spec[spec_nr].mass = FLOAT(value)
            'charge'    : spec[spec_nr].charge = FIX(value)
            'temp'      : spec[spec_nr].temp = FLOAT(value)
            'dens'      : spec[spec_nr].dens = FLOAT(value)
	    'n0'    	: spec[spec_nr].dens = FLOAT(value)
            ELSE :
          ENDCASE
        ENDIF
        IF old_spec_standard EQ 1 THEN BEGIN
          CASE key OF
            'etg'       : old_etg = FIX(value)
            'omne'      : old_omne = FLOAT(value)
            'omni'      : old_omni = FLOAT(value)
            'omte'      : old_omte = FLOAT(value)
            'omt'       : old_omte = FLOAT(value)
            'omti'      : old_omti = FLOAT(value)
            'omi'       : old_omti = FLOAT(value)
            'memi'      : old_memi = FLOAT(value)
            'teti'      : old_teti = FLOAT(value)
            ELSE :
          ENDCASE
        ENDIF
        CASE key OF
          '&species'    : IF spec_nr LT (par.n_spec - 1) THEN BEGIN
                            in_spec_block = 1
                            spec_nr += 1
                          END
          '/'           : in_spec_block = 0
          ; geometry namelist has been introduced with new normalization
          ; which is why it is used as an indicator here
	  '&geometry'   : par.norm08 = 1
          ; --- integer values ---
          'ntimesteps'  : IF STRPOS(value,'*') GT -1 THEN $
                            par.ntimesteps = 100000L ELSE $
                            par.ntimesteps = LONG(value)
          'write_h5'    : par.write_h5 = FIX(value)
          'chpt_h5'     : par.chpt_h5 = FIX(value)
          'istep_field' : istep_field[run] = FIX(value)
          'istep_mom'   : istep_mom[run] = FIX(value)
          'istep_nrg'   : istep_nrg[run] = FIX(value)
          'istep_energy': istep_energy[run] = FIX(value)
          'istep_energy3d': istep_energy3d[run] = FIX(value)
          'istep_nlt'   : istep_nlt[run] = FIX(value)
          'num_nlt_modes'  : par.num_nlt_modes = FIX(value)
          'num_nlt_pod_modes'  : par.num_nlt_pod_modes = FIX(value)
          'nlt_pod'     : par.nlt_pod = FIX(value)
          'nlp_gdt'     : par.nlp_gdt = FIX(value)
          'nlp_kxind'   : par.nlp_kxind = FIX(value)
          'nlp_kyind'   : par.nlp_kyind = FIX(value)
          'nlt_old'     : par.nlt_old = FIX(value)
          'kx_nlt_ind'  : kx_nlt_ind[*] = FIX(STRSPLIT(value,' ',/EXTRACT))
          'ky_nlt_ind'  : ky_nlt_ind[*] = FIX(STRSPLIT(value,' ',/EXTRACT))
          'n_ev'        : par.n_ev = FIX(value)
          'nx0'         : par.nkx0 = FIX(value)
          'nky0'        : par.nky0 = FIX(value)
          'n0_global'   : par.n0_global = FIX(value)
          'nz0'         : par.nz0 = FIX(value)
          'nv0'         : par.nv0 = FIX(value)
          'nw0'         : par.nw0 = FIX(value)
          'n_procs_y'   : par.n_procs_y = FIX(value)
          'n_procs_z'   : par.n_procs_z = FIX(value)
          'n_procs_v'   : par.n_procs_v = FIX(value)
          'n_procs_w'   : par.n_procs_w = FIX(value)
          'ky0_ind'     : par.ky0_ind = FIX(value)
          'kx0_ind'     : par.kx0_ind = FIX(value)
          'n_fields'    : par.n_fields = FIX(value)
          'n_moms'      : par.n_moms = FIX(value)
          'n_energies'  : par.n_energies = FIX(value)
          'n_srcmoms'   : par.n_srcmoms = FIX(value)
          'nrgcols'     : par.nrgcols = FIX(value)
	  'n_pol'   	: par.n_pol = FIX(value)
          'nonlinear'   : par.nonlinear = FIX(value)
          'arakawa_cons_bc': par.arakawa_cons_bc = FIX(value)
          'arakawa_zv': par.arakawa_zv = FIX(value)
          'delzonal'    : par.delzonal = FIX(value)
	  'x_local'     : par.x_local = FIX(value)
	  'y_local' 	: par.y_local = FIX(value)
          'lilo'        : par.lilo = FIX(value)
	  'mag_prof'	: par.mag_prof = FIX(value)
          'rad_bc_type' : par.rad_bc_type = FIX(value)
          'diag_trap_levels' : par.diag_trap_levels = FIX(value)
          'shifted_metric' : par.shifted_metric = FIX(value)
          'norm_flux_projection' : par.norm_flux_projection = FIX(value)
          ; --- float values ---
          'dt_max'      : par.dt_max = FLOAT(value)
          'kx_center'   : par.kx_center = FLOAT(value)
          'lx'          : par.lx = FLOAT(value)
          'kymin'       : par.kymin = FLOAT(value)
          'lv'          : par.lv = FLOAT(value)
          'lz'          : par.lv = FLOAT(value)
          'lw'          : par.lw = FLOAT(value)
          'beta'        : par.beta = FLOAT(value)
          'debye'       : par.debye = FLOAT(value)
          'debye2'      : par.debye = SQRT(FLOAT(value))
          'coll'        : par.coll = FLOAT(value)
          'ExBrate'     : par.ExBrate = FLOAT(value)
          'ExB_stime'   : par.ExB_stime = FLOAT(value)
          'pfsrate'     : par.pfsrate = FLOAT(value)
          'ehat'        : par.ehat = FLOAT(value)
          'q0'          : par.q0 = FLOAT(value)
          'shat'        : par.shat = FLOAT(value)
          'curv'        : par.curv = FLOAT(value)
          'major_R'     : par.major_R = FLOAT(value)
          'major_Z'     : par.major_Z = FLOAT(value)
          'minor_r'     : par.minor_r = FLOAT(value)
          'parscale'    : par.parscale = FLOAT(value)
          'amhd'        : par.amhd = FLOAT(value)
          'trpeps'      : par.trpeps = FLOAT(value)
          'kappa'       : par.kappa = FLOAT(value)
          's_kappa'     : par.s_kappa = FLOAT(value)
          'delta'       : par.delta = FLOAT(value)
          's_delta'     : par.s_delta = FLOAT(value)
          'zeta'        : par.zeta = FLOAT(value)
          's_zeta'      : par.s_zeta = FLOAT(value)
          'drR'         : par.drR = FLOAT(value)
          'drZ'         : par.drZ = FLOAT(value)
          'hyp_x'       : par.hyp_x = FLOAT(value)
          'hyp_y'       : par.hyp_y = FLOAT(value)
          'hyp_z'       : par.hyp_z = FLOAT(value)
          'hyp_v'       : par.hyp_v = FLOAT(value)
          'hyp_perp'    : par.hyp_perp = FLOAT(value)
          'rhostar'     : par.rhostar = FLOAT(value)
          'reset_limit' : par.reset_limit = FLOAT(value)
          'x0'          : par.x0 = FLOAT(value)
	  'Tref'    	: par.Tref = FLOAT(value) ; in keV
	  'Lref'        : par.Lref = FLOAT(value) ; in m
	  'nref'    	: par.nref = FLOAT(value) ; in 1e19/m^3
	  'Bref'    	: par.Bref = FLOAT(value) ; in T
	  'Qref'    	: par.Qref = FLOAT(value) ; in C
          'mref'        : par.mref = FLOAT(value) * 1.6726231e-27 ; in kg
          'omegatorref' : par.omegatorref = FLOAT(value) ; in 1/s
          'vtorref'     : par.omegatorref = FLOAT(value) ; in 1/s ;backward-comp.
        ; --- string values ---
          'mu_grid_type': par.mu_grid_type = $
                            STRSPLIT(value,"'",/EXTRACT)
          'comp_type'   : par.comp_type = value
          'collision_op': par.collision_op = $
                            STRSPLIT(value,"'",/EXTRACT)
	  'magn_geometry' : par.magn_geometry = $
                            STRSPLIT(value,"'",/EXTRACT)
	  'geomfile'    : par.geomfile = $
                            STRSPLIT(value,"'",/EXTRACT)
          'which_ev'    : which_ev = value
          'TYPE'        : IF value EQ 'GLOBAL' THEN par.x_local = 0
          'SVN_REV'     : par.svn_rev = value
          'PRECISION'   : IF value EQ 'SINGLE' THEN $
                            par.prec_single = 1
          'ENDIANNESS'  : endianness = value EQ 'LITTLE' ? 1 : 0
        ; --- compatibility with old par files ---
          'istep_kin'   : BEGIN
                            istep_field[run] = FIX(value)
                            istep_mom[run] = FIX(value)
                          END
          'nkx0'        : par.nkx0 = FIX(value)
          'ns0'         : par.nz0 = FIX(value)
          'nz'          : par.nv0 = 2 * FIX(value)
          'nw'          : par.nw0 = FIX(value)
          'totalnw'     : par.nw0 = FIX(value)
          'n_pesy'      : par.n_procs_y = FIX(value)
          'n_pess'      : par.n_procs_z = FIX(value)
          'n_pesz'      : par.n_procs_v = FIX(value)
          'n_pesw'      : par.n_procs_w = FIX(value)
          'nfields'     : nfields = FIX(value)
    	  'nexc'    	: nexc = FIX(value)
          'stellarator' : IF FIX(value) THEN par.magn_geometry = 'stell.dat'
          'NONLINEAR'   : par.nonlinear = FIX(value)
          'arakawa_cons_bc': par.arakawa_cons_bc = FIX(value)
          'arakawa_zv': par.arakawa_zv = FIX(value)
          'PREC_SINGLE' : par.prec_single = FIX(value)
          'DELZONAL'    : par.delzonal = FIX(value)
          'DSMOOTH_X'   : dsmooth_x = value
          'DSMOOTH_Z'   : dsmooth_z = value
          'hyp_vp'      : par.hyp_v = FLOAT(value)
	  'lasttime'	: gene10_1_standard = 1
          ELSE :
        ENDCASE
;      ENDIF
      l += 1
    ENDWHILE

    IF old_spec_standard EQ 1 THEN BEGIN
;     beta used to be 4*pi*p/B^2 in older Gene versions,
;     now it is 8*pi*p/B^2
      par.beta = 2.0 * par.beta
      IF NOT gene10_1_standard THEN BEGIN
        IF old_etg EQ 0 THEN BEGIN
          name = ['i','e']
          omn = [old_omni,old_omne]
          omt = [old_omti,old_omte]
          mass = [1.0,old_memi]
          temp = [1.0/old_teti,1.0]
          charge = [1,-1]
        ENDIF ELSE BEGIN
          name = ['e','i']
          omn = [old_omne,old_omni]
          omt = [old_omte,old_omti]
          mass = [1.0,1.0/old_memi]
          temp = [1.0,1.0/old_teti]
          charge = [-1,1]
        ENDELSE
        FOR sp = 0, par.n_spec - 1 DO spec[sp] = $
          {name:name[sp],passive:0,omn:omn[sp],omt:omt[sp],mass:mass[sp],$
          charge:charge[sp],temp:temp[sp],dens:1.0,prof:PTR_NEW()}
        IF dsmooth_x EQ '' THEN par.hyp_x = 0
        IF dsmooth_z EQ '' THEN par.hyp_z = 0
      ENDIF ELSE par.gene_version = 10
      IF (par.n_fields EQ 0) AND (par.n_moms EQ 0) THEN BEGIN
        par.n_fields = 2
        par.n_moms = nfields
      ENDIF
    ENDIF

    IF run GT 0 THEN $ ; compare par and spec with previous run
      IF (oldpar.nkx0 NE par.nkx0) OR (oldpar.nky0 NE par.nky0) OR (oldpar.nz0 NE par.nz0) $
      OR (oldpar.lx NE par.lx) OR (oldpar.kymin NE par.kymin) $
      OR (oldpar.n_fields NE par.n_fields) OR (oldpar.n_moms NE par.n_moms) $
      OR (TOTAL((oldspec[0:par.n_spec-1].omt LT 0.999 * spec[0:par.n_spec-1].omt) OR $
      (oldspec[0:par.n_spec-1].omt GT 1.001 * spec[0:par.n_spec-1].omt)) NE 0) $
      OR (TOTAL((oldspec[0:par.n_spec-1].omn LT 0.999 * spec[0:par.n_spec-1].omn) OR $
      (oldspec[0:par.n_spec-1].omn GT 1.001 * spec[0:par.n_spec-1].omn)) NE 0) OR $
      (oldpar.comp_type NE par.comp_type) THEN BEGIN

      printerror, 'serious run parameter mismatch', /popup
;      RETURN
    ENDIF

    IF istep_field[run] EQ -1 THEN istep_field[run] = istep_mom[run]
    IF (istep_energy[run] GT 0) AND (istep_energy3d[run] EQ 0) THEN $
       istep_energy3d[run] = istep_energy[run]

  ENDFOR

  IF (par.nkx0 EQ 0) OR (par.nky0 EQ 0) OR (par.nz0 EQ 0) OR $
    (par.n_fields EQ 0) OR (par.n_spec EQ 0) THEN BEGIN
    printerror, $
      'no nkx0, nky0, nz0, n_spec, n_fields, or n_moms in parameters', $
      /popup
;    RETURN
  ENDIF

  ; for eigenvalue runs, the istep values are irrelevant
  IF par.comp_type EQ 'EV' THEN BEGIN
    par.istep_field = 1
    par.istep_mom = 1
    par.istep_nrg = 1
    par.istep_energy = 1
    par.istep_energy3d = 1
    par.istep_nlt = 1
    par.istep_kx_nlt_ind = 1
    par.istep_ky_nlt_ind = 1
  ENDIF

  IF PTR_VALID(par.istep_field) THEN PTR_FREE, par.istep_field
  IF PTR_VALID(par.istep_mom) THEN PTR_FREE, par.istep_mom
  IF PTR_VALID(par.istep_nrg) THEN PTR_FREE, par.istep_nrg
  IF PTR_VALID(par.istep_energy) THEN PTR_FREE, par.istep_energy
  IF PTR_VALID(par.istep_energy3d) THEN PTR_FREE, par.istep_energy3d
  IF PTR_VALID(par.istep_nlt) THEN PTR_FREE, par.istep_nlt
  IF PTR_VALID(par.kx_nlt_ind) THEN PTR_FREE, par.kx_nlt_ind
  IF PTR_VALID(par.ky_nlt_ind) THEN PTR_FREE, par.ky_nlt_ind

  par.istep_field = PTR_NEW(istep_field)
  par.istep_mom = PTR_NEW(istep_mom)
  par.istep_nrg = PTR_NEW(istep_nrg)
  par.istep_energy = PTR_NEW(istep_energy)
  par.istep_energy3d = PTR_NEW(istep_energy3d)
  par.istep_nlt = PTR_NEW(istep_nlt)
  par.kx_nlt_ind = PTR_NEW(kx_nlt_ind)
  par.ky_nlt_ind = PTR_NEW(ky_nlt_ind)

  IF (which_ev EQ 'none') OR (which_ev EQ '') THEN par.n_ev = 0
  IF par.n_srcmoms EQ 0 THEN par.n_srcmoms = 3

  IF par.x_local THEN par.shifted_metric = 0

  ; --- calculate derived par elements ----

; Warning for all developers!
; par.nky0 is modified if ky0_ind > 0 since the dummy modes ky=[0,kymin,..,(ky0_ind-1)*kymin]
; are added to allow for ky->y FFT
; Therefore, par.nky0 is not necessarily identical to the parameters value
  IF par.ky0_ind GT 0 THEN BEGIN
    par.nky0 += par.ky0_ind
;    PRINT, 'IMPORTANT INFORMATION: the dummy modes ky=[0,kymin,..,(ky0_ind-1)*kymin] are '+$
;      'added to allow for ky->y FFT if ky0_ind > 0. Please consider this modification '+$
;      'while using ky indices and/or self-written diagnostics'
  ENDIF

  IF (par.y_local) THEN BEGIN
     par.ny0 = 2 * par.nky0 
  ENDIF ELSE BEGIN
     par.ny0 = par.nky0
     par.nky0 = par.ny0/2

     par.nkx0 = 2 * par.nkx0   ;the diagnostics tools always
                                ;considers negative and positive 
                                ;kx values, i.e. sykx data is
                                ;ny0*nkx0/2
  ENDELSE
  par.nx0 = par.nkx0

  IF par.kymin GT 1e-6 THEN par.ly = 2.0 * !DPI / par.kymin
  IF par.x_local THEN par.dx = DOUBLE(par.lx) / par.nx0 $
    ELSE par.dx = DOUBLE(par.lx) / (par.nx0 - 1)
  par.dy = DOUBLE(par.ly) / par.ny0
  par.dz = 2.0D * !DPI * par.n_pol / par.nz0

  IF (par.x0 LT 0) THEN BEGIN
     IF (par.x_local) THEN BEGIN
        IF (par.magn_geometry EQ 'circular' OR par.magn_geometry EQ 's_alpha')$
           AND par.trpeps GE 0 AND par.major_R GE 0 THEN $
              par.x0 = par.trpeps*par.major_R/par.minor_r
     ENDIF ELSE BEGIN
        IF (par.magn_geometry EQ 'circular') THEN $
           par.x0 = 0.5 ;set default to 0.5 for old runs
     ENDELSE
  ENDIF ELSE BEGIN
     IF (par.trpeps LE 0) THEN par.trpeps = par.x0*par.minor_r/par.major_R
  ENDELSE

  IF PTR_VALID(par.ky) THEN PTR_FREE, par.ky
  IF PTR_VALID(par.kx) THEN PTR_FREE, par.kx
  IF PTR_VALID(par.lxarr) THEN PTR_FREE, par.lxarr
  IF PTR_VALID(par.z) THEN PTR_FREE, par.z

  par.ky = PTR_NEW(par.kymin*INDGEN(par.nky0))
  par.lxarr = PTR_NEW(DBLARR(par.nky0,/NOZERO))
  par.kx = PTR_NEW(DBLARR(par.nkx0,par.nky0,/NOZERO))
  par.z = PTR_NEW(!DPI*(-par.n_pol+2.0D*par.n_pol/par.nz0*INDGEN(par.nz0)))
  IF ((par.nz0 MOD 2) NE 0) THEN (*par.z) += !DPI*par.n_pol/par.nz0
  

  IF (par.lx EQ 0) AND (par.nky0-par.ky0_ind EQ 1) THEN $
    IF nexc NE 0 THEN par.lx = ABS(nexc) / (par.kymin * ABS(par.shat)) $
    ELSE par.lx = 1.0 / (par.kymin * ABS(par.shat))

  par.adapt_lx = (par.lx EQ 0) AND (nexc EQ 0)

  IF par.adapt_lx THEN par.n_conn = 1

  IF ABS(par.shat) GT 1e-3 THEN BEGIN
    lxmin = 1.0 / (ABS(par.shat) * par.kymin)
    IF par.adapt_lx THEN BEGIN
      (*par.lxarr)[0] = par.lx
      IF par.nky0 GT 1 THEN $
        (*par.lxarr)[1:*] = 1.0 / (ABS(par.shat) * (*par.ky)[1:*])
    ENDIF ELSE (*par.lxarr)[*] = par.lx
  ENDIF ELSE (*par.lxarr)[*] = par.lx

  IF par.nkx0 GT 1 THEN BEGIN
    FOR ikx = 0, par.nkx0 / 2 DO BEGIN
      (*par.kx)[ikx,*] = par.kx_center
      lxarr_fin = WHERE(*par.lxarr NE 0.0)
      IF lxarr_fin[0] NE -1 THEN (*par.kx)[ikx,lxarr_fin] = $
        2.0 * !DPI / (*par.lxarr)[lxarr_fin] * ikx + par.kx_center
    ENDFOR
    FOR ikx = par.nkx0 / 2 + 1, par.nkx0 - 1 DO $
      (*par.kx)[ikx,*] = - (*par.kx)[par.nkx0-ikx,*] + 2.0 * par.kx_center
  ENDIF ELSE IF par.kx0_ind GT 0 THEN $
   (*par.kx)[*,*] = 2.0 * !DPI / (*par.lxarr) * par.kx0_ind $
   ELSE (*par.kx)[*,*] = par.kx_center

  IF (par.major_R EQ 0) AND (par.curv NE 0) THEN par.major_R = 2.0 / par.curv
  IF (par.curv EQ 0) AND (par.major_R NE 0) THEN par.curv = 2.0 / par.major_R
  IF par.q0 EQ 0 THEN par.q0 = par.curv / 2.0 * SQRT(par.ehat)
  IF (par.ehat EQ 0) AND (par.q0 NE 0) AND (par.curv NE 0) THEN $
    par.ehat = (2.0 * par.q0 / par.curv)^2

  series.swap_endian = 0
  IF (endianness NE -1) AND ((BYTE(1,0,1))[0] NE endianness) THEN $
    series.swap_endian = 1

  IF par.collision_op EQ 'none' THEN par.coll = 0.0

  IF par.magn_geometry EQ '' THEN BEGIN
    par.magn_geometry = 's_alpha'
    PRINT, 'using internal s_alpha geometry coefficients'
;    RETURN
  ENDIF

  IF (par.norm_flux_projection EQ -1) THEN BEGIN
     par.norm_flux_projection = 0
     IF (par.svn_rev EQ 'exported') THEN BEGIN
        print, "couldn't determine flux projection; assuming norm_flux_projection=F"
     ENDIF ELSE IF (FIX(par.svn_rev) LE 0) THEN BEGIN
        print, "couldn't determine flux projection; assuming norm_flux_projection=F"
     ENDIF ELSE par.norm_flux_projection = ((FIX(par.svn_rev) GE 3190) AND $
         (FIX(par.svn_rev) LT 3846))
  ENDIF

  par.in_data = 'kxky'
  IF ((par.gene_version EQ 10) OR (par.x_local EQ 0)) THEN par.in_data = 'sxky'
  IF (par.y_local EQ 0) THEN par.in_data = 'sykx'

  IF (NOT par.x_local) AND (par.minor_r LE 0) AND (par.trpeps GT 0) THEN $
    par.minor_r = 2.0 * par.trpeps

  ; change units of experiment data
  par.Tref = par.Tref * 1e3 ; keV -> eV
  par.nref = par.nref * 1e19

  ; derived experimental reference values
  par.cref = SQRT(par.Tref*par.Qref/par.mref) ; in m/s
  par.rhoref = SQRT(par.mref*par.Tref/par.Qref) / par.Bref ; in m

  file_par.exist = 1
  series.par_updated = 1

  IF KEYWORD_SET(archive) THEN $
    file_par = old_file_par ; return to originally requested data paths

END

;##########################################################################

PRO set_series_lengths

  COMMON global_vars

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref ; T in eV

  IF PTR_VALID(series.ky) THEN PTR_FREE, series.ky
  IF PTR_VALID(series.kx) THEN PTR_FREE, series.kx
  IF PTR_VALID(series.lxarr) THEN PTR_FREE, series.lxarr

  series.kymin = par.kymin / rho_ref
  series.ly = par.ly * rho_ref
  series.lx = par.lx * rho_ref

  IF par.x_local THEN series.dx = DOUBLE(series.lx) / par.nx0 $
    ELSE series.dx = DOUBLE(series.lx) / (par.nx0 - 1)
  series.dy = DOUBLE(series.ly) / par.ny0
  series.ky = PTR_NEW(series.kymin*INDGEN(par.nky0))
  series.lxarr = PTR_NEW(DBLARR(par.nky0,/NOZERO))
  series.kx = PTR_NEW(DBLARR(par.nkx0, par.nky0,/NOZERO))

  IF (series.lx EQ 0) AND (par.nky0-par.ky0_ind EQ 1) THEN $
    IF nexc NE 0 THEN series.lx = ABS(nexc) / (series.kymin * ABS(par.shat)) $
    ELSE series.lx = 1.0 / (series.kymin * ABS(par.shat))

  IF ABS(par.shat) GT 1e-3 THEN BEGIN
    lxmin = 1.0 / (ABS(par.shat) * series.kymin)
    IF par.adapt_lx THEN BEGIN
      (*series.lxarr)[0] = par.lx
      IF par.nky0 GT 1 THEN $
        (*series.lxarr)[1:*] = 1.0 / (ABS(par.shat) * (*series.ky)[1:*])
    ENDIF ELSE (*series.lxarr)[*] = series.lx
  ENDIF ELSE (*series.lxarr)[*] = series.lx

  IF (par.nkx0 GT 1) THEN BEGIN
    FOR ikx = 0, par.nkx0 / 2 DO BEGIN
      (*series.kx)[ikx,*] = par.kx_center / rho_ref
      lxarr_fin = WHERE(*series.lxarr NE 0.0)
      IF lxarr_fin[0] NE -1 THEN (*series.kx)[ikx,lxarr_fin] = $
        2.0 * !DPI / (*series.lxarr)[lxarr_fin] * ikx + par.kx_center / rho_ref
    ENDFOR
    FOR ikx = par.nkx0 / 2 + 1, par.nkx0 - 1 DO $
      (*series.kx)[ikx,*] = - (*series.kx)[par.nkx0-ikx,*] + 2.0 * par.kx_center / rho_ref
  ENDIF ELSE IF par.kx0_ind GT 0 THEN $
   (*series.kx)[*,*] = 2.0 * !DPI / (*series.lxarr) * par.kx0_ind $
   ELSE (*series.kx)[*,*] = par.kx_center / rho_ref

END

;##########################################################################

FUNCTION check_gridpoints, geo_lun, geo_file

  COMMON global_vars

  IF par.write_h5 THEN BEGIN
    gridp_dset = H5D_OPEN(geo_lun,'parameters/gridpoints')
    gridpoints = H5D_READ(gridp_dset)
    H5D_CLOSE, gridp_dset
  ENDIF ELSE BEGIN
    line = ''
    READF, geo_lun, line
    ; try to read namelist containing gridpoints
    IF STRMID(line,0,1) EQ '&' THEN BEGIN
       WHILE NOT EOF(geo_lun) AND STRPOS(line,'gridpoints') LT 0 DO $
         READF, geo_lun, line
       IF EOF(geo_lun) THEN BEGIN
         printerror,'corrupted geometry file: no gridpoint information found in ' + $
           geo_file + ', using s-alpha instead', /popup
         RETURN, 0
       ENDIF ELSE BEGIN
         gridpoints = FIX((STRSPLIT(line,'=',/EXTRACT))[1])
         ; go to end of namelist
         WHILE NOT EOF(geo_lun) AND (STRMID(line,0,1) NE '/') AND $
           (STRMID(line,0,1) NE '\') DO READF, geo_lun, line
       ENDELSE
    ENDIF ELSE BEGIN ; no namelist found
      POINT_LUN, geo_lun, 0 ; rewind file
      gridpoints = (FSTAT(geo_lun)).SIZE / 161 ; 161 = 8*20 chars + end of line
    ENDELSE
  ENDELSE

  IF FIX(gridpoints) NE par.nz0 THEN BEGIN
    printerror, 'gridpoints (' + rm0es(gridpoints) + ') in ' + geo_file + $
      ' do not match nz0 (' + rm0es(par.nz0) + '), using s-alpha'
    RETURN, 0
  ENDIF

  RETURN, 1

END

; ############################################################################

PRO set_s_alpha, geom_struct

  COMMON global_vars

  theta = - !DPI * par.n_pol + INDGEN(par.nz0) * par.dz
  arg = theta * par.shat - par.amhd * SIN(theta)

  geom_struct.gxx = 1.0
  geom_struct.gxz = 0.0
  IF par.trpeps NE 0 THEN BEGIN
    geom_struct.gyz = 1.0 / (par.trpeps * par.major_R)
    geom_struct.gzz = 1.0 / (par.trpeps * par.major_R)^2
  ENDIF ELSE BEGIN
    geom_struct.gyz = 0
    geom_struct.gzz = 0
  ENDELSE
  geom_struct.dBdy = 0.0

  IF (SIZE(geom_struct.gxy,/N_DIMENSIONS) EQ 1) THEN BEGIN
    geom_struct.gxy = arg
    geom_struct.gyy = arg^2 + 1.0
    geom_struct.Bfield = 1.0 / (1.0 + par.trpeps * COS(theta))
    geom_struct.dBdx = - COS(theta) * geom_struct.Bfield^2 / par.major_R
    geom_struct.dBdz = par.trpeps * $
      SIN(theta) / ((par.trpeps * COS(theta))^2 + 1.0)
    ; dBdz is former dBfield
  ENDIF ELSE BEGIN
    ; does not make real sense and is only implemented for testing
    FOR i = 0, par.nx0 - 1 DO BEGIN
      geom_struct.gxz[i,*] = arg
      geom_struct.gyy[i,*] = arg^2+1.0
      geom_struct.Bfield[i,*] = 1.0 / (1.0 + par.trpeps * COS(theta))
      geom_struct.dBdx[i,*] = - COS(theta) * geom_struct.Bfield^2 / par.major_R
      geom_struct.dBdz[i,*] = par.trpeps * $
        SIN(theta) / ((par.trpeps * COS(theta))^2 + 1.0)
    ENDFOR
  ENDELSE

  geom_struct.jacobian = 1.0 / geom_struct.Bfield
  IF par.norm08 THEN $
    IF par.magn_geometry EQ 'slab' THEN $
      geom_struct.jacobian *= par.parscale $
    ELSE geom_struct.jacobian *= par.q0

END

; ############################################################################

PRO read_stell_dat, geo_lun, geom
; reads stell.dat files (the oldest TRACER-Gene interface)

  COMMON global_vars

  nz0 = 0
  READF, geo_lun, nz0
  IF nz0 NE par.nz0 THEN BEGIN
     printerror, 'stell nz0 does not match par nz0, using s-alpha'
     set_s_alpha, geom_struct
     RETURN
  ENDIF

  ; values arrange in lines
  temp_arr = DBLARR(nz0,/NOZERO)
  READF, geo_lun, temp_arr & geom.gxx = temp_arr
  READF, geo_lun, temp_arr & geom.gxy = temp_arr
  READF, geo_lun, temp_arr & geom.gyy = temp_arr
  READF, geo_lun, temp_arr & geom.Bfield = temp_arr
  READF, geo_lun, temp_arr ;& geom.wcvx = par.curv * temp_arr
  READF, geo_lun, temp_arr ;& geom.wcvy = par.curv * temp_arr
  READF, geo_lun, temp_arr & geom.jacobian = temp_arr

  ;calculate dBds
;    FOR k=1,par.nz0-2 DO $
;      (*series.geom).dBdz[k]=((*series.geom).Bfield[k+1]-$
;        (*series.geom).Bfield[k-1])/(2.0*par.dz)
;    (*series.geom).dBdz[0]=0.0 & (*series.geom).dBdz[par.nz0-1]=0.0

END

; ############################################################################

PRO read_tracer_ff, geo_lun, geom

  COMMON global_vars

  IF par.norm08 THEN BEGIN
    tmp_geo_data = DBLARR(16,par.nz0,/NOZERO)

    IF par.write_h5 THEN BEGIN
      h5_fields = ['g_xx','g_xy','g_xz','g_yy','g_yz','g_zz',$
        'Bfield','dBdx','dBdy','dBdz','Jacobian',$
        'g_R','g_phi','g_Z','c_1','c_2']
      h5_fields[0:5] = 'metric/' + h5_fields[0:5]
      h5_fields[6:10] = 'Bfield_terms/' + h5_fields[6:10]
      h5_fields[11:15] = 'shape/' + h5_fields[11:15]

      bft_members = H5G_GET_NMEMBERS(geo_lun,'Bfield_terms/')

      FOR j = 0, 15 DO BEGIN
        IF ((bft_members NE 4) OR (j NE 8)) THEN BEGIN
          geo_dset = H5D_OPEN(geo_lun,h5_fields[j])
          geo_field = H5D_READ(geo_dset)
          tmp_geo_data[j,*] = geo_field
          H5D_CLOSE, geo_dset
        ENDIF
      ENDFOR
    ENDIF ELSE BEGIN
       ;rewind file and try to read namelist
       POINT_LUN, geo_lun, 0
       line = ''
       WHILE NOT EOF(geo_lun) AND (STRMID(line,0,1) NE '/') AND $
           (STRMID(line,0,1) NE '\') DO BEGIN
          READF, geo_lun, line
          temp = STRSPLIT(line,'=',/EXTRACT)
          IF (STRTRIM(temp[0]) EQ 'Cy') THEN geom.C_y = float(temp[1])
          IF (STRTRIM(temp[0]) EQ 'Cxy') THEN geom.C_xy = float(temp[1])
       ENDWHILE
       READF, geo_lun, tmp_geo_data
    ENDELSE

    geom.gxx = tmp_geo_data[0,*]
    geom.gxy = tmp_geo_data[1,*]
    geom.gxz = tmp_geo_data[2,*]
    geom.gyy = tmp_geo_data[3,*]
    geom.gyz = tmp_geo_data[4,*]
    geom.gzz = tmp_geo_data[5,*]
    geom.Bfield = tmp_geo_data[6,*]
    geom.dBdx = tmp_geo_data[7,*]
    geom.dBdy = tmp_geo_data[8,*]
    geom.dBdz = tmp_geo_data[9,*]
    geom.jacobian = tmp_geo_data[10,*]
    geom.R = tmp_geo_data[11,*]
    geom.phi = tmp_geo_data[12,*]
    geom.Z = tmp_geo_data[13,*]
    geom.c1 = tmp_geo_data[14,*]
    geom.c2 = tmp_geo_data[15,*]

  ENDIF ELSE BEGIN ; old tracer interface files
    tmp_geo_data = DBLARR(8,par.nz0,/NOZERO)

    READF, geo_lun, tmp_geo_data

    geom.gxx = tmp_geo_data[0,*]
    geom.gxy = tmp_geo_data[1,*]
    geom.Bfield = tmp_geo_data[3,*]
    ;tmp_geo_data[7,*] is dBds
    geom.jacobian = tmp_geo_data[6,*]
  ENDELSE

END

; ############################################################################

FUNCTION read_datakey, geo_lun, data_str, data
; it is not possible to return data as a parameter into a structure
; therefore we need this 'strange' function construction

  COMMON global_vars

  IF par.write_h5 THEN BEGIN
    geo_dset = H5D_OPEN(geo_lun,data_str)
    data = H5D_READ(geo_dset)
    H5D_CLOSE, geo_dset
  ENDIF ELSE BEGIN
    line = ''
    WHILE NOT EOF(geo_lun) AND STRPOS(line,data_str) LT 0 DO $
      READF, geo_lun, line
    IF EOF(geo_lun) THEN BEGIN
      printerror, 'Error while reading ' + data_str + '; using 0.0 instead'
      POINT_LUN, geo_lun, 0 ;rewind
      RETURN, 0
    ENDIF ELSE READF, geo_lun, data ;, FORMAT='(16E20.10)'
  ENDELSE

  RETURN, data

END

; ############################################################################

PRO read_tracer_df, geo_lun, geom

  COMMON global_vars

  IF par.write_h5 THEN BEGIN
    geom.q = read_datakey(geo_lun,'profile/q_prof')
    geom.gxx = read_datakey(geo_lun,'metric/g_xx')
    geom.gxy = read_datakey(geo_lun,'metric/g_xy')
    geom.gxz = read_datakey(geo_lun,'metric/g_xz')
    geom.gyy = read_datakey(geo_lun,'metric/g_yy')
    geom.gyz = read_datakey(geo_lun,'metric/g_yz')
    geom.gzz = read_datakey(geo_lun,'metric/g_zz')
    geom.C_y = read_datakey(geo_lun,'metric/C_y')
    geom.C_xy = read_datakey(geo_lun,'metric/C_xy')
    geom.Bfield = read_datakey(geo_lun,'Bfield_terms/Bfield')
    geom.dBdx = read_datakey(geo_lun,'Bfield_terms/dBdx')
    geom.dBdy = read_datakey(geo_lun,'Bfield_terms/dBdy')
    geom.dBdz = read_datakey(geo_lun,'Bfield_terms/dBdz')
    geom.jacobian = read_datakey(geo_lun,'Bfield_terms/Jacobian')
    IF (par.magn_geometry NE 'circular') THEN BEGIN
      geom.R = read_datakey(geo_lun,'shape/geo_R')
      geom.Z = read_datakey(geo_lun,'shape/geo_Z')
      geom.c1 = read_datakey(geo_lun,'shape/c_1')
      geom.c2 = read_datakey(geo_lun,'shape/c_2')
    ENDIF
  ENDIF ELSE BEGIN
    geom.q = read_datakey(geo_lun,'q',geom.q)
    geom.gxx = read_datakey(geo_lun,'gxx',geom.gxx)
    geom.gxy = read_datakey(geo_lun,'gxy',geom.gxy)
    geom.gxz = read_datakey(geo_lun,'gxz',geom.gxz)
    geom.gyy = read_datakey(geo_lun,'gyy',geom.gyy)
    geom.gyz = read_datakey(geo_lun,'gyz',geom.gyz)
    geom.gzz = read_datakey(geo_lun,'gzz',geom.gzz)
    geom.Bfield = read_datakey(geo_lun,'Bfield',geom.Bfield)
    geom.dBdx = read_datakey(geo_lun,'dBdx',geom.dBdx)
    geom.dBdy = read_datakey(geo_lun,'dBdy',geom.dBdy)
    geom.dBdz = read_datakey(geo_lun,'dBdz',geom.dBdz)
    geom.jacobian = read_datakey(geo_lun,'jacobian',geom.jacobian)
    geom.C_y = read_datakey(geo_lun,'C_y',geom.C_y)
    geom.C_xy = read_datakey(geo_lun,'C_xy',geom.C_xy)
    IF (par.magn_geometry NE 'circular') THEN BEGIN
      geom.R = read_datakey(geo_lun,'geo_R',geom.R)
      geom.Z = read_datakey(geo_lun,'geo_Z',geom.Z)
      geom.c1 = read_datakey(geo_lun,'geo_c1',geom.c1)
      geom.c2 = read_datakey(geo_lun,'geo_c2',geom.c2)
    ENDIF
  ENDELSE

END

; ############################################################################

PRO read_tracer, geo_lun, geom_struct
  COMMON global_vars
  
  CASE par.in_data OF
     'kxky' : read_tracer_ff, geo_lun, geom_struct
     'sxky' : read_tracer_df, geo_lun, geom_struct
     'sykx' : read_tracer_df, geo_lun, geom_struct
     ELSE : printerror, 'ERROR in read_tracer'
  ENDCASE
END

; ############################################################################

PRO read_geometry

  COMMON global_vars

  IF PTR_VALID(series.geom) THEN PTR_FREE, series.geom

  IF par.x_local AND par.y_local THEN geom_struct = {$
    q        : par.q0,$
    gxx      : DBLARR(par.nz0),$
    gxy      : DBLARR(par.nz0),$
    gxz      : DBLARR(par.nz0),$
    gyy      : DBLARR(par.nz0),$
    gyz      : DBLARR(par.nz0),$
    gzz      : DBLARR(par.nz0),$
    Bfield   : DBLARR(par.nz0),$
    dBdx     : DBLARR(par.nz0),$
    dBdy     : DBLARR(par.nz0),$
    dBdz     : DBLARR(par.nz0),$
    jacobian : DBLARR(par.nz0),$
    jac_norm : DBLARR(par.nz0),$
    C_y      : (par.minor_r*par.x0 GT 0)?$
                par.x0*par.minor_r/par.q0:1.0,$
    C_xy     : 1.0,$
    R	     : DBLARR(par.nz0),$
    z	     : DBLARR(par.nz0),$
    phi      : DBLARR(par.nz0),$
    c1	     : DBLARR(par.nz0),$
    c2	     : DBLARR(par.nz0)} $
  ELSE IF (par.x_local EQ 0) THEN $
     geom_struct = {$
      q        : DBLARR(par.nx0)+par.q0,$
      gxx      : DBLARR(par.nx0,par.nz0),$
      gxy      : DBLARR(par.nx0,par.nz0),$
      gxz      : DBLARR(par.nx0,par.nz0),$
      gyy      : DBLARR(par.nx0,par.nz0),$
      gyz      : DBLARR(par.nx0,par.nz0),$
      gzz      : DBLARR(par.nx0,par.nz0),$
      Bfield   : DBLARR(par.nx0,par.nz0),$
      dBdx     : DBLARR(par.nx0,par.nz0),$
      dBdy     : DBLARR(par.nx0,par.nz0),$
      dBdz     : DBLARR(par.nx0,par.nz0),$
      jacobian : DBLARR(par.nx0,par.nz0),$
      jac_norm : DBLARR(par.nx0,par.nz0),$
      C_y      : DBLARR(par.nx0),$
      C_xy     : DBLARR(par.nx0),$
      R	       : DBLARR(par.nx0,par.nz0),$
      z	       : DBLARR(par.nx0,par.nz0),$
      phi      : DBLARR(par.nx0,par.nz0),$
      c1       : DBLARR(par.nx0,par.nz0),$
      c2       : DBLARR(par.nx0,par.nz0)} $
  ELSE IF (par.y_local EQ 0) THEN BEGIN
     ind1 = par.ny0
     geom_struct = {$
      q        : par.q0,$
      gxx      : DBLARR(ind1,par.nz0),$
      gxy      : DBLARR(ind1,par.nz0),$
      gxz      : DBLARR(ind1,par.nz0),$
      gyy      : DBLARR(ind1,par.nz0),$
      gyz      : DBLARR(ind1,par.nz0),$
      gzz      : DBLARR(ind1,par.nz0),$
      Bfield   : DBLARR(ind1,par.nz0),$
      dBdx     : DBLARR(ind1,par.nz0),$
      dBdy     : DBLARR(ind1,par.nz0),$
      dBdz     : DBLARR(ind1,par.nz0),$
      jacobian : DBLARR(ind1,par.nz0),$
      jac_norm : DBLARR(ind1,par.nz0),$
      C_y      : 1.0,$
      C_xy     : 1.0,$
      R	       : DBLARR(ind1,par.nz0),$
      z	       : DBLARR(ind1,par.nz0),$
      phi      : DBLARR(ind1,par.nz0),$
      c1       : DBLARR(ind1,par.nz0),$
      c2       : DBLARR(ind1,par.nz0)}
  ENDIF

  use_s_alpha = 0
  stell_dat_format = 0

  IF par.write_h5 THEN BEGIN
    IF par.magn_geometry EQ 'tracer' THEN BEGIN
      geo_file = get_file_path((STRSPLIT(par.geomfile,'.dat',$
        /EXTRACT,/REGEX))[0],0)
      IF NOT geo_file.exist THEN BEGIN
        ; try without extension
        geo_file.path = gui.out.data_path + par.geomfile + '.h5'
        geo_file.exist = FILE_TEST(geo_file.path)
        IF NOT geo_file.exist THEN geo_file = $
          get_file_path((STRSPLIT(par.magn_geometry,'.dat',$
          /EXTRACT,/REGEX))[0],0)
      ENDIF
    ENDIF ELSE geo_file = get_file_path((STRSPLIT(par.magn_geometry,'.dat',$
      /EXTRACT,/REGEX))[0],0)
    err = 1 - H5F_IS_HDF5(geo_file.path[0])
    IF NOT err THEN BEGIN
      IF !QUIET NE 1 THEN PRINT, 'using ' + geo_file.path[0]
      geo_lun = H5F_OPEN(geo_file.path[0])
    ENDIF ELSE BEGIN
      IF (par.magn_geometry NE 's_alpha') AND (par.magn_geometry NE 's-alpha') THEN $
        printerror, geo_file.path[0] + ' not HDF5, using s-alpha instead'
      use_s_alpha = 1
    ENDELSE

    IF use_s_alpha THEN set_s_alpha, geom_struct ELSE BEGIN
      IF check_gridpoints(geo_lun,geo_file.path[0]) EQ 1 THEN BEGIN
         read_tracer, geo_lun, geom_struct
      ENDIF ELSE set_s_alpha, geom_struct
      H5F_CLOSE, geo_lun
    ENDELSE
  ENDIF ELSE BEGIN
    IF par.magn_geometry EQ 'stell.dat' THEN BEGIN ; very old stell.dat format
      run_labels = STRSPLIT(series.run_labels,',',/EXTRACT)
      file_ext = (run_labels[0] NE 'act') AND $
        (run_labels[0] NE 'dat') ? '_' + run_labels[0] : '.dat'
      OPENR, geo_lun, gui.out.data_path + 'stell_' + rm0es(par.nz0), $
        /GET_LUN, ERROR=ierr
      IF ierr NE 0 THEN BEGIN
        OPENR, geo_lun, gui.out.data_path + 'stell' + file_ext, $
          /GET_LUN, ERROR=ierr
        IF ierr NE 0 THEN BEGIN
          printerror, 'neither stell_' + rm0es(par.nz0) + ' nor stell' + $
            file_ext + ' file found, using s-alpha'
          use_s_alpha = 1
        ENDIF ELSE PRINT, 'using stell' + file_ext
      ENDIF ELSE PRINT, 'using stell_' + rm0es(par.nz0)
      line = ''
      READF, geo_lun, line
      IF STRMID(line,0,1) NE '&' THEN stell_dat_format = 1
      POINT_LUN, geo_lun, 0
    ENDIF ELSE BEGIN
      IF par.magn_geometry EQ 'tracer' THEN BEGIN
        geo_file = get_file_path((STRSPLIT(par.geomfile,'.dat',$
          /EXTRACT,/REGEX))[0],0)
        IF NOT geo_file.exist THEN BEGIN
          ; try without extension
          geo_file.path = gui.out.data_path + par.geomfile
          geo_file.exist = FILE_TEST(geo_file.path)
          IF NOT geo_file.exist THEN geo_file = $
            get_file_path((STRSPLIT(par.magn_geometry,'.dat',$
            /EXTRACT,/REGEX))[0],0)
        ENDIF
      ENDIF ELSE geo_file = get_file_path((STRSPLIT(par.magn_geometry,'.dat',$
        /EXTRACT,/REGEX))[0],0)
      OPENR, geo_lun,geo_file.path[0],/GET_LUN, ERROR=ierr
      IF ierr EQ 0 THEN PRINT, 'using ' + geo_file.path[0] ELSE BEGIN
        IF (par.magn_geometry NE 's_alpha') AND (par.magn_geometry NE 's-alpha') THEN $
          printerror, geo_file.path + ' not found, using s-alpha instead'
        use_s_alpha = 1
      ENDELSE
    ENDELSE

    IF use_s_alpha THEN set_s_alpha, geom_struct ELSE BEGIN
      IF stell_dat_format THEN read_stell_dat, geo_lun, geom_struct $
        ELSE BEGIN ; use newer format
        IF check_gridpoints(geo_lun,geo_file.path[0]) EQ 1 THEN BEGIN
           read_tracer, geo_lun, geom_struct
        ENDIF ELSE set_s_alpha, geom_struct
        FREE_LUN, geo_lun
      ENDELSE
    ENDELSE
  ENDELSE

  geom_struct.jac_norm = geom_struct.jacobian / $
    TOTAL(geom_struct.jacobian) * N_ELEMENTS(geom_struct.jacobian)
  ; add Bref since B enters e.g. in the transport calculation
  geom_struct.Bfield = geom_struct.Bfield * series.Bref
  series.geom = PTR_NEW(geom_struct)

END

; ############################################################################

PRO read_profiles, run=run

  COMMON global_vars

  FOR n = 0, par.n_spec - 1 DO BEGIN
    IF PTR_VALID(spec[n].prof) THEN PTR_FREE, spec[n].prof

    prof_struct = {$
      temp : DBLARR(par.nx0,/NOZERO),$
      dens : DBLARR(par.nx0,/NOZERO),$
      omt  : DBLARR(par.nx0,/NOZERO),$
      omn  : DBLARR(par.nx0,/NOZERO),$
      x    : DBLARR(par.nx0,/NOZERO)}

    prof_struct.temp = spec[n].temp
    prof_struct.dens = spec[n].dens
    prof_struct.omt = spec[n].omt
    prof_struct.omn = spec[n].omn

    IF KEYWORD_SET(run) THEN $
      file_prof = get_file_path('profiles_'+spec[n].name,run) $
    ELSE file_prof = get_file_path('profiles_'+spec[n].name)

    IF NOT file_prof.exist THEN BEGIN
      printerror, file_prof.path[0] + ' not found; using constant profiles'
    ENDIF ELSE BEGIN
      IF par.write_h5 THEN BEGIN
        ierr = 1 - H5F_IS_HDF5(file_prof.path[0])
        IF ierr EQ 0 THEN PRINT, 'using ' + file_prof.path[0] ELSE $
          printerror, file_prof.path[0] + ' not found; using constant profiles'
        prof_lun = H5F_OPEN(file_prof.path[0])

        prof_dset = H5D_OPEN(prof_lun,'/position/x_o_rho_ref')
        prof_struct.x = H5D_READ(prof_dset)
        H5D_CLOSE, prof_dset
        prof_dset = H5D_OPEN(prof_lun,'/density/n')
        prof_struct.dens = H5D_READ(prof_dset)
        H5D_CLOSE, prof_dset
        prof_dset = H5D_OPEN(prof_lun,'/density/omn')
        prof_struct.omn = H5D_READ(prof_dset)
        H5D_CLOSE, prof_dset
        prof_dset = H5D_OPEN(prof_lun,'/temp/T')
        prof_struct.temp = H5D_READ(prof_dset)
        H5D_CLOSE, prof_dset
        prof_dset = H5D_OPEN(prof_lun,'/temp/omt')
        prof_struct.omt = H5D_READ(prof_dset)
        H5D_CLOSE, prof_dset
      ENDIF ELSE BEGIN
        OPENR, prof_lun, file_prof.path[0], /GET_LUN, ERROR=ierr
        IF ierr EQ 0 THEN PRINT, 'using ' + file_prof.path[0] ELSE $
          printerror, file_prof.path[0] + ' not found; using constant profiles'

        dummy_str = ''
        READF, prof_lun, dummy_str

        dummy = STRSPLIT(dummy_str,' ',COUNT=n_cols)

        POINT_LUN, - prof_lun, pos
        READF, prof_lun, dummy_str
        IF STRMID(dummy_str,0,1) NE '#' THEN POINT_LUN, prof_lun, pos

        temp_data = DBLARR(n_cols-1,par.nx0)

        READF, prof_lun, temp_data

        IF NOT EOF(prof_lun) THEN printerror, 'possible error while reading ' + $
          file_prof.path[0]

        FREE_LUN, prof_lun

        IF (n_cols EQ 6) THEN BEGIN
          prof_struct.x = temp_data[0,*] / par.rhostar
          prof_struct.temp = temp_data[1,*]
          prof_struct.dens = temp_data[2,*]
          prof_struct.omt  = temp_data[3,*]
          prof_struct.omn  = temp_data[4,*]
        ENDIF ELSE BEGIN
          prof_struct.x = temp_data[1,*]
          prof_struct.temp = temp_data[2,*]
          prof_struct.dens = temp_data[3,*]
          prof_struct.omt  = temp_data[4,*]
          prof_struct.omn  = temp_data[5,*]
        ENDELSE
      ENDELSE
    ENDELSE

    ; check whether profiles are normalized to reference position x0
    ; and correct if necessary
    mindist = MIN(ABS(prof_struct.x*par.rhostar-par.x0),x0_ind)
    allowed_diff=prof_struct.omt[x0_ind]*prof_struct.temp[x0_ind]*mindist
    ;with resets, the profiles are indeed not normalized to x0
    IF (par.reset_limit LT 1e-5) then begin
       IF ((prof_struct.temp[x0_ind] GT 1.0+allowed_diff) OR $
           (prof_struct.temp[x0_ind] LT 1.0-allowed_diff)) THEN $
              prof_struct.temp *= (1E3/par.Tref/spec[n].temp)
       
       IF ((prof_struct.dens[x0_ind] GT 1.0+allowed_diff) OR $
           (prof_struct.dens[x0_ind] LT 1.0-allowed_diff)) THEN $
              prof_struct.dens *= (1E19/par.nref/spec[n].dens)
    endif 
    spec[n].prof = PTR_NEW(prof_struct)
  ENDFOR ;--- n loop
END
