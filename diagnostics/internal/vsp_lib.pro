;##########################################################################
;# this library contains almost all vsp related internal functions        #
;##########################################################################

FUNCTION vsp_renorm, var, Lref=Lref

  COMMON global_vars

  IF NOT KEYWORD_SET(Lref) THEN Lref = series.Lref
  v_ref = SQRT(series.Tref*series.Qref/series.mref) ; T in eV

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref
  D_gb = v_ref * rho_ref^2 / Lref

  ; #### THIS NEEDS TO BE VERIFIED ####
  CASE var MOD 8 OF ; renormalization is the same for each species
    0 : RETURN, 1.0D/(D_gb*series.nref/Lref)             ; G_es
    1 : RETURN, 1.0D/(D_gb*series.nref/Lref)             ; G_em
    2 : RETURN, 1.0D/(D_gb*series.nref*series.Tref/Lref) ; Q_es
    3 : RETURN, 1.0D/(D_gb*series.nref*series.Tref/Lref) ; Q_em
    4 : RETURN, 1.0D/(series.nref/v_ref^3*rho_ref/Lref)  ; SQRT(<|f_|^2>)
    ELSE : printerror, 'error in vsp_renorm'
  ENDCASE

END

;#############################################################################

FUNCTION get_vsp_string, var_index, fancy=fancy, units=units, ounit = ounit, $
  time=time, flux=flux, normstr=normstr, rhostr=rhostr
; takes an index or an array of indices and returns the corresponding
; label(s); if /units is specified, fancy versions are returned
; ounit: returns only the unit in fancy style

  COMMON global_vars

  rho_sc = 0
  isp = 0
  WHILE (isp LT par.n_spec) AND (rho_sc EQ 0) DO BEGIN
    IF ((spec[isp].charge EQ -1) AND $
      (spec[isp].temp/series.Tref EQ 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN rho_sc = 3
    IF ((spec[isp].charge EQ 1) AND $
      (spec[isp].temp/series.Tref EQ 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN rho_sc = 2
    IF ((spec[isp].charge EQ 1) AND $
      (spec[isp].temp/series.Tref NE 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN rho_sc = 1
    isp = isp + 1
  ENDWHILE

  ref_str = 'ref'
  CASE rho_sc OF
    1 : BEGIN
          vel_str = 'c!Ds!N'
          sp_str = 'i'
          rho_str = '!7q!6!Ds!N'
        END
    2 : BEGIN
          vel_str = 'v!Dti!N'
          sp_str = 'i'
          rho_str = '!7q!6!Di!N'
        END
    3 : BEGIN
          vel_str = 'v!Dte!N'
          sp_str = 'e'
          rho_str = '!7q!6!De!N'
        END
    ELSE : BEGIN
             vel_str = 'v!Dt,' + ref_str + '!N'
             sp_str = ''
             rho_str = '!7q!6!D' + ref_str + '!N'
           END
  ENDCASE

  Lref_str = series.Lref_str

  v = STRARR(5)
  vf = STRARR(5)
  vu = STRARR(5)
  v[0] = '<G_es>'
  v[1] = '<G_em>'
  v[2] = '<Q_es>'
  v[3] = '<Q_em>'
  v[4] = '<|f|^2>^(1/2)'
  vf[0] = '!7C!6!Des!N'
  vf[1] = '!7C!6!Dem!N'
  vf[2] = '!6Q!Des!N'
  vf[3] = '!6Q!Dem!N'
  vf[4] = '!9!!!6f!9!!!6!U2!N'

  vf = '!12<' + vf + '!12>!6!N'

  vf[4] = vf[4] + '!U1/2!N'

  sp_str = 'j' ; n and T normalized to n0j * nref not to nref
  vu[0:1] = '' ; vel_str+rho_str+'!U2!N'+'n!De0!N/'+Lref_str+'!U2!N'
  vu[2:3] = '' ; vel_str+rho_str+'!U2!N'+'p!Dref0!N/'+Lref_str+'!U2!N'
  vu[4] = ''

  IF KEYWORD_SET(ounit) THEN RETURN, vu[var_index]
  IF KEYWORD_SET(units) THEN RETURN, vf[var_index] + ' / (' + vu[var_index] + ')'
  IF KEYWORD_SET(fancy) THEN RETURN, vf[var_index] ELSE $
    RETURN, v[var_index]

END

;#############################################################################

FUNCTION get_vsp_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  start_time = gui.out.start_t

  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
  file_vsp = get_file_path('vsp',run)

  IF NOT file_vsp.exist THEN RETURN, 0.0 ;-1L

  IF par.write_h5 THEN BEGIN
    err = 1 - H5F_IS_HDF5(file_vsp.path)
    IF NOT err THEN vsp_lun = H5F_OPEN(file_vsp.path)
  ENDIF ELSE $
    OPENR, vsp_lun, file_vsp.path, /GET_LUN, ERROR=err, /F77_UNFORMATTED, $
    SWAP_ENDIAN=series.swap_endian

  IF err NE 0 THEN BEGIN
    printerror, 'file error occured in get_vsp_time', /popup
    RETURN, 0L
  ENDIF

  time = par.prec_single ? 0.0 : 0.0D
  time_first_step = time
  time_second_step = time
  time_last_step = time

  IF par.write_h5 THEN BEGIN
    h5_dir_str = 'vsp/time'
    time_dset = H5D_OPEN(vsp_lun,h5_dir_str)
    time_arr = H5D_READ(time_dset)
    time_arr *= time_renorm(0)

    n_step_attrib = H5A_OPEN_NAME(time_dset,'n_steps')
    n_steps = H5A_READ(n_step_attrib)
    H5A_CLOSE, n_step_attrib
    H5D_CLOSE, time_dset

    IF (SIZE(time_arr))[1] NE n_steps THEN BEGIN
      printerror, 'invalid time array for HDF5 get_vsp_time'
      RETURN, 0L
    ENDIF

    series.h5_nst_vsp = n_steps

    IF KEYWORD_SET(get_n_steps) THEN BEGIN
      H5F_CLOSE, vsp_lun
      RETURN, n_steps
    ENDIF

    time_first_step[0] = time_arr[0]
    time_last_step[0] = time_arr[n_steps-1]

    IF KEYWORD_SET(first) THEN BEGIN
      H5F_CLOSE, vsp_lun
      RETURN, time_first_step
    ENDIF

    IF KEYWORD_SET(last) THEN BEGIN
      H5F_CLOSE, vsp_lun
      RETURN, time_last_step
    ENDIF

    IF time_first_step GT start_time THEN BEGIN
      H5F_CLOSE, vsp_lun
      RETURN, 0L
    ENDIF

    IF time_last_step LT start_time THEN BEGIN
      IF !QUIET NE 1 THEN PRINT, 'run ', run_label, ' not in time range'
      H5F_CLOSE, vsp_lun
      RETURN, -2L
    ENDIF

    IF time_last_step EQ time_first_step THEN BEGIN
      PRINT, 'warning: run ', run_label, ' has only one time step'
      H5F_CLOSE, vsp_lun
      RETURN, 0L
    ENDIF

    step_ind = (WHERE(FLOAT(time_arr) GE FLOAT(start_time)))[0] > 0L
    step_time = time_arr[step_ind]

    H5F_CLOSE, vsp_lun

    RETURN, step_ind
  ENDIF ELSE BEGIN
    vsp_size = (FSTAT(vsp_lun)).SIZE

    content_byte_size = 4LL * (2 - par.prec_single) * par.nz0 * par.nv0 * $
      par.nw0 * par.n_spec * 5 + 8
    step_byte_size = content_byte_size + 16 - 4 * par.prec_single

    pos_last_step = (FSTAT(vsp_lun)).SIZE - step_byte_size
    n_steps = ROUND((FSTAT(vsp_lun)).SIZE/step_byte_size)
    IF KEYWORD_SET(get_n_steps) THEN BEGIN
      FREE_LUN, vsp_lun
      RETURN, n_steps
    ENDIF
    READU, vsp_lun, time_first_step
    time_first_step *= time_renorm(0)

    IF KEYWORD_SET(first) THEN BEGIN
      FREE_LUN, vsp_lun
      RETURN, time_first_step
    ENDIF

    POINT_LUN, - vsp_lun, position_before
    position_after = position_before + content_byte_size
    POINT_LUN, vsp_lun, position_after

    IF NOT EOF(vsp_lun) THEN BEGIN
      READU, vsp_lun, time_second_step
      time_second_step *= time_renorm(0)
    ENDIF

    POINT_LUN, vsp_lun, pos_last_step
    READU, vsp_lun, time_last_step
    time_last_step *= time_renorm(0)

    IF KEYWORD_SET(last) THEN BEGIN
      FREE_LUN, vsp_lun
      RETURN, time_last_step
    ENDIF
    avg_step_diff = (time_last_step - time_first_step) / FLOAT(n_steps)

    IF time_first_step GT start_time THEN BEGIN
      FREE_LUN, vsp_lun
      RETURN, 0L
    ENDIF

    IF time_last_step LT start_time THEN BEGIN
      IF !QUIET NE 1 THEN PRINT, 'run ', run_label, ' not in time range'
      FREE_LUN, vsp_lun
      RETURN, -2L
    ENDIF

    IF time_last_step EQ time_first_step THEN BEGIN
      PRINT, 'warning: run ', run_label, ' has only one time step'
      FREE_LUN, vsp_lun
      RETURN, 0L
    ENDIF

    IF time_second_step GT start_time THEN BEGIN
      FREE_LUN, vsp_lun
      RETURN, 0L
    ENDIF

    ; get next step:
    estimated_step = ROUND((start_time-time_first_step)/avg_step_diff)
    estimated_pos_init = step_byte_size * estimated_step - $
      step_byte_size < (pos_last_step - step_byte_size)

    found = 0
    iterations = 0
    max_iterations = 20
    estimated_pos = estimated_pos_init

    WHILE (found EQ 0) AND (iterations LE max_iterations) AND (estimated_pos GE 0) DO BEGIN
      IF estimated_pos GE 0 THEN POINT_LUN, vsp_lun, estimated_pos ELSE BEGIN
        PRINT, 'estimated_pos LT 0, file position error!'
        FREE_LUN, vsp_lun
        RETURN, -1L
      ENDELSE
      time_prev = par.prec_single ? 0.0 : 0.0D
      READU, vsp_lun, time_prev
      time_prev = time_prev  * time_renorm(0)
      POINT_LUN, - vsp_lun, position_before
      position_after = position_before + content_byte_size
      POINT_LUN, vsp_lun, position_after
      time_act = par.prec_single ? 0.0 : 0.0D
      READU, vsp_lun, time_act
      time_act = time_act  * time_renorm(0)
      time_diff_act = time_act - time_prev
      IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
      new_estimated_pos = step_byte_size * $
        ROUND((start_time-time_act)/time_diff_act) + estimated_pos
      IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
        new_estimated_pos = new_estimated_pos + step_byte_size
      estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)
      IF estimated_pos LT 0 THEN estimated_pos = 0L

      iterations = iterations + 1
    ENDWHILE
    IF found EQ 0 THEN BEGIN
      PRINT, 'using simplified search pattern'
      estimated_pos = estimated_pos_init
      WHILE found EQ 0 DO BEGIN
        IF estimated_pos GE 0 THEN POINT_LUN, vsp_lun, estimated_pos ELSE BEGIN
          PRINT, 'estimated_pos LT 0, file position error!'
          FREE_LUN, vsp_lun
          RETURN, -1L
        ENDELSE
        READU, vsp_lun, time_prev
        time_prev = time_prev  * time_renorm(0)
          POINT_LUN, - vsp_lun, position_before
          position_after = position_before + content_byte_size
          POINT_LUN, vsp_lun, position_after
        READU, vsp_lun, time_act
        time_act = time_act  * time_renorm(0)
        time_diff_act = time_act - time_prev
        IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
        new_estimated_pos = step_byte_size * $
          ROUND((start_time-time_act)/time_diff_act) / ABS(ROUND((start_time-time_act)/$
          time_diff_act)) + estimated_pos
        IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
          new_estimated_pos = new_estimated_pos + step_byte_size
        estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)
        IF estimated_pos LT 0 THEN estimated_pos = 0L
        IF found EQ 0 THEN PRINT, 'new estimate:', estimated_pos, (pos_last_step - step_byte_size)
      ENDWHILE
    ENDIF ;ELSE PRINT, 'number of iterations:', iterations

    FREE_LUN, vsp_lun

    step_pos = ROUND(estimated_pos/(1.0*step_byte_size))

    RETURN, step_pos
  ENDELSE

END

;##########################################################################

PRO jump_to_vsp_step, jump_step, vsp_lun
; jumps to the desired step of the vsp file

  COMMON global_vars

  IF jump_step EQ -1 THEN BEGIN
    PRINT, 'omitting run, not in time range'
    RETURN
  ENDIF

  IF par.write_h5 THEN BEGIN
    series.h5_count_vsp = jump_step
    RETURN
  ENDIF

  ; get step byte size, jump to position
  vsp_step_bytes = 0
  time = par.prec_single ? 0.0 : 0.0D
  content = par.prec_single ? $
    FLTARR(par.nz0,par.nv0,par.nw0,gui.out.n_spec_sel,5,/NOZERO) : $
    DBLARR(par.nz0,par.nv0,par.nw0,gui.out.n_spec_sel,5,/NOZERO)
  READU, vsp_lun, time
  READU, vsp_lun, content
  POINT_LUN, - vsp_lun, vsp_step_bytes
  jump_byte = LONG64(vsp_step_bytes) * jump_step
  POINT_LUN, vsp_lun, jump_byte

END

;##########################################################################

PRO read_vsp_step, field_lun, mom_luns

END

;##########################################################################

PRO vsp_loop, diag0

  COMMON global_vars

  run_labels = STRSPLIT(series.run_labels,',',/EXTRACT)
  file_vsp = get_file_path('vsp')

  IF (file_vsp.exist EQ 0) THEN BEGIN
    printerror, 'Skipping vsp loop due to missing vsp file(s)'
    RETURN
  ENDIF

  content_byte_size = 4LL * (2 - par.prec_single) * $
    par.nz0 * par.nv0 * par.nw0 * par.n_spec * 5 + 8

  ; mapping field name to number
  IF par.write_h5 THEN vsp_map = ['G_es','G_em','Q_es','Q_em','<f_>']

  FOR run = 0, N_ELEMENTS(run_labels) - 1 DO BEGIN
    run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]

    IF par.write_h5 THEN BEGIN
      vsp_err = 1 - H5F_IS_HDF5(file_vsp.path[run])
      IF NOT vsp_err THEN vsp_lun = H5F_OPEN(file_vsp.path[run])

      series.h5_count_field = 0
    ENDIF ELSE OPENR, vsp_lun, file_vsp.path[run], /GET_LUN, $
      ERROR=vsp_err, /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian

    IF vsp_err NE 0 THEN BEGIN
      printerror, file_vsp.path, ' does not exist'
      RETURN
    ENDIF

    ; the slowest file with requested data determines jump_step
    jump_step = get_vsp_time(run)
    IF jump_step LT 0 THEN BEGIN
      end_loop = 1
      IF !QUIET NE 1 THEN PRINT, 'omitting run ', run_label
    ENDIF ELSE BEGIN
      jump_to_vsp_step, jump_step, vsp_lun
      IF jump_step GT 0 THEN $
        PRINT, '- starting vsp loop at step ', rm0es(jump_step), ' -'
      IF jump_step EQ 0 THEN PRINT, '- starting vsp loop -'
      IF !QUIET EQ 1 THEN PRINT, 'reading run ', run_label, ': vsp' + $
        (par.write_h5 ? ' (HDF5)' : '')
      end_loop = 0
    ENDELSE

    vsp_time = par.prec_single ? 1.0 : 1.0D
    end_t = par.prec_single ? FLOAT(gui.out.end_t) : gui.out.end_t
    vsp_time *= end_t ; to enter WHILE loop

    eof_vsp = par.write_h5 ? $
      (series.h5_nst_vsp LE series.h5_count_vsp) : EOF(vsp_lun)
    IF eof_vsp THEN BEGIN
      IF par.write_h5 THEN H5F_CLOSE, vsp_lun ELSE FREE_LUN, vsp_lun
      vsp_lun = 0
    ENDIF

    WHILE NOT eof_vsp AND (vsp_time LE end_t) $
      AND NOT end_loop DO BEGIN

      IF par.write_h5 THEN BEGIN
        ; take time from first field
        step_nr = STRING(series.h5_count_vsp,FORMAT='(I010)')
        step_loc = '/vsp/' + vsp_map[0] + '/' + step_nr + '/'
        step_dset = H5D_OPEN(vsp_lun,step_loc)
        time_attrib = H5A_OPEN_NAME(step_dset,'time')
        vsp_time = H5A_READ(time_attrib)
        H5A_CLOSE, time_attrib
        H5D_CLOSE, step_dset
      ENDIF ELSE BEGIN
        READU, vsp_lun, vsp_time
        vsp_time *= time_renorm(0)
      ENDELSE

      IF (vsp_time GE gui.out.start_t - 1e-6) AND $
        (vsp_time LE gui.out.end_t + 1e-6) THEN BEGIN

        PRINT, 'time: ', rm0es(vsp_time)

        series.step_count += 1L

        vsp_data = par.prec_single ? $
          FLTARR(par.nz0,par.nv0,par.nw0,par.n_spec,5,/NOZERO) : $
          DBLARR(par.nz0,par.nv0,par.nw0,par.n_spec,5,/NOZERO)

        IF par.write_h5 THEN BEGIN
          FOR v = 0, 4 DO BEGIN
            step_nr = STRING(series.h5_count_vsp,FORMAT='(I010)')
            step_loc = '/vsp/' + vsp_map[v] + '/' + step_nr + '/'
            vsp_dset = H5D_OPEN(vsp_lun,step_loc)

            vsp_data_field = H5D_READ(vsp_dset)
            H5D_CLOSE, vsp_dset

            vsp_data[0,0,0,0,v] = TEMPORARY(vsp_data_field)
          ENDFOR

          series.h5_count_vsp += 1
        ENDIF ELSE READU, vsp_lun, vsp_data

        vsp_data = (TEMPORARY(vsp_data))[*,*,*,*gui.out.spec_select,*]

        IF cfg.kpar_filter GE 0 THEN BEGIN
          IF (cfg.kpar_filter LT par.nz0 / 2) AND (par.nz0 GT 1) THEN BEGIN
            vsp_data_kpar = FFT(vsp_data,DIMENSION=1)
            IF cfg.kpar_filter EQ 0 THEN BEGIN
              vsp_data_kpar[1:*,*,*,*,*] = 0
            ENDIF ELSE BEGIN
              neg_kpar = par.nz0 - cfg.kpar_filter
              vsp_data_kpar[0:cfg.kpar_filter-1,*,*,*,*] = 0
              vsp_data_kpar[cfg.kpar_filter+1:neg_kpar-1,*,*,*,*] = 0
              IF cfg.kpar_filter NE 1 THEN vsp_data_kpar[neg_kpar+1:*,*,*,*,*] = 0
            ENDELSE
            vsp_data[0,0,0,0,0] = FFT(vsp_data_kpar,DIMENSION=1,/INVERSE)
          ENDIF
        ENDIF

	FOR v = 0, 3 DO vsp_data[0,0,0,0,v] = $
          vsp_data[*,*,*,*,v] * vsp_renorm(v)
    	vsp = PTR_NEW(vsp_data,/NO_COPY)

    	call_diags, diag0, 'loop'

        PTR_FREE, vsp
      ENDIF ELSE BEGIN
        IF par.write_h5 THEN BEGIN
          series.h5_count_vsp += 1
        ENDIF ELSE BEGIN
          POINT_LUN, - vsp_lun, position_before
          position_after = position_before + content_byte_size
          POINT_LUN, vsp_lun, position_after
        ENDELSE
      ENDELSE

      eof_vsp = par.write_h5 ? $
        (series.h5_nst_vsp LE series.h5_count_vsp) : EOF(vsp_lun)
      IF eof_vsp THEN BEGIN
        IF par.write_h5 THEN H5F_CLOSE, vsp_lun ELSE FREE_LUN, vsp_lun
        vsp_lun = 0
      ENDIF
    ENDWHILE

;    IF par.write_h5 THEN H5F_CLOSE, vsp_lun ELSE FREE_LUN, vsp_lun
  ENDFOR

END

;##########################################################################

PRO GetMuWeightsAndKnots, muweights, muknots, lw=lw, nw0=nw0,$
  nw_per_block=nw_per_block,mu_grid_type=mu_grid_type
  COMMON global_vars
  
  IF NOT KEYWORD_SET(lw) THEN lw = DOUBLE(par.lw)
  IF NOT KEYWORD_SET(nw0) THEN nw0 = par.nw0
  IF NOT KEYWORD_SET(nw_per_block) THEN nw_per_block = 4
  IF NOT KEYWORD_SET(mu_grid_type) THEN $
     mu_grid_type = STRTRIM(par.mu_grid_type)
  
  muweights = DBLARR(nw0)
  muknots = DBLARR(nw0)

  dbl_lb = 0.0d0
  CASE mu_grid_type OF
     'equidist' : BEGIN
        deltamu=lw/nw0
        FOR i=0,nw0-1 DO BEGIN
           muknots[i]=(i+0.5)*deltamu
           muweights[i]=deltamu
        ENDFOR
     END
    'eq_vperp' : BEGIN
        deltamu=lw/(nw0*nw0)
        FOR i=0,nw0-1 DO BEGIN
          muknots[i]=(i-0.5)*(i-0.5)*deltamu
          muweights[i]= (2*i-1)*deltamu
        ENDFOR
     END
    'gau_leg' : BEGIN
       gauleg,dbl_lb,lw,muknots,muweights,nw0
       ; rest of the integration is added to the last weight
       muweights[nw0-1]+= exp(muknots[nw0-1]-lw)
    END
    'gau_leg_block' : BEGIN
       ;lw subblocks with equal widths containing nw_per_block
       ;grid points
       ndomain = nw0/nw_per_block
       if ((nw0 MOD nw_per_block) NE 0) THEN $
          printerror,'nw_in not divisible by nw_per_block'
       dbl_ub = 0.0d0
       FOR i=0,ndomain-1 DO BEGIN
          dbl_lb=dbl_ub
          dbl_ub=lw/ndomain*(i+1)
          gauleg, dbl_lb, dbl_ub, dbl_knots, dbl_weights, nw_per_block
          muweights[i*nw_per_block:(i+1)*nw_per_block-1]=dbl_weights
          muknots[i*nw_per_block:(i+1)*nw_per_block-1]=dbl_knots
       ENDFOR
    END
    'gau_leg_block_lin' : BEGIN
       ;lw subblocks with lin. increasing widths containing 
       ;nw_per_block grid points
       ndomain = nw0/nw_per_block
       if ((nw0 MOD nw_per_block) NE 0) THEN $
          printerror,'nw_in not divisible by nw_per_block'
       sub_lw_norm = 0.0d0
       FOR i=0,ndomain-1 DO sub_lw_norm += DOUBLE(i+1)
       dbl_ub = 0.0d0
       FOR i=0,ndomain-1 DO BEGIN
          dbl_lb=dbl_ub
          dbl_ub=dbl_ub+lw*DOUBLE(i+1)/sub_lw_norm
          gauleg, dbl_lb, dbl_ub, dbl_knots, dbl_weights, nw_per_block
          muweights[i*nw_per_block:(i+1)*nw_per_block-1]=dbl_weights
          muknots[i*nw_per_block:(i+1)*nw_per_block-1]=dbl_knots
       ENDFOR
       ; rest of the integration is added to the last weight
       muweights[nw0-1]+= exp(muknots[nw0-1]-lw)
    END
    'gau_lag' : BEGIN
       if (nw0 GE 124) THEN printerror, $
          'Gauss-Laguerre integration is currently not possible for nw0>124.'
       dbl_ub=lw
       gaulag, dbl_lb, dbl_ub, muknots, muweights, nw0
    END
    ELSE : printerror, 'Unkown mu_grid_type: '+mu_grid_type+'!'
 ENDCASE
END

;#############################################################################

;>Compute Gauss-Legendre abscissas and weights in interval [x1, x2]
;!Original Fortran program written by T-M Tran (CRPP, Lausanne)
;!\param n size of weights/knots array
;!\param x1 lower bound
;!\param x2 upper bound
;!\param x  Gauss-Legendre knots
;!\param w  Gauss-Legendre weights
PRO gauleg, x1,x2,x,w,n
  x = DBLARR(n)
  w = DBLARR(n)
  m=(n+1)/2
  xm=0.5d0*(DOUBLE(x2)+DOUBLE(x1))
  xl=0.5d0*(DOUBLE(x2)-DOUBLE(x1))
  eps = (machar(/DOUBLE)).EPS
  FOR i=1,m DO BEGIN
     z=COS(3.14159265358979323846d0*(i-0.25d0)/(n+0.5d0))
     REPEAT BEGIN
        p1=1.d0
        p2=0.d0
        FOR j=1,n DO BEGIN
           p3=p2
           p2=p1
           p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        ENDFOR
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
     ENDREP UNTIL ( ABS(z-z1) LT EPS )

     x[i-1]=xm-xl*z
     x[n+1-i-1]=xm+xl*z
     w[i-1]=2.d0*xl/((1.d0-z*z)*pp*pp)
     w[n+1-i-1]=w[i-1]
  ENDFOR
END

;#############################################################################
;> This routine is adapted from a F77 specfun routine.
PRO gaulag, x1,x2,X,W,N
;       =========================================================
;       Purpose : Compute the zeros of Laguerre polynomial Ln(x)
;                 in the interval [0,oo], and the corresponding
;                 weighting coefficients for Gauss-Laguerre
;                 integration
;       Input :   n    --- Order of the Laguerre polynomial
;                 x1   --- lower bound
;                 x2   --- upper bound
;                 X(n) --- Zeros of the Laguerre polynomial
;                 W(n) --- Corresponding weighting coefficients
;       =========================================================
  X = DBLARR(n)
  W = DBLARR(n)
    
  HN=1.0D0/n
  PF=0.0D0
  PD=0.0D0
  FOR  nr=1,n DO BEGIN
     Z=HN
     if (nr GT 1) THEN Z=X[nr-2]+HN*nr^1.27D0
     it=0
     REPEAT BEGIN
        it+=1
        Z0=Z
        P=1.0D0
        FOR i=1,nr-1 DO P*=(Z-X[i-1])
        F0=1.0D0
        F1=1.0D0-Z
        FOR k=2,n DO BEGIN
           PF=((2.0D0*k-1.0D0-Z)*F1-(k-1.0D0)*F0)/k
           PD=k/Z*(PF-F1)
           F0=F1
           F1=PF
        ENDFOR
        FD=PF/P
        Q=0.0D0
        FOR i=1,nr-1 DO BEGIN
           WP=1.0D0
           FOR j=1,nr-1 DO BEGIN
              if (j NE i) THEN WP*=(Z-X[j-1])
           ENDFOR
           Q=Q+WP
        ENDFOR
        GD=(PD-Q*FD)/P
        Z=Z-FD/GD
     ENDREP UNTIL ((it GT 40) OR (ABS((Z-Z0)/Z) LE 1.0D-15))
     X[nr-1]=Z
     W[nr-1]=1.0D0/(Z*PD*PD)
  ENDFOR
  w=w*exp(x)
  fac=x2/TOTAL(w,/DOUBLE)
  x=x*fac
  w=w*fac
END

;#############################################################################

FUNCTION get_fm, v_par, mu
; returns the Maxwellian distribution function (currently, only for the
; local code)

  COMMON global_vars

  v_abs = REBIN(REFORM(mu,[1,1,par.nw0]),[par.nz0,par.nv0,par.nw0]) * $
    REBIN((*series.geom).Bfield,[par.nz0,par.nv0,par.nw0]) + $
    REBIN(REFORM(v_par,[1,par.nv0]),[par.nz0,par.nv0,par.nw0])
  RETURN, (!DPI)^(-1.5) * EXP(-vabs^2)

END
