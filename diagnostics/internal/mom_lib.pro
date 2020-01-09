;##########################################################################
;# this library contains almost all mom/field related internal functions  #
;##########################################################################

FUNCTION mom_renorm, var, Lref=Lref

  COMMON global_vars

  IF NOT KEYWORD_SET(Lref) THEN Lref = series.Lref

  IF (par.diag_trap_levels GT 0) THEN BEGIN
     IF ((var GT par.n_fields + par.n_moms/3) AND (var LT 100)) THEN $
        var = ((var - par.n_fields) MOD par.n_moms/3) + par.n_fields     
  ENDIF

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref
  c_ref = SQRT(series.Tref*series.Qref/series.mref)
  p_ref = series.nref * series.Tref * series.Qref

  nf = par.n_fields
  CASE var OF
    199  : RETURN, (series.nref*series.Qref)                        ; sum_j q_j n_j
    104  : RETURN, (series.nref/ Lref)                              ; -d/dx n_1
    103  : RETURN, (series.Tref/ Lref)                              ; -d/dx T_1
    102  : RETURN, (series.Tref * (series.Qref / series.Qref) / $   ; d^2/dx^2 phi
             (rho_ref * Lref))
    101  : RETURN, (p_ref * rho_ref / Lref)                         ; p_j0
    100  : RETURN, (series.Tref * rho_ref / Lref)                   ; T_tot
    nf+8 : RETURN, (p_ref * c_ref^2 / series.Bref) * rho_ref / Lref ; TperpI1
    nf+7 : RETURN, (p_ref * c_ref^2 / series.Bref) * rho_ref / Lref ; TparI1
    nf+6 : RETURN, (p_ref / series.Bref) * rho_ref / Lref           ; densI1
    nf+5 : RETURN, par.norm08 ? (c_ref * rho_ref / Lref) : $        ; upar
             (c_ref * rho_ref / Lref) / SQRT(par.ehat)
    nf+4 : RETURN, par.norm08 ? (p_ref * c_ref * rho_ref / Lref) : $; qper
    	     (p_ref * c_ref * rho_ref / Lref) / SQRT(par.ehat)
    nf+3 : RETURN, par.norm08 ? (p_ref * c_ref * rho_ref / Lref) : $; qpar
    	     (p_ref * c_ref * rho_ref / Lref) / SQRT(par.ehat)
    nf+2 : RETURN, (series.Tref * rho_ref / Lref)                   ; Tperp
    nf+1 : RETURN, (series.Tref * rho_ref / Lref)                   ; Tpar
    nf   : RETURN, (series.nref * rho_ref / Lref)                   ; n
    2	 : RETURN, (series.Bref * rho_ref / Lref)                   ; Bpar
    1	 : RETURN, par.norm08 ? (series.Bref * rho_ref^2 / Lref) : $; Apar
    	     0.5 * par.beta * SQRT(par.ehat) * $
             (series.Bref * rho_ref^2 / Lref)
    0	 : RETURN, (series.Tref * (series.Qref / series.Qref) * $   ; phi
             rho_ref / Lref)
    ELSE : printerror, 'error in mom_renorm (var ' + rm0es(var) + $
             ' not found)', /popup
  ENDCASE

END

;#############################################################################

PRO create_mom_struct

  COMMON global_vars

  IF TOTAL(series.request_mom) EQ 0 THEN RETURN
  
  request_dat = series.request_mom  ;structure tags cannot be passed & modified
  check_request_dat, request_dat
  series.request_mom = request_dat

  IF !QUIET NE 1 THEN BEGIN
    mem_use = 16L * par.nky0 * par.nz0 / 1048576.0 * (0L + $
      TOTAL(request_dat[*,0]) * par.nkx0 + $
      TOTAL(request_dat[*,1]) * par.nkx0 + $
      TOTAL(request_dat[*,2]) * (par.nkx0 / 2 + 1) * 2 + $
      TOTAL(request_dat[*,3]) * par.nkx0) * $
      N_ELEMENTS(*gui.out.spec_select)
      PRINT, '    creating data structures (required memory: ' + $
         rm0es(mem_use) + ' MB)'
      PRINT, '    currently in use: ' + mem_usage()
      stime = SYSTIME(1)
  ENDIF

  mom = REPLICATE({kxky:PTR_NEW(),kxsy:PTR_NEW(),$
    sxky:PTR_NEW(),sxsy:PTR_NEW()},$
    N_ELEMENTS(*gui.out.spec_select),200) ;par.n_fields+par.n_moms)

  FOR isp = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO BEGIN
     tmpptr = REFORM(mom[isp,*])
     create_3darr_struct,tmpptr,series.request_mom
     mom[isp,*] = tmpptr
  ENDFOR

  IF !QUIET NE 1 THEN BEGIN
    PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'
    PRINT, '    current memory usage: ' + mem_usage()
  ENDIF

END

;##########################################################################

PRO destroy_mom_struct

  COMMON global_vars

  PTR_FREE, mom.kxky, mom.sxky, mom.kxsy, mom.sxsy

END

;##########################################################################

PRO fft_format, kxky=kxky, sxky=sxky, kxsy=kxsy, sxsy=sxsy, $
  doppler_freq=doppler_freq
; Called by the _init part of a diagnostic to indicate that it needs
; fourier or real space data.
; request_mom: 21,4 array (field,transformation type), 1 = reqd, 0 = not reqd
; Optionally, the fourier data might be modified with an additional
; phase factor correcting for a Doppler shift due to ExB shear flows

  COMMON global_vars

  FOR i = 0, N_ELEMENTS(kxky) - 1 DO series.request_mom[kxky[i],0] = 1
  FOR i = 0, N_ELEMENTS(sxky) - 1 DO series.request_mom[sxky[i],1] = 1
  FOR i = 0, N_ELEMENTS(kxsy) - 1 DO series.request_mom[kxsy[i],2] = 1
  FOR i = 0, N_ELEMENTS(sxsy) - 1 DO series.request_mom[sxsy[i],3] = 1

  ;derived variables
  resolve_derived_vars_deps

  ;further modifications
  IF KEYWORD_SET(doppler_freq) THEN BEGIN
    printerror, 'Info: Doppler correction is activated; omega_doppler = '+$
                rm0es(doppler_freq)
    print, 'initial ExB_stime = ', par.ExB_stime
    series.doppler_corr_mom = doppler_freq
  ENDIF ELSE series.doppler_corr_mom = 0.0

END

;##########################################################################

PRO calc_mom_ffts
; Calculates sxky, kxsy, sxsy from kxky in every time step for all
; requested fields, for all species.
; sxsy: standard via sxky, only if kxsy (not sxky) requested via kxsy
; NOTE: [0,0,0] on the lhs is used as a faster alternative to [*,*,*]

  COMMON global_vars

  FOR s = 0, gui.out.n_spec_sel - 1 DO BEGIN
     tmpptr = REFORM(mom[s,*])
     calc_ffts_3darr,tmpptr,series.request_mom
     mom[s,*] = tmpptr
  ENDFOR

END

;##########################################################################

FUNCTION get_mom_time, run, get_n_steps=get_n_steps, $
  first=first, last=last, step_time=step_time, coarse=coarse

  COMMON global_vars

  start_time = gui.out.start_t

  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
  file_field = get_file_path('field',run)
  file_mom = get_file_path('mom_'+spec[*].name,run)

  IF (file_field.exist + TOTAL(file_mom.exist)) EQ 0 THEN RETURN, 0.0

  request_field = TOTAL(series.request_mom[INDGEN(par.n_fields),*],2)
  request_mom = TOTAL(series.request_mom[INDGEN(par.n_moms)+par.n_fields,*],2)

  IF KEYWORD_SET(last) AND NOT KEYWORD_SET(coarse) THEN BEGIN
    istep_field = - (*par.istep_field)[run]
    istep_mom = - (*par.istep_mom)[run]
  ENDIF ELSE BEGIN ; for /last, use finer scale file, else use coarser
    istep_field = (*par.istep_field)[run]
    istep_mom = (*par.istep_mom)[run]
  ENDELSE

  IF file_field.exist AND TOTAL(request_mom) EQ 0 THEN file_choice = 0 $
    ELSE IF (TOTAL(file_mom[*].exist) NE 0) AND (TOTAL(file_field.exist) EQ 0) $
    THEN BEGIN

    sp = (WHERE(file_mom[*].exist))[0] ; first species with existing mom file
    file_choice = 1
  ENDIF ELSE IF file_field.exist AND (TOTAL(file_mom[*].exist) EQ 0) OR $
    (istep_field GE istep_mom) THEN $ ; if field file exists and is to be used
    file_choice = 0 $
  ELSE IF TOTAL(file_mom[*].exist) GT 0 THEN BEGIN
    sp = (WHERE(file_mom[*].exist))[0] ; first species with existing mom file
    file_choice = 1
  ENDIF

  IF par.write_h5 THEN BEGIN
    err = 1 - H5F_IS_HDF5(file_choice ? file_mom[sp].path : file_field.path)
    IF NOT err THEN tdata_id = $
    H5F_OPEN(file_choice ? file_mom[sp].path : file_field.path)
  ENDIF ELSE $
    OPENR, tdata_lun, file_choice ? file_mom[sp].path : file_field.path, $
    /GET_LUN, ERROR=err, /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian

  IF err NE 0 THEN BEGIN
    printerror, 'Please notify the developers about this error', /popup
    RETURN, 0L
  ENDIF

  IF (NOT KEYWORD_SET(last) AND (TOTAL(request_field) GT 0) AND $
    (TOTAL(request_mom) GT 0)) THEN BEGIN
    lcm_field_mom = lcm((*par.istep_mom)[run],(*par.istep_field)[run])
    interleave_factor = lcm_field_mom / $
      (file_choice ? (*par.istep_mom)[run] : (*par.istep_field)[run]) > 1L
  ENDIF ELSE interleave_factor = 1L

  time = par.prec_single ? 0.0 : 0.0D
  time_first_step = time
  time_second_step = time
  time_last_step = time

  IF par.write_h5 THEN BEGIN
    h5_dir_str = file_choice ? $
      'mom_' + spec[sp].name + '/time' : 'field/time'
    time_dset = H5D_OPEN(tdata_id,h5_dir_str)
    time_arr = H5D_READ(time_dset)
    time_arr *= time_renorm(0)

    n_step_attrib = H5A_OPEN_NAME(time_dset,'n_steps')
    n_steps = H5A_READ(n_step_attrib)
    H5A_CLOSE, n_step_attrib
    H5D_CLOSE, time_dset

    IF (SIZE(time_arr))[1] NE n_steps THEN BEGIN
      printerror, 'invalid time array for HDF5 get_mom_time'
      RETURN, 0L
    ENDIF

    IF file_choice THEN BEGIN
      series.h5_nst_mom = n_steps
      series.h5_nst_field = ROUND(n_steps*FLOAT((*par.istep_mom)[run])/FLOAT((*par.istep_field)[run]))
    ENDIF ELSE BEGIN
      series.h5_nst_field = n_steps
      series.h5_nst_mom = ROUND(n_steps*FLOAT((*par.istep_field)[run])/FLOAT((*par.istep_mom)[run]))
    ENDELSE

    IF KEYWORD_SET(get_n_steps) THEN BEGIN
      H5F_CLOSE, tdata_id
      RETURN, n_steps
    ENDIF

    time_first_step[0] = time_arr[0]
    time_last_step[0] = time_arr[n_steps-1]

    IF KEYWORD_SET(first) THEN BEGIN
      H5F_CLOSE, tdata_id
      RETURN, time_first_step
    ENDIF

    IF KEYWORD_SET(last) THEN BEGIN
      H5F_CLOSE, tdata_id
      RETURN, time_last_step
    ENDIF

    IF time_first_step GT start_time THEN BEGIN
      H5F_CLOSE, tdata_id
      RETURN, 0L
    ENDIF

    IF time_last_step LT start_time THEN BEGIN
      IF !QUIET NE 1 THEN PRINT, 'run ', run_label, ' not in time range'
      H5F_CLOSE, tdata_id
      RETURN, -2L
    ENDIF

    IF time_last_step EQ time_first_step THEN BEGIN
      PRINT, 'warning: run ', run_label, ' has only one time step'
      H5F_CLOSE, tdata_id
      RETURN, 0L
    ENDIF

    step_ind = (WHERE(FLOAT(time_arr) GE FLOAT(start_time)))[0] > 0L
    step_time = time_arr[step_ind]

    H5F_CLOSE, tdata_id

    RETURN, step_ind
  ENDIF ELSE BEGIN
    IF par.in_data NE 'sykx' THEN $
       content_byte_size = (8LL * (2 - par.prec_single) * par.nkx0 * (par.nky0 - $
          par.ky0_ind) * par.nz0 + 8) * (file_choice ? par.n_moms : par.n_fields) $
    ELSE content_byte_size = (8LL * (2 - par.prec_single) * (par.nkx0/2) * par.ny0 $
          * par.nz0 + 8) * (file_choice ? par.n_moms : par.n_fields)

    ; now data block + time entry
    full_entry_size = content_byte_size + 16 - 4 * par.prec_single

    step_byte_size = full_entry_size * interleave_factor

    n_steps = ROUND((FSTAT(tdata_lun)).SIZE/FLOAT(step_byte_size))
    pos_last_step = (n_steps - 1) * step_byte_size

    IF KEYWORD_SET(get_n_steps) THEN BEGIN
      FREE_LUN, tdata_lun
      RETURN, n_steps
    ENDIF

    READU, tdata_lun, time_first_step
    time_first_step *= time_renorm(0)

    IF KEYWORD_SET(first) THEN BEGIN
      FREE_LUN, tdata_lun
      RETURN, time_first_step
    ENDIF

    POINT_LUN, - tdata_lun, position_before
    position_after = position_before + content_byte_size
    POINT_LUN, tdata_lun, position_after

    IF NOT EOF(tdata_lun) THEN BEGIN
      READU, tdata_lun, time_second_step
      time_second_step *= time_renorm(0)
    ENDIF

    POINT_LUN, tdata_lun, pos_last_step
    READU, tdata_lun, time_last_step
    time_last_step *= time_renorm(0)
    step_time = time_last_step

    IF KEYWORD_SET(last) THEN BEGIN
      FREE_LUN, tdata_lun
      RETURN, time_last_step
    ENDIF

    avg_step_diff = (time_last_step - time_first_step) / FLOAT(n_steps)

    IF time_first_step GT start_time THEN BEGIN
      FREE_LUN, tdata_lun
      RETURN, 0L
    ENDIF

    IF time_last_step LT start_time THEN BEGIN
      IF !QUIET NE 1 THEN PRINT, 'run ', run_label, ' not in time range'
      FREE_LUN, tdata_lun
      RETURN, -2L
    ENDIF

    IF time_last_step EQ time_first_step THEN BEGIN
      PRINT, 'warning: run ', run_label, ' has only one time step'
      FREE_LUN, tdata_lun
      RETURN, 0L
    ENDIF

    IF time_second_step GT start_time THEN BEGIN
      FREE_LUN, tdata_lun
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

    WHILE (found EQ 0) AND (iterations LE max_iterations) AND $
      (estimated_pos GE 0) DO BEGIN

      IF estimated_pos GE 0 THEN $
        POINT_LUN, tdata_lun, estimated_pos ELSE BEGIN

        PRINT, 'estimated_pos < 0, file position error!'
        FREE_LUN, tdata_lun
        RETURN, -1L
      ENDELSE
      time_prev = par.prec_single ? 0.0 : 0.0D
      READU, tdata_lun, time_prev
      time_prev *= time_renorm(0)
      POINT_LUN, - tdata_lun, position_before
      position_after = position_before + content_byte_size + $
        full_entry_size * (interleave_factor - 1L)
      POINT_LUN, tdata_lun, position_after
      time_act = par.prec_single ? 0.0 : 0.0D
      READU, tdata_lun, time_act
      time_act *= time_renorm(0)
      time_diff_act = time_act - time_prev
      IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
;IF found THEN PRINT, estimated_pos
;IF found THEN PRINT, time_prev, start_time, time_act, ROUND((start_time-time_act)/time_diff_act)

      new_estimated_pos = step_byte_size * $
        ROUND((start_time-time_act)/time_diff_act) + estimated_pos

      IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
        new_estimated_pos = new_estimated_pos + step_byte_size

      estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)

      IF estimated_pos LT 0 THEN estimated_pos = 0L

      iterations += 1
    ENDWHILE
    IF found EQ 0 THEN BEGIN
      PRINT, 'using simplified search pattern'
      estimated_pos = estimated_pos_init
      WHILE found EQ 0 DO BEGIN
        IF estimated_pos GE 0 THEN POINT_LUN, tdata_lun, estimated_pos ELSE BEGIN
          PRINT, 'estimated_pos LT 0, file position error!'
          FREE_LUN, tdata_lun
          RETURN, -1L
        ENDELSE
        READU, tdata_lun, time_prev
        time_prev *= time_renorm(0)
        POINT_LUN, - tdata_lun, position_before
        position_after = position_before + content_byte_size + $
          full_entry_size*(interleave_factor-1L)
        POINT_LUN, tdata_lun, position_after
        READU, tdata_lun, time_act
        time_act *= time_renorm(0)
        time_diff_act = time_act - time_prev
        IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
        new_estimated_pos = step_byte_size * $
          ROUND((start_time-time_act)/time_diff_act) / ABS(ROUND((start_time-time_act)/$
          time_diff_act)) + estimated_pos
        IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
          new_estimated_pos = new_estimated_pos + step_byte_size
        estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)
        IF estimated_pos LT 0 THEN estimated_pos = 0L
      ENDWHILE
    ENDIF ; ELSE PRINT, 'number of iterations:', iterations

    step_pos = ROUND(estimated_pos/DOUBLE(step_byte_size))
    estimated_pos = step_pos * step_byte_size

    POINT_LUN, tdata_lun, estimated_pos
    step_time = par.prec_single ? 0.0 : 0.0D
    READU, tdata_lun, step_time
    step_time *= time_renorm(0)

    step_pos = estimated_pos / step_byte_size

    FREE_LUN, tdata_lun

    RETURN, step_pos
  ENDELSE

END

;##########################################################################

PRO jump_to_mom_step, run, jump_step, mom_luns, field_lun
; jumps to the desired step of mom and/or field file

  COMMON global_vars

  IF jump_step EQ -1 THEN BEGIN
    PRINT, 'omitting run, not in time range'
    RETURN
  ENDIF

  IF N_ELEMENTS(field_lun) LT 1 THEN field_lun = 0

  IF (field_lun GT 0) AND (TOTAL(mom_luns) GT 0) THEN BEGIN
    lcm_field_mom = lcm((*par.istep_mom)[run],(*par.istep_field)[run])
    interleave_factor_field = $
      lcm_field_mom / (*par.istep_field)[run] > 1L
    interleave_factor_mom = $
      lcm_field_mom / (*par.istep_mom)[run] > 1L
  ENDIF ELSE BEGIN
    interleave_factor_field = 1L
    interleave_factor_mom = 1L
  ENDELSE

  IF par.write_h5 THEN BEGIN
    series.h5_count_field = jump_step * interleave_factor_field
    series.h5_count_mom = jump_step * interleave_factor_mom
    RETURN
  ENDIF

  ; get step byte size, jump to position
  field_step_bytes = 0
  FOR sp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    IF TOTAL(series.request_mom[INDGEN(par.n_moms)+par.n_fields,*]) GT 0 THEN BEGIN
      time = par.prec_single ? 0.0 : 0.0D
      IF par.in_data NE 'sykx' THEN content = par.prec_single ? $
        COMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) : $
        DCOMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) $
        ELSE content = par.prec_single ? $
        COMPLEXARR(par.nkx0/2,par.ny0,par.nz0,/NOZERO) : $
        DCOMPLEXARR(par.nkx0/2,par.ny0,par.nz0,/NOZERO)
      READU, mom_luns[sp], time
      FOR mom_nr = 0, par.n_moms - 1 DO READU, mom_luns[sp], content
      POINT_LUN, - mom_luns[sp], mom_step_bytes
      jump_byte = LONG64(mom_step_bytes) * jump_step * interleave_factor_mom
      POINT_LUN, mom_luns[sp], jump_byte
    ENDIF
  ENDFOR

  IF TOTAL(series.request_mom[INDGEN(par.n_fields),*]) GT 0 THEN BEGIN
    time = par.prec_single ? 0.0 : 0.0D
    IF par.in_data NE 'sykx' THEN content = par.prec_single ? $
      COMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) : $
      DCOMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) $
      ELSE content = par.prec_single ? $
      COMPLEXARR(par.nkx0/2,par.ny0,par.nz0,/NOZERO) : $
      DCOMPLEXARR(par.nkx0/2,par.ny0,par.nz0,/NOZERO)
    READU, field_lun, time
    FOR field_nr = 0, par.n_fields - 1 DO READU, field_lun, content
    POINT_LUN, - field_lun, field_step_bytes
    jump_byte = LONG64(field_step_bytes) * jump_step * interleave_factor_field
    POINT_LUN, field_lun, jump_byte
  ENDIF

END

;##########################################################################

PRO read_mom_step, field_lun, mom_luns

  COMMON global_vars

  mom_time = par.prec_single ? 0.0 : 0.0D

  request_field = TOTAL(series.request_mom[INDGEN(par.n_fields),*],2)
  request_mom = TOTAL(series.request_mom[INDGEN(par.n_moms)+par.n_fields,*],2)

  IF par.in_data NE 'sykx' THEN $
    content_byte_size = 8LL * (2 - par.prec_single) * par.nkx0 * $
    (par.nky0 - par.ky0_ind) * par.nz0 + 8 ELSE $
    content_byte_size = 8LL * (2 - par.prec_single) * (par.nkx0/2) * $
    par.ny0 * par.nz0 + 8

  IF par.write_h5 THEN BEGIN ; mapping field name to number
    field_map = ['phi','A_par','B_par']
    mom_map = par.n_moms LE 6 ? $
      ['dens','T_par','T_perp','q_par','q_perp','u_par'] : $
      ['dens_pass','T_par_pass','T_perp_pass','q_par_pass',$
      'q_perp_pass','u_par_pass','dens_trap','T_par_trap',$
      'T_perp_trap','q_par_trap','q_perp_trap','u_par_trap',$
      'dens_FLR','T_par_FLR','T_perp_FLR','q_par_FLR',$
      'q_perp_FLR','u_par_FLR']
  ENDIF

  IF field_lun GT 0 THEN BEGIN
    eof_field = par.write_h5 ? $
      (series.h5_nst_field LE series.h5_count_field) : EOF(field_lun)
    IF eof_field THEN BEGIN
      IF par.write_h5 THEN H5F_CLOSE, field_lun ELSE FREE_LUN, field_lun
      field_lun = 0
    ENDIF ELSE BEGIN
      time = par.prec_single ? 0.0 : 0.0D

      IF par.in_data NE 'sykx' THEN BEGIN
        content = par.prec_single ? $
          COMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) : $
          DCOMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) 
      ENDIF ELSE BEGIN
        content = par.prec_single ? $
          COMPLEXARR(par.ny0,par.nkx0/2,par.nz0,/NOZERO) : $
          DCOMPLEXARR(par.ny0,par.nkx0/2,par.nz0,/NOZERO)            
      ENDELSE


      IF par.write_h5 THEN BEGIN
        ; take time from first field
        step_nr = STRING(series.h5_count_field,FORMAT='(I010)')
        step_loc = '/field/' + field_map[0] + '/' + step_nr + '/'
        step_dset = H5D_OPEN(field_lun,step_loc)
        time_attrib = H5A_OPEN_NAME(step_dset,'time')
        time = H5A_READ(time_attrib)
        H5A_CLOSE, time_attrib
        H5D_CLOSE, step_dset
      ENDIF ELSE READU, field_lun, time

      time *= time_renorm(0)

      FOR field_nr = 0, par.n_fields - 1 DO BEGIN
        IF request_field[field_nr] GE 1 THEN BEGIN
          IF !QUIET NE 1 THEN BEGIN
            PRINT, '    reading ' + rm0es(content_byte_size/1048576.0) + $
	      ' MB from field file (field time: ' + rm0es(time) + ')'
	    stime = SYSTIME(1)
          ENDIF

          IF par.write_h5 THEN BEGIN
            step_nr = STRING(series.h5_count_field,FORMAT='(I010)')
            step_loc = '/field/' + field_map[field_nr] + '/' + step_nr + '/'
            field_dset = H5D_OPEN(field_lun,step_loc)

            content_reim = H5D_READ(field_dset)
            H5D_CLOSE, field_dset

            content[0,0,0] = COMPLEX(content_reim.real,$
              content_reim.imaginary,DOUBLE=(1-par.prec_single))
            content_reim = 0
          ENDIF ELSE READU, field_lun, content

          IF !QUIET NE 1 THEN $
            PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'

          IF cfg.kpar_filter GE 0 THEN BEGIN
            IF (cfg.kpar_filter LT par.nz0 / 2) AND (par.nz0 GT 1) THEN BEGIN
              content_kpar = FFT(content,DIMENSION=3)
              IF cfg.kpar_filter EQ 0 THEN BEGIN
                content_kpar[*,*,1:*] = 0
              ENDIF ELSE BEGIN
                neg_kpar = par.nz0 - cfg.kpar_filter
                content_kpar[*,*,0:cfg.kpar_filter-1] = 0
                content_kpar[*,*,cfg.kpar_filter+1:neg_kpar-1] = 0
                IF cfg.kpar_filter NE 1 THEN content_kpar[*,*,neg_kpar+1:*] = 0
              ENDELSE
              content[0,0,0] = FFT(content_kpar,DIMENSION=3,/INVERSE)
            ENDIF
          ENDIF

          FOR j = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO BEGIN
	    IF par.in_data EQ 'sxky' THEN BEGIN
    	      IF par.ky0_ind GT 0 THEN $
                (*mom[j,field_nr].sxky)[*,0:par.ky0_ind-1,*] = 0
	      (*mom[j,field_nr].sxky)[0,par.ky0_ind,0] = $
                TEMPORARY(content*mom_renorm(field_nr))
	    ENDIF ELSE IF par.in_data EQ 'sykx' THEN BEGIN
               FOR ik=0,par.nz0-1 DO BEGIN
                  (*mom[j,field_nr].kxsy)[0:par.nkx0/2-1,*,ik] = $
                     TRANSPOSE(content[*,*,ik]*$
                     mom_renorm(field_nr))
                  (*mom[j,field_nr].kxsy)[par.nkx0/2:*,*,ik] = 0.0
               ENDFOR
            ENDIF ELSE BEGIN
    	      IF par.ky0_ind GT 0 THEN BEGIN
                (*mom[j,field_nr].kxky)[*,0:par.ky0_ind-1,*] = 0
                (*mom[j,field_nr].kxky)[0,par.ky0_ind,0] = $
                  TEMPORARY(content*mom_renorm(field_nr))
              ENDIF ELSE (*mom[j,field_nr].kxky) = $
                TEMPORARY(content*mom_renorm(field_nr))
	    ENDELSE
          ENDFOR
        ENDIF ELSE IF NOT par.write_h5 THEN BEGIN ; if data not requested, jump to next field
          POINT_LUN, - field_lun, position_before
          position_after = position_before + content_byte_size
          POINT_LUN, field_lun, position_after
        ENDIF
      ENDFOR

      mom_time = time
    ENDELSE
  ENDIF

  FOR j = 0, gui.out.n_spec_sel - 1 DO BEGIN
    IF mom_luns[j] GT 0 THEN BEGIN
      eof_mom = par.write_h5 ? $
        (series.h5_nst_mom LE series.h5_count_mom) : EOF(mom_luns[j])
      IF eof_mom THEN BEGIN
        IF par.write_h5 THEN H5F_CLOSE, mom_luns[j] ELSE FREE_LUN, mom_luns[j]
	mom_luns[j] = 0
      ENDIF ELSE BEGIN
        time = par.prec_single ? 0.0 : 0.0D

        IF par.in_data NE 'sykx' THEN BEGIN
          content = par.prec_single ? $
            COMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) : $
            DCOMPLEXARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) 
        ENDIF ELSE BEGIN
          content = par.prec_single ? $
            COMPLEXARR(par.ny0,par.nkx0/2,par.nz0,/NOZERO) : $
            DCOMPLEXARR(par.ny0,par.nkx0/2,par.nz0,/NOZERO)            
        ENDELSE

        IF par.write_h5 THEN BEGIN
          ; take time from first moment
          step_nr = STRING(series.h5_count_mom,FORMAT='(I010)')
          sp = (*gui.out.spec_select)[j]
          step_loc = '/mom_' + spec[sp].name + '/' + mom_map[0] + '/' + step_nr + '/'
          step_dset = H5D_OPEN(mom_luns[j],step_loc)
          time_attrib = H5A_OPEN_NAME(step_dset,'time')
          time = H5A_READ(time_attrib)
          H5A_CLOSE, time_attrib
          H5D_CLOSE, step_dset
        ENDIF ELSE READU, mom_luns[j], time

        time *= time_renorm(0)

        FOR mom_nr = 0, par.n_moms - 1 DO BEGIN
          IF request_mom[mom_nr] GE 1 THEN BEGIN
            IF !QUIET NE 1 THEN BEGIN
              PRINT, '    reading ' + rm0es(content_byte_size/1048576.0) + $
                ' MB from mom file (mom time: ' + rm0es(time) + ')'
              stime = SYSTIME(1)
	    ENDIF

            IF par.write_h5 THEN BEGIN
              step_nr = STRING(series.h5_count_mom,FORMAT='(I010)')
              step_loc = '/mom_' + spec[sp].name + '/' + mom_map[mom_nr] + '/' + step_nr + '/'
              mom_dset = H5D_OPEN(mom_luns[j],step_loc)

              content_reim = H5D_READ(mom_dset)
              H5D_CLOSE, mom_dset

              content[0,0,0] = COMPLEX(content_reim.real,$
                content_reim.imaginary,DOUBLE=(1-par.prec_single))
              content_reim = 0
            ENDIF ELSE READU, mom_luns[j], content

            IF cfg.kpar_filter GE 0 THEN BEGIN
              IF (cfg.kpar_filter LT par.nz0 / 2) AND (par.nz0 GT 1) THEN BEGIN
                content_kpar = FFT(content,DIMENSION=3)
                IF cfg.kpar_filter EQ 0 THEN BEGIN
                  content_kpar[*,*,1:*] = 0
                ENDIF ELSE BEGIN
                  neg_kpar = par.nz0 - cfg.kpar_filter
                  content_kpar[*,*,0:cfg.kpar_filter-1] = 0
                  content_kpar[*,*,cfg.kpar_filter+1:neg_kpar-1] = 0
                  IF cfg.kpar_filter NE 1 THEN content_kpar[*,*,neg_kpar+1:*] = 0
                ENDELSE
                content[0,0,0] = FFT(content_kpar,DIMENSION=3,/INVERSE)
              ENDIF
            ENDIF

            IF par.in_data EQ 'sxky' THEN BEGIN
	      IF par.ky0_ind GT 0 THEN $
                (*mom[j,mom_nr+par.n_fields].sxky)[*,0:par.ky0_ind-1,*] = 0
	      (*mom[j,mom_nr+par.n_fields].sxky)[0,par.ky0_ind,0] = $
	      	TEMPORARY(content*mom_renorm(mom_nr+par.n_fields))
	    ENDIF ELSE IF par.in_data EQ 'sykx' THEN BEGIN
              FOR k = 0, par.nz0 - 1 DO BEGIN
                (*mom[j,mom_nr+par.n_fields].kxsy)[0:par.nkx0/2-1,*,k] = $
                  TRANSPOSE(content[*,*,k]*$
                  mom_renorm(mom_nr+par.n_fields))
                (*mom[j,mom_nr+par.n_fields].kxsy)[par.nkx0/2:*,*,k] = 0.0
              ENDFOR
            ENDIF ELSE BEGIN
              IF par.ky0_ind GT 0 THEN BEGIN
                (*mom[j,mom_nr+par.n_fields].kxky)[*,0:par.ky0_ind-1,*] = 0
                (*mom[j,mom_nr+par.n_fields].kxky)[0,par.ky0_ind,0] = $
                  TEMPORARY(content*mom_renorm(mom_nr+par.n_fields))
              ENDIF ELSE *mom[j,mom_nr+par.n_fields].kxky = $
                TEMPORARY(content*mom_renorm(mom_nr+par.n_fields))
	    ENDELSE

            IF !QUIET NE 1 THEN PRINT, '    elapsed time: ' + $
              rm0es(SYSTIME(1)-stime) + ' s'
          ENDIF ELSE IF NOT par.write_h5 THEN BEGIN ; if data not requested, jump to next field
            POINT_LUN, - mom_luns[j], position_before
            position_after = position_before + content_byte_size
            POINT_LUN, mom_luns[j], position_after
          ENDIF
        ENDFOR

        mom_compute_derived_vars, j

        IF mom_time NE 0.0 THEN BEGIN
           IF !QUIET NE 1 THEN PRINT, 'mom time = ',rm0es(time), $
             ', field time = ', rm0es(field_time)

           IF mom_time NE time THEN BEGIN
              IF (j EQ 0) THEN $
                 PRINT, 'warning: field and mom files have asynchronous time stamps' $
              ELSE PRINT, 'warning: species have asynchronous mom time stamps'
           ENDIF
        ENDIF ELSE mom_time = time
     ENDELSE ; eof
   ENDIF ; mom_luns[j] GT 0? 
 ENDFOR ; j loop (species)

END

;##########################################################################

PRO mom_loop, diag0

  COMMON global_vars

  IF TOTAL(series.request_mom) EQ 0 THEN BEGIN
    printerror, 'no mom or field data requested in mom/field loop'
    RETURN
  ENDIF

  last_time = -1.0D * 1e20

  file_field = get_file_path('field')
  file_mom = get_file_path('mom_'+spec[(*gui.out.spec_select)[*]].name)

  request_field = TOTAL(series.request_mom[INDGEN(par.n_fields),*])
  request_mom = TOTAL(series.request_mom[INDGEN(par.n_moms)+par.n_fields,*])

  IF request_field AND (file_field.exist EQ 0) THEN BEGIN
    printerror, 'Skipping mom loop due to missing field file(s)'
    RETURN
  ENDIF

  IF request_mom AND (TOTAL(file_mom.exist) NE gui.out.n_spec_sel) THEN BEGIN
    printerror, 'Skipping mom loop due to missing mom file(s)'
    RETURN
  ENDIF

  IF gui.out.sparsefac GT 1 THEN PRINT, 'Sparse mom loop: reading only every ' + $
    append_th(gui.out.sparsefac) + ' step'

  mom_luns = LONARR(par.n_spec)
  field_lun = 0
  create_mom_struct

  IF par.in_data NE 'sykx' THEN BEGIN
     field_content_byte_size = (8LL * (2 - par.prec_single) * par.nkx0 * $
        (par.nky0 - par.ky0_ind) * par.nz0 + 8) * par.n_fields + $
         16 - 4 * par.prec_single
     mom_content_byte_size = (8LL * (2 - par.prec_single) * par.nkx0 * $
        (par.nky0 - par.ky0_ind) * par.nz0 + 8) * par.n_moms + $
         16 - 4 * par.prec_single
  ENDIF ELSE BEGIN
     field_content_byte_size = (8LL * (2 - par.prec_single) * (par.nkx0/2) * $
        par.ny0 * par.nz0 + 8) * par.n_fields + $
        16 - 4 * par.prec_single
     mom_content_byte_size = (8LL * (2 - par.prec_single) * (par.nkx0/2) * $
        par.ny0 * par.nz0 + 8) * par.n_moms + $
        16 - 4 * par.prec_single     
  ENDELSE

  manual_abort = 0

  FOR run = 0, series.n_runs - 1 DO BEGIN
    run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]

    interleave_field_factor = LONG64(gui.out.sparsefac)
    interleave_mom_factor = LONG64(gui.out.sparsefac)

    IF ((*par.istep_field)[run] GT 0) AND ((*par.istep_mom)[run] GT 0) AND $
      (request_field GT 0) AND (request_mom GT 0) THEN BEGIN

      lcm_field_mom = lcm((*par.istep_mom)[run],(*par.istep_field)[run])
      interleave_field = lcm_field_mom / (*par.istep_field)[run]
      interleave_mom = lcm_field_mom / (*par.istep_mom)[run]
      interleave_field_factor *= LONG64(interleave_field)
      interleave_mom_factor *= LONG64(interleave_mom)

      IF (interleave_field NE 1) THEN printerror, $
        'using just every '+append_th(interleave_field)+$
        ' field step to synchronize with mom file'
      IF (interleave_mom NE 1) THEN printerror, $
        'using just every '+append_th(interleave_mom)+$
        ' mom step to synchronize with field file'
      IF ((interleave_field NE 1) AND (interleave_mom NE 1)) THEN $
        printerror, 'Interleaving mom AND field files! Try to use consistent '+$
        'istep_mom/field settings'
    ENDIF

    IF request_field GT 0 THEN BEGIN
      IF par.write_h5 THEN BEGIN
        err = 1 - H5F_IS_HDF5(file_field.path[run])
        IF NOT err THEN field_lun = H5F_OPEN(file_field.path[run])

        series.h5_count_field = 0
      ENDIF ELSE OPENR, field_lun, file_field.path[run], /GET_LUN, $
        ERROR=err, /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian

      IF err NE 0 THEN BEGIN
        printerror, file_field.path[run] + ' does not exist', /popup
        RETURN
      ENDIF ELSE IF !QUIET NE 1 THEN PRINT, 'reading ', file_field.path[run]
    ENDIF

    IF request_mom GT 0 THEN BEGIN
      FOR j = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO BEGIN
        IF par.write_h5 THEN BEGIN
          err = 1 - H5F_IS_HDF5(file_mom[j].path[run])
          IF NOT err THEN temp_lun = H5F_OPEN(file_mom[j].path[run])

          IF j EQ 0 THEN series.h5_count_mom = 0
        ENDIF ELSE OPENR, temp_lun, file_mom[j].path[run], /GET_LUN, $
          ERROR=err, /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian

        IF err NE 0 THEN BEGIN
          printerror, file_mom[j].path[run] + ' does not exist', /popup
          RETURN
        ENDIF ELSE BEGIN
          IF !QUIET NE 1 THEN PRINT, 'reading ', file_mom[j].path[run]
          mom_luns[j] = temp_lun
        ENDELSE
      ENDFOR
    ENDIF

    IF (series.doppler_corr_mom NE 0.0) THEN BEGIN
       run_start_time = get_mom_time(run,/first) 
       IF (run_start_time GE par.ExB_stime) THEN BEGIN
          print, 'adjusting ExB_stime to '+rm0es(run_start_time)
          par.ExB_stime = run_start_time
       ENDIF
    ENDIF

    jump_step = get_mom_time(run,step_time=step_time) ; determine start step 
    IF (jump_step LT 0) OR manual_abort THEN BEGIN
      end_loop = 1
      IF !QUIET NE 1 THEN PRINT, 'omitting run ', run_label
    ENDIF ELSE BEGIN
      jump_to_mom_step, run, jump_step, mom_luns, field_lun

      IF jump_step GT 0 THEN $
        PRINT, '- starting field/mom loop at step ', rm0es(jump_step), ' -'
      IF jump_step EQ 0 THEN PRINT, '- starting field/mom loop -'
      IF ((TOTAL(file_mom.exist) + file_field.exist) GE 1) AND $
        (!QUIET EQ 1) THEN PRINT, 'reading run ', run_label, ': field/mom' + $
        (par.write_h5 ? ' (HDF5)' : '')
      end_loop = 0
    ENDELSE

    ; for unstable runs which have follow-ups from secure checkpoints
    ; (overlapping time windows):
    IF run LT series.n_runs - 1 THEN $
      start_next = get_mom_time(run+1,/first) ELSE $
      start_next = 2.0 * ABS(gui.out.end_t)

    IF NOT par.x_local THEN read_profiles, run=run

    ; --- start mom loop ---

    printed_time_to_pwin = 0
    WHILE NOT end_loop DO BEGIN
      read_mom_step, field_lun, mom_luns

      IF par.write_h5 THEN BEGIN
        IF field_lun GT 0 THEN $
          series.h5_count_field += interleave_field_factor
        IF TOTAL(mom_luns) GT 0 THEN $
          series.h5_count_mom += interleave_mom_factor
      ENDIF ELSE BEGIN ; set file pointers to next entry
        IF field_lun GT 0 THEN BEGIN
          POINT_LUN, - field_lun, position_before
          position_after = position_before + field_content_byte_size * $
            LONG64(interleave_field_factor-1)
          POINT_LUN, field_lun, position_after
        ENDIF

        FOR j = 0, gui.out.n_spec_sel - 1 DO BEGIN
          IF mom_luns[j] GT 0 THEN BEGIN
            POINT_LUN, - mom_luns[j], position_before
            position_after = position_before + mom_content_byte_size * $
              LONG64(interleave_mom_factor-1)
            POINT_LUN, mom_luns[j], position_after
          ENDIF
        ENDFOR
      ENDELSE

      IF last_time LT mom_time THEN BEGIN
        IF mom_time GT (gui.out.end_t * (1.0D + 1e-7)) THEN end_loop = 1 $
          ELSE IF (mom_time GT start_next) THEN BEGIN

          PRINT, 'overlapping time windows; continuing with follow-up'
          end_loop = 1
        ENDIF ELSE IF mom_time GE gui.out.start_t * (1.0D - 1e-7) THEN BEGIN
          series.step_count += 1

          ; erase previous time in plot window, print new time
          XYOUTS, 0.75, 0.95, '!6t = ' + rm0es(last_time), COLOR=0, /NORMAL, $
            CHARSIZE=1.3
          XYOUTS, 0.75, 0.95, '!6t = ' + rm0es(mom_time), COLOR=1, /NORMAL, $
            CHARSIZE=1.3
          printed_time_to_pwin = 1

          PRINT, '  time: ', rm0es(mom_time)
          last_time = mom_time

          calc_mom_ffts
          call_diags, diag0, 'loop'

          ; break loop by clicking and holding left in plot window
          IF gui.window.main NE 0 THEN BEGIN ; omit for command line version
            CURSOR, xtemp, ytemp, /NOWAIT
            IF !MOUSE.button EQ 1 THEN BEGIN
              PRINT, 'loop aborted by user'
              end_loop = 1
              manual_abort = 1
            ENDIF  ; mouse click
          ENDIF ; gui.window.main
        ENDIF ;mom_time GT start?
      ENDIF ; last time LT data time

      end_loop = end_loop OR (request_field AND (field_lun EQ 0)) $
        OR (request_mom AND (TOTAL(mom_luns) EQ 0))
    ENDWHILE ; current field/mom loop

    IF printed_time_to_pwin THEN XYOUTS, 0.75, 0.95, $
      '!6t = ' + rm0es(last_time), COLOR=0, /NORMAL, CHARSIZE=1.3

    ; close current mom and field files
    IF request_mom GT 0 THEN FOR j = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO $
      IF mom_luns[j] GT 0 THEN IF par.write_h5 THEN $
      H5F_CLOSE, mom_luns[j] ELSE FREE_LUN, mom_luns[j]
    IF request_field GT 0 AND field_lun GT 0 THEN IF par.write_h5 THEN $
      H5F_CLOSE, field_lun ELSE FREE_LUN, field_lun
  ENDFOR

  destroy_mom_struct

  IF jump_step GE 0 THEN PRINT, '- finished field/mom loop -'

END

;##########################################################################

PRO resolve_derived_vars_deps
; Resolve dependencies of derived variables, ie. request those
; mom entries which are required for construction

  COMMON global_vars
  
  nf = par.n_fields
 
  ; check quasineutrality (over all selected species)
  IF TOTAL(series.request_mom[199,*]) GT 0 THEN $
    series.request_mom[nf,par.in_data EQ 'sxky'] = 1   

  ; total density gradient
  IF TOTAL(series.request_mom[104,*]) GT 0 THEN $
    series.request_mom[nf,par.in_data EQ 'sxky'] = 1 

  ; total temperature gradient
  IF TOTAL(series.request_mom[103,*]) GT 0 THEN $
    series.request_mom[100,par.in_data EQ 'sxky'] = 1 

  ; shearing rate
  IF TOTAL(series.request_mom[102,*]) GT 0 THEN BEGIN
    series.request_mom[0,par.in_data EQ 'sxky'] = 1
    series.request_mom[nf,par.in_data EQ 'sxky'] = 1
    ; dummy read for now to have mom_compute_derived_vars executed
  ENDIF
     
  ; (partial, i.e., per species) pressure 
  IF TOTAL(series.request_mom[101,*]) GT 0 THEN BEGIN
    series.request_mom[nf,par.in_data EQ 'sxky'] = 1
    series.request_mom[100,par.in_data EQ 'sxky'] = 1
  ENDIF

  ; total temperature
  IF TOTAL(series.request_mom[100,*]) GT 0 THEN $
    series.request_mom[nf+1:nf+2,par.in_data EQ 'sxky'] = 1

  ; --- transport quantities ---

  ; Gamma_es = - n i ky Phi
  IF TOTAL(series.request_mom[150,*]) GT 0 THEN BEGIN
    series.request_mom[0,par.in_data EQ 'sxky'] = 1
    series.request_mom[nf,par.in_data EQ 'sxky'] = 1
  ENDIF

  ; Q_es = - (1/2 Tpar + Tperp + 3/2 n) i ky Phi
  IF TOTAL(series.request_mom[151,*]) GT 0 THEN BEGIN
    series.request_mom[0,par.in_data EQ 'sxky'] = 1
    series.request_mom[nf:nf+2,par.in_data EQ 'sxky'] = 1
  ENDIF

  ; P_es = - upar i ky Phi
  IF TOTAL(series.request_mom[152,*]) GT 0 THEN BEGIN
    series.request_mom[0,par.in_data EQ 'sxky'] = 1
    series.request_mom[nf+5,par.in_data EQ 'sxky'] = 1
  ENDIF

  ; Gamma_em = upar i ky Apar + densI1 i ky Bpar
  IF TOTAL(series.request_mom[153,*]) GT 0 THEN BEGIN
    IF nf GE 1 THEN BEGIN
      series.request_mom[1,par.in_data EQ 'sxky'] = 1
      series.request_mom[nf+5,par.in_data EQ 'sxky'] = 1

      IF (nf GE 2) AND (par.n_moms GT 6) THEN BEGIN ; add Bpar term only if available
        series.request_mom[2,par.in_data EQ 'sxky'] = 1
        series.request_mom[nf+6,par.in_data EQ 'sxky'] = 1
      ENDIF
    ENDIF ELSE series.request_mom[153,*] = 0
  ENDIF

  ; Q_em = (qpar + qperp + 5/2 upar) i ky Apar - (TparI1 + TperpI1) i ky Bpar
  IF TOTAL(series.request_mom[154,*]) GT 0 THEN BEGIN
    IF nf GE 1 THEN BEGIN
      series.request_mom[1,par.in_data EQ 'sxky'] = 1
      series.request_mom[nf+3:nf+4,par.in_data EQ 'sxky'] = 1

      IF (nf GE 2) AND (par.n_moms GT 6) THEN BEGIN ; add Bpar term only if available
        series.request_mom[2,par.in_data EQ 'sxky'] = 1
        series.request_mom[nf+7:nf+8,par.in_data EQ 'sxky'] = 1
      ENDIF
    ENDIF ELSE series.request_mom[154,*] = 0
  ENDIF

  ; P_em = (Tpar + n) i ky Apar
  IF TOTAL(series.request_mom[155,*]) GT 0 THEN BEGIN
    IF nf GE 1 THEN BEGIN
      series.request_mom[1,par.in_data EQ 'sxky'] = 1
      series.request_mom[nf+1,par.in_data EQ 'sxky'] = 1
    ENDIF ELSE series.request_mom[155,*] = 0
  ENDIF

END

;##########################################################################

PRO mom_compute_derived_vars, j
; Compute derived variables, such as B_x, B_y

  COMMON global_vars

  nf = par.n_fields
  sp = (*gui.out.spec_select)[j]

  in_data = par.in_data EQ 'sykx' ? 'kxsy' : par.in_data
  momtag = get_tag_index(mom[0,0],in_data)

  FOR mom_nr = 100, N_ELEMENTS(series.request_mom[*,0]) - 1 DO BEGIN
    IF TOTAL(series.request_mom[mom_nr,*]) GE 1 THEN BEGIN
      CASE mom_nr OF
        100 : *mom[j,mom_nr].(momtag) = $ ; total temperature
                ((*mom[j,nf+1].(momtag)) + 2.0 * (*mom[j,nf+2].(momtag))) / 3.0
        101 : BEGIN ; (partial) pressure = n_1 T_0 + n_0 T_1
                press_sp = spec[sp].temp*spec[sp].dens
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'profiles need to be implemented in pressure computation'
                *mom[j,mom_nr].(momtag) = ((*mom[j,nf].(momtag)) + (*mom[j,100].(momtag))) * press_sp
              END
        102 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'shearing rate not implemented for global code yet'
                *mom[j,mom_nr].(momtag) = - REBIN(REFORM((*series.kx)^2,$
                  [par.nkx0,par.nky0,1]),[par.nkx0,par.nky0,par.nz0]) * (*mom[j,0].(momtag))
              END
        103 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'radial temperature gradient not implemented for global code yet'
                *mom[j,mom_nr].(momtag) = - COMPLEX(0,1) * $
                  REBIN(REFORM((*series.kx),[par.nkx0,par.nky0,1]),$
                  [par.nkx0,par.nky0,par.nz0]) * (*mom[j,100].(momtag))
              END
        104 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'radial density gradient not implemented for global code yet'
                *mom[j,mom_nr].(momtag) = - COMPLEX(0,1) * $
                  REBIN(REFORM((*series.kx),[par.nkx0,par.nky0,1]),$
                  [par.nkx0,par.nky0,par.nz0]) * (*mom[j,nf].(momtag))
              END
        150 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'Gamma_es not yet implemented for global code'
                *mom[j,mom_nr].(momtag) = - COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,0].(momtag)) * CONJ(*mom[j,nf].(momtag)) * $
                  (spec[j].dens / series.Bref) / REBIN(REFORM(SQRT((*series.geom).gxx),[1,1,par.nz0]),$
                  [par.nkx0,par.nky0,par.nz0])
              END
        151 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'Q_es not yet implemented for global code'
                *mom[j,mom_nr].(momtag) = - COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,0].(momtag)) * CONJ((0.5*(*mom[j,nf+1].(momtag))*series.nref+$
                  (*mom[j,nf+2].(momtag))*series.nref+1.5*(*mom[j,nf].(momtag))*series.Tref)) * $
                  (spec[j].dens * spec[j].temp * series.Qref / series.Bref) / $
                  REBIN(REFORM(SQRT((*series.geom).gxx),[1,1,par.nz0]),$
                  [par.nkx0,par.nky0,par.nz0])
              END
        152 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'P_es not yet implemented for global code'
                *mom[j,mom_nr].(momtag) = - COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,0].(momtag)) * CONJ(*mom[j,nf+5].(momtag)) * $
                  (spec[j].dens * spec[j].mass / series.Bref) / $
                  REBIN(REFORM(SQRT((*series.geom).gxx),[1,1,par.nz0]),$
                  [par.nkx0,par.nky0,par.nz0])
              END
        153 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'Gamma_em not yet implemented for global code'
                *mom[j,mom_nr].(momtag) = COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,1].(momtag)) * CONJ(*mom[j,nf].(momtag)) * $
                  (spec[j].dens / series.Bref)
                IF (nf GE 2) AND (par.n_moms GT 6) THEN BEGIN
                  *mom[j,mom_nr].(momtag) -= COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,2].(momtag)) * CONJ(*mom[j,nf+6].(momtag)) * $
                  (spec[j].dens * spec[j].temp / (spec[j].charge * series.Bref)) * $
                  REBIN(REFORM(1.0/(*series.geom).Bfield,[1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0])
                ENDIF
                *mom[j,mom_nr].(momtag) /= REBIN(REFORM(SQRT((*series.geom).gxx),$
                  [1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0])
              END
        154 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'Q_em not yet implemented for global code'
                *mom[j,mom_nr].(momtag) = COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,1].(momtag)) * CONJ((*mom[j,nf+3].(momtag))+(*mom[j,nf+4].(momtag))) * $
                  (spec[j].dens * spec[j].temp / series.Bref)
                IF (nf GE 2) AND (par.n_moms GT 6) THEN BEGIN
                  *mom[j,mom_nr].(momtag) -= COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,2].(momtag)) * CONJ((*mom[j,nf+7].(momtag))+(*mom[j,nf+8].(momtag))) * $
                  (spec[j].dens * spec[j].temp^2 / (spec[j].charge * series.Bref)) * $
                  REBIN(REFORM(1.0/(*series.geom).Bfield,[1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0])
                ENDIF
                *mom[j,mom_nr].(momtag) /= REBIN(REFORM(SQRT((*series.geom).gxx),$
                  [1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0])
              END
        155 : BEGIN
                IF par.in_data EQ 'sxky' THEN $
                  printerror, 'P_em not yet implemented for global code'
                *mom[j,mom_nr].(momtag) = COMPLEX(0,1) * $
                  REBIN(REFORM(*series.ky,[1,par.nky0]),[par.nkx0,par.nky0,par.nz0]) * $
                  (*mom[j,1].(momtag)) * CONJ((*mom[j,nf].(momtag))*spec[j].temp+$
                  (*mom[j,nf+1].(momtag))*spec[j].dens) * $
                  (spec[j].dens * spec[j].temp / series.Bref) / $
                  REBIN(REFORM(SQRT((*series.geom).gxx),[1,1,par.nz0]),$
                  [par.nkx0,par.nky0,par.nz0])
              END
        199 : BEGIN ; quasineutrality (over selected species)
                *mom[j,mom_nr].(momtag) = (*mom[j,nf].(momtag)) * spec[sp].charge * spec[sp].dens
                ;sum if last species
                IF j EQ gui.out.n_spec_sel - 1 THEN BEGIN
                  sumspec = *mom[0,mom_nr].(momtag)
                  FOR n = 1, gui.out.n_spec_sel - 1 DO sumspec += *mom[n,mom_nr].(momtag)
                  FOR n = 0, gui.out.n_spec_sel - 1 DO *mom[n,mom_nr].(momtag) = sumspec
                ENDIF
              END
        ELSE : printerror, 'error in mom_compute_derived_vars: var ' + $
                 rm0es(mom_nr) + ' not found', /popup
      ENDCASE
    ENDIF
  ENDFOR

END
