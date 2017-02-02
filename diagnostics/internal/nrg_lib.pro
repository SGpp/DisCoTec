;##########################################################################
;# this library contains almost all nrg related internal functions        #
;##########################################################################

FUNCTION nrg_renorm, var, time=time

  COMMON global_vars

  Lref = series.Lref
  mref = series.mref
  Tref = series.Tref ; in eV
  nref = series.nref
  Qref = series.Qref
  Bref = series.Bref
  v_ref = SQRT(Tref*Qref/mref)

  k_B = 1.0D

  ; T is given in eV, the elementary charge will get into it via Qref

  IF KEYWORD_SET(time) THEN RETURN, Lref / v_ref

  rho_ref = SQRT(mref*Tref/Qref) / Bref
  D_gb = v_ref * rho_ref^2 / Lref

  CASE (var MOD 10) OF ; same for each species
    0 : RETURN, (nref*rho_ref/Lref)^2          ; n^2
    1 : RETURN, (v_ref*rho_ref/Lref)^2         ; upar^2
    2 : RETURN, (Tref*Qref/k_B*rho_ref/Lref)^2 ; Tpar^2
    3 : RETURN, (Tref*Qref/k_B*rho_ref/Lref)^2 ; Tperp^2
    4 : RETURN, (D_gb*nref/Lref)           ; G_es
    5 : RETURN, (D_gb*nref/Lref)           ; G_em
    6 : RETURN, (D_gb*nref*Tref*Qref/Lref) ; Q_es
    7 : RETURN, (D_gb*nref*Tref*Qref/Lref) ; Q_em
    8 : RETURN, (D_gb*nref*mref*v_ref/Lref) ; P_es
    9 : RETURN, (D_gb*nref*mref*v_ref/Lref) ; P_em
    ELSE : printerror, 'error in nrg_renorm'
  ENDCASE

END

;##########################################################################

FUNCTION get_nrg_string, var_index, fancy=fancy, units=units, ounit=ounit,$   
  time=time, flux=flux, normstr=normstr, rhostr=rho_str, vel_str=vel_str
; takes an index or an array of indices and returns the corresponding
; label(s); if /units is specified, fancy versions are returned
; ounit: returns only the unit in fancy style

  COMMON global_vars

  rho_sc = 0
  isp = 0

  WHILE (isp LT par.n_spec) AND (rho_sc EQ 0) DO BEGIN       ;look for current normalization requested
    IF ((spec[isp].charge EQ -1) AND $
      (spec[isp].temp/series.Tref EQ 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN BEGIN
       rho_sc = 3               ;normalized on electron mass and temperature
       mref_str='m!De!N'
    ENDIF
    IF ((spec[isp].charge EQ 1) AND $
      (spec[isp].temp/series.Tref EQ 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN BEGIN
       rho_sc = 2               ;normalized on ion  mass and temperature
       mref_str='m!Di!N'
    ENDIF
    IF ((spec[isp].charge EQ 1) AND $
      (spec[isp].temp/series.Tref NE 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN BEGIN
       rho_sc = 1
       mref_str='m!Di!N'
    ENDIF
    IF ((spec[isp].charge EQ -1) AND $
      (spec[isp].temp/series.Tref NE 1.0) AND $
      (spec[isp].mass/series.mref EQ 1.0)) THEN BEGIN
       rho_sc = 4
       mref_str='m!De!N'
    ENDIF
    isp += 1
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
    4 : BEGIN  
          vel_str = '(T!Di!N/m!De!N)!U(1/2)!N'
          sp_str = 'e'
          rho_str = '!7q!6!De!N'
         
    END
    ELSE: BEGIN
          vel_str = 'v!Dt,' + ref_str + '!N'
          sp_str = ''
          rho_str = '!7q!6!D' + ref_str + '!N'
          mref_str = 'm!Dref!N'
    END
  ENDCASE

  Lref_str = series.Lref_str

  v = STRARR(14)
  vf = STRARR(14)
  vu = STRARR(14)

  v[0] = '<|n|^2>'
  v[1] = '<|u|^2>'
  v[2] = '<|T_par|^2>'
  v[3] = '<|T_perp|^2>'
  v[4] = 'G_es'
  v[5] = 'G_em'
  v[6] = 'Q_es'
  v[7] = 'Q_em'
  v[8] = 'P_es'
  v[9] = 'P_em'
  v[10] = 'D_es'
  v[11] = 'D_em'
  v[12] = 'Chi_es'
  v[13] = 'Chi_em'
  vf[0] = '!6n'                             
  vf[1] = '!6u!D!9#!6!N'  ;more beautiful written strings
  vf[2] = '!6T!D!9#!6!N'
  vf[3] = '!6T!D!9x!6!N'
  vf[0:3] = '!9!!' + vf[0:3] + '!9!!!6!U2!N'
  vf[4] = '!7C!6!Des!N'
  vf[5] = '!7C!6!Dem!N'
  vf[6] = '!6Q!Des!N'
  vf[7] = '!6Q!Dem!N'
  vf[8] = '!7P!6!Des!N'
  vf[9] = '!7P!6!Dem!N'
  vf[10] = '!6D!Des!N'
  vf[11] = '!6D!Dem!N'
  vf[12] = '!7v!6!Des!N'
  vf[13] = '!7v!6!Dem!N'
  vf = '!12<' + vf + '!12>!6!N'

  sp_str = 'j' ;as n, T and upar are normalized to current species!!!
  vu[0] = 'n!D'+sp_str+'0!N'+rho_str+'/'+Lref_str
  vu[1] = 'v!DT'+sp_str+'!N'+rho_str+'/'+Lref_str
  vu[2] = 'T!D'+sp_str+'0!N'+rho_str+'/'+Lref_str
  vu[3] = 'T!D'+sp_str+'0!N'+rho_str+'/'+Lref_str
  vu[0:3] = '('+vu[0:3]+')!U2!N'
  vu[4:5] = vel_str+rho_str+'!U2!N'+'n!De0!N/'+Lref_str+'!U2!N'
  vu[6:7] = vel_str+rho_str+'!U2!N'+'p!Dref0!N/'+Lref_str+'!U2!N'
  vu[8:9] = mref_str+vel_str+'!U2!N'+rho_str+'!U2!N'+'n!Dref0!N/'+Lref_str+'!U2!N'
  vu[10:13]= vel_str+rho_str+'!U2!N/'+Lref_str+'!N'    

  IF KEYWORD_SET(ounit) THEN RETURN, vu[var_index]
  IF KEYWORD_SET(units) THEN RETURN, vf[var_index] + $
    ' / (' + vu[var_index] + ')'
  IF KEYWORD_SET(fancy) THEN RETURN, vf[var_index] ELSE $
    RETURN, v[var_index]

END

;#############################################################################

FUNCTION get_nrg_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  file_nrg = get_file_path('nrg',run)

  IF NOT file_nrg.exist THEN RETURN, 0.0 ;-1L

  OPENR, nrg_lun, file_nrg.path, /GET_LUN, ERROR=err
  IF err NE 0 THEN printerror, 'strange error occured in get_nrg_time'
  IF EOF(nrg_lun) THEN BEGIN
    PRINT, 'empty nrg file'
    RETURN, 0.0
  ENDIF
  nrg_size = (FSTAT(nrg_lun)).SIZE
  ; due to a mysterious error on PSI this line needs to be repeated:
  nrg_size = (FSTAT(nrg_lun)).SIZE
  time = 0.0D
  line = ''
  n_spec = 0

  IF KEYWORD_SET(first) THEN BEGIN
    READF, nrg_lun, time
    FREE_LUN, nrg_lun
    time *= time_renorm(0)
    RETURN, time
  ENDIF

  READF, nrg_lun, line
  time_length = STRLEN(line)
  line_length = time_length + 1
  content_byte_size = 0
  WHILE NOT EOF(nrg_lun) AND (line_length GT time_length) DO BEGIN
    READF, nrg_lun, line
    line_length = STRLEN(line)
    IF line_length GT time_length THEN n_spec += 1
    content_byte_size += line_length + 1 ; + 1 because of the newline character
  ENDWHILE

  IF (par.n_spec GT 0) AND (n_spec NE par.n_spec) THEN $
    printerror, 'Warning: number of species in parameter file (' + $
    rm0es(par.n_spec) + ') does not match number of species in nrg file (' + $
    rm0es(n_spec) + ')!'

  nrg_steps = FLOOR(nrg_size/content_byte_size)

  IF KEYWORD_SET(get_n_steps) THEN BEGIN
     FREE_LUN, nrg_lun
     RETURN, nrg_steps
  ENDIF

  IF KEYWORD_SET(last) THEN BEGIN
    POINT_LUN, nrg_lun, nrg_size - content_byte_size
    READF, nrg_lun, time
    FREE_LUN, nrg_lun
    time *= time_renorm(0)
    RETURN, time
  ENDIF

END

;##########################################################################

PRO read_nrg, run_label=run_label

  COMMON global_vars

  IF KEYWORD_SET(run_label) THEN BEGIN
    n_runs = 1
    run_labels = [run_label]
  ENDIF ELSE BEGIN
    n_runs = series.n_runs
    run_labels = STRSPLIT(series.run_labels,',',/EXTRACT)
  ENDELSE

  file_nrg = get_file_path('nrg',$
    WHERE(STRSPLIT(series.run_labels,',',/EXTRACT) EQ run_labels))

  IF NOT file_nrg.exist THEN BEGIN
    printerror, 'no nrg file found'
    RETURN
  ENDIF

  IF PTR_VALID(nrg) THEN PTR_FREE, nrg

  ; prepare sparse reading
  nrg_steps_all = 0L
  FOR run = 0, n_runs - 1 DO $
    nrg_steps_all += get_nrg_time(run,/get_n_steps)
  n_steps_sparse_all = 2000L
  interleave_nrg = ROUND(nrg_steps_all/FLOAT(n_steps_sparse_all))

  FOR run = 0, n_runs - 1 DO BEGIN
    nrg_steps = get_nrg_time(run,/get_n_steps)

    ; create an nrg_steps-dimensional struct array with time and data
    ; (the nrg file contains blocks with one time entry and 8 (or 10) fields
    ; of moments per species)
    IF nrg_steps_all LT 20000L THEN BEGIN ; read every step
      new_nrg_data = REPLICATE($
        {time:0.0,data:FLTARR(par.nrgcols,par.n_spec)},nrg_steps)

      OPENR, nrg_lun, file_nrg.path[run], /GET_LUN
      READF, nrg_lun, new_nrg_data
      FREE_LUN, nrg_lun
    ENDIF ELSE BEGIN ; for large files, read only a subset of steps
      IF (!QUIET NE 1) AND (run EQ 0) THEN PRINT, 'reading ' + $
        rm0es(n_steps_sparse_all) + ' out of ' + $
        rm0es(nrg_steps_all) + ' steps available from the nrg files'

      n_steps_sparse = FLOOR(nrg_steps/interleave_nrg)

      single_step = {time:0.0,data:FLTARR(par.nrgcols,par.n_spec)}

      OPENR, nrg_lun, file_nrg.path[run], /GET_LUN
      POINT_LUN, - nrg_lun, position_before
      READF, nrg_lun, single_step
      POINT_LUN, - nrg_lun, position_after
      step_size = position_after - position_before
      FREE_LUN, nrg_lun

      new_nrg_data = REPLICATE($
        {time:0.0,data:FLTARR(par.nrgcols,par.n_spec)},n_steps_sparse)

      count = 0L
      OPENR, nrg_lun, file_nrg.path[run], /GET_LUN
      POINT_LUN, - nrg_lun, position_before
      WHILE NOT EOF(nrg_lun) AND (count LT n_steps_sparse) DO BEGIN
        READF, nrg_lun, single_step
        new_nrg_data[count] = single_step
        count += 1

        position_after = position_before + 1LL * step_size * $
          count * interleave_nrg
        POINT_LUN, nrg_lun, position_after
      ENDWHILE
      FREE_LUN, nrg_lun

      (*par.istep_nrg)[run] *= interleave_nrg
      nrg_steps = n_steps_sparse
    ENDELSE

    IF run EQ 0 THEN nrg_data = TEMPORARY(new_nrg_data) ELSE BEGIN
      IF new_nrg_data[0].time GT nrg_data[N_ELEMENTS(nrg_data.time)-1].time THEN $
        nrg_data = [TEMPORARY(nrg_data),new_nrg_data] $
      ELSE IF new_nrg_data[0].time EQ nrg_data[N_ELEMENTS(nrg_data.time)-1].time THEN $
        nrg_data = [TEMPORARY(nrg_data),new_nrg_data[1:nrg_steps-1]] $
      ELSE BEGIN
        IF !QUIET NE 1 THEN PRINT, 'overlapping nrg time, skipping older data'
        IF nrg_data[0].time GE new_nrg_data[0].time THEN $
          nrg_data = TEMPORARY(new_nrg_data) ELSE nrg_data = $
          [nrg_data[WHERE(nrg_data.time LT new_nrg_data[0].time)],new_nrg_data]
      ENDELSE
    ENDELSE
  ENDFOR

  ; renormalize
  nrg_data.time *= time_renorm(0)
  FOR ivar = 0, par.nrgcols-1 DO $
    nrg_data.data[ivar,*] *= nrg_renorm(ivar)

  nrg = PTR_NEW(nrg_data) ; store in global vars

  file_nrg.exist = 1

END

;##########################################################################

PRO jump_to_nrg_step, jump_step, nrg_lun

;  COMMON global_vars

END

;##########################################################################

PRO read_nrg_step, field_lun, mom_luns

;  COMMON global_vars

END

;##########################################################################

PRO nrg_loop, diag0

  COMMON global_vars

  series.step_count = 1
  read_nrg

END

;##########################################################################

PRO plot_nrg_data, diag=diag, tcorr=tcorr, corr_field=corr_field, $
  show_SI=show_SI

  COMMON global_vars

  IF NOT PTR_VALID(nrg) THEN BEGIN
    PRINTERROR, 'no nrg data found'
    RETURN
  ENDIF
  nrg_steps = N_ELEMENTS((*nrg).time)
  IF nrg_steps LE 1 THEN BEGIN
    printerror, 'too few time points in nrg file'
    RETURN
  ENDIF
  nrgcols = par.nrgcols
  nfluxes = nrgcols - 4

  rel_fit_area = [0.0,1.0]

  ; fields: n^2, u_par^2, T_par^2,
  ; T_perp^2, Gamma_es, Gamma_em, Q_es,
  ; Q_em; P_es, P_em

  IF file_par.exist EQ 1 THEN BEGIN
    n_spec = par.n_spec
    log_plot = 1 - par.nonlinear
  ENDIF ELSE BEGIN
    n_spec = (SIZE((*nrg)[0].data))[2] / nrgcols
    log_plot = 1
    PRINT, 'par file(s) not available, selecting logarithmic y axis'
  ENDELSE
  
  fields = INDGEN((nrgcols EQ 8) OR log_plot ? 4 : 6)
  IF NOT log_plot THEN fields += 4
  IF KEYWORD_SET(diag) THEN fields = INDGEN(nrgcols)
  field_dim = (SIZE(fields))[1]

  ; get x range in time window; exclude first step
  IF KEYWORD_SET(diag) THEN BEGIN
    t = 0L
    WHILE (t LE (nrg_steps - 2)) AND ((*nrg)[t].time LT gui.out.start_t) DO t += 1
    step_start = t
    IF (*nrg)[nrg_steps-1].time LE gui.out.end_t THEN step_end = nrg_steps - 1 $
      ELSE BEGIN

      WHILE (t LE (nrg_steps - 2)) AND ((*nrg)[t].time LT gui.out.end_t) DO t += 1
      step_end = t
    ENDELSE
  ENDIF ELSE BEGIN
    step_start = 1
    step_end = nrg_steps - 1
  ENDELSE
  
  IF step_start EQ step_end THEN BEGIN
    printerror, 'chosen nrg time window too small or out of range'
    RETURN
  ENDIF

  IF KEYWORD_SET(diag) THEN BEGIN
    step_range = step_end - step_start
    fit_start = FLOOR(rel_fit_area[0]*step_range) + step_start < nrg_steps
    fit_end = FLOOR(rel_fit_area[1]*step_range) + step_start < nrg_steps
    IF fit_start GT fit_end THEN BEGIN
      temp_fit_start = fit_start
      fit_start = fit_end
      fit_end = temp_fit_start
    ENDIF

    gr_vals = DBLARR(n_spec*nrgcols)     ; growth rates
    sat_vals = DBLARR(n_spec*nrgcols)    ; saturation averages 
    sat_err_vals = DBLARR(n_spec*nrgcols)
    gr_string = STRARR(n_spec*nrgcols)
    avg_string = STRARR(n_spec*nrgcols)

    ; averages of D_es, D_em, Chi_es, and Chi_em
    sat_vals2 = DBLARR(n_spec*4)
    sat_err_vals2 = DBLARR(n_spec*4)
   
    request_err = tcorr NE 0

    ; compute correlation time self-consistently from all nrg columns
    IF tcorr EQ -1 THEN BEGIN
      cfield = corr_field[0] EQ -1 ? INDGEN(n_spec*nrgcols) : $
        corr_field[*,0] + nrgcols * corr_field[*,1]
      ncf = N_ELEMENTS(cfield)

      ; if too many nrg steps, use sparse data for correlations
      nrg_data = TRANSPOSE((*nrg)[step_start:step_end].data[cfield])
      nrg_time = (*nrg)[step_start:step_end].time

      ; if too many nrg steps, use sparse data for correlations
      IF step_range GE 10000 THEN BEGIN
        PRINT, 'using sparse nrg data to compute t_corr'

        n_steps = step_range + 1L
        WHILE n_steps GT 10000 DO n_steps /= 2
        selinds = LINDGEN(n_steps) * ((step_range + 1) / n_steps)
        nrg_data = nrg_data[selinds,*]
        nrg_time = nrg_time[selinds]
      ENDIF ELSE n_steps = step_range + 1

      nrg_data -= REBIN(REFORM(TOTAL(nrg_data,1)/n_steps,[1,ncf]),[n_steps,ncf])

      corr = pf_arr([n_steps,ncf])
      FOR t = 0, n_steps - 1 DO corr[t,*] = TOTAL(nrg_data*$
        nrg_data[(INDGEN(n_steps)+t) MOD n_steps,*],1)
      corr /= REBIN(REFORM(TOTAL(nrg_data^2,1),[1,ncf]),[n_steps,ncf])

      tcorr = 0.0
      n_tcorr_valid = 0
      FOR spec_field = 0, ncf - 1 DO BEGIN
        tind = (WHERE(corr[*,spec_field] LE 0.36788))[0]
        IF tind GT 0 THEN BEGIN
          tind_eth = tind - (0.36788 - corr[tind,spec_field]) / $
            (corr[tind-1,spec_field] - corr[tind,spec_field])
          dt = nrg_time[tind] - nrg_time[tind-1]

          IF !QUIET NE 1 THEN PRINT, 'nrg correlation time: ' + $
            rm0es(tind_eth*dt) + ' in time units for nrg/spec field ' + $
            rm0es(spec_field)

          tcorr += tind_eth * dt
          n_tcorr_valid += 1
        ENDIF ELSE IF !QUIET NE 1 THEN PRINT, $
          'nrg correlation time calculation failed, nrg/spec field ' + $
          rm0es(spec_field)
      ENDFOR

      tcorr = n_tcorr_valid GT 0 ? tcorr / n_tcorr_valid : 0
      IF tcorr EQ 0 THEN $
        PRINT, 'could not extract correlation time from nrg data'
    ENDIF

    FOR spec_field = 0, n_spec * nrgcols - 1 DO BEGIN
      ; growth rate, gr_fit[1]: log(y) = Ax + B, [A]
      fit_data = (*nrg)[fit_start:fit_end].data[spec_field]
      IF (WHERE(ABS(fit_data) GT 1e-15))[0] EQ -1 THEN gr_fit = [0,0] ELSE BEGIN
        fit_data = TEMPORARY(fit_data[WHERE(ABS(fit_data) GT 1e-15)])
        fit_time = ((*nrg)[fit_start:fit_end].time)$
          [WHERE(ABS(fit_data) GT 1e-15)]
        gr_fit = LINFIT(fit_time,ALOG(fit_data))
      ENDELSE
      gr_vals[spec_field] = gr_fit[1]

      ; compute saturated values with statistical errors
      stat_avg_err, (*nrg)[fit_start:fit_end].time, $
        (*nrg)[fit_start:fit_end].data[spec_field], $
        sat_val, sat_err_val, tcorr=tcorr, wwidth=wwidth
      sat_vals[spec_field] = sat_val
      sat_err_vals[spec_field] = sat_err_val

      gr_string[spec_field] = spec[FIX(spec_field/float(nrgcols))].name + $
        '!6 growth rate ' + get_nrg_string(spec_field MOD nrgcols,/fancy) + ': ' + $
        rm0es(gr_vals[spec_field]) + ' ' + get_var_string(-1,/time,/ounit)
      avg_string[spec_field] = spec[FIX(spec_field/float(nrgcols))].name + $
        '!6 saturation ' + get_nrg_string(spec_field MOD nrgcols,/fancy) + ': ' + $
        rm0es(sat_vals[spec_field])
    ENDFOR
   
    gr_max = MAX(gr_vals[0:3],MIN=gr_min)
    IF gr_min / gr_max GT 0.98 THEN $
      gr_str_tot = '!6series growth rate: ' + rm0es(gr_vals[0]/2.0) + $
      ' ' + get_var_string(-1,/time,/ounit) $
      ELSE gr_str_tot = '!6series growth rate: not converged'

    counter = 0.0D
    spec_nr = 0
    IF NOT par.x_local THEN printerror, $
      'Warning: the nrg diagnostic is not adapted to global simulations'
    ; note: need to consider gxx(x,z) as well as density and temperature
    ; profiles
    IF NOT par.y_local THEN printerror, $
      'Warning: the nrg diagnostic is not adapted to global simulations'
    ; note: need to consider gxx(y,z) 

    ; g_xx (Q*grad x) or sqrt(gxx) (Q*grad x/|grad x|) is averaged over a flux surface
    read_geometry
    IF par.norm_flux_projection THEN BEGIN
      geo_fac = TOTAL(sqrt((*series.geom).gxx)*(*series.geom).jac_norm) / $
        N_ELEMENTS((*series.geom).jac_norm)
    ENDIF ELSE BEGIN
      geo_fac = TOTAL((*series.geom).gxx*(*series.geom).jac_norm) / $
        N_ELEMENTS((*series.geom).jac_norm)
    ENDELSE

    FOR spec_field = 0, n_spec * nrgcols - 1 DO BEGIN
      spec_nr = (spec_field - (spec_field MOD nrgcols)) / nrgcols ; to get number of current species
      IF (spec_field MOD nrgcols) EQ 4 THEN BEGIN ; G_es
        diffusion_factor = series.Lref / (spec[spec_nr].omn * $
          spec[spec_nr].dens * series.nref * geo_fac)
        sat_vals2[counter] = spec[spec_nr].omn NE 0.0 ? $
          sat_vals[spec_field] * diffusion_factor : !VALUES.F_NAN
        sat_err_vals2[counter] = spec[spec_nr].omn NE 0.0 ? $
          sat_err_vals[spec_field] * diffusion_factor : !VALUES.F_NAN          
        counter += 1.0D
      ENDIF
      IF (spec_field MOD nrgcols) EQ 5 THEN BEGIN ; G_em
        diffusion_factor = series.Lref / (spec[spec_nr].omn * $
          spec[spec_nr].dens * series.nref * geo_fac)
        sat_vals2[counter] = spec[spec_nr].omn NE 0.0 ? $
          sat_vals[spec_field] * diffusion_factor : !VALUES.F_NAN
        sat_err_vals2[counter] = spec[spec_nr].omn NE 0.0 ? $
          sat_err_vals[spec_field] * diffusion_factor : !VALUES.F_NAN          
        counter += 1.0D        
        ; note: values of series.*ref are always '1' except
        ; if SI is requested; exception: Lref
      ENDIF
      IF (spec_field MOD nrgcols) EQ 6 THEN BEGIN ; Q_es
        diffusion_factor = series.Lref / (spec[spec_nr].omt * $
          spec[spec_nr].dens * series.nref * spec[spec_nr].temp * series.Tref * $
          series.Qref * geo_fac)      
        sat_vals2[counter] = spec[spec_nr].omt NE 0.0 ? $
          sat_vals[spec_field] * diffusion_factor : !VALUES.F_NAN
        sat_err_vals2[counter] = spec[spec_nr].omt NE 0.0 ? $
          sat_err_vals[spec_field] * diffusion_factor : !VALUES.F_NAN 
        counter += 1.0D
      ENDIF
      IF (spec_field MOD nrgcols) EQ 7 THEN BEGIN ; Q_em
        diffusion_factor = series.Lref / (spec[spec_nr].omt * $
          spec[spec_nr].dens * series.nref * spec[spec_nr].temp * series.Tref * $
          series.Qref * geo_fac)      
        sat_vals2[counter] = spec[spec_nr].omt NE 0.0 ? $
          sat_vals[spec_field] * diffusion_factor : !VALUES.F_NAN
        sat_err_vals2[counter] = spec[spec_nr].omt NE 0.0 ? $
          sat_err_vals[spec_field] * diffusion_factor : !VALUES.F_NAN 
        counter += 1.0D         
      ENDIF
    ENDFOR

    set_output, diag, /ps, multi=[0,1,1], xsize=25, ysize=17
  ENDIF ELSE set_output, /reset, multi=(n_spec NE 4 ? $
    [0,n_spec<3,1+(n_spec-1)/3] : [0,2,2])

  FOR sp = 0, n_spec - 1 DO BEGIN
    spec_title = spec[sp].name
    n_pages   = 1 + (field_dim GT nrgcols-4)
    fld_pp    = n_pages GT 1 ? [4,field_dim-4] : [field_dim]
    fld_start = n_pages GT 1 ? [0,4] : [0]
    fld_end   = n_pages GT 1 ? [3,nrgcols-1] : [field_dim-1]

    FOR page = 0, n_pages - 1 DO BEGIN ; 1st page: moments, 2nd page: fluxes
      PLOT, (*nrg)[step_start:step_end].time, FINDGEN(step_end-step_start+1), $
      /NODATA, YLOG=log_plot, XTITLE=get_var_string(0,/time,/units), $
      YRANGE=[MIN((*nrg)[step_start:step_end].data[fields[fld_start[page]:fld_end[page]],sp]),$
      MAX((*nrg)[step_start:step_end].data[fields[fld_start[page]:fld_end[page]],sp])], $
      XMARGIN=[10,2], YMARGIN=[3,2], /XSTYLE, /YSTYLE, $
      TITLE=spec_title, COLOR=1, CHARSIZE=1.0

      FOR f = fld_start[page], fld_end[page] DO BEGIN
        CASE fields[f] OF
          8: col = 6
          9: col = 8
          ELSE: col = f - fld_start[page] + 1
        ENDCASE
           
        OPLOT, (*nrg)[step_start:step_end].time, $
          (*nrg)[step_start:step_end].data[fields[f],sp], COLOR=col
      ENDFOR

      ; plot cursor
      yrange = !Y.CRANGE
      IF NOT par.nonlinear THEN yrange = 10^yrange
      OPLOT, gui.out.start_t * [1,1], yrange, COLOR=8, THICK=2
      OPLOT, gui.out.end_t * [1,1], yrange, COLOR=8, THICK=2

      ; plot base line
      IF PRODUCT(yrange) LT 0 THEN OPLOT, !X.CRANGE,[0,0], COLOR=1

      ; plot legend at ordinate
      T3D, /RESET, ROT=[0,0,90], TRANS=[0,-1,0]
      PLOT, [0,1], [0,1], COLOR=1, /NOERASE, /NORMAL, XSTYLE=12, YSTYLE=12, $
        /T3D, POSITION=[0,0,1,1], /NODATA
      FOR f = 0, fld_pp[page] - 1 DO BEGIN
        CASE f OF
          4: col = 6
          5: col = 8
          ELSE: col = f + 1
        ENDCASE
        OPLOT, [0.12+f/FLOAT(fld_pp[page])*0.8,0.12+(f+0.25)/FLOAT(fld_pp[page])*0.8], $
          [0.97,0.97], COLOR=col, /T3D
        XYOUTS, 0.14 + (0.8 * f + 0.2) / FLOAT(fld_pp[page]), 0.965, $
          get_nrg_string(fields[f+fld_start[page]],/fancy), COLOR=1, /NORMAL, /T3D
      ENDFOR

      T3D, /RESET
    ENDFOR ; --- page loop
  ENDFOR ; --- species loop

  IF KEYWORD_SET(diag) THEN BEGIN
    FOR isp = 0,gui.out.n_spec_sel - 1 DO BEGIN
      IF nrgcols EQ 10 THEN BEGIN
        set_output, diag, (*gui.out.spec_select)[isp], $
          header=['t','G_tot','Q_tot','P_tot'], dat=[[((*nrg).time-(*nrg)[0].time)],$
          [(*nrg).data[isp*nrgcols+4]+(*nrg).data[isp*nrgcols+5]],$
          [(*nrg).data[isp*nrgcols+6]+(*nrg).data[isp*nrgcols+7]],$
          [(*nrg).data[isp*nrgcols+8]+(*nrg).data[isp*nrgcols+9]]]
      ENDIF ELSE BEGIN
        set_output, diag, (*gui.out.spec_select)[isp], $
          header=['t','G_tot','Q_tot'], dat=[[((*nrg).time-(*nrg)[0].time)],$
          [(*nrg).data[isp*nrgcols+4]+(*nrg).data[isp*nrgcols+5]],$
          [(*nrg).data[isp*nrgcols+6]+(*nrg).data[isp*nrgcols+7]]]
      ENDELSE
    ENDFOR

    x0 = 0.03        ; x offset left
    x1 = 0.03        ; x offset right
    y0 = 0.18        ; y offset down (including space for additional lines)
    y1 = 0.01        ; y offset up
    csize = 1.0      ; char size
    lh = (1.0 - y0 - y1) / 21.0 ; (0.5 - 2.0 * y0) / 9.5 ; line height
    spp = n_spec < 4 ; species_per_page
    first_colw = 25.0 * !D.X_CH_SIZE / !D.X_SIZE ; width of first column

    ;---- table with linear growth rates ----
    FOR page = 0, CEIL(1.0*n_spec/spp) - 1 DO BEGIN
      PLOT, [0,1], [0,1], XSTYLE=12, YSTYLE=12, $
        XMARGIN=[0,0], YMARGIN=[0,0], /NODATA

      ; title and lines
      XYOUTS, x0, 1.0 - y1 - lh, '!7c!6 / (' + get_var_string(-1,/time,/ounit) + $
        ')', CHARSIZE=csize, COLOR=1
      OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - 1.3 * lh), COLOR=1

      nsp = ((page + 1) * spp < n_spec) - page * spp
      colw = (1.0 - x0 - x1 - first_colw) / nsp ; column width
      FOR isp = 0, nsp - 1 DO BEGIN             ; build species columns
        XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh, $
          spec[isp+page*spp].name, CHARSIZE=csize, COLOR=1
        OPLOT, [1,1] * x0+first_colw+(isp-0.05)*colw, $
          [1.0-y1-(nrgcols+2)*lh,1.0-y1], COLOR=1
      ENDFOR

      ; fill table
      FOR j = 0, nrgcols-1 DO BEGIN ; rows
        XYOUTS, x0, 1.0 - y1 - lh * (j + 2), '!7c!6(' + $
          get_nrg_string(j,/fancy) + ')', CHARSIZE=csize, COLOR=1
        FOR isp = 0, nsp - 1 DO $ ; columns
          XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
          rm0es(gr_vals[j+nrgcols*(isp+page*spp)]), CHARSIZE=csize, COLOR=1
      ENDFOR
      XYOUTS, x0, 1.0-y1-(nrgcols+4)*lh, gr_str_tot, CHARSIZE=csize, COLOR=2
      XYOUTS, x0, 0.01,  '!6calculated in time range [' + $
        rm0es((*nrg)[fit_start].time) + ',' + rm0es((*nrg)[fit_end].time) + $
        '] '+ get_var_string(0,/time,/ounit) + ' (' + $
	rm0es(N_ELEMENTS((*nrg)[fit_start:fit_end].time)) + $
        ' nrg steps)', CHARSIZE=csize, COLOR=4
      IF STRCMP(!D.NAME,'PS') THEN XYOUTS, 0.55 + x0, 0.01, $
        gui.out.data_path, CHARSIZE=csize, COLOR=4
    ENDFOR

    ;---- table with saturation values in code units ----
    FOR page = 0, CEIL(1.0*n_spec/spp) - 1 DO BEGIN

      PLOT, [0,1], [0,1], XSTYLE=12, YSTYLE=12, $
        XMARGIN=[0,0], YMARGIN=[0,0], /NODATA

      ; title and lines
      XYOUTS, x0, 1.0 - y1 - lh, 'average', CHARSIZE=csize, COLOR=1
      OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - 1.3 * lh), COLOR=1
      OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - 5.3 * lh), COLOR=1
      OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - (5.3+nfluxes) * lh), COLOR=1

      nsp = ((page + 1) * spp < n_spec) - page * spp
      colw = (1.0 - x0 - x1 - first_colw) / nsp ; column width
      FOR isp = 0, nsp - 1 DO BEGIN             ; build species columns
        XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh, $
          spec[isp+page*spp].name, CHARSIZE=csize, COLOR=1
        OPLOT, [1,1] * x0+first_colw+(isp-0.05)*colw, $
          [1.0-y1-(nrgcols+6)*lh,1.0-y1], COLOR=1
      ENDFOR

      ; --- moments, fluxes, diffusivities ---
      FOR j = 0, nrgcols-1 DO BEGIN ; rows
        XYOUTS, x0, 1.0 - y1 - lh * (j + 2), $
          get_nrg_string(j,/units), CHARSIZE=csize, COLOR=1
        FOR isp = 0, nsp - 1 DO BEGIN ; columns
          IF sat_err_vals[j+nrgcols*(isp+page*spp)] NE 0 THEN BEGIN 
            errorout = '!9+!6' + rm0es(sat_err_vals[j+nrgcols*(isp+page*spp)],prec=2)
          ENDIF ELSE BEGIN
            errorout = ''
          ENDELSE   
          XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
            rm0es(sat_vals[j+nrgcols*(isp+page*spp)],prec=2) + errorout, $
            CHARSIZE=csize, COLOR=1
        ENDFOR
      ENDFOR

      FOR j = nrgcols, nrgcols + 3 DO BEGIN ; D and chi values
        XYOUTS, x0, 1.0 - y1 - lh * (j + 2), $
          get_nrg_string(10+j-nrgcols,/units), CHARSIZE=csize, COLOR=1
        FOR isp = 0, nsp - 1 DO BEGIN
          IF sat_err_vals2[(j-nrgcols)+4*(isp+page*spp)] NE 0 THEN BEGIN 
            errorout = '!9+!6' + $
              rm0es(sat_err_vals2[(j-nrgcols)+4*(isp+page*spp)],prec=2)
          ENDIF ELSE BEGIN
            errorout = ''
          ENDELSE
          XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
            rm0es(sat_vals2[(j-nrgcols)+4*(isp+page*spp)],prec=2) + errorout, $
            CHARSIZE=csize, COLOR=1
        ENDFOR  
      ENDFOR

      msg = ''
      IF TOTAL([sat_err_vals,sat_err_vals2],/NAN) GT 0 THEN $
        msg = '!9+!6 ' + (tcorr EQ -2 ? $
        'simple standard deviation of unweighted data' : $
        'standard deviation of window averages (window width ' + $
        rm0es(wwidth) + ')') ELSE IF request_err THEN $
        msg = '!6could not obtain errors'
      XYOUTS, x0, 0.05, msg, CHARSIZE=csize, COLOR=1
     
      XYOUTS, x0, 0.01, '!6calculated in time range [' + $
        rm0es((*nrg)[fit_start].time) + ',' + rm0es((*nrg)[fit_end].time) + $
        '] '+ get_var_string(0,/time,/ounit) + ' (' + $
	rm0es(N_ELEMENTS((*nrg)[fit_start:fit_end].time)) + ' nrg steps)', $
	CHARSIZE=csize, COLOR=4

      IF STRCMP(!D.NAME,'PS') THEN XYOUTS, 0.55 + x0, 0.01, $
        gui.out.data_path, CHARSIZE=csize, COLOR=4
    ENDFOR

    IF KEYWORD_SET(show_SI) THEN BEGIN ; --- SI fluxes and diffusivities ---
      IF par.mref EQ 0.0 THEN BEGIN
        mref = WHERE(spec.charge EQ -1) NE -1 ? $ ; identify electrons by charge
          1.0D / spec[WHERE(spec.charge EQ -1)].mass * 9.1e-31 : $
          2.0D * 1.6726231e-27 ; m_D
      ENDIF ELSE mref = par.mref
      series_mref = series.mref EQ 0.0 ? 1.0 : series.mref

      area_ref = (2.0 * !PI)^2 * par.n_pol * (par.Lref / series.Lref)^2
      IF par.x_local THEN BEGIN
        IF par.norm_flux_projection THEN $
          area_ref *= TOTAL((*series.geom).jacobian*SQRT((*series.geom).gxx)) $
        ELSE area_ref *= TOTAL((*series.geom).jacobian)
        area_ref *= (*series.geom).C_y / par.nz0
      ENDIF ELSE BEGIN
        IF par.norm_flux_projection THEN $
          area_ref *= TOTAL(TOTAL((*series.geom).jacobian[par.nx0/2,*]*$
            SQRT((*series.geom).gxx[par.nx0/2,*]),2)) $
         ELSE area_ref *= TOTAL(TOTAL((*series.geom).jacobian[par.nx0/2,*]),2)
         area_ref *= (*series.geom).C_y[par.nx0/2] / par.nz0
      ENDELSE
      IF par.norm_flux_projection THEN area_str = 'A!Dref!N!6' $
        ELSE area_str = "V'"

      v_norm = SQRT(par.Tref*par.Qref/mref) / $
        SQRT(series.Tref*series.Qref/series_mref)
      rho_norm = SQRT(mref*par.Tref/par.Qref) / par.Bref / $
        SQRT(series_mref*series.Tref/series.Qref) * series.Bref
      rho_Lref = rho_norm / (par.Lref / series.Lref)

      D_gb_norm = v_norm * rho_norm^2 / (par.Lref / series.Lref)
      norm_Gamma = D_gb_norm * par.nref / series.nref * series.Lref / par.Lref
      norm_D = norm_Gamma / (par.nref / series.nref * series.Lref / par.Lref)
      norm_Gamma *= par.Qref/series.Qref  ;translates to (W/(eV m^2))
      norm_Q = D_gb_norm * par.nref / series.nref * series.Lref / par.Lref * $
        par.Tref/series.Tref*par.Qref/series.Qref 
      norm_chi = norm_Q / (par.nref / series.nref * $
        series.Lref / par.Lref * par.Tref / series.Tref * par.Qref / series.Qref)
      norm_Pi = par.nref / series.nref * par.Tref / series.Tref * $
        par.Qref / series.Qref * rho_Lref^2

      si_units_str = ['n!Dj0!N','v!DTj0!N','T!Dj0!N','T!Dj0!N',$
        '(W/(eV m!U2!N))','(W/(eV m!U2!N))',$
        '(kW/m!U2!N)','(kW/m!U2!N)','(N/cm)','(N/cm)',$
        '(MW/keV)','(MW/keV)','MW','MW','Nm','Nm',$
        '(m!U2!N/s)','(m!U2!N/s)','(m!U2!N/s)','(m!U2!N/s)']
      si_units_str = ' / ' + si_units_str
      si_units_str[0:3] = '!U1/2!N' + si_units_str[0:3]
      si_units_str[10:15] = ' ' + area_str + si_units_str[10:15]

      si_units_fac = [rho_Lref,rho_Lref,rho_Lref,rho_Lref,$
        norm_Gamma,norm_Gamma,norm_Q*1e-3,norm_Q*1e-3,norm_Pi*100,norm_Pi*100,$
        norm_Gamma*area_ref*1e-3,norm_Gamma*area_ref*1e-3,norm_Q*area_ref*1e-6,$
        norm_Q*area_ref*1e-6,norm_Pi*area_ref,norm_Pi*area_ref,$
        norm_D,norm_D,norm_Chi,norm_Chi]

      ;---- table with saturation values in SI units ----
      FOR page = 0, CEIL(1.0*n_spec/spp) - 1 DO BEGIN

        PLOT, [0,1], [0,1], XSTYLE=12, YSTYLE=12, $
          XMARGIN=[0,0], YMARGIN=[0,0], /NODATA

        ; title and lines
        XYOUTS, x0, 1.0 - y1 - lh, 'average', CHARSIZE=csize, COLOR=1
        OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - 1.3 * lh), COLOR=1
        OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - 5.3 * lh), COLOR=1
        OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - (5.3 + nfluxes) * lh), COLOR=1
        OPLOT, [x0-0.01,1.0-x1], [1,1] * (1.0 - y1 - (5.3 + 2 * nfluxes) * lh), COLOR=1

        nsp = ((page + 1) * spp < n_spec) - page * spp
        colw = (1.0 - x0 - x1 - first_colw) / nsp ; column width
        FOR isp = 0, nsp - 1 DO BEGIN             ; build species columns
          XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh, $
            spec[isp+page*spp].name, CHARSIZE=csize, COLOR=1
          OPLOT, [1,1] * x0 + first_colw + (isp - 0.05) * colw, $
            [1.0-y1-(nrgcols+nfluxes+6)*lh,1.0-y1], COLOR=1
        ENDFOR

        ; --- moments ---
        FOR j = 0, 3 DO BEGIN ; rows
          XYOUTS, x0, 1.0 - y1 - lh * (j + 2), $
            get_nrg_string(j,/fancy) + si_units_str[j], CHARSIZE=csize, COLOR=1
          FOR isp = 0, nsp - 1 DO BEGIN ; columns
            IF sat_err_vals[j+nrgcols*(isp+page*spp)] NE 0 THEN BEGIN 
              errorout = '!9+!6' + $
                rm0es(SQRT(sat_err_vals[j+nrgcols*(isp+page*spp)])*si_units_fac[j],prec=3)
            ENDIF ELSE BEGIN
              errorout = ''
            ENDELSE   
            XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
              rm0es(SQRT(sat_vals[j+nrgcols*(isp+page*spp)])*si_units_fac[j],prec=3)+$
              errorout, CHARSIZE=csize, COLOR=1
          ENDFOR
        ENDFOR

        ; --- fluxes ----
        FOR j = 4, nrgcols-1 DO BEGIN ; rows
          XYOUTS, x0, 1.0 - y1 - lh * (j + 2), $
            get_nrg_string(j,/fancy) + si_units_str[j], CHARSIZE=csize, COLOR=1
          FOR isp = 0, nsp - 1 DO BEGIN ; columns
            nrg_j = j + nrgcols * (isp + page * spp)
            IF sat_err_vals[nrg_j] NE 0 THEN BEGIN 
              errorout = '!9+!6' + rm0es(sat_err_vals[nrg_j]*si_units_fac[j],prec=2)
            ENDIF ELSE BEGIN
              errorout = ''
            ENDELSE   
            XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
              rm0es(sat_vals[nrg_j]*si_units_fac[j],prec=2) + errorout, $
              CHARSIZE=csize, COLOR=1
          ENDFOR
        ENDFOR

        ; --- fluxes * area ----
        FOR j = nrgcols, nrgcols+nfluxes-1 DO BEGIN ; rows
          fl_j = j-nrgcols+10
          XYOUTS, x0, 1.0 - y1 - lh * (j + 2), $
            get_nrg_string(j-nrgcols+4,/fancy) + si_units_str[fl_j], CHARSIZE=csize, COLOR=1
          FOR isp = 0, nsp - 1 DO BEGIN ; columns
            nrg_j = j - nrgcols + 4 + nrgcols * (isp + page * spp)
            IF sat_err_vals[nrg_j] NE 0 THEN BEGIN 
              errorout = '!9+!6' + rm0es(sat_err_vals[nrg_j]*si_units_fac[fl_j],prec=2)
            ENDIF ELSE BEGIN
              errorout = ''
            ENDELSE   
            XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
              rm0es(sat_vals[nrg_j]*si_units_fac[fl_j],prec=2) + $
              errorout, CHARSIZE=csize, COLOR=1
          ENDFOR
        ENDFOR

        ; --- diffusivities ----
        FOR j = nrgcols + nfluxes, nrgcols + nfluxes + 3 DO BEGIN
          diff_j = j - nrgcols - nfluxes + 16
          XYOUTS, x0, 1.0 - y1 - lh * (j + 2), $
            get_nrg_string(diff_j-6,/fancy) + si_units_str[diff_j], $
            CHARSIZE=csize, COLOR=1
          FOR isp = 0, nsp - 1 DO BEGIN 
            nrg_j = j - (nrgcols + nfluxes) + 4 * (isp + page * spp)
            IF sat_err_vals2[nrg_j] NE 0 THEN BEGIN 
              errorout = '!9+!6' + rm0es(sat_err_vals2[nrg_j]*si_units_fac[diff_j],prec=2)
            ENDIF ELSE BEGIN
              errorout = ''
            ENDELSE 
             XYOUTS, x0 + first_colw + isp * colw, 1.0 - y1 - lh * (j + 2), $
               rm0es(sat_vals2[nrg_j]*si_units_fac[diff_j],prec=2) + errorout, $
               CHARSIZE=csize, COLOR=1
          ENDFOR
        ENDFOR

        m_str = rm0es(mref,prec=2)
        base_exp = STRSPLIT(m_str,'e',/EXTRACT)
        IF N_ELEMENTS(base_exp) EQ 2 THEN $
          m_str = base_exp[0] + '!9X!610!U' + base_exp[1] + '!N'
        T_str = rm0es(par.Tref,prec=2)
        base_exp = STRSPLIT(T_str,'e',/EXTRACT)
        IF N_ELEMENTS(base_exp) EQ 2 THEN $
          T_str = base_exp[0] + '!9X!610!U' + base_exp[1] + '!N'
        Q_str = rm0es(par.Qref,prec=2)
        base_exp = STRSPLIT(Q_str,'e',/EXTRACT)
        IF N_ELEMENTS(base_exp) EQ 2 THEN $
          Q_str = base_exp[0] + '!9X!610!U' + base_exp[1] + '!N'
        B_str = rm0es(par.Bref,prec=2)
        base_exp = STRSPLIT(B_str,'e',/EXTRACT)
        IF N_ELEMENTS(base_exp) EQ 2 THEN $
          B_str = base_exp[0] + '!9X!610!U' + base_exp[1] + '!N'
        L_str = rm0es(par.Lref,prec=3)
        base_exp = STRSPLIT(L_str,'e',/EXTRACT)
        IF N_ELEMENTS(base_exp) EQ 2 THEN $
          L_str = base_exp[0] + '!9X!610!U' + base_exp[1] + '!N'
        n_str = rm0es(par.nref,prec=2)
        base_exp = STRSPLIT(n_str,'e',/EXTRACT)
        IF N_ELEMENTS(base_exp) EQ 2 THEN $
          n_str = base_exp[0] + '!9X!610!U' + base_exp[1] + '!N'

        XYOUTS, x0, 0.13,  '!6SI ref. values: m = ' + $
          m_str + ' kg, T = ' + T_str + ' eV, Q = ' + Q_str + ' C, B = ' + $
          B_str + ' T, L = ' + L_str + ' m, n = ' + n_str + ' m!U-3!N', $
          CHARSIZE=csize, COLOR=2

        tmp = get_nrg_string(0,rhostr=rhostr,vel_str=vel_str)

        XYOUTS, x0, 0.09, '!6derived ref. values:  ' + rhostr + $
          ' = ' + rm0es(1000*rho_norm,prec=3) + ' mm, ' + vel_str + $
          ' = ' + rm0es(0.001*v_norm,prec=3) + ' km/s, ' + 'Q!Igb!N' + $
          ' = ' + rm0es(norm_Q*1E-3, prec=3) + ' kW/m!U2!N, ' + area_str + $
          ' = ' + rm0es(area_ref, prec=3) + ' m!U2!N', CHARSIZE=csize, COLOR=2

        msg = ''
        IF TOTAL([sat_err_vals,sat_err_vals2],/NAN) GT 0 THEN $
          msg = '!9+!6 ' + (tcorr EQ -2 ? $
          'simple standard deviation of unweighted data' : $
          'standard deviation of window averages (window width ' + $
          rm0es(wwidth) + ')') ELSE IF request_err THEN $
          msg = '!6could not obtain errors'
        XYOUTS, x0, 0.05, msg, CHARSIZE=csize, COLOR=1

        XYOUTS, x0, 0.01,  '!6calculated in time range [' + $
          rm0es((*nrg)[fit_start].time) + ',' + rm0es((*nrg)[fit_end].time) + $
          '] '+ get_var_string(0,/time,/ounit) + ' (' + $
          rm0es(N_ELEMENTS((*nrg)[fit_start:fit_end].time)) + ' nrg steps)', $
          CHARSIZE=csize, COLOR=4

        IF STRCMP(!D.NAME,'PS') THEN XYOUTS, 0.55 + x0, 0.01, $
          gui.out.data_path, CHARSIZE=csize, COLOR=4
      ENDFOR
    ENDIF
  ENDIF

  set_output, diag, /reset

END
