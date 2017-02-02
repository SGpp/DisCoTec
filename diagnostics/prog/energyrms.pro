FUNCTION energyrms_info

  RETURN, {$
    type      : 'energy',$
    title     : 'Root mean square',$
    help_text : ['Plots rms of energy quantities.'],$
    ext_vars  : [['vars','0','variables to be plotted; default: [0,1]']]}

END

;######################################################################

PRO energyrms_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.vars) LT 1 THEN (*e.vars) = [0,1]
  valid = WHERE(*e.vars LT 5)
  IF valid[0] GE 0 THEN *e.vars = (*e.vars)[valid]

  i = set_internal_vars(diag,{$
    n_vars  : N_ELEMENTS(*e.vars),$
    vars    : *e.vars,$
    rms_id  : PTR_NEW(),$
    time_id : PTR_NEW()})

  IF par.x_local THEN fft_format_en, kxky=(*e.vars) $
    ELSE fft_format_en, sxky=(*e.vars)

END

;######################################################################

PRO energyrms_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  rmsall = FLTARR((*i).n_vars)

  FOR f = 0, (*i).n_vars - 1 DO BEGIN
    v = (*i).vars[f]

    IF par.x_local THEN BEGIN
      rms = 2.0 * TOTAL(TOTAL(ABS(*energy[v].kxky)^2,1),1) - $
        TOTAL(ABS((*energy[v].kxky)[*,0,*])^2,1)
      rms *= (*series.geom).jac_norm / par.nz0
    ENDIF ELSE BEGIN
      rms = 2.0 * TOTAL(ABS(*energy[v].kxky)^2,2) - $
        ABS((*energy[v].kxky)[*,0,*])^2
      rms *= (*series.geom).jac_norm / (par.nx0 * par.nz0)
    ENDELSE
    rmsall[f] = SQRT(TOTAL(rms))
  ENDFOR

  (*i).rms_id = store_step((*i).rms_id,rmsall)
  (*i).time_id = store_step((*i).time_id,energy_time)

END

;######################################################################

PRO energyrms_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  rms = store_step((*i).rms_id,/get,/reg_array)
  time = store_step((*i).time_id,/get,/reg_array)

  IF N_ELEMENTS(time) LE 1 THEN BEGIN
    PRINT, (*diag).name + ': less than two time steps, skipping output'
    RETURN
  ENDIF

  en_var_string = '!6' + ['E!Dtot!N','dE/dt!9!!!6!Dnonconserve!N',$
    'Q (= dE/dt!9!!!6!Ddrive!N)','C (= dE/dt!9!!!6!Dcoll!N)']

  set_output, diag, /ps

  FOR f = 0, (*i).n_vars - 1 DO $
    PLOT, time, rms[f,*], COLOR=1, /XSTYLE, $
    XTITLE=get_var_string(1,/time,/fancy,/unit), TITLE='!6(!13<!9#!6'+$
    en_var_string(f)+'!9#!6!U2!N!13>!6!Dxyz!N)!U1/2!N(t)'

  set_output, diag, header=['t',get_var_string((*i).vars[0])], $
    dat=[[time],[REFORM(rms[0,*])]]
  set_output, diag, /reset

END
