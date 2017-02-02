FUNCTION nrgtstep_info

  RETURN, {$
    type      : 'nrg',$
    title     : 'time step',$
    help_text : ['Displays the (possibly adapted) time step as obtained '+$
                 'from the nrg file together with the dt_max taken '+$
		 'from the parameters file (note: currently only '+$
		 'the last dt_max is displayed)'],$
    ext_vars  : [['var','0','add time trace of an nrg var; default: -1'+$
                  '(none)'],$
                 ['dt_max','1','add dt_max to plot'],$
                 ['log','1','use logarithmic y axis']]}

END

;######################################################################

PRO nrgtstep_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var) NE 1 THEN *e.var = -1
  show_dt_max = KEYWORD_SET(*e.dt_max)
  log = KEYWORD_SET(*e.log)

  i = set_internal_vars(diag,{$
    var         : *e.var,$
    show_dt_max : show_dt_max,$
    log         : log})

END

;######################################################################

PRO nrgtstep_loop, diag

END

;######################################################################

PRO nrgtstep_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  time_nrg_full = (*nrg).time
  t_window_ind = WHERE((time_nrg_full GE gui.out.start_t) AND $
    (time_nrg_full LE gui.out.end_t))
  IF N_ELEMENTS(t_window_ind) LT 2 THEN BEGIN
    PRINT, (*diag).name + " error: insufficiently large time window"
    RETURN
  ENDIF
  time_nrg = time_nrg_full[t_window_ind]

  n_steps = N_ELEMENTS(time_nrg)
  dt_arr = FLTARR(n_steps,/NOZERO)
  dt_arr[0:n_steps-2] = (time_nrg[1:n_steps-1] - time_nrg[0:n_steps-2]) / $
    ABS((*par.istep_nrg)[0])
  dt_arr[n_steps-1] = dt_arr[n_steps-2]
  dt_max = par.dt_max*nrg_renorm(0,/time)

  FOR isp = 0, ((*i).var NE -1) * (gui.out.n_spec_sel - 1) DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    IF (*i).var NE -1 THEN set_output, diag, sp, /ps, multi=[0,1,2] $
      ELSE set_output, diag, /ps

    yrange_x = (*i).log ? MIN(dt_arr) : 0
    yrange = (*i).show_dt_max ? $
      [yrange_x,dt_max>MAX(dt_arr)] : [yrange_x,MAX(dt_arr)]

    PLOT, time_nrg, dt_arr, COLOR=1, YRANGE=yrange, $
      /XSTYLE, XTITLE="!6"+get_var_string(0,/time,/fancy), $
      YTITLE="!7D!6t / ("+get_var_string(0,/time,/ounit)+")", $
      TITLE="!6nrg time step adaptation", YLOG=(*i).log
    IF (*i).show_dt_max THEN OPLOT, $
      [time_nrg[0],time_nrg[n_steps-1]], dt_max * [1,1], $
      COLOR=2, LINE=2

    IF (*i).var NE -1 THEN BEGIN
      var_nrg = (*nrg)[t_window_ind].data[sp*8+(*i).var]

      PLOT, time_nrg, var_nrg, COLOR=1, /XSTYLE, $
        XTITLE="!6"+get_var_string(0,/time,/fancy), $
        YTITLE="!6"+get_nrg_string((*i).var,/fancy), $
        TITLE="!6nrg time trace", YLOG=(*i).log

      set_output, diag, sp, /reset
    ENDIF ELSE set_output, diag, /reset
  ENDFOR
  
END
