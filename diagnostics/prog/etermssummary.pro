FUNCTION etermssummary_info

  RETURN, {$
    type      : 'energy',$
    title     : 'Energy terms summary',$
    help_text : ['Plots fields from the energy file.'],$
    ext_vars  : [['vars','0','variables to be plotted; default: '+$
                '0:10']]}

END

;######################################################################

PRO etermssummary_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.vars) LT 1 THEN (*e.vars) = INDGEN(11)

  i = set_internal_vars(diag,{vars : *e.vars})

  fft_format_en, /terms

END

;######################################################################

PRO etermssummary_loop, diag

END

;######################################################################

PRO etermssummary_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF N_ELEMENTS(*eterms_time) LE 1 THEN BEGIN
    PRINT, (*diag).name + ': less than two time steps, skipping output'
    RETURN
  ENDIF

  tinds = WHERE((*eterms_time GE gui.out.start_t) AND $
    (*eterms_time LE gui.out.end_t),tcount)
  IF tcount LT 2 THEN BEGIN
    PRINT, (*diag).name + ': < 2 time steps, skipping output'
    RETURN
  ENDIF

  n_vars = N_ELEMENTS((*i).vars)

  en_var_string = '!6' + ['E!Dtot!N','dE/dt!9!!!6!Dtot!N',$
    'Q (= dE/dt!9!!!6!Ddrive!N)','Q!Dsource!N (nonloc.)',$
    'C (= dE/dt!9!!!6!Dcoll!N)',$
    'dE/dt!9!!!6!Dhyp_v!N','dE/dt!9!!!6!Dhyp_z!N',$
    'dE/dt!9!!!6!Dhyp_xy!N','nonlin','Poisson zv','curv+rest',$
    'test','dE/dt_perstep']

  set_output, diag, /ps

  time = (*eterms_time)[tinds]
  uniq_inds = UNIQ(time)
  time = time[uniq_inds]
  eterms_dat = (*eterms)[tinds,*]
  eterms_dat = eterms_dat[uniq_inds,*]
  tcount = N_ELEMENTS(uniq_inds)

  FOR f = 0, n_vars - 1 DO BEGIN
    v = (*i).vars[f]
    IF f MOD 4 EQ 0 THEN BEGIN
      varr = v
      IF f + 1 LT n_vars THEN varr = [varr,(*i).vars[f+1]]
      IF f + 2 LT n_vars THEN varr = [varr,(*i).vars[f+2]]
      IF f + 3 LT n_vars THEN varr = [varr,(*i).vars[f+3]]
      yrange = [MIN(eterms_dat[*,varr]),MAX(eterms_dat[*,varr])]
    ENDIF

    IF f MOD 4 EQ 0 THEN BEGIN
      PLOT, time, eterms_dat[*,v], COLOR=1, $
        /XSTYLE, XTITLE=get_var_string(1,/time,/fancy,/unit), $
        YRANGE=yrange, POSITION=[0.15,0.25,0.95,0.95]
    ENDIF ELSE OPLOT, time, eterms_dat[*,v], COLOR=((f MOD 4)+1)

    IF (f MOD 4 EQ 3) OR (f EQ n_vars - 1) THEN $
      plot_legend, INDGEN(N_ELEMENTS(varr)) + 1, $
      en_var_string[varr], x_offset=[0.1,0.0,0.95,0.12], per_line=2, $
      csize=1.5+INTARR(N_ELEMENTS(varr))

    avg = INT_TABULATED(time,eterms_dat[*,v],/DOUBLE) / $
      (time[tcount-1] - time[0])
    PRINT, en_var_string[f] + ' average: ' + rm0es(avg)
  ENDFOR

  set_output, diag, header=['t','var '+rm0es((*i).vars[0])], $
    dat=[[time],[REFORM(eterms_dat[*,(*i).vars[0]])]]
  set_output, diag, /reset

END
