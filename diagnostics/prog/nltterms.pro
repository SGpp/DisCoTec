FUNCTION nltterms_info

  RETURN, {$
    type      : 'nlt',$
    title     : 'NL and LIN time traces',$
    help_text : ['Plots information from nlt_info.dat.'],$
    ext_vars  : [['vars','0','variables to be plotted; default: '+$
                '0:20'],$
                 ['trange','0','time window to plot; default: all'],$
                 ['yrange','0','yrange of the plot; default: min,max']$
                 ]}

END

;######################################################################

PRO nltterms_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.vars) LT 1 THEN (*e.vars) = INDGEN(21)
  IF N_ELEMENTS(*e.trange) LT 1 THEN (*e.trange)=[-1,-1]
  IF N_ELEMENTS(*e.yrange) LT 1 THEN (*e.yrange)=[-1,-1]
  ;print,"trange",*e.trange,(*e.trange)[1]

  i = set_internal_vars(diag,{$
      vars   : *e.vars,$
      trange : *e.trange,$
      yrange : *e.yrange})

END

;######################################################################

PRO nltterms_loop, diag

END

;######################################################################

PRO nltterms_output, diag

  COMMON global_vars

  ;print,'!!!!!!!Start of nltterms_output'
  e = (*diag).external_vars

  IF N_ELEMENTS(*nltterms_time) LE 1 THEN BEGIN
    PRINT, (*diag).name + ': less than two time steps, skipping output'
    RETURN
  ENDIF

  ;n_vars = N_ELEMENTS((*i).vars)
;  en_var_string = '!6' + ['E!Dtot!N','dE/dt!9!!!6!Dtot!N',$
;    'Q (= dE/dt!9!!!6!Ddrive!N)','Q!Dsource!N (nonloc.)',$
;    'C (= dE/dt!9!!!6!Dcoll!N)',$
;    'dE/dt!9!!!6!Dhyp_v!N','dE/dt!9!!!6!Dhyp_z!N',$
;    'dE/dt!9!!!6!Dhyp_xy!N','nonlin','Poisson zv','curv+rest',$
;    'test','dE/dt_perstep']

  set_output, diag, /ps

  kxmin=2.0*!PI/par.lx
  IF par.nlp_gdt THEN BEGIN
    kx_ind=par.nlp_kxind 
    ky_ind=par.nlp_kyind 
    nlt_loop_max=par.num_nlt_pod_modes
  ENDIF ELSE BEGIN
    kx_ind=*par.kx_nlt_ind
    ky_ind=*par.ky_nlt_ind
    nlt_loop_max=par.num_nlt_modes-1
  ENDELSE
  ;print,kx_ind
  ;print,ky_ind
  xrange=FLTARR(2)
  IF (*e.trange)[0] NE -1 THEN BEGIN
   xrange=[(*e.trange)[0],(*e.trange)[1]] 
  ENDIF
  IF (*e.yrange)[0] NE -1 THEN BEGIN
   yrange=[(*e.yrange)[0],(*e.yrange)[1]] 
  ENDIF

  FOR f = 0, nlt_loop_max DO BEGIN

    kxstring=STRING(kx_ind[f]*kxmin,FORMAT='(F5.2)')
    IF kx_ind[f] GT par.nx0/2 THEN kxstring=STRING(-1.0*(par.nx0-kx_ind[f])*kxmin,FORMAT='(F5.2)')
    IF (*e.yrange)[0] EQ -1 THEN yrange = [MIN((*nltterms)[*,f*4:f*4+2]),MAX((*nltterms)[*,f*4:f*4+2])]
    IF par.nlp_gdt THEN BEGIN
        IF f LT par.num_nlt_pod_modes THEN BEGIN
          ptitle='k!Dx!N='+kxstring+'k!Dy!N='$
             +STRING(ky_ind*par.kymin,FORMAT='(F5.2)')+'(POD n='+STRTRIM(f+1,2)+')'
        ENDIF ELSE BEGIN
          ptitle='k!Dx!N='+STRING(kx_ind*kxmin,FORMAT='(F5.2)')+'k!Dy!N='$
             +STRING(ky_ind*par.kymin,FORMAT='(F5.2)')+'(POD n=res.)'
        ENDELSE
    ENDIF ELSE BEGIN
      ptitle='k!Dx!N='+kxstring+'k!Dy!N='$
           +STRING(ky_ind[f]*par.kymin,FORMAT='(F5.2)')
    ENDELSE
    PLOT, *nltterms_time, (*nltterms)[*,f*4+2], COLOR=1, $
         TITLE=ptitle,$
        /XSTYLE, XTITLE=get_var_string(1,/time,/fancy,/unit), $
        XRANGE=xrange,YRANGE=yrange, POSITION=[0.15,0.25,0.95,0.95]
    OPLOT, *nltterms_time, (*nltterms)[*,f*4], COLOR=4
    OPLOT, *nltterms_time, (*nltterms)[*,f*4+1], COLOR=2
    ;OPLOT, *nltterms_time, (*nltterms)[*,0]-(*nltterms)[*,0], COLOR=1
    plot_legend, [1,4,2],['E!Dk!N','dE/dt!D NL !N','dE/dt!D lin !N']  

  ENDFOR

  set_output, diag, /reset

END
