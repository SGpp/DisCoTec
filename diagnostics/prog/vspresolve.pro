FUNCTION vspresolve_info

  RETURN, {$
    type      : 'vsp',$
    title     : 'V-space resolution check',$
    help_text : ['Plots Fourier transformed v_parallel and mu space, '+$
                 'v-space integrals to determine convergence.',$
                 'The former may be used to determine whether '+$
                 'unphysical high-k_[vpar,mu] modes exist, '+$
                 'requiring higher resolutions. The latter can tell '+$
                 'the user when a simulation is overresolved (only '+$
                 'from an integration point of view).',$
                 'Note that 1st: simulations have to fulfill all '+$
                 'those requirements at the same time; 2nd: '+$
                 'in the lower plots, apparent convergence at a '+$
                 'certain lower resolution does not mean that one '+$
                 'may use that value, only that the current value '+$
                 '(nv0) should be sufficient for the integrals; '+$
                 'and 3rd: the mu FFT is done on an equidistant grid.'],$
    ext_vars  : [['zind','0','index of parallel slice(s); '+$
                  'default: nz0/2; -1 for average'],$
                 ['log','1','logarithmic z axis for contours']]}

END

;######################################################################

PRO vspresolve_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.zind) NE 1 THEN *e.zind = par.nz0 / 2
  IF *e.zind EQ -1 THEN *e.zind = INDGEN(par.nz0)
  IF N_ELEMENTS(*e.log) NE 1 THEN *e.log = 0

  i = set_internal_vars(diag,{$
    zind       : *e.zind,$
    log        : *e.log,$
    data_id    : PTR_NEW()})

END

;######################################################################

PRO vspresolve_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  data = N_ELEMENTS((*i).zind) GT 1 ? TOTAL((*vsp)[(*i).zind,*,*,*,4]*$
    REBIN((*series.geom).jac_norm/par.nz0,$
    [par.nz0,par.nv0,par.nw0,gui.out.n_spec_sel]),1) : REFORM($
    (*vsp)[(*i).zind,*,*,*,4],[par.nv0,par.nw0,gui.out.n_spec_sel])

  (*i).data_id = $
    time_avg((*i).data_id,data,vsp_time,fwd=gui.out.res_steps)

END

;######################################################################

PRO vspresolve_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  c_levels = 40

  data = time_avg((*i).data_id,/avg,fwd=gui.out.res_steps,tarr=time)
  data = REFORM(data,[par.nv0,par.nw0,gui.out.n_spec_sel,$
    1+gui.out.res_steps*(series.step_count-1)],/OVERWRITE)

  vpar = - par.lv + INDGEN(par.nv0) * 2.0 * par.lv / (par.nv0 - 1.0)
  GetMuWeightsAndKnots, muweight, mu
  mumod = FLTARR(2*par.nw0,/NOZERO)
  mumod[0] = - REVERSE(mu)
  mumod[par.nw0] = mu

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, coltable=33, multi=[0,2,2]

    FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
      IF gui.out.res_steps THEN BEGIN
        IF time[n] GE 0.0 THEN BEGIN
          timeinfo = '!6t = ' + rm0es(time[n]) + ' ' + $
            get_var_string(1,/time,/ounit)
        ENDIF ELSE BEGIN
          gamma_omega = rm0es(get_eigenvalue(-time[n]),prec=3)
          gamom_str = '(' + gamma_omega[0] + ',' + gamma_omega[1] + ')'
          timeinfo = '!6EV ' + rm0es(-time[n]) + ', (!7c,x!6) = ' + $
            gamom_str
        ENDELSE
      ENDIF ELSE timeinfo = '!6t = ' + rm0es(gui.out.start_t) + $
        ' to ' + rm0es(gui.out.end_t) + ' ' + $
        get_var_string(1,/time,/ounit)
      timeinfo += '        z/qR = ' + (N_ELEMENTS((*i).zind) GT 1 ? $
        'avg.' : rm0es((*par.z)[(*i).zind],prec=3))

      ; compute and plot v_par and v_perp FFTs
      vfft = SHIFT(ABS(FFT(data[*,*,isp],DIMENSION=1)),par.nv0/2-1,0)
      wfft = FLTARR(par.nv0,2*par.nw0,/NOZERO)
      wfft[0,0] = REVERSE(data[*,*,isp],2)
      wfft[0,par.nw0] = data[*,*,isp]
      wfft = SHIFT(ABS(FFT(wfft,DIMENSION=2,/OVERWRITE)),0,par.nw0-1)

      lev_col = contour_levels(vfft,c_levels,log=(*i).log)
        
      CONTOUR, vfft, vpar, mu, LEVELS=lev_col[*,0], $
        C_COLORS=lev_col[*,1], /XSTYLE, /YSTYLE, /FILL, $
        XTITLE='!6FFT(v!D!9#!6!N)', YTITLE='!7l!6', $
        TITLE='!6v!9!D#!N!6 FFT'
          
      XYOUTS, 0.05, 0.005, timeinfo, /NORMAL, CHARSIZE=1.0

      lev_col = contour_levels(wfft,c_levels,log=(*i).log)

      CONTOUR, wfft, vpar, mumod, LEVELS=lev_col[*,0], $
        C_COLORS=lev_col[*,1], /XSTYLE, /YSTYLE, /FILL, $
        XTITLE='!6v!D!9#!6!N', YTITLE='!6FFT(!7l!6)', $
        TITLE='!7l!6 FFT'
          
      XYOUTS, 0.05, 0.005, timeinfo, /NORMAL, CHARSIZE=1.0

      ; compute and plot less-resolved integrals
      ; of v_par and v_perp space
      vint_full = TOTAL(data[*,*,isp],1) / par.nv0
      vint_red = FLTARR(par.nv0/2-2,par.nw0,/NOZERO)

      FOR w = 0, par.nw0 - 1 DO BEGIN
        FOR v = 0, par.nv0 / 2 - 3 DO BEGIN
          vmod = 2 * v + 4
          vint_red[v,w] = TOTAL(INTERPOL(data[*,w,isp],vmod),1) / vmod
        ENDFOR
      ENDFOR

      LOADCT, 41, FILE='internal/colortable.tbl'
      PLOT, mu, vint_full, COLOR=1, /XSTYLE, /NODATA, $
        YRANGE=[0,MAX(vint_full)>MAX(vint_red)], $
        TITLE='!6v!9!D#!N!6 integr.'
      FOR v = 0, par.nv0 / 2 - 3 DO BEGIN
        OPLOT, mu, vint_red[v,*], COLOR=v+2
          IF v LT 10 THEN XYOUTS, $
          !X.CRANGE[0] + 0.7 * (!X.CRANGE[1] - !X.CRANGE[0]), $
          (0.94 - 0.08 * v) * (!Y.CRANGE[1] - !Y.CRANGE[0]) + $
          !Y.CRANGE[0], '!6N!Dv!N=' + rm0es(2*v+4), $
          COLOR=v+2, CHARSIZE=1.0
      ENDFOR
      OPLOT, mu, vint_full, COLOR=1, LINE=1
      XYOUTS, !X.CRANGE[0] + 0.7 * (!X.CRANGE[1] - !X.CRANGE[0]), $
        (0.94 - 0.08 * (v < 10)) * (!Y.CRANGE[1] - !Y.CRANGE[0]) + $
        !Y.CRANGE[0], '!6N!7!Dl!N!6=' + rm0es(par.nv0), $
        COLOR=1, CHARSIZE=1.0

      nwarr = [4,6,8,10,16,24,36]
      nwarr = nwarr[WHERE(nwarr LT par.nw0)]
      nnw = N_ELEMENTS(nwarr)

      IF nnw GT 0 THEN BEGIN
        wint_full = $
          TOTAL(data[*,*,isp]*REBIN(REFORM(muweight,[1,par.nw0]),$
          [par.nv0,par.nw0]),2)
        wint_red = FLTARR(par.nv0,nnw,/NOZERO)

        FOR w = 0, nnw - 1 DO BEGIN
          wmod = nwarr[w]
          GetMuWeightsAndKnots,w_mu,k_mu,nw0=wmod,lw=DOUBLE(par.lw)
          FOR v = 0, par.nv0 - 1 DO wint_red[v,w] = $
            TOTAL(REFORM(INTERPOL(data[v,*,isp],mu,k_mu))*w_mu,1)
        ENDFOR

        PLOT, vpar, wint_full, COLOR=1, /XSTYLE, /NODATA, $
          YRANGE=[0,MAX(wint_full)>MAX(wint_red)], $
          TITLE='!7l!6 integr.'
        FOR w = 0, nnw - 1 DO BEGIN
          OPLOT, vpar, wint_red[*,w], COLOR=w+2
          IF w LT 10 THEN XYOUTS, $
            !X.CRANGE[0] + 0.7 * (!X.CRANGE[1] - !X.CRANGE[0]), $
            (0.94 - 0.08 * w) * (!Y.CRANGE[1] - !Y.CRANGE[0]) + $
            !Y.CRANGE[0], '!6N!7!Dl!N!6=' + rm0es(nwarr[w]), $
            COLOR=w+2, CHARSIZE=1.0
        ENDFOR
        OPLOT, vpar, wint_full, COLOR=1, LINE=1
        XYOUTS, !X.CRANGE[0] + 0.7 * (!X.CRANGE[1] - !X.CRANGE[0]), $
          (0.94 - 0.08 * (w < 10)) * (!Y.CRANGE[1] - !Y.CRANGE[0]) + $
          !Y.CRANGE[0], '!6N!7!Dl!N!6=' + rm0es(par.nw0), $
          COLOR=1, CHARSIZE=1.0
      ENDIF

      LOADCT, 33, FILE='internal/colortable.tbl'
    ENDFOR

    set_output, diag, sp, /reset
  ENDFOR ; --- isp loop 

END
