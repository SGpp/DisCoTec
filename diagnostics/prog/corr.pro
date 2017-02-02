FUNCTION corr_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'Correlations (x/y,t)',$
    help_text : ['Plots spatial (x/y) or temporal cross-correlation',$
                 'functions for pairs of fluctuating quantities.'],$
    ext_vars  : $
      [['var1','0','first variable (see variable list; -1 for B_x, '+$
       '-2 for B_y, -3 for E_x, -4 for E_y); default: phi'],$
      ['var2','0','(array of) second variable(s) (same options '+$
       'as var1); default: density'],$
      ['zind','0','index of z plane, -1 for average; default: nz0/2'],$
      ['tcorr','1','time correlation; sets zind = -1'],$
      ['delzonal','1','remove zonal component before correlation']]}

END

;######################################################################

PRO corr_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var1) NE 1 THEN *e.var1 = 0
  IF N_ELEMENTS(*e.var2) LT 1 THEN *e.var2 = par.n_fields
  IF (TOTAL(WHERE([*e.var1,*e.var2] EQ -1)) + $
    TOTAL(WHERE([*e.var1,*e.var2] EQ -2)) GT 0) AND $
    (par.n_fields LT 2) THEN BEGIN
    (*diag).requested = 0
    printerror, 'Skipping ' + (*diag).name + ': nonexistent magnetic ' + $
      'fields requested'
    RETURN
  ENDIF
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = par.nz0 / 2
  IF NOT KEYWORD_SET(*e.tcorr) THEN *e.tcorr = 0 ELSE *e.zind = -1
  IF NOT KEYWORD_SET(*e.delzonal) THEN *e.delzonal = 0
  IF (*e.zind)[0] EQ -1 THEN *e.zind = INDGEN(par.nz0)

  i = set_internal_vars(diag,{$
    var1        : *e.var1,$
    var2        : *e.var2,$
    zind        : *e.zind,$
    tcorr       : *e.tcorr,$
    delzonal    : *e.delzonal,$
    corr_xy_id  : PTR_NEW(),$
    dat1_t_id   : PTR_NEW(),$
    dat2_t_id   : PTRARR(N_ELEMENTS(*e.var2),gui.out.n_spec_sel),$
    t_arr_id    : PTR_NEW()})

  IF ((WHERE([(*i).var1,(*i).var2] EQ -1))[0] NE -1) OR $
    ((WHERE([(*i).var1,(*i).var2] EQ -2))[0] NE -1) THEN var_req_bxy = 1
  IF (WHERE([(*i).var1,(*i).var2] LE -3))[0] NE -1 THEN var_req_exy = 0
  IF (WHERE([(*i).var1,(*i).var2] GE 0))[0] NE -1 THEN var_req_reg = $
    ([(*i).var1,(*i).var2])[WHERE([(*i).var1,(*i).var2] GE 0)]

  IF (*i).delzonal THEN fft_format, kxky=var_req_reg $
    ELSE fft_format, sxsy=var_req_reg
  fft_format, kxky=var_req_bxy
  fft_format, kxky=var_req_exy

END

;######################################################################

PRO corr_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nvar2 = N_ELEMENTS((*i).var2)
  nz = N_ELEMENTS((*i).zind)

  IF NOT (*i).tcorr THEN corr_xy_all = $
    FLTARR(par.nx0,par.ny0,nvar2,gui.out.n_spec_sel,/NOZERO) ELSE $
    dat1_spec = FLTARR(par.nx0,par.ny0,nz,gui.out.n_spec_sel,/NOZERO)

  jac_norm = NOT (*i).tcorr AND (nz GT 1) ? $
    REBIN(REFORM((*series.geom).jac_norm[(*i).zind],[1,1,nz]),$
    [par.nx0,par.ny0,nz]) : 1

  data1 = FLTARR(par.nkx0,par.ny0,nz,/NOZERO)
  data2 = FLTARR(par.nkx0,par.ny0,nz,/NOZERO)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    IF TOTAL(WHERE([(*i).var1,(*i).var2] EQ -1) GE 0) GT 0 THEN BEGIN
      B_x_kxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      B_x_sxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      temp_data = COMPLEXARR(par.nkx0,2*par.nky0,nz,/NOZERO)

      B_x_kxky[0,0,0] = REBIN(REFORM(*series.ky,[1,par.nky0,1]),$
        [par.nkx0,par.nky0,nz]) * COMPLEX(0,1) * (*mom[0,1].kxky)[*,*,(*i).zind]

      IF (*i).delzonal THEN B_x_kxky[*,0,*] -= $
        TOTAL(REFORM(REFORM(B_x_kxky[*,0,*])*$
        REBIN(REFORM((*series.geom).jac_norm[(*i).zind]/nz,$
        [1,nz]),[par.nkx0,nz]),[par.nkx0,1,nz]),2)

      B_x_sxky[0,0,0] = FFT((B_x_kxky),DIMENSION=1,/INVERSE)
      temp_data[0,0,0] = B_x_sxky
      temp_data[*,par.nky0,*] = 0
      FOR y = par.nky0 + 1, 2 * par.nky0 - 1 DO $
        temp_data[0,y,0] = CONJ(B_x_sxky[*,2*par.nky0-y,*])
      B_x_sxsy = FFT(temp_data,DIMENSION=2,/INVERSE)
    ENDIF
    IF TOTAL(WHERE([(*i).var1,(*i).var2] EQ -2) GE 0) GT 0 THEN BEGIN
      B_y_kxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      B_y_sxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      temp_data = COMPLEXARR(par.nkx0,2*par.nky0,nz,/NOZERO)

      B_y_kxky[0,0,0] = - REBIN(*series.kx,[par.nkx0,par.nky0,nz]) * $
        COMPLEX(0,1) * (*mom[0,1].kxky)[*,*,(*i).zind]

      IF (*i).delzonal THEN B_y_kxky[*,0,*] -= $
        TOTAL(REFORM(REFORM(B_y_kxky[*,0,*])*$
        REBIN(REFORM((*series.geom).jac_norm[(*i).zind]/nz,$
        [1,nz]),[par.nkx0,nz]),[par.nkx0,1,nz]),2)

      B_y_sxky[0,0,0] = FFT((B_y_kxky),DIMENSION=1,/INVERSE)
      temp_data[0,0,0] = B_y_sxky
      temp_data[*,par.nky0,*] = 0
      FOR y = par.nky0 + 1, 2 * par.nky0 - 1 DO $
        temp_data[0,y,0] = CONJ(B_y_sxky[*,2*par.nky0-y,*])
      B_y_sxsy = FFT(temp_data,DIMENSION=2,/INVERSE)
    ENDIF
    IF TOTAL(WHERE([(*i).var1,(*i).var2] EQ -3) GE 0) GT 0 THEN BEGIN
      E_x_kxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      E_x_sxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      temp_data = COMPLEXARR(par.nkx0,2*par.nky0,nz,/NOZERO)

      E_x_kxky[0,0,0] = REBIN(*series.kx,[par.nkx0,par.nky0,nz]) * $
        COMPLEX(0,1) * (*mom[0,0].kxky)[*,*,(*i).zind]

      IF (*i).delzonal THEN E_x_kxky[*,0,*] -= $
        TOTAL(REFORM(REFORM(E_x_kxky[*,0,*])*$
        REBIN(REFORM((*series.geom).jac_norm[(*i).zind]/nz,$
        [1,nz]),[par.nkx0,nz]),[par.nkx0,1,nz]),2)

      E_x_sxky[0,0,0] = FFT((E_x_kxky),DIMENSION=1,/INVERSE)
      temp_data[0,0,0] = E_x_sxky
      temp_data[*,par.nky0,*] = 0
      FOR y = par.nky0 + 1, 2 * par.nky0 - 1 DO $
        temp_data[0,y,0] = CONJ(E_x_sxky[*,2*par.nky0-y,*])
      E_x_sxsy = FFT(temp_data,DIMENSION=2,/INVERSE)
    ENDIF
    IF TOTAL(WHERE([(*i).var1,(*i).var2] EQ -4) GE 0) GT 0 THEN BEGIN
      E_y_kxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      E_y_sxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
      temp_data = COMPLEXARR(par.nkx0,2*par.nky0,nz,/NOZERO)

      E_y_kxky[0,0,0] = REBIN(REFORM(*series.ky,[1,par.nky0,1]),$
        [par.nkx0,par.nky0,nz]) * COMPLEX(0,1) * (*mom[0,0].kxky)[*,*,(*i).zind]

      IF (*i).delzonal THEN E_y_kxky[*,0,*] -= $
        TOTAL(REFORM(REFORM(E_y_kxky[*,0,*])*$
        REBIN(REFORM((*series.geom).jac_norm[(*i).zind]/nz,$
        [1,nz]),[par.nkx0,nz]),[par.nkx0,1,nz]),2)

      E_y_sxky[0,0,0] = FFT((E_y_kxky),DIMENSION=1,/INVERSE)
      temp_data[0,0,0] = E_y_sxky
      temp_data[*,par.nky0,*] = 0
      FOR y = par.nky0 + 1, 2 * par.nky0 - 1 DO $
        temp_data[0,y,0] = CONJ(E_y_sxky[*,2*par.nky0-y,*])
      E_y_sxsy = FFT(temp_data,DIMENSION=2,/INVERSE)
    ENDIF
    IF (*i).var1 GE 0 THEN BEGIN
      IF NOT (*i).delzonal THEN BEGIN
        mom_sxsy = (*mom[isp,(*i).var1].sxsy)[*,*,(*i).zind]
      ENDIF ELSE BEGIN
        mom_kxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
        mom_sxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
        temp_data = COMPLEXARR(par.nkx0,2*par.nky0,nz,/NOZERO)

        mom_kxky[0,0,0] = (*mom[isp,(*i).var1].kxky)[*,*,(*i).zind]

        ; remove zonal component
        FOR z = 0, nz - 1 DO mom_kxky[*,0,z] -= $
          TOTAL(REFORM(REFORM(mom_kxky[*,0,*])*$
          REBIN(REFORM((*series.geom).jac_norm[(*i).zind]/nz,$
          [1,nz]),[par.nkx0,nz]),[par.nkx0,nz]),2)

        mom_sxky[0,0,0] = FFT((mom_kxky),DIMENSION=1,/INVERSE)
        temp_data[0,0,0] = mom_sxky
        temp_data[*,par.nky0,*] = 0
        FOR y = par.nky0 + 1, 2 * par.nky0 - 1 DO $
          temp_data[0,y,0] = CONJ(mom_sxky[*,2*par.nky0-y,*])
        mom_sxsy = FFT(temp_data,DIMENSION=2,/INVERSE)
      ENDELSE
    ENDIF

    CASE (*i).var1 OF
      -1   : data1[0,0,0] = B_x_sxsy
      -2   : data1[0,0,0] = B_y_sxsy
      -3   : data1[0,0,0] = E_x_sxsy
      -4   : data1[0,0,0] = E_y_sxsy
      ELSE : data1[0,0,0] = mom_sxsy
    ENDCASE

    IF (*i).tcorr THEN dat1_spec[0,0,0,isp] = data1

    FOR ivar = 0, nvar2 - 1 DO BEGIN
      ; ### for spatial correlations, no Jacobian ###
      ; ### is currently used for the z averaging ###
      ; ###  However: delzonal uses the Jacobian  ###

      IF (*i).var2[ivar] GE 0 THEN BEGIN
        IF NOT (*i).delzonal THEN BEGIN
          mom_sxsy = (*mom[isp,(*i).var2[ivar]].sxsy)[*,*,(*i).zind]
        ENDIF ELSE BEGIN
          mom_kxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
          mom_sxky = COMPLEXARR(par.nkx0,par.nky0,nz,/NOZERO)
          temp_data = COMPLEXARR(par.nkx0,2*par.nky0,nz,/NOZERO)

          mom_kxky[0,0,0] = (*mom[isp,(*i).var2[ivar]].kxky)[*,*,(*i).zind]

          ; remove zonal component
          FOR z = 0, nz - 1 DO mom_kxky[*,0,z] -= $
            TOTAL(REFORM(REFORM(mom_kxky[*,0,*])*$
            REBIN(REFORM((*series.geom).jac_norm[(*i).zind]/nz,$
            [1,nz]),[par.nkx0,nz]),[par.nkx0,nz]),2)

          mom_sxky[0,0,0] = FFT((mom_kxky),DIMENSION=1,/INVERSE)
          temp_data[0,0,0] = mom_sxky
          temp_data[*,par.nky0,*] = 0
          FOR y = par.nky0 + 1, 2 * par.nky0 - 1 DO $
            temp_data[0,y,0] = CONJ(mom_sxky[*,2*par.nky0-y,*])
          mom_sxsy = FFT(temp_data,DIMENSION=2,/INVERSE)
        ENDELSE
      ENDIF

      CASE (*i).var2[ivar] OF
        -1   : data2[0,0,0] = B_x_sxsy
        -2   : data2[0,0,0] = B_y_sxsy
        -3   : data2[0,0,0] = E_x_sxsy
        -4   : data2[0,0,0] = E_y_sxsy
        ELSE : data2[0,0,0] = mom_sxsy
      ENDCASE

      IF NOT (*i).tcorr THEN BEGIN ; spatial correlation
        IF (*i).var2[ivar] NE (*i).var1 THEN corr_norm_inv = $
          1.0 / SQRT(TOTAL(data1^2*jac_norm)*TOTAL(data2^2*jac_norm)) ELSE $
          corr_norm_inv = 1.0 / TOTAL(data1^2*jac_norm)

        corr_xy = FLTARR(par.nx0,par.ny0,/NOZERO)

        IF nz GT 1 THEN BEGIN
          data1[0,0,0] = data1 * jac_norm
          FOR y = 0, par.ny0 - 1 DO FOR x = 0, par.nx0 - 1 DO corr_xy[x,y] = $
            TOTAL(data1*SHIFT(data2,-x-par.nx0/2,-y-par.ny0/2,0))
        ENDIF ELSE BEGIN
          FOR y = 0, par.ny0 - 1 DO FOR x = 0, par.nx0 - 1 DO corr_xy[x,y] = $
            TOTAL(data1*SHIFT(data2,-x-par.nx0/2,-y-par.ny0/2))
        ENDELSE

        corr_xy_all[0,0,ivar,isp] = TEMPORARY(corr_xy) * corr_norm_inv

      ENDIF ELSE (*i).dat2_t_id[ivar,isp] = $
        store_step((*i).dat2_t_id[ivar,isp],data2)
    ENDFOR
  ENDFOR ; --- species loop

  IF (*i).tcorr THEN BEGIN
    (*i).dat1_t_id = store_step((*i).dat1_t_id,dat1_spec)
    (*i).t_arr_id = store_step((*i).t_arr_id,mom_time)
  ENDIF ELSE $
    (*i).corr_xy_id = time_avg((*i).corr_xy_id,corr_xy_all,mom_time)

END

;######################################################################

PRO corr_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nvar2 = N_ELEMENTS((*i).var2)

  IF NOT (*i).tcorr THEN BEGIN ; --- spatial correlation
    xaxis = (- par.nx0 / 2 + INDGEN(par.nx0)) * par.dx
    yaxis = (- par.ny0 / 2 + INDGEN(par.ny0)) * par.dy

    corr_xy = time_avg((*i).corr_xy_id,/avg)

    corrnormstr = '(!12<!9!!!6f!9!!!6!U2!N!12><!9!!!6g!9!!!6!U2!N!12>!6)!U1/2!N'
    IF N_ELEMENTS((*i).zind) EQ par.nz0 THEN zstr = 'z' $
      ELSE zstr = rm0es((*par.z)[(*i).zind])
    corrstr = '!8C!6(!7D!6x,!7D!6y) = !12<!6f(x,y,' + zstr + $
      ')g(x+!7D!6x,y+' + '!7D!6y,' + zstr + ')!12>!6 / ' + corrnormstr

    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      sp = (*gui.out.spec_select)[isp]

      ; spatial auto-correlations: correlation length
      lambda_xy = FLTARR(nvar2,2)
      FOR ivar = 0, nvar2 - 1 DO BEGIN
        IF (*i).var2[ivar] EQ (*i).var1 THEN BEGIN
          xind = 0
          x = par.nx0 / 2
          WHILE (xind EQ 0) AND (x LT par.nx0) DO BEGIN
            IF corr_xy[x,par.ny0/2,ivar,isp] LE 0.36788 THEN xind = x
            x += 1
          ENDWHILE
          yind = 0
          y = par.ny0 / 2
          WHILE (yind EQ 0) AND (y LT par.ny0) DO BEGIN
            IF corr_xy[par.nx0/2,y,ivar,isp] LE 0.36788 THEN yind = y
            y += 1
          ENDWHILE

          ; linear interpolation to find e-th part point
          IF xind NE 0 THEN BEGIN
            xind_eth = xind - par.nx0 / 2 - $
              (0.36788 - corr_xy[xind,par.ny0/2,ivar,isp]) / $
              (corr_xy[xind-1,par.ny0/2,ivar,isp] - $
              corr_xy[xind,par.ny0/2,ivar,isp])
          ENDIF ELSE xind_eth = !VALUES.F_INFINITY
          IF yind NE 0 THEN BEGIN
            yind_eth = yind - par.ny0 / 2 - $
              (0.36788 - corr_xy[par.nx0/2,yind,ivar,isp]) / $
              (corr_xy[par.nx0/2,yind-1,ivar,isp] - $
              corr_xy[par.nx0/2,yind,ivar,isp])
          ENDIF ELSE yind_eth = !VALUES.F_INFINITY

          PRINT, 'lambda_x,y = ' + rm0es(xind_eth*par.lx/par.nx0) + $
            ',' + rm0es(yind_eth*par.ly/par.ny0) + ' rho_ref'
          lambda_xy[ivar,*] = $
            [xind_eth*par.lx/par.nx0,yind_eth*par.ly/par.ny0]
        ENDIF
      ENDFOR

      set_output, diag, sp, /ps, multi=[0,1,2]

      maxy = MAX(corr_xy[*,par.ny0/2,*,isp],MIN=miny)
      PLOT, xaxis, corr_xy[*,par.ny0/2,0,isp], COLOR=1, CHARSIZE=1.41, $
        POSITION=[0.15,0.6,0.95,0.95], /XSTYLE, XTITLE='!7D!6x', $
        YTITLE='!8C!6(!7D!6x,0)', YRANGE=[miny,maxy]
      PLOTS, !X.CRANGE, [0,0], COLOR=1
      PLOTS, [0,0], !Y.CRANGE, COLOR=1
      OPLOT, [0,lambda_xy[0,0]], $
        !Y.CRANGE[1] * 0.36788 * [1,1], COLOR=1, LINE=1
      OPLOT, lambda_xy[0,0] * [1,1], $
        !Y.CRANGE[1] * 0.36788 * [0,1], COLOR=1, LINE=1
      FOR ivar = 1, nvar2 - 1 DO BEGIN
        OPLOT, xaxis, corr_xy[*,par.ny0/2,ivar,isp], COLOR=ivar+1
        OPLOT, [0,lambda_xy[ivar,0]], $
          !Y.CRANGE[1] * 0.36788 * [1,1], COLOR=ivar+1, LINE=1
        OPLOT, lambda_xy[ivar,0] * [1,1], $
          !Y.CRANGE[1] * 0.36788 * [0,1], COLOR=ivar+1, LINE=1
      ENDFOR

      maxy = MAX(corr_xy[par.nx0/2,*,*,isp],MIN=miny)
      PLOT, yaxis, corr_xy[par.nx0/2,*,0,isp], COLOR=1, CHARSIZE=1.41, $
        POSITION=[0.15,0.2,0.95,0.5], /XSTYLE, XTITLE='!7D!6y', $
        YTITLE='!8C!6(0,!7D!6y)', YRANGE=[miny,maxy]
      PLOTS, !X.CRANGE, [0,0], COLOR=1
      PLOTS, [0,0], !Y.CRANGE, COLOR=1
      OPLOT, [0,lambda_xy[0,1]], $
        !Y.CRANGE[1] * 0.36788 * [1,1], COLOR=1, LINE=1
      OPLOT, lambda_xy[0,1] * [1,1], $
        !Y.CRANGE[1] * 0.36788 * [0,1], COLOR=1, LINE=1
      FOR ivar = 1, nvar2 - 1 DO BEGIN
        OPLOT, yaxis, corr_xy[par.nx0/2,*,ivar,isp], COLOR=ivar+1
        OPLOT, [0,lambda_xy[ivar,1]], $
          !Y.CRANGE[1] * 0.36788 * [1,1], COLOR=ivar+1, LINE=1
        OPLOT, lambda_xy[ivar,1] * [1,1], $
          !Y.CRANGE[1] * 0.36788 * [0,1], COLOR=ivar+1, LINE=1
      ENDFOR

      XYOUTS, 0.1, 0.975, /NORMAL, corrstr, COLOR=1
      plot_info_str, diag

      IF (*i).var1 GE 0 THEN corr_title1 = get_var_string((*i).var1,/fancy)
      IF (*i).var1 EQ -1 THEN corr_title1 = '!6B!Dx!N'
      IF (*i).var1 EQ -2 THEN corr_title1 = '!6B!Dy!N'
      IF (*i).var1 EQ -3 THEN corr_title1 = '!6E!Dx!N'
      IF (*i).var1 EQ -4 THEN corr_title1 = '!6E!Dy!N'
      corr_title2 = STRARR(nvar2)
      FOR ivar = 0, nvar2 - 1 DO BEGIN
        IF (*i).var2[ivar] GE 0 THEN corr_title2[ivar] = $
          get_var_string((*i).var2[ivar],/fancy)
        IF (*i).var2[ivar] EQ -1 THEN corr_title2[ivar] = '!6B!Dx!N'
        IF (*i).var2[ivar] EQ -2 THEN corr_title2[ivar] = '!6B!Dy!N'
        IF (*i).var2[ivar] EQ -3 THEN corr_title2[ivar] = '!6E!Dx!N'
        IF (*i).var2[ivar] EQ -4 THEN corr_title2[ivar] = '!6E!Dy!N'
      ENDFOR

      plot_legend, INDGEN(nvar2) + 1, corr_title1 + ' x ' + corr_title2, $
        x_offset=[0.15,0.0,0.95,0.1], per_line=3

      LOADCT, 33, FILE='internal/colortable.tbl'
      FOR ivar = 0, nvar2 - 1 DO BEGIN
        corr_title = corr_title1 + ' x ' + corr_title2[ivar]

        lev_col = contour_levels(corr_xy[*,*,ivar,isp],20)
        CONTOUR, corr_xy[*,*,ivar,isp], xaxis, yaxis, LEVELS=lev_col[*,0], $
          C_COLORS=lev_col[*,1], /FILL, COLOR=1, /XSTYLE, /YSTYLE, $
          /ISOTROPIC, XTITLE='!7D!6x', YTITLE='!7D!6y', TITLE=corr_title
      ENDFOR

      FOR ivar = 0, nvar2 - 1 DO set_output, diag, sp, $
        dat=[[xaxis],[corr_xy[*,*,ivar,isp]]], append = (ivar GT 0), $
        header=['corr/y=',rm0es(yaxis)],$ 
        commentline=['corr_xy: '+get_var_string((*i).var1)+' x '+ $
        get_var_string((*i).var2[ivar])]

      set_output, diag, sp, /reset
    ENDFOR ; --species loop
  ENDIF ELSE BEGIN ; --- temporal correlation
    nz = N_ELEMENTS((*i).zind)

    dat1_t = store_step((*i).dat1_t_id,/get,/reg_array)
    dat1_t = REFORM(dat1_t,[par.nx0,par.ny0,nz,$
      gui.out.n_spec_sel,series.step_count],/OVERWRITE)
    t_arr = store_step((*i).t_arr_id,/get,/reg_array)
    corr_t = DBLARR(series.step_count,/NOZERO)

    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      sp = (*gui.out.spec_select)[isp]

      set_output, diag, sp, /ps

      FOR ivar = 0, nvar2 - 1 DO BEGIN
        dat2_t = store_step((*i).dat2_t_id[ivar,isp],/get,/reg_array)

        jac_rebin = REBIN((*series.geom).jac_norm[(*i).zind]/nz,$
          [nz,series.step_count])

        FOR t = 0, series.step_count - 1 DO corr_t[t] = $
          TOTAL(TOTAL(TOTAL(TOTAL(REFORM(dat1_t[*,*,*,isp,*])*$
          dat2_t[*,*,*,(INDGEN(series.step_count)+t) $
          MOD series.step_count],1),1)*jac_rebin,1),1)

        corr_norm = SQRT(TOTAL(TOTAL(TOTAL(TOTAL($
          REFORM(dat1_t[*,*,*,isp,*])^2,1),1)*jac_rebin,1),1)*$
          TOTAL(TOTAL(TOTAL(TOTAL(dat2_t^2,1),1)*jac_rebin,1),1))
        corr_t /= corr_norm

        IF ivar EQ 0 THEN PLOT, t_arr[*], corr_t[*], COLOR=1, CHARSIZE=1.41, $
          POSITION=[0.15,0.6,0.95,0.95], /XSTYLE, XTITLE='!7D!6t', $
          YTITLE='!8C!6(!7D!6t)', YRANGE=[0,1] ELSE $
          OPLOT, t_arr, corr_t, COLOR=ivar+1
        PLOTS, !X.CRANGE, [0,0], COLOR=1

        IF (*i).var2[ivar] EQ (*i).var1 THEN BEGIN
          tind = 0
          t = 0
          WHILE (tind EQ 0) AND (t LT series.step_count) DO BEGIN
            IF corr_t[t] LE 0.36788 THEN tind = t
            t += 1
          ENDWHILE
          IF tind NE 0 THEN BEGIN
            tind_eth = tind - (0.36788 - corr_t[tind]) / $
              (corr_t[tind-1] - corr_t[tind])
            dt = t_arr[tind] - t_arr[tind-1]
            PRINT, 'The correlation function has decayed to ' + $
            'the e-th part at t = ' + rm0es(tind_eth) + ', corresponding to ' + $
            rm0es(tind_eth*dt) + ' in time units.'

            OPLOT, [0,tind_eth*dt] + !X.CRANGE[0], $
              !Y.CRANGE[1] * 0.36788 * [1,1], COLOR=ivar+1, LINE=1
            OPLOT, tind_eth * dt + !X.CRANGE[0] * [1,1], $
              !Y.CRANGE[1] * 0.36788 * [0,1], COLOR=ivar+1, LINE=1
          ENDIF ELSE PRINT, $
            'correlation time calculation failed'
        ENDIF
      ENDFOR

      plot_info_str, diag

      IF (*i).var1 GE 0 THEN corr_title1 = get_var_string((*i).var1,/fancy)
      IF (*i).var1 EQ -1 THEN corr_title1 = '!6B!Dx!N'
      IF (*i).var1 EQ -2 THEN corr_title1 = '!6B!Dy!N'
      IF (*i).var1 EQ -3 THEN corr_title1 = '!6E!Dx!N'
      IF (*i).var1 EQ -4 THEN corr_title1 = '!6E!Dy!N'
      corr_title2 = STRARR(nvar2)
      FOR ivar = 0, nvar2 - 1 DO BEGIN
        IF (*i).var2[ivar] GE 0 THEN corr_title2[ivar] = $
          get_var_string((*i).var2[ivar],/fancy)
        IF (*i).var2[ivar] EQ -1 THEN corr_title2[ivar] = '!6B!Dx!N'
        IF (*i).var2[ivar] EQ -2 THEN corr_title2[ivar] = '!6B!Dy!N'
        IF (*i).var2[ivar] EQ -3 THEN corr_title2[ivar] = '!6E!Dx!N'
        IF (*i).var2[ivar] EQ -4 THEN corr_title2[ivar] = '!6E!Dy!N'
      ENDFOR

      plot_legend, INDGEN(nvar2) + 1, corr_title1 + ' x ' + corr_title2, $
        x_offset=[0.15,0.0,0.95,0.1], per_line=3

      set_output, diag, sp, dat=[t_arr[*],corr_t[*]], append=(isp NE 0)

      set_output, diag, sp, /reset
    ENDFOR
  ENDELSE

END
