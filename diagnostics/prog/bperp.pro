FUNCTION bperp_info

  RETURN, {$
    type      : 'mom',$
    title     : 'B_perp',$
    help_text : ['Shows time trace of perpendicular magnetic '+$
                 'field fluctuations.'],$
    ext_vars  : [['normalize','0','0: rms of B_perp (default), '+$
                  '1: rms of B_perp/B_0, 2: B_perp_rms/B_0_rms, '+$
                  '3: divide data by B_y,max(t)'],$
                 ['gamma_kx','0','indices of kx for growth '+$
                  'rate fit; default: off/avg. (if ky is specified; -1); '+$
                  'need to specify only positive index, negative '+$
                  'is added automatically'],$
                 ['gamma_ky','0','indices of ky for growth '+$
                  'rate fit; default: off/1st finite '+$
                  '(if kx is specified)'],$
                 ['add_Exy','1','also display the electric field'],$
                 ['add_Bpar','1','display B_parallel (if available)']]}

END

;######################################################################

PRO bperp_init, diag

  COMMON global_vars

  IF (par.n_fields LE 1) OR (par.nonlinear NE 1) THEN BEGIN
    PRINT, (*diag).name + ': need nonlinear electromagnetic runs'
    (*diag).selected = 0
    RETURN
  ENDIF

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.normalize) LT 1 THEN *e.normalize = 0
  gfit = 0
  kx_avg = 0
  IF (N_ELEMENTS(*e.gamma_kx) LT 1) AND (N_ELEMENTS(*e.gamma_ky) GE 1) $
    THEN *e.gamma_kx = -1
  IF NOT KEYWORD_SET(*e.add_Exy) THEN *e.add_Exy = 0
  *e.add_Bpar = KEYWORD_SET(*e.add_Bpar) AND (par.n_fields GE 2)

  IF N_ELEMENTS(*e.gamma_kx) GE 1 THEN BEGIN
    gfit = 1

    IF *e.gamma_kx[0] EQ -1 THEN BEGIN
      kx = INDGEN(par.nkx0)
      kx_avg = 1
    ENDIF ELSE BEGIN
      ; filter out higher numbers
      pos_kx = WHERE(*e.gamma_kx LE par.nkx0-par.nkx0/2-1)
      kx = pos_kx[0] NE -1 ? (*e.gamma_kx)[pos_kx] : 1 < (par.nkx0 - 1)
      ; add corresponding negative kx (not for zero mode)
      nozero_kx = WHERE(kx NE 0)
      IF nozero_kx[0] NE -1 THEN kx = [kx,par.nkx0-kx[nozero_kx]]
    ENDELSE
    IF N_ELEMENTS(*e.gamma_ky) LT 1 THEN ky = 1
  ENDIF
  IF N_ELEMENTS(*e.gamma_ky) GE 1 THEN BEGIN
    gfit = 1
    ky = *e.gamma_ky
    IF N_ELEMENTS(*e.gamma_kx) LT 1 THEN kx = 1 < (par.nkx0 - 1)
  ENDIF

  i = set_internal_vars(diag,{$
    norm         : *e.normalize,$
    gfit         : gfit,$
    kx_avg       : kx_avg,$
    kx           : gfit ? kx : 0,$
    ky           : gfit ? ky : 0,$
    add_Exy      : *e.add_Exy,$
    add_Bpar     : *e.add_Bpar,$
    B_0_inv      : *e.normalize EQ 1 ? 1.0/(*series.geom).Bfield : 0,$
    B_0_rms_inv  : *e.normalize EQ 2 ? $
                   1.0/SQRT(TOTAL((*series.geom).Bfield^2*$
                   (*series.geom).jac_norm)/par.nz0) : 0,$
    B_ymax_id    : PTR_NEW(),$
    time_id      : PTR_NEW(),$
    B_x_rms_id   : PTR_NEW(),$
    B_y_rms_id   : PTR_NEW(),$
    B_x_avg_id   : PTR_NEW(),$
    B_y_avg_id   : PTR_NEW(),$
    E_x_rms_id   : PTR_NEW(),$
    E_y_rms_id   : PTR_NEW(),$
    E_x_avg_id   : PTR_NEW(),$
    E_y_avg_id   : PTR_NEW(),$
    B_par_rms_id : PTR_NEW(),$
    B_par_avg_id : PTR_NEW()})

  fft_format, kxky=((*i).add_Exy ? [0,1] : 1)
  IF (*i).add_Bpar THEN fft_format, kxky=2

END

;######################################################################

PRO bperp_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  ; note: the sign is omitted since rms is taken later
  IF NOT (*i).gfit THEN BEGIN
    B_x_sqd = TOTAL(ABS(REBIN(REFORM((*series.ky),$
      [1,par.nky0,1]),[par.nkx0,par.nky0,par.nz0])*(*mom[0,1].kxky))^2,1)
    B_y_sqd = TOTAL(ABS(REBIN((*series.kx)[*,0],$
      [par.nkx0,par.nky0,par.nz0])*(*mom[0,1].kxky))^2,1)

    B_x_rms_z = 2.0 * TOTAL(B_x_sqd[1:par.nky0-1,*],1) + B_x_sqd[0,*]
    B_y_rms_z = 2.0 * TOTAL(B_y_sqd[1:par.nky0-1,*],1) + B_y_sqd[0,*]

    IF (*i).add_Exy THEN BEGIN
      E_x_sqd = TOTAL(ABS(REBIN((*series.kx)[*,0],$
        [par.nkx0,par.nky0,par.nz0])*(*mom[0,0].kxky))^2,1)
      E_y_sqd = TOTAL(ABS(REBIN(REFORM((*series.ky),$
        [1,par.nky0,1]),[par.nkx0,par.nky0,par.nz0])*(*mom[0,0].kxky))^2,1)

      E_x_rms_z = 2.0 * TOTAL(E_x_sqd[1:par.nky0-1,*],1) + E_x_sqd[0,*]
      E_y_rms_z = 2.0 * TOTAL(E_y_sqd[1:par.nky0-1,*],1) + E_y_sqd[0,*]
    ENDIF

    IF (*i).add_Bpar THEN BEGIN
      B_par_sqd = TOTAL(ABS(*mom[0,2].kxky)^2,1)

      B_par_rms_z = 2.0 * TOTAL(B_par_sqd[1:par.nky0-1,*],1) + B_par_sqd[0,*]
    ENDIF
  ENDIF ELSE BEGIN
    nkx = N_ELEMENTS((*i).kx)
    nky = N_ELEMENTS((*i).ky)

    B_x_rms_z = ABS(REBIN(REFORM((*series.ky)[(*i).ky],[1,nky,1]),$
      [nkx,nky,par.nz0])*(*mom[0,1].kxky)[(*i).kx,(*i).ky,*])^2
    B_y_rms_z = ABS(REBIN((*series.kx)[(*i).kx,0],$
      [nkx,nky,par.nz0])*(*mom[0,1].kxky)[(*i).kx,(*i).ky,*])^2

    IF (*i).add_Exy THEN BEGIN
      E_x_rms_z = ABS(REBIN((*series.kx)[(*i).kx,0],$
        [nkx,nky,par.nz0])*(*mom[0,0].kxky)[(*i).kx,(*i).ky,*])^2
      E_y_rms_z = ABS(REBIN(REFORM((*series.ky)[(*i).ky],[1,nky,1]),$
        [nkx,nky,par.nz0])*(*mom[0,0].kxky)[(*i).kx,(*i).ky,*])^2
    ENDIF

    IF (*i).add_Bpar THEN $
      B_par_rms_z = ABS((*mom[0,2].kxky)[(*i).kx,(*i).ky,*])^2
  ENDELSE

  IF (*i).norm EQ 1 THEN BEGIN
    B_x_rms_z *= (*i).gfit ? REBIN(REFORM((*i).B_0_inv,[1,1,par.nz0]),$
      [nkx,nky,par.nz0]) : (*i).B_0_inv
    B_y_rms_z *= (*i).gfit ? REBIN(REFORM((*i).B_0_inv,[1,1,par.nz0]),$
      [nkx,nky,par.nz0]) : (*i).B_0_inv

    IF (*i).add_Exy THEN BEGIN
      E_x_rms_z *= (*i).gfit ? REBIN(REFORM((*i).B_0_inv,[1,1,par.nz0]),$
        [nkx,nky,par.nz0]) : (*i).B_0_inv
      E_y_rms_z *= (*i).gfit ? REBIN(REFORM((*i).B_0_inv,[1,1,par.nz0]),$
        [nkx,nky,par.nz0]) : (*i).B_0_inv
    ENDIF

    IF (*i).add_Bpar THEN B_par_rms_z *= $
      (*i).gfit ? REBIN(REFORM((*i).B_0_inv,[1,1,par.nz0]),$
        [nkx,nky,par.nz0]) : (*i).B_0_inv
  ENDIF

  jac_norm = (*i).gfit ? REBIN(REFORM((*series.geom).jac_norm,$
    [1,1,par.nz0]),[nkx,nky,par.nz0]) : (*series.geom).jac_norm
  ; below, an extra REFORM is inserted for cases with Nz = 1
  B_x_rms = SQRT(TOTAL(REFORM(B_x_rms_z*jac_norm,$
    (*i).gfit ? [nkx,nky,par.nz0] : par.nz0),1+2*(*i).gfit)/par.nz0)
  B_y_rms = SQRT(TOTAL(REFORM(B_y_rms_z*jac_norm,$
    (*i).gfit ? [nkx,nky,par.nz0] : par.nz0),1+2*(*i).gfit)/par.nz0)

  IF (*i).add_Exy THEN BEGIN
    E_x_rms = SQRT(TOTAL(REFORM(E_x_rms_z*jac_norm,$
      (*i).gfit ? [nkx,nky,par.nz0] : par.nz0),1+2*(*i).gfit)/par.nz0)
    E_y_rms = SQRT(TOTAL(REFORM(E_y_rms_z*jac_norm,$
      (*i).gfit ? [nkx,nky,par.nz0] : par.nz0),1+2*(*i).gfit)/par.nz0)
  ENDIF

  IF (*i).add_Bpar THEN $
    B_par_rms = SQRT(TOTAL(REFORM(B_par_rms_z*jac_norm,$
      (*i).gfit ? [nkx,nky,par.nz0] : par.nz0),1+2*(*i).gfit)/par.nz0)

  IF (*i).norm EQ 2 THEN BEGIN
    B_x_rms *= (*i).B_0_rms_inv
    B_y_rms *= (*i).B_0_rms_inv

    IF (*i).add_Exy THEN BEGIN
      E_x_rms *= (*i).B_0_rms_inv
      E_y_rms *= (*i).B_0_rms_inv
    ENDIF

    IF (*i).add_Bpar THEN B_par_rms *= (*i).B_0_rms_inv
  ENDIF

  IF (*i).norm EQ 3 THEN BEGIN
    B_y_kxky = COMPLEXARR(par.nx0,par.nky0,par.nz0,/NOZERO)
    temp_data = COMPLEXARR(par.nx0,par.ny0,par.nz0,/NOZERO)
    B_y_kxky[0,0,0] = $
      - COMPLEX(0,1) * REBIN((*series.kx)[0:par.nkx0/2-1,0],$
      [par.nkx0/2,par.nky0,par.nz0]) * (*mom[0,1].kxky)[0:par.nkx0/2-1,*,*]
    B_y_kxky[par.nkx0/2,*,*] = 0
    B_y_kxky[par.nkx0/2+1,0,0] = $
      - COMPLEX(0,1) * REBIN((*series.kx)[par.nkx0/2+1:par.nkx0-1,0],$
      [par.nkx0-par.nkx0/2-1,par.nky0,par.nz0]) * $
      (*mom[0,1].kxky)[par.nkx0/2+1:par.nkx0-1,*,*]
    B_y_sxky = FFT(TEMPORARY(B_y_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)
    temp_data[0,0,0] = B_y_sxky[*,*,*]
    temp_data[*,par.nky0,*] = 0
    FOR y = 1, par.nky0 - 1 DO temp_data[0,par.ny0-y,0] = CONJ(B_y_sxky[*,y,*])
    B_y_sxky = 0
    B_y = FFT(TEMPORARY(temp_data),DIMENSION=2,/INVERSE,DOUBLE=0)

    B_ymax = MAX(ABS(FLOAT(TEMPORARY(B_y))))
    (*i).B_ymax_id = store_step((*i).B_ymax_id,B_ymax)
  ENDIF

  IF (*i).kx_avg THEN BEGIN
    B_x_rms = TOTAL(B_x_rms,1)
    B_y_rms = TOTAL(B_y_rms,1)

    IF (*i).add_Exy THEN BEGIN
      E_x_rms = TOTAL(E_x_rms,1)
      E_y_rms = TOTAL(E_y_rms,1)
    ENDIF

    IF (*i).add_Bpar THEN B_par_rms = TOTAL(B_x_rms,1)
  ENDIF

  (*i).time_id = store_step((*i).time_id,mom_time)
  (*i).B_x_rms_id = store_step((*i).B_x_rms_id,B_x_rms)
  (*i).B_y_rms_id = store_step((*i).B_y_rms_id,B_y_rms)
  (*i).B_x_avg_id = time_avg((*i).B_x_avg_id,B_x_rms,mom_time)
  (*i).B_y_avg_id = time_avg((*i).B_y_avg_id,B_y_rms,mom_time)
  IF (*i).add_Exy THEN BEGIN
    (*i).E_x_rms_id = store_step((*i).E_x_rms_id,E_x_rms)
    (*i).E_y_rms_id = store_step((*i).E_y_rms_id,E_y_rms)
    (*i).E_x_avg_id = time_avg((*i).E_x_avg_id,E_x_rms,mom_time)
    (*i).E_y_avg_id = time_avg((*i).E_y_avg_id,E_y_rms,mom_time)
  ENDIF
  IF (*i).add_Bpar THEN BEGIN
    (*i).B_par_rms_id = store_step((*i).B_par_rms_id,B_par_rms)
    (*i).B_par_avg_id = time_avg((*i).B_par_avg_id,B_par_rms,mom_time)
  ENDIF

END

;######################################################################

PRO bperp_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  set_output, diag, /ps

  time = store_step((*i).time_id,/get,/reg_array)
  B_x_rms = store_step((*i).B_x_rms_id,/get,/reg_array)
  B_y_rms = store_step((*i).B_y_rms_id,/get,/reg_array)
  B_x_avg = time_avg((*i).B_x_avg_id,/avg)
  B_y_avg = time_avg((*i).B_y_avg_id,/avg)
  IF (*i).add_Exy THEN BEGIN
    E_x_rms = store_step((*i).E_x_rms_id,/get,/reg_array)
    E_y_rms = store_step((*i).E_y_rms_id,/get,/reg_array)
    E_x_avg = time_avg((*i).E_x_avg_id,/avg)
    E_y_avg = time_avg((*i).E_y_avg_id,/avg)
  ENDIF
  IF (*i).add_Bpar THEN BEGIN
    B_par_rms = store_step((*i).B_par_rms_id,/get,/reg_array)
    B_par_avg = time_avg((*i).B_par_avg_id,/avg)
  ENDIF

  ; set axes and buffer strings to adjust the position of the y index
  time_str = get_var_string(0,/time,/fancy,/units)
  CASE (*i).norm OF
    1 : BEGIN
      norm_str_B = '<(B[x,y]/B_0)^2>^(1/2) = '
      norm_str_E = '<(E[x,y]/B_0)^2>^(1/2) = '
      norm_str_Bp = '<(Bpar/B_0)^2>^(1/2) = '
      norm_fstr_B = '!6<(B!Dx, !N/B!D0!N)!U2!N>!U1/2!N'
      norm_fstr_E = '!6<(E!Dx, !N/B!D0!N)!U2!N>!U1/2!N'
      norm_fstr_Bp = '!6<(B!D!9#!6!N/B!D0!N)!U2!N>!U1/2!N'
      buf_str = '     !I !N'
    END
    2 : BEGIN
      norm_str_B = '(<B[x,y]^2>/<B_0^2>)^(1/2) = '
      norm_str_E = '(<E[x,y]^2>/<B_0^2>)^(1/2) = '
      norm_str_Bp = '(<Bpar^2>/<B_0^2>)^(1/2) = '
      norm_fstr_B = '!6<B!S!Dx, !R!U2  !N>!U1/2!N/<B!S!D0!R!U2!N>!U1/2!N'
      norm_fstr_E = '!6<E!S!Dx, !R!U2  !N>!U1/2!N/<B!S!D0!R!U2!N>!U1/2!N'
      norm_fstr_Bp = '!6<B!S!D!9#!6!R!U2!N>!U1/2!N/<B!S!D0!R!U2!N>!U1/2!N'
      buf_str = '          '
    END
    3 : BEGIN
      norm_str_B = '<B[x,y]^2>^(1/2) = '
      norm_str_E = '<E[x,y]^2>^(1/2) = '
      norm_str_Bp = '<Bpar^2>^(1/2) = '
      norm_fstr_B = '!6<B!S!Dx, !R!U2  !N>!U1/2!N'
      norm_fstr_E = '!6<E!S!Dx, !R!U2  !N>!U1/2!N'
      norm_fstr_Bp = '!6<B!S!D!9#!6!R!U2!N>!U1/2!N'
      buf_str = ' '

      B_ymax = store_step((*i).B_ymax_id,/get,/reg_array)
      temp_time = time
      FOR t = 0, series.step_count - 1 DO temp_time[t] = $
        temp_time[t-1>0] + (time[t] - time[t-1>0]) * (B_ymax[t] / B_ymax[0])
      time = TEMPORARY(temp_time)

      time_str = '!S!9!A!EA!N!6!R' + time_str
      PRINT, (*diag).name + ': time axis adjusted step-wise by B_ymax[t]/B_ymax[0]'
      PRINT, '  with B_ymax[0] = ' + rm0es(B_ymax[0])
      PRINT, '  rel. B_ymax decay: ' + $
        rm0es((B_ymax[0]-B_ymax[series.step_count-1])/B_ymax[0],prec=5)
    END
    ELSE : BEGIN
      norm_str_B = '<B[x,y]^2>^(1/2) = '
      norm_str_E = '<E[x,y]^2>^(1/2) = '
      norm_str_Bp = '<Bpar^2>^(1/2) = '
      norm_fstr_B = '!6<B!S!Dx, !R!U2  !N>!U1/2!N'
      norm_fstr_E = '!6<E!S!Dx, !R!U2  !N>!U1/2!N'
      norm_fstr_Bp = '!6<B!S!D!9#!6!R!U2!N>!U1/2!N'
      buf_str = ' '
    END
  ENDCASE

  FOR j = 0, (*i).add_Exy + (*i).add_Bpar DO BEGIN
    IF (j EQ 1) AND NOT (*i).add_Exy AND (*i).add_Bpar THEN j = 2

    CASE j OF
      0    : BEGIN
               BE_x_rms = TEMPORARY(B_x_rms)
               BE_y_rms = TEMPORARY(B_y_rms)
               BE_x_avg = TEMPORARY(B_x_avg)
               BE_y_avg = TEMPORARY(B_y_avg)

               norm_str = norm_str_B
               norm_fstr = norm_fstr_B
             END
      1    : BEGIN
               BE_x_rms = TEMPORARY(E_x_rms)
               BE_y_rms = TEMPORARY(E_y_rms)
               BE_x_avg = TEMPORARY(E_x_avg)
               BE_y_avg = TEMPORARY(E_y_avg)

               norm_str = norm_str_E
               norm_fstr = norm_fstr_E
             END
      ELSE : BEGIN
               BE_x_rms = TEMPORARY(B_par_rms)
               BE_y_rms = BE_x_rms
               BE_x_avg = TEMPORARY(B_par_avg)
               BE_y_avg = BE_x_avg

               norm_str = norm_str_Bp
               norm_fstr = norm_fstr_Bp
             END
    ENDCASE

    IF NOT (*i).gfit THEN BEGIN
      PRINT, (*diag).name + ': ' + norm_str + rm0es(BE_x_avg,prec=4) + $
        (j EQ 2 ? '' : (',' + rm0es(BE_y_avg,prec=4)))

      yrange = [MIN(BE_x_rms)<MIN(BE_y_rms),MAX(BE_x_rms)>MAX(BE_y_rms)]

      PLOT, time, BE_x_rms, /XSTYLE, /YSTYLE, COLOR=1, YRANGE=yrange, $
        XTITLE=time_str, YTITLE=norm_fstr

      IF j NE 2 THEN BEGIN
        AXIS, YAXIS=0, /YSTYLE, COLOR=2, YTITLE='!6!D  y!N'+buf_str, $
          YRANGE=yrange

        PLOT, time, BE_x_rms, /XSTYLE, /YSTYLE, COLOR=1, YRANGE=yrange, $
          XTITLE=time_str, YTITLE=norm_fstr, $
          /NOERASE
        OPLOT, time, BE_y_rms, COLOR=2

        OPLOT, [time[0],time[N_ELEMENTS(time)-1]], [1,1] * BE_y_avg, $
          COLOR=2, LINE=1
      ENDIF

      OPLOT, [time[0],time[N_ELEMENTS(time)-1]], [1,1] * BE_x_avg, $
        COLOR=1, LINE=1
;      set_output, diag, dat=[[time],[BE_x_rms],[BE_y_rms]]
    ENDIF ELSE BEGIN
      nkx = (*i).kx_avg ? 1 : N_ELEMENTS((*i).kx)
      BE_x_rms = REFORM(BE_x_rms,[nkx,$
        N_ELEMENTS((*i).ky),series.step_count],/OVERWRITE)
      BE_y_rms = REFORM(BE_y_rms,[nkx,$
        N_ELEMENTS((*i).ky),series.step_count],/OVERWRITE)

      logBx = ALOG(BE_x_rms)
      logBy = ALOG(BE_y_rms)

      FOR y = 0, N_ELEMENTS((*i).ky) - 1 DO BEGIN
        FOR x = 0, nkx - 1 DO BEGIN
          lfx = LINFIT(time,logBx[x,y,*])
          lfy = LINFIT(time,logBy[x,y,*])

          kx_str = (*i).kx_avg ? 'avg.' : rm0es((*i).kx[x])
          PRINT, 'gamma(kx=' + kx_str + ',ky=' + $
            rm0es((*i).ky[y]) + ') for Bx, By: ' + $
            rm0es(lfx[1],prec=6) + ', ' + rm0es(lfy[1],prec=6)

          yrange = [MIN(BE_x_rms[x,y,*])<MIN(BE_y_rms[x,y,*]),$
            MAX(BE_x_rms[x,y,*])>MAX(BE_y_rms[x,y,*])]

          kx_str = (*i).kx_avg ? $
            'avg.' : rm0es((*series.kx)[(*i).kx[x],0],prec=5)

          PLOT, time, BE_x_rms[x,y,*], /XSTYLE, /YSTYLE, COLOR=1, $
            YRANGE=yrange, XTITLE=time_str, $
            YTITLE=norm_fstr, /YLOG, TITLE='!6k!Dx!N='+$
            kx_str+', k!Dy!N='+rm0es((*series.ky)[(*i).ky[y]],prec=5)

          IF j NE 2 THEN BEGIN
            AXIS, YAXIS=0, /YSTYLE, COLOR=2, YTITLE='!6!D  y!N'+buf_str, $
              YRANGE=yrange

            PLOT, time, BE_x_rms[x,y,*], /XSTYLE, /YSTYLE, COLOR=1, $
              YRANGE=yrange, XTITLE=time_str, $
              YTITLE=norm_fstr, /NOERASE, /YLOG
            OPLOT, time, BE_y_rms[x,y,*], COLOR=2

            OPLOT, time, EXP(lfy[0]+time*lfy[1]), COLOR=2, LINE=2
            OPLOT, [time[0],time[N_ELEMENTS(time)-1]], $
              [1,1] * BE_y_avg[x,y,*], COLOR=2, LINE=1
          ENDIF

          OPLOT, time, EXP(lfx[0]+time*lfx[1]), COLOR=1, LINE=2
          OPLOT, [time[0],time[N_ELEMENTS(time)-1]], $
            [1,1] * BE_x_avg[x,y,*], COLOR=1, LINE=1
        ENDFOR
      ENDFOR
    ENDELSE

    plot_info_str, diag
  ENDFOR

  set_output, diag, /reset

END
