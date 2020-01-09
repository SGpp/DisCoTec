FUNCTION frequency_info

  RETURN, {$
    type      : 'mom', $;'mom_uni',$
    title     : 'Frequencies',$
    help_text : ['In standard mode, frequency and growth rate are '+$
                 'evaluated using dphase = log(phi(t)/phi(t-dt)). '+$
                 'Error bars denote one standard deviation. '+$
                 'In FFT mode (e.g., for nonlinear runs), a direct '+$
                 'as well as a windowed FFT are used to extract '+$
                 'dominant/subdominant frequencies. The results are '+$
                 'plotted as a function of ky.'],$
    ext_vars  : [['kxind','0','x modes to be analyzed; default: -1 '+$
                  '(all connected); -2 for all (may yield smoother '+$
                  'curves); -3 for x-resolved plotting of all modes'],$
                 ['kyind','0','y modes to be plotted; default: -1 '+$
                  '(all non-zero)'],$
                 ['Gaussian','0','use with FFT mode; set to scalar '+$
                  'for plotting Gaussian; set to [min,max] to fit '+$
                  'Gaussian in that frequency range; automatically '+$
                  'sets xrange, ylog; default: off'],$
                 ['fftmode','1','Fourier transform of the time '+$
                  'coordinate instead of phase evaluation'],$
                 ['weight','1','weight by the absolute of the '+$
                  'amplitude (ignored in FFT mode)']]}

END

;######################################################################

PRO frequency_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.kxind) LT 1 THEN *e.kxind = -1
  IF N_ELEMENTS(*e.kyind) LT 1 THEN *e.kyind = -1
  IF NOT KEYWORD_SET(*e.fftmode) THEN *e.fftmode = 0
  IF NOT KEYWORD_SET(*e.weight) THEN *e.weight = 0

  gmode = N_ELEMENTS(*e.Gaussian)
  IF (gmode LT 1) OR (gmode GT 2) THEN gmode = 0
  gwidth = gmode EQ 1 ? *e.Gaussian : 0.0
  grange = gmode EQ 2 ? *e.Gaussian : [0.0,0.0]

  IF (*e.kyind)[0] LT 0 THEN *e.kyind = $
    (par.ky0_ind > 1) + INDGEN(par.nky0-(par.ky0_ind>1))

  resolve_kx = 0
  IF (*e.kxind)[0] EQ -3 THEN BEGIN
    IF NOT *e.fftmode THEN BEGIN
      *e.kxind = par.x_local ? [par.nkx0/2+1+INDGEN(par.nkx0/2-1),$
        INDGEN(par.nkx0/2)] : INDGEN(par.nkx0)
      resolve_kx = 1
    ENDIF ELSE *e.kxind = -1
  ENDIF

  manual_kx = (WHERE(*e.kxind LT 0))[0] EQ -1

  connect = 0
  IF ((*e.kxind)[0] EQ -1) AND (N_ELEMENTS(*e.kyind) GT 1) THEN BEGIN
    IF par.adapt_lx THEN BEGIN
      PRINT, (*diag).name + ': cannot plot kx/ky contours if ' + $
        'adapt_lx is set; averaging over connected modes'
      *e.kxind = -1
    ENDIF ELSE IF par.x_local THEN connect = 1
  ENDIF

  IF ((*e.kxind)[0] EQ -1) OR ((*e.kxind)[0] EQ -2) THEN $
    *e.kxind = INDGEN(par.nkx0)

  i = set_internal_vars(diag,{$
    connect    : connect,$
    gmode      : gmode,$
    gwidth     : gwidth,$
    grange     : grange,$
    resolve_kx : resolve_kx,$
    manual_kx  : manual_kx,$
    kxind      : *e.kxind,$
    kyind      : *e.kyind,$
    fftmode    : *e.fftmode,$
    weight     : *e.weight,$
    time_id    : PTR_NEW(),$
    data_id    : PTR_NEW()})
    
  IF par.x_local THEN fft_format, kxky=[0] ELSE fft_format, sxky=[0]

END    

;######################################################################

PRO frequency_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  data = par.x_local ? (*mom[0,0].kxky)[(*i).kxind,(*i).kyind,*] : $
    (*mom[0,0].sxky)[(*i).kxind,(*i).kyind,*]
  data = COMPLEX(TEMPORARY(data))

  (*i).time_id = store_step((*i).time_id,FLOAT(mom_time))
  (*i).data_id = store_step((*i).data_id,data)

END

;######################################################################

PRO frequency_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  with_splines = 1

  IF series.step_count LE 2 THEN BEGIN
    printerror, (*diag).name + ': < 3 time points is insufficient ' + $
      'for calculation of time derivatives'
    PTR_FREE, (*i).time_id, (*i).data_id
    RETURN
  ENDIF

  nkxind = N_ELEMENTS((*i).kxind)
  nkyind = N_ELEMENTS((*i).kyind)

  data = store_step((*i).data_id,/get,/reg_array)
  data = REFORM(data,[nkxind,nkyind,par.nz0,series.step_count],$
    /OVERWRITE)
  time = store_step((*i).time_id,/get,/reg_array)

  IF (*i).connect THEN BEGIN
    IF par.n_conn GT 0 THEN bigN = par.n_conn ELSE $
      IF ABS(par.shat) GT 1e-3 THEN bigN = ROUND(ABS(par.shat)*$
            par.lx*par.kymin*(*series.geom).C_y*par.q0/par.x0) $
      ELSE bigN = 1
  ENDIF

  IF NOT (*i).fftmode THEN BEGIN ; --- phase mode ---
    error_bars = 1
    kylog = 0
    xmargin = [12,3]

    n_steps = series.step_count - 1

    IF (*i).resolve_kx THEN BEGIN
      dphase_dt = COMPLEXARR(nkxind,n_steps,/NOZERO)
      dphase_dt = REFORM(dphase_dt,[nkxind,n_steps],/OVERWRITE)
      omega_avg = FLTARR(nkxind,nkyind,/NOZERO)
      gamma_avg = FLTARR(nkxind,nkyind,/NOZERO)
    ENDIF ELSE BEGIN
      dphase_dt = COMPLEXARR(n_steps,/NOZERO)
      omega = FLTARR(n_steps,nkyind,/NOZERO)
      gamma = FLTARR(n_steps,nkyind,/NOZERO)
      omega_avg = FLTARR(nkyind,/NOZERO)
      gamma_avg = FLTARR(nkyind,/NOZERO)
      IF error_bars THEN BEGIN
        omega_err = FLTARR(nkyind,/NOZERO)
        gamma_err = FLTARR(nkyind,/NOZERO)
      ENDIF
    ENDELSE

    FOR y = 0, nkyind - 1 DO BEGIN
      IF (*i).connect THEN BEGIN
        kxind_pos = INDGEN((par.nkx0-par.nkx0/2-1)/$
          (ABS(bigN)*(*i).kyind[y])+1) * (ABS(bigN) * (*i).kyind[y])
        IF N_ELEMENTS(kxind_pos) GE 2 THEN BEGIN
          kxind_neg = par.nkx0 - kxind_pos[1:*]
          kxind = [kxind_pos,REVERSE(kxind_neg)]
        ENDIF ELSE kxind = kxind_pos
        nkxind = N_ELEMENTS(kxind)
      ENDIF ELSE IF NOT (*i).resolve_kx THEN BEGIN
        kxind = (*i).kxind
        IF NOT (par.nkx0 MOD 2) AND $
          ((WHERE(kxind EQ par.nkx0/2))[0] NE -1) THEN $
          kxind = [kxind[0:par.nkx0/2-1],kxind[par.nkx0/2+1:*]]
        nkxind = N_ELEMENTS(kxind)
      ENDIF ELSE kxind = par.x_local ? $
        SHIFT((*i).kxind,par.nkx0/2-1) : (*i).kxind

      IF (*i).manual_kx THEN kxind = INDGEN(nkxind)

      data_ky = REFORM(data[kxind,y,*,*],[nkxind,par.nz0,n_steps+1])

      IF (*i).weight THEN BEGIN
        data_old = data_ky[*,*,0:n_steps-1]
        IF (*i).resolve_kx THEN BEGIN
          dphase_dt[0,0] = TOTAL(ALOG(data_ky[*,*,1:*]/data_old)*$
            ABS(data_old),2) / TOTAL(ABS(data_old),2)
          dphase_dt /= REBIN(REFORM(time[1:*]-time[0:n_steps-1],$
            [1,n_steps]),[nkxind,n_steps])
        ENDIF ELSE BEGIN
          dphase_dt[0] = TOTAL(TOTAL(ALOG(data_ky[*,*,1:*]/data_old)*$
            ABS(data_old),1),1) / TOTAL(TOTAL(ABS(data_old),1),1)
          dphase_dt /= time[1:*] - time[0:n_steps-1]
        ENDELSE
      ENDIF ELSE BEGIN
        IF par.x_local THEN BEGIN
          dphase_dt[0,0] = (*i).resolve_kx ? TOTAL(ALOG($
            data_ky[*,*,1:*]/data_ky[*,*,0:n_steps-1])*$
            REBIN(REFORM((*series.geom).jac_norm/par.nz0,$
            [1,par.nz0,1]),[nkxind,par.nz0,n_steps]),2,/NAN) / $
            REBIN(REFORM((time[1:*]-time[0:n_steps-1]),$
            [1,n_steps]),[nkxind,n_steps]) : $
            TOTAL(TOTAL(ALOG(data_ky[*,*,1:*]/data_ky[*,*,0:n_steps-1]),1)*$
            REBIN((*series.geom).jac_norm/par.nz0,[par.nz0,n_steps]),$
            1,/NAN) / ((time[1:*] - time[0:n_steps-1]) * nkxind)
        ENDIF ELSE BEGIN
          dphase_dt[0,0] = (*i).resolve_kx ? TOTAL(ALOG($
            data_ky[*,*,1:*]/data_ky[*,*,0:n_steps-1])*REBIN($
            (*series.geom).jac_norm/par.nz0,[nkxind,par.nz0,n_steps]),$
            2,/NAN) / REBIN(REFORM((time[1:*]-time[0:n_steps-1]),$
            [1,n_steps]),[nkxind,n_steps]) : $
            TOTAL(TOTAL(ALOG(data_ky[*,*,1:*]/data_ky[*,*,0:n_steps-1])*$
            REBIN((*series.geom).jac_norm[kxind,*]/par.nz0,$
            [nkxind,par.nz0,n_steps]),1),1,/NAN) / $
            ((time[1:*] - time[0:n_steps-1]) * nkxind)
        ENDELSE
      ENDELSE ; --- no weight

      IF (*i).resolve_kx THEN BEGIN
        omega_avg[0,y] = MEDIAN(IMAGINARY(dphase_dt),DIMENSION=2)
        gamma_avg[0,y] = MEDIAN(FLOAT(dphase_dt),DIMENSION=2)
      ENDIF ELSE BEGIN
        omega[0,y] = IMAGINARY(dphase_dt)
        gamma[0,y] = FLOAT(dphase_dt)

        ; skip amplitude rescalings in linear runs for gamma computation
        resc_inds = $
          WHERE(gamma[*,y]*(time[1:*]-time[0:n_steps-1]) LT -10,resc_count)
        IF resc_count GT 0 THEN BEGIN
          PRINT, (*diag).name + ': skipped rescaled amplitudes for gamma'
          gamma[resc_inds,y] = !VALUES.F_NAN
        ENDIF

        omega_avg[y] = MEDIAN(omega[*,y])
        gamma_avg[y] = MEDIAN(gamma[*,y])
        IF error_bars THEN BEGIN
          omega_err[y] = STDDEV(omega[*,y]) / SQRT(par.nz0*nkxind)
          gamma_err[y] = STDDEV(gamma[*,y],/NAN) / SQRT(par.nz0*nkxind)
        ENDIF
      ENDELSE
    ENDFOR

    ytitle = '!7x!6 / (' + get_var_string(-1,/time,/ounit) + ')'
    dt_max = MAX(time[1:*]-time[0:n_steps-1])
    limit = !PI / dt_max

    set_output, diag, /ps, multi=[0,1,2], charsize=1.2, ysize=22

    ; --- plot frequency time trace ---
    IF NOT (*i).resolve_kx THEN BEGIN
      omax = MAX(omega,MIN=omin)
      PLOT, time[[0,n_steps-1]], [0,1], COLOR=1, /XSTYLE, $
        YRANGE=[omin,omax], XTITLE=get_var_string(0,/time,/units), $
        YTITLE=ytitle, TITLE=title, /NODATA, XMARGIN=xmargin
      FOR y = 0, nkyind - 1 DO $
        OPLOT, time[0:n_steps-1], omega[*,y], COLOR=1+y

      OPLOT, !X.CRANGE, [0,0], COLOR=1
      OPLOT, !X.CRANGE, [1,1] * limit, COLOR=2, THICK=6
      OPLOT, !X.CRANGE, [-1,-1] * limit, COLOR=2, THICK=6
    ENDIF

    ytitle = '!S!7x!6!R!A-!N / (' + $
      get_var_string(-1,/time,/ounit) + ')'
    IF (*i).resolve_kx THEN BEGIN
      IF par.x_local THEN BEGIN
        xaxis = (*series.kx)[0:2*par.nkx0-par.nkx0/2*2-2,(*i).kyind[0]]
        xtitle = '!6k!Dx!N ' + get_var_string(/rhostr)
      ENDIF ELSE BEGIN
        xaxis = INDGEN(par.nkx0) * series.dx - 0.5 * series.lx
        xtitle = 'x / ' + get_var_string(/rhostr)
      ENDELSE

      xaxis = xaxis[SORT(xaxis)]
      omega_avg = omega_avg[SORT(xaxis),*]

      IF nkyind EQ 1 THEN PLOT, xaxis, omega_avg, COLOR=1, $
        /XSTYLE, XTITLE=xtitle, TITLE=title ELSE BEGIN

        LOADCT, 33, FILE='internal/colortable.tbl'
        lev_col = contour_levels(omega_avg,20)
        CONTOUR, omega_avg, xaxis, (*series.ky)[(*i).kyind],$
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
          /XSTYLE, /YSTYLE, XTITLE=xtitle, TITLE=ytitle, $
          YTITLE='!6k!Dy!N'+get_var_string(/rhostr), YLOG=kylog, $
          POSITION=[0.1,0.6,0.885,0.95], XMARGIN=xmargin
        plot_colorbar, lev_col, POSITION=[0.9,0.6,0.925,0.95], prec=2
        LOADCT, 41, FILE='internal/colortable.tbl'
      ENDELSE
    ENDIF ELSE BEGIN
      IF nkyind GT 1 THEN BEGIN
        PLOT, (*series.ky)[(*i).kyind], omega_avg, COLOR=1, /XSTYLE, $
          XTITLE='!6k!Dy!N'+get_var_string(/rhostr), YTITLE=ytitle, $
          TITLE=title, XMARGIN=xmargin, XLOG=kylog
        FOR y = 0, nkyind - 1 DO OPLOT, [1,1] * $
          (*series.ky)[(*i).kyind[y]], [1,1] * omega_avg[y], $
          COLOR=1+y, PSYM=8
        IF error_bars THEN ERRPLOT, (*series.ky)[(*i).kyind], $
          omega_avg - omega_err, omega_avg + omega_err, COLOR=1
      ENDIF

      xaxis = kylog AND (nkyind GT 1) ? 10^!X.CRANGE : !X.CRANGE
      OPLOT, xaxis, [0,0], COLOR=1
      OPLOT, xaxis, [1,1] * limit, COLOR=2, THICK=6
      OPLOT, xaxis, [-1,-1] * limit, COLOR=2, THICK=6
    ENDELSE

    plot_info_str, diag

    ; --- plot growth rate time trace ---
    ytitle = '!7c!6 / (' + get_var_string(-1,/time,/ounit) + ')'

    IF NOT (*i).resolve_kx THEN BEGIN
      gmax = MAX(gamma,MIN=gmin)
      PLOT, time[[0,n_steps-1]], [0,1], COLOR=1, /XSTYLE, $
        YRANGE=[gmin,gmax], XTITLE=get_var_string(0,/time,/units), $
        YTITLE=ytitle, TITLE=title, /NODATA, XMARGIN=xmargin
      FOR y = 0, nkyind - 1 DO $
        OPLOT, time[0:n_steps-1], gamma[*,y], COLOR=1+y
  
      OPLOT, !X.CRANGE, [0,0], COLOR=1
      OPLOT, !X.CRANGE, [1,1] * limit, COLOR=2, THICK=6
      OPLOT, !X.CRANGE, [-1,-1] * limit, COLOR=2, THICK=6
    ENDIF

    ytitle = '!S!7c!6!R!A-!N / (' + $
      get_var_string(-1,/time,/ounit) + ')'
    IF (*i).resolve_kx THEN BEGIN
      IF par.x_local THEN BEGIN
        xaxis = (*series.kx)[0:2*par.nkx0-par.nkx0/2*2-2,(*i).kyind[0]]
        xtitle = '!6k!Dx!N ' + get_var_string(/rhostr)
      ENDIF ELSE BEGIN
        xaxis = INDGEN(par.nkx0) * series.dx - 0.5 * series.lx
        xtitle = 'x / ' + get_var_string(/rhostr)
      ENDELSE

      xaxis = xaxis[SORT(xaxis)]
      gamma_avg = gamma_avg[SORT(xaxis),*]

      IF nkyind EQ 1 THEN PLOT, xaxis, gamma_avg, COLOR=1, $
        /XSTYLE, XTITLE=xtitle, TITLE=title ELSE BEGIN

        LOADCT, 33, FILE='internal/colortable.tbl'
        lev_col = contour_levels(gamma_avg,20)
        CONTOUR, gamma_avg, xaxis, (*series.ky)[(*i).kyind],$
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
          /XSTYLE, /YSTYLE, XTITLE=xtitle, TITLE=ytitle, $
          YTITLE='!6k!Dy!N'+get_var_string(/rhostr), YLOG=kylog, $
          POSITION=[0.1,0.1,0.885,0.45], XMARGIN=xmargin
        plot_colorbar, lev_col, POSITION=[0.9,0.1,0.925,0.45], prec=2
        LOADCT, 41, FILE='internal/colortable.tbl'
      ENDELSE
    ENDIF ELSE BEGIN
      IF nkyind GT 1 THEN BEGIN
        PLOT, (*series.ky)[(*i).kyind], gamma_avg, COLOR=1, /XSTYLE, $
          XTITLE='!6k!Dy!N'+get_var_string(/rhostr), YTITLE=ytitle, $
          TITLE=title, XMARGIN=xmargin, XLOG=kylog
        FOR y = 0, nkyind - 1 DO OPLOT, [1,1] * $
          (*series.ky)[(*i).kyind[y]], [1,1] * gamma_avg[y], $
          COLOR=1+y, PSYM=8
        IF error_bars THEN ERRPLOT, (*series.ky)[(*i).kyind], $
          gamma_avg - gamma_err, gamma_avg + gamma_err, COLOR=1
      ENDIF

      xaxis = kylog AND (nkyind GT 1) ? 10^!X.CRANGE : !X.CRANGE
      OPLOT, xaxis, [0,0], COLOR=1
      OPLOT, xaxis, [1,1] * limit, COLOR=2, THICK=6
      OPLOT, xaxis, [-1,-1] * limit, COLOR=2, THICK=6
    ENDELSE

    plot_info_str, diag

    ; --- print average values, data output ---
    commentline = 'Nyquist frequency: ' + rm0es(limit)
    IF (*i).weight THEN commentline += ', amplitude weighed'

    IF NOT (*i).resolve_kx THEN BEGIN
      PRINT, (*diag).name + ': averages (gamma,omega):'
      PRINT, STRING('ky',FORMAT='(A-7)')+$
             STRING('gamma','omega',FORMAT='(2(1X,A9))')
      FOR y = 0, nkyind - 1 DO PRINT, $
        STRING((*series.ky)[(*i).kyind[y]],FORMAT='(F7.3)') + $
        STRING(gamma_avg[y],omega_avg[y],FORMAT='(2(1X,F9.3))')

      set_output, diag, sp, header=['ky','gamma','omega'], $
        commentline=commentline, dat=[[(*series.ky)[(*i).kyind]],$
        [gamma_avg],[omega_avg]]
    ENDIF ELSE set_output, diag, sp, header=['ky','gamma','omega'], $
      commentline=commentline, dat=[[(*series.ky)[(*i).kyind]],$
      [TOTAL(gamma_avg,1)],[TOTAL(omega_avg,1)]]
  ENDIF ELSE BEGIN ; --- FFT mode ---
    n_freqs = 2 ; maximum number of freqs. (dominant+subdominant)
    twin = [0.25,0.5] ; windowed FFT parameters

    ytitle = get_var_string(0,/fancy) + '(k!Dx!N,k!Dy!N,t)'

    ; prepare data interpolation for equidistant grid
    time_range = time[series.step_count-1] - time[0]
    ; choose an interpolation factor to achieve high FFT speeds
    n_steps_new = ROUND(1.5*series.step_count)
    base2 = 2L
    WHILE n_steps_new GT base2 DO base2 *= 2
    n_steps_new = base2
    time_new = time[0] + INDGEN(n_steps_new) * $
      (time_range / (n_steps_new - 1.0))
    n_steps_win = FIX(n_steps_new*twin[0])
    window_range = time_new[n_steps_win-1] - time_new[0]
    ; time array for irregular-to-regular interpolation
    equid_time = INTERPOL(FINDGEN(N_ELEMENTS(time)),time,time_new)

    omega = 2 * !PI / time_range * $
      (- n_steps_new / 2 + 1 + INDGEN(n_steps_new))
    omega_win = 2 * !PI / window_range * $
      (- n_steps_win / 2 + 1 + INDGEN(n_steps_win))

    PRINT, (*diag).name + ': min/max omega for direct FFT: ' + $
      rm0es(omega[n_steps_new/2]) + ', ' + $
      rm0es(omega[n_steps_new-1])
    PRINT, (*diag).name + ': min/max omega for windowed FFT: ' + $
      rm0es(omega_win[n_steps_win/2]) + ', ' + $
      rm0es(omega_win[n_steps_win-1])

    set_output, diag, /ps, multi=[0,1,1], charsize=1.4, ysize=10

    omega_max = FLTARR(2,n_freqs,nkyind)
    omega_moms = FLTARR(nkyind)

    fftdata_win = FLTARR(nkyind,n_steps_win)

    FOR iky = 0, nkyind - 1 DO BEGIN
      IF (*i).connect THEN BEGIN
        kxind_pos = INDGEN((par.nkx0-par.nkx0/2-1)/$
          (ABS(bigN)*(*i).kyind[iky])+1) * (ABS(bigN) * (*i).kyind[iky])
        IF N_ELEMENTS(kxind_pos) GE 2 THEN BEGIN
          kxind_neg = par.nkx0 - kxind_pos[1:*]
          kxind = [kxind_pos,REVERSE(kxind_neg)]
        ENDIF ELSE kxind = kxind_pos
        nkxind = N_ELEMENTS(kxind)

        data_ip = INTERPOLATE(REFORM(data[kxind,iky,*,*],$
          [nkxind,par.nz0,series.step_count]),$
          INDGEN(nkxind),INDGEN(par.nz0),equid_time,/GRID)
      ENDIF ELSE BEGIN
        data_ip = INTERPOLATE(REFORM(data[*,iky,*,*],$
          [nkxind,par.nz0,series.step_count]),$
          INDGEN(nkxind),INDGEN(par.nz0),equid_time,/GRID)
      ENDELSE

      fftdata_temp = par.x_local ? $
        TOTAL(ABS(FFT(data_ip,DIMENSION=3))^2,1) * $
        REBIN((*series.geom).jac_norm,[par.nz0,n_steps_new]) : $
        TOTAL(ABS(FFT(data_ip,DIMENSION=3))^2*$
        REBIN((*series.geom).jac_norm[(*i).kxind,*],$
        [nkxind,par.nz0,n_steps_new]),1)
      fftdata = SHIFT(TOTAL(TEMPORARY(fftdata_temp),1),n_steps_new/2-1)
      fftdata /= nkxind * par.nz0

      ; loop over twindows
      step_start_win = 0
      step_end_win = FIX(n_steps_new*twin[0]) - 1
      n_windows = 0
      IF step_end_win LE n_steps_new - 1 THEN win_hanning = $
        REBIN(REFORM(HANNING(n_steps_win),$
        [1,1,n_steps_win]),[nkxind,par.nz0,n_steps_win])
      WHILE step_end_win LE n_steps_new - 1 DO BEGIN
        data_filtered = REFORM(win_hanning*data_ip[*,*,$
          step_start_win:step_end_win],[nkxind,par.nz0,n_steps_win])
        fftdata_win_temp = par.x_local ? $
          TOTAL(ABS(FFT(data_filtered,DIMENSION=3))^2,1) * $
          REBIN((*series.geom).jac_norm,[par.nz0,n_steps_win]) : $
          TOTAL(ABS(FFT(data_filtered,DIMENSION=3))^2*$
          REBIN((*series.geom).jac_norm[(*i).kxind,*],$
          [nkxind,par.nz0,n_steps_win]),1)
        fftdata_win[iky,*] += $
          SHIFT(TOTAL(TEMPORARY(fftdata_win_temp),1),n_steps_win/2-1)

        step_start_win += FIX(n_steps_win*(1.0-twin[1]))
        step_end_win += FIX(n_steps_win*(1.0-twin[1]))
        n_windows += 1
      ENDWHILE

      fftdata_win[iky,*] /= nkxind * par.nz0 * n_windows

      title = 'k!Dy!N ' + get_var_string(/rhostr) + ' = ' + $
        rm0es((*series.ky)[(*i).kyind[iky]],prec=3)

      IF (with_splines) THEN BEGIN
         SPLINE_P, omega_win, fftdata_win[iky,*], omega_spline, fftdata_spline

         ; search for n_freqs local maxima
         n_steps_spline = N_ELEMENTS(fftdata_spline)
         fftdata_spline_dt1 = fftdata_spline[1:n_steps_spline-2] - $
                              fftdata_spline[0:n_steps_spline-3]
         fftdata_spline_dt2 = fftdata_spline[1:n_steps_spline-2] - $
                              fftdata_spline[2:n_steps_spline-1]
         max_loc = WHERE((fftdata_spline_dt1 GT 0) AND $
                         (fftdata_spline_dt2 GT 0)) + 1
         fftdata_spline_max = fftdata_spline[max_loc]
         omega_spline_max = omega_spline[max_loc]
         ; filter out irrelevant maxima
         max_loc = WHERE(fftdata_spline_max GE 0.05*MAX(fftdata_spline_max))
         fftdata_spline_max = fftdata_spline_max[max_loc]
         omega_spline_max = omega_spline_max[max_loc]
         ; sort by size
         max_loc = REVERSE(SORT(fftdata_spline_max))
         fftdata_spline_max = fftdata_spline_max[max_loc]
         omega_spline_max = omega_spline_max[max_loc]
         n_max = N_ELEMENTS(fftdata_spline_max) < n_freqs
         omega_max[0,0:n_max-1,iky] = omega_spline_max[0:n_max-1]
         omega_max[1,0:n_max-1,iky] = fftdata_spline_max[0:n_max-1]

         IF !QUIET NE 1 THEN FOR j = 0, n_freqs - 1 DO PRINT, $
            '   ' + rm0es(j+1) + '. frequency (windowed FFT, ky = ' + $
            rm0es((*series.ky)[(*i).kyind[iky]],prec=3) + '): ' + $
            rm0es(omega_max[0,j,iky])
         ; calculate first moment
         omega_moms[iky] = INT_TABULATED(omega_spline,omega_spline*$
            fftdata_spline) / INT_TABULATED(omega_spline,fftdata_spline)
         xrange = 3.0 * MAX(ABS(omega_max[0,*,iky])) * [-1,1]
      ENDIF ELSE BEGIN
         xr_max = 3
         xrange = xr_max*[-1,1]
      ENDELSE

      IF (*i).gmode EQ 1 THEN BEGIN
        omega_axis = (FINDGEN(200) / 100.0 - 1.0) * 10.0
        fftmax = MAX(fftdata,maxind)
        omega_fftmax = omega[maxind]
        gaussfunc = EXP(-((omega_axis-omega_fftmax)/(*i).gwidth)^2) * fftmax

        xrange = [-10.0,10.0]
      ENDIF
      IF (*i).gmode EQ 2 THEN BEGIN
        ; reduce x axis to user specified range
        min_ind = (WHERE(omega GE (*i).grange[0]))[0]
        IF min_ind EQ -1 THEN min_ind = 0
        max_ind = WHERE(omega LE (*i).grange[1],maxcount)
        max_ind = (maxcount EQ -1) OR (maxcount EQ 0) ? 0 : max_ind[maxcount-1]
        omega_axis = omega[min_ind:max_ind]

        gaussfunc = GAUSSFIT(omega_axis,fftdata[min_ind:max_ind],$
          coeffs,NTERMS=3)
        PRINT, 'width of Gaussian (ky=' + $
          rm0es((*series.ky)[(*i).kyind[iky]],prec=3) + $
          '): ' + rm0es(coeffs[2])
        PRINT, 'full width at half maximum (ky=' + $
          rm0es((*series.ky)[(*i).kyind[iky]],prec=3) + $
          '): ' + rm0es(2*SQRT(2*ALOG(2))*coeffs[2])

        xrange = (*i).grange
      ENDIF

      PLOT, omega, fftdata, COLOR=1, XTITLE='!7x!6 / ('+$
        get_var_string(-1,/time,/ounit)+')', YTITLE='!9!!!6'+ytitle+$
        '!9!!!6!U2!N', /XSTYLE, /NODATA, POSITION=[0.15,0.4,0.9,0.9], $
        TITLE=title, XRANGE=((*i).gmode GT 0 ? xrange : undef), $
        YLOG=((*i).gmode GT 0)
      OPLOT, omega, fftdata, COLOR=2
      OPLOT, omega_win, fftdata_win[iky,*], COLOR=3
      OPLOT, [0,0], !Y.CRANGE, COLOR=1

      IF (*i).gmode NE 0 THEN $
        OPLOT, omega_axis, gaussfunc, COLOR=1, THICK=3.0 ELSE $
        dummy = TEMPORARY(undef)

      IF (with_splines) THEN BEGIN
         OPLOT, omega_spline, fftdata_spline, COLOR=4

         ; maxima
         FOR imax = 0, n_freqs - 1 DO OPLOT, $
            [omega_max[0,imax,iky],omega_max[0,imax,iky]], !Y.CRANGE, $
            COLOR=4, LINESTYLE=2
                                ; first moment
         OPLOT, [omega_moms[iky],omega_moms[iky]], !Y.CRANGE, COLOR=5
      ENDIF
      plot_legend, [2,3,4], $
        ['FFT','Windowed FFT','Windowed FFT (splines)'], per_line=2

      plot_info_str, diag
    ENDFOR ; --- ky loop

    IF nkyind GT 1 THEN BEGIN
      xtitle = '!6k!Dy!N ' + get_var_string(/rhostr)
      ytitle = '!7x!6 / (' + get_var_string(-1,/time,/ounit) + ')'

      omax = MAX(omega_max[0,*,*],MIN=omin)
      PLOT, (*series.ky)[(*i).kyind], omega_max[0,0,*], COLOR=1, $
        YTITLE=ytitle, XTITLE=xtitle, /XSTYLE, /NODATA, $
        YRANGE=[omin,omax], POSITION=[0.15,0.4,0.9,0.9]
      OPLOT, !X.CRANGE, [0,0], COLOR=1

      FOR imax = 0, n_freqs - 1 DO OPLOT, $
        (*series.ky)[(*i).kyind], omega_max[0,imax,*], COLOR=imax+2
      OPLOT, (*series.ky)[(*i).kyind], omega_moms, COLOR=n_freqs+2

      plot_info_str, diag

      plot_legend, INDGEN(n_freqs+1) + 2, $
        [rm0es(INDGEN(n_freqs)+1)+'. mode','first moment'], per_line=2

      set_output, diag, header=['ky',rm0es(INDGEN(n_freqs)+1)+$
        '. mode','first moment'], commentline='kxind = '+$
        rm0es((*i).kxind), dat=[[(*series.ky)[(*i).kyind]],$
        [TRANSPOSE(omega_max[0,*,*])],[omega_moms]]

      FOR iky = 0, nkyind - 1 DO BEGIN
        win_max = MAX(fftdata_win[iky,*])
        IF win_max GT 0 THEN fftdata_win[iky,*] = $
          TEMPORARY(fftdata_win[iky,*]) / win_max
      ENDFOR

      kyaxis = (*series.ky)[(*i).kyind]
      kymin = MIN(kyaxis[WHERE(kyaxis GT 0)],MAX=kymax)

      yrange = 3.0 * MAX(ABS(omega_max[0,*,*])) * [-1,1]

      LOADCT, 33, FILE='internal/colortable.tbl'
      lev_col = contour_levels([MIN(fftdata_win),MAX(fftdata_win)],20)
      lev_col = contour_levels([0,1],20)
      CONTOUR, fftdata_win, (*series.ky)[(*i).kyind], omega_win, $
        LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
        /XSTYLE, /YSTYLE, /XLOG, xrange=[kymin,kymax], XTITLE=xtitle, $
        YTITLE=ytitle, CHARSIZE=csize, CHARTHICK=cthick, YRANGE=yrange, $
        TITLE='!12<!9!!!6'+get_var_string(0,/fancy)+'!9!!!6!U2!N!12>!6'
      PLOTS, 10^!X.CRANGE, [0,0], COLOR=255
      LOADCT, 41, FILE='internal/colortable.tbl'

      plot_info_str, diag
    ENDIF
  ENDELSE ; --- FFT mode

  set_output, diag, /reset

END
