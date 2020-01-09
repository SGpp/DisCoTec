FUNCTION svd_time_info

  RETURN, {$
    type      : 'extended',$
    title     : 'SVD time traces',$
    help_text : ['Plots information from SVD time traces.'],$
    ext_vars  : [['mnums','0','number of SVD modes to plot; default: '+$
                '10'],$
                 ['xrange','0','time window to plot; default: all'],$
                 ['yrange','0','yrange of the plot; default: min,max'],$
                 ['num_time','0','number of time steps and number of modes: must enter value'],$
                 ['fft_on','1','plot frequency spectrum; default: off']$
                 ]}

END

;######################################################################

PRO svd_time_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.mnums) LT 1 THEN (*e.mnums) = 10
  IF N_ELEMENTS(*e.xrange) LT 1 THEN (*e.xrange)=[-1,-1]
  IF N_ELEMENTS(*e.yrange) LT 1 THEN (*e.yrange)=[-1,-1]
  IF N_ELEMENTS(*e.num_time) LT 1 THEN print, "Error! Must enter num_time" 
  IF N_ELEMENTS(*e.fft_on) LT 1 THEN *e.fft_on = 0

  i = set_internal_vars(diag,{$
      mnums   : *e.mnums,$
      fft_on  : *e.fft_on,$
      xrange : *e.xrange,$
      num_time : *e.num_time,$
      yrange : *e.yrange})

END

;######################################################################

PRO svd_time_loop, diag

  COMMON global_vars

    i = (*diag).internal_vars

    svd_file_string='c_tot_svd'

    file_svd = get_file_path(svd_file_string)

    ;print,"file_svd",file_svd

    if file_svd.exist eq 0 then begin
      printerror, 'Error: c_tot file(s) not found'
      return
    endif

    last_time = -1.0d * 1e20

    svd_lun = 0l

    ;print,"Here we are in svd_time_loop!"
    ;print,"Reading c_tot_svd file."
    ;print,"Number of pod modes and number of time steps:", (*i).num_time

    if series.n_runs gt 1 then begin
      print,"n_runs must be 1."
      stop 
    endif

    ;time=FLTARR((*i).num_time)
    ;c_tot_in=COMPLEXARR((*i).num_time)
    ;c_modes=COMPLEXARR((*i).mnums,(*i).num_time)

    time_in = par.prec_single ? 0.0 : 0.0D
    time = par.prec_single ? $
      FLTARR((*i).num_time,/NOZERO) : $
      DBLARR((*i).num_time,/NOZERO)
    c_tot_in = par.prec_single ? $
      COMPLEXARR((*i).num_time,/NOZERO) : $
      DCOMPLEXARR((*i).num_time,/NOZERO)
    c_modes = par.prec_single ? $
      COMPLEXARR((*i).num_time,(*i).mnums,/NOZERO) : $
      DCOMPLEXARR((*i).num_time,(*i).mnums,/NOZERO)

    for run = 0, series.n_runs - 1 do begin
      run_label = (strsplit(series.run_labels,',',/extract))[run]
      ;print,"run_label",run_label
      ;print,"file_svd.path[run]",file_svd.path[run]

      ;print,"swap_endian",series.swap_endian
      openr, svd_lun, file_svd.path[run], /get_lun, error=err, $
        /f77_unformatted, swap_endian=series.swap_endian
      if err ne 0 then begin
        printerror, file_svd.path[run] + ' does not exist', /popup
	return
      endif else if !quiet ne 1 then print, 'reading ', $
        file_svd.path[run]

      print,'info time',size(time)
      for itime=0,(*i).num_time-1 do begin
         READU, svd_lun, time_in
         time[itime]=time_in
         print,itime,time[itime]
         READU, svd_lun, c_tot_in
         c_modes[itime,0:(*i).mnums-1]=c_tot_in[0:(*i).mnums-1]
      endfor

      ;plot, REAL_PART(c_modes[*,(*i).mnums-1])
      ;oplot, IMAGINARY(c_modes[*,(*i).mnums-1])

    endfor

  IF (*i).fft_on EQ 0 THEN BEGIN

    set_output, diag, /ps, suffix='_trace'

    FOR f = 0, (*i).mnums-1 DO BEGIN

      IF ((*i).yrange)[0] EQ -1 THEN BEGIN
        r_ct_min=MIN(REAL_PART(c_modes[*,f]))
        i_ct_min=MIN(IMAGINARY(c_modes[*,f]))
        r_ct_max=MAX(REAL_PART(c_modes[*,f]))
        i_ct_max=MAX(IMAGINARY(c_modes[*,f]))
        yrange = [MIN([r_ct_min,i_ct_min]),MAX([r_ct_max,i_ct_max])]
      ENDIF
      IF ((*i).xrange)[0] EQ -1 THEN BEGIN
        xrange=[MIN(time),MAX(time)]
      ENDIF ELSE BEGIN
        xrange=(*i).xrange
      ENDELSE

      ptitle="Time Amplitude, POD n="+STRTRIM(f,2)
      PLOT, time,REAL_PART(c_modes[*,f]), COLOR=1, $
           TITLE=ptitle,$
          /XSTYLE, XTITLE=get_var_string(1,/time,/fancy,/unit), $
          XRANGE=xrange,YRANGE=yrange, POSITION=[0.15,0.25,0.95,0.95]
      OPLOT, time,IMAGINARY(c_modes[*,f]), COLOR=2
      plot_legend, [1,2],['Re[h(t)]','Im[h(t)]']  

    ENDFOR

    set_output, diag, /reset, suffix='_trace'
  ENDIF ELSE BEGIN ;fft_on
    ;print,"Now do fft."

  set_output, diag, /ps, multi=[0,1,1], charsize=1.4, ysize=10,$
        suffix='_freq'
  FOR f=0, (*i).mnums-1 DO BEGIN
;;;;;;;;;;;;;Directly from frequency.pro
;;;;;;;;;;;;;Directly from frequency.pro
;;;;;;;;;;;;;Directly from frequency.pro
;;;;;;;;;Local Modifications
;;;;;;;;;Local Modifications
;series.step_count ===> (*i).num_time
;par.nkx0 ==> loc_nkx0
    loc_nkx0=1
    nkyind=1 ;nkyind ==> 1
    nkxind=1 ;nkxind ==> 1
; (*i).kyind[0] ==> 1 
    data = par.prec_single ? $
      COMPLEXARR(1,1,1,(*i).num_time,/NOZERO) : $
      DCOMPLEXARR(1,1,1,(*i).num_time,/NOZERO)
    data[0,0,0,*]=c_modes[*,f]
    with_splines=0
;;;;;;;;;Local Modifications Done
;;;;;;;;;Local Modifications Done

    n_freqs = 2 ; maximum number of freqs. (dominant+subdominant)
    twin = [0.25,0.5] ; windowed FFT parameters

    ytitle = 'h(t)' 

    ; prepare data interpolation for equidistant grid
    time_range = time[(*i).num_time-1] - time[0]
    ; choose an interpolation factor to achieve high FFT speeds
    n_steps_new = ROUND(1.5*(*i).num_time)
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


    omega_max = FLTARR(2,n_freqs,nkyind)
    omega_moms = FLTARR(nkyind)

    fftdata_win = FLTARR(nkyind,n_steps_win)

    FOR iky = 0, nkyind - 1 DO BEGIN
      data_ip = INTERPOLATE(REFORM(data[*,iky,*,*],$
          [1,1,(*i).num_time]),$
          INDGEN(1),INDGEN(1),equid_time,/GRID)

      fftdata_temp = TOTAL(ABS(FFT(data_ip,DIMENSION=3))^2,1) * $
        REBIN([1],[1,n_steps_new]) 
      fftdata = SHIFT(TOTAL(TEMPORARY(fftdata_temp),1),n_steps_new/2-1)
      fftdata /= nkxind * 1

      ; loop over twindows
      step_start_win = 0
      step_end_win = FIX(n_steps_new*twin[0]) - 1
      n_windows = 0
      IF step_end_win LE n_steps_new - 1 THEN win_hanning = $
        REBIN(REFORM(HANNING(n_steps_win),$
        [1,1,n_steps_win]),[nkxind,1,n_steps_win])
      WHILE step_end_win LE n_steps_new - 1 DO BEGIN
        data_filtered = REFORM(win_hanning*data_ip[*,*,$
          step_start_win:step_end_win],[nkxind,1,n_steps_win])
        fftdata_win_temp = TOTAL(ABS(FFT(data_filtered,DIMENSION=3))^2,1) * $
          REBIN([1],[1,n_steps_win]) 
        fftdata_win[iky,*] += $
          SHIFT(TOTAL(TEMPORARY(fftdata_win_temp),1),n_steps_win/2-1)

        step_start_win += FIX(n_steps_win*(1.0-twin[1]))
        step_end_win += FIX(n_steps_win*(1.0-twin[1]))
        n_windows += 1
      ENDWHILE

      fftdata_win[iky,*] /= nkxind * 1 * n_windows

      title = 'Frequency Spectrum (POD n='+STRTRIM(f+1,2) +')'
      ;title = 'k!Dy!N ' + get_var_string(/rhostr) + ' = ' + $
      ;  rm0es((*series.ky)[(*i).kyind[iky]],prec=3)

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

         ;IF !QUIET NE 1 THEN FOR j = 0, n_freqs - 1 DO PRINT, $
         ;   '   ' + rm0es(j+1) + '. frequency (windowed FFT, ky = ' + $
         ;   rm0es((*series.ky)[(*i).kyind[iky]],prec=3) + '): ' + $
         ;   rm0es(omega_max[0,j,iky])
         ; calculate first moment
         omega_moms[iky] = INT_TABULATED(omega_spline,omega_spline*$
            fftdata_spline) / INT_TABULATED(omega_spline,fftdata_spline)
         IF ((*i).xrange)[0] EQ -1 THEN BEGIN
          xrange = 3.0 * MAX(ABS(omega_max[0,*,iky])) * [-1,1]
         ENDIF ELSE BEGIN
          xrange=(*i).xrange
         ENDELSE
      ENDIF ELSE BEGIN
         xr_max = 3
         IF ((*i).xrange)[0] EQ -1 THEN BEGIN
          xrange = xr_max*[-1,1]
         ENDIF ELSE BEGIN
          xrange=(*i).xrange
         ENDELSE
      ENDELSE

      ;IF (*i).gmode EQ 1 THEN BEGIN
      ;  omega_axis = (FINDGEN(200) / 100.0 - 1.0) * 10.0
      ;  fftmax = MAX(fftdata,maxind)
      ;  omega_fftmax = omega[maxind]
      ;  gaussfunc = EXP(-((omega_axis-omega_fftmax)/(*i).gwidth)^2) * fftmax

;        xrange = [-10.0,10.0]
;      ENDIF
      ;IF (*i).gmode EQ 2 THEN BEGIN
      ;  ; reduce x axis to user specified range
      ;  min_ind = (WHERE(omega GE (*i).grange[0]))[0]
      ;  IF min_ind EQ -1 THEN min_ind = 0
      ;  max_ind = WHERE(omega LE (*i).grange[1],maxcount)
      ;  max_ind = (maxcount EQ -1) OR (maxcount EQ 0) ? 0 : max_ind[maxcount-1]
      ;  omega_axis = omega[min_ind:max_ind]
;
;        gaussfunc = GAUSSFIT(omega_axis,fftdata[min_ind:max_ind],$
;          coeffs,NTERMS=3)
        ;PRINT, 'width of Gaussian (ky=' + $
        ;  rm0es((*series.ky)[(*i).kyind[iky]],prec=3) + $
        ;  '): ' + rm0es(coeffs[2])
        ;PRINT, 'full width at half maximum (ky=' + $
        ;  rm0es((*series.ky)[(*i).kyind[iky]],prec=3) + $
        ;  '): ' + rm0es(2*SQRT(2*ALOG(2))*coeffs[2])

;        xrange = (*i).grange
;      ENDIF

      PLOT, omega, fftdata, COLOR=1, XTITLE='!7x!6 / ('+$
        get_var_string(-1,/time,/ounit)+')', YTITLE='!9!!!6'+ytitle+$
        '!9!!!6!U2!N', /XSTYLE, /NODATA, POSITION=[0.15,0.4,0.9,0.9], $
        TITLE=title, XRANGE=xrange, $
        YLOG=0
      OPLOT, omega, fftdata, COLOR=2
      OPLOT, omega_win, fftdata_win[iky,*], COLOR=3
      OPLOT, [0,0], !Y.CRANGE, COLOR=1

      ;PLOT, time,REAL_PART(c_modes[*,f]), COLOR=1, $
      ;     TITLE=ptitle,$
      ;    /XSTYLE, XTITLE=get_var_string(1,/time,/fancy,/unit), $
      ;    XRANGE=xrange,YRANGE=yrange, POSITION=[0.15,0.25,0.95,0.95]

      ;IF (*i).gmode NE 0 THEN $
      ;  OPLOT, omega_axis, gaussfunc, COLOR=1, THICK=3.0 ELSE $
      ;  dummy = TEMPORARY(undef)

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
  
   ENDFOR ;pod mode loop

   set_output, diag, /reset, suffix='_freq'
  ENDELSE ; --- FFT mode


   

END

;######################################################################

PRO svd_time_output, diag

  COMMON global_vars

;  e = (*diag).external_vars
;
;  set_output, diag, /ps
;
;  FOR f = 0, nlt_loop_max DO BEGIN
;
;    IF (*e.yrange)[0] EQ -1 THEN yrange = [MIN((*nltterms)[*,f*4:f*4+2]),MAX((*nltterms)[*,f*4:f*4+2])]
;    IF par.nlp_gdt THEN BEGIN
;        IF f LT par.num_nlt_pod_modes THEN BEGIN
;          ptitle='k!Dx!N='+STRING(kx_ind*kxmin,FORMAT='(F5.2)')+'k!Dy!N='$
;             +STRING(ky_ind*par.kymin,FORMAT='(F5.2)')+'(POD n='+STRTRIM(f+1,2)+')'
;        ENDIF ELSE BEGIN
;          ptitle='k!Dx!N='+STRING(kx_ind*kxmin,FORMAT='(F5.2)')+'k!Dy!N='$
;             +STRING(ky_ind*par.kymin,FORMAT='(F5.2)')+'(POD n=res.)'
;        ENDELSE
;    ENDIF ELSE BEGIN
;      ptitle='k!Dx!N='+STRING(kx_ind[f]*kxmin,FORMAT='(F5.2)')+'k!Dy!N='$
;           +STRING(ky_ind[f]*par.kymin,FORMAT='(F5.2)')
;    ENDELSE
;    PLOT, *nltterms_time, (*nltterms)[*,f*4+2], COLOR=1, $
;         TITLE=ptitle,$
;        /XSTYLE, XTITLE=get_var_string(1,/time,/fancy,/unit), $
;        XRANGE=xrange,YRANGE=yrange, POSITION=[0.15,0.25,0.95,0.95]
;    OPLOT, *nltterms_time, (*nltterms)[*,f*4], COLOR=4
;    OPLOT, *nltterms_time, (*nltterms)[*,f*4+1], COLOR=2
;    plot_legend, [1,4,2],['E!Dk!N','dE/dt!D NL !N','dE/dt!D lin !N']  
;
;  ENDFOR
;
;  set_output, diag, /reset

END
