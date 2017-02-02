FUNCTION zonal_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'Zonal flow diagnostic',$
    help_text : ['Plots time traces of shearing rate and, if beta>0, of '+$
        'magn. shear fluctuations (the averages shown in the plots are '+$
        'usually to be compared to the maximum linear growth rate). '+$
        'Also, kx spectra of zonal potential/field and of shearing/magn. '+$
        'shear fluct. rates are shown. If an fft_range is provided, '+$
        'a contour of FFT(phi_zonal)(kx,omega) is presented (e.g. to '+$
        'identify GAM freqs.). By inserting a decorrelation time, '+$
        'a correction to the shearing rate '+$
        '(cmp. Phys. Plasmas, Vol. 6, p.923) will be calculated'],$
    ext_vars  : [$
        ['var','0','variable to be analyzed; default: -1 (phi and '+$
         'additionally shearing rate if beta is finite'],$
        ['kxind','0','array of phi_zonal kx-modes to be plotted; '+$
         'default: [1,..,nx0/2-1]'],$
        ['zind','0','index of specific z position; '+$
         'default: -1 (averaging over z for zonal mode)'],$
        ['fft_range','0','2d array containing neg. and pos. freq. limit for '+$
         '(windowed) FFT data; default: 0 (no FFT will be shown); -1 for automatic'],$
        ['corrtime','0','insert a decorrelation time of the underlying '+$
         'turbulence to calculate a correction; default: -1 (off)']]}

END

;######################################################################

PRO zonal_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  apar = 0
  IF N_ELEMENTS(*e.var) NE 1 THEN *e.var = -1
  IF *e.var LT 1 THEN BEGIN
    IF par.n_fields GE 2 THEN apar = 1
    *e.var = 0
  ENDIF

  IF par.lx EQ 0 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + ': ky dependent lx not allowed'
    (*diag).selected = 0
    RETURN
  ENDIF

  valid = 1
  FOR ikx = 0, N_ELEMENTS(*e.kxind)-1 DO $
    valid = valid AND ((*e.kxind)[ikx] LT par.nkx0)

  IF (N_ELEMENTS(*e.kxind) LT 1) OR NOT valid THEN BEGIN
    IF par.nkx0 GT 4 THEN *e.kxind = indgen(par.nkx0/2-1)+1 $
      ELSE *e.kxind = 1
  ENDIF
  IF (N_ELEMENTS(*e.zind) LT 1) THEN *e.zind = -1

  IF N_ELEMENTS(*e.fft_range) LT 1 THEN fft_range = [0,0] $
  ELSE fft_range=*e.fft_range
  IF N_ELEMENTS(fft_range) EQ 1 THEN $
    IF fft_range EQ -1 THEN fft_range = [-1,-1] ELSE fft_range = [0,0]

  IF N_ELEMENTS(*e.corrtime) NE 1 THEN *e.corrtime = -1

  IF *e.corrtime GE 0 THEN BEGIN
    printerror, 'Need to use all kx modes for shearing rate correction' + $
     ' (will be done automatically)'
    *e.kxind = [par.nkx0/2+1+INDGEN(par.nkx0/2-1),INDGEN(par.nkx0/2)]
;    *e.kxind = SHIFT(INDGEN(par.nkx0),par.nkx0/2-1)
  ENDIF

  i = set_internal_vars(diag,{$
    var           : *e.var,$
    kxind         : *e.kxind,$
    zind          : *e.zind,$
    fft_range	  : fft_range,$
    corrtime      : *e.corrtime,$
    shear_id      : PTR_NEW(),$
    shear_avg_id  : PTR_NEW(),$
    sh_fl_id      : apar ? PTR_NEW() : 0,$
    sh_fl_avg_id  : apar ? PTR_NEW() : 0,$
    time_id       : PTR_NEW(),$
    phi_zfl_id    : PTR_NEW(),$
    apar_zfl_id   : PTR_NEW()})

  vars = apar ? [*e.var,1] : *e.var
  fft_format, kxky=vars

END

;######################################################################

PRO zonal_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkxind = N_ELEMENTS((*i).kxind)
  phi_zfl_temp = DCOMPLEXARR(par.nkx0,gui.out.n_spec_sel) ; phi_zonal
  phi_zfl = DCOMPLEXARR(nkxind,gui.out.n_spec_sel)
  IF gui.out.n_spec_sel LT 2 THEN phi_zfl = $
    REFORM(phi_zfl,[nkxind,gui.out.n_spec_sel],/OVERWRITE)
  shear = DBLARR(gui.out.n_spec_sel) ; phi_zonal

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    ; Calculate flux surface averages with Jacobi det
    ; phi_zonal(kx) = <phi(kx,0,z)>_z
    IF ((*i).zind[0] GT -1) THEN phi_zfl_temp[0,isp] = $
       REFORM((*mom[isp,(*i).var].kxky)[*,0,(*i).zind[0]]) $
    ELSE phi_zfl_temp[0,isp] = TOTAL(REFORM((*mom[isp,(*i).var].kxky)[*,0,*],$
      [par.nkx0,par.nz0])*REBIN(REFORM((*series.geom).jac_norm/par.nz0,$
      [1,par.nz0]),[par.nkx0,par.nz0]),2)

    Omega = - (*series.kx)[*,0]^2/(*series.geom).C_xy * phi_zfl_temp[*,isp]
    phi_zfl[0,isp] = phi_zfl_temp[(*i).kxind,isp]
    shear[isp] = SQRT(TOTAL(ABS(Omega)^2,/DOUBLE))
  ENDFOR
  (*i).phi_zfl_id = store_step((*i).phi_zfl_id,phi_zfl)

  ; Calculate ExB shearing rate using all kx modes
  ; Omega = <d/dx (i kx <phi(kx,0,z)>_z
  (*i).shear_id = store_step((*i).shear_id,shear)
  (*i).shear_avg_id = time_avg((*i).shear_avg_id,shear,mom_time)

  ; Calculate magnetic shear fluctuations
  ; shfl = q R d/dx B_y (see e.g. PoP 9, 4103)
  IF (par.n_fields GE 2) AND ((*i).var EQ 0) THEN BEGIN
    IF ((*i).zind[0] GT -1) THEN apar_zfl = $
       REFORM((*mom[0,1].kxky)[*,0,(*i).zind[0]]) $
    ELSE apar_zfl = TOTAL(REFORM((*mom[0,1].kxky)[*,0,*],[par.nkx0,par.nz0])*$
       REBIN(REFORM((*series.geom).jac_norm/par.nz0,[1,par.nz0]),$
             [par.nkx0,par.nz0]),2)

    (*i).apar_zfl_id = store_step((*i).apar_zfl_id,apar_zfl[(*i).kxind])
    shfl = (*series.kx)[*,0]^2 /(*series.geom).C_xy * apar_zfl
    shear_fluc = par.q0 * SQRT(TOTAL(ABS(shfl)^2,/DOUBLE))
    (*i).sh_fl_id = store_step((*i).sh_fl_id,shear_fluc)
    (*i).sh_fl_avg_id = time_avg((*i).sh_fl_avg_id,shear_fluc,mom_time)
  ENDIF

  (*i).time_id = store_step((*i).time_id,mom_time)

END

;######################################################################

PRO zonal_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkxind = N_ELEMENTS((*i).kxind)
  IF (par.n_fields GE 2) AND ((*i).var EQ 0) THEN $
    apar = 1 ELSE apar = 0

  time = store_step((*i).time_id,/get,/reg_array)
  T_ = time[series.step_count-1] - time[0]

  phi_zfl = store_step((*i).phi_zfl_id,/get,/reg_array)
  shear_rate = store_step((*i).shear_id,/get,/reg_array)
  shear_avg = time_avg((*i).shear_avg_id,/avg)
  IF apar THEN BEGIN
    apar_zfl = store_step((*i).apar_zfl_id,/get,/reg_array)
    shear_fluc = store_step((*i).sh_fl_id,/get,/reg_array)
    shear_fluc_avg = time_avg((*i).sh_fl_avg_id,/avg)
  ENDIF

  rho_str = ' ' + get_var_string(/rhostr)
  ddx_str='!6!S!E !N!Ud!R!S-!R!Ddx!X!N'
  shear_str = '!7x!6!DE!N'
  om_kx_av_str = '!12<!9!!' + shear_str + $
    '!9!!!6!U2!N!12>!6!S!Dx!N!R!U1/2!N'
  shfl_str = '!6!Ss!R!E-!N!6'
  sh_fl_av_str = '!12<!9!!' + shfl_str + $
    '!9!!!6!U2!N!12>!6!S!Dx!N!R!U1/2!N'
  zstr = ((*i).zind EQ -1)?'k!D!9#!6!N=0':$
          'z='+rm0es((*par.z)[(*i).zind])

  shcol = 4

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, multi=[0,1,2]
    ; --- plot time traces ---
    PLOT, time, shear_rate[isp,*], COLOR=1, /XSTYLE, $
      TITLE=get_var_string((*i).var,/fancy)+' shearing rate '+shear_str+'='+$
      ddx_str+'v!DE,y!N(k!Dy!N=0,'+zstr+')', $
      YTITLE=om_kx_av_str, YMARGIN=[3,4],$
      XTITLE=get_var_string(0,/time,/units), charsize=1.4
    OPLOT, !X.CRANGE,[shear_avg[isp],shear_avg[isp]], COLOR=shcol,$
      LINESTYLE=1
    deltax = !X.CRANGE[1] - !X.CRANGE[0]
    deltay = !Y.CRANGE[1] - !Y.CRANGE[0]

    XYOUTS, 0.1*deltax+!X.CRANGE[0],0.1*deltay+!Y.CRANGE[0],color=shcol,charsize=1.5,$
      '!12<!6'+om_kx_av_str+'!12>!6!Dt!N = '+rm0es(shear_avg[isp],prec=3)+' '+$
      get_var_string(-1,/time,/ounit)
      
    IF apar THEN BEGIN
      PLOT, time, shear_fluc, COLOR=1, /XSTYLE, $
        TITLE='magn. shear fluct. '+shfl_str+'=qR'+ddx_str+$
        'B!Dy!N(k!Dy!N=0,'+zstr+')', $
        YTITLE=sh_fl_av_str, YMARGIN=[4,3], $
        XTITLE=get_var_string(0,/time,/units), CHARSIZE=1.4
      OPLOT, !X.CRANGE,[shear_fluc_avg,shear_fluc_avg], COLOR=shcol,$
        LINESTYLE=1
      deltax = !X.CRANGE[1] - !X.CRANGE[0]
      deltay = !Y.CRANGE[1] - !Y.CRANGE[0]

      XYOUTS, 0.1 * deltax + !X.CRANGE[0], 0.1 * deltay + !Y.CRANGE[0], $
        COLOR=shcol, CHARSIZE=1.5, '!12<!6'+sh_fl_av_str+'!12>!6!Dt!N = '+$
        rm0es(shear_fluc_avg,prec=3)+' '+get_var_string(-1,/time,/ounit)

      set_output, diag, sp, header = ['t','shear. rate','shear fluct.'],$
        commentlines=['avg. values:',get_var_string((*i).var)+$
        ' shearing rate: '+rm0es(shear_avg[isp]),$
        'magn. shear fluct: '+rm0es(shear_fluc_avg)],dat=[[time],[REFORM(shear_rate[isp,*])],[shear_fluc]]
    ENDIF ELSE set_output, diag, sp, header=['shear. rate'], $
      commentlines=['avg. values:','shearing rate: '+rm0es(shear_avg[isp])], $
      dat=[shear_rate[isp,*]]

    XYOUTS, 0.01,0.975, 'All k!Dx!N modes included:',/NORMAL,color=1,charsize=1.5

    plot_info_str, diag

    PRINT, 'Zonal: shearing rate: ', rm0es(shear_avg[isp]), $
      (apar ? ', magn. shearing fluctuation: ' + rm0es(shear_fluc_avg) : '')

    ; --- plot spectra ---
    kxaxis = (*series.kx)[(*i).kxind,0]

    xlog = ((WHERE(kxaxis LE 0))[0] EQ -1) 
  
    phi_zfl_kx = DBLARR(nkxind,/NOZERO)

    FOR ikx = 0, nkxind - 1 DO $
      phi_zfl_kx[ikx] = INT_TABULATED(time[*],ABS(phi_zfl[ikx,isp,*])) / T_
    ymax = MAX(phi_zfl_kx[WHERE(kxaxis NE 0.0)],min=ymin)
    legend_str = get_var_string((*i).var,/fancy) + '!6!Dzonal!N'
    ytitle='!6zonal potential'
    IF apar THEN BEGIN
      apar_zfl_kx = DBLARR(nkxind,/NOZERO)  
      FOR ikx=0, nkxind-1 DO $
        apar_zfl_kx[ikx] = INT_TABULATED(time,ABS(apar_zfl[ikx,*]))/T_
      ymax = MAX([ymax,apar_zfl_kx[WHERE(kxaxis NE 0.0)]],min=mymin)
      ymin=MIN([ymin,mymin])    
      legend_str = [legend_str,'!6A!D!9#!6,zonal!N']
      ytitle += '/field'
    ENDIF    
    PLOT, kxaxis, phi_zfl_kx, YMARGIN=[6,2], charsize=1.4,$
      COLOR=1, /YLOG, /XSTYLE, XTITLE='!6k!Dx!N'+rho_str, $
      YTITLE=ytitle, XLOG=xlog, YRANGE=[ymin,ymax]
  ;    title='Considering only user specified k!Dx!N modes'
    IF apar THEN OPLOT, kxaxis, apar_zfl_kx, COLOR=shcol
    legend_str = '!12<!9!!!6'+legend_str+'!9!!!12>!6!Dt!N'  
    plot_legend,INDGEN(apar+1)*(shcol-1)+1,legend_str,per_line = 2, $
      x_offset=[0.25,apar*0.5,0.75,apar*0.5+0.05]  
    XYOUTS, 0.01,apar*0.5+0.475, 'Only user specified k!Dx!N:',/NORMAL,color=1,charsize=1.5
  
    ;plot shearing rate
    shear_kx = DBLARR(nkxind,/NOZERO)  
    temp = DBLARR(series.step_count,/NOZERO)  
    FOR ikx=0, nkxind -1 DO BEGIN
      FOR tind=0, series.step_count-1 DO $
        temp[tind] = ABS((*series.kx)[(*i).kxind[ikx],0]^2/(*series.geom).C_xy*$
                         phi_zfl[ikx,isp,tind])
      shear_kx[ikx] = INT_TABULATED(time,temp)/T_
    ENDFOR
    ymax = MAX(shear_kx[WHERE(kxaxis NE 0.0)],min=ymin)
    legend_str = shear_str
    IF apar THEN BEGIN
      shfl_kx = DBLARR(nkxind,/NOZERO)  
      temp = DBLARR(series.step_count,/NOZERO)  
      FOR ikx=0, nkxind -1 DO BEGIN
        FOR tind=0, series.step_count-1 DO $
          temp[tind] = ABS((*series.kx)[(*i).kxind[ikx],0]^2*par.q0*apar_zfl[ikx,tind])
        shfl_kx[ikx] = INT_TABULATED(time,temp)/T_
      ENDFOR
      ymax=MAX([ymax,shfl_kx[WHERE(kxaxis NE 0.0)]],min=mymin)
      ymin=MIN([ymin,mymin])
      legend_str = [legend_str,shfl_str]
      set_output, diag, sp, header = ['kx','shear. rate','shear fluct.'],$
        commentlines=['averaged over time'],dat=[[kxaxis],[shear_kx],[shfl_kx]],/append    
    ENDIF ELSE set_output, diag, sp, header = ['kx','shear. rate'],$
      commentlines=['averaged over time'],dat=[[kxaxis],[shear_kx]],/append

    PLOT, kxaxis, shear_kx, charsize=1.4,$
      COLOR=1, XLOG=xlog, /YLOG, /XSTYLE, XTITLE='!6k!Dx!N'+rho_str, $
      YTITLE='rates / ('+get_var_string(-1,/time,/ounit)+')', $
      YRANGE=[ymin,ymax],YMARGIN=[6,2]
    IF apar THEN OPLOT, kxaxis, shfl_kx, COLOR=shcol
    legend_str = '!12<!9!!'+legend_str+'!9!!!12>!6!Dt!N'  
    plot_legend,INDGEN(apar+1)*(shcol-1)+1,legend_str,per_line = 2, $
      x_offset=[0.25,(1-apar)*0.5,0.75,(1-apar)*0.5+0.05]

    print, 'rms(shear_tavg) = ',SQRT(TOTAL(ABS(shear_kx)^2))

    ; --- plot time FFT ---
    IF (TOTAL(ABS((*i).fft_range)) NE 0) THEN BEGIN ; frequency spectrum of zonal potential
      temp = time_fft(time,phi_zfl[0,isp,*],/SHIFT,omega=my_om_arr)    
      phi_zfl_fft = DBLARR(nkxind,N_ELEMENTS(temp))  
      om4corr = DBLARR(nkxind)
      phi_zfl_fft[0,*] = temp
      FOR ikx=1, nkxind-1 DO BEGIN  
        phi_zfl_fft[ikx,*] = time_fft(time,phi_zfl[ikx,isp,*],/SHIFT)
        om4corr[ikx] = INT_TABULATED(my_om_arr,phi_zfl_fft[ikx,*]*my_om_arr)/$
          INT_TABULATED(my_om_arr,(phi_zfl_fft[ikx,*]))	
        phi_zfl_fft[ikx,*] *= 1.0/MAX(phi_zfl_fft[ikx,*])
      ENDFOR 

      LOADCT, 33      
      lev_col = contour_levels(phi_zfl_fft,40,/LOG)

      yrange = [MIN(my_om_arr),MAX(my_om_arr)]
      IF ((*i).fft_range[0] NE -1) THEN yrange[0] = (*i).fft_range[0]
      IF ((*i).fft_range[1] NE -1) THEN yrange[1] = (*i).fft_range[1]
    
      range = WHERE(kxaxis NE 0.0)
      IF range[0] EQ -1 THEN range=INDGEN(nkxind)
    
      CONTOUR, phi_zfl_fft[range,*], kxaxis[range], my_om_arr, LEVELS=lev_col[*,0], $
        C_COLORS=lev_col[*,1], /FILL, COLOR=1, /XSTYLE, /YSTYLE, $
        YRANGE=yrange,$
        TITLE='!6norm. frequency spectrum of ' + get_var_string((*i).var,/fancy) + '!6!Dzonal!N', $
        XTITLE='!6k!Dx!N'+rho_str,YTITLE='!7x!D!6re!N '+get_var_string(0,/time,/ounit)
      LOADCT, 41, FILE='internal/colortable.tbl'
    ENDIF

    ; --- compute corrections ---  
    IF (*i).corrtime GE 0 THEN BEGIN
      ;"measure" frequency of Re[phi_zonal] by searching for extrema    
      om_f_av = FLTARR(nkxind)
      FOR ikx=0, nkxind-1 DO BEGIN
        temp = FLOAT(phi_zfl[ikx,isp,*])
        datder=DERIV(time,temp)
        n_extr = 0
        t_av = 0
        FOR tind=1, series.step_count -3 DO BEGIN
          IF ((datder[tind] GE 0) AND (datder[tind+1] LT 0)) OR $
          ((datder[tind] LE 0) AND (datder[tind+1] GT 0)) THEN BEGIN
            n_extr = n_extr + 1
            IF n_extr GT 1 THEN t_av = t_av + time[tind] - t_last
            t_last=time[tind]    
          ENDIF
        ENDFOR
        IF n_extr GT 1 THEN om_f_av[ikx] = !PI/((t_av)/(n_extr-1)) $
        ELSE om_f_av[ikx]=0
      ENDFOR   
      
      ; plot of "measured" frequencies
      PLOT, kxaxis, om_f_av,YMARGIN=[6,2],$
        COLOR=1, PSYM=8, XTITLE='!6k!Dx!N'+rho_str, charsize=1.5,$
        YTITLE='!7x!6!Df!N / ('+get_var_string(-1,/time,/ounit)+')',$
        TITLE='"Measured" freq. for '+shear_str+' correction'
      PLOTS, !X.CRANGE,[0,0],color=1 ; draw base line    
      deltax = !X.CRANGE[1] - !X.CRANGE[0]
      deltay = !Y.CRANGE[1] - !Y.CRANGE[0]
      XYOUTS, 0.01*deltax+!X.CRANGE[0],0.025*deltay+!Y.CRANGE[0],color=1,charsize=1.5,$
      '!6t!Icorr!N = '+rm0es((*i).corrtime,prec=3)+' '+$
      get_var_string(0,/time,/ounit) 
  ;    PLOTS, !X.CRANGE,[om_re_min,om_re_min],color=1, LINESTYLE=2 ; draw min freq line   

     ; calculate shearing rate correction
     ; cmp. Phys. Plasmas, Vol. 6, p.923
      shear_corr_kx = DBLARR(nkxind,/NOZERO)
      temp = DBLARR(series.step_count,/NOZERO)
      FOR ikx = 0, nkxind - 1 DO BEGIN    
        F = (om_f_av[ikx] * (*i).corrtime)^2
        corr_fac = ((1 + 3 * F)^2 + 4 * F^3)^0.25 / ((1 + F) * (1 + 4 * F)^0.5)    
        FOR tind = 0, series.step_count - 1 DO $
          temp[tind] = ABS((*series.kx)[(*i).kxind[ikx],0]^2/$
                           (*series.geom).C_xy * phi_zfl[ikx,isp,tind] * $
                           corr_fac)
        shear_corr_kx[ikx] = INT_TABULATED(time,temp)/T_
      ENDFOR

      ymax = MAX(shear_kx[WHERE(kxaxis NE 0.0)],min=ymin)
      ymax2 = MAX(shear_corr_kx[WHERE(kxaxis NE 0.0)],min=ymin2)
      ymax = ymax > ymax2
      ymin = ymin < ymin2
      PLOT, kxaxis, shear_kx, $
        COLOR=1, XTITLE='!6k!Dx!N'+rho_str, /XSTYLE,$
        YTITLE='!9!!!7x!9!!!6 / ('+get_var_string(-1,/time,/ounit)+')',$
        POSITION=[0.15,0.15,0.95,0.5],charsize=1.5,YMARGIN=[6,2],$
        TITLE='corrected and uncorrected '+shear_str, /YLOG, YRANGE=[ymin,ymax]
      OPLOT, kxaxis, shear_corr_kx, COLOR=2
      PLOTS, !X.CRANGE,[0,0],color=1 ; draw base line

      shear_corr = DBLARR(series.step_count,/NOZERO)
      temp = DBLARR(nkxind,/NOZERO)
      FOR tind = 0, series.step_count - 1 DO BEGIN
        FOR ikx = 0, nkxind - 1 DO BEGIN    
          F = (om_f_av[ikx] * (*i).corrtime)^2
          corr_fac = ((1 + 3 * F)^2 + 4 * F^3)^0.25 / ((1 + F) * (1 + 4 * F)^0.5)
          temp[ikx] = ABS((*series.kx)[(*i).kxind[ikx],0]^2/$
              (*series.geom).C_xy * phi_zfl[ikx,isp,tind] * corr_fac)^2
        ENDFOR	
        shear_corr[tind] = SQRT(TOTAL(temp))
      ENDFOR
      shear_corr_avg = INT_TABULATED(time,shear_corr)/T_
    
      PRINT, 'shearing rate: ', rm0es(shear_avg[isp])
      PRINT, 'Corrected shearing rate: ', rm0es(shear_corr_avg) 

      deltax = !X.CRANGE[1] - !X.CRANGE[0]
      deltay = !Y.CRANGE[1] - !Y.CRANGE[0]
      XYOUTS, 0.01 * deltax + !X.CRANGE[0], 10^(0.025 * deltay + !Y.CRANGE[0]), $
        COLOR=2, CHARSIZE=1.5, '!12<!6'+om_kx_av_str+'!12>!6!Dt!Ucorr!N = '+$
        rm0es(shear_corr_avg,prec=3)+' '+get_var_string(-1,/time,/ounit)

      ; print legend
      plot_legend, [1,2], ['!9!!'+shear_str+'!9!!!6(k!Dx!N)',$
        '!9!!!7x!6!DE,corr!N!9!!!6(k!Dx!N)'],$
        x_offset=[0.15,0.0,0.95,0.075], per_line = 2
    
      set_output, diag, sp, header = ['kx','shear. corr.'],$
        commentlines=['Used corr. time of '+rm0es((*i).corrtime)+' for correction',$
        'New corrected avg. shearing rate:'+rm0es(shear_corr_avg)],$
        dat=[[kxaxis],[shear_corr_kx]],/append
    ENDIF ;corrtime

    plot_info_str, diag

    set_output, diag, sp, /reset
  ENDFOR

END
