FUNCTION timetrace_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'Time trace',$
    help_text : ['Draws a linear-linear and a log-linear plot of the time trace '+$'
    	'of the real, imaginary and absolute value of selected modes and variable. '+$
    	'The corresponding kx, ky and z values can be found in the legend. '+$
	'Data output will only be done for first choice of var, kxind, kyind and zind.'],$    
    ext_vars  : [$		 
    	['var','0','index of the fluctuating quantity to be plotted'+$
	 ' (see variable list)'],$
    	['kxind','0','index of kx modes to be plotted (integer '+$
	 'between 0 and nkx0 - 1; default: 0)'],$
    	['kyind','0','index of ky modes to be plotted (integer '+$
	 'between 0 and nky0 - 1; default: 1)'],$
    	['zind','0','index of z-value to be plotted (integer '+$
	 'between 0 and nz0 - 1; default: -1 (average))'],$
    	['wsmooth','0','insert smooth window width; default: off']$    
    	]}
 		 
END

;######################################################################

PRO timetrace_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var) LT 1 THEN *e.var = 0
  IF N_ELEMENTS(*e.kxind) LT 1 THEN *e.kxind = 0
  IF N_ELEMENTS(*e.kyind) LT 1 THEN *e.kyind = par.ky0_ind EQ 0
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = -1
  IF N_ELEMENTS(*e.wsmooth) LT 1 THEN *e.wsmooth = 0

  i = set_internal_vars(diag,{$
    vars    : (*e.var),$
    kxind   : *e.kxind,$
    kyind   : *e.kyind,$
    zind    : *e.zind,$
    wsmooth : *e.wsmooth,$
    time_id : PTR_NEW(),$
    data_id : PTR_NEW()})

  IF par.x_local THEN fft_format, kxky=(*e.var) $
    ELSE fft_format, sxky=(*e.var)

END

;######################################################################

PRO timetrace_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  n_kxind = N_ELEMENTS((*i).kxind)
  n_kyind = N_ELEMENTS((*i).kyind)
  n_zind = N_ELEMENTS((*i).zind)
  n_vars = N_ELEMENTS((*i).vars)

  tt_data = REFORM(DCOMPLEXARR($
    n_kxind,n_kyind,n_zind,n_vars,gui.out.n_spec_sel,/NOZERO),$
    [n_kxind,n_kyind,n_zind,n_vars,gui.out.n_spec_sel])
  
  IF par.x_local THEN BEGIN
   FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     FOR ivar = 0, n_vars - 1 DO BEGIN
       FOR ikx = 0, n_kxind - 1 DO BEGIN
         FOR iky = 0, n_kyind - 1 DO BEGIN
           FOR ins = 0, n_zind - 1 DO BEGIN
             IF (*i).zind[ins] EQ - 1 THEN BEGIN
               tt_data[ikx,iky,ins,ivar,isp] = COMPLEX(0.0,0.0)
               FOR nsi = 0, par.nz0 - 1 DO BEGIN
                 tt_data[ikx,iky,ins,ivar,isp] += $
                   (*mom[isp,(*i).vars[ivar]].kxky)[(*i).kxind[ikx],$
                   (*i).kyind[iky],nsi]*(*series.geom).jac_norm[nsi]
               ENDFOR
               tt_data[ikx,iky,ins,ivar,isp] /= par.nz0
             ENDIF ELSE tt_data[ikx,iky,ins,ivar,isp] = $
               (*mom[isp,(*i).vars[ivar]].kxky)[(*i).kxind[ikx],$
               (*i).kyind[iky],(*i).zind[ins]]
           ENDFOR
         ENDFOR
       ENDFOR
     ENDFOR
   ENDFOR
  ENDIF ELSE BEGIN
   FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     FOR ivar = 0, n_vars - 1 DO BEGIN
       FOR ikx = 0, n_kxind - 1 DO BEGIN
         FOR iky = 0, n_kyind - 1 DO BEGIN
           FOR ins = 0, n_zind - 1 DO BEGIN
             IF (*i).zind[ins] EQ - 1 THEN BEGIN
               tt_data[ikx,iky,ins,ivar,isp] = COMPLEX(0.0,0.0)
               FOR nsi = 0, par.nz0 - 1 DO BEGIN
                 tt_data[ikx,iky,ins,ivar,isp] += $
                   (*mom[isp,(*i).vars[ivar]].sxky)[(*i).kxind[ikx],$
                   (*i).kyind[iky],nsi]*(*series.geom).jac_norm[nsi]
               ENDFOR
               tt_data[ikx,iky,ins,ivar,isp] /= par.nz0
             ENDIF ELSE tt_data[ikx,iky,ins,ivar,isp] = $
               (*mom[isp,(*i).vars[ivar]].sxky)[(*i).kxind[ikx],$
               (*i).kyind[iky],(*i).zind[ins]]
          ENDFOR
         ENDFOR
       ENDFOR
     ENDFOR
   ENDFOR  
  ENDELSE

  (*i).time_id = store_step((*i).time_id,mom_time)
  (*i).data_id = store_step((*i).data_id,tt_data)

END

;######################################################################

PRO timetrace_plot_data, diag, time, tt_data, title, legend_str, text_tit

  COMMON global_vars
  i = (*diag).internal_vars

  n_kxind = N_ELEMENTS(tt_data[*,0,0,0,0])
  n_kyind = N_ELEMENTS(tt_data[0,*,0,0,0])
  n_zind  = N_ELEMENTS(tt_data[0,0,*,0,0])
  n_vars  = N_ELEMENTS(tt_data[0,0,0,*,0])
  
  maxval = MAX(tt_data, min=minval)  
    
; lin - lin plot
  PLOT, time, tt_data[0,0,0,0,*], COLOR = 1, $
    /XSTYLE, YRANGE=[minval,maxval], $
    TITLE=title, CHARSIZE=csize,$
    POSITION=[0.15,0.6,0.95,0.9615]

  count = 0

  FOR ivar=0, n_vars-1 DO BEGIN
   FOR ikx=0, n_kxind-1 DO BEGIN
    FOR iky=0, n_kyind-1 DO BEGIN
     FOR ins=0, n_zind-1 DO BEGIN   
      OPLOT, time, tt_data[ikx,iky,ins,ivar,*], COLOR = 1 + count      
; try to fit GAMS
IF 0 THEN BEGIN       
     Y=REFORM(tt_data[ikx,iky,ins,ivar,*]/MAX(tt_data[ikx,iky,ins,ivar,*]))
     weights = FLTARR(series.step_count)+1.0; /REFORM(tt_data[ikx,iky,ins,ivar,*])
     A = [0.001,0.02,2.*!PI/20.,0.1]
     gamfunct2, time, A, ygam
     OPLOT, time, ygam*MAX(tt_data[ikx,iky,ins,ivar,*]), COLOR=1 + count, $
        LINESTYLE=3
	
;     A = [0.0001,0.01,-0.05]          
     yfit = CURVEFIT(time, Y, weights, A, SIGMA, FUNCTION_NAME='gamfunct2')     
     PRINT, 'A= '
     PRINT, A
     OPLOT, time, yfit*MAX(tt_data[ikx,iky,ins,ivar,*]), COLOR=1 + count, $
        LINESTYLE=2     
ENDIF	
     count = count + 1 	
    ENDFOR
   ENDFOR   
  ENDFOR
 ENDFOR  
  ; add base line  

  IF !Y.CRANGE[0] LE 0 AND !Y.CRANGE[1] GE 0 THEN $
    oplot, !X.CRANGE, [0,0], color=1

; lin -log plot  
  IF minval LT 1e-3 THEN minval=1e-3
  PLOT, time, tt_data[0,0,0,0,*], COLOR = 1, $
    /XSTYLE, YRANGE=[minval,maxval], /YLOG,$
    XTITLE=get_var_string(0,/time,/units),$
    CHARSIZE=csize,POSITION=[0.15,0.2,0.95,0.55]

  count = 0
  PRINT,'timetrace: '+text_tit
  FOR ivar=0, n_vars-1 DO BEGIN
   FOR ikx=0, n_kxind-1 DO BEGIN
    FOR iky=0, n_kyind-1 DO BEGIN
     FOR ins=0, n_zind-1 DO BEGIN
      IF par.x_local THEN line = 'kx='+rm0es((*series.kx)[(*i).kxind[ikx],(*i).kyind[iky]],prec=3)+$
       ', ky= '+rm0es((*series.ky)[(*i).kyind[iky]],prec=3) $
      ELSE line = 'x='+rm0es((*series.lxarr)[(*i).kyind[iky]]*(-0.5+(*i).kxind[ikx]/par.nx0),prec=3)+$
       ', ky= '+rm0es((*series.ky)[(*i).kyind[iky]],prec=3)
      IF (*i).zind[ins] EQ -1 THEN line = line + ', z avg.' $
      ELSE line = line + ',z= '+rm0es((*par.z)[(*i).zind[ins]])
     
      OPLOT, time, tt_data[ikx,iky,ins,ivar,*], COLOR = 1 + count

      res = LINFIT(time,ALOG(REFORM(tt_data[ikx,iky,ins,ivar,*])))
      PRINT, line+': fit f(x) = exp('+rm0es(res[1])+'x + ' +rm0es(res[0])+')'
      PRINT, line+': median = '+rm0es(MEDIAN(tt_data[ikx,iky,ins,ivar,*]),prec=3)+' ('$
       +rm0es(MEDIAN(tt_data[ikx,iky,ins,ivar,*])/$
              MAX(tt_data[ikx,iky,ins,ivar,*]),prec=3)+' norm. to max.)'
      PRINT, line + ': avg = ' + rm0es(TOTAL(tt_data[ikx,iky,ins,ivar,*]/$
        N_ELEMENTS(tt_data[ikx,iky,ins,ivar,*])),prec=3)
      count = count + 1
     ENDFOR
    ENDFOR
   ENDFOR
  ENDFOR

  ; print legend
  plot_legend, INDGEN(n_kxind*n_kyind*n_zind*n_vars) + 1, $
      legend_str, x_offset=[0.15,0.0,0.95,0.12]

END ; plot time_trace

;######################################################################

PRO timetrace_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
  n_kxind = N_ELEMENTS((*i).kxind)
  n_kyind = N_ELEMENTS((*i).kyind)
  n_zind  = N_ELEMENTS((*i).zind)
  n_vars  = N_ELEMENTS((*i).vars)

  time = store_step((*i).time_id,/get,/reg_array)
  tt_data = REFORM(store_step((*i).data_id,/get,/reg_array),[n_kxind,$
    n_kyind,n_zind,n_vars,gui.out.n_spec_sel,series.step_count],/OVERWRITE)

  IF (*i).wsmooth GT 1 THEN $
    tt_data = smooth(tt_data,[1,1,1,1,1,(*i).wsmooth])
    
  csize = 1.0  

 ; build legend_str  
  legend_str = STRARR(n_kxind*n_kyind*n_zind*n_vars)   
  elem=0
  FOR ivar = 0, n_vars-1 DO BEGIN
   FOR ikx = 0, n_kxind-1 DO BEGIN
    FOR iky = 0, n_kyind-1 DO BEGIN
     FOR ins = 0, n_zind-1 DO BEGIN
      legend_str[elem] = ''
      IF (n_vars GT 1) THEN legend_str[elem] = $
         get_var_string((*i).vars[ivar],/fancy)+','
      IF par.x_local THEN legend_str[elem] += $
        '('+rm0es((*series.kx)[(*i).kxind[ikx],(*i).kyind[iky]],prec=3)+$
        ','+rm0es((*series.ky)[(*i).kyind[iky]],prec=3) $
      ELSE legend_str[elem] += $
        '('+rm0es((*series.lxarr)[(*i).kyind[iky]]*(-0.5+(*i).kxind[ikx]/$
        par.nx0),prec=3)+','+rm0es((*series.ky)[(*i).kyind[iky]],prec=3)	
      IF (*i).zind[ins] EQ -1 THEN $
        legend_str[elem] = legend_str[elem] + ')' $
      ELSE legend_str[elem] = legend_str[elem] + ',' + $
        rm0es((*par.z)[(*i).zind[ins]],prec=2)+')'
      elem = elem + 1
     ENDFOR
    ENDFOR
   ENDFOR  
  ENDFOR

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
   sp=(*gui.out.spec_select)[isp]
   
   IF series.step_count GT 1 THEN BEGIN
    set_output, diag, sp, /ps, multi=[0,1,2], charsize=csize
   
;   %%%%%%%%%% Real part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF (n_vars EQ 1) THEN BEGIN
      title    = 'Re [' + get_var_string((*i).vars,/fancy) + ']'
      text_tit = 'Re [' + get_var_string((*i).vars) + ']'
    ENDIF ELSE BEGIN
      title    = 'Re [ ... ]'
      text_tit = 'Re [ ... ]'
    ENDELSE

    timetrace_plot_data, diag, time, TOTAL(DOUBLE(tt_data[*,*,*,*,isp,*]),5), $
      title, legend_str, text_tit

    plot_info_str, diag

;   %%%%%%%%%% Imaginary part %%%%%%%%%%%%%%%%%%%%%%%%%
    IF (n_vars EQ 1) THEN BEGIN
      title    = 'Im [' + get_var_string((*i).vars,/fancy) + ']'
      text_tit = 'Im [' + get_var_string((*i).vars) + ']'
    ENDIF ELSE BEGIN
      title    = 'Im [ ... ]'
      text_tit = 'Im [ ... ]'
    ENDELSE  

    timetrace_plot_data, diag, time, TOTAL(IMAGINARY(tt_data[*,*,*,*,isp,*]),5),$
      title, legend_str, text_tit

    plot_info_str, diag

;   %%%%%%%%%% Absolute value %%%%%%%%%%%%%%%%%%%%%%%%%
    IF (n_vars EQ 1) THEN BEGIN    
      title    = '!9!!!6' + get_var_string((*i).vars,/fancy) + '!9!!!6'
      text_tit = 'Abs [' + get_var_string((*i).vars) + ']'  
    ENDIF ELSE BEGIN
      title    = '!9!!!6 ... !9!!!6'
      text_tit = '| ... |'
    ENDELSE

    timetrace_plot_data, diag, time, TOTAL(ABS(tt_data[*,*,*,*,isp,*]),5), $
      title, legend_str, text_tit

    plot_info_str, diag

    set_output, diag, sp, /reset
   ENDIF ELSE BEGIN
     print, 'timetrace: cannot create plots with one time step'
   ENDELSE
   
;  %%%%%%%%%% DATA output only for first var, kx, ky, z %%%%%%%%%%%%%%%%%%%%%%
   real_temp = REFORM(FLOAT(tt_data[0,0,0,0,isp,*]))
   imag_temp = REFORM(IMAGINARY(tt_data[0,0,0,0,isp,*]))
   set_output, diag, sp, header = ['t','Re['+get_var_string((*i).vars[0])+']',$
     'Im['+get_var_string((*i).vars[0])+']'], dat=[[time],[real_temp],[imag_temp]]
  ENDFOR ;-- isp
END


;######################################################################

PRO gamfunct, X, A, F, pder
  
  f1x = EXP(- A[1] * X)
  f2x = COS( A[2] * X)
  F = A[0] + (1 - A[0]) * f1x * f2x
  
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF N_PARAMS() GE 4 THEN $  
    pder = [[1.-f1x*f2x], [( A[0] - 1.) * X * f1x*f2x], [( A[0] - 1.) * X * SIN( A[2] * X) *f1x ]]  
END  

;######################################################################

PRO gamfunct2, X, A, F, pder
  
  f1x = EXP( - A[1] * X)
  f2x = COS( + A[2] * X)
  f3x = EXP( - A[3] * X)
  F = A[0] + (1 - A[0]) * ( 0.1 * f2x*f3x + 0.9) * f1x
  
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF N_PARAMS() GE 4 THEN $  
    pder = [[1.-( 0.1 * f2x*f3x + 0.9) * f1x], [(A[0]-1.) * X *( 0.1 * f2x*f3x + 0.9) * f1x], $
     [( A[0] - 1.) * 0.1 * X * SIN( A[2] * X) * f3x * f1x],[(A[0]-1.) *0.1*X*f2x*f3x*f1x]]
END
