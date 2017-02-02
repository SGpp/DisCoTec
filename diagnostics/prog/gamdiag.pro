FUNCTION gamdiag_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'GAM diag',$
    help_text : ['Draws a linear-linear and a log-linear plot of the time trace '+$'
    	'of the real, imaginary and absolute value of selected modes and variable. '+$
    	'The corresponding kx, ky and z values can be found in the legend. '+$
	'Data output will only be done for first choice of var, kxind and zind.'],$    
    ext_vars  : [$		 
    	['var','0','index of the fluctuating quantity to be plotted'+$
	 ' (see variable list)'],$
    	['kxind','0','index of kx modes to be plotted (integer '+$
	 'between 0 and nkx0 - 1; default: 1 if kx_center=0 else 0)'],$
    	['zind','0','index of z-value to be plotted (integer '+$
	 'between 0 and nz0 - 1; default: -1 (average))'],$
        ['norm','1','switch on/off normalization to first time step']]}
 		 
END

;######################################################################

PRO gamdiag_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var)   LT 1 THEN *e.var   = 0
  IF N_ELEMENTS(*e.kxind) LT 1 THEN BEGIN
     IF (par.kx_center NE 0.0) THEN $
        *e.kxind = 0 ELSE *e.kxind = 1
  ENDIF
  IF N_ELEMENTS(*e.zind)  LT 1 THEN *e.zind  = -1 ; averaging
  IF N_ELEMENTS(*e.norm) LT 1 THEN *e.norm = 0

  ind = (WHERE(spec[*].charge LT 0))[0]
  Te = (ind GT -1)?spec[ind].temp:1.0
  ind = (WHERE(spec[*].charge GT 0))[0]
  Ti = (ind GT -1)?spec[ind].temp:1.0  

  i = set_internal_vars(diag,{$
    vars    : (*e.var)[0],$
    kxind   : *e.kxind,$
    zind    : *e.zind,$
    norm    : *e.norm,$
    Te      : Te,$
    Ti      : Ti,$
    time_id : PTR_NEW(),$
    data_id : PTR_NEW()})
    
  IF (par.x_local AND par.y_local) THEN fft_format, kxky=(*e.var)[0] $
  ELSE IF par.y_local THEN fft_format, sxky=(*e.var)[0] $
  ELSE fft_format, kxsy=(*e.var)[0]

END

;######################################################################

PRO gamdiag_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  n_kxind = N_ELEMENTS((*i).kxind)
  n_zind = N_ELEMENTS((*i).zind)

  tt_data = DCOMPLEXARR(n_kxind,n_zind,gui.out.n_spec_sel,/NOZERO)
  
  IF (par.x_local AND par.y_local) THEN BEGIN
   FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     FOR ikx = 0, n_kxind - 1 DO BEGIN
        FOR ins = 0, n_zind - 1 DO BEGIN
           IF (*i).zind[ins] EQ - 1 THEN BEGIN
             tt_data[ikx,ins,isp] = COMPLEX(0.0,0.0)
             FOR nsi = 0, par.nz0 - 1 DO BEGIN
               tt_data[ikx,ins,isp] += $
                 (*mom[isp,(*i).vars].kxky)[(*i).kxind[ikx],$
                   0,nsi]*(*series.geom).jac_norm[nsi]
             ENDFOR
             tt_data[ikx,ins,isp] = tt_data[ikx,ins,isp]/par.nz0
           ENDIF ELSE tt_data[ikx,ins,isp] = $
             (*mom[isp,(*i).vars].kxky)[(*i).kxind[ikx],0,(*i).zind[ins]]
         ENDFOR
       ENDFOR
   ENDFOR
  ENDIF ELSE BEGIN
   IF (par.y_local) THEN BEGIN
    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     FOR ikx = 0, n_kxind - 1 DO BEGIN
        FOR ins = 0, n_zind - 1 DO BEGIN
           IF (*i).zind[ins] EQ - 1 THEN BEGIN
             tt_data[ikx,ins,isp] = COMPLEX(0.0,0.0)
             FOR nsi = 0, par.nz0 - 1 DO BEGIN
               tt_data[ikx,ins,isp] += $
                 (*mom[isp,(*i).vars].sxky)[(*i).kxind[ikx],$
                 0,nsi]*(*series.geom).jacobian[ikx,nsi]
             ENDFOR
             tt_data[ikx,ins,isp] = tt_data[ikx,ins,isp]/TOTAL((*series.geom).jacobian[ikx,*])
           ENDIF ELSE tt_data[ikx,ins,isp] = $
             (*mom[isp,(*i).vars].sxky)[(*i).kxind[ikx],0,(*i).zind[ins]]
         ENDFOR
       ENDFOR
    ENDFOR  
   ENDIF ELSE BEGIN
    tt_data = COMPLEX(0.0,0.0)
    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     FOR ikx = 0, n_kxind - 1 DO BEGIN
        FOR ins = 0, n_zind - 1 DO BEGIN
           IF (*i).zind[ins] EQ -1 THEN BEGIN
              tt_data[ikx,ins,isp] += $
                 TOTAL((*mom[isp,(*i).vars].kxsy)[(*i).kxind[ikx],*,*]*$
                       (*series.geom).jacobian)/TOTAL((*series.geom).jacobian)
           ENDIF ELSE BEGIN
              FOR iy = 0, par.ny0 - 1 DO BEGIN
                 tt_data[ikx,ins,isp] += $
                    (*mom[isp,(*i).vars].kxsy)[(*i).kxind[ikx],iy,(*i).zind[ins]]/par.ny0
              ENDFOR
           ENDELSE
        ENDFOR
       ENDFOR
    ENDFOR  
   ENDELSE
  ENDELSE

  (*i).time_id = store_step((*i).time_id,mom_time)
  (*i).data_id = store_step((*i).data_id,tt_data)

END

;######################################################################
PRO gamdiag_plot_residual, h, A_R

  maxval = MAX(A_R, min=minval)
  maxxval = MAX(h)
    
; lin - lin plot
  PLOT, h, A_R, COLOR = 1, $
    /XSTYLE, YRANGE=[0,maxval], $
    XRANGE=[0,1.1*maxxval],PSYM=4,$
    TITLE=title, CHARSIZE=csize,$
    xtitle='h', ytitle='A!IR!N',$
    POSITION=[0.15,0.2,0.95,0.9615]

  h = INDGEN(20)/20.*1.1*maxxval
  A_R = 1./(1.+1.6/h)
  OPLOT, h, A_R, COLOR=2
      
END                          ; plot residual level vs. h=q^2/sqrt(r/R)
;######################################################################

PRO gamdiag_plot_result_vs_x, xaxis, fit_res, theo_res

  title = ['A!IR!N','!7c!6!IG!N','!7x!6!IG!N']

  FOR p=0, N_ELEMENTS(fit_res[*,0])-1 DO BEGIN
      
      maxval = MAX([fit_res[p,*],theo_res[p,*]], min=minval)
      
      ; lin - lin plot
      PLOT, xaxis, fit_res[p,*], COLOR = 1, $
        /XSTYLE, YRANGE=[0,maxval], PSYM=4,$
        YTITLE=title[p], CHARSIZE=csize,$
        xtitle='x/a', title='',$
        POSITION=[0.15,0.2,0.95,0.9615]
      
      OPLOT, xaxis,theo_res[p,*], COLOR=2
      
      plot_legend, [1,2],['GENE','theory'], per_line=2,$
        x_offset=[0.15,0.0,0.95,0.12], psym=[4,0]
      
  ENDFOR
  
END

;###########################################################################

PRO gamdiag_plot_data, diag, time, tt_data, title, legend_str, text_tit,$
                       A_Sugama=A_Sugama, A_GaoZhou=A_GaoZhou, B=B

  COMMON global_vars
  i = (*diag).internal_vars

  n_kxind = N_ELEMENTS(tt_data[*,0,0])
  n_zind  = N_ELEMENTS(tt_data[0,*,0])
  n_tind  = N_ELEMENTS(tt_data[0,0,*])

  maxval = MAX(tt_data, min=minval)
  show_theory = (n_kxind*n_zind EQ 1)
    
; lin - lin plot
  PLOT, time, tt_data[0,0,*], COLOR = 1, $
    /XSTYLE, YRANGE=[minval,maxval], $
    TITLE=title,$ ;CHARSIZE=csize,$
    POSITION=[0.15,0.2,0.95,0.9615]
      
  count = 0

  FOR ikx=0, n_kxind-1 DO BEGIN
    FOR ins=0, n_zind-1 DO BEGIN   
     OPLOT, time, tt_data[ikx,ins,*], COLOR = 1 + 3*count
     IF (show_theory) THEN BEGIN
        IF KEYWORD_SET(A_Sugama) THEN BEGIN
           gamfunct,time,A_Sugama[*,count],F
           OPLOT, time, F, COLOR = 2 + 3*count, linestyle=2
        ENDIF
        IF KEYWORD_SET(A_GaoZhou) THEN BEGIN
           gamfunct,time,A_GaoZhou[*,count],F
           OPLOT, time, F, COLOR = 3 + 3*count, linestyle=2
        ENDIF
        IF KEYWORD_SET(B) THEN BEGIN
           gamfunct_damp,time,B[*,count],F
           OPLOT, time, F, COLOR = 4 + 3*count, linestyle=2     
        ENDIF
     ENDIF
     
;     omega_f = sqrt2/par.major_R
;     gamma_f = omega_f*exp(-q_^2)
;     OPLOT,time,(1-A_R)*exp(-gamma_f*time)*cos(omega_f*time)+A_R, $
;           COLOR = 5 + 4*count,linestyle=2

     count = count + 1 	
   ENDFOR   
  ENDFOR
  
  ; add base line  

  IF !Y.CRANGE[0] LE 0 AND !Y.CRANGE[1] GE 0 THEN $
    oplot, !X.CRANGE, [0,0], color=1

  ; print legend
  IF (show_theory) THEN BEGIN
     plot_legend, [1,4,2,3],['Sim.','Fit', '[Sugama06]','[Gao10Zhou11]'], $
      x_offset=[0.15,0.0,0.95,0.12], per_line = 2
  ENDIF ELSE plot_legend, INDGEN(n_kxind*n_zind) + 1, $
      legend_str, x_offset=[0.15,0.0,0.95,0.12]
  plot_info_str, diag

END ; plot time_trace

;######################################################################

PRO gamdiag_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF series.step_count LT 2 THEN RETURN
  
  n_kxind = N_ELEMENTS((*i).kxind)
  n_zind = N_ELEMENTS((*i).zind)
  
  time = store_step((*i).time_id,/get,/reg_array)
  full_data = REFORM(store_step((*i).data_id,/get,/reg_array),[n_kxind,$
    n_zind,gui.out.n_spec_sel,series.step_count],/OVERWRITE)

  tt_data = DOUBLE(full_data)

  csize = 1.42

 ; build legend_str  
  legend_str = STRARR(n_kxind*n_zind)   
  
  fit_res = DBLARR(4,n_kxind*n_zind,/NOZERO)
  Sugama_res = DBLARR(3,n_kxind*n_zind,/NOZERO)
  GaoZhou_res = DBLARR(3,n_kxind*n_zind,/NOZERO)
  coord_arr = DBLARR(2,n_kxind*n_zind,/NOZERO)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     sp=(*gui.out.spec_select)[isp]

    count = 0
    FOR ikx = 0, n_kxind-1 DO BEGIN
      FOR ins = 0, n_zind-1 DO BEGIN

      ;normalize to t=0 value
      IF (*i).norm THEN $
         tt_data[ikx,ins,isp,*] /= tt_data[ikx,ins,isp,0]

      ;construct some values
      IF (par.x_local) THEN BEGIN
         x_R = par.trpeps
         q_ = par.q0
         kx=(*series.kx)[(*i).kxind[ikx],0]
      ENDIF ELSE BEGIN
         x_R = (0.5 + (-series.lx/2.+(*i).kxind[ikx]*series.dx)*par.rhostar)*$
               par.minor_r/par.major_R
         q_ = (*series.geom).q((*i).kxind[ikx])
         kx=2.0*!PI/par.lx
      ENDELSE

      ;get analytical result
      gamdiag_get_sugama_res,q_,x_R*par.major_R,par.major_R,kx,(*i).ti,(*i).te,A_Sugama
      gamdiag_get_GaoZhou_res,q_,x_R*par.major_R,par.major_R,kx,(*i).ti,(*i).te,A_GaoZhou

      ;use analytic result as initial cond. for fit
      B = [A_GaoZhou,0.0]

      IF (series.step_count GT 5) THEN BEGIN
         weights = FLTARR(series.step_count)+1.0
         yfit = CURVEFIT(time,REFORM(tt_data[ikx,ins,isp,*]),weights,$
                         B, SIGMA, FUNCTION_NAME='gamfunct_damp')
      ENDIF

      ; print comparison
      print, 'x/R=',rm0es(x_R), ', q=',rm0es(q_)
      gamdiag_print_result, A_Sugama, A_GaoZhou, B, TOTAL(tt_data[ikx,ins,isp,*])/series.step_count,$
        MEDIAN(tt_data[ikx,ins,isp,*])

      ;store fit results
      fit_res[*,count] = B
      coord_arr[0,count] = q_
      coord_arr[1,count] = x_R
      Sugama_res[*,count] = A_Sugama
      GaoZhou_res[*,count] = A_GaoZhou

      IF par.x_local THEN legend_str[count] = $
        '('+rm0es((*series.kx)[(*i).kxind[ikx],0],prec=3) $
      ELSE legend_str[count] = $
        '('+rm0es((*series.lxarr)[0]*(-0.5+(*i).kxind[ikx]/$
         par.nx0),prec=3)
	
      IF (*i).zind[ins] EQ -1 THEN $
        legend_str[count] = legend_str[count] + ')' $
      ELSE legend_str[count] = legend_str[count] + ',' + $
        rm0es((*par.z)[(*i).zind[ins]],prec=2)+')'

      count += 1
   ENDFOR  ;ins
  ENDFOR  ;ikx
     
  set_output, diag, sp, /ps, multi=[0,1,1], charsize=csize
   
;   %%%%%%%%%% Real part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  title = get_var_string((*i).vars,/fancy)+'(t)'
  text_tit = get_var_string((*i).vars)+'(t)'
  IF (*i).norm THEN BEGIN
     title += '/'+get_var_string((*i).vars,/fancy)+'(t=0)' 
     text_tit += '/'+get_var_string((*i).vars)+'(t=0)' 
  ENDIF
  
  gamdiag_plot_data, diag, time, TOTAL(tt_data[*,*,isp,*],3),$
    title, legend_str, text_tit,$
    A_Sugama=Sugama_res,A_GaoZhou=GaoZhou_res, B=fit_res

  IF (n_kxind*n_zind GT 1) THEN BEGIN
     gamdiag_plot_residual, (sqrt(coord_arr[1,*])/$
       coord_arr[0,*]^2), fit_res[0,*]
     gamdiag_plot_result_vs_x, coord_arr[1,*]*par.major_R/par.minor_r,$
       fit_res, Sugama_res
      
  ENDIF
  set_output, diag, sp, /reset

; %%%%%%%%%% DATA output only for one set of var, kx, ky, z %%%%%%%%%%%%%%%%%%%%%%
  xind = n_kxind/2
  real_temp = REFORM(DOUBLE(full_data[xind,0,isp,*]));/$
;   MAX(ABS(full_data[xind,0,isp,*]))
  imag_temp = REFORM(IMAGINARY(full_data[xind,0,isp,*]));/$
;   MAX(ABS(full_data[xind,0,isp,*]))

  set_output, diag, sp, header = ['t','Re['+get_var_string((*i).vars[0])+']',$
    'Im['+get_var_string((*i).vars[0])+']'], dat=[[time],[real_temp],[imag_temp]],$
    commentline=['xind= '+rm0es(xind)+', x/a= '+rm0es(coord_arr[1,xind]*par.major_R/par.minor_r)+$
    ', q= '+rm0es(coord_arr[0,xind]),$
    'Fit:    '+rm0es(fit_res[0,xind])+','+rm0es(fit_res[1,xind])+','$
     +rm0es(fit_res[2,xind]),$
    'Analytic (Gao, Zhou/Yu): '+rm0es(GaoZhou_res[0,xind])+','+rm0es(GaoZhou_res[1,xind])+','$
     +rm0es(GaoZhou_res[2,xind])]

  IF (n_kxind*n_zind GT 1) THEN BEGIN
     set_output, diag, sp, header = ['x/a','q','A_R','gamma','omega','A_R(theo)',$
       'gam(theo)','om(theo)'], dat=[[REFORM(coord_arr[1,*])*par.major_R/par.minor_r],$
       [REFORM(coord_arr[0,*])],$                              
       [REFORM(fit_res[0,*])],[REFORM(fit_res[1,*])],[REFORM(fit_res[2,*])],$
       [REFORM(GaoZhou_res[0,*])],[REFORM(GaoZhou_res[1,*])],[REFORM(GaoZhou_res[2,*])]],/append
  ENDIF

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

PRO gamfunct_damp, X, A, F, pder
  
  f1x = EXP(- A[1] * X)
  f2x = COS( A[2] * X)
  f3x = A[0]-A[3]*X
  F = f3x + (1 - f3x) * f1x * f2x
  
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF N_PARAMS() GE 4 THEN $  
    pder = [[1.-f1x*f2x], [( f3x - 1.) * X * f1x*f2x], [( f3x - 1.) * X * SIN( A[2] * X) *f1x ], [-X*(1.-f1x*f2x)]]  
END  

;######################################################################

PRO gamdiag_get_sugama_res, q_, x, R0, kx, Ti, Te, res
; see J. Plasma Phys, 72 (2006), pp.825-828
; and Erratum: J. Plasma Phys. 
; input: q, x, R0, kx, Ti, Te
; output: array res=[residual, gamma, omega]

  erratum = 0  ;the non-corrected version seems to fit better???
  L_par = q_*R0
  sqrt2 = sqrt(2.0)
  TeTi = Te/Ti

  ;renormalize from rho_ref to v_Ti/Omega*q
  khat = kx*q_*sqrt(2.0*Ti)

  help_fac = 1 + 2.*(23.+16.*TeTi+4.*TeTi^2)/(q_*(7+4*TeTi))^2
  omega_hat = 0.5*sqrt(7.+4.*TeTi)*q_*sqrt(help_fac)

  if (erratum) then begin
     gamma_hat = sqrt(!PI)/2.*q_^2/help_fac*$
               (exp(-omega_hat^2)*omega_hat^2*(omega_hat^2+1+2*TeTi)+$
                0.25*(sqrt2*kx*q_*sqrt(TeTi))^2*exp(-0.25*omega_hat^2)*omega_hat^2*$
                (omega_hat^4/128.+(1.+TeTi)/16.*omega_hat^2+$
                 (3./8.+7./16.*TeTi+5./32.*TeTi^2)))
;; maybe, the erratum is missing a factor of 2.0 in the last major bracket??
  endif else begin
     gamma_hat = sqrt(!PI)/2.*q_^2/help_fac*$
               (exp(-omega_hat^2)*omega_hat^2*(omega_hat^2+1+2*TeTi)+$
                0.25*(khat)^2*exp(-0.25*omega_hat^2)*omega_hat^2*$
                (omega_hat^4/64.+(1.+3.*TeTi/8.)*(omega_hat^2/8.+3./4.)))
  endelse

  ;translate to from v_Ti/L_|| to c_ref/Lref
  omega_f = omega_hat*sqrt(2.0*Ti)/L_par
  gamma_f = gamma_hat*sqrt(2.0*Ti)/L_par

  A_R = 1./(1.+1.6*q_^2/sqrt(x/R0))

  res=[A_R,gamma_f,omega_f]

END

;######################################################################

PRO gamdiag_get_GaoZhou_res, q, x, R0, kx, Ti, Te, res
; should be more accurate for non-circular geometries
; input: q, x, R0, kx, Ti, Te
; output: array res=[residual, gamma, omega]
; gamma and omega are based on a prediction by [Z. Gao, Phys. Plasmas 17, 092503 (2010)];
; the residual level is taken from [D. Zhou and W. Yu, Phys. Plasmas 18, 052505 (2011)]

  COMMON global_vars

  k = par.kappa
  s_k = par.s_kappa
  eps = x/R0 ;par.trpeps
  drR = par.drR
  tau = Te/Ti


  ;normalize from rho_s to v_Ti/Omega
  khat = kx * sqrt(2.0*Ti)

  Qterm = (7./4.+tau)^(-2)*(k^2+1.)/2.*(747./32.+481./32.*tau+35./8.*tau^2+0.5*tau^3)$
          -(7./4.+tau)^(-1)*((39.+13.*k^2)/16.+(25.-k^2)/8.*tau+(5.-k^2)/4.*tau^2)$
          +((9.+6.*k^2+9*k^4)/(16.*k^2*(k^2+1.))+(k^2-1)^2/(4*k^2*(k^2+1))*tau)

  omega_hat = sqrt((7./4.+tau)*(2./(k^2+1)))*$
          (1.-s_k*k^2/(4*k^2+4)-eps^2*(9*k^2+3.)/(8*k^2+8.)$
           -drR^2*k^2/(4*k^2+4.)+eps*drR*(4*k^2+1.)/(4*k^2+4))$
          *(1+(k^2+1.)/(4.*q^2)*(23./8.+2*tau+0.5*tau^2)$
            *(7./4.+tau)^(-2)+khat^2/(4.*k^2)*Qterm)

  gamma_hat = 4.*k^2*sqrt(7./4.+tau)/(khat^2*(k^2+1.)^1.5)*$
              (1.+(2*k^2+5.)/(4*k^2+4.)*s_k-(27*k^2+9.)/(8*k^2+8.)*eps^2 $
               -(7.*k^2+4)/(4*k^2+4.)*drR^2+9*k^2/(4*k^2+4.)*eps*drR)$
              *EXP(-sqrt(7./4.+tau)/khat * sqrt(2*k^2/(k^2+1))*(1.+(3*k^2+4.)/(4*k^2+4.)*s_k $
               -(9*k^2+3.)/(8*k^2+8.)*eps^2-(3*k^2+2.)/(4*k^2+4.)*drR^2+$
               (4*k^2+1.)/(4*k^2+4.)*eps*drR))

  ;normalize from v_Ti/R0 to cref/Lref
  gamma_f = gamma_hat * sqrt(2.0*Ti) / par.major_R
  omega_f = omega_hat * sqrt(2.0*Ti) / par.major_R 


  ;[D. Zhou and W. Yu, Phys. Plasmas 18, 052505 (2011)]
  delt = par.delta ; = - COS(!PI/2.+ASIN(par.delta))
  nth = 100
  th = FINDGEN(nth)*2.0*!PI/(nth-1)

;  check equivalence of flux surface as Zhou&Yu use a different R,Z definition
;  R = R0*par.Lref + eps/par.Lref * (cos(th)-delt*sin(th)^2)
;  Z = k * eps/par.Lref * sin(th)
;  PLOT, (*series.geom).R[*], (*series.geom).Z[*], $
;      /ISOTROPIC, COLOR=1, /XSTYLE, /YSTYLE, XTITLE='!6R', YTITLE='!6Z',$
;      YMARGIN=[4,4],XMARGIN=[8,2]
;  OPLOT, R, Z, COLOR=2, THICK=2

  grhrh = (k^2*cos(th)^2+sin(th)^2+2*delt*sin(th)*sin(2*th)+(delt*sin(2*th))^2)/$
        (k^2*(1+delt*sin(th)^2*cos(th))^2)

  I0 = INT_TABULATED(th,((1.+delt*sin(th)^2*cos(th))*grhrh))/(2.0*!PI)
  I1 = INT_TABULATED(th,(1.+delt*sin(th)^2*cos(th))*(cos(th)-delt*sin(th)^2)*grhrh)/(2.0*!PI)

  Sfac = ((25./16.-53./256.*delt)+0.5*sqrt(eps)-$
          (3./64.-93./256.*delt+(9*k^2)/(8.*q^2)*(0.75*I0+I1))*eps)/$
         (k^2*(1.+3./8.*delt*eps)^2*I0)

  A_R = 1./(1.+Sfac*q^2/sqrt(eps))

  res=[A_R,gamma_f,omega_f]

END


;######################################################################

PRO gamdiag_print_result, A1, A2, B, avg, median
    print, '        Fit      Sugama   Gao-Zhou Gao-Zhou/Fit Average  Median'
    print, 'A_R     ',STRING(B[0],FORMAT='(F7.4,5X)'),'  ',$
           STRING(A1[0],FORMAT='(F7.4)'),'  ',$
           STRING(A2[0],FORMAT='(F7.4)'),'  ',$
           STRING(A2[0]/B[0],FORMAT='(F7.4)')+'      ',$
           STRING(avg,FORMAT='(F7.4)')+'  ',$
           STRING(median,FORMAT='(F7.4)')
    print, 'gamma g ', STRING(B[1],FORMAT='(F7.4)'),'  ',$
           STRING(A1[1],FORMAT='(F7.4)'),'  ',$
           STRING(A2[1],FORMAT='(F7.4)'),'  ',$           
           STRING(A2[1]/B[1],FORMAT='(F7.4)')
    print, 'omega o ',STRING(B[2],FORMAT='(F7.4)'),'  ',$
           STRING(A1[2],FORMAT='(F7.4)'),'  ',$
           STRING(A2[2],FORMAT='(F7.4)'),'  ',$
           STRING(A2[2]/B[2],FORMAT='(F7.4)')
    IF N_ELEMENTS(B) GT 3 THEN BEGIN
       print, 'damp. d ', STRING(B[3],FORMAT='(F7.4)'),'  ',$
              STRING(0.0,FORMAT='(F7.4)'),'  ',$
              STRING(0.0,FORMAT='(F7.4)')       
    ENDIF
    print, 'fit function: f(t) = (A_R-d t) + (1-(A_R-d t))*COS(o t)*EXP(-g t)'
    print, 'Sugama:   [H. Sugama and T.-H. Watanabe, J. Plasma Phys, 72 (2006), pp.825-828]
    print, 'Gao-Zhou: [Z. Gao, Phys. Plasmas 17, 092504 (2010)] for gamma/omega'
    print, '          [D. Zhou and W. Yu, Phys. Plasmas 18, 052505 (2011)] for residual'
END
