FUNCTION fluxprof2D_info

  RETURN, {$
    type      : 'mom global',$
    title     : 'Flux Profile (2D)',$
    help_text : ['Displays radial fluxes vs. x and t'],$
    ext_vars  : [['xind','0','indices of x grid points to be plotted; '+$
                  'default: 0:nx0-1'],$
                ['momentum','1','add parallel momentum fluxes; default:off'],$
                ['with_area','1','also display the fluxes multiplied by the '+$
                  'flux-surface area (for nonlocal simulations); default:off']]}

END

;######################################################################

PRO fluxprof2D_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = INDGEN(par.nx0)
  nxind  = N_ELEMENTS(*e.xind)

  xgrid = (- series.lx / 2.0 + (*e.xind) * series.dx)

  IF NOT KEYWORD_SET(*e.momentum) THEN *e.momentum = 0
  IF NOT KEYWORD_SET(*e.with_area) THEN *e.with_area = 0
  IF (par.x_local) THEN *e.with_area = 0

  ;comment the following two lines to switch the axis to rho_ref units
  IF (NOT par.x_local) THEN $
    xgrid = xgrid*par.rhostar + par.x0

  nf = par.n_fields
  with_em = (nf GT 1) AND (par.beta GT 1e-5)
  IF with_em THEN BEGIN
    ; fluxes: G_es, Q_es, G_em, Q_em[, P_es, P_em]
    n_fluxes = (4 + 2 * (*e.momentum)) * (par.n_moms / 6)
    vars = INDGEN(nf+par.n_moms)
  ENDIF ELSE BEGIN
    ; G_es, Q_es[,P_es] (passing, trapped, FLR)
    n_fluxes = (2 + (*e.momentum)) * (par.n_moms / 6)
    IF (par.n_moms EQ 6) THEN BEGIN
       vars = [0,nf,nf+1,nf+2]
       IF *e.momentum THEN vars = [vars,nf+5]
    ENDIF ELSE BEGIN
      vars = 0
      FOR j = 0, par.n_moms / 6 - 1 DO vars = *e.momentum ? $
        [vars,nf+6*j,nf+1+6*j,nf+2+6*j,nf+5+6*j] : $
        [vars,nf+6*j,nf+1+6*j,nf+2+6*j]
    ENDELSE
  ENDELSE

  IF (par.x_local) THEN BEGIN
      jac_fac = DBLARR(par.nx0,par.nz0,/NOZERO)
      FOR i=0, par.nx0 - 1 DO BEGIN
          jac_fac[i,*] = (*series.geom).jacobian
      ENDFOR
      jac_fac /= TOTAL((*series.geom).jacobian)
; the following doesn't work for unkown reason in two spec cases
;        jac_fac = REBIN((*series.geom).jacobian,$
;       [par.nx0,par.nz0])/TOTAL((*series.geom).jacobian)
      area = TOTAL((*series.geom).jacobian*sqrt((*series.geom).gxx))/par.nz0
   ENDIF ELSE BEGIN
      jac_fac = (*series.geom).jacobian * $
                REBIN(1.0/TOTAL((*series.geom).jacobian,2),[par.nx0,par.nz0])
      area = total((*series.geom).jacobian[*,*]*$
                   sqrt((*series.geom).gxx[*,*]),2)/par.nz0
  ENDELSE
  area *= (2.0*!PI)^2*par.n_pol*(*series.geom).C_y

  ; for comparison with nrg (volume average):
  ;*e.xind = INDGEN(par.nx0)
  ;nxind = par.nx0
  ;jac_fac = (*series.geom).jacobian / TOTAL((*series.geom).jacobian) * nxind

  i = set_internal_vars(diag,{$
    xind        : *e.xind,$
    xgrid       : xgrid,$
    n_fluxes    : n_fluxes,$
    momentum    : *e.momentum,$
    with_area   : *e.with_area,$
    jac_fac     : jac_fac,$
    area        : series.ly*par.rhostar*par.minor_r*2.0*!PI*par.n0_global*$
                    TOTAL((*series.geom).jacobian[*,*]*$
                          sqrt((*series.geom).gxx[*,*]),$
                          2-par.x_local)/par.nz0,$
    profx       : FLTARR(nxind,n_fluxes,/NOZERO),$
    tmpprofx    : FLTARR(nxind,par.nz0,/NOZERO),$
    time_id     : PTR_NEW(),$
    profx_id    : PTR_NEW()})

  fft_format, sxky=vars

END

;######################################################################

FUNCTION fluxprof2D_calc_x_prof, diag, var1, var2

  COMMON global_vars

  i = (*diag).internal_vars

  nxind = N_ELEMENTS((*i).xind)

  (*i).tmpprofx = $
    REBIN(FLOAT(var1[(*i).xind,0,*]*CONJ(var2[(*i).xind,0,*])),$
    [nxind,par.nz0]) + 2.0 * TOTAL(FLOAT(CONJ($
    var1[(*i).xind,1:par.nky0-1,*])*var2[(*i).xind,1:par.nky0-1,*]),2)

  (*i).tmpprofx *= (*i).jac_fac[(*i).xind,*]

  RETURN, TOTAL((*i).tmpprofx,2)

END

;######################################################################


PRO fluxprof2D_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nxind = N_ELEMENTS((*i).xind)

  tmp_profx = FLTARR(nxind,(*i).n_fluxes,gui.out.n_spec_sel,/NOZERO)

  nf = par.n_fields
  rho_ref = SQRT(series.mref*series.Tref) / (series.Qref * series.Bref)
  with_em = (nf GT 1) AND (par.beta GT 1e-5)

  ; Calculate vx, Bx, grad q_par
  ; ve_x = - dphi / dy , B_x = dApar / dy
  IF (par.x_local) THEN BEGIN
    sqrtgxx_inv = REBIN(REFORM(1.0/sqrt((*series.geom).gxx),[1,1,par.nz0]),$
                   [par.nx0,par.nky0,par.nz0])
    ve_x = - COMPLEX(0,1) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0,par.nky0,par.nz0]) * (*mom[0,0].sxky) / $
      (series.Bref*(*series.geom).C_xy)* $
      sqrtgxx_inv
    IF with_em THEN B_x = COMPLEX(0.0,1.0) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0,par.nky0,par.nz0]) * (*mom[0,1].sxky) / $
      (series.Bref*(*series.geom).C_xy) * $
      sqrtgxx_inv
  ENDIF ELSE BEGIN
    sqrtgxx_inv = REBIN(REFORM(1.0/sqrt((*series.geom).gxx[*,*]),$
                               [par.nx0,1,par.nz0]),[par.nx0,par.nky0,par.nz0])
    ve_x = - COMPLEX(0,1) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0,par.nky0,par.nz0]) * (*mom[0,0].sxky) / $
      REBIN(series.Bref*(*series.geom).C_xy,[par.nx0,par.nky0,par.nz0])* $
      sqrtgxx_inv
    IF with_em THEN B_x = COMPLEX(0.0,1.0) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0,par.nky0,par.nz0]) * (*mom[0,1].sxky) / $
      REBIN(series.Bref*(*series.geom).C_xy,[par.nx0,par.nky0,par.nz0])* $
      sqrtgxx_inv
  ENDELSE

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]
    n_notrap = 2 + 2 * with_em ; number of fluxes without trapdiag

    FOR j = 0, par.n_moms / 6 - 1 DO BEGIN ; passing, trapped, FLR
      ; Gamma_es = <n * ve_x>
      (*i).profx[*,n_notrap*j+0] = fluxprof2D_calc_x_prof(diag,$
        *mom[isp,nf+6*j].sxky,ve_x) * spec[sp].dens

      ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
      IF (par.x_local) THEN BEGIN
        var1 = (0.5 * (*mom[isp,nf+6*j+1].sxky) + (*mom[isp,nf+6*j+2].sxky))/series.nref + $
	  1.5 * (*mom[isp,nf+6*j].sxky)/(series.Tref*series.Qref)
      ENDIF ELSE BEGIN
         var1 = (0.5 * (*mom[isp,nf+6*j+1].sxky) + (*mom[isp,nf+6*j+2].sxky)) * $
          REBIN((*spec[sp].prof).dens/series.nref,[par.nx0,par.nky0,par.nz0]) + $
	  1.5 * (*mom[isp,nf+6*j].sxky) * REBIN((*spec[sp].prof).temp/$
          (series.Tref*series.Qref),[par.nx0,par.nky0,par.nz0])         
      ENDELSE

      (*i).profx[*,n_notrap*j+1] = fluxprof2D_calc_x_prof(diag,var1,ve_x) * $
        (spec[sp].dens * spec[sp].temp)


      IF (*i).momentum THEN BEGIN ; P_es = <upar ve_x>
         (*i).profx[*,n_notrap*j+2*(1+with_em)] = fluxprof2D_calc_x_prof(diag,$
           *mom[isp,nf+6*j+5].sxky,ve_x) * spec[sp].dens * spec[sp].mass
      ENDIF


      IF with_em THEN BEGIN
        ; Gamma_em = <upar * B_x>
        IF (NOT par.x_local) THEN BEGIN
          var1 = (*mom[isp,nf+6*j+5].sxky) * $
            REBIN((*spec[sp].prof).dens,[par.nx0,par.nky0,par.nz0])
        ENDIF ELSE var1 = *mom[isp,nf+6*j+5].sxky
        (*i).profx[*,n_notrap*j+2] = fluxprof2D_calc_x_prof(diag,$
           var1,B_x)

        ; Q_em = (qpar + qperp + 5/2 upar) B_x
        ; mom(5) und mom(6) are NOT identical to qpar and qperp,
        ; but rather to qpar + 3/2 upar and qperp + upar, resp.
        (*i).profx[*,n_notrap*j+3] = fluxprof2D_calc_x_prof(diag, $1
             ((*mom[isp,nf+6*j+3].sxky)+(*mom[isp,nf+6*j+4].sxky)), B_x)*$
              spec[sp].dens*spec[sp].temp

        ; P_em = (Tpar + n) B_x
        IF (*i).momentum THEN BEGIN 
           (*i).profx[*,n_notrap*j+5] = fluxprof2D_calc_x_prof(diag,$
             ((*mom[isp,nf+6*j+1].sxky)*series.nref +$
              (*mom[isp,nf+6*j].sxky))*series.Tref, B_x)*$
              (spec[sp].dens *spec[sp].temp)
        ENDIF
     ENDIF
     ENDFOR

     tmp_profx[0,0,isp] = (*i).profx
   ENDFOR ; species loop

;print, 'Q = ', TOTAL((*i).profx[*,1])/par.nx0, MIN((*i).profx[*,1]), MAX((*i).profx[*,1])
   (*i).time_id = store_step((*i).time_id,mom_time)
   (*i).profx_id = store_step((*i).profx_id,tmp_profx)

END

;######################################################################

PRO fluxprof2D_plot2D, data, xaxis, yaxis, xtitle=xtitle, $
  ytitle=ytitle, ztitle=ztitle, title=title

  IF NOT KEYWORD_SET(xtitle) THEN xtitle = ''
  IF NOT KEYWORD_SET(ytitle) THEN ytitle = ''
  IF NOT KEYWORD_SET(ztitle) THEN ztitle = ''
  IF NOT KEYWORD_SET(title) THEN title=''

  LOADCT, 33, FILE='internal/colortable.tbl'

  log=0

  c_levels = 128
  max = MAX(data,min=min)
  lev_col = contour_levels([min,max],c_levels,/no_fix_zero,log=log)

  csize = 1.4
  cthick = 3.0

  CONTOUR, data, xaxis, yaxis, $
    LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL, $
    /NORMAL, /XSTYLE, /YSTYLE, XTITLE=xtitle, YTITLE=ytitle, $
    YMARGIN=[4,6], CHARTHICK=cthick

  plot_colorbar, lev_col, POSITION=[0.165,0.85,0.95,0.925], $
    CHARSIZE=0.75*csize, log=log
  xpos = 0.5-!D.X_CH_SIZE*0.5*STRLEN(title)/!D.X_SIZE
  XYOUTS, xpos, 0.95, title, COLOR=1, /NORMAL, $
    CHARSIZE=2*csize, CHARTHICK=cthick

END

;######################################################################

PRO fluxprof2D_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  n_xind = N_ELEMENTS((*i).xind)

  eps = 0 ;1

  time = store_step((*i).time_id,/get,/reg_array)
  profx = REFORM(store_step((*i).profx_id,/get,/reg_array),$
    [n_xind,(*i).n_fluxes,gui.out.n_spec_sel,series.step_count],/OVERWRITE)

  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  n_notrap = (2 + (*i).momentum) * (1 + with_em)

  var_names = STRARR(n_notrap)
  var_names2 = STRARR(n_notrap)
  var_unit = STRARR(n_notrap)
  var_names[0:1] = ['!7C!6!Ies','!6Q!Ies']
  var_names2[0:1] = ['G_es','Q_es']
  var_unit[0] = ['!7C!6!IGB!N']
  var_unit[1] = ['!6Q!IGB!N']
  IF (*i).momentum THEN BEGIN
    var_names[2*(1+with_em)] = '!7P!6!Des'
    var_names2[2*(1+with_em)] = 'P_es'
    var_unit[2*(1+with_em)] = ['!7P!6!IGB!N']
  ENDIF

  IF with_em THEN BEGIN
    var_names[2:3]=['!7C!6!Iem','!6Q!Iem']
    var_names2[2:3]=['G_em','Q_em']
    var_unit[2] = ['!7C!6!IGB!N']
    var_unit[3] = ['!6Q!IGB!N']
    IF (*i).momentum THEN BEGIN
      var_names[5] = '!7P!6!Dem'
      var_names2[5] = 'P_em'
      var_unit[5] = ['!7P!6!IGB!N']
    ENDIF
  ENDIF
    
  nbflux = (2 + (*i).momentum) * (1 + with_em)

  rho_str = get_var_string(/rhostr)
  IF (par.x_local) THEN xtitle = '!6x / '+rho_str $
  ELSE xtitle = '!6x / a'
  
  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN

     sp=(*gui.out.spec_select)[isp]

     set_output, diag, sp, header = ['x/a'],$
                 dat=[[(*i).xgrid]]
     set_output, diag, sp, header = ['t'],$
                 dat=[[time]],/append

     set_output, diag, sp, ps=(eps EQ 0), eps=eps,multi=[0,1,1], $
                 ysize=14.4, xsize=20
     
     FOR iflux = 0, nbflux-1 DO BEGIN
        
        timeline = '!6t='+rm0es(gui.out.start_t)+'-'+$
                   rm0es(gui.out.end_t)+' '+get_var_string(1,/time,/ounit)+$
                   ', ' + STRTRIM(STRING(series.step_count),2) + $
                   ', ' + spec[sp].name

        title = var_names[iflux]+','+spec[sp].name+'!N / '+$
                var_unit[iflux]

        data = reform(profx[*,iflux,isp,*])

        IF (*i).n_fluxes EQ n_notrap THEN BEGIN 
           fluxprof2D_plot2D, data, (*i).xgrid, time,$
                              xtitle=xtitle, ytitle='t / ['+get_var_string(0,/time,/ounit)+']',$
                              title=title

           IF ((*i).with_area) THEN BEGIN
              FOR ix=0, N_ELEMENTS(data[*,0])-1 DO $
                 data[ix,*] = data[ix,*]*(*i).area[(*i).xind[ix]]
              
              fluxprof2D_plot2D, data, (*i).xgrid, time,$
                                 xtitle=xtitle, ytitle='t / ['+get_var_string(0,/time,/ounit)+']',$
                                 title=title+' * A!IFS!N / '+ series.Lref_str + '!E2!N'
           ENDIF

        ENDIF ELSE BEGIN
           printerror, 'trapdiag currently not supported'
        ENDELSE 


        set_output, diag, sp, header = [var_names2[iflux]+','+spec[sp].name],$
                    dat=[[data]],/append

     ENDFOR                     ;-- iflux

  ENDFOR                        ;-- isp

  set_output, diag, sp, /reset


END
