FUNCTION fluxgrad_info

  RETURN, {$
    type      : 'mom global',$
    title     : 'Flux vs. gradient (nonlocal)',$
    help_text : ['to be filled'],$
    ext_vars  : [['xind','0','indices of x grid points for averaging '+$
                 '(integer between 0 and nx0 - 1; default: 2*nx0/5:3*nx0/5)']]}

END

;######################################################################

PRO fluxgrad_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.xind) LT 1 THEN $
    *e.xind = 2 * par.nx0 / 5 + INDGEN(par.nx0/5+1)
  nx = N_ELEMENTS(*e.xind)

  nf = par.n_fields
  with_em = (nf GT 1) AND (par.beta GT 1e-5)
  IF with_em THEN BEGIN
    n_fluxes = 4 * (par.n_moms / 6) ; G_es, Q_es, G_em, Q_em
    vars = INDGEN(nf+par.n_moms)
  ENDIF ELSE BEGIN
    n_fluxes = 2 * (par.n_moms / 6) ; G_es, Q_es (passing, trapped, FLR)
    IF (par.n_moms EQ 6) THEN vars = [0,nf,nf+1,nf+2] ELSE BEGIN $
      vars = 0
      FOR j = 0, par.n_moms / 6 - 1 DO $
        vars = [vars,nf+6*j,nf+1+6*j,nf+2+6*j]
    ENDELSE
  ENDELSE

  jac_fac = (*series.geom).jacobian
  FOR x = 0, par.nx0 - 1 DO $
    jac_fac[x,*] /= TOTAL(jac_fac[x,*])

  ; for comparison with nrg (volume average):
  ;*e.xind = INDGEN(par.nx0)
  ;nx = par.nx0
  ;jac_fac = (*series.geom).jacobian / $
  ;  (TOTAL((*series.geom).jacobian) * nx)

  sqrtgxx_FS = TOTAL((*series.geom).gxx*jac_fac,2)

  i = set_internal_vars(diag,{$
    with_em     : with_em,$
    xind        : *e.xind,$
    xvals       : par.x0+(-0.5*series.lx+INDGEN(par.nx0)*series.dx)*par.rhostar,$
    rsc         : par.rhostar*par.minor_r,$
    jac_fac     : jac_fac,$
    sqrtgxx_FS  : sqrtgxx_FS[(*e.xind)],$
    n_fluxes    : n_fluxes, $
    flux        : DBLARR(nx,n_fluxes,/NOZERO),$ ; defined here to improve performance
    tmpflux     : DBLARR(nx,par.nz0,/NOZERO),$ ; defined here to improve performance
    time_id     : PTR_NEW(),$
    grad_id     : PTR_NEW(),$
    flux_id     : PTR_NEW()})

  fft_format, sxky=vars

END

;######################################################################

FUNCTION fluxgrad_calc_flux, diag, var1, var2

  COMMON global_vars

  i = (*diag).internal_vars

  nx = N_ELEMENTS((*i).xind)
  nky = par.nky0

  FOR x = 0, nx - 1 DO BEGIN
    (*i).tmpflux[x,*] = $
      DOUBLE(var1[x,0]*CONJ(var2[x,0])) + 2.0D * $
      TOTAL(DOUBLE(CONJ(var1[x,1:nky-1,*])*var2[x,1:nky-1,*]),2,/DOUBLE)
  ENDFOR

  (*i).tmpflux *= (*i).jac_fac[(*i).xind,*]

  RETURN, TOTAL((*i).tmpflux,2,/DOUBLE)

END

;######################################################################

PRO fluxgrad_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nx = n_elements((*i).xind)
  nf = par.n_fields
  rho_ref = SQRT(series.mref*series.Tref) / (series.Qref * series.Bref)

  avflux = DBLARR((*i).n_fluxes,gui.out.n_spec_sel,/NOZERO)
  avgrad = DBLARR(2,gui.out.n_spec_sel,/NOZERO)

  sqrtgxx_inv = 1.0/REBIN(REFORM(sqrt((*series.geom).gxx[(*i).xind,*]),$
       [nx,1,par.nz0]),[nx,par.nky0,par.nz0])

  ; calculate vx, Bx, grad q_par
  ; ve_x = - d phi / dy
  ve_x = - DCOMPLEX(0.0,1.0) * $
    REBIN(REBIN(REFORM(*series.ky,[1,par.nky0]),[nx,par.nky0])/$
    REBIN((series.Bref*(*series.geom).C_xy[(*i).xind]),[nx,par.nky0]),$
    [nx,par.nky0,par.nz0]) * (*mom[0,0].sxky)[(*i).xind,*,*] * $
    sqrtgxx_inv

  ; B_x = dApar/dy
  IF (*i).with_em THEN B_x = DCOMPLEX(0.0,1.0) * $
    REBIN(REBIN(REFORM(*series.ky,[1,par.nky0]),[nx,par.nky0])/$
    REBIN((series.Bref*(*series.geom).C_xy[(*i).xind]),[nx,par.nky0]),$
    [nx,par.nky0,par.nz0]) * (*mom[0,1].sxky)[(*i).xind,*,*] * $
    sqrtgxx_inv

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    n_notrap = 2 + 2 * (*i).with_em ; number of fluxes without trapdiag

    FOR j = 0, par.n_moms / 6 - 1 DO BEGIN ; passing, trapped, FLR
      ; Gamma_es = <n ve_x>
      (*i).flux[*,n_notrap*j+0] = fluxgrad_calc_flux(diag,$
        (*mom[isp,nf+6*j].sxky)[(*i).xind,*,*],ve_x) * spec[sp].dens

      ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
      var1 = ((0.5 * (*mom[isp,nf+6*j+1].sxky)[(*i).xind,*,*] + $
        (*mom[isp,nf+6*j+2].sxky)[(*i).xind,*,*]) * $
        REBIN((*spec[sp].prof).dens[(*i).xind]/series.nref,$
        [nx,par.nky0,par.nz0]) + $
        1.5 * (*mom[isp,nf+6*j].sxky)[(*i).xind,*,*] * $
        REBIN((*spec[sp].prof).temp[(*i).xind]/(series.Tref*series.Qref),$
        [nx,par.nky0,par.nz0]))
      (*i).flux[*,n_notrap*j+1] = fluxgrad_calc_flux(diag,var1,ve_x) * $
	(spec[sp].dens * spec[sp].temp)

      IF (*i).with_em THEN BEGIN
        var1 = ((*mom[isp,nf+6*j+5].sxky)[(*i).xind,*,*] * $
          REBIN((*spec[sp].prof).dens[(*i).xind],[nx,par.nky0,par.nz0]))
        ; Gamma_em = <upar * B_x>
        (*i).flux[*,n_notrap*j+2] = $
          fluxgrad_calc_flux(diag,var1,B_x) * spec[sp].dens
        ; Q_em = (qpar + qperp + 5/2 upar) B_x
        ; note: mom(5) und mom(6) are NOT identical to qpar and qperp,
        ; but to qpar + 1.5 upar and qperp + upar
        (*i).flux[*,n_notrap*j+3] = fluxgrad_calc_flux(diag,$
          ((*mom[isp,nf+6*j+3].sxky)[(*i).xind,*,*] + $
          (*mom[isp,nf+6*j+4].sxky)[(*i).xind,*,*]),B_x) * $
          (spec[sp].dens * spec[sp].temp)
      ENDIF ; with_em
    ENDFOR

    output = NOT !QUIET ; switch on for comparison with nrg file
    IF output THEN BEGIN
      PRINT, '--- fluxgrad ---'
      PRINT, spec[sp].name
      PRINT, 'G_es = ', $
        STRING(TOTAL((*i).flux[*,0],1)/nx,FORMAT='(E12.4)')
      PRINT, 'Q_es = ', $
        STRING(TOTAL((*i).flux[*,1],1)/nx,FORMAT='(E12.4)')
      IF (*i).with_em THEN BEGIN
        PRINT, 'G_em = ', $
          STRING(TOTAL((*i).flux[*,2],1)/nx,FORMAT='(E12.4)')
        PRINT, 'Q_em = ', $
          STRING(TOTAL((*i).flux[*,3],1)/nx,FORMAT='(E12.4)')
      ENDIF
    ENDIF

    ; flux surface averaged temperatures and densities
    ; taking all radial grid points due to derivatives
    temp = TOTAL((REFORM((*mom[isp,par.n_fields+1].sxky)[*,0,*])/3.0+$
      REFORM((*mom[isp,par.n_fields+2].sxky)[*,0,*])*(2.0/3.0))*$
      (*i).jac_fac[*,*],2)
    dens = TOTAL(REFORM((*mom[isp,par.n_fields].sxky)[*,0,*])*$
      (*i).jac_fac[*,*],2)

    total_temp = (*spec[sp].prof).temp + temp * (*i).rsc
    total_dens = (*spec[sp].prof).dens + dens * (*i).rsc
    omt = - DERIV((*i).xvals,ALOG((*spec[sp].prof).temp+temp*(*i).rsc)) / $
      par.minor_r
    omn = - DERIV((*i).xvals,ALOG((*spec[sp].prof).dens+dens*(*i).rsc)) / $
      par.minor_r

    avgrad[0,isp] = TOTAL(omt[(*i).xind],/DOUBLE) / nx
    avgrad[1,isp] = TOTAL(omn[(*i).xind],/DOUBLE) / nx

    FOR n = 0, (*i).n_fluxes - 1 DO avflux[n,isp] = $
      TOTAL((*i).flux[*,n]/(spec[sp].dens*total_dens[(*i).xind]*$
      spec[sp].temp*total_temp[(*i).xind]*omt[(*i).xind]*(*i).sqrtgxx_FS),$
      1,/DOUBLE) / nx
  ENDFOR ; species loop

  (*i).time_id = store_step((*i).time_id,mom_time)
  (*i).flux_id = store_step((*i).flux_id,avflux)
  (*i).grad_id = store_step((*i).grad_id,avgrad)

END

;######################################################################

PRO fluxgrad_plot, xval, yval, xtitle, ytitle, timeline, $
  titleline=titleline, noerase=noerase

  COMMON global_vars

  IF NOT KEYWORD_SET(titleline) THEN titleline = ''

  ymax = MAX(yval,MIN=ymin)
  yrange = [ymin,ymax]

  xmax = MAX(xval,MIN=xmin)
  xrange = [xmin<0,1.1*xmax]

  PLOT, xval, yval, COLOR=1, XTITLE=xtitle, $
    /XSTYLE, YMARGIN=[6,2], $; POSITION=[0.15,0.4,0.9,0.9], $
    YRANGE=yrange, YTITLE=ytitle,$
    TITLE=titleline, XRANGE=xrange, $
    NOERASE=noerase, /NODATA
  OPLOT, xval, yval, color=4, PSYM=4
  OPLOT, xval, yval, color=4

  XYOUTS, 0.05, 0.01, /NORMAL, timeline, COLOR=1

;  plot_legend, colorarr[0:np-1], var_names[0:np-1], per_line=2

END

;######################################################################

PRO fluxgrad_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  time = store_step((*i).time_id,/get,/reg_array)
  avflux = REFORM(store_step((*i).flux_id,/get,/reg_array),$
    [(*i).n_fluxes,gui.out.n_spec_sel,series.step_count])
  avgrad = REFORM(store_step((*i).grad_id,/get,/reg_array),$
    [2,gui.out.n_spec_sel,series.step_count])

  n_notrap = 2 + 2 * (*i).with_em

  var_names = STRARR(n_notrap)
  var_names2 = STRARR(n_notrap)
  var_names[0:1] = ['!7C!6!Des','!6Q!Des']
  var_names2[0:1] = ['G_es','Q_es']
  IF (*i).with_em THEN BEGIN
    var_names[2:3] = ['!7C!6!Dem','!6Q!Dem']
    var_names2[2:3] = ['G_em','Q_em']
  ENDIF

  header = STRARR(par.n_moms/6*n_notrap)
  header[0:n_notrap-1] = var_names2 + ',p'
  IF par.n_moms / 6 GT 1 THEN BEGIN
    FOR j = 1, par.n_moms / 6 - 2 DO $
      header[j*n_notrap:(j+1)*n_notrap-1] = $
      var_names2 + ',t' + rm0es(par.n_moms/6-2-j)
    header[j*n_notrap:(j+1)*n_notrap-1] = var_names2 + ',flr'
  ENDIF

  rho_str = get_var_string(/rhostr)

  title = 'avg. over x/a = ' + $
    rm0es((*i).xvals[(*i).xind[0]],prec=2) + '-' + $
    rm0es((*i).xvals[(*i).xind[N_ELEMENTS((*i).xind)-1]],prec=2)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    timeline = '!6t=' + rm0es(gui.out.start_t) + '-' + $
      rm0es(gui.out.end_t) + ' ' + get_var_string(1,/time,/ounit) + $
      ', ' + rm0es(series.step_count) + ', ' + spec[sp].name

    set_output, diag, sp, /ps, multi=[0,1,1], $
      ysize=20, xsize=14.4

    fluxgrad_plot, time, avflux[1,isp,*], $
      't,' , '<!7v!6!D' + spec[sp].name + '!N>', $
      timeline, titleline=title

    fluxgrad_plot, time, avgrad[0,isp,*], $
      't', '<R/L!IT' + spec[sp].name + '!N>', $
      timeline, titleline=title

    fluxgrad_plot, avgrad[0,isp,*], avflux[1,isp,*], $
      '<R/L!DT' + spec[sp].name + '!N>', $
      '<!7v!6!D' + spec[sp].name + '!N>', $
      timeline, titleline=title

    chi_xt_avg = INT_TABULATED(time, avflux[1,isp,*])/$
                 (time[series.step_count-1]-time[0])
    omt_xt_avg = INT_TABULATED(time, avgrad[0,isp,*])/$
                 (time[series.step_count-1]-time[0])   
    
    print, 'Chi'+spec[sp].name+'; '+title+': ', rm0es(chi_xt_avg)+' Chi_gb'
    print, 'omt_'+spec[sp].name+'; '+title+': ', rm0es(omt_xt_avg)

    ; uncomment for Dimits curve:
    max_grad = CEIL(MAX(avgrad[0,isp,*]))
    grad = 6.0 + INDGEN(10) * ((max_grad > 6.5 - 6.0) / max_grad)
    chi_dimits = 33.88 * (1.0 - 6.0 / grad)
    OPLOT, grad, chi_dimits, COLOR=2, LINESTYLE=2

    set_output, diag, sp, header=['t','<R/LT>','chi'], $
      dat=[[time],[REFORM(avgrad[0,isp,*])],[REFORM(avflux[1,isp,*])]]

    set_output, diag, sp, /reset
  ENDFOR ;-- isp

END
