FUNCTION fluxspectra_y_info

  RETURN, {$
    type      : 'mom global',$ ;'mom_uni',$
    title     : 'Flux spectra (ky/kx) (y-global)',$
    help_text : ['Plots flux spectra in ky/kx space, averaged over '+$
    	    	 'the remaining coordinates and time. Dotted lines '+$
                 'in the log/log plot indicate negative values.'],$
    ext_vars  : [['pages','0','select plot pages; default: all '+$
                  '([1,1,1,1,1,1]); -1: [0,0,1,0,0,1]'],$
                 ['renorm','1','switch on/off normalization due to '+$
                  'total flux'],$
                 ['k_mult','1','show (rescaled) multiplication with '+$
                  'k as dashed line in log plots (gives better '+$
                  'impression of transport portions)'],$
                 ['momentum','1','add parallel momentum']]}

END

;######################################################################

PRO fluxspectra_y_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.pages) LT 1 THEN *e.pages = [1,1,1,1,1,1]
  IF (*e.pages)[0] EQ -1 THEN *e.pages = [0,0,1,0,0,1]
  IF N_ELEMENTS(*e.pages) NE 6 THEN BEGIN
    n_pages_old = N_ELEMENTS(*e.pages)
    pages_new = INTARR(6)
    pages_new[0:(n_pages_old-1)<5] = (*e.pages)[0:(n_pages_old-1)<5]
    *e.pages = pages_new
  ENDIF
  IF NOT KEYWORD_SET((*e.renorm)) THEN (*e.renorm) = 0
  IF NOT KEYWORD_SET(*e.k_mult) THEN *e.k_mult = 0
  IF NOT KEYWORD_SET(*e.momentum) THEN *e.momentum = 0
  IF par.lx EQ 0 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + $
      ': ky dependent lx currently not allowed'
    (*diag).selected = 0
    RETURN
  ENDIF

  nf = par.n_fields
  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  IF with_em THEN BEGIN
    ; fluxes: G_es, Q_es, G_em, Q_em[, P_es, P_em]
    n_fluxes = (4 + 2 * (*e.momentum)) * (par.n_moms / 6)
    vars = INDGEN(par.n_fields+par.n_moms)
  ENDIF ELSE BEGIN
    ; G_es, Q_es[,P_es] (passing, trapped, FLR)
    n_fluxes = (2 + (*e.momentum)) * (par.n_moms / 6)
    IF par.n_moms EQ 6 THEN BEGIN
      vars = [0,nf,nf+1,nf+2]
      IF *e.momentum THEN vars = [vars,nf+5]
    ENDIF ELSE BEGIN
      vars = 0
      FOR j = 0, par.n_moms / 6 - 1 DO vars = *e.momentum ? $
        [vars,nf+6*j,nf+1+6*j,nf+2+6*j,nf+5+6*j] : $
        [vars,nf+6*j,nf+1+6*j,nf+2+6*j]
    ENDELSE
  ENDELSE

  i = set_internal_vars(diag,{$
    kyind          : (*series.ky),$
    kxind          : 2.0*!PI*(INDGEN(par.nkx0/2+1))/series.lx,$
; defined here to improve performance:
    spectrumx      : FLTARR(par.nkx0/2+1,n_fluxes,/NOZERO),$ 
    tempspecx      : FLTARR(par.nkx0/2+1,par.nz0,/NOZERO),$
    spectrumy      : FLTARR(par.nky0,n_fluxes,/NOZERO),$
    tempspecy      : FLTARR(par.nky0,par.nz0,/NOZERO),$
; ---
    sumspectx_id   : PTR_NEW(),$
    sumspecty_id   : PTR_NEW(),$
    pages          : *e.pages,$
    renorm         : *e.renorm,$
    k_mult         : *e.k_mult,$
    momentum       : *e.momentum,$
    n_fluxes       : n_fluxes})

  fft_format, kxsy=vars

END

;######################################################################

FUNCTION fluxspectra_y_CalcXspectrum, diag, var

  COMMON global_vars

  i = (*diag).internal_vars

  nkx0o2 = par.nkx0 / 2 
  nky0 = par.nky0

  FOR x = 0, nkx0o2 - 1 DO BEGIN
  var[x,*,*] *= (*series.geom).jac_norm[*,*]
  ENDFOR

  (*i).tempspecx[0,*] = TOTAL(var[0,*,*],2)

  FOR x = 1, nkx0o2 - 1 DO (*i).tempspecx[x,*] = $

     2.0 * TOTAL(var[x,*,*],2)
  (*i).tempspecx[nkx0o2,*] = TOTAL(var[nkx0o2,*,*],2)

  RETURN, TOTAL((*i).tempspecx,2) / (par.nz0)

END

;######################################################################

FUNCTION fluxspectra_y_CalcYspectrum, diag, var

  COMMON global_vars

  i = (*diag).internal_vars

  nkx0o2 = par.nkx0 / 2 + 1
  nky0 = par.nky0
  ny0 = par.ny0

  yfftdata = FLTARR(nkx0o2,2*nky0,par.nz0,/NOZERO)  

  yfftdata[0:nkx0o2-1,*,0:par.nz0-1] = FFT(var[0:nkx0o2-1,*,0:par.nz0-1],DIMENSION=2)

  IF par.x_local THEN BEGIN
    (*i).tempspecy[0,*] = yfftdata[0,0,*] +  $
      2.0 * TOTAL(yfftdata[1:nkx0o2-1,0,*],1)

    (*i).tempspecy[1:nky0-1,*] = 2.0 * TOTAL(yfftdata[*,1:nky0-1,*],1)

    FOR iy = 1,nky0-1 DO BEGIN
;       (*i).tempspecy[iy,*] += yfftdata[*,ny0-iy,*]
    ENDFOR

  ENDIF ELSE BEGIN
    (*i).tempspecy[0,*] = yfftdata[0,0,*] + yfftdata[nkx0o2,0,*] + $
      2.0 * TOTAL(yfftdata[1:nkx0o2-1,0,*],1) 

    (*i).tempspecy[1:nky0-1,*] = 2.0 * $
      TOTAL(yfftdata[*,1:nky0-1,*],1) 
  ENDELSE

  RETURN, TOTAL((*i).tempspecy,2) / (par.nz0)

END

;######################################################################

PRO fluxspectra_y_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  yaxis = -series.ly * .5 + indgen(par.ny0)*series.dy
  nx0 = par.nkx0 / 2 + 1

  sumspectx = REFORM(FLTARR(par.nkx0/2+1,(*i).n_fluxes,$
    gui.out.n_spec_sel,/NOZERO),$
    [par.nkx0/2+1,(*i).n_fluxes,gui.out.n_spec_sel])
  sumspecty = REFORM(FLTARR(par.nky0,(*i).n_fluxes,$
    gui.out.n_spec_sel,/NOZERO),$
    [par.nky0,(*i).n_fluxes,gui.out.n_spec_sel])

  nf = par.n_fields

  c_ref = SQRT(series.Tref*series.Qref/series.mref)
  p_ref = series.nref * series.Tref * series.Qref

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref

  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  with_Bpar_moms = with_em AND (par.n_fields GT 2) AND (par.n_moms GT 6)

  sqrtgxx_inv = 1.0/REBIN(REFORM(sqrt((*series.geom).gxx),$
                [1,par.ny0,par.nz0]),[nx0,par.nky0,par.nz0])

  ; calculate vEx, Bx
  ; ve_x = - d phi / dy
  ve_x = FLTARR(nx0,par.ny0,par.nz0,/NOZERO)
;  ve_x[0,0,0] = COMPLEX(0.0,1.0) * $
;    REBIN(REFORM(-1.0/series.Bref*(*series.ky),[1,par.nky0,1]),$
;    [par.nkx0,par.nky0,par.nz0]) * (*mom[0,0].kxsy)
 
  FOR ix = 0 , nx0 - 1 DO BEGIN
    FOR iz = 0 , par.nz0 - 1 DO BEGIN

       ve_x[ix,0:par.ny0-1,iz] = DERIV(yaxis,(*mom[0,0].kxsy)[ix,0:par.ny0-1,iz])

    ENDFOR
  ENDFOR

;  ve_x *= sqrtgxx_inv           
  ve_x /= (*series.geom).C_xy

  IF with_em THEN BEGIN
  
    B_x = FLTARR(nx0,par.ny0,par.nz0,/NOZERO)
    B_x[0:nx0-1,*,0:par.nz0-1] =DERIV(yaxis,(*mom[0,1].kxsy)[0:nx0-1,*,0:par.nz0-1])
;    B_x *= sqrtgxx_inv
    B_x /= (*series.geom).C_xy

    IF with_Bpar_moms THEN BEGIN
    
      dBpar_dy = FLTARR(nx0,par.ny0,par.nz0,/NOZERO)
  
      dBpar_dy[0:nx0-1,*,0:par.nz0-1] =DERIV(yaxis,(*mom[0,2].kxsy)[0:nx0-1,*,0:par.nz0-1])
  

      ; include prefactor: 1 / B0(z)
      dBpar_dy *= REBIN(REFORM(1.0/(*series.geom).Bfield,$
        [1,1,par.nz0]),[nx0,par.nky0,par.nz0])*sqrtgxx_inv
    ENDIF
  ENDIF

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    ; number of fluxes without trapdiag
    n_notrap = (2 + (*i).momentum) * (1 + with_em)

    FOR j = 0, par.n_moms / 6 - 1 DO BEGIN ; passing, trapped, FLR
      jmod = nf + 6 * j
      ; Gamma_es = <n * ve_x>
      (*i).spectrumy[*,n_notrap*j] = fluxspectra_y_CalcYspectrum(diag,$
        CONJ(*mom[isp,jmod].kxsy)*ve_x) * spec[sp].dens
      (*i).spectrumx[*,n_notrap*j] = fluxspectra_y_CalcXspectrum(diag,$
        CONJ(*mom[isp,jmod].kxsy)*ve_x) * spec[sp].dens
      ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
      (*i).spectrumy[*,n_notrap*j+1] = fluxspectra_y_CalcYspectrum(diag,$
        CONJ((0.5*(*mom[isp,jmod+1].kxsy)[*,*,*]*series.nref+$
        (*mom[isp,jmod+2].kxsy)*series.nref+$
        1.5*(*mom[isp,jmod].kxsy)*series.Tref))*ve_x) * $
        (spec[sp].dens * spec[sp].temp * series.Qref)
      (*i).spectrumx[*,n_notrap*j+1] = fluxspectra_y_CalcXspectrum(diag,$
        CONJ((0.5*(*mom[isp,jmod+1].kxsy*series.nref)+$
        (*mom[isp,jmod+2].kxsy*series.nref)+$
        1.5*(*mom[isp,jmod].kxsy)*series.Tref))*ve_x) * $
        (spec[sp].dens * spec[sp].temp * series.Qref)

      IF (*i).momentum THEN BEGIN ; P_es = <upar ve_x>
        (*i).spectrumy[*,n_notrap*j+2*(1+with_em)] = fluxspectra_y_CalcYspectrum(diag,$
          CONJ(*mom[isp,jmod+5].kxsy)*ve_x) * $
          (spec[sp].dens * spec[sp].mass)
        (*i).spectrumx[*,n_notrap*j+2*(1+with_em)] = fluxspectra_y_CalcXspectrum(diag,$
          CONJ(*mom[isp,jmod+5].kxsy)*ve_x) * $
          (spec[sp].dens * spec[sp].mass)
      ENDIF

      IF with_em THEN BEGIN
        IF NOT with_Bpar_moms THEN BEGIN
          ; Gamma_em = <upar * B_x>
          (*i).spectrumy[*,n_notrap*j+2] = fluxspectra_y_CalcYspectrum(diag,$
            CONJ(*mom[isp,jmod+5].kxsy)*B_x) * spec[sp].dens
          (*i).spectrumx[*,n_notrap*j+2] = fluxspectra_y_CalcXspectrum(diag,$
            CONJ(*mom[isp,jmod+5].kxsy)*B_x) * spec[sp].dens
          ; Q_em = (qpar + qperp + 5/2 upar) B_x
          ; mom(5) und mom(6) are NOT identical to qpar and qperp
          ; but to qpar+1.5upar and qperp+upar
          (*i).spectrumy[*,n_notrap*j+3] = fluxspectra_y_CalcYspectrum(diag,$
            CONJ((*mom[isp,jmod+3].kxsy)+(*mom[isp,jmod+4].kxsy))*B_x) * $
            (spec[sp].dens * spec[sp].temp)
          (*i).spectrumx[*,n_notrap*j+3] = fluxspectra_y_CalcXspectrum(diag,$
            CONJ((*mom[isp,jmod+3].kxsy)+(*mom[isp,jmod+4].kxsy))*B_x) * $
            (spec[sp].dens * spec[sp].temp)
        ENDIF ELSE BEGIN
          ; Gamma_em = <upar * B_x + densI1 - d/dy B_par>
          (*i).spectrumy[*,n_notrap*j+2] = fluxspectra_y_CalcYspectrum(diag,$
            CONJ(*mom[isp,jmod+5].kxsy)*B_x*spec[sp].dens-$
            CONJ(*mom[isp,jmod+6].kxsy)*dBpar_dy*$
            (spec[sp].temp*spec[sp].dens/spec[sp].charge))
          (*i).spectrumx[*,n_notrap*j+2] = fluxspectra_y_CalcXspectrum(diag,$
            CONJ(*mom[isp,jmod+5].kxsy)*B_x*spec[sp].dens-$
            CONJ(*mom[isp,jmod+6].kxsy)*dBpar_dy*$
            (spec[sp].temp*spec[sp].dens/spec[sp].charge))
          ; Q_em = (qpar + qperp + 5/2 upar) B_x - (TparI1 + TperpI1) d/dy B_par
          ; mom(5) und mom(6) are NOT identical to qpar and qperp
          ; but to qpar+1.5upar and qperp+upar
          (*i).spectrumy[*,n_notrap*j+3] = fluxspectra_y_CalcYspectrum(diag,$
            CONJ((*mom[isp,jmod+3].kxsy)+(*mom[isp,jmod+4].kxsy))*B_x*$
            (spec[sp].dens*spec[sp].temp)-$
            CONJ((*mom[isp,jmod+7].kxsy)+(*mom[isp,jmod+8].kxsy))*dBpar_dy*$
            (spec[sp].temp^2*spec[sp].dens/spec[sp].charge))
          (*i).spectrumx[*,n_notrap*j+3] = fluxspectra_y_CalcXspectrum(diag,$
            CONJ((*mom[isp,jmod+3].kxsy)+(*mom[isp,jmod+4].kxsy))*B_x*$
            (spec[sp].dens*spec[sp].temp)-$
            CONJ((*mom[isp,jmod+7].kxsy)+(*mom[isp,jmod+8].kxsy))*dBpar_dy*$
            (spec[sp].temp^2*spec[sp].dens/spec[sp].charge))
        ENDELSE

        IF (*i).momentum THEN BEGIN ; P_em = (Tpar + n) B_x
          (*i).spectrumy[*,n_notrap*j+5] = fluxspectra_y_CalcYspectrum(diag,$
            CONJ((*mom[isp,jmod+1].kxsy)*series.nref+$
            (*mom[isp,jmod].kxsy)*series.Tref)*B_x) * $
            (spec[sp].dens * spec[sp].mass)
          (*i).spectrumx[*,n_notrap*j+5] = fluxspectra_y_CalcXspectrum(diag,$
            CONJ((*mom[isp,jmod+1].kxsy)*series.nref+$
            (*mom[isp,jmod].kxsy)*series.Tref)*B_x) * $
            (spec[sp].dens * spec[sp].mass)
        ENDIF
      ENDIF
    ENDFOR

    sumspectx[0,0,isp] = (*i).spectrumx
    sumspecty[0,0,isp] = (*i).spectrumy

    IF 0 THEN BEGIN ; for comparison with nrg file
      PRINT, '--- fluxspectra ---'
      IF (*i).momentum THEN PRINT, $
        'warning: momentum switch needs to be off!'
      G_es = 0.0
      G_em = 0.0
      Q_es = 0.0
      Q_em = 0.0
      FOR j = 0, par.n_moms / 6 - 1 DO BEGIN
        G_es += TOTAL((*i).spectrumy[*,0+n_notrap*j],1)
        Q_es += TOTAL((*i).spectrumy[*,1+n_notrap*j],1)
        IF with_em THEN BEGIN
          G_em += TOTAL((*i).spectrumy[*,2+n_notrap*j],1)
          Q_em += TOTAL((*i).spectrumy[*,3+n_notrap*j],1)
        ENDIF
      ENDFOR
      PRINT, spec[sp].name
      PRINT, 'G_es = ', rm0es(G_es,prec=5)
      PRINT, 'Q_es = ', rm0es(Q_es,prec=5)
      PRINT, 'G_em = ', rm0es(G_em,prec=5)
      ;PRINT, 'from x spectrum: Q_es = ', rm0es(TOTAL(sumspectx[*,2,0],1))
      PRINT, 'Q_em = ', rm0es(Q_em,prec=5)
    ENDIF
  ENDFOR ; --- species loop

  (*i).sumspectx_id = time_avg((*i).sumspectx_id,sumspectx,mom_time,$
    fwd=gui.out.res_steps)
  (*i).sumspecty_id = time_avg((*i).sumspecty_id,sumspecty,mom_time,$
    fwd=gui.out.res_steps)

END

;######################################################################

PRO fluxspectra_y_plot, diag, xind, spectra, xtitle, var_names, $
  titleline=titleline, k_mult=k_mult, pages=pages, time=time
  
  nspectra = N_ELEMENTS(spectra[0,*])
  maxy = MAX(spectra,MIN=miny)
  colorarr = [1,3,2,4,6,8] ; to be color consistent with older plots
  IF NOT KEYWORD_SET(titleline) THEN titleline = ''
  IF NOT KEYWORD_SET(k_mult) THEN k_mult = 0
  IF N_ELEMENTS(pages) NE 3 THEN $
    STOP, 'fluxspectra_y_plot error: no pages specified'

  FOR j = 1, 3 DO BEGIN ; log-log/log-lin/lin-lin
    IF pages[j-1] THEN BEGIN
      xrange = [xind[0],xind[N_ELEMENTS(xind)-1]]

      IF j EQ 1 THEN BEGIN
        ytickformat = 'betterticks'
        yrange = [1e-5,maxy]
      ENDIF ELSE BEGIN
        ytickformat = ''
        yrange=[miny,maxy]
      ENDELSE

      IF j NE 3 THEN BEGIN
        xtickformat = 'betterticks'
        IF xind[0] EQ 0 THEN xrange[0] = xind[1]
      ENDIF ELSE xtickformat = ''
    
      PLOT, xind, spectra[*,0], COLOR=colorarr[0], XTITLE=xtitle, $
        XLOG=(j/3+1), XTICKFORMAT=xtickformat, /XSTYLE, $
        XTICKLEN=1.0, POSITION=[0.15,0.4,0.9,0.9], $
        YLOG=(1/j), YTICKFORMAT=ytickformat, YRANGE=yrange, $
        TITLE=titleline, XRANGE=xrange, /NODATA

      FOR var = 0, nspectra - 1 DO BEGIN
        IF j EQ 1 THEN $
          OPLOT, xind, ABS(spectra[*,var]), COLOR=colorarr[var], LINE=1
        OPLOT, xind, spectra[*,var], COLOR=colorarr[var]
      ENDFOR

      IF k_mult AND (2 / j) THEN BEGIN
        FOR var = 0, nspectra - 1 DO BEGIN
          scalfac = MAX(ABS(spectra[*,var])) / MAX(ABS(xind*spectra[*,var]))
          OPLOT, xind, scalfac * xind * spectra[*,var], $
            COLOR=colorarr[var], LINESTYLE=2
        ENDFOR
      ENDIF

      IF (j EQ 2) THEN PLOTS, 10^!X.CRANGE, [0,0], COLOR=1
      IF (j EQ 3) THEN PLOTS, !X.CRANGE, [0,0], COLOR=1

      ; print run data (one line)
      plot_info_str, diag, time=time

      plot_legend, colorarr[0:nspectra-1], var_names[0:nspectra-1], $
        per_line=2
    ENDIF
  ENDFOR

END

;######################################################################

PRO fluxspectra_y_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  sumspectx = time_avg((*i).sumspectx_id,/avg,fwd=gui.out.res_steps,$
    tarr=tarr)
  sumspecty = time_avg((*i).sumspecty_id,/avg,fwd=gui.out.res_steps,$
    tarr=tarr)

  IF (*i).renorm THEN BEGIN
    total_flux = TOTAL(sumspecty,1)
    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      FOR ivar = 0, (*i).n_fluxes - 1 DO BEGIN
        sumspecty[*,ivar,isp] = sumspecty[*,ivar,isp] / $
          total_flux[ivar,isp]
        sumspectx[*,ivar,isp] = sumspectx[*,ivar,isp] / $
          total_flux[ivar,isp]
      ENDFOR
    ENDFOR
  ENDIF

  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  n_notrap = (2 + (*i).momentum) * (1 + with_em)
  var_names = STRARR(n_notrap)
  var_names2 = STRARR(n_notrap)
  var_names[0:1] = ['!7C!6!Des','!6Q!Des']
  var_names2[0:1] = ['G_es','Q_es']
  IF (*i).momentum THEN BEGIN
    var_names[2*(1+with_em)] = '!7P!6!Des'
    var_names2[2*(1+with_em)] = 'P_es'
  ENDIF
  IF with_em THEN BEGIN
    var_names[2:3]=['!7C!6!Dem','!6Q!Dem']
    var_names2[2:3]=['G_em','Q_em']
    IF (*i).momentum THEN BEGIN
      var_names[5] = '!7P!6!Dem'
      var_names2[5] = 'P_em'
    ENDIF
  ENDIF

  header = STRARR(par.n_moms/6*n_notrap)
  header[0:n_notrap-1] = var_names2 + ',p'
  IF par.n_moms / 6 GT 1 THEN BEGIN
    FOR j = 1, par.n_moms / 6 - 2 DO $
      header[j*n_notrap:(j+1)*n_notrap-1] = var_names2 + ',t' + $
      rm0es(par.n_moms/6-2-j)
    header[j*n_notrap:(j+1)*n_notrap-1] = var_names2 + ',flr'
  ENDIF

  rho_str = get_var_string(/rhostr)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, multi=[0,1,1], ysize=20, xsize=14.4

    FOR n = 0, (series.step_count - 1) * gui.out.res_steps DO BEGIN
      IF gui.out.res_steps THEN time = tarr[n]

      IF (*i).n_fluxes EQ n_notrap THEN BEGIN
        fluxspectra_y_plot, diag, (*i).kyind,sumspecty[*,*,isp,n], '!6k!Dy!N' + $
          rho_str, var_names + '!N', k_mult=(*i).k_mult, $
          pages=(*i).pages[0:2], time=time
        fluxspectra_y_plot, diag, (*i).kxind,sumspectx[*,*,isp,n], '!6k!Dx!N' + $
          rho_str, var_names + '!N', k_mult=(*i).k_mult, $
          pages=(*i).pages[3:5], time=time
      ENDIF ELSE BEGIN
        fluxtotalx = sumspectx[*,0:n_notrap-1,isp,n]
        FOR j = 1, par.n_moms / 6 - 1 DO fluxtotalx += $
          sumspectx[*,j*n_notrap:(j+1)*n_notrap-1,isp,n]

        vnames = STRARR(1+par.n_moms/6) ; total, pass, t1, ..., tn, flr
        data_y = FLTARR(par.nky0,1+par.n_moms/6,/NOZERO)
        data_x = FLTARR(par.nkx0/2+1,1+par.n_moms/6,/NOZERO)

        PRINT, 'fluxspectra: transport fractions (' + spec[sp].name + ')'
        FOR v = 0, n_notrap - 1 DO BEGIN
          vnames[0] = var_names[v]
          vnames[1] = var_names[v] + ',pass!N'
          data_y[*,0] = 0.0
          data_x[*,0] = 0.0
          FOR j = 0, par.n_moms / 6 - 1 DO BEGIN
            data_y[*,0] += sumspecty[*,j*n_notrap+v,isp,n]
            data_x[*,0] += sumspectx[*,j*n_notrap+v,isp,n]
          ENDFOR
          data_y[*,1] = sumspecty[*,v,isp,n]
          IF TOTAL(data_y[*,0],1) GT 1e-16 THEN PRINT, var_names2[v] + $
            ',pass' + '/' + var_names2[v] + ',total = ' + $       
            rm0es(TOTAL(data_y[*,1],1)/TOTAL(data_y[*,0],1))
          data_x[*,1] = sumspectx[*,v,isp,n]
          FOR j = 1, par.n_moms / 6 - 2 DO BEGIN
            vnames[j+1] = var_names[v] + ',trp' + $
              rm0es(par.n_moms/6-2-j) + '!N'
            data_y[*,j+1] = sumspecty[*,j*n_notrap+v,isp,n]
          IF TOTAL(data_y[*,0],1) GT 1e-16 THEN PRINT, var_names2[v] + $
            ',trp' + rm0es(par.n_moms/6-2-j) + '/' + var_names2[v] + $
            ',total = ' + rm0es(TOTAL(data_y[*,j+1],1)/TOTAL(data_y[*,0],1))
            data_x[*,j+1] = sumspectx[*,j*n_notrap+v,isp,n]
          ENDFOR
          vnames[par.n_moms/6] = var_names[v] + ',flr!N'
          data_y[*,par.n_moms/6] = sumspecty[*,(par.n_moms/6-1)*n_notrap+v,isp,n]
          IF TOTAL(data_y[*,0],1) GT 1e-16 THEN PRINT, var_names2[v] + $
            ',flr' + '/' + var_names2[v] + ',total = ' + $
            rm0es(TOTAL(data_y[*,par.n_moms/6],1)/TOTAL(data_y[*,0],1)) 
          data_x[*,par.n_moms/6] = sumspectx[*,(par.n_moms/6-1)*n_notrap+v,isp,n]
          fluxspectra_y_plot, diag, (*i).kyind, data_y, '!6k!Iy!N' + rho_str, $
            vnames, k_mult=(*i).k_mult, pages=(*i).pages[0:2], time=time
          fluxspectra_y_plot, diag, (*i).kxind, data_x, '!6k!Ix!N' + rho_str, $
            vnames, k_mult=(*i).k_mult, pages=(*i).pages[3:5], time=time
        ENDFOR
      ENDELSE

      set_output, diag, sp, header=['kx',header],$
        dat=[[(*i).kxind],[sumspectx[*,*,isp,n]]], append=(n GT 0)
      set_output, diag, sp, header=['ky',header],$
        dat=[[(*i).kyind],[sumspecty[*,*,isp,n]]], /append
    ENDFOR

    with_chi = 0
    IF gui.out.res_steps OR (*i).momentum THEN with_chi = 0
    IF with_chi THEN BEGIN
      gxx_fac = TOTAL((*series.geom).gxx*(*series.geom).jacobian) / $
        TOTAL((*series.geom).jacobian)
      inv_fac = 1.0 / (series.nref * series.Tref * series.Qref * $
        spec[sp].temp * spec[sp].omt * gxx_fac)
      PRINT, 'Q_es=', rm0es(TOTAL(sumspecty[*,1,isp])), ', ', $
        rm0es(TOTAL(sumspectx[*,1,isp]))
      PRINT, 'G_es=', rm0es(TOTAL(sumspecty[*,0,isp])), ', ', $
        rm0es(TOTAL(sumspectx[*,0,isp]))
      PRINT, 'Chi = ', rm0es(TOTAL(sumspecty[*,1,isp])*inv_fac)
      set_output, diag, sp, header=['ky','chi'], $
        dat=[[(*i).kyind],[REFORM(sumspecty[*,1,isp])*inv_fac]], /append
    ENDIF

    set_output, diag, sp, /reset
  ENDFOR ; --- species loop

END
