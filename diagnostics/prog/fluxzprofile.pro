FUNCTION fluxzprofile_info

  RETURN, {$
    type      : 'mom', $ ;'mom_uni', $
    title     : 'Parallel flux profiles', $
    help_text : ['Plots parallel (z) profiles of transport fluxes, ' + $
    	    	 'averaged over the remaining spatial coordinates and time.'],$
    ext_vars  : ['momentum','1','toggle momentum transport calculation']}

END

;######################################################################

PRO fluxzprofile_init, diag

  COMMON global_vars

  e = (*diag).external_vars
  IF NOT KEYWORD_SET(*e.momentum) THEN *e.momentum = 0

  IF par.kx_center NE 0 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + ': kx_center shifts ' + $
      'currently not possible'
    (*diag).selected = 0
    RETURN
  ENDIF

  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  with_bpar = (par.n_fields GT 2) AND (par.beta GT 1e-5) AND $
    (par.n_moms GT 6)
  nf = par.n_fields
  nm0 = with_bpar ? 9 : 6 ; n_moms without trapdiag
  IF with_em THEN BEGIN
    ; G_es, Q_es, G_em, Q_em, P_es, P_em
    n_fluxes = (4 + 2 * (*e.momentum) + 2 * with_bpar) * (par.n_moms / nm0)
    vars = INDGEN(par.n_fields+par.n_moms)
  ENDIF ELSE BEGIN
    ; G_es, Q_es, P_es (passing, trapped, FLR)
    n_fluxes = (2+(*e.momentum)) * (par.n_moms / 6)
    IF par.n_moms EQ 6 THEN BEGIN
      vars = [0,nf,nf+1,nf+2]
      IF *e.momentum THEN vars = [vars,nf+5]
    ENDIF ELSE BEGIN
      vars = 0
      FOR j = 0, par.n_moms / 6 - 1 DO vars = *e.momentum ? $
        [vars,nf+6*j,nf+1+6*j,nf+2+6*j,nf+5+6*j] : $
        [vars,nf+6*j,nf+1+6*j,nf+2+6*j]
    ENDELSE
;    vars = par.n_moms EQ 6 ? [0,nf,nf+1,nf+2] : $
;      [0,nf,nf+1,nf+2,nf+6,nf+7,nf+8,nf+12,nf+13,nf+14]
  ENDELSE

  i = set_internal_vars(diag,{$
    zind        : *par.z,$
    with_em     : with_em,$
    with_bpar   : with_bpar,$
    nm0         : nm0,$
    zprofile    : FLTARR(par.nz0,n_fluxes,/NOZERO),$
    sumzprof_id : PTR_NEW(),$
    momentum    : *e.momentum,$
    n_fluxes    : n_fluxes})

  fft_format, kxky=vars

END

;######################################################################

FUNCTION CalcZflux, diag, var1, var2

  RETURN, 2.0 * TOTAL(TOTAL(N_ELEMENTS(var2) LT 1 ? var1 : $
    CONJ(var1)*var2,1),1)

END

;######################################################################

PRO fluxzprofile_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nf = par.n_fields

  sumzprof = FLTARR(par.nz0,(*i).n_fluxes,gui.out.n_spec_sel,/NOZERO)
  IF gui.out.n_spec_sel LT 2 THEN sumzprof = REFORM(sumzprof,$
    [par.nz0,(*i).n_fluxes,gui.out.n_spec_sel],/OVERWRITE)

  sqrtgxx_inv = par.x_local?1.0/REBIN(REFORM(sqrt((*series.geom).gxx),$
                            [1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0]):$
                1.0/REBIN(REFORM(sqrt((*series.geom).gxx),$
                         [par.nx0,1,par.nz0]),[par.nx0,par.nky0,par.nz0])

  ; ve_x = - d phi / dy
  ve_x = COMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO)
  ve_x[0,0,0] = COMPLEX(0.0,1.0) * (par.x_local ? -1.0 : $
    REBIN(-1.0/(*series.geom).C_xy,[par.nkx0,par.nky0,par.nz0])) * $
    REBIN(REFORM((*series.ky),[1,par.nky0,1]),$
    [par.nkx0,par.nky0,par.nz0]) * (*mom[0,0].kxky) * sqrtgxx_inv
  ; B_x = d A_par / dy
  IF (*i).with_em THEN BEGIN
    B_x = COMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO)
    B_x[0,0,0] = COMPLEX(0.0,1.0) * (par.x_local ? -1.0 : $
      REBIN(1.0/(*series.geom).C_xy,[par.nkx0,par.nky0,par.nz0])) * $
      REBIN(REFORM((*series.ky),[1,par.nky0,1]),$
      [par.nkx0,par.nky0,par.nz0]) * (*mom[0,1].kxky) * sqrtgxx_inv

    IF (*i).with_Bpar THEN BEGIN
      ; (d B_par / dy) / B0(z)
      dBpar_dy = COMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO)
      dBpar_dy[0,0,0] = COMPLEX(0.0,1.0) * $
        REBIN(REFORM(1.0/series.Bref*(*series.ky),[1,par.nky0,1]),$
        [par.nkx0,par.nky0,par.nz0]) * (*mom[0,2].kxky) * sqrtgxx_inv

      ; include prefactor: 1 / B0(z)
      dBpar_dy *= REBIN(REFORM(1.0/(*series.geom).Bfield,$
        [1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0])
    ENDIF
  ENDIF

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    ; number of fluxes without trapdiag
    n_notrap = 2 + 2 * ((*i).with_em + (*i).with_bpar) + (1+(*i).with_em)*(*i).momentum
    FOR j = 0, par.n_moms / (*i).nm0 - 1 DO BEGIN ; passing, trapped, FLR
      jmod = nf + (*i).nm0 * j

      ; Gamma_es = <n * ve_x>
      (*i).zprofile[*,n_notrap*j+0] = CalcZflux(diag,$
        (*mom[isp,jmod].kxky),ve_x) * spec[sp].dens
      ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
      (*i).zprofile[*,n_notrap*j+1] = CalcZflux(diag,$
        (0.5*(*mom[isp,jmod+1].kxky)/series.nref+$
         (*mom[isp,jmod+2].kxky)/series.nref+$
         1.5*(*mom[isp,jmod].kxky)/series.Tref),ve_x) * $
         spec[sp].dens * spec[sp].temp
      
      ; P_es = m * upar * ve_x
      IF (*i).momentum THEN BEGIN
         (*i).zprofile[*,n_notrap*j+2+2*(*i).with_em+2*(*i).with_bpar] = CalcZflux(diag,$
         (*mom[isp,jmod+5].kxky),ve_x) * spec[sp].dens * spec[sp].mass
      ENDIF
 


      IF (*i).with_em THEN BEGIN
        IF NOT (*i).with_bpar THEN BEGIN
          ; Gamma_em = <upar * B_x>
          (*i).zprofile[*,n_notrap*j+2] = CalcZflux(diag,$
            (*mom[isp,jmod+5].kxky),B_x) * spec[sp].dens
          ; Q_em = (qpar + qperp + 5/2 upar) B_x
          ;        (upar contained in qpar/qperp)
          (*i).zprofile[*,n_notrap*j+3] = CalcZflux(diag,$
            ((*mom[isp,jmod+3].kxky)+(*mom[isp,jmod+4].kxky)),$
            B_x) * spec[sp].dens * spec[sp].temp
        ENDIF ELSE BEGIN

          ; Gamma_em = <upar * B_x - densI1 * d/dy B_par>
          (*i).zprofile[*,n_notrap*j+2] = CalcZflux(diag,$
            CONJ(*mom[isp,jmod+5].kxky)*B_x*spec[sp].dens)
          (*i).zprofile[*,n_notrap*j+4] = CalcZflux(diag,$
            -CONJ(*mom[isp,jmod+6].kxky)*dBpar_dy*$
            (spec[sp].temp*spec[sp].dens/spec[sp].charge))
          ; Q_em = (qpar + qperp + 5/2 upar) B_x - (TparI1 + TperpI1) d/dy B_par
          (*i).zprofile[*,n_notrap*j+3] = CalcZflux(diag,$
            CONJ((*mom[isp,jmod+3].kxky)+(*mom[isp,jmod+4].kxky))*B_x*$
            (spec[sp].dens*spec[sp].temp))
          (*i).zprofile[*,n_notrap*j+5] = CalcZflux(diag,$
            -CONJ((*mom[isp,jmod+7].kxky)+(*mom[isp,jmod+8].kxky))*dBpar_dy*$
            (spec[sp].temp^2*spec[sp].dens/spec[sp].charge))
        ENDELSE
        IF (*i).momentum THEN BEGIN
           (*i).zprofile[*,n_notrap*j+3+2*(*i).with_em+2*(*i).with_bpar] = CalcZflux(diag,$
           (*mom[isp,jmod+1].kxky)*series.nref+$
            (*mom[isp,jmod].kxky)*series.Tref,B_x) * spec[sp].dens * spec[sp].temp
        ENDIF
      ENDIF
    ENDFOR

    sumzprof[0,0,isp] = (*i).zprofile

    IF 0 THEN BEGIN ; for comparison with nrg file
      combine_em_bpar =1

      IF isp EQ 0 THEN PRINT, '--- fluxzprofile ---'

      G_es = 0.0
      G_em = 0.0
      G_bp = 0.0
      Q_es = 0.0
      Q_em = 0.0
      Q_bp = 0.0
      FOR j = 0, par.n_moms / (*i).nm0 - 1 DO BEGIN
        FOR z = 0, par.nz0 - 1 DO BEGIN
          G_es += (*i).zprofile[z,0+n_notrap*j] * (*series.geom).jac_norm[z]
          Q_es += (*i).zprofile[z,1+n_notrap*j] * (*series.geom).jac_norm[z]
          IF (*i).with_em THEN BEGIN
            G_em += (*i).zprofile[z,2+n_notrap*j] * (*series.geom).jac_norm[z]
            Q_em += (*i).zprofile[z,3+n_notrap*j] * (*series.geom).jac_norm[z]
            IF (*i).with_bpar THEN BEGIN
              G_bp += (*i).zprofile[z,4+n_notrap*j] * (*series.geom).jac_norm[z]
              Q_bp += (*i).zprofile[z,5+n_notrap*j] * (*series.geom).jac_norm[z]
            ENDIF
          ENDIF
        ENDFOR
      ENDFOR
      PRINT, spec[sp].name + ':'
      PRINT, 'G_es = ', G_es / par.nz0, FORMAT='(A,E11.4)'
      PRINT, 'Q_es = ', Q_es / par.nz0, FORMAT='(A,E11.4)'
      IF combine_em_bpar THEN BEGIN
        PRINT, 'G_em = ', (G_em + G_bp) / par.nz0, FORMAT='(A,E11.4)'
        PRINT, 'Q_em = ', (Q_em + Q_bp) / par.nz0, FORMAT='(A,E11.4)'
      ENDIF ELSE BEGIN
        PRINT, 'G_em = ', G_em / par.nz0, FORMAT='(A,E11.4)'
        PRINT, 'Q_em = ', Q_em / par.nz0, FORMAT='(A,E11.4)'
        PRINT, 'G_bp = ', G_bp / par.nz0, FORMAT='(A,E11.4)'
        PRINT, 'Q_bp = ', Q_bp / par.nz0, FORMAT='(A,E11.4)'
      ENDELSE
    ENDIF
  ENDFOR ; --- species loop

  (*i).sumzprof_id = time_avg((*i).sumzprof_id,sumzprof,mom_time,$
    fwd=gui.out.res_steps)

END

;######################################################################

PRO fluxzprofile_plot, diag, xind, Fprofile, var_names

  COMMON global_vars
  i = (*diag).internal_vars

  nprofiles = N_ELEMENTS(Fprofile[0,*])
  IF NOT KEYWORD_SET(titleline) THEN titleline = ''
  maxy = MAX(Fprofile,MIN=miny)
  IF maxy LT 0 THEN maxy = 0 ; always show the base line
  IF miny GT 0 THEN miny = 0
  ;colors adjusted such that G_es, G_em,
  ;Q_es, Q_em, P_es, P_em are consistent between fluxzprof and fluxspectra
  if (*i).with_bpar THEN BEGIN
     colorarr = [1,3,2,4,5,7,6,8] 
  ENDIF ELSE BEGIN
     colorarr = [1,3,2,4,6,8]
  ENDELSE

  xrange = [-!PI*par.n_pol,!PI*par.n_pol]
  xtickv = 0.5 * (INDGEN(4*par.n_pol+1) - 2 * par.n_pol)
  xtickname = STRARR(4*par.n_pol+1)
  xtickname = rm0es(xtickv) + '!7p!6'
  xtickname[2*par.n_pol] = '0'
  xtickv *= !PI
  IF par.n_pol GT 1 THEN xtickname[2*INDGEN(par.n_pol+2)+1] = ' '

  PLOT, xind, Fprofile[*,0], COLOR=colorarr[0], XTITLE='!6z/qR', $
    /XSTYLE, XRANGE=xrange, XTICKS=N_ELEMENTS(xtickv)-1, XMINOR=2, $
    XTICKV=xtickv, XTICKNAME=xtickname, XTICKLEN=1.0, TITLE='!6Fluxes ', $
    POSITION=[0.15,0.4,0.9,0.9], YRANGE=[miny,maxy], /NODATA
  PLOTS, !X.CRANGE, [0,0], COLOR=1

  FOR ivar = 0, nprofiles - 1 DO OPLOT, xind, Fprofile[*,ivar], $
    COLOR=colorarr[ivar]

  plot_info_str, diag

  plot_legend, colorarr[0:nprofiles-1], var_names[0:nprofiles-1], $
    per_line=2

END

;######################################################################

PRO fluxzprofile_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  sumzprof = time_avg((*i).sumzprof_id,/avg,fwd=gui.out.res_steps,tarr=time)

  n_notrap = 2 + 2 * ((*i).with_em + (*i).with_bpar) + (1+(*i).with_em)*(*i).momentum
  var_names = STRARR(n_notrap)
  var_names[0:1] = ['!7C!6!Des','!6Q!Des']
  IF (*i).with_em THEN var_names[2:3] = ['!7C!6!Dem','!6Q!Dem']
  IF (*i).with_bpar THEN var_names[4:5] = ['!7C!6!S!Dem!N!R!UBpar!N!D',$
    '!6Q!S!Dem!N!R!UBpar!N!D']
  IF (*i).momentum THEN  var_names[2*(1+(*i).with_em+(*i).with_bpar)] = '!7P!6!Des'
  IF (*i).momentum and (*i).with_em THEN  var_names[2*(1+(*i).with_em+(*i).with_bpar)+1] = '!7P!6!Dem'

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, multi=[0,1,1], ysize=20, xsize=14.4

    FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
      IF (*i).n_fluxes EQ n_notrap THEN BEGIN ; no trapdiag
        fluxzprofile_plot, diag, (*i).zind, sumzprof[*,*,isp,n], $
          var_names + '!N'
      ENDIF ELSE BEGIN
        fluxzprofile_plot, diag, (*i).zind, sumzprof[*,0:n_notrap-1,isp,n], $
          var_names + ',pass!N'
        fluxzprofile_plot, diag, (*i).zind, sumzprof[*,n_notrap:2*n_notrap-1,isp,n], $
          var_names + ',trp!N'
        fluxzprofile_plot, diag, (*i).zind, sumzprof[*,0:n_notrap-1,isp,n] + $
          sumzprof[*,n_notrap:2*n_notrap-1,isp,n] + $
          sumzprof[*,2*n_notrap:3*n_notrap-1,isp,n], var_names + ',total!N'
      ENDELSE
    ENDFOR

    set_output, diag, sp, header=['z',var_names], $
      dat=[[(*i).zind],[sumzprof[*,*,isp,0]]]

    set_output, diag, sp, /reset
  ENDFOR ;--- species loop

END
