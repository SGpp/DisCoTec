FUNCTION checkpoint_info

  RETURN, {$
    type      : 'chpt',$
    title     : 'Plot |g|^2',$
    help_text : ['Prints scalar value or draws 1D or 2D plots of the '+$
                 'spatial mean square of the distribution '+$
                 'function g1, depending on the choice of indices.',$
                 'Data is read from "checkpoint_<run>" or '+$
                 '"checkpoint".','Note that velocity space averages '+$
                 'are not performed as mean squares (!)'],$
    ext_vars  : [['xind','0','kx/x index (-2 for all, -3 for all '+$
                  'with FFT); default: -1 (average)'],$
                 ['yind','0','ky index (-2 for all, -3 for all '+$
                  'with FFT ky2y); default: -1 (average)'],$
                 ['zind','0','z index (-1 for average, -2 for all); '+$
                  'default: nz0/2'],$
                 ['vind','0','v_par index (-2 for all); '+$
                  'default: -1 (average)'],$
                 ['wind','0','mu index (-2 for all); '+$
                  'default: -1 (average)'],$
                 ['g2f','1','plot f1 instead of g1 (local runs only)'],$
                 ['log','1','switch on/off logarithmic scaling']]}

; need to implement time avg., velocity space moments

END

;######################################################################

PRO checkpoint_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  plot_apar = 0 ; plot A_par (only for local runs)

  IF par.adapt_lx THEN printerror, 'Warning: kx labels may be ' + $
    'incorrect if adapt_lx is set'

  xfft = 0
  yfft = 0
  IF N_ELEMENTS(*e.xind) LT 1 THEN xind = -1 ELSE $
    IF (*e.xind)[0] LE -2 THEN BEGIN
  
    xind = INDGEN(par.nx0)
    IF (*e.xind)[0] EQ -3 THEN xfft = 1
    IF (par.x_local AND (xfft EQ 0)) THEN $
      xind = SHIFT(xind,par.nkx0/2-1)
  ENDIF ELSE xind = *e.xind
  IF N_ELEMENTS(*e.yind) LT 1 THEN yind = -1 ELSE $
    IF (*e.yind)[0] LE -2 THEN BEGIN

    yind = INDGEN(par.nky0)
    IF (*e.yind)[0] EQ -3 THEN yfft = 1
  ENDIF ELSE yind = *e.yind
  IF (yind[0] EQ -1) AND (par.ky0_ind GT 0) THEN yind = 0
  IF N_ELEMENTS(*e.zind) LT 1 THEN zind = par.nz0 / 2 ELSE $
    IF (*e.zind)[0] EQ -2 THEN zind = INDGEN(par.nz0) $
    ELSE zind = *e.zind
  IF N_ELEMENTS(*e.vind) LT 1 THEN vind = -1 ELSE $
    IF (*e.vind)[0] EQ -2 THEN vind = INDGEN(par.nv0) $
    ELSE vind = (*e.vind)
  IF N_ELEMENTS(*e.wind) LT 1 THEN wind = -1 ELSE $
    IF (*e.wind)[0] EQ -2 THEN wind = INDGEN(par.nw0) $
    ELSE wind = (*e.wind)

  getMuWeightsAndKnots, mu_weights, mu_knots

  i = set_internal_vars(diag,{$
    xavg       : xind[0] EQ -1,$
    yavg       : yind[0] EQ -1,$
    zavg       : zind[0] EQ -1,$
    vavg       : vind[0] EQ -1,$
    wavg       : wind[0] EQ -1,$
    xind       : xind[0] EQ -1 ? INDGEN(par.nx0) : xind,$
    yind       : yind[0] EQ -1 ? INDGEN(par.nky0) : yind,$
    zind       : zind[0] EQ -1 ? INDGEN(par.nz0) : zind,$
    vind       : vind[0] EQ -1 ? INDGEN(par.nv0) : vind,$
    wind       : wind[0] EQ -1 ? INDGEN(par.nw0) : wind,$
    nxind      : N_ELEMENTS(xind),$
    nyind      : N_ELEMENTS(yind) * (1 + yfft),$
    nzind      : N_ELEMENTS(zind),$
    nvind      : N_ELEMENTS(vind),$
    nwind      : N_ELEMENTS(wind),$
    xfft       : xfft,$
    yfft       : yfft,$
    log        : KEYWORD_SET(*e.log),$
    g2f        : KEYWORD_SET(*e.g2f),$
    plot_apar  : plot_apar,$
    mu_weights : mu_weights,$
    mu_knots   : mu_knots,$
    skip_time  : 0,$
    apar       : plot_apar ? $
                 DCOMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO) : 0,$
    time_id    : PTR_NEW(),$
    data_id    : PTR_NEW()})

  ndims_notime = ((*i).nxind GT 1) + ((*i).nyind GT 1) + $
    ((*i).nzind GT 1) + ((*i).nvind GT 1) + ((*i).nwind GT 1)
  IF ndims_notime EQ 3 THEN BEGIN
    PRINT, (*diag).name + ' warning: three phase space dimensions ' + $
      'selected, taking only first time step'
    (*i).skip_time = 1
  ENDIF
  IF ndims_notime GT 3 THEN BEGIN

    printerror, 'skipping ' + (*diag).name + $
      ': cannot plot more than three dimensions'
    (*diag).selected = 0
  ENDIF

END

;######################################################################

FUNCTION compute_f, diag, g1, calc_apar=calc_apar

  COMMON global_vars

  i = (*diag).internal_vars

  nkx = par.nkx0
  nky = par.nky0
  nz = par.nz0
  nv = par.nv0
  nw = par.nw0
  ns = gui.out.n_spec_sel

  IF ns NE par.n_spec THEN BEGIN
    PRINT, (*diag).name + $
      ' error: must select all species to get physical result!'
    RETURN, -1
  ENDIF

  v_par = - par.lv + INDGEN(nv) * 2.0D * par.lv / (nv - 1.0D)
  dv = v_par[1] - v_par[0]
  mu = (*i).mu_knots
  B0 = (*series.geom).Bfield
  vTj = SQRT(2.0D*spec.temp/spec.mass)

  kx = 2.0D * !DPI / par.lx * INDGEN(nkx)
  kx[nkx/2+1] = - kx[nkx-nkx/2-1-INDGEN(nkx-nkx/2-1)]
  ky = par.kymin * INDGEN(nky)
  zcoord = !DPI * (2 * INDGEN(nz) / DOUBLE(nz) - 1)
  s_alpha_fac = REBIN(REFORM(par.shat*zcoord-par.amhd*SIN(zcoord),$
    [1,1,nz]),[nkx,nky,nz])
  k_perp = SQRT(REBIN(kx^2,[nkx,nky,nz])+$
    2*s_alpha_fac*REBIN(kx,[nkx,nky,nz])*$
    REBIN(REFORM(ky,[1,nky]),[nkx,nky,nz])+(1+s_alpha_fac^2)*$
    REBIN(REFORM(ky^2,[1,nky]),[nkx,nky,nz]))

  lambda_j = REBIN(REFORM(REBIN(REFORM(vTj*spec.mass/spec.charge,$
    [1,1,ns]),[nz,nw,ns])*SQRT(REBIN(REFORM(mu,[1,nw]),[nz,nw,ns])/$
    REBIN(B0,[nz,nw,ns])),[1,1,nz,nw,ns]),[nkx,nky,nz,nw,ns]) * $
    REBIN(k_perp,[nkx,nky,nz,nw,ns])
  Bessel_J0_lambda_j = $
    BESELJ(TEMPORARY(lambda_j),0,DOUBLE=(1-par.prec_single))

  F0_arg = REBIN(REFORM(-v_par^2,[1,nv]),[nz,nv,nw]) - $
    REBIN(REFORM(REBIN(REFORM(mu,[1,nw]),[nz,nw])*$
    REBIN(B0,[nz,nw]),[nz,1,nw]),[nz,nv,nw])

  F0 = !DPI^(-1.5D) * EXP(TEMPORARY(F0_arg))

  F0_term = TOTAL(TOTAL(REBIN(REFORM(REBIN(REFORM(dv*v_par^2,$
    [1,nv]),[nz,nv,nw])*F0,[1,1,nz,nv,nw]),$
    [nkx,nky,nz,nv,nw,ns])*REBIN(REFORM($
    Bessel_J0_lambda_j^2,[nkx,nky,nz,1,nw,ns]),$
    [nkx,nky,nz,nv,nw,ns]),4)*REBIN(REFORM($
    (*i).mu_weights,[1,1,1,nw]),[nkx,nky,nz,nw,ns]),4)

  A_par_denom = k_perp^2
  FOR n = 0, ns - 1 DO A_par_denom += REBIN(REFORM($
    par.beta*spec[n].charge^2/spec[n].mass*spec[n].dens*!DPI*$
    B0,[1,1,nz]),[nkx,nky,nz]) * F0_term[*,*,*,n]
  F0_term = 0

  curr_term = TOTAL(TOTAL(REBIN(REFORM(dv*v_par,[1,1,1,nv]),$
    [nkx,nky,nz,nv,nw,ns])*g1,4)*Bessel_J0_lambda_j*$
    REBIN(REFORM((*i).mu_weights,[1,1,1,nw]),$
    [nkx,nky,nz,nw,ns]),4)

  A_par_numer = 0.0D
  FOR n = 0, ns - 1 DO A_par_numer += REBIN(REFORM($
    REBIN([0.5D*par.beta*!DPI*spec[n].charge*spec[n].dens*vTj[n]],$
    nz)*B0,[1,1,nz]),[nkx,nky,nz]) * curr_term[*,*,*,n]
  curr_term = 0

  A_par = TEMPORARY(A_par_numer) / TEMPORARY(A_par_denom)

  IF KEYWORD_SET(calc_apar) THEN RETURN, A_par

  A_par_re = par.prec_single ? FLOAT(A_par) : DOUBLE(A_par)
  A_par_im = IMAGINARY(TEMPORARY(A_par))
  A_par_rebin = COMPLEX(REBIN(TEMPORARY(A_par_re),$
    [nkx,nky,nz,nv,nw,ns]),REBIN(TEMPORARY(A_par_im),$
    [nkx,nky,nz,nv,nw,ns]),DOUBLE=(1-par.prec_single))

  A_par_Bessel = TEMPORARY(A_par_rebin) * $
    REBIN(REFORM(TEMPORARY(Bessel_J0_lambda_j),$
    [nkx,nky,nz,1,nw,ns]),[nkx,nky,nz,nv,nw,ns])

  gf_diff = REBIN(REFORM(2.0D*spec.charge/(spec.mass*vTj),[1,1,1,ns]),$
    [nz,nv,nw,ns]) * REBIN(REFORM(v_par,[1,nv]),[nz,nv,nw,ns]) * $
    REBIN(F0,[nz,nv,nw,ns])
  gf_diff = REBIN(REFORM(TEMPORARY(gf_diff),[1,1,nz,nv,nw,ns]),$
    [nkx,nky,nz,nv,nw,ns]) * TEMPORARY(A_par_Bessel)

  RETURN, g1 - TEMPORARY(gf_diff)

END

;######################################################################

PRO checkpoint_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).skip_time AND (series.step_count NE 0) THEN RETURN

  IF (*i).plot_apar THEN BEGIN
    (*i).apar = compute_f(diag,*chpt,/calc_apar)
    RETURN
  ENDIF

  data_ptr = (*i).g2f ? PTR_NEW(compute_f(diag,*chpt)) : chpt

  data_spec = REFORM(FLTARR((*i).nxind,(*i).nyind,(*i).nzind,$
    (*i).nvind,(*i).nwind,gui.out.n_spec_sel),$
    [(*i).nxind,(*i).nyind,(*i).nzind,$
    (*i).nvind,(*i).nwind,gui.out.n_spec_sel],/OVERWRITE)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    xdim = (*i).xavg ? par.nx0 : (*i).nxind
    ydim = (*i).yavg ? par.nky0 : (*i).nyind / (1 + (*i).yfft)
    zdim = (*i).zavg ? par.nz0 : (*i).nzind
    vdim = (*i).vavg ? par.nv0 : (*i).nvind
    wdim = (*i).wavg ? par.nw0 : (*i).nwind

    ; select only required data
    data = REFORM((*data_ptr)$
      [(*i).xind,(*i).yind,(*i).zind,(*i).vind,(*i).wind,isp],$
      [xdim,ydim,zdim,vdim,wdim],/OVERWRITE)

    IF (*i).wavg THEN BEGIN
      FOR w = 0, par.nw0 - 1 DO $
        data[0,0,0,0,w] = data[*,*,*,*,w] * (*i).mu_weights[w]
      data = REFORM(TOTAL(TEMPORARY(data),5),$
        [xdim,ydim,zdim,vdim,1],/OVERWRITE)
      wdim = 1
    ENDIF

    IF (*i).vavg THEN BEGIN
      data = REFORM(TOTAL(TEMPORARY(data),4),$
        [xdim,ydim,zdim,1,wdim],/OVERWRITE)
      vdim = 1
    ENDIF

    IF (*i).yfft THEN BEGIN
      ydim *= 2
      data_ext = par.prec_single ? $
        COMPLEXARR(xdim,ydim,zdim,vdim,wdim,/NOZERO) : $
        DCOMPLEXARR(xdim,ydim,zdim,vdim,wdim,/NOZERO)
      data_ext[0,0,0,0,0] = data[0:xdim/2,*,*,*,*]
      data_ext[*,ydim/2,*,*,*] = 0
      data_ext[0,ydim/2+1,0,0,0] = $
        CONJ(data[0,ydim/2-1-INDGEN(ydim/2-1),*,*,*])
      data_ext[1,ydim/2+1,0,0,0] = $
        CONJ((TEMPORARY(data))[xdim-1-INDGEN(xdim-xdim/2),$
        ydim/2-1-INDGEN(ydim/2-1),*,*,*])
      data = TEMPORARY(data_ext)
      data = FFT(data,DIMENSION=2,$
        /INVERSE,DOUBLE=(1-par.prec_single),/OVERWRITE)
    ENDIF

    IF (*i).xfft THEN BEGIN
      IF (*i).yfft THEN BEGIN
        data_temp = par.prec_single ? $
          COMPLEXARR(xdim,ydim,zdim,vdim,wdim,/NOZERO) : $
          DCOMPLEXARR(xdim,ydim,zdim,vdim,wdim,/NOZERO)
        data_temp[0,0,0,0,0] = data[0:xdim/2-1,*,*,*,*]
        data_temp[xdim/2,*,*,*,*] = 0
        data_temp[xdim/2+1,0,0,0,0] = CONJ((TEMPORARY($
          data))[xdim/2-1-INDGEN(xdim/2-1),*,*,*,*])
        data = TEMPORARY(data_temp)
        data = FFT(data,DIMENSION=1,$
          /INVERSE,DOUBLE=(1-par.prec_single),/OVERWRITE)
      ENDIF ELSE data = FFT(data,DIMENSION=1,$
        /INVERSE,DOUBLE=(1-par.prec_single),/OVERWRITE)
    ENDIF

    ; take absolute square
    data = ABS(TEMPORARY(data))^2
    IF (*i).vavg THEN data /= vdim

    ; x-z Jacobian for global runs
    IF NOT par.x_local AND ((*i).zavg OR (*i).xavg) THEN $
      data = TEMPORARY(data) * REBIN(REFORM($
      (*series.geom).jac_norm[(*i).xind,(*i).zind]/((*i).xind*(*i).zind),$
      [xdim,1,zdim,1,1],/OVERWRITE),$
      [xdim,ydim,zdim,vdim,wdim])

    IF (*i).zavg THEN BEGIN
      ; z Jacobian for local runs
      IF par.x_local THEN data = TEMPORARY(data) * REBIN(REFORM($
        (*series.geom).jac_norm[(*i).zind]/(*i).nzind,$
        [1,1,zdim,1,1],/OVERWRITE),[xdim,ydim,zdim,vdim,wdim])

      data = REFORM(TOTAL(TEMPORARY(data),3),$
        [xdim,ydim,1,vdim,wdim],/OVERWRITE)
      zdim = 1
    ENDIF

    IF (*i).yavg THEN BEGIN
      data[*,0,*,*,*] *= 0.5

      data = 2 * REFORM(TOTAL(TEMPORARY(data),2),$
        [xdim,1,zdim,vdim,wdim],/OVERWRITE)
      ydim = 1
    ENDIF

    IF (*i).xavg THEN BEGIN
      data = REFORM(TOTAL(TEMPORARY(data),1),$
        [1,ydim,zdim,vdim,wdim],/OVERWRITE)
      xdim = 1
    ENDIF

    data_spec[0,0,0,0,0,isp] = TEMPORARY(data)
  ENDFOR ; --- species loop

  IF (*i).g2f THEN PTR_FREE, data_ptr

  (*i).time_id = store_step((*i).time_id,chpt_time)
  (*i).data_id = store_step((*i).data_id,data_spec)

  PRINT, 'loop max: ' + mem_usage(/highw)

END

;######################################################################

PRO checkpoint_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).plot_apar THEN BEGIN ; plot A_par as computed from g1
    temp_sxky = FFT((*i).apar,DIMENSION=1,/INVERSE,$
      DOUBLE=(1-par.prec_single))

    temp_data = par.prec_single ? $
      COMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO) : $
      DCOMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO)

    temp_data[0,0,0] = temp_sxky
    temp_data[0,par.nky0,0] = COMPLEXARR(par.nkx0,1,par.nz0)
    temp_data[0,par.nky0+1,0] = $
      CONJ(TEMPORARY(temp_sxky[*,par.nky0-1-INDGEN(par.nky0-1),*]))

    outdata = TEMPORARY(FFT(temp_data,DIMENSION=2,$
      /INVERSE,DOUBLE=(1-par.prec_single),/OVERWRITE))
    outdata = FLOAT(TEMPORARY(outdata))

    xtitle = '!6x/!7q!6!Dref!N'
    xaxis = ((*i).xind / FLOAT(par.nx0-1) - 0.5) * par.lx
    ytitle = '!6y/!7q!6!Dref!N'
    yaxis = FINDGEN(par.ny0) / par.ny0 * par.ly - par.ly / 2

    PRINT, 'A_par min/max: ', rm0es(MIN(outdata[*,*,par.nz0/2])), $
      rm0es(MAX(outdata[*,*,par.nz0/2]))

    set_output, diag, /ps, coltable=33, ysize=10

    lev_col = CONTOUR_LEVELS(outdata[*,*,par.nz0/2],20)
    CONTOUR, outdata[*,*,par.nz0/2], xaxis, yaxis, COLOR=1, $
      /XSTYLE, /YSTYLE, /ISOTROPIC, /FILL, $
      LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], $
      XTITLE=xtitle, YTITLE=ytitle

    set_output, diag, /reset

    RETURN
  ENDIF

  dims = BYTARR((((*i).nxind GT 1)+((*i).nyind GT 1)+$
    ((*i).nzind GT 1)+((*i).nvind GT 1)+((*i).nwind GT 1)+$
    (series.step_count GT 1))>1)

  ndims = 0
  IF (*i).nxind GT 1 THEN BEGIN
    dims[ndims] = 1
    dimarr = (*i).nxind
    ndims += 1
  ENDIF
  IF (*i).nyind GT 1 THEN BEGIN
    dims[ndims] = 2
    IF ndims NE 0 THEN dimarr = [dimarr,(*i).nyind] ELSE $
      dimarr = (*i).nyind
    ndims += 1
  ENDIF
  IF (*i).nzind GT 1 THEN BEGIN
    dims[ndims] = 3
    IF ndims NE 0 THEN dimarr = [dimarr,(*i).nzind] ELSE $
      dimarr = (*i).nzind
    ndims += 1
  ENDIF
  IF (*i).nvind GT 1 THEN BEGIN
    dims[ndims] = 4
    IF ndims NE 0 THEN dimarr = [dimarr,(*i).nvind] ELSE $
      dimarr = (*i).nvind
    ndims += 1
  ENDIF
  IF (*i).nwind GT 1 THEN BEGIN
    dims[ndims] = 5
    IF ndims NE 0 THEN dimarr = [dimarr,(*i).nwind] ELSE $
      dimarr = (*i).nwind
    ndims += 1
  ENDIF
  IF NOT (*i).skip_time AND (series.step_count GT 1) THEN BEGIN
    dims[ndims] = 6
    IF ndims NE 0 THEN dimarr = [dimarr,series.step_count] ELSE $
      dimarr = series.step_count
    ndims += 1
  ENDIF

  IF ndims NE N_ELEMENTS(dims) THEN IF ndims NE 0 THEN BEGIN
    PRINT, (*diag).name + ' error: dimension mismatch, aborting'
    RETURN
  ENDIF

  IF ndims EQ 0 THEN BEGIN
    time = store_step((*i).time_id,/get,/reg_array)
    data = store_step((*i).data_id,/get,/reg_array)
    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      sp = (*gui.out.spec_select)[isp]
      PRINT, 'Total ' + rm0es(spec[sp].name) + $
        ' average (x,y,z,vpar,mu,time): ' + $
        rm0es(TOTAL(data[*,*,*,*,*,isp,*])/series.step_count)
    ENDFOR

    RETURN
  ENDIF

  t_is_axis = (WHERE(dims EQ 6) NE -1) AND (ndims NE 3)
  t_is_label = (WHERE(dims EQ 6) NE -1) AND NOT t_is_axis

  ; for three dimensions (or 1 + time), create multiple plots
  n_plots = (((ndims < 2) - 1) * (dimarr[ndims-1] - 1)) + 1
  IF NOT t_is_label AND (ndims EQ 2) THEN n_plots = 1

  ; line versus contour plotting
  plot_dim = 1 + ((ndims GT 1) AND NOT t_is_axis)

  CASE dims[0] OF
    1 : BEGIN
          xtitle = par.x_local + (*i).xfft EQ 1 ? $
            '!6k!Dx!N!7q!6!Dref!N' : '!6x/!7q!6!Dref!N'
          xaxis = par.x_local + (*i).xfft EQ 1 ? $
            (*par.kx)[(*i).xind,0] : $
            ((*i).xind / FLOAT(par.nx0-1) - 0.5) * par.lx
        END
    2 : BEGIN
          xtitle = (*i).yfft ? '!6y/!7q!6!Dref!N' : '!6k!Dy!N!7q!6!Dref!N'
          xaxis = (*i).yfft ? FINDGEN(par.ny0) / par.ny0 * $
            par.ly - par.ly / 2 : (*par.ky)[(*i).yind+par.ky0_ind]
        END
    3 : BEGIN
          xtitle = '!6z/L!Dref!N'
          xaxis = (*par.z)[(*i).zind]
        END
    4 : BEGIN
          xtitle = '!6v!D!9#!6!N/v!DTj!N'
          xaxis = (*i).vind * (2.0 * par.lv / (par.nv0 - 1)) - par.lv
        END
    5 : BEGIN
          xtitle = '!7l!6/(T!Dj0!N/B!Dref!N)'
          xaxis = (*i).mu_knots[(*i).wind]
        END
    6 : BEGIN
          xtitle = '!6t/(L!Dref!N/c!Dref!N)'
          xaxis = LINDGEN(series.step_count)
        END
    ELSE :
  ENDCASE

  IF plot_dim EQ 2 THEN BEGIN
    CASE dims[1] OF
      2 : BEGIN
            ytitle = (*i).yfft ? '!6y/!7q!6!Dref!N' : '!6k!Dy!N!7q!6!Dref!N'
            yaxis = (*i).yfft ? FINDGEN(par.ny0) / par.ny0 * $
              par.ly - par.ly / 2 : (*par.ky)[(*i).yind+par.ky0_ind]
          END
      3 : BEGIN
            ytitle = '!6z/L!Dref!N'
            yaxis = (*par.z)[(*i).zind]
          END
      4 : BEGIN
            ytitle = '!6v!D!9#!6!N/v!DTj!N'
            yaxis = (*i).vind * (2.0 * par.lv / (par.nv0 - 1)) - par.lv
          END
      5 : BEGIN
            ytitle = '!7l!6/(T!Dj0!N/B!Dref!N)'
            yaxis = (*i).mu_knots[(*i).wind]
          END
      ELSE :
    ENDCASE
  ENDIF

  time = store_step((*i).time_id,/get,reg_array=t_is_axis)
  data = store_step((*i).data_id,/get,reg_array=t_is_axis)

  IF t_is_axis THEN data_plot = TEMPORARY(data)

  title = '!9!!!6' + ((*i).g2f ? 'f' : 'g') + '!9!!!6!U2!N'
  IF (*i).xavg THEN title += ', x avg.' ELSE IF (*i).nxind EQ 1 THEN $
    title += (par.x_local + (*i).xfft EQ 1) ? $
    ', k!Dx!N='+rm0es((*par.kx)[(*i).xind,0],prec=3) : $ 
    ', x='+rm0es(((*i).xind / FLOAT(par.nx0-1) - 0.5) * par.lx,prec=3)
  IF (*i).yavg THEN title += ', y avg.' ELSE IF (*i).nyind EQ 1 THEN $
    title += (*i).yfft ? rm0es(FINDGEN(par.ny0) / par.ny0 * $
    par.ly - par.ly / 2,prec=3) : $
    ', k!Dy!N=' + rm0es((*par.ky)[(*i).yind+par.ky0_ind],prec=3)
  IF (*i).zavg THEN title += ', z avg.' ELSE IF (*i).nzind EQ 1 THEN $
    title += ', z=' + rm0es((*par.z)[(*i).zind],prec=3)
  IF (*i).vavg THEN title += ', !6v!D!9#!6!N avg.' ELSE $
    IF (*i).nvind EQ 1 THEN title += ', !6v!D!9#!6!N=' + $
    rm0es((*i).vind*(2.0*par.lv/(par.nv0+1))-par.lv,prec=3)
  IF (*i).wavg THEN title += ', !7l!6 avg.' ELSE $
    IF (*i).nwind EQ 1 THEN title += ', !7l!6=' + $
    rm0es((*i).mu_knots[(*i).wind],prec=3)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    IF plot_dim EQ 1 THEN set_output, diag, sp, /ps ELSE $
      set_output, diag, sp, /ps, coltable=33

    FOR plotnr = 0, n_plots - 1 DO BEGIN
      IF NOT t_is_axis THEN BEGIN
        IF t_is_label THEN data_plot = $
          REFORM((*data[plotnr])[*,*,*,*,*,isp],$
          dimarr[0:ndims-2],/OVERWRITE) ELSE BEGIN

          IF n_plots GT 1 THEN data_plot = $
            REFORM((REFORM((*data[0])[*,*,*,*,*,isp],dimarr,$
            /OVERWRITE))[(PRODUCT(dimarr[0:(ndims-2)>0])*plotnr):$
            (PRODUCT(dimarr[0:(ndims-2)>0])*(plotnr+1)-1)],$
            [dimarr[0:(ndims-2)>0]]) ELSE $
          data_plot = REFORM((*data[0])[*,*,*,*,*,isp],dimarr,$
            /OVERWRITE)
        ENDELSE
      ENDIF

      IF plot_dim EQ 1 THEN BEGIN
        PLOT, xaxis, data_plot, COLOR=1, /XSTYLE, XTITLE=xtitle, $
          TITLE=title, YLOG=(*i).log
      ENDIF ELSE BEGIN
        lev_col = contour_levels([MIN(data_plot),MAX(data_plot)],$
          30,log=(*i).log)
        
        CONTOUR, data_plot, xaxis, yaxis, COLOR=1, /XSTYLE, /YSTYLE, $
          XTITLE=xtitle, YTITLE=ytitle, TITLE=title, /FILL, $
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1],$
          YMARGIN=[10,3]
        
        plot_colorbar, lev_col, prec=3, log=(*i).log, $
          POSITION=[0.2,0.075,0.95,0.125]
      ENDELSE

      IF dims[ndims-1] NE 6 THEN $
        plot_info_str, diag, time=(*time[0])
    ENDFOR

    set_output, diag, sp, /reset

    IF plot_dim EQ 1 THEN BEGIN ; data output in 1D case
      dat_xtitles = [par.x_local+(*i).xfft EQ 1 ? 'kx' : 'x',$
        (*i).yfft ? 'y' : 'ky','z','vpar','mu','t']
      datxtitle = dat_xtitles[dims[0]-1]

      set_output, diag, sp, $
        header=[datxtitle,'|'+((*i).g2f ? 'f' : 'g')+'|^2'], $
        dat=[[xaxis],[data_plot]]
    ENDIF
  ENDFOR ; --- species loop

  IF NOT t_is_axis THEN PTR_FREE, data
  PTR_FREE, time

END
