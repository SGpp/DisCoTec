FUNCTION spectra_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'Amplitude spectra (ky/kx)',$
    help_text : ['Plots amplitude spectra in ky/kx space, '+$
        'averaged over the remaining spatial coordinates and time.'],$
    ext_vars  : [$		 
        ['vars','0','index of the fluctuating quantity to be plotted '+$
         '(see variable list); default: all; negative values for '+$
         'multiplication of variable (-var-1) with k_[x,y]'],$
        ['kxind','0','kx index at which to evaluate ky spectrum; '+$
         'default: average (-1)'],$
        ['kyind','0','ky index at which to evaluate kx spectrum; '+$
         'default: average (-1)'],$
        ['zind','0','z index at which to evaluate spectra; '+$
         'default: average (-1); note: additional ktheta and k2 '+$
         'spectra will be plotted if a single zind is given'],$
;        ['minmax','0','plot range; default: [1e-4,1e12]'],$
        ['fit','0','relative fit range (2d array); default: no fit']]}

END

;######################################################################

PRO spectra_init, diag

  COMMON global_vars

  e = (*diag).external_vars

;  IF NOT KEYWORD_SET(*e.minmax) THEN *e.minmax = [1e-4,1e12]
  IF N_ELEMENTS(*e.vars) LT 1 THEN $
    IF par.beta GT 0 THEN *e.vars = INDGEN(par.n_moms+par.n_fields) $
    ELSE *e.vars = [0,INDGEN(par.n_moms)+par.n_fields]

  ; look for negative variables, set k multiplication
  neg_var = WHERE(*e.vars LT 0)
  kmult_var = BYTARR(N_ELEMENTS(*e.vars))
  IF neg_var[0] NE -1 THEN BEGIN
    (*e.vars)[neg_var] = - (*e.vars)[neg_var] - 1
    kmult_var[neg_var] = 1B
  ENDIF

  IF N_ELEMENTS(*e.kxind) LT 1 THEN $
    IF (par.kx0_ind EQ 0) THEN *e.kxind = -1 $
    ELSE *e.kxind = 0
  IF N_ELEMENTS(*e.kyind) LT 1 THEN $
    IF par.adapt_lx THEN *e.kyind = par.ky0_ind $
    ELSE *e.kyind = -1
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = -1
  
  IF ((NOT par.x_local) AND ((*e.zind) NE -1)) THEN BEGIN
     printerror, 'zind NE -1 not supported for x global; '+$
                 'switching to zind=-1'
     (*e.zind) = -1
  ENDIF

  IF par.adapt_lx THEN BEGIN
    PRINT, (*diag).name + ' warning: using only kx=0 mode ' + $
      '(without conn.) for ky spectrum'
    *e.kxind = 0
  ENDIF

  IF N_ELEMENTS(*e.fit) NE 2 THEN fit = [0.0,0.0] ELSE fit = *e.fit

  i = set_internal_vars(diag,{$
    vars        : *e.vars,$
    kmult_var   : kmult_var,$
    kxind       : *e.kxind,$
    kyind       : *e.kyind,$
    zind        : *e.zind[0],$
    xsumspec_id : PTR_NEW(),$
    ysumspec_id : PTR_NEW(),$
    fit     	: fit,$
    minmax      : [1e-4,1e12]}) ;*e.minmax})

  fft_format, kxky=*e.vars

END

;######################################################################

PRO spectra_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkxind = N_ELEMENTS((*i).kxind)
  nkyind = N_ELEMENTS((*i).kyind)
  nkxpos = par.nkx0 / 2 > 1
  nvars = N_ELEMENTS((*i).vars)
  nzind = ((*i).zind EQ -1)?par.nz0:1
  zind = ((*i).zind EQ -1)?INDGEN(par.nz0):(*i).zind

  xsumspec = FLTARR(nkxpos,nkyind,nvars,gui.out.n_spec_sel,/NOZERO)
  ysumspec = FLTARR(par.nky0,nkxind,nvars,gui.out.n_spec_sel,/NOZERO)
  xtempdata = REFORM(FLTARR(nkxpos,nzind,/NOZERO),[nkxpos,nzind])
  ytempdata = REFORM(FLTARR(par.nky0,nzind,/NOZERO),[par.nky0,nzind])

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    FOR ivar = 0, nvars - 1 DO BEGIN
      FOR iky = 0, nkyind - 1 DO BEGIN
        IF (*i).kyind[iky] EQ -1 THEN BEGIN ; x spectrum averaged over ky
          IF par.x_local THEN BEGIN
            xtempdata[0,*] = 2.0 * TOTAL(ABS($
              (*mom[isp,(*i).vars[ivar]].kxky)[0,1:par.nky0-1,zind])^2,2) + $
               REFORM(TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxky)[0,0,zind])^2,2),[1,nzind])
            xtempdata[1:nkxpos-1,*] = TOTAL(ABS($
              (*mom[isp,(*i).vars[ivar]].kxky)[1:nkxpos-1,1:par.nky0-1,zind])^2,2) + $
              REVERSE(TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxky)$
              [par.nkx0-nkxpos+1:par.nkx0-1,1:par.nky0-1,zind])^2,2),1) + $
              TOTAL(REFORM(ABS((*mom[isp,(*i).vars[ivar]].kxky)[1:nkxpos-1,0,zind])^2,$
              [nkxpos-1,1,nzind]),2)
          ENDIF ELSE BEGIN
            xtempdata[0,*] = 2.0 * (TOTAL(ABS($
              (*mom[isp,(*i).vars[ivar]].kxky)[0,1:par.nky0-1,zind])^2,2) + $
              TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxky)[0,0,zind])^2,2)) * $
              (*series.geom).jac_norm[0,zind]
            xtempdata[1:nkxpos-1,*] = TOTAL(ABS($
              (*mom[isp,(*i).vars[ivar]].kxky)[1:nkxpos-1,1:par.nky0-1,zind])^2*$
              REBIN(REFORM((*series.geom).jac_norm[1:nkxpos-1,zind],$
              [nkxpos-1,1,par.nz0]),[nkxpos-1,par.nky0-1,par.nz0]),2) + $
              REVERSE(TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxky)$
              [par.nkx0-nkxpos+1:par.nkx0-1,1:par.nky0-1,zind])^2*$
              REBIN(REFORM((*series.geom).jac_norm[par.nkx0-nkxpos+1:par.nkx0-1,zind],$
              [nkxpos-1,1,par.nz0]),[nkxpos-1,par.nky0-1,par.nz0]),2),1) + $
              TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxky)[1:nkxpos-1,0,zind])^2*$
              (*series.geom).jac_norm[1:nkxpos-1,zind],2)
          ENDELSE
        ENDIF ELSE BEGIN ; x spectrum for certain ky modes
          IF (par.ky0_ind EQ 0) AND ((*i).kyind[iky] EQ 0) THEN BEGIN ; --- zero ky mode
            xtempdata[0,*] = ABS((*mom[isp,(*i).vars[ivar]].kxky)[0,0,zind])^2
            IF nkxpos GT 1 THEN xtempdata[1:nkxpos-1,zind] = $
              ABS((*mom[isp,(*i).vars[ivar]].kxky)[1:nkxpos-1,0,zind])^2
          ENDIF ELSE BEGIN
            xtempdata[0,*] = ABS((*mom[isp,(*i).vars[ivar]].kxky)[0,(*i).kyind[iky],zind])^2
            IF nkxpos GT 1 THEN xtempdata[1:nkxpos-1,zind] = $
              ABS((*mom[isp,(*i).vars[ivar]].kxky)[1:nkxpos-1,(*i).kyind[iky],zind])^2 + $
              ABS(REVERSE((*mom[isp,(*i).vars[ivar]].kxky)$
              [par.nkx0-nkxpos+1:par.nkx0-1,(*i).kyind[iky],zind],1))^2
          ENDELSE
        ENDELSE

        IF (par.x_local AND (nzind EQ par.nz0)) THEN xtempdata *= REBIN(REFORM($
          (*series.geom).jac_norm,[1,nzind]),[nkxpos,nzind])

        xsumspec[*,iky,ivar,isp] = nzind GT 1 ? $
          TOTAL(xtempdata,2) / par.nz0 : xtempdata
      ENDFOR ; --- iky loop
      FOR ikx = 0, N_ELEMENTS((*i).kxind) - 1 DO BEGIN
        IF (*i).kxind[ikx] EQ -1 THEN ytempdata[0,0] = $
          TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxky)[*,*,zind])^2,1) ELSE $
          ytempdata[0,0] = $
          ABS(REFORM((*mom[isp,(*i).vars[ivar]].kxky)[(*i).kxind[ikx],*,zind]))^2

        IF (par.x_local AND (nzind EQ par.nz0)) THEN ytempdata *= REBIN(REFORM($
          (*series.geom).jac_norm,[1,nzind]),[par.nky0,nzind])

        ysumspec[*,ikx,ivar,isp] = nzind GT 1 ? $
          TOTAL(ytempdata,2) / par.nz0 : ytempdata
      ENDFOR ; --- ikx loop
    ENDFOR ; --- ivar
  ENDFOR ; --- species loop

  (*i).xsumspec_id = time_avg((*i).xsumspec_id,xsumspec,mom_time,$
                              fwd=gui.out.res_steps)
  (*i).ysumspec_id = time_avg((*i).ysumspec_id,ysumspec,mom_time,$
                              fwd=gui.out.res_steps)

END

;######################################################################

PRO spectra_plot, xdata, ydata, log=log, position=position, $
  title=title, xtitle=xtitle, ytitle=ytitle, fit=fit, $
  colorarr=colorarr, linestyle=linestyle, minmax=minmax, vars=vars, $
  kmult=kmult, label=label, xname=xname, warn=warn,unitstr=unitstr

  IF NOT KEYWORD_SET(label) THEN label = xtitle
  IF N_ELEMENTS(warn) LT 1 THEN warn = 0

  n_vars = N_ELEMENTS(ydata[0,*])
  IF N_ELEMENTS(kmult) NE n_vars THEN kmult = BYTARR(n_vars)
  IF N_ELEMENTS(unitstr) LT 1 THEN unitstr = ''

  ndata = N_ELEMENTS(xdata)

  mymax = MAX(ydata, min=mymin)
  IF mymax GT minmax[1] THEN BEGIN
    mymax = minmax[1]
    IF warn THEN PRINT, 'suppressed values > ' + rm0es(mymax) + $
      ' in ' + xname + ' spectra'
  ENDIF

  IF log THEN BEGIN
    posv = WHERE(xdata GT 0.0)
    tickformat = 'betterticks'
    IF mymin LT minmax[0] THEN BEGIN
      mymin = minmax[0]
      IF warn THEN PRINT, 'suppressed values < ' + rm0es(mymin) + $
        ' in ' + xname + ' spectra'
    ENDIF
  ENDIF ELSE BEGIN
    posv = INDGEN(ndata)
    tickformat = ''
  ENDELSE
  
  PLOT, xdata[posv], ydata[posv,0], COLOR=1, $
    XLOG=log, XTICKFORMAT=tickformat, XTITLE=xtitle+' '+unitstr,$
    YLOG=log, YTICKFORMAT=tickformat, YTITLE=ytitle,$
    YMARGIN=[0,1], YRANGE=[mymin,mymax], XTICKLEN=1.0, $ ;/YSTYLE, $
    /XSTYLE, TITLE=title, POSITION=position, $
    CHARSIZE=1.0,/NODATA

  FOR ivar = 0, n_vars - 1 DO BEGIN
    ydata_mod = (kmult[ivar] ? xdata[posv]^2 : 1) * ydata[posv,ivar]
    OPLOT, xdata[posv], ydata_mod, COLOR=colorarr[ivar], $
      LINESTYLE=linestyle[ivar], PSYM=-2*kmult[ivar]
  ENDFOR

  IF log AND (TOTAL(fit) NE 0) THEN BEGIN
    PRINT, 'fit results for ' + rm_idl_fmt(title) + ' in range ' + $
      xname + '=[' + rm0es(xdata[fit[0]*ndata]) + ',' + $
      rm0es(xdata[fit[1]*ndata]) + ']:'
    FOR ivar = 0, n_vars - 1 DO BEGIN
      fitres = LINFIT(ALOG10(xdata[fit[0]*ndata:fit[1]*ndata]),$
        ALOG10(ydata[fit[0]*ndata:fit[1]*ndata,ivar]))
      OPLOT, xdata[0.75*fit[0]*ndata:(1.25*fit[1]<1.0)*ndata],$
        10^fitres[0]*xdata[0.75*fit[0]*ndata:(1.25*fit[1]<1.0)*ndata]^fitres[1],$
        COLOR=colorarr[ivar], linestyle=3
      xmid =10^(0.5*(ALOG10(xdata[fit[1]*ndata])+ALOG10(xdata[fit[0]*ndata])))
      XYOUTS, xmid, 10^(0.5*fitres[0])*xmid^fitres[1], '!9'+STRING(63B)+' !6'+$
        label+'!E'+rm0es(fitres[1],prec=2)+'!N',color=colorarr[ivar]

      PRINT, '<|' + get_var_string(vars[ivar]) + '|^2> ~ ' + $
        xname + '^' + rm0es(fitres[1],prec=2)
    ENDFOR
  ENDIF

END

;######################################################################

PRO spectra_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkxind = N_ELEMENTS((*i).kxind)
  nkyind = N_ELEMENTS((*i).kyind)
  nvars = N_ELEMENTS((*i).vars)
  nsp = gui.out.n_spec_sel  
  nkxpos = par.nkx0 / 2 > 1
  n_res_steps = gui.out.res_steps * (series.step_count - 1) + 1

  xsumspec = REFORM(time_avg((*i).xsumspec_id,/avg,fwd=gui.out.res_steps,$
    tarr=time),[nkxpos,nkyind,nvars,nsp,n_res_steps])
  ysumspec = REFORM(time_avg((*i).ysumspec_id,/avg,fwd=gui.out.res_steps,$
    tarr=time),[par.nky0,nkxind,nvars,nsp,n_res_steps])

  ; this part ensures that variables have always the same color
  ; (if par.n_fields = const.)
  colorarr = INDGEN(par.n_moms+par.n_fields) + 1
  linestyle = INTARR(par.n_moms+par.n_fields)
  FOR var = 6, par.n_moms - 1 DO BEGIN
     colorarr[var+par.n_fields] = colorarr[par.n_fields+(var MOD 6)]
     linestyle[var+par.n_fields] = var / 6
  ENDFOR

  colorarr = colorarr[(*i).vars]
  linestyle = linestyle[(*i).vars]
  psym = INTARR(nvars)
  kx_str = STRARR(nvars)
  ky_str = STRARR(nvars)

  IF TOTAL((*i).kmult_var) GT 0 THEN BEGIN
    kmult_inds = WHERE((*i).kmult_var)
    psym[kmult_inds] = -2
    kx_str[kmult_inds] = '!6k!Dx!N'
    ky_str[kmult_inds] = '!6k!Dy!N'
  ENDIF

  rho_str = ' ' + get_var_string(/rhostr)
  zstr=((*i).zind EQ -1)?'z avg.':'z = '+rm0es((*par.z)[(*i).zind],prec=3)

  ky_axes = (*series.ky)
  ky_title = '!6k!Dy!N'
  ky_unit = rho_str
  kx = (*series.kx)[0:par.nkx0/2-1,0]
  kx_unit = rho_str
  
  IF (((*i).zind NE -1) AND par.x_local) THEN BEGIN
     rho_cm = SQRT(par.mref*par.Tref/par.Qref) / par.Bref / $
        SQRT(series.mref*series.Tref/series.Qref) * series.Bref*100.

     k2_cm = get_k2_fac((*i).zind)*(*series.ky)/rho_cm
     ktheta_cm = get_ktheta_fac((*i).zind)*(*series.ky)/rho_cm
     ky_axes = [[ky_axes],[k2_cm],[ktheta_cm]]
     ky_title = [ky_title,'!6k!D2!N','!6k!D!7h!6!N']
     ky_unit = [ky_unit,'/ cm!U-1!N','/ cm!U-1!N']
  ENDIF

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, multi=[0,1,2], coltable=41

    FOR ikx = 0, N_ELEMENTS((*i).kxind) - 1 DO BEGIN
      IF (*i).kxind[ikx] EQ -1 THEN kx_str = 'k!Dx!N avg.' $
      ELSE kx_str = 'k!Dx!N' + get_var_string(/rhostr) + '=' + $
                    rm0es(kx[(*i).kxind[ikx]],prec=3)

      maxy = MAX(ysumspec[*,ikx,*,isp,*],MIN=miny)
      logminy = miny
      IF maxy GT (*i).minmax[1] THEN BEGIN
        maxy = (*i).minmax[1]
        PRINT, 'suppressed values > ' + rm0es(maxy) + ' in ky spectra'
      ENDIF
      IF miny LT (*i).minmax[0] THEN BEGIN
        logminy = (*i).minmax[0]
        PRINT, 'suppressed values < ' + rm0es(logminy) + ' in ky spectra'
      ENDIF
 
      ; --- ky spectra ---
      FOR ky_type = 0, N_ELEMENTS(ky_axes[0,*])-1 DO BEGIN
       ; plot k_y spectrum (logarithmic)
       title = ky_title[ky_type]+' spectrum (' + kx_str + ','+zstr+')'
   
       FOR n = 0, n_res_steps - 1 DO BEGIN 
        IF par.nky0 GT (1 + (par.ky0_ind EQ 0)) THEN BEGIN    
          spectra_plot, ky_axes[*,ky_type], REFORM(ysumspec[*,ikx,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle=ky_title[ky_type], minmax=(*i).minmax, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='ky', $
            position=[0.15,0.596,0.95,0.9615], /log, /warn, $
            unitstr=ky_unit[ky_type]

          linplot_end = 0.538
          title = ''
        ENDIF ELSE linplot_end = 0.9615

        ; plot k_y spectrum (linear)
        IF par.nky0 GT 1 THEN BEGIN
          spectra_plot, ky_axes[*,ky_type], REFORM(ysumspec[*,ikx,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle=ky_title[ky_type], minmax=(*i).minmax, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='ky', $
            position=[0.15,0.231,0.95,linplot_end], log=0, $
            unitstr=ky_unit[ky_type]

          plot_legend, colorarr, ky_str + get_var_string((*i).vars,/fancy), $
            x_offset=[0.15,0.0,0.95,0.12], linestyles=linestyle, psym=psym
        ENDIF

        plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)

        IF (ky_type EQ 0) THEN BEGIN
         commentline = 'kxind = '+$
           rm0es((*i).kxind[ikx])+', '+zstr+(gui.out.res_steps ? ', time = '+rm0es(time[n]) : '')
         IF ((*i).zind NE -1) THEN commentline=[commentline,'k2fac_cm = '+$
            rm0es(get_k2_fac((*i).zind)/rho_cm)+', '+'kthetafac_cm = '+$
            rm0es(get_ktheta_fac((*i).zind)/rho_cm)]
         IF N_ELEMENTS(ky) EQ 1 THEN BEGIN
           set_output, diag, sp, header=['ky',get_var_string((*i).vars)],append=(ikx+n GT 0),$
             dat=[ky_axes[*,0],REFORM(ysumspec[*,ikx,*,isp,n])], commentline=commentline
         ENDIF ELSE BEGIN
           set_output, diag, sp, header=['ky',get_var_string((*i).vars)],append=(ikx+n GT 0), $
             dat=[[ky_axes[*,0]],[REFORM(ysumspec[*,ikx,*,isp,n])]], commentline=commentline
         ENDELSE
        ENDIF
       ENDFOR ; ky type
      ENDFOR ; time loop
    ENDFOR ; kx loop
  
    ; --- kx spectra ---
    FOR iky = 0, N_ELEMENTS((*i).kyind) - 1 DO BEGIN
      kx = (*series.kx)[0:par.nkx0/2-1,iky]

      IF (*i).kyind[iky] EQ -1 THEN ky_str = 'k!Dy!N avg.' $
        ELSE ky_str = 'k!Dy!N' + get_var_string(/rhostr) + '=' + $
        rm0es(ky[(*i).kyind[iky]],prec=3)

      ; plot k_x spectrum (logarithmic)
      title = '!6k!Dx!N spectrum (' + ky_str + ',' + zstr + ')'

      FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
        IF par.nkx0 / 2 GT 1 THEN BEGIN
          spectra_plot, kx, REFORM(xsumspec[*,iky,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle='!6k!Dx!N', minmax=(*i).minmax, unitstr=kx_unit, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='kx', $
            position=[0.15,0.596,0.95,0.9615], /log, /warn

          linplot_end = 0.538
          title = ''
        ENDIF ELSE linplot_end = 0.9615

        ; plot time k_x spectrum (linear)
        IF par.nkx0 / 2 GT 1 THEN BEGIN
          spectra_plot, kx, REFORM(xsumspec[*,iky,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle='!6k!Dx!N', minmax=(*i).minmax, unitstr=kx_unit, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='kx', $
            position=[0.15,0.231,0.95,linplot_end], log=0

          plot_legend, colorarr, get_var_string((*i).vars,/fancy), $
            x_offset=[0.15,0.0,0.95,0.12], linestyles=linestyle, psym=psym
        ENDIF

        plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)

        set_output, diag, sp, header=['kx',get_var_string((*i).vars)], $
                    dat=[[kx],[REFORM(xsumspec[*,iky,*,isp,n])]], $
                    commentline='kyind = '+rm0es((*i).kyind[iky])+','+$
                    zstr+(gui.out.res_steps ? ', time = '+rm0es(time[n]) : ''),$
                    /append
        ENDFOR ; time loop
    ENDFOR ; ky loop

    set_output, diag, sp, /reset
  ENDFOR ; species loop

END
