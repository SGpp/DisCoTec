FUNCTION spectra_y_info

  RETURN, {$
    type      : 'mom global',$ ;'mom_uni',$
    title     : 'Amplitude spectra (ky/kx) (y-global)',$
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
        ['minmax','0','plot range; default: [1e-4,1e12]'],$
        ['fit','0','relative fit range (2d array); default: no fit']]}

END

;######################################################################

PRO spectra_y_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF NOT KEYWORD_SET(*e.minmax) THEN *e.minmax = [1e-4,1e12]
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
    xsumspec_id : PTR_NEW(),$
    ysumspec_id : PTR_NEW(),$
    fit     	: fit,$
    minmax      : *e.minmax})

  fft_format, kxsy=*e.vars

END

;######################################################################

PRO spectra_y_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkxind = N_ELEMENTS((*i).kxind)
  nkyind = N_ELEMENTS((*i).kyind)
  nkxpos = par.nkx0 / 2 + 1 > 1
  nvars = N_ELEMENTS((*i).vars)

  xsumspec = FLTARR(nkxpos,nkyind,nvars,gui.out.n_spec_sel,/NOZERO)
  ysumspec = FLTARR(par.nky0,nkxind,nvars,gui.out.n_spec_sel,/NOZERO)
  yfftdata = FLTARR(nkxpos,par.ny0,par.nz0,/NOZERO)
  xtempdata = FLTARR(nkxpos,par.nz0,/NOZERO)
  ytempdata = FLTARR(par.nky0,/NOZERO)

                                ; Computes 'y-spectra'.  Note that the
                                ; flux surface average performed here
                                ; is not correct as the jacobian is
                                ; not included.  Still useful as a
                                ; rough check on the spectrum in y.

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    FOR ivar = 0, nvars - 1 DO BEGIN
      FOR ikx = 0, nkxind - 1 DO BEGIN

           yfftdata[0:nkxpos-1,*,0:par.nz0-1] = FFT((*mom[isp,(*i).vars[ivar]].kxsy),DIMENSION=2)

            ytempdata[0] = 2.0 * TOTAL(TOTAL(ABS($
              yfftdata[0,1:par.nky0-1,*])^2,3)) + $
              TOTAL(ABS(yfftdata[0,0,*])^2,3)

            FOR iy = 1 , par.nky0-1 DO BEGIN

            ytempdata[iy] = TOTAL(TOTAL(ABS($
              yfftdata[*,iy,*])^2,1),2) + $
           TOTAL(TOTAL(ABS(yfftdata[*,par.ny0-iy,*])^2,1),2)

           ENDFOR

        ysumspec[*,ikx,ivar,isp] = ytempdata[*]


      ENDFOR 

      ; Calculates kx spectrum.  Here, the flux surface average
      ; is performed correctly as the jacobian is available in the 
      ; format required.
    
      FOR iky = 0, N_ELEMENTS((*i).kyind) - 1 DO BEGIN
        IF (*i).kyind[iky] EQ -1 THEN xtempdata[0,0] = $
          TOTAL(ABS((*mom[isp,(*i).vars[ivar]].kxsy))^2,2) ELSE $
          xtempdata[0,0] = $
          ABS(REFORM((*mom[isp,(*i).vars[ivar]].kxsy)[*,(*i).kyind[iky],*]))^2

          xsumspec[*,iky,ivar,isp] = TOTAL(xtempdata,2) / par.nz0
          
     ENDFOR                     ; --- iky loop
    ENDFOR ; --- ivar
  ENDFOR ; --- species loop

  (*i).xsumspec_id = time_avg((*i).xsumspec_id,xsumspec,mom_time,$
                              fwd=gui.out.res_steps)
  (*i).ysumspec_id = time_avg((*i).ysumspec_id,ysumspec,mom_time,$
                              fwd=gui.out.res_steps)

END

;######################################################################

PRO spectra_y_plot, xdata, ydata, log=log, position=position, $
  title=title, xtitle=xtitle, ytitle=ytitle, fit=fit, $
  colorarr=colorarr, linestyle=linestyle, minmax=minmax, vars=vars, $
  kmult=kmult, label=label, xname=xname, warn=warn

  IF NOT KEYWORD_SET(label) THEN label = xtitle
  IF N_ELEMENTS(warn) LT 1 THEN warn = 0

  n_vars = N_ELEMENTS(ydata[0,*])
  IF N_ELEMENTS(kmult) NE n_vars THEN kmult = BYTARR(n_vars)

  rho_str = ' ' + get_var_string(/rhostr)

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
    XLOG=log, XTICKFORMAT=tickformat, XTITLE=xtitle+rho_str,$
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

PRO spectra_y_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkxind = N_ELEMENTS((*i).kxind)
  nkyind = N_ELEMENTS((*i).kyind)
  nvars = N_ELEMENTS((*i).vars)
  nsp = gui.out.n_spec_sel  
  nkxpos = par.nkx0 / 2 + 1 > 1
  n_res_steps = gui.out.res_steps * (series.step_count - 1) + 1

 ; HELP, time_avg((*i).xsumspec_id,/avg,fwd=gui.out.res_steps,$
 ;   tarr=time)

  print, nkxind, nkyind, nvars, nsp, n_res_steps

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
    ky_str[kmult_inds] = '!6k!Dx!N'
  ENDIF

  rho_str = ' ' + get_var_string(/rhostr)

  ky = (*series.ky)
  kx = (*series.kx)[0:par.nkx0/2-1,0]

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
      ; plot time averaged k_y spectrum (logarithmic)
      title = '!6k!Dy!N spectrum (' + kx_str + ')'
   
      FOR n = 0, n_res_steps - 1 DO BEGIN 
        IF par.nky0 GT (1 + (par.ky0_ind EQ 0)) THEN BEGIN    
          spectra_y_plot, ky, REFORM(ysumspec[*,ikx,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle='!6k!Dy!N', minmax=(*i).minmax, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='ky', $
            position=[0.15,0.596,0.95,0.9615], /log, /warn

          linplot_end = 0.538
          title = ''
        ENDIF ELSE linplot_end = 0.9615

        ; plot time averaged k_y spectrum (linear)
        IF par.nky0 GT 1 THEN BEGIN
          spectra_y_plot, ky, REFORM(ysumspec[*,ikx,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle='!6k!Dy!N', minmax=(*i).minmax, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='ky', $
            position=[0.15,0.231,0.95,linplot_end], log=0

          plot_legend, colorarr, ky_str + get_var_string((*i).vars,/fancy), $
            x_offset=[0.15,0.0,0.95,0.12], linestyles=linestyle, psym=psym
        ENDIF

        plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
      ENDFOR ; time loop
    ENDFOR ; kx loop
  
    ; --- kx spectra ---
    FOR iky = 0, N_ELEMENTS((*i).kyind) - 1 DO BEGIN
      kx = (*series.kx)[0:par.nkx0/2-1,iky]

      IF (*i).kyind[iky] EQ -1 THEN ky_str = 'k!Dy!N avg.' $
        ELSE ky_str = 'k!Dy!N' + get_var_string(/rhostr) + '=' + $
        rm0es(ky[(*i).kyind[iky]],prec=3)

      ; plot time averaged k_x spectrum (logarithmic)
      title = '!6k!Dx!N spectrum (' + ky_str + ')'

      FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
        IF par.nkx0 / 2 GT 1 THEN BEGIN
          spectra_y_plot, kx, REFORM(xsumspec[*,iky,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle='!6k!Dx!N', minmax=(*i).minmax, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='kx', $
            position=[0.15,0.596,0.95,0.9615], /log, /warn

          linplot_end = 0.538
          title = ''
        ENDIF ELSE linplot_end = 0.9615

        ; plot time averaged k_x spectrum (linear)
        IF par.nkx0 / 2 GT 1 THEN BEGIN
          spectra_y_plot, kx, REFORM(xsumspec[*,iky,*,isp,n]), $
            title=title, ytitle='!12<!9!!!6A!9!!!6!U2!N!12>!6', $
            xtitle='!6k!Dx!N', minmax=(*i).minmax, $
            fit=(*i).fit, colorarr=colorarr, linestyle=linestyle, $
            vars=(*i).vars, kmult=(*i).kmult_var, xname='kx', $
            position=[0.15,0.231,0.95,linplot_end], log=0

          plot_legend, colorarr, kx_str + get_var_string((*i).vars,/fancy), $
            x_offset=[0.15,0.0,0.95,0.12], linestyles=linestyle, psym=psym
        ENDIF

        plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
      ENDFOR ; time loop
    ENDFOR ; ky loop

    n = 0
 ;   FOR iky = 0, N_ELEMENTS((*i).kyind) - 1 DO $
 ;     set_output, diag, sp, header=['kx',get_var_string((*i).vars)], $
 ;     dat=[[kx],[REFORM(xsumspec[*,iky,*,isp,n])]], commentline='kyind = '+$
 ;     rm0es((*i).kyind[iky]), append=(iky GT 0)
 ;   IF N_ELEMENTS(ky) EQ 1 THEN FOR ikx = 0, N_ELEMENTS((*i).kxind) - 1 DO $
 ;     set_output, diag, sp, header=['ky',get_var_string((*i).vars)], /append, $
 ;     dat=[ky,REFORM(ysumspec[*,ikx,*,isp,n])], commentline='kxind = '+$
 ;     rm0es((*i).kxind[ikx]) ELSE $
 ;     FOR ikx = 0, N_ELEMENTS((*i).kxind) - 1 DO $
 ;     set_output, diag, sp, header=['ky',get_var_string((*i).vars)], /append, $
 ;     dat=[[ky],[REFORM(ysumspec[*,ikx,*,isp,n])]], commentline='kxind = '+$
 ;     rm0es((*i).kxind[ikx])

    set_output, diag, sp, /reset
  ENDFOR ; species loop

END
