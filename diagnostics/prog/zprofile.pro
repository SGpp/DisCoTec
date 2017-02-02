FUNCTION zprofile_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'Amplitude profiles (z)',$
    help_text : ['Plots the time averaged parallel (z) profiles. '+$
      'Given kxind AND kyind are undefined, the ouput will contain '+$
      'the absolute amplitude, the squared amplitude and variance '+$
      'if spatial averages over kx and ky are possible. Otherwise, '+$
      'only the absolute value is plotted for each given kxind/kyind '+$
      'combination (if kx and ky are undefined and no averaging is '+$
      'possible, the default is kxind=0, kyind=1). '+$
      'Concerning the data output: at the moment only the first chosen '+$
      'variable will be written to the data file.'],$
    ext_vars  : [$
      ['vars','0','index of the fluctuating quantity to be '+$
       'plotted (see variable list); default: all'],$
      ['kxind','0','kx mode index (integer '+$
       'between 0 and nkx0 - 1; default: undefined)'],$
      ['kyind','0','ky mode index (integer '+$
       'between 0 and nky0 - 1; default: undefined)'],$
      ['log','1','switch on/off logarithmic ordinate'],$
      ['calc_av','1','switch to average plot (default: autodetect)']]}

END

;#######################################################################

PRO zprofile_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars
 
  IF N_ELEMENTS(*e.vars) LT 1 THEN $
    IF par.beta GT 0 THEN *e.vars = INDGEN(par.n_moms+par.n_fields) $
    ELSE *e.vars = [0,INDGEN(par.n_moms)+par.n_fields]

  calc_av = (N_ELEMENTS(*e.kxind) LT 1) AND (N_ELEMENTS(*e.kyind) LT 1) AND $
    (par.kx0_ind EQ 0) AND (par.lx NE 0)

  IF KEYWORD_SET(*e.calc_av) THEN calc_av = 1

  calc_av = calc_av AND ((par.nkx0 GT 1) OR (par.nky0 GT 1))

  IF N_ELEMENTS(*e.kxind) LT 1 THEN *e.kxind = 0
  IF N_ELEMENTS(*e.kyind) LT 1 THEN *e.kyind = ((par.ky0_ind EQ 0) AND $
                                                par.nky0 GT 1)
  IF N_ELEMENTS(*e.log) LT 1 THEN *e.log = 0

  i = set_internal_vars(diag,{$
    vars         : *e.vars,$
    kxind   	 : *e.kxind,$
    kyind   	 : *e.kyind,$
    log     	 : *e.log,$
    calc_av 	 : calc_av,$
    z_profile_id : PTR_NEW()})

  IF par.x_local THEN fft_format, kxky=*e.vars ELSE $
    fft_format, sxky=*e.vars

END

;#######################################################################

PRO zprofile_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF par.x_local THEN facnx = 1.0 ELSE facnx = 1.0 / par.nx0

  ; <|A|>, <|A|^2>^(1/2), <|A|^2-<|A|>^2>^(1/2)
  IF (*i).calc_av THEN n_profs = 3 ELSE $
    n_profs = N_ELEMENTS((*i).kxind) * N_ELEMENTS((*i).kyind)
  z_profile = FLTARR(par.nz0,N_ELEMENTS((*i).vars),gui.out.n_spec_sel,$
    n_profs,/NOZERO)
  IF n_profs LT 2 THEN z_profile = $
    REFORM(z_profile,[par.nz0,N_ELEMENTS((*i).vars),gui.out.n_spec_sel,$
    n_profs],/OVERWRITE)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    FOR ivar = 0, N_ELEMENTS((*i).vars) - 1 DO BEGIN
      IF par.x_local THEN datptr = mom[isp,(*i).vars[ivar]].kxky ELSE $
        datptr = mom[isp,(*i).vars[ivar]].sxky

      IF (*i).calc_av THEN BEGIN ; spatial averaged z profile
        ; --- <|A|> ---
        IF par.x_local THEN z_profile[*,ivar,isp,0] = ABS((*datptr)[0,0,*]) $
          ELSE z_profile[*,ivar,isp,0] = TOTAL(ABS((*datptr)[*,0,*]),1) * facnx

        ; --- <|A|^2>^(1/2) ---
        zprof1_sqd = REFORM(TOTAL(ABS((*datptr)[*,0,*])^2,1))
        IF par.nky0 GT 1 THEN zprof1_sqd += 2.0 * TOTAL(TOTAL(ABS($
          (*datptr)[*,1:par.nky0-1,*])^2,1),1)
        zprof1_sqd *= facnx

        ; --- <|A|^2-<|A|>^2>^(1/2) ---
        z_profile[*,ivar,isp,2] = SQRT(zprof1_sqd-z_profile[*,ivar,isp,0]^2)
        z_profile[*,ivar,isp,1] = SQRT(zprof1_sqd)
      ENDIF ELSE BEGIN ; z profile at certain points --- <|A|> ---
        iprof = 0
        FOR ikx = 0, N_ELEMENTS((*i).kxind) - 1 DO BEGIN
          FOR iky = 0, N_ELEMENTS((*i).kyind) - 1 DO BEGIN
            z_profile[*,ivar,isp,iprof] = $
              ABS((*datptr)[(*i).kxind[ikx],(*i).kyind[iky],*])
            iprof += 1
          ENDFOR
        ENDFOR
      ENDELSE
    ENDFOR ; ivar
  ENDFOR ; --- species loop

  (*i).z_profile_id = time_avg((*i).z_profile_id,z_profile,mom_time,$
    fwd=gui.out.res_steps)

END

;######################################################################

PRO zprofile_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  z_profile = time_avg((*i).z_profile_id,/avg,fwd=gui.out.res_steps,tarr=time)

  ; this part ensures that variables have always the same color
  ; (if par.n_fields = const.)
  
  keep_colors = par.n_moms EQ 6

  IF keep_colors THEN BEGIN
    colorarr = INDGEN(par.n_moms+par.n_fields) + 1
    linestyle = INTARR(par.n_moms+par.n_fields)

    FOR var = 6, par.n_moms - 1 DO BEGIN
      colorarr[var+par.n_fields] = colorarr[par.n_fields+(var MOD 6)]
      linestyle[var+par.n_fields] = var / 6
    ENDFOR
   
    colorarr = colorarr[(*i).vars]
    linestyle = linestyle[(*i).vars]
  ENDIF ELSE BEGIN
    colorarr = INDGEN(N_ELEMENTS((*i).vars)) + 1
    linestyle = INTARR(N_ELEMENTS((*i).vars))
  ENDELSE

  xrange = [-!PI*par.n_pol,!PI*par.n_pol]
  fac=4.0
  num=fix(fac)*par.n_pol
  ;xtickname may be at most 60 entries long
  while num GT 59 do begin
     fac=fac/2.0
     num=fix(fac)*par.n_pol
  endwhile
  xtickv = 0.5 * (4.0/fac) * (INDGEN(fix(fac)*par.n_pol+1)-fix(fac)/2.0*par.n_pol)
  xtickname = STRARR(fix(fac)*par.n_pol+1)
  xtickname = rm0es(xtickv) + '!7p!6'
  xtickname[fix(fac)/2.0*par.n_pol] = '0'  
  xtickv *= !PI
  IF par.n_pol GT 1 THEN xtickname[2*INDGEN(par.n_pol+2)+1] = ' '

  ytitle = ['!12<!9!!!6A(z)!9!!!6!N!12>!6',$
    	    '!12<!9!!!6A(z)!9!!!6!E2!N!12>!6!E1/2!N',$
    	    '!12<!9!!!6A(z)-!12<!6A(z)!12>!6!9!!!6!E2!N!12>!6!E1/2!N']
  title = ''
  IF (*i).log THEN ytickformat = 'betterticks' ELSE ytickformat = ''

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]
   
    set_output, diag, sp, /ps, xsize=25, ysize=17

    FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
      IF (*i).calc_av THEN BEGIN ; plot averages
        med0 = MEDIAN(z_profile[*,*,isp,0,n])
        med2 = MEDIAN(z_profile[*,*,isp,2,n])
        j0 = med0 GE 1e-3 * med2 ? 0 : 2

        FOR j = j0, 2 DO BEGIN ; loop over profile
          maxy = MAX(z_profile[*,*,isp,j,n],MIN=miny)
          IF (*i).log AND (miny LE 1e-3) THEN BEGIN
            logminy = 1e-3
            PRINT, 'suppressed values < 1e-3 in zprofile (log plot)'
          ENDIF ELSE logminy = miny

          PLOT, *par.z, z_profile[*,0,isp,j,n], COLOR=1, XTITLE='!6z/q!I0!NR', $
            YTITLE=ytitle[j], TITLE=title, YLOG=(*i).log, /NODATA, /XSTYLE, $
            XRANGE=xrange, XTICKS=N_ELEMENTS(xtickv)-1, XMINOR=2,$
            XTICKV=xtickv, XTICKNAME=xtickname, XTICKLEN=1.0, CHARSIZE=1, $
            POSITION=[0.1,0.25,0.95,0.95], YTICKFORMAT=ytickformat, $
            YRANGE=[logminy,maxy]

          FOR ivar = 0, N_ELEMENTS((*i).vars) - 1 DO $
            OPLOT, (*par.z), z_profile[*,ivar,isp,j,n], COLOR=colorarr[ivar], $
            LINESTYLE=linestyle[ivar] 

          IF (!Y.CRANGE[0] LE 0) AND (!Y.CRANGE[1] GE 0) THEN $
            OPLOT, !X.CRANGE, [0,0], COLOR=1

          plot_legend, colorarr, get_var_string((*i).vars,/fancy), $
            x_offset=[0.1,0.0,0.95,3.0/25.0], per_line=3, linestyles=linestyle

          plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
        ENDFOR

        dtitle = ['<|A|>','<|A|²>^½','<|A²|-|A|²>^½']
        FOR ivar = 0, N_ELEMENTS((*i).vars) - 1 DO $
          set_output, diag, sp, header=['z/qR',[dtitle[*]]], $
          dat=[[(*par.z)],[REFORM(z_profile[*,ivar,isp,*])]], $
          commentline='variable A = '+get_var_string((*i).vars[ivar]), $
          append=(ivar GT 0)
      ENDIF ELSE BEGIN ; --- plot z profile for certain kx/ky modes ---
        j = 0
        FOR ikx = 0, N_ELEMENTS((*i).kxind) - 1 DO BEGIN
          FOR iky = 0, N_ELEMENTS((*i).kyind) - 1 DO BEGIN
            kystr = rm0es((*series.ky)[(*i).kyind[iky]],prec=3)

            IF par.x_local THEN BEGIN
              kxstr = rm0es((*series.kx)[(*i).kxind[ikx],(*i).kyind[iky]],prec=3)
              title = 'k!Ix!N='+kxstr+', k!Iy!N='+kystr
              IF j EQ 0 THEN dtitle = 'kx=' + kxstr + ',ky=' + kystr ELSE $
                dtitle = [dtitle,'kx='+kxstr+',ky='+kystr]      
            ENDIF ELSE BEGIN
              xstr = rm0es((*i).kxind[ikx]*series.dx-series.lx/2,prec=3)
              title = 'x=' + xstr + ', k!Iy!N=' + kystr
              IF j EQ 0 THEN dtitle = 'x=' + xstr + ',ky=' + kystr ELSE $
                dtitle = [dtitle,'x='+xstr+',ky='+kystr]  
            ENDELSE

            maxy = MAX(z_profile[*,*,isp,j,n],MIN=miny)  
            IF (*i).log AND (miny LE 1e-3) THEN logminy = 1e-3 ELSE $
              logminy=miny

            PLOT, (*par.z), z_profile[*,0,isp,j,n], COLOR=1, $
              XTITLE='!6z/q!I0!NR', YTITLE=ytitle[0], TITLE=title, $
              YLOG=(*i).log, /XSTYLE, XRANGE=xrange, XMINOR=2, $
              XTICKS=N_ELEMENTS(xtickv)-1, XTICKV=xtickv, XTICKNAME=xtickname, $
              XTICKLEN=1.0, POSITION=[0.10,0.25,0.95,0.95], $
              YTICKFORMAT=ytickformat, YRANGE=[logminy,maxy], /NODATA

            FOR ivar = 0,N_ELEMENTS((*i).vars) - 1 DO $
              OPLOT, (*par.z), z_profile[*,ivar,isp,j,n], $
              COLOR=colorarr[ivar], LINESTYLE=linestyle[ivar]

            IF (!Y.CRANGE[0] LE 0) AND (!Y.CRANGE[1] GE 0) THEN $
              OPLOT, !X.CRANGE, [0,0], COLOR=1

            plot_legend, colorarr, get_var_string((*i).vars,/fancy), $
              x_offset=[0.10,0.0,0.95,3.0/25.0], per_line=3, $
              linestyles=linestyle

            plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)

            j += 1
          ENDFOR
        ENDFOR

        FOR ivar = 0,N_ELEMENTS((*i).vars) - 1 DO $
          set_output, diag, sp, header=['z/qR',[dtitle[*]]], $
        dat=[[(*par.z)],[REFORM(z_profile[*,ivar,isp,*,n])]], $
        commentline='variable A = '+get_var_string((*i).vars[ivar]), $
        append=(ivar GT 0)
      ENDELSE
    ENDFOR ; --- steps loop

    set_output, diag, sp, /reset
  ENDFOR ; --- species loop

END
