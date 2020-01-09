FUNCTION energyspect_info

  RETURN, {$
    type      : 'energy',$
    title     : 'Energy spectra',$
    help_text : ['Plots k_x and k_y spectra for the free energy, '+$
                 'energy drive, dissipation, etc.'],$
    ext_vars  : [['vars','0','energy quantities to be analyzed; '+$
                  'default: all'],$
                 ['kx_range','0','plot range, [min_ind,max_ind], '+$
                  'of kx indices; default: all '+$
                  '(without zero mode for log-log plot)'],$
                 ['ky_range','0','plot range, [min_ind,max_ind], '+$
                  'of ky indices; default: all '+$
                  '(without zero mode for log-log plot)'],$
                 ['kperp_range','0','[kperp_min,kperp_max]; if '+$
                  'specified, prints the relative contribution '+$
                  'of the range; unlike kxy_range, assumes '+$
                  'physical values; default: off'],$
                 ['add_transf','1','add spectra of nonlinear '+$
                  'transfer rate']]}

END

;######################################################################

PRO energyspect_init, diag

  COMMON global_vars

  e = (*diag).external_vars
  vars = N_ELEMENTS(*e.vars) GE 1 ? *e.vars : INDGEN(6)
  kx_range = N_ELEMENTS(*e.kx_range) EQ 2 ? *e.kx_range : $
    [0,par.nx0-par.nx0/2-1]
  ky_range = N_ELEMENTS(*e.ky_range) EQ 2 ? *e.ky_range : [0,par.nky0-1]
  en_ratios = N_ELEMENTS(*e.kperp_range) EQ 2
  add_transf = KEYWORD_SET(*e.add_transf)

  IF add_transf THEN BEGIN
    dummy = WHERE(vars EQ 0,vars_have_0)
    dummy = WHERE(vars EQ 5,vars_have_5)

    IF vars_have_0 LT 1 THEN vars = [vars,0]
    IF vars_have_5 LT 1 THEN vars = [vars,5]
  ENDIF

  IF en_ratios THEN BEGIN
    kperp_range = *e.kperp_range

    k_perp2 = pf_arr([par.nx0-par.nx0/2,par.nky0,par.nz0])
    kx_rebin = REBIN((*series.kx)[0:par.nx0-par.nx0/2-1,0],$
      [par.nx0-par.nx0/2,par.nky0])
    ky_rebin = REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0-par.nx0/2,par.nky0])

    FOR z = 0, par.nz0 - 1 DO k_perp2[0,0,z] = $
      (*series.geom).gxx[z] * kx_rebin^2 + $
      2.0*(*series.geom).gxy[z] * kx_rebin * ky_rebin + $
      (*series.geom).gyy[z] * ky_rebin^2

    kperp_inds = WHERE((k_perp2 GE kperp_range[0]^2) AND $
      (k_perp2 LE kperp_range[1]^2),count)
    IF count LT 1 THEN en_ratios = 0
  ENDIF ELSE kperp_inds = 0

  i = set_internal_vars(diag,{$
    vars          : vars,$
    n_vars        : N_ELEMENTS(vars),$
    kx_range      : kx_range,$
    ky_range      : ky_range,$
    add_transf    : add_transf,$
    en_ratios     : en_ratios,$
    kperp_inds    : kperp_inds,$
    data_kx_id    : PTR_NEW(),$
    data_ky_id    : PTR_NEW(),$
    data_kperp_id : PTR_NEW()})

  ;fft_format_en, kxky=*e.vars

END

;######################################################################

PRO energyspect_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  data_kx = pf_arr([par.nkx0-par.nkx0/2,(*i).n_vars])
  data_ky = pf_arr([par.nky0,(*i).n_vars])
  data_kperp = pf_arr((*i).n_vars)

  FOR v = 0, (*i).n_vars - 1 DO BEGIN
    data_zavg = TOTAL(REFORM((*energy[v].kxky)[*,*,*]*REBIN(REFORM($
      (*series.geom).jac_norm/par.nz0,[1,1,par.nz0]),$
      [par.nkx0,par.nky0,par.nz0]),[par.nkx0,par.nky0,par.nz0]),3)
    ; add positive and negative kx
    temp_kxky = FLTARR(par.nkx0-par.nkx0/2,par.nky0)
    temp_kxky[1,0] = data_zavg[par.nkx0-1-INDGEN(par.nkx0-par.nkx0/2-1),*]
    temp_kxky += data_zavg[0:par.nkx0-par.nkx0/2-1,*]

    data_kx[0,v] = 2 * TOTAL(temp_kxky,2) - temp_kxky[*,0]
    data_ky[0,v] = TOTAL(temp_kxky,1)

    IF (*i).en_ratios THEN BEGIN
      ; add positive and negative kx
      temp_0pos = FLTARR(par.nkx0-par.nkx0/2,par.nky0,par.nz0)
      temp_0pos[1,0,0] = $
        (*energy[v].kxky)[par.nkx0-1-INDGEN(par.nkx0-par.nkx0/2-1),*,*]
      temp_0pos += (*energy[v].kxky)[0:par.nkx0-par.nkx0/2-1,*,*]

      ; apply Jacobian even though not all z points may be included!
      temp_0pos *= REBIN(REFORM((*series.geom).jac_norm/par.nz0,$
        [1,1,par.nz0]),[par.nkx0-par.nkx0/2,par.nky0,par.nz0])
      temp_0pos[0,1,0] = 2 * temp_0pos[*,1:*,*]

      data_kperp[v] = TOTAL(temp_0pos[(*i).kperp_inds])
    ENDIF
  ENDFOR

  (*i).data_kx_id = time_avg((*i).data_kx_id,data_kx,energy_time,$
    fwd=gui.out.res_steps)
  (*i).data_ky_id = time_avg((*i).data_ky_id,data_ky,energy_time,$
    fwd=gui.out.res_steps)
  IF (*i).en_ratios THEN (*i).data_kperp_id = time_avg($
    (*i).data_kperp_id,data_kperp,energy_time,fwd=gui.out.res_steps)

END

;######################################################################

PRO energyspect_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  show_power_law_fits = 0 ; fits straight lines to log-log data

  kxinds_0pos = INDGEN((*i).kx_range[1]-(*i).kx_range[0]+1) + $
    (*i).kx_range[0]
  kxind_neg_min = par.nkx0 - ((*i).kx_range[0] > 1)
  kxind_neg_max = par.nkx0 - (*i).kx_range[1]
  kxinds_neg = INDGEN(kxind_neg_min-kxind_neg_max+1) + kxind_neg_max

  kyinds = INDGEN((*i).ky_range[1]-(*i).ky_range[0]+1) + (*i).ky_range[0]

  n_kxinds = N_ELEMENTS(kxinds_0pos)
  n_kyinds = N_ELEMENTS(kyinds)

  x_has_0 = kxinds_0pos[0] EQ 0
  y_has_0 = kyinds[0] EQ 0

  data_kx = time_avg((*i).data_kx_id,/avg,fwd=gui.out.res_steps,tarr=tarr)
  data_ky = time_avg((*i).data_ky_id,/avg,fwd=gui.out.res_steps)
  IF (*i).en_ratios THEN data_kperp = $
    time_avg((*i).data_kperp_id,/avg,fwd=gui.out.res_steps)

  en_var_string = '!6' + ['E!Dtot!N',$
    'dE/dt!9!!!6!Dnoncons!N (w/ z,v poisson)',$
    'Q (= dE/dt!9!!!6!Ddrive!N)','C (= dE/dt!9!!!6!Dcoll!N)',$
    'D (hyp_z,v,x,y) (w/o dissipative BC)',$
    'dE/dt!9!!!6!DNL!N','1/E dE/dt!9!!!6!DNL!N']
  IF par.arakawa_zv AND NOT par.arakawa_cons_bc THEN $
    en_var_string[4] += '(w/o dissipative BC)'

  kx_axis = (*series.kx)[kxinds_0pos,0]
  ky_axis = (*series.ky)[kyinds]

  set_output, diag, /ps, multi=[0,1,2]

  delivered_neg_warning = 0

  FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
    FOR v = 0, (*i).n_vars - 1 DO BEGIN ; plot y spectra
      plot_data = data_ky[kyinds,v,n]

      PLOT, ky_axis, plot_data, COLOR=1, PSYM=-7, /XSTYLE, $
        XTITLE='!6k!Dy!N'+get_var_string(/rhostr), $
        TITLE=en_var_string[(*i).vars[v]]

      PLOT, ky_axis[y_has_0:*], ABS(plot_data[y_has_0:*]), COLOR=1, PSYM=-7, $
        XTITLE='!6k!Dy!N'+get_var_string(/rhostr), $
        TITLE=en_var_string[(*i).vars[v]], /XSTYLE, /XLOG, /YLOG

      neg_inds = WHERE(plot_data LT 0.0,neg_count)
      IF neg_count GE 1 THEN BEGIN
        OPLOT, ky_axis[neg_inds], -plot_data[neg_inds], COLOR=2, PSYM=7
        IF NOT delivered_neg_warning THEN PRINT, (*diag).name + $
          ': negative values indicated by red symbols in log-log plot'
        delivered_neg_warning = 1
      ENDIF

      plot_info_str, diag, time=(gui.out.res_steps ? tarr[n] : undef)

      IF show_power_law_fits THEN BEGIN
        fit = LINFIT(ALOG10(ky_axis[y_has_0:*]),$
          ALOG10(ABS(plot_data[y_has_0:*])))
        OPLOT, ky_axis[[y_has_0,n_kyinds-1]], $
          10^(fit[0]+fit[1]*ALOG10(ky_axis[[y_has_0,n_kyinds-1]])), $
          COLOR=4, LINE=1
        PRINT, ' ' + en_var_string[(*i).vars[v]] + ' ky slope:' + rm0es(fit[1])
      ENDIF
    ENDFOR

    IF (*i).add_transf THEN BEGIN ; nonlinear transfer rate y spectrum
      ind0 = WHERE((*i).vars EQ 0)
      ind5 = WHERE((*i).vars EQ 5)

      plot_data = (data_ky[kyinds,*,n])[*,ind5] / (data_ky[kyinds,*,n])[*,ind0]
      PLOT, ky_axis, plot_data, COLOR=1, PSYM=-7, /XSTYLE, $
        XTITLE='!6k!Dy!N'+get_var_string(/rhostr), TITLE=en_var_string[6]

      PLOT, ky_axis[y_has_0:*], ABS(plot_data[y_has_0:*]), COLOR=1, PSYM=-7, $
        XTITLE='!6k!Dy!N'+get_var_string(/rhostr), TITLE=en_var_string[6], $
        /XSTYLE, /XLOG, /YLOG
      neg_inds = WHERE(plot_data LT 0.0,neg_count)
      IF neg_count GE 1 THEN BEGIN
        OPLOT, ky_axis[neg_inds], -plot_data[neg_inds], COLOR=2, PSYM=7
        IF NOT delivered_neg_warning THEN PRINT, (*diag).name + $
          ': negative values indicated by red symbols in log-log plot'
        delivered_neg_warning = 1
      ENDIF

      plot_info_str, diag, time=(gui.out.res_steps ? tarr[n] : undef)
    ENDIF

    FOR v = 0, (*i).n_vars - 1 DO BEGIN ; plot x spectra
      plot_data = data_kx[kxinds_0pos,v,n]

      PLOT, kx_axis, plot_data, COLOR=1, PSYM=-7, /XSTYLE, $
        XTITLE='!6k!Dx!N'+get_var_string(/rhostr), $
        TITLE=en_var_string[(*i).vars[v]]

      PLOT, kx_axis[x_has_0:*], ABS(plot_data[x_has_0:*]), COLOR=1, PSYM=-7, $
        XTITLE='!6k!Dx!N'+get_var_string(/rhostr), $
        TITLE=en_var_string[(*i).vars[v]], /XSTYLE, /XLOG, /YLOG

      neg_inds = WHERE(plot_data LT 0.0,neg_count)
      IF neg_count GE 1 THEN BEGIN
        OPLOT, kx_axis[neg_inds], -plot_data[neg_inds], COLOR=2, PSYM=7
        IF NOT delivered_neg_warning THEN PRINT, (*diag).name + $
          ': negative values indicated by red symbols in log-log plot'
        delivered_neg_warning = 1
      ENDIF

      plot_info_str, diag, time=(gui.out.res_steps ? tarr[n] : undef)

      IF show_power_law_fits THEN BEGIN
        fit = LINFIT(ALOG10(kx_axis[x_has_0:*]),$
          ALOG10(ABS(plot_data[x_has_0:*])))
        OPLOT, kx_axis[[x_has_0,n_kxinds-1]], $
          10^(fit[0]+fit[1]*ALOG10(kx_axis[[x_has_0,n_kxinds-1]])), $
          COLOR=4, LINE=1
        PRINT, ' ' + en_var_string[(*i).vars[v]] + ' kx slope:' + rm0es(fit[1])
      ENDIF
    ENDFOR

    IF (*i).add_transf THEN BEGIN ; nonlinear transfer rate y spectrum
      plot_data = (data_kx[kxinds_0pos,*,n])[*,ind5] / $
        (data_ky[kxinds_0pos,*,n])[*,ind0]
      PLOT, kx_axis, plot_data, COLOR=1, PSYM=-7, /XSTYLE, $
        XTITLE='!6k!Dx!N'+get_var_string(/rhostr), TITLE=en_var_string[6]

      PLOT, kx_axis[x_has_0:*], ABS(plot_data[x_has_0:*]), COLOR=1, PSYM=-7, $
        XTITLE='!6k!Dx!N'+get_var_string(/rhostr), TITLE=en_var_string[6], $
        /XSTYLE, /XLOG, /YLOG
      neg_inds = WHERE(plot_data LT 0.0,neg_count)
      IF neg_count GE 1 THEN BEGIN
        OPLOT, kx_axis[neg_inds], -plot_data[neg_inds], COLOR=2, PSYM=7
        IF NOT delivered_neg_warning THEN PRINT, (*diag).name + $
          ': negative values indicated by red symbols in log-log plot'
        delivered_neg_warning = 1
      ENDIF

      plot_info_str, diag, time=(gui.out.res_steps ? tarr[n] : undef)
    ENDIF

    IF (*i).en_ratios THEN BEGIN
      FOR v = 0, (*i).n_vars - 1 DO BEGIN
        en_in_range = data_kperp[v,n]
        en_total = TOTAL(data_kx[*,v,n])
        PRINT, en_var_string[(*i).vars[v]] + '_(range,tot,ratio) = ' + $
          rm0es(en_in_range) + ',' + $
          rm0es(en_total) + ',' + rm0es(en_in_range/en_total)
      ENDFOR
    ENDIF
  ENDFOR

  IF NOT gui.out.res_steps THEN BEGIN
    temp_kx = pf_arr([par.nkx0-par.nkx0/2,(*i).n_vars+1])
    temp_kx[0,0] = kx_axis
    temp_kx[0,1] = data_kx
    temp_ky = pf_arr([par.nky0,(*i).n_vars+1])
    temp_ky[0,0] = ky_axis
    temp_ky[0,1] = data_ky

    set_output, diag, dat=temp_ky, $
      commentlines='ky, energy vars: '+STRJOIN(en_var_string[(*i).vars]+', ')
    set_output, diag, dat=temp_kx, /append, $
      commentlines='kx, energy vars: '+STRJOIN(en_var_string[(*i).vars]+', ')
  ENDIF ELSE PRINT, (*diag).name + $
    ': no data file output for resolved time steps'

  set_output, diag, /reset

END
