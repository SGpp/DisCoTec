FUNCTION vspfpdiff_info

  RETURN, {$
    type      : 'vsp',$
    title     : 'Fast particle diffusivity',$
    help_text : ['Plots the particle diffusivity (electrostatic '+$
                 'and electromagnetic) to determine its scalings '+$
                 'with energy and pitch angle.'],$
    ext_vars  : [['zind','0','z indices for average; default: -1 (all)'],$
                 ['eta','0','pitch angle array; default: [0.3,0.9] '+$
                  '(corresponding to trapped/passing particles)'],$
                 ['plot_range','0','plot range (x axis min/max, '+$
                  'y axis min/max; default: '+$
                  '[-1,-1,-1,-1] ([E_min,E_max,D_min,D_max])'],$
                 ['posneg','1','separately plot positive/negative '+$
                  'v_par contributions (otherwise, average)']]}

END

;######################################################################

PRO vspfpdiff_init, diag

  COMMON global_vars

  e = (*diag).external_vars
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = -1
  zind = (*e.zind)[0] EQ -1 ? INDGEN(par.nz0) : *e.zind
  IF N_ELEMENTS(*e.plot_range) EQ 2 THEN $
    plot_range = [*e.plot_range,-1,-1] ELSE $
    plot_range = N_ELEMENTS(*e.plot_range) EQ 4 ? $
    *e.plot_range : [-1,-1,-1,-1]
  eta = N_ELEMENTS(*e.eta) LT 1 ? [0.3,0.9] : *e.eta
  posneg = KEYWORD_SET(*e.posneg)

  nz = N_ELEMENTS(zind)
  nv = par.nv0 / 2
  n_eta = N_ELEMENTS(eta)
  B_field = [(*series.geom).Bfield[zind]]

  dv = 2 * par.lv / (par.nv0 - 1.0)
  v_par = (INDGEN(nv) + 0.5) * dv
  GetMuWeightsAndKnots, mu_weights, mu

  ; interpolation from vpar-mu grid to eta-E grid
  ; eta = vpar / vtot = sqrt(1 - vperp^2/vtot^2)
  ;     = vpar / sqrt(vpar^2 + mu B)
  ; E = vpar^2 + mu B = mu B / (1 - eta^2)
  ; vpar = eta * sqrt(mu B / (1 - eta^2)) = sqrt(mu B / (eta^-2 - 1))
  ; mu B = vpar^2 (eta^-2 - 1)
  ; E = vpar^2 / eta^2

  mu_mod = pf_arr([nz,nv,n_eta])
  mu_ind_lo = pf_arr([nz,nv,n_eta])
  mu_ind_hi = pf_arr([nz,nv,n_eta])
  dmu_lo = pf_arr([nz,nv,n_eta])
  dmu_hi = pf_arr([nz,nv,n_eta])

  FOR z = 0, nz - 1 DO mu_mod[z,0,0] = $
    REFORM(REBIN(v_par^2/B_field[z],[nv,n_eta])*$
    REBIN(REFORM(eta^(-2)-1,[1,n_eta]),[nv,n_eta]),[1,nv,n_eta])

  mu_ind_lo[0,0,0] = (VALUE_LOCATE(mu,mu_mod[*,*,*]) > 0) < (nv - 2)
  mu_ind_hi[0,0,0] = mu_ind_lo[*,*,*] + 1
  dmu_inv = 1.0 / (mu[mu_ind_hi] - mu[mu_ind_lo])
  dmu_lo[0,0,0] = (mu_mod - mu[mu_ind_lo]) * dmu_inv
  dmu_hi[0,0,0] = (mu[mu_ind_hi] - TEMPORARY(mu_mod)) * dmu_inv

  ; inverse Maxwellian for normalization, includes 0.5 for nv0->nv0/2
  f_maxwell_inv = pf_arr([nz,nv,par.nw0])
  vpar_sqd = REBIN(REFORM(v_par^2,[1,nv]),[nz,nv])
  FOR w = 0, par.nw0 - 1 DO f_maxwell_inv[0,0,w] = $
    EXP(vpar_sqd+REBIN(mu[w]*B_field,[nz,nv]))
  f_maxwell_inv = REBIN((0.5*(1+posneg)*!PI^(1.5))*$
    TEMPORARY(f_maxwell_inv),[nz,nv,par.nw0,2])

  ; add factor of 1/d3v to f_maxwell_inv (use exact dv from GENE)
  dv = FLTARR(nv) + par.lv / (nv - 1.0)
  dv[nv-4:nv-1] *= [49,43,59,17] / 48.0
  d3v_inv = REBIN(REFORM(1.0/mu_weights,[1,1,par.nw0]),[nz,nv,par.nw0]) * $
    REBIN(REFORM(1.0/dv,[1,nv,1]),[nz,nv,par.nw0]) * $
    REBIN(1.0/(!PI*(*series.geom).Bfield[zind]),[nz,nv,par.nw0])
  f_maxwell_inv = TEMPORARY(f_maxwell_inv) * d3v_inv

  i = set_internal_vars(diag,{$
    zind          : zind,$
    nz            : nz,$
    plot_range    : plot_range,$
    posneg        : posneg,$
    eta           : eta,$
    n_eta         : n_eta,$
    v_par         : v_par,$
    mu            : mu,$
    mu_ind_lo     : mu_ind_lo,$
    mu_ind_hi     : mu_ind_hi,$
    dmu_lo        : dmu_lo,$
    dmu_hi        : dmu_hi,$
    f_maxwell_inv : f_maxwell_inv,$
    vdat_id       : PTR_NEW(),$
    vneg_id       : PTR_NEW()})

END

;######################################################################

PRO vspfpdiff_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nv = par.nv0 / 2

  IF (*i).posneg THEN BEGIN
    (*i).vdat_id = time_avg((*i).vdat_id,$
      (*vsp)[(*i).zind,nv+INDGEN(nv),*,*,0:1],$
      vsp_time,fwd=gui.out.res_steps)
    (*i).vneg_id = time_avg((*i).vneg_id,$
      (*vsp)[(*i).zind,nv-1-INDGEN(nv),*,*,0:1],$
      vsp_time,fwd=gui.out.res_steps)
  ENDIF ELSE (*i).vdat_id = time_avg((*i).vdat_id,$
    (*vsp)[(*i).zind,nv-1-INDGEN(nv),*,*,0:1]+$
    (*vsp)[(*i).zind,nv+INDGEN(nv),*,*,0:1],$
    vsp_time,fwd=gui.out.res_steps)

END

;######################################################################

PRO vspfpdiff_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  ymarg = !Y.MARGIN
  !Y.MARGIN = [ymarg[0],6]

  nv = par.nv0 / 2
  GetMuWeightsAndKnots, mu_weights, mu

  diffusivity_es = pf_arr([(*i).nz,nv,par.nw0])
  diffusivity_em = pf_arr([(*i).nz,nv,par.nw0])
  diff_E_z = pf_arr([(*i).nz,nv,(*i).n_eta,2])
  diff_E = pf_arr([nv,(*i).n_eta,2])

  vsp_data = REFORM(time_avg((*i).vdat_id,/avg,fwd=gui.out.res_steps,tarr=time),$
    [(*i).nz,nv,par.nw0,gui.out.n_spec_sel,2,1+gui.out.res_steps*(series.step_count-1)])
  IF (*i).posneg THEN BEGIN
    vneg_data = REFORM(time_avg((*i).vneg_id,/avg,fwd=gui.out.res_steps,tarr=time),$
      [(*i).nz,nv,par.nw0,gui.out.n_spec_sel,2,1+gui.out.res_steps*(series.step_count-1)])
    diff_E_neg = pf_arr([nv,(*i).n_eta,2])
  ENDIF

  gxx_inv = REBIN([1.0/(*series.geom).gxx[(*i).zind]],$
    [(*i).nz,nv,(*i).n_eta])

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]
    omn_inv = 1.0 / spec[sp].omn

    set_output, diag, sp, /ps

    FOR n = 0, (series.step_count - 1) * gui.out.res_steps DO BEGIN
      diffusivity_es[0,0,0] = omn_inv * (*i).f_maxwell_inv * $
        REFORM(vsp_data[*,*,*,isp,0,n],[(*i).nz,nv,par.nw0])
      diffusivity_em[0,0,0] = omn_inv * (*i).f_maxwell_inv * $
        REFORM(vsp_data[*,*,*,isp,1,n],[(*i).nz,nv,par.nw0])

      ; coordinate transformation
      temp_z = REBIN(INDGEN((*i).nz),[(*i).nz,nv,(*i).n_eta])
      temp_v = REBIN(REFORM(INDGEN(nv),[1,nv]),$
        [(*i).nz,nv,(*i).n_eta])
      diff_temp0 = $
        (*i).dmu_hi * diffusivity_es[temp_z,temp_v,(*i).mu_ind_lo] + $
        (*i).dmu_lo * diffusivity_es[temp_z,temp_v,(*i).mu_ind_hi]
      diff_temp1 = $
        (*i).dmu_hi * diffusivity_em[temp_z,temp_v,(*i).mu_ind_lo] + $
        (*i).dmu_lo * diffusivity_em[temp_z,temp_v,(*i).mu_ind_hi]

      diff_E_z[0,0,0,0] = diff_temp0 * gxx_inv
      diff_E_z[0,0,0,1] = diff_temp1 * gxx_inv

      diff_E[0,0,0] = REFORM(TOTAL(diff_E_z*$
        REBIN([(*series.geom).jac_norm[(*i).zind]/(*i).nz],$
        [(*i).nz,nv,(*i).n_eta,2]),1),[nv,(*i).n_eta,2])

      IF (*i).posneg THEN BEGIN
        ; normalization to (rho^2 v_th / R) R / rho
        diffusivity_es[0,0,0] = omn_inv * (*i).f_maxwell_inv * $
          REFORM(vneg_data[*,*,*,isp,0,n],[(*i).nz,nv,par.nw0])
        diffusivity_em[0,0,0] = omn_inv * (*i).f_maxwell_inv * $
          REFORM(vneg_data[*,*,*,isp,1,n],[(*i).nz,nv,par.nw0])

        ; coordinate transformation
        diff_temp0 = $
          (*i).dmu_hi * diffusivity_es[temp_z,temp_v,(*i).mu_ind_lo] + $
          (*i).dmu_lo * diffusivity_es[temp_z,temp_v,(*i).mu_ind_hi]
        diff_temp1 = $
          (*i).dmu_hi * diffusivity_em[temp_z,temp_v,(*i).mu_ind_lo] + $
          (*i).dmu_lo * diffusivity_em[temp_z,temp_v,(*i).mu_ind_hi]

        diff_E_z[0,0,0,0] = diff_temp0 * gxx_inv
        diff_E_z[0,0,0,1] = diff_temp1 * gxx_inv

        diff_E_neg[0,0,0] = REFORM(TOTAL(diff_E_z*$
          REBIN([(*series.geom).jac_norm[(*i).zind]/(*i).nz],$
          [(*i).nz,nv,(*i).n_eta,2]),1),[nv,(*i).n_eta,2])
      ENDIF

      FOR p = 0, (*i).n_eta - 1 DO BEGIN
        title = '!6D(!7g!6=' + rm0es((*i).eta[p]) + ')'

        E_tot = (*i).v_par^2 * (spec[sp].temp / (*i).eta[p]^2)

        xrange = [E_tot[0],E_tot[nv-1]]
        yrange = [MIN(ABS(diff_E[*,p,*])),MAX(ABS(diff_E[*,p,*]))]
        IF (*i).plot_range[0] NE -1 THEN xrange[0] = (*i).plot_range[0]
        IF (*i).plot_range[1] NE -1 THEN xrange[1] = (*i).plot_range[1]
        IF (*i).plot_range[2] NE -1 THEN yrange[0] = (*i).plot_range[2]
        IF (*i).plot_range[3] NE -1 THEN yrange[1] = (*i).plot_range[3]

        PLOT, E_tot, ABS(diff_E[*,p,0]), COLOR=1, TITLE=title, $
          XTITLE='!6E / T!Dref!N', YTITLE='!6D(E)', $
          /XLOG, /YLOG, /XSTYLE, /YSTYLE, XRANGE=xrange, YRANGE=yrange
        OPLOT, E_tot, ABS(diff_E[*,p,1]), COLOR=1, LINE=2

        IF (*i).posneg THEN BEGIN
          OPLOT, E_tot, ABS(diff_E_neg[*,p,0]), COLOR=2
          OPLOT, E_tot, ABS(diff_E_neg[*,p,1]), COLOR=2, LINE=2
        ENDIF

        ; mark points resulting from negative values
        neg_range_es = WHERE(diff_E[*,p,0] LT 0.0,count_neg_es)
        neg_range_em = WHERE(diff_E[*,p,1] LT 0.0,count_neg_em)
        IF count_neg_es EQ 1 THEN neg_range_es = [1,1] * neg_range_es
        IF count_neg_em EQ 1 THEN neg_range_em = [1,1] * neg_range_em
        IF count_neg_es GT 0 THEN OPLOT, E_tot[neg_range_es], $
          - diff_E[neg_range_es,p,0], COLOR=8, PSYM=4
        IF count_neg_em GT 0 THEN OPLOT, E_tot[neg_range_em], $
          - diff_E[neg_range_em,p,1], COLOR=8, PSYM=4
        IF (*i).posneg THEN BEGIN
          neg_range_es = WHERE(diff_E_neg[*,p,0] LT 0.0,count_neg_es)
          neg_range_em = WHERE(diff_E_neg[*,p,1] LT 0.0,count_neg_em)
          IF count_neg_es EQ 1 THEN neg_range_es = [1,1] * neg_range_es
          IF count_neg_em EQ 1 THEN neg_range_em = [1,1] * neg_range_em
          IF count_neg_es GT 0 THEN OPLOT, E_tot[neg_range_es], $
            - diff_E_neg[neg_range_es,p,0], COLOR=8, PSYM=4
          IF count_neg_em GT 0 THEN OPLOT, E_tot[neg_range_em], $
            - diff_E_neg[neg_range_em,p,1], COLOR=8, PSYM=4
        ENDIF

        ; lines of relevant slopes
        plotmax = MAX(ABS(diff_E[*,p,*]))
        E_norm_inv = E_tot[0] / E_tot
        OPLOT, E_tot, plotmax * E_norm_inv, COLOR=2, LINE=1
        OPLOT, E_tot, plotmax * (E_norm_inv)^(1.5), COLOR=4, LINE=1
        OPLOT, E_tot, plotmax * SQRT(E_norm_inv), COLOR=6, LINE=1

        ; legend
        IF count_neg_es + count_neg_em GT 0 THEN BEGIN
          PLOT, [1,1] * 0.075, [1,1] * 0.23, COLOR=8, PSYM=4, $
            /NOERASE, XSTYLE=5, YSTYLE=5, XRANGE=[0,1], YRANGE=[0,1]
          XYOUTS, 0.15, 0.22, '!6D(E) < 0', COLOR=8
        ENDIF
        PLOT, [1,2] * 0.05, [1,1] * 0.19, COLOR=1, $
          /NOERASE, XSTYLE=5, YSTYLE=5, XRANGE=[0,1], YRANGE=[0,1]
        XYOUTS, 0.15, 0.18, '!6D!S!D' + spec[sp].name + $
          '!N!R!Ues!N(E' + ((*i).posneg ? ',v!D!9#!6!N!U  !N0)' : $
          ')'), COLOR=1
        OPLOT, [1,2] * 0.05, [1,1] * 0.15, COLOR=1, LINE=2
        XYOUTS, 0.15, 0.14, '!6D!S!D' + spec[sp].name + $
          '!N!R!Uem!N(E' + ((*i).posneg ? ',v!D!9#!6!N!U  !N0)' : $
          ')'), COLOR=1
        IF (*i).posneg THEN BEGIN
          XYOUTS, 0.315, 0.175, '!6!U>!N', COLOR=1
          XYOUTS, 0.33, 0.135, '!6!U>!N', COLOR=1
          XYOUTS, 0.315, 0.185, '!6!D<!N', COLOR=2
          XYOUTS, 0.33, 0.145, '!6!D<!N', COLOR=2
        ENDIF
        OPLOT, [1,2] * 0.05, [1,1] * 0.11, COLOR=2, LINE=1
        XYOUTS, 0.15, 0.10, '!6E!U-1!N', COLOR=2
        OPLOT, [1,2] * 0.05, [1,1] * 0.08, COLOR=4, LINE=1
        XYOUTS, 0.15, 0.07, '!6E!U-3/2!N', COLOR=4
        OPLOT, [1,2] * 0.05, [1,1] * 0.05, COLOR=6, LINE=1
        XYOUTS, 0.15, 0.04, '!6E!U-1/2!N', COLOR=6

        plot_info_str, diag
      ENDFOR
    ENDFOR

    set_output, diag, sp, /reset
  ENDFOR

  !Y.MARGIN = ymarg

END
