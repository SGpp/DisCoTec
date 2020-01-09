FUNCTION itmhdf5_info

  RETURN, {$
    type      : 'mom',$
    title     : 'ITM HDF5 output',$
    help_text : ['Creates an HDF5 output file according to the '+$
                 'ITM specifications.'],$
    ext_vars  : ''}

END

;######################################################################

PRO itmhdf5_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  found_e = 0
  found_i = 0
  FOR s = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[s]

    IF spec[sp].charge EQ -1 THEN BEGIN
      found_e = 1
      sp_e = sp
      s_e = s
    ENDIF
    IF spec[sp].charge GT 0 THEN BEGIN
      found_i = 1
      sp_i = sp
      s_i = s
    ENDIF
  ENDFOR
  IF found_e + found_i LT 2 THEN BEGIN
    PRINT, 'no ions and/or electrons found'
    RETURN
  ENDIF

  read_geometry

  i = set_internal_vars(diag,{$
    sp_e                   : sp_e,$
    sp_i                   : sp_i,$
    s_e                    : s_e,$
    s_i                    : s_i,$
    time_id                : PTR_NEW(),$
    energy_exb_id          : PTR_NEW(),$
    energy_mag_id          : PTR_NEW(),$
    energy_electron_th_id  : PTR_NEW(),$
    energy_electron_par_id : PTR_NEW(),$
    energy_ion_par_id      : PTR_NEW(),$
    energy_ion_th_id       : PTR_NEW(),$
    energy_id              : PTR_NEW(),$
    flux_particle_id       : PTR_NEW(),$
    flux_heat_el_id        : PTR_NEW(),$
    flux_ion_id            : PTR_NEW(),$
    flux_heat_ion_id       : PTR_NEW(),$
    flux_mag_particle_id   : PTR_NEW(),$
    flux_mag_ion_id        : PTR_NEW(),$
    flux_mag_heat_el_id    : PTR_NEW(),$
    flux_mag_heat_ion_id   : PTR_NEW(),$
    zonal_ne_id            : PTR_NEW(),$
    zonal_ni_id            : PTR_NEW(),$
    zonal_Jpl_id           : PTR_NEW(),$
    zonal_phi_id           : PTR_NEW(),$
    zonal_Er_id            : PTR_NEW(),$
    zonal_vor_id           : PTR_NEW(),$
    zonal_Apl_id           : PTR_NEW(),$
    zonal_te_id            : PTR_NEW(),$
    zonal_ti_id            : PTR_NEW(),$
    zonal_ui_id            : PTR_NEW(),$
    spec_ne_id             : PTR_NEW(),$
    spec_te_id             : PTR_NEW(),$
    spec_ti_id             : PTR_NEW(),$
    spec_phi_id            : PTR_NEW(),$
    spec_vor_id            : PTR_NEW(),$
    spec_Jpl_id            : PTR_NEW(),$
    spec_Fe_id             : PTR_NEW(),$
    spec_Qe_id             : PTR_NEW(),$
    spec_Qi_id             : PTR_NEW(),$
    spec_Me_id             : PTR_NEW(),$
    spec_Mi_id             : PTR_NEW(),$
    spec_B_id              : PTR_NEW(),$
    env_ne_id              : PTR_NEW(),$
    env_he_id              : PTR_NEW(),$
    env_te_id              : PTR_NEW(),$
    env_Fe_id              : PTR_NEW(),$
    env_Qe_id              : PTR_NEW(),$
    env_Me_id              : PTR_NEW(),$
    env_ni_id              : PTR_NEW(),$
    env_ti_id              : PTR_NEW(),$
    env_ui_id              : PTR_NEW(),$
    env_Qi_id              : PTR_NEW(),$
    env_Mi_id              : PTR_NEW(),$
    env_vor_id             : PTR_NEW(),$
    env_Jpl_id             : PTR_NEW(),$
    env_phi_id             : PTR_NEW(),$
    axisym_phi_id          : PTR_NEW(),$
    axisym_Apl_id          : PTR_NEW(),$
    axisym_ne_id           : PTR_NEW(),$
    axisym_te_id           : PTR_NEW(),$
    axisym_ni_id           : PTR_NEW(),$
    axisym_ti_id           : PTR_NEW(),$
    axisym_ui_id           : PTR_NEW(),$
    axisym_vor_id          : PTR_NEW(),$
    axisym_Jpl_id          : PTR_NEW(),$
    phi_id                 : PTR_NEW(),$
    ne_id                  : PTR_NEW(),$
    Jpl_id                 : PTR_NEW(),$
    vor_id                 : PTR_NEW()})

  nf = par.n_fields
  fft_format, $
    kxky=(nf GT 1 ? [0,1,nf,nf+1,nf+2,nf+3,nf+4,nf+5] : [0,nf,nf+1,nf+2]), $
    sxky=[0,1,nf,nf+1,nf+2,nf+5], sxsy=[0,1,nf,nf+1,nf+2,nf+5]

END

;######################################################################

FUNCTION avg, data_in, all=all, zon=zon, kys=kys, env=env, axi=axi
; averages certain dimensions; input formats: all=kxky, zon=sxky,
; kys=kxky, env=kxky, axi=sxsy

  COMMON global_vars

  IF KEYWORD_SET(all) THEN RETURN, TOTAL(2.0*data_in*REBIN(REFORM($
    (*series.geom).jac_norm/par.nz0,[1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0]))

  IF KEYWORD_SET(zon) THEN RETURN, TOTAL(REFORM(data_in[*,0,*])*REBIN(REFORM($
    (*series.geom).jac_norm/par.nz0,[1,par.nz0]),[par.nx0,par.nz0]),2)

  IF KEYWORD_SET(kys) THEN RETURN, TOTAL(TOTAL(data_in,1)*REBIN(REFORM($
    (*series.geom).jac_norm/par.nz0,[1,par.nz0]),[par.nky0,par.nz0]),2)

  IF KEYWORD_SET(env) THEN BEGIN
    xsum = TOTAL(data_in,1)
    RETURN, 2.0 * TOTAL(xsum,1) - REFORM(xsum[0,*])
  ENDIF

  IF KEYWORD_SET(axi) THEN RETURN, data_in[*,*,par.nz0/2]

END

;######################################################################

PRO itmhdf5_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  s_e = (*i).s_e
  sp_e = (*i).sp_e
  s_i = (*i).s_i
  sp_i = (*i).sp_i
  nf = par.n_fields

  par_mref = par.mref EQ 0.0 ? 1.0 / spec[sp_e].mass * 9.1e-31 : par.mref
  series_mref = series.mref EQ 0.0 ? 1.0 : series.mref

  n_ref = par.nref / series.nref
  L_ref = par.Lref / series.Lref
  m_ref = par_mref / series_mref
  T_ref = par.Tref / series.Tref
  Q_ref = par.Qref / series.Qref
  B_ref = par.Bref / series.Bref
  rho_ref = SQRT(m_ref*T_ref/Q_ref) / B_ref
  c_ref = SQRT(T_ref*Q_ref/m_ref)

  D_gb_norm = c_ref * rho_ref^2 / L_ref
  Gamma_norm = D_gb_norm * n_ref / L_ref
  Q_norm = D_gb_norm * n_ref / L_ref * T_ref * Q_ref

  Gamma_norm_i = Gamma_norm
  Gamma_norm_e = Gamma_norm
  Q_norm_i = Q_norm
  Q_norm_e = Q_norm
  E_norm = T_ref / L_ref
  curr_norm = c_ref * Q_ref * n_ref * rho_ref / L_ref
  temp_norm = T_ref * rho_ref / (3.0 * L_ref) ; already includes factor 3

  ; --- compute derived quantities ---

  xarr = (FINDGEN(par.nx0) / par.nx0 - 0.5) * par.lx
  ky = par.kymin * INDGEN(par.nky0)
  ; ve_x = - d phi / dy
  ve_x = - COMPLEX(0.0,1.0) * REBIN(REFORM(ky/series.Bref,[1,par.nky0,1]),$
    [par.nkx0,par.nky0,par.nz0]) * (*mom[0,0].kxky)
  ; B_x = ehat * beta * dApar/dy
  B_x = nf GT 1 ? COMPLEX(0.0,1.0) * REBIN(REFORM(ky/series.Bref,$
    [1,par.nky0,1]),[par.nkx0,par.nky0,par.nz0]) * (*mom[0,1].kxky) : 0.0

  Gamma_es = CONJ(*mom[0,nf].kxky) * ve_x * spec[0].dens
  Gamma_em = nf GT 1 ? CONJ(*mom[0,nf+5].kxky) * B_x * spec[0].dens : $
    FLTARR(par.nkx0,par.nky0,par.nz0)
  ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
  Q_e_es = CONJ(0.5*(*mom[s_e,nf+1].kxky)+(*mom[s_e,nf+2].kxky)+$
    1.5*(*mom[s_e,nf].kxky)) * ve_x * (spec[sp_e].dens * spec[sp_e].temp)
  ; Q_em = (qpar + qperp + 5/2 upar) B_x   
  Q_e_em = nf GT 1 ? CONJ((*mom[s_e,nf+3].kxky)+(*mom[s_e,nf+4].kxky)) * B_x * $
    (spec[sp_e].dens * spec[sp_e].temp) : FLTARR(par.nkx0,par.nky0,par.nz0)
  ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
  Q_i_es = CONJ(0.5*(*mom[s_i,nf+1].kxky)+(*mom[s_i,nf+2].kxky)+$
    1.5*(*mom[s_i,nf].kxky)) * ve_x * (spec[sp_i].dens * spec[sp_i].temp)
  ; Q_em = (qpar + qperp + 5/2 upar) B_x   
  Q_i_em = nf GT 1 ? CONJ((*mom[s_i,nf+3].kxky)+(*mom[s_i,nf+4].kxky)) * B_x * $
    (spec[sp_i].dens * spec[sp_i].temp) : FLTARR(par.nkx0,par.nky0,par.nz0)
  ; real part for kx averaging, normalization
  Gamma_e_es = FLOAT(Gamma_es) * Gamma_norm_e
  Gamma_i_es = FLOAT(TEMPORARY(Gamma_es)) * Gamma_norm_i
  Gamma_e_em = FLOAT(Gamma_em) * Gamma_norm_e
  Gamma_i_em = FLOAT(TEMPORARY(Gamma_em)) * Gamma_norm_i
  Q_e_es = FLOAT(TEMPORARY(Q_e_es)) * Q_norm_e
  Q_i_es = FLOAT(TEMPORARY(Q_i_es)) * Q_norm_i
  Q_e_em = FLOAT(TEMPORARY(Q_e_em)) * Q_norm_e
  Q_i_em = FLOAT(TEMPORARY(Q_i_em)) * Q_norm_i

  radE_kxky = COMPLEXARR(par.nkx0,par.nky0,par.nz0)
  radE_kxky[0,0,0] = COMPLEX(0,1) * REBIN((*series.kx)[0:par.nkx0/2-1,0],$
    [par.nkx0/2,par.nky0,par.nz0]) * (*mom[s_e,0].kxky)[0:par.nkx0/2-1,*,*]
  radE_kxky[par.nkx0/2+1,0,0] = $
    COMPLEX(0,1) * REBIN((*series.kx)[par.nkx0/2+1:par.nkx0-1,0],$
    [par.nkx0-par.nkx0/2-1,par.nky0,par.nz0]) * $
    (*mom[s_e,0].kxky)[par.nkx0/2+1:par.nkx0-1,*,*]
  radialE_sxky = FFT(TEMPORARY(radE_kxky),DIMENSION=1,/INVERSE) * E_norm

  parcurrent = (*mom[s_e,nf+5].kxky) * (spec[sp_e].charge * curr_norm) + $
    (*mom[s_i,nf+5].kxky) * (spec[sp_i].charge * curr_norm)
  parcurrent_sxky = (*mom[s_e,nf+5].sxky) * (spec[sp_e].charge * curr_norm) + $
    (*mom[s_i,nf+5].sxky) * (spec[sp_i].charge * curr_norm)
  parcurrent_sxsy = (*mom[s_e,nf+5].sxsy) * (spec[sp_e].charge * curr_norm) + $
    (*mom[s_i,nf+5].sxsy) * (spec[sp_i].charge * curr_norm)

  temp_e = ((*mom[s_e,nf+1].kxky) + 2.0 * (*mom[s_e,nf+2].kxky)) * temp_norm
  temp_e_sxky = ((*mom[s_e,nf+1].sxky) + 2.0 * (*mom[s_e,nf+2].sxky)) * temp_norm
  temp_e_sxsy = ((*mom[s_e,nf+1].sxsy) + 2.0 * (*mom[s_e,nf+2].sxsy)) * temp_norm
  temp_i = ((*mom[s_i,nf+1].kxky) + 2.0 * (*mom[s_i,nf+2].kxky)) * temp_norm
  temp_i_sxky = ((*mom[s_i,nf+1].sxky) + 2.0 * (*mom[s_i,nf+2].sxky)) * temp_norm
  temp_i_sxsy = ((*mom[s_i,nf+1].sxsy) + 2.0 * (*mom[s_i,nf+2].sxsy)) * temp_norm

  vorticity = - (spec[sp_i].charge * curr_norm) * (*mom[s_i,nf].kxky) - $
    (spec[sp_e].charge * curr_norm) * (*mom[s_e,nf].kxky)
  vorticity_sxky = - (spec[sp_i].charge * curr_norm) * (*mom[s_i,nf].sxky) - $
    (spec[sp_e].charge * curr_norm) * (*mom[s_e,nf].sxky)
  vorticity_sxsy = - (spec[sp_i].charge * curr_norm) * (*mom[s_i,nf].sxsy) - $
    (spec[sp_e].charge * curr_norm) * (*mom[s_e,nf].sxsy)

  ; --- store data ---

  (*i).time_id = store_step((*i).time_id,mom_time)
  ; --------------------------
  ; not available, set to zero
  energy_exb = 0.0
  (*i).energy_exb_id = store_step((*i).energy_exb_id,energy_exb)
  energy_mag = 0.0
  (*i).energy_mag_id = store_step((*i).energy_mag_id,energy_mag)
  energy_electron_th = 0.0
  (*i).energy_electron_th_id = store_step((*i).energy_electron_th_id,energy_electron_th)
  energy_ion_th = 0.0
  (*i).energy_ion_th_id = store_step((*i).energy_ion_th_id,energy_ion_th)
  energy = 0.0
  (*i).energy_id = store_step((*i).energy_id,energy)
  energy_electron_par = 0.0
  (*i).energy_electron_par_id = store_step((*i).energy_electron_par_id,energy_electron_par)
  energy_ion_par = 0.0
  (*i).energy_ion_par_id = store_step((*i).energy_ion_par_id,energy_ion_par)
  ; --------------------------
  (*i).flux_particle_id = store_step((*i).flux_particle_id,avg(Gamma_e_es,/all))
  (*i).flux_heat_el_id = store_step((*i).flux_heat_el_id,avg(Q_e_es,/all))
  (*i).flux_ion_id = store_step((*i).flux_ion_id,avg(Gamma_i_es,/all))
  (*i).flux_heat_ion_id = store_step((*i).flux_heat_ion_id,avg(Q_i_es,/all))
  (*i).flux_mag_particle_id = store_step((*i).flux_mag_particle_id,avg(Gamma_e_em,/all))
  (*i).flux_mag_ion_id = store_step((*i).flux_mag_ion_id,avg(Gamma_i_em,/all))
  (*i).flux_mag_heat_el_id = store_step((*i).flux_mag_heat_el_id,avg(Q_e_em,/all))
  (*i).flux_mag_heat_ion_id = store_step((*i).flux_mag_heat_ion_id,avg(Q_i_em,/all))
  (*i).zonal_ne_id = time_avg((*i).zonal_ne_id,avg(*mom[s_e,nf].sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_ni_id = time_avg((*i).zonal_ni_id,avg(*mom[s_i,nf].sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_Jpl_id = time_avg((*i).zonal_Jpl_id,avg(parcurrent_sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_phi_id = time_avg((*i).zonal_phi_id,avg(*mom[s_e,0].sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_Er_id = time_avg((*i).zonal_Er_id,avg(radialE_sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_vor_id = time_avg((*i).zonal_vor_id,avg(vorticity_sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_Apl_id = time_avg((*i).zonal_Apl_id,avg(*mom[s_e,1].sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_te_id = time_avg((*i).zonal_te_id,avg(temp_e_sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_ti_id = time_avg((*i).zonal_ti_id,avg(temp_i_sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).zonal_ui_id = time_avg((*i).zonal_ui_id,avg(*mom[s_i,nf+5].sxky,/zon),mom_time,fwd=gui.out.res_steps)
  (*i).spec_ne_id = time_avg((*i).spec_ne_id,avg(ABS(*mom[s_e,nf].kxky)^2,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_te_id = time_avg((*i).spec_te_id,avg(ABS(temp_e)^2,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_ti_id = time_avg((*i).spec_ti_id,avg(ABS(temp_i)^2,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_phi_id = time_avg((*i).spec_phi_id,avg(ABS(*mom[s_e,0].kxky)^2,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_vor_id = time_avg((*i).spec_vor_id,avg(ABS(vorticity)^2,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_Jpl_id = time_avg((*i).spec_Jpl_id,avg(ABS(parcurrent)^2,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_Fe_id = time_avg((*i).spec_Fe_id,avg(Gamma_e_es,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_Qe_id = time_avg((*i).spec_Qe_id,avg(Q_e_es,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_Qi_id = time_avg((*i).spec_Qi_id,avg(Q_i_es,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_Me_id = time_avg((*i).spec_Me_id,avg(Q_e_em,/kys),mom_time,fwd=gui.out.res_steps)
  (*i).spec_Mi_id = time_avg((*i).spec_Mi_id,avg(Q_i_em,/kys),mom_time,fwd=gui.out.res_steps)
  ; --------------------------
  ; not available, set to zero
  spec_B = FLTARR(par.nx0,par.nky0,par.nz0)
  (*i).spec_B_id = time_avg((*i).spec_B_id,avg(ABS(spec_B)^2,/kys),mom_time,fwd=gui.out.res_steps)
  he = FLTARR(par.nx0,par.nky0,par.nz0)
  (*i).env_he_id = time_avg((*i).env_he_id,avg(ABS(he)^2,/env),mom_time,fwd=gui.out.res_steps)
  ; --------------------------
  (*i).env_ne_id = time_avg((*i).env_ne_id,avg(ABS(*mom[s_e,nf].kxky)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_te_id = time_avg((*i).env_te_id,avg(ABS(temp_e)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_Fe_id = time_avg((*i).env_Fe_id,avg(Gamma_e_es,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_Qe_id = time_avg((*i).env_Qe_id,avg(Q_e_es,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_Me_id = time_avg((*i).env_Me_id,avg(Q_e_em,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_ni_id = time_avg((*i).env_ni_id,avg(ABS(*mom[s_i,nf].kxky)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_ti_id = time_avg((*i).env_ti_id,avg(ABS(temp_i)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_ui_id = time_avg((*i).env_ui_id,avg(ABS(*mom[s_i,nf+5].kxky)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_Qi_id = time_avg((*i).env_Qi_id,avg(Q_i_es,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_Mi_id = time_avg((*i).env_Mi_id,avg(Q_i_em,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_vor_id = time_avg((*i).env_vor_id,avg(ABS(vorticity)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_Jpl_id = time_avg((*i).env_Jpl_id,avg(ABS(parcurrent)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).env_phi_id = time_avg((*i).env_phi_id,avg(ABS(*mom[s_e,0].kxky)^2,/env),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_phi_id = time_avg((*i).axisym_phi_id,avg(*mom[s_e,0].sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_Apl_id = time_avg((*i).axisym_Apl_id,avg(*mom[s_e,1].sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_ne_id = time_avg((*i).axisym_ne_id,avg(*mom[s_e,nf].sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_te_id = time_avg((*i).axisym_te_id,avg(temp_e_sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_ni_id = time_avg((*i).axisym_ni_id,avg(*mom[s_i,nf].sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_ti_id = time_avg((*i).axisym_ti_id,avg(temp_i_sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_ui_id = time_avg((*i).axisym_ui_id,avg(*mom[s_i,nf+5].sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_vor_id = time_avg((*i).axisym_vor_id,avg(vorticity_sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).axisym_Jpl_id = time_avg((*i).axisym_Jpl_id,avg(parcurrent_sxsy,/axi),mom_time,fwd=gui.out.res_steps)
  (*i).phi_id = time_avg((*i).phi_id,1.0*(*mom[s_e,0].sxsy),mom_time,fwd=gui.out.res_steps)
  (*i).ne_id = time_avg((*i).ne_id,1.0*(*mom[s_e,nf].sxsy),mom_time,fwd=gui.out.res_steps)
  (*i).Jpl_id = time_avg((*i).Jpl_id,parcurrent_sxsy,mom_time,fwd=gui.out.res_steps)
  (*i).vor_id = time_avg((*i).vor_id,vorticity_sxsy,mom_time,fwd=gui.out.res_steps)

END

;######################################################################

PRO store_var, var, varname, comment, vartype
; stores a variable var in a dataset named varname and adds a comment

  group_id = vartype[0]
  type_id = vartype[1]
  space_id = vartype[2]

  var_id = H5D_CREATE(group_id,varname,type_id,space_id)

  comment_type = H5T_IDL_CREATE(STRJOIN(REPLICATE(' ',128),/SINGLE))
  comment_space = H5S_CREATE_SCALAR()
  comment_id = H5A_CREATE(var_id,'comment',comment_type,comment_space)
  H5A_WRITE, comment_id, comment
  H5A_CLOSE, comment_id
  H5S_CLOSE, comment_space
  H5T_CLOSE, comment_type

  H5D_WRITE, var_id, var
  H5D_CLOSE, var_id

END

;######################################################################

PRO itmhdf5_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  sp_e = (*i).sp_e
  sp_i = (*i).sp_i

  n_ref = series.nref * par.nref
  L_ref = series.Lref * par.Lref
  m_ref = series.mref * par.mref
  T_ref = series.Tref * par.Tref
  Q_ref = series.Qref * par.Qref
  B_ref = series.Bref * par.Bref
  rho_ref = SQRT(m_ref*T_ref/Q_ref) / B_ref
  c_ref = SQRT(T_ref*Q_ref/m_ref)
  time_norm = L_ref / c_ref
  dens_norm = n_ref * rho_ref / L_ref
  phi_norm = T_ref * rho_ref / L_ref
  Apar_norm = B_ref * rho_ref^2 / L_ref
  u_norm = c_ref * rho_ref / L_ref

  ; --- get time averages/time resolved arrays of simulation data ---

  time = store_step((*i).time_id,/get,/reg_array) * time_norm
  energy_exb = store_step((*i).energy_exb_id,/get,/reg_array)
  energy_mag = store_step((*i).energy_mag_id,/get,/reg_array)
  energy_electron_th = store_step((*i).energy_electron_th_id,/get,/reg_array)
  energy_ion_th = store_step((*i).energy_ion_th_id,/get,/reg_array)
  energy_electron_par = store_step((*i).energy_electron_par_id,/get,/reg_array)
  energy_ion_par = store_step((*i).energy_ion_par_id,/get,/reg_array)
  energy = store_step((*i).energy_id,/get,/reg_array)
  flux_particle = store_step((*i).flux_particle_id,/get,/reg_array)
  flux_heat_el = store_step((*i).flux_heat_el_id,/get,/reg_array)
  flux_ion = store_step((*i).flux_ion_id,/get,/reg_array)
  flux_heat_ion = store_step((*i).flux_heat_ion_id,/get,/reg_array)
  flux_mag_particle = store_step((*i).flux_mag_particle_id,/get,/reg_array)
  flux_mag_ion = store_step((*i).flux_mag_ion_id,/get,/reg_array)
  flux_mag_heat_el = store_step((*i).flux_mag_heat_el_id,/get,/reg_array)
  flux_mag_heat_ion = store_step((*i).flux_mag_heat_ion_id,/get,/reg_array)
  zonal_ne = time_avg((*i).zonal_ne_id,/avg,fwd=gui.out.res_steps) * dens_norm
  zonal_ni = time_avg((*i).zonal_ni_id,/avg,fwd=gui.out.res_steps) * dens_norm
  zonal_Jpl = time_avg((*i).zonal_Jpl_id,/avg,fwd=gui.out.res_steps)
  zonal_phi = time_avg((*i).zonal_phi_id,/avg,fwd=gui.out.res_steps) * phi_norm
  zonal_Er = time_avg((*i).zonal_Er_id,/avg,fwd=gui.out.res_steps)
  zonal_vor = time_avg((*i).zonal_vor_id,/avg,fwd=gui.out.res_steps)
  zonal_Apl = time_avg((*i).zonal_Apl_id,/avg,fwd=gui.out.res_steps) * Apar_norm
  zonal_te = time_avg((*i).zonal_te_id,/avg,fwd=gui.out.res_steps)
  zonal_ti = time_avg((*i).zonal_ti_id,/avg,fwd=gui.out.res_steps)
  zonal_ui = time_avg((*i).zonal_ui_id,/avg,fwd=gui.out.res_steps) * u_norm
  spec_ne = time_avg((*i).spec_ne_id,/avg,fwd=gui.out.res_steps) * dens_norm^2
  spec_te = time_avg((*i).spec_te_id,/avg,fwd=gui.out.res_steps)
  spec_ti = time_avg((*i).spec_ti_id,/avg,fwd=gui.out.res_steps)
  spec_phi = time_avg((*i).spec_phi_id,/avg,fwd=gui.out.res_steps) * phi_norm^2
  spec_vor = time_avg((*i).spec_vor_id,/avg,fwd=gui.out.res_steps)
  spec_Jpl = time_avg((*i).spec_Jpl_id,/avg,fwd=gui.out.res_steps)
  spec_Fe = time_avg((*i).spec_Fe_id,/avg,fwd=gui.out.res_steps)
  spec_Qe = time_avg((*i).spec_Qe_id,/avg,fwd=gui.out.res_steps)
  spec_Qi = time_avg((*i).spec_Qi_id,/avg,fwd=gui.out.res_steps)
  spec_Me = time_avg((*i).spec_Me_id,/avg,fwd=gui.out.res_steps)
  spec_Mi = time_avg((*i).spec_Mi_id,/avg,fwd=gui.out.res_steps)
  spec_B = time_avg((*i).spec_B_id,/avg,fwd=gui.out.res_steps)
  env_ne = time_avg((*i).env_ne_id,/avg,fwd=gui.out.res_steps) * dens_norm^2
  env_he = time_avg((*i).env_he_id,/avg,fwd=gui.out.res_steps)
  env_te = time_avg((*i).env_te_id,/avg,fwd=gui.out.res_steps)
  env_Fe = time_avg((*i).env_Fe_id,/avg,fwd=gui.out.res_steps)
  env_Qe = time_avg((*i).env_Qe_id,/avg,fwd=gui.out.res_steps)
  env_Me = time_avg((*i).env_Me_id,/avg,fwd=gui.out.res_steps)
  env_ni = time_avg((*i).env_ni_id,/avg,fwd=gui.out.res_steps) * dens_norm^2
  env_ti = time_avg((*i).env_ti_id,/avg,fwd=gui.out.res_steps)
  env_ui = time_avg((*i).env_ui_id,/avg,fwd=gui.out.res_steps) * u_norm^2
  env_Qi = time_avg((*i).env_Qi_id,/avg,fwd=gui.out.res_steps)
  env_Mi = time_avg((*i).env_Mi_id,/avg,fwd=gui.out.res_steps)
  env_vor = time_avg((*i).env_vor_id,/avg,fwd=gui.out.res_steps)
  env_Jpl = time_avg((*i).env_Jpl_id,/avg,fwd=gui.out.res_steps)
  env_phi = time_avg((*i).env_phi_id,/avg,fwd=gui.out.res_steps) * phi_norm^2
  axisym_phi = time_avg((*i).axisym_phi_id,/avg,fwd=gui.out.res_steps) * phi_norm
  axisym_Apl = time_avg((*i).axisym_Apl_id,/avg,fwd=gui.out.res_steps) * Apar_norm
  axisym_ne = time_avg((*i).axisym_ne_id,/avg,fwd=gui.out.res_steps) * dens_norm
  axisym_te = time_avg((*i).axisym_te_id,/avg,fwd=gui.out.res_steps)
  axisym_ni = time_avg((*i).axisym_ni_id,/avg,fwd=gui.out.res_steps) * dens_norm
  axisym_ti = time_avg((*i).axisym_ti_id,/avg,fwd=gui.out.res_steps)
  axisym_ui = time_avg((*i).axisym_ui_id,/avg,fwd=gui.out.res_steps) * u_norm
  axisym_vor = time_avg((*i).axisym_vor_id,/avg,fwd=gui.out.res_steps)
  axisym_Jpl = time_avg((*i).axisym_Jpl_id,/avg,fwd=gui.out.res_steps)
  phi = time_avg((*i).phi_id,/avg,fwd=gui.out.res_steps) * phi_norm
  xne = time_avg((*i).ne_id,/avg,fwd=gui.out.res_steps) * dens_norm
  Jpl = time_avg((*i).Jpl_id,/avg,fwd=gui.out.res_steps)
  vor = time_avg((*i).vor_id,/avg,fwd=gui.out.res_steps)

  ; --- create HDF5 file and its directory structure ---

  h5file = gui.out.data_path + 'GENE.h5'
  PRINT, 'writing data to file ' + h5file
  FILE_DELETE, h5file, /ALLOW_NONEXISTENT
  file_id = H5F_CREATE(h5file)

  data_group_id = H5G_CREATE(file_id,'data')
  var0d_group_id = H5G_CREATE(data_group_id,'var0d')
  var1d_group_id = H5G_CREATE(data_group_id,'var1d')
  var2d_group_id = H5G_CREATE(data_group_id,'var2d')
  var3d_group_id = H5G_CREATE(data_group_id,'var3d')
  env1d_group_id = H5G_CREATE(data_group_id,'env1d')
  spect1d_group_id = H5G_CREATE(data_group_id,'spec1d')
  coord_sys_group_id = H5G_CREATE(data_group_id,'coord_sys')
  turbgrid_group_id = H5G_CREATE(coord_sys_group_id,'turbgrid')
  documentation_group_id = H5G_CREATE(file_id,'documentation')

  ; --- define data types and data spaces ---

  step_count = gui.out.res_steps ? series.step_count : 1

  float_type_id = H5T_IDL_CREATE(0.0)
  scalar_space_id = H5S_CREATE_SCALAR()
  t_space_id = H5S_CREATE_SIMPLE(series.step_count)
  x_space_id = H5S_CREATE_SIMPLE(par.nx0)
  xt_space_id = H5S_CREATE_SIMPLE([par.nx0,step_count])
  y_space_id = H5S_CREATE_SIMPLE(par.nky0)
  yt_space_id = H5S_CREATE_SIMPLE([par.nky0,step_count])
  z_space_id = H5S_CREATE_SIMPLE(par.nz0)
  zt_space_id = H5S_CREATE_SIMPLE([par.nz0,step_count])
  v_space_id = H5S_CREATE_SIMPLE(par.nv0)
  w_space_id = H5S_CREATE_SIMPLE(par.nw0)
  xy_space_id = H5S_CREATE_SIMPLE([par.nx0,par.ny0])
  xyt_space_id = H5S_CREATE_SIMPLE([par.nx0,par.ny0,step_count])
  xyz_space_id = H5S_CREATE_SIMPLE([par.nx0,par.ny0,par.nz0])
  xyzt_space_id = H5S_CREATE_SIMPLE([par.nx0,par.ny0,par.nz0,step_count])

  ; --- convert and store documentation and parameter file ---

  doc_path = '../doc/gene.pdf'
  doc_exists = FILE_TEST(doc_path)
  IF doc_exists THEN BEGIN
    n_lines = FILE_LINES(doc_path)
    OPENR, doc_lun, doc_path, /GET_LUN
    content_doc = STRARR(n_lines)
    READF, doc_lun, content_doc
    FREE_LUN, doc_lun

    content_doc[0:n_lines-2] += '\n'
    content_doc = STRJOIN(TEMPORARY(content_doc),/SINGLE)
    ; currently, no method exists for preserving/re-extracting a pdf

    doc_type_id = H5T_IDL_CREATE(content_doc)
    vartype = [documentation_group_id,doc_type_id,scalar_space_id]
    store_var, content_doc, 'Code_Paper', 'GENE documentation', vartype
    H5T_CLOSE, doc_type_id
  ENDIF ELSE PRINT, 'no pdf documentation found, skipping'

  par_path = file_par.path[0]
  n_lines = FILE_LINES(par_path)
  OPENR, par_lun, par_path, /GET_LUN
  content_par = STRARR(n_lines)
  READF, par_lun, content_par
  FREE_LUN, par_lun

  content_par[0:n_lines-2] += '\n'
  content_par = STRJOIN(TEMPORARY(content_par),/SINGLE)
  ; note: this way, one can obtain the original file by importing the
  ; string into IDL and typing
  ; SPAWN, 'echo "' + Parameter_InputFile._DATA + '" > parfile'

  par_type_id = H5T_IDL_CREATE(content_par)
  vartype = [documentation_group_id,par_type_id,scalar_space_id]
  store_var, content_par, 'Parameter_InputFile', 'GENE parameter file', vartype
  H5T_CLOSE, par_type_id

  ; --- store simulation data ---

  ; grid
  grid_type = 'local flux tube, field aligned'
  string_type_id = H5T_IDL_CREATE(grid_type)
  vartype = [coord_sys_group_id,string_type_id,scalar_space_id]
  store_var, grid_type, 'grid_type', 'type of coordinate system', vartype
  H5T_CLOSE, string_type_id

  vartype = [turbgrid_group_id,float_type_id,x_space_id]
  store_var, (FINDGEN(par.nx0) / par.nx0 - 0.5) * par.lx, 'dim1', 'first dimension of grid', vartype
  vartype = [turbgrid_group_id,float_type_id,y_space_id]
  store_var, *series.ky, 'dim2', 'second dimension of grid', vartype
  vartype = [turbgrid_group_id,float_type_id,z_space_id]
  store_var, (2 * FINDGEN(par.nz0) / par.nz0 - 1) * !PI, 'dim3', 'third dimension of grid', vartype
  vartype = [turbgrid_group_id,float_type_id,v_space_id]
  store_var, (2 * FINDGEN(par.nv0) / (par.nv0 - 1) - 1) * par.lv, 'dim_v1', 'first dimension of velocity space grid', vartype
  vartype = [turbgrid_group_id,float_type_id,w_space_id]
  getMuWeightsAndKnots, mu_weights, mu_knots
  store_var, mu_knots, 'dim_v2', 'second dimension of velocity space grid', vartype
  vartype = [turbgrid_group_id,float_type_id,y_space_id]
  store_var, *series.ky, 'dim_spec', 'dimension of spectrum wavenumber grid', vartype
  vartype = [turbgrid_group_id,float_type_id,z_space_id]
  store_var, (2 * FINDGEN(par.nz0) / par.nz0 - 1) * !PI, 'dim_env', 'dimension of parallel envelope grid', vartype

  ; time
  vartype = [var0d_group_id,float_type_id,t_space_id]
  store_var, time, 'time', 'code time for this segment [s]', vartype

  ; 0d data
  vartype = [var0d_group_id,float_type_id,t_space_id]
  store_var, energy_exb, 'energy_exb', 'volume average ExB energy density in [J/m^3]', vartype
  store_var, energy_mag, 'energy_mag', 'volume average magnetic energy density in [J/m^3]', vartype
  store_var, energy_electron_th, 'energy_electron_th', $
    'volume average electron thermal free energy density in [J/m^3]', vartype
  store_var, energy_ion_th, 'energy_ion_th', 'volume average ion thermal free energy density in [J/m^3]', vartype
  store_var, energy_electron_par, 'energy_electron_par', 'volume average electron parallel energy in [J/m^3]', vartype
  store_var, energy_ion_par, 'energy_ion_par', 'volume average ion parallel energy in [J/m^3]', vartype
  store_var, energy, 'energy', 'volume average total free energy density in [J/m^3]', vartype
  store_var, flux_particle, 'flux_particle', 'volume average electron particle flux in [m^-2/s]', vartype
  store_var, flux_heat_el, 'flux_heat_el', 'volume average electron conductive heat flux in [eV m^-2/s]', vartype
  store_var, flux_ion, 'flux_ion', 'volume average ion particle flux in [m^-2/s]', vartype
  store_var, flux_heat_ion, 'flux_heat_ion', 'volume average ion conductive heat flux in [eV m^-2/s]', vartype
  store_var, flux_mag_particle, 'flux_mag_particle', 'volume average magnetic flutter particle flux in [m^-2/s]', vartype
  store_var, flux_mag_heat_el, 'flux_mag_heat_el', $
    'volume average magnetic flutter electron conductive heat flux in [eV m^-2/s]', vartype
  store_var, flux_mag_ion, 'flux_mag_ion', 'volume average magnetic flutter ion flux in [m^-2/s]', vartype
  store_var, flux_mag_heat_ion, 'flux_mag_heat_ion', $
    'volume average magnetic flutter ion conductive heat flux in [eV m^-2/s]', vartype

  ; 1d x data
  var1d = [var1d_group_id,float_type_id,xt_space_id]
  store_var, zonal_phi, 'zonal_phi', 'zonal electric potential in [V]', var1d
  store_var, zonal_Er, 'zonal_Er', 'zonal radial electric field in [V/m]', var1d
  store_var, zonal_vor, 'zonal_vor', 'zonal vorticity in [coulombs/m^3]', var1d
  store_var, zonal_Apl, 'zonal_Apl', 'zonal parallel magnetic potential in [T m]', var1d
  store_var, zonal_Jpl, 'zonal_Jpl', 'zonal parallel current divided by B in [A/m^2 per T]', var1d
  store_var, zonal_ne, 'zonal_ne', 'zonal electron density in [m^-3]', var1d
  store_var, zonal_te, 'zonal_te', 'zonal electron temperature in [eV]', var1d
  store_var, zonal_ni, 'zonal_ni', 'zonal ion density in [m^-3]', var1d
  store_var, zonal_ti, 'zonal_ti', 'zonal ion temperature in [eV]', var1d
  store_var, zonal_ui, 'zonal_ui', 'zonal parallel ion velocity divided by B in [m/s per T]', var1d

  ; 1d ky data
  vartype = [spect1d_group_id,float_type_id,yt_space_id]
  store_var, spec_phi, 'spec_phi', 'electrostatic potential', vartype
  store_var, spec_vor, 'spec_vor', 'vorticity', vartype
  store_var, spec_B, 'spec_B', 'magnetic energy', vartype
  store_var, spec_Jpl, 'spec_Jpl', 'parallel current', vartype
  store_var, spec_ne, 'spec_ne', 'electron density', vartype
  store_var, spec_te, 'spec_te', 'electron temperature', vartype
  store_var, spec_Fe, 'spec_Fe', 'particle flux', vartype
  store_var, spec_Qe, 'spec_Qe', 'electron heat flux', vartype
  store_var, spec_Me, 'spec_Me', 'electron magnetic heat flux', vartype
  store_var, spec_ti, 'spec_ti', 'ion temperature', vartype
  store_var, spec_Qi, 'spec_Qi', 'ion heat flux', vartype
  store_var, spec_Mi, 'spec_Mi', 'ion magnetic heat flux', vartype

  ; 1d z data
  vartype = [env1d_group_id,float_type_id,zt_space_id]
  store_var, env_phi, 'env_phi', 'parallel envelope electrostatic potential in [V^2]', vartype
  store_var, env_vor, 'env_vor', 'parallel envelope vorticity in [coulomb^2/m^6]', vartype
  store_var, env_Jpl, 'env_Jpl', 'parallel envelope parallel current in [A^2/m^4]', vartype
  store_var, env_ne, 'env_ne', 'parallel envelope electron density in [m^-6]', vartype
  store_var, env_he , 'env_he', 'parallel envelope nonadiabatic electron density in [m^-6]', vartype
  store_var, env_te , 'env_te', 'parallel envelope electron temperature in [eV^2]', vartype
  store_var, env_Fe , 'env_Fe', 'parallel envelope electron flux in [m^-2/s]', vartype
  store_var, env_Qe , 'env_Qe', 'parallel envelope electron heat flux in [eV m^-2/s]', vartype
  store_var, env_Me , 'env_Me', 'parallel envelope magnetic electron heat flux in [eV m^-2/s]', vartype
  store_var, env_ni , 'env_ni', 'parallel envelope ion density in [m^-6]', vartype
  store_var, env_ti , 'env_ti', 'parallel envelope ion temperature in [eV^2]', vartype
  store_var, env_ui , 'env_ui', 'parallel envelope ion parallel velocity in [m^2/s^2]', vartype
  store_var, env_Qi , 'env_Qi', 'parallel envelope ion heat flux in [eV m^-2/s]', vartype
  store_var, env_Mi , 'env_Mi', 'parallel envelope magnetic ion heat flux in [eV m^-2/s]', vartype

  ; 2d xy data
  vartype = [var2d_group_id,float_type_id,xyt_space_id]
  store_var, axisym_phi , 'axisym_phi', 'axisymmetric component of electrostatic potential in [V]', vartype
  store_var, axisym_Apl , 'axisym_Apl', 'axisymmetric component of parallel magnetic potential in [T m]', vartype
  store_var, axisym_Jpl , 'axisym_Jpl', 'axisymmetric component of parallel current in [A/m^2]', vartype
  store_var, axisym_vor, 'axisym_vor', 'axisymmetric component of vorticity in [coulomb/m^3]', vartype
  store_var, axisym_ne, 'axisym_ne', 'axisymmetric component of electron density in [m^-3]', vartype
  store_var, axisym_te , 'axisym_te', 'axisymmetric component of electron temperature in [eV]', vartype
  store_var, axisym_ni , 'axisym_ni', 'axisymmetric component of ion density in [m^-3]', vartype
  store_var, axisym_ti , 'axisym_ti', 'axisymmetric component of ion temperature in [eV]', vartype
  store_var, axisym_ui , 'axisym_ui', 'axisymmetric component of ion parallel velocity in [m/s]', vartype

  ; 3d xyz data
  vartype = [var3d_group_id,float_type_id,xyzt_space_id]
  store_var, phi, 'phi', 'electrostatic potential in [V]', vartype
  store_var, vor, 'vor', 'vorticity in [coulombs/m^3]', vartype
  store_var, Jpl, 'Jpl', 'parallel current in [A/m^2]', vartype
  store_var, xne, 'ne', 'electron density in [m^-3]', vartype

  ; --- HDF5 data type and space cleanup ---

  H5S_CLOSE, scalar_space_id
  H5S_CLOSE, t_space_id
  H5S_CLOSE, x_space_id
  H5S_CLOSE, xt_space_id
  H5S_CLOSE, y_space_id
  H5S_CLOSE, yt_space_id
  H5S_CLOSE, z_space_id
  H5S_CLOSE, zt_space_id
  H5S_CLOSE, v_space_id
  H5S_CLOSE, w_space_id
  H5S_CLOSE, xy_space_id
  H5S_CLOSE, xyt_space_id
  H5S_CLOSE, xyz_space_id
  H5S_CLOSE, xyzt_space_id
  H5T_CLOSE, float_type_id
  H5F_CLOSE, file_id

  gui.misc.recent_h5 = gui.misc.recent_h5 + ' ' + h5file

  set_output, diag, /reset

END
