FUNCTION torus_info

  RETURN, {$
    type      : 'mom',$ ;'mom_nl',$
    title     : 'Toroidal representation',$
    help_text : ['Transforms data from the field aligned grid '+$
                 'onto the torus.'],$
    ext_vars  : [$
      ['mode','0','mode of operation; 0: const. phi slice (default), '+$
       '1: partial torus, 2: PR textures, 3: VK data'],$
      ['var','0','variable to display; default: 0'],$
      ['phi_cut','0','toroidal position of cut(s) in units of 2 pi; '+$
       'default: 0'],$
      ['res','0','user defined resolutions: nz[, nx[, ny[,nphi[,tintrp]]]]; '+$
       'default : 400, nx0, ny0, 400, 0 where tinterp>0 enforces '+$
       'equidistant time slices with the number of time steps multiplied '+$
       'by tinterp'],$
      ['ps2png','1','create HQ png files (modes 0 and 1); '+$
       'requires postscript output to be set']]}

END

;######################################################################

PRO torus_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  ; -------------------------------------------------------------------
  ; free parameters which are not included in the external vars
  torR = 1.0      ; major R normalization, default: 1.0 or GENE value
  phi_range = 1.0 ; size of slice relative to full torus, default: 1.0
  glob_x_red  = 0 ; floor(0.8*par.nx0/2) ; omit inner/outer flux surfaces, default: 0
  show_minmax = 0 ; print global min/max of data for normalization
  fix_color_zero = 0 ; fix color for zero amplitude to center of color scheme
  contourfill = 1 ; fill contours or contours lines only; default: 1
  use_sxky = (contourfill EQ 0) ; use sxky data instead of sxsy; default: 0
                                ; note: sxky is more accurate and may
                                ; be required for contour lines; sxsy
                                ; should, however, be more efficient
  ; -------------------------------------------------------------------
  ; uncomment this line to set the normalization manually
  ; man_minmax = [-50,50]
  ; -------------------------------------------------------------------

  IF par.lx EQ 0 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + ': ky dependent lx'
    (*diag).selected = 0
    RETURN
  ENDIF

  ; geometry types: 0 is s-alpha, 1 is circular, 2 is TRACER
  IF (WHERE(ABS((*series.geom).R) GT 1e-5))[0] GE 0 THEN $
     geom_type = 2 ELSE IF (par.magn_geometry EQ 'circular') THEN $
        geom_type = 1 ELSE geom_type = 0
  IF (par.trpeps LT 1e-5) AND (geom_type NE 2) AND (par.x_local) $
    THEN BEGIN

    printerror, 'Skipping ' + (*diag).name + ': trpeps not finite'
    (*diag).selected = 0
    RETURN
  ENDIF

  IF NOT KEYWORD_SET(*e.mode) THEN *e.mode = 0
  IF (*e.mode EQ 3) AND gui.out.n_spec_sel NE 1 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + $
      ': mode 3 accepts only one species'
    (*diag).selected = 0
    RETURN
  ENDIF

  IF NOT KEYWORD_SET(*e.var) THEN *e.var = 0
  res = [par.nx0,par.ny0,400]
  phi_res = 400 ; number of phi points for the hull, default: 400
  tintrpfac = 0  ; time interpolation factor
  CASE N_ELEMENTS(*e.res) OF
    1    : res[2] = *e.res
    2    : res[*] = [(*e.res)[1],res[1],(*e.res)[0]]
    3    : res[*] = [(*e.res)[1],(*e.res)[2],(*e.res)[0]]
    4    : BEGIN
             res[*] = [(*e.res)[1],(*e.res)[2],(*e.res)[0]]
             phi_res = (*e.res)[3]
           END
    5    : BEGIN
             res[*] = [(*e.res)[1],(*e.res)[2],(*e.res)[0]]
             phi_res = (*e.res)[3]
             tintrpfac = (*e.res)[4]
           END
    ELSE :
  ENDCASE
  IF NOT KEYWORD_SET(*e.ps2png) THEN *e.ps2png = 0
  ps2png = (*e.ps2png OR (*e.mode NE 0)) AND (*e.mode NE 3)

  IF ((geom_type NE 1) AND (geom_type NE 2)) OR $
    (*e.mode EQ 0) OR (*e.mode EQ 3) THEN glob_x_red = 0

  IF N_ELEMENTS(*e.phi_cut) LT 1 THEN BEGIN
     IF *e.mode EQ 2 THEN *e.phi_cut = [0.0,0.5] $
     ELSE *e.phi_cut = 0.0
  ENDIF
  IF *e.mode EQ 1 THEN *e.phi_cut = (*e.phi_cut)[0]
  n_cuts = N_ELEMENTS(*e.phi_cut)
  *e.phi_cut *= 2 * !PI

  use_man_minmax = N_ELEMENTS(man_minmax) EQ 2

  get_torus_pars, geom_type, torR, torZ, r0, $
    round_q0, bigM, rho

  ;apply correction factor if C_y NE ;r0/q0; this is currently only
  ;necessary for Miller geometry, and therefore only considered in the code
  ;below for local simulations
  Cyq0_r0=1.0
  if par.magn_geometry eq 'miller' then Cyq0_r0=(*series.geom).C_y*par.q0/r0
  
  IF par.x_local AND (par.n0_global LT 0) THEN BEGIN
     IF ((bigM * round_q0) MOD 1 GT 1e-3) THEN BEGIN
        printerror, 'Skipping ' + (*diag).name + ': no proper bigM'
        (*diag).selected = 0
        RETURN
     ENDIF
  ENDIF

  i = set_internal_vars(diag,{$
    mode           : *e.mode,$
    var            : *e.var,$
    geom_type      : geom_type,$
    x_res          : res[0],$
    y_res          : res[1],$
    z_res          : res[2],$
    x_res_red      : res[0] - 2 * glob_x_red,$
    phi_res        : phi_res,$
    tintrpfac      : tintrpfac,$
    phi_cut        : *e.phi_cut MOD (2.0 * !PI),$
    n_cuts         : n_cuts,$
    phi_range      : phi_range,$
    torR           : torR,$
    torZ           : torZ,$
    r0      	   : r0,$
    dx             : 0.0,$
    dy             : 0.0,$
    q0             : round_q0,$
    bigM           : bigM,$
    y_shift        : FLTARR(res[0],/NOZERO),$
    glob_x_red     : glob_x_red,$
;    use_sh_met     : 0,$
;    shmet_corr     : FLTARR(res[0],res[2]+1,/NOZERO),$
    closest_y_cut  : *e.mode NE 3 ? pf_arr([2,res[2],n_cuts,res[0]]) : 0,$
    dy_cut         : *e.mode NE 3 ? pf_arr([2,res[2],n_cuts,res[0]]) : 0,$
    closest_y_hull : (*e.mode NE 3) AND (*e.mode NE 0) ? $
                     pf_arr([2,res[2],phi_res,res[0]]) : 0,$
    dy_hull        : (*e.mode NE 3) AND (*e.mode NE 0) ? $
                     pf_arr([2,res[2],phi_res,res[0]]) : 0,$
    space_r        : geom_type EQ 0 ? $
                     FLTARR(res[0],/NOZERO) : $
                     FLTARR(res[2],res[0]-2*glob_x_red,/NOZERO),$
    space_theta    : geom_type EQ 0 ? $
                     FLTARR(res[2],/NOZERO) : $
                     FLTARR(res[2],res[0]-2*glob_x_red,/NOZERO),$
    data_minmax    : FLTARR(2,gui.out.n_spec_sel),$
    show_minmax    : show_minmax,$
    use_man_minmax : use_man_minmax,$
    man_minmax     : use_man_minmax ? man_minmax : 0,$
    cut_data_id    : PTR_NEW(),$
    hull_data_id   : PTR_NEW(),$
    vk_data_id     : PTR_NEW(),$
    vk_data_lun    : 0,$
    rms            : 0,$ ;switch time avg. to rms
    png_out        : ps2png,$
    norm_local     : 1 - par.nonlinear * (ps2png OR (*e.mode EQ 3)),$
    transparent    : 1,$
    use_sxky       : use_sxky,$
    fix_color_zero : fix_color_zero,$
    contourfill    : contourfill})


  IF (*i).norm_local AND (*i).show_minmax THEN BEGIN
    PRINT, (*diag).name + $
      ': using per-step normalization, show_minmax = 0'
    (*i).show_minmax = 0
  ENDIF
  IF (*i).norm_local AND (*i).use_man_minmax THEN BEGIN
    PRINT, (*diag).name + $
      ': using per-step normalization, use_man_minmax = 0'
    (*i).use_man_minmax = 0
  ENDIF

  IF (*i).use_man_minmax THEN BEGIN
    (*i).man_minmax *= 1.0

    val_min = (*i).man_minmax[0]
    val_max = (*i).man_minmax[1]

    IF ABS(val_max-val_min) LT 1e-8 THEN BEGIN
      PRINT, (*diag).name + $
        ': manual min = max, use_man_minmax = 0'
      (*i).use_man_minmax = 0
    ENDIF

    IF val_min GT val_max THEN (*i).man_minmax = [val_max,val_min]
  ENDIF

  IF (*i).x_res_red LT 2 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + $
      ': reduced radial range to < 2 points'
    (*diag).selected = 0
    RETURN
  ENDIF

  (*i).dx = rho * par.lx / (*i).x_res ; normalized to Lref
  (*i).dy = rho * par.ly / (*i).y_res ; normalized to Lref

  use_sh_met = 0

  (*i).y_shift = get_y_shift((use_sxky?par.nx0:(*i).x_res), $
       (*i).dx, (*i).dy, Cyq0_r0)

  IF NOT par.x_local THEN BEGIN
    C_y = INTERPOL((*series.geom).C_y[*],(*i).x_res)
    q_prof = INTERPOL((*series.geom).q[*],(*i).x_res)

    IF par.shifted_metric THEN BEGIN
      ; prepare shifted metric correction
      gxx = REFORM((*series.geom).gxx,[par.nx0,1,par.nz0])
      interp_3d_data, gxx, (*i).x_res, 1, (*i).z_res, FLTARR((*i).x_res)
      gxx = REFORM(gxx,/OVERWRITE)
      gxz = REFORM((*series.geom).gxz,[par.nx0,1,par.nz0])
      interp_3d_data, gxz, (*i).x_res, 1, (*i).z_res, FLTARR((*i).x_res)
      gxz = REFORM(gxz,/OVERWRITE)

      use_sh_met = 1
    ENDIF
  ENDIF

  (*i).data_minmax[0,*] = 1e20
  (*i).data_minmax[1,*] = -1e20

  ; all cases except VK (mode 3): cut geometry
  ; for partial torus and PR: hull geometry
  zres_inv = 1.0 / (*i).z_res
  IF (*i).mode LT 3 THEN BEGIN
    ; calculating cut coordinates
    IF geom_type NE 2 THEN BEGIN
      r = (*i).r0 + (INDGEN((*i).x_res_red) - (*i).x_res_red / 2) * (*i).dx
      IF (par.Lref GT 0.0) THEN r*=par.Lref
      chi = (2 * INDGEN((*i).z_res) * zres_inv - 1.0) * !PI
      IF geom_type EQ 0 THEN BEGIN
        (*i).space_r[*] = r[*]
        (*i).space_theta[*] = chi[*]
      ENDIF ELSE BEGIN
        FOR z = 0, (*i).z_res - 1 DO BEGIN
          (*i).space_r[z,*] = r[*]
          (*i).space_theta[z,*] = $
            2.0 * ATAN(SQRT((1.0+r[*])/(1.0-r[*]))*TAN(chi[z]/2.0))
       ENDFOR
       ; IDL decides that ATAN(TAN(-pi/2)) is pi/2 which is why:
        (*i).space_theta[0,*]-= 2.0*!PI
       
      ENDELSE
    ENDIF ELSE BEGIN ; now: TRACER 2D case
      IF (par.x_local) THEN BEGIN
        Z_extd = [(*series.geom).Z,(*series.geom).Z[0]]
        Z_interp = (INTERPOL(Z_extd,(*i).z_res+1,/SPLINE))[0:(*i).z_res-1]
         
        R_extd = [(*series.geom).R,(*series.geom).R[0]]
        R_interp = (INTERPOL(R_extd,(*i).z_res+1,/SPLINE))[0:(*i).z_res-1]
         
        c1_extd = [(*series.geom).c1/(*series.geom).gxx,$
          (*series.geom).c1[0]/(*series.geom).gxx[0]]
        c1_interp = (INTERPOL(c1_extd,(*i).z_res+1,/SPLINE))[0:(*i).z_res-1]
        c2_extd = [(*series.geom).c2/(*series.geom).gxx,$
          (*series.geom).c2[0]/(*series.geom).gxx[0]]
        c2_interp = (INTERPOL(c2_extd,(*i).z_res+1,/SPLINE))[0:(*i).z_res-1]
         
        FOR x = 0, (*i).x_res_red - 1 DO BEGIN
          x_phys = (x - (*i).x_res_red / 2) * (*i).dx * par.Lref
          Z_pos = Z_interp+c2_interp*x_phys
          R_pos = R_interp+c1_interp*x_phys

          theta_interp = ATAN(Z_interp-(*i).torZ,R_interp-(*i).torR)
          (*i).space_theta[*,x] = ATAN(Z_pos-(*i).torZ,R_pos-(*i).torR)
;          (*i).space_r[*,x] = SQRT((Z_pos)^2+(R_pos-(*i).torR)^2)

          ; in a local simulation using the realistic rhostar,
          ; the radial box length might intersect the origin
          ; and might thus show weird results;
          ; hence, check whether radius vector intersects the origin
          ; and set r to zero if so
          FOR z = 0, (*i).z_res - 1 DO BEGIN
            IF ((ABS(theta_interp[z]-(*i).space_theta[z,x]) GT !PI/4.0) $
            AND (ABS(theta_interp[z]-(*i).space_theta[z,x]) LT 1.5*!PI)) $
            THEN (*i).space_r[z,x] = 0.0 ELSE $
              (*i).space_r[z,x] = SQRT((Z_pos[z]-(*i).torZ)^2+(R_pos[z]-(*i).torR)^2)
          ENDFOR
        ENDFOR

      ENDIF ELSE BEGIN ; nonlocal realistic geometry
        Z_interp = FLTARR((*i).z_res,par.nx0,/NOZERO)
        R_interp = FLTARR((*i).z_res,par.nx0,/NOZERO)
        FOR x = 0, par.nx0 - 1 DO BEGIN
          Z_extd = [REFORM((*series.geom).Z[x,*]),$
            REFORM((*series.geom).Z[x,0])]
          Z_interp[*,x] = (INTERPOL(Z_extd,(*i).z_res+1,/SPLINE))[0:(*i).z_res-1]
          R_extd = [REFORM((*series.geom).R[x,*]),$
            REFORM((*series.geom).R[x,0])]
          R_interp[*,x] = (INTERPOL(R_extd,(*i).z_res+1,/SPLINE))[0:(*i).z_res-1]            
        ENDFOR

        FOR x = 0, (*i).x_res_red - 1 DO BEGIN
          (*i).space_r[*,x] = SQRT((Z_interp[*,x+glob_x_red]-(*i).torZ)^2+$
                                   (R_interp[*,x+glob_x_red]-(*i).torR)^2)
          (*i).space_theta[*,x] = ATAN(Z_interp[*,x+glob_x_red]-(*i).torZ,$
                                       R_interp[*,x+glob_x_red]-(*i).torR)
        ENDFOR         
      ENDELSE
    ENDELSE

    ; calculating y interpolation for cuts
    two_o_z_res = 2.0 / FLOAT((*i).z_res)
    dy_inv = 1.0 / (*i).dy
    yres_o_2 = (*i).y_res / 2

    FOR x = 0, (*i).x_res - 1 DO BEGIN
      IF par.x_local THEN BEGIN
        phi_auxterm = REBIN(REFORM((*i).phi_cut*(*series.geom).C_y,$
          [1,(*i).n_cuts]),[(*i).z_res,(*i).n_cuts])
        qprof_auxterm = !PI * ((*i).r0 + par.shat * $
          (x - (*i).x_res / 2) * (*i).dx)*Cyq0_r0
      ENDIF ELSE BEGIN
        phi_auxterm = REBIN(REFORM((*i).phi_cut*C_y[x],$
          [1,(*i).n_cuts]),[(*i).z_res,(*i).n_cuts])
        qprof_auxterm = !PI * q_prof[x] * C_y[x]

        IF use_sh_met THEN BEGIN
          xc = par.nx0 / 2
          nx = ABS(x-xc) + 1
          q0 = (q_prof[xc] + q_prof[xc-1])*0.5

          IF x GE xc THEN BEGIN
            gxx_x = gxx[xc:x,*]
            gxz_x = gxz[xc:x,*]
            qprof_x = q_prof[xc:x]
            shmet_int = TOTAL(REBIN(par.rhostar*par.minor_r*qprof_x,$
              [nx,(*i).z_res+1])*gxz_x/gxx_x,1)+ $
              !PI*(two_o_z_res*INDGEN((*i).z_res+1)-1.0)*$
              (q_prof[x]-q0)
          ENDIF ELSE BEGIN
            gxx_x = gxx[x:xc-1,*]
            gxz_x = gxz[x:xc-1,*]
            qprof_x = q_prof[x:xc-1]
            shmet_int = - TOTAL(REBIN(par.rhostar*par.minor_r*qprof_x,$
              [nx-1,(*i).z_res+1])*gxz_x/gxx_x,1)+ $
              !PI*(two_o_z_res*INDGEN((*i).z_res+1)-1.0)*$
              (q_prof[x]-q0)
          ENDELSE
          shmet_corr = (C_y[x] * series.dx) * TEMPORARY(shmet_int)
          (*i).y_shift[x] = round(2.0 * !PI / (*i).dy*(C_y[xc]*q0))
          shmet_corr = (TEMPORARY(shmet_corr))[0:(*i).z_res-1]
        ENDIF
      ENDELSE

      slice_y = qprof_auxterm * (two_o_z_res * REBIN(INDGEN((*i).z_res),$
        [(*i).z_res,(*i).n_cuts]) - 1.0) - phi_auxterm
      IF use_sh_met THEN slice_y -= $
        REBIN(shmet_corr,[(*i).z_res,(*i).n_cuts])
      slice_y = TEMPORARY(slice_y) * dy_inv + yres_o_2
      slice_y = ((slice_y MOD (*i).y_res) + (*i).y_res) MOD (*i).y_res

      closest_y = FIX(slice_y)
      dy_lo = slice_y - closest_y
      dy_hi = 1.0 - dy_lo
      closest_y_hi = (((closest_y + 1) MOD (*i).y_res) + $
        (*i).y_res) MOD (*i).y_res
      closest_y_lo = ((closest_y MOD (*i).y_res) + $
        (*i).y_res) MOD (*i).y_res

      (*i).closest_y_cut[0,*,*,x] = closest_y_lo
      (*i).closest_y_cut[1,*,*,x] = closest_y_hi
      (*i).dy_cut[0,*,*,x] = dy_hi
      (*i).dy_cut[1,*,*,x] = dy_lo
    ENDFOR

    ; calculating y interpolation for hulls
    IF (*i).mode NE 0 THEN BEGIN
      xind = [0,(*i).x_res_red-1]+(*i).glob_x_red
      FOR x = 0, 1 DO BEGIN
        IF par.x_local THEN $
          qprof_auxterm = !PI * ((*i).r0 + par.shat * $
          (xind[x] - (*i).x_res / 2) * (*i).dx) * Cyq0_r0 ELSE $
          qprof_auxterm = !PI * q_prof[xind[x]] * C_y[xind[x]]

        FOR p = 0, phi_res - 1 DO BEGIN
          phi_cut_1 = (*i).phi_cut[0]
          phi_auxterm = par.x_local ? $
            (2.0 * !PI * phi_range * p / phi_res + phi_cut_1) * $
            (*series.geom).C_y : $
            ((2.0*!PI*phi_range*p/phi_res+phi_cut_1) MOD (2.0*!PI))*$
            C_y[xind[x]]

          estim_y = (qprof_auxterm * (two_o_z_res * INDGEN((*i).z_res) - $
            1.0) - phi_auxterm) * dy_inv + yres_o_2
          estim_y = ((estim_y MOD (*i).y_res) + (*i).y_res) MOD (*i).y_res

          closest_y_lo = FIX(estim_y)
          dy_lo = estim_y - closest_y_lo
          dy_hi = 1.0 - dy_lo
          closest_y_hi = (closest_y_lo + 1) MOD (*i).y_res

          (*i).closest_y_hull[0,*,p,x] = closest_y_lo
          (*i).closest_y_hull[1,*,p,x] = closest_y_hi
          (*i).dy_hull[0,*,p,x] = dy_hi
          (*i).dy_hull[1,*,p,x] = dy_lo
        ENDFOR
      ENDFOR
    ENDIF
  ENDIF ELSE BEGIN ; full torus geometry for VK
    PRINT, 'calculating coordinates for VK data set'

    run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
    filename_run = run_array[0]
    IF N_ELEMENTS(run_array) GT 1 THEN $
      filename_run += '_' + run_array[N_ELEMENTS(run_array)-1]

    coord_path = gui.out.data_path + 'vis_coord_' + filename_run
    OPENW, coord_lun, coord_path, /GET_LUN

    coord_loc_xyz = $
      INTARR(6,1L*(*i).x_res*(*i).bigM*(*i).y_res*(*i).z_res,/NOZERO)
    n_points = 0L

    xres_o_2 = (*i).x_res / 2
    yres_o_2 = (*i).y_res / 2
    space_range_inv_32767 = 32767.0 / $
      (torR + (*i).r0 + ((*i).x_res / 2) * (*i).dx)

    y_res_full = (*i).y_res * (*i).bigM
    y_arr_full = INDGEN(y_res_full)
    y_arr_full_mod = y_arr_full MOD (*i).y_res

    FOR x = 0, (*i).x_res - 1 DO BEGIN
      in_x = (x - xres_o_2) * (*i).dx
      r = (*i).r0 + in_x
      x_rebin = REBIN([x],y_res_full)

      IF par.x_local THEN BEGIN
        phi_auxterm = 1.0 / (*series.geom).C_y
        qprof_auxterm = ((*i).r0 + in_x * par.shat) / (*series.geom).C_y
      ENDIF ELSE BEGIN
        phi_auxterm = 1.0 / C_y[x]
        qprof_auxterm = q_prof[x]
      ENDELSE

      FOR z = 0, (*i).z_res - 1 DO BEGIN
        theta = !PI * (2.0 * z * zres_inv - 1.0)
        r_auxterm = space_range_inv_32767 * (torR + r * COS(theta))
        r_sin_theta = REBIN([space_range_inv_32767*(r*SIN(theta))],$
          y_res_full)

        in_y = (y_arr_full - yres_o_2) * (*i).dy
        phi = - phi_auxterm * in_y + qprof_auxterm * theta

        coord_loc_xyz[*,n_points+y_arr_full] = TRANSPOSE(REFORM($
          [r_auxterm*COS(phi),r_auxterm*SIN(phi),r_sin_theta,x_rebin,$
          y_arr_full_mod,REBIN([z],y_res_full)],[y_res_full,6]))

        n_points += y_res_full
      ENDFOR
    ENDFOR

    PRINT, 'writing coordinates to file'
    WRITEU, coord_lun, coord_loc_xyz
    FREE_LUN, coord_lun

    data_out_path = gui.out.data_path + 'vis_data_' + filename_run
    OPENW, vk_data_lun_temp, data_out_path, /GET_LUN
    (*i).vk_data_lun = vk_data_lun_temp
  ENDELSE

  fft_format, sxky=use_sxky*(*e.var), sxsy=(use_sxky EQ 0)*(*e.var) ;,doppler_freq=par.omegatorref*par.Lref/par.cref

END

;######################################################################

PRO torus_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).mode LT 3 THEN cut_data = pf_arr([(*i).z_res+((*i).mode EQ 1),$
    (*i).x_res,(*i).n_cuts,gui.out.n_spec_sel])
  IF ((*i).mode EQ 1) OR ((*i).mode EQ 2) THEN hull_data = $
    pf_arr([(*i).z_res+1,(*i).phi_res,2,gui.out.n_spec_sel])

  FOR sp = 0, gui.out.n_spec_sel - 1 DO BEGIN

    IF NOT (*i).use_sxky THEN BEGIN ; should be faster, less memory required(?)
       data = FLOAT((*mom[sp,(*i).var].sxsy)[*,*,*])
;  uncomment to check coordinate lines (correct alignment of cuts and hulls)
;  data[*,par.ny0/2-1:par.ny0/2+1,*] = 10.0 * MAX(ABS(data))
       interp_3d_data, data, (*i).x_res, (*i).y_res, (*i).z_res, (*i).y_shift
    ENDIF ELSE BEGIN ; more accurate
       data = (*mom[sp,(*i).var].sxky)[*,*,*]
       interp_3d_data, data, (*i).x_res, (*i).y_res, (*i).z_res, $
            (par.x_local?(*i).y_shift/(*i).y_res:(*i).bigM*(*series.geom).q),$
             infft=[0,1,0]
    ENDELSE

    ; storing data
    IF (*i).mode LT 3 THEN BEGIN ; no VK
      ; storing cut data
      xind = REBIN(REFORM(INDGEN((*i).x_res),[1,(*i).x_res]),$
        [(*i).z_res,(*i).x_res])
      zind = REBIN(INDGEN((*i).z_res),[(*i).z_res,(*i).x_res])

      FOR c = 0, (*i).n_cuts - 1 DO BEGIN
        cyc0 = REFORM((*i).closest_y_cut[0,*,c,*])
        cyc1 = REFORM((*i).closest_y_cut[1,*,c,*])
        
        cut_data[0,0,c,sp] = REFORM((*i).dy_cut[0,*,c,*]*$
          data[xind,cyc0,zind]+(*i).dy_cut[1,*,c,*]*data[xind,cyc1,zind])
      ENDFOR

      IF (*i).mode GE 1 THEN BEGIN
        ; storing hull data (in z-phi)
        zind = REBIN(INDGEN((*i).z_res),[(*i).z_res,(*i).phi_res])
        FOR x = 0, 1 DO BEGIN
          xind = REBIN([(*i).glob_x_red+x*((*i).x_res_red-1)],$
            [(*i).z_res,(*i).phi_res])

          hull_data[0,0,x,sp] = REFORM((*i).dy_hull[0,*,*,x]*$
            data[xind,REFORM((*i).closest_y_hull[0,*,*,x]),zind]+$
            (*i).dy_hull[1,*,*,x]*$
            data[xind,REFORM((*i).closest_y_hull[1,*,*,x]),zind])
        ENDFOR

;        IF (*i).mode EQ 1 THEN $
;          cut_data[(*i).z_res,*,sp] = cut_data[0,*,sp]
        hull_data[(*i).z_res,*,*,sp] = hull_data[0,*,*,sp]
      ENDIF
    ENDIF ELSE BEGIN ; VK
      IF (*i).norm_local THEN BEGIN
        IF data_min NE data_max THEN data_range_inv = 1.0 / $
          (data_max - data_min) ELSE data_range_inv = 1.0

        data_norm = UINT(ROUND(((data[*,*,*]-data_min)*$
          data_range_inv)*65535.0)<65535>0)
        WRITEU, (*i).vk_data_lun, data_norm
      ENDIF ELSE (*i).vk_data_id = time_avg((*i).vk_data_id,data,$
                                     fwd=gui.out.res_steps)
    ENDELSE
    
    IF NOT (*i).norm_local THEN BEGIN
      data_min = MIN((*i).mode GT 0 ? data : cut_data[*,*,*,sp],$
        MAX=data_max)
      (*i).data_minmax[*,sp] = [(*i).data_minmax[0,sp] < data_min,$
        (*i).data_minmax[1,sp] > data_max]
    ENDIF
  ENDFOR

  IF (*i).mode LT 3 THEN BEGIN ; no VK
    (*i).cut_data_id = time_avg((*i).cut_data_id,cut_data,mom_time,$
                                  fwd=gui.out.res_steps,rms=(*i).rms)
    IF (*i).mode GT 0 THEN (*i).hull_data_id = $
       time_avg((*i).hull_data_id,hull_data,$
                mom_time,fwd=gui.out.res_steps,rms=(*i).rms)
  ENDIF

END

;######################################################################

PRO torus_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  coltable = 33 ; rainbow colors
;  coltable = 3 ; orange plasma colors
;  coltable = 99 ; red-yellow
;  coltable = 0 ; black-white (use with contourfill=0?)

  csize = 1.4
  cthick = 5.0

  n_res_steps = gui.out.res_steps * (((*i).tintrpfac>1)*$
                series.step_count - 1) + 1

  eps = (n_res_steps EQ 1) AND ((*i).n_cuts EQ 1)

  IF (*i).mode EQ 3 THEN BEGIN
    vk_data = time_avg((*i).vk_data_id,/get,fwd=gui.out.res_steps,$
                      tarr=time,intrp=(*i).tintrpfac)
    vk_data = REFORM(vk_data,[(*i).x_res, (*i).y_res, (*i).z_res,$
              gui.out.n_spec_sel,n_res_steps],/OVERWRITE)

    PRINT, '  writing VK data file'
    FOR n = 0, n_res_steps - 1 DO BEGIN
      IF (*i).data_minmax[1,0] - (*i).data_minmax[0,0] EQ 0.0 THEN $
        (*i).data_minmax[1,0] = (*i).data_minmax[0,0] + 1.0
      data_range_inv = $
        1.0 / ((*i).data_minmax[1,0] - (*i).data_minmax[0,0])

      data_norm = UINT(ROUND(((vk_data[*,*,*,n]-$
        (*i).data_minmax[0,0])*data_range_inv)*65535.0)<65535>0)
      WRITEU, (*i).vk_data_lun, data_norm
    ENDFOR
    FREE_LUN, (*i).vk_data_lun

    run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
    filename_run = run_array[0]
    IF N_ELEMENTS(run_array) GT 1 THEN $
      filename_run += '_' + run_array[N_ELEMENTS(run_array)-1]

    vis_par_path = gui.out.data_path + 'vis_par_' + filename_run
    OPENW, vis_par_lun, vis_par_path, /GET_LUN
    PRINTF, vis_par_lun, $
      'coord:  signed 2B binary; x y z resolution (space):', $
      rm0es([(*i).x_res,(*i).bigM*(*i).y_res,(*i).z_res])
;      rm0es([(*i).x_res,n_points,1])
    PRINTF, vis_par_lun, $
      'data: unsigned 2B binary; x y z resolution (local):', $
      rm0es([(*i).x_res,(*i).y_res,(*i).z_res])
    PRINTF, vis_par_lun, 'time = '+(gui.out.res_steps ? $
           ', time = '+rm0es(time[n]) : 't avg.')
;    y_shift_base = 2.0 * !PI * par.shat * (*i).dx / (*i).dy
;    PRINTF, vis_par_lun, 'y_shift: ', rm0es(y_shift_base)
    FREE_LUN, vis_par_lun
  ENDIF ELSE BEGIN
    cut_data = time_avg((*i).cut_data_id,/avg,fwd=gui.out.res_steps,$
      tarr=time,intrp=(*i).tintrpfac,rms=(*i).rms)
    cut_data = REFORM(cut_data,[(*i).z_res+((*i).mode EQ 1),$
      (*i).x_res,(*i).n_cuts,gui.out.n_spec_sel,n_res_steps],/OVERWRITE)

    IF (*i).mode GT 0 THEN BEGIN
      hull_data = time_avg((*i).hull_data_id,/avg,rms=(*i).rms,$
        fwd=gui.out.res_steps,tarr=time,intrp=(*i).tintrpfac)
      hull_data = REFORM(hull_data,[(*i).z_res+1,(*i).phi_res,2,$
        gui.out.n_spec_sel,n_res_steps],/OVERWRITE)
    ENDIF
    IF (*i).mode EQ 1 THEN BEGIN ; partial torus
      n_points_hull = ((*i).z_res + 1L) * (*i).phi_res
      n_poly_hull = (*i).z_res * ((*i).phi_res - 1L)
      n_points_cut = ((*i).z_res + 1L) * (*i).x_res_red
      n_poly_cut = (*i).z_res * ((*i).x_res_red - 1L)

      vert_hull = FLTARR(3,n_points_hull,2,/NOZERO)
      poly_hull = LONARR(5L*n_poly_hull,/NOZERO)
      vert_cut = FLTARR(3,n_points_cut,/NOZERO)
      poly_cut = LONARR(5L*n_poly_cut,/NOZERO)

      ; create hull vertices and polygons
      two_o_z_res = 2.0 / FLOAT((*i).z_res)
      FOR x = 0, 1 DO BEGIN
        xmod = x * ((*i).x_res_red - 1)

        FOR p = 0, (*i).phi_res - 1 DO BEGIN
          phi = 2.0 * !PI * ((*i).phi_res - p) * $
            (*i).phi_range / (*i).phi_res

          IF par.x_local THEN BEGIN
            vert_hull[0,((*i).z_res+1L)*p+INDGEN((*i).z_res),x] = $
              ((*i).torR + (*i).space_r[xmod] * $
              COS((*i).space_theta[*])) * COS(phi)
            vert_hull[1,((*i).z_res+1L)*p+INDGEN((*i).z_res),x] = $
              - (*i).space_r[xmod] * SIN((*i).space_theta[*])
            vert_hull[2,((*i).z_res+1L)*p+INDGEN((*i).z_res),x] = $
              ((*i).torR + (*i).space_r[xmod] * $
              COS((*i).space_theta[*])) * SIN(phi)
          ENDIF ELSE BEGIN
            vert_hull[0,((*i).z_res+1L)*p+INDGEN((*i).z_res),x] = $
              ((*i).torR + (*i).space_r[*,xmod] * $
              COS((*i).space_theta[*,x])) * COS(phi)
            vert_hull[1,((*i).z_res+1L)*p+INDGEN((*i).z_res),x] = $
              - (*i).space_r[*,xmod] * SIN((*i).space_theta[*,x])
            vert_hull[2,((*i).z_res+1L)*p+INDGEN((*i).z_res),x] = $
              ((*i).torR + (*i).space_r[*,xmod] * $
              COS((*i).space_theta[*,x])) * SIN(phi)
          ENDELSE

          vert_hull[*,(*i).z_res+((*i).z_res+1L)*p,x] = $
            vert_hull[*,((*i).z_res+1L)*p,x]
        ENDFOR
      ENDFOR

      point = 0L
      FOR p = 0, (*i).phi_res - 2 DO BEGIN
        FOR z = 0, (*i).z_res - 1 DO BEGIN
          zy_ind = 5L * (z + 1L * (*i).z_res * p)
          poly_hull[zy_ind:zy_ind+4L] = $
            [4,point,point+1L,point+(*i).z_res+2L,point+(*i).z_res+1L]
          point = point + 1L
        ENDFOR
        point = point + 1L
      ENDFOR

      ; create cut vertices and polygons
      FOR x = 0, (*i).x_res_red - 1 DO BEGIN
        IF par.x_local THEN BEGIN
          vert_cut[0,((*i).z_res+1L)*x+INDGEN((*i).z_res)] = $
            (*i).torR+(*i).space_r[x]*COS((*i).space_theta[*])
          vert_cut[1,((*i).z_res+1L)*x+INDGEN((*i).z_res)] = $
            -(*i).space_r[x]*SIN((*i).space_theta[*])
          vert_cut[2,((*i).z_res+1L)*x+INDGEN((*i).z_res)] = 0.0
        ENDIF ELSE BEGIN
          vert_cut[0,((*i).z_res+1L)*x+INDGEN((*i).z_res)] = $
            (*i).torR+(*i).space_r[*,x]*COS((*i).space_theta[*,x])
          vert_cut[1,((*i).z_res+1L)*x+INDGEN((*i).z_res)] = $
            -(*i).space_r[*,x]*SIN((*i).space_theta[*,x])
          vert_cut[2,((*i).z_res+1L)*x+INDGEN((*i).z_res)] = 0.0
        ENDELSE

        vert_cut[*,(*i).z_res+((*i).z_res+1L)*x] = $
          vert_cut[*,((*i).z_res+1L)*x]
      ENDFOR

      point = 0L
      FOR x = 0, (*i).x_res_red - 2 DO BEGIN
        FOR z = 0, (*i).z_res - 1 DO BEGIN
          zy_ind = 5L * (z + 1L * (*i).z_res * x)
          poly_cut[zy_ind:zy_ind+4L] = $
            [4,point,point+1L,point+(*i).z_res+2L,point+(*i).z_res+1L]
          point = point + 1L
        ENDFOR
        point = point + 1L
      ENDFOR

      T3D, /RESET
      T3D, ROTATE=[0.0,-22.5,0.0]
      vert_cut = VERT_T3D(vert_cut)

      T3D, /RESET
      T3D, ROTATE=[20.0,0.0,0.0]
      T3D, ROTATE=[0.0,0.0,0.0]
      T3D, ROTATE=[0.0,0.0,0.0]
      T3D, TRANSLATE=[-0.48,0.58,0.0]
      FOR x = 0, 1 DO vert_hull[*,*,x] = VERT_T3D(vert_hull[*,*,x])
      vert_cut = VERT_T3D(vert_cut)

      cut_data_norm = BYTARR((*i).z_res+1,(*i).x_res_red,/NOZERO)
      hull_data_norm = BYTARR((*i).z_res+1,(*i).phi_res,2,/NOZERO)
    ENDIF ; mode EQ 1 (partial torus with IDL) 

    IF (*i).png_out THEN BEGIN
      ; note: typically, the image will have N x N pixels, but for
      ; general geometry, it may be adapted to have N x M, M < N pixels
      png_resolution = 800
      n_levels = 255
  
      xmarg = !X.MARGIN
      !X.MARGIN = 0
      ymarg = !Y.MARGIN
      !Y.MARGIN = 0

      bg_color = !P.BACKGROUND
      !P.BACKGROUND = 16777215L
    ENDIF ELSE n_levels = 30

    IF (NOT (*i).contourfill) THEN n_levels = 15

    IF (*i).mode NE 1 THEN BEGIN ; not partial torus
      IF (*i).png_out THEN BEGIN ; prepare direct graphics output
        ; Note that this is done to speed up the plotting process.
        ; While technically, one might use this technique for the
        ; postscript output, as well, it produces pixel-based
        ; graphics. Additionally, the speedup is larger when the
        ; full color range is required, as is the case for png output.
        IF (*i).geom_type NE 2 THEN BEGIN
          n_pixels = png_resolution
          png_resolution = [n_pixels,n_pixels]

          xy_arr = 2.0 * MAX((*i).space_r) * $
            (FINDGEN(n_pixels) / n_pixels - 0.5)
          img_theta = ATAN(REBIN(REFORM(xy_arr,[1,n_pixels]),$
            [n_pixels,n_pixels]),REBIN(xy_arr,[n_pixels,n_pixels]))
          img_r = SQRT(REBIN(xy_arr^2,[n_pixels,n_pixels])+$
            REBIN(REFORM(xy_arr^2,[1,n_pixels]),[n_pixels,n_pixels]))

          img_rselect = WHERE((img_r GE MIN((*i).space_r)) AND $
            (img_r LE MAX((*i).space_r)))
          img_theta = img_theta[img_rselect]
          img_r = img_r[img_rselect]

          n_img_sel = N_ELEMENTS(img_rselect)

          IF (*i).geom_type EQ 0 THEN BEGIN ; s-alpha
            ; linear interpolation to get pixel value
            theta_inds_0 = VALUE_LOCATE((*i).space_theta,img_theta)
            theta_inds_1 = (theta_inds_0 + 1) MOD (*i).z_res
            r_inds_0 = VALUE_LOCATE((*i).space_r,img_r)
            r_inds_1 = (r_inds_0 + 1) < ((*i).x_res - 1)

            dr0 = img_r - (*i).space_r[r_inds_0]
            dr1 = (*i).space_r[r_inds_1] - img_r
            dt0 = img_theta - (*i).space_theta[theta_inds_0]
            dt1 = (*i).space_theta[theta_inds_1] - img_theta

            ; correction at the boundary
            dt1_lt_negpi = WHERE(dt1 LT -!PI)
            IF dt1_lt_negpi[0] NE -1 THEN dt1[dt1_lt_negpi] += 2.0 * !PI

            drdt_temp = (dr0 + dr1) * (dt0 + dt1)
            drdt_zeroinds = WHERE(drdt_temp EQ 0)
            IF drdt_zeroinds[0] NE -1 THEN drdt_temp[drdt_zeroinds] = 1.0
            drdt_norm =  1.0 / TEMPORARY(drdt_temp)
          ENDIF ELSE BEGIN
            r_inds_0 = LONARR(n_img_sel,/NOZERO)
            r_inds_1 = LONARR(n_img_sel,/NOZERO)
            theta_inds_0_0 = LONARR(n_img_sel,/NOZERO)
            theta_inds_0_1 = LONARR(n_img_sel,/NOZERO)
            theta_inds_1_0 = LONARR(n_img_sel,/NOZERO)
            theta_inds_1_1 = LONARR(n_img_sel,/NOZERO)
            dr0 = FLTARR(n_img_sel,/NOZERO)
            dt0 = FLTARR(n_img_sel,/NOZERO)
            dt1 = FLTARR(n_img_sel,/NOZERO)

            FOR ir = 0L, n_img_sel - 1L DO BEGIN
              r_inds_0[ir] = VALUE_LOCATE((*i).space_r[0,*],img_r[ir])
              r_inds_1[ir] = (r_inds_0[ir] + 1) < ((*i).x_res_red - 1)
              dr0[ir] = (img_r[ir] - (*i).space_r[0,r_inds_0[ir]])

              dr = (*i).space_r[0,r_inds_1[ir]] - (*i).space_r[0,r_inds_0[ir]]
              IF dr NE 0 THEN dr0[ir] /= dr

              theta_inds_0_0[ir] = $
                VALUE_LOCATE((*i).space_theta[*,r_inds_0[ir]],img_theta[ir])
              theta_inds_0_1[ir] = (theta_inds_0_0[ir] + 1) MOD (*i).z_res
              dt = (*i).space_theta[theta_inds_0_1[ir],r_inds_0[ir]] - $
                (*i).space_theta[theta_inds_0_0[ir],r_inds_0[ir]]
              IF dt LT -!PI THEN dt += 2.0 * !PI
              dt0[ir] = (img_theta[ir] - $
                (*i).space_theta[theta_inds_0_0[ir],r_inds_0[ir]]) / dt < 1.0

              theta_inds_1_0[ir] = $
                VALUE_LOCATE((*i).space_theta[*,r_inds_1[ir]],img_theta[ir])
              theta_inds_1_1[ir] = (theta_inds_1_0[ir] + 1) MOD (*i).z_res
              dt = (*i).space_theta[theta_inds_1_1[ir],r_inds_1[ir]] - $
                (*i).space_theta[theta_inds_1_0[ir],r_inds_1[ir]]
              IF dt LT -!PI THEN dt += 2.0 * !PI
              dt1[ir] = (img_theta[ir] - $
                (*i).space_theta[theta_inds_1_0[ir],r_inds_1[ir]]) / dt < 1.0
            ENDFOR
          ENDELSE
        ENDIF ELSE BEGIN ; geom_type = 2 (tracer, chease, ...)
          plot_height_1 = MAX((*i).space_r*SIN((*i).space_theta),$
            MIN=plot_height_0)
          plot_width_1 = MAX((*i).space_r*COS((*i).space_theta),$
            MIN=plot_width_0)
          plot_height = plot_height_1 - plot_height_0
          plot_width = plot_width_1 - plot_width_0

          IF plot_height GE plot_width THEN BEGIN
            png_res_x = ROUND(png_resolution*plot_width/plot_height)
            png_res_y = png_resolution
          ENDIF ELSE BEGIN
            png_res_x = png_resolution
            png_res_y = ROUND(png_resolution*plot_height/plot_width)
          ENDELSE
          png_resolution = [png_res_x,png_res_y]

          x_arr = FINDGEN(png_res_x) * $
            (plot_width / png_res_x) + plot_width_0
          y_arr = FINDGEN(png_res_y) * $
            (plot_height / png_res_y) + plot_height_0

          img_theta = ATAN(REBIN(REFORM(y_arr,[1,png_res_y]),$
            [png_res_x,png_res_y]),REBIN(x_arr,[png_res_x,png_res_y]))
          img_r = SQRT(REBIN(x_arr^2,[png_res_x,png_res_y])+$
            REBIN(REFORM(y_arr^2,[1,png_res_y]),[png_res_x,png_res_y]))

          inner_border_r = (*i).space_r[*,0]
          outer_border_r = (*i).space_r[*,(*i).x_res_red-1]
          inner_border_theta = (*i).space_theta[*,0]
          outer_border_theta = (*i).space_theta[*,(*i).x_res_red-1]

          FOR z = 0, (*i).z_res - 1 DO BEGIN
            z_plus_1 = (z + 1) MOD (*i).z_res

            use_boundary_inner = ABS(inner_border_theta[z]-$
               inner_border_theta[z_plus_1]) GT !PI
;            use_boundary_outer = outer_border_theta[z] GT $
;                                 outer_border_theta[z_plus_1]
            use_boundary_outer = ABS(outer_border_theta[z]-$
                outer_border_theta[z_plus_1]) GT !PI

            theta_range_inner = NOT use_boundary_inner ? $
              WHERE((img_theta GE inner_border_theta[z]) AND $
              (img_theta LT inner_border_theta[z_plus_1])) : $
              WHERE(((img_theta GE inner_border_theta[z]) AND $
              (img_theta LT inner_border_theta[z_plus_1] + 2 * !PI)) OR $
              ((img_theta GE inner_border_theta[z] - 2 * !PI) AND $
              (img_theta LT inner_border_theta[z_plus_1])))
            theta_range_outer = NOT use_boundary_outer ? $
              WHERE((img_theta GE outer_border_theta[z]) AND $
              (img_theta LT outer_border_theta[z_plus_1])) : $
              WHERE(((img_theta GE outer_border_theta[z]) AND $
              (img_theta LT outer_border_theta[z_plus_1] + 2 * !PI)) OR $
              ((img_theta GE outer_border_theta[z] - 2 * !PI) AND $
              (img_theta LT outer_border_theta[z_plus_1])))

            IF theta_range_inner[0] NE -1 THEN BEGIN
              img_r_inthetarange = img_r[theta_range_inner]
              img_theta_inthetarange = img_theta[theta_range_inner]

              img_rselect_z_inner = WHERE((img_r_inthetarange - $
                ((inner_border_theta[z_plus_1] - $
                img_theta_inthetarange) * inner_border_r[z] + $
                (img_theta_inthetarange - inner_border_theta[z]) * $
                inner_border_r[z_plus_1]) / $
                (inner_border_theta[z_plus_1] - inner_border_theta[z])) GT 0.0)
              IF (img_rselect_z_inner[0] NE -1) THEN $
                img_rselect_inner = N_ELEMENTS(img_rselect_inner) LT 1 ? $
                theta_range_inner[img_rselect_z_inner] : $
                [img_rselect_inner,theta_range_inner[img_rselect_z_inner]]
            ENDIF
            IF theta_range_outer[0] NE -1 THEN BEGIN
              img_r_inthetarange = img_r[theta_range_outer]
              img_theta_inthetarange = img_theta[theta_range_outer]

              img_rselect_z_outer = WHERE((img_r_inthetarange - $
                ((outer_border_theta[z_plus_1] - $
                img_theta_inthetarange) * outer_border_r[z] + $
                (img_theta_inthetarange - outer_border_theta[z]) * $
                outer_border_r[z_plus_1]) / $
                (outer_border_theta[z_plus_1] - outer_border_theta[z])) LT 0.0)
              IF (img_rselect_z_outer[0] NE -1) THEN $
                img_rselect_outer = N_ELEMENTS(img_rselect_outer) LT 1 ? $
                theta_range_outer[img_rselect_z_outer] : $
                [img_rselect_outer,theta_range_outer[img_rselect_z_outer]]
            ENDIF
          ENDFOR

          min_inner = MIN(img_rselect_inner,MAX=max_inner)
          min_outer = MIN(img_rselect_outer,MAX=max_outer)
          min_all = min_inner > min_outer
          max_all = max_inner < max_outer

          img_rselect = WHERE((HISTOGRAM(img_rselect_inner,$
            MIN=min_all,MAX=max_all) NE 0) AND $
            (HISTOGRAM(img_rselect_outer,$
            MIN=min_all,MAX=max_all) NE 0)) + min_all
        ENDELSE
      ENDIF ELSE BEGIN ; ps output, triangulation for contour plot
        IF (*i).geom_type EQ 0 THEN BEGIN ; s-alpha
          ; central white circle
          theta_bg = (*i).space_theta
          ; plot data geometry
          r_mod = REBIN(REFORM((*i).space_r,[1,(*i).x_res]),$
            [(*i).z_res,(*i).x_res])
          theta_mod = REBIN(REFORM((*i).space_theta,[(*i).z_res,1]),$
            [(*i).z_res,(*i).x_res])
          plot_x = r_mod * COS(theta_mod)
          plot_y = r_mod * SIN(theta_mod)
        ENDIF ELSE BEGIN
          ; central white circle
          theta_bg = REFORM((*i).space_theta[*,0])
          IF (*i).geom_type EQ 2 THEN r_bg = REFORM((*i).space_r[*,0])

          ; plot data geometry
          plot_x = (*i).torR + (*i).space_r * COS((*i).space_theta)
          plot_y = (*i).torZ + (*i).space_r * SIN((*i).space_theta)
        ENDELSE

        ; prepare central white circle
        lev_bg = [0,1]
        dat_bg = REBIN(REFORM([0.0,1.0],[1,2]),[(*i).z_res,2])
        IF (*i).geom_type NE 2 THEN BEGIN
          r_bg = [0.0,(*i).r0-((*i).x_res_red+2)/2*(*i).dx]
          r_bg_mod = REBIN(REFORM(r_bg,[1,2]),[(*i).z_res,2])
          IF (par.Lref GT 0.0) THEN r_bg_mod *= par.Lref
        ENDIF ELSE BEGIN
          r_bg_mod = FLTARR((*i).z_res,2)
          r_bg_mod[0,1] = r_bg
        ENDELSE
        theta_bg_mod = REBIN(REFORM(theta_bg,[(*i).z_res,1]),[(*i).z_res,2])
        plot_x_bg = ((*i).geom_type EQ 0 ? 0.0 : (*i).torR) + $
          r_bg_mod * COS(theta_bg_mod)
        plot_y_bg = ((*i).geom_type EQ 0 ? 0.0 : (*i).torZ) + $
          r_bg_mod * SIN(theta_bg_mod)
  
        TRIANGULATE, plot_x_bg, plot_y_bg, triangulation_bg, TOLERANCE=0
        TRIANGULATE, plot_x, plot_y, triangulation, TOLERANCE=0
      ENDELSE
    ENDIF

    FOR sp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      sp_id = (*gui.out.spec_select)[sp]

      IF (*i).png_out THEN BEGIN
        set_output, diag, sp_id, coltable=coltable, /reset

        IF N_ELEMENTS(png_resolution) EQ 1 THEN $
          png_resolution = [1,1] * png_resolution

        SET_PLOT, 'Z'
        DEVICE, SET_RESOLUTION=png_resolution
        image = BYTARR(4,png_resolution[0],png_resolution[1])
        frame = INTARR(png_resolution[0],png_resolution[1])
        framenr = 0
        n_png_files = series.step_count * ((*i).n_cuts + 2 * ((*i).mode EQ 2))

        run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
        filename_run = spec[sp_id].name + '_' + run_array[0]
        IF N_ELEMENTS(run_array) GT 1 THEN $
          filename_run += '_' + run_array[N_ELEMENTS(run_array)-1]

        png_name_base = gui.out.out_path + (*diag).name + filename_run + '_'
        IF (*i).mode NE 0 THEN png_name = png_name_base + '[cut/hull]_' ELSE $
          png_name = png_name_base
        IF (*i).n_cuts GT 1 THEN png_name += '[cut]_'
        framenr_str = '[frame]'
        png_name += framenr_str + '.png'

        PRINT, 'writing ' + rm0es(n_png_files) + ' png file' + $
          (n_png_files EQ 1 ? ':' : 's:')
        PRINT, png_name

        step_digits = FIX(ALOG10(series.step_count)) + 1
      ENDIF ELSE set_output, diag, sp_id, ps=1-eps, $
	eps=eps, coltable=coltable, xsize=19 ;, ysize = 20

      IF (*i).data_minmax[1,sp] - (*i).data_minmax[0,sp] EQ 0.0 THEN $
        (*i).data_minmax[1,sp] = (*i).data_minmax[0,sp] + 1.0

      IF NOT (*i).norm_local THEN BEGIN
        IF (*i).use_man_minmax THEN $
          (*i).data_minmax[*,sp] = (*i).man_minmax

        IF (*i).show_minmax THEN PRINT, (*diag).name + $
          ': global min/max for species ' + spec[sp_id].name + ': ' + $
          rm0es((*i).data_minmax[0,sp]) + '/' + rm0es((*i).data_minmax[1,sp])
      ENDIF

      FOR n = 0, n_res_steps - 1 DO BEGIN ; reading data and plotting
        IF (*i).png_out THEN IF $
          N_ELEMENTS(image) NE 1L * png_resolution[0] * png_resolution[1] $
          THEN image = BYTARR(4,png_resolution[0],png_resolution[1])

        lev_col = contour_levels((*i).norm_local ? $
          cut_data[*,*,*,sp,n] : (*i).data_minmax[*,sp],n_levels,$
          no_fix_zero=(NOT (*i).fix_color_zero))

        FOR c = 0, (*i).n_cuts - 1 DO BEGIN
          IF (*i).mode NE 1 THEN BEGIN
            IF (*i).png_out THEN BEGIN
              img = BYTARR(png_resolution)
              cdat = REFORM(cut_data[*,(*i).glob_x_red:(*i).glob_x_red+$
                     (*i).x_res_red-1,c,sp,n],[(*i).z_res,(*i).x_res_red])
              IF (*i).geom_type EQ 0 THEN BEGIN
                img_unnorm = (dr1 * (dt1 * cdat[theta_inds_0,r_inds_0] + $
                  dt0 * cdat[theta_inds_1,r_inds_0]) + $
                  dr0 * (dt1 * cdat[theta_inds_0,r_inds_1] + $
                  dt0 * cdat[theta_inds_1,r_inds_1])) * drdt_norm
              ENDIF ELSE IF (*i).geom_type EQ 1 THEN BEGIN
                img_unnorm = FLTARR(n_img_sel,/NOZERO)
                FOR ip = 0L, n_img_sel - 1 DO BEGIN
                  img_unnorm[ip] = (1.0 - dr0[ip]) * $
                    ((1.0 - dt0[ip]) * cdat[theta_inds_0_0[ip],r_inds_0[ip]] + $
                    dt0[ip] * cdat[theta_inds_0_1[ip],r_inds_0[ip]]) + $
                    dr0[ip] * ((1.0 - dt1[ip]) * cdat[theta_inds_1_0[ip],r_inds_1[ip]] + $
                    dt1[ip] * cdat[theta_inds_1_1[ip],r_inds_1[ip]])
                ENDFOR
              ENDIF ELSE BEGIN
                plot_x = (*i).space_r * COS((*i).space_theta)
                plot_y = (*i).space_r * SIN((*i).space_theta)

                TRIANGULATE, plot_x, plot_y, triangulation, boundary, TOLERANCE=0
                ; known bug: EXTRAPOLATE sometimes fails!
                ; switching off this option however yields 
                ; an ill defined (constant?) outer boundary line
                img_unnorm = (TRIGRID(plot_x,plot_y,cdat,$
                  triangulation, EXTRAPOLATE=boundary,$
                  NX=png_resolution[0],NY=png_resolution[1],$
                  XOUT=xarr,YOUT=yarr))[img_rselect]
              ENDELSE
;              img_min = MIN(img_unnorm,MAX=img_max)
;              IF img_min LT (*i).data_minmax[0,sp] THEN $
;                PRINT, framenr, ' img min: ', img_min, (*i).data_minmax[0,sp]
;              IF img_max GT (*i).data_minmax[1,sp] THEN $
;                PRINT, framenr, ' img max: ', img_max, (*i).data_minmax[1,sp]
;              img_min = img_min < (*i).data_minmax[0,sp]
;              img_max = img_max > (*i).data_minmax[1,sp]

              img[img_rselect] = contour_levels(img_unnorm,n_levels,$
                min=(*i).data_minmax[0,sp],max=(*i).data_minmax[1,sp],/img,$
                no_fix_zero=(NOT (*i).fix_color_zero))

              TV, img
            ENDIF ELSE BEGIN ; ps output for modes 0, 2
              IF par.x_local THEN BEGIN
                CONTOUR, cut_data[*,*,c,sp,n], plot_x, plot_y, $
                  TRIANGULATION=triangulation, C_COLORS=lev_col[*,1], $
                  LEVELS=lev_col[*,0], FILL=(*i).contourfill, /ISOTROPIC, $
                  XSTYLE=1+4*(*i).png_out, YSTYLE=1+4*(*i).png_out;, $
                ; display axis by uncommenting the following lines
;                  XMARGIN=[8,1], YMARGIN=[3,1], $
;                  XTITLE='R / m', YTITLE='Z / m', CHARTHICK=cthick, $
;                  CHARSIZE=csize

                set_output, diag, sp_id, dat=cut_data[*,*,c,sp,n], $
                  append=(c+n GT 0), commentline='phi = '+$
                  rm0es((*i).phi_cut[c]/!PI)+((*i).phi_cut[c] EQ 0 ? '' : $
                  ' pi')+', t = '+(n_res_steps GT 1 ? rm0es(time[n]) : $
                  rm0es(gui.out.start_t) + '-' + rm0es(gui.out.end_t))
              ENDIF ELSE BEGIN
                CONTOUR, cut_data[*,(*i).glob_x_red:(*i).glob_x_red+$
                  (*i).x_res_red-1,c,sp,n], plot_x, plot_y, $
                  TRIANGULATION=triangulation, C_COLORS=lev_col[*,1], $
                  LEVELS=lev_col[*,0], FILL=(*i).contourfill, /ISOTROPIC, $
                  XSTYLE=1+4*(*i).png_out, YSTYLE=1+4*(*i).png_out ;,$
                ; display axis by uncommenting the following lines
;                  XMARGIN=[8,1], YMARGIN=[3,2], $
;                  XTITLE='R / m', YTITLE='Z / m', CHARTHICK=cthick, $
;                  CHARSIZE=csize
                
                set_output, diag, sp_id, $
                  dat=cut_data[*,(*i).glob_x_red:(*i).glob_x_red+$
                  (*i).x_res_red-1,c,sp,n], append=(c+n GT 0)
              ENDELSE

              ; overplot central white circle
              LOADCT, 0, FILE='internal/colortable.tbl'
              col_bg = !P.BACKGROUND * [1,1]
              CONTOUR, dat_bg, plot_x_bg, plot_y_bg, $
                TRIANGULATION=triangulation_bg, C_COLORS=col_bg, $
                LEVELS=lev_bg, /FILL, /OVERPLOT
              LOADCT, coltable, FILE='internal/colortable.tbl'

              ;uncomment to plot colorbar
;              plot_colorbar, lev_col, position=[0.925,0.975], $
;                 charsize=0.75*csize, orientation=1
            ENDELSE
          ENDIF ELSE BEGIN ; mode EQ 1 (partial torus)
            cut_data_norm[0,0] = $
              ROUND((cut_data[*,(*i).glob_x_red:(*i).glob_x_red+$
              (*i).x_res_red-1,sp,n]-(*i).data_minmax[0,sp])/$
              ((*i).data_minmax[1,sp]-(*i).data_minmax[0,sp])*255)
            hull_data_norm[0,0,0] = $
              ROUND((hull_data[*,*,*,sp,n]-(*i).data_minmax[0,sp])/$
              ((*i).data_minmax[1,sp]-(*i).data_minmax[0,sp])*255)

            PLOT, [0.2,0.8], [0.2,0.8], /NODATA, XSTYLE=5, YSTYLE=5, COLOR=0
            SET_SHADING, LIGHT=[0.0,0.0,2.0], REJECT=0

            FOR x = 0, 1 DO dummy = POLYSHADE(vert_hull[*,*,x],$
              poly_hull,SHADE=hull_data_norm[*,*,x],/DATA)
            dummy = POLYSHADE(vert_cut,poly_cut,SHADE=cut_data_norm,/DATA)
          ENDELSE

          IF (*i).png_out THEN BEGIN
            frame = REVERSE(TVRD(),2)
            TVLCT, red, green, blue, /GET

            ; set rgb color channels
            image[0,*,*] = red[frame[*,*]]
            image[1,*,*] = green[frame[*,*]]
            image[2,*,*] = blue[frame[*,*]]

            ; set alpha channel (255B = opaque, 0B = transparent)
            image[3,*,*] = 255B ; opaque as default
            IF (*i).transparent THEN BEGIN
              tind = WHERE(frame EQ 0) ; find black pixels
              IF tind[0] NE -1 THEN BEGIN
                tind = TEMPORARY(tind) * 4 + 3
                image[tind] = 0B ; set to transparent
              ENDIF
            ENDIF

            run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
            filename_run = spec[sp_id].name + '_' + run_array[0]
            IF N_ELEMENTS(run_array) GT 1 THEN $
              filename_run += '_' + run_array[N_ELEMENTS(run_array)-1]

            framenr_str = rm0es(framenr)
            WHILE STRLEN(framenr_str) LT step_digits DO $
              framenr_str = '0' + framenr_str
            png_name_base = $
              gui.out.out_path + (*diag).name + filename_run + '_'
            IF (*i).mode NE 0 THEN png_name = png_name_base + 'cut' ELSE $
              png_name = png_name_base
            IF (*i).n_cuts GT 1 THEN BEGIN
              cut_str = rm0es(c)
              png_name += cut_str + '_'
            ENDIF
            png_name += framenr_str + '.png'

            WRITE_PNG, png_name, image, /ORDER 
            ; the TRANSPARENT keyword does not work as desired,
            ; hence the forth dimension of image is used for the alpha
            ; channel (255B = opaque, 0B = transparent)

            IF N_ELEMENTS(n_pixels) NE 1 THEN n_pixels = MAX(png_resolution)
          ENDIF
        ENDFOR ; loop over cuts

        IF (*i).png_out THEN BEGIN
          IF (*i).mode EQ 2 THEN BEGIN
            new_theta = INDGEN(n_pixels)*2.0*!PI/(n_pixels-1)
            new_phi = INDGEN(n_pixels)*2.0*!PI/(n_pixels-1)
            
            old_theta = FLTARR((*i).z_res+1,(*i).phi_res,/NOZERO)
            old_phi = FLTARR((*i).z_res+1,(*i).phi_res,/NOZERO)

            FOR x = 0, 1 DO BEGIN
              hdat_norm = REFORM(contour_levels(hull_data[*,*,x,sp,n],$
                n_levels,min=(*i).data_minmax[0,sp],max=(*i).data_minmax[1,sp],/img,$
                no_fix_zero=(NOT (*i).fix_color_zero)),[(*i).z_res+1,(*i).phi_res])

              interpolate_to_theta = (*i).geom_type NE 0
              IF interpolate_to_theta THEN BEGIN
                ; interpolate from regular z grid 
                ; (irregular theta grid) to high res., regular theta
                old_theta = FLTARR((*i).z_res+1,(*i).phi_res,/NOZERO)
                FOR p = 0, (*i).phi_res - 1 DO BEGIN
                  old_theta[*,p] = [(*i).space_theta[*,$
                    x*((*i).x_res_red-1)],(*i).space_theta[0,$
                    x*((*i).x_res_red-1)]]
                  old_phi[*,p] = 2.0 * !PI * p * $
                    (*i).phi_range / (*i).phi_res
                ENDFOR

                TRIANGULATE, old_theta, old_phi, triangulation, TOLERANCE=0

                hdat_ip = TRIGRID(old_theta,old_phi,hdat_norm,$
                  triangulation,$;EXTRAPOLATE=boundary,$
                  NX=n_pixels,NY=n_pixels,XOUT=xarr1,YOUT=yarr1)
              ENDIF ELSE BEGIN
                hdat_thetaip = INTERPOLATE(TEMPORARY(hdat_norm),$
                  INDGEN(n_pixels)*((*i).z_res/FLOAT(n_pixels-1)),$
                  INDGEN((*i).phi_res),/GRID)
                hdat_ip = INTERPOLATE(TEMPORARY(hdat_thetaip),$
                  INDGEN(n_pixels)*(((*i).phi_res-1)/FLOAT(n_pixels-1)))
              ENDELSE

              ; using a transpose to have a phi-z plot (i.e., phi as x axis)
              img = BYTARR(n_pixels,n_pixels)
              img[0,0] = TRANSPOSE(hdat_ip)

              frame = REVERSE(img,2)

              IF N_ELEMENTS(image) NE 1L * n_pixels * n_pixels THEN $
                image = BYTARR(3,n_pixels,n_pixels)

              TVLCT, red, green, blue, /GET
              image[0,*,*] = red[frame[*,*]]
              image[1,*,*] = green[frame[*,*]]
              image[2,*,*] = blue[frame[*,*]]
;              image[3,*,*] = 255B

              png_name = png_name_base + 'hull' + rm0es(x) + '_' + $
                framenr_str + '.png'

              WRITE_PNG, png_name, image, /ORDER
           ENDFOR
           IF (framenr EQ 0) THEN write_inc_file, diag, filename_run, framenr_str
          ENDIF

          framenr += 1
        ENDIF ELSE plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
      ENDFOR

      IF (*i).png_out THEN SET_PLOT, 'X'
      set_output, diag, sp_id, /reset, eps=eps
    ENDFOR
  ENDELSE

  IF (*i).png_out THEN BEGIN
    !X.MARGIN = xmarg
    !Y.MARGIN = ymarg
  ENDIF

  threshold = 100
  mem = mem_usage(/highw,threshold=threshold)
  IF (!QUIET EQ 1) AND threshold THEN PRINT, $
    'maximum memory used by ' + (*diag).name + ': ' + mem

END

;######################################################################

PRO write_inc_file, diag, filename_run, framenr_str

  COMMON global_vars

  i = (*diag).internal_vars

  ; Important note:
  ; povray requires a torus with major radius 1 to employ
  ; the map_type 5 (torus mapping), thus we normalize to torR
  ; and correct for/translate by torZ/torR *afterwards*
  lref = (*i).torR

  filename = gui.out.out_path+(*diag).name + filename_run + '_'+$
           framenr_str + '.inc'
  png_name_base = gui.out.out_path + (*diag).name + filename_run + '_'

  print, 'writing '+filename
  OPENW,lun,filename,/GET_LUN

  PRINTF, lun, '// finish of the torus'
  PRINTF, lun, '#declare tor_amb = 0.25;'
  PRINTF, lun, '#declare tor_dif = 0.40;'
  PRINTF, lun, '#declare tor_phg = 0.05;'
  PRINTF, lun, ''

  PRINTF, lun, ''
  PRINTF, lun, '#declare major_R = '+rm0es((*i).torR)+';'
  PRINTF, lun, '#declare major_Z = '+rm0es((*i).torZ)+';'
  PRINTF, lun, '#declare Lref    = '+rm0es(lref)+';'
  PRINTF, lun, ''
  PRINTF, lun, '//rescale torus coordinates by huge factor'  
  PRINTF, lun, '//to avoid floating point errors, then' 
  PRINTF, lun, '//scale to torus with unit major radius'   
  PRINTF, lun, '//to apply the textures and finally normalize'
  PRINTF, lun, '//to desired scale'  
  PRINTF, lun, '#declare rsc = 1000.0;'
  PRINTF, lun, ''

  hull_name = ['inner','outer']
  spline_type = ['linear','linear'] ;['linear','quadratic']
  FOR ihull=0,1 DO BEGIN
     ix = ihull*((*i).x_res_red-1)
     PRINTF, lun, '#declare '+hull_name[ihull]+'_hull = lathe{'
     PRINTF, lun, ' '+spline_type[ihull]+'_spline // Spline_Type'
     PRINTF, lun, '  '+rm0es((*i).z_res+1)+','
     min_r = 1000
     max_r = 0
     min_z = 1000
     max_z = 0  

     FOR z=0,(*i).z_res-1 DO BEGIN
        IF ((*i).geom_type EQ 0) THEN BEGIN
           last_r = ((*i).torR + (*i).space_r[ix]*COS((*i).space_theta[z]))
           last_z = ((*i).torZ + (*i).space_r[ix]*SIN((*i).space_theta[z]))
        ENDIF ELSE BEGIN
           last_r = ((*i).torR + (*i).space_r[z,ix]*COS((*i).space_theta[z,ix]))
           last_z = ((*i).space_r[z,ix]*SIN((*i).space_theta[z,ix]));+(*i).torZ
        ENDELSE
; maybe, it's better to use the R/Z boundaries being computed
; in the previous code parts instead of "measuring" them here
; (cases exist where the cuts do not perfectly fit)
        if (last_r lt min_r) then min_r=last_r
        if (last_r gt max_r) then max_r=last_r
        if (last_z lt min_z) then min_z=last_z
        if (last_z gt max_z) then max_z=last_z     
        PRINTF, lun, '<'+rm0es(last_r)+'*rsc, '+rm0es(last_z)+'*rsc>'
     ENDFOR
     ; repeat inner z point to close boundary curve
     z =0 
     IF ((*i).geom_type EQ 0) THEN BEGIN
        last_r = ((*i).torR+(*i).space_r[ix]*COS((*i).space_theta[z]))
        last_z = ((*i).torZ+(*i).space_r[ix]*SIN((*i).space_theta[z]))
     ENDIF ELSE BEGIN
        last_r = ((*i).torR+(*i).space_r[z,ix]*COS((*i).space_theta[z,ix]))
        last_z = (*i).space_r[z,ix]*SIN((*i).space_theta[z,ix]) ;+(*i).torZ
     ENDELSE
  
     PRINTF, lun, '<'+rm0es(last_r)+'*rsc, '+rm0es(last_z)+'*rsc>'     
     PRINTF, lun, '  // sturm'
     PRINTF, lun, ''
     PRINTF, lun, '//scale to "unit" torus'
     PRINTF, lun, '  scale<1./(rsc*major_R),1./(rsc*major_R),1./(rsc*major_R)>'
     PRINTF, lun, ''
     PRINTF, lun, '  texture{'
     png_name = png_name_base + 'hull'+rm0es(ihull)+'_' + framenr_str + '.png'
     PRINTF, lun, '    pigment{'
     PRINTF, lun, '      image_map{'
     PRINTF, lun, '        png "'+png_name+'"'
     PRINTF, lun, '        map_type 5'
     PRINTF, lun, '        interpolate 2'
     PRINTF, lun, '      }'
     PRINTF, lun, '    }'
     PRINTF, lun, '  finish { ambient tor_amb phong tor_phg diffuse tor_dif}'
  ;  PRINTF, lun, '  rotate<0,90,0>'
     PRINTF, lun, '  } // end of texture'
     PRINTF, lun, ''
     PRINTF, lun, '  scale<major_R/Lref,major_R/Lref,major_R/Lref>'     
     PRINTF, lun, '  rotate <0.0,'+rm0es(-(*i).phi_cut[0]*180.0/!PI)+',0.0>'
     PRINTF, lun, '  translate <0.0,major_Z/Lref,0.0>'
     PRINTF, lun, '} // end of lathe object'
  ENDFOR

  PRINTF, lun, ''
  PRINTF, lun, '#declare offset1 = 0.00001;'
  PRINTF, lun, '#declare offset2 = 0.001;'
  PRINTF, lun, '#declare min_r = '+rm0es(min_r)+'/Lref-offset2;'
  PRINTF, lun, '#declare max_r = '+rm0es(max_r)+'/Lref+offset2;'
  PRINTF, lun, '#declare min_z = '+rm0es(min_z)+'/Lref-offset2;'
  PRINTF, lun, '#declare max_z = '+rm0es(max_z)+'/Lref+offset2;'
  PRINTF, lun, ''

  FOR icut = 0, (*i).n_cuts-1 DO BEGIN
     IF (((*i).phi_cut[icut] GT 0.0) AND $
         ((*i).phi_cut[icut] LE !PI)) THEN $
            sign = '+' ELSE sign = '-'
     
     PRINTF, lun, '// A plane, which gets cut '+rm0es(icut)+$
             ' projected onto.'
     PRINTF, lun, '#declare cut'+rm0es(icut)+' ='
     PRINTF, lun, 'plane { <0, 0, '+sign+'1>, '+sign+'offset1'
     PRINTF, lun, '  texture{'
     PRINTF, lun, '    pigment {'
     PRINTF, lun, '      image_map{'
     png_name = png_name_base + 'cut' + rm0es(icut) + '_' + framenr_str + '.png'
     PRINTF, lun, '        png "'+png_name+'"'
     PRINTF, lun, '        map_type    0  // planar projection'
     PRINTF, lun, '        interpolate 2'
     PRINTF, lun, '        once           // do not checker the texture'
     PRINTF, lun, '      }'
     PRINTF, lun, '    }'
     PRINTF, lun, '    finish {ambient tor_amb diffuse tor_dif phong tor_phg }'
     PRINTF, lun, '  }'
     PRINTF, lun, '  scale <max_r-min_r,max_z-min_z,1.0>'
     PRINTF, lun, '  translate <min_r,min_z+major_Z/Lref,0.0>'
     PRINTF, lun, '  rotate <0.0,'+rm0es(-(*i).phi_cut[icut]*180.0/!PI)+',0.0>'
     PRINTF, lun, '}'
     PRINTF, lun, '' 
  ENDFOR

  PRINTF, lun, ''
  PRINTF, lun, '#declare plasma = union{'
  PRINTF, lun, ' intersection {'
  IF ((*i).n_cuts EQ 2) THEN BEGIN
     PRINTF, lun, '    difference {'
     PRINTF, lun, '        lathe{ outer_hull }'
  ENDIF
  PRINTF, lun, '	lathe{ inner_hull }'
  IF ((*i).n_cuts EQ 2) THEN PRINTF, lun, '    }'

  diff_phi = ((*i).phi_cut[(*i).n_cuts-1]-(*i).phi_cut[0])
  IF diff_phi GT !PI THEN $
     PRINTF, lun, '    union {'
  FOR iicut=0,(diff_phi NE !PI) DO BEGIN
     icut = iicut*((*i).n_cuts-1)
     PRINTF, lun, '      box { <-max_r,min_z,offset1>, <max_r,max_z,max_r>'
     PRINTF, lun, '            translate <0.0,major_Z/Lref,0.0>'
     phi_box = (-(*i).phi_cut[icut]-(icut GT 0)*!PI)*180.0/!PI
     PRINTF, lun, '            rotate <0.0,'+rm0es(phi_box)+',0.0>'
     PRINTF, lun, '      }'
  ENDFOR
  IF diff_phi GT !PI THEN $
     PRINTF, lun, '    }'  

;  PRINTF, lun, '    box{ plasma_bounding_box }'
  PRINTF, lun, ' }'
  FOR icut = 0, (*i).n_cuts-1 DO $
     PRINTF, lun, ' plane{ cut'+rm0es(icut)+' }'
  PRINTF, lun, '}'

  FREE_LUN, lun

END
