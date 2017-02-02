FUNCTION magfline_info

  RETURN, {$
    type      : 'mom',$
    title     : 'Magnetic field lines',$
    help_text : ['Computes Poincare sections of perturbed magnetic '+$
                 'field lines.'],$
    ext_vars  : [['n_lin','0','number of lines in x[/y] direction; '+$
                  'default: [1,1]; optionally, add initial '+$
                  'positions: [Nx[,Ny],x1,y1,...,xn,yn]'],$
                 ['n_pol','0','number of poloidal cycles per field '+$
                  'line; default: 100'],$
                 ['kyarr','0','considered y modes; default: -1 (all)'],$
                 ['plot_mode','0','bitwise: output per time step; '+$
                  '1: Poincare section; 2: field line diffusivity; '+$
                  '4: B_xy at nz0/2; 8: zoom-to-fit radial box; '+$
                  'default: 1'],$
                 ['res','0','tracing resolution: Nz[,Nx[,Ny]]; use '+$
                  'small prime factors for Nx, Ny for higher speed; '+$
                  'default: [nz0>128,nx0>128,ny0>128]']]}

END

;######################################################################

PRO magfline_init, diag

  COMMON global_vars

  IF par.n_fields LE 1 THEN BEGIN
    PRINT, (*diag).name + ' error: need electromagnetic runs'
    (*diag).selected = 0
    RETURN
  ENDIF
  IF NOT par.nonlinear THEN PRINT, (*diag).name + $
    ' warning: linear runs only with one finite ky or no adapt_lx!'

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.n_lin) LT 1 THEN *e.n_lin = [1,1]
  IF N_ELEMENTS(*e.n_lin) EQ 1 THEN *e.n_lin = [*e.n_lin,1]
  ne_nlin = N_ELEMENTS(*e.n_lin)
  IF (ne_nlin GT 2) AND ((ne_nlin MOD 2) NE 0) THEN $
    *e.n_lin = [(*e.n_lin)[0],1,(*e.n_lin)[1:ne_nlin-1]]
  FOR j = 0, 1 DO IF (*e.n_lin)[j] LT 1 THEN (*e.n_lin)[j] = 1
  IF (N_ELEMENTS(*e.n_lin) - 2) / 2 NE (*e.n_lin)[0] * (*e.n_lin)[1] $
    THEN *e.n_lin = (*e.n_lin)[0:1]
  n_lines = (*e.n_lin)[0] * (*e.n_lin)[1]

  IF N_ELEMENTS(*e.n_pol) NE 1 THEN *e.n_pol = 100L
  IF N_ELEMENTS(*e.kyarr) LT 1 THEN *e.kyarr = -1
  IF (*e.kyarr)[0] EQ -1 THEN *e.kyarr = INDGEN(par.nky0)
  *e.kyarr = (*e.kyarr)[WHERE(*e.kyarr LT par.nky0)]
  IF N_ELEMENTS(*e.plot_mode) NE 1 THEN *e.plot_mode = 1

  xyz_res = [par.nx0,par.ny0,par.nz0] > 128
  CASE N_ELEMENTS(*e.res) OF
    3    : xyz_res = [(*e.res)[2],(*e.res)[0:1]]
    2    : xyz_res[[0,2]] = [(*e.res)[1],(*e.res)[0]]
    1    : xyz_res[2] = (*e.res)[0]
    ELSE :
  ENDCASE

  plot_mode = BYTARR(4)
  IF *e.plot_mode MOD 2 GE 1 THEN BEGIN
    plot_mode[0] = 1 ; plot section
    *e.plot_mode -= 1
  ENDIF
  IF *e.plot_mode MOD 4 GE 1 THEN BEGIN
    plot_mode[1] = 1 ; plot diffusivity
    *e.plot_mode -= 2
  ENDIF
  IF *e.plot_mode MOD 8 GE 1 THEN BEGIN
    plot_mode[2] = 1 ; plot field
    *e.plot_mode -= 4
  ENDIF
  IF *e.plot_mode EQ 8 THEN plot_mode[3] = 1 ; zoom box
  IF TOTAL(plot_mode[0:2]) LE 0 THEN BEGIN
    PRINT, (*diag).name + ' error: no valid plot mode selected'
    (*diag).selected = 0
    RETURN
  ENDIF
  IF TOTAL(plot_mode[0:1]) EQ 0 THEN *e.n_pol = 0

  IF *e.n_pol LE 0 THEN xyz_res = [par.nx0,par.ny0,par.nz0]
  x_res = xyz_res[0] > (par.nx0 / 2 * 2)
  y_res = xyz_res[1] > (par.ny0 / 2 * 2)
  z_res = xyz_res[2] > (par.nz0 / 2 * 2)

  ; set initial (seed) field line positions
  center_seed = 1 ; radial seed line positioning
  fline_init = FLTARR(n_lines,2,/NOZERO)
  IF N_ELEMENTS(*e.n_lin) EQ 2 THEN BEGIN
    FOR x = 0, (*e.n_lin)[0] - 1 DO BEGIN
      FOR y = 0, (*e.n_lin)[1] - 1 DO $
        fline_init[y+x*(*e.n_lin)[1],*] = $
        [x_res*(x+center_seed*0.5)/FLOAT((*e.n_lin)[0]),$
        y_res*(y+0.5)/FLOAT((*e.n_lin)[1])]
    ENDFOR
  ENDIF ELSE FOR l = 0, n_lines - 1 DO $
    fline_init[l,*] = [(*e.n_lin)[2+2*l],(*e.n_lin)[3+2*l]]

  verbose = (1L * x_res * y_res * z_res GT 50000000L) OR $
    (1L * n_lines * *e.n_pol * z_res GE 5000000L)

  i = set_internal_vars(diag,{$
  verbose     : verbose,$
  n_lines     : n_lines,$
  n_lines_x   : (*e.n_lin)[0],$
  n_lines_y   : (*e.n_lin)[1],$
  n_pol       : *e.n_pol,$
  kyarr       : *e.kyarr,$
  x_res       : x_res,$
  y_res       : y_res,$
  z_res       : z_res,$
  fline_init  : fline_init,$
  poincare_id : PTR_NEW(),$
  time_id     : PTR_NEW(),$
  rad_bc_id   : PTR_NEW(),$
  plot_mode   : plot_mode,$
  B_xy_id     : PTR_NEW()})

  fft_format, kxky=1

END

;######################################################################

PRO magfline_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars  

  nky = N_ELEMENTS((*i).kyarr)

  ; auxiliary quantities, defined for speed gain
  n_lines = (*i).n_lines
  x_res = (*i).x_res
  y_res = (*i).y_res
  z_res = (*i).z_res
  xres_o_2 = x_res / 2
  yres_o_2 = y_res / 2
  zres_inv = 1.0 / z_res
  dx_inv = x_res / par.lx
  dy_inv = y_res / par.ly
  y_shift_base = 2.0 * !PI * par.shat * dy_inv / dx_inv
  B0_inv_extd = 1.0 / [(*series.geom).Bfield,(*series.geom).Bfield[0]]

  ; --- calculate perpendicular magnetic field fluctuations ---

  ; extend magnetic field fluctuations parallely (quasi-periodic BC)
  bigN = par.shat NE 0.0 ? ROUND(2*!PI*par.shat*par.lx/par.ly) : 0
  IF bigN LT 0 THEN PRINT, (*diag).name + $
    ' warning: negative mode mapping (BC) not tested'
  kx_mod = LONARR(par.nkx0,par.nky0) - 1
  FOR y = 0, par.nky0 - 1 DO kx_mod[0:par.nkx0-par.nkx0/2-1,y] = $
    INDGEN(par.nkx0-par.nkx0/2) + bigN * y
  IF (WHERE(kx_mod GE par.nkx0-par.nkx0/2))[0] NE -1 THEN $
    kx_mod[WHERE(kx_mod GE par.nkx0-par.nkx0/2)] = - 2 * par.nkx0
  FOR y = 0, par.nky0 - 1 DO kx_mod[par.nkx0/2+1:*,y] = $
    - (par.nkx0 - par.nkx0 / 2 - 1 - $
    INDGEN(par.nkx0-par.nkx0/2-1)) + bigN * y + par.nkx0
  IF (WHERE(kx_mod GE par.nkx0))[0] NE -1 THEN $
    kx_mod[WHERE(kx_mod GE par.nkx0)] = - 2 * par.nkx0

  data_extd = COMPLEXARR(par.nkx0,par.nky0,par.nz0+1)
  data_extd[0,0,0] = *mom[0,1].kxky
  FOR x = 0, par.nkx0 - 1 DO BEGIN
    kx_mod_applies = WHERE(kx_mod[x,*] GE 0)

    IF kx_mod_applies[0] GE 0 THEN $
      data_extd[x,kx_mod_applies,par.nz0] = $
      (*mom[0,1].kxky)[kx_mod[x,kx_mod_applies],kx_mod_applies,$
      REBIN([0],N_ELEMENTS(kx_mod_applies))]
  ENDFOR

  ; divide A_par by B_0
  data_extd *= REBIN(REFORM(B0_inv_extd,[1,1,par.nz0+1]),$
    [par.nkx0,par.nky0,par.nz0+1])

  IF bigN MOD 2 THEN BEGIN
    phase_fac = 1 - 2 * ((INDGEN(par.nky0)*bigN) MOD 2)
    data_extd[*,*,par.nz0] *= $
      REBIN(REFORM(phase_fac,[1,par.nky0]),[par.nkx0,par.nky0])
  ENDIF

  IF (*i).verbose THEN PRINT, '    ' + (*diag).name + $
    ': computing field fluctuations'

  ; interpolate perpendicular magnetic field fluctuations
  apar_ip = z_res NE par.nz0 ? INTERPOLATE(TEMPORARY(data_extd),$
    INDGEN(z_res+1)*(zres_inv*par.nz0)) : TEMPORARY(data_extd)

  B_xy_kxky = COMPLEXARR(x_res,nky,z_res+1)
  temp_data = COMPLEXARR(x_res,y_res,z_res+1)

  ; obtain real space B_x (with Fourier interpolation)
  B_xy_kxky[0,0,0] = $
    COMPLEX(0,1) * REBIN(REFORM((*series.ky)[(*i).kyarr],[1,nky]),$
    [par.nkx0/2,nky,z_res+1]) * apar_ip[0:par.nkx0/2-1,(*i).kyarr,*]
  B_xy_kxky[x_res-par.nkx0/2+1,0,0] = $
    COMPLEX(0,1) * REBIN(REFORM((*series.ky)[(*i).kyarr],[1,nky]),$
    [par.nkx0/2-1,nky,z_res+1]) * $
    apar_ip[par.nkx0/2+1:par.nkx0-1,(*i).kyarr,*]
  B_xy_sxky = FFT(B_xy_kxky,DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, nky - 1 DO temp_data[0,(*i).kyarr[y],0] = B_xy_sxky[*,y,*]
  FOR y = 0, nky - 1 DO IF (*i).kyarr[y] NE 0 THEN $
    temp_data[0,y_res-(*i).kyarr[y],0] = CONJ(B_xy_sxky[*,y,*])
  B_x = FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0)
  B_x = FLOAT(TEMPORARY(B_x))

  ; obtain real space B_y (with Fourier interpolation)
  B_xy_kxky[0,0,0] = $
    - COMPLEX(0,1) * REBIN((*series.kx)[0:par.nkx0/2-1,0],$
    [par.nkx0/2,nky,z_res+1]) * apar_ip[0:par.nkx0/2-1,(*i).kyarr,*]
  B_xy_kxky[x_res-par.nkx0+par.nkx0/2+1,0,0] = $
    - COMPLEX(0,1) * REBIN((*series.kx)[par.nkx0/2+1:par.nkx0-1,0],$
    [par.nkx0-par.nkx0/2-1,nky,z_res+1]) * $
    apar_ip[par.nkx0/2+1:par.nkx0-1,(*i).kyarr,*]
  B_xy_sxky = FFT(TEMPORARY(B_xy_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, nky - 1 DO temp_data[0,(*i).kyarr[y],0] = B_xy_sxky[*,y,*]
  FOR y = 0, nky - 1 DO IF (*i).kyarr[y] NE 0 THEN $
    temp_data[0,y_res-(*i).kyarr[y],0] = CONJ(B_xy_sxky[*,y,*])
  B_xy_sxky = 0
  B_y = FFT(TEMPORARY(temp_data),DIMENSION=2,/INVERSE,DOUBLE=0)
  B_y = FLOAT(TEMPORARY(B_y))

  IF (*i).verbose THEN PRINT, '    ' + (*diag).name + $
    ' memory info: ' + mem_usage()

  ; rescale magnetic field for linear runs
  IF NOT par.nonlinear THEN BEGIN
    B_xy_sqd = TOTAL(B_x^2) + TOTAL(B_y^2)
    PRINT, '    ' + (*diag).name + ' info: rescaling B_xy^2 from ' + $
      rm0es(B_xy_sqd) + ' to ', FORMAT='(A,$)'
    scale_fac = 10.0 / SQRT(B_xy_sqd)
    B_x *= scale_fac
    B_y *= scale_fac

    PRINT, rm0es(TOTAL(B_x^2)+TOTAL(B_y^2))
  ENDIF

  (*i).time_id = store_step((*i).time_id,mom_time)

  ; store field at outboard midplane for contour plot
  IF (*i).plot_mode[2] THEN BEGIN
    store_B_xy = FLTARR(x_res,y_res,2,/NOZERO)
    store_B_xy[0,0,0] = B_x[*,*,z_res/2]
    store_B_xy[0,0,1] = B_y[*,*,z_res/2]

    (*i).B_xy_id = store_step((*i).B_xy_id,store_B_xy)

    IF TOTAL((*i).plot_mode[0:1]) EQ 0 THEN RETURN
  ENDIF

  ; include factors for calculating field line shift
  shift_fac = 2 * !PI * par.q0 * zres_inv
  B_x *= dx_inv * shift_fac
  B_y *= dy_inv * shift_fac

  ; --- follow field lines ---

  IF (*i).verbose THEN PRINT, '    ' + (*diag).name + $
    ': starting field line integration'

  poincare = FLTARR(n_lines,(*i).n_pol,2,/NOZERO)
  fline_x = (*i).fline_init[*,0]
  fline_y = (*i).fline_init[*,1]
  rad_bc = INTARR(n_lines,(*i).n_pol)

  FOR p = 0, (*i).n_pol - 1 DO BEGIN
    FOR fline_z = 0, z_res - 1 DO BEGIN
      fline_z_0 = REBIN([fline_z],n_lines)

      ; nearest neighbors on grid, relative field line position
      pos_x1 = FIX(fline_x)
      pos_x2 = (pos_x1 + 1) MOD x_res
      pos_y1 = FIX(fline_y)
      pos_y2 = (pos_y1 + 1) MOD y_res
      dx1 = fline_x - pos_x1
      dx2 = 1.0 - dx1
      dy1 = fline_y - pos_y1
      dy2 = 1.0 - dy1

      ; field fluctuation at origin
      B_interp_x = dy2 * (dx2 * B_x[pos_x1,pos_y1,fline_z_0] + $
        dx1 * B_x[pos_x2,pos_y1,fline_z_0]) + $
        dy1 * (dx2 * B_x[pos_x1,pos_y2,fline_z_0] + $
        dx1 * B_x[pos_x2,pos_y2,fline_z_0])
      B_interp_y = dy2 * (dx2 * B_y[pos_x1,pos_y1,fline_z_0] + $
        dx1 * B_y[pos_x2,pos_y1,fline_z_0]) + $
        dy1 * (dx2 * B_y[pos_x1,pos_y2,fline_z_0] + $
        dx1 * B_y[pos_x2,pos_y2,fline_z_0])

      ; field line position at half way
      fline_x_hw = fline_x + 0.5 * B_interp_x
      fline_y_hw = fline_y + 0.5 * B_interp_y

      ; apply x/y boundary conditions
      fline_x_hw = (fline_x_hw + x_res) MOD x_res
      fline_y_hw = ((fline_y_hw MOD y_res) + y_res) MOD y_res

      ; nearest half way neighbors, relative field line position
      pos_x1 = FIX(fline_x_hw)
      pos_x2 = (pos_x1 + 1) MOD x_res
      pos_y1 = FIX(fline_y_hw)
      pos_y2 = (pos_y1 + 1) MOD y_res
      dx1 = fline_x_hw - pos_x1
      dx2 = 1.0 - dx1
      dy1 = fline_y_hw - pos_y1
      dy2 = 1.0 - dy1
      fline_z_fw = fline_z_0 + 1

      ; field fluctuation at half way
      B_interp_x = $
        0.5 * (dy2 * (dx2 * B_x[pos_x1,pos_y1,fline_z_0] + $
        dx1 * B_x[pos_x2,pos_y1,fline_z_0]) + $
        dy1 * (dx2 * B_x[pos_x1,pos_y2,fline_z_0] + $
        dx1 * B_x[pos_x2,pos_y2,fline_z_0])) + $
        0.5 * (dy2 * (dx2 * B_x[pos_x1,pos_y1,fline_z_fw] + $
        dx1 * B_x[pos_x2,pos_y1,fline_z_fw]) + $
        dy1 * (dx2 * B_x[pos_x1,pos_y2,fline_z_fw] + $
        dx1 * B_x[pos_x2,pos_y2,fline_z_fw]))
      B_interp_y = $
        0.5 * (dy2 * (dx2 * B_y[pos_x1,pos_y1,fline_z_0] + $
        dx1 * B_y[pos_x2,pos_y1,fline_z_0]) + $
        dy1 * (dx2 * B_y[pos_x1,pos_y2,fline_z_0] + $
        dx1 * B_y[pos_x2,pos_y2,fline_z_0])) + $
        0.5 * (dy2 * (dx2 * B_y[pos_x1,pos_y1,fline_z_fw] + $
        dx1 * B_y[pos_x2,pos_y1,fline_z_fw]) + $
        dy1 * (dx2 * B_y[pos_x1,pos_y2,fline_z_fw] + $
        dx1 * B_y[pos_x2,pos_y2,fline_z_fw]))

      ; field line position at full way
      fline_x += B_interp_x
      fline_y += B_interp_y

      ; apply x boundary condition
      fline_x_orig = fline_x
      fline_x = (fline_x + x_res) MOD x_res
      ; store radial box exceedances
      dfline_x_orig = fline_x_orig - fline_x
      out_of_box_R = WHERE(dfline_x_orig GT 1)
      IF out_of_box_R[0] NE -1 THEN rad_bc[out_of_box_R,p:*] += 1
      out_of_box_L = WHERE(dfline_x_orig LT -1)
      IF out_of_box_L[0] NE -1 THEN rad_bc[out_of_box_L,p:*] -= 1

      ; apply y/z boundary conditions
      IF fline_z EQ z_res - 1 THEN fline_y -= $
        y_shift_base * (fline_x - xres_o_2)
      fline_y = ((fline_y MOD y_res) + y_res) MOD y_res
    ENDFOR

    ; store Poincare data
    poincare[0,p,0] = fline_x
    poincare[0,p,1] = fline_y

    IF (*i).verbose THEN print_progress_cl, p + 1, (*i).n_pol, $
      first=(p EQ 0), /erase_at_100, header_str='      '
  ENDFOR

  (*i).poincare_id = store_step((*i).poincare_id,poincare)
  (*i).rad_bc_id = store_step((*i).rad_bc_id,rad_bc)

END

;######################################################################

PRO magfline_output, diag

  COMMON global_vars

  ; -------------------------------------------------------
  ; -------------------- plot settings ---------------------
  ; --------------------------------------------------------
  color_by_line = 1     ; else: color by poloidal turn
  show_init = 0         ; mark initial positions
  add_hist = 0          ; show Poincare histogram
                        ; (only with Poincare section)
  add_rad_disp = 0      ; show radial displacement
                        ; (only with field line diffusivity)
  diff_plot_simple = -1 ; plot only the line average to
                        ; avoid clutter and slow display
                        ; default: -1 (simple >1000 lines)
  add_time_vs_diff = 1  ; show time resolved diffusivity
                        ; (only with field line diffusivity,
                        ; multiple time steps)
  add_vec_plot = 1      ; show magnetic field vector plot
                        ; (only with B_xy plotting on)
  arrow_seed = [20,20]  ; number of vector plot arrows (x/y)
  ; --------------------------------------------------------

  i = (*diag).internal_vars

  time = store_step((*i).time_id,/get)
  IF TOTAL((*i).plot_mode[0:1]) GT 0 THEN BEGIN
    poincare = store_step((*i).poincare_id,/get)
    rad_bc = store_step((*i).rad_bc_id,/get)
  ENDIF
  IF (*i).plot_mode[2] THEN B_xy = store_step((*i).B_xy_id,/get)

  plot_tvdiff = 0
  IF add_time_vs_diff THEN BEGIN
    fl_diff_t = FLTARR(series.step_count,2,/NOZERO)
    tarr = FLTARR(series.step_count,/NOZERO)
  ENDIF

  ysize = 17 * par.ly / par.lx
  IF (*i).plot_mode[2] THEN ysize = ysize * (add_vec_plot ? 1 : 2)
  ysize = ysize < 25 > 12

  set_output, diag, /ps, coltable=33, ysize=ysize
  ymarg = !Y.MARGIN
  !Y.MARGIN = [ymarg[0],5]

  FOR n = 0, series.step_count - 1 DO BEGIN
    ; --- plot Poincare sections ---
    IF (*i).plot_mode[0] THEN BEGIN
      ybox = [-0.5*par.ly,0.5*par.ly]
      IF (*i).plot_mode[3] THEN xbox = $
        ([MIN((*poincare[n])[*,*,0]),MAX((*poincare[n])[*,*,0])] / $
        (*i).x_res - 0.5) * par.lx ELSE xbox = [-0.5*par.lx,0.5*par.lx]

      PLOT, xbox, ybox, /NODATA, COLOR=1, /XSTYLE, /YSTYLE, $
        TICKLEN=-0.5*!P.TICKLEN, ISOTROPIC=(1-(*i).plot_mode[3]), $
        XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
        TITLE='!6Poincare section, t='+rm0es(*time[n])

      IF color_by_line THEN BEGIN
        carr = (contour_levels([0,((*i).n_lines>2)-1],(*i).n_lines>2))[*,1]
        FOR l = 0L, (*i).n_lines - 1L DO BEGIN
          poincare_x = $
            ((*poincare[n])[l,*,0] / (*i).x_res - 0.5) * par.lx
          poincare_y = $
            ((*poincare[n])[l,*,1] / (*i).y_res - 0.5) * par.ly
          OPLOT, poincare_x, poincare_y, COLOR=carr[l], PSYM=3
        ENDFOR
      ENDIF ELSE BEGIN
        carr = (contour_levels([0,((*i).n_pol>2)-1],(*i).n_pol>2))[*,1]
        FOR p = 0L, (*i).n_pol - 1L DO BEGIN
          poincare_x = $
            ((*poincare[n])[*,p,0] / (*i).x_res - 0.5) * par.lx
          poincare_y = $
            ((*poincare[n])[*,p,1] / (*i).y_res - 0.5) * par.ly
          OPLOT, poincare_x, poincare_y, COLOR=carr[p], PSYM=3
        ENDFOR
      ENDELSE

      ; mark initial (seed) positions
      IF show_init THEN BEGIN
        fline_init_xy = (*i).fline_init * $
          [REBIN([par.lx/(*i).x_res],(*i).n_lines),$
          REBIN([par.ly/(*i).y_res],(*i).n_lines)] - 0.5 * $
          [REBIN([par.lx],(*i).n_lines),REBIN([par.ly],(*i).n_lines)]
        OPLOT, fline_init_xy[*,0], fline_init_xy[*,1], COLOR=1, PSYM=6
      ENDIF

      ; plot Poincare histograms
      IF add_hist THEN BEGIN
        ; create histograms on a plot-ready grid
        cut_res = [200,200]
        hist_x = (FINDGEN(cut_res[0]) / cut_res[0] - 0.5) * par.lx
        hist_y = (FINDGEN(cut_res[1]) / cut_res[1] - 0.5) * par.ly
        xy_grid = LONARR(cut_res[0],cut_res[1])

        poincare_x = ROUND((*poincare[n])[*,*,0]*$
          (cut_res[0]/FLOAT((*i).x_res))) MOD cut_res[0]
        poincare_y = ROUND((*poincare[n])[*,*,1]*$
          (cut_res[1]/FLOAT((*i).y_res))) MOD cut_res[1]

        FOR p = 0L, (*i).n_pol - 1L DO FOR l = 0L, (*i).n_lines - 1L DO $
          xy_grid[poincare_x[l,p],poincare_y[l,p]]++

        lev_col = contour_levels(xy_grid,40)
        CONTOUR, xy_grid, hist_x, hist_y, $
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /XSTYLE, /YSTYLE, $
          XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
          TITLE='Poincare histogram, t='+rm0es(*time[n]), $
          TICKLEN=-0.5*!P.TICKLEN, /FILL, /ISOTROPIC

        plot_colorbar, lev_col, orientation=1, $
          position=[0.19,0.95,0.94,0.97]
      ENDIF
    ENDIF

    ; --- plot field line diffusivity ---
    IF (*i).plot_mode[1] THEN BEGIN
      IF diff_plot_simple LT 0 THEN $
        diff_plot_simple = (*i).n_lines_x GT 1000

      ; reverse radial boundary condition
      ; (i.e., remove L_x jumps from the radial displacement)
      poincare_mod = ((*poincare[n])[*,*,0] / (*i).x_res - 0.5) * par.lx
      poincare_mod += (*rad_bc[n]) * par.lx
      ; radial displacement
      dr = poincare_mod[*,*] - REBIN(((*i).fline_init[*,0]/$
        (*i).x_res-0.5)*par.lx,[(*i).n_lines,(*i).n_pol])

      ; compute field line diffusivity
      dr_sqd = dr^2 / (2 * !PI * par.q0 * (*i).n_lines_y)
      fl_diff_ysum = FLTARR((*i).n_lines_x,(*i).n_pol,/NOZERO)
      FOR l = 0, (*i).n_lines_x - 1 DO fl_diff_ysum[l,0] = REFORM($
        TOTAL(dr_sqd[l*(*i).n_lines_y+LINDGEN((*i).n_lines_y),*],1),$
        [1,(*i).n_pol])
      ; divide by poloidal turns
      fl_diff = TEMPORARY(fl_diff_ysum) * REBIN(REFORM($
        1.0/(LINDGEN((*i).n_pol)+1),[1,(*i).n_pol]),[(*i).n_lines_x,(*i).n_pol])

      fl_diff_xavg = TOTAL(fl_diff,1) / (*i).n_lines_x
      diff_min = diff_plot_simple ? MIN(fl_diff_xavg,MAX=diff_max) : $
        MIN(fl_diff,MAX=diff_max)
      PLOT, [1,(*i).n_pol], [diff_min,diff_max], COLOR=1, $
        /NODATA, /XSTYLE, YSTYLE=(1-diff_plot_simple), $
        XTITLE='!6poloidal turns', YTITLE='!6D!Dfl!N'

      IF NOT diff_plot_simple THEN BEGIN
        carr = (*i).n_lines_x EQ 1 ? 1 : $
          (contour_levels([0,(*i).n_lines_x-1],(*i).n_lines_x))[*,1]
        FOR l = 0, (*i).n_lines_x - 1 DO $
          OPLOT, LINDGEN((*i).n_pol) + 1, fl_diff[l,*], COLOR=carr[l]
      ENDIF

      IF (*i).n_lines_x NE 1 THEN BEGIN
        npol_axis = LINDGEN((*i).n_pol) + 1

        LOADCT, 41, FILE='internal/colortable.tbl'

        IF diff_plot_simple THEN BEGIN
          OPLOT, npol_axis, fl_diff_xavg, COLOR=1
        ENDIF ELSE BEGIN
          lthick = 1 + 7 * ((*i).n_lines_x GT 5)
          OPLOT, npol_axis, fl_diff_xavg, COLOR=0, THICK=lthick
          OPLOT, npol_axis, fl_diff_xavg, COLOR=1, LINE=1, THICK=lthick
        ENDELSE

        LOADCT, 33, FILE='internal/colortable.tbl'
      ENDIF

      fl_diff_avg = TOTAL(fl_diff_xavg) / (*i).n_pol
      PRINT, '<D_fl>: ' + rm0es(fl_diff_avg) + ', D_fl(n_pol=' + $
        rm0es((*i).n_pol) + '): ' + rm0es(fl_diff_xavg[(*i).n_pol-1])

      IF add_time_vs_diff AND (series.step_count GT 1) THEN BEGIN
        plot_tvdiff = 1
        fl_diff_t[n,0] = fl_diff_avg
        fl_diff_t[n,1] = fl_diff_xavg[(*i).n_pol-1]
        tarr[n] = *time[n]
      ENDIF

      IF add_rad_disp THEN BEGIN
        dr_min = MIN(dr,MAX=dr_max)
        PLOT, [1,(*i).n_pol], [dr_min,dr_max], COLOR=1, $
          /NODATA, /XSTYLE, /YSTYLE, XTITLE='!6poloidal turns', $
          YTITLE='!7D!6r'

        carr = (*i).n_lines_x EQ 1 ? 1 : $
          (contour_levels([0,(*i).n_lines_x-1],(*i).n_lines_x))[*,1]
        FOR l = 0, (*i).n_lines_x - 1 DO OPLOT, $
          LINDGEN((*i).n_pol) + 1, dr[l,*], COLOR=carr[l]
      ENDIF
    ENDIF

    ; --- plot field fluctuations ---
    IF (*i).plot_mode[2] THEN BEGIN
      !P.MULTI = [0,1+add_vec_plot,2]
      IF add_vec_plot THEN BEGIN
        xmarg = !X.MARGIN
        !X.MARGIN -= [3.8,1.2]
      ENDIF

      time_str = ', t=' + rm0es(*time[n])

      FOR j = 0, 1 DO BEGIN
        lev_col = contour_levels((*B_xy[n])[*,*,j],40)
        CONTOUR, (*B_xy[n])[*,*,j], (FINDGEN((*i).x_res) / $
          (*i).x_res - 0.5) * par.lx, (FINDGEN((*i).y_res) / $
          (*i).y_res - 0.5) * par.ly, /FILL, /XSTYLE, /YSTYLE, $
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], $
          XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
          TITLE='!6B!D'+(j EQ 0?'x':'y')+'!N(x,y)/B!D0!N'+$
          (add_vec_plot ? '' : time_str)

        plot_colorbar, lev_col, orientation=1, $
          position=(add_vec_plot ? $
          [0.118+0.5*j,0.95,0.463+0.5*j,0.965] : $
          [0.19,0.95-0.5*j,0.94,0.97-0.5*j])

        PRINT, 'B_' + (j EQ 0 ? 'x' : 'y') + $
          ',max: ' + rm0es(MAX((*B_xy[n])[*,*,j]))
      ENDFOR

      ; add vector plot of the perpendicular field fluctuations
      IF add_vec_plot THEN BEGIN
        !P.MULTI=[1,1,2]
        arrows = vector_plot_arrows((*B_xy[n])[*,*,0],$
          (*B_xy[n])[*,*,1],arrow_seed=arrow_seed)

        ; plot arrows
        LOADCT, 41, FILE='internal/colortable.tbl'
        PLOT, par.lx * 0.5 * [-1,1,1,-1,-1], $
          par.ly * 0.5 * [-1,-1,1,1,-1], /XSTYLE, /YSTYLE, $
          XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
          TITLE='!6B(x,y)'+time_str, COLOR=1
        FOR v = 0, arrow_seed[0] * arrow_seed[1] - 1 DO $
          PLOTS, arrows[v,*,0], arrows[v,*,1], COLOR=1
        LOADCT, 33, FILE='internal/colortable.tbl'

        !X.MARGIN = xmarg
      ENDIF

      !P.MULTI = [0,1,1]
    ENDIF
  ENDFOR

  ; time vs. field line diffusivity
  IF plot_tvdiff THEN BEGIN
    dmin = MIN(fl_diff_t,MAX=dmax)

    LOADCT, 41, FILE='internal/colortable.tbl'

    PLOT, tarr, fl_diff_t[*,0], COLOR=1, /XSTYLE, /YSTYLE, $
      YRANGE=[dmin,dmax], XTITLE=get_var_string(0,/time,/units), $
      YTITLE='!6<D!Dfl!N>,                   '
    AXIS, YAXIS=0, YRANGE=[dmin,dmax], COLOR=2, /YSTYLE, $
      YTITLE='!6             D!Dfl!N(n!Dpol!N='+rm0es((*i).n_pol)+')'
    PLOT, tarr, fl_diff_t[*,0], COLOR=1, /XSTYLE, /YSTYLE, $
      YRANGE=[dmin,dmax], XTITLE=get_var_string(0,/time,/units), $
      YTITLE='!6<D!Dfl!N>,                   ', /NOERASE
    OPLOT, tarr, fl_diff_t[*,1], COLOR=2

    ; compute and plot averages
    d_npolavg_id = PTR_NEW()
    d_lastpol_id = PTR_NEW()
    FOR n = 0, series.step_count - 1 DO BEGIN
      d_npolavg_id = time_avg(d_npolavg_id,fl_diff_t[n,0],tarr[n])
      d_lastpol_id = time_avg(d_lastpol_id,fl_diff_t[n,1],tarr[n])
    ENDFOR
    d_npolavg = time_avg(d_npolavg_id,/avg)
    d_lastpol = time_avg(d_lastpol_id,/avg)
    OPLOT, [tarr[0],tarr[series.step_count-1]], $
      [1,1] * d_npolavg, COLOR=1, LINE=1
    OPLOT, [tarr[0],tarr[series.step_count-1]], $
      [1,1] * d_lastpol, COLOR=1, LINE=2
    PRINT, '<D_fl(pol. turn avg.)>_t = ' + rm0es(d_npolavg)
    PRINT, '<D_fl(last pol. turn)>_t = ' + rm0es(d_lastpol)

    LOADCT, 33, FILE='internal/colortable.tbl'
  ENDIF

  PTR_FREE, time
  IF TOTAL((*i).plot_mode[0:1]) GT 0 THEN PTR_FREE, poincare, rad_bc
  IF (*i).plot_mode[2] THEN BEGIN
    set_output, diag, dat=*B_xy[0]
    PTR_FREE, B_xy
  ENDIF

  !Y.MARGIN = ymarg
  set_output, diag, /reset

END
