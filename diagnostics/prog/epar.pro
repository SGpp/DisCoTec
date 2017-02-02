FUNCTION epar_info

  RETURN, {$
    type      : 'mom',$
    title     : 'Parallel E field',$
    help_text : ['Computes E_par as created by E_perp fluctuations '+$
                 'along perturbed (B_perp) field lines, along with '+$
                 'the d/dt A_par term.'],$
    ext_vars  : [['n_levels','0','set to number of contour levels '+$
                  'for 2D contour plots; 0 for Epar_min,max plot; '+$
                  'default: -1 (Epar_rms)'],$
                 ['sparse_plot','0','plot only every n-th step'],$
                 ['damping','0','[spectral damping exponent,'+$
                  'starting k_xy]; default: [-1,-1(off)]'],$
                 ['stat_dyn','1','plot the static (phi) and dynamic '+$
                  '(A_par) contributions separately'],$
                 ['jpar_mult','1','add plot with E_par * j_par']]}

END

;######################################################################

PRO epar_init, diag

  COMMON global_vars

  IF par.n_fields LE 1 THEN BEGIN
    PRINT, (*diag).name + ' error: need electromagnetic runs'
    (*diag).selected = 0
    RETURN
  ENDIF
  IF NOT par.nonlinear THEN PRINT, (*diag).name + $
    ' warning: linear runs only with one finite ky or no adapt_lx!'

  e = (*diag).external_vars
  IF N_ELEMENTS(*e.n_levels) NE 1 THEN *e.n_levels = -1
  minmax_plot = *e.n_levels EQ 0
  rms_plot = *e.n_levels EQ -1
  IF N_ELEMENTS(*e.sparse_plot) NE 1 THEN *e.sparse_plot = 1
  IF NOT KEYWORD_SET(*e.stat_dyn) THEN *e.stat_dyn = 0
  IF NOT KEYWORD_SET(*e.jpar_mult) THEN *e.jpar_mult = 0

  IF N_ELEMENTS(*e.damping) NE 2 THEN *e.damping = [-1,-1]
  damping = (*e.damping)[1] GE 0
  damp_kx = 0
  damp_ky = 0
  damp_kx_fac = 0
  damp_ky_fac = 0
  IF damping THEN BEGIN
    damp_exp = (*e.damping)[0]
    damp_kx_start = (WHERE((*series.kx)[*,0] GE (*e.damping)[1]))[0]
    IF damp_kx_start GE 0 THEN BEGIN
      damp_kx = 1
      damp_kx_fac = FLTARR(par.nkx0) + 1.0
      damp_kx_fac[damp_kx_start] = $
        ((*series.kx)[damp_kx_start:par.nkx0/2-1,0]/$
        (*series.kx)[damp_kx_start,0])^damp_exp
      damp_kx_fac[par.nkx0/2] = ABS($
        ((*series.kx)[par.nkx0/2:(par.nkx0-(damp_kx_start>1)),0]/$
        (*series.kx)[damp_kx_start,0])^damp_exp)
      damp_kx_fac = REBIN(damp_kx_fac,[par.nkx0,par.nky0,par.nz0])
    ENDIF
    damp_ky_start = (WHERE((*series.ky)[*] GE (*e.damping)[1]))[0]
    IF damp_ky_start GE 0 THEN BEGIN
      damp_ky = 1
      damp_ky_fac = FLTARR(par.nky0) + 1.0
      damp_ky_fac[damp_ky_start] = ((*series.ky)[damp_ky_start:*]/$
        (*series.ky)[damp_ky_start])^damp_exp
      damp_ky_fac = REBIN(REFORM(damp_ky_fac,[1,par.nky0]),$
        [par.nkx0,par.nky0,par.nz0])
    ENDIF
  ENDIF

  i = set_internal_vars(diag,{$
    n_levels    : *e.n_levels,$
    sparse_plot : *e.sparse_plot,$
    stat_dyn    : ROUND(*e.stat_dyn),$
    jpar_mult   : ROUND(*e.jpar_mult),$
    damp_kx     : damp_kx,$
    damp_ky     : damp_ky,$
    damp_kx_fac : damp_kx_fac,$
    damp_ky_fac : damp_ky_fac,$
    time_id     : PTR_NEW(),$
    A_par_id    : PTR_NEW(),$
    E_par_id    : PTR_NEW(),$
    j_par_id    : PTR_NEW(),$
    minmax_plot : minmax_plot,$
    rms_plot    : rms_plot})

  fft_format, kxky=((*i).jpar_mult ? [0,1,par.n_fields+5] : [0,1])

END

;######################################################################

PRO epar_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  ; divide A_par by B_0
  B0_inv = 1.0 / (*series.geom).Bfield
  A_par_norm = *mom[0,1].kxky * $
    REBIN(REFORM(B0_inv,[1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0])

  temp_data = COMPLEXARR(par.nx0,par.ny0,par.nz0,/NOZERO)
  IF NOT (par.nky0 MOD 2) THEN temp_data[*,par.nky0,*] = 0

  ; obtain real space B_x
  B_x_kxky = COMPLEX(0,1) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
    [par.nx0,par.nky0,par.nz0]) * A_par_norm
  IF (*i).damp_kx THEN B_x_kxky *= (*i).damp_kx_fac
  IF (*i).damp_ky THEN B_x_kxky *= (*i).damp_ky_fac
  B_x_sxky = FFT(TEMPORARY(B_x_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, par.nky0 - 1 DO temp_data[0,y,0] = B_x_sxky[*,y,*]
  FOR y = 1, par.nky0 - 1 DO $
    temp_data[0,par.ny0-y,0] = CONJ(B_x_sxky[*,y,*])
  B_x_sxky = 0
  B_x = FLOAT(FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0))

  ; obtain real space B_y
  B_y_kxky = - COMPLEX(0,1) * REBIN((*series.kx)[*,0],$
    [par.nx0,par.nky0,par.nz0]) * A_par_norm
  IF (*i).damp_kx THEN B_y_kxky *= (*i).damp_kx_fac
  IF (*i).damp_ky THEN B_y_kxky *= (*i).damp_ky_fac
  B_y_sxky = FFT(TEMPORARY(B_y_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, par.nky0 - 1 DO temp_data[0,y,0] = B_y_sxky[*,y,*]
  FOR y = 1, par.nky0 - 1 DO $
    temp_data[0,par.ny0-y,0] = CONJ(B_y_sxky[*,y,*])
  B_y_sxky = 0
  B_y = FLOAT(FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0))

  ; obtain real space E_y
  E_y_kxky = COMPLEX(0,1) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
    [par.nx0,par.nky0,par.nz0]) * (*mom[0,0].kxky)
  IF (*i).damp_kx THEN E_y_kxky *= (*i).damp_kx_fac
  IF (*i).damp_ky THEN E_y_kxky *= (*i).damp_ky_fac
  E_y_sxky = FFT(TEMPORARY(E_y_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, par.nky0 - 1 DO temp_data[0,y,0] = E_y_sxky[*,y,*]
  FOR y = 1, par.nky0 - 1 DO $
    temp_data[0,par.ny0-y,0] = CONJ(E_y_sxky[*,y,*])
  E_y_sxky = 0
  E_y = FLOAT(FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0))

  ; obtain real space E_x
  E_x_kxky = - COMPLEX(0,1) * REBIN((*series.kx)[*,0],$
    [par.nx0,par.nky0,par.nz0]) * (*mom[0,0].kxky)
  IF (*i).damp_kx THEN E_x_kxky *= (*i).damp_kx_fac
  IF (*i).damp_ky THEN E_x_kxky *= (*i).damp_ky_fac
  E_x_sxky = FFT(TEMPORARY(E_x_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, par.nky0 - 1 DO temp_data[0,y,0] = E_x_sxky[*,y,*]
  FOR y = 1, par.nky0 - 1 DO $
    temp_data[0,par.ny0-y,0] = CONJ(E_x_sxky[*,y,*])
  E_x_sxky = 0
  E_x = FLOAT(FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0))

  E_par = TEMPORARY(B_x) * TEMPORARY(E_x) + $
    TEMPORARY(B_y) * TEMPORARY(E_y)

  ; obtain real space A_par
  A_par_kxky = *mom[0,1].kxky
  IF (*i).damp_kx THEN A_par_kxky *= (*i).damp_kx_fac
  IF (*i).damp_ky THEN A_par_kxky *= (*i).damp_ky_fac
  A_par_sxky = FFT(TEMPORARY(A_par_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

  FOR y = 0, par.nky0 - 1 DO temp_data[0,y,0] = A_par_sxky[*,y,*]
  FOR y = 1, par.nky0 - 1 DO $
    temp_data[0,par.ny0-y,0] = CONJ(A_par_sxky[*,y,*])
  A_par_sxky = 0
  A_par = FLOAT(FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0))

  IF (*i).jpar_mult THEN BEGIN ; obtain real space j_par
    j_par = pf_arr([par.nx0,par.ny0,par.nz0],/zero)

    FOR n = 0, gui.out.n_spec_sel - 1 DO BEGIN
      u_par_kxky = *mom[n,par.n_fields+5].kxky
      IF (*i).damp_kx THEN u_par_kxky *= (*i).damp_kx_fac
      IF (*i).damp_ky THEN u_par_kxky *= (*i).damp_ky_fac
      u_par_sxky = FFT(TEMPORARY(u_par_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)

      FOR y = 0, par.nky0 - 1 DO temp_data[0,y,0] = u_par_sxky[*,y,*]
      FOR y = 1, par.nky0 - 1 DO $
        temp_data[0,par.ny0-y,0] = CONJ(u_par_sxky[*,y,*])
      u_par_sxky = 0
      u_par = FLOAT(FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=0))

      j_par += spec[(*gui.out.spec_select)[n]].charge * $
        spec[(*gui.out.spec_select)[n]].dens * TEMPORARY(u_par)
    ENDFOR

    (*i).j_par_id = store_step((*i).j_par_id,j_par[*,*,par.nz0/2])
  ENDIF

  (*i).A_par_id = store_step((*i).A_par_id,A_par[*,*,par.nz0/2])
  (*i).E_par_id = store_step((*i).E_par_id,E_par[*,*,par.nz0/2])
  (*i).time_id = store_step((*i).time_id,mom_time)

END

;######################################################################

PRO epar_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  A_par = store_step((*i).A_par_id,/get,/reg_array)
  E_par = store_step((*i).E_par_id,/get)
  time = store_step((*i).time_id,/get,/reg_array)
  IF (*i).jpar_mult THEN BEGIN
    j_par = store_step((*i).j_par_id,/get)
    j_spec = STRJOIN(spec[*gui.out.spec_select].name,',')
  ENDIF

  n_steps = series.step_count

  IF n_steps LT 2 THEN BEGIN
    PRINT, (*diag).name + ' error: must analyze more than one time step'
    RETURN
  ENDIF

  dApardt = FLTARR(par.nx0,par.ny0,n_steps,/NOZERO)
  dApardt[0,0,1] = (A_par[*,*,1:n_steps-1] - A_par[*,*,0:n_steps-2]) * $
    REBIN(REFORM(1.0/(time[1:n_steps-1]-time[0:n_steps-2]),$
    [1,1,n_steps-1]),[par.nx0,par.ny0,n_steps-1])
  dApardt[0,0,0] = (A_par[*,*,1] - A_par[*,*,0]) * $
    REBIN(REFORM(1.0/(time[1]-time[0]),[1]),[par.nx0,par.ny0,1])
  A_par = 0

  xaxis = (FINDGEN(par.nx0) / par.nx0 - 0.5) * par.lx
  yaxis = (FINDGEN(par.ny0) / par.ny0 - 0.5) * par.ly

  ymarg = !Y.MARGIN
  !Y.MARGIN = [ymarg[0],5]
  ysize = 17.5 * par.ly / par.lx < 25 > 12
  n_plots = (*i).stat_dyn AND NOT ((*i).minmax_plot OR (*i).rms_plot) ? 3 : 1
  IF (*i).jpar_mult THEN BEGIN
    n_plots += 1
    n_plots_x = 2
  ENDIF ELSE n_plots_x = n_plots
  set_output, diag, /ps, coltable=(((*i).minmax_plot OR (*i).rms_plot) $
    ? 41 : 47), ysize=ysize, multi=[n_plots,n_plots_x,n_plots/n_plots_x], $
    charsize=1.0

  IF (*i).minmax_plot THEN BEGIN
    step = 0L
    n_steps_sparse = CEIL(n_steps/FLOAT((*i).sparse_plot))
    min_data = FLTARR(n_steps_sparse,1+2*(*i).stat_dyn,/NOZERO)
    max_data = FLTARR(n_steps_sparse,1+2*(*i).stat_dyn,/NOZERO)
    minmax_time = FLTARR(n_steps_sparse,/NOZERO)
  ENDIF
  IF (*i).rms_plot THEN BEGIN
    step = 0L
    n_steps_sparse = CEIL(n_steps/FLOAT((*i).sparse_plot))
    rms_data = FLTARR(n_steps_sparse,1+2*(*i).stat_dyn,/NOZERO)
    rms_time = FLTARR(n_steps_sparse,/NOZERO)

    IF (*i).jpar_mult THEN rms_Epar_jpar = FLTARR(n_steps_sparse,/NOZERO)
  ENDIF

  plottitle = '!6E!D!9#!6' + $
    [(*i).stat_dyn ? ',total' : '',',flutter',',inductive'] + '!N'

  IF (*i).rms_plot THEN rms_avg_id = PTRARR(1+2*(*i).stat_dyn+(*i).jpar_mult)

  FOR n = 0, n_steps - 1 DO BEGIN
    IF n MOD (*i).sparse_plot EQ 0 THEN BEGIN
      E_par_tot = (*E_par[n]) - dApardt[*,*,n]

      FOR j = 0, 2 * (*i).stat_dyn DO BEGIN
        CASE j OF
          0    : E_data = TEMPORARY(E_par_tot)
          1    : E_data = *E_par[n]
          ELSE : E_data = - dApardt[*,*,n]
        ENDCASE

        IF (*i).minmax_plot THEN BEGIN
          min_data[step,j] = ABS(MIN(E_data,MAX=maxval))
          max_data[step,j] = ABS(TEMPORARY(maxval))
          minmax_time[step] = time[n]
          IF j EQ 2 * (*i).stat_dyn THEN step += 1
        ENDIF ELSE IF (*i).rms_plot THEN BEGIN
          rms_data[step,j] = SQRT(TOTAL(TOTAL(TOTAL((ABS(E_data))^2,1),1)*$
            (*series.geom).jac_norm)/(1L*par.nx0*par.ny0*par.nz0))
          rms_avg_id[j] = time_avg(rms_avg_id[j],rms_data[step,j],time[n])

          IF (j EQ 0) AND ((*i).jpar_mult) THEN BEGIN
            Epar_jpar = E_data * (*j_par[n])

            rms_Epar_jpar[step] = SQRT(TOTAL(TOTAL(TOTAL((ABS(Epar_jpar))^2,1),1)*$
            (*series.geom).jac_norm)/(1L*par.nx0*par.ny0*par.nz0))
            rms_avg_id[2*(*i).stat_dyn+1] = $
              time_avg(rms_avg_id[2*(*i).stat_dyn+1],rms_Epar_jpar[step],time[n])
          ENDIF
          rms_time[step] = time[n]
          IF j EQ 2 * (*i).stat_dyn THEN step += 1
        ENDIF ELSE BEGIN
          lev_col = contour_levels(E_data,(*i).n_levels)
          CONTOUR, E_data, xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
          plot_colorbar, lev_col, orientation=1

          store_colors, store_rgb, new_ct=41
          IF (*i).stat_dyn THEN !P.MULTI[0] = !P.MULTI[0] MOD n_plots + 1
          PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, $
            NOERASE=1-(*i).stat_dyn, /NODATA, /ISOTROPIC, $
            XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
            TITLE=plottitle[j]
          store_colors, store_rgb, /restore

          plot_info_str, diag, time=time[n]

          IF (j EQ 0) AND ((*i).jpar_mult) THEN BEGIN
            jpar_Epar = E_data * (*j_par[n])

            lev_col = contour_levels(jpar_Epar,(*i).n_levels)
            CONTOUR, jpar_Epar, xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
              LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
            plot_colorbar, lev_col, orientation=1

            store_colors, store_rgb, new_ct=41
            IF (*i).stat_dyn THEN !P.MULTI[0] = !P.MULTI[0] MOD n_plots + 1
            PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, $
              NOERASE=1-(*i).stat_dyn, /NODATA, /ISOTROPIC, $
              XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
              TITLE='!6E!D!9#!6,total!N!9X!6j!D!9#!6'+j_spec+'!N'
            store_colors, store_rgb, /restore

            plot_info_str, diag, time=time[n]
          ENDIF
        ENDELSE
      ENDFOR
    ENDIF
  ENDFOR

  IF (*i).minmax_plot THEN BEGIN
    FOR j = 0, 2 * (*i).stat_dyn DO BEGIN
      yrange = [MIN(min_data[*,j])<MIN(max_data[*,j]),$
        MAX(min_data[*,j])>MAX(max_data[*,j])]

      PLOT, minmax_time, min_data[*,j], /NODATA, COLOR=1, /XSTYLE, /YSTYLE, $
        XTITLE='!6'+get_var_string(0,/time,/units), YTITLE=plottitle[j], $
        YRANGE=yrange
      OPLOT, minmax_time, min_data[*,j], COLOR=4
      OPLOT, minmax_time, max_data[*,j], COLOR=2
      XYOUTS, 0.05, 0.1, '!9!!!6E!D!9#!6,min!N!9!!!6', COLOR=4, /NORMAL
      XYOUTS, 0.05, 0.15, '!9!!!6E!D!9#!6,max!N!9!!!6', COLOR=2, /NORMAL

      plot_info_str, diag
    ENDFOR
  ENDIF
  IF (*i).rms_plot THEN BEGIN
    FOR j = 0, 2 * (*i).stat_dyn DO BEGIN
      rms_avg = time_avg(rms_avg_id[j],/avg)

      yrange = [MIN(rms_data[*,j]),MAX(rms_data[*,j])]

      PLOT, rms_time, rms_data[*,j], COLOR=1, /XSTYLE, /YSTYLE, $
        XTITLE='!6'+get_var_string(0,/time,/units), YTITLE=plottitle[j], $
        YRANGE=yrange
      OPLOT, !X.CRANGE, [1,1] * rms_avg, COLOR=1, LINE=1

      PRINT, (*diag).name + ': <E_par,' + (STRSPLIT(plottitle[j],$
        '!,',/EXTRACT))[4] + '^2>^(1/2) = ' + rm0es(rms_avg,prec=4)

      plot_info_str, diag

      IF (j EQ 0) AND (*i).jpar_mult THEN BEGIN
        rms_avg = time_avg(rms_avg_id[2*(*i).stat_dyn+1],/avg)

        yrange = [MIN(rms_Epar_jpar),MAX(rms_Epar_jpar)]

        PLOT, rms_time, rms_Epar_jpar, COLOR=1, /XSTYLE, /YSTYLE, $
          XTITLE='!6'+get_var_string(0,/time,/units), $
          YTITLE='!6E!D!9#!6,total!N!9X!6j!D!9#!6'+j_spec+'!N', $
          YRANGE=yrange
        OPLOT, !X.CRANGE, [1,1] * rms_avg, COLOR=1, LINE=1

        PRINT, (*diag).name + ': <E_par,total*j_par,' + j_spec + $
          '^2>^(1/2) = ' + rm0es(rms_avg,prec=4)

        plot_info_str, diag
      ENDIF
    ENDFOR
  ENDIF

  set_output, diag, /reset
  !Y.MARGIN = ymarg

  PTR_FREE, E_par

END
