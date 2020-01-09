FUNCTION nrgphas_info

  RETURN, {$
    type      : 'nrg',$
    title     : 'nrg phase space',$
    help_text : ['nrg phase space'],$
    ext_vars  : [['vars','0','nrg vars; default: [1,6]'],$
                 ['cont','1','add a contour plot']]}
 		 
END

;######################################################################

PRO nrgphas_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.vars) LT 1 THEN *e.vars = [1,6]
  IF NOT KEYWORD_SET(*e.cont) THEN *e.cont = 0

  i = set_internal_vars(diag,{$
    vars : *e.vars,$
    cont : *e.cont})

END

;######################################################################

PRO nrgphas_loop, diag

END

;######################################################################

PRO nrgphas_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  time_nrg_full = (*nrg).time
  t_window_ind = WHERE((time_nrg_full GE gui.out.start_t) AND $
    (time_nrg_full LE gui.out.end_t))
  IF N_ELEMENTS(t_window_ind) LT 2 THEN BEGIN
    PRINT, (*diag).name + " error: insufficiently large time window"
    RETURN
  ENDIF
  time_nrg = time_nrg_full[t_window_ind]

  n_colors = 50
  carr = ROUND(FINDGEN(n_colors)/n_colors*255.0)
  dat_per_color = ROUND(N_ELEMENTS(time_nrg)/FLOAT(n_colors)+0.5)
  dat_colors = FINDGEN(dat_per_color)
  dat_index = LONARR(dat_per_color,n_colors,/NOZERO)
  FOR cind = 0, n_colors - 1 DO dat_index[*,cind] = $
    dat_colors[*] + cind * dat_per_color < N_ELEMENTS(time_nrg) - 1

  n_steps_red = -1
  IF (*i).cont THEN BEGIN
    n_steps = N_ELEMENTS(time_nrg)
    n_steps_max = 10000L
    IF n_steps GT n_steps_max THEN BEGIN
      PRINT, (*diag).name + " warning: reducing number of steps for contour"

      n_steps_red = n_steps_max
      time_nrg_cont = INTERPOL(time_nrg,n_steps_red)
    ENDIF ELSE time_nrg_cont = time_nrg
      
    grid_xres = 200
    grid_yres = 200
    xy_grid = LONARR(grid_xres,grid_yres)
  ENDIF

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, coltable=33, multi=[0,1,2+(*i).cont]

    FOR vind = 0, N_ELEMENTS((*i).vars) - 1 DO BEGIN
      var = (*i).vars[vind]
      data_nrg = (*nrg)[t_window_ind].data[sp*8+var]
      deriv_nrg = DERIV(time_nrg,data_nrg)
      IF n_steps_red NE -1 THEN $
        data_nrg_cont = INTERPOL(data_nrg,n_steps_red) ELSE $
        data_nrg_cont = data_nrg
      deriv_nrg_cont = DERIV(time_nrg_cont,data_nrg_cont)

      IF (*i).cont THEN BEGIN
        data_min = MIN(data_nrg_cont,MAX=data_max)
        deriv_min = MIN(deriv_nrg_cont,MAX=deriv_max)
        IF (ABS(data_max-data_min) LT 1e-8) OR (ABS(deriv_max-deriv_min) LT 1e-8) $
          THEN BEGIN
          PRINT, (*diag).name + " error: constant data or derivative"
          RETURN
        ENDIF
        data_range_inv = 1.0 / (data_max - data_min)
        deriv_range_inv = 1.0 / (deriv_max - deriv_min)

        xy_grid[*,*] = 0
        FOR t = 0, N_ELEMENTS(time_nrg_cont) - 1 DO $
          xy_grid[ROUND((data_nrg_cont[t]-data_min)*data_range_inv*(grid_xres-1)),$
          ROUND((deriv_nrg_cont[t]-deriv_min)*deriv_range_inv*(grid_yres-1))] += 1

        lev_col = contour_levels([MIN(xy_grid),MAX(xy_grid)],n_levels)
        CONTOUR, xy_grid, FINDGEN(grid_xres) / (grid_xres * data_range_inv) + data_min, $
          FINDGEN(grid_yres) / (grid_yres * deriv_range_inv) + deriv_min, /FILL, $
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], COLOR=1, $
          /XSTYLE, /YSTYLE, XTITLE="!6"+get_nrg_string(var,/fancy), $
          YTITLE="!6d"+get_nrg_string(var,/fancy)+"!6/dt"
      ENDIF

      PLOT, [MIN(data_nrg),MAX(data_nrg)], [MIN(deriv_nrg),MAX(deriv_nrg)], $
        /NODATA, COLOR=1, XTITLE="!6"+get_nrg_string(var,/fancy), $
        YTITLE="!6d"+get_nrg_string(var,/fancy)+"!6/dt"
      FOR cind = 0, n_colors - 1 DO PLOTS, data_nrg[dat_index[*,cind]], $
        deriv_nrg[dat_index[*,cind]], COLOR=carr[cind], PSYM=3

      PLOT, [MIN(time_nrg),MAX(time_nrg)], [MIN(data_nrg),MAX(data_nrg)], $
        /NODATA, COLOR=1, XTITLE="!6"+get_var_string(0,/time,/fancy), $
        YTITLE="!6"+get_nrg_string(var,/fancy), /XSTYLE
      FOR cind = 0, n_colors - 1 DO OPLOT, time_nrg[dat_index[*,cind]], $
        data_nrg[dat_index[*,cind]], COLOR=carr[cind]

      set_output, diag, sp, suffix=get_nrg_string(var), dat=[data_nrg,deriv_nrg]
    ENDFOR

    set_output, diag, sp, /reset
  ENDFOR

END
