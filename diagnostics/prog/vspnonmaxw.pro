FUNCTION vspnonmaxw_info

  RETURN, {$
    type      : 'vsp',$
    title     : 'Non-Maxwellian v-space',$
    help_text : ['Plots contours of the Non-Maxwellian part of f_'+$
                 'as a function of v_parallel and mu'],$
    ext_vars  : [['zind','0','index of parallel slice(s); '+$
                  'default: nz0/2'],$
                 ['log','1','logarithmic z axis']]}

; NOTE: cmp. ArXiv:0907.4413v1, Fig. 7

END

;######################################################################

PRO vspnonmaxw_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = par.nz0 / 2
  IF N_ELEMENTS(*e.log) NE 1 THEN *e.log = 0

  i = set_internal_vars(diag,{$
    zind       : *e.zind,$
    log        : *e.log,$
    data_id    : PTR_NEW()})

END

;######################################################################

PRO vspnonmaxw_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  (*i).data_id = time_avg((*i).data_id,(*vsp)[(*i).zind,*,*,*,4],$
    vsp_time,fwd=gui.out.res_steps)

END

;######################################################################

PRO vspnonmaxw_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  c_levels = 40
  nzind = N_ELEMENTS((*i).zind)

  data = time_avg((*i).data_id,/avg,fwd=gui.out.res_steps,tarr=time)
  data = REFORM(data,[nzind,par.nv0,par.nw0,gui.out.n_spec_sel,$
    1+gui.out.res_steps*(series.step_count-1)],/OVERWRITE)

  vpar = - par.lv + INDGEN(par.nv0) * 2.0 * par.lv / (par.nv0 - 1.0)
  vpar_hires = - par.lv + INDGEN(500) * 2.0 * par.lv / 499.0
  GetMuWeightsAndKnots, mu_weights, mu

  ; Maxwellian: 1/pi^(3/2) * EXP(-vpar^2-mu*B)
  maxwellian = !PI^(-1.5) * EXP(REBIN(REFORM(-vpar^2,[1,par.nv0]),$
    [par.nz0,par.nv0,par.nw0])-REBIN(REFORM(REBIN(REFORM(mu,$
    [1,par.nw0]),[par.nz0,par.nw0])*REBIN((*series.geom).Bfield,$
    [par.nz0,par.nw0]),[par.nz0,1,par.nw0]),[par.nz0,par.nv0,par.nw0]))

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, coltable=33

    first = 1
    FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
      IF gui.out.res_steps THEN BEGIN
        IF time[n] GE 0.0 THEN BEGIN
          timeinfo = '!6t = ' + rm0es(time[n]) + ' ' + $
            get_var_string(1,/time,/ounit)
        ENDIF ELSE BEGIN
          gamma_omega = rm0es(get_eigenvalue(-time[n]),prec=3)
          gamom_str = '(' + gamma_omega[0] + ',' + gamma_omega[1] + ')'
          timeinfo = '!6EV ' + rm0es(-time[n]) + ', (!7c,x!6) = ' + $
            gamom_str
        ENDELSE
      ENDIF ELSE timeinfo = '!6t = ' + rm0es(gui.out.start_t) + ' to ' + $
        rm0es(gui.out.end_t) + ' ' + get_var_string(1,/time,/ounit)

      FOR z = 0, nzind - 1 DO BEGIN
        pos = [0.1,0.1+0.9/nzind*(nzind-1-z),0.925,1.0-0.9/nzind*z]
        dpos_x = pos[2] - pos[0]
        dpos_y = pos[3] - pos[1]
        dy = 0.15

        zdata = REFORM(data[z,*,*,isp,n],[par.nv0,par.nw0])
        zdmax = MAX(zdata,MIN=zdmin)
        IF zdmax NE zdmin THEN zdata = 2 * (zdata - zdmin) / $
          (zdmax - zdmin) - 1 ELSE printerror, (*diag).name + $
          ' warning: zdata contains uniform values!'
        zmwmax = MAX(maxwellian[z,*,*],MIN=zmwmin)
        zmaxwell = 2 * (maxwellian[z,*,*] - zmwmin) / (zmwmax - zmwmin) - 1
        zdata -= zmaxwell
        IF (*i).log THEN zdata = ABS(TEMPORARY(zdata))

        vperp = SQRT(mu*(*series.geom).Bfield[z])

        IF (WHERE(zdata NE 0.0))[0] NE -1 THEN BEGIN
          title = ((*i).log ? '!9#' : '') + '!6f[-1,1] - f!DM!N[-1,1]' + $
            ((*i).log ? '!9#!6' : '') + ' @ z/qR = ' + $
            rm0es((*par.z)[(*i).zind[z]],prec=3)

          lev_col = contour_levels(zdata,c_levels,log=(*i).log)
	  
          CONTOUR, zdata, vpar, vperp, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], $
            /XSTYLE, /YSTYLE, /FILL, XTITLE='!6v!D!9#!6!N', $
            YTITLE='!6v!D!9x!6!N', TITLE=title, NOERASE=(z GT 0),$
            POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,pos[0]+0.9*dpos_x,$
	    pos[3]-dy*dpos_y]
	    
          ; separation line of trapped and passing particles
          IF (par.magn_geometry EQ '') OR (par.magn_geometry EQ 's-alpha') $
            OR (par.magn_geometry EQ 's_alpha') AND par.n_pol EQ 1 THEN $
            OPLOT, vpar_hires, SQRT((*series.geom).Bfield[z]*(1.0-par.trpeps)*$
            vpar_hires^2/(par.trpeps*(1.0+COS(-!PI+(*i).zind[z]*par.dz)))), COLOR=3

          plot_colorbar, lev_col, prec=3, log=(*i).log, $
	    POSITION=[pos[0]+0.95*dpos_x,pos[1]+dy*dpos_y,pos[2],pos[3]-dy*dpos_y]

          XYOUTS, 0.05, 0.005, timeinfo, /NORMAL, CHARSIZE=1.0

          set_output, diag, sp, dat=[zdata], $
	    append=(first EQ 0), commentline=timeinfo+$
	    ', var = '+get_vsp_string(4)+' , z = '+rm0es((*par.z)[(*i).zind[z]])
        ENDIF ELSE BEGIN
          PRINT, 'vsp: empty plot'
          lev_col = contour_levels(zdata)	  
          CONTOUR, zdata, vpar, vperp, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], NOERASE=(z GT 0),$
            /XSTYLE, /YSTYLE, /NODATA, XTITLE='!6v!D!9#!6!N', $
            YTITLE='!7l!6', TITLE='z/qR = '+rm0es((*par.z)[(*i).zind[z]],prec=3)	  
        ENDELSE
	first = 0
      ENDFOR
    ENDFOR

    set_output, diag, sp, /reset
  ENDFOR ; --- isp loop 

END
