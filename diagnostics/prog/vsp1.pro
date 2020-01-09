FUNCTION vsp1_info

  RETURN, {$
    type      : 'vsp',$
    title     : 'velocity space',$
    help_text : ['vsp diag'],$
    ext_vars  : [['var','0','variable to be plotted '+$
    	    	  '0 = G_es, 1 = G_em, 2 = Q_es, 3 = Q_em, 4 = abs(f_); '+$
		  'default: 4'],$
    	    	['zind','0','index of parallel slices to be plotted; '+$
		 'default: [0,nz0/2,nz0-1]'],$
		['log','1','switch on/off logarithmic color encoding']]}

END

;######################################################################

PRO vsp1_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars
  
  IF N_ELEMENTS(*e.var) NE 1 THEN *e.var = 4 ; 0 = G_es, 1 = G_em, 2 = Q_es, 3 = Q_em, 4 = abs(f_)
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = [0,par.nz0/2,par.nz0-1]
  IF N_ELEMENTS(*e.log) NE 1 THEN *e.log = 0

  i = set_internal_vars(diag,{$
    var        : *e.var,$
    zind       : *e.zind,$
    log        : *e.log,$
    time_id    : PTR_NEW(),$
    data_id    : PTR_NEW()})

END

;######################################################################

PRO vsp1_loop, diag

  COMMON global_vars
  
  i = (*diag).internal_vars  

  (*i).data_id =  time_avg((*i).data_id,(*vsp)[(*i).zind,*,*,*,$
        (*i).var],vsp_time,fwd=gui.out.res_steps)
      
END

;######################################################################

PRO vsp1_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nzind = N_ELEMENTS((*i).zind)
  nsp = gui.out.n_spec_sel
  nvar = N_ELEMENTS((*i).var)
  n_res_steps = gui.out.res_steps * (series.step_count - 1) + 1

  data = REFORM(time_avg((*i).data_id,/avg,fwd=gui.out.res_steps,$
    tarr=time),[nzind,par.nv0,par.nw0,nsp,nvar,n_res_steps])

  vpar = - par.lv + INDGEN(par.nv0) * 2.0 * par.lv / (par.nv0 - 1.0)
  GetMuWeightsAndKnots, mu_weights, mu
  
  c_levels = 40


  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]
    
    set_output, diag, sp, /ps, coltable=33, multi=[0,1,1]

    first = 1
    FOR t = 0,  n_res_steps - 1 DO BEGIN
      FOR iz = 0, nzind - 1 DO BEGIN
        pos = [0.1,0.1+0.9/nzind*(nzind-1-iz),0.925,1.0-0.9/nzind*iz]
        dpos_x = pos[2] - pos[0]
        dpos_y = pos[3] - pos[1]
        dy = 0.15
      
        IF (WHERE(data[iz,*,*,isp,0,t] NE 0.0))[0] NE -1 THEN BEGIN
          lev_col = contour_levels(data[iz,*,*,isp,0,t],c_levels,log=(*i).log)
	  
          CONTOUR, data[iz,*,*,isp,0,t], vpar, mu, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], $
            /XSTYLE, /YSTYLE, /FILL, XTITLE='!6v!I!9#!6!N', $
            YTITLE='!7l!6', TITLE=get_vsp_string((*i).var,/fancy)+' @ z/qR = '+$
            rm0es((*par.z)[(*i).zind[iz]],prec=3), NOERASE=(iz GT 0),$
            POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,pos[0]+0.9*dpos_x,$
	    pos[3]-dy*dpos_y]
          
          plot_trapped_passing_bnd, (*i).zind[iz], vpar, mu
;          plot_j0_minmax, sp, (*i).zind[iz], mu

          plot_colorbar, lev_col, prec=3, log=(*i).log, $
            POSITION=[pos[0]+0.95*dpos_x,pos[1]+dy*dpos_y,pos[2],pos[3]-dy*dpos_y]

        ENDIF ELSE BEGIN
           PRINT, 'vsp: empty plot'
           lev_col = contour_levels(data[iz,*,*,isp,0,t])	  
           CONTOUR, data[iz,*,*,isp,0,t], vpar, mu, $
             LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], NOERASE=(iz GT 0),$
             /XSTYLE, /YSTYLE, /NODATA, XTITLE='!6v!I!9#!6!N', $
             YTITLE='!7l!6', TITLE='z/qR = '+rm0es((*par.z)[(*i).zind[iz]],prec=3),$
             POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,pos[0]+0.9*dpos_x,$
	     pos[3]-dy*dpos_y]
        ENDELSE

        plot_info_str, diag, time=(gui.out.res_steps ? time[t] : undef)

        set_output, diag, sp, dat=[[REFORM(data[iz,*,*,isp,0,t])]], $
	   append=(first EQ 0), commentline=(gui.out.res_steps ? $
           ', time = '+rm0es(time[t]) : 't avg.')+$
	   ', var = '+get_vsp_string((*i).var)+' , z = '+rm0es((*par.z)[(*i).zind[iz]])
     ENDFOR ;iz
   ENDFOR ;t

    set_output, diag, sp, /reset
  ENDFOR ; --- isp loop 


END

;################################################################################

PRO plot_trapped_passing_bnd, zind, vpar, mu, Blev=Blev
  COMMON global_vars

 

  vpar_int = INTERPOL(vpar,512>par.nv0)
  mu_int = mu ;INTERPOL(mu,64>par.nw0)

  show_grid = 0
  IF (show_grid) THEN BEGIN
    FOR ivp=0, N_ELEMENTS(vpar)-1 DO BEGIN
      FOR imu = 0, N_ELEMENTS(mu)-1 DO BEGIN
         PLOTS, vpar[ivp], mu[imu], PSYM=4, color=10, SYMSIZE=0.5
      ENDFOR
    ENDFOR
  ENDIF
  
  nbf = 2-(par.x_local AND par.y_local)
  lBfield = FLTARR(nbf,/NOZERO) ;local Bfield
  lBlev= FLTARR(nbf,/NOZERO)    ;local Blev
  IF (par.x_local AND par.y_local) THEN BEGIN
     lBfield = (*series.geom).Bfield[zind]
     lBlev = max((*series.geom).Bfield)
     maxrat = lBfield/lBlev
  ENDIF ELSE BEGIN
     ;evaluate min/max trapped/passing boundary for x/y global
     minrat = 100
     maxrat = -1
     FOR i = 0, N_ELEMENTS((*series.geom).Bfield[*,0])-1 DO BEGIN
        ratio = (*series.geom).Bfield[i,zind]/MAX((*series.geom).Bfield[i,*])
        IF (ratio LT minrat) THEN BEGIN
           minrat = ratio
           lBlev[0] = MAX((*series.geom).Bfield[i,*])
           lBfield[0] = (*series.geom).Bfield[i,zind]
        ENDIF
        IF (ratio GT maxrat) THEN BEGIN
           maxrat = ratio
           lBlev[1] = MAX((*series.geom).Bfield[i,*])
           lBfield[1] = (*series.geom).Bfield[i,zind]
        ENDIF
     ENDFOR
  ENDELSE

  IF KEYWORD_SET(Blev) THEN lBlev = Blev

  FOR bf = 0, nbf-1 DO BEGIN
     FOR imu = 0, N_ELEMENTS(mu_int)-1 DO BEGIN
        energy_old= vpar_int[0]^2 + mu_int[imu]*lBfield[bf]
        FOR ivp = 1, N_ELEMENTS(vpar_int)-1 DO BEGIN
           energy = vpar_int[ivp]^2 + mu_int[imu]*lBfield[bf]
           IF ((energy_old-mu_int[imu]*lBlev[bf])*(energy-mu_int[imu]*lBlev[bf]) LE 0) THEN $
              PLOTS, vpar_int[ivp], mu_int[imu], PSYM=4, color=4, SYMSIZE=0.5
           energy_old = energy
        ENDFOR
        IF (maxrat EQ 1.0) THEN $ ;plot vertical 'line' for zero trapping
           PLOTS,0., mu_int[imu], PSYM=4, color=4, SYMSIZE=0.5
     ENDFOR
  ENDFOR


  ; separation line of trapped and passing particles (analytic result)
;  IF (par.magn_geometry EQ '') OR (par.magn_geometry EQ 's-alpha') $
;     OR (par.magn_geometry EQ 's_alpha') AND par.n_pol EQ 1 THEN $
;        OPLOT, vpar_int, (1.0 - par.trpeps) * vpar_int^2 / (par.trpeps * $
;        (1.0 + COS(-!PI+zind*par.dz))), COLOR=3
END

PRO plot_j0_minmax, sp, zind, mu
  COMMON global_vars

  N = 64>par.nw0
  mu_int = INTERPOL(mu,N)
  j0 = FLTARR(N,/NOZERO)

  k_perp = par.kx_center

  FOR imu=0, N-1 DO BEGIN
     rho = sqrt(2.0*spec[sp].mass*spec[sp].temp*mu_int[imu]/$
                spec[sp].charge^2/(*series.geom).Bfield[zind])
     j0[imu] = BESELJ(k_perp*rho,0)
  ENDFOR
  OPLOT, j0, mu_int, color=127
END
