FUNCTION betaprofile_info

  RETURN, {$
    type      : 'mom global',$
    title     : 'Beta profile (nonlocal)',$
    help_text : ['Plots beta and beta_crit profile'],$    
    ext_vars  : [$		 
    	['xind','0','index of x grid points to be plotted (integer '+$
	 'between 0 and nkx0 - 1; default: 0)']$
    	]}
 		 
END

;######################################################################

PRO betaprofile_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars

  IF (par.n_spec NE gui.out.n_spec_sel) THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + $
      ': requires all species to be selected'
    (*diag).selected = 0
    RETURN
  ENDIF

  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = INDGEN(par.nx0)

  jac_fac = (*series.geom).jacobian

  FOR ix=0,par.nx0-1 DO $
    jac_fac[ix,*] = jac_fac[ix,*]/TOTAL(jac_fac[ix,*])

  x_a = par.x0 + par.rhostar * (INDGEN(par.nx0) * par.dx - 0.5 * par.lx)
  shat = x_a / (*series.geom).q * DERIV(x_a,(*series.geom).q)
  eps = x_a * par.minor_r/par.major_R

  bfield_flxavg = TOTAL((*series.geom).Bfield*jac_fac,2)

  i = set_internal_vars(diag,{$
    xind    : *e.xind,$
    xaxis   : par.x0 + (-series.lx/2.+(*e.xind)*series.dx)*par.rhostar,$
    rsc     : par.rhostar * par.minor_r,$
    bcrit_fac : 0.6*shat/(*series.geom).q^2/(1.+eps),$
    beta_fac : par.beta/2 / bfield_flxavg,$  ; Factor of 2 ???
    jac_fac : jac_fac,$
    initial_prof : DBLARR(par.nx0,4,gui.out.n_spec_sel),$
    beta_id : PTR_NEW(),$
    beta_crit_id : PTR_NEW()})

  IF (*i).rsc EQ 0.0 THEN printerror, 'Warning: rescaling of fluctuations equals zero'
    
  fft_format, sxky=[par.n_fields,par.n_fields+1,par.n_fields+2]

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN  
    sp=(*gui.out.spec_select)[isp]
    (*i).initial_prof[*,0,isp] = (*spec[sp].prof).temp
    (*i).initial_prof[*,1,isp] = (*spec[sp].prof).dens
    (*i).initial_prof[*,2,isp] = (*spec[sp].prof).omt
    (*i).initial_prof[*,3,isp] = (*spec[sp].prof).omn
  ENDFOR

END

;######################################################################

PRO betaprofile_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
  n_xind = N_ELEMENTS((*i).xind)

  temp = DBLARR(n_xind,/NOZERO)
  dens = DBLARR(n_xind,/NOZERO)

  pressure = DBLARR(n_xind)
  grad_sum = DBLARR(n_xind)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN  
    sp = (*gui.out.spec_select)[isp]
    temp = 0.0
    dens = 0.0
    FOR iz = 0, par.nz0 - 1 DO BEGIN
      temp += DOUBLE((*mom[isp,par.n_fields+1].sxky)[(*i).xind,0,iz]/3.+$
      	(*mom[isp,par.n_fields+2].sxky)[(*i).xind,0,iz]*2./3.)*$
      	(*i).jac_fac[*,iz]
      dens += DOUBLE((*mom[isp,par.n_fields].sxky)[(*i).xind,0,iz])*$
      	(*i).jac_fac[*,iz]
    ENDFOR

;TO DO: need to check for spec[sp]%temp and dens factors!!
    temp = (*spec[sp].prof).temp + temp*(*i).rsc
    dens = (*spec[sp].prof).dens + dens*(*i).rsc

    omt = -DERIV((*i).xaxis,ALOG(temp))/par.minor_r
    omn = -DERIV((*i).xaxis,ALOG(dens))/par.minor_r

    pressure += dens*temp
    grad_sum += (omn+omt)
  ENDFOR


  beta_crit = (*i).bcrit_fac[(*i).xind]/grad_sum
  beta = (*i).beta_fac[(*i).xind] * pressure

  IF (gui.misc.n_sel_diags EQ 1) THEN BEGIN
     oldmulti = !P.MULTI
     !P.MULTI = [0,1,1]
     betaprofile_plot, diag, beta, beta_crit
     !P.MULTI = oldmulti
 
     XYOUTS, 0.01, 0.01, 't='+rm0es(mom_time,prec=3),/NORMAL, COLOR=1 
     WAIT, 0.1
  ENDIF
    
  (*i).beta_id = time_avg((*i).beta_id,beta,mom_time)
  (*i).beta_crit_id = time_avg((*i).beta_crit_id,beta_crit,mom_time)  
      
END

;######################################################################
PRO betaprofile_plot, diag, beta, beta_crit
   COMMON global_vars

   i = (*diag).internal_vars

   ymax = 1.5*MAX(beta,min=ymin)
   ymin *= 1.5

   PLOT, (*i).xaxis, beta, color=1,/XSTYLE,$
     xtitle='x / a', ytitle='beta(x), beta_crit', YRANGE=[ymin<0,ymax]
   OPLOT,(*i).xaxis, beta_crit,color=2,linestyle=1
   IF (PRODUCT(!Y.CRANGE) LT 0) THEN OPLOT, !X.CRANGE, [0,0], color = 1
  
END
;######################################################################

PRO betaprofile_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  beta = time_avg((*i).beta_id,/avg)
  beta_crit = time_avg((*i).beta_crit_id,/avg)  
  
  n_xind = N_ELEMENTS((*i).xind)

  set_output, diag, /eps, multi=[0,1,1], charsize=csize

  betaprofile_plot, diag, beta, beta_crit
  
  XYOUTS, 0.01, 0.01, 'avg. between t='+rm0es(gui.out.start_t,prec=3)+$
     '-'+rm0es(gui.out.end_t,prec=3)+' '+get_var_string(0,/time,/ounit),$
     /NORMAL, COLOR=1 

  set_output, diag, sp, header = ['x/a','beta','beta_crit'],$
       dat=[[(*i).xaxis],[beta],[beta_crit]]
  
  set_output, diag, sp, /reset, /eps


END
