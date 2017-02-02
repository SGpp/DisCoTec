FUNCTION profile_info

  RETURN, {$
    type      : 'mom global',$
    title     : 'Profile (nonlocal)',$
    help_text : ['Plots T, n and gradient profile for each species'],$    
    ext_vars  : [$		 
    	['xind','0','index of x grid points to be plotted (integer '+$
	 'between 0 and nkx0 - 1; default: 0)']$
    	]}
 		 
END

;######################################################################

PRO profile_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars

  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = INDGEN(par.nx0)

  jac_fac = (*series.geom).jacobian

  FOR ix=0,par.nx0-1 DO $
    jac_fac[ix,*] = jac_fac[ix,*]/TOTAL(jac_fac[ix,*])
  
  i = set_internal_vars(diag,{$
    xind    : *e.xind,$
    xaxis   : par.x0 + (-series.lx/2.+(*e.xind)*series.dx)*par.rhostar,$
    rsc     : par.rhostar * par.minor_r,$
    jac_fac : jac_fac,$
    initial_prof : DBLARR(par.nx0,4,gui.out.n_spec_sel),$
    temp_id : PTR_NEW(),$
    dens_id : PTR_NEW()})

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

PRO profile_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
  n_xind = N_ELEMENTS((*i).xind)

  temp = DBLARR(n_xind,gui.out.n_spec_sel)  
  dens = DBLARR(n_xind,gui.out.n_spec_sel)
  
  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN  
    FOR iz = 0, par.nz0 - 1 DO BEGIN
      temp[*,isp] += ((*mom[isp,par.n_fields+1].sxky)[(*i).xind,0,iz]/3.+$
      	(*mom[isp,par.n_fields+2].sxky)[(*i).xind,0,iz]*2./3.)*$
      	(*i).jac_fac[*,iz]
      dens[*,isp] += (*mom[isp,par.n_fields].sxky)[(*i).xind,0,iz]*$
      	(*i).jac_fac[*,iz]
    ENDFOR
  ENDFOR

  IF (gui.misc.n_sel_diags EQ 1) THEN BEGIN
;only 0'th species shown here
     sp = (*gui.out.spec_select)[0]
     total_temp = (*spec[sp].prof).temp + temp[*,0]*(*i).rsc
     total_dens = (*spec[sp].prof).dens + dens[*,0]*(*i).rsc
     
     total_omt = -DERIV((*i).xaxis,ALOG(total_temp))/par.minor_r
     total_omn = -DERIV((*i).xaxis,ALOG(total_dens))/par.minor_r

     IF ((par.x_local) OR (par.lilo)) THEN BEGIN
        total_omt += spec[sp].omt
        total_omn += spec[sp].omn
     ENDIF

     oldmulti = !P.MULTI
     !P.MULTI = [0,2,2]
     profile_plot, diag, 0,total_temp, total_dens, total_omt, total_omn
     !P.MULTI = oldmulti
     
     XYOUTS, 0.01, 0.01, 't='+rm0es(mom_time,prec=3),/NORMAL, COLOR=1 
     WAIT, 0.1
  ENDIF
    
  (*i).temp_id = time_avg((*i).temp_id,temp,mom_time)
  (*i).dens_id = time_avg((*i).dens_id,dens,mom_time)  
      
END

;######################################################################
PRO profile_plot, diag, isp, total_temp, total_dens, total_omt, total_omn
   COMMON global_vars

   i = (*diag).internal_vars
   sp=(*gui.out.spec_select)[isp]

   ymax = MAX([total_temp, (*spec[sp].prof).temp],min=ymin)
   PLOT, (*i).xaxis, total_temp, color=1,/XSTYLE,$
     xtitle='x / a', ytitle='T(x)', YRANGE=[ymin<0,ymax],thick=2,charthick=1.42
   OPLOT,(*i).xaxis, (*spec[sp].prof).temp,color=1,linestyle=1,thick=2
   OPLOT,(*i).xaxis, (*i).initial_prof[*,0,isp],color=2,linestyle=1,thick=2
   IF (PRODUCT(!Y.CRANGE) LT 0) THEN OPLOT, !X.CRANGE, [0,0], color = 1

   ymax = MAX([total_omt, (*spec[sp].prof).omt],min=ymin)
   PLOT, (*i).xaxis, total_omt, color=1,/XSTYLE,$
     xtitle='x / a', ytitle='omt(x)', YRANGE=[ymin<0,ymax],thick=2
   OPLOT,(*i).xaxis, (*spec[sp].prof).omt,color=2,linestyle=0,thick=2
   OPLOT,(*i).xaxis, (*i).initial_prof[*,2,isp],color=2,linestyle=0,thick=2
   IF (PRODUCT(!Y.CRANGE) LT 0) THEN OPLOT, !X.CRANGE, [0,0], color = 1

   ymax = MAX([total_dens, (*spec[sp].prof).dens],min=ymin)
   PLOT, (*i).xaxis, total_dens, color=1,/XSTYLE,$
     xtitle='x / a', ytitle='n(x)', YRANGE=[ymin<0,ymax],thick=2
   OPLOT,(*i).xaxis, (*spec[sp].prof).dens,color=1,linestyle=1,thick=2
   OPLOT,(*i).xaxis, (*i).initial_prof[*,1,isp],color=2,linestyle=1,thick=2
   IF (PRODUCT(!Y.CRANGE) LT 0) THEN OPLOT, !X.CRANGE, [0,0], color = 1
 
   ymax = MAX([total_omn, (*spec[sp].prof).omn],min=ymin)   
   PLOT, (*i).xaxis, total_omn, color=1,/XSTYLE,$
     xtitle='x / a', ytitle='omn(x)', YRANGE=[ymin<0,ymax],thick=2
   OPLOT,(*i).xaxis, (*spec[sp].prof).omn,color=2,linestyle=0,thick=2 
   OPLOT,(*i).xaxis, (*i).initial_prof[*,3,isp],color=2,linestyle=0,thick=2
   IF (PRODUCT(!Y.CRANGE) LT 0) THEN OPLOT, !X.CRANGE, [0,0], color = 1
END
;######################################################################

PRO profile_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
  temp_fluc = time_avg((*i).temp_id,/avg)
  dens_fluc = time_avg((*i).dens_id,/avg)  
  
  n_xind = N_ELEMENTS((*i).xind)

  nsp = gui.out.n_spec_sel
   
  FOR isp = 0, nsp - 1 DO BEGIN
  
   sp=(*gui.out.spec_select)[isp]
   
   set_output, diag, sp, /ps, multi=[0,1,2], charsize=csize

   total_temp = (*spec[sp].prof).temp + temp_fluc[*,isp]*(*i).rsc
   total_dens = (*spec[sp].prof).dens + dens_fluc[*,isp]*(*i).rsc
   
   total_omt = -DERIV((*i).xaxis,ALOG(total_temp))/par.minor_r
   total_omn = -DERIV((*i).xaxis,ALOG(total_dens))/par.minor_r

   IF ((par.x_local) OR (par.lilo)) THEN BEGIN
      total_omt += spec[sp].omt
      total_omn += spec[sp].omn
   ENDIF

   profile_plot, diag, isp, total_temp, total_dens, total_omt, total_omn 

   XYOUTS, 0.01, 0.01, 'avg. between t='+rm0es(gui.out.start_t,prec=3)+$
     '-'+rm0es(gui.out.end_t,prec=3)+' '+get_var_string(0,/time,/ounit),$
     /NORMAL, COLOR=1 

   set_output, diag, sp, header = ['x/a','T','omt','n','omn','T_ini','omt_ini','n_ini','omn_ini'],$
     dat=[[(*i).xaxis],[total_temp],[total_omt],[total_dens],[total_omn],$
          [(*spec[sp].prof).temp[(*i).xind]],[(*spec[sp].prof).omt[(*i).xind]],$
          [(*spec[sp].prof).dens[(*i).xind]],[(*spec[sp].prof).omn[(*i).xind]]]

   set_output, diag, sp, /reset

  ENDFOR ;-- isp
END
