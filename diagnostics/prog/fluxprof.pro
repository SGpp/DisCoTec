FUNCTION fluxprof_info

  RETURN, {$
    type      : 'mom global',$
    title     : 'Flux Profile (nonlocal)',$
    help_text : ['Plot radial profiles of particle and heat fluxes'],$    
    ext_vars  : [$		 
    	['xind','0','indices of x grid points to be plotted (integer '+$
	 'between 0 and nx0 - 1; default: 0:nx0-1)'],$
        ['norm_projection','1','switch on/off radial projection '+$
         'using the normalized contravariant radial unit vector ' +$
         '(grad x is taken otherwise)']$
    	]}
 		 
END

;######################################################################

PRO fluxprof_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars

  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = INDGEN(par.nx0)
  n_xi  = n_elements(*e.xind)
  IF NOT KEYWORD_SET(*e.norm_projection) THEN *e.norm_projection = 0


  IF (par.x_local OR par.lilo) THEN xgrid = -series.lx/2.+(*e.xind)*series.dx $
  ELSE xgrid = par.x0 + (-series.lx/2.+(*e.xind)*series.dx)*par.rhostar

  nf = par.n_fields
  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  IF with_em THEN BEGIN
    n_fluxes = 4*(par.n_moms/6) ; G_es, Q_es, G_em, Q_em
    vars = INDGEN(par.n_fields+par.n_moms)
  ENDIF ELSE BEGIN
    n_fluxes = 2*(par.n_moms/6) ; G_es, Q_es (passing, trapped, FLR)
    IF (par.n_moms EQ 6) THEN vars = [0,nf,nf+1,nf+2] ELSE BEGIN $
      vars = 0
      FOR j = 0, par.n_moms/6-1 DO vars = [vars,nf+6*j,nf+1+6*j,nf+2+6*j]
    ENDELSE
  ENDELSE   
  
  IF (*e.norm_projection) THEN BEGIN
     proj_fac = par.x_local?1.0/REBIN(REFORM(sqrt((*series.geom).gxx),$
         [1,1,par.nz0]),[par.nkx0,par.nky0,par.nz0]):$
          1.0/REBIN(REFORM(sqrt((*series.geom).gxx),$
         [par.nx0,1,par.nz0]),[par.nx0,par.nky0,par.nz0])
  ENDIF ELSE proj_fac = 1.0
  
  proj_fac *= par.x_local?1.0/(*series.geom).C_xy:$
              1.0/REBIN(REFORM((*series.geom).C_xy,$
              [par.nx0,1,1]),[par.nkx0,par.nky0,par.nz0])

  jac_fac = (*series.geom).jacobian
  area = (2.0*!PI)^2*par.n_pol
  IF (par.x_local) THEN BEGIN
    jac_fac[*] = jac_fac[*]/TOTAL(jac_fac)
    IF (*e.norm_projection) THEN $
       area *= TOTAL((*series.geom).jacobian*sqrt((*series.geom).gxx))/par.nz0 $
    ELSE area *= TOTAL((*series.geom).jacobian)/par.nz0  ;use dVdx instead
  ENDIF ELSE BEGIN
    FOR ix=0,par.nx0-1 DO $ 
      jac_fac[ix,*] = jac_fac[ix,*]/TOTAL(jac_fac[ix,*])
    IF (*e.norm_projection) THEN area *= total((*series.geom).jacobian[*,*]*$
            sqrt((*series.geom).gxx[*,*]),2)/par.nz0 $
    ELSE area *= total((*series.geom).jacobian[*,*],2)/par.nz0 ;use dVdx instead
  ENDELSE
  area *= (*series.geom).C_y

  ; for comparison with nrg (volume average):
  ;(*e.xind) = INDGEN(par.nx0)  
  ;n_xi=par.nx0
  ;jac_fac = (*series.geom).jacobian/TOTAL((*series.geom).jacobian)*n_xi    
  
  i = set_internal_vars(diag,{$
    kyind       : (*series.ky),$
    xind        : *e.xind,$
    norm_proj   : *e.norm_projection,$
    xgrid       : xgrid ,$                         
    k_mult      : 0 ,$     ; not using k_mult for the moment
    n_fluxes    : n_fluxes, $
    jac_fac     : jac_fac,$
    proj_fac    : proj_fac,$
    area        : area,$	
    profx       : DBLARR(n_xi,n_fluxes,/NOZERO),$ ; defined here to improve performance
    tmpprofx    : DBLARR(n_xi,par.nz0,/NOZERO),$ ; defined here to improve performance     
    spectrumy   : DBLARR(par.nky0,n_fluxes,/NOZERO),$ ; defined here to improve performance
    tempspeky   : DBLARR(par.nky0,par.nz0,/NOZERO),$ ; defined here to improve performance
    avprofx_id  : PTR_NEW(),$
    avspecty_id : PTR_NEW()}) 

  fft_format, sxky=vars

END

;######################################################################

;######################################################################

FUNCTION fluxprof_calc_x_prof, diag, var1, var2

  COMMON global_vars

  i = (*diag).internal_vars

  n_xi = n_elements((*i).xind) ; to make source code more concise
  nx0 = par.nx0
  nky0 = par.nky0
  
  FOR nx = 0, n_xi-1 DO BEGIN
      (*i).tmpprofx[nx,*]=$
       DOUBLE(var1[(*i).xind[nx],0]*CONJ(var2[(*i).xind[nx],0])) + $
        2.0D*TOTAL(DOUBLE(CONJ(var1[(*i).xind[nx],1:nky0-1,*])*var2[(*i).xind[nx],1:nky0-1,*]),2,/DOUBLE)
  ENDFOR

  IF par.x_local THEN BEGIN
    FOR nzi = 0, par.nz0 - 1 DO (*i).tmpprofx[*,nzi] = $
      (*i).tmpprofx[*,nzi] * (*i).jac_fac[nzi]
  ENDIF ELSE BEGIN
    FOR nzi = 0, par.nz0 - 1 DO (*i).tmpprofx[*,nzi] = $
      (*i).tmpprofx[*,nzi] * (*i).jac_fac[(*i).xind,nzi]      
  ENDELSE

  RETURN, TOTAL((*i).tmpprofx,2,/DOUBLE)

END

;######################################################################

FUNCTION fluxprof_calc_y_spect, diag, var1, var2

  COMMON global_vars

  i = (*diag).internal_vars

  n_xi = n_elements((*i).xind) ; to make source code more concise
  nky0 = par.nky0
  lx = par.lx
  
  ximin=(*i).xind[0]
  ximax=(*i).xind[n_xi-1]
  
  tempspeky = DCOMPLEXARR(par.nky0,/NOZERO)
  IF (par.x_local) THEN BEGIN
    FOR iky = 0, nky0 - 1 DO tempspeky[iky] = $
      TOTAL(DOUBLE(CONJ(var1[ximin:ximax,iky,*])*var2[ximin:ximax,iky,*])*$
            (*i).jac_fac[*])	
  ENDIF ELSE BEGIN
    FOR iky = 0, nky0 - 1 DO tempspeky[iky] = $
      TOTAL(DOUBLE(CONJ(var1[ximin:ximax,iky,*])*var2[ximin:ximax,iky,*])*$
            (*i).jac_fac[ximin:ximax,*])      
  ENDELSE

  RETURN, 2.0D*tempspeky / n_xi

END

;######################################################################

PRO fluxprof_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
  n_xi = n_elements((*i).xind) 
  
  var1 = DCOMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO)

  avprofx =DBLARR(n_xi,(*i).n_fluxes,gui.out.n_spec_sel,/NOZERO)
  avspecty=DBLARR(par.nky0,(*i).n_fluxes,gui.out.n_spec_sel,/NOZERO)

  nf = par.n_fields
  rho_ref = sqrt(series.mref*series.Tref)/(series.Qref*series.Bref)
  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)

  ; Calculate vx, Bx, grad q_par
  ve_x = DCOMPLEXARR(par.nx0,par.nky0,par.nz0,/NOZERO)
  IF with_em THEN B_x = DCOMPLEXARR(par.nx0,par.nky0,par.nz0,/NOZERO)

  ; ve_x = - d phi / dy
  ve_x = -DCOMPLEX(0.0,1.0)*REBIN(REFORM((*series.ky),$
         [1,par.nky0,1]),[par.nkx0,par.nky0,par.nz0])* $
          (*mom[0,0].sxky) / series.Bref * (*i).proj_fac
  ; B_x = dApar/dy
  IF with_em THEN $
     B_x = DCOMPLEX(0.0,1.0)*REBIN(REFORM((*series.ky),$
       [1,par.nky0,1]),[par.nkx0,par.nky0,par.nz0])* $
        (*mom[0,1].sxky) / series.Bref * (*i).proj_fac

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN  
    sp=(*gui.out.spec_select)[isp]
    n_notrap = 2+2*with_em ; number of fluxes without trapdiag
    FOR j = 0,par.n_moms/6-1 DO BEGIN ; passing, trapped, FLR
      ; Gamma_es = <n * ve_x>
       (*i).profx[*,n_notrap*j+0]    =fluxprof_calc_x_prof(diag, $
          (*mom[isp,nf+6*j].sxky), ve_x)*spec[sp].dens
           ;
       (*i).spectrumy[*,n_notrap*j+0]=fluxprof_calc_y_spect(diag, $
         (*mom[isp,nf+6*j].sxky), ve_x)*spec[sp].dens
	
       ; Q_es = (1/2 Tpar + Tperp + 3/2 n) ve_x
       IF par.x_local THEN var1 = ((0.5*(*mom[isp,nf+6*j+1].sxky)+$
          (*mom[isp,nf+6*j+2].sxky))/series.nref+$
	  1.5*(*mom[isp,nf+6*j].sxky)/(series.Tref*series.Qref)) $
       ELSE BEGIN 
          FOR ix=0,par.nkx0-1 DO BEGIN
             var1[ix,*,*] = ((0.5*(*mom[isp,nf+6*j+1].sxky)[ix,*,*]+(*mom[isp,nf+6*j+2].sxky)[ix,*,*])*$
	      (*spec[sp].prof).dens[ix]/series.nref+1.5*(*mom[isp,nf+6*j].sxky)[ix,*,*]*$
              (*spec[sp].prof).temp[ix]/(series.Tref*series.Qref))
          ENDFOR
       ENDELSE
       (*i).profx[*,n_notrap*j+1] = fluxprof_calc_x_prof(diag,var1,ve_x)$
	    *spec[sp].dens*spec[sp].temp
       ;
       (*i).spectrumy[*,n_notrap*j+1] = fluxprof_calc_y_spect(diag,var1,ve_x)$
	    *spec[sp].dens*spec[sp].temp
	 
       IF with_em THEN BEGIN
          IF par.x_local THEN var1=(*mom[isp,nf+6*j+5].sxky) $
          ELSE FOR ix=0,par.nkx0-1 DO $
              var1[ix,*,*] = ((*mom[isp,nf+6*j+5].sxky)[ix,*,*]*(*spec[sp].prof).dens[ix])
          ; Gamma_em = <upar * B_x>
          (*i).spectrumy[*,n_notrap*j+2] = fluxprof_calc_y_spect(diag, $
             var1,B_x)*spec[sp].dens
          (*i).profx[*,n_notrap*j+2] = fluxprof_calc_x_prof(diag, $
             var1,B_x)*spec[sp].dens
          ; Q_em = (qpar + qperp + 5/2 upar) B_x
          ; mom(5) und mom(6) are NOT identical to qpar and qperp, but qpar+1.5upar and qperp+upar
          (*i).spectrumy[*,n_notrap*j+3] = fluxprof_calc_y_spect(diag, $
           ((*mom[isp,nf+6*j+3].sxky)+(*mom[isp,nf+6*j+4].sxky)), B_x)$
             *spec[sp].dens*spec[sp].temp
          (*i).profx[*,n_notrap*j+3] = fluxprof_calc_x_prof(diag, $
           ((*mom[isp,nf+6*j+3].sxky)+(*mom[isp,nf+6*j+4].sxky)), B_x)$
             *spec[sp].dens*spec[sp].temp
       ENDIF ; em effect
	           
    ENDFOR
    avprofx[*,*,isp]  =  (*i).profx
    avspecty[*,*,isp] =  (*i).spectrumy
    
    output = NOT !QUIET         ; switch on for comparison with nrg - file
    IF output THEN BEGIN 
       print, '--- fluxprof ---'
       G_es = 0 & G_em = 0 & Q_es = 0 & Q_em = 0  
       FOR j= 0, par.n_moms/6-1 DO BEGIN
          G_es = G_es + TOTAL((*i).spectrumy[*,0+n_notrap*j],1)
          Q_es = Q_es + TOTAL((*i).spectrumy[*,1+n_notrap*j],1)
          IF with_em THEN BEGIN
             G_em = G_em + TOTAL((*i).spectrumy[*,2+n_notrap*j],1)
             Q_em = Q_em + TOTAL((*i).spectrumy[*,3+n_notrap*j],1)
          ENDIF      
       ENDFOR 
       print, spec[sp].name
       print, 'G_es = ', string(G_es,format='(E12.4)'), $
              string(total((*i).profx[*,0],1)/par.nx0,format='(E12.4)')
       print, 'Q_es = ', string(Q_es, format='(E12.4)'), $
              string(total((*i).profx[*,1],1)/par.nx0,format='(E12.4)')
       IF with_em THEN BEGIN
          print, 'G_em = ', string(G_em, format='(E12.4)'), $
                 string(total((*i).profx[*,2],1)/par.nx0,format='(E12.4)')
          print, 'Q_em = ', string(Q_em, format='(E12.4)'), $
                 string(total((*i).profx[*,3],1)/par.nx0,format='(E12.4)')
       ENDIF
    ENDIF                       ; output    
    
    IF (gui.misc.n_sel_diags EQ 1) THEN BEGIN
       fluxprof_plot, (*i).xgrid,avprofx[*,*,isp], '!6x/a',$
                      'v', rm0es(mom_time), plot_type=1,/no_legend
       WAIT, 0.05
    ENDIF
 ENDFOR                         ; species loop
  
  (*i).avprofx_id  = time_avg((*i).avprofx_id,avprofx,mom_time)
  (*i).avspecty_id = time_avg((*i).avspecty_id,avspecty,mom_time)
  
  
END

;######################################################################

PRO fluxprof_plot, xind, arrin, xtitle, var_names, timeline, $
  titleline=titleline, k_mult=k_mult, plot_type=plot_type,$
  no_legend=no_legend

 ; plot_type = 1 -> x profile
 ;             2 -> ky spectra
  
  np = n_elements(arrin[0,*])
  maxy = max(arrin, min = miny)
  colorarr = [1,3,2,4,5,6] ; to be color consistent with former plots
  IF NOT KEYWORD_SET(titleline) THEN titleline=''
  IF N_ELEMENTS(k_mult) EQ 0 THEN k_mult = 0

  
  n_start=1
  if (plot_type eq 1) then n_start=3
  FOR j = n_start,3 DO BEGIN ; log-log/log-lin/lin-lin - Plot
    xrange = [xind[0],xind[N_ELEMENTS(xind)-1]]
    
    IF (1/j) THEN BEGIN
      ytickformat='betterticks'
      yrange = [1e-5, maxy]
    ENDIF ELSE BEGIN
      ytickformat=''
      yrange=[miny, maxy]
    ENDELSE

    IF FIX(j/3+1) THEN BEGIN
      xtickformat='betterticks' 
      IF xind[0] EQ 0 THEN xrange[0]=xind[1]
    ENDIF ELSE xtickformat=''
    
    plot, xind, arrin[*,0],color=colorarr[0],xtitle=xtitle,$
       xlog=j/3+1,xtickformat=xtickformat,xstyle=1,$
       xticklen=1.0,position=[0.15,0.4,0.9,0.9],$
       ylog=1/j,ytickformat=ytickformat,yrange=yrange,$
       title=titleline, xrange=xrange
       
    FOR var=1, np-1 DO $
      oplot, xind, arrin[*,var],color=colorarr[var]      
      
    IF k_mult AND (2/j) THEN BEGIN
      FOR var=0, np-1 DO BEGIN
        scalfac = max(ABS(arrin[*,var]))/$
	  max(ABS(xind*arrin[*,var]))
        oplot, xind, scalfac*xind*arrin[*,var],$
          color=colorarr[var],linestyle=2
      ENDFOR
    ENDIF              

    IF (j EQ 2) THEN PLOTS, 10^!X.CRANGE, [0,0], COLOR=1
    IF (j EQ 3) THEN PLOTS, !X.CRANGE, [0,0], COLOR=1
   
    ; Laufdaten ausgeben (eine Zeile)
    xyouts,0.15,0.225,/normal, timeline,color=1

    ; Legende ausgeben
    IF NOT KEYWORD_SET(no_legend) THEN $
      plot_legend, colorarr[0:np-1], var_names[0:np-1], $
        per_line = 2
  ENDFOR
END


;######################################################################

PRO fluxprof_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).norm_proj THEN proj_str = 'A' $
  ELSE proj_str = 'dVdx'

  v_ref = SQRT(par.Tref*par.Qref/par.mref)
  rho_ref = SQRT(par.mref*par.Tref/par.Qref) / par.Bref
  D_gb = v_ref * rho_ref^2 / par.Lref
  G_gb = D_gb*par.nref/par.Lref
  Q_gb = D_gb*par.nref*par.Tref*par.Qref/par.Lref

  avprofx  = time_avg((*i).avprofx_id,/avg)
  avspecty = time_avg((*i).avspecty_id,/avg)

  n_xind = N_ELEMENTS((*i).xind)

;PLOT, (*i).xgrid, (*i).area, color=1

  with_em = (par.n_fields GT 1) AND (par.beta GT 1e-5)
  n_notrap = 2+2*with_em    

  var_names = STRARR(n_notrap)
  var_names2 = STRARR(n_notrap)
  var_names[0:1] = ['!7C!6!Les','!6Q!Les']
  var_names2[0:1] = ['G_es','Q_es']

  IF with_em THEN BEGIN
    var_names[2:3]=['!7C!6!Lem','!6Q!Lem']
    var_names2[2:3]=['G_em','Q_em']
  ENDIF
  header = STRARR(par.n_moms/6*n_notrap)
  header[0:n_notrap-1] = var_names2+',p'
  IF par.n_moms/6 GT 1 THEN BEGIN
    FOR j = 1, par.n_moms/6 - 2 DO $
      header[j*n_notrap:(j+1)*n_notrap-1] = var_names2+',t'+rm0es(par.n_moms/6-2-j)    
    header[j*n_notrap:(j+1)*n_notrap-1] = var_names2+',flr'
  ENDIF

  rho_str = get_var_string(/rhostr)

  avprofx_per_surf = avprofx[*,*,0]

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
  
   sp=(*gui.out.spec_select)[isp]

   timeline = '!6t='+rm0es(gui.out.start_t)+'-'+$
     rm0es(gui.out.end_t)+' '+get_var_string(1,/time,/ounit)+$
     ', ' + STRTRIM(STRING(series.step_count),2) + $
     ', ' + spec[sp].name

   set_output, diag, sp, /ps, multi=[0,1,1], $
     ysize=20, xsize=14.4
     
   IF (*i).n_fluxes EQ n_notrap THEN BEGIN 
     fluxprof_plot, (*i).kyind,avspecty[*,*,isp], '!6k!Iy!N'+rho_str,$
       var_names+'!N', timeline, k_mult=(*i).k_mult, plot_type=2         
     fluxprof_plot, (*i).xgrid,avprofx[*,*,isp], '!6x/a',$
       var_names+'!N', timeline, k_mult=(*i).k_mult, plot_type=1
;     fluxprof_plot, (*i).xgrid,SHIFT(ABS(FFT(avprofx[*,*,isp],/INVERSE,DIMENSION=1)),par.nkx0/2), '!6x/a',$
;       var_names+'!N', timeline, k_mult=(*i).k_mult, plot_type=1
     FOR ivar=0,N_ELEMENTS(avprofx[0,*])-1 DO $
       avprofx_per_surf[*,ivar]=avprofx[*,ivar,isp]*(*i).area
     fluxprof_plot, (*i).xgrid,avprofx_per_surf[*,*], '!6x/a',$
       var_names+'!N', timeline, k_mult=(*i).k_mult, plot_type=1,$
       titleline = 'Flux * '+proj_str
     si_heatflux = avprofx_per_surf[*,1]
     IF with_em THEN si_heatflux += avprofx_per_surf[*,3]
     si_heatflux *= Q_gb*par.Lref^2*1E-6

     fluxprof_plot, (*i).xgrid,si_heatflux, $
       '!6x/a', 'Q!Itot!N', timeline, k_mult=(*i).k_mult, plot_type=1,$
       titleline = 'Heat flux * '+proj_str+' in MW',/no_legend

     q_tot_xavg = TOTAL(avprofx[*,1,isp])/n_xind
     IF (with_em) THEN q_tot_xavg += TOTAL(avprofx[*,3,isp])/n_xind

     print, 'Q'+spec[sp].name+', avgd in x = ['+$
            rm0es((*i).xgrid[0],prec=2)+','+$
            rm0es((*i).xgrid[n_xind-1],prec=2)+']: ', $
            rm0es(q_tot_xavg)+' Q_gb'

     print, 'Q'+spec[sp].name+'*'+proj_str+' avgd in x = ['+$
            rm0es((*i).xgrid[0],prec=2)+','+$
            rm0es((*i).xgrid[n_xind-1],prec=2)+']: ', $
            rm0es(TOTAL(si_heatflux)/n_xind)+' MW'
     

   ENDIF ELSE BEGIN
     fluxtotalx = avprofx[*,0:n_notrap-1,isp]     
     FOR j = 1, par.n_moms/6 - 1 DO BEGIN
        fluxtotalx += avprofx[*,j*n_notrap:(j+1)*n_notrap-1,isp]	
     ENDFOR

     vnames = STRARR(1+par.n_moms/6) ;total, pass, t1, ..., tn, flr
     data_y = DBLARR(par.nky0,1+par.n_moms/6,/NOZERO)
     data_x = DBLARR(FIX(par.nx0/2.0)+1,1+par.n_moms/6,/NOZERO)     

     print, 'flux: transport fractions ('+spec[sp].name+')'	
     FOR v=0, n_notrap-1 DO BEGIN
       vnames[0] = var_names[v]
       vnames[1] = var_names[v]+',pass!N'
       data_y[*,0] = 0.0
       data_x[*,0] = 0.0
       FOR j = 0, par.n_moms/6 - 1 DO BEGIN
         data_y[*,0] += avspecty[*,j*n_notrap+v,isp]       
	 data_x[*,0] += avprofx[*,j*n_notrap+v,isp]       
       ENDFOR
       data_y[*,1] = avspecty[*,v,isp]
       IF TOTAL(data_y[*,0],1) GT 1e-16 THEN print, var_names2[v]+',pass'+ $
    	  '/'+var_names2[v]+',total = '+$       
	  rm0es(TOTAL(data_y[*,1],1)/TOTAL(data_y[*,0],1))     
       data_x[*,1] = avprofx[*,v,isp]	  
       FOR j = 1, par.n_moms/6 - 2 DO BEGIN
        vnames[j+1]=var_names[v]+',trp'+rm0es(par.n_moms/6-2-j)+'!N'
	data_y[*,j+1] = avspecty[*,j*n_notrap+v,isp]
       IF TOTAL(data_y[*,0],1) GT 1e-16 THEN print, var_names2[v]+',trp'+$
         rm0es(par.n_moms/6-2-j)+'/'+var_names2[v]+',total = '+$
	 rm0es(TOTAL(data_y[*,j+1],1)/TOTAL(data_y[*,0],1))
        data_x[*,j+1] = avprofx[*,j*n_notrap+v,isp]
       ENDFOR
       vnames[par.n_moms/6] = var_names[v]+',flr!N'
       data_y[*,par.n_moms/6] = avspecty[*,(par.n_moms/6-1)*$
         n_notrap+v,isp]
       IF TOTAL(data_y[*,0],1) GT 1e-16 THEN print, var_names2[v]+',flr'+ $
    	  '/'+var_names2[v]+',total = '+$              
	  rm0es(TOTAL(data_y[*,par.n_moms/6],1)/TOTAL(data_y[*,0],1)) 
       data_x[*,par.n_moms/6] = avprofx[*,(par.n_moms/6-1)*$
         n_notrap+v,isp]	 
       fluxprof_plot, (*i).kyind,data_y, '!6k!Iy!N'+rho_str,$
         vnames, timeline, k_mult=(*i).k_mult
       fluxprof_plot, (*i).xgrid, data_x,$
        '!6x/a', vnames, timeline, k_mult=(*i).k_mult 
       FOR ivar=0,par.n_moms/6-1 DO data_x = data_x*(*i).area
       fluxprof_plot, (*i).xgrid, data_x,$
        '!6x/a', vnames, timeline, k_mult=(*i).k_mult,titleline='Flux * total flux surface area'
       fluxprof_plot, (*i).xgrid, data_x*Q_gb*par.Lref^2/1E6,$
        '!6x/a', vnames, timeline, k_mult=(*i).k_mult,titleline='Flux * A in MW'
     ENDFOR     
   ENDELSE 

   set_output, diag, sp, header = ['x/a',[header,'Q*'+proj_str+'/MW']],$
     dat=[[(*i).xgrid],[avprofx[*,*,isp]],[si_heatflux]]
   set_output, diag, sp, header = ['ky',header],$
     dat=[[(*i).kyind],[avspecty[*,*,isp]]],/append

     
   set_output, diag, sp, /reset   

 

  ENDFOR ;-- isp
END
