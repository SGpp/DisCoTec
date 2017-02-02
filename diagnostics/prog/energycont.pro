FUNCTION energycont_info

  RETURN, {$
    type      : 'energy',$
    title     : 'Contour plots',$
    help_text : ['Draws a time series of color-coded contour plots '+$
                 'of an energy quantity in a specified x-y plane.'],$
    ext_vars  : [['vars','0','quantities to be plotted; '+$
                  '0: energy, 1: dedt_nc, 2: drive, 3: collisions, 4: dissipation (hyp_z,v,x,y),5: Nonlinearity'],$
    	    	 ['zind','0','parallel (z) index of x-y plane; -1=z-average.  '+$
                  'default: z-average'],$
                 ['kx2x','1','toggle x direction: Fourier to real '+$
                  'space !!!Leave as kx,ky!!!'],$
                 ['ky2y','1','toggle y direction: Fourier to real '+$
                  'space !!!Leave as kx,ky!!!']]}

END

;######################################################################

PRO energycont_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.vars) LT 1 THEN (*e.vars) = [0,1]
  ;IF NOT KEYWORD_SET(*e.zind) THEN (*e.zind) = par.nz0 / 2
  ;IF (*e.zind) GE par.nz0 THEN (*e.zind) = par.nz0 / 2
  IF NOT KEYWORD_SET(*e.zind) THEN (*e.zind) = -1
  IF (*e.zind) GE par.nz0 THEN (*e.zind) = -1
  IF NOT KEYWORD_SET(*e.kx2x) THEN (*e.kx2x) = 0
  IF NOT KEYWORD_SET(*e.ky2y) THEN (*e.ky2y) = 0
  IF par.lx EQ 0 THEN BEGIN
;    printerror, 'Skipping ' + (*diag).name + $
;      ': ky dependent lx not allowed'
;    (*diag).selected = 0
;    RETURN
PRINT, (*diag).name + ' warning: ky dependent lx, results may be unphysical'
lxorig = series.lx
series.lx = 100.0
dxorig = series.dx
series.dx = series.lx / par.nx0
  ENDIF

  ; x direction 0/+1 in real space (local/global BC),
  ; 0/-1 in Fourier (odd/even nkx0)
  add_x_pt = 0
  IF NOT *e.kx2x AND NOT (par.nkx0 MOD 2) THEN add_x_pt = -1
  IF *e.kx2x AND par.x_local AND NOT (par.nkx0 MOD 2) THEN add_x_pt = 1

  ; y direction: +1 in real space (BC), -1 in Fourier (zero mode)
  IF *e.ky2y THEN add_y_pt = 1 ELSE add_y_pt = -1

  i = set_internal_vars(diag,{$
    vars      : *e.vars,$
    add_x_pt  : add_x_pt,$
    add_y_pt  : add_y_pt,$
    zind      : *e.zind,$
    xind      : FINDGEN(par.nkx0+add_x_pt),$
    yind      : FINDGEN(2*par.nky0+add_y_pt),$
    fft_mode  : 3-*e.kx2x-2*(*e.ky2y),$
    time_id   : PTR_NEW(),$
    data_id   : PTR_NEW(),$
    data_t    : PTR_NEW()})

  CASE (*i).fft_mode OF
    0 : BEGIN ; x, y
      fft_format_en, sxsy=*e.vars
      (*i).xind = (*i).xind * series.dx - 0.5 * series.lx
      (*i).yind = (*i).yind * series.dy - 0.5 * series.ly
    END
    1 : BEGIN ; kx, y
      fft_format_en, kxsy=*e.vars
      (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
        (2.0 * !PI / series.lx)
      (*i).yind = (*i).yind * series.dy - 0.5 * series.ly
    END
    2 : BEGIN ; x, ky
      fft_format_en, sxky=*e.vars
      (*i).xind = (*i).xind * series.dx - 0.5 * series.lx
      (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)
    END
    3 : BEGIN ; kx, ky
      fft_format_en, kxky=*e.vars
      (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
        (2.0 * !PI / series.lx)
      (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)
    END
    ELSE:
  ENDCASE

IF par.lx EQ 0 THEN BEGIN
 series.lx = lxorig
 series.dx = dxorig
ENDIF

END

;######################################################################

PRO energycont_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkx0 = par.nkx0
  nky0 = par.nky0

  data = FLTARR(nkx0+(*i).add_x_pt,2*nky0+(*i).add_y_pt,$
    N_ELEMENTS((*i).vars),/NOZERO)

  FOR var = 0, N_ELEMENTS((*i).vars) - 1 DO BEGIN
    CASE (*i).fft_mode OF
      0 : BEGIN ; x, y
        nk_energy=SIZE((*energy[0].sxsy))
        data_temp = FLTARR(nk_energy[1],nk_energy[2],$
          N_ELEMENTS((*i).vars),/NOZERO)
        IF (*i).zind NE -1 THEN BEGIN 
          data_temp[0,0,var] = (*energy[(*i).vars[var]].sxsy)[*,*,(*i).zind]
        ENDIF ELSE BEGIN
         FOR xi=0,nk_energy[1]-1 DO BEGIN
          FOR yi=0,nk_energy[2]-1 DO BEGIN
           data_temp[xi,yi,var]=TOTAL( (*energy[(*i).vars[var]].sxsy)[xi,yi,*]$
		*(*series.geom).jacobian[*])/TOTAL((*series.geom).jacobian)
          ENDFOR
         ENDFOR
        ENDELSE

        ;data[0,0,var] = (*energy[(*i).vars[var]].sxsy)[*,*,(*i).zind]
        data[0,0,var] = data_temp[*,*,var]

        IF (*i).add_x_pt THEN data[nkx0,*,var] = data[0,*,var]
        data[0,2*nky0,var] = data[*,0,var]
      END

      1 : BEGIN ; kx, y

        nk_energy=SIZE((*energy[0].kxsy))
        data_temp = FLTARR(nk_energy[1],nk_energy[2],$
          N_ELEMENTS((*i).vars),/NOZERO)
        IF (*i).zind NE -1 THEN BEGIN 
          data_temp[0,0,var] = (*energy[(*i).vars[var]].kxsy)[*,*,(*i).zind]
        ENDIF ELSE BEGIN
         FOR xi=0,nk_energy[1]-1 DO BEGIN
          FOR yi=0,nk_energy[2]-1 DO BEGIN
           data_temp[xi,yi,var]=TOTAL( (*energy[(*i).vars[var]].kxsy)[xi,yi,*]$
		*(*series.geom).jacobian[*])/TOTAL((*series.geom).jacobian)
          ENDFOR
         ENDFOR
        ENDELSE

        ; zero and positive kx
        data[nkx0-nkx0/2-1,0,var] = $
           data_temp[0:nkx0-nkx0/2-1,*,var]
          ;(*energy[(*i).vars[var]].kxsy)[0:nkx0-nkx0/2-1,*,(*i).zind]
        ; negative kx
        data[0,0,var] = REVERSE(data[nkx0-nkx0/2:*,*],1)

        data[0,2*nky0,var] = data[*,0,var]
      END

      2 : BEGIN ; x, ky

        nk_energy=SIZE((*energy[0].sxky))
        data_temp = FLTARR(nk_energy[1],nk_energy[2],$
          N_ELEMENTS((*i).vars),/NOZERO)
        IF (*i).zind NE -1 THEN BEGIN 
          data_temp[0,0,var] = (*energy[(*i).vars[var]].sxky)[*,*,(*i).zind]
        ENDIF ELSE BEGIN
         FOR xi=0,nk_energy[1]-1 DO BEGIN
          FOR yi=0,nk_energy[2]-1 DO BEGIN
           data_temp[xi,yi,var]=TOTAL( (*energy[(*i).vars[var]].sxky)[xi,yi,*]$
		*(*series.geom).jacobian[*])/TOTAL((*series.geom).jacobian)
          ENDFOR
         ENDFOR
        ENDELSE

        ; zero and positive ky
        data[0,nky0-1,var] = (data_temp[*,*,var])
          ;((*energy[(*i).vars[var]].sxky)[*,*,(*i).zind])
        ; negative ky; the loop is faster than REVERSE
        FOR y = 0, nky0 - 2 DO $
          data[0,y,var] = data[*,2*nky0-2-y,var]

        IF (*i).add_x_pt THEN data[nkx0,*,var] = data[0,*,var]
      END

      3 : BEGIN ; kx, ky

        nk_energy=SIZE((*energy[0].kxky))
        data_temp = FLTARR(nk_energy[1],nk_energy[2],$
          N_ELEMENTS((*i).vars),/NOZERO)
        IF (*i).zind NE -1 THEN BEGIN 
          data_temp[0,0,var] = (*energy[(*i).vars[var]].kxky)[*,*,(*i).zind]
        ENDIF ELSE BEGIN
         FOR xi=0,nk_energy[1]-1 DO BEGIN
          FOR yi=0,nk_energy[2]-1 DO BEGIN
           data_temp[xi,yi,var]=TOTAL( (*energy[(*i).vars[var]].kxky)[xi,yi,*]$
		*(*series.geom).jacobian[*])/TOTAL((*series.geom).jacobian)
          ENDFOR
         ENDFOR
        ENDELSE
        ; zero and positive kx, zero and positive ky
        data[nkx0-nkx0/2-1,nky0-1,var] = (data_temp[0:nkx0-nkx0/2-1,*,var])
          ;((*energy[(*i).vars[var]].kxky)[0:nkx0-nkx0/2-1,*,(*i).zind])
        ; negative kx, positive ky
        ;data[0,nky0-1,var] = ((*energy[(*i).vars[var]].kxky)$
        ;  [nkx0/2+1:nkx0-1,*,(*i).zind])
        data[0,nky0-1,var] = data_temp[nkx0/2+1:nkx0-1,*,var]
        ; negative ky
        data[0,0,var] = REVERSE(REVERSE(data[*,nky0-1:*,var],2),1)
      END

      ELSE:
    ENDCASE
  ENDFOR

  (*i).time_id = store_step((*i).time_id,energy_time)
  (*i).data_id = store_step((*i).data_id,data)
  (*i).data_t = time_avg((*i).data_t,data,energy_time)

END

;######################################################################

FUNCTION get_pos_encont, aspect_ratio, ivar, n_vars, x_as_y

  pos = FLTARR(10)  

  IF aspect_ratio LT 1.0 THEN BEGIN ; vertical
    lmargin = 6*!D.Y_CH_SIZE & rmargin = 8*!D.X_CH_SIZE
    bmargin = 4*!D.Y_CH_SIZE & tmargin = 1*!D.Y_CH_SIZE
    between = 0.5*!D.Y_CH_SIZE
    height=(!D.Y_SIZE-tmargin-bmargin-(n_vars-1)*between)/n_vars
    width = !D.X_SIZE-rmargin-lmargin
    IF x_as_y THEN $
      IF (aspect_ratio LE FLOAT(height)/width) THEN height = width*aspect_ratio $
      ELSE width = height/aspect_ratio

    res = convert_coord([lmargin,lmargin+width],$
                        [!D.Y_SIZE-tmargin-ivar*between-(ivar+1.0)*height,$
                         !D.Y_SIZE-tmargin-ivar*(height+between)],$
                        /device,/to_normal)
    pos[0] = res[0,0] & pos[1] = res[1,0] ; CONTOUR left, bottom
    pos[2] = res[0,1] & pos[3] = res[1,1] ; CONTOUR right, top  
    res = convert_coord([2*!D.X_CH_SIZE],$
                        [!D.Y_SIZE-tmargin-ivar*(height+between)-!D.Y_CH_SIZE],$
                        /device,/to_normal)
    pos[4] = res[0,0] & pos[5] = res[1,0] ; Var_names
    res = convert_coord([lmargin+width+1*!D.X_CH_SIZE, lmargin+width+3*!D.X_CH_SIZE],$
                        [!D.Y_SIZE-tmargin-ivar*between-(ivar+1.0)*height,$
                         !D.Y_SIZE-tmargin-ivar*(height+between)],$
                        /device,/to_normal)
    pos[6] = res[0,0] & pos[7] = res[1,0] ; COLORBAR left, bottom
    pos[8] = res[0,1] & pos[9] = res[1,1] ; COLORBAR right, top         
    
  ENDIF ELSE BEGIN
    lmargin = 6*!D.Y_CH_SIZE & rmargin = 3*!D.X_CH_SIZE
    bmargin = 4*!D.Y_CH_SIZE & tmargin = 5*!D.Y_CH_SIZE  
    between = 8*!D.X_CH_SIZE
    height=!D.Y_SIZE-tmargin-bmargin
    width =(!D.X_SIZE-rmargin-lmargin-(n_vars-1)*between)/n_vars
    
    IF x_as_y THEN $
      IF (aspect_ratio LE FLOAT(height)/width) THEN height = FIX(width*aspect_ratio) $
      ELSE width = FIX(height/aspect_ratio)/n_vars
    
    res = convert_coord([lmargin+ivar*(width+between),$
                         lmargin+ivar*(width+between)+width],$
                        [!D.Y_SIZE-tmargin-height, !D.Y_SIZE-tmargin],$
                        /device,/to_normal)
    
    pos[0] = res[0,0] & pos[1] = res[1,0] ; CONTOUR left, bottom
    pos[2] = res[0,1] & pos[3] = res[1,1] ; CONTOUR right, top  
    res = convert_coord([lmargin-3*!D.X_CH_SIZE+ivar*(width+between)],$
                        [!D.Y_SIZE-2.5*!D.Y_CH_SIZE],$
                        /device,/to_normal)
    pos[4] = res[0,0] & pos[5] = res[1,0] ; Var_names           
    res = convert_coord([lmargin+ivar*(width+between),$
                         lmargin+ivar*(width+between)+width],$
                        [!D.Y_SIZE-3*!D.Y_CH_SIZE,$
                         !D.Y_SIZE-!D.Y_CH_SIZE],$
                        /device,/to_normal) 
    pos[6] = res[0,0] & pos[7] = res[1,0] ; COLORBAR left, bottom
    pos[8] = res[0,1] & pos[9] = res[1,1] ; COLORBAR right, top 
  ENDELSE
  
  RETURN, pos
    
END

;######################################################################

PRO energycont_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF par.arakawa_zv AND  NOT par.arakawa_cons_bc THEN BEGIN
    en_var_string = '!6' + ['E!Dtot!N','dE/dt!9!!!6!D NC!N (Warning:z,v,poisson)',$
     'Q (= dE/dt!9!!!6!Ddrive!N)','C (= dE/dt!9!!!6!Dcoll!N)','D (hyp_z,v,x,y)(No diss b.c)',$
     'dE/dt!D NL !N','1/E dE/dt !D NL !N']
  ENDIF ELSE BEGIN
    en_var_string = '!6' + ['E!Dtot!N','dE/dt!9!!!6!D NC!N (Warning:z,v,poisson)',$
     'Q (= dE/dt!9!!!6!Ddrive!N)','C (= dE/dt!9!!!6!Dcoll!N)','D (hyp_z,v,x,y)',$
     'dE/dt!D NL !N','1/E dE/dt !D NL !N']
  ENDELSE
  
  n_vars = N_ELEMENTS((*i).vars)
  count = 0

  x_as_y = 0
  aspect_ratio = (*i).yind[N_ELEMENTS((*i).yind)-1] / $
    (*i).xind[N_ELEMENTS((*i).xind)-1]
    
  ; identical scaling of x and y axes for sxsy, kxky
  IF ((*i).fft_mode EQ 0) OR ((*i).fft_mode EQ 3) THEN x_as_y = 1
      
  IF (aspect_ratio*n_vars LT 1.0) THEN BEGIN    
    xysize = [25,17]
  ENDIF ELSE BEGIN  
    xysize = [17,25]
  ENDELSE

  x_axis = 'x/'
  y_axis = 'y/'
  IF (*i).fft_mode MOD 2 EQ 1 THEN x_axis = 'k!Dx!N '
  IF (*i).fft_mode GT 1 THEN y_axis = 'k!Dy!N '
  
  xtitle = '!6' + x_axis + get_var_string(/rhostr)
  ytitle = '!6' + y_axis + get_var_string(/rhostr)

  csize = 1.41 - (n_vars < 4) * 0.125
  cthick = 3.0

  c_levels = 20
  cont_minmax = FLTARR(2,n_vars)
  
  time = store_step((*i).time_id,/get)
  data = store_step((*i).data_id,/get)
  data_t = time_avg((*i).data_t,/avg,fwd=0, tarr=tarr)

  set_output, diag, coltable=33, $
    xsize=xysize[0], ysize=xysize[1], charsize=csize,$
    ps=(series.step_count GT 1), eps=(series.step_count EQ 1)

  first = 1
   
  FOR n = 0, series.step_count - 1 DO BEGIN
    tstring='t='+STRING(*time[n],format='(f6.2)')
    FOR var = 0, n_vars - 1 DO BEGIN    ; loop over all variables to plot

      lev_col = contour_levels([MIN((*data[n])[*,*,var]),$
        MAX((*data[n])[*,*,var])],c_levels)
      
      pos = get_pos_encont(aspect_ratio,var,n_vars,x_as_y)
    
      CONTOUR, (*data[n])[*,*,var], (*i).xind, (*i).yind, $
        LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
        /XSTYLE, /YSTYLE, $
        XTITLE=xtitle, YTITLE=ytitle, POSITION=pos[0:3],$
        COLOR=((*i).fft_mode EQ 3)*55, ISOTROPIC=x_as_y, NOERASE=(var NE 0), $
        CHARSIZE=csize, CHARTHICK=cthick, /NORMAL, YTICK_GET=YTICKV, $
	XTICK_GET=XTICKV

      IF ((*i).fft_mode EQ 3) THEN BEGIN
        AXIS, 0,0, XAX=0,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0
        AXIS, 0,0, XAX=1,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0	  
        AXIS, 0,0, YAX=0,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0
        AXIS, 0,0, YAX=1,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0	  
      ENDIF

      XYOUTS, pos[4], pos[5],/NORMAL, COLOR=1,$
        en_var_string[(*i).vars[var]], CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
      XYOUTS, 0.8,0.01,/NORMAL, COLOR=1,$
        tstring, CHARSIZE=1.5*csize,$
        CHARTHICK=cthick

      plot_colorbar, lev_col, position=pos[6:9], charsize=0.75*csize
	
      ; data output	
      set_output, diag, dat=[[REFORM((*data[n])[*,*,var])]],$
        commentlines='time: '+rm0es(*time[n])+ ', variable: '+$
        rm_idl_fmt(en_var_string[(*i).vars[var]])+'('+rm_idl_fmt(xtitle)+','+$
        rm_idl_fmt(ytitle)+')',append = (first EQ 0)
      first =0 	
    ENDFOR

IF 0 THEN BEGIN    
    XYOUTS, 0, 0, 't='+rm0es(*time[n])+' '+$
      get_var_string(0,/time,/ounit), /NORMAL, COLOR=1,CHARSIZE=csize
ENDIF
  ENDFOR

    ;Now output time averaged contour plots
    ;Now output time averaged contour plots
    ;Now output time averaged contour plots
  IF par.arakawa_zv AND  NOT par.arakawa_cons_bc THEN BEGIN
   en_var_string = '!6' + ['<E!Dtot!N>!Dt!N','<dE/dt!9!!!6!Dnonconserve!N>!Dt!N',$
    '<Q (= dE/dt!9!!!6!Ddrive!N)>!Dt!N','<C (= dE/dt!9!!!6!Dcoll!N)>!Dt!N',$
    '<D=hyp_z,v,x,y>!Dt!N(Warning:no diss b.c)',$
     '<dE/dt!D NL !N>!Dt!N','<1/E dE/dt !D NL !N>!Dt!N']
  ENDIF ELSE BEGIN
   en_var_string = '!6' + ['<E!Dtot!N>!Dt!N','<dE/dt!9!!!6!Dnonconserve!N>!Dt!N',$
    '<Q (= dE/dt!9!!!6!Ddrive!N)>!Dt!N','<C (= dE/dt!9!!!6!Dcoll!N)>!Dt!N',$
    '<D=hyp_z,v,x,y>!Dt!N',$
     '<dE/dt!D NL !N>!Dt!N','<1/E dE/dt !D NL !N>!Dt!N']
  ENDELSE


    FOR var = 0, n_vars - 1 DO BEGIN    ; loop over all variables to plot

      lev_col = contour_levels([MIN(data_t[*,*,var]),$
        MAX(data_t[*,*,var])],c_levels)
      
      pos = get_pos_encont(aspect_ratio,var,n_vars,x_as_y)
      CONTOUR, data_t[*,*,var], (*i).xind, (*i).yind, $
        LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
        /XSTYLE, /YSTYLE, $
        XTITLE=xtitle, YTITLE=ytitle, POSITION=pos[0:3],$
        COLOR=((*i).fft_mode EQ 3)*55, ISOTROPIC=x_as_y, NOERASE=(var NE 0), $
        CHARSIZE=csize, CHARTHICK=cthick, /NORMAL, YTICK_GET=YTICKV, $
	XTICK_GET=XTICKV

      IF ((*i).fft_mode EQ 3) THEN BEGIN
        AXIS, 0,0, XAX=0,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0
        AXIS, 0,0, XAX=1,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0	  
        AXIS, 0,0, YAX=0,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0
        AXIS, 0,0, YAX=1,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0	  
      ENDIF

      XYOUTS, pos[4], pos[5],/NORMAL, COLOR=1,$
        en_var_string[(*i).vars[var]], CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
       ntime=N_ELEMENTS(time)
       tmin=*time[0]
       tmax=*time[ntime-1]
       trange='t=['+STRING(tmin,format='(f6.2)')+','+STRING(tmax,format='(f6.2)')+']'
      XYOUTS, 0.7,0.01,/NORMAL, COLOR=1,$
        trange, CHARSIZE=1.5*csize,$
        CHARTHICK=cthick

      plot_colorbar, lev_col, position=pos[6:9], charsize=0.75*csize
	
      ; data output	
      ;set_output, diag, dat=[[REFORM((*data[n])[*,*,var])]],$
      ;  commentlines='time: '+rm0es(*time[n])+ ', variable: '+$
      ;  rm_idl_fmt(en_var_string[(*i).vars[var]])+'('+rm_idl_fmt(xtitle)+','+$
      ;  rm_idl_fmt(ytitle)+')',append = (first EQ 0)
      ;first =0 	
    ENDFOR
; IF n_vars GT 2 THEN BEGIN
; Now plot time averaged dissipation
; diss=data_t[*,*,1]-data_t[*,*,2]
; lev_col = contour_levels([MIN(diss[*,*]),$
;        MAX(diss[*,*])],c_levels)
;     
;     pos = get_pos_encont(aspect_ratio,0,n_vars,x_as_y)
;     CONTOUR, diss[*,*], (*i).xind, (*i).yind, $
;       LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
;       /XSTYLE, /YSTYLE, $
;       XTITLE=xtitle, YTITLE=ytitle, POSITION=pos[0:3],$
;       COLOR=((*i).fft_mode EQ 3)*55, ISOTROPIC=x_as_y,  $
;       CHARSIZE=csize, CHARTHICK=cthick, /NORMAL, YTICK_GET=YTICKV, $
;XTICK_GET=XTICKV
;
;     IF ((*i).fft_mode EQ 3) THEN BEGIN
;       AXIS, 0,0, XAX=0,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
;  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0
;       AXIS, 0,0, XAX=1,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
;  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0	  
;       AXIS, 0,0, YAX=0,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
;  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0
;       AXIS, 0,0, YAX=1,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
;  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0	  
;     ENDIF
;
;     XYOUTS, pos[4], pos[5],/NORMAL, COLOR=1,$
;       en_var_string[4], CHARSIZE=1.5*csize,$
;       CHARTHICK=cthick
;
;     plot_colorbar, lev_col, position=pos[6:9], charsize=0.75*csize
;ENDIF

  set_output, diag, /reset, eps=(series.step_count EQ 1)

  PTR_FREE, time, data

END
