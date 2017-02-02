FUNCTION cont_info

  RETURN, {$
    type      : 'mom',$
    title     : 'Contour plots',$
    help_text : ['Draws a time series of color-coded contour plots '+$
                 'of a fluctuating quantity in a chosen x-y plane.'],$
    ext_vars  : [['vars','0','quantities to be plotted (see '+$
    	    	  'variable list)'],$
    	    	 ['zind','0','parallel (z) index of x-y plane; '+$
                  'default: outboard midplane (nz0/2); set to '+$
                  '-1 for average'],$
                 ['kx2x','1','toggle x direction: Fourier to real '+$
                  'space'],$
                 ['ky2y','1','toggle y direction: Fourier to real '+$
                  'space'],$
                 ['ps2mpg','1','toggle mpg output instead of ps']]}

END

;######################################################################

PRO cont_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  vars = N_ELEMENTS(*e.vars) GE 1 ? *e.vars : [0,par.n_fields]
  zind = N_ELEMENTS(*e.zind) EQ 1 ? *e.zind : par.nz0 / 2
  IF zind EQ -1 THEN BEGIN
    IF NOT par.x_local OR NOT par.y_local THEN BEGIN
      printerror, 'Skipping ' + (*diag).name + $
      ': z average implemented only for local runs'
      (*diag).selected = 0
      RETURN
    ENDIF

    nz = par.nz0
    zind = INDGEN(nz)
  ENDIF ELSE nz = 1
  IF NOT KEYWORD_SET(*e.kx2x) THEN *e.kx2x = 0
  IF NOT KEYWORD_SET(*e.ky2y) THEN *e.ky2y = 0
  IF NOT KEYWORD_SET(*e.ps2mpg) THEN *e.ps2mpg = 0
  IF par.lx EQ 0 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + $
      ': ky dependent lx not allowed'
    (*diag).selected = 0
    RETURN
  ENDIF

  ; x direction 0/+1 in real space (local/global BC),
  ; 0/-1 in Fourier (odd/even nkx0)
  add_x_pt = 0
  IF NOT *e.kx2x AND NOT (par.nkx0 MOD 2) THEN add_x_pt = -1
  IF *e.kx2x AND par.x_local AND NOT (par.nkx0 MOD 2) THEN add_x_pt = 1

  ; y direction: +1 in real space (BC), -1 in Fourier (zero mode)
  IF *e.ky2y THEN add_y_pt = 1 ELSE add_y_pt = -1

  i = set_internal_vars(diag,{$
    vars      : vars,$
    add_x_pt  : add_x_pt,$
    add_y_pt  : add_y_pt,$
    nz        : nz,$
    zind      : zind,$
    xind      : FINDGEN(par.nkx0+add_x_pt),$
    yind      : FINDGEN(2*par.nky0+add_y_pt),$
    fft_mode  : 3-*e.kx2x-2*(*e.ky2y),$
    data_id   : PTR_NEW(),$
    mpg_out   : *e.ps2mpg})

  CASE (*i).fft_mode OF
    0 : BEGIN ; x, y
      fft_format, sxsy=vars
      (*i).xind = (*i).xind * series.dx - 0.5 * series.lx
      (*i).yind = (*i).yind * series.dy - 0.5 * series.ly
    END
    1 : BEGIN ; kx, y
      fft_format, kxsy=vars
      (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
        (2.0 * !PI / series.lx)
      (*i).yind = (*i).yind * series.dy - 0.5 * series.ly
    END
    2 : BEGIN ; x, ky
      fft_format, sxky=vars
      (*i).xind = (*i).xind * series.dx - 0.5 * series.lx
      (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)
    END
    3 : BEGIN ; kx, ky
      fft_format, kxky=vars
      (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
        (2.0 * !PI / series.lx)
      (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)
    END
    ELSE:
  ENDCASE

END

;######################################################################

PRO cont_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkx0 = par.nkx0
  nky0 = par.nky0

  IF (*i).nz GT 1 THEN jac_norm = REBIN(REFORM($
    (*series.geom).jac_norm/par.nz0,[1,1,par.nz0]),$
    [par.nx0,par.ny0,par.nz0])

  data = FLTARR(nkx0+(*i).add_x_pt,2*nky0+(*i).add_y_pt,$
    N_ELEMENTS((*i).vars),gui.out.n_spec_sel,/NOZERO)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    FOR var = 0, N_ELEMENTS((*i).vars) - 1 DO BEGIN
      CASE (*i).fft_mode OF
        0 : BEGIN ; x, y
          data[0,0,var,isp] = (*i).nz GT 1 ? TOTAL($
            (*mom[isp,(*i).vars[var]].sxsy)*jac_norm,3) : $
            (*mom[isp,(*i).vars[var]].sxsy)[*,*,(*i).zind]

          IF (*i).add_x_pt THEN data[nkx0,*,var,isp] = data[0,*,var,isp]
          data[0,2*nky0,var,isp] = data[*,0,var,isp]
        END

        1 : BEGIN ; kx, y
          ; zero and positive kx
          data[nkx0-nkx0/2-1,0,var,isp] = (*i).nz GT 1 ? ABS(TOTAL($
            (*mom[isp,(*i).vars[var]].kxsy)[0:nkx0-nkx0/2-1,*,*]*jac_norm,3)) : $
            ABS((*mom[isp,(*i).vars[var]].kxsy)[0:nkx0-nkx0/2-1,*,(*i).zind])
          ; negative kx
          data[0,0,var,isp] = REVERSE(data[nkx0-nkx0/2:*,*],1)

          data[0,2*nky0,var,isp] = data[*,0,var,isp]
        END

        2 : BEGIN ; x, ky
          ; zero and positive ky
          data[0,nky0-1,var,isp] = (*i).nz GT 1 ? ABS(TOTAL($
            (*mom[isp,(*i).vars[var]].sxky)*jac_norm,3)) : $
            ABS((*mom[isp,(*i).vars[var]].sxky)[*,*,(*i).zind])
          ; negative ky; the loop is faster than REVERSE
          FOR y = 0, nky0 - 2 DO $
            data[0,y,var,isp] = data[*,2*nky0-2-y,var,isp]

          IF (*i).add_x_pt THEN data[nkx0,*,var,isp] = data[0,*,var,isp]
        END

        3 : BEGIN ; kx, ky
          ; zero and positive kx, zero and positive ky
          data[nkx0-nkx0/2-1,nky0-1,var,isp] = (*i).nz GT 1 ? ABS(TOTAL($
            (*mom[isp,(*i).vars[var]].kxky)[0:nkx0-nkx0/2-1,*,*]*jac_norm,3)) : $
            ABS((*mom[isp,(*i).vars[var]].kxky)[0:nkx0-nkx0/2-1,*,(*i).zind])
          ; negative kx, positive ky
          data[0,nky0-1,var,isp] = (*i).nz GT 1 ? ABS(TOTAL($
            (*mom[isp,(*i).vars[var]].kxky)[nkx0/2+1:nkx0-1,*,*]*jac_norm,3)) : $
            ABS((*mom[isp,(*i).vars[var]].kxky)[nkx0/2+1:nkx0-1,*,(*i).zind])
          ; negative ky
          data[0,0,var,isp] = REVERSE(REVERSE(data[*,nky0-1:*,var,isp],2),1)
        END

        ELSE :
      ENDCASE
    ENDFOR
  ENDFOR ; --- species loop

  (*i).data_id = time_avg((*i).data_id,data,mom_time,$
    fwd=gui.out.res_steps)

END

;######################################################################

FUNCTION get_pos, aspect_ratio, ivar, n_vars, x_as_y, mpgoff

  pos = FLTARR(10)  

  IF aspect_ratio LT 1.0 THEN BEGIN ; vertical
    lmargin = 6*!D.Y_CH_SIZE*mpgoff & rmargin = 8*!D.X_CH_SIZE*mpgoff
    bmargin = 4*!D.Y_CH_SIZE*mpgoff & tmargin = 1*!D.Y_CH_SIZE*mpgoff
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
                        [!D.Y_SIZE-tmargin-ivar*(height+between)-(2-mpgoff)*!D.Y_CH_SIZE],$
                        /device,/to_normal)
    pos[4] = res[0,0] & pos[5] = res[1,0] ; Var_names
    res = convert_coord([lmargin+width+1*!D.X_CH_SIZE, lmargin+width+3*!D.X_CH_SIZE],$
                        [!D.Y_SIZE-tmargin-ivar*between-(ivar+1.0)*height,$
                         !D.Y_SIZE-tmargin-ivar*(height+between)],$
                        /device,/to_normal)
    pos[6] = res[0,0] & pos[7] = res[1,0] ; COLORBAR left, bottom
    pos[8] = res[0,1] & pos[9] = res[1,1] ; COLORBAR right, top         
    
  ENDIF ELSE BEGIN
    lmargin = 6*!D.Y_CH_SIZE*mpgoff & rmargin = 3*!D.X_CH_SIZE*mpgoff
    bmargin = 4*!D.Y_CH_SIZE*mpgoff & tmargin = 5*!D.Y_CH_SIZE*mpgoff  
    between = 8*!D.X_CH_SIZE*mpgoff+(1-mpgoff)*!D.X_CH_SIZE
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
    res = convert_coord([lmargin-3*!D.X_CH_SIZE*mpgoff+ivar*(width+between)],$
                        [!D.Y_SIZE-2.5*(2-mpgoff)*!D.Y_CH_SIZE],$
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

PRO cont_output, diag

  COMMON global_vars

  i = (*diag).internal_vars  
  
  n_vars = N_ELEMENTS((*i).vars)
  count = 0

  aspect_ratio = (*i).yind[N_ELEMENTS((*i).yind)-1] / $
    (*i).xind[N_ELEMENTS((*i).xind)-1]
    
  ; identical scaling of x and y axes for sxsy, kxky
  x_as_y = ((*i).fft_mode EQ 0) OR ((*i).fft_mode EQ 3)
      
  npixel = 800>(par.nx0>par.ny0) ;minimum 800 pixels for movies in larger dim.
  IF aspect_ratio * n_vars LT 1 THEN BEGIN    
    xysize = [25,17]
    mpg_resolution = [npixel,ROUND(npixel*aspect_ratio)]
  ENDIF ELSE BEGIN  
    xysize = [17,25]
    mpg_resolution = [ROUND(npixel/aspect_ratio),npixel]
  ENDELSE

  x_axis = 'x/'
  y_axis = 'y/'
  IF (*i).fft_mode MOD 2 EQ 1 THEN x_axis = 'k!Dx!N '
  IF (*i).fft_mode GT 1 THEN y_axis = 'k!Dy!N '
  
  xtitle = '!6' + x_axis + get_var_string(/rhostr)
  ytitle = '!6' + y_axis + get_var_string(/rhostr)

  csize = 1.41 - (n_vars < 4) * 0.125
  cthick = 3.0

  IF (*i).mpg_out THEN c_levels = 50 $  ;use more levels for movies
  ELSE c_levels = 20

  ;switch on/off global min/max (in time) for color levels here:
  use_global_minmax = ((*i).mpg_out AND par.nonlinear)
  cont_minmax = FLTARR(2,n_vars)

  tintrpfac = ((*i).mpg_out AND par.nonlinear) ;interpol. to equidist. time steps for movies
  n_res_steps = gui.out.res_steps * ((tintrpfac>1)*$
                series.step_count - 1) + 1
  
  data = time_avg((*i).data_id,/avg,fwd=gui.out.res_steps,tarr=time,intrp=tintrpfac)
  data = REFORM(data,[par.nkx0+(*i).add_x_pt,2*par.nky0+(*i).add_y_pt,$
    N_ELEMENTS((*i).vars),gui.out.n_spec_sel,n_res_steps],/OVERWRITE)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    IF (*i).mpg_out THEN BEGIN
      set_output, diag, sp, coltable=33, $
        xsize=xysize[0], ysize=xysize[1], charsize=csize, /reset

      IF use_global_minmax THEN BEGIN
         FOR n = 0, n_res_steps - 1 DO BEGIN
            FOR var = 0, n_vars - 1 DO BEGIN
               min = MIN(data[*,*,var,isp,n],MAX=max)
               IF min LT cont_minmax[0,var] THEN $
                  cont_minmax[0,var] = min
               IF max GT cont_minmax[1,var] THEN $
                  cont_minmax[1,var] = max
            ENDFOR
         ENDFOR
      ENDIF

      mpgfile = 1
    ENDIF 
   
    set_output, diag, sp,coltable=33, $
      xsize=xysize[0], ysize=xysize[1], charsize=csize,$
      ps=(n_res_steps GT 1), eps=(n_res_steps EQ 1), $
      resolution=mpg_resolution, mpgfile=mpgfile

    first = 1
    framenr = 1
   
    FOR n = 0, n_res_steps - 1 DO BEGIN
      FOR var = 0, n_vars - 1 DO BEGIN ; loop over all variables to plot
        IF (use_global_minmax) THEN lev_col = $
           contour_levels(cont_minmax[*,var],c_levels) $
        ELSE lev_col = contour_levels([MIN(data[*,*,var,isp,n]),$
           MAX(data[*,*,var,isp,n])],c_levels)

        pos = get_pos(aspect_ratio,var,n_vars,x_as_y,1-(*i).mpg_out)
    
        CONTOUR, data[*,*,var,isp,n], (*i).xind, (*i).yind, $
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL, $
          XSTYLE=1+4*(*i).mpg_out, YSTYLE=1+4*(*i).mpg_out, $
          XTITLE=xtitle, YTITLE=ytitle, POSITION=pos[0:3], $
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
        XYOUTS, pos[4], pos[5], /NORMAL, COLOR=1, $
          get_var_string((*i).vars[var],/fancy), CHARSIZE=1.5*csize, $
          CHARTHICK=cthick

        IF NOT (*i).mpg_out THEN $
          plot_colorbar, lev_col, position=pos[6:9], charsize=0.75*csize
	
        ; data output	
        set_output, diag, sp, dat=[[REFORM(data[*,*,var,isp,n])]], $
          commentlines=(gui.out.res_steps ? ', time = '+rm0es(time[n]):'avg')+$
          ', variable: '+get_var_string((*i).vars[var])+$
          '('+rm_idl_fmt(xtitle)+','+$
          rm_idl_fmt(ytitle)+')', append=(first EQ 0)
        first = 0 	
      ENDFOR

      IF NOT (*i).mpg_out THEN plot_info_str, diag, $
        time=(gui.out.res_steps ? time[n] : undef)

      set_output, diag, sp, write_frame=framenr, $
        mpgfile=mpgfile, resolution=mpg_resolution
      framenr += 1
    ENDFOR

    set_output, diag, sp, /reset, eps=(n_res_steps EQ 1), $
      write_frame=framenr, mpgfile=mpgfile, resolution=mpg_resolution
  ENDFOR ; --- species loop

END
