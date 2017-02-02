FUNCTION nltspect_info

  RETURN, {$
    type      : 'nlt',$
    title     : 'NLT spectra',$
    help_text : ['Plots nonlinear energy transfer functions '],$
    ext_vars  : [['kx_data','0','  1000: sum over kx.  '+$
                  'Otherwise: slice at kx index (negative indices allowed).'+$
                  '  (Default: sum).'],$
                 ['ky_data','0','  1000: sum over ky.  '+$
                  'Otherwise: slice at ky index (negative indices allowed).'+$
                  '  (Default: sum).'],$
    	    	 ['kx_range','0', '  Specifies plot range in kx e.g., [-0.5,0.5] (Default: max,min)' ],$
    	    	 ['ky_range','0', '  Specifies plot range in ky e.g., [-1.5,1.5] (Default: max,min)' ]$
                                                                                              ]}





END

;######################################################################

PRO nltspect_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF NOT KEYWORD_SET(*e.kx_range) THEN (*e.kx_range)=[999,999]
  IF NOT KEYWORD_SET(*e.ky_range) THEN (*e.ky_range)=[999,999]
  ;IF *e.kx_data EQ !NULL THEN (*e.kx_data)=1000
  ;IF *e.ky_data EQ !NULL THEN (*e.ky_data)=1000
  IF N_ELEMENTS(*e.kx_data) EQ 0 THEN (*e.kx_data)=1000
  IF N_ELEMENTS(*e.ky_data) EQ 0 THEN (*e.ky_data)=1000

  IF par.lx EQ 0 THEN BEGIN
    STOP, 'lx cannot be 0'
  ENDIF

  i = set_internal_vars(diag,{$
    kx_data   : *e.kx_data,$
    ky_data   : *e.ky_data,$
    kx_range   : *e.kx_range,$
    ky_range   : *e.ky_range,$
    yind      : FINDGEN(2*par.nky0-1),$
    xind      : FINDGEN(par.nkx0),$
    time_id   : PTR_NEW(),$
    data_kx   : PTR_NEW(),$
    data_ky   : PTR_NEW(),$
    data_id   : PTR_NEW()})

    (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
      (2.0 * !PI / series.lx)
    (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)

END

;######################################################################

PRO nltspect_loop, diag

  COMMON global_vars

  ;jac_fac = (*series.geom).jacobian

  i = (*diag).internal_vars

  nkx0 = par.nkx0
  nky0 = par.nky0

  num_nlt_out = par.num_nlt_modes
  IF par.nlp_gdt THEN num_nlt_out = par.num_nlt_pod_modes+1

  data_kx = FLTARR(nkx0,num_nlt_out,/NOZERO)
  data_ky = FLTARR(2*nky0-1,num_nlt_out,/NOZERO)
  ;data_full = FLTARR(nkx0,2*nky0-1,par.num_nlt_modes,/NOZERO)

  FOR var = 0, num_nlt_out - 1 DO BEGIN
   data_kx[*,var]=TOTAL((*nlt2d)[*,*,var],2)
   data_ky[*,var]=TOTAL((*nlt2d)[*,*,var],1)
  ENDFOR

  data=*nlt2d

  (*i).time_id = store_step((*i).time_id,nlt_time)
  (*i).data_kx = time_avg((*i).data_kx,data_kx,nlt_time,fwd=0)
  (*i).data_ky = time_avg((*i).data_ky,data_ky,nlt_time,fwd=0)
  (*i).data_id = time_avg((*i).data_id,data,nlt_time,fwd=0)

END

;######################################################################

FUNCTION get_pos_enspect, aspect_ratio, ivar, n_vars, x_as_y

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

PRO nltspect_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  out_suffix=''
  IF par.nlp_gdt THEN out_suffix='_POD'

   data_full = time_avg((*i).data_id,/avg,fwd=0, tarr=tarr)
   xind=(*i).xind
   yind=(*i).yind
   num_nlt_out = par.num_nlt_modes
   IF par.nlp_gdt THEN num_nlt_out = par.num_nlt_pod_modes + 1
   en_var_string_kx=STRARR(num_nlt_out)
   en_var_string_ky=STRARR(num_nlt_out)
   comments_kx=STRARR(num_nlt_out)
   comments_ky=STRARR(num_nlt_out)
   kxmin=2.0*!PI/par.lx
   kxstring=STRARR(num_nlt_out)
   kxind0=FLTARR(num_nlt_out)

   IF (*i).kx_data EQ 1000 THEN BEGIN
    data_ky = time_avg((*i).data_ky,/avg,fwd=0, tarr=tarr)

  FOR var=0,num_nlt_out-1 DO BEGIN
     IF par.nlp_gdt THEN BEGIN
      kxind0[var]=kxmin*par.nlp_kxind
      IF par.nlp_kxind GT par.nx0/2 THEN BEGIN
       kxind0[var]=-1.0*(par.nx0-par.nlp_kxind)*kxmin
      ENDIF
      kxstring[var]=STRING(kxind0[var],FORMAT='(F5.2)')
      mode_num=STRTRIM(var+1,2)
      IF var EQ num_nlt_out-1 THEN mode_num='res.'

      en_var_string_ky[var]="<NLT>!Dkx'!N(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
      ;comments_kx[var]="kx spectrum: <NLT>!Dkx'!N(k!Dx!N="+STRING(par.nlp_kxind*kxmin,FORMAT='(F5.2)')+'k!Dy!N='$
      ;       +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
      comments_ky[var]="ky spectrum: <NLT>_kx'(k_x="+kxstring[var]+'k_y='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
     ENDIF ELSE BEGIN
      kxind0[var]=kxmin*(*par.kx_nlt_ind)[var]
      IF (*par.kx_nlt_ind)[var] GT par.nx0/2 THEN BEGIN
        kxind0[var]=-1.0*(par.nx0-(*par.kx_nlt_ind)[var])*kxmin
      ENDIF
      kxstring[var]=STRING(kxind0[var],FORMAT='(F5.2)')
      en_var_string_ky[var]="<NLT>!Dkx'!N(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
      comments_ky[var]="<NLT>kx'(k_x="+kxstring[var]+'ky='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
     ENDELSE
  ENDFOR

   ENDIF ELSE BEGIN
    kx_index=(*i).kx_data
    FOR var=0,num_nlt_out-1 DO BEGIN

     IF par.nlp_gdt THEN BEGIN
      mode_num=STRTRIM(var+1,2)
      IF var EQ num_nlt_out-1 THEN mode_num='res.'

      en_var_string_ky[var]="NLT(k'!Dx!N="+STRING(kxmin*kx_index,'(F5.2)')+")(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
      comments_ky[var]="ky spectrum: NLT(k'_x="+STRING(kxmin*kx_index,'(F5.2)')+")(k_x="+kxstring[var]+'k_y='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
     ENDIF ELSE BEGIN
      en_var_string_ky[var]="NLT(k'!Dx!N="+STRING(kxmin*kx_index,'(F5.2)')+")(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
      comments_ky[var]="NLT(k'x="+STRING(kxmin*kx_index,'(F5.2)')+")(kx="+kxstring[var]+'ky='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
     ENDELSE
    ENDFOR
    kx_index=kx_index+(par.nkx0/2-1)
    data_ky=FLTARR(2*par.nky0-1,num_nlt_out)
    data_ky[*,*]=data_full[kx_index,*,*]
   ENDELSE  

   IF (*i).ky_data EQ 1000 THEN BEGIN
    data_kx = time_avg((*i).data_kx,/avg,fwd=0, tarr=tarr)
    FOR var=0,num_nlt_out-1 DO BEGIN

     IF par.nlp_gdt THEN BEGIN
      mode_num=STRTRIM(var+1,2)
      IF var EQ num_nlt_out-1 THEN mode_num='res.'

      en_var_string_kx[var]="<NLT>!Dky'!N(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
      comments_kx[var]="<NLT>ky'(kx="+kxstring[var]+'ky='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
     ENDIF ELSE BEGIN
      en_var_string_kx[var]="<NLT>!Dky'!N(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
      comments_kx[var]="<NLT>ky'(kx="+kxstring[var]+'ky='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
     ENDELSE

    ENDFOR
   ENDIF ELSE BEGIN

    ky_index=(*i).ky_data
    FOR var=0,num_nlt_out-1 DO BEGIN

     IF par.nlp_gdt THEN BEGIN
      mode_num=STRTRIM(var+1,2)
      IF var EQ num_nlt_out-1 THEN mode_num='res.'

      en_var_string_kx[var]="NLT(k'!Dy!N="+STRING(par.kymin*ky_index,'(F5.2)')+")(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
      comments_kx[var]="NLT(k'y="+STRING(par.kymin*ky_index,'(F5.2)')+")(kx="+kxstring[var]+'ky='$
             +STRING(par.nlp_kyind*par.kymin,FORMAT='(F5.2)')+')'+'(POD n='+mode_num+')'
     ENDIF ELSE BEGIN
      en_var_string_kx[var]="NLT(k'!Dy!N="+STRING(par.kymin*ky_index,'(F5.2)')+")(k!Dx!N="+kxstring[var]+'k!Dy!N='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
      comments_kx[var]="NLT(k'y="+STRING(par.kymin*ky_index,'(F5.2)')+")(kx="+kxstring[var]+'ky='$
             +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
     ENDELSE

    ENDFOR
    ky_index=ky_index+(2*par.nky0-1)/2
    data_kx=FLTARR(par.nkx0,num_nlt_out)
    print,"data_full",size(data_full)
    print,"data_kx",size(data_kx)
    data_kx[*,*]=data_full[*,ky_index,*]

   ENDELSE  

    count = 0
  
    x_as_y = 0
      
    ; identical scaling of x and y axes for sxsy, kxky
    x_as_y = 1
  
    x_axis = "k'!Dx!N "
    y_axis = "k'!Dy!N "
    
    xtitle = '!6' + x_axis + get_var_string(/rhostr)
    ytitle = '!6' + y_axis + get_var_string(/rhostr)
  
    csize = 1.41 - (num_nlt_out < 4) * 0.125
    cthick = 3.0
     
    colorarr = [1,3,2,4,6,8] ; to be color consistent with older plots
    set_output, diag, /ps, multi=[0,1,1], ysize=20, xsize=14.4,$
          suffix=out_suffix
  
    time=store_step((*i).time_id,/get)
    ntime=N_ELEMENTS(time)
    ;print,"time",*time[0],*time[ntime-1]
    tmin=*time[0]
    tmax=*time[ntime-1]
    trange='t=['+STRING(tmin,format='(f6.2)')+','+STRING(tmax,format='(f6.2)')+']'
  
  ;Output plot file here.


    IF (*i).ky_range[0] EQ 999 THEN BEGIN
      xrange=[MIN(yind),MAX(yind)]
    ENDIF ELSE BEGIN
      xrange=(*i).ky_range
    ENDELSE

    FOR var = 0, num_nlt_out - 1 DO BEGIN    ; loop over all variables to plot
      PLOT, yind, data_ky[*,var], COLOR=colorarr[0], XTITLE=ytitle, $
          XLOG=0,  $
          XTICKLEN=1.0, POSITION=[0.15,0.4,0.9,0.9], $
          YLOG=0,  YRANGE=yrange, $
          TITLE=en_var_string_ky[var], XRANGE=xrange 

      XYOUTS, 0.15,0.01,/NORMAL, COLOR=1,$
        trange, CHARSIZE=1.5*csize,$
        CHARTHICK=cthick

    ENDFOR
  
    IF (*i).kx_range[0] EQ 999 THEN BEGIN
      xrange=[MIN(xind),MAX(xind)]
    ENDIF ELSE BEGIN
      xrange=(*i).kx_range
    ENDELSE
    FOR var = 0, num_nlt_out - 1 DO BEGIN    ; loop over all variables to plot
  
      PLOT, xind, data_kx[*,var], COLOR=colorarr[0], XTITLE=xtitle, $
          XLOG=0,  $
          XTICKLEN=1.0, POSITION=[0.15,0.4,0.9,0.9], $
          YLOG=0,  YRANGE=yrange, $
          TITLE=en_var_string_kx[var], XRANGE=xrange 

      XYOUTS, 0.15,0.01,/NORMAL, COLOR=1,$
        trange, CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
  
    ENDFOR
  
    ; Output ascii file for kx spectra
    FOR var = 0, num_nlt_out - 1 DO BEGIN     ; loop over all variables to plot
      ;data_now=data_kx[*,var]
      set_output, diag,$
          dat=[[(*i).xind],[data_kx[*,var]]],append=(var GT 0),$
          suffix=out_suffix,commentlines=comments_kx[var]
    ENDFOR
  
    ; Output ascii file for ky spectra
    FOR var = 0, num_nlt_out - 1 DO BEGIN     ; loop over all variables to plot
      set_output, diag,$
          dat=[[(*i).yind],[data_ky[*,var]]], /append,$
          suffix=out_suffix,commentlines=comments_ky[var]
    ENDFOR
  
  
    set_output, diag, /reset, eps=(series.step_count EQ 1),$
          suffix=out_suffix
  
    ;Output data file
    ;PTR_FREE, time, data
  
    END
