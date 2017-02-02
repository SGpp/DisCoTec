FUNCTION nltcont_info

  RETURN, {$
    type      : 'nlt',$
    title     : 'NLT contour plots',$
    help_text : ['Plots contour plots for nonlinear energy transfer functions (to do: vrange option) '],$
    ext_vars  : [ ['kx_range','0', '  Specifies plot range in kx e.g., [-0.5,0.5] (Default: max,min)' ],$
    	    	 ['ky_range','0', '  Specifies plot range in ky e.g., [-1.5,1.5] (Default: max,min)' ],$
                 ['tout','0','plot time-resolved nlts']$
                                                                                              ]}


END

;######################################################################

PRO nltcont_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF NOT KEYWORD_SET(*e.kx_range) THEN BEGIN
    print, 'kx_range not set'
    (*e.kx_range)=[999,999]
  ENDIF
  IF NOT KEYWORD_SET(*e.ky_range) THEN BEGIN
    print, 'ky_range not set'
    (*e.ky_range)=[999,999]
  ENDIF
  IF NOT KEYWORD_SET(*e.tout) THEN BEGIN
    (*e.tout)=0
  ENDIF

  IF par.lx EQ 0 THEN BEGIN
    STOP, 'lx cannot be 0'
  ENDIF

  i = set_internal_vars(diag,{$
    kx_range   : *e.kx_range,$
    ky_range   : *e.ky_range,$
    tout      : *e.tout,$
    yind      : FINDGEN(2*par.nky0-1),$
    xind      : FINDGEN(par.nkx0),$
    time_id   : PTR_NEW(),$
    data_tavg : PTR_NEW(),$
    data_id   : PTR_NEW()})

    (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
      (2.0 * !PI / series.lx)
    (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)

END

;######################################################################

PRO nltcont_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkx0 = par.nkx0
  nky0 = par.nky0

  data=*nlt2d
  (*i).time_id = store_step((*i).time_id,nlt_time)
  (*i).data_id = store_step((*i).data_id,data)
  (*i).data_tavg = time_avg((*i).data_tavg,data,nlt_time)

END

;######################################################################

PRO nltcont_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

    IF (*i).ky_range[0] EQ 999 THEN BEGIN
      yrange=[MIN((*i).yind),MAX((*i).yind)]
    ENDIF ELSE BEGIN
      yrange=(*i).ky_range
    ENDELSE
    IF (*i).kx_range[0] EQ 999 THEN BEGIN
      xrange=[MIN((*i).xind),MAX((*i).xind)]
    ENDIF ELSE BEGIN
      xrange=(*i).kx_range
    ENDELSE

  kxmin=2.0*!PI/par.lx
  num_nlt_out = par.num_nlt_modes
  IF par.nlp_gdt THEN num_nlt_out = par.num_nlt_pod_modes+1

  en_var_string=STRARR(num_nlt_out)
  kxstring=STRARR(num_nlt_out)
  kxind0=FLTARR(num_nlt_out)
  
  FOR var=0,num_nlt_out-1 DO BEGIN
     IF par.nlp_gdt THEN BEGIN
      mode_num=STRTRIM(var+1,2)
      kxind0[var]=kxmin*par.nlp_kxind
      IF par.nlp_kxind GT par.nx0/2 THEN BEGIN
       kxind0[var]=-1.0*(par.nx0-par.nlp_kxind)*kxmin
      ENDIF
      kxstring[var]=STRING(kxind0[var],FORMAT='(F5.2)')
      ;IF kx_ind[f] GT par.nx0/2 THEN kxstring=STRING(-1.0*(par.nx0-kx_ind[f])*kxmin,FORMAT='(F5.2)')
      ;kxstring=STRING(kx_ind[f]*kxmin,FORMAT='(F5.2)')
      IF var EQ num_nlt_out-1 THEN mode_num='res.'
      en_var_string[var]='T!DNL!N(k!Dx!N='+$
         kxstring[var]+',k!Dy!N='+STRING(par.kymin*par.nlp_kyind,FORMAT='(F5.2)')+$
         ')'+'(POD n='+mode_num+')'
     ENDIF ELSE BEGIN
       kxind0[var]=kxmin*(*par.kx_nlt_ind)[var]
       IF (*par.kx_nlt_ind)[var] GT par.nx0/2 THEN BEGIN
         kxind0[var]=-1.0*(par.nx0-(*par.kx_nlt_ind)[var])*kxmin
       ENDIF
       kxstring[var]=STRING(kxind0[var],FORMAT='(F5.2)')
       en_var_string[var]="T!D NL !N(k!Dx!N="+kxstring[var]+'k!Dy!N='$
            +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
     ENDELSE
  ENDFOR

  n_vars = num_nlt_out ;par.num_nlt_modes
  count = 0

  aspect_ratio = (*i).yind[N_ELEMENTS((*i).yind)-1] / $
    (*i).xind[N_ELEMENTS((*i).xind)-1]
    
  xysize = [17,25]

  x_axis = "k'!Dx!N "
  y_axis = "k'!Dy!N "
  
  xtitle = '!6' + x_axis + get_var_string(/rhostr)
  ytitle = '!6' + y_axis + get_var_string(/rhostr)

  csize = 1.41 
  cthick = 3.0

  c_levels = 20
  cont_minmax = FLTARR(2,n_vars)
  
  time = store_step((*i).time_id,/get)
  ;help,time
  ;print,*time[0]
  data = store_step((*i).data_id,/get)

  data_tavg = time_avg((*i).data_tavg,/avg,fwd=0, tarr=tarr)
   
 IF (*i).tout THEN BEGIN
 FOR var = 0, n_vars - 1 DO BEGIN    ; loop over all variables to plot
  first = 1
  IF par.nlp_gdt THEN BEGIN
    out_suffix = 'POD'+STRING(var,FORMAT='(I2.2)')
  ENDIF ELSE BEGIN
    out_suffix = STRING(var,FORMAT='(I2.2)')
  ENDELSE
  set_output, diag, coltable=33, $
    xsize=xysize[0], ysize=xysize[1],suffix=out_suffix, charsize=csize,$
    ps=(series.step_count GT 1), eps=(series.step_count EQ 1)
    FOR n = 0, series.step_count - 2 DO BEGIN
      timechar='t='+STRING(*time[n],'(F8.4)')

      lev_col = contour_levels([MIN((*data[n])[*,*,var]),$
        MAX((*data[n])[*,*,var])],c_levels)
      
      CONTOUR, (*data[n])[*,*,var], (*i).xind, (*i).yind, $
        LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
        /XSTYLE, /YSTYLE, $
        XTITLE=xtitle, YTITLE=ytitle, $
        POSITION=[0.1,0.5,0.9,0.9], $
        CHARSIZE=csize, CHARTHICK=cthick, /NORMAL, YTICK_GET=YTICKV, $
	XTICK_GET=XTICKV,XRANGE=xrange,YRANGE=yrange

        AXIS, 0,0, XAX=0,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0
        AXIS, 0,0, XAX=1,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0	  
        AXIS, 0,0, YAX=0,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0
        AXIS, 0,0, YAX=1,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0	  

      XYOUTS, 0.1,0.92,/NORMAL, COLOR=1,$
        en_var_string[var], CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
      XYOUTS, 0.7,0.45,/NORMAL, COLOR=1,$
        timechar, CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
      ;XYOUTS, (*par.kx_nlt_ind)[var]*kxmin,(*par.ky_nlt_ind)[var]*par.kymin,/DATA, COLOR=1,$
      ;  'X', CHARSIZE=1.5*csize,$
      ;  CHARTHICK=cthick
      IF par.nlp_gdt THEN BEGIN
        arrow, kxstring[var]-5*kxmin,par.nlp_kyind*par.kymin-5*par.kymin,$
                kxtring[var],par.nlp_kyind*par.kymin,/DATA,/FILL
      ENDIF ELSE BEGIN
        arrow, kxstring[var]-5*kxmin,(*par.ky_nlt_ind)[var]*par.kymin-5*par.kymin,$
                kxstring[var],(*par.ky_nlt_ind)[var]*par.kymin,/DATA,/FILL
      ENDELSE
      ;arrow, (*par.kx_nlt_ind)[var]*kxmin-5*kxmin,(*par.ky_nlt_ind)[var]*par.kymin-5*par.kymin,$
      ;          (*par.kx_nlt_ind)[var]*kxmin,(*par.ky_nlt_ind)[var]*par.kymin,/DATA,/FILL

      plot_colorbar, lev_col, position=[0.91,0.5,0.93,0.9], charsize=0.75*csize
	
      ; data output	
      ERASE
      IF par.nlp_gdt THEN BEGIN
        out_suffix = 'POD'+STRING(var,FORMAT='(I2.2)')
      ENDIF ELSE BEGIN
        out_suffix = STRING(var,FORMAT='(I2.2)')
      ENDELSE
      set_output, diag, dat=[[REFORM((*data[n])[*,*,var])]],$
        commentlines='time: '+rm0es(*time[n])+ ', variable: '+$
        rm_idl_fmt(en_var_string[var])+'('+rm_idl_fmt(xtitle)+','+$
        rm_idl_fmt(ytitle)+')',append = (first EQ 0),suffix=out_suffix
      first =0 	
    ENDFOR ;time loop
    ;set_output, diag, /reset, eps=(series.step_count EQ 1),suffix=STRING(var,FORMAT='(I2.2)')
    set_output, diag, /reset, eps=(series.step_count EQ 1),suffix=out_suffix

  ENDFOR ;var loop
  ENDIF

    ;Now output time averaged contour plots
    ;Now output time averaged contour plots
    ;Now output time averaged contour plots
  ntime=N_ELEMENTS(time)
  ;print,"time",*time[0],*time[ntime-1]
  tmin=*time[0]
  tmax=*time[ntime-1]
  trange='t=['+STRING(tmin,format='(f6.2)')+','+STRING(tmax,format='(f6.2)')+']'
  first=1
  IF par.nlp_gdt THEN BEGIN
    out_suffix = 'POD_tavg'
  ENDIF ELSE BEGIN
    out_suffix = '_tavg'
  ENDELSE
  set_output, diag, coltable=33, $
    xsize=xysize[0], ysize=xysize[1],suffix=out_suffix, charsize=csize,$
    ps=(series.step_count GT 1), eps=(series.step_count EQ 1)

  FOR var=0,num_nlt_out-1 DO BEGIN
     IF par.nlp_gdt THEN BEGIN
      mode_num=STRTRIM(var+1,2)
      IF var EQ num_nlt_out-1 THEN mode_num='res.'
      en_var_string[var]='<T!DNL!N>!Dt!N(k!Dx!N='+kxstring[var]$
         +',k!Dy!N='+STRING(par.kymin*par.nlp_kyind,FORMAT='(F5.2)')+$
         ')'+'(POD n='+mode_num+')'
     ENDIF ELSE BEGIN
       en_var_string[var]="<T!D NL!N>!Dt!N(k!Dx!N="+kxstring[var]+'k!Dy!N='$
            +STRING((*par.ky_nlt_ind)[var]*par.kymin,FORMAT='(F5.2)')+')'
     ENDELSE
  ENDFOR

    FOR var = 0, n_vars - 1 DO BEGIN    ; loop over all variables to plot

      lev_col = contour_levels([MIN(data_tavg[*,*,var]),$
        MAX(data_tavg[*,*,var])],c_levels)
      
      CONTOUR, data_tavg[*,*,var], (*i).xind, (*i).yind, $
        LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
        /XSTYLE, /YSTYLE, $
        XTITLE=xtitle, YTITLE=ytitle, $
        POSITION=[0.1,0.5,0.9,0.9], $
        CHARSIZE=csize, CHARTHICK=cthick, /NORMAL, YTICK_GET=YTICKV, $
	XTICK_GET=XTICKV,XRANGE=xrange,YRANGE=yrange

        AXIS, 0,0, XAX=0,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0
        AXIS, 0,0, XAX=1,COLOR=55, XTICKV=XTICKV, XTICKS=N_ELEMENTS(XTICKV)-1, $
	  XTICKNAME = REPLICATE(' ',N_ELEMENTS(XTICKV)), XMINOR=0	  
        AXIS, 0,0, YAX=0,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0
        AXIS, 0,0, YAX=1,COLOR=55, YTICKV=YTICKV, YTICKS=N_ELEMENTS(YTICKV)-1, $
	  YTICKNAME = REPLICATE(' ',N_ELEMENTS(YTICKV)), YMINOR=0	  

      XYOUTS, 0.1,0.92,/NORMAL, COLOR=1,$
        en_var_string[var], CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
      XYOUTS, 0.4,0.01,/NORMAL, COLOR=1,$
        trange, CHARSIZE=1.5*csize,$
        CHARTHICK=cthick
      ;XYOUTS, (*par.kx_nlt_ind)[var]*kxmin,(*par.ky_nlt_ind)[var]*par.kymin,/DATA, COLOR=1,$
      ;  'X', CHARSIZE=1.5*csize,$
      ;  CHARTHICK=cthick
      IF par.nlp_gdt THEN BEGIN
        arrow, kxstring[var]-5*kxmin,par.nlp_kyind*par.kymin-5*par.kymin,$
                kxstring[var],par.nlp_kyind*par.kymin,/DATA,/FILL
      ENDIF ELSE BEGIN
        arrow, kxstring[var]-5*kxmin,(*par.ky_nlt_ind)[var]*par.kymin-5*par.kymin,$
                kxstring[var],(*par.ky_nlt_ind)[var]*par.kymin,/DATA,/FILL
      ENDELSE
;	pro arrow, xt, yt, xh, yh, device=device, data=data, norm=norm, $
;          help=hlp, color=color, linestyle=linestyle, thick=thick, $
;	  fill=fill, length=alen0, width=awid0, shaft=ash0

      plot_colorbar, lev_col, position=[0.91,0.5,0.93,0.9], charsize=0.75*csize

      ERASE
	
      ; data output	

      IF par.nlp_gdt THEN BEGIN
        out_suffix = 'POD_tavg'
      ENDIF ELSE BEGIN
        out_suffix = '_tavg'
      ENDELSE
      set_output, diag, dat=[[REFORM(data_tavg[*,*,var])]],$
        commentlines='time average '+ ', variable: '+$
        rm_idl_fmt(en_var_string[var])+'('+rm_idl_fmt(xtitle)+','+$
        rm_idl_fmt(ytitle)+')',append = (first EQ 0),suffix='_tavg'
      first =0 	
    ENDFOR
    set_output, diag, /reset, eps=(series.step_count EQ 1),suffix=out_suffix


  PTR_FREE, time, data
  
    END

;-------------------------------------------------------------
;+
; NAME:
;       ARROW
; PURPOSE:
;       Draw arrows on screen.
; CATEGORY:
; CALLING SEQUENCE:
;       arrow, xt, yt, xh, yh
; INPUTS:
;       xt, yt = x,y of arrow tail.     in
;       xh, yh = x,y of arrow head.     in
; KEYWORD PARAMETERS:
;       Keywords:
;         /DATA means use data coordinates (def).
;         /NORM means use normalized coordinates.
;         /DEVICE means use device coordinates.
;         COLOR=c arrow outline color (def=255).
;           Make same as FILL for no outline.
;         THICK=t arrow outline thickness (def=0).
;         LINESTYLE=s arrow outline line style (def=0).
;         FILL=f set arrow fill color (0-255, def = -1 = no fill).
;         LENGTH=l  Arrow head length in % plot width (def = 3).
;         WIDTH=w  Arrow head width in % plot width (def = 2).
;         SHAFT=s  Arrow shaft width in % plot width (def = 1).
; OUTPUTS:
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner, 21 Sep, 1990
;
; Copyright (C) 1990, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	pro arrow, xt, yt, xh, yh, device=device, data=data, norm=norm, $
          help=hlp, color=color, linestyle=linestyle, thick=thick, $
	  fill=fill, length=alen0, width=awid0, shaft=ash0
 
	if (n_params(0) lt 4) or keyword_set(hlp) then begin
	  print,' Draw arrows on screen.'
	  print,' arrow, xt, yt, xh, yh'
	  print,'   xt, yt = x,y of arrow tail.     in'
	  print,'   xh, yh = x,y of arrow head.     in'
	  print,' Keywords:'
	  print,'   /DATA means use data coordinates (def).'
	  print,'   /NORM means use normalized coordinates.'
	  print,'   /DEVICE means use device coordinates.'
	  print,'   COLOR=c arrow outline color (def=255).'
	  print,'     Make same as FILL for no outline.'
	  print,'   THICK=t arrow outline thickness (def=0).'
	  print,'   LINESTYLE=s arrow outline line style (def=0).'
	  print,'   FILL=f set arrow fill color (0-255, def = -1 = no fill).'
	  print,'   LENGTH=l  Arrow head length in % plot width (def = 3).'
	  print,'   WIDTH=w  Arrow head width in % plot width (def = 2).'
	  print,'   SHAFT=s  Arrow shaft width in % plot width (def = 1).'
	  return
	endif
 
	;---------  Make sure parameters set  ---------
	if n_elements(device) eq 0 then device = 0
	if n_elements(data) eq 0 then data = 0
	if n_elements(norm) eq 0 then norm = 0
	if (device+data+norm) eq 0 then data = 1
	if n_elements(alen0) eq 0 then alen0 = 3.0
	if n_elements(awid0) eq 0 then awid0 = 2.0
	if n_elements(ash0) eq 0 then ash0 = 1.0
	if n_elements(color) eq 0 then color = !p.color
	if n_elements(linestyle) eq 0 then linestyle = !p.linestyle
	if n_elements(thick) eq 0 then thick = !p.thick
	if n_elements(fill) eq 0 then fill = -1
 
	alen = alen0/100.
	awid2 = awid0/100./2.
	ash2 = ash0/100./2.
	;  Arrows have a problem in normalized coordinates.  A length of 1 unit
	;  will in general be different in X and Y.  This may be handled easily
	;  by converting to an isotropic coordinate system, computing the arrow,
	;  then converted back to norm. to plot.  The factor ff is the x to y
	;  shape ratio that does this conversion.
	ff = float(!d.x_size)/float(!d.y_size)	; Fudge fact to cor norm coord.
 
 
	;---------  Convert arrow to normalized coordinates  --------
	if keyword_set(data) then begin
	  t = convert_coord([xt,xh],[yt,yh],/data,/to_norm)
	  x1 = (t(0,*))(0)
	  x2 = (t(0,*))(1)
	  y1 = (t(1,*))(0)
	  y2 = (t(1,*))(1)
	endif
	if keyword_set(device) then begin
	  t = convert_coord([xt,xh],[yt,yh],/dev,/to_norm)
	  x1 = (t(0,*))(0)
	  x2 = (t(0,*))(1)
	  y1 = (t(1,*))(0)
	  y2 = (t(1,*))(1)
	endif
	if keyword_set(norm) then begin
	  x1 = xt
	  x2 = xh
	  y1 = yt
	  y2 = yh
	endif
	x1 = x1*ff
	x2 = x2*ff
 
	plots, [x1,x1]/ff, [y1,y1], /norm	; Plot one.
	dx = x2 - x1  & dy = y2 - y1	; Set up arrow.
	m1 = sqrt(dx^2 + dy^2)>.1
	u1x = dx/m1  & u1y = dy/m1		; Unit vector along arrow.
	u2x = -u1y  & u2y = u1x		; Unit vector across arrow.
	x2b = x2 - alen*u1x  & y2b = y2 - alen*u1y	  ; Midpt back of head.
	hx1 = x2b + ash2*u2x  & hy1 = y2b + ash2*u2y     ; ARROW HEAD BACK.
	hx2 = x2b - ash2*u2x  & hy2 = y2b - ash2*u2y     ; ARROW HEAD BACK.
	hx3 = x2b + awid2*u2x  & hy3 = y2b + awid2*u2y   ; ARROW HEAD BACK.
	hx4 = x2b - awid2*u2x  & hy4 = y2b - awid2*u2y   ; ARROW HEAD BACK.
	tx1 = x1 + ash2*u2x  & ty1 = y1 + ash2*u2y	     ; ARROW TAIL.
	tx2 = x1 - ash2*u2x  & ty2 = y1 - ash2*u2y	     ; ARROW TAIL.
	xp = [tx1, hx1, hx3, x2, hx4, hx2, tx2, tx1]
	yp = [ty1, hy1, hy3, y2, hy4, hy2, ty2, ty1]
	plots, xp/ff, yp, /norm		; Plot temp arrow.
 
	if fill ne -1 then polyfill, /norm, xp/ff, yp, color=fill
	plots, xp/ff, yp, /norm, color=color, linestyle=linestyle, $
	  thick=thick
 
	end




