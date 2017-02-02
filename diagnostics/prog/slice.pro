FUNCTION slice_info

  RETURN, {$
    type      : 'mom',$
    title     : 'Slices',$
    help_text : ['Averages over squared absolute of a variable and '+$
                 'plots the result in 1D or 2D (or prints the '+$
                 'resulting 0D value)' ],$
    ext_vars  : [['var','0','variable index (see variable list); '+$
                  'default: 0'],$
		 ['xind','0','kx or x index; default: -1 (average); '+$
                  '-2 for all'],$
    	    	 ['yind','0','ky or y index; default: -1 (average); '+$
                  '-2 for all'],$
    	    	 ['zind','0','z index; default: nz0/2; -1 for '+$
                  'average; -2 for all'],$
		 ['kx2x','1','toggle x: real space/Fourier space'],$
                 ['ky2y','1','toggle y: real space/Fourier space']]}

END

;######################################################################

PRO slice_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var) LT 1 THEN *e.var = 0
  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = -1
  IF N_ELEMENTS(*e.yind) LT 1 THEN *e.yind = -1
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = par.nz0 / 2
  IF N_ELEMENTS(*e.kx2x) NE 1 THEN *e.kx2x = 0
  IF N_ELEMENTS(*e.ky2y) NE 1 THEN *e.ky2y = 0

  IF NOT *e.kx2x AND (par.x_local EQ 0) THEN printerror, $
    'Warning: radial FFT is not recommended for nonlocal version'
  
  IF (NOT par.y_local) AND ((*e.yind)[0] EQ -1 OR (*e.zind)[0] EQ -1) $
       AND (NOT (*e.ky2y)) THEN print,$
       'WARNING: FFT to ky in y global: y-z averages are performed without Jacobian!'

  if (not par.y_local) and *e.ky2y then begin
     nkx0 = par.nx0/2+par.nx0/2 mod 2
     ;in fourier mode we have the positive kx only
     IF (*e.xind)[0] EQ -2 THEN *e.xind = *e.kx2x ? INDGEN(par.nx0) : $
       INDGEN(nkx0)
  endif else begin
     ;in fourier mode we have positive and negative kx
     IF (*e.xind)[0] EQ -2 THEN *e.xind = *e.kx2x ? INDGEN(par.nx0) : $
       SHIFT(INDGEN(par.nx0),par.nx0/2-1)
  endelse

  IF (*e.yind)[0] EQ -2 THEN *e.yind = *e.ky2y ? INDGEN(par.ny0) : $
    INDGEN(par.nky0)
  IF (*e.zind)[0] EQ -2 THEN *e.zind = INDGEN(par.nz0)

  i = set_internal_vars(diag,{$
    var      : *e.var,$
    xind     : *e.xind,$
    yind     : *e.yind,$
    zind     : *e.zind,$
    kx2x     : *e.kx2x,$
    ky2y     : *e.ky2y,$
    modifier : 3,$ ; default: 3 (abs^2)
    time_id  : PTR_NEW(),$
    data_id  : PTR_NEW()})

  CASE 2 * (*i).kx2x + (*i).ky2y OF
    0 : fft_format, kxky=(*i).var
    1 : fft_format, kxsy=(*i).var
    2 : fft_format, sxky=(*i).var
    3 : fft_format, sxsy=(*i).var
  ENDCASE

END

;######################################################################

FUNCTION data_modifier, diag, data, legend=legend

  COMMON global_vars

  i = (*diag).internal_vars

  IF KEYWORD_SET(legend) THEN BEGIN
    CASE (*i).modifier OF
      0 : RETURN, '!9!!!6' + data + '!9!!!6'
      1 : RETURN, '!6Re(' + data + ')'
      2 : RETURN, '!6Im(' + data + ')'
      3 : RETURN, '!9!!!6' + data + '!9!!!6!U2!N'
      4 : RETURN, '!9!!!6' + data + '!9!!!6'
      ELSE : RETURN, data
    ENDCASE
  ENDIF ELSE BEGIN
    CASE (*i).modifier OF
      0 : RETURN, ABS(TEMPORARY(data))
      1 : RETURN, FLOAT(TEMPORARY(data))
      2 : RETURN, IMAGINARY(TEMPORARY(data))
      3 : RETURN, ABS(TEMPORARY(data))^2
      4 : RETURN, ABS(TEMPORARY(data))
      ELSE : RETURN, TEMPORARY(data)
    ENDCASE
  ENDELSE

END

;######################################################################

PRO slice_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  fft_stat = 2 * (*i).kx2x + (*i).ky2y

  xavg = (*i).xind[0] EQ -1
  yavg = (*i).yind[0] EQ -1
  zavg = (*i).zind[0] EQ -1
  add_kx = par.nx0/2 mod 2
  xind = xavg ? (((not par.y_local) and (fft_stat eq 1)) ? INDGEN(par.nx0/2+add_kx) :$
       INDGEN(par.nx0)) : (*i).xind

  yind = yavg ? (((*i).ky2y ) ? INDGEN(par.ny0) : $
       INDGEN(par.nky0)) : (*i).yind

  zind = zavg ? INDGEN(par.nz0) : (*i).zind

  n_vars = N_ELEMENTS((*i).var)

  data = PTRARR(n_vars,gui.out.n_spec_sel)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    FOR ivar = 0, n_vars - 1 DO BEGIN

       xdim = N_ELEMENTS(xind)
       ydim = N_ELEMENTS(yind)
       zdim = N_ELEMENTS(zind)
       ;ydim = yavg ? ((*i).ky2y ? par.ny0 : par.nky0) : $
       ;       N_ELEMENTS((*i).yind)
       ;zdim = zavg ? par.nz0 : N_ELEMENTS((*i).zind)
       
       CASE fft_stat OF
          0 : datptr = mom[isp,(*i).var[ivar]].kxky
          1 : datptr = mom[isp,(*i).var[ivar]].kxsy
          2 : datptr = mom[isp,(*i).var[ivar]].sxky
          3 : datptr = mom[isp,(*i).var[ivar]].sxsy
       ENDCASE
     
       ; select only required data
       data_sp = REFORM((*datptr)[xind,yind,zind,0],$
                        [xdim,ydim,zdim],/OVERWRITE)
       
       ; now take absolute square or real/imaginary part
       ; note: Fourier averages are of quadratic type
       ;       regardless of the modifier choice
       data_sp = data_modifier(diag,data_sp)

       ; x-z Jacobian for x global runs
       IF NOT par.x_local AND (zavg OR xavg) THEN BEGIN
          data_sp = TEMPORARY(data_sp) * REBIN(REFORM($
                    (*series.geom).jac_norm[xind,zind,0],$
                    [xdim,1,zdim],/OVERWRITE),[xdim,ydim,zdim])
       ENDIF
     
       ; y-z Jacobian for y global runs
       IF NOT par.y_local AND (zavg OR yavg) THEN BEGIN
          if ((*i).ky2y) then begin
             data_sp = TEMPORARY(data_sp) * REBIN(REFORM($
                       (*series.geom).jac_norm[yind,zind,0],$
                       [1,ydim,zdim],/OVERWRITE),[xdim,ydim,zdim])
          ;else: warning printed
          endif
       ENDIF
       
       IF zavg THEN BEGIN
          ; z Jacobian for local runs
          IF (par.x_local and par.y_local) THEN data_sp = TEMPORARY(data_sp) * $
               REBIN(REFORM((*series.geom).jac_norm[zind],$
               [1,1,zdim],/OVERWRITE),[xdim,ydim,zdim])
          
          data_sp = REFORM(TOTAL(TEMPORARY(data_sp),3),$
               [xdim,ydim,1],/OVERWRITE) / par.nz0
          zdim = 1
       ENDIF

       IF yavg THEN BEGIN
          IF (*i).ky2y THEN BEGIN
             data_sp = REFORM(TOTAL(TEMPORARY(data_sp),2),$
                [xdim,1,zdim],/OVERWRITE) / par.ny0
          ENDIF ELSE BEGIN
             data_sp[*,0,*] *= 0.5
             data_sp = 2 * REFORM(TOTAL(TEMPORARY(data_sp),2),$
                                  [xdim,1,zdim],/OVERWRITE)
          ENDELSE
          ydim = 1
       ENDIF

       IF xavg THEN BEGIN
          if (*i).kx2x then begin
             ;x real space integration
             data_sp = REFORM(TOTAL(TEMPORARY(data_sp),1)/par.nx0,$
                  [1,ydim,zdim],/OVERWRITE)
          endif else begin
             if (not par.y_local) and (*i).ky2y then begin
                ;mirror positive kx to negative kx
                data_sp[0,*,*] *= 0.5
                data_sp = REFORM(TOTAL(2.0*TEMPORARY(data_sp),1),$
                      [1,ydim,zdim],/OVERWRITE)
             endif else begin
                ;sum up amplitudes for negative and positive kx
                data_sp = REFORM(TOTAL(TEMPORARY(data_sp),1),$
                      [1,ydim,zdim],/OVERWRITE)
             endelse
          endelse
       ;   data_sp = REFORM(TOTAL(TEMPORARY(data_sp),1),$
       ;                    [1,ydim,zdim],/OVERWRITE)
       ;   IF (*i).kx2x THEN data_sp /= par.nx0
       ;   
          xdim = 1
       ENDIF

       data[ivar,isp] = PTR_NEW(data_sp)

    ENDFOR  ; --- ivar loop
  ENDFOR    ; --- species loop

  (*i).time_id = store_step((*i).time_id,mom_time)
  (*i).data_id = store_step((*i).data_id,data)

END

;######################################################################

PRO slice_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nxind = N_ELEMENTS((*i).xind)
  nyind = N_ELEMENTS((*i).yind)
  nzind = N_ELEMENTS((*i).zind)

  IF (*i).kx2x THEN x_axis = INDGEN(par.nx0)*series.dx-series.lx/2 $
    ELSE x_axis = (*series.kx)[*,0]

  IF (*i).ky2y THEN y_axis = -series.ly/2+INDGEN(par.ny0)*series.dy $
    ELSE y_axis = (*series.ky)

  tmparr = [nxind,nyind,nzind]

  ndims = TOTAL(tmparr GT 1)

  IF ndims GT 0 THEN BEGIN
    choice = WHERE(tmparr GT 1)

    my_label = STRARR(ndims)
    my_axis = PTRARR(ndims)

    FOR j = 0, ndims - 1 DO BEGIN
      CASE choice[j] OF
        0 : BEGIN
              IF ((*i).kx2x) THEN my_label[j] = 'x / !7q!6!Dref!N' $
                ELSE my_label[j] = 'k!Ix!N!7q!6!Dref!N'
              my_axis[j] = PTR_NEW(x_axis[(*i).xind])
;              my_axis[j] = PTR_NEW((*series.geom).q[(*i).xind])
;              my_label[j] = 'q'
            END
        1 : BEGIN
              IF ((*i).ky2y) THEN my_label[j] = 'y / !7q!6!Dref!N' $
                ELSE my_label[j] = 'k!Iy!N!7q!6!Dref!N' 
              my_axis[j] = PTR_NEW(y_axis[(*i).yind])
            END
        2 : BEGIN
              my_axis[j] = PTR_NEW((*par.z)[(*i).zind])
              my_label[j] = 'z/L!Dref!N'
            END
      ENDCASE
    ENDFOR
  ENDIF

  info = (*i).kx2x ? 'x' : 'k!Dx!N'
  IF ((*i).xind)[0] EQ -1 THEN info += ' av.' ELSE IF nxind EQ 1 THEN $
    info += '='+rm0es(x_axis[(*i).xind],prec=3)
  IF (*i).ky2y THEN info += ', y' ELSE info += ', k!Dy!N'
  IF ((*i).yind)[0] EQ -1 THEN info += ' av.' ELSE IF nyind EQ 1 THEN $
    info += '='+rm0es(y_axis[(*i).yind],prec=3)
  IF ((*i).zind)[0] EQ -1 THEN info += ', z avg.' ELSE IF nzind EQ 1 THEN $
    info += ', z='+rm0es((*par.z)[(*i).zind],prec=3)

  IF ndims NE 0 THEN set_output, diag, /ps, coltable=33, charsize=csize
  nsp = gui.out.n_spec_sel
  nvars = N_ELEMENTS((*i).var)
  first = 1

  time = store_step((*i).time_id,/get)
  data = store_step((*i).data_id,/get)

  FOR n = 0, series.step_count - 1 DO BEGIN
   FOR isp = 0, nsp - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    FOR ivar= 0 , nvars - 1 DO BEGIN

      mytitle = data_modifier(diag,get_var_string((*i).var[ivar],/fancy)+$
                 ((*i).var[ivar] LT par.n_fields ? '' : ',' + spec[sp].name),/legend)
        
      CASE ndims OF
         0 : BEGIN
            IF (isp+ivar) EQ 0 THEN PRINT, 'time: ' + rm0es(*time[n])
            namestr = rm_idl_fmt(data_modifier(diag,get_var_string((*i).var[ivar]),/legend))+': '
            PRINT, STRING(namestr,format='(A-16)') + rm0es(*(*data[n])[ivar,isp])
         END
         1 : BEGIN
            IF (ivar EQ 0) THEN BEGIN
               collect_dat = REFORM(*(*data[n])[ivar,isp])
               collect_names = mytitle ;rm_idl_fmt(mytitle)
            ENDIF ELSE BEGIN
               collect_dat = [[collect_dat], [REFORM(*(*data[n])[ivar,isp])]]
               collect_names = [[collect_names],[mytitle]] ;[rm_idl_fmt(mytitle)]]
            ENDELSE

            IF (ivar EQ nvars-1) THEN BEGIN
               slice_plot, collect_dat, *my_axis[0],$
                  xtitle = my_label[0], info=info, title=collect_names,$
                  pos = [0.125,0.1+0.9/nsp*(nsp-1-isp),0.925,1.0-0.9/nsp*isp], $
                  no_erase=(isp GT 0)

               set_output, diag, commentlines=['species: '+spec[sp].name+', time: '+$
                   rm0es(*time[n])+', '+rm_idl_fmt(info)], header=[[rm_idl_fmt(my_label[0])],$
                   [rm_idl_fmt(collect_names)]], $
                   append=(n GT 0), dat=[[*my_axis[0]],[collect_dat]]
               plot_info_str, diag, time=(*time[n])
            ENDIF
	  END
      2 : BEGIN
            slice_plot, REFORM(*(*data[n])[ivar,isp]), *my_axis[0], $
              *my_axis[1], xtitle = my_label[0], ytitle = my_label[1],$
              pos = [0.125,0.1+0.9/(nvars*nsp)*(nvars*nsp-1-nvars*isp-ivar),0.925,1.0-0.9/(nvars*nsp)*(nvars*isp+ivar)], $
              no_erase=(isp+ivar GT 0), info=info, title=mytitle
            set_output, diag, commentlines=['species: '+spec[sp].name+', time: '+$
              rm0es(*time[n]), '2D array '+rm_idl_fmt(mytitle)+'('+rm_idl_fmt(my_label[0])+$
              ','+rm_idl_fmt(my_label[1])+') at '+info],$
              append=(first EQ 0), dat=[[REFORM(*(*data[n])[ivar,isp])]]
            plot_info_str, diag, time=(*time[n])            
	  END
      3 : BEGIN
            slice_plot, REFORM(*(*data[n])[ivar,isp]), *my_axis[0], $
              *my_axis[1], *my_axis[2], xtitle=my_label[0], ytitle=my_label[1],$
              info=info, title=mytitle
            set_output, diag, commentlines=['species: '+spec[sp].name+', time: '+$
              rm0es(*time[n]), '3D array '+rm_idl_fmt(mytitle)+'('+rm_idl_fmt(my_label[0])+$
              ','+rm_idl_fmt(my_label[1])+','+rm_idl_fmt(my_label[2])+') at '+info],$
              append=(first EQ 0), dat=[[REFORM(*(*data[n])[ivar,isp])]]
            plot_info_str, diag, time=(*time[n])
	  END
      ELSE : PRINT, ndims, ' dimensions not supported'
     ENDCASE

     first = 0
     PTR_FREE, (*data[n])[ivar,isp]

     ENDFOR  ; --- ivar loop
   ENDFOR ; --- species loop
  ENDFOR

  PTR_FREE, time, data

  IF ndims NE 0 THEN BEGIN
    PTR_FREE, my_axis

    set_output, diag, /reset
  ENDIF

END

;######################################################################

PRO slice_plot, dat, xaxis, yaxis, zaxis, xtitle=xtitle, ytitle=ytitle,$
  title=title, info=info, pos=pos, no_erase=no_erase, oplot=oplot, $
  norm = norm


  IF NOT KEYWORD_SET(xtitle) THEN xtitle = ''
  IF NOT KEYWORD_SET(ytitle) THEN ytitle = ''
  IF NOT KEYWORD_SET(title) THEN title = ''
  IF NOT KEYWORD_SET(info) THEN info = ''
  IF NOT KEYWORD_SET(pos) THEN pos = [0.0,0.0,1.0,1.0]
  IF N_ELEMENTS(no_erase) NE 1 THEN no_erase = 1

  COMMON global_vars

  datmax = MAX(dat,min=datmin)
  dpos_x = pos[2] - pos[0]
  dpos_y = pos[3] - pos[1]
  dy = 0.15

  color = INDGEN(128)+1
  color[0] = 4
  color[5] = 1

  IF N_ELEMENTS(dat) EQ 1 THEN BEGIN
    ; do nothing
  ENDIF ELSE IF N_ELEMENTS(yaxis) LT 1 THEN BEGIN
    LOADCT, 41, FILE='internal/colortable.tbl'

    ndims = SIZE(dat,/n_dimensions)
    IF ndims GT 1 THEN nvars = N_ELEMENTS(dat[0,*]) $
    ELSE nvars = 1

    mydat = FLTARR(N_ELEMENTS(xaxis),ndims)
    mydat[*,*] = dat ;copy data to allow for interference free modifications
   
    IF (ndims GT 1) THEN BEGIN
       mytitle = 'norm.: '
       FOR ivar = 0, nvars-1 DO BEGIN
          mydat[*,ivar] /= MAX(ABS(mydat[*,ivar]))
          mytitle += title[ivar]+' '
       ENDFOR
    ENDIF ELSE BEGIN
       nvars = 1
       mytitle = title
    ENDELSE

    datmax = MAX(mydat,min=datmin)

    PLOT, xaxis, mydat[*,0],$
          XSTYLE=1,YSTYLE=1, XTITLE=xtitle,YTITLE=ytitle,$
          YRANGE=[datmin,datmax],$
          TITLE=mytitle, color=1, /NODATA, NOERASE=no_erase,$
          POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,$
                    pos[0]+0.95*dpos_x,pos[3]-dy*dpos_y]

    FOR ivar = 0, nvars - 1 DO $
       OPLOT, xaxis, mydat[*,ivar], color=color[ivar]

    IF PRODUCT(!Y.CRANGE) LT 0 THEN OPLOT, !X.CRANGE, [0,0],color=1
  ENDIF ELSE IF N_ELEMENTS(zaxis) LT 1 THEN BEGIN
    LOADCT, 33, FILE='internal/colortable.tbl'  
    lev_col = contour_levels([datmin,datmax],c_levels)            
    c_levels = 20
    CONTOUR,dat[*,*], xaxis,yaxis, NOERASE=no_erase,$
      LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
      XSTYLE=1,YSTYLE=1, XTITLE=xtitle,YTITLE=ytitle,$
      POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,$
        pos[0]+0.9*dpos_x,pos[3]-dy*dpos_y],$
      title=title
    plot_colorbar, lev_col, charsize=0.75*!P.CHARSIZE,prec=2,$
      POSITION=[pos[0]+0.95*dpos_x,pos[1]+dy*dpos_y,$
        pos[2],pos[3]-dy*dpos_y]

 ;  LOADCT, 41, FILE='internal/colortable.tbl'

 ;  SURFACE, dat[*,*],xaxis, yaxis, $
 ;    COLOR=1, /XSTYLE,/YSTYLE, $ ;/LEGO, $
 ;    POSITION=[0.1,0.1,0.95,0.475],$
 ;    YTITLE=ytitle, XTITLE=xtitle ;, TITLE=title
  ENDIF ELSE BEGIN
    LOADCT, 33, FILE='internal/colortable.tbl'
    print, 'starting slicer ...'
    datptr = PTR_NEW(dat)
    SLICER3,datptr,/MODAL
    PTR_FREE, datptr
  ENDELSE

  XYOUTS, 0.005, 0.025, info, COLOR=1, CHARSIZE=1.0, /NORMAL

END
