FUNCTION Derivprof2D_info

  RETURN, {$
    type      : 'mom global',$
    title     : 'Deriv. profile (2D)',$
    help_text : ['Displays the 0th, 1st and 2nd radial derivative '+$
                 'of the flux surface average of the chosen variable'+$
                 'vs. x and t'],$
    ext_vars  : [['var','0','index of the variable; default: 0'],$
                 ['xind','0','indices of x grid points to be plotted; '+$
                  'default: 0:nx0-1']]}

END

;######################################################################

PRO Derivprof2D_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var) LT 1 THEN *e.var = 0
  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = INDGEN(par.nx0)
  nxind  = N_ELEMENTS(*e.xind)

  derivxgrid = (- series.lx / 2.0 + (*e.xind) * series.dx) 
  xgrid = derivxgrid
  IF (NOT par.x_local) THEN $
    xgrid = xgrid*par.rhostar + par.x0

  IF (par.x_local) THEN BEGIN
     jac_fac = REBIN(REFORM((*series.geom).jacobian,[1,par.nz0]),$
        [par.nx0,par.nz0])/TOTAL((*series.geom).jacobian) 
  ENDIF ELSE jac_fac = (*series.geom).jacobian * $
    REBIN(1.0/TOTAL((*series.geom).jacobian,2),[par.nx0,par.nz0])

  ; for comparison with nrg (volume average):
  ;*e.xind = INDGEN(par.nx0)
  ;nxind = par.nx0
  ;jac_fac = (*series.geom).jacobian / TOTAL((*series.geom).jacobian) * nxind

  i = set_internal_vars(diag,{$
    xind        : *e.xind,$
    var         : *e.var,$
    derivxgrid  : derivxgrid,$
    xgrid       : xgrid,$
    jac_fac     : jac_fac,$
    time_id     : PTR_NEW(),$
    profx_id    : PTR_NEW()})

  fft_format, sxky=(*i).var

END

;######################################################################

PRO Derivprof2D_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nxind = N_ELEMENTS((*i).xind)
  profx = FLTARR(nxind,gui.out.n_spec_sel,3,/NOZERO)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]
    
    ; calculate flux surface average
    fluxsrfavg = TOTAL((REFORM((*mom[isp,(*i).var].sxky)[(*i).xind,0,*])*$
                     (*i).jac_fac[(*i).xind,*]),2)
    ; 0'th derivative
    profx[*,isp,0] = fluxsrfavg
    ; first (radial) derivative
    profx[*,isp,1] = -DERIV((*i).derivxgrid,fluxsrfavg)
    ; second (radial) derivative
    profx[*,isp,2] = DERIV((*i).derivxgrid,profx[*,isp,1])    
  ENDFOR                          ; species loop

  (*i).time_id = store_step((*i).time_id,mom_time)
  (*i).profx_id = store_step((*i).profx_id,profx)

END

;######################################################################

PRO Derivprof2D_plot2D, data, xaxis_in, yaxis_in, xtitle=xtitle_in, $
  ytitle=ytitle_in, title=title

  IF NOT KEYWORD_SET(xtitle) THEN xtitle = ''
  IF NOT KEYWORD_SET(ytitle) THEN ytitle = ''
  IF NOT KEYWORD_SET(ztitle) THEN ztitle = ''
  IF NOT KEYWORD_SET(title)  THEN title = ''

  LOADCT, 33, FILE='internal/colortable.tbl'
;  LOADCT, 0

  c_levels = 128
  lev_col = contour_levels([MIN(data),MAX(data)],c_levels,/no_fix_zero)

  csize = 1.4
  cthick = 3.0

  transpose = 0 ;1

  IF (transpose) THEN BEGIN
    data = TRANSPOSE(data)
    xaxis = yaxis_in
    yaxis = xaxis_in
    xtitle = ytitle_in
    ytitle = xtitle_in
  ENDIF ELSE BEGIN
    xaxis = xaxis_in
    yaxis = yaxis_in
    xtitle = xtitle_in
    ytitle = ytitle_in
  ENDELSE

  CONTOUR, data, xaxis, yaxis, $
    LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL, $
    /NORMAL, /XSTYLE, /YSTYLE, XTITLE=xtitle, YTITLE=ytitle, $
    YMARGIN=[4,6], CHARTHICK=cthick


  plot_colorbar, lev_col, POSITION=[0.165,0.85,0.95,0.925], $
    CHARSIZE=0.75*csize
  XYOUTS, 0.4, 0.95, title, COLOR=1, /NORMAL, $
    CHARSIZE=2*csize, CHARTHICK=cthick

END

;######################################################################

PRO Derivprof2D_plot1D, data, xaxis, yaxis, xtitle=xtitle, $
  title=title, ytitle=ytitle

  IF NOT KEYWORD_SET(xtitle) THEN xtitle = ''
  IF NOT KEYWORD_SET(title)  THEN title = ''

  csize = 1.4
  cthick = 1.4

  FOR t=0,N_ELEMENTS(yaxis)-1 DO BEGIN

    PLOT, xaxis, data[*,0,t], $
      /XSTYLE, /YSTYLE, XTITLE=xtitle, YTITLE=title, $
      CHARTHICK=cthick, color=1, /NODATA
    OPLOT, !X.CRANGE, [0.0,0.0], color=1

    OPLOT, xaxis, data[*,0,t], color=4
    OPLOT, xaxis, data[*,1,t], color=3
    OPLOT, xaxis, data[*,2,t], color=2
  ENDFOR

  XYOUTS, 0.45, 0.95, title, COLOR=1, /NORMAL, $
    CHARSIZE=2*csize, CHARTHICK=cthick

END

 ;######################################################################

PRO Derivprof2D_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  eps = 0

  time = store_step((*i).time_id,/get,/reg_array)
  profx = store_step((*i).profx_id,/get,/reg_array)

  n_xind = N_ELEMENTS((*i).xind)

  rho_str = get_var_string(/rhostr)
  IF (par.x_local) THEN xtitle = '!6x / '+rho_str $
  ELSE xtitle = '!6x / a'


  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
  
   sp=(*gui.out.spec_select)[isp]

   timeline = '!6t='+rm0es(gui.out.start_t)+'-'+$
     rm0es(gui.out.end_t)+' '+get_var_string(1,/time,/ounit)+$
     ', ' + STRTRIM(STRING(series.step_count),2) + $
     ', ' + spec[sp].name

   set_output, diag, sp, ps=(eps EQ 0), eps=eps,multi=[0,1,1], $
     ysize=14.4, xsize=20, coltable=41
   
   IF gui.out.res_steps THEN BEGIN
      title = '!N <'+get_var_string((*i).var,/fancy)+'>!IFS!N, '
      title += '-!9d!6!Ix!N <'+get_var_string((*i).var,/fancy)+'>!IFS!N, '
      title += '-!9d!6!Ix!E2!N <'+get_var_string((*i).var,/fancy)+'>!IFS!N'
      Derivprof2D_plot1D, reform(profx[*,isp,*,*]), (*i).xgrid, time,$
         xtitle=xtitle, ytitle='t / ['+get_var_string(0,/time,/ounit)+']',$
         title=title
   ENDIF ELSE BEGIN
       Derivprof2D_plot2D, reform(profx[*,isp,0,*]), (*i).xgrid, time,$
         xtitle=xtitle, ytitle='t / ['+get_var_string(0,/time,/ounit)+']',$
         title='!N <'+get_var_string((*i).var,/fancy)+'>!IFS!N'

       Derivprof2D_plot2D, reform(profx[*,isp,1,*]), (*i).xgrid, time,$
         xtitle=xtitle, ytitle='t / ['+get_var_string(0,/time,/ounit)+']',$
         title='-!9d!6!Ix!N <'+get_var_string((*i).var,/fancy)+'>!IFS!N'

       Derivprof2D_plot2D, reform(profx[*,isp,2,*]), (*i).xgrid, time,$
         xtitle=xtitle, ytitle='t / ['+get_var_string(0,/time,/ounit)+']',$
         title='-!9d!6!Ix!E2!N <'+get_var_string((*i).var,/fancy)+'>!IFS!N'
   ENDELSE

;   set_output, diag, sp, header = ['x/a',header],$
;     dat=[[(*i).xgrid],[avprofx[*,*,isp]]]
     
   set_output, diag, sp, /reset,eps=eps

  ENDFOR ;-- isp

END
