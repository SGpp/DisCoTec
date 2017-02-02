FUNCTION srcmom1_info

  RETURN, {$
    type      : 'srcmom',$
    title     : 'sources moments',$
    help_text : ['velocity space moments of sources'],$
    ext_vars  : [['vars','0','variables to be plotted '+$
    	    	  '(0=ck_heat, 1=ck_part, 2=heat src); '+$
                  'default: all'],$
                 ['moms','0','moments to be plotted '+$
    	    	  '(0=0th moment, 1=vpar moment, 2=v^2 moment)'],$
                 ['xind','0','x indices to be plotted; default: all']$
                 ]}

END

;######################################################################

PRO srcmom1_init, diag

  COMMON global_vars
  
  e = (*diag).external_vars
  
  IF N_ELEMENTS(*e.vars) LT 1 THEN *e.vars = INDGEN(3)
  IF N_ELEMENTS(*e.moms) LT 1 THEN *e.moms = INDGEN(3)
  IF N_ELEMENTS(*e.xind) LT 1 THEN *e.xind = INDGEN(par.nx0)

  var_names = ['ck_heat', 'ck_part','heat src']
  mom_names = ['0th', 'v!I!9#!6!N', 'v!E2!N']

  i = set_internal_vars(diag,{$
    vars       : *e.vars,$
    moms       : *e.moms,$
    xind       : *e.xind,$
    var_names  : var_names[*e.vars],$
    mom_names  : mom_names[*e.moms],$
    xaxis      : par.x0 + (-series.lx/2.+(*e.xind)*series.dx)*par.rhostar,$
    average_id : PTR_NEW()})

  series.request_srcmom[*e.vars] = 1

END

;######################################################################

PRO srcmom1_loop, diag

  COMMON global_vars
  
  i = (*diag).internal_vars

  n_vars = N_ELEMENTS((*i).vars)
  n_moms = N_ELEMENTS((*i).moms)
  n_xind = N_ELEMENTS((*i).xind)
  n_specs = gui.out.n_spec_sel

  data = FLTARR(n_xind,n_moms,n_vars,n_specs,/NOZERO)
  data = REFORM(data,[n_xind,n_moms,n_vars,n_specs],/OVERWRITE)

  FOR isp = 0, n_specs - 1 DO BEGIN
    FOR ivar = 0, N_ELEMENTS((*i).vars) - 1 DO BEGIN
      data[*,*,ivar,isp] = (*srcmom[isp,(*i).vars[ivar]])$
        [(*i).xind,(*i).moms,0]
    ENDFOR
    srcmom1_plot, diag, data[*,*,*,isp]
  ENDFOR

  (*i).average_id = time_avg((*i).average_id,data,srcmom_time)

END

;######################################################################

PRO srcmom1_plot, diag, dat, title=title, xtitle=xtitle, ytitle=ytitle

  COMMON global_vars

  i = (*diag).internal_vars

  dat = REFORM(dat,[N_ELEMENTS((*i).xind),N_ELEMENTS((*i).moms),$
                    N_ELEMENTS((*i).vars)],/OVERWRITE)

  n_moms = N_ELEMENTS(dat[0,*,0])

  ymin = MIN(dat,max=ymax)

  IF NOT(KEYWORD_SET(title)) THEN title = ''
  IF NOT(KEYWORD_SET(xtitle)) THEN xtitle = 'x/a'

  count = 0
  
  FOR ivar = 0, N_ELEMENTS((*i).vars) - 1 DO BEGIN
    FOR imom = 0, n_moms - 1 DO BEGIN
      IF count EQ 0 THEN BEGIN
        PLOT, (*i).xaxis,dat[*,imom,ivar],/XSTYLE,$
          YRANGE=[ymin,ymax],/YSTYLE,color = count+1,$
          POSITION=[0.15,0.25,0.95,0.95],xtitle=xtitle,$
          title=title
        legendstr = (*i).var_names[ivar]+$
                   ', ('+(*i).mom_names[imom]+')'
      ENDIF ELSE BEGIN
        OPLOT, (*i).xaxis, dat[*,imom,ivar], color = count+1
        legendstr = [legendstr,(*i).var_names[ivar]+$
                   ', ('+(*i).mom_names[imom]+')']
      ENDELSE
      count +=1
    ENDFOR
  ENDFOR

  IF PRODUCT(!Y.CRANGE) LT 0.0 THEN $
     OPLOT, !X.CRANGE, [0.0,0.0], color=1

  plot_legend, INDGEN(count)+1, legendstr, $
         x_offset=[0.15,0.0,0.95,0.12], per_line = 3

  timeline = ''
  XYOUTS, 0.15, 0.125, /NORMAL, timeline, COLOR=1 

END 

;######################################################################

PRO srcmom1_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  n_vars = N_ELEMENTS((*i).vars)
  n_moms = N_ELEMENTS((*i).moms)
  n_xind = N_ELEMENTS((*i).xind)
  n_specs = gui.out.n_spec_sel
  
  data = REFORM(time_avg((*i).average_id,/avg),$
                [n_xind,n_moms,n_vars,n_specs])

  FOR isp = 0, n_specs - 1 DO BEGIN
     sp=(*gui.out.spec_select)[isp]
   
     set_output, diag, sp, /eps, multi=[0,1,2], coltable=41
     srcmom1_plot, diag, data[*,*,*,isp], title="Time avg'd source/sink moments"

     set_output, diag, sp, /reset, /eps
     
  ENDFOR 

  FOR isp = 0, n_specs - 1 DO BEGIN
     sp=(*gui.out.spec_select)[isp]
     set_output, diag, sp, header=['source_profile'],$
        dat=[[data[*,*,0,isp]],[data[*,*,1,isp]]], append=(isp GT 0)
  ENDFOR 



END
