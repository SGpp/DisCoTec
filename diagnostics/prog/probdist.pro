FUNCTION probdist_info

  RETURN, {$
    type      : 'mom',$ ;'mom_uni',$
    title     : 'Probability distributions',$
    help_text : ['Plots the probability distribution function of '+$
                 'either a single variable or of a pair of '+$
                 'quantities'],$
    ext_vars  : [['var1','0','first variable; default: phi'],$
                 ['var2','0','second variable; default: none'],$
                 ['range1','0','amplitude range of the first '+$
                  'variable ([min,max]); default: [-30,30]'],$
                 ['range2','0','amplitude range of the second '+$
                  'variable ([min,max]); default: [-30,30]']]}

END

;######################################################################

PRO probdist_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var1) NE 1 THEN *e.var1 = 0
  use_2d = N_ELEMENTS(*e.var2) NE 1 ? 0 : 1

  IF N_ELEMENTS(*e.range1) NE 2 THEN *e.range1 = [-30.0,30.0]
  IF N_ELEMENTS(*e.range2) NE 2 THEN *e.range2 = [-30.0,30.0]

  ; bin size and axes
  bin_size_1 = 1.0
  bin_size_2 = 1.0
  bin_count_1 = ((*e.range1)[1] - (*e.range1)[0]) / bin_size_1 + 1
  bin_count_2 = ((*e.range2)[1] - (*e.range2)[0]) / bin_size_2 + 1
  x_axis = (*e.range1)[0] + ((*e.range1)[1] - (*e.range1)[0]) * $
    FINDGEN(bin_count_1) / (bin_count_1 - 1)
  y_axis = (*e.range2)[0] + ((*e.range2)[1] - (*e.range2)[0]) * $
    FINDGEN(bin_count_2) / (bin_count_2 - 1)

  i = set_internal_vars(diag,{$
    use_2d      : use_2d,$
    var1        : *e.var1,$
    var2        : use_2d ? *e.var2 : 0,$
    bin_size_1  : bin_size_1,$
    bin_size_2  : bin_size_2,$
    bin_count_1 : bin_count_1,$
    bin_count_2 : bin_count_2,$
    x_axis      : x_axis,$
    y_axis      : y_axis,$  
    def_min_1   : (*e.range1)[0],$
    def_max_1   : (*e.range1)[1],$
    def_min_2   : (*e.range2)[0],$
    def_max_2   : (*e.range2)[1],$
    z_values_id : PTR_NEW(),$
    abs_min_1   : LONARR(gui.out.n_spec_sel),$
    abs_max_1   : LONARR(gui.out.n_spec_sel),$
    abs_min_2   : LONARR(gui.out.n_spec_sel),$
    abs_max_2   : LONARR(gui.out.n_spec_sel)})

  var_req = use_2d ? [(*i).var1,(*i).var2] : (*i).var1
  fft_format, sxsy=var_req

END

;######################################################################

PRO probdist_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).use_2d THEN BEGIN ; 2D histogram
    z_values = FLTARR((*i).bin_count_1,(*i).bin_count_2,gui.out.n_spec_sel)

    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      data_1 = REFORM(*mom[isp,(*i).var1].sxsy)
      data_2 = REFORM(*mom[isp,(*i).var2].sxsy)
      (*i).abs_min_1 = MIN(data_1,MAX=data_max)
      (*i).abs_max_1 = data_max
      (*i).abs_min_2 = MIN(data_2,MAX=data_max)
      (*i).abs_max_2 = data_max

      z_values[*,*,isp] += HIST_2D(data_1,data_2,BIN1=(*i).bin_size_1,$
        BIN2=(*i).bin_size_2,MIN1=(*i).def_min_1,MIN2=(*i).def_min_2,$
        MAX1=(*i).def_max_1,MAX2=(*i).def_max_2)
    ENDFOR
  ENDIF ELSE BEGIN ; 1D histogram
    z_values = FLTARR((*i).bin_count_1,gui.out.n_spec_sel)

    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      data = REFORM(*mom[isp,(*i).var1].sxsy)
      (*i).abs_min_1[isp] = MIN(data,MAX=abs_max)
      (*i).abs_max_1[isp] = abs_max
      z_values[*,isp] += HISTOGRAM(data,BINSIZE=(*i).bin_size_1,$
        MIN=(*i).def_min_1,NBINS=(*i).bin_count_1)
    ENDFOR
  ENDELSE

  (*i).z_values_id = time_avg((*i).z_values_id,z_values,mom_time)

END

;######################################################################

PRO probdist_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  z_values_avg = time_avg((*i).z_values_id,/avg)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps

    IF (*i).use_2d THEN BEGIN
      SURFACE, z_values_avg[*,*,isp], (*i).x_axis, (*i).y_axis, $
        TITLE='!6Density function: ', /LEGO, ZTITLE='!6count', $
        XTITLE=get_var_string((*i).var1,/fancy), $
        YTITLE=get_var_string((*i).var2,/fancy), $
        COLOR=1, XCHARSIZE=2.0, YCHARSIZE=2.0, /XSTYLE, /YSTYLE

      PRINT, (*diag).name + ' info: (species ' + spec[sp].name + ')'
      PRINT, 'Absolute amplitude range ' + get_var_string((*i).var1) + $
        ': [' + rm0es((*i).abs_min_1[isp]) + ',' + $
        rm0es((*i).abs_max_1[isp]) + ']'
      PRINT, 'Absolute amplitude range ' + get_var_string((*i).var2) + $
        ': [' + rm0es((*i).abs_min_2[isp]) + ',' + $
        rm0es((*i).abs_max_2[isp]) + ']'

      set_output, diag, dat=[REFORM(z_values_avg[*,*,isp])]
    ENDIF ELSE BEGIN
      PLOT, (*i).x_axis, z_values_avg[*,isp], $
        TITLE='!6density function of '+get_var_string((*i).var1,/fancy), $
        YTITLE='!6count', XTITLE=get_var_string((*i).var1,/fancy), $
        /XSTYLE, COLOR=1, PSYM=10

      PRINT, (*diag).name + ': species ' + spec[sp].name + $
        ', absolute amplitude range ' + get_var_string((*i).var1) + ': [' + $
        rm0es((*i).abs_min_1[isp]) + ',' + rm0es((*i).abs_max_1[isp]) + ']'

      set_output, diag, sp, header=[get_var_string((*i).var1),'count'], $
        dat=[[(*i).x_axis],[z_values_avg[*,isp]]]
    ENDELSE

    plot_info_str, diag

    set_output, diag, sp, /reset
  ENDFOR ; --- species loop

END
