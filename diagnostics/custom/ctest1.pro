PRO ctest1_init, diag, var=var, sw=sw

  COMMON global_vars

  IF NOT KEYWORD_SET(var) THEN var = 0
  IF NOT KEYWORD_SET(sw) THEN sw = 0

  i = set_internal_vars(diag,{$
    var    : var,$
    iota   : 0L,$
    rms : 0.0})

  fft_format, kxky=0

END

;######################################################################

PRO ctest1_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  (*i).iota = (*i).iota + 1

  rms = TOTAL(TOTAL(ABS((*mom[0,0].kxky)[*,*,*])^2,3),1)
  rms = SQRT(2.0 * TOTAL(rms[1:*]) + rms[0])

  (*i).rms = rms

END

;######################################################################

PRO ctest1_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  set_output, diag, 0, /ps

  PRINT, '-----'
  PRINT, rm0es((*i).iota), ' steps handled by ctest1_loop'
  PRINT, 'compare series.step_count: ', rm0es(series.step_count)
  PRINT, 'rms: ', rm0es((*i).rms)
  PRINT, '-----'
  
  xarr = [0,1]
  yarr = [-1,1]
  PLOT, xarr, yarr, COLOR=1

  set_output, diag, 0, dat=[xarr,yarr]
  set_output, diag, 0, /reset

END

PRO ctest1 & END



