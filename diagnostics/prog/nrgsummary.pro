FUNCTION nrgsummary_info

  RETURN, {$
    type      : 'nrg',$
    title     : 'summary',$
    help_text : ['Plots all nrg data, lists growth rates and'+$
                 'saturation levels'],$
    ext_vars  : [['tcorr','0','correlation time for error '+$
                 'computation; default: -1 (automatic calculation) '+$
                 'for nonlinear, 0 (no errors) for linear runs; '+$
                 'use -2 for simple standard deviation; set '+$
                 'to nonzero positive value for manual setting'],$
                 ['corr_field','0','nrg field/species for which '+$
                 'to compute tcorr (if tcorr=-1); default: '+$
                 '[[0,6],[0,0]] (n^2, Q_es for first species); '+$
                 'set to -1 for all fields/species'],$
                 ['show_SI','1','include table with SI values']]}

END

;######################################################################

PRO nrgsummary_init, diag

END

;######################################################################

PRO nrgsummary_loop, diag

END

;######################################################################

PRO nrgsummary_output, diag

  COMMON global_vars

  e = (*diag).external_vars
  tcorr = N_ELEMENTS(*e.tcorr) NE 1 ? -3 : *e.tcorr
  corr_field = N_ELEMENTS(*e.corr_field) GE 1 ? *e.corr_field : $
    [[0,6],[0,0]]
  IF (N_ELEMENTS(corr_field) EQ 1) AND (corr_field[0] GE 0) THEN $
    corr_field = [[corr_field],[0]]
  show_SI = KEYWORD_SET(*e.show_SI)
  IF (tcorr LT 0) AND (tcorr NE -1) AND (tcorr NE -2) THEN $
    tcorr = par.nonlinear ? -1 : 0

  plot_nrg_data, diag=diag, tcorr=tcorr, corr_field=corr_field, $
    show_SI=show_SI
  
END
