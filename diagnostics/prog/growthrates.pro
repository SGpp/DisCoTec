FUNCTION growthrates_info

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN BEGIN
    n_vars = N_ELEMENTS((*scanlog).dimarr)
    IF n_vars GT 0 THEN BEGIN
      ext_vars = STRARR(3,n_vars)
      ext_vars[0,0] = 'eigenvalues'
      ext_vars[0,1:n_vars-1] = (*scanlog).dimarr[1:n_vars-1]
      ext_vars[1,*] = '0'
      ext_vars[2,*] = 'type index between '+rm0es(0)+' and '+$
        rm0es(STRING((*scanlog).dsize[*]-1))+' meaning between '+$
	rm0es(STRING([0,(*scanlog).axxis[*,0]]))+' and '+rm0es(STRING([(*scanlog).dsize[0]-1,(*scanlog).rd[*]]))
    ENDIF ELSE ext_vars = ''
  ENDIF ELSE ext_vars = ''
  
  RETURN, {$
    type      : 'scan',$
    title     : 'Growth rates / Freq.',$
    help_text : ['plots out the growth rates over chosen scanparams'],$
    ext_vars  : ext_vars}
          
END

;######################################################################

PRO growthrates_init, diag

  COMMON global_vars
  
  IF (*diag).table_entry[0,0] EQ '' THEN BEGIN
    PRINT, 'skipping growthrates diag'
    (*diag).selected = 0
    RETURN
  ENDIF

  n_vars = N_ELEMENTS((*diag).table_entry[3,*])
  IF n_vars GE 1 THEN eigenvalues = (*diag).table_entry[3,0]
  plotparA = (*diag).table_entry[3,1:N_ELEMENTS((*diag).table_entry[3,*])-1]
  reimabs = ''

  i = set_internal_vars(diag,{$
    eigenvalues : eigenvalues,$  
    plotparA    : plotparA,$  
    reimabs     : reimabs})

END

;######################################################################

PRO growthrates_loop, diag

END

;######################################################################

PRO growthrates_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  PRINT, 'growthrates_output'

  scan_plot, diag, ' modestructure', (*scanlog).eigenfield, $
    (*i).reimabs, 0, 1, 3, 0.01, 0.06, 0.99, 0.16, 2

END











