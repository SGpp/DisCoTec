FUNCTION showscan_info

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN BEGIN
    n_vars = N_ELEMENTS((*scanlog).dimarr)
    IF n_vars GT 0 THEN BEGIN
      ext_vars = STRARR(3,n_vars)
      ext_vars[0,0] = 'showdate'
      ext_vars[0,1:n_vars-1] = (*scanlog).dimarr[1:n_vars-1]
      ext_vars[1,*] = '0'
      ext_vars[2,1:n_vars-1] = 'type index between '+rm0es(0)+' and '+$
        STRING((*scanlog).dsize[1:n_vars-1]-1)+' meaning between '+$
	rm0es(STRING((*scanlog).axxis[*,0]))+' and '+rm0es(STRING((*scanlog).rd[*]))
    ENDIF ELSE ext_vars = ''
  ENDIF ELSE ext_vars = ''  

  RETURN, {$
          type      : 'scan',$
          title     : 'show data scan diagnostic',$
          help_text : ['Plots the chosen date'+$
                       ' over the specified parameters. '+$
                       'The dates from 0 to 3 are: '+$
                       '0 : scalarproduct of the two vectors, '+$
                       '1 : time for eigenvalue calculation, '+$
                       '2 : iterations for eigenvalue calculation, '+$
                       '3 : lx, '+$
                       '4 : ly, '+$                       
                       '5 : dt_max. '+$
                       '6 : total_time'+$
                       '7 : Number of iterations'],$
          ext_vars  : ext_vars}
  
END

;#####################################################################

PRO showscan_init,diag

COMMON global_vars

  IF (*diag).table_entry[0,0] EQ '' THEN BEGIN
      PRINT, 'skipping showscan diag'
      (*diag).selected = 0
      RETURN
  ENDIF
  
  n_vars = N_ELEMENTS((*diag).table_entry[3,*])
  date2plot = 0
  IF n_vars GE 1 THEN date2plot = FIX((*diag).table_entry[3,0])
  plotparA = (*diag).table_entry[3,1:n_vars-1]
  
  eigenvalues = '*'
  runnum = 0
  
  CASE date2plot OF
      0 : TITLE='scalarproduct of the two vectors '
      1 : TITLE='time for eigenvalue calculation'
      2 : TITLE='iterations for eigenvalue calculation'
      3 : TITLE='lx'
      4 : TITLE='ly'
      5 : TITLE='dt_max'
      6 : TITLE='total_time'
      7 : TITLE='Number of iterations'
      ELSE : TITLE = '???'
  ENDCASE
  
  i = set_internal_vars(diag,{$
      title : title,$
      date2plot : date2plot,$
      eigenvalues  : eigenvalues,$  
      runnum : runnum,$
      plotparA  : plotparA})
;fft_format, kxky=vars
END

;######################################################################

PRO showscan_loop, diag

END

;######################################################################

PRO showscan_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

;first call of finddate for datafield generation
  abc = (*scanlog).dsize
  run_label = STRING(1,FORMAT='(I04)')
  plotn=N_ELEMENTS(finddate(diag,run_label,(*i).date2plot,(*i).title))
  abc[0]=plotn
  datafield = MAKE_ARRAY(abc,/DOUBLE) 
  FOR ix = 0, N_ELEMENTS(datafield)/plotn - 1 DO BEGIN
      run_label = STRING(ix+1,FORMAT='(I04)')  
      print,finddate(diag,run_label,(*i).date2plot,(*i).title) 
      datafield[ix*plotn:(ix+1)*plotn-1] = finddate(diag,run_label,(*i).date2plot,(*i).title)
  ENDFOR

  PRINT, 'showscan_output'
  PRINT, datafield[*]

  dsizesave = (*scanlog).dsize
  (*scanlog).dsize[0] =   abc[0]
  scan_plot, diag, (*i).title, datafield, $
    'non', 0, 1, 2, 0.01, 0.06, 0.99, 0.16, 2
  (*scanlog).dsize = dsizesave

END
