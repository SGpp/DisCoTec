FUNCTION nrgscan_info

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN BEGIN
    n_vars = N_ELEMENTS((*scanlog).dimarr) + 2
    IF n_vars GT 0 THEN BEGIN
      ext_vars = STRARR(3,n_vars)
      ext_vars[0,0] = 'mode'
      ext_vars[0,1:n_vars-3] = (*scanlog).dimarr[1:n_vars-3]
      ext_vars[0,n_vars-1] = 'column'
      ext_vars[0,n_vars-2] = 'norm column'
      ext_vars[1,*] = '0'
      ext_vars[2,0:n_vars-3] = 'type index between ' + rm0es(0) + $
        ' and ' + rm0es((*scanlog).dsize[*]-1)+' meaning between ' + $
	rm0es([0,(*scanlog).axxis[*,0]]) + ' and ' + $
        rm0es([(*scanlog).dsize[0]-1,(*scanlog).rd[*]])
      ext_vars[2,n_vars-1] = 'type columnnumber of nrg value to plot.' + $
        ' (DEFAULT=4)'
      ext_vars[2,n_vars-2] = 'type columnnumber normalization nrg value' + $
        ' with. (DEFAULT=0)'
    ENDIF ELSE ext_vars = ''
  ENDIF ELSE ext_vars = ''

  RETURN, {$
    type      : 'scan',$
    title     : 'nrg scan diagnostic',$
    help_text : ['Plots the structure of the chosen column of the nrg'+$
                 ' over the specified parameters. '+$
                 'The columns from 0 to 7 are: '+$
                 get_nrg_string(0)+', '+get_nrg_string(1)+', '+$
                 get_nrg_string(2)+', '+get_nrg_string(3)+', '+$
                 get_nrg_string(4)+', '+get_nrg_string(5)+', '+$
                 get_nrg_string(6)+', '+get_nrg_string(7)+$
                 '     there is also an 8th column for the field phi^2'+$
                 '     and a 9th for phi*n'],$
    ext_vars  : ext_vars}
  
END

;#####################################################################

PRO nrgscan_init, diag

  COMMON global_vars
  
  IF (*diag).table_entry[0,0] EQ '' THEN BEGIN
    PRINT, 'skipping nrgscan diag'
    (*diag).selected = 0
    RETURN
  ENDIF

  n_vars = N_ELEMENTS((*diag).table_entry[3,*])
  IF n_vars GE 1 THEN eigenvalues = (*diag).table_entry[3,0]
  eigenvalues='*'
  plotparA=REFORM((*diag).table_entry[3,1:n_vars-3])
  nrgpar=(*diag).table_entry[3,n_vars-1]

  IF STRLEN(nrgpar) EQ 0 THEN nrgpar = 4
  nrgnormpar = (*diag).table_entry[3,n_vars-2]
  IF STRLEN(nrgnormpar) EQ 0 THEN nrgnormpar = 0

  nrgfield = MAKE_ARRAY((*scanlog).dsize,/DOUBLE,/NOZERO)

  i = set_internal_vars(diag,{$
    eigenvalues : eigenvalues,$
    plotparA    : plotparA,$
    nrgpar      : FIX(nrgpar) MOD 8,$
    nrgnormpar  : FIX(nrgnormpar) MOD 10,$
    init        : REPLICATE({nrgfield : nrgfield},gui.out.n_spec_sel),$
    eigen 	  : REPLICATE({nrgfield : nrgfield},gui.out.n_spec_sel)})

  fft_format, kxky = 0
  series.request_nrg = 1

END



;######################################################################

PRO nrgscan_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars
  runnr=FIX(series.run_labels[0])

  ; calculate <|phi|^2> if necessary
  IF (*i).nrgnormpar GT 7 AND $
    scan_plotnec(diag,[(*i).eigenvalues,(*i).plotparA],runnr)$
    EQ 1 THEN BEGIN
      phi_sq_av = TOTAL(TOTAL(TOTAL(ABS((*mom[0,0].kxky)$
      [*,(par.ky0_ind EQ 0):*,*])^2,1),1)*(*series.geom).jac_norm)/par.nz0
      IF par.ky0_ind EQ 0 THEN $
        phi_sq_av += TOTAL(TOTAL(TOTAL(ABS((*mom[isp,(*i).var].kxky)$
	  [*,0,*])^2,1),1)*(*series.geom).jac_norm)/par.nz0
  ENDIF
  
  ; set nrgindex
  IF mom_time GT 0 THEN nrgindex = N_ELEMENTS((*nrg).time) - 1 $
  ELSE nrgindex = WHERE((*nrg).time EQ mom_time)  

  ; loop over selected species
  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN 
    sp=(*gui.out.spec_select)[isp] ; sp is the "real" species index

    CASE (*i).nrgnormpar OF
      8 : divisor = phi_sq_av
      9 : divisor = sqrt((*nrg)[nrgindex].data[0,sp])*sqrt(phi_sq_av)
      ELSE : divisor=(*nrg)[nrgindex].data[(*i).nrgnormpar MOD 8,sp]
    ENDCASE      
    
    IF mom_time GT 0 THEN BEGIN
      (*i).init[isp].nrgfield[(runnr-1)*par.nky0:runnr*par.nky0-1] $
        = (*nrg)[nrgindex].data[(*i).nrgpar,sp]/divisor
    ENDIF ELSE BEGIN
      found=0
      FOR isort = 0, (*scanlog).dsize[0]-1 DO BEGIN
          IF ((*scanlog).sortfield[(runnr-1)*(*scanlog).dsize[0]+isort]) $
	  EQ FIX(mom_time) THEN BEGIN
            (*i).eigen[isp].nrgfield[(runnr-1)*(*scanlog).dsize[0]+isort] $
                = (*nrg)[nrgindex].data[(*i).nrgpar,sp]/divisor
             found+=1
          ENDIF
      ENDFOR
      IF found NE 1 THEN print,'not all values found'
    ENDELSE    
    
  ENDFOR ; species loop

END

;######################################################################

PRO nrgscan_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
;  print, 'nrgscan_output'
  
  nrg_field_names = [get_nrg_string(INDGEN(8),/fancy),$
    '!12<!9!!!7U!9!!!6!E2!N!12>!6',get_nrg_string(0,/fancy)+$
    '!E1/2!N!12<!9!!!7U!9!!!6!E2!N!12>!6!E1/2!N']
    
  n_names = N_ELEMENTS(nrg_field_names)
    
  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp=(*gui.out.spec_select)[isp] ; sp is the "real" species index
    
    ; some abbreviations to avoid long titles
    CASE spec[sp].name OF
      'electrons' : spec_name = 'e'
      'ions'	  : spec_name = 'i'
      ELSE  	  : spec_name = spec[sp].name
    ENDCASE
    
    title = '('+nrg_field_names[(*i).nrgpar]+'!I'+spec_name+ $
      '!N) / ('+nrg_field_names[(*i).nrgnormpar]+'!I'+spec_name+'!N)'
    
    IF par.ntimesteps GT 0 THEN BEGIN
      scan_plot, diag,title,(*i).init[isp].nrgfield,$
        'non', 0, 1, 2, 0.01, 0.06, 0.99, 0.16, 2, sp=sp
    ENDIF ELSE BEGIN
      print, (*i).eigen[isp].nrgfield
      scan_plot, diag,title,(*i).eigen[isp].nrgfield,$
        'non', 0, 1, 2, 0.01, 0.06, 0.99, 0.16, 2, sp = sp
    ENDELSE
  ENDFOR ; species loop
  
END
