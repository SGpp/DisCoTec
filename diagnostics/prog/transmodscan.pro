FUNCTION transmodscan_info

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN BEGIN
    n_vars = N_ELEMENTS((*scanlog).dimarr)+4
    IF n_vars GT 0 THEN BEGIN
      ext_vars = STRARR(3,n_vars)
      ext_vars[0,0] = 'eigenvalues'
      ext_vars[0,n_vars-1] = 'exponent'
      ext_vars[0,n_vars-2] = 'xmodes'
      ext_vars[0,n_vars-3] = 'const'
      ext_vars[0,n_vars-4] = 'refvalues'
      ext_vars[0,1:n_vars-5] = (*scanlog).dimarr[1:n_vars-5]
      ext_vars[1,*] = '0'
      ext_vars[2,0:n_vars-5] = 'type index between '+rm0es(0)+' and '+$
        rm0es((*scanlog).dsize[*]-1)+' (corresponds to values from '+$
	rm0es(STRING([1,(*scanlog).axxis[*,0]]))+' to '+$
        rm0es(STRING([(*scanlog).dsize[0]-1,(*scanlog).rd[*]]))+')'
      ext_vars[2,n_vars-1]='choose the weighting function exponent; default 2'
      ext_vars[2,n_vars-2]='choose number of accounted k_x modes (similar to nkx0); default: all'
      ext_vars[2,n_vars-3]='set constant to multiply transportvalues with, default: 1.0'
      ext_vars[2,n_vars-4]='transport normalization value taken from nrg file; default: 14 (indices start at 0 and run over all species)'
    ENDIF ELSE ext_vars = ''
  ENDIF ELSE ext_vars = ''  

  RETURN, {$
      type      : 'scan',$
      title     : 'Transport model',$
      help_text : ['Plots quasilinear transport model values'+$
      	    	   ' over chosen scanparams'],$
      ext_vars  : ext_vars}

END

;#####################################################################


PRO transmodscan_init,diag

  COMMON global_vars

  e = (*diag).external_vars
  IF NOT KEYWORD_SET(*e.exponent) THEN (*e.exponent)=2
  IF NOT KEYWORD_SET(*e.xmodes) THEN (*e.xmodes)=par.nx0
  IF NOT KEYWORD_SET(*e.const) THEN (*e.const)=1
  IF NOT KEYWORD_SET(*e.refvalues) THEN (*e.refvalues)=14
  IF (*e.xmodes) GT par.nx0 THEN (*e.xmodes)=par.nx0
  IF par.nx0 MOD 2 NE (*e.xmodes) MOD 2 THEN $
    IF (*e.xmodes) MOD 2 EQ 0 THEN (*e.xmodes)-=1 ELSE (*e.xmodes)+=1
  

  IF (*diag).table_entry[0,0] EQ '' THEN RETURN  
  
  n_vars = N_ELEMENTS((*diag).table_entry[3,*])
  IF n_vars GE 1 THEN eigenvalues = (*diag).table_entry[3,0]
  plotparA=REFORM((*diag).table_entry[3,1:n_vars-2])
  IF STRLEN(eigenvalues) EQ 0 THEN eigenvalues='*'

  initgrowthfield=MAKE_ARRAY((*scanlog).dsize,/DOUBLE) 
  eigengrowthfield=MAKE_ARRAY((*scanlog).dsize,/DOUBLE) 
  kyval = MAKE_ARRAY((*scanlog).dsize,/DOUBLE)
  kperp_sq = MAKE_ARRAY((*scanlog).dsize,/DOUBLE)
  ratiofield=REFORM(MAKE_ARRAY([N_ELEMENTS(*e.refvalues),(*scanlog).dsize],/DOUBLE) )

  maxsize=(*scanlog).dsize

  ind=1
  IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 THEN BEGIN
      maxsize(ind)=1
  ENDIF ELSE BEGIN 
  ENDELSE

  max_ky=MAKE_ARRAY(maxsize,/DOUBLE)
  ky_max=MAKE_ARRAY(maxsize,/DOUBLE)

  o_nky0 = par.nky0-par.ky0_ind ;original nky0
  
  i = set_internal_vars(diag,{$
      refvalues   : *e.refvalues,$
      xmodes      : *e.xmodes,$
      const      : *e.const,$
      eigenvalues  : eigenvalues,$  
      plotparA  : plotparA,$  
      ky_max    : ky_max,$
      max_ky    : max_ky,$
      o_nky0    : o_nky0,$
      initrunnum  : 1,$
      eigenrunnum  : 1,$
      gamofky : DBLARR(o_nky0),$
      ratiofield : ratiofield,$
      initgrowthfield : initgrowthfield,$
      eigengrowthfield : eigengrowthfield,$
      kperp_sq : kperp_sq,$
      kyval : kyval,$
      exponent : *e.exponent,$
      max_arr : DBLARR(o_nky0),$
      maxsize : maxsize})    


                                ; request mom/field variables
  series.request_nrg = 1
  fft_format, kxky=0 

                                ; read omega.dat
END

;######################################################################
FUNCTION transmodscan_calcweight, diag,ix,iz
  COMMON global_vars
  
  i = (*diag).internal_vars

  iy = par.ky0_ind

  phi_exp = (ABS((*mom[0,0].kxky)[ix,iy,iz]))^(*i).exponent*$
            (*series.geom).jac_norm[iz]
print, (*par.kx)[ix,iy]
  k_sq_weighed = (((*par.kx)[ix,iy])^2*$
              (*series.geom).gxx[iz]+(*par.ky)[iy]*(*par.kx)[ix,iy]*$
              2*(*series.geom).gxy[iz]+((*par.ky)[iy])^2*$
              (*series.geom).gyy[iz])*phi_exp

  return, [k_sq_weighed, phi_exp]
END

;######################################################################
PRO transmodscan_calc, diag, isort, dontcalc
  COMMON global_vars

  i = (*diag).internal_vars

  print, 'calc_transport'


;  OPENR,file_lun,gui.out.data_path+'omega_'+$
;        (STRSPLIT(series.run_labels,',',/EXTRACT))[0],/GET_LUN,ERROR=err 
  IF dontcalc EQ 0 THEN BEGIN
      k_perp_sq = DBLARR((*i).o_nky0)
      phi_sum = DBLARR((*i).o_nky0)
      IF (par.nx0 MOD 2) EQ 0 THEN BEGIN
          FOR ix = 0, ((*i).xmodes/2-1) DO BEGIN
              FOR iz = 0, par.nz0-1 DO BEGIN
                  retval=transmodscan_calcweight(diag, ix,iz)
                  k_perp_sq=k_perp_sq + retval[0]
                  phi_sum=phi_sum +retval[1]
              ENDFOR
          ENDFOR
          FOR ix = par.nx0-1, (par.nx0-(*i).xmodes/2)+1,-1 DO BEGIN
              FOR iz = 0, par.nz0-1 DO BEGIN
                  retval=transmodscan_calcweight(diag, ix,iz)
                  k_perp_sq=k_perp_sq + retval[0]
                  phi_sum=phi_sum +retval[1]
              ENDFOR
          ENDFOR
      ENDIF ELSE BEGIN
; ix == 1
          FOR ix = 0, FIX((*i).xmodes/2) DO BEGIN
              FOR iz = 0, par.nz0-1 DO BEGIN 
                  retval=transmodscan_calcweight(diag,ix,iz)
                  k_perp_sq=k_perp_sq + retval[0]
                  phi_sum=phi_sum +retval[1]
              ENDFOR
          ENDFOR
          FOR ix = par.nx0-1, par.nx0-FIX((*i).xmodes/2)+1,-1  DO BEGIN
              FOR iz = 0, par.nz0-1 DO BEGIN 
                  retval=transmodscan_calcweight(diag,ix,iz)
                  k_perp_sq=k_perp_sq + retval[0]
                  phi_sum=phi_sum +retval[1]
              ENDFOR
          ENDFOR
      ENDELSE
      k_perp_sq=k_perp_sq/phi_sum

      IF FIX(mom_time) GT 0  THEN BEGIN
         s_ind = (FIX(series.run_labels)-1)*(*i).o_nky0
         e_ind = FIX(series.run_labels)*(*i).o_nky0-1
         (*i).initgrowthfield[s_ind:e_ind]=$
            REAL_PART((*scanlog).eigenfield[s_ind:e_ind])/k_perp_sq
         (*i).kyval[s_ind:e_ind] = (*par.ky)[par.ky0_ind:*]
         (*i).kperp_sq[s_ind:e_ind] = k_perp_sq
      ENDIF ELSE BEGIN
         ind_ev = ((*i).eigenrunnum-1)*(*scanlog).dsize[0]+isort
         (*i).eigengrowthfield[ind_ev]=$
            REAL_PART((*scanlog).eigenfield[ind_ev])/k_perp_sq
         (*i).kperp_sq[ind_ev] = k_perp_sq
         (*i).kyval[s_ind:e_ind] = (*par.ky)[par.ky0_ind:*]
      ENDELSE

  ENDIF Else print, 'out of plotrange'  
;  FREE_LUN, file_lun
END



;######################################################################

PRO transmodscan_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  (*i).eigenrunnum=FIX(series.run_labels)

; calc value for plot
  dontcalc=1
  IF scan_plotnec(diag,[(*i).eigenvalues,(*i).plotparA],(*i).eigenrunnum)$
    EQ 1 THEN BEGIN
;      fft_format, kxky=0 
      dontcalc=0
  ENDIF

; calc nrgindex for run and set ratiofield

  IF mom_time GT 0 THEN BEGIN
      nrgindex = N_ELEMENTS((*nrg).time) - 1
      FOR ind=0,N_ELEMENTS((*i).refvalues)-1 DO BEGIN
          (*i).ratiofield[((*i).eigenrunnum-1)*N_ELEMENTS((*i).refvalues)+ind]=$
           (*nrg)[nrgindex].data[(*i).refvalues[ind]]
      ENDFOR
  ENDIF ELSE BEGIN
      nrgindex = WHERE((*nrg).time EQ mom_time)  
      FOR isort = 0, (*scanlog).dsize[0]-1 DO BEGIN
          IF ((*scanlog).sortfield[((*i).eigenrunnum-1)*$
            (*scanlog).dsize[0]+isort]) EQ FIX(mom_time) THEN BEGIN
              FOR ind=0,N_ELEMENTS((*i).refvalues)-1 DO BEGIN
                  (*i).ratiofield[(((*i).eigenrunnum-1)*(*scanlog).dsize[0]+isort)*N_ELEMENTS((*i).refvalues)+ind]=$
                   (*nrg)[nrgindex].data[(*i).refvalues[ind]]
              ENDFOR
          ENDIF
      ENDFOR    
  ENDELSE
  
; initvalues
  IF mom_time GT 0 THEN BEGIN
      transmodscan_calc,diag,0,dontcalc
;     if scan over ky --> max after scanrange
      IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 THEN BEGIN
          IF (FIX(series.run_labels) MOD ((*scanlog).dsize[1]/(*i).o_nky0)$
              EQ 0) AND (FIX(series.run_labels) GT 0) THEN BEGIN
              (*i).max_ky[(*i).initrunnum-1]=MAX((*i).$
              initgrowthfield[((*i).initrunnum-1)*(*i).o_nky0*$
                (*scanlog).dsize[1]:(*i).initrunnum*(*i).o_nky0*$
                (*scanlog).dsize[1]-1],Max_Subscript,/NAN)
              (*i).ky_max[(*i).initrunnum-1]=(*scanlog).rd[0]-$
                ((*scanlog).dsize[1]-Max_Subscript-1)*(*scanlog).deltad[0]
              (*i).initrunnum+=1
          ENDIF
;     max of ky if nky > 1
      ENDIF ELSE BEGIN
          (*i).max_ky[(*i).initrunnum-1]=$
            MAX((*i).initgrowthfield[(*i).initrunnum*(*i).o_nky0-1:$
            (*i).initrunnum*(*i).o_nky0+(*i).o_nky0-2],Max_Subscript,/NAN)
          (*i).ky_max[(*i).initrunnum-1]=(*scanlog).rd[0]-$
            ((*scanlog).dsize[1]-Max_Subscript-1)*(*scanlog).deltad[0]
          (*i).initrunnum+=1    
      ENDELSE
; eigenvalues
  ENDIF ELSE BEGIN
      found=0
;     calc sorted values     
      FOR isort = 0, (*scanlog).dsize[0]-1 DO BEGIN
          IF ((*scanlog).sortfield[((*i).eigenrunnum-1)*$
             (*scanlog).dsize[0]+isort]) EQ FIX(mom_time) THEN BEGIN
              transmodscan_calc,diag,isort,dontcalc
              found+=1
          ENDIF
      ENDFOR
      IF found NE 1 THEN print,'not all values found'
;     if las value of 1 run      
      IF FIX(mom_time) EQ -1 THEN BEGIN
;         if over ky is scanned
          IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 THEN BEGIN
;             if end of ky-direction (and not the first run)        
              IF (*i).eigenrunnum MOD (*scanlog).dsize[1] EQ 0 $
                AND (*i).eigenrunnum GT 1 THEN BEGIN
;                 max of ky-range
                  FOR iev=0, (*scanlog).dsize[0]-1 DO BEGIN
                      (*i).max_ky[((*i).eigenrunnum/(*scanlog).$
                        dsize[1]-1)*(*scanlog).dsize[0]+iev]=$
                        MAX((*i).eigengrowthfield[(((*i).eigenrunnum-$
                        (*scanlog).dsize[1])*(*scanlog).dsize[0]+iev):$
                        (*i).eigenrunnum*(*scanlog).dsize[0]-1:$
                        ((*scanlog).dsize[0])],Max_Subscript,/NAN)
                      (*i).ky_max[(*i).eigenrunnum/(*scanlog).$
                        dsize[1]-1]=(*scanlog).rd[0]-$
                        ((*scanlog).dsize[1]-Max_Subscript-1)*$
                        (*scanlog).deltad[0]
                  ENDFOR
              ENDIF
          ENDIF ELSE BEGIN
              FOR iev=0, (*scanlog).dsize[0]-1 DO BEGIN
                  (*i).max_ky[((*i).eigenrunnum-1)*(*scanlog).dsize[0]+iev]=$
                    (*i).eigengrowthfield[((*i).eigenrunnum-1)*(*scanlog).dsize[0]+iev]
              ENDFOR
          ENDELSE
      ENDIF
  ENDELSE
END



;######################################################################

PRO transmodscan_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  
;  print, 'transmodscan_output'
  IF par.ntimesteps GT 0 THEN BEGIN
      growthfield=(*i).initgrowthfield
  ENDIF ELSE BEGIN
      growthfield=(*i).eigengrowthfield
  ENDELSE

  set_output, diag, header = ['ky','kperp_sq_av', 'growthfield'],$
              dat=[[TRANSPOSE((*i).kyval)],[TRANSPOSE((*i).kperp_sq)],$
                    [TRANSPOSE(growthfield)]]


;*************calculate quasilinear ratios***************************
  FOR ind=N_ELEMENTS((*i).refvalues)-1,0,-1 DO BEGIN
      (*i).ratiofield[ind:*:N_ELEMENTS((*i).refvalues)]=$
         (*i).ratiofield[ind:*:N_ELEMENTS((*i).refvalues)]/(*i).ratiofield[0:*:N_ELEMENTS((*i).refvalues)]
  ENDFOR
  (*i).ratiofield=(*i).ratiofield*(*i).const

;************************plot all ratio fields***********************
  FOR ind=0,N_ELEMENTS((*i).refvalues)-1 DO BEGIN
      scan_plot, diag, 'quasilinear transport',growthfield*$
       (*i).ratiofield[ind:*:N_ELEMENTS((*i).refvalues)],'abs', 0, 1, 2, 0.01, 0.06, 0.99, 0.16, 2,$
       suffix='-rval_'+rm0es(STRING(ind))
      SPAWN,"echo "+"'"+'set ylabel "nrg '+rm0es(STRING((*i).refvalues[ind]))+$
            '"; splot '+'"'+gui.out.out_path+'transmodscan-rval_'+rm0es(STRING(ind))+$
            '_0001.dat"'+" u 1:2:3 w lp' | gnuplot -persist" 
  ENDFOR
  
;  print, (*i).ratiofield
;************************plot maximazion over ky********************* 

  IF NOT ((*scanlog).dsize[1] EQ (*i).maxsize[1]) THEN BEGIN
      IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 $
        THEN plotparA1='0'     
      dummy=(*scanlog).dsize
      (*scanlog).dsize=(*i).maxsize
      
      maxplotfield=(*i).max_ky
;printjob for all given refvalues
      FOR ind=0,N_ELEMENTS((*i).refvalues)-1 DO BEGIN
          FOR subs=0,N_ELEMENTS((*i).max_ky)-1 DO BEGIN
              adummy=(*i).max_ky[subs]*$
                     (*i).ratiofield[[ind+N_ELEMENTS((*i).refvalues)*WHERE(growthfield EQ (*i).max_ky[subs])]]
              maxplotfield[subs]=adummy[0]
          ENDFOR
          scan_plot, diag, 'max of quasilinear transport',maxplotfield,'abs',$
                     0, 1, 3, 0.01, 0.06, 0.99, 0.16, 2, suffix='m-rval_'+rm0es(STRING(ind))
          SPAWN,"echo "+"'"+'set ylabel "nrg '+rm0es(STRING((*i).refvalues[ind]))+$
                '"; plot '+'"'+gui.out.out_path+'transmodscanm-rval_'+$
                rm0es(STRING(ind))+'_0001.dat"'+" w lp' | gnuplot -persist" 
      ENDFOR
      scan_plot, diag, 'ky_max',(*i).ky_max,'abs',$
                 0, 1, 3, 0.01, 0.06, 0.99, 0.16, 2, suffix='kymax'
      (*scanlog).dsize=dummy
  ENDIF
  


END
