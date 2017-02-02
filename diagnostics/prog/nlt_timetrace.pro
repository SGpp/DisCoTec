FUNCTION nlt_timetrace_info

  RETURN, {$
    type      : 'nlt',$
    title     : 'NLT time traces',$
    help_text : ['Plots time traces for subsets of nonlinear transfer functions '],$
    ext_vars  : [['kx_data','0','  1000: sum over kx.  '+$
                  'Otherwise: slice at kx index (negative indices allowed).'+$
                  '  (Default: sum).'],$
                 ['ky_data','0','  1000: sum over ky.  '+$
                  'Otherwise: slice at ky index (negative indices allowed).'+$
                  '  (Default: sum).']$
                                                                                              ]}





END

;######################################################################

PRO nlt_timetrace_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.kx_data) EQ 0 THEN (*e.kx_data)=1000
  IF N_ELEMENTS(*e.ky_data) EQ 0 THEN (*e.ky_data)=1000

  IF par.lx EQ 0 THEN BEGIN
    STOP, 'lx cannot be 0'
  ENDIF

  i = set_internal_vars(diag,{$
    kx_data   : *e.kx_data,$
    ky_data   : *e.ky_data,$
    yind      : FINDGEN(2*par.nky0-1),$
    xind      : FINDGEN(par.nkx0),$
    time_id   : PTR_NEW(),$
    data_id   : PTR_NEW()})

    (*i).xind = ((*i).xind - par.nkx0 + par.nkx0 / 2 + 1) * $
      (2.0 * !PI / series.lx)
    (*i).yind = ((*i).yind - par.nky0 + 1) * (2.0 * !PI / series.ly)

END

;######################################################################

PRO nlt_timetrace_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nkx0 = par.nkx0
  nky0 = par.nky0

  data=*nlt2d

  (*i).time_id = store_step((*i).time_id,nlt_time)
  (*i).data_id = store_step((*i).data_id,data)

END

;######################################################################

PRO nlt_timetrace_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  out_suffix=''
  IF par.nlp_gdt THEN out_suffix='_POD'

   IF par.nlp_gdt THEN BEGIN
     num_nlt_modes_out=par.num_nlt_pod_modes+1
   ENDIF ELSE BEGIN
     num_nlt_modes_out=par.num_nlt_modes
   ENDELSE

   data = store_step((*i).data_id,/get)
   time_id = store_step((*i).time_id,/get)
   xind=(*i).xind
   yind=(*i).yind
   en_var_string=STRARR(num_nlt_modes_out)
   comments=STRARR(num_nlt_modes_out)
   kxmin=2.0*!PI/par.lx
   
   data_kx=FLTARR(series.step_count-1,par.nkx0,num_nlt_modes_out)
   time_trace=FLTARR(series.step_count-1,num_nlt_modes_out)
   time=FLTARR(series.step_count-1)

   FOR n=0,series.step_count-2 DO BEGIN
    ;print,n,*time_id[n],series.step_count-2
    time[n]=*time_id[n] 
   ENDFOR

   IF (*i).ky_data EQ 1000 THEN BEGIN
     FOR var=0,num_nlt_modes_out-1 DO BEGIN
      FOR n=0,series.step_count-2 DO BEGIN
        data_kx[n,*,var]=TOTAL((*data[n])[*,*,var],2)
      ENDFOR
     ENDFOR
   ENDIF ELSE BEGIN
    ky_index=(*i).ky_data
    ky_index=ky_index+(2*par.nky0-1)/2
    FOR var=0,num_nlt_modes_out-1 DO BEGIN
     FOR n=0,series.step_count-2 DO BEGIN
       data_kx[n,*,var]=(*data[n])[*,ky_index,var]
     ENDFOR
    ENDFOR
   ENDELSE  

   IF (*i).kx_data EQ 1000 THEN BEGIN
     FOR var=0,num_nlt_modes_out-1 DO BEGIN
       time_trace[*,var]=TOTAL(data_kx[*,*,var],2) 
     ENDFOR
   ENDIF ELSE BEGIN
     kx_index=(*i).kx_data
     kx_index=kx_index+(par.nkx0/2-1)
     FOR var=0,num_nlt_modes_out-1 DO BEGIN
       time_trace[*,var]=data_kx[*,kx_index,var]
     ENDFOR
   ENDELSE  

   ;set titles

   kxmin=2.0*!PI/par.lx
   kxstring=STRARR(num_nlt_modes_out)
   kxind0=FLTARR(num_nlt_modes_out)
   FOR var=0,num_nlt_modes_out-1 DO BEGIN

    IF par.nlp_gdt THEN BEGIN
      mode_num=STRTRIM(var+1,2)
      kxind0[var]=kxmin*par.nlp_kxind
      IF par.nlp_kxind GT par.nx0/2 THEN BEGIN
       kxind0[var]=-1.0*(par.nx0-par.nlp_kxind)*kxmin
      ENDIF
      kxstring[var]=STRING(kxind0[var],FORMAT='(F5.2)')
      IF var EQ num_nlt_modes_out-1 THEN mode_num='res.'
      en_var_string[var]='T!DNL!N(k!Dx!N='+kxstring[var]$
         +',k!Dy!N='+STRING(par.kymin*par.nlp_kyind,FORMAT='(F5.2)')+$
         ')'+'(POD n='+mode_num+')'
      comments[var]='T_NL(k_x='+kxstring[var]$
         +',k_y='+STRING(par.kymin*par.nlp_kyind,FORMAT='(F5.2)')+$
         ')'+'(POD n='+mode_num+')'
    ENDIF ELSE BEGIN
      kxind0[var]=kxmin*(*par.kx_nlt_ind)[var]
      IF (*par.kx_nlt_ind)[var] GT par.nx0/2 THEN BEGIN
        kxind0[var]=-1.0*(par.nx0-(*par.kx_nlt_ind)[var])*kxmin
      ENDIF
      kxstring[var]=STRING(kxind0[var],FORMAT='(F5.2)')
      en_var_string[var]='T!DNL!N(k!Dx!N='+kxstring[var]$
         +',k!Dy!N='+STRING(par.kymin*(*par.ky_nlt_ind)[var],FORMAT='(F5.2)')+')'
      comments[var]='T_NL(k_x='+kxstring[var]$
         +',k_y='+STRING(par.kymin*(*par.ky_nlt_ind)[var],FORMAT='(F5.2)')+')'
    ENDELSE

     IF (*i).kx_data EQ 1000 THEN BEGIN
       IF (*i).ky_data EQ 1000 THEN BEGIN
         en_var_string[var]='<'+en_var_string[var]+">!Dkx',ky'!N"       
         comments[var]='<'+comments[var]+">_kx',ky'"       
       ENDIF ELSE BEGIN ;ky_data
         en_var_string[var]='<'+en_var_string[var]+">!Dkx'!N[k'!Dy!N="+$
              STRING(par.kymin*(*i).ky_data,'(F5.2)')+']'       
         comments[var]='<'+comments[var]+">_kx'[k'y="+$
              STRING(par.kymin*(*i).ky_data,'(F5.2)')+']'       
       ENDELSE ;ky_data
     ENDIF ELSE BEGIN ;kx_data
       IF (*i).ky_data EQ 1000 THEN BEGIN
         en_var_string[var]='<'+en_var_string[var]+">_ky'[k'_x="+$
              STRING(kxmin*(*i).kx_data,'(F5.2)')+']'       
         comments[var]='<'+comments[var]+">_ky'[k'_x="+$
              STRING(kxmin*(*i).kx_data,'(F5.2)')+']'       
       ENDIF ELSE BEGIN ;ky_data
         en_var_string[var]=en_var_string[var]+"[k'_x="+$
              STRING(kxmin*(*i).kx_data,'(F5.2)')+$
              ",k'_y="+STRING(par.kymin*(*i).ky_data,'(F5.2)')+']'   
         comments[var]=comments[var]+"[k'_x="+$
              STRING(kxmin*(*i).kx_data,'(F5.2)')+$
              ",k'_y="+STRING(par.kymin*(*i).ky_data,'(F5.2)')+']'   
       ENDELSE ;ky_data
     ENDELSE ;kx_data

   ENDFOR
  
   csize = 1.41 
   cthick = 3.0
   set_output, diag, /ps  ,suffix=out_suffix
  ;Output plot file here.

    FOR var = 0, num_nlt_modes_out - 1 DO BEGIN    ; loop over all variables to plot
      yrange=[MIN(time_trace[*,var]),MAX(time_trace[*,var])]
      PLOT, time, time_trace[*,var], COLOR=1, XTITLE=get_var_string(1,/time,/fancy,/unit), $
           POSITION=[0.15,0.4,0.9,0.9], $
            YRANGE=yrange, $
          TITLE=en_var_string[var] 
    ENDFOR
  
  
    FOR var = 0, num_nlt_modes_out - 1 DO BEGIN     ; loop over all variables to plot
      set_output, diag,$
          dat=[[time],[time_trace[*,var]]],append=(var GT 0),suffix=out_suffix,$
          commentlines=comments[var]
    ENDFOR
  
      set_output, diag, /reset, eps=(series.step_count EQ 1),suffix=out_suffix
  
    PTR_FREE, time_id, data

    END



