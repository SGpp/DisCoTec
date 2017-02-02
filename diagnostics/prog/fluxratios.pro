; Sep 9 2013 Ken Liao changes.
;   changed some help text
;   removed hardcoded limit on 70 entries in scan
;   changed phi to phi^2 in numerator and denominator
;   fix calculation of ratio_local2
;   fix ytitle and color for 2-dimensional surface plots
;   fix species label in comment line in ascii output
;   change ascii append file behavior
;   It is wrong to try to average over last 10 lines of nrg, since the nrg file is organized such that each time point takes several lines (depending on number of species). I have changed it to take last data point, or the closest point to mom_time if need_phi = true. This is only fixed for initial value problem solver.
FUNCTION fluxratios_info

  COMMON global_vars

  RETURN, {$
    type      : 'scan',$
    title     : 'Flux ratio scan',$
    help_text : ['Plots the ratios between the selected variables. '+$
                 'One ratio is calculated per pair. '+$
                 'Variables 0 to 7 correspond to the nrg variables, '+$
                 '8 to D_es, 9 to D_em, 10 to chi_es, 11 to chi_em, '+$
                 '12 to phi^2, and 13 to phi*n'],$
    ext_vars  : [['vars','0','array containing variable pairs, '+$
                  'e.g.: [4,0] gives <Gamma_es> / <|n|^2> '+$
                  '(default: [8,10])'],$
                 ['species','0','array containing pairs of '+$
                  'species for the variable pairs, where 1 is the first species and 0 is replaced by the selected species (default: 0,0)']]}

END

;######################################################################

 PRO fluxratios_init, diag

  COMMON global_vars

  IF NOT PTR_VALID(scanlog) THEN BEGIN
    PRINT, 'skipping ' + (*diag).name + ': no scanfile found'
    RETURN
  ENDIF

  ; set default values for vars and species
  e = (*diag).external_vars
  IF N_ELEMENTS(*e.vars) LT 1 THEN *e.vars = [8,10]
  IF N_ELEMENTS(*e.species) LT 1 THEN *e.species = [0,0]
  ; [0,0]: D/chi will be calculated for all species selected in the gui

  ; check if eigenvalue solver was used
  ; and if yes, how many eigenvalues are requested to be found
  ev_solv = par.n_ev GT 0 ; zero for initial value solver, 1 for eigenvalue solver
  n_ev = par.n_ev ; numbers of requested eigenvalues, zero for initial value solver

  variables = FIX(*e.vars)
  var_pairs = INTARR(N_ELEMENTS(variables)/2,2) ; for odd numbers, skip last

  ; exception handling for false input:
  w1 = WHERE(variables GT 13)
  w2 = WHERE(variables LT 0)
  IF (w1[0] NE -1) OR (w2[0] NE -1) THEN BEGIN
    printerror, 'Incorrect input in vars'
    (*diag).selected = 0
    RETURN
  ENDIF

  FOR v = 0, N_ELEMENTS(variables) / 2 - 1 DO var_pairs[v,*] = variables[2*v:2*v+1]

  species = INTARR(N_ELEMENTS(var_pairs))
  FOR v = 0, N_ELEMENTS(species) - 1 DO $
    IF v LE N_ELEMENTS(*e.species) - 1 THEN species[v] = $
    FIX((*e.species)[v]) ELSE species[v] = 0
  ; for more species than vars, skip additional species; for less, add zeroes

  ; exception handling for false input:
  w1 = WHERE(species GT gui.out.n_spec_sel)
  w2 = WHERE(species LT 0)
  IF (w1[0] NE -1) OR (w2[0] NE -1) THEN BEGIN
    printerror, 'Incorrect input in spec'
    (*diag).selected = 0
    RETURN
  ENDIF

  sp0 = WHERE(species EQ 0)
  IF sp0[0] NE -1 THEN species[sp0] = -1
  FOR isp = 0, N_ELEMENTS(species) - 1 DO $
    IF species[isp] NE -1 THEN species[isp] = (*gui.out.spec_select)[species[isp]-1]
  ; for each pair of variable one pair of species
  spec_pairs = INTARR(N_ELEMENTS(var_pairs[*,0]),2) 
  FOR v = 0, N_ELEMENTS(species) / 2 - 1 DO spec_pairs[v,*] = species[2*v:2*v+1]
  ; put the input species numbers as pairs into the 2D array 'var_pairs'

  need_phi = 0
  w1 = WHERE(variables EQ 12)
  w2 = WHERE(variables EQ 13)
  IF (w1[0] NE -1) OR (w2[0] NE -1) THEN BEGIN
    fft_format, kxky=0
    need_phi = 1
  ENDIF

  i = set_internal_vars(diag,{$
    var_pairs  : var_pairs,$
    spec_pairs : spec_pairs,$
    ev_solv    : ev_solv,$
    n_ev       : n_ev,$
    need_phi   : need_phi,$
    phi_id     : PTR_NEW(),$
    mom_times  : PTR_NEW()})

END

;######################################################################

PRO fluxratios_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).need_phi THEN BEGIN
    IF NOT (*i).ev_solv OR ((ABS(FIX(mom_time)) NE 0) AND $
      (ABS(FIX(mom_time)) LE (*i).n_ev)) THEN BEGIN

      ; averaging over kx and ky
      phi_rms_z = TOTAL((2.0*TOTAL(ABS($
        (*mom[0,0].kxky)[*,1:par.nky0-1,*])^2,2)+ABS($
        (*mom[0,0].kxky)[*,0,*])^2),1)
      ; z averaging
      phi_rms_z = TEMPORARY(phi_rms_z) * (*series.geom).jac_norm
      phi_rms = SQRT(TOTAL(phi_rms_z)/par.nz0)

      (*i).phi_id = store_step((*i).phi_id,phi_rms)
      (*i).mom_times = store_step((*i).mom_times,mom_time)
    ENDIF
  ENDIF

END

;######################################################################

 PRO fluxratios_output, diag

  COMMON global_vars
 
  i = (*diag).internal_vars

  ; --- get abscissa, i.e.: open scan.log and read e.g. ky values ---

  scan_file = gui.out.data_path + 'scan.log'
  OPENR, scan_lun, scan_file, /GET_LUN, ERROR=err

  IF err NE 0 THEN BEGIN
    PRINT, 'No scan.log found in data path, aborting'
    RETURN
  ENDIF
  arr_size = (*scanlog).n_scans

  ; read header line and get scan dimension
  line = ''
  READF, scan_lun, line
  line = STRTRIM(STRSPLIT(line, '/',/EXTRACT))
  line = STRTRIM(STRSPLIT(line[0],'|',/EXTRACT))
  ; line is array of '#RUN' and all scan values with species number
  dimension = N_ELEMENTS(line) - 1
  x_sp = STRARR(dimension,2)
  ; 0th entry of x_sp is x name, 1st entry is species number
  ; start with '1' (line[0] is '#RUN')
  FOR d = 1, dimension DO x_sp[d-1,*] = STRTRIM(STRSPLIT(line[d],' ',/EXTRACT))

  xaxis = FLTARR(arr_size,dimension)
  n = 0L
  line = ''
  WHILE (not EOF(scan_lun)) DO BEGIN
    READF, scan_lun, line
    IF line NE '' THEN BEGIN
      scanvalue = FLOAT(STRTRIM(STRSPLIT(line,'|',/EXTRACT)))
      ; note: no separator for last column
      FOR d = 0, dimension - 1 DO xaxis[n,d] = scanvalue[d+1]
      n += 1
    ENDIF
  ENDWHILE

  FREE_LUN, scan_lun

  ; --- get data for normalization (for all species marked in the gui) ---

  IF gui.out.n_spec_sel EQ 0 THEN BEGIN
    PRINT, 'No species selected!'
    RETURN
  ENDIF

  ; check number of ratios requested for each pair of variables
  ; '0'th-second array element is for denominator, '1' th is for enumerator:
  ; if all marked species are to be used, n_spec is negative
  ; negative values: traverse marked species; positive: only this species
  ; otherwise one could not distinguish the case of '1 requested' from
  ; '0 requested and only one species marked'
  n_spec = INTARR(N_ELEMENTS((*i).var_pairs[*,0]),2)
  FOR v = 0, N_ELEMENTS((*i).var_pairs[*,0]) - 1 DO BEGIN
    IF (*i).spec_pairs[v,0] EQ -1 THEN $
      n_spec[v,0] = - gui.out.n_spec_sel ELSE n_spec[v,0] = 1
    IF (*i).spec_pairs[v,1] EQ -1 THEN $
      n_spec[v,1] = - gui.out.n_spec_sel ELSE n_spec[v,1] = 1
  ENDFOR
   
  ; get ordinate data
  file_path = STRING(gui.out.data_path,'/nrg_*')                                                           
 
  nrg_files = FILE_SEARCH(file_path,/TEST_REGULAR)
  nrg_files = nrg_files[WHERE(1-FILE_TEST(nrg_files,/ZERO_LENGTH))]
  IF nrg_files[0] EQ '' THEN BEGIN
    PRINT,'No nrg file found in data path! -> EXIT'
    RETURN
  ENDIF

  phi_data = (*i).need_phi ? store_step((*i).phi_id,/get) : PTR_NEW()
  mom_times  = (*i).need_phi ? store_step((*i).mom_times,/get) : PTR_NEW()

;-----------------------------------------------------------------------
;get ordinate data: out of .nrg

  IF (*i).n_ev eq 0 THEN ratio = DBLARR(arr_size,N_ELEMENTS((*i).var_pairs[*,1]),MAX(MAX(ABS(n_spec[*,0]))>MAX(ABS(n_spec[*,1]))),1) $
  ;store averaged data (ratio) from the j-th nrg-file (of arr_size possible)
  ;the outer MAX is necessary in the case of equality -> '>' gives an array
              ELSE  ratio = DBLARR(arr_size,N_ELEMENTS((*i).var_pairs[*,0]),MAX(MAX(ABS(n_spec[*,0]))>MAX(ABS(n_spec[*,1]))),(*i).n_ev)  
              ;data for each eigenvalue for ratio from the j-th nrg-file (of arr_size possible)


   nr_files = N_ELEMENTS(nrg_files)
   FOR j = 0, nr_files - 1 DO BEGIN    ;do for each found files -> j is filenumber
     OPENR, nrg_lun, nrg_files[j], /GET_LUN, ERROR=err
  
     nr_lines = FILE_LINES(nrg_files[j])    ;numbers of lines in current file

    IF (*i).ev_solv eq 0 THEN BEGIN
      IF nr_lines le 0 THEN BEGIN ;special case: j-th nrg-file is completely empty
        nrg_data = MAKE_ARRAY(10,par.n_spec,value='!VALUES.D_NAN',1)
        nrg_time = MAKE_ARRAY(1,value=!VALUES.D_NAN)
      ENDIF ELSE BEGIN
        nr_points = nr_lines/(par.n_spec+1)
        nrg_data = DBLARR(10,nr_points,par.n_spec)
        nrg_time = DBLARR(nr_points)
        tline = 0D
        sline = DBLARR(10)
        FOR k = 0L, nr_points - 1 DO BEGIN
          READF, nrg_lun, tline
          nrg_time[k] = tline
          FOR l = 0, par.n_spec-1 DO BEGIN
            READF, nrg_lun, sline
            nrg_data[*,k,l] = sline
          ENDFOR
        ENDFOR
      ENDELSE
    ENDIF ELSE BEGIN
      IF nr_lines le 0 THEN BEGIN ;special case: j-th nrg-file is completely empty
        tabular = DBLARR(8,1) 
        tabular[*,*] = !VALUES.F_NAN
      ENDIF ELSE tabular = DBLARR(8,nr_lines) ;the tabular to read the values of nrg.dat-file 

      tline = 0.0
      sline = FLTARR(8)
  
      FOR k = 0L, nr_lines - 1 DO BEGIN 
        IF (k MOD (par.n_spec+1)) EQ 0 THEN BEGIN
          READF, nrg_lun, tline
          tabular[0,k] = tline
          tabular[1:7,k] = 0.0
        ENDIF ELSE BEGIN
          READF, nrg_lun, sline
          tabular[*,k] = sline
        ENDELSE
      ENDFOR
    ENDELSE


  IF (*i).need_phi THEN BEGIN
    minval = MIN(ABS((*mom_times[j]) - nrg_time),minind) ;find the nrg time point that matches the mom time
    nr_points = minind + 1
    n_average = 1
  ENDIF ELSE n_average = 3

;----------------------------------------------------------------------------------
 ;calculation

 ;If need_phi is true, then need to be smarter about which nrg lines to average
 ;in order to match up to mom_time, so extra noise isn't added

FOR v=0, (N_ELEMENTS((*i).var_pairs[*,0])-1) DO BEGIN  ; v-Loop
;v is the run index over each pair of requested variables the ratio
;has to be calculated from

 IF (*i).ev_solv eq 0 THEN BEGIN ; ---------------------------------------------
  dividend = DBLARR(MAX(ABS(n_spec[*,0])),1)
  divisor = DBLARR(MAX(ABS(n_spec[*,1])),1)
  ; If a species is not specified, then use the same species

  SWITCH (*i).var_pairs[v,0] OF   ;average the data (enumerator)
    0:
    1:
    2:
    3:
    4:
    5:
    6:
    7: BEGIN
        FOR isp=0,ABS(n_spec[v,0])-1 DO BEGIN
          speci = (*i).spec_pairs[v,0] GE 0 ? (*i).spec_pairs[v,0] : (*gui.out.spec_select)[isp]
          dividend[isp] = TOTAL(nrg_data[(*i).var_pairs[v,0],nr_points-n_average:nr_points-1,speci],2)/n_average
        ENDFOR
        BREAK
       END
     8:
     9:
    10:
    11: BEGIN
        FOR isp=0,ABS(n_spec[v,0])-1 DO BEGIN
          speci = (*i).spec_pairs[v,0] GE 0 ? (*i).spec_pairs[v,0] : (*gui.out.spec_select)[isp]
          dividend[isp] = TOTAL(nrg_data[(*i).var_pairs[v,0]-4,nr_points-n_average:nr_points-1,speci],2)/n_average
        ENDFOR
         BREAK
        END
    12: BEGIN
         dividend[*,*] = (*phi_data[j])^2
         BREAK
        END
    13: BEGIN
         FOR isp=0,ABS(n_spec[v,0])-1 DO BEGIN
           speci = (*i).spec_pairs[v,0] GE 0 ? (*i).spec_pairs[v,0] : (*gui.out.spec_select)[isp]
           dividend[isp] = TOTAL(nrg_data[0,nr_points-n_average:nr_points-1,speci],2)/n_average
         ENDFOR
         dividend = (*phi_data[j])*SQRT(dividend)
         BREAK
        END
    ELSE: PRINT, 'Invalid enumerator number!', RETURN
  
  ENDSWITCH
  
  ; If a species is not specified, then use the same species
  SWITCH (*i).var_pairs[v,1] OF   ;average the data (denominator)
    0:
    1:
    2:
    3:
    4:
    5:
    6:
    7: BEGIN
        FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN
          speci = (*i).spec_pairs[v,1] GE 0 ? (*i).spec_pairs[v,1] : (*gui.out.spec_select)[isp]
          divisor[isp] = TOTAL(nrg_data[(*i).var_pairs[v,1],nr_points-n_average:nr_points-1,speci],2)/n_average
        ENDFOR
        BREAK
       END
     8:
     9:
    10:
    11: BEGIN
        FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN
          speci = (*i).spec_pairs[v,1] GE 0 ? (*i).spec_pairs[v,1] : (*gui.out.spec_select)[isp]
          divisor[isp] = TOTAL(nrg_data[(*i).var_pairs[v,1]-4,nr_points-n_average:nr_points-1,speci],2)/n_average
        ENDFOR
        BREAK
       END
    12: BEGIN
         divisor[*,*] = (*phi_data[j])^2
         BREAK
        END
    13: BEGIN
         FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN
           speci = (*i).spec_pairs[v,1] GE 0 ? (*i).spec_pairs[v,1] : (*gui.out.spec_select)[isp]
           divisor[isp] = TOTAL(nrg_data[0,nr_points-n_average:nr_points-1,(*i).spec_pairs[v,1]],2)/n_average
         ENDFOR
         divisor = (*phi_data[j])*SQRT(divisor)
         BREAK
        END
    ELSE: PRINT, 'Invalid denominator number!', RETURN 
 ENDSWITCH
ENDIF ;no eigenvalue solver ------------------------------------------------


IF (*i).ev_solv eq 1 THEN BEGIN ;eigenvalue solver ----------------------------

   dividend = DBLARR(MAX(ABS(n_spec[*,0])),(*i).n_ev)  
   divisor = DBLARR(MAX(ABS(n_spec[*,1])),(*i).n_ev)
   ;now 2D-arrays for each eigenvalue one column, too

  SWITCH (*i).var_pairs[v,0] OF ;get the data sorted by species and eigenvalue (enumerator)
    0:
    1:
    2:
    3:
    4:
    5:
    6:
    7: BEGIN
        counter = 0    ;now: number of the eigenvalue
        FOR n=0,nr_lines-1 DO BEGIN
           IF ((n mod (par.n_spec+1)) eq 0) THEN counter = BYTE(ABS(tabular[0,n]))
           IF counter gt (*i).n_ev THEN BREAK
           IF n_spec[v,0] eq 1 THEN BEGIN
           ;only one species is required
             isp=0
             IF ((n mod (par.n_spec+1)) eq ((*i).spec_pairs[v,0]+1)) THEN $
               IF counter ne 0 THEN dividend[isp,(counter-1)]=tabular[(*i).var_pairs[v,0],n]
               ;the ev with number '0' is initialization for initial value solver
           ENDIF ELSE BEGIN
           ;all marked species are required
             FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN 
               IF ((n mod (par.n_spec+1)) eq ((*gui.out.spec_select)[isp]+1)) THEN $
                 IF counter ne 0 THEN dividend[isp,(counter-1)]=tabular[(*i).var_pairs[v,0],n]
             ENDFOR
           ENDELSE
        ENDFOR
        BREAK
       END
    8:
    9:
   10:
   11: BEGIN
        counter = 0    ;now: number of the eigenvalue
        FOR n=0,nr_lines-1 DO BEGIN
           IF ((n mod (par.n_spec+1)) eq 0) THEN counter = BYTE(ABS(tabular[0,n]))
           IF counter gt (*i).n_ev THEN BREAK
           IF n_spec[v,0] eq 1 THEN BEGIN
           ;only one species is required
             isp=0
             IF ((n mod (par.n_spec+1)) eq ((*i).spec_pairs[v,0]+1)) THEN $
               IF counter ne 0 THEN dividend[isp,(counter-1)]=tabular[(*i).var_pairs[v,0]-4,n]
               ;the ev with number '0' is initialization for initial value solver
            ENDIF ELSE BEGIN
           ;all marked species are required
             FOR isp=0,ABS(n_spec[v,0])-1 DO BEGIN 
               IF ((n mod (par.n_spec+1)) eq ((*gui.out.spec_select)[isp]+1)) THEN $
                 IF counter ne 0 THEN dividend[isp,(counter-1)]=tabular[(*i).var_pairs[v,0]-4,n]
             ENDFOR
           ENDELSE
        ENDFOR
        BREAK
       END
    12:BEGIN
        FOR isp=0,MAX(ABS(n_spec[*,0]))-1 DO dividend[isp,*] = (*phi_data[j]) ;phi is all the same for each species!
        BREAK
        END
    13:BEGIN
        counter = 0    ;now: number of the eigenvalue
        FOR n=0,nr_lines-1 DO BEGIN
           IF ((n mod (par.n_spec+1)) eq 0) THEN counter = BYTE(ABS(tabular[0,n]))
           IF counter gt (*i).n_ev THEN BREAK
           IF n_spec[v,0] eq 1 THEN BEGIN
           ;only one species is required
             isp=0
             IF ((n mod (par.n_spec+1)) eq ((*i).spec_pairs[v,0]+1)) THEN $
               IF counter ne 0 THEN dividend[isp,(counter-1)]=tabular[0,n]
               ;the ev with number '0' is initialization for initial value solver
           ENDIF ELSE BEGIN
           ;all marked species are required
             FOR isp=0,ABS(n_spec[v,0])-1 DO BEGIN 
               IF ((n mod (par.n_spec+1)) eq ((*gui.out.spec_select)[isp]+1)) THEN $
                 IF counter ne 0 THEN dividend[isp,(counter-1)]=tabular[0,n]
             ENDFOR
           ENDELSE
       ENDFOR
        FOR isp=0,MAX(ABS(n_spec[*,0]))-1 DO dividend[isp,*] = (*phi_data[j])*SQRT(dividend[isp,*])
        ;phi is all the same for each species!
        BREAK
       END
    ELSE: PRINT, 'Invalid enumerator number!', RETURN 
 ENDSWITCH

 SWITCH (*i).var_pairs[v,1] OF ;get the data sorted by species and eigenvalue (denominator)
    0:
    1:
    2:
    3:
    4:
    5:
    6:
    7: BEGIN
        counter = 0    ;now: number of the eigenvalue
        FOR n=0,nr_lines-1 DO BEGIN
           IF ((n mod (par.n_spec+1)) eq 0) THEN counter = BYTE(ABS(tabular[0,n]))
           IF counter gt (*i).n_ev THEN BREAK
           IF n_spec[v,0] eq 1 THEN BEGIN
           ;only one species is required
             isp=0
             IF ((n mod (par.n_spec+1)) eq ((*i).spec_pairs[v,1]+1)) THEN $
               IF counter ne 0 THEN divisor[isp,(counter-1)]=tabular[(*i).var_pairs[v,1],n]
               ;the ev with number '0' is initialization for initial value solver
            ENDIF ELSE BEGIN
           ;all marked species are required
             FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN 
               IF ((n mod (par.n_spec+1)) eq ((*gui.out.spec_select)[isp]+1)) THEN $
                 IF counter ne 0 THEN divisor[isp,(counter-1)]=tabular[(*i).var_pairs[v,1],n]
             ENDFOR
           ENDELSE
        ENDFOR
        BREAK
       END
    8:
    9:
   10:
   11: BEGIN
        counter = 0    ;now: number of the eigenvalue
        FOR n=0,nr_lines-1 DO BEGIN
           IF ((n mod (par.n_spec+1)) eq 0) THEN counter = BYTE(ABS(tabular[0,n]))
           IF counter gt (*i).n_ev THEN BREAK
           IF n_spec[v,1] eq 1 THEN BEGIN
           ;only one species is required
             isp=0
             IF ((n mod (par.n_spec+1)) eq ((*i).spec_pairs[v,1]+1)) THEN $
               IF counter ne 0 THEN divisor[isp,(counter-1)]=tabular[(*i).var_pairs[v,1]-4,n]
               ;the ev with number '0' is initialization for initial value solver
           ENDIF ELSE BEGIN
           ;all marked species are required
             FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN 
               IF ((n mod (par.n_spec+1)) eq ((*gui.out.spec_select)[isp]+1)) THEN $
                 IF counter ne 0 THEN divisor[isp,(counter-1)]=tabular[(*i).var_pairs[v,1]-4,n]
             ENDFOR
           ENDELSE
        ENDFOR
        BREAK
       END
    12:BEGIN
        FOR isp = 0, MAX(ABS(n_spec[*,1])) - 1 DO divisor[isp,*] = *phi_data[j]
        BREAK
        END
    13:BEGIN
        counter = 0    ;now: number of the eigenvalue
        FOR n=0,nr_lines-1 DO BEGIN
           IF ((n mod (par.n_spec+1)) eq 0) THEN counter = BYTE(ABS(tabular[0,n]))
           IF counter gt (*i).n_ev THEN BREAK
           IF n_spec[v,1] eq 1 THEN BEGIN
           ;only one species is required
             isp=0
             IF ((n mod (par.n_spec+1)) eq ((*i).spec_pairs[v,1]+1)) THEN $
               IF counter ne 0 THEN divisor[0,(counter-1)]=tabular[0,n]
               ;the ev with number '0' is initialization for initial value solver
           ENDIF ELSE BEGIN
           ;all marked species are required
             FOR isp=0,ABS(n_spec[v,1])-1 DO BEGIN 
               IF ((n mod (par.n_spec+1)) eq ((*gui.out.spec_select)[isp]+1)) THEN $
                 IF counter ne 0 THEN divisor[isp,(counter-1)]=tabular[0,n]
             ENDFOR
           ENDELSE
        ENDFOR
        FOR isp = 0, MAX(ABS(n_spec[*,1])) - 1 DO $
          divisor[isp,*] = (*phi_data[j]) * SQRT(divisor[isp,*])
        ;phi is all the same for each species!
        BREAK
       END 
    ELSE: PRINT, 'Invalid denominator number!', RETURN 
 ENDSWITCH

ENDIF ; eigenvalue solver used --------------------------------------------

  SWITCH (*i).var_pairs[v,0] OF      ;now transform D from Gamma and chi from Q (enumerator)
    8:
    9: BEGIN 
         IF (*i).spec_pairs[v,0] ne -1 THEN dividend[0,*]=dividend[0,*]/spec[(*i).spec_pairs[v,0]].omn $
            ELSE FOR isp=0,ABS(n_spec[v,0])-1 DO dividend[isp,*]=dividend[isp,*]/spec[(*gui.out.spec_select)[isp]].omn
         BREAK
       END
   10:
   11: BEGIN
        IF (*i).spec_pairs[v,0] ne -1 THEN dividend[0,*]= $
          dividend[0,*]/(spec[(*i).spec_pairs[v,0]].omt*spec[(*i).spec_pairs[v,0]].dens*spec[(*i).spec_pairs[v,0]].temp) 
        FOR isp=0,ABS(n_spec[v,0])-1 DO dividend[isp,*]= $
          dividend[isp,*]/(spec[(*gui.out.spec_select)[isp]].omt* $
          spec[(*gui.out.spec_select)[isp]].dens*spec[(*gui.out.spec_select)[isp]].temp)
        BREAK
       END
   ELSE: BREAK
  ENDSWITCH

  SWITCH (*i).var_pairs[v,1] OF    ;now transform D from Gamma and chi from Q (denominator)
    8:
    9: BEGIN 
         IF (*i).spec_pairs[v,1] ne -1 THEN divisor[0,*]=divisor[0,*]/spec[(*i).spec_pairs[v,1]].omn $
            ELSE FOR isp=0,ABS(n_spec[v,1])-1 DO divisor[isp,*]=divisor[isp,*]/spec[(*gui.out.spec_select)[isp]].omn
         BREAK
       END
   10:
   11: BEGIN
        IF (*i).spec_pairs[v,1] ne -1 THEN divisor[0,*]= $
          divisor[0,*]/(spec[(*i).spec_pairs[v,1]].omt*spec[(*i).spec_pairs[v,1]].dens*spec[(*i).spec_pairs[v,1]].temp) 
        FOR isp=0,ABS(n_spec[v,1])-1 DO divisor[isp,*]= $
          divisor[isp,*]/(spec[(*gui.out.spec_select)[isp]].omt* $
          spec[(*gui.out.spec_select)[isp]].dens*spec[(*gui.out.spec_select)[isp]].temp)
        BREAK
       END
   ELSE: BREAK
  ENDSWITCH

  FOR isp=0,MAX(MAX(ABS(n_spec[v,0]))>MAX(ABS(n_spec[v,1])))-1 DO BEGIN
   ;take care of possibly different  length of dividend and divisor:
   IF isp gt ABS(n_spec[v,0])-1 THEN ratio[j,v,isp,*]=dividend[0,*]/divisor[isp,*] $
   ELSE IF isp gt ABS(n_spec[v,1])-1 THEN ratio[j,v,isp,*]=dividend[isp,*]/divisor[0,*] $
   ELSE ratio[j,v,isp,*]=dividend[isp,*]/divisor[isp,*]
  ENDFOR

 ENDFOR ;end of v-Loop (traverse variable pairs)

 FREE_LUN, nrg_lun

ENDFOR ;end of j-Loop (traverse files for data of different x-axis-values)

;--------------------------------------------------------------------------------
;output in ps-file

    xaxis_local = xaxis
    ratio_local = ratio


FOR v=0, N_ELEMENTS((*i).var_pairs[*,0])-1 DO BEGIN ;second v-Loop (now for output)

     icon=''
     IF N_ELEMENTS((*i).var_pairs[*,0]) eq 1 THEN BEGIN 
     ;only in case of only one pair of variables,
     ;the file can be named with information about variables
        CASE (*i).var_pairs[v,0] OF    ;enumerator
         8 : CASE (*i).var_pairs[v,1] OF  ;denominator
              8 : icon='D_es-D_es'
              9 : icon='D_es-D_em'
              10: icon='D_es-Chi_es'
              11: icon='D_es-Chi_em'
              12: icon='D_es-phi_2'
              13: icon='D_es-phi_n'
            ELSE: icon='D_es-'+get_nrg_string((*i).var_pairs[v,1]) 
           ENDCASE
         9 : CASE (*i).var_pairs[v,1] OF  ;denominator ...
              8 : icon='D_em-D_es'
              9 : icon='D_em-D_em'
              10: icon='D_em-Chi_es'
              11: icon='D_em-Chi_em'
              12: icon='D_em-phi_2'
              13: icon='D_em-phi_n'
            ELSE: icon='D_em-'+get_nrg_string((*i).var_pairs[v,1]) 
            ENDCASE
         11: CASE (*i).var_pairs[v,1] OF 
              8 : icon='Chi_es-D_es'
              9 : icon='Chi_es-D_es'
              10: icon='Chi_es-D_es'
              11: icon='Chi_es-D_es'
              12: icon='Chi_es-phi_2'
              13: icon='Chi_es-phi_n'
            ELSE: icon='Chi_es-'+get_nrg_string((*i).var_pairs[v,1]) 
           ENDCASE
         10: CASE (*i).var_pairs[v,1] OF 
              8 : icon='Chi_em-D_es' 
              9 : icon='Chi_em-D_em' 
              10: icon='Chi_em-Chi_es' 
              11: icon='Chi_em-Chi_em' 
              12: icon='Chi_em-phi_2'
              13: icon='Chi_em-phi_n' 
            ELSE: icon='Chi_em-'+get_nrg_string((*i).var_pairs[v,1])  
           ENDCASE
         12: CASE (*i).var_pairs[v,1] OF 
              8 : icon='phi_2-D_es' 
              9 : icon='phi_2-D_em' 
              10: icon='phi_2-Chi_es' 
              11: icon='phi_2-Chi_em' 
              12: icon='phi_2-phi_2'
              13: icon='phi_2-phi_n' 
            ELSE: icon='phi_2-'+get_nrg_string((*i).var_pairs[v,1])  
           ENDCASE
         13: CASE (*i).var_pairs[v,1] OF 
              8 : icon='phi_n-D_es' 
              9 : icon='phi_n-D_em' 
              10: icon='phi_n-Chi_es' 
              11: icon='phi_n-Chi_em' 
              12: icon='phi_n-phi_2'
              13: icon='phi_n-phi_n' 
            ELSE: icon='phi_n-'+get_nrg_string((*i).var_pairs[v,1])  
           ENDCASE
         ELSE: CASE (*i).var_pairs[v,1] OF 
              8 : icon= get_nrg_string((*i).var_pairs[v,0])+'-D_es' 
              9 : icon=get_nrg_string((*i).var_pairs[v,0])+'-D_em' 
              10: icon=get_nrg_string((*i).var_pairs[v,0])+'-Chi_es' 
              11: icon=get_nrg_string((*i).var_pairs[v,0])+'-Chi_em'
              12: icon=get_nrg_string((*i).var_pairs[v,0])+'-phi_2'
              13: icon=get_nrg_string((*i).var_pairs[v,0])+'-phi_n' 
            ELSE: icon=get_nrg_string((*i).var_pairs[v,0])+'-'+get_nrg_string((*i).var_pairs[v,1])  
            ENDCASE
        ENDCASE
     ENDIF 

    IF v eq 0 THEN BEGIN ;open ps-file and define some title only once (at the begining of v-loop)
        ;ps_file=''
        ;IF TOTAL(ABS(n_spec[*,0])>ABS(n_spec[*,1])) eq 1 THEN ps_file = STRING('fluxratio') $
        ;;to find the number of necessary plots: number of vars*number ofspecies of the variable
        ;                            ELSE 
        ps_file = 'fluxratio'
        IF (*i).ev_solv eq 0 THEN BEGIN ;layout for initial value solver
          IF TOTAL(ABS(n_spec[*,0])>ABS(n_spec[*,1])) eq 1 THEN $ ;layout depends on number of plots
             set_output, ps_file, /ps, charsize=1.5, xsize=25,ysize=17, multi=[0,1,1]$
          ELSE IF TOTAL(ABS(n_spec[*,0])>ABS(n_spec[*,1])) eq 2 THEN $
             set_output, ps_file, /ps, charsize=1.5, xsize=25,ysize=17, multi=[0,2,1] $
          ELSE set_output, ps_file, /ps, charsize=1.5, xsize=25,ysize=17, multi=[0,2,2]
        ENDIF ELSE BEGIN ;eigenvalue solver used (need space for legend)
                         ;and: if 2-dimensional plot:it is not possible any more to put  all ev into one plot! 
          IF dimension eq 1 THEN BEGIN
             IF TOTAL(ABS(n_spec[*,0])>ABS(n_spec[*,1])) eq 1 THEN $ ;layout depends on number of plots
                set_output, ps_file, /ps, charsize=1.5, xsize=25,ysize=17, multi=[0,1,1] $
             ELSE set_output, ps_file, /ps, charsize=1.5, xsize=25,ysize=17, multi=[0,1,1]     
          ENDIF ELSE BEGIN ;2-dimensional plot it is not possible any more to put  all ev into one plot!
             IF TOTAL(ABS(n_spec[*,0])>ABS(n_spec[*,1]))*(*i).n_ev eq 1 THEN $ ;layout depends on number of plots
                set_output, ps_file, /ps, charsize=2.0, xsize=25,ysize=17, multi=[0,1,1] $
             ELSE set_output, ps_file, /ps, charsize=2.0, xsize=25,ysize=17, multi=[0,2,1]
          ENDELSE  
        ENDELSE

        ;title for x axis:
        x_title=STRARR(dimension)
        FOR d=0,dimension-1 DO BEGIN
          IF STRCMP(x_sp[d,0], 'ky') THEN x_title[d] = '!7q!6!Ds!N * k!Dy!N' ELSE x_title[d] = x_sp[d,0]
          IF STRCMP(x_sp[d,0], 'omn')THEN x_title[d]= '!7x!6!Dn!D'+' '+ spec[BYTE(x_sp[d,1]-1)].name+'!N'
          IF STRCMP(x_sp[d,0], 'omt')THEN x_title[d] = '!7x!6!DT!D'+' '+ spec[BYTE(x_sp[d,1]-1)].name+'!N'
          IF STRCMP(x_sp[d,0], 'beta')THEN x_title[d] = '!7b'
        ENDFOR


        ;title for y-axis and for whole plot title:
        y_title = [get_nrg_string(INDGEN(12), /fancy),$                ;casual variables
            '!12<!9!7U!12>!U2!N!6','!12<!9!7U!12>!12<!6n!12>!6']       ;phi^2 und phi*n

        y_title2 = [get_nrg_string(INDGEN(12)),$        ;casual variables
            '<phi^2>','<phi*n>']                       ;phi^2 und phi*n

           ;'!12<!6D!9!6!Ies!N!12>!6', '!12<!6D!9!6!Iem!N!12>!6',$     ;D_es, D_em
           ;'!12<!9!7V!9!6!Ies!N!12>!6', '!12<!9!7V!9!6!Iem!N!12>!6',$ ;Chi_es, Chi_em

   ENDIF
   
   IF dimension eq 1 THEN BEGIN  ;only one-dimensional plot----------------------
    IF (*i).ev_solv eq 0 THEN BEGIN ;initial value solver
       IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN ;'loop' over species_pairs
           PLOT, xaxis_local, ratio_local[*,v,0,0],$
             TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
             '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name,$
             XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
             /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,0,*],/NAN) , MAX(ratio_local[*,*],/NAN)]
           OPLOT,  xaxis_local, ratio_local[*,v,0,0], COLOR=2+v, PSYM=-7 
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN
           FOR isp=0,(ABS(n_spec[v,0]))-1 DO BEGIN  ;loop over species_pairs
             PLOT, xaxis_local, ratio_local[*,v,isp,0],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name,$
               XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             OPLOT,  xaxis_local, ratio_local[*,v,isp,0], COLOR=2+v, PSYM=-7 
         ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN
           FOR isp=0,(ABS(n_spec[v,1]))-1 DO BEGIN
             PLOT, xaxis_local, ratio_local[*,v,isp,0],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name,$
               XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             OPLOT, xaxis_local, ratio_local[*,v,isp,0], COLOR=2+v, PSYM=-7 
           ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN
           FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN
             PLOT, xaxis_local, ratio_local[*,v,isp,0],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name,$
               XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             OPLOT,  xaxis_local, ratio_local[*,v,isp,0], COLOR=2+v, PSYM=-7 
           ENDFOR
       ENDIF
       ;now: legend plotting (only one legend per page!)
       ;IF (v mod (v<4)) eq ((v<4)-1) THEN plot_legend, color_arr,legend_arr, per_line = 2
    ENDIF ELSE BEGIN  ;eigenvalue solver  
       !Y.MARGIN=[9,2] ;to give space for legend
       !X.MARGIN=[11,2]
       IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN ;no loop over species
           PLOT, xaxis_local, ratio_local[*,v,0,0],$
             TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
             '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name,$
             XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
             /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,0,*],/NAN) , MAX(ratio_local[*,v,0,*],/NAN)]
           FOR counter = 0, (*i).n_ev-1 DO $  ;loop over eigenvalues
             ;only show the requested eigenvalues, not supplementary foundones,as these data are often -NaN
             OPLOT,  xaxis_local, ratio_local[*,v,0,counter], COLOR=2+counter, PSYM=-7   
        ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN ;loop over species of enum
           FOR isp=0,(ABS(n_spec[v,0]))-1 DO BEGIN
             IF MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1])) ge 3 THEN !y.margin=[4,2]
             IF ((isp mod 4) eq 2) or ((isp mod 4) eq 3) THEN !y.margin=[10,2] 
             ;only the both pictures on the bottom have another ymargin
             PLOT, xaxis_local, ratio_local[*,v,isp,0],POS=[0.1,0.35,0.95,0.9],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name,$
               XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             FOR counter = 0, (*i).n_ev-1 DO $
             ;only show the requested eigenvalues, not supplementary foundones,as these data are often -NaN
               OPLOT,  xaxis_local, ratio_local[*,v,isp,counter], COLOR=2+counter, PSYM=-7       
         ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN ;loop over species of denom
           FOR isp=0,(ABS(n_spec[v,1]))-1 DO BEGIN
             IF MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1])) ge 3 THEN !y.margin=[4,2]
             IF ((isp mod 4) eq 2) or ((isp mod 4) eq 3) THEN !y.margin=[10,2] 
             ;only the both pictures on the bottom have another ymargin
             PLOT, xaxis_local, ratio_local[*,v,isp,0],POS=[0.1,0.35,0.95,0.9],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name,$
               XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             FOR counter = 0, (*i).n_ev-1 DO $
             ;only show the requested eigenvalues, not supplementary foundones,as these data are often -NaN
              OPLOT,  xaxis_local, ratio_local[*,v,isp,counter], COLOR=2+counter, PSYM=-7      
           ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN ;loop over species
           FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN
             IF MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1])) ge 3 THEN !y.margin=[4,2]
             IF ((isp mod 4) eq 2) or ((isp mod 4) eq 3) THEN !y.margin=[10,2] 
             ;only the both pictures on the bottom have another ymargin
             PLOT, xaxis_local, ratio_local[*,v,isp,0],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name,$
               XTITLE= x_title,  YTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               /NODATA, COLOR=1,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
               FOR counter = 0, (*i).n_ev-1 DO $
             ;only show the requested eigenvalues, not supplementary foundones,as these data are often -NaN
               OPLOT,  xaxis_local, ratio_local[*,v,isp,counter], COLOR=2+counter, PSYM=-7        
           ENDFOR
       ENDIF
     ;now for legend plotting:
     color_arr=INDGEN((*i).n_ev)+2
     modenr_arr=INDGEN((*i).n_ev)
     legend_arr='mode '+STRING(modenr_arr)
     plot_legend,color_arr,legend_arr, per_line = 2
     !Y.MARGIN=[4,2]  ;reset margins!
     !X.MARGIN=[10,3]
 ENDELSE ;end eigenvalue solver
ENDIF ;end one dimensional plot

IF dimension eq 2 THEN BEGIN ;two dimensional plot-----------------------------
    ;create the data structure required by SURFACE---------------

    xcount=1 ;get number of values to be read in in x and y direction (i.e. for first and second scan-parameter)
    j=1
    WHILE xaxis_local[j,0] ne xaxis_local[0,0] DO BEGIN 
        xcount=xcount+1
        j=j+1
        if j eq arr_size -1 THEN break ;just for security reasons
    ENDWHILE
    ycount=arr_size/xcount ;this is because arr_size always has to be xcount*ycount!

    ;create data structures
    IF (*i).ev_solv eq 0 THEN begin
      ;ratio_local2=DBLARR(xcount,ycount,N_ELEMENTS((*i).var_pairs[*,0]),MAX(MAX(ABS(n_spec[*,0]))>(ABS(n_spec[*,1]))),1)
      ratio_local2=REFORM(ratio_local,xcount,ycount,N_ELEMENTS((*i).var_pairs[*,0]),MAX(MAX(ABS(n_spec[*,0]))>(ABS(n_spec[*,1]))),1)
    ENDIF ELSE BEGIN
      ratio_local2=DBLARR(xcount,ycount,N_ELEMENTS((*i).var_pairs[*,0]),MAX(MAX(ABS(n_spec[*,0]))>(ABS(n_spec[*,1]))),(*i).n_ev)
      FOR counter = 0, (*i).n_ev-1 DO BEGIN ;loop over eigenvalues
        FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN      
          j=0                     ;fill z-structure
          FOR j1=0,ycount-1 DO BEGIN
              FOR j0=0,xcount-1 DO BEGIN
                  ratio_local2[j0,j1,v,isp,counter]=ratio_local[j,v,isp,counter]
                  j=j+1
              ENDFOR
          ENDFOR
        ENDFOR
      ENDFOR
    ENDELSE

    x=DBLARR(xcount)            ;fill x-structure
    x=xaxis_local[0:xcount-1,0]

    y=DBLARR(ycount)            ;fill y-structure
    y[0]=xaxis_local[0,1]       ;take the first
    j2=0
    FOR j=1, arr_size-1 DO BEGIN
        IF xaxis_local[j,1] ne xaxis_local[j-1,1] THEN BEGIN
            j2=j2+1
            y[j2]=xaxis_local[j,1] ;take only if new value
        ENDIF
    ENDFOR
    
;----------------------------------------------------------------------------

    IF (*i).ev_solv eq 0 THEN BEGIN ;initial value solver
       IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN ;no loop over species
          FOR isp=0,(ABS(n_spec[v,0]))-1 DO BEGIN  ;loop over species_pairs
            SURFACE, ratio_local2[*,*,v,0,0], x,y,$
             TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
             '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name,$
             XTITLE= x_sp[0,0],YTITLE=x_sp[1,0], ZTITLE= y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
             COLOR=2,ZRANGE=[MIN(ratio_local[*,0,*],/NAN) , MAX(ratio_local[*,*],/NAN)]
          ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN
           FOR isp=0,(ABS(n_spec[v,0]))-1 DO BEGIN  ;loop over species_pairs
              SURFACE, ratio_local2[*,*,v,isp,0], x,y,$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name,$
               XTITLE= x_title[0], YTITLE=x_title[1], ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               COLOR=2;,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             ;OPLOT,  xaxis_local, ratio_local[*,v,isp,0], COLOR=2+v, PSYM=-7 
         ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN
           FOR isp=0,(ABS(n_spec[v,1]))-1 DO BEGIN  
             SURFACE, ratio_local2[*,*,v,isp,0], x,y,$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name,$
               XTITLE= x_title[0], YTITLE=x_title[1], ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               COLOR=2 ;,YRANGE=[MIN(ratio_local[*,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             ;OPLOT,  xaxis_local, ratio_local[*,v,isp,0], COLOR=2+v, PSYM=-7 
           ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN
           FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN
              SURFACE, ratio_local2[*,*,v,isp,0], x,y,$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name,$
               XTITLE=x_title[0], YTITLE=x_title[1],  ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               COLOR=2;,YRANGE=[MIN(ratio_local[*,v,isp,*],/NAN) , MAX(ratio_local[*,v,isp,*],/NAN)]
             ;OPLOT,  xaxis_local, ratio_local[*,v,isp,0], COLOR=2+v, PSYM=-7 
           ENDFOR
       ENDIF
       ;now: legend plotting (only one legend per page!)
       ;IF (v mod (v<4)) eq ((v<4)-1) THEN plot_legend, color_arr,legend_arr, per_line = 2
    ENDIF ELSE BEGIN  ;eigenvalue solver  -----------------------------------------------
       IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN ;no loop over species
          FOR isp=0,(ABS(n_spec[v,0]))-1 DO BEGIN  ;loop over species_pairs
           FOR counter = 0, (*i).n_ev-1 DO BEGIN  ;loop over eigenvalues
             ;T3D, /RESET,  TRANSLATE=[-0.5,-0.5,-0.5],ROTATE=[-40,40,20]
             ;T3D, TRANSLATE=[0.5,0.5,0.5]
             SURFACE, ratio_local2[*,*,v,isp,counter], x,y,$
              TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
              '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name+'    mode '+STRING(counter),$
              XTITLE= x_title[0], YTITLE=x_title[1], ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
              COLOR=counter+2, XRANGE=[MIN(x,/NAN) , MAX(x,/NAN)], YRANGE=[MIN(y,/NAN), MAX(y,/NAN)], $
              ZRANGE=[MIN(ratio_local[*,v,isp,counter],/NAN) , MAX(ratio_local[*,v,isp,counter],/NAN)]
          ENDFOR
         ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] ne -1) THEN BEGIN ;loop over species of enum
          FOR isp=0,(ABS(n_spec[v,0]))-1 DO BEGIN
           FOR counter = 0, (*i).n_ev-1 DO BEGIN  ;loop over eigenvalues
               SURFACE, ratio_local2[*,*,v,isp,counter], x,y,POS=[0.1,0.35,0.95,0.9],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*i).spec_pairs[v,1]].name+'    mode '+STRING(counter),$
               XTITLE= x_title[0], YTITLE=x_title[1], ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               COLOR=counter+2
          ENDFOR 
         ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] ne -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN ;loop over species of denom
          FOR isp=0,(ABS(n_spec[v,1]))-1 DO BEGIN
           FOR counter = 0, (*i).n_ev-1 DO BEGIN  ;loop over eigenvalues
              SURFACE, ratio_local2[*,*,v,isp,counter],x,y,POS=[0.1,0.35,0.95,0.9],$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*i).spec_pairs[v,0]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name+'    mode '+STRING(counter),$
               XTITLE= x_title[0], YTITLE=x_title[1], ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               COLOR=counter+2

          ENDFOR
         ENDFOR
       ENDIF ELSE IF ((*i).spec_pairs[v,0] eq -1) and ((*i).spec_pairs[v,1] eq -1) THEN BEGIN ;loop over species
          FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN
           FOR counter = 0, (*i).n_ev-1 DO BEGIN  ;loop over eigenvalues
              SURFACE, ratio_local2[*,*,v,isp,counter],x,y,$
               TITLE= y_title[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
               '/'+y_title[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name+'    mode '+STRING(counter),$
               XTITLE= x_title[0], YTITLE=x_title[1], ZTITLE = y_title[(*i).var_pairs[v,0]]+'/'+y_title[(*i).var_pairs[v,1]],$
               COLOR=counter+2
          ENDFOR
         ENDFOR
       ENDIF
     ENDELSE ;end eigenvalue solver
 ENDIF ;end two dimensional plot
ENDFOR ;second v-Loop

;-------------------------------------------------------------------------------------
;output in dat-file

IF (*i).ev_solv eq 0 THEN BEGIN ;initial value solver
  ;create header array for set_output
  header_arr=STRARR(dimension+1)
  FOR d=0, dimension-1 DO BEGIN
      header_arr[d]=x_sp[d,0]
      IF (x_sp[d,1] ne x_sp[d,0]) THEN header_arr[d] += '_'+spec[BYTE(x_sp[d,1]-1)].name   ;not all scan-parameters do have a species number
  ENDFOR 
  header_arr[dimension]='ratio'

  FOR v=0, (N_ELEMENTS((*i).var_pairs[*,0])-1) DO BEGIN  ; v-Loop
    FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN ;isp-Loop


         ;create array for values for set_output
          dat_arr = DBLARR(arr_size,dimension+1)
          FOR d=0, dimension-1 DO BEGIN
              dat_arr[*,d]=xaxis_local[*,d]
          ENDFOR
          dat_arr[*,dimension]= ratio_local[*,v,isp,0]
 
          spec_pair = reform((*i).spec_pairs[v,*])
          if spec_pair[0] eq -1 then spec_pair[0] = isp
          if spec_pair[1] eq -1 then spec_pair[1] = isp
          set_output,diag, isp, header= header_arr,$
                   commentline= ' RATIO=   '+y_title2[(*i).var_pairs[v,0]]+'_'+spec[spec_pair[0]].name+$
                                '/'+y_title2[(*i).var_pairs[v,1]]+'_'+spec[spec_pair[1]].name,$
                   dat=dat_arr ,append=(n_elements(var_pairs) gt 2)
    ENDFOR ;species-loop (isp)
  ENDFOR ;variable-pairs-loop (v)

ENDIF

IF (*i).ev_solv ne 0 THEN BEGIN ;eigen value solver
  ;create header array for set_output
  header_arr=STRARR(dimension+1)
  FOR d=0, dimension-1 DO BEGIN
      header_arr[d]=x_sp[d,0]
      IF (x_sp[d,1] ne x_sp[d,0]) THEN header_arr[d] += '_'+spec[BYTE(x_sp[d,1]-1)].name   ;not all scan-parameters do have a species number
  ENDFOR 
  header_arr[dimension]='ratio'

  FOR v=0, (N_ELEMENTS((*i).var_pairs[*,0])-1) DO BEGIN  ; v-Loop
    FOR isp=0,MAX(ABS(n_spec[v,0])>ABS(n_spec[v,1]))-1 DO BEGIN ;isp-Loop
      FOR counter = 0, (*i).n_ev-1 DO BEGIN  ;loop over eigenvalues 

         ;create array for values for set_output
          dat_arr = DBLARR(arr_size,dimension+1)
          FOR d=0, dimension-1 DO BEGIN
              dat_arr[*,d]=xaxis_local[*,d]
          ENDFOR
          dat_arr[dimension]= ratio_local[*,v,isp,counter]
 
          set_output,diag, isp, header= header_arr,$
                   commentline= ' RATIO=   '+y_title2[(*i).var_pairs[v,0]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
                                '/'+y_title2[(*i).var_pairs[v,1]]+'_'+spec[(*gui.out.spec_select)[isp]].name+$
                                '    number of eigenvalue = '+ STRING(counter),$
                   dat=dat_arr ,append=((counter gt 0) or (isp gt 0) or (v gt 0)) 
      ENDFOR ;eigenvalue-loop (counter)
    ENDFOR ;species-loop (isp)
  ENDFOR ;variable-pairs-loop (v)
ENDIF

  set_output, 'fluxratio',/reset

  FOR j = 0, N_ELEMENTS(phi_data) - 1 DO $
    IF PTR_VALID(phi_data[j]) THEN PTR_FREE, phi_data[j]

END
