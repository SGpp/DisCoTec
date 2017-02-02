;#############################################################################
;#    this library contains almost all energy related internal functions     #
;#############################################################################

FUNCTION energy_renorm, var, Lref=Lref
;The normalization for the 3D quantities should be as follows:
;(n_ref^2 m_ref T_ref) (rho_ref/L_ref)^2 (n_0j^2 B0 T_0j m_j) |g_j|^2/
;F_0j dv_|| d_mu
;The total energy is the same thing summed over kx,ky and integrated
;over z.The time derivative quantities are the same thing with units of a time
;derivative.

;Note that n_0j, T_0j, m_j etc. are _NOT_ included in the energy arrays
;and -- technically -- still need to be considered (par.Xref)

  COMMON global_vars

  IF NOT KEYWORD_SET(Lref) THEN Lref = series.Lref

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref
  c_ref = SQRT(series.Tref*series.Qref/series.mref)

  E_ref = series.nref^2 * series.mref * series.Qref * series.Tref * $
    (rho_ref/series.Lref)^2

  CASE var OF
    5	 : RETURN, E_ref * c_ref / series.Lref ; NL
    4	 : RETURN, E_ref * c_ref / series.Lref ; dE/dt_hyp_z,v,x,y
    3	 : RETURN, E_ref * c_ref / series.Lref ; dE/dt_coll
    2	 : RETURN, E_ref * c_ref / series.Lref ; dE/dt_drive
    1	 : RETURN, E_ref * c_ref / series.Lref ; dE/dt_noncons
    0	 : RETURN, E_ref ; E
    ELSE : printerror, 'error in energy_renorm (var ' + rm0es(var) + ' not found)'
  ENDCASE

END

;#############################################################################

PRO create_energy_struct

  COMMON global_vars
  
  IF TOTAL(series.request_energy) EQ 0 THEN RETURN
  request_dat = series.request_energy ;struct. tags cannot be passed & modified
  check_request_dat, request_dat
  series.request_energy = request_dat

  energy = REPLICATE({kxky:PTR_NEW(),kxsy:PTR_NEW(),$
    sxky:PTR_NEW(),sxsy:PTR_NEW()},6)        ;num_energies

  create_3darr_struct,energy,series.request_energy

  IF !QUIET NE 1 THEN BEGIN
    PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'
    PRINT, '    current memory usage: ' + mem_usage()
  ENDIF

END

;##########################################################################

PRO destroy_energy_struct, free_eterms=free_eterms

  COMMON global_vars

  IF KEYWORD_SET(free_eterms) THEN BEGIN
    IF PTR_VALID(eterms) THEN PTR_FREE, eterms
    IF PTR_VALID(eterms_time) THEN PTR_FREE, eterms_time
    RETURN
  ENDIF

  PTR_FREE, energy.kxky, energy.sxky, energy.kxsy, energy.sxsy

END

;##########################################################################

PRO fft_format_en, kxky=kxky, sxky=sxky, kxsy=kxsy, sxsy=sxsy, terms=terms
; Called by the _init part of a diagnostic to indicate that it needs
; fourier or real space data.
; request_energy: 4,4 array (field,transformation type), 1 = reqd, 0 = not reqd

  COMMON global_vars

  FOR i = 0, N_ELEMENTS(kxky) - 1 DO series.request_energy[kxky[i],0] = 1
  FOR i = 0, N_ELEMENTS(sxky) - 1 DO series.request_energy[sxky[i],1] = 1
  FOR i = 0, N_ELEMENTS(kxsy) - 1 DO series.request_energy[kxsy[i],2] = 1
  FOR i = 0, N_ELEMENTS(sxsy) - 1 DO series.request_energy[sxsy[i],3] = 1

  IF KEYWORD_SET(terms) THEN series.request_eterms = 1

END

;##########################################################################

PRO calc_energy_ffts
; Calculates sxky, kxsy, sxsy from kxky in every time step for all
; requested fields
; sxsy: standard via sxky, only if kxsy (not sxky) requested via kxsy
; NOTE: [0,0,0] on the lhs is used as a faster alternative to [*,*,*]

  COMMON global_vars

  calc_ffts_3darr,energy,series.request_energy

END

;##########################################################################

FUNCTION get_energy_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  start_time = gui.out.start_t

  IF KEYWORD_SET(last) THEN BEGIN
    istep_energy3d = - (*par.istep_energy3d)[run]
  ENDIF ELSE BEGIN ; for /last, use finer scale file, else use coarser
    istep_energy3d = (*par.istep_energy3d)[run]
  ENDELSE

  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
  file_energy = get_file_path('energy3d',run)

  IF NOT file_energy.exist THEN RETURN, 0.0

  request_energy = TOTAL(series.request_energy[INDGEN(par.n_energies),*],2) ;n_energies

  OPENR, tdata_lun, file_energy.path, /GET_LUN, ERROR=err, /F77_UNFORMATTED, $
    SWAP_ENDIAN=series.swap_endian

  IF err NE 0 THEN BEGIN
    printerror, 'Please notify the developers about this error'
    RETURN, 0L
  ENDIF

  IF par.in_data NE 'sykx' THEN $
     content_byte_size = (4LL * (2 - par.prec_single) * par.nkx0 * (par.nky0 - $
          par.ky0_ind) * par.nz0 + 8) * par.n_energies $
    ELSE content_byte_size = (4LL * (2 - par.prec_single) * (par.nkx0/2) * par.ny0 $
          * par.nz0 + 8) * par.n_energies

  step_byte_size = content_byte_size + 16 - 4 * par.prec_single ; data block + time entry

  time = par.prec_single ? 0.0 : 0.0D
  time_first_step = time
  time_second_step = time
  time_last_step = time

  pos_last_step = (FSTAT(tdata_lun)).SIZE - step_byte_size
  n_steps = ROUND((FSTAT(tdata_lun)).SIZE/step_byte_size)
  IF KEYWORD_SET(get_n_steps) THEN BEGIN
    FREE_LUN, tdata_lun
    RETURN, n_steps
  ENDIF    
  READU, tdata_lun, time_first_step
  time_first_step *= time_renorm(0)

  IF KEYWORD_SET(first) THEN BEGIN
    FREE_LUN, tdata_lun  
    RETURN, time_first_step
  ENDIF    
  
  POINT_LUN, - tdata_lun, position_before
  position_after = position_before + content_byte_size
  POINT_LUN, tdata_lun, position_after
  
  IF NOT EOF(tdata_lun) THEN BEGIN
    READU, tdata_lun, time_second_step
    time_second_step *= time_renorm(0)
  ENDIF
  
  POINT_LUN, tdata_lun, pos_last_step
  READU, tdata_lun, time_last_step
  time_last_step *= time_renorm(0)
  
  IF KEYWORD_SET(last) THEN BEGIN
    FREE_LUN, tdata_lun
    RETURN, time_last_step
  ENDIF    
  avg_step_diff = (time_last_step - time_first_step) / FLOAT(n_steps)

  IF time_first_step GT start_time THEN BEGIN 
    FREE_LUN, tdata_lun
    RETURN, 0L
  ENDIF
  
  IF time_last_step LT start_time THEN BEGIN
    IF !QUIET NE 1 THEN PRINT, 'run ', run_label, ' not in time range'
    FREE_LUN, tdata_lun
    RETURN, -2L
  ENDIF

  IF time_last_step EQ time_first_step THEN BEGIN
    PRINT, 'warning: run ', run_label, ' has only one time step'
    FREE_LUN, tdata_lun
    RETURN, 0L
  ENDIF

  IF time_second_step GT start_time THEN BEGIN 
    FREE_LUN, tdata_lun
    RETURN, 0L
  ENDIF
  
  ; get next step:
  estimated_step = ROUND((start_time-time_first_step)/avg_step_diff)
  estimated_pos_init = step_byte_size * estimated_step - $
    step_byte_size < (pos_last_step - step_byte_size)
;PRINT, 'estimated:', estimated_step, estimated_pos_init

  found = 0
  iterations = 0
  max_iterations = 20
  estimated_pos = estimated_pos_init

  WHILE (found EQ 0) AND (iterations LE max_iterations) AND (estimated_pos GE 0) DO BEGIN
    IF estimated_pos GE 0 THEN POINT_LUN, tdata_lun, estimated_pos ELSE BEGIN
      PRINT, 'estimated_pos LT 0, file position error!'
      FREE_LUN, tdata_lun
      RETURN, -1L
    ENDELSE
    time_prev = par.prec_single ? 0.0 : 0.0D
    READU, tdata_lun, time_prev
    time_prev = time_prev  * time_renorm(0)
    POINT_LUN, - tdata_lun, position_before
    position_after = position_before + content_byte_size
    POINT_LUN, tdata_lun, position_after
    time_act = par.prec_single ? 0.0 : 0.0D
    READU, tdata_lun, time_act
    time_act = time_act  * time_renorm(0)
    time_diff_act = time_act - time_prev
    IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
    new_estimated_pos = step_byte_size * $
      ROUND((start_time - time_act) / time_diff_act) + estimated_pos
    IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
      new_estimated_pos = new_estimated_pos + step_byte_size
    estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)
    IF estimated_pos LT 0 THEN estimated_pos = 0L
;    IF found EQ 0 THEN PRINT, 'new estimate:', estimated_pos, (pos_last_step - step_byte_size)
    iterations = iterations + 1
  ENDWHILE
  IF found EQ 0 THEN BEGIN
    PRINT, 'using simplified search pattern'
    estimated_pos = estimated_pos_init
    WHILE found EQ 0 DO BEGIN
      IF estimated_pos GE 0 THEN POINT_LUN, tdata_lun, estimated_pos ELSE BEGIN
        PRINT, 'estimated_pos LT 0, file position error!'
        FREE_LUN, tdata_lun
        RETURN, -1L
      ENDELSE
      READU, tdata_lun, time_prev
      time_prev = time_prev  * time_renorm(0)
        POINT_LUN, - tdata_lun, position_before
        position_after = position_before + content_byte_size
        POINT_LUN, tdata_lun, position_after
      READU, tdata_lun, time_act
      time_act = time_act  * time_renorm(0)
      time_diff_act = time_act - time_prev
      IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
      new_estimated_pos = step_byte_size * $
        ROUND((start_time - time_act) / time_diff_act) / ABS(ROUND((start_time - time_act) /$
	time_diff_act)) + estimated_pos
      IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
        new_estimated_pos = new_estimated_pos + step_byte_size
      estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)
      IF estimated_pos LT 0 THEN estimated_pos = 0L
      IF found EQ 0 THEN PRINT, 'new estimate:', estimated_pos, (pos_last_step - step_byte_size)
    ENDWHILE
  ENDIF ;ELSE PRINT, 'number of iterations:', iterations

  FREE_LUN, tdata_lun

  step_pos = ROUND(estimated_pos/(1.0*step_byte_size))

  RETURN, step_pos

END

;##########################################################################

PRO jump_to_energy_step, run, jump_step, energy_lun
; jumps to the desired step of energy file

  COMMON global_vars

  IF jump_step EQ -1 THEN BEGIN
    PRINT, 'omitting run, not in time range'
    RETURN
  ENDIF

  IF N_ELEMENTS(energy_lun) LT 1 THEN energy_lun = 0

  ; get step byte size, jump to position
  energy_step_bytes = 0
  IF (TOTAL(series.request_energy[INDGEN(par.n_energies),*]) GT 0) THEN BEGIN ;n_energies
    time = par.prec_single ? 0.0 : 0.0D

    IF par.in_data NE 'sykx' THEN content = par.prec_single ? $
       FLTARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) : $
       DBLARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) $
    ELSE content = par.prec_single ? $
       FLTARR(par.nkx0/2,par.ny0,par.nz0,/NOZERO) : $
       DBLARR(par.nkx0/2,par.ny0,par.nz0,/NOZERO)

    READU, energy_lun, time
    FOR energy_nr = 0, par.n_energies - 1 DO READU, energy_lun, content ;n_energies
    POINT_LUN, - energy_lun, energy_step_bytes
    jump_byte = LONG64(energy_step_bytes) * jump_step
    POINT_LUN, energy_lun, jump_byte
  ENDIF

END

;##########################################################################

PRO read_energy_step, energy_lun

  COMMON global_vars

  energy_time = par.prec_single ? 0.0 : 0.0D

  request_energy = TOTAL(series.request_energy[INDGEN(par.n_energies),*],2) ;n_energies

  IF par.in_data NE 'sykx' THEN $
    content_byte_size = 4LL * (2 - par.prec_single) * par.nkx0 * (par.nky0 - $
    par.ky0_ind) * par.nz0 + 8 ELSE $
    content_byte_size = 4LL * (2 - par.prec_single) * (par.nkx0/2) * $
    par.ny0 * par.nz0 + 8


  IF energy_lun GT 0 THEN BEGIN
    IF EOF(energy_lun) THEN BEGIN
      FREE_LUN, energy_lun
      energy_lun = 0
    ENDIF ELSE BEGIN
      time = par.prec_single ? 0.0 : 0.0D

      IF par.in_data NE 'sykx' THEN BEGIN
         content = par.prec_single ? $
          FLTARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO) : $
          DBLARR(par.nkx0,par.nky0-par.ky0_ind,par.nz0,/NOZERO)
      ENDIF ELSE BEGIN
         content = par.prec_single ? $
          FLTARR(par.ny0,par.nkx0/2,par.nz0,/NOZERO) : $
          DBLARR(par.ny0,par.nkx0/2,par.nz0,/NOZERO)      
      ENDELSE


      READU, energy_lun, time
      time = time * time_renorm(0)

      FOR energy_nr = 0, par.n_energies - 1 DO BEGIN  ;n_energies
        IF request_energy[energy_nr] GE 1 THEN BEGIN
          IF !QUIET NE 1 THEN BEGIN
            PRINT, '    reading ' + rm0es(content_byte_size/1048576.0) + $
	      ' MB from energy file (energy time: ' + rm0es(time) + ')'
	    stime = SYSTIME(1)
          ENDIF

          READU, energy_lun, content

          IF cfg.kpar_filter GE 0 THEN BEGIN
            IF (cfg.kpar_filter LT par.nz0 / 2) AND (par.nz0 GT 1) THEN BEGIN
              content_kpar = FFT(content,DIMENSION=3)
              IF cfg.kpar_filter EQ 0 THEN BEGIN
                content_kpar[*,*,1:*] = 0
              ENDIF ELSE BEGIN
                neg_kpar = par.nz0 - cfg.kpar_filter
                content_kpar[*,*,0:cfg.kpar_filter-1] = 0
                content_kpar[*,*,cfg.kpar_filter+1:neg_kpar-1] = 0
                IF cfg.kpar_filter NE 1 THEN content_kpar[*,*,neg_kpar+1:*] = 0
              ENDELSE
              content[0,0,0] = FFT(content_kpar,DIMENSION=3,/INVERSE)
            ENDIF
          ENDIF

          IF !QUIET NE 1 THEN $
            PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'

          IF par.in_data EQ 'sxky' THEN BEGIN
            IF par.ky0_ind GT 0 THEN $
              (*energy[energy_nr].sxky)[*,0:par.ky0_ind-1,*] = 0
            (*energy[energy_nr].sxky)[0,par.ky0_ind,0] = $
              TEMPORARY(content*energy_renorm(energy_nr))
          ENDIF ELSE IF par.in_data EQ 'sykx' THEN BEGIN
               FOR ik=0,par.nz0-1 DO BEGIN
                  (*energy[energy_nr].kxsy)[0:par.nkx0/2-1,*,ik] = $
                     TRANSPOSE(content[*,*,ik]*$
                     energy_renorm(energy_nr))
                  (*energy[energy_nr].kxsy)[par.nkx0/2:*,*,ik] = 0.0
               ENDFOR
          ENDIF ELSE BEGIN
            IF par.ky0_ind GT 0 THEN BEGIN
              (*energy[energy_nr].kxky)[*,0:par.ky0_ind-1,*] = 0
              (*energy[energy_nr].kxky)[0,par.ky0_ind,0] = $
                TEMPORARY(content*energy_renorm(energy_nr))
            ENDIF ELSE (*energy[energy_nr].kxky) = $
              TEMPORARY(content*energy_renorm(energy_nr))
          ENDELSE
        ENDIF ELSE BEGIN ; if data not requested, jump to next field
          POINT_LUN, - energy_lun, position_before
          position_after = position_before + content_byte_size
          POINT_LUN, energy_lun, position_after
        ENDELSE
      ENDFOR

      energy_time = time
    ENDELSE
  ENDIF

END

;##########################################################################

PRO energy_loop, diag0

  COMMON global_vars

  ;IF 5 LT 1 THEN BEGIN
  ;  printerror, 'no energy data available'
  ;  RETURN
  ;ENDIF

  IF (TOTAL(series.request_energy) EQ 0) AND NOT series.request_eterms THEN BEGIN
    series.request_energy[*,3]=[1,1,1,1,1,1]
    ;printerror, 'no energy data requested in energy loop'
    ;RETURN
  ENDIF

  IF series.request_eterms THEN get_eterms_data

  request_energy = TOTAL(series.request_energy[INDGEN(par.n_energies),*]) ;n_energies

  IF request_energy GT 0 THEN BEGIN
    file_energy = get_file_path('energy3d')

    IF file_energy.exist EQ 0 THEN BEGIN
      printerror, 'Skipping energy loop due to missing energy file(s)'
      RETURN
    ENDIF

    last_time = -1.0D * 1e20

    IF gui.out.sparsefac GT 1 THEN PRINT, 'Sparse energy loop: reading only every ' + $
      append_th(gui.out.sparsefac) + ' step'

    energy_lun = 0L
    create_energy_struct

    IF ((WHERE(*par.istep_energy3d EQ 0) EQ -1) AND $
      (request_energy GT 0)) OR $
      (gui.out.sparsefac GT 1) THEN BEGIN

      energy_step_byte_size = (4LL * (2 - par.prec_single) * par.nkx0 * $
        (par.nky0 - par.ky0_ind) * par.nz0 + 8) * par.n_energies + $   ;n_energies
        16 - 4 * par.prec_single
    ENDIF

    FOR run = 0, series.n_runs - 1 DO BEGIN
      run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]

      OPENR, energy_lun, file_energy.path[run], /GET_LUN, ERROR=err, $
        /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian
      IF err NE 0 THEN BEGIN
        printerror, file_energy.path[run] + ' does not exist', /popup
	RETURN
      ENDIF ELSE IF !QUIET NE 1 THEN PRINT, 'reading ', $
        file_energy.path[run]

      jump_step = get_energy_time(run)
      IF jump_step LT 0 THEN BEGIN
        end_loop = 1
        IF !QUIET NE 1 THEN PRINT, 'omitting run ', run_label
      ENDIF ELSE BEGIN
        jump_to_energy_step, run, jump_step, energy_lun
        IF jump_step GT 0 THEN $
          PRINT, '- starting energy loop at step ', rm0es(jump_step), ' -'
        IF jump_step EQ 0 THEN PRINT, '- starting energy loop -'
        IF (file_energy.exist GE 1) AND $
          (!QUIET EQ 1) THEN PRINT, 'reading run ', run_label, ': energy'
        end_loop = 0
      ENDELSE

      ; for unstable runs which have follow-ups from secure checkpoints
      ; (overlapping time windows):
      IF run LT series.n_runs - 1 THEN start_next = get_energy_time(run+1,/first) $
        ELSE start_next = 2.0 * ABS(gui.out.end_t)

      sparse_count = 0L

      ; --- start energy loop ---

      WHILE NOT end_loop DO BEGIN
        IF sparse_count MOD gui.out.sparsefac EQ 0 THEN $
          read_energy_step, energy_lun ELSE BEGIN

          IF energy_lun GT 0 THEN BEGIN
            POINT_LUN, - energy_lun, position_before
            position_after = position_before + energy_step_byte_size
            POINT_LUN, energy_lun, position_after
          ENDIF
        ENDELSE

        IF last_time LT energy_time THEN BEGIN
          IF energy_time GT (gui.out.end_t + 1e-5) THEN end_loop = 1 $
            ELSE IF (energy_time GT start_next) THEN BEGIN
            PRINT, 'overlapping time windows; continuing with follow-up'
            end_loop = 1
          ENDIF ELSE IF energy_time GE gui.out.start_t - 1e-5 THEN BEGIN
            IF sparse_count MOD gui.out.sparsefac EQ 0 THEN BEGIN
              sparse_count = 0L

              last_time = energy_time
              series.step_count += 1
              PRINT, '  time: ', rm0es(energy_time)

              calc_energy_ffts
              call_diags, diag0, 'loop'
            ENDIF

            ; break loop by clicking and holding left in plot window
            CURSOR, xtemp, ytemp, /NOWAIT
            IF !MOUSE.button EQ 1 THEN BEGIN
              IF energy_lun GT 0 THEN FREE_LUN, energy_lun

              destroy_energy_struct

              PRINT, 'loop aborted by user'
              RETURN
            ENDIF ; mouse click
          ENDIF
        ENDIF ; last time LT data time

        end_loop = end_loop OR (request_energy AND (energy_lun EQ 0))

        sparse_count += 1
      ENDWHILE ; current energy loop

      ; close current energy files
      IF request_energy GT 0 AND energy_lun GT 0 THEN FREE_LUN, energy_lun
    ENDFOR

    destroy_energy_struct

    IF jump_step GE 0 THEN PRINT, '- finished energy loop -'
  ; make sure the output part is called if only eterms are requested
  ENDIF ELSE IF series.request_eterms THEN series.step_count = 1L

END

;##########################################################################

PRO read_eterms, run, time=time, data=data, $
  get_steps=get_steps, n_steps=n_steps, n_en=n_en
; Reads the energy_[run] file and makes the corresponding data available
; as a time array and etermdat[fiels,time].
; Alternatively, set /get_steps to return the step number in n_steps.

  COMMON global_vars

  n_steps = 0

  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
  eterms_file = get_file_path('energy',run)

  exist = FILE_TEST(eterms_file.path,/READ,/REGULAR)
  IF NOT exist THEN BEGIN
    printerror, eterms_file.path + ' does not exist, skipping eterms'
    RETURN
  ENDIF

  n_lines = FILE_LINES(eterms_file.path)
  IF n_lines LE 14 THEN BEGIN
    printerror, eterms_file.path + ' has wrong format, skipping eterms'
    RETURN
  ENDIF

  OPENR, eterms_lun, eterms_file.path, /GET_LUN, ERROR=err
  IF err NE 0 THEN BEGIN
    printerror, eterms_file.path + ' could not be read, skipping eterms'
    RETURN
  ENDIF

  header = STRARR(14)
  READF, eterms_lun, header
  n_en = STRMID(header[13],0,6) EQ ' #  14' ? 14 : 13

  IF n_en EQ 14 THEN BEGIN
    line = ''
    READF, eterms_lun, line
  ENDIF

  n_steps = n_lines - n_en - 1
  IF KEYWORD_SET(get_steps) THEN BEGIN
    FREE_LUN, eterms_lun
    RETURN
  ENDIF

  rawdata = FLTARR(n_en,n_steps,/NOZERO)
  READF, eterms_lun, rawdata

  FREE_LUN, eterms_lun

  time = REFORM(rawdata[0,*])
  data = TRANSPOSE(rawdata[1:*,*])

END

;##########################################################################

PRO get_eterms_data

  COMMON global_vars

  n_steps_all = LONARR(series.n_runs)
  ;print,"series.n_runs",series.n_runs
  ;print,"n_steps_all",n_steps_all
  FOR run = 0, series.n_runs - 1 DO BEGIN
    read_eterms, run, /get_steps, n_steps=n_steps, n_en=n_en
    n_steps_all[run] = n_steps
  ENDFOR

  n_steps_tot = TOTAL(n_steps_all)
  IF n_steps_tot LT 1 THEN RETURN
  start_step = series.n_runs GT 1 ? [0L,n_steps_all[0:series.n_runs-2]] : 0L

  eterms_time = PTR_NEW(FLTARR(n_steps_tot,/NOZERO),/NO_COPY)
  eterms = PTR_NEW(FLTARR(n_steps_tot,n_en,/NOZERO),/NO_COPY)

  FOR run = 0, series.n_runs - 1 DO BEGIN
    read_eterms, run, time=time, data=data

    (*eterms_time)[start_step[run]] = time
    (*eterms)[start_step[run],0] = data
  ENDFOR

END
