;#############################################################################
;#    this library contains almost all nlt related internal functions     #
;#############################################################################

FUNCTION nlt_renorm, var, Lref=Lref
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

    RETURN, E_ref * c_ref / series.Lref ; units of E/t

END

;#############################################################################

PRO create_nlt_struct

  COMMON global_vars
  

  IF !QUIET NE 1 THEN BEGIN
    mem_use = 16L * (2*par.nky0-1)*par.nkx0 / 1048576.0 
    PRINT, '  creating data structures (required memory: ' + $
      rm0es(mem_use) + ' MB)'
    PRINT, '    currently in use: ' + mem_usage()
    stime = SYSTIME(1)
  ENDIF

  num_nlt_use = par.num_nlt_modes
  IF par.nlp_gdt THEN num_nlt_use = par.num_nlt_pod_modes+1
  IF par.prec_single THEN BEGIN
    nlt2D = PTR_NEW(FLTARR(par.nkx0,2*par.nky0-1,num_nlt_use))
  ENDIF ELSE BEGIN
    nlt2D = PTR_NEW(DBLARR(par.nkx0,2*par.nky0-1,num_nlt_use))
  ENDELSE 
  
  IF !QUIET NE 1 THEN BEGIN
    PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'
    PRINT, '    current memory usage: ' + mem_usage()
  ENDIF

END

;##########################################################################

PRO destroy_nlt_struct, free_nltterms=free_nltterms

  COMMON global_vars

  IF KEYWORD_SET(free_nltterms) THEN BEGIN
    IF PTR_VALID(nltterms) THEN PTR_FREE, nltterms
    IF PTR_VALID(nltterms_time) THEN PTR_FREE, nltterms_time
    RETURN
  ENDIF

  PTR_FREE, nlt2d

END

;##########################################################################

FUNCTION get_nlt_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  start_time = gui.out.start_t

  IF KEYWORD_SET(last) THEN BEGIN
    istep_nlt = - (*par.istep_nlt)[run]
  ENDIF ELSE BEGIN ; for /last, use finer scale file, else use coarser
    istep_nlt = (*par.istep_nlt)[run]
  ENDELSE

  IF par.nlp_gdt THEN BEGIN
    nlp_file_string='nlp_ky'+STRING(par.nlp_kyind,FORMAT='(i4.4)')+$
     'kx'+STRING(par.nlp_kxind,FORMAT='(i4.4)')
    content_byte_size = (4LL * (2 - par.prec_single) * $
      par.nkx0 * (2*par.nky0 - 1) + 8) * (par.num_nlt_pod_modes+1)
  ENDIF ELSE BEGIN
    nlp_file_string='nlt'    
    content_byte_size = (4LL * (2 - par.prec_single) * $
      par.nkx0 * (2*par.nky0 - 1) + 8) * par.num_nlt_modes
  ENDELSE
  print,"Using file:",nlp_file_string
  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
  file_nlt = get_file_path(nlp_file_string,run)

  IF NOT file_nlt.exist THEN RETURN, 0.0

  OPENR, tdata_lun, file_nlt.path, /GET_LUN, ERROR=err, /F77_UNFORMATTED, $
    SWAP_ENDIAN=series.swap_endian

  IF err NE 0 THEN BEGIN
    printerror, 'Please notify the developers about this error'
    RETURN, 0L
  ENDIF

  ;IF par.nlt_old THEN BEGIN
  ;  content_byte_size = (4LL * (2 - par.prec_single) * $
  ;    par.nkx0 * (2*par.nky0 - 1) + 8) * 2*par.num_nlt_modes
  ;ENDIF
  step_byte_size = content_byte_size + 16 - 4 * par.prec_single ; data block + time entry

  time = par.prec_single ? 0.0 : 0.0D
  time_first_step = time
  time_second_step = time
  time_last_step = time

  ;pos_last_step = (FSTAT(tdata_lun)).SIZE - step_byte_size
  n_steps = ROUND((FSTAT(tdata_lun)).SIZE/FLOAT(step_byte_size))
  n_steps = n_steps - 2 ;otherwise it doesn't work . . . ?
  pos_last_step = (n_steps - 1) * step_byte_size
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

PRO jump_to_nlt_step, run, jump_step, nlt_lun
; jumps to the desired step of nlt file

  COMMON global_vars

  IF jump_step EQ -1 THEN BEGIN
    PRINT, 'omitting run, not in time range'
    RETURN
  ENDIF

  num_nlt_read=par.num_nlt_modes
  IF par.nlp_gdt THEN num_nlt_read = par.num_nlt_pod_modes+1

  IF N_ELEMENTS(nlt_lun) LT 1 THEN nlt_lun = 0

  ; get step byte size, jump to position
  nlt_step_bytes = 0
  ;IF (TOTAL(series.request_nlt[INDGEN(par.n_energies),*]) GT 0) THEN BEGIN
    time = par.prec_single ? 0.0 : 0.0D
    content = par.prec_single ? $
      FLTARR(par.nkx0,2*par.nky0-1,/NOZERO) : $
      DBLARR(par.nkx0,2*par.nky0-1,/NOZERO)
    READU, nlt_lun, time
    print,'time',time
    FOR nlt_nr = 0, num_nlt_read - 1 DO READU, nlt_lun, content
    POINT_LUN, - nlt_lun, nlt_step_bytes
    jump_byte = LONG64(nlt_step_bytes) * jump_step
    POINT_LUN, nlt_lun, jump_byte
  ;ENDIF
END

;##########################################################################

PRO read_nlt_step, nlt_lun

  COMMON global_vars


  nlt_time = par.prec_single ? 0.0 : 0.0D

  ;request_nlt = TOTAL(series.request_nlt[INDGEN(par.n_energies),*],2)

  content_byte_size = 4LL * (2 - par.prec_single) * par.nkx0 * (2*par.nky0-1) + 8

  IF nlt_lun GT 0 THEN BEGIN
    IF EOF(nlt_lun) THEN BEGIN
      FREE_LUN, nlt_lun
      nlt_lun = 0
    ENDIF ELSE BEGIN
      time = par.prec_single ? 0.0 : 0.0D
      content = par.prec_single ? $
        FLTARR(par.nkx0,2*par.nky0-1,/NOZERO) : $
        DBLARR(par.nkx0,2*par.nky0-1,/NOZERO)
      content2 = par.prec_single ? $
        FLTARR(par.nkx0,2*par.nky0-1,/NOZERO) : $
        DBLARR(par.nkx0,2*par.nky0-1,/NOZERO)

      READU, nlt_lun, time
      time = time * time_renorm(0)

      num_nlt_read = par.num_nlt_modes
      IF par.nlp_gdt THEN num_nlt_read = par.num_nlt_pod_modes+1

      FOR nlt_nr = 0, num_nlt_read  - 1 DO BEGIN
       ;IF request_nlt[energy_nr] GE 1 THEN BEGIN
          IF !QUIET NE 1 THEN BEGIN
            PRINT, '    reading ' + rm0es(content_byte_size/1048576.0) + $
	      ' MB from nlt file (nlt time: ' + rm0es(time) + ')'
	    stime = SYSTIME(1)
          ENDIF

          READU, nlt_lun, content
          ;IF par.nlt_old THEN BEGIN
          ;  ;print,"Reading second part of nlt.",par.nlt_old
          ;  READU, nlt_lun, content2
          ;  content+=content2
          ;ENDIF
          (*nlt2d)[par.nkx0/2-1:par.nkx0-1,*,nlt_nr] = content[0:par.nkx0/2,*]*nlt_renorm(nlt_nr)
          (*nlt2d)[0:par.nkx0/2-2,*,nlt_nr] = content[par.nkx0/2+1:par.nkx0-1,*]*nlt_renorm(nlt_nr)

          IF !QUIET NE 1 THEN $
            PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'


      ENDFOR

      nlt_time = time
    ENDELSE
  ENDIF

END

;##########################################################################

PRO nlt_loop, diag0

  COMMON global_vars
  IF par.num_nlt_modes EQ 0 THEN BEGIN
    printerror, 'no nlt data available'
    return
  ENDIF

  get_nltterms_data
  series.step_count=1

  IF par.nlp_gdt THEN BEGIN
    nlp_file_string='nlp_ky'+STRING(par.nlp_kyind,FORMAT='(i4.4)')+$
     'kx'+STRING(par.nlp_kxind,FORMAT='(i4.4)')
  ENDIF ELSE BEGIN
    nlp_file_string='nlt'    
  ENDELSE

    file_nlt = get_file_path(nlp_file_string)
    if file_nlt.exist eq 0 then begin
      printerror, 'skipping nlt loop due to missing nlt file(s)'
      return
    endif

    last_time = -1.0d * 1e20

    if gui.out.sparsefac gt 1 then print, 'sparse nlt loop: reading only every ' + $
      append_th(gui.out.sparsefac) + ' step'

    nlt_lun = 0l
    create_nlt_struct

    if ((where(*par.istep_nlt eq 0) eq -1) ) or $
      (gui.out.sparsefac gt 1) then begin

      nlt_step_byte_size = (4ll * (2 - par.prec_single) * par.nkx0 * $
        (2*par.nky0 - 1) + 8) * par.num_nlt_modes + $
        16 - 4 * par.prec_single
    endif

    for run = 0, series.n_runs - 1 do begin
      run_label = (strsplit(series.run_labels,',',/extract))[run]

      openr, nlt_lun, file_nlt.path[run], /get_lun, error=err, $
        /f77_unformatted, swap_endian=series.swap_endian
      if err ne 0 then begin
        printerror, file_nlt.path[run] + ' does not exist', /popup
	return
      endif else if !quiet ne 1 then print, 'reading ', $
        file_nlt.path[run]

      jump_step = get_nlt_time(run)
      if jump_step lt 0 then begin
        end_loop = 1
        if !quiet ne 1 then print, 'omitting run ', run_label
      endif else begin
        jump_to_nlt_step, run, jump_step, nlt_lun
        if jump_step gt 0 then $
          print, '- starting nlt loop at step ', rm0es(jump_step), ' -'
        if jump_step eq 0 then print, '- starting nlt loop -'
        if (file_nlt.exist ge 1) and $
          (!quiet eq 1) then print, 'reading run ', run_label, ': nlt'
        end_loop = 0
      endelse

      ; for unstable runs which have follow-ups from secure checkpoints
      ; (overlapping time windows):
      if run lt series.n_runs - 1 then start_next = get_nlt_time(run+1,/first) $
        else start_next = 2.0 * abs(gui.out.end_t)

      sparse_count = 0l

      ; --- start nlt loop ---

      while not end_loop do begin
        if sparse_count mod gui.out.sparsefac eq 0 then $
          read_nlt_step, nlt_lun else begin

          if nlt_lun gt 0 then begin
            point_lun, - nlt_lun, position_before
            position_after = position_before + nlt_step_byte_size
            point_lun, nlt_lun, position_after
          endif
        endelse

        if last_time lt nlt_time then begin
          if nlt_time gt (gui.out.end_t + 1e-5) then end_loop = 1 $
            else if (nlt_time gt start_next) then begin
            print, 'overlapping time windows; continuing with follow-up'
            end_loop = 1
          endif else if nlt_time ge gui.out.start_t - 1e-5 then begin
            if sparse_count mod gui.out.sparsefac eq 0 then begin
              sparse_count = 0l

              last_time = nlt_time
              series.step_count += 1
              print, '  time: ', rm0es(nlt_time)

              call_diags, diag0, 'loop'
            endif

            ; break loop by clicking and holding left in plot window
            cursor, xtemp, ytemp, /nowait
            if !mouse.button eq 1 then begin
              if nlt_lun gt 0 then free_lun, nlt_lun

              destroy_nlt_struct

              print, 'loop aborted by user'
              return
            endif ; mouse click
          endif
        endif ; last time lt data time

        end_loop = end_loop or ( nlt_lun eq 0 )

        sparse_count += 1
      endwhile ; current nlt loop

      ; close current nlt files
      if nlt_lun gt 0 then free_lun, nlt_lun
    endfor

    destroy_nlt_struct

    if jump_step ge 0 then print, '- finished nlt loop (nlt_loop) -'
  ; make sure the output part is called if only nltterms are requested

end

;##########################################################################

pro read_nltterms, run, time=time, data=data, $
  get_steps=get_steps, n_steps=n_steps
; Reads the nlt_info_[run] file and makes the corresponding data available
; as a time array and etermdat[fiels,time].
; Alternatively, set /get_steps to return the step number in n_steps.

  COMMON global_vars
  IF par.nlp_gdt THEN BEGIN
    num_nltterms=4*(par.num_nlt_pod_modes+1)
  ENDIF ELSE BEGIN
    num_nltterms=240
  ENDELSE
  IF par.nlt_old THEN num_nltterms=80
  print,"num_nltterms",num_nltterms
 
  n_steps = 0

  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]

  IF par.nlp_gdt THEN BEGIN
    nlp_file_string='nlp_info_ky'+STRING(par.nlp_kyind,FORMAT='(i4.4)')+$
     'kx'+STRING(par.nlp_kxind,FORMAT='(i4.4)')
  ENDIF ELSE BEGIN
    nlp_file_string='nlt_info'    
  ENDELSE
  nltterms_file = get_file_path(nlp_file_string,run)

  exist = FILE_TEST(nltterms_file.path,/READ,/REGULAR)
  IF NOT exist THEN BEGIN
    printerror, nltterms_file.path + ' does not exist, skipping nltterms'
    RETURN
  ENDIF

  n_lines = FILE_LINES(nltterms_file.path)
  IF n_lines LE 1 THEN BEGIN
    printerror, nltterms_file.path + ' has wrong format, skipping nltterms'
    RETURN
  ENDIF

  OPENR, nltterms_lun, nltterms_file.path, /GET_LUN, ERROR=err
  IF err NE 0 THEN BEGIN
    printerror, nltterms_file.path + ' could not be read, skipping nltterms'
    RETURN
  ENDIF

  header = STRARR(1)
  READF, nltterms_lun, header

  n_steps = n_lines - 1
  IF KEYWORD_SET(get_steps) THEN BEGIN
    FREE_LUN, nltterms_lun
    RETURN
  ENDIF

  rawdata = FLTARR(num_nltterms+1,n_steps,/NOZERO)
  READF, nltterms_lun, rawdata

  FREE_LUN, nltterms_lun

  time = REFORM(rawdata[0,*])
  data = TRANSPOSE(rawdata[1:*,*])

END

;##########################################################################

PRO get_nltterms_data

  COMMON global_vars

  num_nltterms=240
  IF par.nlt_old THEN num_nltterms=80
  print,"num_nltterms",num_nltterms

  n_steps_all = LONARR(series.n_runs)
  ;print,"series.n_runs",series.n_runs
  ;print,"n_steps_all",n_steps_all
  FOR run = 0, series.n_runs - 1 DO BEGIN
    read_nltterms, run, /get_steps, n_steps=n_steps
    n_steps_all[run] = n_steps
  ENDFOR

  n_steps_tot = TOTAL(n_steps_all)
  IF n_steps_tot LT 1 THEN RETURN
  start_step = series.n_runs GT 1 ? [0L,n_steps_all[0:series.n_runs-2]] : 0L

  nltterms_time = PTR_NEW(FLTARR(n_steps_tot,/NOZERO),/NO_COPY)
  nltterms = PTR_NEW(FLTARR(n_steps_tot,num_nltterms,/NOZERO),/NO_COPY)
  ;print,'n_steps_tot',n_steps_tot

  FOR run = 0, series.n_runs - 1 DO BEGIN
    read_nltterms, run, time=time, data=data

    (*nltterms_time)[start_step[run]] = time
    ;print,'run',run
    ;print,'start_step[run]',start_step[run]
    ;help,*nltterms
    ;help,*nltterms_time
    (*nltterms)[start_step[run],0] = data
  ENDFOR

END
