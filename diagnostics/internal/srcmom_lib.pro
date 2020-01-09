;##########################################################################
;# this library contains almost all mom/field related internal functions  #
;##########################################################################

FUNCTION srcmom_renorm, var, Lref=Lref

  COMMON global_vars

  IF NOT KEYWORD_SET(Lref) THEN Lref = series.Lref

  IF var GT par.n_fields + 5 THEN var = ((var-par.n_fields) MOD 6) + $
    par.n_fields

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref
  c_ref = SQRT(series.Tref*series.Qref/series.mref)
  p_ref = series.nref * series.Tref * series.Qref

  CASE var OF
    2    : RETURN, 1.0
    1    : RETURN, 1.0
    0    : RETURN, 1.0 
    ELSE : printerror, 'error in srcmom_renorm (var ' + rm0es(var) + ' not found)'
  ENDCASE


END

;#############################################################################

PRO create_srcmom_struct

  COMMON global_vars
  
  IF TOTAL(series.request_srcmom) EQ 0 THEN RETURN


  IF !QUIET NE 1 THEN BEGIN
    mem_use = 8L * par.nx0 / 1048576.0 * (0L + $
      TOTAL(series.request_srcmom)) * $
      N_ELEMENTS(*gui.out.spec_select)
    PRINT, '    creating data structures (required memory: ' + $
      rm0es(mem_use) + ' MB)'
    PRINT, '    currently in use: ' + mem_usage()
    stime = SYSTIME(1)
  ENDIF

  srcmom = REPLICATE(PTR_NEW(),$
    N_ELEMENTS(*gui.out.spec_select),par.n_srcmoms)

  ind_arr = WHERE(series.request_srcmom[*,0] GT 0)
  FOR isp = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO BEGIN
    IF ind_arr[0] NE -1 THEN FOR ivar = 0, N_ELEMENTS(ind_arr) - 1 DO $
      srcmom[isp,ind_arr[ivar]] = par.prec_single ? $
      PTR_NEW(FLTARR(par.nx0,3,/NOZERO),/NO_COPY) : $
      PTR_NEW(DBLARR(par.nx0,3,/NOZERO),/NO_COPY)
  ENDFOR

  IF !QUIET NE 1 THEN BEGIN
    PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'
    PRINT, '    current memory usage: ' + mem_usage()
  ENDIF

END

;##########################################################################

PRO destroy_srcmom_struct

  COMMON global_vars

  PTR_FREE, srcmom

END

;##########################################################################

FUNCTION get_srcmom_time, run, get_n_steps=get_n_steps, first=first, last=last,$
                       step_time=step_time

  COMMON global_vars

  start_time = gui.out.start_t

  run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
  file_srcmom = get_file_path('srcmom_'+spec[*].name,run)

  IF TOTAL(file_srcmom.exist) EQ 0 THEN RETURN, 0.0
  sp = (WHERE(file_srcmom[*].exist))[0]

  request_srcmom = TOTAL(series.request_srcmom)  

  OPENR, tdata_lun, file_srcmom[sp].path, /GET_LUN, ERROR=err, $
      /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian

  content_byte_size = LONG64(4L*(2-par.prec_single)*par.nx0*3L+8L) * 3L
  full_entry_size = content_byte_size + 16L - 4L * par.prec_single ; data block + time entry

  interleave_factor = 1L
  step_byte_size = full_entry_size*interleave_factor

  time = par.prec_single ? 0.0 : 0.0D
  time_first_step = time
  time_second_step = time
  time_last_step = time

  n_steps = ROUND((FSTAT(tdata_lun)).SIZE/FLOAT(step_byte_size))
  pos_last_step = (n_steps-1)*step_byte_size

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
  step_time = time_last_step
  
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

  WHILE (found EQ 0) AND (iterations LE max_iterations) AND $
     (estimated_pos GE 0) DO BEGIN
    IF estimated_pos GE 0 THEN POINT_LUN, tdata_lun, estimated_pos ELSE BEGIN
      PRINT, 'estimated_pos LT 0, file position error!'
      FREE_LUN, tdata_lun
      RETURN, -1L
    ENDELSE
    time_prev = par.prec_single ? 0.0 : 0.0D
    READU, tdata_lun, time_prev
    time_prev = time_prev  * time_renorm(0)
    POINT_LUN, - tdata_lun, position_before
    position_after = position_before + content_byte_size + $
       full_entry_size*(interleave_factor-1L)
    POINT_LUN, tdata_lun, position_after
    time_act = par.prec_single ? 0.0 : 0.0D
    READU, tdata_lun, time_act
print, 'time_act = ', time_act
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
      time_prev *= time_renorm(0)
      POINT_LUN, - tdata_lun, position_before
      position_after = position_before + content_byte_size + $
        full_entry_size*(interleave_factor-1L)
      POINT_LUN, tdata_lun, position_after
      READU, tdata_lun, time_act
      time_act *= time_renorm(0)
      time_diff_act = time_act - time_prev
      IF (start_time GT time_prev) AND (start_time LE time_act) THEN found = 1
      new_estimated_pos = step_byte_size * $
        ROUND((start_time-time_act)/time_diff_act) / ABS(ROUND((start_time-time_act)/$
	time_diff_act)) + estimated_pos
      IF (new_estimated_pos EQ estimated_pos) AND (found EQ 0) THEN $
        new_estimated_pos = new_estimated_pos + step_byte_size
      estimated_pos = new_estimated_pos < (pos_last_step - step_byte_size)
      IF estimated_pos LT 0 THEN estimated_pos = 0L
;      IF found EQ 0 THEN PRINT, 'new estimate:', estimated_pos, (pos_last_step - step_byte_size)
    ENDWHILE
  ENDIF ;ELSE PRINT, 'number of iterations:', iterations

  step_pos = ROUND(estimated_pos/(1.0*step_byte_size))
  estimated_pos = step_pos*step_byte_size

  POINT_LUN, tdata_lun, estimated_pos
  step_time = par.prec_single ? 0.0 : 0.0D
  READU, tdata_lun, step_time
  step_time *= time_renorm(0)

  ; now set the file pointer exactly to the first step in the time window
  WHILE ((step_time LT start_time) AND NOT EOF(tdata_lun)) DO BEGIN
     estimated_pos += step_byte_size
     POINT_LUN, tdata_lun, estimated_pos
     READU, tdata_lun, step_time
     step_time *= time_renorm(0)  
  ENDWHILE
  
  step_pos = estimated_pos / step_byte_size
 
  FREE_LUN, tdata_lun  

  RETURN, step_pos

END

;##########################################################################

PRO jump_to_srcmom_step, run, jump_step, srcmom_luns
; jumps to the desired step of mom and/or field file

  COMMON global_vars

  IF jump_step EQ -1 THEN BEGIN
    PRINT, 'omitting run, not in time range'
    RETURN
  ENDIF

  interleave_factor_srcmom = 1L

  ; get step byte size, jump to position
  FOR sp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     IF (TOTAL(series.request_srcmom) GT 0) THEN BEGIN
      time = par.prec_single ? 0.0 : 0.0D
      content = par.prec_single ? $
        FLTARR(par.nx0,3,/NOZERO) : $
      	DBLARR(par.nx0,3,/NOZERO)
      READU, srcmom_luns[sp], time
      FOR mom_nr = 0, par.n_srcmoms - 1 DO READU, srcmom_luns[sp], content
      POINT_LUN, - srcmom_luns[sp], srcmom_step_bytes
      jump_byte = LONG64(srcmom_step_bytes) * jump_step * $
                  interleave_factor_srcmom
      POINT_LUN, srcmom_luns[sp], jump_byte
    ENDIF
  ENDFOR

END

;##########################################################################

PRO read_srcmom_step, srcmom_luns

  COMMON global_vars

  srcmom_time = par.prec_single ? 0.0 : 0.0D

  request_srcmom = TOTAL(series.request_srcmom)

  content_byte_size = 4LL * (2 - par.prec_single) * par.nkx0 * 3L + 8L

  FOR j = 0, gui.out.n_spec_sel - 1 DO BEGIN
    IF srcmom_luns[j] GT 0 THEN BEGIN
      IF EOF(srcmom_luns[j]) THEN BEGIN
        FREE_LUN, srcmom_luns[j]
	srcmom_luns[j] = 0
      ENDIF ELSE BEGIN
        time = par.prec_single ? 0.0 : 0.0D
        content = par.prec_single ? $
          FLTARR(par.nkx0,3,/NOZERO) : $
	  DBLARR(par.nkx0,3,/NOZERO)

        READU, srcmom_luns[j], time
        time = time * time_renorm(0)

        FOR mom_nr = 0, par.n_srcmoms - 1 DO BEGIN
          IF series.request_srcmom[mom_nr] GE 1 THEN BEGIN
            IF !QUIET NE 1 THEN BEGIN
              PRINT, '    reading ' + rm0es(content_byte_size/1048576.0) + $
                ' MB from srcmom file (srcmom time: ' + rm0es(time) + ')'
                stime = SYSTIME(1)
            ENDIF
            
            READU, srcmom_luns[j], content
            (*srcmom[j,mom_nr]) = $
	      	TEMPORARY(content*srcmom_renorm(mom_nr))

            IF !QUIET NE 1 THEN PRINT, '    elapsed time: ' + $
              rm0es(SYSTIME(1)-stime) + ' s'
          ENDIF ELSE BEGIN ; if data not requested, jump to next field
            POINT_LUN, - srcmom_luns[j], position_before
            position_after = position_before + content_byte_size
            POINT_LUN, srcmom_luns[j], position_after
          ENDELSE
        ENDFOR
        srcmom_time = time
      ENDELSE
    ENDIF
  ENDFOR

END

;##########################################################################

PRO srcmom_loop, diag0

  COMMON global_vars

  IF TOTAL(series.request_srcmom) EQ 0 THEN BEGIN
    printerror, 'no srcmom data requested in srcmom loop'
    RETURN
  ENDIF

  last_time = -1.0D * 1e20

  file_srcmom = get_file_path('srcmom_'+spec[(*gui.out.spec_select)[*]].name)

  request_srcmom = TOTAL(series.request_srcmom)

  IF request_srcmom AND (TOTAL(file_srcmom.exist) NE $
                         gui.out.n_spec_sel) THEN BEGIN
    printerror, 'Skipping srcmom loop due to missing srcmom file(s)'
    RETURN
  ENDIF

  IF gui.out.sparsefac GT 1 THEN PRINT, $
     'Sparse srcmom loop: reading only every ' + $
     append_th(gui.out.sparsefac) + ' step'

  srcmom_luns = LONARR(par.n_spec)
  create_srcmom_struct

  srcmom_content_byte_size = (4LL * (2 - par.prec_single) * par.nx0 * 3L + $
    8L) * 3L + 16 - 4 * par.prec_single

  FOR run = 0, series.n_runs - 1 DO BEGIN
     run_label = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
     
     interleave_srcmom_factor = LONG64(gui.out.sparsefac)
     
     IF request_srcmom GT 0 THEN BEGIN
        FOR j = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO BEGIN
           OPENR, temp_lun, file_srcmom[j].path[run], /GET_LUN, ERROR=err, $
                  /F77_UNFORMATTED, SWAP_ENDIAN=series.swap_endian
           IF err NE 0 THEN BEGIN
              printerror, file_srcmom[j].path[run] + ' does not exist', /popup
              RETURN
           ENDIF ELSE BEGIN
              IF !QUIET NE 1 THEN PRINT, 'reading ', file_srcmom[j].path[run]
              srcmom_luns[j] = temp_lun
           ENDELSE
        ENDFOR
     ENDIF
     
     jump_step = get_srcmom_time(run,step_time=step_time) ;determine start step
     IF jump_step LT 0 THEN BEGIN
        end_loop = 1
        IF !QUIET NE 1 THEN PRINT, 'omitting run ', run_label
     ENDIF ELSE BEGIN
        jump_to_srcmom_step, run, jump_step, srcmom_luns
        IF jump_step GT 0 THEN $
           PRINT, '- starting srcmom loop at step ', rm0es(jump_step), ' -'
        IF jump_step EQ 0 THEN PRINT, '- starting srcmom loop -'
        IF (TOTAL(file_srcmom.exist) GE 1) AND $
           (!QUIET EQ 1) THEN PRINT, 'reading run ', run_label, ': srcmom'
        end_loop = 0
     ENDELSE
     
     ; for unstable runs which have follow-ups from secure checkpoints
     ; (overlapping time windows):
     IF run LT series.n_runs - 1 THEN start_next = get_srcmom_time(run+1,/first) $
     ELSE start_next = 2.0 * ABS(gui.out.end_t)

     IF (NOT par.x_local) THEN read_profiles, run=run
     
     ; --- start mom loop ---

     WHILE NOT end_loop DO BEGIN
        read_srcmom_step, srcmom_luns 

        ;set file pointers to next entry
        FOR j = 0, gui.out.n_spec_sel - 1 DO BEGIN
           IF srcmom_luns[j] GT 0 THEN BEGIN
              POINT_LUN, - srcmom_luns[j], position_before
              position_after = position_before + srcmom_content_byte_size*$
                               LONG64(interleave_srcmom_factor-1)
              POINT_LUN, srcmom_luns[j], position_after
           ENDIF
        ENDFOR
        
        IF last_time LT srcmom_time THEN BEGIN
           IF srcmom_time GT (gui.out.end_t * (1.0D + 1e-7)) THEN end_loop = 1 $
           ELSE IF (srcmom_time GT start_next) THEN BEGIN
              PRINT, 'overlapping time windows; continuing with follow-up'
              end_loop = 1
           ENDIF ELSE IF srcmom_time GE gui.out.start_t * (1.0D - 1e-7) THEN BEGIN
              last_time = srcmom_time
              series.step_count += 1
              PRINT, '  time: ', rm0es(srcmom_time)
              
              call_diags, diag0, 'loop'
              
              ; break loop by clicking and holding left in plot window
              IF gui.window.main NE 0 THEN BEGIN ; omit for command line version
                 CURSOR, xtemp, ytemp, /NOWAIT
                 IF !MOUSE.button EQ 1 THEN BEGIN
                    FOR j = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO $
                       IF srcmom_luns[j] GT 0 THEN FREE_LUN, srcmom_luns[j]

                    destroy_srcmom_struct
                    
                    PRINT, 'loop aborted by user'
                    RETURN
                 ENDIF  ; mouse click
              ENDIF ; gui.window.main
           ENDIF ;srcmom_time GT start?
        ENDIF ; last time LT data time
        
        end_loop = end_loop OR (request_srcmom AND (TOTAL(srcmom_luns) EQ 0))
     ENDWHILE                   ; current srcmom loop
    
                                ; close current mom and field files
     IF request_srcmom GT 0 THEN BEGIN
        FOR j = 0, N_ELEMENTS(*gui.out.spec_select) - 1 DO $
           IF srcmom_luns[j] GT 0 THEN FREE_LUN, srcmom_luns[j]
     ENDIF
  ENDFOR
  
  destroy_srcmom_struct
  
  IF jump_step GE 0 THEN PRINT, '- finished srcmom loop -'

END
