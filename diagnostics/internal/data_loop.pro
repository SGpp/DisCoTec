PRO call_diags, diag0, stage
; stage can be 'init', 'loop', or 'output'

  COMMON global_vars

  ; custom diagnostics
  cdiag = cdiag0
  WHILE PTR_VALID(cdiag) DO BEGIN
    IF (*cdiag).selected EQ 1 THEN BEGIN
      IF !QUIET NE 1 THEN BEGIN
        PRINT, '    starting ',(*cdiag).name, '_' + stage
        stime = SYSTIME(1)
      ENDIF
      CALL_PROCEDURE, (*cdiag).name + '_' + stage, cdiag
      IF !QUIET NE 1 THEN $
        PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'
    ENDIF
    cdiag = (*cdiag).next
  ENDWHILE

  ; regular diagnostics
  diag = diag0
  WHILE PTR_VALID(diag) DO BEGIN
    IF (*diag).selected THEN BEGIN
      IF stage EQ 'init' THEN set_external_vars, diag
      diag_call = (*diag).name + '_' + stage
      IF !QUIET NE 1 THEN BEGIN
         PRINT, '    starting ', (*diag).name, '_' + stage
         stime = SYSTIME(1)
      ENDIF
      CALL_PROCEDURE, diag_call, diag
      IF !QUIET NE 1 THEN $
        PRINT, '    elapsed time: ' + rm0es(SYSTIME(1)-stime) + ' s'
      IF stage EQ 'output' THEN set_external_vars, diag, /destroy
    ENDIF
    diag = (*diag).next
  ENDWHILE

END

;##########################################################################

PRO data_loop, tab=tab

  COMMON global_vars

  IF NOT vm_sav THEN recompile_cdiags

  series.request_mom[*,*] = 0 ; reset requesting of fields and fft
  series.request_nrg = 0
  series.request_energy[*,*] = 0
  series.request_srcmom = 0

  IF NOT file_par.exist THEN BEGIN
    printerror, 'skipping time loop (missing parameter file)'
    RETURN
  ENDIF ELSE BEGIN
    read_par
    series.Lref = set_Lref(gui.out.Lref_sel)
    read_geometry
    set_series_lengths
    IF NOT par.x_local THEN read_profiles
  ENDELSE

  ; init stage of diagnostics
  diag0 = (*gui.tab.arr)[tab].diag0
  gui.misc.n_sel_diags = get_sel_diags(diag0)
  call_diags, diag0, 'init'

  series.step_count = 0L

  CALL_PROCEDURE, (STRSPLIT((*gui.tab.arr)[tab].name,' ',/EXTRACT))[0] + $
    '_loop', diag0

  IF !QUIET NE 1 THEN PRINT, '- number of steps: ', $
    rm0es(series.step_count), ' -'

  ; output stage of diagnostics
  IF series.step_count GE 1 THEN $
    call_diags, diag0, 'output'

  ; free pointers
  PTR_FREE, series.geom, series.kx, series.ky, series.lxarr
  destroy_energy_struct, /free_eterms
  free_internal_vars, diag0, cdiag0

END
