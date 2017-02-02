PRO create_diag_struct
; Reads all available diagnostics and converts
; information to the internal diag structure (diag).

  COMMON global_vars

  ; --- include regular diags ---
  IF vm_sav THEN BEGIN
    HELP, /FUNCTIONS, NAMES='*_info', OUTPUT=diag_files
    IF cfg.vm_lowcase THEN diag_files = STRLOWCASE(diag_files)
  ENDIF ELSE diag_files = FILE_SEARCH('prog/*.pro',/TEST_REGULAR)

  FOR i = 0, N_ELEMENTS(diag_files) - 1 - vm_sav DO BEGIN
    IF vm_sav THEN BEGIN
      name_base = (STRSPLIT(diag_files[i+1],'_',/EXTRACT))[0]
      diag_info = CALL_FUNCTION(STRTRIM(diag_files[i+1]))
    ENDIF ELSE BEGIN
      name_base = (STRSPLIT(diag_files[i],'./\',/EXTRACT))[1]
      diag_info = CALL_FUNCTION(name_base+'_info')
    ENDELSE

    CASE diag_info.type OF
      'mom_lin' : str_type = 'mom lin'
      'mom_nl'  : str_type = 'mom nl'
      'mom_uni' : str_type = ['mom lin','mom nl']
      ELSE      : str_type = diag_info.type
    ENDCASE

; currently, only global code developers can see global diags:
global_tabs = ['mom_global','srcmom']

IF ((TOTAL(diag_info.type EQ global_tabs) EQ 0) OR show_global) THEN BEGIN

    IF diag_info.ext_vars[0,0] EQ '' THEN BEGIN
      table_entry = ''
      external_vars = ''
      n_vars = 0L
    ENDIF ELSE BEGIN
      n_vars = N_ELEMENTS(diag_info.ext_vars[0,*])
      table_entry = STRARR(4,n_vars)
      FOR j = 0, n_vars - 1 DO BEGIN
        table_entry[0:2,j] = diag_info.ext_vars[*,j]
        IF j EQ 0 THEN external_vars = $
          CREATE_STRUCT(diag_info.ext_vars[0,j],PTR_NEW()) $
        ELSE external_vars = $
          CREATE_STRUCT(external_vars, diag_info.ext_vars[0,j],$
          PTR_NEW())
      ENDFOR
    ENDELSE

    diag_struc = {$
      str_type	    : str_type,$
      name          : name_base,$
      title         : diag_info.title,$
      help_text     : diag_info.help_text,$
      table_entry   : table_entry,$
      external_vars : external_vars,$
      internal_vars : PTR_NEW(),$
      selected      : 0,$
      next          : PTR_NEW()}

    IF i EQ 0 THEN BEGIN
      table_arr = {$
        name        : str_type[0],$
        id          : 0,$
        table_id    : 0,$
        n_params    : n_vars,$
        n_diags     : 1,$
        diag0       : PTR_NEW(diag_struc)}
      FOR j = 1, N_ELEMENTS(str_type) - 1 DO BEGIN
        table_struc = {$
          name 	    : str_type[j],$
          id   	    : 0,$
          table_id  : 0,$
          n_params  : n_vars,$
          n_diags   : 1,$
          diag0	    : PTR_NEW(diag_struc)}
        table_arr = [table_arr,table_struc]
      ENDFOR
    ENDIF ELSE BEGIN
      FOR j = 0, N_ELEMENTS(str_type) - 1 DO BEGIN
        index = WHERE(table_arr.name EQ str_type[j])
        IF index EQ -1 THEN BEGIN
          table_struc = {$
            name      : str_type[j],$
            id        : 0,$
            table_id  : 0,$
  	    n_params  : n_vars,$
            n_diags   : 1,$
            diag0     : PTR_NEW(diag_struc)}
          table_arr = [table_arr,table_struc]
        ENDIF ELSE BEGIN

    	  diag_ptr = table_arr[index].diag0

          WHILE PTR_VALID((*diag_ptr).next) DO $
            diag_ptr = (*diag_ptr).next
          (*diag_ptr).next = PTR_NEW(diag_struc)
          IF table_arr[index].n_params LT n_vars THEN $
    	    table_arr[index].n_params = n_vars
          table_arr[index].n_diags += 1
        ENDELSE
      ENDFOR
    ENDELSE
ENDIF ; --- special global diag treatment
  ENDFOR

  gui.tab.arr = PTR_NEW(table_arr)

  show = 0
  IF show THEN BEGIN
   FOR n = 0, N_ELEMENTS(table_arr) - 1 DO BEGIN
     print, 'name = ', table_arr[n].name
     print, 'n_params = ', table_arr[n].n_params
     print, 'n_diags = ', table_arr[n].n_diags
     diag = table_arr[n].diag0
     WHILE PTR_VALID(diag) DO BEGIN
       print, '-', (*diag).name
       diag = (*diag).next
     ENDWHILE
   ENDFOR
  ENDIF

  ; --- include custom diags ---

  IF vm_sav EQ 0 THEN BEGIN
    cdiag_files = $
      FILE_SEARCH('custom/*.pro',/TEST_REGULAR,COUNT=count)

    first_line = ''
    cdiag_count = 0
    GET_LUN, cdiag_lun
    FOR i = 0, count - 1 DO BEGIN
      OPENR, cdiag_lun, cdiag_files[i]
      READF, cdiag_lun, first_line
      CLOSE, cdiag_lun

      init_prog = (STRSPLIT(first_line,',',/EXTRACT))[0]
      uscore_pos = STRSPLIT(init_prog,'_')
      name = STRMID(init_prog,0,uscore_pos[N_ELEMENTS(uscore_pos)-1]-1)
      name = (STRSPLIT(name,' ',/EXTRACT))[1]

      IF name NE '' THEN BEGIN
        new_cdiag = PTR_NEW({$
          cdiag_nr      : cdiag_count,$
          name          : name,$
          selected      : 0,$
          internal_vars : PTR_NEW(),$
          next          : PTR_NEW()})

        IF NOT PTR_VALID(cdiag_struct) THEN BEGIN   ; first diag
          cdiag0 = new_cdiag
          cdiag_struct = new_cdiag
        ENDIF ELSE BEGIN                            ; append to list
          (*cdiag_struct).next = new_cdiag
          cdiag_struct = new_cdiag
        ENDELSE

        cdiag_count += 1
      ENDIF
    ENDFOR
    FREE_LUN, cdiag_lun
  ENDIF ELSE cdiag_count = 0

  gui.misc.n_cdiags = cdiag_count

END

;#########################################################################

PRO reload_diag_var_names, tab

  COMMON global_vars

  diag = (*gui.tab.arr)[tab].diag0
  (*gui.tab.arr)[tab].n_params = 0

  WHILE PTR_VALID(diag) DO BEGIN
    diag_info = CALL_FUNCTION((*diag).name+'_info')
    IF diag_info.ext_vars[0,0] EQ '' THEN BEGIN
      table_entry = ''
      external_vars = ''
    ENDIF ELSE BEGIN
      n_vars = N_ELEMENTS(diag_info.ext_vars[0,*])
      IF (*gui.tab.arr)[tab].n_params LT n_vars THEN $
    	 (*gui.tab.arr)[tab].n_params = n_vars
      table_entry = STRARR(4,n_vars)
      FOR j = 0, n_vars - 1 DO BEGIN
        table_entry[0:2,j] = diag_info.ext_vars[*,j]
	IF j EQ 0 THEN external_vars = $
  	  CREATE_STRUCT(diag_info.ext_vars[0,j],PTR_NEW()) $
	ELSE external_vars = $
	  CREATE_STRUCT(external_vars, diag_info.ext_vars[0,j],$
          PTR_NEW())
      ENDFOR
    ENDELSE

    diag_struc = {$
      str_type      : (*diag).str_type,$
      name          : (*diag).name,$
      title         : diag_info.title,$
      help_text     : diag_info.help_text,$
      table_entry   : table_entry,$
      external_vars : external_vars,$
      internal_vars : (*diag).internal_vars,$
      selected      : (*diag).selected,$
      next          : (*diag).next}

    (*diag) = diag_struc
    diag = (*diag).next
  ENDWHILE

END

;#########################################################################

FUNCTION get_diag_ptr, ptr0, n
; gets i'th pointer of the diag structure

  COMMON global_vars
  diag = ptr0

  i = 0
  WHILE PTR_VALID(diag) AND (i LT n) DO BEGIN
    diag = (*diag).next
    i = i + 1
  ENDWHILE

  IF i EQ n THEN RETURN, diag $
  ELSE RETURN, PTR_NEW()

END

;#########################################################################

FUNCTION get_cdiag_pointer, cdiag_nr
; gets a pointer to the cdiag structure

  COMMON global_vars

  cdiag = cdiag0
  WHILE PTR_VALID(cdiag) DO BEGIN
    IF (*cdiag).cdiag_nr EQ cdiag_nr THEN RETURN, cdiag
    cdiag = (*cdiag).next
  ENDWHILE
  RETURN, PTR_NEW()

END

;#############################################################################

PRO recompile_cdiags
; This procedure makes use of RESOLVE_ROUTINE, which can only recompile
; routines the names of which are identical to the file names in the '.' dir.
; Therefore, custom diags have to have 'PRO [diagname] & END'.

  COMMON global_vars

  CD, 'custom'

  cdiag = cdiag0
  WHILE PTR_VALID(cdiag) DO BEGIN
    IF (*cdiag).selected EQ 1 THEN BEGIN
      PRINT, 'recompiling ', (*cdiag).name
      RESOLVE_ROUTINE, (*cdiag).name, /COMPILE_FULL_FILE
    ENDIF
    cdiag = (*cdiag).next
  ENDWHILE

  CD, '..'

END

;##########################################################################

PRO set_external_vars, diag, destroy=destroy, dummy = dummy
; set external var pointers to table entries
; needs to be done BEFORE calling any diag_init
; set destroy if all external vars pointers shall be freed

  IF (*diag).table_entry[0,0] EQ '' THEN n_vars = 0 $
    ELSE n_vars = N_TAGS((*diag).external_vars)

  FOR v = 0, n_vars - 1 DO BEGIN
    PTR_FREE, (*diag).external_vars.(v)

    IF NOT KEYWORD_SET(destroy) THEN $
      (*diag).external_vars.(v) = $
      convert_table_entry((*diag).table_entry[3,v],/ptr)
  ENDFOR

END

;##########################################################################


FUNCTION set_internal_vars, diag, internal_var_set
; Enables internal variables for user defined diagnostics.
; Input: diag pointer, internal variable set for this diag.
; This function is to be used for custom diags, as well.

  IF PTR_VALID(diag) THEN (*diag).internal_vars = $
    PTR_NEW(internal_var_set) ELSE BEGIN

    PRINT, 'error in set_internal_vars:'
    PRINT, '  diag pointer not valid, not registering internal variables'
  ENDELSE

  int_var_pointer = (*diag).internal_vars
  RETURN, int_var_pointer
END

;##########################################################################

FUNCTION get_sel_diags, diag0

  n_sel_diags = 0

  diag = diag0
  WHILE PTR_VALID(diag) DO BEGIN
     IF (*diag).selected THEN n_sel_diags += 1
     diag = (*diag).next
  ENDWHILE

  RETURN, n_sel_diags
END

;##########################################################################

PRO free_internal_vars, diag0, cdiag0
; frees internal vars only, not diag or var structures

  diag = diag0
  WHILE PTR_VALID(diag) DO BEGIN
    IF PTR_VALID((*diag).internal_vars) THEN BEGIN
      ; -- look for pointers in internal var structure and destroy them
      FOR tag = 0, N_TAGS((*(*diag).internal_vars))-1 DO BEGIN
        IF TOTAL(PTR_VALID((*(*diag).internal_vars).(tag))) NE 0 THEN $
	  PTR_FREE, (*(*diag).internal_vars).(tag)
      ENDFOR
      PTR_FREE, (*diag).internal_vars
    ENDIF
    diag = (*diag).next
  ENDWHILE

  cdiag = cdiag0
  WHILE PTR_VALID(cdiag) DO BEGIN
    IF PTR_VALID((*cdiag).internal_vars) THEN $
      PTR_FREE, (*cdiag).internal_vars
    cdiag = (*cdiag).next
  ENDWHILE

END

;##########################################################################

PRO free_diag_memory

  COMMON global_vars

  IF PTR_VALID(gui.tab.arr) THEN $
   FOR itab=0, N_ELEMENTS((*gui.tab.arr))-1 DO BEGIN
    diag = (*gui.tab.arr)[itab].diag0
    WHILE PTR_VALID(diag) DO BEGIN
      diag_del = diag
      diag = (*diag).next
      PTR_FREE, diag_del
    ENDWHILE
   ENDFOR

  IF PTR_VALID(cdiag0) THEN BEGIN
   cdiag = cdiag0
   WHILE PTR_VALID(cdiag) DO BEGIN
     cdiag_del = cdiag
     cdiag = (*cdiag).next
     PTR_FREE, cdiag_del
   ENDWHILE
  ENDIF
END
