FUNCTION rm0es, flt_in, prec=prec, int_digits=int_digits
; This function takes a float/integer/string and returns a string
; with blanks and leading/trailing zeroes removed. If arrays of
; floats/integers/strings are provided, string arrays of same
; dimension are returned. If prec is set to an integer value, prec
; digits after the decimal point are returned. Note that setting
; prec=0 does not perform a ROUND conversion, use that function
; instead to get the nearest integer. For floats with exponent,
; the number in front of the exponent is rm0esed.

  IF NOT KEYWORD_SET(int_digits) THEN int_digits = 0

  str_out_arr = STRARR(N_ELEMENTS(flt_in))

  ; catch byte-type input
  IF SIZE(flt_in,/TYPE) EQ 1 THEN flt_in = LONG(flt_in)

  FOR k = 1, N_ELEMENTS(flt_in) DO BEGIN
    str_temp = STRTRIM(STRING(flt_in[k-1]),2)
    str_temp_length = STRLEN(str_temp)

    intdetect = 1B
    FOR i = 0, str_temp_length - 1 DO $
      IF STRMID(str_temp,i,1) EQ '.' THEN intdetect = 0B

    IF intdetect EQ 1B THEN BEGIN
      str_out = str_temp

      IF int_digits THEN WHILE $
        STRLEN(str_out) LT int_digits DO str_out = '0' + str_out
    ENDIF ELSE BEGIN
      exponent = ''
      exponent_split = STRSPLIT(str_temp,'eE',/EXTRACT)
      IF N_ELEMENTS(exponent_split) GT 1 THEN BEGIN
        str_temp = exponent_split[0]
        exponent = 'e' + exponent_split[1]
        str_temp_length = STRLEN(str_temp)
      ENDIF

      IF KEYWORD_SET(prec) THEN BEGIN
        flt_mod = FLOAT(str_temp)
        flt_mod = FLOAT(ROUND(flt_mod*10.0^prec)) * 10.0^(-prec)
        str_temp = STRTRIM(STRING(flt_mod),2)
        str_temp_length = STRLEN(str_temp)
      ENDIF

      act_char = '0'
      i = 1
      WHILE (act_char EQ '0') AND (i LE str_temp_length) DO BEGIN
        j = str_temp_length - i
        act_char = STRMID(str_temp,j,1)
        i = i + 1
      ENDWHILE
      str_out = STRMID(str_temp,0,j+1)
      str_out_length = STRLEN(str_out)
      IF STRMID(str_out,str_out_length-1,1) EQ '.' THEN $
        str_out = STRMID(str_out,0,str_out_length-1)

      IF int_digits THEN WHILE $
        STRLEN((STRSPLIT(str_out,'.',/EXTRACT))[0]) LT $
        int_digits DO str_out = '0' + str_out

      str_out = str_out + exponent
    ENDELSE

    str_out_arr[k-1] = str_out
  ENDFOR

  IF N_ELEMENTS(flt_in) EQ 1 THEN RETURN, str_out_arr[0]
  RETURN, str_out_arr

END

;##########################################################################

FUNCTION append_th, long_in, fancy=fancy
; appends st, nd, rd, or th to a number
; fancy : set for fancy style with exponent

  IF N_ELEMENTS(long_in) LT 1 THEN RETURN, 'error'
  long_in = LONG(long_in)
  sign = 1 * (long_in GE 0) - (long_in LT 0)
  long_in = ABS(long_in)

  append = REPLICATE('th',N_ELEMENTS(long_in))
  long_in_1 = WHERE((long_in MOD 10) EQ 1)
  IF long_in_1[0] NE -1 THEN append[long_in_1] = 'st'
  long_in_2 = WHERE((long_in MOD 10) EQ 2)
  IF long_in_2[0] NE -1 THEN append[long_in_2] = 'nd'
  long_in_3 = WHERE((long_in MOD 10) EQ 3)
  IF long_in_3[0] NE -1 THEN append[long_in_3] = 'rd'
  long_in_11 = WHERE(long_in EQ 11)
  IF long_in_11[0] NE -1 THEN append[long_in_11] = 'th'
  long_in_12 = WHERE(long_in EQ 12)
  IF long_in_12[0] NE -1 THEN append[long_in_12] = 'th'
  long_in_13 = WHERE(long_in EQ 13)
  IF long_in_13[0] NE -1 THEN append[long_in_13] = 'th'

  IF KEYWORD_SET(fancy) THEN append = '!U' + append + '!N'

  RETURN, rm0es(sign*long_in) + append

END

;##########################################################################

FUNCTION rm_idl_fmt, in_str
; removes IDL format substring

  out_str = in_str

  FOR j = 0, N_ELEMENTS(out_str) - 1 DO BEGIN
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!D',/REGEX,/EXTRACT),'_')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!I',/REGEX,/EXTRACT),'_')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!9#',/REGEX,/EXTRACT),'par')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!9x',/REGEX,/EXTRACT),'perp')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!7q',/REGEX,/EXTRACT),'rho')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!7l',/REGEX,/EXTRACT),'mu')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!7U',/REGEX,/EXTRACT),'phi')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!N',/REGEX,/EXTRACT),'')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!12',/REGEX,/EXTRACT),'')
    out_str[j] = STRJOIN(STRSPLIT(out_str[j],'!6',/REGEX,/EXTRACT),'')
  ENDFOR

  RETURN, out_str

END

;##########################################################################

FUNCTION get_tag_index, in_struct, in_tag
; get index of tag in structure
; to be used like: struct.(tag_index)

  my_tag = STRUPCASE(in_tag)
  tagnames = TAG_NAMES(in_struct)

  RETURN, (WHERE(tagnames EQ my_tag))

END

;##########################################################################

FUNCTION time_avg, id, value, time, avg=avg, fwd=fwd, tarr=tarr, intrp=intrp,$
                   rms=rms
; Time average calculation
; This function sums up value (which can be a number or a
; (multi-dimensional) array) weighed with the time interval. The id
; pointer that identifies a specific summation is returned. The variable
; time corresponds to the data step time. If /avg is set, no further
; summation is performed (value and time are ignored), and the result is
; returned.
; If /fwd is specified (usually via gui.out.res_steps), the store_step
; function is used instead to supply time-resolved data. Additionally,
; a time array is constructed and returned via tarr.
; The time-resolved data will be interpolated to an equidistant time
; array with intrp*number of original time steps if intrp is specified
; and >0
; If /rms is specified, the root-mean-square will be computed instead
; of a simple time average; usage of /avg for final output is still required

  FORWARD_FUNCTION store_step

  IF KEYWORD_SET(fwd) THEN BEGIN
    IF NOT KEYWORD_SET(avg) THEN BEGIN
      IF PTR_VALID(id) THEN BEGIN
        time_id = (*id)[0]
        data_id = (*id)[1]
      ENDIF
      time_id = store_step(time_id,time)
      data_id = store_step(data_id,value)

      res = PTR_NEW([time_id,data_id])
    ENDIF ELSE BEGIN
      time_id = (*id)[0]
      data_id = (*id)[1]

      IF KEYWORD_SET(intrp) THEN BEGIN
         oldtime = store_step(time_id,/get,/reg_array)
         tdim = N_ELEMENTS(oldtime)
         dt = (oldtime[tdim-1]-oldtime[0])/(intrp*tdim-1)
         tarr = oldtime[0] + FINDGEN(intrp*tdim)*dt
         nointrp = 0
         IF (intrp EQ 1) THEN IF (TOTAL(oldtime-tarr) EQ 0) THEN nointrp = 1
         IF (nointrp) THEN BEGIN
            res = store_step(data_id,/get,/reg_array)
         ENDIF ELSE BEGIN
            print, 'interpolating to equidistant time trace with factor '+$
                   rm0es(intrp)
            oldres = store_step(data_id,/get,/reg_array)
            oldsize = SIZE(oldres)
            rdim = oldsize[oldsize[0]+2]/tdim
            oldres = REFORM(oldres,[rdim,tdim])
            IF (2 EQ 1) THEN BEGIN ;saves memory but is awfully slow 
                                   ;(left here for mem. limited cases)
               tmparr = INTERPOL(oldres[0,*],oldtime,tarr)
               oldres = oldres[1:*,*]
               res = TRANSPOSE(tmparr)

               FOR ir = 1, rdim - 1 DO BEGIN
                  tmparr = INTERPOL(oldres[0,*],oldtime,tarr)
                  IF (ir LT (rdim-1)) THEN oldres = oldres[1:*,*]
                  res = [[res],TRANSPOSE(tmparr)]
               ENDFOR
            ENDIF ELSE BEGIN
               res = FLTARR(rdim, N_ELEMENTS(tarr),/NOZERO)
               FOR ir = 0, rdim - 1 DO res[ir,*] = $
                  INTERPOL(oldres[ir,*],oldtime,tarr)
            ENDELSE
            res = REFORM(res,[oldsize[1:oldsize[0]-1],N_ELEMENTS(tarr)])
         ENDELSE
      ENDIF ELSE BEGIN
         tarr = store_step(time_id,/get,/reg_array)
         res = store_step(data_id,/get,/reg_array)
      ENDELSE
    ENDELSE

    PTR_FREE, id

    RETURN, res
  ENDIF

  IF NOT PTR_VALID(id) THEN BEGIN
    IF KEYWORD_SET(avg) THEN BEGIN
      PRINT, 'time_avg error: average requested without id pointer'
      RETURN, 0.0
    ENDIF
    IF (N_ELEMENTS(value) LT 1) OR (N_ELEMENTS(time) LT 1) THEN BEGIN
      PRINT, 'time_avg error: need to specify value and time'
      RETURN, 0.0
    ENDIF

    IF KEYWORD_SET(rms) THEN value = TEMPORARY(value)^2

    RETURN, PTR_NEW({value_sum:value,time_sav:-time,time_first:time})
  ENDIF

  IF KEYWORD_SET(avg) THEN BEGIN
    IF (*id).time_sav LE 0.0 THEN RETURN, (*id).value_sum ; only one step
    timesize = SIZE((*id).time_sav)
    prec_single = 1
    IF N_ELEMENTS(timesize) EQ 3 THEN $
      IF timesize[1] EQ 5 THEN prec_single = 0

    time_range_inv = (prec_single ? 1.0 : 1.0D) / $
      ((*id).time_sav - (*id).time_first)

    value_avg = (*id).value_sum * time_range_inv
    PTR_FREE, id

    IF KEYWORD_SET(rms) THEN value_avg = SQRT(value_avg)

    RETURN, value_avg
  ENDIF

  IF (*id).time_sav LE 0.0 THEN BEGIN ; first step: use second dt
    (*id).time_sav = - (*id).time_sav
    (*id).value_sum = TEMPORARY((*id).value_sum) * (time - (*id).time_sav)
    (*id).time_first = (*id).time_first - (time - (*id).time_sav)
  ENDIF

  IF KEYWORD_SET(rms) THEN value = ABS(TEMPORARY(value))^2 * (time - (*id).time_sav) $
  ELSE value = TEMPORARY(value) * (time - (*id).time_sav)

  (*id).value_sum = TEMPORARY((*id).value_sum) + TEMPORARY(value)

  (*id).time_sav = time

  RETURN, id

END

;##########################################################################

FUNCTION store_step, id, step, get=get, reg_array=reg_array
; Step vault
; Stores consecutive single or multi dimensional arrays and returns
; a pointer id that is used to reference the data. Upon completion
; of the step storage, one can retrieve the data via the /get keyword.
; If /reg_array is specified, the data is returned as an array
; of one more dimension than the input steps; otherwise, the return
; value is an array with as many elements as steps were stored, with
; each element being a pointer to the corresponding step data.
; Note that trailing degenerated dimensions are preserved if present
; in the incoming step.

  IF NOT PTR_VALID(id) THEN BEGIN
    IF KEYWORD_SET(get) THEN BEGIN
      PRINT, 'store_step error: data requested without id pointer'
      RETURN, 0.0
    ENDIF
    IF N_ELEMENTS(step) LT 1 THEN BEGIN
      PRINT, 'store_step error: need to specify step data'
      RETURN, 0.0
    ENDIF

    step_size = SIZE(step)
    n_dim = step_size[0]
    IF n_dim EQ 0 THEN dimensions = 1 ELSE dimensions = step_size[1:n_dim]

    RETURN, PTR_NEW([PTR_NEW(dimensions),PTR_NEW(step)])
  ENDIF

  IF KEYWORD_SET(get) THEN BEGIN
    IF KEYWORD_SET(reg_array) THEN BEGIN
      steps = (*id)[1:*]
      n_steps = N_ELEMENTS(steps)
      step_size = SIZE(*(steps[0]))
      n_dim = step_size[0]

      dimensions = (*id)[0]
      PTR_FREE, id

      elements = step_size[n_dim+2]
      type = step_size[n_dim+1]
      data = MAKE_ARRAY(DIMENSION=[elements,n_steps],TYPE=type,/NOZERO)

      FOR n = 0, n_steps - 1 DO BEGIN
        data[0,n] = REFORM(TEMPORARY(*steps[n]),elements,/OVERWRITE)

        PTR_FREE, steps[n]
      ENDFOR

      data = n_dim EQ 0 ? REFORM(data,/OVERWRITE) : $
        REFORM(data,[*dimensions,n_steps],/OVERWRITE)

      PTR_FREE, dimensions

      RETURN, data
    ENDIF

    PTR_FREE, (*id)[0]
    data = (*id)[1:*]
    PTR_FREE, id
    RETURN, data
  ENDIF

  data = *id
  PTR_FREE, id

  RETURN, PTR_NEW([data,PTR_NEW(step)])

END

;##########################################################################

PRO plot_info_str, diag, time=time
; plots strings to the bottom of the plot pages

  FORWARD_FUNCTION get_var_string

  COMMON global_vars

  IF cfg.info_str EQ 0 THEN RETURN

  ; obtain diag variable info string
  IF (cfg.info_str GT 1) AND PTR_VALID(diag) THEN BEGIN
    vinfo_str = '!3'

    diag_info = CALL_FUNCTION((*diag).name+'_info')
    n_vars = diag_info.ext_vars[0,0] EQ '' ? 0 : $
      N_ELEMENTS(diag_info.ext_vars[0,*])
    e = (*diag).external_vars

    max_neval = 8 ; maximal number of array elements to be shown

    FOR v = 0, n_vars - 1 DO BEGIN
      ; syntax *e.(v) means addressing a structure like an array
      val = N_ELEMENTS(*e.(v)) LT 1 ? 'n/s' : *e.(v)
      IF diag_info.ext_vars[1,v] EQ '1' THEN BEGIN
        IF N_ELEMENTS(val) NE 1 THEN val = 0
      ENDIF ELSE BEGIN
        ne_val = N_ELEMENTS(val)
        IF ne_val GT 1 THEN BEGIN
          val_temp = '['
          n_max = ne_val < max_neval
          FOR n = 0, n_max - 1 DO val_temp += $
            ((n EQ n_max - 1) AND (max_neval LT ne_val) ? $
            '...' : rm0es(val[n])) + (n EQ n_max - 1 ? ']' : ',')
          val = val_temp
        ENDIF
      ENDELSE

      vinfo_str += (v EQ 0 ? '' : ',') + $
        diag_info.ext_vars[0,v] + '=' + rm0es(val)
    ENDFOR

    cs = STRLEN(vinfo_str) GT 60 ? 0.7 : 1.0
  ENDIF

  ; save current colors, load new table
  COMMON COLORS, R_ORIG, G_ORIG, B_ORIG, R_CURR, G_CURR, B_CURR
  r_store = R_CURR
  g_store = G_CURR
  b_store = B_CURR
  LOADCT, 41, FILE='internal/colortable.tbl'

  IF (cfg.info_str GT 1) AND PTR_VALID(diag) THEN XYOUTS, 0.995, 0.005, $
    vinfo_str, /NORMAL, COLOR=9, CHARSIZE=cs, ALIGNMENT=1, /NOCLIP

  ; obtain time info string
  IF N_ELEMENTS(time) EQ 1 THEN BEGIN
    IF time LT 0.0 THEN BEGIN
      gamma_omega = rm0es(get_eigenvalue(-time),prec=3)
      tinfo_str = '!3EV' + rm0es(-time) + ': !4c,x!3=' + $
        gamma_omega[0] + ',' + gamma_omega[1]
    ENDIF ELSE tinfo_str = $
      '!3t=' + rm0es(time) + get_var_string(1,/time,/ounit)
  ENDIF ELSE tinfo_str = '!3t!9e!3[' + rm0es(gui.out.start_t) + ',' + $
    rm0es(gui.out.end_t) + ']' + get_var_string(1,/time,/ounit) + $
    ',f!Dsparse!N=' + rm0es(gui.out.sparsefac)

  IF N_ELEMENTS(cs) NE 1 THEN cs = 1.0
  XYOUTS, 0.005, 0.005, tinfo_str, /NORMAL, COLOR=1, CHARSIZE=cs, /NOCLIP

  ; restore old colors
  TVLCT, r_store, g_store, b_store
  R_CURR = r_store
  G_CURR = g_store
  B_CURR = b_store

END

;##########################################################################

PRO print_progress_cl, currval, maxval, first=first, $
  erase_at_100=erase_at_100, header_str=header_str
; Prints the progress percentage of currval relative to maxval; note: the
; routine recognizes the last call by (currval GE maxval), thus it is
; advisible to ensure this point is reached.
; In order to ensure correct display of the percentages, it is recommended
; that only a single external routine call the present routine at a time
; (i.e., until the progress of that external routine reaches 100%).
; first:        needs to be set at first step
; erase_at_100: remove printed output upon reaching 100%
; header_str:   put a particular string in front of the 'Progress:' string

  IF (N_ELEMENTS(currval) NE 1) OR (N_ELEMENTS(maxval) NE 1) THEN BEGIN
    PRINT, 'print_progress_cl: invalid input parameters'
    RETURN
  ENDIF

  IF N_ELEMENTS(header_str) NE 1 THEN header_str = ''

  currnorm = STRTRIM(STRING(100.0*currval/maxval),2)
  IF currval LT 0.1 * maxval THEN currnorm = ' ' + currnorm

  n_digits = 2 ; number of digits for percentage (does not apply to '100%')
  currnorm = STRMID(currnorm,0,n_digits) + '%'

  IF KEYWORD_SET(first) THEN BEGIN
    PRINT, FORMAT='(A,$)', header_str + 'Progress: ' + currnorm
    RETURN
  ENDIF

  ; string with (n_digits + 1) backspace characters to erase
  ; previous percentage
  bspaces = STRJOIN(REPLICATE(STRING(8B),n_digits+1))
  IF currval LT maxval THEN BEGIN
    PRINT, FORMAT='(A,$)', bspaces + currnorm
    RETURN
  ENDIF

  ; completion (100%)
  bspaces = STRJOIN(REPLICATE(STRING(8B),n_digits+1+10+STRLEN(header_str)))
  IF NOT KEYWORD_SET(erase_at_100) THEN PRINT, FORMAT='(A)', $
    bspaces + '100%' ELSE PRINT, FORMAT='(A,$)', bspaces

END

;##########################################################################

FUNCTION get_eigenvalue, n, run

  COMMON global_vars

  IF N_ELEMENTS(n) NE 1 THEN n = 1
  IF N_ELEMENTS(run) NE 1 THEN run = series.run_labels[0]

  OPENR, ev_lun, gui.out.data_path + 'eigenvalues_' + rm0es(run), $
    /GET_LUN, ERR=err
  IF err NE 0 THEN BEGIN
    PRINT, 'error opening eigenvalues_' + rm0es(run)
    RETURN, [0.0,0.0]
  ENDIF

  line = ''
  found_ev_line = 0
  WHILE NOT found_ev_line AND NOT EOF(ev_lun) DO BEGIN
    READF, ev_lun, line
    IF (STRLEN(line) GE 11) THEN $
       found_ev_line = (STRPOS(line,'eigenvalues') GE 0)
  ENDWHILE
  IF NOT found_ev_line THEN BEGIN
    PRINT, 'error: no eigenvalues in eigenvalues_' + rm0es(run)
    RETURN, [0.0,0.0]
  ENDIF

  ev_count = 1
  ev_value = FLTARR(2)
  WHILE NOT EOF(ev_lun) AND (ev_count LE n) DO BEGIN
    READF, ev_lun, line

    IF ev_count EQ n THEN BEGIN
      linesep = STRSPLIT(STRTRIM(line,2),' ',/EXTRACT)
      IF N_ELEMENTS(linesep) LT 2 THEN BEGIN
        PRINT, 'error: wrong eigenvalues format in eigenvalues_' + rm0es(run)
        RETURN, [0.0,0.0]
      ENDIF
      ev_value[*] = FLOAT(linesep[0:1])
    ENDIF

    ev_count += 1
  ENDWHILE

  IF ev_count LE n THEN PRINT, 'error: not enough eigenvalues found in ' + $
    'eigenvalues_' + rm0es(run)

  RETURN, ev_value

END

;##########################################################################

PRO undefine, varname

  dummy = SIZE(TEMPORARY(varname))

END

;##########################################################################

PRO add2arr, arr, elem

  arr = N_ELEMENTS(arr) LT 1 ? elem : [arr,elem]

END

;##########################################################################

FUNCTION pf_arr, dims, zero=zero, dbl=dbl, cmplx=cmplx
; proper fast array allocation
; returns an FLTARR with /NOZERO and preserved trailing hollow dimensions,
;   e.g.: [24,16,1] rather than [24,16] as would FLTARR
; dims: array of dimensions
; zero: set in order not to use /NOZERO, i.e.: zero all elements
; dbl:  use double precision, i.e.: DBLARR rather than FLTARR

  ndims = N_ELEMENTS(dims)
  IF ndims LT 1 OR ndims GT 8 THEN BEGIN
    PRINT, 'pf_arr: too few/many dimensions'
    RETURN, -1
  ENDIF

  IF (WHERE(dims LE 0))[0] NE -1 THEN BEGIN
    PRINT, 'pf_arr: zero or negative dimension sizes'
    RETURN, -2
  ENDIF

  nozero = NOT KEYWORD_SET(zero) ? 1 : 0

  IF NOT KEYWORD_SET(cmplx) THEN BEGIN
    RETURN, NOT KEYWORD_SET(dbl) ? $
      REFORM(FLTARR(dims,NOZERO=nozero),dims,/OVERWRITE) : $
      REFORM(DBLARR(dims,NOZERO=nozero),dims,/OVERWRITE)
  ENDIF ELSE BEGIN
    RETURN, NOT KEYWORD_SET(dbl) ? $
      REFORM(COMPLEXARR(dims,NOZERO=nozero),dims,/OVERWRITE) : $
      REFORM(DCOMPLEXARR(dims,NOZERO=nozero),dims,/OVERWRITE)
  ENDELSE

END

;##########################################################################

FUNCTION str2fltarr, str_var
; converts input string to float array
; comma, colon, basic math operators, and some par values are recognized
; note: math operators will be treated in serial order (e.g., '1+4*5' equals 25)
; brackets will be removed and ignored

  COMMON global_vars

  IF STRPOS(str_var,'[') EQ 0 THEN $
    str_var = STRMID(str_var,1,STRLEN(str_var)-2)
  IF STRPOS(str_var,'(') EQ 0 THEN $
    str_var = STRMID(str_var,1,STRLEN(str_var)-2)
  IF STRPOS(str_var,',') GE 0 THEN BEGIN
    comma_split = STRSPLIT(str_var,',',/EXTRACT,COUNT=n_comma)
    FOR j = 0, n_comma - 1 DO $
      add2arr, out_arr, str2fltarr(comma_split[j])
    RETURN, out_arr
  ENDIF
  IF STRPOS(str_var,':') GE 0 THEN BEGIN
    colon_split = STRSPLIT(str_var,':',/EXTRACT,COUNT=n_colon)
    IF n_colon GT 2 THEN BEGIN
      printerror, 'invalid string'
      RETURN, 0
    ENDIF ELSE IF n_colon EQ 2 THEN BEGIN
      add2arr, out_arr, (FLOAT(str2fltarr(colon_split[0]))+$
	INDGEN(FIX(str2fltarr(colon_split[1]))-$
	FIX(str2fltarr(colon_split[0]))+1))
    ENDIF ELSE IF n_colon EQ 1 THEN $
      add2arr, out_arr, FLOAT(str2fltarr(colon_split[0]))
    RETURN, out_arr
  ENDIF
  math_op = ['+','-','/','*']
  FOR j = 0, 3 DO BEGIN
    IF STRPOS(str_var,math_op[j]) GE 0 THEN BEGIN
      math_split = STRSPLIT(str_var,math_op[j],/EXTRACT,COUNT=n_op)
      IF n_op EQ 2 THEN BEGIN
        num1 = str2fltarr(math_split[0])
	num2 = str2fltarr(math_split[1])
	CASE math_op[j] OF
          '*' : add2arr, out_arr, num1*num2
	  '/' : add2arr, out_arr, num1/num2
	  '+' : add2arr, out_arr, num1+num2
	  '-' : add2arr, out_arr, num1-num2
	ENDCASE
	RETURN, out_arr
      ENDIF
      IF STRSPLIT(str_var,math_op[j]) EQ 1 THEN BEGIN
        num1 = str2fltarr(math_split[0])
        CASE math_op[j] OF
          '*' : add2arr, out_arr, num1
	  '/' : add2arr, out_arr, 1/num1
	  '+' : add2arr, out_arr, num1
	  '-' : add2arr, out_arr, -num1
	ENDCASE
	RETURN, out_arr
      ENDIF

      printerror, 'invalid string'
      RETURN, 0
    ENDIF
  ENDFOR

  IF (str_var EQ 'par.nkx0') OR (str_var EQ 'nkx0') THEN RETURN, par.nkx0
  IF (str_var EQ 'par.nx0') OR (str_var EQ 'nx0') THEN RETURN, par.nx0
  IF (str_var EQ 'par.nky0') OR (str_var EQ 'nky0') THEN RETURN, par.nky0
  IF (str_var EQ 'par.ny0') OR (str_var EQ 'ny0') THEN RETURN, par.ny0
  IF (str_var EQ 'par.nz0') OR (str_var EQ 'nz0') THEN RETURN, par.nz0
  IF (str_var EQ 'par.nw0') OR (str_var EQ 'nw0') THEN RETURN, par.nw0
  IF (str_var EQ 'par.nv0') OR (str_var EQ 'nv0') THEN RETURN, par.nv0
  IF (str_var EQ 'par.n_fields') OR (str_var EQ 'nf') THEN RETURN, par.n_fields
  IF (str_var EQ 'par.n_moms') OR (str_var EQ 'nm') THEN RETURN, par.n_moms

  RETURN, FLOAT(str_var)

END

;##########################################################################

FUNCTION convert_table_entry, str_var, dummy=dummy, ptr=ptr, int=int

  IF str_var NE '' THEN BEGIN
    IF KEYWORD_SET(ptr) THEN RETURN, PTR_NEW(FLOAT(str2fltarr(str_var)))
    IF KEYWORD_SET(int) THEN RETURN, FIX(str2fltarr(str_var))
    RETURN, str2fltarr(str_var)
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(ptr) THEN RETURN, PTR_NEW(dummy)
    RETURN, dummy
  ENDELSE

END

;##########################################################################

FUNCTION get_var_string, var_index, fancy=fancy, units=units, ounit=ounit,$
  time=time, normstr=normstr, rhostr=rhostr, sp = sp
; takes an index or an array of indices and returns the corresponding
; label(s); if /units is specified, fancy versions are returned
; ounit: returns only the unit in fancy style

  COMMON global_vars

  isp = 0

  Tref_is_Te = 0
  Tref_is_Ti = 0
  mref_is_me = 0
  mref_is_mi = 0
  Lref_str = series.Lref_str

  WHILE (isp LT par.n_spec) DO BEGIN
    Tref_is_Te = Tref_is_Te OR (spec[isp].charge EQ -1 $
      AND ((spec[isp].temp*series.Tref GT 0.9999) AND $
      (spec[isp].temp*series.Tref LT 1.0001)))
    Tref_is_Ti = Tref_is_Ti OR (spec[isp].charge EQ 1 $
      AND ((spec[isp].temp*series.Tref GT 0.9999) AND $
      (spec[isp].temp*series.Tref LT 1.0001)) AND $
      (spec[isp].dens EQ 1.0))
    mref_is_me = mref_is_me OR (spec[isp].charge EQ -1 $
      AND ((spec[isp].mass*series.mref GT 0.9999) AND $
      (spec[isp].mass*series.mref LT 1.0001)))
    mref_is_mi = mref_is_mi OR (spec[isp].charge EQ 1 $
      AND ((spec[isp].mass*series.mref GT 0.9999) AND $
      (spec[isp].mass*series.mref LT 1.0001)) AND $
      (spec[isp].dens EQ 1.0))
    isp += 1
  ENDWHILE

  IF series.Tref EQ 1.0 AND series.mref EQ 1 THEN ref_str = 'ref' $
  ELSE ref_str = ''

  IF Tref_is_Ti THEN T_ref_str = '!6T!Di0!N'
  IF Tref_is_Te THEN T_ref_str = '!6T!De0!N'

  IF Tref_is_Te AND mref_is_me THEN BEGIN
    v_ref_str = '!6v!Dte!N'
    sp_ref_str = 'e'
    rho_ref_str = '!7q!6!De!N'
  ENDIF ELSE IF Tref_is_Te AND mref_is_mi THEN BEGIN
    v_ref_str = '!6c!Ds!N'
    sp_ref_str = 'i'
    rho_ref_str = '!7q!6!Ds!N'
  ENDIF ELSE IF Tref_is_Ti AND mref_is_mi THEN BEGIN
    v_ref_str = '!6v!Dti!N'
    sp_ref_str = 'i'
    rho_ref_str = '!7q!6!Di!N'
  ENDIF ELSE BEGIN
    v_ref_str = '!6v!Dt,'+ref_str+'!N'
    sp_ref_str = ''
    rho_ref_str = '!7q!6!D'+ref_str+'!N'
  ENDELSE

  IF KEYWORD_SET(rhostr) THEN  RETURN, rho_ref_str
  IF KEYWORD_SET(normstr) THEN RETURN, Lref_str

  IF KEYWORD_SET(time) THEN BEGIN
    IF var_index EQ -1 THEN BEGIN
      IF KEYWORD_SET(ounit) THEN RETURN, v_ref_str + '/' + Lref_str
      IF KEYWORD_SET(units) THEN RETURN, '!7c!6 / ('+v_ref_str+'/'+Lref_str+')'
      IF KEYWORD_SET(fancy) THEN RETURN, '!7c!6' $
        ELSE RETURN, 'gamma'
    ENDIF ELSE BEGIN
      IF KEYWORD_SET(ounit) THEN RETURN, Lref_str+ '/'+ v_ref_str
      IF KEYWORD_SET(units) THEN RETURN, '!6t / ('+Lref_str+ '/'+ v_ref_str+')'
      RETURN, 't'
    ENDELSE
  ENDIF

  IF N_ELEMENTS(sp) NE 1 THEN sp_str = '' $
  ELSE IF sp LT par.n_spec THEN BEGIN
    IF spec[sp].name EQ 'electrons' THEN sp_str = 'e' $
    ELSE IF spec[sp].name EQ 'ions' THEN sp_str = 'i' $
    ELSE sp_str = spec[sp].name
  ENDIF

  nf = par.n_fields
  v = STRARR(200) ;3+par.n_moms)
  vf = STRARR(200) ;3+par.n_moms)
  vu = STRARR(200) ;3+par.n_moms)

  v[0] = 'phi'
  v[1] = 'Apar'
  v[2] = 'Bpar'
  v[nf] = 'n'
  v[nf+1] = 'Tpar'
  v[nf+2] = 'Tperp'
  v[nf+3] = 'qpar'
  v[nf+4] = 'qperp'
  v[nf+5] = 'upar'
  IF par.n_moms GT 6 THEN BEGIN
    v[nf+6] = 'N00'
    v[nf+7] = 'N20'
    v[nf+8] = 'N02'
  ENDIF

  vf[0] = '!7U!6'
  vf[1] = '!6A!D!9#!N!6'
  vf[2] = '!6B!D!9#!N!6'
  vf[nf] = '!6n'
  vf[nf+1] = '!6T!D!9#!6!N'
  vf[nf+2] = '!6T!D!9x!6!N'
  vf[nf+3] = '!6q!D!9#!6!N'
  vf[nf+4] = '!6q!D!9x!6!N'
  vf[nf+5] = '!6u!D!9#!6!N'
  IF par.n_moms GT 6 THEN BEGIN
    vf[nf+6] = '!6N!D00!N'
    vf[nf+7] = '!6N!D20!N'
    vf[nf+8] = '!6N!D02!N'
  ENDIF

  IF sp_str NE '' THEN BEGIN
    v[nf:nf+5] += '_' + sp_str
    vf[nf:nf+5] += '!D' + sp_str + '!N'
  ENDIF

  n_moms_base = par.n_moms EQ 9 ? 9 : 6

  FOR j = 1, par.n_moms / n_moms_base - 2 DO BEGIN
    v[nf+j*n_moms_base:nf+j*n_moms_base+n_moms_base-1] = $
      v[nf:nf+n_moms_base-1] + ',trap' + rm0es(par.n_moms/n_moms_base-2-j)
    vf[nf+j*n_moms_base:nf+j*n_moms_base+n_moms_base-1] = $
      vf[nf:nf+n_moms_base-1] + '!D,t' + $
      rm0es(par.n_moms/n_moms_base-2-j) + '!N'
  ENDFOR



  IF par.n_moms GT n_moms_base THEN BEGIN
    v[nf+(par.n_moms/n_moms_base-1)*n_moms_base:$
      nf+(par.n_moms/n_moms_base-1)*n_moms_base+n_moms_base-1] = $
      v[nf:nf+n_moms_base-1] + ',flr'
    vf[nf+(par.n_moms/n_moms_base-1)*n_moms_base:$
      nf+(par.n_moms/n_moms_base-1)*n_moms_base+n_moms_base-1] = $
      vf[nf:nf+n_moms_base-1] + '!D,flr!N'
    v[nf:nf+n_moms_base-1] = v[nf:nf+n_moms_base-1] + ',pass'
    vf[nf:nf+n_moms_base-1] = vf[nf:nf+n_moms_base-1] + '!D,p!N'
  ENDIF

  FOR j = 0, par.n_moms / n_moms_base - 1 DO BEGIN
    vf[nf+j*n_moms_base+3] += '+1.5p!Dj0!N' + vf[nf+n_moms_base-1]
    vf[nf+j*n_moms_base+4] += '+p!Dj0!N' + vf[nf+n_moms_base-1]
  ENDFOR

  ;derived variables
  v[100] = 'T_tot'
  v[101] = 'p_j'
  v[102] = 'd^2/dx^2 '+v[0]
  v[103] = '-d/dx '+v[100]
  v[104] = '-d/dx '+v[nf]
  v[199] = 'sum q_j n_j'

  vf[100] = '!6T!Dtot!N'
  vf[101] = 'p!Dj1!N'
  vf[102] = '!6!S!E !N!Ud!E2!N!R!S-!R!Ddx!E2!N!X!N'+vf[0]
  vf[103] = '!6!S!E !N!Ud!R!S-!R!Ddx!X!N'+vf[100]
  vf[104] = '!6!S!E !N!Ud!R!S-!R!Ddx!X!N'+vf[nf]
  vf[199] = '!7R!6!Dj!Nq!Dj!Nn!Dj!N'

  sp_ref_str = '' ; n and T normalized to (n0j * nref) not to nref ???
  vu[0] = T_ref_str + rho_ref_str + '/(e' + Lref_str + ')'
  vu[1] = 'B!D' + ref_str + '!N' + rho_ref_str + '!U2!N/' + Lref_str
  vu[2] = 'B!D' + ref_str + '!N' + rho_ref_str + '/' + Lref_str
  vu[nf] = 'n!Dj0!N' + rho_ref_str + '/' + Lref_str
  vu[nf+1] = 'T!Dj0!N' + rho_ref_str + '/' + Lref_str
  vu[nf+2] = 'T!Dj0!N' + rho_ref_str + '/' + Lref_str
  vu[nf+3] = 'p!Dj0!N' + v_ref_str + rho_ref_str + '/' + Lref_str
  vu[nf+4] = 'p!Dj0!N' + v_ref_str + rho_ref_str + '/' + Lref_str
  vu[nf+5] = v_ref_str + rho_ref_str + '/' + Lref_str
  IF par.n_moms GT 6 THEN BEGIN
    vu[nf+6] = '(p!Dj0!N/B!D' + ref_str + '!N) (' + $
      rho_ref_str + '/' + Lref_str + ')'
    vu[nf+7] = '(p!Dj0!N' + v_ref_str + '!U2!N' + '/B!D' + ref_str + $
      '!N) (' + rho_ref_str + '/' + Lref_str + ')'
    vu[nf+8] = '(p!Dj0!N' + v_ref_str + '!U2!N' + '/B!D' + ref_str + $
      '!N) (' + rho_ref_str + '/' + Lref_str + ')'
  ENDIF
  FOR j = 1, par.n_moms / n_moms_base - 1 DO $
    vu[nf+j*n_moms_base:nf+j*n_moms_base+n_moms_base-1] = $
    vu[nf:nf+n_moms_base-1]

  ;derived variables
  vu[100] = vu[nf+1]
  vu[101] = 'p!Dref!N' + rho_ref_str + '/' + Lref_str
  vu[102] = 'T!Dref0!N/(e' + rho_ref_str + Lref_str + ')'
  vu[103] = 'T!Dj0!N/' + Lref_str
  vu[104] = 'n!Dj0!N/' + Lref_str
  vu[199] = 'n!Dref!N e'

  IF KEYWORD_SET(ounit) THEN RETURN, vu[var_index]
  IF KEYWORD_SET(units) THEN RETURN, vf[var_index] + ' / (' + $
    vu[var_index] + ')'
  IF KEYWORD_SET(fancy) THEN RETURN, vf[var_index] ELSE $
    RETURN, v[var_index]

END

;######################################################################

FUNCTION betterticks, axis, index, value

  RETURN, STRING(FORMAT='(A,I0,A)','!610!U',ROUND(ALOG10(value)),'!N')

END

;######################################################################

FUNCTION get_file_path, name, run, set_fmt=set_fmt
; input parameters:
; name : string (array) which contains the file type, e.g.: 'field'
; run  : integer (array) which contains index of run_label
; returns file path/exist structure

  COMMON global_vars

  IF N_ELEMENTS(run) LT 1 OR KEYWORD_SET(set_fmt) THEN BEGIN
    runstr = (STRSPLIT(series.run_labels,',',/EXTRACT,COUNT=n_runs))
    run = INDGEN(n_runs)
  ENDIF ELSE BEGIN
    runstr = (STRSPLIT(series.run_labels,',',/EXTRACT))[run]
    n_runs = N_ELEMENTS(runstr)
  ENDELSE

  data_path = gui.out.data_path

  IF KEYWORD_SET(set_fmt) THEN BEGIN
    PTR_FREE, series.filename_fmt
    file_struct = {exist:0,path:STRARR(N_ELEMENTS(runstr))}
    filename_fmt = -1
    FOR j = 0, n_runs - 1 DO BEGIN
      IF (runstr[j] NE 'act') AND (runstr[j] NE 'dat') THEN BEGIN
        IF FILE_TEST(data_path+'parameters_'+runstr[j]) AND NOT $
          FILE_TEST(data_path+'parameters_'+runstr[j],/ZERO_LENGTH) THEN BEGIN

          file_struct.path[j] = data_path + 'parameters_' + runstr[j]
          filename_fmt = 1
        ENDIF ELSE IF FILE_TEST(data_path+runstr[j]+'_parameters') AND NOT $
          FILE_TEST(data_path+'_parameters'+runstr[j],/ZERO_LENGTH) THEN BEGIN

          file_struct.path[j] = data_path + runstr[j] + '_parameters'
          filename_fmt = 2
        ENDIF ELSE IF FILE_TEST(data_path+'parnew_' + runstr[j]) AND NOT $
          FILE_TEST(data_path+'parnew_'+runstr[j],/ZERO_LENGTH) THEN BEGIN

          file_struct.path[j] = data_path + 'parnew_' + runstr[j]
          filename_fmt = 3
        ENDIF ELSE IF FILE_TEST(data_path+'run_'+runstr[j]+$
          '/parameters.dat') AND NOT FILE_TEST(data_path+'run_'+runstr[j]+$
          '/parameters.dat',/ZERO_LENGTH) THEN BEGIN

          filename_fmt = 4
          file_struct.path[j] = data_path + 'run_' + runstr[j] + '/parameters.dat'
        ENDIF
      ENDIF ELSE IF FILE_TEST(data_path+'parameters.dat') AND NOT $
        FILE_TEST(data_path+'parameters.dat'+runstr[j],/ZERO_LENGTH) THEN BEGIN

        file_struct.path[j] = data_path + 'parameters.dat'
	filename_fmt = 0
      ENDIF ELSE IF FILE_TEST(data_path+'parnew.dat') AND NOT $
        FILE_TEST(data_path+'parnew.dat'+runstr[j],/ZERO_LENGTH) THEN BEGIN

        file_struct.path[j] = data_path + 'parnew.dat'
	filename_fmt = 5
      ENDIF
      IF j EQ 0 THEN filename_fmt_arr = filename_fmt ELSE $
        filename_fmt_arr = [filename_fmt_arr, filename_fmt]
    ENDFOR
    file_struct.exist = TOTAL(FILE_TEST(file_struct.path)) EQ n_runs
    series.filename_fmt = PTR_NEW(filename_fmt_arr)
  ENDIF ELSE BEGIN
    ; only use HDF5 for field, mom, vsp, checkpoint, geometry
    gfile = STRMID(par.geomfile,0,3)
    mgfile = STRMID(par.magn_geometry,0,3)
    IF (WHERE(['fie','mom','vsp','che','pro',gfile,mgfile] EQ $
      STRMID(name,0,3)))[0] NE -1 THEN $
      h5ext = par.write_h5 ? '.h5' : '' ELSE h5ext = ''

    file_struct = REPLICATE({exist:0,path:STRARR(N_ELEMENTS(runstr))},$
      N_ELEMENTS(name))

    FOR n = 0, N_ELEMENTS(name) - 1 DO BEGIN
      FOR j = 0, N_ELEMENTS(runstr) - 1 DO BEGIN
        CASE (*series.filename_fmt)[run[j]] OF
          -1 : file_struct[n].path[j] = ''
          0 : file_struct[n].path[j] = data_path + name[n] + '.dat' + h5ext
          1 : file_struct[n].path[j] = data_path + name[n] + '_' + runstr[j] + h5ext
          2 : file_struct[n].path[j] = data_path + runstr[j] + '_' + name[n] + h5ext
          3 : BEGIN ; for compatibility with old mom file format
      	    name_split = STRSPLIT(name[n],'_',/EXTRACT)
	    IF N_ELEMENTS(name_split) EQ 2 THEN $
	      file_struct[n].path[j] = data_path + name_split[0] + $
              name_split[1] + '_' + runstr[j] + h5ext $
              ELSE file_struct[n].path[j] = data_path + name[n] + '_' + runstr[j] + h5ext
	  END
          4 : file_struct[n].path[j] = data_path + 'run_' + runstr[j] + $
      	    '/' + name + '.dat' + h5ext
          5 : BEGIN ; for compatibility with gene10 file format
      	    name_split = STRSPLIT(name[n],'_',/EXTRACT)
	    IF N_ELEMENTS(name_split) EQ 2 THEN $
	      file_struct[n].path[j] = data_path + name_split[0] + $
              name_split[1] + '.dat' + h5ext $
              ELSE file_struct[n].path[j] = data_path + name[n] + '.dat' + h5ext
	  END
          ELSE : file_struct[n].path = data_path + name[n] + '.dat' + h5ext
        ENDCASE
      ENDFOR
      file_struct[n].exist = $
        (TOTAL(FILE_TEST(file_struct[n].path)) EQ n_runs) AND $
        (TOTAL(FILE_TEST(file_struct[n].path,/ZERO_LENGTH)) EQ 0)
    ENDFOR
  ENDELSE

  RETURN, file_struct

END

;##########################################################################

FUNCTION convert_run_string, run_string_in
; replace dash ranges with comma chains, return new run_string
; 'dat' or 'act' is the label of an unrenamed run (.dat files)

  COMMON global_vars

  comma_split = STRSPLIT(run_string_in,',',/EXTRACT)
  comma_pos = STRSPLIT(run_string_in,',',LENGTH=cspl_len)
  run_string_out = ''
  FOR i = 0, N_ELEMENTS(comma_split) - 1 DO BEGIN
    dash_split = STRSPLIT(comma_split[i],'-',/EXTRACT)

    ; '-' split: insert all number in range <number1>-<number2>
    IF N_ELEMENTS(dash_split) GT 1 THEN BEGIN
      comma_split[i] = dash_split[0]
      ch_end = STRMID(dash_split[1],0,1)
      byte_end = BYTE(STRMID(dash_split[1],0,1))
      IF (ch_end GE 'A' AND ch_end LE 'Z') OR $
        (ch_end GE 'a' AND ch_end LE 'z') THEN BEGIN

        ch_start = STRMID(dash_split[0],STRLEN(dash_split[0])-1,1)
        IF (ch_start GE 'A' AND ch_start LE 'Z') OR $
          (ch_start GE 'a' AND ch_start LE 'z') THEN BEGIN
	  base = STRMID(comma_split[i],0,STRLEN(comma_split[i])-1)
          byte_start = BYTE(ch_start) + 1B
        ENDIF ELSE BEGIN
          base = comma_split[i]
          byte_start = 97B
	ENDELSE
        FOR j = byte_start[0], byte_end[0] DO BEGIN
	  IF j EQ 91B THEN j = 97B ; jump from 'Z' to 'a'
          comma_split[i] += ','+base+STRING(BYTE(j))
	ENDFOR
      ENDIF ELSE BEGIN
        IF ch_end EQ '0' THEN BEGIN ; for cases with leading zeroes
          len = STRLEN(dash_split[0])
          fmt = "(I0" + rm0es(len) + ")"
          FOR j = FIX(dash_split[0]) + 1, FIX(dash_split[1]) DO $
            comma_split[i] += ',' + STRING(j,FORMAT=fmt)
        ENDIF ELSE BEGIN ; regular case (no leading zeroes)
          FOR j = FIX(dash_split[0]) + 1, FIX(dash_split[1]) DO $
            comma_split[i] += ',' + rm0es(j)
        ENDELSE
      ENDELSE
    ENDIF
    ; '=' split: insert numbers starting from <number1> until last number found
    IF STRPOS(comma_split[i],'=',STRLEN(comma_split[i])-1) GT 0 THEN BEGIN
      comma_split[i] = STRMID(comma_split[i],0, STRLEN(comma_split[i])-1) ; run_label without '='
      len = STRLEN(comma_split[i])
      fmt="(I0"+rm0es(len)+")" ;for input numbers starting with 000...
      j = FIX(comma_split[i])+1
      filename = gui.out.data_path + 'nrg_' + STRING(j,FORMAT=fmt)
      WHILE FILE_TEST(filename) DO BEGIN
        comma_split[i] += ',' + STRING(j,FORMAT=fmt)
        j += 1
        filename = gui.out.data_path + 'nrg_' + STRING(j,FORMAT=fmt)
      ENDWHILE
    ENDIF
    ; '+' split: add letters in alphabetical order to run string
    ; until files don't exist
    IF STRPOS(comma_split[i],'+',STRLEN(comma_split[i])-1) GT 0 THEN BEGIN
      base = STRMID(comma_split[i],0, STRLEN(comma_split[i])-1) ; run_label without '+'
      comma_split[i] = base ; first file without letter
      ch = 97B              ; 97B = 'a'
      filename = gui.out.data_path + 'nrg_' + base + STRING(ch)
      WHILE FILE_TEST(filename) DO BEGIN
        comma_split[i] += ',' + base + STRING(ch)
        ch = ch + 1B
        filename = gui.out.data_path + 'nrg_' + base + STRING(ch)
      ENDWHILE
    ENDIF

    run_string_out += comma_split[i] + $
      STRMID(run_string_in,comma_pos[i]+cspl_len[i],1)
  ENDFOR

  RETURN, run_string_out

END

;##########################################################################

PRO plot_legend, coltable, captions, per_line=per_line, x_offset=x_offset,$
  linestyles=linestyles, psym=psym, thick=thick, csize=csize, cthick=cthick
; coltable: intarray with color indices of the graphs
; captions: strarray with captions to the graphs
; per_line: number of captions per line
; x_offset: array with x offsets
; linestyle: array with linestyles for each line
; psym: array containing symbol definition for each line; default: 0
; thick: array with line thickness
; csize: array of char size
; cthick: array of char thickness

  IF NOT KEYWORD_SET(per_line) THEN per_line = 4
  IF NOT KEYWORD_SET(x_offset) THEN x_offset = [0.15,0.0,0.9,0.2]
  IF N_ELEMENTS(coltable) EQ 0 THEN BEGIN
    PRINT, 'error, expected call: plot_legend, coltable, captions'
    RETURN
  ENDIF
  IF N_ELEMENTS(coltable) NE N_ELEMENTS(captions) THEN BEGIN
    PRINT, 'error: coltable and captions have to be of same dimension'
    RETURN
  ENDIF
  IF NOT KEYWORD_SET(linestyles) THEN BEGIN
    linestyles = INTARR(N_ELEMENTS(coltable))
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(linestyles) NE N_ELEMENTS(coltable) THEN BEGIN
      PRINT, 'error: coltable and linestyles have to be of same dimension'
      RETURN
    ENDIF
  ENDELSE
  IF NOT KEYWORD_SET(psym) THEN BEGIN
    psym = INTARR(N_ELEMENTS(coltable))
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(psym) NE N_ELEMENTS(coltable) THEN BEGIN
      PRINT, 'error: coltable and psym have to be of same dimension'
      RETURN
    ENDIF
  ENDELSE
  IF NOT KEYWORD_SET(thick) THEN BEGIN
    thick = INTARR(N_ELEMENTS(coltable))
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(thick) NE N_ELEMENTS(coltable) THEN BEGIN
      PRINT, 'error: coltable and thick have to be of same dimension'
      RETURN
    ENDIF
  ENDELSE
  IF NOT KEYWORD_SET(csize) THEN BEGIN
    csize = INTARR(N_ELEMENTS(coltable))
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(csize) NE N_ELEMENTS(coltable) THEN BEGIN
      PRINT, 'error: coltable and csize have to be of same dimension'
      RETURN
    ENDIF
  ENDELSE
  IF NOT KEYWORD_SET(cthick) THEN BEGIN
    cthick = INTARR(N_ELEMENTS(coltable))
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(cthick) NE N_ELEMENTS(coltable) THEN BEGIN
      PRINT, 'error: coltable and cthick have to be of same dimension'
      RETURN
    ENDIF
  ENDELSE
  n_vars = N_ELEMENTS(coltable)

  width = 14 / per_line
  n_lines = n_vars / per_line

  IF (n_vars MOD per_line) NE 0 THEN n_lines = n_lines + 1

  emptyticks = STRARR(8)
  emptyticks[*] = ' '

  IF x_offset[1] LT 0.05 THEN x_offset[1] = 0.03

  PLOT, [0,14], [0,n_lines+1], /NODATA, TICKLEN=0.0, $
    XTICKNAME=emptyticks, YTICKNAME=emptyticks,$
    YMARGIN=[0,0], /XSTYLE, /YSTYLE, COLOR=1, $
    POSITION=x_offset, /NOERASE

  FOR var = 0, n_vars - 1 DO BEGIN
    linewidth = width * (var MOD per_line)
    y_pos = n_lines - FIX(var/per_line)

    PLOTS, [linewidth+1,linewidth+2], [y_pos,y_pos], COLOR=coltable[var], $
      linestyle = linestyles[var], psym=psym[var], thick=thick[var]

    XYOUTS, 2.2 + linewidth, y_pos, captions[var], COLOR=1,$
      charsize=csize[var],charthick=cthick[var]
  ENDFOR

END

;##########################################################################

FUNCTION contour_levels, data, n_levels, log=log, no_fix_zero=no_fix_zero,$
                         img=img,min=min,max=max
; computes the levels and colors for contour plots
; returns level/color array
; log: logarithmic color distribution
; no_fix_zero : suppress fixing zero to green color
; img : return color for each data point (e.g., for png img plotting)

  IF N_ELEMENTS(data) LT 2 THEN BEGIN
    PRINT, 'contour_levels error: wrong data format'
    RETURN, [[0.0,1.0],[0.0,1.0]] ; black and white
  ENDIF
  IF N_ELEMENTS(n_levels) EQ 0 THEN BEGIN
     IF NOT KEYWORD_SET(img) THEN n_levels = 20 $
     ELSE n_levels = 255
  ENDIF
  IF n_levels LT 2 THEN n_levels = 20
  IF NOT KEYWORD_SET(log) THEN log = 0
  IF NOT KEYWORD_SET(min) THEN min = MIN(data)
  IF NOT KEYWORD_SET(max) THEN max = MAX(data)

  IF max - min LT 1e-14 THEN BEGIN
    PRINT, 'contour_levels warning: data (near) constant'
    max = max + min + 1e-12
  ENDIF
  absmax = ABS(max) > ABS(min)

  IF NOT KEYWORD_SET(img) THEN BEGIN
     levels = FLTARR(n_levels)
     colors = INTARR(n_levels)
  ENDIF ELSE BEGIN
     levels = data
     colors = INTARR(N_ELEMENTS(data))
  ENDELSE

  IF log EQ 1 THEN BEGIN
    IF min LT 0 THEN BEGIN
      PRINT, 'contour_levels warning: log scale only for positive values'
      log = 0
    ENDIF
    IF min EQ 0 THEN min = 1e-5 ;absmax / n_levels
    min = ALOG10(min)
    max = ALOG10(max)
  ENDIF

  IF NOT KEYWORD_SET(img) THEN BEGIN
     levels = min + (max - min) * INDGEN(n_levels) / FLOAT(n_levels-1)
     levels[0] = levels[0] - ABS(min) * 0.001 ; suppress out-of-range
     levels[n_levels-1] = levels[n_levels-1] + ABS(max) * 0.001
  ENDIF

  IF log EQ 0 THEN BEGIN     ; linear scale
    IF NOT KEYWORD_SET(no_fix_zero) THEN BEGIN
     IF min * max LT 0  THEN BEGIN ; green is fixed to zero
      IF KEYWORD_SET(img) THEN BEGIN
         FOR j = 0L, N_ELEMENTS(colors)-1 DO BEGIN
            IF levels[j] LE 0 THEN $
               colors[j] = 127.0 * (1.0 + levels[j]>min / ABS(min)) $
            ELSE $
               colors[j] = 127.0 * (1.0 + levels[j]<max / ABS(max))
         ENDFOR
      ENDIF ELSE BEGIN
         temp = MIN(ABS(levels),ind_zero)
         FOR j = 0L, ind_zero DO $
            colors[j] = 127.0 * (1.0 + levels[j] / ABS(min))
         FOR j = 1L + ind_zero, N_ELEMENTS(colors) - 1L DO $
            colors[j] = 127.0 * (1.0 + levels[j] / ABS(max))
      ENDELSE
     ENDIF ELSE BEGIN          ; fully negative or positive scale
      IF min GE 0 THEN colors = levels / absmax * 255.0 ELSE $
        colors = - levels / absmax * 255.0
     ENDELSE
    ENDIF ELSE BEGIN
       IF NOT KEYWORD_SET(img) THEN $
          colors = INDGEN(n_levels)*255.0/n_levels+1 $
       ELSE $
          colors = ROUND(n_levels/(max-min)*(levels-min))
    ENDELSE
  ENDIF ELSE BEGIN     ; logarithmic scale
    colors = INDGEN(n_levels)*255.0/n_levels+1 ; avoid color=0
    levels = 10^levels
    IF KEYWORD_SET(img) THEN BEGIN
       print, 'logarithmic scale not implemented for img option'
       STOP
    ENDIF
  ENDELSE

  IF NOT KEYWORD_SET(img) THEN BEGIN
     level_colors = FLTARR(n_levels,2)
     level_colors[*,0] = levels[*]
     level_colors[*,1] = colors[*] > 0.0
  ENDIF ELSE BEGIN
     level_colors = TEMPORARY(colors)
  ENDELSE

  RETURN, level_colors

END

;##########################################################################

PRO plot_colorbar, lev_col, position=pos, orientation=orientation,$
  charsize=charsize, prec=prec, log=log, xmargin=xmargin, ymargin=ymargin
; lev_col: contour levels
; position: array with relative position coordinates
;   2 elements: x positions for orientation=0, y position for orientation=1
;   4 elements: x and y positions
; orientation: 0 = vertical, 1 = horizontal (default: automatic)
; charsize: char size of color bar axis labels
; prec: precision of color bar axis labels
; log: logarithmic color distribution (should agree with log switch in
;   contour_levels)
; xmargin: graphics keyword (see IDL help)
; ymargin: graphics keyword (see IDL help)

  IF N_ELEMENTS(orientation) NE 1 THEN BEGIN
    IF N_ELEMENTS(pos) EQ 4 THEN orientation = $
      (pos[3] - pos[1]) / (pos[2] - pos[0]) LT 1 ELSE orientation = 0
  ENDIF
  IF N_ELEMENTS(pos) NE 4 THEN BEGIN
    IF N_ELEMENTS(pos) EQ 2 THEN BEGIN
      pos = orientation ? [!X.WINDOW[0],pos[0],!X.WINDOW[1],pos[1]] : $
        [pos[0],!Y.WINDOW[0],pos[1],!Y.WINDOW[1]]
    ENDIF ELSE BEGIN
      pos = orientation ? [!X.WINDOW[0],0.95,!X.WINDOW[1],0.97] : $
        [0.95,!Y.WINDOW[0],0.97,!Y.WINDOW[1]]

      IF TOTAL(!P.MULTI[1:2]) GT 2 THEN BEGIN ; 2+ color bars per page
        hor_scale = 1.0 ; 1.0 / !P.MULTI[1]
        vert_scale = 1.0 / !P.MULTI[2]
        pos *= [hor_scale,vert_scale,hor_scale,vert_scale]

        hor_shift = 0.0 ; (!P.MULTI[1] - 1 - !P.MULTI[0] MOD !P.MULTI[1]) / FLOAT(!P.MULTI[1])
        vert_shift = (!P.MULTI[0] / !P.MULTI[1]) / FLOAT(!P.MULTI[2])
        pos += [hor_shift,vert_shift,hor_shift,vert_shift]
      ENDIF
    ENDELSE
  ENDIF
  IF NOT KEYWORD_SET(charsize) THEN charsize = 0.707
  IF NOT KEYWORD_SET(prec) THEN prec = 1
  IF NOT KEYWORD_SET(log) THEN log = 0
  IF N_ELEMENTS(xmargin) NE 2 THEN xmargin = [10,3] ; IDL default
  IF N_ELEMENTS(ymargin) NE 2 THEN ymargin = [4,2] ; IDL default

  c_levels = N_ELEMENTS(lev_col[*,0])

  IF NOT orientation THEN BEGIN ; vertical
    color_bar = FLTARR(2,c_levels)
    FOR j = 0, c_levels - 1 DO color_bar[*,j] = lev_col[j,0]

    IF TOTAL(pos) NE 0 THEN BEGIN
      CONTOUR, color_bar, [0,1], lev_col[*,0], $
        LEVELS=lev_col[*,0], CHARSIZE=charsize, $
        C_COLORS=lev_col[*,1], YLOG=log, $
        YRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], $
        XSTYLE=5, /FILL, /NOERASE, POSITION=pos, YSTYLE=5, $
        XMARGIN=xmargin, YMARGIN=ymargin

      store_colors, store_rgb, new_ct=41
      PLOT, [0,1], lev_col[*,0], COLOR=1, YLOG=log, $
        YRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], XMARGIN=xmargin, $
        YMARGIN=ymargin, /XSTYLE, YSTYLE=5, CHARSIZE=charsize, $
        /NODATA, /NOERASE, XTICKLEN=0, XTICKNAME=[' ',' '], $
        YTICKLEN=0, XTICKS=1, XMINOR=1, YTICKS=0, YMINOR=1, $
        POSITION=pos
      AXIS, YAXIS=0, YTICKFORMAT='(A1)', COLOR=1, /YSTYLE
      AXIS, YAXIS=1, COLOR=1, /YSTYLE, CHARSIZE=charsize
    ENDIF ELSE BEGIN
      CONTOUR, color_bar, [0,1], lev_col[*,0], $
        LEVELS=lev_col[*,0], CHARSIZE=charsize, $
        C_COLORS=lev_col[*,1], YLOG=log, $
        YRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], $
        XSTYLE=5, /FILL, /NOERASE, YSTYLE=5, $
        XMARGIN=xmargin, YMARGIN=ymargin

      store_colors, store_rgb, new_ct=41
      PLOT, [0,1], lev_col[*,0], COLOR=1, YLOG=log, $
        YRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], XMARGIN=xmargin, $
        YMARGIN=ymargin, /XSTYLE, /YSTYLE, CHARSIZE=charsize, $
        /NODATA, /NOERASE, XTICKLEN=0, XTICKNAME=[' ',' '], $
        YTICKLEN=0, XTICKS=1, XMINOR=1, YTICKS=0, YMINOR=1
    ENDELSE
  ENDIF ELSE BEGIN ; horizontal
    color_bar = FLTARR(c_levels,2)
    FOR j = 0, c_levels - 1 DO color_bar[j,*] = lev_col[j,0]

    IF TOTAL(pos) NE 0 THEN BEGIN
      CONTOUR, color_bar, lev_col[*,0], [0,1], $
        LEVELS=lev_col[*,0], CHARSIZE=charsize, $
        C_COLORS=lev_col[*,1], XLOG=log, $
        XRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], $
        XSTYLE=5, /FILL, /NOERASE, POSITION=pos, YSTYLE=5, $
        XMARGIN=xmargin, YMARGIN=ymargin

      store_colors, store_rgb, new_ct=41
      PLOT, lev_col[*,0], [0,1], COLOR=1, XLOG=log, $
        XRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], XMARGIN=xmargin, $
        YMARGIN=ymargin, /XSTYLE, /YSTYLE, CHARSIZE=charsize, $
        /NODATA, /NOERASE, XTICKLEN=0, YTICKNAME=[' ',' '], $
        YTICKLEN=0, XTICKS=0, XMINOR=1, YTICKS=1, YMINOR=1, $
        POSITION=pos
    ENDIF ELSE BEGIN
      CONTOUR, color_bar, lev_col[*,0], [0,1], $
        LEVELS=lev_col[*,0], CHARSIZE=charsize, $
        C_COLORS=lev_col[*,1], XLOG=log, $
        XRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], $
        XSTYLE=5, /FILL, /NOERASE, YSTYLE=5, $
        XMARGIN=xmargin, YMARGIN=ymargin

      store_colors, store_rgb, new_ct=41
      PLOT, lev_col[*,0], [0,1], COLOR=1, XLOG=log, $
        XRANGE=[lev_col[0,0],lev_col[c_levels-1,0]], XMARGIN=xmargin, $
        YMARGIN=ymargin, /XSTYLE, /YSTYLE, CHARSIZE=charsize, $
        /NODATA, /NOERASE, XTICKLEN=0, YTICKNAME=[' ',' '], $
        YTICKLEN=0, XTICKS=0, XMINOR=1, YTICKS=1, YMINOR=1
    ENDELSE
  ENDELSE

  store_colors, store_rgb, /restore

END

;##########################################################################

PRO store_colors, store_rgb, restore=restore, new_ct=new_ct
; stores the current color set unter store_rgb and loads a new color table
; if new_ct is specified; call with /restore to revert to the previous set

  IF NOT KEYWORD_SET(restore) THEN BEGIN
    ; save current colors
    COMMON COLORS, R_ORIG, G_ORIG, B_ORIG, R_CURR, G_CURR, B_CURR
    store_rgb = [[R_CURR],[G_CURR],[B_CURR]]

    ; load new color table
    IF N_ELEMENTS(new_ct) EQ 1 THEN $
      LOADCT, new_ct, FILE='internal/colortable.tbl'

    RETURN
  ENDIF

  IF N_ELEMENTS(store_rgb) NE 768 THEN BEGIN
    PRINT, 'store_colors error: invalid store_rgb format'
    RETURN
  ENDIF

  ; restore old colors
  TVLCT, store_rgb[*,0], store_rgb[*,1], store_rgb[*,2]
  R_CURR = store_rgb[*,0]
  G_CURR = store_rgb[*,1]
  B_CURR = store_rgb[*,2]

END

;##########################################################################

FUNCTION vector_plot_arrows, field_x, field_y, comp_deriv=comp_deriv, $
  arrow_seed=arrow_seed, lx=lx, ly=ly
; Creates arrows for a 2D vector plot using two 2D fields for the x
; and y directions. The arrow seed field is an equidistant grid. The
; default for the box is par.lx by par.ly.
; Instead of providing the x and y fields, one can specify only a scalar
; field. In this case, comp_deriv=1 produces the gradient field from it,
; while comp_deriv=2 produces the perpendicular vector product field.
; In those cases, one can specify dx and dy to obtain physical
; derivatives (default: dx = dy = 1).
; The output can be plotted as follows:
; 1. create (reduced) plot box, e.g.:
;   lxm = 2 * (x_res - 1.0) / x_res - 1
;   lym = 2 * (y_res - 1.0) / y_res - 1
;   PLOT, lx * 0.5 * [-1,lxm,lxm,-1,-1], $
;     ly * 0.5 * [-1,-1,lym,lym,-1], /XSTYLE, /YSTYLE
; 2. plot arrows:
; FOR v = 0, arrow_seed[0] * arrow_seed[1] - 1 DO $
;   PLOTS, arrows[v,*,0], arrows[v,*,1], $
;   CLIP=[-0.5*lx,-0.5*ly,0.5*lx*lxm,0.5*ly*lym], NOCLIP=0

  COMMON global_vars

  IF N_ELEMENTS(arrow_seed) NE 2 THEN arrow_seed = [10,10]
  n_arrows = arrow_seed[0] * arrow_seed[1]

  lx = N_ELEMENTS(lx) GT 0 ? lx : par.lx
  ly = N_ELEMENTS(ly) GT 0 ? ly : par.ly

  IF NOT KEYWORD_SET(comp_deriv) THEN BEGIN
    IF (N_ELEMENTS(field_x) LT 1) OR (N_ELEMENTS(field_y) LT 1) THEN BEGIN
      PRINT, 'vector_plot_arrows error: no valid fields'
      RETURN, -1
    ENDIF
    res_fx = SIZE(field_x)
    res_fy = SIZE(field_y)
    IF (res_fx[0] LT 2) OR (res_fy[0] LT 2) THEN BEGIN
      PRINT, 'vector_plot_arrows error: too few dimensions'
      RETURN, -2
    ENDIF
    x_res = res_fx[1]
    y_res = res_fx[2]
    IF (x_res NE res_fy[1]) OR (y_res NE res_fy[2]) THEN BEGIN
      PRINT, 'vector_plot_arrows error: incompatible dimensions'
      RETURN, -3
    ENDIF
  ENDIF ELSE BEGIN
    res_f = SIZE(field_x)
    IF res_f[0] LT 2 THEN BEGIN
      PRINT, 'vector_plot_arrows error: too few dimensions'
      RETURN, -4
    ENDIF
    x_res = res_f[1]
    y_res = res_f[2]

    ; add boundary points for derivatives
    field_extd = FLTARR(x_res+2,y_res+2,/NOZERO)
    field_extd[1,1] = field_x
    field_extd[*,0] = [field_x[x_res-1,y_res-1],field_x[*,y_res-1],field_x[0,y_res-1]]
    field_extd[*,y_res+1] = [field_x[x_res-1,0],field_x[*,0],field_x[0,0]]
    field_extd[0,*] = [field_x[x_res-1,y_res-1],REFORM(field_x[x_res-1,*]),field_x[x_res-1,0]]
    field_extd[x_res+1,*] = [field_x[0,y_res-1],REFORM(field_x[0,*]),field_x[0,0]]

    dxt2_inv = 0.5 * x_res / lx
    dyt2_inv = 0.5 * y_res / ly

    IF comp_deriv EQ 1 THEN BEGIN ; gradient
      field_x[0,0] = dxt2_inv * $
        (field_extd[2:x_res+1,1:y_res] - field_extd[0:x_res-1,1:y_res])
      field_y = dyt2_inv * $
        (field_extd[1:x_res,2:y_res+1] - field_extd[1:x_res,0:y_res-1])
    ENDIF ELSE BEGIN ; vector product
      field_x[0,0] = dyt2_inv * $
        (field_extd[1:x_res,2:y_res+1] - field_extd[1:x_res,0:y_res-1])
      field_y = - dxt2_inv * $
        (field_extd[2:x_res+1,1:y_res] - field_extd[0:x_res-1,1:y_res])
    ENDELSE

    field_extd = 0
  ENDELSE

  ; number of segments per arrow
  n_segments = 8
  ; maximum arrow length (for a box width of 1.0)
  max_length_norm = 0.1

  max_length = SQRT(MAX(field_x^2)+MAX(field_y^2)) / (lx + ly) * 2.0
  seg_norm_fac_x = max_length_norm / (max_length * n_segments * lx)
  seg_norm_fac_y = max_length_norm / (max_length * n_segments * ly)

  ; regular seed field
  x_pos = REFORM(REBIN((FINDGEN(arrow_seed[0])+0.5)/arrow_seed[0],$
    [arrow_seed[0],arrow_seed[1]]),n_arrows)
  y_pos = REFORM(REBIN(REFORM((FINDGEN(arrow_seed[1])+0.5)/arrow_seed[1],$
    [1,arrow_seed[1]]),[arrow_seed[0],arrow_seed[1]]),n_arrows)

  ; store seed positions
  arrows = FLTARR(n_arrows,n_segments+3,2,/NOZERO)
  arrows[0,0,0] = x_pos
  arrows[0,0,1] = y_pos

  ; interpolate the field linearly to advance the
  ; arrays segment by segment
  FOR s = 1, n_segments - 1 DO BEGIN
    x_pos[0] = (x_res - 1) * arrows[*,s-1,0]
    y_pos[0] = (y_res - 1) * arrows[*,s-1,1]

    x0 = FLOOR(x_pos)
    x1 = x0 + 1
    y0 = FLOOR(y_pos)
    y1 = y0 + 1

    dx0 = x_pos - x0
    dx1 = 1.0 - dx0
    dy0 = y_pos - y0
    dy1 = 1.0 - dy0

    x_field_interp = dx1 * dy1 * field_x[x0,y0] + $
                     dx0 * dy1 * field_x[x1,y0] + $
                     dx1 * dy0 * field_x[x0,y1] + $
                     dx0 * dy0 * field_x[x1,y1]
    y_field_interp = dx1 * dy1 * field_y[x0,y0] + $
                     dx0 * dy1 * field_y[x1,y0] + $
                     dx1 * dy0 * field_y[x0,y1] + $
                     dx0 * dy0 * field_y[x1,y1]

    arrows[0,s,0] = arrows[*,s-1,0] + x_field_interp * seg_norm_fac_x
    arrows[0,s,1] = arrows[*,s-1,1] + y_field_interp * seg_norm_fac_y
  ENDFOR

  ; create arrowhead
  ; angle between the last segment and the arrowhead side
  arhead_angle = 0.125 * !PI
  ; length of the arrowhead side in units of the segment length
  arheadlen_o_seglen = 2.0
  sin_angle = SIN(arhead_angle)
  cos_angle = COS(arhead_angle)

  dx_norm = arheadlen_o_seglen * $
    (arrows[*,n_segments-1,0] - arrows[*,n_segments-2,0])
  dy_norm = arheadlen_o_seglen * $
    (arrows[*,n_segments-1,1] - arrows[*,n_segments-2,1])

  arrows[0,n_segments,0] = arrows[*,n_segments-1,0] - $
    (cos_angle * dx_norm - sin_angle * dy_norm)
  arrows[0,n_segments,1] = arrows[*,n_segments-1,1] - $
    (cos_angle * dy_norm + sin_angle * dx_norm)
  arrows[0,n_segments+1,0] = arrows[*,n_segments-1,*]
  arrows[0,n_segments+2,0] = arrows[*,n_segments-1,0] - $
    (cos_angle * dx_norm + sin_angle * dy_norm)
  arrows[0,n_segments+2,1] = arrows[*,n_segments-1,1] - $
    (cos_angle * dy_norm - sin_angle * dx_norm)

  arrows = TEMPORARY(arrows) < 1.0 > 0.0

  ; renormalize arrows to physical box
  arrows[0,0,0] = lx * (arrows[*,*,0] - 0.5)
  arrows[0,0,1] = ly * (arrows[*,*,1] - 0.5)

  RETURN, arrows

END

;##########################################################################

FUNCTION set_Lref, droplist_status
; returns Lref renormalization factor with respect to
; the current droplist status

  COMMON global_vars

  IF gui.droplist.mref_base NE 0 THEN BEGIN          ;if droplist for m does not exist
    WIDGET_CONTROL, gui.droplist.mref_base, MAP = 1  ;make visible in case it was invisible
  ENDIF

  IF gui.droplist.Tref_base NE 0 THEN BEGIN          ;if droplist for T does not exist
    WIDGET_CONTROL, gui.droplist.Tref_base, MAP = 1  ;make visible in case it was invisible
  ENDIF

  Lref_arr = '!6L!Dref!N'
  droplist_arr = 1.0D ; the order of the array entries is: L_ref, R_0, omts, omns
                                ; caution: the number of entries is
                                ; adjusted dynamically
  IF par.magn_geometry EQ 's_alpha' THEN BEGIN
    droplist_arr = [droplist_arr,1.0D/par.major_R]
    Lref_arr = [Lref_arr,'R!D0!N']
  ENDIF ELSE IF ((par.magn_geometry EQ 'circular') OR $
                 (STRPOS(par.magn_geometry,'miller') GE 0)) THEN BEGIN
    droplist_arr = [droplist_arr,1.0D/par.major_R,1.0D/par.minor_r]
    Lref_arr = [Lref_arr,'R!D0!N','a']
  ENDIF

  FOR isp = 0, par.n_spec-1 DO BEGIN
    IF spec[isp].omt NE 0.0 THEN BEGIN     ;test if normalization is non-zero
     droplist_arr = [droplist_arr,spec[isp].omt]
     Lref_arr = [Lref_arr,'L!DT,'+spec[isp].name+'!N']
    ENDIF
    IF spec[isp].omn NE 0.0 THEN BEGIN
     droplist_arr = [droplist_arr,spec[isp].omn]
     Lref_arr = [Lref_arr,'L!Dn,'+spec[isp].name+'!N']
    ENDIF
  ENDFOR

  series.Lref_str = Lref_arr[droplist_status]

  RETURN, droplist_arr[droplist_status]

END

;##########################################################################

FUNCTION time_renorm, var

  COMMON global_vars

  v_ref = SQRT(series.Tref*series.Qref/series.mref)  ; T in eV

  CASE var OF
    0  : RETURN, par.prec_single ? $ ; time
           FLOAT(series.Lref/v_ref) : series.Lref / v_ref
    -1 : RETURN, par.prec_single ? $ ; gamma, omega
           FLOAT(v_ref/series.Lref) : v_ref / series.Lref
    ELSE : printerror, 'error in time_renorm'
  ENDCASE

  PRINT, 'series.Lref/v_ref=' + rm0es(series.Lref/v_ref)

END

;##########################################################################

FUNCTION length_renorm, var

  COMMON global_vars

  rho_ref = SQRT(series.mref*series.Tref) / series.Qref / series.Bref

  CASE var OF
    0  : RETURN, 1.0D / (rho_ref)    ; length
    -1 : RETURN, rho_ref             ; fourier modes
    ELSE : printerror, 'error in length_renorm'
  ENDCASE

END

;##########################################################################

PRO printerror, errormsg, popup=popup

  COMMON global_vars

  IF N_ELEMENTS(errormsg) LT 1 THEN $
    errormsg = 'an error occured, please contact gene@ipp.mpg.de'
  IF errormsg NE '' THEN PRINT, errormsg
  IF KEYWORD_SET(popup) THEN message = $
    DIALOG_MESSAGE(errormsg,/ERROR,DIALOG_PARENT=gui.window.main) ELSE BEGIN

    WIDGET_CONTROL, gui.text.message, SET_VALUE=errormsg
  ENDELSE

END

;############################################################################

PRO plot_geom_coeff, coeff, coeffstr, xstr=xstr, cut1=cut1, cut2=cut2, $
                     show_yz_conn=show_yz_conn

  COMMON global_vars

  IF NOT KEYWORD_SET(xstr) THEN xstr='z'

  csize = 2.0

  YMARGIN=[4,4]
  XMARGIN=[8,2]

  coeffmax = MAX(coeff,MIN=coeffmin)
  IF (coeffmax EQ coeffmin) THEN range = [0.0,1.1*coeffmax] $
    ELSE range=[coeffmin<0,coeffmax]

  IF SIZE(coeff,/N_DIMENSIONS) EQ 1 THEN BEGIN
    IF xstr EQ 'x' THEN BEGIN
      xtitle = 'x / ' + get_var_string(/rhostr)
      xaxis = - 0.5 * series.lx + INDGEN(par.nx0) * series.dx

      PLOT, xaxis, coeff, YRANGE=range, /YSTYLE, COLOR=1, $
        XTITLE=xtitle, YTITLE=coeffstr, /XSTYLE, CHARSIZE=csize, $
        XMARGIN=XMARGIN, YMARGIN=YMARGIN
    ENDIF ELSE IF xstr EQ 'z' THEN BEGIN
       xtitle = 'z'
       xaxis = (*par.z)
       xtickv = 0.5 * (INDGEN(4*par.n_pol+1) - 2 * par.n_pol)
       xtickname = STRARR(4*par.n_pol+1)
       xtickname = rm0es(xtickv) + '!7p!6'
       xtickname[2*par.n_pol] = '0'
       xtickv *= !PI
       IF par.n_pol GT 1 THEN xtickname[2*INDGEN(par.n_pol+2)+1] = ' '

       PLOT, xaxis, coeff, YRANGE=range, /YSTYLE, COLOR=1, $
         XTITLE=xtitle, YTITLE=coeffstr, /XSTYLE, XMINOR=2,$
         XRANGE=[-par.n_pol*!PI,par.n_pol*!PI], XTICKS=N_ELEMENTS(xtickv)-1, $
         XTICKV=xtickv, CHARSIZE=csize, XTICKNAME=xtickname,$
         XMARGIN=XMARGIN, YMARGIN=YMARGIN
    ENDIF
  ENDIF ELSE IF SIZE(coeff,/N_DIMENSIONS) EQ 2 THEN BEGIN
    IF NOT par.x_local THEN BEGIN
       xtitle = 'x / ' + get_var_string(/rhostr)
       xaxis = - 0.5 * series.lx + INDGEN(par.nx0) * series.dx
    ENDIF ELSE IF NOT par.y_local THEN BEGIN
       xtitle = 'y / ' + get_var_string(/rhostr)
       xaxis = - 0.5 * series.ly + INDGEN(par.ny0) * series.dy
    ENDIF       
    ytitle = 'z'
    ytickv = 0.5 * (INDGEN(4*par.n_pol+1) - 2 * par.n_pol)
    ytickname = STRARR(4*par.n_pol+1)
    ytickname = rm0es(ytickv) + '!7p!6'
    ytickname[2*par.n_pol] = '0'
    ytickv *= !PI

    SURFACE, coeff, xaxis, (*par.z), COLOR=1, CHARSIZE=csize, $
      XTITLE=xtitle, YTITLE=ytitle, /XSTYLE, /YSTYLE, ZTITLE=coeffstr, $
      MIN_VALUE=range[0], MAX_VALUE=range[1], YTICKS=N_ELEMENTS(ytickv)-1, $
      YRANGE=[-par.n_pol*!PI,par.n_pol*!PI], YMINOR=2, YTICKV=ytickv, $
      YTICKNAME=ytickname, XMARGIN=XMARGIN, YMARGIN=YMARGIN

    IF KEYWORD_SET(show_yz_conn) THEN BEGIN
       dat_fft = FFT(coeff,dimension=1)
       dat_ext_fft = COMPLEXARR(par.ny0,par.nz0)
       
       sign = 1.0               ;1.0/(par.n0_global) ;1 ;-1
       FOR j=0, par.ny0/2 DO BEGIN
          dat_ext_fft[j,*] = dat_fft[j,*]*EXP(-2.0*sign*!PI*COMPLEX(0.0,1.0)*j*par.n0_global*par.q0)
       ENDFOR
       FOR j=par.ny0/2+1,par.ny0-1 DO BEGIN
          dat_ext_fft[j,*] = dat_fft[j,*]*EXP(-2.0*sign*!PI*COMPLEX(0.0,1.0)*(j-par.ny0)*par.n0_global*par.q0)
       ENDFOR
       dat_ext_right = FLOAT(FFT(dat_ext_fft,dimension=1,/INVERSE))
       FOR j=0, par.ny0/2 DO BEGIN
          dat_ext_fft[j,*] = dat_fft[j,*]*EXP(2.0*sign*!PI*COMPLEX(0.0,1.0)*j*par.n0_global*par.q0)
       ENDFOR
       FOR j=par.ny0/2+1,par.ny0-1 DO BEGIN
          dat_ext_fft[j,*] = dat_fft[j,*]*EXP(2.0*sign*!PI*COMPLEX(0.0,1.0)*(j-par.ny0)*par.n0_global*par.q0)
       ENDFOR
       dat_ext_left = FLOAT(FFT(dat_ext_fft,dimension=1,/INVERSE))
       
       dat_all = FLTARR(3*par.ny0,3*par.nz0)
       dat_all[par.ny0:2*par.ny0-1,par.nz0:2*par.nz0-1] = coeff
       dat_all[2*par.ny0:3*par.ny0-1,par.nz0:2*par.nz0-1] = coeff
       dat_all[0:par.ny0-1,par.nz0:2*par.nz0-1] = coeff
       
       dat_all[par.ny0:2*par.ny0-1,2*par.nz0:3*par.nz0-1] = dat_ext_right
       dat_all[par.ny0:2*par.ny0-1,0:par.nz0-1] = dat_ext_left
       
       dat_all[0:par.ny0-1,2*par.nz0:3*par.nz0-1] = dat_ext_right
       dat_all[2*par.ny0:3*par.ny0-1,2*par.nz0:3*par.nz0-1] = dat_ext_right
       dat_all[0:par.ny0-1,0:par.nz0-1] = dat_ext_left
       dat_all[2*par.ny0:3*par.ny0-1,0:par.nz0-1] = dat_ext_left
       
       c_levels = 20
       lev_col = contour_levels([MIN(dat_all),MAX(dat_all)],c_levels)
       
       CONTOUR, dat_all,-par.ly+INDGEN(3*par.ny0)*par.ly/par.ny0, $
                -3.0*!PI+INDGEN(3*par.nz0)*par.dz,$
                LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL, $
                XSTYLE=1, YSTYLE=1, COLOR = 1,$
                XTITLE=xtitle, YTITLE=ytitle, title=coeffstr,$
                ISOTROPIC=x_as_y, $
                CHARSIZE=csize, CHARTHICK=cthick, /NORMAL ;$
       lstyle = 2
       OPLOT, [0.0,0.0],!Y.CRANGE, color=0, LINESTYLE=lstyle
       OPLOT, [par.ly,par.ly],!Y.CRANGE, color=0, LINESTYLE=lstyle
       
       OPLOT, !X.CRANGE, [-!PI,-!PI],color=0, LINESTYLE=lstyle
       OPLOT, !X.CRANGE, [!PI,!PI],color=0, LINESTYLE=lstyle
    ENDIF

    IF N_ELEMENTS(cut1) GT 0 THEN BEGIN
       FOR icut=0, N_ELEMENTS(cut1)-1 DO BEGIN
          IF (cut1[icut] EQ -1) THEN dat = REFORM(TOTAL(coeff[*,*],1))/par.ny0 $
          ELSE dat = REFORM(coeff[cut1[icut],*])
          PLOT, (*par.z), dat, CHARSIZE=csize, color=1,$
                XTITLE=ytitle, YTITLE=coeffstr, /XSTYLE, /YSTYLE, $
                XTICKS=N_ELEMENTS(ytickv)-1, TITLE='cut1 = '+rm0es(cut1[icut]),$
                XRANGE=[-3.0*par.n_pol*!PI,3.0*par.n_pol*!PI], XMINOR=2, $
                XTICKV=ytickv, XTICKNAME=ytickname, THICK = 2
           IF (cut1[icut] GE 0) THEN $
              OPLOT, -3.0*!PI+INDGEN(3*par.nz0)*par.dz, color = 2, $
                     [REFORM(dat_ext_left[cut1[icut],*]),dat,REFORM(dat_ext_right[cut1[icut],*])] $
           ELSE OPLOT, -3.0*!PI+INDGEN(3*par.nz0)*par.dz, [dat,dat,dat], color=2
       ENDFOR
    ENDIF
    IF N_ELEMENTS(cut2) GT 0 THEN BEGIN
       FOR icut=0, N_ELEMENTS(cut2)-1 DO BEGIN
          PLOT, xaxis, coeff[*,cut2[icut]], CHARSIZE=csize, color=1,$
                XRANGE=[-0.5*par.ly,1.5*par.ly],$
                XTITLE=xtitle, YTITLE=coeffstr, /XSTYLE, /YSTYLE, $
                TITLE='cut2 = '+rm0es(cut2[icut])
          OPLOT, xaxis+par.ly, coeff[*,cut2[icut]], color=1
       ENDFOR
    ENDIF
  ENDIF

END

;##########################################################################
FUNCTION get_ktheta_fac, zind
;ktheta_fac : conversion factor to be applied to k_y in order to
;  get poloidal wave number (assuming k_z << k_y)

  COMMON global_vars

  ga1 = (*series.geom).gxx * (*series.geom).gyy - $
        (*series.geom).gxy^2
  ;covariant metric elements
  g_yy = ((*series.geom).gxx * (*series.geom).gzz - $
          (*series.geom).gxz^2)*(*series.geom).jacobian^2
  g_yz = ((*series.geom).gxy * (*series.geom).gxz - $
          (*series.geom).gxx * (*series.geom).gyz)*$
         (*series.geom).jacobian^2
  g_zz = ga1*(*series.geom).jacobian^2
  
  IF (par.x_local) THEN BEGIN
    ;norm. e_theta
     norm_e_theta = sqrt(((*series.geom).q*(*series.geom).C_y)^2*$
       g_yy+2*(*series.geom).q*(*series.geom).C_y*g_yz+g_zz)

    ;k_theta_fac
     ktheta_fac = (*series.geom).q*(*series.geom).C_y/$
                  norm_e_theta


  ENDIF ELSE BEGIN
     norm_e_theta = FLTARR(par.nx0,par.nz0,/NOZERO)
     ktheta_fac = FLTARR(par.nx0,par.nz0,/NOZERO)
        
     FOR ii = 0, par.nx0 - 1 DO BEGIN
        norm_e_theta[ii,*] = sqrt(((*series.geom).q[ii]*$
          (*series.geom).C_y[ii])^2*g_yy[ii,*]+$
          2*(*series.geom).q[ii]*(*series.geom).C_y[ii]*$
          g_yz[ii,*]+g_zz[ii,*])
        ktheta_fac[ii,*] = (*series.geom).q[ii]*(*series.geom).C_y[ii]/$
          norm_e_theta[ii,*]
     ENDFOR
  ENDELSE
  
  IF N_ELEMENTS(zind) LT 1 THEN RETURN, ktheta_fac $
  ELSE IF (par.x_local) THEN RETURN, ktheta_fac[zind] $
  ELSE RETURN, ktheta_fac[*,zind]


END

;##########################################################################

FUNCTION get_k2_fac, zind
;k2_fac : conversion factor to be applied to k_y in order to get
;  wave number in direction perp. to B field and radial
;  direction; relevant k value for most cmp. with exp.!)
  COMMON global_vars
  k2_fac = sqrt(((*series.geom).gxx*(*series.geom).gyy-$
      ((*series.geom).gxy)^2)/(*series.geom).gxx) 
  IF N_ELEMENTS(zind) LT 1 THEN RETURN, k2_fac $
  ELSE IF (par.x_local) THEN RETURN, k2_fac[zind] $
  ELSE RETURN, k2_fac[*,zind]

END

;############################################################################

PRO show_geometry_coeffs

  COMMON global_vars

  IF gui.out.out_format[2] EQ 1 THEN multi = [0,1,1] $
  ELSE multi = [0,2+2*par.norm08,2+1*par.norm08]

  eps = 0
  show_yz_conn = ((NOT par.y_local) AND (multi[2] EQ 1))
  coltable = show_yz_conn*33
  set_output, 'geometry', eps=eps, xsize=25, ysize=17, multi=multi, $
              coltable=coltable

  plot_geom_coeff, (*series.geom).gxx, 'g!Uxx!N'
;    cut1=[-1,0,par.ny0/2,par.ny0-1], cut2=[0, par.nz0/2, par.nz0-1]
  plot_geom_coeff, (*series.geom).gxy, 'g!Uxy!N'


  IF par.norm08 THEN BEGIN
    plot_geom_coeff, (*series.geom).gxz, 'g!Uxz!N'
    plot_geom_coeff, (*series.geom).gyy, 'g!Uyy!N'
    plot_geom_coeff, (*series.geom).gyz, 'g!Uyz!N'
    plot_geom_coeff, (*series.geom).gzz, 'g!Uzz!N'
  ENDIF

  plot_geom_coeff, (*series.geom).Bfield * series.Bref, 'B(z)', show_yz_conn=show_yz_conn

  IF par.norm08 THEN BEGIN
    plot_geom_coeff, (*series.geom).dBdx, 'dBdx'
    plot_geom_coeff, (*series.geom).dBdy, 'dBdy'
    plot_geom_coeff, (*series.geom).dBdz, 'dBdz'
    plot_geom_coeff, (*series.geom).dBdz, 'dBdz'

    IF (multi[2] EQ 1) THEN BEGIN
      ga1 = (*series.geom).gxx * (*series.geom).gyy - $
        (*series.geom).gxy^2
      ga2 = (*series.geom).gxx * (*series.geom).gyz - $
        (*series.geom).gxy * (*series.geom).gxz
      ga3 = (*series.geom).gxy * (*series.geom).gyz - $
        (*series.geom).gyy * (*series.geom).gxz

      IF (STRTRIM(par.magn_geometry,2) eq 'gist') THEN BEGIN
         K_x = - (*series.geom).dBdy
         K_y = (*series.geom).dBdx
      ENDIF ELSE BEGIN
         K_x = - ((*series.geom).dBdy + ga2 / ga1 * (*series.geom).dBdz)
         K_y = ((*series.geom).dBdx-ga3 / ga1 * (*series.geom).dBdz)         
      ENDELSE

      K_z = ((*series.geom).dBdy+ga2 / ga1 * (*series.geom).dBdx)

      IF (par.x_local) THEN BEGIN
        K_x /= (*series.geom).C_xy
        K_y /= (*series.geom).C_xy
        K_z /= (*series.geom).C_xy
      ENDIF ELSE BEGIN
        FOR ii = 0, par.nx0 - 1 DO BEGIN
          K_x[ii,*] /= (*series.geom).C_xy[ii]
          K_y[ii,*] /= (*series.geom).C_xy[ii]
          K_z[ii,*] /= (*series.geom).C_xy[ii]
        ENDFOR
      ENDELSE

      plot_geom_coeff,K_x,'K_x', show_yz_conn=show_yz_conn
      plot_geom_coeff,K_y,'K_y', show_yz_conn=show_yz_conn
      plot_geom_coeff,K_z,'K_z'

      IF ((par.magn_geometry NE 's_alpha') AND $
          (par.magn_geometry NE 's_alpha_B')) THEN BEGIN
         ktheta_fac = get_ktheta_fac()
         k2_fac = get_k2_fac()
         plot_geom_coeff,ktheta_fac,'k!D!4h!6!N factor'
         plot_geom_coeff,k2_fac,'k!D2!6!N factor'
      ENDIF


; for debugging consistency checks
;      shat_test = (*series.geom).gxy[*,par.nz0-1]/(*series.geom).gxx[*,par.nz0-1]-$
;              (*series.geom).gxy[*,0]/(*series.geom).gxx[*,0]  
;      shat_test /= 2.0*!PI*(*series.geom).C_y*par.q0/par.x0
;      shat_test_avg = TOTAL(shat_test)/par.ny0
;      PLOT,- 0.5 * series.ly + INDGEN(par.ny0) * series.dy, shat_test, title='shat test',$
;         xtitle='y', color=1
;      OPLOT, !X.CRANGE, [par.shat,par.shat], color=2
;      OPLOT, !X.CRANGE, [shat_test_avg,shat_test_avg], color=3
;      print, 'shat/shat_test: ', par.shat/shat_test_avg
;
;      jac_control = 1./sqrt($
;        (*series.geom).gxx*(*series.geom).gyy*(*series.geom).gzz+$
;        (*series.geom).gxy*(*series.geom).gyz*(*series.geom).gxz+$
;        (*series.geom).gxz*(*series.geom).gxy*(*series.geom).gyz-$
;        (*series.geom).gxz*(*series.geom).gyy*(*series.geom).gxz-$
;        (*series.geom).gyz*(*series.geom).gyz*(*series.geom).gxx-$
;        (*series.geom).gzz*(*series.geom).gxy*(*series.geom).gxy)
;      plot_geom_coeff,jac_control,'jac control'
    ENDIF
  ENDIF

  plot_geom_coeff, (*series.geom).jacobian, 'jacobian'

  IF NOT par.x_local THEN BEGIN
    plot_geom_coeff, (*series.geom).q, 'q', xstr='x'
    IF (multi[2] EQ 1) THEN BEGIN
      x_a = par.x0 + par.rhostar * (INDGEN(par.nx0) * par.dx - 0.5 * par.lx)
      shat = x_a / (*series.geom).q * DERIV(x_a,(*series.geom).q)
      plot_geom_coeff, shat, '!S!5s!R!7!U^!N', xstr='x'

      set_output, 'geometry', header=['x/a','q','shat'],$
                  dat=[[x_a],[(*series.geom).q],[shat]]

      IF (par.magn_geometry NE 'circular') THEN BEGIN
        PLOT, (*series.geom).R[par.nx0-1,*], (*series.geom).Z[par.nx0-1,*], $
          /ISOTROPIC, COLOR=1, /XSTYLE, /YSTYLE, XTITLE='!6R', YTITLE='!6Z',$
          YMARGIN=[4,4],XMARGIN=[8,2]
        FOR ix = 0, par.nx0-2 DO $
           OPLOT, (*series.geom).R[ix,*], (*series.geom).Z[ix,*], COLOR=ix+2
      ENDIF
    ENDIF
  ENDIF ELSE BEGIN
    IF ((par.magn_geometry NE 's_alpha') AND $
      (par.magn_geometry NE 'circular')) THEN $
         IF (par.y_local) THEN PLOT, (*series.geom).R[*], (*series.geom).Z[*], $
          /ISOTROPIC, COLOR=1, /XSTYLE, /YSTYLE, XTITLE='!6R', YTITLE='!6Z',$
          YMARGIN=[4,4],XMARGIN=[8,2]
  ENDELSE

  XYOUTS, 0.01, 0.96, '!6geometry coefficients (' + par.magn_geometry + $
    ')', COLOR=1, CHARSIZE=1.75, /NORMAL, CHARTHICK=1.5

  set_output, 'geometry', /reset, eps=eps

END

;##########################################################################

PRO check_request_dat, request_dat
  COMMON global_vars

  IF TOTAL(request_dat) EQ 0 THEN RETURN

  IF par.in_data EQ 'sxky' THEN BEGIN
    ; --- request all sxky needed for other formats ---
    request_dat[WHERE((request_dat[*,0] GT 0) OR $
    (request_dat[*,1] GT 0) OR (request_dat[*,2] GT 0) $
      OR (request_dat[*,3] GT 0)),1] = 1
    ; use kxky if kxsy is requested but neither kxky nor sxsy
    ind_arr = WHERE((request_dat[*,2] GT 0) AND $
      (request_dat[*,0] EQ 0) AND (request_dat[*,3] EQ 0))
    IF ind_arr[0] NE -1 THEN request_dat[ind_arr,0] = 1
  ENDIF ELSE IF par.in_data EQ 'sykx' THEN BEGIN
    ; --- request all sxky needed for other formats ---
    request_dat[WHERE((request_dat[*,0] GT 0) OR $
    (request_dat[*,1] GT 0) OR (request_dat[*,2] GT 0) $
      OR (request_dat[*,3] GT 0)),2] = 1

    ; use sxsy if sxsy, sxky, kxky is requested
    ind_arr = WHERE((request_dat[*,0] GT 0) OR $
      (request_dat[*,1] GT 0) OR (request_dat[*,3] GT 0))
    IF ind_arr[0] NE -1 THEN request_dat[ind_arr,3] = 1

    ; use sxky if sxky, kxky is requested
    ind_arr = WHERE((request_dat[*,0] GT 0) OR $
      (request_dat[*,1] GT 0))
    IF ind_arr[0] NE -1 THEN request_dat[ind_arr,1] = 1    
  ENDIF ELSE BEGIN               ; --- gene11 local
    ; --- request all kxky needed for other formats ---
    request_dat[WHERE((request_dat[*,0] GT 0) OR $
    (request_dat[*,1] GT 0) OR (request_dat[*,2] GT 0) $
      OR (request_dat[*,3] GT 0)),0] = 1
  ENDELSE


END

;##########################################################################

PRO create_3darr_struct, datptr, request_dat
  COMMON global_vars

  ; 0 - kxky
  ind_arr = WHERE(request_dat[*,0] GT 0)
  IF ind_arr[0] NE -1 THEN FOR ivar = 0, N_ELEMENTS(ind_arr) - 1 DO $
     datptr[ind_arr[ivar]].kxky = par.prec_single ? $
       PTR_NEW(COMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO),/NO_COPY) : $
       PTR_NEW(DCOMPLEXARR(par.nkx0,par.nky0,par.nz0,/NOZERO),/NO_COPY)
  ; 1 - sxky
  ind_arr = WHERE(request_dat[*,1] GT 0)
  IF ind_arr[0] NE -1 THEN FOR ivar = 0, N_ELEMENTS(ind_arr) - 1 DO $
     datptr[ind_arr[ivar]].sxky = par.prec_single ? $
       PTR_NEW(COMPLEXARR(par.nx0,par.nky0,par.nz0,/NOZERO),/NO_COPY) : $
       PTR_NEW(DCOMPLEXARR(par.nx0,par.nky0,par.nz0,/NOZERO),/NO_COPY)
  ; 2 - kxsy
  ind_arr = WHERE(request_dat[*,2] GT 0)
  IF ind_arr[0] NE -1 THEN FOR ivar = 0, N_ELEMENTS(ind_arr) - 1 DO $
     datptr[ind_arr[ivar]].kxsy = par.prec_single ? $
       PTR_NEW(COMPLEXARR(par.nkx0/2+1,par.ny0,par.nz0,/NOZERO),/NO_COPY) : $
       PTR_NEW(DCOMPLEXARR(par.nkx0/2+1,par.ny0,par.nz0,/NOZERO),/NO_COPY)
  ; 3 - sxsy
  ind_arr = WHERE(request_dat[*,3] GT 0)
  IF ind_arr[0] NE -1 THEN FOR ivar = 0, N_ELEMENTS(ind_arr) - 1 DO $
     datptr[ind_arr[ivar]].sxsy = par.prec_single ? $
       PTR_NEW(FLTARR(par.nx0,par.ny0,par.nz0,/NOZERO),/NO_COPY) : $
       PTR_NEW(DBLARR(par.nx0,par.ny0,par.nz0,/NOZERO),/NO_COPY)

END


;##########################################################################

PRO calc_ffts_3darr, datptr, request_dat
; Calculates sxky, kxsy, sxsy from kxky in every time step for all
; requested fields, for all species.
; sxsy: standard via sxky, only if kxsy (not sxky) requested via kxsy
; NOTE: [0,0,0] on the lhs is used as a faster alternative to [*,*,*]

  COMMON global_vars

  nkx0o2 = par.nkx0 / 2

  FOR f = 0, N_ELEMENTS(request_dat[*,0]) - 1 DO BEGIN
     IF par.in_data EQ 'sxky' THEN BEGIN ; --- (input data is sxky)
       ; calculate kxky data
        IF request_dat[f,0] EQ 1 THEN (*datptr[f].kxky)[0,0,0] = $
           FFT((*datptr[f].sxky),DIMENSION=1,DOUBLE=(1-par.prec_single))

        ; calculate sxsy data
        IF request_dat[f,3] EQ 1 THEN BEGIN
           temp_data = DCOMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO)
           temp_data[0,0,0] = $
              (*datptr[f].sxky)[*,0:par.nky0-1,*]
           temp_data[*,par.nky0,*] = 0
           FOR j = par.nky0 + 1, 2 * par.nky0 - 1 DO temp_data[0,j,0] = $
              CONJ((*datptr[f].sxky)[*,2*par.nky0-j,*])
           (*datptr[f].sxsy)[0,0,0] = $
              FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=(1-par.prec_single))
        ENDIF
        
        ; calculate kxsy data
        IF request_dat[f,2] EQ 1 THEN BEGIN
           IF (request_dat[f,3] EQ 1 AND $
               request_dat[f,0] EQ 0) THEN BEGIN
              temp_data = FFT((*datptr[f].sxsy),DIMENSION=1,/INVERSE,$
                              DOUBLE=(1-par.prec_single))
              (*datptr[f].kxsy)[0,0,0] = temp_data[0:nkx0o2,*,*]
           ENDIF ELSE BEGIN
              IF request_dat[f,0] EQ 0 THEN (*datptr[f].kxky)[0,0,0] = $
                 FFT((*datptr[f].sxky),DIMENSION=1,DOUBLE=(1-par.prec_single))
              temp_data = DCOMPLEXARR(nkx0o2+1,2*par.nky0,par.nz0,/NOZERO)
              temp_data[0,0,0] = (*datptr[f].kxky)[0:nkx0o2,*,*]
              temp_data[0:nkx0o2,par.nky0,*] = 0
              temp_data[0,par.nky0+1,0] = $
                 CONJ((*datptr[f].kxky)[0,par.nky0-1-INDGEN(par.nky0-1),*])
              temp_data[1,par.nky0+1,0] = $
                 CONJ((*datptr[f].kxky)[nkx0o2+1+INDGEN(nkx0o2-1),$
                                       par.nky0-1-INDGEN(par.nky0-1),*])
              (*datptr[f].kxsy)[0,0,0] = $
                 FFT(temp_data,DIMENSION=2,/INVERSE,DOUBLE=(1-par.prec_single))
           ENDELSE
        ENDIF
     ENDIF ELSE IF par.in_data EQ 'sykx' THEN BEGIN 
        ; --- input data is sykx, i.e. kxsy is filled
        
        ; calculate sxsy data
        IF request_dat[f,3] EQ 1 THEN BEGIN
           temp_data = DCOMPLEXARR(par.nx0,par.ny0,par.nz0,/NOZERO)
           temp_data[0:nkx0o2,*,*] = (*datptr[f].kxsy)[0:nkx0o2,*,*]
           FOR i = nkx0o2 + 1, par.nkx0 - 1 DO temp_data[i,*,*] = $
              CONJ((*datptr[f].kxsy)[par.nx0-i,*,*])
          (*datptr[f].sxsy)[0,0,0] = $
             FFT(temp_data,DIMENSION=1,/INVERSE,DOUBLE=(1-par.prec_single))
       ENDIF
        
       ; calculate sxky data
        IF request_dat[f,1] EQ 1 THEN BEGIN
           temp_data = FFT((*datptr[f].sxsy),DIMENSION=2,DOUBLE=(1-par.prec_single))
          (*datptr[f].sxky)[0,0,0] = temp_data[*,0:par.nky0-1,*]
        ENDIF

        ; calculate kxky data
        IF request_dat[f,0] EQ 1 THEN $
          (*datptr[f].kxky)[0,0,0] = $
          FFT((*datptr[f].sxky),DIMENSION=1,DOUBLE=(1-par.prec_single))
      ENDIF ELSE BEGIN           ; --- gene11 (output data is kxky)
;        IF request_dat[f,0] THEN BEGIN ;seed freq. of 80kHz for testing (occasionally needed by tbg)
;           (*datptr[f].kxky)[*,*,*] = 0.0 ;0.1*SIN(2.0*!PI*80E3*par.Lref/par.cref*mom_time)
;           (*datptr[f].kxky)[0,2,*] = 10.0*SIN(2.0*!PI*80E3*par.Lref/par.cref*mom_time)
;        ENDIF
        ; Doppler correction
        IF request_dat[f,0] AND (series.doppler_corr_mom NE 0.0) THEN BEGIN
          IF mom_time GT par.ExB_stime THEN BEGIN
            FOR iky = 0, par.nky0 - 1 DO (*datptr[f].kxky)[*,iky,*] *= $
               EXP(-DCOMPLEX(0.0,1.0)*(*par.ky)[iky]*(*series.geom).C_y*$
                   par.Lref/par.rhoref*series.doppler_corr_mom*$
                   (mom_time-par.ExB_stime))
          ENDIF
        ENDIF

        ; calculate sxky data
        IF request_dat[f,1] THEN (*datptr[f].sxky)[0,0,0] = $
          FFT((*datptr[f].kxky),DIMENSION=1,/INVERSE,DOUBLE=(1-par.prec_single))

        ; calculate kxsy data
        IF request_dat[f,2] THEN BEGIN
          temp_data = par.prec_single ? $
            COMPLEXARR(nkx0o2+1,2*par.nky0,par.nz0,/NOZERO) : $
            DCOMPLEXARR(nkx0o2+1,2*par.nky0,par.nz0,/NOZERO)

          temp_data[0,0,0] = (*datptr[f].kxky)[0:nkx0o2,*,*]
          temp_data[0,par.nky0,0] = COMPLEXARR(nkx0o2+1,1,par.nz0)
          temp_data[0,par.nky0+1,0] = $
            CONJ((*datptr[f].kxky)[0,par.nky0-1-INDGEN(par.nky0-1),*])
          temp_data[1,par.nky0+1,0] = $
            CONJ((*datptr[f].kxky)[par.nkx0-1-INDGEN(par.nkx0-nkx0o2-1),$
            par.nky0-1-INDGEN(par.nky0-1),*])

          (*datptr[f].kxsy) = TEMPORARY(FFT(temp_data,DIMENSION=2,$
            /INVERSE,DOUBLE=(1-par.prec_single),/OVERWRITE))
        ENDIF

        ; calculate sxsy data
        IF request_dat[f,3] THEN BEGIN
          IF request_dat[f,2] AND NOT request_dat[f,1] THEN BEGIN
            temp_data = par.prec_single ? $
              COMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO) : $
              DCOMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO)

            temp_data[0,0,0] = (*datptr[f].kxsy)[0:nkx0o2-1,*,*]
            temp_data[nkx0o2,0,0] = COMPLEXARR(1,2*par.nky0,par.nz0)
            temp_data[nkx0o2+1,0,0] = $
              CONJ((*datptr[f].kxsy)[nkx0o2-1-INDGEN(nkx0o2-1),*,*])

            *datptr[f].sxsy = TEMPORARY(FFT(temp_data,$
              DIMENSION=1,/INVERSE,/OVERWRITE))
          ENDIF ELSE BEGIN
            IF NOT request_dat[f,1] THEN BEGIN
              temp_sxky = FFT(*datptr[f].kxky,DIMENSION=1,/INVERSE,$
                DOUBLE=(1-par.prec_single))

              temp_data = par.prec_single ? $
                COMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO) : $
                DCOMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO)

              temp_data[0,0,0] = temp_sxky
              temp_data[0,par.nky0,0] = COMPLEXARR(par.nkx0,1,par.nz0)
              temp_data[0,par.nky0+1,0] = $
                CONJ(TEMPORARY(temp_sxky[*,par.nky0-1-INDGEN(par.nky0-1),*]))
            ENDIF ELSE BEGIN
              temp_data = par.prec_single ? $
                COMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO) : $
                DCOMPLEXARR(par.nkx0,2*par.nky0,par.nz0,/NOZERO)

              temp_data[0,0,0] = *datptr[f].sxky
              temp_data[0,par.nky0,0] = COMPLEXARR(par.nkx0,1,par.nz0)
              temp_data[0,par.nky0+1,0] = $
                CONJ((*datptr[f].sxky)[*,par.nky0-1-INDGEN(par.nky0-1),*])
            ENDELSE

            (*datptr[f].sxsy)[0,0,0] = TEMPORARY(FFT(temp_data,DIMENSION=2,$
              /INVERSE,DOUBLE=(1-par.prec_single),/OVERWRITE))
          ENDELSE
        ENDIF
      ENDELSE
   ENDFOR
END


;##########################################################################

FUNCTION time_fft, time, arr, interp_fac=interp_fac, width=width, $
  overlap=overlap, shift=myshift, filter_func=filter_func, omega=omega
; returns abs. value of windowed fft of time series
; time : array containing time values
; arr  : array with time dependent information
; optional arguments:
; interp_fac (scalar): use interp_fac*(number of input time steps) for
;                      interpolation (default: 2.0)
; width (scalar):      relative width of time windows (default: 0.25)
; overlap (scalar):    relative overlap of time windows (default: 0.5)
; shift (logical):     the usual FFT output [0,pos.freq.,neg. freq.]
;                      will be shifted to [neg. freq.,0,pos. freq]
; filter_func:         filter funtion applied on each time window (default: HANNING)
; omega:               returns array with omega values

  IF N_ELEMENTS(interp_fac) NE 1 THEN interp_fac = 2.0
  IF N_ELEMENTS(width) NE 1 THEN width = 0.25
  IF N_ELEMENTS(overlap) NE 1 THEN overlap = 0.5
  IF N_ELEMENTS(filter_func) NE 1 THEN filter_func = 'HANNING'
  IF N_ELEMENTS(myshift) NE 1 THEN myshift = 0

  arr = REFORM(arr)
  IF (N_ELEMENTS(time) NE N_ELEMENTS(arr)) THEN BEGIN
    printerror, 'error in time_fft function: array sizes are not conform'
    RETURN, -1
  ENDIF

  ; first of all we need to create an equidistant time trace by interpolation
  n_tsteps = N_ELEMENTS(time)          ; # time steps in original arrays
  n_interp = FIX(interp_fac*n_tsteps)  ; # time steps in interpol'd arrays

  newdt = (time[n_tsteps-1] - time[0]) / (n_interp - 1)
  newtime = time[0] + INDGEN(n_interp) * newdt

  data_ip = INTERPOL(arr,time,newtime) ; interpolated data

  ; now we loop over the time windows
  t_ind_start = 0
  t_ind_end = FIX(n_interp*width) - 1
  n_tind = t_ind_end - t_ind_start + 1

  T_window = newtime[t_ind_end]-newtime[t_ind_start]
  omega = 2.*!PI/T_window * (-n_tind/2+1+INDGEN(n_tind))

  fftdata_win = DCOMPLEXARR(n_tind)

  n_wins = 0
  WHILE t_ind_end LE n_interp - 1 DO BEGIN
    data_filtered = CALL_FUNCTION(filter_func,n_tind) * $
      data_ip[t_ind_start:t_ind_end]
    fftdata_win += ABS(FFT(data_filtered))
    t_ind_start += FIX(n_tind*(1.0-overlap))
    t_ind_end += FIX(n_tind*(1.0-overlap))
    n_wins += 1
  ENDWHILE
  fftdata_win = fftdata_win / n_wins

  IF myshift THEN fftdata_win = SHIFT(fftdata_win,n_tind/2-1) $
    ELSE omega = SHIFT(omega,n_tind/2-1)

  RETURN, fftdata_win

END

;##########################################################################

FUNCTION str_h5proof, str_in
; replace all instances of / and . with , in order to be able to use a
; string with HDF5

  IF N_ELEMENTS(str_in) LT 1 THEN RETURN, ''

  str_out = str_in
  FOR j = 0, N_ELEMENTS(str_in) - 1 DO BEGIN
    str_in_sep = STRSPLIT(str_in[j],"/.",/EXTRACT)

    str_out_temp = ''
    FOR k = 0, N_ELEMENTS(str_in_sep) - 1 DO $
      str_out_temp += str_in_sep[k] + ','
    str_out_temp = STRMID(str_out_temp,0,STRLEN(str_out_temp)-1)

    str_out[j] = str_out_temp
  ENDFOR

  RETURN, str_out[*]

END

;##########################################################################

FUNCTION mem_usage, highw=highw, threshold=threshold
; returns the current memory usage in MBytes; set /highw for highwater mark
; if threshold (in MB) is set, it is reset to 1 if exceeded, 0 otherwise

  IF NOT KEYWORD_SET(highw) THEN highw = 0

  mem = MEMORY(CURRENT=(1-highw),HIGHWATER=highw) / 1048576.0
  IF KEYWORD_SET(threshold) THEN threshold = threshold LT mem ? 1 : 0

  RETURN, rm0es(mem) + ' MB'

END

;##########################################################################

PRO set_widget_colors

  COMMON global_vars

  xdef_path = '$HOME/.Xdefaults'
  bg_default = 'gray'
  fg_default = 'black'

  xdef_exists = FILE_TEST(xdef_path,/READ,/WRITE)

  ; if defaults are set and no Xdefaults exists, do nothing
  IF (cfg.bg_color EQ bg_default) AND (cfg.fg_color EQ fg_default) AND $
    NOT xdef_exists THEN RETURN

  IF xdef_exists THEN BEGIN
  ; if Xdefaults exists, check whether the current values are already set
    n_lines = FILE_LINES(xdef_path)
    lines = STRARR(n_lines)

    OPENR, xdef_lun, xdef_path, /GET_LUN
    READF, xdef_lun, lines
    FREE_LUN, xdef_lun

    pos_bg = WHERE(STRMID(lines,0,TRANSPOSE(INTARR(n_lines)+14)) EQ $
      'Idl*background')
    IF pos_bg[0] NE -1 THEN BEGIN
      pos_bg_mod = N_ELEMENTS(pos_bg) NE 1 ? $
        pos_bg[N_ELEMENTS(pos_bg)-1] : pos_bg
      bg_col_xdef = rm0es(STRSPLIT(lines[pos_bg_mod],':',/EXTRACT))
      bg_col_xdef = N_ELEMENTS(bg_col_xdef) EQ 2 ? bg_col_xdef[1] : ''
    ENDIF ELSE bg_col_xdef = cfg.bg_color EQ bg_default ? bg_default : ''

    pos_fg = WHERE(STRMID(lines,0,TRANSPOSE(INTARR(n_lines)+14)) EQ $
      'Idl*foreground')
    IF pos_fg[0] NE -1 THEN BEGIN
      pos_fg_mod = N_ELEMENTS(pos_fg) NE 1 ? $
        pos_fg[N_ELEMENTS(pos_fg)-1] : pos_fg
      fg_col_xdef = rm0es(STRSPLIT(lines[pos_fg_mod],':',/EXTRACT))
      fg_col_xdef = N_ELEMENTS(fg_col_xdef) EQ 2 ? fg_col_xdef[1] : ''
    ENDIF ELSE fg_col_xdef = cfg.fg_color EQ fg_default ? fg_default : ''

    ; if the Xdefault settings are different, set the new ones
    IF bg_col_xdef NE cfg.bg_color THEN BEGIN
      IF pos_bg[0] EQ -1 THEN $
        lines = [lines,'Idl*background: '+cfg.bg_color] ELSE $
        lines[pos_bg[N_ELEMENTS(pos_bg)-1]] = $
          'Idl*background: ' + cfg.bg_color
      gui.misc.wcolors_set = 1
    ENDIF
    IF fg_col_xdef NE cfg.fg_color THEN BEGIN
      IF pos_fg[0] EQ -1 THEN $
        lines = [lines,'Idl*foreground: '+cfg.fg_color] ELSE $
        lines[pos_fg[N_ELEMENTS(pos_fg)-1]] = $
          'Idl*foreground: ' + cfg.fg_color
      gui.misc.wcolors_set = 1
    ENDIF

    ; remove duplicate entries in Xdefaults
    ind_del_nonuniq_bg = N_ELEMENTS(pos_bg) GT 1 ? $
      WHERE(INDGEN(n_lines) NE pos_bg[0:N_ELEMENTS(pos_bg)-2]) : $
      INDGEN(n_lines)
    ind_del_nonuniq = N_ELEMENTS(pos_fg) GT 1 ? $
      WHERE(ind_del_nonuniq_bg NE pos_fg[0:N_ELEMENTS(pos_fg)-2]) : $
      ind_del_nonuniq_bg
    lines = lines[ind_del_nonuniq]
  ENDIF ELSE BEGIN
  ; if Xdefaults does not exist, create it and write settings
    lines = ['Idl*background: '+cfg.bg_color,$
      'Idl*foreground: '+cfg.fg_color]
    gui.misc.wcolors_set = 1
  ENDELSE

  ; write new Xdefaults content
  IF gui.misc.wcolors_set THEN BEGIN
    OPENW, xdef_lun, xdef_path, /GET_LUN
    FOR j = 0, N_ELEMENTS(lines) - 1 DO PRINTF, xdef_lun, lines[j]
    FREE_LUN, xdef_lun
  ENDIF

END

;############################################################################

PRO set_size_thick, stored_vals, reset=reset

  COMMON global_vars

  IF NOT KEYWORD_SET(reset) THEN BEGIN
    stored_vals = [!P.CHARSIZE,!X.CHARSIZE,!Y.CHARSIZE,!Z.CHARSIZE,$
      !P.THICK,!X.THICK,!Y.THICK,!Z.THICK,!P.CHARTHICK]

    !P.CHARSIZE = cfg.charsize
    !X.CHARSIZE = cfg.charsize
    !Y.CHARSIZE = cfg.charsize
    !Z.CHARSIZE = cfg.charsize
    !P.THICK = cfg.linethick
    !X.THICK = cfg.linethick
    !Y.THICK = cfg.linethick
    !Z.THICK = cfg.linethick
    !P.CHARTHICK = cfg.charthick
    RETURN
  ENDIF

  IF N_ELEMENTS(stored_vals) NE 9 THEN BEGIN
    printerror, 'failed to reset font size and line thickness'
    RETURN
  ENDIF

  !P.CHARSIZE = stored_vals[0]
  !X.CHARSIZE = stored_vals[1]
  !Y.CHARSIZE = stored_vals[2]
  !Z.CHARSIZE = stored_vals[3]
  !P.THICK = stored_vals[4]
  !X.THICK = stored_vals[5]
  !Y.THICK = stored_vals[6]
  !Z.THICK = stored_vals[7]
  !P.CHARTHICK = stored_vals[8]

END

;##########################################################################

PRO read_cfg

  COMMON global_vars

  OPENR, cfg_lun, 'diag.cfg', /GET_LUN, ERROR=err
  IF !QUIET NE 1 THEN IF err THEN PRINT, 'no diag.cfg found' ELSE $
    PRINT, 'reading diag.cfg'
  IF err THEN RETURN

  n_lines = FILE_LINES('diag.cfg')
  lines = STRARR(n_lines)
  READF, cfg_lun, lines

  FOR j = 0, n_lines - 1 DO BEGIN
    key_val = rm0es(STRSPLIT(lines[j],'=',/EXTRACT))
    IF N_ELEMENTS(key_val) EQ 2 THEN BEGIN
      key = key_val[0]
      val = key_val[1]
      CASE key OF
        'ps_viewer'      : cfg.ps_viewer = val
        'vm_lowcase'     : cfg.vm_lowcase = val EQ 'yes'
        'short_gui'      : cfg.short_gui = val EQ 'yes'
        'speedup_gui'    : cfg.speedup_gui = val EQ 'yes'
        'always_refresh' : cfg.always_refresh = val EQ 'yes'
        'scanlog_repair' : cfg.scanlog_repair = val EQ 'yes'
        'info_str'       : cfg.info_str = val
        'charsize'       : cfg.charsize = val
        'charthick'      : cfg.charthick = val
        'linethick'      : cfg.linethick = val
        'bg_color'       : cfg.bg_color = val
        'fg_color'       : cfg.fg_color = val
        'kpar_filter'    : cfg.kpar_filter = val
        ELSE :
      ENDCASE
    ENDIF
  ENDFOR

  FREE_LUN, cfg_lun

  IF cfg.ps_viewer NE '' THEN BEGIN
    viewlen = STRLEN(cfg.ps_viewer)
    IF STRMID(cfg.ps_viewer,viewlen-1) NE ' ' THEN cfg.ps_viewer += ' '
  ENDIF

  set_widget_colors

END

;############################################################################

FUNCTION lcm, number1, number2
; returns lowest common multiple

  sum1 = number1
  sum2 = number2

  WHILE sum1 NE sum2 DO $
    IF sum1 LT sum2 THEN sum1 += number1 ELSE sum2 += number2

  RETURN, sum1

END

;############################################################################

PRO get_torus_pars, geom_type, torR, torZ, r0, round_q0, bigM, rho
; returns some parameters required for a full or partial torus visualization
; input: geom_type (0: s-alpha, 1: circular, 2: tracer/chease etc.)

  COMMON global_vars

  ; torR defines the torus major radius which should be used 
  ; to normalize the data to a torus with major radius one
  ; or to simplify (r, theta operations)
  ; As major_R might have been accidentally set to one, the
  ; average geom.R is considered first, then the par.major_R
  IF geom_type EQ 2 THEN BEGIN ;tracer/chease geometry
    torR = par.x_local ? (TOTAL((*series.geom).R*$
      (*series.geom).jacobian) / TOTAL((*series.geom).jacobian)) : $
      TOTAL((*series.geom).R[0,*]) / par.nz0
  ENDIF ELSE IF par.major_R GT 0.0 THEN torR = par.major_R * par.Lref

  ; torZ defines the average height (default: 0.0)
  IF geom_type EQ 2 THEN torZ = $
    TOTAL((*series.geom).Z*(*series.geom).jacobian) / $
    TOTAL((*series.geom).jacobian) ELSE torZ = 0.0

  ; r0 defines the position of the flux surface
  ; in the chosen flux surface label (normalized to Lref)
  ; note: might need to be set by hand
  ; as old GENE sims did not store the rho_tor value when
  ; tracer/chease had been run
  IF geom_type EQ 2 THEN BEGIN
    IF par.x0 LE 0 THEN BEGIN
      ; estimate flux surface label by average radius of 
      ; flux surface (might significantly differ from, e.g., rho_tor)
      ; r0 is normalized to Lref
      r0 = TOTAL(SQRT(((*series.geom).R-torR)^2+((*series.geom).Z-torZ)^2)*$
        (*series.geom).jacobian) / (TOTAL((*series.geom).jacobian) * par.Lref)
     ENDIF ELSE r0 = par.x0 ;assuming now that x0 is normalized to Lref in generalized geometries
                            ;before: r0 = minor_r GT 0 ? par.x0 * par.minor_r : par.x0
  ENDIF ELSE BEGIN ; s-alpha or circular
    IF (par.minor_R EQ 0.0) OR (par.x0 LT 0) THEN $
      r0 = par.trpeps * par.major_R ELSE $
      r0 = par.x0 * par.minor_r
  ENDELSE

  IF (geom_type EQ 2) AND (par.x0 LT 0) THEN BEGIN
    ; provide simulation box center by hand if necessary:
    ; r0 = 0.8
    ; PRINT, 'replaced r0 by ', rm0es(r0)

    PRINT, 'Estimated r0 (in units of Lref): ', rm0es(r0)
    PRINT, 'Compare with the actual flux surface label!'
  ENDIF ELSE PRINT, 'Using r0 (in units of Lref): ',rm0es(r0)
    
  rho_o_Lref = (SQRT(par.mref*par.Tref/par.Qref) / par.Bref/ par.Lref)

  bigM = par.n0_global

  IF bigM LE 0 THEN BEGIN
    round_q0 = 0.1 * ROUND(10.0*par.q0)

    if (par.magn_geometry EQ 'miller') then $
       bigM = FLOOR(2.0*!PI/(rho_o_Lref*par.ly)*(*series.geom).C_y) $
    else bigM = FLOOR(2.0*!PI*r0/(rho_o_Lref*par.ly*par.q0)) 
    PRINT, 'floor(bigM) using reference vals: ', rm0es(bigM)
    bigM_min = bigM

    IF ((bigM * round_q0) MOD 1) GT 1e-3 THEN BEGIN
      WHILE (((bigM_min * round_q0) MOD 1) GT 1e-3) $
        AND (bigM_min LT 100) DO bigM_min += 1
      bigM = bigM_min
    ENDIF
  ENDIF ELSE round_q0 = par.q0

  IF par.x_local THEN BEGIN
     if (par.magn_geometry EQ 'miller') then $
        rho=par.kymin*(*series.geom).C_y*par.q0/(bigM*round_q0) $
     else rho = par.kymin * r0 / (bigM * round_q0) 
    PRINT, 'rho/Lref from reference values: ', rho_o_Lref
    PRINT, 'value chosen for visualization: ', rho
    PRINT, 'number of fluxtube repetitions to cover the torus: ', bigM
  ENDIF ELSE rho = par.rhostar * par.minor_r
END

;############################################################################

FUNCTION get_y_shift, x_res, dx, dy, Cyq0_r0

  COMMON global_vars

  IF par.x_local THEN BEGIN
     y_shift_base = 2.0 * !PI * par.shat * dx / dy * Cyq0_r0
     y_shift = y_shift_base*(INDGEN(x_res)-x_res/2)
  ENDIF ELSE BEGIN
     y_shift_base = 2.0 * !PI / dy
     y_shift = y_shift_base*INTERPOL((*series.geom).C_y*(*series.geom).q,x_res)
  ENDELSE

  RETURN, y_shift

END

;############################################################################

PRO interp_3d_data, data, x_res, y_res, z_res, y_shift, infft=infft
; Interpolates data onto a new grid with resolution (x_res,y_res,z_res).
; The output is always in position space, input can be fourier-transformed,
; set according entry in 3d-logical-arr to true
;
; For z interpolation, the y_shift from the boundary condition
; can be computed on the new x grid via get_y_shift; also, the boundary
; point in z is left attached, resulting in z_res + 1 points.

  IF N_ELEMENTS(infft) LT 1 THEN infft = [0,0,0]

  IF (TOTAL(infft) EQ 0) THEN BEGIN
     nx0 = N_ELEMENTS(data[*,0,0])
     ny0 = N_ELEMENTS(data[0,*,0])
     nz0 = N_ELEMENTS(data[0,0,*])

     IF x_res NE nx0 THEN BEGIN
        data_extd = FLTARR(nx0+1,ny0,nz0,/NOZERO)
        data_extd[0,0,0] = data ; note: same as [*,*,*] but faster
        data_extd[nx0,0,0] = data[0,*,*]
        
        data = (INTERPOLATE(TEMPORARY(data_extd),$
           INDGEN(x_res+1)*(nx0/FLOAT(x_res)),$
           INDGEN(ny0),INDGEN(nz0),/GRID))[0:x_res-1,*,*]
     ENDIF

     IF y_res NE ny0 THEN BEGIN
        data_extd = FLTARR(x_res,ny0+1,nz0,/NOZERO)
        data_extd[0,0,0] = data
        data_extd[0,ny0,0] = data[*,0,*]

        data = (INTERPOLATE(TEMPORARY(data_extd),INDGEN(x_res),$
          INDGEN(y_res+1)*(ny0/FLOAT(y_res)),$
          INDGEN(nz0),/GRID))[*,0:y_res-1,*]
     ENDIF

     IF z_res NE nz0 THEN BEGIN
        y_mod = FLTARR(x_res,y_res,/NOZERO)
        FOR y = 0, y_res - 1 DO y_mod[*,y] = $
           (y - ROUND(y_shift[*])) MOD y_res
        y_mod_neg = WHERE(y_mod LT 0)
        IF y_mod_neg[0] NE -1 THEN y_mod[y_mod_neg] += y_res
    
        data_extd = FLTARR(x_res,y_res,nz0+1,/NOZERO)
        data_extd[0,0,0] = data
        data_extd[0,0,nz0] = data[REBIN(INDGEN(x_res),[x_res,y_res]),$
          y_mod[*,*],REBIN([0],[x_res,y_res])]

        data = INTERPOLATE(TEMPORARY(data_extd),$
                           INDGEN(z_res+1)*(nz0/FLOAT(z_res)))
     ENDIF

  ENDIF ELSE IF ((infft[0] EQ 0) AND (infft[1] NE 0) AND (infft[2] EQ 0)) THEN BEGIN
     nx0 = N_ELEMENTS(data[*,0,0])
     nky0 = N_ELEMENTS(data[0,*,0])
     ny0 = 2*nky0
     nz0 = N_ELEMENTS(data[0,0,*])

     ; interpolation of the z coordinate
     IF z_res NE nz0 THEN BEGIN
        data_extd = DCOMPLEXARR(nx0,nky0,nz0+1,/NOZERO)
        data_extd[0,0,0] = data ; note: same as [*,*,*] but faster

        ; apply parallel b.c.
        FOR j = 0, nky0 -1 DO BEGIN
           data_extd[*,j,nz0] = data[*,j,0]*$
                 EXP(-2.0*!PI*complex(0.0,1.0)*y_shift[*]*j)
        ENDFOR

        ; calculate sxsy data
        temp_data = DCOMPLEXARR(nx0,ny0,nz0+1,/NOZERO)
        FOR ix = 0, nx0 - 1 DO BEGIN
           FOR j = 0, nky0 - 1 DO temp_data[ix,j,*] = data_extd[ix,j,*]
           temp_data[ix,nky0,*] = COMPLEX(0.0,0.0)
           FOR j = nky0 + 1, 2 * nky0 - 1 DO temp_data[ix,j,*] = $
              CONJ(data_extd[ix,2*nky0-j,*])
        ENDFOR
        data_extd = FFT(temp_data,DIMENSION=2,/INVERSE)  

        data = INTERPOLATE(TEMPORARY(data_extd),$
          INDGEN(z_res)*(nz0/FLOAT(z_res)))

     ENDIF ELSE BEGIN
        ; calculate sxsy data
        temp_data = DCOMPLEXARR(nx0,2*nky0,z_res,/NOZERO)
        FOR ix = 0, nx0 - 1 DO BEGIN
           FOR j = 0, nky0 - 1 DO temp_data[ix,j,*] = data[ix,j,*]
           temp_data[ix,par.nky0,*] = COMPLEX(0.0,0.0)
           FOR j = nky0 + 1, 2 * nky0 - 1 DO temp_data[ix,j,*] = $
              CONJ(data[ix,2*nky0-j,*])
        ENDFOR
        data = FFT(temp_data,DIMENSION=2,/INVERSE)  
     ENDELSE

     ; interpolation of the x coordinate
     IF x_res NE nx0 THEN BEGIN
        data_extd = FLTARR(nx0+1,2*nky0,z_res,/NOZERO)
        data_extd[0,0,0] = data ; note: same as [*,*,*] but faster
        data_extd[nx0,*,*] = data[0,*,*]

        data = (INTERPOLATE(TEMPORARY(data_extd),$
           INDGEN(x_res+1)*(nx0/FLOAT(x_res)),$
           INDGEN(2*nky0),INDGEN(z_res),/GRID))[0:x_res-1,*,*]
     ENDIF
  
     ; interpolation of the y coordinate
     IF y_res NE ny0 THEN BEGIN
        data_extd = FLTARR((*i).x_res,ny0+1,(*i).z_res,/NOZERO)
        data_extd[0,0,0] = data
        data_extd[*,ny0,*] = data[*,0,*]
        
        data = (INTERPOLATE(TEMPORARY(data_extd),INDGEN(x_res),$
          INDGEN(y_res+1)*(ny0/FLOAT(y_res)),$
          INDGEN(z_res),/GRID))[*,0:y_res-1,*]
     ENDIF
  ENDIF ELSE BEGIN
     printerror, 'infft = ', infft, ' not implemented yet'
     stop
  ENDELSE
END

;############################################################################

FUNCTION get_RZ_coords, z_res, x_res, x_axis, torR=torR, torZ=torZ
; returns R[z,x] and Z[z,x] coordinates
; for this purpose, the local code requires the x axis on the same grid
; cylindrical coordinates can be optionally be returned as well
; where torR and torZ can be used to define the central point

  COMMON global_vars

  IF NOT KEYWORD_SET(torR) THEN torR = 0.0
  IF NOT KEYWORD_SET(torZ) THEN torZ = 0.0

  R_pos = FLTARR(z_res,x_res,/NOZERO)
  Z_pos = FLTARR(z_res,x_res,/NOZERO)  
  space_r = FLTARR(z_res,x_res,/NOZERO)
  space_theta = FLTARR(z_res,x_res,/NOZERO)  

  IF par.x_local THEN BEGIN
    Z_extd = [(*series.geom).Z,(*series.geom).Z[0]]
    Z_interp = (INTERPOL(Z_extd,z_res+1,/SPLINE))[0:z_res-1]

    R_extd = [(*series.geom).R,(*series.geom).R[0]]
    R_interp = (INTERPOL(R_extd,z_res+1,/SPLINE))[0:z_res-1]

    c1_extd = [(*series.geom).c1/(*series.geom).gxx,$
      (*series.geom).c1[0]/(*series.geom).gxx[0]]

    c1_interp = (INTERPOL(c1_extd,z_res+1,/SPLINE))[0:z_res-1]

    c2_extd = [(*series.geom).c2/(*series.geom).gxx,$
      (*series.geom).c2[0]/(*series.geom).gxx[0]]

    c2_interp = (INTERPOL(c2_extd,z_res+1,/SPLINE))[0:z_res-1]
     
    FOR x = 0, x_res - 1 DO BEGIN
      R_pos[*,x] = R_interp+c1_interp*x_axis[x]*par.Lref
      Z_pos[*,x] = Z_interp+c2_interp*x_axis[x]*par.Lref

      theta_interp = ATAN(Z_interp-torZ,R_interp-torR)
      space_theta[*,x] = ATAN(Z_pos[*,x]-torZ,R_pos[*,x]-torR)

      ; in a local simulation using the realistic rhostar,
      ; the radial box length might intersect the origin
      ; and might thus show weird results;
      ; hence, check whether radius vector intersects the origin
      ; and set r to zero if so
      FOR z = 0, z_res - 1 DO BEGIN
        IF ((ABS(theta_interp[z]-space_theta[z,x]) $
          GT !PI/4.0) AND (ABS(theta_interp[z]-space_theta[z,x]) $
          LT 1.5*!PI)) THEN space_r[z,x] = 0.0 ELSE $
          space_r[z,x] = SQRT((R_pos[z,x]-torR)^2+(Z_pos[z,x]-torZ)^2)
      ENDFOR
    ENDFOR
  ENDIF ELSE BEGIN ; nonlocal realistic geometry
    IF x_res NE par.nx0 THEN BEGIN
      printerror, 'only x resolution fixed to nx0 is implemented'
      STOP
    ENDIF
    FOR x = 0, x_res - 1 DO BEGIN
      Z_extd = [REFORM((*series.geom).Z[x,*]),$
        REFORM((*series.geom).Z[x,0])]
      Z_pos[*,x] = (INTERPOL(Z_extd,z_res+1,/SPLINE))[0:z_res-1]
      R_extd = [REFORM((*series.geom).R[x,*]),$
        REFORM((*series.geom).R[x,0])]
      R_pos[*,x] = (INTERPOL(R_extd,z_res+1,/SPLINE))[0:z_res-1]
    ENDFOR
    space_r = SQRT((R_pos-torR)^2+(Z_pos-torZ)^2)
    space_theta = ATAN((Z_pos-torZ),(R_pos-torR))
  ENDELSE

  RETURN, {R:R_pos,Z:Z_pos,space_r:space_r,space_theta:space_theta}

END



;############################################################################
; Lagrange interpolation routines as implemented in GENE
; > Third order lagrange interpolation
PRO lag3interp, y_in,x_in,n_in,y_out,x_out,n_out
  IF (x_in[n_in-1] GT x_in[0]) THEN BEGIN
     jstart=3
     jfirst=1
     jlast=n_out
     jstep=1
  ENDIF ELSE BEGIN
     jstart=n_in-2
     jfirst=n_out
     jlast=1
     jstep=-1
  ENDELSE

  j1=jstart
  FOR j=jfirst,jlast,jstep DO BEGIN
     x=x_out[j-1]
     WHILE ((x GE x_in[j1-1]) AND (j1 LT n_in-1) AND (j1 GT 2)) DO $
        j1=j1+jstep

     j2=j1+jstep
     j0=j1-jstep
     jm=j1-2*jstep

     ;...  extrapolate inside or outside
     x2=x_in[j2-1]
     x1=x_in[j1-1]
     x0=x_in[j0-1]
     xm=x_in[jm-1]

     aintm=(x-x0)*(x-x1)*(x-x2)/((xm-x0)*(xm-x1)*(xm-x2))
     aint0=(x-xm)*(x-x1)*(x-x2)/((x0-xm)*(x0-x1)*(x0-x2))
     aint1=(x-xm)*(x-x0)*(x-x2)/((x1-xm)*(x1-x0)*(x1-x2))
     aint2=(x-xm)*(x-x0)*(x-x1)/((x2-xm)*(x2-x0)*(x2-x1))
     
     y_out[j-1]=aintm*y_in[jm-1]+aint0*y_in[j0-1] $
                +aint1*y_in[j1-1]+aint2*y_in[j2-1]
  
  ENDFOR

END

;############################################################################
;>Returns Derivative based on a 3rd order lagrange interpolation
PRO lag3deriv,y_in,x_in,n_in,dydx_out,x_out,n_out

  IF (x_in[n_in-1] GT x_in[0]) THEN BEGIN
     jstart=3
     jfirst=1
     jlast=n_out
     jstep=1
  ENDIF ELSE BEGIN
     jstart=n_in-2
     jfirst=n_out
     jlast=1
     jstep=-1
  ENDELSE

  j1=jstart
  FOR j=jfirst,jlast,jstep DO BEGIN
     x=x_out[j-1]
     WHILE ((x GE x_in[j1-1]) AND (j1 LT n_in-1) AND (j1 GT 2)) DO $
        j1=j1+jstep

     j2=j1+jstep
     j0=j1-jstep
     jm=j1-2*jstep

     ;...  extrapolate inside or outside
     x2=x_in[j2-1]
     x1=x_in[j1-1]
     x0=x_in[j0-1]
     xm=x_in[jm-1]

     aintm=((x-x1)*(x-x2)+(x-x0)*(x-x2)+(x-x0)*(x-x1)) $
           /((xm-x0)*(xm-x1)*(xm-x2))
     aint0=((x-x1)*(x-x2)+(x-xm)*(x-x2)+(x-xm)*(x-x1)) $
           /((x0-xm)*(x0-x1)*(x0-x2))
     aint1=((x-x0)*(x-x2)+(x-xm)*(x-x2)+(x-xm)*(x-x0)) $
           /((x1-xm)*(x1-x0)*(x1-x2))
     aint2=((x-x0)*(x-x1)+(x-xm)*(x-x1)+(x-xm)*(x-x0)) $
           /((x2-xm)*(x2-x0)*(x2-x1))

     dydx_out[j-1]=aintm*y_in[jm-1]+aint0*y_in[j0-1] $
                 +aint1*y_in[j1-1]+aint2*y_in[j2-1]

  ENDFOR
  
END

;############################################################################

PRO stat_avg_err, time, data, avg, err, tcorr=tcorr, wwidth=wwidth
; time: array with simulation times (input)
; data: array with data at the simulation times (input)
; avg: rms of data values (return value)
; err: standard deviation of avg (return value)
; tcorr: correlation time, basis for time windows for error computation;
;   0: sets err = 0 (default); >0.0: manual value; -2: standard deviation
; wwidth: will contain window width

  n_steps = N_ELEMENTS(time)
  IF n_steps NE N_ELEMENTS(data) THEN BEGIN
    printerror, 'stat_avg_err error: incompatibe time, data arrays', /popup
    STOP
  ENDIF

  IF (N_ELEMENTS(tcorr) NE 1) OR (n_steps LE 2) THEN BEGIN
    tcorr = 0
    PRINT, 'stat_avg_err warning: skipping error computation'
  ENDIF

  err = 0

  ; --- compute average ---

  dt_full = time[n_steps-1] - time[0]
  avg = INT_TABULATED(time,data,/DOUBLE) / dt_full
  IF tcorr EQ 0 THEN RETURN

  ; --- use window statistics for errors ---

  IF tcorr GT 0 THEN BEGIN
    ; start with window width = 5 tcorr; if fewer than 10, reduce width, min 2
    IF dt_full GE 50 * tcorr THEN BEGIN
      n_win = FLOOR(dt_full/(5*tcorr))
      wwidth = dt_full / n_win
    ENDIF ELSE BEGIN
      IF dt_full GE 20 * tcorr THEN BEGIN
        n_win = 10
        wwidth = dt_full / n_win
      ENDIF ELSE BEGIN
        IF !QUIET NE 1 THEN PRINT, $
          'stat_avg_err warning: insufficient time range, skipping error computation'
        tcorr = 0
        RETURN
      ENDELSE
    ENDELSE

    win_inds = LONARR(n_win+1)
    FOR n = 0, n_win - 1 DO win_inds[n] = (WHERE(time-time[0] GE n*wwidth))[0]
    win_inds[n_win] = n_steps - 1
    IF (WHERE(win_inds LT 0))[0] GE 0 THEN BEGIN
      printerror, 'stat_avg_err error: could not detect windows', /popup
      STOP
    ENDIF

    ; ensure window has 2+ data points
    IF (WHERE((win_inds[1:n_win]-win_inds[0:n_win-1]) LE 2))[0] GE 0 THEN BEGIN
      PRINT, 'stat_avg_err warning: windows with <= 2 steps, skipping error computation'
      tcorr = 0
      RETURN
    ENDIF

    win_avg = FLTARR(n_win)
    FOR n = 0, n_win - 1 DO BEGIN
      win_ind_0 = win_inds[n]
      win_ind_1 = win_inds[n+1] - 1
      win_avg[n] = INT_TABULATED(time[win_ind_0:win_ind_1],$
        data[win_ind_0:win_ind_1],/DOUBLE) / (time[win_ind_1]-time[win_ind_0])
    ENDFOR

    err = STDDEV(win_avg,/DOUBLE)
  ENDIF

  ; --- return simple standard deviation ---

  IF tcorr EQ -2 THEN err = STDDEV(data,/DOUBLE)

END
