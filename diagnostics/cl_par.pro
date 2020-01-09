PRO read_cl_params

  COMMON global_vars

  OPENR, parlun, 'cl_params', /GET_LUN

  line = ''

  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  gui.out.data_path = val
  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  gui.out.out_path = val

  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  gui.out.start_t = DOUBLE(val)
  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  gui.out.end_t = DOUBLE(val)

  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  gui.out.sparsefac = LONG(val)

  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  runs_input = val
  run_labels = convert_run_string(runs_input)
  series.run_labels = run_labels
  run_label_arr = STRSPLIT(run_labels,',',/EXTRACT)
  series.n_runs = N_ELEMENTS(run_label_arr)

  file_par = get_file_path('parameters',/set_fm)

  IF PTR_VALID(nrg) THEN PTR_FREE, nrg
  IF PTR_VALID(series.geom) THEN PTR_FREE, series.geom
  read_par
  set_series_lengths

  READF, parlun, line
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  spec_sel = FIX(val)
  gui.out.spec_select = PTR_NEW(spec_sel)

  FREE_LUN, parlun

END

;#########################################################################

PRO cl_setdiags, momtab=momtab

  COMMON global_vars

  OPENR, parlun, 'cl_params', /GET_LUN

  line = ''

  FOR j = 0, 6 DO READF, parlun, line
  READF, parlun, line
  FREE_LUN, parlun
  val = rm0es((STRSPLIT(line,';',/EXTRACT))[0])
  diagline = val

  table_array = *gui.tab.arr
  momloc = WHERE(table_array.name EQ 'mom')
  IF momloc[0] EQ -1 THEN BEGIN
    momtab = -1
    RETURN
  ENDIF
  momtab = momloc[0]
  table_array = table_array[momloc]

  diag = get_diag_ptr(table_array.diag0,0)

  found = 0
  WHILE (found EQ 0) AND PTR_VALID(diag) DO BEGIN
    IF (*diag).name EQ diagline[0] THEN BEGIN
      found = 1
      (*diag).selected = 1
    ENDIF

    diag = (*diag).next
  ENDWHILE

END
