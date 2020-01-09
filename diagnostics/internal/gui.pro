PRO gui_main

  CLOSE, /ALL
  COMMON global_vars

  read_cfg
  create_diag_struct

  set_output, /reset, dname='dummy' ; for drawing white plot window, see out_device.pro
  load_gui_elements                 ; those functions are in diag_struct.pro
  create_diag_table
  load_form, /last_form, /preserve_msg
  create_diag_table, /refresh
  navhist, /add_entry

  set_size_thick, stored_size_thick

  WIDGET_CONTROL, gui.window.main, /REALIZE
  load_gui_logo
  XMANAGER, 'gui_main', gui.window.main

  free_diag_memory

  PTR_FREE, nrg, vsp, scanlog, $
    gui.out.spec_select, gui.tab.arr,$
    par.lxarr, par.ky, par.kx, par.z, $
    par.istep_field, par.istep_mom, par.istep_nrg, par.istep_energy, $
    series.filename_fmt, series.geom, par.istep_nlt, par.istep_energy3d, $
    series.ky, series.lxarr, series.kx, spec[*].prof, $
    par.kx_nlt_ind, par.ky_nlt_ind

  IF !QUIET NE 1 THEN HELP, /HEAP_VARIABLES

  HEAP_GC
  IF !QUIET NE 1 THEN BEGIN
    PRINT, 'Maximum memory use since last mem_usage call: ' + $
      mem_usage(/highw)
    PRINT, 'Memory still in use: ' + mem_usage()
  ENDIF

  set_size_thick, stored_size_thick, /reset

END

;#########################################################################

PRO load_gui_elements

  COMMON global_vars

  gui.window.main = WIDGET_BASE($
    TITLE     = 'GENE diagnostics',$
    SCR_XSIZE = 1000,$
    SCR_YSIZE = cfg.short_gui ? 720 : 800,$
    XOFFSET   = 200,$
    YOFFSET   = 10,$
    EVENT_PRO = 'gui_main_event')

  plotwindow = WIDGET_DRAW($
    gui.window.main,$
    XSIZE    = 480,$
    YSIZE    = 380,$
    XOFFSET  = 10,$
    FRAME    = 1,$
    YOFFSET  = 10,$
    RETAIN   = 2)

  gui.tab.base = WIDGET_TAB($
    gui.window.main,$
    XOFFSET   = 10,$
    YOFFSET   = 410,$
    SCR_XSIZE = 820,$
    SCR_YSIZE = cfg.short_gui ? 300 : 380)

  len = MAX(STRLEN((*gui.tab.arr).name))

  FOR j = 0, N_ELEMENTS(*gui.tab.arr) - 1 DO $
    (*gui.tab.arr)[j].id = WIDGET_BASE($
    gui.tab.base,$
    TITLE  = STRING((*gui.tab.arr)[j].name,$
    FORMAT = '(A-'+rm0es(len)+')'))

  gui.text.data_path = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '',$
    XOFFSET   = 600,$
    YOFFSET   = 10,$
    XSIZE     = 36,$
    /EDITABLE)

  gui.button.path_bbnav = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = '<<',$
    XOFFSET   = 840,$
    YOFFSET   = 10,$
    SCR_XSIZE = 20,$
    SCR_YSIZE = 30)

  gui.button.path_backnav = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = '<',$
    XOFFSET   = 860,$
    YOFFSET   = 10,$
    SCR_XSIZE = 20,$
    SCR_YSIZE = 30)

  gui.button.path_forwnav = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = '>',$
    XOFFSET   = 880,$
    YOFFSET   = 10,$
    SCR_XSIZE = 20,$
    SCR_YSIZE = 30)

  gui.button.path_ffnav = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = '>>',$
    XOFFSET   = 900,$
    YOFFSET   = 10,$
    SCR_XSIZE = 20,$
    SCR_YSIZE = 30)

  gui.button.tmp1_path = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'tmp1',$
    XOFFSET   = 920,$
    YOFFSET   = 10,$
    SCR_XSIZE = 35,$
    SCR_YSIZE = 30)

  gui.button.tmp2_path = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'tmp2',$
    XOFFSET   = 955,$
    YOFFSET   = 10,$
    SCR_XSIZE = 35,$
    SCR_YSIZE = 30)

  gui.button.data_path = WIDGET_BUTTON($
    gui.window.main,$
    VALUE     = 'data path:',$
    XOFFSET   = 510,$
    YOFFSET   = 10,$
    SCR_XSIZE = 85,$
    SCR_YSIZE = 30)

  gui.button.out_path = WIDGET_BUTTON($
    gui.window.main,$
    VALUE     = 'output path:',$
    XOFFSET   = 510,$
    YOFFSET   = 40,$
    SCR_XSIZE = 85,$
    SCR_YSIZE = 30)

  gui.text.out_path = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '',$
    XOFFSET   = 600,$
    YOFFSET   = 40,$
    XSIZE     = 36,$
    /EDITABLE)

  widgetlabel = WIDGET_LABEL($
    gui.window.main,$
    VALUE     = 'default: output',$
    XOFFSET   = 840,$
    YOFFSET   = 46)

  gui.button.runs = WIDGET_BUTTON($
    gui.window.main,$
    VALUE     = 'runs:',$
    XOFFSET   = 510,$
    YOFFSET   = 80,$
    SCR_XSIZE = 85,$
    SCR_YSIZE = 30)

  gui.text.runs = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '',$
    XOFFSET   = 600,$
    YOFFSET   = 80,$
    XSIZE     = 36,$
    /EDITABLE)

  gui.button.series_info = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'series info',$
    XOFFSET   = 840,$
    YOFFSET   = 80,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  gui.button.nrg_data = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'nrg data',$
    XOFFSET   = 840,$
    YOFFSET   = 110,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

 gui.button.refvalues_info = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'reference values',$
    XOFFSET   = 840,$
    YOFFSET   = 140,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  gui.button.geometry = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'geometry',$
    XOFFSET   = 840,$
    YOFFSET   = 180,$
    SCR_XSIZE = 150 - show_global * 76,$
    SCR_YSIZE = 30)

  IF show_global THEN $
    gui.button.profiles = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'profiles',$
    XOFFSET   = 916,$
    YOFFSET   = 180,$
    SCR_XSIZE = 74,$
    SCR_YSIZE = 30)

  widgetlabel = WIDGET_LABEL($
    gui.window.main,$
    VALUE     = 'start time:',$
    XOFFSET   = 510,$
    YOFFSET   = 116)

  gui.text.start_t = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '',$
    XOFFSET   = 600,$
    YOFFSET   = 110,$
    XSIZE     = 9,$
    /EDITABLE)

  gui.button.first_step = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'first',$
    XOFFSET   = 670,$
    YOFFSET   = 110,$
    SCR_XSIZE = 50,$
    SCR_YSIZE = 30)

  widgetlabel = WIDGET_LABEL($
    gui.window.main,$
    VALUE     = 'end time:',$
    XOFFSET   = 510,$
    YOFFSET   = 146)

  gui.text.end_t = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '',$
    XOFFSET   = 600,$
    YOFFSET   = 140,$
    XSIZE     = 9,$
    /EDITABLE)

  gui.button.last_step = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'last',$
    XOFFSET   = 670,$
    YOFFSET   = 140,$
    SCR_XSIZE = 50,$
    SCR_YSIZE = 30)

  gui.bgroup.resolve_steps = CW_BGROUP($
    gui.window.main,$
    ['on','off'],$
    LABEL_LEFT = 'time'+STRING(10B)+'average',$
    /EXCLUSIVE,$
    /FRAME,$
    COLUMN    = 1,$
    SPACE     = -4,$
    XPAD      = -1,$
    YPAD      = -1,$
    XOFFSET   = 727,$
    YOFFSET   = 110,$
    XSIZE     = 50,$
    YSIZE     = 50)

  widgetlabel = WIDGET_LABEL($
    gui.window.main,$
    VALUE     = 'sort method:',$
    XOFFSET   = 510,$
    YOFFSET   = 186)

  gui.text.sortmeth = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '0',$
    XOFFSET   = 600,$
    YOFFSET   = 180,$
    XSIZE     = !D.X_CH_SIZE,$
    /EDITABLE)

  widgetlabel = WIDGET_LABEL($
    gui.window.main,$
    VALUE     = 'sparse factor:',$
    XOFFSET   = 680,$
    YOFFSET   = 186)

  gui.text.sparsefac = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '1',$
    XOFFSET   = 780,$
    YOFFSET   = 180,$
    XSIZE     = !D.X_CH_SIZE,$
    /EDITABLE)

  widgetlabel = WIDGET_LABEL($
    gui.window.main,$
    VALUE     = 'normalize to',$
    XOFFSET   = 675,$
    YOFFSET   = 223)

  xsize = !D.X_CH_SIZE * (5 + MAX(STRLEN(gui.out.out_format_str)))
  ysize = 10 + ROUND(2.1*!D.Y_CH_SIZE*N_ELEMENTS(gui.out.out_format))
  gui.bgroup.out_format = CW_BGROUP($
    gui.window.main,$
    gui.out.out_format_str,$
    /NONEXCLUSIVE,$
    LABEL_TOP = 'output format',$
    /FRAME,$
    SET_VALUE = gui.out.out_format,$
    XOFFSET   = 840,$
    YOFFSET   = 220,$
    XSIZE     = xsize < 150,$
    YSIZE     = ysize < 121,$
    COLUMN    = 1,$
    SPACE     = -4,$
    SCROLL    = (xsize GE 150) OR (ysize GE 121),$
    X_SCROLL_SIZE = 80 * ((xsize GE 150) OR (ysize GE 121)),$
    Y_SCROLL_SIZE = 80 * ((xsize GE 150) OR (ysize GE 121)))

  gui.text.message = WIDGET_TEXT($
    gui.window.main,$
    VALUE     = '(no message)',$
    XOFFSET   = 510,$
    YOFFSET   = 365,$
    XSIZE     = 77)

  gui.button.start = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'start',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 420 : 430,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  gui.button.load = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'load form',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 460 : 480,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  gui.button.save = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'save form',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 490 : 510,$
    SCR_XSIZE = 89,$
    SCR_YSIZE = 30)

  gui.button.save_as = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = '... as',$
    XOFFSET   = 931,$
    YOFFSET   = cfg.short_gui ? 490 : 510,$
    SCR_XSIZE = 59,$
    SCR_YSIZE = 30)

  gui.button.clear = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'clear form',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 520 : 540,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  gui.button.var_list = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'variable list',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 550 : 590,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  IF vm_sav EQ 0 THEN gui.button.cdiags = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'custom diags',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 580 : 620,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

  gui.button.recent_h5 = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'recent H5',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 610 : 670,$
    SCR_XSIZE = 74,$
    SCR_YSIZE = 30)

  gui.button.h5_files = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'H5 files',$
    XOFFSET   = 916,$
    YOFFSET   = cfg.short_gui ? 610 : 670,$
    SCR_XSIZE = 74,$
    SCR_YSIZE = 30)

  gui.button.recent_ps = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'recent ps',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 640 : 700,$
    SCR_XSIZE = 74,$
    SCR_YSIZE = 30)

  gui.button.ps_files = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'ps files',$
    XOFFSET   = 916,$
    YOFFSET   = cfg.short_gui ? 640 : 700,$
    SCR_XSIZE = 74,$
    SCR_YSIZE = 30)

  gui.button.exit = WIDGET_BUTTON($
    gui.window.main,$
    VALUE = 'exit',$
    XOFFSET   = 840,$
    YOFFSET   = cfg.short_gui ? 680 : 750,$
    SCR_XSIZE = 150,$
    SCR_YSIZE = 30)

END

;#########################################################################

PRO load_gui_logo

  COMMON global_vars

  logo_path = 'internal/gui_logo.png'

  png_exists = QUERY_PNG(logo_path)
  IF NOT png_exists THEN BEGIN
    PRINT, 'error loading logo file ' + logo_path
    RETURN
  ENDIF

  logo = READ_PNG(logo_path)
  TVSCL, logo, /TRUE

END

;#########################################################################

PRO create_Lref_droplist
; droplist with all possible normalizations

  COMMON global_vars

  dl_string = ['L_ref']

  IF par.magn_geometry EQ 's_alpha' THEN dl_string = [dl_string,'R_0'] $
    ELSE IF ((par.magn_geometry EQ 'circular') OR $
                 (STRPOS(par.magn_geometry,'miller') GE 0)) THEN $
    dl_string = [dl_string,'R_0','a']

  FOR isp = 0, par.n_spec-1 DO BEGIN
    IF spec[isp].omt NE 0.0 THEN $
      dl_string = [dl_string,'L_T_'+spec[isp].name]
    IF spec[isp].omn NE 0.0 THEN $
      dl_string = [dl_string,'L_n_'+spec[isp].name]
  ENDFOR

  IF gui.droplist.Lref NE 0 THEN BEGIN
    WIDGET_CONTROL, gui.droplist.Lref, SET_VALUE=dl_string
  ENDIF ELSE BEGIN
    gui.droplist.Lref_base = WIDGET_BASE($
      gui.window.main,$
      SCR_XSIZE = 110,$
      SCR_YSIZE = 30,$
      XOFFSET   = 675,$
      YOFFSET   = 240)

    gui.droplist.Lref = WIDGET_DROPLIST($
      gui.droplist.Lref_base,$
      VALUE     = dl_string,$
      XOFFSET   = 0,$
      YOFFSET   = 0)
  ENDELSE

  gui.out.Lref_sel = 0
  series.Lref = set_Lref(gui.out.Lref_sel)

END

;#########################################################################

PRO create_mref_droplist

  COMMON global_vars

  IF add_norm THEN BEGIN ; this droplist is only built for nonstandard use!
    dl_string = 'm_' + spec.name
    value = (WHERE(spec.mass EQ 1.0))[0]
    IF value EQ -1 THEN BEGIN
      dl_string = ['m_ref',dl_string]
      value = 0
    ENDIF

    IF gui.droplist.mref NE 0 THEN BEGIN ; if droplist exists
      WIDGET_CONTROL, gui.droplist.mref, SET_VALUE=dl_string
    ENDIF ELSE BEGIN
      gui.droplist.mref_base = WIDGET_BASE($
        gui.window.main,$
        SCR_XSIZE = 110,$
        SCR_YSIZE = 30,$
        XOFFSET   = 675,$
        YOFFSET   = 270)

      gui.droplist.mref = WIDGET_DROPLIST($ ; if droplist does not exist
        gui.droplist.mref_base,$
        VALUE     = dl_string,$
        XOFFSET   = 0,$
        YOFFSET   = 0)
    ENDELSE

    WIDGET_CONTROL, gui.droplist.mref, SET_DROPLIST_SELECT=value
  ENDIF

  series.mref = 1.0D

END

;#########################################################################

PRO create_Tref_droplist

  COMMON global_vars

  IF add_norm THEN BEGIN ; this droplist is only built for nonstandard use!
    dl_string = 'T_' + spec.name
    value = (WHERE(spec.temp EQ 1.0))[0]
    IF value EQ -1 THEN BEGIN
      dl_string = ['T_ref',dl_string]
      value = 0
    ENDIF

    IF gui.droplist.Tref NE 0 THEN WIDGET_CONTROL, gui.droplist.Tref, $
      SET_VALUE=dl_string ELSE BEGIN

      gui.droplist.Tref_base = WIDGET_BASE($
        gui.window.main,$
        SCR_XSIZE = 110,$
        SCR_YSIZE = 30,$
        XOFFSET   = 675,$
        YOFFSET   = 300)

      gui.droplist.Tref = WIDGET_DROPLIST($
        gui.droplist.Tref_base,$
        VALUE = dl_string,$
        XOFFSET   = 0,$
        YOFFSET   = 0)
    ENDELSE

    WIDGET_CONTROL, gui.droplist.Tref, SET_DROPLIST_SELECT=value
  ENDIF

  series.Tref = 1.0D

END

;#########################################################################

PRO create_species_bgroup

  COMMON global_vars

  IF cfg.speedup_gui THEN WIDGET_CONTROL, gui.window.main, UPDATE=0

  IF gui.bgroup.species NE 0 THEN $
    WIDGET_CONTROL, gui.bgroup.species, /DESTROY

  xsize = !D.X_CH_SIZE * (5 + MAX(STRLEN(spec.name))) > 60
  ysize = (10 + ROUND(2.1*!D.Y_CH_SIZE*par.n_spec)) > 25

  gui.bgroup.species = CW_BGROUP($
    gui.window.main,$
    spec.name,$
    /NONEXCLUSIVE,$
    LABEL_TOP = 'particle species',$
    /FRAME,$
    SET_VALUE = INTARR(par.n_spec > 1),$
    XOFFSET   = 510,$
    YOFFSET   = 220,$
    XSIZE     = xsize < 150,$
    YSIZE     = ysize < 121,$
    COLUMN    = 1,$
    SPACE     = -4,$
    SCROLL    = (xsize GE 150) OR (ysize GE 121),$
    X_SCROLL_SIZE = 80 * ((xsize GE 150) OR (ysize GE 121)),$
    Y_SCROLL_SIZE = 80 * ((xsize GE 150) OR (ysize GE 121)))

  IF cfg.speedup_gui THEN WIDGET_CONTROL, gui.window.main, UPDATE=1

END

;#########################################################################

PRO create_diag_table, tab=tab, refresh=refresh, rebuild=rebuild

  COMMON global_vars

  IF KEYWORD_SET(rebuild) THEN $
    WIDGET_CONTROL, (*gui.tab.arr)[tab].table_id,/DESTROY

  IF NOT KEYWORD_SET(tab) THEN BEGIN
    tstart=0
    tend = N_ELEMENTS(*gui.tab.arr)-1
  ENDIF ELSE BEGIN
    tstart = tab
    tend = tab
  ENDELSE

  FOR tab=tstart,tend DO BEGIN
   IF (*gui.tab.arr)[tab].n_diags GT 0 THEN BEGIN
    n_params = (*gui.tab.arr)[tab].n_params
    n_diags = (*gui.tab.arr)[tab].n_diags
    cell_width = INTARR(n_params+3)         ; nr of columns
    cell_width[*] = 110                              ; set cell size
    cell_width[0:2] = [18,180,18]
    tab_value = STRARR(n_params+3,n_diags)
    tab_editable = INTARR(n_params+3,n_diags)

    column_init = STRARR(n_params+3)        ; initialize column layout
    FOR i = 0, n_params - 1 DO column_init[i+3] = '--'
    i = 0
    tab_diag = (*gui.tab.arr)[tab].diag0           ; read diag content
    WHILE PTR_VALID(tab_diag) DO BEGIN
      tab_value[0,i] = '?'
      tab_value[1,i] = (*tab_diag).title
      IF (*tab_diag).selected THEN tab_value[2,i] = 'X'
      vars = (*tab_diag).table_entry
      FOR l = 0, n_params - 1 DO BEGIN
        IF (vars[0,0] NE '') AND (l LE N_ELEMENTS(vars[0,*]) - 1) THEN BEGIN
          IF vars[1,l] EQ '0' THEN BEGIN
            tab_value[3+l,i] = vars[3,l]
            tab_editable[3+l,i] = 1                ; allow cell edit
          ENDIF
          IF vars[1,l] EQ '1' THEN BEGIN
            IF STRTRIM(vars[3,l],2) EQ '1' THEN $
              tab_value[3+l,i] = '  on' ELSE tab_value[3+l,i] = '  off'
          ENDIF
        ENDIF ELSE tab_value[3+l,i] = '--'
      ENDFOR ;
      i = i + 1
      tab_diag = (*tab_diag).next
    ENDWHILE

    IF NOT KEYWORD_SET(refresh) THEN BEGIN         ; create table widget
      (*gui.tab.arr)[tab].table_id = WIDGET_TABLE($
        (*gui.tab.arr)[tab].id,$
        COLUMN_WIDTHS = cell_width,$
        /SCROLL,$
        SCR_XSIZE = 810,$
        SCR_YSIZE = cfg.short_gui ? 270 : 350,$
        X_SCROLL_SIZE = 810,$
        XSIZE     = (*gui.tab.arr)[tab].n_params+3,$
        YSIZE     = (*gui.tab.arr)[tab].n_diags,$
        VALUE     = tab_value,$
        COLUMN_LABELS = column_init,$
        /ALL_EVENTS,$
        EDITABLE=tab_editable,$
        /NO_ROW_HEADERS)
    ENDIF ELSE BEGIN
      WIDGET_CONTROL, (*gui.tab.arr)[tab].table_id, SET_VALUE=tab_value,$
        SET_TABLE_SELECT=[-1,-1,-1,-1]
    ENDELSE

   ENDIF ; --- n_diags > 0 ??
  ENDFOR ; --- loop over tabs
  tab = tend

END

;############################################################################

PRO apply_form
; checks for changes in main gui elements and
; transfers GUI information to gui.out
; calls read_par

  COMMON global_vars

  WIDGET_CONTROL, gui.text.data_path, GET_VALUE=data_path
  WIDGET_CONTROL, gui.text.out_path, GET_VALUE=out_path
  WIDGET_CONTROL, gui.text.runs, GET_VALUE=runs_input
  WIDGET_CONTROL, gui.text.start_t, GET_VALUE=start_t
  WIDGET_CONTROL, gui.text.end_t, GET_VALUE=end_t
  WIDGET_CONTROL, gui.text.sparsefac, GET_VALUE=sparsefac
  WIDGET_CONTROL, gui.text.sortmeth, GET_VALUE=sortmeth
  tab = WIDGET_INFO(gui.tab.base,/TAB_CURRENT)

  path_length = STRLEN(data_path)
  IF STRMID(data_path,path_length-1,1) NE '/' THEN data_path += '/'
  IF gui.out.data_path NE data_path THEN BEGIN
    gui.out.data_path = data_path
    series.par_updated = 0
  ENDIF

  path_length = STRLEN(out_path)
  IF STRMID(out_path,path_length-1,1) NE '/' THEN out_path += '/'
  IF path_length EQ 0 THEN out_path = 'output/'
  gui.out.out_path = out_path
  IF start_t EQ '' THEN start_t = '0.0'
  gui.out.start_t = DOUBLE(start_t)
  IF end_t EQ '' THEN end_t = '1e10'
  gui.out.end_t = DOUBLE(end_t)

  run_labels = convert_run_string(runs_input)
  run_label_arr = STRSPLIT(run_labels,',',/EXTRACT)

  IF run_labels NE series.run_labels THEN BEGIN
    series.run_labels = run_labels
    series.par_updated = 0
  ENDIF ELSE IF (STRPOS(run_labels,'act') GE 0) OR $
    (STRPOS(run_labels,'dat') GE 0) OR cfg.always_refresh OR $
    (STRPOS(data_path,'scanfiles') GE 0) THEN series.par_updated = 0

  series.n_runs = N_ELEMENTS(run_label_arr)
  file_par = get_file_path('parameters',/set_fm)

  update_scanlog = 0

  IF series.par_updated EQ 0 THEN BEGIN
    IF PTR_VALID(nrg) THEN PTR_FREE, nrg
    IF PTR_VALID(series.geom) THEN PTR_FREE, series.geom
    oldpar = par
    oldspec = spec
    read_par
    IF series.par_updated THEN BEGIN
      set_series_lengths
      IF (par.ky0_ind GT 0) AND !QUIET NE 1 THEN PRINT, $
        'ALERT: adding dummy modes ky=[0,kymin,..,(ky0_ind-1)*kymin] ' + $
        'to allow for ky->y FFT if ky0_ind > 0.'
    ENDIF

    new_species = (par.n_spec NE oldpar.n_spec)
    IF NOT new_species THEN $
      new_species = (TOTAL(spec[*].name NE oldspec[*].name) NE 0)
    IF new_species THEN BEGIN
      create_species_bgroup
      create_mref_droplist
      create_Tref_droplist
    ENDIF

    new_Lref_norm = (par.curv NE oldpar.curv) OR new_species
    IF NOT new_Lref_norm THEN $
      new_Lref_norm = (TOTAL(spec[*].omn NE oldspec[*].omn) NE 0) OR $
        (TOTAL(spec[*].omt NE oldspec[*].omt) NE 0)
    IF new_Lref_norm THEN create_Lref_droplist

    update_scanlog = 1
  ENDIF

  IF sparsefac EQ '' THEN sparsefac = '1'
  gui.out.sparsefac = LONG(sparsefac) > 1L

  IF (*gui.tab.arr)[tab].name EQ 'scan' THEN BEGIN
   IF update_scanlog OR gui.out.sortmeth NE FIX(sortmeth) THEN BEGIN
     IF PTR_VALID(scanlog) THEN PTR_FREE, scanlog
      gui.out.sortmeth = FIX(sortmeth)
      read_scanlog
   ENDIF
  ENDIF

  check_bgroups, /no_reset

END

;############################################################################

PRO load_form, last_form=last_form, preserve_msg=preserve_msg

  COMMON global_vars

  IF NOT KEYWORD_SET(preserve_msg) THEN printerror, '' ELSE $
    IF gui.misc.wcolors_set THEN $
    printerror, 'widget colors were set, restart IDL to refresh'

  create_species_bgroup
  create_Lref_droplist
  create_mref_droplist
  create_Tref_droplist

;----
  IF NOT(KEYWORD_SET(last_form)) THEN BEGIN
  ;take the commented lines if you want to put the guiform-files in 'output'
    ;WIDGET_CONTROL, gui.text.out_path, GET_VALUE=out_path
    ;path_length = STRLEN(out_path)
    ;IF STRMID(out_path,path_length-1,1) NE '/' THEN out_path = out_path + '/'
    ;IF path_length EQ 0 THEN out_path = 'output/'
    ;gui_filter = out_path + '*.gui'
    gui_filter ='*.gui'
    gui_files = DIALOG_PICKFILE(/READ,PATH='internal/guiforms',FILTER=gui_filter)
    IF N_ELEMENTS(gui_files) ne 1 THEN RETURN ELSE $
       IF gui_files EQ '' THEN RETURN
    gui.out.current_gui = gui_files[0]
  ENDIF ELSE BEGIN
    ;WIDGET_CONTROL, gui.text.out_path, GET_VALUE=out_path
    ;path_length = STRLEN(out_path)
    ;IF STRMID(out_path,path_length-1,1) NE '/' THEN out_path = out_path + '/'
    ;IF path_length EQ 0 THEN out_path = 'output/'
    ;gui.out.current_gui = out_path+'last_form.gui'
    gui.out.current_gui = 'internal/guiforms/'+'last_form.gui'
  ENDELSE
;----

  OPENR, guiform_lun, gui.out.current_gui, /GET_LUN, ERROR=err

  IF err NE 0 THEN BEGIN
    printerror, 'No guiform found'
    RETURN
  ENDIF

  IF !QUIET EQ 1 THEN BEGIN
    form_split = STRSPLIT(gui.out.current_gui,'/',/EXTRACT)
    form_name = form_split[N_ELEMENTS(form_split)-1]
  ENDIF ELSE form_name = gui.out.current_gui
  PRINT, 'loading ' + form_name

  tmp1_path = ''
  tmp2_path = ''
  data_path = ''
  out_path = ''
  runs_input = ''
  start_t = ''
  end_t = ''
  sparsefac = ''
  spec_str = ''
  res_steps = gui.out.res_steps
  out_format = gui.out.out_format

  READF, guiform_lun, tmp1_path
  READF, guiform_lun, tmp2_path
  READF, guiform_lun, data_path
  READF, guiform_lun, out_path
  READF, guiform_lun, runs_input
  READF, guiform_lun, start_t, end_t
  READF, guiform_lun, sparsefac
  READF, guiform_lun, tab

  IF STRMID(data_path,STRLEN(data_path)-1,1) NE '/' THEN data_path += '/'

  gui.misc.tmp1_path = tmp1_path
  gui.misc.tmp2_path = tmp2_path
  WIDGET_CONTROL, gui.text.data_path, SET_VALUE=data_path
  WIDGET_CONTROL, gui.text.out_path, SET_VALUE=out_path
  WIDGET_CONTROL, gui.text.runs, SET_VALUE=runs_input
  WIDGET_CONTROL, gui.text.start_t, SET_VALUE=start_t
  WIDGET_CONTROL, gui.text.end_t, SET_VALUE=end_t
  WIDGET_CONTROL, gui.text.sparsefac, SET_VALUE=sparsefac
  WIDGET_CONTROL, gui.droplist.Lref, SET_DROPLIST_SELECT=0
  WIDGET_CONTROL, gui.tab.base, SET_TAB_CURRENT=tab

  series.run_labels = runs_input
  gui.out.data_path = data_path

  IF runs_input NE '' THEN BEGIN
    series.par_updated = 0
    apply_form
    species = INTARR(par.n_spec > 1)
    READF, guiform_lun, species
  ENDIF ELSE BEGIN
    READF, guiform_lun, spec_str
    species = 0
  ENDELSE

  ; call diag_info again AFTER apply_form
  IF ((*gui.tab.arr)[tab].name EQ 'scan') AND $
    (FILE_TEST(data_path+'scan.log')) THEN BEGIN
    reload_diag_var_names, tab
    create_diag_table, tab=tab,/rebuild
  ENDIF ; ELSE printerror, 'no scan.log file found'

  READF, guiform_lun, res_steps
  WIDGET_CONTROL, gui.bgroup.resolve_steps, SET_VALUE=res_steps

  READF, guiform_lun, out_format
  WIDGET_CONTROL, gui.bgroup.out_format, SET_VALUE=out_format

  create_species_bgroup
  WIDGET_CONTROL, gui.bgroup.species, SET_VALUE=species

  series.Lref = 1.0D

  create_Lref_droplist
  series.mref = 1.0D
  create_mref_droplist
  series.Tref = 1.0D
  create_Tref_droplist

  ; read tables
  mode = 0
  n_vars = 0
  itab = -1
  line = ''

  WHILE NOT EOF(guiform_lun) DO BEGIN     ; loop over lines until end of file
    READF, guiform_lun, line
    CASE mode OF
      0: BEGIN                            ; diag header, selected, var count
         temp = STRSPLIT(line,'#',/EXTRACT)
         IF N_ELEMENTS(temp) NE 4 THEN BEGIN
           printerror,'form could not be restored correctly: diag tables'
           CLOSE, guiform_lun
;	   SPAWN, ['rm','-f','internal/guiform'], /NOSHELL
           clear_form, /only_tables
           RETURN
         ENDIF
         IF temp[0] NE itab THEN BEGIN
            itab = temp[0]
	    IF itab GE N_ELEMENTS((*gui.tab.arr)) THEN BEGIN
	      printerror,'form could not be restored correctly: tabs'
	      CLOSE, guiform_lun
;              SPAWN, ['rm','-f','internal/guiform'], /NOSHELL
              clear_form, /only_tables
              RETURN
            ENDIF
            diag = (*gui.tab.arr)[itab].diag0
         ENDIF
         IF PTR_VALID(diag) THEN BEGIN     ; set selected
            (*diag).selected = temp[2]
            n_vars = temp[3]
            vars = (*diag).table_entry
            v = 0
            mode = 1                        ; now variables
            IF n_vars EQ 0 THEN BEGIN       ; if no variables, go to next diag
              mode = 0
              diag = (*diag).next
            ENDIF
          ENDIF ELSE BEGIN
            printerror,'form could not be restored correctly: diags'
            CLOSE, guiform_lun
;            SPAWN, ['rm','-f','internal/guiforms/guiform'], /NOSHELL
            clear_form, /only_tables
            RETURN
          ENDELSE
        END
        1: BEGIN                                         ; variable values
          IF v LE N_ELEMENTS(vars[0,*]) - 1 THEN BEGIN
            IF vars[0,0] NE '' THEN vars[3,v] = line     ; set variable value
            n_vars = n_vars - 1                          ; count variables
            v = v + 1
            IF n_vars EQ 0 THEN BEGIN
              (*diag).table_entry = vars
              diag = (*diag).next
              mode = 0
            ENDIF
          ENDIF ELSE BEGIN
            ;printerror, 'form could not be restored correctly: vars';!!!
            print, 'form could not be restored correctly: vars'
            CLOSE, guiform_lun
;            SPAWN, ['rm','-f','internal/guiforms/guiform'], /NOSHELL
            clear_form, /only_tables
            RETURN
          ENDELSE
        END
    ENDCASE
  ENDWHILE

  FREE_LUN, guiform_lun

  RETURN

END

;###########################################################################

PRO navhist, bb=bb, back=back, forw=forw, ff=ff, add_entry=add_entry
; procedure responsible for the path/run history navigation
; always needs to be called with exactly one keyword:
; bb : fast backward (first), b : backward, f : forward,
; ff : fast forward (last), add_entry : extend list by most recent values

  COMMON global_vars

  nmax_hist = N_ELEMENTS(gui.hist.store)
  IF nmax_hist LT 3 THEN $
    printerror, 'path_navhist: too small nmax_hist', /popup
  nm_m1 = nmax_hist - 1
  nm_m2 = nmax_hist - 2

  IF KEYWORD_SET(add_entry) THEN BEGIN
    ; get rid of items where stored and current item are identical
    curr_dpath = gui.out.data_path
    curr_runs = series.run_labels
    curr_stime = gui.out.start_t
    curr_etime = gui.out.end_t
    uniq_ind = $
      WHERE((gui.hist.store.dpath NE curr_dpath) OR $
      (gui.hist.store.runs NE curr_runs) OR $
      (gui.hist.store.stime NE curr_stime) OR $
      (gui.hist.store.etime NE curr_etime))
    n_uniq = uniq_ind[0] NE -1 ? N_ELEMENTS(uniq_ind) : nmax_hist

    IF (n_uniq NE 0) AND (n_uniq NE nmax_hist) THEN BEGIN
      gui.hist.store[1:n_uniq] = gui.hist.store[uniq_ind]

      IF n_uniq NE nm_m1 THEN BEGIN
        gui.hist.store[n_uniq+1:nm_m1].dpath = ''
        gui.hist.store[n_uniq+1:nm_m1].runs = ''
        gui.hist.store[n_uniq+1:nm_m1].stime = 0.0D
        gui.hist.store[n_uniq+1:nm_m1].etime = 0.0D
      ENDIF
    ENDIF ELSE BEGIN

      IF gui.hist.nhist GT nm_m1 THEN BEGIN
        ; limit number of stored items: remove last
        gui.hist.store[1:nm_m1] = gui.hist.store[0:nm_m2]
      ENDIF ELSE BEGIN
        ; add new item
        IF gui.hist.nhist NE 0 THEN gui.hist.store[1:gui.hist.nhist] = $
          gui.hist.store[0:gui.hist.nhist-1]
        gui.hist.nhist += 1
      ENDELSE
    ENDELSE

    ; store most recent item
    gui.hist.store[0].dpath = gui.out.data_path
    gui.hist.store[0].runs = series.run_labels
    gui.hist.store[0].stime = gui.out.start_t
    gui.hist.store[0].etime = gui.out.end_t
    gui.hist.pos = 0

    RETURN
  ENDIF

  IF gui.hist.pos EQ 0 THEN BEGIN ; store current setting
    WIDGET_CONTROL, gui.text.data_path, GET_VALUE=data_path
    WIDGET_CONTROL, gui.text.runs, GET_VALUE=run_labels
    WIDGET_CONTROL, gui.text.start_t, GET_VALUE=start_t
    WIDGET_CONTROL, gui.text.end_t, GET_VALUE=end_t
    gui.out.data_path = data_path
    series.run_labels = run_labels
    gui.out.start_t = rm0es(start_t)
    gui.out.end_t = rm0es(end_t)
  ENDIF

  histpos = gui.hist.pos
  IF KEYWORD_SET(bb) THEN histpos = gui.hist.nhist - 1
  IF KEYWORD_SET(back) THEN histpos = (histpos + 1) < (gui.hist.nhist - 1)
  IF KEYWORD_SET(forw) THEN histpos = (histpos - 1) > 0
  IF KEYWORD_SET(ff) THEN histpos = 0

  navhist = gui.hist.store[histpos]
  gui.hist.pos = histpos

  gui.out.data_path = navhist.dpath
  series.run_labels = navhist.runs
  gui.out.start_t = navhist.stime
  gui.out.end_t = navhist.etime
  WIDGET_CONTROL, gui.text.data_path, SET_VALUE=navhist.dpath
  WIDGET_CONTROL, gui.text.runs, SET_VALUE=navhist.runs
  WIDGET_CONTROL, gui.text.start_t, SET_VALUE=rm0es(navhist.stime)
  WIDGET_CONTROL, gui.text.end_t, SET_VALUE=rm0es(navhist.etime)

  series.par_updated = 0
  apply_form

END

;###########################################################################

PRO show_profiles

  COMMON global_vars

  xtitle = 'x / ' + get_var_string(/rhostr)
  csize = 1.41

  set_output, 'profiles', /ps, xsize=25, ysize=17, multi=[0,par.n_spec,2]

  FOR isp = 0, par.n_spec - 1 DO BEGIN
    mymax = max([(*spec[isp].prof).temp,(*spec[isp].prof).omt])
    PLOT, (*spec[isp].prof).x, (*spec[isp].prof).temp, COLOR=1, $
      /XSTYLE, CHARSIZE=csize, YRANGE=[0,mymax], $
      TITLE=spec[isp].name, YTITLE='T, omt', XTITLE=xtitle
    OPLOT, (*spec[isp].prof).x, (*spec[isp].prof).omt, COLOR=2

    mymax = max([(*spec[isp].prof).dens,(*spec[isp].prof).omn])
    PLOT, (*spec[isp].prof).x, (*spec[isp].prof).dens, COLOR=1, $
      /XSTYLE, CHARSIZE=csize, YRANGE=[0,mymax], $
      TITLE=spec[isp].name, YTITLE='n, omn', XTITLE=xtitle
    OPLOT, (*spec[isp].prof).x, (*spec[isp].prof).omn, COLOR=2
  ENDFOR

  XYOUTS, 0.01, 0.96, '!6profiles', COLOR=1, CHARSIZE=1.75, $
    /NORMAL, CHARTHICK=1.5

  stref = (WHERE(spec[*].temp EQ 1))[0]
  snref = (WHERE(spec[*].dens EQ 1))[0]

  IF gui.out.out_format[2] EQ 1 THEN BEGIN
    rhostar_x = (SQRT(par.mref*par.Tref/par.Qref) / par.Bref) * $
      SQRT((*spec[stref].prof).temp[*]) / par.Lref
    PLOT, (*spec[0].prof).x, 1.0 / rhostar_x, COLOR=1, $
      /XSTYLE, CHARSIZE=csize, $ ;YRANGE=[0,mymax], $
      TITLE='!6', YTITLE='!61/!7q!6*', XTITLE=xtitle
    OPLOT, !X.CRANGE, 1.0 / par.rhostar * [1,1], COLOR=2

    beta_x = par.beta * (*spec[snref].prof).dens * (*spec[stref].prof).temp
    PLOT, (*spec[0].prof).x, beta_x, COLOR=1,$
      /XSTYLE, CHARSIZE=csize, $ ;YRANGE=[0,mymax], $
      TITLE='!6', YTITLE='!7b!6(x)', XTITLE=xtitle
    OPLOT, !X.CRANGE, par.beta * [1,1], COLOR=2

    coll_x = par.coll * (*spec[snref].prof).dens / (*spec[stref].prof).temp^2
    PLOT, (*spec[0].prof).x, coll_x, COLOR=1,$
      /XSTYLE, CHARSIZE=csize, $ ;YRANGE=[0,mymax], $
      TITLE='!6', YTITLE='!6coll(x)', XTITLE=xtitle
    OPLOT, !X.CRANGE, par.coll * [1,1], COLOR=2

    set_output, 'profiles', /reset, /ps
  ENDIF

END

;############################################################################

PRO series_info, summary=summary

  COMMON global_vars

  IF file_par.exist EQ 0 THEN RETURN

  simple_geom = (par.magn_geometry EQ '') OR $
    (par.magn_geometry EQ 's-alpha') OR (par.magn_geometry EQ 's_alpha_B') OR $
    (par.magn_geometry EQ 'slab') OR (par.magn_geometry EQ 'circular')

  PLOT, [0,1], [0,1], XSTYLE=12, YSTYLE=12, XMARGIN=[0,0], YMARGIN=[0,0], /NODATA
  PLOT, [0,0.93], [0.93,0.93], XSTYLE=12, YSTYLE=12, XMARGIN=[0,0], YMARGIN=[0,0]
  XYOUTS, 0.01, 0.95, '!5runs ' + series.run_labels, COLOR=1, CHARSIZE=1.7

  ch_size = 1.9

  ; --- species table
  listentry = STRARR(1+par.n_spec,6)
  listentry[0,*] = ['!7x!5!Dn!N','!7x!5!DT!N','!5m','!5   q','!5T','!5n']
  spec_str = '!5species'

  FOR isp = 0, par.n_spec-1 DO BEGIN
    listentry[isp+1,*] = rm0es([spec[isp].omn,spec[isp].omt,spec[isp].mass,$
      spec[isp].charge,spec[isp].temp,spec[isp].dens],prec=3)
    listentry[isp+1,2] = rm0es(STRING(spec[isp].mass,FORMAT='(G12.2)'))
    listentry[isp+1,3] = '   ' + listentry[isp+1,3]
    spec_str = [spec_str,spec[isp].name]
  ENDFOR

  n_columns = N_ELEMENTS(listentry[0,*])
  n_rows = N_ELEMENTS(listentry[*,0])
  ymax = (par.n_spec + 1) * 0.06 < 0.3
  ymin = 0.01
  xmin = 0.2
  xmax = 1.0
  xstep = (xmax - xmin) / n_columns
  ystep = (ymax - ymin) / (n_rows - 1) < 0.06

  FOR row=0, n_rows-1 DO BEGIN
    ycoord = ymax - ystep * row
    XYOUTS, 0.01, ycoord, spec_str[row],CHARSIZE=ch_size, COLOR=1
    FOR column=0, n_columns-1 DO BEGIN
      xcoord = xmin + ((xmax - xmin) / n_columns) * column
      XYOUTS, xcoord, ycoord, $
        listentry[row,column], CHARSIZE=ch_size, COLOR=1
    ENDFOR
  ENDFOR

  ; --- general parameters
  listentry = STRARR(2,18)
  listentry[0,*] = ['!5geometry','!5type','!5n!Dkx0!N,!5n!Dky0!N,!5n!Dz0!N',$
    '!5n!Dv!N,!5n!Dw!N','!5L!Dx!N','!5k!Dx,min!N',$
    '!5L!Dy!N', '!5k!Dy,min!N','!5L!Dv!N,!5L!Dw!N',$
    simple_geom?'!5r/R!D0!N':'!5x!D0!N','!5q!D0!N',$
    par.amhd GT 0 ? '!S!5s!R!7!U^!N,!7a!5!DMHD!N' : '!S!5s!R!7!U^!N',$
    '!5R!D0!N/L!Dref!N','!51/!7q!5!U*!N','!7b','!7m!5!Dcoll!N',$
    '!7c!5!DExB!N','!5D!D(x,y,z,v!D!9#!5)!N']
  listentry[1,*] = rm0es([0,0,par.nkx0,par.nv0,$
    par.lx,0,par.ly,par.kymin,par.lv,$
    simple_geom?par.trpeps:par.x0,par.q0,$
    par.shat,par.major_R,1./par.rhostar,$
    par.beta,par.coll,par.ExBrate,0],prec=3)

  IF (par.magn_geometry EQ '') OR (par.magn_geometry EQ 's_alpha') THEN $
    listentry[1,0] = '!5s-!7a!5' ELSE $
    listentry[1,0] = (STRSPLIT(par.magn_geometry,'.dat',/EXTRACT,/REGEX))[0]
  listentry[1,1] = par.nonlinear ? 'nonlinear' : 'linear'
  IF par.delzonal THEN listentry[1,1] += ', !sZF!r--!N'
  listentry[1,2] = rm0es(par.nkx0) + ',' + rm0es(par.nky0-par.ky0_ind)
  IF par.ky0_ind GT 0 THEN listentry[1,2] += '(!8+' + rm0es(par.ky0_ind) + '!5)'
  listentry[1,2] += ',' + rm0es(par.nz0)
  listentry[1,3] = rm0es(par.nv0) + ',' + rm0es(par.nw0)
  IF par.lx EQ 0 THEN BEGIN
    listentry[1,4] = '1/(!S!5s!R!7!U^!N!5 k!Dy!N)'
    listentry[1,5] = '2!7p!5/L!Dx'
  ENDIF ELSE listentry[1,5] = rm0es(2.0*!PI/par.lx,prec=4)
  listentry[1,8] = rm0es(par.lv) + ',' + rm0es(par.lw)
  IF par.n0_global GT 0 THEN $
    listentry[*,7] += [' !5(n!I0!N)',' ('+rm0es(par.n0_global)+')']
  IF (par.major_R EQ 0.0) AND (par.q0 EQ 0) THEN BEGIN
     IF (par.ehat NE 0.0) THEN BEGIN
        listentry[0,10] = '!S!7e!R!U^!N'
        listentry[1,10] = rm0es(par.ehat,prec=3)
     ENDIF ELSE listentry[0:1,10] = ''
  ENDIF
  IF par.amhd GT 0 THEN listentry[1,11] += rm0es(par.amhd,prec=3)
  IF (par.major_R GT 0) AND (par.minor_R GT 0) THEN BEGIN
    listentry[0,12] = '(!5R!D0!N,!7a!5)/L!Dref!N'
    listentry[1,12] = rm0es(par.major_R) + ',' + rm0es(par.minor_r)
  ENDIF ELSE IF par.minor_r GT 0 THEN BEGIN
    listentry[0,12] = '!7a!5/L!Dref!N'
    listentry[1,12] = rm0es(par.minor_r)
  ENDIF ELSE IF par.major_R EQ 0 THEN listentry[0:1,12] = ''
  IF (par.pfsrate GT 0) THEN $
    listentry[*,16] += [',!7c!5!Dpfs!N',','+rm0es(par.pfsrate,prec=3)]
  listentry[1,17] = rm0es(par.hyp_x) + ',' + rm0es(par.hyp_y) + $
    ',' + rm0es(par.hyp_z) + ',' + rm0es(par.hyp_v)

  IF NOT par.x_local THEN BEGIN
     listentry[*,1] += [' [b.c.]',' ['+rm0es(par.rad_bc_type)+']']
  ENDIF

  n_entries = N_ELEMENTS(listentry[0,*])
  ymin = 0.4
  ymax = 0.85
  FOR column = 0, 1 DO BEGIN
    FOR row = 0, n_entries / 2 - 1 DO BEGIN
      XYOUTS, 0.01 + column * 0.55, $
        ymax - row * (ymax - ymin) / (n_entries / 2 - 1), $
        listentry[0,column*n_entries/2+row], CHARSIZE=ch_size, COLOR=1
      XYOUTS, 0.25 + column * 0.5, $
        ymax - row * (ymax - ymin) / (n_entries / 2 - 1), $
        listentry[1,column*n_entries/2+row], CHARSIZE=ch_size, COLOR=1
    ENDFOR
  ENDFOR

END

;###########################################################################

PRO refvalues_info

  COMMON global_vars

  IF file_par.exist EQ 0 THEN RETURN

  PLOT, [0,1], [0,1], XSTYLE=12, YSTYLE=12, XMARGIN=[0,0], YMARGIN=[0,0], /NODATA
  PLOT, [0,0.93], [0.93,0.93], XSTYLE=12, YSTYLE=12, XMARGIN=[0,0], YMARGIN=[0,0]
  XYOUTS, 0.01, 0.95, '!5runs ' + series.run_labels, COLOR=1, CHARSIZE=1.7

  nvals = 8

  list_entry = STRARR(2,nvals)
  list_entry[0,*] = ['L!Dref!N / m   :','B!Dref!N / T    :','Q!Dref!N / C    :',$
    'm!Dref!N / kg  :','T!Dref!N / eV   :','n!Dref!N        :','v!Dref!N / (m/s) :',$
    '!7q!6!Dref!N / m    :']

  list_entry[1,*] = [par.Lref,par.Bref,par.Qref,par.mref,par.Tref,par.nref,$
    SQRT(par.Tref*par.Qref/par.mref),SQRT(par.mref*par.Tref/par.Qref)/par.Bref]

  chsize = 2.0
  FOR j = 0, 7 DO BEGIN
    XYOUTS, 0.1, 0.85 - j * (0.85 / nvals), list_entry[0,j], $
      COLOR=1, CHARSIZE=chsize
    XYOUTS, 0.3, 0.85 - j * (0.85 / nvals), list_entry[1,j], $
      COLOR=1, CHARSIZE=chsize
  ENDFOR

END

;############################################################################

PRO display_var_list

  COMMON global_vars

  list_entry = STRARR(2,par.n_fields+par.n_moms)

  list_entry[0,*] = STRTRIM(STRING(INDGEN(par.n_fields+par.n_moms)),2)
  list_entry[1,*] = get_var_string(INDGEN(par.n_fields+par.n_moms),/units)

  chsize = 1.8

  PLOT, [0,1], [0,1], XSTYLE=12, YSTYLE=12, XMARGIN=[0,0], YMARGIN=[0,0], /NODATA
  FOR j = 0, par.n_fields + par.n_moms - 1 DO BEGIN
    XYOUTS, 0.1, 0.95 - j * (0.95 / (par.n_fields + par.n_moms)), $
      list_entry[0,j], COLOR=1, CHARSIZE=chsize
    XYOUTS, 0.3, 0.95 - j * (0.95 / (par.n_fields + par.n_moms)), $
      list_entry[1,j], COLOR=1, CHARSIZE=chsize
  ENDFOR

END

;############################################################################

PRO display_help, diag
; reads the list with the available diagnostics and converts
; information to the internal data structure (name, help, variables)

  COMMON global_vars

  gui.window.help = WIDGET_BASE($
    GROUP_LEADER = gui.window.main,$
    TITLE='diagnostic help: '+(*diag).name,$
    SCR_XSIZE = 550,$
    SCR_YSIZE = 340,$
    XOFFSET   = 210,$
    YOFFSET   = 20,$
    EVENT_PRO = 'gui_main_event')

  temp = WIDGET_TEXT($
    gui.window.help,$
    SCR_XSIZE = 540,$
    SCR_YSIZE = 150,$
    XOFFSET   = 10,$
    YOFFSET   = 10,$
    /SCROLL,$
    /WRAP,$
    VALUE     = (*diag).help_text)

  temp = WIDGET_LABEL($
    gui.window.help,$
    VALUE     = 'parameters:',$
    XOFFSET   = 10,$
    YOFFSET   = 160)

  vars = (*diag).table_entry
  n_vars = N_ELEMENTS(vars[0,*])
  var_buffer = STRARR(n_vars)
  max_length = 0
  IF vars[0,0] NE '' THEN BEGIN
    FOR v = 0, n_vars - 1 DO BEGIN
      name_length = STRLEN(vars[0,v])
      max_length = name_length > max_length
    ENDFOR
    FOR v = 0, n_vars - 1 DO BEGIN
      var_str = vars[0,v]
      WHILE STRLEN(var_str) LT max_length DO var_str = var_str + ' '
      var_str = var_str + ':' + vars[2,v]
      var_buffer[v] = var_str
    ENDFOR

    temp = WIDGET_TEXT($
      gui.window.help,$
      SCR_XSIZE = 540,$
      SCR_YSIZE = 150,$
      XOFFSET   = 10,$
      YOFFSET   = 190,$
      /SCROLL,$
      VALUE     = var_buffer)
  ENDIF

  WIDGET_CONTROL, gui.window.help, /REALIZE

END

;############################################################################

PRO custom_diag_list

  COMMON global_vars

  IF gui.misc.n_cdiags LT 1 THEN BEGIN
    printerror, 'no custom diags found'
    RETURN
  ENDIF

  gui.window.cdiags = WIDGET_BASE($
    GROUP_LEADER = gui.window.main,$
    TITLE='custom diagnostics',$
    SCR_XSIZE = 550,$
    SCR_YSIZE = 340,$
    XOFFSET   = 210,$
    YOFFSET   = 20,$
    EVENT_PRO = 'gui_main_event')

  tab_value = STRARR(2,gui.misc.n_cdiags)

  FOR cdiag_nr = 0, gui.misc.n_cdiags - 1 DO BEGIN
    cdiag = get_cdiag_pointer(cdiag_nr)
    tab_value[0,cdiag_nr] = (*cdiag).name
    IF (*cdiag).selected EQ 1 THEN tab_value[1,cdiag_nr] = 'X' $
      ELSE tab_value[1,cdiag_nr] = ''
  ENDFOR

  ; create table widget
  gui.table.cdiags = WIDGET_TABLE($
    gui.window.cdiags,$
    COLUMN_WIDTHS = [180,18],$
    XOFFSET   = 10,$
    YOFFSET   = 10,$
    SCR_XSIZE = 530,$
    SCR_YSIZE = 320,$
    XSIZE     = 2,$
    YSIZE     = gui.misc.n_cdiags,$
    VALUE     = tab_value,$
    COLUMN_LABELS = ['',''],$
    /ALL_EVENTS,$
    EDITABLE=[0,0],$
    /NO_HEADERS)

  WIDGET_CONTROL, gui.window.cdiags, /REALIZE
  XMANAGER, 'custom_diag_list', gui.window.cdiags

END

;############################################################################

PRO custom_diag_list_event, ev

  COMMON global_vars

  CASE ev.id OF

    ; --- events of the custom diag table, toggle custom diags on/off ---
    gui.table.cdiags : IF (ev.type EQ 4) AND (ev.sel_left EQ 1) THEN BEGIN
      WIDGET_CONTROL, gui.table.cdiags, GET_VALUE=diag_tab
      cdiag = get_cdiag_pointer(ev.sel_top)
      IF(PTR_VALID(cdiag)) THEN BEGIN
        IF diag_tab[1,ev.sel_top] EQ 'X' THEN BEGIN
          diag_tab[1,ev.sel_top]=' '
          (*cdiag).selected = 0
        ENDIF ELSE BEGIN
          diag_tab[1,ev.sel_top]='X'
          (*cdiag).selected = 1
        ENDELSE
      ENDIF
      WIDGET_CONTROL, gui.table.cdiags, SET_VALUE=diag_tab
    ENDIF

    ELSE:
  ENDCASE

END

;############################################################################

PRO clear_form,  only_tables=only_tables

  COMMON global_vars

  printerror, ''

  FOR itab = 0, N_ELEMENTS(*gui.tab.arr) - 1 DO BEGIN
    diag = (*gui.tab.arr)[itab].diag0
    WHILE PTR_VALID(diag) DO BEGIN        ; loop over all diagnostics
      (*diag).selected = 0                ; deselect diagnostic
      vars = (*diag).table_entry        ; loop over all parameters
      IF vars[0,0] NE '' THEN BEGIN
        FOR v = 0, N_ELEMENTS(vars[0,*]) - 1 DO $
          IF vars[1,v] EQ 1 THEN vars[3,v] = '0' ELSE vars[3,v] = ''
      ENDIF
      diag = (*diag).next
    ENDWHILE
  ENDFOR
  create_diag_table, /refresh

  IF NOT KEYWORD_SET(only_tables) THEN BEGIN
    WIDGET_CONTROL, gui.text.data_path, SET_VALUE=''    ; clear input boxes
    WIDGET_CONTROL, gui.text.out_path, SET_VALUE=''
    WIDGET_CONTROL, gui.text.runs, SET_VALUE=''
    WIDGET_CONTROL, gui.text.start_t, SET_VALUE=''
    WIDGET_CONTROL, gui.text.end_t, SET_VALUE=''
    WIDGET_CONTROL, gui.text.sparsefac, SET_VALUE='1'
    WIDGET_CONTROL, gui.bgroup.species, SET_VALUE=0     ; reset buttons
    WIDGET_CONTROL, gui.bgroup.resolve_steps, SET_VALUE=0
    WIDGET_CONTROL, gui.bgroup.out_format, SET_VALUE=[1,0]
    WIDGET_CONTROL, gui.droplist.Lref, SET_DROPLIST_SELECT=0
    WIDGET_CONTROL, gui.tab.base, SET_TAB_CURRENT=0

    gui.misc.tmp1_path = ''
    gui.misc.tmp2_path = ''
    series.Lref = 1.0D ; required because of above statement
  ENDIF

END

;############################################################################

PRO save_form, save_as=save_as, last_form=last_form

  COMMON global_vars

  IF FILE_SEARCH('internal/guiforms', /TEST_DIRECTORY) EQ '' THEN $
    SPAWN, 'mkdir' + ' ' + 'internal/guiforms'

  IF KEYWORD_SET(save_as) THEN BEGIN
    gui_filter = '*.gui'
    gui_files = DIALOG_PICKFILE(/READ,PATH='internal/guiforms',FILTER=gui_filter)
    IF (N_ELEMENTS(gui_files) NE 1) OR (gui_files EQ '') THEN RETURN
    gui.out.current_gui = gui_files[0]
    IF (STRMID(gui.out.current_gui,STRLEN(gui.out.current_gui)-4,4) NE '.gui') $
      THEN gui.out.current_gui += '.gui'
   ENDIF

  IF KEYWORD_SET(last_form) THEN gui.out.current_gui = 'internal/guiforms/last_form.gui'

  OPENW, guiform_lun, gui.out.current_gui, /GET_LUN, ERROR=err
  IF err NE 0 THEN BEGIN
    message = DIALOG_MESSAGE('/guiform_lun could not be created.',$
      /ERROR,DIALOG_PARENT=gui.window.main)
    RETURN
  ENDIF

  PRINT, 'Saving guiform in '+ gui.out.current_gui

  WIDGET_CONTROL, gui.text.data_path, GET_VALUE=data_path
  IF STRMID(data_path,STRLEN(data_path)-1,1) NE '/' THEN $
    data_path += '/'
  WIDGET_CONTROL, gui.text.out_path, GET_VALUE=out_path
  WIDGET_CONTROL, gui.text.runs, GET_VALUE=runs_input
  WIDGET_CONTROL, gui.text.start_t, GET_VALUE=start_t
  WIDGET_CONTROL, gui.text.end_t, GET_VALUE=end_t
  WIDGET_CONTROL, gui.text.sparsefac, GET_VALUE=sparsefac
  WIDGET_CONTROL, gui.bgroup.species, GET_VALUE=species
  WIDGET_CONTROL, gui.bgroup.resolve_steps, GET_VALUE=res_steps
  WIDGET_CONTROL, gui.bgroup.out_format, GET_VALUE=out_format

  ; return to L_ref normalization
  start_t = start_t / time_renorm(0)
  end_t = end_t / time_renorm(0)

  PRINTF, guiform_lun, gui.misc.tmp1_path                 ; write to file
  PRINTF, guiform_lun, gui.misc.tmp2_path
  PRINTF, guiform_lun, data_path
  PRINTF, guiform_lun, out_path
  PRINTF, guiform_lun, runs_input
  PRINTF, guiform_lun, STRTRIM(STRING(start_t),2),STRTRIM(STRING(end_t),2)
  PRINTF, guiform_lun, sparsefac
  PRINTF, guiform_lun, WIDGET_INFO(gui.tab.base,/TAB_CURRENT)
  PRINTF, guiform_lun, species
  PRINTF, guiform_lun, res_steps
  PRINTF, guiform_lun, out_format

  FOR itab = 0, N_ELEMENTS((*gui.tab.arr)) - 1 DO BEGIN  ; save tables
    ; reset scan tab if no scanlog found
    IF itab EQ WIDGET_INFO(gui.tab.base,/TAB_CURRENT) THEN BEGIN
      IF ((*gui.tab.arr)[itab].name EQ 'scan') AND $
        (NOT FILE_TEST(data_path+'scan.log')) AND $
        ((*gui.tab.arr)[itab].n_params GT 0) THEN BEGIN

        reload_diag_var_names, itab
        create_diag_table, tab=itab, /rebuild
      ENDIF
    ENDIF

    diag = (*gui.tab.arr)[itab].diag0
    diag_count = 0
    WHILE PTR_VALID(diag) DO BEGIN
      vars = (*diag).table_entry
      IF vars[0,0] EQ '' THEN n_vars = 0 ELSE n_vars = N_ELEMENTS(vars[0,*])

      ; itab # diag number # 1 if selected # number of variables
      PRINTF, guiform_lun, itab, ' # ',diag_count, ' # ', (*diag).selected, $
        ' # ', n_vars

      FOR v = 0, n_vars - 1 DO PRINTF, guiform_lun, vars[3,v] ; save var value
      diag = (*diag).next
      diag_count = diag_count + 1
    ENDWHILE
  ENDFOR

  FREE_LUN, guiform_lun

END

;#############################################################################

PRO show_h5_files

  COMMON global_vars

  WIDGET_CONTROL, gui.text.out_path, GET_VALUE=out_path
  path_length = STRLEN(out_path)
  IF STRMID(out_path,path_length-1,1) NE '/' THEN out_path = out_path + '/'
  IF path_length EQ 0 THEN out_path = 'output/'
  h5_filter = out_path + '*.h5'
  h5_files = $
    DIALOG_PICKFILE(/READ,FILTER=h5_filter,/MULTIPLE_FILE)
  IF N_ELEMENTS(h5_files) EQ 1 THEN IF h5_files EQ '' THEN RETURN
  FOR f = 0, N_ELEMENTS(h5_files) - 1 DO dummy = H5_BROWSER(h5_files[f])

END

;#############################################################################

PRO show_ps_files

  COMMON global_vars

  WIDGET_CONTROL, gui.text.out_path, GET_VALUE=out_path
  path_length = STRLEN(out_path)
  IF STRMID(out_path,path_length-1,1) NE '/' THEN out_path = out_path + '/'
  IF path_length EQ 0 THEN out_path = 'output/'
  ps_filter = out_path + '*.ps;*.eps'
  ps_files = $
    DIALOG_PICKFILE(/READ,FILTER=ps_filter,/MULTIPLE_FILES)
  IF N_ELEMENTS(ps_files) EQ 1 THEN IF ps_files EQ '' THEN RETURN
  FOR f = 0, N_ELEMENTS(ps_files) - 1 DO $
    SPAWN, cfg.ps_viewer + ps_files[f] + '&'

END

;#############################################################################

PRO check_bgroups, no_reset=no_reset

  COMMON global_vars

  WIDGET_CONTROL, gui.bgroup.species, GET_VALUE=species
  IF NOT KEYWORD_SET(no_reset) THEN BEGIN ; if no species selected, select all
    IF TOTAL(species) EQ 0 THEN species = species + 1
    WIDGET_CONTROL, gui.bgroup.species, SET_VALUE=species
  ENDIF

  spec_select = INTARR(TOTAL(species)>1)
  k = 0
  FOR j = 0, N_ELEMENTS(species) - 1 DO IF species[j] EQ 1 THEN BEGIN
    spec_select[k] = j
    k = k + 1
  ENDIF
  IF PTR_VALID(gui.out.spec_select) THEN PTR_FREE, gui.out.spec_select
  gui.out.spec_select = PTR_NEW(spec_select)
  gui.out.n_spec_sel = TOTAL(species)

  WIDGET_CONTROL, gui.bgroup.resolve_steps, GET_VALUE=res_steps
  gui.out.res_steps = res_steps

  WIDGET_CONTROL, gui.bgroup.out_format, GET_VALUE=out_format
  IF NOT KEYWORD_SET(no_reset) AND (TOTAL(out_format) EQ 0) THEN $
    PRINT, 'WARNING: no output format selected'
  gui.out.out_format = out_format

END

;#############################################################################

PRO gui_main_event, ev, dummy=dummy
; event procedure for widget application
; ev: event passed by the widget system

  COMMON global_vars

  CASE ev.id OF
    gui.button.nrg_data : BEGIN              ; click on show nrg data button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
      IF par.n_spec GT 0 THEN BEGIN
        IF NOT PTR_VALID(nrg) THEN read_nrg
        plot_nrg_data
      ENDIF ELSE BEGIN
        PTR_FREE, nrg
	printerror, 'cannot show nrg data without parameters'
      ENDELSE
    END

    ; data path and run label navigation
    gui.button.data_path : BEGIN
        WIDGET_CONTROL, gui.text.data_path, GET_VALUE=start_path
        data_path = DIALOG_PICKFILE(/DIRECTORY,/MUST_EXIST,PATH=start_path,$
          TITLE='Please select a data path')
        IF (data_path NE '') THEN WIDGET_CONTROL, $
          gui.text.data_path, SET_VALUE=data_path
      END

    gui.button.path_bbnav : navhist, /bb

    gui.button.path_backnav : navhist, /back

    gui.button.path_forwnav : navhist, /forw

    gui.button.path_ffnav : navhist, /ff

    ; get user-defined text in text-fields
    gui.button.tmp1_path : IF gui.misc.tmp1_path NE '' THEN $
      WIDGET_CONTROL, gui.text.data_path, SET_VALUE=gui.misc.tmp1_path $
      ELSE gui.misc.tmp1_path = $
      DIALOG_PICKFILE(/DIRECTORY,TITLE='Please select a default for button tmp1')

    gui.button.tmp2_path : IF gui.misc.tmp2_path NE '' THEN $
      WIDGET_CONTROL, gui.text.data_path, SET_VALUE=gui.misc.tmp2_path $
      ELSE gui.misc.tmp2_path = $
      DIALOG_PICKFILE(/DIRECTORY,TITLE='Please select a default for button tmp2')

    gui.button.out_path : BEGIN
        WIDGET_CONTROL, gui.text.out_path, GET_VALUE=start_path
        IF STRTRIM(start_path) EQ '' THEN start_path=GETENV('PWD')+'/output'
        out_path = DIALOG_PICKFILE(/DIRECTORY,PATH=start_path,$
          TITLE='Please select an output directory')
        IF (out_path NE '') THEN $
          WIDGET_CONTROL, gui.text.out_path, SET_VALUE=out_path
      END

    gui.button.runs : BEGIN
      WIDGET_CONTROL, gui.text.data_path, GET_VALUE=start_path
      runstr_in = DIALOG_PICKFILE(PATH=start_path,FILTER='parameters*',$
        /FIX_FILTER,/MUST_EXIST,/MULTIPLE_FILES,$
        TITLE='Please select a run')
      IF runstr_in[0] NE '' THEN BEGIN
        runstr_out = ''
        ; check for path changes in run pick dialog
        cur_path = (STRSPLIT(runstr_in[0],'parameters',/REGEX,/EXTRACT))[0]
        IF (cur_path NE start_path) THEN $
          WIDGET_CONTROL, gui.text.data_path, SET_VALUE=cur_path

        ; erase path/parameters
        add_act = 0
        FOR istr = 0, N_ELEMENTS(runstr_in) - 1 DO BEGIN
          runstr_in[istr] = (STRSPLIT(runstr_in[istr],'parameters',$
            /REGEX,/EXTRACT,COUNT=nsplit))[nsplit-1]
          IF runstr_in[istr] EQ '.dat' THEN add_act = 1 ELSE BEGIN
            runstr_out += STRMID(runstr_in[istr],1)
            IF istr NE (N_ELEMENTS(runstr_in) - 1) THEN runstr_out += ','
          ENDELSE
        ENDFOR
        IF add_act THEN runstr_out += $
           ((N_ELEMENTS(runstr_in) GT 1)?',':'')+'act'
        WIDGET_CONTROL, gui.text.runs, SET_VALUE=runstr_out
        apply_form
      ENDIF
    END

    gui.text.runs : IF ev.type EQ 0 THEN apply_form

    gui.button.first_step : BEGIN      ; click on first button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
      func_name = 'get_'+(STRSPLIT((*gui.tab.arr)[WIDGET_INFO(gui.tab.base,$
        /TAB_CURRENT)].name,' ',/EXTRACT))[0]+'_time'
      t_min_arr = DBLARR(series.n_runs)
      FOR run = 0, series.n_runs - 1 DO $
        t_min_arr[run] = CALL_FUNCTION(func_name,run,/first)
      gui.out.start_t = MIN(t_min_arr)
      WIDGET_CONTROL, gui.text.start_t, SET_VALUE=rm0es(gui.out.start_t)
      IF NOT PTR_VALID(nrg) THEN read_nrg ;due to cursor
      plot_nrg_data ;due to cursor
    END

    gui.button.last_step : BEGIN      ; click on last button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
      func_name = 'get_'+(STRSPLIT((*gui.tab.arr)[WIDGET_INFO(gui.tab.base,$
        /TAB_CURRENT)].name,' ',/EXTRACT))[0]+'_time'
      t_max_arr = DBLARR(series.n_runs)
      FOR run = 0, series.n_runs - 1 DO $
        t_max_arr[run] = CALL_FUNCTION(func_name,run,/last)
      gui.out.end_t = MAX(t_max_arr) ;> 0.0
      WIDGET_CONTROL, gui.text.end_t, SET_VALUE=rm0es(gui.out.end_t)
      IF NOT PTR_VALID(nrg) THEN read_nrg ;due to cursor
      plot_nrg_data ;due to cursor
    END

    gui.bgroup.species :

    gui.bgroup.resolve_steps :

    gui.bgroup.out_format :

    gui.droplist.Lref : BEGIN          ; select in Lref droplist
      apply_form
      ; reset nrg data to L_ref
      IF PTR_VALID(nrg) THEN BEGIN
        (*nrg).time *= time_renorm(-1)
        FOR fspec = 0, 8 * par.n_spec - 1 DO $
          (*nrg)[*].data[fspec] /= nrg_renorm(fspec)
      ENDIF

      gui.out.start_t *= time_renorm(-1)
      gui.out.end_t *= time_renorm(-1)
      ; get the chosen reference
      gui.out.Lref_sel = WIDGET_INFO(gui.droplist.Lref,/DROPLIST_SELECT)
      series.Lref = set_Lref(gui.out.Lref_sel)

      ; now the actual renormalization
      IF PTR_VALID(nrg) THEN BEGIN
        (*nrg).time *= nrg_renorm(0,/time)
        FOR fspec = 0, 8 * par.n_spec - 1 DO $
          (*nrg)[*].data[fspec] *= nrg_renorm(fspec)
      ENDIF
      WIDGET_CONTROL, gui.text.start_t, $
        SET_VALUE=rm0es(gui.out.start_t*time_renorm(0))
      WIDGET_CONTROL, gui.text.end_t, $
        SET_VALUE=rm0es(gui.out.end_t*time_renorm(0))
    END

    gui.text.start_t : BEGIN    ; required because cursor with 'show nrg data' image
                                ; should be plotted just in time with
                                ; setting of new cursor position
      start_t = 0.0D
      WIDGET_CONTROL, gui.text.start_t, GET_VALUE=start_t
      gui.out.start_t = rm0es(start_t)
      WIDGET_CONTROL, gui.text.start_t, SET_VALUE=rm0es(start_t)
      IF NOT PTR_VALID(nrg) THEN read_nrg
      plot_nrg_data
    END

    gui.text.end_t : BEGIN      ; required because cursor with 'show nrg data' image
                                ; should be plotted just in time with
                                ; setting of new cursor position
      end_t = 0.0D
      WIDGET_CONTROL, gui.text.end_t, GET_VALUE=end_t
      gui.out.end_t = rm0es(end_t)
      WIDGET_CONTROL, gui.text.end_t, SET_VALUE=rm0es(end_t)
      IF NOT PTR_VALID(nrg) THEN read_nrg
      plot_nrg_data
    END

    gui.text.sparsefac :

    gui.droplist.mref : BEGIN
      apply_form
      ; reset nrg data to L_perp
      IF PTR_VALID(nrg) THEN BEGIN
        (*nrg).time = (*nrg).time / time_renorm(0)
	FOR fspec = 0, 8 * par.n_spec - 1 DO (*nrg)[*].data[fspec] = $
          (*nrg)[*].data[fspec] / nrg_renorm(fspec)
      ENDIF
      gui.out.start_t = gui.out.start_t / time_renorm(0)
      gui.out.end_t = gui.out.end_t / time_renorm(0)
      mref_sel = WIDGET_INFO(gui.droplist.mref,/DROPLIST_SELECT)
      mref_found = ((WHERE(spec.mass EQ 1.0))[0] NE -1)
      IF mref_found THEN BEGIN
        series.mref = 1.0 / spec[mref_sel].mass
	series.mref_str = 'm!D' + spec[mref_sel].name + '!N'
      ENDIF ELSE IF mref_sel EQ 0 THEN BEGIN
        series.mref=1.0D
	series.mref_str = 'm!Dref!N'
      ENDIF ELSE BEGIN
        series.mref = 1.0 / spec[mref_sel-1].mass
	series.mref_str = 'm!D' + spec[mref_sel-1].name + '!N'
      ENDELSE
      IF PTR_VALID(nrg) THEN BEGIN
        (*nrg).time = (*nrg).time * nrg_renorm(0,/time)
        FOR fspec = 0, 8 * par.n_spec - 1 DO $
          (*nrg)[*].data[fspec] = (*nrg)[*].data[fspec] * nrg_renorm(fspec)
      ENDIF
      WIDGET_CONTROL, gui.text.start_t, $
        SET_VALUE=rm0es(gui.out.start_t*time_renorm(0))
      WIDGET_CONTROL, gui.text.end_t, $
        SET_VALUE=rm0es(gui.out.end_t*time_renorm(0))
    END

    gui.droplist.Tref : BEGIN
      apply_form
      ; reset nrg data to Tref
      IF PTR_VALID(nrg) THEN BEGIN
        (*nrg).time = (*nrg).time / time_renorm(0)
	FOR fspec = 0, 8 * par.n_spec - 1 DO (*nrg)[*].data[fspec] = $
          (*nrg)[*].data[fspec] / nrg_renorm(fspec)
      ENDIF
      gui.out.start_t = gui.out.start_t / time_renorm(0)
      gui.out.end_t = gui.out.end_t / time_renorm(0)
      Tref_sel = WIDGET_INFO(gui.droplist.Tref,/DROPLIST_SELECT)
      Tref_found = ((WHERE(spec.temp EQ 1.0))[0] NE -1)
      IF Tref_found THEN BEGIN
        series.Tref = 1.0 / spec[Tref_sel].temp
	series.Tref_str = 'T!D' + spec[Tref_sel].name + '!N'
      ENDIF ELSE IF Tref_sel EQ 0 THEN BEGIN
        series.Tref = 1.0D
	series.Tref_str = 'T!Dref!N'
      ENDIF ELSE BEGIN
        series.Tref = 1.0 / spec[Tref_sel-1].temp
	series.Tef_str = 'T!D' + spec[Tref_sel-1] + '!N'
      ENDELSE
      IF PTR_VALID(nrg) THEN BEGIN
        (*nrg).time = (*nrg).time * nrg_renorm(0,/time)
        FOR fspec = 0, 8 * par.n_spec - 1 DO (*nrg)[*].data[fspec] = $
          (*nrg)[*].data[fspec] * nrg_renorm(fspec)
      ENDIF
      WIDGET_CONTROL, gui.text.start_t, $
        SET_VALUE=rm0es(gui.out.start_t*time_renorm(0))
      WIDGET_CONTROL, gui.text.end_t, $
        SET_VALUE=rm0es(gui.out.end_t*time_renorm(0))
    END

;    gui.button.nrg_diag : BEGIN                 ; click on nrg diag button
;      WIDGET_CONTROL, /HOURGLASS
;      apply_form
;      IF par.n_spec GT 0 THEN BEGIN
;        gui.misc.recent_h5 = ''
;        gui.misc.recent_ps = ''
;        IF (series.par_updated EQ 0) OR NOT PTR_VALID(nrg) THEN BEGIN
;          read_nrg
;        ENDIF
;        check_bgroups
;        plot_nrg_data, /diag
;      ENDIF ELSE BEGIN
;      	PTR_FREE, nrg
;	printerror, 'cannot run nrg diagnostic without parameters'
;       ENDELSE
;      gui.out.start_t = gui.out.start_t/time_renorm(0) ; to return to the default renormalized time values,
;      gui.out.end_t = gui.out.end_t/time_renorm(0)     ; which could be destroyed by apply_form
;    END

    gui.button.series_info : BEGIN               ; click on series_info button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
;      IF NOT series.par_updated THEN read_par
      series_info
    END

    gui.button.refvalues_info : BEGIN            ; click on ref-values_info button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
      refvalues_info
    END

    gui.button.geometry : BEGIN                  ; click on geometry button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
      gui.misc.recent_h5 = ''
      gui.misc.recent_ps = ''
      IF (series.par_updated EQ 0) OR NOT PTR_VALID(series.geom) THEN $
        read_geometry
      show_geometry_coeffs
    END

    gui.button.profiles : BEGIN                  ; click on profile button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
      gui.misc.recent_h5 = ''
      gui.misc.recent_ps = ''
      read_profiles
      show_profiles
    END

    gui.button.var_list : BEGIN                  ; click on varlist button
      WIDGET_CONTROL, /HOURGLASS
      apply_form
;      IF series.par_updated EQ 0 THEN read_par
      display_var_list
    END

    gui.button.start : BEGIN ; click on start button -- saving guiform
      WIDGET_CONTROL, /HOURGLASS
      printerror, ''
      current_tab = WIDGET_INFO(gui.tab.base,/TAB_CURRENT)
      save_form, /last_form
      apply_form
      navhist, /add_entry
      IF par.n_spec GT 0 THEN BEGIN
        gui.misc.recent_h5 = ''
        gui.misc.recent_ps = ''
        check_bgroups
        XYOUTS, 0.3, 0.0, '!3click and hold in this window to abort', $
          /NORMAL, CHARSIZE=1.0, COLOR=1
        data_loop, tab=current_tab
        IF PTR_VALID(series.geom) THEN PTR_FREE, series.geom
        IF TOTAL(gui.out.out_format) GE 1 THEN $
          XYOUTS, 0.3, 0.0, '!3click and hold in this window to abort', $
          /NORMAL, CHARSIZE=1.0, COLOR=0
        create_diag_table, tab=current_tab, /refresh
        message = DIALOG_MESSAGE('finished diagnosing',/INFORMATION,$
          DIALOG_PARENT=gui.window.main)
      ENDIF ELSE printerror, 'cannot run data loop without parameter file!'
    END

    gui.button.load : BEGIN                           ; click on load button
      load_form
      create_diag_table, /refresh
    END

    gui.button.save_as : save_form, /save_as          ; click on save as button

    gui.button.save : save_form, /last_form           ; click on save button

    gui.button.clear : clear_form                     ; click on clear form button

    gui.button.cdiags : BEGIN                         ; click on custom diags button
      IF NOT WIDGET_INFO(gui.window.cdiags,/VALID_ID) THEN BEGIN
        WIDGET_CONTROL, gui.button.cdiags, SET_VALUE='hide custom diags'
        custom_diag_list
      ENDIF ELSE BEGIN
        WIDGET_CONTROL, gui.window.cdiags, /DESTROY
        WIDGET_CONTROL, gui.button.cdiags, SET_VALUE='show custom diags'
      ENDELSE
    END

                                                ; click on recent H5 button
    gui.button.recent_h5 : IF gui.misc.recent_h5 NE '' THEN BEGIN
      recent_h5_files = STRSPLIT(gui.misc.recent_h5,/EXTRACT)
      window_sav = !D.WINDOW
      FOR j = 0, N_ELEMENTS(recent_h5_files) - 1 DO $
        dummy = H5_BROWSER(recent_h5_files[j])
      WSET, window_sav
    ENDIF

    gui.button.h5_files : BEGIN                 ; click on H5 files button
      window_sav = !D.WINDOW
      show_h5_files
      WSET, window_sav
    END

                                                ; click on recent ps button
    gui.button.recent_ps : IF gui.misc.recent_ps NE '' THEN BEGIN
      recent_ps_files = STRSPLIT(gui.misc.recent_ps,/EXTRACT)
      FOR j = 0, (SIZE(recent_ps_files))[1] - 1 DO $
        SPAWN, cfg.ps_viewer + recent_ps_files[j] + $
        ' &' ; /NOSHELL would block the gui !
    ENDIF

    gui.button.ps_files : show_ps_files         ; click on ps files button

                                                ; click on exit button
    gui.button.exit : WIDGET_CONTROL, gui.window.main, /DESTROY

    ; --- events of tabs
    gui.tab.base : BEGIN
      tab = WIDGET_INFO(gui.tab.base,/TAB_CURRENT)
      create_diag_table, tab=tab, /refresh
      CASE (*gui.tab.arr)[tab].name OF
        'scan' : BEGIN
;          WIDGET_CONTROL, gui.text.sortmeth, /EDITABLE
          WIDGET_CONTROL, gui.text.data_path, GET_VALUE=data_path
          IF STRMID(data_path,STRLEN(data_path)-1,1) NE '/' THEN $
            data_path += '/'
          IF FILE_TEST(data_path+'scan.log') THEN BEGIN
            series.par_updated = 0
            WIDGET_CONTROL, gui.text.runs, SET_VALUE='0001'
            apply_form
;            read_scanlog
            reload_diag_var_names, tab
            create_diag_table, tab=tab, /rebuild
          ENDIF ELSE printerror, 'no scan.log file found'
        END
        ELSE :
      ENDCASE
    END

    ; --- events of the diag tables ---
    (*gui.tab.arr)[WIDGET_INFO(gui.tab.base,/TAB_CURRENT)].table_id : BEGIN
      tab = WIDGET_INFO(gui.tab.base,/TAB_CURRENT)
      myid = (*gui.tab.arr)[tab].table_id
      diag_tab = STRARR((*gui.tab.arr)[tab].n_params+3,(*gui.tab.arr)[tab].n_diags)

      IF ev.type EQ 4 THEN BEGIN                ; selection change
        ; --- toggle all diags on/off ---
        ; Note: Clicking on the on/off header cell a second time
        ; results in a double -1,-1,-1,-1 event, one instance of
        ; which can be avoided doing the SET_TABLE_SELECT below.
        ; However, after toggling all diags, one needs to click
        ; one additional time in the header cell before it will
        ; accept the next all diag toggle. Currently, no remedy
        ; for this issue is known.
        IF (ev.sel_left EQ 2) AND (ev.sel_top EQ 0) AND (ev.sel_bottom GT 0) THEN BEGIN
          WIDGET_CONTROL, myid, GET_VALUE=diag_tab
          test_alloff = 1

          row = ev.sel_top
          WHILE test_alloff AND (row LE ev.sel_bottom) DO BEGIN
            diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,row)
            IF PTR_VALID(diag) THEN IF diag_tab[2,row] EQ 'X' THEN $
              test_alloff = 0
            row += 1
          ENDWHILE

          FOR row = ev.sel_top, ev.sel_bottom DO BEGIN
            diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,row)
            IF PTR_VALID(diag) THEN BEGIN
              diag_tab[2,row] = test_alloff ? 'X' : ' '
              (*diag).selected = test_alloff
            ENDIF
          ENDFOR

          WIDGET_CONTROL, myid, SET_VALUE=diag_tab, $
            SET_TABLE_SELECT=[-1,-1,-1,-1]
          WIDGET_CONTROL, myid, SET_TABLE_VIEW=viewport
        ENDIF

        IF (ev.sel_top EQ ev.sel_bottom) THEN BEGIN  ; no multiple selections

          viewport = WIDGET_INFO(myid,/TABLE_VIEW)

          ; --- help buttons ---
          IF ev.sel_left EQ 0 THEN BEGIN
            IF NOT WIDGET_INFO(gui.window.help,/VALID_ID) THEN BEGIN
              diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,ev.sel_top)
              display_help, diag
            ENDIF ELSE WIDGET_CONTROL, gui.window.help, /DESTROY
            WIDGET_CONTROL, myid, SET_TABLE_SELECT=[-1,-1,-1,-1]
            WIDGET_CONTROL, myid, SET_TABLE_VIEW=viewport
          ENDIF

          ; --- toggle diags on/off ---
          IF ev.sel_left EQ 2 THEN BEGIN
            WIDGET_CONTROL, myid, GET_VALUE=diag_tab
            diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,ev.sel_top)
            IF PTR_VALID(diag) THEN BEGIN
              IF diag_tab[2,ev.sel_top] EQ 'X' THEN BEGIN
                diag_tab[2,ev.sel_top] = ' '
                (*diag).selected = 0
              ENDIF ELSE BEGIN
                diag_tab[2,ev.sel_top] = 'X'
                (*diag).selected = 1
              ENDELSE
            ENDIF
            WIDGET_CONTROL, myid, SET_VALUE=diag_tab, $
              SET_TABLE_SELECT=[-1,-1,-1,-1]
            WIDGET_CONTROL, myid, SET_TABLE_VIEW=viewport
          ENDIF

          ; --- update diag table labels ---
          diag_tab[0,0] = ' '
          diag_tab[1,0] = ' '
          diag_tab[2,0] = ' '
          diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,ev.sel_top)     ; get selected diagnostic
          IF PTR_VALID(diag) THEN BEGIN
            vars = (*diag).table_entry
            IF vars[0,0] EQ '' THEN n_vars = 0 ELSE n_vars = N_ELEMENTS(vars[0,*])
            v = 0
            FOR i = 0, (*gui.tab.arr)[tab].n_params - 1 DO BEGIN
              IF v LT n_vars THEN BEGIN        ; write name to table header
                diag_tab[i+3,0] = vars[0,v]
                v = v + 1
              ENDIF ELSE diag_tab[i+3,0] = '--'
            ENDFOR
            WIDGET_CONTROL, myid, COLUMN_LABELS=diag_tab[*,0]
          ENDIF

         ; --- toggle switch parameters of diags ---
          IF (ev.sel_left GE 3) AND (ev.sel_left LT 3 + (*gui.tab.arr)[tab].n_params) $
            THEN BEGIN

            WIDGET_CONTROL, myid, GET_VALUE=diag_tab
            diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,ev.sel_top)
            IF(PTR_VALID(diag)) THEN BEGIN
	      n_vars = N_ELEMENTS((*diag).table_entry[0,*])
	      IF ((*diag).table_entry[0,0] NE '') AND (n_vars GT ev.sel_left-3) THEN BEGIN
	        IF (*diag).table_entry[1,ev.sel_left-3] EQ '1' THEN BEGIN
		  IF diag_tab[ev.sel_left,ev.sel_top] EQ '  off' THEN BEGIN
            	    diag_tab[ev.sel_left,ev.sel_top] = '  on'
		    (*diag).table_entry[3,ev.sel_left-3] = '1'
		  ENDIF ELSE BEGIN
            	    diag_tab[ev.sel_left,ev.sel_top] = '  off'
		    (*diag).table_entry[3,ev.sel_left-3] = '0'
		  ENDELSE
		ENDIF
              ENDIF
            ENDIF ; PTR_VALID ?
            WIDGET_CONTROL, myid, SET_VALUE=diag_tab, $
              SET_TABLE_SELECT=[-1,-1,-1,-1]
            WIDGET_CONTROL, myid, SET_TABLE_VIEW=viewport
          ENDIF ; change within table (horizontal borders)?
        ENDIF ; change within table (vertical borders)?
      ENDIF ; selection change event ?

      ; --- set regular parameters of diags ---
      IF ev.type EQ 0 THEN BEGIN       ; cell value change
        WIDGET_CONTROL, myid, GET_VALUE=diag_tab
        diag = get_diag_ptr((*gui.tab.arr)[tab].diag0,ev.y)
        IF PTR_VALID(diag) THEN BEGIN
	  n_vars = N_ELEMENTS((*diag).table_entry[0,*])
	  IF ((*diag).table_entry[0,0] NE '') AND (n_vars GT ev.x-3) THEN $
	    (*diag).table_entry[3,ev.x-3] = diag_tab[ev.x,ev.y]
        ENDIF ; PTR_VALID ?
      ENDIF ; cell value change event
    END

    ELSE:
  ENDCASE

END
