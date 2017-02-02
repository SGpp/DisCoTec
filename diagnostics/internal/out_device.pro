PRO set_output, diag, species, suffix=suffix, ps=ps, reset=reset, dat=dat, $
  coltable=coltable, multi=multi, xsize=xsize, ysize=ysize, charsize=charsize, $
  eps=eps, header=header, resolution=resolution, dname=dname, $
  commentlines=commentlines, append=append, write_frame=write_frame, mpgfile=mpgfile

; diag  : pointer to current diagnostic structure (obligatory argument)
; species : (optional) current species number to mark output files with species label
; suffix : (optional) additional suffix to file names
; ps    : keyword to switch on output to PS file
; reset : keyword that has to be used after all plotting events to close files
; dat   : set dat to arrays to be stored in a data file (ASCII or HDF5)
; coltable : (optional) keyword to select a certain color table; default: 41
; multi : (optional) keyword to set multi option (cmp. IDL plot command, default [0,1,1])
; xsize : (optional) keyword to set x size in cm; default: 17
; ysize : (optional) keyword to set y size in cm; default: 25
; charsize : (optional) keyword to set device char size; default: 1.5
; eps   : keyword to switch on EPS output -- only use for single page output !!
; !! the eps keyword needs to be set also when resetting at the end of the diag !!
; resolution : two element array with x and y resolution; default: [800,600]
; header : (optional) string or array of strings to be written as first line to data files
; commentlines : (optional) add comment lines (array for > 1 lines) to data file
; append : keyword to append data to data file (two blank lines will be written in between);
;  	   default: no append
; dname : (optional) string containing device name; default: !D.NAME
; write_frame: (optional) keyword providing a frame number for PNG output
; mpgfile : (optional) MPEG output requires a variable to store the
;           file ID, hence initialize a variable with 1 and provide
;           this as a parameter to set_output

  COMMON global_vars

  IF N_ELEMENTS(diag) EQ 0 THEN diag = ''
  IF PTR_VALID(diag) THEN name = (*diag).name ELSE name = diag
  spec_str = ''
  FOR j = 0, N_ELEMENTS(species) - 1 DO spec_str += spec[species[j]].name
  IF NOT KEYWORD_SET(suffix) THEN suffix=''
  IF N_ELEMENTS(coltable) LT 1 THEN coltable = 41
  IF NOT KEYWORD_SET(multi) THEN multi = [0,1,1]
  IF NOT KEYWORD_SET(xsize) THEN xsize = 17 ; A4 minus 2x2cm border
  IF NOT KEYWORD_SET(ysize) THEN ysize = 25
  IF NOT KEYWORD_SET(charsize) THEN charsize = 1.5
  IF N_ELEMENTS(resolution) NE 2 THEN resolution=[800,600]
  IF NOT KEYWORD_SET(dname) THEN dname = !D.NAME
  IF KEYWORD_SET(ps) AND KEYWORD_SET(eps) THEN eps = 0 ; cannot do both, switching to PS
  IF N_ELEMENTS(append) NE 1 THEN append = 0
  framenr_str = ''
  IF (N_ELEMENTS(write_frame) GT 0) AND (series.step_count GT 1) THEN BEGIN
     framenr = write_frame
     fmt = "(I0"+STRING((FIX(ALOG10(series.step_count))+1)>4,FORMAT="(I1)")+")"
     framenr_str = "_"+STRING(framenr,FORMAT=fmt)
  ENDIF ELSE framenr = 0

  IF N_ELEMENTS(mpgfile) EQ 0 THEN mpgfile = 0

  run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
  file_ext = spec_str + suffix + '_' + run_array[0]
  IF N_ELEMENTS(run_array) GT 1 THEN $
    file_ext += '_' + run_array[N_ELEMENTS(run_array)-1]

  IF (NOT FILE_TEST(gui.out.out_path,/DIRECTORY)) AND (KEYWORD_SET(ps) OR $
    KEYWORD_SET(eps) OR KEYWORD_SET(dat)) THEN $
    SPAWN, 'mkdir ' + gui.out.out_path

  kpar_ext = cfg.kpar_filter LT 0 ? '' : '_kpar' + rm0es(cfg.kpar_filter)

  filename = gui.out.out_path + name + file_ext + kpar_ext

  IF N_ELEMENTS(dat) GT 0 THEN BEGIN
    IF gui.out.out_format[0] THEN BEGIN ; ASCII output
      datfile = filename + '.dat'
      OPENW, data_lun, datfile, /GET_LUN, ERROR=err, append=append
      IF err NE 0 THEN BEGIN
        PRINT, 'unable to write data file'
        RETURN
      ENDIF

      IF KEYWORD_SET(append) THEN BEGIN
        PRINTF, data_lun, ''
        PRINTF, data_lun, ''
      ENDIF ELSE BEGIN
        PRINT, 'writing ASCII data file ', datfile
        PRINTF, data_lun, '# CREATED BY GENE DIAGNOSTICS TOOL on '+SYSTIME()
        run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
        line = '# diag: '+name+', run(s): '+run_array[0]
        IF N_ELEMENTS(run_array) GT 1 THEN $
           line += '-' + run_array[N_ELEMENTS(run_array)-1]
        line += ', time window: '+rm0es(gui.out.start_t)+'-'+rm0es(gui.out.end_t)
        line += '('+rm0es(series.step_count)+' steps)'
        IF spec_str NE '' THEN line +=', species: '+spec_str
        PRINTF, data_lun, line
        line = '# original data path: '+gui.out.data_path
        PRINTF, data_lun, line
        line = '# normalization: '+rm_idl_fmt(series.Lref_str)
        IF (series.Lref EQ 1.0) THEN line += '='+rm0es(par.Lref,prec=4)+'m'
        IF add_norm THEN BEGIN
           line += ', '+rm_idl_fmt(series.mref_str)
           IF (series.mref EQ 1.0) THEN line += '='+rm0es(par.mref,prec=4)+'kg'
           line += ', '+rm_idl_fmt(series.Tref_str)
           IF (series.Tref EQ 1.0) THEN line += '='+rm0es(par.Tref,prec=4)+'eV'
           IF (series.nref EQ 1.0) THEN line += ', n_ref='+rm0es(par.nref,prec=4)+'/m^3'
           IF (par.omegatorref NE 0.0) THEN line += ', omega_tor,ref='+rm0es(par.omegatorref,prec=4)+'/(rad/s)'
        ENDIF
        PRINTF, data_lun, line
        IF (series.doppler_corr_mom NE 0.0) THEN PRINTF, data_lun, $
           '# mom data shifted by omega_tor '+rm0es(series.doppler_corr_mom)+' cref/Lref'
        PRINTF, data_lun, '#'
      ENDELSE

      IF KEYWORD_SET(commentlines) THEN BEGIN
        FOR nlines = 0, N_ELEMENTS(commentlines) - 1 DO $
          PRINTF, data_lun, '# ' + commentlines[nlines]
      ENDIF

      n_col_str = rm0es(N_ELEMENTS(dat[0,*]))
      IF KEYWORD_SET(header) THEN BEGIN
        header='  '+header
        header[0]='#'+header[0]
        PRINTF, data_lun, header, FORMAT='('+n_col_str+'A-12)'
      ENDIF

      IF N_ELEMENTS(dat) EQ 1 THEN $
        PRINTF, data_lun, dat, FORMAT='('+n_col_str+'E12.4)' $
      ELSE PRINTF, data_lun, TRANSPOSE(dat), FORMAT='('+n_col_str+'E12.4)'
      FREE_LUN, data_lun
  ENDIF
    IF gui.out.out_format[1] THEN BEGIN ; HDF5 output
      h5file = filename + '.h5'

      IF NOT append THEN BEGIN
        FILE_DELETE, h5file, /ALLOW_NONEXISTENT

        file_id = H5F_CREATE(h5file)
        PRINT, 'writing HDF5 data file ', h5file

        ; general information block
        run_array = STRSPLIT(series.run_labels,',',/EXTRACT)
        run_str = run_array[0]
        IF N_ELEMENTS(run_array) GT 1 THEN run_str += $
          '-' + run_array[N_ELEMENTS(run_array)-1]
        info_type = ['created by','creation time','diagnostic name','run(s)',$
          'time window','step count','species','normalization L_ref',$
          'Tref/eV', 'nref/m^-3', 'mref/kg', 'Lref/m', 'Bref/T','omegatorref/(rad/s)']
        info_content = ['GENE Diagnostics Tool',SYSTIME(),name,run_str,$
          rm0es(gui.out.start_t)+'-'+rm0es(gui.out.end_t),$
          rm0es(series.step_count),spec_str,rm_idl_fmt(series.Lref_str),$
          rm0es(par.Tref), rm0es(par.nref), rm0es(par.mref), rm0es(par.Lref), $
          rm0es(par.Bref),rm0es(par.omegatorref)]
        IF (series.doppler_corr_mom NE 0.0) THEN info_content = $
           [[info_content],'mom data shifted by omega_tor '+$
            rm0es(series.doppler_corr_mom)+' cref/Lref']
        info_type = str_h5proof(info_type)
        info_content = str_h5proof(info_content)

        info_group_id = H5G_CREATE(file_id,'general information')
        dataspace_id = H5S_CREATE_SCALAR()
        FOR j = 0, N_ELEMENTS(info_type) - 1 DO BEGIN
          datatype_id = H5T_IDL_CREATE(info_content[j])
          info_content_current = info_content[j]
          dataset_id = $
            H5D_CREATE(info_group_id,info_type[j],datatype_id,dataspace_id)
          H5D_WRITE, dataset_id, info_content_current
          H5D_CLOSE, dataset_id
          H5T_CLOSE, datatype_id
        ENDFOR
        H5S_CLOSE, dataspace_id
        H5G_CLOSE, info_group_id

        ; comment block
        IF N_ELEMENTS(commentlines) GT 0 THEN BEGIN
          datatype_id = H5T_IDL_CREATE(commentlines)
          dataspace_id = H5S_CREATE_SIMPLE(N_ELEMENTS(commentlines))
          dataset_id = H5D_CREATE(file_id,'comments',datatype_id,dataspace_id)
          H5D_WRITE, dataset_id, commentlines
          H5D_CLOSE, dataset_id
          H5S_CLOSE, dataspace_id
          H5T_CLOSE, datatype_id
        ENDIF
      ENDIF ELSE file_id = H5F_OPEN(h5file,/WRITE)

      ; data block, using header info for naming
      IF NOT KEYWORD_SET(header) THEN BEGIN
        IF ((NOT append) OR (framenr GT 0)) THEN BEGIN
          datatype_id = H5T_IDL_CREATE(dat)
          dataspace_id = H5S_CREATE_SIMPLE(SIZE(dat,/DIMENSIONS))
          dataset_id = H5D_CREATE(file_id,'data'+framenr_str,$
                                  datatype_id,dataspace_id)
          H5D_WRITE, dataset_id, dat
          H5D_CLOSE, dataset_id
          H5S_CLOSE, dataspace_id
          H5T_CLOSE, datatype_id
        ENDIF
      ENDIF ELSE BEGIN
        IF N_ELEMENTS(header) EQ 1 THEN BEGIN
          datatype_id = H5T_IDL_CREATE(dat)
          dataspace_id = H5S_CREATE_SIMPLE(SIZE(dat,/DIMENSIONS))
          dataset_id = $
            H5D_CREATE(file_id,str_h5proof(header+framenr_str),$
                       datatype_id,dataspace_id)
          H5D_WRITE, dataset_id, dat
          H5D_CLOSE, dataset_id
          H5S_CLOSE, dataspace_id
          H5T_CLOSE, datatype_id
        ENDIF ELSE BEGIN
          datatype_id = H5T_IDL_CREATE(dat[*,0])
          dataspace_id = H5S_CREATE_SIMPLE(SIZE(dat[*,0],/DIMENSIONS))
          FOR j = 0, N_ELEMENTS(header) - 1 DO BEGIN
            data_name = j EQ 0 ? str_h5proof(header[j]+framenr_str) : $
              str_h5proof(header[j]+framenr_str) + ',' + $
                        str_h5proof(header[0]+framenr_str)
            dataset_id = $
              H5D_CREATE(file_id,data_name,datatype_id,dataspace_id)
            H5D_WRITE, dataset_id, dat[*,j]
            H5D_CLOSE, dataset_id
          ENDFOR
          H5S_CLOSE, dataspace_id
          H5T_CLOSE, datatype_id
        ENDELSE
      ENDELSE

      H5F_CLOSE, file_id

      gui.misc.recent_h5 = gui.misc.recent_h5 + ' ' + h5file
    ENDIF

    RETURN
  ENDIF

  !P.MULTI = multi
  !P.CHARSIZE = charsize

  active = WHERE(gui.out.out_format[2:*] EQ 1)
  FOR j = 0, N_ELEMENTS(active) - 1 DO $
    IF active[j] GE 0 THEN active[j] += 2

  IF NOT KEYWORD_SET(reset) AND N_ELEMENTS(write_frame) EQ 0 THEN BEGIN
    IF TOTAL(active) LE 0 THEN BEGIN
      IF coltable EQ 99 THEN set_redyellow ELSE $
      LOADCT, coltable, FILE='internal/colortable.tbl'
      RETURN
    ENDIF
    ps_nr = WHERE(gui.out.out_format_str EQ 'postscript')
    ps_device = (gui.out.out_format[ps_nr] EQ 1) AND $
                (mpgfile NE 1)

    IF ps_device AND (TOTAL(gui.out.out_format[2:*]) GT 1) THEN BEGIN
      dat_status = gui.out.out_format[0:1]
      gui.out.out_format = 0
      gui.out.out_format[0:1] = dat_status
      gui.out.out_format[ps_nr] = 1
      printerror, 'postscript output incompatible with other ' + $
        'graphic devices, thus they will be ignored'
    ENDIF
    IF ps_device THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, /COLOR, XSIZE=xsize, YSIZE=ysize, XOFFSET=2, YOFFSET=2,$
      ENCAPSULATED=KEYWORD_SET(eps)
;      LOADCT, coltable, FILE='internal/colortable.tbl'
      filename += '.'
      IF KEYWORD_SET(eps) THEN filename += 'e'
      filename += 'ps'
      PRINT, 'writing file ', filename
      DEVICE, FILENAME=filename
    ENDIF ELSE BEGIN
      SET_PLOT, 'Z'
      ERASE
      DEVICE, SET_RESOLUTION=resolution, SET_PIXEL_DEPTH=24, /DECOMPOSED

      FOR j = 0, N_ELEMENTS(active) - 1 DO BEGIN
         out_format = gui.out.out_format_str[active[j]]
         IF (active[j] EQ 2) THEN out_format = 'mpg'

         IF (out_format EQ 'mpg') THEN $
            mpgfile = MPEG_OPEN(resolution,QUALITY=100)
      ENDFOR

;      LOADCT, coltable, FILE='internal/colortable.tbl'
    ENDELSE
 
    IF coltable EQ 99 THEN set_redyellow ELSE $
    LOADCT, coltable, FILE='internal/colortable.tbl'
  ENDIF ELSE IF KEYWORD_SET(reset) THEN BEGIN ; reset
    CASE dname OF
      'PS' : BEGIN
        filename += '.'
        IF KEYWORD_SET(eps) THEN filename += 'e'
        filename += 'ps'
        DEVICE, /CLOSE

        OPENU, ps_lun, filename, /GET_LUN
        line = ''
        file_pos = 0L
        title_replaced = 0
        WHILE (NOT EOF(ps_lun)) AND (title_replaced EQ 0) DO BEGIN
          POINT_LUN, - ps_lun, file_pos
          READF, ps_lun, line
          IF STRMID(line,2,31) EQ 'Title: Graphics produced by IDL' THEN BEGIN
            old_length = STRLEN(line)
            line = '%%Title: ' + name + file_ext
            IF STRLEN(line) GT 32 THEN line = STRMID(line,0,32) + ' +'
            new_length = STRLEN(line)
            FOR i = 1, old_length - new_length DO line += ' '
            POINT_LUN, ps_lun, file_pos
            PRINTF, ps_lun, line
            title_replaced = 1
          ENDIF
        ENDWHILE
        FREE_LUN, ps_lun
        gui.misc.recent_ps = gui.misc.recent_ps + ' ' + filename
      END
      'Z' : BEGIN
         FOR j = 0, N_ELEMENTS(active) - 1 DO BEGIN
            out_format = gui.out.out_format_str[active[j]]
            IF (active[j] EQ 2) THEN out_format = 'mpg'

            IF (framenr EQ 0) THEN BEGIN ;don't write last entry twice
               write_to_z_device, filename, resolution, framenr, $
                    framenr_str, out_format, mpgfile
            ENDIF
            
            IF (out_format NE 'mpg') THEN $
               DEVICE,/CLOSE $
            ELSE BEGIN
               myfilename = filename + '.'+out_format
               print, 'wrote file: '+myfilename
               MPEG_SAVE, mpgfile, FILENAME=myfilename
               MPEG_CLOSE, mpgfile
               ERASE
            ENDELSE
         ENDFOR
      END
      ELSE :
    ENDCASE

    USERSYM, [-1,0,1,0], [0,1,0,-1], /FILL ; create plotting symbol
    IF coltable EQ 99 THEN set_redyellow ELSE $
    LOADCT, coltable, FILE='internal/colortable.tbl'

    IF FILE_TEST('win32') EQ 0 THEN BEGIN
      SET_PLOT, 'X'
      DEVICE, DECOMPOSED=0, TRUE_COLOR=24
    ENDIF ELSE BEGIN
      SET_PLOT, 'WIN'
      DEVICE, DECOMPOSED=0
    ENDELSE
  ENDIF ELSE IF (framenr GT 0) THEN BEGIN
     CASE dname OF
      'Z' : BEGIN
        FOR j = 0, N_ELEMENTS(active) - 1 DO BEGIN
           out_format = gui.out.out_format_str[active[j]]
           IF (active[j] EQ 2) THEN out_format = 'mpg'
              
           write_to_z_device, filename, resolution, framenr, $
                framenr_str, out_format, mpgfile
        ENDFOR
      END
      ELSE :
    ENDCASE
  ENDIF

END


PRO write_to_z_device, filename, resolution, framenr, framenr_str, out_format,$
                       mpgfile

  frame = INTARR(resolution[0],resolution[1],/NOZERO)
  IF (out_format NE 'mpg') THEN BEGIN
     frame = TVRD()
     myfilename = filename + framenr_str+'.'+out_format 
  ENDIF ELSE BEGIN
     frame = REVERSE(TVRD(),2)
     myfilename = filename + '.'+out_format
  ENDELSE

  TVLCT, red, green, blue, /GET
  image = BYTARR(3,resolution[0],resolution[1])

  image[0,*,*] = red[frame[*,*]]
  image[1,*,*] = green[frame[*,*]]
  image[2,*,*] = blue[frame[*,*]]

  IF (out_format NE 'mpg') THEN BEGIN
     WRITE_IMAGE, myfilename, out_format, image, /TRUE
     PRINT, 'wrote file: ' + myfilename
  ENDIF ELSE BEGIN
     print, 'writing frame '+rm0es(framenr)+' to '+myfilename
     MPEG_PUT, mpgfile, FRAME=framenr, IMAGE=image[*,*,*]
  ENDELSE
END



PRO set_redyellow
   steps = 255
   redVector = REPLICATE(255, steps)
   blueVector = REPLICATE(0, steps)
;The green vector (according to the algorithm) is:

   scaleFactor = FINDGEN(steps) / (steps - 1)
   beginNum = 255
   endNum = 0
   greenVector = beginNum + (endNum - beginNum) * scaleFactor

;Alright, now load these color vectors, and there you have it, a color table smoothly progressing from yellow to red!

   TVLCT, redVector, greenVector, blueVector
END
