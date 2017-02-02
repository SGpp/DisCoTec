FUNCTION read_iterdb_1d_data, file, var_name_in, nr_var=nr_var
  ; returns 1D data array for variable var_name_in 
  ; from 1D iterdb file
  ; nr_var : variable containing the lengths of the arrays
  ;          default: nj

  IF NOT KEYWORD_SET(nr_var) THEN nr_var='nj'

  IF NOT FILE_TEST(file) THEN BEGIN
      print, "Error: file "+file+" does not exist!"
      RETURN, -1
  ENDIF

  OPENR,lun,file,err=err,/get_lun
  found = 0

  ;check for nj: the size of the vectors printed in the file
  WHILE ((NOT (EOF(lun))) AND (NOT found)) DO BEGIN
     ; read header
     line = ''
     READF, lun, line
     IF (STRPOS(line, nr_var) GT 0) THEN BEGIN
        found = 1
        nj = 0
        READF, lun, nj
     ENDIF
  ENDWHILE

  IF (NOT found) THEN BEGIN
     print, "couldn't read vector size nj in "+file
     return, -1
  ENDIF

  ;rewind
  POINT_LUN, lun, 0
  found = 0

  WHILE ((NOT (EOF(lun))) AND (NOT found)) DO BEGIN
     READF, lun, line
     IF (STRPOS(line, var_name_in) GT 0) THEN BEGIN
        found = 1
        data = FLTARR(nj,/NOZERO)
        READF, lun, data
     ENDIF
  ENDWHILE

  FREE_LUN, lun

  IF (found) THEN RETURN, data $
  ELSE RETURN, -1

END

;########################################################################

FUNCTION read_iterdb_2d_data, file, var_name_in
  ; returns 2D data array for variable var_name_in 
  ; from 2D iterdb file


  IF NOT FILE_TEST(file) THEN BEGIN
      print, "Error: file "+file+" does not exist!"
      RETURN, -1
  ENDIF

  OPENR,lun,file,err=err,/get_lun
  found = 0

  WHILE ((NOT (EOF(lun))) AND (NOT found)) DO BEGIN
     ; read header
     line = ''
     FOR i = 1, 5 DO BEGIN
         READF, lun, line
;         print, line
     ENDFOR

     line = ''
     READF, lun, line
     var_name = STRTRIM(STRMID(line,1,15))
     
     READF, lun, line
     READF, lun, line 
     nx = FIX(STRTRIM(STRMID(line,1,15)))
     READF, lun, line 
     ny = FIX(STRTRIM(STRMID(line,1,15)))
     xaxis = FLTARR(nx,/NOZERO)
     READF, lun, xaxis
     time = 0.0
     READF, lun, time
     ydata = FLTARR(nx,/NOZERO)
     READF, lun, ydata
     IF (var_name EQ var_name_in) THEN found = 1

     ;comment lines
     FOR i=1,3 DO BEGIN
         READF, lun, line
;         print, line
     ENDFOR
  ENDWHILE

  FREE_LUN, lun

  IF (found) THEN BEGIN
     data = FLTARR(nx,2,/NOZERO)
     data[*,0] = xaxis
     data[*,1] = ydata
     RETURN, data
  ENDIF ELSE RETURN, -1

END

;########################################################################
  
