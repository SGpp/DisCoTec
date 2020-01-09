;##########################################################################
;# this library contains almost all scan related internal functions       #
;##########################################################################

FUNCTION scan_renorm, var

  CASE var OF
    0 : RETURN, 1.0D
    1 : RETURN, 1.0D
    2 : RETURN, 1.0D
    3 : RETURN, 1.0D
    4 : RETURN, 1.0D
    ELSE : printerror, 'error in scan_renorm'
  ENDCASE

END

;#############################################################################

FUNCTION get_scan_time, run, get_n_steps=get_n_steps, $
  first=first, last=last, step_time=step_time, coarse=coarse

  RETURN, get_mom_time(run,get_n_steps=get_n_steps,$
    first=first, last=last, step_time=step_time, coarse=coarse)

END

;#############################################################################

FUNCTION compvalues, meth, temptemp, comptemp

  arr1 = temptemp
  arr2 = comptemp
;chosen will be vglval=MINIMUM
  CASE meth OF
    ; smallest distance in the complex plane
    1 : BEGIN
          arr1_real = FLOAT(arr1)
          arr1_imag = IMAGINARY(arr1)
          inf_inds = WHERE(NOT FINITE(arr1_real),count)
          IF count GT 0 THEN arr1_real[inf_inds] = 0
          inf_inds = WHERE(NOT FINITE(arr1_imag),count)
          IF count GT 0 THEN arr1_imag[inf_inds] = 0
          arr1 = COMPLEX(arr1_real,arr1_imag)

          arr2_real = FLOAT(arr2)
          arr2_imag = IMAGINARY(arr2)
          inf_inds = WHERE(NOT FINITE(arr2_real),count)
          IF count GT 0 THEN arr2_real[inf_inds] = 0
          inf_inds = WHERE(NOT FINITE(arr2_imag),count)
          IF count GT 0 THEN arr2_imag[inf_inds] = 0
          arr2 = COMPLEX(arr2_real,arr2_imag)

          vglval = TOTAL(ABS(arr1-arr2))
        END
    2 : BEGIN
          vglval=N_ELEMENTS(arr2)                  
          FOR j=1,N_ELEMENTS(arr2)-1 DO BEGIN ;sort imaginary part by size
            IF FINITE(IMAGINARY(arr2(j))) EQ 1 AND FINITE(IMAGINARY(arr2(j-1))) THEN BEGIN
              IF FINITE(IMAGINARY(arr1(j))) EQ 1 AND FINITE(IMAGINARY(arr1(j-1))) THEN BEGIN
                IF IMAGINARY(arr2(j)) GT IMAGINARY(arr2(j-1)) THEN BEGIN
                  IF IMAGINARY(arr1(j)) GT IMAGINARY(arr1(j-1)) THEN vglval -= 1
                ENDIF ELSE BEGIN
                  IF IMAGINARY(arr1(j)) LE IMAGINARY(arr1(j-1)) THEN vglval -= 1
                ENDELSE
              ENDIF
            ENDIF
          ENDFOR
        END
    3 : BEGIN
                                ;hybrid
          FOR i=0,N_ELEMENTS(arr1)-1 DO BEGIN
              exc1=REAL_PART(arr1(i))
              exc2=IMAGINARY(arr1(i))
              IF FINITE(exc1) NE 1 THEN exc1=FLOAT(0.)
              IF FINITE(exc2) NE 1 THEN exc2=FLOAT(0.)
              arr1(i)=COMPLEX(exc1,exc2)
          ENDFOR
          FOR i=0,N_ELEMENTS(arr2)-1 DO BEGIN
              exc1=REAL_PART(arr2(i))
              exc2=IMAGINARY(arr2(i))
              IF FINITE(exc1) NE 1 THEN exc1=FLOAT(0.)
              IF FINITE(exc2) NE 1 THEN exc2=FLOAT(0.)
              arr2(i)=COMPLEX(exc1,exc2)
          ENDFOR
          vglval=TOTAL(ABS(arr1-arr2))    
          FOR j=1,N_ELEMENTS(arr2)-1 DO BEGIN ;sort imaginary part by size
              IF FINITE(IMAGINARY(arr2(j))) EQ 1 AND FINITE(IMAGINARY(arr2(j-1))) THEN BEGIN
                  IF FINITE(IMAGINARY(arr1(j))) EQ 1 AND FINITE(IMAGINARY(arr1(j-1))) THEN BEGIN
                      IF IMAGINARY(arr2(j)) GT IMAGINARY(arr2(j-1)) THEN BEGIN
                          IF IMAGINARY(arr1(j)) GT IMAGINARY(arr1(j-1)) THEN BEGIN
                              vglval=vglval*0.9
                          ENDIF
                      ENDIF ELSE BEGIN
                          IF IMAGINARY(arr1(j)) LE IMAGINARY(arr1(j-1)) THEN BEGIN
                              vglval=vglval*0.9
                          ENDIF
                      ENDELSE
                  ENDIF
              ENDIF
          ENDFOR
      END
      4 : BEGIN
                                ;squared method1
          FOR i=0,N_ELEMENTS(arr1)-1 DO BEGIN
              exc1=REAL_PART(arr1(i))
              exc2=IMAGINARY(arr1(i))
              IF FINITE(exc1) NE 1 THEN exc1=FLOAT(0.)
              IF FINITE(exc2) NE 1 THEN exc2=FLOAT(0.)
              arr1(i)=COMPLEX(exc1,exc2)
          ENDFOR
          FOR i=0,N_ELEMENTS(arr2)-1 DO BEGIN
              exc1=REAL_PART(arr2(i))
              exc2=IMAGINARY(arr2(i))
              IF FINITE(exc1) NE 1 THEN exc1=FLOAT(0.)
              IF FINITE(exc2) NE 1 THEN exc2=FLOAT(0.)
              arr2(i)=COMPLEX(exc1,exc2)
          ENDFOR
          vglval=TOTAL((REAL_PART(arr1-arr2))^4+(IMAGINARY(arr1-arr2))^4)
      END
      ELSE :
  ENDCASE

  RETURN, vglval

END

;##########################################################################

FUNCTION permute, meth, valutemp, sorttemp, comptemp
; mintemp: vglvalue (to be minimal)
; mintemptemp=DCOMPLEXARR(arrsize)

  arrsize=N_ELEMENTS(sorttemp)
  mintemp=100000000LL
  valumintemptemp=valutemp
  sortmintemptemp=sorttemp
  IF(arrsize GT 1) THEN BEGIN
      FOR m=0,arrsize-1 DO BEGIN
          CASE m OF
              0: BEGIN
                  valurest=valutemp(1:arrsize-1)
                  sortrest=sorttemp(1:arrsize-1)
                  restarr=permute(meth,valurest,sortrest,comptemp(1:arrsize-1))
              END
              (arrsize-1): BEGIN
                  valurest=valutemp(0:m-1)
                  sortrest=sorttemp(0:m-1)
                  restarr=permute(meth,valurest,sortrest,comptemp(1:arrsize-1))
              END
              ELSE: BEGIN
                  valurest=[valutemp(0:m-1),valutemp(m+1:arrsize-1)]
                  sortrest=[sorttemp(0:m-1),sorttemp(m+1:arrsize-1)]
                  restarr=permute(meth,valurest,sortrest,comptemp(1:arrsize-1))
              END
          ENDCASE
          valutemptemp=[valutemp(m),valurest]
          sorttemptemp=[sorttemp(m),sortrest]
          vglvalue=compvalues(meth,valutemptemp,comptemp)
          IF mintemp GT vglvalue THEN BEGIN
              valumintemptemp=valutemptemp
              sortmintemptemp=sorttemptemp
              mintemp=vglvalue
          ENDIF
      ENDFOR
  ENDIF
  valutemp=valumintemptemp
  sorttemp=sortmintemptemp

  return,valumintemptemp

END

;##########################################################################


function sortevs,meth,valutemp,sorttemp,eigenfield,dsize,index
; poss: comparability (Integer)
; steps: distance of Neighbour to compare with in given dimension (Integer)
; valstosort: # finite values in valutemp
; dim: dimensionindex of dsize to compare with
; dimstep: distance in onedimensional subscript of the elements to compare
; compindex: first index (1-d) of the elements to compare
; k: Number of comparable elements of the comparepart 
; comptemp: comparearray (array(dsize(0)))
  steps=1
  poss=1
  valstosort=0
  FOR i=0,dsize(0)-1 DO BEGIN
      valstosort+=FINITE(valutemp(i))
  ENDFOR
  WHILE poss NE 0 DO BEGIN
      poss=0 
      FOR dim = 1,N_ELEMENTS(dsize)-1 DO BEGIN
          dimstep=1
          FOR i=0,(dim-1) DO BEGIN
              dimstep=dimstep*dsize(i)
          ENDFOR    
          compindex=index-dimstep*steps
          IF (compindex GE 0) THEN BEGIN
              poss=1
              k=0
              for i=compindex,((compindex+dsize(0))-1) DO BEGIN
                  k+=FINITE(eigenfield(i))
              ENDFOR
              IF k GE valstosort THEN BEGIN ;comparable Eigenvalues
                  comptemp=eigenfield(compindex : compindex+dsize(0)-1)
                  valutemp=permute(meth,valutemp,sorttemp,comptemp)
                  return, valutemp
              ENDIF
          ENDIF ELSE BEGIN
              
          ENDELSE
      ENDFOR
      steps=steps+1
  ENDWHILE
  comptemp=eigenfield(index-dsize(0):index-1)
  valutemp=permute(meth,valutemp,sorttemp,comptemp)
  return, valutemp
end

;#############################################################################

PRO read_scanlog

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN PTR_FREE, scanlog

  ftest_scan = FILE_TEST(gui.out.data_path+'scan.log')
  fline_scan = ftest_scan ? FILE_LINES(gui.out.data_path+'scan.log') : 0

  IF (fline_scan LT 3) AND cfg.scanlog_repair AND $
    FILE_TEST(gui.out.data_path+'parameters') THEN BEGIN

    IF !QUIET NE 1 THEN PRINT, 'attempting to repair scanlog...'
    sslink = '../tools/perl/scanscript'
    IF FILE_TEST("scanscript",/SYMLINK) THEN FILE_DELETE, "scanscript"
    FILE_LINK, sslink, "scanscript", /NOEXPAND_PATH
    SPAWN, './scanscript --mks --o=' + gui.out.data_path
    FILE_DELETE, "scanscript", /ALLOW_NONEXISTENT
    IF !QUIET NE 1 THEN PRINT, '...' + $
      (FILE_LINES(gui.out.data_path+'scan.log') GE 2 ? 'success' : 'failed')
  ENDIF

  OPENR, scanlog_lun, gui.out.data_path + 'scan.log', /GET_LUN, ERROR=err

  IF err NE 0 THEN BEGIN
    printerror, 'error: no scan.log file found in ' + gui.out.data_path
    RETURN
  ENDIF

  line=''
  READF,scanlog_lun,line
  dimarr=STRSPLIT(line,'|',COUNT=dimensions,/EXTRACT)
  dimensions=dimensions-1
  eigenarr=STRSPLIT(dimarr[dimensions],'/',COUNT=eigenvals,/EXTRACT)
  eigenvals = eigenvals -1
  dimarr[dimensions]=eigenarr[0]
  FOR indm=1,dimensions DO BEGIN
     splitnames=STRSPLIT(dimarr[indm],' ',/EXTRACT)
     IF (MAX(STRCMP(splitnames[0],['omn','omt','mass','charge','temp','dens'])) EQ 1) THEN BEGIN
        dimarr[indm]=splitnames[0]+' '+spec[splitnames[1]-1].name
     ENDIF ELSE BEGIN
        dimarr[indm]=splitnames[0]
     ENDELSE
  ENDFOR
  dsize=REPLICATE(1,dimensions+1) ;resolution of dimension (number of steps)&num of eigenvalues
  deltad=MAKE_ARRAY(dimensions,/FLOAT)    ;stepwidth of parameters
                                ;determine how many elements per file by going through file
  line=''
  READF,scanlog_lun,line
  el=FLOAT(STRTRIM(STRSPLIT(line,'|',/EXTRACT))) ;Run&Pars&EVs
  rd=DINDGEN(dimensions)
  FOR i=1,dimensions DO rd[i-1]=el[i] ;end of Interval
  WHILE NOT EOF(scanlog_lun) DO BEGIN
      line=''
      READF,scanlog_lun,line
      el=FLOAT(STRTRIM(STRSPLIT(line,'|',COUNT=linelen,/EXTRACT)))
      FOR i=1,dimensions DO BEGIN
          IF deltad[i-1] GT 0 THEN BEGIN
              IF el[i] GT rd[i-1] THEN BEGIN
                  dsize[i]=dsize[i]+1 ;size of Dimensions
                  deltad[i-1]=el[i]-rd[i-1] ;stepwidth of Dimensions
                  rd[i-1]=el[i] ;largest parameters of intervals
              ENDIF
          ENDIF ELSE BEGIN
              IF deltad[i-1] LT 0 THEN BEGIN
                  IF el[i] LT rd[i-1] THEN BEGIN
                      dsize[i]=dsize[i]+1 ;size of Dimensions
                      deltad[i-1]=el[i]-rd[i-1] ;stepwidth of Dimensions
                      rd[i-1]=el[i] ;largest parameters of intervals
                  ENDIF    
              ENDIF ELSE BEGIN
                  IF el[i] NE rd[i-1] THEN BEGIN
                      dsize[i]=dsize[i]+1 ;size of Dimensions
                      deltad[i-1]=el[i]-rd[i-1] ;stepwidth of Dimensions
                      rd[i-1]=el[i] ;largest parameters of intervals
                  ENDIF 
              ENDELSE 
          ENDELSE
      ENDFOR 
  ENDWHILE
  FREE_LUN,scanlog_lun
  dsize(0) = eigenvals
  sortfield=MAKE_ARRAY(dsize,/INTEGER) 
  eigenfield=MAKE_ARRAY(dsize,/COMPLEX) ;,VALUE=DCOMPLEX('NaN','NaN'))

  axxis=MAKE_ARRAY(dimensions, MAX(dsize), /FLOAT)
  axxindex=MAKE_ARRAY(dimensions, /INTEGER)

  sortmeth = gui.out.sortmeth

  OPENR,scanlog_lun,gui.out.data_path+'scan.log',/GET_LUN,ERROR=er
  IF err NE 0 THEN PRINT, 'not found'
  READF,scanlog_lun,line
  OPENR,scansort_lun,gui.out.data_path+'scansort.log',/GET_LUN,ERROR=err   
  i=0
  IF err EQ 0 THEN BEGIN
      WHILE (NOT EOF(scansort_lun)) AND $
        (NOT N_ELEMENTS(sortfield) LT i+dsize[0]-1) DO BEGIN
          line=''
          READF,scansort_lun,line
          estr=STRTRIM(STRSPLIT(line,' ',COUNT=sortentries,/EXTRACT))
          IF sortentries EQ dsize[0] THEN BEGIN
              sortfield[i:i+dsize[0]-1]=estr
              i=i+dsize[0]
          ENDIF
      ENDWHILE
      FREE_LUN,scansort_lun
      IF i NE N_ELEMENTS(sortfield) THEN BEGIN
          err=-1
          print, 'scansort.log is not valid'
      ENDIF
  ENDIF
  IF (sortmeth EQ 0) AND (err EQ 0) THEN BEGIN
      evnum=0
      WHILE NOT EOF(scanlog_lun) DO BEGIN
          READF,scanlog_lun,line
          estr=STRTRIM(STRSPLIT(line,'|',COUNT=sortentries,/EXTRACT))
          IF axxindex[0] EQ 0 THEN BEGIN
              FOR i=1, dimensions DO BEGIN
                  axxis[i-1, 0]=FLOAT(estr[i])
                  axxindex[*]=1
              ENDFOR       
          ENDIF ELSE BEGIN
              FOR i=1, dimensions DO BEGIN
                  IF (deltad[i-1])*(FLOAT(estr[i])-axxis[i-1, axxindex[i-1]-1]) $
                    GT 0 THEN BEGIN
                      axxis[i-1, axxindex[i-1]]=FLOAT(estr[i])
                      axxindex[i-1]+=1
                  ENDIF        
              ENDFOR             
          ENDELSE
          FOR i=0,dsize(0)-1 DO BEGIN
              IF sortfield(evnum) LT 0 THEN BEGIN
                  estrev=STRTRIM(STRSPLIT(estr[dimensions-sortfield(evnum)],' ',/EXTRACT))
                  IF STRCMP(estrev[0],'NaN',3) OR sortfield(evnum) EQ 0 THEN BEGIN
                      eigenfield[evnum]=COMPLEX('NaN','NaN')
                  ENDIF ELSE BEGIN
                      IF estrev[0] LT 0 THEN BEGIN
                          eigenfield[evnum]=COMPLEX(FLOAT(0.0),'NaN')
                      ENDIF ELSE BEGIN
                          eigenfield[evnum]=COMPLEX(FLOAT(estrev[0]),FLOAT(estrev[1]))
                      ENDELSE
                  ENDELSE
              ENDIF ELSE BEGIN
                  IF sortfield(evnum) EQ 0 THEN BEGIN
                      eigenfield[evnum]=COMPLEX(FLOAT(0.0),'NaN')
                  ENDIF ELSE BEGIN 
                      eigenfield[evnum]=COMPLEX('NaN','NaN')
                  ENDELSE
              ENDELSE
              evnum+=1
          ENDFOR
      ENDWHILE
                                ;read scansort.log and fill sortfield
  ENDIF ELSE BEGIN
                                ;read eigenfield and sort eigenvalues and print scansort.log
                                ;read EV's
      j=0
      WHILE NOT EOF(scanlog_lun) DO BEGIN
          line=''
          READF,scanlog_lun,line
          estr=STRTRIM(STRSPLIT(line,'|',COUNT=entries,/EXTRACT)) ;Run&Pars&EVs
          IF axxindex[0] EQ 0 THEN BEGIN
              FOR i=1, dimensions DO BEGIN
                  axxis[i-1, 0]=FLOAT(estr[i])
                  axxindex[*]=1
              ENDFOR       
          ENDIF ELSE BEGIN
              FOR i=1, dimensions DO BEGIN
                  IF (deltad[i-1])*(FLOAT(estr[i])-axxis[i-1, axxindex[i-1]-1]) $
                    GT 0 THEN BEGIN
                      axxis[i-1, axxindex[i-1]]=FLOAT(estr[i])
                      axxindex[i-1]+=1
                  ENDIF
              ENDFOR             
          ENDELSE
                                ;read one line of eigenvalues
          entries=entries-dimensions-1
                                ;initialize sortfield
          IF sortmeth EQ 0 THEN entries=dsize(0)
          FOR i=1,entries DO BEGIN
              estrev=STRTRIM(STRSPLIT(estr[dimensions+i],' ',/EXTRACT))
              IF i LE dsize(0) THEN BEGIN
                  IF STRCMP(estrev[0],'NaN',3) THEN BEGIN
                      eigenfield[j]=COMPLEX('NaN','NaN')
                      sortfield[j]=1
                  ENDIF ELSE BEGIN
                      IF estrev[0] LT 0.001 THEN BEGIN
                          eigenfield[j]=COMPLEX(FLOAT(0.0),'NaN')
                          sortfield[j]=0
                      ENDIF ELSE BEGIN                        
                          eigenfield[j]=COMPLEX(FLOAT(estrev[0]),FLOAT(estrev[1]))
                          sortfield[j]=-i
                      ENDELSE
                  ENDELSE
                  j=j+1
              ENDIF ELSE BEGIN
                  IF NOT STRCMP(estrev[0],'NaN',3) THEN BEGIN
                      IF (sortmeth GT 0) && (estrev[0] GT MIN(REAL_PART(eigenfield[j-dsize(0):j-1]),minsubs,/NAN)) THEN BEGIN
                          eigenfield[j-dsize(0)+minsubs]=COMPLEX(FLOAT(estrev[0]),FLOAT(estrev[1]))
                          sortfield[j-dsize(0)+minsubs]=-i
                      ENDIF                        
                  ENDIF
              ENDELSE
          ENDFOR   
          IF sortmeth GT 0 THEN BEGIN
              IF dsize(0) GE 2 THEN BEGIN
                                ;order set of eigenvalues by realpartsize
                  valutemp1=eigenfield[j-dsize(0):j-1]
                  sorttemp1=sortfield[j-dsize(0):j-1]
                  valutemp2=MAKE_ARRAY(dsize(0),/COMPLEX,VALUE=COMPLEX('NaN','NaN'))
                  sorttemp2=MAKE_ARRAY(dsize(0),/INTEGER,VALUE=1) 
                  i=0
                  WHILE MIN(sorttemp1) LT 1 DO BEGIN
                      hello=MAX(REAL_PART(valutemp1),maxsubscript,/NAN)
                      valutemp2(i)=valutemp1(maxsubscript)
                      valutemp1(maxsubscript)=COMPLEX(-1.,0.)
                      sorttemp2(i)=sorttemp1(maxsubscript)
                      sorttemp1(maxsubscript)=1
                      i+=1
                  ENDWHILE
                  
                                ;order eigenvalues by neighbours
                  k=0
                  FOR i=(j-dsize(0)),j-1 DO BEGIN
                      k=k+FINITE(eigenfield(i))
                  ENDFOR
                  IF (j GT dsize(0)) THEN BEGIN    
                      valutemp2=sortevs(sortmeth,valutemp2,sorttemp2,eigenfield,dsize,j-dsize(0))
                  ENDIF ELSE BEGIN ;else there's nothing to compare with
                  ENDELSE
                  FOR i=1,dsize(0) DO BEGIN
                      eigenfield[j-dsize(0)-1+i]=valutemp2[i-1]
                      sortfield[j-dsize(0)-1+i] =sorttemp2[i-1]
                  ENDFOR
              ENDIF
          ENDIF ELSE BEGIN
                                ; don't sort anything
          ENDELSE
      ENDWHILE
      OPENW,scansort_lun,gui.out.data_path+'scansort.log',/GET_LUN,ERROR=err 
      IF err NE 0 THEN BEGIN
          
          print,'not possible to write in file'
      ENDIF ELSE BEGIN
          print, 'written sortdata into scansort.log; feel free to manipulate'
          PRINTF,scansort_lun,sortfield
          FREE_LUN,scansort_lun    
      ENDELSE
  ENDELSE
  FREE_LUN,scanlog_lun
  ;return data
   scanlog = PTR_NEW({$
    dimarr     : dimarr,$     ; dimarr: names of the scanned parameters; starting with index 1
    dsize      : dsize,$      ; dsize[0]: number of eigenvalues,dsize[1]: number of parameters in every dimension
    deltad     : deltad,$     ; delta: size of parameter steps 
    rd         : rd,$         ; rd: largest parameter of every dimension
    n_scans    : FIX(el[0]),$ ; n_scans: number of scans/files
    eigenfield : eigenfield,$ ; eigenfield: array, containing the eigenvalues
    sortfield  : sortfield,$  ; sortfield: array, containing the sortinformation 0 means (0.0,NaN) (negative eigenvalue found, 1 means (NaN,NaN) (no eigenvalue found)
    axxis      : axxis})

END

;##########################################################################

PRO scan_loop, diag0

  COMMON global_vars

;  IF NOT FILE_TEST(gui.out.data_path+'scan.log') THEN BEGIN
;    printerror, 'skipping scan loop due to lack of scan.log file'
;    RETURN
;  ENDIF  

  IF NOT PTR_VALID(scanlog) THEN BEGIN
    printerror, 'No scan log structure loaded. Did you define a sort method?'
    RETURN
  ENDIF
  
  IF series.Lref NE 1.0 OR series.mref NE 1.0 OR series.Tref NE 1.0 THEN BEGIN
    series.Lref = set_Lref(0) ; reset to L_perp
    series.mref = 1.0D
    series.Tref = 1.0D
    series.Qref = 1.0D
    series.Bref = 1.0D    
    WIDGET_CONTROL, gui.droplist.Lref, SET_DROPLIST_SELECT=0
    WIDGET_CONTROL, gui.droplist.Tref, SET_DROPLIST_SELECT=0   
    WIDGET_CONTROL, gui.droplist.mref, SET_DROPLIST_SELECT=0
    print, 'INFO: scan loop does not support renormalization yet!'
  ENDIF

  steps = 0  
  IF series.request_nrg OR (TOTAL(series.request_mom) GT 0) THEN BEGIN  
   oldseries = series
  
   FOR run=1, (*scanlog).n_scans DO BEGIN ;do this for each file of the scan!
    series.run_labels = STRING(run, FORMAT='(I04)')
    file_par = get_file_path('parameters',/set_fm)
    IF NOT file_par.exist THEN BEGIN
       PRINTERROR, 'error in scan loop: missing parameter files'
       RETURN
    ENDIF
    
    read_par
    read_geometry
    set_series_lengths 
    IF series.request_nrg THEN read_nrg
    
    IF (TOTAL(series.request_mom) GT 0) THEN BEGIN ; read mom/field data    
     IF par.ntimesteps GT 0 THEN BEGIN ;initial value solver
      steps += 1
;      print, 'reading last time step'
      last_time = get_mom_time(0,/last,/coarse)

      gui.out.start_t = last_time * 0.99999
      gui.out.end_t =   last_time * 1.00001

      mom_loop, diag0
     ENDIF ELSE BEGIN ;eigen value solver
   
;    print, 'reading negative time steps'
       gui.out.start_t = -100.0
       ; distinguish old (neg. times) runs from newer (with comp_type)
       gui.out.end_t = par.comp_type EQ 'EV' ? last_time * 1.000001 : -0.5
       steps += 1
       mom_loop, diag0
     ENDELSE  
    ENDIF ; reading mom/field data
    
    PTR_FREE, series.geom     
   ENDFOR  
   series = oldseries  
  ENDIF ; request_nrg and request_mom
  
  series.step_count = steps>1
  
END

;######################################################################

FUNCTION finddate,diag,run_label,date2plot,title

  COMMON global_vars
  
  retval=0
  IF (date2plot EQ 0) THEN BEGIN
      OPENR,scandate_lun,gui.out.data_path+'eigenvalues_'+$
            run_label,/GET_LUN,ERROR=err
      IF err NE 0 THEN BEGIN
          print, 'could not open file'
          retval=FLOAT('NaN')
      ENDIF ELSE BEGIN
          line=' '    
          READF,scandate_lun,line
          IF(STRCMP(line,'scalar product',14,/FOLD_CASE)EQ 1)THEN BEGIN
              READF, scandate_lun,line
              retarr=DOUBLE(STRTRIM(STRSPLIT(line,' ',/EXTRACT)))
              retval=NORM( retarr, /DOUBLE, LNORM=2)
              READF, scandate_lun,line
              WHILE(STRCMP(line,'eigenvalues',11,/FOLD_CASE)EQ 0) DO BEGIN
                  line=STRTRIM(line)
                  retarr=DOUBLE(STRTRIM(STRSPLIT(line,' ',/EXTRACT)))
                  retval=[retval,NORM( retarr, /DOUBLE, LNORM=2) ]
                  READF, scandate_lun,line
              ENDWHILE
          ENDIF ELSE BEGIN
              print, 'no scalar products calculated'
              retval=FLOAT('NaN')
          ENDELSE
          FREE_LUN, scandate_lun
      ENDELSE
  ENDIF ELSE BEGIN
      OPENR,scandate_lun,gui.out.data_path+'parameters_'+$
            run_label,/GET_LUN,ERROR=err
      IF err NE 0 THEN BEGIN
          print, 'could not open file'
          retval=FLOAT('NaN')
      ENDIF ELSE BEGIN
          line=' '    
          WHILE NOT EOF(scandate_lun) DO BEGIN
              READF,scandate_lun,line
              line=STRTRIM(line)
              linarr=STRTRIM(STRSPLIT(line,'=',COUNT=found,/EXTRACT))
              IF (found NE 2) THEN BEGIN
                  linarr=STRTRIM(STRSPLIT(line,':',COUNT=found,/EXTRACT))
                  IF (found GE 2) THEN BEGIN
                      line=STRTRIM(STRSPLIT(linarr[1],' ',/EXTRACT))
                      linarr[1]=line[0]
                  ENDIF
              ENDIF
              
              IF (STRCMP(linarr[0],title, /FOLD_CASE) EQ 1) THEN BEGIN
                  retval=DOUBLE(STRTRIM(linarr[1]))
              ENDIF
          ENDWHILE
          FREE_LUN, scandate_lun
      ENDELSE
  ENDELSE      
  RETURN, retval
END

;##########################################################################
;##########################################################################

PRO scan_plotprep, diag, k, xax, yax, zax, plotorderr, parindex, mytitlex,$
  mytitley, mytitlez, mytitletot

  COMMON global_vars
  i = (*diag).internal_vars

  FOR m = 1,N_ELEMENTS((*scanlog).dsize) - 1 DO BEGIN
    IF (*scanlog).dsize[m] EQ 1 THEN BEGIN
        (*i).plotparA[parindex]='0:0'
    ENDIF
    IF STRLEN((*i).plotparA[parindex]) EQ 0  THEN (*i).plotparA[parindex]='*'
    IF (*i).plotparA[parindex] EQ '*' THEN BEGIN
        (*i).plotparA[parindex]='0:'+$
          STRING(rm0es((*scanlog).dsize[parindex+1]-1))
    ENDIF
    interval=$
      FIX(STRTRIM(STRSPLIT((*i).plotparA[parindex],COUNT=inter,':',/EXTRACT)))
    IF inter EQ 1 THEN BEGIN
        (*i).plotparA[parindex]=$
          rm0es(STRING(interval[0]))+':'+rm0es(STRING(interval[0]))
        k=k-1
        mytitletot=' ('+(*scanlog).dimarr[m]+'='+$
                    rm0es(STRING((*scanlog).rd[m-1]-((*scanlog).dsize[m]-$
                    interval[0]-1)*(*scanlog).deltad[m-1])) +')'+mytitletot
    ENDIF ELSE BEGIN
        IF interval[0] EQ interval[1] THEN BEGIN
            (*i).plotparA[parindex]=rm0es(STRING(interval[0]))+':'+$
              rm0es(STRING(interval[0]))
            k=k-1
            mytitletot=' ('+(*scanlog).dimarr[m]+'='+$
              rm0es(STRING((*scanlog).rd[m-1]-$
              ((*scanlog).dsize[m]-interval[0]-1)*(*scanlog).deltad[m-1])) +$
              ')'+mytitletot
        ENDIF ELSE BEGIN
            CASE k OF
                0: BEGIN
                    mytitlex=(*scanlog).dimarr[m]
                    xax=(*scanlog).axxis[m-1,interval[0]:interval[1]]
                END
                1: BEGIN
                    mytitley=(*scanlog).dimarr[m]
                    yax=(*scanlog).axxis[m-1,interval[0]:interval[1]]
                END
                2: BEGIN
                    mytitlez=(*scanlog).dimarr[m]
                    zax=(*scanlog).axxis[m-1,interval[0]:interval[1]]
                END
                3: BEGIN
                    print,'too much parameters, set to 0',k
                    (*i).plotparA[parindex]='0:0'
                    mytitletot=' ('+(*scanlog).dimarr[m]+'='+$
                      rm0es(STRING((*scanlog).rd[m-1]-((*scanlog).dsize[m]-1)*$
                      (*scanlog).deltad[m-1])) +')'+mytitletot
                    k=k-1
                END
            ENDCASE
        ENDELSE
    ENDELSE
    k=k+1
    plotorderr=STRING(plotorderr)+STRING(',')+STRING(((*i).plotparA[parindex]))
    parindex+=1
  ENDFOR

END

;######################################################################

PRO scan_plotsetup, diag, eigeninterval, again

  COMMON global_vars

  i = (*diag).internal_vars

  IF (*i).eigenvalues EQ '*' THEN BEGIN
    again = (*scanlog).dsize(0) + 1
    eigeninterval[0] = 0
    eigeninterval[1] = (*scanlog).dsize(0) - 1
  ENDIF ELSE IF STRPOS((*i).eigenvalues,':') NE -1 THEN BEGIN
    eigeninterval = FIX(STRTRIM(STRSPLIT((*i).eigenvalues,$
      COUNT=eigeninter,':',/EXTRACT)))
    IF eigeninter EQ 2 THEN $
      again = eigeninterval[1] - eigeninterval[0] + 2
  ENDIF

END

;######################################################################

PRO plotiplot, xax, yax, plotarr, start

  CASE start OF
      0: IPLOT, xax, yax, plotarr, OVERPLOT=0,/NO_SAVEPROMPT, COLOR=[30,200,50], $
        NAME='RIEMANNSURFACE', TITLE=' !7c!6 '
      2: IPLOT, xax, yax, plotarr, OVERPLOT=1, /NO_SAVEPROMPT, COLOR=[30,200,50], $
        NAME='RIEMANNSURFACE', TITLE=' !7c!6 ';, /VIEW_NEXT
      ELSE: PRINT, 'error with iplot'
  ENDCASE

END

;######################################################################

FUNCTION riemsort, sortarr, refarr

  IF NORM(refarr-sortarr) GT NORM(refarr-[sortarr[1],sortarr[0]]) THEN BEGIN
    riemx = sortarr[1]
    sortarr[1] = sortarr[0]
    sortarr[0] = riemx
  ENDIF

  RETURN, sortarr

END

;######################################################################

PRO riemplot, diag, riemeigenfield, arrsize

  COMMON global_vars
;  i = (*diag).internal_vars

  sind=0
  dsizepos=MAKE_ARRAY(2,/BYTE)
  FOR ind=1, N_ELEMENTS(arrsize) -1 DO BEGIN
      IF arrsize[ind] GT 1 THEN BEGIN
          dsizepos[sind]=ind
          sind+=1
      ENDIF
  ENDFOR
  riemarr=riemeigenfield
  OPENW,scanriem_lun,gui.out.out_path+'riemann.log',/GET_LUN,ERROR=err
  IF err NE 0 THEN print, 'could not open file'
  PRINTF,scanriem_lun, '# '+(*scanlog).dimarr[dsizepos[0]]+' '+(*scanlog).dimarr[dsizepos[1]]+' remod1 remod2 immod1 immod2'
  sizearr=0
  sizearr=SIZE(riemarr,/DIMENSIONS)
  withiplot=0
     IF withiplot EQ 1 THEN BEGIN
  firstplot=1
      ENDIF
;first dimension mode 0
  FOR ri=0, sizearr[1]-1 DO BEGIN
      IF (ri GT 0) THEN riemarr[*,ri,0]=riemsort(riemarr[*,ri,0], riemarr[*,ri-1,0])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,sizearr[2]],/DOUBLE)
      ENDIF
      FOR rj=0, sizearr[2]-1 DO BEGIN
          IF (rj GT 0) THEN riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri,rj-1])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,rj]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      IF firstplot EQ 1 THEN BEGIN
          firstplot=0
          plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),0
          plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),1
      ENDIF ELSE BEGIN
          plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
          plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDELSE
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
  ENDFOR
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
;first dimension mode 1
  FOR ri=0, sizearr[1]-1 DO BEGIN
      IF (ri GT 0) THEN riemarr[*,ri,0]=riemsort(riemarr[*,ri,0], riemarr[*,ri-1,0])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,sizearr[2]],/DOUBLE)
      ENDIF
      FOR rj=0, sizearr[2]-1 DO BEGIN
          IF (rj GT 0) THEN riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri,rj-1])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,rj]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
  ENDFOR
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

;second dimension mode 0
  FOR rj=0, sizearr[2]-1 DO BEGIN
      IF (rj GT 0) THEN riemarr[*,0,rj]=riemsort(riemarr[*,0,rj], riemarr[*,0,rj-1])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,sizearr[1]],/DOUBLE)
      ENDIF
      FOR ri=0, sizearr[1]-1 DO BEGIN
          IF (ri GT 0) THEN riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri-1,rj])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,ri]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
  ENDFOR

  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

;second dimension mode 1
  FOR rj=0, sizearr[2]-1 DO BEGIN
      IF (rj GT 0) THEN riemarr[*,0,rj]=riemsort(riemarr[*,0,rj], riemarr[*,0,rj-1])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,sizearr[1]],/DOUBLE)
      ENDIF
      FOR ri=0, sizearr[1]-1 DO BEGIN
          IF (ri GT 0) THEN riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri-1,rj])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,ri]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
  ENDFOR

  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

;first dimension diag mode 0
;(first sort to xrj=sizearr[2]-1)
  FOR rj=0, sizearr[2]-1 DO BEGIN
      IF (rj GT 0) THEN riemarr[*,0,rj]=riemsort(riemarr[*,0,rj], riemarr[*,0,rj-1])
  ENDFOR
;start
  xri=0
  xrj=sizearr[2]-1
  WHILE(xri LE sizearr[1]-1) DO BEGIN
      PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[0,xri,xrj]), IMAGINARY(riemarr[0,xri,xrj])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,MIN([sizearr[2]-1-xrj,sizearr[1]-1-xri])+1],/DOUBLE)
      iplotarr[*,0]=[(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[0,xri,xrj]), IMAGINARY(riemarr[0,xri,xrj])]
      ENDIF
      FOR xi=1, MIN([sizearr[2]-1-xrj,sizearr[1]-1-xri]) DO BEGIN
          ri=xri+xi
          rj=xrj+xi
          riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri-1,rj-1])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,xi]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
      IF (xrj EQ 0) THEN BEGIN
          xri=xri+1
          IF (xri LE sizearr[1]-1) THEN riemarr[*,xri,0]=riemsort(riemarr[*,xri,0], riemarr[*,xri-1,0])
      ENDIF ELSE BEGIN
          xrj=xrj-1
      ENDELSE
  ENDWHILE

  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

;first dimension diag mode 1
;(first sort to xrj=sizearr[2]-1)
  FOR rj=0, sizearr[2]-1 DO BEGIN
      IF (rj GT 0) THEN riemarr[*,0,rj]=riemsort(riemarr[*,0,rj], riemarr[*,0,rj-1])
  ENDFOR
;start
  xri=0
  xrj=sizearr[2]-1
  WHILE(xri LE sizearr[1]-1) DO BEGIN
      PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[1,xri,xrj]), IMAGINARY(riemarr[1,xri,xrj])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,MIN([sizearr[2]-1-xrj,sizearr[1]-1-xri])+1],/DOUBLE)
      iplotarr[*,0]=[(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[1,xri,xrj]), IMAGINARY(riemarr[1,xri,xrj])]
      ENDIF
      FOR xi=1, MIN([sizearr[2]-1-xrj,sizearr[1]-1-xri]) DO BEGIN
          ri=xri+xi
          rj=xrj+xi
          riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri-1,rj-1])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,xi]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
      IF (xrj EQ 0) THEN BEGIN
          xri=xri+1
          IF (xri LE sizearr[1]-1) THEN riemarr[*,xri,0]=riemsort(riemarr[*,xri,0], riemarr[*,xri-1,0])
      ENDIF ELSE BEGIN
          xrj=xrj-1
      ENDELSE
  ENDWHILE

  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

;second dimension diag mode 0
;(first sort to xrj=sizearr[2]-1)
  FOR rj=0, sizearr[2]-1 DO BEGIN
      IF (rj GT 0) THEN riemarr[*,0,rj]=riemsort(riemarr[*,0,rj], riemarr[*,0,rj-1])
  ENDFOR
;start
  xri=0
  xrj=0
  WHILE(xrj LT sizearr[2]-1) DO BEGIN
      PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[0,xri,xrj]), IMAGINARY(riemarr[0,xri,xrj])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,MIN([sizearr[2]-1-xrj,xri])+1],/DOUBLE)
      iplotarr[*,0]=[(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[0,xri,xrj]), IMAGINARY(riemarr[0,xri,xrj])]
      ENDIF
      FOR xi=1, MIN([sizearr[2]-1-xrj,xri]) DO BEGIN
          ri=xri-xi
          rj=xrj+xi
          riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri+1,rj-1])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,xi]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[0,ri,rj]), IMAGINARY(riemarr[0,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
      IF (xri LT sizearr[1]-1) THEN BEGIN
          xri=xri+1
          riemarr[*,xri,0]=riemsort(riemarr[*,xri,0], riemarr[*,xri-1,0])
      ENDIF ELSE BEGIN
          xrj=xrj+1
          riemarr[*,sizearr[1]-1,xrj]=riemsort(riemarr[*,sizearr[1]-1,xrj], riemarr[*,sizearr[1]-1,xrj-1])
      ENDELSE
  ENDWHILE

  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

;second dimension diag mode 1
  xri=0
  xrj=0
  WHILE(xrj LT sizearr[2]-1) DO BEGIN
      PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[1,xri,xrj]), IMAGINARY(riemarr[1,xri,xrj])
     IF withiplot EQ 1 THEN BEGIN
      iplotarr=MAKE_ARRAY([4,MIN([sizearr[2]-1-xrj,xri])+1],/DOUBLE)
      iplotarr[*,0]=[(*scanlog).axxis[dsizepos[0]-1,xri],(*scanlog).axxis[dsizepos[1]-1,xrj], FLOAT(riemarr[1,xri,xrj]), IMAGINARY(riemarr[1,xri,xrj])]
      ENDIF
      FOR xi=1, MIN([sizearr[2]-1-xrj,xri]) DO BEGIN
          ri=xri-xi
          rj=xrj+xi
          riemarr[*,ri,rj]=riemsort(riemarr[*,ri,rj], riemarr[*,ri+1,rj-1])
          PRINTF, scanriem_lun,(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])
     IF withiplot EQ 1 THEN BEGIN
          iplotarr[*,xi]=[(*scanlog).axxis[dsizepos[0]-1,ri],(*scanlog).axxis[dsizepos[1]-1,rj], FLOAT(riemarr[1,ri,rj]), IMAGINARY(riemarr[1,ri,rj])]
      ENDIF
      ENDFOR
     IF withiplot EQ 1 THEN BEGIN
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[2,*]),2
      plotiplot,REFORM(iplotarr[0,*]),REFORM(iplotarr[1,*]),REFORM(iplotarr[3,*]),3
      ENDIF
      PRINTF,scanriem_lun, ' '
      PRINTF,scanriem_lun, ' '
      IF (xri LT sizearr[1]-1) THEN BEGIN
          xri=xri+1
          riemarr[*,xri,0]=riemsort(riemarr[*,xri,0], riemarr[*,xri-1,0])
      ENDIF ELSE BEGIN
          xrj=xrj+1
          riemarr[*,sizearr[1]-1,xrj]=riemsort(riemarr[*,sizearr[1]-1,xrj], riemarr[*,sizearr[1]-1,xrj-1])
      ENDELSE
  ENDWHILE

  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, '#'
  PRINTF,scanriem_lun, ' '
  PRINTF,scanriem_lun, ' '

  FREE_LUN, scanriem_lun

END

;######################################################################

FUNCTION scan_genplotfield,field2plot,plotorder,dsize,test

  COMMON global_vars

  field2plot2=field2plot
  arrsize=dsize
  plotorderarr=STRSPLIT(plotorder,',',COUNT=dims,/EXTRACT)
  FOR index=dims-1,0,-1 DO BEGIN
      multiplicl=1
      multiplicr=1
      IF index GT 0 THEN BEGIN
          FOR jndex=0,index-1 DO BEGIN
              multiplicl=multiplicl*arrsize[jndex]
          ENDFOR
      ENDIF ELSE BEGIN
          IF test EQ 0 THEN BEGIN
          ;Riemanntest-----------------------------------------
              riemarr=REFORM(field2plot2,arrsize,/OVERWRITE)
              riemarr=REFORM(riemarr)
              IF (SIZE(riemarr,/N_DIMENSIONS) EQ 3)  THEN BEGIN
                  sizearr=SIZE(riemarr,/DIMENSIONS)
                  ri=0
                  rj=0
                  WHILE(ri LT sizearr[1]-1)DO BEGIN
                      riemarr[*,ri+1,0]=riemsort(riemarr[*,ri+1,0], riemarr[*,ri,0])
                      ri+=1
                  ENDWHILE
                  WHILE(rj LT sizearr[2]-1)DO BEGIN
                      riemarr[*,ri,rj+1]=riemsort(riemarr[*,ri,rj+1], riemarr[*,ri,rj])
                      rj+=1
                  ENDWHILE
                  WHILE(ri GT 0)DO BEGIN
                      riemarr[*,ri-1,rj]=riemsort(riemarr[*,ri-1,rj], riemarr[*,ri,rj])
                      ri-=1
                  ENDWHILE
                  WHILE(rj GT 1)DO BEGIN
                      riemarr[*,ri,rj-1]=riemsort(riemarr[*,ri,rj-1], riemarr[*,ri,rj])
                      rj-=1
                  ENDWHILE
                  IF (NORM(riemarr[*,ri,rj-1]-riemarr[*,ri,rj]) GT $
                    NORM(riemarr[*,ri,rj-1]-[riemarr[1,ri,rj],riemarr[0,ri,rj]])) THEN BEGIN

                      print, 'OUTPUT WARNING: probably a Riemann surface; look at "riemann.log"!'
                      riemplot,diag,riemarr,arrsize
                      startstr=' '

                      str="echo "+"'"+'set pm3d; splot '+'"'+gui.out.out_path+'riemann.log"'+$
                        " u 1:2:3  w l lw 1 palette' | gnuplot -persist"
                      SPAWN, str
                      str="echo "+"'"+'set pm3d; splot '+'"'+gui.out.out_path+'riemann.log"'+$
                        " u 1:2:4  w l lw 1 palette' | gnuplot -persist"
                      SPAWN, str
                  ENDIF ELSE BEGIN
                  ENDELSE
              ENDIF
          ENDIF
      ;--------------------------------------------
      ENDELSE
      IF index LT dims-1 THEN BEGIN
          FOR jndex=index+1,dims-1 DO BEGIN
              multiplicr=multiplicr*arrsize[jndex]
          ENDFOR
      ENDIF
      range=STRSPLIT(plotorderarr[index],':',COUNT=ddot,/EXTRACT)
      IF ddot EQ 1 THEN range=[range,range]
      range=FIX(rm0es(range))
      field2plot1=DCOMPLEXARR(multiplicl*multiplicr*(range[1]-range[0]+1))
      FOR kndex=0,multiplicr-1 DO BEGIN
          field2plot1[kndex*multiplicl*(range[1]-range[0]+1):(kndex+1)*$
            multiplicl*(range[1]-range[0]+1)-1]=field2plot2[multiplicl*$
            (range[0]+kndex*arrsize[index]):multiplicl*$
            (kndex*arrsize[index]+range[1]+1)-1]
      ENDFOR
      field2plot2=field2plot1
      arrsize[index]=range[1]-range[0]+1
  ENDFOR
  plotfield=REFORM(field2plot2,arrsize,/OVERWRITE)
  return, plotfield
END

;######################################################################

PRO scan_plot, diag, totaltitle, field2plot, reimabs, outa, outb, outc,$
               lega, legb, legc, legd, outset, sp = sp, suffix=suffix

;totaltitle
;outset: 2=only one call of plot ; 1/-1/0  first/last/no call of series
;of plot calls
;sp: species index (for output file names)

  COMMON global_vars

  i = (*diag).internal_vars

  IF NOT KEYWORD_SET(suffix) THEN suffix=''
  IF (*i).eigenvalues EQ '' THEN (*i).eigenvalues='*'
  IF (*scanlog).dsize[0]EQ 1 THEN (*i).eigenvalues='0'

  againsave=(*i).eigenvalues

  coltablearr=MAKE_ARRAY(1,/BYTE,VALUE=0)
  captionsarr=MAKE_ARRAY(1,/STRING,VALUE='')
  again=1
  eigeninterval=[0,0]
  scan_plotsetup, diag, eigeninterval, again

  reim=(1+eigeninterval[1]-eigeninterval[0])*2
  outputset=0

  filearr=0.

  riemanntest=0

  WHILE reim GT 0 DO BEGIN
      mytitletot=totaltitle
      str=''
      mytitle=''
      mytitlex=''
      mytitley=''
      mytitlez=''
      again=1
      colrun=0
      eigeninterval=[0,0]
      scan_plotsetup,diag,  eigeninterval, again

      REPEAT BEGIN
          plotorderl=''
          IF again GT 1 THEN BEGIN
              againsave=STRING(eigeninterval[1]-again+2)
          ENDIF
          plotorderl+=againsave
          mytitle='mode '+STRTRIM(STRING(againsave),2)
          IF (again EQ 1)OR(again EQ (*scanlog).dsize(0)+1) OR $
            ((eigeninterval[1] NE 0)AND(again EQ $
            eigeninterval[1]-eigeninterval[0]+2)) THEN BEGIN
              k=0
              xax=1
              yax=1
              zax=1
              plotorderr=''
              parindex=0
              scan_plotprep, diag, k, xax, yax, zax, plotorderr,$
                parindex, mytitlex, mytitley, mytitlez,mytitletot
          ENDIF
          plotorder=plotorderl+plotorderr
          plotfield=scan_genplotfield(field2plot,plotorder,$
                    (*scanlog).dsize,riemanntest)
          riemanntest=1
          IF STRCMP(reimabs,'abs',3,/FOLD_CASE) THEN BEGIN
              plotfield=ABS(plotfield)
              plotmin=0
              plotmax=MAX(ABS(field2plot))
              IF (again EQ 1)OR(again EQ (*scanlog).dsize(0)+1) OR $
                ((eigeninterval[1] NE 0)AND $
                (again EQ eigeninterval[1]-eigeninterval[0]+2))THEN $
                     amplitude=' abs(!7x!6!Dcomplex!N) '
          ENDIF ELSE BEGIN
              IF STRCMP(reimabs,'non',3,/FOLD_CASE) THEN BEGIN
		  plotfield=FLOAT(plotfield)
                  plotmin=MIN(field2plot)
                  plotmax=MAX(field2plot)
              ENDIF ELSE BEGIN
                  IF reim GT (eigeninterval[1]-eigeninterval[0]+1) THEN BEGIN
                      IF (again EQ 1)OR(again EQ (*scanlog).dsize(0)+1) OR $
                        ((eigeninterval[1] NE 0)AND(again EQ eigeninterval[1]$
                        -eigeninterval[0]+2))THEN amplitude= ' !7c!6 '
                      plotfield=FLOAT(plotfield)
                      plotmin=0
                      plotmax=MAX(FLOAT(field2plot),/NAN)
                  ENDIF ELSE BEGIN
                      IF (again EQ 1)OR(again EQ (*scanlog).dsize(0)+1) OR $
                        ((eigeninterval[1] NE 0)AND(again EQ eigeninterval[1]$
                        -eigeninterval[0]+2))THEN amplitude=' !7x!6!r!a-!N '
                      plotfield=IMAGINARY(plotfield)
                      plotmin=MIN(IMAGINARY(field2plot),/NAN)
                      plotmax=MAX(IMAGINARY(field2plot),/NAN)
                  ENDELSE
              ENDELSE
          ENDELSE
          plotfield=REFORM(plotfield) ;discard dipensable dimensions
          xax=REFORM(xax) ;discard dipensable dimensions
          yax=REFORM(yax) ;discard dipensable dimensions
          zax=REFORM(zax) ;discard dipensable dimensions
          same=1

          IF (again EQ 1)OR(again EQ (*scanlog).dsize(0)+1) OR $
            ((eigeninterval[1] NE 0)AND(again EQ eigeninterval[1]- $
            eigeninterval[0]+2))THEN BEGIN
              same=0
              colrun=0
          ENDIF ELSE BEGIN
              same = 1
          ENDELSE

          case colrun OF
              0: color=[225,184,0]
              1: color=[0,225,184]
              2: color=[184,0,225]
              3: color=[184,225,0]
              4: color=[0,184,225]
              5: color=[225,0,184]
              ELSE : color=[0,0,255]
          ENDCASE
          CASE k OF
              0: BEGIN
                  IF STRCMP(reimabs,'abs',3,/FOLD_CASE)$
                    OR STRCMP(reimabs,'non',3,/FOLD_CASE) THEN BEGIN
                      reim-=2
                  ENDIF ELSE BEGIN
                      reim-=1
                  ENDELSE
              END
              1: BEGIN
                  IF STRCMP(reimabs,'abs',3,/FOLD_CASE)$
                    OR STRCMP(reimabs,'non',3,/FOLD_CASE) THEN BEGIN
                      IF coltablearr[0] EQ 0 THEN BEGIN
                          color1=2
                          coltablearr[0]=2
                          captionsarr[0]=mytitle
                      ENDIF ELSE BEGIN
                          color1=N_ELEMENTS(coltablearr)+2
                          coltablearr=[coltablearr,color1]
                          captionsarr=[captionsarr,mytitle]
                      ENDELSE
                      ;                   IF again NE 0 THEN BEGIN
                  ENDIF ELSE BEGIN
                      IF reim GT 1+eigeninterval[1]-eigeninterval[0] THEN BEGIN
                          IF coltablearr[0] EQ 0 THEN BEGIN
                              coltablearr[0]=(reim MOD (*scanlog).dsize(0)) +2
                              captionsarr[0]=mytitle
                          ENDIF ELSE BEGIN
                              coltablearr=[coltablearr,(reim MOD $
                                                        (*scanlog).dsize(0)) +2]
                              captionsarr=[captionsarr,mytitle]
                          ENDELSE
                      ENDIF
                      color1=(reim MOD (*scanlog).dsize[0]) +2
                  ENDELSE

                  IF outputset EQ 0 THEN BEGIN
                    IF (outset EQ 1) OR (outset EQ 2) THEN BEGIN
                       set_output, diag, sp, /ps, multi=[outa, outb, outc], suffix=suffix
                      ENDIF
                      outputset=1
                  ENDIF
                  IF SIZE(filearr,/N_ELEMENTS) EQ 1 THEN BEGIN
	              filearr=[[xax],[plotfield]]
                      filenames=[mytitlex,mytitle+STRING(reim MOD (*scanlog).dsize[0]) ]
                  ENDIF ELSE BEGIN
                      filearr=[[filearr],[plotfield]]
                      filenames=[filenames,mytitle+STRING(reim MOD (*scanlog).dsize[0])]
                  ENDELSE
                  IF same EQ 1 THEN BEGIN
                      OPLOT, xax, plotfield,COLOR=color1,PSYM=-color1*2+2
                  ENDIF ELSE BEGIN
                      PLOT, xax,plotfield,/NODATA,/XSTYLE,$
                            XTITLE=mytitlex, YTITLE=amplitude,COLOR=1,$
                            TITLE=mytitletot,YRANGE=[plotmin,plotmax]
                      OPLOT, !X.CRANGE, [0,0],COLOR=1
                      OPLOT, xax, plotfield,COLOR=color1,PSYM=-color1*2+2
                  ENDELSE
                  IF STRCMP(reimabs,'abs',3,/FOLD_CASE)$
                    OR STRCMP(reimabs,'non',3,/FOLD_CASE) THEN BEGIN
                      reim-=2
                  ENDIF ELSE BEGIN
                      reim-=1
                  ENDELSE
              END
              2: BEGIN
                  IF STRCMP(reimabs,'abs',3,/FOLD_CASE)$
                    OR STRCMP(reimabs,'non',3,/FOLD_CASE) THEN BEGIN
                      reim-=2
                  ENDIF ELSE BEGIN
                      reim-=1
                  ENDELSE
                  IF SIZE(filearr,/N_ELEMENTS) EQ 1 THEN BEGIN
                      filearr=MAKE_ARRAY([N_ELEMENTS(xax)*N_ELEMENTS(yax),3])
                      FOR index = 0, N_ELEMENTS(xax)-1 DO BEGIN
                          FOR jndex = 0, N_ELEMENTS(yax)-1 DO BEGIN
                              filearr[index+jndex*N_ELEMENTS(xax),0:2]$
                                =[xax[index],yax[jndex],plotfield[jndex*N_ELEMENTS(xax)+index]]
                              filenames=[mytitlex,mytitley,mytitle+STRING(reim MOD (*scanlog).dsize[0]) ]
                          ENDFOR
                      ENDFOR
                  ENDIF ELSE BEGIN
                      filearr=[[filearr],[plotfield[0:N_ELEMENTS(plotfield)-1]]]
                      filenames=[filenames,mytitle+STRING(reim MOD (*scanlog).dsize[0])]
                  ENDELSE

                  isurface,plotfield,xax,yax,NAME=mytitle,COLOR=color,$
                            INSERT_LEGEND=[-1,1],XTITLE=mytitlex,$
                            ZTITLE=amplitude,YTITLE=mytitley,$
                            TITLE=mytitletot,TRANSPARENCY=20,$
                            /USE_TRIANGLES,IDENTIFIER=ident,OVERPLOT=same,$
			    /NO_SAVEPROMPT,/DISABLE_SPLASH_SCREEN,$
			    XTICKFONT_INDEX=4,YTICKFONT_INDEX=4,$
			    ZTICKFONT_INDEX=4
               END
              3: BEGIN
                  IF STRCMP(reimabs,'abs',3,/FOLD_CASE) OR $
                    STRCMP(reimabs,'non',3,/FOLD_CASE) THEN BEGIN
                      reim-=2
                  ENDIF ELSE BEGIN
                      reim-=1
                  ENDELSE
                  IF reim EQ (*scanlog).dsize(0) THEN reim = 0
                  maxval=MAX(plotfield,MIN=minval,/NAN)
                  plotfieldcol=(plotfield-minval)*255/(maxval-minval)
                  plotfieldcol=BYTE(plotfieldcol)
                  plotfieldcol=(plotfieldcol+1) mod 256

                  plotptr=PTR_NEW(plotfieldcol,/NO_COPY)
                  SLICER3,plotptr,DATA_NAMES=mytitle,/MODAL
               END
          ENDCASE
          again=again-1
          colrun=colrun+1
      ENDREP UNTIL again LE 1

  ENDWHILE
  IF outputset EQ 1 THEN BEGIN
      plot_legend, coltablearr, captionsarr,per_line=(*scanlog).dsize[0],$
                   x_offset=[lega, legb, legc, legd]
      IF (outset EQ -1) OR (outset EQ 2) THEN BEGIN
          set_output, diag, sp, /reset , suffix=suffix
      ENDIF
  ENDIF

  set_output, diag, sp, header = filenames,$
              dat=[filearr], suffix=suffix
END

;######################################################################

FUNCTION scan_plotnec, diag, plotpars, eigenrunnum
; returns 1, if data are relevant for plotting

  COMMON global_vars
  i = (*diag).internal_vars

  number=eigenrunnum-1
  FOR dim=0, N_ELEMENTS((*scanlog).dsize)-1 DO BEGIN
      IF NOT (((*scanlog).dsize[dim] EQ 0)OR(plotpars[dim] EQ '*')) THEN BEGIN
          IF (STRLEN(plotpars[dim]) GT 0) THEN BEGIN
              test=number MOD (*scanlog).dsize(dim)
              border=FIX(STRTRIM(STRSPLIT(plotpars[dim],COUNT=inter,':',/EXTRACT)))
              IF inter EQ 1 THEN border =[border,border]
              IF ((test LT border[0])OR(test GT border[1])) THEN BEGIN
                  return,0
              ENDIF
          ENDIF
          number=FLOOR(number / (*scanlog).dsize(dim))
      ENDIF
  ENDFOR
  return,1

END

;##########################################################################
