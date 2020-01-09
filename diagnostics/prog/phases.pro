FUNCTION phases_info

  RETURN, {$
    type      : 'mom', $;'mom_uni',$
    title     : 'Cross phases',$
    help_text : ['Plots the PDFs of the cross phases between pairs '+$
                 'of variables in ky space.'],$
    ext_vars  : [['vars','0','variable pairs for cross phase '+$
                  'calculation, e.g. [0,nf,0,nf+1] for phi x n and '+$
                  'phi x Tpar cross phases; default: depends on '+$
                  'n_fields and n_moms'],$
                 ['xind','0','x indices for cross phase '+$
                  'calculation; default: -1 (average)'],$
                 ['zind','0','z point for cross phase '+$
                  'calculation; default: -1 (average)'],$
                 ['log','1','toggle logarithmic scale in ky'],$
                 ['weight','1','toggle amplitude weighting in '+$
                  'ky direction']]}

END

;######################################################################

PRO phases_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF NOT KEYWORD_SET(*e.weight) THEN *e.weight = 0
  IF NOT KEYWORD_SET(*e.log) THEN *e.log = 0
  IF N_ELEMENTS(*e.xind) LT 1 THEN xind = INDGEN(par.nkx0) ELSE $
    CASE (*e.xind)[0] OF
      -1 : xind = INDGEN(par.nkx0)
      ELSE : xind = *e.xind
    ENDCASE

  IF (par.nky0 EQ 1) AND (par.ky0_ind EQ 0) THEN BEGIN
    printerror, 'Skipping + ' + (*diag).name + ': no non-zero ky modes'
    (*diag).selected = 0
    RETURN
  ENDIF ELSE kyind = par.ky0_ind > 1 + INDGEN(par.nky0-(par.ky0_ind>1))

  IF N_ELEMENTS(*e.zind) LT 1 THEN zinds = INDGEN(par.nz0) ELSE $
    CASE (*e.zind)[0] OF
      -1 : zinds = INDGEN(par.nz0)
      -2 : zinds = par.nz0/2
      ELSE : zinds = *e.zind
    ENDCASE  

  IF !QUIET NE 1 THEN IF *e.weight THEN PRINT, (*diag).name + $
    ': weighing by amplitude' ELSE PRINT, (*diag).name + ': no weighing'

  binsize = 0.1
  nbins = FIX(2*!PI/binsize)

  nf = par.n_fields
  IF N_ELEMENTS(*e.vars) LE 1 THEN BEGIN
    IF par.n_moms GT 6 THEN BEGIN
      CASE nf OF
        1 : BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2],[0,nf+6],[0,nf+6+1],[0,nf+6+2]]
          vars = [0,nf,nf+1,nf+2,nf+6,nf+6+1,nf+6+2]
        END
        2 : BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2],[0,nf+6],[0,nf+6+1],[0,nf+6+2],$
            [1,nf+5],[1,nf+3],[1,nf+6+4],[1,nf+6+5],[1,nf+6+3],[1,nf+6+4]]
          vars = [0,1,nf,nf+1,nf+2,nf+3,nf+5,nf+6,nf+6+1,nf+6+2,nf+6+3,$
            nf+6+4,nf+6+5]    
        END
        3 : BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2],[1,nf+5],[1,nf+3],[1,nf+4],$
            [2,1],[2,nf+3],[2,nf+4]]
          vars = [0,1,2,nf,nf+1,nf+2,nf+3,nf+4,nf+5]
        END
      ENDCASE
    ENDIF ELSE BEGIN
      CASE nf OF
        1 : BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2]]
          vars = [0,nf,nf+1,nf+2]
        END
        2 : IF par.beta GT 1e-5 THEN BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2],[1,nf+5],[1,nf+3],[1,nf+4]]
          vars = [0,1,nf,nf+1,nf+2,nf+3,nf+4,nf+5]
        ENDIF ELSE BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2]]
          vars = [0,nf,nf+1,nf+2]
        ENDELSE
        3 : BEGIN
          corr = [[0,nf],[0,nf+1],[0,nf+2],[1,nf+5],[1,nf+3],[1,nf+4],$
            [2,1],[2,nf+3],[2,nf+4]]
          vars = [0,1,2,nf,nf+1,nf+2,nf+3,nf+4,nf+5]           
        END          
      ENDCASE
    ENDELSE
  ENDIF ELSE BEGIN
    corr = INTARR(2,N_ELEMENTS(*e.vars)/2)
    vars = FIX(*e.vars)

    FOR v = 0, N_ELEMENTS(*e.vars) / 2 - 1 DO $
      corr[*,v] = vars[2*v:2*v+1]
  ENDELSE

  i = set_internal_vars(diag,{$
    weight     : *e.weight,$
    log        : *e.log,$
    xind       : xind,$
    kyind      : kyind,$
    zind       : zinds,$
    nbins      : nbins,$
    binsize    : binsize,$
    sumhist_id : PTR_NEW(),$
    corr       : corr})

  PRINT, (*diag).name + ': using binsize = ', rm0es(binsize), $
    ', resulting in ', rm0es(nbins), ' bins'

  fft_format, sxky=vars

END

;######################################################################

PRO phases_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  n_corrs = N_ELEMENTS((*i).corr) / 2
  nxind = N_ELEMENTS((*i).xind)
  nkyind = N_ELEMENTS((*i).kyind)
  nzind = N_ELEMENTS((*i).zind)

  hist = FLTARR((*i).nbins+1,nkyind,/NOZERO)
  sumhist = FLTARR((*i).nbins+1,nkyind,n_corrs,gui.out.n_spec_sel,/NOZERO)
  IF gui.out.n_spec_sel LT 2 THEN sumhist = $
    REFORM(sumhist,[(*i).nbins+1,nkyind,n_corrs,gui.out.n_spec_sel],/OVERWRITE)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    FOR cv = 0, n_corrs - 1 DO BEGIN          
      var1 = (*i).corr[0,cv]
      var2 = (*i).corr[1,cv]

      tempdata1 = REFORM(COMPLEX((*mom[isp,var1].sxky)$
        [(*i).xind,(*i).kyind,(*i).zind,0]),[nxind,nkyind,nzind])
      tempdata2 = REFORM(COMPLEX((*mom[isp,var2].sxky)$
        [(*i).xind,(*i).kyind,(*i).zind,0]),[nxind,nkyind,nzind])

      IF (*i).weight NE 0 THEN BEGIN
        angle = ATAN(tempdata1/tempdata2,/PHASE)

        ; take product of abs values as weighting factor
        abstempdata1 = REFORM(ABS(TEMPORARY(tempdata1)),[nxind,nkyind,nzind])
        abstempdata2 = REFORM(ABS(TEMPORARY(tempdata2)),[nxind,nkyind,nzind])

        absval = REFORM(abstempdata1/REBIN(REFORM(TOTAL(TOTAL(abstempdata1,$
          3),1),[1,nkyind,1,1]),[nxind,nkyind,nzind,1])*abstempdata2/$
          REBIN(REFORM(TOTAL(TOTAL(abstempdata2,3),1),[1,nkyind,1,1]),$
          [nxind,nkyind,nzind,1]),[nxind,nkyind,nzind])
        FOR j = 0, nkyind - 1 DO BEGIN ; weighting with abs value of mode
          index = REFORM(FIX((angle[*,j,*]+!PI)/(*i).binsize))
          dummy = HISTOGRAM(index,MIN=0,MAX=FIX(2*!PI/(*i).binsize),$
            REVERSE_INDICES=ri)

          FOR b = 0, (*i).nbins DO hist[b,j] = (ri[b] EQ ri[b+1]) ? $
            0.0 : TOTAL((REFORM(absval[*,j,*]))[ri[ri[b]:ri[b+1]-1]])
        ENDFOR
      ENDIF ELSE BEGIN
        angle = ATAN(TEMPORARY(tempdata1)/TEMPORARY(tempdata2),/PHASE)

        FOR j = 0, nkyind - 1 DO hist[0,j] = HISTOGRAM(angle[*,j,*],$
          MIN=-!PI,MAX=!PI,/NAN,BINSIZE=(*i).binsize)
      ENDELSE

      hist = TEMPORARY(hist) / FLOAT(nxind*nzind)
      sumhist[0,0,cv,isp] = hist
    ENDFOR ; --- cv loop
  ENDFOR ; --- species loop

  (*i).sumhist_id = time_avg((*i).sumhist_id,sumhist,mom_time,$
    fwd=gui.out.res_steps)

END

;######################################################################

PRO phases_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  sumhist = time_avg((*i).sumhist_id,/avg,fwd=gui.out.res_steps,tarr=time)

  n_corrs = N_ELEMENTS((*i).corr) / 2
  cpp = 3 ; cross phases per page
  nkyind = N_ELEMENTS((*i).kyind)
  winkind = - !PI + INDGEN((*i).nbins+1) * 2.0 * !PI / (*i).nbins

  rho_str = ' ' + get_var_string(/rhostr)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]
   
    set_output, diag, sp, /ps, multi=[0,cpp<n_corrs,(nkyind GT 1)+1], $
      ysize=20, charsize=2.0

    FOR n = 0, gui.out.res_steps * (series.step_count - 1) DO BEGIN
      FOR page = 0, n_corrs / cpp + ((n_corrs MOD cpp) GT 0) - 1 DO BEGIN
        FOR cv = page * cpp, ((page + 1) * cpp - 1) < (n_corrs - 1) DO BEGIN
          var1 = (*i).corr[0,cv]
          var2 = (*i).corr[1,cv]
          maxz = MAX(sumhist[*,*,cv,isp])

          IF nkyind GT 1 THEN BEGIN
            SURFACE, sumhist[*,*,cv,isp], winkind, (*series.ky)[(*i).kyind], $
              /LEGO, COLOR=1, YLOG=(*i).log, XRANGE=[-!pi,!PI], /XSTYLE, $
              /YSTYLE, ZRANGE=[0,maxz], TITLE=get_var_string(var1,/fancy) + $
              ' x ' + get_var_string(var2,/fancy), YTITLE='!6k!Dy!N'+rho_str, $
              XTITLE='!7a!6', XTICKS=4, XTICKV=[-!PI,-!PI/2.0,0.0,!PI/2,!PI], $
              XTICKNAME=['-!7p!6','-!7p!6/2','0','!7p!6/2','!7p!6']
          ENDIF ELSE BEGIN
            PLOT, winkind, sumhist[*,0,cv,isp], XRANGE=[-!PI,!PI], /XSTYLE, $
              YTITLE=get_var_string(var1,/fancy) + ' x ' + $
              get_var_string(var2,/fancy) + ' at !6k!Dy!N' + rho_str + '=' + $
              rm0es((*par.ky)[(*i).kyind[0]]), XTITLE = '!7a!6', XTICKS=4, $
              XTICKV=[-!PI,-!PI/2.0,0.0,!PI/2,!PI], $
              XTICKNAME=['-!7p!6','-!7p!6/2','0','!7p!6/2','!7p!6'], COLOR=1

            avsum = 0.0
            avang = winkind[0]
            ; initial value
            FOR jav = 0, N_ELEMENTS(winkind) - 1 DO avsum += $
              sumhist[jav,0,cv,isp,n] * MIN([ABS(winkind[0]-winkind[jav]),$
              2*(!PI)-ABS(winkind[0]-winkind[jav])])

            ; minimum
            FOR iav = 0,N_ELEMENTS(winkind) - 1 DO BEGIN
              sum = 0.0
              FOR jav = 0, N_ELEMENTS(winkind) - 1 DO sum += $
                sumhist[jav,0,cv,isp,n] * MIN([(ABS(winkind[iav]-winkind[jav])),$
                (2*(!PI)-ABS(winkind[iav]-winkind[jav]))])
              IF avsum GT sum THEN BEGIN
                avsum = sum
                avang = winkind[iav]
              ENDIF
            ENDFOR
;            avsum = ASIN(TOTAL(SIN(sumhist[*,0,cv,isp]))/((*i).nbins+1))
            PLOTS, avang, maxz, COLOR=2
            PLOTS, avang, 0, COLOR=2, /CONTINUE
          ENDELSE

          set_output, diag, sp, header=['phase/ky=',$
            rm0es((*series.ky)[(*i).kyind])], commentlines='cross phase: '+$
            get_var_string(var1)+' x '+ get_var_string(var2), $
            dat=[[winkind],[sumhist[*,*,cv,isp,n]]], append=(cv NE 0)
        ENDFOR ; cv

        IF nkyind GT 1 THEN BEGIN
          LOADCT, 33, FILE='internal/colortable.tbl'

          FOR cv = page * cpp, ((page + 1) * cpp - 1) < (n_corrs - 1) DO BEGIN   
            var1 = (*i).corr[0,cv]
            var2 = (*i).corr[1,cv]

            lev_col = contour_levels(sumhist[*,*,cv,isp,n],20)

            IF NOT (*i).log THEN BEGIN
              CONTOUR, sumhist[*,*,cv,isp,n], winkind, (*series.ky)[(*i).kyind], $
                ZVALUE=1.0, /T3D, /NOCLIP, LEVELS=lev_col[*,0], $
                C_COLOR=lev_col[*,1], COLOR=0, XTICKLEN=1.0, YTICKLEN=1.0, $
                XRANGE=[-!PI,!pi], /XSTYLE, /YSTYLE, $
                TITLE=get_var_string(var1,/fancy)+' x '+get_var_string(var2,/fancy), $
                /FILL, YTITLE='!6k!Dy!N'+rho_str, XTITLE='!7a!6', XTICKS=4, $
                XTICKV=[-!PI,-!PI/2.0,0.0,!PI/2,!PI], $
                XTICKNAME=['-!7p!6',' ','0',' ','!7p!6']

              FOR xind = 1, 11 DO BEGIN
                xval = - !PI + xind / 6.0 * !PI
                PLOTS, [xval,xval], !Y.CRANGE, COLOR=125, LINE=1
              ENDFOR
              PLOTS, [0,0], !Y.CRANGE, COLOR=125, LINE=0
            ENDIF ELSE BEGIN
              CONTOUR, sumhist[*,*,cv,isp,n], winkind, $
                (*series.ky)[(*i).kyind], /YLOG, LEVELS=lev_col[*,0], $
                C_COLOR=lev_col[*,1], YRANGE=yrange, COLOR=0, /YSTYLE, $
                XTICKLEN=0, YTICKLEN=-0.02, XRANGE=[-!PI,!pi], $
                TITLE=get_var_string(var1,/fancy)+' x '+get_var_string(var2,/fancy), $
                /FILL, YTITLE='!6k!Dy!N'+rho_str, XTITLE='!7a!6', XTICKS=4, $
                XTICKV=[-!PI,-!PI/2.0, 0.0, !PI/2, !PI], $
                XTICKNAME=['-!7p!6', ' ', '0', ' ', '!7p!6']

              FOR xind = 1, 11 DO BEGIN
                xval = - !PI + xind / 6.0 * !PI
                PLOTS, [xval,xval], 10^!Y.CRANGE, COLOR=125, LINE=1
              ENDFOR
              PLOTS, [0,0], 10^!Y.CRANGE, COLOR=125, LINE=0
            ENDELSE
          ENDFOR ; cv
        ENDIF ; nkyind GT 1

        plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
      ENDFOR ; page
    ENDFOR ; steps

    set_output, diag, sp, /reset
  ENDFOR ; --- species loop

END
