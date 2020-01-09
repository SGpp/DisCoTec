FUNCTION ball_info

  RETURN, {$
    type      : 'mom',$
    title     : 'Ballooning modes',$
    help_text : ['Plots modes in the ballooning representation '+$
                 '(stretching out beyond the parallel box).'],$
    ext_vars  : [['vars','0','variables to be investigated; default: 0'],$
                 ['range','0','relative zoom-in range, e.g., [0.4,0.6]; default: off'],$
                 ['kyind','0','i_ky or [i_ky1,...,i_kyn]; '+$
                  'default: lowest non-zero mode'],$
                 ['n_pod','0','set to nonzero integer to plot '+$
                  'Proper Orthogonal Decomposition data for n_pod '+$
                  'modes; set to -1 for maximal n_pod'],$
                 ['norm','1','switch on/off normalization to zero '+$
                  'ballooning angle which is often used for '+$
                  'benchmarks but which is also erasing phase '+$
                  'information'],$
;                 ['reim','1','plot real and imaginary part in '+$
;                  'addition to the absolute value'],$
                 ['var_ratios','1','show pairwise ratios of plot '+$
                  'variables; must specify even number of vars']]}

END

;######################################################################

PRO ball_init, diag

  COMMON global_vars

  IF ABS(par.shat) LT 1e-5 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + $
      ': shat not finite, use zprofile diagnostic instead'
    (*diag).selected = 0
    RETURN
  ENDIF

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.vars) LT 1 THEN *e.vars = 0
  IF NOT KEYWORD_SET(*e.range) THEN *e.range = [-1,-1]
  IF (*e.range)[0] GE (*e.range)[1] THEN *e.range = [-1,-1]
  IF N_ELEMENTS(*e.kyind) LT 1 THEN $
    *e.kyind = (par.ky0_ind GT 0) ? par.ky0_ind : 1
  IF N_ELEMENTS(*e.n_pod) NE 1 THEN *e.n_pod = 0
  use_pod = *e.n_pod NE 0
  IF NOT KEYWORD_SET(*e.norm) THEN *e.norm = 0
; old input parameter to switch between abs only and abs/re/im;
; now, hard-wired to the latter
;  IF NOT KEYWORD_SET(*e.reim) THEN *e.reim = 0
  reim = 1
  IF NOT KEYWORD_SET(*e.var_ratios) THEN *e.var_ratios = 0

  IF par.ky0_ind GT 0 THEN BEGIN
    IF TOTAL(WHERE(*e.kyind LT par.ky0_ind)) NE -1 THEN BEGIN
      PRINT, (*diag).name + ': skipping analysis of ky=0 or dummy modes'
      IF TOTAL(WHERE(*e.kyind GE par.ky0_ind)) NE -1 THEN $
        *e.kyind = (*e.kyind)[WHERE(*e.kyind GE par.ky0_ind)] $
        ELSE (*diag).selected = 0
    ENDIF
  ENDIF

  IF *e.var_ratios AND ((N_ELEMENTS(*e.vars) MOD 2) NE 0) THEN BEGIN
    PRINT, (*diag).name + ': odd number of variables; not plotting ratios'
    *e.var_ratios = 0
  ENDIF

  Cyq0_x0 = (par.magn_geometry NE 'circular') AND $
     ((par.x0 GT 0.0) AND (par.x0 LT 1e10)) ? $
     (*series.geom).C_y * par.q0 / par.x0 : 1.0


  ; number of meaningful (i.e., connected) x modes
  nxmod = par.x_local ? ((par.nx0 - 1) / 2) * 2 + 1 : par.nx0

  ; store modified index of the radial direction
  jbar = 0
  aux_jbar_fac = 1
  IF nxmod GT 1 THEN BEGIN
    jbar = INTARR(nxmod)
    jbar[0] = par.nx0 / 2 + 1 + INDGEN(nxmod/2)
    jbar[nxmod/2] = INDGEN(nxmod-nxmod/2)
    IF par.shat LT 0 THEN jbar = REVERSE(jbar)
    shear_sign = par.shat LT 0 ? -1 : 1

    nexc = ROUND(par.kymin*ABS(par.shat)*par.lx*Cyq0_x0)

    IF reim THEN BEGIN
      twopii = 2.0 * !PI * COMPLEX(0.0,1.0)
      nky = N_ELEMENTS(*e.kyind)
      shift = nexc * FIX(*e.kyind)
      phasefac = (-1)^shift * EXP(-twopii*par.n0_global*par.q0*(*e.kyind))

      aux_jbar_fac = COMPLEXARR(nxmod,nky,par.nz0,/NOZERO)
      FOR k = 0, par.nz0 -1 DO BEGIN
        FOR j = 0, nky - 1 DO BEGIN
          aux_jbar_fac[nxmod/2+1:nxmod-1,j,k] = shear_sign * $
            phasefac[j]^[jbar[nxmod/2+1:nxmod-1]/shift[j]]
          aux_jbar_fac[0:nxmod/2-1,j,k] = shear_sign * $
            CONJ(phasefac[j])^ABS([(jbar[0:nxmod/2-1]-par.nx0)/shift[j]])
          aux_jbar_fac[nxmod/2,j,k] = 1
        ENDFOR
      ENDFOR
    ENDIF
  ENDIF

  i = set_internal_vars(diag,{$
    n_vars        : N_ELEMENTS(*e.vars),$
    vars          : *e.vars,$
    range         : *e.range,$
    kyind         : *e.kyind,$
    norm          : *e.norm,$
    reim          : reim,$
    var_ratios    : *e.var_ratios,$
    use_pod       : use_pod,$
    n_pod         : *e.n_pod,$
    varb_id       : PTR_NEW(),$
    nxmod         : nxmod,$
    jbar          : jbar,$
    Cyq0_x0       : Cyq0_x0,$
    aux_jbar_fac  : aux_jbar_fac})

  IF par.x_local THEN fft_format, kxky=(*e.vars) $
    ELSE fft_format, sxky=(*e.vars)

END

;######################################################################

PRO ball_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nky = N_ELEMENTS((*i).kyind)

  varb_temp = COMPLEXARR((*i).nxmod*par.nz0,nky,(*i).n_vars,$
    gui.out.n_spec_sel,/NOZERO)
  IF gui.out.n_spec_sel LT 2 THEN varb_temp = $
    REFORM(varb_temp,[(*i).nxmod*par.nz0,nky,(*i).n_vars,$
    gui.out.n_spec_sel],/OVERWRITE)

  yind_arr = REBIN(REFORM(INDGEN(nky),[1,nky,1]),$
    [(*i).nxmod,nky,par.nz0])
  xzind_arr = REBIN(REFORM(INDGEN(par.nz0),[1,1,par.nz0]),$
    [(*i).nxmod,nky,par.nz0])+REBIN(INDGEN((*i).nxmod)*par.nz0,$
    [(*i).nxmod,nky,par.nz0])

  IF par.x_local THEN BEGIN
    FOR isp = 0, gui.out.n_spec_sel - 1 DO FOR v = 0, (*i).n_vars - 1 DO $
      varb_temp[xzind_arr,yind_arr,REBIN([v],[(*i).nxmod,nky,par.nz0]),$
      REBIN([isp],[(*i).nxmod,nky,par.nz0])] = (*i).aux_jbar_fac * $
      (*mom[isp,(*i).vars[v]].kxky)[(*i).jbar,(*i).kyind,*]
  ENDIF ELSE BEGIN
    IF (par.n0_global NE 0) THEN printerror, $
       'n0_global not implemented in global ball diag'
    nx0_inv = 1.0 / par.nx0
    xvar = - 0.5 * par.lx + INDGEN(par.nx0) * par.dx
    FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
      FOR v = 0, (*i).n_vars - 1 DO BEGIN
        FOR modeind = 0, nky - 1 DO BEGIN ; -- start mode loop
          FOR j = 0, (*i).nxmod - 1 DO BEGIN
            exp_fac = EXP(-2.0D0*!DPI/par.lx*COMPLEX(0.0,1.0)*(*i).jbar[j]*xvar)
            FOR k = 0, par.nz0 - 1 DO varb_temp[k+j*par.nz0,modeind,isp] = $
              nx0_inv * TOTAL((*mom[isp,(*i).vars[v]].sxky)$
              [INDGEN((*i).nxmod),(*i).kyind[modeind],k]*exp_fac)
          ENDFOR
        ENDFOR ; --- end mode loop
      ENDFOR ; --- end vars loop
    ENDFOR ; --- end species loop
  ENDELSE

  (*i).varb_id = time_avg((*i).varb_id,varb_temp,mom_time,$
    fwd=(gui.out.res_steps OR (*i).use_pod))

END

;######################################################################

PRO ball_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  nky = N_ELEMENTS((*i).kyind)

  range_detail_lo = (*i).range[0]
  range_detail_hi = (*i).range[1]

  varb = time_avg((*i).varb_id,/avg,$
    fwd=(gui.out.res_steps OR (*i).use_pod),tarr=time)

  IF (*i).use_pod THEN BEGIN
    IF series.step_count LT 2 THEN BEGIN
      printerror, (*diag).name + $
        ' error: need multiple time steps for POD analysis'
      RETURN
    ENDIF

    IF (*i).n_pod GT series.step_count THEN BEGIN
      (*i).n_pod = series.step_count

      PRINT, (*diag).name + ' warning: reducing n_pod to ' + $
        rm0es((*i).n_pod) + ' (number of time steps)'
    ENDIF

    IF (*i).n_pod EQ -1 THEN (*i).n_pod = series.step_count
  ENDIF

  n_varplots = (*i).var_ratios ? (*i).n_vars / 2 : (*i).n_vars

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN ; --- species loop
    sp = (*gui.out.spec_select)[isp]

    set_output, diag, sp, /ps, multi=[0,1,2]

    FOR n = 0, (gui.out.res_steps AND NOT (*i).use_pod) * $
      (series.step_count - 1) DO BEGIN

      IF gui.out.res_steps AND NOT (*i).use_pod THEN BEGIN
        IF time[n] GE 0.0 THEN BEGIN
          timeinfo = '!C t = ' + rm0es(time[n]) + ' ' + $
            get_var_string(1,/time,/ounit)
        ENDIF ELSE BEGIN
          gamma_omega = rm0es(get_eigenvalue(-time[n]),prec=3)
          gamom_str = '(' + gamma_omega[0] + ',' + gamma_omega[1] + ')'
          timeinfo = '!C EV ' + rm0es(-time[n]) + ', (!7c,x!6) = ' + $
            gamom_str
        ENDELSE
      ENDIF ELSE timeinfo = '!C t=' + rm0es(gui.out.start_t) + '-' + $
        rm0es(gui.out.end_t) + get_var_string(1,/time,/ounit)

      FOR v = 0, n_varplots - 1 DO BEGIN
        var_string = (*i).var_ratios ? $
          '(' + get_var_string((*i).vars[v]) + '/' + $
          get_var_string((*i).vars[v+1]) + ')' : $
          get_var_string((*i).vars[v])

        FOR modeind = 0, nky - 1 DO BEGIN ; --- y mode loop
          IF par.adapt_lx THEN bigN = 1 ELSE bigN = ROUND(ABS(par.shat)*$
            series.lx*(*series.ky)[(*i).kyind[modeind]]*(*i).Cyq0_x0) * par.n_pol

          ; plot indices: select 2/4/... pi ranges to plot
          ; pinddim: number of points in the final curve
          pinddim = ((((par.nx0 - 1) / 2) / bigN) * 2 + 1) * par.nz0
          ; pindrange: number of wing 2/4/... pi ranges per sign for curve
          pindrange = (pinddim / par.nz0 - 1) / 2
          ; origrange: number of wing 2/4/... pi ranges per sign in phib
          origrange = (par.nx0 - 1) / 2
          ; pind: x axis plot index
          pind = INTARR(pinddim)
          ; reloffset: first 2/4/... pi range on wing that is relevant
          reloffset = ROUND((origrange/FLOAT(bigN)-$
            FLOAT(FIX(origrange/FLOAT(bigN)+0.001)))*bigN)

          FOR p = 0, 2 * pindrange DO BEGIN
            pind[p*par.nz0+INDGEN(par.nz0)] = $
              (reloffset + p * bigN) * par.nz0 + INDGEN(par.nz0)
          ENDFOR

          ; x axis plot variable
          tppr = (2 * FINDGEN(pinddim) - pinddim) / par.nz0 * par.n_pol

          kyreal = (*series.ky)[(*i).kyind[modeind]]
          kyreal = ROUND(1e3*kyreal) * 1e-3

          fulltitle_base = '!N!6k!Dy!N = ' + rm0es(kyreal)
          var_string = (*i).var_ratios ? $
            '!N!6(' + get_var_string((*i).vars[v],/fancy) + '/' + $
            get_var_string((*i).vars[v+1],/fancy) + ')' : $
            get_var_string((*i).vars[v],/fancy)
          var_string2 = (*i).var_ratios ? $
            '!N!6(' + get_var_string((*i).vars[v]) + '/' + $
            get_var_string((*i).vars[v+1]) + ')' : $
            get_var_string((*i).vars[v])
          xtitle = '!N!7h!6!Dp!N / !7p!6'
          ytitle = '!N!9#' + var_string + '!D!6B!N!9#!6'

          IF NOT (*i).use_pod THEN BEGIN
            plotdata = (*i).var_ratios ? $
              varb[pind,modeind,v,isp,n] / varb[pind,modeind,v+1,isp,n] : $
              varb[pind,modeind,v,isp,n]
          ENDIF ELSE BEGIN
            pod_data_in = (*i).var_ratios ? $
              REFORM(varb[pind,modeind,v,isp,*]/varb[pind,modeind,v+1,isp,*],$
              [pinddim,series.step_count]) : $
              REFORM(varb[pind,modeind,v,isp,*],[pinddim,series.step_count])

            LA_SVD, pod_data_in, sing_vals_t, U_matrix_tt, V_matrix_tz

;            plotdata_ext = (TRANSPOSE(TEMPORARY(V_matrix_tz)))[*,0:(*i).n_pod-1]
            plotdata_ext = TRANSPOSE((TEMPORARY(V_matrix_tz))[0:(*i).n_pod-1,*])

            IF (*i).n_pod GT 1 THEN BEGIN ; plot singular values
              PLOT, INDGEN((*i).n_pod) + 1, sing_vals_t, COLOR=1, $
                /XSTYLE, TITLE=fulltitle_base, XTITLE='!6POD #', $
                YTITLE='!6S!DPOD!N('+var_string+'!D!6B!N)'
              FOR p = 0, (*i).n_pod - 1 DO OPLOT, $
                [p+1,p+1], [1,1] * sing_vals_t[p], COLOR=((p+1)<255), PSYM=4

              tt_min = MIN(ABS(U_matrix_tt),MAX=tt_max)
              PLOT, time[[0,series.step_count-1]], [tt_min,tt_max], COLOR=1, $
                /NODATA, /XSTYLE, /YSTYLE, XTITLE=get_var_string(0,/time,/units), $
                YTITLE=ytitle, TITLE=fulltitle_base
              FOR p = 0, (*i).n_pod - 1 DO OPLOT, $
                time, ABS(U_matrix_tt[p,*]), COLOR=((p+1)<255)

              plot_info_str, diag
            ENDIF
          ENDELSE

          FOR p = 0, ((*i).n_pod > 1) - 1 DO BEGIN 
            IF (*i).use_pod THEN plotdata = plotdata_ext[*,p]

            IF (*i).norm THEN plotdata /= plotdata[pinddim/2] 

            plotmin = MIN(ABS(plotdata),MAX=plotmax) > 1e-6
            absplotmin = plotmin
            IF (*i).reim THEN BEGIN
              plotmin = plotmin < MIN(FLOAT(plotdata),$
                MAX=plotmax1) < MIN(IMAGINARY(plotdata),$
                MAX=plotmax2)
              plotmax = plotmax > plotmax1 > plotmax2
            ENDIF

            fulltitle = fulltitle_base
            IF (*i).use_pod THEN fulltitle += ': !6POD ' + $
              rm0es(p+1) + '/' + rm0es((*i).n_pod)

            FOR log = 0, 1 DO BEGIN
              yrange = [log?(absplotmin>1e-6):plotmin,plotmax]

              !P.MULTI[0] = log

              PLOT, tppr, ABS(plotdata), /NODATA, $
                /XSTYLE, XTITLE=xtitle, XRANGE=[tppr[0],tppr[pinddim-1]], $
                YTITLE=((*i).reim ? '' : ytitle), YRANGE=yrange, $
                TITLE=fulltitle, COLOR=1, YLOG=log

              IF (*i).reim THEN BEGIN ; display colored axis annotations
                AXIS, YAXIS=0, YRANGE=yrange, YLOG=log, $
                  YTITLE='!6,!15Re!6('+var_string+$
                  '!D!6B!N)    ', COLOR=2
                AXIS, YAXIS=0, YRANGE=yrange, YLOG=log, $
                  YTITLE='!6               ,!15Im!6('+$
                  var_string+'!D!6B!N)    ', COLOR=4

                !P.MULTI[0] = log
                PLOT, tppr, ABS(plotdata), /NODATA, $
                  /XSTYLE, XTITLE=xtitle, XRANGE=[tppr[0],tppr[pinddim-1]], $
                  YTITLE=ytitle+((*i).reim ? '               ' : ''), $
                  YRANGE=yrange, TITLE=fulltitle, COLOR=1, YLOG=log, /NOERASE
              ENDIF

              IF (*i).reim THEN BEGIN
                OPLOT, tppr, IMAGINARY(plotdata), COLOR=4
                OPLOT, tppr, FLOAT(plotdata), COLOR=2
                OPLOT, !X.CRANGE, [0,0], COLOR=1, LINE=1
              ENDIF
              OPLOT, tppr, ABS(plotdata), COLOR=1
            ENDFOR ;log=0,1

            IF (*i).range[0] GE 0.0 THEN BEGIN ; plot zoomed region
              !P.MULTI[0] = 0
              detail_range = ROUND(range_detail_lo*(pinddim-1)) + $
                             INDGEN(ROUND(range_detail_hi*(pinddim-1))-$
                                    ROUND(range_detail_lo*(pinddim-1)))
              plotmin = MIN(ABS(plotdata[detail_range]),$
                            MAX=plotmax) > 1e-6
              IF (*i).reim THEN BEGIN
                 plotmin = plotmin < MIN(FLOAT(plotdata[detail_range]),$
                    MAX=plotmax1) < MIN(IMAGINARY(plotdata[detail_range]),$
                    MAX=plotmax2)
                 plotmax = plotmax > plotmax1 > plotmax2
              ENDIF
              rangetitle = '!N!6 range [' + rm0es(range_detail_lo) + $
                ',' + rm0es(range_detail_hi) + '], k!Dy!N = ' + rm0es(kyreal)
              IF (*i).use_pod THEN rangetitle += ': !6POD ' + $
                rm0es(p+1) + '/' + rm0es((*i).n_pod)

              FOR log = 0, 1 DO BEGIN
                PLOT, tppr[detail_range], ABS(plotdata[detail_range]), /NODATA, $
                  /XSTYLE, XTITLE='!N!7h!6!Dp!N / !7p!6', $
                  YRANGE=[log?(plotmin>1e-6):plotmin,plotmax], /YSTYLE, $
                  YTITLE=ytitle, TITLE=rangetitle, COLOR=1, YLOG=log

                IF (*i).reim THEN BEGIN
                  OPLOT, tppr[detail_range], IMAGINARY(plotdata[detail_range]), COLOR=4
                  OPLOT, tppr[detail_range], FLOAT(plotdata[detail_range]), COLOR=2
                  OPLOT, !X.CRANGE, [0,0], COLOR=1, LINE=1
                ENDIF
                OPLOT, tppr[detail_range], ABS(plotdata[detail_range]), COLOR=1
              ENDFOR
            ENDIF

            plot_info_str, diag, $
              time=(gui.out.res_steps ? time[n] : undef)

            IF (*i).use_pod THEN BEGIN
              set_output, diag, sp, header=['theta','Re('+var_string2+')',$
                'Im('+var_string2+')'], commentline='ky='+rm0es(kyreal)+$
                ', POD='+rm0es(p+1)+'/'+rm0es((*i).n_pod)+$
                ' with SV='+rm0es(sing_vals_t[p]), append=(v+p+modeind GT 0)+$
                ((*i).norm?'; *** DATA NORMALIZED to theta=0 value ***':''),$
                dat=[[tppr[*]],[FLOAT(plotdata)],[IMAGINARY(plotdata)]]
            ENDIF

          ENDFOR ; --- POD number loop
        ENDFOR ; --- y mode loop

        IF NOT (*i).use_pod THEN $
          set_output, diag, sp, header=['theta','Re('+var_string2+')',$
          'Im('+var_string2+')'], commentline='ky='+rm0es(kyreal)+$
          ((*i).norm?'; *** DATA NORMALIZED to theta=0 value ***':''),$
          dat=[[tppr[*]],[FLOAT(plotdata)],[IMAGINARY(plotdata)]], append=(v GT 0)
      ENDFOR ; --- var loop
    ENDFOR

    set_output, diag, sp, /reset
  ENDFOR ; --- species loop

END
