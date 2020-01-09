FUNCTION svd_quasilinear_info

  RETURN, {$
    type      : 'extended',$
    title     : 'SVD quasilinear test',$
    help_text : ['Calculates quasilinear expectations from SVD modes.'],$
    ext_vars  : ['mnums','0','number of SVD modes to plot; default: 1' ]}

END

;######################################################################

PRO svd_quasilinear_init, diag

  COMMON global_vars

END

;######################################################################

PRO svd_quasilinear_loop, diag

  COMMON global_vars


    svd_file_string='field'
    file_svd = get_file_path(svd_file_string)
    print,'file_svd',file_svd.path[0]
    nrg_file_string='nrg'
    file_nrg = get_file_path(nrg_file_string)
    print,'file_nrg',file_nrg.path[0]

    ;print,"file_svd",file_svd

    if file_svd.exist eq 0 then begin
      printerror, 'Error: c_tot file(s) not found'
      return
    endif

    last_time = -1.0d * 1e20

    svd_lun = 0l
    nrg_lun = 0l

    ;print,"Here we are in svd_time_loop!"
    ;print,"Reading c_tot_svd file."
    ;print,"Number of pod modes and number of time steps:", (*i).num_time

    if series.n_runs gt 1 then begin
      print,"n_runs must be 1."
      stop 
    endif

    ;time=FLTARR((*i).num_time)
    ;c_tot_in=COMPLEXARR((*i).num_time)
    ;c_modes=COMPLEXARR((*i).mnums,(*i).num_time)

    time_in = par.prec_single ? 0.0 : 0.0D
    nmodes=get_nrg_time(0,/get_n_steps)
    print,"nmodes",nmodes
    print,par.n_spec

    ;time = par.prec_single ? $
    ;  FLTARR((*i).num_time,/NOZERO) : $
    ;  DBLARR((*i).num_time,/NOZERO)
    ;print,"nx0",par.nx0
    ;print,"nky0",par.nky0
    ;print,"nz0",par.nz0
    phi_in = par.prec_single ? $
      COMPLEXARR(par.nx0,par.nz0,/NOZERO) : $
      DCOMPLEXARR(par.nx0,par.nz0,/NOZERO)
    apar_in = par.prec_single ? $
      COMPLEXARR(par.nx0,par.nz0,/NOZERO) : $
      DCOMPLEXARR(par.nx0,par.nz0,/NOZERO)
    phisq_tot = par.prec_single ? $
      FLTARR(nmodes) : $
      DBLARR(nmodes)
    nrg_modes = par.prec_single ? $
      FLTARR(8,par.n_spec,nmodes) : $
      DBLARR(8,par.n_spec,nmodes)
    nrg_in = par.prec_single ? $
      FLTARR(8) : $
      DBLARR(8)
    ;new_nrg_data = REPLICATE($
    ;  {time:0.0,data:FLTARR(8,par.n_spec)},nrg_steps)
    ;print,phisq_tot

    for run = 0, series.n_runs - 1 do begin
      openr, nrg_lun, file_nrg.path[run], /get_lun, error=err
      openr, svd_lun, file_svd.path[run], /get_lun, error=err, $
        /f77_unformatted, swap_endian=series.swap_endian
      if err ne 0 then begin
        printerror, file_svd.path[run] + ' does not exist', /popup
	return
      endif else if !quiet ne 1 then print, 'reading ', $
        file_svd.path[run]

      print,'info time',size(time)
      for itime=0,nmodes-1 do begin
         READU, svd_lun, time_in
         ;print,itime,time_in
         READU, svd_lun, phi_in
         READU, svd_lun, apar_in
         READF, nrg_lun, time_in
         FOR species_num=0,par.n_spec-1 DO BEGIN
           READF, nrg_lun, nrg_in
           nrg_modes[*,species_num,itime]=nrg_in
         ENDFOR
         ;print,time_in,nrg_modes[6,0,itime]

         FOR kxind=0, par.nx0-1 DO BEGIN
           phisq_tot(itime)=phisq_tot(itime)+ABS(TOTAL(phi_in[kxind,*]$
		   *(*series.geom).jacobian[*])/TOTAL((*series.geom).jacobian))^2
         ENDFOR

      endfor
      ;plot, phisq_tot,/ylog
      plot, nrg_modes[6,0,*]
      print,nrg_modes[6,0,0]

    endfor
    free_lun,nrg_lun
    free_lun,svd_lun

    phisq_tot0=TOTAL(phisq_tot)
    phisq_1=phisq_tot[0]
    phisq_res=TOTAL(phisq_tot[1:nmodes-1])
    QiES_tot=TOTAL(nrg_modes[6,0,*])
    QiES_1=nrg_modes[6,0,0]
    QiES_res=TOTAL(nrg_modes[6,0,1:nmodes-1])
    svd_quasi_file = get_file_path('svd_quasilinear')
    print,svd_quasi_file.path[0]
    R1=QiES_1/phisq_1
    R_res=QiES_res/phisq_res
    openw,svd_lun,svd_quasi_file.path[0],/get_lun,error=err
    printf,svd_lun,"Total phi_sq:     ",phisq_tot0
    printf,svd_lun,"phi_sq 1:         ",phisq_1
    printf,svd_lun,"phi_sq res.:      ",phisq_res
    printf,svd_lun,"Total Q_i_ES:     ",QiES_tot
    printf,svd_lun,"Q_i_ES 1:         ",QiES_1
    printf,svd_lun,"Q_i_ES res.:      ",QiES_res
    printf,svd_lun,"R1:               ",QiES_1/phisq_1
    printf,svd_lun,"R_res:            ",QiES_res/phisq_res
    printf,svd_lun,"alpha:            ",phisq_res/phisq_1
    printf,svd_lun,"beta:             ",R_res/R1
    printf,svd_lun,"Q_ql/Q_true       ",R1*(phisq_tot0)/QiES_tot
    free_lun,svd_lun


END

;######################################################################




