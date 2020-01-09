;This library is dedicated for synthetic diagnostics
;for the DIII-D tokamak. The filter routines etc. 
;are based on source code being kindly provided 
;by A.E. White at APS 2010

;#############################################################################
PRO BES_filter, GENE_data, GENE_R, GENE_Z, R00=R00, Z00=Z00, $
                BES_data=BES_data, BES_xloc=BES_xloc, $
                BES_zloc=BES_zloc, get_synth_BES_pos=get_synth_BES_pos, $
                filter=filter, raw_sig=raw_sig, plot=plot
  COMMON global_vars

  IF N_ELEMENTS(R00) LT 1 THEN BEGIN
     R00 = 1.6995E+00 + 0.5*0.61747
     ;R00 = major_radius of bes channel 9 is at middle of shot
     
     R00 = ((*series.geom).R[par.nz0/2])
     ;[Holland et al., PoP 16, 052301 (2009)]:
     ;synthetic BES array centered in the middle of each simulation
  ENDIF
  IF N_ELEMENTS(Z00) LT 1 THEN $
     Z00 = -0.04 ; BES is 4 cm below midplane for 128913

  NR_BES = 5
  NZ_BES = 6
  N_BES  = NR_BES*NZ_BES

  dx = 0.009 ; 0.9cm spacing
  dz = 0.012 ; 1.2cm spacing

  ;taking M. Shafer et al., Rev. Sci. Instrum 77, 10F110 (2006)
  IF N_ELEMENTS(rho) LT 1 THEN rho = par.x0
  exptrho = [0.5,0.75]
  nextrho = min(abs(rho-exptrho),minidx)
  CASE minidx OF
     0 : BES_lr = 0.022  ;rad. FWHM estimated from Fig. 4 at R~202.5 cm
     1 : BES_lr = 0.0115 ;stated in the text for R~220 cm
  ENDCASE

  ; proper translation of FWHM and normalization for Gaussian 
  BES_lr /= (2.0*sqrt(2.0*ALOG(2.0)))
  BES_ltheta = 0.013/(2.0*sqrt(2.0*ALOG(2.0))) ;1.3cm FWHM in  poloidal direction

  ;radial response function (slightly shifted and asymmetric Gaussian)
  BES_rasymm = 0.1 ;0.46 ;inner FWHM slightly less
  BES_rshift = 0.004 ;rad. response function slightly shifted

  BES_xloc = FLTARR(NR_BES, NZ_BES)
  BES_zloc = FLTARR(NR_BES, NZ_BES)
  
  FOR ir = 0, NR_BES - 1 DO FOR iz = 0, NZ_BES - 1 DO BEGIN
     BES_xloc[ir,iz] = R00 + (ir-NR_BES/2) * dx
     BES_zloc[ir,iz] = Z00 - iz * dz
  ENDFOR

  BES_xloc += BES_rshift

  IF N_ELEMENTS(GENE_data) LE 1 THEN RETURN

  ;create local high res. rectangular grid
  hr_lr = dx*NR_BES+2*(3.0*2.0*MAX([BES_rasymm,1-BES_rasymm])*BES_lr>dx)
  hr_lz = dz*(NZ_BES-1)+2*(3.0*BES_ltheta>dz)
  nind1 = 128
  nind2 = 128
  hr_r = FLTARR(nind1,/NOZERO)
  hr_z = FLTARR(nind2,/NOZERO)
  
  hr_r = REBIN(R00+(-0.5+FINDGEN(nind1)/(nind1-1))*hr_lr,[nind1,nind2])
  hr_z = REBIN(REFORM(Z00+(3.0*BES_ltheta>dz)-(FINDGEN(nind2)/(nind2-1))*$
           hr_lz,[1,nind2]),[nind1,nind2])

  TRIANGULATE, GENE_R, GENE_Z, tria, TOLERANCE=0
  hr_dat = TRIGRID(GENE_R,GENE_Z,GENE_data,tria,NX=nind1,NY=nind2,$
                   XOUT=hr_r[*,0],YOUT=hr_z[0,*])

  BES_xloc2 = FLTARR(NR_BES, NZ_BES)
  BES_zloc2 = FLTARR(NR_BES, NZ_BES)

  BES_hr_idx = LONARR(NR_BES, NZ_BES)

  BES_spot = FLTARR(nind1,nind2,NR_BES,NZ_BES)
  BES_data  = FLTARR(NR_BES,NZ_BES)


  ;create BES array
  FOR ir = 0, NR_BES - 1 DO FOR iz = 0, NZ_BES - 1 DO BEGIN
     mind = MIN((hr_r-BES_xloc[ir,iz])^2+(hr_z-BES_zloc[ir,iz])^2,$
                minidx)
     BES_hr_idx[ir,iz] = minidx
     BES_xloc2[ir,iz] = hr_R[minidx]
     BES_zloc2[ir,iz] = hr_Z[minidx]

     FOR ind2 = 0, nind2-1 DO BEGIN
        FOR ind1 = 0, nind1-1 DO BEGIN
           rdiff = (hr_r[ind1,ind2]-BES_xloc[ir,iz])
           rdiff /= 2.0*((rdiff LT 0)?BES_rasymm:(1.0-BES_rasymm))
           zdiff = (hr_z[ind1,ind2]-BES_zloc[ir,iz])
           BES_spot[ind1,ind2,ir,iz] = EXP(-0.5*((rdiff/BES_lr)^2+$
                 (zdiff/BES_ltheta)^2)) ;improve: map theta to r/z
        ENDFOR
     ENDFOR

;    circular spots
;     BES_spot[*,*,ir,iz] = EXP(-0.5*(((hr_r-BES_xloc[ir,iz])/BES_lr*)^2+$
;         ((hr_z-BES_zloc[ir,iz])/BES_ltheta)^2)) ;improve: map theta to r/z
     BES_spot_norm = TOTAL(BES_spot[*,*,ir,iz])
     
     IF ABS(BES_spot_norm) GT 0 THEN BEGIN
        BES_spot[*,*,ir,iz] /= BES_spot_norm
        BES_data[ir,iz] = TOTAL(hr_dat*BES_spot[*,*,ir,iz])
     ENDIF ELSE print, 'found empty BES spots'
  ENDFOR

  raw_sig = hr_dat[BES_hr_idx]

  IF KEYWORD_SET(filter) THEN $
     GENE_data[*] *= TOTAL(TOTAL(BES_spot,4),3)/(N_BES)

  IF KEYWORD_SET(plot) THEN BEGIN
     psym = 3 ;dot
     psym = 7 ;X
     lthick = 1.5 ;contour lines
     lev_col = contour_levels(GENE_data[*,*],64)
     xmin = MIN(hr_r,max=xmax)
     zmin = MIN(hr_z,max=zmax)
     
     CONTOUR, GENE_data[*,*],GENE_R, GENE_Z,$
              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
              /FILL, /ISOTROPIC, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
              YRANGE=[zmin,zmax], POSITION=[0.0,0.5,0.5,1.0]
     OPLOT, GENE_R, GENE_Z, PSYM=3, COLOR=1
     
     CONTOUR, hr_dat,hr_r,hr_z,$
              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
              /FILL, /ISOTROPIC, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
              YRANGE=[zmin,zmax], POSITION=[0.5,0.5,1.0,1.0],/NOERASE
     CONTOUR, BES_spot[*,*,0,0],hr_r,hr_z,c_colors=[123,1,254],$
       levels=[0.1,0.5,0.9]*MAX(BES_spot[*,*,0,0]),/OVERPLOT, THICK=lthick
     CONTOUR, BES_spot[*,*,0,1],hr_r,hr_z,c_colors=[123,1,254],$
       levels=[0.1,0.5,0.9]*MAX(BES_spot[*,*,0,1]),/OVERPLOT, THICK=lthick
     CONTOUR, BES_spot[*,*,NR_BES-1, NZ_BES-1],hr_r,hr_z,c_colors=[123,1,254],$
       levels=[0.1,0.5,0.9]*MAX(BES_spot[*,*,NR_BES-1,NZ_BES-1]),/OVERPLOT, THICK=lthick
     OPLOT, BES_xloc, BES_zloc, PSYM=psym, COLOR=1, SYMSIZE=0.5
     OPLOT, BES_xloc2, BES_zloc2, PSYM=3, COLOR=1
     
     CONTOUR, raw_sig,BES_xloc,BES_zloc,$
              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
              /FILL, /ISOTROPIC, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
              YRANGE=[zmin,zmax], POSITION=[0.0,0.0,0.5,0.5],/NOERASE
     OPLOT, BES_xloc, BES_zloc, PSYM=psym, COLOR=1, SYMSIZE=0.5
     
     CONTOUR, BES_data,BES_xloc,BES_zloc,$
              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
              /FILL, /ISOTROPIC, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
              YRANGE=[zmin,zmax], POSITION=[0.5,0.0,1.0,0.5],/NOERASE
     OPLOT, BES_xloc, BES_zloc, PSYM=psym, COLOR=1, SYMSIZE=0.5
  ENDIF
  
  IF KEYWORD_SET(get_synth_BES_pos) THEN BEGIN
     BES_xloc = BES_xloc2
     BES_zloc = BES_zloc2
  ENDIF

END

;#############################################################################

PRO CECE_filter, GENE_data, GENE_R, GENE_Z, R00=R00, Z00=Z00, $
                CECE_sig=CECE_sig, CECE_xloc=CECE_xloc, $
                CECE_zloc=CECE_zloc, get_synth_CECE_pos=get_synth_CECE_pos,$
                CECE_spot=CECE_spot, filter=filter, plot=plot, $
                raw_sig = raw_sig, rho=rho

  COMMON global_vars

  IF N_ELEMENTS(R00) LT 1 THEN $
     R00 = 1.6995E+00
  IF N_ELEMENTS(Z00) LT 1 THEN $
     Z00 = 0.049 ; CECE is 4.9cm above midplane (Terry, GENRAY)

  nind1 = N_ELEMENTS(GENE_R[*,0])
  nind2 = N_ELEMENTS(GENE_R[0,*])

  ;create N_CECE *pairs* of CECE channels (for statistics)
  ;same locations as BES, but above rather than below midplane
  ;spot = exp(-0.5*((R-R0)^2/Lr^2 + (Z-Z0)^2/Lz^2))
  ;when R-R0=2*Lr, or Z-Z0=2*Lz, spot = 1/e**2
  ;radial 1/e**2 diameter=4*L is 1cm radially, 3.8cm vertically

  IF N_ELEMENTS(rho) LT 1 THEN rho = par.x0
  exptrho = [0.55,0.65,0.70]
  nextrho = min(abs(rho-exptrho),minidx)
;  print, 'using CECE array for rho = ', exptrho[minidx]
  CASE minidx OF
     0 : CECE_r0 = [0.47,0.52,0.57,0.62,0.67] ; for rho = 0.55
     1 : CECE_r0 = [0.47,0.52,0.57,0.62,0.67]+0.15 ; for rho = 0.65
     2 : CECE_r0 = [0.75,0.775,0.8,0.825,0.85] ; for rho = 0.70
  ENDCASE
  CECE_r0 *= 0.59  ;translate from minor r units to m?
  CECE_r0 += R00   ;shift by major R
;print, '1st choice: ',CECE_r0

;  ;symmetric about domain center with 2.5cm step size as above
;  CECE_r0 = ((*series.geom).R[par.nz0/2])+(FINDGEN(5)*0.025-0.05)
;print, '2nd choice: ',CECE_r0
  N_CECE = N_ELEMENTS(CECE_r0)
  CECE_xloc = FLTARR(2*N_CECE)
  CECE_zloc = Z00+FLTARR(2*N_CECE)

  CECE_lr = 0.01/4.
  CECE_lz = 0.038/4. ;(in PoP17, (2010): 0.032/4 ?)

  ;taking C. Holland et al. Phys. Plasmas 16,052301
;  CECE_lr = 0.00416*0.01
;  CECE_lz = 0.014*0.01   ;???? is that correct????


  CECE_dx = 0.005 ; 0.5cm spacing (called spot_sep1 in Anne's code)

  ;create CECE array
  FOR ir = 0, N_CECE - 1 DO BEGIN
     CECE_xloc[2*ir]   = CECE_r0[ir]
     CECE_xloc[2*ir+1] = CECE_r0[ir]+CECE_dx
  ENDFOR
;print, CECE_xloc
  IF N_ELEMENTS(GENE_data) LE 1 THEN RETURN 

  ;create local high res. rectangular grid
  hr_lr = (MAX(CECE_xloc)-MIN(CECE_xloc)+2.0*(3.0*CECE_lr>CECE_dx))>0.10
  hr_lz = 0.10 ;2.0*20.0*CECE_lz
  nind1 = 128
  nind2 = 128
  hr_r = FLTARR(nind1,/NOZERO)
  hr_z = FLTARR(nind2,/NOZERO)
  
  hr_r = (MIN(CECE_xloc)+MAX(CECE_xloc))/2.+REBIN((-0.5+FINDGEN(nind1)/(nind1-1))*hr_lr,[nind1,nind2])
;  hr_z = REBIN(REFORM(Z00+(20.0*CECE_lz)-(FINDGEN(nind2)/(nind2-1))*$
;           hr_lz,[1,nind2]),[nind1,nind2])
  hr_z = REBIN(REFORM(Z00+(0.05)-(FINDGEN(nind2)/(nind2-1))*$
           hr_lz,[1,nind2]),[nind1,nind2])

  TRIANGULATE, GENE_R, GENE_Z, tria, TOLERANCE=0
  hr_dat = TRIGRID(GENE_R,GENE_Z,GENE_data,tria,NX=nind1,NY=nind2,$
                   XOUT=hr_r[*,0],YOUT=hr_z[0,*])

  CECE_xloc2 = FLTARR(2*N_CECE)
  CECE_zloc2 = FLTARR(2*N_CECE)
  CECE_hr_idx = LONARR(2*N_CECE)

  CECE_spot = FLTARR(nind1,nind2,2*N_CECE)
  CECE_spot_norm = FLTARR(nind1,nind2,2*N_CECE)
  CECE_sig  = FLTARR(2*N_CECE)

  FOR ir = 0, 2*N_CECE-1 DO BEGIN
     mind = MIN((hr_r-CECE_xloc[ir])^2+(hr_z-CECE_zloc[ir])^2,minidx)

     CECE_hr_idx[ir] = minidx
     CECE_xloc2[ir] = hr_R[minidx]
     CECE_zloc2[ir] = hr_Z[minidx]
     CECE_spot[*,*,ir] = EXP(-0.5*(((hr_r-CECE_xloc[ir])/CECE_lr)^2+$
                               ((hr_z-CECE_zloc[ir])/CECE_lz)^2))
     CECE_spot_norm = TOTAL(CECE_spot[*,*,ir])

     IF CECE_spot_norm GT 0 THEN BEGIN
        CECE_spot[*,*,ir] /= CECE_spot_norm
        CECE_sig[ir] = TOTAL(hr_dat*CECE_spot[*,*,ir])
     ENDIF ELSE print, 'found empty CECE spots'
  ENDFOR

  raw_sig = hr_dat[CECE_hr_idx]

  IF KEYWORD_SET(filter) THEN GENE_data[*] *= TOTAL(CECE_spot,3)/(2*N_CECE)

  IF KEYWORD_SET(plot) THEN BEGIN
     LOADCT, 33
     lev_col = contour_levels(GENE_data[*,*],64)
     xmin = MIN(hr_r,max=xmax)
     zmin = MIN(hr_z,max=zmax)
     isotrp = 1
     CONTOUR, GENE_data[*,*],GENE_R, GENE_Z,$
              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
              /FILL, ISOTROPIC=isotrp, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
              YRANGE=[zmin,zmax], POSITION=[0.0,0.5,0.5,1.0]
     OPLOT, GENE_R, GENE_Z, PSYM=3, COLOR=1
     
     CONTOUR, hr_dat,hr_r,hr_z,$
              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
              /FILL, ISOTROPIC=isotrp, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
              YRANGE=[zmin,zmax], POSITION=[0.5,0.5,1.0,1.0],/NOERASE
     CONTOUR, CECE_spot[*,*,0],hr_r,hr_z,c_colors=[123,1],$
              levels=[1./EXP(2.0),0.5]*MAX(CECE_spot[*,*,0]),/OVERPLOT, THICK=0.1
     CONTOUR, CECE_spot[*,*,N_CECE-1],hr_r,hr_z,c_colors=[123,1],$
       levels=[1./EXP(2.0),0.5]*MAX(CECE_spot[*,*,N_CECE-1]),/OVERPLOT, THICK=0.1
     CONTOUR, CECE_spot[*,*,N_CECE],hr_r,hr_z,c_colors=[123,1],$
       levels=[1./EXP(2.0),0.5]*MAX(CECE_spot[*,*,N_CECE]),/OVERPLOT, THICK=0.1
     CONTOUR, CECE_spot[*,*,2*N_CECE-1],hr_r,hr_z,c_colors=[123,1],$
       levels=[1./EXP(2.0),0.5]*MAX(CECE_spot[*,*,2*N_CECE-1]),/OVERPLOT, THICK=0.1
     OPLOT, CECE_xloc, CECE_zloc, PSYM=3, COLOR=1
     OPLOT, CECE_xloc2, CECE_zloc2, PSYM=3, COLOR=1
     
;     CONTOUR, raw_sig,CECE_xloc,CECE_zloc,$
;              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
;              /FILL, /ISOTROPIC, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
;              YRANGE=[zmin,zmax], POSITION=[0.0,0.0,0.5,0.5],/NOERASE
;     OPLOT, CECE_xloc, CECE_zloc, PSYM=3, COLOR=1
     
;     CONTOUR, CECE_sig,CECE_xloc,CECE_zloc,$
;              C_COLORS=lev_col[*,1], LEVELS=lev_col[*,0], $
;              /FILL, /ISOTROPIC, XSTYLE=1, YSTYLE=1,XRANGE=[xmin,xmax],$
;              YRANGE=[zmin,zmax], POSITION=[0.5,0.0,1.0,0.5],/NOERASE
;     OPLOT, CECE_xloc, CECE_zloc, PSYM=3, COLOR=1
  ENDIF
  
  IF KEYWORD_SET(get_synth_CECE_pos) THEN BEGIN
     CECE_xloc = CECE_xloc2
     CECE_zloc = CECE_zloc2
  ENDIF

END

;#############################################################################

PRO REFL_filter, GENE_data, GENE_R, GENE_Z, R00=R00, Z00=Z00, $
                REFL_sig=REFL_sig, REFL_xloc=REFL_xloc, $
                REFL_zloc=REFL_zloc, get_synth_REFL_pos=get_synth_REFL_pos,$
                REFL_spot=REFL_spot, filter=filter, plot=plot,$
                raw_sig = raw_sig

  COMMON global_vars

  IF N_ELEMENTS(R00) LT 1 THEN $
     R00 = 1.6995E+00
  IF N_ELEMENTS(Z00) LT 1 THEN $
     Z00 = 0.052 ; original comments: REFL is at same position as CECE;
                 ; approximately; zloc given by Terry,
                 ; GENRAY is 5.2cm above
                 ; PoP 17, 056103: z_cece-z_refl = 0.3cm
                 ; so it's probably okay (tbg)

  nind1 = N_ELEMENTS(GENE_R[*,0])
  nind2 = N_ELEMENTS(GENE_R[0,*])

  ;create N_REFL *pairs* of REFL channels (for statistics)
  ;same locations as BES, but above rather than below midplane
  ;spot = exp(-0.5*((R-R0)^2/Lr^2 + (Z-Z0)^2/Lz^2))
  ;when R-R0=2*Lr, or Z-Z0=2*Lz, spot = 1/e**2
  ;radial 1/e**2 diameter=4*L is XX cm radially, YY cm vertically

  IF N_ELEMENTS(rho) LT 1 THEN rho = par.x0
  exptrho = [0.55,0.65,0.70] ;???? same as for CECE????????
  nextrho = min(abs(rho-exptrho),minidx)
;  print, 'using CECE array for rho = ', exptrho[minidx]
  CASE minidx OF
     0 : REFL_r0 = [0.47,0.52,0.57,0.62,0.67] ; for rho = 0.55
     1 : REFL_r0 = [0.47,0.52,0.57,0.62,0.67]+0.15 ; for rho = 0.65
     2 : REFL_r0 = [0.75,0.775,0.8,0.825,0.85] ; for rho = 0.70
  ENDCASE
  REFL_r0 *= 0.59  ;translate from minor r units to m?
  REFL_r0 += R00   ;shift by major R

  N_REFL = N_ELEMENTS(REFL_r0)
  REFL_xloc = FLTARR(2*N_REFL)
  REFL_zloc = Z00+FLTARR(2*N_REFL)

  REFL_lr = 0.01/4.
  REFL_lz = 0.035/4.

  REFL_dx = 0.005 ; 0.5cm spacing (called spot_sep1 in Anne's code)

  ;create REFL array
  FOR ir = 0, N_REFL - 1 DO BEGIN
     REFL_xloc[2*ir]   = REFL_r0[ir]
     REFL_xloc[2*ir+1] = REFL_r0[ir]+REFL_dx
  ENDFOR
  IF N_ELEMENTS(GENE_data) LE 1 THEN RETURN 

  ;create local high res. rectangular grid
  hr_lr = (MAX(REFL_xloc)-MIN(REFL_xloc)+2.0*(3.0*REFL_lr>REFL_dx))>0.10
  hr_lz = 0.10 ;2.0*20.0*CECE_lz
  nind1 = 128
  nind2 = 128
  hr_r = FLTARR(nind1,/NOZERO)
  hr_z = FLTARR(nind2,/NOZERO)
  
  hr_r = (MIN(REFL_xloc)+MAX(REFL_xloc))/2.+REBIN((-0.5+FINDGEN(nind1)/(nind1-1))*hr_lr,[nind1,nind2])
;  hr_z = REBIN(REFORM(Z00+(20.0*CECE_lz)-(FINDGEN(nind2)/(nind2-1))*$
;           hr_lz,[1,nind2]),[nind1,nind2])
  hr_z = REBIN(REFORM(Z00+(0.05)-(FINDGEN(nind2)/(nind2-1))*$
           hr_lz,[1,nind2]),[nind1,nind2])

  TRIANGULATE, GENE_R, GENE_Z, tria, TOLERANCE=0
  hr_dat = TRIGRID(GENE_R,GENE_Z,GENE_data,tria,NX=nind1,NY=nind2,$
                   XOUT=hr_r[*,0],YOUT=hr_z[0,*])

  REFL_xloc2 = FLTARR(2*N_REFL)
  REFL_zloc2 = FLTARR(2*N_REFL)
  REFL_hr_idx = LONARR(2*N_REFL)

  REFL_spot = FLTARR(nind1,nind2,2*N_REFL)
  REFL_spot_norm = FLTARR(nind1,nind2,2*N_REFL)
  REFL_sig  = FLTARR(2*N_REFL)

  FOR ir = 0, 2*N_REFL-1 DO BEGIN     
     mind = MIN((hr_r-REFL_xloc[ir])^2+(hr_z-REFL_zloc[ir])^2,minidx)

     REFL_hr_idx[ir] = minidx
     REFL_xloc2[ir] = hr_R[minidx]
     REFL_zloc2[ir] = hr_Z[minidx]
     REFL_spot[*,*,ir] = EXP(-0.5*(((hr_r-REFL_xloc[ir])/REFL_lr)^2+$
                               ((hr_z-REFL_zloc[ir])/REFL_lz)^2))
     REFL_spot_norm = TOTAL(REFL_spot[*,*,ir])

     IF REFL_spot_norm GT 0 THEN BEGIN
        REFL_spot[*,*,ir] /= REFL_spot_norm
        REFL_sig[ir] = TOTAL(hr_dat*REFL_spot[*,*,ir])
     ENDIF ELSE print, 'found empty REFL spots'
  ENDFOR

  raw_sig = hr_dat[REFL_hr_idx]

  IF KEYWORD_SET(filter) THEN GENE_data[*] *= TOTAL(REFL_spot,3)/(2*N_REFL)

  IF KEYWORD_SET(plot) THEN printerror, 'plot not implemented yet in REFL'

  IF KEYWORD_SET(get_synth_REFL_pos) THEN BEGIN
     REFL_xloc = REFL_xloc2
     REFL_zloc = REFL_zloc2
  ENDIF

END

;#############################################################################

PRO plot_D3D_separatrix, iterdbfile=iterdbfile, oplot=oplot

  IF NOT KEYWORD_SET(iterdbfile) THEN $
     iterdb_file = GETENV('HOME')+'/gene11/profdata/'+$
                   'iterdb.138038_APS2009_Lmode_FINAL'
  
  Rbound=read_iterdb_1d_data(iterdb_file,'r points for plasma boundary',$
                             nr_var='nplasbdry')
  Zbound=read_iterdb_1d_data(iterdb_file,'z points for plasma boundary',$
                             nr_var='nplasbdry')
  
  IF NOT KEYWORD_SET(oplot) THEN $
     PLOT, 0.01*Rbound, 0.01*Zbound,$
           COLOR=1,/XSTYLE, /YSTYLE, /ISOTROPIC,$
           XTITLE='!6R / m', YTITLE='!6Z / m' $
  ELSE OPLOT, 0.01*Rbound, 0.01*Zbound,COLOR=1
  
END


;#############################################################################



FUNCTION calc_rfspect1d, x, y, NFFT=NFFT, FREQ=freq, NO_HANNING=no_hanning, $
  ERR=err, COHERENCY=coherency, PHASE=PHASE, APHAS_ERR=phase_err, ANGLE1=angle1

; x (y) : time traces of signal(s)
; NFFT  : number of frequencies in FFT; the time trace will be split
;         and Hanning-filtered if > NFFT elements are available
; FREQ  : returns the dimensionless array of (positive) frequencies
; NO_HANNING : enforce running without Hanning filters
; ERR   : error estimate for auto/cross spectrum
; COHERENCY : coherency of the signals
; PHASE : phase between the signals
; APHASE_ERR : associated error

;
; T Goerler 08-02-2013: rerenamed calc_rfspect1d_phase_one to
; calc_rfspect1d and switched inverse FFT back to FFT
;
; AE White 07-13-09
; hunt down any differences
; in the SIGN of the phase angle that may be showing up
; Analysis must be checked because measurements of 
; lemo polarity at digitizer showed no switching in hardware.
;
; AE White 04-03-08
; Added the cross phase output for CECE REFL studies
;
; C. Holland, UCSD
; v1.0: 2-6-2008
;
;
; Calculates magnitude of cross-power Sxy between two signals x and y, 
; assumed to be real, which is properly normalized such that if x=y,
; MEAN(x) = 0, NFFT=NX, and no hanning windo is used, then the
; mean-square value of x, 
;
; x_msv = (sum_i=0,NX-1 x_i^2)/NX 
;
; is equal to sum_i=0,NX/2 Sxy_i
;
; Thus if x and y were voltages sampled every second, the units of Sxy
; would be volt**2/sec
;
; 3-13-2008: fixed bug in when abs. value of S taken
;

  NX = N_ELEMENTS(x)
  IF N_ELEMENTS(y) EQ 0 THEN y = x
  IF (N_ELEMENTS(y) NE N_ELEMENTS(x)) THEN BEGIN
     PRINT, 'x and y not same size, returning 0'
     freq = FINDGEN(NX)
     RETURN, FLTARR(NX)
  ENDIF
  IF NOT KEYWORD_SET(NFFT) THEN NFFT=NX

  freq = FINDGEN(NFFT/2+1)
  S = COMPLEXARR(NFFT)
  Px = FLTARR(NFFT)
  Py = FLTARR(NFFT)
  NM = NX/NFFT  ;number of Hanning windows
  ;print, 'number of Hannig windows in calc_rfspect1d:', NM

  IF KEYWORD_SET(no_hanning) THEN hanW = FLTARR(NFFT) + 1 $
  ELSE hanW = HANNING(NFFT)
  hw_norm = NFFT/TOTAL(hanW^2)
  FOR im = 0, NM-1 DO BEGIN
     idx = im*NFFT + LINDGEN(NFFT)
     xhat = (x[idx]-MEAN(x[idx]))*hanW
     yhat = (y[idx]-MEAN(y[idx]))*hanW
     xk = FFT(xhat)  ;, -1) !REVERT?
     yk = FFT(yhat)  ;, -1) !REVERT?
     S += CONJ(xk)*yk
     
     Px += ABS(xk)^2
     Py += ABS(yk)^2
  ENDFOR
  ;Bendat and Piersol Random Data 2000 pp 445 Eqn. 11.149
  ;Gxy_k= 2/(N delta t) x [X_k* x Y_k] k= 0...N/2 
  ;Gxy=  S[0:NFFT/2]
  Gxy=S

  ;check this:
  angle1 = ATAN(Gxy, /phase)
  phase = ATAN(Gxy[0:NFFT/2], /phase)
  ;phase = angle1
  ;;angle1 = angle1[0:NFFT/2]*hw_norm/NM
  ;;phase = phase[0:NFFT/2]*hw_norm/NM    

  S = ABS(S[0:NFFT/2])*hw_norm/NM
  S[1:NFFT/2-1] *= 2.0
  
  Px = Px[0:NFFT/2]*hw_norm/NM
  Py = Py[0:NFFT/2]*hw_norm/NM
  Px[1:NFFT/2-1] *= 2.0
  Py[1:NFFT/2-1] *= 2.0
  coherency = S/SQRT(Px*Py)
  err = SQRT(Px*Py*(1. - coherency^2)/NM)
  phase_err = sqrt(1. - coherency^2)/(coherency * 2 * sqrt(NM))
 
  RETURN, S
END ;calc_rfspect1d
