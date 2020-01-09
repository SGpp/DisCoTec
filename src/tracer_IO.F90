!>Reads equilibrium data from an EFIT g-eqdsk file.
module tracer_IO
  use discretization
  use lagrange_interpolation
  use spline_interp
  use par_mod, only: pi,psi_o_twopi,q_scalefac,print_ini_msg, lilo
  use file_io, only: get_unit_nr
  use communications
  use mtrandom

  implicit none

  public:: initialize_efit, get_from_efit, initialize_efit_small
  public:: get_c11, finalize_efit
  public:: get_q_from_efit, get_rhotor_from_efit, get_fpol_from_efit
  public:: get_psi_from_efit, get_psirz_from_efit, get_pres_dpdx_from_efit

  private
  integer:: EFITFILE
  real, parameter :: per_phi=2*pi
  real:: rmin, rmax, zmin, zmax
  !for splines:
  real, dimension(:,:,:), allocatable :: splpsi, splrho
  real:: dx, dy
  integer, dimension(:), allocatable  :: imi,ima,jmi,jma
  integer, dimension(:,:), allocatable:: ipoint
  integer :: icp,i,j,ierr


  integer :: nwEQD,nhEQD !Grid dimensions treated as global parameters
  character(len=100) :: geomdir
  character(len=80)  :: geomfile
  real, dimension(:), allocatable     :: rad,zet,rhotor,linpsi
  real, dimension(:), allocatable     :: qpsi,fpol,ffprime,pres
  real, dimension(:,:), allocatable   :: psi,rhotor2D
  logical :: fpol_present = .true.

  real    :: dummy,hrad,hzet,btf,rtf,psib,psia

contains

  !> Reads EFIT file and sets up a spline representation of the poloidal and toroidal flux.
  subroutine initialize_efit(ext_geomdir,ext_geomfile, xval_a, rzstart, Lref, Bref, zmin, zmax, rmag, zmag)
    character(len=100), intent(in) :: ext_geomdir
    character(len=80), intent(in)  :: ext_geomfile
    real, dimension(pi1gl:pi2gl),intent(in):: xval_a
    real, intent(out):: Lref, Bref, zmin, zmax, rmag, zmag
    real, dimension(pi1gl:pi2gl,2), intent(out):: rzstart
    logical:: file_exists

    geomdir=ext_geomdir
    geomfile=ext_geomfile

    if (trim(geomfile).eq.'') then
       WRITE(*,"(4A)") 'Empty geomfile string! Please provide a valid file name'
       stop
    endif

    Inquire(file=trim(geomdir)//'/'//trim(geomfile),exist=file_exists)
    if (.not.file_exists) then
       WRITE(*,"(4A)") trim(geomdir),'/',trim(geomfile),' does not exist'
       stop
    endif

    call get_unit_nr(EFITFILE)
    !Get grid dimensions 
    call dimeq(nwEQD,nhEQD)

    allocate(rad(nwEQD), zet(nhEQD), rhotor(nwEQD), linpsi(nwEQD))
    allocate(qpsi(nwEQD), fpol(nwEQD), ffprime(nwEQD), pres(nwEQD))
    allocate(psi(nwEQD,nhEQD))
!    allocate(rhotor2D(nwEQD,nhEQD))

    !Read equilibrium file
    call eqfile(psib,psia,btf,rtf,rmag,zmag,rad,zet,psi,qpsi,pres)

    call init_spline_coeffs(nwEQD,nhEQD,rad,zet)
    call setup_psi_phi_grids(rmag,zmag,Lref,psia,psib,btf,rhotor)

    !rhotor2D is presently not used
!    do j=1,nhEQD
!       do i=1,nwEQD
!          call lag3interp(rhotor,linpsi,nwEQD,rhotor2D(i,j),psi(i,j))
!       enddo
!    enddo

    zmin=zet(1)
    zmax=zet(nhEQD)

    Bref=abs(btf)

!    call initialize_rho_spline

    rzstart = find_startpoint(xval_a,rmag)

  end subroutine initialize_efit

  !>This routine does a short initialization just to determine reference values
  !!for an automatic calculation of rhostar.
  subroutine initialize_efit_small(ext_geomdir,ext_geomfile, Lref_out, Bref_out, rmag,tracer_initialized)
    character(len=100), intent(in) :: ext_geomdir
    character(len=80), intent(in)  :: ext_geomfile
    real, intent(out):: Lref_out, Bref_out, rmag
    real:: zmag
    real, save:: Bref, Lref
    logical, intent(in):: tracer_initialized

    if (.not.tracer_initialized) then
       geomdir=ext_geomdir
       geomfile=ext_geomfile

       call get_unit_nr(EFITFILE)
       !Get grid dimensions 
       call dimeq(nwEQD,nhEQD)

       allocate(rad(nwEQD), zet(nhEQD), rhotor(nwEQD), linpsi(nwEQD))
       allocate(qpsi(nwEQD), fpol(nwEQD), ffprime(nwEQD), pres(nwEQD))
       allocate(psi(nwEQD,nhEQD))
       !Read equilibrium file
       call eqfile(psib,psia,btf,rtf,rmag,zmag,rad,zet,psi,qpsi,pres)
       Bref=abs(btf)
       call init_spline_coeffs(nwEQD,nhEQD,rad,zet)
       call setup_psi_phi_grids(rmag,zmag,Lref,psia,psib,btf,rhotor)

       deallocate(rad,zet,rhotor,linpsi,qpsi,fpol,ffprime,psi,pres)
       call finalize_spline_coeffs
    endif

    Bref_out=Bref
    Lref_out=Lref
  end subroutine initialize_efit_small

  subroutine setup_psi_phi_grids(rmag,zmag,Lref,psiaxis,psisep,Bref,rhotor)
    real, intent(out):: Lref
    real, dimension(nwEQD), intent(out) :: rhotor
    real, intent(inout):: psiaxis, psisep
    real, intent(in):: Bref,rmag,zmag
    real :: dpsi, dpsi_fine, xmax,xmin
    real:: psif, dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2, fpol_int, phiedge
    real:: phiedge_old, dw, dh, sumhits_real
    real, dimension(:), allocatable :: qpsi_fine, phi_fine, &
         linpsi_fine, rhotor_fine
    integer :: i, nw_fine, ntot, it, c, n_it=40, xind, yind
    integer:: count
    integer:: xmid, ymid, nw, nh, ind, acc_count
    integer,dimension(:,:),allocatable:: hits
    real, dimension(:),allocatable:: x, y, xwidth
!!$    character(len=2):: it_str
    
    !we do not want to over-resolve to avoid artificial 'bridges' between the 
    !plasma and the private flux region
    nw_fine = max(500,2*max(nwEQD,nhEQD)) 

    !Generate toroidal flux vs. poloidal flux array
    dpsi = (psisep-psiaxis)/(nwEQD-1)
    do i=1,nwEQD
       linpsi(i)=psiaxis+(i-1)*dpsi
    enddo
    allocate(phi_fine(nw_fine))
    allocate(linpsi_fine(nw_fine))
    dpsi_fine = (psisep-psiaxis)/(nw_fine-1)
    do i=1,nw_fine
       linpsi_fine(i)=psiaxis+(i-1)*dpsi_fine
    enddo

    if (fpol_present) then
       !compute edge toroidal flux using Monte-Carlo method
       !(should be more robust than integration of the q-profile)
       !each process uses a different seed, which yields faster convergence
       !if more processors are used
       call sgrnd(42*nint((mype+1)**(pi/4)))
       ntot=nw_fine**2
       nw=nw_fine
       nh=nw_fine
       dw=(rad(nwEQD)-rad(1))/(nw_fine-1)
       dh=(zet(nhEQD)-zet(1))/(nw_fine-1)
       acc_count=0
       allocate(hits(nw,nh))
       allocate(x(nw),y(nh),xwidth(nh))
       phiedge=0.
       hits=0
       do xind=1,nw
          x(xind)=rad(1)+(xind-1)*dw
       enddo
       do yind=1,nh
          y(yind)=zet(1)+(yind-1)*dh
       enddo
       !we iterate until an accuracy of 5*10^-4 is reached in three consecutive iterations
       do it=1,n_it
          ntot=ntot*1.5
          phiedge_old=phiedge
          phiedge=0.
          !reset hits to have independent iterations (at the price of slower, but more meaningful
          !convergence)
          hits=0
          !sum is real to avoid integer overflow
          sumhits_real=0.
          !compute psi at ntot random positions, count hits in each grid cell,
          !if psiaxis<psi<psiedge (we don't know the shape of the plasma in advance)
          do c=1,ntot
             xind=nint(grnd()*(nw-1))+1
             yind=nint(grnd()*(nh-1))+1
             call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
                  &ipoint,x(xind),y(yind),&
                  &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
             !need to allow some freedom for psiaxis due to limited interpolation accuracy
             if (psif.ge.(psiaxis-1e-4*abs(psiaxis)).and.psif.le.psisep) then
                hits(xind,yind)=hits(xind,yind)+1
             endif
          enddo
          !collect results from all processes
          call my_sum_to_all_int(hits,nw*nh,MY_MPI_COMM_WORLD)

          !delete all hits outside the plasma, starting from the magnetic axis
          xmid=minloc(abs(x-rmag),1)
          ymid=minloc(abs(y-zmag),1)
          !top half
          do yind=ymid,nh
             xmin=xmid
             xmax=xmid
             call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
                  &ipoint,x(xmid),y(yind),&
                  &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
             if (hits(xmid,yind)==0.and.(psif.lt.psiaxis.or.psif.gt.psisep)) then
                hits(:,yind:)=0
                exit
             endif
             do xind=xmid,nw
                if (hits(xind,yind)>0) then
                   xmax=xind
                else
                   hits(xind:,yind)=0
                   exit
                endif
             enddo
             do xind=xmid-1,1,-1
                if (hits(xind,yind)>0) then
                   xmin=xind
                else
                   hits(1:xind,yind)=0
                   exit
                endif
             enddo
             xmid=(xmax+xmin)/2
             xwidth(yind)=(xmax-xmin)
          enddo
          xmid=minloc(abs(x-rmag),1)
          !lower half
          do yind=ymid-1,1,-1
             xmin=xmid
             xmax=xmid
             call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
                  &ipoint,x(xmid),y(yind),&
                  &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
             if (hits(xmid,yind)==0.and.(psif.lt.psiaxis.or.psif.gt.psisep)) then
                hits(:,1:yind)=0
                exit
             endif
             do xind=xmid,nw
                if (hits(xind,yind)>0) then
                   xmax=xind
                else
                   hits(xind:,yind)=0
                   exit
                endif
             enddo
             do xind=xmid-1,1,-1
                if (hits(xind,yind)>0) then
                   xmin=xind
                else
                   hits(1:xind,yind)=0
                   exit
                endif
             enddo
             xmid=(xmax+xmin)/2
             xwidth(yind)=(xmax-xmin)
          enddo
          !delete everything above and below the narrowest regions of the plasma (=X-points)
          hits(:,1:ymid+1-minloc(xwidth(ymid:1:-1),1))=0
          hits(:,ymid+minloc(xwidth(ymid:)-1,1):)=0
   
          !'count' will contain the number of grid cells with hits
          count=0
          !compute edge toroidal flux using all remaining hits
          do yind=1,nh
             do xind=1,nw
                if (hits(xind,yind).gt.0) then
                   count=count+1
                   call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
                        &ipoint,x(xind),y(yind),&
                        &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
                   call lag3interp(fpol,linpsi,nwEQD,&
                        &fpol_int,psif)
                   phiedge=phiedge+hits(xind,yind)*abs(fpol_int)/x(xind)*dw*dh
                   sumhits_real=sumhits_real+real(hits(xind,yind))
                endif
             enddo
          enddo
!!$ file output for debugging
!!$          if (mype==0) then
!!$             write(it_str,'(I2)') it
!!$             open(450,file='hits'//trim(adjustl(it_str))//'.dat',status='replace')
!!$             do yind=1,nh
!!$                do xind=1,nw
!!$                   write(450,'(2ES16.8,I12)') x(xind),y(yind),hits(xind,yind)
!!$                enddo
!!$                write(450,*)
!!$             enddo
!!$             close(450)
!!$          endif
          phiedge=phiedge/(sumhits_real/count)/2./pi
          if (abs((phiedge-phiedge_old))/phiedge.lt.5.e-4) then
             acc_count=acc_count+1
             !we exit if three consecutive iterations reach an accuracy better than 5e-4 
             if (acc_count.ge.3) exit
          else
             acc_count=0
          endif
          if (it.eq.n_it.and.acc_count.lt.3) &
             stop 'Could not compute toroidal flux to sufficient accuracy, check equilibrium!'
       enddo

       !construct phi-grid from data of final iteration
       phi_fine=0.
       phi_fine(nw_fine)=phiedge
       !successively, peel off hits from outside to inside if they exceed the examined psi value,
       !calculate phi(psi) using the remaining hits
       do ind=nw_fine-1,2,-1
          count=0
          sumhits_real=0.
          do yind=1,nh
             do xind=1,nw
                if (hits(xind,yind)==0) cycle
                call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
                     &ipoint,x(xind),y(yind),&
                     &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
                if (psif.gt.linpsi_fine(ind)) then
                   hits(xind,yind)=0.
                   cycle
                else
                   count=count+1
                endif
                call lag3interp(fpol,linpsi,nwEQD,&
                     &fpol_int,psif)
                phi_fine(ind)=phi_fine(ind)+hits(xind,yind)*abs(fpol_int)/x(xind)*dw*dh
                sumhits_real=sumhits_real+real(hits(xind,yind))
             enddo
          enddo
          phi_fine(ind)=phi_fine(ind)/(sumhits_real/count)/2./pi
       enddo
       !set phi at axis manually
       phi_fine(1)=0.

       !interpolate to EQSDK radial resolution
       allocate(rhotor_fine(nw_fine))
       rhotor_fine=sqrt(phi_fine/phiedge)
       call lag3interp(rhotor_fine,linpsi_fine,nw_fine,&
            &rhotor,linpsi,nwEQD)    
       deallocate(rhotor_fine)
       deallocate(hits,x,y,xwidth)
       Lref=sqrt(2.*abs(phiedge/Bref))
    endif
    deallocate(linpsi_fine,phi_fine)


    !for comparison (and cases without fpol): q-profile integration 
    nw_fine=max(2000,nwEQD)
    allocate(phi_fine(nw_fine),qpsi_fine(nw_fine),linpsi_fine(nw_fine))
    dpsi_fine = (psisep-psiaxis)/(nw_fine-1)
    do i=1,nw_fine
       linpsi_fine(i)=psiaxis+(i-1)*dpsi_fine
    enddo
    call lag3interp(qpsi,linpsi,nwEQD,&
         &qpsi_fine,linpsi_fine,nw_fine)
    phi_fine(1)=0.
    do i=2,nw_fine
       phi_fine(i)=phi_fine(i-1)+0.5*(qpsi_fine(i-1)+qpsi_fine(i))*dpsi_fine
    enddo
    deallocate(qpsi_fine)


    if (.not.fpol_present) then
           allocate(rhotor_fine(nw_fine))
       rhotor_fine = sqrt(phi_fine/phi_fine(nw_fine))
       Lref=sqrt(2.*abs(phi_fine(nw_fine)/Bref))
       call lag3interp(rhotor_fine,linpsi_fine,nw_fine,&
            &rhotor,linpsi,nwEQD)    
       deallocate(rhotor_fine)
    else
       ! If a difference between the two evaluations is found, either:
       ! 1) the q-profile at the edge prevents a proper integration (which does not affect the MC result)
       ! 2) the MC integration fails to recognize the separatrix -- this can be checked by uncommenting 
       !    the lines above marked with '!!$' comments and plotting the output
       ! So far, the MC scheme has been tested with EFIT files containing data from 
       ! AUG, NSTX, C-Mod, DIII-D, JET.
       if (abs(phiedge-phi_fine(nw_fine))/phiedge.gt.0.01.and.mype==0) then
          write(*,'(A)') "WARNING: Results from toroidal flux evaluations differ by more than 1%"
          write(*,'(A,ES16.8,A,ES16.8)') 'MC: ',phiedge,'; q-profile: ',phi_fine(nw_fine)
       endif
    endif

    deallocate(phi_fine,linpsi_fine)

  end subroutine setup_psi_phi_grids


  subroutine init_spline_coeffs(nx, ny, rad, zet)
    integer, intent(in):: nx, ny
    real, dimension(nx):: rad
    real, dimension(ny):: zet
    integer:: i
    dx = rad(2) - rad(1)
    dy = zet(2) - zet(1)

    !------------------------------ Initialize 2D spline -----------------------------------!

    !Rectangular domain

    allocate(imi(ny),ima(ny),jmi(nx),jma(nx))    

    imi = 1
    ima = nx
    jmi = 1
    jma = ny

    !Number of data in splpsi
    icp=0
    do i=1,ny
       if ( imi(i) > 0 .and. ima(i) > 0 ) then
          icp = icp + ima(i) - imi(i) + 1
       endif
    enddo

    allocate(splpsi(6,6,icp),ipoint(nx,ny))

    call s2dcut(nx,ny,dx,dy,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)

  end subroutine init_spline_coeffs

  subroutine finalize_spline_coeffs
    deallocate(imi,ima,jmi,jma,splpsi,ipoint)
  end subroutine finalize_spline_coeffs


  subroutine finalize_efit
    deallocate(rad,zet, rhotor, linpsi, psi)
!    deallocate(rhotor2D)
    deallocate(qpsi,fpol,ffprime,pres)
    call finalize_spline_coeffs
!    call finalize_rho_spline
  end subroutine finalize_efit

  !> For a given x0 (taken to be rho_tor), we find the starting point 
  !! for the tracing on the R, Z grid.
  function find_startpoint(xval_a,rmag)
    real,dimension(pi1gl:pi2gl),intent(in):: xval_a
    real,dimension(pi1gl:pi2gl,2):: find_startpoint
    real,intent(in):: rmag
    real :: rstart,zstart,Deltar,Deltaz,dr,dz,maxr
    !    real:: rhof,drhodr,drhodz,d2rhodr2,d2rhodrdz,d2rhodz2
    real:: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
    real:: relpsif, psi_acc,psit, relpsit !t=target
    real:: rad_fine, zet_fine, zstart_nextit
    integer:: r_ind, z_ind, iteration
    integer:: nr=300, nz=300, i1, i2

    nr = max(nr,size(rad))
    nz = max(nz,size(zet))

    if (lilo) then
       i1=pi1gl; i2=pi1gl; rstart=0.; zstart=0.
    else
       i1=pi1gl; i2=pi2gl
    endif

    do i=i1,i2

       !we only look outside of the magnetic axis
       rstart=(rmax+rmag)*0.5
       zstart=(zmax+zmin)*0.5
       zstart_nextit=zstart

       Deltar = 2.0*(rmax-rstart)
       !we examine only half the z-domain, as some EQDSKs contain also 
       !coils, which lead to false starting points
       !this assumes the plasma center is not in the top/bottom quarter of the domain
       Deltaz = 1.0*min((zstart-zmin),(zmax-zstart))
       call get_psi_from_efit(xval_a(i), psit, relpsit)
       do iteration=0,16
          dr=Deltar/(real(nr)*2.**iteration)
          dz=Deltaz/(real(nz)*2.**iteration)

          !find the outermost point of each flux surface xval_a(i)
          maxr = 0.

          do z_ind=0,nz
             zet_fine = zstart+(z_ind-nz/2)*dz
             if ((zet_fine.gt.zmax).or.(zet_fine.lt.zmin)) cycle
             do r_ind=0,nr
                rad_fine = rstart+(r_ind-nr/2)*dr
                if (rad_fine.gt.rmax) cycle
                !we use psi instead of rhotor2D as the quality of 
                !the latter is too bad (it's based on an interpolation
                !from psi and thus we can/should use psi, anyway)
                call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
                     &ipoint,rad_fine,zet_fine,&
                     &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
                !we filter out radial positions which have dpsidr<0, 
                !since this may sometimes confuse this algorithm, 
                !and is not desired anyway, since we're looking for 
                !starting positions on the low-field side
                if (abs(psif-psit).le.1.5*dpsidr*dr.and.dpsidr.ge.0.) then
                   if (rad_fine.gt.maxr) then
                      maxr=rad_fine
                      zstart_nextit=zet_fine
                   endif
                endif
             enddo
          enddo
          rstart=maxr
          zstart=zstart_nextit

          !check 
          call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
               &ipoint,rstart,zstart,&
               &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
          relpsif = (psif-psia)/(psib-psia)
          psi_acc = abs((relpsif-relpsit)/relpsit)
          if (psi_acc.lt.1E-5) exit
       enddo
       if (psi_acc.gt.1E-5) then
          if (mype.eq.0) then
             print*, 'Warning: rzstart for tracer not exactly at chosen rhotor'
             print*, 'rel. psi target/real: ', relpsit, relpsif
          endif
       endif
       find_startpoint(i,1)=rstart
       find_startpoint(i,2)=zstart
    enddo

  end function find_startpoint

  !> Compute the B-field components at an arbitrary R, Z position. 
  subroutine get_from_efit(rrr, zzz, Brad, Bphi, Bzet, dBraddr, dBraddp, dBraddz,&             
       &dBphidr, dBphidp, dBphidz, dBzetdr, dBzetdp, dBzetdz)


    real, intent(IN) :: rrr,zzz
    real, intent(OUT):: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
    real, intent(OUT):: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz
    real    :: d2psidz2,d2psidrdz,dpsidz,dpsidr,psif,d2psidr2
    real    :: fpolf,ffprimef,inv_rrr

    !---------------------------------------------------------------------------!

    call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,ipoint,rrr,zzz,&
         &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

    inv_rrr = 1.0/rrr
    !B components
    Brad = -inv_rrr*dpsidz
    Bzet = inv_rrr*dpsidr
    dBraddp = 0.
    dBphidp = 0.
    dBzetdp = 0.
    dBraddr = -inv_rrr*(d2psidrdz+Brad)
    dBzetdz = inv_rrr*d2psidrdz
    dBraddz = -inv_rrr*d2psidz2
    dBzetdr = inv_rrr*(d2psidr2-Bzet)

    if (fpol_present) then !interpolate & use fpol, ffprime
       call lag3interp(fpol,linpsi,nwEQD,fpolf,psif)
       call lag3interp(ffprime,linpsi,nwEQD,ffprimef,psif)
       Bphi    = inv_rrr*fpolf
       dBphidr = -inv_rrr*Bphi + ffprimef/fpolf * inv_rrr* dpsidr
       dBphidz = ffprimef*inv_rrr*dpsidz/fpolf
    else
       Bphi    = -btf*rtf/rrr
       dBphidr = -Bphi/rrr
       dBphidz = 0.0
    endif

  end subroutine get_from_efit

  !>Get safety factor from EFIT
  subroutine get_q_from_efit(rhotor_in,q_out,shat_out)
    real, intent(in) :: rhotor_in
    real, intent(out) :: q_out, shat_out
    real:: dqdx

    call lag3interp(qpsi,rhotor,nwEQD,q_out,rhotor_in)
    call lag3deriv(qpsi,rhotor,nwEQD,dqdx,rhotor_in)
    shat_out=rhotor_in/q_out*dqdx

  end subroutine get_q_from_efit

  !>Get pressure gradient from EFIT
  subroutine get_pres_dpdx_from_efit(rhotor_in,pres_out,dpdx_out)
    real, dimension(pi1gl:pi2gl), intent(in) :: rhotor_in
    real, dimension(pi1gl:pi2gl), intent(out) :: pres_out, dpdx_out

    call lag3interp(pres,rhotor,nwEQD,pres_out,rhotor_in,px0)
    call lag3deriv(pres,rhotor,nwEQD,dpdx_out,rhotor_in,px0)

  end subroutine get_pres_dpdx_from_efit

  subroutine get_fpol_from_efit(rhotor_in,fpol_out,ffprime_out)
    real, intent(in) :: rhotor_in
    real, intent(out) :: fpol_out, ffprime_out

    call lag3interp(fpol,rhotor,nwEQD,fpol_out,rhotor_in)
    call lag3interp(ffprime,rhotor,nwEQD,ffprime_out,rhotor_in)

  end subroutine get_fpol_from_efit

  subroutine get_psi_from_efit(rhotor_in,psi_out,relpsi_out)
    real, intent(in) :: rhotor_in
    real, intent(out) :: psi_out, relpsi_out

    call lag3interp(linpsi,rhotor,nwEQD,psi_out,rhotor_in)
    relpsi_out = (psi_out-psia)/(psib-psia)

  end subroutine get_psi_from_efit

  subroutine get_psirz_from_efit(r0,z0,psi_out,relpsi_out)
    real, intent(in) :: r0,z0
    real, intent(out) :: psi_out, relpsi_out
    real:: dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
    integer :: ierr

    call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,ipoint,r0,z0,&
         &psi_out,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

    relpsi_out = (psi_out-psia)/(psib-psia)

  end subroutine get_psirz_from_efit

  !>Get rhotor from EFIT
  subroutine get_rhotor_from_efit(r0,z0,rhof)
    real, intent(in) :: r0,z0
    real, intent(out) :: rhof
    real:: psif, dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
    integer :: ierr

    call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,&
         &ipoint,r0,z0,&
         &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
    call lag3interp(rhotor,linpsi,nwEQD,rhof,psif)

  end subroutine get_rhotor_from_efit

  !> Get the resolution in the poloidal plane from the EFIT file.
  subroutine dimeq(nw,nh)

    integer, intent(OUT) :: nw,nh

    integer :: idum,i
    character(10), dimension(6) ::  cdum

    open(EFITFILE, file=trim(geomdir)//'/'//trim(geomfile), status="old", action="read")

    read(EFITFILE,"(6A8,3I4)")(cdum(i),i=1,6),idum,nw,nh
    close(EFITFILE)


    return

  end subroutine dimeq

  !> Read EFIT input and generate a rho_tor grid from the poloidal flux and the safety
  !! factor profile.
  subroutine eqfile(psisep,psiaxis,bt0,rzero,rmaxis,zmaxis,&
       &rad,zet,psiRZ,qpsi,pres)
    real, intent(OUT) :: psisep, psiaxis, bt0, rzero
    real, intent(OUT) :: rmaxis, zmaxis
    real, dimension(nwEQD), intent(OUT) :: rad
    real, dimension(nhEQD), intent(OUT) :: zet
    real, dimension(nwEQD,nhEQD), intent(OUT) :: psiRZ
    real, dimension(nwEQD), intent(out) :: qpsi, pres
    real :: psiAxis2, psiSep2
    real :: xdim,zdim,rcentr,rleft,zmid,xdum
    real :: bcentr,plas_cur
    real, dimension(nwEQD) :: pprime !,phi
    !    real, dimension(nwEQD) :: fpol,ffprim,qpsi
    !    real, dimension(:), allocatable :: LCFS, limEQD
    character(10), dimension(6) :: cdum
    !    integer :: n_bndyxy,nlimEQD
    integer :: i,j, idum
    character(10) :: realfmt="(5E16.9)"


    call get_unit_nr(EFITFILE)
    open(EFITFILE, file=trim(geomdir)//'/'//trim(geomfile), status="old", action="read")

    read(EFITFILE,"(6A8,3I4)")(cdum(i),i=1,6),idum,idum,idum
    read(EFITFILE,realfmt) xdim,zdim,rcentr,rleft,zmid
    read(EFITFILE,realfmt) rmaxis,zmaxis,psiAxis,psiSep,bcentr
    read(EFITFILE,realfmt) plas_cur,psiAxis2,xdum,rmaxis,xdum
    read(EFITFILE,realfmt) zmaxis,xdum,psiSep2,xdum,xdum
    read(EFITFILE,realfmt) (fpol(i),i=1,nwEQD)
    read(EFITFILE,realfmt) (pres(i),i=1,nwEQD)
    read(EFITFILE,realfmt) (ffprime(i),i=1,nwEQD)
    read(EFITFILE,realfmt) (pprime(i),i=1,nwEQD)
    read(EFITFILE,realfmt) ((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
    read(EFITFILE,realfmt) (qpsi(i),i=1,nwEQD)

    !    !Boundary Data, currently not used
    !    read(EFITFILE,*) n_bndyxy,nlimEQD    

    !    allocate(LCFS(2*n_bndyxy),limEQD(2*nlimEQD))

    !    read(EFITFILE,"(5E16.9)") (LCFS(i),i=1,2*n_bndyxy)
    !    read(EFITFILE,"(5E16.9)") (limEQD(i),i=1,2*nlimEQD)

    close(EFITFILE)

    call set_eqcoords(xdim,zdim,rleft,zmid,rad,zet)           

    !switching to 2pi division for ITER equilibria 
    psi_o_twopi = psi_o_twopi.or.(index(cdum(1),'FEATF4E').gt.0)

    if ((mype.eq.0)) then
       if (psiAxis.ne.psiAxis2) write (*,"(A)") &
            &'Warning: Two different psi-axis values found in efit file header - taking 1st'
       if (psiSep.ne.psiSep2) write (*,"(A)") &
            &'Warning: Two different psi-separatrix values found in efit file header - taking 1st'
       if (q_scalefac.ne.1.0) write(*,"(A,F6.3)") &
            &'Warning: Poloidal flux is being rescaled by a factor of ', q_scalefac
    endif

    !optionally rescale the poloidal flux for a better match of the qpsi profile
    psi=psi*q_scalefac
    psiaxis=psiaxis*q_scalefac
    psisep=psisep*q_scalefac
    ffprime=ffprime/q_scalefac
    !pprime=pprime/q_scalefac

    !optionally divide psi by 2pi to account for different input definitions 
    if (psi_o_twopi) then
       if (mype==0) write (*,"(A)") &
            'Poloidal flux is being divided by 2pi'
       psirz=psirz/(2.0*pi)
       psiaxis=psiaxis/(2.0*pi)
       psisep=psisep/(2.0*pi)
       !pprime=pprime*2*pi
       ffprime=ffprime*(2.0*pi)
    endif

    fpol_present = (abs(sum(fpol)).gt.0).and.(abs(sum(ffprime)).gt.0)

    if (psiaxis.gt.psisep) then
       if (mype.eq.0) write(*,*) &
            &'Warning: flipped sign of psi in tracer'
       psirz = - psirz
       psiaxis = -psiaxis
       psisep = -psisep
       ffprime = -ffprime
    endif

    if (fpol_present) then
       if (any(fpol.lt.0.)) then
          if (mype.eq.0) write(*,*) &
               &'Warning: flipped sign of fpol in tracer'
          fpol = -fpol
       endif
    endif

    if (fpol_present) then
       !take reference magnetic field at magnetic axis
       bt0   = fpol(1)/rmaxis
       rzero = rmaxis
    else
       !take reference magnetic field at center
       bt0   = bcentr
       rzero = rcentr
       if ((rcentr.ne.rmaxis).and.(mype==0)) write(*,'(A)') &
            &"Warning: normalizing to efit BCENTR instead of magnetic axis value"
    endif

    if (bt0.gt.0) bt0=-bt0


  end subroutine eqfile

  !> Set up R, Z grids for the spline representation. 
  subroutine set_eqcoords(xdim,zdim,rleft,zmid,rad,zet)

    real, intent(IN) :: xdim,zdim,rleft,zmid
    real, dimension(nwEQD), intent(OUT) :: rad
    real, dimension(nhEQD), intent(OUT) :: zet

    real:: z1
    integer :: j,k
    intrinsic real

    do j=1,nwEQD
       rad(j) = rleft + (j-1)*(xdim/real(nwEQD-1))
    enddo

    z1 = zmid - zdim/2.

    do k=1,nhEQD
       zet(k) = z1 + (k-1)*(zdim/real(nhEQD-1))
    enddo

    rmin=minval(rad) ; rmax=maxval(rad)
    zmin=minval(zet) ; zmax=maxval(zet)

  end subroutine set_eqcoords

  !> Generate spline representation of rho_tor.
  subroutine initialize_rho_spline
    !------------------------------ Initialize 2D spline -----------------------------------!
    !we operate on the same grid as for the psi spline, therefore we don't initialize separate
    !arrays for imi, ima, ipoint, etc.

    allocate(splrho(6,6,icp))

    call s2dcut(nwEQD,nhEQD,dx,dy,rhotor2D,imi,ima,jmi,jma,icp,splrho,ipoint)

  end subroutine initialize_rho_spline

  !> Get the initial condition for c11=drho_tor/dR.
  subroutine get_c11(r,z,c11_int,Lref,Bref)
    real:: psi,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2, rho, q, shat
    integer:: ierr
    real, intent(in):: r,z
    real, intent(out):: c11_int
    real, intent(in):: Lref, Bref


    call spline(nwEQD,nhEQD,rad,zet,dx,dy,icp,splpsi,ipoint,r,z,&
         &psi,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

    call lag3interp(rhotor,linpsi,nwEQD,&
         &rho,psi)    
    call get_q_from_efit(rho,q,shat)

    c11_int=dpsidr*q/rho/Lref/Bref

  end subroutine get_c11

  subroutine finalize_rho_spline
    deallocate(splrho)
  end subroutine finalize_rho_spline

end module tracer_IO
