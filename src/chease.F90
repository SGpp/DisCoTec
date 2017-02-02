#include "redef.h"
!>Interface to the magnetic equilibrium code CHEASE
MODULE chease_mod
  use discretization
  use coordinates
  use par_geom, only: geomtype
  use lagrange_interpolation
  implicit none
  public:: init_chease 
  public :: force_x
  public :: force_z

#ifdef WITHFUTILS
  public:: get_chease, get_chease_refvals
#endif

  private

  LOGICAL :: force_x !used to force to the nearest flux surface, only local
  LOGICAL :: force_z !used to force using the z grid of CHEASE, only local

  REAL:: dqdpsi0,r0,dxdpsi0,p0     ! some ref quantities in the middle of the box 
  CHARACTER(len=80):: x_def           ! x definitions :
  ! - arho_p=a*rho_t=a*sqrt(Psi/Psi_e), Psi poloidal flux 
  ! - arho_t=a*rho_p=a*sqrt(Phi/Phi_e), Phi toroidal flux
  ! - C_psi = q0/(r0 B0)*Psi (old definition) 
  REAL:: xoa0                        ! center of the simulation domain in global simulations
  ! of the reference flux surface x/a (between 0 and 1)
  REAL:: flux_pos                  ! only for the local case : position
  ! of the reference flux surface x/a (between 0 and 1)
  INTEGER:: NPSI,NCHI, NPSI0,NPSI1, NCHI0 ! indices for the CHEASE mesh
  integer:: nx0int,nchiint
  REAL:: R0EXP, B0EXP
  REAL:: p_spec
  LOGICAL:: pressure_term
  LOGICAL:: print_ini_msg
  REAL :: rhostar                    ! rho_ref/a
  CHARACTER(len=80):: geomfile
  CHARACTER(len=100):: geomdir
  REAL, DIMENSION(:), ALLOCATABLE   :: dxdpsi
  REAL, DIMENSION(:), ALLOCATABLE   :: Rgeom
  REAL, DIMENSION(:), ALLOCATABLE   :: ageom
  REAL, DIMENSION(:), ALLOCATABLE   :: V
  REAL, DIMENSION(:),   ALLOCATABLE   :: dVdpsi
  REAL, DIMENSION(:),   ALLOCATABLE   :: q                            ! safety factor
  REAL, DIMENSION(:),   ALLOCATABLE   :: dqdpsi                      ! dqdpsisi
  REAL, DIMENSION(:),   ALLOCATABLE   :: p                           ! pressure
  REAL, DIMENSION(:),   ALLOCATABLE   :: dpdpsi                      ! dpressure/dpsi
  REAL, DIMENSION(:),   ALLOCATABLE   :: rho_t                       ! sqrt(PHI_t/B0)
  REAL, DIMENSION(:),   ALLOCATABLE   :: dpsidrhotor
  REAL, DIMENSION(:),   ALLOCATABLE   :: tmp1d
  REAL, DIMENSION(:),   ALLOCATABLE   :: PSI,CHI                        
  REAL, DIMENSION(:,:), ALLOCATABLE   :: g11,g12,g22,g33           ! metric (PSI,CHI,PHI)
  REAL, DIMENSION(:,:), ALLOCATABLE   :: B_2d,dBdPSI_2d,dBdCHI_2d,J_2d  
  REAL, DIMENSION(:,:), ALLOCATABLE   :: R, Z, dPsidr, dPsidz
  real, dimension(:), allocatable:: dxdpsi_int, Rgeom_int, ageom_int, V_int, dVdpsi_int, q_int
  real, dimension(:), allocatable:: dqdpsi_int, psi_int, chi_int, p_int, tmp1d_int, dpdpsi_int
  real, dimension(:), allocatable:: rho_t_int, dpsidrhotor_int
  real, dimension(:,:), allocatable:: g11_int, g12_int, g22_int, g33_int, B_2d_int
  real, dimension(:,:), allocatable:: dBdpsi_2d_int, dBdchi_2d_int, J_2d_int, R_int, Z_int
  real, dimension(:,:), allocatable:: dPsidr_int, dPsidz_int
 
 
  
  ! SUBROUTINES DEFINITIONS
CONTAINS
  
  !> Initialize chease module
  SUBROUTINE init_chease(px_def,pxoa0,pflux_pos,pp_spec,ppressure_term,&
       prhostar,pprint_ini_msg,pgeomfile,pgeomdir,pforce_x,pforce_z)
    REAL    :: pflux_pos, pp_spec, pxoa0,prhostar
    CHARACTER(len=80) ::px_def,pgeomfile,pgeomdir
    LOGICAL :: ppressure_term,pprint_ini_msg,pforce_x,pforce_z
#ifdef WITHFUTILS
    flux_pos=pflux_pos
    x_def=px_def; xoa0=pxoa0
    p_spec=pp_spec; 
    pressure_term=ppressure_term
    print_ini_msg=pprint_ini_msg
    rhostar=prhostar
    geomfile=pgeomfile;geomdir=pgeomdir
    if (pforce_x) then
       force_x=pforce_x
    else 
       force_x=.False.
    endif 
    if (pforce_z) then
       force_z=pforce_z
    else
       force_z=.False.
    end if 
    
#else
    print*,"You have to compile with FUTILS switch set to yes to use CHEASE geometry."
    stop
#endif
  END SUBROUTINE init_chease
  
#ifdef WITHFUTILS
  subroutine get_chease_refvals(Br,Lr,major_R,minor_r,geomfile_in,geomdir_in)
    real, intent(inout):: Br, Lr, major_R, minor_r
    character(len=80),intent(in):: geomfile_in
    character(len=100),intent(in):: geomdir_in

    geomfile=geomfile_in
    geomdir=geomdir_in
    call read_CHEASE_input_HDF5_small
    ! set reference field and equilibrium normalization length
    ! here, major_R is an input that determines the normalization used
    Br       = B0EXP         !  B at magnetic axis
    Lr       = R0EXP/major_R ! 
    minor_r=maxval(rho_t)/R0EXP*major_R
    deallocate(rho_t)

  end subroutine get_chease_refvals

  subroutine read_CHEASE_input_HDF5_small
    
    USE futils
    
    INTEGER :: hdf5_ioutgyro
    
    IF (geomfile.eq.'') geomfile='ogyropsi.h5'
    
    CALL openf(trim(geomdir)//'/'//trim(geomfile), hdf5_ioutgyro)
    CALL getatt(hdf5_ioutgyro, "/data", "NPSI", NPSI)
    CALL getatt(hdf5_ioutgyro, "/data", "R0EXP", R0EXP)
    CALL getatt(hdf5_ioutgyro, "/data", "B0EXP", B0EXP)
    
    allocate(rho_t(1:NPSI))
    CALL getarr(hdf5_ioutgyro, "/data/var1d/rho_tor", rho_t)
    CALL closef(hdf5_ioutgyro)

  end subroutine read_CHEASE_input_HDF5_small

  SUBROUTINE get_chease(geom,q0,&
       &trpeps,shat,beta,dpdx_pm,lx,x0,C_y,C_xy,q_prof,dqdx_prof,Lr,Br,major_R,minor_r,&
       &nx0)
    TYPE(geomtype),INTENT(INOUT) :: geom
    integer, intent(in) ::nx0
    REAL, DIMENSION(pi1gl:pi2gl), INTENT(OUT) :: dpdx_pm,C_y,C_xy,q_prof,dqdx_prof
    REAL, INTENT(OUT) :: q0,trpeps,shat,beta
    real, intent(INOUT) :: lx,x0,Lr,Br,major_R,minor_r
    INTEGER :: k, kind, kpi, i,ix0,nx0int,nchiint 
    REAL    :: CPI, mu0
    REAL    :: pref
    REAL :: Cy
    LOGICAL:: interpol
    REAL, DIMENSION(pi1gl:pi2gl,pj1:pj2,lk1:lk2):: gxx, gxz, gyy, gyz, dBdx

    
    CPI = 4.*atan(1.)
    !*******************************
    ! Read file from CHEASE
    call read_CHEASE_input_HDF5
    
    if (.not.force_x .and. .not.force_z) then
       interpol=.true.
    else if  (.not.force_x .and. force_z) then
       if (.not.xy_local) then
          print*, 'For global simulations must interpolate in both directions, cannot use the force flags'
          stop
       else if (xy_local) then
          IF (NCHI.NE.nz0) THEN 
             print*, 'NCHI in geomfile should have the same value as nz0 in the parameter file'
             stop
          END IF
          interpol=.true.
       endif
    else if (force_x .and. .not.force_z) then
       if (.not.xy_local) then
          print*, 'For global simulations must interpolate in both directions, cannot use the force flags'
          stop
       else
          interpol=.true.
       endif
    else if (force_x .and. force_z) then
       interpol=.false.
       IF (NCHI.NE.nz0) THEN 
          print*, 'NCHI in geomfile should have the same value as nz0 in the parameter file'
          stop
       END IF
    endif

    !*******************************
    ! set reference field and equilibrium normalization length
    Br       = B0EXP         !  B at magnetic axis
    Lr       = R0EXP/major_R ! 

    !*******************************
    !define the dimensions and interpolate
    if (.not.force_x .and. .not.force_z) then
       if (.not.xy_local) then
          nx0int=nx0
          nchiint=nchi+1
       else
          nx0int=1
          nchiint=nchi+1
       end if
    else if (.not.force_x .and.  force_z)then
       nx0int=1 
       nchiint=nchi+1
    else if (force_x .and. .not.force_z) then
       nx0int=NPSI
       nchiint=nchi+1
    else if (force_x .and. force_z) then
       nx0int=NPSI
       nchiint=nchi
    endif

    call interpolate_geometry(NPSI,nx0,NCHI,nz0,lx,x0,nx0int,nchiint,interpol) 

    !*******************************
    ! set x definition, find center of simulation domain, 
    ! set lx and x0 in the global case

    call set_xbox(x0,lx,q0,interpol)

    !*******************************
    ! chi angle check
    nchi0=nz0/2
    IF (MODULO(nz0,2).eq.1) THEN
       if (mype.eq.0) print*, '!!!!!  nz0 must be even   !!!!!!!!!!'
       STOP
    END IF
   
    ! change CHI so that values are going from -pi to pi
    chi_int(nz0/2+1:nz0)=chi_int(nz0/2+1:nz0)-2*CPI

    Cy=r0/q0
    ! setting shear and trpeps (for local code, or for info in global code)
    shat     = Cy/dxdpsi0*dqdpsi0 
    trpeps   = r0/R0EXP

    !set minor_r to appropriate value for global code
    minor_r=maxval(rho_t)/R0EXP*major_R
    
    do k=0, nz0-1
       !shift from CHEASE to GENE
       kpi=modulo(nz0/2+k,nz0)+1
       zval(k)=chi_int(kpi)
    end do
    dz = (2.0 * CPI) / nz0
   !******************************
    Do i=pi1gl,pi2gl
       
       if (xy_local) then
           ix0=NPSI0
       else
           ix0=NPSI1+i-pi1gl   ! for x-parallelisation this pi1gl should be kept
       end if
   
       !*******************************
       ! Calculation of values needed in GENE, 
       ! note all output quantities are normalized
       DO k=lk1,lk2
          kind = k-lk1+1 + lk0*my_pez
          !indices are moved so that g(1) is calculated at chi=-pi
          kpi=MODULO(nz0/2+kind-1,nz0)+1
          ! Calculation of the metric (normalized to Lnorm)
          gxx(i,:,k) = dxdpsi_int(ix0)**2 * g11_int(ix0,kpi)
          geom%gij(i,:,k) = dxdpsi_int(ix0)*Cy * &
               &(dqdpsi_int(ix0) * chi_int(kpi)*g11_int(ix0,kpi)+ &
               &q_int(ix0)*g12_int(ix0,kpi))
          gyy(i,:,k) = Cy**2 * (dqdpsi_int(ix0)**2 * chi_int(kpi)**2 * g11_int(ix0,kpi)   + &
               &2*q_int(ix0)*dqdpsi_int(ix0)*chi_int(kpi)*g12_int(ix0,kpi)+q_int(ix0)**2*g22_int(ix0,kpi) + &
               &g33_int(ix0,kpi) )
          gxz(i,:,k) = Lr * dxdpsi_int(ix0) * g12_int(ix0,kpi)
          gyz(i,:,k) = Lr * Cy * (dqdpsi_int(ix0) * chi_int(kpi) * g12_int(ix0,kpi) + &
               &q_int(ix0) * g22_int(ix0,kpi) )
          geom%gzz(i,:,k) = Lr**2 * g22_int(ix0,kpi)
          
          geom%Bfield(i,:,k)   = B_2d_int(ix0,kpi)/Br
          dBdx(i,:,k) = Lr/(Br*dxdpsi_int(ix0))*dBdPSI_2d_int(ix0,kpi)
          geom%dBdz(i,:,k)  = dBdCHI_2d_int(ix0,kpi)/Br
          geom%jacobian(i,:,k) = 1.0/(dxdpsi_int(ix0)*Cy)*J_2d_int(ix0,kpi)/Lr    
          geom%R(i,k) = R_int(ix0,kpi)
          geom%Z(i,k) = Z_int(ix0,kpi)
          geom%R_hat(i,k) = R_int(ix0,kpi)/Lr
          geom%Z_hat(i,k) = Z_int(ix0,kpi)/Lr
          geom%dxdR(i,k) = dxdpsi_int(ix0)*dPsidr_int(ix0,kpi) !dxdr
          geom%dxdZ(i,k) = dxdpsi_int(ix0)*dPsidz_int(ix0,kpi) !dxdz
       END DO ! k loop
       if (yx_order) then
          geom%gii=gyy
          geom%giz=gyz
          geom%gjj=gxx
          geom%gjz=gxz
          geom%dBdi=0.
          geom%dBdj=dBdx
       else
          geom%gii=gxx    
          geom%giz=gxz
          geom%gjj=gyy
          geom%gjz=gyz
          geom%dBdi=dBdx
          geom%dBdj=0.
       end if

       q_prof(i)=q_int(ix0)
       dqdx_prof(i)=Lr*dqdpsi_int(ix0)/dxdpsi_int(ix0)
       C_y(i)= Cy/Lr                             ! normalized C_y
       C_xy(i)=1/(2*cpi*dxdpsi_int(ix0)*Cy*Br)             ! Normalized C_xy
       
       IF (pressure_term) THEN
           mu0=4*CPI*1e-7
           pref=p0/p_spec
           beta=2*mu0*pref/Br**2
           dpdx_pm(i)=-Lr*(1.0/dxdpsi_int(ix0))*dpdpsi_int(ix0)*2*mu0/Br**2
       ELSE
           dpdx_pm(i)=0
       END IF
        
    end Do  ! i loop
    !*******************************
    !deallocation
    call CLOSE_GEOM(interpol)
    
    IF ((mype.eq.0).and.(print_ini_msg)) THEN 
       ! Print infos
       if (interpol) then
          if (force_x .and. .not.force_z) then
             WRITE(*,'(A)') "  "
             write(*,'(A)')  ' Interpolating in z and forcing to nearest flux surface in x'
             WRITE(*,'(A)') "  "
          else if (.not.force_x .and. force_z) then
             WRITE(*,'(A)') "  "
             write(*,'(A)')  ' Interpolating in x and same z grid as CHEASE'
             WRITE(*,'(A)') "  "
          end if
          if (xy_local) then
             WRITE(*,'(A,i4)') ' CHEASE PSI grid = ', NPSI
             WRITE(*,'(A,i4)') ' GENE x grid = ', nx0       
             WRITE(*,'(A,i4)') ' CHEASE CHI grid = ', nchi
             WRITE(*,'(A,i4)') ' GENE z grid = ', nz0
             WRITE(*,'(A)') "  "
          end if
       else
          write(*,'(A)')  ' Forcing to nearest flux surface in x and same z grid as CHEASE'
          WRITE(*,'(A)') "  "
       endif
       print*, 'Some  reference values : '
       WRITE(*,'(A,F8.4)')  ' Bref = ', Br
       write(*,'(A,F8.4)')  ' Lref = ', Lr
       write(*,'(A,F8.4)')  ' q0   = ', q0
       write(*,'(A,F10.6)') ' shat = ', shat
       write(*,'(A,F10.6)') ' trpeps = ', trpeps
       write(*,'(A,F10.6)') ' minor_r = ', minor_r
       IF (.not.pressure_term) WRITE(*,'(A)') " Beta from parameter file (not consistent with CHEASE) : "
       write(*,'(A,F8.4)')  ' beta = ', beta
       write(*,'(A,F8.4)')  ' dpdx_pm0 = ', dpdx_pm((pi1gl+pi2gl)/2)
       WRITE(*,'(A)') "  "
       WRITE(*,'(A)') "***********************  "
       WRITE(*,'(A)') "  "
    END IF
       
  END SUBROUTINE get_chease

  ! READ CHEASE 

  SUBROUTINE read_CHEASE_input_HDF5
    
    USE futils
    
    INTEGER :: hdf5_ioutgyro
    
    IF (geomfile.eq.'') geomfile='ogyropsi.h5'
    
    CALL openf(trim(geomdir)//'/'//trim(geomfile), hdf5_ioutgyro)

    CALL getatt(hdf5_ioutgyro, "/data", "NPSI", NPSI)
    CALL getatt(hdf5_ioutgyro, "/data", "NCHI",  NCHI)
    CALL getatt(hdf5_ioutgyro, "/data", "R0EXP", R0EXP)
    CALL getatt(hdf5_ioutgyro, "/data", "B0EXP", B0EXP)
    
    IF ((mype.eq.0).and.(print_ini_msg)) THEN
       WRITE(*,'(A)') '  '
       WRITE(*,'(3A)') '*********** reading CHEASE file ', trim(geomdir)//'/'//trim(geomfile) , ' *********'
       WRITE(*,'(A)') "  "
    END IF    
    
    ALLOCATE(dxdpsi(1:NPSI),Rgeom(1:NPSI), ageom(1:NPSI), &
         V(1:NPSI),dVdpsi(1:NPSI), q(1:NPSI),     dqdpsi(1:NPSI),  &
         PSI(1:NPSI),   CHI(1:NCHI),  &
         p(1:NPSI),tmp1d(1:NPSI), dpdpsi(1:NPSI), &
         rho_t(1:NPSI),dpsidrhotor(1:NPSI))
    
    ALLOCATE(g11(1:NPSI,1:NCHI),g12(1:NPSI,1:NCHI),g22(1:NPSI,1:NCHI),g33(1:NPSI,1:NCHI))
    ALLOCATE(B_2d(1:NPSI,1:NCHI),dBdPSI_2d(1:NPSI,1:NCHI),&
         dBdCHI_2d(1:NPSI,1:NCHI),J_2d(1:NPSI,1:NCHI))
    
    ALLOCATE(R(1:NPSI,1:NCHI),Z(1:NPSI,1:NCHI),&
         dPsidr(1:NPSI,1:NCHI),dPsidz(1:NPSI,1:NCHI))

    ! 1d quantities
    CALL getarr(hdf5_ioutgyro, "/data/grid/PSI", PSI)
    CALL getarr(hdf5_ioutgyro, "/data/grid/CHI"  , CHI)
    !
    CALL getarr(hdf5_ioutgyro, "/data/var1d/Rgeom" ,   Rgeom)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/ageom" ,   ageom)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/Volume" ,   V)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/dVdpsi" ,   dVdpsi)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/q", q)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/dqdpsi", dqdpsi)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/rho_tor", rho_t)
    CALL getarr(hdf5_ioutgyro, "/data/var1d/dpsidrhotor",dpsidrhotor)
    IF (pressure_term) THEN
       CALL getarr(hdf5_ioutgyro, "/data/var1d/p", p)
       CALL getarr(hdf5_ioutgyro, "/data/var1d/dpdpsi", dpdpsi)
    END IF



   ! 2-dim quantities

    CALL getarr(hdf5_ioutgyro, "/data/var2d/g11", g11)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/g12", g12)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/g22", g22)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/g33", g33)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/B",   B_2d)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/dBdpsi", dBdPSI_2d)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/dBdchi", dBdchi_2d)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/Jacobian",    J_2d)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/R",   R)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/Z",   Z)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/dPsidR", dPsidr)
    CALL getarr(hdf5_ioutgyro, "/data/var2d/dPsidZ", dPsidz)  

    ! close HDF file
    CALL closef(hdf5_ioutgyro)

  END SUBROUTINE read_CHEASE_input_HDF5

  !-----------------------------------------------------------
  SUBROUTINE CLOSE_GEOM(interpol)
    logical, intent(in):: interpol

    DEALLOCATE(dxdpsi,Rgeom,ageom,q,dqdpsi,PSI,CHI,tmp1d,V,dVdpsi,p,dpdpsi,&
         rho_t,dpsidrhotor)
    DEALLOCATE(g11,g12,g22,g33,B_2d,dBdPSI_2d,dBdCHI_2d,J_2d)
    DEALLOCATE(R,Z,dPsidr,dPsidz)
!    if (interpol) then
       DEALLOCATE(dxdpsi_int,Rgeom_int,ageom_int,q_int,dqdpsi_int,PSI_int,CHI_int,&
            tmp1d_int,V_int,dVdpsi_int,p_int,dpdpsi_int,&
            rho_t_int,dpsidrhotor_int)
       DEALLOCATE(g11_int,g12_int,g22_int,g33_int,B_2d_int,dBdPSI_2d_int,dBdCHI_2d_int,J_2d_int)
       DEALLOCATE(R_int,Z_int,dPsidr_int,dPsidz_int)
 !   endif
  END SUBROUTINE CLOSE_GEOM
  !-----------------------------------------------------------

  !>  set x definition, find center of the simulation domain in CHEASE mesh (NPSI0)
  ! check compatibility between the number of x points and CHEASE mesh
  ! for global simulation
  SUBROUTINE set_xbox(x0,lx,q0,interpol)
    REAL, INTENT(OUT)      :: q0
    REAL, INTENT(INOUT)    :: x0,lx
    INTEGER                :: i
    real                   :: min_diff, check_mesh, delta_x, cpi
    REAL, DIMENSION(:), allocatable:: x_o_a         ! =x/a
    CHARACTER(len=80) :: x_str
    logical:: interpol
    
    cpi=4.*atan(1.)
    if (.not.force_x) then
       if (.not.xy_local) then
          ALLOCATE(x_o_a(1:nx0))
       else
          ALLOCATE(x_o_a(1))
       endif
    else 
       ALLOCATE(x_o_a(1:NPSI))
    end if
   
    NPSI0=NPSI/2
    
    if (.not.xy_local) then
       flux_pos = xoa0 - 1.0e-6       ! - 1e-6 : when using even nx point x0 
       ! might be exactly in the middle of two mesh points
    end if
    
    select case (x_def)
       case('C_psi')
       x_str = 'q0/(r0 B0) Psi'
       x_o_a(:)=q_int(:)/(ageom_int(:)**2*B0EXP)*PSI_int(:)
       dxdpsi_int(:)=q_int(:)/(ageom_int(:)*B0EXP)
    case('arho_t')
       x_str = 'a rho_t'
       x_o_a(:)=rho_t_int(:)/maxval(rho_t)
       dxdpsi_int(:)=1./dpsidrhotor_int(:)
    case('arho_p')
       x_str ='a rho_p'
       x_o_a=sqrt(PSI_int/PSI(NPSI))
       dxdpsi_int(:)=ageom(NPSI)/(2.0*sqrt(PSI_int(:)/PSI(NPSI)))/PSI(NPSI)
    case('arho_v')
       x_str ='a rho_v'
       x_o_a=sqrt(V_int/V(NPSI))
       dxdpsi_int(:)=ageom(NPSI)*dVdpsi_int/(2.0*sqrt(V(NPSI)*V_int(:)))
    case default
       if (mype.eq.0) print*, 'x_def has to be set to C_psi, arho_t, arho_p, arho_v'
       stop
    end select
    min_diff = abs(x_o_a(1)-flux_pos)
    check_mesh=0
    !delta_x=x_o_a(2)-x_o_a(1)
    
    if (.not.xy_local) then
       delta_x=x_o_a(2)-x_o_a(1)
    end if
    
    if (interpol .and. .not.force_x) then
       if (.not.xy_local ) then
          DO i=1,nx0
             IF (abs(x_o_a(i)-flux_pos)<=min_diff) then
                min_diff = abs(x_o_a(i)-flux_pos)
                NPSI0=i
             END IF
             !for the moment we leave out the outer points to avoid the deviations there
             if (i<nx0-5.and.i>5) check_mesh=max(check_mesh,&
                  &abs((x_o_a(i+1)-x_o_a(i))-delta_x))
          END DO
       else if (xy_local) then 
          NPSI0=1
       end if       
    else
       DO i=1,NPSI
          IF (abs(x_o_a(i)-flux_pos)<=min_diff) then
             min_diff = abs(x_o_a(i)-flux_pos)
             NPSI0=i
          END IF
       END DO
    end if
     
    IF (.not.xy_local ) THEN
       if (.not.interpol) then
          ! check that the grid is equidistant in x_def
          check_mesh=check_mesh/delta_x
          IF (check_mesh>1e-2) THEN
             if (mype.eq.0) print*, 'The mesh in the chease file is not equidistant in x, change x_def value', check_mesh
             stop
          END if
       endif
       
       IF (MODULO(nx0,2).eq.1) THEN         
          !x0=x_o_a(NPSI0)
          NPSI1=1
          !lx=(x_o_a(NPSI1+nx0-1)-x_o_a(NPSI1))/rhostar
          ! some ref quantities
          r0=rho_t_int(NPSI0)
          q0=q_int(NPSI0)
          p0=p_int(NPSI0)
          dxdpsi0=dxdpsi_int(NPSI0)
          dqdpsi0=dqdpsi_int(NPSI0)
       ELSE
          !x0=(x_o_a(NPSI0+1)+x_o_a(NPSI0))/2.0
          NPSI1=1
          !lx=(x_o_a(NPSI1+nx0-1)-x_o_a(NPSI1))/rhostar
          ! some ref quantities
          r0=(rho_t_int(NPSI0+1)+rho_t_int(NPSI0))/2.0
          q0=(q_int(NPSI0+1)+q_int(NPSI0))/2.0
          p0=(p_int(NPSI0+1)+p_int(NPSI0))/2.
          dxdpsi0=(dxdpsi_int(NPSI0+1)+dxdpsi_int(NPSI0))/2.0 
          dqdpsi0=(dqdpsi_int(NPSI0+1)+dqdpsi_int(NPSI0))/2.0  
       END IF
       if (mype==0.and.print_ini_msg) print *,'lowest x=',&
            &x_o_a(1),', highest x=',x_o_a(nx0)

    ELSE !local
       r0=rho_t_int(NPSI0)
       q0=q_int(NPSI0)
       p0=p_int(NPSI0)
       dxdpsi0=dxdpsi_int(NPSI0)
       dqdpsi0=dqdpsi_int(NPSI0)
    END IF
   
    IF ((mype.eq.0).and.(print_ini_msg)) THEN
       write(*,'(2A)')      ' x definition : ', x_str 
       write(*,'(A,F8.4)')  ' center of the simulation domain at x/a =', x_o_a(NPSI0)
       IF (.not.xy_local) THEN
          write(*,'(A,F8.4,A)') ' set lx  = ', lx , ' (in rho_ref unit)'
          write(*,'(A,F8.4)')   ' box size lx/a  = ', lx*rhostar 
       END IF
    END IF
    deallocate(x_o_a)
    
  END SUBROUTINE set_xbox

  SUBROUTINE interpolate_geometry(npsi,nx0,nchi,nz0,lx,x0,nx0int,nchiint,interpol)
    integer, intent(in):: npsi, nx0, nchi, nz0,nx0int,nchiint
    real, intent(in):: lx, x0
    integer:: i,k
    real:: dx, dz, cpi
    real,dimension(1:npsi):: x_in
    real,dimension(1:nx0int):: x_out
    real,dimension(1:nchiint):: z_in
    real,dimension(1:nz0):: z_out
    real,dimension(1:npsi,1:nchi+1):: tmparr
    logical:: interpol
    
    
    CPI = 4.*atan(1.)
    !temporary arrays for interpolation
    allocate(dxdpsi_int(1:nx0int),Rgeom_int(1:nx0int),ageom_int(1:nx0int),V_int(1:nx0int),dVdpsi_int(1:nx0int),&
         &q_int(1:nx0int), dqdpsi_int(1:nx0int), psi_int(1:nx0int), chi_int(1:nz0), p_int(1:nx0int),&
         &tmp1d_int(1:nx0int), dpdpsi_int(1:nx0int), rho_t_int(1:nx0int), dpsidrhotor_int(1:nx0int),&
         &g11_int(1:nx0int,1:nz0), g12_int(1:nx0int,1:nz0), g22_int(1:nx0int,1:nz0), g33_int(1:nx0int,1:nz0),&
         &B_2d_int(1:nx0int,1:nz0), dBdpsi_2d_int(1:nx0int,1:nz0),dBdchi_2d_int(1:nx0int,1:nz0),&
         &J_2d_int(1:nx0int,1:nz0), R_int(1:nx0int,1:nz0), Z_int(1:nx0int,1:nz0), dPsidr_int(1:nx0int,1:nz0),&
         &dPsidz_int(1:nx0int,1:nz0))
    select case (x_def)
    case('C_psi')
       x_in=q/(ageom**2*B0EXP)*psi
    case('arho_t')
       x_in=rho_t(:)/maxval(rho_t)
    case('arho_p')
       x_in=sqrt(PSI/PSI(NPSI))
    case('arho_v')
       x_in=sqrt(V/V(NPSI))
    end select
 IF (.not.xy_local) THEN
       dx=lx*rhostar/(nx0-1)
       do i=1,nx0int
          x_out(i)=x0-lx*rhostar/2+(i-1)*dx
       enddo
    ELSE
       if (.not.force_x) then
          x_out=flux_pos
       else 
          x_out=x_in
       end if
    END IF
    !enter this if interpolating in x direction. not force_x ensures that will work for the global
    if (interpol .and. .not.force_x) then
       call lag3interp(q,x_in,npsi,q_int,x_out,nx0int)
       call lag3interp(dxdpsi,x_in,npsi,dxdpsi_int,x_out,nx0int)
       call lag3interp(Rgeom,x_in,npsi,Rgeom_int,x_out,nx0int)
       call lag3interp(ageom,x_in,npsi,ageom_int,x_out,nx0int)
       call lag3interp(V,x_in,npsi,V_int,x_out,nx0int)
       call lag3interp(dVdpsi,x_in,npsi,dVdpsi_int,x_out,nx0int)
       call lag3interp(dqdpsi,x_in,npsi,dqdpsi_int,x_out,nx0int)
       call lag3interp(rho_t,x_in,npsi,rho_t_int,x_out,nx0int)
       call lag3interp(dpsidrhotor,x_in,npsi,dpsidrhotor_int,x_out,nx0int)
       call lag3interp(psi,x_in,npsi,psi_int,x_out,nx0int)
       call lag3interp(p,x_in,npsi,p_int,x_out,nx0int)
       call lag3interp(dpdpsi,x_in,npsi,dpdpsi_int,x_out,nx0int)
    else
       q_int(:)=q(:)
       dxdpsi_int(:)=dxdpsi(:)
       Rgeom_int(:)=Rgeom(:)
       ageom_int(:)=ageom(:)
       V_int(:)=V(:)
       dVdpsi_int(:)=dVdpsi(:)
       dqdpsi_int(:)=dqdpsi(:)
       rho_t_int(:)=rho_t(:)
       dpsidrhotor_int(:)=dpsidrhotor(:)
       if (dpsidrhotor_int(1) .le. 1e-6) then
         dpsidrhotor_int(1)=1e-6 !when not interpolating and local this might be 0.0 in chease output, need to correct 
       endif
       psi_int(:)=psi(:)
       p_int(:)=p(:)
       dpdpsi_int(:)=dpdpsi(:)
    endif
    
    !enter this if interp is true, will interpolate in z direction again def value takes into acount for global case
    if (interpol) then
       if (.not.force_z) then
          dz=2*cpi/nz0
          do k=1,nchi+1
             z_in(k)=(k-1)*2*cpi/nchi
          enddo
          do k=1,nz0
             z_out(k)=(k-1)*dz
          enddo
       else if (force_z) then
          dz=2*cpi/nz0
          do k=1,nchi+1
             z_in(k)=(k-1)*2*cpi/nchi
          enddo
          z_out=z_in(1:NCHI)
       endif
   
       if (.not.force_z) then
          call lag3interp(chi,z_in,nchi,chi_int,z_out,nz0)
       else
          chi_int=chi
       end if
       !we have to create a temporary array here to enforce chi-periodicity
       tmparr(1:npsi,1:nchi)=g11
       tmparr(1:npsi,nchi+1)=g11(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,g11_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=g12
       tmparr(1:npsi,nchi+1)=g12(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,g12_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=g22
       tmparr(1:npsi,nchi+1)=g22(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,g22_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=g33
       tmparr(1:npsi,nchi+1)=g33(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,g33_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=B_2d
       tmparr(1:npsi,nchi+1)=B_2d(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,B_2d_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=dBdpsi_2d
       tmparr(1:npsi,nchi+1)=dBdpsi_2d(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,dBdpsi_2d_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=dBdchi_2d
       tmparr(1:npsi,nchi+1)=dBdchi_2d(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,dBdchi_2d_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=J_2d
       tmparr(1:npsi,nchi+1)=J_2d(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,J_2d_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=R
       tmparr(1:npsi,nchi+1)=R(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,R_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=Z
       tmparr(1:npsi,nchi+1)=Z(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,Z_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=dPsidr
       tmparr(1:npsi,nchi+1)=dPsidr(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,dPsidr_int,x_out,nx0int,z_out,nz0)
       tmparr(1:npsi,1:nchi)=dPsidz
       tmparr(1:npsi,nchi+1)=dPsidz(1:npsi,1)
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,dPsidz_int,x_out,nx0int,z_out,nz0)
    else  
       !we don't have to interpolate
       chi_int(:)=chi(:)
       g11_int(:,:)=g11(:,:)
       g12_int(:,:)=g12(:,:)
       g22_int(:,:)=g22(:,:)
       g33_int(:,:)=g33(:,:)
       B_2d_int(:,:)=B_2d(:,:)
       dBdpsi_2d_int(:,:)=dBdpsi_2d(:,:)
       dBdchi_2d_int(:,:)=dBdchi_2d(:,:)
       J_2d_int(:,:)=J_2d(:,:)
       R_int(:,:)=R(:,:)
       Z_int(:,:)=Z(:,:)
       dPsidr_int(:,:)=dPsidr(:,:)
       dPsidz_int(:,:)=dPsidz(:,:)
    end if

  END SUBROUTINE interpolate_geometry

#endif

END MODULE chease_mod
