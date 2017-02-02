!>This module provides external contributions
!!
!!Provides external contributions, e.g., of the electrostatic
!!potential, and contains functions to add those contributions
!!to the Vlasov equations (we currently consider only weak
!!external terms thus they do not enter the field equations or
!!the Maxwellian)
#include "redef.h"
#include "intrinsic_sizes.h"
MODULE external_contr
  use par_mod
  use discretization_adptv_module
  use coordinates, only: kxmin, xval, vp
  use fourier, only: initialize_fourier_x_1d, to_real_x_1d, finalize_fourier_x_1d,&
       &initialize_fourier_boundary, to_fourier_boundary, to_real_boundary,&
       &finalize_fourier_boundary
  use geometry, only: rhostar, Lref, q0, C_y, minor_r
  use communications, ONLY: MPI_COMPLEX_TYPE,mpi_comm_x
  use profiles_mod, only: write_tempdens_profiles, cref, omegator, omegatorprime
  use file_io, only: get_unit_nr
 
  IMPLICIT NONE

  !add single kx mode
  public:: init_external_contributions, finalize_external_contributions, check_external_contr
  public:: add_phi_ext,build_ext_contr_prefac,add_external_contr, mem_est_external_contr
  public:: phi0_ext, phase_phi_ext, omt0_ext, phase_omt_ext, omn0_ext, phase_omn_ext, &
       kxind_phi_ext, kxind_omt_ext, kxind_omn_ext, apar0_ext, phase_apar_ext, kxind_apar_ext
  ! routines for adaptive grids
  public:: apply_shear_flow_adptv
  !flow shear as described by G. Hammett et al.
  public:: initialize_shear_flow, apply_shear_flow, finalize_shear_flow
  public:: ExBrate, pfsrate, ExB_stime, set_ExBrate, lilo_w_full_omegator_profile
  public:: with_external, set_ext_defaults

  private

  !input parameters
  !amplitudes and phases
  Real :: phi0_ext=0.0D0, phase_phi_ext=0.0D0
  Real :: apar0_ext=0.0D0, phase_apar_ext=0.0D0
  Real :: omt0_ext=0.0D0, phase_omt_ext=0.0D0
  Real :: omn0_ext=0.0D0, phase_omn_ext=0.0D0
  Real :: ExBrate=0.0D0, pfsrate=0.0D0
  Real, dimension(:), allocatable :: pfsrate_arr
  Complex, dimension(:,:), allocatable :: Omega_ExB
  Real :: ExB_stime = 0.0 !time where ExB shear flow is switched on
  
  !kx mode to be initialized 
  !(negative counterpart will be initialized as well!)
  Integer :: kxind_phi_ext=-1, kxind_omt_ext=-1, kxind_omn_ext=-1, kxind_apar_ext=-1

  !prefactors
  Complex, Dimension(:,:), Allocatable, Private :: p_phi_ext
  Complex, Dimension(:,:,:,:), Allocatable, Private :: p_apar_ext
  Complex, Dimension(:,:,:,:,:,:), Allocatable, Private :: p_omt_ext,p_omn_ext
  Integer, Dimension(:), Allocatable, Private :: kxind_ordered

  !Private variables
  Logical:: write_pe
  Logical:: with_external

  Integer:: init_status_ExB=0 !(for Hammett model)
  Real, Dimension(:), Allocatable:: kxExB

  !This flag allows use of a full omegator profile (i.e. nonlocal) with LILO
  Logical :: lilo_w_full_omegator_profile = .false.   
  
CONTAINS

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_external_contr(mem_req_in)
    real:: mem_req_in
    real:: mem_loc

    mem_loc=0.
    if ((apar0_ext.gt.0.0).and.(kxind_apar_ext.gt.0).and.(kxind_apar_ext.lt.nx0/2))&
         mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*ll0*ln0
    if ((phi0_ext.gt.0.0).and.(kxind_phi_ext.gt.0).and.(kxind_phi_ext.lt.nx0/2))&
         mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0
    if ((omt0_ext.gt.0.0).and.(kxind_omt_ext.gt.0).and.(kxind_omt_ext.lt.nx0/2)) &
         mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lijklmn0
    if ((omn0_ext.gt.0.0).and.(kxind_omn_ext.gt.0).and.(kxind_omn_ext.lt.nx0/2)) &
         mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lijklmn0
    if (ExB) then
       if (x_local) then
          mem_loc=mem_loc+SIZE_OF_REAL_MB*lg0
       else
          mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*pi0*pj0
       endif
    endif

    mem_est_external_contr=mem_req_in+mem_loc
  End Function mem_est_external_contr

  Subroutine check_external_contr
    with_phi_ext=((phi0_ext.gt.0.0).and.(kxind_phi_ext.gt.0).and.&
         &(kxind_phi_ext.lt.nx0/2))
    with_apar_ext=((apar0_ext.gt.0.0).and.(kxind_apar_ext.gt.0).and.&
         &(kxind_apar_ext.lt.nx0/2))
    with_omt_ext=((omt0_ext.gt.0.0).and.(kxind_omt_ext.gt.0).and.&
         &(kxind_omt_ext.lt.nx0/2))
    with_omn_ext=((omn0_ext.gt.0.0).and.(kxind_omn_ext.gt.0).and.&
         &(kxind_omn_ext.lt.nx0/2))
    
    with_external=with_phi_ext.or.with_omt_ext.or.with_omn_ext.or.with_apar_ext

  End Subroutine check_external_contr

  Subroutine set_ext_defaults
    phi0_ext=0.0D0
    phase_phi_ext=0.0D0
    apar0_ext=0.0D0
    phase_apar_ext=0.0D0
    omt0_ext=0.0D0
    phase_omt_ext=0.0D0
    omn0_ext=0.0D0
    phase_omn_ext=0.0D0
    ExBrate=0.0D0
    pfsrate=0.0D0
    ExB_stime = 0.0 !time where ExB shear flow is switched on
  End Subroutine set_ext_defaults


  Subroutine init_external_contributions
    !COMPLEX, DIMENSION(0:nx0-1) :: tmp, phi_ext, omt_ext, omn_ext
    COMPLEX, DIMENSION(li1:li2) :: tmp, phi_ext, omt_ext, omn_ext, apar_ext, tmp2
    COMPLEX,DIMENSION(0:nx0-1) :: global_x_temp
    INTEGER:: i,j,k,l,m,n, relevant_pex, ierr
    real::vT

    write_pe = ((mype==0).AND.(print_ini_msg))
    
    if (with_external) then
       if (x_local) then !transformation from Fortran FFT order to natural order
          !(required in convolution)
          allocate(kxind_ordered(0:nx0-1))
          if (nx0.gt.1) then
             do i=0,hkx+evenx
                kxind_ordered(hkx+(i)) = i
             enddo
             do i=lkx,nx0-1
                kxind_ordered(i-lkx)=i
             enddo
          else
             kxind_ordered = li1
          endif
       else !initialize fourier transforms
          !PRINT*,"Calling initialize_fourier_x_1d from external_contr.",li1,li2
          Call initialize_fourier_x_1d
       endif
       
       !initialize external phi
       if (with_phi_ext) then
          ALLOCATE(p_phi_ext(li1:li2,lj1:lj2))
          if (write_pe) then
             WRITE(*,'(A)') 'Using external modifications of electrostatic potential with'
             WRITE(*,'(A,ES10.3,A,I3,A,ES10.3)') 'amplitude=', phi0_ext,', kxind=', &
                  & kxind_phi_ext,', phase=', phase_phi_ext
          endif
          
          if (x_local) then
             phi_ext = CMPLX(0.0,0.0)       
             phi_ext(li1+kxind_phi_ext) = 0.5*phi0_ext*CMPLX(0.0,-1.0)*EXP(imag*phase_phi_ext)
             phi_ext(li2-kxind_phi_ext+1) = 0.5*phi0_ext*CMPLX(0.0,1.0)*EXP(-imag*phase_phi_ext)
             !p_phi_ext = - d/dx phi_ext * d/dy = -imag*kx*phi_ext* imag*ky =
             !          = kx*phi_ext *ky
             do j=lj1,lj2
                p_phi_ext(li1:li2,j) = kx(li1:li2)*phi_ext(li1:li2)*ky(j)
             enddo
          else
             phi_ext = CMPLX(0.0,0.0)
             relevant_pex = kxind_phi_ext/li0
             IF (my_pex.EQ.relevant_pex) THEN
                phi_ext(kxind_phi_ext) = 0.5*phi0_ext*CMPLX(0.0,-1.0)*EXP(imag*phase_phi_ext)
             END IF
             relevant_pex = (nx0-kxind_phi_ext)/li0
             IF (my_pex.EQ.relevant_pex) THEN
                phi_ext(nx0-kxind_phi_ext) = 0.5*phi0_ext*CMPLX(0.0,1.0)*EXP(-imag*phase_phi_ext)
             END IF
             
             tmp = -imag*kx(li1:li2)*phi_ext
             Call to_real_x_1d(tmp,tmp2)
             do j=lj1,lj2
                p_phi_ext(:,j)=tmp2(:)
             enddo
          endif
       endif
       
       !initialize external apar
       if (with_apar_ext) then
          ALLOCATE(p_apar_ext(li1:li2,lj1:lj2,ll1:ll2,ln1:ln2))
          call set_vp_coordinate_vars  !repeated in vel_space, but needed here, too

          if (write_pe) then
             WRITE(*,'(A)') 'Using external modifications of magnetic potential with'
             WRITE(*,'(A,ES10.3,A,I3,A,ES10.3)') 'amplitude=', apar0_ext,', kxind=', &
                  & kxind_apar_ext,', phase=', phase_apar_ext
          endif
          
          if (x_local) then
             apar_ext = CMPLX(0.0,0.0)       
             apar_ext(li1+kxind_apar_ext) = 0.5*apar0_ext*CMPLX(0.0,-1.0)*EXP(imag*phase_apar_ext)
             apar_ext(li2-kxind_apar_ext+1) = 0.5*apar0_ext*CMPLX(0.0,1.0)*EXP(-imag*phase_apar_ext)
             do n=ln1,ln2
                vT = sqrt(2.0*spec(n)%temp/spec(n)%mass)
                do l=ll1,ll2
                   do j=lj1,lj2 
                      !sign change with respect to phi due to definition of chi
                      p_apar_ext(li1:li2,j,l,n) = -kx(li1:li2)*apar_ext(li1:li2)*ky(j)*vp(l)*vT
                   enddo
                enddo
             enddo
          else
             apar_ext = CMPLX(0.0,0.0)
             relevant_pex = kxind_apar_ext/li0
             IF (my_pex.EQ.relevant_pex) THEN
                apar_ext(kxind_phi_ext) = 0.5*apar0_ext*CMPLX(0.0,-1.0)*EXP(imag*phase_apar_ext)
             END IF
             relevant_pex = (nx0-kxind_apar_ext)/li0
             IF (my_pex.EQ.relevant_pex) THEN
                apar_ext(nx0-kxind_apar_ext) = 0.5*apar0_ext*CMPLX(0.0,1.0)*EXP(-imag*phase_apar_ext)
             END IF
             
             tmp = imag*kx(li1:li2)*apar_ext !sign change with respect to phi due to definition of chi
             Call to_real_x_1d(tmp,tmp2)
             do n=ln1,ln2
                vT = sqrt(2.0*spec(n)%temp/spec(n)%mass)
                do l=ll1,ll2
                   do j=lj1,lj2
                      p_apar_ext(:,j,l,n)=tmp2*vp(l)*vT
                   enddo
                enddo
             enddo
          endif
       endif
       
       !initialize external temperature variation (required only for tertiary inst. analysis)
       if (with_omt_ext) then
          if (write_pe) then
             WRITE(*,'(A)') 'Using external temperature gradient modifications with'
             WRITE(*,'(A,ES10.3,A,I3,A,ES10.3)') 'amplitude=', omt0_ext,', kxind=', &
                  & kxind_omt_ext,', phase=', phase_omt_ext
          endif
          
          omt_ext = cmplx(0.0,0.0)
          
          if (x_local) then
             omt_ext(li1+kxind_omt_ext) = 0.5*omt0_ext*CMPLX(0.0,-1.0)*EXP(imag*phase_omt_ext)
             omt_ext(li2-kxind_omt_ext+1) = 0.5*omt0_ext*CMPLX(0.0,1.0)*EXP(-imag*phase_omt_ext)
             ALLOCATE(p_omt_ext(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
             DO n=ln1,ln2
                do m=lm1,lm2
                   do l=ll1,ll2
                      do k=lk1,lk2
                         do j=lj1,lj2
                            p_omt_ext(:,j,k,l,m,n) = omt_ext(li1:li2)
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          else
             relevant_pex = kxind_omt_ext/li0
             IF (my_pex.EQ.relevant_pex) THEN
                omt_ext(li1+kxind_omt_ext) = 0.5*omt0_ext*CMPLX(0.0,-1.0)*EXP(imag*phase_omt_ext)
             END IF
             relevant_pex = (nx0-kxind_omt_ext)/li0
             IF (my_pex.EQ.relevant_pex) THEN
                omt_ext(li2-kxind_omt_ext+1) = 0.5*omt0_ext*CMPLX(0.0,1.0)*EXP(-imag*phase_omt_ext)
             END IF
             tmp = omt_ext
             Call to_real_x_1d(tmp,omt_ext)
             CALL mpi_allgather(omt_ext,li0,MPI_COMPLEX_TYPE,global_x_temp,li0,MPI_COMPLEX_TYPE,mpi_comm_x,ierr)
             do n=0, n_spec-1 
                spec(n)%omt_prof(0:nx0-1) = spec(n)%omt_prof(0:nx0-1) &
                     & + global_x_temp(0:nx0-1)
             enddo
          endif
       endif
       
       !initialize external density variation (required only for tertiary inst. analysis)
       if (with_omn_ext) then
          if (write_pe) then
             WRITE(*,'(A)') 'Using external density gradient modifications with'
             WRITE(*,'(A,ES10.3,A,I3,A,ES10.3)') 'amplitude=', omn0_ext,', kxind=', & 
                  & kxind_omn_ext, ', phase=', phase_omn_ext
          endif
          
          omn_ext = cmplx(0.0,0.0)
          
          if (x_local) then
             omn_ext(kxind_omn_ext) = 0.5*omn0_ext*CMPLX(0.0,-1.0)*&
                  & EXP(imag*phase_omn_ext)
             omn_ext(nx0-kxind_omn_ext) = 0.5*omn0_ext*CMPLX(0.0,1.0)*&
                  & exp(-imag*phase_omn_ext)
             allocate(p_omn_ext(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))   
             do n=ln1,ln2
                do m=lm1,lm2
                   do l=ll1,ll2
                      do k=lk1,lk2
                         do j=lj1,lj2
                            p_omn_ext(:,j,k,l,m,n) = omn_ext(li1:li2)
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          else
             relevant_pex = kxind_omn_ext/li0
             IF (my_pex.EQ.relevant_pex) THEN
                omn_ext(kxind_omn_ext) = 0.5*omn0_ext*CMPLX(0.0,-1.0)*&
                     & EXP(imag*phase_omn_ext)
             END IF
             relevant_pex = (nx0-kxind_omn_ext)/li0
             IF (my_pex.EQ.relevant_pex) THEN
                omn_ext(nx0-kxind_omn_ext) = 0.5*omn0_ext*CMPLX(0.0,1.0)*&
                     & EXP(-imag*phase_omn_ext)
             END IF
             tmp = omn_ext
             Call to_real_x_1d(tmp,omn_ext)
             CALL mpi_allgather(omn_ext,li0,MPI_COMPLEX_TYPE,&
                  &global_x_temp,li0,MPI_COMPLEX_TYPE,mpi_comm_x,ierr)
             do n=0, n_spec-1 
                spec(n)%omn_prof(0:nx0-1) = spec(n)%omn_prof(0:nx0-1) &
                     & + global_x_temp
             enddo
          endif
       endif
       
       if(.not.x_local) call finalize_fourier_x_1d
       
       if ((write_pe).and.(.not.x_local)) then
          do n=0,n_spec-1
             Call write_tempdens_profiles(spec(n),.false.)
          end do
       endif
    endif


    if (ExB) then
       if (.not.x_local) then
          allocate(Omega_ExB(pi1:pi2,lj1:lj2))
!          Omega_ExB = 0.0
       endif
       allocate(pfsrate_arr(pi1:pi2))
       call set_ExBrate
    endif

  end subroutine init_external_contributions


  !> Adds the vel. space prefactors to p_omt_ext & p_omn_ext
  !! \param fac_ind index of factor (1 = p_omt_ext, 2 = p_omn_ext)
  !! \param k,l,m,n currently addressed element in remaining dimensions
  !! \param vel_fac velocity space prefactor to be added to p_*_ext
  Subroutine build_ext_contr_prefac(fac_ind,vel_fac,k,l,m,n)
    integer, intent(in) :: fac_ind, k, l, m, n
    real, intent(in)    :: vel_fac
    integer :: j

    If (x_local) then
       If (fac_ind.eq.1) then !external temperature contributions
          Do j=lj1,lj2
             p_omt_ext(li1:li2,j,k,l,m,n) = p_omt_ext(li1:li2,j,k,l,m,n)*&
                  &vel_fac*imag*ky(j)
          End Do
       Elseif (fac_ind.eq.2) then !external density contributions
          Do j=lj1,lj2
             p_omn_ext(li1:li2,j,k,l,m,n) = p_omn_ext(li1:li2,j,k,l,m,n)*&
                  &vel_fac*imag*ky(j)
          End Do
       Else
          STOP 'Unknown fac_ind in build_ext_contr_prefac'
       End If
    End If

  end Subroutine build_ext_contr_prefac


  !> Adds extra (radial) equilibrium contributions to 
  !! (local) right hand side of the Vlasov equation.
  !! In the local code this is done by a convolution 
  !! with one certain mode (and its negative counterpart) 
  !! which is given by the index kxind_ext.
  !! In the global code the external contributions is
  !! absorbed in the usual prefactors (e.g. pdg1dy)
  !! \param localrhs localrhs as two-dimensional array
  !! \param fac_ind index of factor (0 = p_phi_ext, 1 = p_omt_ext, 2 = p_omn_ext,3 = p_apar_ext)
  !! \param g1_or_chi first two dimensions of chi or g1
  !! \param kxind_ext kx mode of external equilibrium contributio
  !! \param k,l,m,n indices in remaining dimensions
  Subroutine add_external_contr(localrhs, fac_ind, g1_or_chi, kxind_ext,k,l,m,n)
    Complex,Dimension(li1:li2,lj1:lj2), intent(inout) :: localrhs
    Complex,Dimension(li1:li2,lj1:lj2), intent(in) :: g1_or_chi
    Integer, intent(in) :: fac_ind, kxind_ext, k, l, m, n

    Complex, Dimension(li1:li2,lj1:lj2) :: prefac_ext
    Integer :: i, j

    If (fac_ind.eq.0) then
       prefac_ext = p_phi_ext(li1:li2,lj1:lj2)
    elseif (fac_ind.eq.1) then
       prefac_ext = p_omt_ext(li1:li2,lj1:lj2,k,l,m,n)
    elseif (fac_ind.eq.2) then
       prefac_ext = p_omn_ext(li1:li2,lj1:lj2,k,l,m,n)
    elseif (fac_ind.eq.3) then
       prefac_ext = p_apar_ext(li1:li2,lj1:lj2,l,n)
    else
       STOP 'Unknown fac_ind in add_external_contr'
    endif

    do j=lj1,lj2
       do i=li2-kxind_ext,li2-evenx
          localrhs(kxind_ordered(i),j) = localrhs(kxind_ordered(i),j) + &
               &prefac_ext(kxind_ext,j)*g1_or_chi(kxind_ordered(i-kxind_ext),j)
       enddo
       do i=li1+kxind_ext,li2-kxind_ext
          localrhs(kxind_ordered(i),j) = localrhs(kxind_ordered(i),j) + &
               & prefac_ext(kxind_ext,j)*g1_or_chi(kxind_ordered(i-kxind_ext),j)
       enddo
       do i=li1+kxind_ext,li2-kxind_ext-evenx
          localrhs(kxind_ordered(i),j) = localrhs(kxind_ordered(i),j) + &
               &prefac_ext(li2-kxind_ext+1,j)*g1_or_chi(kxind_ordered(i+kxind_ext),j)
       enddo
       do i=li1,li1+kxind_ext
          localrhs(kxind_ordered(i),j) = localrhs(kxind_ordered(i),j) + &
               & prefac_ext(li2-kxind_ext+1,j)*g1_or_chi(kxind_ordered(i+kxind_ext),j)
       enddo
    enddo

  End Subroutine add_external_contr


  !> Add external phi and/or apar contribution to right hand side in global code
  Subroutine add_phi_ext(m_pdg1dy,p_fac,l,n)
    Real, Dimension(li1:li2),intent(inout) :: m_pdg1dy !<pdg1dy prefactor
    Real, Dimension(pi1:pi2),intent(in) :: p_fac !<additional factor in front of p_phi_ext
    integer,intent(in)::l,n
    if (with_phi_ext)&
         &m_pdg1dy = m_pdg1dy + p_fac*real(p_phi_ext(li1:li2,lj1))
    if (with_apar_ext)&
         &m_pdg1dy = m_pdg1dy + p_fac*real(p_apar_ext(li1:li2,lj1,l,n))
  End Subroutine add_phi_ext


  !>Deallocates all array used for external contributions
  Subroutine finalize_external_contributions
    if (with_external) then
       if (allocated(kxind_ordered)) deallocate(kxind_ordered)
       if (with_phi_ext) deallocate(p_phi_ext)
       if (with_apar_ext) deallocate(p_apar_ext)
       if (allocated(p_omt_ext)) deallocate(p_omt_ext)   
       if (allocated(p_omn_ext)) deallocate(p_omn_ext)
    endif
    if (ExB) then
       if (.not.x_local) deallocate(Omega_ExB)
       deallocate(pfsrate_arr)
    endif
  end Subroutine finalize_external_contributions

  !>Determine ExBrate from VROT profile
  !!Needs to be called before prefactor computation 
  !!as pfsrate might change
  subroutine set_ExBrate
    integer :: ix0, j

    if (ExBrate.eq.-1111.0.or.lilo_w_full_omegator_profile) then
       if (abs(cref).lt.epsilon(abs(cref))) stop &
            &'Cannot evaluate ExBrate without properly set reference values'
       if (.not.allocated(omegator)) stop &
            &'No toroidal angular velocity available'

       !use the center value even in x_global cases as reference
       ix0 = (pi1gl+pi2gl)/2
       if (x_local.or.lilo) ExBrate = -omegatorprime(ix0)*&
            &xval_a(ix0)/q0*Lref/cref

       !cannot do this nonlocally at the moment due to auto-parallelization
       !where ExBrate will be <>-1111.0 at the next call and will therefore
       !break the Omega_ExB defintion

       if ((.not.x_local).and.((.not.lilo).or.lilo_w_full_omegator_profile)) then
          do j=lj1,lj2
!             Omega_ExB(pi1:pi2,j) = -imag*(omegator(pi1:pi2)-omegator(ix0))*C_y(pi1:pi2)*ky(j)*&
!                  &Lref/cref*1.0/(rhostar*minor_r)
             Omega_ExB(pi1:pi2,j) = -imag*omegator(pi1:pi2)*C_y(pi1:pi2)*ky(j)*&
                  &Lref/cref*1.0/(rhostar*minor_r)
          enddo
          ! if(mype == 0) write(*,*) "Omega_ExB(:)",Omega_ExB
          if (pfsrate.eq.-1111.0) then
             if (mype==0) print*, 'parallel flow shear not implemented for x global'
             pfsrate = 0.0
             pfsrate_arr = 0.0
          endif
       endif
    endif

    if (x_local.or.(lilo.and..not.lilo_w_full_omegator_profile)) then
       if (pfsrate.eq.-1111.0) pfsrate=ExBrate
       pfsrate_arr = pfsrate
       if (lilo) then
          do j=lj1,lj2 
             Omega_ExB(pi1:pi2,j) = -imag*ky(j)*ExBrate*xval(pi1:pi2)
          enddo
       endif
    endif

  end subroutine set_ExBrate

  !>Initializes the Hammett ExB shear model
  !!Called in time_scheme.F90 as this model is only available in 
  !!certain (explicit) solvers
  subroutine initialize_shear_flow
    !kxExB is used to monitor the time evolution of the 
    !kx grid for the Hammett ExB shear model
    if(ExB) then
       if (init_status_ExB==0) then
          if (x_local) allocate(kxExB(lg1:lg2))
          init_status_ExB=1
       endif
       if (x_local) kxExB=0.
    end if
       
  end subroutine initialize_shear_flow

  !> Models a radially constant ExB shearing rate by a discrete shearing of the simulation domain in time 
  !! (the distribution function is shifted by one i-index). In the nonlocal code, a phase factor is applied
  !! at every timestep
  subroutine apply_shear_flow(p_g1,dt)
    complex, dimension(li1:li2,lj1:lj2,1:lklmn0), intent(inout) :: p_g1 
    real, intent(in):: dt
    complex, dimension(li1:li2,lj1:lj2) :: tmparr
    integer:: i,j, ljs, n_shift, shifts, klmn

    if (x_local) then
       if (yx_order) then
          !time dependent kx grid keeps track of the effect of ExB shear
          kxExB(li1:li2)=kxExB(li1:li2)-ky(li1:li2)*ExBrate*dt
          if (any(abs(kxExB(li1:li2)-kx(0)).gt.0.5*kxmin)) then
             if (.not.y_local) then
                do klmn=1,lklmn0
                   call to_fourier_boundary(p_g1(:,:,klmn),tmparr)
                   call apply_kx_shift_yx(tmparr,(klmn.eq.lklmn0))
                   call to_real_boundary(tmparr,p_g1(:,:,klmn))
                enddo
             else
                do klmn=1,lklmn0
                   call apply_kx_shift_yx(p_g1(:,:,klmn),(klmn.eq.lklmn0))
                enddo
             endif
          endif
       else 
          do j=lj1,lj2
             !time dependent kx grid keeps track of the effect of ExB shear
             kxExB(j)=kxExB(j)-ky(j)*ExBrate*dt
             shifts=int(abs((kxExB(j)-kx(0))/(0.5*kxmin)))
             !if kxExB(j) deviates by more than half a grid spacing 
             !in kx from the initial grid, copy g_1 to the next kx index
             if ((kxExB(j)-kx(0)).lt.-0.5*kxmin) then
                !downward shift
                do n_shift=1,shifts
                   do i=lkx,nx0-1-1
                      p_g1(i,j,:)=p_g1(i+1,j,:)
                   enddo
                   p_g1(nx0-1,j,:)=p_g1(0,j,:)
                   do i=0,hkx-1
                      p_g1(i,j,:)=p_g1(i+1,j,:)
                   enddo
                   do i=hkx,hkx+evenx
                      p_g1(i,j,:)=0.
                   enddo
                   kxExB(j)=kxExB(j)+kxmin
                enddo
             elseif ((kxExB(j)-kx(0)).gt.0.5*kxmin) then
                !upward shift
                do n_shift=1,shifts
                   p_g1(hkx+evenx,j,:)=0.
                   do i=hkx,1,-1
                      p_g1(i,j,:)=p_g1(i-1,j,:)
                   enddo
                   p_g1(0,j,:)=p_g1(nx0-1,j,:)
                   do i=nx0-1,lkx+1,-1
                      p_g1(i,j,:)=p_g1(i-1,j,:)
                   enddo
                   p_g1(lkx,j,:)=0.
                   kxExB(j)=kxExB(j)-kxmin
                enddo
             endif
          enddo
       end if
    else ! not x_local
       ljs=lj1
       if (p_has_0_mode) ljs = lj1+1
       do j=ljs,lj2
          do i=li1,li2
             p_g1(i,j,:)=p_g1(i,j,:)*exp(Omega_ExB(i,j)*dt)
          enddo
       enddo
    endif

  end subroutine apply_shear_flow

  subroutine apply_kx_shift_yx(g2d,last)
    complex, dimension(li1:li2,lj1:lj2), intent(inout) :: g2d 
    logical, intent(in) :: last
    integer:: i,j,n_shift, shifts, minus_y
    
    if (y_local) then
       do i=li1,li2
          !time dependent kx grid keeps track of the effect of ExB shear
          shifts=int(abs((kxExB(i)-kx(0))/(0.5*kxmin)))
          !if kxExB(j) deviates by more than half a grid spacing 
          !in kx from the initial grid, copy g_1 to the next kx index
          if ((kxExB(i)-kx(0)).lt.-0.5*kxmin) then
             !downward shift
             do n_shift=1,shifts
                do j=lkx,nx0-2
                   g2d(i,j)=g2d(i,j+1)
                enddo
                g2d(i,nx0-1)=g2d(i,0)
                do j=0,hkx-1
                   g2d(i,j)=g2d(i,j+1)
                enddo
                do j=hkx,hkx+evenx
                   g2d(i,j)=0.
                enddo
                if (last) kxExB(i)=kxExB(i)+kxmin
             enddo
          elseif ((kxExB(i)-kx(0)).gt.0.5*kxmin) then
             !upward shift
             do n_shift=1,shifts
                g2d(i,hkx+evenx)=0.
                do j=hkx,1,-1
                   g2d(i,j)=g2d(i,j-1)
                enddo
                g2d(i,0)=g2d(i,nx0-1)
                do j=nx0-1,lkx+1,-1
                   g2d(i,j)=g2d(i,j-1)
                enddo
                g2d(i,lkx)=0.
                if (last) kxExB(i)=kxExB(i)-kxmin
             enddo
          endif
       enddo
    else !y global
       do i=li1,li2
          shifts=int(abs((kxExB(i)-kx(0))/(0.5*kxmin)))
          !if kxExB(j) deviates by more than half a grid spacing 
          !in kx from the initial grid, copy g_1 to the next kx index
          if ((kxExB(i)-kx(0)).lt.-0.5*kxmin) then
             !downward shift
             do n_shift=1,shifts
                !shift down
                do j=0,nx0-2
                   g2d(i,j)=g2d(i,j+1)
                enddo
                !set highest kx mode to zero
                g2d(i,nx0-1)=0.
                if (last) kxExB(i)=kxExB(i)+kxmin
             enddo             
          elseif ((kxExB(i)-kx(0)).gt.0.5*kxmin) then
             !upward shift
             do n_shift=1,shifts
                do j=lj2,lj1+1,-1
                   g2d(i,j)=g2d(i,j-1)
                enddo
                minus_y=mod(2*nky0-i,nky0)
                if (i.le.(nky0-1)/2) then
                   !fill kx=0 zero mode using reality condition 
                   !and next larger kx if ky>0
                   g2d(i,lj1)=conjg(g2d(minus_y,lj1+1))
                else
                   !fill kx=0 zero mode by copying positive ky value at
                   !kx=0 as all positive ky values have already been shifted
                   !(FFT order: ky=0,...,kmax,-kmin,...,-1)
                   g2d(i,lj1)=conjg(g2d(minus_y,lj1))
                endif
                if (last) kxExB(i)=kxExB(i)-kxmin
             enddo
          endif
       enddo
    endif

  end subroutine apply_kx_shift_yx

  !>Deallocates the variables used for the Hammett ExB shear model
  subroutine finalize_shear_flow
    integer :: omega_handle, i, j
    
    if (ExB.and.init_status_ExB==1) then 
       if (x_local) deallocate (kxExB)
       init_status_ExB=0
    endif

    if (ExB .and. mype == 0 .and. lilo) then
        call get_unit_nr(omega_handle)
        open(unit = omega_handle, file = trim(diagdir)//'/Omega_ExB_final'&
             //trim(file_extension), form = 'formatted')
        write(omega_handle,*) "# xval_a, xval, real(Omega_ExB), aimag(Omega_ExB)"
        do j = lj1, lj2
         write(omega_handle,*) "kyind",j
         do i = pi1,pi2
           write(omega_handle,"(4G12.4)") xval_a(i),xval(i),real(Omega_ExB(i,j)),aimag(Omega_ExB(i,j))
         end do
        end do
        close(omega_handle)
    end if


  end subroutine finalize_shear_flow


  ! --- modified routines for adaptive grids -----

  subroutine apply_shear_flow_adptv(dt)
    real, intent(in):: dt
    integer:: i,j, ljs, l1, l2, m1, m2

    if (x_local) then
       stop "no adaptive grids for x local version!"
    else ! not x_local
       ljs=lj1
       if (p_has_0_mode) ljs = lj1+1
       do i=li1,li2
          do j=ljs,lj2
             l1 = ll1_adp(i)
             l2 = ll2_adp(i)
             m1 = lm1_adp(i)
             m2 = lm2_adp(i)
             g_1(i,j,:,l1:l2,m1:m2,:)=g_1(i,j,:,l1:l2,m1:m2,:)*&
                  &exp(Omega_ExB(i,j)*dt)
          enddo
       enddo
    endif

  end subroutine apply_shear_flow_adptv

END MODULE external_contr
