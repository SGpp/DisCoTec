#include "redef.h"
#include "intrinsic_sizes.h"
!>Local variant of the gyroaverage module
!!
!!Here, we just multiply with jfac and go back to the calling function.
MODULE gyro_average_ff_mod
  USE discretization
  USE coordinates
  USE par_in, only: spec
  USE par_other,only: n_fields
  USE aux_func, ONLY: Bessel_J1_over_x, j0
  USE boundaries, only: exchange_z, ubexc, pb_phase
  USE geometry, only: geom

  IMPLICIT NONE
  PUBLIC:: gyro_average_ff, initialize_gyro_average_ff, jfac, get_jfac, &
       & I1_factor, gyro_average_ff_bpar, mem_est_gyro_average_ff, &
       & finalize_gyro_average_ff
  PRIVATE

  ! jfac is declared here, as it belongs to the gyro-averaging
  ! it is not private to be used also directly in block1.F90
  ! for performance reasons.
  REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: jfac, I1_factor

  INTERFACE gyro_average_ff
     MODULE PROCEDURE gyro_average_scalar,gyro_average_1D, gyro_average_2D,&
          &gyro_average_3D, gyro_average_4D, gyro_average_5D, gyro_average_fields
  END INTERFACE

CONTAINS

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_gyro_average_ff(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !jfac
    mem_loc=SIZE_OF_REAL_MB*li0*lj0*lz0*lm0*ln0
    !jfac_cmplx
    mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*li0*lj0*lz0*lm0*ln0
    IF (n_fields .GT. 2) THEN
      !ifac
      mem_loc=SIZE_OF_REAL_MB*li0*lj0*lz0*lm0*ln0
      !ifac_cmplx
      mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*li0*lj0*lz0*lm0*ln0
    END IF

    mem_est_gyro_average_ff=mem_req_in+mem_loc
  End Function mem_est_gyro_average_ff 


  subroutine initialize_gyro_average_ff
    integer:: i,j,k,m,n
    real:: k_perp2, jfactor, ifactor, rho
    complex, dimension(:,:,:,:,:), ALLOCATABLE :: jfac_cmplx, ifac_cmplx

    ALLOCATE(jfac(li1:li2, lj1:lj2, lbz:ubz, lm1:lm2, ln1:ln2))
    ALLOCATE(jfac_cmplx(li1:li2, lj1:lj2, lbz:ubz, lm1:lm2, ln1:ln2))
    IF (n_fields .GT. 2) THEN
      ALLOCATE(I1_factor(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2))
      ALLOCATE(ifac_cmplx(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2))
    END IF

    do n=ln1,ln2
       do m = lm1,lm2
          do k=lk1,lk2
             rho = sqrt(2.0*spec(n)%mass*spec(n)%temp*mu(m)/&
                  &(spec(n)%charge**2*geom%Bfield(pi1,pj1,k)))
             do j = lj1,lj2
                do i = li1,li2
                   k_perp2 = geom%gii(pi1,pj1,k)*ki(i)**2 + &
                        & 2.0*geom%gij(pi1,pj1,k)*ki(i)*kj(j) + &
                        & geom%gjj(pi1,pj1,k)*kj(j)**2
                   jfactor = j0(sqrt(k_perp2)*rho)
                   jfac_cmplx(i,j,k,m,n)=cmplx(jfactor,0.0)
                   IF (n_fields .GT. 2) THEN
                     ifactor = 2.0 * Bessel_J1_over_x(SQRT(k_perp2)*rho)
                     ifac_cmplx(i,j,k,m,n) = CMPLX(ifactor,0.0)
                   END IF
                end do
             end do
          end do
       end do
    end do

    if (evenx.eq.1) then
       !the unphysical highest positive kx-mode is set to zero
       if (yx_order) then
          jfac_cmplx(:,hkx+1,:,:,:)=cmplx(0.)
          IF (n_fields .GT. 2) ifac_cmplx(:,hkx+1,:,:,:) = CMPLX(0.0)
       else
          jfac_cmplx(hkx+1,:,:,:,:)=cmplx(0.)
          IF (n_fields .GT. 2) ifac_cmplx(hkx+1,:,:,:,:) = CMPLX(0.0)
       end if
    endif
    
    do n=ln1,ln2
       do m=lm1,lm2
          call exchange_z(jfac_cmplx(:,:,:,m,n))
          IF (n_fields .GT. 2) CALL exchange_z(ifac_cmplx(:,:,:,m,n))
       enddo
    enddo

    call remove_pb_phase(jfac_cmplx)
    IF (n_fields .GT. 2) call remove_pb_phase(ifac_cmplx)

    jfac=real(jfac_cmplx)
    DEALLOCATE(jfac_cmplx)
    IF (n_fields .GT. 2) THEN
      I1_factor = REAL(ifac_cmplx)
      DEALLOCATE(ifac_cmplx)
    END IF

  end subroutine initialize_gyro_average_ff

  !! the phase factor of exchange_z is *removed*, 
  !! to retain it in the product jfac*phi  or ifac*bpar
  subroutine remove_pb_phase(fac_cmplx)
    complex,intent(inout)::fac_cmplx(li1:li2, lj1:lj2, lbz:ubz, lm1:lm2, ln1:ln2)
    integer::i,j,k,m,n
    if (ubexc.ge.lj1) then
       do n=ln1,ln2
          do m=lm1,lm2
             do k=1,nzb        
                if (yx_order) then
                   do j=lj1,lj2
                      do i=li1,ubexc
                         if (my_pez.eq.0) fac_cmplx(i,j,lbz+k-1,m,n)= &
                              fac_cmplx(i,j,lbz+k-1,m,n)*pb_phase(i)
                         if (my_pez.eq.n_procs_z-1) fac_cmplx(i,j,lk2+k,m,n)= &
                              fac_cmplx(i,j,lk2+k,m,n)*conjg(pb_phase(i))
                      enddo
                   end do
                else
                   do j=lj1,ubexc
                      if (my_pez.eq.0) fac_cmplx(:,j,lbz+k-1,m,n)= &
                           fac_cmplx(:,j,lbz+k-1,m,n)*pb_phase(j)
                      if (my_pez.eq.n_procs_z-1) fac_cmplx(:,j,lk2+k,m,n)= &
                           fac_cmplx(:,j,lk2+k,m,n)*conjg(pb_phase(j))
                   enddo
                end if
             enddo
          enddo
       enddo
    endif
  end subroutine remove_pb_phase
    

SUBROUTINE gyro_average_scalar(func,barfunc,i,j,k,m,n)
  COMPLEX, intent(IN) :: func
  COMPLEX, intent(OUT) :: barfunc
  INTEGER :: i,j,k,m,n

  barfunc = jfac(i,j,k,m,n)*func
END SUBROUTINE gyro_average_scalar

SUBROUTINE gyro_average_1D(func, barfunc,j,k,m,n)
  COMPLEX, DIMENSION(li1:li2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2), INTENT(OUT) :: barfunc
  INTEGER :: j,k,m,n

  barfunc = jfac(:,j,k,m,n)*func(li1:li2)

END SUBROUTINE gyro_average_1D

SUBROUTINE gyro_average_2D(func, barfunc,k,m,n)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(OUT) :: barfunc
  INTEGER :: k,m,n

  barfunc = jfac(:,:,k,m,n)*func
  
END SUBROUTINE gyro_average_2D

SUBROUTINE gyro_average_ff_bpar(func,barfunc,k,m,n)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(OUT) :: barfunc
  INTEGER :: k, m, n
  
  barfunc = I1_factor(:,:,k,m,n) * func
END SUBROUTINE gyro_average_ff_bpar

SUBROUTINE gyro_average_3D(func, barfunc, m,n)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(OUT) :: barfunc
  INTEGER :: m,n

  barfunc = jfac(:,:,:,m,n)*func
END SUBROUTINE gyro_average_3D

SUBROUTINE gyro_average_4D(func, barfunc, n)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2), INTENT(OUT) :: barfunc
  INTEGER :: n

  barfunc = jfac(:,:,:,:,n)*func

END SUBROUTINE gyro_average_4D

SUBROUTINE gyro_average_5D(func, barfunc)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2), INTENT(OUT) :: barfunc

  barfunc = jfac*func
END SUBROUTINE gyro_average_5D

  SUBROUTINE gyro_average_fields(func,barfunc)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2),INTENT(OUT) :: barfunc
    INTEGER :: n,m

    DO n=ln1,ln2
       DO m=lm1,lm2
          barfunc(:,:,:,m,n) = jfac(:,:,:,m,n)*func
       END DO
    END DO
    
  END SUBROUTINE gyro_average_fields

  FUNCTION get_jfac(i,j,k,m,n)
    INTEGER :: i,j,k,m,n
    real :: get_jfac

    get_jfac = jfac(i,j,k,m,n)
  END FUNCTION get_jfac


  subroutine finalize_gyro_average_ff

    DEALLOCATE(jfac)
    IF (n_fields.gt.2) DEALLOCATE(I1_factor)

  end subroutine finalize_gyro_average_ff

END MODULE gyro_average_ff_mod
