!>Provides the velocity grid and some velocity space related functions
!!
!!
#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
!define to include profile matrices in the gyroaverage [like (B0F0/T0(i)) in the field equations]
#define GY_WITH_MAT
module vel_space
  use par_in,only: arakawa_zv, perf_vec, trap_pass, coll_f_fm_on, &
       &n0_global, Omega0_tor
  use par_other,only: n_fields
  use discretization
  use boundary_exchange_z, only: exchange_z_3d_equil
  use coordinates
  use gyro_average_df_mod
  use gyro_average_ff_mod
  use lagrange_interpolation
  use geometry
  use axpy
  use hybrid, only: vp_weight_trap, trap, vp_star, Bfieldmax,&
       &initialize_hyb_vars, finalize_hyb_vars
  use equilibrium_fields, only: with_comoving_other, phi0, dens_co

  Implicit none

  PUBLIC :: mat_00, mat_00_trap, mat_10, mat_001, fm
  PUBLIC :: initialize_vel_space, finalize_vel_space, mem_est_vel_space
  PUBLIC :: calc_moments, calc_vsp_moment
  PUBLIC :: initialize_calc_moments_2, calc_moments_2, finalize_calc_moments_2
  public :: mem_est_calc_moments_2
  public :: gyro_op_wrapper

  PRIVATE

  Real, Dimension(:,:,:,:,:,:), Allocatable:: fm
  !background rotation velocity
  Real, Dimension(:,:,:,:,:), Allocatable:: tmp_arr1, tmp_arr2


  ! velocity space integration elements
  Real, Dimension(:,:,:,:,:), Allocatable:: mat_00, mat_00_trap, mat_001
  Real, Dimension(:,:,:,:,:,:), Allocatable:: mat_10
  
#ifdef G_DAGGER
  CHARACTER :: G_op = 'C'
#else
  CHARACTER :: G_op = 'N'
#endif
  

contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_vel_space(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    logical::fm_has_vzbounds
    fm_has_vzbounds = (arakawa_zv.or.(collision_op.ne.'none'.and.coll_f_fm_on))

    !mat_00
    mem_loc=SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0
    !mat_00_trap
    if(trap_pass) mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0
    !mat_10
    mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lklmn0
    !mat_001
    IF (n_fields .GT. 2) mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0
    !fm
    if (fm_has_vzbounds) then
       mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lz0*lv0*lm0*pn0
    else
       mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*pn0
    endif

    !check for use of calc_moments_2
    if (perf_vec(9).eq.2) then
       !calc_mom_temp
!       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lijk0*llm0
    endif
    !!tmpvar,tmpvar2, 
    !mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*(li0+li0)*lj0*lk0

    mem_est_vel_space=mem_req_in+mem_loc

  End Function mem_est_vel_space

  !!\todo The calculation of Blocal (including exchange) should be united with
  !!the identical one in prefactors.F90 -- possibly Bfield could be redefined 
  !!with boundaries.
  subroutine initialize_vel_space
    Implicit None
    integer:: i, j, k, l, m, n, pni, k1, k2, l1, l2
    Real, Dimension(pi1:pi2) :: vabs2
    real, dimension(:,:,:), allocatable :: Blocal
    real, dimension(:,:,:,:), allocatable ::fac_co
    logical::fm_has_vzbounds

    fm_has_vzbounds = (arakawa_zv.or.(collision_op.ne.'none'.and.coll_f_fm_on))

    Allocate(mat_00(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2))
    if(trap_pass) allocate(mat_00_trap(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2)) 
    Allocate(mat_10(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    IF (n_fields .GT. 2) ALLOCATE(mat_001(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))
    if (fm_has_vzbounds) then
       k1=lbz
       k2=ubz
       l1=lbv
       l2=ubv
    else
       k1=lk1
       k2=lk2
       l1=ll1
       l2=ll2
    endif
       
    Allocate(Blocal(pi1:pi2, pj1:pj2, k1:k2))
    Allocate(fm(pi1:pi2, pj1:pj2, k1:k2, l1:l2, lm1:lm2, pn1:pn2))
    if (with_comoving_other) Allocate(fac_co(pi1:pi2,pj1:pj2,k1:k2, ln1:ln2))

#ifndef GY_WITH_MAT
    !if(mype.eq.0)print*,'velocity space moments: gyroaverage does not include mats'
#else
    !if(mype.eq.0)print*,'velocity space moments: gyroaverage includes mats'
#endif

    call set_vp_coordinate_vars
    call set_mu_coordinate_vars
    
    if(trap_pass) call initialize_hyb_vars

    !... init v space integration factors

    ! higher order: fourth order with open formula on left
    ! side and closed formula on right side
    If (nv0.LE.4) Then
      write(*,"(a)") "Error: nv0 must be greater than 4"
      stop
    End If



    ! the perpendicular velocity space integral element is 2*pi v_perp dv_perp = pi \hat{B} dmu.
    ! 
    ! mat_nm is v_parallel^n * mu^m * pi \hat{B} dv_parallel dmu
    ! Note: some mats are renormalized

    Do m=lm1,lm2
       Do l=ll1,ll2
          Do k=lk1,lk2
             Do j=pj1,pj2
                if(trap_pass) then
                   mat_00_trap(pi1:pi2,j,k,l,m) = pi * mu_weight(m)*&
                        &geom%Bfield(pi1:pi2,j,k)* vp_weight_trap(pi1:pi2,k,l,m)
                end if
                mat_00(pi1:pi2,j,k,l,m) = pi * mu_weight(m)*&
                     &geom%Bfield(pi1:pi2,j,k)* vp_weight(l)
             End Do
          End Do
       End Do
    End Do

    do n=ln1,ln2
       Do k=lk1,lk2
          ! higher moments of velocity space
          Do l=ll1,ll2
             Do m=lm1,lm2
                ! parallel velocity
                mat_10(:,:,k,l,m,n) = sqrt(2.*spec(n)%temp/spec(n)%mass) * &
                     vp(l) * mat_00(:,:,k,l,m)
             End Do
          End Do
       End Do
    End Do
    
    IF (n_fields .GT. 2) THEN
      DO k = lk1, lk2
        DO l = ll1, ll2
          DO m = lm1, lm2
            mat_001(:,:,k,l,m) = mu(m) * mat_00(:,:,k,l,m)
          END DO
        END DO
      END DO
    END IF

    !compute bfield with boundaries, if necessary
    Blocal(pi1:pi2,pj1:pj2,lk1:lk2)=geom%Bfield(pi1:pi2,pj1:pj2,lk1:lk2)
    if (fm_has_vzbounds) call exchange_z_3d_equil(Blocal,n0_global*q0)
    if (with_comoving_other) then
       do n=ln1,ln2
          do k=lk1,lk2
             do j=pj1,pj2
                do i=pi1,pi2
                   fac_co(i,j,k,n)=dens_co(i,k,n)
                enddo
             enddo
          enddo
          if (fm_has_vzbounds) call exchange_z_3d_equil(fac_co(:,:,:,n),n0_global*q0)
       enddo
    endif


    pni = pn1
    do n=ln1,ln2
       if (pn0.gt.1) pni=n
       do m=lm1,lm2
          do l=l1,l2
             if (l.lt.0.or.l.gt.(nv0-1)) then
                fm(:,:,:,l,m,pni)=0
             else
                do k=k1,k2
                   do j=pj1,pj2
                      !Squared velocity (temporary)
                      vabs2(:) = (mu(m)*Blocal(pi1:pi2,j,k)+vp(l)**2)/&
                           & spec(n)%temp_prof(pi1:pi2)
                      !Equilibrium distribution function
                      fm(:,j,k,l,m,pni) = spec(n)%dens_prof(pi1:pi2) *&
                           & (pi*spec(n)%temp_prof(pi1:pi2))**(-1.5) *&
                           & exp(-vabs2(:))
                      if (with_comoving_other) fm(:,j,k,l,m,pni) = &
                           & fm(:,j,k,l,m,pni) * fac_co(:,j,k,pni)
                   enddo
                enddo
             endif
          enddo
       enddo
    enddo

    deallocate(Blocal)
    if (with_comoving_other) deallocate(fac_co)

  End Subroutine initialize_vel_space
  
  !>Deallocates velocity space related arrays
  Subroutine finalize_vel_space
    Implicit none

    Deallocate(fm)
    Deallocate(mat_00, mat_10)
    if(trap_pass) then
       call finalize_hyb_vars
       deallocate(mat_00_trap)
    end if
    IF (n_fields .GT. 2) DEALLOCATE(mat_001)

  End Subroutine finalize_vel_space


  !> Calculates velocity space moments (integrations) of J0 gyroaveraged distribution functions
  !!
  !! \param n_mom number of moments to be calculated
  !! \param withbounds indicates whether p_dist includes ghost cells
  !! \param p_dist distribution function
  !! \param p_mat remaining integrand * variables of integration
  !! \param p_mom resulting moment (on the *local* vel. space process)
  !! \param p_gy_op optional parameter defining the gyroaverage operator to be used
  !!        default: 0 (0=J0, 1=I1)
  !!
  !! \note This routine is not written totally general. At the moment it is only possible to
  !! call it in a local context without boundaries and in a global context with boundaries (but
  !! with nyb=0). 
  SUBROUTINE calc_moments(n_mom,withbounds,p_dist,p_mat,p_mom,p_gy_op)
    Implicit none
    INTEGER, INTENT(in):: n_mom
    LOGICAL, INTENT(in):: withbounds
    COMPLEX, DIMENSION(:,:,:,:,:,:),INTENT(in):: p_dist
    REAL, DIMENSION(pi1:pi2,pj1:pj2,lk1:lk2,1:n_mom,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(in):: p_mat 
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_mom,ln1:ln2),INTENT(out):: p_mom
    INTEGER, DIMENSION(1:n_mom), optional :: p_gy_op

    ! Local variables
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: g_av
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: dist_mat
    INTEGER, DIMENSION(1:n_mom) :: gy_op
    INTEGER:: j,k,l,m,n,i_mom,joff,koff,loff,moff,noff

    PERFON('calc_moments')
    IF (present(p_gy_op)) THEN
       gy_op = p_gy_op
    ELSE
       gy_op = 0
    ENDIF

    p_mom = cmplx(0.0,0.0)

    ! the passed array p_dist is running locally from 1 to the end
    ! dependent on if we have the boundaries included, the calculation of
    ! the inner point indices is different
    ! the first two dimensions are always without boundaries, that means
    ! running over 1:li0 1:lj0
    ! 
    joff = 1-lj1
    if (withbounds) then
       koff = 1-lk1+nzb
       loff = 1-ll1+nvb
       moff = 1-lm1+nwb
    else
       koff = 1-lk1
       loff = 1-ll1
       moff = 1-lm1
    endif
    noff=1-ln1
    
    DO n=ln1,ln2
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2

                IF (xy_local) THEN
                   DO i_mom = 1, n_mom
                      !gyro average
                      IF (i_mom.EQ.1) THEN
                         CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff,n+noff),g_av,k,m,n,gy_op(i_mom))
                      ELSEIF (gy_op(i_mom).NE.gy_op(i_mom-1)) THEN
                         CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff,n+noff),g_av,k,m,n,gy_op(i_mom))
                      ENDIF

                      !calculation of the moments
                      ! p_mom(:,:,k,i_mom,n) = p_mom(:,:,k,i_mom,n) &
                      !     & + p_mat(pi1,lj1,k,i_mom,l,m,n) * g_av
                      CALL axpy_ij(lij0,p_mat(pi1,pj1,k,i_mom,l,m,n),g_av,&
                           & p_mom(:,:,k,i_mom,n))
                   ENDDO
                ELSE
                   ! in the following call, we assume that nyb = 0. If that is not the
                   ! case, one has to rewrite the gyro_average method
                   DO i_mom = 1, n_mom
#ifndef GY_WITH_MAT
                      IF (i_mom.EQ.1) THEN
                         CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff,n+noff),g_av,k,m,n,gy_op(i_mom))
                      ELSEIF (gy_op(i_mom).NE.gy_op(i_mom-1)) THEN
                         CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff,n+noff),g_av,k,m,n,gy_op(i_mom))
                      ENDIF
                      DO j=lj1,lj2
                         p_mom(:,j,k,i_mom,n) = p_mom(:,j,k,i_mom,n) &
                              & + p_mat(:,pj1,k,i_mom,l,m,n)*g_av(:,j)
                      END DO
#else
                      do j=lj1,lj2
                         dist_mat(:,j)= p_dist(:,j+joff,k+koff,l+loff,m+moff,n+noff)&
                              &* p_mat(:,pj1,k,i_mom,l,m,n)
                      enddo
                      CALL gyro_op_wrapper(dist_mat,g_av,k,m,n,gy_op(i_mom))
                      p_mom(:,:,k,i_mom,n) = p_mom(:,:,k,i_mom,n) + g_av
#endif
                   END DO
                END IF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! In p_mom we have now in each process the moments over 
    ! the process local velocity space
    ! Note: Sum over all MPI processes needs to be done

    PERFOFF

  END SUBROUTINE calc_moments

  subroutine initialize_calc_moments_2
!    allocate(calc_mom_temp(lijk0,llm0))
  end subroutine initialize_calc_moments_2

  !>Give an estimate of the memory requirements
  Real Function mem_est_calc_moments_2()
    !calc_mom_temp
    mem_est_calc_moments_2 = 0
!    mem_est_calc_moments_2 = SIZE_OF_COMPLEX_MB*lijk0*llm0
  End Function mem_est_calc_moments_2

  subroutine finalize_calc_moments_2
!    deallocate(calc_mom_temp)
  end subroutine finalize_calc_moments_2

  !> Alternative calculation of velocity space moments (integrations)
  !! Gyroaverage operations need to be done in advance or included in
  !! the p_mat prefactors
  !!
  !! \param momnum number of moments to be calculated
  !! \param sp_dim space dimensions
  !! \param vsp_dim velocity space dimensions
  !! \param p_g_1 distribution function
  !! \param p_mat remaining integrand * variables of integration
  !! \param p_mom resulting moment (on the *local* vel. space process)
  SUBROUTINE calc_moments_2(sp_dim,vsp_dim,momnum,p_g_1,p_mat,p_mom)
    Implicit none
    INTEGER, INTENT(in):: sp_dim,vsp_dim,momnum
    COMPLEX, DIMENSION(sp_dim,vsp_dim,ln1:ln2), INTENT(in):: p_g_1
    REAL, DIMENSION(sp_dim,vsp_dim,1:momnum,ln1:ln2), INTENT(in):: p_mat 
    COMPLEX, DIMENSION(sp_dim,momnum,ln1:ln2), INTENT(out):: p_mom

    ! Local variables
    INTEGER:: vsp, n, mnum    

    !PRINT*,mype,': sp_dim = ',sp_dim,', vsp_dim = ',vsp_dim,' ==> ',sp_dim*vsp_dim*16,' Bytes'
    PERFON('calcmom2')    
    do n=ln1,ln2
       do mnum=1,momnum
          p_mom(:,mnum,n)=p_mat(:,1,mnum,n)*p_g_1(:,1,n)
          do vsp=2,vsp_dim
             p_mom(:,mnum,n)=p_mom(:,mnum,n)+p_mat(:,vsp,mnum,n)*p_g_1(:,vsp,n)
          enddo
!          calc_mom_temp=p_mat(:,:,mnum,n)*p_g_1(:,:,n)
!          p_mom(:,mnum,n)=sum(calc_mom_temp,2)
       enddo
    enddo
    PERFOFF

  END SUBROUTINE calc_moments_2


!!************************************************************************************
!!A full velocity space moment of the particle coordinate distribution function involves some 
!!further (FLR) terms which depend on phi and Bpar (Apar drops out for local Maxwellians)
!!The prefactors of phi and Bpar which include mu-Integrations of the gyroaverage operators
!!can be computed during the initialization

  !>Computes velocity space moments of the particle coordinate distribution function (thus including FLR effects)
  !!
  !!This routine computes velocity space moments
  !!of the particle coordinate distribution function 
  !!which includes velocity space averaged and gyroaveraged
  !!electrostatic potentials and parallel magnetic fields
  !!
  !! \param n_mom number of moments to be calculated
  !! \param p_dist distribution function
  !! \param p_emfields eletrostatic(magnetic) fields
  !! \param p_mat remaining integrand * variables of integration
  !! \param p_mom resulting moment (on the *local* vel. space process)
  !!
  !! \note This routine is not written totally general:
  !!       no boundary points in x and y are allowed
  SUBROUTINE calc_vsp_moment(n_mom,p_dist,withbounds,p_emfields,p_mat,p_mom,nc_corr)
    Implicit none
    INTEGER, INTENT(in):: n_mom
    COMPLEX, DIMENSION(:,:,:,:,:,:), INTENT(in) :: p_dist
    logical, intent(in)::withbounds, nc_corr
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in)   :: p_emfields
    REAL, DIMENSION(pi1:pi2,pj1:pj2,lk1:lk2,1:n_mom,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(in):: p_mat 
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_mom,ln1:ln2),INTENT(out):: p_mom
    Complex, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: tmpvar
    Complex, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: tmpvar2
    Complex, DIMENSION(li1:li2) :: gy_tmpvar, pullback_dist_func
    Real, DIMENSION(pi1:pi2) :: invtemp

    Integer :: j,k,l,m,n,o,pni
    Integer :: ioff,joff,koff,loff,moff,noff

    PERFON('calc_vsp_mom')
    
    ioff = 1-li1
    joff = 1-lj1
    if (withbounds) then
       koff = 1-lk1+nzb
       loff = 1-ll1+nvb
       moff = 1-lm1+nwb
    else
       koff = 1-lk1
       loff = 1-ll1
       moff = 1-lm1
    endif
    noff = 1-ln1

    
    p_mom = cmplx(0.0,0.0)

    pni=pn1
    Do n=ln1,ln2
       if (pn0.gt.1) pni=n
       invtemp = 1.0/(spec(n)%temp*spec(n)%temp_prof(pi1:pi2))
       DO m=lm1,lm2
          tmpvar(li1:li2,lj1:lj2,lk1:lk2) = spec(n)%charge*p_emfields(li1:li2,lj1:lj2,lk1:lk2,1)
          if (xy_local) then
             tmpvar2(:,:,:)=jfac(:,:,lk1:lk2,m,n)*tmpvar(:,:,:)*invtemp(pi1)
             IF (n_fields .GT. 2) tmpvar2 = tmpvar2 + spec(n)%temp * mu(m) * &
                p_emfields(li1:li2,lj1:lj2,lk1:lk2,3) * &
                I1_factor(li1:li2,lj1:lj2,lk1:lk2,m,n) * invtemp(pi1)
          else
             call gyro_average_df(tmpvar,tmpvar2,m,n,transA_in='N') !fields don't use "daggered" operators

             Do k=lk1,lk2
                Do j=lj1,lj2
                   tmpvar2(li1:li2,j,k)=tmpvar2(li1:li2,j,k)*invtemp(li1:li2)
                End Do
             End Do
          endif
          Do l=ll1,ll2
             Do k=lk1,lk2
                Do j=lj1,lj2
                   !!Gyroaverage of F1+(q <phi> + T(x0) mu <Bpar>) F0/T(x)
                   if (xy_local) then

                      tmpvar(li1:li2,j,k) = p_dist(li1+ioff:li2+ioff,j+joff,k+koff,l+loff,m+moff,n+noff)&
                           &+tmpvar2(li1:li2,j,k)*fm(pi1,pj1,k,l,m,pni)
                      call gyro_average_ff(tmpvar(:,j,k),gy_tmpvar,j,k,m,n)

                      !! T* F = <F1 + (q <phi> + T(x0) mu <Bpar>) F0/T(x)> - q phi F0/T0
                      pullback_dist_func = gy_tmpvar - &
                        spec(n)%charge * p_emfields(li1:li2,j,k,1) * &
                        fm(pi1,pj1,k,l,m,pni) * invtemp(pi1)

                      !Now add p_mat which contains e.g. pi Bfield dvpar dmu for the 00-moment
                      do o=1,n_mom
                         p_mom(li1:li2,j,k,o,n) = p_mom(li1:li2,j,k,o,n)+pullback_dist_func(li1:li2)*p_mat(pi1,pj1,k,o,l,m,n)
                      enddo

                   else
#ifndef GY_WITH_MAT
                      tmpvar(li1:li2,j,k) = p_dist(li1+ioff:li2+ioff,j+joff,k+koff,l+loff,m+moff,n+noff)&
                           &+tmpvar2(li1:li2,j,k)*fm(li1:li2,pj1,k,l,m,n)
                      call gyro_average_df(tmpvar(:,j,k),gy_tmpvar,j,k,m,n,transA_in=G_op)
                      !! particle coordinate distribution function
                      !! T* F = <F1 + [q <phi> + T(x0) mu <Bpar>] F0/T0(x)> - q phi F0/T0(x)
                      pullback_dist_func(li1:li2) = gy_tmpvar(li1:li2)-spec(n)%charge*p_emfields(li1:li2,j,k,1)*&
                           &fm(li1:li2,pj1,k,l,m,n)*invtemp(li1:li2)

                      do o=1,n_mom
                         p_mom(li1:li2,j,k,o,n) = p_mom(li1:li2,j,k,o,n)+pullback_dist_func(li1:li2)*&
                              &p_mat(li1:li2,pj1,k,o,l,m,n)
                      enddo
#else
                      do o=1,n_mom
                         !include mat in gyroaverage (like in the field equations)
                         tmpvar(li1:li2,j,k) = (p_dist(li1+ioff:li2+ioff,j+joff,k+koff,l+loff,m+moff,n+noff)&
                              &+tmpvar2(li1:li2,j,k)*fm(li1:li2,pj1,k,l,m,n))*p_mat(li1:li2,pj1,k,o,l,m,n)
                         call gyro_average_df(tmpvar(:,j,k),gy_tmpvar,j,k,m,n,transA_in=G_op)
                         !! particle coordinate distribution function
                         !! T* F = <F1 + [q <phi> + T(x0) mu <Bpar>] F0/T0(x)> - q phi F0/T0(x)
                         !add phi term without gyroaverage of F0/T0*mat (the same as in field equations)
                         pullback_dist_func(li1:li2) = gy_tmpvar(li1:li2)-spec(n)%charge*p_emfields(li1:li2,j,k,1)*&
                              &fm(li1:li2,pj1,k,l,m,n)*invtemp(li1:li2)*p_mat(li1:li2,pj1,k,o,l,m,n)
                         if (nc_corr) then
                            pullback_dist_func(li1:li2) = pullback_dist_func(li1:li2) + tmpvar2(li1:li2,j,k)*&
                              & fm(li1:li2,pj1,k,l,m,n)*p_mat(li1:li2,pj1,k,o,l,m,n)
                         endif
                         p_mom(li1:li2,j,k,o,n) = p_mom(li1:li2,j,k,o,n)+pullback_dist_func(li1:li2)
                      enddo
#endif
                   endif
                End Do
             End Do
          End Do
       End Do
    End Do

    PERFOFF

  END SUBROUTINE calc_vsp_moment

  !!computes the gyroaverage of a distribution function type variable
  !! \param p_dist distribution function variable to be gyro-averaged
  !! \param g_av result
  !! \param k,m,n indices
  !! \param gy_op switch controlling the use of I1 or J0 -type gyroaverage
  SUBROUTINE gyro_op_wrapper(p_dist,g_av,k,m,n,gy_op)
    Implicit none
    COMPLEX, DIMENSION(:,:),INTENT(in):: p_dist
    COMPLEX, DIMENSION(li1:li2,lj1:lj2),INTENT(OUT) :: g_av
    INTEGER:: gy_op,k,m,n

    IF (xy_local) THEN
       IF (gy_op.EQ.1) THEN 
          !For Bpar applications, we need I1 instead of jfac for gyroaverage.
          CALL gyro_average_ff_bpar(p_dist,g_av,k,m,n)
       ELSE 
          !J0 gyroaverage
          CALL gyro_average_ff(p_dist,g_av,k,m,n) 
       ENDIF
    ELSE
       IF (gy_op.EQ.1) THEN
          stop 'Bpar not implemented in the global code'
       ELSE
          !ifdef G_DAGGER : use !G^\dagger gyroaverage operators
          !(G_op='C') to obtain gyroaverage in particle space
          !(not to be used for emfields)
          CALL gyro_average_df(p_dist,g_av,k,m,n,transA_in=G_op) 
       ENDIF
    ENDIF

  END SUBROUTINE gyro_op_wrapper
  
  
  
end module vel_space

