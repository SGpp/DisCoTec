#include "redef.h"
#include "intrinsic_sizes.h"
!>Computes the electromagnetic fields from the modified distribution function
!!\todo clarify naming of routines ??
Module fieldsolver_ff
  Use par_mod
  Use communications
  USE coordinates
  use aux_func
  USE gyro_average_ff_mod
  use vel_space, only: fm, mat_00
  use geometry, ONLY: geom, parscale, major_R
  use charge_curr_dens
  use adiabatic_response
  use hybrid, only: trap
  use equilibrium_fields, only: with_comoving_other, dens_co
  use antenna
  use spatial_averages

  implicit none

  public:: initialize_field_solve_ff, field_solve_ff, field_solve_ff_f, &
       finalize_field_solve_ff, mem_est_field_solve_ff, Apar_pre_antenna,&
       antenna_contrib, a11det_inv_antenna

  private
  
  Real, Dimension(:,:,:), Allocatable:: a11det_inv, a22det_inv, a11det_inv_f
  real, dimension(:,:,:,:), allocatable:: a11det_inv_antenna
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: C1_C1C3C22, C2_C1C3C22, C3_C1C3C22


contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_field_solve_ff(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

 
    !a11det_inv, a11det_inv
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*lijk0
    !a11det_inv_antenna
    IF (antenna_type.ne.0) mem_loc = mem_loc + SIZE_OF_REAL_MB * lijk0
    IF (n_fields .GT. 2) THEN
      !CxC1C3C22
      mem_loc=mem_loc+3*SIZE_OF_REAL_MB*lijk0
    ELSE
      !a22det_inv
      mem_loc=mem_loc+SIZE_OF_REAL_MB*lijk0
    END IF

    mem_est_field_solve_ff=mem_req_in+mem_loc+mem_est_antenna(mem_req_in)
  End Function mem_est_field_solve_ff

  !>Computes the electromagnetic fields from the modified distribution function
  Subroutine field_solve_ff(moments,p_emfields)
    ! Arguments
    !>modified distribution function of type g_1
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields), Intent(IN) :: moments
    !>array of the electromagnetic fields
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(out) :: p_emfields
    complex, dimension(li1:li2):: phi_fs_avg
    INTEGER :: k
    
    DEBUG(1,"start fieldsolve")    
    PERFON('FldSolvesf')

    IF (n_fields .GT. 2) THEN
      ! solve coupled Phi-Bpar system (Apar remains decoupled)
      p_emfields(li1:li2,:,lk1:lk2,1) = moments(:,:,:,1) * C3_C1C3C22 - moments(:,:,:,3) * C2_C1C3C22
      p_emfields(li1:li2,:,lk1:lk2,3) = moments(:,:,:,3) * C1_C1C3C22 - moments(:,:,:,1) * C2_C1C3C22
    ELSE
      ! only independent field equations to solve
      p_emfields(li1:li2,:,lk1:lk2,1) = moments(:,:,:,1) * a22det_inv
    END IF

    if (n_fields.gt.1) p_emfields(li1:li2,:,lk1:lk2,2) =  moments(:,:,:,2) * a11det_inv

    if (antenna_type.ne.0.and.antenna_contrib) &
         call add_Apar_antenna(p_emfields,a11det_inv_antenna)

    if (del_phi) p_emfields(li1:li2,:,lk1:lk2,1) = 0.0
  
    If (delzonal) Then
       Call del_zonal_phi(p_emfields(:,:,:,1))
    End If
    If ((delzonal_fields) .and. (n_fields.ge.2)) Then
       Call del_zonal_phi(p_emfields(:,:,:,2))
    End If
    if (only_Er) then
       if (p_has_0_mode) then
          call flux_surface_average(p_emfields(:,:,:,1),.true.,phi_fs_avg)
          if (yx_order) then
             do k=lbz,ubz
                p_emfields(lj1,:,k,1)=phi_fs_avg(:)
             enddo
             if (nky0.gt.1) p_emfields(lj1+1:lj2,:,:,1)=0.0
          else
             do k=lbz,ubz
                p_emfields(:,lj1,k,1)=phi_fs_avg(:)
             enddo
             if (nky0.gt.1) p_emfields(:,lj1+1:lj2,:,1)=0.0             
          endif
       else
          p_emfields(:,:,:,1) = 0.0
       endif
    endif
    If (abs(add_zonal_phi).gt.epsilon(0.0)) Then
       if (p_has_0_mode) then
          if (yx_order) then
             p_emfields(:,lj1,:,1) = p_emfields(:,lj1,:,1) + add_zonal_phi
          else
             p_emfields(lj1,:,:,1) = p_emfields(lj1,:,:,1) + add_zonal_phi
          endif
       endif
    Endif

    !> set field of 00 mode to zero, for the neoclassical limit
    if (p_has_00_mode) then
      p_emfields(li1,lj1,:,:)=0.0
    endif

    !PERFOFF

    PERFOFF

    DEBUG(1,"==== END   field_solve ====")
  End Subroutine field_solve_ff

  !>Computes the electromagnetic fields from the distribution function f (given without boundaries!)
  Subroutine field_solve_ff_f(p_f_1,p_emfields)
    ! Arguments
    !>modified distribution function of type g_1
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), Intent(IN) :: p_f_1
    !>array of the electromagnetic fields
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(out) :: p_emfields
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields):: moments
    
    DEBUG(1,"start fieldsolve")    
    PERFON('FldSolvesf')

    call calc_charge_curr_dens(p_f_1,moments)
        
    !PERFON('fieldeq')
    
    IF (n_fields .GT. 2) THEN
      ! solve coupled Phi-Bpar system (Apar remains decoupled)
      p_emfields(li1:li2,:,lk1:lk2,1) = moments(:,:,:,1) * C3_C1C3C22 - moments(:,:,:,3) * C2_C1C3C22
      p_emfields(li1:li2,:,lk1:lk2,3) = moments(:,:,:,3) * C1_C1C3C22 - moments(:,:,:,1) * C2_C1C3C22
    ELSE
      ! only independent field equations to solve
      p_emfields(li1:li2,:,lk1:lk2,1) = moments(:,:,:,1) * a22det_inv
    END IF
    
    if (n_fields.gt.1) p_emfields(li1:li2,:,lk1:lk2,2) =  moments(:,:,:,2) * a11det_inv_f

    If (delzonal) Then
       Call del_zonal_phi(p_emfields(:,:,:,1))
    End If
    If ((delzonal_fields) .and. (n_fields.ge.2)) Then
       Call del_zonal_phi(p_emfields(:,:,:,2))
    End If
    !PERFOFF

    PERFOFF

    DEBUG(1,"==== END   field_solve ====")
  End Subroutine field_solve_ff_f

  
  !>Delete the zonal component of phi if delzonal is switched on
  !!Delete the zonal component of A_par if delzonal_fields is switched on
  Subroutine del_zonal_phi(p_phi)
    
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz),intent(inout):: p_phi
    Complex, Dimension(li1:li2) :: zon_phi
    Integer :: i,k
    if (p_has_0_mode) then
       
       PERFON('Zonaldel')
#ifdef WITHOMP_ZONAL
       !$omp parallel do default(none)&
       !$omp shared(pi1, pj1, li1,li2,lj1,lk1,lk2)&
       !$omp shared(zon_phi, p_phi,nz0,geom)
#endif
       Do i=li1,li2
          zon_phi(i) = Sum(p_phi(i,lj1,lk1:lk2)*geom%jacobian(pi1,pj1,lk1:lk2)) / (nz0*geom%avg_jaco)
       Enddo
#ifdef WITHOMP_ZONAL
       !$omp end parallel do
#endif       
       Call sum_to_all_complex(zon_phi, Size(zon_phi), mpi_comm_z)
       
#ifdef WITHOMP_ZONAL
       !$omp parallel do default(none)&
       !$omp shared(lj1,lk1,lk2)&
       !$omp shared(zon_phi, p_phi,delzonal_factor)
#endif
       Do k=lk1,lk2
          p_phi(:,lj1,k) = p_phi(:,lj1,k) - delzonal_factor*zon_phi
       End Do
#ifdef WITHOMP_ZONAL
       !$omp end parallel do
#endif
       PERFOFF
       
    endif
  End Subroutine del_zonal_phi

  !>Initializes the field solver
  Subroutine initialize_field_solve_ff
    real :: k_perp2
    real :: sum0, sum1, sum2, sum3, sum4, sum5, C1C3C22_inv
    real :: zint_mu, zarg
    real, dimension(0:n_spec-1) :: b2spec
    real, dimension(li1:li2, lj1:lj2, 0:n_spec-1):: VIntCorr
    integer::   i, j, k, l, m, n, pni, ierr
#ifdef EXTERNAL_ERF
  EXTERNAL erf
#endif
    real :: erf

    Allocate(a11det_inv(li1:li2, lj1:lj2, lk1:lk2))
    Allocate(a11det_inv_f(li1:li2, lj1:lj2, lk1:lk2))
    IF (antenna_type.ne.0) THEN
       call initialize_antenna
       ALLOCATE(a11det_inv_antenna(li1:li2,lj1:lj2,lk1:lk2,1:2))
    END IF
    IF (n_fields .GT. 2) THEN
       ALLOCATE(C1_C1C3C22(li1:li2,lj1:lj2,lk1:lk2))
       ALLOCATE(C2_C1C3C22(li1:li2,lj1:lj2,lk1:lk2))
       ALLOCATE(C3_C1C3C22(li1:li2,lj1:lj2,lk1:lk2))
    ELSE
       ALLOCATE(a22det_inv(li1:li2,lj1:lj2,lk1:lk2))
    END IF

    pni=pn1
    Do k=lk1,lk2
       ! Loop over all perpendicular wave numbers.
       Do i = li1,li2
          Do j = lj1,lj2
             Do n=ln1,ln2
                if (pn0.gt.1) pni=n
                sum0=0.0
                Do m = lm1,lm2
                   Do l = ll1,ll2
                      ! Calculating the correction factor for amperes law
                      sum0 = sum0 + mat_00(pi1,pj1,k,l,m)*vp(l)**2*fm(pi1,pj1,k,l,m,pni) &
                           & *get_jfac(i,j,k,m,n)*get_jfac(i,j,k,m,n)
                   End Do
                End Do
                ! sum0 has the local sum over velocity space, now calculate the global sum
                sum0=sum_vw(sum0)
                VIntCorr(i,j,n) = sum0
             End Do
          End Do
       End Do

       ! broadcast result to all other species
       Do n=0, n_procs_s-1
          Call mpi_bcast(&
               VIntCorr(li1,lj1,n*ln0), li0*lj0*ln0,&
               MPI_REAL_TYPE, n, mpi_comm_spec, ierr)
       End Do

       Do i = li1,li2
          Do j = lj1,lj2
             k_perp2 = geom%gii(pi1,pj1,k)*ki(i)**2+2.0*geom%gij(pi1,pj1,k)*&
                  &ki(i)*kj(j)+geom%gjj(pi1,pj1,k)*kj(j)**2
             ! remove Bfield for benchmark with Dimits
             do n=0,n_spec-1
                b2spec(n)= k_perp2*spec(n)%mass*spec(n)%temp/&
                     &(spec(n)%charge*geom%Bfield(pi1,pj1,k))**2
             enddo
             
             if (.not.(adiabatic_electrons.or.trap_pass)) then
                ! here all cases with full electron dynamics are covered, 
                ! including adiabatic ions if n_spec=1
                
                sum1 = 0.0
                sum2 = debye2 * k_perp2
                sum4 = 0.0
                sum5 = 0.0
                do n=0,n_spec-1
                   if (.not.spec(n)%passive) then
                      sum1 = sum1 + beta/2*spec(n)%charge**2*spec(n)%dens/spec(n)%mass*VIntCorr(i,j,n)
                      sum2 = sum2 + spec(n)%charge**2*spec(n)%dens/spec(n)%temp * &
                           & (1.0-gamma0(b2spec(n)))*dens_co(pi1,k,n)

                      IF (n_fields .GT. 2) THEN
                         sum4 = sum4 + 2.0 * spec(n)%temp * spec(n)%dens / geom%Bfield(pi1,pj1,k)**2 * &
                              (gamma0(b2spec(n)) - gamma1(b2spec(n)))*dens_co(pi1,k,n)
                         sum5 = sum5 - spec(n)%charge * spec(n)%dens / geom%Bfield(pi1,pj1,k) * &
                              (gamma0(b2spec(n)) - gamma1(b2spec(n)))*dens_co(pi1,k,n)
                      END IF
                   endif
                enddo
                
                IF (n_spec.EQ.1) THEN
                   ! in the one species/adiabatic ions case, the tau parameter represents the
                   ! charges, densities and temperatures of the adiabatic ion species,
                   ! where tau=Zeff*Te/Ti
                   sum2 = sum2 + tau  !ADD DENS_CO????  Adiabatic species and rotation possible???
                END IF
                
                
                If ((ki(i).Ne.0.).Or.(kj(j).Ne.0.)) Then
                   IF (n_fields .GT. 2) THEN
                      C1C3C22_inv = 1.0 / (sum2 * (- 2.0 / beta - sum4) - sum5**2)
                      
                      C1_C1C3C22(i,j,k) = sum2 * C1C3C22_inv
                      C2_C1C3C22(i,j,k) = sum5 * C1C3C22_inv
                      C3_C1C3C22(i,j,k) = (- 2.0 / beta - sum4) * C1C3C22_inv
                   ELSE
                      a22det_inv(i,j,k) = 1./sum2
                   END IF
                   a11det_inv_f(i,j,k) =  1./(k_perp2)
                   a11det_inv(i,j,k)=1./(k_perp2 + 2.0*sum1)
                   
                   IF (antenna_type.ne.0) then
                      a11det_inv_antenna(i,j,k,1) = k_perp2 / (k_perp2 + 2.0 * sum1)
                      a11det_inv_antenna(i,j,k,2)= -2.*sum1/(k_perp2+2.*sum1)
                   endif
                ELSE
                   IF (n_fields .GT. 2) THEN
                      C1_C1C3C22(i,j,k) = 0.0
                      C2_C1C3C22(i,j,k) = 0.0
                      C3_C1C3C22(i,j,k) = 0.0
                   ELSE
                      a22det_inv(i,j,k) = 0.0
                   END IF
                   a11det_inv(i,j,k) = 0.0
                   a11det_inv_f(i,j,k) = 0.0
                   IF (antenna_type.ne.0) a11det_inv_antenna(i,j,k,1:2) = 0.0
                End If
             elseif (adiabatic_electrons) then
                if (with_comoving_other) stop 'adiabatic electrons not implemented in with_comoving_other frame'
                ! ============================================================
                ! Now adiabatic electrons.
                ! ============================================================
                If ((ki(i).Ne.0).Or.(kj(j).Ne.0)) Then
                   ! In the adiabatic electron case, the mass ratio is taken as big, so that
                   ! one can approximate the gamma0-function by 0
                   sum3 = 0.0
                   DO n=0,n_spec-1
                      if (.not.spec(n)%passive) then                      
                         sum3 = sum3 + spec(n)%charge**2 * spec(n)%dens / spec(n)%temp* &
                              & (1.0 - gamma0(b2spec(n)))
                      endif
                   END DO
                   ! === adiabatic electrons ===
                   a22det_inv(i,j,k) = 1.0 / ( sum3  + 1.0 )
                else
                   a22det_inv(i,j,k) = 0.0
                end If
             elseif (no_electron_response) then
                if (with_comoving_other) stop 'no electron response model not implemented in with_comoving_other frame'
                if ((ki(i).ne.0).or.(kj(j).ne.0)) then
                   sum3 = 0.0
                   DO n=0,n_spec-1
                      if (.not.spec(n)%passive) then                      
                         sum3 = sum3 + spec(n)%charge**2 * spec(n)%dens / spec(n)%temp* &
                              & (1.0 - gamma0(b2spec(n)))
                      endif
                   END DO
                   ! === no electron response ===
                   a22det_inv(i,j,k) = 1.0 / sum3 
                else
                   a22det_inv(i,j,k) = 0.0
                end If
             else 
                if (with_comoving_other) stop 'hybrid model not implemented in with_comoving_other frame'
                ! ===============================================================
                ! Hybrid model : adiabatic passing and trapped kinetic electrons 
                ! ===============================================================
                If ((ki(i).Ne.0).Or.(kj(j).Ne.0)) Then
                   sum2 = debye2 * k_perp2
                   do n=ln1,ln2
                      if (.not.spec(n)%passive) then
                         if(spec(n)%charge.ne.-1) then
                            ! contribution from ion species
                            sum2 = sum2 + spec(n)%charge**2 * spec(n)%dens/spec(n)%temp * &
                                 (1.0-gamma0(b2spec(n)))
                         else 
                            ! contribution from electrons :
                            ! Compute integration over mu for polarization drift term ...
                            ! ... on local mu-domain
                            zint_mu = 0.0
                            do m = lm1, lm2
                               zarg = sqrt(mu(m)*geom%Bfield(pi1,pj1,k)) * &
                                    trap(pi1,k) / sqrt(1.0 - trap(pi1,k)**2)
                               
                               zint_mu = zint_mu + mu_weight(m) * &
                                    get_jfac(i,j,k,m,n)**2 * erf(zarg) * &
                                    EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)))
                               
                            end do
                            ! ... now the global sum over mu-domains
                            CALL my_real_sum_w_0d(zint_mu)
                            
                            zint_mu = geom%Bfield(pi1,pj1,k) * zint_mu
                            
                            ! 
                            sum2 = sum2 + spec(n)%charge**2 * spec(n)%dens/spec(n)%temp * &
                                 ( (trap(pi1,k) - zint_mu) + (1.0 - trap(pi1,k)) )
                         end if
                      endif
                   enddo
                   call my_sum_to_all_real_0d(sum2,MPI_COMM_SPEC)
                   a22det_inv(i,j,k) = 1.0/sum2
                   
                Else ! when k_perp = 0, lambda_j = 0 and b_j = 0
                   
                   do n=0,n_spec-1
                      if (.not.spec(n)%passive) then
                         if(spec(n)%charge.eq.-1) then
                            a22det_inv(i,j,k) = 1.0/((1.0 - trap(pi1,k)) * &
                                 &spec(n)%charge**2 * spec(n)%dens/spec(n)%temp)
                         end if
                      end if
                   end do
                End If
             End if ! trap_pass
             
          End Do ! j 
       End Do ! i
    End Do ! k
    
  End Subroutine initialize_field_solve_ff
  
  
  !>Cleanup of the field solver module
  Subroutine finalize_field_solve_ff

    Deallocate(a11det_inv, a11det_inv_f)
    if (antenna_type.ne.0) then
       call finalize_antenna
       DEALLOCATE(a11det_inv_antenna)
    endif
    IF (n_fields .GT. 2) THEN
      DEALLOCATE(C1_C1C3C22,C2_C1C3C22,C3_C1C3C22)
    ELSE
      DEALLOCATE(a22det_inv)
    END IF
    
  End Subroutine finalize_field_solve_ff




end Module fieldsolver_ff
