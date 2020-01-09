#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
!> This modules contains routines for post-processing of 
!! energy terms
Module diagnostics_energy
  Use par_mod 
  Use file_io, only: get_unit_nr
  Use vel_space, only: fm, mat_00
  use calc_rhs, only: calc_rhs_only,ptr_bar_emfields, ptr_barchi, &
       ptr_dbarchidxy, ptr_dgdxy, rhs_nl, chi_block, rhs_f0,this_nonlinear_term
  use geometry
  use communications
  use aux_fields, only: calc_aux_fields
  use fieldsolver_ff, only: a11det_inv_antenna
  Use gyro_average_ff_mod, only: gyro_average_ff,gyro_average_ff_bpar, get_jfac
  Use gyro_average_df_mod, only: gyro_average_df
  use collisions, only: equ_collisions
  use dfdzv_terms !, only: add_dfdz, add_dfdv, add_hypv, add_hypz
  use dgdxy_terms
  use dzv_terms
  use dfdxy_terms, only: add_dfdxy
  use dchidz_term, only: add_dchidz
  use dchidxy_terms, only: add_dchidxy, add_dchidxy_orig
  use sources_mod, only: add_f1_sources_block, explicit_buffer, add_kBuffer_explicit
  !USE nonlinearity, ONLY: add_nonlinearity
  use blockindex
  use axpy
  use prefactors
  use numerical_damping
  use boundaries, only: pb_xshift, x_boundary_block, exchange_x
  use f0_term, only: add_f0_term
  use antenna, only: antenna_type, n_antennae, get_Apar_antenna, get_dApar_dt_antenna, Apar_pre_antenna, &
       add_dApar_dt_antenna, antenna_contrib
  use spatial_averages
#ifdef WITHFUTILS
  use futils
  use hashtable
#endif
  Implicit None

  PUBLIC :: initialize_all_diags_energy, exec_all_diags_energy, &
       &finalize_all_diags_energy, mem_est_diag_energy,get_energy_total,&
       &check_diag_energy, get_cfgamma_h, get_cfgamma_f, get_cfgamma_g, get_cfgamma_ant,&
       &get_energy_terms, energy_integral, scalar_product, block_mult, diag_energy_exchange
       

  PRIVATE 

  INTEGER:: ENERGYFILE,ETERMSFILE,ENEXCHANGEFILE
  real :: last_energy
  real,dimension(:,:,:),allocatable:: kperp2
#ifdef WITHFUTILS
  integer :: fidenergy_h5, fidenergy3d_h5
  integer :: isnap_energy3d = 0, isnap_energy = 0
  integer, parameter :: bufsize = 20
  TYPE(BUFFER_TYPE),SAVE :: hbufenergy_h5
  ! IBM compilers demand:
  ! A variable declared in the scope of a module, hbufenergy_h5, 
  ! that is of a derived type with default initialization, must have the SAVE attribute.
  integer, dimension(2) :: dims
  integer :: rank
  character(len=FILENAME_MAX) :: dset_name_energy3d
#endif

Contains

!!!******************************************************************!!!

  SUBROUTINE check_diag_energy
    if (turbdeal) then
       stop 'The energy diagnostic is not implemented for turbo dealiasing yet'
       !this would require some extra treatment of the nonlinearity
       !where several stage need to be called
    endif
    if(istep_energy_exchange.gt.0) then
       if(.not.xy_local.and.mype==0) stop "istep_energy_exchange not ready for global." 
       if(n_fields.gt.2) then
          stop "istep_energy_exchange not ready for bpar (need to verify that bpar is in chi)."
       end if
    end if
  END SUBROUTINE check_diag_energy

  !> Estimates the memory required by energy diags
  Real Function mem_est_diag_energy(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc = 0.0

    !in sum_int_3d: v3d_tmp
    mem_loc = mem_loc + SIZE_OF_REAL_MB*lijk0

    !in get_cfgamma: phi_bar
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0

    !in get_energy_terms:
    !rhs,cfgamma,g_temp (temporary!)
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*3*lijklmn0
    !f_in, h_in (temporary!)
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*lzvw0*ln0
    if (arakawa_zv) mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*lzvw0*ln0

    !fields_in,dE_terms
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*(n_fields+11)

    !in energy_write_3d:
    !ar3dr
    mem_loc = mem_loc + SIZE_OF_REAL_MB*lijk0
    !energyx,energyf,dest3d
    mem_loc = mem_loc + SIZE_OF_REAL_MB*nx0*(lj0*lk0+nky0*(lk0+nz0))

    !in get_energy_3d:
    !cfgamma
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*1*lijklmn0
    !f_in (temporary!)
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*lzvw0*ln0    
    !fields_in
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*n_fields    

    !in get_energy_total: energy
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lijk0

    mem_est_diag_energy=mem_req_in+mem_loc
  End Function mem_est_diag_energy

!!!******************************************************************!!!

  SUBROUTINE initialize_all_diags_energy

    if (istep_energy.gt.0) call initialize_diag_rhs_energy
    if (istep_energy3d.gt.0) call initialize_diag_energy
    if (istep_energy_exchange.gt.0) call initialize_diag_energy_exchange
    last_energy=0.0

  END SUBROUTINE initialize_all_diags_energy

!!!******************************************************************!!!

  SUBROUTINE exec_all_diags_energy(itime,istep_energy,istep_energy3d)
    integer, intent(in) :: itime, istep_energy, istep_energy3d
    logical :: write1d,write3d

    write1d = .false.
    write3d = .false.

    IF (do_energy) THEN
       IF (istep_energy.GT.0) write1d = (MODULO(itime,istep_energy).eq.0)
       IF (istep_energy3d.GT.0) write3d = (MODULO(itime,istep_energy3d).eq.0)
       IF (write1d.or.write3d) CALL diag_energy(write1d,write3d)
    END IF

  END SUBROUTINE exec_all_diags_energy

!!!******************************************************************!!!

  SUBROUTINE finalize_all_diags_energy

    if (istep_energy.gt.0) call finalize_diag_rhs_energy
    if (istep_energy3d.gt.0) call finalize_diag_energy
    if (istep_energy_exchange.gt.0) call finalize_diag_energy_exchange

  END SUBROUTINE finalize_all_diags_energy
!!!******************************************************************!!!



!!!******************************************************************!!!
!!!******************************************************************!!!
  Subroutine initialize_diag_energy

    if(write_std) then
       call get_unit_nr(ENERGYFILE)
       
       if(mype.eq.0) open(ENERGYFILE, file=trim(diagdir)//'/energy3d'//&
            trim(file_extension),form='unformatted',status='unknown')
    end if
    
#ifdef WITHFUTILS
    if(write_h5) then
       if ((my_pev+my_pew+my_pespec).eq.0) then
          call creatf(trim(diagdir)//'/energy3d'//trim(file_extension)//'.h5', &
               fidenergy3d_h5, "energy3d", 'd', mpi_comm_xyz)
          call creatg(fidenergy3d_h5, '/energy3d')
          call creatg(fidenergy3d_h5, '/energy3d/var1')
          call creatg(fidenergy3d_h5, '/energy3d/var2')
          call creatg(fidenergy3d_h5, '/energy3d/var3')
          call creatg(fidenergy3d_h5, '/energy3d/var4')
          
          rank = 0
          call creatd(fidenergy3d_h5, rank, dims, "/energy3d/time", "time")
       end if
    end if
#endif
  End Subroutine Initialize_diag_energy


  SUBROUTINE finalize_diag_energy
    if(write_std) then
       if (mype.eq.0) close(ENERGYFILE)
    end if
#ifdef WITHFUTILS
    if(write_h5) then
       if ((my_pev+my_pew+my_pespec).eq.0) then
          call closef(fidenergy3d_h5)
       end if
    end if
#endif
  END SUBROUTINE finalize_diag_energy

  Subroutine initialize_diag_energy_exchange

    if(write_std) then
       call get_unit_nr(ENEXCHANGEFILE)
       if(mype.eq.0) open(ENEXCHANGEFILE, file=trim(diagdir)//'/energy_exchange'//&
            trim(file_extension),status='unknown')
       write(*,*) "ENEXCHANGEFILE",ENEXCHANGEFILE
    end if
    
!Needs to be modified for hdf5
!#ifdef WITHFUTILS
!    if(write_h5) then
!       if ((my_pev+my_pew+my_pespec).eq.0) then
!          call creatf(trim(diagdir)//'/energy3d'//trim(file_extension)//'.h5', &
!               fidenergy3d_h5, "energy3d", 'd', mpi_comm_xyz)
!          call creatg(fidenergy3d_h5, '/energy3d')
!          call creatg(fidenergy3d_h5, '/energy3d/var1')
!          call creatg(fidenergy3d_h5, '/energy3d/var2')
!          call creatg(fidenergy3d_h5, '/energy3d/var3')
!          call creatg(fidenergy3d_h5, '/energy3d/var4')
!          
!          rank = 0
!          call creatd(fidenergy3d_h5, rank, dims, "/energy3d/time", "time")
!       end if
!    end if
!#endif
  End Subroutine initialize_diag_energy_exchange

  SUBROUTINE finalize_diag_energy_exchange
    if(write_std) then
       if (mype.eq.0) close(ENEXCHANGEFILE)
    end if
!Needs to be modified for hdf5
!#ifdef WITHFUTILS
!    if(write_h5) then
!       if ((my_pev+my_pew+my_pespec).eq.0) then
!          call closef(fidenergy3d_h5)
!       end if
!    end if
!#endif
  END SUBROUTINE finalize_diag_energy_exchange




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_diag_rhs_energy
    implicit none
    integer:: i,j,k

    if(write_std) then
       if (mype==0) then
          call get_unit_nr(ETERMSFILE)    
          open(ETERMSFILE, file=trim(diagdir)//'/energy'//&
               &trim(file_extension),status='unknown')
          write(ETERMSFILE,*) "#  1.  time"
          write(ETERMSFILE,*) "#  2.  Total Energy" 
          write(ETERMSFILE,*) "#  3.  Total dE/dt"
          if (only_neo) then
             write(ETERMSFILE,*) "#  4.  Drive 0 term (neoclassical computation)"
          else
             write(ETERMSFILE,*) "#  4.  Drive 1 term"
          endif
          write(ETERMSFILE,*) "#  5.  Heat sources (nonlocal)"
          write(ETERMSFILE,*) "#  6.  Collisional Dissipation"
          write(ETERMSFILE,*) "#  7.  hyp_v dissipation"
          write(ETERMSFILE,*) "#  8.  hyp_z dissipation"
          write(ETERMSFILE,*) "#  9.  hyp_x/y dissipation"
          write(ETERMSFILE,*) "#  10. nonlinearity"
          write(ETERMSFILE,*) "#  11. [Poisson err]_z,vpar"
          write(ETERMSFILE,*) "#  12. curvature + remaining terms (if any)"
          write(ETERMSFILE,*) "#  13. (abs(10)+abs(11)+abs(12))/abs(4): test for conservation"
          write(ETERMSFILE,*) "#  14. change in energy: (E(t)-E(t_last_energy_calc))/(dt*istep_energy)"
100       format ('#',"        1",14I12)
          write(ETERMSFILE,100) 2,3,4,5,6,7,8,9,10,11,12,13,14
       end if
    end if

#ifdef WITHFUTILS
    if(write_h5) then
       if ((mype).eq.0) then
          call creatf(trim(diagdir)//'/energy'//trim(file_extension)//'.h5', &
               fidenergy_h5, "energy", 'd')
          call creatg(fidenergy_h5, '/energy')
          call htable_init(hbufenergy_h5, bufsize)
          call set_htable_fileid(hbufenergy_h5, fidenergy_h5, '/energy')
       end if
    end if
#endif

    if (any(antenna_type.eq.(/2,3/))) then
       allocate(kperp2(li1:li2,lj1:lj2,lk1:lk2))
       do k=lk1,lk2
          do j=lj1,lj2
             do i=li1,li2
                kperp2(i,j,k)=geom%gii(pi1gl,pj1,k)*ki(i)**2+2.*geom%gij(pi1gl,pj1,k)*ki(i)*kj(j)+&
                     geom%gjj(pi1gl,pj1,k)*kj(j)**2
             enddo
          enddo
       enddo
    endif
  end subroutine initialize_diag_rhs_energy

  Subroutine finalize_diag_rhs_energy
    if(write_std) then
       if (mype.eq.0) then
          close(ETERMSFILE)
       endif
    end if

#ifdef WITHFUTILS
    if(write_h5) then
       if(mype.eq.0) then
          call htable_hdf5_flush(hbufenergy_h5)
          call closef(fidenergy_h5)
       end if
    end if
#endif
    if (any(antenna_type.eq.(/2,3/))) deallocate(kperp2)
  End Subroutine finalize_diag_rhs_energy


  !>This subroutine performs the mat_00 integral and sum over species in order to reduce from 6D to 3D
  subroutine energy_integral(v6d,v3d)
    !This subroutine performs the mat_00 integral and sum over species

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: v6d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(out) :: v3d
!    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d_temp
    integer :: j,k,l,m,n

    v3d=(0.0,0.0)
    if (xy_local) then
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   call axpy_ij(lij0,mat_00(pi1,pj1,k,l,m),v6d(:,:,k,l,m,n),&
                        &v3d(:,:,k))
                enddo
             enddo
          enddo
       enddo
    else
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      v3d(:,j,k) = v3d(:,j,k) + v6d(:,j,k,l,m,n)*mat_00(li1:li2,pj1,k,l,m)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    call my_complex_sum_vwspec(v3d,size(v3d)) 
!    v3d=v3d+v3d_temp

  end subroutine energy_integral

  subroutine energy_integral_vw(v5d,v3d)
    !This subroutine performs the mat_00 integral and sum over species

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2), intent(in) :: v5d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(out) :: v3d
    integer :: j,k,l,m

    v3d=(0.0,0.0)
    if (xy_local) then
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   call axpy_ij(lij0,mat_00(pi1,pj1,k,l,m),v5d(:,:,k,l,m),&
                        &v3d(:,:,k))
                enddo
             enddo
          enddo
    else
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      v3d(:,j,k) = v3d(:,j,k) + v5d(:,j,k,l,m)*mat_00(li1:li2,pj1,k,l,m)
                   enddo
                enddo
             enddo
          enddo
    endif
    call my_complex_sum_vw(v3d(:,:,:),size(v3d(:,:,:))) 

  end subroutine energy_integral_vw


  !>This subroutine performs the mat_00 integral and sum over species in order to reduce from 6D to 3D
  !!(strip mining version)
  subroutine add_energy_integral_block(v6d,v3d,lb1,lb2)
    implicit none
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(in) :: v6d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(inout) :: v3d
    integer :: i,j,klmn

    if (xy_local) then
       do klmn=lb1,lb2
          call axpy_ij(lij0,mat_00(pi1,pj1,sk(klmn),sl(klmn),sm(klmn)),&
               &v6d(:,:,klmn),v3d(:,:,sk(klmn)))
       enddo
    else
       do klmn=lb1,lb2
          do j=lj1,lj2
             do i=li1,li2
                v3d(i,j,sk(klmn)) = v3d(i,j,sk(klmn)) + &
                     &v6d(i,j,klmn)*mat_00(i,pj1,sk(klmn),sl(klmn),sm(klmn))
             enddo
          enddo
       enddo
    endif

!    call my_complex_sum_vwspec(v3d,size(v3d)) 

  end subroutine add_energy_integral_block


  subroutine scalar_product(v1,v2,sp)

    implicit none
    complex, intent(in), dimension(li0*lj0*lk0*ll0*lm0*ln0) :: v1,v2
    complex, intent(out) :: sp

    sp=sum(conjg(v1)*v2)

  end subroutine scalar_product


  !>This subroutine returns Cj/F0*((fj*) + qj phi_bar_j*F0/Tj +mu bpar_bar_j*F0)
  !!Alternatively, one could use Cj/F0*((gj*) + (qj chi_bar_j*F0/Tj))
  SUBROUTINE get_cfgamma_f(f_in,fields_in,cfgamma)
    Implicit None
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: fields_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: f_in    
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: cfgamma
    complex, dimension(li1:li2, lj1:lj2):: bpar_bar
    complex, dimension(li1:li2, lj1:lj2):: psi_term
    integer :: j,k,l,m,n,pni
    
    pni=pn1
    
    if (xy_local) then 
       DO n=ln1,ln2
          if (pn0.gt.1) pni=n
          DO m=lm1,lm2
             DO l=ll1,ll2
                DO k=lk1,lk2
                !psi_term=qj/Tj phi_bar (+ mu bpar_bar)
                CALL gyro_average_ff(fields_in(:,:,k,1),psi_term,k,m,n)
                if (n_fields.gt.2) then
                   !phi and bpar contribution
                   CALL gyro_average_ff_bpar(fields_in(:,:,k,3),bpar_bar,k,m,n)
                   psi_term = psi_term*spec(n)%charge/spec(n)%temp+mu(m)*bpar_bar
                else
                   !only phi contribution
                   psi_term=psi_term*spec(n)%charge/spec(n)%temp
                endif
                cfgamma(:,:,k,l,m,n)=spec(n)%temp*spec(n)%dens/fm(pi1,pj1,k,l,m,pni)*&
                     &(CONJG(f_in(:,:,k,l,m,n))+(fm(pi1,pj1,k,l,m,pni)*CONJG(psi_term)))
                END DO
             END DO
          END DO
       END DO
    else
       DO n=ln1,ln2
          DO m=lm1,lm2
             DO l=ll1,ll2
                DO k=lk1,lk2
                   !here psi_term = phi_bar (!)   for practical reasons
                   call gyro_average_df(fields_in(:,:,k,1),psi_term,k,m,n) !fields don't use "daggered" operators
                   DO j=lj1,lj2
                      cfgamma(:,j,k,l,m,n)=spec(n)%temp*spec(n)%temp_prof(li1:li2)*spec(n)%dens/&
                           & fm(li1:li2,pj1,k,l,m,n)*(CONJG(f_in(:,j,k,l,m,n)))+&
                           & spec(n)%dens*spec(n)%charge*CONJG(psi_term(:,j))
                   END DO
                END DO
             END DO
          END DO
       END DO       
    endif

  END SUBROUTINE get_cfgamma_f

  SUBROUTINE get_cfgamma_g(g_in,cfgamma_g)
    Implicit None
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: g_in    
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: cfgamma_g
    integer :: j,k,l,m,n,pni
    pni=pn1

    if (xy_local) then 
       DO n=ln1,ln2
          if (pn0.gt.1) pni=n
          DO m=lm1,lm2
             DO l=ll1,ll2
                DO k=lk1,lk2
                   cfgamma_g(:,:,k,l,m,n)=spec(n)%temp*spec(n)%dens/fm(pi1,pj1,k,l,m,pni)*&
                        &(CONJG(g_in(:,:,k,l,m,n)))
                END DO
             END DO
          END DO
       END DO
    else
       DO n=ln1,ln2
          DO m=lm1,lm2
             DO l=ll1,ll2
                DO k=lk1,lk2
                   DO j=lj1,lj2
                      cfgamma_g(:,j,k,l,m,n)=spec(n)%temp*spec(n)%temp_prof(li1:li2)*spec(n)%dens/&
                           & fm(li1:li2,pj1,k,l,m,n)*(CONJG(g_in(:,j,k,l,m,n)))
                   END DO
                END DO
             END DO
          END DO
       END DO       
    endif

  END SUBROUTINE get_cfgamma_g


  !>This subroutine returns Cj/F0*(hj*)=cfgamma_h
  subroutine get_cfgamma_h(h_in,cfgamma_h)
    Implicit None
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: h_in    
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: cfgamma_h
    integer :: j,k,l,m,n,pni
    pni=pn1

    if (xy_local) then
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   cfgamma_h(:,:,k,l,m,n)=spec(n)%temp*spec(n)%dens/fm(pi1,pj1,k,l,m,pni)*&
                        conjg(h_in(:,:,k,l,m,n))
                end do
             end do
          end do
       end do
    else
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      cfgamma_h(:,j,k,l,m,n)=spec(n)%temp*spec(n)%temp_prof(li1:li2)*spec(n)%dens/&
                           fm(li1:li2,pj1,k,l,m,n)*conjg(h_in(:,j,k,l,m,n))
                   end do
                end do
             end do
          end do
       end do
    endif
  end subroutine get_cfgamma_h

  !>This subroutine gets chi from g and h
  subroutine get_chi_from_gh(g_in,h_in,chi_out)
    Implicit None
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in) :: h_in    
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: chi_out
    integer :: k,l,m,n,pni
    pni=pn1

    if (xy_local) then
       do n=ln1,ln2
          !if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   chi_out(:,:,k,l,m,n)=(spec(n)%temp)/(spec(n)%charge*fm(pi1,pj1,k,l,m,n))*&
                        (h_in(:,:,k,l,m,n)-g_in(:,:,k,l,m,n))
                end do
             end do
          end do
       end do
    else
       stop "get_chi_from_gh not ready for global!"
    endif
  end subroutine get_chi_from_gh

  !>This subroutine returns the additional prefactor to the RHS which is required for antenna runs
  subroutine get_cfgamma_ant(cfgamma_ant)
    Implicit None
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: cfgamma_ant
    complex, dimension(lk1:lk2):: A_ant
    real:: denom, sum0
    integer :: i,j,k,l,m,n,pni,n2,p,ierr
    real, dimension(li1:li2, lj1:lj2, lk1:lk2, 0:n_spec-1):: VIntCorr
    pni=pn1

    if (xy_local) then
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
                   VIntCorr(i,j,k,n) = sum0
                End Do
             End Do
          End Do
       enddo
       
       ! broadcast result to all other species
       Do n=0, n_procs_s-1
          Call mpi_bcast(&
               VIntCorr(li1,lj1,lk1,n*ln0), li0*lj0*lk0*ln0,&
               MPI_REAL_TYPE, n, mpi_comm_spec, ierr)
       End Do
       cfgamma_ant=(0.,0.)
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=ll1,ll2
                do p=1,n_antennae
                   call get_Apar_antenna(p,A_ant,i,j)
                   if (all((/j.ge.lj1,j.le.lj2/))) then
                      do k=lk1,lk2
                         denom=kperp2(i,j,k)
                         do n2=0,n_spec-1
                            denom=denom+beta*spec(n2)%charge**2*spec(n2)%dens/spec(n2)%mass*VIntCorr(i,j,k,n2)!
                         enddo
                         cfgamma_ant(i,j,k,l,m,n)=cfgamma_ant(i,j,k,l,m,n)+&
                              spec(n)%charge*spec(n)%dens*sqrt(2.*spec(n)%temp/spec(n)%mass)*vp(l)*conjg(A_ant(k))&
                              *get_jfac(i,j,k,m,n)*a11det_inv_antenna(i,j,k,1)
                         if (j+ky0_ind.eq.0) then
                            denom=kperp2(nx0-i,j,k)
                            do n2=0,n_spec-1
                               denom=denom+beta*spec(n2)%charge**2*spec(n2)%dens/spec(n2)%mass*VIntCorr(nx0-i,j,k,n2)!
                            enddo                            
                            cfgamma_ant(nx0-i,j,k,l,m,n)=cfgamma_ant(nx0-i,j,k,l,m,n)+&
                                 spec(n)%charge*spec(n)%dens*sqrt(2.*spec(n)%temp/spec(n)%mass)*vp(l)*A_ant(k)&
                                 *get_jfac(nx0-i,j,k,m,n)*a11det_inv_antenna(nx0-i,j,k,1)
                         endif
                      enddo
                   endif
                end do
             enddo
          end do
       end do
    else
       stop 'antenna not available for global simulations'
    endif
  end subroutine get_cfgamma_ant



  !>This subroutine gets the terms of interest in the energy equation.
  subroutine get_energy_terms(g_1_in,e_tot,dE_terms)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: g_1_in
    real :: e_tot(11)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms

    ! Local variables
    complex, allocatable, dimension(:,:,:,:,:,:) :: g_temp,rhs,cfgamma,f_in,h_in,cfgamma_g,cfgamma_ant
    complex, allocatable, dimension(:,:,:) :: rhs_block,g_block
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: fields_in
    complex,dimension(lk1:lk2):: dApar_dt, A_ant
    logical :: calc_this_term(11)
    integer :: lbg1,lbg2,klmn,iblock,iterm,i,j,k,p

!    if (nblocks.ne.1) stop 'get_energy_terms is currently only working with nblocks=1'
    lbg0 = lklmn0 / nblocks

    !    write(*,*) "Initialization"
    allocate(rhs(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(cfgamma(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(cfgamma_g(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    if (any(antenna_type.eq.(/2,3/))) &
         allocate(cfgamma_ant(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(f_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    if (arakawa_zv) then 
       allocate(h_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    else
       !dummy array that is passed to routines, but not used
       allocate(h_in(1:1,1:1,1:1,1:1,1:1,1:1))
    endif

    e_tot= 0.0 
    dE_terms = (0.0,0.0)
    calc_this_term(1:2) = .false.
   
    if (arakawa_zv) then 
       call calc_aux_fields(g_1_in,fields_in,f_in,.true.,h_in) !compute f and h
       call get_cfgamma_h(h_in,cfgamma)
    else
       call calc_aux_fields(g_1_in,fields_in,f_in,.false.) !compute f
       call get_cfgamma_f(f_in,fields_in,cfgamma)
    endif

    !for calculation of the drive term - independent of arakawa, etc.
    call get_cfgamma_g(g_1_in,cfgamma_g)

    ! ---------------------------- TOTAL ENERGY --------------------------
    ! IMPORTANT: Don't change the order!
    ! calling calc_rhs_only will update some pointers (e.g., ptr_bar_emfields)
    ! which are required for subsequent calls to, for instance, the nonlinearity

    !Energy
!    call get_energy_3d(g_1_in,dE_terms(:,:,:,1))
!    call get_energy_total(g_1_in,e_tot(1))
!!  the following should be faster (less memory demanding, respectively)
    if (any(antenna_type.eq.(/2,3/))) then
       !0.5 because cfgamma contains a factor 2 to account for the complex conjugate for the RHS terms
       rhs = 0.5*cfgamma*f_in(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    else
       rhs = 0.5*cfgamma*g_1_in   
    endif
    call energy_integral(rhs,dE_terms(li1:li2,lj1:lj2,lk1:lk2,1))

    !------ only for antenna ------!
    if (any(antenna_type.eq.(/2,3/))) then
       !add magnetic field contribution (including antenna!)
       dE_terms(:,:,:,1) = dE_terms(:,:,:,1)+kperp2/beta*abs(fields_in(:,:,lk1:lk2,2))**2
       !from here on, all cfgamma usage will be modified to include cfgamma_ant!
       call get_cfgamma_ant(cfgamma_ant)
       cfgamma=cfgamma+cfgamma_ant
    endif
    !------------------------------!

    !Total dE/dt
    !need g_temp because some compilers don't like sending an input variable 
    ! (g_1_in) to a subroutine that might modify it (calc_rhs_only)
    allocate(g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    g_temp=g_1_in
    rhs=(0.0,0.0)

    call calc_rhs_only(f_in,g_temp,fields_in,rhs,2,h_in)
    !collisions are taken out of calc_rhs_only in case of operator splitting
    !note that collisions exchange mu ghost cells
    if (coll_split) call equ_collisions(f_in,rhs,replace_rhs=.false.)

    deallocate(g_temp)

    rhs = cfgamma*rhs
    call energy_integral(rhs,dE_terms(:,:,:,2))
    select case(antenna_type)
       !NOTE: With antenna_type=2, there is the following external term to be taken into account in dE/dt.
       !As this term does not stem from the RHS, it appears to break the usual energy balance as written out in the file.
       !This is not actually the case, though, as it is simply an external modification of the energy evolution. 
       case(2)
       do p=1,n_antennae
          call get_Apar_antenna(p,A_ant,i,j)
          call get_dApar_dt_antenna(p,dApar_dt,i,j)
          if (all((/j.ge.lj1,j.le.lj2/))) then
             do k=lk1,lk2
                !factor 2 to account for complex conjugate
                dE_terms(i,j,k,2)=dE_terms(i,j,k,2)+2.*kperp2(i,j,k)/beta*conjg(A_ant(k))*dApar_dt(k)*&
                     a11det_inv_antenna(i,j,k,1)
                if (j+ky0_ind.eq.0.and.i.ne.0) &
                     dE_terms(nx0-i,j,k,2)=dE_terms(nx0-i,j,k,2)+2.*kperp2(nx0-i,j,k)/beta*conjg(A_ant(k))&
                     *conjg(dApar_dt(k))*a11det_inv_antenna(nx0-i,j,k,1)
             enddo
          endif
       enddo
       !NOTE: With antenna_type=3, the same comments as above apply, though the form of the additional term differs.
       case(3)
       do p=1,n_antennae
          call get_dApar_dt_antenna(p,dApar_dt,i,j)
          if (all((/j.ge.lj1,j.le.lj2/))) then
             do k=lk1,lk2
                !factor 2 to account for complex conjugate
                dE_terms(i,j,k,2)=dE_terms(i,j,k,2)+2.*kperp2(i,j,k)/beta*conjg(fields_in(i,j,k,2))*dApar_dt(k)
                if (j+ky0_ind.eq.0.and.i.ne.0) &
                     dE_terms(nx0-i,j,k,2)=dE_terms(nx0-i,j,k,2)+2.*kperp2(nx0-i,j,k)/beta*&
                     conjg(fields_in(nx0-i,j,k,2))*conjg(dApar_dt(k))
             enddo
          endif
       enddo
       case default
       end select

    !------------------ INDIVIDUAL TERMS -------------------------------
    ! ---------------------- nonblocked terms ---------------------------------
    !collisions
    if (collision_op.ne.'none') then
       call equ_collisions(f_in,rhs,replace_rhs=.true.)
       rhs = cfgamma*rhs
       call energy_integral(rhs,dE_terms(:,:,:,5))
    endif
    deallocate(rhs)

PERFON('en_block')
    !------------------- strip mining terms ---------------------------
    allocate(rhs_block(li1:li2,lj1:lj2,1:lbg0),g_block(lbi:ubi,lj1:lj2,1:lbg0))
    rhs_block = (0.0,0.0)
    Do iblock=1,nblocks
       lbg1 =(iblock-1)*lbg0+1
       lbg2 =lbg1+lbg0-1

       if (.not.xy_local) then
          if (n_fields.gt.1) then
             do klmn=lbg1,lbg2
                chi_block(li1:li2,lj1:lj2,klmn-lbg1+1) = &
                     &ptr_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1) &
                     &- vTvpar(sl(klmn),sn(klmn)) * &
                     & ptr_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),2) 
             end do
             call exchange_x(x_boundary_block,chi_block)
          else
             ! here chi is independent of v_par, therefore one exchange for each j,k,m,n
             ! is enough
             do klmn=lbg1,lbg2
                chi_block(li1:li2,lj1:lj2,klmn-lbg1+1) = &
                     &ptr_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1)
             end do
             call exchange_x(x_boundary_block,chi_block)
          end if
       end if

       do iterm=3,10
          calc_this_term(iterm) = .true.
          select case(iterm) 
             ! ---------------- SOURCE/DRIVE TERMS -------------------
          case(3)
             if (.not.only_neo) then
                if (xy_local) then
                   call add_dchidxy_orig(fields_in,ptr_bar_emfields,&
                   &rhs_block,ptr_barchi,ptr_dbarchidxy,lbg1,lbg2,.false.)
                else
                   call add_dchidxy(chi_block,rhs_block,ptr_dbarchidxy,lbg1,lbg2,.false.)
                end if
             else
                if (rhs_f0) call add_f0_term(rhs_block,lbg1,lbg2,time)
             endif
             if (antenna_type.eq.3.and.antenna_contrib) &
                  call add_dApar_dt_antenna(rhs_block,lbg1,lbg2)
          case(4)
             if (.not.xy_local) then
                call add_f1_sources_block(rhs_block,lbg1,lbg2,0)
                if (explicit_buffer) call add_kBuffer_explicit(f_in,rhs_block,lbg1,lbg2)
             else
                calc_this_term(iterm) = .false.
             endif
             ! ---------------- SINK TERMS -------------------
          case(5)
              ! (note: collisions are not included here but later)  
             calc_this_term(iterm) = .false.
          case(6)
             !hyp_v
             if (arakawa_zv) then 
                if (hyp_on_h) then
                   call add_hypv_ak(h_in,rhs_block,lbg1,lbg2)
                else
                   call add_hypv_ak(f_in,rhs_block,lbg1,lbg2)
                endif
             else
                call add_hypv_block(f_in,rhs_block,lbg1,lbg2)
             endif
          case(7)
             !hyp_z
             if (arakawa_zv) then
                if (hyp_on_h) then
                   call add_hypz_ak(h_in,rhs_block,lbg1,lbg2)  
                   if (hypz_compensation) call equ_comp_hypz(fields_in,ptr_bar_emfields,rhs_block,lbg1,lbg2)
                else
                   call add_hypz_ak(f_in,rhs_block,lbg1,lbg2) 
                endif
             else
                call add_hypz_block(f_in,rhs_block,lbg1,lbg2)
             endif
          case(8)
             if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then
                if (arakawa_zv) then
                   if (hyp_on_h) then
                      call add_dfdxy(h_in,rhs_block,lbg1,lbg2) 
                   else
                      call add_dfdxy(f_in,rhs_block,lbg1,lbg2) 
                   endif                
                else
                    call add_dfdxy(f_in,rhs_block,lbg1,lbg2) 
                endif
             else
                calc_this_term(iterm) = .false.
             endif
             ! ----------------- Analytically vanishing terms -----------------
          case(9)
             if (rhs_nl) then
                if (.not.xy_local) then
                   ptr_barchi(:,:,:) = chi_block(li1:li2,lj1:lj2,:)
                   do klmn=lbg1,lbg2
                      g_block(li1:li2,:,klmn-lbg1+1)=g_1_in(:,:,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
                   enddo
                   call exchange_x(x_boundary_block,g_block)
                   call add_dgdxy(g_block, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
                   rhs_block = 0.0
                   
                   CALL this_nonlinear_term%add(g_block, ptr_dgdxy, fields_in,&
                        &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, 2)

                   !CALL add_nonlinearity(g_block, ptr_dgdxy, fields_in, &
                   !     &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, 2)
                else
                   if (.not.nonlin_h) then
                      CALL this_nonlinear_term%add(g_1_in, ptr_dgdxy, fields_in,&
                           &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, 2)
                   else
                      CALL this_nonlinear_term%add(h_in, ptr_dgdxy, fields_in,&
                           &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, 2)
                   endif
                endif
             else
                calc_this_term(iterm) = .false.
             endif
          case(10)
             ! Poisson bracket in z,vpar, i.e.
             ! cfgamma * (dfdzv+dchidzv) 
             !!this includes the (z,vpar) hyperdiffusion if hyp_on_h=T
             if (arakawa_zv) then
                call equ_dzv(h_in,rhs_block,lbg1,lbg2)
                if (hyp_on_h) then
                   if (hypz_compensation) call equ_comp_hypz(fields_in,ptr_bar_emfields,rhs_block,lbg1,lbg2)
                endif
             else
                call equ_dfdzv(f_in,rhs_block,lbg1,lbg2)
                call add_dchidz(fields_in, ptr_bar_emfields, rhs_block, lbg1, lbg2)
             end if
          case default
             stop 'wrong index in blockterm loop'
          end select

          if (calc_this_term(iterm)) then
             if(iterm==3.and.(.not.only_neo).and.(.not.any(antenna_type.eq.(/2,3/)))) then 
                !if this is the drive 1 term then we need only cfgamma_g instead of all cfgamma
               call block_mult(rhs_block,cfgamma_g,lbg1,lbg2)
             else
               call block_mult(rhs_block,cfgamma,lbg1,lbg2)
             end if
             call add_energy_integral_block(rhs_block,dE_terms(:,:,:,iterm),lbg1,lbg2)
             rhs_block = (0.0,0.0)
          endif
       enddo
    enddo
    deallocate(rhs_block,g_block)
PERFOFF

    do iterm=1,10
       if (calc_this_term(iterm)) &
            &call my_complex_sum_vwspec(dE_terms(:,:,:,iterm),size(dE_terms(:,:,:,iterm)))
       !------------ direct space integration --------------------       
       if ((p_has_00_mode).and.(xy_local).and.(.not.only_neo)) dE_terms(li1,lj1,lk1:lk2,iterm)=(0.0,0.0)
       call sum_int_3d(dE_terms(:,:,:,iterm),e_tot(iterm))
       !compute dEdt_NC
       ! Note: nonlinearity should not be included
       !Warning, z,v poisson bracket does not cancel k by k and can also include dissipative boundary condition.
       if (iterm.gt.2.and.iterm.ne.9) dE_terms(:,:,:,11)=dE_terms(:,:,:,11)+dE_terms(:,:,:,iterm)
    enddo

    if ((arakawa_zv.and.hyp_on_h).or.(.not.arakawa_zv)) then
       e_tot(10)= e_tot(10)-e_tot(6)-e_tot(7) !eliminated hyperdiffusions
       dE_terms(:,:,:,11)=dE_terms(:,:,:,11)-dE_terms(:,:,:,6)-dE_terms(:,:,:,7) !eliminated hyperdiffusions
       !dE_terms(:,:,:,10)=dE_terms(:,:,:,10)-dE_terms(:,:,:,6)-dE_terms(:,:,:,7) !eliminated hyperdiffusions
    !else !hyperdiffusions are not included in parallel dynamics scheme
    endif
    !rest
    e_tot(11)=e_tot(2)-sum(e_tot(3:10))

#if 0
!now the other way round (just for testing)
    allocate(rhs(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    g_temp=g_1_in
    rhs=(0.0,0.0)   
    call calc_rhs_only(f_in,g_temp,fields_in,rhs,2,h_in)
    !collisions are taken out of calc_rhs_only in case of operator splitting
    if(coll_split) call equ_collisions(f_in,rhs,replace_rhs=.false.)

    deallocate(g_temp)

    if (arakawa_zv) then
       call calc_aux_fields(rhs,fields_in,f_in,.true.,h_in)
       call get_cfgamma_h(h_in,cfgamma)
    else
       call calc_aux_fields(rhs,fields_in,f_in,.false.)
       call get_cfgamma_f(f_in,fields_in,cfgamma)
    endif

!linear in dt

    call energy_integral(cfgamma*g_1_in,erhs1)
    if((p_has_00_mode).and.(xy_local).and.(.not.only_neo)) erhs1(li1,lj1,:)=(0.0,0.0)
    call sum_int_3d(erhs1,erg)
!quadratic in dt

    call energy_integral(cfgamma*rhs,erhs1)
    if((p_has_00_mode).and.(xy_local).and.(.not.only_neo)) erhs1(li1,lj1,:)=(0.0,0.0)
    call sum_int_3d(erhs1,erg2)

if(mype.eq.0) then
   write (*, '(4(ES14.6,X))') (erg-e_tot(2))/e_tot(2),(erg+e_tot(2))*dt, &
        &(erg+e_tot(2))*dt+erg2*dt**2, e_tot(1)+(erg+e_tot(2))*dt+erg2*dt**2
endif
    deallocate(rhs)
#endif

    deallocate(cfgamma,f_in)
    if (any(antenna_type.eq.(/2,3/))) deallocate(cfgamma_ant)
    deallocate(h_in)

  end subroutine get_energy_terms


  subroutine block_mult(rhs_block, cfgamma, lb1, lb2)
    integer, intent(in) :: lb1, lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout) :: rhs_block
    complex, dimension(li1:li2,lj1:lj2,1:lklmn0), intent(in) :: cfgamma

    rhs_block(:,:,lb1:lb2) = rhs_block(:,:,lb1:lb2) * cfgamma(:,:,lb1:lb2)

  end subroutine block_mult


  !>Driver for energy diagnostics 
  !!
  subroutine diag_energy(write_1d,write_3d)
    logical :: write_1d,write_3d
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2) :: dE_diss
    real :: e_tot(11)
    real :: cancel_check

PERFON('dia_en')

    call get_energy_terms(g_1,e_tot,dE_terms)

    if (write_3d) then
       if(write_std) then
          dE_diss = dE_terms(:,:,:,6)+dE_terms(:,:,:,7)+dE_terms(:,:,:,8)
          call energy_write_3d(dE_terms(:,:,:,1),dE_terms(:,:,:,11),dE_terms(:,:,:,3),&
               &dE_terms(:,:,:,5),dE_diss,&
               dE_terms(:,:,:,9),ENERGYFILE) !add dE_heat ??
       end if

#ifdef WITHFUTILS
       if(write_h5) then
          isnap_energy3d = itime/istep_energy3d + 1 
          if ((my_pev+my_pew+my_pespec).eq.0) then
             call append(fidenergy3d_h5, "/energy3d/time", real(time,8))
             call attach(fidenergy3d_h5, "/energy3d/time", "n_steps", isnap_energy3d)
             
             write(dset_name_energy3d, "(A, '/', i10.10)") "/energy3d/var1", isnap_energy3d
             call putarrnd(fidenergy3d_h5, dset_name_energy3d, dE_terms(li1:li2,lj1:lj2,lk1:lk2,1), (/3, 2, 1/))
             call attach(fidenergy3d_h5, dset_name_energy3d, "time", time)

             write(dset_name_energy3d, "(A, '/', i10.10)") "/energy3d/var2", isnap_energy3d
             call putarrnd(fidenergy3d_h5, dset_name_energy3d, dE_terms(li1:li2,lj1:lj2,lk1:lk2,11), (/3, 2, 1/))
             call attach(fidenergy3d_h5, dset_name_energy3d, "time", time)
             
             write(dset_name_energy3d, "(A, '/', i10.10)") "/energy3d/var3", isnap_energy3d
             call putarrnd(fidenergy3d_h5, dset_name_energy3d, dE_terms(li1:li2,lj1:lj2,lk1:lk2,3), (/3, 2, 1/))
             call attach(fidenergy3d_h5, dset_name_energy3d, "time", time)
             
             write(dset_name_energy3d, "(A, '/', i10.10)") "/energy3d/var4", isnap_energy3d
             call putarrnd(fidenergy3d_h5, dset_name_energy3d, dE_terms(li1:li2,lj1:lj2,lk1:lk2,5), (/3, 2, 1/))
             call attach(fidenergy3d_h5, dset_name_energy3d, "time", time)
             
             isnap_energy3d = isnap_energy3d+1
             call flushh5(fidenergy3d_h5)
          end if
       end if
#endif
    endif

    if (write_1d) then
       if(abs(e_tot(3)).gt.1.0e-12) then
          cancel_check=(abs(e_tot(9))+abs(e_tot(10))+&
               abs(e_tot(11)))/abs(e_tot(3)) !non-conservation (normalized to drive)
       else
          cancel_check=0.0
       end if

       if(write_std) then
          if(mype==0) then
             write(ETERMSFILE,"(14ES12.4)") time, &
                  e_tot(1),&  !total energy
                  e_tot(2),&  !total dE/dt
                  e_tot(3),&  !Drive
                  e_tot(4),&  !Heat sources (only nonlocal)
                  e_tot(5),&  !Collisions 
                  e_tot(6),&  !hyp_v
                  e_tot(7),&  !hyp_z
                  e_tot(8),&  !hyp_x/y
                  e_tot(9),&  !nonlinearity
                  e_tot(10),& ![Poisson error]_z,v
                  e_tot(11),& !curvature (+rest, if any)
                  cancel_check,& !Test of cancellation
                  (e_tot(1)-last_energy)/(dt*istep_energy) !Actual exact change in energy (2.0 because of c.c)
             call flush(ETERMSFILE)
          end if
       end if
       
#ifdef WITHFUTILS
       if(write_h5) then
          isnap_energy = itime/istep_energy + 1 
          if ((mype).eq.0) then
             call attach(fidenergy_h5, "/energy", "n_steps", isnap_energy)    
             call add_record(hbufenergy_h5, "time", "simulation time", real(time,8))
             call add_record(hbufenergy_h5, "Total energy", "E_{tot}", real(e_tot(1),8))
             call add_record(hbufenergy_h5, "Total dE_dt", "dE_dt", real(e_tot(2),8))
             call add_record(hbufenergy_h5, "Drive 1 term", "T_\par", real(e_tot(3),8))
             call add_record(hbufenergy_h5, "Heat sources (nonlocal)", "Q_{sources}", real(e_tot(4),8))
             call add_record(hbufenergy_h5, "Collisional Dissipation", "coll", real(e_tot(5),8))
             call add_record(hbufenergy_h5, "hyp_v dissipation", "hyp_v", real(e_tot(6),8))
             call add_record(hbufenergy_h5, "hyp_z dissipation", "hyp_z", real(e_tot(7),8))
             call add_record(hbufenergy_h5, "hyp_x_y dissipation", "hyp_{x_y}", real(e_tot(8),8))
             call add_record(hbufenergy_h5, "nonlinearity", "hyp_z", real(e_tot(9),8))
             call add_record(hbufenergy_h5, "[Poisson err]_z,vpar", "[Poisson err]_z,v", real(e_tot(10),8))
             call add_record(hbufenergy_h5, "curvature + remaining terms", "curvature + rem terms", real(e_tot(11),8))
             call htable_endstep(hbufenergy_h5)
             call htable_hdf5_flush(hbufenergy_h5)
             call flushh5(fidenergy_h5)
          end if
       end if
#endif
    endif
    last_energy=e_tot(1)

PERFOFF
  end subroutine diag_energy

  !>energy_exchange diagnostics version 2
  !! See, e.g., J. Candy PHYSICS OF PLASMAS 20, 082503 (2013)
  subroutine diag_energy_exchange
    ! Local variables
    complex, allocatable, dimension(:,:,:,:,:,:) :: g_temp,rhs,chi_rhs,chi,f_temp,h_temp
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: fields
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: fields_rhs
    integer :: ierr,n
    real :: energy_exchange_loc(0:n_spec-1),energy_exchange(0:n_spec-1)
    complex :: energy_exchange_3d(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2)


    allocate(rhs(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(f_temp(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    allocate(h_temp(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    allocate(chi_rhs(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    energy_exchange_loc=0.0
    energy_exchange=0.0
     
    g_temp = g_1
    !get h
    call calc_aux_fields(g_temp,fields,f_temp,.true.,h_temp) 
    !get chi
    call get_chi_from_gh(g_1,h_temp,chi)
    !get dh/dt
    call calc_rhs_only(f_temp,g_temp,fields,rhs,2,h_temp)
    !collisions are taken out of calc_rhs_only in case of operator splitting
    !note that collisions exchange mu ghost cells
    if (coll_split) call equ_collisions(f_temp,rhs,replace_rhs=.false.)
    !now get dchi/dt
    call calc_aux_fields(rhs,fields_rhs,f_temp,.true.,h_temp) !note: use chi for h
    call get_chi_from_gh(rhs,h_temp,chi_rhs)

    !put it together in h_temp    
    do n=ln1,ln2
      h_temp(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n)=0.5*spec(n)%charge*( &
           -rhs(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n)*conjg(chi(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n)) &
           +conjg(chi_rhs(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n))*g_1(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n))
    end do

    do n=ln1,ln2
      call energy_integral_vw(h_temp(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n),energy_exchange_3d(:,:,:,n))
      call sum_int_3d(energy_exchange_3d(:,:,:,n),energy_exchange_loc(n))
    end do

    call mpi_allreduce(energy_exchange_loc,energy_exchange,n_spec,MPI_REAL_TYPE,MPI_SUM,mpi_comm_spec,ierr)

    !put it together in h_temp    
    !energy_exchange_loc=0.0
    !energy_exchange2=0.0
    !do n=ln1,ln2
    !  h_temp(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n)=0.5*spec(n)%charge*( &
    !       +chi_rhs(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n)*conjg(g_1(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n)))
    !end do

    !do n=ln1,ln2
    !  call energy_integral_vw(h_temp(:,:,lk1:lk2,ll1:ll2,lm1:lm2,n),energy_exchange_3d(:,:,:,n))
    !  call sum_int_3d(energy_exchange_3d(:,:,:,n),energy_exchange_loc(n))
    !end do

    !call mpi_allreduce(energy_exchange_loc,energy_exchange2,n_spec,MPI_REAL_TYPE,MPI_SUM,mpi_comm_spec,ierr)

    !if(write_std) then
    !   if(mype==0) then
    !      write(ENEXCHANGEFILE2,"(5ES12.4)") time, energy_exchange, energy_exchange2
    !      call flush(ENEXCHANGEFILE2)
    !   end if
    !end if

    if(write_std) then
       if(mype==0) then
          write(ENEXCHANGEFILE,"(3ES12.4)") time, energy_exchange 
          call flush(ENEXCHANGEFILE)
       end if
    end if

    deallocate(rhs)
    deallocate(g_temp)
    deallocate(f_temp)
    deallocate(h_temp)
    deallocate(chi_rhs)
    deallocate(chi)

  end subroutine diag_energy_exchange

  !>Gathers an input 3D array and outputs as unformatted data to file determined by filehandle.
  !!
  subroutine energy_write_3d(energy,dedt_nc,eg_flux,eg_coll,eg_diss,eg_nl,filehandle)
    !****
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2),intent(in) :: energy,dedt_nc,eg_flux,&
                 eg_coll,eg_diss,eg_nl
    integer :: filehandle

    ! Local variables
    real :: ar3dr(li1:li2,lj1:lj2,lk1:lk2)
    integer:: ierr 
    real, Dimension(0:ni0-1, lj1:lj2,  lk1:lk2):: energyx
    real, Dimension(0:ni0-1, 0:nj0-1, lk1:lk2):: energyf
    real, Dimension(0:ni0-1, 0:nj0-1, 0:nz0-1)::  dest3d
    Integer :: j,k,write_pe,which_eg=1,num_out
    !COMPLEX :: ar3dc(li1:li2,lj1:lj2,lk1:lk2),e_tot(16)

    num_out=6

    write_pe=0

    IF (mype.eq.0) write(filehandle) time

    do which_eg=1,num_out

       energyx=0.0
       energyf=0.0
       dest3d=0.0

       select case(which_eg) 
       case(1)
          ar3dr=real(energy)
       case(2)
          ar3dr=real(dedt_nc)
       case(3)
          ar3dr = real(eg_flux)
       case(4)
          ar3dr = real(eg_coll)
       case(5)
          ar3dr = real(eg_diss)
       case(6)
          ar3dr = real(eg_nl)
       end select

       if ((my_pev+my_pew+my_pespec).eq.0) then
          ! Gather energy over x-Distribution on pex==0.
          Do k = lk1, lk2
             Do j = lj1, lj2
                Call mpi_gather(ar3dr(li1,j,k), li0, MPI_REAL_TYPE,&
                     energyx(0,j,k), li0, MPI_REAL_TYPE,&
                     0, mpi_comm_x, ierr)
             Enddo
          Enddo

          ! Gather energy over y-Distribution on pey==0.
          Do k = lk1, lk2
             Call mpi_gather(energyx(0,lj1,k), Size(energyx(:,:,k)), MPI_REAL_TYPE,&
                  energyf(0,0,k), Size(energyx(:,:,k)), MPI_REAL_TYPE,&
                  0, mpi_comm_y, ierr)
          Enddo

          ! Gather energyf over z-Distribution on pez=0
          CALL mpi_gather(energyf(0,0,lk1),SIZE(energyf(:,:,:)), MPI_REAL_TYPE,&
               dest3d(0,0,0),SIZE(energyf),MPI_REAL_TYPE,&
               0, mpi_comm_z, ierr)

          IF (mype.eq.write_pe) write(filehandle) dest3d

       endif
    end do !which_eg
    !     IF (mype.eq.0) call flush(filehandle)

  end subroutine energy_write_3d

  !>This subroutine calculates the energy as a function of kx,ky,z
  !!Note: the kx=0,ky=0 mode is zeroed
  subroutine get_energy_3d(g_1_in,energy)
    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(out) :: energy
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields) :: fields_in
    complex, allocatable, dimension(:,:,:,:,:,:) :: f_in,temp

    allocate(f_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    allocate(temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))

    energy=(0.0,0.0)
    CALL calc_aux_fields(g_1_in,fields_in,f_in,.false.)

    !Get energy
    call get_cfgamma_f(f_in,fields_in,temp)
    temp=temp*g_1_in
    call energy_integral(temp,energy)
    if((p_has_00_mode).and.(xy_local).and.(.not.only_neo)) energy(li1,lj1,:)=(0.0,0.0)

    deallocate(f_in,temp)

  end subroutine get_energy_3d

  !>This subroutine calculates the total energy
  !!by simply calling get_energy_3d and summing / integrating over kx,ky,z
  subroutine get_energy_total(g_1_in,e_tot)
    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: energy
    real, intent(out) :: e_tot

    energy=(0.0,0.0)
    e_tot= 0.0

    !Get energy
    call get_energy_3d(g_1_in,energy)
    call sum_int_3d(energy,e_tot)

  end subroutine get_energy_total

End Module diagnostics_energy
