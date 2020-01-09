#include "redef.h"
#include "intrinsic_sizes.h"
#undef VOL_AVG
!>Contains special diagnostics for the radially nonlocal code
Module diagnostics_df
  Use par_mod
  Use file_io
  Use vel_space
  Use communications
  Use boundaries
  Use aux_fields
  use geometry
  use sources_mod
  use f0_term, only: add_f0_term_df_1d
  use coordinates
  use lagrange_interpolation, only: lag3deriv
  use diagnostics_neoclass, only: set_mats_diag_neoclass,istep_neoclass
  use profile_smoothing

  Implicit None

  !The basic diag interface
  PUBLIC :: initialize_all_diags_df, exec_all_diags_df, &
       &finalize_all_diags_df,mem_est_all_diags_df, &
       &check_diag_df, set_mats_diag_profile,&
       &istep_prof,avgprof_stime, avgprof,&
       &istep_srcmom

  PRIVATE 

  REAL, DIMENSION(:,:,:), ALLOCATABLE,public :: new_profs
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: avgprof
  REAL :: avgprof_stime = -1.0

  Character(Len=8), PRIVATE:: filestat='replace', filepos='rewind'

  INTEGER, DIMENSION(:),allocatable,private:: PROFFILE, SRCMOMFILE
  Logical :: write_prof_pe, write_srcmom_pe

  REAL, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE,private :: prof_mat, srcmom_mat
  integer::n_prof_mats = 6
  integer::n_nc_mats = 4
  REAL :: old_prof_time = -1.0
 
  Integer, private :: block_nr
  Integer :: istep_prof=-1, istep_srcmom=-1
  character(len=11):: srcmom_format

Contains

  Subroutine check_diag_df

    if (.not.xy_local.and.istep_neoclass.gt.0.and.istep_prof.le.0) then
       istep_prof = istep_neoclass
       if (mype==0) Write(*,"(A)") &
            &'Switch on diag_df (istep_prof = istep_neoclass) for profile information'
    endif

    if (istep_prof.gt.0) then
       if (xy_local) then
          istep_prof = 0
          if (mype==0) print*, &
               &'istep_prof > 0 not possible in local simulations'
       elseif (ky0_ind.gt.0) then
          istep_prof = 0
          if (mype==0) print*, &
               &'istep_prof > 0 not possible in linear, ky>0 simulations'
       elseif (.not.y_local) then
          istep_prof = 0
          if (mype==0) print*, &
               &'istep_prof > 0 not implemented yet for y global simulations'
       endif
    endif

  End Subroutine check_diag_df

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_all_diags_df(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc = 0.0
    if (istep_prof.gt.0) mem_loc = mem_est_diag_profile(mem_loc)
    if (istep_srcmom.gt.0) mem_loc = mem_est_diag_srcmom(mem_loc)

    mem_est_all_diags_df=mem_req_in+mem_loc
  End Function mem_est_all_diags_df

!*****************************************************************

  Subroutine initialize_all_diags_df
    if (istep_prof.gt.0) call initialize_diag_profile
    if (istep_srcmom.gt.0) call initialize_diag_srcmom
  End Subroutine Initialize_all_diags_df

  SUBROUTINE exec_all_diags_df(itime,time,reset)
    integer, intent(in) :: itime
    real, intent(in) :: time
    logical, intent(inout) :: reset
    !logical, intent(inout) :: overflow_exit, underflow_exit, reset
    
    IF (istep_prof.GT.0) THEN
       IF (MODULO(itime,istep_prof) == 0) THEN
          call calc_aux_fields(g_1,emfields,f_,.false.)
          CALL diag_profile(reset)
       ENDIF
    END IF
    IF (istep_srcmom.GT.0) THEN
       IF (MODULO(itime,istep_srcmom) == 0) THEN
          CALL diag_srcmom(time)
       ENDIF
    END IF
  END SUBROUTINE exec_all_diags_df

  Subroutine finalize_all_diags_df
    if (istep_prof.gt.0) call finalize_diag_profile
    if (istep_srcmom.gt.0) call finalize_diag_srcmom
  End Subroutine finalize_all_diags_df

!*****************************************************************
  !>Give an estimate of the memory requirements of this diag
  Real Function mem_est_diag_profile(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    
    !new_profs
    mem_loc = nx0*ln0*4*SIZE_OF_REAL_MB
    !prof_mat
    mem_loc = mem_loc + pi0*pj0*lk0*(n_prof_mats+n_nc_mats)*ll0*lm0*ln0*SIZE_OF_REAL_MB
    !local vars
    !momc
    mem_loc = mem_loc + lijk0*n_prof_mats*ln0*SIZE_OF_COMPLEX_MB
    !momc_neo
    mem_loc = mem_loc + lijk0*n_nc_mats*ln0*SIZE_OF_COMPLEX_MB
    !mom
    mem_loc = mem_loc + lijk0*3*SIZE_OF_REAL_MB
    !mom_neo
    mem_loc = mem_loc + 4*pi0*lk0 * SIZE_OF_REAL_MB
    !dfields_dy, tmp_x, tmp_dx
    mem_loc = mem_loc + (lijk0*n_fields+lx0+li0)*SIZE_OF_COMPLEX_MB

    mem_est_diag_profile=mem_req_in+mem_loc
  End Function mem_est_diag_profile

!*****************************************************************

  SUBROUTINE initialize_diag_profile
    Implicit None
    Integer :: n

    ALLOCATE(PROFFILE(ln1:ln2))
    allocate(new_profs(0:nx0-1,0:n_spec-1,1:4))

    block_nr=0

    write_prof_pe = (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec))

    if (write_prof_pe) then
       ! Moment files for each species in a different file
       DO n=ln1,ln2
          call get_unit_nr(PROFFILE(n))
          OPEN(PROFFILE(n), file=trim(diagdir)//'/profile_'//&
               &trim(spec(n)%name)//''//trim(file_extension),&
               form='formatted', status=filestat, position=filepos)  
          ! write file header
          write(PROFFILE(n),"(3A)") "#   x/a             x/rho_ref       T/Tref          ",&
                "n/nref             omt             omn            Gamma            Q         Pi        ",&
                "Gamma_neo          Q_neo        Pi_neo         j_boot"
       END DO
    endif

    Allocate(prof_mat(pi1:pi2,pj1:pj2,lk1:lk2,1:(n_prof_mats+n_nc_mats),ll1:ll2,lm1:lm2,ln1:ln2))

    IF (.NOT.ALLOCATED(avgprof)) ALLOCATE(avgprof(0:nx0-1,0:n_spec-1,0:7))
    avgprof=0.0
    old_prof_time = -1

    Call set_mats_diag_profile

    !Set default smoothing width for resets here, as it depends on rhostar
    if (smooth_reset .and. smooth_reset_width .LE. 0) then
        smooth_reset_width = 5*rhostar
    endif

  END SUBROUTINE initialize_diag_profile


!*************************************************

  SUBROUTINE set_mats_diag_profile
    Implicit None
    Integer :: k,l,m,n
    real, dimension(:,:,:,:,:,:,:), allocatable::neoclass_mat

    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO n=ln1,ln2
                !mat: first index vpar, second vperp (!)
                !mat_00
                prof_mat(:,:,k,1,l,m,n) = mat_00(:,:,k,l,m)
                ! mat_20 + mat_02
                prof_mat(:,:,k,2,l,m,n) = (vp(l)*vp(l) + mu(m)*geom%Bfield(pi1:pi2,:,k)) &
                     & * mat_00(:,:,k,l,m)
                !mat_10 
                prof_mat(:,:,k,3,l,m,n) = vp(l) * mat_00(:,:,k,l,m)
                if (n_fields.gt.1) then
                   !v_T(x_0) * mat_10
                   prof_mat(:,:,k,4,l,m,n) = sqrt(2.0*spec(n)%temp/spec(n)%mass) * vp(l) * &
                        &mat_00(:,:,k,l,m)
                   !v_T(x_0) * (mat_30 + mat_12)
                   prof_mat(:,:,k,5,l,m,n) = (vp(l)*vp(l) + mu(m)*geom%Bfield(pi1:pi2,:,k)) * &
                        &prof_mat(:,:,k,4,l,m,n)
                   !v_T(x0) * (mat_20)
                   prof_mat(:,:,k,6,l,m,n) = vp(l) * prof_mat(:,:,k,4,l,m,n)
                endif
             END do
          END DO
       END DO
    END DO    
    allocate(neoclass_mat(pi1:pi2,pj1:pj2,lk1:lk2,1:4,ll1:ll2,lm1:lm2,ln1:ln2))
    call set_mats_diag_neoclass(neoclass_mat,4)
    prof_mat(:,:,:,7:10,:,:,:)=neoclass_mat

    deallocate(neoclass_mat)

  END SUBROUTINE set_mats_diag_profile

  !> Profile diagnostic: Calculates temperature and density profiles
  !! 1.) total n / n_j0
  !! 2.) total T / T_j0
  !! 3.) L_ref/L_n_j
  !! 4.) L_ref/L_T_j
  !! 5.) Particle Flux (x) / (c_ref n_ref (rho_ref/Lref)^2)
  !! 6.) Heat flux (x) / (c_ref p_ref (rho_ref/Lref)^2)
  !! 7.) Momentum flux (x) / (c_ref^2 n_ref m_ref (rho_ref/Lref)^2)
  !! 8.) neoclassical Particle Flux (x) / (c_ref n_ref (rho_ref/Lref)^2)
  !! 9.) neoclassical Heat flux (x) / (c_ref p_ref (rho_ref/Lref)^2)
  !! 10.) neoclassical momentum flux (x) / (c_ref^2 n_ref m_ref (rho_ref/Lref)^2)
  !! 11.) neoclassical bootstrap current (x) / (c_ref Bref (rho_ref/Lref))
  Subroutine diag_profile(reset)
    Logical, intent(out) :: reset
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_prof_mats,ln1:ln2):: momc
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_nc_mats,ln1:ln2):: momc_neo
    Real,    Dimension(li1:li2,lj1:lj2,lk1:lk2,3):: mom
    real, Dimension(pi1:pi2,lk1:lk2,1:4):: mom_neo
    real, Dimension(li1:li2,1:4):: neoflux
    real, Dimension(0:nx0-1,1:4):: neoflux_out
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields) :: dfields_dy
    Real, Dimension(li1:li2,1:7) :: var
    Real, Dimension(0:nx0-1,1:7) :: out_arr
    Real, Dimension(0:nx0-1) :: temp_arr
    real,dimension(0:nx0-1,1:4):: tmpprof
    Integer:: i, j, k, n, o, ierr
    REAL, DIMENSION(li1:li2) :: fnorm
    logical:: loc_reset

    PERFON('d_prof')

    PERFON('d_prof1')

    reset = .false.

    momc=cmplx(0,0)
    momc_neo = cmplx(0,0)
    
    ! f_ has boundary points in all dimensions
    Call calc_vsp_moment(n_prof_mats,f_,.true.,emfields,prof_mat(:,:,:,1:6,:,:,:),momc,.false.)
    Call calc_vsp_moment(n_nc_mats,f_,.true.,emfields,prof_mat(:,:,:,7:10,:,:,:),momc_neo,nc_exb_corr) 

    PERFOFF
    
    do o=1,n_fields
       do k=lk1,lk2
          do j=lj1,lj2
             do i=li1,li2
                ! Definition of bx is NOT identical to B_{1,x} in the script, as
                ! the Bfield is part of bx but NOT part of B_{1,x} ! 
                dfields_dy(i,j,k,o)=imag*ky(j)*emfields(i,j,k,o)*flux_geomfac(i,pj1,k)
             enddo
          enddo
       enddo
    enddo
    
#ifndef VOL_AVG
    ! normalization
    fnorm = 1.0 / (REAL(nz0)*geom%avg_jaco_yz(li1:li2))
#else
    !switch on for volume averages (don't forget sum over x)
    fnorm = 1.0 / (REAL(nz0*nx0)*geom%avg_jaco)
#endif

    do n=ln1,ln2
       var = 0.0

       !Calculate flux surface average of n_1 and T_1
       if (p_has_0_mode) then
          do k=lk1,lk2
             var(li1:li2,1) = var(li1:li2,1)+&
                  &real(momc(li1:li2,lj1,k,1,n))*geom%jacobian(li1:li2,pj1,k)
             var(li1:li2,2) = var(li1:li2,2)+(real(2.0/3.0*momc(li1:li2,lj1,k,2,n))-&
                  &spec(n)%temp_prof(li1:li2)*real(momc(li1:li2,lj1,k,1,n)))*&
                  &geom%jacobian(li1:li2,pj1,k)/(spec(n)%dens_prof(li1:li2))
          enddo
       endif

       do o=1,2
          Call my_sum_to_0_real(var(:,o), li0, mpi_comm_vw)
          Call my_sum_to_0_real(var(:,o), li0, mpi_comm_z)
       enddo
       
       ! total_dens = n_0 + rhostar*<n_1>_FS
       var(:,1) = var(:,1)*rhostar*minor_r*fnorm+spec(n)%dens_prof(li1:li2)
       ! total_temp = T_0 + rhostar*<T_1>_FS
       var(:,2) = var(:,2)*rhostar*minor_r*fnorm+spec(n)%temp_prof(li1:li2)


       !!****************** Fluxes **********************! 

       ! particle flux
       mom(:,:,:,1) = -REAL(CONJG(momc(:,:,:,1,n)) * dfields_dy(:,:,:,1))
       ! heat flux
       mom(:,:,:,2) = -Real(Conjg(momc(:,:,:,2,n)) * dfields_dy(:,:,:,1))
       ! momentum flux
       mom(:,:,:,3) = -real(conjg(momc(:,:,:,3,n)) *  dfields_dy(:,:,:,1))
       !add electromagnetic contribution
       if (n_fields.gt.1) then
          mom(:,:,:,1) = mom(:,:,:,1) + &
               & REAL(CONJG(momc(:,:,:,4,n)) * dfields_dy(:,:,:,2))
          mom(:,:,:,2) = mom(:,:,:,2) + &
               & Real(Conjg(momc(:,:,:,5,n)) * dfields_dy(:,:,:,2))
          mom(:,:,:,3) = mom(:,:,:,3) + &
               & real(conjg(momc(:,:,:,6,n)) * dfields_dy(:,:,:,2))
       endif

       
       ! take into account negative ky, except for 0-mode
       mom = 2. * mom
       if (p_has_0_mode) mom(:,lj1,:,:) = 0.5 * mom(:,lj1,:,:)
       
       !! Fourier space averaging
       do k=lk1,lk2
          do j=lj1,lj2
             ! particle flux
             var(:,5) = var(:,5) + mom(:,j,k,1)*geom%jacobian(li1:li2,pj1,k)
             ! heat flux
             var(:,6) = var(:,6) + mom(:,j,k,2)*geom%jacobian(li1:li2,pj1,k)
             !momentum flux
             var(:,7) = var(:,7) + mom(:,j,k,3)*geom%jacobian(li1:li2,pj1,k)
          enddo
       enddo
       var(li1:li2,5) = var(li1:li2,5)*spec(n)%dens*fnorm
       var(li1:li2,6) = var(li1:li2,6)*spec(n)%dens*spec(n)%temp*fnorm
       var(li1:li2,7) = var(li1:li2,7)*spec(n)%dens*&
                       &sqrt(2.0*spec(n)%temp*spec(n)%mass)*fnorm
       
       do o=5,7
          Call my_sum_to_0_real(var(:,o), li0, mpi_comm_vw)
          Call my_sum_to_0_real(var(:,o), li0, mpi_comm_z)
       enddo
       
       !neoclassical moments: no yx_order allowed so far
       !each processor takes only first ky mode, later select ky=0
       !sum over velocity space
       mom_neo = real(momc_neo(:,lj1,:,:,n))
       Call my_sum_to_0_real(mom_neo, Size(mom_neo), mpi_comm_vw)
       !flux_surface_average (only ky=0 mode contributes):
       neoflux = 0.0
       if (p_has_0_mode) then
          do o=1,4
             do k=lk1,lk2
                neoflux(pi1:pi2,o) = neoflux(pi1:pi2,o) + &
                    &mom_neo(pi1:pi2,k,o)*geom%Jacobian(pi1:pi2,pj1,k)
             enddo
             Call my_sum_to_0_real(neoflux(:,o), Size(neoflux(:,o)), mpi_comm_z)
             neoflux(:,o) = neoflux(:,o)*fnorm
             !call my_real_gather_to_0(neoflux(:,o),pi1+1,pi0,nx0,mpi_comm_x)
             if ((my_pey+my_pez+my_pev+my_pew).eq.0) then 
                Call mpi_gather(neoflux(li1,o), li0, MPI_REAL_TYPE,&
                     neoflux_out(0,o), li0, MPI_REAL_TYPE,&
                     0, mpi_comm_x, ierr)
             endif
          enddo
       endif
       
       if ((my_pey+my_pez+my_pev+my_pew).eq.0) then 

          do o=1,2
             Call mpi_gather(var(li1,o), li0, MPI_REAL_TYPE,&
                  out_arr(0,o), li0, MPI_REAL_TYPE,&
                  0, mpi_comm_x, ierr)
             !careful: LOG only where collected data makes sense..
             if (write_prof_pe) then
                temp_arr = -LOG(out_arr(:,o))
                CALL lag3deriv(temp_arr, xval, nx0, &
                     &out_arr(:,o+2), xval, nx0)
!               call calc_gradient(out_arr(:,o),out_arr(:,o+2))
                out_arr(:,o+2) = out_arr(:,o+2)/(rhostar*minor_r)
             endif
          enddo

          do o=5,7
             Call mpi_gather(var(li1,o), li0, MPI_REAL_TYPE,&
                  out_arr(0,o), li0, MPI_REAL_TYPE,&
                  0, mpi_comm_x, ierr)
          enddo
         
          IF (write_prof_pe) THEN
             
#ifdef VOL_AVG
             write(*,"(A,F14.6)") spec(n)%name,time
             write(*,"(A,ES12.4)") 'G=', SUM(out_arr(:,5))
             write(*,"(A,ES12.4)") 'Q=', SUM(out_arr(:,6))
             write(*,"(A,ES12.4)") 'P=', SUM(out_arr(:,7))
             print*, ''
#endif

             !for now, the output is GNUPLOT compatible
             !in future, we might change to binary or HDF5 output
             WRITE(PROFFILE(n), "(A,F14.6,1X,I5)") '#',time, block_nr
             do i=0,nx0-1
                write(PROFFILE(n), "(13ES16.6)") xval_a(i), xval(i), &
                     &spec(n)%temp*out_arr(i,2),spec(n)%dens*out_arr(i,1), &
                     &out_arr(i,4),out_arr(i,3),out_arr(i,5),out_arr(i,6),&
                     &out_arr(i,7),neoflux_out(i,1),neoflux_out(i,2)&
                     &,neoflux_out(i,3),neoflux_out(i,4)
             enddo
             
             WRITE(PROFFILE(n), "(A)") ''
             WRITE(PROFFILE(n), "(A)") ''

            call flush(PROFFILE(n))
          ENDIF

       endif !my_pey, my_pez, my_pev, my_pew == 0?

       !now distribute the result to all other processes
       call mpi_bcast(out_arr(0,1),7*nx0,MPI_REAL_TYPE,0,mpi_comm_xyzvw,ierr) !works for species parall.??

       !calculate running average
       if ((avgprof_stime.ge.0.0).and.(avgprof_stime.le.time)) then
          if (avgprof(nx0/2,n,0).eq.0.0) then !first avg. step
             avgprof(:,n,0) = xval_a
             if (comp_type.eq.'NC') then
                do o=1,7
                   avgprof(:,n,o) = out_arr(:,o)
                enddo
             endif
             do o=1,4
                new_profs(:,n,o) = out_arr(:,o)
             enddo
             !set avgprof_stime to actually used starting time
             if (n==ln2) avgprof_stime = time 
          else
             do o=1,7
                avgprof(:,n,o) = avgprof(:,n,o)+&
                     & out_arr(:,o)*(time-old_prof_time)
             enddo
             do o=1,2
                new_profs(:,n,o) = avgprof(:,n,o)/(time-avgprof_stime)
                temp_arr = -LOG(new_profs(:,n,o))
                CALL lag3deriv(temp_arr, xval, nx0, &
                     &new_profs(:,n,o+2), xval, nx0)
!                call calc_gradient(new_profs(:,n,o),new_profs(:,n,o+2))
                new_profs(:,n,o+2) = new_profs(:,n,o+2)/(rhostar*minor_r)
             enddo
          endif
       else
          do o=1,4
             new_profs(:,n,o) = out_arr(:,o)
          enddo
       endif

       !--- check for gyroorder consistency ---
       !(i.e., adapt profiles if deviations from initial profile are 
       !larger than reset_limit)
       if (smooth_reset) then
          call smooth_profiles(new_profs(:,n,1:4),tmpprof,smooth_reset_width)
          reset = (reset.or.(ANY(ABS(tmpprof(:,1)/spec(n)%dens_prof-1.0)&
               .GT.reset_limit)).or.(ANY(ABS(tmpprof(:,2)/spec(n)%temp_prof-1.0)&
               .GT.reset_limit)))
        else
          reset = (reset.OR.(ANY(ABS(new_profs(:,n,1)/spec(n)%dens_prof-1.0)&
               .GT.reset_limit)).OR.(ANY(ABS(new_profs(:,n,2)/spec(n)%temp_prof-1.0)&
               .GT.reset_limit)))
       endif

    enddo !n loop
    do o=1,4
       call mpi_allgather(MPI_IN_PLACE,nx0*ln0,MPI_REAL_TYPE,new_profs(0,0,o),nx0*ln0,&
            MPI_REAL_TYPE,mpi_comm_spec,ierr)
    enddo
    
    loc_reset=reset
    call mpi_allreduce(loc_reset,reset,1,MPI_LOGICAL,MPI_LOR,mpi_comm_spec,ierr)

    if (avgprof_stime.ge.0.0) then
       do o=1,7
          call my_real_gather_to_0(avgprof(:,:,o),px0*ln1+1,ln0*px0,&
               &px0*n_spec,mpi_comm_spec)
       enddo
    endif

    old_prof_time = time

    block_nr = block_nr+1

    PERFOFF
    
  End Subroutine diag_profile

  Subroutine finalize_diag_profile
    Integer :: i, n, o, ierr
    logical :: dummy

    if (old_prof_time.ne.time) call diag_profile(dummy)

    Deallocate(prof_mat, new_profs)

    if ((nonlinear).and.(avgprof_stime.ge.0.0)) then        
       do o = 1,7
          avgprof(:,:,o)=avgprof(:,:,o)/(time-avgprof_stime)
       enddo

       if (write_prof_pe) then
          do n=ln1,ln2
             WRITE(PROFFILE(n), "(A,I3,A,F14.6,A,F14.6)") '#time averaged profiles, block nr ', block_nr, &
                  &', time window: ', avgprof_stime, ' - ', time
             do i=0,nx0-1
                write(PROFFILE(n), "(8ES16.6)") xval_a(i), xval(i), &
                     &spec(n)%temp*avgprof(i,n,2),spec(n)%dens*avgprof(i,n,1), &
                     &avgprof(i,n,4),avgprof(i,n,3),avgprof(i,n,5:7)
!                write(PROFFILE(n), "(6ES12.4)") avgprof(i,n,1:6)
             enddo
          enddo
       endif
    endif

    if (par_in_dir.ne.'skip_parfile') then
       DEALLOCATE(avgprof)
    else
       call mpi_bcast(avgprof,SIZE(avgprof),MPI_REAL_TYPE,&
              &0,MY_MPI_COMM_WORLD,ierr)
    endif


    If (write_prof_pe) THEN
       do n=ln1,ln2
          CLOSE(PROFFILE(n))
       enddo
    endif

    DEALLOCATE(PROFFILE)

  End Subroutine finalize_diag_profile
  
!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of this diag
  Real Function mem_est_diag_srcmom(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    
    mem_loc = px0*ln0*9*SIZE_OF_REAL_MB

    mem_est_diag_srcmom=mem_req_in+mem_loc
  End Function mem_est_diag_srcmom

!*****************************************************************

  Subroutine initialize_diag_srcmom
    Implicit None

    Integer :: n

    ALLOCATE(SRCMOMFILE(ln1:ln2))

    write_srcmom_pe = (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec))

    if (write_srcmom_pe) then
       ! Moment files for each species in a different file
       DO n=ln1,ln2
          call get_unit_nr(SRCMOMFILE(n))
          OPEN(SRCMOMFILE(n), file=trim(diagdir)//'/srcmom_'//&
               &trim(spec(n)%name)//''//trim(file_extension),&
               form='unformatted', status=filestat, position=filepos)  
       END DO
    endif

    write(srcmom_format,'(A,I3,A)') '(',px0,'ES14.6)'

  End Subroutine initialize_diag_srcmom

  Subroutine diag_srcmom(time)
    real, intent(in) :: time
    integer :: i,k,l,m,n,o,s,ierr
    complex,dimension(li1:li2) :: tmp
    real,dimension(li1:li2,ln1:ln2,3,3) :: src_mom,src_mom_tmp
    real,dimension(0:nx0-1,3) :: out_arr
    
    PERFON('src_moms')
    
    src_mom = 0.0
    
    if (p_has_0_mode) then
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do s=1,3
                      tmp = 0.0
                      if ((s==1).and.heatsrc_on) &
                           & call add_ck_heat(l,m,n,tmp)
                      if ((s==2).and.partsrc_on.and.psource_type==1) &
                           & call add_ck_part_1(l,m,n,tmp)
                      if ((s==2).and.partsrc_on.and.psource_type==2) &
                           & call add_ck_part_2(l,m,n,tmp)
                      if ((s==3).and.(f0_heat_src_active.or.include_f0_contr)) &
                           & call add_f0_term_df_1d(&
                           & time,k,l,m,n,tmp)
                      tmp = tmp*geom%jacobian(li1:li2,lj1,k)
                      src_mom(:,n,s,1) = src_mom(:,n,s,1)+&
                           &tmp*mat_00(li1:li2,lj1,k,l,m)
                      src_mom(:,n,s,2) = src_mom(:,n,s,2)+&
                           &tmp*mat_10(li1:li2,lj1,k,l,m,n)
                      src_mom(:,n,s,3) = src_mom(:,n,s,3)+&
                           &tmp*(mu(m)*geom%Bfield(li1:li2,lj1,k)+vp(l)**2)*&
                           &mat_00(li1:li2,lj1,k,l,m)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    
    !sum v,w,z
    Call mpi_allreduce(src_mom,src_mom_tmp,li0*ln0*3*3,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_zvw, ierr)
    
    do i=li1,li2
       src_mom(i,:,:,:) = src_mom_tmp(i,:,:,:)/real(nz0*geom%avg_jaco_yz(i))
    enddo

    !now collect the radial data on the 0'th process of each species block
    do n=ln1,ln2
!       if (write_srcmom_pe) write(SRCMOMFILE(n), "(F13.6)") time 
       if (write_srcmom_pe) write(SRCMOMFILE(n)) time
       do s=1,3
          out_arr(li1:li2,1:3) = src_mom(li1:li2,n,s,1:3)
          do o=1,3
             call my_real_gather_to_0(out_arr(:,o),li1+1,li0,&
                  &px0,mpi_comm_x)
          enddo
          if (write_srcmom_pe) then
             !& write(SRCMOMFILE(n),srcmom_format) src_mom(:,n,s,o)
             write(SRCMOMFILE(n)) out_arr
          endif
       enddo
       if (write_srcmom_pe) call flush(SRCMOMFILE(n))
    enddo
             
    PERFOFF   

  End Subroutine diag_srcmom


  Subroutine finalize_diag_srcmom
    integer :: n

    if (write_srcmom_pe) then
       do n=ln1,ln2
          CLOSE(SRCMOMFILE(n))
       enddo
    endif
    
    DEALLOCATE(SRCMOMFILE)

  End Subroutine finalize_diag_srcmom

End Module diagnostics_df
