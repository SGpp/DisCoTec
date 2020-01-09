#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
!> This modules contains routines for specialized diagnostics 
Module diagnostics_extended
  Use par_mod 
  Use file_io, only: get_unit_nr
  Use vel_space, only: fm, mat_00
  use calc_rhs, only: calc_rhs_only,ptr_bar_emfields, ptr_barchi, &
       ptr_dbarchidxy, ptr_dgdxy, rhs_nl, this_nonlinear_term
  use geometry
  use communications
  use aux_fields, only: calc_aux_fields
  Use gyro_average_ff_mod, only: gyro_average_ff
  Use gyro_average_df_mod, only: gyro_average_df
  use collisions, only: equ_collisions
  use dfdzv_terms !, only: add_dfdz, add_dfdv, add_hypv, add_hypz
  use dgdxy_terms
  use dzv_terms
  use dfdxy_terms, only: add_dfdxy
  use dchidz_term, only: add_dchidz
  use dchidxy_terms, only: add_dchidxy, add_dchidxy_orig
  use sources_mod, only: add_f1_sources_block
  use blockindex
  use axpy
  use prefactors
  use numerical_damping
  USE boundaries, only: pb_xshift
  use diagnostics_energy, only: get_cfgamma_h, get_cfgamma_g, get_cfgamma_f,&
      get_energy_terms,scalar_product, block_mult
#ifdef WITHFUTILS
  use futils
  use hashtable
#endif
  Implicit None

  PUBLIC :: initialize_all_diags_nlt, exec_all_diags_nlt, &
       &finalize_all_diags_nlt 

  PRIVATE 

  INTEGER:: nlt_info_handle,nlt_out,nlt_modes=60
  complex, dimension(:,:,:,:,:,:,:), allocatable ::  nlt_pod_mode,nlt_pod_mode_wcons
  integer, dimension(:), allocatable :: nlt_pod_out,nlt_info_pod
  real :: last_energy,last_energy_nlt(60)
  real, dimension(:,:), allocatable :: last_energy_nlt_pod
  real, dimension(:), allocatable :: time_storage

Contains

!!!******************************************************************!!!

  SUBROUTINE initialize_all_diags_nlt

    if (istep_nlt.gt.0) then
#ifdef with_extended_diags
       call initialize_diag_nlt 
       if(nlt_pod) call initialize_diag_nlt_pod
#else
       if (mype==0) Write(*,'(a)') 'ERROR: with_extended_diags must be set in switches.h for using diag_nlev'
       stop
#endif
    end if
    last_energy_nlt=0.0

  END SUBROUTINE initialize_all_diags_nlt

!!!******************************************************************!!!

  SUBROUTINE exec_all_diags_nlt(itime,istep_nlt)
    integer, intent(in) :: itime, istep_nlt

#ifdef with_extended_diags
    IF(istep_nlt.GT.0) THEN
       IF (MOD(itime,istep_nlt).eq.0) then
         CALL diag_nlt
         if(nlt_pod) call diag_nlt_pod
       END IF
    END IF
#endif

  END SUBROUTINE exec_all_diags_nlt

!!!******************************************************************!!!

  SUBROUTINE finalize_all_diags_nlt

#ifdef with_extended_diags
    if (istep_nlt.gt.0) then
       call finalize_diag_nlt
       if(nlt_pod) call finalize_diag_nlt_pod
    end if
#endif

  END SUBROUTINE finalize_all_diags_nlt
!!!******************************************************************!!!

#ifdef with_extended_diags
  subroutine initialize_diag_nlt
    
    integer :: i
    

      do i=1,num_nlt_modes 
        if(mype==0) write(*,*) "kx_index",kx_nlt_ind(i),"ky_index",ky_nlt_ind(i)
        if((kx_nlt_ind(i).gt.nx0-1).or.(ky_nlt_ind(i).gt.nky0-1)) stop "Selected wavenumbers out of range for diag_nlt!" 
      end do


      if(mype==0) then
        call get_unit_nr(nlt_out)
        open(unit=nlt_out,file=trim(diagdir)//'/nlt.dat',form='unformatted',status='unknown')
        call get_unit_nr(nlt_info_handle)
        open(unit=nlt_info_handle,file=trim(diagdir)//'/nlt_info.dat',status='unknown')
!        write(nlt_info_handle,*) "#Nonlinear transfer functions analyzed for the following wavenumbers:"
!        write(nlt_info_handle,*) "#kx index, ky index, kx, ky" 
!        do i=1,num_nlt_modes
!          write(nlt_info_handle,*) "#",kx_nlt_ind(i),ky_nlt_ind(i),kx(kx_nlt_ind(i)),ky(ky_nlt_ind(i))
!        end do
        write(nlt_info_handle,*) "#time,nlt_tot_k1,rhs_lin_k1,energy_k1,dedt_k1, . . . k2, . . . k3, . . ."
      end if

  end subroutine initialize_diag_nlt

  subroutine finalize_diag_nlt
  
   if(mype==0) then
    close(nlt_info_handle)
    close(nlt_out)
   end if

  end subroutine finalize_diag_nlt

  subroutine initialize_diag_nlt_pod
    
    integer :: i,ierr,pod_io,dummy_i,j
    character(len=4) :: kx_ind,ky_ind
    integer :: max_cons(1)
    complex, dimension(:,:,:,:,:,:), allocatable :: nlt_pod_in 
    complex :: g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: norm
    !real :: e_tot(11)
    !complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms

      allocate(last_energy_nlt_pod(20,num_nlt_pod_modes+1))
      last_energy_nlt_pod=0.0
      allocate(nlt_pod_out(num_nlt_modes))   !These are the i/o handles for the nlt.dat files
      allocate(nlt_info_pod(num_nlt_modes))  !These are the i/o handles for the nlt_info.dat files

      do i=1,num_nlt_modes 
        if((kx_nlt_ind(i).gt.nx0-1).or.(ky_nlt_ind(i).gt.nky0-1)) stop "Selected wavenumbers out of range for diag_nlt!" 
      end do

      !calculate the maximum number of connections needed
      max_cons=maxval(nx0_nlt_pod)
      !nlt_pod_mode: This is the pod mode at only kx=kx_center
      allocate(nlt_pod_mode(1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,num_nlt_modes,num_nlt_pod_modes))   
      !nlt_pod_mode_wcons:  This is the pod mode with all connections   
      allocate(nlt_pod_mode_wcons(0:max_cons(1)-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,num_nlt_modes,num_nlt_pod_modes))
      nlt_pod_mode=cmplx(0.0,0.0)
      nlt_pod_mode_wcons=cmplx(0.0,0.0)

      do i=1,num_nlt_modes  !Cycle through wavenumbers 
        if(mype==0) then
          write(kx_ind,'(i4.4)') kx_nlt_ind(i)
          write(ky_ind,'(i4.4)') ky_nlt_ind(i)
          call get_unit_nr(nlt_pod_out(i))
          open(unit=nlt_pod_out(i),file=trim(diagdir)//'/nlp_ky'//ky_ind//'kx'//kx_ind//'.dat',form='unformatted',status='unknown')
          call get_unit_nr(nlt_info_pod(i))
          open(unit=nlt_info_pod(i),file=trim(diagdir)//'/nlp_info_ky'//ky_ind//'kx'//kx_ind//'.dat',status='unknown')
          write(nlt_info_pod(i),*) "#time,nlt_tot_pod1,nlt_lin_pod1,energy_pod1,dedt_pod1,...pod2,...pod3"
         end if
  
          allocate(nlt_pod_in(nx0_nlt_pod(i),1,0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1))
          call get_unit_nr(pod_io)
          !Open files:  POD modes from earlier analysis
          if(mype==0) open(unit=pod_io,file=trim(SVD_df_file_path)//'/SVD_df_ky'//ky_ind//'kx'//kx_ind,&
                            form='unformatted',status='unknown')
          !Cycle through pod modes - take the first num_nlt_pod_modes in the file
          do j=1,num_nlt_pod_modes
            if(mype==0) then
              read(pod_io) dummy_i
              write(*,*) "Reading mode", dummy_i,j
              read(pod_io) nlt_pod_in
            end if
            !Send the mode to all processors
            call MPI_BCAST(nlt_pod_in,nx0_nlt_pod(i)*nz0*nv0*nw0*n_spec,MPI_COMPLEX_TYPE,0,my_mpi_comm_world,ierr)
           !!!!!!!!!!!!!!!!Tests
           ! all_sum=sum(sum(sum(sum(sum(sum(nlt_pod_in,1),1),1),1),1),1)
           ! call MPI_ALLREDUCE(all_sum,all_sum2,1,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
           ! write(*,*) "(in)mype,sum,i",mype,sqrt(real(conjg(all_sum)*all_sum)),i
           ! if(mype==0) write(*,*) "sum of all procs (wcons):",all_sum2,i
           ! all_sum=sum(sum(sum(sum(sum(sum(nlt_pod_in(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),1),1),1),1),1),1)
           ! write(*,*) "mype,loc_sum",sqrt(real(conjg(all_sum)*all_sum))
           !!!!!!!!!!!!!!!Tests
            !Now give each local array the right data
            nlt_pod_mode_wcons(0:nx0_nlt_pod(i)-1,:,:,:,:,i,j)=nlt_pod_in(:,1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)          
            nlt_pod_mode(1,:,:,:,:,i,j)=nlt_pod_in(1,1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)          

            !Normalize modes
            call fill_connections(nlt_pod_mode_wcons(0:nx0_nlt_pod(i)-1,:,:,:,:,i,j),&
                 &g_temp,kx_nlt_ind(i),ky_nlt_ind(i),nx0_nlt_pod(i))
            !call spro_of_choice(1,g_temp,g_temp,li1,li2,lj1,lj2,norm)
            call spro_of_choice(2,g_temp,g_temp,li1,li2,lj1,lj2,norm)
            nlt_pod_mode_wcons(:,:,:,:,:,i,j)=nlt_pod_mode_wcons(:,:,:,:,:,i,j)/sqrt(norm)
            nlt_pod_mode(:,:,:,:,:,i,j)=nlt_pod_mode(:,:,:,:,:,i,j)/sqrt(norm)
           !!!!!!!!!!!!!!!!Tests
           ! all_sum=sum(sum(sum(sum(sum(nlt_pod_mode_wcons(:,:,:,:,:,i,j),1),1),1),1),1)
           ! call MPI_ALLREDUCE(all_sum,all_sum2,1,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
           ! write(*,*) "(wcon)mype,sum,i",mype,sqrt(real(conjg(all_sum)*all_sum)),i
           ! if(mype==0) write(*,*) "sum of all procs (wcons):",all_sum2,i
           !!!!!!!!!!!!!!!!Tests
            
          end do !num_nlt_modes
          deallocate(nlt_pod_in)
  
          if(mype==0) close(pod_io)
  
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tests
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tests
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tests
      !open(unit=888,file=trim(diagdir)//'/nlp_tests.dat',status='unknown')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tests
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tests
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tests


  end subroutine initialize_diag_nlt_pod


  subroutine finalize_diag_nlt_pod
  
    integer :: i

    do i=1,num_nlt_modes
      if(mype==0) close(nlt_info_pod(i))
      if(mype==0) close(nlt_pod_out(i))
    end do
    deallocate(nlt_pod_mode)
    deallocate(nlt_pod_mode_wcons)
    deallocate(nlt_pod_out)
    deallocate(nlt_info_pod)
    deallocate(last_energy_nlt_pod)

  end subroutine finalize_diag_nlt_pod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


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

  !>This subroutine performs the mat_00 integral and sum over species in order to reduce from 6D to 3D
  !!(strip mining version)
  subroutine add_energy_integral_block(v6d,v3d,lb1,lb2)
    implicit none
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(in) :: v6d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(inout) :: v3d
    integer :: j,klmn

    if (xy_local) then
       do klmn=lb1,lb2
          call axpy_ij(lij0,mat_00(pi1,pj1,sk(klmn),sl(klmn),sm(klmn)),&
               &v6d(:,:,klmn),v3d(:,:,sk(klmn)))
       enddo
    else
       do klmn=lb1,lb2
          do j=lj1,lj2
             v3d(:,j,sk(klmn)) = v3d(:,j,sk(klmn)) + &
                  &v6d(:,j,klmn)*mat_00(li1:li2,pj1,sk(klmn),sl(klmn),sm(klmn))
          enddo
       enddo
    endif

!    call my_complex_sum_vwspec(v3d,size(v3d)) 

  end subroutine add_energy_integral_block

  !>Sums / integrates over kx (x),ky,z
  subroutine sum_int_3d_comp(v3d,v1d)
    !This subroutine performs the z integral and sum over kx,ky

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(in) :: v3d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d_tmp

    complex, intent(out) :: v1d
    complex :: v1d_loc
    integer :: ierr
      
    v3d_tmp = v3d
    if ((xy_local).and.(evenx.eq.1)) v3d_tmp(hkx+1,:,:) = 0.0
    if (p_has_0_mode) v3d_tmp(:,lj1,:) = 0.5*v3d_tmp(:,lj1,:)

    if (xy_local) then
       v1d_loc = 2.0*sum(sum(sum(v3d_tmp,1),1)*geom%jacobian(pi1,pj1,:))/&
            (real(nz0)*geom%avg_jaco)
    else
       v1d_loc = 2.0*sum(sum(v3d_tmp,2)*geom%jacobian(:,pj1,:))/&
            (real(nx0*nz0)*geom%avg_jaco)
    endif

    call mpi_allreduce(v1d_loc,v1d,1,MPI_COMPLEX_TYPE,MPI_SUM, mpi_comm_xyz, ierr)

  end subroutine sum_int_3d_comp


  !>This subroutine calculates nonlinear transfer functions
  subroutine diag_nlt

    implicit none
    integer :: i,imode,ierr
!    complex :: nlt_rhs(li1:li2,-1*(nky0-1):nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) !The full nonlinear term
    real :: nlt_rhs_2d(li1:li2,-1*(nky0-1):nky0-1)                                 !Integrated over z,v,mu
    real :: nlt_tot(nlt_modes)
    complex :: cfgamma_g(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: cfgamma_chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_xy_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_xy_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: g_1_xy1(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: g_1_xy2(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: chi_xy1(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: chi_xy2(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    !real :: kxp,kyp,kxpp,kypp !kxp=kx_prime, kxpp = kx_double_prime = kx-kxp, etc.
    !integer :: ip,jp,ipp,jpp !Indices corresponding to above  
    !real :: ckkp1,ckkp2              !coupling coefficient
    real :: kx0,ky0           !kx and ky values for which to calculate coupling coefficient
    integer :: i0,j0          !Indices corresponding to above
    !integer :: ip0,jp0        !Loop indices for kp
    integer :: k,l,m,n
    integer :: xfersize
    !logical :: conjg_p,conjg_pp
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields) :: fields_1
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2) :: f_1    
    complex, allocatable, dimension(:,:,:,:,:,:) :: apar_bar,h_1
    complex, allocatable, dimension(:,:,:,:,:,:) :: cfgamma_globy_chi,chi_globy,g_1_globy,dummy_globy !Global in ky
    complex, allocatable, dimension(:,:,:,:,:,:) :: cfgamma_globy_g
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms
    real :: dE_z(lk1:lk2)
    real :: e_tot(11),rhs_temp,r_temp,rhs_tot(nlt_modes)
    real :: rhs_tot_temp(nlt_modes),en_k(nlt_modes),dendt_k(nlt_modes) 

    if(mype==0) write(nlt_out) time
  
    allocate(h_1(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    !allocate(temp(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    if(n_fields.gt.1) allocate(apar_bar(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))

    if(n_procs_y.gt.1) then
      allocate(cfgamma_globy_chi(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
      cfgamma_globy_chi=cmplx(0.0,0.0)
      allocate(cfgamma_globy_g(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
      cfgamma_globy_g=cmplx(0.0,0.0)
      allocate(chi_globy(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
      chi_globy=cmplx(0.0,0.0)
      allocate(g_1_globy(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
      g_1_globy=cmplx(0.0,0.0)
      allocate(dummy_globy(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
      dummy_globy=cmplx(0.0,0.0)
    end if

!Get cfgamma (same as in get_energy_terms)
    if (arakawa_zv) then 
       call calc_aux_fields(g_1,fields_1,f_1,.true.,h_1) 
       call get_cfgamma_h(h_1,cfgamma_chi)
    else
       call calc_aux_fields(g_1,fields_1,f_1,.false.) 
       call get_cfgamma_f(f_1,fields_1,cfgamma_chi)
    endif

    !now get cfgamma_g - i.e. only g part
    call get_cfgamma_g(g_1,cfgamma_g)
    !now get cfgamma_chi - i.e. total - cfgamma_g
    cfgamma_chi=cfgamma_chi-cfgamma_g


    !construct chi
    do k=lk1,lk2
      do l=ll1,ll2
        do m=lm1,lm2
          do n=ln1,ln2
            call gyro_average_ff(fields_1(:,:,k,1),chi(:,:,k,l,m,n),k,m,n)
            if(n_fields.gt.1) call gyro_average_ff(fields_1(:,:,k,2),apar_bar(:,:,k,l,m,n),k,m,n)
          end do
        end do
      end do
    end do    

   if(n_fields.gt.1) then
     do l=ll1,ll2
      do n=ln1,ln2
          chi(:,:,lk1:lk2,l,lm1:lm2,n)= chi(:,:,lk1:lk2,l,lm1:lm2,n)-&
              (2.0*spec(n)%temp/spec(n)%mass)**0.5*vp(l)*apar_bar(:,:,lk1:lk2,l,lm1:lm2,n)
      end do
     end do
   end if

   if(allocated(apar_bar)) deallocate(apar_bar)

   !get y modes on all processors (this can be improved)
   if(n_procs_y.gt.1) then
    xfersize=lg0*nky0*lk0*ll0*lm0*ln0
    g_1_globy(:,lj1:lj2,:,:,:,:)=g_1(:,:,:,:,:,:)
    cfgamma_globy_chi(:,lj1:lj2,:,:,:,:)=cfgamma_chi(:,:,:,:,:,:)
    cfgamma_globy_g(:,lj1:lj2,:,:,:,:)=cfgamma_g(:,:,:,:,:,:)
    chi_globy(:,lj1:lj2,:,:,:,:)=chi(:,:,:,:,:,:)
    call MPI_ALLREDUCE(g_1_globy,dummy_globy,xfersize,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_y,ierr)
    g_1_globy=dummy_globy
    call MPI_ALLREDUCE(chi_globy,dummy_globy,xfersize,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_y,ierr)
    chi_globy=dummy_globy
    call MPI_ALLREDUCE(cfgamma_globy_chi,dummy_globy,xfersize,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_y,ierr)
    cfgamma_globy_chi=dummy_globy
    call MPI_ALLREDUCE(cfgamma_globy_g,dummy_globy,xfersize,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_y,ierr)
    cfgamma_globy_g=dummy_globy
    deallocate(dummy_globy)
   end if

   !loop over the modes of interest (input in parameters file, max=10 at this time).

   !For getting entire RHS later
   call get_energy_terms(g_1,e_tot,dE_terms)

   do i=1,num_nlt_modes   

     kx0=kx(kx_nlt_ind(i))
     ky0=ky(ky_nlt_ind(i))
     i0=kx_nlt_ind(i)
     j0=ky_nlt_ind(i)
     !nlt_rhs=cmplx(0.0,0.0)

     if(n_procs_y.eq.1) then
       call construct_nlt(i,i0,j0,kx0,ky0,cfgamma_g,cfgamma_chi,chi,nlt_rhs_2d)
     else
       call construct_nlt_ypar(i,i0,j0,kx0,ky0,g_1_globy,cfgamma_globy_g,cfgamma_globy_chi,chi_globy,nlt_rhs_2d)
     end if
     
     nlt_tot(i)=sum(sum(nlt_rhs_2d,1),1)

     !output
     if(mype==0) write(nlt_out) nlt_rhs_2d

     !This outputs the entire linear RHS for energy equation for kx0,ky0.  
     !This should balance (on time-average in saturated state) with the sum 
     !over all k of the nonlinear transfer function (nlt_tot)
      rhs_tot(i)=0.0
      en_k(i)=0.0
      dendt_k(i)=0.0
      if(j0.ge.lj1.and.j0.le.lj2) then
         dE_z(:)=real(dE_terms(i0,j0,:,2)-dE_terms(i0,j0,:,9))
         rhs_temp=sum(dE_z*geom%jacobian(pi1,pj1,:)/(real(nz0)*geom%avg_jaco),1)  !total linear  RHS
         call mpi_allreduce(rhs_temp,r_temp,1,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)
         rhs_tot(i)=r_temp
         rhs_temp=sum(real(dE_terms(i0,j0,:,1))*geom%jacobian(pi1,pj1,:)/(real(nz0)*geom%avg_jaco),1)  !energy
         call mpi_allreduce(rhs_temp,r_temp,1,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)
         en_k(i)=r_temp
         dendt_k(i)=(en_k(i)-last_energy_nlt(i))/real(istep_nlt*dt)
      end if

   end do !num_nlt_modes

   if(n_procs_y.gt.1) then
     call mpi_allreduce(rhs_tot,rhs_tot_temp,nlt_modes,MPI_REAL_TYPE,MPI_SUM,&
          &mpi_comm_y,ierr)  
     rhs_tot=rhs_tot_temp
     call mpi_allreduce(en_k,rhs_tot_temp,nlt_modes,MPI_REAL_TYPE,MPI_SUM,&
          &mpi_comm_y,ierr)  
     en_k=rhs_tot_temp
     call mpi_allreduce(dendt_k,rhs_tot_temp,nlt_modes,MPI_REAL_TYPE,MPI_SUM,&
          &mpi_comm_y,ierr)  
     dendt_k=rhs_tot_temp
   end if
   last_energy_nlt=en_k

   !Output time traces of 1) total nonlinear, 2) total linear, 3) total energy, 4) time derivative of energy
   if(mype==0) then
      write(nlt_info_handle,'(1es16.6)', advance='no') time
      do imode=1,nlt_modes-1
         write(nlt_info_handle,'(4es16.6)', advance='no') &
              &nlt_tot(imode),rhs_tot(imode),en_k(imode),dendt_k(imode)
      enddo
      write(nlt_info_handle,'(4es16.6)') &
           &nlt_tot(nlt_modes),rhs_tot(nlt_modes),en_k(nlt_modes),dendt_k(nlt_modes)      
      If (mype == 0) call flush(nlt_info_handle)
   endif

   If (mype == 0) call flush(nlt_out)

   if(n_procs_y.gt.1) then
      deallocate(cfgamma_globy_chi)
      deallocate(cfgamma_globy_g)
      deallocate(chi_globy)
      deallocate(g_1_globy)
    end if

  end subroutine diag_nlt

  subroutine construct_nlt(i,i0,j0,kx0,ky0,cfgamma_g,cfgamma_chi,chi,nlt_rhs_2d)

    complex, intent(in) :: cfgamma_g(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex, intent(in) :: cfgamma_chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex, intent(in) :: chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    integer, intent(in) :: i0,j0          !Indices corresponding to above
    real, intent(out) :: nlt_rhs_2d(li1:li2,-1*(nky0-1):nky0-1)                                 !Integrated over z,v,mu
    complex :: nlt_rhs(li1:li2,-1*(nky0-1):nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) !The full nonlinear term
    real,intent(in) :: kx0,ky0           !kx and ky values for which to calculate coupling coefficient
    integer,intent(in) :: i
    complex :: cfgamma_xy_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: cfgamma_xy_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: g_1_xy1(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: g_1_xy2(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: chi_xy1(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: chi_xy2(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    real :: kxp,kyp,kxpp,kypp !kxp=kx_prime, kxpp = kx_double_prime = kx-kxp, etc.
    integer :: ip,jp,ipp,jpp !Indices corresponding to above  
    real :: ckkp1,ckkp2              !coupling coefficient
    integer :: ip0,jp0        !Loop indices for kp
    logical :: conjg_p,conjg_pp


    nlt_rhs=cmplx(0.0,0.0)
     !loop over k-prime
     do ip0=0,nx0-1
       do jp0=-(nky0-1),nky0-1
         if((kx_nlt_ind(i).ne.ip0.or.ky_nlt_ind(i).ne.jp0).and.&   !don't calculate for kx0,ky0
            (ip0.ne.0.or.jp0.ne.0)) then       !don't calculate for (0,0)

           kxp=kx(ip0)   !kx_prime
           kyp=kymin*jp0 !ky_prime
           kxpp=kx0-kxp  !kx-kx_prime
           kypp=ky0-kyp  !ky-ky_prime
           ckkp1=kxp*ky0-kx0*kyp    !coupling coefficient for T_chi
           ckkp2=kxpp*ky0-kx0*kypp    !coupling coefficient for T_g

           if((kxpp.le.maxval(kx).and.kxpp.ge.minval(kx)).and.&
              (abs(kypp).le.maxval(ky))) then    !Make sure that kpp is in the box
             call get_k_indices(kxp,kyp,ip,jp,conjg_p)        !get indices corresponding to kp
                                                              !must take complex conjugate for ky>0 modes (indicated by conjg_p)
             call get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp)   !get indices corresponding to kpp

               cfgamma_xy_g(:,:,:,:)=cfgamma_g(i0,j0,:,:,:,:)
               cfgamma_xy_chi(:,:,:,:)=cfgamma_chi(i0,j0,:,:,:,:)

               if(conjg_p) then
                 chi_xy1(:,:,:,:)=conjg(chi(ip,jp,:,:,:,:))
                 g_1_xy2(:,:,:,:)=conjg(g_1(ip,jp,:,:,:,:))
               else
                 chi_xy1(:,:,:,:)=chi(ip,jp,:,:,:,:)
                 g_1_xy2(:,:,:,:)=g_1(ip,jp,:,:,:,:)
               end if

               if(conjg_pp) then
                 g_1_xy1(:,:,:,:)=conjg(g_1(ipp,jpp,:,:,:,:))
                 chi_xy2(:,:,:,:)=conjg(chi(ipp,jpp,:,:,:,:))
               else
                 g_1_xy1(:,:,:,:)=g_1(ipp,jpp,:,:,:,:)
                 chi_xy2(:,:,:,:)=chi(ipp,jpp,:,:,:,:)
               end if

             !Note all complex conjugations have already been made
             nlt_rhs(ip0,jp0,:,:,:,:)=ckkp1*cfgamma_xy_chi*chi_xy1*g_1_xy1+ckkp2*cfgamma_xy_g*chi_xy2*g_1_xy2
             if(nlt_symmetrize) then
                 nlt_rhs(ip0,jp0,:,:,:,:)=nlt_rhs(ip0,jp0,:,:,:,:)-&
                        ckkp1*cfgamma_xy_chi*chi_xy2*g_1_xy2-ckkp2*cfgamma_xy_g*chi_xy1*g_1_xy1
             end if
           end if !max kx ky check

         end if !if(kx_nlt_ind(i).ge.ip.or.ky_nlt_ind(j).ne.jp) then
       end do !ip loop
     end do !jp loop

     call nlt_integral(nlt_rhs,nlt_rhs_2d)

  end subroutine construct_nlt


  subroutine construct_nlt_ypar(i,i0,j0,kx0,ky0,g_1_globy,cfgamma_globy_g,cfgamma_globy_chi,chi_globy,nlt_rhs_2d)

    complex, intent(in) :: cfgamma_globy_g(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex, intent(in) :: cfgamma_globy_chi(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex, intent(in) :: chi_globy(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex, intent(in) :: g_1_globy(li1:li2,0:nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    integer, intent(in) :: i0,j0          !Indices corresponding to above
    real, intent(out) :: nlt_rhs_2d(li1:li2,-1*(nky0-1):nky0-1)                                 !Integrated over z,v,mu
    integer, intent(in) :: i
    complex :: nlt_rhs(li1:li2,-1*(nky0-1):nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) !The full nonlinear term
    real, intent(in) :: kx0,ky0           !kx and ky values for which to calculate coupling coefficient
    complex :: cfgamma_xy_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: cfgamma_xy_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: g_1_xy1(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: g_1_xy2(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: chi_xy1(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: chi_xy2(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
!    complex :: cfgamma_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    real :: kxp,kyp,kxpp,kypp !kxp=kx_prime, kxpp = kx_double_prime = kx-kxp, etc.
    integer :: ip,jp,ipp,jpp !Indices corresponding to above  
    real :: ckkp1,ckkp2              !coupling coefficient
    integer :: ip0,jp0        !Loop indices for kp
    logical :: conjg_p, conjg_pp

    nlt_rhs=cmplx(0.0,0.0)

     !loop over k-prime
     do ip0=0,nx0-1
       do jp0=-(nky0-1),nky0-1
         if((kx_nlt_ind(i).ne.ip0.or.ky_nlt_ind(i).ne.jp0).and.&   !don't calculate for kx0,ky0
            (ip0.ne.0.or.jp0.ne.0)) then       !don't calculate for (0,0)

           kxp=kx(ip0)   !kx_prime
           kyp=kymin*jp0 !ky_prime
           kxpp=kx0-kxp  !kx-kx_prime
           kypp=ky0-kyp  !ky-ky_prime
           ckkp1=kxp*ky0-kx0*kyp    !coupling coefficient for T_chi
           ckkp2=kxpp*ky0-kx0*kypp    !coupling coefficient for T_g

           if((kxpp.le.maxval(kx).and.kxpp.ge.minval(kx)).and.&
              (abs(kypp).le.maxval(ky))) then    !Make sure that kpp is in the box
             call get_k_indices(kxp,kyp,ip,jp,conjg_p)        !get indices corresponding to kp
                                                              !must take complex conjugate for ky>0 modes (indicated by conjg_p)
             call get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp)   !get indices corresponding to kpp
               cfgamma_xy_chi(:,:,:,:)=cfgamma_globy_chi(i0,j0,:,:,:,:)
               cfgamma_xy_g(:,:,:,:)=cfgamma_globy_g(i0,j0,:,:,:,:)

               if(conjg_p) then
                 chi_xy1(:,:,:,:)=conjg(chi_globy(ip,jp,:,:,:,:))
                 g_1_xy2(:,:,:,:)=conjg(g_1_globy(ip,jp,:,:,:,:))
               else
                 chi_xy1(:,:,:,:)=chi_globy(ip,jp,:,:,:,:)
                 g_1_xy2(:,:,:,:)=g_1_globy(ip,jp,:,:,:,:)
               end if

               if(conjg_pp) then
                 g_1_xy1(:,:,:,:)=conjg(g_1_globy(ipp,jpp,:,:,:,:))
                 chi_xy2(:,:,:,:)=conjg(chi_globy(ipp,jpp,:,:,:,:))
               else
                 g_1_xy1(:,:,:,:)=g_1_globy(ipp,jpp,:,:,:,:)
                 chi_xy2(:,:,:,:)=chi_globy(ipp,jpp,:,:,:,:)
               end if

             !Note all complex conjugations have already been made
             nlt_rhs(ip0,jp0,:,:,:,:)=ckkp1*cfgamma_xy_chi*chi_xy1*g_1_xy1+ckkp2*cfgamma_xy_g*chi_xy2*g_1_xy2
             if(nlt_symmetrize) then
                 nlt_rhs(ip0,jp0,:,:,:,:)=nlt_rhs(ip0,jp0,:,:,:,:)-&
                        ckkp1*cfgamma_xy_chi*chi_xy2*g_1_xy2-ckkp2*cfgamma_xy_g*chi_xy1*g_1_xy1
             end if

           end if !max kx ky check

         end if !if(kx_nlt_ind(i).ge.ip.or.ky_nlt_ind(j).ne.jp) then
       end do !ip loop
     end do !jp loop

     call nlt_integral(nlt_rhs,nlt_rhs_2d)

  end subroutine construct_nlt_ypar


  !>This subroutine gets the correct indices for a given kx and ky (including negative).  If ky is negative, take_conjg=T
  subroutine get_k_indices(kx_in,ky_in,i,j,take_conjg)
    implicit none

    real, intent(in) :: kx_in,ky_in
    integer, intent(out) :: i,j
    logical, intent(out) :: take_conjg
    real :: kx0,ky0
             !call get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp) 

    take_conjg=.false.
    kx0=kx_in
    ky0=ky_in 

    j=nint(ky0/kymin)
    if(ky0.lt.0.0) then
      take_conjg=.true.
      j=-1*j
      kx0=-1.0*kx0
    end if

    if(kx0.ge.0.0) then
      i=nint(kx0/kxmin)
    else
      i=nint(kx0/kxmin+nx0)
    end if

  end subroutine get_k_indices

  !>This subroutine performs the mat_00 integral, z integral and sum over species in order to reduce from 6D to 2D
  subroutine nlt_integral(v6d,v2d)
    !This subroutine performs the mat_00 integral and sum over species

    implicit none
    complex, dimension(0:nx0-1,-1*(nky0-1):nky0-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: v6d
    real, dimension(0:nx0-1,-1*(nky0-1):nky0-1), intent(out) :: v2d
    real :: v2d_temp
    complex, dimension(0:nx0-1,-1*(nky0-1):nky0-1,lk1:lk2) :: v3d
    integer :: i,j,k,l,m,n,ierr

    v3d=(0.0,0.0)
     do n=ln1,ln2
        do m=lm1,lm2
           do l=ll1,ll2
              do k=lk1,lk2
                 call axpy_ij(nx0*(2*nky0-1),mat_00(pi1,pj1,k,l,m),v6d(:,:,k,l,m,n),&
                      &v3d(:,:,k))
              enddo
           enddo
        enddo
     enddo

    call my_complex_sum_vwspec(v3d,size(v3d)) 


    do j=-(nky0-1),nky0-1
       do i=0,nx0-1
        v2d(i,j)=real(sum(v3d(i,j,:)*geom%jacobian(pi1,pj1,:)/(real(nz0)*geom%avg_jaco),1) )
        call mpi_allreduce(v2d(i,j),v2d_temp,1,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)
        v2d(i,j)=v2d_temp
      end do
    end do

  end subroutine nlt_integral

  !>This subroutine calculates nonlinear transfer functions for POD modes
  subroutine diag_nlt_pod

    implicit none
    integer :: i,ierr,j,jj
    real :: nlt_rhs_2d(li1:li2,-1*(nky0-1):nky0-1)                                 !Integrated over z,v,mu
    real :: nlt_tot(20,num_nlt_pod_modes+1)  !+1 is for residual
    complex :: cfgamma_g(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: cfgamma_chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    real :: kx0,ky0           !kx and ky values for which to calculate coupling coefficient
    integer :: i0,j0          !Indices corresponding to above
    integer :: k,l,m,n
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields) :: fields_1
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2) :: f_1    
    complex, allocatable, dimension(:,:,:,:,:,:) :: h_1,apar_bar
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms
    real :: dE_z(lk1:lk2)
    real :: rhs_temp,r_temp,rhs_tot(20,num_nlt_pod_modes+1),&
            en_k(20,num_nlt_pod_modes+1),&
            dendt_k(20,num_nlt_pod_modes+1),tot_out(20,(num_nlt_pod_modes+1)*4)
    complex :: e_tot(11)
    complex :: g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: g_nlt_pod(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: hn_t,hn_all(num_nlt_pod_modes)
    character(len=3) :: char_n_pod

    write(char_n_pod,"(i3.3)") 4*(num_nlt_pod_modes+1)+1
    allocate(h_1(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    if(n_fields.gt.1) allocate(apar_bar(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))

!Get only fields right now from g_1 - need cfgamma for modes
    if (arakawa_zv) then
       call calc_aux_fields(g_1,fields_1,f_1,.true.,h_1) 
    else
       call calc_aux_fields(g_1,fields_1,f_1,.false.) 
    endif

    !construct chi
    do k=lk1,lk2
      do l=ll1,ll2
        do m=lm1,lm2
          do n=ln1,ln2
            call gyro_average_ff(fields_1(:,:,k,1),chi(:,:,k,l,m,n),k,m,n)
            if(n_fields.gt.1) call gyro_average_ff(fields_1(:,:,k,2),apar_bar(:,:,k,l,m,n),k,m,n)
          end do
        end do
      end do
    end do    

   if(n_fields.gt.1) then
     do l=ll1,ll2
      do n=ln1,ln2
          chi(:,:,lk1:lk2,l,lm1:lm2,n)= chi(:,:,lk1:lk2,l,lm1:lm2,n)-&
              (2.0*spec(n)%temp/spec(n)%mass)**0.5*vp(l)*apar_bar(:,:,lk1:lk2,l,lm1:lm2,n)
      end do
     end do
   end if

   if(allocated(apar_bar)) deallocate(apar_bar)

    do i=1,num_nlt_modes  !Cycle through wavenumbers 
      if(mype==0) write(nlt_pod_out(i)) time
       kx0=kx(kx_nlt_ind(i))
       ky0=ky(ky_nlt_ind(i))
       i0=kx_nlt_ind(i)
       j0=ky_nlt_ind(i)

      do j=1,num_nlt_pod_modes+1 !Cycle through POD modes for each wavenumber, '+1' is for residual

        if(j.ne.num_nlt_pod_modes+1) then
          g_nlt_pod=cmplx(0.0,0.0)  !This one has only the kx_center mode
          g_nlt_pod(kx_nlt_ind(i),ky_nlt_ind(i),:,:,:,:)=nlt_pod_mode(1,:,:,:,:,i,j)
          !Need g_temp with all connections used
          call fill_connections(nlt_pod_mode_wcons(0:nx0_nlt_pod(i)-1,:,:,:,:,i,j),&
               &g_temp,kx_nlt_ind(i),ky_nlt_ind(i),nx0_nlt_pod(i))
          !Project out time coefficient of mode of interest (note: modes are already normalized for this scalar product)
          !call spro_of_choice(1,g_temp,g_1,li1,li2,lj1,lj2,hn_t)
          call spro_of_choice(2,g_temp,g_1,li1,li2,lj1,lj2,hn_t)
          !if(mype==0) write(*,*) "i,j,hn_t",i,j,hn_t
          g_temp=g_nlt_pod*hn_t  !Now g_temp has the mode of interest with its appropriate time weight at kx_center
          hn_all(j)=hn_t
       else !Now get nlt for g_1-sum(hn*gn) i.e. the residual distribution function
          g_nlt_pod=cmplx(0.0,0.0)  !This one has only the kx_center mode
          do jj=1,num_nlt_pod_modes 
            g_nlt_pod(i0,j0,:,:,:,:)=g_nlt_pod(i0,j0,:,:,:,:)+hn_all(jj)*nlt_pod_mode(1,:,:,:,:,i,jj)
          end do
          !The following four lines shouldn't be necessary since I only use the i0,j0 modes
          !g_temp=cmplx(0.0,0.0)
          !g_temp(i0,j0,:,:,:,:)=cmplx(1.0,0.0)
          !g_temp=g_temp*g_1 !Get g_1 only at the right wavenumber
          !g_temp=g_temp-g_nlt_pod 
          g_temp=g_1-g_nlt_pod
       end if

        !Now get cfgamma for modes
   !Get cfgamma (same as in get_energy_terms)
   if (arakawa_zv) then 
      call calc_aux_fields(g_1,fields_1,f_1,.true.,h_1) 
      call get_cfgamma_h(h_1,cfgamma_chi)
   else 
      call calc_aux_fields(g_1,fields_1,f_1,.false.) 
      call get_cfgamma_f(f_1,fields_1,cfgamma_chi)
   endif

       !now get cfgamma_g - i.e. only g part
       call get_cfgamma_g(g_temp,cfgamma_g)
       !now get cfgamma_chi - i.e. total - cfgamma_g
       cfgamma_chi=cfgamma_chi-cfgamma_g

      !Now calculate nlt's
       !Get the nlt
       call construct_nlt(i,i0,j0,kx0,ky0,cfgamma_g,cfgamma_chi,chi,nlt_rhs_2d)
       
       nlt_tot(i,j)=sum(sum(nlt_rhs_2d,1),1)

       !output
       if(mype==0) write(nlt_pod_out(i)) nlt_rhs_2d

       !The following outputs the entire linear RHS for energy equation for POD modes at kx0,ky0.  
       !This should balance (on time-average in saturated state) with the sum 
       !over all k of the nonlinear transfer function (nlt_tot)
       if(j.ne.num_nlt_pod_modes+1) then
         !Fill g_temp with connections again (need to use connections since parallel boundary condition (for hyp_z, etc.) uses other kx)
         call fill_connections(nlt_pod_mode_wcons(0:nx0_nlt_pod(i)-1,:,:,:,:,i,j),g_temp,kx_nlt_ind(i),ky_nlt_ind(i),nx0_nlt_pod(i))
         !if(mype==0) write(*,*) "i,j,hn_t^2,e_tot(1)",i,j,real(conjg(hn_t)*hn_t),real(e_tot(1))
         !This gets the dE terms that we want - however we only want kx_center (see below)
         call get_energy_terms_sp(hn_t*g_temp,g_1,e_tot,dE_terms)
       else  !Get residual
         g_nlt_pod=cmplx(0.0,0.0)
         do jj=1,num_nlt_pod_modes
           !Fill g_temp with connections again
           call fill_connections(nlt_pod_mode_wcons(0:nx0_nlt_pod(i)-1,:,:,:,:,i,jj),&
                &g_temp,kx_nlt_ind(i),ky_nlt_ind(i),nx0_nlt_pod(i))
           g_nlt_pod=g_nlt_pod+hn_all(jj)*g_temp 
         end do
         g_temp=g_1-g_nlt_pod 
         call get_energy_terms_sp(g_temp,g_1,e_tot,dE_terms)
       end if

       rhs_tot(i,j) = 0.0
       en_k(i,j)    = 0.0
       dendt_k(i,j) = 0.0
      
       !This gets total RHS at only kx_center
       dE_z(:)=real(dE_terms(i0,j0,:,2)-dE_terms(i0,j0,:,9))
       rhs_temp=real(sum(dE_z*geom%jacobian(pi1,pj1,:)/(real(nz0)*geom%avg_jaco),1))  !total linear  RHS
       call mpi_allreduce(rhs_temp,r_temp,1,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)
       rhs_tot(i,j)=r_temp
       rhs_temp=sum(real(dE_terms(i0,j0,:,1))*geom%jacobian(pi1,pj1,:)/(real(nz0)*geom%avg_jaco),1)  !energy
       call mpi_allreduce(rhs_temp,r_temp,1,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)
       en_k(i,j)=r_temp
       dendt_k(i,j)=(en_k(i,j)-last_energy_nlt_pod(i,j))/real(istep_nlt*dt)
       tot_out(i,4*(j-1)+1:4*(j-1)+4) =(/nlt_tot(i,j),rhs_tot(i,j),en_k(i,j),dendt_k(i,j)/)

     end do !num_nlt_pod_modes
     !Output time traces of 1) total nonlinear, 2) total linear, 3) total energy, 4) time derivative of energy
     if(mype==0) write(nlt_info_pod(i),'('//char_n_pod//'es16.6)') time,tot_out(i,:)

   end do !num_nlt_modes

   last_energy_nlt_pod=en_k

  end subroutine diag_nlt_pod

  subroutine fill_connections(g_kxlin,g_kxnl,i_in,j_in,ncon)
    integer, intent(IN) :: ncon
    complex, intent(in) :: g_kxlin(0:ncon-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex, intent(out) :: g_kxnl(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    integer, intent(in) :: i_in,j_in
    integer :: i,kx_nl_ind
    
    !open(unit=100,file=trim(diagdir)//'/kx_indices.dat',status='unknown')
    g_kxnl=cmplx(0.0,0.0)
    !kx>0
    do i=0,ncon/2
      kx_nl_ind=i_in+i*pb_xshift(j_in) 
      if(kx(kx_nl_ind).lt.0.0) stop "Error in fill_connections."
      g_kxnl(kx_nl_ind,j_in,:,:,:,:)=g_kxlin(i,:,:,:,:)
    !  if(time.lt.1.0.and.mype==0) write(100,*) "kx<0:kx(ind)",kx(kx_nl_ind)
    !  if(time.lt.1.0.and.mype==0) write(100,*) "kx>0:i_in,j_in,nl_index,lin_index"
    !  if(time.lt.1.0.and.mype==0) write(100,*) i_in,j_in,kx(kx_nl_ind),i 
    end do
    !now kx<0
    do i=1,ncon/2
      kx_nl_ind=i_in-i*pb_xshift(j_in)  
      kx_nl_ind=nx0+kx_nl_ind
      if(kx(kx_nl_ind).ge.0.0) stop "Error in fill_connections."
      g_kxnl(kx_nl_ind,j_in,:,:,:,:)=g_kxlin(ncon-i,:,:,:,:)
    !  if(time.lt.1.0.and.mype==0) write(100,*) "kx<0:kx(ind)",kx(kx_nl_ind)
    !  if(time.lt.1.0.and.mype==0) write(100,*) "kx<0:i_in,j_in,nl_index,lin_index"
    !  if(time.lt.1.0.and.mype==0) write(100,*) i_in,j_in,kx_nl_ind,ncon-i
    end do
    !close(100)
    
  end subroutine fill_connections

  !>This subroutine gets the terms of interest in the energy equation as a scalar product.
  subroutine get_energy_terms_sp(g_1_in,g_2_in,e_tot,dE_terms)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: g_2_in
    complex, allocatable, dimension(:,:,:,:,:,:) :: g_temp,rhs,cfgamma,f_1_in,f_2_in,cfgamma_g
    complex, allocatable, dimension(:,:,:,:,:,:) :: h_2_in
    complex, allocatable, dimension(:,:,:) :: rhs_block
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: fields_1_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: fields_2_in
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms

    !complex,dimension(li1:li2,lj1:lj2,lk1:lk2):: erhs1
    complex :: e_tot(11)
    logical :: calc_this_term(11)
    integer :: lbg1, lbg2, iblock, iterm

!    if (nblocks.ne.1) stop 'get_energy_terms is currently only working with nblocks=1'
    lbg0 = lklmn0 / nblocks

    !    write(*,*) "Initialization"
    allocate(rhs(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(cfgamma(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(cfgamma_g(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(f_1_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    allocate(f_2_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    if (arakawa_zv) then
       allocate(h_2_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    endif

    e_tot= 0.0 
    dE_terms = (0.0,0.0)
    calc_this_term(1:2) = .false.
    
    call calc_aux_fields(g_1_in,fields_1_in,f_1_in,.false.)
    call get_cfgamma_f(f_1_in,fields_1_in,cfgamma)

    if (arakawa_zv) then 
       call calc_aux_fields(g_2_in,fields_2_in,f_2_in,.true.,h_2_in)
    else
       call calc_aux_fields(g_2_in,fields_2_in,f_2_in,.false.)
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
    rhs = 0.5*cfgamma*g_2_in   !To account for the complex conjugate for the RHS terms
    call energy_integral(rhs,dE_terms(:,:,:,1))

    !Total dE/dt
    !need g_temp because some compilers don't like sending an input variable 
    ! (g_1_in) to a subroutine that might modify it (calc_rhs_only)
    allocate(g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    g_temp=g_2_in
    rhs=(0.0,0.0)
    if(arakawa_zv)then
       call calc_rhs_only(f_2_in,g_temp,fields_2_in,rhs,2,h_2_in)
    else
       call calc_rhs_only(f_2_in,g_temp,fields_2_in,rhs,2)
    endif
    !collisions are taken out of calc_rhs_only in case of operator splitting
    if(coll_split) then
       call equ_collisions(f_2_in,rhs,replace_rhs=.false.)
    endif
    deallocate(g_temp)
    rhs = cfgamma*rhs
    call energy_integral(rhs,dE_terms(:,:,:,2))

    !------------------ INDIVIDUAL TERMS -------------------------------
    ! ---------------------- nonblocked terms ---------------------------------
    !collisions
    if (collision_op.ne.'none') then
       call equ_collisions(f_2_in,rhs,.true.)
       rhs = cfgamma*rhs
       call energy_integral(rhs,dE_terms(:,:,:,5))
    endif
    deallocate(rhs)

PERFON('en_block')
    !------------------- strip mining terms ---------------------------
    allocate(rhs_block(li1:li2,lj1:lj2,1:lbg0))
    rhs_block = (0.0,0.0)
    Do iblock=1,nblocks
       lbg1 =(iblock-1)*lbg0+1
       lbg2 =lbg1+lbg0-1

       do iterm=3,10
          calc_this_term(iterm) = .true.
          select case(iterm) 
             ! ---------------- SOURCE/DRIVE TERMS -------------------
          case(3)
             call add_dchidxy_orig(fields_2_in,ptr_bar_emfields,&
                  &rhs_block,ptr_barchi,ptr_dbarchidxy,lbg1,lbg2,.false.)
          case(4)
             if (.not.xy_local) then
                call add_f1_sources_block(rhs_block,lbg1,lbg2,0)
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
                   call add_hypv_ak(h_2_in,rhs_block,lbg1,lbg2)
                else
                   call add_hypv_ak(f_2_in,rhs_block,lbg1,lbg2)
                endif
             else
                call add_hypv_block(f_2_in,rhs_block,lbg1,lbg2)
             endif
          case(7)
             !hyp_z
             if (arakawa_zv) then
                if (hyp_on_h) then
                   call add_hypz_ak(h_2_in,rhs_block,lbg1,lbg2)  
                   if (hypz_compensation) call equ_comp_hypz(fields_2_in,ptr_bar_emfields,rhs_block,lbg1,lbg2)
                else
                   call add_hypz_ak(f_2_in,rhs_block,lbg1,lbg2) 
                endif
             else
                call add_hypz_block(f_2_in,rhs_block,lbg1,lbg2)
             endif
          case(8)
             if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then
                if (arakawa_zv) then
                   if (hyp_on_h) then
                      call add_dfdxy(h_2_in,rhs_block,lbg1,lbg2) !f_in is h
                   else
                      call add_dfdxy(f_2_in,rhs_block,lbg1,lbg2) !f_in is h
                   endif                
                else
                    call add_dfdxy(f_2_in,rhs_block,lbg1,lbg2) !f_in is f
                endif
             else
                calc_this_term(iterm) = .false.
             endif
             ! ----------------- Analytically vanishing terms -----------------
          case(9)
             if (rhs_nl) then
                if ((.not.xy_local).and.(nblocks.gt.1)) then
                   call add_dgdxy(g_2_in, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
                   rhs_block = 0.0
                endif
                
                CALL this_nonlinear_term%add(g_2_in, ptr_dgdxy, fields_2_in,&
                     &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, 2)
             else
                calc_this_term(iterm) = .false.
             endif
          case(10)
             ! Poisson bracket in z,vpar, i.e.
             ! cfgamma * (dfdzv+dchidzv) 
             !!this includes the (z,vpar) hyperdiffusion
             if (arakawa_zv) then
                call equ_dzv(h_2_in,rhs_block,lbg1,lbg2)
                if (hyp_on_h) then
                   if (hypz_compensation) call equ_comp_hypz(fields_2_in,ptr_bar_emfields,rhs_block,lbg1,lbg2)
                endif
             else
                call equ_dfdzv(f_2_in,rhs_block,lbg1,lbg2)
                call add_dchidz(fields_2_in, ptr_bar_emfields, rhs_block, lbg1, lbg2)
             end if
          case default
             stop 'wrong index in blockterm loop'
          end select

          if (calc_this_term(iterm)) then
             if(iterm==3) then !if this is the drive term then we need only cfgamma_g instead of all cfgamma
               call block_mult(rhs_block,cfgamma_g,lbg1,lbg2)
             else
               call block_mult(rhs_block,cfgamma,lbg1,lbg2)
             end if
             call add_energy_integral_block(rhs_block,dE_terms(:,:,:,iterm),lbg1,lbg2)
             rhs_block = (0.0,0.0)
          endif
       enddo
    enddo
    deallocate(rhs_block)
PERFOFF

    do iterm=1,10
       if (calc_this_term(iterm)) &
            &call my_complex_sum_vwspec(dE_terms(:,:,:,iterm),size(dE_terms(:,:,:,iterm)))
       !------------ direct space integration --------------------       
       if ((p_has_00_mode).and.(xy_local)) dE_terms(li1,lj1,lk1:lk2,iterm)=(0.0,0.0)
       call sum_int_3d_comp(dE_terms(:,:,:,iterm),e_tot(iterm))
       !compute dEdt_NC
       ! Note: nonlinearity should not be included
       !Warning, z,v poisson bracket does not cancel k by k and also includes dissipative boundary condition.
       if (iterm.gt.2.and.iterm.ne.9) dE_terms(:,:,:,11)=dE_terms(:,:,:,11)+dE_terms(:,:,:,iterm)
    enddo

if (hyp_on_h) then
    e_tot(10)= e_tot(10)-e_tot(6)-e_tot(7) !eliminated hyperdiffusions
    dE_terms(:,:,:,11)=dE_terms(:,:,:,11)-dE_terms(:,:,:,6)-dE_terms(:,:,:,7) !eliminated hyperdiffusions
else
    if (.not.arakawa_zv) e_tot(10)= e_tot(10)-e_tot(6)-e_tot(7) !eliminated hyperdiffusions
    if (.not.arakawa_zv) dE_terms(:,:,:,11)=dE_terms(:,:,:,11)-dE_terms(:,:,:,6)-dE_terms(:,:,:,7)
endif
    !rest
    e_tot(11)=e_tot(2)-sum(e_tot(3:10))

    deallocate(cfgamma,f_1_in,f_2_in)
    if (arakawa_zv) deallocate(h_2_in)

  end subroutine get_energy_terms_sp

  subroutine spro_of_choice(wsp,vec1,vec2,vli1,vli2,vlj1,vlj2,sp)
    integer, intent(in) :: wsp,vli1,vli2,vlj1,vlj2
    complex, intent(in), dimension(vli1:vli2,vlj1:vlj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)  :: vec1,vec2
    complex, intent(out) :: sp
    complex, dimension(vli1:vli2,vlj1:vlj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)  :: vec1_temp,vec2_temp
    complex :: e_tot(11)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,1:11) :: dE_terms
    integer :: vsize


    select case (wsp)
      case(1) !integral
        vsize=(vli2-vli1+1)*(vlj2-vlj1+1)*lk0*ll0*lm0*ln0
        call integral_weight(vec1,vec1_temp,vli1,vli2,vlj1,vlj2)
        call integral_weight(vec2,vec2_temp,vli1,vli2,vlj1,vlj2)
        call scalar_product_general(vec1_temp,vec2_temp,vsize,sp)
      case(2) !energy
        if(vli1.ne.li1.or.vli2.ne.li2.or.vlj1.ne.lj1.or.vlj2.ne.lj2) stop "Error in spro_of_choice."
        call get_energy_terms_sp(vec1,vec2,e_tot,dE_terms) 
        sp=e_tot(1)
        !e_temp=cmplx(0.0,0.0)    
        !do i=vli1,vli2
        !  do j=vlj1,vlj2 
        !    e_temp=e_temp+sum(dE_terms(i,j,:,1)*geom%jacobian(pi1,pj1,:)/(real(nz0)*geom%avg_jaco),1) !energy
        !  end do
        !end do
        !call mpi_allreduce(e_temp,sp,1,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_z,ierr)
    end select

  end subroutine spro_of_choice

  subroutine integral_weight(g_1_in,g_1_out,vli1,vli2,vlj1,vlj2)
    implicit none
    integer, intent(in) :: vli1,vli2,vlj1,vlj2
    complex, dimension(vli1:vli2,vlj1:vlj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_1_in
    complex, dimension(vli1:vli2,vlj1:vlj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: g_1_out
    integer :: k,l,m,n,pni
    complex, dimension(vli1:vli2,vlj1:vlj2) :: one_array
    logical :: f0_flag
    pni=pn1
    f0_flag=SVD_f0_flag

    one_array=cmplx(1.0,0.0) 
    !For now g_1_out will be used to store the integral weighting
    !if(evenx.ne.1) then
    do n=ln1,ln2
      if (pn0.gt.1) pni=n
      do k=lk1,lk2
         do m=lm1,lm2
            do l=ll1,ll2
              g_1_out(:,:,k,l,m,n)=one_array(:,:)*sqrt(mat_00(pi1,pj1,k,l,m)*geom%jacobian(pi1gl,pj1,k))
              if(f0_flag) then
                g_1_out(:,:,k,l,m,n)=g_1_out(:,:,k,l,m,n)/sqrt(fm(pi1,pj1,k,l,m,pni))
              end if
            enddo
         enddo
      enddo
    enddo
    !else
    !  write(*,*) "Error.  evenx==1."
    !  stop
    !end if

    g_1_out=g_1_out/sqrt((real(nz0)*geom%avg_jaco))
    g_1_out=g_1_in*g_1_out


  end subroutine integral_weight

  subroutine integral_unweight(g_1_in,g_1_out,vli1,vli2,vlj1,vlj2)
    implicit none
    integer, intent(in) :: vli1,vli2,vlj1,vlj2
    complex, dimension(vli1:vli2,vlj1:vlj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_1_in
    complex, dimension(vli1:vli2,vlj1:vlj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: g_1_out
    integer :: k,l,m,n,pni
    complex, dimension(vli1:vli2,vlj1:vlj2) :: one_array
    logical :: f0_flag

    pni=pn1

    f0_flag=SVD_f0_flag
    one_array=cmplx(1.0,0.0) 
    !For now g_1_out will be used to store the integral weighting
    !if(evenx.ne.1) then
    do n=ln1,ln2
      if (pn0.gt.1) pni=n
      do k=lk1,lk2
         do m=lm1,lm2
            do l=ll1,ll2
              g_1_out(:,:,k,l,m,n)=one_array(:,:)/(sqrt(mat_00(pi1,pj1,k,l,m)*geom%jacobian(pi1gl,pj1,k)))
             !var(o,n) = sum( sum( Sum(mom(:,:,:,o,n),1),1 )*geom%jacobian(pi1,pj1,:) )  
              if(f0_flag) then
                g_1_out(:,:,k,l,m,n)=g_1_out(:,:,k,l,m,n)*sqrt(fm(pi1,pj1,k,l,m,pni))
              end if
            enddo
         enddo
      enddo
    enddo
    !else
    !  write(*,*) "Error.  evenx==1."
    !  stop
    !end if

    g_1_out=g_1_out*sqrt((real(nz0)*geom%avg_jaco))
    g_1_out=g_1_in*g_1_out


  end subroutine integral_unweight

  subroutine scalar_product_general(v1,v2,vsize,sp)

    implicit none
    integer, intent(in) :: vsize
    complex, intent(in), dimension(vsize) :: v1,v2
    complex, intent(out) :: sp
    complex :: sptemp
    integer :: ierr

    sptemp=sum(conjg(v1)*v2)
    call mpi_allreduce(sptemp,sp,1,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)

  end subroutine scalar_product_general
#endif
 

End Module diagnostics_extended
