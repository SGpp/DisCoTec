#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
! Free energy spectral analysis (kx,ky,kperp ...)
Module diag_Gyro_LES_spec2d
  
  Use par_mod 
  Use file_io, only: get_unit_nr
  use communications
  use antenna, only: antenna_type
  use collisions, only: equ_collisions
  use diagnostics_energy, only: get_cfgamma_f,get_cfgamma_h,get_cfgamma_ant,&
       &energy_integral
  use dfdzv_terms 
  use dgdxy_terms
  use dzv_terms
  use dfdxy_terms, only: add_dfdxy
  use dchidz_term, only: add_dchidz
  use dchidxy_terms, only: add_dchidxy, add_dchidxy_orig
  use spatial_averages, only: sum_int_3d, sum_int_z
  use prefactors
  use numerical_damping
  USE calc_rhs, only: this_nonlinear_term
  use diag_Gyro_LES_common

  Implicit None
  
  public:: initialize_diag_spectral_twoD, finalize_diag_spectral_twoD, diag_spectral_twoD

  !****** Quantities for 2D spectral diags
  private
  Character(Len=8) :: filestat='replace', filepos='rewind'
  integer, dimension(8) :: fnumber ! Will store the file id's for all quantity/species couples
  integer, dimension(:,:),allocatable:: fnumber_spec ! Will store the file id's for all quantity/species couples
  integer, dimension(:), allocatable:: twoD_counter ! Stores the number of windows for average
  integer, dimension(:,:), allocatable:: twoD_counter_spec ! Stores the number of windows for average
  real, dimension(:,:,:), allocatable:: cumul_kxky
  real, dimension(:,:,:,:), allocatable:: cumul_kxky_spec
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: temp_rhs, temp_fe6d

contains

  !******************************************************************************************
  !************************************ entropy/electrostatic balance ***********************
  !******************************************************************************************

  !******************************************************************************************************
  !**************************************  Diag 2D ******************************************************
  !*******************************************************************************************************

  Subroutine initialize_diag_spectral_twoD

#ifdef with_extended_diags
    Character(len=20), dimension(8) :: label
    integer :: n,m,i,j,nx_index
    real, dimension(0:nx0-1) :: kxglob

    ! Arrays to cumulate for time average
    allocate(cumul_kxky(li1:li2,lj1:lj2,1:8))
    allocate(twoD_counter(8))
    allocate(cumul_kxky_spec(li1:li2,lj1:lj2,ln1:ln2,1:24))
    allocate(fnumber_spec(ln1:ln2,1:24))
    allocate(twoD_counter_spec(ln1:ln2,1:24))

    fnumber=0
    fnumber_spec=0
    twoD_counter=0
    twoD_counter_spec=0
    cumul_kxky=0
    cumul_kxky_spec=0


    ! We initialize the label to be appended to the name of file
    label(1)='fe'
    label(2)='drive'
    label(3)='hyp_v'
    label(4)='hyp_z'
    label(5)='hyp_kperp'
    label(6)='coll'
    label(7)='parallel'
    label(8)='curvature'

    ! We will write kx and ky in the header of the file 	    
    ! We create kxglob without hkx+1 mode and invert it	    

    if (evenx.eq.1) then
    kxglob(0:nx0/2-2)=kx(nx0/2+1:nx0-1)
    kxglob(nx0/2-1:nx0-1)=kx(0:nx0/2)
    nx_index =  nx0-2
    
    else
    kxglob(0:nx0/2-1)=kx(nx0/2+1:nx0-1)
    kxglob(nx0/2:nx0-1)=kx(0:nx0/2)
    nx_index =  nx0-1
    endif

    if (mype.eq.0) then
       do m=1,8
          call get_unit_nr(fnumber(m))
          OPEN(fnumber(m), file=trim(diagdir)//'/Spectral_2D_'&
               //trim(adjustl(label(m)))//trim(file_extension),form = 'FORMATTED', &
               status=filestat, position=filepos)         

          ! We can directly write preliminary numbers as headers	    
          
          Do i=0,nx_index
             write(fnumber(m),"(ES12.4)") kxglob(i)
          enddo
          Do  j=0,nky0-1     
             write(fnumber(m),"(ES12.4)") ky(j)
          enddo
       enddo

    endif

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do n=ln1,ln2
            do m=1,8
                call get_unit_nr(fnumber_spec(n,m))
                OPEN(fnumber_spec(n,m), file=trim(diagdir)//'/Spectral_2D_W_'&
                    //trim(adjustl(label(m)))//'_'//trim(spec(n)%name)//trim(file_extension),form = 'FORMATTED', &
                    status=filestat, position=filepos)         

                ! We can directly write preliminary numbers as headers	    
                Do i=0,nx_index
                    write(fnumber_spec(n,m),"(ES12.4)") kxglob(i)
                enddo
                Do  j=0,nky0-1     
                    write(fnumber_spec(n,m),"(ES12.4)") ky(j)
                enddo
            enddo
        enddo
        do n=ln1,ln2
            do m=9,16
                call get_unit_nr(fnumber_spec(n,m))
                OPEN(fnumber_spec(n,m), file=trim(diagdir)//'/Spectral_2D_E_'&
                    //trim(adjustl(label(m-8)))//'_'//trim(spec(n)%name)//trim(file_extension),form = 'FORMATTED', &
                    status=filestat, position=filepos)         

                ! We can directly write preliminary numbers as headers	    
                Do i=0,nx_index
                    write(fnumber_spec(n,m),"(ES12.4)") kxglob(i)
                enddo
                Do  j=0,nky0-1     
                    write(fnumber_spec(n,m),"(ES12.4)") ky(j)
                enddo
            enddo
        enddo
        do n=ln1,ln2
            do m=17,24
                call get_unit_nr(fnumber_spec(n,m))
                OPEN(fnumber_spec(n,m), file=trim(diagdir)//'/Spectral_2D_FE_'&
                    //trim(adjustl(label(m-16)))//'_'//trim(spec(n)%name)//trim(file_extension),form = 'FORMATTED', &
                    status=filestat, position=filepos)         

                ! We can directly write preliminary numbers as headers	    
                Do i=0,nx_index
                    write(fnumber_spec(n,m),"(ES12.4)") kxglob(i)
                enddo
                Do  j=0,nky0-1     
                    write(fnumber_spec(n,m),"(ES12.4)") ky(j)
                enddo
            enddo
        enddo
    endif
#endif

  End Subroutine initialize_diag_spectral_twoD

  Subroutine finalize_diag_spectral_twoD

#ifdef with_extended_diags
    Implicit None
    Integer:: n,m

    if (mype.eq.0) then
       do m=1,8
          close(fnumber(m))
       enddo
    endif

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do n=ln1,ln2
            do m=1,24
                close(fnumber_spec(n,m))
            enddo
        enddo
    endif

    deallocate(fnumber_spec)
    deallocate(twoD_counter)
    deallocate(twoD_counter_spec)
    deallocate(cumul_kxky)
    deallocate(cumul_kxky_spec)
#endif

  End Subroutine finalize_diag_spectral_twoD


!!!******************************************************************************************
!!! Subroutine that goes from 6d to 2d     ***************************************************
!!!*******************************************************************************************

  Subroutine diag_spectral_twoD

#ifdef with_extended_diags
    implicit none

    real, dimension(li1:li2,lj1:lj2,1:8) :: v3d
    real, dimension(li1:li2,lj1:lj2,ln1:ln2,1:24) :: v4d
    integer :: m,n,lbg1,lbg2
    complex, dimension(:,:,:),pointer :: ptr2
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4

    Allocate(g_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(temp_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(temp_fe6d(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma3(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
       

    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1 
    lbg2 = lklmn0  
    v4d=0.0
    v3d=0.0

    if (arakawa_zv) then
       call get_cfgamma_h(h_,cfgamma)
       call get_cfgamma(f_,emfields,cfgamma1,cfgamma2)
    else
       call get_cfgamma_f(f_,emfields,cfgamma) !because we compute f
       call get_cfgamma(f_,emfields,cfgamma1,cfgamma2)
    endif
    !with antenna, the free energy as well as the free energy per species will contain the 
    !antenna contribution -- however, with electromagnetics the entropy and electrostatic parts
    !are not fully separated -- the entropy contains picks up some Bperp terms, the electrostatic part 
    !picks up some Bpar terms!
    if (any(antenna_type.eq.(/2,3/))) then
       allocate(cfgamma_ant(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       call get_cfgamma_ant(cfgamma_ant)
       !add antenna term to cfgamma for calculation of total energy
       cfgamma=cfgamma+cfgamma_ant
       !add antenna term to cfgamma1 (as it already contains some electromagnetic contributions)
       !to have it also in the species-dependent energy
       cfgamma1=cfgamma1+cfgamma_ant
       deallocate(cfgamma_ant)
    endif

    !------------------------------------------------------------------
    !Dissipation  v / 3
    temp_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then 
       if (hyp_on_h) then 
          call add_hypv_ak(h_,temp_rhs,lbg1,lbg2)
       else
          call add_hypv_ak(f_,temp_rhs,lbg1,lbg2)
       endif
    else
        call add_hypv_block(f_,temp_rhs,lbg1,lbg2)
    endif

    temp_fe6d = cfgamma*temp_rhs
    call twoD_array(temp_fe6d,v3d(:,:,3))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,19))
    enddo
    temp_fe6d = cfgamma1*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,3))
    enddo
    temp_fe6d = cfgamma2*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,11))
    enddo

    !--------------------------------------------------------------------------------
    !Dissipation  z / 4
    temp_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then
       if (hyp_on_h) then 
          call add_hypz_ak(h_,temp_rhs,lbg1,lbg2)
          if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,temp_rhs,lbg1,lbg2)
       else  
          call add_hypz_ak(f_,temp_rhs,lbg1,lbg2)
       endif  
    else
        call add_hypz_block(f_,temp_rhs,lbg1,lbg2)
    endif

    temp_fe6d = cfgamma*temp_rhs
    call twoD_array(temp_fe6d,v3d(:,:,4))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,20))
    enddo
    temp_fe6d = cfgamma1*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,4))
    enddo
    temp_fe6d = cfgamma2*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,12))
    enddo

    !--------------------------------------------------------------------------------
    !Dissipation  kperp / 5

    temp_rhs = cmplx(0.0,0.0)
    if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then
       !/todo:check if hyp_on_h should affect perpendicular hyperdiffusions!
       !if (hyp_on_h) then
       !   call add_dfdxy(h_,temp_rhs,lbg1,lbg2) !f_in is f
       !else
          call add_dfdxy(f_,temp_rhs,lbg1,lbg2) !f_in is f
       !endif
    endif

    temp_fe6d = cfgamma*temp_rhs
    call twoD_array(temp_fe6d,v3d(:,:,5))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,21))
    enddo
    temp_fe6d = cfgamma1*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,5))
    enddo
    temp_fe6d = cfgamma2*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,13))
    enddo


    !--------------------------------------------------------------------------------
    !Collisions / 6
    if (collision_op.ne.'none') then
       call equ_collisions(f_,temp_rhs,replace_rhs=.true.)
    endif

    temp_fe6d = cfgamma*temp_rhs
    call twoD_array(temp_fe6d,v3d(:,:,6))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,22))
    enddo
    temp_fe6d = cfgamma1*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,6))
    enddo
    temp_fe6d = cfgamma2*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,14))
    enddo

    !--------------------------------------------------------------------------------
    !dgdxy and dchidxy terms with both curvature and gradient contributions (.true.)
    temp_rhs = cmplx(0.0,0.0)
    call add_dgdxy(g_1, temp_rhs, ptr1, pdg1di, pdg1dj,lbg1,lbg2)
    call add_dchidxy_orig(emfields,ptr4, temp_rhs,ptr2,ptr3,lbg1,lbg2,.true.)

    !exclude curvature (.false.)   
    g_rhs  = cmplx(0.0,0.0)
    call add_dchidxy_orig(emfields,ptr4, g_rhs,ptr2,ptr3,lbg1,lbg2,.false.)

    !output arrays for gradient term
    temp_fe6d = cfgamma*g_rhs
    call twoD_array(temp_fe6d,v3d(:,:,2))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,18))
    enddo
    temp_fe6d = cfgamma1*g_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,2))
    enddo
    temp_fe6d = cfgamma2*g_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,10))
    enddo

    !Curvature / 8
    !substract the gradient contribution  
    temp_rhs = temp_rhs - g_rhs   


    !output arrays for curvature term
    temp_fe6d = cfgamma*temp_rhs
    call twoD_array(temp_fe6d,v3d(:,:,8))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,24))
    enddo
    temp_fe6d = cfgamma1*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,8))
    enddo
    temp_fe6d = cfgamma2*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,16))
    enddo

    !--------------------------------------------------------------------------------

    !Parallel / 7
    temp_rhs = cmplx(0.0,0.0)

    !Parallel term goes always with h
    if (arakawa_zv) then
       call equ_dzv(h_,temp_rhs,lbg1,lbg2)
       if (hyp_on_h) then
          if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,temp_rhs,lbg1,lbg2)
       endif
    else
       call equ_dfdzv(f_,temp_rhs,lbg1,lbg2)
       call add_dchidz(emfields, ptr4, temp_rhs, lbg1, lbg2)
    end if

    temp_fe6d = cfgamma*temp_rhs
    call twoD_array(temp_fe6d,v3d(:,:,7))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,23))
    enddo
    temp_fe6d = cfgamma1*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,7))
    enddo
    temp_fe6d = cfgamma2*temp_rhs
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,15))
    enddo

    !--------------------------------------------------------------------------------

    !Free energy / 1
    temp_fe6d(:,:,:,:,:,:) = cfgamma*g_1/2.
    call twoD_array(temp_fe6d,v3d(:,:,1))
    do n=ln1,ln2
       call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,1))
    enddo

    if (n_spec.gt.1) then
       call get_electrostatic(emfields,cfgamma3)
       call sum_int_z(real(cfgamma3(:,:,:,n)),v4d(:,:,n,9))
       do n=ln1,ln2
          v4d(:,:,n,17) = v4d(:,:,n,1) + v4d(:,:,n,9)
       enddo
    else
       temp_fe6d(:,:,:,:,:,:) = cfgamma2*g_1/2.
       do n=ln1,ln2
          call twoD_array_spec(temp_fe6d(:,:,:,:,:,n),v4d(:,:,n,9))
          v4d(:,:,n,17) = v4d(:,:,n,1)
       enddo
    endif

    deallocate(g_rhs,temp_rhs,temp_fe6d)
    deallocate(cfgamma)
    deallocate(cfgamma1)
    deallocate(cfgamma2)

    !At this point, all output is collected in v3d and v4d arrays
    do m=1,8
       ! We do proceed with cumulating if condition is verified
       if (itime.gt.0 .and. .not.( mod(itime,istep_fe_twoD*avg_window) .eq. 0) ) then
          cumul_kxky(:,:,m)=cumul_kxky(:,:,m)+v3d(:,:,m)     
          twoD_counter(m) = twoD_counter(m) +1   
       endif
       ! Then we do output files
       if  (itime.gt.0 .and. mod(itime,istep_fe_twoD*avg_window) .eq. 0 ) then
          cumul_kxky(:,:,m)=cumul_kxky(:,:,m)+v3d(:,:,m)     
          twoD_counter(m) = twoD_counter(m) +1   
          v3d(:,:,m)=cumul_kxky(:,:,m)/real(twoD_counter(m))
          call write2d(v3d(:,:,m),m)
          cumul_kxky(:,:,m) = 0
          twoD_counter(m)=0
       endif ! itime > 0 / twoD_counter...
    enddo ! loop over quant

    do n=ln1,ln2
        do m=1,24
            ! We do proceed with cumulating if condition is verified
            if (itime.gt.0 .and. .not.( mod(itime,istep_fe_twoD*avg_window) .eq. 0) ) then
            cumul_kxky_spec(:,:,n,m)=cumul_kxky_spec(:,:,n,m)+v4d(:,:,n,m)     
            twoD_counter_spec(n,m) = twoD_counter_spec(n,m) +1   
            endif
            ! Then we do output files
            if  (itime.gt.0 .and. mod(itime,istep_fe_twoD*avg_window) .eq. 0 ) then
                cumul_kxky_spec(:,:,n,m)=cumul_kxky_spec(:,:,n,m)+v4d(:,:,n,m)     
                twoD_counter_spec(n,m) = twoD_counter_spec(n,m) +1   
                v4d(:,:,n,m)=cumul_kxky_spec(:,:,n,m)/real(twoD_counter_spec(n,m))
                call write2d_spec(v4d(:,:,n,m),n,m)
                cumul_kxky_spec(:,:,n,m) = 0
                twoD_counter_spec(n,m)=0
            endif ! itime > 0 / twoD_counter...
        enddo ! loop over quant
    enddo ! loop over species

    deallocate(cfgamma3)
#endif
 
  End subroutine diag_spectral_twoD

!!!******************************************************************************************
!!!*********************** Write 2d  ************************
#ifdef with_extended_diags

  subroutine write2d(arr2d,m)

    integer, intent(in) :: m
    Real, Dimension(li1:li2,lj1:lj2), intent(in):: arr2d
    Real, Dimension(0:nx0-1, 0:nky0-1)  ::  dest2d, dest2d_linear
    Real, Dimension(0:nx0-1,lj1:lj2):: arr2dx
    Real :: dest1d_linear
    Integer :: i,j,ierr


    if ((my_pev+my_pew+my_pespec).eq.0) then
       ! Gather arr2d over x-Distribution on pex==0.
       do j = lj1, lj2
          Call mpi_gather(arr2d(li1,j), li0, MPI_REAL_TYPE,&
               arr2dx(0,j), li0, MPI_REAL_TYPE,&
               0, mpi_comm_x, ierr)
       enddo


       ! Gather arr2d over y-Distribution on pey==0.
       Call mpi_gather(arr2dx(0,lj1),size(arr2dx),MPI_REAL_TYPE,&
            dest2d(0,0),size(arr2dx),MPI_REAL_TYPE,0,mpi_comm_y,ierr)

       if (mype.eq.0) then

          write(fnumber(m),"(ES12.4)") time
          
          if (nonlinear) then
            write(fnumber(m),"(ES12.4)") 2*sum(sum(dest2d(0:nx0-1,1:nky0-1),1),1) + sum(dest2d(:,0),1)
          else
            dest2d_linear = dest2d  
            if ((x_local).and.(evenx.eq.1)) dest2d_linear(hkx+1,:) = 0.0
            if (p_has_0_mode) dest2d_linear(:,lj1) = 0.5*dest2d_linear(:,lj1)
            dest1d_linear = 2.0*sum(sum(dest2d_linear,1),1)
            write(fnumber(m),"(ES12.4)") dest1d_linear
          endif

          if (evenx.eq.1) then

            do i=nx0/2+1,nx0-1
                do j = 0,nky0-1
                    write(fnumber(m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo

            do i=0,nx0/2-1
                do j=0,nky0-1
                    write(fnumber(m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo
          else

            do i=nx0/2+1,nx0-1
                do j = 0,nky0-1
                    write(fnumber(m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo

            do i=0,nx0/2
                do j=0,nky0-1
                    write(fnumber(m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo
          endif 
         call flush(fnumber(m))
       endif
    endif

    call my_barrier()

  End Subroutine write2d

!!!*********************** Write 2d for species dependent  ************************

  subroutine write2d_spec(arr2d,n,m)

    integer, intent(in) :: n,m
    Real, Dimension(li1:li2,lj1:lj2), intent(in):: arr2d
    Real, Dimension(0:nx0-1, 0:nky0-1)  ::  dest2d, dest2d_linear
    real :: dest1d_linear
    Real, Dimension(0:nx0-1,lj1:lj2):: arr2dx
    Integer :: i,j,ierr


    if ((my_pev+my_pew).eq.0) then
       ! Gather arr2d over x-Distribution on pex==0.
       do j = lj1, lj2
          Call mpi_gather(arr2d(li1,j), li0, MPI_REAL_TYPE,&
               arr2dx(0,j), li0, MPI_REAL_TYPE,&
               0, mpi_comm_x, ierr)
       enddo


       ! Gather arr2d over y-Distribution on pey==0.
       Call mpi_gather(arr2dx(0,lj1),size(arr2dx),MPI_REAL_TYPE,&
            dest2d(0,0),size(arr2dx),MPI_REAL_TYPE,0,mpi_comm_y,ierr)

       if (mype.eq. pexyzvwspec(0,0,0,0,0,my_pespec) ) then
          write(fnumber_spec(n,m),"(ES12.4)") time

          if (nonlinear) then
            write(fnumber_spec(n,m),"(ES12.4)") 2*sum(sum(dest2d(0:nx0-1,1:nky0-1),1),1) + sum(dest2d(:,0),1)
          else
            dest2d_linear = dest2d  
            if ((x_local).and.(evenx.eq.1)) dest2d_linear(hkx+1,:) = 0.0
            if (p_has_0_mode) dest2d_linear(:,lj1) = 0.5*dest2d_linear(:,lj1)
            dest1d_linear = 2.0*sum(sum(dest2d_linear,1),1)
            write(fnumber_spec(n,m),"(ES12.4)") dest1d_linear
          endif

          if (evenx.eq.1) then

            do i=nx0/2+1,nx0-1
                do j = 0,nky0-1
                    write(fnumber_spec(n,m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo
            do i=0,nx0/2-1
                do j=0,nky0-1
                    write(fnumber_spec(n,m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo

        else
            do i=nx0/2,nx0-1
                do j = 0,nky0-1
                    write(fnumber_spec(n,m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo
            do i=0,nx0/2-1
                do j=0,nky0-1
                    write(fnumber_spec(n,m),"(ES12.4)") dest2d(i,j)
                enddo
            enddo
        endif
        call flush(fnumber_spec(n,m))
       endif
    endif

    call my_barrier()

  End Subroutine write2d_spec

!!!******************************************************************************************

  subroutine twoD_array(v6d,v2d_out)

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: v6d
    real, dimension(li1:li2,lj1:lj2), intent(out) :: v2d_out
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d
    real, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d_bis

    v3d = cmplx(0.0)
    v3d_bis = 0.0    
    v2d_out = 0.0

    call energy_integral(v6d,v3d)

    v3d_bis=real(v3d)
    if ((x_local).and.(evenx.eq.1)) v3d_bis(hkx+1,:,:) = 0.0
    if (p_has_0_mode) v3d_bis(li1,lj1,:) = 0.0

    call sum_int_z(v3d_bis,v2d_out)

  end subroutine twoD_array

!!!******************************************************************************************

  subroutine twoD_array_spec(v5d,v2d_out)

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2), intent(in) :: v5d
    real, dimension(li1:li2,lj1:lj2), intent(out) :: v2d_out
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d
    real, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d_bis

    v3d = cmplx(0.0)
    v3d_bis = 0.0    
    v2d_out = 0.0

    call integral_vpw(v5d,v3d)

    v3d_bis=real(v3d)
    if ((x_local).and.(evenx.eq.1)) v3d_bis(hkx+1,:,:) = 0.0
    if (p_has_0_mode) v3d_bis(li1,lj1,:) = 0.0

    call sum_int_z(v3d_bis,v2d_out)

  end subroutine twoD_array_spec
#endif


end Module diag_Gyro_LES_spec2d
