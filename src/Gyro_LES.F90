#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
! module containing all the FGS group simulations routines
! dynamic procedure routine
! Free energy spectral analysis (kx,ky,kperp ...)
! Free energy consevation analysis vs time
! Transfer function
Module Gyro_LES

Use par_mod
Use file_io, only: get_unit_nr
Use vel_space, only:  mat_00, fm
Use geometry
Use communications
use aux_fields
Use discretization
Use blockindex
Use prefactors
use axpy
use numerical_damping
USE calc_rhs,only: this_nonlinear_term

    Implicit None
    Public :: initialize_GyroLES, exec_GyroLES, finalize_GyroLES

    PRIVATE
    
    Character(Len=8) :: filestat='replace', filepos='rewind'
    Integer, Dimension(:), Allocatable  :: GyroLES_FILE_spec
    Complex, Dimension(:,:,:,:,:,:), Allocatable  :: filter_g_1
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: filter_f_ , filter_h_
    Complex, Dimension(:,:,:,:), Allocatable::  filter_emfields
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: filter_nlt    
    Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_filter_g_1 
    Complex, Dimension(:,:,:,:,:,:), Allocatable:: filter_nlt_filter_g_1  
    Complex, Dimension(:,:,:,:,:,:), Allocatable  :: t_fgs
    complex, dimension(:,:,:,:,:,:), allocatable :: m_x, m_y, tmp_submod
    complex, dimension(:,:,:,:), allocatable :: tmp_4D_dp


Contains

subroutine initialize_GyroLES
  
    implicit none
    
    integer :: n

    allocate(GyroLES_FILE_spec(ln1:ln2))
    
    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do n=ln1,ln2
            call get_unit_nr(GyroLES_FILE_spec(n))
            open(GyroLES_FILE_spec(n), file=trim(diagdir)//&
            &'/GyroLES'//'_'//trim(spec(n)%name)//trim(file_extension), form='formatted', &
            &status=filestat, position=filepos)
            write(GyroLES_FILE_spec(n),*) "#  1. time"
            write(GyroLES_FILE_spec(n),*) "#  2. hyp_x amplitude" 
            write(GyroLES_FILE_spec(n),*) "#  3. hyp_y amplitude"
            300     format ('#',"        1",3I12)
            write(GyroLES_FILE_spec(n),300) 2,3
        enddo
    endif
    
       adapt_istep_GyroLES = int(1.0/dt_max) 
       istep_GyroLES = min(istep_GyroLES,adapt_istep_GyroLES)

end subroutine initialize_GyroLES

subroutine exec_GyroLES

    Complex, dimension(:,:,:), pointer :: ptr2
    Complex, dimension(:,:,:,:), pointer :: ptr1, ptr3
    real,dimension(ln1:ln2) :: mymy, mxmx, mxmy, Tmx, Tmy, c_x, c_y
    real :: eps_alpha 
    Integer :: i, j, n, lbg1, lbg2, stage

    ptr1 => NULL()
    ptr2 => NULL()
    ptr3 => NULL()
    lbg1 = 1
    lbg2 = lklmn0
    stage = 0

    Allocate(filter_g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(filter_f_(lbi:ubi, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2))
    if (nonlin_h) Allocate(filter_h_(lbi:ubi, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2))

    Allocate(filter_emfields(lbi:ubi,lj1:lj2,lbz:ubz,1:n_fields))
    Allocate(nlt(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    Allocate(filter_nlt(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))     
    Allocate(nlt_filter_g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))  
    Allocate(filter_nlt_filter_g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))   
    Allocate(t_fgs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(m_x(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    Allocate(m_y(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(tmp_submod(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(tmp_4D_dp(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))

    !Apply the filter to g_1, compute the filtered fields
    Call filter(g_1, filter_g_1)

    if (.not.nonlin_h) then
       call calc_aux_fields(filter_g_1,filter_emfields,filter_f_,.false.)
       
       ! compute the T sub grid term
       nlt = cmplx(0.0,0.0)
       Call this_nonlinear_term%add(g_1, ptr1 ,emfields, ptr2, ptr3, nlt, lbg1, lbg2, stage)
    
       nlt_filter_g_1 = cmplx(0.0,0.0)
       Call this_nonlinear_term%add(filter_g_1, ptr1, filter_emfields, ptr2, ptr3, nlt_filter_g_1, lbg1, lbg2, stage)
    else
       call calc_aux_fields(filter_g_1,filter_emfields,filter_f_,.true.,h_out=filter_h_)
       
       ! compute the T sub grid term
       nlt = cmplx(0.0,0.0)
       Call this_nonlinear_term%add(h_, ptr1 ,emfields, ptr2, ptr3, nlt, lbg1, lbg2, stage)
    
       nlt_filter_g_1 = cmplx(0.0,0.0)
       Call this_nonlinear_term%add(filter_h_, ptr1, filter_emfields, ptr2, ptr3, nlt_filter_g_1, lbg1, lbg2, stage)
    endif
    Call filter(nlt, filter_nlt)
    Call filter(nlt_filter_g_1, filter_nlt_filter_g_1)     
    t_fgs = filter_nlt - filter_nlt_filter_g_1

    ! compute hyper-diff as models, with amplitude set to 1
    
    m_x = cmplx(0.0,0.0)
    m_y = cmplx(0.0,0.0)

    do j=lj1,lj2
        m_y(:,j,:,:,:,:) = djdiff(j)*filter_f_(li1:li2,j,lk1:lk2,ll1:ll2,lm1:lm2,:)
        do i=li1,li2
            m_x(i,j,:,:,:,:) = didiff(i)*filter_f_(i,j,lk1:lk2,ll1:ll2,lm1:lm2,:)
        enddo
    enddo

    ! terms needed for coefficients of the dynamic procedure
    tmp_submod = cmplx(0.0,0.0)
    tmp_4D_dp = cmplx(0.0,0.0)
    mxmx = 0.0
    tmp_submod = m_x*conjg(m_x)
    do n=ln1,ln2
        call integral_vpw(tmp_submod(:,:,:,:,:,n),tmp_4D_dp(:,:,:,n))
        if(p_has_00_mode) tmp_4D_dp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(tmp_4D_dp(:,:,:,n),mxmx(n))
    enddo
 
    tmp_submod = cmplx(0.0,0.0)
    tmp_4D_dp = cmplx(0.0,0.0)
    mxmy = 0.0
    tmp_submod = m_x*conjg(m_y)
    do n=ln1,ln2
        call integral_vpw(tmp_submod(:,:,:,:,:,n),tmp_4D_dp(:,:,:,n))
        if(p_has_00_mode) tmp_4D_dp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(tmp_4D_dp(:,:,:,n),mxmy(n))
    enddo
    
    tmp_submod = cmplx(0.0,0.0)
    tmp_4D_dp = cmplx(0.0,0.0)
    mymy = 0.0
    tmp_submod = m_y*conjg(m_y)
    do n=ln1,ln2
        call integral_vpw(tmp_submod(:,:,:,:,:,n),tmp_4D_dp(:,:,:,n))
        if(p_has_00_mode) tmp_4D_dp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(tmp_4D_dp(:,:,:,n),mymy(n))
    enddo
  
    tmp_submod = cmplx(0.0,0.0)
    tmp_4D_dp = cmplx(0.0,0.0)
    Tmx  = 0.0
    tmp_submod = t_fgs*conjg(m_x)
    do n=ln1,ln2
        call integral_vpw(tmp_submod(:,:,:,:,:,n),tmp_4D_dp(:,:,:,n))
        if(p_has_00_mode) tmp_4D_dp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(tmp_4D_dp(:,:,:,n),Tmx(n))
    enddo
 
    tmp_submod = cmplx(0.0,0.0)
    tmp_4D_dp = cmplx(0.0,0.0)
    Tmy = 0.0
    tmp_submod = t_fgs*conjg(m_y)
    do n=ln1,ln2
        call integral_vpw(tmp_submod(:,:,:,:,:,n),tmp_4D_dp(:,:,:,n))
        if(p_has_00_mode) tmp_4D_dp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(tmp_4D_dp(:,:,:,n),Tmy(n))
    enddo
   
    eps_alpha = 1.0 - (sqrt(fracx*fracy)**(-(hyp_x_order - 2./3.)))
   
    ! Calculation of the amplitude 
    c_x = 0.0
    c_y = 0.0
    do n=ln1,ln2
        c_x(n) = (mymy(n)*Tmx(n) - mxmy(n)*Tmy(n))/(mxmy(n)*mxmy(n) - mymy(n)*mxmx(n))/eps_alpha
        c_y(n) = (mxmx(n)*Tmy(n) - mxmy(n)*Tmx(n))/(mxmy(n)*mxmy(n) - mymy(n)*mxmx(n))/eps_alpha
    enddo


    !Update the values of the amplitude
    Do n=ln1,ln2
        if (c_x(n).gt.0.0) then
            hyp_x_spec(n) = c_x(n)
        else
            hyp_x_spec(n) = 0.0001
        endif
    
        if (c_y(n).gt.0.0) then
            hyp_y_spec(n) = c_y(n)
        else
            hyp_y_spec(n) = 0.0001
        endif
    enddo

    call my_barrier()
    
    !Write the amplitude to the file
    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do n=ln1,ln2
            write(GyroLES_FILE_spec(n),"(3ES12.4)") time, hyp_x_spec(n), hyp_y_spec(n)
           call flush(GyroLES_FILE_spec(n))
        enddo
    endif       
    
    Deallocate(filter_g_1)
    Deallocate(filter_f_)
    if (nonlin_h) Deallocate(filter_h_)
    Deallocate(filter_emfields)
    Deallocate(nlt)
    Deallocate(filter_nlt)
    Deallocate(nlt_filter_g_1)
    Deallocate(filter_nlt_filter_g_1)
    Deallocate(t_fgs)
    Deallocate(m_x)
    Deallocate(m_y)
    Deallocate(tmp_submod)
    deallocate(tmp_4D_dp)
    

End subroutine exec_GyroLES

subroutine finalize_GyroLES

    integer :: n
    
        IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then    
            Do n=ln1,ln2 
                Close (GyroLES_FILE_spec(n))
            enddo
        endif
        
    deallocate(GyroLES_FILE_spec)

end subroutine finalize_GyroLES


!!!**************************************************************************!!!
!!!********************* Subroutine that the module needs ************************!!! 

!!!**************************************************************************!!!
!!!*********************** Calculate the filter  ************************
subroutine filter(f_in,f_out)

    implicit none
    
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: f_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: f_out
    Integer :: fkx, fky
    
    ! Filter in the kx direction, it gives you the indice of the maximum value of the kx mode you keep
        fkx = int((fracx*hkx))
        ! Filter in the ky direction, it gives you the indice of the maximum value of the ky mode you keep 
        !fky = int((fracy*lj2))
        fky = int(fracy*(nky0-1))

        !Copy the input function to the output function FGS
        f_out = f_in
        !Filter FGS in the kx direction for modes higher than fkx and lower than -fkx 
 
        f_out( li1 + fkx + 1 : li2 - fkx,:,:,:,:,:) = CMPLX(0.,0.)
   
        IF (fky.ge.lj1 .and. fky.ge.lj2) then
        ELSE IF (fky.ge.lj1.and.fky.lt.lj2) then
            f_out(:,fky + 1 : lj2,:,:,:,:) = CMPLX(0.,0.)
        ELSE IF (fky.ge.lj1.and.fky.eq.lj2) then
        ELSE
            f_out(:,:,:,:,:,:)= CMPLX(0.,0.)
        ENDIF

End Subroutine filter

!!!***************************************************************************************
subroutine integral_vpw(v5d,v3d)
!This subroutine performs the mat_00 integral and doesnt have explicit species dependence

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2), intent(in) :: v5d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2),intent(out) :: v3d
    integer :: k,m,l
    
    v3d=cmplx(0.0,0.0)
        do k=lk1,lk2
            do m=lm1,lm2
                do l=ll1,ll2
                    call saxpy(2*li0*lj0,mat_00(pi1,pj1,k,l,m),v5d(li1,lj1,k,l,m),1,v3d(li1,lj1,k),1)
                enddo
            enddo
        enddo

    call my_complex_sum_vw(v3d,size(v3d)) 
 

end subroutine integral_vpw  

!!!******************************************************************************************
subroutine sum_int_3d(v3d,v1d)
    !This subroutine performs the z integral and sum over kx,ky

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(in) :: v3d
    real, intent(out) :: v1d
    real, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d_tmp
    real :: v1d_loc
    integer :: ierr

    v3d_tmp = 0.0
    v1d = 0.0
    v1d_loc = 0.0

    v3d_tmp = real(v3d)
    if ((x_local).and.(evenx.eq.1)) v3d_tmp(hkx+1,:,:) = 0.0
    if (p_has_0_mode) v3d_tmp(:,lj1,:) = 0.5*v3d_tmp(:,lj1,:)

    v1d_loc = 2.0*sum(sum(sum(v3d_tmp,1),1)*geom%jacobian(pi1,pj1,:))/&
            (real(nz0)*geom%avg_jaco)

    call mpi_allreduce(v1d_loc,v1d,1,MPI_REAL_TYPE,MPI_SUM, mpi_comm_xyz, ierr)

  end subroutine sum_int_3d
!!!******************************************************************************************

End Module Gyro_LES
