#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
#undef WITH_FLR
!> This modules contains routines for post-processing of 
!! fsa_moments terms
Module diagnostics_fsa_moments
  Use par_mod 
  Use file_io, only: get_unit_nr
  Use vel_space, only: fm, mat_00, calc_vsp_moment
  use calc_rhs, only: calc_rhs_only,ptr_bar_emfields, ptr_barchi, &
       ptr_dbarchidxy, ptr_dgdxy, rhs_nl, chi_block, this_nonlinear_term
  use geometry
  use communications
  use aux_fields, only: calc_aux_fields
  use collisions, only: equ_collisions
  use dfdzv_terms !, only: add_dfdz, add_dfdv, add_hypv, add_hypz
  use dgdxy_terms
  use dzv_terms
  use dfdxy_terms, only: add_dfdxy
  use dchidz_term, only: add_dchidz
  use dchidxy_terms, only: add_dchidxy, add_dchidxy_orig
  use f0_term 
  use sources_mod, only: add_f1_sources,add_kBuffer_explicit,explicit_buffer
  use spatial_averages
  use blockindex
  use axpy
  use prefactors
  use numerical_damping
  use boundaries, only: x_boundary_block, exchange_x
#ifdef WITHFUTILS
  use futils
  use hashtable
#endif
  Implicit None

  PUBLIC :: initialize_all_diags_fsa_moments, exec_all_diags_fsa_moments, &
       &finalize_all_diags_fsa_moments, mem_est_diag_fsa_moments,&
       &check_diag_fsa_moments

  PRIVATE 

  INTEGER,dimension(:),allocatable :: FSAMOMFILE
  Logical :: write_prof_pe
  integer::block_nr = 0
  integer::n_fsa_mats = 1
  integer::n_fsa_terms = 12
  REAL, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE,private :: fsa_moments_mat

Contains

!!!*******************public routines *******************************!!!

  SUBROUTINE check_diag_fsa_moments
    if (turbdeal.and.istep_fsa_moments.gt.0) then
       stop 'The fsa_moments diagnostic is not implemented for turbo dealiasing yet'
       !this would require some extra treatment of the nonlinearity
       !where several stage need to be called
    endif
    if (nonlin_h.and.istep_fsa_moments.gt.0) &
         stop 'fsa_moments diagnostic is not implemented for nonlinearity using "h" distribution!'
#ifndef with_extended_diags
    if (istep_fsa_moments.gt.0) then
       Write (*,'(A)') 'ERROR: set with_extended_diags in switches.h to use fsa_moments diagnostics'
       stop
    endif
#endif
  END SUBROUTINE check_diag_fsa_moments

  !> Estimates the memory required by fsa_moments diags
  !> todo complete this!
  Real Function mem_est_diag_fsa_moments(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc = 0.0

#ifdef with_extended_diags
    !rhs,g_temp (temporary!)
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*2*lijklmn0
    !f_in, h_in (temporary!)
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*lzvw0*ln0
    if (arakawa_zv) mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*lzvw0*ln0

    !fields_in,fsamom_terms
    mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*lij0*(n_fields+11)
#endif

    mem_est_diag_fsa_moments=mem_req_in+mem_loc

  End Function mem_est_diag_fsa_moments

!!!******************************************************************!!!

  SUBROUTINE initialize_all_diags_fsa_moments

#ifdef with_extended_diags
    if (istep_fsa_moments.gt.0) then
       call initialize_diag_fsa_moments
    end if
#endif

  END SUBROUTINE initialize_all_diags_fsa_moments

!!!******************************************************************!!!

  SUBROUTINE exec_all_diags_fsa_moments(itime,time)
    integer, intent(in) :: itime
    real, intent(in) :: time

#ifdef with_extended_diags
    IF(istep_fsa_moments.GT.0) THEN
       IF (MOD(itime,istep_fsa_moments).eq.0) CALL diag_fsa_moments
    END IF
#endif

  END SUBROUTINE exec_all_diags_fsa_moments

!!!******************************************************************!!!

  SUBROUTINE finalize_all_diags_fsa_moments

#ifdef with_extended_diags
    if (istep_fsa_moments.gt.0) then
       call finalize_diag_fsa_moments
    end if
#endif

  END SUBROUTINE finalize_all_diags_fsa_moments
!!!******************************************************************!!!

#ifdef with_extended_diags

!!!***********private subroutines ***********************************!!!
  subroutine initialize_diag_fsa_moments
    implicit none
    Integer :: n

    ALLOCATE(FSAMOMFILE(ln1:ln2))

    block_nr=0

    write_prof_pe = (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec))

    if (write_prof_pe) then
       ! Moment files for each species in a different file
       DO n=ln1,ln2
          call get_unit_nr(FSAMOMFILE(n))
          OPEN(FSAMOMFILE(n), file=trim(diagdir)//'/fsamom_'//&
               &trim(spec(n)%name)//''//trim(file_extension),&
               form='formatted', status='replace', position='rewind')  
          ! write file header
          write(FSAMOMFILE(n),"(4A)")      "#   x/a             (2)tot Q        ",&
               "(3)tot dQ       (4)coll         (5)dchidxy      (6)f1_sources   ",&
               "(7)f1_buffer    (8)f0 term      (9)v_D f_1      (10)hyp_v       ",&
               "(11)hyp_z       (12)hyp_xy      (13)zv poisson"!  (14)nonlinearity"
       END DO
    endif

    Allocate(fsa_moments_mat(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2,1:n_fsa_mats,ln1:ln2))

    Call set_mats_diag_fsa_moments

  end subroutine initialize_diag_fsa_moments

  !------------------------------

  Subroutine finalize_diag_fsa_moments
    Integer :: n

    if(write_std) then
       do n=ln1,ln2
          if (mype.eq.0) then
             close(FSAMOMFILE(n))
          endif
       enddo
       deallocate(FSAMOMFILE)
    end if
    deallocate(fsa_moments_mat)
  End Subroutine finalize_diag_fsa_moments
  
  !------------------------------

  SUBROUTINE set_mats_diag_fsa_moments
    Implicit None
    Integer :: k,l,m,n

    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO n=ln1,ln2
                !mat: first index vpar, second vperp (!)
                !mat_00 -- density moment            check normalization
                !fsa_moments_mat(:,:,k,l,m,3,n) = mat_00(:,:,k,l,m)
                !mat_00 -- parallel momentum moment  check normalization
                !fsa_moments_mat(:,:,k,l,m,2,n) = vp(l)*mat_00(:,:,k,l,m)
                !mat_20 + mat_02 -- energy moment    check normalization!!
                fsa_moments_mat(:,:,k,l,m,1,n) = (vp(l)*vp(l) + mu(m)*geom%Bfield(pi1:pi2,:,k)) &
                     & * mat_00(:,:,k,l,m)
             END do
          END DO
       END DO
    END DO    

  END SUBROUTINE set_mats_diag_fsa_moments
  

  subroutine scalar_product(v1,v2,sp)

    implicit none
    complex, intent(in), dimension(li0*lj0*lk0*ll0*lm0*ln0) :: v1,v2
    complex, intent(out) :: sp

    sp=sum(conjg(v1)*v2)

  end subroutine scalar_product

  subroutine get_fsa_moments(rhs_in,withbounds,emfields_in,n_fsa_mats,p_fsamom_terms)
    complex,dimension(:,:,:,:,:,:),intent(in)::rhs_in
    logical,intent(in)::withbounds
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in):: emfields_in
    integer,intent(in)::n_fsa_mats
    real,dimension(li1:li2,ln1:ln2,1:n_fsa_mats),intent(out)::p_fsamom_terms
    
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2,1:n_fsa_mats):: momc
    real, Dimension(li1:li2,ln1:ln2):: fsavg
    integer::o
    momc=cmplx(0,0)

#ifdef WITH_FLR
    Call calc_vsp_moment(n_fsa_mats,rhs_in,withbounds,emfields_in,fsa_moments_mat,momc,.false.)
#else
    !for testing: skip gyroaverage on F
    Call my_calc_vsp_moment(n_fsa_mats,rhs_in,withbounds,fsa_moments_mat,momc)
#endif
    do o=1,n_fsa_mats
       call flux_surface_average(momc(:,:,:,:,o),fsavg)
       Call my_real_sum_vw(fsavg,size(fsavg))
       p_fsamom_terms(:,:,o)=fsavg
    enddo

 end subroutine get_fsa_moments

  !>This subroutine gets the terms of interest in the fsa_moments equation.
  !!\todo Simplify structure for arakawa_zv
  subroutine get_fsamom_terms(g_1_in,fsamom_terms)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: g_1_in
    real,dimension(li1:li2,ln1:ln2,1:n_fsa_terms) :: fsamom_terms

    ! Local variables
    complex, allocatable, dimension(:,:,:,:,:,:) :: h_in,f_in,rhs,g_temp
    complex, allocatable, dimension(:,:,:) :: rhs_block
    complex, allocatable, dimension(:,:,:) :: g_block
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: fields_in
    logical :: calc_this_term(n_fsa_terms)
    integer :: lbg1, lbg2,klmn, iblock, iterm

    lbg0 = lklmn0 / nblocks

    allocate( rhs(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    allocate(f_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    if (arakawa_zv) then
       allocate(h_in(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    endif

    fsamom_terms = 0.0
    calc_this_term(1:3) = .false.
    if (arakawa_zv) then
       call calc_aux_fields(g_1_in,fields_in,f_in,.true.,h_in) !compute f and h
    else
       call calc_aux_fields(g_1_in,fields_in,f_in,.false.) !compute f
    endif

    ! ---------------------------- TOTAL fsa_moments --------------------------
    ! IMPORTANT: Don't change the order!
    ! calling calc_rhs_only will update some pointers (e.g., ptr_bar_emfields)
    ! which are required for subsequent calls to, for instance, the nonlinearity

    ! direct moments of f_ : density, momentum, energy
    ! f_ has boundary points in all dimensions
    call get_fsa_moments(f_in,.true.,fields_in,n_fsa_mats,fsamom_terms(:,:,1))
    
    !Total change of dens,mom,energy

    !need g_temp because some compilers don't like sending an input variable 
    allocate(g_temp(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    g_temp=g_1_in
    rhs=(0.0,0.0)
    if (arakawa_zv) then
       call calc_rhs_only(f_in,g_temp,fields_in,rhs,2,h_in)
    else
       call calc_rhs_only(f_in,g_temp,fields_in,rhs,2)
    endif
    !collisions are taken out of calc_rhs_only in case of operator splitting: add here
    if(coll_split) then
       call equ_collisions(f_in,rhs,replace_rhs=.false.)
    endif
    deallocate(g_temp)
    
    call get_fsa_moments(rhs,.false.,fields_in,n_fsa_mats,fsamom_terms(:,:,2))

    !------------------ INDIVIDUAL TERMS -------------------------------
    !---------------------- nonblocked terms ---------------------------------
    !collisions
    rhs=(0.0,0.0)
    if (collision_op.ne.'none') then
       call equ_collisions(f_in,rhs,.true.)
       call get_fsa_moments(rhs,.false.,fields_in,n_fsa_mats,fsamom_terms(:,:,3))
    endif

PERFON('fsa_block')
    !------------------- strip mining terms ---------------------------
    allocate(rhs_block(li1:li2,lj1:lj2,1:lbg0))
    allocate(g_block(lbi:ubi,lj1:lj2,1:lbg0))
    do iterm=4,n_fsa_terms
       calc_this_term(iterm) = .true.
       rhs = (0.0,0.0)
       Do iblock=1,nblocks
          rhs_block = (0.0,0.0)
          
          lbg1 = (iblock-1)*lbg0+1
          lbg2 = lbg1+lbg0-1
          select case(iterm) 
             ! ---------------- SOURCE/DRIVE TERMS -------------------
          case(4)
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

             if (xy_local) then
                call add_dchidxy_orig(fields_in,ptr_bar_emfields,&
                &rhs_block,ptr_barchi,ptr_dbarchidxy,lbg1,lbg2,.true.)
                !with_curv=T....full dchidxy terms
             else
                call add_dchidxy(chi_block,rhs_block,ptr_dbarchidxy,lbg1,lbg2,.true.)
                !with_curv=T....full dchidxy terms
             end if
          case(5)
          ! krook particle and haet sources and/or profiled heat source
             if (.not.xy_local) then
                call add_f1_sources(rhs_block,lbg1,lbg2,0)
             else
                calc_this_term(iterm) = .false.
             endif
          case(6)
          ! standard explicit_buffer=.false.
             if (.not.x_local.and.explicit_buffer) then 
                call add_kBuffer_explicit(f_in,rhs_block,lbg1,lbg2)
             else
                calc_this_term(iterm) = .false.
             endif
          case(7)
             ! f0_source term
             if (include_f0_contr) then
                call add_f0_term(rhs_block,lbg1,lbg2,time)
             else
                calc_this_term(iterm) = .false.
             endif
          case(8)
             !dgdxz terms  (this is the vD grad f_1 nonlocal term in neoclassics)
             if (.not.xy_local) then
                call get_g_block(g_1_in,g_block,iblock)
                call exchange_x(x_boundary_block,g_block)
                call add_dgdxy(g_block, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
             else
                call add_dgdxy(g_1_in, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
             endif
             ! ---------------- SINK TERMS -------------------
          case(9)
             !hyp_v
             if (hyp_v.gt.0) then
                if (arakawa_zv) then 
                   if (hyp_on_h) then
                      call add_hypv_ak(h_in,rhs_block,lbg1,lbg2)
                   else
                      call add_hypv_ak(f_in,rhs_block,lbg1,lbg2) 
                   endif
                else
                   call add_hypv_block(f_in,rhs_block,lbg1,lbg2)
                endif
             else
                calc_this_term(iterm) = .false.
             endif
          case(10)
             !hyp_z
             if (hyp_z.gt.0) then
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
             else
                calc_this_term(iterm) = .false.
             endif
          case(11)
             !hyp_xy
             if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then
                if (arakawa_zv) then
if (hyp_on_h) then
                    call add_dfdxy(h_in,rhs_block,lbg1,lbg2) !f_in is h
else
                    call add_dfdxy(f_in,rhs_block,lbg1,lbg2) !f_in is h
endif                
                else
                    call add_dfdxy(f_in,rhs_block,lbg1,lbg2) !f_in is f
                endif
             else
                calc_this_term(iterm) = .false.
             endif
             ! ----------------- Analytically vanishing terms -----------------
          case(12)
             ! Poisson bracket in z,vpar, i.e.
             ! (dfdzv+dchidzv) 
             !!this includes the (z,vpar) hyperdiffusion
             if (arakawa_zv) then
                call equ_dzv(h_in,rhs_block,lbg1,lbg2)
                if (hyp_on_h) then
                   if (hypz_compensation) call equ_comp_hypz(fields_in,ptr_bar_emfields,rhs_block,lbg1,lbg2)
                endif
             else
                call equ_dfdzv(f_in,rhs_block,lbg1,lbg2)
                call add_dchidz(fields_in, ptr_bar_emfields, rhs_block, lbg1, lbg2)
             end if
          case(13)
             !nonlinearity
             if (rhs_nl) then
                if ((.not.xy_local).and.(nblocks.gt.1)) then
                   call add_dgdxy(g_1_in, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
                   rhs_block = 0.0
                endif
                CALL this_nonlinear_term%add(g_1_in, ptr_dgdxy, fields_in,&
                     &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, 2)
             else
                calc_this_term(iterm) = .false.
             endif
          case default
             stop 'wrong index in blockterm loop'
          end select

          if (calc_this_term(iterm)) &
               call add_rhs_block(rhs,rhs_block,iblock)

       enddo  !close block loop
       
       !rhs for iterm is completed. we can use the fsa_moments routine
       if (calc_this_term(iterm)) &
            call get_fsa_moments(rhs,.false.,fields_in,n_fsa_mats,fsamom_terms(:,:,iterm))

    enddo !close iterm loop

    deallocate(rhs_block)
    deallocate(g_block)
    deallocate(rhs)
    if (arakawa_zv) deallocate(h_in)

PERFOFF


  end subroutine get_fsamom_terms


  subroutine add_rhs_block(rhs_inout,rhs_block_in,iblock)
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),intent(inout) :: rhs_inout
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0),intent(in) :: rhs_block_in
    integer, intent(in) :: iblock

    rhs_inout(:,:,:,iblock)=rhs_block_in

  end subroutine add_rhs_block
  
  subroutine get_g_block(g_1_in,g_1_block_out,iblock)
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),intent(in) :: g_1_in
    complex, dimension(lbi:ubi, lj1:lj2, 1:lbg0),intent(out) :: g_1_block_out
    integer, intent(in) :: iblock

    g_1_block_out(li1:li2,lj1:lj2,1:lbg0)=g_1_in(li1:li2,lj1:lj2,1:lbg0,iblock)

  end subroutine get_g_block


  !>Driver for fsa_moments diagnostics 
  !!
  subroutine diag_fsa_moments
    real,dimension(li1:li2,ln1:ln2,1:n_fsa_terms) :: fsamom_terms

PERFON('dia_en')

    call get_fsamom_terms(g_1,fsamom_terms)
    if(write_std) then
       call fsa_moments_write(fsamom_terms)
    end if


PERFOFF
  end subroutine diag_fsa_moments

  !>Gathers an input 3D array and outputs as unformatted data to file determined by filehandle.
  !!
  !> todo possibly switch to binary output instead of gnuplot compatible!!
  subroutine fsa_moments_write(fsamom_terms)
    !****
    implicit none
    real,dimension(li1:li2,ln1:ln2,1:n_fsa_terms),intent(in) :: fsamom_terms
    Real, Dimension(0:nx0-1,1:n_fsa_terms) :: out_arr

    integer::i, n, o, ierr

    ! Local variables
       
    do n=ln1,ln2
       if ((my_pey+my_pez+my_pev+my_pew).eq.0) then 
          do o=1,n_fsa_terms
             Call mpi_gather(fsamom_terms(li1,n,o), li0, MPI_REAL_TYPE,&
                  out_arr(0,o), li0, MPI_REAL_TYPE,&
                  0, mpi_comm_x, ierr)
          enddo

          IF (write_prof_pe) THEN
             
             !for now, the output is GNUPLOT compatible
             !in future, we might change to binary or HDF5 output
             WRITE(FSAMOMFILE(n), "(A,F14.6,1X,I5)") '#',time, block_nr
             
             do i=0,nx0-1
                write(FSAMOMFILE(n), "(13ES16.6)") xval_a(i), out_arr(i,:)
             enddo
             
             WRITE(FSAMOMFILE(n), "(A)") ''
             WRITE(FSAMOMFILE(n), "(A)") ''

             block_nr = block_nr + 1

            call flush(FSAMOMFILE(n))
          ENDIF
       endif
    enddo

  end subroutine fsa_moments_write


  SUBROUTINE my_calc_vsp_moment(n_mom,p_dist,withbounds,p_mat,p_mom)
    Implicit none
    INTEGER, INTENT(in):: n_mom
    COMPLEX, DIMENSION(:,:,:,:,:,:), INTENT(in) :: p_dist
    logical, intent(in)::withbounds
    REAL, DIMENSION(pi1:pi2,pj1:pj2,lk1:lk2,1:n_mom,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(in):: p_mat 
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_mom,ln1:ln2),INTENT(out):: p_mom
    Complex, DIMENSION(li1:li2) :: pullback_dist_func

    Integer :: j,k,l,m,n,o
    Integer :: ioff,joff,koff,loff,moff,noff
    
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

    Do n=ln1,ln2
       DO m=lm1,lm2
          Do l=ll1,ll2
             Do k=lk1,lk2
                Do j=lj1,lj2
                   !!skipGyroaverage of F1+(q <phi> + T(x0) mu <Bpar>) F0/T(x)
                   pullback_dist_func(li1:li2) = p_dist(li1+ioff:li2+ioff,j+joff,k+koff,l+loff,m+moff,n+noff)
                   do o=1,n_mom
                      p_mom(li1:li2,j,k,o,n) = p_mom(li1:li2,j,k,o,n)+pullback_dist_func(li1:li2)*&
                           &p_mat(li1:li2,pj1,k,o,l,m,n)
                   enddo
                End Do
             End Do
          End Do
       End Do
    End Do

  END SUBROUTINE my_calc_vsp_moment

#endif


End Module diagnostics_fsa_moments
