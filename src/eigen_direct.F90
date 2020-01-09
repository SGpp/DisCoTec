#include "redef.h"
!>Direct eigenvalue solver to compute all eigenvalues/vectors of the linear gyrokinetic operator
!!
!!The computation is fully parallel and based on the Mario Sempf's wrappers for the SCALAPACK 
!!library. Note that direct eigenvalue solvers scale like the dimensionality of the matrix cubed,
!!so that this becomes very slow very quickly. Furthermore, the memory requirements are quite 
!!significant
module eigen_direct

  use eigen_parameters
  use file_io, only: get_unit_nr
  use fullmatrix_aux
  use communications
  use diagnostics
  use diagnostics_energy, only: get_energy_total
  use scalapack_aux
  use par_mod, only: time_ev
  implicit none   
  public:: ev_direct
  private
  
  complex,dimension(:,:),allocatable:: sl_mat

contains
  !>Computes all eigenvalues and optionally left/right eigenvectors and writes them to
  !!eigenvalues.dat / analyzes them
  subroutine ev_direct
    double precision :: time_1, time_2
    complex,dimension(vlen*n_procs_sim):: eigval
    real,dimension(vlen*n_procs_sim):: growth_rates
    integer,dimension(vlen*n_procs_sim):: eval_sort
    complex,dimension(:,:),allocatable:: bl_ev_left, bl_ev_right,bl_q
    integer:: evfile, ind, ierr
    integer,dimension(1):: mloc
    integer :: ev_handle
    real,allocatable,dimension(:) :: scale_ev
    logical :: dummy

    PERFON('ev_direct')

    vlen_gl=vlen*n_procs_sim
    allocate(scale_ev(0:vlen_gl-1))

    !generate explicit representation of the linear operator
    call sl_init(sl_con,n_procs_sim,1)
    call descinit(sl_desc,vlen_gl,vlen_gl,vlen,vlen_gl,0,0,sl_con,vlen,ierr)  
    allocate(sl_mat(0:vlen-1,0:vlen_gl-1))
    call get_expl_sl_mat(sl_mat)

    if(mype.eq.0) then
       write(*,"(a)") '******************************************************'
       write(*,"(a)") 'Starting eigenvalue computation with SCALAPACK'
       write(*,"(a)") '******************************************************'
    endif

    !transform linear operator to 2D block cyclic representation
    call initialize_bl
    call sl_2_bl(sl_mat,bl_mat)

    allocate(bl_ev_left(loc_rows,loc_cols),bl_ev_right(loc_rows,loc_cols),&
         bl_q(loc_rows,loc_cols))

    !initialize bl_q with unit matrix
    sl_mat=0.
    do ind=0,vlen-1
       sl_mat(ind,ind+mype*vlen)=1.
    enddo
    call sl_2_bl(sl_mat,bl_q)

    call get_systime(time_1)
    call c_eigensolve(ev_left,ev_right,.false.,vlen_gl,bl_mat, &
         bl_desc,eigval,bl_ev_left,bl_desc,bl_ev_right, &
         bl_desc,bl_q,bl_desc)

    !output to file and standard out
    call get_unit_nr(evfile)
    if(vec_out_limit) then 
       growth_rates=real(eigval)
       do ind=1,vlen_gl
          mloc=maxloc(growth_rates(1:vlen_gl))
          eval_sort(ind)=mloc(1)
          growth_rates(mloc(1))=-1.0e10 
       enddo
       call get_unit_nr(ev_handle)
       open(unit=ev_handle,file=trim(diagdir)//'/ev_sort.dat',status='unknown')
       do ind=1,vlen_gl
          write(ev_handle,*) ind,eval_sort(ind)
       end do
       close(ev_handle)
    end if

    if (mype.eq.0) then
       open(evfile, file=trim(diagdir)//'/eigenvalues'//trim(file_extension),&
            &form='formatted', status='replace')             
       write(*,"(a)") 'computed eigenvalues: '
       write(evfile, "(a,F7.3)") '#eigenvalues (re/im part) at ky=', ky(lg1)
       if(vec_out_limit) then 
          do ind=1,n_ev
             write(evfile, "(8ES20.8)") eigval(eval_sort(ind))
          enddo
       else !vec_out_limit
          do ind=1,vlen_gl
             write(evfile, "(8ES20.8)") eigval(ind)
          enddo
       end if !vec_out_limit
       close(evfile)
    endif

    call get_systime(time_2)
    if (mype.eq.0) write(*,"(A,F8.2)") 'time for eigenvalue computation: ', time_2-time_1
    time_ev=time_2-time_1

    call initialize_all_diags

    if (ev_left.or.ev_right) then
       PERFON('ev_out')
       call get_systime(time_1)
       
       !output right eigenvectors
       if (ev_right) then  
          if(ev_out) call ev_out_open(ev_handle,1)
          call bl_2_sl(bl_ev_right,sl_mat)  
         if(vec_out_limit) then
          do ind=1,n_ev
             call ccopy(vlen,sl_mat(0,eval_sort(ind)-1),1,g_1(li1,lj1,lk1,ll1,lm1,ln1),1)
             time=ind
             call get_energy_total(g_1,scale_ev(1))
             g_1 = g_1 / SQRT(scale_ev(1))
             call exec_all_diags(0,time,dummy,dummy,dummy)
             if(ev_out) call ev_out_write(g_1,ind,vlen_gl,ev_handle)
          enddo

         else !vec_out_limit
          do ind=0,vlen_gl-1
             call ccopy(vlen,sl_mat(0,ind),1,g_1(li1,lj1,lk1,ll1,lm1,ln1),1)
             time=ind
             call get_energy_total(g_1,scale_ev(ind))
             g_1 = g_1 / SQRT(scale_ev(ind))
             call exec_all_diags(0,time,dummy,dummy,dummy)
             if(ev_out) call ev_out_write(g_1,ind+1,vlen_gl,ev_handle)
          enddo
         end if !vec_out_limit
          if(ev_out) call ev_out_close(ev_handle)
          !use most unstable mode as initial value for testing
          mloc=maxloc(real(eigval))-1
          if(mype.eq.0) print*,'maximum growth rate (check): ',eigval(mloc(1)+1) 
!          call ccopy(vlen,sl_mat(0,mloc(1)),1,g_1(li1,lj1,lk1,ll1,lm1,ln1),1)
       else
          mloc=maxloc(real(eigval))-1
          if(mype.eq.0) print*,'maximum growth rate (check): ',eigval(mloc+1) 
       endif

       !output left eigenvectors (reuses sl_mat)       
       if (ev_left) then 
          call ev_out_open(ev_handle,2)
          call bl_2_sl(bl_ev_left,sl_mat)

         if(vec_out_limit) then
          vec_out_limit=.false.  !ev_out_write uses this as a flag
          do ind=1,n_ev
             call ccopy(vlen,sl_mat(0,eval_sort(ind)-1),1,g_1(li1,lj1,lk1,ll1,lm1,ln1),1)
             time=ind
!             call exec_all_diags(0,time,dummy,dummy,dummy)
             call get_energy_total(g_1,scale_ev(1))
             g_1 = g_1 / SQRT(scale_ev(1))
             if(ev_out) call ev_out_write(g_1,ind,vlen_gl,ev_handle)
          enddo
         else ! vec_out_limit
          do ind=0,vlen_gl-1
             call ccopy(vlen,sl_mat(0,ind),1,g_1(li1,lj1,lk1,ll1,lm1,ln1),1)
             time=ind
             !rescale left eigenvectors to preserve orthonormality 
             if(ev_right) g_1=g_1*sqrt(scale_ev(ind))
             call ev_out_write(g_1,ind+1,vlen_gl,ev_handle)
          enddo
         end if
          call ev_out_close(ev_handle)
       endif

       PERFOFF
       call get_systime(time_2)
       if (mype.eq.0) write(*,"(A,F8.2)") 'time for eigenvector output: ', time_2-time_1
    endif

    call finalize_all_diags
    
    deallocate(sl_mat,bl_mat,bl_ev_left,bl_ev_right)

    call blacs_gridexit(bl_con)
    call blacs_gridexit(sl_con)

    PERFOFF

  end subroutine ev_direct

!include Mario's library
#include "mpl_routines.F90"

end module eigen_direct

