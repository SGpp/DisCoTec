#include "redef.h"
#ifdef WITHSLEPC
#include "petscversion.h"
#endif
!>Parameters for eigenvalue computations
module eigen_parameters
  
  use discretization, only: mype
  use par_in, only: comp_type, nonlinear
  use coordinates, only: kymax, ky

  implicit none

  public

  character(len=20) :: which_ev='none'

  complex:: ev_shift=(10.,0.)
#ifdef COMBI
  complex:: ev_shift_ref=(0.0,0.0)
#endif
  
  logical:: ev_left=.false., ev_right=.true.
  
  !parameters for slepc
  integer:: n_ev
  !the following params are set to approp. defaults in slepc_aux.F90 if <0
  integer:: ev_max_it=-1, ev_n_test=-1
  integer:: ev_n_test_sr=16,n_ev_sr = 1, ev_max_it_sr=12  !for extreme smallest_real eigenvalues
  real:: ev_prec=1e-4
  integer:: ksp_max_it=0
  integer:: it_ev=-1

contains

  subroutine check_eigenvalue_comp
    logical:: write_pe=.false.

    if (mype.eq.0) write_pe=.true.

    !ky>1 modes can exhibit frequencies being significantly larger 
    !than 10 (for instance, ETG) and hence ev_shift is adapted
    if (kymax.gt.1) ev_shift=cmplx(max(10.,2.*kymax),0.0)

    if (comp_type.eq.'EV') then
       if (nonlinear) then
          if (write_pe) write(*,"(A)") "nonlinear eigenvalue computations not possible"
          stop
       endif
       if(which_ev.eq.'none') then
#if (PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR==0)
          if (write_pe) write(*,"(A)") "no which_ev given in parameters, using harmonic"
          which_ev='harmonic'
#else
          if (write_pe) write(*,"(A)") "no which_ev given in parameters, using jd"
          which_ev='jd'
#endif         
       end if
       if (n_ev.eq.0) n_ev=1
    endif

    select case (which_ev)
    case('jd')
       if(ev_n_test.eq.0) ev_n_test=max(2*n_ev,n_ev+30)
    case('gd')
       if(ev_n_test.eq.0) ev_n_test=max(2*n_ev,n_ev+25)
    case('shift_invert')
       if(ev_n_test.eq.0) ev_n_test=max(2*n_ev,n_ev+15)
    case('harmonic')
       if(ev_n_test.eq.0) ev_n_test=max(2*n_ev,n_ev+25)
    case default
       if(ev_n_test.eq.0) ev_n_test=max(2*n_ev,n_ev+15)
       !IE1p
    end select
    
  end subroutine check_eigenvalue_comp

end module eigen_parameters
