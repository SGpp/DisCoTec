#include "redef.h"
#include "intrinsic_sizes.h"

module time_averages
  use discretization
  use communications
  use par_mod
  use spatial_averages
  use file_io

  implicit none

  public:: set_avgfluxes, avgfluxes, avgflux_stime, reset_avgflux
  public:: initialize_avgfluxes, finalize_avgfluxes
  public:: avgflux_type, write_flux_final
  public:: avgflux_alpha

  private

  Real:: avgflux_stime=-1.0, old_avgflux_time = -1.0
  Integer :: write_flux_final = 0 !0: don't write, 1: write w/o mod., 2:norm to ky |phi|^2
  Integer :: avgflux_type = 0 !0: upper sum, 1: exp. avg.
  Real :: avgflux_alpha = 0.008 !decay rate for exp. averaging, default: 1/125
  REAL,DIMENSION(:,:),ALLOCATABLE :: avgfluxes

contains

  !>Give an estimate of the memory requirements
  Real Function mem_est_avgfluxes(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !avgfluxes
    mem_loc = mem_loc + 2.*SIZE_OF_COMPLEX_MB*ln0*2

    mem_est_avgfluxes = mem_req_in + mem_loc
  end Function mem_est_avgfluxes

  subroutine initialize_avgfluxes
    IF (.NOT.ALLOCATED(avgfluxes)) ALLOCATE(avgfluxes(0:n_spec-1,2))

    avgfluxes=0.0
    old_avgflux_time = -1

  end subroutine initialize_avgfluxes

  !>compute simple linear time average (upper sum)
  !!of nrg fluxes or return last time step values (linear runs)
  !! \param var the var array used in diag_nrg
  !! \param emfields the em. fields for optional normalization
  !! \param time the current time step
  subroutine set_avgfluxes(var, emfields, time,nrgcols)
    integer, intent(in):: nrgcols
    Real, Dimension(nrgcols,ln1:ln2),intent(in) :: var
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in)::emfields
    real, intent(in) :: time

    real :: dt_, normval, fac
    integer :: n, o

    if ((avgflux_stime.ge.0.0).and.(avgflux_stime.le.time)) then  

       select case (write_flux_final)
       case(2)
          !take ky |phi|^2 as normalization
          call volume_average_abssq(emfields(:,:,:,1),.true.,normval)
          normval = normval * kymin
       case default
          normval = 1.0
       end select

       fac=1./normval

       if (nonlinear) then
          if (old_avgflux_time.lt.0) then  !first call to this routine
             old_avgflux_time = time
             avgflux_stime = time !set avgflux_stime to actually used starting time

             do n=ln1,ln2
                avgfluxes(n,1) = sum(var(5:6,n))*fac
                avgfluxes(n,2) = sum(var(7:8,n))*fac
             enddo
          else
             dt_ = time-old_avgflux_time
             fac = dt_*fac

             select case (avgflux_type)
             case(0) !simple upper sum integration

                if (old_avgflux_time.eq.avgflux_stime) then
                   !reset avgfluxes
                   avgfluxes = 0.0
                endif
                do n=ln1,ln2
                   avgfluxes(n,1) = avgfluxes(n,1)+sum(var(5:6,n))*fac
                   avgfluxes(n,2) = avgfluxes(n,2)+sum(var(7:8,n))*fac
                enddo
             case(1) !exp. average
                fac = avgflux_alpha*fac
                do n=ln1,ln2
                   avgfluxes(n,1) = (1-fac)*avgfluxes(n,1)+sum(var(5:6,n))*fac
                   avgfluxes(n,2) = (1-fac)*avgfluxes(n,2)+sum(var(7:8,n))*fac
                enddo
             case default
                stop 'unknown avgflux_type'
             end select
             old_avgflux_time = time
          endif
       else !linear (update to latest time step)
          do n=ln1,ln2
             avgfluxes(n,1) = sum(var(5:6,n))*fac
             avgfluxes(n,2) = sum(var(7:8,n))*fac
          enddo
       endif

       do o=1,2
          call my_real_gather_to_0(avgfluxes(:,o),ln1+1,ln0,n_spec,&
               & mpi_comm_spec)
       enddo
       
    endif

  end subroutine set_avgfluxes
  
  Subroutine reset_avgflux

    old_avgflux_time = -1.0
    avgfluxes = 0.0

  End Subroutine reset_avgflux

  subroutine finalize_avgfluxes(time,istep_nrg)
    real, intent(in) :: time
    integer, intent(in) :: istep_nrg 

    integer :: n, ierr, FLUXFINALFILE
    real :: delta_T

    if (avgflux_stime.ge.0.0) then 
       if (nonlinear.and.(avgflux_type.ne.1)) then
          delta_T = old_avgflux_time-avgflux_stime
          if (delta_T.gt.0) avgfluxes=avgfluxes/delta_T
       endif

       if (mype==0) then 
          Write (*,*)
          if (my_sim.ge.0) then
             Write (*,'(A,I4)') 'Time averaged fluxes, my_sim = ', my_sim
          else
             Write (*,'(A)') 'Time averaged fluxes'
          endif
          Write (*,'(A,F6.2,A,I3,A)') '(starting from t=',avgflux_stime,&
               '; every ',istep_nrg,' time steps)'
          Write (*,'(A)') 'Species   | particle flux | heat flux'
          
          DO n=0, n_spec-1
             WRITE (*,"(2A,ES12.4,A,ES12.4)") spec(n)%name, '| ',&
                  &avgfluxes(n,1), '  | ',avgfluxes(n,2)
          END DO

          if (write_flux_final.gt.0) then
             call get_unit_nr(FLUXFINALFILE)   
             OPEN(FLUXFINALFILE, file=trim(diagdir)//'/fluxfinal'//&
                  &trim(file_extension), &
                  &form='formatted', status='replace', position='rewind')
             DO n=0, n_spec-1
                WRITE (FLUXFINALFILE,"(2(ES12.4))") &
                     &avgfluxes(n,1), avgfluxes(n,2)
             END DO
             CLOSE(FLUXFINALFILE)
          endif
       endif
       
       call my_barrier
       
       if (mype==0) then
          Write (*,*)
       endif
    endif
    
    if (par_in_dir.ne.'skip_parfile') then
       DEALLOCATE(avgfluxes)
    else
       call mpi_bcast(avgfluxes,SIZE(avgfluxes),MPI_REAL_TYPE,&
            &0,MY_MPI_COMM_WORLD,ierr)
    endif
  end subroutine finalize_avgfluxes


end module time_averages
