#include "redef.h"
#include "intrinsic_sizes.h"
  
module convergence_monitoring
  use discretization
  use coordinates
  use communications
  use par_in
  use par_other
  use geometry
  use aux_fields
  use RK_coefficients, only: implicit_scheme
  use file_io
  use quasilinear_model

#ifdef WITHFUTILS
  use futils
#endif

  implicit none

  public:: istep_omega, omega_prec, convergence_exit
  public:: set_cm_defaults, mem_est_convergence_monitoring
  public:: initialize_diag_convergence, diag_convergence, finalize_diag_convergence
  public:: set_phi_old, rescale_phi_old, check_convergence_monitoring
  public:: om_conv

  private

  logical:: convergence_exit, use_gamma=.false.
  real:: omega_prec, ga_av


  integer:: istep_omega=-1, OMEGAFILE, GAMMAFILE
  integer:: lie, lje
  integer::ndigits
  complex:: om_conv
  complex, dimension(:,:,:), allocatable:: phi_old

  !diag_gamma
  real:: oldval = 1.
  integer:: count=0, galen=20
  real,dimension(:),allocatable:: gamma
  real, dimension(:,:,:),allocatable:: phisq
  logical:: loc_del_phi

contains

  subroutine check_convergence_monitoring
    !switch off omega diagnostic for nonlinear runs, GAM investigations (ky=0)
    !or nky0>1 as the lx-adapation per ky mode has been removed some time ago
    if ((nonlinear).or.((nky0.eq.1).and.((ky0_ind.eq.0).or. &
         &((abs(kymin).lt.1e-5).and.(n0_global.eq.0.or.n0_global.eq.-1111))))&
         &.or.((nky0.gt.1).and.y_local)) istep_omega=0

    if (istep_omega .eq. -1) istep_omega = 20

    ! diag_omega will only converge if there is just 1 independent mode per ky.
    ! For x_local=.f. this is always true, for x_local=.true. only with adapt_lx, 
    ! for y_local=.f. never.
    ! For these cases, the simpler diag_gamma is used instead
    if ((.not.y_local).or.((x_local).and.(.not.adapt_lx).and.(nx0.ne.1))) use_gamma=.true.

    if ((ExB).and.(istep_omega.gt.0)) then
       if (mype==0) &
            write(*,'(A)') "WARNING: Omega diagnostic deactivated due to incompatibility with ExB shear flow model"
       istep_omega = 0    
    endif

    !determine, if we have phi=0
    if (only_Er .and. only_zonal) then
       !write(*,"(A)") "for ky=0 computations, only_Er to phi=<phi>_fs=0"
       loc_del_phi = .true.
    endif
    if (del_fields) loc_del_phi = .true.

    if (only_neo.and.istep_omega.gt.0) then
       istep_omega = 0
       if (mype==0) write(*,"(A)") "switching off omega diagnostics for pure neoclassical run"
    elseif ((loc_del_phi.or.del_phi).and.(istep_omega.gt.0)) then
       if (mype==0) &
            write(*,'(A)') "WARNING: automatic convergence monitoring deactivated for del_phi mode"
       istep_omega = 0    
    endif

  end subroutine check_convergence_monitoring

  subroutine set_cm_defaults

    omega_prec=1e-3
    convergence_exit = .false.
    !use_gamma=.false.
    loc_del_phi=.false.

  end subroutine set_cm_defaults

  !>Give an estimate of the memory requirements of diag_omega
  Function mem_est_convergence_monitoring(mem_req_in)
    real:: mem_est_convergence_monitoring, mem_req_in
    real:: mem_loc=0

    !phi_old
    mem_loc = lijk0*SIZE_OF_COMPLEX_MB

    mem_est_convergence_monitoring = mem_loc+mem_req_in

  End Function mem_est_convergence_monitoring

  subroutine initialize_diag_convergence

    !formatted output precision
    ndigits= int(log10(1.0/omega_prec))
    !if (ndigits<4)ndigits=4

    if (use_gamma) then
       call initialize_diag_gamma
    else
       call initialize_diag_omega
    end if

  end subroutine initialize_diag_convergence

  subroutine diag_convergence
    
    if (use_gamma) then
       call diag_gamma
    else
       call diag_omega
    end if

  end subroutine diag_convergence


  subroutine finalize_diag_convergence

    if (use_gamma) then
       call finalize_diag_gamma
    else
       call finalize_diag_omega
    end if


  end subroutine finalize_diag_convergence


!!!******************************************************************!!!

  SUBROUTINE initialize_diag_omega

    if (xy_local) then
       if (yx_order) then
          lie = li2
          lje = hkx !take only pos. modes
       else
          lie = hkx !take only pos. modes
          lje = lj2
       end if
    else
       lie = li2 !take all x and y
       lje = lj2
    endif
    ALLOCATE(phi_old(li1:lie,lj1:lje,lk1:lk2))
    om_conv = 0.0

    call get_unit_nr(OMEGAFILE)   
    If(mype==0) OPEN(OMEGAFILE, file=trim(diagdir)//'/omega'//&
         &trim(file_extension), &
         &form='formatted', status='replace', position='rewind')
  END SUBROUTINE initialize_diag_omega


  !> Omega diagnostic: Calculates growth rate and real frequency
  SUBROUTINE diag_omega
    real :: stddev_loc2,stddev_glob2
    real :: weight_loc,weight_glob
    complex, dimension(li1:lie,lj1:lje,lk1:lk2) :: om_arr
    complex :: om_av, fac_loc
    integer:: ierr
    character(len=32)  :: fmt
#ifdef COMBI_MGR
    complex :: om_out
    real :: tmp1, tmp2
#endif

    !PERFON('diag_omega')

    if (om_conv.eq.(0.0,0.0)) then
       om_arr=log(emfields(li1:lie,lj1:lje,lk1:lk2,1)/phi_old)/(dt*istep_omega)
       weight_loc=sum(abs(phi_old))
       Call mpi_allreduce(weight_loc, weight_glob, 1,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_xyz, ierr)
       fac_loc=sum(om_arr*abs(phi_old))
       Call mpi_allreduce(fac_loc, om_av, 1,&
            MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_xyz, ierr)
       om_av = om_av/weight_glob
       stddev_loc2=sum(abs(om_arr-om_av)**2*abs(phi_old))
       Call mpi_allreduce(stddev_loc2,stddev_glob2, 1,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_xyz, ierr)

       ! calculate 99% confidence interval (factor*std_dev/sqrt(n)) 
       ! assuming normal distribution of omega(kx, z)
       if (stddev_glob2*6.7/((nz0*n_procs_x*(lie-li1+1)*(lje-lj1+1)-1)*&
            weight_glob).lt.(0.1*omega_prec)**2) then
          if (implicit_scheme) om_av=(1-exp(-om_av*dt))/dt
          om_conv = om_av
          if (mype/n_procs_y.eq.0) then
             write(fmt,"(a,4(i2.2,a))") &
                  "(a,f8.5,a,f",ndigits+4,".",ndigits,&
                  ",a,f",ndigits+4,".",ndigits,")"
             write(*,fmt) 'ky=', ky,' gamma = ', real(om_av),&
                  ' omega= ',aimag(om_av)
             if(implicit_scheme) write (*,"(a)") &
                  &'(error induced by implicit Euler has been corrected)'
          endif
       endif
    endif

#ifdef COMBI_MGR
    if (implicit_scheme) om_av=(1-exp(-om_av*dt))/dt
    om_out = om_av
    if (mype==0) then
        tmp1 = real(om_out)
        tmp2 = aimag(om_out)
        call write_omega( itime, tmp1, tmp2 )
    endif
#endif

    convergence_exit = om_conv.ne.(0.0,0.0)

    phi_old=emfields(li1:lie,lj1:lje,lk1:lk2,1)

    !PERFOFF

  END SUBROUTINE diag_omega

  subroutine reset_om_diag
    om_conv=0.0
  end subroutine reset_om_diag

  subroutine set_phi_old
    if (.not.use_gamma) phi_old=emfields(li1:lie,lj1:lje,lk1:lk2,1)
  end subroutine set_phi_old

  subroutine rescale_phi_old(globsum)
    real:: globsum

    if (use_gamma) then
       oldval=oldval/globsum/globsum
    else
       phi_old=phi_old/globsum
    end if
  end subroutine rescale_phi_old


  Subroutine finalize_diag_omega
    character(len=28)  :: fmt
#ifdef WITHFUTILS
    integer :: fidomega_h5
    integer :: rank
    integer, dimension(2) :: dims
#endif
    if (mype==0) then
       write(fmt,"(a,2(i2.2,a))") "(f7.3,2(1X,f",ndigits+6,".",ndigits,"))"
       write(OMEGAFILE,fmt) kymin,real(om_conv), aimag(om_conv)
       close(OMEGAFILE)
    endif

    if (quasilin_model.gt.0) call set_quasilin(real(om_conv))

    DEALLOCATE(phi_old)

#ifdef WITHFUTILS
    if(write_h5) then
       if((mype).eq.0) then
          ! OMEGA - GAMMA
          call creatf(trim(diagdir)//'/omega'//trim(file_extension)//'.h5', &
               fidomega_h5, "omega-gamma diagnostics", 'd')
          rank = 0
          call creatd(fidomega_h5, rank, dims, "/ky", "ky")
          call creatd(fidomega_h5, rank, dims, "/gamma", "gamma")
          call creatd(fidomega_h5, rank, dims, "/omega", "omega")

          call append(fidomega_h5, '/ky', real(ky,8))
          call append(fidomega_h5, '/gamma', real(om_conv,8))
          call append(fidomega_h5, '/omega', real(aimag(om_conv),8))

          call closef(fidomega_h5)
       end if
    end if
#endif

  End Subroutine finalize_diag_omega


  subroutine initialize_diag_gamma

    allocate(phisq(li1:li2,lj1:lj2,lk1:lk2))
    allocate(gamma(0:galen-1))
    gamma=0.
    count=0
    call get_unit_nr(GAMMAFILE)   
    If(mype==0) OPEN(GAMMAFILE, file=trim(diagdir)//'/omega'//&
         &trim(file_extension), &
         &form='formatted', status='replace', position='rewind')
  end subroutine initialize_diag_gamma

  subroutine diag_gamma
    real, dimension(li1:li2,lj1:lj2,lk1:lk2):: phisq
    real:: newval_loc, newval_glob, stddev2
    integer:: i, j, k, ierr, pii
    character(len=30)  :: fmt

    ! take into account negative kj, except for 0-mode
    phisq = 2. * real(emfields(:,:,lk1:lk2,1) * conjg(emfields(:,:,lk1:lk2,1)))
    if (p_has_0_mode)  then
       if (xy_local.and.yx_order) then
          phisq(li1,:,:) = 0.5 * phisq(li1,:,:)
       else
          phisq(:,lj1,:) = 0.5 * phisq(:,lj1,:)
       end if
    end if
    
    newval_loc=0.
    do k=lk1,lk2
       do j=lj1,lj2
          pii=li1
          do i=li1,li2
             if (pmi0.gt.1) pii=i
             newval_loc=newval_loc+phisq(i,j,k)*geom%jacobian(pii,pj1,k)
          enddo
       end do
    end do
    Call mpi_allreduce(newval_loc, newval_glob, 1,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_xyz, ierr)

    gamma(mod(count,galen))=0.5*log(newval_glob/oldval)/(istep_omega*dt)
    ga_av=sum(gamma)/galen
    stddev2=sum(abs(gamma-ga_av)**2)

    ! calculate 99% confidence interval 
    ! assuming normal distribution
    if ((stddev2*6.7/galen.lt.(0.1*omega_prec)**2).and.(count.ge.galen)) then
       om_conv = ga_av
       if (mype/n_procs_y.eq.0) then
          write(fmt,"(a,2(i2.2,a))") &
               "(a,f",ndigits+4,".",ndigits,")"
          write(*,fmt) 'growth rate converged: gamma = ', ga_av
       endif
       convergence_exit=.true.
    endif

    !remember value
    oldval=newval_glob
    count=count+1

  end subroutine diag_gamma

  subroutine finalize_diag_gamma
    character(len=28)  :: fmt

    if (mype==0) then
       write(fmt,"(a,4(i2.2,a))") &
            "(f7.3,f",ndigits+4,".",ndigits,",f",ndigits+4,".",ndigits,")"
       write(GAMMAFILE,fmt) kymin,ga_av,0.
       close(GAMMAFILE)
    endif

    if (quasilin_model.gt.0) call set_quasilin(ga_av)

    deallocate(phisq,gamma)
  end subroutine finalize_diag_gamma

end module convergence_monitoring
