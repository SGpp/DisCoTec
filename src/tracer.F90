!> Interface module for field line tracing of EFIT equilibria. This module borrows heavily from
!! the TRACER code, courtesy of P. Xanthopoulos. The underlying algorithm is explained in 
!! PoP 13, 092301 (2006) and PoP 16, 082303 (2009).
module tracer_mod
  use discretization
  use lagrange_interpolation
  use tracer_aux
  use tracer_rk5
  use tracer_IO
  use par_other, only: pi,print_ini_msg
  use par_in, only: lilo
  use par_geom, only: geomtype
  use communications

  implicit none

  public:: init_tracer, get_tracer, finalize_tracer, get_efit_refvals, tracer_initialized
  private 
  real:: Lref, Bref, x0, rmag, zmag, psi_out, relpsi_out
  real, dimension(:), allocatable:: xval_a
  !for edge_opt, we generate a well-resolved grid and interpolate down to the actual resolution
  integer, parameter:: edge_opt_res=1024
  logical:: tracer_initialized=.false., write_pe
contains

  !> Setup and tracing of EFIT equilibria.
  subroutine init_tracer(geomdir, geomfile, x0_in, lx, rhostar, minor_r, major_R, edge_opt, trpeps)
    character(len=100), intent(in) :: geomdir
    character(len=80), intent(in)  :: geomfile
    real, intent(in):: x0_in, lx, rhostar, edge_opt
    real, intent(out):: major_R, minor_r, trpeps
    real:: dx, c11_int, r, p, zsym, hn1, err, dphi, zmin, zmax
    real:: q, last_q, dqdx, qefit, shat_efit, distance
    real:: z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
               &dBpdz,dBzdr,dBzdp,dBzdz, rhorz
    real:: psirz, relpsirz, psif, relpsif
    integer, parameter :: neqn=8 !Number of ODEs
    real, dimension(pi1gl:pi2gl):: q_prof, dqdx_prof, acc
    real, dimension(pi1gl:pi2gl,2):: rzstart
    real, dimension(neqn):: y 
    integer:: ierr,i,k,kk, n_runs, resolution, resol_b, resol_f, n_iterations, tres
    integer:: i1,i2
    write_pe = (mype.eq.0)

    n_iterations=10 !floor(10-0.75*log(nx0))
    allocate(xval_a(pi1gl:pi2gl))
    x0=x0_in

    !define radial grid (at this point, the coordinates module has not yet been initialized)
    if (x_local) then
       xval_a(pi1)=x0_in
    else
       !at the moment, for lilo we trace the same flux surface nx0 times
       if (lilo) then
          xval_a=x0_in
       else
          dx=lx*rhostar/(px0-1)
          do i=pi1gl,pi2gl
             xval_a(i)=x0_in-lx*rhostar/2+i*dx
          enddo
       endif
    endif

    if (edge_opt.ne.0.) then
       tres=edge_opt_res
    else
       tres=nz0
    endif

    if (write_pe) then
       write(*,*)
       write(*,"(A)") 'Tracing EFIT geometry...'
    endif

    call initialize_efit(geomdir,geomfile, xval_a, rzstart, &
         &Lref, Bref, zmin, zmax, rmag, zmag)

    if (write_pe) write(*,"(A,F6.3,A,F6.3,A)") 'Setting Lref =',Lref,&
         &'m and Bref =',Bref,'T.'
    
    !for lilo, we only trace one flux surface and copy the result
    if (lilo) then
       i1=pi1gl; i2=pi1gl
    else
       i1=pi1gl; i2=pi2gl
    endif

    do i=i1,i2
       n_runs=1
       
       !get and save psi value for this flux surface 
       !(for testing)
       call get_psi_from_efit(xval_a(i),psif,relpsif)

       !perform iterations with increasing tracing distance 
       !and resolution to find correct safety factor
       do while (n_runs.le.n_iterations+1)
          if (n_runs.le.n_iterations) then
             distance=20.*pi*10**(0.4*n_runs)
             resolution=4.*distance*10**(0.2*n_runs)
          else
             distance=2.*q*pi
             resolution=tres
          endif
          call set_resolution(resolution,resol_b,resol_f)

          call initialize_tracer_aux(n_runs,resolution,tres,i)
          !tracing for every flux surface in the radial grid
          !starting position
          r=rzstart(i,1)
          p=0.
          z=rzstart(i,2)

          !find local symmetry plane, which has Br=0
          call symm_plane(r,zmin,zmax,zsym,ierr)
          if (ierr == 0)  then
             z=zsym
             call get_from_efit(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,&
                  &dBpdr,dBpdp,dBpdz,dBzdr,dBzdp,dBzdz)   
          else
             stop 'Error getting symmetry plane'
          endif

          call get_psirz_from_efit(r,z,psirz,relpsirz)
          if ((abs(relpsif-relpsirz)/relpsirz).gt.1E-2) then
             call get_rhotor_from_efit(r,z,rhorz)
             if (mype.eq.0) then
                write(*,'(A,2(F12.6,1X))') 'target: rhotor, rel. psi = ', &
                     &xval_a(i), relpsif
                write(*,'(A,2(F12.6,1X))') 'found:  rhotor, rel. psi = ', &
                     &rhorz, relpsirz
             endif
             stop 'missed flux surface in tracer'
          endif

          !initial condition for dx/dR at the starting position
          call get_c11(r,z,c11_int,Lref,Bref)

          dphi=distance/real(resolution-mod(resolution,2))
          hn1=0.1*dphi !guessed step for integrator

          !we start in the middle of our parallel grid and trace backwards
          k=resol_b
          
          !calculate quantities at this point
          call set_matrix_c_initial(r,z,Bref,c11_int)
          call calc_matrix_d
          call calc_metric_coefficients(r)
          call get_from_efit(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
               &dBpdz,dBzdr,dBzdp,dBzdz)   
          call calc_quantities(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,&
               &dBpdp,dBpdz,dBzdr,dBzdp,dBzdz,rmag,zmag)
          call prep_output(i,k,r,p,z,n_runs==n_iterations+1)
          !'err' determines the accuracy of the integrator and will 
          !have a strong influence on tracing speed
          err=epsilon(1.0)**0.6

          call set_rk_init_cond(r,z,Bp,Bz,Bref,c11_int,y)

          !backward tracing loop
          do kk=1,resol_b
             k=resol_b-kk
             p=p-dphi

             call odeint(y,p+dphi,p,err,hn1,0.) 
             
             call update_c_matrix(y,r,z)
             call calc_matrix_d
             call calc_metric_coefficients(r)
             call get_from_efit(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
                  &dBpdz,dBzdr,dBzdp,dBzdz)   
             call calc_quantities(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,&
                  &dBpdp,dBpdz,dBzdr,dBzdp,dBzdz,rmag,zmag)
             call prep_output(i,k,r,p,z,n_runs==n_iterations+1)
          enddo

          !return to starting position
          r=rzstart(i,1)
          p=0.
          z=zsym

          call get_from_efit(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
               &dBpdz,dBzdr,dBzdp,dBzdz) 
          call set_rk_init_cond(r,z,Bp,Bz,Bref,c11_int,y)

          !forward tracing loop
          do kk=1,resol_f
             k=resol_b+kk
             p=p+dphi

             call odeint(y,p-dphi,p,err,hn1,0.) 

             call update_c_matrix(y,r,z)
             call calc_matrix_d
             call calc_metric_coefficients(r)

             call get_from_efit(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
                  &dBpdz,dBzdr,dBzdp,dBzdz)   
             call calc_quantities(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,&
                  &dBpdp,dBpdz,dBzdr,dBzdp,dBzdz,rmag,zmag)
             call prep_output(i,k,r,p,z,n_runs==n_iterations+1)
          enddo

          !in the first iteration, we compute the safety factor and 
          !its radial derivative and compare with the data base
          if (n_runs.lt.n_iterations+1) then
             call calc_q(xval_a(i), q, dqdx, resolution,i)
             if (dqdx==-100.) call calc_dqdx(xval_a(i),resolution,dqdx,i)
             q_prof(i)=q
             dqdx_prof(i)=dqdx
          endif
          if (n_runs.gt.1.and.n_runs.ne.n_iterations+1) then
             acc(i)=abs(q-last_q)/q
             if (lilo) acc(pi1gl:pi2gl)=acc(pi1gl)
             !we require an accuracy of 0.01%
             if (acc(i).gt.1e-4) then
                n_runs=n_runs+1
                last_q=q
             else
                n_runs=n_iterations+1
             endif
          else
             n_runs=n_runs+1
             last_q=q
          endif
       enddo
       call get_q_from_efit(xval_a(i),qefit,shat_efit)
       if (abs((q-qefit)/qefit).gt.0.005) then
          !if (mype.eq.0) then
          if (write_pe) then
             write(*,'(A,F6.3)') 'WARNING: rel. q difference larger than 0.5% at x=',&
                  &xval_a(i)
             write(*,'(2(A,F7.4,1X),A,F4.1,A)') '         q_tracer = ',&
                  &q,'q_efit = ',qefit,'rel. diff: ',&
                  &abs((q-qefit)/qefit)*100.,'%'
          endif
       endif
       if (abs((dqdx*xval_a(i)/q-shat_efit)/shat_efit).gt.0.03.and.shat_efit.gt.0.1) then
          !if (mype.eq.0) then
          if (write_pe) then
             write(*,'(A,F6.3)') 'WARNING: rel. shat difference larger than 3.0% at x=',&
                  &xval_a(i)
             write(*,'(2(A,F7.4,1X),A,F4.1,A)') '         shat_tracer = ',&
                  &dqdx*xval_a(i)/q,'shat_efit = ',shat_efit,'rel. diff: ',&
                  &abs((dqdx*xval_a(i)/q-shat_efit)/shat_efit)*100.,'%'
          endif
       endif
       if (write_pe.and.(pi0.eq.1)) then
          write(*,"(2(A,F7.4,1X))") 'Tracing finished at rho_tor = ',xval_a(i)
          write(*,"(3(A,F7.4),A)") '(psi = ',psif,', rel. psi = ', &
               &relpsif,', rho_pol = ', sqrt(relpsif),')'
       !printing relative psi for comparison with other efit plotters
       endif
    enddo
    if (any(acc.gt.1e-4)) then
       write(*,"(A,G11.4,A)") 'Warning: The accuracy of the q-profile determination is limited to <',&
            maxval(acc)*100.,'%.'
    endif
    if (write_pe.and.pi0.gt.1) write(*,"(A)") 'Tracing finished.'
    
    major_R=rmag/Lref
    minor_r=1.0

    !for (approximate) computation of nustar
    trpeps = x0*minor_r/major_R

    tracer_initialized=.true.

  end  subroutine init_tracer

  !> Extracts geometry data from the tracer module.
  subroutine get_tracer(geom,&
       &C_y_out,C_xy_out,q_prof_out,dqdx_prof_out,dpdx_prof_out,&
       &q0_out,shat_out, Lref_out, Bref_out,edge_opt)  
    type(geomtype),intent(inout):: geom
    real, dimension(pi1gl:pi2gl),intent(out):: C_y_out,C_xy_out,q_prof_out,dqdx_prof_out,dpdx_prof_out
    real, dimension(pi1gl:pi2gl,lk1:lk2):: gxx,gxy,gxz,gyy,gyz,gzz,&
         Bfield,dBdx,dBdy,dBdz,jacobian,geo_R,geo_p,geo_Z,geo_c1,geo_c2
    real, dimension(pi1gl:pi2gl):: C_y,C_xy,q_prof,dqdx_prof, pres_prof, dpdx_prof
    real, intent(out):: q0_out,shat_out
    real:: q0, shat, Bref_out, Lref_out, edge_opt !, pref

    if (edge_opt.ne.0.) then
       call get_from_tracer_aux(gxx,gxy,gxz,gyy,gyz,gzz,Bfield,dBdx,dBdy,dBdz,jacobian,geo_R,geo_p,geo_Z,&
            &geo_c1,geo_c2,C_xy,C_y,q_prof,dqdx_prof,q0,shat,Bref,Lref,xval_a,x0,edge_opt,edge_opt_res,rmag,zmag)
    else
       call get_from_tracer_aux(gxx,gxy,gxz,gyy,gyz,gzz,Bfield,dBdx,dBdy,dBdz,jacobian,geo_R,geo_p,geo_Z,&
            &geo_c1,geo_c2,C_xy,C_y,q_prof,dqdx_prof,q0,shat,Bref,Lref,xval_a,x0,edge_opt,nz0,rmag,zmag)
    endif

    call get_pres_dpdx_from_efit(xval_a,pres_prof,dpdx_prof)

    geom%gij(pi1gl:pi2gl,pj1,lk1:lk2)=gxy(pi1gl:pi2gl,lk1:lk2)
    geom%gzz(pi1gl:pi2gl,pj1,lk1:lk2)=gzz(pi1gl:pi2gl,lk1:lk2)
    if (yx_order) then
       geom%gii(pi1gl:pi2gl,pj1,lk1:lk2)=gyy(pi1gl:pi2gl,lk1:lk2)
       geom%giz(pi1gl:pi2gl,pj1,lk1:lk2)=gyz(pi1gl:pi2gl,lk1:lk2)
       geom%gjj(pi1gl:pi2gl,pj1,lk1:lk2)=gxx(pi1gl:pi2gl,lk1:lk2)
       geom%gjz(pi1gl:pi2gl,pj1,lk1:lk2)=gxz(pi1gl:pi2gl,lk1:lk2)
       geom%dBdj(pi1gl:pi2gl,pj1,lk1:lk2)=dBdx(pi1gl:pi2gl,lk1:lk2)
       geom%dBdi(pi1gl:pi2gl,pj1,lk1:lk2)=dBdy(pi1gl:pi2gl,lk1:lk2)
    else
       geom%gii(pi1gl:pi2gl,pj1,lk1:lk2)=gxx(pi1gl:pi2gl,lk1:lk2)
       geom%giz(pi1gl:pi2gl,pj1,lk1:lk2)=gxz(pi1gl:pi2gl,lk1:lk2)
       geom%gjj(pi1gl:pi2gl,pj1,lk1:lk2)=gyy(pi1gl:pi2gl,lk1:lk2)
       geom%gjz(pi1gl:pi2gl,pj1,lk1:lk2)=gyz(pi1gl:pi2gl,lk1:lk2)
       geom%dBdi(pi1gl:pi2gl,pj1,lk1:lk2)=dBdx(pi1gl:pi2gl,lk1:lk2)
       geom%dBdj(pi1gl:pi2gl,pj1,lk1:lk2)=dBdy(pi1gl:pi2gl,lk1:lk2)
    end if

    geom%Bfield(pi1gl:pi2gl,pj1,lk1:lk2)  =Bfield(pi1gl:pi2gl,lk1:lk2)
    geom%dBdz(pi1gl:pi2gl,pj1,lk1:lk2)    =dBdz(pi1gl:pi2gl,lk1:lk2)
    geom%jacobian(pi1gl:pi2gl,pj1,lk1:lk2)=jacobian(pi1gl:pi2gl,lk1:lk2)
    geom%R(pi1gl:pi2gl,lk1:lk2)           =geo_R(pi1gl:pi2gl,lk1:lk2)
    geom%PHI(pi1gl:pi2gl,lk1:lk2)         =geo_p(pi1gl:pi2gl,lk1:lk2)
    geom%Z(pi1gl:pi2gl,lk1:lk2)           =geo_Z(pi1gl:pi2gl,lk1:lk2)
    geom%R_hat(pi1gl:pi2gl,lk1:lk2)       =geo_R(pi1gl:pi2gl,lk1:lk2)/Lref
    geom%Z_hat(pi1gl:pi2gl,lk1:lk2)       =geo_Z(pi1gl:pi2gl,lk1:lk2)/Lref
    geom%dxdR(pi1gl:pi2gl,lk1:lk2)        =geo_c1(pi1gl:pi2gl,lk1:lk2)
    geom%dxdZ(pi1gl:pi2gl,lk1:lk2)        =geo_c2(pi1gl:pi2gl,lk1:lk2)
    C_y_out=C_y
    C_xy_out=C_xy
    q_prof_out=q_prof
    dqdx_prof_out=dqdx_prof
    q0_out=q0
    shat_out=shat
    Lref_out=Lref
    Bref_out=Bref
!    if (x_local.or.lilo) then
!       pref = pres_prof(pi1gl)
!    else
!       call lag3interp(pres_prof,xval_a,px0,pref,x0)
!    endif


    if (mype==0.and.print_ini_msg) then
       write(*,"(A,F6.3,A,F6.3,A)") 'Tracing results: q0 =', q0,', shat = ', shat,'.'
       write(*,*)
!       write(*,'(A)') 'from efit:'
!       write(*,'(A)') ' x0           pres         beta_p       dpdx_p'
!       write(*,'(4(ES12.5,1X))') x0, pref, 8.0E-7*pi*pref/Bref**2, -dpdx_prof/pref
    endif

    dpdx_prof_out = -dpdx_prof*8.0E-7*pi/Bref**2.0

  end subroutine get_tracer

  subroutine finalize_tracer
    deallocate(xval_a)
    call finalize_efit
    call finalize_tracer_aux(.true.)
    tracer_initialized=.false.
  end subroutine finalize_tracer

  subroutine get_efit_refvals(Bref,Lref,major_R,minor_r,geomdir,geomfile)
    real, intent(out):: Bref, Lref, major_R, minor_r
    character(len=100), intent(in) :: geomdir
    character(len=80), intent(in)  :: geomfile
    
    call initialize_efit_small(geomdir,geomfile, Lref, Bref, rmag,tracer_initialized)
    major_R=rmag/Lref
    minor_r=1.0

  end subroutine get_efit_refvals
  
end module tracer_mod
