!>Auxiliaries for the tracer module.
module tracer_aux
  use discretization
  use tracer_IO
  use lagrange_interpolation
  use file_io, only: get_unit_nr
  use par_mod, only: lilo

  implicit none
  public:: symm_plane, set_resolution, set_matrix_c_initial, calc_quantities, calc_metric_coefficients,&
       calc_matrix_d, set_rk_init_cond, update_c_matrix, initialize_tracer_aux, finalize_tracer_aux,&
       prep_output, calc_q, calc_dqdx, get_from_tracer_aux, write_add_info
  private
  real, dimension(3,3) :: c,d
  real, dimension(3,3) :: gij,g_ij
  real:: Bfield, dBdv1, dBdv2, dBdtau, jac, Bphi, Btheta, q0, shat
  real, dimension(:,:),allocatable:: gxx, gxy, gxz, gyy, gyz, gzz, jac_array, Bfield_array, Bphi_arr, Btheta_arr
  real, dimension(:,:),allocatable:: geo_R, geo_p, geo_Z, geo_c1, geo_c2, dBdx_array, dBdy_array, dBdz_array
  real, dimension(:),allocatable:: gxx_n, gxy_n, gxz_n, gyy_n, gyz_n, gzz_n, jac_array_n
  real, dimension(:),allocatable:: Bfield_array_n, Bphi_arr_n, Btheta_arr_n, geo_R_n, geo_p_n, geo_Z_n
  real, dimension(:),allocatable:: geo_c1_n, geo_c2_n, dBdx_array_n, dBdy_array_n, dBdz_array_n
  real, dimension(:),allocatable:: q_prof_mod, dqdx_prof_mod, C_xy, C_y
  
  !additional quantities  
  real:: area, fsa_sgxx, avgshat, avg_sgxx
  logical:: initialized=.false.
contains

  subroutine write_add_info(rmag,zmag,Lref,x0,res)
    real, intent(in):: rmag, zmag,Lref,x0
    real:: pi,area, fsa_sgxx
    integer:: kmin=-1, kmax=-2,k,fileid,res

    if (mype==0) then
       pi=4.*atan(1.)
       call get_unit_nr(fileid)
       open(fileid,file='geom_info',status='replace')
       area=0
       fsa_sgxx=0
       do k=0,res-1
          area=area+(2*pi)**2*sqrt(gxx(pi1gl,k))*jac_array(pi1gl,k)*x0/q0
          fsa_sgxx=fsa_sgxx+sqrt(gxx(pi1gl,k))*jac_array(pi1gl,k)
       enddo
       write(fileid,*) 'Flux surface area:', area*Lref**2/res
       k=0
       do while (atan2((geo_Z(pi1gl,k)-zmag),(geo_R(pi1gl,k)-rmag)).lt.-pi*7/18) 
          k=k+1
          kmin=k
       enddo
       k=0
       do while (atan2((geo_Z(pi1gl,k)-zmag),(geo_R(pi1gl,k)-rmag)).lt.pi*7/18) 
          kmax=k
          k=k+1
       enddo
       avg_sgxx=0
       do k=kmin,kmax
          avg_sgxx=avg_sgxx+sqrt(gxx(pi1gl,k))*jac_array(pi1gl,k)/sum(jac_array(pi1gl,kmin:kmax))
       enddo
       write(fileid,*) 'R=1/2*(Rout+Rin)=', 0.5*(maxval(geo_R)+minval(geo_R))
       write(fileid,*) 'r=1/2*(Rout-Rin)',0.5*(maxval(geo_R)-minval(geo_R))
       write(fileid,*) 'dr/drho',0.5*(1./sqrt(gxx(pi1gl,maxloc(geo_R(pi1gl,:))-1))+&
            1./sqrt(gxx(pi1gl,minloc(geo_R(pi1gl,:))-1))) 
        write (fileid,*) 'r_outb=',geo_R(pi1gl,nz0/2)-rmag
       write (fileid,*) 'r_avg=',sum(sqrt((geo_R(pi1gl,:)-rmag)**2+(geo_Z(pi1gl,:)-zmag)**2))/nz0
       if (kmax.ge.kmin) &
            write(fileid,*) 'Avg. outboard shear:', &
            (gxy(pi1gl,kmax)/gxx(pi1gl,kmax)-gxy(pi1gl,kmin)/gxx(pi1gl,kmin))/(2.*pi)*res/(kmax-kmin+1)
       write(fileid,*) 'FSA sqrt(gxx):', fsa_sgxx/sum(jac_array(pi1gl,:))
       write(fileid,*) 'Avg. outboard sqrt(gxx):', avg_sgxx
       
       close(fileid)
    endif

  end subroutine write_add_info

  subroutine calc_q(xval_a, q, dqdx, res,i)
    real,intent(in):: xval_a
    real,intent(out):: q, dqdx
    real:: q_comp, C_y_comp
    integer, intent(in):: res, i
    real,dimension(0:res):: deviation
    real :: ftor_pr, fpol_pr, deriv, pi !, fpol_efit, ffprime_efit, dfpol_efit, test
    logical:: force_periodicity=.true.
    integer:: k, count, count2, jump, position, negderiv, posderiv, i_ind
!    real :: rcentr, zcentr, ang, ang0, ang_last

    pi=4.*atan(1.)

    !q0=dF_tor/dF_pol
    ftor_pr=sum(jac_array_n*Bphi_arr_n/geo_R_n)
    fpol_pr=sum(jac_array_n*Btheta_arr_n)
    q_comp=abs(ftor_pr/fpol_pr)
    C_y_comp=xval_a/q_comp

!    rcentr = SUM(geo_R_n*jac_array_n)/SUM(jac_array_n)
!    zcentr = SUM(geo_Z_n*jac_array_n)/SUM(jac_array_n)

    if (force_periodicity) then
       deviation=0.
       negderiv=0
       posderiv=0
       count=0
       count2 = 0
       jump = 0
       position=0

!       ang0 = atan((geo_Z_n(0)-zcentr)/(geo_R_n(0)-rcentr))
!       ang_last = 0.0
       do k=1,res-1
!          ang = atan((geo_Z_n(k)-zcentr)/(geo_R_n(k)-rcentr))-ang0+pi*jump
!          if (ang.lt.ang_last) then
!             ang = ang+pi
!             jump = jump+1
!          endif
!          ang_last = ang

          deviation(k)= &
               !(Bfield_array(i,k)-Bfield_array(i,0))**2+(jac_array(i,k)-jac_array(i,0))**2+&
               !(gxx(i,k)-gxx(i,0))**2+(dBdz_array(i,k)-dBdz_array(i,0))**2+&
                  !(gxz(i,k)-gxz(i,0))**2+(gzz(i,k)-gzz(i,0))**2+
               (geo_R_n(k)-geo_R_n(0))**2+&
               (geo_Z_n(k)-geo_Z_n(0))**2
          deriv=(deviation(k)-deviation(k-1))/(geo_p_n(k)-geo_p_n(k-1))
          !keep track of long time behavior of derivative
          if (deriv.lt.0) then
             negderiv=negderiv+1
             if (negderiv.gt.5) posderiv=0
          endif
          !we count minima in the deviation only if deviation is small and had been decreasing for some time
          if (deriv.gt.0.and.negderiv.gt.10.and.deviation(k-1).lt.1e-2*xval_a**2) then
             count=count+1
             position=k-1
             negderiv=0
          endif
          !keep track of long time behavior of derivative
          if (deriv.gt.0) then
             posderiv=posderiv+1
             if (posderiv.gt.5) negderiv=0
          endif
       enddo
!       count2 = INT(ang_last/(2.*pi))

       if (count.gt.0) then
          q=(geo_p_n(position)-geo_p_n(0))/(2.*pi)/count
          C_y(i)=xval_a/q
          dqdx=q/xval_a/(2.*count*pi)*(gxy_n(position)/gxx_n(position)-gxy_n(0)/gxx_n(0))
          if (abs((q_comp-q)/q_comp).gt.0.05) then
             !if (mype==0) write (*,'(A)') 'Using standard q profile.'
             q=q_comp
             C_y(i)=C_y_comp
             dqdx=-100.
          else
             
             !if (mype==0) write (*,'(A)') 'Enforcing periodicity on geometry coefficients.'
          endif
       else
          print*, 'Warning: count=0 in force_periodicity'
          q=q_comp
          C_y(i)=C_y_comp
          dqdx=-100.          
       endif
    else
       q=q_comp
       C_y(i)=C_y_comp
       dqdx=-100.        
    endif
    q_prof_mod(i)=q
    dqdx_prof_mod(i)=dqdx

    !for lilo runs, we copy q etc. to all radial indices
    if (lilo) then
       do i_ind=pi1gl,pi2gl
          q_prof_mod(i_ind)=q
          dqdx_prof_mod(i_ind)=dqdx
          C_y(i_ind)=C_y(i)
       enddo
    endif

  end subroutine calc_q

  subroutine calc_dqdx(xval_a,resolution, dqdx,i)
    real, intent(in):: xval_a
    real,intent(out):: dqdx    
    integer, intent(in):: resolution, i
    real, dimension(resolution) :: y
    real :: intercept, slope
    
    y(:)=gxy_n(:)/gxx_n(:)
    call least_sq(geo_p_n(:),y,slope,intercept)
    dqdx=slope*q_prof_mod(i)**2/xval_a
    dqdx_prof_mod(i)=dqdx

  contains
    
    subroutine least_sq(x,y,m,c)
      real, dimension(:), intent(IN) :: x, y
      real, intent(OUT) :: m, c
      
      integer :: n
      real :: sum_x, sum_xsq, sum_xy, sum_y
      
      n=size(x)
      
      !Evaluate sums
      
      sum_x=sum(x)
      sum_y=sum(y)
      sum_xy=dot_product(x,y)
      sum_xsq=dot_product(x,x)
      
      !Evaluate coefficients m & c
      
      m=(n*sum_xy-sum_x*sum_y)/  &
           (n*sum_xsq-sum_x*sum_x)
      c=(sum_y-m*sum_x)/n
    end subroutine least_sq
    
  end subroutine calc_dqdx

  subroutine coordinate_transformations(res,Lref,Bref,x0,edge_opt)
    integer:: res,k
    real:: Lref, Bref, x0, edge_opt
    real, dimension(pi1gl:pi2gl):: gthetaphi

    !first we transform to a radially independent y coordinate
    C_xy=C_y*q0/x0
    do k=0,res-1
       gxy(:,k)=gxy(:,k)/C_xy
       gyz(:,k)=gyz(:,k)/C_xy
       gyy(:,k)=gyy(:,k)/C_xy**2
       dBdy_array(:,k)=dBdy_array(:,k)*C_xy
       jac_array(:,k)=1./sqrt(gxx(:,k)*(gyy(:,k)*gzz(:,k)-gyz(:,k)**2)-gxy(:,k)*(gxy(:,k)*gzz(:,k)-&
       gyz(:,k)*gxz(:,k))+gxz(:,k)*(gxy(:,k)*gyz(:,k)-gxz(:,k)*gyy(:,k)))
    enddo
    C_y=x0/q0

    !we normalize the metrics
    do k=0,res-1
       gxx(:,k)=gxx(:,k)
       gxy(:,k)=gxy(:,k)
       gxz(:,k)=gxz(:,k)*Lref/q_prof_mod
       gyy(:,k)=gyy(:,k)
       gyz(:,k)=gyz(:,k)*Lref/q_prof_mod
       gzz(:,k)=gzz(:,k)*(Lref/q_prof_mod)**2
       jac_array(:,k)=jac_array(:,k)/Lref*q_prof_mod
       Bfield_array(:,k)=Bfield_array(:,k)/Bref
       dBdx_array(:,k)=dBdx_array(:,k)/Bref*Lref
       dBdy_array(:,k)=dBdy_array(:,k)/Bref*Lref
       dBdz_array(:,k)=dBdz_array(:,k)/Bref*q_prof_mod
    enddo

    !Now we convert the normalized tracer output in order to use the straight field line angle as parallel coordinate,
    !We do this by writing x, y and z as $x=\rho$, $y=\frac{\rho_0}{q_0}(q\theta-\phi)$ and $z=\frac{\phi}{q}$,
    !We transform to x'=x, y'=y, $z'=\theta$ by writing the metrics in terms of $\rho$, $\theta$ and $\phi$,
    !dBdy will be zero after this transformation, since the e_y basis vector then points in the toroidal direction.
    !dBdx also has to be modified, since it formerly contained a term proportional to (the finite) dBdy.
    !The jacobian should actually be the same afterwards, but we recalculate it anyway here.
    do k=0,res-1
       gthetaphi=q0/x0*gyz(:,k)+q_prof_mod*gzz(:,k)
       gxz(:,k)=1/x0*q0/q_prof_mod*gxy(:,k)-1./q_prof_mod**2*geo_p(:,k)*dqdx_prof_mod*gxx(:,k)
       gzz(:,k)=1./x0**2*q0**2/q_prof_mod**2*gyy(:,k)-1/q_prof_mod**4*dqdx_prof_mod**2*(geo_p(:,k))**2*&
         gxx(:,k)-gzz(:,k)-2./q_prof_mod**2*gxz(:,k)*geo_p(:,k)*dqdx_prof_mod+2./q_prof_mod*gthetaphi
       gyz(:,k)=x0*q_prof_mod/q0*gzz(:,k)+x0/q0/q_prof_mod*dqdx_prof_mod*geo_p(:,k)*gxz(:,k)-x0/q0*gthetaphi
       dBdx_array(:,k)=dBdx_array(:,k)+dBdy_array(:,k)*x0/q0/q_prof_mod*geo_p(:,k)*dqdx_prof_mod
       dBdy_array(:,k)=0.
       jac_array(:,k)=1./sqrt(gxx(:,k)*(gyy(:,k)*gzz(:,k)-gyz(:,k)**2)-gxy(:,k)*(gxy(:,k)*gzz(:,k)&
            -gyz(:,k)*gxz(:,k))+gxz(:,k)*(gxy(:,k)*gyz(:,k)-gxz(:,k)*gyy(:,k)))
    enddo
    if (edge_opt==0.) initialized=.true.


  end subroutine coordinate_transformations

  !>Normalize and transfer metrics to the tracer main module, where they are retrieved by the geometry module.
  !! In addition, we also do some conversions to make the y-coordinate radially independent, and to change the 
  !! definition of the parallel coordinate from z=phi/q to z=theta.
  subroutine get_from_tracer_aux(gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,Bfield_out,&
       dBdx_out,dBdy_out,dBdz_out,jacobian_out,geo_R_out,geo_p_out,geo_Z_out,geo_c1_out,geo_c2_out,&
       C_xy_out,C_y_out,q_prof_out,dqdx_prof_out,q0_out,shat_out,Bref,Lref,xval_a,x0,edge_opt,res,rmag, zmag)
    real, dimension(pi1gl:pi2gl,lk1:lk2),intent(out):: gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,&
         Bfield_out,dBdx_out,dBdy_out,dBdz_out,jacobian_out,geo_R_out,geo_p_out,geo_Z_out,geo_c1_out,geo_c2_out
    real, dimension(pi1gl:pi2gl), intent(out):: C_xy_out, C_y_out,q_prof_out,dqdx_prof_out
    real, dimension(pi1gl:pi2gl), intent(in):: xval_a
    real, intent(out):: q0_out,shat_out
    real, intent(in):: Bref, Lref,x0, rmag, zmag, edge_opt
    integer, intent(in):: res
    !the following variables are for the edge_opt case
    real, dimension(0:res-1):: zval, dzprimedz
    real, dimension(lk1:lk2):: zprime
    real, dimension(0:res-1):: tmp
    real:: pi
    integer:: i,k

    pi=4.*atan(1.)

    if (x_local) then
       q0=q_prof_mod(pi1gl)
       shat=dqdx_prof_mod(pi1gl)*xval_a(pi1gl)/q_prof_mod(pi1gl)
    else
       if (mod((pi1gl+pi2gl),2).eq.0) then !odd number of radial points
          q0=q_prof_mod((pi2gl+pi1gl)/2)
          shat=dqdx_prof_mod((pi2gl+pi1gl)/2)*xval_a((pi2gl+pi1gl)/2)/q_prof_mod((pi2gl+pi1gl)/2)
       else !even number of radial points
          !estimate center by linear interpolation of neighboring grid points
          q0=0.5*(q_prof_mod((pi1gl+pi2gl)/2)+q_prof_mod((pi1gl+pi2gl+1)/2))
          shat=0.5*(dqdx_prof_mod((pi1gl+pi2gl)/2)*xval_a((pi2gl+pi1gl)/2)/q_prof_mod((pi2gl+pi1gl)/2)+&
               dqdx_prof_mod((pi1gl+pi2gl+1)/2)*xval_a((pi2gl+pi1gl+1)/2)/q_prof_mod((pi2gl+pi1gl+1)/2))
       end if
    endif

    if (.not.initialized) call coordinate_transformations(res,Lref,Bref,x0,edge_opt)
    if (x_local.and.edge_opt==0) call write_add_info(rmag,zmag,Lref,x0,res)
    if (edge_opt.ne.0.) then
       if (.not.initialized) then
          if (x_local) call write_add_info(rmag,zmag,Lref,x0,res)
          do k=0,res-1
             zval(k)=-pi+k*2.*pi/res
             dzprimedz(k)=edge_opt*pi/log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/&
                  sqrt(edge_opt**2*zval(k)**2+1)
             gxz(:,k)=gxz(:,k)*dzprimedz(k)
             gyz(:,k)=gyz(:,k)*dzprimedz(k)
             gzz(:,k)=gzz(:,k)*dzprimedz(k)**2
             jac_array(:,k)=jac_array(:,k)/dzprimedz(k)
             dBdz_array(:,k)=dBdz_array(:,k)/dzprimedz(k)
          enddo
          initialized=.true.
       endif
       do k=0,res-1
          zval(k)=-pi+k*2.*pi/res
       enddo
       do k=lk1,lk2
          zprime(k)=sinh((-pi+k*2.*pi/nz0)*log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/pi)/edge_opt
       enddo
       do i=pi1gl,pi2gl
          tmp=gxx(i,0:res-1)
          call lag3interp(tmp,zval,res,gxx_out(i,lk1:lk2),zprime,lk0)
          tmp=gxy(i,0:res-1)
          call lag3interp(tmp,zval,res,gxy_out(i,lk1:lk2),zprime,lk0)
          tmp=gxz(i,0:res-1)
          call lag3interp(tmp,zval,res,gxz_out(i,lk1:lk2),zprime,lk0)
          tmp=gyy(i,0:res-1)
          call lag3interp(tmp,zval,res,gyy_out(i,lk1:lk2),zprime,lk0)
          tmp=gyz(i,0:res-1)
          call lag3interp(tmp,zval,res,gyz_out(i,lk1:lk2),zprime,lk0)
          tmp=gzz(i,0:res-1)
          call lag3interp(tmp,zval,res,gzz_out(i,lk1:lk2),zprime,lk0)
          tmp=dBdx_array(i,0:res-1)
          call lag3interp(tmp,zval,res,dBdx_out(i,lk1:lk2),zprime,lk0)
          tmp=dBdy_array(i,0:res-1)
          call lag3interp(tmp,zval,res,dBdy_out(i,lk1:lk2),zprime,lk0)
          tmp=dBdz_array(i,0:res-1)
          call lag3interp(tmp,zval,res,dBdz_out(i,lk1:lk2),zprime,lk0)
          tmp=jac_array(i,0:res-1)
          call lag3interp(tmp,zval,res,jacobian_out(i,lk1:lk2),zprime,lk0)
          tmp=Bfield_array(i,0:res-1)
          call lag3interp(tmp,zval,res,Bfield_out(i,lk1:lk2),zprime,lk0)
          tmp=geo_R(i,0:res-1)
          call lag3interp(tmp,zval,res,geo_R_out(i,lk1:lk2),zprime,lk0)
          tmp=geo_Z(i,0:res-1)
          call lag3interp(tmp,zval,res,geo_Z_out(i,lk1:lk2),zprime,lk0)
          tmp=geo_p(i,0:res-1)
          call lag3interp(tmp,zval,res,geo_p_out(i,lk1:lk2),zprime,lk0)
          tmp=geo_c1(i,0:res-1)
          call lag3interp(tmp,zval,res,geo_c1_out(i,lk1:lk2),zprime,lk0)
          tmp=geo_c2(i,0:res-1)
          call lag3interp(tmp,zval,res,geo_c2_out(i,lk1:lk2),zprime,lk0)
       enddo
    else
       do k=lk1,lk2
          gxx_out(:,k)=gxx(:,k)
          gxy_out(:,k)=gxy(:,k)
          gxz_out(:,k)=gxz(:,k)
          gyy_out(:,k)=gyy(:,k)
          gyz_out(:,k)=gyz(:,k)
          gzz_out(:,k)=gzz(:,k)
          jacobian_out(:,k)=jac_array(:,k)
          Bfield_out(:,k)=Bfield_array(:,k)
          dBdx_out(:,k)=dBdx_array(:,k)
          dBdy_out(:,k)=dBdy_array(:,k)
          dBdz_out(:,k)=dBdz_array(:,k)
          geo_R_out(:,k)=geo_R(:,k)
          geo_p_out(:,k)=geo_p(:,k)
          geo_Z_out(:,k)=geo_Z(:,k)
          geo_c1_out(:,k)=geo_c1(:,k)
          geo_c2_out(:,k)=geo_c2(:,k)
       enddo
    endif

    C_xy_out=C_xy
    C_y_out=C_y
    q_prof_out=q_prof_mod
    dqdx_prof_out=dqdx_prof_mod
    q0_out=q0
    shat_out=shat
  end subroutine get_from_tracer_aux

  !>Allocates any arrays required in the tracer_aux module. This routine is called twice, since we need to do
  !! two separate tracing runs with different resolutions.
  subroutine initialize_tracer_aux(n_runs,resolution,tres,i)
    integer, intent(in):: n_runs,resolution,tres,i
    if (i==pi1gl.and.n_runs==1) then
       allocate(q_prof_mod(pi1gl:pi2gl),dqdx_prof_mod(pi1gl:pi2gl),C_y(pi1gl:pi2gl),C_xy(pi1gl:pi2gl))
       allocate(gxx(pi1gl:pi2gl,0:tres-1),gxy(pi1gl:pi2gl,0:tres-1),gxz(pi1gl:pi2gl,0:tres-1))
       allocate(gyy(pi1gl:pi2gl,0:tres-1),gyz(pi1gl:pi2gl,0:tres-1),gzz(pi1gl:pi2gl,0:tres-1))
       allocate(jac_array(pi1gl:pi2gl,0:tres-1),Bfield_array(pi1gl:pi2gl,0:tres-1))
       allocate(dBdx_array(pi1gl:pi2gl,0:tres-1),dBdy_array(pi1gl:pi2gl,0:tres-1),dBdz_array(pi1gl:pi2gl,0:tres-1))
       allocate(geo_R(pi1gl:pi2gl,0:tres-1),geo_p(pi1gl:pi2gl,0:tres-1),geo_Z(pi1gl:pi2gl,0:tres-1))
       allocate(geo_c1(pi1gl:pi2gl,0:tres-1),geo_c2(pi1gl:pi2gl,0:tres-1))
       allocate(Bphi_arr(pi1gl:pi2gl,0:tres-1),Btheta_arr(pi1gl:pi2gl,0:tres-1))
    else
       call finalize_tracer_aux(.false.)
    endif
    allocate(gxx_n(0:resolution-1),gxy_n(0:resolution-1),gxz_n(0:resolution-1))
    allocate(gyy_n(0:resolution-1),gyz_n(0:resolution-1),gzz_n(0:resolution-1))
    allocate(jac_array_n(0:resolution-1),Bfield_array_n(0:resolution-1))
    allocate(dBdx_array_n(0:resolution-1),dBdy_array_n(0:resolution-1),dBdz_array_n(0:resolution-1))
    allocate(geo_R_n(0:resolution-1),geo_p_n(0:resolution-1),geo_Z_n(0:resolution-1))
    allocate(geo_c1_n(0:resolution-1),geo_c2_n(0:resolution-1))
    allocate(Bphi_arr_n(0:resolution-1),Btheta_arr_n(0:resolution-1))

  end subroutine initialize_tracer_aux

  subroutine finalize_tracer_aux(all)
    logical::  all
    if (all) then
       deallocate(gxx,gxy,gxz,gyy,gyz,gzz,jac_array,Bfield_array,dBdx_array,dBdy_array,dBdz_array)
       deallocate(geo_R,geo_p,geo_Z,geo_c1,geo_c2,Bphi_arr,Btheta_arr)
       deallocate(gxx_n,gxy_n,gxz_n,gyy_n,gyz_n,gzz_n,jac_array_n,Bfield_array_n,dBdx_array_n,dBdy_array_n,dBdz_array_n)
       deallocate(geo_R_n,geo_p_n,geo_Z_n,geo_c1_n,geo_c2_n,Bphi_arr_n,Btheta_arr_n)
       deallocate(q_prof_mod,dqdx_prof_mod,C_xy,C_y)
       initialized=.false.
    else 
       deallocate(gxx_n,gxy_n,gxz_n,gyy_n,gyz_n,gzz_n,jac_array_n,Bfield_array_n,dBdx_array_n,dBdy_array_n,dBdz_array_n)
       deallocate(geo_R_n,geo_p_n,geo_Z_n,geo_c1_n,geo_c2_n,Bphi_arr_n,Btheta_arr_n)
    endif
  end subroutine finalize_tracer_aux

  !>Copies data to output arrays.
  subroutine prep_output(i,k,r,p,z,final)
    integer, intent(in):: i,k
    real, intent(in):: r,p,z
    logical, intent(in):: final
    integer:: i_ind
    
    gxx_n(k)=gij(1,1)
    gxy_n(k)=gij(1,2)
    gxz_n(k)=gij(1,3)
    gyy_n(k)=gij(2,2)
    gyz_n(k)=gij(2,3)
    gzz_n(k)=gij(3,3)
    jac_array_n(k)=jac
    Bfield_array_n(k)=Bfield
    dBdx_array_n(k)=dBdv1
    dBdy_array_n(k)=dBdv2
    dBdz_array_n(k)=dBdtau
    Bphi_arr_n(k)=Bphi
    Btheta_arr_n(k)=Btheta
    geo_R_n(k)=r
    geo_p_n(k)=p
    geo_Z_n(k)=z
    geo_c1_n(k)=c(1,1)
    geo_c2_n(k)=c(1,2)

    if (final) then
       !for lilo runs, we copy the geometry to every radial index
       if (lilo) then
          do i_ind=pi1gl,pi2gl
             gxx(i_ind,k)=gij(1,1)
             gxy(i_ind,k)=gij(1,2)
             gxz(i_ind,k)=gij(1,3)
             gyy(i_ind,k)=gij(2,2)
             gyz(i_ind,k)=gij(2,3)
             gzz(i_ind,k)=gij(3,3)
             jac_array(i_ind,k)=jac
             Bfield_array(i_ind,k)=Bfield
             dBdx_array(i_ind,k)=dBdv1
             dBdy_array(i_ind,k)=dBdv2
             dBdz_array(i_ind,k)=dBdtau
             Bphi_arr(i_ind,k)=Bphi
             Btheta_arr(i_ind,k)=Btheta
             geo_R(i_ind,k)=r
             geo_p(i_ind,k)=p
             geo_Z(i_ind,k)=z
             geo_c1(i_ind,k)=c(1,1)
             geo_c2(i_ind,k)=c(1,2)
          enddo
       else
          gxx(i,k)=gij(1,1)
          gxy(i,k)=gij(1,2)
          gxz(i,k)=gij(1,3)
          gyy(i,k)=gij(2,2)
          gyz(i,k)=gij(2,3)
          gzz(i,k)=gij(3,3)
          jac_array(i,k)=jac
          Bfield_array(i,k)=Bfield
          dBdx_array(i,k)=dBdv1
          dBdy_array(i,k)=dBdv2
          dBdz_array(i,k)=dBdtau
          Bphi_arr(i,k)=Bphi
          Btheta_arr(i,k)=Btheta
          geo_R(i,k)=r
          geo_p(i,k)=p
          geo_Z(i,k)=z
          geo_c1(i,k)=c(1,1)
          geo_c2(i,k)=c(1,2)
       endif
    endif
  end subroutine prep_output

  !> Updates the matrix C^i_j=du^i/dx^j after a Runge-Kutta step.
  subroutine update_c_matrix(y,r,z)
    real, dimension(8), intent(in):: y
    real, intent(out):: r,z

    r = y(1)
    z = y(2)
    c(1,1) = y(3)
    c(1,2) = y(4)
    c(1,3) = y(5)
    c(2,1) = y(6)
    c(2,2) = y(7)
    c(2,3) = y(8)

  end subroutine update_c_matrix

  !> Initial condition for Runge-Kutta scheme.
  subroutine set_rk_init_cond(r,z,Bp,Bz,Bref,c11_0,y)
    real, dimension(8), intent(out):: y
    real, intent(in):: r,z,Bref,c11_0,Bp,Bz
    y(1)=r
    y(2)=z
    y(3)=c11_0
    y(4)=0.
    y(5)=0.
    y(6)=0.
    y(7)=+Bp/c11_0/Bref
    y(8)=-r*Bz/c11_0/Bref

  end subroutine set_rk_init_cond

  !> Compute derivatives and the jacobian. Some computations done in the original tracer code have been
  !! commented out, since they are not required in GENE, but might be of interest for special applications.
  subroutine calc_quantities(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,dBpdz,dBzdr,dBzdp,dBzdz,rmag,z0)
    real, intent(in):: r, z, Br, Bp, Bz, dBrdr, dBrdp, dBrdz, dBpdr, dBpdp, dBpdz, dBzdr, dBzdp, dBzdz, rmag, z0
    real:: divB, b3, dBdr, dBdp, dBdz, inv
    !Field and derivatives 

    Bfield=sqrt(Br**2+Bz**2+Bp**2)
    divB=Br/r+dBrdr+dBpdp/r+dBzdz

    if (abs(divB).gt.4.0*epsilon(divB)) then
       !for unknown reason, divB is slightly less accurate for
       !TRANSP equilibria and therefore 4*eps instead of 
       !1*eps is employed as upper boundary here (should be fine)
!       if (mype==0) &
!          & 
       print*, 'div B/epsilon = ', divB/epsilon(divB)
       stop 'magnetic field not divergence free in tracer'
    endif

    b3=Bp/(Bfield*r)  !b^/tau
    Btheta=Bz*(r-rmag)/((r-rmag)**2+(z-z0)**2)-Br*(z-z0)/((r-rmag)**2+(z-z0)**2) !B^theta
    Bphi=Bp

    !Derivs wrt cylindrical system

    dBdr=Br*dBrdr/Bfield+Bz*dBzdr/Bfield+Bp*dBpdr/Bfield
    dBdp=Br*dBrdp/Bfield+Bz*dBzdp/Bfield+Bp*dBpdp/Bfield
    dBdz=Br*dBrdz/Bfield+Bz*dBzdz/Bfield+Bp*dBpdz/Bfield
    
    !Derivs wrt Clebsch system
    
    dBdv1=d(1,1)*dBdr+d(2,1)*dBdz
    dBdv2=d(1,2)*dBdr+d(2,2)*dBdz 
    dBdtau=d(1,3)*dBdr+d(2,3)*dBdz+dBdp 
    
    !Jacobian (using 3 methods)
    
    !    jac_chain=-r/detc !Chain rule *gives also sign*
    !    jac_contravar = SQRT(1./detgij(gij)) !Via contravariant elements
    !    jac_covar = SQRT(detgij(g_ij)) !Via covariant elements
    jac=sqrt(1./detgij(gij))
    
    !Stream function JB^phi
    
    inv=(jac*Bfield*b3)
    
    !Cartesian coordinates
    
    !    cart_x=r*COS(phi)
    !    cart_y=r*SIN(phi)
    
    !Poloidal angle
    
    !    IF (z > z0 .AND. r > rmag)   theta=ATAN((z-z0)/(r-rmag))
    !    IF (z > z0 .AND. r < rmag)   theta=Pi-ATAN((z-z0)/(rmag-r))
    !    IF (z < z0 .AND. r > rmag)   theta=ATAN((z-z0)/(r-rmag))
    !    IF (z <= z0 .AND. r < rmag)   theta=-Pi+ATAN((z-z0)/(r-rmag))
    
    !----------------------------- Curvature ----------------------------------!
    
    !Auxilliaries
    
    !    fac1=(gij(1,3)*gij(2,2)-gij(1,2)*gij(2,3))/&
    !         &(gij(1,1)*gij(2,2)-gij(1,2)*gij(1,2))
    !    fac2=(gij(1,1)*gij(2,3)-gij(1,2)*gij(1,3))/&
    !         &(gij(1,1)*gij(2,2)-gij(1,2)*gij(1,2))
    !    fac3=SQRT(gij(1,1))*Bfield*jac*b3
    
    !Covariant components
    
    !    k1=dBdv1/Bfield+fac1*dBdtau/Bfield
    !    k2=dBdv2/Bfield+fac2*dBdtau/Bfield
    
    !Normal curvature
    
    !    k_norm=(gij(1,1)*dBdv1+gij(1,2)*dBdv2&
    !         &+gij(1,3)*dBdtau)/SQRT(gij(1,1))/Bfield
    
    !Geodesic curvature
    
    !    k_geo=-(dBdv2+fac2*dBdtau)/fac3              
    
    !square of Total curvature (explicit expression)
    
    !    kappa_sq=(dBdv1**2*gij(1,1)+dBdv2**2*gij(2,2)&
    !         &+dBdtau**2*(gij(3,3)-b3**2)&
    !         &+2*dBdv1*dBdv2*gij(1,2)+2*dBdv1*dBdtau*gij(1,3)&
    !         &+2*dBdv2*dBdtau*gij(2,3))/Bfield**2  
    
    !square of Total curvature (implicit expression) 
    
    !    kappa_sq_imp=k_norm**2+k_geo**2          
    
    !-------------------------- Consistency checks ----------------------------!
    
    !Curvature
    !    zero_kappa=ABS(kappa_sq-kappa_sq_imp)
    
    !Clebsch system
    
    !    zero_cl1=Br*c(1,1)+Bz*c(1,2)+Bp/r*c(1,3)
    !    zero_cl2=Br*c(2,1)+Bz*c(2,2)+Bp/r*c(2,3)
    !    zero_cl3=ABS(Bfield**2-inv**2*(gij(1,1)*gij(2,2)-gij(1,2)**2))
    
    !---------------------- Determinant of metrics -----------------------------!
  contains
    !---------------------------------------------------------------------------!
    
    function detgij(x)
      
      real :: detgij
      real, dimension(3,3) :: x
      
      detgij=x(1,1)*x(2,2)*x(3,3)+2*x(1,2)*x(2,3)*x(1,3)&
           &-x(2,2)*x(1,3)**2-x(1,1)*x(2,3)**2-x(3,3)*x(1,2)**2 
      
    end function detgij
    
    !----------------------------------------------------------------------------!
    
    
  end subroutine calc_quantities
  
  !> Compute co- and contravariant metrics from the C^i_j matrix.
  subroutine calc_metric_coefficients(r)
    real, intent(in):: r

    gij(1,1) = c(1,1)*c(1,1) + c(1,2)*c(1,2) + c(1,3)*c(1,3)/r**2
    gij(1,2) = c(1,1)*c(2,1) + c(1,2)*c(2,2) + c(1,3)*c(2,3)/r**2
    gij(1,3) = c(1,1)*c(3,1) + c(1,2)*c(3,2) + c(1,3)*c(3,3)/r**2
    gij(2,1) = gij(1,2) 
    gij(2,2) = c(2,1)*c(2,1) + c(2,2)*c(2,2) + c(2,3)*c(2,3)/r**2
    gij(2,3) = c(2,1)*c(3,1) + c(2,2)*c(3,2) + c(2,3)*c(3,3)/r**2
    gij(3,1) = gij(1,3)
    gij(3,2) = gij(2,3)
    gij(3,3) = c(3,1)*c(3,1) + c(3,2)*c(3,2) + c(3,3)*c(3,3)/r**2


    g_ij(1,1) = d(1,1)*d(1,1) + d(2,1)*d(2,1) + d(3,1)*d(3,1)*r**2
    g_ij(1,2) = d(1,1)*d(1,2) + d(2,1)*d(2,2) + d(3,1)*d(3,2)*r**2
    g_ij(1,3) = d(1,1)*d(1,3) + d(2,1)*d(2,3) + d(3,1)*d(3,3)*r**2
    g_ij(2,1) = g_ij(1,2) 
    g_ij(2,2) = d(1,2)*d(1,2) + d(2,2)*d(2,2) + d(3,2)*d(3,2)*r**2
    g_ij(2,3) = d(1,2)*d(1,3) + d(2,2)*d(2,3) + d(3,2)*d(3,3)*r**2
    g_ij(3,1) = g_ij(1,3)
    g_ij(3,2) = g_ij(2,3)
    g_ij(3,3) = d(1,3)*d(1,3) + d(2,3)*d(2,3) + d(3,3)*d(3,3)*r**2

  end subroutine calc_metric_coefficients


  subroutine set_matrix_c_initial(r,z,Bref,c11_0)
    real, intent(in):: r,z,Bref
    real, intent(in):: c11_0
    real::Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,dBpdz,dBzdr,dBzdp,dBzdz
    
    call get_from_efit(r,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
         &dBpdz,dBzdr,dBzdp,dBzdz)
    c=0.
    c(1,1)=c11_0
    c(2,2)=+Bp/c11_0/Bref
    c(2,3)=-r*Bz/c11_0/Bref
    c(3,3)=1.

  end subroutine set_matrix_c_initial

  !> Invert the C^i_j matrix to find the covariant conversion rules.
  subroutine calc_matrix_d
    real:: detc
    intrinsic TINY

    detc=c(1,1)*c(2,2)-c(1,2)*c(2,1)
    if (detc < tiny(1.0)) stop &
         '!Fatal error: (r0,z0) out of bounds or line hit boundary!'
    
    d(1,1)=c(2,2)/detc
    d(1,2)=-c(1,2)/detc
    d(1,3)=(c(1,2)*c(2,3)-c(1,3)*c(2,2))/detc
    d(2,1)=-c(2,1)/detc
    d(2,2)=c(1,1)/detc
    d(2,3)=(c(1,3)*c(2,1)-c(1,1)*c(2,3))/detc
    d(3,1)=0.
    d(3,2)=0.
    d(3,3)=1.
    
  end subroutine calc_matrix_d

  subroutine set_resolution(resol,resol_b,resol_f)
    integer, intent(in):: resol
    integer, intent(out):: resol_b, resol_f
    integer :: res_mod

    res_mod=mod(resol,2)
    resol_b=(resol-res_mod)/2
    resol_f=resol_b - 1 + res_mod
    
  end subroutine set_resolution


  !>Locates the symmetry plane for a given radial (R) position, where Br==0.
  subroutine symm_plane(r0,zmin,zmax,z0,ierr)

    real, intent(IN) :: r0,zmin,zmax 
    real, intent(OUT) :: z0 
    integer, intent(OUT) :: ierr
    real :: Br_min,Br_max 
    real :: left,right,Br,Bp,Bz,dBrdr,dBrdz,dBpdr,dBpdz,dBzdr,dBzdz
    real :: dBrdp,dBpdp,dBzdp,ztry_min,ztry_max,zdiff

   ztry_min=0.3*zmin ; ztry_max=0.3*zmax
   do
       call get_from_efit(r0,ztry_min,Br,Bp,Bz,dBrdr,dBrdp,&
            &dBrdz,dBpdr,dBpdp,dBpdz,dBzdr,dBzdp,dBzdz) 
       Br_min=Br

       call get_from_efit(r0,ztry_max,Br,Bp,Bz,dBrdr,dBrdp,&
            &dBrdz,dBpdr,dBpdp,dBpdz,dBzdr,dBzdp,dBzdz) 
       Br_max=Br
       if (Br_min < 0 .and. Br_max > 0) then
          left=ztry_min ; right=ztry_max
!          left = zmin; right = zmax
          exit
       elseif (Br_min > 0 .and. Br_max < 0) then
          left=ztry_max ; right=ztry_min
!          left = zmax; right=zmin
          exit
       else
          zdiff=abs(ztry_min-ztry_max)
          if (zdiff > 2.*spacing(zdiff)) then
             ztry_min=ztry_min/2. ; ztry_max=ztry_max/2.
          else
             print *, 'Symmetry plane could not be located!'
             ierr=1
             return
          endif
       endif
    enddo
    
    !bisect
    
    do
       z0=(left+right)/2. 

       call get_from_efit(r0,z0,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
            dBpdz,dBzdr,dBzdp,dBzdz) 

       if (Br < 0.) then
          left=z0
       elseif (Br > 0.) then
          right=z0
       else
!          print*, 'root exactly found'
          exit !root found
       endif
       
       if (abs(right-left) < 2.*spacing(z0)) then
!          print*, 'root found with best possible accuracy'
          exit
       endif
    enddo
    
    !Success
    
    ierr=0
    
  end subroutine symm_plane


end module tracer_aux
