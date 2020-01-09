#include "redef.h"
#include "intrinsic_sizes.h"
#undef GYRO_HERE
#undef WITH_NEO_F0
!> The diag_neoclass module is used for preparation and output
!! of neoclassical transport data
Module diagnostics_neoclass
  Use par_mod
  use aux_fields, only: f_, emfields
  Use communications
  Use file_io
  use geometry
  use spatial_averages
  Use vel_space !, only: mat_00, mat_10, fm

Implicit None
  PUBLIC :: initialize_all_diags_neoclass, initialize_all_diags_neoclass2D, exec_all_diags_neoclass, exec_all_diags_neoclass2D, & 
    & finalize_all_diags_neoclass,finalize_all_diags_neoclass2D, mem_est_diag_neoclass, istep_neoclass, istep_neoclass2D,& 
    & check_diag_neoclass,check_diag_neoclass2D, neoflux, set_mats_diag_neoclass, write_neoclassfile 


  PRIVATE 
  Integer:: istep_neoclass = -1, istep_neoclass2D = -1  
  CHARACTER(LEN=8) :: filestat='replace', filepos='rewind'
  Integer:: NEOCLASSFILE, NEOCLASSFILE2D 
  REAL, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE:: neoclass_mat
  Real, Dimension(:,:,:),ALLOCATABLE :: neoflux,neofluxf0
  Real, Dimension(:,:,:,:),ALLOCATABLE :: neoflux2D 
  Real, Dimension(:), ALLOCATABLE :: fnorm

  INTEGER :: n_neoclass_mats = 4

#ifdef G_DAGGER
  CHARACTER :: G_op = 'C'
#else
  CHARACTER :: G_op = 'N'
#endif

CONTAINS

!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Give an estimate of the memory requirements of diag_neoclass
  Real Function mem_est_diag_neoclass(mem_req_in)
    real, parameter:: mem_cn=SIZE_OF_COMPLEX/(1024.)**2,&
         & mem_rn=SIZE_OF_REAL/(1024.)**2
    real:: mem_req_in
    real:: mem_loc=0

    !neoclass_mats
    mem_loc = n_neoclass_mats*mem_cn*pi0*pj0*lklmn0
    
#ifdef WITH_NEO_F0
    !neofluxf0
    mem_loc = mem_loc + 2*px0*n_spec*mem_rn
#endif
    !neoflux
    mem_loc = mem_loc + n_neoclass_mats*nx0*n_spec*mem_rn

    !local variables
    !neo_prefs
    mem_loc = mem_loc + pi0*pj0*lk0*mem_rn
    !momc_neo
    mem_loc = mem_loc + lijk0*n_neoclass_mats*ln0*mem_cn
    !mom_neo
    mem_loc = mem_loc + lijk0*n_neoclass_mats*ln0*mem_rn
    !2D diagnostic
    mem_loc = mem_loc +  lijk0*n_neoclass_mats*ln0*mem_rn 
    mem_est_diag_neoclass = mem_loc+mem_req_in

  End Function mem_est_diag_neoclass

!!!******************************************************************!!!

  subroutine check_diag_neoclass
    logical :: write_pe

    write_pe = ((mype.le.0).AND.(print_ini_msg))
       
    if (xy_local.and.yx_order.and.(istep_neoclass.gt.0)) then
       if (write_pe) write(*,"(A)") "diag_neoclass not implemented for xy_local+yx_order --setting istep_neoclass=0"
       istep_neoclass=0
    endif
    
    if (istep_neoclass .gt. 0) then
       if (xy_local) then
          if (ky0_ind .eq. 0) then
             if (write_pe) write(*,"(A)") "with diag_neoclass: fields are zeroed out for ky=kx=0 in local code"
          else
             if (write_pe) write(*,"(A)") "diag_neoclass switched off because kx=ky=0 mode is not in the system"
             istep_neoclass = 0
          endif      
       end if
       
       if (.not.only_neo.and.include_f0_contr) then
          if (write_pe) write(*,"(A)") "WARNING: possible time scale separation between turbulence and neoclassics"
          if ((xy_local).and.(write_pe)) &
               &write(*,"(A)") "for neoclassical fluxes, computation with kx=ky=0 only is recommended"
       end if
    end if
#ifdef WITH_NEO_F0
    if (.not.y_local) stop "neoclassical f0 flux is only implemented for y_local"
#endif

  end subroutine check_diag_neoclass
!!!******************************************************************!!!
!!  - Review!!!
 subroutine check_diag_neoclass2D
    logical :: write_pe
    write_pe = ((mype.le.0).AND.(print_ini_msg))

   if (xy_local.and.yx_order.and.(istep_neoclass2D.gt.0)) then
       if (write_pe) write(*,"(A)") "diag_neoclass not implemented for xy_local+yx_order --setting istep_neoclass2D=0"
       istep_neoclass2D=0
    endif
    
    if (istep_neoclass2D .gt. 0) then
       if (xy_local) then
          if (ky0_ind .eq. 0) then
             if (write_pe) write(*,"(A)") "with diag_neoclass2D: fields are zeroed out for ky=kx=0 in local code"
          else
             if (write_pe) write(*,"(A)") "diag_neoclass2D switched off because kx=ky=0 mode is not in the system"
             istep_neoclass2D = 0
          endif      
       end if
       
       if (.not.only_neo.and.include_f0_contr) then
          if (write_pe) write(*,"(A)") "WARNING: possible time scale separation between turbulence and neoclassics"
          if ((xy_local).and.(write_pe)) &
               &write(*,"(A)") "for neoclassical fluxes, computation with kx=ky=0 only is recommended"
       end if
    end if
#ifdef WITH_NEO_F0
    if (.not.y_local) stop "neoclassical f0 flux is only implemented for y_local"
#endif

 end subroutine check_diag_neoclass2D


!!!******************************************************************!!!

  SUBROUTINE initialize_all_diags_neoclass

    call get_unit_nr(NEOCLASSFILE)
    if (mype==0) then
       OPEN(NEOCLASSFILE, file=trim(diagdir)//'/neoclass'&
            &//trim(file_extension), form='formatted',&
            status=filestat, position=filepos)
       if (x_local) then
          WRITE(NEOCLASSFILE,'(A)') '# <G_neo>     <Q_neo>       <Pi_neo>          <j_bootstrp>'
       else
          WRITE(NEOCLASSFILE,'(A)') '# <G_neo(x0)>  <Q_neo(x0)>  <Pi_neo(x0)>      <j_bootstrp(x0)>'
       endif
    endif

    ALLOCATE(neoclass_mat(pi1:pi2,pj1:pj2,lk1:lk2,1:n_neoclass_mats,ll1:ll2,lm1:lm2,ln1:ln2))
    if (.not.allocated(neoflux)) allocate(neoflux(0:nx0-1,0:n_spec-1,1:n_neoclass_mats))

    call set_mats_diag_neoclass(neoclass_mat,n_neoclass_mats)

#ifdef WITH_NEO_F0
    ! normalization (flux surface average)
    ALLOCATE(fnorm(pg1gl:pg2gl))
    if (xy_local) then
       fnorm = 1.0 / (Real(nz0)*geom%avg_jaco)
    elseif(.not.x_local) then
       fnorm = 1.0 / (REAL(nz0)*geom%avg_jaco_yz)
    else !.not.y_local
       fnorm = 1.0 / (REAL(nky0*nz0)*geom%avg_jaco_yz)
    endif
    
    IF (.NOT.ALLOCATED(neofluxf0)) ALLOCATE(neofluxf0(pg1gl:pg2gl,0:n_spec-1,0:1))
          
    call get_neofluxf0(neofluxf0)
    if (xy_local) then
       neoflux_f0_scalar=neofluxf0(pi1,:,:)
    elseif (.not.x_local)then
       neoflux_f0_scalar=(neofluxf0(nx0/2,:,:)+neofluxf0(nx0/2-1,:,:))/2.0
    else
       neoflux_f0_scalar=neofluxf0(pi1,:,:)
    endif
#endif


  END SUBROUTINE initialize_all_diags_neoclass
!!!******************************************************************!!!

SUBROUTINE initialize_all_diags_neoclass2D 

  call get_unit_nr(NEOCLASSFILE2D)
  
  if(.not.allocated(neoclass_mat)) ALLOCATE(neoclass_mat(pi1:pi2,pj1:pj2,lk1:lk2,1:n_neoclass_mats,ll1:ll2,lm1:lm2,ln1:ln2))
  if (.not.allocated(neoflux2D)) allocate(neoflux2D( 0:nky0-1, 0:nz0-1, 1:n_neoclass_mats, 0:n_spec-1))

  if (mype==0) then
     OPEN(NEOCLASSFILE2D, file=trim(diagdir)//'/neo2D'&
          &//trim(file_extension), form='unformatted',&
          status=filestat, position=filepos,ACCESS="STREAM")
     write(NEOCLASSFILE2D) n_spec
     write(NEOCLASSFILE2D) nky0 
     write(NEOCLASSFILE2D) ly
     write(NEOCLASSFILE2D) nz0
     write(NEOCLASSFILE2D) zval
  endif  

END SUBROUTINE INITIALIZE_ALL_DIAGS_NEOCLASS2D

!!!******************************************************************!!!
  !! sets the integration matrices for the neoclassical fluxes and bootstrap current
  SUBROUTINE set_mats_diag_neoclass(loc_neo_mat,n_nc_mats)
    Implicit None
    integer,intent(in)::n_nc_mats
    real,dimension(pi1:pi2,pj1:pj2,lk1:lk2,1:n_nc_mats,ll1:ll2,lm1:lm2,ln1:ln2),intent(out)::loc_neo_mat
    Integer :: j, k,l,m,n
    real, dimension(pi1:pi2,pj1:pj2,lk1:lk2):: Q_neo_pref,G_neo_pref,P_neo_pref
    DO n=ln1,ln2
       !!Q_neo_pref = 2 K_x n_0j(x0) T_0j(x0)^2/qj Bfield
       !!G_neo_pref = 2 K_x n_0j(x0) T_0j(x0)  /qj Bfield
       !!P_neo_pref = 2 K_x n_0j(x0) T_0j(x0)^2/qj Bfield v_th(x0)
       if (yx_order) then
          G_neo_pref(pi1:pi2,pj1:pj2,:) =  geom%K_j(pi1:pi2,pj1:pj2,lk1:lk2)&
               &*spec(n)%temp*spec(n)%dens&
               &/(geom%Bfield(pi1:pi2,pj1:pj2,lk1:lk2)*spec(n)%charge)!/sqrt(geom%gjj(pi1:pi2,pj1:pj2,:))
          if (norm_flux_projection) G_neo_pref(pi1:pi2,pj1:pj2,:) = &
                &G_neo_pref(pi1:pi2,pj1:pj2,:)/sqrt(geom%gjj(pi1:pi2,pj1:pj2,lk1:lk2))
       else
          G_neo_pref(pi1:pi2,pj1:pj2,:) =  geom%K_i(pi1:pi2,pj1:pj2,lk1:lk2)&
               &*spec(n)%temp*spec(n)%dens&
               &/(geom%Bfield(pi1:pi2,pj1:pj2,lk1:lk2)*spec(n)%charge)!/sqrt(geom%gii(pi1:pi2,pj1:pj2,:))
          if (norm_flux_projection) G_neo_pref(pi1:pi2,pj1:pj2,:) = &
                &G_neo_pref(pi1:pi2,pj1:pj2,:)/sqrt(geom%gii(pi1:pi2,pj1:pj2,lk1:lk2))
       end if

       Q_neo_pref(pi1:pi2,pj1:pj2,lk1:lk2) = G_neo_pref(pi1:pi2,pj1:pj2,lk1:lk2)*spec(n)%temp
       P_neo_pref(pi1:pi2,pj1:pj2,lk1:lk2) = Q_neo_pref(pi1:pi2,pj1:pj2,lk1:lk2)/sqrt(2*spec(n)%temp/spec(n)%mass)

       DO m=lm1,lm2
          DO l=ll1,ll2    
             DO k=lk1,lk2
                DO j=pj1,pj2
                   !mat_20+1/2 mat_01  for G_neo
                   loc_neo_mat(:,j,k,1,l,m,n) = G_neo_pref(:,j,k)&
                                               &*(2.*vp(l)*vp(l)+mu(m)*geom%Bfield(pi1:pi2,j,k))*mat_00(pi1:pi2,j,k,l,m)
                   !mat_40+3/2 mat_21+1/2 mat_02  for Q_neo
                   loc_neo_mat(:,j,k,2,l,m,n) = Q_neo_pref(:,j,k)&
                                               &*(2.* vp(l)*vp(l)*vp(l)*vp(l)&
                                               &+3.*vp(l)*vp(l)*mu(m)*geom%Bfield(pi1:pi2,j,k)&
                                               &+(mu(m)*geom%Bfield(pi1:pi2,j,k))**2)&
                                               &*mat_00(pi1:pi2,j,k,l,m)
                   !1/2 mat_11 + mat_30 for P_neo
                   loc_neo_mat(:,j,k,3,l,m,n) = P_neo_pref(:,j,k)&
                                                &*(vp(l)*mu(m)*geom%Bfield(pi1:pi2,j,k) + 2*vp(l)**3)&
                                                &*mat_00(pi1:pi2,j,k,l,m)
                   !B0*mat_10 essentially (bootstrap current)                            
                   loc_neo_mat(:,j,k,4,l,m,n) = spec(n)%dens*sqrt(2*spec(n)%temp/spec(n)%mass)&
                                               &*geom%Bfield(pi1:pi2,j,k)*vp(l)*mat_00(pi1:pi2,j,k,l,m)
                END DO
             END DO
          END DO
       END DO
    END DO    
  END SUBROUTINE set_mats_diag_neoclass
  
#ifdef WITH_NEO_F0
  !!computes f0 contribution to neoclassical fluxes (parallelized in i,j)
  !!f0 is integrated analytically
  !!the normalization by rhostar makes sense in global computations only..
  !!>todo y-global yx_order, flux surface average
  subroutine get_neofluxf0(loc_neofluxf0)
    implicit none
    real, dimension(pg1gl:pg2gl,0:n_spec-1,0:1),intent(out):: loc_neofluxf0
    Integer :: j, k, n, o
    real, dimension(pi1:pi2,pj1:pj2):: pQneo0,pGneo0

    !! zeroth order fluxes, hkd thesis, only for finite rhostar
    loc_neofluxf0 = 0.0
    if (rhostar.gt.0) then
       do n=ln1,ln2
          do j=pj1,pj2
             pGneo0(pi1:pi2,j)=spec(n)%dens*spec(n)%dens_prof(pi1:pi2)/spec(n)%charge&
                  &*spec(n)%temp*spec(n)%temp_prof(pi1:pi2)/rhostar
             pQneo0(pi1:pi2,j)=2.5*pGneo0(pi1:pi2,j)*spec(n)%temp*spec(n)%temp_prof(pi1:pi2)
          enddo

          !!perform flux surface average of the geometry factors
          !!>todo:y average in yx order?
          !!>todo: mpi sum y in yx order / y global?  right now: ignore pjj:pj2 2nd dimension
          if (p_has_0_mode) then
             if (yx_order) then
                do k=lk1,lk2
                    loc_neofluxf0(pi1:pi2,n,0)= loc_neofluxf0(pi1:pi2,n,0)&
                         &+geom%K_j(pi1:pi2,pj1,k)/sqrt(geom%gjj(pi1:pi2,pj1,k))&
                         &/geom%Bfield(pi1:pi2,pj1,k)*geom%jacobian(pi1:pi2,pj1,k)
                enddo
             else
                do k=lk1,lk2
                    loc_neofluxf0(pi1:pi2,n,0)= loc_neofluxf0(pi1:pi2,n,0)&
                         &+geom%K_i(pi1:pi2,pj1,k)/sqrt(geom%gii(pi1:pi2,pj1,k))&
                         &/geom%Bfield(pi1:pi2,pj1,k)*geom%jacobian(pi1:pi2,pj1,k)
                enddo
             endif
          endif
          

          !!mpi sum z (only my_pez=0 gets result) 
          Call my_sum_to_0_real(loc_neofluxf0(pi1:pi2,n,0),Size(loc_neofluxf0(pi1:pi2,n,0)),mpi_comm_z)

          !!Q and Gamma only differ in prefactor 
          !Q
          loc_neofluxf0(pi1:pi2,n,1)=pQneo0(pi1:pi2,pj1)*loc_neofluxf0(pi1:pi2,n,0)
          !Gamma
          loc_neofluxf0(pi1:pi2,n,0)=pGneo0(pi1:pi2,pj1)*loc_neofluxf0(pi1:pi2,n,0)
          
          !gather all x profile points
          do o=0,1 
             call my_real_gather_to_0(loc_neofluxf0(:,n,o),pi1+1,pi0,px0,mpi_comm_x)
             !normalization for yz average
             loc_neofluxf0(:,n,o) = loc_neofluxf0(:,n,o)*fnorm
          enddo

       enddo !n loop

       !gather species points
       do o=0,1
          call my_real_gather_to_0(loc_neofluxf0(:,:,o),px0*ln1+1,px0*ln0,px0*n_spec,&
               & mpi_comm_spec)
       enddo
    endif    !rhostar.le.0 : nothing
  end subroutine get_neofluxf0
#endif


!!!******************************************************************!!!
  !> neoclass diagnostic: Calculates flux surface averaged velocity space moments:
  !! 1.) G_neo / (c_ref n_ref (rho_ref/Lref)^2)
  !! 2.) Q_neo / (c_ref p_ref (rho_ref/Lref)^2)
  !! 3.) P_neo / (c_ref^2 n_ref m_ref (rho_ref/Lref)^2)
  !! 4.) j_B   / (c_ref n_ref B_ref (rho_ref/Lref)) 
  Subroutine exec_all_diags_neoclass
    real, dimension(0:n_spec-1,1:n_neoclass_mats)::neoflux_scalar
    Integer:: n, o
#ifdef WITH_NEO_F0
    real, dimension(0:n_spec-1,0:1)::neofluxf0_scalar
#endif

    call get_neoflux(neoflux)

    !select scalar flux (drop profile information) for writing into neoclass file
    do o=1,n_neoclass_mats !... <G_neo>, <Q_neo>, <P_neo>, <bootstrap current>
       do n=ln1,ln2
          if (xy_local) then
             !xy_local : select kx=0 mode (x is not parallelized)
             if(p_has_00_mode) neoflux_scalar(n,o)=neoflux(li1,n,o)
          elseif(y_local) then
             !x_global: profile in x: select x0 position
             neoflux_scalar(n,o)=(neoflux(nx0/2,n,o)+neoflux(nx0/2-evenx,n,o))/2.0
          elseif(x_local) then
             !y_global:  select kx=0 mode (x is not parallelized)
             neoflux_scalar(n,o) = neoflux(lg1,n,o)
          else
             if (mype.eq.0) stop "x y global not supported"
          endif
       enddo
    enddo

    !gather species information to 0 process
    do o=1,n_neoclass_mats
       call my_real_gather_to_0(neoflux_scalar(:,o),ln1+1,ln0,n_spec,&
            & mpi_comm_spec)
    enddo

#ifdef WITH_NEO_F0
    call write_neoclassfile(neoflux_scalar,n_neoclass_mats,neofluxf0_scalar)     
#else
    call write_neoclassfile(neoflux_scalar,n_neoclass_mats)     
#endif

  End Subroutine exec_all_diags_neoclass
  
 Subroutine exec_all_diags_neoclass2D 
    real, dimension(0:n_spec-1,1:n_neoclass_mats)::neoflux_scalar

#ifdef WITH_NEO_F0
    real, dimension(0:n_spec-1,0:1)::neofluxf0_scalar
#endif

    call get_neoflux2D(neoflux2D) 
    call write_neoclassfile2D(neoflux2D)

  End Subroutine exec_all_diags_neoclass2D




  !computes the neoclassical fluxes
  Subroutine get_neoflux(neoflux)
    real, dimension(0:nx0-1,0:n_spec-1,1:n_neoclass_mats),intent(out) :: neoflux
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_neoclass_mats,ln1:ln2):: momc_neo
    real, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_neoclass_mats,ln1:ln2):: mom_neo
    Integer:: n, o

    momc_neo=(0.0)

    !\ATTENTION: neoclassical flux contains (numerical) FLR corrections if emfields are considered  
    !(needs fields WITH x and y coordinate: extended memory need, but simpler code)
    Call calc_vsp_moment(n_neoclass_mats,f_,.true.,emfields,neoclass_mat,momc_neo,nc_exb_corr)

    mom_neo(:,:,:,:,:) = Real(momc_neo(:,:,:,:,:))
    
    !-------------- Sum over velocity space  -----------------
    Call my_sum_to_0_real(mom_neo(:,:,:,:,:), Size(mom_neo(:,:,:,:,:)), mpi_comm_vw)
    
    !-------------- Flux surface averaging -------------------
    do o=1,n_neoclass_mats !... <G_neo>, <Q_neo>, <P_neo>, <bootstrap current>
       do n=ln1,ln2
          call flux_surface_average(mom_neo(:,:,:,o,n),.false.,neoflux(lg1:lg2,n,o))
          !gather profile information to 0 process
          if (.not.x_local) call my_real_gather_to_0(neoflux(:,n,o),pi1+1,pi0,nx0,mpi_comm_x)
       enddo
    enddo

  End Subroutine get_neoflux
  

 Subroutine get_neoflux2D(flux2D)
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_neoclass_mats,ln1:ln2):: momc_neo
    real, Dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_neoclass_mats,ln1:ln2):: mom_neo
    real, dimension(0:nky0-1, 0:nz0-1, 1:n_neoclass_mats, 0:n_spec-1),intent(out) :: flux2D
    Integer:: i 
    Integer:: n, o

    momc_neo=(0.0)

    !\ATTENTION: neoclassical flux contains (numerical) FLR corrections if emfields are considered  
    !(needs fields WITH x and y coordinate: extended memory need, but simpler code)
    Call calc_vsp_moment(n_neoclass_mats,f_,.true.,emfields,neoclass_mat,momc_neo,nc_exb_corr)

    mom_neo(:,:,:,:,:) = Real(momc_neo(:,:,:,:,:))
    
    !-------------- Sum over velocity space  -----------------
    Call my_sum_to_0_real(mom_neo(:,:,:,:,:), Size(mom_neo(:,:,:,:,:)), mpi_comm_vw)
    
    !-------------- Flux surface averaging -------------------
    do o=1,n_neoclass_mats !... <G_neo>, <Q_neo>, <P_neo>, <bootstrap current>
       do n=ln1,ln2
          flux2D(li1:li2,lk1:lk2,o, n) = mom_neo(li1:li2,lj1,lk1:lk2,o,n)
          do i=li1, li2
             call my_real_gather_to_0(flux2D(i,:,o,n), lk1+1, lk0, nz0, mpi_comm_z)
          end do 
       enddo
    enddo

  End Subroutine get_neoflux2D
  



  !writes neoclass data into one file containing all species
  subroutine write_neoclassfile(loc_neoflux,n_neoclass_mats,loc_neofluxf0)
    integer::n_neoclass_mats,n
    real, dimension(0:n_spec-1,1:n_neoclass_mats),intent(in) :: loc_neoflux
    real, dimension(0:n_spec-1,0:1)  ,intent(in),optional :: loc_neofluxf0
    real,dimension(1:n_neoclass_mats):: oarr
#ifdef WITH_NEO_F0
    real,dimension(0:1):: oarrf0
#endif

    if (mype == 0) then
       write(NEOCLASSFILE, "(F15.6)") time
       do n=0,n_spec-1
          oarr = loc_neoflux(n,:)
#ifdef WITH_NEO_F0
          oarrf0 = loc_neofluxf0(n,:)
          write(NEOCLASSFILE, "(6ES12.4)") oarr,oarrf0
#else
          write(NEOCLASSFILE, "(4ES12.4)") oarr
#endif
       enddo
      call flush(NEOCLASSFILE)
    endif
  end subroutine write_neoclassfile


subroutine write_neoclassfile2D(flux2D) 
!writes 2D neoclassical flux
    real, dimension(0:nky0-1, 0:nz0-1, 1:n_neoclass_mats, 0:n_spec-1),intent(in) :: flux2D
    integer:: j, k, n, o

    if (mype == 0) then
       write(NEOCLASSFILE2D) time

       DO n = 0, n_spec-1
          DO o =1,n_neoclass_mats
             DO j = 0, nky0-1
                DO k = 0, nz0-1
                   WRITE(NEOCLASSFILE2D) flux2D(j, k , o, n)
                END DO
             END DO
          END DO
       END DO
       
       call flush(NEOCLASSFILE2D)
    endif
end subroutine write_neoclassfile2D


  Subroutine finalize_all_diags_neoclass
    Implicit none
    Integer :: ierr
    
    If (mype==0) CLOSE(NEOCLASSFILE)
    Deallocate(neoclass_mat)

    if (par_in_dir.ne.'skip_parfile') then
       DEALLOCATE(neoflux)
    else !make neoflux available on all cores for external usage (e.g., trinity)
       call mpi_bcast(neoflux,SIZE(neoflux),MPI_REAL_TYPE,&
            &0,MY_MPI_COMM_WORLD,ierr)
    endif
#ifdef WITH_NEO_F0 
    Deallocate(fnorm)
    if (par_in_dir.ne.'skip_parfile') then
       DEALLOCATE(neofluxf0)
    else !make neoflux available on all cores for external usage (e.g., trinity)
       call mpi_bcast(neofluxf0,SIZE(neofluxf0),MPI_REAL_TYPE,&
            &0,MY_MPI_COMM_WORLD,ierr)
    endif
#endif

  End Subroutine finalize_all_diags_neoclass

  Subroutine finalize_all_diags_neoclass2D 

    Implicit none
    deallocate(neoflux2D)
    If (mype==0) CLOSE(NEOCLASSFILE2D) 
    
  End Subroutine finalize_all_diags_neoclass2D





End Module diagnostics_neoclass
