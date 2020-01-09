#include "redef.h"

module exchange_z_ff
  use discretization
  use par_in
  use par_other
  use fourier
  use communications

  implicit none
  
  public:: init_parallel_boundary_ff, parallel_boundary_ff, finalize_parallel_boundary_ff  
  public:: exchange_5df_ff_1, exchange_5df_ff_2, initialize_exch_z_cyc, finalize_exch_z_cyc
  public:: delete_entry_lower, delete_entry_upper, pb_xshift, pb_phase, ubexc
  private
  
  integer, dimension(:), allocatable:: pb_xshift !parallel b.c. x shift
  complex, dimension(:), allocatable:: pb_phase !phase factor
  integer, dimension(:), allocatable:: pb_lind, pb_uind, pb_lind2, pb_uind2
  complex, dimension(:,:,:), allocatable:: zeroarray

  logical, dimension(:,:), allocatable:: delete_entry_lower, delete_entry_upper
  logical:: external_zbc
  integer:: ubexc
  integer:: init_status_exch_z_cyc = 0

  integer:: vectype_zexc

contains
  !****m* boundaries/init_parallel_boundary
  ! DESCRIPTION
  ! Initializes the boundary conditions.
  ! This function has to be called first before any of the z exchange
  ! functions.
  ! 
  ! SYNOPSIS
  subroutine init_parallel_boundary_ff(q0, shat)
    implicit none
    real, intent(in) :: q0, shat

    integer :: ubglob, i, j
    integer :: py_num

    
    ! ubglob is the highest ky-mode which will get nonzero boundaries in exchange_z, i.e. 
    ! only modes lh1:ubexc have to be calculated. 
    ! ubexc is also a marker: -1 means no connection, -2 means periodic connection
    !                         -3 means no exchange at all
    if ((abs(shat).gt.1e-3).and.(.not.only_zonal)) then
       IF(y_local) THEN
          IF (nx0.GT.1) THEN
             ! up to which ky mode do we make a connection? It either the last
             ! ky mode in the system (nky0-1+ky0_ind) or it is the last ky mode
             ! for which a connection exists (that means the matching x mode is still in the
             ! system. Whatever leads to a lower index is chosen.
             ubglob = min(nky0-1,(nx0-evenx)/abs(nexc)-ky0_ind) 
             py_num = ubglob/lh0
             if (my_pey.lt.py_num) ubexc=lh2    ! all lower processes are involved
             if (my_pey.eq.py_num) ubexc=ubglob ! actual process only up to ubglob
             if (my_pey.gt.py_num) ubexc=-1     ! higher processes are marked as not involved
             IF (ubexc.GE.lh1) THEN
                allocate(pb_xshift(lh1:ubexc),pb_lind(lh1:ubexc),pb_uind(lh1:ubexc),pb_phase(lh1:ubexc))
                do j=lh1,ubexc
                   pb_xshift(j)=nexc*(j+ky0_ind)
                enddo
                CALL set_pb_indices
             END IF
          ELSE !(nx0.eq.1)
             !for 1 x-mode linear runs:
             ubexc=-1
          END IF
       ELSE
          !global y
          if (nx0.gt.1) then
             ubglob = nky0-1 
          else
             !only ky=0 mode has nonzero boundary condition
             ubglob = 0 
             !this will set ubexc = 0
             !note that values of ubglob other than 0 or nky0-1 make no sense
             !because negative ky are included after fourier transform.
          endif
          py_num = ubglob/lh0
          if (my_pey.lt.py_num) ubexc=lh2    ! all lower processes are involved
          if (my_pey.eq.py_num) ubexc=ubglob ! actual process only up to ubglob
          if (my_pey.gt.py_num) ubexc=-1     ! higher processes are marked as not involved
          IF (ubexc.GE.lh1) THEN
             allocate(pb_xshift(lh1:lh2),pb_lind(lh1:lh2),pb_uind(lh1:lh2),pb_phase(lh1:lh2))
             allocate(pb_lind2(lh1:lh2),pb_uind2(lh1:lh2))
             do i=lh1,nky0/2
                pb_xshift(i)=nexc*i
             enddo
             do i=nky0/2+1,lh2
                pb_xshift(i)=-nexc*(nky0-i)
             enddo
             CALL set_pb_indices
          END IF
       END IF
    else
       if (.not.only_zonal.and.(n0_global.ne.0 .and. n0_global.ne.-1111)) then 
          !shat=0 : n0_global phaseshift remains (if n0_global is defined)
          ubexc = lh2
          ALLOCATE(pb_xshift(lh1:lh2),pb_lind(lh1:lh2),pb_uind(lh1:lh2),pb_phase(lh1:lh2))
          if (.not.y_local) allocate(pb_lind2(lh1:lh2),pb_uind2(lh1:lh2))
          pb_xshift=0
          CALL set_pb_indices
       else
          !zonal modes or shat==0 and n0_global==0:
          ubexc=-2  !periodic boundary condition
       endif
    end if

    if (nzb.le.0) then
      ubexc=-3  !no exchange at all
    end if
    
    !internal z boundaries (ghost cells)
    if ((my_pez.eq.0).or.(my_pez.eq.(n_procs_z-1))) then
       external_zbc=.true.
    else
       external_zbc=.false.
    endif

    if (allocated(pb_phase)) then
       pb_phase= (-1)**(pb_xshift) !this factor has to be retained in products like J0*phi!

       if ((n0_global.ne.0.and.n0_global.ne.-1111).and.(ubexc.ge.lh1)) then
          if(y_local) then
             do j=lh1,ubexc
                pb_phase(j) = pb_phase(j)*exp(-imag*(j+ky0_ind)*2.0*pi*n0_global*q0)
             end do
          else
             do j=lh1,nky0/2
                pb_phase(j) = pb_phase(j)*exp(-imag*j*2.0*pi*n0_global*q0)
             end do
             do j=nky0/2+1,lh2
                pb_phase(j)= pb_phase(j)*exp(imag*(nky0-j)*2.0*pi*n0_global*q0)
             end do
          end if
       end if
    end if

    !allocate array with zeros for later usage in e.g. exchange_z
    ALLOCATE(zeroarray(li1:li2,lj1:lj2,1:nzb))
    zeroarray = (0.,0.)
    
    !print some information
    if ((mype.eq.0).and.(print_ini_msg)) then 
       write(*,"(a,i3)") 'N for parallel boundary condition: ',nexc
       if (ubexc.eq.-2) write(*,"(a)") "Using periodic boundary conditions in parallel direction"
    endif

    allocate(delete_entry_lower(li1:li2,lj1:lj2), delete_entry_upper(li1:li2,lj1:lj2))
    if (ubexc.gt.-2) then
       delete_entry_lower=.true.
       delete_entry_upper=.true.
       do j=lh1,ubexc
          do i=pb_lind(j), pb_uind(j)
             if (yx_order) then
                if(y_local) then
                   delete_entry_lower(j,lg1+mod(i-lg1+pb_xshift(j),nx0))=.false.
                   delete_entry_upper(j,lg1+mod(i-lg1,nx0))=.false.
                end if
             else
                delete_entry_lower(lg1+mod(i-lg1+pb_xshift(j),nx0),j)=.false.
                delete_entry_upper(lg1+mod(i-lg1,nx0),j)=.false.
             endif
          enddo
       enddo
    else
       !periodic
       delete_entry_lower=.false.
       delete_entry_upper=.false.
    endif

  end subroutine init_parallel_boundary_ff

  subroutine exchange_5df_ff_1(u)
    implicit none
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),Intent(Inout):: u
    integer:: m,n
 
    do n=ln1,ln2
       do m=lm1,lm2
          call parallel_boundary_ff(u(:,:,:,:,m,n),ll0,nvb)
       enddo
    enddo

  end subroutine exchange_5df_ff_1

  !>Exchange ghost cells and apply boundary condition in z direction
  !!\param u at least 3D array where the z direction has to be the third index
  !!\param l0 number of elements on indices > 3
  !!\param b0 number of boundary cells on indices > 0 which can be skipped
  Subroutine parallel_boundary_ff(u, l0, b0)
    implicit none
    Integer, Intent(In):: l0, b0
    COMPLEX,INTENT(INOUT):: u(li1:li2, lj1:lj2, lbz:ubz, -b0:l0+b0-1)
    Integer:: i, j, k, l, tag, dest_pes, recv_pes
    Integer:: stat(MPI_STATUS_SIZE), ierr
    Complex, dimension(li1:li2, lj1:lj2, 1:nzb):: tmps_r
    integer:: kls,klb,krs,krb,s_x,minus_y
    complex,dimension(li1:li2,0:lj0-1):: temp1, temp2


    if (ubexc .eq. -3) return !nzb=0, no points to exchange

    If (n_procs_z.Eq.1) Then       
       ! no MPI parallelization over z coordinate
       Do l = 0,l0-1
          if (ubexc.eq.-2) then
             call ccopy(lij0*nzb,u(li1,lj1,lk2-nzb+1,l),1,u(li1,lj1,lbz,l),1)
             call ccopy(lij0*nzb,u(li1,lj1,lk1,l),1,u(li1,lj1,lk2+1,l),1)
          else 
             if (y_local) then
                call ccopy(lij0*nzb,zeroarray(li1,lj1,1),1,u(li1,lj1,lbz,l),1)
                call ccopy(lij0*nzb,zeroarray(li1,lj1,1),1,u(li1,lj1,lk2+1,l),1)
             end if
          endif
          if (ubexc.ge.0) then
             do k=1,nzb  
                !compute left and right (l,r) k index in ghost cell and simulation domain (g,s)
                kls=lk1+k-1
                klb=lbz+k-1
                krs=lk2-nzb+k
                krb=lk2+k
                ! shift of the x-modes according to boundary condition
                ! kx(i) with i > (nx0-1)/2 are assumed to be 0 
                if (yx_order) then
                   if (y_local) then
                      do i=li1,ubexc
                         do j=pb_lind(i), pb_uind(i)
                            !f(z-lz) = f(z)*conjg(phase_factor)
                            u(i,lj1+mod(j-lj1+pb_xshift(i),nx0),klb,l)=u(i,lj1+mod(j-lj1,nx0),krs,l)*CONJG(pb_phase(i))
                            !f(z+lz) = f(z)*phase_factor
                            u(i,lj1+mod(j-lj1,nx0),krb,l)=u(i,lj1+mod(j-lj1+pb_xshift(i),nx0),kls,l)*pb_phase(i)
                         enddo
                      enddo
                   else !global y:
                      !f(z+lz) = f(z)*phase_factor
                      call to_fourier_boundary(u(:,:,kls,l),temp1)     
                      temp2=(0.,0.)
                      do i=li1,ubexc
                         if ((even_ni0.eq.1).and.(i.eq.nky0/2)) cycle !no conn. for highest (unbalanced) ky mode
                         minus_y=mod(2*nky0-i,nky0)
                         do j=pb_lind(i), pb_uind(i)
                            s_x=j+pb_xshift(i)

                            if (s_x.ge.0) then
                               temp2(i,j)=temp1(i,s_x)*pb_phase(i)
                            else
                               !kx mode is given implicitly via reality condition
                               temp2(i,j)=conjg(temp1(minus_y,-s_x))*pb_phase(i)
                            end if
                         end do
                      end do
                      call to_real_boundary(temp2,u(:,:,krb,l))

                      !f(z-lz) = f(z)*conjg(phase_factor)
                      call to_fourier_boundary(u(:,:,krs,l),temp1)
                      temp2=(0.,0.)
                      do i=li1,ubexc
                         if ((even_ni0.eq.1).and.(i.eq.nky0/2)) cycle  !no conn. for highest (unbalanced) ky mode
                         minus_y=mod(2*nky0-i,nky0)
                         do j=pb_lind2(i), pb_uind2(i)
                            s_x=j-pb_xshift(i)
                            if (s_x.ge.0) then
                               temp2(i,j)=temp1(i,s_x)*conjg(pb_phase(i))
                            else
                               !kx mode is given implicitly via reality condition
                               temp2(i,j)=conjg(temp1(minus_y,-s_x))*conjg(pb_phase(i))
                            end if
                         end do
                      end do
                      call to_real_boundary(temp2,u(:,:,klb,l))
                   end if
                else   !not yx_order==F
                   do j=lj1,ubexc
                      do i=pb_lind(j), pb_uind(j)
                         !f(z-lz) = f(z)*conjg(phase_factor)
                         u(li1+mod(i-li1+pb_xshift(j),nx0),j,klb,l)=u(li1+mod(i-li1,nx0),j,krs,l)*CONJG(pb_phase(j))
                         !f(z+lz) = f(z)*phase_factor
                         u(li1+mod(i-li1,nx0),j,krb,l)=u(li1+mod(i-li1+pb_xshift(j),nx0),j,kls,l)*pb_phase(j)
                      enddo
                   enddo
                end if
             enddo
          endif
       enddo
    Else 
          !
          !     Send Boundary to next pe in z-Direction:
          !
       Do l = 0, l0-1
          dest_pes = mod(my_pez+1,n_procs_z)
          recv_pes = mod(my_pez-1+n_procs_z,n_procs_z)
          tag = 333+l
          Call mpi_sendrecv(&
               u(li1,lj1,lk2-nzb+1,l), lij0*nzb, &
               &MPI_COMPLEX_TYPE, dest_pes, tag,&
               &tmps_r(li1,lj1,1), lij0*nzb, MPI_COMPLEX_TYPE, recv_pes, tag,&
               &mpi_comm_z, stat, ierr)
          If (my_pez > 0) Then
             u(:,:,lbz:lk1-1,l) = tmps_r
          Else
             If (ubexc.eq.-1) Then
                u(:,:,lbz:lk1-1,l)=(0.,0.)
             elseif (ubexc.eq.-2) then
                u(:,:,lbz:lk1-1,l) = tmps_r
             else
                Do k = 1, nzb
                   klb=lbz+k-1
                   if (yx_order) then
                      if (y_local) then
                         u(:,:,klb,l)=(0.,0.) 
                         do i=li1,ubexc
                            do j=pb_lind(i), pb_uind(i)
                               u(i,lj1+mod(j-lj1+pb_xshift(i),nx0),klb,l)=tmps_r(i,lj1+mod(j-lj1,nx0),k)*conjg(pb_phase(i))
                            enddo
                         enddo
                      else
                         call to_fourier_boundary(tmps_r(:,:,k),temp1)
                         temp2=(0.,0.)
                         do i=li1,ubexc
                            if ((even_ni0.eq.1).and.(i.eq.nky0/2)) cycle  !no conn. for highest (unbalanced) ky mode
                            minus_y=mod(2*nky0-i,nky0)
                            do j=pb_lind2(i), pb_uind2(i)
                               s_x=j-pb_xshift(i)
                               if (s_x.ge.0) then
                                  temp2(i,j)=temp1(i,s_x)*conjg(pb_phase(i))
                               else
                                  !kx mode is given implicitly via reality condition
                                  temp2(i,j)=conjg(temp1(minus_y,-s_x))*conjg(pb_phase(i))
                               end if
                            end do
                         end do
                         call to_real_boundary(temp2,u(:,:,klb,l))
                      end if
                   else  
                      u(:,:,klb,l)=(0.,0.)                
                      do j=lj1,ubexc
                         do i=pb_lind(j), pb_uind(j)
                            u(li1+mod(i-li1+pb_xshift(j),nx0),j,klb,l)=tmps_r(li1+mod(i-li1,nx0),j,k)*conjg(pb_phase(j))
                         enddo
                      enddo
                   end if
                Enddo
             endif
          Endif
          !
          !     Send Boundary to previous pe in z-Direction:
          !
          dest_pes = mod(my_pez-1+n_procs_z,n_procs_z)
          recv_pes = mod(my_pez+1,n_procs_z)
          tag = 340+l0+l
          Call mpi_sendrecv(&
               u(li1,lj1,lk1,l), Size(u(:,:,lk1:lk1+nzb-1,l)), MPI_COMPLEX_TYPE, dest_pes, tag,&
               tmps_r(li1,lj1,1), Size(tmps_r), MPI_COMPLEX_TYPE, recv_pes, tag,&
               mpi_comm_z, stat, ierr)
          If (my_pez < n_procs_z-1) Then
             u(:,:,lk2+1:ubz,l) = tmps_r
          Else
             if (ubexc.eq.-1) then
                u(:,:,lk2+1:ubz,l)=cmplx(0,0)
             elseif (ubexc.eq.-2) then
                u(:,:,lk2+1:ubz,l) = tmps_r
             else
                Do k = 1, nzb
                   krb=lk2+k
                   if (yx_order) then
                      if (y_local) then
                         u(:,:,krb,l)=(0.,0.)
                         do i=li1,ubexc
                            do j=pb_lind(i), pb_uind(i)
                               u(i,lj1+mod(j-lj1,nx0),krb,l)= tmps_r(i,lj1+mod(j-lj1+pb_xshift(i),nx0),k)*pb_phase(i)
                            enddo
                         enddo
                      else
                         !global y
                         call to_fourier_boundary(tmps_r(:,:,k),temp1)     
                         temp2=(0.,0.)
                         do i=li1,ubexc
                            if ((even_ni0.eq.1).and.(i.eq.nky0/2)) cycle  !no conn. for highest (unbalanced) ky mode
                            do j=pb_lind(i), pb_uind(i)
                               s_x=j+pb_xshift(i)    
                               if (s_x.ge.0) then
                                  temp2(i,j)=temp1(i,s_x)*pb_phase(i)
                               else
                                  !kx mode is given implicitly via reality condition
                                  minus_y=mod(2*nky0-i,nky0)
                                  temp2(i,j)=conjg(temp1(minus_y,-s_x))*pb_phase(i)
                               end if
                            end do
                         end do
                         call to_real_boundary(temp2,u(:,:,krb,l))
                      end if
                   else
                      u(:,:,krb,l)=(0.,0.)
                      do j=lj1,ubexc
                         do i=pb_lind(j), pb_uind(j)
                            u(li1+mod(i-li1,nx0),j,krb,l)= tmps_r(li1+mod(i-li1+pb_xshift(j),nx0),j,k)*pb_phase(j)
                         enddo
                      enddo
                   end if
                Enddo
             endif
          Endif
       enddo
    End If

  End Subroutine parallel_boundary_ff


!-----------------------------------SECOND VERSION---------------------------------------------------

  subroutine exchange_5df_ff_2(u)
    implicit none
    COMPLEX,DIMENSION(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),INTENT(inout):: u
    Integer:: i, j, k, l, m, n, kls, klb, krs, krb

    if (ubexc .eq. -3) return !nzb=0, no points to exchange

    if (n_procs_z.gt.1) then
       call exch_z_cyc(u)
    endif

    if (external_zbc) then
       !shifts in kx
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                if ((ubexc.eq.-2).and.(n_procs_z.eq.1)) then 
                   !periodic
                   call ccopy(li0*lj0*nzb,u(li1,lj1,lk2-nzb+1,l,m,n),1,u(li1,lj1,lbz,l,m,n),1)
                   call ccopy(li0*lj0*nzb,u(li1,lj1,lk1,l,m,n),1,u(li1,lj1,lk2+1,l,m,n),1)
                elseif (ubexc.eq.-1) then
                   !zero boundary condition	
                   if(my_pez.eq.0) call ccopy(li0*lj0*nzb,zeroarray(li1,lj1,1),1,u(li1,lj1,lbz,l,m,n),1)
                   if(my_pez.eq.(n_procs_z-1)) call ccopy(li0*lj0*nzb,zeroarray(li1,lj1,1),1,u(li1,lj1,lk2+1,l,m,n),1)
                elseif (ubexc.ge.0) then
                   do k=1,nzb  
                      kls=lk1+k-1
                      klb=lbz+k-1
                      krs=lk2-nzb+k
                      krb=lk2+k
                      if (n_procs_z.eq.1) then
                         call ccopy(li0*(ubexc-lj1+1),u(li1,lj1,krs,l,m,n),1,u(li1,lj1,klb,l,m,n),1)
                         call ccopy(li0*(ubexc-lj1+1),u(li1,lj1,kls,l,m,n),1,u(li1,lj1,krb,l,m,n),1) 
                      endif
                      ! shift of the x-modes according to boundary condition
                      ! kx(i) with i > (nx0-1)/2 are assumed to be 0 
                      if (yx_order) then
                         if (y_local) then
                            if (my_pez.eq.0) then
                               do i=li1,ubexc
                                  u(i,:,klb,l,m,n)=cshift(u(i,:,klb,l,m,n),-pb_xshift(i))*conjg(pb_phase(i))
                               enddo
                               where (delete_entry_lower)
                                  u(:,:,klb,l,m,n)=(0.,0.)
                               end where
                            end if
                            if (my_pez.eq.(n_procs_z-1)) then
                               do i=li1,ubexc
                                  u(i,:,krb,l,m,n)=cshift(u(i,:,krb,l,m,n),pb_xshift(i))*pb_phase(i)
                               enddo
                               where (delete_entry_upper)
                                  u(:,:,krb,l,m,n)=(0.,0.)
                               end where
                            endif
                         else
                            stop 'perf_vec(7)=2 is not implemented yet for global y'
                         end if
                      else
                         if (my_pez.eq.0) then
                            do j=lj1,ubexc
                               u(:,j,klb,l,m,n)=cshift(u(:,j,klb,l,m,n),-pb_xshift(j))*conjg(pb_phase(j))
                            enddo
                            where (delete_entry_lower)
                               u(:,:,klb,l,m,n)=(0.,0.)
                            end where
                         end if
                         if(my_pez.eq.(n_procs_z-1)) then
                            do j=lj1,ubexc
                               u(:,j,krb,l,m,n)=cshift(u(:,j,krb,l,m,n),pb_xshift(j))*pb_phase(j)
                            enddo
                            where (delete_entry_upper)
                               u(:,:,krb,l,m,n)=(0.,0.)
                            endwhere
                         endif
                      endif
                   enddo
                endif
             enddo
          enddo
       enddo
    endif

  End Subroutine exchange_5df_ff_2


  subroutine finalize_parallel_boundary_ff
    implicit none

    If (allocated(pb_xshift)) deallocate(pb_xshift)
    If (allocated(pb_lind)) deallocate(pb_lind)
    If (allocated(pb_uind)) deallocate(pb_uind)
    If (allocated(pb_lind2)) deallocate(pb_lind2)
    If (allocated(pb_uind2)) deallocate(pb_uind2)
    If (allocated(pb_phase)) deallocate(pb_phase)
    
    deallocate(zeroarray)
    deallocate(delete_entry_lower, delete_entry_upper)

  end subroutine finalize_parallel_boundary_ff


  subroutine set_pb_indices
    integer:: i

    if (y_local) then
       if(nexc.gt.0) then
          pb_lind=lkx
          pb_uind=lkx+nx0-1-evenx-pb_xshift
       else
          pb_lind=lkx-pb_xshift
          pb_uind=lkx+nx0-1-evenx
       end if
    else
       do i=li1,li2
          if(pb_xshift(i).gt.0) then
             pb_lind(i)=0
             pb_uind(i)=nx0-1-pb_xshift(i)
             pb_lind2(i)=max(-(nx0-1)+pb_xshift(i),0)
             pb_uind2(i)=nx0-1
          else
             pb_lind(i)=max(-(nx0-1)-pb_xshift(i),0)
             pb_uind(i)=nx0-1
             pb_lind2(i)=0
             pb_uind2(i)=nx0-1+pb_xshift(i)
          end if
       end do
    end if
  end subroutine set_pb_indices

  !auxiliary routines used by local and global version
  subroutine initialize_exch_z_cyc
    integer:: bl_num, bl_len, bl_str, ierr

    if (init_status_exch_z_cyc==1) return

    bl_num=lv0*lw0*ln0-2*lv0*nwb-2*nvb
    bl_len=li0*ly0*nzb
    bl_str=li0*ly0*lz0

    call mpi_type_vector(bl_num,bl_len,bl_str,MPI_COMPLEX_TYPE,vectype_zexc,ierr)
    call mpi_type_commit(vectype_zexc,ierr)
    init_status_exch_z_cyc = 1

  end subroutine initialize_exch_z_cyc

  subroutine exch_z_cyc(u)
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(inout):: u
    integer:: tag, dest_pez, recv_pez, stat(MPI_STATUS_SIZE), ierr

    PERFON('mpi_z_exc')
    dest_pez = mod(my_pez-1+n_procs_z,n_procs_z)
    recv_pez = mod(my_pez+1,n_procs_z)
    tag=345
    Call mpi_sendrecv(&
         & u(li1,lj1,lk1,ll1,lm1,ln1), 1, vectype_zexc, dest_pez, tag,&
         & u(li1,lj1,lk2+1,ll1,lm1,ln1), 1, vectype_zexc, recv_pez, tag,&
         & mpi_comm_z, stat, ierr)
    tag=444
    
    dest_pez = mod(my_pez+1,n_procs_z)
    recv_pez = mod(my_pez-1+n_procs_z,n_procs_z)
    Call mpi_sendrecv(&
         & u(li1,lj1,lk2-nzb+1,ll1,lm1,ln1), 1, vectype_zexc, dest_pez, tag,&
         & u(li1,lj1,lbz,ll1,lm1,ln1), 1, vectype_zexc, recv_pez, tag,&
         & mpi_comm_z, stat, ierr)
    
    PERFOFF

  end subroutine exch_z_cyc

  Subroutine finalize_exch_z_cyc
    integer :: ierr

    call mpi_type_free(vectype_zexc,ierr)
    if (ierr.ne.0) stop 'mpi_type_free failed in finalize_exch_z_cyc'

    init_status_exch_z_cyc = 0

  End Subroutine finalize_exch_z_cyc


!-----------------------------Exchanges required for equilibrium quantities---------------------------------



end module exchange_z_ff
