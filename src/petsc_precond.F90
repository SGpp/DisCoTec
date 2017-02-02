#include "redef.h"
#include "switches.h"
module petsc_precond
  use par_other, only: precond_approx, prec, dt, p_has_0_mode
  use par_in
  use eigen_parameters
  use discretization
  use calc_rhs
  use communications, only: MY_MPI_COMM_WORLD, get_systime
  use geometry,only: shat
#ifdef WITHSLEPC 
  use petsc_aux
  use petscmat
  use petscpc
#endif

  implicit none

  public:: pc_type, pc_sub_it, pc_factor_level, sub_pc_type, pc_blocks, order_rcm, ksp_prec
  public:: check_petsc_precond

#ifdef WITHSLEPC
  public :: initialize_L_g, finalize_L_g
  public :: pc_obj, set_pc_obj
  public :: L_g_mat, L_g_initialized
#endif

  private

#ifdef WITHSLEPC 
#include "finclude/petscvecdef.h"   
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "petscversion.h"

  Mat L_g_mat
  PC pc_obj
  PetscInt n_local,first_local
#endif

  integer:: pc_sub_it=500
  integer:: pc_factor_level=-1
  integer:: pc_blocks=0
  character(len=20) :: pc_type='default', sub_pc_type='default'
  
  logical:: L_g_initialized
  integer:: x_np, v_np, w_np, z_np, x_bl_ub, v_bl_ub, w_bl_ub, n_dep
  integer,dimension(:),allocatable:: kxzcol, kxzoff
  complex:: pre_shift
  logical:: periodic_z=.true., order_rcm=.false.
  real:: ksp_prec=1e-6

  character(len=50):: runtime_option_string

contains

  subroutine check_petsc_precond

    if (pc_type.eq.'default') then
       select case (comp_type)
       case('EV')
          select case (which_ev)
          case('jd','gd','shift_invert','jd_as')
             pc_type='asm'
          case default
             pc_type='none'
          end select
       case('NC')
          pc_type='asm'
       case('IV')
          if(precomp_nc)then
             pc_type='asm'
          else
             pc_type='none'
          end if
       end select
    end if
    
    if (sub_pc_type.eq.'default') sub_pc_type='lu'

    select case (comp_type)
    case('EV')
       select case (which_ev)
       case('jd','gd','jd_as')
          pre_shift=ev_shift
       case('shift_invert')
          pre_shift=(30.,0.)
       case default
          !IE1p
          pre_shift=ev_shift
       end select
    case('NC','IV')
       pre_shift=0.
    end select

  end subroutine check_petsc_precond

#ifdef WITHSLEPC
 
  subroutine set_pc_obj
    call set_pc_defaults
    call init_pc_obj
  end subroutine set_pc_obj

  subroutine set_pc_defaults
    if (pc_blocks.eq.0) then
       pc_blocks=ln0
#ifdef prec_coll
       if (collision_op.eq.'none') pc_blocks=pc_blocks*lm0
#else
       pc_blocks=pc_blocks*lm0
#endif
    end if

    if(pc_factor_level.eq.-1) then
       !set resolution dependend default
       if (comp_type.eq.'EV') then
          if(xy_local) then
             pc_factor_level=2+ceiling(sqrt(real(ni0*nz0*nv0))/4)
          else 
             pc_factor_level=2+ceiling(sqrt(real(nz0*nv0))/4)
          end if
       else
          if(xy_local) then
             pc_factor_level=2+ceiling(sqrt(real(nv0))*real(nz0)/8)
          else
             pc_factor_level=2+ceiling(sqrt(real(nv0*ni0*nz0)/8))
          end if
       end if
       
    end if

  end subroutine set_pc_defaults

  subroutine init_pc_obj
    KSP sub_ksp_obj(0:pc_blocks-1)
    integer,dimension(pc_blocks):: block_sizes

    select case(pc_type)
    case('cl')
       call PCSetFromOptions(pc_obj,globerr)
       call PCSetUp(pc_obj,globerr)
    case('none')
       call PCSetType(pc_obj,PCNONE,globerr)
       call PCSetUp(pc_obj,globerr)
    case('jacobi')
       call PCSetType(pc_obj,PCJACOBI,globerr)
    case('bjacobi')
       if(mype.eq.0) write(*,'(A,I3,A)') 'using block Jacobi preconditioning with ', pc_blocks, ' blocks per processor'
       call PCSetType(pc_obj,PCBJACOBI,globerr)

       !define block sizes
       block_sizes=vlen/pc_blocks
       block_sizes(pc_blocks)=vlen-sum(block_sizes(1:pc_blocks-1))

       call PCBJacobiSetLocalBlocks(pc_obj,pc_blocks,block_sizes,globerr)
       call PCSetFromOptions(pc_obj,globerr)
       call PCSetUP(pc_obj,globerr)
       call PCBJacobiGetSubKSP(pc_obj,n_local,first_local,sub_ksp_obj,globerr)
       call set_sub_ksp(sub_ksp_obj)

    case('asm')
       if(mype.eq.0) write(*,'(A,I3,A)') 'using additive Schwarz preconditioning with ', pc_blocks, ' blocks per processor'

       call PCSetType(pc_obj,PCASM,globerr)
       call PCASMSetOverlap(pc_obj,2,globerr)
       if (comp_type.eq.'NC') then
          write(runtime_option_string,'(a,I5)') '-pc_asm_blocks ',pc_blocks*n_procs_sim 
       else
          write(runtime_option_string,'(a,I5)') '-st_pc_asm_blocks ',pc_blocks*n_procs_sim 
       end if
       call PetscOptionsInsertString(trim(runtime_option_string),globerr)
       call PCSetFromOptions(pc_obj,globerr)
       call PCASMSetType(pc_obj,PC_ASM_RESTRICT,globerr)
       call PCSetUP(pc_obj,globerr)
       call PCASMGetSubKSP(pc_obj,n_local,first_local,sub_ksp_obj,globerr)
       call set_sub_ksp(sub_ksp_obj)

#if PETSC_HAVE_PARMS
    case('parms_bj')
       if(mype.eq.0) write(*,'(A)') 'using parms Block Jacobi preconditioning'
       call PCSetType(pc_obj,PCPARMS,globerr)
       call PCPARMSSetGlobal(pc_obj,PC_PARMS_GLOBAL_BJ,globerr)
       call set_parms_local

   case('parms_ras')
       if(mype.eq.0) write(*,'(A)') 'using parms restricted additive Schwarz preconditioning'
       call PCSetType(pc_obj,PCPARMS,globerr)
       call PCPARMSSetGlobal(pc_obj,PC_PARMS_GLOBAL_RAS,globerr)
       call set_parms_local

       !case('parms_schur') would be the most interesting, however, even with a lot of tuning it is not
       !always possible to find a working setting for each test case, let alone a robust setting for
       !all cases.
#endif

    case default
       stop 'pc_type is not valid'
    end select
    CHKERRQ(globerr)    

  end subroutine init_pc_obj

#if PETSC_HAVE_PARMS
  subroutine set_parms_local

    select case (sub_pc_type)
    case('parms_ilut')
       call PCPARMSSetLocal(pc_obj,PC_PARMS_LOCAL_ILUT,globerr)
       call PCPARMSSetFill(pc_obj,10000,0,0,globerr)
       if (comp_type.eq.'NC') then
          write(runtime_option_string,'(a,es9.2)') '-pc_parms_droptol_factors ',ksp_prec 
          call PetscOptionsInsertString(trim(runtime_option_string),globerr)
       else
          write(runtime_option_string,'(a,es9.2)') '-st_pc_parms_droptol_factors ',10.*ev_prec 
          call PetscOptionsInsertString(trim(runtime_option_string),globerr)
       end if

    case('parms_arms')
       call PCPARMSSetLocal(pc_obj,PC_PARMS_LOCAL_ARMS,globerr)
       !I have only managed to find a robust parameter setting for parms_levels=0. 
       !For more levels, the other droptol/fill values have to be specified as well.
       call PCPARMSSetFill(pc_obj,0,0,10000,globerr)
       !tolerances have to be lowered for nlevels>0
       if (comp_type.eq.'NC') then
          write(runtime_option_string,'(a,es9.2)') '-pc_parms_droptol_last_schur ',ksp_prec 
          call PetscOptionsInsertString(trim(runtime_option_string),globerr) 
          call PetscOptionsInsertString('-pc_parms_levels 0',globerr)
       else
          write(runtime_option_string,'(a,es9.2)') '-st_pc_parms_droptol_last_schur ',100.*ev_prec 
          call PetscOptionsInsertString(trim(runtime_option_string),globerr)
          call PetscOptionsInsertString('-st_pc_parms_levels 0',globerr)
       end if

    case default
       stop 'please select a valid local parms preconditioner for sub_pc_type'
    end select
    call PCSetFromOptions(pc_obj,globerr)
    call PCSetUP(pc_obj,globerr)

  end subroutine set_parms_local
#endif


  subroutine set_sub_ksp(sub_ksp_obj)
    KSP sub_ksp_obj(0:n_local-1)
    KSP block_ksp_obj
    PC block_pc_obj
    integer:: i

    PERFON('block_pc_setup')
    do i = 0,n_local-1
       block_ksp_obj = sub_ksp_obj(i)
#if (PETSC_VERSION_MAJOR>3) || (PETSC_VERSION_MINOR>4)
       call KSPSetTolerances(block_ksp_obj,PETSC_DEFAULT_REAL,&
            &PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,pc_sub_it,globerr)
#else
       call KSPSetTolerances(block_ksp_obj,PETSC_DEFAULT_DOUBLE_PRECISION,&
            &PETSC_DEFAULT_DOUBLE_PRECISION,&
            &PETSC_DEFAULT_DOUBLE_PRECISION,pc_sub_it,globerr)
#endif
       call KSPGetPC(block_ksp_obj,block_pc_obj,globerr)
       call set_sub_pc(block_pc_obj)
    end do
    PERFOFF
  end subroutine set_sub_ksp

  subroutine set_sub_pc(sub_pc_obj)
    PC sub_pc_obj

    select case(sub_pc_type)
    case('ilu')
       call PCSetType(sub_pc_obj,PCILU,globerr)
#if ((PETSC_VERSION_MINOR<2) && (PETSC_VERSION_RELEASE==1))
       if(order_rcm) then
          call PCFactorSetMatOrderingType(sub_pc_obj,MATORDERING_RCM,globerr)
       else
          call PCFactorSetMatOrderingType(sub_pc_obj,MATORDERING_QMD,globerr)
       end if
#else
       if(order_rcm) then
          call PCFactorSetMatOrderingType(sub_pc_obj,MATORDERINGRCM,globerr)
       else
          call PCFactorSetMatOrderingType(sub_pc_obj,MATORDERINGQMD,globerr)
       end if
#endif
       if(mype.eq.0) write(*,'(A,I3)') 'setting pc_factor_level to ',pc_factor_level
       call PCFactorSetLevels(sub_pc_obj,pc_factor_level,globerr)

    case('lu') 
       call PCSetType(sub_pc_obj,PCLU,globerr)
#ifdef PETSC_HAVE_UMFPACK
       if (.not.xy_local) then
          call PCFactorSetMatSolverPackage(sub_pc_obj,MATSOLVERUMFPACK,globerr)
          if(which_ev.eq.'EV') call PetscOptionsInsertString('-mat_umfpack_droptol 1.e-7',globerr)
       end if
#endif
    end select

    call PCSetFromOptions(sub_pc_obj,globerr)
    call PCSetUp(sub_pc_obj,globerr)

  end subroutine set_sub_pc
    
  !>Initializes the L_g_mat matrix, which is a PETSc matrix in MPIAIJ format that contains the 
  !!purely differential part of the linear operator, neglecting the integral terms from the fields
  !!or energy conserving terms in the collision operator
  !!
  !!CalFullRhs has to be initialized before this routine is called
  subroutine initialize_L_g
    complex,dimension(:),allocatable:: g_loc, rhs_loc
    double precision :: time_1, time_2
    integer:: col, memory
    integer:: n_entries
    integer,dimension(2,0:vlen-1):: coords
    integer,dimension(:),allocatable:: d_nnz, o_nnz
    integer:: pair, nbands
    
    PERFON('build_Lg')
    
    nbands=9

#ifdef prec_ara
    if (arakawa_zv) nbands=nbands+4
#endif
    
#ifdef prec_x
    if(.not.xy_local) nbands=nbands+4
#endif
    
#ifdef prec_coll
    if(collision_op.ne.'none') nbands=nbands+6
#endif

    memory=8
    if (prec.eq.'DOUBLE') memory=2*memory
    memory=int(memory*vlen*nbands/1024.**2)

    if (mype.eq.0) then
       write(*,"(A)") '--------------------------------------------------------------'
       write(*,"(A)") ''
       write(*,"(A)") 'Creating explicit representation of the linear operator without the field part'
       write(*,"(a,i8)") 'Rank of the linear operator: ',vlen_gl
       write(*,"(A,i5,A,F8.2,A)") 'Memory requirements: ',memory,' MB per core'
    endif
    
    call get_systime(time_1)
    
    allocate(g_loc(0:vlen-1), rhs_loc(0:vlen-1))

    if(x_local.and.(abs(shat).gt.epsilon(shat)).and.(.not.p_has_0_mode)) periodic_z=.false.
    !maximum number of nonzeros per row in diagonal/off-diagonal portion of local submatrix 
    !for z (-kx) and vpar directions

    call initialize_set_g
    allocate(d_nnz(0:vlen-1),o_nnz(0:vlen-1))
    d_nnz=0
    o_nnz=0
    
    !this loop is used to compute the number of entries in L_g_mat
    do col=0,n_dep-1
       call set_g(col,g_loc,coords,n_entries)
       do pair=0,n_entries-1
          if(coords(1,pair)/vlen.eq.mype) then
             d_nnz(coords(2,pair))=d_nnz(coords(2,pair))+1
          else
             o_nnz(coords(2,pair))=o_nnz(coords(2,pair))+1
          end if
       end do
    enddo

#if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR<3))
    call MatCreateMPIAIJ(MY_MPI_COMM_WORLD,vlen,vlen,vlen_gl,&
         vlen_gl,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,L_g_mat,globerr)
#else
    call MatCreateAIJ(MY_MPI_COMM_WORLD,vlen,vlen,vlen_gl,&
         vlen_gl,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,L_g_mat,globerr)
    call MatMPIAIJSetPreallocation(L_g_mat,PETSC_NULL_INTEGER,d_nnz,&
         PETSC_NULL_INTEGER,o_nnz,globerr)
#endif

    deallocate(d_nnz,o_nnz)

    precond_approx=.true.

    do col=0,n_dep-1
       call set_g(col,g_loc,coords,n_entries)
       call CalFullRhs(g_loc,rhs_loc,0)
       call extract_nonzeros(rhs_loc,coords,n_entries)
    enddo

    precond_approx=.false.
    call finalize_set_g

    call MatAssemblyBegin(L_g_mat, MAT_FINAL_ASSEMBLY, globerr)
    call MatAssemblyEnd(L_g_mat, MAT_FINAL_ASSEMBLY, globerr)

    deallocate(g_loc,rhs_loc)
    
    call get_systime(time_2)
    if (mype.eq.0) then
       write(*,"(A,F8.2)") 'Time to build explicit L_g: ', time_2-time_1  
       write(*,"(A)") ''
       write(*,"(A)") '--------------------------------------------------------------'
       write(*,"(A)") ''
    endif
    PERFOFF
    L_g_initialized=.true.

  end subroutine initialize_L_g

  subroutine finalize_L_g

    call MatDestroy(L_g_mat,globerr)
    L_g_initialized=.false.

  end subroutine finalize_L_g

  subroutine extract_nonzeros(loc_rhs,rescoord,n_entries)
    complex,dimension(0:vlen-1),intent(inout):: loc_rhs !<the distribution function g
    integer,intent(inout)::  n_entries
    integer,dimension(2,0:vlen-1):: rescoord
    integer:: ind, col(1), row(1)
    complex:: val(1)

    do ind=0,n_entries-1
       col=rescoord(1,ind)
       row=rescoord(2,ind)+vlen*mype
       if (impl_shift) then
          !L_g=(1-dt*L)
          val=-dt*loc_rhs(rescoord(2,ind)) 
          if (row(1).eq.col(1))  val=val+(1.0,0.)
       else
          !L_g=L-1*ev_shift
          val=loc_rhs(rescoord(2,ind))
          if (row(1).eq.col(1)) then 
             val=val-pre_shift
          endif
       endif
!       if(mype.eq.0) print*,rescoord(:,ind),val
       call MatSetValues(L_g_mat,1,row,1,col,val,INSERT_VALUES,globerr)
    enddo

  end subroutine extract_nonzeros

  subroutine initialize_set_g
    integer:: zkx_i

    if(.not.xy_local) then
       !coupled in x by finite differences
#ifdef prec_x
       x_np=5
#else
       x_np=1
#endif
       x_bl_ub=(ni0-1)/x_np
    else
       x_np=1
       x_bl_ub=0
    end if

    !coupled in z(-kx)
    if(.not.periodic_z) then
       z_np=5
       allocate(kxzcol(0:4),kxzoff(0:4))
       kxzcol=(nz0*ni0-1)/5
       do zkx_i=0,z_np-1
          if ((kxzcol(zkx_i)*5+zkx_i).ge.nz0*ni0) kxzcol(zkx_i)=kxzcol(zkx_i)-1
          kxzoff(zkx_i)=zkx_i
       end do
    else
       !z boundary points affected by the periodicity have to be computed separately
       z_np=9
       allocate(kxzcol(0:8),kxzoff(0:8))
       !4 boundary points are computed individually
       kxzcol(0:1)=0
       kxzoff(0)=0
       kxzoff(1)=1
       kxzcol(7:8)=0
       kxzoff(7)=nz0-2
       kxzoff(8)=nz0-1

       !(nz0-4) inner points need a distance of 5 points to be independent
       kxzcol(2:6)=(nz0-5)/5
       do zkx_i=2,z_np-3
          if ((kxzcol(zkx_i)*5+zkx_i).ge.(nz0-2)) kxzcol(zkx_i)=kxzcol(zkx_i)-1
          kxzoff(zkx_i)=zkx_i
       end do
    end if

    !coupled in v
    v_np=5
    v_bl_ub=(nv0-1)/v_np

    !coupled in mu only with collisions
    if (collision_op.eq.'none') then
       w_np=1
    else
       w_np=3
    end if
    w_bl_ub=(nw0-1)/w_np

    n_dep=x_np*z_np*v_np*w_np

  end subroutine initialize_set_g

  !>Auxiliary routine that sets all independent entries of g_1 belonging to the specified block
  !!to 1, rest zero, to compute a set of independent columns of the linear operator simultaneously  
  subroutine set_g(dep_ind,loc_g,rescoord,n_entries)
    complex,dimension(0:vlen-1),intent(inout):: loc_g !<the distribution function g
    integer,intent(inout)::  n_entries
    integer,intent(in):: dep_ind
    integer,dimension(2,0:vlen-1):: rescoord

    integer:: x_i,zkx_i,v_i,w_i
    integer:: x_bl,zkx_bl,v_bl,w_bl,xind,zkxind,vind,wind,n
    integer:: ind,res_ind,sten,proc, rest

    !remove old entries
    loc_g=(0.,0.)

    rescoord=0
    n_entries=0

    !all coordinates with the same v_i,zkx_i indices are not connected via 
    !the purely differential operators and can therefore be computed simultaneously 
    w_i=dep_ind/(z_np*x_np*v_np)
    rest=dep_ind - w_i*z_np*x_np*v_np
    v_i=rest/(z_np*x_np)
    rest=rest - v_i*z_np*x_np
    zkx_i=rest/x_np
    x_i=rest - zkx_i*x_np

    do n=ln1,ln2
       do w_bl=0,w_bl_ub
          wind=w_bl*w_np+w_i
          do v_bl=0,v_bl_ub
             vind=v_bl*v_np+v_i
             do zkx_bl=0,kxzcol(zkx_i)
                zkxind=zkx_bl*5+kxzoff(zkx_i)
                do x_bl=0,x_bl_ub
                   xind=x_bl*x_np+x_i
                   call map_to_local(xind,zkxind,vind,wind,n,ind,proc)
                   !skip if not valid
                   if(proc.ge.0) then
                      !initialize g
                      if(mype.eq.proc) then
                         loc_g(ind)=(1.,0.)
                      end if
                      ind=ind+proc*vlen

                      do sten=-2,2
                         call map_to_local(xind,zkxind+sten,vind,wind,n,res_ind,proc)
                         if(mype.eq.proc) then
                            rescoord(:,n_entries)=(/ind, res_ind/)
                            n_entries=n_entries+1
                         end if
                      end do
                      do sten=-(v_np-1)/2,(v_np-1)/2
                         if(sten.ne.0) then
                            call map_to_local(xind,zkxind,vind+sten,wind,n,res_ind,proc)
                            if(mype.eq.proc) then
                               rescoord(:,n_entries)=(/ind, res_ind/)
                               n_entries=n_entries+1
                            end if
                         end if
                      end do
#ifdef prec_ara
                      if (arakawa_zv) then
                         !4 more points in stencil
                         call map_to_local(xind,zkxind+1,vind+1,wind,n,res_ind,proc)
                         if(mype.eq.proc) then
                            rescoord(:,n_entries)=(/ind, res_ind/)
                            n_entries=n_entries+1
                         end if
                         call map_to_local(xind,zkxind+1,vind-1,wind,n,res_ind,proc)
                         if(mype.eq.proc) then
                            rescoord(:,n_entries)=(/ind, res_ind/)
                            n_entries=n_entries+1
                         end if
                         call map_to_local(xind,zkxind-1,vind+1,wind,n,res_ind,proc)
                         if(mype.eq.proc) then
                            rescoord(:,n_entries)=(/ind, res_ind/)
                            n_entries=n_entries+1
                         end if
                         call map_to_local(xind,zkxind-1,vind-1,wind,n,res_ind,proc)
                         if(mype.eq.proc) then
                            rescoord(:,n_entries)=(/ind, res_ind/)
                            n_entries=n_entries+1
                         end if
                      end if
#endif

#ifdef prec_x
                      !compute the positions of the result
                      if(.not.xy_local) then
                         do sten=-2,2
                            if(sten.ne.0) then
                               call map_to_local(xind+sten,zkxind,vind,wind,n,res_ind,proc)
                               if(mype.eq.proc) then
                                  rescoord(:,n_entries)=(/ind, res_ind/)
                                  n_entries=n_entries+1
                               end if
                            end if
                         end do
                      end if
#endif

#ifdef prec_coll
                      if(collision_op.ne.'none') then
                         !collisions have a 3x3 stencil in velocity space. The delta mu=0 row has already been considered, 
                         !two rows are remaining
                         do sten=-1,1
                            call map_to_local(xind,zkxind,vind+sten,wind-1,n,res_ind,proc)
                            if(mype.eq.proc) then
                               rescoord(:,n_entries)=(/ind, res_ind/)
                               n_entries=n_entries+1
                            end if
                         end do
                         do sten=-1,1
                            call map_to_local(xind,zkxind,vind+sten,wind+1,n,res_ind,proc)
                            if(mype.eq.proc) then
                               rescoord(:,n_entries)=(/ind, res_ind/)
                               n_entries=n_entries+1
                            end if
                         end do
                      end if
#endif
                   end if
                end do
             end do
          enddo
       enddo
    enddo

  end subroutine set_g

  subroutine finalize_set_g
    deallocate(kxzoff,kxzcol)
  end subroutine finalize_set_g
  
  
  !>Determines the local processor number and coordinate in a g-like array from 
  !!global zkx,v,w,n indices
  subroutine map_to_local(xind,zkxind,vind,wind,n,loc_ind,loc_proc)
    integer,intent(in):: xind,zkxind,vind,wind,n
    integer,intent(out):: loc_ind,loc_proc

    integer:: zind,kxind,px,pz,pv,pw,lvind,lwind,ps,pind

    !default values
    loc_ind=-1
    loc_proc=-1

    !check validity
    if (wind.ge.nw0) return
    if (wind.lt.0) return
    if (vind.ge.nv0) return
    if (vind.lt.0) return   
    if (.not.periodic_z) then
       if (zkxind.lt.0) return
       if (zkxind.ge.nz0*ni0) return          
    end if

    if (.not.xy_local) then
       if (xind.lt.0) return
       if (xind.ge.ni0) return 
    end if
    
    !map to kx and z indices
    if(xy_local) then
       if(.not.periodic_z) then
          kxind=zkxind/nz0
          zind=zkxind-kxind*nz0
          !change to usual ordering of the kx
          kxind=mod(kxind-(ni0-1)/2+ni0,ni0)
       else
          zind=zkxind
          kxind=0
       endif
       px=0
    else
       zind=zkxind
       px=xind/li0
       kxind=xind-li0*px
    end if
    
    !map to local coordinates (i.e. 0:vlen-1)
    ps=n/ln0
    pind=n-ln0*ps
    
    if (periodic_z) then
       !periodic in z
       zind=mod(zind+nz0,nz0)
    end if
    pz=zind/lk0
    zind=zind-lk0*pz

    pv=vind/ll0
    lvind=vind-ll0*pv

    !lm1 is 0 on each proc
    pw=wind/lm0
    lwind=wind-lm0*pw
    
    loc_ind=kxind+zind*li0+lvind*li0*lk0+lwind*li0*lk0*ll0+pind*li0*lk0*ll0*lm0

    !compute the processor that holds this index
    loc_proc=px+pz*n_procs_x+pv*n_procs_x*n_procs_z+pw*n_procs_x*n_procs_z*n_procs_v+ps*n_procs_x*n_procs_z*n_procs_v*n_procs_w

  end subroutine map_to_local
#endif
end module petsc_precond
