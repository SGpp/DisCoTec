#include "switches.h"
#include "redef.h"
#include "intrinsic_sizes.h"
#undef FLR_NUMERIC
!> The diagnostics module is used for preparation and output
!! of data for post-processing, e.g. moments of the
!! distribution function are calculated and written to output
!! files
Module diagnostics
  use mpi
  Use par_mod
  Use file_io
  Use vel_space, only: mat_00, mat_10, fm, calc_moments, gyro_op_wrapper, calc_vsp_moment
!  Use prefactors
  Use communications
  use aux_fields
  !the gyroaverage will be needed for nonlocal diag_vsp!
  Use gyro_average_ff_mod
  USE gyro_average_df_mod
  USE antenna, ONLY: Apar_pre_antenna,antenna_type, amp, n_antennae, lv_antenna_modes
  use flr_corr
  use aux_func
  use MatrixModule
  use BandedMatrixModule, only: BandedMatrix, get_number_of_bands, &
       &dot_multiply,set_value,add_value
  USE VectorModule, ONLY: Vector, static_size_of_vector
  USE boundaries, only: pb_xshift,ubexc
  use diagnostics_df
  use diagnostics_fsa_moments
  use diagnostics_energy
  use diagnostics_neoclass
  use convergence_monitoring
  use time_averages
  USE diag_finit
  use checkpoint
  use eigen_parameters, only: which_ev, n_ev, ev_prec
  USE geometry, ONLY: c_xy, write_geometry, flux_geomfac
  use Gyro_LES
  use diag_Gyro_LES
#ifdef WITHFUTILS
  use futils 
  use hashtable
#endif
  use axpy
  use x_derivatives, only: x_deriv_exc
  use diagnostics_extended


  Implicit None


  !The basic diag interface
  PUBLIC :: check_diag, initialize_all_diags, exec_all_diags, &
       finalize_all_diags, mem_est_diag, diag_nrg

  !Variables that determine how often routines are called
  PUBLIC :: istep_field, istep_mom, istep_nrg, istep_vsp, istep_g1, &
       istep_neoclass, istep_neoclass2D, istep_gav, & 
       istep_schpt,istep_dfout, write_geom, istep_antenna

  !Variables for running time averages
  PUBLIC :: gav_stime

  !Toggle momentum flux calculation
  PUBLIC :: momentum_flux

  !Routines for eigenvalue analysis (-> move to new module?)
  PUBLIC :: ev_out_open, ev_out_write, ev_out_close

  !Misc (move to new/other module?)
  PUBLIC :: rescale_g_1, set_mats, set_mats_diag_nrg, &
       diag_Blev, diag_trap_levels, cat_output, NRGFILE

  PRIVATE 

  Real :: gav_stime=-1.0
  
  Integer:: istep_field=-1, istep_mom=-1, istep_nrg=-1, istep_vsp=0
  Integer:: istep_dfout=0, istep_antenna=0
   
  Character(Len=8):: filestat='replace', filepos='rewind'

  INTEGER:: NRGFILE, VSPFILE, FIELDFILE, G1FILE, antenna_file
  integer:: nrgcols=10
  
  INTEGER, DIMENSION(:),allocatable:: MOMFILE
  INTEGER, DIMENSION(:),allocatable:: DF_EV_FILE
  character(len=FILENAME_MAX), dimension(:), allocatable :: df_file_out
  integer, dimension(:,:), allocatable :: kx_indices
  integer :: dfinfo !MPI_Info object for tuning IO 
  Integer, allocatable, dimension(:,:) :: dfarr_of_sizes, dfarr_of_subsizes, dfarr_of_starts
  INTEGER(KIND=MPI_OFFSET_KIND), dimension(:), allocatable :: df_offset
  integer, allocatable, dimension(:) :: dfarr
  INTEGER, DIMENSION(:),allocatable:: nkx_keep_array
  integer, dimension(:), allocatable :: df_mpi_handle
  INTEGER:: par_handle
  REAL, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE:: nrg_mat
  INTEGER :: n_nrg_mats = 6, n_moms_Bpar
  INTEGER, DIMENSION(:), allocatable :: nrg_mat_op
  REAL, DIMENSION(:,:), ALLOCATABLE :: pTj_qjB0
  REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE:: mat_20, mat_01, nmat_00, nmat_20, nmat_01
  REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE:: mat_30, mat_11
  COMPLEX, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: gav
  INTEGER:: g1arr
  INTEGER(KIND=MPI_OFFSET_KIND):: g1_offset


  !running average variables
  REAL:: old_gav_time=-1.0

  !administration
  REAL:: last_exec_diag_time

  Integer:: diag_trap_levels=0
  real,dimension(0:9):: diag_Blev= 0.
  character(len=12):: nrgformat
  logical:: write_geom=.true., cat_output=.false., momentum_flux=.true.

  !hdf5 output
#ifdef WITHFUTILS
  ! files, timestep identifiers, and coordinates variables 
  integer, parameter :: bufsize = 20
  integer :: fidnrg_h5, fidfield_h5, fidvsp_h5, fidcoord_h5
  integer, allocatable, dimension(:) :: fidmom_h5
  integer :: isnap_field = 0, isnap_mom = 0, isnap_vsp = 0
  TYPE(BUFFER_TYPE), allocatable, dimension(:) :: hbufnrg_h5
  character (len=6), dimension(1:9) :: mom_label =(/'dens  ', &
       &'T_par ', 'T_perp','q_par ','q_perp','u_par ','densI1','TparI1','TppI1 '/)
  character (len=8), dimension(1:10) :: nrg_label1 =(/'n2      ', 'u_par2  ', &
       &'T_par   ', 'T_per   ', 'Gamma_es', 'Gamma_em', 'Q_es    ','Q_em    ',&
        'P_es    ','P_em    '/)
  character (len=10), dimension(1:10) :: nrg_label2  !see initialize_diag_nrg
  character (len=5), dimension(1:3) :: field_label =(/'phi  ', &
       &'A_par', 'B_par'/)
  character (len=5), dimension(1:5) :: vsp_label =(/'G_es', 'G_em',&
       & 'Q_es', 'Q_em', '<f_>'/)
  character (len=5), allocatable, dimension(:) :: momtrp_label
#endif

CONTAINS

  SUBROUTINE check_diag

    !set derived global variables    
    IF (istep_nrg.LT.0) istep_nrg = 10
    
    IF (istep_field .EQ. -1) THEN
       IF (istep_mom .EQ. -1) THEN
          istep_field = 0
          istep_mom = 0
       ELSE 
          istep_field = istep_mom
       END IF
    ELSE
       IF (istep_mom .EQ. -1) istep_mom = istep_field
    END IF

    IF (istep_field .LT. 1) write_fielddiff = .false.
  
    if(comp_type.eq.'EV') istep_schpt=0

    call check_diag_df
    
    call check_diag_fsa_moments

    if (istep_energy.lt.0) then
       if (xy_local) then
          istep_energy = 100
       else
          istep_energy = 0
       endif
    endif
    if (istep_energy3d.lt.0) istep_energy3d = 0
    if ((istep_energy.gt.0).or.(istep_energy3d.gt.0).or.&
         (istep_energy_exchange.gt.0) ) &
         &call check_diag_energy
    
    if (istep_neoclass .gt. 0) call check_diag_neoclass
    if (istep_neoclass2D .gt. 0) call check_diag_neoclass2D 

    if (istep_antenna.gt.0.and.antenna_type.ne.2.and.antenna_type.ne.3) &
         &istep_antenna=0

    if (.not.momentum_flux) nrgcols=8
    if (istep_nlt.gt.0.and.nonlin_h) stop &
         'diag_nlt diagnostic currently is not implemented for nonlinearity using "h" distribution!'


    if (omega_prec.lt.1.e-15) then
       omega_prec=1.e-15
       if (mype.eq.0) write(*,"(A,E14.6)") "convergence monitoring: setting omega_prec = ",omega_prec
    endif

  END SUBROUTINE check_diag


!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_diag(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc = 0.0

    !module variables
    !mat_20, mat_01
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0
    !mat_30, mat_11
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*pi0*pj0*lklmn0

    if (istep_nrg.gt.0)   mem_loc=mem_est_diag_nrg(mem_loc)
    if (istep_field.gt.0) mem_loc=mem_est_diag_field(mem_loc)
    if (istep_mom.gt.0)   mem_loc=mem_est_diag_mom(mem_loc)
    if (istep_vsp.gt.0)   mem_loc=mem_est_diag_vsp(mem_loc)
    if (istep_gav.gt.0)   mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*lklmn0


    !flr corrections
    mem_loc = mem_est_flr_corrs(mem_loc)

    !nonlocal diags
    if (.not.xy_local) mem_loc=mem_est_all_diags_df(mem_loc)
    
    !neoclassic diag
    if (istep_neoclass.gt.0) mem_loc=mem_est_diag_neoclass(mem_loc)    
    
    if (diag_GyroLES) mem_loc = mem_est_diag_GyroLES(mem_loc)    
    
    mem_est_diag=mem_req_in+mem_loc
  End Function mem_est_diag 

!!!******************************************************************!!!

  !>Initializes GENE internal diagnostics
  !!
  !!Each individual GENE diagnostic initialization should be called 
  !!by this routine
  SUBROUTINE initialize_all_diags
    INTEGER :: k, n

    PERFON('inidiag')
    if (comp_type.ne.'IV') istep_omega=0

    if (istep_nrg.gt.0) call initialize_diag_nrg
    if (istep_field.gt.0) call initialize_diag_field
    if (istep_mom.gt.0) call initialize_diag_mom
    if (istep_vsp.gt.0) call initialize_diag_vsp
    if (istep_omega.gt.0) call initialize_diag_convergence
!    if (istep_g1.gt.0) call initialize_diag_g1
    if (istep_gav.gt.0) call initialize_diag_gav
    if (write_checkpoint.or.(istep_schpt.gt.0).or.(istep_g1.gt.0).or. &
         &(reset_limit.ne.1000).or.(istep_gav.gt.0)) then
       call initialize_checkpoint_write
    end if
    call initialize_all_diags_energy
    call initialize_all_diags_nlt

    if (istep_antenna.gt.0) call initialize_diag_antenna
    
    if (diag_GyroLES) call initialize_all_diags_GyroLES
    if (GyroLES) call initialize_GyroLES
    
    IF (write_fielddiff) CALL initialize_all_diags_finit(istep_field)

    if(istep_dfout.gt.0) call initialize_diag_df_ev

    if (.not.xy_local) call initialize_all_diags_df
    
    if (.not.xy_local) call initialize_all_diags_fsa_moments

    if (istep_neoclass.gt.0) call initialize_all_diags_neoclass 
    if (istep_neoclass2D.gt.0) call initialize_all_diags_neoclass2D 

    
    !Some arrays which are shared by several diagnostics
    !In a next step, they should be encapsulated in a new structure,
    !routine, or whatever
    ALLOCATE(mat_20(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2))
    ALLOCATE(mat_30(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    ALLOCATE(mat_01(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2))
    ALLOCATE(mat_11(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    IF (n_fields .GT. 2) THEN
      ALLOCATE(nmat_00(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))
      ALLOCATE(nmat_20(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))
      ALLOCATE(nmat_01(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))

      ALLOCATE(pTj_qjB0(lk1:lk2,ln1:ln2))
      DO n = ln1, ln2
        DO k = lk1, lk2
          pTj_qjB0(k,n) = spec(n)%temp / (spec(n)%charge * geom%Bfield(pi1,pj1,k))
        END DO
      END DO
    ENDIF

    call set_mats

    call initialize_flr_corrs

#ifdef WITHFUTILS
    if ((write_h5).and.(mype==0)) then
       !coordinates
       call creatf(trim(diagdir)//'/coord'//trim(file_extension)//'.h5',fidcoord_h5, "coordinates", 'd')
       call creatg(fidcoord_h5, '/coord')
       
       if(xy_local) then
          call putarr(fidcoord_h5, "/coord/kx", kx   , "k_x")
       else
          call putarr(fidcoord_h5, "/coord/x" , xval , "x"  )
       end if
       call putarr(fidcoord_h5, "/coord/ky", ky   , "k_y")
       call putarr(fidcoord_h5, "/coord/z" , zval , "z"  )
       call putarr(fidcoord_h5, "/coord/vp", vp   , "vp"  )
       call putarr(fidcoord_h5, "/coord/mu", mu   , "mu" )
       call putarr(fidcoord_h5, "/coord/vp_weight", vp_weight, "vp_weight" )
       call putarr(fidcoord_h5, "/coord/mu_weight", mu_weight, "mu_weight" )
       call closef(fidcoord_h5)
    end if ! write_h5 
#endif

    !write geometric quantities which are partially needed for idl diagnostics
    if(write_geom.and.(.not.cat_output)) call write_geometry

    !reset for new problem
    last_exec_diag_time = -3.14159265

    PERFOFF

  END SUBROUTINE initialize_all_diags

!!!******************************************************************!!!
  SUBROUTINE set_mats
    INTEGER :: k, l, m, n

    DO m=lm1,lm2
       DO l=ll1,ll2    
          DO k=lk1,lk2
             ! temperatures
             mat_20(:,:,k,l,m) = vp(l)*vp(l) * mat_00(:,:,k,l,m)
             mat_01(:,:,k,l,m) = mu(m)*geom%Bfield(pi1:pi2,:,k) * mat_00(:,:,k,l,m)
             DO n=ln1,ln2
                mat_30(:,:,k,l,m,n) = vp(l)*vp(l) * mat_10(:,:,k,l,m,n) !q_||,||
                mat_11(:,:,k,l,m,n) = mu(m)*geom%Bfield(pi1:pi2,:,k) *&
                     mat_10(:,:,k,l,m,n) !q_||,perp
             END DO

             IF (n_fields .GT. 2) THEN
               nmat_00(:,:,k,l,m) = & ! densI1
                 mu(m) * geom%Bfield(pi1:pi2,:,k) * mat_00(:,:,k,l,m)
               nmat_20(:,:,k,l,m) = & ! TparI1
                 mu(m) * geom%Bfield(pi1:pi2,:,k) * mat_20(:,:,k,l,m)
               nmat_01(:,:,k,l,m) = & ! TperpI1
                 mu(m) * geom%Bfield(pi1:pi2,:,k) * mat_01(:,:,k,l,m)
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE set_mats

!!!******************************************************************!!!
  SUBROUTINE exec_all_diags(itime,time,overflow_exit,underflow_exit,reset)
    integer, intent(in) :: itime
    real, intent(in) :: time
    logical, intent(inout) :: overflow_exit, underflow_exit, reset
    logical:: no_aux_yet
    !For nrg call
    real, allocatable, dimension(:,:) :: nrg_data

    !avoid double entries if exec_all_diags is, e.g., called during
    !convergence exit
    IF (time.eq.last_exec_diag_time) THEN
       RETURN
    ELSE
       last_exec_diag_time=time
       !keep track whether calc_aux_fields has already been called
       !(should thus be provided to 'call_diag' if the sub-diagnostic
       ! relies on either emfields or f):
       no_aux_yet=.true.
    ENDIF

    PERFON('diag_exec')
    
    IF (call_diag(itime,istep_neoclass,no_aux_yet)) THEN
       CALL exec_all_diags_neoclass
    END IF

    IF (call_diag(itime,istep_neoclass2D,no_aux_yet)) THEN 
       CALL exec_all_diags_neoclass2D
    END IF


    allocate(nrg_data(nrgcols,ln1:ln2))
    IF (call_diag(itime,istep_nrg,no_aux_yet)) &
       &CALL diag_nrg(overflow_exit,underflow_exit,nrg_data,.true.)
    deallocate(nrg_data)

    IF ((istep_omega.GT.0).AND.(.not.convergence_exit)) THEN
       IF (MODULO(itime,istep_omega) == 0) THEN
          if (no_aux_yet) then
             call calc_aux_fields(g_1,emfields,f_,.false.)
             no_aux_yet=.false.
          endif
          IF (itime.GT.0) then
             CALL diag_convergence
          ELSE
             call set_phi_old
          ENDIF
       END IF
    END IF

    IF (call_diag(itime,istep_field,no_aux_yet)) &
         &CALL diag_field

    IF (call_diag(itime,istep_mom,no_aux_yet)) &
         &CALL diag_mom

    IF (call_diag(itime,istep_dfout,no_aux_yet)) THEN
       if(dfout_mpiio) then
          CALL diag_df_ev_mpi
       else
          CALL diag_df_ev
       end if
    END IF

    IF (call_diag(itime,istep_energy,no_aux_yet).OR.&
         call_diag(itime,istep_energy3d,no_aux_yet)) &
         &call exec_all_diags_energy(itime, istep_energy, &
         &istep_energy3d)

    IF (istep_energy_exchange.gt.0) THEN
       IF(MODULO(itime,istep_energy_exchange).eq.0) THEN   
         call diag_energy_exchange
       END IF
    END IF

    IF (call_diag(itime,istep_nlt,no_aux_yet)) &
         &call exec_all_diags_nlt(itime, istep_nlt)

    IF (call_diag(itime,istep_vsp,no_aux_yet)) &
         & CALL diag_vsp

    IF (write_fielddiff) CALL exec_all_diags_finit(itime,time)

    if (.not.xy_local) call exec_all_diags_df(itime,time,reset)
    
    if (.not.xy_local) call exec_all_diags_fsa_moments(itime,time)

    IF (call_diag(itime,istep_g1)) call checkpoint_write(g_1,'g1')

    IF (call_diag(itime,istep_gav)) CALL diag_gav

    if (call_diag(itime,istep_antenna)) call diag_antenna

    ! write secure checkpoint
    IF (call_diag(itime,istep_schpt)) &
         &call checkpoint_write(g_1,'s_checkpoint')

    !for eigenvalue runs
    if((comp_type.eq.'EV').and.write_checkpoint) then
         call checkpoint_write(g_1,'checkpoint')
    end if

    !for neoclassics runs
    if((comp_type.eq.'NC').and.write_checkpoint) &
         &call checkpoint_write(g_1,'checkpoint')
    
    !For gyroLES routine 
    IF (GyroLES) THEN
       IF (MODULO(itime,istep_GyroLES) .EQ. 0) THEN
          IF (no_aux_yet) THEN
             if (.not.nonlin_h) then
                call calc_aux_fields(g_1,emfields,f_,.false.)  !only emfields required
             else
                call calc_aux_fields(g_1,emfields,f_,.true.,h_out=h_)
             endif
             no_aux_yet=.false.
          ENDIF
          CALL exec_GyroLES
       ENDIF
    ENDIF
    !For gyroLES diagnostics
    !Warning, after this routine, f_ and h_ might be interchanged
    !due to calling calc_rhs_only
    IF (diag_GyroLES) call exec_all_diags_GyroLES(itime,time)    

   
    PERFOFF
    
  END SUBROUTINE exec_all_diags


!!!******************************************************************!!!
  !>Finalizes all GENE internal diagnostics
  SUBROUTINE finalize_all_diags

#ifdef COMBI_MGR
    !mh pass checkpoint data to interface
    if ((comp_type.eq.'IV').and.write_checkpoint) then
      read_checkpoint = .true.
      !call checkpoint_write(g_1,'checkpoint')
      call checkpoint_write_memory(g_1,time,dt,li1, li2, lj1, lj2, lk1, lk2, ll1, &
                                   ll2, lm1, lm2, ln1, ln2,ni0,nj0,nz0,nv0,nw0, &
                                   n_spec, comm_cart )
    end if
#else
    if ((comp_type.eq.'IV').and.write_checkpoint) call checkpoint_write(g_1,'checkpoint')
#endif

    IF (istep_omega.GT.0) call finalize_diag_convergence
    IF (istep_nrg.GT.0) call finalize_diag_nrg
    IF (istep_field.GT.0) call finalize_diag_field
    IF (istep_mom.GT.0) call finalize_diag_mom
    IF (istep_vsp.GT.0) call finalize_diag_vsp
!    if (istep_g1.GT.0) call finalize_diag_g1
    if (istep_gav.GT.0) call finalize_diag_gav
    if (write_checkpoint.or.(istep_schpt.gt.0).or.(istep_g1.gt.0).or. &
         &(reset_limit.ne.1000).or.(istep_gav.gt.0)) call finalize_checkpoint_write
    if (istep_antenna.gt.0) call finalize_diag_antenna
    
    if (.not.xy_local) call finalize_all_diags_df
    
    if (.not.xy_local) call finalize_all_diags_fsa_moments
    
    if (istep_neoclass.GT.0) call finalize_all_diags_neoclass
    if (istep_neoclass2D.GT.0) call finalize_all_diags_neoclass2D

    call finalize_all_diags_energy
    call finalize_all_diags_nlt

    IF (diag_GyroLES) call finalize_all_diags_GyroLES
    IF (GyroLES) call finalize_GyroLES
    
    IF (write_fielddiff) CALL finalize_all_diags_finit

    if(istep_dfout.gt.0) call finalize_diag_df_ev

    DEALLOCATE(mat_20,mat_30,mat_01,mat_11)
    IF (n_fields .GT. 2) DEALLOCATE(nmat_00,nmat_20,nmat_01,pTj_qjB0)

    call finalize_flr_corrs

  END SUBROUTINE finalize_all_diags

!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Give an estimate of the memory requirements of diag_nrg
  Real Function mem_est_diag_nrg(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !nrg_mats
    mem_loc = 6.*SIZE_OF_COMPLEX_MB*pi0*pj0*lklmn0
    IF (n_fields .GT. 2) mem_loc = mem_loc + SIZE_OF_COMPLEX_MB * pi0 * pj0 * lklmn0
    !local variables
    !momc
    mem_loc = mem_loc + lijk0*6*ln0*SIZE_OF_COMPLEX_MB
    IF (n_fields .GT. 2) mem_loc = mem_loc + lijk0 * ln0 * SIZE_OF_COMPLEX_MB
    !mom
    mem_loc = mem_loc + lijk0*8*ln0*SIZE_OF_REAL_MB
    !dfields_dy
    mem_loc = mem_loc + lijk0*n_fields*SIZE_OF_COMPLEX_MB
    !corrections (target?)
    mem_loc = mem_loc + lijk0*ln0*3*SIZE_OF_COMPLEX_MB
    IF (n_fields .GT. 2) mem_loc = mem_loc + lijk0 * ln0 * 2 * SIZE_OF_COMPLEX_MB
    IF (equil_par_curr) mem_loc = mem_loc + lijk0 * ln0 * 3 * SIZE_OF_COMPLEX_MB

    mem_est_diag_nrg = mem_loc+mem_req_in

  End Function mem_est_diag_nrg

!!!******************************************************************!!!

  SUBROUTINE initialize_diag_nrg
    integer:: ndigits   
#ifdef WITHFUTILS
    integer :: n

    !BlueGene compiler has problems with array definition
    nrg_label2(1)='n^2'
    nrg_label2(2)='u_\\par^2'
    nrg_label2(3)='T_\\par'
    nrg_label2(4)='T_\\perp'
    nrg_label2(5)='Gamma_{es}'
    nrg_label2(6)='Gamma_{em}'
    nrg_label2(7)='Q_{es}'
    nrg_label2(8)='Q_{em}'
    nrg_label2(9)='P_{es}'
    nrg_label2(10)='P_{em}'
#endif

    if (mype==0) then
       if (.not.cat_output) then
          call get_unit_nr(NRGFILE)
          OPEN(NRGFILE, file=trim(diagdir)//'/nrg'&
               &//trim(file_extension), form='formatted',&
               status=filestat, position=filepos)
       endif
#ifdef WITHFUTILS
       if(write_h5) then
          call creatf(trim(diagdir)//'/nrg'//trim(file_extension)//'.h5', &
               fidnrg_h5, "phase space averages", 'd')
          call creatg(fidnrg_h5, '/nrg')
          allocate(hbufnrg_h5(0:n_spec-1))
          do n=0, n_spec-1
             call creatg(fidnrg_h5, '/nrg'//trim(spec(n)%name))
             call htable_init(hbufnrg_h5(n), bufsize)
             call set_htable_fileid(hbufnrg_h5(n), fidnrg_h5, '/nrg'//trim(spec(n)%name))
          end do
       end if
#endif
    endif

    IF (n_fields .GT. 2) n_nrg_mats = 8

    ALLOCATE(nrg_mat(pi1:pi2,pj1:pj2,lk1:lk2,1:n_nrg_mats,ll1:ll2,lm1:lm2,ln1:ln2))
    ALLOCATE(nrg_mat_op(1:n_nrg_mats))

    call set_mats_diag_nrg

    !set the number of digits for output
    ndigits=int(log10(1.0/omega_prec))
    if (ndigits<4) ndigits=4
    write(nrgformat,'(A,I2.2,A,I2.2,A,I2.2,A)') '(',nrgcols,'ES',8+ndigits,'.',ndigits,')'
    
    call initialize_avgfluxes

  END SUBROUTINE initialize_diag_nrg


!!!******************************************************************!!!
  SUBROUTINE set_mats_diag_nrg

    Integer :: j, k,l,m,n
    DO n=ln1,ln2
       DO m=lm1,lm2
          DO l=ll1,ll2    
             DO k=lk1,lk2
                DO j=pj1,pj2
                   nrg_mat(:,j,k,1,l,m,n) = mat_00(:,j,k,l,m)
                   nrg_mat(:,j,k,2,l,m,n) = mat_10(:,j,k,l,m,n)
                   nrg_mat(:,j,k,3,l,m,n) = vp(l)*vp(l) * mat_00(:,j,k,l,m) !mat_20
                   nrg_mat(:,j,k,4,l,m,n) = mu(m)*geom%Bfield(pi1:pi2,j,k) * mat_00(:,j,k,l,m) !mat_01
                   nrg_mat(:,j,k,5,l,m,n) = nrg_mat(:,j,k,3,l,m,n)+nrg_mat(:,j,k,4,l,m,n) !mat_20+mat_01
                   nrg_mat(:,j,k,6,l,m,n) = (vp(l)*vp(l) + mu(m)*geom%Bfield(pi1:pi2,j,k))*&
                        mat_10(:,j,k,l,m,n) !mat_30 + mat_11
                   IF (n_fields .GT. 2) THEN ! (nrg_mats to be computed with I1)
                      nrg_mat(:,j,k,7,l,m,n) = mu(m) * geom%Bfield(pi1:pi2,j,k) * mat_00(:,j,k,l,m)
                      nrg_mat(:,j,k,8,l,m,n) = mu(m) * geom%Bfield(pi1:pi2,j,k) * nrg_mat(:,j,k,5,l,m,n)
                   ENDIF
                END DO
             END DO
          END DO
       END DO
    END DO

    nrg_mat_op = 0
    IF (n_fields.GT.2) nrg_mat_op(7:8) = 1

  END SUBROUTINE set_mats_diag_nrg

!!!******************************************************************!!!
  !> NRG diagnostic: Calculates volume averages of velocity space moments:
  !! 1.) <|n|^2> / (n_j0 rho_ref/L_ref)^2
  !! 2.) <|u_par|^2> / (v_Tj rho_ref/Lref)^2
  !! 3.) <|T_par|^2> / (T_j0 rho_ref/Lref)^2
  !! 4.) <|T_perp|^2> / (T_j0 rho_ref/Lref)^2
  !! 5.) G_es / (c_ref n_ref (rho_ref/Lref)^2)
  !! 6.) G_em / (c_ref n_ref (rho_ref/Lref)^2)
  !! 7.) Q_es / (c_ref p_ref (rho_ref/Lref)^2) 
  !! 8.) Q_em / (c_ref p_ref (rho_ref/Lref)^2)
  !! 9.) P_es / (c_ref^2 m_ref n_ref (rho_ref/Lref)^2)
  !! 10.) P_em / (c_ref^2 m_ref n_ref (rho_ref/Lref)^2)
  Subroutine diag_nrg(overflow_exit,underflow_exit,var,write_nrg)
    Logical, intent(out) :: overflow_exit, underflow_exit
    logical, intent(in) :: write_nrg
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2,n_nrg_mats,ln1:ln2):: momc
    Real,    Dimension(li1:li2,lj1:lj2,lk1:lk2,nrgcols,ln1:ln2):: mom
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields) :: dfields_dy
    Real, Dimension(nrgcols,ln1:ln2), intent(out) :: var
    Real, Dimension(nrgcols,ln1:ln2) :: rat
    Integer:: i, k, n, o, v, ierr, pes
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2,1:n_flr_corrs), TARGET :: corrections
    logical :: OUTPUT=.false.
    real ::local_sum, global_sum
#ifdef WITHFUTILS
    integer:: ng, isnap_nrg
#endif

    !PERFON('d_nrg')

    !PERFON('d_nrg1')    
    momc=cmplx(0,0)

    ! f_ has boundary points in all dimensions
#ifndef FLR_NUMERIC
    IF (OUTPUT) THEN
       CALL calculate_test_sum(f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "f_ = ",global_sum
    END IF

    call calc_moments(n_nrg_mats,.true.,f_,nrg_mat,momc,p_gy_op=nrg_mat_op)
    IF (OUTPUT) THEN
       CALL calculate_test_sum(momc(:,:,:,1,ln1),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "momc = ",global_sum
    END IF
#else
    if (n_fields.gt.2) stop 'calc_vsp_moment not implemented yet for bpar'
    call calc_vsp_moment(n_nrg_mats,f_,.true.,emfields,nrg_mat,momc,.false.)
#endif
    
    !>00 mode is NOT zeroed out for nrg data
    !if (p_has_00_mode) momc(li1,lj1,:,:,:)=0.0

    !PERFOFF
    do o=1,n_fields
       do k=lk1,lk2
          call calc_dfields_dy(dfields_dy(:,:,k,o),k,o)
       enddo
    enddo

    do n=ln1,ln2
       ! particle flux (5 = es, 6 = em)
       mom(:,:,:,5,n) = -REAL(CONJG(momc(:,:,:,1,n)) *  dfields_dy(:,:,:,1))
       ! heat flux (7 = es, 8 = em)
       mom(:,:,:,7,n) = -Real(Conjg(momc(:,:,:,5,n)) *  dfields_dy(:,:,:,1))
       ! momentum flux (9 =es, 10 = em)
       if (momentum_flux) &
            mom(:,:,:,9,n) = -real(conjg(momc(:,:,:,2,n)) *  dfields_dy(:,:,:,1))
       !same for electromagnetic contribution
       if (n_fields.gt.1) then
          mom(:,:,:,6,n) = REAL(CONJG(momc(:,:,:,2,n)) *  dfields_dy(:,:,:,2))
          mom(:,:,:,8,n) = Real(Conjg(momc(:,:,:,6,n)) *  dfields_dy(:,:,:,2))
          if (momentum_flux) &
               mom(:,:,:,10,n) = real(conjg(momc(:,:,:,3,n)) * dfields_dy(:,:,:,2))
          IF (n_fields .GT. 2) THEN
            DO k = lk1, lk2
              mom(:,:,k,6,n) = mom(:,:,k,6,n) - pTj_qjB0(k,n) * REAL(CONJG(momc(:,:,k,7,n))*dfields_dy(:,:,k,3))
              mom(:,:,k,8,n) = mom(:,:,k,8,n) - pTj_qjB0(k,n) * REAL(CONJG(momc(:,:,k,8,n))*dfields_dy(:,:,k,3))
            END DO
          END IF
       else
          mom(:,:,:,6,n) = 0.
          mom(:,:,:,8,n) = 0.
          if (momentum_flux) mom(:,:,:,10,n) = 0.
       endif
    enddo

    !-------------- Sum over velocity space  -----------------
    Call my_sum_to_0_complex(momc, Size(momc), mpi_comm_vw)

    do n=ln1,ln2
       do v=5,nrgcols
          Call my_sum_to_0_real(mom(:,:,:,v,n), Size(mom(:,:,:,v,n)), mpi_comm_vw)
       enddo
    enddo

#ifndef FLR_NUMERIC
    ! now calculate some additional FLR correction terms
    ! corrections(:,:,:,:,1) : correction for dens
    ! corrections(:,:,:,:,2) : correction for tpar
    ! corrections(:,:,:,:,3) : correction for tperp
    ! corrections(:,:,:,:,4:6) : only for Bpar_off = .f. case
    ! corrections(:,:,:,:,n_flr_corrs-2:n_flr_corrs) : only for j_0parallel case
    IF (n_fields .LE. 2) THEN
      CALL correct_for_FLR(emfields(:,:,:,1), corrections)
    ELSE
      CALL correct_for_FLR_bpar(emfields(:,:,:,1),emfields(:,:,:,3),corrections)
    END IF
#else
    corrections = 0.0
#endif
    if (x_local) then !used to preserve omp parallel workspace for x_local versions
       do n=ln1,ln2
          ! Add polarization density
          momc(:,:,:,1,n) = momc(:,:,:,1,n) +corrections(:,:,:,n,1)

          IF (equil_par_curr) momc(:,:,:,2,n) = momc(:,:,:,2,n) + &
            corrections(:,:,:,n,n_flr_corrs-2)
          
          ! particle flux should also be corrected
          mom(:,:,:,5,n) = mom(:,:,:,5,n)-REAL( CONJG(corrections(:,:,:,n,1))*dfields_dy(:,:,:,1))

          if (n_fields.gt.1.and.momentum_flux) then
             !correction for electromagnetic momentum flux
             mom(:,:,:,10,n) = mom(:,:,:,10,n)-real(conjg(corrections(:,:,:,n,2))*dfields_dy(:,:,:,2))
          endif
          IF ((equil_par_curr) .AND. (n_fields .GT. 1)) mom(:,:,:,6,n) = mom(:,:,:,6,n) + &
            SQRT(2.0*spec(n)%temp/spec(n)%mass) * &
            REAL(CONJG(corrections(:,:,:,n,n_flr_corrs-2))*dfields_dy(:,:,:,2))
          IF (n_fields .GT. 2) THEN
            DO k = lk1, lk2
              mom(:,:,k,6,n) = mom(:,:,k,6,n) - pTj_qjB0(k,n) * &
                REAL(CONJG(corrections(:,:,k,n,4))*dfields_dy(:,:,k,3))
            END DO
          END IF
          
          ! Add FLR Correction terms to the temperatures and heat fluxes
          momc(:,:,:,3,n) = 2*(momc(:,:,:,3,n)+corrections(:,:,:,n,2)) &
               &-momc(:,:,:,1,n)
          momc(:,:,:,4,n) = momc(:,:,:,4,n)+corrections(:,:,:,n,3)-momc(:,:,:,1,n)
          mom(:,:,:,7,n)  = mom(:,:,:,7,n)-&
               &REAL( CONJG((corrections(:,:,:,n,2)+corrections(:,:,:,n,3)))*dfields_dy(:,:,:,1))
          IF ((equil_par_curr) .AND. (n_fields .GT. 1)) mom(:,:,:,8,n) = mom(:,:,:,8,n) + &
            SQRT(2.0*spec(n)%temp/spec(n)%mass) * &
            REAL(CONJG(corrections(:,:,:,n,n_flr_corrs-1)+corrections(:,:,:,n,n_flr_corrs))*dfields_dy(:,:,:,2))
          IF (n_fields .GT. 2) THEN
            DO k = lk1, lk2
              mom(:,:,k,8,n) = mom(:,:,k,8,n) - pTj_qjB0(k,n) * &
                REAL(CONJG(corrections(:,:,k,n,5)+corrections(:,:,k,n,6))*dfields_dy(:,:,k,3))
            END DO
          END IF
       enddo
    else  !with temp and dens profile information
       do n=ln1,ln2
          ! Add polarization density
          momc(:,:,:,1,n) = momc(:,:,:,1,n)+corrections(:,:,:,n,1)
          if (n_fields.gt.1.and.momentum_flux) then
             !correction for electromagnetic momentum flux
             mom(:,:,:,10,n) = mom(:,:,:,10,n)-real(conjg(corrections(:,:,:,n,2))*dfields_dy(:,:,:,2))
          endif

          ! Add profile information to u_par and FLR Correction terms 
          ! to temperatures 
          do i=li1,li2
             momc(i,:,:,2,n) = momc(i,:,:,2,n)/spec(n)%dens_prof(i)
             momc(i,:,:,3,n) = (2.0D0*(momc(i,:,:,3,n)+corrections(i,:,:,n,2)) &
                  &-spec(n)%temp_prof(i)*momc(i,:,:,1,n))/spec(n)%dens_prof(i)
             momc(i,:,:,4,n) = (momc(i,:,:,4,n)+corrections(i,:,:,n,3)-&
                  &spec(n)%temp_prof(i)*momc(i,:,:,1,n))/spec(n)%dens_prof(i)
          enddo
          ! electrostatic particle and heat fluxes should also be corrected
          mom(:,:,:,5,n) = mom(:,:,:,5,n)-REAL( CONJG(corrections(:,:,:,n,1))*dfields_dy(:,:,:,1))
          mom(:,:,:,7,n) = mom(:,:,:,7,n)-&
               &REAL( CONJG((corrections(:,:,:,n,2)+corrections(:,:,:,n,3)))*dfields_dy(:,:,:,1))
       enddo
    endif

    mom(:,:,:,1:4,:) = Real(momc(:,:,:,1:4,:)*Conjg(momc(:,:,:,1:4,:)))
   
    ! take into account negative kj, except for 0-mode
    mom = 2. * mom
    if (p_has_0_mode)  then
       if (xy_local.and.yx_order) then
          mom(li1,:,:,:,:) = 0.5 * mom(li1,:,:,:,:)
       else
          mom(:,lj1,:,:,:) = 0.5 * mom(:,lj1,:,:,:)
       end if
    end if

    !space averaging
    call sum_3d(nrgcols,mom,var)
    
    !change normalization (to species independent)
    do n=ln1,ln2
       var(2,n) = var(2,n) * spec(n)%mass/spec(n)%temp
       var(5,n) = var(5,n) * spec(n)%dens
       var(6,n) = var(6,n) * spec(n)%dens
       var(7,n) = var(7,n) * spec(n)%dens*spec(n)%temp
       var(8,n) = var(8,n) * spec(n)%dens*spec(n)%temp
       if (momentum_flux) then
          var(9,n) = var(9,n) * spec(n)%dens*spec(n)%mass
          var(10,n) = var(10,n) * spec(n)%dens*2*spec(n)%temp
       endif
    enddo

    ! data of all species is written in one file       
    If (mype == 0.and.write_nrg) Write(NRGFILE, "(F13.6)") time
    Do pes = 0, n_procs_s-1
       ! send to process 0 
       Call my_real_get(rat, var, Size(var),&
            0, pexyzvwspec(0, 0, 0, 0, 0, pes))
       ! check for overflows
       IF (mype==0) THEN
          IF ((MAXVAL(rat).GT.overflow_limit) &
               .or. (rat(1,ln1).ne.rat(1,ln1))) THEN
             overflow_exit=.TRUE.
          ENDIF
          IF ((.not.nonlinear).and. &
               & (MAXVAL(rat).LT.underflow_limit)) &
               underflow_exit=.TRUE.          
       END IF
       ! write data for all species on pes to file
       If (mype == 0.and.write_nrg) then
          do n=ln1,ln2
             Write(NRGFILE, nrgformat) rat(:,n)
          enddo
#ifdef COMBI_MGR
          call set_nrg( time, rat(1,ln1) )
#endif
       end If
       
#ifdef WITHFUTILS
      if(write_nrg) then
       if ((write_h5).and.(mype).eq.0) then
          isnap_nrg = itime/istep_nrg + 1
          call attach(fidnrg_h5, "/nrg", "n_steps", isnap_nrg)    
          do n=ln1,ln2
             ng = pes*ln0 + n-ln1
             call add_record(hbufnrg_h5(ng), "time", "simulation time", real(time,8))
             do o=1,8
                call add_record(hbufnrg_h5(ng), trim(nrg_label1(o)), &
                     &trim(nrg_label2(o)),real(rat(o,n),8))
             enddo
             call htable_endstep(hbufnrg_h5(ng))
             call htable_hdf5_flush(hbufnrg_h5(ng))
          end do
          call flushh5(fidnrg_h5)
          
       end if
      end if !write_nrg
#endif
    End Do

    !calculate running average
    call set_avgfluxes(var,emfields,time,nrgcols)

    ! if any of the rat arrays of the different species triggered an overflow
    ! send it to all other processes.
    CALL MPI_bcast(overflow_exit, 1, MPI_LOGICAL, 0, MY_MPI_COMM_WORLD, ierr)
    CALL MPI_bcast(underflow_exit, 1, MPI_LOGICAL, 0, MY_MPI_COMM_WORLD, ierr)
    If ((mype == 0).and.(.not.cat_output)) call flush(NRGFILE)
    !PERFOFF

  End Subroutine diag_nrg


  Subroutine finalize_diag_nrg

      if (mype==0) then
         if (.not.cat_output) CLOSE(NRGFILE)
#ifdef WITHFUTILS
         if(write_h5) then
            call closef(fidnrg_h5)
            deallocate(hbufnrg_h5)
         end if
#endif
      endif

      Deallocate(nrg_mat,nrg_mat_op)

      call finalize_avgfluxes(time, istep_nrg)
      
  End Subroutine finalize_diag_nrg

!!!******************************************************************!!!
!!!******************************************************************!!!
  
  !>Give an estimate of the memory requirements of diag_field
  Real Function mem_est_diag_field(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc = lijk0*SIZE_OF_COMPLEX_MB
    mem_est_diag_field = mem_loc+mem_req_in

  End Function mem_est_diag_field

!!!******************************************************************!!!

  Subroutine initialize_diag_field
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: rank, o
#endif

    if(write_std) then
       call get_unit_nr(FIELDFILE)    

       if(mype.eq.0) open(FIELDFILE, file=trim(diagdir)//'/field'//&
            &trim(file_extension),&
            form='unformatted',status=filestat,position=filepos)
    end if

#ifdef WITHFUTILS
    if (write_h5.and.((my_pev+my_pew+my_pespec).eq.0)) then
       isnap_field = 0
       call creatf(trim(diagdir)//'/field'//trim(file_extension)//'.h5', &
            fidfield_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidfield_h5, '/field')
       rank = 0
       call creatd(fidfield_h5, rank, dims, "/field/time", "time")
       
       do o=1,n_fields
          call creatg(fidfield_h5, '/field/'//trim(field_label(o)))
       enddo
    end if
#endif
    
  End Subroutine initialize_diag_field


  !*****************************************************************
  !>FIELD diagnostic
  !!
  !!Writes 3D data (kx,ky,z) of electrostatic field
  !!phi and electromagnetic field Apar (if beta > 0)
  !!and Bpar (if beta > 0 and BPAR precompiler switch active)
  !!to field.dat file
  SUBROUTINE diag_field
    integer:: o
    Complex, Dimension(li1:li2, lj1:lj2,lk1:lk2) :: tmp_field
#ifdef WITHFUTILS
    character(len=FILENAME_MAX) :: dset_name_field
#endif

    !PERFON('d_field')

    if ((my_pev+my_pew+my_pespec).eq.0) then
       if(write_std) then
          IF (mype.eq.0) WRITE(FIELDFILE) time
          do o=1,n_fields
             !to avoid temporary array in call to write_3d:
             tmp_field = emfields(li1:li2,lj1:lj2,lk1:lk2,o)
             IF ((o .EQ. 2) .AND. (antenna_type.eq.1)) tmp_field = Apar_pre_antenna(:,:,lk1:lk2)
             call write3d(FIELDFILE,tmp_field,0)
          enddo
       end if

#ifdef WITHFUTILS
      if(write_h5) then
         call append(fidfield_h5, "/field/time", real(time,8))
         call attach(fidfield_h5, "/field/time", "n_steps", isnap_field+1)
         
         do o=1,n_fields
            write(dset_name_field, "(A, '/', i10.10)") "/field/"//trim(field_label(o)), isnap_field
            !to avoid temporary array in call to putarrnd:
            tmp_field = emfields(li1:li2,lj1:lj2,lk1:lk2,o)
            IF ((o .EQ. 2) .AND. (antenna_type.eq.1)) tmp_field = Apar_pre_antenna(:,:,lk1:lk2)
            call putarrnd(fidfield_h5, dset_name_field, tmp_field, (/3, 2, 1/))
            call attach(fidfield_h5, dset_name_field, "time", time)      
         enddo

         isnap_field = isnap_field+1
         call flushh5(fidfield_h5)
      end if
#endif

      !write
      if(write_std) then
         IF (mype.eq.0) call flush(FIELDFILE)
      end if
   endif

    !PERFOFF

  END SUBROUTINE diag_field

  Subroutine finalize_diag_field

    if(write_std) then
       if (mype.eq.0) close(FIELDFILE)
    end if

#ifdef WITHFUTILS
    if(write_h5) then
       if ((my_pev+my_pew+my_pespec).eq.0) then

          call closef(fidfield_h5)
       end if
    end if
#endif

  End Subroutine finalize_diag_field


!!!******************************************************************!!!
!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of diag_mom
  Real Function mem_est_diag_mom(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !vmom
    mem_loc = 6.*lijk0*diag_trap_levels*SIZE_OF_COMPLEX_MB
    !vmom_corr
    mem_loc = mem_loc + 4.*lijk0*SIZE_OF_COMPLEX_MB
    !upar_corr, qpar_corr, qperp_corr, NOT dummy_corr
    IF (equil_par_curr) mem_loc = mem_loc + 2.0 * lijk0 * SIZE_OF_COMPLEX_MB
    !corrections
    mem_loc = mem_loc + 3.*lijk0*ln0*SIZE_OF_COMPLEX_MB
    IF (n_fields .GT. 2) THEN
      mem_loc = mem_loc + 2.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB ! FLR corrs
      mem_loc = mem_loc + 3.0 * lijk0 * SIZE_OF_COMPLEX_MB ! additional moments
      mem_loc = mem_loc + 3.0 * lijk0 * SIZE_OF_COMPLEX_MB ! additional moment vmom_corrs
      mem_loc = mem_loc + 3.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB ! additional moment corrs
    END IF
    IF (equil_par_curr) mem_loc = mem_loc + 3.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB

    mem_est_diag_mom = mem_loc+mem_req_in

  End Function mem_est_diag_mom

!!!******************************************************************!!!

  SUBROUTINE initialize_diag_mom
    Integer :: n, ierr
    Real :: maxB_loc
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: rank, o, trap_level
#endif

    if(diag_trap_levels.eq.0) then
       n_moms = 6
    elseif(diag_trap_levels.eq.1) then
       n_moms  = (diag_trap_levels+2)*6
!WARNING: probably we need to modify the trapped/passing boundary for each x,y position
       if (sum(diag_Blev).eq.0.) then
          maxB_loc = maxval(geom%Bfield)
          call mpi_allreduce(maxB_loc,diag_Blev(0),1,MPI_REAL_TYPE,MPI_MAX, mpi_comm_z, ierr)
       endif
    else
       n_moms  = (diag_trap_levels+2)*6
    endif

    IF (n_fields .GT. 2) THEN
      n_moms = n_moms + 3 ! two additional moments for B_par contributions to fluxes
      n_moms_Bpar = 3
    ELSE
      n_moms_Bpar = 0
    END IF

    if(write_std) then
       ALLOCATE(MOMFILE(ln1:ln2))

       if (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
          ! Moment files for each species in a different file
          DO n=ln1,ln2
             call get_unit_nr(MOMFILE(n))

             OPEN(MOMFILE(n), file=trim(diagdir)//'/mom_'//&
                  &trim(spec(n)%name)//''//trim(file_extension),&
                  form='unformatted', status=filestat, position=filepos)  
          END DO
       endif
    end if

#ifdef WITHFUTILS
    if ((write_h5).and.(my_pev+my_pew+my_pespec).eq.0) then
       isnap_mom = 0
       allocate(fidmom_h5(0:n_spec-1))
       allocate(momtrp_label(0:diag_trap_levels+1))
       if (diag_trap_levels.gt.0) THEN
          momtrp_label(:) = '_trap'
          momtrp_label(diag_trap_levels) = '_pass'
          momtrp_label(diag_trap_levels+1) = '_FLR'
       else
          momtrp_label(:) = ''
       endif

       do n=0,n_spec-1
          call creatf(trim(diagdir)//'/mom_'//trim(spec(n)%name)//trim(file_extension)//'.h5', &
               fidmom_h5(n), "moments diagnostics", 'd', mpi_comm_xyz)
          call creatg(fidmom_h5(n), '/mom_'//trim(spec(n)%name))
          rank = 0
          call creatd(fidmom_h5(n), rank, dims, "/mom_"//trim(spec(n)%name)//"/time", "time")

          do o = 1, 6 + n_moms_Bpar
             do trap_level=0,diag_trap_levels
                call creatg(fidmom_h5(n), '/mom_'//trim(spec(n)%name)//'/'//&
                     &trim(mom_label(o))//trim(momtrp_label(trap_level)))
             enddo
             if (diag_trap_levels.gt.0) call creatg(fidmom_h5(n), '/mom_'//&
                  &trim(spec(n)%name)//'/'//trim(mom_label(o))//&
                  &trim(momtrp_label(diag_trap_levels+1)))
          enddo
        end do
    endif
#endif       

  END SUBROUTINE initialize_diag_mom

  !> MOM (moments) diagnostic: Writes 3D data (kx,ky,z) of (combinations of) velocity space moments to mom_<species name>.dat:
  !! 1.) n     / n_j0 rho_ref / L_ref
  !! 2.) Tpar  / T_j0 rho_ref / L_ref
  !! 3.) Tperp / T_j0 rho_ref / L_ref
  !! 4.) qpar+1.5upar / p_j0 c_ref rho_ref / L_ref
  !! 5.) qperp+1.upar / p_j0 c_ref rho_ref / L_ref
  !! 6.) upar  / c_ref rho_ref / L_ref
  SUBROUTINE diag_mom
    INTEGER :: k, l, m, n

    ! density etc. are 3D arrays, the last index separates the different trapped and 
    ! passing particles

    COMPLEX, DIMENSION(:,:,:,:,:),allocatable:: vmom
    COMPLEX, DIMENSION(:,:,:,:),allocatable :: vmom_corr !f1 independent part of the full vel. sp. moment
    COMPLEX, DIMENSION(:,:,:),allocatable :: dummy_corr
    REAL :: energy
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: f_av
    INTEGER :: i, j, o, n_corrs, write_pe, trap_level
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2,1:n_flr_corrs) :: corrections
#ifdef WITHFUTILS
    complex, dimension(:,:,:,:,:),allocatable:: tmp_buf,tmpcorr_buf
    character(len=FILENAME_MAX) :: dset_name_mom
    integer:: ng, pes, ierr
#endif

    !PERFON('diag_mom')

    if (equil_par_curr) THEN
       ! still not implemented for B_parallel!
       n_corrs = 6 !all vmoms need a correction
    ELSE
       ALLOCATE(dummy_corr(li1:li2,lj1:lj2,lk1:lk2))
       dummy_corr = 0.0
       n_corrs = 3 !the last 3 vmoms are zero anyway
       IF (n_fields .GT. 2) n_corrs = 6
    ENDIF

    write_pe = pexyzvwspec(0,0,0,0,0,my_pespec)

    allocate(vmom(li1:li2, lj1:lj2, lk1:lk2, 6+n_moms_Bpar, 0:diag_trap_levels))
    allocate(vmom_corr(li1:li2, lj1:lj2, lk1:lk2, n_corrs))


#ifdef WITHFUTILS
    if(write_h5) then
       allocate(tmp_buf(li1:li2,lj1:lj2,lk1:lk2,6+n_moms_Bpar,0:n_procs_s-1))
       allocate(tmpcorr_buf(li1:li2,lj1:lj2,lk1:lk2,n_corrs,0:n_procs_s-1))
    end if
#endif

    do n=ln1,ln2
       if(write_std) then
          IF (mype.eq.write_pe) WRITE(MOMFILE(n)) time
       end if
       vmom = 0.0
       if (xy_local) then
#if WITHOMP_DIAG
          !$omp parallel do default(none) &
          !$omp private(m,l,k,trap_level,f_av,energy) &
          !$omp shared(f_,jfac,ll1,ll2,lk1,lk2,lm1,lm2,mu,geom) &
          !$omp shared(diag_trap_levels,diag_Blev,vp,li1,li2,lj1,pi1,pj1,n) &
          !$omp shared(mat_00,mat_10,mat_20,mat_30,mat_01,mat_11,lij0) &
          !$omp shared(nmat_00,nmat_20,nmat_01,n_fields) &
          !$omp shared(vmom)
#endif
          do k=lk1,lk2
             do m=lm1,lm2
                do l=ll1,ll2
                   energy = vp(l)*vp(l)+mu(m)*geom%Bfield(pi1,pj1,k)
                   !gyroaveraging
                   CALL gyro_op_wrapper(f_(:,:,k,l,m,n),f_av,k,m,n,0)
                   !>00mode is not zeroed out for mom data
                   !if ((istep_neoclass .gt. 0).and.(p_has_00_mode))&
                   !   & f_av(li1,li2)=0.0
                   do trap_level=0,diag_trap_levels
                      if ((energy.Lt.mu(m)*diag_Blev(trap_level)).or.(trap_level.eq.diag_trap_levels)) then
                         call axpy_ij(lij0,mat_00(pi1,pj1,k,l,m),  f_av,vmom(:,:,k,1,trap_level)) !dens
                         call axpy_ij(lij0,mat_20(pi1,pj1,k,l,m),  f_av,vmom(:,:,k,2,trap_level)) !Tpar
                         call axpy_ij(lij0,mat_01(pi1,pj1,k,l,m),  f_av,vmom(:,:,k,3,trap_level)) !Tperp
                         call axpy_ij(lij0,mat_30(pi1,pj1,k,l,m,n),f_av,vmom(:,:,k,4,trap_level)) !qpar
                         call axpy_ij(lij0,mat_11(pi1,pj1,k,l,m,n),f_av,vmom(:,:,k,5,trap_level)) !qperp
                         call axpy_ij(lij0,mat_10(pi1,pj1,k,l,m,n),f_av,vmom(:,:,k,6,trap_level)) !upar
                         exit
                      endif
                   enddo
                   IF (n_fields .GT. 2) THEN
                     !gyroaveraging for Bpar cases
                     CALL gyro_op_wrapper(f_(:,:,k,l,m,n),f_av,k,m,n,1)
                     do trap_level=0,diag_trap_levels
                       if ((energy.Lt.mu(m)*diag_Blev(trap_level)).or.(trap_level.eq.diag_trap_levels)) then
                         CALL axpy_ij(lij0,nmat_00(pi1,pj1,k,l,m),f_av,vmom(:,:,k,7,trap_level)) !densI1
                         CALL axpy_ij(lij0,nmat_20(pi1,pj1,k,l,m),f_av,vmom(:,:,k,8,trap_level)) !TparI1
                         CALL axpy_ij(lij0,nmat_01(pi1,pj1,k,l,m),f_av,vmom(:,:,k,9,trap_level)) !TperpI1
                         exit
                       endif
                     enddo
                   END IF
                enddo
             enddo
          enddo
#if WITHOMP_DIAG
          !$omp end parallel do
#endif
       else !with geometry and/or temp,dens profile information
#if WITHOMP_DIAG
          !$omp parallel do default(none) &
          !$omp private(m,l,k,j,i,trap_level,f_av,energy) &
          !$omp shared(f_,jfac,ll1,ll2,lk1,lk2,lm1,lm2,mu,geom,spec) &
          !$omp shared(diag_trap_levels,diag_Blev,vp,li1,li2,lj1,lj2,pi1,pj1,n) &
          !$omp shared(mat_00,mat_10,mat_20,mat_30,mat_01,mat_11,li0,lj0) &
          !$omp shared(nmat_00,nmat_20,nmat_01,n_fields) &
          !$omp shared(vmom)
#endif
          do k=lk1,lk2
             do m=lm1,lm2
                do l=ll1,ll2
                   !gyroaveraging
                   CALL gyro_op_wrapper(f_(:,:,k,l,m,n),f_av,k,m,n,0)
                   
                   do trap_level=0,diag_trap_levels
                      do i=li1,li2
                         !\todo: check the following energy x-dependence (temperature missing?)
                         energy = vp(l)*vp(l)+mu(m)*geom%Bfield(i,pj1,k)
                         if ((energy.Lt.mu(m)*diag_Blev(trap_level)).or.(trap_level.eq.diag_trap_levels)) then
                            do j=lj1,lj2
                               vmom(i,j,k,1,trap_level) = vmom(i,j,k,1,trap_level) + & !dens
                                    & mat_00(i,pj1,k,l,m) * f_av(i,j)
                               vmom(i,j,k,2,trap_level) = vmom(i,j,k,2,trap_level) + & !Tpar
                                    & mat_20(i,pj1,k,l,m) * f_av(i,j)
                               vmom(i,j,k,3,trap_level) = vmom(i,j,k,3,trap_level) + & !Tperp
                                    & mat_01(i,pj1,k,l,m) * f_av(i,j)
                               vmom(i,j,k,4,trap_level) = vmom(i,j,k,4,trap_level) + & !qpar
                                    & mat_30(i,pj1,k,l,m,n) * f_av(i,j)
                               vmom(i,j,k,5,trap_level) = vmom(i,j,k,5,trap_level) + & !qperp
                                    & mat_11(i,pj1,k,l,m,n) * f_av(i,j)
                               vmom(i,j,k,6,trap_level) = vmom(i,j,k,6,trap_level) + & !upar
                                    & mat_10(i,pj1,k,l,m,n) * f_av(i,j)/spec(n)%dens_prof(i)
                            enddo
                            !exit ???
                         endif
                      enddo
                   enddo
                   !B1_parallel is not implemented for nonlocal code versions
                enddo
             enddo
          enddo
#if WITHOMP_DIAG
          !$omp end parallel do          
#endif
       endif
       
       CALL my_complex_sum_vw(vmom, SIZE(vmom))
       
       if ((my_pev+my_pew).eq.0) then   
          ! calculating and writing FLR corrections
          IF (n_fields .LE. 2) THEN
             CALL correct_for_FLR(emfields(:,:,:,1), corrections)
          ELSE
             CALL correct_for_FLR_bpar(emfields(:,:,:,1),emfields(:,:,:,3),corrections)
          END IF

          ! dens_corr
          vmom_corr(:,:,:,1) = corrections(:,:,:,n,1)
          
          ! parallel & perpendicular temperature
          if (x_local) then
             !Tpar (wo FLR)
             vmom(:,:,:,2,:)  = 2*vmom(:,:,:,2,:) - vmom(:,:,:,1,:)
             vmom_corr(:,:,:,2) = 2*corrections(:,:,:,n,2)-corrections(:,:,:,n,1)
             !Tperp (wo FLR)  
             vmom(:,:,:,3,:) = vmom(:,:,:,3,:) - vmom(:,:,:,1,:)
             vmom_corr(:,:,:,3) = corrections(:,:,:,n,3)-corrections(:,:,:,n,1)

             IF (equil_par_curr) THEN
                !upar
                vmom_corr(:,:,:,6) = corrections(:,:,:,n,n_flr_corrs-2)
                !qpar_corr
                vmom_corr(:,:,:,4) = corrections(:,:,:,n,n_flr_corrs-1)
                !qperp_corr
                vmom_corr(:,:,:,5) = corrections(:,:,:,n,n_flr_corrs)
             END IF
          else ! with temp,dens profile information
             do i=li1,li2
                !Tpar (wo FLR)
                vmom(i,:,:,2,:) = (2*vmom(i,:,:,2,:)-spec(n)%temp_prof(i)*vmom(i,:,:,1,:))/&
                     &spec(n)%dens_prof(i)
                vmom_corr(i,:,:,2) = (2*corrections(i,:,:,n,2)-&
                     &spec(n)%temp_prof(i)*corrections(i,:,:,n,1))/spec(n)%dens_prof(i)
                !Tperp (wo FLR)
                vmom(i,:,:,3,:) = (vmom(i,:,:,3,:)-spec(n)%temp_prof(i)*vmom(i,:,:,1,:))/&
                     &spec(n)%dens_prof(i)
                vmom_corr(i,:,:,3) = (corrections(i,:,:,n,3)-&
                     &spec(n)%temp_prof(i)*corrections(i,:,:,n,1))/spec(n)%dens_prof(i)
             enddo
          endif

          IF (n_fields .GT. 2) THEN
            vmom_corr(:,:,:,4) = corrections(:,:,:,n,4)
            vmom_corr(:,:,:,5) = corrections(:,:,:,n,5)
            vmom_corr(:,:,:,6) = corrections(:,:,:,n,6)
          END IF

          if (diag_trap_levels.eq.0) then
             !add FLR corrections to moments
             do o=1,n_corrs
                IF (o .LE. 3) THEN
                  vmom(:,:,:,o,0) = vmom(:,:,:,o,0) + vmom_corr(:,:,:,o)
                ELSE
                  IF (n_fields .GT. 2) THEN
                    vmom(:,:,:,o+3,0) = vmom(:,:,:,o+3,0) + vmom_corr(:,:,:,o)
                  ELSE
                    vmom(:,:,:,o,0) = vmom(:,:,:,o,0) + vmom_corr(:,:,:,o)
                  END IF
                END IF
             enddo
          endif

          if(write_std) then
             do trap_level=diag_trap_levels,0,-1
                do o = 1, 6 + n_moms_Bpar
                   CALL write3d(MOMFILE(n),vmom(:,:,:,o,trap_level),write_pe)
                enddo
             enddo

             if (diag_trap_levels.gt.0) then !add FLR
                do o=1,n_corrs
                   CALL write3d(MOMFILE(n),vmom_corr(:,:,:,o),write_pe)
                enddo
                do o=n_corrs+1,6
                   ! no flr corrections, writing zeros for compability with IDL diagnostic
                   call write3d(MOMFILE(n),dummy_corr,write_pe)
                enddo
             endif
          endif
       end if

       if(write_std) then
          IF (mype.eq.write_pe) call flush(MOMFILE(n))
       end if

#ifdef WITHFUTILS
       if (write_h5.and.((my_pev+my_pew) .eq. 0))  then
          if(my_pespec .eq. 0) then
             do pes = 0, n_procs_s-1
                ng = pes*ln0 + n-ln1
                call append(fidmom_h5(ng), "/mom_"//trim(spec(ng)%name)//"/time", real(time,8))
                call attach(fidmom_h5(ng), "/mom_"//trim(spec(ng)%name)//"/time", "n_steps", isnap_mom+1)
             end do
          endif

          do trap_level=0,diag_trap_levels
             call mpi_gather(vmom(:,:,:,:,trap_level), (6+n_moms_Bpar)*lijk0, MPI_COMPLEX_TYPE, tmp_buf,&
                  &(6+n_moms_Bpar)*lijk0, MPI_COMPLEX_TYPE, 0, mpi_comm_spec, ierr)

             if(my_pespec .eq. 0) then
                do pes = 0, n_procs_s-1
                   ng = pes*ln0 + n-ln1
                   do o = 1, 6 + n_moms_Bpar
                      write(dset_name_mom, "(A, '/', i10.10)") "/mom_"&
                           &//trim(spec(ng)%name)//"/"//trim(mom_label(o))//&
                           &trim(momtrp_label(trap_level)), isnap_mom
                      call putarrnd(fidmom_h5(ng), trim(dset_name_mom),&
                           & tmp_buf(:,:,:,o,pes), (/3, 2, 1/))
                      call attach(fidmom_h5(ng), trim(dset_name_mom), "time", time)

                   enddo
                   call flushh5(fidmom_h5(ng))
                 end do
             end if
          enddo

          if (diag_trap_levels .gt. 0) then
             call mpi_gather(vmom_corr(:,:,:,:), n_corrs*lijk0, MPI_COMPLEX_TYPE, tmpcorr_buf,&
                  &n_corrs*lijk0, MPI_COMPLEX_TYPE, 0, mpi_comm_spec, ierr)
             
             if(my_pespec .eq. 0) then
                do pes = 0, n_procs_s-1
                   ng = pes*ln0 + n-ln1

                   do o = 1, n_corrs

                      write(dset_name_mom, "(A, '/', i10.10)") "/mom_"&
                           &//trim(spec(ng)%name)//"/"//trim(mom_label(o))//&
                           &trim(momtrp_label(diag_trap_levels+1)), isnap_mom
                      call putarrnd(fidmom_h5(ng), trim(dset_name_mom), tmpcorr_buf(:,:,:,o,pes), (/3, 2, 1/))
                      call attach(fidmom_h5(ng), trim(dset_name_mom), "time", time)
                   enddo

                   do o=n_corrs+1,6
                      ! no flr corrections, writing zeros for compability with IDL diagnostic
                      write(dset_name_mom, "(A, '/', i10.10)") "/mom_"&
                           &//trim(spec(ng)%name)//"/"//trim(mom_label(o))//&
                           &trim(momtrp_label(diag_trap_levels+1)), isnap_mom
                      call putarrnd(fidmom_h5(ng), trim(dset_name_mom), dummy_corr, (/3, 2, 1/))
                      call attach(fidmom_h5(ng), trim(dset_name_mom), "time", time)
                   enddo

                   call flushh5(fidmom_h5(ng))
                end do
             end if
          end if
       endif
#endif
    enddo
    
    deallocate(vmom, vmom_corr)
    IF (.not.equil_par_curr) DEALLOCATE(dummy_corr)
    
#ifdef WITHFUTILS
    if(write_h5) then
       isnap_mom = isnap_mom+1
       deallocate(tmp_buf,tmpcorr_buf)
    end if
#endif
    
    !PERFOFF
    
  END SUBROUTINE diag_mom

  Subroutine finalize_diag_mom
    integer :: n

    if(write_std) then
       if ((mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec))) then
          do n=ln1,ln2
             close(MOMFILE(n))
          enddo
       endif
       DEALLOCATE(MOMFILE)
    end if

#ifdef WITHFUTILS
    if (write_h5.and.(my_pev+my_pew+my_pespec).eq.0) then
       do n=0, n_spec-1
          call closef(fidmom_h5(n))
       end do
       deallocate(fidmom_h5,momtrp_label)
    end if
#endif
  End Subroutine finalize_diag_mom


!!!******************************************************************!!!
!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of diag_vsp
  Real Function mem_est_diag_vsp(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !local variables
    !var
    mem_loc = lklmn0*5.*SIZE_OF_REAL_MB
    !fullarr
    mem_loc = mem_loc + nz0*nv0*nw0*n_spec*5*SIZE_OF_REAL_MB
    !dfields_dy
    mem_loc = mem_loc + lij0*n_fields*SIZE_OF_COMPLEX_MB
    !phi_bar, temporary, dens_with_corr
    mem_loc = mem_loc + (2.*li0+li0)*lj0*SIZE_OF_COMPLEX_MB

    mem_est_diag_vsp = mem_loc+mem_req_in

  End Function mem_est_diag_vsp

!!!******************************************************************!!!

  Subroutine initialize_diag_vsp
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: o, rank
#endif

    if (mype==0) then
       if(write_std) then    
          call get_unit_nr(VSPFILE)
          OPEN(VSPFILE, file=trim(diagdir)//&
               &'/vsp'//trim(file_extension), form='unformatted', &
               status=filestat, position=filepos)
       END IF

#ifdef WITHFUTILS
       if (write_h5) then
          isnap_vsp = 0
          call creatf(trim(diagdir)//'/vsp'//trim(file_extension)//'.h5', &
               fidvsp_h5, "Velocity space diagnostics", 'd')
          call creatg(fidvsp_h5, '/vsp')
          do o = 1,5
             call creatg(fidvsp_h5, '/vsp/'//trim(vsp_label(o)))
          enddo
          rank = 0
          call creatd(fidvsp_h5, rank, dims, "/vsp/time", "time")
       end if
#endif
    endif

  End Subroutine Initialize_diag_vsp

  !>VSP (velocity space) diagnostic
  !!
  !!Calculates guiding-center-fluxes G_es/G_em and Q_es/em
  !!and sqrt(<|f_|^2>) averaged over kx and ky
  Subroutine diag_vsp
    Real:: var(lk1:lk2, ll1:ll2, lm1:lm2,ln1:ln2, 5)
    Real:: fullarr(0:nz0-1, 0:nv0-1, 0:nw0-1, 0:n_spec-1, 5)
    Integer :: k, l, m, n, o, i, ii, pni, ierr
    Real :: fnorm
    Complex, Dimension(li1:li2, lj1:lj2,1:n_fields):: dfields_dy
    COMPLEX, DIMENSION(li1:li2, lj1:lj2):: phi_bar,temporary
    ! dens_with_corr: density with FLR correction
    Complex, Dimension(li1:li2, lj1:lj2):: dens_with_corr
#ifdef WITHFUTILS
    character(len=FILENAME_MAX):: dset_name_vsp
#endif

    !PERFON('diag_vsp')
    pni=pn1
    do n=ln1,ln2
       if (pn0.gt.1) pni=n
       do m=lm1,lm2
          do l=ll1,ll2
             do k = lk1, lk2
                !!************************************************************!
                !!******************* calculate bx and vx ********************!
                !!************************************************************!
                do o=1,n_fields
                   call calc_dfields_dy(dfields_dy(:,:,o),k,o)
                enddo
                
                CALL gyro_op_wrapper(emfields(:,:,k,1),phi_bar,k,m,n,0)
                
                if (xy_local) then
                   temporary = f_(:,:,k,l,m,n)+spec(n)%charge/spec(n)%temp*&
                        & phi_bar*fm(pi1,pj1,k,l,m,pni)
                else
                   do i=li1,li2
                      temporary(i,:) = f_(i,:,k,l,m,n)+spec(n)%charge/spec(n)%temp*&
                           & phi_bar(i,:)*fm(i,pj1,k,l,m,n)                
                   enddo
                endif
                
                Call gyro_op_wrapper(temporary,dens_with_corr,k,m,n,0)
                
                if (xy_local) then
                   call axpy_ij(lij0,-spec(n)%charge/spec(n)%temp*fm(pi1,pj1,k,l,m,pni),&
                        emfields(:,:,k,1),dens_with_corr)
                else
                   do i=li1,li2
                      dens_with_corr(i,:) = dens_with_corr(i,:) - &
                           & spec(n)%charge/spec(n)%temp*fm(i,pj1,k,l,m,n)*&
                           & emfields(i,:,k,1)
                   enddo
                endif
                
                !             the statements above calculate:
                !             dens_with_corr = jfac(:,:,k,m,n)*f_(li1:li2,:,k,l,m,n)-&
                !                fm(k,l,m)*(1-jfac(:,:,k,m,n)**2)*phi(li1:li2,:,k)*&
                !                spec(n)%charge/spec(n)%temp
                
                !!************************************************************!
                !!******************* calculate transport ********************!
                !!************************************************************!
                !... G_es = <vx ( n + FLR corr)>
                var(k,l,m,n,1) = -2.0*Sum(&
                     Real(Conjg(dfields_dy(:,:,1))*dens_with_corr))
                
                !... G_em = <bx u_parallel>
                if (n_fields.gt.1) then
                   var(k,l,m,n,2) = 2.0*Sum(&
                        Real(Conjg(dfields_dy(:,:,2))*dens_with_corr))
                else
                   var(k,l,m,n,2) = 0
                endif
                
                !... <f_>
                var(k,l,m,n,5) = 2.0*Sum(&
                     Real(Conjg(f_(li1:li2,lj1:lj2,k,l,m,n))*f_(li1:li2,lj1:lj2,k,l,m,n)))
                
                !substract kj=0 mode (it should not be counted twice)
                if (p_has_0_mode) then
                   if (xy_local.and.yx_order) then !here ky=0 is at i=li1
                      var(k,l,m,n,5) = var(k,l,m,n,5)&
                           -Real(Sum(Conjg(f_(li1,:,k,l,m,n)) * f_(li1,:,k,l,m,n)))
                   else 
                      var(k,l,m,n,5) = var(k,l,m,n,5)&
                           -Real(Sum(Conjg(f_(li1:li2,lj1,k,l,m,n)) * f_(li1:li2,lj1,k,l,m,n)))
                   endif
                endif
             End Do
          End Do
       End Do
    End Do
       
    !... Q_es = <(3/2) p vx>
    var(:,:,:,:,3) = var(:,:,:,:,1)
    
    !... Q_em = <[(5/2) vpl + qpl + qpp] bx>
    var(:,:,:,:,4) = var(:,:,:,:,2)

    !   Calculate sum over n_procs_x, leave result in pex = 0.
    Call my_sum_to_0_real(var, Size(var), mpi_comm_x)
    !   Calculate sum over n_procs_y, leave result in pey = 0.
    Call my_sum_to_0_real(var, Size(var), mpi_comm_y)
!WARNING need to implement (x,y) dep. in the following part :
    do n=ln1,ln2
       var(:,:,:,n,1) = var(:,:,:,n,1) * mat_00(pi1,pj1,:,:,:) * spec(n)%dens
       var(:,:,:,n,2) = var(:,:,:,n,2) * mat_10(pi1,pj1,:,:,:,n) * spec(n)%dens
       var(:,:,:,n,3) = var(:,:,:,n,3) * (mat_20(pi1,pj1,:,:,:) + &
            &mat_01(pi1,pj1,:,:,:)) * spec(n)%dens * spec(n)%temp
       var(:,:,:,n,4) = var(:,:,:,n,4) * (mat_30(pi1,pj1,:,:,:,n) + &
            &mat_11(pi1,pj1,:,:,:,n)) * spec(n)%dens * spec(n)%temp
       var(:,:,:,n,5) = Sqrt( var(:,:,:,n,5) )
    enddo
    
    fullarr = 0.0

    Do ii = 1,5
       Do n = ln1,ln2
          Do m = lm1,lm2
             Do l = ll1,ll2
                CALL mpi_gather(var(lk1,l,m,n,ii),lk0, MPI_REAL_TYPE,&
                     fullarr(0,l-ll1,m-lm1,n-ln1,ii),lk0,MPI_REAL_TYPE,&
                     0, mpi_comm_z, ierr)
             Enddo
             CALL my_real_gather_to_0(fullarr(:,:,m-lm1,n-ln1,ii),&
                  &1,nz0*ll0,nz0*nv0,mpi_comm_v)
          Enddo
          CALL my_real_gather_to_0(fullarr(:,:,:,n-ln1,ii),&
               &1,nz0*nv0*lm0,nz0*nv0*nw0,mpi_comm_w)
       Enddo
       CALL my_real_gather_to_0(fullarr(:,:,:,:,ii),&
            &1,nz0*nv0*nw0*ln0,nz0*nv0*nw0*n_spec,mpi_comm_spec)
    Enddo

    if (xy_local) then
       fnorm = 1.0
    else
       fnorm = 1.0 / REAL(ni0)
    endif 
    If (mype == 0) Then
       fullarr = fullarr * fnorm
       if(write_std) then
          Write(VSPFILE) time
          Write(VSPFILE) fullarr
         call flush(VSPFILE)
       end if
    End If

#ifdef WITHFUTILS
    if(write_h5) then
      if(mype.eq.0) then          
          call append(fidvsp_h5, "/vsp/time", real(time,8))
          call attach(fidvsp_h5, "/vsp/time", "n_steps", isnap_vsp+1)
          do o = 1,5
             write(dset_name_vsp, "(A, '/', i10.10)") "/vsp/"//&
                  &trim(vsp_label(o)), isnap_vsp
             call putarr(fidvsp_h5, dset_name_vsp, fullarr(:,:,:,:,o)*fnorm)
             call attach(fidvsp_h5, dset_name_vsp, "time", time)             
          enddo

          call flushh5(fidvsp_h5)
       end if
      isnap_vsp = isnap_vsp+1
   end if
#endif

    Call my_barrier()
    !PERFOFF
  End Subroutine diag_vsp

  Subroutine finalize_diag_vsp
    If (mype==0) then
       if (write_std) CLOSE(VSPFILE)
#ifdef WITHFUTILS
       if (write_h5) call closef(fidvsp_h5)
#endif
    endif
  End Subroutine finalize_diag_vsp


!!!******************************************************************!!!
!!!***************** additional subroutines  **** *******************!!!
!!!******************************************************************!!!
  
  
  !> Calculates the y-derivative of the electromagnetic fields
  !! for computing the ExB drift velocity that enters the fluxes
  Subroutine calc_dfields_dy(dfields_dy,k,o)
    complex, dimension(li1:li2,lj1:lj2),intent(out):: dfields_dy
    integer,intent(in):: o, k
    complex, dimension(lbi:ubi,lj1:lj2):: field_wb
    integer:: j

      if (yx_order) then
         if (y_local) then
            do j=lj1,lj2
               dfields_dy(:,j)=imag*ki(:)*emfields(:,j,k,o)*flux_geomfac(pi1,pj1,k)
            end do
         else
            !for .not. y_local: compute the y derivative, then multiply geom_fac
            field_wb(li1:li2,:) = emfields(:,:,k,o)
            call x_deriv_exc(field_wb,dfields_dy) 
            do j=lj1,lj2
               dfields_dy(:,j)=dfields_dy(:,j)*flux_geomfac(:,pj1,k)
            enddo
         end if
      else
         do j=lj1,lj2
            dfields_dy(:,j)=imag*kj(j)*emfields(:,j,k,o)
            if (x_local) then
               dfields_dy(:,j) = dfields_dy(:,j)*flux_geomfac(pi1,pj1,k)
            else
               dfields_dy(:,j) = dfields_dy(:,j)*flux_geomfac(:,pj1,k)
            endif
         end do
      end if

  End Subroutine calc_dfields_dy


!!!******************************************************************!!!
!!!******************************************************************!!!

  !> Rescales the distribution function in order to avoid overflows;
  !! only used in linear mode
  subroutine rescale_g_1(overflow_exit)
    Logical, intent(inout) :: overflow_exit
    integer::ierr
    real:: locsum,globsum

    locsum=sum(abs(g_1))

    Call mpi_allreduce(locsum,globsum,1,MPI_REAL_TYPE,&
         MPI_SUM, mpi_comm_vwspec, ierr) 
    locsum=globsum
    Call mpi_allreduce(locsum,globsum,1,MPI_REAL_TYPE,&
         MPI_SUM, mpi_comm_z, ierr) 
    locsum=globsum
    Call mpi_allreduce(locsum,globsum,1,MPI_REAL_TYPE,&
         MPI_SUM, mpi_comm_x, ierr) 

    
    
    overflow_exit = (globsum.NE.globsum) ! check for NaN, NaNQ

    IF (.NOT.overflow_exit) THEN
       g_1=g_1/globsum
       if (istep_omega.ne.0) call rescale_phi_old(globsum)
       
       if ((mype.eq.0).and.print_ini_msg) &
            & write(*,"(a)") 'rescaled amplitudes'
       call calc_aux_fields(g_1,emfields,f_,.false.)
    ENDIF
  end subroutine rescale_g_1


!!!******************************************************************!!!
!!!******************************************************************!!!

  !>Initializes diag_df_ev
  subroutine initialize_diag_df_ev
    implicit none
    integer :: n,info_handle,ierr
    character(len=4) :: file_ky_num,file_kx_num
    integer :: i,j,idkx,ind,ncon_tot,ncon_left,current_shift
    integer :: mloc(1)
    integer :: kx_shift_index
    integer, allocatable, dimension(:) :: ind_temp
    integer :: pb_xshift_allky(0:nky0-1),pbtemp(0:nky0-1)
    character(len=FILENAME_MAX) :: df_path

    ALLOCATE(DF_EV_FILE(num_ky_modes))
    pb_xshift_allky=0
    if(ubexc.gt.lh1) pb_xshift_allky(lh1:ubexc)=pb_xshift(lh1:ubexc) 
    call MPI_ALLREDUCE(pb_xshift_allky,pbtemp,nky0,MPI_INTEGER,MPI_SUM,mpi_comm_y,ierr) 
    pb_xshift_allky=pbtemp


    if(mype==0) then
      call get_unit_nr(info_handle)
      open(unit=info_handle,file=trim(diagdir)//'/dfout.info',status='unknown')
    end if

    allocate(nkx_keep_array(num_ky_modes))  
    nkx_keep_array=0
    !allocate with nx0 so that it has space for maximum
    allocate(kx_indices(nx0,num_ky_modes))  
    kx_indices=0

    do i=1,num_ky_modes

      !idkx=nexc*which_ky(i) !number kx indices to skip for parallel boundary condition 
      idkx=pb_xshift_allky(which_ky(i)) 
      ncon_tot=1

      ind=which_kx_center(i)+idkx
      do while((ind.lt.nx0/2.and.ind.gt.-nx0/2).and.which_ky(i).ne.0.and.idkx.ne.0)
      !Add one to the number of kx connections until kx_max is reached
        ncon_tot=ncon_tot+1
        ind=ind+idkx 
      end do
      ind=which_kx_center(i)-idkx
      ncon_left=0
      do while((ind.gt.-nx0/2.and.ind.lt.nx0/2).and.which_ky(i).ne.0.and.idkx.ne.0)
      !Do the same in the negative direction
        ncon_tot=ncon_tot+1
        ind=ind-idkx 
        ncon_left=ncon_left+1
      end do

      !Take the maximum number of modes such that there are the same number
      !on each side of kx_center
      ncon_left=min(ncon_left,ncon_tot-1-ncon_left)
      ncon_tot=ncon_left*2+1
      if(which_ky(i)==0) then
         ncon_tot=1      !If this is a zonal mode keep only one kx
      end if
      if(.not.nonlinear.and.adapt_lx) ncon_tot=nx0    !If this is a linear run, keep all kx
      nkx_keep_array(i)=min(ncon_tot,num_kx_modes)
      !if(allocated(kx_indices)) deallocate(kx_indices)
      if(allocated(ind_temp)) deallocate(ind_temp)
      allocate(ind_temp(nkx_keep_array(i)))

      if(mype==0) then
         write(info_handle,*) "************************************************************"
         write(info_handle,*) "************************************************************"
         write(info_handle,*) "ky index",  which_ky(i)
         write(info_handle,*) "kx center", which_kx_center(i)
         write(info_handle,*) "kx shift",  pb_xshift_allky(which_ky(i))
         write(info_handle,*) "nkx_keep",  nkx_keep_array(i)
      end if

      if(mype==0) write(info_handle,*) "kx modes:"

      if(nkx_keep_array(i)==1) then
        kx_indices(1,i)=which_kx_center(i)
        if(mype==0) write(info_handle,*) kx_indices(1,i)
      else
        current_shift=pb_xshift_allky(which_ky(i)) 
        do j=-(nkx_keep_array(i)-1)/2,(nkx_keep_array(i)-1)/2
          kx_shift_index=which_kx_center(i)+j*current_shift
           if(kx_shift_index.lt.0) then
             kx_shift_index=nx0+kx_shift_index
           end if 
           if(.not.nonlinear.and.adapt_lx) kx_shift_index=j+(nkx_keep_array(i)-1)/2
           kx_indices(j+(nkx_keep_array(i)-1)/2+1,i)=kx_shift_index
        end do
        ind_temp(:)=kx_indices(1:nkx_keep_array(i),i)
        do j=1,nkx_keep_array(i)
          mloc=minloc(ind_temp)
          kx_indices(j,i)=ind_temp(mloc(1))
          ind_temp(mloc(1))=nx0+100
          if(mype==0) write(info_handle,*) kx_indices(j,i)
        end do
      end if
    end do

    if(mype==0) close(info_handle)

    if(dfout_mpiio) then

      allocate(df_offset(num_ky_modes))
      df_offset=0
      allocate(dfarr_of_sizes(6,num_ky_modes))
      allocate(dfarr_of_subsizes(6,num_ky_modes))
      allocate(dfarr_of_starts(6,num_ky_modes))
      allocate(dfarr(num_ky_modes))
      call mpi_info_create(dfinfo,ierr)
      call mpi_info_set(dfinfo,"romio_ds_write", "disable",ierr)
      call mpi_info_set(dfinfo,"romio_ds_read", "disable",ierr)
      do i=1, num_ky_modes
        if(which_ky(i).ge.lj1.and.which_ky(i).le.lj2) then
          dfarr_of_sizes(:,i)=(/nkx_keep_array(i),1,nz0,nv0,nw0,n_spec/)
          dfarr_of_subsizes(:,i)=(/nkx_keep_array(i),1,lk0,ll0,lm0,ln0/)
          dfarr_of_starts(:,i)=(/0,0,my_pez*lk0,my_pev*ll0,&
               &my_pew*lm0,my_pespec*ln0/)
          call mpi_type_create_subarray(6,dfarr_of_sizes(:,i),dfarr_of_subsizes(:,i),&
              &dfarr_of_starts(:,i),MPI_ORDER_FORTRAN,MPI_COMPLEX_TYPE,dfarr(i),ierr)
          call mpi_type_commit(dfarr(i),ierr)
  
            write(file_ky_num,"(i4.4)") which_ky(i)  
            write(file_kx_num,"(i4.4)") abs(which_kx_center(i)) 
            if(which_kx_center(i).ge.0) then
               df_path=trim(diagdir)//'/'//'df_ky'//&
                    &file_ky_num//'kx'//file_kx_num//'.dat'
            else
               df_path=trim(diagdir)//'/'//'df_ky'//&
                    &file_ky_num//'kxn'//file_kx_num//'.dat'
            end if
          call mpi_file_open(mpi_comm_zvwspec,df_path,MPI_MODE_WRONLY+MPI_MODE_CREATE,dfinfo,DF_EV_FILE(i),ierr)
        end if !correct ky process
      end do !num_ky_modes
    else !not dfout_mpiio
      do n=1,num_ky_modes
       call get_unit_nr(DF_EV_FILE(n))
          write(file_ky_num,"(i4.4)") which_ky(n)  
          write(file_kx_num,"(i4.4)") abs(which_kx_center(n)) 
          if(which_kx_center(n).ge.0) then
             open(unit=DF_EV_FILE(n),file=trim(diagdir)//'/'//'df_ky'//&
                  &file_ky_num//'kx'//file_kx_num//'.dat',form='unformatted',&
                  &status='unknown')    
          else
             open(unit=DF_EV_FILE(n),file=trim(diagdir)//'/'//'df_ky'//&
                  &file_ky_num//'kxn'//file_kx_num//'.dat',form='unformatted',&
                  &status='unknown')    
          end if
      end do
    end if  !dfout_mpiio

  end subroutine initialize_diag_df_ev

!> This subroutine outputs certain kx/ky's of the distribution function.  Useful when outputting the entire 
!!distribution function is not feasible due to memory constraints.  Also used in conjunction with extended
!!diagnostics
  subroutine diag_df_ev
    
    implicit none
    
    integer :: n,m,l,ierr,i,j
    
    complex, allocatable, dimension(:,:,:,:,:):: g_full
    complex, dimension(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2):: g_loc
    complex, dimension(0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1) :: g_glob
!    integer, allocatable, dimension(:) :: kx_indices,ind_temp
!    integer :: current_shift
    

    do i=1,num_ky_modes

       if((which_ky(i).ge.lj1).and.(which_ky(i).le.lj2)) then
          
         ! write(*,*) "nkx_keep_array(i)",nkx_keep_array(i)
          allocate(g_full(nkx_keep_array(i),nz0,nv0,nw0,n_spec))
          
          do j=1,nkx_keep_array(i)
             g_loc=g_1(kx_indices(j,i),which_ky(i),:,:,:,:)
             
             Do n = ln1,ln2
                Do m = lm1,lm2
                   Do l = ll1,ll2
                      CALL mpi_gather(g_loc(lk1,l,m,n),lk0, MPI_COMPLEX_TYPE,&
                           g_glob(0,l-ll1,m-lm1,n-ln1),lk0,MPI_COMPLEX_TYPE,&
                           0, mpi_comm_z, ierr)
                   Enddo

                   CALL my_complex_gather_to_0(g_glob(:,:,m-lm1,n-ln1),&
                        &1,nz0*ll0,nz0*nv0,mpi_comm_v)
                Enddo

                CALL my_complex_gather_to_0(g_glob(:,:,:,n-ln1),&
                     &1,nz0*nv0*lm0,nz0*nv0*nw0,mpi_comm_w)
             Enddo
             
             CALL my_complex_gather_to_0(g_glob(:,:,:,:),&
                  &1,nz0*nv0*nw0*ln0,nz0*nv0*nw0*n_spec,mpi_comm_spec)

             g_full(j,:,:,:,:)=g_glob 
             
             
          end do !end kx_index loop
          if((my_pez==0).and.(my_pev==0).and.(my_pew==0).and.(my_pespec==0)) &
               write(DF_EV_FILE(i)) time
          
          if((my_pez==0).and.(my_pev==0).and.(my_pew==0).and.(my_pespec==0)) &
               write(DF_EV_FILE(i)) g_full
          
          deallocate(g_full)
       end if !ky if statement  
       call my_barrier
       
    end do !ky loop
    
  end subroutine diag_df_ev


!> This subroutine outputs certain kx/ky's of the distribution function.  Useful when outputting the entire 
!!distribution function is not feasible due to memory constraints.  Also used in conjunction with extended
!!diagnostics.  This version uses MPI-IO
  subroutine diag_df_ev_mpi
    
    implicit none
    complex, allocatable,dimension(:,:,:,:,:,:) :: g_
    integer :: i,j,ierr,n_procs_dfout
    
    n_procs_dfout=n_procs_z*n_procs_v*n_procs_w*n_procs_s

    do i=1,num_ky_modes

       if((which_ky(i).ge.lj1).and.(which_ky(i).le.lj2)) then

         allocate(g_(0:nkx_keep_array(i)-1,0:0,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
          
          do j=0,nkx_keep_array(i)-1
             g_(j,0,:,:,:,:)=g_1(kx_indices(j+1,i),which_ky(i),:,:,:,:)
          end do !end kx_index loop

         !write header
         call mpi_file_set_view(DF_EV_FILE(i),df_offset(i),MPI_BYTE,MPI_BYTE,"native",dfinfo,ierr)
         if (mype.eq.0) then
           call mpi_file_write(DF_EV_FILE(i),time,1,MPI_REAL_TYPE,MPI_STATUS_IGNORE,ierr)
         endif

         df_offset(i)=df_offset(i)+4
         if (prec.eq.'DOUBLE') df_offset(i)=df_offset(i)+4
  
         call mpi_file_set_view(DF_EV_FILE(i),df_offset(i),MPI_BYTE,dfarr(i),"native",dfinfo,ierr)
         call mpi_file_write_all(DF_EV_FILE(i),g_(0,0,lk1,ll1,lm1,ln1),nkx_keep_array(i)*1*lk0*ll0*lm0*ln0,&
              MPI_COMPLEX_TYPE,MPI_STATUS_IGNORE,ierr)
         if (prec.eq.'DOUBLE') then
            df_offset(i)=df_offset(i)+nkx_keep_array(i)*1*lk0*ll0*lm0*ln0*n_procs_dfout*16
         else
            df_offset(i)=df_offset(i)+nkx_keep_array(i)*1*lk0*ll0*lm0*ln0*n_procs_dfout*8
         end if

         deallocate(g_)

       end if !ky if statement  
       
    end do !ky loop
    
  end subroutine diag_df_ev_mpi

  subroutine finalize_diag_df_ev
    integer :: n
    
    if(.not.dfout_mpiio) then
      do n=1,num_ky_modes
        close(DF_EV_FILE(n))
      end do
    end if

    DEALLOCATE(DF_EV_FILE)
    deallocate(nkx_keep_array)
    deallocate(kx_indices)
    !reuse DF_EV_FILE(1) handle
    !open(unit=DF_EV_FILE(1),file='dfout.info',status='old',position='append')
    !write(DF_EV_FILE(1),"(A,G12.4)") "df_n_time=",itime/istep_dfout 
    !close(DF_EV_FILE(1))

    if(dfout_mpiio) then
      deallocate(df_offset)
      deallocate(dfarr_of_sizes)
      deallocate(dfarr_of_subsizes)
      deallocate(dfarr_of_starts)
    end if

  end subroutine finalize_diag_df_ev

  subroutine ev_out_open(ev_handle,r_or_l)
   integer, intent(out) :: ev_handle
   integer, intent(in) :: r_or_l
   character(len=FILENAME_MAX) :: filename , fullpath
 
     call get_unit_nr(ev_handle)
     filename='eigenvectors'//trim(file_extension)
     if(r_or_l.eq.1) then
       filename='eigenvectors_r'//trim(file_extension)
     end if
     if(r_or_l.eq.2) then
       filename='eigenvectors_l'//trim(file_extension)
     end if
       fullpath=trim(diagdir)//'/'//trim(filename)
     open(unit=ev_handle,file=fullpath,form='unformatted',status='unknown')
     if(vec_out_limit) then
       call get_unit_nr(par_handle)
       open(unit=par_handle,file=trim(diagdir)//'/ev_z_structure.dat',status='unknown')
     end if
     
  end subroutine ev_out_open
  
  subroutine ev_out_write(g_,ev_number,numb_ev,ev_handle)
    integer, intent(in) :: ev_number, numb_ev,ev_handle
    integer :: n,m,l,ierr,i
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), intent(in):: g_
    complex:: g_full(0:nx0-1,0:nky0-1,0:nz0-1, 0:nv0-1, 0:nw0-1, 0:n_spec-1)
    complex:: g_glob(0:nky0-1,0:nz0-1, 0:nv0-1, 0:nw0-1, 0:n_spec-1)
    complex, dimension(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2):: g_loc
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields):: ev_field
    complex,dimension(lk1:lk2,1:n_fields):: ev_field_loc
    complex,dimension(0:nz0-1,1:n_fields):: ev_field_glob
    complex,dimension(li1:li2,lj1:lj2,0:nz0-1,1:n_fields):: ev_field_full
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2) :: f_in
    
    integer :: z_ind,x_ind
    
    do i=0,nx0-1
       g_loc=g_(i,lj1,:,:,:,:)
       Do n = ln1,ln2
          Do m = lm1,lm2
             Do l = ll1,ll2
                CALL mpi_gather(g_loc(lk1,l,m,n),lk0, MPI_COMPLEX_TYPE,&
                     g_glob(0,0,l-ll1,m-lm1,n-ln1),lk0,MPI_COMPLEX_TYPE,&
                     0, mpi_comm_z, ierr)
             Enddo

             CALL my_complex_gather_to_0(g_glob(:,:,:,m-lm1,n-ln1),&
                  &1,nky0*nz0*ll0,nky0*nz0*nv0,mpi_comm_v)

          Enddo
          
          CALL my_complex_gather_to_0(g_glob(:,:,:,:,n-ln1),&
                  &1,nky0*nz0*nv0*lm0,nky0*nz0*nv0*nw0,mpi_comm_w)
       Enddo
       
       CALL my_complex_gather_to_0(g_glob(:,:,:,:,:),&
            &1,nky0*nz0*nv0*nw0*ln0,nky0*nz0*nv0*nw0*n_spec,mpi_comm_spec)
       
       g_full(i,:,:,:,:,:)=g_glob
    end do
    
    if(mype==0) write(ev_handle) ev_number
    if((mype==0).and.(mod(ev_number,100)==0)) write(*,*) "Eigenvector written.", ev_number, " of ", numb_ev
    if(mype==0) write(ev_handle) g_full
    call my_barrier()
    
    if(vec_out_limit.and.(ev_number.le.n_ev.or.which_ev.ne.'all_mpl')) then 
       !output parallel mode structure of top ten modes in gnuplotable format
       call calc_aux_fields(g_,ev_field,f_in,.false.)
       do i=0,nx0-1
          ev_field_loc=ev_field(i,lj1,lk1:lk2,:)
          do n=1,n_fields
             call mpi_gather(ev_field_loc(lk1,n),lk0,MPI_COMPLEX_TYPE,&
                  &ev_field_glob(0,n),lk0,MPI_COMPLEX_TYPE,&
                  &0,mpi_comm_z,ierr)
          end do
          ev_field_full(i,lj1,:,:)=ev_field_glob
       end do
       
       if(mype==0) then
          write(par_handle,*) ""
          write(par_handle,*) ""
          do x_ind=(nx0-1)/2+1,nx0-1
             do z_ind=0,nz0-1
                if(n_fields==1) then
                   write(par_handle,'(3es12.4)') real(conjg(ev_field_full(&
                        &x_ind,lj1,z_ind,1))*ev_field_full(x_ind,lj1,z_ind,1)),&
                        &real(ev_field_full(x_ind,lj1,z_ind,1)),&
                        &aimag(ev_field_full(x_ind,lj1,z_ind,1))
                else
                   write(par_handle,'(6es12.4)') real(conjg(ev_field_full(&
                        &x_ind,lj1,z_ind,1))*ev_field_full(x_ind,lj1,z_ind,1)),&
                        &real(ev_field_full(x_ind,lj1,z_ind,1)),&
                        &aimag(ev_field_full(x_ind,lj1,z_ind,1)),&
                        &real(conjg(ev_field_full(x_ind,lj1,z_ind,2))*&
                        &ev_field_full(x_ind,lj1,z_ind,2)),&
                        &real(ev_field_full(x_ind,lj1,z_ind,2)),&
                        &aimag(ev_field_full(x_ind,lj1,z_ind,2))
                end if
             end do
          end do
          do x_ind=0,(nx0-1)/2
             do z_ind=0,nz0-1
                if(n_fields==1) then
                   write(par_handle,'(3es12.4)') real(conjg(ev_field_full(&
                        &x_ind,lj1,z_ind,1))*ev_field_full(x_ind,lj1,z_ind,1)),&
                        &real(ev_field_full(x_ind,lj1,z_ind,1)),&
                        &aimag(ev_field_full(x_ind,lj1,z_ind,1))
                else
                   write(par_handle,'(6es12.4)') real(conjg(ev_field_full(&
                        &x_ind,lj1,z_ind,1))*ev_field_full(x_ind,lj1,z_ind,1)),&
                        &real(ev_field_full(x_ind,lj1,z_ind,1)),&
                        &aimag(ev_field_full(x_ind,lj1,z_ind,1)),&
                        &real(conjg(ev_field_full(x_ind,lj1,z_ind,2))*&
                        &ev_field_full(x_ind,lj1,z_ind,2)),&
                        &real(ev_field_full(x_ind,lj1,z_ind,2)),&
                        &aimag(ev_field_full(x_ind,lj1,z_ind,2))
                end if
             end do
          end do
       end if !mype==0
    end if
    
  end subroutine ev_out_write
  
  subroutine ev_out_close(ev_handle)
    integer, intent(in) :: ev_handle
    
    close(ev_handle)
    if(vec_out_limit) close(par_handle) 

  end subroutine ev_out_close

#if 0
!*************************************************************************************
  subroutine initialize_diag_g1

    Integer :: ierr
    
    call get_unit_nr(G1FILE)
    if (mype == 0) then
       !delete old file
       open(G1FILE, file=trim(diagdir)//'/g1'//trim(file_extension), &
            &status='replace')
       close(G1FILE,status='delete')
    endif
    call my_barrier()
    
    call mpi_type_create_subarray(5,(/nx0*nky0,nz0,nv0,nw0,n_spec/),&
         &(/li0*lj0,lk0,ll0,lm0,ln0/),&
         &(/my_pey*lj0*li0,my_pez*lk0,my_pev*ll0,my_pew*lm0,my_pespec*ln0/),&
         &MPI_ORDER_FORTRAN,MPI_COMPLEX_TYPE,g1arr,ierr)
    call mpi_type_commit(g1arr,ierr)
    
    call mpi_file_open(MY_MPI_COMM_WORLD,trim(diagdir)//'/g1'//&
         &trim(file_extension),MPI_MODE_WRONLY+MPI_MODE_CREATE,&
         &MPI_INFO_NULL,G1FILE,ierr)
    g1_offset=0
    
  end subroutine initialize_diag_g1

  !>Writes out the whole g1 distribution function
  subroutine diag_g1
    integer:: ierr

    !PERFON('diag_g1')
    !write header
    if (mype.eq.0) then
       call mpi_file_write(G1FILE,time,1,MPI_REAL_TYPE,MPI_STATUS_IGNORE,ierr)
    endif

    if (prec.eq.'DOUBLE') then 
       g1_offset=g1_offset+8
    else
       g1_offset=g1_offset+4
    endif
        
    call mpi_file_set_view(G1FILE,g1_offset,MPI_BYTE,g1arr,"native",&
         &MPI_INFO_NULL,ierr)
    call mpi_file_write_all(G1FILE,g_1(li1,lj1,lk1,ll1,lm1,ln1),li0*lj0*lk0*&
         &ll0*lm0*ln0,MPI_COMPLEX_TYPE,MPI_STATUS_IGNORE,ierr)

    if (prec.eq.'DOUBLE') then 
       g1_offset=g1_offset+nx0*nky0*nz0*nv0*nw0*n_spec*16
    else
       g1_offset=g1_offset+nx0*nky0*nz0*nv0*nw0*n_spec*8
    endif
    !PERFOFF

  end subroutine diag_g1

  subroutine finalize_diag_g1
    integer :: ierr
    
    call mpi_file_close(G1FILE,ierr)
    call mpi_type_free(g1arr,ierr)

  end subroutine finalize_diag_g1
#endif
!*************************************************************************************
!*************************************************************************************

  subroutine initialize_diag_gav

    if (gav_stime.ge.0.) &
         ALLOCATE(gav(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

    gav = 0.0

    old_gav_time = -1.0

  end subroutine initialize_diag_gav

  !>Computes time averaged distribution function
  subroutine diag_gav
    if ((gav_stime.ge.0.0).and.(gav_stime.le.time)) then
       if (old_gav_time.lt.0.0) then
          gav_stime = time !set gav_stime to actually used starting time
       else
          gav = gav+g_1*(time-old_gav_time)
       endif
       old_gav_time = time
    endif
  end subroutine diag_gav


  !>Writes out and finalizes time averaged distribution function
  subroutine finalize_diag_gav

    if (gav_stime.ge.0) then

       if (gav_stime.lt.old_gav_time) then
          gav = gav/(old_gav_time-gav_stime)
       else
          gav = g_1
       endif

       call checkpoint_write(gav,'gav')

       DEALLOCATE(gav)
    endif

  end subroutine finalize_diag_gav

  subroutine initialize_diag_antenna
    character(len=250) :: strarr !<Character array
    integer:: p
  
    If (mype.eq.0) then
       call get_unit_nr(antenna_file)
       open(antenna_file, file=trim(diagdir)//&
            &'/antenna'//trim(file_extension), form='formatted', &
            &status=filestat, position=filepos)
       Write(strarr,"(A)") "# time"
       Write(strarr,"(2A,3I3,A,3I3,A)") trim(strarr),"           ",lv_antenna_modes(1:3)," (real) ",&
            lv_antenna_modes(1:3)," (imag)"
       do p=2,n_antennae
          Write(strarr,"(2A,3I3,A,3I3,A)")  trim(strarr)," ",lv_antenna_modes(p*3-2:p*3),&
               " (real) ",lv_antenna_modes(p*3-2:p*3)," (imag)"
       enddo

       write(antenna_file,"(A)") trim(strarr)
    endif

  end subroutine initialize_diag_antenna
  
  subroutine diag_antenna
    character(len=250) :: strarr !<Character array
    integer:: p
   
    if (mype==0) then
       Write(strarr,"(ES16.8)") time
       do p=1,n_antennae
          Write(strarr,"(2A,ES16.8,A,ES16.8)") trim(strarr),' ',real(amp(p)),' ',aimag(amp(p))
       enddo
       
       write(antenna_file,"(A)") trim(strarr)
    endif

  end subroutine diag_antenna

  subroutine finalize_diag_antenna
    if (mype==0) close(antenna_file)
  end subroutine finalize_diag_antenna

!*************************************************************************************

  !>Small routine abbreviating the check for diagnostic
  !!calls and calc_aux_fields
  Logical function call_diag(itime,istep,no_aux_yet)
    Integer, intent(in) :: itime, istep
    Logical, intent(inout),optional :: no_aux_yet
    
    IF (istep.gt.0) THEN
       call_diag = (MODULO(itime,istep).eq.0) 
    ELSE
       call_diag = .false.
    ENDIF

    IF (call_diag) THEN
       IF (present(no_aux_yet)) THEN
          if (no_aux_yet) then
             call calc_aux_fields(g_1,emfields,f_,.false.)
             no_aux_yet=.false.
          endif
       ENDIF
    ENDIF

  END function call_diag


  !>Computes x,y,z sum including jacobian
  !!\param n_moms Number of moments/variables
  !!\param mom velocity space moments to be summed
  !!\param var resulting arrays
  subroutine sum_3d(n_moms,mom,var)
    integer, intent(in) :: n_moms
    Real,Dimension(li1:li2,lj1:lj2,lk1:lk2,n_moms,ln1:ln2),intent(in):: mom
    Real,Dimension(n_moms,ln1:ln2) :: var

    integer :: j, k, n, o
    Real :: fnorm

    !!************************************************************!
    !!************ Fourier space averaging ******* ***************!
    !!************************************************************!
    if (xy_local) then
       fnorm = 1.0 / (Real(nz0)*geom%avg_jaco)
       do n=ln1,ln2
#if WITHOMP_DIAG
          !$omp parallel do default(none) &
          !$omp private(o) &
          !$omp shared(var, mom, geom,pi1,pj1,n,n_moms)
#endif
          do o = 1, n_moms
             var(o,n) = sum( sum( Sum(mom(:,:,:,o,n),1),1 )*geom%jacobian(pi1,pj1,:) )  
          enddo
#if WITHOMP_DIAG
          !$omp end parallel do
#endif
       enddo
    else !nonlocal
       fnorm = 1.0 / (REAL(nz0)*real(ni0)*geom%avg_jaco)
       var = 0.
       do n=ln1,ln2
          do o = 1, n_moms
             do k=lk1,lk2
                do j=lj1,lj2
                   var(o,n) = var(o,n) + sum(mom(:,j,k,o,n)*geom%jacobian(pi1:pi2,pj1,k))
                enddo
             enddo
          enddo
       enddo
    endif

    Call my_sum_to_0_real(var, Size(var), mpi_comm_xyz)

    var=var*fnorm

  end subroutine sum_3d

End Module diagnostics
