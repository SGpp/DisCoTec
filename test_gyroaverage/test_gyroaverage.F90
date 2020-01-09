#include "redef.h"

PROGRAM test_gyroaverage_df
  USE localpolynombase_mod
  USE communications, ONLY: MY_MPI_COMM_WORLD,COMM_X, communicators,&
       &mpi_comm_x, mpi_comm_spec,my_barrier, calculate_test_sum,&
       &mpi_comm_v, initialize_comm_sim
  USE discretization
  USE coordinates
  USE par_other
  USE par_in
  USE geometry
  USE BoundaryDescriptionModule
  use MatrixModule
  USE BandedMatrixModule
  use VectorModule
  USE mpi
  USE DerivativeMatrixModule
  use Grid1DModule	
  USE file_io
  USE gyro_average_df_mod
  USE gyro_average_dd_mod
  USE test_pmm
  USE tga_io
  USE phys_ini

  IMPLICIT NONE

  INTEGER :: iProc

  PERFINIT
  PERFON('all')

  ! For now, all of the typical GENE inputs are just hardcoded into the
  ! test_gyroaverage_io module.
 
  !call tga_init
  call read_tga_parameters(par_in_dir)
  diagdir='./'
  parscheme='c4th'

  !Initializations   

  nblocks=1

  call mpi_init(ierr)
  CALL initialize_comm_sim(MPI_COMM_WORLD)

  call initialize_matrix_module(comm_cart)
  CALL initialize_BandedMatrix_module
  call initialize_Vector_module

  CALL mpi_comm_size(comm_cart,n_procs, ierr)
  CALL mpi_comm_rank(comm_cart,rank, ierr)
  write(*,"(A,I4,A)") "Using ",n_procs," processors."
  
  call initialize_current_parall

  ! Initialize FFTs for "wrapper" routines in gyro_average_dd
  CALL initialize_fourier_x_1d

  call initialize_gyro_average_df
  !call initialize_gyro_average_dd

  WRITE(*,*) "Initialized gyro matrix"

  CALL print_gyro_matrix_df

  PERFOFF

  call mpi_finalize(ierr)
  
END PROGRAM test_gyroaverage_df 

    
