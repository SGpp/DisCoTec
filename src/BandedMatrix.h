!> Just knows the mathematical extend of the matrix, that is the 
!! number of rows and columns. It also has a data object embedded
!! which contains the real data in form of a StoreBandedMatrixObject.
  TYPE BandedMatrix
     INTEGER :: NCols
     INTEGER :: NRows
     TYPE(StoreBandedMatrixObject) :: DATA
  END TYPE BandedMatrix

  TYPE BandedMatrixReal
    INTEGER :: NCols
    INTEGER :: NRows
    TYPE(StoreBandedMatrixObjectReal) :: DATA
  END TYPE BandedMatrixReal

  INTERFACE initialize
    MODULE PROCEDURE mp_initialize_matrix, mp_initialize_vector, mp_initialize_matrix_real, mp_initialize_vector_real
  END INTERFACE

  INTERFACE initialize_BandedMatrix_module
    MODULE PROCEDURE mp_initialize_BandedMatrix_module
  END INTERFACE

  INTERFACE allocate
    MODULE PROCEDURE mp_allocate, mp_allocate_real
  END INTERFACE

  INTERFACE finalize
    MODULE PROCEDURE mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE

  INTERFACE finalize_BandedMatrix_module
     MODULE PROCEDURE mp_finalize_BandedMatrix_module
  END INTERFACE

  INTERFACE set_value
    MODULE PROCEDURE mp_set_complex_value, mp_set_real_value, mp_set_real_value_real
  END INTERFACE

  INTERFACE add_value
    MODULE PROCEDURE mp_add_complex_value, mp_add_real_value, mp_add_real_value_real
  END INTERFACE

  INTERFACE commit_values
    module procedure mp_commit_values, mp_commit_values_real
  END INTERFACE

  INTERFACE set_zero
    MODULE PROCEDURE mp_set_zero, mp_set_zero_real
  END INTERFACE

  INTERFACE get_global_matrix_locally
    module procedure mp_get_global_matrix_locally, mp_get_global_matrix_locally_real
  END INTERFACE

  INTERFACE convert_Banded_to_Full
    module procedure mp_convert_Banded_to_Full, mp_convert_Banded_to_Full_real
  END INTERFACE

  INTERFACE convert_Full_to_Banded
    module procedure mp_convert_Full_to_Banded, mp_convert_Full_to_Banded_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
    module procedure mp_get_local_abs_square_sum, mp_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE mat_get_value
    MODULE PROCEDURE mp_get_value, mp_get_value_real
  END INTERFACE

  INTERFACE mat_get_row_pointer
     MODULE PROCEDURE bm_get_row_pointer
  END INTERFACE

  interface autotune
     module procedure bm_autotune
  end interface

  INTERFACE dot_multiply
    MODULE PROCEDURE mp_dot_multiply, bm_matmat_dot_multiply, mp_dot_multiply_real
    module procedure bm_dot_multiply_Banded_with_Full, bm_matmat_dot_multiply_real
    module procedure bm_dot_multiply_Banded_with_Full_real
  END INTERFACE

  INTERFACE square_matrix
     module procedure bm_square_matrix
  END INTERFACE

  INTERFACE show
    MODULE PROCEDURE mp_show_on_screen, mp_show_in_file, mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE add_matrix
    MODULE PROCEDURE mp_add_to_matrix, mp_add_matrix!, mp_add_Matrix_to_BandedMatrix
    MODULE PROCEDURE mp_add_to_matrix_real
  END INTERFACE

  interface subtract_matrix
  !INTERFACE OPERATOR(-)
     MODULE PROCEDURE mp_subtract_matrix, mp_subtract_from_matrix
    MODULE PROCEDURE mp_subtract_matrix_real, mp_subtract_from_matrix_real
  END INTERFACE

  interface multiply_matrix_with_scalar
  !INTERFACE OPERATOR(*)
    module procedure mp_multiply_matrix_with_real, mp_multiply_matrix_with_real_real
    module procedure mp_scale_matrix_by_real, mp_scale_matrix_by_real_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
    module procedure mp_assign_matrix, mp_assign_matrix_real
  END INTERFACE

  INTERFACE output_data
    module procedure mp_output_data_matrix, mp_output_data_matrix_real
  END INTERFACE

  INTERFACE isInitialized
    module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

  INTERFACE print_storage_details
     MODULE PROCEDURE bm_print_storage_details
  END INTERFACE

  INTERFACE get_number_of_bands
     module procedure bm_get_number_of_bands
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE bm_sum_generic, bm_sum_2D, bm_sum_3D, bm_sum_0D
  END INTERFACE

  INTERFACE row_axpy
     module procedure bm_row_axpy
  END INTERFACE

  interface transpose_and_conjugate
     module procedure bm_transpose_and_conjugate
  end interface

  interface transpose_storage
     module procedure bm_transpose_storage
  end interface

  INTERFACE LU_factor
     module procedure bm_LU_factor
  END INTERFACE

  INTERFACE LU_factor_ok
     module procedure bm_LU_factor_ok
  END INTERFACE

  INTERFACE LU_solve
     module procedure bm_LU_solve
  END INTERFACE

#ifdef ALL_ROUTINES
  INTERFACE invert
    module procedure mp_invert_matrix, mp_invert_matrix_real
  END INTERFACE


#endif
