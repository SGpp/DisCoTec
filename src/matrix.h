  TYPE,public :: Matrix
     INTEGER :: NCols
     INTEGER :: NRows
     TYPE(StoreFullMatrixObject) :: DATA
  END TYPE matrix

  TYPE,public :: MatrixReal
     INTEGER :: NCols
     INTEGER :: NRows
     TYPE(StoreFullMatrixObjectReal) :: DATA
  END TYPE matrixReal

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_matrix, mp_initialize_vector
     MODULE PROCEDURE mp_initialize_matrix_real, mp_initialize_vector_real
  END INTERFACE

  INTERFACE initialize_matrix_module
     MODULE PROCEDURE mp_initialize_matrix_module
  END INTERFACE

  INTERFACE allocate
    MODULE PROCEDURE mp_allocate, mp_allocate_real
  END INTERFACE

  INTERFACE attach
    MODULE PROCEDURE mp_attach_vector, mp_attach_matrix, mp_attach_vector_real, mp_attach_matrix_real
  END INTERFACE

  INTERFACE finalize
    MODULE PROCEDURE mp_finalize_matrix, mp_finalize_matrix_real
  END INTERFACE

  INTERFACE finalize_matrix_module
     MODULE PROCEDURE mp_finalize_module
  END INTERFACE

  INTERFACE set_value
    MODULE PROCEDURE mp_set_complex_value, mp_set_real_value, mp_set_real_value_real
  END INTERFACE

  INTERFACE add_value
    MODULE PROCEDURE mp_add_complex_value, mp_add_real_value, mp_add_real_value_real
  END INTERFACE

  interface add_to_row
    module procedure mat_add_to_row
  end interface

  INTERFACE commit_values
    module procedure mp_commit_values, mp_commit_values_real
  END INTERFACE

  INTERFACE set_zero
    MODULE PROCEDURE mp_set_zero, mp_set_zero_real
  END INTERFACE

  INTERFACE get_global_matrix_locally
    module procedure mp_get_global_matrix_locally, mp_get_global_matrix_locally_real
  END INTERFACE

  INTERFACE get_local_abs_square_sum
    module procedure mp_get_local_abs_square_sum, mp_get_local_abs_square_sum_real
  END INTERFACE

  INTERFACE mat_get_value
    MODULE PROCEDURE mp_get_value, mp_get_value_real
  END INTERFACE

  INTERFACE dot_multiply
    MODULE PROCEDURE mp_dot_multiply_MatMat, mp_dot_multiply_MatMat_real
  END INTERFACE

  INTERFACE show
    MODULE PROCEDURE mp_show_on_screen, mp_show_in_file, mp_show_on_screen_real, mp_show_in_file_real
  END INTERFACE

  INTERFACE add_matrix
    MODULE PROCEDURE mp_add_to_matrix, mp_add_matrix, mp_add_to_matrix_real, mp_add_matrix_real
  END INTERFACE

  interface subtract_matrix
  !INTERFACE OPERATOR(-)
    MODULE PROCEDURE mp_subtract_matrix, mp_subtract_from_matrix, mp_subtract_matrix_real, mp_subtract_from_matrix_real
  END INTERFACE

  interface multiply_matrix_with_scalar
  !INTERFACE OPERATOR(*)
    module procedure mp_multiply_matrix_with_real, mp_multiply_matrix_with_real_real
    module procedure mp_scale_matrix_by_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
    module procedure mp_assign_matrix, mp_assign_matrix_real
  END INTERFACE

  INTERFACE invert
    module procedure mp_invert_matrix, mp_invert_matrix_real
  END INTERFACE

  INTERFACE output_data
    module procedure mp_output_data_matrix, mp_output_data_matrix_real
  END INTERFACE

  INTERFACE isInitialized
    module procedure mp_isInitialized, mp_isInitialized_real
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE my_matrix_sum_generic, my_matrix_sum_2D, my_matrix_sum_3D, my_matrix_sum_0D
    MODULE PROCEDURE my_matrix_sum_generic_real, my_matrix_sum_2D_real, my_matrix_sum_3D_real, my_matrix_sum_0D_real
  END INTERFACE

  INTERFACE LU_factor
    MODULE PROCEDURE mp_LU_factor, mp_LU_factor_real
  END INTERFACE

  INTERFACE LU_solve
    MODULE PROCEDURE mp_LU_solve, mp_LU_solve_real
  END INTERFACE
