program test_MatrixModule
	use MatrixModule
	implicit none

	logical :: passed_all_tests=.true.
	if (test_initialize()) then
		write(*,"(A40,A10)") "test_initialize","passed"
	else
		write(*,"(A40,A10)") "test_initialize","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_initialize_matrix_module()) then
		write(*,"(A40,A10)") "test_initialize_matrix_module","passed"
	else
		write(*,"(A40,A10)") "test_initialize_matrix_module","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_allocate()) then
		write(*,"(A40,A10)") "test_allocate","passed"
	else
		write(*,"(A40,A10)") "test_allocate","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_attach()) then
		write(*,"(A40,A10)") "test_attach","passed"
	else
		write(*,"(A40,A10)") "test_attach","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_finalize()) then
		write(*,"(A40,A10)") "test_finalize","passed"
	else
		write(*,"(A40,A10)") "test_finalize","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_set_value()) then
		write(*,"(A40,A10)") "test_set_value","passed"
	else
		write(*,"(A40,A10)") "test_set_value","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_add_value()) then
		write(*,"(A40,A10)") "test_add_value","passed"
	else
		write(*,"(A40,A10)") "test_add_value","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_commit_values()) then
		write(*,"(A40,A10)") "test_commit_values","passed"
	else
		write(*,"(A40,A10)") "test_commit_values","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_set_zero()) then
		write(*,"(A40,A10)") "test_set_zero","passed"
	else
		write(*,"(A40,A10)") "test_set_zero","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_get_global_matrix_locally()) then
		write(*,"(A40,A10)") "test_get_global_matrix_locally","passed"
	else
		write(*,"(A40,A10)") "test_get_global_matrix_locally","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_get_local_abs_square_sum()) then
		write(*,"(A40,A10)") "test_get_local_abs_square_sum","passed"
	else
		write(*,"(A40,A10)") "test_get_local_abs_square_sum","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_mat_get_value()) then
		write(*,"(A40,A10)") "test_mat_get_value","passed"
	else
		write(*,"(A40,A10)") "test_mat_get_value","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_dot_multiply()) then
		write(*,"(A40,A10)") "test_dot_multiply","passed"
	else
		write(*,"(A40,A10)") "test_dot_multiply","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_show()) then
		write(*,"(A40,A10)") "test_show","passed"
	else
		write(*,"(A40,A10)") "test_show","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_add_matrix()) then
		write(*,"(A40,A10)") "test_add_matrix","passed"
	else
		write(*,"(A40,A10)") "test_add_matrix","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_subtract_matrix()) then
		write(*,"(A40,A10)") "test_subtract_matrix","passed"
	else
		write(*,"(A40,A10)") "test_subtract_matrix","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_multiply_matrix_with_scalar()) then
		write(*,"(A40,A10)") "test_multiply_matrix_with_scalar","passed"
	else
		write(*,"(A40,A10)") "test_multiply_matrix_with_scalar","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_ASSIGNMENT_equal()) then
		write(*,"(A40,A10)") "test_ASSIGNMENT_equal","passed"
	else
		write(*,"(A40,A10)") "test_ASSIGNMENT_equal","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_invert()) then
		write(*,"(A40,A10)") "test_invert","passed"
	else
		write(*,"(A40,A10)") "test_invert","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_output_data()) then
		write(*,"(A40,A10)") "test_output_data","passed"
	else
		write(*,"(A40,A10)") "test_output_data","  FAILED"
		passed_all_tests=.false.
	end if

	if (test_my_sum()) then
		write(*,"(A40,A10)") "test_my_sum","passed"
	else
		write(*,"(A40,A10)") "test_my_sum","  FAILED"
		passed_all_tests=.false.
	end if

	if (passed_all_tests) then
		print*,'All tests passed.'
	else
		print*,'Some tests FAILED.'
	end if
contains

function test_initialize() result(passed)
	logical :: passed

	passed=.true.
end function test_initialize

function test_initialize_matrix_module() result(passed)
	logical :: passed

	passed=.true.
end function test_initialize_matrix_module

function test_allocate() result(passed)
	logical :: passed

	passed=.true.
end function test_allocate

function test_attach() result(passed)
	logical :: passed

	passed=.true.
end function test_attach

function test_finalize() result(passed)
	logical :: passed

	passed=.true.
end function test_finalize

function test_set_value() result(passed)
	logical :: passed

	passed=.true.
end function test_set_value

function test_add_value() result(passed)
	logical :: passed

	passed=.true.
end function test_add_value

function test_commit_values() result(passed)
	logical :: passed

	passed=.true.
end function test_commit_values

function test_set_zero() result(passed)
	logical :: passed

	passed=.true.
end function test_set_zero

function test_get_global_matrix_locally() result(passed)
	logical :: passed

	passed=.true.
end function test_get_global_matrix_locally

function test_get_local_abs_square_sum() result(passed)
	logical :: passed

	passed=.true.
end function test_get_local_abs_square_sum

function test_mat_get_value() result(passed)
	logical :: passed

	passed=.true.
end function test_mat_get_value

function test_dot_multiply() result(passed)
	logical :: passed

	passed=.true.
end function test_dot_multiply

function test_show() result(passed)
	logical :: passed

	passed=.true.
end function test_show

function test_add_matrix() result(passed)
	logical :: passed

	passed=.true.
end function test_add_matrix

function test_subtract_matrix() result(passed)
	logical :: passed

	passed=.true.
end function test_subtract_matrix

function test_multiply_matrix_with_scalar() result(passed)
	logical :: passed

	passed=.true.
end function test_multiply_matrix_with_scalar

function test_ASSIGNMENT_equal() result(passed)
	logical :: passed

	passed=.true.
end function test_ASSIGNMENT_equal

function test_invert() result(passed)
	logical :: passed

	passed=.true.
end function test_invert

function test_output_data() result(passed)
	logical :: passed

	passed=.true.
end function test_output_data

function test_my_sum() result(passed)
	logical :: passed

	passed=.true.
end function test_my_sum
end program
