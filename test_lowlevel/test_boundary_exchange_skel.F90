program test_boundary_exchange
	use boundary_exchange
	implicit none

	logical :: passed_all_tests=.true.

	!-----------------------------------------
	! Testing interface exchange_x
	!-----------------------------------------
	if (test_exchange_x_1D()) then
		write(*,"(A40,I3,A10)") "test_exchange_x_1D on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_exchange_x_1D on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_exchange_x_2D()) then
		write(*,"(A40,I3,A10)") "test_exchange_x_2D on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_exchange_x_2D on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_exchange_x_3D()) then
		write(*,"(A40,I3,A10)") "test_exchange_x_3D on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_exchange_x_3D on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_exchange_x_4D()) then
		write(*,"(A40,I3,A10)") "test_exchange_x_4D on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_exchange_x_4D on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_exchange_x_5D()) then
		write(*,"(A40,I3,A10)") "test_exchange_x_5D on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_exchange_x_5D on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	if (test_exchange_x_6D()) then
		write(*,"(A40,I3,A10)") "test_exchange_x_6D on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_exchange_x_6D on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface initialize_type
	!-----------------------------------------
	if (test_mp_initialize_type()) then
		write(*,"(A40,I3,A10)") "test_mp_initialize_type on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_initialize_type on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface finalize_type
	!-----------------------------------------
	if (test_mp_finalize_type()) then
		write(*,"(A40,I3,A10)") "test_mp_finalize_type on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_finalize_type on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	!-----------------------------------------
	! Testing interface set_mpi_type
	!-----------------------------------------
	if (test_mp_set_mpi_type()) then
		write(*,"(A40,I3,A10)") "test_mp_set_mpi_type on process ",mype," passed"
	else
		write(*,"(A40,I3,A10)") "test_mp_set_mpi_type on process ",mype,"  FAILED"
		passed_all_tests=.false.
	end if

	if (passed_all_tests) then
		print*,'All tests on process ',mype,' passed.'
	else
		print*,'Some tests on process ',mype,' FAILED.'
	end if
contains

function test_exchange_x_1D() result(passed)
	logical :: passed

	passed=.true.
end function test_exchange_x_1D

function test_exchange_x_2D() result(passed)
	logical :: passed

	passed=.true.
end function test_exchange_x_2D

function test_exchange_x_3D() result(passed)
	logical :: passed

	passed=.true.
end function test_exchange_x_3D

function test_exchange_x_4D() result(passed)
	logical :: passed

	passed=.true.
end function test_exchange_x_4D

function test_exchange_x_5D() result(passed)
	logical :: passed

	passed=.true.
end function test_exchange_x_5D

function test_exchange_x_6D() result(passed)
	logical :: passed

	passed=.true.
end function test_exchange_x_6D

function test_mp_initialize_type() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_initialize_type

function test_mp_finalize_type() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_finalize_type

function test_mp_set_mpi_type() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_set_mpi_type
end program
