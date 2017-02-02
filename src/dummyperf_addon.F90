SUBROUTINE perf_get(name,inc_time, inc_MFlops)
  CHARACTER(len=*),intent(IN) :: name
  REAL(8),intent(OUT) :: inc_time, inc_MFlops

  inc_time = 0.0
  inc_MFlops = 0.0
END SUBROUTINE perf_get

SUBROUTINE perf_reset(name)
  CHARACTER(len=*), intent(IN) :: name

END SUBROUTINE perf_reset

SUBROUTINE perf_context_start(name)
  CHARACTER(len=*), INTENT(IN) :: name
END SUBROUTINE perf_context_start

SUBROUTINE perf_context_end()
END SUBROUTINE perf_context_end
