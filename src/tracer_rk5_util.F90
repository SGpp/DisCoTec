module tracer_rk5_util

  implicit none
  
  public:: nrerror, assert_eq
  interface assert_eq
     module procedure assert_eq2,assert_eq3,assert_eq4,assert_eqn
  end interface
  private
  !Subprograms
  
  
contains
  
  function assert_eq2(n1,n2,string)
    character(LEN=*), intent(IN) :: string
    integer, intent(IN) :: n1,n2
    integer :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2
  !BL
  function assert_eq3(n1,n2,n3,string)
    character(LEN=*), intent(IN) :: string
    integer, intent(IN) :: n1,n2,n3
    integer :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq3'
    end if
  end function assert_eq3
  !BL
  function assert_eq4(n1,n2,n3,n4,string)
    character(LEN=*), intent(IN) :: string
    integer, intent(IN) :: n1,n2,n3,n4
    integer :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq4'
    end if
  end function assert_eq4
  !BL
  function assert_eqn(nn,string)
    character(LEN=*), intent(IN) :: string
    integer, dimension(:), intent(IN) :: nn
    integer :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eqn'
    end if
  end function assert_eqn
  !BL
  subroutine nrerror(string)
    character(LEN=*), intent(IN) :: string
    write (*,*) 'nrerror: ',string
    stop 'program terminated by nrerror'
  end subroutine nrerror
  
  
end module tracer_rk5_util
