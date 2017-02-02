PROGRAM test_ListObject
  USE ListModule
  IMPLICIT NONE

  type(List) :: mylist

  type(Node) :: mynode

  mynode%coordinates = (/1,2/)
  mynode%value = 1.0

  call initialize(mylist)
  call push(mylist, mynode)
  call push(mylist, mynode)
  call push(mylist, mynode)
  call push(mylist, mynode)

  CALL finalize(mylist)
END PROGRAM test_ListObject
