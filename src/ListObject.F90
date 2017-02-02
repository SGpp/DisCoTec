#if 0
#define CONCAT(x,y) x ## y
#define MODULENAME(a,b) CONCAT(a,b)
#define TYPENAME(a) CONCAT(a,NodeData)
#define LISTNAME(a) CONCAT(a,List)
!This is an alternative (clever) way to determine all the precompiler 
!variables below. Unfortunately, it doesn't work on every compiler, 
!e.g., on pathscale, where the ## concatentation is not implemented yet
!(see bug #10229)
#endif

#define VALUE_NAME MatrixEntry
#define MODULENODENAME MatrixEntryNodeModule
#define MODULELISTNAME MatrixEntryListModule
#define TYPENAME MatrixEntryNodeData
#define LISTNAME MatrixEntryList
#define TYPENAMEREAL MatrixEntryNodeDataReal

MODULE MODULENODENAME
  TYPE TYPENAME
     INTEGER :: coord(2) !< row, column in the final gyromatrix
     COMPLEX :: value !< value of the entry
  END TYPE TYPENAME

  TYPE TYPENAMEREAL
       INTEGER :: coord(2)
       REAL :: value
  END TYPE TYPENAMEREAL

  INTERFACE output
     module procedure mp_output
  END INTERFACE

CONTAINS
  SUBROUTINE mp_output(this)
    type(TYPENAME) :: this

    WRITE(*,"(2(A,I4),A,2ES20.10)",advance='no') "(",this%coord(1),",",this%coord(2),")=",this%value
  END SUBROUTINE mp_output
END MODULE MODULENODENAME
#include "ListObjectTemplate.f90"

#undef VALUE_NAME
#undef MODULENODENAME
#undef MODULELISTNAME
#undef TYPENAME
#undef LISTNAME

#define VALUE_NAME Integer
#define MODULENODENAME IntegerNodeModule
#define MODULELISTNAME IntegerListModule
#define TYPENAME IntegerNodeData
#define LISTNAME IntegerList

MODULE MODULENODENAME
  TYPE TYPENAME
     INTEGER :: value
  END TYPE TYPENAME

  INTERFACE output
     module procedure mp_output
  END INTERFACE
CONTAINS
  SUBROUTINE mp_output(this)
    type(TYPENAME) :: this

    WRITE(*,"(I10)", advance='no') this%value
  END SUBROUTINE mp_output
END MODULE MODULENODENAME
#include "ListObjectTemplate.f90"

#undef VALUE_NAME
#undef MODULENODENAME
#undef MODULELISTNAME
#undef TYPENAME
#undef LISTNAME

#define VALUE_NAME MatrixEntry
#define MODULENODENAME MatrixEntryNodeModuleReal
#define MODULELISTNAME MatrixEntryListModuleReal
#define TYPENAME MatrixEntryNodeDataReal
#define LISTNAME MatrixEntryListReal

MODULE MODULENODENAME
  TYPE TYPENAME
     INTEGER :: coord(2) !< row, column in the final gyromatrix
     REAL :: value !< value of the entry
  END TYPE TYPENAME

  INTERFACE output
     module procedure mp_output
  END INTERFACE

CONTAINS
  SUBROUTINE mp_output(this)
    type(TYPENAME) :: this

    WRITE(*,"(2(A,I4),A,2ES20.10)",advance='no') "(",this%coord(1),",",this%coord(2),")=",this%value
  END SUBROUTINE mp_output
END MODULE MODULENODENAME
#include "ListObjectTemplate.f90"

