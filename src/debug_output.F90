#include "intrinsic_sizes.h"
MODULE debug_output
  IMPLICIT NONE
  
  INTEGER,PARAMETER :: MAXDIMS=6
  INTERFACE output_data
     MODULE PROCEDURE output_data_1D,output_data_2D,output_data_3D,output_data_4D,&
          &output_data_5D,output_data_6D
  END INTERFACE

CONTAINS
  SUBROUTINE output_data_1D(fn_main,mat)
    CHARACTER(*) :: fn_main
    COMPLEX :: mat(:)

    INTEGER,DIMENSION(1:MAXDIMS) :: dims
    
    dims(1) = SIZE(mat,1)
    CALL output_data_generic(fn_main,mat,1,dims)
  END SUBROUTINE output_data_1D

  SUBROUTINE output_data_2D(fn_main,mat)
    CHARACTER(*) :: fn_main
    COMPLEX :: mat(:,:)

    INTEGER,DIMENSION(1:MAXDIMS) :: dims
    
    dims(1) = SIZE(mat,1)
    dims(2) = SIZE(mat,2)
    CALL output_data_generic(fn_main,mat,2,dims)
  END SUBROUTINE output_data_2D

  SUBROUTINE output_data_3D(fn_main,mat)
    CHARACTER(*) :: fn_main
    COMPLEX :: mat(:,:,:)

    INTEGER,DIMENSION(1:MAXDIMS) :: dims
    
    dims(1) = SIZE(mat,1)
    dims(2) = SIZE(mat,2)
    dims(3) = SIZE(mat,3)
    CALL output_data_generic(fn_main,mat,3,dims)
  END SUBROUTINE output_data_3D

  SUBROUTINE output_data_4D(fn_main,mat)
    CHARACTER(*) :: fn_main
    COMPLEX :: mat(:,:,:,:)

    INTEGER,DIMENSION(1:MAXDIMS) :: dims
    
    dims(1) = SIZE(mat,1)
    dims(2) = SIZE(mat,2)
    dims(3) = SIZE(mat,3)
    dims(4) = SIZE(mat,4)
    CALL output_data_generic(fn_main,mat,4,dims)
  END SUBROUTINE output_data_4D

  SUBROUTINE output_data_5D(fn_main,mat)
    CHARACTER(*) :: fn_main
    COMPLEX :: mat(:,:,:,:,:)

    INTEGER,DIMENSION(1:MAXDIMS) :: dims
    
    dims(1) = SIZE(mat,1)
    dims(2) = SIZE(mat,2)
    dims(3) = SIZE(mat,3)
    dims(4) = SIZE(mat,4)
    dims(5) = SIZE(mat,5)
    CALL output_data_generic(fn_main,mat,5,dims)
  END SUBROUTINE output_data_5D

  SUBROUTINE output_data_6D(fn_main,mat)
    CHARACTER(*) :: fn_main
    COMPLEX :: mat(:,:,:,:,:,:)

    INTEGER,DIMENSION(1:MAXDIMS) :: dims
    
    dims(1) = SIZE(mat,1)
    dims(2) = SIZE(mat,2)
    dims(3) = SIZE(mat,3)
    dims(4) = SIZE(mat,4)
    dims(5) = SIZE(mat,5)
    dims(6) = SIZE(mat,6)
    CALL output_data_generic(fn_main,mat,6,dims)
  END SUBROUTINE output_data_6D
  
SUBROUTINE output_data_generic(fn_main,mat,ndim,dims)
  CHARACTER(*) :: fn_main
  COMPLEX :: mat(*)
  INTEGER :: ndim
  INTEGER,DIMENSION(1:MAXDIMS) :: dims
  
  INTEGER :: thisunit,idim,totaldim
  LOGICAL :: op
  CHARACTER(len=FILENAME_MAX) :: filename

  thisunit=30
  DO 
     INQUIRE(thisunit,opened=op)
     IF (op) THEN
        thisunit = thisunit+1
     ELSE 
        EXIT
     END IF
  END DO

  WRITE(filename,"(3A)") "./",TRIM(fn_main),".dat"
  PRINT*,"Writing to file ",trim(filename)
  OPEN(thisunit,file=TRIM(filename))
  totaldim = 1.0
  DO idim=1,ndim
     WRITE(thisunit,"(I3,1X)",advance='no') dims(idim)
     totaldim = totaldim * dims(idim)
  END DO
  WRITE(thisunit,"(A)") ""
  WRITE(thisunit,"(4ES20.10)") mat(1:totaldim)
  CLOSE(thisunit)
END SUBROUTINE output_data_generic

END MODULE debug_output
