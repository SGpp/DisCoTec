MODULE Matrix_Vector_module
  use MatrixModule
  USE VectorModule

  INTERFACE dot_multiply
     module procedure MatVec_dot_multiply, MatVec_dot_multiply_real
  END INTERFACE

CONTAINS
  SUBROUTINE MatVec_dot_multiply(mat,vec, res)
    TYPE(Matrix),intent(IN) :: mat
    TYPE(Vector),intent(IN) :: vec
    TYPE(Vector),intent(INOUT) :: res

    !PRINT*,vec%NCols,vec%NRows,res%NCols,res%NRows
    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       CALL dot_multiply(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"Matrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF
  END SUBROUTINE MatVec_dot_multiply

  SUBROUTINE MatVec_dot_multiply_real(mat,vec, res)
    TYPE(MatrixReal),intent(IN) :: mat
    TYPE(VectorReal),intent(IN) :: vec
    TYPE(VectorReal),intent(INOUT) :: res

    !PRINT*,vec%NCols,vec%NRows,res%NCols,res%NRows
    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       CALL dot_multiply(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"Matrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF
  END SUBROUTINE MatVec_dot_multiply_real
END MODULE Matrix_Vector_module
