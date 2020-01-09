MODULE all_rhs_terms
  USE nonlinear_term_mod,only: nonlinear_term_t
  IMPLICIT NONE

  PRIVATE

  CLASS(nonlinear_term_t),POINTER :: this_nonlinear_term
  !DEC$ ATTRIBUTES ALIGN:64 :: this_nonlinear_term
  
  PUBLIC :: this_nonlinear_term
END MODULE all_rhs_terms
