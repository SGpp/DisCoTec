;#############################################################################
;#    this library contains almost all internal functions for extended diags     #
;#    Note: data input from files with differen names - will be done in diag_loop
;#    Note: Procedure 'shells' left here for now
;#############################################################################

;#############################################################################

PRO create_extended_struct

  COMMON global_vars

END

;##########################################################################

PRO destroy_extended_struct

  COMMON global_vars

END

;##########################################################################

FUNCTION get_extended_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  first=0.0
  last=0.0
  get_n_steps = 1

END

;##########################################################################

PRO jump_to_extended_step
; jumps to the desired step of nlt file

  COMMON global_vars

END

;##########################################################################

PRO read_extended_step 

  COMMON global_vars

end

;##########################################################################

PRO extended_loop, diag0

  COMMON global_vars

  ;print,'Here we are in ext_loop!'

  call_diags, diag0, 'loop'

end

;##########################################################################
