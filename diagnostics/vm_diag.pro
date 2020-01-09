PRO vm_diag

  JOURNAL, 'vm_diag.log'

  !QUIET = 1

  PRINT, '--- starting GENE diagnostics, VM version ---'
  PRINT, ''

  IF !QUIET NE 1 THEN !EXCEPT = 2 ELSE !EXCEPT = 1

  CLOSE, /ALL
  global_vars
  COMMON global_vars
  vm_sav = 1
  gui_main

  PRINT, ''
  PRINT, '--- exiting GENE diagnostics, VM version ---'

  JOURNAL

END
