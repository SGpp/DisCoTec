diagfiles = FILE_SEARCH('prog/*.pro',/TEST_REGULAR)
cdiagfiles = FILE_SEARCH('custom/*.pro',/TEST_REGULAR)

FILE_DELETE, 'inc_files.pro', /ALLOW_NONEXISTENT

OPENW, inc_file, 'inc_files.pro', /GET_LUN

FOR i = 0, N_ELEMENTS(diagfiles) - 1 DO $
  PRINTF, inc_file, '.r ' + diagfiles[i]
FOR i = 0, N_ELEMENTS(cdiagfiles) - 1 DO $
  PRINTF, inc_file, '.r ' + cdiagfiles[i]

CLOSE, inc_file
FREE_LUN, inc_file

@inc_files.pro

FILE_DELETE, 'inc_files.pro', /ALLOW_NONEXISTENT
