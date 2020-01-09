!QUIET = 1

PRINT, '--- starting GENE diagnostics ---'
PRINT, ''

IF !QUIET NE 1 THEN !EXCEPT = 2 ELSE !EXCEPT = 1

.r internal/global_vars.pro
.r internal/tools.pro
.r internal/diag_struct.pro
.r internal/data_struct.pro
.r internal/data_loop.pro
.r internal/mom_lib.pro
.r internal/nrg_lib.pro
.r internal/vsp_lib.pro
.r internal/scan_lib.pro
.r internal/srcmom_lib.pro
.r internal/chpt_lib.pro
.r internal/energy_lib.pro
.r internal/out_device.pro
.r internal/gui.pro
.r internal/iterdb_lib.pro
.r internal/synth_diag_D3D.pro
.r internal/nlt_lib.pro
.r internal/extended_lib.pro

@internal/inc_diags.pro

CLOSE, /ALL
global_vars
gui_main

PRINT, ''
PRINT, '--- exiting GENE diagnostics ---'
