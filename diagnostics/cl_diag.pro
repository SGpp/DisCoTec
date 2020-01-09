!QUIET = 1

PRINT, '--- starting GENE diagnostics, command line version ---'

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
.r internal/chpt_lib.pro
.r internal/energy_lib.pro
.r internal/out_device.pro
.r internal/gui.pro
.r cl_par.pro

@internal/inc_diags.pro

CLOSE, /ALL
global_vars

read_cfg
read_cl_params
create_diag_struct
cl_setdiags, momtab=momtab

COMMON global_vars

ps_format = WHERE(gui.out.out_format_str EQ 'postscript')
IF ps_format[0] NE -1 THEN gui.out.out_format[ps_format] = 1

IF par.n_spec GT 0 THEN data_loop, tab=momtab

IF PTR_VALID(series.geom) THEN PTR_FREE, series.geom
free_diag_memory

PTR_FREE, nrg, vsp, scanlog, $
  gui.out.spec_select, gui.tab.arr,$
  par.lxarr, par.ky, par.kx, par.z, $
  par.istep_field, par.istep_mom, par.istep_nrg, par.istep_energy, $
  series.filename_fmt, series.geom,$
  series.ky, series.lxarr, series.kx, spec[*].prof

IF !QUIET NE 1 THEN HELP, /HEAP_VARIABLES
HEAP_GC
IF !QUIET NE 1 THEN PRINT, 'Maximum memory use since last mem_usage call: ' + mem_usage(/highw)
IF !QUIET NE 1 THEN PRINT, 'Memory still in use: ' + mem_usage()

PRINT, '--- exiting GENE diagnostics, command line version ---'
