# subspace writer

to get distributed sparse grid subspace size files for a scenario

## usage

call like a DisCoTec example executable, i.e.

```bash
export OMP_NUM_THREADS=1
mpiexec -np $nprocs ./subspace_writer $ctparam
```

with `$nprocs` the number of processes in one process group and `$ctparam` the
combination scheme (that specifies the ctscheme json file to read).

If you want the sizes of a conjoint grid, you also need to specify the other
scheme as

```ini
[subspace]
otherctscheme = ....json
```

in $ctparam.
