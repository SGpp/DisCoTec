What is this project about?
---------------------------
[![Build Status](https://simsgs.informatik.uni-stuttgart.de/jenkins/buildStatus/icon?job=DisCoTec-main-spack)](https://simsgs.informatik.uni-stuttgart.de/jenkins/view/DisCoTec/job/DisCoTec-main-spack/)

This project contains __DisCoTec__, a code for running the distributed sparse grid combination technique with MPI parallelization. While it originates from the excellent [SGpp project](https://github.com/SGpp/SGpp), all the parallelization makes it a very different code, such that it is its own project now.

Contributing
---------
* Please develop new features in a new branch and then create a pull request to
  notify other users about the changes. So everyone can check whether there are
  side-effects to their work.
* New features will only be merged to the main branch if they are sufficiently tested.

Requirements
--------------
cmake >= (3.24.2),
libmpich-dev (>= 3.2-6), or other MPI library

Additional (optional) dependencies (can be obtained with cmake):
- OpenMP
- libboost-all-dev (>= 1.60)
- HDF5
- HighFive 

Installation instructions: 
--------------------------
#### Complete build:
```bash
git clone https://github.com/SGpp/DisCoTec.git DisCoTec
cmake -S DisCoTec -B DisCoTec/build
cmake --build DisCoTec/build -j
```

#### Advanced build:

All examples and the tools can be build on their own. To build a specific example, run the specific cmake command in the example folder. For example, run the following commands
```bash
git clone https://github.com/SGpp/DisCoTec.git DisCoTec
cd DisCoTec/examples/combi_example
cmake -S . -B build
cmake --build build -j
```
to build the combi_example.


For building only the DisCoTec library, run cmake with the `src` folder as source folder.

#### Optional CMake Options
- `DISCOTEC_TEST=**ON**|OFF` - Build tests if you build the complete project.
- `DISCOTEC_BUILD_MISSING_DEPS=**ON**|OFF`- First order dependencies that are not found are built automatically (glpk is always built).
- `DISCOTEC_TIMING=**ON**|OFF` - Enables internal timing
- `DISCOTEC_USE_HDF5=**ON**|OFF`
- `DISCOTEC_USE_HIGHFIVE=**ON**|OFF` - Enables HDF5 support via HighFive. If `DISCOTEC_USE_HIGHFIVE=ON`, `DISCOTEC_USE_HDF5` has also to be `ON`.
- `DISCOTEC_UNIFORMDECOMPOSITION=**ON **|OFF` - Enables the uniform decomposition of the grid.
- `DISCOTEC_GENE=ON|**OFF**` - Currently GEne is not supported with CMake!
- `DISCOTEC_OPENMP=ON|**OFF**` - Enables OpenMP support.
- `DISCOTEC_ENABLEFT=ON|**OFF**` - Enables the use of the FT library.
- `DISCOTEC_USE_LTO=**ON**|OFF` - Enables link time optimization if the compiler supports it.
- `DISCOTEC_OMITREADYSIGNAL=ON|**OFF**` - Omit the ready signal in the MPI communication. This can be used to reduce the communication overhead.
- `DISCOTEC_USENONBLOCKINGMPICOLLECTIVE=ON|**OFF**` - TODO: Add description


To run the compiled tests, go to folder `tests` and run
```bash
mpiexec -np 9 ./test_distributedcombigrid_boost
```
where you can use all the parameters of the boost test suite.
Or you can run the tests with `ctest` in the build folder.


GENE  submodules as dependencies for GENE examples
----------------
_Warning: The CMake Integration is currently not adappted to use GENE and SeLaLib!_

There are gene versions as submodules: a linear one in the gene_mgr folder, and 
a nonlinear one in gene-non-linear. To get them, you need access to their 
respective repos. Then go into the folder and

``` bash
git submodule init
git submodule update
```
or use the `--recursive` flag when cloning the repo.
and then you just hope you can `make` it by adapting the local library flags.
and then you just hope you can `make` it by adapting the local library flags ;)

Third-Level-Manager
-------------------
The third-level-manager handles the synchronization and data transfer between
the systems. The following steps will guide you from compilation to execution:

1. To compile first navigate to third_level_manager.
   The folder contains a Makefile.template which you can adjust to make it
   compatible with your maschine.

2. The manager takes a .ini file as an input parameter. The same folder as in 1.
   holds the example file `example.ini` where the number of systems and the port
   on which the manager listens can be adjusted. Currently, only 2 systems are
   supported.

3. Use the run script to execute the manager. You can pass your own parameter
   file to the script whereas by default it reads the example.ini.

On Hazel Hen
------------

Warning: outdated not testet with the new CMake build scripts!!


### Execution

Let's assume we want to run the example under
combi/distributedcombigrid/examples/distributed_third_level/ and distribute our
combischeme to **HLRS** and **helium**, while the third-level-manager is running
on helium at port 9999. The following steps are necessary:

1. Since data connections to the HLRS are not possible without using ssh
   tunnels, we set them up in advance. Run `ssh -R  9998:localhost:9999
   username@eslogin002.hww.de` from helium to log in to HLRS while creating an
   ssh tunnel.

2. Adjust the parameter files on HLRS and helium to fit the simulation. Use
   hostname=eslogin002 on HLRS and hostname=localhost on helium. Set
   dataport=9999 on both systems.

3. Run the third-level-manger on helium.

4. Connect to eslogin002 in a different terminal and run the forwarding
   script `forward.sh 9999 9998 pipe1`. This will forward the port 9998 to 9999
   on eslogin002. (We only need the local forwarding because the configuration
   of the ssh server on the HLRS does not allow us to access the ssh tunnel
   from a different host than eslogin002. Since our application runs on the
   compute nodes (for now) this detour is necessary.)

5. Start the simulation. The example directory holds a separate run file
   `run.sh` which needs to be modified to fit HLRS and helium. Also set the
   corresponding boost library location.
