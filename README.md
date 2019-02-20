What is this project about?
---------------------------
This project contains the __distributed combigrid module__. 

Important:
----------
Please do use feature branches in development.

Guidelines
---------
*  Please develop new features in a new branch and then create a merge request 
to notify other users about the changes. So everyone can check whether there are 
side-effects to other branches.
* Of course, you are still very welcome to directly fix obvious bugs in the 
master branch.
* Before merging new features to the master branch, please make sure that they
are sufficiently commented. 
* Although the distributed combigrid module is independent of the other modules
in SG++, it will remain a part of this project. To ensure compability please
make sure that you only change files in the distributedcombigrid folder. 
* In the future the automated testing and code style checking on commits that is 
done in SG++ will also be done for this project.

Installation instructions: 
--------------------------
Compile with e.g.
```
scons -j4 SG_JAVA=0 COMPILE_BOOST_TESTS=1 SG_ALL=0 SG_DISTRIBUTEDCOMBIGRID=1 VERBOSE=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 COMPILER=mpich OPT=1
``` 
and for debugging, add
```
CPPFLAGS=-g
``` 

There are gene versions as submodules: a linear one in the gene_mgr folder, and 
a nonlinear one in gene-non-linear. To get them, you need access to their 
respective repos. Then go into the folder and

```
git submodule init
git submodule update
```
and then you just hope you can `make` it by adapting the local library flags ;)


On Hazel Hen
--------------
load modules: PrgEnv-gnu, scons, python 2.7

Unfortunately, the boost module on Hazel Hen does not work with our code any more.
Hence, you have to install boost yourself.
Do not forget to set boost paths in SConfigure, i.e. BOOST_INCLUDE_PATH, 
BOOST_LIBRARY_PATH

compile with
scons -j 16 SG_ALL=0 SG_DISTRIBUTEDCOMBIGRID=1 VERBOSE=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CXX=CC OPT=1

(the linking of the boost tests might fail. however this is not a problem, the
sg++ libraries should be there)



