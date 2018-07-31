What is this project about?
---------------------------
This project contains the __distributed combigrid module__. It also contains the 
other SG++ modules from the release 2.0 version, because the build system 
expects some of them. I think this will simplify future merging.
Note, however, the only module which is to be further developed and maintained
in this project is the distributedcombigrid module. If you want to work on
the "classical" SG++ modules, you should get access to the corresponding project.

Important:
----------
Please do only push changes to the master branch which are of general nature.
Please put anything that is specific to a certain application, e.g., GENE or
DUNE, into an own branch. Furthermore, please create new branches to develop new
features, e.g. fault-tolerance functionality.

Guidelines
---------
* The master branch is supposed to be kept clean (I know that the current state
is far from "clean", but the goal should be to change this some day) Please 
develop new features in a new branch and then create a merge request to notify 
other users about the changes. So everyone can check whether there are 
side-effects to other branches.
* Of course, you are still very welcome to directly fix obvious bugs in the 
master branch.
* Before merging new features to the master branch, please make sure that they
are sufficiently commented. The current state of the master branch is not a good
example for good comment style, so please look at the other modules of SG++ to
get an impression how the comments should look like.
* Although the distributed combigrid module is independent of the other modules
in SG++, it will remain a part of this project. To ensure compability please
make sure that you only change files in the distributedcombigrid folder. 
* For the above reason, the other SG++ modules and also the build system will be
updated on a regular basis to the newest version. So, if you have to modify the
build system, please try to figure out a way which allows easy merging.
* In future the automated testing and code style checking on commits that is 
done in SG++ will also be done for this project.

Installation instructions: 
--------------------------
Compile with  
scons -j 4 SG_ALL=0 SG_DISTRIBUTEDCOMBIGRID=1 VERBOSE=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CXX=mpicxx.mpich OPT=1 

(The distributedcombigridmodule is completely independent from the other modules 
now. It is not necessary any more to compile or link the combigrid module.)

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



