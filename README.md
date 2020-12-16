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

Requirements
--------------
scons (Python2, >= 2.4.1-1)
libboost-all-dev (>= 1.60)
libmpich-dev (>= 3.2-6)


Installation instructions: 
--------------------------
Compile by running
```
./compile.sh
``` 
and for debugging, add
```
CPPFLAGS=-g
``` 
to the scons command. To clean the compiled files use

```
scons -c
```


Gene submodules:
----------------
There are gene versions as submodules: a linear one in the gene_mgr folder, and 
a nonlinear one in gene-non-linear. To get them, you need access to their 
respective repos. Then go into the folder and

```
git submodule init
git submodule update
```
and then you just hope you can `make` it by adapting the local library flags ;)
