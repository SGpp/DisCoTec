# hyper.deal for discontinous Galerkin Combination Technique

## Authors

* Peter Munch @ga69jux
* Katharina Kormann @ne65nib
* Martin Kronbichler @ne96pad

## Installation

```bash
# get the matrixfree-branch of deal.II
git clone https://gitlab.lrz.de/ga69jux/matrixfree.git
cd matrixfree 
git checkout for_hyperdeal
cd ..

# get and build metis
git clone https://github.com/scibuilder/metis.git
cd metis
cmake .
make 
cd ..

# get and build p4est
wget http://p4est.github.io/release/p4est-2.0.tar.gz
matrixfree/doc/external-libs/p4est-setup.sh p4est-2.0.tar.gz `pwd`

cd matrixfree
git checkout 095600820d999db798a43530662d384425bd4db0
cd ../matrixfree-build/
cmake -D CMAKE_CXX_FLAGS="-std=c++17 -march=native"     -D DEAL_II_CXX_FLAGS_RELEASE="-std=c++17 -O3"     -D CMAKE_BUILD_TYPE="DebugRelease"     -D CMAKE_INSTALL_PREFIX="../matrixfree-install/"     -D DEAL_II_WITH_MPI="ON"   -D DEAL_II_WITH_LAPACK="ON"   -D DEAL_II_WITH_THREADS="OFF"   -D DEAL_II_WITH_P4EST="ON"   -D P4EST_DIR="../FAST"   -D DEAL_II_WITH_METIS:BOOL="ON"   -D METIS_DIR:FILEPATH="../metis" -D  ../matrixfree
make install
cd ../hyper.deal.combi
git checkout d91fc0731e6ce855660b7a5d0a8092cc6ce63461
git patch /home/sgsscratch/hyperdeal.patch
cd ../hyper.deal.combi-build
cmake ../hyper.deal.combi -D DEAL_II_DIR=../matrixfree-install -D CMAKE_CXX_FLAGS="-std=c++17 -O3" -D CMAKE_CXX_COMPILER="mpicxx"
make release
cd applications/advection_reference_dealii
make 
ins compile.sh noch DOC=0 hinzufügen
Discotec branch projektarbeit
. ./compile.sh
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/libboost_serialization.so:/home/sgsscratch/goehrims/matrixfree-install/lib:/home/sgsscratch/goehrims/DisCoTec/lib/sgpp:/home/sgsscratch/goehrims/DisCoTec/glpk/lib
cd distributedcombigrid/examples/dealii_example
make
cd ../../tools/
make