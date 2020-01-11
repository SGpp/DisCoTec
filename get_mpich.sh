if [ -f mpich/lib/libmpich.so ]; then
  echo "libmpich.so found -- nothing to build."
else
  echo "Downloading mpich source."
  wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
  tar xfz mpich-3.2.tar.gz
  rm mpich-3.2.tar.gz
  echo "configuring and building mpich."
  cd mpich-3.2
  ./configure \
          --prefix=`pwd`/../mpich \
          --enable-static=false \
          --enable-alloca=true \
          --disable-long-double \
          --enable-threads=single \
          --enable-fortran=no \
          --enable-fast=all \
          --enable-g=none \
          --enable-timing=none
  make -j4
  make install
  cd -
  rm -rf mpich-3.2
fi


echo "Downloading boost from source."
wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
tar -zxvf boost_1_69_0.tar.gz
rm boost_1_69_0.tar.gz
cd boost_1_69_0
./bootstrap.sh
./b2
cd -
rm -rf boost_1_69_0
