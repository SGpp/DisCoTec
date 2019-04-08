Building the SimpleAmqpClient C++ Library
=========================================

The library building instructions are described in the official git repo
(https://github.com/alanxz/SimpleAmqpClient). However we need to make some
adjustments since installing the dependency rabbitmq-c
(https://github.com/alanxz/rabbitmq-c) systemwide is not possible.

The Building steps are the following:

1. Clone and build rabbitmq-c as described in the git repo.
2. Clone the SimpleAmqpClient Library, create a build folder in the main
   directory and enter it. (`mkdir build && cd build`)
4. Run `cmake ..` once. This will throw errors since it cannot find the
   rabbitmq-c (maybe also boost) libs and headers).
3. Use `ccmake ..` to set the corrseponding paths
   (use c to save your configuration). For rabbitmq-c the includes are under
   rabbitmq-c/librabbitmq and the library is under
   rabbitmq-c/build/librabbitmq/librabbitmq.so
4. Use `cmake ..` again, this time there should not be any error.
5. Use `make` to build the library.

4. (Optional) The tests are not build by default. You have to clone gtest from
   https://github.com/abseil/googletest/tree/ec44c6c1675c25b9827aacd08c02433cccde7780
   into SimpleAmqpClient/third-party/ and pass `-DENABLE_TESTING=ON` to cmake.
   But I haven't yet figured out how to successfully run the tests.


