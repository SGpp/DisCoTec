# Contributing

## Report a Bug

If you encounter a bug, please write an [issue](https://github.com/SGpp/DisCoTec/issues),
describing the system and environment where you encountered the error.
If you are using Spack, please include the output of `spack spec discotec`
(to your installed discotec version) and/or `spack load --list` if you are
using spack modules.

### Share Compilation Pain

If you find yourself in the very common (but sometimes seemingly endless)
pain of compilation on HPC systems, feel free to open an issue as well.
If it is a problem in the DisCoTec Spack package or CMake setup, our
contributors may be able to solve it, otherwise they may be able to give
pointers to common solutions.

## Solve an Issue or Contribute a Feature

To find a good place to start coding, have a look at the
[issue tracker](https://github.com/SGpp/DisCoTec/issues) and comment if
you start working on an issue.

If you want to implement a new feature not mentioned in the issues,
please also create an issue to prevent duplicate work.

DisCoTec uses the Google bracket style (2 spaces) and camelCase naming
convention; member variables are appended with `_`.

If you are not affiliated with the SGpp Team, we typically expect you
develop the code on your fork and create a pull request for review.
Please include (unit) tests in the [./tests](/tests) folder.
Reviews will be based on function, understandability, and testability /
test coverage.

## Something else? Get in touch!

If you want to contribute in another way or collaborate with The SGpp Team, get in touch at <theresa@sparsegrids.org>.
