# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import distutils.sysconfig
import os
import subprocess
import re

import Helper

def doConfigure(env, moduleFolders, languageWrapperFolders):
  print()
  print("Checking programs and libraries:")

  config = env.Configure(custom_tests={"CheckExec" : Helper.CheckExec,
                                       "CheckJNI" : Helper.CheckJNI,
                                       "CheckFlag" : Helper.CheckFlag})

  # now set up all further environment settings that should never fail
  # compiler setup should be always after checking headers and flags,
  # as they can make the checks invalid, e.g., by setting "-Werror"

  # TODO(pfandedd): discuss with julian
#   if env["PLATFORM"] not in ["cygwin", "win32"]:
  if env["OPT"] == True:
    env.Append(CPPFLAGS=["-DNDEBUG", "-O3"])
  else:
    env.Append(CPPFLAGS=["-g", "-O0"])

  # make settings case-insensitive
  env["COMPILER"] = env["COMPILER"].lower()
  env["ARCH"] = env["ARCH"].lower()

  if env["COMPILER"] in ("gnu", "openmpi", "mpich"):
    configureGNUCompiler(config)
  elif env["COMPILER"] == "clang":
    configureClangCompiler(config)
  elif env["COMPILER"] in ("intel", "intel.mpi"):
    configureIntelCompiler(config)
  else:
    Helper.printErrorAndExit("You must specify a valid value for Compiler.",
                             "Available configurations are:",
                             "gnu, clang, intel, openmpi, mpich, intel.mpi")

  # special treatment for different platforms
  if env["PLATFORM"] == "darwin":
    # the "-undefined dynamic_lookup"-switch is required to actually build a shared library
    # in OSX. "-dynamiclib" alone results in static linking of
    # all further dependent shared libraries
    # beware: if symbols are missing that are actually required
    # (because the symbols don't reside in a shared library),
    # there will be no error during compilation
    # the python binding (pysgpp) requires lpython and a flat namespace
    # also for the python binding, the library must be suffixed with '*.so' even
    # though it is a dynamiclib and not a bundle (see SConscript in src/pysgpp)
    env.AppendUnique(LINKFLAGS=["-flat_namespace", "-undefined", "dynamic_lookup", "-lpython"])
    # The GNU assembler (GAS) is not supported in Mac OS X.
    # A solution that fixed this problem is by adding -Wa,-q to the compiler flags.
    # From the man pages for as (version 1.38):
    # -q Use the clang(1) integrated assembler instead of the GNU based system assembler.
    # Note that the CPPFLAG is exactly "-Wa,-q", where -Wa passes flags to the assembler and
    # -q is the relevant flag to make it use integrated assembler
    env.AppendUnique(CPPFLAGS=["-Wa,-q"])
    env.AppendUnique(CPPPATH="/usr/local/include")
    env.AppendUnique(LIBPATH="/usr/local/lib")
    env["SHLIBSUFFIX"] = ".dylib"
  elif env["PLATFORM"] == "cygwin":
    # required to find the static libraries compiled before the shared libraries
    # the static libraries are required as the linker on windows cannot ignore undefined symbols
    # (as is done on linux automatically and is done on OSX with the settings above)
    #env.Append(LIBPATH=[BUILD_DIR])
    pass

  # will lead to a warning on cygwin and win32
  # ("all code is position independent")
  if env["PLATFORM"] not in ["cygwin", "win32"]:
    env.AppendUnique(CPPFLAGS=["-fPIC"])

  # setup the include base folder
  # env.Append(CPPPATH=["#/src/sgpp"])
  for moduleFolder in moduleFolders:
    if (moduleFolder in languageWrapperFolders) or (not env["SG_" + moduleFolder.upper()]):
      continue
    env.AppendUnique(CPPPATH=["#/" + moduleFolder + "/src/"])

  # detour compiler output
  env["PRINT_CMD_LINE_FUNC"] = Helper.printCommand

  checkCpp17(config)

  # glpk library
  config.env.AppendUnique(CPPPATH=[config.env['GLPK_INCLUDE_PATH']])
  config.env.AppendUnique(LIBPATH=[config.env['GLPK_LIBRARY_PATH']])

  checkBoostTests(config)
  env = config.Finish()

  print("Configuration done.")
  print()

def checkCpp17(config):
  # check C++17 support
  if not config.CheckFlag(""):
    Helper.printErrorAndExit("The compiler doesn't seem to receive sensible CPPFLAGS.\
      Maybe you forgot to set OPT? Abort!")
    Exit(1)
  if not config.CheckFlag("-std=c++17"):
    Helper.printErrorAndExit("The compiler doesn't seem to support the C++17 standard. Abort!")
    Exit(1)

  config.env.PrependUnique(CPPFLAGS="-std=c++17")
def checkBoostTests(config):
  config.env.AppendUnique(CPPPATH=[config.env["BOOST_INCLUDE_PATH"]])
  config.env.AppendUnique(LIBPATH=[config.env["BOOST_LIBRARY_PATH"]])
  # Check the availability of the boost unit test dependencies
  if config.env["COMPILE_BOOST_TESTS"]:
    config.env.AppendUnique(CPPPATH=[config.env["BOOST_INCLUDE_PATH"]])
    config.env.AppendUnique(LIBPATH=[config.env["BOOST_LIBRARY_PATH"]])

    if not config.CheckHeader(os.path.join("boost", "test", "unit_test.hpp"), language="c++"):
      config.env["COMPILE_BOOST_TESTS"] = False
      Helper.printWarning("****************************************************",
                          "No Boost Unit Test Headers found. Omitting Boost unit tests.",
                          "Please install the corresponding package, e.g., on Ubuntu",
                          "> sudo apt-get install libboost-test-dev",
                          "****************************************************")

    if not config.CheckLib("boost_unit_test_framework", autoadd=0, language="c++"):
      config.env["COMPILE_BOOST_TESTS"] = False
      Helper.printWarning("****************************************************",
                          "No Boost Unit Test library found. Omitting Boost unit tests.",
                          "Please install the corresponding package, e.g., on Ubuntu",
                          "> sudo apt-get install libboost-test-dev",
                          "****************************************************")

def configureGNUCompiler(config):
  if config.env["COMPILER"] == "openmpi":
    config.env["CC"] = ("mpicc.openmpi")
    config.env["LINK"] = ("mpic++.openmpi")
    config.env["CXX"] = ("mpic++.openmpi")
    config.env["FC"] = ("mpifort.openmpi")
    config.env["CPPDEFINES"]["USE_MPI"] = "1"
    Helper.printInfo("Using openmpi.")
  elif config.env["COMPILER"] == "mpich":
    config.env["CC"] = ("mpicc.mpich")
    config.env["LINK"] = ("mpicxx.mpich")
    config.env["CXX"] = ("mpicxx.mpich")
    config.env["FC"] = ("mpifort.mpich")
    config.env["CPPDEFINES"]["USE_MPI"] = "1"
    Helper.printInfo("Using mpich.")

  versionString = subprocess.check_output([config.env["CXX"], "-dumpversion"]).strip()
  Helper.printInfo("Using {} ({})".format(config.env["CXX"], versionString))

  if not config.CheckExec(config.env["CXX"]) or not config.CheckExec(config.env["CC"]) or \
      not config.CheckExec(config.env["LINK"]) :
    Helper.printErrorAndExit("Compiler not found!")

  allWarnings = \
      "-Wall -pedantic -Wextra \
      -Wconversion -Wformat=2 \
      -Wformat-nonliteral -Wformat-security -Winit-self  \
      -Wmissing-format-attribute \
      -Wmissing-include-dirs -Wpacked -Wredundant-decls \
      -Wswitch-default -Wswitch-enum -Wunreachable-code -Wunused \
      -Wno-unused-parameter".split(" ")

  # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or
  #     http://www.swig.org/Release/CHANGES, 03/02/2006
  #     "If you are going to use optimizations turned on with gcc > 4.0 (for example -O2),
  #     ensure you also compile with -fno-strict-aliasing"
  config.env.Append(CPPFLAGS=allWarnings + [
      "-fno-strict-aliasing",
      "-funroll-loops", "-mfpmath=sse",
      # "-fsanitize=address",
      # "-DDEBUG_OUTPUT=1",
      ])
  config.env.Append(CPPFLAGS=["-fopenmp"])
  config.env.Append(LINKFLAGS=[
    "-fopenmp",
    # "-fsanitize=leak",
    # "-fsanitize=address"#,"-lasan","-lubsan" ## ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH=$(shell which llvm-symbolizer)# cf https://gist.github.com/kwk/4171e37f4bcdf7705329
    ])

  # required for profiling
  config.env.Append(CPPFLAGS=["-fno-omit-frame-pointer"])

  print("config: "+ str(config.env))

  if config.env["DEBUG_OUTPUT"]:
    config.env.Append(CPPFLAGS=["-DDEBUG_OUTPUT"])
  if config.env["TIMING"]:
    config.env.Append(CPPFLAGS=["-DTIMING"])
  if config.env["BUILD_STATICLIB"]:
    config.env.Append(CPPFLAGS=["-D_BUILD_STATICLIB"])

  if config.env["USENONBLOCKINGMPICOLLECTIVE"]:
    config.env.Append(CPPFLAGS=["-DUSENONBLOCKINGMPICOLLECTIVE"])
  if config.env["OMITREADYSIGNAL"]:
    config.env.Append(CPPFLAGS=["-DOMITREADYSIGNAL"])
  if config.env["UNIFORMDECOMPOSITION"]:
    config.env.Append(CPPFLAGS=["-DUNIFORMDECOMPOSITION"])
  if config.env["USE_HDF5"]:
    config.env.Append(CPPFLAGS=["-DHAVE_HIGHFIVE"])
  if config.env["ENABLEFT"]:
    config.env.Append(CPPFLAGS=["-DENABLEFT"])
  if config.env["ISGENE"]:
    config.env.Append(CPPFLAGS=["-DISGENE"])
  if config.env["USE_VTK"]:
    config.env.Append(CPPFLAGS=["-DUSE_VTK"])
    config.env.Append(LINKFLAGS=[
      "-lvtkCommonCore",
      "-lvtkIOXML",
      "-lvtkIOParallelXML"
    ])
    config.env.AppendUnique(CPPPATH=[config.env['VTK_INCLUDE_PATH']])

  if config.env["ARCH"] == "sse3":
    config.env.AppendUnique(CPPFLAGS=["-msse3"])
  elif config.env["ARCH"] == "sse42":
    config.env.AppendUnique(CPPFLAGS=["-msse4.2"])
  elif config.env["ARCH"] == "avx":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
  elif config.env["ARCH"] == "fma4":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
    config.env.AppendUnique(CPPFLAGS=["-mfma4"])
  elif config.env["ARCH"] == "avx2":
    config.env.AppendUnique(CPPFLAGS=["-mavx2"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  elif config.env["ARCH"] == "avx512":
    config.env.AppendUnique(CPPFLAGS=["-mavx512f"])
    config.env.AppendUnique(CPPFLAGS=["-mavx512cd"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  else:
    Helper.printErrorAndExit("You must specify a valid ARCH value for gnu.",
                             "Available configurations are: sse3, sse42, avx, fma4, avx2, avx512")

  # check if using MinGW (g++ on win32)
  if config.env["PLATFORM"] == "win32":
    # disable warnings which occur when including Boost in the tests
    config.env.Append(CPPFLAGS=["-Wno-switch-enum", "-Wno-deprecated-declarations"])
    # also use "lib" prefix on MinGW for consistency with Linux (default is no prefix)
    config.env["SHLIBPREFIX"] = "lib"

def configureClangCompiler(config):
  config.env["CC"] = ("clang")
  config.env["LINK"] = ("clang++")
  config.env["CXX"] = ("clang++")

  versionString = subprocess.check_output([config.env["CXX"], "--version"])
  try:
    versionString = re.match(r"^.*version ([^ ]+)", versionString).group(1)
  except AttributeError:
    versionString = "(unknown version)"
  Helper.printInfo("Using {} {}".format(config.env["CXX"], versionString))

  if not config.CheckExec(config.env["CXX"]) or not config.CheckExec(config.env["CC"]) or \
      not config.CheckExec(config.env["LINK"]) :
    Helper.printErrorAndExit("Compiler not found!")

  allWarnings = "-Wall -Wextra -Wno-unused-parameter".split(" ")

  # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or
  #     http://www.swig.org/Release/CHANGES, 03/02/2006
  #    "If you are going to use optimisations turned on with gcc > 4.0 (for example -O2),
  #     ensure you also compile with -fno-strict-aliasing"

  config.env.Append(CPPFLAGS=allWarnings + [
      "-DDEFAULT_RES_THRESHOLD=-1.0", "-DTASKS_PARALLEL_UPDOWN=4"])
  config.env.Append(CPPFLAGS=["-fopenmp=libiomp5"])
  config.env.Append(LINKFLAGS=["-fopenmp=libiomp5"])

  if config.env["DEBUG_OUTPUT"]:
    config.env.Append(CPPFLAGS=["-DDEBUG_OUTPUT"])


  if config.env["BUILD_STATICLIB"]:
    config.env.Append(CPPFLAGS=["-D_BUILD_STATICLIB"])

  if config.env["ARCH"] == "sse3":
    config.env.AppendUnique(CPPFLAGS=["-msse3"])
  elif config.env["ARCH"] == "sse42":
    config.env.AppendUnique(CPPFLAGS=["-msse4.2"])
  elif config.env["ARCH"] == "avx":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
  elif config.env["ARCH"] == "fma4":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
    config.env.AppendUnique(CPPFLAGS=["-mfma4"])
  elif config.env["ARCH"] == "avx2":
    config.env.AppendUnique(CPPFLAGS=["-mavx2"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  elif config.env["ARCH"] == "avx512":
    config.env.AppendUnique(CPPFLAGS=["-mavx512f"])
    config.env.AppendUnique(CPPFLAGS=["-mavx512cd"])
    config.env.AppendUnique(CPPFLAGS=["-mfma"])
  else:
    Helper.printErrorAndExit("You must specify a valid ARCH value for clang.",
                             "Available configurations are: sse3, sse4.2, avx, fma4, avx2, avx512")

def configureIntelCompiler(config):
  config.env.AppendUnique(CPPFLAGS=["-Wall", "-ansi", "-Wno-deprecated", "-wd1125",
                                    "-fno-strict-aliasing",
                                    "-ip", "-ipo", "-funroll-loops",
                                    "-ansi-alias", "-fp-speculation=safe",
                                    "-DDEFAULT_RES_THRESHOLD=-1.0", "-DTASKS_PARALLEL_UPDOWN=4",
                                    "-no-offload"])
  if config.env["COMPILER"] == "intel.mpi":
    config.env["CC"] = ("mpiicc")
    config.env["LINK"] = ("mpiicpc")
    config.env["CXX"] = ("mpiicpc")
    config.env["CPPDEFINES"]["USE_MPI"] = "1"
  else:
    config.env["CC"] = ("icc")
    config.env["LINK"] = ("icpc")
    config.env["CXX"] = ("icpc")

  versionString = subprocess.check_output([config.env["CXX"], "-dumpversion"]).strip()
  Helper.printInfo("Using {} {}".format(config.env["CXX"], versionString))

  if not config.CheckExec(config.env["CXX"]) or not config.CheckExec(config.env["CC"]) or \
      not config.CheckExec(config.env["LINK"]) :
    Helper.printErrorAndExit("Compiler not found!")

  config.env.AppendUnique(CPPFLAGS=["-qopenmp"])
  config.env.AppendUnique(LINKFLAGS=["-qopenmp"])

  if config.env["DEBUG_OUTPUT"]:
    config.env.Append(CPPFLAGS=["-DDEBUG_OUTPUT"])


  if config.env["BUILD_STATICLIB"]:
    config.env.AppendUnique(CPPFLAGS=["-D_BUILD_STATICLIB"])

  if config.env["ARCH"] == "sse3":
    config.env.AppendUnique(CPPFLAGS=["-msse3"])
  elif config.env["ARCH"] == "sse42":
    config.env.AppendUnique(CPPFLAGS=["-msse4.2"])
  elif config.env["ARCH"] == "avx":
    config.env.AppendUnique(CPPFLAGS=["-mavx"])
  elif config.env["ARCH"] == "avx2":
    config.env.AppendUnique(CPPFLAGS=["-xCORE-AVX2"])
    config.env.AppendUnique(CPPFLAGS=["-fma"])
  elif config.env["ARCH"] == "avx512":
    config.env.AppendUnique(CPPFLAGS=["-xCOMMON-AVX512"])
    config.env.AppendUnique(CPPFLAGS=["-fma"])
  elif config.env["ARCH"] == "mic":
    config.env.AppendUnique(CPPFLAGS=["-mmic"])
    config.env.AppendUnique(LINKFLAGS=["-mmic"])
    config.env["CPPDEFINES"]["USEMIC"] = "1"
  else:
    Helper.printErrorAndExit("You must specify a valid ARCH value for intel.",
                             "Available configurations are: sse3, sse4.2, avx, avx2, avx512, mic")

  config.env.AppendUnique(CPPPATH=[distutils.sysconfig.get_python_inc()])
