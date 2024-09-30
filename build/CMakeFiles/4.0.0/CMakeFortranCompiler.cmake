set(CMAKE_Fortran_COMPILER "/usr/local/bin/mpif90")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "14.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")


set(CMAKE_Fortran_COMPILER_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk/usr")
set(CMAKE_COMPILER_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk/usr")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/local/bin/gcc-ar-14")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_Fortran_COMPILER_LINKER "/usr/bin/ld")
set(CMAKE_Fortran_COMPILER_LINKER_ID "AppleClang")
set(CMAKE_Fortran_COMPILER_LINKER_VERSION 820.1)
set(CMAKE_Fortran_COMPILER_LINKER_FRONTEND_VARIANT GNU)
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/local/bin/gcc-ranlib-14")
set(CMAKE_TAPI "/Library/Developer/CommandLineTools/usr/bin/tapi")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED )
set(CMAKE_LINKER_PUSHPOP_STATE_SUPPORTED )
set(CMAKE_Fortran_LINKER_PUSHPOP_STATE_SUPPORTED )
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()


set(CMAKE_Fortran_SYSROOT_FLAG "-isysroot")
set(CMAKE_Fortran_OSX_DEPLOYMENT_TARGET_FLAG "-mmacosx-version-min=")

set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/usr/local/Cellar/open-mpi/5.0.7/include;/usr/local/Cellar/open-mpi/5.0.7/lib;/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/x86_64-apple-darwin21/14/finclude;/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/x86_64-apple-darwin21/14/include;/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/x86_64-apple-darwin21/14/include-fixed;/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "emutls_w;heapt_w;mpi_usempif08;mpi_usempi_ignore_tkr;mpi_mpifh;mpi;gfortran;gcc_s.1.1;gcc;quadmath")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/usr/local/Cellar/open-mpi/5.0.7/lib;/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc/x86_64-apple-darwin21/14;/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current/gcc;/Users/amaricci/opt/dmft_tools/gnu/2.7.3-1-g788688a/lib;/Users/amaricci/opt/scifor/gnu/4.20.2/lib;/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/current;/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk/System/Library/Frameworks")
