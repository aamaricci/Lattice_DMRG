# Install script for directory: /Users/amaricci/projects_src/Lattice_DMRG

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/include/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/include" TYPE DIRECTORY FILES "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/include/")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND "/usr/local/bin/cmake" -E rm -f )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/lib/libdmrg.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/lib" TYPE STATIC_LIBRARY FILES "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/libdmrg.a")
  if(EXISTS "$ENV{DESTDIR}/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/lib/libdmrg.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/lib/libdmrg.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/lib/libdmrg.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/CMakeFiles/dmrg.dir/install-cxx-module-bmi-RELEASE.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/etc/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/etc" TYPE DIRECTORY FILES "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/etc/")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/bin/dmrg_config_user.sh")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/bin" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/etc/dmrg_config_user.sh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/bin/dmrg_config_global.sh")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/bin" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/etc/dmrg_config_global.sh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/version")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/version")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/.pkgconfig.d/dmrg.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/.pkgconfig.d" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/etc/dmrg.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/amaricci/.modules.d/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/amaricci/.modules.d" TYPE DIRECTORY FILES "/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/etc/modules/")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE(
"
ADD LIBRARY TO YOUR SYSTEM: 
Pick ONE method below [or add it in your bash profile, e.g. ~/.bashrc]:
[33mMethod 1: use the provided dmrg environment module[m:
   $ module use $HOME/.modules.d
   $ module load dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d

[33mMethod 2: source the config script[m:
   $ source /Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/bin/dmrg_config_user.sh

[33mMethod 3: use pkg-config with the provided dmrg.pc[m:
   $ export PKG_CONFIG_PATH=/Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/etc/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs dmrg

[33mMethod ADMIN: Add this line to the system shell configuration file, e.g. /etc/bash.bashrc[m
   $ source /Users/amaricci/opt/dmrg/gnu/devel_mpi/dble/3.2.1-11-gc8aed3d/bin/dmrg_config_global.sh
")

endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/src/MATRIX/cmake_install.cmake")
  include("/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/src/LIST/cmake_install.cmake")
  include("/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/src/FOCK/cmake_install.cmake")
  include("/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/src/BLOCKS/cmake_install.cmake")
  include("/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/src/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/amaricci/projects_src/Lattice_DMRG/build_dev_mpi_release/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
