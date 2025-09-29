cmake_minimum_required(VERSION 3.12)

# Create main project.
project(pdaf VERSION 1.16.0 LANGUAGES Fortran)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Not used now - KB
if(USE_GIT_VERSION_INFO)
   set(GIT_COMMIT_ID ${pdaf_VERSION})
   set(GIT_BRANCH_NAME unknown)
   option(PDAF_EMBED_VERSION "Embed PDAF version information" OFF)
   if(PDAF_EMBED_VERSION)
      add_custom_target(version
         ${CMAKE_COMMAND} -DINFILE=${PROJECT_SOURCE_DIR}/util/pdaf_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/pdaf_version.F90 -P "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Retrieving description of last PDAF commit..."
        SOURCES ${PROJECT_SOURCE_DIR}/util/pdaf_version.F90.in "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
        VERBATIM
      )
      set_property(TARGET version PROPERTY FOLDER pdaf)
   endif()
endif(USE_GIT_VERSION_INFO)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   # To find MKL, see https://cmake.org/cmake/help/latest/module/FindBLAS.html#intel-mkl
   enable_language(C)
   set(BLA_VENDOR Intel10_64lp_seq)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

option(PDAF_USE_MPI "Enable MPI parallelization" ON)
if(PDAF_USE_MPI AND NOT ${MPI_FOUND})
   find_package(MPI REQUIRED)
   message(STATUS ${MPI_FOUND})
else()
   message(STATUS "PDAF: MPI is already configured - continuing")
endif(PDAF_USE_MPI AND NOT ${MPI_FOUND})

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

add_subdirectory(src)
add_subdirectory(external)
#KBadd_subdirectory(testsuite/src/dummymodel_1D)
#KBadd_subdirectory(test)
