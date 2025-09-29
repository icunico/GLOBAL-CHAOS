#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
add_library(optim
            CG+_mpi/cgfam.f
            CG+_mpi/cgsearch.f
            CG+/cgfam.f
            CG+/cgsearch.f
            LBFGS/lbfgsb.f
            LBFGS/linpack.f
            LBFGS/timer.f
           )
#KBtarget_include_directories(optim
#KB    PUBLIC
#KB        $<INSTALL_INTERFACE:include>
#KB        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
#KB)

# In principle the code should be fixed - but here we provide a working fix
# Set -fallow-argument-mismatch for major version >= 10
#KBif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#KBstring(REPLACE "." ";" VERSION_LIST ${CMAKE_Fortran_COMPILER_VERSION})
#KBlist(GET VERSION_LIST 0 VERSION_MAJOR)
#KBif(${VERSION_MAJOR} GREATER "9")
#KBtarget_compile_options(pdaf-d PRIVATE "-fallow-argument-mismatch")
#KBendif()
#KBendif()
