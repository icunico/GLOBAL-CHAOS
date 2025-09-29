set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include )

add_library(adj STATIC
   bioptimod_memory.f90 fm34.f90 lidata_test.F90 SLAE.F90 adj_new.f90
)

add_executable(adj_exe
   main_adj.f90
)
set_property(TARGET adj_exe PROPERTY OUTPUT_NAME adj.xx)

target_link_libraries(adj_exe PRIVATE adj)

add_executable(forward_exe
   main_forward.f90
)
set_property(TARGET forward_exe PROPERTY OUTPUT_NAME forward.xx)

target_link_libraries(forward_exe PRIVATE adj)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(adj PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
