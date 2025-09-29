cmake_minimum_required(VERSION 3.10)

project(adj VERSION 1.0.0 LANGUAGES Fortran)

add_subdirectory(src)


enable_testing()

# define tests
add_test(
  NAME check_output
  COMMAND bash ${CMAKE_SOURCE_DIR}/tests/test01.sh $<TARGET_FILE:adj_exe>
)

