cmake_minimum_required (VERSION 2.8)
project (coalescence Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()


set (CMAKE_NO_SYSTEM_FROM_IMPORTED yes)

set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -O3 -mcmodel=medium")

set (CMAKE_Fortran_MODULE_DIRECTORY mod)
file (GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.f*)
add_executable(coalesce ${_sources})


# List of all libraries

find_package(FSON REQUIRED)
include_directories(${FSON_INCLUDE_DIRS})

include_directories($ENV{HOME}/include)
find_library(FSON fson PATHS $ENV{HOME}/lib)

message(${_sources})
message(${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
    coalesce
    ${FSON_LIBRARIES}
  )
