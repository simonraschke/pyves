cmake_minimum_required(VERSION 3.10)
project(_pves)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

set(CMAKE_CXX_STANDARD 20)
set(PYBIND11_CPP_STANDARD -std=c++20)

message("C   COMPILER " ${CMAKE_C_COMPILER})
message("CXX COMPILER " ${CMAKE_CXX_COMPILER})

find_package(Python3 COMPONENTS Interpreter Development  REQUIRED)
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

SET(SOURCE_DIR "pyves/src")
include_directories(${SOURCE_DIR})

# file(GLOB SOURCES ${SOURCE_DIR}/*.cpp)

set(TEST_DIR "test")
file(GLOB TEST_FILES ${TEST_DIR}/*.cpp)

add_subdirectory(lib/catch)
include_directories(lib/catch/single_include/catch2)
add_executable("${PROJECT_NAME}_test" ${TEST_FILES})
target_compile_options("${PROJECT_NAME}_test" PRIVATE -fmax-errors=3)
# target_link_options("${PROJECT_NAME}_test" PRIVATE -W -Wno-undef)
target_link_libraries("${PROJECT_NAME}_test" Catch2::Catch2)

add_subdirectory(lib/tbb)
add_subdirectory(lib/pybind11)

include_directories(pyves/src)
include_directories(lib/eigen)
include_directories(lib/tbb/include)
include_directories(lib/pybind11/include)

pybind11_add_module(_pyves pyves/src/bindings.cpp)
target_compile_options(_pyves PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3)
# target_link_libraries(_pyves PRIVATE tbb)