cmake_minimum_required(VERSION 3.10)
project(_pves LANGUAGES C CXX)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

set(CMAKE_CXX_STANDARD 17)
set(PYBIND11_CPP_STANDARD -std=c++17)

message("C   COMPILER " ${CMAKE_C_COMPILER})
message("CXX COMPILER " ${CMAKE_CXX_COMPILER})

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

SET(SOURCE_DIR "pyves/src")
include_directories(${SOURCE_DIR})

set(TEST_DIR "test")
file(GLOB TEST_FILES ${TEST_DIR}/*.cpp)

add_subdirectory(${CMAKE_SOURCE_DIR}/lib/catch)
include_directories(${CMAKE_SOURCE_DIR}/lib/catch/single_include/catch2)

add_executable("${PROJECT_NAME}_test" ${TEST_FILES} pyves/src/random.cpp)
target_compile_options("${PROJECT_NAME}_test"
    PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CXX_COMPILER_ID:Clang>:-ggdb -g -Wall -Wextra -Wpedantic -ferror-limit=2>
        $<$<CXX_COMPILER_ID:GNU>:-ggdb -g -Wall -Wextra -Wpedantic -fmax-errors=2>
)
target_link_libraries("${PROJECT_NAME}_test" Catch2::Catch2 tbb pthread)

add_subdirectory(${CMAKE_SOURCE_DIR}/lib/pybind11)

include_directories(${CMAKE_SOURCE_DIR}/lib/eigen)
include_directories(${CMAKE_SOURCE_DIR}/lib/tbb/include)
include_directories(${CMAKE_SOURCE_DIR}/lib/pybind11/include)

file(GLOB SOURCES ${SOURCE_DIR}/*.cpp)
pybind11_add_module(_pyves ${SOURCES})
target_compile_options(_pyves
    PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CXX_COMPILER_ID:Clang>:-ggdb -g -Wall -Wextra -Wpedantic -ferror-limit=2>
        $<$<CXX_COMPILER_ID:GNU>:-ggdb -g -Wall -Wextra -Wpedantic -fmax-errors=2>
)
target_link_libraries(_pyves PRIVATE tbb pthread)