cmake_minimum_required(VERSION 3.10)
project(_pves LANGUAGES C CXX)

set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_CXX_STANDARD 17)
set(PYBIND11_CPP_STANDARD -std=c++17)

message("C   COMPILER " ${CMAKE_C_COMPILER})
message("CXX COMPILER " ${CMAKE_CXX_COMPILER})



# FIND MODULES
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})



#SET VARIABLES
SET(TEST_DIR test)
SET(SRC_DIR pyves/src)
file(GLOB SRC_FILES ${SRC_DIR}/*.cpp ) #maybe remove bindings.cpp
message(${SRC_FILES})
string(REPLACE "${SRC_DIR}/bindings.cpp" " " SRC_FILES_NO_BIND "${SRC_FILES}")
message(${SRC_FILES_NO_BIND})
file(GLOB TEST_FILES ${SRC_FILES} ${TEST_DIR}/*.cpp)



# ADD ALL SUBDIRS
add_subdirectory(${CMAKE_SOURCE_DIR}/lib/pybind11)



# INCLUDE LIB DIRS
include_directories(${SRC_DIR})
include_directories(${CMAKE_SOURCE_DIR}/lib/eigen)
include_directories(${CMAKE_SOURCE_DIR}/lib/tbb/include)
include_directories(${CMAKE_SOURCE_DIR}/lib/pybind11/include)



# MAKE TEST EXECUTABLE
add_executable("${PROJECT_NAME}_test" ${TEST_FILES})
target_include_directories(${PROJECT_NAME}_test PRIVATE ${CMAKE_SOURCE_DIR}/lib/catch/single_include/catch2)
target_compile_options(${PROJECT_NAME}_test
    PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CXX_COMPILER_ID:Clang>:-ggdb -g -Wall -Wextra -Wpedantic -ferror-limit=2>
        $<$<CXX_COMPILER_ID:GNU>:-ggdb -g -Wall -Wextra -Wpedantic -fmax-errors=2>
)
execute_process(COMMAND python3-config --embed --ldflags OUTPUT_VARIABLE PY_EMBED_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
string(STRIP "${PY_EMBED_LIBS}" PY_EMBED_LIBS)
target_link_options(${PROJECT_NAME}_test PRIVATE -W -Wno-undef)
target_link_libraries(${PROJECT_NAME}_test PRIVATE tbb pthread ${PY_EMBED_LIBS})



# MAKE PYTHON SHARED LIBRARY
pybind11_add_module(_pyves ${SRC_FILES})
target_compile_options(_pyves
    PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CXX_COMPILER_ID:Clang>:-ggdb -g -Wall -Wextra -Wpedantic -ferror-limit=2>
        $<$<CXX_COMPILER_ID:GNU>:-ggdb -g -Wall -Wextra -Wpedantic -fmax-errors=2>
)
target_link_options(_pyves PRIVATE -W -Wno-undef)
target_link_libraries(_pyves PRIVATE tbb pthread)