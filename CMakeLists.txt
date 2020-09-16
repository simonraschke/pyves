cmake_minimum_required(VERSION 3.14)
project(_pves LANGUAGES C CXX)

set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_CXX_STANDARD 17)
set(PYBIND11_CPP_STANDARD -std=c++17)

message("C   COMPILER " ${CMAKE_C_COMPILER})
message("CXX COMPILER " ${CMAKE_CXX_COMPILER})



# FIND MODULES
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "PYTHON_LIBRARIES    = ${Python3_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE   = ${Python3_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON3_INCLUDE_DIRS}")
include_directories(${Python3_INCLUDE_DIRS})

set(EIGEN_DEFINES "")
set(EIGEN_LINKING_EXTENDED "")
# find_package(BLAS QUIET)
# if(BLAS_FOUND)
#     message(STATUS "found BLAS")
#     set(EIGEN_DEFINES "EIGEN_USE_BLAS")
#     set(EIGEN_LINKING_EXTENDED ${BLAS_LIBRARIES})
# endif()

message(STATUS "EIGEN_DEFINES = ${EIGEN_DEFINES}")
message(STATUS "EIGEN_LINKING_EXTENDED = ${EIGEN_LINKING_EXTENDED}")
add_compile_definitions(${EIGEN_DEFINES})


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
include_directories(${CMAKE_SOURCE_DIR}/lib/pybind11/include)
include_directories(${CMAKE_SOURCE_DIR}/lib/taskflow)



# MAKE TEST EXECUTABLE
# if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#     message("debug mode. building tests")
    add_executable("${PROJECT_NAME}_test" ${TEST_FILES})
    target_include_directories(${PROJECT_NAME}_test PRIVATE ${CMAKE_SOURCE_DIR}/lib/catch/single_include/catch2)
    target_compile_options(${PROJECT_NAME}_test
        PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/W3>
            $<$<CXX_COMPILER_ID:Clang>:-ggdb -g -Wall -Wextra -Wpedantic -ferror-limit=2>
            $<$<CXX_COMPILER_ID:GNU>:-ggdb -g -Wall -Wextra -Wpedantic -fmax-errors=2>
    )
    message(STATUS ${Python_VERSION})
    message("Python Major Version " ${Python3_VERSION_MAJOR})
    message("Python Minor Version " ${Python3_VERSION_MINOR})
    # if(Python_MINOR GREATER_EQUAL 8)
    #     message(">=8")
    #     execute_process(COMMAND python3-config --embed --ldflags OUTPUT_VARIABLE PY_EMBED_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    # else()
    #     message("<8")
    #     execute_process(COMMAND python3-config --ldflags OUTPUT_VARIABLE PY_EMBED_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    # endif()
    # string(STRIP "${PY_EMBED_LIBS}" PY_EMBED_LIBS)
    target_link_options(${PROJECT_NAME}_test PRIVATE -W -Wno-undef)
    target_link_libraries(${PROJECT_NAME}_test PRIVATE pthread ${Python3_LIBRARIES} ${EIGEN_LINKING_EXTENDED})
# endif()


# MAKE PYTHON SHARED LIBRARY
pybind11_add_module(_pyves ${SRC_FILES})
set_target_properties(_pyves PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_compile_options(_pyves
    PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CXX_COMPILER_ID:Clang>:-ggdb -g -Wall -Wextra -Wpedantic -ferror-limit=2>
        $<$<CXX_COMPILER_ID:GNU>:-ggdb -g -Wall -Wextra -Wpedantic -fmax-errors=2>
)
target_link_options(_pyves PRIVATE -W -Wno-undef)
target_link_libraries(_pyves PRIVATE pthread ${EIGEN_LINKING_EXTENDED})