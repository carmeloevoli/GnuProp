cmake_minimum_required(VERSION 3.5)
project(gammaprop CXX)
set(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(GAMMAPROP_EXTRA_SOURCES)
set(GAMMAPROP_EXTRA_INCLUDES)
set(GAMMAPROP_EXTRA_LIBRARIES)

set(CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -pedantic -ftree-vectorize -O3")
set(GAMMAPROP_INCLUDES)
set(GAMMAPROP_LIBRARIES)

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
    #set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# PLOG (provided)
include_directories(external/plog/include)

# NamedType (provided)
include_directories(external/NamedType/include)

# GSL (required)
find_package(GSL REQUIRED)
list(APPEND GAMMAPROP_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND GAMMAPROP_EXTRA_LIBRARIES ${GSL_LIBRARIES})

# BENIAMINO (required)
find_package(SIMPROP REQUIRED)
list(APPEND GAMMAPROP_EXTRA_INCLUDES ${SIMPROP_INCLUDE_DIR})
list(APPEND GAMMAPROP_EXTRA_LIBRARIES ${SIMPROP_LIBRARIES})

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
# option(ENABLE_TESTING "Build tests and enable test target" OFF)
# if(ENABLE_TESTING)
#     set(GOOGLETEST_VERSION 1.10.0)
#     include_directories(external/gtest/googletest/include)
#     add_subdirectory(external/gtest/googletest)
#     if(APPLE)
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
#     endif(APPLE)
# endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
include_directories(include ${GAMMAPROP_EXTRA_INCLUDES})

add_library(gammaprop SHARED
    src/gammaprop.cpp
    #src/cosmoneutrinos.cpp
    src/KelnerAharonian2008.cpp
    src/losses.cpp
    )
    target_link_libraries(gammaprop ${GAMMAPROP_EXTRA_LIBRARIES})

# make library
#include_directories(include ${GAMMAPROP_INCLUDES})
#add_library(GAMMAPROP_LIB ${GAMMAPROP_SRCS} ${GAMMAPROP_INCLUDES})

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# copy data files in build
file(GLOB GAMMAPROP_DATA "data/*.txt")
file(COPY ${GAMMAPROP_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# make apps
add_executable(protonlosses apps/printProtonlosses.cpp)
target_link_libraries (protonlosses gammaprop ${GAMMAPROP_EXTRA_LIBRARIES})

add_executable(nuxsecs apps/printNuxsecs.cpp)
target_link_libraries (nuxsecs gammaprop ${GAMMAPROP_EXTRA_LIBRARIES})

add_executable(run apps/run.cpp)
target_link_libraries (run gammaprop ${GAMMAPROP_EXTRA_LIBRARIES})

# add_executable(generationeenergy apps/generationEnergy.cpp)
# target_link_libraries (generationeenergy beniamino ${GAMMAPROP_EXTRA_LIBRARIES})

# add_executable(uhecr apps/uhecr.cpp)
# target_link_libraries (uhecr beniamino ${GAMMAPROP_EXTRA_LIBRARIES})

# add_executable(nuflux apps/nuflux.cpp)
# target_link_libraries (nuflux beniamino ${GAMMAPROP_EXTRA_LIBRARIES})

# # testing
# if(ENABLE_TESTING)
#     enable_testing()
#     add_executable(test_common test/testCommon.cpp)
#     target_link_libraries(test_common BENIAMINO gtest gtest_main ${GAMMAPROP_EXTRA_LIBRARIES})
#     add_test(test_common test_common)

#     add_executable(test_rng test/testRng.cpp)
#     target_link_libraries(test_rng BENIAMINO gtest gtest_main ${GAMMAPROP_EXTRA_LIBRARIES})
#     add_test(test_rng test_rng)

#     add_executable(test_pid test/testPid.cpp)
#     target_link_libraries(test_pid BENIAMINO gtest gtest_main ${GAMMAPROP_EXTRA_LIBRARIES})
#     add_test(test_pid test_pid)

#     add_executable(test_cosmology test/testCosmology.cpp)
#     target_link_libraries(test_cosmology BENIAMINO gtest gtest_main ${GAMMAPROP_EXTRA_LIBRARIES})
#     add_test(test_cosmology test_cosmology)

#     add_executable(test_photonFields test/testPhotonFields.cpp)
#     target_link_libraries(test_photonFields BENIAMINO gtest gtest_main ${GAMMAPROP_EXTRA_LIBRARIES})
#     add_test(test_photonFields test_photonFields)
# endif(ENABLE_TESTING)

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "GAMMAPROP_INCLUDES ...... = ${GAMMAPROP_EXTRA_INCLUDES}")
message (STATUS "GAMMAPROP_LIBRARIES ..... = ${GAMMAPROP_EXTRA_LIBRARIES}")
