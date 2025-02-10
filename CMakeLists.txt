cmake_minimum_required(VERSION 3.10)
project(gnuprop CXX)
set(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(GNUPROP_EXTRA_SOURCES)
set(GNUPROP_EXTRA_INCLUDES)
set(GNUPROP_EXTRA_LIBRARIES)

set(CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -pedantic -ftree-vectorize -O3")
set(GNUPROP_INCLUDES)
set(GNUPROP_LIBRARIES)

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
list(APPEND GNUPROP_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND GNUPROP_EXTRA_LIBRARIES ${GSL_LIBRARIES})

# BENIAMINO (required)
find_package(SIMPROP REQUIRED)
list(APPEND GNUPROP_EXTRA_INCLUDES ${SIMPROP_INCLUDE_DIR})
list(APPEND GNUPROP_EXTRA_LIBRARIES ${SIMPROP_LIBRARIES})

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
include_directories(include ${GNUPROP_EXTRA_INCLUDES})

add_library(gammaprop SHARED
    src/gnuprop.cpp
    src/cached/cached.cpp
    src/cached/load.cpp
    src/interactions/GammaPairProduction.cpp
    src/interactions/InverseCompton.cpp
    src/interactions/PhotoPion.cpp
    #src/losses.cpp
    src/rates.cpp
    #src/cosmoneutrinos.cpp
    src/tridiag.cpp
    )
    target_link_libraries(gammaprop ${GNUPROP_EXTRA_LIBRARIES})

# make library
#include_directories(include ${GNUPROP_INCLUDES})
#add_library(GNUPROP_LIB ${GNUPROP_SRCS} ${GNUPROP_INCLUDES})

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# copy data files in build
file(GLOB GNUPROP_DATA "data/*.txt")
file(COPY ${GNUPROP_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# make examples
add_executable(protonlosses examples/dumpProtonLosses.cpp)
target_link_libraries (protonlosses gammaprop ${GNUPROP_EXTRA_LIBRARIES})

add_executable(photopion examples/dumpPhotoPionRates.cpp)
target_link_libraries (photopion gammaprop ${GNUPROP_EXTRA_LIBRARIES})

add_executable(gammaabsorption examples/dumpGammaAbsorption.cpp)
target_link_libraries (gammaabsorption gammaprop ${GNUPROP_EXTRA_LIBRARIES})

add_executable(pairabsorption examples/dumpPairAbsorption.cpp)
target_link_libraries (pairabsorption gammaprop ${GNUPROP_EXTRA_LIBRARIES})

# add_executable(opticaldepth examples/computeOpticalDepth.cpp)
# target_link_libraries (opticaldepth gammaprop ${GNUPROP_EXTRA_LIBRARIES})

# make apps 
add_executable(compute_tables apps/compute_tables.cpp)
target_link_libraries (compute_tables gammaprop ${GNUPROP_EXTRA_LIBRARIES})

add_executable(run apps/run.cpp)
target_link_libraries (run gammaprop ${GNUPROP_EXTRA_LIBRARIES})

# add_executable(generationeenergy apps/generationEnergy.cpp)
# target_link_libraries (generationeenergy beniamino ${GNUPROP_EXTRA_LIBRARIES})

# add_executable(uhecr apps/uhecr.cpp)
# target_link_libraries (uhecr beniamino ${GNUPROP_EXTRA_LIBRARIES})

# add_executable(nuflux apps/nuflux.cpp)
# target_link_libraries (nuflux beniamino ${GNUPROP_EXTRA_LIBRARIES})

# # testing
# if(ENABLE_TESTING)
#     enable_testing()
#     add_executable(test_common test/testCommon.cpp)
#     target_link_libraries(test_common BENIAMINO gtest gtest_main ${GNUPROP_EXTRA_LIBRARIES})
#     add_test(test_common test_common)

#     add_executable(test_rng test/testRng.cpp)
#     target_link_libraries(test_rng BENIAMINO gtest gtest_main ${GNUPROP_EXTRA_LIBRARIES})
#     add_test(test_rng test_rng)

#     add_executable(test_pid test/testPid.cpp)
#     target_link_libraries(test_pid BENIAMINO gtest gtest_main ${GNUPROP_EXTRA_LIBRARIES})
#     add_test(test_pid test_pid)

#     add_executable(test_cosmology test/testCosmology.cpp)
#     target_link_libraries(test_cosmology BENIAMINO gtest gtest_main ${GNUPROP_EXTRA_LIBRARIES})
#     add_test(test_cosmology test_cosmology)

#     add_executable(test_photonFields test/testPhotonFields.cpp)
#     target_link_libraries(test_photonFields BENIAMINO gtest gtest_main ${GNUPROP_EXTRA_LIBRARIES})
#     add_test(test_photonFields test_photonFields)
# endif(ENABLE_TESTING)

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "GNUPROP_INCLUDES ...... = ${GNUPROP_EXTRA_INCLUDES}")
message (STATUS "GNUPROP_LIBRARIES ..... = ${GNUPROP_EXTRA_LIBRARIES}")
