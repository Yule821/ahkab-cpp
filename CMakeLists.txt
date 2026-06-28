cmake_minimum_required(VERSION 3.16)
project(ahkab-cpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(AHKAB_BUILD_TESTS "Build test suite" ON)
option(AHKAB_BUILD_CLI "Build CLI executable" ON)

include(FetchContent)

find_package(Eigen3 REQUIRED)

if(AHKAB_BUILD_TESTS)
  find_package(GTest REQUIRED)
  enable_testing()
endif()

add_subdirectory(libahkab)

if(AHKAB_BUILD_CLI)
  add_subdirectory(src/cli)
endif()

if(AHKAB_BUILD_TESTS)
  add_subdirectory(tests)
endif()
