cmake_minimum_required(VERSION 3.15)
project(sat_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math")

find_package(pybind11 REQUIRED)
pybind11_add_module(sat_engine src/engine.cpp src/sgp4.cpp src/deepspace.cpp)
target_include_directories(sat_engine PRIVATE src)

# OpenMP lets the batch propagation loop use every core; optional so the build
# still succeeds on toolchains without it.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(sat_engine PRIVATE OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP enabled - parallel propagation")
else()
  message(STATUS "OpenMP not found - single-threaded propagation")
endif()

install(TARGETS sat_engine DESTINATION .)
