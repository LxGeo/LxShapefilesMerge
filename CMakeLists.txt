cmake_minimum_required(VERSION 3.12)
project(PROJECTNAME VERSION 0.1.0)

if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++17 -O3")

include(CTest)
enable_testing()

add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/" build)
add_subdirectory(app)
add_subdirectory(lib)