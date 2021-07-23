cmake_minimum_required(VERSION 3.12)
project(LxShapefilesMerge)

#set(EXECUTABLE_OUTPUT_PATH ../bin/${CMAKE_BUILD_TYPE})


find_package(fmt CONFIG REQUIRED)

find_package(Boost REQUIRED COMPONENTS filesystem graph log)
if(Boost_FOUND)
	set(Boost_USE_STATIC_LIBS ON)
	set(Boost_USE_MULTITHREAD OFF)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
else()
  message(STATUS "NOTICE: This exe requires Boost and will not be compiled.")
endif()

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
  include(${CGAL_USE_FILE})
else()
  message(STATUS "NOTICE: This exe requires CGAL and will not be compiled.")
endif()
 
find_package(GDAL REQUIRED)
if (GDAL_FOUND)
  include_directories(${GDAL_INCLUDE_DIRS})
else()
  message(STATUS "NOTICE: This exe requires GDAL and will not be compiled.")
endif()

find_package(OpenMP)
if (OpenMP_FOUND)
else()
  message(STATUS "NOTICE: This exe uses OpenMP for better performances (not found).")
endif()

include_directories(include)
include_directories("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/shared/include/")
include_directories("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/IO_DATA/include/")

add_definitions(-D_CONSOLE -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)

file(GLOB_RECURSE LX_SHPMERGE_SRCS include/* src/*)
add_executable(LxShapefilesMerge ${LX_SHPMERGE_SRCS})
add_dependencies(${PROJECT_NAME} LxGeoFactoryShared)
add_dependencies(${PROJECT_NAME} IO_DATA)

target_link_libraries(LxShapefilesMerge ${GDAL_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_LOG_LIBRARY} ${OpenCV_LIBS} fmt::fmt fmt::fmt-header-only)

target_link_libraries(${PROJECT_NAME} debug ../out/Debug/LxGeoFactoryShared)
target_link_libraries(${PROJECT_NAME} optimized ../out/Release/LxGeoFactoryShared)

target_link_libraries(${PROJECT_NAME} debug ../out/Debug/IO_DATA)
target_link_libraries(${PROJECT_NAME} optimized ../out/Release/IO_DATA)

target_compile_options(LxShapefilesMerge PRIVATE ${OpenMP_CXX_FLAGS})