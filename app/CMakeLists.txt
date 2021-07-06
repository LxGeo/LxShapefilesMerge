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

add_definitions(-D_CONSOLE -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)

file(GLOB_RECURSE LX_SHPMERGE_SRCS include/* src/*)
add_executable(LxShapefilesMerge ${LX_SHPMERGE_SRCS})

#add_dependencies(LxShapefilesMerge LxGeoFactoryShared)

target_link_libraries(LxShapefilesMerge ${GDAL_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_LOG_LIBRARY} ${OpenCV_LIBS} fmt::fmt fmt::fmt-header-only)


target_compile_options(LxShapefilesMerge PRIVATE ${OpenMP_CXX_FLAGS})