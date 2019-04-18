cmake_minimum_required(VERSION 2.8.3)
project(GNSS_utils)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(GTest)
find_package(yaml-cpp REQUIRED)
find_package(catkin COMPONENTS
	roscpp
	nav_msgs
	geometry_msgs
	sensor_msgs
	inertial_sense)


if (NOT TARGET geometry)
	add_subdirectory(lib/geometry)
endif()

add_library(GNSS_utils
	src/gtime.cpp
	src/datetime.cpp
	src/satellite.cpp
	)
target_include_directories(GNSS_utils PUBLIC 
	include
	${EIGEN3_INCLUDE_DIRS}
	${YAML_CPP_INCLUDE_DIR}
	lib/geometry/include)
target_link_libraries(GNSS_utils
	${YAML_CPP_LIBRARIES}
	stdc++fs
	geometry)

if (${catkin_FOUND})
	message(STATUS "Found Catkin, building GNSS_utils ROS demo")
	add_executable(point_position_demo
		src/positioning_demo.cpp)
	target_include_directories(point_position_demo PUBLIC include ${catkin_INCLUDE_DIRS})
	target_link_libraries(point_position_demo GNSS_utils)
	catkin_package(
		INCLUDE_DIRS include
		LIBRARIES GNSS_utils)
endif()

if (${GTEST_FOUND})
	message(STATUS "Found Gtest, building GNSS_utils tests")
	add_definitions(-DGNSS_UTILS_DIR="${CMAKE_CURRENT_LIST_DIR}")
	add_executable(GNSS_utils_test
	    src/test/test_time.cpp
	    src/test/test_wgs84.cpp
	    src/test/test_satellite.cpp
	    src/reference/reference_algorithms.cpp)
	target_link_libraries(GNSS_utils_test 
		${GTEST_LIBRARIES}
		gtest_main
		gtest
		pthread
		GNSS_utils)
	target_include_directories(GNSS_utils_test PUBLIC
		include
		${GTEST_INCLUDE_DIRS})
endif()
