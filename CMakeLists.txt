cmake_minimum_required(VERSION 2.8.3)
project(gnss_utils)

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

add_library(gnss_utils
    src/gtime.cpp
    src/datetime.cpp
    src/satellite.cpp
    )
target_include_directories(gnss_utils PUBLIC 
    include
    ${EIGEN3_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    lib/geometry/include)
target_link_libraries(gnss_utils
    ${YAML_CPP_LIBRARIES}
    stdc++fs
    geometry)

if (${GTEST_FOUND})
    message(STATUS "Found Gtest, building gnss_utils tests")
    add_definitions(-DGNSS_UTILS_DIR="${CMAKE_CURRENT_LIST_DIR}")
    add_executable(gnss_utils_test
        src/test/test_time.cpp
        src/test/test_wgs84.cpp
        src/test/test_satellite.cpp
        src/reference/reference_algorithms.cpp)
    target_link_libraries(gnss_utils_test
            ${GTEST_LIBRARIES}
            gtest_main
            gtest
            pthread
            gnss_utils)
    target_include_directories(gnss_utils_test PUBLIC
            include
            ${GTEST_INCLUDE_DIRS})
endif()

if (${catkin_FOUND})
    if(${BUILD_GNSS_UTILS_ROS})
    message(STATUS "Found Catkin, building gnss_utils ROS demo")
    add_message_files(FILES SatInfo.msg)
    generate_messages(DEPENDENCIES std_msgs geometry_msgs inertial_sense)
    catkin_package(
        INCLUDE_DIRS include
        CATKIN_DEPENDS geometry_msgs nav_msgs inertial_sense sensor_msgs
        LIBRARIES gnss_utils)
    add_executable(point_position_demo
        src/positioning_demo.cpp)
    target_include_directories(point_position_demo PUBLIC
        include
        ${catkin_INCLUDE_DIRS}
        lib/logger/include)
    target_link_libraries(point_position_demo gnss_utils ${catkin_LIBRARIES})
    add_dependencies(point_position_demo
        gnss_utils_generate_messages_cpp)
    else()
    message(STATUS "Disabling ROS Node Build")
    endif()
endif()
