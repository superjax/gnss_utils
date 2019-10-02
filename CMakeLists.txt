cmake_minimum_required(VERSION 2.8.3)
project(gnss_utils)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)

function(find_deps package url)
    find_package(${package} QUIET)
    if (${package}_FOUND)
        message(STATUS "${PROJECT_NAME} found package ${package}")
    elseif(TARGET ${package})
        if(EXISTS ${CMAKE_SOURCE_DIR}/lib/${package})
            message(STATUS 
                "${PROJECT_NAME} found ${package} in ${CMAKE_SOURCE_DIR}/lib/")
            set(${package}_INCLUDE_DIRS
                ${CMAKE_SOURCE_DIR}/lib/${package}/include PARENT_SCOPE)
        else()
            message(STATUS 
                "${PROJECT_NAME} could not find existing ${package} TARGET")
        endif()
    elseif(EXISTS ${CMAKE_CURRENT_LIST_DIR}/lib/${package})
        message(STATUS 
            "${PROJECT_NAME} found ${package} in ${PROJECT_NAME}/lib/")
        add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/lib/${package})
        set(${package}_INCLUDE_DIRS 
            ${CMAKE_CURRENT_LIST_DIR}/lib/${package}/include PARENT_SCOPE)
    else()
        message(STATUS "${PROJECT_NAME} did not find ${package}")
        message(STATUS "cloning ${package} into local ${PROJECT_NAME}/lib/")
        # Clone and create the proper variables
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
        execute_process(COMMAND 
            git clone ${url} ${CMAKE_CURRENT_LIST_DIR}/lib/${package})
        add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/lib/${package})
        set(${package}_INCLUDE_DIRS 
            ${CMAKE_CURRENT_LIST_DIR}/lib/${package}/include PARENT_SCOPE)
    endif()
endfunction()

find_deps(geometry https://github.com/superjax/geometry)

include_directories(gnss_utils
    ${EIGEN3_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    ${geometry_INCLUDE_DIRS}
    include)
add_library(gnss_utils
    src/gtime.cpp
    src/datetime.cpp
    src/satellite.cpp
    )
# target_include_directories(gnss_utils PUBLIC 
#     ${EIGEN3_INCLUDE_DIRS}
#     ${YAML_CPP_INCLUDE_DIR}
#     ${geometry_INCLUDE_DIRS}
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>  
#     $<INSTALL_INTERFACE:include/>)
target_link_libraries(gnss_utils
    ${YAML_CPP_LIBRARIES}
    stdc++fs)

option(BUILD_TESTS "Build Tests" OFF)
if (BUILD_TESTS)
    find_package(GTest REQUIRED)
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
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  
        $<INSTALL_INTERFACE:include/>
        ${GTEST_INCLUDE_DIRS})
endif()

option(BUILD_ROS_EXAMPLE "Build ROS Example" OFF)
if (BUILD_ROS_EXAMPLE)
    find_package(catkin REQUIRED COMPONENTS
    roscpp
    nav_msgs
    geometry_msgs
    sensor_msgs
    inertial_sense)
    message(STATUS "Building gnss_utils ROS demo")
    add_message_files(FILES SatInfo.msg)
    generate_messages(DEPENDENCIES std_msgs geometry_msgs inertial_sense)
    catkin_package(
        INCLUDE_DIRS include
        CATKIN_DEPENDS geometry_msgs nav_msgs inertial_sense sensor_msgs
        LIBRARIES gnss_utils)
    add_executable(point_position_demo
        src/positioning_demo.cpp)
    # target_include_directories(point_position_demo PUBLIC
    #     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  
    #     $<INSTALL_INTERFACE:include/>
    #     ${catkin_INCLUDE_DIRS})
    target_link_libraries(point_position_demo gnss_utils ${catkin_LIBRARIES})
    add_dependencies(point_position_demo
        gnss_utils_generate_messages_cpp)

    # install package.xml for ROS release
    install(FILES package.xml DESTINATION share/${PROJECT_NAME})
else()
    message(STATUS "Disabling ROS Example Node Build")
endif()


# install
set(LIB_DEST lib/${PROJECT_NAME})
set(INCLUDE_DEST include)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION ${LIB_DEST}
  LIBRARY DESTINATION ${LIB_DEST}
)
install(DIRECTORY include/${PROJECT_NAME}
  DESTINATION ${INCLUDE_DEST}
  FILES_MATCHING PATTERN "*.h"
)

# install CMake package configuration
install(EXPORT ${PROJECT_NAME}-targets DESTINATION ${LIB_DEST})
install(FILES ${PROJECT_NAME}-config.cmake DESTINATION ${LIB_DEST})

