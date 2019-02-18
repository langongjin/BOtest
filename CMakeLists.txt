cmake_minimum_required(VERSION 3.8)
project(BOtest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(
        Boost REQUIRED COMPONENTS
        system
        filesystem
        thread
        regex
        unit_test_framework)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    find_library(Boost boost PATHS /opt/local/lib)
    include_directories(${Boost_LIBRARY_PATH})
endif()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(NLOpt REQUIRED nlopt>=2.4)
include_directories(${NLOpt_INCLUDE_DIRS})
link_directories(${NLOpt_LIBRARY_DIRS})

include_directories(
        limbo/
        external/)

#message("Hello ${NLOpt_INCLUDE_DIRS}")
set(SOURCE_FILES main.cpp)

add_executable(BOtest ${SOURCE_FILES})
target_compile_definitions(BOtest PUBLIC USE_NLOPT=1)
target_link_libraries(
        BOtest
#        ${TBB_LIBRARIES}
        ${Boost_LIBRARIES}
        ${NLOpt_LIBRARIES})
