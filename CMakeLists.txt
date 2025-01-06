cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(MultiwireProportionalCounter LANGUAGES CXX)

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(Garfield REQUIRED)

set(SRC "${CMAKE_CURRENT_LIST_DIR}/src")
set(LIB "${CMAKE_CURRENT_LIST_DIR}/lib")

include_directories("${CMAKE_CURRENT_LIST_DIR}/include")

# ---Build executables-----------------------------------------------------------
add_compile_options("-std=c++17")

add_executable(generate_gas ${SRC}/generate_gas.cpp)
target_compile_options(generate_gas PRIVATE "-O2")
target_link_libraries(generate_gas Garfield::Garfield ROOT::RIO)

add_executable(test ${SRC}/test.cpp)
target_compile_options(test PRIVATE "-O2")
target_link_libraries(test Garfield::Garfield ROOT::RIO)

add_executable(view_field ${SRC}/view_field.cpp)
target_compile_options(view_field PRIVATE "-O2")
target_link_libraries(view_field Garfield::Garfield ROOT::RIO)

add_executable(mwpc ${SRC}/mwpc.cpp)
target_compile_options(mwpc PRIVATE "-fopenmp;-O3")
target_link_libraries(mwpc ${LIB}/libCustomContainers.so Garfield::Garfield ROOT::RIO ROOT::Tree)