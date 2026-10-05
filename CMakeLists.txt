cmake_minimum_required(VERSION 3.16)
project(claude_vrp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(claude_vrp src/main.cc src/instance.cc src/solution.cc src/construction.cc src/local_search.cc src/route_elimination.cc src/iterated_local_search.cc src/destroy_repair.cc src/adaptive_lns.cc src/fleet_minimization.cc src/ejection_search.cc src/genetic_search.cc)
target_include_directories(claude_vrp PRIVATE include)

find_package(Threads REQUIRED)
target_link_libraries(claude_vrp PRIVATE Threads::Threads)
target_compile_options(claude_vrp PRIVATE -Wall -Wextra -Wpedantic)
