cmake_minimum_required(VERSION 3.16)
project(claude_vrp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Threads REQUIRED)

# Solver library, shared by the executable and the tests.
add_library(claude_vrp_lib STATIC
    src/instance.cc src/solution.cc src/construction.cc src/local_search.cc
    src/route_elimination.cc src/iterated_local_search.cc
    src/destroy_repair.cc src/adaptive_lns.cc src/fleet_minimization.cc
    src/ejection_search.cc src/genetic_search.cc)
target_include_directories(claude_vrp_lib PUBLIC include)
target_link_libraries(claude_vrp_lib PUBLIC Threads::Threads)
target_compile_options(claude_vrp_lib PRIVATE -Wall -Wextra -Wpedantic)

add_executable(claude_vrp src/main.cc)
target_link_libraries(claude_vrp PRIVATE claude_vrp_lib)
target_compile_options(claude_vrp PRIVATE -Wall -Wextra -Wpedantic)

option(CLAUDE_VRP_BUILD_TESTS "Build the tests (run with ctest)" ON)
if(CLAUDE_VRP_BUILD_TESTS)
    enable_testing()
    add_executable(claude_vrp_tests tests/solver_tests.cc)
    target_link_libraries(claude_vrp_tests PRIVATE claude_vrp_lib)
    target_compile_options(claude_vrp_tests PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_definitions(claude_vrp_tests PRIVATE
        CLAUDE_VRP_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

    # One ctest entry per unit test (see allTests() in tests/solver_tests.cc).
    foreach(test instance segments insertion construction validate
            local_search soft_local_search iterated_local_search
            destroy_repair search_phases write_routes)
        add_test(NAME ${test} COMMAND claude_vrp_tests ${test})
    endforeach()

    # End to end: solve an instance, write the routes, check them with the
    # independent Python verifier.
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_FOUND)
        set(solution_file ${CMAKE_CURRENT_BINARY_DIR}/end_to_end_C1_10_1.sol)
        add_test(NAME end_to_end_solve
            COMMAND claude_vrp ${CMAKE_CURRENT_SOURCE_DIR}/data/C1_10_1.TXT 5
                    --threads 2 --output ${solution_file})
        set_tests_properties(end_to_end_solve PROPERTIES
            FIXTURES_SETUP end_to_end)
        add_test(NAME end_to_end_verify
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/experiments/verify_solutions.py
                    ${solution_file}
                    --data ${CMAKE_CURRENT_SOURCE_DIR}/data)
        set_tests_properties(end_to_end_verify PROPERTIES
            FIXTURES_REQUIRED end_to_end)
    endif()
endif()
