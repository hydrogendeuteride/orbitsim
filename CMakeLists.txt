cmake_minimum_required(VERSION 3.25)
project(orbitsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_library(orbitsim INTERFACE
        include/orbitsim/lambert.hpp)
target_include_directories(orbitsim
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/glm
)
target_compile_features(orbitsim INTERFACE cxx_std_20)

add_executable(orbitsim_example main.cpp)
target_link_libraries(orbitsim_example PRIVATE orbitsim)

option(ORBITSIM_ENABLE_GTEST "Build GoogleTest-based tests (in addition to assert-based tests)" OFF)

if (ORBITSIM_ENABLE_GTEST)
    include(FetchContent)

    # Prefer vendored googletest (no network).
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest/CMakeLists.txt")
        add_subdirectory(third_party/googletest)
    else()
        # Otherwise try system package, then FetchContent (requires network at configure time).
        find_package(GTest CONFIG QUIET)
        if (NOT GTest_FOUND)
            FetchContent_Declare(
                    googletest
                    GIT_REPOSITORY https://github.com/google/googletest.git
                    GIT_TAG v1.14.0
            )
            set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
            FetchContent_MakeAvailable(googletest)
        endif ()
    endif ()

    add_executable(orbitsim_gtests
            tests/math_gtest.cpp
            tests/game_simulation_gtest.cpp
            tests/spacecraft_state_cache_gtest.cpp
            tests/trajectories_gtest.cpp
            tests/soi_rails_gtest.cpp
            tests/frames_gtest.cpp
            tests/maneuvers_gtest.cpp
            tests/geodesy_gtest.cpp
            tests/nodes_gtest.cpp
            tests/coordinate_frames_gtest.cpp
            tests/lambert_gtest.cpp
    )
    target_link_libraries(orbitsim_gtests PRIVATE orbitsim GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(orbitsim_gtests)
endif ()
