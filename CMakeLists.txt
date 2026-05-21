cmake_minimum_required(VERSION 3.25)
project(orbitsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

if (TARGET glm::glm)
    set(ORBITSIM_GLM_TARGET glm::glm)
elseif (TARGET glm)
    set(ORBITSIM_GLM_TARGET glm)
elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../glm/CMakeLists.txt")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../glm" "${CMAKE_CURRENT_BINARY_DIR}/glm" EXCLUDE_FROM_ALL)
    if (TARGET glm::glm)
        set(ORBITSIM_GLM_TARGET glm::glm)
    elseif (TARGET glm)
        set(ORBITSIM_GLM_TARGET glm)
    endif ()
else ()
    find_package(glm CONFIG REQUIRED)
    set(ORBITSIM_GLM_TARGET glm::glm)
endif ()

add_library(orbitsim INTERFACE
        include/orbitsim/kepler.hpp
        include/orbitsim/kepler_events.hpp
        include/orbitsim/kepler_trajectory.hpp
        include/orbitsim/kepler_maneuver.hpp
        include/orbitsim/lambert.hpp)
target_include_directories(orbitsim
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(orbitsim INTERFACE cxx_std_20)
if (ORBITSIM_GLM_TARGET)
    target_link_libraries(orbitsim INTERFACE ${ORBITSIM_GLM_TARGET})
endif ()

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
            tests/kepler_gtest.cpp
    )
    target_link_libraries(orbitsim_gtests PRIVATE orbitsim GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(orbitsim_gtests)
endif ()
