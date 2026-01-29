cmake_minimum_required(VERSION 3.25)
project(orbitsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_library(orbitsim INTERFACE)
target_include_directories(orbitsim
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/glm
)
target_compile_features(orbitsim INTERFACE cxx_std_20)

add_executable(orbitsim_example main.cpp)
target_link_libraries(orbitsim_example PRIVATE orbitsim)

enable_testing()
add_executable(orbitsim_tests tests/orbitsim_test.cpp)
target_link_libraries(orbitsim_tests PRIVATE orbitsim)
add_test(NAME orbitsim_tests COMMAND orbitsim_tests)
