cmake_minimum_required(VERSION 3.12.4)

if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif ()

option(USE_OPENBLAS "Use OpenBLAS" ON)
option(BUILD_TESTS "Build test programs" ON)

include("cmake/HunterGate.cmake")
HunterGate(
        URL "https://github.com/ruslo/hunter/archive/v0.23.147.tar.gz"
        SHA1 "feada7fa5f9220b4ff35ac111ddb637cafa13b05"
)

project(CppNNet2)
set(CMAKE_CXX_STANDARD 14)

if (BUILD_TESTS)
    hunter_add_package(GTest)
    find_package(GTest CONFIG REQUIRED)
endif (BUILD_TESTS)

if (USE_OPENBLAS)
    hunter_add_package(OpenBLAS)
    find_package(OpenBLAS CONFIG REQUIRED)
    add_compile_definitions(USES_OPENBLAS)
endif (USE_OPENBLAS)

include_directories(include)
add_library(CppNNet2 INTERFACE)

if (BUILD_TESTS)
    enable_testing()
    add_executable(Tests tests/Tests.cpp)
    if (USE_OPENBLAS)
        target_link_libraries(Tests OpenBLAS::OpenBLAS GTest::main)
    else (USE_OPENBLAS)
        target_link_libraries(Tests GTest::main)
    endif (USE_OPENBLAS)
endif (BUILD_TESTS)