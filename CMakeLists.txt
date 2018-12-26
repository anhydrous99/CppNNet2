cmake_minimum_required(VERSION 3.10)

if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif ()

include("cmake/HunterGate.cmake")
HunterGate(
        URL "https://github.com/ruslo/hunter/archive/v0.23.83.tar.gz"
        SHA1 "12dec078717539eb7b03e6d2a17797cba9be9ba9"
)

project(CppNNet2)
set(CMAKE_CXX_STANDARD 14)

hunter_add_package(ZLIB)
hunter_add_package(CURL)
hunter_add_package(nlohmann_json)
hunter_add_package(SDL_mixer)
hunter_add_package(GTest)
hunter_add_package(OpenBLAS)
find_package(ZLIB CONFIG REQUIRED)
find_package(CURL CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(SDL_mixer CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)
find_package(OpenBLAS CONFIG REQUIRED)

include_directories(include)

add_library(CppNNet2 INTERFACE)

enable_testing()

add_executable(Matrix_Test tests/Matrix_Tests.cpp)

target_link_libraries(Matrix_Test CppNNet2 GTest::main OpenBLAS::OpenBLAS)