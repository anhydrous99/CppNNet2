cmake_minimum_required(VERSION 3.10)

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
find_package(ZLIB CONFIG REQUIRED)
find_package(CURL CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(SDL_mixer CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)

add_library(CppNNet2 main.cpp)