cmake_minimum_required(VERSION 3.9.6)
project(photon)

file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
file(GLOB_RECURSE HEADER_FILES src/*.hpp)

message("Source files: ${SOURCE_FILES}")
message("Header files: ${HEADER_FILES}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}                             \
        -DDEBUG                                                     \
        -g                                                          \
        -O0                                                         \
    ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
        -Weverything                    \
    ")

add_executable(photon ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(photon PUBLIC ${CMAKE_SOURCE_DIR})
