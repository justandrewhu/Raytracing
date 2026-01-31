cmake_minimum_required(VERSION 3.16)
project(raytracer VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optimization flags for release
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external)

# Source files
set(SOURCES
    src/main.cpp
)

# Create executable
add_executable(raytracer ${SOURCES})

# Link math library on Unix
if(UNIX)
    target_link_libraries(raytracer m)
endif()

# Copy scenes to build directory
file(COPY ${CMAKE_SOURCE_DIR}/scenes DESTINATION ${CMAKE_BINARY_DIR})

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
