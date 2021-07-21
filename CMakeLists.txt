cmake_minimum_required(VERSION 3.19)
project(1605094)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

set(OpenGlLinkers -lglut -lGLU -lGL)
add_executable(1605094
        ray_tracing_header.h
        ray_tracing_main.cpp)
target_link_libraries(1605094 ${OpenGlLinkers})
