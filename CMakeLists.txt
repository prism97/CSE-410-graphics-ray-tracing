cmake_minimum_required(VERSION 3.19)
project(1605094)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

set(OpenGlLinkers -lglut -lGLU -lGL)
add_executable(1605094
        1605094_Header.h
        1605094_Main.cpp)
target_link_libraries(1605094 ${OpenGlLinkers})
