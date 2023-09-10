cmake_minimum_required(VERSION 3.10)

project(2dv)

set(CMAKE_CXX_STANDARD 11)



add_executable(
    ${PROJECT_NAME}
    src/main.cpp
    src/parameter.cpp
    src/force.cpp
    src/division.cpp
    src/geo.cpp
    src/wangMath.cpp
    src/IOV.cpp
    src/wangSystem.cpp
    src/2dvInitial.cpp

    src/2dvOrgan.cpp
    src/2dvCell.cpp
    src/2dvLine.cpp
    src/2dvVertex.cpp
)