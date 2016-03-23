cmake_minimum_required(VERSION 3.4)
project(PIC__)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    complex.cpp
    complex.h
    constants.h
    correlation.cpp
    divergencecleanup.cpp
    fields.cpp
    initialization.cpp
    input.cpp
    input.h
    main.cpp
    matrix3d.cpp
    matrix3d.h
    matrixelement.cpp
    matrixElement.h
    mover.cpp
    output.cpp
    output.h
    particle.cpp
    particle.h
    random.cpp
    random.h
    simulation.cpp
    simulation.h
    specialmath.cpp
    specialmath.h
    util.cpp
    util.h
    vector3d.cpp
    vector3d.h)

add_executable(PIC__ ${SOURCE_FILES})