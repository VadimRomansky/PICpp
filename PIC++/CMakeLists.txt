cmake_minimum_required(VERSION 3.4)
project(PIC__)
set(MPICH2_LIB_PATH "-L/home/vadim/mpich2-install/lib")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_COMPILER "/home/vadim/mpich2-install/bin/mpicxx")
set(CMAKE_C_COMPILER "/home/vadim/mpich2-install/bin/mpicc")

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
    vector3d.h mpi_util.cpp mpi_util.h dichotomousSolver.h dichotomousSolver.cpp)

link_libraries(${MPICH2_LIB_PATH} -lmpe)

add_executable(PIC__ ${SOURCE_FILES})