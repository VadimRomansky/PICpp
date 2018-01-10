cmake_minimum_required(VERSION 3.4)
project(PIC__)
set(MPICH2_LIB_PATH "-L/home/vadim/mpich2-install/lib")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_COMPILER "/home/vadim/mpich2-install/bin/mpicxx")
set(CMAKE_C_COMPILER "/home/vadim/mpich2-install/bin/mpicc")

set(SOURCE_FILES
        boundaryFieldEvaluator.h
        boundaryFieldEvaluator.cpp
    complex.cpp
    complex.h
    constants.h
    correlation.cpp
        correlationMap.cpp
        dichotomousSolver.h
        dichotomousSolver.cpp
    divergencecleanup.cpp
        divergencecleanup_buneman.cpp
    fields.cpp
        fields_exchange.cpp
        fields_filtering.cpp
        fourier.h
        fourier.cpp
    initialization.cpp
    input.cpp
    input.h
        largeVectorBasis.h
        largeVectorBasis.cpp
    main.cpp
        massMatrix.h
    matrix3d.cpp
    matrix3d.h
    matrixelement.cpp
    matrixElement.h
        memory_debug.h
        memory_util.h
        memory_util.cpp
        momentum_collecting.cpp
        momentum_exchange.cpp
    mover.cpp
    output.cpp
    output.h
    particle.cpp
    particle.h
        paths.h
    random.cpp
    random.h
        rightPartEvaluator.h
        rightPartEvaluator.cpp
    simulation.cpp
    simulation.h
        smoothing.cpp
    specialmath.cpp
    specialmath.h
    util.cpp
    util.h
    vector3d.cpp
    vector3d.h
        mpi_util.cpp
        mpi_util.h)

link_libraries(${MPICH2_LIB_PATH} -lmpe)

add_executable(PIC__ ${SOURCE_FILES})