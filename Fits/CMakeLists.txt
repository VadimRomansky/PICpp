cmake_minimum_required(VERSION 3.8)
project(Fits)
set(FITS_LIB_PATH "-L/home/vadim/cfitsio/lib/")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)

link_libraries(${FITS_LIB_PATH} -lcfitsio)
#target_link_libraries(Fits libcfitsio.a -lcfitsio)
add_executable(Fits ${SOURCE_FILES})