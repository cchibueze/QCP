cmake_minimum_required(VERSION 3.12)

project(my_project VERSION 0.1
        DESCRIPTION "My Fortran program with optional MPI"
        LANGUAGES Fortran)

option(USE_MPI "Enable MPI support for the Fortran program" OFF)

set(CMAKE_Fortran_STANDARD GNU)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -g")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -Wall -O2")

if(USE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

find_package(LAPACK)

add_subdirectory(src)

install(TARGETS my_exe DESTINATION bin)