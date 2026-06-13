cmake_minimum_required(VERSION 3.12)

project(my_project VERSION 0.1
        DESCRIPTION "My Fortran program with optional MPI"
        LANGUAGES Fortran)

#ONLY ALTER ON/OFF IN THE TWO LINES BELOW (BOTH ON NOT PERMITTED)
option(USE_MPI "Enable MPI    support for the Fortran program" OFF)
option(USE_OMP "Enable OpenMP support for the Fortran program" ON)

#set(CMAKE_Fortran_STANDARD F95)
set(FVERSION "-std=gnu")
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)
if(USE_OMP)
  set(OMP "-fopenmp")
else()
  set(OMP "")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OMP} -Wall -g")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OMP} -Wall -O2")
set(CMAKE_BUILD_TYPE Debug)

if(USE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

find_package(LAPACK)
if(LAPACK_FOUND)
   set(lapackblas_libraries "${LAPACK_LIBRARIES}")
endif()

# Source code
add_subdirectory(src)

install(TARGETS my_exe DESTINATION bin)