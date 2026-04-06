cmake_minimum_required(VERSION 3.12)

project(my_project VERSION 0.1
        DESCRIPTION "My Fortran program"
        LANGUAGES Fortran)
enable_language(Fortran C)

# Currently setting the Fortran compiler to use -std=gnu, change this if you
# want a specific standard
set(FVERSION "-std=gnu")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -g")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -Wall -march=native -O3")
set(CMAKE_BUILD_TYPE Release)

find_package(LAPACK)
if(LAPACK_FOUND)
   set(lapackblas_libraries "${LAPACK_LIBRARIES}")
endif()

# Source code
add_subdirectory(src)

install(TARGETS my_exe DESTINATION "bin")
