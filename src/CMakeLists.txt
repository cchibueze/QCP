set(MY_MODS AOData.f90 AtomData.f90 CalcData.f90 CC.f90 CI.f90
            GeometryOptimizations.f90 InputProcessor.f90 Integrals.f90 IntegralTensors.f90 IntegralTransformer.f90
            MatVecOps.f90 MolData.f90 MP2.f90 NucAOsGenerator.f90 NucData.f90 NuclearDerivatives.f90 NuclearVibration.f90
            Parallel.f90 PrintModule.f90 Relativistic.f90 ResultData.f90 SCFMatBuilder.f90 SCFRoutines.f90 SpecialFunctions.f90 Timing.f90)

add_library(mylib ${MY_MODS})
add_executable(my_exe QCP.f90)
target_link_libraries(my_exe PRIVATE mylib)

target_compile_options(mylib PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
target_compile_options(my_exe PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)

target_compile_definitions(mylib PRIVATE $<$<BOOL:${USE_MPI}>:USE_MPI>)
target_compile_definitions(my_exe PRIVATE $<$<BOOL:${USE_MPI}>:USE_MPI>)

if(USE_MPI)
  target_link_libraries(mylib PUBLIC MPI::MPI_Fortran)
  target_link_libraries(my_exe PRIVATE MPI::MPI_Fortran)
endif()

if(LAPACK_FOUND)
  target_link_libraries(my_exe PRIVATE ${LAPACK_LIBRARIES})
endif()
