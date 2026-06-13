!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                              PARALLELIZATION MODULE                                                   !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
!#######################################################################################################################!   
module Parallel
implicit none
save
private
public :: set_omp_variables, omp, num_threads

!=======================================================================================================================!
!                                                    MODULE VARIABLES DELCARATION                                       !
!=======================================================================================================================!
! openmp variables
logical :: omp
integer :: num_threads

!#######################################################################################################################!  
    contains


    subroutine set_omp_variables(num_threads_in)
        use Print_module
        !$ use omp_lib
        integer, intent(in) :: num_threads_in

            omp = .false.
            num_threads = num_threads_in
            !$ omp = .true.
            if (omp) then
                call printtext('OpenMP parallelization is enabled!!')
            else
                call printtext('OpenMP parallelization is disabled!!')
            endif
            if (omp) then
                !$ call omp_set_num_threads(num_threads)
                !$ print *, 'Number of threads to be used for parallelization:', omp_get_num_threads()
            endif
            
    end subroutine set_omp_variables

end module parallel
!#######################################################################################################################!type ao                        !atomic orbital type declaration
