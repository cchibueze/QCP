module Print_module
implicit none
save
contains
    subroutine bigblockheader(string)
    character (len=*) :: string
        write(77,*)'#######################################################################'
        write(77,*)'#                                                                     #'
        write(77,*)'#                                                                     #'
        write(77,*)'# ',string    
        write(77,*)'#                                                                     #'
        write(77,*)'#                                                                     #'
        write(77,*)'#######################################################################'
    end subroutine bigblockheader
    
    subroutine smallblockheader(string)
    character (len=*) :: string
        write(77,*)'#######################################################################'
        write(77,*)'#                                                                     #'
        write(77,*)'# ',string      
        write(77,*)'#                                                                     #'
        write(77,*)'#######################################################################'
    end subroutine smallblockheader
    
    subroutine onelineheader(string)
    character (len=*) :: string
        write(77,*) string
        write(77,*)'======================================================================='
    end subroutine onelineheader
    
    subroutine twolinesfooter(string)
    character (len=*) :: string
        write(77,*) string
        write(77,*)'======================================================================='
        write(77,*)'======================================================================='
    end subroutine twolinesfooter
    
    subroutine writelines(n)
    integer , intent(in) :: n
    integer :: i
        do i=1,n
            write(77,*)
        enddo
    end subroutine writelines

    subroutine printtext(string, int, r)
#ifdef USE_MPI
    use mpi
#endif
    real(8), intent(in), optional :: r
    integer, intent(in), optional :: int
    character (len=*), intent(in) :: string
#ifdef USE_MPI
    integer :: ierr, rank
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank  , ierr)
        if (rank == 0) then
#endif
        if    (present(r)) then
            print *, string, r
        else if(present(int)) then
            print *, string, int
        else 
            print *, string
        endif
#ifdef USE_MPI
        else
            continue
        endif
#endif
    end subroutine

end module Print_module
    