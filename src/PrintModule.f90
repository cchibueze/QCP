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
    
end module Print_module
    