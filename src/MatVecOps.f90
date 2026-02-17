!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                           MATRIX AND VECTOR OPERATIONS                                                !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!    
module matvec_operations
implicit none
save
!=======================================================================================================================!
!                                                      MODULE SUBROUTINES                                               !
!=======================================================================================================================!

contains
    
    subroutine inv(m,A,Ainv)
        integer :: m
        real (kind=8),intent(in) :: A(m,m)
        real (kind=8),intent(out):: Ainv(size(A,1),size(A,2))
        real (kind=8)            :: work(size(A,1))            ! work array for LAPACK
        integer         :: n,info,ipiv(size(A,1))     ! pivot indices

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
    
        ! SGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call dGETRF(n,n,Ainv,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        ! SGETRI computes the inverse of a matrix using the LU factorization
        ! computed by SGETRF.
        call DGETRI(n,Ainv,n,ipiv,work,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end subroutine inv

    !#######################################################################################################################! 
    
    
    subroutine gdiag(m,A,B,w,V)
    integer , intent(in) :: m
    real (kind = 8) , dimension(1:m,1:m) , intent(in) :: A,B
    real (kind = 8) , dimension(1:m,1:m) :: C,D
    real (kind = 8) , dimension(1:m) ,  intent(out) :: w
    real (kind = 8) , dimension(1:m,1:m) , intent(out) :: V
    integer :: LWORK,nb,ilaenv
    integer , allocatable , dimension(:) :: WORK
    integer :: INFO
        C = A
        D = B
        LWORK = -1
        allocate(WORK(3*(m**2)))
        call dsygv( 2, 'V', 'U', m, C, m, D, m, w, WORK, LWORK, INFO )
        !print *, 'ab', info
        LWORK = max( 3*(m**2), int(WORK(1)) )
        call dsygv( 2, 'V', 'U', m, C, m, D, m, w, WORK, LWORK, INFO )
        IF( INFO /= 0 ) THEN
            print *, "info: ",info 
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF
        V = C
    end subroutine gdiag
    
    !#######################################################################################################################! 
    
    
    subroutine diag(m,A,w,V)
    integer , intent(in) :: m
    real (kind = 8) , dimension(1:m,1:m) , intent(in) :: A
    real (kind = 8) , dimension(1:m,1:m) :: B
    real (kind = 8) , dimension(1:m) ,  intent(out) :: w
    real (kind = 8) , dimension(1:m,1:m) , intent(out) :: V
    integer :: LWORK,nb,ilaenv
    integer , allocatable , dimension(:) :: WORK
    integer :: INFO
        B = A
        LWORK = -1
        !nb = ilaenv( 1, 'DSYTRD', 'L', m, -1, -1, -1 )
        !lwork = max((nb+2)*m, 3*m)
        !print *, lwork
        !allocate(work(lwork))
        allocate(WORK(3*(m**2)))
        call dsyev( 'V', 'U', m, B, m, w, WORK, LWORK, INFO )
        !print *, 'ab', info
        LWORK = max( 3*(m**2), int(WORK(1)) )
        call dsyev( 'V', 'U', m, B, m, w, WORK, LWORK, INFO )
        IF( INFO /= 0 ) THEN
            print *, "info: ",info 
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF
        V = B
    end subroutine diag

    !#######################################################################################################################! 
    
    
    subroutine matprint(m,A)
    integer , intent(in) :: m
    real (kind = 8) , dimension(1:m,1:m) , intent(in) :: A
    integer :: i,j
        write(77,*)''
        do i=1,m
            write(77,9837), (A(i,j) , j=1,m)
            9837 format(12f10.6)
        enddo
        write(77,*)''
    end subroutine matprint

    
    !#######################################################################################################################! 
    
    
    subroutine vecprint(m,v)
    integer , intent(in) :: m
    real (kind = 8) , dimension(1:m) :: v
    integer :: i
        write(77,*)''
        do i=1,m
            write(77,9837) v(i)
            9837 format(f15.6)
        enddo
        write(77,*)''
    end subroutine vecprint
    
 
    
    
    !#######################################################################################################################! 
    
    
    subroutine norm(d,v,n)
    integer , intent(in) :: d
    real (kind = 8) :: n2
    real (kind = 8) , intent(out) :: n
    real (kind = 8) , dimension(1:d) :: v
    integer :: i
        n2 = dot_product(v,v)
        n = sqrt(n2)
    end subroutine norm
end module matvec_operations    
!#######################################################################################################################!