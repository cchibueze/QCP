!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                    INTEGRAL TRANSFORMATION                                            !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!   
module Integral_Transformation
use Print_module
!use SCF_matrix_builder 
implicit none
save
    contains
    
    subroutine two_index_transform(m,Cmat,h,h_trans)
    integer , intent(in) :: m
    real (kind = 8) , dimension(m,m), intent(in) :: Cmat
    real (kind = 8) , dimension(m,m), intent(in) :: h
    real (kind = 8) , dimension(m,m), intent(out) :: h_trans
        !call smallblockheader('2 index transformation running')
        !call writelines(2)
        h_trans = matmul( transpose(Cmat), matmul(H,Cmat) ) !H MATRIX IN MO BASIS
    end subroutine
    
    subroutine four_index_transform(m,Cmat,J_ee,J_MO)
    integer , intent(in) :: m
    real (kind = 8) , dimension(m,m), intent(in) :: Cmat
    real (kind = 8) , dimension(m**4), intent(in) :: J_ee
    real (kind = 8) , dimension(m,m,m,m) , intent(out) :: J_MO
    real (kind = 8) , dimension(m,m,m,m) :: iqrs_tmp
    real (kind = 8) , dimension(m,m,m,m) :: ijrs_tmp
    real (kind = 8) , dimension(m,m,m,m) :: ijks_tmp
    integer :: i,j,k,l,p,q,r,s
    
        iqrs_tmp = 0
        ijrs_tmp = 0
        ijks_tmp = 0
        J_MO = 0
        do s=1,m
            do r=1,m
                do q=1,m
                    do p=1,m
                        do i=1,m
                            iqrs_tmp(i,q,r,s) = iqrs_tmp(i,q,r,s) + Cmat(p,i)*J_ee(p + m*( (q-1) + m*( (r-1) + m*(s-1) ) ))
                        enddo
                    enddo
                enddo
            enddo
        enddo


        do s=1,m
            do r=1,m
                do q=1,m
                    do j=1,m
                        do i=1,m
                            ijrs_tmp(i,j,r,s) = ijrs_tmp(i,j,r,s) + Cmat(q,j)*iqrs_tmp(i,q,r,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        

        do s=1,m
            do r=1,m
                do k=1,m
                    do j=1,m
                        do i=1,m
                            ijks_tmp(i,j,k,s) = ijks_tmp(i,j,k,s) + Cmat(r,k)*ijrs_tmp(i,j,r,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    
        do s=1,m
            do l=1,m
                do k=1,m
                    do j=1,m
                        do i=1,m
                            J_MO(i,j,k,l) = J_MO(i,j,k,l) + Cmat(s,l)*ijks_tmp(i,j,k,s) ! J_TENSOR IN MO BASIS
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine
end module
    
!#######################################################################################################################!