!#######################################################################################################################!   
module Relativistic
use Print_module
use SCF_matrix_builder
use molecular_data
implicit none
save
private :: Massvelocity_term_calc,Darwin_term_calc
contains

    subroutine Relativistic_Correction()
    use result_data
    real (kind=8) :: MV_term,D_term
        call smallblockheader('SCALAR RELATIVISTIC CORRECTION CALCULATION')
        call writelines(2)
        call onelineheader('Starting SRC routine')
        call writelines(2)
        call Massvelocity_term_calc(Cmo,MV_term)
        call Darwin_term_calc(Cmo,D_term)
        Esrc = MV_term + D_term
        call onelineheader('Scalar Relativistic Perturbation Energy calculated!')
        write(77,*)''
        write(77,*) 'Mass-Velocity Term: ',MV_term  
        write(77,*) 'Darwin Term: ',D_term
        write(77,*) 'Total Correction: ',Esrc
        write(77,*) 'Total Scalar Relativistic HF Energy: ',EHF+Esrc
        call writelines(2)
        write(77,*)
        call twolinesfooter("SCALAR RELATIVISTIC CORRECTION CALCULATION DONE")
        call writelines(8)
    end subroutine Relativistic_Correction

    !#######################################################################################################################! 
    
    subroutine Massvelocity_term_calc(C,MV_term)
    real (kind=8) , intent(in) , dimension(aotot,aotot) :: C
    real (kind=8) , intent(out) :: MV_term
    integer :: i,p,q
    real (kind=8) , dimension(aotot,aotot) :: MV
        call MV_calc(MV)
        MV_term = 0
        do i=1,int(eltot/2)
            do p=1,aotot
                do q=1,aotot
                    MV_term = MV_term + 2*C(p,i)*c(q,i)*MV(p,q)
                enddo
            enddo
        enddo       
    end subroutine Massvelocity_term_calc

    
    !#######################################################################################################################! 
    
    subroutine Darwin_term_calc(C,D_term)
    real (kind=8) , intent(in) , dimension(aotot,aotot) :: C
    real (kind=8) , intent(out) :: D_term
    integer :: i,p,q
    real (kind=8) , dimension(aotot,aotot) :: D
        call D_calc(D)
        D_term = 0
        do i=1,int(eltot/2)
            do p=1,aotot
                do q=1,aotot
                    D_term = D_term + 2*C(p,i)*c(q,i)*D(p,q)
                enddo
            enddo
        enddo    
    end subroutine Darwin_term_calc

end module Relativistic
!#######################################################################################################################! 