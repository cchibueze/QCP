!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                            NUCLEAR VIBRATION CALCULATIONS                                             !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
        
!#######################################################################################################################!       
module nuclear_vibration
use molecular_data
use nuclei_data
use Print_module
use nuclear_derivatives
use matvec_operations
implicit none
save
private :: sqrt_mass_calc
contains
    
    subroutine sqrt_mass_calc(M_sqrt)
    integer :: i,a,l
    real (kind = 8) :: ml
    real (kind = 8) , intent(out) , dimension(3*nuctot,3*nuctot) :: M_sqrt
        M_sqrt = 0.0
        do i=0,nuctot-1
            do a=1,3
              l=3*i+a
              ml = nuclei(i+1)%mass
              M_sqrt(l,l) = sqrt(ml)
            enddo
        enddo
    end subroutine sqrt_mass_calc
    
    
    !#######################################################################################################################!
    

    subroutine  freq_analysis(h)
    integer :: i,a,c
    real (kind=8) , dimension(3*nuctot,3*nuctot) :: Hess,M_12,M_21,M_21T,M_21T_Hess,G,U,S
    real (kind=8) , dimension(3*nuctot) :: labdas,W,Ab
    real (kind=8) , intent(in) :: h
    real (kind=8) :: lab,wi,lf = 5.29177210903E-11,edtm_a,edtm2,pi=3.14159265359
        call smallblockheader('FREQUENCY CALCULATION')
        call writelines(2)
        call onelineheader('Starting frequency calculation routine')
        call writelines(3)
        
        call Hessian_calc(h,Hess)
        
        call sqrt_mass_calc(M_12)
        M_21 = 0.0
        do i =1,3*nuctot
            M_21(i,i) = 1/(M_12(i,i))
        enddo
        M_21T = transpose(M_21)
        
        M_21T_Hess = matmul(M_21T,Hess)
        G = matmul(M_21T_Hess,M_21)
        call diag(3*nuctot,G,labdas,U)
        S = matmul(M_21,U)
        
        W = 0.0
        Ab = 0.0
        do i=1,3*nuctot
            W(i) = (1/(200*pi*137)) * sqrt(abs(labdas(i))) / lf
            
            wi = sqrt(abs(labdas(i)))
            
            edtm2 = 0
            do a=1,3
                edtm_a = 0
                
                do c=0,nuctot-1
                    edtm_a = edtm_a + nuclei(c+1)%charge * S( (3*c+a),i )
                enddo
                
                edtm2 = edtm2 + edtm_a**2
            enddo
            
            Ab(i) = (6.02E23 * pi / (3 * 137)) * edtm2 * (lf**2) * 10 !L/mol/cm
        enddo
        
        
        
        call onelineheader('MODES AND FREQUENCIES (FIRST 5/6 MODES ARE TRANSLATIONAL/ROTATIONAL) ')
        write(77,*)''
        write(77,16), 'mode','frequency (cm^-1)','intensity (L/mol/cm)'
        write(77,*)''
16      format(A7,2x,A18,2x,A20)
        do i=1,3*nuctot
            write(77,17), i, W(i),Ab(i)
17          format(I5,6x,f8.2,10x,f10.5)
        enddo
        call writelines(3)
        call onelineheader('EIGENMODE VECTORS (ORDERED FROM LEFT TO RIGHT BY MODE NUMBER)')
        write(77,*)''
        call matprint(3*nuctot,U)
        call writelines(3)
        call twolinesfooter('FREQUENCY CALCULATION DONE')
        call writelines(8)
    end subroutine freq_analysis
end module nuclear_vibration
!#######################################################################################################################!  