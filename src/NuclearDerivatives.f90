!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                  NUCLEAR DERIVATIVES                                                  !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!       
module nuclear_derivatives
use calculation_data
use molecular_data
use SCF_routines
use MP2_correlation
implicit none
save
   
contains  
    
    subroutine inf_disp(dx,i,a)
    use nuclei_data
    use ao_data
    integer :: j
    integer , intent(in) :: i, a
    real (kind = 8) , intent(in) :: dx
        nuclei(i)%position(a) = nuclei(i)%position(a) + dx
        do j=1,aotot
            if (aos(j)%nucno == i) then
                aos(j)%position = nuclei(i)%position
            endif
        enddo
    end subroutine inf_disp
    

    
    
    
    !#######################################################################################################################! 
    
    
    
    subroutine grad_calc(dx,g)
    real (kind = 8) , intent(in) :: dx
    real (kind = 8) :: Ei,Ef,dE,gia
    real (kind = 8) , intent(out) , dimension(3*nuctot) :: g
    integer :: i,a
        do i=0,nuctot-1
            do a=1,3
                call inf_disp(dx,i+1,a) !!!NEW
                call HFC(0,Ef)
                if (go_wf == 'mp2') then
                    call MP2_calc(0,Ef)
                endif
                
                call inf_disp(-2*dx,i+1,a) !!!NEW
                
                call HFC(0,Ei) !!!NEW
                if (go_wf == 'mp2') then
                    call MP2_calc(0,Ei)
                endif
                
                dE = Ef - Ei
                !print *, dE
                
                gia = dE/ (2*dx) !!!NEW
                
                
                g(3*i+a)=gia
                call inf_disp(dx,i+1,a)!!!NEW
            enddo
        enddo
        end subroutine grad_calc

    
    
    !#######################################################################################################################!  

    
    subroutine Hessian_calc(dx,Hess)
    integer :: i1,i2,j1,j2,i,j
    real (kind=8) , intent(out) , dimension(3*nuctot,3*nuctot) :: Hess
    real (kind=8) , dimension(3*nuctot,3*nuctot) :: HessT
    real (kind=8) , intent(in) :: dx
    real (kind=8) ::E11,E00,E10,E01,d2Epot,d2Epot_dxidxj
        !call HFC(E0,0)
        Hess = 0.0
        do i1=0,nuctot-1
            do i2=1,3
                do j1=0,nuctot-1
                    do j2=1,3
                        i = 3*i1+i2
                        j = 3*j1+j2
                    
                        call inf_disp(dx,i1+1,i2)
                        call inf_disp(dx,j1+1,j2)
                        call HFC(0,E11)
                        if (go_wf == 'mp2') then
                            call MP2_calc(0,E11)
                        endif
                
                        call inf_disp(-2*dx,i1+1,i2)
                        call HFC(0,E01)
                        if (go_wf == 'mp2') then
                            call MP2_calc(0,E01)
                        endif
                
                        call inf_disp(-2*dx,j1+1,j2)
                        call HFC(0,E00)
                        if (go_wf == 'mp2') then
                            call MP2_calc(0,E00)
                        endif      
                
                        call inf_disp(2*dx,i1+1,i2)
                        call HFC(0,E10)
                        if (go_wf == 'mp2') then
                            call MP2_calc(0,E10)
                        endif     
                
                        call inf_disp(-dx,i1+1,i2)
                        call inf_disp(dx,j1+1,j2)
                        
                        d2Epot = E11 + E00 - E01 - E10
                        d2Epot_dxidxj = d2Epot/ (4*(dx**2))           

                        Hess(i,j) = d2Epot_dxidxj
                    enddo
                enddo
            enddo
        enddo
        HessT = transpose(Hess)
        Hess = (Hess + HessT)/2
    end subroutine Hessian_calc 
    
    
    
    
    !#######################################################################################################################! 
    
    
    
    subroutine grad_calc_forwarddiff(dx,g)
    real (kind = 8) , intent(in) :: dx
    real (kind = 8) :: Ei,Ef,dE,gia
    real (kind = 8) , intent(out) , dimension(3*nuctot) :: g
    integer :: i,a
        call HFC(0,Ei)
        do i=0,nuctot-1
            do a=1,3
                call inf_disp(dx,i+1,a) 
                call HFC(0,Ef)
                
                dE = Ef - Ei
                
                gia = dE/ dx
                
                g(3*i+a)=gia
                call inf_disp(-dx,i+1,a)
            enddo
        enddo
        end subroutine grad_calc_forwarddiff

    
    
    !#######################################################################################################################!  

    
    subroutine Hessian_calc_forwarddiff(dx,Hess)
    integer :: i1,i2,j1,j2,i,j
    real (kind=8) , intent(out) , dimension(3*nuctot,3*nuctot) :: Hess
    real (kind=8) , dimension(3*nuctot,3*nuctot) :: HessT
    real (kind=8) , intent(in) :: dx
    real (kind=8) ::E0,Ei,Ej,Eij,d2Epot,d2Epot_dxidxj
        call HFC(0,E0)
        Hess = 0.0
        do i1=0,nuctot-1
            do i2=1,3
                do j1=0,nuctot-1
                    do j2=1,3
                        i = 3*i1+i2
                        j = 3*j1+j2
                    
                        call inf_disp(dx,i1+1,i2)
                        call HFC(0,Ei)
                    
                        call inf_disp(dx,j1+1,j2)
                        call HFC(0,Eij)

                        call inf_disp(-dx,i1+1,i2)
                        call HFC(0,Ej)
                    
                        call inf_disp(-dx,j1+1,j2)
                    
                        d2Epot = E0 + Eij - Ei - Ej
                        d2Epot_dxidxj = d2Epot/(dx**2)            

                        Hess(i,j) = d2Epot_dxidxj
                    enddo
                enddo
            enddo
        enddo
        HessT = transpose(Hess)
        Hess = (Hess + HessT)/2
    end subroutine Hessian_calc_forwarddiff
    
end module nuclear_derivatives
!#######################################################################################################################!