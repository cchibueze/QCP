!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                 GEOMETRY OPTIMIZATIONS                                                !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
        
!#######################################################################################################################!     
module geometry_optimization
use molecular_data
use nuclei_data
use calculation_data
use Print_module
use SCF_routines
use MP2_correlation
use nuclear_derivatives
use matvec_operations
implicit none
save
contains
    
    subroutine NGC(h,Ef)
    integer :: i,a,ii,pp,ll,kk
    real (kind = 8) , intent(in) :: h
    real (kind = 8) :: E0,En,Ecur,labda,g0n,gkn,Ef
    real (kind = 8) , dimension(3*nuctot) :: g0,d0,gk_1,dk_1,gk,dk,dg,bk
    real (kind = 8) , dimension(10) :: labdas=(/1.0E0,1E-1,2E-1,3E-1,4E-1,5E-1,6E-1,7E-7,8E-1,3E-2/),E_list


    call smallblockheader('NON-LINEAR CONJUGENT GRADIENT GEOMETRY OPTIMIZATION')
    call writelines(2)

    ! STARTING NON LINEAR GRADIENT OPTIMIZATION
    !=======================================================================================================================!

    call onelineheader('Starting non-linear conjugent gradient (NGC) geometry optimization')
    
    ! CALCULATING INITIAL ELECTRONIC ENERGY, GRADIENT AND GRADIENT NORM
    !=======================================================================================================================
        call HFC(0,E0)
        if (go_wf == 'mp2') then
            call MP2_calc(0,E0)
        endif
        
        call grad_calc(h,g0)

        call norm(3*nuctot,g0,g0n)

    
    ! CHECK FOR A MINIMUM BEFORE A NGC OPTIMIZATION IS CARRIED OUT
    !=======================================================================================================================!
    
    
        if (g0n < mg_ngc) then
            call writelines(2)
            call onelineheader("Already at desired minimum")
            Ef = E0
            write(77,*) "grad:",g0n
            write(77,*) "E final:",Ef

        
    ! IN CASE OF NO INITIAL MINIMUM; CARRY ON WITH AN NGC OPTIMIZATION
    !=======================================================================================================================!
        else

        
        
    ! CARRY OUT THE FIRST LINE SEARCH BEFORE THE NGC LOOP
    !=======================================================================================================================!
            d0 = -g0
            
            ll = size(labdas)
            call writelines(2)
            write(77,*) "E initial:",E0
            write(77,*)''
            write(77,*) 'Step: 1'
            E_list = 0.0
            do i=1,ll
                labda = labdas(i)
                pp = 1
                do ii=0,nuctot-1
                    do a=1,3
                        call inf_disp(d0(pp)*labda,ii+1,a)
                        pp = pp + 1
                    enddo
                enddo
                
                call HFC(0,En)
                if (go_wf == 'mp2') then
                    call MP2_calc(0,En)
                endif
                
                E_list(i) = En
                pp = 1
                do ii=0,nuctot-1
                    do a=1,3
                        call inf_disp(-d0(pp)*labda,ii+1,a)
                        pp = pp + 1
                    enddo
                enddo
            enddo
    
            labda = labdas(minloc(E_list,dim=1))
            Ecur = E_list(minloc(E_list,dim=1))
            write(77,*) 'Current Energy:',Ecur
            pp = 1
            do ii=0,nuctot-1
                do a=1,3
                    call inf_disp(d0(pp)*labda,ii+1,a)
                    pp = pp + 1
                enddo
            enddo
    
            gk_1 = g0
            dk_1 = d0
 
        
    ! STARTING NON LINEAR GRADIENT OPTIMIZATION
    !=======================================================================================================================!
        
            do kk=1,ngc_max_iter-1
                call grad_calc(h,gk)
                call norm(3*nuctot,gk,gkn)
                write(77,*) "gradnorm:", gkn
                print *, gkn
                write(77,*)''
                if (gkn < mg_ngc) then
                    Ef = Ecur
                    write(77,*) ''
                    call onelineheader("Total Energy Converged!")
                    write(77,*)''
                    write(77,*) "NGC iteration steps:",kk
                    write(77,*) "E final:", Ef
                    exit
                    
                else
                    dg = gk - gk_1
        
                    !bk = dot_product(gk,gk)/dot_product(gk_1,gk_1)!FR
                    !bk = dot_product(gk,dg)/dot_product(gk_1,gk_1) !PR
                    bk = dot_product(gk,dg)/dot_product(dk_1,dg)  !HS
        
                    dk = -gk + bk * dk_1
        
                    E_list = 0.0
                    do i=1,ll
                        labda = labdas(i)
                        pp = 1
                        do ii=0,nuctot-1
                            do a=1,3
                                call inf_disp(dk(pp)*labda,ii+1,a)
                                pp = pp + 1
                            enddo
                        enddo
                        
                        call HFC(0,En)
                        if (go_wf == 'mp2') then
                            call MP2_calc(0,En)
                        endif
                        
                        E_list(i) = En
                        pp = 1
                        do ii=0,nuctot-1
                            do a=1,3
                                call inf_disp(-dk(pp)*labda,ii+1,a)
                                pp = pp + 1
                            enddo
                        enddo
                    enddo
    
                    labda = labdas(minloc(E_list,dim=1))         
    !               write(77,*) minloc(E_list,dim=1)
                    Ecur = E_list(minloc(E_list,dim=1))
                    write(77,*) 'Step: ', kk+1
                    write(77,*) 'Current Energy:',Ecur
                    
                    pp = 1
                    do ii=0,nuctot-1
                        do a=1,3
                            call inf_disp(dk(pp)*labda,ii+1,a)
                            pp = pp + 1
                        enddo
                    enddo
    
                    gk_1 = gk
                    dk_1 = dk
                endif
            enddo
        endif
        write(77,*)''
        write(77,*)''

    ! FINAL OPTIMIZATION GEOMETRY OUTPUT
    !=======================================================================================================================!
    
        call onelineheader('Final Geometry after NGC (Atomic Units)')
        write(77,*)''
        do i=1,nuctot
            write(77,10) nuclei(i)%name,nuclei(i)%position
    10      format(1x,A9,3f14.5)        
        
        enddo
        
        call writelines(2)
        write(77,*) 
        call twolinesfooter('NGC DONE')
        call writelines(7)
    end subroutine NGC




    !#######################################################################################################################!     
    
    

    
    subroutine Newton_Raphson(h,Ef)
    integer :: i,j,a
    real (kind=8) , dimension(3*nuctot,3*nuctot) :: Hess,Hess_inv
    real (kind=8) , dimension(3*nuctot) :: g,Hess_inv_g
    real (kind=8) , intent(in) :: h
    real (kind=8) ::E0,Ecur,Ef,gn

        call smallblockheader('NEWTON-RAPHSON GEOMETRY OPTIMIZATION')
        call writelines(2)
        call onelineheader("Starting Newton-Raphson (NR) Optimization")
        call writelines(2)
        
        call HFC(0,E0)
        write(77,*)'Initial Energy:',E0
        call writelines(2)
        write(77,*)'Step: 0'
        do i=1,NR_max_iter
            call grad_calc(h,g)
            call norm(3*nuctot,g,gn)
            write(77,*) 'gradient:',gn 
            write(77,*)''
            if (gn < mg_nr) then
                call HFC(0,Ecur)
                Ef = Ecur
                write(77,*) ''
                write(77,*) "Total Energy Converged!"
                write(77,*)'======================================================================='
                write(77,*)''
                write(77,*) "NR iteration steps:",i-1
                write(77,*) "gradnorm:",gn
                write(77,*) "E final:", Ef
                exit
            else
                write(77,*) 'Step:',i
                call Hessian_calc(h,Hess)
                call inv(3*nuctot,Hess,Hess_inv)
                Hess_inv_g = matmul(Hess_inv,g)
                do j=0,nuctot-1
                    do a=1,3
                        call inf_disp(-yNR*Hess_inv_g(3*j+a),j+1,a)
                    enddo
                enddo
            endif
        enddo
        call writelines(3)

        call onelineheader('Final Geometry after NR (Atomic Units)')
        write(77,*)''
        do i=1,nuctot
            write(77,10) nuclei(i)%name,nuclei(i)%position
    10      format(1x,A9,3f14.5)        
        enddo   
        call writelines(2)
        call twolinesfooter('Newton-Raphson Optimization DONE')
        call writelines(8)
    end subroutine Newton_Raphson
end module geometry_optimization    
!#######################################################################################################################! 