!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                     CC CALCULATION                                                    !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!   
module CC_correlation
use Print_module
use matvec_operations
use molecular_data
use integral_tensors
use result_data
use Integral_Transformation
use Timing

implicit none
save
contains


    !#######################################################################################################################! 
    subroutine closed_shell_CCD_calc(output)
    integer , intent(in) :: output
    
    integer :: i,j,k,l,a,b,c,d,ccc,O,V
    real (kind = 8) :: ei,ej,ea,eb,A2,B2,C2,D2,E2,dt2_amps_aibj
    real (kind = 8) :: ts,te,ti,tf
    
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO
    real (kind = 8) , allocatable, dimension(:,:,:,:) :: t2_amps
    real (kind = 8) , allocatable, dimension(:,:) :: E2oo, E2vv 
    real (kind = 8) , allocatable, dimension(:,:,:,:) :: Uvovo, B2oooo, C2vovo, D2vovo
    
        call timer(ts)

        O = int(eltot/2) 
        V = aotot - O 
        if (output == 1) then
            call smallblockheader('CLOSED SHELL CCSD CORRELATION CALCULATION')
            call writelines(2)
            call onelineheader('Starting CCSD routine')
            call writelines(2)
        endif
        
  
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        allocate(J_MO(aotot,aotot,aotot,aotot))
        J_MO = 0
        !N^5
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                       GENERATING AUXILLIARY INTEGRAL TENSORS                                          !
    !=======================================================================================================================!
        

        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                             OBTAINING CCSD CLUSTER AMPLITUDES                                         !
    !=======================================================================================================================!  

        allocate(t2_amps(V,O,V,O))
        t2_amps = 0
        
        !!!! GENERATE FIRST ORDER (MP2) AMPLITUDES (t0 = 0) !!!!
        
        !O^2 * V^2
        do a=1,V
            ea = Eig(O+a)
            do i=1,O
                ei = Eig(i)
                do b=1,V
                    eb = Eig(O+b)
                    do j=1,O
                        ej = Eig(j)
                        t2_amps(a,i,b,j) = J_MO(O+a,i,O+b,j)/(ei+ej-ea-eb)
                    enddo
                enddo
            enddo
        enddo

        !!!! GENERATE CCD AMPLITUDES VIA PERTUBATION-BASED QUASI-NEWTON METHOD
        
        do ccc=1,15
            !!!! CALCULATING AUXILLIARY TENSORS
            allocate(Uvovo(V,O,V,O),B2oooo(O,O,O,O),C2vovo(V,O,V,O))
            allocate(D2vovo(V,O,V,O),E2oo(O,O),E2vv(V,V))
            Uvovo = 0
            B2oooo = 0
            C2vovo = 0
            D2vovo = 0
            E2oo = 0
            E2vv = 0
            
            !Uvovo O^2 * V^2
            do a=1,V
                do i=1,O
                    do b=1,V
                        do j=1,O
                            Uvovo(a,i,b,j) = 2* t2_amps(a,i,b,j) - t2_amps(a,j,b,i)
                        enddo
                    enddo
                enddo
            enddo
            !B2oooo O^4 * V^2
            do j=1,O
                do i=1,O
                    do l=1,O
                        do k=1,O
                            B2oooo(k,l,i,j) = B2oooo(k,l,i,j) + J_MO(k,i,l,j)
                            do d=1,V
                                do c=1,V
                                    B2oooo(k,l,i,j) = B2oooo(k,l,i,j) + t2_amps(c,i,d,j) * J_MO(k,O+c,l,O+d)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            !C2vovo and D2vovo O^3 * V^3
            do i=1,O
                do a=1,V
                    do k=1,O
                        do c=1,V
                            C2vovo(c,k,a,i) = C2vovo(c,k,a,i) + J_MO(k,i,O+a,O+c)
                            D2vovo(c,k,a,i) = D2vovo(c,k,a,i) + (2*J_MO(O+a,i,k,O+c) - J_MO(O+a,O+c,k,i))
                            do l=1,O
                                do d=1,V
                                    C2vovo(c,k,a,i) = C2vovo(c,k,a,i) - t2_amps(a,l,d,i) * J_MO(k,O+d,l,O+c) / 2
                                    D2vovo(c,k,a,i) = D2vovo(c,k,a,i) + Uvovo(a,i,d,l) * (2*J_MO(l,O+d,k,O+c) - J_MO(l,O+c,k,O+d)) / 2
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            !E2oo O^4 * V
            do j=1,O
                E2oo(j,j) = E2oo(j,j) + Eig(j)
            enddo
            do j=1,O
                do k=1,O
                    do l=1,O
                        do c=1,V
                            do d=1,V
                                E2oo(k,j) = E2oo(k,j) + Uvovo(c,l,d,j) * J_MO(k,O+d,l,O+c)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            !E2vv O * V^4
            do b=1,V
                E2vv(b,b) = E2vv(b,b) + Eig(O+b)
            enddo
            do b=1,V
                do c=1,V
                    do k=1,O
                        do l=1,O
                            do d=1,V
                                E2vv(c,b) = E2vv(c,b) - Uvovo(b,k,d,l) * J_MO(l,O+d,k,O+c)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
        
            !!!AMPLITUDE INCREMENTS!!!!!
            do a=1,V
                ea = Eig(O+a)
                do i=1,O
                    ei = Eig(i)
                    do b=1,V
                        eb = Eig(O+b)
                        do j=1,O
                            ej = Eig(j)
                            
                            !!!!ADDING T2 AMPLITUDE INCREMENTS!!!!
                            A2=0
                            B2=0
                            C2=0
                            D2=0
                            E2=0
                            
                            !O^2 * V^4
                            do d=1,V
                                do c=1,V
                                    A2 = A2 + t2_amps(c,i,d,j)*J_MO(O+a,O+c,O+b,O+d)
                                enddo
                            enddo
                            !O^4 * V^2
                            do l=1,O
                                do k=1,O
                                    B2 = B2 + t2_amps(a,k,b,l) * B2oooo(k,l,i,j)
                                enddo
                            enddo
                            !O^3 * V^3
                            do k=1,O
                                do c=1,V
                                    C2 = C2 + t2_amps(b,k,c,j) * C2vovo(c,k,a,i) + t2_amps(a,k,c,i) * C2vovo(c,k,b,j)
                                enddo
                            enddo
                            C2 = C2 * -0.5d0
                            do k=1,O
                                do c=1,V
                                    C2 = C2 - t2_amps(b,k,c,i) * C2vovo(c,k,a,j) - t2_amps(a,k,c,j) * C2vovo(c,k,b,i)
                                    D2 = D2 +   Uvovo(b,j,c,k) * D2vovo(c,k,a,i) +   Uvovo(a,i,c,k) * D2vovo(c,k,b,j)
                                enddo
                            enddo
                            D2 = D2 * 0.5d0
                            !O^3 * V^2
                            do k=1,O
                                E2 = E2 + t2_amps(a,i,b,k) * E2oo(k,j) + t2_amps(b,j,a,k) * E2oo(k,i)
                            enddo
                            E2 = -E2
                            !O^2 * V^3
                            do c=1,V
                                E2 = E2 + t2_amps(a,i,c,j) * E2vv(c,b) + t2_amps(b,j,c,i) * E2vv(c,a)
                            enddo
                            !O^2 * V^2
                            dt2_amps_aibj = ( J_MO(O+a,i,O+b,j) + (A2 + B2 + C2 + D2 + E2) )/(ei + ej - ea - eb)
                            
                            !print *, dt2_amps_aibj
                            t2_amps(a,i,b,j) = t2_amps(a,i,b,j) + dt2_amps_aibj
                        enddo
                    enddo
                enddo
            enddo
            
            !deallocate auxilliary tensors!!!
            deallocate(Uvovo,B2oooo,C2vovo,D2vovo,E2oo,E2vv)
        !'ccc loop cycle'    
        enddo      
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                                     OBTAINING CCSD ENERGY                                             !
    !=======================================================================================================================! 
        !O^2 * V^2
        Eccd = 0
        do a=1,V
            do i=1,O
                do b=1,V
                    do j=1,O
                        Eccd = Eccd + t2_amps(a,i,b,j) * (2*J_MO(i,O+a,j,O+b) - J_MO(i,O+b,j,O+a) )
                    enddo
                enddo
            enddo
        enddo
        
        
    !#######################################################################################################################!
        
        
        if (output == 1) then
            call onelineheader("Correlation Energy Calculated!")
            write(77,*)''
            write(77,*) 'Hatree-Fock Energy: ',EHF
            write(77,*) 'CCD Energy: ',EHF+Eccd
            write(77,*) 'CCD Correlation Energy: ',Eccd
            call writelines(2)
            call twolinesfooter("CCD CALCULATION DONE")
            call writelines(8)
        endif

        call timer(te)
        print *, 'Time taken for CCD calculation: ', te-ts, ' seconds'
        
    end subroutine closed_shell_CCD_calc
    !#######################################################################################################################!
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !#######################################################################################################################! 
    subroutine closed_shell_CCSD_calc(output)
    integer , intent(in) :: output
    
    integer :: i,j,k,l,a,b,c,d,p,q,r,s,calcF,ccc,O,V,trs=0
    integer , dimension(aotot) :: mo_indices
    
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: iqrs_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ijrs_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ijks_tmp
    
    real (kind = 8) :: ei,ej,ea,eb,P1,P2,A2,B2,C2,D2,E2,dt1_amps_ai,dt2_amps_aibj
    real (kind = 8) , dimension(aotot,aotot) :: H_MO,FI_MO,H_T1,FI_T1,t1,x,y
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO,J_T1!,iajb_array2
    
    
    real (kind = 8) , allocatable, dimension(:,:) :: t1_amps
    real (kind = 8) , allocatable, dimension(:,:,:,:) :: t2_amps
    real (kind = 8) , allocatable, dimension(:,:) :: E2oo, E2vv 
    real (kind = 8) , allocatable, dimension(:,:,:,:) :: Uvovo, B2oooo, C2vovo, D2vovo
    
        O = int(eltot/2) 
        V = aotot - O 
        if (output == 1) then
            call smallblockheader('CLOSED SHELL CCSD CORRELATION CALCULATION')
            call writelines(2)
            call onelineheader('Starting CCSD routine')
            call writelines(2)
        endif
        
  
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        H_MO = matmul( transpose(Cmo), matmul(H,Cmo) ) !H MATRIX IN MO BASIS
        allocate(J_MO(aotot,aotot,aotot,aotot),J_T1(aotot,aotot,aotot,aotot))
        
        H_MO = 0
        J_MO = 0
        !N^3
        call two_index_transform(aotot,Cmo,H,H_MO)         !H_MO = matmul( transpose(Cmo), matmul(H,Cmo) ) !H MATRIX IN MO BASIS
        !N^5
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
        
        !N^2 * O
        FI_MO = 0
        do p=1,aotot
            do q=1,aotot
                FI_MO(p,q) = H_MO(p,q)
                do i=1,O
                    FI_MO(p,q) = FI_MO(p,q) + 2*J_MO(p,q,i,i) - J_MO(p,i,i,q)
                enddo
            enddo
        enddo
        

        
    !#######################################################################################################################!   
        
        
        
        
        
        
        
        
        
        
        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                             OBTAINING CCSD CLUSTER AMPLITUDES                                         !
    !=======================================================================================================================!  

        allocate(t1_amps(V,O))
        allocate(t2_amps(V,O,V,O))
        
        t1_amps = 0
        t2_amps = 0
        !!!! GENERATE FIRST ORDER (MP2) AMPLITUDES (t0 = 0) !!!!
        !!!! NOTE THERE IS NO t=0 SINGLES FIRST ORDER CONTRIBUTION DUE TO A DIAGONAL FOCK MATRIX IN THE MO BASIS!!!!
        
        !O^2 * V^2
        do a=1,V
            ea = Eig(O+a)
            do i=1,O
                ei = Eig(i)
                do b=1,V
                    eb = Eig(O+b)
                    do j=1,O
                        ej = Eig(j)
                        t2_amps(a,i,b,j) = J_MO(O+a,i,O+b,j)/(ei+ej-ea-eb)
                    enddo
                enddo
            enddo
        enddo
        

        !!!! GENERATE CCD AMPLITUDES VIA PERTUBATION-BASED QUASI-NEWTON METHOD
        
        do ccc=1,30
            !print *, 'kk'
            
            !!!! INVOKING T1 TRANSFORMATION !!!!        
            if (ccc == 1) then
                H_T1 = H_MO
                J_T1 = J_MO
                FI_T1 = FI_MO
            else

                !!!! CREATE T1, X AND Y MATRICES !!!!
                t1 = 0
                do a=1,V
                    do i=1,O
                        t1(a+O,i) = t1_amps(a,i)
                    enddo
                enddo
                

                
                x=0
                y=0
                
                do i=1,aotot
                    x(i,i) = 1
                    y(i,i) = 1
                enddo
                
                x = x - t1
                y = y + transpose(t1)
                
                !print *, t1

            !=======================================================================================================================!
            !                                          INTEGRAL TRANSFORMATION (MO ==> T1)                                          !
            !=======================================================================================================================!
                !N^3
                H_T1 = matmul( matmul(x,H_MO),transpose(y) ) !H MATRIX IN T1 BASIS
                
                !N^5
                allocate(iqrs_tmp(aotot,aotot,aotot,aotot),ijrs_tmp(aotot,aotot,aotot,aotot),ijks_tmp(aotot,aotot,aotot,aotot))
                iqrs_tmp = 0 ! MAKE J TENSOR IN T1 BASIS, ijkl are general indices here (unlike convinient occupied incices!!!)
                ijrs_tmp = 0
                ijks_tmp = 0
                J_T1 = 0
                do s=1,aotot
                    do r=1,aotot
                        do q=1,aotot
                            do p=1,aotot
                                do i=1,aotot
                                    iqrs_tmp(i,q,r,s) = iqrs_tmp(i,q,r,s) + x(i,p)*J_MO(p,q,r,s)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo


                do s=1,aotot
                    do r=1,aotot
                        do q=1,aotot
                            do j=1,aotot
                                do i=1,aotot
                                    ijrs_tmp(i,j,r,s) = ijrs_tmp(i,j,r,s) + y(j,q)*iqrs_tmp(i,q,r,s)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
        

                do s=1,aotot
                    do r=1,aotot
                        do k=1,aotot
                            do j=1,aotot
                                do i=1,aotot
                                    ijks_tmp(i,j,k,s) = ijks_tmp(i,j,k,s) + x(k,r)*ijrs_tmp(i,j,r,s)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
    
                do s=1,aotot
                    do l=1,aotot
                        do k=1,aotot
                            do j=1,aotot
                                do i=1,aotot
                                    J_T1(i,j,k,l) = J_T1(i,j,k,l) + y(l,s)*ijks_tmp(i,j,k,s) ! J_TENSOR IN T1 BASIS
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
                deallocate(iqrs_tmp,ijrs_tmp,ijks_tmp)
                
                !!!! CALCULATING FOCK MATRIX IN T1 BASIS !!!!
                !N^2 * O
                do p=1,aotot
                    do q=1,aotot
                        FI_T1(p,q) = H_T1(p,q)
                        do i=1,O
                            FI_T1(p,q) = FI_T1(p,q) + 2*J_T1(p,q,i,i) - J_T1(p,i,i,q)
                        enddo
                    enddo
                enddo
                
            endif
             
            !!!! DONE WITH MO ==> T1 INTEGRAL TRANSFORMATION
            
            
            
            
            !!!! CALCULATING AUXILLIARY TENSORS
            allocate(Uvovo(V,O,V,O),B2oooo(O,O,O,O),C2vovo(V,O,V,O))
            allocate(D2vovo(V,O,V,O),E2oo(O,O),E2vv(V,V))
            Uvovo = 0
            B2oooo = 0
            C2vovo = 0
            D2vovo = 0
            E2oo = 0
            E2vv = 0
            
            !Uvovo O^2 * V^2
            do a=1,V
                do i=1,O
                    do b=1,V
                        do j=1,O
                            Uvovo(a,i,b,j) = 2* t2_amps(a,i,b,j) - t2_amps(a,j,b,i)
                        enddo
                    enddo
                enddo
            enddo
            !B2oooo O^4 * V^2
            do j=1,O
                do i=1,O
                    do l=1,O
                        do k=1,O
                            B2oooo(k,l,i,j) = B2oooo(k,l,i,j) + J_T1(k,i,l,j)
                            do d=1,V
                                do c=1,V
                                    B2oooo(k,l,i,j) = B2oooo(k,l,i,j) + t2_amps(c,i,d,j) * J_T1(k,O+c,l,O+d)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            !C2vovo and D2vovo O^3 * V^3
            do i=1,O
                do a=1,V
                    do k=1,O
                        do c=1,V
                            C2vovo(c,k,a,i) = C2vovo(c,k,a,i) + J_T1(k,i,O+a,O+c)
                            D2vovo(c,k,a,i) = D2vovo(c,k,a,i) + (2*J_T1(O+a,i,k,O+c) - J_T1(O+a,O+c,k,i))
                            do l=1,O
                                do d=1,V
                                    C2vovo(c,k,a,i) = C2vovo(c,k,a,i) - t2_amps(a,l,d,i) * J_T1(k,O+d,l,O+c) / 2
                                    D2vovo(c,k,a,i) = D2vovo(c,k,a,i) + Uvovo(a,i,d,l) * (2*J_T1(l,O+d,k,O+c) - J_T1(l,O+c,k,O+d)) / 2
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            !E2oo O^4 * V
            do j=1,O
                do k=1,O
                    E2oo(k,j) = E2oo(k,j) + FI_T1(k,j)
                    do l=1,O
                        do c=1,V
                            do d=1,V
                                E2oo(k,j) = E2oo(k,j) + Uvovo(c,l,d,j) * J_T1(k,O+d,l,O+c)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            !E2vv O * V^4
            do b=1,V
                do c=1,V
                    E2vv(c,b) = E2vv(c,b) + FI_T1(O+b,O+c)
                    do k=1,O
                        do l=1,O
                            do d=1,V
                                E2vv(c,b) = E2vv(c,b) - Uvovo(b,k,d,l) * J_T1(l,O+d,k,O+c)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            
            
            
            
            
            
            !!!AMPLITUDE INCREMENTS!!!!!
            do a=1,V
                ea = Eig(O+a)
                do i=1,O
                    ei = Eig(i)
        
                    !!!!ADDING T1 AMPLITUDE INCREMENTS!!!!
                    
                    P1=0
                    P2=0
                    
                    !O^2 * V^2
                    do d=1,V
                        do l=1,O
                            P1 = P1 + t2_amps(a,i,d,l) * FI_T1(l,O+d)
                        enddo
                    enddo
                    
                    P1 = 2.0d0 * P1
                                   
                    do d=1,V
                        do k=1,O
                            P2 = P2 - t2_amps(a,k,d,i) * FI_T1(k,O+d)
                            !O^3 * V^2
                            do l=1,O
                                P2 = P2 - t2_amps(a,k,d,l) * ( 2*J_T1(k,i,l,O+d) - J_T1(k,O+d,l,i) )
                            enddo
                            !O^2 * V^3
                            do c=1,V
                                P2 = P2 + t2_amps(c,k,d,i) * ( 2*J_T1(O+a,O+d,k,O+c) - J_T1(O+a,O+c,k,O+d) )
                            enddo
                            
                        enddo
                    enddo
                        
                    dt1_amps_ai = ( FI_T1(O+a,i) + (P1 + P2) )/(ei - ea)
                    !print *, dt1_amps_ai
                    t1_amps(a,i) = t1_amps(a,i) + dt1_amps_ai
                    
                    
                    
                    
                    do b=1,V
                        eb = Eig(O+b)
                        do j=1,O
                            ej = Eig(j)
                            
                            !!!!ADDING T2 AMPLITUDE INCREMENTS!!!!
                            A2=0
                            B2=0
                            C2=0
                            D2=0
                            E2=0
                            
                            !O^2 * V^4
                            do d=1,V
                                do c=1,V
                                    A2 = A2 + t2_amps(c,i,d,j)*J_T1(O+a,O+c,O+b,O+d)
                                enddo
                            enddo
                            !O^4 * V^2
                            do l=1,O
                                do k=1,O
                                    B2 = B2 + t2_amps(a,k,b,l) * B2oooo(k,l,i,j)
                                enddo
                            enddo
                            !O^3 * V^3
                            do k=1,O
                                do c=1,V
                                    C2 = C2 + t2_amps(b,k,c,j) * C2vovo(c,k,a,i) + t2_amps(a,k,c,i) * C2vovo(c,k,b,j)
                                enddo
                            enddo
                            C2 = C2 * -0.5d0
                            do k=1,O
                                do c=1,V
                                    C2 = C2 - t2_amps(b,k,c,i) * C2vovo(c,k,a,j) - t2_amps(a,k,c,j) * C2vovo(c,k,b,i)
                                    D2 = D2 +   Uvovo(b,j,c,k) * D2vovo(c,k,a,i) +   Uvovo(a,i,c,k) * D2vovo(c,k,b,j)
                                enddo
                            enddo
                            D2 = D2 * 0.5d0
                            !O^3 * V^2
                            do k=1,O
                                E2 = E2 + t2_amps(a,i,b,k) * E2oo(k,j) + t2_amps(b,j,a,k) * E2oo(k,i)
                            enddo
                            E2 = -E2
                            !O^2 * V^3
                            do c=1,V
                                E2 = E2 + t2_amps(a,i,c,j) * E2vv(c,b) + t2_amps(b,j,c,i) * E2vv(c,a)
                            enddo
                            !O^2 * V^2
                            dt2_amps_aibj = ( J_T1(O+a,i,O+b,j) + (A2 + B2 + C2 + D2 + E2) )/(ei + ej - ea - eb)
                            !print *, dt2_amps_aibj
                            t2_amps(a,i,b,j) = t2_amps(a,i,b,j) + dt2_amps_aibj
                        enddo
                    enddo
                enddo
            enddo
            
            !deallocate auxilliary tensors!!!
            deallocate(Uvovo,B2oooo,C2vovo,D2vovo,E2oo,E2vv)
        !'ccc loop cycle'    
        enddo
        
 
    !#######################################################################################################################!   
        
        
        
        
        
        
        
        
        
        
        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                                     OBTAINING CCSD ENERGY                                             !
    !=======================================================================================================================! 
        !O^2 * V^2
        Eccsd = 0
        do a=1,V
            do i=1,O
                do b=1,V
                    do j=1,O
                        Eccsd = Eccsd + ( t1_amps(a,i) * t1_amps(b,j) + t2_amps(a,i,b,j) ) * (2*J_MO(i,O+a,j,O+b) - J_MO(i,O+b,j,O+a) )
                    enddo
                enddo
            enddo
        enddo
        
        
    !#######################################################################################################################!
        
        
        if (output == 1) then
            call onelineheader("Correlation Energy Calculated!")
            write(77,*)''
            write(77,*) 'Hatree-Fock Energy: ',EHF
            write(77,*) 'CCSD Energy: ',EHF+Eccsd
            write(77,*) 'CCSD Correlation Energy: ',Eccsd
            call writelines(2)
            call twolinesfooter("CCSD CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine closed_shell_CCSD_calc
    !#######################################################################################################################!
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    subroutine closed_shell_CCD_calc_slow(output)
    integer , intent(in) :: output
    
    integer :: i,j,k,l,a,b,c,d,p,q,r,s,calcf,cc,O,V
    integer , dimension(aotot) :: mo_indices
    
    real (kind = 8) :: ei,ej,ea,eb,Q1,Q2,R1,R2,R3,dt_amps_aibj
    real (kind = 8) , dimension(aotot,aotot) :: H_MO,FI
    real (kind = 8) , dimension(aotot,aotot,aotot,aotot) :: J_MO!,iajb_array2
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: t_amps
    logical :: nott=0
    
    O = int(eltot/2)
    V = aotot - O

        if (output == 1) then
            call smallblockheader('CLOSED SHELL CCD CORRELATION CALCULATION')
            call writelines(2)
            call onelineheader('Starting CCD routine')
            call writelines(2)
        endif
        
        Eccd = 0
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        call two_index_transform(aotot,Cmo,H,H_MO)         !H_MO = matmul( transpose(Cmo), matmul(H,Cmo) ) !H MATRIX IN MO BASIS
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS

    !#######################################################################################################################!
        

        
        
        
        
        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                             OBTAINING CCD CLUSTER AMPLITUDES                                          !
    !=======================================================================================================================!  
        allocate(t_amps(V,O,V,O))
        t_amps = 0
        !!!! GENERATE FIRST ORDER (MP2) AMPLITUDES (t0 = 0) !!!!
        do a=1,V
            ea = Eig(O+a)
            do i=1,O
                ei = Eig(i)
                do b=1,V
                    eb = Eig(O+b)
                    do j=1,O
                        ej = Eig(j)
                        t_amps(a,i,b,j) = J_MO(O+a,i,O+b,j)/(ei+ej-ea-eb)
                    enddo
                enddo
            enddo
        enddo
        

        !!!! GENERATE CCD AMPLITUDES VIA PERTUBATION-BASED QUASI-NEWTON METHOD
        do cc=1,20
            
            do a=1,V
                ea = Eig(O+a)
                do i=1,O
                    ei = Eig(i)
                    do b=1,V
                        eb = Eig(O+b)
                        do j=1,O
                            ej = Eig(j)
                            Q1=0
                            Q2=0
                            R1=0
                            R2=0
                            R3=0
                            
                            Q1 =  t_amps(a,i,b,j) * (eb - ej) + t_amps(b,j,a,i) * (ea - ei)
                            do c=1,V
                                do k=1,O
                                    Q1 = Q1 + t_amps(a,i,c,k)*(2*J_MO(O+b,j,k,O+c) - J_MO(O+b,O+c,k,j)) + t_amps(b,j,c,k)*(2*J_MO(O+a,i,k,O+c) - J_MO(O+a,O+c,k,i))
                                enddo
                            enddo
                            
                            
                            
                            
         
                            
                            do c=1,V
                                do k=1,O
                                    Q2 = Q2 + t_amps(a,k,c,j)*J_MO(O+b,O+c,k,i) + t_amps(a,k,c,i)*J_MO(O+b,j,k,O+c) + t_amps(b,k,c,i)*J_MO(O+a,O+c,k,j) + t_amps(b,k,c,j)*J_MO(O+a,i,k,O+c)
                                enddo
                            enddo
                            Q2 = -Q2
                            
                            do c=1,V
                                do d=1,V
                                    Q2 = Q2 + t_amps(c,i,d,j)*J_MO(O+a,O+c,O+b,O+d)
                                enddo
                            enddo
                            
                            do k=1,O
                                do l=1,O
                                    Q2 = Q2 + t_amps(a,k,b,l)*J_MO(k,i,l,j)
                                enddo
                            enddo
                            
                            
                  
                            
                            do c=1,V
                                do d=1,V
                                    do k=1,O
                                        do l=1,O
                                            R1 = R1 + ( t_amps(a,i,c,k)*t_amps(b,j,d,l) + t_amps(b,j,c,k)*t_amps(a,i,d,l) ) * ( 2*J_MO(k,O+c,l,O+d) - J_MO(k,O+d,l,O+c) )
                                            
                                            R2 = R2 + ( t_amps(a,i,b,k)*t_amps(c,l,d,j) + t_amps(b,j,a,k)*t_amps(c,l,d,i) ) * ( 2*J_MO(l,O+c,k,O+d) - J_MO(l,O+d,k,O+c) )
                                            R2 = R2 + ( t_amps(a,i,c,k)*t_amps(b,l,d,j) + t_amps(b,j,c,k)*t_amps(a,l,d,i) ) * ( 2*J_MO(k,O+c,l,O+d) - J_MO(k,O+d,l,O+c) )
                                            R2 = R2 + ( t_amps(a,i,c,j)*t_amps(b,k,d,l) + t_amps(b,j,c,i)*t_amps(a,k,d,l) ) * ( 2*J_MO(k,O+c,l,O+d) - J_MO(k,O+d,l,O+c) )
                                            
                                            R3 = R3 + ( t_amps(a,k,b,l)*t_amps(c,i,d,j) + t_amps(b,k,a,l)*t_amps(c,j,d,i) ) * J_MO(k,O+c,l,O+d)
                                            R3 = R3 + ( t_amps(a,k,c,i)*t_amps(b,l,d,j) + t_amps(b,k,c,j)*t_amps(a,l,d,i) ) * J_MO(k,O+c,l,O+d)
                                            R3 = R3 + ( t_amps(a,k,d,j)*t_amps(b,l,c,i) + t_amps(b,k,d,i)*t_amps(a,l,c,j) ) * J_MO(k,O+c,l,O+d)
                                        
                                        enddo
                                    enddo
                                enddo
                            enddo
                            
                            R2 = -R2
                            R3 = R3/2
                            
                            
         
                            
                            dt_amps_aibj = ( J_MO(O+a,i,O+b,j) + (Q1 + Q2 + R1 + R2 + R3) )/(ei + ej - ea - eb)
                            t_amps(a,i,b,j) = t_amps(a,i,b,j) + dt_amps_aibj

                        enddo
                    enddo
                enddo
            enddo
            
            
        enddo
        
 
    !#######################################################################################################################!   
        
        
        
        
        
        
        
        
        
  
        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                                     OBTAINING CCD ENERGY                                              !
    !=======================================================================================================================! 
        Eccd = 0
        do a=1,V
            do i=1,O
                do b=1,V
                    do j=1,O
                        Eccd = Eccd + t_amps(a,i,b,j) * (2*J_MO(i,O+a,j,O+b) - J_MO(i,O+b,j,O+a) )
                    enddo
                enddo
            enddo
        enddo
        

        
    !#######################################################################################################################!
        
        
        if (output == 1) then
            call onelineheader("Correlation Energy Calculated!")
            write(77,*)''
            write(77,*) 'Hatree-Fock Energy: ',EHF
            write(77,*) 'CCD Energy: ',EHF+Eccd
            write(77,*) 'CCD Correlation Energy: ',Eccd
            call writelines(2)
            call twolinesfooter("CCD CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine closed_shell_CCD_calc_slow
    !#######################################################################################################################!
    
    
    
    
    
    
    
    
    
    !#######################################################################################################################!  
    subroutine closed_shell_CCSD_calc_slow(output)
    integer , intent(in) :: output
    
    integer :: i,j,k,l,a,b,c,d,p,q,r,s,calcF,ccc,O,V,trs=0
    integer , dimension(aotot) :: mo_indices
    
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: iqrs_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ijrs_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ijks_tmp
    
    real (kind = 8) :: ei,ej,ea,eb,P1,P2,Q1,Q2,R1,R2,R3,dt1_amps_ai,dt2_amps_aibj
    real (kind = 8) , dimension(aotot,aotot) :: H_MO,FI_MO,H_T1,FI_T1,t1,x,y
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO,J_T1!,iajb_array2
    
    
    real (kind = 8) , allocatable, dimension(:,:) :: t1_amps
    real (kind = 8) , allocatable, dimension(:,:,:,:) :: t2_amps
        O = int(eltot/2) 
        V = aotot - O 
        if (output == 1) then
            call smallblockheader('CLOSED SHELL CCSD CORRELATION CALCULATION')
            call writelines(2)
            call onelineheader('Starting CCSD routine')
            call writelines(2)
        endif
        
  
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        H_MO = matmul( transpose(Cmo), matmul(H,Cmo) ) !H MATRIX IN MO BASIS
        allocate(J_MO(aotot,aotot,aotot,aotot),J_T1(aotot,aotot,aotot,aotot))
        
        if (trs/=0) then
        allocate(iqrs_tmp(aotot,aotot,aotot,aotot),ijrs_tmp(aotot,aotot,aotot,aotot),ijks_tmp(aotot,aotot,aotot,aotot))
        
        iqrs_tmp = 0 !MAKE J TENSOR IN MO BASIS, ijkl are general indices here (unlike convinient occupied incices!!!)
        ijrs_tmp = 0
        ijks_tmp = 0
        J_MO = 0
        do s=1,aotot
            do r=1,aotot
                do q=1,aotot
                    do p=1,aotot
                        do i=1,aotot
                            iqrs_tmp(i,q,r,s) = iqrs_tmp(i,q,r,s) + Cmo(p,i)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                        enddo
                    enddo
                enddo
            enddo
        enddo


        do s=1,aotot
            do r=1,aotot
                do q=1,aotot
                    do j=1,aotot
                        do i=1,aotot
                            ijrs_tmp(i,j,r,s) = ijrs_tmp(i,j,r,s) + Cmo(q,j)*iqrs_tmp(i,q,r,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        

        do s=1,aotot
            do r=1,aotot
                do k=1,aotot
                    do j=1,aotot
                        do i=1,aotot
                            ijks_tmp(i,j,k,s) = ijks_tmp(i,j,k,s) + Cmo(r,k)*ijrs_tmp(i,j,r,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    
        do s=1,aotot
            do l=1,aotot
                do k=1,aotot
                    do j=1,aotot
                        do i=1,aotot
                            J_MO(i,j,k,l) = J_MO(i,j,k,l) + Cmo(s,l)*ijks_tmp(i,j,k,s) ! J_TENSOR IN MO BASIS
                        enddo
                    enddo
                enddo
            enddo
        enddo
        endif
        H_MO = 0
        J_MO = 0
        call two_index_transform(aotot,Cmo,H,H_MO)         !H_MO = matmul( transpose(Cmo), matmul(H,Cmo) ) !H MATRIX IN MO BASIS
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
        
        FI_MO = 0
        do p=1,aotot
            do q=1,aotot
                FI_MO(p,q) = H_MO(p,q)
                do i=1,O
                    FI_MO(p,q) = FI_MO(p,q) + 2*J_MO(p,q,i,i) - J_MO(p,i,i,q)
                enddo
            enddo
        enddo
        

        
    !#######################################################################################################################!   
        
        
        
        
        
        
        
        
        
        
        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                             OBTAINING CCSD CLUSTER AMPLITUDES                                         !
    !=======================================================================================================================!  

        allocate(t1_amps(V,O))
        allocate(t2_amps(V,O,V,O))
        allocate(iqrs_tmp(aotot,aotot,aotot,aotot),ijrs_tmp(aotot,aotot,aotot,aotot),ijks_tmp(aotot,aotot,aotot,aotot))
        t1_amps = 0
        t2_amps = 0
        !!!! GENERATE FIRST ORDER (MP2) AMPLITUDES (t0 = 0) !!!!
        !!!! NOTE THERE IS NO t=0 SINGLES FIRST ORDER CONTRIBUTION DUE TO A DIAGONAL FOCK MATRIX IN THE MO BASIS!!!!
        do a=1,V
            ea = Eig(O+a)
            do i=1,O
                ei = Eig(i)
                do b=1,V
                    eb = Eig(O+b)
                    do j=1,O
                        ej = Eig(j)
                        t2_amps(a,i,b,j) = J_MO(O+a,i,O+b,j)/(ei+ej-ea-eb)
                    enddo
                enddo
            enddo
        enddo
        

        !!!! GENERATE CCD AMPLITUDES VIA PERTUBATION-BASED QUASI-NEWTON METHOD
        
        do ccc=1,10
            !print *, 'kk'
            
            !!!! INVOKING T1 TRANSFORMATION !!!!        
            if (ccc == 1) then
                H_T1 = H_MO
                J_T1 = J_MO
                FI_T1 = FI_MO
            else

                !!!! CREATE T1, X AND Y MATRICES !!!!
                t1 = 0
                do a=1,V
                    do i=1,O
                        t1(a+O,i) = t1_amps(a,i)
                    enddo
                enddo
                

                
                x=0
                y=0
                
                do i=1,aotot
                    x(i,i) = 1
                    y(i,i) = 1
                enddo
                
                x = x - t1
                y = y + transpose(t1)
                
                !print *, t1

            !=======================================================================================================================!
            !                                          INTEGRAL TRANSFORMATION (MO ==> T1)                                          !
            !=======================================================================================================================!
        
                H_T1 = matmul( matmul(x,H_MO),transpose(y) ) !H MATRIX IN T1 BASIS
        
                iqrs_tmp = 0 ! MAKE J TENSOR IN T1 BASIS, ijkl are general indices here (unlike convinient occupied incices!!!)
                ijrs_tmp = 0
                ijks_tmp = 0
                J_T1 = 0
                do s=1,aotot
                    do r=1,aotot
                        do q=1,aotot
                            do p=1,aotot
                                do i=1,aotot
                                    iqrs_tmp(i,q,r,s) = iqrs_tmp(i,q,r,s) + x(i,p)*J_MO(p,q,r,s)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo


                do s=1,aotot
                    do r=1,aotot
                        do q=1,aotot
                            do j=1,aotot
                                do i=1,aotot
                                    ijrs_tmp(i,j,r,s) = ijrs_tmp(i,j,r,s) + y(j,q)*iqrs_tmp(i,q,r,s)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
        

                do s=1,aotot
                    do r=1,aotot
                        do k=1,aotot
                            do j=1,aotot
                                do i=1,aotot
                                    ijks_tmp(i,j,k,s) = ijks_tmp(i,j,k,s) + x(k,r)*ijrs_tmp(i,j,r,s)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
    
                do s=1,aotot
                    do l=1,aotot
                        do k=1,aotot
                            do j=1,aotot
                                do i=1,aotot
                                    J_T1(i,j,k,l) = J_T1(i,j,k,l) + y(l,s)*ijks_tmp(i,j,k,s) ! J_TENSOR IN T1 BASIS
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
                
                
                !!!! CALCULATING FOCK MATRIX IN T1 BASIS !!!!
        
                do p=1,aotot
                    do q=1,aotot
                        FI_T1(p,q) = H_T1(p,q)
                        do i=1,O
                            FI_T1(p,q) = FI_T1(p,q) + 2*J_T1(p,q,i,i) - J_T1(p,i,i,q)
                        enddo
                    enddo
                enddo
                
            endif
             
            !!!! DONE WITH MO ==> T1 INTEGRAL TRANSFORMATION
            
            
            do a=1,V
                ea = Eig(O+a)
                do i=1,O
                    ei = Eig(i)

                    
                    !!!!ADDING T1 AMPLITUDE INCREMENTS!!!!
                    
                    P1=0
                    P2=0
                    
                    do d=1,V
                        do l=1,O
                            P1 = P1 + t2_amps(a,i,d,l) * FI_T1(l,O+d)
                        enddo
                    enddo
                    
                    P1 = 2 * P1
                    
                    
                    do d=1,V
                        do k=1,O
                            
                            P2 = P2 - t2_amps(a,k,d,i) * FI_T1(k,O+d)
                            
                            do l=1,O
                                P2 = P2 - t2_amps(a,k,d,l) * ( 2*J_T1(k,i,l,O+d) - J_T1(k,O+d,l,i) )
                            enddo
                            
                            do c=1,V
                                P2 = P2 + t2_amps(c,k,d,i) * ( 2*J_T1(O+a,O+d,k,O+c) - J_T1(O+a,O+c,k,O+d) )
                            enddo
                            
                        enddo
                    enddo
                    
                    
                    
                    dt1_amps_ai = ( FI_T1(O+a,i) + (P1 + P2) )/(ei - ea)
                    !print *, dt1_amps_ai
                    t1_amps(a,i) = t1_amps(a,i) + dt1_amps_ai
                    
                    do b=1,V
                        eb = Eig(O+b)
                        do j=1,O
                            ej = Eig(j)
                            
                            !!!!ADDING T2 AMPLITUDE INCREMENTS!!!!
                            
                            Q1=0
                            Q2=0
                            R1=0
                            R2=0
                            R3=0
                            
                            do c=1,V
                                Q1 = Q1 + ( t2_amps(a,i,c,j) * FI_T1(O+b,O+c) + t2_amps(b,j,c,i) * FI_T1(O+a,O+c) )
                            enddo
                            
                            do k=1,O
                                Q1 = Q1 - t2_amps(a,i,b,k) * FI_T1(k,j)
                                Q1 = Q1 - t2_amps(b,j,a,k) * FI_T1(k,i) 
                            enddo
                            
                            do c=1,V
                                do k=1,O
                                    Q1 = Q1 + t2_amps(a,i,c,k)*(2*J_T1(O+b,j,k,O+c) - J_T1(O+b,O+c,k,j)) 
                                    Q1 = Q1 + t2_amps(b,j,c,k)*(2*J_T1(O+a,i,k,O+c) - J_T1(O+a,O+c,k,i))
                                enddo
                            enddo
                            
    
                            
                            do c=1,V
                                do k=1,O
                                    Q2 = Q2 + t2_amps(a,k,c,j)*J_T1(O+b,O+c,k,i) + t2_amps(a,k,c,i)*J_T1(O+b,j,k,O+c) 
                                    Q2 = Q2 + t2_amps(b,k,c,i)*J_T1(O+a,O+c,k,j) + t2_amps(b,k,c,j)*J_T1(O+a,i,k,O+c)
                                enddo
                            enddo
                            Q2 = -Q2
                            
                            do c=1,V
                                do d=1,V
                                    Q2 = Q2 + t2_amps(c,i,d,j)*J_T1(O+a,O+c,O+b,O+d)
                                enddo
                            enddo
                            
                            do k=1,O
                                do l=1,O
                                    Q2 = Q2 + t2_amps(a,k,b,l)*J_T1(k,i,l,j)
                                enddo
                            enddo
                            
                            
                            
                            
                            
                            
                            
                            do c=1,V
                                do d=1,V
                                    do k=1,O
                                        do l=1,O
                                            R1 = R1 + ( t2_amps(a,i,c,k)*t2_amps(b,j,d,l) + t2_amps(b,j,c,k)*t2_amps(a,i,d,l) ) * ( 2*J_T1(k,O+c,l,O+d) - J_T1(k,O+d,l,O+c) )
                                            
                                            
                                            R2 = R2 + ( t2_amps(a,i,b,k)*t2_amps(c,l,d,j) + t2_amps(b,j,a,k)*t2_amps(c,l,d,i) ) * ( 2*J_T1(l,O+c,k,O+d) - J_T1(l,O+d,k,O+c) )
                                            R2 = R2 + ( t2_amps(a,i,c,k)*t2_amps(b,l,d,j) + t2_amps(b,j,c,k)*t2_amps(a,l,d,i) ) * ( 2*J_T1(k,O+c,l,O+d) - J_T1(k,O+d,l,O+c) )
                                            R2 = R2 + ( t2_amps(a,i,c,j)*t2_amps(b,k,d,l) + t2_amps(b,j,c,i)*t2_amps(a,k,d,l) ) * ( 2*J_T1(k,O+c,l,O+d) - J_T1(k,O+d,l,O+c) )
                                            
                                            R3 = R3 + ( t2_amps(a,k,b,l)*t2_amps(c,i,d,j) + t2_amps(b,k,a,l)*t2_amps(c,j,d,i) ) * J_T1(k,O+c,l,O+d)
                                            R3 = R3 + ( t2_amps(a,k,c,i)*t2_amps(b,l,d,j) + t2_amps(b,k,c,j)*t2_amps(a,l,d,i) ) * J_T1(k,O+c,l,O+d)
                                            R3 = R3 + ( t2_amps(a,k,d,j)*t2_amps(b,l,c,i) + t2_amps(b,k,d,i)*t2_amps(a,l,c,j) ) * J_T1(k,O+c,l,O+d)
                                        
                                        enddo
                                    enddo
                                enddo
                            enddo
                            
                            R2 = -R2
                            R3 = R3/2
                            
                            
                            
                            
                            
                            
                            
                            
                            dt2_amps_aibj = ( J_T1(O+a,i,O+b,j) + (Q1 + Q2 + R1 + R2 + R3) )/(ei + ej - ea - eb)
                            !print *, dt2_amps_aibj
                            t2_amps(a,i,b,j) = t2_amps(a,i,b,j) + dt2_amps_aibj
                        enddo
                    enddo
                enddo
            enddo

            
        enddo
        
 
    !#######################################################################################################################!   
        
        
        
        
        
        
        
        
        
        
        
        
    !#######################################################################################################################!
          
    !=======================================================================================================================!
    !                                                     OBTAINING CCSD ENERGY                                              !
    !=======================================================================================================================! 
        Eccsd = 0
        do a=1,V
            do i=1,O
                do b=1,V
                    do j=1,O
                        Eccsd = Eccsd + ( t1_amps(a,i) * t1_amps(b,j) + t2_amps(a,i,b,j) ) * (2*J_MO(i,O+a,j,O+b) - J_MO(i,O+b,j,O+a) )
                    enddo
                enddo
            enddo
        enddo
        
        
    !#######################################################################################################################!
        
        
        if (output == 1) then
            call onelineheader("Correlation Energy Calculated!")
            write(77,*)''
            write(77,*) 'Hatree-Fock Energy: ',EHF
            write(77,*) 'CCSD Energy: ',EHF+Eccsd
            write(77,*) 'CCSD Correlation Energy: ',Eccsd
            call writelines(2)
            call twolinesfooter("CCSD CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine closed_shell_CCSD_calc_slow
end module CC_correlation
!#######################################################################################################################!