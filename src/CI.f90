!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                     CI CALCULATION                                                    !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!   
module CI_correlation
use Print_module
use molecular_data
use result_data
use integral_tensors
use matvec_operations
use SCF_matrix_builder 
use Integral_Transformation
implicit none
save
contains
    subroutine closed_shell_CID_calc(output)
    integer , intent(in) :: output
    
    integer :: i,j,k,l,p,q,r,s,t,detno,d,a,O,V
    integer , dimension(aotot) :: mo_indices
    
    real (kind = 8) :: Ecorr
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: iqrs_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ijrs_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ijks_tmp
    
    real (kind = 8) , allocatable , dimension(:,:) :: H_ij
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_ijkl!,iajb_array2
    real (kind = 8) , allocatable , dimension(:,:) :: HCI,C_CI
    real (kind = 8) , allocatable , dimension(:) :: E_CI
    
        O = int(eltot/2)
        V = aotot - O
        
        allocate(E_CI(1+O*V))
        allocate(HCI(1+O*V,1+O*V),C_CI(1+O*V,1+O*V))
        allocate(H_ij(aotot,aotot))
        allocate(iqrs_tmp(aotot,aotot,aotot,aotot),ijrs_tmp(aotot,aotot,aotot,aotot))
        allocate(ijks_tmp(aotot,aotot,aotot,aotot),J_ijkl(aotot,aotot,aotot,aotot))
        if (output == 1) then
            call smallblockheader('CLOSED SHELL CID CORRELATION CALCULATION')
            call writelines(2)
            call onelineheader('Starting CI routine')
            call writelines(2)
        endif
        
        detno = 1 + O*V
        
        Ecid = 0
        
        H_ij = matmul( transpose(Cmo), matmul(H,Cmo) ) !H MATRIX IN MO BASIS
        
        iqrs_tmp = 0
        ijrs_tmp = 0
        ijks_tmp = 0
        J_ijkl = 0
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
                            J_ijkl(i,j,k,l) = J_ijkl(i,j,k,l) + Cmo(s,l)*ijks_tmp(i,j,k,s) ! J_TENSOR IN MO BASIS
                        enddo
                    enddo
                enddo
            enddo
        enddo
        deallocate(iqrs_tmp,ijrs_tmp,ijks_tmp)
        HCI = 0
        HCI(1,1) = Eel ! HARTREE-FOCK ENERGY <0|H|0> AT 1,1 POSITION OF CI MATRIX
        
        do i=1,aotot
            mo_indices(i) = i
        enddo
        
        do i=1,O
            do a=1,V
                HCI( 1 , V*(i-1)+a+1 ) = J_ijkl(i,a+O,i,a+O) ! OFF-DIAGIONAL ELEMENTS <0|H|D> IN FIRST ROW OF CI MATRIX
                HCI( V*(i-1)+a+1  , 1) = J_ijkl(i,a+O,i,a+O) ! OFF-DIAGIONAL ELEMENTS <D|H|0> IN FIRST COLUMN OF CI MATRIX
            enddo
        enddo
        
        do i=1,O
            do a=1,V
                mo_indices(i) = a+O         ! MAKING A DOUBLE EXCITATIONS i ==> a (both alpha and beta!)
                d = V*(i-1)+a+1
                do s=1,O
                                            ! DIAGONAL ELEMENTS <D|H|D> OF CI MATRIX
                    k = mo_indices(s)
                    HCI(d,d) = HCI(d,d) + 2*H_ij(k,k)
                    do t=1,O
                        l = mo_indices(t)
                        HCI(d,d) = HCI(d,d) + 2*J_ijkl(k,k,l,l) - J_ijkl(k,l,l,k)
                    enddo
                    
                enddo
                
                
                mo_indices(i) = i           ! RELAXATION a ==> i (both alpha and beta!)
            enddo
        enddo
        deallocate(J_ijkl,H_ij)
        do i=1,(O*V + 1)
            HCI(i,i) = HCI(i,i) + Vnn
        enddo
            
        call diag(detno,HCI,E_CI,C_CI)
        
        Ecid = E_CI(1) - (EHF)  
        
        if (output == 1) then
            call writelines(2)
            call onelineheader("CID matrix")
            call matprint(detno,HCI)
            call writelines(2)
        
            call writelines(2)
            call onelineheader("CID energies")
            call vecprint(detno,E_CI)
            call writelines(2)
            call onelineheader("CID coefficient matrix")
            call matprint(detno,C_CI)
            call writelines(2)
          
            call onelineheader("Correlation Energy Calculated!")
            write(77,*)''
            write(77,*) 'Hatree-Fock Energy: ',EHF
            write(77,*) 'CI Energy: ',EHF+Ecid
            write(77,*) 'CI Correlation Energy: ',Ecid
            call writelines(2)
            call twolinesfooter("CI CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine closed_shell_CID_calc
    
!#######################################################################################################################!     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################! 
    
    subroutine closed_shell_CISES_calc(output)
    use nuclei_data
    integer , intent(in) :: output
    
    integer :: i,j,k,l,p,q,r,s,t,detno,d,a,b,O,V
    integer , dimension(aotot) :: mo_indices
    
    
    real (kind = 8) , dimension(3) :: DP,edtm
    real (kind = 8) , allocatable , dimension(:,:) :: DPx,DPy,DPz,DPx_MO,DPy_MO,DPz_MO
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO!,iajb_array2
    real (kind = 8) , allocatable , dimension(:,:) :: HsCI,Cs_CI,HtCI,Ct_CI
    real (kind = 8) , allocatable , dimension(:) :: Ws_CI,Wt_CI
    real (kind = 8) :: pi=3.14159265359, fn, fntot
    
        O = int(eltot/2)
        V = aotot - O
        
        allocate(Ws_CI(O*V),Wt_CI(O*V))
        allocate(HsCI(O*V,O*V),HtCI(O*V,O*V),Cs_CI(O*V,O*V),Ct_CI(O*V,O*V))
        allocate(DPx(aotot,aotot),DPy(aotot,aotot),DPz(aotot,aotot))
        allocate(DPx_MO(aotot,aotot),DPy_MO(aotot,aotot),DPz_MO(aotot,aotot))
        allocate(J_MO(aotot,aotot,aotot,aotot))
        if (output == 1) then
            call smallblockheader('CLOSED SHELL CIS CALCULATION')
            call writelines(2)
            call onelineheader('Starting CIS routine')
            call writelines(2)
        endif
        
        detno = O*V

    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
    !#######################################################################################################################!   
        
        
    !=======================================================================================================================!
    !                                          BUILDING HCI MATRIX AND PERFORM CIS                                          !
    !=======================================================================================================================!  
          
        HsCI = 0
        HtCI = 0
        do i=1,O
            do a=1,V
                HsCI( V*(i-1)+a , V*(i-1)+a ) = HsCI( V*(i-1)+a , V*(i-1)+a ) + ( Eig(O+a) - Eig(i) )
                HtCI( V*(i-1)+a , V*(i-1)+a ) = HtCI( V*(i-1)+a , V*(i-1)+a ) + ( Eig(O+a) - Eig(i) )
                do j=1,O
                    do b=1,V
                        HsCI(V*(i-1)+a,V*(j-1)+b) = HsCI(V*(i-1)+a,V*(j-1)+b) + 2*J_MO(a+O,i,j,b+O) - J_MO(a+O,b+O,j,i)
                        HtCI(V*(i-1)+a,V*(j-1)+b) = HtCI(V*(i-1)+a,V*(j-1)+b) + - J_MO(a+O,b+O,j,i)
                    enddo
                enddo
                
            enddo
        enddo
        
        deallocate(J_MO)
        call diag(detno,HsCI,Ws_CI,Cs_CI)
        call diag(detno,HtCI,Wt_CI,Ct_CI)
        

    !#######################################################################################################################!
        
        
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                           OBTAINING TRANSITION PROBABILITIES                                          !
    !=======================================================================================================================!
      
    call DP_calc(DPx,DPy,DPz) 
    DPx_MO = matmul( transpose(Cmo), matmul(DPx,Cmo) ) !DP MATRIX IN MO BASIS
    DPy_MO = matmul( transpose(Cmo), matmul(DPy,Cmo) ) !DP MATRIX IN MO BASIS
    DPz_MO = matmul( transpose(Cmo), matmul(DPz,Cmo) ) !DP MATRIX IN MO BASIS
    deallocate(DPx,DPy,DPz)
    DP = 0
    do i=1,O
        DP(1) =  DP(1) + 2 * DPx_MO(i,i)
        DP(2) =  DP(2) + 2 * DPy_MO(i,i)
        DP(3) =  DP(3) + 2 * DPz_MO(i,i)
    enddo 
    DP = -DP   
    do j=1,nuctot
        DP =  DP + nuclei(j)%charge * nuclei(j)%position
    enddo        
    write(77,*) 'dipole moment:', ( (DP(1)**2 + DP(2)**2 + DP(3)**2)**0.5 ) * 2.541746473, 'D'
        
    call onelineheader("Obtaining Singlet Transition Probabilities")
    fntot = 0
    do j = 1,detno
        edtm = 0
        do i=1,O
            do a=1,V
                edtm(1) = edtm(1) + Cs_CI( V*(i-1)+a , j ) * DPx_MO(i,a+O)
                edtm(2) = edtm(2) + Cs_CI( V*(i-1)+a , j ) * DPy_MO(i,a+O)
                edtm(3) = edtm(3) + Cs_CI( V*(i-1)+a , j ) * DPz_MO(i,a+O)
            enddo
        enddo
        
        edtm = edtm * -2**0.5
        fn = (2*Ws_CI(j) / 3) * ( (edtm(1)**2 + edtm(2)**2 + edtm(3)**2) )
        write(77,*) 'f:',  fn , 'at' ,( (137*2*pi)/Ws_CI(j) ) * (5.29177210903 / 100) , 'nm'! , Ws_CI(j) * 27, 'eV)' !* 2.541746473, 'D'
        fntot = fntot + fn
    enddo
    deallocate(DPx_MO,DPy_MO,DPz_MO)    
    write(77,*) 'Total sum of oscillator strengths:', fntot    
    !#######################################################################################################################!
        

        if (output == 1) then
            !call writelines(2)
            !call onelineheader("CISES matrix")
            !call matprint(detno,HCI)
            !call writelines(2)
        
            call writelines(2)
            call onelineheader("Singlet CIS Excitation Energies")
            call vecprint(detno,Ws_CI*27.2114079527)
            call writelines(2)
            call onelineheader("Triplet CIS Excitation Energies")
            call vecprint(detno,Wt_CI*27.2114079527)
            call writelines(2)
            !call onelineheader("CISES coefficient matrix")
            !call matprint(detno,C_CI)
            !call writelines(2)
          
            call onelineheader("CISES SPECTRA Calculated!")
            write(77,*)''
            call writelines(2)
            call twolinesfooter("CISES CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine closed_shell_CISES_calc
    
    
    
!#######################################################################################################################! 
    
    subroutine CIS_lowest_davidson_calc(output)
    use nuclei_data
    integer , intent(in) :: output
    
    integer :: i,a,j,b,O,V,OV,rt,idx,l,k
    integer :: nrts ,ss
    
    real(kind = 8) :: tol1=1E-6, tol2 = 1E-9, theta, tvecpn, resn, d
    
    integer , allocatable , dimension(:) :: rtidx
    
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO!,iajb_array2
    
    real (kind = 8) , allocatable , dimension(:,:) :: X,AX,M,P,S
    
    real (kind = 8) , allocatable , dimension(:) :: SPA,W,Ws,tvec,tvecp,res,Au,u
    
    real (kind = 8) :: pi=3.14159265359
    
        O = int(eltot/2)
        V = aotot - O
        OV = O*V
        nrts = OV/2
        ss = 10
        
        allocate(SPA(OV),tvec(OV),tvecp(OV),Au(OV),u(OV))
        allocate(X(OV,ss),AX(OV,ss))
        allocate(W(nrts))
        allocate(J_MO(aotot,aotot,aotot,aotot))
        if (output == 1) then
            call smallblockheader('CLOSED SHELL CIS CALCULATION')
            call writelines(2)
            call onelineheader('Starting CIS routine')
            call writelines(2)
        endif
        


    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
    !#######################################################################################################################!   
        
        do i=1,O
            do a=1,V
                SPA( V*(i-1)+a) = ( Eig(O+a) - Eig(i) ) + 2*J_MO(a+O,i,i,a+O) - J_MO(a+O,a+O,i,i)
            enddo
        enddo

        do i=1,nrts 
            rtidx(i) = minloc(SPA,dim=1,mask=(SPA > 0))
            SPA(rtidx(i)) = SPA(rtidx(i)) - 100000
            print *, rtidx(i)
        enddo

        do i=1,nrts
            SPA(rtidx(i)) = SPA(rtidx(i)) + 100000
        enddo

    !=======================================================================================================================!
    !                                           PERFORM CIS VIA THE DAVIDSON ALGORITHM                                      !
    !=======================================================================================================================!  
        allocate(P(OV,OV))
        do rt=1,1!nrts 
            !LOOPING OVER N ROOTS
            idx = rtidx(rt)
            
            !SETTING DAVIDSON MATRICES ZERO BEFORE CALCULATING  ROOT
            X = 0
            AX = 0  
            X(idx,1) = 1
            
            !RETRIEVING DIRECTION VECTORS FOR FIRST GUESS
            do i=1,O
                do a=1,V
                    
                    AX(V*(i-1)+a,1) = AX(V*(i-1)+a,1) + ( Eig(O+a) - Eig(i) ) *  X(V*(i-1)+a,1)

                    do j=1,O
                        do b=1,V
                            
                            AX(V*(i-1)+a,1) = AX(V*(i-1)+a,1) + (2*J_MO(a+O,i,j,b+O) - J_MO(a+O,b+O,j,i)) * X(V*(j-1)+b,1) 
                            
                        enddo
                    enddo
                enddo
            enddo
            
            theta = 0
            do i=1,OV
                theta = theta + X(i,1) * AX(i,1)
            enddo
            
            res = AX(:,1) - theta * X(:,1)
                
            do l=1,OV
                d = theta - SPA(l)
                if (abs(d) < tol1) then
                    tvec(l) = 0
                else
                    tvec(l) = res(l)/d
                endif
            enddo
            
            !CARRYING OUT WITH THE DAVIDSON ALGORITHM FOR A SPECIFIC ROOT
            do k=2,ss-1
                
                print *, 'k:',k

                P = 0
                P = -matmul(X(:,1:k-1),transpose(X(:,1:k-1)))
                do l=1,OV
                    P(l,l) = P(l,l) + 1
                enddo
                
                tvecp = matmul(P,tvec)
                
                tvecpn = 0
                do l=1,OV
                    tvecpn = tvecpn + tvecp(l)**2
                enddo
                tvecpn = tvecpn**0.5
                
                X(:,k) = tvecp/tvecpn
                
                do i=1,O
                    do a=1,V
                    
                        AX(V*(i-1)+a,k) = AX(V*(i-1)+a,k) + ( Eig(O+a) - Eig(i) ) *  X(V*(i-1)+a,k)

                        do j=1,O
                            do b=1,V
                            
                                AX(V*(i-1)+a,k) = AX(V*(i-1)+a,k) + (2*J_MO(a+O,i,j,b+O) - J_MO(a+O,b+O,j,i)) * X(V*(j-1)+b,k) 
                            
                            enddo
                        enddo
                    enddo
                enddo
                
                allocate(M(k,k))
                allocate(S(k,k))
                allocate(Ws(k))
                M=0
                M = matmul(transpose(X(:,1:k)),AX(:,1:k))
            
                call diag(k,M,Ws,S)
                
                u = matmul(X(:,1:k),S(:,1))
                Au = matmul(AX(:,1:k),S(:,1))
                
                theta = Ws(1)
                
                res = Au - theta*u
                
                resn = 0
                do l=1,OV
                    resn = resn + res(l)**2
                enddo
                resn = resn**0.5
                
                if (resn < tol2) then
                    exit
                endif      
                
                do l=1,OV
                    d = theta - SPA(l)
                    if (abs(d) < tol1) then
                        tvec(l) = 0
                    else
                        tvec(l) = res(l)/d
                    endif
                enddo
        
                
                deallocate(M,S,Ws)
            enddo

            
            
        enddo  
        deallocate(J_MO)

    !#######################################################################################################################!
        

        if (output == 1) then
        
            call writelines(2)
            call onelineheader("Singlet CIS Excitation Energies")
            
            print *, 'theta:', theta*27
          
            call onelineheader("CISES SPECTRA Calculated!")
            write(77,*)''
            call writelines(2)
            call twolinesfooter("CISES CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine CIS_lowest_davidson_calc
    
    
    
    
    
    
!#######################################################################################################################! 
    
    
    
    
    
    
    
    
    
    
    subroutine CIS_multiple_davidson_calc(output)
    use nuclei_data
    integer , intent(in) :: output
    
    integer :: i,a,j,b,   rt,idx,   l,k,iter
    integer :: O, V, OV, nrts ,ss, smax
    
    real(kind = 8) :: tol1=1E-6, tol2 = 1E-9, theta, tvecpn, resn, d
    
    integer , allocatable , dimension(:) :: rtidx
    
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO!,iajb_array2
    
    real (kind = 8) , allocatable , dimension(:,:) :: X,AX,M,P,S,AU,U
    
    real (kind = 8) , allocatable , dimension(:) :: SPA,W,Ws,tvec,tvecp,res
    
    real (kind = 8) :: pi=3.14159265359
    
        O = int(eltot/2)
        V = aotot - O
        OV = O*V
        nrts = 2
        ss = 1*nrts
        smax = int((OV-ss-1)/ss)
        print *, 'ov,smax:', ov, smax
        allocate(rtidx(ss))
        allocate(SPA(OV),tvec(OV),tvecp(OV))
        allocate(U(OV,ss),AU(OV,ss))
        allocate(X(OV,OV),AX(OV,OV))
        allocate(W(nrts))

        if (output == 1) then
            call smallblockheader('CLOSED SHELL CIS CALCULATION')
            call writelines(2)
            call onelineheader('Starting CIS routine')
            call writelines(2)
        endif
        


    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        
        allocate(J_MO(aotot,aotot,aotot,aotot))
        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
    !#######################################################################################################################!   
        
        do i=1,O
            do a=1,V
                SPA( V*(i-1)+a) = ( Eig(O+a) - Eig(i) ) + 2*J_MO(a+O,i,i,a+O) - J_MO(a+O,a+O,i,i)
            enddo
        enddo

        do i=1,ss
            rtidx(i) = minloc(SPA,dim=1,mask=(SPA > 0))
            SPA(rtidx(i)) = SPA(rtidx(i)) - 100000
            print *, rtidx(i)
        enddo

        do i=1,ss
            SPA(rtidx(i)) = SPA(rtidx(i)) + 100000
        enddo
        
        print *, SPA
        print *, rtidx
    !=======================================================================================================================!
    !                                           PERFORM CIS VIA THE DAVIDSON ALGORITHM                                      !
    !=======================================================================================================================!  
        allocate(P(OV,OV))


        !SETTING DAVIDSON MATRICES ZERO BEFORE CALCULATING  ROOT
        X = 0
        AX = 0 
        
        
        !RETRIEVING DIRECTION VECTORS FOR FIRST GUESS
        do rt=1,ss
            
            idx = rtidx(rt)
            
            X(idx,rt) = 1

            do i=1,O
                do a=1,V
                    
                    AX(V*(i-1)+a,rt) = AX(V*(i-1)+a,rt) + ( Eig(O+a) - Eig(i) ) *  X(V*(i-1)+a,rt)

                    do j=1,O
                        do b=1,V
                            
                            AX(V*(i-1)+a,rt) = AX(V*(i-1)+a,rt) + (2*J_MO(a+O,i,j,b+O) - J_MO(a+O,b+O,j,i)) * X(V*(j-1)+b,rt) 
                            
                        enddo
                    enddo
                enddo
            enddo
            
        enddo   

            
        !CARRYING OUT WITH THE DAVIDSON ALGORITHM FOR A SPECIFIC ROOT
        do iter=1,smax
                
            print *, 'iter:',iter
            
            k = iter*ss
            
            print *, 'k:',k

            allocate(M(k,k))
            allocate(S(k,k))
            allocate(Ws(k))

            
            M=0
            M = matmul(transpose(X(:,1:k)),AX(:,1:k))
            
            call diag(k,M,Ws,S)
                
            U = matmul(X(:,1:k),S(:,1:ss))
            AU = matmul(AX(:,1:k),S(:,1:ss))
               
            print *, 'Ws'
            print *, Ws
            resn = 0
            do rt=1,nrts
                res = AU(:,rt) - Ws(rt) * U(:,rt)
                do l=1,OV
                    resn = resn + res(l)**2
                enddo
            enddo
            resn = resn**0.5
                
            if (resn < tol2) then
                exit
            endif 
            
            do rt=1,ss
                
                theta = Ws(rt)
                res = AU(:,rt) - theta * U(:,rt)
                tvec = 0
                do l=1,OV
                    d = theta - SPA(l)
                    if (abs(d) < tol1) then
                        tvec(l) = 0
                    else
                        tvec(l) = res(l)/d
                    endif
                enddo
                
                P = 0
                P = -matmul(X(:,1:k+rt-1),transpose(X(:,1:k+rt-1)))
                do l=1,OV
                    P(l,l) = P(l,l) + 1
                enddo
                
                tvecp = matmul(P,tvec)
                
                tvecpn = 0
                do l=1,OV
                    tvecpn = tvecpn + tvecp(l)**2
                enddo
                tvecpn = tvecpn**0.5
                
                X(:,k+rt) = tvecp/tvecpn
                
                do i=1,O
                    do a=1,V
                    
                        AX(V*(i-1)+a,k+rt) = AX(V*(i-1)+a,k+rt) + ( Eig(O+a) - Eig(i) ) *  X(V*(i-1)+a,k+rt)

                        do j=1,O
                            do b=1,V
                            
                                AX(V*(i-1)+a,k+rt) = AX(V*(i-1)+a,k+rt) + (2*J_MO(a+O,i,j,b+O) - J_MO(a+O,b+O,j,i)) * X(V*(j-1)+b,k+rt) 
                            
                            enddo
                        enddo
                    enddo
                enddo
                
            enddo
                

        
                
            deallocate(M,S,Ws)
        enddo

            
            

        deallocate(J_MO)

    !#######################################################################################################################!
        

        if (output == 1) then
        
            call writelines(2)
            call onelineheader("Singlet CIS Excitation Energies")
            
            print *, 'theta:', theta*27
          
            call onelineheader("CISES SPECTRA Calculated!")
            write(77,*)''
            call writelines(2)
            call twolinesfooter("CISES CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine CIS_multiple_davidson_calc

    
    
    
    
    
    !#######################################################################################################################!
    
    
    
    
    
    
    
    subroutine closed_shell_TDHF_calc(output)
    use nuclei_data
    ! Written by Chima Chibueze
    ! Solving the TD-HF equations
    ! Literature: Dreeuw, A.; Head-Gordon, M. Chem. Rev. 2005, 105, 4009?4037
    
    integer , intent(in) :: output
    
    integer :: i,j,k,l,p,q,r,s,t,detno,d,a,b,O,V
    real (kind = 8) :: pi=3.14159265359
    real (kind = 8) , dimension(3) :: DP,edtm
    real (kind = 8) :: elli
    real (kind = 8) , allocatable , dimension(:,:) :: DPx,DPy,DPz,DPx_MO,DPy_MO,DPz_MO
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: J_MO!,iajb_array2
    real (kind = 8) , allocatable , dimension(:,:) :: Amat,Bmat,Dmat,Zmat,DZvec
    real (kind = 8) , allocatable , dimension(:,:) :: AmBmat, ApBmat,orbi
    real (kind = 8) , allocatable , dimension(:) :: W
    
        if (output == 1) then
            call smallblockheader('CLOSED SHELL TD-HF CALCULATION')
            call writelines(2)
            call onelineheader('Starting TD-HF routine')
            call writelines(2)
        endif
        
        O = int(eltot/2)
        V = aotot - O
        detno = O*V
        allocate(J_MO(aotot,aotot,aotot,aotot))


    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                          INTEGRAL TRANSFORMATION (AO ==> MO)                                          !
    !=======================================================================================================================!
        

        call four_index_transform(aotot,Cmo,J_ee,J_MO)        !J_AO ==> J_MO MATRIX IN MO BASIS
    !#######################################################################################################################!
        
        
        
        
        
        
        


    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                      BUILDING A,B AND Z MATRICES AND PERFORM CIS                                      !
    !=======================================================================================================================!  
        
        allocate(Amat(detno,detno),Bmat(detno,detno),Dmat(detno,detno),Zmat(detno,detno))
        allocate(AmBmat(detno,detno),ApBmat(detno,detno),orbi(detno,detno))
        allocate(W(detno)) 
        Amat = 0
        Bmat = 0
        ApBmat = 0
        AmBmat = 0
        Dmat = 0
        orbi = 0
        do i=1,O
            do a=1,V
                Amat( V*(i-1)+a , V*(i-1)+a ) = Amat( V*(i-1)+a , V*(i-1)+a ) + ( Eig(O+a) - Eig(i) )
                orbi( V*(i-1)+a , V*(i-1)+a ) = orbi( V*(i-1)+a , V*(i-1)+a ) + ( Eig(O+a) - Eig(i) )
                
                do j=1,O
                    do b=1,V
                        Amat(V*(i-1)+a,V*(j-1)+b) = Amat(V*(i-1)+a,V*(j-1)+b) + 2*J_MO(a+O,i,j,b+O) - J_MO(a+O,b+O,j,i)
                        Bmat(V*(i-1)+a,V*(j-1)+b) = Bmat(V*(i-1)+a,V*(j-1)+b) + 2*J_MO(a+O,i,b+O,j) - J_MO(a+O,j,b+O,i)
                    enddo
                enddo
                
            enddo
        enddo

        deallocate(J_MO)
        
        AmBmat = Amat - Bmat
        ApBmat = Amat + Bmat
!        deallocate(Amat,Bmat)

        
        Dmat = matmul(AmBmat,ApBmat)
        !deallocate(AmBmat,ApBmat)
        
        call gdiag(detno,AmBmat,ApBmat,W,Zmat)

        do i=1,detno
            W(i) = (W(i))**0.5*27
        enddo
        

    !#######################################################################################################################!
        
        
        
        
        
        
        
        
        
        
        
        
    !#######################################################################################################################!
        
    !=======================================================================================================================!
    !                                           OBTAINING TRANSITION PROBABILITIES                                          !
    !=======================================================================================================================!
        allocate(DPx(aotot,aotot),DPy(aotot,aotot),DPz(aotot,aotot))
        allocate(DPx_MO(aotot,aotot),DPy_MO(aotot,aotot),DPz_MO(aotot,aotot))   
        call DP_calc(DPx,DPy,DPz) 
        DPx_MO = matmul( transpose(Cmo), matmul(DPx,Cmo) ) !DP MATRIX IN MO BASIS
        DPy_MO = matmul( transpose(Cmo), matmul(DPy,Cmo) ) !DP MATRIX IN MO BASIS
        DPz_MO = matmul( transpose(Cmo), matmul(DPz,Cmo) ) !DP MATRIX IN MO BASIS
        deallocate(DPx,DPy,DPz)
        DP = 0
        do i=1,O
            DP(1) =  DP(1) + 2 * DPx_MO(i,i)
            DP(2) =  DP(2) + 2 * DPy_MO(i,i)
            DP(3) =  DP(3) + 2 * DPz_MO(i,i)
        enddo
    
        DP = -DP
    
        do j=1,nuctot
            DP =  DP + nuclei(j)%charge * nuclei(j)%position
        enddo
        
        write(77,*) 'dipole moment:', ( (DP(1)**2 + DP(2)**2 + DP(3)**2)**0.5 ) * 2.541746473, 'D'

        
         elli = 0
         do j = 1,detno
            edtm = 0
            do i=1,O
                do a=1,V
                    edtm(1) = edtm(1) + Zmat( V*(i-1)+a , j ) * DPx_MO(i,a+O)
                    edtm(2) = edtm(2) + Zmat( V*(i-1)+a , j ) * DPy_MO(i,a+O)
                    edtm(3) = edtm(3) + Zmat( V*(i-1)+a , j ) * DPz_MO(i,a+O)
                enddo
            enddo
            edtm = edtm * 2
            write(77,*) 'f:',  (2*W(j) / 3) * ( (edtm(1)**2 + edtm(2)**2 + edtm(3)**2) ) , 'at' ,( (137*2*pi)/W(j) ) * (5.29177210903 / 100) , 'nm'! , W_CI(j) * 27, 'eV)' !* 2.541746473, 'D'
            elli = elli + (2*W(j) / 3) * ( (edtm(1)**2 + edtm(2)**2 + edtm(3)**2) )
        enddo
        write(77,*) 'sum of oscillator strengths:', elli
        
    !#######################################################################################################################!
        
        
        
        
        
        
        
        if (output == 1) then
            call writelines(2)
            !call onelineheader("delta matrix")
            !call matprint(detno,orbi)
            !call onelineheader("A matrix")
            !call matprint(detno,Amat-orbi)
            !call onelineheader("B matrix")
            !call matprint(detno,Bmat)
            !call onelineheader("D matrix")
            !call matprint(detno,Dmat)
            !call writelines(2)
        
            call writelines(2)
            call onelineheader("TD-HF energies")
            call vecprint(detno,W)
            call writelines(2)
            call onelineheader("TD-HF coefficient matrix")
            call matprint(detno,Zmat)
            call writelines(2)
          
            call onelineheader("TD-HF SPECTRA Calculated!")
            write(77,*)''
            call writelines(2)
            call twolinesfooter("TD-HF CALCULATION DONE")
            call writelines(8)
        endif
        
    end subroutine closed_shell_TDHF_calc
    
end module CI_correlation
!#######################################################################################################################!