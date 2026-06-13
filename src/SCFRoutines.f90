!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                      SCF ROUTINES                                                     !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
!#######################################################################################################################!
module SCF_routines
use parallel
!$ use omp_lib
#if USE_MPI
use MPI
#endif
use Timing
use result_data
use integral_tensors
use molecular_data
use calculation_data
use Print_module
use SCF_matrix_builder
use matvec_operations
implicit none
save
contains
    
    !=======================================================================================================================!
    !                                               RESTRICTED HARTREE-FOCK                                                 !
    !=======================================================================================================================!
    
    subroutine HFC(output,Eout)
    integer , intent(in) :: output    
    real (kind = 8) , optional , intent(out) :: Eout
    integer :: a,i,p,q, ld
    real (kind = 8) :: Eel_prev,ts,tf
    real (kind = 8) , dimension(aotot) :: lab
    real (kind = 8) , dimension(aotot,aotot) :: S,V,VT,VTS,ss,ss_12,V_ss_12,U,UT
    real (kind = 8) , dimension(aotot,aotot) :: T,Vne,J,K,G,F,F_,C_,D,UTF
    logical :: do_ld = .false.

        Cmo = 0.0
        call Vnn_calc(Vnn)
        call S_calc(S)
        

        call T_calc(T)
        call Vne_calc(Vne)
        H = T + Vne

        call diag(aotot,S,lab,V)
        VT = transpose(V)
        
        ss_12  = 0.0
        !lin. dep. check!
        ld = 0

        do p=1,aotot
            if (do_ld) then
            if (lab(p) < 1E-5) then
                ld = ld + 1
            else
                ss_12(p,p) = 1/((lab(p)**0.5))
            endif
            else
                ss_12(p,p) = 1/((lab(p)**0.5))
            endif
        enddo
        !print *, lab 
        !print *, ld
        !print *, ''
    
        
        !SYMMETRIC LOWDIN ORTHOGONALIZATION
        V_ss_12 = matmul(V,ss_12)
        U = matmul(V_ss_12,VT)
        UT = transpose(U)   
        !print *, U
        

        call timer(ts)
        call J_ee_calc()
        call timer(tf)

        call printtext('Done calculate two-electron integrals in (s):', r=(tf - ts))

        ! bare nuclei Hamiltonian as initial Fock matrix
        
        F = H + 0
        
        UTF = matmul(UT,F)
        F_ = matmul(UTF,U)
        
        !prevent lin. dep.
        
        do p=1,ld
            F_(p,1:aotot) = 0
            F_(1:aotot,p) = 0
            F_(p,p) = 100*p**2
        enddo
        
        call matprint(aotot,F_)
        call diag(aotot,F_,Eig,C_)

        Cmo = matmul(U,C_)
        D = 0
        do i=1,int(eltot/2)
            do p=1,aotot
                do q=1,aotot
                    D(q,p) = D(q,p) + Cmo(p,i) * Cmo(q,i)
                enddo
            enddo
        enddo
        Eel_prev = 0

        if (output==1) then
            call smallblockheader('SINGLE-POINT HARTREE-FOCK CALCULATION')
            call writelines(2)
            call onelineheader('Starting Hartree-Fock (HF) routine')
            write(77,*)''
        endif
    
    
        do a=1,max_iter
            
            call J_calc(D,J)
            call K_calc(D,K)
            G = 2*J - K
            F = H + G
            UTF = matmul(UT,F)
            F_ = matmul(UTF,U)

            !prevent lin. dep.
            do p=1,ld
                F_(p,1:aotot) = 0
                F_(1:aotot,p) = 0
                F_(p,p) = 100*p**2
            enddo
        
            call diag(aotot,F_,Eig,C_)
            Cmo = matmul(U,C_)
            D = 0
            do i=1,int(eltot/2)
                do p=1,aotot
                    do q=1,aotot
                        D(q,p) = D(q,p) + Cmo(p,i) * Cmo(q,i)
                    enddo
                enddo
            enddo
            Eel = 0
            do p=1,aotot
                do q=1,aotot
                    Eel = Eel + D(q,p)*(H(p,q)+F(p,q))
                enddo
            enddo
        
            if (output==1) then
                write(77,*) ''
                write(77,*) "Cycle",a
                write(77,*) "Electronic Energy:", Eel
                write(77,*)''        
            endif
        
            if (abs(Eel_prev - Eel) < eps) then
                
                EHF = Eel + Vnn
                if (present(Eout)) then
                    Eout = EHF
                endif
                
                if (output == 1) then
                    write(77,*)''
                    call onelineheader("Energy Converged!")
                    write(77,*)''
                    write(77,*) "SCF steps:",a
                    write(77,*)''
                    write(77,*) "Orbital Energies:"
                    write(77,10), (Eig(i) , i=1,aotot)
        10          format(12f10.3)
                    write(77,*)''
                    write(77,*) "Orbital Coefficients:"
                    call matprint(aotot,Cmo)
                    write(77,*)''
                    write(77,*) "Electronic Energy:",Eel
                    write(77,*)''
                    write(77,*) "Nuclear Energy:",Vnn
                    write(77,*)''
                    write(77,*) "Total Energy:",EHF
                    call writelines(2)
                    write(77,*) 
                    call twolinesfooter("HARTREE-FOCK CALCULATION DONE")
                    call writelines(7)
                endif
            
            
                exit
            
            elseif (a == max_iter - 1) then
                write(77,*) ''
                write(77,*) "Warning"
                write(77,*) "SFC not Converged; Max. SCF steps reached"
                write(77,*) "DONE"
            
            endif
            Eel_prev = Eel
        enddo

    end subroutine HFC
    
    !=======================================================================================================================!
    !                                             UNRESTRICTED HARTREE-FOCK                                                 !
    !=======================================================================================================================!
    
    
    subroutine UHFC(output)
    integer , intent(in) :: output
    character (len=30) :: state
    integer :: a,i,k,l,p,q
    real (kind = 8) :: Eel_prev,ts,tf,sz,s2_theoretical,s2_actual,s2aux
    real (kind = 8) , dimension(aotot) :: lab
    real (kind = 8) , dimension(aotot,aotot) :: S,V,VT,VTS,ss,ss_12,V_ss_12,U,UT
    real (kind = 8) , dimension(aotot,aotot) :: T,Vne,J,Ka,Kb,dK,Ga,Gb,Fa,Fb,Fa_,Fb_
    real (kind = 8) , dimension(aotot,aotot) :: Da,Db,D,Ca_,Cb_,UTFa,UTFb,Sab

        
    ! ONLY ACCOUNT FOR SINGLETS, DOUBLETS, TRIPLETS AND QUARTETS FOR NOW!!!
        if (multiplicity == 2) then
            state = 'doublet'
            na = ceiling(real(eltot)/2)
            nb = floor(real(eltot)/2)
        else if (multiplicity == 3) then
            state = 'triplet'
            na = eltot/2 + 1
            nb = eltot/2 - 1
        else if (multiplicity == 4) then
            state = 'quartet'
            na = ceiling(real(eltot)/2) + 1
            nb = floor(real(eltot)/2) - 1
        else
            state = 'singlet'
            na = eltot/2
            nb = eltot/2
        endif
        
        sz = (real(na) - real(nb))/2
        s2_theoretical = sz * (sz + 1)
        s2_actual = s2_theoretical
        
        Ca = 0
        Cb = 0
        
        call Vnn_calc(Vnn)
        call S_calc(S)
        call T_calc(T)
        call Vne_calc(Vne)
        H = T + Vne
        
        call diag(aotot,S,lab,V)
        VT = transpose(V)

        ss_12  = 0.0
        do a=1,aotot
            ss_12(a,a) = 1/((lab(a)**0.5))
        enddo
    
        !SYMMETRIC LOWDIN ORTHOGONALIZATION
        V_ss_12 = matmul(V,ss_12)
        U = matmul(V_ss_12,VT)
        UT = transpose(U)   

        call J_ee_calc()
        
        Fa = H
        Fb = H 
        
        UTFa = matmul(UT,Fa)
        UTFb = matmul(UT,Fb)
        
        Fa_ = matmul(UTFa,U)
        Fb_ = matmul(UTFb,U)
         
        call diag(aotot,Fa_,Eiga,Ca_)
        call diag(aotot,Fb_,Eigb,Cb_)
        
        Ca = matmul(U,Ca_)
        Cb = matmul(U,Cb_)
        
        Eel_prev = 0
    
        if (output==1) then
            call smallblockheader('SINGLE-POINT UNRESTRICTED HARTREE-FOCK CALCULATION')
            call writelines(2)
            call onelineheader('Starting Unrestricted Hartree-Fock (UHF) routine')
            write(77,*)''
            write(77,13)'Carrying out UHF for a ',state,'state, with multiplicity ',multiplicity,'.'
13          format(A23,A7,A26,I1,A1)          
            write(77,*)'Number of alpha electrons:',na
            write(77,*)'Number of beta electrons: ',nb         
        endif
    
    
        do a=1,max_iter
            Da = 0
            Db = 0
            D = 0
            do p=1,aotot
                do q=1,aotot
                    do i=1,na
                        Da(q,p) = Da(q,p) + Ca(p,i)*Ca(q,i)
                    enddo
                    do i = 1,nb
                        Db(q,p) = Db(q,p) + Cb(p,i)*Cb(q,i)
                    enddo
                enddo
            enddo
            D = Da + Db
            call Juhf_calc(D,J)
            call Kuhf_calc(Da,Ka)
            call Kuhf_calc(Db,Kb)
            dK = (Ka - Kb) / 2
            
            Ga = J - Ka
            Gb = J - Kb
        
            Fa = H + Ga
            Fb = H + Gb
            
            Fa(1:nb,na+1:aotot) = Fa(1:nb,na+1:aotot) + dK(1:nb,na+1:aotot)
            Fa(na+1:aotot,1:nb) = Fa(na+1:aotot,1:nb) + dK(na+1:aotot,1:nb)
            
            Fb(1:nb,na+1:aotot) = Fb(1:nb,na+1:aotot) - dK(1:nb,na+1:aotot)
            Fb(na+1:aotot,1:nb) = Fb(na+1:aotot,1:nb) - dK(na+1:aotot,1:nb)
            UTFa = matmul(UT,Fa)
            UTFb = matmul(UT,Fb)
        
            Fa_ = matmul(UTFa,U)
            Fb_ = matmul(UTFb,U)
         
            call diag(aotot,Fa_,Eiga,Ca_)
            call diag(aotot,Fb_,Eigb,Cb_)
        
            Ca = matmul(U,Ca_)
            Cb = matmul(U,Cb_)
            
            Eel = 0
            do p=1,aotot
                do q=1,aotot
                    Eel = Eel + Da(q,p)*(H(p,q)+Fa(p,q))
                    Eel = Eel + Db(q,p)*(H(p,q)+Fb(p,q))
                enddo
            enddo

            Eel=Eel/2
        
            if (output==1) then
                write(77,*) ''
                write(77,*) "Cycle",a
                write(77,*) "Electronic Energy:", Eel
                write(77,*)''        
            endif
        
            !if (a > 48) then
            if (abs(Eel_prev - Eel) < eps) then
                
                Sab = matmul(transpose(Ca),matmul(S,Cb))
                s2_actual = s2_actual + nb
                do k = 1, na
                    do l = 1, nb
                        s2_actual = s2_actual - Sab(k,l) * Sab(k,l)
                        !s2aux = 0
                        !do p = 1, aotot
                        !    do q = 1, aotot
                        !        s2aux = s2aux + Ca(p,k) * Cb(q,l) * S(p,q)
                        !    enddo
                        !enddo
                        !s2_actual = s2_actual - s2aux**2
                    enddo
                enddo
                
                !call matprint(aotot, matmul( transpose(Sab),Sab ) )
                
                EHF = Eel + Vnn
                if (output == 1) then
                    write(77,*)''
                    call onelineheader("Energy Converged!")
                    write(77,*)''
                    write(77,*) "SCF steps:",a
                    write(77,*)''
                    write(77,*) "Alpha Orbital Energies:"
                    write(77,10), (Eiga(i) , i=1,aotot)
        10          format(12f10.3)
                    write(77,*) "Beta Orbital Energies:"
                    write(77,16), (Eigb(i) , i=1,aotot)
        16          format(12f10.3)
                    write(77,*)''
                    write(77,*) "Alpha Orbital Coefficients:"
                    call matprint(aotot,Ca)
                    write(77,*)''
                    write(77,*) "Beta 0rbital Coefficients:"
                    call matprint(aotot,Cb)
                    write(77,*)''
                    write(77,*) "Electronic Energy:",Eel
                    write(77,*)''
                    write(77,*) "Nuclear Energy:",Vnn
                    write(77,*)''
                    write(77,*) "Total Energy:",EHF
                    write(77,*)
                    write(77,*) "z-component of spin operator S_z:", sz
                    write(77,*)
                    write(77,*) "Theoretical expectation value of the squared spin operator, <S^2>_exact:", s2_theoretical
                    write(77,*)
                    write(77,*) "Actual <S^2>_UHF value for the obtained UHF wavefunction:", s2_actual
                    write(77,*)
                    write(77,*) "Difference between <S^2>_UHF and <S^2>_exact:", s2_actual - s2_theoretical
                    call writelines(2)
                    write(77,*) 
                    call twolinesfooter("UNRESTRICTED HARTREE-FOCK CALCULATION DONE")
                    call writelines(7)
                endif
            
            
                exit
            
            elseif (a == max_iter - 1) then
                write(77,*) ''
                write(77,*) "Warning"
                write(77,*) "SFC not Converged; Max. SCF steps reached"
                write(77,*) "DONE"
            
            endif
            Eel_prev = Eel
        enddo

    end subroutine UHFC

    
    

    !=======================================================================================================================!
    !                                                  GENERAL HARTREE-FOCK                                                 !
    !=======================================================================================================================!
    
    
    subroutine GHFC(output)
    integer , intent(in) :: output
    character (len=30) :: state
    integer :: a,i,k,l,p,q
    real (kind = 8) :: Eel_prev,ts,tf,sz_theoretical,sz,s2_theoretical,s2,Naa,Nbb,S2rohf,dS2uhf,dS2z,dS2xy
    real (kind = 8) , dimension(aotot*2) :: lab
    real (kind = 8) , dimension(eltot,eltot) :: Saa,Sbb,Sab
    real (kind = 8) , dimension(aotot,aotot) :: S_block,V,VT,VTS,ss,ss_12,V_ss_12,U_block
    real (kind = 8) , dimension(aotot,aotot) :: T,Vne,J,Kaa,Kbb,Kab,Kba,Faa,Fbb,Daa,Dab,Dba,Dbb,Dtot
    real (kind = 8) , dimension(aotot*2,aotot*2) :: U,UT,F,F_,C_,UTF,S,D,DF
    real (kind = 8) , dimension(aotot,eltot) :: Caa,Cbb

    
    ! ONLY ACCOUNT FOR SINGLETS, DOUBLETS, TRIPLETS AND QUARTETS FOR NOW!!!
        if (multiplicity == 2) then
            state = 'doublet'
            na = ceiling(real(eltot)/2)
            nb = floor(real(eltot)/2)
        else if (multiplicity == 3) then
            state = 'triplet'
            na = eltot/2 + 1
            nb = eltot/2 - 1
        else if (multiplicity == 4) then
            state = 'quartet'
            na = ceiling(real(eltot)/2) + 1
            nb = floor(real(eltot)/2) - 1
        else
            state = 'singlet'
            na = eltot/2
            nb = eltot/2
        endif
        
        sz_theoretical = (real(na) - real(nb))/2
        s2_theoretical = sz_theoretical * (sz_theoretical + 1)
        
        Cmo=0
        
        call Vnn_calc(Vnn)
        call S_calc(S_block)
        call T_calc(T)
        call Vne_calc(Vne)
        H = T + Vne
        call diag(aotot,S_block,lab,V)
        VT = transpose(V)

        ss_12  = 0.0
        do a=1,aotot
            ss_12(a,a) = 1/((lab(a)**0.5))
        enddo
    
        !SYMMETRIC LOWDIN ORTHOGONALIZATION
        V_ss_12 = matmul(V,ss_12)
        U_block = matmul(V_ss_12,VT)
        U = 0
        U(1:aotot,1:aotot) = U_block
        U((aotot+1):(aotot*2),(aotot+1):(aotot*2)) = U_block
        UT = transpose(U)
            
        S = 0
        S(1:aotot,1:aotot) = S_block
        S(aotot+1:2*aotot,aotot+1:2*aotot) = S_block

        call J_ee_calc()

        F = 0
        !H�ckel bare nucleus Hamiltonian
        F(1:aotot,1:aotot) = H !+ J                         - Kaa
        F((aotot+1):(aotot*2),(aotot+1):(aotot*2)) = H !+ J - Kbb

        !Trying to mix some alpha and beta spins
        Kab = 0.0
        Kba = Kab
        F(1:aotot,(aotot+1):(aotot*2)) =                  0 - Kab
        F((aotot+1):(aotot*2),1:aotot) =                  0 - Kba

        UTF = matmul(UT,F)
        F_ = matmul(UTF,U)

        call diag(2*aotot,F_,Eig,C_)

        Cmo = matmul(U,C_)
        Caa = Cmo(1:aotot,1:eltot)
        Cbb = Cmo((aotot+1):(aotot*2),1:eltot)
        
        Daa = 0
        Dab = 0
        Dba = 0
        Dbb = 0   
        do p = 1,aotot
            do q = 1,aotot
                do i=1,eltot
                    Daa(q,p) = Daa(q,p) + Caa(p,i)*Caa(q,i)
                    Dab(q,p) = Dab(q,p) + Cbb(p,i)*Caa(q,i)
                    Dba(q,p) = Dba(q,p) + Caa(p,i)*Cbb(q,i)
                    Dbb(q,p) = Dbb(q,p) + Cbb(p,i)*Cbb(q,i)                        
                enddo
            enddo
        enddo  
        Dtot = Daa + Dbb
        Eel_prev = 0
    
        if (output==1) then
            call smallblockheader('SINGLE-POINT GENERAL HARTREE-FOCK CALCULATION')
            call writelines(2)
            call onelineheader('Starting General Hartree-Fock (GHF) routine')
            write(77,*)''
            write(77,13)'Carrying out GHF for a ',state,'state, with multiplicity ',multiplicity,'.'
13          format(A23,A7,A26,I1,A1)          
            write(77,*)'Number of alpha electrons:',na
            write(77,*)'Number of beta electrons: ',nb         
        endif
    
    
        do a=1,max_iter
            call Jghf_calc(Dtot,J)
            
            call Kghf_calc(Daa,Kaa)
            call Kghf_calc(Dab,Kab)
            call Kghf_calc(Dba,Kba)
            call Kghf_calc(Dbb,Kbb)

            F = 0
            F(1:aotot,1:aotot) = H + J                         - Kaa
            F(1:aotot,(aotot+1):(aotot*2)) =                 0 - Kab
            F((aotot+1):(aotot*2),1:aotot) =                 0 - Kba
            F((aotot+1):(aotot*2),(aotot+1):(aotot*2)) = H + J - Kbb
            
            UTF = matmul(UT,F)
            F_ = matmul(UTF,U)
            
            call diag(aotot*2,F_,Eig,C_)
        
            Cmo = matmul(U,C_)
            Caa = Cmo(1:aotot,1:eltot)
            Cbb = Cmo((aotot+1):(aotot*2),1:eltot)
            
            Daa = 0
            Dab = 0
            Dba = 0
            Dbb = 0
            
            do p = 1,aotot
                do q = 1,aotot
                    do i=1,eltot
                        Daa(q,p) = Daa(q,p) + Caa(p,i)*Caa(q,i)
                        Dab(q,p) = Dab(q,p) + Cbb(p,i)*Caa(q,i)
                        Dba(q,p) = Dba(q,p) + Caa(p,i)*Cbb(q,i)
                        Dbb(q,p) = Dbb(q,p) + Cbb(p,i)*Cbb(q,i)                        
                    enddo
                enddo
            enddo  
            D(1:aotot,1:aotot) =                         Daa
            D(1:aotot,(aotot+1):(aotot*2)) =             Dab
            D((aotot+1):(aotot*2),1:aotot) =             Dba
            D((aotot+1):(aotot*2),(aotot+1):(aotot*2)) = Dbb
            Dtot = Daa + Dbb
            Eel = 0
            F(1:aotot,1:aotot) = F(1:aotot,1:aotot) + H
            F((aotot+1):(aotot*2),(aotot+1):(aotot*2)) = F((aotot+1):(aotot*2),(aotot+1):(aotot*2)) + H
            DF = matmul(D,F)

            do i = 1,aotot*2
                Eel = Eel + DF(i,i)
            enddo
            Eel=Eel/2
        
            if (output==1) then
                write(77,*) ''
                write(77,*) "Cycle",a
                write(77,*) "Electronic Energy:", Eel
                write(77,*)''        
            endif
        
            !if (a > 48) then
            if (abs(Eel_prev - Eel) < eps) then
                Saa = 0
                Sbb = 0
                Sab = 0
                Saa = matmul(transpose(Caa),matmul(S,Caa))
                Sbb = matmul(transpose(Cbb),matmul(S,Cbb))   
                Sab = matmul(transpose(Caa),matmul(S,Cbb))
                
                Naa = 0
                Nbb = 0
                do i=1,eltot
                    Naa = Naa + Saa(i,i)
                    Nbb = Nbb + Sbb(i,i)
                enddo
                sz = (Naa - Nbb)/2
                S2rohf = (Naa - Nbb)*(Naa - Nbb + 2)/4
                dS2uhf = Nbb
                dS2z = eltot
                dS2xy = 0
                do p=1,eltot
                    dS2xy = dS2xy + Sab(p,p)
                    do q=1,eltot
                        dS2uhf = dS2uhf - Sab(p,q)**2
                        dS2z =   dS2z   - (Saa(p,q) - Sbb(p,q))**2
                    enddo
                enddo
                dS2z = dS2z / 4
                dS2xy = dS2xy**2
                s2 = S2rohf + dS2uhf + dS2z + dS2xy
                EHF = Eel + Vnn
                if (output == 1) then
                    write(77,*)''
                    call onelineheader("Energy Converged!")
                    write(77,*)''
                    write(77,*) "SCF steps:",a
                    write(77,*)''
                    write(77,*) "Orbital Energies:"
                    write(77,10), (Eig(i) , i=1,aotot*2)
        10          format(12f10.3)
                    write(77,*)''
                    write(77,*) "Orbital Coefficients:"
                    call matprint(aotot*2,Cmo)
                    write(77,*)''
                    write(77,*) "Electronic Energy:",Eel
                    write(77,*)''
                    write(77,*) "Nuclear Energy:",Vnn
                    write(77,*)''
                    write(77,*) "Total Energy:",EHF
                    write(77,*)
                    write(77,*) "Theoretical number of alpha electrons, Na and Nb:", na, nb
                    write(77,*)
                    write(77,*) "Actual number of alpha and beta electrons, Na and Nb:", Naa, Nbb
                    write(77,*)
                    write(77,*) "Theoretical z-component of spin operator S_z:", sz_theoretical
                    write(77,*)
                    write(77,*) "Theoretical expectation value of the squared spin operator, <S^2>_exact:", s2_theoretical
                    write(77,*)
                    write(77,*) "Actual z-component of spin operator S_z, <S^2>:", sz
                    write(77,*)
                    write(77,*) "Actual squared spin expectation value, <S^2>:", s2
                    write(77,*)
                    write(77,*) "Individual contributions to <S^2> = ROHF + DODS + Z + XY:"
                    write(77,*)
                    write(77,*) "ROHF contribution <S^2>_ROHF = (Na - Nb)(Na - Nb + 2)/4:", S2rohf
                    write(77,*)
                    write(77,*) "DODS contribution <S^2>_DODS = Nb - sum_ij |S^{ab}_{ij}|^2:", dS2uhf
                    write(77,*)
                    write(77,*) "Z-contamination contribution <S^2>_Z = (Ne - sum_ij |S^{aa}_ij - S^{bb}_ij|^2)/4:", dS2z
                    write(77,*)
                    write(77,*) "XY-contamination contribution <S^2>_XY = |sum_i S^{ab}_{ii}|^2:", dS2xy
                    call writelines(2)
                    write(77,*) 
                    call twolinesfooter("GENERAL HARTREE-FOCK CALCULATION DONE")
                    call writelines(7)
                endif
            
            
                exit
            
            elseif (a == max_iter - 1) then
                write(77,*) ''
                write(77,*) "Warning"
                write(77,*) "SFC not Converged; Max. SCF steps reached"
                write(77,*) "DONE"
            
            endif
            Eel_prev = Eel
        enddo

    end subroutine GHFC
    
end module SCF_routines     
!#######################################################################################################################!   