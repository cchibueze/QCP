!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                            SCF MATRICES BUILDING ROUTINES                                             !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
        
!#######################################################################################################################!     
module SCF_matrix_builder
use integrals
use molecular_data
use ao_data
implicit none
save
private :: S_ij_calc,T_ij_calc,Vne_ij_calc,J_ij_calc,K_ij_calc,D_ij_calc,DP_ij_calc,MV_ij_calc
contains

! First the nuclear potential energy has to be calculated, independently from the electronic one (Born-Oppenheimmer)
! This is just a i<j sum of coulomb interaction between the nuclei
! Although treated classicaly, this is relatively accurate since the nuclei are much smaller then the electrons.
    
    
    !=======================================================================================================================!
    !                                                NUCLEAR POTENTIAL CALCULATOR                                           !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Calculates nuclear potential in a classical fashion.                                                                  !
    ! Returns a float.                                                                                                      !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!
    
    subroutine Vnn_calc(Vnn)
    use nuclei_data
    real (kind = 8) :: Vnn,d,q1,q2
    real (kind = 8) , dimension(1:3) :: r1,r2
    integer :: i,j
        Vnn = 0
        do i=1,nuctot
            do j=i+1,nuctot
                r1 = nuclei(i)%position
                r2 = nuclei(j)%position
                q1 = nuclei(i)%charge
                q2 = nuclei(j)%charge
                d = ((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 + (r1(3)-r2(3))**2)**0.5
                Vnn = Vnn + (q1*q2)/d
            enddo
        enddo
    end subroutine Vnn_calc
    
    

    
    ! The following subroutines together generate the Fock matrix, needed to find the optimized mo-coefficients and electronic energies.
    


    
    !=======================================================================================================================!
    !                                               OVERLAP MATRIX ELEMENT CALCULATOR                                       !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates overlap between two contracted Gaussians                                                                    !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!   
    
    subroutine S_ij_calc(i,j,s_ij)
    integer , intent(in) :: i,j
    integer :: k,l
    real (kind = 8) :: s_ij!,overlap
    type (ao) :: ao1,ao2
        ao1 = aos(i)
        ao2 = aos(j)
        s_ij = 0.0
        do k = 1,ao1%num_cont
            do l = 1,ao2%num_cont
                s_ij = s_ij + ao1%coef(k)*ao2%coef(l)* &
                overlap(ao1%exp(k),ao1%shell,ao1%position, &
                        ao2%exp(l),ao2%shell,ao2%position)
            enddo
        enddo
    end subroutine

    

    
    
    

    !=======================================================================================================================!
    !                                               OVERLAP MATRIX  CALCULATOR                                              !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                               Generates the overlap matrix S                                          !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!
    
    subroutine S_calc(S)
    real (kind = 8) :: s_ij
    real (kind = 8) , dimension(aotot,aotot) :: S
    integer :: i,j
        do i = 1,aotot
            do j = i,aotot
                call S_ij_calc(i,j,s_ij)
                S(i,j) = s_ij
                S(j,i) = s_ij
            enddo
        enddo    
    end subroutine S_calc
    
    

    
    
    !=======================================================================================================================!
    !                                              DIPOLE MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates dipole moment for two contracted Gaussians                                                                  !
    ! Returns float.                                                                                                        !                                                      
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine DP_ij_calc(i,j,dpx_ij,dpy_ij,dpz_ij)
    integer , intent(in) :: i,j
    integer :: k,l
    real (kind = 8) :: dpx,dpy,dpz,dpx_ij,dpy_ij,dpz_ij
    type (ao) :: ao1,ao2
    ao1 = aos(i)
    ao2 = aos(j)
        dpx_ij = 0
        dpy_ij = 0
        dpz_ij = 0
        
        do k = 1,ao1%num_cont
            do l = 1,ao2%num_cont
                call dipole(ao1%exp(k),ao1%shell,ao1%position, &
                        ao2%exp(l),ao2%shell,ao2%position,dpx,dpy,dpz)
                
                dpx_ij = dpx_ij + ao1%coef(k)*ao2%coef(l)* dpx
                dpy_ij = dpy_ij + ao1%coef(k)*ao2%coef(l)* dpy
                dpz_ij = dpz_ij + ao1%coef(k)*ao2%coef(l)* dpz
                
            enddo
        enddo
    end subroutine DP_ij_calc


    
    !=======================================================================================================================!
    !                                               DIPOLE MATRIX  CALCULATOR                                               !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                               Generates the dipole matrix D                                           !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine DP_calc(DPx,DPy,DPz)
    real (kind = 8) :: dpx_ij,dpy_ij,dpz_ij
    real (kind = 8) , dimension(aotot,aotot) :: DPx,DPy,DPz
    integer :: i,j
        do i = 1,aotot
            do j = i,aotot
                call DP_ij_calc(i,j,dpx_ij,dpy_ij,dpz_ij)
                DPx(i,j) = dpx_ij
                DPx(j,i) = dpx_ij
                DPy(i,j) = dpy_ij
                DPy(j,i) = dpy_ij
                DPz(i,j) = dpz_ij
                DPz(j,i) = dpz_ij
            enddo
        enddo    
    end subroutine DP_calc
    
    
    
    
    
    
    
    !=======================================================================================================================!
    !                                               KINETIC MATRIX ELEMENT CALCULATOR                                       !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates kinetic energy for two contracted Gaussians                                                                 !
    ! Returns float.                                                                                                        !                                                      
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine T_ij_calc(i,j,t_ij)
    integer , intent(in) :: i,j
    integer :: k,l
    real (kind = 8) :: t_ij!,kinetic
    type (ao) :: ao1,ao2
        ao1 = aos(i)
        ao2 = aos(j)
        t_ij = 0.0
        do k = 1,ao1%num_cont
            do l = 1,ao2%num_cont
                t_ij = t_ij + ao1%coef(k)*ao2%coef(l)* &
                kinetic(ao1%exp(k),ao1%shell,ao1%position, &
                        ao2%exp(l),ao2%shell,ao2%position)
            enddo
        enddo
    end subroutine T_ij_calc


    
    !=======================================================================================================================!
    !                                               KINETIC MATRIX  CALCULATOR                                              !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                               Generates the kinetic matrix T                                          !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine T_calc(T)
    real (kind = 8) :: t_ij
    real (kind = 8) , dimension(aotot,aotot) :: T
    integer :: i,j
        do i = 1,aotot
            do j = i,aotot
                call T_ij_calc(i,j,t_ij)
                T(i,j) = t_ij
                T(j,i) = t_ij
            enddo
        enddo    
    end subroutine T_calc

    
    
    
    !=======================================================================================================================!
    !                                    MASS-VELOCITY TERM MATRIX ELEMENT CALCULATOR                                       !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates mass-velocity term energy for two contracted Gaussians                                                      !
    ! Returns float.                                                                                                        !                                                      
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine MV_ij_calc(i,j,mv_ij)
    integer , intent(in) :: i,j
    integer :: k,l
    real (kind = 8) :: mv_ij!,massvelocity
    type (ao) :: ao1, ao2
        ao1 = aos(i)
        ao2 = aos(j)
        mv_ij = 0.0
        do k = 1,ao1%num_cont
            do l = 1,ao2%num_cont
                mv_ij = mv_ij + ao1%coef(k)*ao2%coef(l)* &
                massvelocity(ao1%exp(k),ao1%shell,ao1%position, &
                        ao2%exp(l),ao2%shell,ao2%position)
            enddo
        enddo
    end subroutine MV_ij_calc


    
    !=======================================================================================================================!
    !                                           MASS-VELOCITY TERM MATRIX  CALCULATOR                                       !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the mass-velocity term matrix MV                                  !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine MV_calc(MV)
    real (kind = 8) :: mv_ij
    real (kind = 8) , dimension(aotot,aotot) :: MV
    integer :: i,j
        do i = 1,aotot
            do j = i,aotot
                call MV_ij_calc(i,j,mv_ij)
                MV(i,j) = mv_ij
                MV(j,i) = mv_ij
            enddo
        enddo    
    end subroutine MV_calc





    !=======================================================================================================================!
    !                                           NUCLEAR ATTRACTION MATRIX ELEMENT CALCULATOR                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates one electron coulombic potential between two contracted Gaussians                                           !
    ! Returns float.                                                                                                        !                                                                                                           !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================! 
    subroutine Vne_ij_calc(i,j,vne_ij)
    use nuclei_data
    integer , intent(in) :: i,j
    integer :: k,l,c
    real (kind = 8) :: vne_ij,val!,nuclear_attraction
    type (ao) :: ao1,ao2
        ao1 = aos(i)
        ao2 = aos(j)
        vne_ij = 0.0
        do k = 1,ao1%num_cont
            do l = 1,ao2%num_cont
                val = 0.0
                do c=1,nuctot
                    val = val + nuclear_attraction(ao1%exp(k),ao1%shell,ao1%position, &
                        ao2%exp(l),ao2%shell,ao2%position,nuclei(c)%position,nuclei(c)%charge)
                enddo
                vne_ij = vne_ij + ao1%coef(k)*ao2%coef(l)*val
            enddo
        enddo
    end subroutine Vne_ij_calc
    


    
    
    
    
    
    
    
    !=======================================================================================================================!
    !                                           NUCLEAR ATTRACTION MATRIX  CALCULATOR                                       !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                       Generates the nuclear attraction matrix Vne                                     !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine Vne_calc(Vne)
    real (kind = 8) :: vne_ij
    real (kind = 8) , dimension(aotot,aotot) :: Vne
    integer :: i,j
        do i = 1,aotot
            do j = i,aotot
                call Vne_ij_calc(i,j,vne_ij)
                Vne(i,j) = vne_ij
                Vne(j,i) = vne_ij
            enddo
        enddo    
    end subroutine Vne_calc



    
    
    !=======================================================================================================================!
    !                                             DARWIN TERM MATRIX ELEMENT CALCULATOR                                     !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates one electron darwin term energy between two contracted Gaussians                                            !
    ! Returns float.                                                                                                        !                                                                                                           !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================! 

    !NOTE : NO INTEGRAL SUBROUTINE NEEDED, SINCE THIS IS AN INTEGRAL WITH A DELTA FUNCTION IN INTEGRAND!!!

    subroutine D_ij_calc(i,j,d_ij)
    use nuclei_data
    integer , intent(in) :: i,j
    real (kind = 8) :: val,d_ij
    type (ao) :: ao1,ao2
    integer :: k,l,c
    ao1 = aos(i)
    ao2 = aos(j)
        d_ij = 0.0
        do k = 1,ao1%num_cont
            do l = 1,ao2%num_cont
                val = 0.0
                do c=1,nuctot
                    val = val + darwin(ao1%exp(k),ao1%shell,ao1%position, &
                        ao2%exp(l),ao2%shell,ao2%position,nuclei(c)%position,nuclei(c)%charge)
                enddo
                d_ij = d_ij + ao1%coef(k)*ao2%coef(l)*val
            enddo
        enddo
    end subroutine D_ij_calc



    

    
    
    
    !=======================================================================================================================!
    !                                               DARWIN TERM MATRIX  CALCULATOR                                          !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                             Generates the darwin term matrix D                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine D_calc(D)
    real (kind = 8) :: d_ij
    real (kind = 8) , dimension(aotot,aotot) :: D
    integer :: i,j
        D = 0.0
        do i = 1,aotot
            do j = i,aotot
                call D_ij_calc(i,j,d_ij)
            
                D(i,j) = d_ij
                D(j,i) = d_ij
            enddo
        enddo    
    end subroutine D_calc
    
    
    

    
    
    
    !=======================================================================================================================!
    !                                              J_ee INTEGRAL TENSOR GENERATOR                                           !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electron electron repulsion and exchange integrals, needed for the J and K matrices                         !
    ! Returns 4D tensors J_ee                                                                                               !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine J_ee_calc()
    use integral_tensors
    type (ao) :: ao1,ao2,ao3,ao4
    integer :: a,b,c,d,ab,cd,k,l,r,s
    real (kind = 8) :: eri,t,tt,ts,tf!,electron_repulsion
        !print(J_ee.shape)
        do b=1,aotot
            do a=b,aotot
                ab = (a*(a-1)/2)+b
                do d=1,aotot
                    do c=d,aotot
                        cd = (c*(c-1)/2)+d
                        if (ab>=cd) then
                            
                            ao1 = aos(a)
                            ao2 = aos(b)
                            ao3 = aos(c)
                            ao4 = aos(d)
                                


                            eri = 0
                            do k=1,ao1%num_cont
                                do l=1,ao2%num_cont
                                    do r=1,ao3%num_cont
                                        do s=1,ao4%num_cont
                                            eri = eri + ao1%coef(k)*ao2%coef(l)*&
                                                        ao3%coef(r)*ao4%coef(s)*&
                                                        electron_repulsion(ao1%exp(k),ao1%shell,ao1%position,&
                                                        ao2%exp(l),ao2%shell,ao2%position,&
                                                        ao3%exp(r),ao3%shell,ao3%position,&
                                                        ao4%exp(s),ao4%shell,ao4%position)
                                        enddo
                                    enddo
                                enddo
                            enddo
                            J_ee(a + aotot*( (b-1) + aotot*( (c-1) + aotot*(d-1) ) )) = eri
                            J_ee(b + aotot*( (a-1) + aotot*( (c-1) + aotot*(d-1) ) )) = eri
                            J_ee(b + aotot*( (a-1) + aotot*( (d-1) + aotot*(c-1) ) )) = eri
                            J_ee(a + aotot*( (b-1) + aotot*( (d-1) + aotot*(c-1) ) )) = eri
                            J_ee(c + aotot*( (d-1) + aotot*( (a-1) + aotot*(b-1) ) )) = eri
                            J_ee(c + aotot*( (d-1) + aotot*( (b-1) + aotot*(a-1) ) )) = eri
                            J_ee(d + aotot*( (c-1) + aotot*( (b-1) + aotot*(a-1) ) )) = eri
                            J_ee(d + aotot*( (c-1) + aotot*( (a-1) + aotot*(b-1) ) )) = eri
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine J_ee_calc

    
    

    !=======================================================================================================================!
    !                                                   J MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electronic repulsion matrix element for two contracted Gaussians                                            !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine J_ij_calc(D,r,s,j_ij)
    use integral_tensors
    integer :: j,p,q
    integer , intent(in) :: r,s
    real (kind = 8) :: j_ij
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: D
        j_ij = 0.0
        do p = 1,aotot
            do q = 1,aotot
                j_ij = j_ij + D(q,p)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
            enddo
        enddo
    end subroutine J_ij_calc


    
    

    
    !=======================================================================================================================!
    !                                                   J MATRIX  CALCULATOR                                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the electron repulsion matrix J                                   !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine J_calc(P,J)
    real (kind = 8) :: j_rs
    real (kind = 8) , dimension(aotot,aotot) :: J
    integer :: r,s
    real (kind = 8) , dimension(aotot,aotot) , intent(in) :: P
        J = 0
        do r = 1,aotot
            do s = r,aotot
                call J_ij_calc(P,r,s,j_rs)
                J(r,s) = j_rs
                J(s,r) = j_rs
            enddo
        enddo    
    end subroutine J_calc





    !=======================================================================================================================!
    !                                                   K MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electronic exchange matrix element for two contracted Gaussians                                             !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine K_ij_calc(D,r,s,k_ij)
    use integral_tensors
    integer :: p,q
    integer , intent(in) :: r,s
    real (kind = 8) :: k_ij
    real (kind = 8) , dimension(aotot,aotot) , intent(in) :: D
        k_ij = 0.0
        do p = 1,aotot
            do q = 1,aotot
                k_ij = k_ij + D(q,p)*J_ee(p + aotot*( (s-1) + aotot*( (r-1) + aotot*(q-1) ) ))
            enddo
        enddo
    end subroutine K_ij_calc


    
    

    
    !=======================================================================================================================!
    !                                                   K MATRIX  CALCULATOR                                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the electron exchange matrix K                                    !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine K_calc(D,K)
    real (kind = 8) :: k_ij
    real (kind = 8) , dimension(aotot,aotot) :: K
    integer :: r,s
    real (kind = 8) , dimension(aotot,aotot) , intent(in) :: D
        do r = 1,aotot
            do s = r,aotot
                call K_ij_calc(D,r,s,k_ij)
                K(r,s) = k_ij
                K(s,r) = k_ij
            enddo
        enddo    
    end subroutine K_calc
    
    
    
    
    
    !=======================================================================================================================!
    !                                               UHF J MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electronic repulsion matrix element for two contracted Gaussians                                            !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine Juhf_ij_calc(D,r,s,j_ij)
    use integral_tensors
    integer :: p,q
    integer , intent(in) :: r,s
    real (kind = 8) :: j_ij
    real (kind = 8) , dimension(aotot,aotot) , intent(in) :: D
        j_ij = 0.0
        do p = 1,aotot
            do q = 1,aotot
                j_ij = j_ij + D(q,p)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
            enddo
        enddo
    end subroutine Juhf_ij_calc


    
    

    
    !=======================================================================================================================!
    !                                               UHF J MATRIX  CALCULATOR                                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the electron repulsion matrix J                                   !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine Juhf_calc(D,J)
    real (kind = 8) :: j_ij
    real (kind = 8) , dimension(1:aotot,1:aotot) :: J
    integer :: r,s
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: D
        do r = 1,aotot
            do s = r,aotot
                call Juhf_ij_calc(D,r,s,j_ij)
                J(r,s) = j_ij
                J(s,r) = j_ij
            enddo
        enddo    
    end subroutine Juhf_calc





    !=======================================================================================================================!
    !                                               UHF K MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electronic exchange matrix element for two contracted Gaussians                                             !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine Kuhf_ij_calc(D,r,s,k_ij)
    use integral_tensors
    integer :: p,q
    integer , intent(in) :: r,s
    real (kind = 8) :: k_ij
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: D
        k_ij = 0.0
        do p = 1,aotot
            do q = 1,aotot
                k_ij = k_ij + D(q,p)*J_ee(p + aotot*( (s-1) + aotot*( (r-1) + aotot*(q-1) ) ))
            enddo
        enddo
    end subroutine Kuhf_ij_calc


    
    

    
    !=======================================================================================================================!
    !                                               UHF K MATRIX  CALCULATOR                                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the electron exchange matrix K                                    !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine Kuhf_calc(D,K)
    real (kind = 8) :: k_ij
    real (kind = 8) , dimension(1:aotot,1:aotot) :: K
    integer :: r,s
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: D
        do r = 1,aotot
            do s = r,aotot
                call Kuhf_ij_calc(D,r,s,k_ij)
                K(r,s) = k_ij
                K(s,r) = k_ij
            enddo
        enddo    
    end subroutine Kuhf_calc

    
    
    
    !=======================================================================================================================!
    !                                               GHF J MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electronic repulsion matrix element for two contracted Gaussians                                            !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine Jghf_ij_calc(D,r,s,j_ij)
    use integral_tensors
    integer :: p,q
    integer , intent(in) :: r,s
    real (kind = 8) :: j_ij
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: D
        j_ij = 0
        do p = 1,aotot
            do q = 1,aotot
                j_ij = j_ij + D(q,p)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
            enddo
        enddo
    end subroutine Jghf_ij_calc


    
    

    
    !=======================================================================================================================!
    !                                               GHF J MATRIX  CALCULATOR                                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the electron repulsion matrix J                                   !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine Jghf_calc(D,J)
    real (kind = 8) :: j_rs
    real (kind = 8) , dimension(1:aotot,1:aotot) :: J
    integer :: r,s,i
    real (kind = 8) , dimension(1:aotot*2,1:aotot*2) , intent(in) :: D
    real (kind = 8) , dimension(1:aotot,1:aotot) :: P
        do r = 1,aotot
            do s = r,aotot
                call Jghf_ij_calc(D,r,s,j_rs)
                J(r,s) = j_rs
                J(s,r) = j_rs
            enddo
        enddo    
    end subroutine Jghf_calc

    
    
    !=======================================================================================================================!
    !                                               GHF K MATRIX ELEMENT CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electronic exchange matrix element for two contracted Gaussians                                             !
    ! Returns float.                                                                                                        !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!  
    
    subroutine Kghf_ij_calc(D,r,s,k_ij)
    use integral_tensors
    integer :: p,q
    integer , intent(in) :: r,s
    real (kind = 8) :: k_ij
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: D
        k_ij = 0.0
        do p = 1,aotot
            do q = 1,aotot
                k_ij = k_ij + D(q,p)*J_ee(p + aotot*( (s-1) + aotot*( (r-1) + aotot*(q-1) ) ))
            enddo
        enddo
    end subroutine Kghf_ij_calc


    
    

    
    !=======================================================================================================================!
    !                                               GHF K MATRIX  CALCULATOR                                                !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                           Generates the electron exchange matrix K                                    !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine Kghf_calc(P,K)
    real (kind = 8) :: k_rs
    real (kind = 8) , dimension(1:aotot,1:aotot) :: K
    integer :: r,s,i
    real (kind = 8) , dimension(1:aotot,1:aotot) , intent(in) :: P
        do r = 1,aotot
            do s = 1,aotot
                call Kghf_ij_calc(P,r,s,k_rs)
                K(r,s) = k_rs
            enddo
        enddo    
    end subroutine Kghf_calc

end module SCF_matrix_builder
!#######################################################################################################################!

    
    
    