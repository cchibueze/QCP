!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                        MATRIX ELEMENT INTEGRAL CALCULATORS                                            !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
        
!#######################################################################################################################! 
module integrals
use special_functions
implicit none
save
private :: E_gen,R_gen
contains
    
    
    ! The 1 and 2 electron (that is, 2 and 4 centered gaussians) integrals are evaluated by the David-McMurchie Scheme.
    ! This scheme makes use of recursive relations between hermite gaussians.
    ! For the theory, see the following references:
    ! T. Helgaker, P. Jorgensen, J. Olsen, T. Helgaker: Molecular electronic-structure theory. Wiley, New York, 2000.
    ! MAJOR CREDITS GO TO DR. J.J. GOINGS FOR THE ANALOGOUS PYTHON IMPLEMENTATION OF THESE INTEGRALS
    
    !=======================================================================================================================!
    !                                                RECURSIVE COEFFICIENT GENERATOR                                        !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Recursive definition of Hermite Gaussian coefficients.                                                                !
    ! Returns a float.                                                                                                      !
    ! a: orbital exponent on Gaussian 'a'                                                           !
    ! b: orbital exponent on Gaussian 'b'                                                            !
    ! i,j: orbital angular momentum number on Gaussian 'a' and 'b'                                                          !
    ! t: number nodes in Hermite (depends on type of integral,                                                              !
    ! e.g. always zero for overlap integrals)                                                                               !
    ! Qx: distance between origins of Gaussian 'a' and 'b'                                                                  !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    
    recursive function E_gen(i,j,t,Qx,a,b) result(En)
    integer , intent(in) :: i,j,t
    real (kind = 8) , intent(in) :: Qx,a,b
    real (kind = 8) :: En,p,q
        p = a + b
        q = (a*b)/p
        if ((t < 0) .or. (t > (i + j))) then
            ! out of bounds for t  
            En = 0.0
        elseif ((i == 0) .and. (j == 0) .and. (t == 0)) then
            ! base case
            En = exp(-q*Qx**2) ! K_AB
        elseif (j == 0) then
            ! decrement index i
            En = (1/(2*p))*E_gen(i-1,j,t-1,Qx,a,b) - &
                   ((q*Qx)/a)*E_gen(i-1,j,t,Qx,a,b)    + &
                   (t+1)*E_gen(i-1,j,t+1,Qx,a,b)
        else
            ! decrement index j
            En = (1/(2*p))*E_gen(i,j-1,t-1,Qx,a,b) + &
                   ((q*Qx)/b)*E_gen(i,j-1,t,Qx,a,b)    + &
                   (t+1)*E_gen(i,j-1,t+1,Qx,a,b)
        endif
    end function E_gen
    
    
    
    
    !=======================================================================================================================!
    !                                   COULOMB AUXILLIARY HERMITE GAUSSIAN INTEGRAL CALCULATOR                             !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Returns the Coulomb auxiliary Hermite integrals                                                                       !
    ! Returns a float.                                                                                                      !
    ! Arguments:                                                                                                            !
    ! t,u,v:   order of Coulomb Hermite derivative in x,y,z                                                                 !
    !             (see defs in Helgaker and Taylor)                                                                         !
    ! n:       order of Boys function                                                                                       !
    ! PCx,y,z: Cartesian vector distance between Gaussian                                                                   !
    !             composite center P and nuclear center C                                                                   !
    ! RPC:     Distance between P and C                                                                                     !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!
    
    recursive function R_gen(t,u,v,n,p,PCx,PCy,PCz,RPC) result(R)

    !implicit none
    integer , intent(in) :: t,u,v,n
    real (kind = 8) , intent(in) :: p,PCx,PCy,PCz,RPC
    real (kind = 8) :: Tb,R!,boys

        Tb = p*RPC*RPC
        R = 0.0
        if (t == 0 .and. u == 0 .and. v == 0) then
            R = R + ((-2*p)**n)*boys(n,Tb)
        elseif (t == 0 .and. u == 0) then
            if (v > 1) then !TAKE IN TO ACCOUNT TUV-2 !!!!
                R = R + (v-1)*R_gen(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
            endif
            R = R + PCz*R_gen(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
        elseif (t == 0) then
            if (u > 1) then
                R = R + (u-1)*R_gen(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
            endif
            R = R + PCy*R_gen(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
        else
            if (t > 1 ) then
                R = R + (t-1)*R_gen(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
            endif
            R = R + PCx*R_gen(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
        endif
        end function R_gen

    
  
    function overlap_norm_sqr(a,lmn) result(norm_sqr)
    real (kind = 8) , parameter :: pi=3.14159265359
    integer , dimension(1:3) , intent(in) :: lmn
    real (kind = 8) , intent(in) :: a
    real (kind = 8) :: f1, f2, f3, norm_sqr
    integer :: l,m,n
        l = lmn(1) ! shell angular momentum on Gaussian 'a'
        m = lmn(2)
        n = lmn(3)
        f1 = real(odd_factorial(2*l-1))
        f2 = real(odd_factorial(2*m-1))
        f3 = real(odd_factorial(2*n-1))
        norm_sqr = (f1*f2*f3) / ( (4*a)**(l+m+n) ) 
        norm_sqr = norm_sqr * (pi / (2*a))**(1.5)
    end function overlap_norm_sqr

    !=======================================================================================================================!
    !                                                   OVERLAP INTEGRAL CALCULATOR                                         !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates overlap integral between two Gaussians                                                                      !
    ! Returns a float.                                                                                                      !
    ! a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)                                                          !
    ! b: orbital exponent on Gaussian 'b' (e.g. beta in the text)                                                           !
    ! lmn1: int array containing orbital angular momentum (e.g. (1,0,0))                                                    !
    !         for Gaussian 'a'                                                                                              !
    ! lmn2: int array containing orbital angular momentum for Gaussian 'b'                                                  !
    ! nA: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]                                                       !
    ! nB: list containing origin of Gaussian 'b'                                                                             !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!
    
    function overlap(a,lmn1,nA,b,lmn2,nB) result(S)
    real (kind = 8) , parameter :: pi=3.14159265359,e=2.71828182846
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2
    integer :: l1,l2,m1,m2,n1,n2
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB
    real (kind = 8) , intent(in) :: a,b
    real (kind = 8) :: S1,S2,S3,S!,E_gen
        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
    
        S1 = E_gen(l1,l2,0,(nA(1)-nB(1)),a,b) ! X
        S2 = E_gen(m1,m2,0,nA(2)-nB(2),a,b) ! Y
        S3 = E_gen(n1,n2,0,nA(3)-nB(3),a,b) ! Z
    
        S =  S1*S2*S3*((pi/(a+b))**1.5)
    end function overlap

    
    !=======================================================================================================================!
    !                                              DIPOLE INTEGRAL CALCULATOR                                               !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates dipole integral between two Gaussians                                                                       !
    ! Returns a float.                                                                                                      !
    ! a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)                                                          !
    ! b: orbital exponent on Gaussian 'b' (e.g. beta in the text)                                                           !
    ! lmn1: int array containing orbital angular momentum (e.g. (1,0,0))                                                    !
    !         for Gaussian 'a'                                                                                              !
    ! lmn2: int array containing orbital angular momentum for Gaussian 'b'                                                  !
    ! nA: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]                                                      !
    ! nB: list containing origin of Gaussian 'b'                                                                            !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    subroutine dipole(a,lmn1,nA,b,lmn2,nB,Dx,Dy,Dz)
    real (kind = 8) , parameter :: pi=3.14159265359,e=2.71828182846
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2
    integer :: l1,l2,m1,m2,n1,n2
    real (kind = 8) :: p
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB
    real (kind = 8) , intent(in) :: a,b
    real (kind = 8) , dimension(1:3) :: nP
    real (kind = 8) :: Dx,Dy,Dz!,overlap
        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
    
        p = a + b
        nP = (a*nA + b*nB)/p ! Gaussian composite center
        
        Dx = ( E_gen(l1,l2,1,nA(1)-nB(1),a,b) + nP(1)*E_gen(l1,l2,0,nA(1)-nB(1),a,b) ) * E_gen(m1,m2,0,nA(2)-nB(2),a,b) * E_gen(n1,n2,0,nA(3)-nB(3),a,b) * ((pi/p)**1.5)
        Dy = E_gen(l1,l2,0,nA(1)-nB(1),a,b) * ( E_gen(m1,m2,1,nA(2)-nB(2),a,b) + nP(2)*E_gen(m1,m2,0,nA(2)-nB(2),a,b) ) * E_gen(n1,n2,0,nA(3)-nB(3),a,b) * ((pi/p)**1.5)
        Dz = E_gen(l1,l2,0,nA(1)-nB(1),a,b) * E_gen(m1,m2,0,nA(2)-nB(2),a,b) * ( E_gen(n1,n2,1,nA(3)-nB(3),a,b) + nP(3)*E_gen(n1,n2,0,nA(3)-nB(3),a,b) ) * ((pi/p)**1.5)
        
    end subroutine dipole
    
    

    !=======================================================================================================================!
    !                                               KINETIC INTEGRAL CALCULATOR                                             !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates kinetic integral between two Gaussians                                                                      !
    ! Returns a float.                                                                                                      !
    ! a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)                                                          !
    ! b: orbital exponent on Gaussian 'b' (e.g. beta in the text)                                                           !
    ! lmn1: int array containing orbital angular momentum (e.g. (1,0,0))                                                    !
    !         for Gaussian 'a'                                                                                              !
    ! lmn2: int array containing orbital angular momentum for Gaussian 'b'                                                  !
    ! nA: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]                                                      !
    ! nB: list containing origin of Gaussian 'b'                                                                            !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    function kinetic(a,lmn1,nA,b,lmn2,nB) result(T)
    real (kind = 8) , parameter :: pi=3.14159265359,e=2.71828182846
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2
    integer :: l1,l2,m1,m2,n1,n2
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB
    real (kind = 8) , intent(in) :: a,b
    real (kind = 8) :: term0,term1,term2,T!,overlap
        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
        
        term0 = -2*(b**2) * overlap(a,(/l1,m1,n1/),nA,b,(/l2+2,m2,n2/),nB) +          &
                b*(2*l2+1) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2,n2/),nB) +         &
                -0.5*(l2*(l2-1)) * overlap(a,(/l1,m1,n1/),nA,b,(/l2-2,m2,n2/),nB)

        term1 = -2*(b**2) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2+2,n2/),nB) +          &
                b*(2*m2+1) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2,n2/),nB) +         &
                -0.5*(m2*(m2-1)) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2-2,n2/),nB)            

        term2 = -2*(b**2) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2,n2+2/),nB) +          &
                b*(2*n2+1) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2,n2/),nB) +         &
                -0.5*(n2*(n2-1)) * overlap(a,(/l1,m1,n1/),nA,b,(/l2,m2,n2-2/),nB)

        T = term0+term1+term2
    end function kinetic
    
    



    
    !=======================================================================================================================!
    !                                         MASS-VELOCITY INTEGRAL CALCULATOR                                             !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates mass-velocity integral between two Gaussians                                                                !
    ! Returns a float.                                                                                                      !
    ! a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)                                                          !
    ! b: orbital exponent on Gaussian 'b' (e.g. beta in the text)                                                           !
    ! lmn1: int array containing orbital angular momentum (e.g. (1,0,0))                                                    !
    !         for Gaussian 'a'                                                                                              !
    ! lmn2: int array containing orbital angular momentum for Gaussian 'b'                                                  !
    ! nA: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]                                                      !
    ! nB: list containing origin of Gaussian 'b'                                                                            !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    function massvelocity(a,lmn1,nA,b,lmn2,nB) result(MV)
    real (kind = 8) , parameter :: pi=3.14159265359,e=2.71828182846
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2
    integer :: l1,l2,m1,m2,n1,n2
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB
    real (kind = 8) , intent(in) :: a,b
    real (kind = 8) :: term1,term2,term3,term4,term5,term6,MV!,overlap
        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
        MV = 0
        term1 = 16 * (a**2) * (b**2) *                             overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2+2,m2  ,n2  /),nB)  +     &
                -8 * (a**2) * b      *               (2*l2+1)    * overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                4  * (a**2) *                        (l2*(l2-1)) * overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2-2,m2  ,n2  /),nB)  +     &
                -8 * a      * (b**2) * (2*l1+1)    *               overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2+2,m2  ,n2  /),nB)  +     &
                4  * a      * b      * (2*l1+1)    * (2*l2+1)    * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                -2 * a      *          (2*l1+1)    * (l2*(l2-1)) * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2-2,m2  ,n2  /),nB)  +     &
                4  *          (b**2) * (l1*(l1-1)) *               overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2+2,m2  ,n2  /),nB)  +     &
                -2 *          b      * (l1*(l1-1)) * (2*l2+1)    * overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                                       (l1*(l1-1)) * (l2*(l2-1)) * overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2-2,m2  ,n2  /),nB)
    
        term2 = 16 * (a**2) * (b**2) *                             overlap(a,(/l1  ,m1+2,n1  /),nA,b,(/l2  ,m2+2,n2  /),nB)  +     &
                -8 * (a**2) * b      *               (2*m2+1)    * overlap(a,(/l1  ,m1+2,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                4  * (a**2) *                        (m2*(m2-1)) * overlap(a,(/l1  ,m1+2,n1  /),nA,b,(/l2  ,m2-2,n2  /),nB)  +     &
                -8 * a      * (b**2) * (2*m1+1)    *               overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2+2,n2  /),nB)  +     &
                4  * a      * b      * (2*m1+1)    * (2*m2+1)    * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                -2 * a      *          (2*m1+1)    * (m2*(m2-1)) * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2-2,n2  /),nB)  +     &
                4  *          (b**2) * (m1*(m1-1)) *               overlap(a,(/l1  ,m1-2,n1  /),nA,b,(/l2  ,m2+2,n2  /),nB)  +     &
                -2 *          b      * (m1*(m1-1)) * (2*m2+1)    * overlap(a,(/l1  ,m1-2,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                                       (m1*(m1-1)) * (m2*(m2-1)) * overlap(a,(/l1  ,m1-2,n1  /),nA,b,(/l2  ,m2-2,n2  /),nB)

        term3 = 16 * (a**2) * (b**2) *                             overlap(a,(/l1  ,m1  ,n1+2/),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                -8 * (a**2) * b      *               (2*n2+1)    * overlap(a,(/l1  ,m1  ,n1+2/),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                4  * (a**2) *                        (n2*(n2-1)) * overlap(a,(/l1  ,m1  ,n1+2/),nA,b,(/l2  ,m2  ,n2-2/),nB)  +     &
                -8 * a      * (b**2) * (2*n1+1)    *               overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                4  * a      * b      * (2*n1+1)    * (2*n2+1)    * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                -2 * a      *          (2*n1+1)    * (n2*(n2-1)) * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)  +     &
                4  *          (b**2) * (n1*(n1-1)) *               overlap(a,(/l1  ,m1  ,n1-2/),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                -2 *          b      * (n1*(n1-1)) * (2*n2+1)    * overlap(a,(/l1  ,m1  ,n1-2/),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                                       (n1*(n1-1)) * (n2*(n2-1)) * overlap(a,(/l1  ,m1  ,n1-2/),nA,b,(/l2  ,m2  ,n2-2/),nB)
    
        term4 = 16 * (a**2) * (b**2) *                             overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2+2,n2  /),nB)  +     &
                -8 * (a**2) * b      *               (2*m2+1)    * overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                4  * (a**2) *                        (m2*(m2-1)) * overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2-2,n2  /),nB)  +     &
                -8 * a      * (b**2) * (2*l1+1)    *               overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2+2,n2  /),nB)  +     &
                4  * a      * b      * (2*l1+1)    * (2*m2+1)    * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                -2 * a      *          (2*l1+1)    * (m2*(m2-1)) * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2-2,n2  /),nB)  +     &
                4  *          (b**2) * (l1*(l1-1)) *               overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2+2,n2  /),nB)  +     &
                -2 *          b      * (l1*(l1-1)) * (2*m2+1)    * overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                                       (l1*(l1-1)) * (m2*(m2-1)) * overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2-2,n2  /),nB)
    
        term5 = 16 * (a**2) * (b**2) *                             overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                -8 * (a**2) * b      *               (2*n2+1)    * overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                4  * (a**2) *                        (n2*(n2-1)) * overlap(a,(/l1+2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)  +     &
                -8 * a      * (b**2) * (2*l1+1)    *               overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                4  * a      * b      * (2*l1+1)    * (2*n2+1)    * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                -2 * a      *          (2*l1+1)    * (n2*(n2-1)) * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)  +     &
                4  *          (b**2) * (l1*(l1-1)) *               overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                -2 *          b      * (l1*(l1-1)) * (2*n2+1)    * overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                                       (l1*(l1-1)) * (n2*(n2-1)) * overlap(a,(/l1-2,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)
    
        term6 = 16 * (a**2) * (b**2) *                             overlap(a,(/l1  ,m1+2,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                -8 * (a**2) * b      *               (2*n2+1)    * overlap(a,(/l1  ,m1+2,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                4  * (a**2) *                        (n2*(n2-1)) * overlap(a,(/l1  ,m1+2,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)  +     &
                -8 * a      * (b**2) * (2*m1+1)    *               overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                4  * a      * b      * (2*m1+1)    * (2*n2+1)    * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                -2 * a      *          (2*m1+1)    * (n2*(n2-1)) * overlap(a,(/l1  ,m1  ,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)  +     &
                4  *          (b**2) * (m1*(m1-1)) *               overlap(a,(/l1  ,m1-2,n1  /),nA,b,(/l2  ,m2  ,n2+2/),nB)  +     &
                -2 *          b      * (m1*(m1-1)) * (2*n2+1)    * overlap(a,(/l1  ,m1-2,n1  /),nA,b,(/l2  ,m2  ,n2  /),nB)  +     &
                                       (m1*(m1-1)) * (n2*(n2-1)) * overlap(a,(/l1  ,m1-2,n1  /),nA,b,(/l2  ,m2  ,n2-2/),nB)
    
        MV = term1  + term2 + term3 + 2*(term4 + term5 + term6)
        MV = MV * (-1) / (8*137*137)
        !print *, l1,m1,n1
        !print *, l2,m2,m2
        !print *, term1,term2,term3
        !print *, term4,term5,term6
        !print *,''
    end function massvelocity  
    



    !=======================================================================================================================!
    !                                           NUCLEAR ATTRACTION  INTEGRAL CALCULATOR                                     !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates nuclear attraction energy integral between two Gaussians                                                    !
    ! Returns a float.                                                                                                      !
    ! a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)                                                       !
    ! b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)                                                        !
    ! lmn1: int array containing orbital angular momentum (e.g. (1,0,0))                                                    !
    !     for Gaussian 'a'                                                                                                  !
    ! lmn2: int array containing orbital angular momentum for Gaussian 'b'                                                  !
    ! nA:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]                                                   !
    ! nB:    list containing origin of Gaussian 'b'                                                                         !
    ! nC:    list containing origin of nuclear center 'C'                                                                   !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!
    
    function nuclear_attraction(a,lmn1,nA,b,lmn2,nB,nC,Zc) result(Vne)
    real (kind = 8) , parameter :: pi=3.14159265359
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2
    integer :: l1,l2,m1,m2,n1,n2,t,u,v
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB,nC
    real (kind = 8) , dimension(1:3) :: nP
    real (kind = 8) , intent(in) :: a,b,Zc
    real (kind = 8) :: nRPC,p,Vne!,R_gen,E_gen
        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
    
        p = a + b
        nP = (a*nA + b*nB)/p ! Gaussian composite center
        nRPC = ((nP(1)-nC(1))**2 + (nP(2)-nC(2))**2 + (nP(3)-nC(3))**2)**0.5

        vne = 0.0
        do t=0,(l1+l2) 
            do u=0,(m1+m2)
                do v=0,(n1+n2)
                    Vne = Vne + E_gen(l1,l2,t,nA(1)-nB(1),a,b) * &
                                E_gen(m1,m2,u,nA(2)-nB(2),a,b) * &
                                E_gen(n1,n2,v,nA(3)-nB(3),a,b) * &
                                R_gen(t,u,v,0,p,nP(1)-nC(1),nP(2)-nC(2),nP(3)-nC(3),nRPC)
                enddo
            enddo
        enddo
        Vne = Vne * -(2*pi/p)*Zc
    end function nuclear_attraction 
    
    
    !=======================================================================================================================!
    !                                            DARWIN TERM  INTEGRAL CALCULATOR                                           !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates Darwin term attraction energy integral between two Gaussians                                                !
    ! Returns a float.                                                                                                      !
    ! a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)                                                       !
    ! b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)                                                        !
    ! lmn1: int array containing orbital angular momentum (e.g. (1,0,0))                                                    !
    !     for Gaussian 'a'                                                                                                  !
    ! lmn2: int array containing orbital angular momentum for Gaussian 'b'                                                  !
    ! nA:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]                                                   !
    ! nB:    list containing origin of Gaussian 'b'                                                                         !
    ! nC:    list containing origin of nuclear center 'C'                                                                   !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!
    function darwin(a,lmn1,nA,b,lmn2,nB,nC,Zc) result(D)
    real (kind = 8) , parameter :: pi=3.14159265359,e=2.71828182846,c=137
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2
    integer :: l1,l2,m1,m2,n1,n2
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB,nC
    real (kind = 8) , intent(in) :: a,b,Zc
    real (kind = 8) , dimension(3) :: nCA,nCB
    real (kind = 8) :: nCA2,nCB2,D!,R_gen,E_gen

        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
    
        nCA = nC - nA
        nCB = nC - nB
        nCA2 = dot_product(nCA,nCA)
        nCB2 = dot_product(nCB,nCB)
        D = Zc * ( e**(-(a*nCA2+b*nCB2))) * (nCA(1)**l1)*(nCA(2)**m1)*(nCA(3)**n1)*(nCB(1)**l2)*(nCB(2)**m2)*(nCB(3)**n2)
        D = D * ( pi / (2*c*c) )
    end function darwin



    !=======================================================================================================================!
    !                                           ELECTRON REPULSION INTEGRAL CALCULATOR                                      !
    !=======================================================================================================================!
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    ! Evaluates electron electron repulsion energy integral between four Gaussians                                          !
    ! Returns a float.                                                                                                      !
    ! a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'                                                               !
    ! lmn1,lmn2                                                                                                             !
    ! lmn3,lmn4: int array containing orbital angular momentum                                                              !
    !            for Gaussian 'a','b','c','d', respectively                                                                 !
    ! nA,nB,nC,nD:   list containing origin of Gaussian 'a','b','c','d'                                                     !
    !                                                                                                                       !
    !                                                                                                                       !
    !                                                                                                                       !
    !=======================================================================================================================!

    function electron_repulsion2(a,lmn1,nA,b,lmn2,nB,c,lmn3,nC,d,lmn4,nD) result(vee)
    real (kind = 8) , parameter :: pi=3.14159265359
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2,lmn3,lmn4
    integer :: l1,l2,m1,m2,n1,n2,l3,m3,n3,l4,m4,n4,t,u,v,tau,nu,phi
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB,nC,nD
    real (kind = 8) , dimension(1:3) :: nP,nQ
    real (kind = 8) , intent(in) :: a,b,c,d
    real (kind = 8) :: nRPQ,p,q,vee,Eijt,Eklu,Emnv,Eijt2,Eklu2!,E_gen,R_gen

        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
    
        l3 = lmn3(1) ! shell angular momentum on Gaussian 'c'
        m3 = lmn3(2)
        n3 = lmn3(3)
    
        l4 = lmn4(1) ! shell angular momentum on Gaussian 'd'
        m4 = lmn4(2)
        n4 = lmn4(3)

        p = a+b ! composite exponent for nP (from Gaussians 'a' and 'b')
        q = c+d ! composite exponent for nQ (from Gaussians 'c' and 'd')
        nP = (a*nA + b*nB)/p ! nA and nB composite center
        nQ = (c*nC + d*nD)/q ! nC and nD composite center
        nRPQ = ((nP(1)-nQ(1))**2 + (nP(2)-nQ(2))**2 + (nP(3)-nQ(3))**2)**0.5

        vee = 0
        do t = 0,l1+l2
            Eijt = E_gen(l1,l2,t,nA(1)-nB(1),a,b)
            do u = 0,m1+m2
                Eklu = E_gen(m1,m2,u,nA(2)-nB(2),a,b)
                do v = 0,n1+n2
                    Emnv = E_gen(n1,n2,v,nA(3)-nB(3),a,b)
                    do tau = 0,l3+l4
                        Eijt2 = E_gen(l3,l4,tau,nC(1)-nD(1),c,d)
                        do nu = 0,m3+m4
                            Eklu2 = E_gen(m3,m4,nu,nC(2)-nD(2),c,d)
                            do phi = 0,n3+n4
                                vee = vee + Eijt * Eklu * Emnv * Eijt2 * Eklu2 * &
                                            E_gen(n3,n4,phi,nC(3)-nD(3),c,d) * &
                                            (-1)**(tau+nu+phi) * &
                                            R_gen(t+tau,u+nu,v+phi,0,p*q/(p+q),nP(1)-nQ(1),nP(2)-nQ(2),nP(3)-nQ(3),nRPQ)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        vee = vee * (2*(pi**2.5))/(p*q*((p+q)**0.5))
    end function electron_repulsion2
    
    function electron_repulsion(a,lmn1,nA,b,lmn2,nB,c,lmn3,nC,d,lmn4,nD) result(vee)
    real (kind = 8) , parameter :: pi=3.14159265359
    integer , dimension(1:3) , intent(in) :: lmn1,lmn2,lmn3,lmn4
    integer :: l1,l2,m1,m2,n1,n2,l3,m3,n3,l4,m4,n4,t,u,v,tau,nu,phi
    real (kind = 8) , dimension(1:3) , intent(in) :: nA,nB,nC,nD
    real (kind = 8) , dimension(1:3) :: nP,nQ
    real (kind = 8) , intent(in) :: a,b,c,d
    real (kind = 8) :: nRPQ,p,q,vee,Emnv,Eklu,Ftuv2,Eijt2,Eklu2!,E_gen,R_gen
    real (kind = 8) , allocatable , dimension(:,:,:) :: gcd

        l1 = lmn1(1) ! shell angular momentum on Gaussian 'a'
        m1 = lmn1(2)
        n1 = lmn1(3)
    
        l2 = lmn2(1) ! shell angular momentum on Gaussian 'b'
        m2 = lmn2(2)
        n2 = lmn2(3)
    
        l3 = lmn3(1) ! shell angular momentum on Gaussian 'c'
        m3 = lmn3(2)
        n3 = lmn3(3)
    
        l4 = lmn4(1) ! shell angular momentum on Gaussian 'd'
        m4 = lmn4(2)
        n4 = lmn4(3)
        
        !print *, 'a'
        allocate( gcd( 0:(l1+l2), 0:(m1+m2), 0:(n1+n2) ) )
        !print *, size(gcd)
        p = a+b ! composite exponent for nP (from Gaussians 'a' and 'b')
        q = c+d ! composite exponent for nQ (from Gaussians 'c' and 'd')
        nP = (a*nA + b*nB)/p ! nA and nB composite center
        nQ = (c*nC + d*nD)/q ! nC and nD composite center
        nRPQ = ((nP(1)-nQ(1))**2 + (nP(2)-nQ(2))**2 + (nP(3)-nQ(3))**2)**0.5
        
        gcd = 0
        do tau = 0,l3+l4
            Eijt2 = E_gen(l3,l4,tau,nC(1)-nD(1),c,d)
            do nu = 0,m3+m4
                Eklu2 = E_gen(m3,m4,nu,nC(2)-nD(2),c,d)
                do phi = 0,n3+n4
                    Ftuv2 = Eijt2 * Eklu2 * E_gen(n3,n4,phi,nC(3)-nD(3),c,d) * (-1)**(tau+nu+phi)
                    
                    do v = 0,n1+n2
                        do u = 0,m1+m2
                            do t = 0,l1+l2
                                gcd(t,u,v) = gcd(t,u,v) + Ftuv2 * &
                                            R_gen(t+tau,u+nu,v+phi,0,p*q/(p+q),nP(1)-nQ(1),nP(2)-nQ(2),nP(3)-nQ(3),nRPQ)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        
        vee = 0
        do v = 0,n1+n2
            Emnv = E_gen(n1,n2,v,nA(3)-nB(3),a,b)
            do u = 0,m1+m2
                Eklu = E_gen(m1,m2,u,nA(2)-nB(2),a,b)
                do t = 0,l1+l2
                    vee = vee + E_gen(l1,l2,t,nA(1)-nB(1),a,b) * gcd(t,u,v) * Eklu * Emnv 
                enddo
            enddo
        enddo
        vee = vee * (2*(pi**2.5))/(p*q*((p+q)**0.5))
        deallocate(gcd)
    end function electron_repulsion
end module integrals   
!#######################################################################################################################!    