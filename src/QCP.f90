
!#######################################################################################################################!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                #####              ########           #####                                          #!
!#                               ##   ##            ####               ##   ##                                         #!
!#                              ##     ##          ##                 ##   ##                                          #!
!#                             ##       ##        ##                 ##   ##                                           #!
!#                            ##       ##        ##                 ######                                             #!
!#                            ##     ###         ##                ##                                                  #!
!#                             ##   #####         ###             ##                                                   #!
!#                              #####   ####       ########      ##                                                    #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#                                                                                                                     #!
!#              NON-)RELATIVISTIC (POST)-HATREE-FOCK/ROOTHAAN-HALL ELECTRONIC STRUCTURE CALCULATIONS                   #!
!#                                                                                                                     #!
!#                                      A QUANTUM CHEMICAL AB INITIO PROGRAM                                           #!
!#                                                                                                                     #!
!#                                          WRITTEN BY CHIMA CHIBUEZE                                                  #!    
!#                                                                                                                     #!    
!#                                                                                                                     #!    
!#######################################################################################################################!
    
    
    
    
    
    
!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                              THE QCP MAIN PROGRAM                                                     !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
program QCP

!initialization + system setup
use input_processor
use calculation_data
use nuclei_aos_generator

!gradient based geometry optimization
use nuclear_derivatives
use geometry_optimization

!console + file output printing
use Print_module

!data storage
use integral_tensors
use result_data 

!workhorses
use SCF_routines

!post HF WFT methods
use MP2_correlation
use CI_correlation
use CC_correlation

!vibrational frequency analysis
use nuclear_vibration

!relativistic post hf corrections
use Relativistic

implicit none

real (kind = 8) :: Ef
real (kind = 8) :: tf,ts
character (len=30) :: filename
integer :: ghf=0

    call CPU_TIME(ts)

    !=======================================================================================================================!
    !                                               INITIALIZE TOPOLOGY                                                     !
    !=======================================================================================================================!
    filename = 'INPUT.dat'
    call inp_reader(filename)
    call get_nuclei()
    call set_eltot()
    call get_aos()
    call normalize_aos()
    call print_mol_data()
    open(unit=77,file='OUTPUT.dat')
    
    !=======================================================================================================================!
    !                                              STARTING A QCPCC CALCULATION                                             !
    !=======================================================================================================================!
     
    allocate(H(aotot,aotot))
    allocate(J_ee(aotot**4))
    if (multiplicity == 1) then
        allocate(Eig(aotot))
        allocate(Cmo(aotot,aotot))
    endif
    max_iter=200
    eps=1E-5
    yNR= 1.0

    
    !=======================================================================================================================!
    !                           PERFORMING A NON-LINEAR CONJUGENT GRADIENT GEOMETRY OPTIMIZATION                            !
    !=======================================================================================================================!
    if (ngc_go == 1) then
        print *, 'Doing Non-Linear Conjugent Gradient Geometry Optimization!!'
        call NGC(dh,Ef)
    endif

    !=======================================================================================================================!
    !                                   PERFORMING A NEWTON-RAPHSON GEOMETRY OPTIMIZATION                                   !
    !=======================================================================================================================!  
    if (nr_go == 1) then
        print *, 'Newton-Raphson Geometry Optimization!!'
        call Newton_Raphson(dhH,Ef)
    endif


    !=======================================================================================================================!
    !                                       PERFORMING AN SINGLE POINT HF CALCULATION                                       !
    !=======================================================================================================================!
    if (sp == 1) then
        if (multiplicity == 1) then
            print *, 'Doing HF!'
            call HFC(1)
        else if (ghf == 1) then
            print *, 'Doing GHF!'
            allocate(Eig(aotot*2))
            allocate(Cmo(aotot*2,aotot*2))
            call GHFC(1)
        else
            print *, 'Doing UHF!'
            allocate(Eiga(aotot),Eigb(aotot))
            allocate(Ca(aotot,aotot),Cb(aotot,aotot))
            call UHFC(1)
        endif
    endif

    !=======================================================================================================================!
    !                                      PERFORMING AN MP2 CORRELATION CALCULATION                                        !
    !=======================================================================================================================!
    if (mp2 == 1) then
        if (multiplicity == 1) then
            print *, 'Doing MP2!!'
            call MP2_calc(1)
        else
            print *, 'Doing UMP2!!'
            call UMP2_calc()
        endif
    endif
    
    !=======================================================================================================================!
    !                                      PERFORMING A CC CORRELATION CALCULATION                                         !
    !=======================================================================================================================!
    if (cc == 1) then
        
        if (cctype == 'D') then
            print *, 'Doing CCD!!'
            call closed_shell_CCD_calc(1)
        else if (cctype == 'SD') then
            print *, 'Doing CCSD!!'
            call closed_shell_CCSD_calc(1)
        endif
        
    endif 
    
    !=======================================================================================================================!
    !                                      PERFORMING A CI CORRELATION CALCULATION                                          !
    !=======================================================================================================================!
    if (cid == 1) then
        print *, 'Doing Closed Shell CID!!'
        call closed_shell_CID_calc(1)
    endif
    
    !=======================================================================================================================!
    !                                             PERFORMING A CISES CALCULATION                                            !
    !=======================================================================================================================!
    if (cises == 1) then
        print *, 'Calculating CIS excitation spectra!!'
        !call CIS_multiple_davidson_calc(1)
        call closed_shell_CISES_calc(1) !
    endif

    !=======================================================================================================================!
    !                                             PERFORMING A TD-HF CALCULATION                                            !
    !=======================================================================================================================!
    if (tdhf == 1) then
        print *, 'Calculating TD-HF excitation spectra!!'
        call closed_shell_TDHF_calc(1)
    endif

    !=======================================================================================================================!
    !                                PERFORMING A SCALAR RELATIVISTIC CORRECTION CALCULATION                                !
    !=======================================================================================================================!
    if (src == 1) then 
        print *, 'Calculating Scalar Relativistic Correction!!'
        call Relativistic_correction()
    endif

    
    call writelines(3)
    write(77,*)'FINAL END OF CALCULATION'
    call cpu_time(tf)
    write(77,32) "total run time:", &
    tf - ts, "seconds"
    32 format(1x,A15,f11.2,A8)
    write(77,*)'======================================================================='
    write(77,*)''
    close(77)
end program QCP
!#######################################################################################################################! 