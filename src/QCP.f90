
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

!console + file output printing
use Print_module

!data storage
use integral_tensors
use result_data 

!workhorses
use SCF_routines

!relativistic post hf corrections
use Relativistic

implicit none

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