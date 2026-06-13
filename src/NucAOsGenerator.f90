!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                             NUCLEI AND AOS GENERATOR                                                  !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
!#######################################################################################################################! 
module nuclei_aos_generator
use molecular_data
implicit none
save
private
public :: get_nuclei, get_aos, normalize_aos, print_mol_data

!#######################################################################################################################!  
    contains

    subroutine get_nuclei()
    use atomic_data
    use nuclei_data
    use calculation_data
    integer :: i
    real (kind = 8) :: x,y,z
    character (len=30) :: label
        call set_nuclearcharge()
        allocate(nuclei(nuctot))
        do i=1,nuctot
            x = xcoords(i)
            y = ycoords(i)
            z = zcoords(i)
            label = labels(i)
            call create_nucleus(i,label,x,y,z)
        enddo
    end subroutine get_nuclei
    
    subroutine create_nucleus(nuc_no,label,x,y,z) !CREATES N'TH NUCLEUS FOR A CERTAIN ATOMIC LABEL AND POSITION
    use nuclei_data
    character (len=30)  :: label
    character (len=30) :: name,count
    real (kind = 8) :: x,y,z,mass,charge
    integer :: nuc_no,ao_on_nuc,nuc_type_no
        call get_nucleus_data(label,mass,charge,nuc_type_no,count)
        name = trim(label)//trim(count)
        nuclei(nuc_no)%name = name
        nuclei(nuc_no)%position = (/x,y,z/)
        nuclei(nuc_no)%mass = mass
        nuclei(nuc_no)%charge = charge
        nuclei(nuc_no)%type_no = nuc_type_no 
    end subroutine create_nucleus

    subroutine get_nucleus_data(label,mass,charge,nuc_type_no,count) !GIVES NUCLEUS PROPERTIES FOR A CERTTAIN ATOMIC LABEL 
    use atomic_data
    character (len=12) , intent(inout) :: label
    character (len=30) , intent(out) :: count
    real (kind = 8) , intent(out) :: mass,charge
    integer , intent(out) :: nuc_type_no
    integer :: j
        do j=1,100 
            if (trim(label)==trim(atomsymbol(j))) then
                mass = atomic_mass(j) * 1836
                charge = nuclear_charge(j)
                nuc_type_no = j
                atomcount(j) = atomcount(j) + 1
                write(count,*) atomcount(j) 
                count = adjustl(count)
            endif
        enddo
    end subroutine get_nucleus_data
    
    !#######################################################################################################################!
    
    subroutine get_aos()
    use calculation_data
    use nuclei_data
    use atomic_data
    use ao_data
    use Print_module
    real (kind = 8) , parameter :: pi=3.14159265359
    integer :: i,j,k,l,ios,nuc_type_no,num_cont,junk
    character (len=30) :: scnt,pcnt,dcnt,fcnt
    character(len=30) :: file_name,nuc_name,atom_name,read_line,angmom
    integer , dimension(4) :: shell_count = (/0,0,0,0/)
    character(len=10) , dimension(1)  :: ssuf=(/'s'/)
    character(len=10) , dimension(3)  :: psuf=(/'px','py','pz'/)
    character(len=10) , dimension(6)  :: dsuf=(/'dxx','dyy','dzz','dxy','dxz','dyz'/)
    character(len=10) , dimension(10) :: fsuf=(/'fxxx','fyyy','fzzz','fxxy','fxxz','fxyy','fyyz','fxzz','fyzz','fxyz'/)
    integer , dimension(3)   :: s_shells =       (/0,0,0/)
    integer , dimension(3,3) :: p_shells=reshape((/1,0,0,  &
                                                   0,1,0,  &
                                                   0,0,1/),&
                                                   shape(p_shells), order=(/2,1/) )
    integer , dimension(6,3) :: d_shells=reshape((/2,0,0,  &
                                                   0,2,0,  &
                                                   0,0,2,  &
                                                   1,1,0,  &
                                                   1,0,1,  &
                                                   0,1,1/),&
                                                   shape(d_shells), order=(/2,1/) )
    integer , dimension(10,3) :: f_shells=reshape((/3,0,0,  &
                                                   0,3,0,  &                                                   
                                                   0,0,3,  &
                                                   2,1,0,  &
                                                   2,0,1,  &
                                                   1,2,0,  &
                                                   0,2,1,  &
                                                   1,0,2,  &
                                                   0,1,2,  &
                                                   1,1,1/),&
                                                   shape(f_shells), order=(/2,1/) )
    real (kind = 8) , dimension(1:3) :: nuc_pos
    real (kind = 8) , allocatable , dimension(:) :: exp_tmp, coef_tmp, coef_tmp2
    type (ao) , allocatable , dimension(:) :: aos_tmp
        aotot = 0
        allocate(aos(aotot))
        file_name = './../../bases/'//trim(basis)//'.bas'

        open(2, file = file_name, iostat=ios,status='old')
        if (ios/=0) then
            write(*,*)'Error opening file',file_name, 'please check file' 
            stop
        endif
        
        do i=1,nuctot
            ! ensuring read line is not an atom name
            read_line = 'not an atomname'
            
            ! rewind file for reading in AOs for the next nucleus
            rewind(2)
            
            nuc_type_no = nuclei(i)%type_no
            nuc_name = nuclei(i)%name
            nuc_pos = nuclei(i)%position
            atom_name = trim(atomname(nuc_type_no))
            shell_count = 0
                        
            do while (read_line /= atom_name) 
                read(2,*,iostat=ios) read_line
                if (ios/=0) then
                    write(*,*)'Error reading basis set file, please check file' 
                    stop
                endif
            enddo

            read(2,*) angmom

            do while (angmom == 'S' .or. angmom == 'P' .or. angmom == 'L' .or. angmom == 'D' .or. angmom == 'F')
                backspace(2)
                read(2,*) angmom, num_cont
                
                allocate(exp_tmp(num_cont),coef_tmp(num_cont))
                if (angmom == 'L') then
                    allocate(coef_tmp2(num_cont))
                    do j = 1,num_cont
                        read(2,*) junk, exp_tmp(j), coef_tmp(j), coef_tmp2(j)
                    enddo
                else
                    do j = 1,num_cont
                        read(2,*) junk, exp_tmp(j), coef_tmp(j)
                    enddo
                endif
                
                select case (angmom)
                    case ('S') 
                        ! ensure each gaussian is normalized w.r.t. read in contraction coeffcients
                        do l=1,num_cont
                            coef_tmp(l) = coef_tmp(l) * 1!((2.0*(exp_tmp(l)/pi))**0.75)
                        enddo
                        
                        shell_count(1) = shell_count(1) + 1
                        write(scnt,*) shell_count(1)
                        scnt = adjustl(scnt)
                        !print *, scnt
                        !aos ==> aos_tmp
                        allocate(aos_tmp(aotot+1))
                        aos_tmp(1:aotot) = aos(1:aotot) 
                        deallocate(aos)
                        allocate(aos(aotot+1))     
                        
                        aos_tmp(aotot+1)%name = trim(nuc_name)//'_'//trim(ssuf(1))//'_'//trim(scnt) 
                        aos_tmp(aotot+1)%position = nuc_pos
                        aos_tmp(aotot+1)%shell = s_shells(1:3)
                        allocate(aos_tmp(aotot+1)%exp(num_cont),aos_tmp(aotot+1)%coef(num_cont)) 
                        aos_tmp(aotot+1)%exp = exp_tmp 
                        aos_tmp(aotot+1)%coef = coef_tmp
                        aos_tmp(aotot+1)%nucno = i
                        aos_tmp(aotot+1)%nuc_type_no = nuc_type_no
                        aos_tmp(aotot+1)%num_cont = num_cont
                        
                        
                        !aos_tmp ==> aos
                        aos = aos_tmp 
                        deallocate(aos_tmp)
                        
                        !update ao count
                        aotot = aotot + 1
                    case ('P') 
                        ! ensure each gaussian is normalized w.r.t. read in contraction coeffcients
                        do l=1,num_cont
                            coef_tmp(l) = coef_tmp(l) * 1!((128*((exp_tmp(l)**5)/(pi**3)))**0.25)
                        enddo
                        shell_count(2) = shell_count(2) + 1
                        write(pcnt,*) shell_count(2)
                        pcnt = adjustl(pcnt)
                        
                        !aos ==> aos_tmp
                        allocate(aos_tmp(aotot+3))
                        aos_tmp(1:aotot) = aos(1:aotot) 
                        deallocate(aos)
                        allocate(aos(aotot+3))                       
                        
                        do k=1,3
                            aos_tmp(aotot+k)%name = trim(nuc_name)//'_'//trim(psuf(k))//'_'//trim(pcnt) 
                            aos_tmp(aotot+k)%position = nuc_pos
                            aos_tmp(aotot+k)%shell = p_shells(k,1:3)
                            allocate(aos_tmp(aotot+k)%exp(num_cont),aos_tmp(aotot+k)%coef(num_cont)) 
                            aos_tmp(aotot+k)%exp = exp_tmp 
                            aos_tmp(aotot+k)%coef = coef_tmp
                            aos_tmp(aotot+k)%nucno = i
                            aos_tmp(aotot+k)%nuc_type_no = nuc_type_no
                            aos_tmp(aotot+k)%num_cont = num_cont
                        enddo

                        !aos_tmp ==> aos
                        aos = aos_tmp 
                        deallocate(aos_tmp)
                        
                        !update ao count
                        aotot = aotot + 3

                    case ('L')
                        ! ensure each gaussian is normalized w.r.t. read in contraction coeffcients
                        do l=1,num_cont
                            coef_tmp(l)  = coef_tmp(l)  * 1!((2.0*( exp_tmp(l)    / pi    ))**0.75)
                            coef_tmp2(l) = coef_tmp2(l) * 1!((128*((exp_tmp(l)**5)/(pi**3)))**0.25)
                        enddo
                        
                        shell_count(1) = shell_count(1) + 1
                        shell_count(2) = shell_count(2) + 1
                        write(scnt,*) shell_count(1)                        
                        write(pcnt,*) shell_count(2)
                        scnt = adjustl(scnt)
                        pcnt = adjustl(pcnt)

                        !aos ==> aos_tmp
                        allocate(aos_tmp(aotot+4))
                        aos_tmp(1:aotot) = aos(1:aotot) 
                        deallocate(aos)
                        allocate(aos(aotot+4))

                        aos_tmp(aotot+1)%name = trim(nuc_name)//'_'//trim(ssuf(1))//'_'//trim(scnt) 
                        aos_tmp(aotot+1)%position = nuc_pos
                        aos_tmp(aotot+1)%shell = s_shells(1:3)
                        allocate(aos_tmp(aotot+1)%exp(num_cont),aos_tmp(aotot+1)%coef(num_cont)) 
                        aos_tmp(aotot+1)%exp = exp_tmp 
                        aos_tmp(aotot+1)%coef = coef_tmp
                        aos_tmp(aotot+1)%nucno = i
                        aos_tmp(aotot+1)%nuc_type_no = nuc_type_no
                        aos_tmp(aotot+1)%num_cont = num_cont                         
                        
                        do k=1,3
                            aos_tmp(aotot+k+1)%name = trim(nuc_name)//'_'//trim(psuf(k))//'_'//trim(pcnt) 
                            aos_tmp(aotot+k+1)%position = nuc_pos
                            aos_tmp(aotot+k+1)%shell = p_shells(k,1:3)
                            allocate(aos_tmp(aotot+k+1)%exp(num_cont),aos_tmp(aotot+k+1)%coef(num_cont)) 
                            aos_tmp(aotot+k+1)%exp = exp_tmp 
                            aos_tmp(aotot+k+1)%coef = coef_tmp2
                            aos_tmp(aotot+k+1)%nucno = i
                            aos_tmp(aotot+k+1)%nuc_type_no = nuc_type_no
                            aos_tmp(aotot+k+1)%num_cont = num_cont
                        enddo

                        !aos_tmp ==> aos
                        aos = aos_tmp 
                        deallocate(aos_tmp)
                        
                        !update ao count
                        aotot = aotot + 4
            
                    case ('D')
                        ! ensure each gaussian is normalized w.r.t. read in contraction coeffcients
                        do l=1,num_cont
                            coef_tmp(l) = coef_tmp(l) * 1!((2.0*(exp_tmp(l)/pi))**0.75)
                        enddo
                        
                        shell_count(3) = shell_count(3) + 1
                        write(dcnt,*) shell_count(3)
                        dcnt = adjustl(dcnt)

                        !aos ==> aos_tmp
                        allocate(aos_tmp(aotot+6))
                        aos_tmp(1:aotot) = aos(1:aotot) 
                        deallocate(aos)
                        allocate(aos(aotot+6))

                        do k=1,6
                            aos_tmp(aotot+k)%name = trim(nuc_name)//'_'//trim(dsuf(k))//'_'//trim(dcnt) 
                            aos_tmp(aotot+k)%position = nuc_pos
                            aos_tmp(aotot+k)%shell = d_shells(k,1:3)
                            allocate(aos_tmp(aotot+k)%exp(num_cont),aos_tmp(aotot+k)%coef(num_cont)) 
                            aos_tmp(aotot+k)%exp = exp_tmp 
                            aos_tmp(aotot+k)%nucno = i
                            aos_tmp(aotot+k)%nuc_type_no = nuc_type_no
                            aos_tmp(aotot+k)%num_cont = num_cont
                        enddo      

                        !aos_tmp ==> aos
                        aos = aos_tmp 
                        deallocate(aos_tmp)
                        
                        !update ao count
                        aotot = aotot + 6
                        
                    case ('F')
                        ! ensure each gaussian is normalized w.r.t. read in contraction coeffcients
                        do l=1,num_cont
                            coef_tmp(l) = coef_tmp(l) * 1!((2.0*(exp_tmp(l)/pi))**0.75)
                        enddo 
                        
                        shell_count(4) = shell_count(4) + 1
                        write(fcnt,*) shell_count(4)
                        fcnt = adjustl(fcnt)
                        
                        !aos ==> aos_tmp
                        allocate(aos_tmp(aotot+10))
                        aos_tmp(1:aotot) = aos(1:aotot) 
                        deallocate(aos)
                        allocate(aos(aotot+10))     
            
                        do k=1,10
                            aos_tmp(aotot+k)%name = trim(nuc_name)//'_'//trim(fsuf(k))//'_'//trim(fcnt) 
                            aos_tmp(aotot+k)%position = nuc_pos
                            aos_tmp(aotot+k)%shell = f_shells(k,1:3)
                            allocate(aos_tmp(aotot+k)%exp(num_cont),aos_tmp(aotot+k)%coef(num_cont)) 
                            aos_tmp(aotot+k)%exp = exp_tmp 
                            aos_tmp(aotot+k)%coef = coef_tmp
                            aos_tmp(aotot+k)%nucno = i
                            aos_tmp(aotot+k)%nuc_type_no = nuc_type_no
                            aos_tmp(aotot+k)%num_cont = num_cont
                        enddo

                        !aos_tmp ==> aos
                        aos = aos_tmp 
                        deallocate(aos_tmp)
                        
                        !update ao count
                        aotot = aotot + 10
                            
                end select
                deallocate(exp_tmp,coef_tmp)
                if (angmom == 'L') deallocate(coef_tmp2)

                ! checking for further read-in of AOs on this nucleus
                read(2,*) angmom
                
            enddo
            
            
        enddo
        
        !close basis set file!
        close(2)
               

        nao = aotot
        nmo = nao
        nocc = int(eltot/2.0d0)
        nvir = nao - nocc

        call printtext('number of AOs: ',int=nao)
        call printtext('number of MOs: ',int=nmo)
        call printtext('number of occupied orbitals: ',int=nocc)
        call printtext('number of virtual orbitals: ',int=nvir)

    end subroutine get_aos
    
    !#######################################################################################################################! 
    
    subroutine normalize_aos()
    use ao_data
    use integrals
    integer :: i,k
    real (kind = 8) :: s_gg
        do i=1,aotot
            do k = 1,aos(i)%num_cont
                s_gg = overlap(aos(i)%exp(k),aos(i)%shell,aos(i)%position, &
                            aos(i)%exp(k),aos(i)%shell,aos(i)%position)
                aos(i)%coef(k) = aos(i)%coef(k) / sqrt(s_gg)
            enddo
        enddo

    end subroutine

    !#######################################################################################################################! 
    
    subroutine print_mol_data()
    use print_module
    use nuclei_data
    use ao_data
    integer :: a,b,no=0
        open(unit=77,file='../../OUTPUT.dat')
        !print *, eltot
        if ( (real(eltot)/2 - ceiling(real(eltot)/2) /= 0) .and. (real(multiplicity)/2 - ceiling(real(multiplicity)/2) /= 0) ) then
            write(77,*) 'Wrongly chosen multiplicity for the number of electrons.'
            write(77,*) 'Odd multiplicity for even number of electrons and vice versa.' 
            stop
        else if ( (real(eltot)/2 - ceiling(real(eltot)/2) == 0) .and. (real(multiplicity)/2 - ceiling(real(multiplicity)/2) == 0) ) then
            write(77,*) 'Wrongly chosen multiplicity for the number of electrons,'
            write(77,*) 'Odd multiplicity for even number of electrons and vice versa.'
            stop
        endif

        
        call printtext('Calculation Running; Check Output file when finished!')
        
        write(77,*) ''
        call bigblockheader('QCP PROGRAM OUTPUT')
        call writelines(10)
        
        call smallblockheader('MOLECULAR TOPOLOGY AND ELECTRONIC DATA')
        call writelines(2)
        
        call onelineheader('AMOUNT OF MOLECULAR SUB-ENTITIES')
        write(77,*) ''
        write(77,*) 'number of nuclei:   ',nuctot
        write(77,*) 'number of AOs:      ',aotot
        write(77,*) 'number of electrons:',eltot
        write(77,*) 'molecular charge:   ',molcharge
        write(77,*) 'spin multiplicity:  ',multiplicity
        call writelines(3)
        
        call onelineheader('NUCLEAR DATA (Atomic Units)')
        write(77,*) ''
        write(77,2) 'name', 'mass', 'charge', 'x-pos', 'y-pos', 'z-pos', 'nuc_type_no'
    2   format(A5,A8,A9,A8,A8,A9,A13) 
        write(77,*)''
        do a = 1,nuctot !PRINT INITIAL TOPOLOGY DETAILS   
            write(77,10) nuclei(a)%name , nuclei(a)%mass, nuclei(a)%charge, (nuclei(a)%position(b) , b=1,3), nuclei(a)%type_no
            10 format(1x,A5,5f8.2,I10,I10)
            write(77,*)''
        enddo
    
        call writelines(2)
        call onelineheader('BASIS SET DATA (Atomic Units)')
        write(77,*)''
        do a = 1,aotot
            write(77,*) (aos(a)%name)
            write(77,* ) , 'nucleus type number: ', aos(a)%nucno
            write(77,'(1x,A13,3f12.5)') , 'position:     ', (aos(a)%position)
            write(77,100) 'shell:        ',(aos(a)%shell)
            100 format(1x,A14,3I12)    
            write(77,'(1x,A14,3I5)') , 'nuclei number:', (aos(a)%nucno)
            write(77,'(1x,A13,10f12.5)') , 'exponents:   ', (aos(a)%exp)
            write(77,'(1x,A13,10f12.5)') , 'coefficients: ', (aos(a)%coef)
            write(77, *) , 'number of contractions: ', size(aos(a)%coef)
            
            write(77,*)''
        enddo
        
        call writelines(2)
        call twolinesfooter('END MOLECULAR TOPOLOGY AND ELECTRONIC DATA')
        call writelines(7)
    end subroutine print_mol_data

end module nuclei_aos_generator