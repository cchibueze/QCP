!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                      INPUT PROCESSOR                                                  !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
!#######################################################################################################################! 
module input_processor
implicit none
save

!#######################################################################################################################!
    contains
 
    subroutine inp_reader(filename)
    use calculation_data
    use molecular_data
    integer :: i,ios
    real (kind = 8) :: x,y,z
    real (kind = 8) :: a2b ! 1/a0 to go from Ångstrom to Bohr units of distance
    character (len=30) , intent(in) :: filename
    character (len=30) :: junk, length_units, atom_symbol

        open(1, file = filename, iostat=ios,status='old')
        if (ios/=0) then
            write(*,*)'Error opening file, please check file' 
            stop
        endif
            print *, 'hey'
        read(1,*) junk, basis
        basis = trim(adjustl(basis))
        read(1,*)
        
        read(1,*) junk, ngc_go
        if (ngc_go == 1) then
            read(1,*) junk, go_wf
            read(1,*) junk, ngc_max_iter
            read(1,*) junk, mg_ngc
        else
            read(1,*)
            read(1,*)
        endif
    
        read(1,*)
    
        read(1,*) junk, nr_go
        
        if (nr_go == 1) then
            read(1,*) junk, NR_max_iter
            read(1,*) junk, mg_NR
        else
            read(1,*)
            read(1,*)
        endif
 
        read(1,*)
    
        read(1,*) junk, sp
        read(1,*)
        read(1,*) junk, mp2
        read(1,*)
        read(1,*) junk, cid
        read(1,*)
        
        read(1,*) junk, cc
        
        if (cc == 1) then
            read(1,*) junk, cctype
        else
            read(1,*)
        endif
        
        read(1,*)
        
        read(1,*) junk, src
        read(1,*)
        read(1,*) junk, fc
        read(1,*)
        read(1,*) junk, cises
        read(1,*)
        read(1,*) junk, tdhf
        read(1,*)
        
        read(1,*) junk, molcharge
        read(1,*)
        read(1,*) junk, multiplicity
        read(1,*)
        read(1,*) junk, length_units
        read(1,*)

        nuctot = 0 
        do
            read(1,*,iostat=ios)
            if (ios/=0) exit
            nuctot = nuctot + 1
        enddo
    
        allocate(xcoords(nuctot))
        allocate(ycoords(nuctot))
        allocate(zcoords(nuctot))
        allocate(labels(nuctot))

        rewind(1)
        do i=1,34
            read(1,*)
        enddo
    
        if (length_units == 'A') then
            a2b = 1.88972612463
        else
            a2b = 1
        endif
    
        do i=1,nuctot
            read(1,*) atom_symbol,x,y,z
            labels(i)=trim(atom_symbol)
            xcoords(i)=x*a2b
            ycoords(i)=y*a2b
            zcoords(i)=z*a2b
        enddo

        close(1)
    end subroutine inp_reader
!#######################################################################################################################! 
    subroutine set_eltot()
    use nuclei_data
    use molecular_data
    integer :: i
        eltot = 0
        do i = 1,nuctot
            eltot = eltot + int(nuclei(i)%charge)
        enddo
        eltot = eltot - molcharge
    end subroutine
    end module input_processor  
!#######################################################################################################################! 