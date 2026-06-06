!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                    MP2 CALCULATION                                                    !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!   
module MP2_correlation
use integral_tensors
use result_data
use molecular_data
use Print_module
use Parallel
use omp_lib
use Timing

implicit none
save
contains
    subroutine MP2_calc(output,Eout)
    integer , intent(in) :: output
    real (kind = 8) , optional , intent(out) :: Eout
    integer :: O,V
    integer :: i,j,a,b,p,q,r,s
    real (kind = 8) :: iajb,ibja, ti, tf
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: uakt_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: uakb_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: iakb_tmp
    real (kind = 8) :: ei,ej,ea,eb
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ivkt_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: ivjt_tmp
    real (kind = 8) , allocatable , dimension(:,:,:,:) :: iajt_tmp

    
        if (output == 1) then
            call smallblockheader('MP2 CORRELATION CALCULATION')
            call writelines(2)
            call onelineheader('Starting MP2 routine')
            call writelines(2)
            call timer(ti)
        endif
        
        O = int(eltot/2)
        V = aotot - O
        Emp2 = 0
        if (O > V) then
            allocate(uakt_tmp(aotot,V,aotot,aotot),uakb_tmp(aotot,V,aotot,V),iakb_tmp(O,V,aotot,V))
            !write(77,*) "O>V"
            uakt_tmp = 0
            uakb_tmp = 0
            iakb_tmp = 0
            do s=1,aotot
                do r=1,aotot
                    do a=1,V
                        do q=1,aotot
                            do p=1,aotot
                                uakt_tmp(p,a,r,s) = uakt_tmp(p,a,r,s) + Cmo(q,a+O)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo


            do b=1,V
                do s=1,aotot
                    do r=1,aotot
                        do a=1,V
                            do p=1,aotot
                                uakb_tmp(p,a,r,b) = uakb_tmp(p,a,r,b) + Cmo(s,b+O)*uakt_tmp(p,a,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        

            do b=1,V
                do r=1,aotot
                    do a=1,V
                        do i=1,O
                            do p=1,aotot            
                                iakb_tmp(i,a,r,b) = iakb_tmp(i,a,r,b) + Cmo(p,i)*uakb_tmp(p,a,r,b)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
    
            do i=1,O
                ei = Eig(i)
                do a=1,V
                    ea = Eig(a+O)
                    do j=1,O
                        ej = Eig(j)
                        do b=1,V
                            eb = Eig(b+O)
                            iajb = 0
                            ibja = 0
                            do r=1,aotot
                                iajb = iajb + Cmo(r,j)*iakb_tmp(i,a,r,b)
                                ibja = ibja + Cmo(r,j)*iakb_tmp(i,b,r,a)
                            enddo
                            Emp2 = Emp2 + (2*(iajb**2) - iajb*ibja)/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(uakt_tmp,uakb_tmp,iakb_tmp)
        else
            !write(77,*) "O =< V"
            allocate(ivkt_tmp(O,aotot,aotot,aotot),ivjt_tmp(O,aotot,O,aotot),iajt_tmp(O,V,O,aotot))
            ivkt_tmp=0
            ivjt_tmp=0
            iajt_tmp=0
            do s=1,aotot
                do r=1,aotot
                    do q=1,aotot
                        do p=1,aotot
                            do i=1,O
                                ivkt_tmp(i,q,r,s) = ivkt_tmp(i,q,r,s) + Cmo(p,i)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,O
                    do r=1,aotot
                        do q=1,aotot
                            do i=1,O
                                ivjt_tmp(i,q,j,s) = ivjt_tmp(i,q,j,s) + Cmo(r,j)*ivkt_tmp(i,q,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,O
                    do a=1,V
                        do q=1,aotot
                            do i=1,O
                                iajt_tmp(i,a,j,s) = iajt_tmp(i,a,j,s) + Cmo(q,a+O)*ivjt_tmp(i,q,j,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
                            
            do b=1,V
                eb = Eig(b+O)
                do j=1,O
                    ej = Eig(j)
                    do a=1,V
                        ea = Eig(a+O)
                        do i=1,O
                            ei = Eig(i)
                            iajb = 0
                            ibja = 0
                            do s=1,aotot
                                iajb = iajb + Cmo(s,b+O)*iajt_tmp(i,a,j,s)
                                ibja = ibja + Cmo(s,b+O)*iajt_tmp(j,a,i,s)
                            enddo
                            Emp2 = Emp2 + (2*(iajb**2) - iajb*ibja)/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(ivkt_tmp,ivjt_tmp,iajt_tmp)
        endif
        
        if (present(Eout)) then
            Eout = EHF + Emp2
        endif
        
        if (output == 1) then
            call onelineheader("Correlation Energy Calculated!")
            write(77,*)''
            write(77,*) 'Hatree-Fock Energy: ',EHF
            write(77,*) 'MP2 Energy: ',EHF+Emp2
            write(77,*) 'MP2 Correlation Energy: ',Emp2
            call writelines(2)
            call twolinesfooter("MP2 CALCULATION DONE")
            call writelines(8)
            call timer(tf)
            print *, "Time taken for MP2 calculation: ",tf-ti," seconds"
        endif
    end subroutine MP2_calc
    
!#######################################################################################################################!
    
    
    
    
    
    
    
    
    
    
    
    
!#######################################################################################################################!
    
    subroutine UMP2_calc()
    integer :: Oa,Ob,Va,Vb
    integer :: i,j,a,b,p,q,r,s
    real (kind = 8) :: iajb,ibja
    real (kind = 8) :: ei,ej,ea,eb
    real (kind = 8) , allocatable :: uakt_tmp(:,:,:,:),uakb_tmp(:,:,:,:),iakb_tmp(:,:,:,:)
    real (kind = 8) , allocatable :: ivkt_tmp(:,:,:,:),ivjt_tmp(:,:,:,:),iajt_tmp(:,:,:,:)

        call smallblockheader('Unrestricted Secoind-Order M�ller-Plesset (UMP2) CORRELATION CALCULATION')
        call writelines(2)
        call onelineheader('Starting the UMP2 routine')
        call writelines(2)
        
        
        Oa = na
        Ob = nb
        Va = aotot - Oa
        Vb = aotot - Ob
        Emp2 = 0
        if (Oa > Va) then
            
            !!!! J_aaaa contribution !!!!!
            
            allocate(uakt_tmp(aotot,Va,aotot,aotot),uakb_tmp(aotot,Va,aotot,Va),iakb_tmp(Oa,Va,aotot,Va))
            uakt_tmp = 0
            uakb_tmp = 0
            iakb_tmp = 0
            do s=1,aotot
                do r=1,aotot
                    do a=1,Va
                        do q=1,aotot
                            do p=1,aotot
                                uakt_tmp(p,a,r,s) = uakt_tmp(p,a,r,s) + Ca(q,a+Oa)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo


            do b=1,Va
                do s=1,aotot
                    do r=1,aotot
                        do a=1,Va
                            do p=1,aotot
                                uakb_tmp(p,a,r,b) = uakb_tmp(p,a,r,b) + Ca(s,b+Oa)*uakt_tmp(p,a,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        

            do b=1,Va
                do r=1,aotot
                    do a=1,Va
                        do i=1,Oa
                            do p=1,aotot            
                                iakb_tmp(i,a,r,b) = iakb_tmp(i,a,r,b) + Ca(p,i)*uakb_tmp(p,a,r,b)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
    
            do i=1,Oa
                ei = Eiga(i)
                do a=1,Va
                    ea = Eiga(a+Oa)
                    do j=1,Oa
                        ej = Eiga(j)
                        do b=1,Va
                            eb = Eiga(b+Oa)
                            iajb = 0
                            ibja = 0
                            do r=1,aotot
                                iajb = iajb + Ca(r,j)*iakb_tmp(i,a,r,b)
                                ibja = ibja + Ca(r,j)*iakb_tmp(i,b,r,a)
                            enddo
                            Emp2 = Emp2 + ((iajb**2) - iajb*ibja)/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(uakt_tmp,uakb_tmp,iakb_tmp)
            

            !!!! J_bbbb contribution !!!!!
            
            allocate(uakt_tmp(aotot,Vb,aotot,aotot),uakb_tmp(aotot,Vb,aotot,Vb),iakb_tmp(Ob,Vb,aotot,Vb))
            uakt_tmp = 0
            uakb_tmp = 0
            iakb_tmp = 0
            do s=1,aotot
                do r=1,aotot
                    do a=1,Vb
                        do q=1,aotot
                            do p=1,aotot
                                uakt_tmp(p,a,r,s) = uakt_tmp(p,a,r,s) + Cb(q,a+Ob)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo


            do b=1,Vb
                do s=1,aotot
                    do r=1,aotot
                        do a=1,Vb
                            do p=1,aotot
                                uakb_tmp(p,a,r,b) = uakb_tmp(p,a,r,b) + Cb(s,b+Ob)*uakt_tmp(p,a,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        

            do b=1,Vb
                do r=1,aotot
                    do a=1,Vb
                        do i=1,Ob
                            do p=1,aotot            
                                iakb_tmp(i,a,r,b) = iakb_tmp(i,a,r,b) + Cb(p,i)*uakb_tmp(p,a,r,b)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
    
            do i=1,Ob
                ei = Eigb(i)
                do a=1,Vb
                    ea = Eigb(a+Ob)
                    do j=1,Ob
                        ej = Eigb(j)
                        do b=1,Vb
                            eb = Eigb(b+Ob)
                            iajb = 0
                            ibja = 0
                            do r=1,aotot
                                iajb = iajb + Cb(r,j)*iakb_tmp(i,a,r,b)
                                ibja = ibja + Cb(r,j)*iakb_tmp(i,b,r,a)
                            enddo
                            Emp2 = Emp2 + ((iajb**2) - iajb*ibja)/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(uakt_tmp,uakb_tmp,iakb_tmp)
            
            !!!! J_aabb (and J_bbaa indirectly) contribution !!!!!
            
            allocate(uakt_tmp(aotot,Va,aotot,aotot),uakb_tmp(aotot,Va,aotot,Vb),iakb_tmp(Oa,Va,aotot,Vb))
            uakt_tmp = 0
            uakb_tmp = 0
            iakb_tmp = 0
            do s=1,aotot
                do r=1,aotot
                    do a=1,Va
                        do q=1,aotot
                            do p=1,aotot
                                uakt_tmp(p,a,r,s) = uakt_tmp(p,a,r,s) + Ca(q,a+Oa)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo


            do b=1,Vb
                do s=1,aotot
                    do r=1,aotot
                        do a=1,Va
                            do p=1,aotot
                                uakb_tmp(p,a,r,b) = uakb_tmp(p,a,r,b) + Cb(s,b+Ob)*uakt_tmp(p,a,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        

            do b=1,Vb
                do r=1,aotot
                    do a=1,Va
                        do i=1,Oa
                            do p=1,aotot            
                                iakb_tmp(i,a,r,b) = iakb_tmp(i,a,r,b) + Ca(p,i)*uakb_tmp(p,a,r,b)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
    
            do i=1,Oa
                ei = Eiga(i)
                do a=1,Va
                    ea = Eiga(a+Oa)
                    do j=1,Ob
                        ej = Eigb(j)
                        do b=1,Vb
                            eb = Eigb(b+Ob)
                            iajb = 0
                            do r=1,aotot
                                iajb = iajb + Cb(r,j)*iakb_tmp(i,a,r,b)
                            enddo
                            Emp2 = Emp2 + (2*(iajb**2))/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(uakt_tmp,uakb_tmp,iakb_tmp)
            !END  
        else

            !!!! J_aaaa contribution !!!!!
            

            allocate(ivkt_tmp(Oa,aotot,aotot,aotot),ivjt_tmp(Oa,aotot,Oa,aotot),iajt_tmp(Oa,Va,Oa,aotot))
            ivkt_tmp=0
            ivjt_tmp=0
            iajt_tmp=0
            do s=1,aotot
                do r=1,aotot
                    do q=1,aotot
                        do p=1,aotot
                            do i=1,Oa
                                ivkt_tmp(i,q,r,s) = ivkt_tmp(i,q,r,s) + Ca(p,i)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,Oa
                    do r=1,aotot
                        do q=1,aotot
                            do i=1,Oa
                                ivjt_tmp(i,q,j,s) = ivjt_tmp(i,q,j,s) + Ca(r,j)*ivkt_tmp(i,q,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,Oa
                    do a=1,Va
                        do q=1,aotot
                            do i=1,Oa
                                iajt_tmp(i,a,j,s) = iajt_tmp(i,a,j,s) + Ca(q,a+Oa)*ivjt_tmp(i,q,j,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
                            
            do b=1,Va
                eb = Eiga(b+Oa)
                do j=1,Oa
                    ej = Eiga(j)
                    do a=1,Va
                        ea = Eiga(a+Oa)
                        do i=1,Oa
                            ei = Eiga(i)
                            iajb = 0
                            ibja = 0
                            do s=1,aotot
                                iajb = iajb + Ca(s,b+Oa)*iajt_tmp(i,a,j,s)
                                ibja = ibja + Ca(s,b+Oa)*iajt_tmp(j,a,i,s)
                            enddo
                            Emp2 = Emp2 + ((iajb**2) - iajb*ibja)/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(ivkt_tmp,ivjt_tmp,iajt_tmp)
            
            !!!! J_bbbb contribution !!!!!
            

            allocate(ivkt_tmp(Ob,aotot,aotot,aotot),ivjt_tmp(Ob,aotot,Ob,aotot),iajt_tmp(Ob,Vb,Ob,aotot))
            ivkt_tmp=0
            ivjt_tmp=0
            iajt_tmp=0
            do s=1,aotot
                do r=1,aotot
                    do q=1,aotot
                        do p=1,aotot
                            do i=1,Ob
                                ivkt_tmp(i,q,r,s) = ivkt_tmp(i,q,r,s) + Cb(p,i)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,Ob
                    do r=1,aotot
                        do q=1,aotot
                            do i=1,Ob
                                ivjt_tmp(i,q,j,s) = ivjt_tmp(i,q,j,s) + Cb(r,j)*ivkt_tmp(i,q,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,Ob
                    do a=1,Vb
                        do q=1,aotot
                            do i=1,Ob
                                iajt_tmp(i,a,j,s) = iajt_tmp(i,a,j,s) + Cb(q,a+Ob)*ivjt_tmp(i,q,j,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
                            
            do b=1,Vb
                eb = Eigb(b+Ob)
                do j=1,Ob
                    ej = Eigb(j)
                    do a=1,Vb
                        ea = Eigb(a+Ob)
                        do i=1,Ob
                            ei = Eigb(i)
                            iajb = 0
                            ibja = 0
                            do s=1,aotot
                                iajb = iajb + Cb(s,b+Ob)*iajt_tmp(i,a,j,s)
                                ibja = ibja + Cb(s,b+Ob)*iajt_tmp(j,a,i,s)
                            enddo
                            Emp2 = Emp2 + ((iajb**2) - iajb*ibja)/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(ivkt_tmp,ivjt_tmp,iajt_tmp)
            
            !!!! J_aabb (and J_bbaa indirectly) contribution !!!!!
            

            allocate(ivkt_tmp(Oa,aotot,aotot,aotot),ivjt_tmp(Oa,aotot,Ob,aotot),iajt_tmp(Oa,Va,Ob,aotot))
            ivkt_tmp=0
            ivjt_tmp=0
            iajt_tmp=0
            do s=1,aotot
                do r=1,aotot
                    do q=1,aotot
                        do p=1,aotot
                            do i=1,Oa
                                ivkt_tmp(i,q,r,s) = ivkt_tmp(i,q,r,s) + Ca(p,i)*J_ee(p + aotot*( (q-1) + aotot*( (r-1) + aotot*(s-1) ) ))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,Ob
                    do r=1,aotot
                        do q=1,aotot
                            do i=1,Oa
                                ivjt_tmp(i,q,j,s) = ivjt_tmp(i,q,j,s) + Cb(r,j)*ivkt_tmp(i,q,r,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
            do s=1,aotot
                do j=1,Ob
                    do a=1,Va
                        do q=1,aotot
                            do i=1,Oa
                                iajt_tmp(i,a,j,s) = iajt_tmp(i,a,j,s) + Ca(q,a+Ob)*ivjt_tmp(i,q,j,s)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
                            
            do b=1,Vb
                eb = Eigb(b+Ob)
                do j=1,Ob
                    ej = Eigb(j)
                    do a=1,Va
                        ea = Eiga(a+Oa)
                        do i=1,Oa
                            ei = Eiga(i)
                            iajb = 0
                            ibja = 0
                            do s=1,aotot
                                iajb = iajb + Cb(s,b+Ob)*iajt_tmp(i,a,j,s)
                            enddo
                            Emp2 = Emp2 + (2*(iajb**2))/(ei+ej-ea-eb)
                        enddo
                    enddo
                enddo
            enddo
            deallocate(ivkt_tmp,ivjt_tmp,iajt_tmp)
            !END
        endif
        call onelineheader("Correlation Energy Calculated!")
        write(77,*)''
        write(77,*) 'Unrestricted Hatree-Fock Energy: ',EHF
        write(77,*) 'UMP2 Energy: ',EHF+Emp2
        write(77,*) 'UMP2 Correlation Energy: ',Emp2
        call writelines(2)
        call twolinesfooter("UMP2 CALCULATION DONE")
        call writelines(8)
    
    end subroutine UMP2_calc
    end module MP2_correlation
!#######################################################################################################################!