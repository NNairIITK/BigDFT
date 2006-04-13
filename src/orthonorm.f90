

	subroutine loewe(iproc,norb,nvctr,nseg,keyg,keyv,psi)
! loewdin orthogonalisation
        implicit real*8 (a-h,o-z)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        real*8, allocatable :: ovrlp(:,:,:),evall(:),tpsi(:)

 if (norb.eq.1) then
      iorb=1
      mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
      mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
      ipsi_c=nvctr(2*iorb-2)+1
      ipsi_f=nvctr(2*iorb-1)+1
      call  wnrm(mvctr_c,mvctr_f,psi(ipsi_c),psi(ipsi_f),scpr) ; scpr=1.d0/sqrt(scpr)
      call wscal(mvctr_c,mvctr_f,scpr,psi(ipsi_c),psi(ipsi_f)) 
 else

        allocate(ovrlp(norb,norb,3),evall(norb),tpsi(nvctr(2*norb)))

        offdiag=0.d0
     do 100,iorb=1,norb
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1
       call wcopy(mavctr_c,mavctr_f,psi(ipsi_c),psi(ipsi_f),tpsi(ipsi_c),tpsi(ipsi_f))

! Full matrix for testing
!         do 100,jorb=1,norb
! Lower triangle
         do 100,jorb=1,iorb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

	call wdot(  & 
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),  & 
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),psi(jpsi_c),psi(jpsi_f),scpr)

        if (iorb.ne.jorb) offdiag=offdiag+abs(scpr)
       ovrlp(iorb,jorb,1)=scpr
100    continue
       if (iproc.eq.0) write(*,*) 'offdiagonal overlap',offdiag
      
! testing
!       tt=0.d0
!       do jorb=1,norb
!       do iorb=1,jorb-1
!       tt=max(tt,abs(ovrlp(jorb,iorb,1)-ovrlp(iorb,jorb,1)))
!!       write(*,'(i3,i3,1x,e21.14,1x,e21.14,1x,e9.2)') jorb,iorb,ovrlp(jorb,iorb,1),ovrlp(iorb,jorb,1),ovrlp(jorb,iorb,1)-ovrlp(iorb,jorb,1)
!       enddo ; enddo
!       write(*,*) 'iproc,overlap: max violation of symmetry:',iproc,tt

!       if (iproc.eq.0) then 
!        write(6,*) 'loewe S'
!        do i=1,norb
!!        write(6,'(14(1x,e6.1))') (abs(ovrlp(i,j,1)),j=1,min(14,norb))
!        write(6,'(i2,i3,14(1x,e21.14))') iproc,i,(abs(ovrlp(i,j,1)),j=1,min(14,norb))
!        enddo
!       endif

! LAPACK
        call DSYEV('V','L',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
        if (info.ne.0) write(6,*) 'info loewe', info
!        if (iproc.eq.0) then 
!          write(6,*) 'overlap eigenvalues'
!77        format(8(1x,e10.3))
!          if (norb.le.16) then
!          write(6,77) evall
!          else
!          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
!          endif
!        endif

! calculate S^{-1/2} ovrlp(*,*,3)
        do 2935,lorb=1,norb
        do 2935,jorb=1,norb
2935    ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
!        do 3985,j=1,norb
!        do 3985,i=1,norb
!        ovrlp(i,j,3)=0.d0
!        do 3985,l=1,norb
!3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
! BLAS:
        call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

! new eigenvectors
     do 200,iorb=1,norb
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1
       call wzero(mavctr_c,mavctr_f,psi(ipsi_c),psi(ipsi_f))

         do 200,jorb=1,norb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

        scpr=ovrlp(jorb,iorb,3)
	call waxpy(  & 
                   scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),tpsi(jpsi_c),tpsi(jpsi_f), &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f))
200     continue

        deallocate(ovrlp,evall,tpsi)

 endif

        return
        end



	subroutine checkortho(iproc,toler,norb,nvctr,nseg,keyg,keyv,psi)
! checks whether loewe worked correctly
        implicit real*8 (a-h,o-z)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))

        write(*,*) 'Start checking orthogonality',iproc
        dev=0.d0
     do 100,iorb=1,norb
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1

         do 100,jorb=1,norb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

	call wdot(  & 
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),  & 
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),psi(jpsi_c),psi(jpsi_f),scpr)
        if (iorb.eq.jorb) then
        dev=dev+(scpr-1.d0)**2
        else
        dev=dev+scpr**2
        endif
!        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
!        if (iorb.ne.jorb .and. abs(scpr).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
100    continue
        write(*,*) 'Deviation from orthogonality ',iproc,dev

        return
	end
