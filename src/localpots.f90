

	subroutine compresstest(iproc,nproc,n1,n2,n3,norb,nbox_c,nseg,nvctr,keyg,keyv,psi,hpsi)
! test compress-uncompress
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nbox_c(2,3,norb),nseg(0:2*norb),nvctr(0:2*norb)
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb))
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        allocatable :: psifscf(:)


! begin test compress-uncompress
       onem=1.d0-eps_mach
       norb_p=int(onem+dble(norb)/dble(nproc))
       do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

        allocate(psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )

	do i=1,mvctr_c
        call random_number(tt)
        psi(ipsi_c+i-1)=tt
        hpsi(ipsi_c+i-1)=tt
        enddo

	do i=1,7*mvctr_f
        call random_number(tt)
        psi(ipsi_f+i-1)=tt
        hpsi(ipsi_f+i-1)=tt
        enddo

        call uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psi(ipsi_c),psi(ipsi_f),psifscf)
        call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psifscf,psi(ipsi_c),psi(ipsi_f))

        tc=0.d0
	do i=1,mvctr_c
        tc=max(tc,abs(psi(ipsi_c+i-1)-hpsi(ipsi_c+i-1)))
        enddo
        if (tc.gt.1.d-10) stop 'coarse compress error'

        tf=0.d0
	do i=1,7*mvctr_f
        tf=max(tf,abs(psi(ipsi_f+i-1)-hpsi(ipsi_f+i-1))) 
        enddo
        write(*,'(a,i4,2(1x,e9.2))') 'COMPRESSION TEST: iorb, coarse, fine error:',iorb,tc,tf
        if (tf.gt.1.d-10) stop 'fine compress error'

        deallocate(psifscf)

     enddo

	return
	end subroutine compresstest


        subroutine compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                            mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                            mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                            psifscf,psi_c,psi_f)
! Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16))
        real*8, allocatable :: psig(:,:,:,:,:,:),ww(:)

        allocate(psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2),  & 
                 ww((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)))
! decompose wavelets into coarse scaling functions and wavelets
	call analyse_shrink(nu1-nl1,nu2-nl2,nu3-nl3,ww,psifscf,psig)

! coarse part
	do iseg=1,mseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
	  do i=i0,i1
            psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
          enddo
        enddo

! fine part
	do iseg=1,mseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
	  do i=i0,i1
            psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
            psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
            psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
            psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
            psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
            psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
            psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
          enddo
        enddo

        deallocate(psig,ww)

	end subroutine compress


        subroutine uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                              mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                              mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                              psi_c,psi_f,psifscf)
! Expands the compressed wavefunction in vector form (psi_c,psi_f) 
! into fine scaling functions (psifscf)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16))
        real*8, allocatable :: psig(:,:,:,:,:,:),ww(:)

        allocate(psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2),  & 
                 ww((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)))

        call zero(8*(nu1-nl1+1)*(nu2-nl2+1)*(nu3-nl3+1),psig)

! coarse part
	do iseg=1,mseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
	  do i=i0,i1
            psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
          enddo
         enddo

! fine part
	do iseg=1,mseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
	  do i=i0,i1
            psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
            psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
            psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
            psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
            psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
            psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
            psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
          enddo
         enddo

! calculate fine scaling functions
	call synthese_grow(nu1-nl1,nu2-nl2,nu3-nl3,ww,psig,psifscf)

        deallocate(psig,ww)

	end subroutine uncompress


        subroutine applylocpotkin(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
                   nseg,nvctr,keyg,keyv,psi,pot,hpsi,epot_sum,ekin_sum)
!  Applies the local potential and kinetic energy operator to a wavefunction
! Input: pot,psi
! Output: hpsi,epot,ekin
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nbox_c(2,3,norb),occup(norb),pot((2*n1+31)*(2*n2+31)*(2*n3+31))
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb))
        real*8, allocatable, dimension(:) :: psifscf,psifscfk,psir
        include 'mpif.h'

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

        hgridh=hgrid*.5d0

! Determine aximal size of work arrays
      nl1=10000 ; nu1=0
      nl2=10000 ; nu2=0
      nl3=10000 ; nu3=0
      do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
        nl1=min(nl1,nbox_c(1,1,iorb)) ; nu1=max(nu1,nbox_c(2,1,iorb))
        nl2=min(nl2,nbox_c(1,2,iorb)) ; nu2=max(nu2,nbox_c(2,2,iorb))
        nl3=min(nl3,nbox_c(1,3,iorb)) ; nu3=max(nu3,nbox_c(2,3,iorb))
      enddo
! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )
        allocate(psifscfk((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )
! Wavefunction in real space
        allocate(psir((2*(nu1-nl1)+31)*(2*(nu2-nl2)+31)*(2*(nu3-nl3)+31)))


        ekin_sum=0.d0
        epot_sum=0.d0
     do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

        call uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  &
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   &
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   &
                    psi(ipsi_c),psi(ipsi_f),psifscf)

        call convolut_kinetic(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,hgridh,psifscf,psifscfk)

        ekin=0.d0 
        do i=1,(2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)
          ekin=ekin+psifscf(i)*psifscfk(i)
        enddo
        ekin_sum=ekin_sum+occup(iorb)*ekin 

        call convolut_magic_n(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psifscf,psir) 

        call  potloc(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,pot,psir,epot)
        epot_sum=epot_sum+occup(iorb)*epot
        write(*,'(a,i5,2(1x,e17.10))') 'iorb,ekin,epot',iorb,ekin,epot

        call convolut_magic_t(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psir,psifscf)

        do i=1,(2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)
          psifscf(i)=psifscf(i)+psifscfk(i)
        enddo

        call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psifscf,hpsi(ipsi_c),hpsi(ipsi_f))
     enddo

        deallocate(psifscf,psifscfk,psir)
   
        return
	end subroutine applylocpotkin


        subroutine potloc(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,pot,psir,epot)
! applies the potential to a localized orbital in its local box
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
        dimension psir(-14+2*nl1:2*nu1+16,-14+2*nl2:2*nu2+16,-14+2*nl3:2*nu3+16)

        epot=0.d0 
        do i3=-14+2*nl3,2*nu3+16
        do i2=-14+2*nl2,2*nu2+16
        do i1=-14+2*nl1,2*nu1+16
          tt=pot(i1,i2,i3)*psir(i1,i2,i3)
          epot=epot+tt*psir(i1,i2,i3)
          psir(i1,i2,i3)=tt
        enddo ; enddo ; enddo

!        do i3=-14+2*nl3,2*nu3+16
!        i1=i3 ; i2=i3 
!        write(77,*) i3,pot(i1,i2,i3),psir(i1,i2,i3)
!        enddo

	return
        end subroutine potloc



        subroutine convolut_magic_n(n1,n2,n3,x,y)
! Applies the magic filter matrix ( no transposition) ; data set grows
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
        dimension x(0:n1,0:n2,0:n3),y(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil)
        real*8, allocatable :: ww(:,:,:)

        allocate(ww(-lupfil:n1-lowfil,-lupfil:n2-lowfil,0:n3))
 
!  (i1,i2*i3) -> (i2*i3,I1)
        ndat=(n2+1)*(n3+1)
        call convrot_grow(n1,ndat,x,y)
!  (i2,i3*I1) -> (i3*i1,I2)
        ndat=(n3+1)*(n1+1+lupfil-lowfil)
        call convrot_grow(n2,ndat,y,ww)
!  (i3,I1*I2) -> (iI*I2,I3)
        ndat=(n1+1+lupfil-lowfil)*(n2+1+lupfil-lowfil)
        call convrot_grow(n3,ndat,ww,y)

        deallocate(ww)

	return
	end subroutine convolut_magic_n


        subroutine convolut_magic_t(n1,n2,n3,x,y)
! Applies the magic filter matrix transposed ; data set shrinks
! The input array x is overwritten
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
        dimension x(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil),y(0:n1,0:n2,0:n3)
        real*8, allocatable :: ww(:,:,:)

        allocate(ww(0:n1,-lupfil:n2-lowfil,-lupfil:n3-lowfil))

!  (I1,I2*I3) -> (I2*I3,i1)
        ndat=(n2+1+lupfil-lowfil)*(n3+1+lupfil-lowfil)
        call convrot_shrink(n1,ndat,x,ww)
!  (I2,I3*i1) -> (I3*i1,i2)
        ndat=(n3+1+lupfil-lowfil)*(n1+1)
        call convrot_shrink(n2,ndat,ww,x)
!  (I3,i1*i2) -> (i1*i2,i3)
        ndat=(n1+1)*(n2+1)
        call convrot_shrink(n3,ndat,x,y)

        deallocate(ww)

	return
	end subroutine convolut_magic_t


	subroutine synthese_grow(n1,n2,n3,ww,x,y)
! A synthesis wavelet transformation where the size of the data is allowed to grow
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension x(0:n1,2,0:n2,2,0:n3,2)
        dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
        dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

! i1,i2,i3 -> i2,i3,I1
        nt=(2*n2+2)*(2*n3+2)
        call  syn_rot_grow(n1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(2*n3+2)*(2*n1+16)
        call  syn_rot_grow(n2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(2*n1+16)*(2*n2+16)
        call  syn_rot_grow(n3,nt,ww,y)

	return
        end subroutine synthese_grow


        subroutine syn_rot_grow(n,nt,x,y)
        implicit real*8 (a-h,o-z)
!  Daubechies_16 symetric 
!        REAL*8 CH(-M:M)/0.d0,&
!        -0.0033824159510050025955D0,-0.00054213233180001068935D0,&
!        0.031695087811525991431D0,0.0076074873249766081919D0,&
!        -0.14329423835127266284D0,-0.061273359067811077843D0,&
!        0.48135965125905339159D0,0.77718575169962802862D0,0.36444189483617893676D0,&
!        -0.051945838107881800736D0,-0.027219029917103486322D0,&
!        0.049137179673730286787D0,0.0038087520138944894631D0,&
!        -0.014952258337062199118D0,-0.00030292051472413308126D0,&
!        0.0018899503327676891843D0&

        dimension x(0:2*n+1,nt),y(nt,-7:2*n+8)
        include 'sym_16.inc'

        do 200,l=1,nt

       DO I=-7,2*n+8
         y(l,I)=0.D0
       enddo

       DO I=0,n
         I2=2*I
         DO J=-M+1,M
           y(l,J+I2)=y(l,J+I2)+CH(J)*x(I,l)+CG(J)*x(n+1+I,l)
         ENDDO
       ENDDO

200     continue

        return
        end subroutine syn_rot_grow


	subroutine analyse_shrink(n1,n2,n3,ww,y,x)
! A analysis wavelet transformation where the size of the data is forced to shrink
! The input array y is overwritten
        implicit real*8 (a-h,o-z)
        dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
        dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
        dimension x(0:n1,2,0:n2,2,0:n3,2)

! I1,I2,I3 -> I2,I3,i1
        nt=(2*n2+16)*(2*n3+16)
        call  ana_rot_shrink(n1,nt,y,ww)
! I2,I3,i1 -> I3,i1,i2
        nt=(2*n3+16)*(2*n1+2)
        call  ana_rot_shrink(n2,nt,ww,y)
! I3,i1,i2 -> i1,i2,i3
        nt=(2*n1+2)*(2*n2+2)
        call  ana_rot_shrink(n3,nt,y,x)

	return
        end subroutine analyse_shrink



        subroutine ana_rot_shrink(n,nt,x,y)
        implicit real*8 (a-h,o-z)
!  Daubechies_8 symetric 
!        parameter(chm3=.230377813308896501d0, chm2=.714846570552915647d0, & 
!                  chm1=.630880767929858908d0, ch00=-.279837694168598542d-1, & 
!                  chp1=-.187034811719093084d0, chp2=.308413818355607636d-1, & 
!                  chp3=.328830116668851997d-1, chp4=-.105974017850690321d-1)


       dimension x(-7:2*n+8,nt),y(nt,0:2*n+1)
       INCLUDE 'sym_16.inc'

       do 200,l=1,nt

       DO I=0,n
         I2=2*I
         CI=0.D0
         DI=0.D0
         DO J=-M+1,M
           CI=CI+CHT(J)*x(J+I2,l)
           DI=DI+CGT(J)*x(J+I2,l)
         ENDDO
         y(l,I)=CI
         y(l,n+1+I)=DI
       ENDDO
       
200     continue

        return
        end subroutine ana_rot_shrink

