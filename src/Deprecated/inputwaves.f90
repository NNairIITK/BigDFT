

        subroutine crtinpwave(iproc,nproc,n1,n2,n3,nbox_c,nbox_f,norb,  & 
                              wan_par,nseg,nvctr,keyg,keyv,hgrid,psi)
! calls the routine that generates the input guess for one orbital
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        dimension nbox_c(2,3,norb),nbox_f(2,3,norb),wan_par(0:3,norb)

       onem=1.d0-eps_mach
       norb_p=int(onem+dble(norb)/dble(nproc))
    do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)

         rx=wan_par(1,iorb)
         ry=wan_par(2,iorb)
         rz=wan_par(3,iorb)
         gau_a=wan_par(0,iorb)

         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1_c=nbox_c(1,1,iorb) ; nu1_c=nbox_c(2,1,iorb)
         nl2_c=nbox_c(1,2,iorb) ; nu2_c=nbox_c(2,2,iorb)
         nl3_c=nbox_c(1,3,iorb) ; nu3_c=nbox_c(2,3,iorb)
         nl1_f=nbox_f(1,1,iorb) ; nu1_f=nbox_f(2,1,iorb)
         nl2_f=nbox_f(1,2,iorb) ; nu2_f=nbox_f(2,2,iorb)
         nl3_f=nbox_f(1,3,iorb) ; nu3_f=nbox_f(2,3,iorb)

          call crtonewave(n1,n2,n3, & 
          nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
          mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),  &
          hgrid,gau_a,rx,ry,rz,psi(ipsi_c),psi(ipsi_f))
     enddo

        return
        end subroutine crtinpwave


        subroutine crtonewave(n1,n2,n3, & 
                   nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                   mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  &
                   hgrid,gau_a,rx,ry,rz,psi_c,psi_f)
! returns an input guess orbital that is a Gaussian centered at a Wannier center
! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays psi_c, psi_f
        implicit real*8 (a-h,o-z)
        parameter(nw=16000)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work
        real*8, allocatable :: psig_c(:,:,:), psig_f(:,:,:,:)

        allocate(wprojx(0:n1,2),wprojy(0:n2,2),wprojz(0:n3,2),work(0:nw,2))
        allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c))
        allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f))

        n_gau=0
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

! First term: coarse projector components
          do i3=nl3_c,nu3_c
          do i2=nl2_c,nu2_c
          do i1=nl1_c,nu1_c
            psig_c(i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
          enddo ; enddo ; enddo

! First term: fine projector components
          do i3=nl3_f,nu3_f
          do i2=nl2_f,nu2_f
          do i1=nl1_f,nu1_f
            psig_f(1,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
            psig_f(2,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(3,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(4,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(5,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(6,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
            psig_f(7,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
          enddo ; enddo ; enddo

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
            psi_c(i-i0+jj)=psig_c(i,i2,i3)
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
            psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
            psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
            psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
            psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
            psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
            psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
            psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
          enddo
        enddo
  
          deallocate(wprojx,wprojy,wprojz,work,psig_c,psig_f)

	return
	end subroutine crtonewave
