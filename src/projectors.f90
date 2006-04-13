	

        subroutine numb_proj(ityp,ntypes,psppar,mproj)
! Determines the number of projectors
        implicit real*8 (a-h,o-z)
        dimension psppar(0:2,0:4,ntypes)

        mproj=0
! ONLY GTH PSP form (not HGH)
          do i=1,2
          do j=1,2
            if (psppar(i,j,ityp).ne.0.d0) mproj=mproj+2*i-1
          enddo
          enddo

        return
        end


        subroutine applyprojectors(iproc,nproc,ntypes,nat,iatype,psppar,occup, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norb,nseg,keyg,keyv,nvctr,psi,hpsi,eproj_sum)
! Applies all the projectors onto a wavefunction
! Input: psi_c,psi_f
! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension psppar(0:2,0:4,ntypes),iatype(nat)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb))
        dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
        dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
        dimension proj(nprojel),occup(norb)

  eproj_sum=0.d0
  onem=1.d0-eps_mach
  norb_p=int(onem+dble(norb)/dble(nproc))
! loop over all my orbitals
  do 100,iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
    maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
    maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
    iseg_c=nseg(2*iorb-2)+1
    iseg_f=nseg(2*iorb-1)+1
    mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
    mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
    ipsi_c=nvctr(2*iorb-2)+1
    ipsi_f=nvctr(2*iorb-1)+1

! loop over all projectors
    iproj=0
    eproj=0.d0
    istart_c=1
    do iat=1,nat
        mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
        mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
        jseg_c=nseg_p(2*iat-2)+1
        jseg_f=nseg_p(2*iat-1)+1
        mbvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
        mbvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
     ityp=iatype(iat)
! ONLY GTH PSP form (not HGH)
     do i=1,2
     do j=1,2
     if (psppar(i,j,ityp).ne.0.d0) then
     do l=1,2*i-1
        iproj=iproj+1
        istart_f=istart_c+mbvctr_c
        call wdot(  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  &
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),  &
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),scpr)
! test
        call wdot(  &
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  &
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),tcpr)
        if (scpr.ne.tcpr) stop 'projectors: scpr.ne.tcpr'
! testend

        scprp=scpr*psppar(i,j,ityp)
!        write(*,*) 'pdot scpr',scpr,psppar(i,j,ityp)
        eproj=eproj+scprp*scpr

        call waxpy(  &
             scprp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  &
                   keyg(1,iseg_c),keyg(1,iseg_f),hpsi(ipsi_c),hpsi(ipsi_f))

        istart_c=istart_f+7*mbvctr_f
     enddo
     endif
     enddo
     enddo
   enddo
     write(*,*) 'iorb,eproj',iorb,eproj
     eproj_sum=eproj_sum+occup(iorb)*eproj
     if (iproj.ne.nproj) stop '1:applyprojectors'
     if (istart_c-1.ne.nprojel) stop '2:applyprojectors'
100 continue

         return
         end



        subroutine crtproj(iproc,nterm,n1,n2,n3, & 
                     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                     radius_f,cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                     mvctr_c,mvctr_f,proj_c,proj_f)
! returns the compressed form of a Gaussian projector 
! x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays proj_c, proj_f
        implicit real*8 (a-h,o-z)
        parameter(ntermx=3,nw=16000)
        parameter(eps_mach=1.d-12)
        dimension lx(3),ly(3),lz(3)
        dimension proj_c(mvctr_c),proj_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work

        allocate(wprojx(0:n1,2,ntermx),wprojy(0:n2,2,ntermx),wprojz(0:n3,2,ntermx),work(0:nw,2))


        onem=1.d0-eps_mach
        rad_c=radius_f*cpmult
        rad_f=radius_f*fpmult

! make sure that the coefficients returned by CALL GAUSS_TO_DAUB are zero outside [ml:mr] 
        err_norm=0.d0 
      do 100,iterm=1,nterm
        n_gau=lx(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,factor,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
        n_gau=ly(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
        n_gau=lz(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
       if (ml1.gt.min(nl1_c,nl1_f)) write(*,*) 'ERROR: ml1'
       if (ml2.gt.min(nl2_c,nl2_f)) write(*,*) 'ERROR: ml2'
       if (ml3.gt.min(nl3_c,nl3_f)) write(*,*) 'ERROR: ml3'
       if (mu1.lt.max(nu1_c,nu1_f)) write(*,*) 'ERROR: mu1'
       if (mu2.lt.max(nu2_c,nu2_f)) write(*,*) 'ERROR: mu2'
       if (mu3.lt.max(nu3_c,nu3_f)) write(*,*) 'ERROR: mu3'
100   continue
        if (iproc.eq.0) write(*,*) 'max err_norm ',err_norm

! First term: coarse projector components
!          write(*,*) 'rad_c=',rad_c
          mvctr=0
          do i3=nl3_c,nu3_c
          dz2=(i3*hgrid-rz)**2
          do i2=nl2_c,nu2_c
          dy2=(i2*hgrid-ry)**2
          do i1=nl1_c,nu1_c
          dx=i1*hgrid-rx
          if (dx**2+(dy2+dz2).lt.rad_c**2) then
            mvctr=mvctr+1
            proj_c(mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
          endif
          enddo ; enddo ; enddo
          if (mvctr.ne.mvctr_c) then 
            write(*,*)  'mvctr,mvctr_c',mvctr,mvctr_c
            stop 'mvctr >< mvctr_c'
          endif

! First term: fine projector components
          mvctr=0
          do i3=nl3_f,nu3_f
          dz2=(i3*hgrid-rz)**2
          do i2=nl2_f,nu2_f
          dy2=(i2*hgrid-ry)**2
          do i1=nl1_f,nu1_f
          dx=i1*hgrid-rx
          if (dx**2+(dy2+dz2).lt.rad_f**2) then
            mvctr=mvctr+1
            proj_f(1,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
            proj_f(2,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
            proj_f(3,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
            proj_f(4,mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
            proj_f(5,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
            proj_f(6,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
            proj_f(7,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
          endif
          enddo ; enddo ; enddo
          if (mvctr.ne.mvctr_f) stop 'mvctr >< mvctr_f'
  

         do iterm=2,nterm

! Other terms: coarse projector components
         mvctr=0
         do i3=nl3_c,nu3_c
         dz2=(i3*hgrid-rz)**2
         do i2=nl2_c,nu2_c
         dy2=(i2*hgrid-ry)**2
         do i1=nl1_c,nu1_c
         dx=i1*hgrid-rx
         if (dx**2+(dy2+dz2).lt.rad_c**2) then
           mvctr=mvctr+1
           proj_c(mvctr)=proj_c(mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
         endif
         enddo ; enddo ; enddo

! Other terms: fine projector components
         mvctr=0
         do i3=nl3_f,nu3_f
         dz2=(i3*hgrid-rz)**2
         do i2=nl2_f,nu2_f
         dy2=(i2*hgrid-ry)**2
         do i1=nl1_f,nu1_f
         dx=i1*hgrid-rx
         if (dx**2+(dy2+dz2).lt.rad_f**2) then
           mvctr=mvctr+1
           proj_f(1,mvctr)=proj_f(1,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
           proj_f(2,mvctr)=proj_f(2,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
           proj_f(3,mvctr)=proj_f(3,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
           proj_f(4,mvctr)=proj_f(4,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
           proj_f(5,mvctr)=proj_f(5,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
           proj_f(6,mvctr)=proj_f(6,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
           proj_f(7,mvctr)=proj_f(7,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
         endif
         enddo ; enddo ; enddo
          

          enddo
  
          deallocate(wprojx,wprojy,wprojz,work)

	return
	end
