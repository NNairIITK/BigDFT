

        subroutine readmywaves(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,rxyz,nseg,nvctr,keyg,keyv,psi)
! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
        implicit real*8 (a-h,o-z)
        character*30 filename
        character*4 f4
        logical cif1,cif2,cif3
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        dimension xya(-1:1,-1:1),xa(-1:1)
        dimension rxyz(3),rxyz_old(3),nbox_c(2,3,norb)
        allocatable :: psifscf(:,:,:)
        allocatable :: psigold(:,:,:,:,:,:),psifscfold(:,:,:),psifscfoex(:,:,:),wwold(:)

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

     do 100,iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
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


        write(f4,'(i4.4)') iorb
        filename = 'wavefunction.'//f4
        open(unit=99,file=filename,status='unknown')

         read(99,*) iorbold
         if (iorbold.ne.iorb) stop 'readallwaves'
         read(99,*) hgridold
         read(99,*) nold1,nold2,nold3
         read(99,*) nlold1,nuold1,nlold2,nuold2,nlold3,nuold3
         read(99,*) (rxyz_old(j),j=1,3)
         read(99,*) mvctrold_c, mvctrold_f

      if (hgridold.eq. hgrid .and. mvctrold_c.eq.mvctr_c .and. mvctrold_f.eq.mvctr_f  & 
          .and. nold1.eq.n1  .and. nold2.eq.n2 .and. nold3.eq.n3  & 
          .and. nlold1.eq.nl1 .and. nuold1.eq.nu1  & 
          .and. nlold2.eq.nl2 .and. nuold2.eq.nu2  &  
          .and. nlold3.eq.nl3 .and. nuold3.eq.nu3) then

      if (iproc.eq.0) write(*,*) 'wavefunction ',iorb,' needs NO transformation'
	do j=1,mvctrold_c
            read(99,*) i1,i2,i3,tt
            psi(ipsi_c+j-1)=tt
         enddo
	do j=1,7*mvctrold_f-6,7
            read(99,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            psi(ipsi_f+j-1+0)=t1
            psi(ipsi_f+j-1+1)=t2
            psi(ipsi_f+j-1+2)=t3
            psi(ipsi_f+j-1+3)=t4
            psi(ipsi_f+j-1+4)=t5
            psi(ipsi_f+j-1+5)=t6
            psi(ipsi_f+j-1+6)=t7
         enddo

       else
       if (iproc.eq.0) write(*,*) 'wavefunction ',iorb,' is transformed'

         allocate(psifscf(-7+2*nl1:2*nu1+8,-7+2*nl2:2*nu2+8,-7+2*nl3:2*nu3+8))

         allocate(psigold(nlold1:nuold1,2,nlold2:nuold2,2,nlold3:nuold3,2),  & 
               psifscfold(-7+2*nlold1:2*nuold1+8,-7+2*nlold2:2*nuold2+8,-7+2*nlold3:2*nuold3+8), &
               psifscfoex(-8+2*nlold1:2*nuold1+9,-8+2*nlold2:2*nuold2+9,-8+2*nlold3:2*nuold3+9), &
                   wwold((2*(nuold1-nlold1)+16)*(2*(nold2-nlold2)+16)*(2*(nold3-nlold3)+16)))

         call zero(8*(nuold1-nlold1+1)*(nuold2-nlold2+1)*(nuold3-nlold3+1),psigold)

	do iel=1,mvctrold_c
            read(99,*) i1,i2,i3,tt
            psigold(i1,1,i2,1,i3,1)=tt
         enddo
	do iel=1,mvctrold_f
            read(99,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            psigold(i1,2,i2,1,i3,1)=t1
            psigold(i1,1,i2,2,i3,1)=t2
            psigold(i1,2,i2,2,i3,1)=t3
            psigold(i1,1,i2,1,i3,2)=t4
            psigold(i1,2,i2,1,i3,2)=t5
            psigold(i1,1,i2,2,i3,2)=t6
            psigold(i1,2,i2,2,i3,2)=t7
         enddo

! calculate fine scaling functions
	call synthese_grow(nuold1-nlold1,nuold2-nlold2,nuold3-nlold3,wwold,psigold,psifscfold)

         do i3=-7+2*nlold3,2*nuold3+8
         do i2=-7+2*nlold2,2*nuold2+8
           i1=-8+2*nlold1
           psifscfoex(i1,i2,i3)=0.d0
           do i1=-7+2*nlold1,2*nuold1+8
             psifscfoex(i1,i2,i3)=psifscfold(i1,i2,i3)
           enddo
           i1=2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

         i3=-8+2*nlold3
         do i2=-8+2*nlold2,2*nuold2+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo
         i3=2*nuold3+9
         do i2=-8+2*nlold2,2*nuold2+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

         i2=-8+2*nlold2
         do i3=-8+2*nlold3,2*nuold3+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo
         i2=2*nuold2+9
         do i3=-8+2*nlold3,2*nuold3+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

! transform to new structure	
        dx=rxyz(1)-rxyz_old(1)
        dy=rxyz(2)-rxyz_old(2)
        dz=rxyz(3)-rxyz_old(3)
        write(*,*) 'dxyz',dx,dy,dz
        hgridh=.5d0*hgrid
        hgridhold=.5d0*hgridold
        call zero((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16),psifscf)
        do i3=-7+2*nl3,2*nu3+8
        z=i3*hgridh
        j3=nint((z-dz)/hgridhold)
        cif3=(j3.ge.-7+2*nlold3 .and. j3.le.2*nuold3+8)
        do i2=-7+2*nl2,2*nu2+8
        y=i2*hgridh
        j2=nint((y-dy)/hgridhold)
        cif2=(j2.ge.-7+2*nlold2 .and. j2.le.2*nuold2+8)
        do i1=-7+2*nl1,2*nu1+8
        x=i1*hgridh
        j1=nint((x-dx)/hgridhold)
        cif1=(j1.ge.-72*nlold1 .and. j1.le.2*nuold1+8)

!        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfold(j1,j2,j3)
!        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfoex(j1,j2,j3)

        if (cif1 .and. cif2 .and. cif3) then 
        zr = ((z-dz)-j3*hgridhold)/hgridhold
        do l2=-1,1
        do l1=-1,1
        ym1=psifscfoex(j1+l1,j2+l2,j3-1)
        y00=psifscfoex(j1+l1,j2+l2,j3  )
        yp1=psifscfoex(j1+l1,j2+l2,j3+1)
        xya(l1,l2)=ym1 + (1.d0 + zr)*(y00 - ym1 + zr*(.5d0*ym1 - y00  + .5d0*yp1))
        enddo
        enddo

        yr = ((y-dy)-j2*hgridhold)/hgridhold
        do l1=-1,1
        ym1=xya(l1,-1)
        y00=xya(l1,0)
        yp1=xya(l1,1)
        xa(l1)=ym1 + (1.d0 + yr)*(y00 - ym1 + yr*(.5d0*ym1 - y00  + .5d0*yp1))
        enddo

        xr = ((x-dx)-j1*hgridhold)/hgridhold
        ym1=xa(-1)
        y00=xa(0)
        yp1=xa(1)
        psifscf(i1,i2,i3)=ym1 + (1.d0 + xr)*(y00 - ym1 + xr*(.5d0*ym1 - y00  + .5d0*yp1))

        endif

        enddo
        enddo
        enddo


        call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psifscf,psi(ipsi_c),psi(ipsi_f))

         deallocate(psigold,psifscfold,psifscfoex,wwold)
         deallocate(psifscf)

      endif

        close(99)
100  continue


	end


        subroutine writemywaves(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,  & 
                   rxyz,nseg,nvctr,keyg,keyv,psi)
! write all my wavefunctions in files by calling writeonewave
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nbox_c(2,3,norb),rxyz(3)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))

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

       call writeonewave(iorb,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,hgrid,rxyz,  & 
                         mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c)  & 
                        ,mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f), & 
                             psi(ipsi_c),psi(ipsi_f))
       enddo
       return
       end



        subroutine writeonewave(iorb,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,hgrid,rxyz,  & 
                           mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                           mseg_f,mvctr_f,keyg_f,keyv_f, & 
                              psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        character*30 filename
        character*4 f4
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f),rxyz(3)


        write(f4,'(i4.4)') iorb
        filename = 'wavefunction.'//f4
        write(*,*) 'opening ',filename
        open(unit=99,file=filename,status='unknown')
         write(99,*) iorb
         write(99,*) hgrid
         write(99,*) n1,n2,n3
         write(99,*) nl1,nu1,nl2,nu2,nl3,nu3
         write(99,'(3(1x,e24.17))') (rxyz(j),j=1,3)
         write(99,*) mvctr_c, mvctr_f

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
            tt=psi_c(i-i0+jj) 
            write(99,'(3(i4),1x,e19.12)') i,i2,i3,tt
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
            t1=psi_f(1,i-i0+jj)
            t2=psi_f(2,i-i0+jj)
            t3=psi_f(3,i-i0+jj)
            t4=psi_f(4,i-i0+jj)
            t5=psi_f(5,i-i0+jj)
            t6=psi_f(6,i-i0+jj)
            t7=psi_f(7,i-i0+jj)
            write(99,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
          enddo
         enddo

          close(99)

	write(*,*) iorb,'th wavefunction written'


	end
