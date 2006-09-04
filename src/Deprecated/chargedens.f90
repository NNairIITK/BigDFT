

        subroutine sumrho(parallel,iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
                              nseg,nvctr,keyg,keyv,psi,rho)
! Calculates the charge density by summing the square of all orbitals
! Input: psi
! Output: rho
        implicit real*8 (a-h,o-z)
        logical parallel
        integer count1,count2,count_rate,count_max
        parameter(eps_mach=1.d-12)
	dimension rho((2*n1+31)*(2*n2+31)*(2*n3+31)),occup(norb)
        dimension nbox_c(2,3,norb)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        real*8, allocatable, dimension(:) :: psifscf,psir,rho_p
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
! Wavefunction in real space
        allocate(psir((2*(nu1-nl1)+31)*(2*(nu2-nl2)+31)*(2*(nu3-nl3)+31)))

 if (parallel) then
        allocate(rho_p((2*n1+31)*(2*n2+31)*(2*n3+31)))
	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

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

        call convolut_magic_n(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psifscf,psir) 

        const=occup(iorb)/hgridh**3
	call addpartrho(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,const,psir,rho_p)

      enddo

        call cpu_time(tr0)
        call system_clock(count1,count_rate,count_max)
        call MPI_ALLREDUCE(rho_p,rho,(2*n1+31)*(2*n2+31)*(2*n3+31),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call cpu_time(tr1)
        call system_clock(count2,count_rate,count_max)
        tel=dble(count2-count1)/dble(count_rate)
        write(78,*) 'RHO: ALLREDUCE TIME',iproc,tr1-tr0,tel
        write(78,*) '---------------------------------------------'

        deallocate(rho_p)
 else

	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

     do iorb=1,norb
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

        call convolut_magic_n(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psifscf,psir) 

        const=occup(iorb)/hgridh**3
	call addpartrho(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,const,psir,rho)

     enddo
 endif

! Check
        tt=0.d0
        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         tt=tt+rho(i)
        enddo
        tt=tt*hgridh**3
	if (iproc.eq.0) write(*,*) 'Total charge from routine chargedens',tt,iproc


        deallocate(psifscf,psir)

        return
	end subroutine sumrho


	subroutine addpartrho(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,const,psir,rho)
! Adds the contribution of one orbital in its local box
        implicit real*8 (a-h,o-z)
        dimension rho(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
        dimension psir(-14+2*nl1:2*nu1+16,-14+2*nl2:2*nu2+16,-14+2*nl3:2*nu3+16)

        do i3=-14+2*nl3,2*nu3+16
        do i2=-14+2*nl2,2*nu2+16
        do i1=-14+2*nl1,2*nu1+16
         rho(i1,i2,i3)=rho(i1,i2,i3)+const*psir(i1,i2,i3)**2
        enddo ; enddo ; enddo

        return
        end subroutine addpartrho
