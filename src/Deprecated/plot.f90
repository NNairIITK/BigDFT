
        subroutine plot_wf(n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, & 
                                        psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f)
        real*8, allocatable, dimension(:) :: ww,psifscf,psir,psig

! Work arrays dimensioned large enough that they can be used in different contexts
        allocate(ww((2*n1+31)*(2*n2+31)*(2*n3+31)))
! Array that holds both scaling function and wavelet coefficients
        allocate(psig((2*n1+2)*(2*n2+2)*(2*n3+2))) 
! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psifscf((2*n1+16)*(2*n2+16)*(2*n3+16)))
! Wavefunction in real space
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))

        call uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, & 
                        ww,psig,psi_c,psi_f,psifscf)
        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,ww,psifscf,psir) 

        call plot_pot(hgrid,n1,n2,n3,10,psir)

        deallocate(ww,psifscf,psir,psig)
        return
	end




        subroutine plot_pot(hgrid,n1,n2,n3,iounit,pot)
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=.5d0*hgrid
         rx=8.d0
         ry=8.d0
         rz=8.d0
        open(iounit) 
        open(iounit+1) 
        open(iounit+2) 

        i3=nint(rz/hgridh)
        i2=nint(ry/hgridh)
        write(*,*) 'plot_p, i2,i3 ',i2,i3
        do i1=-14,2*n1+16
        write(iounit,*) i1*hgridh,pot(i1,i2,i3)
        enddo

        i1=nint(rx/hgridh)
        i2=nint(ry/hgridh)
        write(*,*) 'plot_p, i1,i2 ',i1,i2
        do i3=-14,2*n3+16
        write(iounit+1,*) i3*hgridh,pot(i1,i2,i3)
        enddo

        i1=nint(rx/hgridh)
        i3=nint(rz/hgridh)
        write(*,*) 'plot_p, i1,i3 ',i1,i3
        do i2=-14,2*n2+16
        write(iounit+2,*) i2*hgridh,pot(i1,i2,i3)
        enddo

        close(iounit) 
        close(iounit+1) 
        close(iounit+2) 

        return
        end



        subroutine plot_psifscf(iunit,hgrid,nl1,nu1,nl2,nu2,nl3,nu3,psifscf)
        implicit real*8 (a-h,o-z)
        dimension psifscf(-7+nl1:2*nu1+8,-7+nl2:2*nu2+8,-7+nl3:2*nu3+8)

	hgridh=.5d0*hgrid

! along x-axis
	i3=n3
	i2=n2
	do i1=-7+nl1:2*nu1+8
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! 111 diagonal
	do i=max(-7+nl1,-7+nl2,-7+nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
        i1=i ; i2=i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! 1-1-1 diagonal
	do i=max(-7+nl1,-7+nl2,-7+nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
        i1=i ; i2=-i ; i3=-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! -11-1 diagonal
	do i=max(-7+nl1,-7+nl2,-7+nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
        i1=-i ; i2=i ; i3=-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! -1-11 diagonal
	do i=max(-7+nl1,-7+nl2,-7+nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
        i1=-i ; i2=-i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

        return
        end

