

        subroutine plot_psifscf(iunit,hgrid,nl1,nu1,nl2,nu2,nl3,nu3,psifscf)
! plots psifscf
        implicit real(kind=8) (a-h,o-z)
        dimension psifscf(-7+2*nl1:2*nu1+8,-7+2*nl2:2*nu2+8,-7+2*nl3:2*nu3+8)

        hgridh=.5d0*hgrid

! 111 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=i ; i2=i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

        nc1=2*(nu1-nl1)+1
        nc2=2*(nu2-nl2)+1
        nc3=2*(nu3-nl3)+1
! 1-1-1 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=i ; i2=nc2-i ; i3=nc3-i
            write(iunit,'(3(1x,e10.3),1x,e12.5))') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

! -11-1 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=nc1-i ; i2=i ; i3=nc3-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

! -1-11 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=nc1-i ; i2=nc2-i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

        return
        end
