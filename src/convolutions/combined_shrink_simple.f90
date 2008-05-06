subroutine comb_rot_shrink_loc(ndat,x,y,icf,nfl,nfu,ib)
! In one dimension,    
! Applies the magic filter transposed, then analysis wavelet transformation.
! The size of the data is forced to shrink
    implicit real(kind=8) (a-h,o-z)
    integer,parameter:: lowfil2=-14,lupfil2=16
    dimension x(lowfil2+2*nfl:2*nfu+lupfil2,ndat),y(ndat,nfl:nfu)
    integer ib(2,ndat)
    include 'v.inc'

!    open(unit=10,file='simple_shrink.flop')
!    nflop=0
!    ! count the flops:
!    do j=1,ndat
!           do i=ib(1,j),ib(2,j)
!             do l=lowfil2+2*i,lupfil2+2*i
!                nflop=nflop+2
!             enddo
!           enddo
!    enddo

    ! the convolution itself:
    call system_clock(ncount0,ncount_rate,ncount_max)
    do j=1,ndat
       do i=ib(1,j),ib(2,j)
         ci=0.d0
         do l=lowfil2+2*i,lupfil2+2*i
           ci=ci+fil2(l-2*i,icf)*x(l,j)
         enddo
         y(j,i)=ci
       enddo
    enddo
       
!    call system_clock(ncount1,ncount_rate,ncount_max)
!    tel=dble(ncount1-ncount0)/dble(ncount_rate)

!    write(10,*) tel, 1.d-6*nflop/tel
  end subroutine comb_rot_shrink_loc

subroutine comb_rot_shrink_loc_1(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
! In one dimension,    
! Applies the magic filter transposed, then analysis wavelet transformation.
! The size of the data is forced to shrink
    implicit real(kind=8) (a-h,o-z)
    integer,parameter:: lowfil2=-14,lupfil2=16
    real(kind=8) x(-14:2*n1+16,-14:2*n2+16,         -14:2*n3+16) ! input
    real(kind=8) y(2,     -14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)! output
    integer ib(2, -14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
    include 'v.inc'

    nflop=0
!    open(unit=20,file='tree_shrink.flop')
!    call system_clock(ncount0,ncount_rate,ncount_max)

    do j2=-14+2*nfl2,2*nfu2+16
        do j3=-14+2*nfl3,2*nfu3+16
           do i=ib(1,j2,j3),ib(2,j2,j3)
             ci1=0.d0
             ci2=0.d0
             !nflop=nflop+(lupfil2-lowfil2+1)*2*2
             do l=lowfil2+2*i,lupfil2+2*i
               ci1=ci1+fil2(l-2*i,1)*x(l,j2,j3)
               ci2=ci2+fil2(l-2*i,2)*x(l,j2,j3)
             enddo
             y(1,j2,j3,i)=ci1
             y(2,j2,j3,i)=ci2
           enddo
        enddo
    enddo

    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)
    !write(20,*) tel, 1.d-6*nflop/tel

end


subroutine comb_rot_shrink_loc_2(ndat,x,y,nfl,nfu,ib)
! In one dimension,    
! Applies the magic filter transposed, then analysis wavelet transformation.
! The size of the data is forced to shrink
    implicit real(kind=8) (a-h,o-z)
    integer,parameter:: lowfil2=-14,lupfil2=16
    dimension x(2,lowfil2+2*nfl:2*nfu+lupfil2,ndat),y(2,2,ndat,nfl:nfu)
    integer ib(2,ndat)
    include 'v.inc'

    nflop=0
    !open(unit=20,file='tree_shrink.flop')
    !call system_clock(ncount0,ncount_rate,ncount_max)

    do j=1,ndat
       do i=ib(1,j),ib(2,j)
         ci11=0.d0
         ci12=0.d0
         ci21=0.d0
         ci22=0.d0

         !nflop=nflop+(lupfil2-lowfil2+1)*2*4
         do l=lowfil2+2*i,lupfil2+2*i
           ci11=ci11+fil2(l-2*i,1)*x(1,l,j)
           ci12=ci12+fil2(l-2*i,2)*x(1,l,j)
           ci21=ci21+fil2(l-2*i,1)*x(2,l,j)
           ci22=ci22+fil2(l-2*i,2)*x(2,l,j)
         enddo
         y(1,1,j,i)=ci11
         y(1,2,j,i)=ci12
         y(2,1,j,i)=ci21
         y(2,2,j,i)=ci22
       enddo
    enddo

    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)
    !write(20,*) tel, 1.d-6*nflop/tel
       
end
