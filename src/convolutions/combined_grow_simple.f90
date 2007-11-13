


subroutine comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
                             ,w2,w1,xc,xf,y,ibyz_c,ibzxx_c,ibxxyy_c,&
                             ibyz_f,ibzxx_f,ibxxyy_f,ibyyzz_r)
implicit none

    integer n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3

    real(kind=8),intent(in)::xc(0:n1,0:n2,0:n3),  xf(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)! input
    real(kind=8),intent(out)::y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)! output

    real(kind=8) w1(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)!work
    real(kind=8) w2(max((n3+1)*(2*n1+31)*(2*n2+31),&
                2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))) ! work

    integer,intent(in)::ibyz_c(2,0:n2,0:n3)
    integer,intent(in)::ibzxx_c(2,        0:n3,-14:2*n1+16)
    integer,intent(in)::ibxxyy_c(2,             -14:2*n1+16,-14:2*n2+16)

    integer,intent(in)::ibyz_f(2,nfl2:nfu2,nfl3:nfu3) 
    integer,intent(in)::ibzxx_f(2,          nfl3:nfu3,2*nfl1-14:2*nfu1+16)
    integer,intent(in)::ibxxyy_f(2,                    2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

    integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)! boundaries of the real space array

    call comb_grow_c(n1,n2,n3,w2,xc,y,ibyz_c,ibzxx_c,ibxxyy_c,ibyyzz_r)

    call comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
                              ,w1,w2,xf,y,ibyz_f,ibzxx_f,ibxxyy_f)                

end subroutine comb_grow_all


        subroutine comb_grow_c(n1,n2,n3,ww,x,y,ibyz,ibzxx,ibxxyy,ibyyzz_r)
! In 3d,            
! Applies synthesis wavelet transformation 
! then convolves with magic filter
!  the size of the data is allowed to grow
! The input array x is not overwritten
! However, the output array y contains nonphysical values
! outside of the localization region
! that remain from the first comb_grow

        implicit real(kind=8) (a-h,o-z)
        real(kind=8) x(0:n1,0:n2,0:n3)
        real(kind=8) ww(         0:n3,-14:2*n1+16,-14:2*n2+16) ! work
        real(kind=8) y(               -14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
        integer ibyz  (    2,0:n2,0:n3)
        integer ibzxx(       2,       0:n3,-14:2*n1+16)
        integer ibxxyy(2       ,-14:2*n1+16,-14:2*n2+16)
    integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)! boundaries of the real space array
        
! i1,i2,i3 -> i2,i3,I1
        call comb_rot_grow_loc_square_1(n1,n2,n3,x,y,ibyz,ibzxx,.true.) 

! i2,i3,I1 -> i3,I1,I2
        call comb_rot_grow_loc_square_1(n2,n3,2*n1+30,y,ww,ibzxx,ibxxyy,.true.) 

! i3,I1,I2  -> I1,I2,I3
        call comb_rot_grow_loc_square_1(n3,2*n1+30,2*n2+30,ww,y,ibxxyy,ibyyzz_r,.true.) 

        END SUBROUTINE



subroutine comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
                             ,w1,w2,x,y,ibyz,ibzxx,ibxxyy)
! In 3d,            
! Applies synthesis wavelet transformation 
! then convolves with magic filter
!  the size of the data is allowed to grow
! The input array x is not overwritten

implicit real(kind=8) (a-h,o-z)

real(kind=8) x(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
real(kind=8) w1(4,           nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)
real(kind=8) w2(2,                      nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
real(kind=8) y(                                      -14:2*n1+16  ,       -14:2*n2+16  ,-14:2*n3+16)
integer ibyz(2,nfl2:nfu2,nfl3:nfu3)
integer ibzxx( 2,        nfl3:nfu3,2*nfl1-14:2*nfu1+16)
integer ibxxyy(2,                  2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

          m1=nfu1-nfl1;    m2=nfu2-nfl2;    m3=nfu3-nfl3

! i1,i2,i3 -> i2,i3,I1
        nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
        call  comb_rot_grow_loc_1(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,w1,ibyz,ibzxx) 

! i2,i3,I1 -> i3,I1,I2
        nt=(nfu3-nfl3+1)*(2*m1+31)
        call  comb_rot_grow_loc_2(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,ibzxx,ibxxyy) 

! i3,I1,I2  -> I1,I2,I3: add the result to y
        call  comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w2,y,ibxxyy)

        END SUBROUTINE



    subroutine  comb_rot_grow_loc_1(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib,ib2)
! In one dimesnion,    
! with optimised cycles
! Applies synthesis wavelet transformation 
! then convolves with magic filter
!  the size of the data is allowed to grow

    implicit none

    integer,intent(in)::nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
    integer,intent(in)::ib(2,         nfl2:nfu2,nfl3:nfu3)
    real(kind=8),intent(in )::x(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)

    integer,intent(in)::ib2(2,                    nfl3:nfu3,-14+2*nfl1:2*nfu1+16)

    
    real(kind=8),intent(out)::y(2,2,        nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)

    integer ncount0,ncount1,ncount_rate,ncount_max,nflop
    integer l2,l3,i,t,l1

    real(kind=8) tel
    real(kind=8) y2i__11,y2i__21,y2i1_11,y2i1_21
    real(kind=8) y2i__12,y2i__22,y2i1_12,y2i1_22

    include 'v_17.f90'


    open(unit=20,file='tree.flop')

    nflop=0
    do l2=nfl2,nfu2
        do l3=nfl3,nfu3
            if (ib(2,l2,l3).ge.ib(1,l2,l3)) nflop=nflop+(ib(2,l2,l3)-ib(1,l2,l3)+1)*31*2*7
        enddo
    enddo

    call system_clock(ncount0,ncount_rate,ncount_max)

!   y=0.d0
    do l1=-14+2*nfl1,2*nfu1+16
        do l3=nfl3,nfu3
            y(:,:,ib2(1,l3,l1):ib2(2,l3,l1),l3,l1)=0.d0
        enddo
    enddo
    
    do l3=nfl3,nfu3
    do l2=nfl2,nfu2

        if (ib(1,l2,l3).le.ib(2,l2,l3)) then

            y(1,1,l2,l3,2*ib(2,l2,l3)+16)=    fil2(16,2)*x(1,ib(2,l2,l3),l2,l3)
            y(2,1,l2,l3,2*ib(2,l2,l3)+16)=&
            fil2(16,1)*x(2,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(3,ib(2,l2,l3),l2,l3)
        
            y(1,2,l2,l3,2*ib(2,l2,l3)+16)=&
            fil2(16,1)*x(4,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(5,ib(2,l2,l3),l2,l3)
            y(2,2,l2,l3,2*ib(2,l2,l3)+16)=&
            fil2(16,1)*x(6,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(7,ib(2,l2,l3),l2,l3)

            do i=ib(1,l2,l3)-7,ib(2,l2,l3)+7 
                y2i__11=0.d0
                y2i__21=0.d0
                y2i__12=0.d0
                y2i__22=0.d0

                y2i1_11=0.d0
                y2i1_21=0.d0
                y2i1_12=0.d0
                y2i1_22=0.d0
                do t=max(i-8,ib(1,l2,l3)),min(i+7,ib(2,l2,l3))
                    y2i__11=y2i__11                               +fil2(2*(i-t)  ,2)*x(1,t,l2,l3)
                    y2i__21=y2i__21+fil2(2*(i-t)  ,1)*x(2,t,l2,l3)+fil2(2*(i-t)  ,2)*x(3,t,l2,l3)
                    y2i__12=y2i__12+fil2(2*(i-t)  ,1)*x(4,t,l2,l3)+fil2(2*(i-t)  ,2)*x(5,t,l2,l3)
                    y2i__22=y2i__22+fil2(2*(i-t)  ,1)*x(6,t,l2,l3)+fil2(2*(i-t)  ,2)*x(7,t,l2,l3)

                    y2i1_11=y2i1_11                               +fil2(2*(i-t)+1,2)*x(1,t,l2,l3)
                    y2i1_21=y2i1_21+fil2(2*(i-t)+1,1)*x(2,t,l2,l3)+fil2(2*(i-t)+1,2)*x(3,t,l2,l3)
                    y2i1_12=y2i1_12+fil2(2*(i-t)+1,1)*x(4,t,l2,l3)+fil2(2*(i-t)+1,2)*x(5,t,l2,l3)
                    y2i1_22=y2i1_22+fil2(2*(i-t)+1,1)*x(6,t,l2,l3)+fil2(2*(i-t)+1,2)*x(7,t,l2,l3)
                enddo
                y(1,1,l2,l3,2*i  )=y2i__11
                y(2,1,l2,l3,2*i  )=y2i__21
                y(1,2,l2,l3,2*i  )=y2i__12
                y(2,2,l2,l3,2*i  )=y2i__22

                y(1,1,l2,l3,2*i+1)=y2i1_11
                y(2,1,l2,l3,2*i+1)=y2i1_21
                y(1,2,l2,l3,2*i+1)=y2i1_12
                y(2,2,l2,l3,2*i+1)=y2i1_22
            enddo
        endif

    enddo
    enddo

    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)

    write(20,*) tel, 1.d-6*nflop/tel
end



    subroutine  comb_rot_grow_loc_2(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib,ib2)
! In one dimesnion,    
! with optimised cycles
! Applies synthesis wavelet transformation 
! then convolves with magic filter
!  the size of the data is allowed to grow

    implicit none

    integer,intent(in)::nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
    integer,intent(in)::ib( 2,          nfl3:nfu3,-14+2*nfl1:2*nfu1+16)
    real(kind=8),intent(in )::x(2,2,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)
    integer,intent(in)::ib2(2,                       -14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
    
    real(kind=8),intent(out)::y(2,            nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
    integer ncount0,ncount1,ncount_rate,ncount_max,nflop
    integer l1,l3,i,t,l2

    real(kind=8) tel,y2i__1,y2i__2,y2i1_1,y2i1_2

    include 'v_17.f90'

    open(unit=20,file='tree.flop')
    nflop=0
    do l3=nfl3,nfu3
        do l1=-14+2*nfl1,2*nfu1+16
            if (ib(2,l3,l1).ge.ib(1,l3,l1)) nflop=nflop+(ib(2,l3,l1)-ib(1,l3,l1)+1)
        enddo
    enddo
    nflop=nflop*31*2*4

    call system_clock(ncount0,ncount_rate,ncount_max)

!     y=0.d0
    do l2=-14+2*nfl2,2*nfu2+16
        do l1=-14+2*nfl1,2*nfu1+16
            y(:,ib2(1,l1,l2):ib2(2,l1,l2),l1,l2)=0.d0
        enddo
    enddo
     
    do l1=-14+2*nfl1,2*nfu1+16
    do l3=nfl3,nfu3

         if (ib(1,l3,l1).le.ib(2,l3,l1)) then
            y(1,l3,l1,2*ib(2,l3,l1)+16)=&
          fil2(16,1)*x(1,1,ib(2,l3,l1),l3,l1)+fil2(16,2)*x(2,1,ib(2,l3,l1),l3,l1)
          y(2,l3,l1,2*ib(2,l3,l1)+16)=&
          fil2(16,1)*x(1,2,ib(2,l3,l1),l3,l1)+fil2(16,2)*x(2,2,ib(2,l3,l1),l3,l1)
      
          do i=ib(1,l3,l1)-7,ib(2,l3,l1)+7 
              y2i__1=0.d0
              y2i__2=0.d0
              y2i1_1=0.d0
              y2i1_2=0.d0
              do t=max(i-8,ib(1,l3,l1)),min(i+7,ib(2,l3,l1))
                  y2i__1=y2i__1+fil2(2*(i-t)  ,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,1,t,l3,l1)
                  y2i__2=y2i__2+fil2(2*(i-t)  ,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,2,t,l3,l1)
                  y2i1_1=y2i1_1+fil2(2*(i-t)+1,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,1,t,l3,l1)
                  y2i1_2=y2i1_2+fil2(2*(i-t)+1,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,2,t,l3,l1)
              enddo
              y(1,l3,l1,2*i  )=y2i__1
              y(2,l3,l1,2*i  )=y2i__2
              y(1,l3,l1,2*i+1)=y2i1_1
              y(2,l3,l1,2*i+1)=y2i1_2
          enddo
        endif

    enddo
    enddo

    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)
    write(20,*) tel, 1.d-6*nflop/tel
end

    subroutine  comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
! In one dimesnion,    
! with optimised cycles
! Applies synthesis wavelet transformation 
! then convolves with magic filter
! then adds the result to y.
! The size of the data is allowed to grow

    implicit none

    integer,intent(in)::nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1,n2,n3
    integer,intent(in)::ib(2,        -14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
    real(kind=8),intent(in)::x(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
    real(kind=8),intent(out)::y(           -14       :2*n1  +16,-14       :2*n2  +16,-14:2*n3+16)

    integer ncount0,ncount1,ncount_rate,ncount_max,nflop
    integer l1,l2,i,t
    integer l1_0,l1_1,ll1

    real(kind=8) tel,y2i__0,y2i__1,y2i1_0,y2i1_1,y2i,y2i1
    include 'v_17.f90'

    open(unit=20,file='tree.flop')

    nflop=0
    do l1=-14+2*nfl1,2*nfu1+16
        do l2=-14+2*nfl2,2*nfu2+16
            if (ib(2,l1,l2).ge.ib(1,l1,l2)) nflop=nflop+(ib(2,l1,l2)-ib(1,l1,l2)+1)*31*2*2
        enddo
    enddo

    call system_clock(ncount0,ncount_rate,ncount_max)

    do l2=-14+2*nfl2,2*nfu2+16
        do l1=-14+2*nfl1,2*nfu1+16
            if (ib(1,l1,l2).le.ib(2,l1,l2)) then
                y(l1,l2,2*ib(2,l1,l2)+16)=y(l1,l2,2*ib(2,l1,l2)+16)+&
                fil2(16,1)*x(1,ib(2,l1,l2),l1,l2)+fil2(16,2)*x(2,ib(2,l1,l2),l1,l2)
            
                do i=ib(1,l1,l2)-7,ib(2,l1,l2)+7 
                    y2i=0.d0
                    y2i1=0.d0
                    do t=max(i-8,ib(1,l1,l2)),min(i+7,ib(2,l1,l2))
                        y2i=y2i+fil2(2*(i-t),1)*x(1,t,l1,l2)+fil2(2*(i-t),2)*x(2,t,l1,l2)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                    enddo
                    y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
                    y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
                enddo
            endif
        enddo
    enddo

    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)

    write(20,*) tel, 1.d-6*nflop/tel
end





    subroutine  comb_rot_grow_loc_square_1(n1,n2,n3,x,y,ib,ib2,loczero)
! In one dimesnion,
! with unoptimized cycles   
! Applies synthesis wavelet transformation 
! then convolves with magic filter
!  the size of the data is allowed to grow
    implicit none
    integer,intent(in)::n1,n2,n3
    real(kind=8),intent(in):: x(0:n1,0:n2,0:n3)
    integer,intent(in)::ib(2,0:n2,0:n3)
    integer,intent(in)::ib2(2,0:n3,-14:2*n1+16)
    logical,intent(in)::loczero
    
    real(kind=8),intent(out)::y(0:n2,0:n3,-14:2*n1+16)
    
    integer i,t,l2,l3,l1
    integer ll1,ll3,l10,l11,l30,l31,ll2,l21
    integer ncount0,ncount1,ncount2,ncount_rate,ncount_max,nflop
    real(kind=8) y2i,y2i1,tel
    real(kind=8) y2i__11, y2i__12, y2i1_11, y2i1_12, y2i__21, y2i__22, y2i1_21, y2i1_22
    real(kind=8) t0,t1
    
    include 'v_17.f90'

    open(unit=10,file='zero.square')

    nflop=0

    do l2=0,n2
        do l3=0,n3
            if (ib(2,l2,l3).ge.ib(1,l2,l3)) nflop=nflop+(ib(2,l2,l3)-ib(1,l2,l3)+1)*31*2
        enddo
    enddo
    call system_clock(ncount0,ncount_rate,ncount_max)

!   initialize the y array with zeroes
!   but only inside the region defined by ib2 array
!   which is the ib array for the next step

    if (loczero) call make_loczero(n1,n2,n3,ib2,y)
    
    call system_clock(ncount1,ncount_rate,ncount_max)
    
    t0=dble(ncount1-ncount0)/dble(ncount_rate)
    
!    the actual convolution
    
    do l3=0,n3
        do l2=0,n2
            
            if (ib(1,l2,l3).le.ib(2,l2,l3)) then
    
                y(l2,l3,2*ib(2,l2,l3)+16)=fil2(16,1)*x(ib(2,l2,l3),l2,l3)
            
                do i=ib(1,l2,l3)-7,ib(2,l2,l3)+7 
                    y2i=0.d0
                    y2i1=0.d0
                    do t=max(i-8,ib(1,l2,l3)),min(i+7,ib(2,l2,l3))
                        y2i=y2i+fil2(2*(i-t),1)*x(t  ,l2,l3)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(t  ,l2,l3)
                    enddo
                    y(l2,l3,2*i)=y2i
                    y(l2,l3,2*i+1)=y2i1
                enddo
        
            endif
        enddo
    enddo
    
    call system_clock(ncount2,ncount_rate,ncount_max)
    t1=dble(ncount2-ncount1)/dble(ncount_rate)
    tel=dble(ncount2-ncount0)/dble(ncount_rate)

    write(10,'(3f10.3,f10.0)') t0,t1,tel, 1.d-6*nflop/tel

end


subroutine make_loczero(n1,n2,n3,ib2,y)
!   initialize the y array with zeroes
!   but only inside the region defined by ib2 array
!   which is the ib array for the next step

    implicit none
    integer,intent(in)::n1,n2,n3
    integer,intent(in)::ib2(2,0:n3,-14:2*n1+16)
    real(kind=8),intent(out)::y(0:n2,0:n3,-14:2*n1+16)

    integer l1,l2,l3
    
        do l1=-14,2*n1+16
            do l3=0,n3
                do l2=ib2(1,l3,l1),ib2(2,l3,l1)
                    y(l2,l3,l1)=0.d0
                enddo
            enddo
        enddo
end subroutine make_loczero
