


subroutine make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
    ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
        ibyz_c,ibzxx_c,ibxxyy_c,ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
!    creates complicated ib arrays    
    implicit none
    integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
    integer i1,i2,i3,nt,m1,m2,m3

    integer,intent(in):: ibyz_c(2,0:n2,0:n3),ibxy_c(2,0:n1,0:n2)
    integer,intent(in):: ibyz_f(2,0:n2,0:n3),ibxy_f(2,0:n1,0:n2)

!    for shrink:    
    integer,intent(out):: ibzzx_c(2,-14:2*n3+16,0:n1) 
    integer,intent(out):: ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

    integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
    integer,intent(out):: ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
    integer,intent(out):: ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

!    for grow:    
    integer,intent(out):: ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
    integer,intent(out):: ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

    integer,intent(out):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
    integer,intent(out):: ibzxx_f(2,          nfl3:nfu3,2*nfl1-14:2*nfu1+16)
    integer,intent(out):: ibxxyy_f(2,                    2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

    integer,allocatable,dimension(:,:,:)::ibyx_c,ibxzz_c,ibzzyy_c
    integer,allocatable,dimension(:,:,:)::ibyx_f,ibxzz_f,ibzzyy_f

    logical,allocatable:: logrid_big(:)
    
!    for real space:
    integer,intent(out):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

    m1=nfu1-nfl1; m2=nfu2-nfl2; m3=nfu3-nfl3

    
!   (0:n3,-14:2*n1+16,-14:2*n2+16) from grow
!   (-14:2*n2+16,-14:2*n3+16,0:n1) from shrink
!   (-14:2*n1+16,-14:2*n2+16,-14:2*n3+16) from real space

    allocate(logrid_big((2*n1+31)*(2*n2+31)*(2*n3+31)))

    allocate(ibyx_c(2,0:n2,0:n1),ibxzz_c(2,-14:2*n3+16,0:n1))
    allocate(ibzzyy_c(2,-14:2*n2+16,-14:2*n3+16))

    allocate(ibyx_f(2,nfl2:nfu2,nfl1:nfu1),ibxzz_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1))
    allocate(ibzzyy_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16))


!    transpose the original array    
    call ib_transpose(ibxy_c,ibyx_c,n1,n2)

    do i2=nfl2,nfu2
        do i1=nfl1,nfu1
            ibyx_f (:,i2,i1)=ibxy_f(:,i1,i2)
            ibxy_ff(:,i1,i2)=ibxy_f(:,i1,i2)
        enddo
    enddo

!    create transposed arrays (dimensions 1 and 3 switched)
    call make_ib(logrid_big,ibyx_c,ibxzz_c,ibzzyy_c,0,n3,0,n2,0,n1)
    call make_ib(logrid_big,ibyx_f,ibxzz_f,ibzzyy_f,nfl3,nfu3,nfl2,nfu2,nfl1,nfu1)
    
!    transpose them back
    call ib_transpose(ibxzz_c,ibzzx_c,n1,2*n3+30)
    call ib_transpose(ibzzyy_c,ibyyzz_c,2*n3+30,2*n2+30)

    call ib_transpose(ibxzz_f,ibzzx_f,m1,2*m3+30)
    call ib_transpose(ibzzyy_f,ibyyzz_f,2*m3+30,2*m2+30)

!    for grow:

    do i2=nfl2,nfu2
        do i3=nfl3,nfu3
            ibyz_ff(:,i2,i3)=ibyz_f(:,i2,i3)
        enddo
    enddo

    call make_ib_c(logrid_big,ibyz_c,ibzxx_c,ibxxyy_c,n1,n2,n3)
    
    call make_ib(logrid_big,ibyz_ff,ibzxx_f,ibxxyy_f,&
    nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
    call squares_1d(ibxxyy_f,2*nfl1-14,2*nfu1+16,2*nfl2-14,2*nfu2+16)
    
!    make real space borders
    nt=(2*n1+31)*(2*n2+31)
    call ib_to_logrid_rot(ibxxyy_c,logrid_big,0,n3,nt)
    
    nt=(2*n2+31)*(2*n3+31)        
!    call ib_from_logrid(ibyyzz_r,logrid_big,-14,2*n1+16,nt)
    call ib_from_logrid(ibyyzz_r,logrid_big,0,2*n1+30,nt)
            
    deallocate(logrid_big,ibyx_c,ibzzyy_c,ibyx_f,ibzzyy_f)

    
end subroutine make_all_ib

subroutine make_ib_c(logrid_big,ibyz,ibzxx,ibxxyy,n1,n2,n3)
!    This subroutine mimics the comb_grow_f one
    implicit none
    integer nt,n1,n2,n3
    integer ibyz(  2,0:n2,0:n3)! input
    integer ibzxx( 2,          0:n3,-14:2*n1+16)!output
    integer ibxxyy(2,                    -14:2*n1+16,-14:2*n2+16)!output
    logical logrid_big(           0:n3,-14:2*n1+16,-14:2*n2+16)! work array

        call squares(ibyz,n2,n3)
! i1,i2,i3 -> i2,i3,I1
        nt=(n2+1)*(n3+1)
        call ib_to_logrid_rot(  ibyz,logrid_big,0,n1,nt)

! i2,i3,I1 -> i3,I1,I2
        nt=(n3+1)*(2*n1+31)
        call ib_from_logrid( ibzxx,logrid_big,0,n2,nt)
        call squares(ibzxx,n3,2*n1+30)
        
        call ib_to_logrid_rot( ibzxx,logrid_big,0,n2,nt)

! i3,I1,I2  -> I1,I2,I3
        nt=(2*n1+31)*(2*n2+31)
        call ib_from_logrid(ibxxyy,logrid_big,0,n3,nt)
        call squares(ibxxyy,2*n1+30,2*n2+30)
    
end subroutine make_ib_c


subroutine make_ib(logrid_big,ibyz,ibzxx,ibxxyy,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!    This subroutine mimics the comb_grow_f one
    implicit none
    integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
    integer ibyz(  2,nfl2:nfu2,nfl3:nfu3)! input
    integer ibzxx( 2,          nfl3:nfu3,2*nfl1-14:2*nfu1+16)!output
    integer ibxxyy(2,                    2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)!output
    logical logrid_big(           nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array

! i1,i2,i3 -> i2,i3,I1
        nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
        call ib_to_logrid_rot(  ibyz,logrid_big,nfl1,nfu1,nt)

! i2,i3,I1 -> i3,I1,I2
        nt=(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)
        call ib_from_logrid( ibzxx,logrid_big,nfl2,nfu2,nt)
        call ib_to_logrid_rot( ibzxx,logrid_big,nfl2,nfu2,nt)

! i3,I1,I2  -> I1,I2,I3
        nt=(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31)
        call ib_from_logrid(ibxxyy,logrid_big,nfl3,nfu3,nt)
    
end subroutine make_ib




subroutine ib_transpose(ib,ib_t,n1,n2)
    implicit none
    integer n1,n2
    integer i1,i2
    integer ib(2,0:n1,0:n2),ib_t(2,0:n2,0:n1)

    do i1=0,n1
        do i2=0,n2
            ib_t(:,i2,i1)=ib(:,i1,i2)
        enddo
    enddo

end subroutine ib_transpose



subroutine make_logrid(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,logrid)
    implicit none
    integer n1,n2,n3
    integer i1,i2,i3 ! counters
    integer j1,j2,j3 ! positions relative to box center
    integer nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
    integer nm1,nm2,nm3 ! medium point of the box
    integer nd1,nd2,nd3 ! dimension of the box
    integer nrad,nrad2  ! internal sphere radius and radius squared
    logical logrid(0:n1,0:n2,0:n3)

    logrid=.false.

    nm1=(nfl1+nfu1)/2
    nm2=(nfl2+nfu2)/2
    nm3=(nfl3+nfu3)/2

    nd1=nm1-nfl1
    nd2=nm2-nfl2
    nd3=nm3-nfl3

    nrad=min(nd1,nd2,nd3)
    nrad2=nrad*nrad

    do i1=nfl1,nfu1        
        do i2=nfl2,nfu2
            do i3=nfl3,nfu3
                j1=i1-nm1
                j2=i2-nm2
                j3=i3-nm3
                if (j1*j1+j2*j2+j3*j3.le.nrad2) logrid(i1,i2,i3)=.true.
            enddo
        enddo
    enddo
    
end subroutine make_logrid



subroutine ib_to_logrid_rot(ib,logrid,nfl,nfu,ndat)
! This one mimics the comb_rot_grow_f_loc
    implicit none
    integer ndat,nfl,nfu,l,i
    integer ib(2,ndat)! input
    logical logrid(ndat,-14+2*nfl:2*nfu+16)! output

    logrid=.false.

    do l=1,ndat
        do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
            logrid(l,i)=.true.
        enddo
    enddo

end subroutine ib_to_logrid_rot



subroutine ibyz_to_logrid(n1,n2,n3,logrid,ibyz)
    implicit none
    integer n1,n2,n3
    logical logrid(0:n1,0:n2,0:n3)
    integer ibyz(2,0:n2,0:n3)

    integer i1,i2,i3
    
    logrid=.false.
    do i3=0,n3
        do i2=0,n2
            do i1=ibyz(1,i2,i3),ibyz(2,i2,i3)
                logrid(i1,i2,i3)=.true.
            enddo
        enddo
    enddo

end subroutine ibyz_to_logrid


subroutine ib_from_logrid(ib,logrid,ml1,mu1,ndat)
    ! mimics the bounds subroutine    
    implicit none
    integer i,i1
    integer ml1,mu1,ndat
    integer ib(2,ndat)
    logical logrid(ml1:mu1,ndat)

    do i=1,ndat
        ib(1,i)= 1000
        ib(2,i)=-1000
        
        inner1:do i1=ml1,mu1
            if (logrid(i1,i)) then 
                ib(1,i)=i1
                exit inner1
            endif
        enddo inner1
        
        inner2:do i1=mu1,ml1,-1
            if (logrid(i1,i)) then 
                ib(2,i)=i1
                exit inner2
            endif
        enddo inner2
    enddo

end subroutine ib_from_logrid


!subroutine bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
!    implicit real(kind=8) (a-h,o-z)
!    logical logrid
!    dimension logrid(0:n1,0:n2,0:n3)
!    dimension ibyz(2,0:n2,0:n3),ibxz(2,0:n1,0:n3),ibxy(2,0:n1,0:n2)
!    
!    
!    do i3=0,n3 
!        do i2=0,n2 
!            ibyz(1,i2,i3)= 1000
!            ibyz(2,i2,i3)=-1000
!            
!            inner1:do i1=0,n1
!             if (logrid(i1,i2,i3)) then 
!                ibyz(1,i2,i3)=i1
!                exit inner1
!             endif
!            enddo inner1
!            
!            inner2:do i1=n1,0,-1
!             if (logrid(i1,i2,i3)) then 
!                ibyz(2,i2,i3)=i1
!                exit inner2
!             endif
!            enddo inner2
!        enddo
!    enddo
!
!
!    do i3=0,n3 
!        do i1=0,n1
!            ibxz(1,i1,i3)= 1000
!            ibxz(2,i1,i3)=-1000
!            
!            inner3:do i2=0,n2 
!             if (logrid(i1,i2,i3)) then 
!                 ibxz(1,i1,i3)=i2
!                 exit inner3
!             endif
!            enddo inner3
!
!            inner4:do i2=n2,0,-1
!             if (logrid(i1,i2,i3)) then 
!                 ibxz(2,i1,i3)=i2
!                 exit inner4
!             endif
!            enddo inner4
!
!        enddo    
!    enddo
!
!
!    do i2=0,n2 
!        do i1=0,n1 
!            ibxy(1,i1,i2)= 1000
!            ibxy(2,i1,i2)=-1000
!            
!            inner5: do i3=0,n3
!             if (logrid(i1,i2,i3)) then 
!                 ibxy(1,i1,i2)=i3
!                 exit inner5
!             endif
!            enddo inner5
!        
!            inner6:do i3=n3,0,-1
!             if (logrid(i1,i2,i3)) then 
!                 ibxy(2,i1,i2)=i3
!                 exit inner6
!             endif
!            enddo inner6
!        
!        enddo
!    enddo
!
!    return
!END SUBROUTINE



subroutine check_ibyz(n1,n2,n3,logrid,logrid1,ibyz)
    implicit none
    integer n1,n2,n3
    logical,dimension(0:n1,0:n2,0:n3)::logrid,logrid1
    integer ibyz(2,0:n2,0:n3)

    integer i1,i2,i3
    logical flag

    call  ibyz_to_logrid(n1,n2,n3,logrid1,ibyz)

    flag=.true.
    do i1=0,n1
        do i2=0,n2
            do i3=0,n3
                flag=(flag.and.(logrid(i1,i2,i3) .eqv. logrid1(i1,i2,i3)))
            enddo
        enddo
    enddo
    write(*,*)'flag=',flag
    stop

end subroutine check_ibyz

subroutine squares_1d(ib,nfl2,nfu2,nfl3,nfu3)
	! modifies the ib array 
	! so that it is made up of blocks of size 2
	! the localization region is enlarged as a result

	! works for even nfl2 only
	implicit none
	integer,intent(in)::nfl2,nfu2,nfl3,nfu3
   	integer,intent(inout)::ib(2,nfl2:nfu2,nfl3:nfu3)

	integer i2,i3,ii2,ii3,ibmin,ibmax

	 do i3=nfl3,nfu3
		 do i2=nfl2/2,(nfu2-1)/2
			ii2=2*i2
			ibmin=min(ib(1,ii2,i3),ib(1,ii2+1,i3))

			ib(1,ii2,i3)=ibmin
			ib(1,ii2+1,i3)=ibmin

			ibmax=max(ib(2,ii2,i3),ib(2,ii2+1,i3))

			ib(2,ii2,i3)=ibmax
			ib(2,ii2+1,i3)=ibmax
		 enddo
	 enddo
end subroutine squares_1d


subroutine squares(ib,n2,n3)
	! modifies the ib array 
	! so that it is made up of squares 2x2
	! the localization region is enlarged as a result
	implicit none
	integer,intent(in)::n2,n3
   	integer,intent(inout)::ib(2,0:n2,0:n3)

	integer i2,i3,ii2,ii3,ibmin,ibmax

	 do i3=0,(n3-1)/2
		 do i2=0,(n2-1)/2
			ii3=2*i3
			ii2=2*i2
			ibmin=min(ib(1,ii2,ii3),ib(1,ii2+1,ii3),&
			ib(1,ii2,ii3+1),ib(1,ii2+1,ii3+1))

			ib(1,ii2,ii3)=ibmin
			ib(1,ii2+1,ii3)=ibmin
			ib(1,ii2,ii3+1)=ibmin
			ib(1,ii2+1,ii3+1)=ibmin

			ibmax=max(ib(2,ii2,ii3),ib(2,ii2+1,ii3),&
			ib(2,ii2,ii3+1),ib(2,ii2+1,ii3+1))

			ib(2,ii2,ii3)=ibmax
			ib(2,ii2+1,ii3)=ibmax
			ib(2,ii2,ii3+1)=ibmax
			ib(2,ii2+1,ii3+1)=ibmax
		 enddo
	 enddo
end subroutine squares
