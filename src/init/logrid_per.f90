subroutine make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f)
  !    creates complicated ib arrays    
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer i1,i2,i3,nt,m1,m2,m3,i_stat,i_all

  integer,intent(in):: ibyz_f(2,0:n2,0:n3+ndebug),ibxy_f(2,0:n1,0:n2+ndebug)

  !    for shrink:    
  integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer,intent(out):: ibzzx_f(2,0:2*n3+1,nfl1:nfu1) 
  integer,intent(out):: ibyyzz_f(2,0:2*n2+1,0:2*n3+1)

  !    for grow:    
  integer,intent(out):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer,intent(out):: ibzxx_f(2,nfl3:nfu3,0:2*n1+1)
  integer,intent(out):: ibxxyy_f(2,0:2*n1+1,0:2*n2+1)

  character(len=*), parameter :: subname=' make_all_ib'
  logical,allocatable:: logrid_big(:)

  allocate(logrid_big((2*n1+2)*(2*n2+2)*(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_big,'logrid_big',subname)

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  !   (0:n3,0:2*n1+1,0:2*n2+1) from grow
  !   (0:2*n2+1,0:2*n3+1,0:n1) from shrink

  !	for shrink:
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        ibxy_ff(:,i1,i2)=ibxy_f(:,i1,i2)
     enddo
  enddo

  call make_ib_inv_per(logrid_big,ibxy_ff,ibzzx_f,ibyyzz_f,&
             n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

  !    for grow:
  do i2=nfl2,nfu2
     do i3=nfl3,nfu3
        ibyz_ff(:,i2,i3)=ibyz_f(:,i2,i3)
!		write(9,*) i2,i3,ibyz_ff(1,i2,i3),ibyz_ff(2,i2,i3)
     enddo
  enddo

  call make_ib_per(logrid_big,ibyz_ff,ibzxx_f,ibxxyy_f,n1,n2,n3,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

  i_all=-product(shape(logrid_big))*kind(logrid_big)
  deallocate(logrid_big,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_big',subname)

end subroutine make_all_ib_per


subroutine make_ib_inv_per(logrid_big,ibxy,ibzzx,ibyyzz,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  !    This subroutine mimics the comb_grow_f one
  implicit none
  integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1,n2,n3
  integer,intent(out):: ibxy(2,nfl1:nfu1,nfl2:nfu2)
  integer,intent(out):: ibzzx(2,0:2*nfu3+1,nfl1:nfu1) 
  integer,intent(out):: ibyyzz(2,0:2*n2+1,0:2*n3+1)
  logical logrid_big(nfl3:nfu3,0:2*n1+1,0:2*n2+1)! work array

  ! I3,i1,i2 -> i1,i2,i3 
  nt=(nfu1-nfl1+1)*(nfu2-nfl2+1)
  call ib_to_logrid_inv_per(ibxy,logrid_big,n3,nt)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+2)*(nfu1-nfl1+1)
  call ib_from_logrid_inv(ibzzx,logrid_big,nfl2,nfu2,nt)
  call ib_to_logrid_inv_per(ibzzx,logrid_big,n2,nt)

  ! I1,I2,I3  -> I2,I3,i1
  nt=(2*n2+2)*(2*n3+2)
  call ib_from_logrid_inv( ibyyzz,logrid_big,nfl1,nfu1,nt)

end subroutine make_ib_inv_per


subroutine ib_to_logrid_inv_per(ib,logrid,n,ndat)
  ! This one mimics the comb_rot_grow_f_loc
  implicit none
  integer ndat,l,i,n,ii
  integer ib(2,ndat)! input
  logical logrid(0:2*n+1,ndat)! output

  logrid=.false.

  do l=1,ndat
     do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
		ii=modulo(i,2*n+2) 
        logrid(ii,l)=.true.
     enddo
  enddo

end subroutine ib_to_logrid_inv_per


subroutine make_ib_per(logrid_big,ibyz,ibzxx,ibxxyy,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  !    This subroutine mimics the comb_grow_f one
  implicit none
  integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1,n2,n3
  integer ibyz(  2,nfl2:nfu2,nfl3:nfu3)! input
  integer ibzxx( 2,          nfl3:nfu3,0:2*n1+1)!output
  integer ibxxyy(2,                    0:2*n1+1,0:2*n2+1)!output
  logical logrid_big(        nfl3:nfu3,0:2*n1+1,0:2*n2+1)! work array

  ! i1,i2,i3 -> i2,i3,I1
  nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call ib_to_logrid_rot_per(  ibyz,logrid_big,n1,nt)

  ! i2,i3,I1 -> i3,I1,I2
  nt=(nfu3-nfl3+1)*(2*n1+2)
  call ib_from_logrid( ibzxx,logrid_big,nfl2,nfu2,nt)
  call ib_to_logrid_rot_per( ibzxx,logrid_big,n2,nt)

  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+2)*(2*n2+2)
  call ib_from_logrid(ibxxyy,logrid_big,nfl3,nfu3,nt)

end subroutine make_ib_per

subroutine ib_to_logrid_rot_per(ib,logrid,n,ndat)
  ! This one mimics the comb_rot_grow_f_loc
  implicit none
  integer ndat,n,l,i,ii
  integer ib(2,ndat)! input
  logical logrid(ndat,0:2*n+1)! output

  logrid=.false.

  do l=1,ndat
     do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
		ii=modulo(i,2*n+2) 
        logrid(l,i)=.true.
     enddo
  enddo

end subroutine ib_to_logrid_rot_per

subroutine make_logrid_f(n1,n2,n3, & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     logrid)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f

  logical,intent(out),dimension(0:n1,0:n2,0:n3)::logrid
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  logrid=.false.
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
		logrid(i,i2,i3)=.true. 
     enddo
  enddo

END SUBROUTINE make_logrid_f
