



subroutine synthese_per_old(nd1,nd2,nd3,x,y,ww)
  ! a periodic synthesis (backward) wavelet transformation
  ! the input array x is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  syn_rot_per_old(nd1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  syn_rot_per_old(nd2,nt,y,ww)
  ! i3,i1,i2  -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  syn_rot_per_old(nd3,nt,ww,y)

end subroutine synthese_per_old

subroutine synthese_per_old_self(nd1,nd2,nd3,x,y,ww)
  ! a periodic synthesis (backward) wavelet transformation
  ! the input array x is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  syn_rot_per_old(nd1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  syn_rot_per_old(nd2,nt,y,ww)
  ! i3,i1,i2  -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  syn_rot_per_old(nd3,nt,ww,x)

end subroutine synthese_per_old_self


subroutine analyse_per_old(nd1,nd2,nd3,y,x,ww)
  ! an analysis (forward) periodic wavelet transformation
  ! the input array y is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: y
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  ana_rot_per_old(nd1,nt,y,x)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  ana_rot_per_old(nd2,nt,x,ww)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  ana_rot_per_old(nd3,nt,ww,x)

end subroutine analyse_per_old

subroutine analyse_per_old_self(nd1,nd2,nd3,y,x,ww)
  ! an analysis (forward) periodic wavelet transformation
  ! the input array y is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: y
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  ana_rot_per_old(nd1,nt,y,x)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  ana_rot_per_old(nd2,nt,x,ww)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  ana_rot_per_old(nd3,nt,ww,y)

end subroutine analyse_per_old_self

subroutine syn_repeated_per(nd1,nd2,nd3,x,num_trans,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  character(len=*), parameter :: subname='syn_repeated_per'
  integer :: nn1,nn2,nn3,i_trans,i_all,i_stat,i1,i2,i3,i
  real(wp), dimension(:), allocatable :: xx,yy,ww

  if (num_trans >= 1)  then

     allocate(yy((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,yy,'yy',subname)

     allocate(xx((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,xx,'xx',subname)

  endif

  if (num_trans >= 2) then

     nn1=(nd1+1)/2-1
     nn2=(nd2+1)/2-1
     nn3=(nd3+1)/2-1

     allocate(ww((nn1+1)*(nn2+1)*(nn3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     do i_trans=1,num_trans-1

        n1=2*(n1+1)-1
        n2=2*(n2+1)-1
        n3=2*(n3+1)-1

        if (n1.gt.nd1) stop 'n1 beyond borders'
        if (n2.gt.nd2) stop 'n2 beyond borders'
        if (n3.gt.nd3) stop 'n3 beyond borders'

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 xx(i)=x(i1,i2,i3)
                 i=i+1
              enddo
           enddo
        enddo

        call synthese_per_old(n1,n2,n3,xx,yy,ww)

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 x(i1,i2,i3)=yy(i)
                 i=i+1
              enddo
           enddo
        enddo

     enddo

     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww',subname)

  endif

  if (num_trans >= 1) then

     n1=2*(n1+1)-1
     n2=2*(n2+1)-1
     n3=2*(n3+1)-1

     call synthese_per_old_self(n1,n2,n3,x,xx,yy)

     i_all=-product(shape(xx))*kind(xx)
     deallocate(xx,stat=i_stat)
     call memocc(i_stat,i_all,'xx',subname)

     i_all=-product(shape(yy))*kind(yy)
     deallocate(yy,stat=i_stat)
     call memocc(i_stat,i_all,'yy',subname)

  endif

end subroutine syn_repeated_per



subroutine ana_repeated_per(nd1,nd2,nd3,x,num_trans,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  character(len=*), parameter :: subname='ana_repeated_per'
  integer :: nn1,nn2,nn3,i_trans,i_all,i_stat,i1,i2,i3,i
  real(wp), dimension(:), allocatable :: xx,yy,ww

  n1=nd1
  n2=nd2
  n3=nd3

  if (num_trans.ge.1)  then

     allocate(yy((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,yy,'yy',subname)

     allocate(xx((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,xx,'xx',subname)

     call analyse_per_old_self(n1,n2,n3,x,yy,xx)

     n1=(n1+1)/2-1
     n2=(n2+1)/2-1
     n3=(n3+1)/2-1

  endif

  if (num_trans.ge.2) then

     allocate(ww((n1+1)*(n2+1)*(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     do i_trans=2,num_trans

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 xx(i)=x(i1,i2,i3)
                 i=i+1
              enddo
           enddo
        enddo

        call analyse_per_old(n1,n2,n3,xx,yy,ww)

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 x(i1,i2,i3)=yy(i)
                 i=i+1
              enddo
           enddo
        enddo

        n1=(n1+1)/2-1
        n2=(n2+1)/2-1
        n3=(n3+1)/2-1

     enddo

     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww',subname)

  endif

  if (num_trans.ge.1) then 
     i_all=-product(shape(xx))*kind(xx)
     deallocate(xx,stat=i_stat)
     call memocc(i_stat,i_all,'xx',subname)

     i_all=-product(shape(yy))*kind(yy)
     deallocate(yy,stat=i_stat)
     call memocc(i_stat,i_all,'yy',subname)
  endif

end subroutine ana_repeated_per