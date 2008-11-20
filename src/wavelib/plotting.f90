subroutine plot_wf(orbname,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,rx,ry,rz,wfd,&
     bounds,psi)
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(in) :: bounds
  character(len=10) :: orbname 
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid,rx,ry,rz
  real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
  !local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: nw1,nw2,i_stat,i_all,i
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:), allocatable :: psir,w1,w2
  real(wp), dimension(:,:,:), allocatable :: x_c
  real(wp), dimension(:,:,:,:), allocatable :: x_f

  nw1=max((n3+1)*(2*n1+31)*(2*n2+31),&
       (n1+1)*(2*n2+31)*(2*n3+31),&
       2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

  nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
       (n1+1)*(n2+1)*(2*n3+31))

  do i=0,3
     scal(i)=1.d0
  enddo

  allocate(x_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,x_c,'x_c',subname)
  allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f,'x_f',subname)
  allocate(w1(nw1+ndebug),stat=i_stat)
  call memocc(i_stat,w1,'w1',subname)
  allocate(w2(nw2+ndebug),stat=i_stat)
  call memocc(i_stat,w2,'w2',subname)
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  call razero((n1+1)*(n2+1)*(n3+1),x_c)
  call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f)

  call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
       wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
       wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
       scal,psi(1),psi(wfd%nvctr_c+1),x_c,x_f)

  call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
       psir,bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,bounds%gb%ibyz_ff,&
       bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,bounds%ibyyzz_r)

  !call plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,psir)
  call plot_pot_full(rx,ry,rz,hgrid,n1,n2,n3,orbname,psir)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(x_c))*kind(x_c)
  deallocate(x_c,stat=i_stat)
  call memocc(i_stat,i_all,'x_c',subname)

  i_all=-product(shape(x_f))*kind(x_f)
  deallocate(x_f,stat=i_stat)
  call memocc(i_stat,i_all,'x_f',subname)

  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1',subname)

  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2',subname)

END SUBROUTINE plot_wf

subroutine plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,pot)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,iounit
  real(gp), intent(in) :: rx,ry,rz,hgrid
  real(dp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3
  real(gp) :: hgridh

  hgridh=.5d0*hgrid
  open(iounit) 
  open(iounit+1) 
  open(iounit+2) 

  i3=nint(rz/hgridh)
  i2=nint(ry/hgridh)
  write(*,*) 'plot_p, i2,i3,n2,n3 ',i2,i3,n2,n3
  do i1=-14,2*n1+16
     write(iounit,*) real(i1,gp)*hgridh,pot(i1,i2,i3)
  enddo

  i1=nint(rx/hgridh)
  i2=nint(ry/hgridh)
  write(*,*) 'plot_p, i1,i2 ',i1,i2
  do i3=-14,2*n3+16
     write(iounit+1,*) real(i3,gp)*hgridh,pot(i1,i2,i3)
  enddo

  i1=nint(rx/hgridh)
  i3=nint(rz/hgridh)
  write(*,*) 'plot_p, i1,i3 ',i1,i3
  do i2=-14,2*n2+16
     write(iounit+2,*) real(i2,gp)*hgridh,pot(i1,i2,i3)
  enddo

  close(iounit) 
  close(iounit+1) 
  close(iounit+2) 

end subroutine plot_pot

subroutine plot_pot_full(rx,ry,rz,hgrid,n1,n2,n3,orbname,pot)
  use module_base
  implicit none
  character(len=10), intent(in) :: orbname
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: rx,ry,rz,hgrid
  real(dp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3
  real(gp) :: hgridh

!virtual orbitals are identified by their name
!  write(*,'(1x,a,i0)')'printing orbital number= ',iounit
  write(*,'(A)')'printing '//orbname
  hgridh=.5d0*hgrid
!  write(orbname,'(i0)')iounit
!  open(unit=22,file='psi'//orbname//'.pot',status='unknown')
  open(unit=22,file=orbname//'.pot',status='unknown')
!  write(22,*)'orbital'//orbname
  write(22,*)orbname
  write(22,*) 2*n1+1,2*n2+1,2*n3+1
  write(22,*) n1*hgrid,' 0. ',n2*hgrid
  write(22,*) ' 0. ',' 0. ',n3*hgrid
  write(22,*)'xyz'
!   the density can easily obtained later, if needed.
!!  open(unit=23,file='den'//orbname//'.pot',status='unknown')
!  open(unit=23,file=orbname//'.dens.pot',status='unknown')
!  write(23,*)orbname
!  write(23,*) 2*n1+1,2*n2+1,2*n3+1
!  write(23,*) n1*hgrid,' 0. ',n2*hgrid
!  write(23,*) ' 0. ',' 0. ',n3*hgrid
!  write(23,*)'xyz'
  do i3=0,2*n3
     do i2=0,2*n2
        do i1=0,2*n1
           write(22,*)pot(i1,i2,i3)
           !write(23,*)pot(i1,i2,i3)**2
        end do
     end do
  end do
  close(unit=22) 
!  close(unit=23) 

end subroutine plot_pot_full

subroutine plot_psifscf(iunit,hgrid,n1,n2,n3,psifscf)
  use module_base
  implicit none
  integer, intent(in) :: iunit,n1,n2,n3
  real(gp), intent(in) :: hgrid
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(in) :: psifscf
  !local variables
  integer :: i1,i2,i3,i
  real(gp) :: hgridh

  hgridh=.5d0*hgrid

  ! along x-axis
  i3=n3
  i2=n2
  do i1=-7,2*n1+8
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
        real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! 111 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=i ; i2=i ; i3=i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! 1-1-1 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=i ; i2=-i ; i3=-i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! -11-1 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=-i ; i2=i ; i3=-i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! -1-11 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=-i ; i2=-i ; i3=i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

end subroutine plot_psifscf

subroutine plot_density(geocode,filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nelec,&
     alat1,alat2,alat3,ngatherarr,rho)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nelec,nproc
  real(gp), intent(in) :: alat1,alat2,alat3
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(dp), dimension(n1i*n2i*n3p), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density'
  integer :: nl1,nl2,nl3,i_all,i_stat,i1,i2,i3,ind,ierr
  real(dp) :: later_avg
  real(dp), dimension(:), pointer :: pot_ion

  if (iproc == 0) then
     open(unit=22,file=filename,status='unknown')
     write(22,*)'normalised density'
     write(22,*) 2*n1+2,2*n2+2,2*n3+2
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     write(22,*)'xyz   periodic' !not true in general but needed in the case
  end if

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if

  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(pot_ion(n1i*n2i*n3i+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)

     call MPI_ALLGATHERV(rho,n1i*n2i*n3p,&
          mpidtypd,pot_ion,ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
  else
     pot_ion => rho
  end if

  if (iproc == 0) then
     do i3=0,2*n3+1
        do i2=0,2*n2+1
           do i1=0,2*n1+1
              ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              write(22,*)pot_ion(ind)/real(nelec,dp)
           end do
        end do
     end do
     close(22)
     !laterally averaged potential following the isolated direction in surfaces case
     !if (geocode == 'S') then
        open(unit=23,file=filename//'_avg'//geocode,status='unknown')
        later_avg=0.0_dp
        do i2=0,2*n2+1
           do i3=0,2*n3+1
              do i1=0,2*n1+1
                 ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                 later_avg=later_avg+pot_ion(ind)
              end do
           end do
           later_avg=later_avg/real((2*n1+2)*(2*n3+2),dp) !2D integration/2D Volume
           write(23,*)i2,alat2/real(2*n2+2,dp)*i2,later_avg
        end do
        close(23)
     !end if
  end if

  if (nproc > 1) then
       i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion',subname)
  end if

end subroutine plot_density

subroutine plot_density_cube(filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
     hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nspin,nproc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(n1i*n2i*n3p,nspin), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density_cube'
  character(len=3) :: advancestring
  integer :: nl1,nl2,nl3,i_all,i_stat,i1,i2,i3,ind,ierr,icount,j,iat
  real(dp), dimension(:,:), pointer :: pot_ion

  if (iproc == 0) then
     open(unit=22,file=filename,status='unknown')
     write(22,*)'CUBE file for charge density'
     write(22,*)'Case for total spin'
     !number of atoms
     write(22,'(i5,3(f12.6))') at%nat,0.0_gp,0.0_gp,0.0_gp
     !grid and grid spacings
     write(22,'(i5,3(f12.6))') 2*n1+2,hxh,0.0_gp,0.0_gp
     write(22,'(i5,3(f12.6))') 2*n2+2,0.0_gp,hyh,0.0_gp
     write(22,'(i5,3(f12.6))') 2*n3+2,0.0_gp,0.0_gp,hzh
     !atomic number and positions
     do iat=1,at%nat
        write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
     end do
  end if

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (at%geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (at%geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if

  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(pot_ion(n1i*n2i*n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)

     call MPI_ALLGATHERV(rho,n1i*n2i*n3p,&
          mpidtypd,pot_ion,ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
  else
     pot_ion => rho
  end if

  if (iproc == 0) then
     !the loop is reverted for a cube file
     !charge normalised to the total charge
     icount=0
     do i1=0,2*n1+1
        do i2=0,2*n2+1
           do i3=0,2*n3+1
              icount=icount+1
              if (icount == 6) then
                 advancestring='yes'
                 icount=0
              else
                 advancestring='no'
              end if

              ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              write(22,'(1x,1pe12.6)',advance=advancestring)pot_ion(ind,1)
           end do
        end do
     end do
     close(22)
  end if

  if (nproc > 1) then
       i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion',subname)
  end if

end subroutine plot_density_cube
