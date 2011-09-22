!> @file
!!  Routines to plot wavefunctions
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Routine to plot wavefunctions (old version)
subroutine plot_wf_old(kindplot,orbname,nexpo,at,lr,hx,hy,hz,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
  character(len=*) :: kindplot
  character(len=10) :: comment
  character(len=11) :: orbname
  integer, intent(in) :: nexpo
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
  !local variables
  character(len=*), parameter :: subname='plot_wf_old'
  integer :: i_stat,i_all
  integer :: nl1,nl2,nl3,n1i,n2i,n3i,n1,n2,n3
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: psir

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i

  call initialize_work_arrays_sumrho(lr,w)
 
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if
 
  call daub_to_isf(lr,w,psi,psir)

  if (lr%geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (lr%geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if

  if (trim(kindplot)=='POT') then
     !call plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,psir)
     call plot_pot_full(nexpo,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
          nl1,nl2,nl3,orbname,psir,comment)
  else if (trim(kindplot)=='CUBE') then
     call plot_cube_full(nexpo,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
          nl1,nl2,nl3,orbname,psir,comment)
  else
     stop 'ERROR: plotting format not recognized'
  end if
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf_old


subroutine plot_wf_cube(orbname,at,lr,hx,hy,hz,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
!HU  character(len=10) :: orbname,comment 
  character(len=10) :: comment
  character(len=11) :: orbname
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
  !local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: nw1,nw2,i_stat,i_all,i,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  integer :: nxc,nxf,nl1,nl2,nl3,n1i,n2i,n3i
  real(gp) :: rx, ry, rz
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:), allocatable :: psir,w1,w2,x_c_psifscf,x_f_psig

  rx=rxyz(1,1)
  ry=rxyz(2,1)
  rz=rxyz(3,1)

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3

  do i=0,3
     scal(i)=1.d0
  enddo

  select case(lr%geocode)
  case('F')
     !dimension of the work arrays
     ! shrink convention: nw1>nw2
     nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))
     nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
  case('S')
     !dimension of the work arrays
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
     nxf=1
  case('P')
     !dimension of the work arrays, fully periodic case
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=1
  end select
  !work arrays
  allocate(x_c_psifscf(nxc+ndebug),stat=i_stat)
  call memocc(i_stat,x_c_psifscf,'x_c_psifscf',subname)
  allocate(x_f_psig(nxf+ndebug),stat=i_stat)
  call memocc(i_stat,x_f_psig,'x_f_psig',subname)
  allocate(w1(nw1+ndebug),stat=i_stat)
  call memocc(i_stat,w1,'w1',subname)
  allocate(w2(nw2+ndebug),stat=i_stat)
  call memocc(i_stat,w2,'w2',subname)

  allocate(psir(n1i*n2i*n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(nxc,x_c_psifscf)
     call razero(nxf,x_f_psig)
     call razero(n1i*n2i*n3i,psir)
  end if
  select case(lr%geocode)
  case('F')
     call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),  & 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1), &
          scal,psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,x_f_psig)

     call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,&
          x_c_psifscf,x_f_psig,  & 
          psir(1),lr%bounds%kb%ibyz_c,lr%bounds%gb%ibzxx_c,lr%bounds%gb%ibxxyy_c,&
          lr%bounds%gb%ibyz_ff,lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)
     
  case('P')
     call uncompress_per(n1,n2,n3,lr%wfd%nseg_c,&
          lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1),&
          psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,psir(1))

     call convolut_magic_n_per_self(2*n1+1,2*n2+1,2*n3+1,&
          x_c_psifscf,psir(1))
     
  case('S')
     call uncompress_slab(n1,n2,n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
          lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keyg(1,lr%wfd%nseg_c+1),&
          lr%wfd%keyv(lr%wfd%nseg_c+1),   &
          psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,psir(1))

     call convolut_magic_n_slab_self(2*n1+1,2*n2+15,2*n3+1,x_c_psifscf,psir(1))
  end select

  i_all=-product(shape(x_c_psifscf))*kind(x_c_psifscf)
  deallocate(x_c_psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'x_c_psifscf',subname)
  i_all=-product(shape(x_f_psig))*kind(x_f_psig)
  deallocate(x_f_psig,stat=i_stat)
  call memocc(i_stat,i_all,'x_f_psig',subname)
  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1',subname)
  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2',subname)

  if (lr%geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (lr%geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if


  !call plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,psir)
!HU  call plot_pot_full(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
!HU       nl1,nl2,nl3,orbname,psir,comment)
  call plot_cube_full(1,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
       nl1,nl2,nl3,orbname,psir,comment)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

END SUBROUTINE plot_wf_cube


subroutine plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,pot)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,iounit,n1i,n2i,n3i,nl1,nl2,nl3
  real(gp), intent(in) :: rx,ry,rz,hx,hy,hz
  real(dp), dimension(*), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3,ind
  real(gp) :: hxh,hyh,hzh

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  open(iounit) 
  open(iounit+1) 
  open(iounit+2) 

  i3=nint(rz/hzh)
  i2=nint(ry/hyh)
  write(*,*) 'plot_p, i2,i3,n2,n3 ',i2,i3,n2,n3
  do i1=0,2*n1
     ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
     write(iounit,*) real(i1,gp)*hxh,pot(ind)
  enddo

  i1=nint(rx/hxh)
  i2=nint(ry/hyh)
  write(*,*) 'plot_p, i1,i2 ',i1,i2
  do i3=0,2*n3
     ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
     write(iounit+1,*) real(i3,gp)*hzh,pot(ind)
  enddo

  i1=nint(rx/hxh)
  i3=nint(rz/hzh)
  write(*,*) 'plot_p, i1,i3 ',i1,i3
  do i2=0,2*n2
     ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
     write(iounit+2,*) real(i2,gp)*hyh,pot(ind)
  enddo

  close(iounit) 
  close(iounit+1) 
  close(iounit+2) 

END SUBROUTINE plot_pot


subroutine plot_pot_full(nexpo,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
     nl1,nl2,nl3,orbname,pot,comment)
  use module_base
  implicit none
  character(len=10), intent(in) :: comment
  character(len=11), intent(in) :: orbname
  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,n1i,n2i,n3i,nexpo
  real(gp), intent(in) :: hx,hy,hz
  real(dp), dimension(*), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3,ind
  real(gp) :: hxh,hyh,hzh

!virtual orbitals are identified by their name
!  write(*,'(1x,a,i0)')'printing orbital number= ',iounit
  write(*,'(A)')'printing '//orbname
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz
!  write(orbname,'(i0)')iounit
!  open(unit=22,file='psi'//orbname//'.pot',status='unknown')
  open(unit=22,file=orbname//'.pot',status='unknown')
!  write(22,*)'orbital'//orbname
  write(22,*)orbname
  write(22,*) 2*n1+2,2*n2+2,2*n3+2
  !here there can be shifts in visualising periodic objects (n1 -> n1+1)
  write(22,*) n1*hx,' 0. ',n2*hy
  write(22,*) ' 0. ',' 0. ',n3*hz
  write(22,*)'xyz periodic '//comment
  do i3=0,2*n3+1
     do i2=0,2*n2+1
        do i1=0,2*n1+1
           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           write(22,*)pot(ind)**nexpo
        end do
     end do
  end do
  close(unit=22) 
!  close(unit=23) 

END SUBROUTINE plot_pot_full


subroutine plot_cube_full(nexpo,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
     nl1,nl2,nl3,orbname,pot,comment)
  use module_base
  use module_types
  implicit none
!HU  character(len=10), intent(in) :: orbname,comment
  character(len=10), intent(in) :: comment
  character(len=11), intent(in) :: orbname
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,n1i,n2i,n3i,nexpo
  real(gp), intent(in) :: hx,hy,hz
  real(dp), dimension(*), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3,ind,iat,j,icount
  real(gp) :: hxh,hyh,hzh
  character(len=3) :: advancestring

!virtual orbitals are identified by their name
!  write(*,'(1x,a,i0)')'printing orbital number= ',iounit
  write(*,'(A)')'printing '//orbname
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  open(unit=22,file=trim(orbname)//'.cube',status='unknown')
  write(22,*) trim(orbname)
  write(22,*)'CUBE file for orbital wavefunction'
  !number of atoms
!  if (at%geocode=='P') then
     write(22,'(i5,3(f12.6),a)') at%nat,0.0_gp,0.0_gp,0.0_gp,' modified origin'
!  else if (at%geocode=='S') then
!     write(22,'(i5,3(f12.6))') at%nat,0.0_gp,-hyh,0.0_gp
!  else if (at%geocode=='F') then
!     write(22,'(i5,3(f12.6))') at%nat,-hxh,-hyh,-hzh
!  end if
  !grid and grid spacings
  write(22,'(i5,3(f19.12))') 2*n1+2,hxh,0.0_gp,0.0_gp
  write(22,'(i5,3(f19.12))') 2*n2+2,0.0_gp,hyh,0.0_gp
  write(22,'(i5,3(f19.12))') 2*n3+2,0.0_gp,0.0_gp,hzh
  !atomic number and positions
  do iat=1,at%nat
     write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
  end do

  !the loop is reverted for a cube file
  !charge normalised to the total charge

  do i1=0,2*n1+1
     do i2=0,2*n2+1
        icount=0
        do i3=0,2*n3+1
           icount=icount+1
           if (icount == 6 .or. i3==2*n3+1) then
              advancestring='yes'
              icount=0
           else
              advancestring='no'
           end if

           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           write(22,'(1x,1pe13.6)',advance=advancestring) pot(ind)**nexpo
        end do
     end do
  end do
  close(22)

END SUBROUTINE plot_cube_full


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

END SUBROUTINE plot_psifscf


subroutine read_potfile(geocode,filename,n1,n2,n3,n1i,n2i,n3i,n3d,i3s,rho)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=*), intent(in) :: filename
  integer, intent(in) :: n1i,n2i,n3i,n3d,n1,n2,n3,i3s
  real(dp), dimension(n1i*n2i*n3d), intent(out) :: rho
  !local variables
  integer :: nl1,nl2,nl3,i1,i2,i3,ind
  real(dp) :: value

  open(unit=22,file=filename,status='unknown')
  read(22,*)!'normalised density'
  read(22,*)!2*n1+2,2*n2+2,2*n3+2
  read(22,*)!alat1,' 0. ',alat2
  read(22,*)!' 0. ',' 0. ',alat3
  read(22,*)!xyz   periodic' !not true in general but needed in the case

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=15
     nl3=15
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
  else
     nl2=15
  end if

  call razero(max(n1i*n2i*n3d,1),rho)

  do i3=0,2*n3+1
     do i2=0,2*n2+1
        do i1=0,2*n1+1
           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-i3s)*n1i*n2i
           read(22,*)value
           if (i3+nl3 >= i3s .and. i3+nl3 <= i3s+n3d-1) then
              rho(ind)=value
           end if
        end do
     end do
  end do
  close(22)

END SUBROUTINE read_potfile


subroutine plot_density_old(geocode,filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
     alat1,alat2,alat3,ngatherarr,rho)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nproc
  real(gp), intent(in) :: alat1,alat2,alat3
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(dp), dimension(n1i*n2i*n3p), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density'
  integer :: nl1,nl2,nl3,i_all,i_stat,i1,i2,i3,ind,ierr,nbxz,nby
  real(dp) :: later_avg
  real(dp), dimension(:), pointer :: pot_ion

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (geocode /= 'F') then
     nl1=1
     nl3=1
     nbxz = 1
  else
     nl1=15
     nl3=15
     nbxz = 0
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if

  if (iproc == 0) then
     open(unit=22,file=filename,status='unknown')
     write(22,*)'normalised density'
     write(22,*) 2*(n1+nbxz),2*(n2+nby),2*(n3+nbxz)
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     
     write(22,*)'xyz   periodic' !not true in general but needed in the case
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
     do i3=0,2*(n3+nbxz)-1
        do i2=0,2*(n2+nby)-1
           do i1=0,2*(n1+nbxz)-1
              ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              write(22,*)pot_ion(ind)
           end do
        end do
     end do
     close(22)
     !laterally averaged potential following the isolated direction in surfaces case
     !if (geocode == 'S') then
        open(unit=23,file=filename//'_avg'//geocode,status='unknown')
        do i2=0,2*n2+1
           later_avg=0.0_dp
           do i3=0,2*n3+1
              do i1=0,2*n1+1
                 ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                 later_avg=later_avg+pot_ion(ind)
              end do
           end do
           later_avg=later_avg/real((2*n1+2)*(2*n3+2),dp) !2D integration/2D Volume
           !problem with periodic/isolated BC
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
END SUBROUTINE plot_density_old


subroutine plot_density_cube_old(geocode,filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
     hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
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
  character(len=5) :: suffix
  character(len=15) :: message
  integer :: nl1,nl2,nl3,i_all,i_stat,i1,i2,i3,ind,ierr,icount,j,iat,ia,ib,nbxz,nby
  real(dp) :: a,b
  real(dp), dimension(:,:), pointer :: pot_ion

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (at%geocode /= 'F') then
     nl1=1
     nl3=1
     nbxz = 1
  else
     nl1=15
     nl3=15
     nbxz = 0
  end if
  !value of the buffer in the y direction
  if (at%geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if

  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(pot_ion(n1i*n2i*n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)

     call MPI_ALLGATHERV(rho(1,1),n1i*n2i*n3p,&
          mpidtypd,pot_ion(1,1),ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)

     !case for npspin==2
     if (nspin==2) then
        call MPI_ALLGATHERV(rho(1,2),n1i*n2i*n3p,&
             mpidtypd,pot_ion(1,2),ngatherarr(0,1),&
             ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
     end if

  else
     pot_ion => rho
  end if

  if (iproc == 0) then

     if (nspin /=2) then
        suffix=''
        message='total spin'
        a=1.0_dp
        ia=1
        b=0.0_dp
        ib=1
        call cubefile_write
     else
        suffix='-up'
        message='spin up'
        a=1.0_dp
        ia=1
        b=0.0_dp
        ib=2
        call cubefile_write

        suffix='-down'
        message='spin down'
        a=0.0_dp
        ia=1
        b=1.0_dp
        ib=2
        call cubefile_write

        suffix=''
        message='total spin'
        a=1.0_dp
        ia=1
        b=1.0_dp
        ib=2
        call cubefile_write

        suffix='-u-d'
        message='spin difference'
        a=1.0_dp
        ia=1
        b=-1.0_dp
        ib=2
        call cubefile_write
     end if
  end if

  if (nproc > 1) then
     i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion',subname)
  end if

contains

  subroutine cubefile_write
        open(unit=22,file=trim(filename)//trim(suffix)//'.cube',status='unknown')
        write(22,*)'CUBE file for charge density'
        write(22,*)'Case for '//trim(message)
        !number of atoms
!        if (geocode=='P') then
           write(22,'(i5,3(f12.6),a)') at%nat,0.0_gp,0.0_gp,0.0_gp,' modified'
!        else if (geocode=='S') then
!           write(22,'(i5,3(f12.6))') at%nat,0.0_gp,-hyh,0.0_gp
!        else if (geocode=='F') then
!           write(22,'(i5,3(f12.6))') at%nat,-hxh,-hyh,-hzh
!        end if
        !grid and grid spacings
        write(22,'(i5,3(f12.6))') 2*(n1+nbxz),hxh,0.0_gp,0.0_gp
        write(22,'(i5,3(f12.6))') 2*(n2+nby) ,0.0_gp,hyh,0.0_gp
        write(22,'(i5,3(f12.6))') 2*(n3+nbxz),0.0_gp,0.0_gp,hzh
        !atomic number and positions
        do iat=1,at%nat
           write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
        end do

        !the loop is reverted for a cube file
        !charge normalised to the total charge
        do i1=0,2*(n1+nbxz)-1
           do i2=0,2*(n2+nby)-1
              icount=0
              do i3=0,2*(n3+nbxz)-1
                 icount=icount+1
                 if (icount == 6 .or. i3==2*(n3+nbxz) - 1) then
                    advancestring='yes'
                    icount=0
                 else
                    advancestring='no'
                 end if

                 ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                 write(22,'(1x,1pe13.6)',advance=advancestring)&
                      a*pot_ion(ind,ia)+b*pot_ion(ind,ib)
              end do
           end do
        end do
        close(22)
  END SUBROUTINE cubefile_write

END SUBROUTINE plot_density_cube_old


subroutine read_density_cube_old(filename, n1i,n2i,n3i, nspin, hxh,hyh,hzh, nat, rxyz,  rho)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(out) ::  n1i,n2i,n3i
  integer, intent(in) :: nspin
  real(gp), intent(out) :: hxh,hyh,hzh
  real(gp), pointer :: rxyz(:,:)
  real(dp), dimension(:), pointer :: rho
  integer, intent(out) ::  nat
 
  !local variables
  character(len=*), parameter :: subname='read_density_cube'
  character(len=5) :: suffix
  character(len=15) :: message
  character(len=3) :: advancestring
  integer :: i_all,i_stat,i1,i2,i3,ind,icount,j,iat,ia

  if (nspin /=2) then
     suffix=''
     message='total spin'
     ia=1
     call cubefile_read
  else
     suffix='-up'
     message='spin up'
     ia=1
     call cubefile_read
     
     suffix='-down'
     message='spin down'
     ia=2
     call cubefile_read
     
  end if

contains

  subroutine cubefile_read
       real(dp) dum1,dum2, dum3
        integer idum
        open(unit=22,file=trim(filename)//trim(suffix)//'.cube',status='old')
        read(22,*)! 'CUBE file for charge density'
        read(22,*)! 'Case for '//trim(message)

        read(22,'(i5,3(f12.6),a)')  nat , dum1, dum2, dum3 ! ,0.0_gp,0.0_gp,0.0_gp,' modified'
           

        read(22,'(i5,3(f12.6))') n1i , hxh,   dum1 ,dum2
        read(22,'(i5,3(f12.6))') n2i ,dum1 , hyh  ,  dum2
        read(22,'(i5,3(f12.6))') n3i ,dum1 , dum2 , hzh
        !atomic number and positions

        if( associated(rxyz) ) then
           i_all=-product(shape(rxyz))*kind(rxyz)
           deallocate(rxyz,stat=i_stat)
           call memocc(i_stat,i_all,'rxyz',subname)
        end if
        
        allocate(rxyz(3,nat+ndebug),stat=i_stat)
        call memocc(i_stat,rxyz,'rxyz',subname)
        

        if( associated(rho ).and. ia==1 ) then
           i_all=-product(shape(rho))*kind(rho)
           deallocate(rho,stat=i_stat)
           call memocc(i_stat,i_all,'rho',subname)
        end if
        if(ia==1) then
           allocate(rho(n1i*n2i*n3i+ndebug) ,stat=i_stat)
           call memocc(i_stat,rho,'rho',subname)
        endif

        do iat=1,nat
           read(22,'(i5,4(f12.6))') idum , dum1 , (rxyz(j,iat),j=1,3)
           ! write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
        end do

        do i1=1,n1i
           do i2=1,n2i
              icount=0
              do i3=1,n3i
                 icount=icount+1
                 if (icount == 6  .or. i3==n3i  ) then
                    advancestring='yes'
                    icount=0
                 else
                    advancestring='no'
                 end if

                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i

                 read(22,'(1x,1pe13.6)',advance=advancestring)  rho(ind) 
              end do
           end do
        end do
        close(22)
  END SUBROUTINE cubefile_read

END SUBROUTINE read_density_cube_old


!>   Write a (sum of two) field in the ISF basis in the cube format
!!   Recent changes by Ali:
!!   1) Filling the 2nd column in atomic coordinates rows by the pseudo-cores charge. I found it standard in few other packages. 
!!      In particular it is needed by the recent charge analysis tool.
!!   2) Outputting the electric-dipole moment is an useful piece of data both for the end-user and for developing step 
!!      as a tool to investigate the consistency  (e.g. for symmetrical directions). 
!!      I already did it in this subroutine, but we can do it as a separate subroutine.
subroutine write_cube_fields(filename,message,at,factor,rxyz,n1,n2,n3,n1i,n2i,n3i,n1s,n2s,n3s,hxh,hyh,hzh,&
     a,x,nexpo,b,y)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename,message
  integer, intent(in) :: n1,n2,n3,n1i,n2i,n3i,n1s,n2s,n3s,nexpo
  real(gp), intent(in) :: hxh,hyh,hzh,a,b,factor
  type(atoms_data), intent(in) :: at
  real(wp), dimension(n1i,n2i,n3i), intent(in) :: x,y
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables
  character(len=3) :: advancestring
  integer :: nl1,nl2,nl3,nbx,nby,nbz,i1,i2,i3,icount,j,iat
  real(dp) :: later_avg
  real(gp) :: dipole_el(3) , dipole_cores(3),q
  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (at%geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
  end if
  !value of the buffer in the y direction
  if (at%geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if

  open(unit=22,file=trim(filename)//'.cube',status='unknown')
  write(22,*)'CUBE file for ISF field'
  write(22,*)'Case for '//trim(message)
  write(22,'(i5,3(f12.6))') at%nat,0.0_gp,0.0_gp,0.0_gp
  !grid and grid spacings
  write(22,'(i5,3(f12.6))') 2*(n1+nbx),hxh,0.0_gp,0.0_gp
  write(22,'(i5,3(f12.6))') 2*(n2+nby),0.0_gp,hyh,0.0_gp
  write(22,'(i5,3(f12.6))') 2*(n3+nbz),0.0_gp,0.0_gp,hzh
  !atomic number and positions
  dipole_el   (1:3)=0_gp
  dipole_cores(1:3)=0_gp
  do iat=1,at%nat
     !write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
     write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)), at%nelpsp(at%iatype(iat))*1. &
          ,(rxyz(j,iat),j=1,3)
     dipole_cores(1:3)=dipole_cores(1:3)+at%nelpsp(at%iatype(iat)) * rxyz(1:3,iat)
  end do


  !the loop is reverted for a cube file
  !charge normalised to the total charge
  do i1=0,2*(n1+nbx) - 1
     do i2=0,2*(n2+nby) - 1
        icount=0
        do i3=0,2*(n3+nbz) - 1
           icount=icount+1
           if (icount == 6 .or. i3==2*(n3+nbz) - 1) then
              advancestring='yes'
              icount=0
           else
              advancestring='no'
           end if
           !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           write(22,'(1x,1pe13.6)',advance=advancestring)&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
           q= ( a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3) )* hxh*hyh*hzh 
           dipole_el(1)=dipole_el(1)+ q* at%alat1/real(2*(n1+nbx),dp)*i1 
           dipole_el(2)=dipole_el(2)+ q* at%alat2/real(2*(n2+nby),dp)*i2
           dipole_el(3)=dipole_el(3)+ q* at%alat3/real(2*(n3+nbz),dp)*i3
        end do
     end do
  end do
  close(22)
  !average in x direction
  open(unit=23,file=trim(filename)//'_avg_x',status='unknown')
  !  do i1=0,2*n1+1
  do i1=0,2*(n1+nbx) - 1
     later_avg=0.0_dp
     do i3=0,2*(n3+nbz) -1
        do i2=0,2*(n2+nby) - 1
           later_avg=later_avg+&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
     later_avg=later_avg/real((2*(n2+nby))*(2*(n3+nbz)),dp) !2D integration/2D Volume
     !to be checked with periodic/isolated BC
     write(23,*)i1+n1s,at%alat1/real(factor*2*(n1+nbx),dp)*(i1+2*n1s),later_avg
  end do
  close(23)
  !average in y direction
  open(unit=23,file=trim(filename)//'_avg_y',status='unknown')
  do i2=0,2*(n2+nby) - 1
     later_avg=0.0_dp
     do i3=0,2*(n3+nbz) - 1
        do i1=0,2*(n1+nbx) -1 
           later_avg=later_avg+&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
     later_avg=later_avg/real((2*(n1+nbx))*(2*(n3+nbz)),dp) !2D integration/2D Volume
     !to be checked with periodic/isolated BC
     write(23,*)i2+n2s,at%alat2/real(factor*2*(n2+nby),dp)*(i2+n2s),later_avg
  end do
  close(23)
  !average in z direction
  open(unit=23,file=trim(filename)//'_avg_z',status='unknown')
  do i3=0,2*(n3+nbz) - 1
     later_avg=0.0_dp
     do i2=0,2*(n2+nby) - 1
        do i1=0,2*(n1+nbx) -1 
           later_avg=later_avg+&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
     later_avg=later_avg/real((2*(n1+nbx))*(2*(n2+nby)),dp) !2D integration/2D Volume
     !to be checked with periodic/isolated BC
     write(23,*)i3+n3s,at%alat3/real(factor*2*(n3+nbz)+2,dp)*(i3+n3s),later_avg
  end do
  close(23)
  if (trim(filename)=='electronic_density') then
     dipole_el=dipole_el        !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     dipole_cores=dipole_cores  !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     open(unit=24,file='dipole',status='unknown')
     write(24,'(a)') " #  Dipole moment of the whole system  (Px, Py, Pz,  |P| [e.bohr])"  ! or [D] 
     write(24,99) "electronic charge: ", dipole_el(1:3) , sqrt(sum(dipole_el**2))
     write(24,99) "pseudo cores:      ", dipole_cores(1:3) , sqrt(sum(dipole_cores**2))
     write(24,99) "Total (cores-el.): ", dipole_cores-dipole_el , sqrt(sum((dipole_cores-dipole_el)**2))
99   format (a20,3f15.7,"    ==> ",f15.5)
     !99 format (a20,4ES15.7)
     close(24)
  endif
END SUBROUTINE write_cube_fields


subroutine plot_density(filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
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
  real(dp), dimension(max(n1i*n2i*n3p,1),nspin), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density'
  character(len=5) :: suffix
  character(len=65) :: message
  integer :: i_all,i_stat,ierr,ia,ib,isuffix,fformat
  real(dp) :: a,b
  real(dp), dimension(:,:), pointer :: pot_ion

  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(pot_ion(n1i*n2i*n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)

     call MPI_ALLGATHERV(rho(1,1),n1i*n2i*n3p,&
          mpidtypd,pot_ion(1,1),ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)

     !case for npspin==2
     if (nspin==2) then
        call MPI_ALLGATHERV(rho(1,2),n1i*n2i*n3p,&
             mpidtypd,pot_ion(1,2),ngatherarr(0,1),&
             ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
     end if

  else
     pot_ion => rho
  end if

  ! Format = 1 -> cube (default)
  ! Format = 2 -> ETSF
  ! ...
  fformat = 1
  isuffix = index(filename, ".cube", back = .true.)
  if (isuffix > 0) then
     isuffix = isuffix - 1
     fformat = 1
  else
     isuffix = index(filename, ".etsf", back = .true.)
     if (isuffix <= 0) isuffix = index(filename, ".etsf.nc", back = .true.)
     if (isuffix > 0) then
        isuffix = isuffix - 1
        fformat = 2
     else
        isuffix = len(filename)
     end if
  end if

  if (iproc == 0) then

     if (nspin /=2) then
        message='total spin'
        if (fformat == 1) then
           suffix=''
           a=1.0_dp
           ia=1
           b=0.0_dp
           ib=1
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1,n2,n3,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
        else
           call write_etsf_density(filename(:isuffix),message,&
                at,rxyz,n1,n2,n3,n1i,n2i,n3i,hxh,hyh,hzh,&
                pot_ion, 1)
        end if
     else
        if (fformat == 1) then
           suffix='-up'
           message='spin up'
           a=1.0_dp
           ia=1
           b=0.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1,n2,n3,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))

           suffix='-down'
           message='spin down'
           a=0.0_dp
           ia=1
           b=1.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1,n2,n3,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))

           suffix=''
           message='total spin'
           a=1.0_dp
           ia=1
           b=1.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1,n2,n3,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))

           suffix='-u-d'
           message='spin difference'
           a=1.0_dp
           ia=1
           b=-1.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1,n2,n3,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
        else
           message = 'spin up, down, total, difference'
           call write_etsf_density(filename(:isuffix),message,&
                at,rxyz,n1,n2,n3,n1i,n2i,n3i,hxh,hyh,hzh,&
                pot_ion, 2)
        end if

     end if

  end if


  if (nproc > 1) then
     i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion',subname)
  end if

END SUBROUTINE plot_density


!>  Read a density file using file format depending on the extension.
subroutine read_density(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz,iatypes, znucl)
  use module_base
  use module_types
  use module_interfaces, except_this_one => read_density
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(out) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:), pointer :: rho
  real(gp), dimension(:,:), pointer, optional :: rxyz
  integer, intent(out), optional ::  nat
  integer, dimension(:), pointer, optional :: iatypes, znucl

  character(len = *), parameter :: subname = "read_density"
  integer :: isuffix,fformat,nat_read, i_stat, i_all
  real(gp), dimension(:,:), pointer :: rxyz_read
  integer, dimension(:), pointer :: iatypes_read, znucl_read

  !check the arguments
  if (.not.(present(rxyz) .and. present(nat) .and. present(iatypes) .and. present(znucl)) .and. &
       & (present(rxyz) .or. present(nat) .or. present(iatypes) .or. present(znucl))) then
     stop 'wrong usage of read_densityt, rxyz, znucl and iatypes should be present'
  end if

  ! Format = 1 -> cube (default)
  ! Format = 2 -> ETSF
  ! ...
  fformat = 1
  isuffix = index(filename, ".cube", back = .true.)
  if (isuffix > 0) then
     isuffix = isuffix - 1
     fformat = 1
  else
     isuffix = index(filename, ".etsf", back = .true.)
     if (isuffix <= 0) isuffix = index(filename, ".etsf.nc", back = .true.)
     if (isuffix > 0) then
        isuffix = isuffix - 1
        fformat = 2
     else
        isuffix = len(filename)
     end if
  end if

  if (fformat == 1) then
     call read_cube(filename(1:isuffix),geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
          nat_read,rxyz_read, iatypes_read, znucl_read)
  else
     call read_etsf(filename(1:isuffix),geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
          nat_read,rxyz_read, iatypes_read, znucl_read)
  end if

  if (present(rxyz) .and. present(nat) .and. present(iatypes) .and. present(znucl)) then
     rxyz => rxyz_read
     iatypes => iatypes_read
     znucl => znucl_read
     nat=nat_read
  else
     i_all=-product(shape(rxyz_read))*kind(rxyz_read)
     deallocate(rxyz_read,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz_read',subname)
     i_all=-product(shape(iatypes_read))*kind(iatypes_read)
     deallocate(iatypes_read,stat=i_stat)
     call memocc(i_stat,i_all,'iatypes_read',subname)
     i_all=-product(shape(znucl_read))*kind(znucl_read)
     deallocate(znucl_read,stat=i_stat)
     call memocc(i_stat,i_all,'znucl_read',subname)
  end if
END SUBROUTINE read_density


subroutine plot_wf(orbname,nexpo,at,factor,lr,hx,hy,hz,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
  character(len=*) :: comment
  character(len=*) :: orbname
  integer, intent(in) :: nexpo
  real(dp), intent(in) :: factor
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: i_stat,i_all
  integer :: n1i,n2i,n3i,n1,n2,n3,n1s,n2s,n3s
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: psir

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1s=lr%ns1
  n2s=lr%ns2
  n3s=lr%ns3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i

  call initialize_work_arrays_sumrho(lr,w)

  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if

  call daub_to_isf(lr,w,psi,psir)

  call write_cube_fields(orbname,' ',&
       at,factor,rxyz,n1,n2,n3,n1i,n2i,n3i,n1s,n2s,n3s,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
       1.0_gp,psir,nexpo,0.0_gp,psir)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf


!> Read density or potential in cube format
subroutine read_cube(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz, iatypes, znucl)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(out) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:), pointer :: rho
  real(gp), dimension(:,:), pointer   :: rxyz
  integer, intent(out)   ::  nat
  integer, dimension(:), pointer   :: iatypes, znucl
  !local variables
  character(len=*), parameter :: subname='read_cube'
  character(len=5) :: suffix
  character(len=15) :: message
  integer :: ia
  logical :: exists

  ! Test if we have up and down densities.
  inquire(file=trim(filename)//"-up.cube",exist=exists)
  if (exists) then
     inquire(file=trim(filename)//"-down.cube",exist=exists)
     if (.not.exists) then
        write(*,*) "WARNING! found a -up.cube file, but no -down.cube..."
        nspin = 1
     else
        nspin = 2
     end if
  else
     nspin = 1
  end if

  if (nspin /=2) then
     suffix=''
     message='total spin'
     ia=1
     !read the header of the file
     call read_cube_header(filename//trim(suffix),geocode,nspin,n1i,n2i,n3i,hxh,hyh,hzh,&
          rho,nat,rxyz, iatypes, znucl) 
     !fill the pointer which was just allocated
     call read_cube_field(filename//trim(suffix),geocode,n1i,n2i,n3i,rho)

  else
     suffix='-up'
     message='spin up'
     ia=1
     !read the header of the file (suppose that it is the same for spin up and down)
     call read_cube_header(filename//trim(suffix),geocode,nspin,n1i,n2i,n3i,hxh,hyh,hzh,&
          rho,nat,rxyz, iatypes, znucl) 
     !fill the pointer which was just allocated
     call read_cube_field(filename//trim(suffix),geocode,n1i,n2i,n3i,rho(1,ia))
    
     suffix='-down'
     message='spin down'
     ia=2
     call read_cube_field(filename//trim(suffix),geocode,n1i,n2i,n3i,rho(1,ia)) 
  end if

contains

  subroutine read_cube_header(filename,geocode,nspin,n1i,n2i,n3i,hxh,hyh,hzh,rho,&
       nat,rxyz, iatypes, znucl)
    use module_base
    use module_types
    implicit none
    character(len=*), intent(in) :: filename
    character(len=1), intent(in) :: geocode
    integer, intent(in) :: nspin
    integer, intent(out) :: n1i,n2i,n3i
    real(gp), intent(out) :: hxh,hyh,hzh
    real(dp), dimension(:,:), pointer :: rho
    integer, intent(out) :: nat
    real(gp), dimension(:,:), pointer :: rxyz
    integer, dimension(:), pointer :: iatypes, znucl
    !local variables
    character(len=*), parameter :: subname='read_density_cube'
    integer :: n1t,n2t,n3t,n1,n2,n3,idum,iat,i_stat,i_all,j
    integer :: nl1,nl2,nl3,nbx,nby,nbz
    real(gp) :: dum1,dum2,dum3
    integer, dimension(:), allocatable :: znucl_

    if (geocode /= 'F') then
       nl1=1
       nl3=1
       nbx = 1
       nbz = 1
    else
       nl1=15
       nl3=15
       nbx = 0
       nbz = 0
    end if
    !value of the buffer in the y direction
    if (geocode == 'P') then
       nl2=1
       nby = 1
    else
       nl2=15
       nby = 0
    end if

    open(unit=22,file=trim(filename)//".cube",status='old')
    read(22,*)! 'CUBE file for charge density'
    read(22,*)! 'Case for '//trim(message)

    read(22,'(i5,3(f12.6),a)')  nat , dum1, dum2, dum3 
    read(22,'(i5,3(f12.6))') n1t , hxh  , dum1 , dum2
    read(22,'(i5,3(f12.6))') n2t , dum1 , hyh  , dum2
    read(22,'(i5,3(f12.6))') n3t , dum1 , dum2 , hzh

    !grid positions
    n1=n1t/2-nbx
    n1i=2*n1+(1-nbx)+2*nl1
    n2=n2t/2-nby
    n2i=2*n2+(1-nby)+2*nl2
    n3=n3t/2-nbz
    n3i=2*n3+(1-nbz)+2*nl3

    !atomic positions
    allocate(rxyz(3,nat+ndebug),stat=i_stat)
    call memocc(i_stat,rxyz,'rxyz',subname)
    allocate(iatypes(nat+ndebug),stat=i_stat)
    call memocc(i_stat,iatypes,'iatypes',subname)
    allocate(znucl_(nat+ndebug),stat=i_stat)
    call memocc(i_stat,znucl_,'znucl_',subname)
    znucl_(:) = -1

    if(associated(rho)) then
       i_all=-product(shape(rho))*kind(rho)
       deallocate(rho,stat=i_stat)
       call memocc(i_stat,i_all,'rho',subname)
    end if

    allocate(rho(n1i*n2i*n3i,nspin+ndebug) ,stat=i_stat)
    call memocc(i_stat,rho,'rho',subname)

    do iat=1,nat
       read(22,'(i5,4(f12.6))') idum , dum1 , (rxyz(j,iat),j=1,3)
       do j = 1, nat, 1
          if (znucl_(j) == idum .or. znucl_(j) == -1) then
             znucl_(j) = idum
             exit
          end if
       end do
       iatypes(iat) = j
       ! write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
    end do

    do j = 1, nat, 1
       if (znucl_(j) == -1) then
          exit
       end if
    end do
    allocate(znucl(j-1+ndebug),stat=i_stat)
    call memocc(i_stat,znucl,'znucl',subname)
    znucl(1:j-1) = znucl_(1:j-1)

    i_all=-product(shape(znucl_))*kind(znucl_)
    deallocate(znucl_,stat=i_stat)
    call memocc(i_stat,i_all,'znucl_',subname)

    close(22)
     
  END SUBROUTINE read_cube_header

END SUBROUTINE read_cube


!>   Read a cube field which have been plotted previously by write_cube_fields
subroutine read_cube_field(filename,geocode,n1i,n2i,n3i,rho)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1i,n2i,n3i
  real(dp), dimension(n1i*n2i*n3i) :: rho
  !local variables
  character(len=*), parameter :: subname='read_density_cube'
  character(len=3) :: advancestring
  integer :: n1t,n2t,n3t,n1,n2,n3,i1,i2,i3,nat,iat
  integer :: nl1,nl2,nl3,nbx,nby,nbz,icount,ind
  real(gp) :: dum1,dum2,dum3,tt

  if (geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if

  open(unit=22,file=trim(filename)//'.cube',status='old')
  read(22,*)! 'CUBE file for charge density'
  read(22,*)! 'Case for '//trim(message)

  read(22,'(i5,3(f12.6),a)')  nat , dum1, dum2, dum3 
  read(22,'(i5,3(f12.6))') n1t , dum3,   dum1 ,dum2
  read(22,'(i5,3(f12.6))') n2t ,dum1 , dum3  ,  dum2
  read(22,'(i5,3(f12.6))') n3t ,dum1 , dum2 , dum3

  !grid positions
  n1=n1t/2-nbx
  if (n1i /= 2*n1+(1-nbx)+2*nl1) stop 'n1i not valid'
  n2=n2t/2-nby
  if (n2i /= 2*n2+(1-nby)+2*nl2) stop 'n2i not valid'
  n3=n3t/2-nbz
  if (n3i /= 2*n3+(1-nbz)+2*nl3) stop 'n3i not valid'

  !zero the pointer
  call razero(n1i*n2i*n3i,rho)

  do iat=1,nat
     !read(22,'(i5,4(f12.6))')! idum , dum1 , (rxyz(j,iat),j=1,3)
     read(22,*)! idum , dum1 , (rxyz(j,iat),j=1,3)
  end do

  !the loop is reverted for a cube file
  !charge normalised to the total charge
  do i1=0,2*(n1+nbx) - 1
     do i2=0,2*(n2+nby) - 1
        icount=0
        do i3=0,2*(n3+nbz) - 1
           icount=icount+1
           if (icount == 6 .or. i3==2*(n3+nbz) - 1) then
              advancestring='yes'
              icount=0
           else
              advancestring='no'
           end if
           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           read(22,'(1x,1pe13.6)',advance=advancestring) tt !rho(ind) 
           rho(ind)=tt
!           write(16,*)i1,i2,i3,ind,rho(ind)
           !read(22,*)',advance=advancestring) rho(ind) 
        end do
     end do
  end do
  close(22)

 ! write(14,*)rho

END SUBROUTINE read_cube_field






subroutine plot_wfSquare_cube(orbname,at,lr,hx,hy,hz,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
!HU  character(len=10) :: orbname,comment 
  character(len=10) :: comment
  character(len=11) :: orbname
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
  !local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: nw1,nw2,i_stat,i_all,i,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  integer :: nxc,nxf,nl1,nl2,nl3,n1i,n2i,n3i
  real(gp) :: rx, ry, rz
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:), allocatable :: psir,w1,w2,x_c_psifscf,x_f_psig

  rx=rxyz(1,1)
  ry=rxyz(2,1)
  rz=rxyz(3,1)

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3

  do i=0,3
     scal(i)=1.d0
  enddo

  select case(lr%geocode)
  case('F')
     !dimension of the work arrays
     ! shrink convention: nw1>nw2
     nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))
     nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
  case('S')
     !dimension of the work arrays
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
     nxf=1
  case('P')
     !dimension of the work arrays, fully periodic case
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=1
  end select
  !work arrays
  allocate(x_c_psifscf(nxc+ndebug),stat=i_stat)
  call memocc(i_stat,x_c_psifscf,'x_c_psifscf',subname)
  allocate(x_f_psig(nxf+ndebug),stat=i_stat)
  call memocc(i_stat,x_f_psig,'x_f_psig',subname)
  allocate(w1(nw1+ndebug),stat=i_stat)
  call memocc(i_stat,w1,'w1',subname)
  allocate(w2(nw2+ndebug),stat=i_stat)
  call memocc(i_stat,w2,'w2',subname)

  allocate(psir(n1i*n2i*n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(nxc,x_c_psifscf)
     call razero(nxf,x_f_psig)
     call razero(n1i*n2i*n3i,psir)
  end if
  select case(lr%geocode)
  case('F')
     call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),  & 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1), &
          scal,psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,x_f_psig)

     call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,&
          x_c_psifscf,x_f_psig,  & 
          psir(1),lr%bounds%kb%ibyz_c,lr%bounds%gb%ibzxx_c,lr%bounds%gb%ibxxyy_c,&
          lr%bounds%gb%ibyz_ff,lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)
     
  case('P')
     call uncompress_per(n1,n2,n3,lr%wfd%nseg_c,&
          lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1),&
          psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,psir(1))

     call convolut_magic_n_per_self(2*n1+1,2*n2+1,2*n3+1,&
          x_c_psifscf,psir(1))
     
  case('S')
     call uncompress_slab(n1,n2,n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
          lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keyg(1,lr%wfd%nseg_c+1),&
          lr%wfd%keyv(lr%wfd%nseg_c+1),   &
          psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,psir(1))

     call convolut_magic_n_slab_self(2*n1+1,2*n2+15,2*n3+1,x_c_psifscf,psir(1))
  end select

  i_all=-product(shape(x_c_psifscf))*kind(x_c_psifscf)
  deallocate(x_c_psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'x_c_psifscf',subname)
  i_all=-product(shape(x_f_psig))*kind(x_f_psig)
  deallocate(x_f_psig,stat=i_stat)
  call memocc(i_stat,i_all,'x_f_psig',subname)
  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1',subname)
  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2',subname)

  if (lr%geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (lr%geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if


  !call plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,psir)
!HU  call plot_pot_full(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
!HU       nl1,nl2,nl3,orbname,psir,comment)
  call plot_cube_full(2,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
       nl1,nl2,nl3,orbname,psir,comment)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

END SUBROUTINE plot_wfSquare_cube
