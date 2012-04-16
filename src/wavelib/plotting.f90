!> @file
!!  Routines to plot wavefunctions
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine plot_density_cube_old(filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
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
  character(len=*), parameter :: subname='plot_density_cube_old'
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
  character(len=*), parameter :: subname='read_density_cube_old'
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


!> Write a (sum of two) field in the ISF basis in the cube format
subroutine write_cube_fields(filename,message,at,factor,rxyz,n1i,n2i,n3i,n1s,n2s,n3s,hxh,hyh,hzh,&
     a,x,nexpo,b,y)
  !n(c) use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename,message
  integer, intent(in) :: n1i,n2i,n3i,n1s,n2s,n3s,nexpo
  real(gp), intent(in) :: hxh,hyh,hzh,a,b,factor
  type(atoms_data), intent(in) :: at
  real(wp), dimension(n1i,n2i,n3i), intent(in) :: x,y
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables
  character(len=3) :: advancestring
  integer :: nl1,nl2,nl3,nbx,nby,nbz,i1,i2,i3,icount,j,iat,nc1,nc2,nc3
  real(dp) :: later_avg
  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (at%geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
     nc1=n1i-31
     nc3=n3i-31
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
     nc1=n1i
     nc3=n3i
  end if
  !value of the buffer in the y direction
  if (at%geocode == 'P') then
     nl2=1
     nby = 1
     nc2=n2i-31
  else
     nl2=15
     nby = 0
     nc2=n2i
  end if

!!$  !cube dimensions
!!$  nc1=2*(n1+nbx)
!!$  nc2=2*(n2+nby)
!!$  nc3=2*(n3+nbz)

  open(unit=22,file=trim(filename)//'.cube',status='unknown')
  write(22,*)'CUBE file for ISF field'
  write(22,*)'Case for '//trim(message)
  write(22,'(i5,3(f12.6))') at%nat,0.0_gp,0.0_gp,0.0_gp
  !grid and grid spacings
  write(22,'(i5,3(f12.6))') nc1,hxh,0.0_gp,0.0_gp
  write(22,'(i5,3(f12.6))') nc2,0.0_gp,hyh,0.0_gp
  write(22,'(i5,3(f12.6))') nc3,0.0_gp,0.0_gp,hzh
  !atomic number and positions
  do iat=1,at%nat
     write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)), at%nelpsp(at%iatype(iat))*1. &
          ,(rxyz(j,iat),j=1,3)
  end do


  !the loop is reverted for a cube file
  !charge normalised to the total charge
  do i1=0,nc1 - 1
     do i2=0,nc2 - 1
        icount=0
        do i3=0,nc3 - 1
           icount=icount+1
           if (icount == 6 .or. i3==nc3 - 1) then
              advancestring='yes'
              icount=0
           else
              advancestring='no'
           end if
           !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           write(22,'(1x,1pe13.6)',advance=advancestring)&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
  end do
  close(22)
  !average in x direction
  open(unit=23,file=trim(filename)//'_avg_x',status='unknown')
  !  do i1=0,2*n1+1
  do i1=0,nc1 - 1
     later_avg=0.0_dp
     do i3=0,nc3 -1
        do i2=0,nc2 - 1
           later_avg=later_avg+&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
     later_avg=later_avg/real(nc2*nc3,dp) !2D integration/2D Volume
     !to be checked with periodic/isolated BC
     write(23,*)i1+n1s,at%alat1/real(factor*nc1,dp)*(i1+2*n1s),later_avg
  end do
  close(23)
  !average in y direction
  open(unit=23,file=trim(filename)//'_avg_y',status='unknown')
  do i2=0,nc2 - 1
     later_avg=0.0_dp
     do i3=0,nc3 - 1
        do i1=0,nc1 -1 
           later_avg=later_avg+&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
     later_avg=later_avg/real(nc1*nc3,dp) !2D integration/2D Volume
     !to be checked with periodic/isolated BC
     write(23,*)i2+n2s,at%alat2/real(factor*nc2,dp)*(i2+n2s),later_avg
  end do
  close(23)
  !average in z direction
  open(unit=23,file=trim(filename)//'_avg_z',status='unknown')
  do i3=0,nc3 - 1
     later_avg=0.0_dp
     do i2=0,nc2 - 1
        do i1=0,nc1 -1 
           later_avg=later_avg+&
                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
        end do
     end do
     later_avg=later_avg/real(nc1*nc2,dp) !2D integration/2D Volume
     !to be checked with periodic/isolated BC
     write(23,*)i3+n3s,at%alat3/real(factor*nc3+2,dp)*(i3+n3s),later_avg
  end do
  close(23)
END SUBROUTINE write_cube_fields


subroutine plot_density(iproc,nproc,filename,at,rxyz,box,nspin,rho)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  type(atoms_data), intent(in) :: at
  type(denspot_distribution), intent(in) :: box
  character(len=*), intent(in) :: filename
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(max(box%ndimpot,1),nspin), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density'
  character(len=5) :: suffix
  character(len=65) :: message
  integer :: i_all,i_stat,ierr,ia,ib,isuffix,fformat,n1i,n2i,n3i
  real(dp) :: a,b
  real(gp) :: hxh,hyh,hzh
  real(dp), dimension(:,:), pointer :: pot_ion

  n1i=box%ndims(1)
  n2i=box%ndims(2)
  n3i=box%ndims(3)

  hxh=box%hgrids(1)
  hyh=box%hgrids(2)
  hzh=box%hgrids(3)

  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(pot_ion(box%ndimgrid,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)

     call MPI_ALLGATHERV(rho(1,1),box%ndimpot,&
          mpidtypd,pot_ion(1,1),box%ngatherarr(0,1),&
          box%ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)

     !case for npspin==2
     if (nspin==2) then
        call MPI_ALLGATHERV(rho(1,2),box%ndimpot,&
             mpidtypd,pot_ion(1,2),box%ngatherarr(0,1),&
             box%ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
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
        isuffix = len(trim(filename))
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
                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
        else
           call write_etsf_density(filename(:isuffix),message,&
                at,rxyz,n1i,n2i,n3i,hxh,hyh,hzh,&
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
                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))

           suffix='-down'
           message='spin down'
           a=0.0_dp
           ia=1
           b=1.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))

           suffix=''
           message='total spin'
           a=1.0_dp
           ia=1
           b=1.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))

           suffix='-u-d'
           message='spin difference'
           a=1.0_dp
           ia=1
           b=-1.0_dp
           ib=2
           call write_cube_fields(filename(:isuffix)//trim(suffix),message,&
                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
        else
           message = 'spin up, down, total, difference'
           call write_etsf_density(filename(:isuffix),message,&
                at,rxyz,n1i,n2i,n3i,hxh,hyh,hzh,&
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


!> Read a density file using file format depending on the extension.
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
  real(dp), dimension(:,:,:,:), pointer :: rho
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
     stop 'wrong usage of read_density, rxyz, znucl and iatypes should be present'
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


subroutine plot_wf(orbname,nexpo,at,factor,lr,hx,hy,hz,rxyz,psi)
  use module_base
  use module_types
  implicit none
  !Arguments
  character(len=*) :: orbname
  integer, intent(in) :: nexpo
  real(dp), intent(in) :: factor
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  !Local variables
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
       at,factor,rxyz,n1i,n2i,n3i,n1s,n2s,n3s,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
       1.0_gp,psir,nexpo,0.0_gp,psir)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf


!> Read the densit and put the values in the rhopot arrays according to the parallelization indicated by
!! nscatterarr array
subroutine read_potential_from_disk(iproc,nproc,filename,geocode,ngatherarr,n1i,n2i,n3i,n3p,nspin,hxh,hyh,hzh,pot)
  use module_base
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,n3p,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(dp), dimension(n1i,n2i,max(n3p,1),nspin), intent(out) :: pot
  !local variables
  character(len=*), parameter :: subname='read_potential_from_disk'
  integer :: n1t,n2t,n3t,nspint,ierror,ierr,i_all,i_stat,ispin
  real(gp) :: hxt,hyt,hzt
  real(dp), dimension(:,:,:,:), pointer :: pot_from_disk

  !only the first processor should read this
  if (iproc == 0) then
     write(*,'(1x,a)')'Reading local potential from file:'//trim(filename)
     call read_density(trim(filename),geocode,&
          n1t,n2t,n3t,nspint,hxt,hyt,hzt,pot_from_disk)
     if (abs(hxt-hxh) <= 1.e-5_gp .and. abs(hyt-hyh) <= 1.e-5_gp .and. abs(hzt-hzh) <= 1.e-5_gp .and. &
          nspint == nspin .and. &
          n1i  == n1t  .and. n2i == n2t .and. n3i == n3t) then
     else
        write(*,*)'ERROR (to be documented): some of the parameters do not coincide'
        write(*,*)hxh,hyh,hzh,hxt,hyt,hzt,nspin,nspint,n1i,n2i,n3i,n1t,n2t,n3t
        call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
     end if
  else
     allocate(pot_from_disk(1,1,1,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot_from_disk,'pot_from_disk',subname)
  end if

  if (nproc > 1) then
     do ispin=1,nspin
        call MPI_SCATTERV(pot_from_disk(1,1,1,ispin),&
             ngatherarr(0,1),ngatherarr(0,2),mpidtypd, &
             pot(1,1,1,ispin),&
             n1i*n2i*n3p,mpidtypd,0,MPI_COMM_WORLD,ierr)
     end do
  else
     call vcopy(n1i*n2i*n3i*nspin,pot_from_disk(1,1,1,1),1,pot(1,1,1,1),1)
  end if

  i_all=-product(shape(pot_from_disk))*kind(pot_from_disk)
  deallocate(pot_from_disk,stat=i_stat)
  call memocc(i_stat,i_all,'pot_from_disk',subname)

end subroutine read_potential_from_disk


!> Read density or potential in cube format
subroutine read_cube(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz, iatypes, znucl)
  !n(c) use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(out) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:,:,:), pointer :: rho
  real(gp), dimension(:,:), pointer   :: rxyz
  integer, intent(out)   ::  nat
  integer, dimension(:), pointer   :: iatypes, znucl
  !local variables
  !n(c) character(len=*), parameter :: subname='read_cube'
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
     call read_cube_field(filename//trim(suffix),geocode,n1i,n2i,n3i,rho(1,1,1,ia))
    
     suffix='-down'
     message='spin down'
     ia=2
     call read_cube_field(filename//trim(suffix),geocode,n1i,n2i,n3i,rho(1,1,1,ia)) 
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
    real(dp), dimension(:,:,:,:), pointer :: rho
    integer, intent(out) :: nat
    real(gp), dimension(:,:), pointer :: rxyz
    integer, dimension(:), pointer :: iatypes, znucl
    !local variables
    character(len=*), parameter :: subname='read_cube_header'
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

    allocate(rho(n1i,n2i,n3i,nspin+ndebug) ,stat=i_stat)
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
  !n(c) use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1i,n2i,n3i
  real(dp), dimension(n1i*n2i*n3i) :: rho
  !local variables
  !n(c) character(len=*), parameter :: subname='read_cube_field'
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


!> Calculate the dipole of a Field given in the rho array.
!! The parallel distribution used is the one of the potential
subroutine calc_dipole(iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin, &
   hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nspin,nproc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(n1i,n2i,max(n3p, 1),nspin), target, intent(in) :: rho
  character(len=*), parameter :: subname='calc_dipole'
  integer :: i_all,i_stat,ierr
  real(gp) :: dipole_el(3) , dipole_cores(3), tmpdip(3),q,qtot
  integer  :: iat,i1,i2,i3,nbx,nby,nbz, nl1,nl2,nl3, ispin
  real(dp), dimension(:,:,:,:), pointer :: ele_rho,rho_buf
  
  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(ele_rho(n1i,n2i,n3i,nspin),stat=i_stat)
     call memocc(i_stat,ele_rho,'ele_rho',subname)

!Commented out, it is enough to allocate the rho at 1
!!$     ! rho_buf is used instead of rho for avoiding the case n3p=0 in 
!!$     ! some procs which makes MPI_ALLGATHERV failed.
!!$     if (n3p.eq.0) then
!!$       allocate(rho_buf(n1i,n2i,n3p+1,nspin),stat=i_stat)
!!$       call memocc(i_stat,rho_buf,'rho_buf',subname)
!!$       rho_buf = 0.0_dp
!!$     else
!!$       allocate(rho_buf(n1i,n2i,n3p,nspin),stat=i_stat)
!!$       call memocc(i_stat,rho_buf,'rho_buf',subname)
!!$       rho_buf = rho
!!$     endif  

     do ispin=1,nspin
        call MPI_ALLGATHERV(rho(1,1,1,ispin),n1i*n2i*n3p,&
             mpidtypd,ele_rho(1,1,1,ispin),ngatherarr(0,1),&
             ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
     end do

  else
     ele_rho => rho
  end if

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

  qtot=0.d0
  dipole_cores(1:3)=0_gp
  do iat=1,at%nat
     dipole_cores(1:3)=dipole_cores(1:3)+at%nelpsp(at%iatype(iat)) * rxyz(1:3,iat)
  end do

  dipole_el   (1:3)=0_gp
  do ispin=1,nspin
     do i1=0,2*(n1+nbx) - 1
        do i2=0,2*(n2+nby) - 1
           do i3=0,2*(n3+nbz) - 1
              !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              !q= ( ele_rho(ind,ispin) ) * hxh*hyh*hzh 
              q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * hxh*hyh*hzh 
              qtot=qtot+q
              dipole_el(1)=dipole_el(1)+ q* at%alat1/real(2*(n1+nbx),dp)*i1 
              dipole_el(2)=dipole_el(2)+ q* at%alat2/real(2*(n2+nby),dp)*i2
              dipole_el(3)=dipole_el(3)+ q* at%alat3/real(2*(n3+nbz),dp)*i3
           end do
        end do
     end do

  end do

  if(iproc==0) then
     !dipole_el=dipole_el        !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     !dipole_cores=dipole_cores  !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     write(*,'(1x,a)')repeat('-',61)//' Electric Dipole Moment'
     tmpdip=dipole_cores+dipole_el
     write(*,96) "|P| = ", sqrt(sum(tmpdip**2)), " (AU)       ", "(Px,Py,Pz)= " , tmpdip(1:3)  
     tmpdip=tmpdip/0.393430307_gp  ! au2debye              
     write(*,96) "|P| = ", sqrt(sum(tmpdip**2)), " (Debye)    ", "(Px,Py,Pz)= " , tmpdip(1:3) 
96   format (a8,Es14.6 ,a,a,3ES13.4)
     !     write(*,'(a)') "  ================= Dipole moment in e.a0    (0.39343 e.a0 = 1 Debye) ================"  ! or [Debye] 
     !     write(*,97) "    Px " ,"     Py ","     Pz ","   |P| " 
     !     write(*,98) "electronic charge: ", dipole_el(1:3) , sqrt(sum(dipole_el**2))
     !     write(*,98) "pseudo cores:      ", dipole_cores(1:3) , sqrt(sum(dipole_cores**2))
     !     write(*,98) "Total (cores-el.): ", dipole_cores+dipole_el , sqrt(sum((dipole_cores+dipole_el)**2))
     !97 format (20x,3a15  ,"    ==> ",a15)
     !98 format (a20,3f15.7,"    ==> ",f15.5)

  endif

  if (nproc > 1) then
     i_all=-product(shape(ele_rho))*kind(ele_rho)
     deallocate(ele_rho,stat=i_stat)
     call memocc(i_stat,i_all,'ele_rho',subname)
!!$     i_all=-product(shape(rho_buf))*kind(rho_buf)
!!$     deallocate(rho_buf,stat=i_stat)
!!$     call memocc(i_stat,i_all,'rho_buf',subname)
  else
     nullify(ele_rho)
!!$     nullify(rho_buf)
  end if

END SUBROUTINE calc_dipole
