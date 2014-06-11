!>     program to convert CARTESIAN COORDINATES TO Z-MATRIX
program carint
  use internal_coordinates
  use dynamic_memory
  use yaml_output
  implicit none

  real(kind=8),parameter :: degree = 57.29578d+00
  !integer,parameter :: nmax=100
  integer :: nat, iat
  character(len=4) :: test
  !!real(kind=8),dimension(3,nmax) :: xyz, geo
  !!integer,dimension(nmax) :: na, nb, nc
  real(kind=8),dimension(:,:),allocatable :: xyz, geo
  integer,dimension(:),allocatable :: na, nb, nc
  real(kind=8) :: x, y, z
  integer :: numat, i, istat
  character(len=64) :: tt

  call f_lib_initialize()

  open(unit=99,file='posinp.xyz')
  read(99,*) nat
  xyz = f_malloc((/ 3, nat /),id='xyz')
  na = f_malloc(nat,id='na')
  nb = f_malloc(nat,id='nb')
  nc = f_malloc(nat,id='nc')
  geo = f_malloc((/ 3, nat /),id='geo')
  read(99,*) tt
  do iat=1,nat
      read(99,*) tt, xyz(1,iat), xyz(2,iat), xyz(3,iat)
  end do
  close(unit=99)


  call yaml_open_sequence('initial coordinates')
  do i=1,nat
     call yaml_sequence(advance='no')
     call yaml_open_map(flow=.true.)
     call yaml_map('atom number',i)
     call yaml_map('positions',xyz(1:3,i))
     call yaml_close_map()
  end do
  call yaml_close_sequence()

  call xyzint(xyz,nat,na,nb,nc,degree,geo)


!!  do i=1,nat
!!     write(6,*)' '
!!  end do
!!  do i=1,nat
!!     if (geo(3,i).gt.180.d0) geo(3,i)=geo(3,i)-360.d0
!!     write(6,103)na(i),geo(1,i),nb(i),geo(2,i),nc(i),geo(3,i)
!!  end do
  call yaml_open_sequence('internal coordinates')
  do i=1,nat
     call yaml_sequence(advance='no')
     call yaml_open_map(flow=.true.)
     call yaml_map('atom number',i)
     call yaml_map('references',(/na(i),nb(i),nc(i)/))
     call yaml_map('values',geo(1:3,i))
     call yaml_close_map()
  end do
  call yaml_close_sequence()
!!  write(6,*)' '
!!110     format(a4,3f15.10)
!!103     format(i5,f12.5,i5,f12.5,i5,f12.5)


  ! subtract 180 degrees
  geo(2:3,1:nat) = geo(2:3,1:nat) - 180.d0

  ! convert to rad
  geo(2:3,1:nat) = geo(2:3,1:nat) / degree
  call internal_to_cartesian(nat, geo, xyz)

  call yaml_open_sequence('final coordinates')
  do i=1,nat
     call yaml_sequence(advance='no')
     call yaml_open_map(flow=.true.)
     call yaml_map('atom number',i)
     call yaml_map('positions',xyz(1:3,i))
     call yaml_close_map()
  end do
  call yaml_close_sequence()

!!  do i=1,nat
!!     write(6,*)' '
!!  end do
!!  do i=1,nat
!!     write(6,'(3f12.5)') xyz(1:3,i)
!!  end do
!!  write(6,*)' '

  call f_free(xyz)
  call f_free(geo)
  call f_free(na)
  call f_free(nb)
  call f_free(nc)


  call f_lib_finalize()

end program carint

