!!****f* BigDFT/test_forces
!! NAME
!!   test_forces
!!
!! FUNCTION
!!    Runs BigDFT and test whether the forces are the 
!!    derivative of the energy.
!!    Performs the integration of the calculated forces over
!!    some random displacement and compare the result with the 
!!    difference of the energy between the final and the initial 
!!    position
!!
!! SYNOPSIS
!! WARNING
!! AUTHOR
!!    Luigi Genovese
!! COPYRIGHT
!!    Copyright (C) 2005 CEA
!! CREATION DATE
!!    09/2006
!!
!! SOURCE
!!
program test_forces
  use libBigDFT
  
  implicit none
  include 'mpif.h'
  integer, parameter :: n=31
  real(kind=8), dimension(:,:), allocatable :: rxyz,fxyz,drxyz,rxyz_old
  real(kind=8), dimension(:), allocatable :: weight
  character*20 tatonam,units
  integer, allocatable, dimension(:) :: iatype
  character*20 :: atomnames(100)
  integer :: nat,nproc,iproc,ntypes,ityp,iat,i,ierr
  real(kind=8) :: energy,energy0,FxdRx,FydRy,FzdRz,path,sumx,sumy,sumz,dx
  logical :: parallel=.true.
  real(kind=8), pointer :: psi(:,:), eval(:)
  integer, pointer :: keyg(:,:), keyv(:)
  integer :: nseg_c, nseg_f, nvctr_c, nvctr_f
  integer :: norb, norbp, n1, n2, n3,i_all,i_stat
  real*8 :: hgrid

  !Body
  !Start MPI in parallel version
  if (parallel) then
     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
     write(6,*) 'mpi started',iproc,nproc
  else
     nproc=1
     iproc=0
  endif
  
! read starting atomic positions
  open(unit=9,file='posinp',status='old')
  read(9,*) nat,units
  if (iproc.eq.0) write(6,*) 'nat=',nat
  allocate(rxyz(3,nat),stat=i_all)
  allocate(iatype(nat),stat=i_stat)
  i_all=i_all+i_stat
  allocate(fxyz(3,nat),stat=i_stat)
  i_all=i_all+i_stat
  allocate(drxyz(3,nat),stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' test_forces: problem of memory allocation'
     stop
  end if
  ntypes=0
  do iat=1,nat
     read(9,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),tatonam
     do ityp=1,ntypes
        if (tatonam.eq.atomnames(ityp)) then
           iatype(iat)=ityp
           goto 200
        endif
     end do
     ntypes=ntypes+1
     if (ntypes.gt.100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     iatype(iat)=ntypes
200  continue
     if (units.eq.'angstroem') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/.529177d0 
        end do
     else if  (units.eq.'atomic' .or. units.eq.'bohr') then
     else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
     end if
  end do
  close(9)
  do ityp=1,ntypes
     if (iproc.eq.0) write(*,*) 'atoms of type ',ityp,' are ',atomnames(ityp)
  enddo

  !calculate the starting point
  if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy0,fxyz, &
             & psi, keyg, keyv, nvctr_c, nvctr_f, nseg_c, nseg_f, norbp, norb, eval, &
             & 0, .false., .false., n1, n2, n3, hgrid, rxyz_old)
  deallocate(psi,stat=i_all)
  deallocate(eval,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(keyg,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(keyv,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' test_forces: problem of memory deallocation'
     stop
  end if

  allocate(weight(n),stat=i_all)
  if (i_all /= 0) then
     write(*,*)' test_forces: problem of memory allocation'
     stop
  end if
  !prepare the array of the correct weights of the iteration steps
  if (mod(n,2).ne.1) stop 'the number of iteration steps has to be odd'
  weight(1)=1.d0/3.d0
  weight(2)=4.d0/3.d0
  do i=3,n-2,2
     weight(i)=2.d0/3.d0
     weight(i+1)=4.d0/3.d0
  enddo
  weight(n)=1.d0/3.d0

  !then start the integration steps
  dx=1.d-2
  path=0.d0

  !calculate the displacement at each integration step
  !(use sin instead of random numbers)
  do iat=1,nat
     drxyz(1,iat)=dx*abs(dsin(real(iat,kind=8)+.2d0))
     drxyz(2,iat)=dx*abs(dsin(real(iat,kind=8)+.4d0))
     drxyz(3,iat)=dx*abs(dsin(real(iat,kind=8)+.7d0))
  end do
 
  do i=1,n


     if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz, &
             & psi, keyg, keyv, nvctr_c, nvctr_f, nseg_c, nseg_f, norbp, norb, eval, &
             & 0, .false., .false., n1, n2, n3, hgrid, rxyz_old)
     deallocate(psi,stat=i_all)
     deallocate(eval,stat=i_stat)
     i_all=i_all+i_stat
     deallocate(keyg,stat=i_stat)
     i_all=i_all+i_stat
     deallocate(keyv,stat=i_stat)
     if (i_all+i_stat /= 0) then
        write(*,*)' test_forces: problem of memory deallocation'
        stop
     end if
     


     FxdRx=0.d0
     FydRy=0.d0
     FzdRz=0.d0
     do iat=1,nat
        !integrate forces*displacement
        FxdRx=FxdRx+fxyz(1,iat)*drxyz(1,iat)
        FydRy=FydRy+fxyz(2,iat)*drxyz(2,iat)
        FzdRz=FzdRz+fxyz(3,iat)*drxyz(3,iat)

        !update atomic positions
        rxyz(1,iat)=rxyz(1,iat)+drxyz(1,iat)
        rxyz(2,iat)=rxyz(2,iat)+drxyz(2,iat)
        rxyz(3,iat)=rxyz(3,iat)+drxyz(3,iat)

     end do
     path=path-weight(i)*(FxdRx+FydRy+FzdRz)


     if (iproc==0) then 
        !verify if the system is still translational invariant by summing the forces
        sumx=0.d0
        sumy=0.d0
        sumz=0.d0
        do iat=1,nat
           sumx=sumx+fxyz(1,iat)
           sumy=sumy+fxyz(2,iat)
           sumz=sumz+fxyz(3,iat)
        end do

        write(67,'(a12,i3,2(a10,3x,e13.5),e13.5,e13.5)')'Iteration',i,'Energy ',energy,&
             'Force sums',sumx,sumy,sumz
        write(67,'(a12,i3,a24,3x,e16.8)')'Max. Iter',n,'Energy Difference',energy-energy0
        write(67,'(a39,3x,e16.8,a12,e16.8)')'Force Integral',path,'Delta',energy-energy0-path
        write(67,'(a39,3x,e16.8,a12,e16.8)')'Error on Delta',sqrt(sumx**2+sumy**2+sumz**2)
     end if

  end do

  deallocate(rxyz,stat=i_all)
  deallocate(iatype,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(fxyz,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(drxyz,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(weight,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' test_forces: problem of memory deallocation'
     stop
  end if
  
  if (parallel) call MPI_FINALIZE(ierr)


end program test_forces
