!create the descriptors for the density and the potential
subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1i,n2i,n3i,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  use Poisson_Solver

  implicit none
  character(len=1), intent(in) :: geocode,datacode
  integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,ixc
  integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
  integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
  !local variables
  integer :: jproc

  if (datacode == 'D') then
     do jproc=0,iproc-1
        call PS_dim4allocation(geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d            !number of planes for the density
        nscatterarr(jproc,2)=n3p            !number of planes for the potential
        nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
        nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
     end do
     do jproc=iproc+1,nproc-1
        call PS_dim4allocation(geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d
        nscatterarr(jproc,2)=n3p
        nscatterarr(jproc,3)=i3s+i3xcsh-1
        nscatterarr(jproc,4)=i3xcsh
     end do
  end if

  call PS_dim4allocation(geocode,datacode,iproc,nproc,n1i,n2i,n3i,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s)
  nscatterarr(iproc,1)=n3d
  nscatterarr(iproc,2)=n3p
  nscatterarr(iproc,3)=i3s+i3xcsh-1
  nscatterarr(iproc,4)=i3xcsh

  ngatherarr(:,1)=n1i*n2i*nscatterarr(:,2)
  ngatherarr(:,2)=n1i*n2i*nscatterarr(:,3)

end subroutine createDensPotDescriptors


!partition the orbitals between processors to ensure load balancing
!the criterion will depend on GPU computation
!and/or on the sizes of the different localisation region
subroutine orbitals_communicators(iproc,nproc,lr,orbs,comms)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(out) :: comms
  !local variables
  integer :: jproc,i,nvctr_tot,j
!!$  real(kind=8), parameter :: eps_mach=1.d-12
!!$  real(kind=8) :: tt

  !calculate the number of elements to be sent to each process
  !and the array of displacements
  !cubic strategy: -the components are equally distributed among the wavefunctions
  !                -each processor has all the orbitals in transposed form
  !                -each wavefunction is equally distributed in its transposed form
  !send buffer

  !check of allocation of important arrays
  if (.not. associated(orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if

  !initialise the array
  do jproc=0,nproc-1
     comms%nvctr_par(jproc)=0 !size 0 nproc-1
  end do

  !balance the components between processors
  !in the most symmetric way
  i=1
  j=1
  loop_components: do 
     jproc=mod(i-1,nproc)
     if (.true.) then !here there is the criterion for filling a processor
        comms%nvctr_par(jproc)=comms%nvctr_par(jproc)+1
        j=j+1
     end if
     if (j > lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) exit loop_components
     i=i+1
  end do loop_components

  !check the distribution
  nvctr_tot=0
  do jproc=0,nproc-1
     nvctr_tot=nvctr_tot+comms%nvctr_par(jproc)
  end do
  if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
     write(*,*)'ERROR: partition of components incorrect'
     stop
  end if


  !if (iproc == 0) print *,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,'nvctrp',comms%nvctr_par(:)

  !calculate the dimension of the wavefunction
  !for the given processor
  !take into account max one k-point per processor
  orbs%npsidim=max((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc),&
       comms%nvctr_par(iproc)*orbs%norb)*orbs%nspinor

  if (iproc == 0) write(*,'(1x,a,i0)') &
       'Wavefunction memory occupation per orbital (Bytes): ',&
       orbs%npsidim*8

  !here the k-points should be taken into account

  do jproc=0,nproc-1
     comms%ncntd(jproc)=comms%nvctr_par(jproc)*orbs%norb_par(iproc)*orbs%nspinor
  end do
  comms%ndspld(0)=0
  do jproc=1,nproc-1
     comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
  end do
  !receive buffer
  do jproc=0,nproc-1
     comms%ncntt(jproc)=comms%nvctr_par(iproc)*orbs%norb_par(jproc)*orbs%nspinor
  end do
  comms%ndsplt(0)=0
  do jproc=1,nproc-1
     comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)
  end do

end subroutine orbitals_communicators
