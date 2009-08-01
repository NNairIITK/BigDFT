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
  character(len=*), parameter :: subname='orbitals_communicators'
  integer :: jproc,i,nvctr_tot,j,ikpts,iorbp,iorb,jorb,norb_tot,ikpt,i_stat,i_all
  integer :: ncomp_res
  integer, dimension(:,:), allocatable :: nvctr_par,norb_par !for all the components and orbitals (with k-pts)
  

  !calculate the number of elements to be sent to each process
  !and the array of displacements
  !cubic strategy: -the components are equally distributed among the wavefunctions
  !                -each processor has all the orbitals in transposed form
  !                -each wavefunction is equally distributed in its transposed form
  !                -this holds for each k-point, which regroups different processors

  !check of allocation of important arrays
  if (.not. associated(orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if

  !initialise the array
  do ikpts=1,orbs%nkptsp
     do jproc=0,nproc-1
        comms%nvctr_par(jproc,ikpts)=0 
     end do
  end do

  allocate(nvctr_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,nvctr_par,'nvctr_par',subname)
  allocate(norb_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)


  !initialise the arrays
  do ikpts=0,orbs%nkpts
     do jproc=0,nproc-1
        nvctr_par(jproc,ikpts)=0 
        norb_par(jproc,ikpts)=0 
     end do
  end do

  !balance the components between processors
  !in the most symmetric way
  !here the components are taken into account for all the k-points
  i=1
  j=1
  loop_components: do 
     jproc=mod(i-1,nproc)
     if (.true.) then !here there is the criterion for filling a processor
        nvctr_par(jproc,0)=nvctr_par(jproc,0)+1
        j=j+1
     end if
     if (j > (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nkpts) exit loop_components
     i=i+1
  end do loop_components

  ikpts=1
  ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
  do jproc=0,nproc-1
     loop_comps: do
        print *,jproc,nvctr_par(jproc,0),ncomp_res 
        if (nvctr_par(jproc,0) >= ncomp_res) then
           nvctr_par(jproc,ikpts)= ncomp_res
           ikpts=ikpts+1
           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
           nvctr_par(jproc,0)=nvctr_par(jproc,0)-ncomp_res
        else
           nvctr_par(jproc,ikpts)= nvctr_par(jproc,0)
           ncomp_res=ncomp_res-nvctr_par(jproc,0)
           nvctr_par(jproc,0)=0
        end if
        if (nvctr_par(jproc,0) == 0 ) then
           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
           exit loop_comps
        end if
     end do loop_comps
  end do

  !some checks
!!$  if (ikpts /= orbs%nkpts ) then
!!$     write(*,*)' ERROR:ikpts not correct:',ikpts,orbs%nkpts
!!$     stop
!!$  end if
  !check the distribution
  do ikpts=1,orbs%nkpts
     nvctr_tot=0
     do jproc=0,nproc-1
        nvctr_tot=nvctr_tot+nvctr_par(jproc,ikpts)
     end do
     if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
        write(*,*)'ERROR: partition of components incorrect, kpoint:',ikpts
        stop
     end if
  end do

  !calculate the same k-point distribution for the orbitals
  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=1
  ikpts=1
  do jproc=0,nproc-1
     do iorbp=1,orbs%norb_par(jproc)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (jorb == orbs%norb) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do
  !some checks
  if (orbs%norb /= 0) then
     if (ikpts /= orbs%nkpts+1 ) then
        write(*,*)' ERROR:ikpts not correct, orbitals:',ikpts,orbs%nkpts
        stop
     end if
     !check the distribution
     do ikpts=1,orbs%nkpts
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+norb_par(jproc,ikpts)
        end do
        if(norb_tot /= orbs%norb) then
           write(*,*)'ERROR: partition of components incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if

!!$  !balance the components between processors
!!$  !in the most symmetric way
!!$  i=1
!!$  j=1
!!$  loop_components: do 
!!$     jproc=mod(i-1,nproc)
!!$     if (.true.) then !here there is the criterion for filling a processor
!!$        comms%nvctr_par(jproc)=comms%nvctr_par(jproc)+1
!!$        j=j+1
!!$     end if
!!$     if (j > lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) exit loop_components
!!$     i=i+1
!!$  end do loop_components
!!$
!!$  !check the distribution
!!$  nvctr_tot=0
!!$  do jproc=0,nproc-1
!!$     nvctr_tot=nvctr_tot+comms%nvctr_par(jproc)
!!$  end do
!!$  if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
!!$     write(*,*)'ERROR: partition of components incorrect'
!!$     stop
!!$  end if


  !if (iproc == 0) print *,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,'nvctrp',comms%nvctr_par(:)

  !assign the partition of the k-points to the communication array
  do ikpts=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikpts
     do jproc=0,nproc-1
        comms%nvctr_par(jproc,ikpts)=nvctr_par(jproc,ikpt) 
     end do
  end do

  !with this distribution the orbitals and the components are ordered following k-points
  !there must be no overlap for the components
  !here we will print out the k-points components distribution, in the transposed and in the direct way

  print *,'iproc,nvctr_par,norb_par',iproc,nvctr_par(:,:),norb_par(:,:)

  do jproc=0,nproc-1
     comms%ncntd(jproc)=0
     do ikpts=1,orbs%nkpts
        comms%ncntd(jproc)=comms%ncntd(jproc)+&
             nvctr_par(jproc,ikpts)*norb_par(iproc,ikpts)*orbs%nspinor
     end do
  end do
  comms%ndspld(0)=0
  do jproc=1,nproc-1
     comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
  end do
  !receive buffer
  do jproc=0,nproc-1
     comms%ncntt(jproc)=0
     do ikpts=1,orbs%nkpts
        comms%ncntt(jproc)=comms%ncntt(jproc)+&
             nvctr_par(iproc,ikpts)*norb_par(jproc,ikpts)*orbs%nspinor
     end do
  end do
  comms%ndsplt(0)=0
  do jproc=1,nproc-1
     comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)
  end do

  print *,'iproc,comms',iproc,comms%ncntd,comms%ndspld,comms%ncntt,comms%ndsplt

!!$  do jproc=0,nproc-1
!!$     comms%ncntd(jproc)=comms%nvctr_par(jproc)*orbs%norb_par(iproc)*orbs%nspinor
!!$  end do
!!$  comms%ndspld(0)=0
!!$  do jproc=1,nproc-1
!!$     comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
!!$  end do
!!$  !receive buffer
!!$  do jproc=0,nproc-1
!!$     comms%ncntt(jproc)=comms%nvctr_par(iproc)*orbs%norb_par(jproc)*orbs%nspinor
!!$  end do
!!$  comms%ndsplt(0)=0
!!$  do jproc=1,nproc-1
!!$     comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)
!!$  end do

  i_all=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_par',subname)
  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)

  !calculate the dimension of the wavefunction
  !for the given processor
  !take into account max one k-point per processor
  orbs%npsidim=max((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc)*orbs%nspinor,&
       sum(comms%ncntt(0:nproc-1)))

  if (iproc == 0) write(*,'(1x,a,i0)') &
       'Wavefunctions memory occupation for root processor (Bytes): ',&
       orbs%npsidim*8



end subroutine orbitals_communicators
