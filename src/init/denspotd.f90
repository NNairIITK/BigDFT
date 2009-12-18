!!****f* BigDFT/createDensPotDescriptors
!! FUNCTION
!!   Create the descriptors for the density and the potential
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1i,n2i,n3i,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  use Poisson_Solver

  implicit none
  !Arguments
  character(len=1), intent(in) :: geocode,datacode
  integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,ixc
  integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
  integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
  !Local variables
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
!!***

!!****f* BigDFT/orbitals_communicators
!! FUNCTION
!!   Partition the orbitals between processors to ensure load balancing
!!   the criterion will depend on GPU computation
!!   and/or on the sizes of the different localisation region
!!
!! SOURCE
!!
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
  integer :: ncomp_res,iskpts,nkptsp,ierr
  logical, dimension(:), allocatable :: GPU_for_comp
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

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines
  allocate(GPU_for_comp(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,GPU_for_comp,'GPU_for_comp',subname)

  if (nproc > 1 .and. .not. GPUshare) then
     call MPI_ALLGATHER(GPUblas,1,MPI_LOGICAL,GPU_for_comp(0),1,MPI_LOGICAL,&
          MPI_COMM_WORLD,ierr)
  else
     GPU_for_comp(0)=GPUblas
  end if

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


  i_all=-product(shape(GPU_for_comp))*kind(GPU_for_comp)
  deallocate(GPU_for_comp,stat=i_stat)
  call memocc(i_stat,i_all,'GPU_for_comp',subname)


  ikpts=1
  ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
  do jproc=0,nproc-1
     loop_comps: do
        if (nvctr_par(jproc,0) >= ncomp_res) then
           nvctr_par(jproc,ikpts)= ncomp_res
           ikpts=ikpts+1
           nvctr_par(jproc,0)=nvctr_par(jproc,0)-ncomp_res
           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        else
           nvctr_par(jproc,ikpts)= nvctr_par(jproc,0)
           ncomp_res=ncomp_res-nvctr_par(jproc,0)
           nvctr_par(jproc,0)=0
           exit loop_comps
        end if
        if (nvctr_par(jproc,0) == 0 ) then
           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
           exit loop_comps
        end if

     end do loop_comps
  end do

  !some checks
!!!  if (ikpts /= orbs%nkpts ) then
!!!     write(*,*)' ERROR:ikpts not correct:',ikpts,orbs%nkpts
!!!     stop
!!!  end if
  !check the distribution
  do ikpts=1,orbs%nkpts
     !print *,'iproc,cpts:',lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nvctr_par(:,ikpts)
     nvctr_tot=0
     do jproc=0,nproc-1
        nvctr_tot=nvctr_tot+nvctr_par(jproc,ikpts)
     end do
     if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
        write(*,*)'ERROR: partition of components incorrect, kpoint:',ikpts
        stop
     end if
  end do

  !calculate the number of k-points treated by each processor in the component distribution
  nkptsp=0
  iskpts=-1
  do ikpts=1,orbs%nkpts
     if (nvctr_par(iproc,ikpts) /= 0) then
        nkptsp=nkptsp+1
        if (iskpts == -1) then
           iskpts=ikpts-1
        end if
     end if
  end do
  
  orbs%nkptsp=nkptsp
  orbs%iskpts=iskpts

  !this function which associates a given k-point to a processor 
  !the association is chosen such that each k-point is associated to only
  !one processor
  !if two processors treat the same k-point the processor which highest rank is chosen
  do ikpts=1,orbs%nkpts
     loop_jproc: do jproc=nproc-1,0,-1
        if (nvctr_par(jproc,ikpts) /= 0) then
           orbs%ikptproc(ikpts)=jproc
           exit loop_jproc
        end if
     end do loop_jproc
  end do
  
  !print *,'check',orbs%ikptproc(:)

  !calculate the same k-point distribution for the orbitals
  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=1
  ikpts=1

  !print *,'here',orbs%norb_par(:)
  do jproc=0,nproc-1
     do iorbp=1,orbs%norb_par(jproc)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (mod(jorb,orbs%norb)==0) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do
  !some checks
  if (orbs%norb /= 0) then
!!!     if (ikpts /= orbs%nkpts) then
!!!        write(*,*)' ERROR:ikpts not correct, orbitals:',ikpts,orbs%nkpts
!!!        stop
!!!     end if
     !check the distribution
     do ikpts=1,orbs%nkpts
        !print *,'partition',ikpts,orbs%nkpts,'ikpts',norb_par(:,ikpts)
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+norb_par(jproc,ikpts)
        end do
        if(norb_tot /= orbs%norb) then
           write(*,*)'ERROR: partition of orbitals incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if

  !allocate communication arrays
  allocate(comms%nvctr_par(0:nproc-1,orbs%nkptsp+ndebug),stat=i_stat)
  call memocc(i_stat,comms%nvctr_par,'nvctr_par',subname)

  allocate(comms%ncntd(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ncntd,'ncntd',subname)

  allocate(comms%ncntt(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ncntt,'ncntt',subname)
  allocate(comms%ndspld(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ndspld,'ndspld',subname)
  allocate(comms%ndsplt(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ndsplt,'ndsplt',subname)


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

  !print *,'iproc,nvctr_par,norb_par',iproc,nvctr_par(:,:),norb_par(:,:)

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

  !print *,'iproc,comms',iproc,comms%ncntd,comms%ndspld,comms%ncntt,comms%ndsplt

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
!!***

!!subroutine print_distribution_schemes(iproc,nproc,nkpts,norb_par,nvctr_par)
!!  use module_base
!!  implicit none
!!  integer, intent(in) :: iproc,nproc,nkpts
!!  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par,nvctr_par
!!  !local variables
!!  
!!end subroutine print_distribution_schemes
