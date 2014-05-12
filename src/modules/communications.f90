!> @file
!!  File defining the structures and the routines for the communication between processes
!! @author
!!    Copyright (C) 2014-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining routines related to communications (mainly transpositions)
module communications

  use module_base
  use communications_base, only: comms_linear, comms_cubic

  implicit none

  private

  public :: transpose_localized
  public :: untranspose_localized
  public :: transpose_switch_psir
  public :: transpose_communicate_psir
  public :: transpose_unswitch_psirt
  public :: start_onesided_communication
  public :: synchronize_onesided_communication
  public :: communicate_locreg_descriptors_basics
  public :: communicate_locreg_descriptors_keys
  public :: transpose_v
  public :: untranspose_v


  contains

    subroutine transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: npsidim_orbs
      type(orbitals_Data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psi
      real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
      type(local_zone_descriptors),intent(in),optional :: lzd
      
      ! Local variables
      integer :: i_tot, i_c, i_f, iorb, iiorb, ilr, i, ind, istat, iall, m, ind7, i7
      real(kind=8),dimension(:),allocatable :: psi_c, psi_f
      character(len=*),parameter :: subname='transpose_switch_psi'
    
      allocate(psi_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psi_c, 'psi_c', subname)
      allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psi_f, 'psi_f', subname)
    
    
      if(present(lzd)) then
      ! split up psi into coarse and fine part
    
      
          i_tot=0
          i_c=0
          i_f=0
    
          do iorb=1,orbs%norbp
             iiorb=orbs%isorb+iorb
             ilr=orbs%inwhichlocreg(iiorb)
    
             call vcopy(lzd%llr(ilr)%wfd%nvctr_c,psi(i_tot+1),1,psi_c(i_c+1),1)
    
             i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
             i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
    
             call vcopy(7*lzd%llr(ilr)%wfd%nvctr_f,psi(i_tot+1),1,psi_f(i_f+1),1)
    
             i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
             i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f
    
          end do
        
    
      else
          ! only coarse part is used...
          call vcopy(collcom%ndimpsi_c, psi(1), 1, psi_c(1), 1)
      end if
    
      ! coarse part
    
      !$omp parallel default(private) &
      !$omp shared(collcom, psi, psiwork_c, psiwork_f, lzd, psi_c,psi_f,m)
    
      m = mod(collcom%ndimpsi_c,7)
      if(m/=0) then
         do i=1,m
            ind = collcom%isendbuf_c(i)
            psiwork_c(ind) = psi_c(i)
         end do
      end if
      !$omp do
      do i = m+1,collcom%ndimpsi_c,7
         psiwork_c(collcom%isendbuf_c(i+0)) = psi_c(i+0)
         psiwork_c(collcom%isendbuf_c(i+1)) = psi_c(i+1)
         psiwork_c(collcom%isendbuf_c(i+2)) = psi_c(i+2)
         psiwork_c(collcom%isendbuf_c(i+3)) = psi_c(i+3)
         psiwork_c(collcom%isendbuf_c(i+4)) = psi_c(i+4)
         psiwork_c(collcom%isendbuf_c(i+5)) = psi_c(i+5)
         psiwork_c(collcom%isendbuf_c(i+6)) = psi_c(i+6)
      end do
      !$omp end do
     
      ! fine part
    
      !$omp do
       do i=1,collcom%ndimpsi_f
          ind=collcom%isendbuf_f(i)
          i7=7*i
          ind7=7*ind
          psiwork_f(ind7-6)=psi_f(i7-6)
          psiwork_f(ind7-5)=psi_f(i7-5)
          psiwork_f(ind7-4)=psi_f(i7-4)
          psiwork_f(ind7-3)=psi_f(i7-3)
          psiwork_f(ind7-2)=psi_f(i7-2)
          psiwork_f(ind7-1)=psi_f(i7-1)
          psiwork_f(ind7-0)=psi_f(i7-0)
      end do
      !$omp end do
      !$omp end parallel
    
    
      iall=-product(shape(psi_c))*kind(psi_c)
      deallocate(psi_c, stat=istat)
      call memocc(istat, iall, 'psi_c', subname)
      iall=-product(shape(psi_f))*kind(psi_f)
      deallocate(psi_f, stat=istat)
      call memocc(istat, iall, 'psi_f', subname)
      
    end subroutine transpose_switch_psi



    subroutine transpose_communicate_psi(iproc, nproc, collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimpsi_c),intent(in) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(in) :: psiwork_f
      !real(kind=8),dimension(sum(collcom%nrecvcounts_c)),intent(out) :: psitwork_c
      !real(kind=8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out) :: psitwork_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psitwork_f
      
      ! Local variables
      integer :: ierr, istat, iall
      !!integer :: iisend, iirecv, ist, ist_c, ist_f, jproc
      real(kind=8),dimension(:),allocatable :: psiwork, psitwork
      integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      character(len=*),parameter :: subname='transpose_communicate_psi'
    
      !call mpi_comm_size(bigdft_mpi%mpi_comm, nproc, ierr)
      !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)
    
      allocate(psiwork(collcom%ndimpsi_c+7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psiwork, 'psiwork', subname)
      allocate(psitwork(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psitwork, 'psitwork', subname)
      allocate(nsendcounts(0:nproc-1), stat=istat)
      call memocc(istat, nsendcounts, 'nsendcounts', subname)
      allocate(nsenddspls(0:nproc-1), stat=istat)
      call memocc(istat, nsenddspls, 'nsenddspls', subname)
      allocate(nrecvcounts(0:nproc-1), stat=istat)
      call memocc(istat, nrecvcounts, 'nrecvcounts', subname)
      allocate(nrecvdspls(0:nproc-1), stat=istat)
      call memocc(istat, nrecvdspls, 'nrecvdspls', subname)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!iisend=0
      !!iirecv=0
      !!do jproc=0,nproc-1
      !!    if(collcom%nsendcounts_c(jproc)>0) call vcopy(collcom%nsendcounts_c(jproc), psiwork_c(ist_c), 1, psiwork(ist), 1)
      !!    ist_c=ist_c+collcom%nsendcounts_c(jproc)
      !!    ist=ist+collcom%nsendcounts_c(jproc)
      !!    if(collcom%nsendcounts_f(jproc)>0) call vcopy(7*collcom%nsendcounts_f(jproc), psiwork_f(ist_f), 1, psiwork(ist), 1)
      !!    ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
      !!    ist=ist+7*collcom%nsendcounts_f(jproc)
      !!    nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
      !!    nsenddspls(jproc)=iisend
      !!    nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
      !!    nrecvdspls(jproc)=iirecv
      !!    iisend=iisend+nsendcounts(jproc)
      !!    iirecv=iirecv+nrecvcounts(jproc)
      !!end do
    
      !write(*,'(a,i4,4x,100i8)') 'iproc, nsendcounts', iproc, nsendcounts
      !write(*,'(a,i4,4x,100i8)') 'iproc, nsenddspls', iproc, nsenddspls
      !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvcounts', iproc, nrecvcounts
      !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvdspls', iproc, nrecvdspls
      
      ! coarse part
      call mpi_alltoallv(psiwork_c, collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, psitwork_c, &
           collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      
      ! fine part
      call mpi_alltoallv(psiwork_f, 7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, psitwork_f, &
           7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      !!call mpi_alltoallv(psiwork, nsendcounts, nsenddspls, mpi_double_precision, psitwork, &
      !!     nrecvcounts, nrecvdspls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!do jproc=0,nproc-1
      !!    if(collcom%nrecvcounts_c(jproc)>0) call vcopy(collcom%nrecvcounts_c(jproc), psitwork(ist), 1, psitwork_c(ist_c), 1)
      !!    ist_c=ist_c+collcom%nrecvcounts_c(jproc)
      !!    ist=ist+collcom%nrecvcounts_c(jproc)
      !!    if(collcom%nrecvcounts_f(jproc)>0) call vcopy(7*collcom%nrecvcounts_f(jproc), psitwork(ist), 1, psitwork_f(ist_f), 1)
      !!    ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
      !!    ist=ist+7*collcom%nrecvcounts_f(jproc)
      !!end do
    
      iall=-product(shape(psiwork))*kind(psiwork)
      deallocate(psiwork, stat=istat)
      call memocc(istat, iall, 'psiwork', subname)
      iall=-product(shape(psitwork))*kind(psitwork)
      deallocate(psitwork, stat=istat)
      call memocc(istat, iall, 'psitwork', subname)
      iall=-product(shape(nsendcounts))*kind(nsendcounts)
      deallocate(nsendcounts, stat=istat)
      call memocc(istat, iall, 'nsendcounts', subname)
      iall=-product(shape(nsenddspls))*kind(nsenddspls)
      deallocate(nsenddspls, stat=istat)
      call memocc(istat, iall, 'nsenddspls', subname)
      iall=-product(shape(nrecvcounts))*kind(nrecvcounts)
      deallocate(nrecvcounts, stat=istat)
      call memocc(istat, iall, 'nrecvcounts', subname)
      iall=-product(shape(nrecvdspls))*kind(nrecvdspls)
      deallocate(nrecvdspls, stat=istat)
      call memocc(istat, iall, 'nrecvdspls', subname)
    
    
    end subroutine transpose_communicate_psi



    subroutine transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
      
      ! Local variables
      integer :: i,ind,sum_c,sum_f,m,i7,ind7
    
      sum_c = sum(collcom%nrecvcounts_c)
      sum_f = sum(collcom%nrecvcounts_f)
    
      !$omp parallel private(i,ind,i7,ind7) &
      !$omp shared(psit_c,psit_f, psitwork_c, psitwork_f,collcom,sum_c,sum_f,m)
    
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i = 1,m
          ind=collcom%iextract_c(i)
          psit_c(ind)=psitwork_c(i)
        end do
      end if
    
      ! coarse part
    
      !$omp do
      do i=m+1, sum_c,7
          psit_c(collcom%iextract_c(i+0))=psitwork_c(i+0)
          psit_c(collcom%iextract_c(i+1))=psitwork_c(i+1)
          psit_c(collcom%iextract_c(i+2))=psitwork_c(i+2)
          psit_c(collcom%iextract_c(i+3))=psitwork_c(i+3)
          psit_c(collcom%iextract_c(i+4))=psitwork_c(i+4)
          psit_c(collcom%iextract_c(i+5))=psitwork_c(i+5)
          psit_c(collcom%iextract_c(i+6))=psitwork_c(i+6)
      end do
      !$omp end do
    
      ! fine part
    
      !$omp do
      do i=1,sum_f
          ind=collcom%iextract_f(i)
          i7=7*i
          ind7=7*ind
          psit_f(ind7-6)=psitwork_f(i7-6)
          psit_f(ind7-5)=psitwork_f(i7-5)
          psit_f(ind7-4)=psitwork_f(i7-4)
          psit_f(ind7-3)=psitwork_f(i7-3)
          psit_f(ind7-2)=psitwork_f(i7-2)
          psit_f(ind7-1)=psitwork_f(i7-1)
          psit_f(ind7-0)=psitwork_f(i7-0)
      end do
      !$omp end do
      
      !$omp end parallel
    
    end subroutine transpose_unswitch_psit



    subroutine transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psitwork_f
      
      ! Local variables
      integer :: i, ind, sum_c,sum_f,m, i7, ind7
    
      sum_c = sum(collcom%nrecvcounts_c)
      sum_f = sum(collcom%nrecvcounts_f)
    
      !$omp parallel default(private) &
      !$omp shared(collcom, psit_c,psit_f, psitwork_c, psitwork_f,sum_c,sum_f,m)
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i=1,m
           ind = collcom%iexpand_c(i)
           psitwork_c(ind) = psit_c(i)
        end do
      end if
    
      ! coarse part
    
      !$omp do
      do i=m+1,sum_c,7
          psitwork_c(collcom%iexpand_c(i+0))=psit_c(i+0)
          psitwork_c(collcom%iexpand_c(i+1))=psit_c(i+1)
          psitwork_c(collcom%iexpand_c(i+2))=psit_c(i+2)
          psitwork_c(collcom%iexpand_c(i+3))=psit_c(i+3)
          psitwork_c(collcom%iexpand_c(i+4))=psit_c(i+4)
          psitwork_c(collcom%iexpand_c(i+5))=psit_c(i+5)
          psitwork_c(collcom%iexpand_c(i+6))=psit_c(i+6)
      end do
      !$omp end do
    
      ! fine part
    
      !$omp do
      do i=1,sum_f
          i7=7*i
          ind7=7*ind
          ind=collcom%iexpand_f(i)
          psitwork_f(ind7-6)=psit_f(i7-6)
          psitwork_f(ind7-5)=psit_f(i7-5)
          psitwork_f(ind7-4)=psit_f(i7-4)
          psitwork_f(ind7-3)=psit_f(i7-3)
          psitwork_f(ind7-2)=psit_f(i7-2)
          psitwork_f(ind7-1)=psit_f(i7-1)
          psitwork_f(ind7-0)=psit_f(i7-0)
      end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_switch_psit



    subroutine transpose_communicate_psit(iproc, nproc, collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
      
      ! Local variables
      integer :: ierr
      !!integer :: iall, ist, ist_c, ist_f, jproc, iisend, iirecv, istat
      !!real(kind=8),dimension(:),allocatable :: psiwork, psitwork
      !!integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      !!character(len=*),parameter :: subname='transpose_communicate_psit'
    
      !call mpi_comm_size(bigdft_mpi%mpi_comm, nproc, ierr)
      !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)
    
      !!allocate(psiwork(collcom%ndimpsi_c+7*collcom%ndimpsi_f), stat=istat)
      !!call memocc(istat, psiwork, 'psiwork', subname)
      !!allocate(psitwork(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f)), stat=istat)
      !!call memocc(istat, psitwork, 'psitwork', subname)
      !!allocate(nsendcounts(0:nproc-1), stat=istat)
      !!call memocc(istat, nsendcounts, 'nsendcounts', subname)
      !!allocate(nsenddspls(0:nproc-1), stat=istat)
      !!call memocc(istat, nsenddspls, 'nsenddspls', subname)
      !!allocate(nrecvcounts(0:nproc-1), stat=istat)
      !!call memocc(istat, nrecvcounts, 'nrecvcounts', subname)
      !!allocate(nrecvdspls(0:nproc-1), stat=istat)
      !!call memocc(istat, nrecvdspls, 'nrecvdspls', subname)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!iisend=0
      !!iirecv=0
      !!do jproc=0,nproc-1
      !!    if(collcom%nrecvcounts_c(jproc)>0) call vcopy(collcom%nrecvcounts_c(jproc), psitwork_c(ist_c), 1, psitwork(ist), 1)
      !!    ist_c=ist_c+collcom%nrecvcounts_c(jproc)
      !!    ist=ist+collcom%nrecvcounts_c(jproc)
      !!    if(collcom%nrecvcounts_f(jproc)>0) call vcopy(7*collcom%nrecvcounts_f(jproc), psitwork_f(ist_f), 1, psitwork(ist), 1)
      !!    ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
      !!    ist=ist+7*collcom%nrecvcounts_f(jproc)
      !!    nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
      !!    nsenddspls(jproc)=iisend
      !!    nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
      !!    nrecvdspls(jproc)=iirecv
      !!    iisend=iisend+nsendcounts(jproc)
      !!    iirecv=iirecv+nrecvcounts(jproc)
      !!end do
    
    
      ! coarse part
       call mpi_alltoallv(psitwork_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, psiwork_c, &
            collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      ! fine part
       call mpi_alltoallv(psitwork_f, 7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, psiwork_f, &
            7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      !!call mpi_alltoallv(psitwork, nrecvcounts, nrecvdspls, mpi_double_precision, psiwork, &
      !!     nsendcounts, nsenddspls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!do jproc=0,nproc-1
      !!    if(collcom%nsendcounts_c(jproc)>0) call vcopy(collcom%nsendcounts_c(jproc), psiwork(ist), 1, psiwork_c(ist_c), 1)
      !!    ist_c=ist_c+collcom%nsendcounts_c(jproc)
      !!    ist=ist+collcom%nsendcounts_c(jproc)
      !!    if(collcom%nsendcounts_f(jproc)>0) call vcopy(7*collcom%nsendcounts_f(jproc), psiwork(ist), 1, psiwork_f(ist_f), 1)
      !!    ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
      !!    ist=ist+7*collcom%nsendcounts_f(jproc)
      !!end do
    
      !!iall=-product(shape(psiwork))*kind(psiwork)
      !!deallocate(psiwork, stat=istat)
      !!call memocc(istat, iall, 'psiwork', subname)
      !!iall=-product(shape(psitwork))*kind(psitwork)
      !!deallocate(psitwork, stat=istat)
      !!call memocc(istat, iall, 'psitwork', subname)
      !!iall=-product(shape(nsendcounts))*kind(nsendcounts)
      !!deallocate(nsendcounts, stat=istat)
      !!call memocc(istat, iall, 'nsendcounts', subname)
      !!iall=-product(shape(nsenddspls))*kind(nsenddspls)
      !!deallocate(nsenddspls, stat=istat)
      !!call memocc(istat, iall, 'nsenddspls', subname)
      !!iall=-product(shape(nrecvcounts))*kind(nrecvcounts)
      !!deallocate(nrecvcounts, stat=istat)
      !!call memocc(istat, iall, 'nrecvcounts', subname)
      !!iall=-product(shape(nrecvdspls))*kind(nrecvdspls)
      !!deallocate(nrecvdspls, stat=istat)
      !!call memocc(istat, iall, 'nrecvdspls', subname)
    
    end subroutine transpose_communicate_psit



    subroutine transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
      use module_base
      use module_types
      implicit none
      
      ! Caling arguments
      integer, intent(in) :: npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimpsi_c),intent(in) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(in) :: psiwork_f
      real(kind=8),dimension(npsidim_orbs),intent(out) :: psi
      type(local_zone_descriptors),intent(in),optional :: lzd
      
      ! Local variables
      integer :: i, ind, iorb, iiorb, ilr, i_tot, i_c, i_f, istat, iall, m, i7, ind7
      real(kind=8),dimension(:),allocatable :: psi_c, psi_f
      character(len=*),parameter :: subname='transpose_unswitch_psi'
      
      
      allocate(psi_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psi_c, 'psi_c', subname)
      allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psi_f, 'psi_f', subname)
      
      !$omp parallel default(private) &
      !$omp shared(collcom, psiwork_c, psi_c,psi_f,psiwork_f,m)
    
      m = mod(collcom%ndimpsi_c,7)
    
      if(m/=0) then
        do i = 1,m
         ind=collcom%irecvbuf_c(i)
         psi_c(ind)=psiwork_c(i) 
        end do
      end if
    
      ! coarse part
    
      !$omp do
        do i=m+1,collcom%ndimpsi_c,7
            psi_c(collcom%irecvbuf_c(i+0))=psiwork_c(i+0)
            psi_c(collcom%irecvbuf_c(i+1))=psiwork_c(i+1)
            psi_c(collcom%irecvbuf_c(i+2))=psiwork_c(i+2)
            psi_c(collcom%irecvbuf_c(i+3))=psiwork_c(i+3)
            psi_c(collcom%irecvbuf_c(i+4))=psiwork_c(i+4)
            psi_c(collcom%irecvbuf_c(i+5))=psiwork_c(i+5)
            psi_c(collcom%irecvbuf_c(i+6))=psiwork_c(i+6)
        end do
      !$omp end do
      
      ! fine part
     
      !$omp do
       do i=1,collcom%ndimpsi_f
            ind=collcom%irecvbuf_f(i)
            i7=7*i
            ind7=7*ind
            psi_f(ind7-6)=psiwork_f(i7-6)
            psi_f(ind7-5)=psiwork_f(i7-5)
            psi_f(ind7-4)=psiwork_f(i7-4)
            psi_f(ind7-3)=psiwork_f(i7-3)
            psi_f(ind7-2)=psiwork_f(i7-2)
            psi_f(ind7-1)=psiwork_f(i7-1)
            psi_f(ind7-0)=psiwork_f(i7-0)
        end do
      !$omp end do
      !$omp end parallel
    
        if(present(lzd)) then
            ! glue together coarse and fine part
    
            i_tot=0
            i_c=0
            i_f=0
            do iorb=1,orbs%norbp
                iiorb=orbs%isorb+iorb
                ilr=orbs%inwhichlocreg(iiorb)
    
                call vcopy(lzd%llr(ilr)%wfd%nvctr_c,psi_c(i_c+1),1,psi(i_tot+1),1)
    
                i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
                i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
                
                call vcopy(7*lzd%llr(ilr)%wfd%nvctr_f,psi_f(i_f+1),1,psi(i_tot+1),1)
    
    
                i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
                i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f
    
    
            end do
        !!$omp end parallel 
    
        else
            call vcopy(collcom%ndimpsi_c, psi_c(1), 1, psi(1), 1)
        end if
      
      iall=-product(shape(psi_c))*kind(psi_c)
      deallocate(psi_c, stat=istat)
      call memocc(istat, iall, 'psi_c', subname)
      iall=-product(shape(psi_f))*kind(psi_f)
      deallocate(psi_f, stat=istat)
      call memocc(istat, iall, 'psi_f', subname)
    
    end subroutine transpose_unswitch_psi



    subroutine transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psi, psit_c, psit_f, lzd)
      use module_base
      use module_types
      !use module_interfaces, except_this_one => transpose_localized
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psi
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
      type(local_zone_descriptors),optional,intent(in) :: lzd
      
      ! Local variables
      real(kind=8),dimension(:),allocatable :: psiwork_c, psiwork_f, psitwork_c, psitwork_f
      integer :: istat, iall
      character(len=*),parameter :: subname='transpose_localized'
      
      allocate(psiwork_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psiwork_c, 'psiwork_c', subname)
      allocate(psiwork_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psiwork_f, 'psiwork_f', subname)
      allocate(psitwork_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psitwork_c, 'psitwork_c', subname)
      allocate(psitwork_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psitwork_f, 'psitwork_f', subname)
      
      call timing(iproc,'Un-TransSwitch','ON')
      if(present(lzd)) then
          call transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
      else
          call transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f)
      end if
      call timing(iproc,'Un-TransSwitch','OF')
    
      call timing(iproc,'Un-TransComm  ','ON')
      if(nproc>1) then
          call transpose_communicate_psi(iproc, nproc, collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
      else
          psitwork_c=psiwork_c
          psitwork_f=psiwork_f
      end if
      call timing(iproc,'Un-TransComm  ','OF')
    
      call timing(iproc,'Un-TransSwitch','ON')
      call transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
      call timing(iproc,'Un-TransSwitch','OF')
      
      iall=-product(shape(psiwork_c))*kind(psiwork_c)
      deallocate(psiwork_c, stat=istat)
      call memocc(istat, iall, 'psiwork_c', subname)
      iall=-product(shape(psiwork_f))*kind(psiwork_f)
      deallocate(psiwork_f, stat=istat)
      call memocc(istat, iall, 'psiwork_f', subname)
      iall=-product(shape(psitwork_c))*kind(psitwork_c)
      deallocate(psitwork_c, stat=istat)
      call memocc(istat, iall, 'psitwork_c', subname)
      iall=-product(shape(psitwork_f))*kind(psitwork_f)
      deallocate(psitwork_f, stat=istat)
      call memocc(istat, iall, 'psitwork_f', subname)
      
    end subroutine transpose_localized



    subroutine untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, psi, lzd)
      use module_base
      use module_types
      !use module_interfaces, except_this_one => untranspose_localized
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
      real(kind=8),dimension(npsidim_orbs),intent(out) :: psi
      type(local_zone_descriptors),optional,intent(in) :: lzd
      
      ! Local variables
      real(kind=8),dimension(:),allocatable :: psiwork_c, psiwork_f, psitwork_c, psitwork_f
      integer :: istat, iall
      character(len=*),parameter :: subname='untranspose_localized'
      
      allocate(psiwork_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psiwork_c, 'psiwork_c', subname)
      allocate(psiwork_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psiwork_f, 'psiwork_f', subname)
      allocate(psitwork_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psitwork_c, 'psitwork_c', subname)
      allocate(psitwork_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psitwork_f, 'psitwork_f', subname)
    
      call timing(iproc,'Un-TransSwitch','ON')
      call transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
      call timing(iproc,'Un-TransSwitch','OF')
    
      call timing(iproc,'Un-TransComm  ','ON')
      if(nproc>1) then
          call transpose_communicate_psit(iproc, nproc, collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
      else
          psiwork_c=psitwork_c
          psiwork_f=psitwork_f
      end if
      call timing(iproc,'Un-TransComm  ','OF')
    
      call timing(iproc,'Un-TransSwitch','ON')
      if(present(lzd)) then
          call transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
      else
          call transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi)
      end if
      call timing(iproc,'Un-TransSwitch','OF')
      
      iall=-product(shape(psiwork_c))*kind(psiwork_c)
      deallocate(psiwork_c, stat=istat)
      call memocc(istat, iall, 'psiwork_c', subname)
      iall=-product(shape(psiwork_f))*kind(psiwork_f)
      deallocate(psiwork_f, stat=istat)
      call memocc(istat, iall, 'psiwork_f', subname)
      iall=-product(shape(psitwork_c))*kind(psitwork_c)
      deallocate(psitwork_c, stat=istat)
      call memocc(istat, iall, 'psitwork_c', subname)
      iall=-product(shape(psitwork_f))*kind(psitwork_f)
      deallocate(psitwork_f, stat=istat)
      call memocc(istat, iall, 'psitwork_f', subname)
      
    end subroutine untranspose_localized






    subroutine transpose_switch_psir(collcom_sr, psir, psirwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psir
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork
    
      ! Local variables
      integer :: i, m, ind
    
      !$omp parallel default(private) &
      !$omp shared(collcom_sr, psir, psirwork, m)
    
      m = mod(collcom_sr%ndimpsi_c,7)
      if(m/=0) then
          do i=1,m
              ind = collcom_sr%isendbuf_c(i)
              psirwork(ind) = psir(i)
          end do
      end if
      !$omp do
      do i = m+1,collcom_sr%ndimpsi_c,7
         psirwork(collcom_sr%isendbuf_c(i+0)) = psir(i+0)
         psirwork(collcom_sr%isendbuf_c(i+1)) = psir(i+1)
         psirwork(collcom_sr%isendbuf_c(i+2)) = psir(i+2)
         psirwork(collcom_sr%isendbuf_c(i+3)) = psir(i+3)
         psirwork(collcom_sr%isendbuf_c(i+4)) = psir(i+4)
         psirwork(collcom_sr%isendbuf_c(i+5)) = psir(i+5)
         psirwork(collcom_sr%isendbuf_c(i+6)) = psir(i+6)
      end do
      !$omp end do
      !$omp end parallel
    
    
    end subroutine transpose_switch_psir
    
    subroutine transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork
    
      ! Local variables
      integer :: ierr
    
    
      if (nproc>1) then
          call mpi_alltoallv(psirwork, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, psirtwork, &
               collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
          call vcopy(collcom_sr%ndimpsi_c, psirwork(1), 1, psirtwork(1), 1)
      end if
    
    
    end subroutine transpose_communicate_psir
    
    subroutine transpose_unswitch_psirt(collcom_sr, psirtwork, psirt)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirt
    
      ! Local variables
      integer :: i, ind, sum_c, m
    
      sum_c = sum(collcom_sr%nrecvcounts_c)
    
      !$omp parallel private(i,ind) &
      !$omp shared(psirt, psirtwork, collcom_sr, sum_c, m)
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i = 1,m
          ind=collcom_sr%iextract_c(i)
          psirt(ind)=psirtwork(i)
        end do
      end if
    
      !$omp do
      do i=m+1, sum_c,7
          psirt(collcom_sr%iextract_c(i+0))=psirtwork(i+0)
          psirt(collcom_sr%iextract_c(i+1))=psirtwork(i+1)
          psirt(collcom_sr%iextract_c(i+2))=psirtwork(i+2)
          psirt(collcom_sr%iextract_c(i+3))=psirtwork(i+3)
          psirt(collcom_sr%iextract_c(i+4))=psirtwork(i+4)
          psirt(collcom_sr%iextract_c(i+5))=psirtwork(i+5)
          psirt(collcom_sr%iextract_c(i+6))=psirtwork(i+6)
      end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_unswitch_psirt
    
    
    



 
    subroutine start_onesided_communication(iproc, nproc, nsendbuf, sendbuf, nrecvbuf, recvbuf, comm, lzd)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer, intent(in):: iproc, nproc, nsendbuf, nrecvbuf
      real(kind=8), dimension(nsendbuf), intent(in):: sendbuf
      real(kind=8), dimension(nrecvbuf), intent(out):: recvbuf
      type(p2pComms), intent(inout):: comm
      type(local_zone_descriptors), intent(in) :: lzd
      
      ! Local variables
      !character(len=*), parameter :: subname='start_onesided_communication'
      integer :: jproc, joverlap, mpisource, istsource, mpidest, istdest, ierr, nit
      integer :: ioffset_send, ist, i2, i3, ist2, ist3, info, nsize, size_of_double
    
    
      call timing(iproc, 'Pot_comm start', 'ON')
    
      if(.not.comm%communication_complete) stop 'ERROR: there is already a p2p communication going on...'
    
      nproc_if: if (nproc>1) then
    
          ! Allocate MPI memory window
          call mpi_type_size(mpi_double_precision, size_of_double, ierr)
          call mpi_info_create(info, ierr)
          call mpi_info_set(info, "no_locks", "true", ierr)
          call mpi_win_create(sendbuf(1), int(nsendbuf*size_of_double,kind=mpi_address_kind), size_of_double, &
               info, bigdft_mpi%mpi_comm, comm%window, ierr)
          call mpi_info_free(info, ierr)
    
          call mpi_win_fence(mpi_mode_noprecede, comm%window, ierr)
          
          do jproc=0,nproc-1
              do joverlap=1,comm%noverlaps(jproc)
                  mpisource=comm%comarr(1,joverlap,jproc)
                  istsource=comm%comarr(2,joverlap,jproc)
                  mpidest=comm%comarr(3,joverlap,jproc)
                  istdest=comm%comarr(4,joverlap,jproc)
                  nit=comm%comarr(5,joverlap,jproc)
                  ioffset_send=comm%comarr(6,joverlap,jproc)
                  call mpi_type_create_hvector(nit, 1, int(size_of_double*ioffset_send,kind=mpi_address_kind), &
                       comm%mpi_datatypes(0,jproc), comm%mpi_datatypes(joverlap,jproc), ierr)
                  call mpi_type_commit(comm%mpi_datatypes(joverlap,jproc), ierr)
                  if (iproc==mpidest) then
                      call mpi_type_size(comm%mpi_datatypes(joverlap,jproc), nsize, ierr)
                      nsize=nsize/size_of_double
                      if(nsize>0) then
                          call mpi_get(recvbuf(istdest), nsize, &
                               mpi_double_precision, mpisource, int((istsource-1),kind=mpi_address_kind), &
                               1, comm%mpi_datatypes(joverlap,jproc), comm%window, ierr)
                      end if
                  end if
              end do
          end do
    
      else nproc_if
    
          ist=1
          do i3=comm%ise(5,iproc),comm%ise(6,iproc)
              ist3=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
              do i2=comm%ise(3,iproc),comm%ise(4,iproc)
                  ist2=(i2-1)*lzd%glr%d%n1i
                  call vcopy(comm%ise(2,iproc)-comm%ise(1,iproc)+1, sendbuf(ist3+ist2+1), 1, recvbuf(ist), 1)
                  ist=ist+comm%ise(2,iproc)-comm%ise(1,iproc)+1
              end do
          end do
    
      end if nproc_if
      
      
      ! Flag indicating whether the communication is complete or not
      if(nproc>1) then
          comm%communication_complete=.false.
      else
          comm%communication_complete=.true.
      end if
    
      call timing(iproc, 'Pot_comm start', 'OF')
    
    end subroutine start_onesided_communication
    
    
    subroutine synchronize_onesided_communication(iproc, nproc, comm)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc
      type(p2pComms),intent(inout):: comm
      
      ! Local variables
      integer:: ierr, jproc, joverlap
      
      
      if(.not.comm%communication_complete) then
          call mpi_win_fence(0, comm%window, ierr)
          do jproc=0,nproc-1
              do joverlap=1,comm%noverlaps(jproc)
                  call mpi_type_free(comm%mpi_datatypes(joverlap,jproc), ierr)
              end do
          end do
          call mpi_win_free(comm%window, ierr)
      end if
    
      ! Flag indicating that the communication is complete
      comm%communication_complete=.true.
    
    end subroutine synchronize_onesided_communication


    !> Locreg communication
    subroutine communicate_locreg_descriptors_basics(iproc, nlr, rootarr, orbs, llr)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nlr
      integer,dimension(nlr),intent(in) :: rootarr
      type(orbitals_data),intent(in) :: orbs
      type(locreg_descriptors),dimension(nlr),intent(inout) :: llr
    
      ! Local variables
      integer:: ierr, istat, iall, ilr, iilr
      character(len=1),dimension(:),allocatable :: worksend_char, workrecv_char
      logical,dimension(:),allocatable :: worksend_log, workrecv_log
      integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
      real(8),dimension(:,:),allocatable :: worksend_dbl, workrecv_dbl
      character(len=*),parameter :: subname='communicate_locreg_descriptors_basics'
    
      allocate(worksend_char(orbs%norbp), stat=istat)
      call memocc(istat, worksend_char, 'worksend_char', subname)
      allocate(worksend_log(orbs%norbp), stat=istat)
      call memocc(istat, worksend_log, 'worksend_log', subname)
      allocate(worksend_int(11,orbs%norbp), stat=istat)
      call memocc(istat, worksend_int, 'worksend_int', subname)
      allocate(worksend_dbl(6,orbs%norbp), stat=istat)
      call memocc(istat, worksend_dbl, 'worksend_dbl', subname)
    
      allocate(workrecv_char(orbs%norb), stat=istat)
      call memocc(istat, workrecv_char, 'workrecv_char', subname)
      allocate(workrecv_log(orbs%norb), stat=istat)
      call memocc(istat, workrecv_log, 'workrecv_log', subname)
      allocate(workrecv_int(11,orbs%norb), stat=istat)
      call memocc(istat, workrecv_int, 'workrecv_int', subname)
      allocate(workrecv_dbl(6,orbs%norb), stat=istat)
      call memocc(istat, workrecv_dbl, 'workrecv_dbl', subname)
    
    
      iilr=0
      do ilr=1,nlr
          if (iproc==rootarr(ilr)) then
              iilr=iilr+1
              worksend_char(iilr)=llr(ilr)%geocode
              worksend_log(iilr)=llr(ilr)%hybrid_on
              worksend_int(1,iilr)=llr(ilr)%ns1
              worksend_int(2,iilr)=llr(ilr)%ns2
              worksend_int(3,iilr)=llr(ilr)%ns3
              worksend_int(4,iilr)=llr(ilr)%nsi1
              worksend_int(5,iilr)=llr(ilr)%nsi2
              worksend_int(6,iilr)=llr(ilr)%nsi3
              worksend_int(7,iilr)=llr(ilr)%localnorb
              worksend_int(8:10,iilr)=llr(ilr)%outofzone(1:3)
              worksend_int(11,iilr)=ilr
              worksend_dbl(1:3,iilr)=llr(ilr)%locregCenter(1:3)
              worksend_dbl(4,iilr)=llr(ilr)%locrad
              worksend_dbl(5,iilr)=llr(ilr)%locrad_kernel
              worksend_dbl(6,iilr)=llr(ilr)%locrad_mult
          end if
      end do
    
      call mpi_allgatherv(worksend_char, orbs%norbp, mpi_character, workrecv_char, orbs%norb_par(:,0), &
           orbs%isorb_par, mpi_character, bigdft_mpi%mpi_comm, ierr)
      call mpi_allgatherv(worksend_log, orbs%norbp, mpi_logical, workrecv_log, orbs%norb_par(:,0), &
           orbs%isorb_par, mpi_logical, bigdft_mpi%mpi_comm, ierr)
      call mpi_allgatherv(worksend_int, 11*orbs%norbp, mpi_integer, workrecv_int, 11*orbs%norb_par(:,0), &
           11*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
      call mpi_allgatherv(worksend_dbl, 6*orbs%norbp, mpi_double_precision, workrecv_dbl, 6*orbs%norb_par(:,0), &
           6*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      do ilr=1,nlr
          iilr=workrecv_int(11,ilr)
          llr(iilr)%geocode=workrecv_char(ilr)
          llr(iilr)%hybrid_on= workrecv_log(ilr)
          llr(iilr)%ns1=workrecv_int(1,ilr)
          llr(iilr)%ns2=workrecv_int(2,ilr)
          llr(iilr)%ns3=workrecv_int(3,ilr)
          llr(iilr)%nsi1=workrecv_int(4,ilr)
          llr(iilr)%nsi2=workrecv_int(5,ilr)
          llr(iilr)%nsi3=workrecv_int(6,ilr)
          llr(iilr)%localnorb=workrecv_int(7,ilr)
          llr(iilr)%outofzone(1:3)=workrecv_int(8:10,ilr)
          llr(iilr)%locregCenter(1:3)=workrecv_dbl(1:3,ilr)
          llr(iilr)%locrad=workrecv_dbl(4,ilr)
          llr(iilr)%locrad_kernel=workrecv_dbl(5,ilr)
          llr(iilr)%locrad_mult=workrecv_dbl(6,ilr)
      end do
    
    
      iall=-product(shape(worksend_int))*kind(worksend_int)
      deallocate(worksend_int,stat=istat)
      call memocc(istat, iall, 'worksend_int', subname)
      iall=-product(shape(workrecv_int))*kind(workrecv_int)
      deallocate(workrecv_int,stat=istat)
      call memocc(istat, iall, 'workrecv_int', subname)
      allocate(worksend_int(13,orbs%norbp), stat=istat)
      call memocc(istat, worksend_int, 'worksend_int', subname)
      allocate(workrecv_int(13,orbs%norb), stat=istat)
      call memocc(istat, workrecv_int, 'workrecv_int', subname)
    
    
      iilr=0
      do ilr=1,nlr
          if (iproc==rootarr(ilr)) then
              iilr=iilr+1
              worksend_int(1,iilr)=llr(ilr)%d%n1
              worksend_int(2,iilr)=llr(ilr)%d%n2
              worksend_int(3,iilr)=llr(ilr)%d%n3
              worksend_int(4,iilr)=llr(ilr)%d%nfl1
              worksend_int(5,iilr)=llr(ilr)%d%nfu1
              worksend_int(6,iilr)=llr(ilr)%d%nfl2
              worksend_int(7,iilr)=llr(ilr)%d%nfu2
              worksend_int(8,iilr)=llr(ilr)%d%nfl3
              worksend_int(9,iilr)=llr(ilr)%d%nfu3
              worksend_int(10,iilr)=llr(ilr)%d%n1i
              worksend_int(11,iilr)=llr(ilr)%d%n2i
              worksend_int(12,iilr)=llr(ilr)%d%n3i
              worksend_int(13,iilr)=ilr
          end if
      end do
    
      call mpi_allgatherv(worksend_int, 13*orbs%norbp, mpi_integer, workrecv_int, 13*orbs%norb_par(:,0), &
           13*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    
      do ilr=1,nlr
          iilr=workrecv_int(13,ilr)
          llr(iilr)%d%n1=workrecv_int(1,ilr)
          llr(iilr)%d%n2=workrecv_int(2,ilr)
          llr(iilr)%d%n3=workrecv_int(3,ilr)
          llr(iilr)%d%nfl1=workrecv_int(4,ilr)
          llr(iilr)%d%nfu1=workrecv_int(5,ilr)
          llr(iilr)%d%nfl2=workrecv_int(6,ilr)
          llr(iilr)%d%nfu2=workrecv_int(7,ilr)
          llr(iilr)%d%nfl3=workrecv_int(8,ilr)
          llr(iilr)%d%nfu3=workrecv_int(9,ilr)
          llr(iilr)%d%n1i=workrecv_int(10,ilr)
          llr(iilr)%d%n2i=workrecv_int(11,ilr)
          llr(iilr)%d%n3i=workrecv_int(12,ilr)
      end do
    
    
      iall=-product(shape(worksend_char))*kind(worksend_char)
      deallocate(worksend_char,stat=istat)
      call memocc(istat, iall, 'worksend_char', subname)
      iall=-product(shape(worksend_log))*kind(worksend_log)
      deallocate(worksend_log,stat=istat)
      call memocc(istat, iall, 'worksend_log', subname)
      iall=-product(shape(worksend_int))*kind(worksend_int)
      deallocate(worksend_int,stat=istat)
      call memocc(istat, iall, 'worksend_int', subname)
      iall=-product(shape(worksend_dbl))*kind(worksend_dbl)
      deallocate(worksend_dbl,stat=istat)
      call memocc(istat, iall, 'worksend_dbl', subname)
    
      iall=-product(shape(workrecv_char))*kind(workrecv_char)
      deallocate(workrecv_char,stat=istat)
      call memocc(istat, iall, 'workrecv_char', subname)
      iall=-product(shape(workrecv_log))*kind(workrecv_log)
      deallocate(workrecv_log,stat=istat)
      call memocc(istat, iall, 'workrecv_log', subname)
      iall=-product(shape(workrecv_int))*kind(workrecv_int)
      deallocate(workrecv_int,stat=istat)
      call memocc(istat, iall, 'workrecv_int', subname)
      iall=-product(shape(workrecv_dbl))*kind(workrecv_dbl)
      deallocate(workrecv_dbl,stat=istat)
      call memocc(istat, iall, 'workrecv_dbl', subname)
    
    end subroutine communicate_locreg_descriptors_basics
    
    
    subroutine communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, rootarr, onwhichmpi)
       use module_base
       use module_types
       use yaml_output
       implicit none
    
       ! Calling arguments
       integer,intent(in):: iproc, nproc, nlr
       type(locreg_descriptors),intent(in) :: glr
       type(locreg_descriptors),dimension(nlr),intent(inout) :: llr
       type(orbitals_data),intent(in) :: orbs
       integer,dimension(nlr),intent(in) :: rootarr
       integer,dimension(orbs%norb),intent(in) :: onwhichmpi
    
       ! Local variables
       integer:: ierr, istat, iall, jorb, ilr, jlr, jtask, root, icomm, nrecv, nalloc, max_sim_comms
       integer :: maxrecvdim, maxsenddim
       logical :: isoverlap
       character(len=*),parameter:: subname='communicate_wavefunctions_descriptors2'
       integer,dimension(:),allocatable :: requests
       integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
       logical,dimension(:,:),allocatable :: covered
       !integer :: total_sent, total_recv
    
       ! This maxval is put out of the allocate to avoid compiler crash with PathScale.
       jorb = maxval(orbs%norb_par(:,0))
       allocate(requests(8*nproc*jorb), stat=istat)
       call memocc(istat, requests, 'requests', subname)
    
       allocate(covered(nlr,0:nproc-1), stat=istat)
       call memocc(istat, covered, 'covered', subname)
    
       allocate(worksend_int(4,nlr), stat=istat)
       call memocc(istat, worksend_int, 'worksend_int', subname)
    
       allocate(workrecv_int(4,nlr), stat=istat)
       call memocc(istat, workrecv_int, 'workrecv_int', subname)
    
       nrecv=0
       !nsend=0
       icomm=0
       maxsenddim=0
       do ilr=1,nlr
           root=rootarr(ilr)
           covered(ilr,:)=.false.
           do jorb=1,orbs%norb
               jlr=orbs%inwhichlocreg(jorb)
               jtask=onwhichmpi(jorb)
               ! check we're on a sending or receiving proc
               if (iproc /= root .and. iproc /= jtask) cycle
               ! don't communicate to ourselves, or if we've already sent this locreg
               if (jtask == root .or. covered(ilr,jtask)) cycle
               call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
               if (isoverlap) then         
                   covered(ilr,jtask)=.true.
                   if (iproc == root) then
                      !write(*,'(5(a,i0))') 'process ',iproc,' sends locreg ',ilr,' to process ',&
                      !    jtask,' with tags ',4*ilr+0,'-',4*ilr+3
                      worksend_int(1,ilr)=llr(ilr)%wfd%nvctr_c
                      worksend_int(2,ilr)=llr(ilr)%wfd%nvctr_f
                      worksend_int(3,ilr)=llr(ilr)%wfd%nseg_c
                      worksend_int(4,ilr)=llr(ilr)%wfd%nseg_f
                      icomm=icomm+1
                      call mpi_isend(worksend_int(1,ilr), 4, mpi_integer, jtask,&
                           itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
                      maxsenddim=max(maxsenddim,llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
                      !nsend=nsend+1
                   else if (iproc == jtask) then
                      !write(*,'(5(a,i0))') 'process ',iproc,' receives locreg ',ilr,' from process ',&
                      !    root,' with tags ',4*ilr+0,'-',4*ilr+3
                      icomm=icomm+1
                      call mpi_irecv(workrecv_int(1,ilr), 4, mpi_integer, root,&
                           itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
                      nrecv=nrecv+1
                   end if
               end if
           end do
       end do
      
       call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
       call mpi_barrier(mpi_comm_world,ierr)
    
       iall=-product(shape(worksend_int))*kind(worksend_int)
       deallocate(worksend_int,stat=istat)
       call memocc(istat, iall, 'worksend_int', subname)
    
       nalloc=0
       maxrecvdim=0
       do jlr=1,nlr 
          if (covered(jlr,iproc)) then
             llr(jlr)%wfd%nvctr_c=workrecv_int(1,jlr)
             llr(jlr)%wfd%nvctr_f=workrecv_int(2,jlr)
             llr(jlr)%wfd%nseg_c=workrecv_int(3,jlr)
             llr(jlr)%wfd%nseg_f=workrecv_int(4,jlr)
    !         call allocate_wfd(llr(jlr)%wfd,subname)
             nalloc=nalloc+1
             maxrecvdim=max(maxrecvdim,llr(jlr)%wfd%nseg_c+llr(jlr)%wfd%nseg_f)
          end if
       end do
       if (f_err_raise(nalloc /= nrecv,'problem in communicate locregs: mismatch in receives '//&
            trim(yaml_toa(nrecv))//' and allocates '//trim(yaml_toa(nalloc))//' for process '//trim(yaml_toa(iproc)),&
            err_name='BIGDFT_RUNTIME_ERROR')) return
    
       iall=-product(shape(workrecv_int))*kind(workrecv_int)
       deallocate(workrecv_int,stat=istat)
       call memocc(istat, iall, 'workrecv_int', subname)
    
       !should reduce memory by not allocating for all llr
       allocate(workrecv_int(6*maxrecvdim,nlr), stat=istat)
       call memocc(istat, workrecv_int, 'workrecv_int', subname)
       allocate(worksend_int(6*maxsenddim,nlr), stat=istat)
       call memocc(istat, worksend_int, 'worksend_int', subname)
    
       ! divide communications into chunks to avoid problems with memory (too many communications)
       ! set maximum number of simultaneous communications
       !total_sent=0
       !total_recv=0
       max_sim_comms=10000
       icomm=0
       do ilr=1,nlr
          root=rootarr(ilr)
          do jtask=0,nproc-1
             if (.not. covered(ilr,jtask)) cycle
             if (iproc == root) then
               !write(*,'(5(a,i0))') 'process ',iproc,' sends locreg ',ilr,' to process ',&
               !     jtask,' with tags ',4*ilr+0,'-',4*ilr+3
               call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyglob(1,1),1,worksend_int(1,ilr),1)
               call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keygloc(1,1),1,&
                    worksend_int(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
               call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyvloc(1),1,&
                    worksend_int(4*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
               call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyvglob(1),1,&
                    worksend_int(5*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
               icomm=icomm+1
               call mpi_isend(worksend_int(1,ilr),6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                    jtask, itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
             else if (iproc == jtask) then
                !write(*,'(5(a,i0))') 'process ',iproc,' receives locreg ',ilr,' from process ',&
                !    root,' with tags ',4*ilr+0,'-',4*ilr+3
                icomm=icomm+1
                call mpi_irecv(workrecv_int(1,ilr),6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                     root, itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
             end if
          end do
          if (mod(ilr,max_sim_comms)==0 .or. ilr==nlr) then
             do jlr=max(ilr-max_sim_comms+1,1),ilr
                if (covered(jlr,iproc))  call allocate_wfd(llr(jlr)%wfd)
             end do
             call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
             if (f_err_raise(ierr /= 0,'problem in communicate locregs: error in mpi_waitall '//&
                  trim(yaml_toa(ierr))//' for process '//trim(yaml_toa(iproc)),&
                  err_name='BIGDFT_RUNTIME_ERROR')) return
             call mpi_barrier(mpi_comm_world,ierr)
             icomm=0
          end if
       end do
    
       iall=-product(shape(worksend_int))*kind(worksend_int)
       deallocate(worksend_int,stat=istat)
       call memocc(istat, iall, 'worksend_int', subname)
    
       do ilr=1,nlr 
          if (covered(ilr,iproc)) then
             call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),workrecv_int(1,ilr),1,llr(ilr)%wfd%keyglob(1,1),1)
             call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
                  workrecv_int(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keygloc(1,1),1)
             call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
                  workrecv_int(4*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keyvloc(1),1)
             call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
                  workrecv_int(5*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keyvglob(1),1)
          end if
       end do
    
       iall=-product(shape(workrecv_int))*kind(workrecv_int)
       deallocate(workrecv_int,stat=istat)
       call memocc(istat, iall, 'workrecv_int', subname)
    
       !print*,'iproc,sent,received,num sent,num received',iproc,total_sent,total_recv,nsend,nrecv
       iall=-product(shape(requests))*kind(requests)
       deallocate(requests,stat=istat)
       call memocc(istat, iall, 'requests', subname)
    
       iall=-product(shape(covered))*kind(covered)
       deallocate(covered,stat=istat)
       call memocc(istat, iall, 'covered', subname)
    
    contains
    
     pure function itag(ilr,recv)
     implicit none
     integer, intent(in) :: ilr,recv
     integer :: itag
    
     itag=ilr+recv*nlr
    
     end function itag
    
    END SUBROUTINE communicate_locreg_descriptors_keys
    
    
    !> Transposition of the arrays, variable version (non homogeneous)
    subroutine transpose_v(iproc,nproc,orbs,wfd,comms,psi_add,work_add,&
               out_add) !optional
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc
      type(orbitals_data), intent(in) :: orbs
      type(wavefunctions_descriptors), intent(in) :: wfd
      !type(local_zone_descriptors),intent(in) :: lzd
      type(comms_cubic), intent(in) :: comms
      !>address of the wavefunction elements (choice)
      !if out_add is absent, it is used for transpose
      real(wp), intent(inout) :: psi_add
      real(wp), intent(inout) :: work_add
      !> size of the buffers, optional.
      real(wp), optional :: out_add
      !local variables
      character(len=*), parameter :: subname='transpose_v'
      integer :: ierr
      external :: switch_waves_v,psitransspi,MPI_ALLTOALLV
    
      call timing(iproc,'Un-TransSwitch','ON')
    
      if (nproc > 1) then
         call switch_waves_v(nproc,orbs,&
              wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par(0,1),psi_add,work_add)
    
         call timing(iproc,'Un-TransSwitch','OF')
         call timing(iproc,'Un-TransComm  ','ON')
         if (present(out_add)) then
            call MPI_ALLTOALLV(work_add,comms%ncntd,comms%ndspld,mpidtypw, &
                 out_add,comms%ncntt,comms%ndsplt,mpidtypw,bigdft_mpi%mpi_comm,ierr)
         else
            call MPI_ALLTOALLV(work_add,comms%ncntd,comms%ndspld,mpidtypw, &
                 psi_add,comms%ncntt,comms%ndsplt,mpidtypw,bigdft_mpi%mpi_comm,ierr)
         end if
         call timing(iproc,'Un-TransComm  ','OF')
         call timing(iproc,'Un-TransSwitch','ON')
      else
         if(orbs%nspinor /= 1) then
            !for only one processor there is no need to transform this
            call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi_add,.true.)
         end if
      end if
    
      call timing(iproc,'Un-TransSwitch','OF')
    
    END SUBROUTINE transpose_v


    
    
    subroutine untranspose_v(iproc,nproc,orbs,wfd,comms,psi_add,&
         work_add,out_add) !optional
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc
      type(orbitals_data), intent(in) :: orbs
      type(wavefunctions_descriptors), intent(in) :: wfd
      type(comms_cubic), intent(in) :: comms
      !real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: psi
      real(wp),intent(inout) :: psi_add
      real(wp),intent(inout) :: work_add
      real(wp),intent(out),optional :: out_add !< Optional argument
      !local variables
      integer :: ierr
      external :: switch_waves_v,psitransspi,MPI_ALLTOALLV
    
    
      call timing(iproc,'Un-TransSwitch','ON')
    
      if (nproc > 1) then
         call timing(iproc,'Un-TransSwitch','OF')
         call timing(iproc,'Un-TransComm  ','ON')
         call MPI_ALLTOALLV(psi_add,comms%ncntt,comms%ndsplt,mpidtypw,  &
              work_add,comms%ncntd,comms%ndspld,mpidtypw,bigdft_mpi%mpi_comm,ierr)
         call timing(iproc,'Un-TransComm  ','OF')
         call timing(iproc,'Un-TransSwitch','ON')
         if (present(out_add)) then
            !!call unswitch_waves_v(nproc,orbs,&
            !!     wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par(0,1),work,outadd)
            call unswitch_waves_v(nproc,orbs,&
                 wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par(0,1),work_add,out_add)
         else
            !!call unswitch_waves_v(nproc,orbs,&
            !!     wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par(0,1),work,psi)
            call unswitch_waves_v(nproc,orbs,&
                 wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par(0,1),work_add,psi_add)
         end if
      else
         if(orbs%nspinor /= 1) then
            call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi_add,.false.)
         end if
      end if
    
      call timing(iproc,'Un-TransSwitch','OF')
    END SUBROUTINE untranspose_v
    


end module communications



subroutine switch_waves_v(nproc,orbs,nvctr,nvctr_par,psi,psiw)
  !n(c) use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,nvctr
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nproc,orbs%nkpts), intent(in) :: nvctr_par
  real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(out) :: psiw
  !local variables
  integer :: iorb,i,j,ij,ijproc,ind,it,it1,it2,it3,it4,ikptsp
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt


  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0
  do ikptsp =1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)
     !if (ikpt < orbs%iokpt(1) .or. ikpt > orbs%iokpt(orbs%norbp)) cycle

     !calculate the number of orbitals belonging to k-point ikptstp
     !calculate to which k-point it belongs
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikpt)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikpt)
                 it=ind+i
                 psiw(it)=psi(ij,1,iorb+isorbp)
                 ij=ij+1
              enddo
              ijproc=ijproc+nvctr_par(j,ikpt)
           enddo
        enddo
     else if (orbs%nspinor == 2) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikpt)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikpt)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 psiw(it1)=psi(ij,1,iorb+isorbp)
                 psiw(it2)=psi(ij,2,iorb+isorbp)
                 ij=ij+1
              enddo
              ijproc=ijproc+nvctr_par(j,ikpt)
           enddo
        enddo
     else if (orbs%nspinor == 4) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikpt)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikpt)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 it3=ind+2*i+2*nvctr_par(j,ikpt)-1
                 it4=ind+2*i+2*nvctr_par(j,ikpt)
                 psiw(it1)=psi(ij,1,iorb+isorbp)
                 psiw(it2)=psi(ij,2,iorb+isorbp)
                 psiw(it3)=psi(ij,3,iorb+isorbp)
                 psiw(it4)=psi(ij,4,iorb+isorbp)
                 ij=ij+1
              enddo
              ijproc=ijproc+nvctr_par(j,ikpt)
           enddo
        enddo
     end if
     !update starting orbitals
     isorb=isorb+norbp_kpt
     isorbp=isorbp+norbp_kpt
     !and starting point for psi
     ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
END SUBROUTINE switch_waves_v


subroutine unswitch_waves_v(nproc,orbs,nvctr,nvctr_par,psiw,psi)
  !n(c) use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,nvctr
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nproc,orbs%nkpts), intent(in) :: nvctr_par
  real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(in) :: psiw
  real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(out) :: psi
  !local variables
  integer :: iorb,i,j,ij,ijproc,ind,it,it1,it2,it3,it4,ikptsp
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt

  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)
     !if (ikpt < orbs%iokpt(1) .or. ikpt > orbs%iokpt(orbs%norbp)) cycle

     !calculate the number of orbitals belonging to k-point ikptstp
     !calculate to which k-point it belongs
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikpt)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikpt)
                 it=ind+i
                 psi(ij,orbs%nspinor,iorb+isorbp)=psiw(it)
                 ij=ij+1
              end do
              ijproc=ijproc+nvctr_par(j,ikpt)
           end do
        end do
     else if (orbs%nspinor == 2) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikpt)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikpt)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 psi(ij,1,iorb+isorbp)=psiw(it1)
                 psi(ij,2,iorb+isorbp)=psiw(it2)
                 ij=ij+1
              end do
              ijproc=ijproc+nvctr_par(j,ikpt)
           end do
        end do
     else if (orbs%nspinor == 4) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikpt)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikpt)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 it3=ind+2*i+2*nvctr_par(j,ikpt)-1
                 it4=ind+2*i+2*nvctr_par(j,ikpt)
                 psi(ij,1,iorb+isorbp)=psiw(it1)
                 psi(ij,2,iorb+isorbp)=psiw(it2)
                 psi(ij,3,iorb+isorbp)=psiw(it3)
                 psi(ij,4,iorb+isorbp)=psiw(it4)
                 ij=ij+1
              end do
              ijproc=ijproc+nvctr_par(j,ikpt)
           end do
        end do
     end if
     !update starting orbitals
     isorb=isorb+norbp_kpt
     isorbp=isorbp+norbp_kpt
     !and starting point for psi
     ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
  
END SUBROUTINE unswitch_waves_v



!> The cubic routines
subroutine psitransspi(nvctrp,orbs,psi,forward)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nvctrp
  logical, intent(in) :: forward
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(orbs%nspinor*nvctrp,orbs%norb,orbs%nkpts), intent(inout) :: psi
  !local variables
  character(len=*), parameter :: subname='psitransspi'
  integer :: i,iorb,isp,i_all,i_stat,ikpts
  real(wp), dimension(:,:,:,:), allocatable :: tpsit

  allocate(tpsit(nvctrp,orbs%nspinor,orbs%norb,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,tpsit,'tpsit',subname)
  if(forward) then
     !we can use vcopy here
     do ikpts=1,orbs%nkpts
        do iorb=1,orbs%norb
           do isp=1,orbs%nspinor
              do i=1,nvctrp
                 tpsit(i,isp,iorb,ikpts)=psi(i+(isp-1)*nvctrp,iorb,ikpts)
              enddo
           enddo
        enddo
     end do
     if (orbs%nspinor == 2) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 psi(2*i-1,iorb,ikpts)=tpsit(i,1,iorb,ikpts)
                 psi(2*i,iorb,ikpts)=tpsit(i,2,iorb,ikpts)
              enddo
           enddo
        end do
     else if (orbs%nspinor == 4) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 psi(2*i-1,iorb,ikpts)=tpsit(i,1,iorb,ikpts)
                 psi(2*i,iorb,ikpts)=tpsit(i,2,iorb,ikpts)
                 psi(2*i+2*nvctrp-1,iorb,ikpts)=tpsit(i,3,iorb,ikpts)
                 psi(2*i+2*nvctrp,iorb,ikpts)=tpsit(i,4,iorb,ikpts)
              enddo
           enddo
        end do
     end if
  else
     if (orbs%nspinor == 2) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 tpsit(i,1,iorb,ikpts)=psi(2*i-1,iorb,ikpts)
                 tpsit(i,2,iorb,ikpts)=psi(2*i,iorb,ikpts)
              enddo
           enddo
        end do
     else if (orbs%nspinor == 4) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 tpsit(i,1,iorb,ikpts)=psi(2*i-1,iorb,ikpts)
                 tpsit(i,2,iorb,ikpts)=psi(2*i,iorb,ikpts)
                 tpsit(i,3,iorb,ikpts)=psi(2*i-1+2*nvctrp,iorb,ikpts)
                 tpsit(i,4,iorb,ikpts)=psi(2*i+2*nvctrp,iorb,ikpts)
              enddo
           enddo
        end do
     end if

     !here we can use vcopy
     do ikpts=1,orbs%nkpts
        do iorb=1,orbs%norb
           do isp=1,orbs%nspinor
              do i=1,nvctrp
                 psi(i+(isp-1)*nvctrp,iorb,ikpts)=tpsit(i,isp,iorb,ikpts)
              enddo
           enddo
        enddo
     end do
  end if

  i_all=-product(shape(tpsit))*kind(tpsit)
  deallocate(tpsit,stat=i_stat)
  call memocc(i_stat,i_all,'tpsit',subname)
END SUBROUTINE psitransspi



!> Transposition of the arrays, variable version (non homogeneous)
subroutine toglobal_and_transpose(iproc,nproc,orbs,Lzd,comms,psi,&
     work,outadd) !optional
  use module_base
  use module_types
  use communications, only: transpose_v
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(comms_cubic), intent(in) :: comms
  real(wp), dimension(:), pointer :: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  character(len=*), parameter :: subname='toglobal_and_transpose'
  integer :: i_all,i_stat
  integer :: psishift1,totshift,iorb,ilr,ldim,Gdim
  real(wp) :: workdum
  real(wp), dimension(:), pointer :: workarr

  call timing(iproc,'Un-TransSwitch','ON')

  !for linear scaling must project the wavefunctions to whole simulation box
  if(Lzd%linear) then
!     if(.not. present(work) .or. .not. associated(work)) stop 'transpose_v needs optional argument work with Linear Scaling'
     psishift1 = 1
     totshift = 1
     Gdim = max((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor,&
           sum(comms%ncntt(0:nproc-1)))
     allocate(workarr(Gdim+ndebug),stat=i_stat)
     call memocc(i_stat,workarr,'workarr',subname)
     call to_zero(Gdim,workarr)
     do iorb=1,orbs%norbp
        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
        ldim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

        !!call Lpsi_to_global(Lzd%Glr,Gdim,Lzd%Llr(ilr),psi(psishift1),&
        !!     ldim,orbs%norbp,orbs%nspinor,orbs%nspin,totshift,workarr)
        call Lpsi_to_global2(iproc, ldim, gdim, orbs%norbp, orbs%nspinor, &
             orbs%nspin, lzd%glr, lzd%llr(ilr), psi(psishift1), workarr(totshift))
        psishift1 = psishift1 + ldim
        totshift = totshift + (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor
     end do

     !reallocate psi to the global dimensions
     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)
     allocate(psi(Gdim+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)
     call vcopy(Gdim,workarr(1),1,psi(1),1) !psi=work
     i_all=-product(shape(workarr))*kind(workarr)
     deallocate(workarr,stat=i_stat)
     call memocc(i_stat,i_all,'workarr',subname)
  end if

  if (nproc > 1 .and. .not. associated(work)) then
     call f_err_throw('The working pointer must be associated',&
          err_name='BIGDFT_RUNTIME_ERROR')
  end if

  if (present(outadd)) then
      call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),work(1),outadd(1))
  else
     if (.not. associated(work)) then
        call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),workdum)
     else
        call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),work(1))
     end if
  end if

  !!if (nproc > 1) then
  !!   !control check
  !!   if (.not. present(work) .or. .not. associated(work)) then
  !!      if(iproc == 0) write(*,'(1x,a)')&
  !!           "ERROR: Unproper work array for transposing in parallel"
  !!      stop
  !!   end if


  !!   !!call switch_waves_v(nproc,orbs,&
  !!   !!     Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,comms%nvctr_par(0,1),psi,work)
  !!   call switch_waves_v(nproc,orbs,&
  !!        Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,comms%nvctr_par,psi,work)

  !!   call timing(iproc,'Un-TransSwitch','OF')
  !!   call timing(iproc,'Un-TransComm  ','ON')
  !!   if (present(outadd)) then
  !!      call MPI_ALLTOALLV(work,comms%ncntd,comms%ndspld,mpidtypw, &
  !!           outadd,comms%ncntt,comms%ndsplt,mpidtypw,bigdft_mpi%mpi_comm,ierr)
  !!   else
  !!      call MPI_ALLTOALLV(work,comms%ncntd,comms%ndspld,mpidtypw, &
  !!           psi,comms%ncntt,comms%ndsplt,mpidtypw,bigdft_mpi%mpi_comm,ierr)
  !!   end if

  !!   call timing(iproc,'Un-TransComm  ','OF')
  !!   call timing(iproc,'Un-TransSwitch','ON')
  !!else
  !!   if(orbs%nspinor /= 1) then
  !!      !for only one processor there is no need to transform this
  !!      call psitransspi(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs,psi,.true.)
  !!   end if
  !!end if

  !!call timing(iproc,'Un-TransSwitch','OF')

END SUBROUTINE toglobal_and_transpose
