!> @file
!!  File defining the structures and the routines for the communication between processes
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining routines related to communications (mainly transpositions)
module communications

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
  public :: toglobal_and_transpose


  contains


    subroutine transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
      use module_types, only: orbitals_data, local_zone_descriptors
      use wrapper_linalg, only: vcopy
      use dynamic_memory
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psi
      real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
      type(local_zone_descriptors),intent(in),optional :: lzd
      
      ! Local variables
      integer :: i_tot, i_c, i_f, iorb, iiorb, ilr, i, ind, m, ind7, i7
      real(kind=8),dimension(:),allocatable :: psi_c, psi_f
      character(len=*),parameter :: subname='transpose_switch_psi'

      call f_routine(id='transpose_switch_psi')
    
      psi_c = f_malloc(collcom%ndimpsi_c,id='psi_c')
      psi_f = f_malloc(7*collcom%ndimpsi_f,id='psi_f')
    
    
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
    
    
      call f_free(psi_c)
      call f_free(psi_f)

      call f_release_routine()
      
    end subroutine transpose_switch_psi



    subroutine transpose_communicate_psi(iproc, nproc, collcom, transpose_action, &
               psiwork_c, psiwork_f, wt, psitwork_c, psitwork_f)
      use module_base
      use communications_base, only: work_transpose, TRANSPOSE_FULL, TRANSPOSE_POST, &
                                     TRANSPOSE_GATHER
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, transpose_action
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimpsi_c),intent(in) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(in) :: psiwork_f
      type(work_transpose),intent(inout) :: wt
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psitwork_f
      
      ! Local variables
      integer :: ist, ist_c, ist_f, iisend, iirecv, jproc
      !integer :: ierr
      !!real(kind=8),dimension(:),allocatable :: psiwork, psitwork
      !!integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      !!character(len=*),parameter :: subname='transpose_communicate_psi'
    
      !call mpi_comm_size(bigdft_mpi%mpi_comm, nproc, ierr)
      !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)

      call f_routine(id='transpose_communicate_psi')
    
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_POST) then
          !!wt%psiwork = f_malloc_ptr(max(collcom%ndimpsi_c+7*collcom%ndimpsi_f,1),id='wt%psiwork')
          !!wt%psitwork = f_malloc_ptr(max(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f),1),id='wt%psitwork')
          !!wt%nsendcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nsendcounts')
          !!wt%nsenddspls = f_malloc_ptr(0.to.nproc-1,id='wt%nsenddspls')
          !!wt%nrecvcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvcounts')
          !!wt%nrecvdspls = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvdspls')
    
          ist=1
          ist_c=1
          ist_f=1
          iisend=0
          iirecv=0
          do jproc=0,nproc-1
              if(collcom%nsendcounts_c(jproc)>0) then
                  call vcopy(collcom%nsendcounts_c(jproc), psiwork_c(ist_c), 1, wt%psiwork(ist), 1)
              end if
              ist_c=ist_c+collcom%nsendcounts_c(jproc)
              ist=ist+collcom%nsendcounts_c(jproc)
              if(collcom%nsendcounts_f(jproc)>0) then
                  call vcopy(7*collcom%nsendcounts_f(jproc), psiwork_f(ist_f), 1, wt%psiwork(ist), 1)
              end if
              ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
              ist=ist+7*collcom%nsendcounts_f(jproc)
              wt%nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
              wt%nsenddspls(jproc)=iisend
              wt%nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
              wt%nrecvdspls(jproc)=iirecv
              iisend=iisend+wt%nsendcounts(jproc)
              iirecv=iirecv+wt%nrecvcounts(jproc)
          end do
    
          !write(*,'(a,i4,4x,100i8)') 'iproc, nsendcounts', iproc, nsendcounts
          !write(*,'(a,i4,4x,100i8)') 'iproc, nsenddspls', iproc, nsenddspls
          !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvcounts', iproc, nrecvcounts
          !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvdspls', iproc, nrecvdspls
          
          !!! coarse part
          !!call mpi_alltoallv(psiwork_c, collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, psitwork_c, &
          !!     collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          !!
          !!! fine part
          !!call mpi_alltoallv(psiwork_f, 7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, psitwork_f, &
          !!     7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          !!call mpi_alltoallv(wt%psiwork, wt%nsendcounts, wt%nsenddspls, mpi_double_precision, wt%psitwork, &
          !!     wt%nrecvcounts, wt%nrecvdspls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          if (nproc>1) then
              call mpiialltoallv(wt%psiwork(1), wt%nsendcounts(0), wt%nsenddspls(0), mpi_double_precision, wt%psitwork(1), &
                   wt%nrecvcounts(0), wt%nrecvdspls(0), mpi_double_precision, bigdft_mpi%mpi_comm, wt%request)
          else
              call vcopy(wt%nsendcounts(0), wt%psiwork(1), 1, wt%psitwork(1), 1)
              wt%request = MPI_REQUEST_NULL
          end if
      end if
    
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_GATHER) then

          if (nproc>1) then
              call mpiwait(wt%request)
          end if

          ist=1
          ist_c=1
          ist_f=1
          do jproc=0,nproc-1
              if(collcom%nrecvcounts_c(jproc)>0) then
                  call vcopy(collcom%nrecvcounts_c(jproc), wt%psitwork(ist), 1, psitwork_c(ist_c), 1)
              end if
              ist_c=ist_c+collcom%nrecvcounts_c(jproc)
              ist=ist+collcom%nrecvcounts_c(jproc)
              if(collcom%nrecvcounts_f(jproc)>0) then
                  call vcopy(7*collcom%nrecvcounts_f(jproc), wt%psitwork(ist), 1, psitwork_f(ist_f), 1)
              end if
              ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
              ist=ist+7*collcom%nrecvcounts_f(jproc)
          end do
    
          !!call f_free_ptr(wt%psiwork)
          !!call f_free_ptr(wt%psitwork)
          !!call f_free_ptr(wt%nsendcounts)
          !!call f_free_ptr(wt%nsenddspls)
          !!call f_free_ptr(wt%nrecvcounts)
          !!call f_free_ptr(wt%nrecvdspls)
      end if
    
      call f_release_routine()
    
    end subroutine transpose_communicate_psi



    subroutine transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
      use module_base
      implicit none
      
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
      
      ! Local variables
      integer :: i,ind,sum_c,sum_f,m,i7,ind7

      call f_routine(id='transpose_unswitch_psit')
    
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

      call f_release_routine()
    
    end subroutine transpose_unswitch_psit



    subroutine transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
      use module_base
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psitwork_f
      
      ! Local variables
      integer :: i, ind, sum_c,sum_f,m, i7, ind7

      call f_routine(id='transpose_switch_psit')
    
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
          ind=collcom%iexpand_f(i)
          ind7=7*ind
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

      call f_release_routine()
    
    end subroutine transpose_switch_psit



    subroutine transpose_communicate_psit(iproc, nproc, collcom, transpose_action, &
               psitwork_c, psitwork_f, wt, psiwork_c, psiwork_f)
      use module_base
      use communications_base, only: work_transpose, TRANSPOSE_FULL, TRANSPOSE_POST, &
                                     TRANSPOSE_GATHER
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, transpose_action
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      type(work_transpose),intent(inout) :: wt
      real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
      
      ! Local variables
      integer :: ist, ist_c, ist_f, jproc, iisend, iirecv
      !integer :: ierr
      !!real(kind=8),dimension(:),allocatable :: psiwork, psitwork
      !!integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      !!character(len=*),parameter :: subname='transpose_communicate_psit'
    
      !call mpi_comm_size(bigdft_mpi%mpi_comm, nproc, ierr)
      !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)

      call f_routine(id='transpose_communicate_psit')
    
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_POST) then
          !!wt%psiwork = f_malloc_ptr(max(collcom%ndimpsi_c+7*collcom%ndimpsi_f,1),id='wt%psiwork')
          !!wt%psitwork = f_malloc_ptr(max(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f),1),id='wt%psitwork')
          !!wt%nsendcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nsendcounts')
          !!wt%nsenddspls = f_malloc_ptr(0.to.nproc-1,id='wt%nsenddspls')
          !!wt%nrecvcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvcounts')
          !!wt%nrecvdspls = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvdspls')
    
          ist=1
          ist_c=1
          ist_f=1
          iisend=0
          iirecv=0
          do jproc=0,nproc-1
              if(collcom%nrecvcounts_c(jproc)>0) then
                  call vcopy(collcom%nrecvcounts_c(jproc), psitwork_c(ist_c), 1, wt%psitwork(ist), 1)
              end if
              ist_c=ist_c+collcom%nrecvcounts_c(jproc)
              ist=ist+collcom%nrecvcounts_c(jproc)
              if(collcom%nrecvcounts_f(jproc)>0) then
                  call vcopy(7*collcom%nrecvcounts_f(jproc), psitwork_f(ist_f), 1, wt%psitwork(ist), 1)
              end if
              ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
              ist=ist+7*collcom%nrecvcounts_f(jproc)
              wt%nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
              wt%nsenddspls(jproc)=iisend
              wt%nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
              wt%nrecvdspls(jproc)=iirecv
              iisend=iisend+wt%nsendcounts(jproc)
              iirecv=iirecv+wt%nrecvcounts(jproc)
          end do
    
          !!! coarse part
          !! call mpi_alltoallv(psitwork_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, psiwork_c, &
          !!      collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
          !!! fine part
          !! call mpi_alltoallv(psitwork_f, 7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, psiwork_f, &
          !!      7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          if (nproc>1) then
              call mpiialltoallv(wt%psitwork(1), wt%nrecvcounts(0), wt%nrecvdspls(0), mpi_double_precision, wt%psiwork(1), &
                   wt%nsendcounts(0), wt%nsenddspls(0), mpi_double_precision, bigdft_mpi%mpi_comm, wt%request)
          else
              call vcopy(wt%nrecvcounts(0), wt%psitwork(1), 1, wt%psiwork(1), 1)
              wt%request = MPI_REQUEST_NULL
          end if
      end if
    
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_GATHER) then

          if (nproc>1) then
              call mpiwait(wt%request)
          end if

          ist=1
          ist_c=1
          ist_f=1
          do jproc=0,nproc-1
              if(collcom%nsendcounts_c(jproc)>0) then
                  call vcopy(collcom%nsendcounts_c(jproc), wt%psiwork(ist), 1, psiwork_c(ist_c), 1)
              end if
              ist_c=ist_c+collcom%nsendcounts_c(jproc)
              ist=ist+collcom%nsendcounts_c(jproc)
              if(collcom%nsendcounts_f(jproc)>0) then
                  call vcopy(7*collcom%nsendcounts_f(jproc), wt%psiwork(ist), 1, psiwork_f(ist_f), 1)
              end if
              ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
              ist=ist+7*collcom%nsendcounts_f(jproc)
          end do
    
          !!call f_free_ptr(wt%psiwork)
          !!call f_free_ptr(wt%psitwork)
          !!call f_free_ptr(wt%nsendcounts)
          !!call f_free_ptr(wt%nsenddspls)
          !!call f_free_ptr(wt%nrecvcounts)
          !!call f_free_ptr(wt%nrecvdspls)
      end if

      call f_release_routine()
    
    end subroutine transpose_communicate_psit



    subroutine transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
      use module_types, only: orbitals_data, local_zone_descriptors
      use wrapper_linalg, only: vcopy
      use dynamic_memory
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
      integer :: i, ind, iorb, iiorb, ilr, i_tot, i_c, i_f, m, i7, ind7
      real(kind=8),dimension(:),allocatable :: psi_c, psi_f
      character(len=*),parameter :: subname='transpose_unswitch_psi'
      
      call f_routine(id='transpose_unswitch_psi')
      
      psi_c = f_malloc(collcom%ndimpsi_c,id='psi_c')
      psi_f = f_malloc(7*collcom%ndimpsi_f,id='psi_f')
      
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
      
      call f_free(psi_c)
      call f_free(psi_f)

      call f_release_routine()
    
    end subroutine transpose_unswitch_psi



    subroutine transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
               transpose_action, psi, psit_c, psit_f, lzd, wt_)
      use module_types, only: orbitals_data, local_zone_descriptors
      use dynamic_memory
      use dictionaries, only: f_err_throw,f_err_define
      use communications_base, only: work_transpose, work_transpose_null, &
                                     TRANSPOSE_FULL, TRANSPOSE_POST, &
                                     TRANSPOSE_GATHER, ERR_LINEAR_TRANSPOSITION
      !use module_interfaces, except_this_one => transpose_localized
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs, transpose_action
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psi
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
      type(local_zone_descriptors),intent(in) :: lzd
      type(work_transpose),intent(inout),target,optional :: wt_
      
      ! Local variables
      !!real(kind=8),dimension(:),allocatable :: psiwork_c, psiwork_f, psitwork_c, psitwork_f
      character(len=*),parameter :: subname='transpose_localized'
      type(work_transpose),pointer :: wt

      call timing(iproc,'Un-TransSwitch','ON')
      call f_routine(id='transpose_localized')

      ! Check the arguments
      if (transpose_action /= TRANSPOSE_FULL .and. &
          transpose_action /= TRANSPOSE_POST .and. &
          transpose_action /= TRANSPOSE_GATHER)  then
          call f_err_throw('transpose_localized was called with errorneous arguments',&
               err_id=ERR_LINEAR_TRANSPOSITION)
      end if

      if (transpose_action == TRANSPOSE_POST .or. &
          transpose_action == TRANSPOSE_GATHER) then
          if (.not.present(wt_)) then
              call f_err_throw('transpose_localized was called with errorneous arguments',&
                   err_id=ERR_LINEAR_TRANSPOSITION)
          end if
      end if

      ! Point to the provided work arrays, otherwise allocate
      if (present(wt_)) then
          wt => wt_
          if (.not.associated(wt%psiwork_c) .or. &
              .not.associated(wt%psiwork_c) .or. &
              .not.associated(wt%psiwork_f) .or. &
              .not.associated(wt%psitwork_c) .or. &
              .not.associated(wt%psitwork_f) .or. &
              .not.associated(wt%psiwork) .or. &
              .not.associated(wt%psitwork) .or. &
              .not.associated(wt%nsendcounts) .or. &
              .not.associated(wt%nsenddspls) .or. &
              .not.associated(wt%nrecvcounts) .or. &
              .not.associated(wt%nrecvdspls)) then
              call f_err_throw('not all workarrays are associated', err_id=ERR_LINEAR_TRANSPOSITION)
          end if
      else
          if (transpose_action == TRANSPOSE_FULL .or. &
              transpose_action == TRANSPOSE_POST) then
              allocate(wt)
              wt = work_transpose_null()
              wt%psiwork_c = f_malloc_ptr(collcom%ndimpsi_c,id='psiwork_c')
              wt%psiwork_f = f_malloc_ptr(7*collcom%ndimpsi_f,id='psiwork_f')
              wt%psitwork_c = f_malloc_ptr(collcom%ndimind_c,id='psitwork_c')
              wt%psitwork_f = f_malloc_ptr(7*collcom%ndimind_f,id='psitwork_f')
              wt%psiwork = f_malloc_ptr(max(collcom%ndimpsi_c+7*collcom%ndimpsi_f,1),id='wt%psiwork')
              wt%psitwork = f_malloc_ptr(max(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f),1),id='wt%psitwork')
              wt%nsendcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nsendcounts')
              wt%nsenddspls = f_malloc_ptr(0.to.nproc-1,id='wt%nsenddspls')
              wt%nrecvcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvcounts')
              wt%nrecvdspls = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvdspls')
          end if
      end if


      
      !!psiwork_c = f_malloc(collcom%ndimpsi_c,id='psiwork_c')
      !!psiwork_f = f_malloc(7*collcom%ndimpsi_f,id='psiwork_f')
      !!psitwork_c = f_malloc(collcom%ndimind_c,id='psitwork_c')
      !!psitwork_f = f_malloc(7*collcom%ndimind_f,id='psitwork_f')
      !##wt%psiwork_c = f_malloc_ptr(collcom%ndimpsi_c,id='psiwork_c')
      !##wt%psiwork_f = f_malloc_ptr(7*collcom%ndimpsi_f,id='psiwork_f')
      !##wt%psitwork_c = f_malloc_ptr(collcom%ndimind_c,id='psitwork_c')
      !##wt%psitwork_f = f_malloc_ptr(7*collcom%ndimind_f,id='psitwork_f')

      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_POST) then
          call transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, wt%psiwork_c, wt%psiwork_f, lzd)
          !#wt%psiwork = f_malloc_ptr(max(collcom%ndimpsi_c+7*collcom%ndimpsi_f,1),id='wt%psiwork')
          !#wt%psitwork = f_malloc_ptr(max(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f),1),id='wt%psitwork')
          !#wt%nsendcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nsendcounts')
          !#wt%nsenddspls = f_malloc_ptr(0.to.nproc-1,id='wt%nsenddspls')
          !#wt%nrecvcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvcounts')
          !#wt%nrecvdspls = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvdspls')
      end if
      call timing(iproc,'Un-TransSwitch','OF')
    
      call timing(iproc,'Un-TransComm  ','ON')
      call transpose_communicate_psi(iproc, nproc, collcom, transpose_action, &
           wt%psiwork_c, wt%psiwork_f, wt, wt%psitwork_c, wt%psitwork_f)
      call timing(iproc,'Un-TransComm  ','OF')
    
      call timing(iproc,'Un-TransSwitch','ON')
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_GATHER) then
          !##call f_free_ptr(wt%psiwork)
          !##call f_free_ptr(wt%psitwork)
          !##call f_free_ptr(wt%nsendcounts)
          !##call f_free_ptr(wt%nsenddspls)
          !##call f_free_ptr(wt%nrecvcounts)
          !##call f_free_ptr(wt%nrecvdspls)
          call transpose_unswitch_psit(collcom, wt%psitwork_c, wt%psitwork_f, psit_c, psit_f)
      end if

      
      !!call f_free(psiwork_c)
      !!call f_free(psiwork_f)
      !!call f_free(psitwork_c)
      !!call f_free(psitwork_f)
      !##call f_free_ptr(wt%psiwork_c)
      !##call f_free_ptr(wt%psiwork_f)
      !##call f_free_ptr(wt%psitwork_c)
      !##call f_free_ptr(wt%psitwork_f)

      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_GATHER) then
          if (.not.present(wt_)) then
              call f_free_ptr(wt%psiwork)
              call f_free_ptr(wt%psitwork)
              call f_free_ptr(wt%nsendcounts)
              call f_free_ptr(wt%nsenddspls)
              call f_free_ptr(wt%nrecvcounts)
              call f_free_ptr(wt%nrecvdspls)
              call f_free_ptr(wt%psiwork_c)
              call f_free_ptr(wt%psiwork_f)
              call f_free_ptr(wt%psitwork_c)
              call f_free_ptr(wt%psitwork_f)
              deallocate(wt)
              nullify(wt)
          end if
      end if

      call f_release_routine()
      call timing(iproc,'Un-TransSwitch','OF')
      
    end subroutine transpose_localized



    subroutine untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
               transpose_action, psit_c, psit_f, psi, lzd, wt_)
      use module_types, only: orbitals_data, local_zone_descriptors
      use dynamic_memory
      use dictionaries, only: f_err_throw,f_err_define
      use communications_base, only: work_transpose, work_transpose_null, &
                                     TRANSPOSE_FULL, TRANSPOSE_POST, &
                                     TRANSPOSE_GATHER, ERR_LINEAR_TRANSPOSITION
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs, transpose_action
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
      real(kind=8),dimension(npsidim_orbs),intent(out) :: psi
      type(local_zone_descriptors),intent(in) :: lzd
      type(work_transpose),intent(inout),target,optional :: wt_
      
      ! Local variables
      type(work_transpose),pointer :: wt

      call f_routine(id='untranspose_localized')

      ! Check the arguments
      if (transpose_action /= TRANSPOSE_FULL .and. &
          transpose_action /= TRANSPOSE_POST .and. &
          transpose_action /= TRANSPOSE_GATHER)  then
          call f_err_throw('untranspose_localized was called with errorneous arguments',&
               err_id=ERR_LINEAR_TRANSPOSITION)
      end if

      if (transpose_action == TRANSPOSE_POST .or. &
          transpose_action == TRANSPOSE_GATHER) then
          if (.not.present(wt_)) then
              call f_err_throw('untranspose_localized was called with errorneous arguments',&
                   err_id=ERR_LINEAR_TRANSPOSITION)
          end if
      end if

      ! Point to the provided work arrays, otherwise allocate
      if (present(wt_)) then
          wt => wt_
          if (.not.associated(wt%psiwork_c) .or. &
              .not.associated(wt%psiwork_c) .or. &
              .not.associated(wt%psiwork_f) .or. &
              .not.associated(wt%psitwork_c) .or. &
              .not.associated(wt%psitwork_f) .or. &
              .not.associated(wt%psiwork) .or. &
              .not.associated(wt%psitwork) .or. &
              .not.associated(wt%nsendcounts) .or. &
              .not.associated(wt%nsenddspls) .or. &
              .not.associated(wt%nrecvcounts) .or. &
              .not.associated(wt%nrecvdspls)) then
              call f_err_throw('not all workarrays are associated', err_id=ERR_LINEAR_TRANSPOSITION)
          end if
      else
          if (transpose_action == TRANSPOSE_FULL .or. &
              transpose_action == TRANSPOSE_POST) then
              allocate(wt)
              wt = work_transpose_null()
              wt%psiwork_c = f_malloc_ptr(collcom%ndimpsi_c,id='psiwork_c')
              wt%psiwork_f = f_malloc_ptr(7*collcom%ndimpsi_f,id='psiwork_f')
              wt%psitwork_c = f_malloc_ptr(collcom%ndimind_c,id='psitwork_c')
              wt%psitwork_f = f_malloc_ptr(7*collcom%ndimind_f,id='psitwork_f')
              wt%psiwork = f_malloc_ptr(max(collcom%ndimpsi_c+7*collcom%ndimpsi_f,1),id='wt%psiwork')
              wt%psitwork = f_malloc_ptr(max(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f),1),id='wt%psitwork')
              wt%nsendcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nsendcounts')
              wt%nsenddspls = f_malloc_ptr(0.to.nproc-1,id='wt%nsenddspls')
              wt%nrecvcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvcounts')
              wt%nrecvdspls = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvdspls')
          end if
      end if

      !##if (transpose_action == TRANSPOSE_FULL .or. &
      !##    transpose_action == TRANSPOSE_POST) then
      !##    wt = work_transpose_null()
      !##end if
      
      !##wt%psiwork_c = f_malloc_ptr(collcom%ndimpsi_c,id='psiwork_c')
      !##wt%psiwork_f = f_malloc_ptr(7*collcom%ndimpsi_f,id='psiwork_f')
      !##wt%psitwork_c = f_malloc_ptr(collcom%ndimind_c,id='psitwork_c')
      !##wt%psitwork_f = f_malloc_ptr(7*collcom%ndimind_f,id='psitwork_f')
    
      call timing(iproc,'Un-TransSwitch','ON')
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_POST) then
          call transpose_switch_psit(collcom, psit_c, psit_f, wt%psitwork_c, wt%psitwork_f)
          !##wt%psiwork = f_malloc_ptr(max(collcom%ndimpsi_c+7*collcom%ndimpsi_f,1),id='wt%psiwork')
          !##wt%psitwork = f_malloc_ptr(max(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f),1),id='wt%psitwork')
          !##wt%nsendcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nsendcounts')
          !##wt%nsenddspls = f_malloc_ptr(0.to.nproc-1,id='wt%nsenddspls')
          !##wt%nrecvcounts = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvcounts')
          !##wt%nrecvdspls = f_malloc_ptr(0.to.nproc-1,id='wt%nrecvdspls')
      end if
      call timing(iproc,'Un-TransSwitch','OF')
    
      call timing(iproc,'Un-TransComm  ','ON')
      call transpose_communicate_psit(iproc, nproc, collcom, transpose_action, &
           wt%psitwork_c, wt%psitwork_f, wt, wt%psiwork_c, wt%psiwork_f)
      call timing(iproc,'Un-TransComm  ','OF')
    
      call timing(iproc,'Un-TransSwitch','ON')
      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_GATHER) then
          !##call f_free_ptr(wt%psiwork)
          !##call f_free_ptr(wt%psitwork)
          !##call f_free_ptr(wt%nsendcounts)
          !##call f_free_ptr(wt%nsenddspls)
          !##call f_free_ptr(wt%nrecvcounts)
          !##call f_free_ptr(wt%nrecvdspls)
          call transpose_unswitch_psi(npsidim_orbs, orbs, collcom, wt%psiwork_c, wt%psiwork_f, psi, lzd)
      end if
      call timing(iproc,'Un-TransSwitch','OF')
      
      !##call f_free_ptr(wt%psiwork_c)
      !##call f_free_ptr(wt%psiwork_f)
      !##call f_free_ptr(wt%psitwork_c)
      !##call f_free_ptr(wt%psitwork_f)

      if (transpose_action == TRANSPOSE_FULL .or. &
          transpose_action == TRANSPOSE_GATHER) then
          if (.not.present(wt_)) then
              call f_free_ptr(wt%psiwork)
              call f_free_ptr(wt%psitwork)
              call f_free_ptr(wt%nsendcounts)
              call f_free_ptr(wt%nsenddspls)
              call f_free_ptr(wt%nrecvcounts)
              call f_free_ptr(wt%nrecvdspls)
              call f_free_ptr(wt%psiwork_c)
              call f_free_ptr(wt%psiwork_f)
              call f_free_ptr(wt%psitwork_c)
              call f_free_ptr(wt%psitwork_f)
              deallocate(wt)
              nullify(wt)
          end if
      end if

      call f_release_routine()
      
    end subroutine untranspose_localized


    subroutine transpose_switch_psir(collcom_sr, psir, psirwork)
      use module_base
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psir
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork
    
      ! Local variables
      integer :: i, m, ind

      call f_routine(id='transpose_switch_psir')
    
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
    
      call f_release_routine()

    end subroutine transpose_switch_psir

    
    subroutine transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)
      use module_base
      use wrapper_linalg, only: vcopy
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork
    
      ! Local variables
      integer :: ierr
    
      call f_routine(id='transpose_communicate_psir')
    
      if (nproc>1) then
          call mpi_alltoallv(psirwork, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, psirtwork, &
               collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         !call vcopy(collcom_sr%ndimpsi_c, psirwork(1), 1, psirtwork(1), 1)
         call f_memcpy(src=psirwork,dest=psirtwork)
      end if
    
      call f_release_routine()
    
    end subroutine transpose_communicate_psir
    
    subroutine transpose_unswitch_psirt(collcom_sr, psirtwork, psirt)
      use module_base
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirt
    
      ! Local variables
      integer :: i, ind, sum_c, m

      call f_routine(id='transpose_unswitch_psirt')
    
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

      call f_release_routine()
    
    end subroutine transpose_unswitch_psirt
    
    
    subroutine start_onesided_communication(iproc, nproc, n1, n2, n3p, sendbuf, nrecvbuf, recvbuf, comm, lzd)
      use module_base
      use module_types, only: local_zone_descriptors
      use communications_base, only: p2pComms, bgq, RMA_SYNC_ACTIVE, RMA_SYNC_PASSIVE, rma_sync, &
                                     TYPES_NESTED, TYPES_SIMPLE, type_strategy
      implicit none
      
      ! Calling arguments
      integer, intent(in):: iproc, nproc, n1, n2, nrecvbuf
      integer,dimension(0:nproc-1),intent(in) :: n3p
      type(p2pComms), intent(inout):: comm
      real(kind=8), dimension(n1*n2*n3p(iproc)*comm%nspin), intent(in):: sendbuf
      real(kind=8), dimension(nrecvbuf), intent(out):: recvbuf
      type(local_zone_descriptors), intent(in) :: lzd
      
      ! Local variables
      !character(len=*), parameter :: subname='start_onesided_communication'
      integer :: joverlap, mpisource, istsource, mpidest, istdest, ierr, nit, ispin, ispin_shift
      integer :: ioffset_send, ist, i2, i3, ist2, ist3, info, nsize, size_of_double, isend_shift
      integer :: i,iel,ii,ind,it,ncount
      !integer :: islices, ilines, ist1, ish1, ish2
      integer,dimension(:),allocatable :: npotarr, blocklengths, types
      integer(kind=mpi_address_kind),dimension(:),allocatable :: displacements
      integer(kind=mpi_address_kind) :: lb, extent


!!integer :: itemsize
!!type(c_ptr) :: baseptr

      !!! if on BG/Q avoid mpi_get etc. as unstable
      !!logical, parameter :: bgq=.false.


      !!do ist=1,nsendbuf
      !!    write(5400,'(a,2i12,es18.7)') 'iproc, ist, sendbuf(ist)', iproc, ist, sendbuf(ist)
      !!end do
      !recvbuf=-123456789.d0

      !write(*,'(a,i12,es16.8)') 'in start_onesided_communication: nsendbuf, sum(sendbuf)', nsendbuf, sum(sendbuf)
    
    
      call f_routine(id='start_onesided_communication')
      call timing(iproc, 'Pot_comm start', 'ON')

      blocklengths = f_malloc(comm%onedtypeovrlp*maxval(comm%comarr(5,:)),id='blocklengths')
      displacements = f_malloc(comm%onedtypeovrlp*maxval(comm%comarr(5,:)),id='displacements')
      types = f_malloc(comm%onedtypeovrlp*maxval(comm%comarr(5,:)),id='types')
      types(:) = mpi_double_precision

!!call mpi_type_extent(mpi_double_precision, itemsize, ierr)
!!!call mpi_type_size(mpi_double_precision, itemsize, ierr)
!!call mpi_alloc_mem(int(n1*n2*n3p(iproc)*comm%nspin*itemsize,kind=mpi_address_kind), MPI_INFO_NULL, baseptr, ierr)
!!call c_f_pointer(baseptr, sendbuf, (/n1*n2*n3p(iproc)*comm%nspin/))

      ! the size of the potential without spin (maybe need to find a better way to determine this...)
      npotarr = f_malloc(0.to.nproc-1,id='npotarr')
      npotarr(0:nproc-1)=n1*n2*n3p(0:nproc-1)
      !!if (nproc>1) then
      !!    call mpiallred(npotarr(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
      !!end if
      !npot=nsendbuf/comm%nspin
    
      if(.not.comm%communication_complete) then
          call f_err_throw('ERROR: there is already a p2p communication going on...')
      end if

      !nproc_if: if (nproc>1) then

          spin_loop: do ispin=1,comm%nspin

              ispin_shift = (ispin-1)*comm%nrecvbuf
    
              ! Allocate MPI memory window. Only necessary in the first iteration.
              if (ispin==1) then
                  if (nproc>1) then
                      call mpi_type_size(mpi_double_precision, size_of_double, ierr)
                  else
                      size_of_double = 8
                  end if
                  if (nproc>1 .and. rma_sync==RMA_SYNC_ACTIVE) then
                      comm%window = mpiwindow(n1*n2*n3p(iproc)*comm%nspin, sendbuf(1), bigdft_mpi%mpi_comm)
                  else if (nproc>1 .and. rma_sync==RMA_SYNC_PASSIVE) then
                      call mpi_win_create(sendbuf(1), int(n1*n2*n3p(iproc)*comm%nspin*size_of_double,kind=mpi_address_kind), &
                           size_of_double, MPI_INFO_NULL, bigdft_mpi%mpi_comm, comm%window, ierr)
                  end if
    
                  !!if (rma_sync==RMA_SYNC_ACTIVE) then
                  !!    call mpi_win_fence(mpi_mode_noprecede, comm%window, ierr)
                  !!end if
                  !!if (rma_sync==RMA_SYNC_PASSIVE) then
                  !!    call mpi_win_lock_all(0, comm%window, ierr)
                  !!end if
              end if
              
              do joverlap=1,comm%noverlaps
                  mpisource=comm%comarr(1,joverlap)
                  istsource=comm%comarr(2,joverlap)
                  mpidest=comm%comarr(3,joverlap)
                  istdest=comm%comarr(4,joverlap)
                  nit=comm%comarr(5,joverlap)
                  ioffset_send=comm%comarr(6,joverlap)
                  isend_shift = (ispin-1)*npotarr(mpisource)
                  ! only create the derived data types in the first iteration, otherwise simply reuse them
                  if (ispin==1) then
                      if (nproc>1 .and. type_strategy==TYPES_NESTED) then
                          call mpi_type_create_hvector(nit, 1, int(size_of_double*ioffset_send,kind=mpi_address_kind), &
                               comm%mpi_datatypes(0), comm%mpi_datatypes(joverlap), ierr)
                          call mpi_type_commit(comm%mpi_datatypes(joverlap), ierr)
                      end if
                      iel=0
                      !nsize = 0
                      ii = 1
                      do it=1,nit
                          do i=1,comm%onedtypeovrlp
                              iel = iel + 1
                              displacements(iel) = int(((it-1)*ioffset_send+comm%onedtypearr(1,i))*&
                                                       size_of_double, &
                                                       kind=mpi_address_kind)
                              blocklengths(iel) = comm%onedtypearr(2,i)
                              !nsize = nsize + blocklengths(iel)
                              !if (iproc==0) write(*,*) 'ist, ncount, ind', &
                              !int(displacements(iel)/size_of_double,kind=8)+1, blocklengths(iel), ii
                              ii = ii + blocklengths(iel)
                          end do
                      end do
                      if (nproc>1 .and. type_strategy==TYPES_SIMPLE) then
                          call mpi_type_create_struct(iel, blocklengths, displacements, types, &
                               comm%mpi_datatypes(joverlap), ierr)
                          call mpi_type_commit(comm%mpi_datatypes(joverlap), ierr)
                      end if
                  end if
                  if (iproc==mpidest) then
                      if (.not.bgq) then
                           if (nproc>1) then
                               !ii = nsize
                               call mpi_type_size(comm%mpi_datatypes(joverlap), nsize, ierr)
                               call mpi_type_get_extent(comm%mpi_datatypes(joverlap), lb, extent, ierr)
                               extent=extent/size_of_double
                               nsize=nsize/size_of_double
                               !if (ii/=nsize) then
                               !    write(*,*) 'ispin, ii, nsize', ispin, ii, nsize
                               !    stop 'ii/=nsize'
                               !end if
                           else
                               nsize = 0
                               do i=1,comm%onedtypeovrlp
                                   nsize = nsize + nit*comm%onedtypearr(2,i)
                               end do
                           end if
                           if(nsize>0) then
                               if (nproc>1) then
                                   if (isend_shift+istsource+extent-1 > npotarr(mpisource)*comm%nspin) then
                                       call f_err_throw('out of window: ist='//trim(yaml_toa(isend_shift+istsource,fmt='(i0)'))//&
                                            &', n='//trim(yaml_toa(int(extent,kind=8),fmt='(i0)'))//&
                                            &', nwin='//trim(yaml_toa(npotarr(mpisource)*comm%nspin,fmt='(i0)')),&
                                            err_name='BIGDFT_RUNTIME_ERROR')
                                   end if
                               end if
                               if (nproc> 1 .and. rma_sync==RMA_SYNC_PASSIVE) then
                                   call mpi_win_lock(MPI_LOCK_EXCLUSIVE, mpisource, 0, comm%window, ierr)
                               end if
                               if (nproc>1) then
                                   call mpi_get(recvbuf(ispin_shift+istdest), nsize, &
                                        mpi_double_precision, mpisource, int((isend_shift+istsource-1),kind=mpi_address_kind), &
                                        1, comm%mpi_datatypes(joverlap), comm%window, ierr)
                               else
                                   ind = 0
                                   do i=1,iel
                                      ist = int(displacements(i)/size_of_double,kind=4) + 1
                                      ncount = blocklengths(i)
                                      !!write(*,*) 'joverlap, ist, ncount, ind', &
                                      !!    joverlap, isend_shift+istsource-1+ist, ncount, ispin_shift+istdest+ind, &
                                      call f_memcpy(n=ncount, src=sendbuf(isend_shift+istsource-1+ist), &
                                           dest=recvbuf(ispin_shift+istdest+ind))
                                      ind = ind + ncount
                                   end do
                               end if
                               if (nproc> 1 .and. rma_sync==RMA_SYNC_PASSIVE) then
                                   call mpi_win_unlock(mpisource, comm%window, ierr)
                               end if
                           end if
                       else
                           call mpi_get(recvbuf(ispin_shift+istdest), nit*lzd%glr%d%n1i*lzd%glr%d%n2i, &
                                mpi_double_precision, mpisource, int((isend_shift+istsource-1),kind=mpi_address_kind), &
                                nit*lzd%glr%d%n1i*lzd%glr%d%n2i, mpi_double_precision, comm%window, ierr)
                       end if
                       !!else
                       !!    call mpi_type_size(comm%mpi_datatypes(joverlap), nsize, ierr)
                       !!    nsize=nsize/size_of_double
                       !!    if(nsize>0) then
                       !!        nsize=nsize/nit
                       !!        ist1=ispin_shift+istdest
                       !!        ish1=isend_shift+istsource-1
                       !!        do islices=1,nit
                       !!            ist2=ist1
                       !!            ish2=ish1
                       !!            if (islices<nit) then
                       !!                call mpi_get(recvbuf(ist1), nsize, &
                       !!                     mpi_double_precision, mpisource,&
                       !!                     int(ish1,kind=mpi_address_kind), &
                       !!                     1, comm%mpi_datatypes(0), comm%window,ierr)
                       !!            else
                       !!                do ilines=1,comm%ise(4)-comm%ise(3)+1
                       !!                    write(*,'(5(a,i0))') 'proc ',iproc,' gets ',comm%ise(2)-comm%ise(1)+1, &
                       !!                        ' elements at position ',ist2,' from position ',ish2+1,' on proc ',mpisource
                       !!                    call mpi_get(recvbuf(ist2), comm%ise(2)-comm%ise(1)+1, &
                       !!                         mpi_double_precision, mpisource,&
                       !!                         int(ish2,kind=mpi_address_kind), &
                       !!                         comm%ise(2)-comm%ise(1)+1, mpi_double_precision, comm%window,ierr)
                       !!                    ist2=ist2+comm%ise(2)-comm%ise(1)+1
                       !!                    ish2=ish2+lzd%glr%d%n1i
                       !!                end do
                       !!            end if
                       !!            ist1=ist1+nsize
                       !!            ish1=ish1+ioffset_send
                       !!        end do
                       !!    end if
                       !!end if
                  end if
              end do
    
          end do spin_loop

          if (nproc>1) then
              comm%communication_complete=.false.
          else
              comm%communication_complete=.true.
          end if

      !else nproc_if

      !    !call vcopy(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i, sendbuf(1), 1, recvbuf(1), 1)
      !    !write(*,*) 'lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i*comm%nspin', lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i*comm%nspin
      !    call f_memcpy(n=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i*comm%nspin, src=sendbuf(1), dest=recvbuf(1))
      !    comm%communication_complete=.true.

      !end if nproc_if
      

      call f_free(npotarr)
      call f_free(blocklengths)
      call f_free(displacements)
      call f_free(types)
    
      call timing(iproc, 'Pot_comm start', 'OF')
      call f_release_routine()
    
    end subroutine start_onesided_communication
    
    
    subroutine synchronize_onesided_communication(iproc, nproc, comm)
      use module_base
      use communications_base, only: p2pComms, bgq, RMA_SYNC_ACTIVE, RMA_SYNC_PASSIVE, rma_sync
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc
      type(p2pComms),intent(inout):: comm
      
      ! Local variables
      integer:: ierr, joverlap
      !integer :: mpidest, mpisource
      
      call timing(iproc, 'Pot_comm start', 'ON')      

      if(.not.comm%communication_complete) then
          if (rma_sync==RMA_SYNC_ACTIVE) then
              !!call mpi_win_fence(mpi_mode_nosucceed, comm%window, ierr)
              !!call mpi_win_free(comm%window, ierr)

              call mpi_fenceandfree(comm%window, mpi_mode_nosucceed)
          else if (rma_sync==RMA_SYNC_PASSIVE) then
            !!  !!call mpi_win_unlock_all(comm%window, ierr)
            !!  do joverlap=1,comm%noverlaps
            !!      mpisource=comm%comarr(1,joverlap)
            !!      mpidest=comm%comarr(3,joverlap)
            !!      if (iproc==mpidest) then
            !!          if (.not.bgq) then
            !!               call mpi_win_unlock(mpisource, comm%window, ierr)
            !!           end if
            !!      end if
            !!  end do
            call mpi_win_free(comm%window, ierr)
          end if
          do joverlap=1,comm%noverlaps
              call mpi_type_free(comm%mpi_datatypes(joverlap), ierr)
          end do
      end if
    
      ! Flag indicating that the communication is complete
      comm%communication_complete=.true.
          
      call timing(iproc, 'Pot_comm start', 'OF')      

    end subroutine synchronize_onesided_communication


    !> Locreg communication
    subroutine communicate_locreg_descriptors_basics(iproc, nproc, nlr, rootarr, orbs, llr)
      use module_base
      use module_types, only: orbitals_data, locreg_descriptors
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nlr
      integer,dimension(nlr),intent(in) :: rootarr
      type(orbitals_data),intent(in) :: orbs
      type(locreg_descriptors),dimension(nlr),intent(inout) :: llr
    
      ! Local variables
      integer:: ierr, ilr, iilr, jproc
      ! integer:: istat, iall
      character(len=1),dimension(:),allocatable :: worksend_char, workrecv_char
      logical,dimension(:),allocatable :: worksend_log, workrecv_log
      integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
      integer,dimension(:),allocatable :: norb_par
      real(8),dimension(:,:),allocatable :: worksend_dbl, workrecv_dbl
      character(len=*),parameter :: subname='communicate_locreg_descriptors_basics'

      call f_routine(id=subname)
    
!!$      allocate(worksend_char(orbs%norbp), stat=istat)
!!$      call memocc(istat, worksend_char, 'worksend_char', subname)
      worksend_char= f_malloc_str(len(worksend_char),orbs%norbp,&
           id='worksend_char')
      worksend_log = f_malloc(orbs%norbp,id='worksend_log')
      worksend_int = f_malloc((/ 27, orbs%norbp /),id='worksend_int')
      worksend_dbl = f_malloc((/ 6, orbs%norbp /),id='worksend_dbl')
    
      workrecv_char= f_malloc_str(len(workrecv_char),orbs%norb,&
           id='workrecv_char')
!!$      allocate(workrecv_char(orbs%norb), stat=istat)
!!$      call memocc(istat, workrecv_char, 'workrecv_char', subname)
      workrecv_log = f_malloc(orbs%norb,id='workrecv_log')
      workrecv_int = f_malloc((/ 27, orbs%norb /),id='workrecv_int')
      workrecv_dbl = f_malloc((/ 6, orbs%norb /),id='workrecv_dbl')
    
    
      iilr=0
      do ilr=1,nlr
          if (iproc==rootarr(ilr)) then
              iilr=iilr+1
              worksend_char(iilr)=llr(ilr)%geocode
              worksend_log(iilr)=llr(ilr)%hybrid_on
              worksend_int(1,iilr)=ilr
              worksend_int(2,iilr)=llr(ilr)%ns1
              worksend_int(3,iilr)=llr(ilr)%ns2
              worksend_int(4,iilr)=llr(ilr)%ns3
              worksend_int(5,iilr)=llr(ilr)%nsi1
              worksend_int(6,iilr)=llr(ilr)%nsi2
              worksend_int(7,iilr)=llr(ilr)%nsi3
              worksend_int(8,iilr)=llr(ilr)%localnorb
              worksend_int(9:11,iilr)=llr(ilr)%outofzone(1:3)
              worksend_int(12,iilr)=llr(ilr)%wfd%nvctr_c
              worksend_int(13,iilr)=llr(ilr)%wfd%nvctr_f
              worksend_int(14,iilr)=llr(ilr)%wfd%nseg_c
              worksend_int(15,iilr)=llr(ilr)%wfd%nseg_f
              worksend_int(16,iilr)=llr(ilr)%d%n1
              worksend_int(17,iilr)=llr(ilr)%d%n2
              worksend_int(18,iilr)=llr(ilr)%d%n3
              worksend_int(19,iilr)=llr(ilr)%d%nfl1
              worksend_int(20,iilr)=llr(ilr)%d%nfu1
              worksend_int(21,iilr)=llr(ilr)%d%nfl2
              worksend_int(22,iilr)=llr(ilr)%d%nfu2
              worksend_int(23,iilr)=llr(ilr)%d%nfl3
              worksend_int(24,iilr)=llr(ilr)%d%nfu3
              worksend_int(25,iilr)=llr(ilr)%d%n1i
              worksend_int(26,iilr)=llr(ilr)%d%n2i
              worksend_int(27,iilr)=llr(ilr)%d%n3i
              worksend_dbl(1:3,iilr)=llr(ilr)%locregCenter(1:3)
              worksend_dbl(4,iilr)=llr(ilr)%locrad
              worksend_dbl(5,iilr)=llr(ilr)%locrad_kernel
              worksend_dbl(6,iilr)=llr(ilr)%locrad_mult
          end if
      end do

      norb_par = f_malloc(0.to.nproc-1,id='norb_par')
      do jproc=0,nproc-1
          norb_par(jproc) = orbs%norb_par(jproc,0)
      end do
    
      call mpi_allgatherv(worksend_char, orbs%norbp, mpi_character, workrecv_char, norb_par, &
           orbs%isorb_par, mpi_character, bigdft_mpi%mpi_comm, ierr)
      call mpi_allgatherv(worksend_log, orbs%norbp, mpi_logical, workrecv_log, norb_par, &
           orbs%isorb_par, mpi_logical, bigdft_mpi%mpi_comm, ierr)
      do jproc=0,nproc-1
          norb_par(jproc) = 27*orbs%norb_par(jproc,0)
      end do
      call mpi_allgatherv(worksend_int, 27*orbs%norbp, mpi_integer, workrecv_int, norb_par, &
           27*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
!!$      call mpiallgather(sendbuf=worksend_int,sendcount=27*orbs%norbp,recvbuf=workrecv_int,recvcounts=27*orbs%norb_par(:,0), &
!!$           displs=27*orbs%isorb_par, comm=bigdft_mpi%mpi_comm)
      do jproc=0,nproc-1
          norb_par(jproc) = 6*orbs%norb_par(jproc,0)
      end do
      call mpi_allgatherv(worksend_dbl, 6*orbs%norbp, mpi_double_precision, workrecv_dbl, norb_par, &
           6*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)

       call f_free(norb_par)
    
      do ilr=1,nlr
          iilr=workrecv_int(1,ilr)
          llr(iilr)%geocode=workrecv_char(ilr)
          llr(iilr)%hybrid_on= workrecv_log(ilr)
          llr(iilr)%ns1=workrecv_int(2,ilr)
          llr(iilr)%ns2=workrecv_int(3,ilr)
          llr(iilr)%ns3=workrecv_int(4,ilr)
          llr(iilr)%nsi1=workrecv_int(5,ilr)
          llr(iilr)%nsi2=workrecv_int(6,ilr)
          llr(iilr)%nsi3=workrecv_int(7,ilr)
          llr(iilr)%localnorb=workrecv_int(8,ilr)
          llr(iilr)%outofzone(1:3)=workrecv_int(9:11,ilr)
          llr(iilr)%wfd%nvctr_c=workrecv_int(12,ilr)
          llr(iilr)%wfd%nvctr_f=workrecv_int(13,ilr)
          llr(iilr)%wfd%nseg_c=workrecv_int(14,ilr)
          llr(iilr)%wfd%nseg_f=workrecv_int(15,ilr)
          llr(iilr)%d%n1=workrecv_int(16,ilr)
          llr(iilr)%d%n2=workrecv_int(17,ilr)
          llr(iilr)%d%n3=workrecv_int(18,ilr)
          llr(iilr)%d%nfl1=workrecv_int(19,ilr)
          llr(iilr)%d%nfu1=workrecv_int(20,ilr)
          llr(iilr)%d%nfl2=workrecv_int(21,ilr)
          llr(iilr)%d%nfu2=workrecv_int(22,ilr)
          llr(iilr)%d%nfl3=workrecv_int(23,ilr)
          llr(iilr)%d%nfu3=workrecv_int(24,ilr)
          llr(iilr)%d%n1i=workrecv_int(25,ilr)
          llr(iilr)%d%n2i=workrecv_int(26,ilr)
          llr(iilr)%d%n3i=workrecv_int(27,ilr)
          llr(iilr)%locregCenter(1:3)=workrecv_dbl(1:3,ilr)
          llr(iilr)%locrad=workrecv_dbl(4,ilr)
          llr(iilr)%locrad_kernel=workrecv_dbl(5,ilr)
          llr(iilr)%locrad_mult=workrecv_dbl(6,ilr)
      end do
    
    
      call f_free(worksend_int)
      call f_free(workrecv_int)
      !!worksend_int = f_malloc((/ 13, orbs%norbp /),id='worksend_int')
      !!workrecv_int = f_malloc((/ 13, orbs%norb /),id='workrecv_int')
    
    
      !!iilr=0
      !!do ilr=1,nlr
      !!    if (iproc==rootarr(ilr)) then
      !!        iilr=iilr+1
      !!        worksend_int(1,iilr)=llr(ilr)%d%n1
      !!        worksend_int(2,iilr)=llr(ilr)%d%n2
      !!        worksend_int(3,iilr)=llr(ilr)%d%n3
      !!        worksend_int(4,iilr)=llr(ilr)%d%nfl1
      !!        worksend_int(5,iilr)=llr(ilr)%d%nfu1
      !!        worksend_int(6,iilr)=llr(ilr)%d%nfl2
      !!        worksend_int(7,iilr)=llr(ilr)%d%nfu2
      !!        worksend_int(8,iilr)=llr(ilr)%d%nfl3
      !!        worksend_int(9,iilr)=llr(ilr)%d%nfu3
      !!        worksend_int(10,iilr)=llr(ilr)%d%n1i
      !!        worksend_int(11,iilr)=llr(ilr)%d%n2i
      !!        worksend_int(12,iilr)=llr(ilr)%d%n3i
      !!        worksend_int(13,iilr)=ilr
      !!    end if
      !!end do
    
      !!call mpi_allgatherv(worksend_int, 13*orbs%norbp, mpi_integer, workrecv_int, 13*orbs%norb_par(:,0), &
      !!     13*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    
      !!do ilr=1,nlr
      !!    iilr=workrecv_int(13,ilr)
      !!    llr(iilr)%d%n1=workrecv_int(1,ilr)
      !!    llr(iilr)%d%n2=workrecv_int(2,ilr)
      !!    llr(iilr)%d%n3=workrecv_int(3,ilr)
      !!    llr(iilr)%d%nfl1=workrecv_int(4,ilr)
      !!    llr(iilr)%d%nfu1=workrecv_int(5,ilr)
      !!    llr(iilr)%d%nfl2=workrecv_int(6,ilr)
      !!    llr(iilr)%d%nfu2=workrecv_int(7,ilr)
      !!    llr(iilr)%d%nfl3=workrecv_int(8,ilr)
      !!    llr(iilr)%d%nfu3=workrecv_int(9,ilr)
      !!    llr(iilr)%d%n1i=workrecv_int(10,ilr)
      !!    llr(iilr)%d%n2i=workrecv_int(11,ilr)
      !!    llr(iilr)%d%n3i=workrecv_int(12,ilr)
      !!end do
    
    
!!$      iall=-product(shape(worksend_char))*kind(worksend_char)
!!$      deallocate(worksend_char,stat=istat)
!!$      call memocc(istat, iall, 'worksend_char', subname)
      call f_free_str(len(worksend_char),worksend_char)
      call f_free(worksend_log)
      !!call f_free(worksend_int)
      call f_free(worksend_dbl)

!!$      iall=-product(shape(workrecv_char))*kind(workrecv_char)
!!$      deallocate(workrecv_char,stat=istat)
!!$      call memocc(istat, iall, 'workrecv_char', subname)
      call f_free_str(len(workrecv_char),workrecv_char)
      call f_free(workrecv_log)
      !!call f_free(workrecv_int)
      call f_free(workrecv_dbl)

      call f_release_routine()
    
    end subroutine communicate_locreg_descriptors_basics
    
    
    subroutine communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, rootarr, onwhichmpi)
       use module_base
       use module_types, only: orbitals_data, locreg_descriptors
       use locregs, only: allocate_wfd, check_overlap_cubic_periodic
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
       integer :: ierr, jorb, ilr, jlr, root, max_sim_comms, norb_max
       !integer :: icomm, ilr_old, jtask, nalloc, nrecv
       integer :: maxrecvdim, maxsenddim, ioffset, window, ist_dest, ist_source, ithread, nthread, mthread
       integer :: iorb, jjorb, ncount, iiorb, size_of_int, info, jproc, ii, iilr, ncover, isproc, nprocp
       logical :: isoverlap
       character(len=*), parameter:: subname='communicate_locreg_descriptors_keys'
       !integer,dimension(:),allocatable :: requests
       !integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
       integer,dimension(:),allocatable :: worksend, workrecv, ncomm, nrecvcounts, types, derived_types!, tmparr_int, istarr
       integer,dimension(:),allocatable :: cover_id, ncount_par, iscount_par
       integer,dimension(:,:),allocatable :: locregs_local
       integer,dimension(:,:),allocatable :: blocklengths
       integer(kind=mpi_address_kind),dimension(:,:),allocatable :: displacements
       logical,dimension(:),allocatable :: covered
       logical,dimension(:,:),allocatable :: mask
       !$ integer :: omp_get_thread_num, omp_get_max_threads

       call f_routine(id=subname)

       ithread = 0
       nthread = 1
       !$ nthread = omp_get_max_threads()

       !@ NEW VESRION #############################################
       ! should be 1D later...
       covered = f_malloc(nlr,id='covered')
       norb_max = maxval(orbs%norb_par(:,0))
       ! In principle bad practice to allocate with nthread, but as norb_max is usually small, it should be ok.
       locregs_local = f_malloc((/1.to.norb_max,0.to.nthread/),id='locregs_local')
       mask = f_malloc((/1.to.norb_max,0.to.nthread/),id='mask')
       ncount_par = f_malloc0(0.to.nproc-1,id='ncount_par')
       iscount_par = f_malloc0(0.to.nproc-1,id='iscount_par')

       do iorb=1,orbs%norbp
           iiorb = orbs%isorb + iorb
           ilr = orbs%inwhichlocreg(iiorb)
           locregs_local(iorb,0:nthread) = ilr
       end do

       ! Determine which locregs process iproc should get.
       ncover = 0
       !$omp parallel default(none) &
       !$omp shared(ncover, nlr, rootarr, covered, orbs, glr, llr, iproc) &
       !$omp private(ilr, root, jorb, jjorb, jlr, isoverlap)
       !$omp do reduction(+: ncover)
       do ilr=1,nlr
           root=rootarr(ilr)
           covered(ilr)=.false.
           do jorb=1,orbs%norbp
               jjorb=orbs%isorb+jorb
               jlr=orbs%inwhichlocreg(jjorb)
               ! don't communicate to ourselves, or if we've already sent this locreg
               if (iproc == root .or. covered(ilr)) cycle
               call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
               if (isoverlap) then         
                   covered(ilr)=.true.
               end if
           end do
           ! For spin polarized calculations, norbup and norbp are different,
           ! therefore one also has to check this distribution.
           do jorb=1,orbs%norbup
               jjorb=orbs%isorbu+jorb
               jlr=orbs%inwhichlocreg(jjorb)
               ! don't communicate to ourselves, or if we've already sent this locreg
               if (iproc == root .or. covered(ilr)) cycle
               call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
               if (isoverlap) then         
                   covered(ilr)=.true.
               end if
           end do
           if (covered(ilr)) then
               ncover = ncover + 1
           end if
       end do
       !$omp end do
       !$omp end parallel

       cover_id = f_malloc(ncover,id='cover_id')
       ii = 0
       do ilr=1,nlr
           if (covered(ilr)) then
               ii = ii + 1
               cover_id(ii) = ilr
           end if
       end do
       if (ii/=ncover) call f_err_throw('ii/=ncover')

       call f_free(covered)


       ! Each process makes its data available in a contiguous workarray.
       maxsenddim=0
       do iorb=1,orbs%norbp
           iiorb=orbs%isorb+iorb
           ilr=orbs%inwhichlocreg(iiorb)
           maxsenddim = maxsenddim + 6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
       end do
       worksend = f_malloc(max(maxsenddim,1),id='worksend_int')

       ! Copy the keys to a contiguous array, ordered by the locreg ID (which
       ! is not necessarily identical to the orbital ID order).
       mask(1:orbs%norbp,0) = .true.
       ioffset=0
       !do iilr=1,nlr
           do iorb=1,orbs%norbp
               !iiorb=orbs%isorb+iorb
               !ilr=orbs%inwhichlocreg(iiorb)
               iiorb = minloc(locregs_local(1:orbs%norbp,0), 1, mask(1:orbs%norbp,0))
               ilr=locregs_local(iiorb,0)
               mask(iiorb,0) = .false.
               !if (ilr==iilr) then
                   !write(*,*) 'COPY: iproc, ilr', iproc, ilr
                   ncount=llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f
                   call vcopy(2*ncount, llr(ilr)%wfd%keygloc(1,1), 1, worksend(ioffset+1), 1)
                   call vcopy(2*ncount, llr(ilr)%wfd%keyglob(1,1), 1, worksend(ioffset+2*ncount+1), 1)
                   call vcopy(ncount, llr(ilr)%wfd%keyvloc(1), 1, worksend(ioffset+4*ncount+1), 1)
                   call vcopy(ncount, llr(ilr)%wfd%keyvglob(1), 1, worksend(ioffset+5*ncount+1), 1)
                   ioffset=ioffset+6*ncount
               !end if
           end do
       !end do

       ! Initialize the MPI window
       !!call mpi_type_size(mpi_integer, size_of_int, ierr)
       !!call mpi_info_create(info, ierr)
       !!call mpi_info_set(info, "no_locks", "true", ierr)
       !!call mpi_win_create(worksend(1), int(maxsenddim*size_of_int,kind=mpi_address_kind), size_of_int, &
       !!     info, bigdft_mpi%mpi_comm, window, ierr)
       !!call mpi_info_free(info, ierr)
       !!call mpi_win_fence(mpi_mode_noprecede, window, ierr)
       call mpi_type_size(mpi_integer, size_of_int, ierr)
       window =  mpiwindow(maxsenddim, worksend(1),  bigdft_mpi%mpi_comm)

       ! Allocate the receive buffer
       maxrecvdim=0
       !do ilr=1,nlr
       do ii=1,ncover
           ilr = cover_id(ii)
           root=rootarr(ilr)
           !if (covered(ilr,iproc)) then
               ncount=6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
               maxrecvdim=maxrecvdim+ncount
           !end if
       end do
       workrecv = f_malloc(maxrecvdim,id='workrecv')

       ! Do the communication. Communicate only once for each pair of MPI tasks using a derived datatype.
       ncomm = f_malloc0(0.to.nproc-1,id='ncomm')
       nrecvcounts = f_malloc0(0.to.nproc-1,id='nrecvcounts')
       blocklengths = f_malloc0((/1.to.maxval(orbs%norb_par(:,0)),0.to.nproc-1/),id='blocklengths')
       displacements = f_malloc0((/1.to.maxval(orbs%norb_par(:,0)),0.to.nproc-1/),id='displacements')

       nthread = 1
       ithread = 0
       !$ nthread = omp_get_max_threads()
       nprocp = nproc/nthread
       mthread = nproc-nthread*nprocp
       !$omp parallel default(none) &
       !$omp shared(mthread, ncover, cover_id, rootarr, ncomm, llr, blocklengths, nlr, norb_max, orbs, locregs_local, mask) &
       !$omp shared(displacements, nrecvcounts, ncount_par, size_of_int) &
       !$omp private(isproc, ilr, root, ii, ncount, ist_source) &
       !$omp firstprivate(ithread, nprocp)
       !$ ithread = omp_get_thread_num()
       isproc = ithread*nprocp + min(mthread,ithread)
       if (ithread<mthread) nprocp = nprocp + 1
       do ii=1,ncover
           ilr = cover_id(ii)
           root=rootarr(ilr)
           if (root>=isproc .and. root<isproc+nprocp) then
               ncomm(root) = ncomm(root) + 1
               ncount=6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
               ncount_par(root) = ncount_par(root) + ncount
               ist_source=get_offset(nlr, llr, orbs, norb_max, root, ilr, locregs_local(:,ithread), mask(:,ithread))
               blocklengths(ncomm(root),root) = ncount
               displacements(ncomm(root),root) = int(ist_source*size_of_int,kind=mpi_address_kind)
               nrecvcounts(root) = nrecvcounts(root) + ncount
           end if
       end do
       !$omp end parallel


       iscount_par(0) = 0
       do jproc=1,nproc-1
           iscount_par(jproc) = iscount_par(jproc-1) + ncount_par(jproc-1)
       end do

       !!! The values in displacements must be sorted...
       !!!tmparr_long = f_malloc(maxval(ncomm(:)),id='tmparr')
       !!!tmparr_int = f_malloc(maxval(ncomm(:)),id='tmparr')
       !!!istarr = f_malloc(maxval(ncomm(:)),id='istarr')
       !!do jproc=0,nproc-1
       !!    !tmparr_long = displacements(1:ncomm(jproc),jproc)
       !!    !tmparr_int = blocklengths(1:ncomm(jproc),jproc)
       !!    ii = 0
       !!    do icomm=1,ncomm(jproc)
       !!        iloc = minloc(tmparr_long(1:ncomm(jproc),1)
       !!        displacements(icomm,jproc) = tmparr_long(iloc)
       !!        blocklengths(icomm,jproc) = tmparr_int(iloc)
       !!        !istarr(icomm) = ii
       !!        ii = ii + tmparr_int(iloc)
       !!    end do
       !!end do

       derived_types = f_malloc(0.to.nproc-1,id='derived_types')
       types = f_malloc(maxval(ncomm),id='types')
       types(:) = mpi_integer
       ist_dest=1
       do jproc=0,nproc-1
           !write(*,'(a,11i8)') 'iproc, jproc, ist_dest, blocklengths(1:ncomm(jproc)), displacements(1:ncomm(jproc))', iproc, jproc, ist_dest, blocklengths(1:ncomm(jproc),jproc), displacements(1:ncomm(jproc),jproc)
           call mpi_type_create_struct(ncomm(jproc), blocklengths(1,jproc), displacements(1,jproc), &
                types, derived_types(jproc), ierr)
           call mpi_type_commit(derived_types(jproc), ierr)
           !call mpi_type_size(types(1), ii, ierr)
           !write(*,*) 'size types', ii
           call mpi_type_size(derived_types(jproc), ii, ierr)
           if (ii/=nrecvcounts(jproc)*size_of_int) then
               !write(*,'(a,5i8)') 'iproc, jproc, ncomm(jproc), ii, nrecvcounts(jproc)*size_of_int', iproc, jproc, ncomm(jproc), ii, nrecvcounts(jproc)*size_of_int
               !write(*,'(a,10i8)') 'iproc, jproc, blocklengths(1:ncomm(jproc)), displacements(1:ncomm(jproc))', iproc, jproc, blocklengths(1:ncomm(jproc),jproc), displacements(1:ncomm(jproc),jproc)
               call f_err_throw('wrong size of derived_types(jproc)',err_name='BIGDFT_RUNTIME_ERROR')
           end if
           if (nrecvcounts(jproc)>0) then
               call mpi_get(workrecv(ist_dest), nrecvcounts(jproc), mpi_integer, jproc, &
                    int(0,kind=mpi_address_kind), 1, derived_types(jproc), window, ierr)
           end if
           ist_dest = ist_dest + nrecvcounts(jproc)
       end do

           !    call mpi_get(workrecv(ist_dest), ncount, mpi_integer, root, &
           !         int(ist_source,kind=mpi_address_kind), ncount, mpi_integer, window, ierr)
           !    ist_dest=ist_dest+ncount
           !end if
       !end do

       ! Synchronize the communication
       call mpi_win_fence(mpi_mode_nosucceed, window, ierr)
       do jproc=0,nproc-1
           call mpi_type_free(derived_types(jproc), ierr)
       end do
       call mpi_win_free(window, ierr)
       call f_free(derived_types)

       call f_free(ncomm)
       call f_free(nrecvcounts)
       call f_free(blocklengths)
       call f_free(displacements)
       call f_free(types)

       !do ii=1,size(workrecv)
       !    write(100+iproc,*) ii, workrecv(ii)
       !end do


       ! Copy the date from the workarrays to the correct locations
       call f_free(worksend)
       !ist_dest=0
       !do jproc=0,nproc-1
           !do ilr=1,nlr
           !    if (.not.covered(ilr,iproc)) cycle
           do ii=1,ncover
               ilr = cover_id(ii)
               root = rootarr(ilr)
               ist_dest = iscount_par(root)
               !ist_dest = 0
               !do jproc=0,root-1
               !    ist_dest = ist_dest + ncount_par(jproc)
               !end do
               !do jorb=1,orbs%norb_par(root,0)
               !    jjorb = orbs%isorb_par(jproc) + jorb
               !    jlr = orbs%inwhichlocreg(jjorb)
               !    if (ilr==jlr) then
                       call allocate_wfd(llr(ilr)%wfd)
                       ncount=llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f
                       call vcopy(2*ncount, workrecv(ist_dest+1), 1, llr(ilr)%wfd%keygloc(1,1), 1)
                       call vcopy(2*ncount, workrecv(ist_dest+2*ncount+1), 1, llr(ilr)%wfd%keyglob(1,1), 1)
                       call vcopy(ncount, workrecv(ist_dest+4*ncount+1), 1, llr(ilr)%wfd%keyvloc(1), 1)
                       call vcopy(ncount, workrecv(ist_dest+5*ncount+1), 1, llr(ilr)%wfd%keyvglob(1), 1)
                       !ist_dest=ist_dest+6*ncount
               !    end if
               !end do
               iscount_par(root) = iscount_par(root) + 6*ncount
           end do
       !end do
       call f_free(workrecv)
       !call f_free(covered)
       call f_free(cover_id)
       call f_free(ncount_par)
       call f_free(iscount_par)

       call f_free(locregs_local)
       call f_free(mask)

       !@ END NEW VESRION ##########################################


    
!!!       ! This maxval is put out of the allocate to avoid compiler crash with PathScale.
!!!       jorb = maxval(orbs%norb_par(:,0))
!!!       requests = f_malloc(8*nproc*jorb,id='requests')
!!!       covered = f_malloc((/ 1.to.nlr, 0.to.nproc-1 /),id='covered')
!!!       worksend_int = f_malloc((/ 4, nlr /),id='worksend_int')
!!!       workrecv_int = f_malloc((/ 4, nlr /),id='workrecv_int')
!!!    
!!!       ! divide communications into chunks to avoid problems with memory (too many communications)
!!!       ! set maximum number of simultaneous communications.
!!!       ! SM: set this value such tag the MPI tag is never greater than 4000000 (otherwise crash on Cray... is the limit maybe 2^22?)
!!!       max_sim_comms=min(nlr,int(4000000.d0/real(nproc,kind=8)))
!!!       !max_sim_comms=min(nlr,2)
!!!
!!!
!!!       nrecv=0
!!!       !nsend=0
!!!       icomm=0
!!!       maxsenddim=0
!!!       do ilr=1,nlr
!!!           root=rootarr(ilr)
!!!           covered(ilr,:)=.false.
!!!           do jorb=1,orbs%norb
!!!               jlr=orbs%inwhichlocreg(jorb)
!!!               jtask=onwhichmpi(jorb)
!!!               ! check we're on a sending or receiving proc
!!!               if (iproc /= root .and. iproc /= jtask) cycle
!!!               ! don't communicate to ourselves, or if we've already sent this locreg
!!!               if (jtask == root .or. covered(ilr,jtask)) cycle
!!!               call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
!!!               if (isoverlap) then         
!!!                   covered(ilr,jtask)=.true.
!!!                   if (iproc == root) then
!!!                   !!   !write(*,'(5(a,i0))') 'process ',iproc,' sends locreg ',ilr,' to process ',&
!!!                   !!   !    jtask,' with tags ',4*ilr+0,'-',4*ilr+3
!!!                   !!   worksend_int(1,ilr)=llr(ilr)%wfd%nvctr_c
!!!                   !!   worksend_int(2,ilr)=llr(ilr)%wfd%nvctr_f
!!!                   !!   worksend_int(3,ilr)=llr(ilr)%wfd%nseg_c
!!!                   !!   worksend_int(4,ilr)=llr(ilr)%wfd%nseg_f
!!!                   !!   icomm=icomm+1
!!!                   !!   call mpi_isend(worksend_int(1,ilr), 4, mpi_integer, jtask,&
!!!                   !!        itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
!!!                      maxsenddim=max(maxsenddim,llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
!!!                   !!   !nsend=nsend+1
!!!                   else if (iproc == jtask) then
!!!                   !!   !write(*,'(5(a,i0))') 'process ',iproc,' receives locreg ',ilr,' from process ',&
!!!                   !!   !    root,' with tags ',4*ilr+0,'-',4*ilr+3
!!!                   !!   icomm=icomm+1
!!!                   !!   call mpi_irecv(workrecv_int(1,ilr), 4, mpi_integer, root,&
!!!                   !!        itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
!!!                      nrecv=nrecv+1
!!!                   end if
!!!               end if
!!!           end do
!!!           !!if (mod(ilr,max_sim_comms)==0 .or. ilr==nlr) then
!!!           !!   call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
!!!           !!   if (f_err_raise(ierr /= 0,'problem in communicate locregs: error in mpi_waitall '//&
!!!           !!        trim(yaml_toa(ierr))//' for process '//trim(yaml_toa(iproc)),&
!!!           !!        err_name='BIGDFT_RUNTIME_ERROR')) return
!!!           !!   call mpi_barrier(mpi_comm_world,ierr)
!!!           !!   icomm=0
!!!           !!end if
!!!       end do
!!!      
!!!       !!call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
!!!       !!call mpi_barrier(mpi_comm_world,ierr)
!!!    
!!!       call f_free(worksend_int)
!!!    
!!!       nalloc=0
!!!       maxrecvdim=0
!!!       do jlr=1,nlr 
!!!          !write(*,*) 'iproc, jlr, covered, nseg_c', iproc, jlr, covered(jlr,iproc), llr(jlr)%wfd%nseg_c
!!!          if (covered(jlr,iproc)) then
!!!             !!llr(jlr)%wfd%nvctr_c=workrecv_int(1,jlr)
!!!             !!llr(jlr)%wfd%nvctr_f=workrecv_int(2,jlr)
!!!             !!llr(jlr)%wfd%nseg_c=workrecv_int(3,jlr)
!!!             !!llr(jlr)%wfd%nseg_f=workrecv_int(4,jlr)
!!!    !         call allocate_wfd(llr(jlr)%wfd,subname)
!!!             nalloc=nalloc+1
!!!             maxrecvdim=max(maxrecvdim,llr(jlr)%wfd%nseg_c+llr(jlr)%wfd%nseg_f)
!!!          end if
!!!       end do
!!!       if (f_err_raise(nalloc /= nrecv,'problem in communicate locregs: mismatch in receives '//&
!!!            trim(yaml_toa(nrecv))//' and allocates '//trim(yaml_toa(nalloc))//' for process '//trim(yaml_toa(iproc)),&
!!!            err_name='BIGDFT_RUNTIME_ERROR')) return
!!!    
!!!       call f_free(workrecv_int)
!!!    
!!!       !should reduce memory by not allocating for all llr
!!!       workrecv_int = f_malloc((/ 6*maxrecvdim, nlr /),id='workrecv_int')
!!!       worksend_int = f_malloc((/ 6*maxsenddim, nlr /),id='worksend_int')
!!!    
!!!       !!! divide communications into chunks to avoid problems with memory (too many communications)
!!!       !!! set maximum number of simultaneous communications
!!!       !!!total_sent=0
!!!       !!!total_recv=0
!!!       !!max_sim_comms=1000
!!!       icomm=0
!!!       ilr_old=0
!!!       do ilr=1,nlr
!!!          root=rootarr(ilr)
!!!          do jtask=0,nproc-1
!!!             if (.not. covered(ilr,jtask)) cycle
!!!             if (iproc == root) then
!!!               !write(*,'(5(a,i0))') 'process ',iproc,' sends locreg ',ilr,' to process ',&
!!!               !     jtask,' with tags ',4*ilr+0,'-',4*ilr+3
!!!               call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyglob(1,1),1,worksend_int(1,ilr),1)
!!!               call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keygloc(1,1),1,&
!!!                    worksend_int(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
!!!               call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyvloc(1),1,&
!!!                    worksend_int(4*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
!!!               call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyvglob(1),1,&
!!!                    worksend_int(5*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
!!!               icomm=icomm+1
!!!               call mpi_isend(worksend_int(1,ilr),6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
!!!                    jtask, itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
!!!             else if (iproc == jtask) then
!!!                !write(*,'(5(a,i0))') 'process ',iproc,' receives locreg ',ilr,' from process ',&
!!!                !    root,' with tags ',4*ilr+0,'-',4*ilr+3
!!!                icomm=icomm+1
!!!                call mpi_irecv(workrecv_int(1,ilr),6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
!!!                     root, itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
!!!             end if
!!!          end do
!!!          if (mod(ilr,max_sim_comms)==0 .or. ilr==nlr) then
!!!             !do jlr=max(ilr-max_sim_comms+1,1),ilr
!!!             do jlr=ilr_old+1,ilr
!!!                !write(*,'(2(a,i0))') 'process ',iproc,' allocates locreg ',jlr
!!!                if (covered(jlr,iproc))  call allocate_wfd(llr(jlr)%wfd)
!!!             end do
!!!             call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
!!!             if (f_err_raise(ierr /= 0,'problem in communicate locregs: error in mpi_waitall '//&
!!!                  trim(yaml_toa(ierr))//' for process '//trim(yaml_toa(iproc)),&
!!!                  err_name='BIGDFT_RUNTIME_ERROR')) return
!!!             call mpi_barrier(mpi_comm_world,ierr)
!!!             icomm=0
!!!             ilr_old=ilr
!!!          end if
!!!       end do
!!!    
!!!       call f_free(worksend_int)
!!!    
!!!       do ilr=1,nlr 
!!!          if (covered(ilr,iproc)) then
!!!             call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),workrecv_int(1,ilr),1,llr(ilr)%wfd%keyglob(1,1),1)
!!!             call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
!!!                  workrecv_int(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keygloc(1,1),1)
!!!             call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
!!!                  workrecv_int(4*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keyvloc(1),1)
!!!             call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
!!!                  workrecv_int(5*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keyvglob(1),1)
!!!          end if
!!!       end do
!!!    
!!!       call f_free(workrecv_int)
!!!    
!!!       !print*,'iproc,sent,received,num sent,num received',iproc,total_sent,total_recv,nsend,nrecv
!!!       call f_free(requests)
!!!       call f_free(covered)

       call f_release_routine()
    
    contains
    
     pure function itag(ilr,recv)
     implicit none
     integer, intent(in) :: ilr,recv
     integer :: itag
    
     !itag=ilr+recv*nlr
     !itag=ilr+recv*max_sim_comms
     itag = mod(ilr-1,max_sim_comms)*nproc + recv + 1
    
     end function itag

    
    END SUBROUTINE communicate_locreg_descriptors_keys


    !> Get the offset of the data of locreg iilr
    function get_offset(nlr, llr, orbs, norb_max, iiproc, iilr, locregs_local, mask)
      use module_base
      use module_types, only: orbitals_data, locreg_descriptors
      implicit none
      integer,intent(in) :: nlr, norb_max, iiproc, iilr
      type(locreg_descriptors),dimension(nlr),intent(inout) :: llr
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(norb_max),intent(inout) :: locregs_local
      logical,dimension(norb_max),intent(inout) :: mask
      integer :: get_offset

      ! Local variables
      integer :: ilr, jorb, jjorb, jlr, ncount, iiorb, iorb

      !call f_routine(id='get_offset')

      !!mask = f_malloc(norb_max,id='norb_max')
      !!locregs_local = f_malloc(norb_max,id='locregs_local')

      !! NEW #########################################
      do iorb=1,orbs%norb_par(iiproc,0)
          iiorb = orbs%isorb_par(iiproc) + iorb
          ilr = orbs%inwhichlocreg(iiorb)
          locregs_local(iorb) = ilr
      end do

      !mask(1:orbs%norb_par(iiproc,0)) = .true.
      mask(:) = .true.
      !!write(*,*) 'iproc, orbs%norb_par(iiproc,0), locregs_local(1:orbs%norb_par(iiproc,0)), mask(1:orbs%norb_par(iiproc,0))', &
      !!            bigdft_mpi%iproc, orbs%norb_par(iiproc,0), locregs_local(1:orbs%norb_par(iiproc,0)), &
      !!            mask(1:orbs%norb_par(iiproc,0))
      get_offset = 0
      do iorb=1,orbs%norb_par(iiproc,0)
          iiorb = minloc(locregs_local(1:orbs%norb_par(iiproc,0)), 1, mask(1:orbs%norb_par(iiproc,0)))
          !!write(*,*) 'iproc, iiorb, locregs_local(1:orbs%norb_par(iiproc,0)), mask(1:orbs%norb_par(iiproc,0))', &
          !!            bigdft_mpi%iproc, iiorb, locregs_local(1:orbs%norb_par(iiproc,0)), mask(1:orbs%norb_par(iiproc,0))
          ilr=locregs_local(iiorb)
          mask(iiorb) = .false.
          if (ilr==iilr) exit !locreg found
          ncount=6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
          get_offset=get_offset+ncount
      end do

      !!call f_free(mask)
      !!call f_free(locregs_local)

      !! OLD #########################################
      !!get_offset=0
      !!ilr_loop: do ilr=1,nlr
      !!    do jorb=1,orbs%norb_par(iiproc,0)
      !!        jjorb=orbs%isorb_par(iiproc)+jorb
      !!        jlr=orbs%inwhichlocreg(jjorb)
      !!        !if (iproc==0) write(*,*) 'ilr, jlr, go', ilr, jlr, get_offset
      !!        if (jlr==ilr) then
      !!            if (jlr==iilr) exit ilr_loop! locreg found
      !!            ncount=6*(llr(jlr)%wfd%nseg_c+llr(jlr)%wfd%nseg_f)
      !!            get_offset=get_offset+ncount
      !!        end if
      !!    end do
      !!end do ilr_loop

      !call f_release_routine()

    end function get_offset
    
    
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
    

    !> Transposition of the arrays, variable version (non homogeneous)
    subroutine toglobal_and_transpose(iproc,nproc,orbs,Lzd,comms,psi,&
         work,outadd) !optional
      use module_base
      use module_types
      use communications_base, only: comms_cubic
      use locreg_operations, only: Lpsi_to_global2
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
         workarr = f_malloc0_ptr(Gdim,id='workarr')
         !call to_zero(Gdim,workarr)
         do iorb=1,orbs%norbp
            ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
            ldim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
    
            !!call Lpsi_to_global(Lzd%Glr,Gdim,Lzd%Llr(ilr),psi(psishift1),&
            !!     ldim,orbs%norbp,orbs%nspinor,orbs%nspin,totshift,workarr)
            call Lpsi_to_global2(iproc, ldim, gdim, orbs%norbp, &
                 orbs%nspin, lzd%glr, lzd%llr(ilr), psi(psishift1:), workarr(totshift:))
            psishift1 = psishift1 + ldim
            totshift = totshift + (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor
         end do
    
         !reallocate psi to the global dimensions
         call f_free_ptr(psi)
         psi = f_malloc_ptr(Gdim,id='psi')
         call vcopy(Gdim,workarr(1),1,psi(1),1) !psi=work
         call f_free_ptr(workarr)
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




end module communications



subroutine switch_waves_v(nproc,orbs,nvctr,nvctr_par,psi,psiw)
  use module_defs, only: wp
  use module_types, only: orbitals_data
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
  use module_defs, only: wp
  use module_types, only: orbitals_data
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
  integer :: i,iorb,isp,ikpts
  real(wp), dimension(:,:,:,:), allocatable :: tpsit

  tpsit = f_malloc((/ nvctrp, orbs%nspinor, orbs%norb, orbs%nkpts /),id='tpsit')
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

  call f_free(tpsit)
END SUBROUTINE psitransspi



