!> @file
!!  File defining the structures and the routines for the communication between processes
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining structures and routines related to basic communications
module communications_base
  use module_base
  implicit none

  private

  !> Contains the information needed for communicating the wavefunctions
  !! between processors for the transposition
  type, public :: comms_cubic
     integer, dimension(:), pointer :: ncntd,ncntt
     integer, dimension(:), pointer :: ndspld,ndsplt
     integer, dimension(:,:), pointer :: nvctr_par
  end type comms_cubic

  !> Contains the information needed for communicating in the linear version
  type,public:: comms_linear
    integer :: nptsp_c, ndimpsi_c, ndimind_c, ndimind_f, nptsp_f, ndimpsi_f
    integer,dimension(:),pointer :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
    integer,dimension(:),pointer :: isendbuf_c, iextract_c, iexpand_c, irecvbuf_c
    integer,dimension(:),pointer :: norb_per_gridpoint_c, indexrecvorbital_c
    integer,dimension(:),pointer :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
    integer,dimension(:),pointer :: isendbuf_f, iextract_f, iexpand_f, irecvbuf_f
    integer,dimension(:),pointer :: norb_per_gridpoint_f, indexrecvorbital_f
    integer,dimension(:),pointer :: isptsp_c, isptsp_f !< starting index of a given gridpoint (basically summation of norb_per_gridpoint_*)
    real(kind=8),dimension(:),pointer :: psit_c, psit_f
    integer,dimension(:),pointer :: nsendcounts_repartitionrho, nrecvcounts_repartitionrho
    integer,dimension(:),pointer :: nsenddspls_repartitionrho, nrecvdspls_repartitionrho
    integer :: ncomms_repartitionrho, window
    integer,dimension(:,:),pointer :: commarr_repartitionrho
  end type comms_linear

  !> Contains all parameters needed for the point to point communication of the potential
  type, public :: p2pComms
    integer, dimension(:), pointer :: noverlaps
    real(kind=8), dimension(:), pointer :: recvBuf
    integer, dimension(:,:,:), pointer :: comarr
    integer :: nrecvBuf
    integer :: window
    integer, dimension(:,:), pointer :: ise !< Starting / ending index of recvBuf in x,y,z dimension after communication (glocal coordinates)
    integer, dimension(:,:), pointer :: mpi_datatypes
    logical :: communication_complete
  end type p2pComms

  interface check_array_consistency
     module procedure check_array_consistency0
     module procedure check_array_consistency1
     module procedure check_array_consistency2
  end interface
       
  !> Public routines
  public :: comms_linear_null
  public :: allocate_MPI_communication_arrays
  public :: allocate_local_comms_cubic
  public :: allocate_MPI_comms_cubic_repartition
  public :: allocate_MPI_comms_cubic_repartitionp2p
  public :: deallocate_comms_linear
  public :: deallocate_MPI_communication_arrays
  public :: deallocate_MPI_comms_cubic_repartition
  public :: deallocate_MPI_comms_cubic_repartitionp2p
  public :: deallocate_p2pComms
  public :: allocateCommunicationsBuffersPotential
  public :: deallocateCommunicationsBuffersPotential

  public :: check_array_consistency

  contains

    !> Creators and destructors
    pure function comms_linear_null() result(comms)
      implicit none
      type(comms_linear) :: comms
      call nullify_comms_linear(comms)
    end function comms_linear_null


    pure subroutine nullify_comms_linear(comms)
      implicit none
      type(comms_linear),intent(inout) :: comms
      nullify(comms%nsendcounts_c)
      nullify(comms%nsenddspls_c)
      nullify(comms%nrecvcounts_c)
      nullify(comms%nrecvdspls_c)
      nullify(comms%isendbuf_c)
      nullify(comms%iextract_c)
      nullify(comms%iexpand_c)
      nullify(comms%irecvbuf_c)
      nullify(comms%norb_per_gridpoint_c)
      nullify(comms%indexrecvorbital_c)
      nullify(comms%isptsp_c)
      nullify(comms%psit_c)
      nullify(comms%nsendcounts_f)
      nullify(comms%nsenddspls_f)
      nullify(comms%nrecvcounts_f)
      nullify(comms%nrecvdspls_f)
      nullify(comms%isendbuf_f)
      nullify(comms%iextract_f)
      nullify(comms%iexpand_f)
      nullify(comms%irecvbuf_f)
      nullify(comms%norb_per_gridpoint_f)
      nullify(comms%indexrecvorbital_f)
      nullify(comms%isptsp_f)
      nullify(comms%psit_f)
      nullify(comms%nsendcounts_repartitionrho)
      nullify(comms%nrecvcounts_repartitionrho)
      nullify(comms%nsenddspls_repartitionrho)
      nullify(comms%nrecvdspls_repartitionrho)
      nullify(comms%commarr_repartitionrho)
    end subroutine nullify_comms_linear


    subroutine allocate_MPI_communication_arrays(nproc, comms, only_coarse)
      implicit none
      integer,intent(in) :: nproc
      type(comms_linear),intent(inout) :: comms
      logical,intent(in),optional :: only_coarse
      logical :: allocate_fine
      if (present(only_coarse)) then
          allocate_fine=.not.only_coarse
      else
          allocate_fine=.true.
      end if
      comms%nsendcounts_c=f_malloc_ptr(0.to.nproc-1,id='comms%nsendcounts_c')
      comms%nsenddspls_c=f_malloc_ptr(0.to.nproc-1,id='comms%nsenddspls_c')
      comms%nrecvcounts_c=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvcounts_c')
      comms%nrecvdspls_c=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvdspls_c')
      if (allocate_fine) then
          comms%nsendcounts_f=f_malloc_ptr(0.to.nproc-1,id='comms%nsendcounts_f')
          comms%nsenddspls_f=f_malloc_ptr(0.to.nproc-1,id='comms%nsenddspls_f')
          comms%nrecvcounts_f=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvcounts_f')
          comms%nrecvdspls_f=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvdspls_f')
      end if
    end subroutine allocate_MPI_communication_arrays


    subroutine allocate_local_comms_cubic(comms, only_coarse)
      implicit none
      type(comms_linear),intent(inout) :: comms
      logical,intent(in),optional :: only_coarse
      logical :: allocate_fine
      if (present(only_coarse)) then
          allocate_fine=.not.only_coarse
      else
          allocate_fine=.true.
      end if
      comms%irecvbuf_c=f_malloc_ptr(comms%ndimpsi_c,id='comms%irecvbuf_c')
      comms%indexrecvorbital_c=f_malloc_ptr(comms%ndimind_c,id='comms%indexrecvorbital_c')
      comms%iextract_c=f_malloc_ptr(comms%ndimind_c,id='comms%iextract_c')
      comms%iexpand_c=f_malloc_ptr(comms%ndimind_c,id='comms%iexpand_c')
      comms%isendbuf_c=f_malloc_ptr(comms%ndimpsi_c,id='comms%isendbuf_c')
      comms%isptsp_c=f_malloc_ptr(max(comms%nptsp_c,1),id='comms%isptsp_c')
      comms%norb_per_gridpoint_c=f_malloc_ptr(comms%nptsp_c,id='comms%norb_per_gridpoint_c')
      if (allocate_fine) then
          comms%irecvbuf_f=f_malloc_ptr(comms%ndimpsi_f,id='comms%irecvbuf_f')
          comms%indexrecvorbital_f=f_malloc_ptr(comms%ndimind_f,id='comms%indexrecvorbital_f')
          comms%iextract_f=f_malloc_ptr(comms%ndimind_f,id='comms%iextract_f')
          comms%iexpand_f=f_malloc_ptr(comms%ndimind_f,id='comms%iexpand_f')
          comms%isendbuf_f=f_malloc_ptr(comms%ndimpsi_f,id='comms%isendbuf_f')
          comms%isptsp_f=f_malloc_ptr(max(comms%nptsp_f,1),id='comms%isptsp_f')
          comms%norb_per_gridpoint_f=f_malloc_ptr(comms%nptsp_f,id='comms%norb_per_gridpoint_f')
      end if
    end subroutine allocate_local_comms_cubic


    subroutine allocate_MPI_comms_cubic_repartition(nproc, comms)
      implicit none
      integer,intent(in) :: nproc
      type(comms_linear),intent(inout) :: comms
      comms%nsendcounts_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nsendcounts_repartitionrho')
      comms%nrecvcounts_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvcounts_repartitionrho')
      comms%nsenddspls_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nsenddspls_repartitionrho')
      comms%nrecvdspls_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvdspls_repartitionrho')
    end subroutine allocate_MPI_comms_cubic_repartition


    subroutine allocate_MPI_comms_cubic_repartitionp2p(ncommunications, commarr_repartitionrho)
      implicit none
      integer,intent(in) :: ncommunications
      integer,dimension(:,:),pointer,intent(inout) :: commarr_repartitionrho
      commarr_repartitionrho=f_malloc_ptr((/4,ncommunications/),id='commarr_repartitionrho')
    end subroutine allocate_MPI_comms_cubic_repartitionp2p


    subroutine deallocate_comms_linear(comms)
      implicit none
      type(comms_linear),intent(inout) :: comms
      call deallocate_MPI_communication_arrays(comms)
      call deallocate_local_comms_cubic(comms)
      call deallocate_MPI_comms_cubic_repartition(comms)
      if (associated(comms%psit_c)) call f_free_ptr(comms%psit_c)
      if (associated(comms%psit_f)) call f_free_ptr(comms%psit_f)
      call deallocate_MPI_comms_cubic_repartitionp2p(comms%commarr_repartitionrho)
    end subroutine deallocate_comms_linear


    subroutine deallocate_MPI_communication_arrays(comms)
      implicit none
      type(comms_linear),intent(inout) :: comms
      if (associated(comms%nsendcounts_c)) call f_free_ptr(comms%nsendcounts_c)
      if (associated(comms%nsenddspls_c)) call f_free_ptr(comms%nsenddspls_c)
      if (associated(comms%nrecvcounts_c)) call f_free_ptr(comms%nrecvcounts_c)
      if (associated(comms%nrecvdspls_c)) call f_free_ptr(comms%nrecvdspls_c)
      if (associated(comms%nsendcounts_f)) call f_free_ptr(comms%nsendcounts_f)
      if (associated(comms%nsenddspls_f)) call f_free_ptr(comms%nsenddspls_f)
      if (associated(comms%nrecvcounts_f)) call f_free_ptr(comms%nrecvcounts_f)
      if (associated(comms%nrecvdspls_f)) call f_free_ptr(comms%nrecvdspls_f)
    end subroutine deallocate_MPI_communication_arrays

    subroutine deallocate_local_comms_cubic(comms)
      implicit none
      type(comms_linear),intent(inout) :: comms
      if (associated(comms%irecvbuf_c)) call f_free_ptr(comms%irecvbuf_c)
      if (associated(comms%indexrecvorbital_c)) call f_free_ptr(comms%indexrecvorbital_c)
      if (associated(comms%iextract_c)) call f_free_ptr(comms%iextract_c)
      if (associated(comms%iexpand_c)) call f_free_ptr(comms%iexpand_c)
      if (associated(comms%isendbuf_c)) call f_free_ptr(comms%isendbuf_c)
      if (associated(comms%irecvbuf_f)) call f_free_ptr(comms%irecvbuf_f)
      if (associated(comms%indexrecvorbital_f)) call f_free_ptr(comms%indexrecvorbital_f)
      if (associated(comms%iextract_f)) call f_free_ptr(comms%iextract_f)
      if (associated(comms%iexpand_f)) call f_free_ptr(comms%iexpand_f)
      if (associated(comms%isendbuf_f)) call f_free_ptr(comms%isendbuf_f)
      if (associated(comms%isptsp_c)) call f_free_ptr(comms%isptsp_c)
      if (associated(comms%isptsp_f)) call f_free_ptr(comms%isptsp_f)
      if (associated(comms%norb_per_gridpoint_c)) call f_free_ptr(comms%norb_per_gridpoint_c)
      if (associated(comms%norb_per_gridpoint_f)) call f_free_ptr(comms%norb_per_gridpoint_f)
    end subroutine deallocate_local_comms_cubic

    subroutine deallocate_MPI_comms_cubic_repartition(comms)
      implicit none
      type(comms_linear),intent(inout) :: comms
      if (associated(comms%nsendcounts_repartitionrho)) call f_free_ptr(comms%nsendcounts_repartitionrho)
      if (associated(comms%nrecvcounts_repartitionrho)) call f_free_ptr(comms%nrecvcounts_repartitionrho)
      if (associated(comms%nsenddspls_repartitionrho)) call f_free_ptr(comms%nsenddspls_repartitionrho)
      if (associated(comms%nrecvdspls_repartitionrho)) call f_free_ptr(comms%nrecvdspls_repartitionrho)
    end subroutine deallocate_MPI_comms_cubic_repartition


    subroutine deallocate_MPI_comms_cubic_repartitionp2p(commarr_repartitionrho)
      implicit none
      integer, dimension(:,:), pointer,intent(inout) :: commarr_repartitionrho
      call f_free_ptr(commarr_repartitionrho)
    end subroutine deallocate_MPI_comms_cubic_repartitionp2p

    subroutine deallocate_p2pComms(p2pcomm)
      implicit none
      ! Calling arguments
      type(p2pComms),intent(inout):: p2pcomm
      ! Local variables
      call f_free_ptr(p2pcomm%noverlaps)
      call f_free_ptr(p2pcomm%recvBuf)
      call f_free_ptr(p2pcomm%comarr)
      call f_free_ptr(p2pcomm%ise)
      if (.not.p2pcomm%communication_complete) then
          stop 'cannot deallocate mpi data types if communication has not completed'
      end if
      call f_free_ptr(p2pcomm%mpi_datatypes)
    end subroutine deallocate_p2pComms


    subroutine allocateCommunicationsBuffersPotential(comgp)
      implicit none
      ! Calling arguments
      type(p2pComms),intent(inout):: comgp
      comgp%recvBuf = f_malloc_ptr(comgp%nrecvBuf,id='comgp%recvBuf')
    end subroutine allocateCommunicationsBuffersPotential
    
    
    subroutine deallocateCommunicationsBuffersPotential(comgp)
      implicit none
      type(p2pComms),intent(inout):: comgp
      call f_free_ptr(comgp%recvBuf)
    end subroutine deallocateCommunicationsBuffersPotential


    !> Check the consistency of arrays after a gather (example: atomic coordinates)
    subroutine check_array_consistency0(maxdiff, nproc, array, ndims, mpi_comm)
      use dynamic_memory
      implicit none
      integer, intent(in) :: mpi_comm
      integer, intent(in) :: ndims, nproc
      real(gp), intent(in) :: array
      real(gp), intent(out) :: maxdiff

      integer :: ierr, jproc, i
      real(gp), dimension(:,:), allocatable :: rxyz_glob

      include 'check_array-inc.f90'

    END SUBROUTINE check_array_consistency0


    !> Check the consistency of arrays after a gather (example: atomic coordinates)
    subroutine check_array_consistency1(maxdiff, nproc, array, mpi_comm)
      use dynamic_memory
      implicit none
      integer, intent(in) :: mpi_comm
      integer, intent(in) :: nproc
      real(gp), dimension(:), intent(in) :: array
      real(gp), intent(out) :: maxdiff

      integer :: ierr, jproc, i, ndims
      real(gp), dimension(:,:), allocatable :: rxyz_glob

      ndims = size(array)

      include 'check_array-inc.f90'

    END SUBROUTINE check_array_consistency1


    !> Check the consistency of arrays after a gather (example: atomic coordinates)
    subroutine check_array_consistency2(maxdiff, nproc, array, mpi_comm)
      use dynamic_memory
      implicit none
      integer, intent(in) :: mpi_comm
      integer, intent(in) :: nproc
      real(gp), dimension(:,:), intent(in) :: array
      real(gp), intent(out) :: maxdiff

      integer :: ierr, jproc, i, ndims
      real(gp), dimension(:,:), allocatable :: rxyz_glob

      ndims = size(array)

      include 'check_array-inc.f90'

    END SUBROUTINE check_array_consistency2

end module communications_base
