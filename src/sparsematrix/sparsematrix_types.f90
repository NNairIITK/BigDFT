module sparsematrix_types
  use f_precisions, only: f_double
  use wrapper_MPI, only: mpi_environment
  implicit none

  ! Basic precision
  integer,parameter,public :: mp=f_double  !< matrix-type precision

  !> Contains the matrices
  type,public :: matrices
      real(kind=mp),dimension(:),pointer :: matrix_compr,matrix_comprp
      real(kind=mp),dimension(:,:,:),pointer :: matrix,matrixp
      real(kind=mp) :: power !< power of the matrix; eg.: 0 => original matrix, -1 => inverse, 0.5 => ^1/2, etc.
  end type matrices

  !> Contains the parameters needed for the sparse matrix matrix multiplication
  type,public :: sparse_matrix_matrix_multiplication
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      integer :: nout, nseq, nseg
      integer :: nfvctrp !< modified number of matrix columns per MPI task for an optimized load balancing during matmul
      integer :: isfvctr !< modified starting column of the matrix for an optimized load balancing during matmul
      integer :: nvctrp_mm !< modified number of compressed matrix elements per MPI task
      integer :: isvctr_mm !< modified starting entry of the compressed matrix elements
      integer :: isseg !< segment containing the first entry (i.e. isvctr+1)
      integer :: ieseg !< segment containing the last entry (i.e. isvctr+nvctrp)
      integer,dimension(:),pointer :: isvctr_mm_par, nvctr_mm_par !<array that contains the values of nvctrp_mm and isvctr_mm of all MPI tasks
      integer,dimension(:),pointer :: ivectorindex, ivectorindex_new, nsegline, istsegline, indices_extract_sequential
      integer,dimension(:,:),pointer :: onedimindices, onedimindices_new, line_and_column_mm, line_and_column
      !!integer,dimension(:,:,:),pointer :: keyg
      integer,dimension(2) :: istartendseg_mm !< starting and ending segments of the matrix subpart which is actually used for the multiplication
                                              !! WARNING: the essential bounds are given by istartend_mm, the segments are used to speed up the code
      integer,dimension(2) :: istartend_mm !< starting and ending indices of the matrix subpart which is actually used for the multiplication
      integer,dimension(2) :: istartend_mm_dj !< starting and ending indices of the submatrices (partitioned in disjoint chunks)
      !!integer :: ncl_smmm !< number of elements for the compress local after a sparse matrix matrix multiplication
      integer :: nccomm_smmm !<number of communications required for the compress distributed after a sparse matrix matrix multiplication
      integer,dimension(:,:),pointer :: luccomm_smmm !<lookup array for the communications required for the compress distributed after a sparse matrix matrix multiplication
      integer :: nvctrp !< number of compressed matrix elements per MPI task
      integer :: isvctr !< starting entry of the compressed matrix elements
      integer,dimension(:),pointer :: isvctr_par, nvctr_par !<array that contains the values of nvctrp and isvctr of all MPI tasks
      integer :: nconsecutive_max !< max number of blocks (i.e. consecutive entries) for the sparse matmul
      integer,dimension(:,:,:),pointer :: consecutive_lookup !< lookup arrays for these blocks
  end type sparse_matrix_matrix_multiplication

  type,public :: sparse_matrix
      integer :: nvctr, nseg, nvctrp, isvctr, parallel_compression, nfvctr, nfvctrp, isfvctr, nspin
      integer :: offset_matrixindex_in_compressed_fortransposed
      integer,dimension(:),pointer :: keyv, nsegline, istsegline, isvctr_par, nvctr_par, isfvctr_par, nfvctr_par
      integer,dimension(:,:,:),pointer :: keyg
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_arr!, orb_from_index
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_fortransposed
      logical :: store_index
      type(sparse_matrix_matrix_multiplication) :: smmm
      integer :: ntaskgroup !< total number of MPI taskgroups to which this task belongs
      integer :: ntaskgroupp !< number of MPI taskgroups to which this task belongs
      !> (1:2,1,:) gives the start and end of the taskgroups (in terms of total matrix indices), 
      !! (1:2,2,:) gives the start and end of the disjoint submatrix handled by the taskgroup
      !! The last rank is of dimension ntaskgroup
      integer,dimension(:,:,:),pointer :: taskgroup_startend
      integer,dimension(:),pointer :: taskgroupid !< dimension ntaskgroupp, gives the ID of the taskgroups to which a task belongs
      integer,dimension(:,:),pointer :: inwhichtaskgroup !< dimension (2,0:nproc-1), tells in which taskgroup a given task is
      type(mpi_environment),dimension(:),pointer :: mpi_groups
      integer,dimension(2) :: istartendseg_t !< starting and ending indices of the matrix subpart which is actually used i
                                             !! for the transposed operation (overlap calculation / orthocontraint)
                                             !! WARNING: the essential bounds are given by istartend_t, the segments are you used to speed up the code
      integer,dimension(2) :: istartend_t !< starting and ending indices of the matrix subpart which is actually used i
                                          !! for the transposed operation (overlap calculation / orthocontraint
      integer :: nvctrp_tg !< size of the taskgroup matrix (per processor)
      integer :: isvctrp_tg !< offset of the taskgroup matrix, given with respect to the non-parallelized matrix
      integer,dimension(2) :: iseseg_tg !< first and last segment of the taskgroup sparse matrix
      integer,dimension(2) :: istartend_local !< first and last element of the sparse matrix which is actually used by a given MPI task
      integer,dimension(2) :: istartendseg_local !< first and last segment of the sparse matrix which is actually used by a given MPI task
      integer,dimension(:,:),pointer :: tgranks !< global task IDs of the tasks in each taskgroup
      integer,dimension(:),pointer :: nranks !< number of task on each taskgroup
      integer :: nccomm !<number of communications required for the compress distributed in the dense parallel format
      integer,dimension(:,:),pointer :: luccomm !<lookup array for the communications required for the compress distributed in the dense parallel format
      logical :: smatmul_initialized !< indicated whether the sparse matmul type has been initialized
  end type sparse_matrix

  type,public :: sparse_matrix_metadata
      character(len=1) :: geocode !< boundary conditions F(ree), W(ire), S(urface), P(eriodic)
      real(kind=mp),dimension(3) :: cell_dim !< dimensions of the simulation cell
      integer :: nfvctr !< size of the matrix
      integer :: nat !< number of atoms
      integer :: ntypes !< number of atoms types
      character(len=20) :: units !< units of the atomic positions 
      integer,dimension(:),pointer :: nzatom !< atomic core charge
      integer,dimension(:),pointer :: nelpsp !< number of electrons
      character(len=20),dimension(:),pointer :: atomnames !< name of the atoms
      integer,dimension(:),pointer :: iatype !< indicates the atoms type
      real(kind=mp),dimension(:,:),pointer :: rxyz !< atomic positions
      integer,dimension(:),pointer :: on_which_atom !< indicates which element of the matrix belong to which atom
  end type sparse_matrix_metadata

end module sparsematrix_types
