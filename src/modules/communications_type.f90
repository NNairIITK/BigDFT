module communications_type
  use module_base
  implicit none

  private

  type,public:: collective_comms
    integer :: nptsp_c, ndimpsi_c, ndimind_c, ndimind_f, nptsp_f, ndimpsi_f
    integer,dimension(:),pointer :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
    integer,dimension(:),pointer :: isendbuf_c, iextract_c, iexpand_c, irecvbuf_c
    integer,dimension(:),pointer :: norb_per_gridpoint_c, indexrecvorbital_c
    integer,dimension(:),pointer :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
    integer,dimension(:),pointer :: isendbuf_f, iextract_f, iexpand_f, irecvbuf_f
    integer,dimension(:),pointer :: norb_per_gridpoint_f, indexrecvorbital_f
    integer,dimension(:),pointer :: isptsp_c, isptsp_f !<starting index of a given gridpoint (basically summation of norb_per_gridpoint_*)
    real(kind=8),dimension(:),pointer :: psit_c, psit_f
    integer,dimension(:),pointer :: nsendcounts_repartitionrho, nrecvcounts_repartitionrho
    integer,dimension(:),pointer :: nsenddspls_repartitionrho, nrecvdspls_repartitionrho
    integer :: ncomms_repartitionrho, window
    integer,dimension(:,:),pointer :: commarr_repartitionrho
  end type collective_comms

end module communications_type


