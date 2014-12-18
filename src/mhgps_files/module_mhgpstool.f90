!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module module_mhgpstool
    use module_base
    implicit none

    private

    !datatypes
    public :: mhgpstool_data

    !routines
    public :: read_folders
    public :: read_data
    public :: count_saddle_points
    public :: init_mhgpstool_data
    public :: finalize_mhgpstool_data

    type mhgpstool_data
        integer  :: nid
        real(gp), allocatable :: fp_arr(:,:,:)
        real(gp), allocatable :: fp_arr_sad(:,:,:)
        real(gp), allocatable :: en_arr(:,:)
        real(gp), allocatable :: en_arr_sad(:,:)

        integer, allocatable  :: sadneighb(:,:,:)

        integer, allocatable  :: sadneighb_merged(:,:)
        integer, allocatable  :: minnumber_merged(:)
        integer, allocatable  :: sadnumber_merged(:)

        integer, allocatable  :: exclude(:)
    end type
    contains
!=====================================================================
subroutine init_mhgpstool_data(nat,nfolder,nsad,mdat)
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    integer, intent(in) :: nfolder
    integer, intent(in) :: nsad(:)
    type(mhgpstool_data), intent(inout) :: mdat
    !local
    integer :: nid
    integer :: nsadmax, nminmax, nsadtot, nmintot
    
    mdat%nid = nat !s-overlap
!    nid = 4*nat !sp-overlap

    nsadmax = maxval(nsad)
    nminmax = 2*nsadmax !worst case
    nsadtot = sum(nsad)
    nmintot = 2*nsadtot !worst case
    mdat%fp_arr     = f_malloc((/nid,nminmax,nfolder/),id='fp_arr')
    mdat%fp_arr_sad = f_malloc((/nid,nsadmax,nfolder/),id='fp_arr_sad')
    mdat%en_arr     = f_malloc((/nminmax,nfolder/),id='en_arr')
    mdat%en_arr_sad = f_malloc((/nsadmax,nfolder/),id='en_arr_sad')

    mdat%sadneighb  = f_malloc((/2,nsadmax,nfolder/),id='sadneighb')

    mdat%sadneighb_merged = f_malloc((/2,nsadtot/),id='sadneighb_merged')
    mdat%minnumber_merged = f_malloc((/nmintot/),id='minnumber_merged')
    mdat%sadnumber_merged = f_malloc((/nsadtot/),id='sadnumber_merged')

    mdat%exclude = f_malloc((/nsadtot/),id='exclude')
end subroutine init_mhgpstool_data
!=====================================================================
subroutine finalize_mhgpstool_data(mdat)
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat


    call f_free(mdat%fp_arr)
    call f_free(mdat%fp_arr_sad)
    call f_free(mdat%en_arr)
    call f_free(mdat%en_arr_sad)

    call f_free(mdat%sadneighb)

    call f_free(mdat%sadneighb_merged)
    call f_free(mdat%minnumber_merged)
    call f_free(mdat%sadnumber_merged)

    call f_free(mdat%exclude)
end subroutine finalize_mhgpstool_data
!=====================================================================
subroutine read_folders(nfolder,folders)
    use module_base
    implicit none
    !parameters
    integer, intent(out) :: nfolder
    character(len=500), allocatable :: folders(:)
    !local
    integer :: u, istat
    integer :: ifolder
    character(len=500) :: line
    u=f_get_free_unit()
    open(u,file='mhgpstool.inp')
    nfolder=0
    do
        read(u,'(a)',iostat=istat)line
        if(istat/=0)exit
        nfolder=nfolder+1
    enddo
    folders = f_malloc_str(500,(/1.to.nfolder/),id='folders')
    rewind(u)
    do ifolder=1,nfolder
        read(u,'(a)')folders(ifolder)
    enddo
    close(u)
end subroutine read_folders
!=====================================================================
subroutine count_saddle_points(nfolder,folders,nsad)
    use yaml_output
    use module_io, only: check_struct_file_exists
    implicit none
    !parameters
    integer, intent(in) :: nfolder
    character(len=500), intent(in) :: folders(:)
    integer, intent(out) :: nsad(nfolder)
    !local
    integer :: ifolder
    integer :: isad
    character(len=5) :: isadc
    character(len=600) :: fsaddle, fminR, fminL
    logical :: fsaddleEx, fminRex, fminLex
    fsaddleEx=.false.; fminRex=.false.; fminLex=.false.
    call yaml_comment('Saddle counts ....',hfill='-')
    do ifolder = 1, nfolder
        isad=0
        do
            write(fsaddle,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
                               '/sad',isad+1,'_finalF'
            call check_struct_file_exists(fsaddle,fsaddleEx)
            write(fminL,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
                  '/sad',isad+1,'_minFinalL'
            call check_struct_file_exists(fminL,fminLex)
            write(fminR,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
                  '/sad',isad+1,'_minFinalR'
            call check_struct_file_exists(fminR,fminRex)
            if((.not. fsaddleEx) .or. (.not. fminLex) .or. (.not. fminRex))then
                exit
            endif
            isad=isad+1
        enddo
        nsad(ifolder)=isad
        call yaml_map(trim(adjustl(folders(ifolder))),isad)
    enddo
    call yaml_map('TOTAL',sum(nsad))
end subroutine count_saddle_points
!=====================================================================
subroutine read_data(folders,mdat)
    use module_base
    implicit none
    !parameters
    character(len=500), intent(in) :: folders(:)
    type(mhgpstool_data), intent(inout) :: mdat
    !local
    integer :: nfolder
    nfolder = size(folders,1)

    
end subroutine read_data
!=====================================================================

end module
