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
    public :: read_and_merge_data
    public :: count_saddle_points
    public :: init_mhgpstool_data
    public :: finalize_mhgpstool_data

    type mhgpstool_data
        integer :: nid
        integer :: nat
        integer :: nsad
        integer :: nmin
        integer :: nsadtot
        integer :: nmintot
        real(gp), allocatable :: fp_arr(:,:)
        real(gp), allocatable :: en_arr(:)
        real(gp), allocatable :: fp_arr_sad(:,:)
        real(gp), allocatable :: en_arr_sad(:)
        integer, allocatable  :: sadneighb(:,:)
        integer, allocatable  :: minnumber(:)
        integer, allocatable  :: sadnumber(:)
        integer, allocatable  :: exclude(:)
        character(len=600), allocatable :: path_sad(:)
        character(len=600), allocatable :: path_min(:)
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
    integer :: nsadtot, nmintot
    
    mdat%nat = nat
    mdat%nid = nat !s-overlap
!    nid = 4*nat !sp-overlap

    mdat%nsad = -1
    mdat%nmin = -1

    mdat%nsadtot = sum(nsad)
    mdat%nmintot = 2*nsadtot !worst case
    mdat%fp_arr     = f_malloc((/nid,nmintot/),id='fp_arr')
    mdat%fp_arr_sad = f_malloc((/nid,nsadtot/),id='fp_arr_sad')
    mdat%en_arr     = f_malloc((/nmintot/),id='en_arr')
    mdat%en_arr_sad = f_malloc((/nsadtot/),id='en_arr_sad')

    mdat%sadneighb  = f_malloc((/2,nsadtot/),id='sadneighb')

    mdat%minnumber = f_malloc((/nmintot/),id='minnumber')
    mdat%sadnumber = f_malloc((/nsadtot/),id='sadnumber')

    mdat%exclude = f_malloc((/nsadtot/),id='exclude')
    
    mdat%path_sad = f_malloc_str(600,(/1.to.nsadtot/),id='path_sad')
    mdat%path_min = f_malloc_str(600,(/1.to.nmintot/),id='path_min')
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

    call f_free(mdat%minnumber)
    call f_free(mdat%sadnumber)

    call f_free(mdat%exclude)
    call f_free_str(600,mdat%path_sad)
    call f_free_str(600,mdat%path_min)
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
subroutine read_and_merge_data(folders,mdat)
    use module_base
    implicit none
    !parameters
    character(len=500), intent(in) :: folders(:)
    type(mhgpstool_data), intent(inout) :: mdat
    !local
    integer :: nfolder
    integer :: ifolder
    nfolder = size(folders,1)

    do ifolder =1, nfolder
        
    enddo
    
end subroutine read_and_merge_data
!=====================================================================
subroutine write_data()
    implicit none
    !parameters
    !local
end subroutine write_data
!=====================================================================
subroutine identical(ndattot,ndat,nid,epot,fp,en_arr,fp_arr,en_delta,&
                    fp_delta,lnew,kid,k_epot)
    use module_base
    use module_fingerprints
    implicit none
    !parameters
    integer, intent(in) :: ndattot
    integer, intent(in) :: ndat
    integer, intent(in) :: nid
    real(gp), intent(in) :: epot
    real(gp), intent(in) :: fp(nid)
    real(gp), intent(in) :: en_arr(ndattot)
    real(gp), intent(in) :: fp_arr(nid,ndattot)
    real(gp), intent(in) :: en_delta
    real(gp), intent(in) :: fp_delta
    logical, intent(out) :: lnew
    integer, intent(out) :: kid
    integer, intent(out) :: k_epot
    !local
    integer :: k, klow, khigh, nsm
    real(gp) :: dmin, d 
    !search in energy array
    call hunt_mt(en_arr,max(1,min(ndat,ndattot)),epot,k_epot)
    lnew=.true.
    
    ! find lowest configuration that might be identical
    klow=k_epot
    do k=k_epot,1,-1
        if (epot-en_arr(k).lt.0.d0) stop 'zeroA'
        if (epot-en_arr(k).gt.en_delta) exit
        klow=k
    enddo
    
    ! find highest  configuration that might be identical
    khigh=k_epot+1
    do k=k_epot+1,ndat
        if (en_arr(k)-epot.lt.0.d0) stop 'zeroB'
        if (en_arr(k)-epot.gt.en_delta) exit
        khigh=k
    enddo
    
    nsm=0
    dmin=huge(1.e0_gp)
    do k=max(1,klow),min(ndat,khigh)
        call fpdistance(nid,fp,fp_arr(1,k),d)
        if (d.lt.fp_delta) then
            lnew=.false.
            nsm=nsm+1
            if (d.lt.dmin) then 
                dmin=d
                kid=k
            endif
        endif
    enddo
    if (nsm.gt.1) write(*,*) '(MH) WARNING: more than one'//&
                             ' identical configuration found'
end subroutine identical
!=====================================================================
subroutine insert_sad(mdat,k_epot,epot,fp,neighb1,neighb2,path)
    !insert at k_epot+1
    use module_base
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat
    integer, intent(in) :: k_epot
    real(gp), intent(in) :: epot
    real(gp), intent(in) :: fp(mdat%nid)
    integer, intent(in) :: neighb1, neighb2
    character(len=600)   :: path
    !local
    integer :: i,k
    if(mdat%nsad+1>mdat%nsadtot)stop 'nsad+1>=nsadtot, out of bounds'

    mdat%nsad=mdat%nsad+1
    do k=mdat%nsad-1,k_epot+1,-1
        mdat%en_arr_sad(k+1)=mdat%en_arr_sad(k)
        mdat%sadnumber(k+1)=mdat%sadnumber(k)
        mdat%path_sad(k+1)=mdat%path_sad(k)
        mdat%sadneighb(1,k+1)=mdat%sadneighb(1,k)
        mdat%sadneighb(2,k+1)=mdat%sadneighb(2,k)
        do i=1,mdat%nid
            mdat%fp_arr_sad(i,k+1)=mdat%fp_arr_sad(i,k)
         enddo
    enddo
    mdat%en_arr_sad(k_epot+1)=epot
    mdat%sadnumber(k_epot+1)=mdat%nsad
    mdat%path_sad(k_epot+1)=path
    mdat%sadneighb(1,k_epot+1)=neighb1
    mdat%sadneighb(1,k_epot+1)=neighb2
    do i=1,mdat%nid
        mdat%fp_arr_sad(i,k+1)=fp(i)
    enddo
end subroutine insert_sad
!=====================================================================
subroutine insert_min(mdat,k_epot,epot,fp,path)
    !insert at k_epot+1
    use module_base
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat
    integer, intent(in) :: k_epot
    real(gp), intent(in) :: epot
    real(gp), intent(in) :: fp(mdat%nid)
    character(len=600)   :: path
    !local
    integer :: i,k
    if(mdat%nmin+1>mdat%nmintot)stop 'nmin+1>=nmintot, out of bounds'

    mdat%nmin=mdat%nmin+1
    do k=mdat%nmin-1,k_epot+1,-1
        mdat%en_arr(k+1)=mdat%en_arr(k)
        mdat%minnumber(k+1)=mdat%minnumber(k)
        mdat%path_min(k+1)=mdat%path_min(k)
        do i=1,mdat%nid
            mdat%fp_arr(i,k+1)=mdat%fp_arr(i,k)
         enddo
    enddo
    mdat%en_arr(k_epot+1)=epot
    mdat%minnumber(k_epot+1)=mdat%nmin
    mdat%path_min(k_epot+1)=path
    do i=1,mdat%nid
        mdat%fp_arr(i,k+1)=fp(i)
    enddo
end subroutine insert_min
!=====================================================================
!> C x is in interval [xx(jlo),xx(jlow+1)
![ ; xx(0)=-Infinity ; xx(n+1) = Infinity
subroutine hunt_mt(xx,n,x,jlo)
  use module_base
  implicit none
  !Arguments
  integer :: jlo,n
  real(gp) :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  if (n.le.0) stop 'hunt_mt'
  if (n == 1) then
     if (x.ge.xx(1)) then
        jlo=1
     else
        jlo=0
     endif
     return
  endif
  ascnd=xx(n).ge.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    continue
     jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    continue
     jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 continue
  if(jhi-jlo == 1)then
     if(x == xx(n))jlo=n
     if(x == xx(1))jlo=1
     return
  endif
  jm=(jhi+jlo)/2
  if(x.ge.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt_mt



end module
