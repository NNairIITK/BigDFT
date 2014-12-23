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
    use module_atoms, only: atomic_structure
    use module_userinput, only: userinput
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
    public :: write_data

    type mhgpstool_data
        integer :: nid
        integer :: nat
        integer :: nsad
        integer :: nmin
        integer :: nsadtot
        integer :: nmintot
        type(atomic_structure) :: astruct
        type(userinput) :: mhgps_uinp
        real(gp), allocatable :: rcov(:)
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
    use module_atoms, only: nullify_atomic_structure
    implicit none
    !parameters
    integer, intent(in) :: nat
    integer, intent(in) :: nfolder
    integer, intent(in) :: nsad(:)
    type(mhgpstool_data), intent(inout) :: mdat
    !local
   
    call nullify_atomic_structure(mdat%astruct)
 
    mdat%nat = nat
    mdat%nid = nat !s-overlap
!    nid = 4*nat !sp-overlap

    mdat%rcov = f_malloc((/nat/),id='rcov')

    mdat%nsad = 0
    mdat%nmin = 0

    mdat%nsadtot = sum(nsad)
    mdat%nmintot = 2*mdat%nsadtot !worst case
    mdat%fp_arr     = f_malloc((/mdat%nid,mdat%nmintot/),id='fp_arr')
    mdat%fp_arr_sad = f_malloc((/mdat%nid,mdat%nsadtot/),id='fp_arr_sad')
    mdat%en_arr     = f_malloc((/mdat%nmintot/),id='en_arr')
    mdat%en_arr_sad = f_malloc((/mdat%nsadtot/),id='en_arr_sad')

    mdat%sadneighb  = f_malloc((/2,mdat%nsadtot/),id='sadneighb')

    mdat%minnumber = f_malloc((/mdat%nmintot/),id='minnumber')
    mdat%sadnumber = f_malloc((/mdat%nsadtot/),id='sadnumber')

    mdat%exclude = f_malloc((/mdat%nsadtot/),id='exclude')
    
    mdat%path_sad = f_malloc_str(600,(/1.to.mdat%nsadtot/),id='path_sad')
    mdat%path_min = f_malloc_str(600,(/1.to.mdat%nmintot/),id='path_min')
end subroutine init_mhgpstool_data
!=====================================================================
subroutine finalize_mhgpstool_data(mdat)
    use module_atoms, only: deallocate_atomic_structure
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat


    call deallocate_atomic_structure(mdat%astruct)
    call f_free(mdat%rcov)
    
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
    integer :: isad, isadfolder
    character(len=5) :: isadc
    character(len=600) :: fsaddle, fminR, fminL
    logical :: fsaddleEx, fminRex, fminLex
    fsaddleEx=.false.; fminRex=.false.; fminLex=.false.
    call yaml_comment('Saddle counts ....',hfill='-')
    isad=0
    do ifolder = 1, nfolder
        isadfolder=0
        do
            call construct_filenames(folders,ifolder,isad+1,fsaddle,&
                 fminL,fminR)
            call check_struct_file_exists(fsaddle,fsaddleEx)
            call check_struct_file_exists(fminL,fminLex)
            call check_struct_file_exists(fminR,fminRex)
            if((.not. fsaddleEx) .or. (.not. fminLex) .or. &
                                                  (.not. fminRex))then
                exit
            endif
            isad=isad+1
            isadfolder=isadfolder+1
        enddo
        nsad(ifolder)=isadfolder
        call yaml_map(trim(adjustl(folders(ifolder))),isadfolder)
    enddo
    call yaml_map('TOTAL',sum(nsad))
end subroutine count_saddle_points
!=====================================================================
subroutine read_and_merge_data(folders,nsad,mdat)
    use module_base
    use yaml_output
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    use module_fingerprints
    implicit none
    !parameters
    character(len=500), intent(in) :: folders(:)
    integer, intent(in) :: nsad(:)
    type(mhgpstool_data), intent(inout) :: mdat
    !local
    integer :: nfolder
    integer :: ifolder
    integer :: isad
    character(len=600) :: fsaddle, fminR, fminL
    real(gp) :: rxyz(3,mdat%nat), fp(mdat%nid), epot
    real(gp) :: en_delta, fp_delta
    real(gp) :: en_delta_sad, fp_delta_sad
    logical  :: lnew
    integer  :: kid
    integer  :: k_epot
    integer  :: id_minleft, id_minright, id_saddle
    integer :: isadcount
integer :: i
    nfolder = size(folders,1)

    en_delta = mdat%mhgps_uinp%en_delta_min
    fp_delta = mdat%mhgps_uinp%fp_delta_min
    en_delta_sad = mdat%mhgps_uinp%en_delta_sad
    fp_delta_sad = mdat%mhgps_uinp%fp_delta_sad
    call yaml_comment('Thresholds ....',hfill='-')
    call yaml_map('en_delta',en_delta)
    call yaml_map('fp_delta',fp_delta)
    call yaml_map('en_delta_sad',en_delta_sad)
    call yaml_map('fp_delta_sad',fp_delta_sad)
    call yaml_comment('Merging ....',hfill='-')
    isadcount=0
    do ifolder =1, nfolder
        do isad =1, nsad(ifolder)
            isadcount=isadcount+1
            call construct_filenames(folders,ifolder,isadcount,fsaddle,&
                 fminL,fminR)
            call deallocate_atomic_structure(mdat%astruct)
            !insert left minimum
            call set_astruct_from_file(trim(fminL),0,mdat%astruct,&
                 energy=epot)
            if (mdat%astruct%nat /= mdat%nat) then
                call f_err_throw('Error in read_and_merge_data:'//&
                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
                     err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
                 fp(1))
            call identical(mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
                 mdat%en_arr,mdat%fp_arr,en_delta,fp_delta,lnew,kid,&
                 k_epot)
            if(lnew)then
                call yaml_comment('Minimum '//trim(adjustl(fminL))//&
                                  ' is new.')
                call insert_min(mdat,k_epot,epot,fp,fminL)
                id_minleft=mdat%minnumber(k_epot+1)
            else
                id_minleft=mdat%minnumber(kid)
                call yaml_comment('Minimum '//trim(adjustl(fminL))//&
                                  ' is identical to minimum '//&
                                  trim(yaml_toa(id_minleft)))
            endif 
            call deallocate_atomic_structure(mdat%astruct)

            !insert right minimum
            call set_astruct_from_file(trim(fminR),0,mdat%astruct,&
                 energy=epot)
            if (mdat%astruct%nat /= mdat%nat) then
                call f_err_throw('Error in read_and_merge_data:'//&
                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
                     err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
                 fp(1))

            call identical(mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
                 mdat%en_arr,mdat%fp_arr,en_delta,fp_delta,lnew,kid,&
                 k_epot)
            if(lnew)then
                call yaml_comment('Minimum '//trim(adjustl(fminR))//&
                                  ' is new.')
                call insert_min(mdat,k_epot,epot,fp,fminR)
                id_minright=mdat%minnumber(k_epot+1)
            else
                id_minright=mdat%minnumber(kid)
                call yaml_comment('Minimum '//trim(adjustl(fminR))//&
                                  ' is identical to minimum '//&
                                  trim(yaml_toa(id_minright)))
            endif 

            !insert saddle
            call deallocate_atomic_structure(mdat%astruct)
            call set_astruct_from_file(trim(fsaddle),0,mdat%astruct,&
                 energy=epot)
            if (mdat%astruct%nat /= mdat%nat) then
                call f_err_throw('Error in read_and_merge_data:'//&
                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
                     err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
                 fp(1))
            call identical(mdat%nsadtot,mdat%nsad,mdat%nid,epot,fp,&
                 mdat%en_arr_sad,mdat%fp_arr_sad,en_delta_sad,&
                 fp_delta_sad,lnew,kid,k_epot)
            if(lnew)then
                call yaml_comment('Saddle '//trim(adjustl(fsaddle))//&
                     ' is new.')
                call insert_sad(mdat,k_epot,epot,fp,id_minleft,&
                     id_minright,fsaddle)
                id_saddle=mdat%sadnumber(k_epot+1)
            else
                id_saddle=mdat%sadnumber(kid)
                call yaml_comment('Saddle '//trim(adjustl(fsaddle))//&
                     ' is identical to saddle '//trim(yaml_toa(id_saddle)))
                if(.not.( ((mdat%sadneighb(1,kid)==id_minleft)&
                         .and.(mdat%sadneighb(2,kid)==id_minright))&
                     &.or.((mdat%sadneighb(2,kid)==id_minleft) &
                         .and.(mdat%sadneighb(1,kid)==id_minright))))then
                    call yaml_warning('following saddle point has'//&
                         ' more than two neighboring minima: '//&
                         trim(yaml_toa(id_saddle)))
!                    nexclude=nexclude+1
!                    exclude(nexclude) = sadnumber_merged(kid_sad)
                endif

            endif 
        enddo
    enddo
    
end subroutine read_and_merge_data
!=====================================================================
subroutine write_data(mdat)
    use module_base
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat
    !local
    integer :: u
    integer :: imin, isad
    integer, allocatable :: mn(:)

    mn = f_malloc((/mdat%nmin/),id='mn')

    !write mdat file for minima
    u=f_get_free_unit()
    open(u,file='mindat')
    do imin = 1,mdat%nmin
        mn(mdat%minnumber(imin)) = imin
        write(u,*)mdat%en_arr(imin)
    enddo 
    close(u)
    
    !write tsdat file for saddle points and connection information
    open(u,file='tsdat')
    do isad=1,mdat%nsad
        write(u,'(es24.17,1x,a,2(1x,i0.0))')mdat%en_arr_sad(isad),'0   0',mn(mdat%sadneighb(1,isad)),mn(mdat%sadneighb(2,isad))
    enddo
    close(u)
    call f_free(mn)

end subroutine write_data
!=====================================================================
subroutine identical(ndattot,ndat,nid,epot,fp,en_arr,fp_arr,en_delta,&
                    fp_delta,lnew,kid,k_epot)
    use module_base
    use yaml_output
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
write(*,*)'dmin',dmin
    if (nsm.gt.1) then
        call yaml_warning('more than one identical configuration'//&
             ' found')
    endif
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
write(*,*)'insert at',k_epot+1
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
    mdat%sadneighb(2,k_epot+1)=neighb2
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
write(*,*)'insert at',k_epot+1
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
        mdat%fp_arr(i,k_epot+1)=fp(i)
    enddo
end subroutine insert_min
!=====================================================================
subroutine construct_filenames(folders,ifolder,isad,fsaddle,fminL,fminR)
    implicit none
    !parameters
    character(len=500), intent(in) :: folders(:)
    integer, intent(in)  ::  ifolder, isad
    character(len=600), intent(out) :: fsaddle, fminR, fminL
    !local
    
    write(fsaddle,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
                       '/sad',isad,'_finalF'
    write(fminL,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
          '/sad',isad,'_minFinalL'
    write(fminR,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
          '/sad',isad,'_minFinalR'
end subroutine
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
