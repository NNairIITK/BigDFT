!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module module_globaltool
    use module_base
    use module_atoms, only: atomic_structure
    implicit none

    private

    !datatypes
    public :: gt_data

    !routines
    public :: read_globaltool_uinp
    public :: write_globaltool_uinp
!!    public :: read_and_merge_data
!!    public :: count_saddle_points
    public :: init_gt_data
    public :: finalize_gt_data
!!    public :: write_data
    public :: getPairId
    public :: unpair

    type gt_uinp
        real(gp) :: en_delta
        real(gp) :: fp_delta
        integer  :: ndict
        character(len=500), allocatable :: directories(:)
    end type

    type gt_data
        integer :: nid
        integer :: nat
        integer :: ntrans
        integer :: nmin
        integer :: ntransax
        integer :: nminmax
        type(atomic_structure) :: astruct
        type(gt_uinp) :: uinp
        real(gp), allocatable :: rcov(:)
        real(gp), allocatable :: fp_arr(:,:)
        real(gp), allocatable :: en_arr(:)
        integer, allocatable  :: transpairs(:)!Use Cantors pairing function
                                               !to identify pairs
        !counts how many distinct neighbored minimum pairs 
        !a saddle has:
        integer, allocatable  :: nneighbpairs(:)
        !counts how often a minimum pair is found:
        integer, allocatable  :: paircounter(:,:)
        integer, allocatable  :: minnumber(:)
        integer, allocatable  :: sadnumber(:)
        integer, allocatable  :: exclude(:)
        character(len=600), allocatable :: path_sad(:)
        character(len=600), allocatable :: path_min(:)
    end type
    contains
!=====================================================================
function getPairId(IDmin1,IDmin2)
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: IDmin1, IDmin2
    integer :: getPairId
    !local
    integer :: k1,k2

    k1= min(IDmin1,Idmin2)
    k2= max(IDmin1,Idmin2)

    !Cantor's pairing function:
    getPairId = (k1+k2)*(k1+k2+1)/2 + k2
end function getPairId
!=====================================================================
subroutine unpair(pairID,IDmin1,IDmin2)
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: pairID
    integer, intent(out) :: IDmin1, IDmin2
    !local
    integer :: w,t
    integer :: k1,k2

    w = floor(0.5_gp*(sqrt(real(8*pairID+1,gp))-1.0_gp))
    t = (w**2 + w)/2
    k1 = pairID - t
    k2 = w - k1
    IDmin1 = min(k1,k2)
    IDmin2 = max(k1,k2)
end subroutine unpair
!=====================================================================
subroutine count_poslocm(gdat)
    use module_base
    use yaml_output
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    character(len=600) :: filename
    integer :: idict,ifile
    character(len=4) :: cifile
    logical :: exists

    call yaml_comment('Counting poslocms ....',hfill='-')
    gdat%nminmax=0
    do idict=1,gdat%uinp%ndict
        ifile=0
        do
            write(filename,'(a,i4.4,a)')trim(adjustl(&
                 gdat%uinp%directories(idict)))//'/poslocm_',&
                ifile+1,'_'
            call check_struct_file_exists(filename,exists=exists)
            if(.not. exists)exit
            ifile=ifile+1
        enddo
        call yaml_map('Number of poslocms in '//trim(adjustl(&
             gdat%uinp%directories(idict))),ifile)
        gdat%nminmax=gdat%nminmax+ifile
    enddo
    call yaml_map('Total poslocms found',gdat%nminmax)
end subroutine
!=====================================================================
subroutine init_nat_rcov(gdat)
    use module_base
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    implicit none
    !parameter
    type(gt_data), intent(inout) :: gdat
    !local
    call deallocate_atomic_structure(gdat%astruct)
    call set_astruct_from_file(trim(adjustl(gdat%uinp%directories(1)))//&
         '/poslocm_0001_',0,gdat%astruct)
    gdat%nat=gdat%astruct%nat
    gdat%rcov = f_malloc((/gdat%nat/),id='rcov')
    call give_rcov(0,gdat%astruct,gdat%rcov)
end subroutine
!=====================================================================
subroutine init_gt_data(gdat)
    use module_base
    use module_atoms, only: nullify_atomic_structure
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local

    call nullify_atomic_structure(gdat%astruct)
     
    call count_poslocm(gdat)

    call init_nat_rcov(gdat)
    gdat%nid = gdat%nat !s-overlap
!    nid = 4*nat !sp-overlap


    gdat%fp_arr = f_malloc((/gdat%nid,gdat%nminmax/),id='fp_arr')
    gdat%en_arr = f_malloc((/gdat%nminmax/),id='en_arr')
end subroutine init_gt_data
!=====================================================================
subroutine finalize_gt_data(gdat)
    use module_atoms, only: deallocate_atomic_structure
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat

    call deallocate_atomic_structure(gdat%astruct)
    call f_free_str(500,gdat%uinp%directories)
    call f_free(gdat%fp_arr)
    call f_free(gdat%en_arr)
    call f_free(gdat%rcov)
end subroutine finalize_gt_data
!=====================================================================
subroutine read_globaltool_uinp(gdat)
    use module_base
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    integer :: u, istat, idict
    character(len=500) :: line
    real(gp) :: rdmy
    u=f_get_free_unit()
    open(u,file='globaltool.inp')
    gdat%uinp%ndict=0
    read(u,*,iostat=istat)gdat%uinp%en_delta,gdat%uinp%fp_delta
    if(istat/=0)then
        call f_err_throw('Error in first line of globaltool.inp',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
    do
        read(u,'(a)',iostat=istat)line
        if(istat/=0)exit
        gdat%uinp%ndict=gdat%uinp%ndict+1
    enddo
    gdat%uinp%directories = f_malloc_str(500,(/1.to.gdat%uinp%ndict/),&
                            id='gdat%uinp%directories')
    rewind(u)
    read(u,*,iostat=istat)gdat%uinp%en_delta,gdat%uinp%fp_delta
    do idict=1,gdat%uinp%ndict
        read(u,'(a)')gdat%uinp%directories(idict)
    enddo
    close(u)
end subroutine read_globaltool_uinp
!=====================================================================
subroutine write_globaltool_uinp(gdat)
    use module_base
    use yaml_output
    implicit none
    !parameters
    type(gt_data), intent(in) :: gdat
    !local
    integer :: idict

    call yaml_comment('User input ....',hfill='-')
    call yaml_mapping_open('thresholds')
    call yaml_map('en_delta',gdat%uinp%en_delta)
    call yaml_map('fp_delta',gdat%uinp%fp_delta)
    call yaml_mapping_close()
    call yaml_map('Number of directories',gdat%uinp%ndict)
    call yaml_mapping_open('directories')
    do idict=1,gdat%uinp%ndict
    call yaml_scalar(trim(adjustl(gdat%uinp%directories(idict))))
    enddo
    call yaml_mapping_close()
end subroutine write_globaltool_uinp
!=====================================================================
!!!subroutine count_saddle_points(nfolder,folders,nsad)
!!!    use yaml_output
!!!    use module_io, only: check_struct_file_exists
!!!    implicit none
!!!    !parameters
!!!    integer, intent(in) :: nfolder
!!!    character(len=500), intent(in) :: folders(:)
!!!    integer, intent(out) :: nsad(nfolder)
!!!    !local
!!!    integer :: ifolder
!!!    integer :: isad, isadfolder
!!!    character(len=5) :: isadc
!!!    character(len=600) :: fsaddle, fminR, fminL
!!!    logical :: fsaddleEx, fminRex, fminLex
!!!    fsaddleEx=.false.; fminRex=.false.; fminLex=.false.
!!!    call yaml_comment('Saddle counts ....',hfill='-')
!!!    isad=0
!!!    do ifolder = 1, nfolder
!!!        isadfolder=0
!!!        do
!!!            call construct_filenames(folders,ifolder,isadfolder+1,fsaddle,&
!!!                 fminL,fminR)
!!!            call check_struct_file_exists(fsaddle,fsaddleEx)
!!!            call check_struct_file_exists(fminL,fminLex)
!!!            call check_struct_file_exists(fminR,fminRex)
!!!            if((.not. fsaddleEx) .or. (.not. fminLex) .or. &
!!!                                                  (.not. fminRex))then
!!!                exit
!!!            endif
!!!            isad=isad+1
!!!            isadfolder=isadfolder+1
!!!        enddo
!!!        nsad(ifolder)=isadfolder
!!!        call yaml_map(trim(adjustl(folders(ifolder))),isadfolder)
!!!    enddo
!!!    call yaml_map('TOTAL',sum(nsad))
!!!end subroutine count_saddle_points
!!!!=====================================================================
!!!subroutine read_and_merge_data(folders,nsad,mdat)
!!!    use module_base
!!!    use yaml_output
!!!    use module_atoms, only: set_astruct_from_file,&
!!!                            deallocate_atomic_structure
!!!    use module_fingerprints
!!!    implicit none
!!!    !parameters
!!!    character(len=500), intent(in) :: folders(:)
!!!    integer, intent(in) :: nsad(:)
!!!    type(gt_data), intent(inout) :: mdat
!!!    !local
!!!    integer :: nfolder
!!!    integer :: ifolder
!!!    integer :: isad
!!!    character(len=600) :: fsaddle, fminR, fminL
!!!    real(gp) :: rxyz(3,mdat%nat), fp(mdat%nid), epot
!!!    real(gp) :: en_delta, fp_delta
!!!    real(gp) :: en_delta_sad, fp_delta_sad
!!!    logical  :: lnew
!!!    integer  :: kid
!!!    integer  :: k_epot
!!!    integer  :: id_minleft, id_minright, id_saddle
!!!    integer :: isadfolder
!!!    nfolder = size(folders,1)
!!!
!!!    en_delta = mdat%mhgps_uinp%en_delta_min
!!!    fp_delta = mdat%mhgps_uinp%fp_delta_min
!!!    en_delta_sad = mdat%mhgps_uinp%en_delta_sad
!!!    fp_delta_sad = mdat%mhgps_uinp%fp_delta_sad
!!!    mdat%nexclude=0
!!!    call yaml_comment('Thresholds ....',hfill='-')
!!!    call yaml_map('en_delta',en_delta)
!!!    call yaml_map('fp_delta',fp_delta)
!!!    call yaml_map('en_delta_sad',en_delta_sad)
!!!    call yaml_map('fp_delta_sad',fp_delta_sad)
!!!    call yaml_comment('Merging ....',hfill='-')
!!!    do ifolder =1, nfolder
!!!        isadfolder=0
!!!        do isad =1, nsad(ifolder)
!!!            isadfolder=isadfolder+1
!!!            call construct_filenames(folders,ifolder,isadfolder,fsaddle,&
!!!                 fminL,fminR)
!!!            call deallocate_atomic_structure(mdat%astruct)
!!!            !insert left minimum
!!!            call set_astruct_from_file(trim(fminL),0,mdat%astruct,&
!!!                 energy=epot)
!!!            if (mdat%astruct%nat /= mdat%nat) then
!!!                call f_err_throw('Error in read_and_merge_data:'//&
!!!                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
!!!                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
!!!                     err_name='BIGDFT_RUNTIME_ERROR')
!!!            end if
!!!            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
!!!                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
!!!                 fp(1))
!!!write(*,*)trim(adjustl(fminL))
!!!write(*,*)'***'
!!!            call identical('min',mdat,mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
!!!                 mdat%en_arr,mdat%fp_arr,en_delta,fp_delta,lnew,kid,&
!!!                 k_epot)
!!!            if(lnew)then
!!!                call yaml_comment('Minimum '//trim(adjustl(fminL))//&
!!!                                  ' is new.')
!!!                call insert_min(mdat,k_epot,epot,fp,fminL)
!!!                id_minleft=mdat%minnumber(k_epot+1)
!!!            else
!!!                id_minleft=mdat%minnumber(kid)
!!!                call yaml_comment('Minimum '//trim(adjustl(fminL))//&
!!!                                  ' is identical to minimum '//&
!!!                                  trim(yaml_toa(id_minleft)))
!!!            endif 
!!!            call deallocate_atomic_structure(mdat%astruct)
!!!
!!!            !insert right minimum
!!!            call set_astruct_from_file(trim(fminR),0,mdat%astruct,&
!!!                 energy=epot)
!!!            if (mdat%astruct%nat /= mdat%nat) then
!!!                call f_err_throw('Error in read_and_merge_data:'//&
!!!                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
!!!                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
!!!                     err_name='BIGDFT_RUNTIME_ERROR')
!!!            end if
!!!            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
!!!                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
!!!                 fp(1))
!!!
!!!write(*,*)trim(adjustl(fminR))
!!!write(*,*)'***'
!!!            call identical('min',mdat,mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
!!!                 mdat%en_arr,mdat%fp_arr,en_delta,fp_delta,lnew,kid,&
!!!                 k_epot)
!!!            if(lnew)then
!!!                call yaml_comment('Minimum '//trim(adjustl(fminR))//&
!!!                                  ' is new.')
!!!                call insert_min(mdat,k_epot,epot,fp,fminR)
!!!                id_minright=mdat%minnumber(k_epot+1)
!!!            else
!!!                id_minright=mdat%minnumber(kid)
!!!                call yaml_comment('Minimum '//trim(adjustl(fminR))//&
!!!                                  ' is identical to minimum '//&
!!!                                  trim(yaml_toa(id_minright)))
!!!            endif 
!!!
!!!            !insert saddle
!!!            call deallocate_atomic_structure(mdat%astruct)
!!!            call set_astruct_from_file(trim(fsaddle),0,mdat%astruct,&
!!!                 energy=epot)
!!!            if (mdat%astruct%nat /= mdat%nat) then
!!!                call f_err_throw('Error in read_and_merge_data:'//&
!!!                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
!!!                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
!!!                     err_name='BIGDFT_RUNTIME_ERROR')
!!!            end if
!!!            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
!!!                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
!!!                 fp(1))
!!!            call identical('sad',mdat,mdat%nsadtot,mdat%nsad,mdat%nid,epot,fp,&
!!!                 mdat%en_arr_sad,mdat%fp_arr_sad,en_delta_sad,&
!!!                 fp_delta_sad,lnew,kid,k_epot)
!!!            if(lnew)then
!!!                call yaml_comment('Saddle '//trim(adjustl(fsaddle))//&
!!!                     ' is new.')
!!!                call insert_sad(mdat,k_epot,epot,fp,id_minleft,&
!!!                     id_minright,fsaddle)
!!!                id_saddle=mdat%sadnumber(k_epot+1)
!!!            else
!!!                id_saddle=mdat%sadnumber(kid)
!!!                call yaml_comment('Saddle '//trim(adjustl(fsaddle))//&
!!!                     ' is identical to saddle '//trim(yaml_toa(id_saddle)))
!!!                call add_neighbors(mdat,kid,id_minleft,id_minright)
!!!!                if(.not.( ((mdat%sadneighb(1,kid)==id_minleft)&
!!!!                         .and.(mdat%sadneighb(2,kid)==id_minright))&
!!!!                     &.or.((mdat%sadneighb(2,kid)==id_minleft) &
!!!!                         .and.(mdat%sadneighb(1,kid)==id_minright))))then
!!!!                    call yaml_warning('following saddle point has'//&
!!!!                         ' more than two neighboring minima: '//&
!!!!                         trim(yaml_toa(id_saddle)))
!!!!                    mdat%nexclude=mdat%nexclude+1
!!!!                    mdat%exclude(mdat%nexclude) = id_saddle
!!!!                endif
!!!
!!!            endif 
!!!        enddo
!!!    enddo
!!!    
!!!end subroutine read_and_merge_data
!!!!=====================================================================
!!!subroutine write_data(mdat)
!!!    use yaml_output
!!!    use module_base
!!!    implicit none
!!!    !parameters
!!!    type(gt_data), intent(inout) :: mdat
!!!    !local
!!!    integer :: u, u2, u3
!!!    integer :: imin, isad
!!!    integer, allocatable :: mn(:)
!!!    integer :: ipair, it
!!!    logical :: exclude
!!!    character(len=5) :: ci
!!!    integer :: isadc, iminc
!!!
!!!    mn = f_malloc((/mdat%nmin/),id='mn')
!!!
!!!    !write mdat file for minima
!!!    u3=f_get_free_unit()
!!!    open(u3,file='copy_configurations.sh')
!!!    write(u3,*)'#!/bin/bash'
!!!    write(u3,*)'mkdir minima'
!!!    write(u3,*)'mkdir saddlepoints'
!!!    u=f_get_free_unit()
!!!    open(u,file='mindat')
!!!    do imin = 1,mdat%nmin
!!!        mn(mdat%minnumber(imin)) = imin
!!!        write(u,*)mdat%en_arr(imin)
!!!        write(ci,'(i5.5)')imin
!!!        write(u3,*)'cp '//trim(adjustl(mdat%path_min(imin)))//&
!!!                   '.EXT minima/min'//ci//'.EXT'
!!!    enddo 
!!!    close(u)
!!!    
!!!    !write tsdat file for saddle points and connection information
!!!    open(u,file='tsdat')
!!!    u2=f_get_free_unit()
!!!    open(u2,file='tsdat_exclude')
!!!    isadc=0
!!!    do isad=1,mdat%nsad
!!!        exclude=.false.
!!!        ipair=maxloc(mdat%paircounter(1:mdat%nneighbpairs(isad),isad),1)
!!!if(mdat%nneighbpairs(isad)>5)exclude=.true.
!!!!        do it = 1, mdat%nneighbpairs(isad)
!!!!            if(it/=ipair)then
!!!!                 if(mdat%paircounter(it,isad)>20000)then
!!!!            call yaml_comment('Saddle '//trim(adjustl(yaml_toa(mdat%sadnumber(isad))))//&
!!!!                 'converged at least twice to another minimum pair.'//&
!!!!                 ' Too ambigous. Will not consider this saddle point.')
!!!!                    exclude=.true.
!!!!                  endif
!!!!            endif
!!!!        enddo
!!!!!write(*,*)mdat%paircounter(:,isad)
!!!write(*,*)'imaxloc',ipair
!!!        if(exclude)then
!!!            write(u2,'(es24.17,1x,a,2(1x,i0.0))')mdat%en_arr_sad(isad),&
!!!                 '0   0',mn(mdat%sadneighb(1,ipair,isad)),&
!!!                  mn(mdat%sadneighb(2,ipair,isad))
!!!        else
!!!            isadc=isadc+1
!!!            write(u,'(es24.17,1x,a,2(1x,i0.0))')mdat%en_arr_sad(isad),&
!!!                 '0   0',mn(mdat%sadneighb(1,ipair,isad)),&
!!!                  mn(mdat%sadneighb(2,ipair,isad))
!!!            write(ci,'(i5.5)')isadc
!!!            write(u3,*)'cp '//trim(adjustl(mdat%path_sad(isad)))//&
!!!                       '.EXT saddlepoints/sad'//ci//'.EXT'
!!!        endif
!!!do ipair=1,mdat%nneighbpairs(isad)
!!!write(*,*)ipair,mdat%paircounter(ipair,isad)
!!!enddo
!!!    
!!!!        if(.not. any(mdat%exclude .eq. mdat%sadnumber(isad)))then
!!!!            write(u,'(es24.17,1x,a,2(1x,i0.0))')mdat%en_arr_sad(isad),&
!!!!                 '0   0',mn(mdat%sadneighb(1,isad)),mn(mdat%sadneighb(2,isad))
!!!!        else
!!!!            write(u2,'(es24.17,1x,a,2(1x,i0.0))')mdat%en_arr_sad(isad),&
!!!!                 '0   0',mn(mdat%sadneighb(1,isad)),mn(mdat%sadneighb(2,isad))
!!!!        endif
!!!    enddo
!!!    close(u)
!!!    close(u2)
!!!    close(u3)
!!!    call f_free(mn)
!!!
!!!end subroutine write_data
!!!!=====================================================================
!!!subroutine identical(cf,mdat,ndattot,ndat,nid,epot,fp,en_arr,fp_arr,en_delta,&
!!!                    fp_delta,lnew,kid,k_epot)
!!!    use module_base
!!!    use yaml_output
!!!    use module_fingerprints
!!!    implicit none
!!!    !parameters
!!!    type(gt_data), intent(in) :: mdat
!!!    integer, intent(in) :: ndattot
!!!    integer, intent(in) :: ndat
!!!    integer, intent(in) :: nid
!!!    real(gp), intent(in) :: epot
!!!    real(gp), intent(in) :: fp(nid)
!!!    real(gp), intent(in) :: en_arr(ndattot)
!!!    real(gp), intent(in) :: fp_arr(nid,ndattot)
!!!    real(gp), intent(in) :: en_delta
!!!    real(gp), intent(in) :: fp_delta
!!!    logical, intent(out) :: lnew
!!!    integer, intent(out) :: kid
!!!    integer, intent(out) :: k_epot
!!!    character(len=3), intent(in) :: cf
!!!    !local
!!!    integer :: k, klow, khigh, nsm
!!!    real(gp) :: dmin, d 
!!!    !search in energy array
!!!    call hunt_mt(en_arr,max(1,min(ndat,ndattot)),epot,k_epot)
!!!    lnew=.true.
!!!    
!!!    ! find lowest configuration that might be identical
!!!    klow=k_epot
!!!    do k=k_epot,1,-1
!!!        if (epot-en_arr(k).lt.0.d0) stop 'zeroA'
!!!        if (epot-en_arr(k).gt.en_delta) exit
!!!        klow=k
!!!    enddo
!!!    
!!!    ! find highest  configuration that might be identical
!!!    khigh=k_epot+1
!!!    do k=k_epot+1,ndat
!!!        if (en_arr(k)-epot.lt.0.d0) stop 'zeroB'
!!!        if (en_arr(k)-epot.gt.en_delta) exit
!!!        khigh=k
!!!    enddo
!!!    
!!!    nsm=0
!!!    dmin=huge(1.e0_gp)
!!!    do k=max(1,klow),min(ndat,khigh)
!!!        call fpdistance(nid,fp,fp_arr(1,k),d)
!!!write(*,*)'fpdist '//cf,abs(en_arr(k)-epot),d
!!!!if(cf=='min')then
!!!!if(d<3.d-3)then
!!!!if(abs(en_arr(k)-epot)>1.d-4)then
!!!!write(*,*)trim(adjustl(mdat%path_min(k)))
!!!!    stop
!!!!endif
!!!!endif
!!!!endif
!!!        if (d.lt.fp_delta) then
!!!            lnew=.false.
!!!            nsm=nsm+1
!!!            if (d.lt.dmin) then 
!!!                dmin=d
!!!                kid=k
!!!            endif
!!!        endif
!!!    enddo
!!!write(*,*)'dmin',dmin
!!!    if (nsm.gt.1) then
!!!        call yaml_warning('more than one identical configuration'//&
!!!             ' found')
!!!    endif
!!!end subroutine identical
!!!!=====================================================================
!!!subroutine add_neighbors(mdat,kid,neighb1,neighb2)
!!!    use module_base
!!!    implicit none
!!!    !parameters
!!!    type(gt_data), intent(inout) :: mdat
!!!    integer, intent(in) :: kid
!!!    integer, intent(in) :: neighb1, neighb2
!!!    !local
!!!    integer :: ipair
!!!    logical :: found
!!!
!!!    !first check, if neihgbor pair is already
!!!    !in list
!!!    found=.false.
!!!write(*,*)
!!!write(*,*)neighb1,neighb2
!!!write(*,*)'---'
!!!    neighbloop: do ipair=1,mdat%nneighbpairs(kid)
!!!        if( ((mdat%sadneighb(1,ipair,kid)==neighb1)&
!!!               .and.(mdat%sadneighb(2,ipair,kid)==neighb2))&
!!!           &.or.((mdat%sadneighb(2,ipair,kid)==neighb1) &
!!!               .and.(mdat%sadneighb(1,ipair,kid)==neighb2)) )then
!!!            mdat%paircounter(ipair,kid) = mdat%paircounter(ipair,kid)+1
!!!            found=.true.
!!!            exit neighbloop
!!!        endif
!!!    enddo neighbloop
!!!
!!!    if(.not. found) then !pair is new, add it to list
!!!        mdat%nneighbpairs(kid) = mdat%nneighbpairs(kid) + 1
!!!        mdat%paircounter(mdat%nneighbpairs(kid),kid)  = 1
!!!        mdat%sadneighb(1,mdat%nneighbpairs(kid),kid) = neighb1
!!!        mdat%sadneighb(2,mdat%nneighbpairs(kid),kid) = neighb2
!!!    endif
!!!do ipair=1,mdat%nneighbpairs(kid)
!!!write(*,*)mdat%sadneighb(1,ipair,kid),mdat%sadneighb(2,ipair,kid),mdat%paircounter(ipair,kid)
!!!enddo
!!!end subroutine add_neighbors
!!!!=====================================================================
!!!subroutine insert_sad(mdat,k_epot,epot,fp,neighb1,neighb2,path)
!!!    !insert at k_epot+1
!!!    use module_base
!!!    implicit none
!!!    !parameters
!!!    type(gt_data), intent(inout) :: mdat
!!!    integer, intent(in) :: k_epot
!!!    real(gp), intent(in) :: epot
!!!    real(gp), intent(in) :: fp(mdat%nid)
!!!    integer, intent(in) :: neighb1, neighb2
!!!    character(len=600)   :: path
!!!    !local
!!!    integer :: i,k
!!!    if(mdat%nsad+1>mdat%nsadtot)stop 'nsad+1>=nsadtot, out of bounds'
!!!
!!!    mdat%nsad=mdat%nsad+1
!!!    do k=mdat%nsad-1,k_epot+1,-1
!!!        mdat%en_arr_sad(k+1)=mdat%en_arr_sad(k)
!!!        mdat%sadnumber(k+1)=mdat%sadnumber(k)
!!!        mdat%path_sad(k+1)=mdat%path_sad(k)
!!!        mdat%nneighbpairs(k+1) = mdat%nneighbpairs(k)
!!!        mdat%paircounter(:,k+1) = mdat%paircounter(:,k)
!!!        mdat%sadneighb(1,:,k+1)=mdat%sadneighb(1,:,k)
!!!        mdat%sadneighb(2,:,k+1)=mdat%sadneighb(2,:,k)
!!!        do i=1,mdat%nid
!!!            mdat%fp_arr_sad(i,k+1)=mdat%fp_arr_sad(i,k)
!!!         enddo
!!!    enddo
!!!    mdat%en_arr_sad(k_epot+1)=epot
!!!    mdat%sadnumber(k_epot+1)=mdat%nsad
!!!    mdat%path_sad(k_epot+1)=path
!!!    mdat%nneighbpairs(k_epot+1) = 1
!!!    mdat%paircounter(1,k_epot+1) = 1
!!!    mdat%sadneighb(1,1,k_epot+1)=neighb1
!!!    mdat%sadneighb(2,1,k_epot+1)=neighb2
!!!    do i=1,mdat%nid
!!!        mdat%fp_arr_sad(i,k+1)=fp(i)
!!!    enddo
!!!end subroutine insert_sad
!!!!=====================================================================
!!!subroutine insert_min(mdat,k_epot,epot,fp,path)
!!!    !insert at k_epot+1
!!!    use module_base
!!!    implicit none
!!!    !parameters
!!!    type(gt_data), intent(inout) :: mdat
!!!    integer, intent(in) :: k_epot
!!!    real(gp), intent(in) :: epot
!!!    real(gp), intent(in) :: fp(mdat%nid)
!!!    character(len=600)   :: path
!!!    !local
!!!    integer :: i,k
!!!write(*,*)'insert at',k_epot+1
!!!    if(mdat%nmin+1>mdat%nmintot)stop 'nmin+1>=nmintot, out of bounds'
!!!
!!!    mdat%nmin=mdat%nmin+1
!!!    do k=mdat%nmin-1,k_epot+1,-1
!!!        mdat%en_arr(k+1)=mdat%en_arr(k)
!!!        mdat%minnumber(k+1)=mdat%minnumber(k)
!!!        mdat%path_min(k+1)=mdat%path_min(k)
!!!        do i=1,mdat%nid
!!!            mdat%fp_arr(i,k+1)=mdat%fp_arr(i,k)
!!!         enddo
!!!    enddo
!!!    mdat%en_arr(k_epot+1)=epot
!!!    mdat%minnumber(k_epot+1)=mdat%nmin
!!!    mdat%path_min(k_epot+1)=path
!!!    do i=1,mdat%nid
!!!        mdat%fp_arr(i,k_epot+1)=fp(i)
!!!    enddo
!!!end subroutine insert_min
!!!!=====================================================================
!!!subroutine construct_filenames(folders,ifolder,isad,fsaddle,fminL,fminR)
!!!    implicit none
!!!    !parameters
!!!    character(len=500), intent(in) :: folders(:)
!!!    integer, intent(in)  ::  ifolder, isad
!!!    character(len=600), intent(out) :: fsaddle, fminR, fminL
!!!    !local
!!!    
!!!    write(fsaddle,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
!!!                       '/sad',isad,'_finalF'
!!!    write(fminL,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
!!!          '/sad',isad,'_minFinalL'
!!!    write(fminR,'(a,i5.5,a)')trim(adjustl(folders(ifolder)))//&
!!!          '/sad',isad,'_minFinalR'
!!!end subroutine
!!!!=====================================================================
!!!!> C x is in interval [xx(jlo),xx(jlow+1)
!!!![ ; xx(0)=-Infinity ; xx(n+1) = Infinity
!!!subroutine hunt_mt(xx,n,x,jlo)
!!!  use module_base
!!!  implicit none
!!!  !Arguments
!!!  integer :: jlo,n
!!!  real(gp) :: x,xx(n)
!!!  !Local variables
!!!  integer :: inc,jhi,jm
!!!  logical :: ascnd
!!!  if (n.le.0) stop 'hunt_mt'
!!!  if (n == 1) then
!!!     if (x.ge.xx(1)) then
!!!        jlo=1
!!!     else
!!!        jlo=0
!!!     endif
!!!     return
!!!  endif
!!!  ascnd=xx(n).ge.xx(1)
!!!  if(jlo.le.0.or.jlo.gt.n)then
!!!     jlo=0
!!!     jhi=n+1
!!!     goto 3
!!!  endif
!!!  inc=1
!!!  if(x.ge.xx(jlo).eqv.ascnd)then
!!!1    continue
!!!     jhi=jlo+inc
!!!     if(jhi.gt.n)then
!!!        jhi=n+1
!!!     else if(x.ge.xx(jhi).eqv.ascnd)then
!!!        jlo=jhi
!!!        inc=inc+inc
!!!        goto 1
!!!     endif
!!!  else
!!!     jhi=jlo
!!!2    continue
!!!     jlo=jhi-inc
!!!     if(jlo.lt.1)then
!!!        jlo=0
!!!     else if(x.lt.xx(jlo).eqv.ascnd)then
!!!        jhi=jlo
!!!        inc=inc+inc
!!!        goto 2
!!!     endif
!!!  endif
!!!3 continue
!!!  if(jhi-jlo == 1)then
!!!     if(x == xx(n))jlo=n
!!!     if(x == xx(1))jlo=1
!!!     return
!!!  endif
!!!  jm=(jhi+jlo)/2
!!!  if(x.ge.xx(jm).eqv.ascnd)then
!!!     jlo=jm
!!!  else
!!!     jhi=jm
!!!  endif
!!!  goto 3
!!!END SUBROUTINE hunt_mt
subroutine check_struct_file_exists(filename,exists)
    use module_base
    implicit none
    !parameter
    character(len=*), intent(in) :: filename
    logical, optional, intent(out) :: exists
    !local
    logical :: xyzexists=.false.,asciiexists=.false.
    integer :: indx,inda

    if(present(exists))then
        exists=.true.
    endif

    indx=index(filename,'.xyz')
    inda=index(filename,'.ascii')
    if(indx==0)then
        inquire(file=trim(adjustl(filename))//'.xyz',exist=xyzexists)
    else
        inquire(file=trim(adjustl(filename)),exist=xyzexists)
    endif
    if(inda==0)then
        inquire(file=trim(adjustl(filename))//'.ascii',&
                exist=asciiexists)
    else
        inquire(file=trim(adjustl(filename)),exist=asciiexists)
    endif
    if(.not. (xyzexists .or. asciiexists))then
        if(present(exists))then
            exists=.false.
        else
            call f_err_throw('File '//trim(adjustl(filename))//&
                             ' does not exist.')
        endif
    endif
end subroutine
subroutine give_rcov(iproc,astruct,rcov)
  use module_base
!  use module_types
  use yaml_output
  use module_atoms
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(atomic_structure), intent(in) :: astruct
  real(kind=gp), dimension(astruct%nat), intent(out) :: rcov
  !Local variables
!  type(atoms_iterator) :: it
  integer :: iat

  do iat=1,astruct%nat
    call covalent_radius(astruct%atomnames(astruct%iatype(iat)),rcov(iat))
    if (iproc == 0) call yaml_map('(MH) RCOV '//trim(astruct%atomnames(astruct%iatype(iat))),rcov(iat))
  end do
  !python metod
!  it=atoms_iter(astruct)
!  do while(atoms_iter_next(it))
!    call covalent_radius(it%name,rcov(it%iat))
!  end do

end subroutine give_rcov




end module
