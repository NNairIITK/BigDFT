!! @file
!! @author Bastian Schaefer
!!    Copyright (C) 2014-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module fo minimum hopping guided path search method
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
    public :: identMHminMHGPSmin
    
    type sadneighb
        integer :: npairx=-1
        integer, allocatable :: neighb(:,:)
        integer, allocatable :: paircounter(:) 
        character(len=600), allocatable :: neighbPath(:,:)
    end type

    type mhgpstool_data
        integer :: nid
        integer :: nat
        integer :: nsad
        integer :: nmin
        integer :: nsadtot
        integer :: nmintot
        integer :: nexclude
        type(atomic_structure) :: astruct
        type(userinput) :: mhgps_uinp
        real(gp), allocatable :: rcov(:)
        real(gp), allocatable :: fp_arr(:,:)
        real(gp), allocatable :: en_arr(:)
        real(gp), allocatable :: fp_arr_sad(:,:)
        real(gp), allocatable :: en_arr_sad(:)
!        integer, allocatable  :: sadneighb(:,:,:)
        type(sadneighb), allocatable :: snghb(:)
        !counts how many distinct neighbored minimum pairs 
        !a saddle has:
        integer, allocatable  :: nneighbpairs(:)
        !counts how often a minimum pair is found:
!        integer, allocatable  :: paircounter(:,:)
        integer, allocatable  :: minnumber(:)
        integer, allocatable  :: sadnumber(:)
        integer, allocatable  :: exclude(:)
        character(len=600), allocatable :: path_sad(:)
        character(len=600), allocatable :: path_min(:)
    end type

    contains
!=====================================================================
subroutine add_sadneighb(snghb,ileft,iright,ipair,fminL,fminR)
    use module_base
    implicit none
    !parameters
    type(sadneighb), intent(inout) :: snghb
    integer, intent(in) :: ileft, iright, ipair
    character(len=*), intent(in) :: fminL, fminR
    !internal
    integer, parameter :: npairxdef=5
    integer :: i
    integer, allocatable :: neighbtmp(:,:)
    integer, allocatable :: counttmp(:)
    character(len=600), allocatable :: neighbPathtmp(:,:)

    if( snghb%npairx < 0)then
!!write(*,*)'hier a'
        allocate(snghb%neighb(2,npairxdef))
        allocate(snghb%neighbPath(2,npairxdef))
        allocate(snghb%paircounter(npairxdef))
!        snghb%neighb = f_malloc((/2,npairxdef/),id='snghb%neighb') 
!        snghb%neighbPath = f_malloc_str(600,(/1.to.2,1.to.npairxdef/),id='snghb%neighbPath')
!        snghb%paircounter = f_malloc((/npairxdef/),id='snghb%paircounter') 
        snghb%paircounter = 0
        snghb%npairx=npairxdef
    else if(ipair>snghb%npairx)then

!!write(*,*)'hier b'
        neighbtmp = f_malloc((/2,snghb%npairx/),id='neighbtmp') 
!        neighbPathtmp = f_malloc_str(600,(/1.to.2,1.to.snghb%npairx/),id='neighbPathtmp')
        allocate(neighbPathtmp(2,snghb%npairx))
        counttmp = f_malloc((/snghb%npairx/),id='counttmp') 
    
        neighbtmp = snghb%neighb
        neighbPathtmp = snghb%neighbPath
        counttmp = snghb%paircounter
       
        deallocate(snghb%neighb)
        deallocate(snghb%neighbPath)
        deallocate(snghb%paircounter)
        !call f_free(snghb%neighb)
        !call f_free(snghb%neighbPath)
        !call f_free(snghb%paircounter)
        snghb%npairx = snghb%npairx + npairxdef
        allocate(snghb%neighb(2,snghb%npairx))
        allocate(snghb%neighbPath(2,snghb%npairx))
        allocate(snghb%paircounter(snghb%npairx))
        !!snghb%neighb = f_malloc((/2,snghb%npairx/),id='snghb%neighb')
        !!snghb%neighbPath = f_malloc_str(600,(/1.to.2,1.to.snghb%npairx/),id='snghb%neighbPath')
        !!snghb%paircounter = f_malloc((/snghb%npairx/),id='snghb%neighb')

        snghb%neighb(:,1:(snghb%npairx-npairxdef)) = neighbtmp(:,1:(snghb%npairx-npairxdef))
        snghb%neighbPath(:,1:(snghb%npairx-npairxdef)) = neighbPathtmp(:,1:(snghb%npairx-npairxdef))
        snghb%paircounter(1:(snghb%npairx-npairxdef)) = counttmp(1:(snghb%npairx-npairxdef))
        snghb%paircounter((snghb%npairx-npairxdef+1):snghb%npairx)=0
        

    
        !deallocate(neighbtmp)
        deallocate(neighbPathtmp)
        !deallocate(counttmp)
        call f_free(neighbtmp)
        !call f_free(neighbPathtmp)
        call f_free(counttmp)
    endif

    snghb%neighb(1,ipair)=ileft
    snghb%neighb(2,ipair)=iright
    snghb%neighbPath(1,ipair)=fminL
    snghb%neighbPath(2,ipair)=fminR
    snghb%paircounter(ipair)=snghb%paircounter(ipair)+1
    
end subroutine
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
    integer :: istat
   
    call nullify_atomic_structure(mdat%astruct)
 
    mdat%nat = nat
    mdat%nid = nat !s-overlap
!    nid = 4*nat !sp-overlap

    mdat%rcov = f_malloc((/nat/),id='rcov')

    mdat%nsad = 0
    mdat%nmin = 0

    mdat%nexclude = 0

    mdat%nsadtot = sum(nsad)
    mdat%nmintot = 2*mdat%nsadtot !worst case
    mdat%fp_arr     = f_malloc((/mdat%nid,mdat%nmintot/),id='fp_arr')
    mdat%fp_arr_sad = f_malloc((/mdat%nid,mdat%nsadtot/),id='fp_arr_sad')
    mdat%en_arr     = f_malloc((/mdat%nmintot/),id='en_arr')
    mdat%en_arr = huge(1.0_gp)
    mdat%en_arr_sad = f_malloc((/mdat%nsadtot/),id='en_arr_sad')
    mdat%en_arr_sad = huge(1.0_gp)


!    mdat%snghb = f_malloc((/mdat%nsadtot/),id='snghb')
    allocate(mdat%snghb(mdat%nsadtot),stat=istat)
    if(istat/=0)stop 'could not allocate mdat%snghb'

!    mdat%sadneighb  = f_malloc((/2,mdat%nsadtot,mdat%nsadtot/),id='sadneighb')
    mdat%nneighbpairs = f_malloc((/mdat%nsadtot/),id='nneighbpairs')
    mdat%nneighbpairs = 0
!    mdat%paircounter   = f_malloc((/mdat%nsadtot,mdat%nsadtot/),id='paircounter')
!    mdat%paircounter  = 0

    mdat%minnumber = f_malloc((/mdat%nmintot/),id='minnumber')
    mdat%sadnumber = f_malloc((/mdat%nsadtot/),id='sadnumber')

    mdat%exclude = f_malloc((/mdat%nsadtot/),id='exclude')
    mdat%exclude = 0
    
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

!    call f_free(mdat%sadneighb)
    call f_free(mdat%nneighbpairs)
!    call f_free(mdat%paircounter)
!!    call f_free(mdat%snghb)
    deallocate(mdat%snghb)

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
    character(len=600) :: fsaddle, fminR, fminL
    logical :: fsaddleEx, fminRex, fminLex

    fsaddleEx=.false.; fminRex=.false.; fminLex=.false.
    call yaml_comment('Saddle counts ....',hfill='-')
    isad=0
    do ifolder = 1, nfolder
        isadfolder=0
        do
            call construct_filenames(folders,ifolder,isadfolder+1,fsaddle,&
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

subroutine identMHminMHGPSmin(MHminPath,mdat,mn)
    use module_base
    use yaml_output
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    use module_fingerprints
!identifies minima from MH database with minima from mhgps database
    implicit none
    !parameters
    character(len=*), intent(in) :: MHminPath
    type(mhgpstool_data), intent(inout) :: mdat
    integer, intent(in) :: mn(mdat%nmin)
    !local
    integer :: imin
    character(len=600) :: filename
    logical :: exists
    real(gp) :: rxyz(3,mdat%nat), fp(mdat%nid), epot
    real(gp) :: en_delta, fp_delta
    logical  :: lnew
    integer  :: kid
    integer  :: k_epot
    en_delta = mdat%mhgps_uinp%en_delta_min
    fp_delta = mdat%mhgps_uinp%fp_delta_min
    imin=0
    call yaml_comment('Mapping of MH minima to MHGPS minima ....',hfill='-')
write(*,*)'ID MH  ---> ID MHGPS '
    do
        imin=imin+1
        write(filename,'(a,i6.6,a)')trim(adjustl(MHminPath))//'/min',imin,'.xyz'
        inquire(file=filename,exist=exists)
        if(.not.exists)exit
        call deallocate_atomic_structure(mdat%astruct)
        !insert left minimum
        call set_astruct_from_file(trim(filename),0,mdat%astruct,&
             energy=epot)
        if (mdat%astruct%nat /= mdat%nat) then
            call f_err_throw('Error in identMHminMHGPSmin:'//&
                 ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                 //' /= '//trim(yaml_toa(mdat%nat))//')',&
                 err_name='BIGDFT_RUNTIME_ERROR')
        end if
        call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
             mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
             fp(1))
        call identical('min',mdat,mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
             mdat%en_arr,mdat%fp_arr,en_delta,fp_delta,lnew,kid,&
             k_epot)
        if(lnew)then
write(*,'(i6.6,1x,a,1x,i6.5)')imin,'-->',-12345
        else
write(*,'(i6.6,1x,a,1x,i6.6)')imin,'-->',mn(mdat%minnumber(kid))
        endif 
    
    enddo
    call deallocate_atomic_structure(mdat%astruct)
end subroutine
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
    real(gp), dimension(mdat%nid) :: fp
    real(gp) :: epot,en_delta, fp_delta
    real(gp) :: en_delta_sad, fp_delta_sad
    logical  :: lnew
    integer  :: kid
    integer  :: k_epot
    integer  :: id_minleft, id_minright, id_saddle
    integer :: isadfolder
    nfolder = size(folders,1)

    en_delta = mdat%mhgps_uinp%en_delta_min
    fp_delta = mdat%mhgps_uinp%fp_delta_min
    en_delta_sad = mdat%mhgps_uinp%en_delta_sad
    fp_delta_sad = mdat%mhgps_uinp%fp_delta_sad
    mdat%nexclude=0
    call yaml_comment('Thresholds ....',hfill='-')
    call yaml_map('en_delta',en_delta)
    call yaml_map('fp_delta',fp_delta)
    call yaml_map('en_delta_sad',en_delta_sad)
    call yaml_map('fp_delta_sad',fp_delta_sad)
    call yaml_comment('Merging ....',hfill='-')
    do ifolder =1, nfolder
        isadfolder=0
        do isad =1, nsad(ifolder)
            isadfolder=isadfolder+1
            call construct_filenames(folders,ifolder,isadfolder,fsaddle,&
                 fminL,fminR)
            call deallocate_atomic_structure(mdat%astruct)
            !insert left minimum
            call set_astruct_from_file(trim(fminL),0,mdat%astruct,&
                 energy=epot)
            if (mdat%astruct%nat /= mdat%nat) then
write(*,*)trim(fminL)
                call f_err_throw('Error in read_and_merge_data:'//&
                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
                     err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
                 fp(1))
write(*,*)trim(adjustl(fminL))
write(*,*)'***'
            call identical('min',mdat,mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
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
write(*,*)trim(fminR)
                call f_err_throw('Error in read_and_merge_data:'//&
                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
                     err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
                 fp(1))

write(*,*)trim(adjustl(fminR))
write(*,*)'***'
            call identical('min',mdat,mdat%nmintot,mdat%nmin,mdat%nid,epot,fp,&
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
write(*,*)trim(fsaddle)
                call f_err_throw('Error in read_and_merge_data:'//&
                     ' wrong size ('//trim(yaml_toa(mdat%astruct%nat))&
                     //' /= '//trim(yaml_toa(mdat%nat))//')',&
                     err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call fingerprint(mdat%nat,mdat%nid,mdat%astruct%cell_dim,&
                 mdat%astruct%geocode,mdat%rcov,mdat%astruct%rxyz,&
                 fp(1))
            call identical('sad',mdat,mdat%nsadtot,mdat%nsad,mdat%nid,epot,fp,&
                 mdat%en_arr_sad,mdat%fp_arr_sad,en_delta_sad,&
                 fp_delta_sad,lnew,kid,k_epot)
            if(lnew)then
                call yaml_comment('Saddle '//trim(adjustl(fsaddle))//&
                     ' is new.')
                call insert_sad(mdat,k_epot,epot,fp,id_minleft,&
                     id_minright,fsaddle,fminL,fminR)
                id_saddle=mdat%sadnumber(k_epot+1)
            else
                id_saddle=mdat%sadnumber(kid)
                call yaml_comment('Saddle '//trim(adjustl(fsaddle))//&
                     ' is identical to saddle '//trim(yaml_toa(id_saddle)))
                call add_neighbors(mdat,kid,id_minleft,id_minright,fminL,fminR)
!                if(.not.( ((mdat%sadneighb(1,kid)==id_minleft)&
!                         .and.(mdat%sadneighb(2,kid)==id_minright))&
!                     &.or.((mdat%sadneighb(2,kid)==id_minleft) &
!                         .and.(mdat%sadneighb(1,kid)==id_minright))))then
!                    call yaml_warning('following saddle point has'//&
!                         ' more than two neighboring minima: '//&
!                         trim(yaml_toa(id_saddle)))
!                    mdat%nexclude=mdat%nexclude+1
!                    mdat%exclude(mdat%nexclude) = id_saddle
!                endif

            endif 
        enddo
    enddo
    
end subroutine read_and_merge_data


subroutine write_data(mdat)
    use yaml_output
    use module_base
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat
    !local
    integer :: u, u2, u3, u4, u5
    integer :: imin, isad
    integer, allocatable :: mn(:)
    integer :: ipair
    logical :: exclude
    character(len=9) :: ci
    integer :: isadc
    integer :: imin_well_aligned

    mn = f_malloc((/mdat%nmin/),id='mn')

    !write mdat file for minima
    u3=f_get_free_unit()
    open(u3,file='copy_configurations.sh')
    u4=f_get_free_unit()
    open(u4,file='copy_well_aligned_configurations.sh')
    write(u3,*)'#!/bin/bash'
    write(u3,*)'mkdir minima'
    write(u3,*)'mkdir saddlepoints'
    write(u4,*)'#!/bin/bash'
    write(u4,*)'mkdir minima_well_aligned'
    write(u4,*)'mkdir saddlepoints'
    u=f_get_free_unit()
    open(u,file='mindat')
    do imin = 1,mdat%nmin
        mn(mdat%minnumber(imin)) = imin
        write(u,*)mdat%en_arr(imin)
        write(ci,'(i9.9)')imin
        write(u3,'(a)')'cp '//trim(adjustl(mdat%path_min(imin)))//&
                   '.EXT minima/min'//ci//'.EXT'
    enddo 
    close(u)
    
    !write tsdat file for saddle points and connection information
    u=f_get_free_unit()
    open(u,file='tsdat')
    u2=f_get_free_unit()
    open(u2,file='tsdat_exclude')
    u5=f_get_free_unit()
    open(u5,file='tsfiles')
    isadc=0
    imin_well_aligned=-1
    do isad=1,mdat%nsad
        imin_well_aligned=imin_well_aligned+2
        exclude=.false.
        ipair=maxloc(mdat%snghb(isad)%paircounter(1:mdat%nneighbpairs(isad)),1)
!!        ipair=maxloc(mdat%paircounter(1:mdat%nneighbpairs(isad),isad),1)
if(mdat%nneighbpairs(isad)>5)exclude=.true.
!        do it = 1, mdat%nneighbpairs(isad)
!            if(it/=ipair)then
!                 if(mdat%paircounter(it,isad)>20000)then
!            call yaml_comment('Saddle '//trim(adjustl(yaml_toa(mdat%sadnumber(isad))))//&
!                 'converged at least twice to another minimum pair.'//&
!                 ' Too ambigous. Will not consider this saddle point.')
!                    exclude=.true.
!                  endif
!            endif
!        enddo
!!write(*,*)mdat%paircounter(:,isad)
write(*,'(a,3(1x,i9.9))')'imaxloc',ipair,min(mn(mdat%snghb(isad)%neighb(1,ipair)),mn(mdat%snghb(isad)%neighb(2,ipair))),&                           
                 max(mn(mdat%snghb(isad)%neighb(1,ipair)),mn(mdat%snghb(isad)%neighb(2,ipair)))
        if(exclude)then
            write(u2,'(es24.17,1x,a,2(1x,i9.9),2x,2(1x,i9.9))')mdat%en_arr_sad(isad),&
                 '0   0',min(mn(mdat%snghb(isad)%neighb(1,ipair)),mn(mdat%snghb(isad)%neighb(2,ipair))),&
                  max(mn(mdat%snghb(isad)%neighb(1,ipair)),mn(mdat%snghb(isad)%neighb(2,ipair))),&
                  imin_well_aligned,imin_well_aligned+1
            write(ci,'(i9.9)')imin_well_aligned
            write(u4,'(a)')'cp '//trim(adjustl(mdat%snghb(isad)%neighbPath(1,ipair)))//&                                                             
                       '.EXT minima_well_aligned/min'//ci//'.EXT'
            write(ci,'(i9.9)')imin_well_aligned+1
            write(u4,'(a)')'cp '//trim(adjustl(mdat%snghb(isad)%neighbPath(2,ipair)))//&                                                             
                       '.EXT minima_well_aligned/min'//ci//'.EXT'
        else
            isadc=isadc+1
            write(u,'(es24.17,1x,a,2(1x,i9.9),2x,2(1x,i9.9))')mdat%en_arr_sad(isad),&
                 '0   0',min(mn(mdat%snghb(isad)%neighb(1,ipair)),mn(mdat%snghb(isad)%neighb(2,ipair))),&
                 max(mn(mdat%snghb(isad)%neighb(1,ipair)),mn(mdat%snghb(isad)%neighb(2,ipair))),&
                 imin_well_aligned,imin_well_aligned+1
            write(ci,'(i9.9)')isadc
            write(u3,'(a)')'cp '//trim(adjustl(mdat%path_sad(isad)))//&
                       '.EXT saddlepoints/sad'//ci//'.EXT'
            write(u4,'(a)')'cp '//trim(adjustl(mdat%path_sad(isad)))//&
                       '.EXT saddlepoints/sad'//ci//'.EXT'
            write(u5,'(a,x,a)')'"'//trim(adjustl(mdat%path_sad(isad)))//'"',ci
            write(ci,'(i9.9)')imin_well_aligned
            write(u4,'(a)')'cp '//trim(adjustl(mdat%snghb(isad)%neighbPath(1,ipair)))//&                                                             
                       '.EXT minima_well_aligned/min'//ci//'.EXT'
            write(ci,'(i9.9)')imin_well_aligned+1
            write(u4,'(a)')'cp '//trim(adjustl(mdat%snghb(isad)%neighbPath(2,ipair)))//&                                                             
                       '.EXT minima_well_aligned/min'//ci//'.EXT'
        endif
do ipair=1,mdat%nneighbpairs(isad)
write(*,*)ipair,mdat%snghb(isad)%paircounter(ipair)
enddo
    
!        if(.not. any(mdat%exclude .eq. mdat%sadnumber(isad)))then
!            write(u,'(es24.17,1x,a,2(1x,i9.9))')mdat%en_arr_sad(isad),&
!                 '0   0',mn(mdat%sadneighb(1,isad)),mn(mdat%sadneighb(2,isad))
!        else
!            write(u2,'(es24.17,1x,a,2(1x,i9.9))')mdat%en_arr_sad(isad),&
!                 '0   0',mn(mdat%sadneighb(1,isad)),mn(mdat%sadneighb(2,isad))
!        endif
    enddo
    close(u)
    close(u2)
    close(u3)
    close(u4)
    close(u5)

    !uncomment following call of identMHminMHGPSmin if
    !identification of MH minima id with minima ID in MHPGS databse
    !is desired
!    call identMHminMHGPSmin('MH_database',mdat,mn)

    call f_free(mn)
end subroutine write_data
!=====================================================================
subroutine identical(cf,mdat,ndattot,ndat,nid,epot,fp,en_arr,fp_arr,en_delta,&
                    fp_delta,lnew,kid,k_epot)
    use module_base
    use yaml_output
    use module_fingerprints
    implicit none
    !parameters
    type(mhgpstool_data), intent(in) :: mdat
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
    character(len=3), intent(in) :: cf
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
        if (abs(epot-en_arr(k)).le.en_delta) then
            call fpdistance(nid,fp,fp_arr(1,k),d)
            write(*,*)'fpdist '//cf,abs(en_arr(k)-epot),d
            if (d.lt.fp_delta) then
                lnew=.false.
                nsm=nsm+1
                if (d.lt.dmin) then 
                    dmin=d
                    kid=k
                endif
            endif
            dmin=min(dmin,d)
        endif
    enddo
    if (nsm.gt.1) then
        call yaml_warning('more than one identical configuration'//&
             ' found')
    endif
end subroutine identical
!=====================================================================
subroutine add_neighbors(mdat,kid,neighb1,neighb2,fminL,fminR)
    use module_base
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat
    integer, intent(in) :: kid
    integer, intent(in) :: neighb1, neighb2
    character(len=*), intent(in) :: fminL, fminR
    !local
    integer :: ipair
    logical :: found

    !first check, if neihgbor pair is already
    !in list
    found=.false.
write(*,*)
write(*,*)neighb1,neighb2
write(*,*)'---'
    neighbloop: do ipair=1,mdat%nneighbpairs(kid)
        if( ((mdat%snghb(kid)%neighb(1,ipair)==neighb1)&
               .and.(mdat%snghb(kid)%neighb(2,ipair)==neighb2))&
           &.or.((mdat%snghb(kid)%neighb(2,ipair)==neighb1) &
               .and.(mdat%snghb(kid)%neighb(1,ipair)==neighb2)) )then
            mdat%snghb(kid)%paircounter(ipair) = mdat%snghb(kid)%paircounter(ipair)+1
            found=.true.
            exit neighbloop
        endif
!!        if( ((mdat%sadneighb(1,ipair,kid)==neighb1)&
!!               .and.(mdat%sadneighb(2,ipair,kid)==neighb2))&
!!           &.or.((mdat%sadneighb(2,ipair,kid)==neighb1) &
!!               .and.(mdat%sadneighb(1,ipair,kid)==neighb2)) )then
!!            mdat%paircounter(ipair,kid) = mdat%paircounter(ipair,kid)+1
!!            found=.true.
!!            exit neighbloop
!!        endif
    enddo neighbloop

    if(.not. found) then !pair is new, add it to list
        mdat%nneighbpairs(kid) = mdat%nneighbpairs(kid) + 1
        call add_sadneighb(mdat%snghb(kid),neighb1,neighb2,mdat%nneighbpairs(kid),fminL,fminR)
!!        mdat%paircounter(mdat%nneighbpairs(kid),kid)  = 1
!!        mdat%sadneighb(1,mdat%nneighbpairs(kid),kid) = neighb1
!!        mdat%sadneighb(2,mdat%nneighbpairs(kid),kid) = neighb2
    endif
do ipair=1,mdat%nneighbpairs(kid)
write(*,*)mdat%snghb(kid)%neighb(1,ipair),mdat%snghb(kid)%neighb(2,ipair),mdat%snghb(kid)%paircounter(ipair)
enddo
end subroutine add_neighbors
!=====================================================================
subroutine insert_sad(mdat,k_epot,epot,fp,neighb1,neighb2,fsad,fminL,fminR)
    !insert at k_epot+1
    use module_base
    implicit none
    !parameters
    type(mhgpstool_data), intent(inout) :: mdat
    integer, intent(in) :: k_epot
    real(gp), intent(in) :: epot
    real(gp), intent(in) :: fp(mdat%nid)
    integer, intent(in) :: neighb1, neighb2
    character(len=600)   :: fsad,fminL,fminR
    !local
    integer :: i,k
    if(mdat%nsad+1>mdat%nsadtot)stop 'nsad+1>=nsadtot, out of bounds'

    mdat%nsad=mdat%nsad+1
    do k=mdat%nsad-1,k_epot+1,-1
        mdat%en_arr_sad(k+1)=mdat%en_arr_sad(k)
        mdat%sadnumber(k+1)=mdat%sadnumber(k)
        mdat%path_sad(k+1)=mdat%path_sad(k)
        mdat%nneighbpairs(k+1) = mdat%nneighbpairs(k)
        mdat%snghb(k+1) = mdat%snghb(k)
!!        mdat%snghb(k+1)%paircounter(:) = mdat%snghb(k)%paircounter(:)
!!        mdat%snghb(k+1)%neighb(1,:)=mdat%snghb(k)%neighb(1,:)
!!        mdat%snghb(k+1)%neighb(2,:)=mdat%snghb(k)%neighb(2,:)
        do i=1,mdat%nid
            mdat%fp_arr_sad(i,k+1)=mdat%fp_arr_sad(i,k)
         enddo
    enddo
    mdat%en_arr_sad(k_epot+1)=epot
    mdat%sadnumber(k_epot+1)=mdat%nsad
    mdat%path_sad(k_epot+1)=fsad
    mdat%nneighbpairs(k_epot+1) = 1

    call add_sadneighb(mdat%snghb(k_epot+1),neighb1,neighb2,mdat%nneighbpairs(k_epot+1),fminL,fminR)
    mdat%snghb(k_epot+1)%paircounter(1) = 1
!    mdat%snghb(k_epot+1)%paircounter(1) = 1
!    mdat%snghb(k_epot+1)%sadneighb(1,1)=neighb1
!    mdat%snghb(k_epot+1)%sadneighb(2,1)=neighb2
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
                       '/sad',isad,'_finalM'
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
