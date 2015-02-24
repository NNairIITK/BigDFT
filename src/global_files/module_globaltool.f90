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
    public :: read_and_merge_data
    public :: init_gt_data
    public :: finalize_gt_data
    public :: write_merged
    public :: write_transitionpairs
    public :: unpair

    type gt_uinp
        real(gp) :: en_delta
        real(gp) :: fp_delta
        integer  :: ndir
        character(len=500), allocatable :: directories(:)
        character(len=3) :: fptype

        integer  :: ntranspairs
        logical  :: search_transpairs
        character(len=600), allocatable :: trans_pairs_paths(:,:)
        real(gp), allocatable :: fp_arr_trans_pairs(:,:,:)
        real(gp), allocatable :: en_arr_trans_pairs(:,:)
    end type

    type gt_data
        logical :: oldfilename
        integer :: nid
        integer :: nat
        integer :: ntrans
        integer :: nmin
        integer :: ntransmax
        integer :: nminmax   !number of poslocs over all directories
        integer, allocatable :: nminpd(:) !number of minima per directory
        integer :: nminmaxpd !maximum entry of nminpd
        type(atomic_structure) :: astruct
        type(gt_uinp) :: uinp
        real(gp), allocatable :: rcov(:)
        logical, allocatable :: input_transpair_found(:)
        character(len=600), allocatable :: trans_pairs_paths_found(:,:)

        real(gp), allocatable :: fp_arr(:,:)
        real(gp), allocatable :: en_arr(:)
        character(len=600), allocatable :: path_min(:)

        integer :: nposlocs
        real(gp), allocatable :: fp_arr_currDir(:,:)
        real(gp), allocatable :: en_arr_currDir(:)
        character(len=600), allocatable :: path_min_currDir(:)

        integer, allocatable  :: transpairs(:)!Use Cantors pairing function
                                              !to identify pairs
        integer, allocatable  :: minnumber(:)
        integer, allocatable  :: sadnumber(:)

        real(gp), allocatable :: gmon_ener(:)
        real(gp), allocatable :: gmon_fp(:,:)
        character(len=600), allocatable :: gmon_path(:)
        character(len=1), allocatable :: gmon_stat(:)
        integer :: gmon_nposlocs !nposlocs according to gmon file
                                 !can be different to nposlocs
                                 !above due to kills of global executable
                                 !(new poslocm might be written, but not
                                 !global.mon)

        integer, allocatable :: mn(:)
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
    do idict=1,gdat%uinp%ndir
        ifile=0
        call check_filename(gdat,idict)
        do
            call construct_filename(gdat,idict,ifile+1,filename)
            call check_struct_file_exists(filename,exists=exists)
            if(.not. exists)exit
            ifile=ifile+1
        enddo
        call yaml_map('Number of poslocms in '//trim(adjustl(&
             gdat%uinp%directories(idict))),ifile)
        gdat%nminpd(idict)=ifile
        gdat%nminmax=gdat%nminmax+ifile
    enddo
    call yaml_map('Total poslocms found',gdat%nminmax)
end subroutine count_poslocm
!=====================================================================
subroutine construct_filename(gdat,idict,ifile,filename)
    implicit none
    !parameters
    type(gt_data), intent(in) :: gdat
    integer, intent(in) :: idict
    integer, intent(in) :: ifile
    character(len=600), intent(out) :: filename

    if(.not.gdat%oldfilename)then
        !for bigdft >= 1.7.6
        write(filename,'(a,i4.4)')trim(adjustl(&
             gdat%uinp%directories(idict)))//'/poslocm_',ifile
    else
        !for bigdft < 1.7.6
        write(filename,'(a,i4.4,a)')trim(adjustl(&
             gdat%uinp%directories(idict)))//'/poslocm_',&
             ifile,'_'
    endif
end subroutine construct_filename
!=====================================================================
subroutine init_nat_rcov(gdat)
    use module_base
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    implicit none
    !parameter
    type(gt_data), intent(inout) :: gdat
    !local
    character(len=600) :: filename
    call check_filename(gdat,1)
    call construct_filename(gdat,1,1,filename)
    call deallocate_atomic_structure(gdat%astruct)
    call set_astruct_from_file(trim(adjustl(filename)),0,gdat%astruct)
    gdat%nat=gdat%astruct%nat
    gdat%rcov = f_malloc((/gdat%nat/),id='rcov')
    call give_rcov(0,gdat%astruct,gdat%rcov)
end subroutine init_nat_rcov
!=====================================================================
subroutine init_gt_data(gdat)
    use module_base
    use module_atoms, only: nullify_atomic_structure
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local

    gdat%ntrans=0
    gdat%nmin=0
    gdat%ntransmax=0
    gdat%nminmax=0
    gdat%nminmaxpd=0
    gdat%oldfilename=.false.

    call nullify_atomic_structure(gdat%astruct)

    gdat%nminpd = f_malloc((/gdat%uinp%ndir/),id='nminpd')

    call count_poslocm(gdat)
    gdat%nminmaxpd = maxval(gdat%nminpd)


    gdat%fp_arr = f_malloc((/gdat%nid,gdat%nminmax/),id='fp_arr')
    gdat%en_arr = f_malloc((/gdat%nminmax/),id='en_arr')
    gdat%en_arr=huge(1.0_gp)
    gdat%path_min =f_malloc_str(600,(/1.to.gdat%nminmax/),id='path_min')
    gdat%fp_arr_currDir = f_malloc((/gdat%nid,gdat%nminmaxpd/),&
                          id='fp_arr_currDir')
    gdat%en_arr_currDir = f_malloc((/gdat%nminmaxpd/),id='en_arr_currDir')
    gdat%en_arr_currDir=huge(1.0_gp)
    gdat%path_min_currDir = f_malloc_str(600,(/1.to.gdat%nminmaxpd/),&
                            id='path_min_currDir')
    gdat%transpairs = f_malloc((/gdat%nminmax/),id='transpairs')
    gdat%transpairs=huge(1)
    gdat%minnumber = f_malloc((/gdat%nminmax/),id='minnumber')
    gdat%sadnumber = f_malloc((/gdat%nminmax/),id='sadnumber')
    gdat%gmon_ener = f_malloc((/gdat%nminmaxpd/),id='gmon_ener')
    gdat%gmon_fp = f_malloc((/gdat%nid,gdat%nminmaxpd/),id='gmon_fp')
    gdat%gmon_path = f_malloc_str(600,(/1.to.gdat%nminmaxpd/),&
                            id='gmon_path')
    gdat%gmon_stat = f_malloc_str(1,(/gdat%nminmaxpd/),id='gmon_stat')
    gdat%gmon_nposlocs=huge(1)
end subroutine init_gt_data
!=====================================================================
subroutine finalize_gt_data(gdat)
    use module_atoms, only: deallocate_atomic_structure
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat

    call deallocate_atomic_structure(gdat%astruct)
    call f_free(gdat%nminpd)
    call f_free_str(500,gdat%uinp%directories)
!    call f_free_str(600,gdat%uinp%trans_pairs_paths)
    if(gdat%uinp%search_transpairs)then
    deallocate(gdat%uinp%trans_pairs_paths)
    deallocate(gdat%trans_pairs_paths_found)
    endif
    call f_free(gdat%input_transpair_found)
!    call f_free(gdat%trans_pairs_paths_found)
    call f_free(gdat%uinp%fp_arr_trans_pairs)
    call f_free(gdat%uinp%en_arr_trans_pairs)
    call f_free(gdat%fp_arr)
    call f_free(gdat%en_arr)
    call f_free_str(600,gdat%path_min)
    call f_free(gdat%fp_arr_currDir)
    call f_free(gdat%en_arr_currDir)
    call f_free_str(600,gdat%path_min_currDir)
    call f_free(gdat%rcov)
    call f_free(gdat%transpairs)
    call f_free(gdat%minnumber)
    call f_free(gdat%sadnumber)
    call f_free(gdat%gmon_ener)
    call f_free(gdat%gmon_fp)
    call f_free_str(600,gdat%gmon_path)
    call f_free(gdat%mn)
    call f_free_str(1,gdat%gmon_stat)
end subroutine finalize_gt_data
!=====================================================================
subroutine read_globaltool_uinp(gdat)
    use module_base
    use yaml_output
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    integer :: u, istat, idict, iline, idx
    character(len=2000) :: line
    character(len=1024) :: comment
    real(gp) :: energy
    real(gp) :: rdmy
    logical :: exists

    inquire(file='globaltool.inp',exist=exists)
    if(.not. exists)then
        call f_err_throw('The input file globaltool.inp does not exist.',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif

    u=f_get_free_unit()
    open(u,file='globaltool.inp')
    read(u,*,iostat=istat)gdat%uinp%en_delta,gdat%uinp%fp_delta,gdat%uinp%fptype
    if(istat/=0)then
        call f_err_throw('Error in first line of globaltool.inp',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
    gdat%uinp%ndir=0
    do
        read(u,'(a)',iostat=istat)line
        if(istat/=0)exit
        gdat%uinp%ndir=gdat%uinp%ndir+1
    enddo
    gdat%uinp%directories = f_malloc_str(500,(/1.to.gdat%uinp%ndir/),&
                            id='gdat%uinp%directories')
    rewind(u)
    read(u,*,iostat=istat)gdat%uinp%en_delta,gdat%uinp%fp_delta
    do idict=1,gdat%uinp%ndir
        read(u,'(a)')gdat%uinp%directories(idict)
    enddo
    close(u)

    call init_nat_rcov(gdat)
    if(trim(adjustl(gdat%uinp%fptype))=='SO')then
        gdat%nid = gdat%nat !s-overlap
    elseif(trim(adjustl(gdat%uinp%fptype))=='SPO')then
        gdat%nid = 4*gdat%nat !sp-overlap
    else
        call f_err_throw('Fingerprint type "'//&
             trim(adjustl(gdat%uinp%fptype))//'" unknown.',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif

    !read transition pairs for which MH aligned version should
    !be found
    gdat%uinp%ntranspairs=0
    inquire(file='transition_pairs.inp',exist=gdat%uinp%search_transpairs)
    if(gdat%uinp%search_transpairs)then
        open(u,file='transition_pairs.inp')
        do
            read(u,'(a)',iostat=istat)line
            if(istat/=0)exit
            gdat%uinp%ntranspairs = gdat%uinp%ntranspairs + 1
        enddo
!        gdat%uinp%trans_pairs_paths = f_malloc_str(600,(/2,1.to.gdat%uinp%ntranspairs/),&
!                            id='gdat%uinp%directories')
        allocate(gdat%uinp%trans_pairs_paths(2,gdat%uinp%ntranspairs))
        gdat%uinp%fp_arr_trans_pairs = f_malloc((/gdat%nid,2,gdat%uinp%ntranspairs/),&
                                       id='fp_arr_trans_pairs')
        gdat%uinp%en_arr_trans_pairs = f_malloc((/2,gdat%uinp%ntranspairs/),&
                                       id='en_arr_trans_pairs')
!        gdat%trans_pairs_paths_found = f_malloc_str(600,(/2,1.to.gdat%uinp%ntranspairs/),&
!                            id='gdat%trans_pairs_paths_found')
        allocate(gdat%trans_pairs_paths_found(2,gdat%uinp%ntranspairs))
        gdat%input_transpair_found = f_malloc((/1.to.gdat%uinp%ntranspairs/),&
                            id='gdat%input_transpair_found')
        call yaml_map('Number of input transition pairs for which MH aligned versions are looked for',gdat%uinp%ntranspairs)
        rewind(u)
        do iline=1,gdat%uinp%ntranspairs
            read(u,'(a)')line
            idx=index(line," ")
            gdat%uinp%trans_pairs_paths(1,iline)=trim(adjustl(line(1:idx)))
            gdat%uinp%trans_pairs_paths(2,iline)=trim(adjustl(line(idx:2000)))
            call deallocate_atomic_structure(gdat%astruct)
            call set_astruct_from_file(trim(adjustl(gdat%uinp%trans_pairs_paths(1,iline))),0,&
                 gdat%astruct,comment=comment,energy=energy)
            if (gdat%astruct%nat /= gdat%nat) then
               call f_err_throw('Error in read_globaltool_uinp:'//&
                ' wrong nat ('//trim(yaml_toa(gdat%astruct%nat))&
                //' /= '//trim(yaml_toa(gdat%nat))//')',&
                err_name='BIGDFT_RUNTIME_ERROR')
            end if
            gdat%uinp%en_arr_trans_pairs(1,iline) = energy
            call fingerprint(gdat%nat,gdat%nid,gdat%astruct%cell_dim,&
                 gdat%astruct%geocode,gdat%rcov,gdat%astruct%rxyz,&
                 gdat%uinp%fp_arr_trans_pairs(1,1,iline))
            call deallocate_atomic_structure(gdat%astruct)
            call set_astruct_from_file(trim(adjustl(gdat%uinp%trans_pairs_paths(2,iline))),0,&
                 gdat%astruct,comment=comment,energy=energy)
            if (gdat%astruct%nat /= gdat%nat) then
               call f_err_throw('Error in read_globaltool_uinp:'//&
                ' wrong nat ('//trim(yaml_toa(gdat%astruct%nat))&
                //' /= '//trim(yaml_toa(gdat%nat))//')',&
                err_name='BIGDFT_RUNTIME_ERROR')
            end if
            gdat%uinp%en_arr_trans_pairs(2,iline) = energy
            call fingerprint(gdat%nat,gdat%nid,gdat%astruct%cell_dim,&
                 gdat%astruct%geocode,gdat%rcov,gdat%astruct%rxyz,&
                 gdat%uinp%fp_arr_trans_pairs(1,2,iline))
            gdat%input_transpair_found(iline)=.false.
        enddo
        close(u)
    endif

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
    call yaml_map('Number of directories',gdat%uinp%ndir)
    call yaml_mapping_open('directories')
    do idict=1,gdat%uinp%ndir
    call yaml_scalar(trim(adjustl(gdat%uinp%directories(idict))))
    enddo
    call yaml_mapping_close()
end subroutine write_globaltool_uinp
!=====================================================================
subroutine read_poslocs(gdat,idict)
    use module_base
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    integer, intent(in) :: idict
    !local
    integer :: iposloc, ifile, indx
    character(len=600) :: filename
    logical :: exists
    character(len=1024) :: comment
    real(gp) :: energy

    ifile=0
    iposloc=0
    do
        call construct_filename(gdat,idict,ifile+1,filename)
        call check_struct_file_exists(filename,exists=exists)
        if(.not. exists)exit
        ifile=ifile+1
        call deallocate_atomic_structure(gdat%astruct)
        call set_astruct_from_file(trim(adjustl(filename)),0,&
             gdat%astruct,comment=comment,energy=energy)
        if (gdat%astruct%nat /= gdat%nat) then
            call f_err_throw('Error in read_poslocs:'//&
                 ' wrong nat ('//trim(yaml_toa(gdat%astruct%nat))&
                 //' /= '//trim(yaml_toa(gdat%nat))//')',&
                 err_name='BIGDFT_RUNTIME_ERROR')
        end if

        indx=index(comment,'fnrm')
        if(indx==0)cycle
        iposloc=iposloc+1
        gdat%en_arr_currDir(iposloc) = energy
        call fingerprint(gdat%nat,gdat%nid,gdat%astruct%cell_dim,&
             gdat%astruct%geocode,gdat%rcov,gdat%astruct%rxyz,&
             gdat%fp_arr_currDir(1,iposloc))
        gdat%path_min_currDir(iposloc)=trim(adjustl(filename))//'.'//&
                          trim(adjustl(gdat%astruct%inputfile_format))
    enddo
    gdat%nposlocs=iposloc
    gdat%ntransmax=gdat%nposlocs

end subroutine read_poslocs
!=====================================================================
subroutine write_merged(gdat)
    use module_base
    use yaml_output
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    integer :: imin
    call yaml_comment('Merged minima ....',hfill='-')
    gdat%mn = f_malloc((/1.to.gdat%nmin/),id='gdat%mn')
    do imin=1,gdat%nmin
        gdat%mn(gdat%minnumber(imin)) = imin
        write(*,'(i6.6,1x,es24.17,1x,a)')gdat%minnumber(imin),gdat%en_arr(imin),&
                                trim(adjustl(gdat%path_min(imin)))
    enddo
end subroutine write_merged
!=====================================================================
subroutine write_transitionpairs(gdat)
    use module_base
    use yaml_output
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(in) :: gdat
    !local
    integer :: itrans
    integer :: IDmin1, IDmin2
    integer :: kIDmin1, kIDmin2
    real(gp) :: fpd
    call yaml_comment('Transition pairs unified ....',hfill='-')
    write(*,'(a)')'  #Trans IDmin1 IDmin2  Ener1                '//&
         '    Ener2                    |DeltaEner|         '//&
         '     FPdist'
    do itrans=1,gdat%ntrans
        call unpair(gdat%transpairs(itrans),IDmin1,IDmin2)
        kIDmin1=gdat%mn(IDmin1)
        kIDmin2=gdat%mn(IDmin2)
        call fpdistance(gdat%nid,gdat%fp_arr(1,kIDmin1),&
             gdat%fp_arr(1,kIDmin2),fpd)
        write(*,'(a,1x,i4.4,3x,i4.4,2x,4(1x,es24.17),1x,i8.8)')'   Trans',&
             IDmin1,IDmin2,gdat%en_arr(kIDmin1),gdat%en_arr(kIDmin2),&
             abs(gdat%en_arr(kIDmin1)-gdat%en_arr(kIDmin2)),fpd,gdat%transpairs(itrans)
    enddo

    if(gdat%uinp%search_transpairs)then
        call yaml_comment('Identified transition pairs (transition_pairs.inp)  ....',hfill='-')
        do itrans=1,gdat%uinp%ntranspairs
            if(.not.gdat%input_transpair_found(itrans))then
                call f_err_throw('Transition pair no. '//yaml_toa(itrans)//&
                     ' not found',err_name='BIGDFT_RUNTIME_ERROR')
            endif
            write(*,'(2(5x,a))')trim(adjustl(gdat%trans_pairs_paths_found(1,itrans))),&
                                trim(adjustl(gdat%trans_pairs_paths_found(2,itrans)))
        enddo
    endif
end subroutine write_transitionpairs
!=====================================================================
subroutine read_and_merge_data(gdat)
    use module_base
    use yaml_output
    use module_atoms, only: set_astruct_from_file,&
                            deallocate_atomic_structure
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    integer :: idict

    call yaml_comment('Merging poslocms ....',hfill='-')
    do idict =1, gdat%uinp%ndir
        call check_filename(gdat,idict)
        call read_poslocs(gdat,idict)
        call add_poslocs_to_database(gdat)
        call read_globalmon(gdat,idict)
        call add_transpairs_to_database(gdat)
    enddo

end subroutine read_and_merge_data
!=====================================================================
subroutine check_filename(gdat,idict)
    use module_base
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    integer, intent(in) :: idict
    !local
    character(len=600) :: filename
    logical :: exists
    !for bigdft >= 1.7.6
    write(filename,'(a,i4.4)')trim(adjustl(&
         gdat%uinp%directories(idict)))//'/poslocm_',1
    call check_struct_file_exists(trim(adjustl(filename)),exists)
    if(exists)then
        gdat%oldfilename=.false.
        return
    endif
    !for bigdft < 1.7.6
    write(filename,'(a,i4.4,a)')trim(adjustl(&
         gdat%uinp%directories(idict)))//'/poslocm_',&
         1,'_'
    call check_struct_file_exists(trim(adjustl(filename)),exists)
    if(exists)then
        gdat%oldfilename=.true.
        return
    endif

    if(.not. exists)then
        call f_err_throw(trim(adjustl(filename))//' does not exist.',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
end subroutine check_filename
!=====================================================================
subroutine add_transpairs_to_database(gdat)
    use module_base
    use yaml_output
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    real(gp) :: ecurr
    real(gp) :: fpcurr(gdat%nid)
    character(len=1) :: statcurr
    integer :: idcurr, idnext
    integer :: kidcurr, kidnext
    integer :: kid, k_epot
    logical :: lnew
    integer :: id_transpair
    integer :: iposloc
    integer :: iposloc_curr, iposloc_next
    integer :: loc_id_transpair
    integer :: i
    real(gp) :: fpd

    call yaml_comment('Reconstructing transition pairs ....',hfill='-')
    write(*,'(a)')'  #trans IDmin1 IDmin2  Ener1                '//&
         '    Ener2                    |DeltaEner|         '//&
         '     FPdist'
    if(gdat%gmon_stat(1)/='P')then
        call f_err_throw('Error in global.mon: Does not start with'//&
             ' P line',err_name='BIGDFT_RUNTIME_ERROR')
    endif
    do iposloc = 1, gdat%gmon_nposlocs
        !check if escaped
        if(gdat%gmon_stat(iposloc)=='S')cycle
        !get ID of minimum
        call identical('min',gdat,gdat%nminmax,gdat%nmin,gdat%nid,&
             gdat%gmon_ener(iposloc),gdat%gmon_fp(1,iposloc),&
             gdat%en_arr,gdat%fp_arr,gdat%uinp%en_delta,&
             gdat%uinp%fp_delta,lnew,kid,k_epot)
        if(lnew)then
            call f_err_throw('Structure does not exist',&
                 err_name='BIGDFT_RUNTIME_ERROR')
        endif
        kidnext=kid
        idnext = gdat%minnumber(kid)
        iposloc_next = iposloc
        !check if restart happened
        if(gdat%gmon_stat(iposloc)=='P')then
            kidcurr = kidnext
            idcurr  = idnext
            iposloc_curr = iposloc
            ecurr = gdat%gmon_ener(iposloc)
            fpcurr(:) = gdat%gmon_fp(:,iposloc)
            statcurr = gdat%gmon_stat(iposloc)
            cycle
        endif
        call check_given_pairs(gdat,iposloc_curr,iposloc_next)
        id_transpair = getPairId(idcurr,idnext)
        call fpdistance(gdat%nid,gdat%fp_arr(1,kidcurr),&
             gdat%fp_arr(1,kidnext),fpd)
        write(*,'(a,1x,i4.4,3x,i4.4,2x,4(1x,es24.17))')'   trans',idcurr,&
             idnext,gdat%en_arr(kidcurr),gdat%en_arr(kidnext),&
             abs(gdat%en_arr(kidcurr)-gdat%en_arr(kidnext)),fpd

        call inthunt_gt(gdat%transpairs,&
             max(1,min(gdat%ntrans,gdat%nminmax)),id_transpair,&
             loc_id_transpair)
        !comment the if query if everey pair should be added to the
        !database, even if it is already in the database
        if(gdat%transpairs(max(1,loc_id_transpair))/=id_transpair)then!add to database
            !shift
            gdat%ntrans=gdat%ntrans+1
            do i=gdat%ntrans-1,loc_id_transpair+1,-1
                gdat%transpairs(i+1)=gdat%transpairs(i)
            enddo
            gdat%transpairs(loc_id_transpair+1) = id_transpair
        endif
        !check if accepted
        if(gdat%gmon_stat(iposloc)=='A')then
            kidcurr = kidnext
            idcurr = idnext
            iposloc_curr = iposloc
            ecurr = gdat%gmon_ener(iposloc)
            fpcurr(:) = gdat%gmon_fp(:,iposloc)
            statcurr = gdat%gmon_stat(iposloc)
        endif
    enddo
end subroutine add_transpairs_to_database
!=====================================================================
subroutine check_given_pairs(gdat,iposloc_curr,iposloc_next)
    use module_base
    use yaml_output
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    integer, intent(in) :: iposloc_curr, iposloc_next
    !loal
    integer :: itrans

    do itrans=1,gdat%uinp%ntranspairs
        if(.not.gdat%input_transpair_found(itrans))then
            !check icurr
            if (equal(gdat%nid,gdat%uinp%en_arr_trans_pairs(1,itrans),&
                gdat%gmon_ener(iposloc_curr),&
                gdat%uinp%fp_arr_trans_pairs(1,1,itrans),&
                gdat%gmon_fp(1,iposloc_curr),&
                gdat%uinp%en_delta,gdat%uinp%fp_delta))then
                if(equal(gdat%nid,gdat%uinp%en_arr_trans_pairs(2,itrans),&
                   gdat%gmon_ener(iposloc_next),&
                   gdat%uinp%fp_arr_trans_pairs(1,2,itrans),&
                   gdat%gmon_fp(1,iposloc_next),&
                   gdat%uinp%en_delta,gdat%uinp%fp_delta))then

                    gdat%trans_pairs_paths_found(1,itrans)=gdat%gmon_path(iposloc_curr)
                    gdat%trans_pairs_paths_found(2,itrans)=gdat%gmon_path(iposloc_next)
                    gdat%input_transpair_found(itrans)=.true.
                endif
            else if((equal(gdat%nid,gdat%uinp%en_arr_trans_pairs(2,itrans),&
                gdat%gmon_ener(iposloc_curr),&
                gdat%uinp%fp_arr_trans_pairs(1,2,itrans),&
                gdat%gmon_fp(1,iposloc_curr),&
                gdat%uinp%en_delta,gdat%uinp%fp_delta))) then
                if(equal(gdat%nid,gdat%uinp%en_arr_trans_pairs(1,itrans),&
                   gdat%gmon_ener(iposloc_next),&
                   gdat%uinp%fp_arr_trans_pairs(1,1,itrans),&
                   gdat%gmon_fp(1,iposloc_next),&
                   gdat%uinp%en_delta,gdat%uinp%fp_delta))then

                    gdat%trans_pairs_paths_found(1,itrans)=gdat%gmon_path(iposloc_next)
                    gdat%trans_pairs_paths_found(2,itrans)=gdat%gmon_path(iposloc_curr)
                    gdat%input_transpair_found(itrans)=.true.
                endif
            endif
        endif
    enddo
end subroutine check_given_pairs
!=====================================================================
subroutine read_globalmon(gdat,idict)
    use module_base
    use yaml_output
    !reads the global.mon file and
    !associates each line with the
    !corresponding structure from the
    !polocm files
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    integer, intent(in) :: idict
    !local
    character(len=600) :: filename
    character(len=200) :: line
    character(len=1)   :: stat
    logical :: found
    integer :: u
    integer :: icount
    integer :: iline
    integer :: istat
    integer :: restartoffset
    integer :: istep
    real(gp) :: ristep
    real(gp) :: energy
    real(gp) :: rdmy
    real(gp) :: fp(gdat%nid)
integer :: itmp


    u=f_get_free_unit()
    filename=trim(adjustl(gdat%uinp%directories(idict)))//'/data/global.mon'
    call yaml_comment('Parsing '//trim(adjustl(filename))//' ....',&
         hfill='-')
    open(u,file=trim(adjustl(filename)))
    icount=0
    iline=0
    restartoffset=0
    do
        read(u,'(a)',iostat=istat)line
        if(istat/=0)exit
        iline=iline+1
        read(line,*)ristep
        istep=nint(ristep)
        if(istep/=0)then
            read(line,*,iostat=istat)ristep,energy,rdmy,rdmy,rdmy,rdmy,rdmy,stat
        else
            read(line,*,iostat=istat)ristep,energy,rdmy,rdmy,stat
            if(icount/=0)restartoffset=icount
        endif
write(*,*)stat
        if(istat/=0)then
            call f_err_throw('Error while parsing '//&
                 trim(adjustl(filename))//', istep='//&
                 trim(yaml_toa(istep)),err_name='BIGDFT_RUNTIME_ERROR')
        endif
        !check possible states
        if(stat/="P" .and. stat/="I" .and. stat/="S" .and.&
           stat/="A" .and. stat/="R")then
            call f_err_throw('Status "'//trim(adjustl(stat))//&
                 '" is unknown. Has format of global.mon been changed?'//&
                 ' If so, update parsing routine',&
                 err_name='BIGDFT_RUNTIME_ERROR')
        endif
        if(stat/="I")then
            if(icount/=istep+restartoffset)then
                call f_err_throw('Error while parsing '//&
                     trim(adjustl(filename))//', istep+restartoffset='//&
                     trim(yaml_toa(istep+restartoffset))//' icount='//&
                     trim(yaml_toa(icount)),err_name='BIGDFT_RUNTIME_ERROR')
            endif
            icount=icount+1
write(*,*)'ins',icount,energy
            gdat%gmon_ener(icount)=energy
            gdat%gmon_stat(icount)=stat
            call gmon_line_to_fp(gdat,icount,iline,found)
            if(.not. found)then
                write(*,*)'to be found:',energy,icount
                do itmp = icount+10 , 1 , -1
                    write(*,*)itmp,gdat%en_arr_currDir(itmp)
                enddo

                call f_err_throw('Could not find istep= '//&
                     trim(yaml_toa(istep))//'of '//trim(adjustl(filename))//&
                     ' among the poslocm files.',err_name='BIGDFT_RUNTIME_ERROR')
            endif
        endif
    enddo
    if(icount/=gdat%nposlocs)then
        !might happen if global was killed after writing the first
        ! poslocm file and before writing global.mon
        ! -> Do not stop, but print warning
!        call f_err_throw('Number of poslocm files is not identical '//&
!             'to number of steps in global.mon',err_name='BIGDFT_RUNTIME_ERROR')
        call yaml_warning('Number of poslocm files in '//&
             trim(adjustl(gdat%uinp%directories(idict)))//&
             ' is not identical to number of steps in '//&
             trim(adjustl(filename)))
    endif
    gdat%gmon_nposlocs=icount
    close(u)

end subroutine read_globalmon
!=====================================================================
subroutine gmon_line_to_fp(gdat,icount,iline,found)
    use module_base
    use yaml_output
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    integer, intent(in) :: icount
    integer, intent(in) :: iline
    logical, intent(out) :: found
    !local
    integer :: iposloc
    !ethresh can be small, since identical energy value should be in
    !global.mon and in the corresponding poslocm file
    real(gp), parameter :: ethresh = 1.d-12

    found=.false.
    !the corresponding poslocm file should always
    !be one of the last files. Therefore, read from behind
    !We start at icount+10, because glbal might have been killed
    !after writing the first poslocm file and before writing global.mon
    do iposloc = min(icount+10,gdat%nminmaxpd) , 1 , -1
write(*,'(i4.4,2(1x,es24.17))')iposloc,gdat%en_arr_currDir(iposloc),abs(gdat%en_arr_currDir(iposloc)-gdat%gmon_ener(icount))
        if(abs(gdat%en_arr_currDir(iposloc)-&
           gdat%gmon_ener(icount))<ethresh)then
            found=.true.
            call yaml_scalar('Line '//trim(yaml_toa(iline))//&
                 ' corresponds to '//&
                 trim(adjustl(gdat%path_min_currDir(iposloc))))
            gdat%gmon_fp(:,icount) = gdat%fp_arr_CurrDir(:,iposloc)
            gdat%gmon_path(icount) = gdat%path_min_currDir(iposloc)
            exit
        endif
    enddo

end subroutine gmon_line_to_fp
!=====================================================================
subroutine add_poslocs_to_database(gdat)
    use module_base
    use yaml_output
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    !local
    integer :: iposloc
    integer :: kid
    integer :: k_epot
    logical :: lnew

    do iposloc=1,gdat%nposlocs
       call identical('minIns',gdat,gdat%nminmax,gdat%nmin,gdat%nid,&
            gdat%en_arr_currDir(iposloc),gdat%fp_arr_currDir(1,iposloc),&
            gdat%en_arr,gdat%fp_arr,gdat%uinp%en_delta,&
            gdat%uinp%fp_delta,lnew,kid,k_epot)
       if(lnew)then
           call yaml_comment('Minimum'//trim(adjustl(&
                gdat%path_min_currDir(iposloc)))//' is new.')
           call insert_min(gdat,k_epot,gdat%en_arr_currDir(iposloc),&
                gdat%fp_arr_currDir(1,iposloc),&
                gdat%path_min_currDir(iposloc))
       else
           call yaml_comment('Minimum '//trim(adjustl(&
                gdat%path_min_currDir(iposloc)))//&
                ' is identical to minimum '//&
                trim(yaml_toa(gdat%minnumber(kid))))
       endif

    enddo
end subroutine add_poslocs_to_database
!=====================================================================
subroutine identical(cf,gdat,ndattot,ndat,nid,epot,fp,en_arr,fp_arr,en_delta,&
                    fp_delta,lnew,kid,k_epot)
    use module_base
    use yaml_output
    use module_fingerprints
    implicit none
    !parameters
    type(gt_data), intent(in) :: gdat
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
    character(len=*), intent(in) :: cf
    !local
    integer :: k, klow, khigh, nsm
    real(gp) :: dmin, d
    !search in energy array
    call hunt_gt(en_arr,max(1,min(ndat,ndattot)),epot,k_epot)
    lnew=.true.

    ! find lowest configuration that might be identical
    klow=k_epot
    do k=k_epot,1,-1
        if (epot-en_arr(k).lt.0.0_gp) stop 'zeroA'
        if (epot-en_arr(k).gt.en_delta) exit
        klow=k
    enddo

    ! find highest  configuration that might be identical
    khigh=k_epot+1
    do k=k_epot+1,ndat
        if (en_arr(k)-epot.lt.0.0_gp) stop 'zeroB'
        if (en_arr(k)-epot.gt.en_delta) exit
        khigh=k
    enddo

    nsm=0
    dmin=huge(1.e0_gp)
    do k=max(1,klow),min(ndat,khigh)
        if (abs(epot-en_arr(k)).le.en_delta) then
            call fpdistance(nid,fp,fp_arr(1,k),d)
            write(*,*)'fpdist '//trim(adjustl(cf)),abs(en_arr(k)-epot),d
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
subroutine insert_min(gdat,k_epot,epot,fp,path)
    !insert at k_epot+1
    use module_base
    implicit none
    !parameters
    type(gt_data), intent(inout) :: gdat
    integer, intent(in) :: k_epot
    real(gp), intent(in) :: epot
    real(gp), intent(in) :: fp(gdat%nid)
    character(len=600)   :: path
    !local
    integer :: i,k
write(*,*)'insert at',k_epot+1
    if(gdat%nmin+1>gdat%nminmax)stop 'nmin+1>=nminmax, out of bounds'

    gdat%nmin=gdat%nmin+1
    do k=gdat%nmin-1,k_epot+1,-1
        gdat%en_arr(k+1)=gdat%en_arr(k)
        gdat%minnumber(k+1)=gdat%minnumber(k)
        gdat%path_min(k+1)=gdat%path_min(k)
        do i=1,gdat%nid
            gdat%fp_arr(i,k+1)=gdat%fp_arr(i,k)
         enddo
    enddo
    gdat%en_arr(k_epot+1)=epot
    gdat%minnumber(k_epot+1)=gdat%nmin
    gdat%path_min(k_epot+1)=path
    do i=1,gdat%nid
        gdat%fp_arr(i,k_epot+1)=fp(i)
    enddo
end subroutine insert_min
!!!!=====================================================================
!!!!> C x is in interval [xx(jlo),xx(jlow+1)
!!!![ ; xx(0)=-Infinity ; xx(n+1) = Infinity
subroutine hunt_gt(xx,n,x,jlo)
  use module_base
  implicit none
  !Arguments
  integer :: jlo,n
  real(gp) :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  include 'hunt-inc.f90'
END SUBROUTINE hunt_gt
!=====================================================================
subroutine inthunt_gt(xx,n,x,jlo)
  use module_base
  implicit none
  !Arguments
  integer :: jlo,n
  integer :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  include 'hunt-inc.f90'
END SUBROUTINE inthunt_gt
!=====================================================================
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
                 ' does not exist.',err_name='BIGDFT_RUNTIME_ERROR')
        endif
    endif
end subroutine check_struct_file_exists
!=====================================================================
subroutine give_rcov(iproc,astruct,rcov)
  use module_base
  use yaml_output
  use module_atoms
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(atomic_structure), intent(in) :: astruct
  real(kind=gp), dimension(astruct%nat), intent(out) :: rcov
  !Local variables
  integer :: iat

  do iat=1,astruct%nat
    call covalent_radius(astruct%atomnames(astruct%iatype(iat)),&
         rcov(iat))
    if (iproc == 0) call yaml_map('(MH) RCOV '//&
                         trim(astruct%atomnames(astruct%iatype(iat))),&
                         rcov(iat))
  end do
  !python metod
!  it=atoms_iter(astruct)
!  do while(atoms_iter_next(it))
!    call covalent_radius(it%name,rcov(it%iat))
!  end do

end subroutine give_rcov
!=====================================================================
logical function equal(nid,ener1,ener2,fp1,fp2,en_delta,fp_delta)
    use module_base
    use module_fingerprints
    implicit none
    !parameter
    integer, intent(in) :: nid
    real(gp), intent(in) :: ener1, ener2
    real(gp), intent(in) :: fp1(nid),fp2(nid)
    real(gp), intent(in) :: en_delta,fp_delta
    !local
    real(gp) :: d

    equal=.false.
    if (abs(ener1-ener2).le.en_delta) then
        call fpdistance(nid,fp1,fp2,d)
        if(d.le.fp_delta)equal=.true.
    endif
end function equal
end module
