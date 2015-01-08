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
!!    public :: count_saddle_points
    public :: init_gt_data
    public :: finalize_gt_data
    public :: write_merged
    public :: getPairId
    public :: unpair

    type gt_uinp
        real(gp) :: en_delta
        real(gp) :: fp_delta
        integer  :: ndir
        character(len=500), allocatable :: directories(:)
    end type

    type gt_data
        integer :: nid
        integer :: nat
        integer :: ntrans
        integer :: nmin
        integer :: ntransax
        integer :: nminmax   !number of poslocs over all directories
        integer, allocatable :: nminpd(:) !number of minima per directory
        integer :: nminmaxpd !maximum entry of nminpd
        type(atomic_structure) :: astruct
        type(gt_uinp) :: uinp
        real(gp), allocatable :: rcov(:)

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
        character(len=1), allocatable :: gmon_stat(:)
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
end subroutine
!=====================================================================
subroutine construct_filename(gdat,idict,ifile,filename)
    implicit none
    !parameters
    type(gt_data), intent(in) :: gdat
    integer, intent(in) :: idict
    integer, intent(in) :: ifile
    character(len=600), intent(out) :: filename

    !for bigdft >= 1.7.6
    write(filename,'(a,i4.4)')trim(adjustl(&
         gdat%uinp%directories(idict)))//'/poslocm_',ifile
    !for bigdft < 1.7.6
!    write(filename,'(a,i4.4,a)')trim(adjustl(&
!         gdat%uinp%directories(idict)))//'/poslocm_',&
!         ifile,'_'
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
    call construct_filename(gdat,1,1,filename)
    call deallocate_atomic_structure(gdat%astruct)
    call set_astruct_from_file(trim(adjustl(filename)),0,gdat%astruct)
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

    gdat%nminpd = f_malloc((/gdat%uinp%ndir/),id='nminpd')
     
    call count_poslocm(gdat)
    gdat%nminmaxpd = maxval(gdat%nminpd)

    call init_nat_rcov(gdat)
    gdat%nid = gdat%nat !s-overlap
!    nid = 4*nat !sp-overlap


    gdat%fp_arr = f_malloc((/gdat%nid,gdat%nminmax/),id='fp_arr')
    gdat%en_arr = f_malloc((/gdat%nminmax/),id='en_arr')
    gdat%path_min =f_malloc_str(600,(/1.to.gdat%nminmax/),id='path_min')
    gdat%fp_arr_currDir = f_malloc((/gdat%nid,gdat%nminmaxpd/),&
                          id='fp_arr_currDir')
    gdat%en_arr_currDir = f_malloc((/gdat%nminmaxpd/),id='en_arr_currDir')
    gdat%path_min_currDir = f_malloc_str(600,(/1.to.gdat%nminmaxpd/),&
                            id='path_min_currDir')
    gdat%transpairs = f_malloc((/gdat%nminmax/),id='transpairs')
    gdat%minnumber = f_malloc((/gdat%nminmax/),id='minnumber')
    gdat%sadnumber = f_malloc((/gdat%nminmax/),id='sadnumber')
    gdat%gmon_ener = f_malloc((/gdat%nminmaxpd/),id='gmon_ener')
    gdat%gmon_fp = f_malloc((/gdat%nid,gdat%nminmaxpd/),id='gmon_fp')
    gdat%gmon_stat = f_malloc_str(100,(/gdat%nminmaxpd/),id='gmon_stat')
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
    call f_free_str(600,gdat%gmon_stat)
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
    gdat%uinp%ndir=0
    read(u,*,iostat=istat)gdat%uinp%en_delta,gdat%uinp%fp_delta
    if(istat/=0)then
        call f_err_throw('Error in first line of globaltool.inp',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
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
     
end subroutine read_poslocs
!=====================================================================
subroutine write_merged(gdat)
    use module_base
    use yaml_output
    implicit none
    !parameters
    type(gt_data), intent(in) :: gdat
    !local
    integer :: imin
    call yaml_comment('Merged minima ....',hfill='-')
    do imin=1,gdat%nmin
        write(*,'(i6.6,1x,es24.17,1x,a)')imin,gdat%en_arr(imin),&
                                trim(adjustl(gdat%path_min(imin)))
    enddo
end subroutine write_merged
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

    call yaml_comment('Merging ....',hfill='-')
    do idict =1, gdat%uinp%ndir
        call read_poslocs(gdat,idict)
        call add_poslocs_to_database(gdat)
        call read_globalmon(gdat,idict)
    enddo
    
end subroutine read_and_merge_data
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
    character(len=1)   :: stat
    integer :: u
    integer :: icount
    integer :: iline
    integer :: istat
    integer :: restartoffset
    real(gp) :: ristep
    integer :: istep
    character(len=200) :: line
    real(gp) :: energy
    real(gp) :: rdmy
    real(gp) :: fp(gdat%nid)
    logical :: found

 
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
            read(line,*,iostat=istat)istep,energy,rdmy,rdmy,rdmy,rdmy,rdmy,stat
        else
            read(line,*,iostat=istat)istep,energy,rdmy,rdmy,stat
            if(icount/=0)restartoffset=icount
        endif
        if(istat/=0)then
            call f_err_throw('Error while parsing '//&
                 trim(adjustl(filename))//', istep='//&
                 trim(yaml_toa(istep)))
        endif
        if(stat/="I")then
            if(icount/=istep+restartoffset)then
                call f_err_throw('Error while parsing '//&
                     trim(adjustl(filename))//', istep+restartoffset='//&
                     trim(yaml_toa(istep+restartoffset))//' icount='//&
                     trim(yaml_toa(icount)))
            endif
            icount=icount+1
            gdat%gmon_ener(icount)=energy
            gdat%gmon_stat(icount)=stat
            call gmon_line_to_fp(gdat,icount,iline,found)
            if(.not. found)then
                call f_err_throw('Could not find istep= '//&
                     trim(yaml_toa(istep))//'of '//trim(adjustl(filename))//&
                     ' among the poslocm files.')
            endif    
        endif
    enddo    
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
    integer :: istep
    !ethresh can be small, since identical energy value should be in
    !global.mon and in the corresponding poslocm file
    real(gp), parameter :: ethresh = 1.d-13

    istep=icount-1
    found=.false.
    !the corresponding poslocm file should always
    !be one of the last files. Therefore, read from behind
    do iposloc = icount , 1 , -1
        if(abs(gdat%en_arr_currDir(iposloc)-&
           gdat%gmon_ener(icount))<ethresh)then
            found=.true.
            call yaml_scalar('Line '//trim(yaml_toa(iline))//&
                 ' corresponds to '//&
                 trim(adjustl(gdat%path_min_currDir(iposloc))))
            gdat%gmon_fp(:,icount) = gdat%fp_arr(:,iposloc)
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
       call identical('min',gdat,gdat%nminmax,gdat%nmin,gdat%nid,&
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
    character(len=3), intent(in) :: cf
    !local
    integer :: k, klow, khigh, nsm
    real(gp) :: dmin, d 
    !search in energy array
    call hunt_gt(en_arr,max(1,min(ndat,ndattot)),epot,k_epot)
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
write(*,*)'fpdist '//cf,abs(en_arr(k)-epot),d
!if(cf=='min')then
!if(d<3.d-3)then
!if(abs(en_arr(k)-epot)>1.d-4)then
!write(*,*)trim(adjustl(gdat%path_min(k)))
!    stop
!endif
!endif
!endif
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
!!!subroutine add_neighbors(gdat,kid,neighb1,neighb2)
!!!    use module_base
!!!    implicit none
!!!    !parameters
!!!    type(gt_data), intent(inout) :: gdat
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
!!!    neighbloop: do ipair=1,gdat%nneighbpairs(kid)
!!!        if( ((gdat%sadneighb(1,ipair,kid)==neighb1)&
!!!               .and.(gdat%sadneighb(2,ipair,kid)==neighb2))&
!!!           &.or.((gdat%sadneighb(2,ipair,kid)==neighb1) &
!!!               .and.(gdat%sadneighb(1,ipair,kid)==neighb2)) )then
!!!            gdat%paircounter(ipair,kid) = gdat%paircounter(ipair,kid)+1
!!!            found=.true.
!!!            exit neighbloop
!!!        endif
!!!    enddo neighbloop
!!!
!!!    if(.not. found) then !pair is new, add it to list
!!!        gdat%nneighbpairs(kid) = gdat%nneighbpairs(kid) + 1
!!!        gdat%paircounter(gdat%nneighbpairs(kid),kid)  = 1
!!!        gdat%sadneighb(1,gdat%nneighbpairs(kid),kid) = neighb1
!!!        gdat%sadneighb(2,gdat%nneighbpairs(kid),kid) = neighb2
!!!    endif
!!!do ipair=1,gdat%nneighbpairs(kid)
!!!write(*,*)gdat%sadneighb(1,ipair,kid),gdat%sadneighb(2,ipair,kid),gdat%paircounter(ipair,kid)
!!!enddo
!!!end subroutine add_neighbors
!=====================================================================
!!subroutine insert_sad(gdat,k_epot,epot,fp,neighb1,neighb2,path)
!!    !insert at k_epot+1
!!    use module_base
!!    implicit none
!!    !parameters
!!    type(gt_data), intent(inout) :: gdat
!!    integer, intent(in) :: k_epot
!!    real(gp), intent(in) :: epot
!!    real(gp), intent(in) :: fp(gdat%nid)
!!    integer, intent(in) :: neighb1, neighb2
!!    character(len=600)   :: path
!!    !local
!!    integer :: i,k
!!    if(gdat%nsad+1>gdat%nsadtot)stop 'nsad+1>=nsadtot, out of bounds'
!!
!!    gdat%nsad=gdat%nsad+1
!!    do k=gdat%nsad-1,k_epot+1,-1
!!        gdat%en_arr_sad(k+1)=gdat%en_arr_sad(k)
!!        gdat%sadnumber(k+1)=gdat%sadnumber(k)
!!        gdat%path_sad(k+1)=gdat%path_sad(k)
!!        gdat%nneighbpairs(k+1) = gdat%nneighbpairs(k)
!!        gdat%paircounter(:,k+1) = gdat%paircounter(:,k)
!!        gdat%sadneighb(1,:,k+1)=gdat%sadneighb(1,:,k)
!!        gdat%sadneighb(2,:,k+1)=gdat%sadneighb(2,:,k)
!!        do i=1,gdat%nid
!!            gdat%fp_arr_sad(i,k+1)=gdat%fp_arr_sad(i,k)
!!         enddo
!!    enddo
!!    gdat%en_arr_sad(k_epot+1)=epot
!!    gdat%sadnumber(k_epot+1)=gdat%nsad
!!    gdat%path_sad(k_epot+1)=path
!!    gdat%nneighbpairs(k_epot+1) = 1
!!    gdat%paircounter(1,k_epot+1) = 1
!!    gdat%sadneighb(1,1,k_epot+1)=neighb1
!!    gdat%sadneighb(2,1,k_epot+1)=neighb2
!!    do i=1,gdat%nid
!!        gdat%fp_arr_sad(i,k+1)=fp(i)
!!    enddo
!!end subroutine insert_sad
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
  if (n.le.0) stop 'hunt_gt'
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
END SUBROUTINE hunt_gt
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
