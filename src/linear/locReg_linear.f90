!> @file
!!  Routines used by the linear scaling version
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS





!> Determine a set of localisation regions from the centers and the radii.
!! cut in cubes the global reference system
subroutine check_linear_inputguess(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,linear)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: nlr
  logical,intent(out) :: linear
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  !local variables
  logical :: warningx,warningy,warningz
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3
  real(gp) :: rx,ry,rz,cutoff
  
  linear = .true.

  !determine the limits of the different localisation regions
  do ilr=1,nlr

     !initialize logicals
     warningx = .false.
     warningy = .false.
     warningz = .false.

     rx=cxyz(1,ilr)
     ry=cxyz(2,ilr)
     rz=cxyz(3,ilr)

     cutoff=locrad(ilr)

     isx=floor((rx-cutoff)/hx)
     isy=floor((ry-cutoff)/hy)
     isz=floor((rz-cutoff)/hz)

     iex=ceiling((rx+cutoff)/hx)
     iey=ceiling((ry+cutoff)/hy)
     iez=ceiling((rz+cutoff)/hz)

     ln1 = iex-isx
     ln2 = iey-isy
     ln3 = iez-isz

     ! First check if localization region fits inside box
        if (iex - isx >= Glr%d%n1 - 14) then
           warningx = .true.
        end if
        if (iey - isy >= Glr%d%n2 - 14) then
           warningy = .true.
        end if
        if (iez - isz >= Glr%d%n3 - 14) then
           warningz = .true.
        end if 

     !If not, then don't use linear input guess (set linear to false)
     if(warningx .and. warningy .and. warningz .and. (Glr%geocode .ne. 'F')) then
       linear = .false.
       if(iproc == 0) then
          write(*,*)'Not using the linear scaling input guess, because localization'
          write(*,*)'region greater or equal to simulation box.'
       end if
       exit 
     end if
  end do
      
end subroutine check_linear_inputguess





!> Tranform wavefunction between localisation region and the global region



!!!> Tranform wavefunction between localisation region and the global region.
!!!! The global region can have any boundary conditions.
!!!! @warning 
!!!! WARNING: Make sure psi is set to zero where Glr does not collide with Llr (or everywhere)
!!subroutine Lpsi_to_global_free_to_any(iproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, lpsi, psi)
!!
!!  use module_base
!!  use module_types
!!
!! implicit none
!!
!!  ! Subroutine Scalar Arguments
!!  integer,intent(in):: iproc
!!  integer :: Gdim          ! dimension of psi 
!!  integer :: Ldim          ! dimension of lpsi
!!  integer :: norb          ! number of orbitals
!!  integer :: nspinor       ! number of spinors
!!  integer :: nspin         ! number of spins 
!!  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!  type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
!!  
!!  !Subroutine Array Arguments
!!  real(wp),dimension(Gdim),intent(inout) :: psi       !Wavefunction (compressed format)
!!  real(wp),dimension(Ldim),intent(in) :: lpsi         !Wavefunction in localization region
!!  
!!  !local variables
!!  integer :: igrid,isegloc,isegG,ix!,iorbs
!!  integer :: lmin,lmax,Gmin,Gmax
!!  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
!!  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
!!  integer :: length      ! Length of the overlap between Lseg and Gseg
!!  integer :: lincrement  ! Increment for writing orbitals in loc_psi
!!  integer :: Gincrement  ! Increment for reading orbitals in psi
!!  integer :: nseg        ! total number of segments in Llr
!!  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
!!  character(len=*), parameter :: subname='Lpsi_to_global'
!!  integer :: i_all
!!  integer :: start,Gstart,Lindex
!!  integer :: lfinc,Gfinc,spinshift,ispin,Gindex,isegstart
!!  integer :: istart
!!  !integer :: i_stat
!!
!!  call f_routine(id=subname)
!!
!!  ! This routine is only intended for conversions between locregs with the same boundary conditions.
!!  if (glr%geocode/= 'F' .or. llr%geocode/='F') then
!!      call f_err_throw('Lpsi_to_global2 can only be used for locregs with free boundary conditions', &
!!           err_name='BIGDFT_RUNTIME_ERROR')
!!  end if
!!
!!  if(nspin/=1) stop 'not fully implemented for nspin/=1!'
!!
!!! Define integers
!!  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
!!  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
!!  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
!!  icheck = 0
!!  spinshift = Gdim / nspin
!! 
!!! Get the keymask: shift for every segment of Llr (with respect to Glr)
!!! allocate(keymask(2,nseg),stat=i_stat)
!!  keymask = f_malloc((/2,nseg/),id='keymask')
!!
!!  call shift_locreg_indexes(Glr,Llr,keymask,nseg)
!!
!!!####################################################
!!! Do coarse region
!!!####################################################
!!  isegstart=1
!!
!! 
!!  !$omp parallel default(private) &
!!  !$omp shared(Glr,Llr, keymask,lpsi,icheck,psi,norb) &
!!  !$omp firstprivate(isegstart,nseg,lincrement,Gincrement,spinshift,nspin) 
!!
!!  !$omp do reduction(+:icheck)
!!  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
!!     lmin = keymask(1,isegloc)
!!     lmax = keymask(2,isegloc)
!!     istart = llr%wfd%keyvloc(isegloc)-1
!!
!!     
!!     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
!!        Gmin = Glr%wfd%keygloc(1,isegG)
!!        Gmax = Glr%wfd%keygloc(2,isegG)
!!
!!        ! For each segment in Llr check if there is a collision with the segment in Glr
!!        !if not, cycle
!!        if(lmin > Gmax) then
!!            isegstart=isegG
!!        end if
!!        if(Gmin > lmax) exit global_loop_c
!!
!!        !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_c
!!        if(lmin > Gmax)  cycle global_loop_c
!!
!!        ! Define the offset between the two segments
!!        offset = lmin - Gmin
!!        if(offset < 0) then
!!           offset = 0
!!        end if
!!
!!        ! Define the length of the two segments
!!        length = min(lmax,Gmax)-max(lmin,Gmin)
!!
!!        !Find the common elements and write them to the new global wavefunction
!!        icheck = icheck + (length + 1)
!!
!!        ! WARNING: index goes from 0 to length because it is the offset of the element
!!
!!        do ix = 0,length     
!!           istart = istart + 1
!!           do ispin=1,nspin
!!              Gindex = Glr%wfd%keyvloc(isegG)+offset+ix+spinshift*(ispin-1)
!!              Lindex = istart+lincrement*norb*(ispin-1)
!!              psi(Gindex) = lpsi(Lindex) 
!!           end do
!!        end do
!!     end do global_loop_c
!!  end do local_loop_c
!!  !$omp end do
!!
!!
!!!##############################################################
!!! Now do fine region
!!!##############################################################
!!
!!  start = Llr%wfd%nvctr_c
!!  Gstart = Glr%wfd%nvctr_c
!!  lfinc  = Llr%wfd%nvctr_f
!!  Gfinc = Glr%wfd%nvctr_f
!!
!!  isegstart=Glr%wfd%nseg_c+1
!!
!!  !$omp do reduction(+:icheck)
!!  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
!!     lmin = keymask(1,isegloc)
!!     lmax = keymask(2,isegloc)
!!     istart = llr%wfd%keyvloc(isegloc)-1
!!
!!     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
!!
!!        Gmin = Glr%wfd%keygloc(1,isegG)
!!        Gmax = Glr%wfd%keygloc(2,isegG)
!!
!!        ! For each segment in Llr check if there is a collision with the segment in Glr
!!        ! if not, cycle
!!        if(lmin > Gmax) then
!!            isegstart=isegG
!!        end if
!!        if(Gmin > lmax)  exit global_loop_f
!!        !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f
!!        if(lmin > Gmax)  cycle global_loop_f
!!
!!        offset = lmin - Gmin
!!        if(offset < 0) offset = 0
!!
!!        length = min(lmax,Gmax)-max(lmin,Gmin)
!!
!!        !Find the common elements and write them to the new global wavefunction
!!        ! First set to zero those elements which are not copied. WARNING: will not work for npsin>1!!
!! 
!!        icheck = icheck + (length + 1)
!!
!!        ! WARNING: index goes from 0 to length because it is the offset of the element
!!        do ix = 0,length
!!        istart = istart + 1
!!           do igrid=1,7
!!              do ispin = 1, nspin
!!                 Gindex = Gstart + (Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid + spinshift*(ispin-1)
!!                 Lindex = start+(istart-1)*7+igrid + lincrement*norb*(ispin-1) 
!!                 psi(Gindex) = lpsi(Lindex) 
!!              end do
!!           end do
!!        end do
!!     end do global_loop_f
!!  end do local_loop_f
!!  !$omp end do
!!
!!  !$omp end parallel
!!
!!  !Check if the number of elements in loc_psi is valid
!!  if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
!!    write(*,*)'There is an error in Lpsi_to_global: sum of fine and coarse points used',icheck
!!    write(*,*)'is not equal to the sum of fine and coarse points in the region',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
!!    stop
!!  end if
!!
!!  call f_free(keymask)
!!
!!  call f_release_routine()
!!
!!END SUBROUTINE Lpsi_to_global_free_to_any
