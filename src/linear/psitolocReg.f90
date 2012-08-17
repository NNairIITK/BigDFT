!> @file
!! Wavefunction put into a localisation region
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Find the shift necessary for the indexes of every segment of Blr
!!   to make them compatible with the indexes of Alr. These shifts are
!!   returned in the array keymask(nseg), where nseg should be the number
!!   of segments in Blr.
!! @warning 
!!   This routine supposes that the region Blr is contained in the region Alr.
!!   This should always be the case, if we concentrate on the overlap between two regions.
subroutine shift_locreg_indexes(Alr,Blr,keymask,nseg)

  use module_base
  use module_types
 
 implicit none

! Arguments
 type(locreg_descriptors),intent(in) :: Alr,Blr   ! The two localization regions
 integer,intent(in) :: nseg
 integer,intent(out) :: keymask(2,nseg)

! Local variable
 integer :: iseg      !integer for the loop
 integer :: Bindex    !starting index of segments in Blr
 integer :: x,y,z     !coordinates of start of segments in Blr 
 integer :: shift(3)  !shift between the beginning of the segment in Blr and the origin of Alr
 integer ::  tmp

!Big loop on all segments
 do iseg=1,nseg

!##########################################
! For the Starting index
    Bindex = Blr%wfd%keygloc(1,iseg)
    tmp = Bindex -1
    z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
    tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
    y   = tmp / (Blr%d%n1+1)
    x   = tmp - y * (Blr%d%n1+1)
 
! Shift between the beginning of the segment and the start of the Alr region
    shift(1) = x + Blr%ns1 - Alr%ns1
    shift(2) = y + Blr%ns2 - Alr%ns2
    shift(3) = z + Blr%ns3 - Alr%ns3

! Write the shift in index form
    keymask(1,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1

!######################################
! For the ending index

    Bindex = Blr%wfd%keygloc(2,iseg)
    tmp = Bindex -1
    z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
    tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
    y   = tmp / (Blr%d%n1+1)
    x   = tmp - y * (Blr%d%n1+1)

! Shift between the beginning of the segment and the start of the Alr region
    shift(1) = x + Blr%ns1 - Alr%ns1
    shift(2) = y + Blr%ns2 - Alr%ns2
    shift(3) = z + Blr%ns3 - Alr%ns3

! Write the shift in index form
    keymask(2,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1
 end do

END SUBROUTINE shift_locreg_indexes


!> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
!! @warning: The quantity must not be stored in a compressed form.
subroutine global_to_local(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho)

 use module_base
 use module_types
 
 implicit none

! Arguments
 type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
 type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
 integer, intent(in) :: size_rho  ! size of rho array
 integer, intent(in) :: size_Lrho ! size of Lrho array
 integer, intent(in) :: nspin  !number of spins
 real(wp),dimension(size_rho),intent(in) :: rho  ! quantity in global region
 real(wp),dimension(size_Lrho),intent(out) :: Lrho ! piece of quantity in local region

! Local variable
 integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
 integer :: indSmall, indSpin, indLarge ! indexes for the arrays
 logical:: z_inside, y_inside, x_inside
 integer:: iz, iy, m
 
! Cut out a piece of the quantity (rho) from the global region (rho) and
! store it in a local region (Lrho).

 if(Glr%geocode == 'F') then
     ! Use loop unrolling here
     indSmall=0
     indSpin=0
     do ispin=1,nspin
         ! WARNING: I added the factors 2.
         do i3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
             iz=(i3-1)*Glr%d%n2i*Glr%d%n1i
             do i2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
                 iy=(i2-1)*Glr%d%n1i
                 m=mod(Llr%d%n1i+Llr%nsi1-Llr%nsi1,4)
                 if(m/=0) then
                     do i1=Llr%nsi1+1,Llr%nsi1+m
                        indSmall=indSmall+1
                        indLarge=iz+iy+i1
                        Lrho(indSmall)=rho(indLarge+indSpin)
                     end do
                  end if
                  do i1=Llr%nsi1+1+m,Llr%d%n1i+Llr%nsi1,4
                     Lrho(indSmall+1)=rho(iz+iy+i1+0+indSpin)
                     Lrho(indSmall+2)=rho(iz+iy+i1+1+indSpin)
                     Lrho(indSmall+3)=rho(iz+iy+i1+2+indSpin)
                     Lrho(indSmall+4)=rho(iz+iy+i1+3+indSpin)
                     indSmall=indSmall+4
                  end do
             end do
         end do
         indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
     end do
 else
     ! General case
     indSmall=0
     indSpin=0
     do ispin=1,nspin
         ! WARNING: I added the factors 2.
         do ii3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
             i3 = mod(ii3-1,Glr%d%n3i)+1
             z_inside = (i3>0 .and. i3<=Glr%d%n3i+1)
             iz=(i3-1)*Glr%d%n2i*Glr%d%n1i
             do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
                 i2 = mod(ii2-1,Glr%d%n2i)+1
                 y_inside = (i2>0 .and. i2<=Glr%d%n2i+1)
                 iy=(i2-1)*Glr%d%n1i
                 do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                     i1 = mod(ii1-1,Glr%d%n1i)+1 
                     x_inside = (i1 > 0 .and. i1 <= Glr%d%n1i+1)
                     ! indSmall is the index in the local localization region
                     indSmall=indSmall+1
                     !!if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                                       !This initializes the buffers of locreg to zeros if outside the simulation box.
                     !!    i3 <= Glr%d%n3i+1 .and. i2 <= Glr%d%n2i+1 .and. i1 <= Glr%d%n1i+1) then       !Should use periodic image instead... MUST FIX THIS.
                     !!   ! indLarge is the index in the global localization region. 
                     !!   indLarge=(i3-1)*Glr%d%n2i*Glr%d%n1i + (i2-1)*Glr%d%n1i + i1
                     if(z_inside .and. y_inside .and. x_inside) then
                        indLarge=iz+iy+i1
                        Lrho(indSmall)=rho(indLarge+indSpin)
                     else
                        Lrho(indSmall)= 0.0_wp
                     end if
                 end do
             end do
         end do
         indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
     end do
 end if

END SUBROUTINE global_to_local


!> Tranform wavefunction between localisation region and the global region
!! @warning 
!!     Psi must be initialized to zero before entering this routine. Each Lpsi is added to the corresponding place in Global.
!!     Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
subroutine Lpsi_to_global(Glr,Gdim,Llr,lpsi,Ldim,norb,nspinor,nspin,shift,psi)

  use module_base
  use module_types

 implicit none

  ! Subroutine Scalar Arguments
  integer :: Gdim          ! dimension of psi 
  integer :: Ldim          ! dimension of lpsi
  integer :: norb          ! number of orbitals
  integer :: nspinor       ! number of spinors
  integer :: nspin         ! number of spins 
  integer :: shift         ! shift to correct place in locreg
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
  
  !Subroutine Array Arguments
  real(wp),dimension(Gdim),intent(inout) :: psi       !Wavefunction (compressed format)
  real(wp),dimension(Ldim),intent(in) :: lpsi         !Wavefunction in localization region
  
  !local variables
  integer :: igrid,isegloc,isegG,ix!,iorbs
  integer :: lmin,lmax,Gmin,Gmax
  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
  integer :: length      ! Length of the overlap between Lseg and Gseg
  integer :: lincrement  ! Increment for writing orbitals in loc_psi
  integer :: Gincrement  ! Increment for reading orbitals in psi
  integer :: nseg        ! total number of segments in Llr
  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
  character(len=*), parameter :: subname='Lpsi_to_global'
  integer :: i_stat,i_all
  integer :: start,Gstart,Lindex
  integer :: lfinc,Gfinc,spinshift,Gindex!,ispin

! Define integers
  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0
  spinshift = Gdim / nspin                  !MUST CHANGE THIS
 
! Get the keymask: shift for every segment of Llr (with respect to Glr)
  allocate(keymask(2,nseg),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  call shift_locreg_indexes(Glr,Llr,keymask,nseg)

!####################################################
! Do coarse region
!####################################################
  do isegloc = 1,Llr%wfd%nseg_c
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
     do isegG = 1,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)
        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle

        ! Define the offset between the two segments
        offset = lmin - Gmin
        if(offset < 0) then
           offset = 0
        end if

        ! Define the length of the two segments
        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new global wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           ! loop over the orbitals
!           do iorbs=1,norb*nspinor
!              do ispin=1,nspin
!                 Gindex = Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+shift+spinshift*(ispin-1)
!                 Lindex = icheck+lincrement*(iorbs-1)+lincrement*norb*(ispin-1)
                 Gindex = Glr%wfd%keyvloc(isegG)+offset+ix+shift!+spinshift*(ispin-1)
                 Lindex = icheck!+lincrement*norb*(ispin-1)
                 psi(Gindex) = psi(Gindex) + lpsi(Lindex)
!              end do
!           end do
        end do
     end do
  end do
! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_c) then
    write(*,*)'There is an error in Lpsi_to_global: number of coarse points used',icheck
    write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
  end if

!##############################################################
! Now do fine region
!##############################################################

  icheck = 0
  start = Llr%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Llr%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  do isegloc = Llr%wfd%nseg_c+1,nseg
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)

! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
     do isegG = Glr%wfd%nseg_c+1,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new global wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           do igrid=0,6
!              do iorbs=1,norb*nspinor
!                do ispin = 1, nspin
!                   Gindex = Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc+&
!                            shift + spinshift*(ispin-1)
!                   Lindex = start+icheck+lincrement*(iorbs-1)+igrid*lfinc + lincrement*norb*(ispin-1) 
                   Gindex=Gstart+Glr%wfd%keyvloc(isegG)+offset+ix+igrid*Gfinc+shift!+spinshift*(ispin-1)
                   Lindex=start+icheck+igrid*lfinc !+ lincrement*norb*(ispin-1)
                   psi(Gindex) = psi(Gindex) + lpsi(Lindex)
!                end do
!              end do
           end do
        end do
     end do
  end do
 ! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_f) then
    write(*,*)'There is an error in Lpsi_to_global: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Llr%wfd%nvctr_f
  end if

  i_all=-product(shape(keymask))*kind(keymask)
  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE Lpsi_to_global


!> Tranform one wavefunction between Global region and localisation region
subroutine psi_to_locreg2(iproc, nproc, ldim, gdim, Llr, Glr, gpsi, lpsi)

  use module_base
  use module_types
 
 implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: iproc                  ! process ID
  integer,intent(in) :: nproc                  ! number of processes
  integer,intent(in) :: ldim          ! dimension of lpsi 
  integer,intent(in) :: gdim          ! dimension of gpsi 
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  
  !Subroutine Array Arguments
  real(wp),dimension(gdim),intent(in) :: gpsi       !Wavefunction (compressed format)
  real(wp),dimension(ldim),intent(out) :: lpsi   !Wavefunction in localization region
  
  !local variables
  integer :: igrid,isegloc,isegG,ix!,iorbs
  integer :: lmin,lmax,Gmin,Gmax
  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
  integer :: length      ! Length of the overlap between Lseg and Gseg
  integer :: lincrement  ! Increment for writing orbitals in loc_psi
  integer :: Gincrement  ! Increment for reading orbitals in psi
  integer :: nseg        ! total number of segments in Llr
  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
  character(len=*), parameter :: subname='psi_to_locreg'
  integer :: i_stat,i_all
  integer :: start,Gstart
  integer :: lfinc,Gfinc,isegstart

! Define integers
  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0

! Initialize loc_psi
  call razero(ldim, lpsi)

! Get the keymask: shift for every segment of Llr (with respect to Glr)
  allocate(keymask(2,nseg),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  call shift_locreg_indexes(Glr,Llr,keymask,nseg)

!####################################################
! Do coarse region
!####################################################
  isegstart=1
  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
 
     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax) cycle local_loop_c
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c
        
        ! Define the offset between the two segments
        offset = lmin - Gmin
        if(offset < 0) then
           offset = 0
        end if
    
        ! Define the length of the two segments
        length = min(lmax,Gmax)-max(lmin,Gmin)
 
        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           lpsi(icheck) = gpsi(Glr%wfd%keyvloc(isegG)+offset+ix)
        end do
     end do global_loop_c
  end do local_loop_c

! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_c) then
    write(*,*)'There is an error in psi_to_locreg2: number of coarse points used',icheck
    write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
  end if

!##############################################################
! Now do fine region
!##############################################################

  icheck = 0
  start = Llr%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Llr%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  isegstart=Glr%wfd%nseg_c+1
  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
 
     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax) cycle local_loop_f
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_f

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           do igrid=1,7
              lpsi(start+(icheck-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid)
           end do
        end do
     end do global_loop_f
  end do local_loop_f
  
 ! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_f) then
    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Llr%wfd%nvctr_f
  end if

  i_all=-product(shape(keymask))*kind(keymask)
  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE psi_to_locreg2


!> Tranform wavefunction between localisation region and the global region
!! @warning 
!!    Psi must be initialized to zero before entering this routine. Each Lpsi is added to the corresponding place in Global.
!!    Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
!subroutine Lpsi_to_global2(Glr,Gdim,Llr,lpsi,Ldim,norb,nspinor,nspin,shift,psi)
subroutine Lpsi_to_global2(iproc, nproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, lpsi, psi)

  use module_base
  use module_types

 implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in):: iproc, nproc
  integer :: Gdim          ! dimension of psi 
  integer :: Ldim          ! dimension of lpsi
  integer :: norb          ! number of orbitals
  integer :: nspinor       ! number of spinors
  integer :: nspin         ! number of spins 
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
  
  !Subroutine Array Arguments
  real(wp),dimension(Gdim),intent(inout) :: psi       !Wavefunction (compressed format)
  real(wp),dimension(Ldim),intent(in) :: lpsi         !Wavefunction in localization region
  
  !local variables
  integer :: igrid,isegloc,isegG,ix!,iorbs
  integer :: lmin,lmax,Gmin,Gmax
  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
  integer :: length      ! Length of the overlap between Lseg and Gseg
  integer :: lincrement  ! Increment for writing orbitals in loc_psi
  integer :: Gincrement  ! Increment for reading orbitals in psi
  integer :: nseg        ! total number of segments in Llr
  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
  character(len=*), parameter :: subname='Lpsi_to_global'
  integer :: i_stat,i_all
  integer :: start,Gstart,Lindex
  integer :: lfinc,Gfinc,spinshift,ispin,Gindex,isegstart

! Define integers
  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0
  spinshift = Gdim / nspin
 
! Get the keymask: shift for every segment of Llr (with respect to Glr)
  allocate(keymask(2,nseg),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  call shift_locreg_indexes(Glr,Llr,keymask,nseg)

!####################################################
! Do coarse region
!####################################################
  isegstart=1
  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)

     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax) cycle local_loop_c
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c

        ! Define the offset between the two segments
        offset = lmin - Gmin
        if(offset < 0) then
           offset = 0
        end if

        ! Define the length of the two segments
        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new global wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           ! loop over the orbitals
           do ispin=1,nspin
              Gindex = Glr%wfd%keyvloc(isegG)+offset+ix+spinshift*(ispin-1)
              Lindex = icheck+lincrement*norb*(ispin-1)
              psi(Gindex) = psi(Gindex) + lpsi(Lindex)
           end do
        end do
     end do global_loop_c
  end do local_loop_c

! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_c) then
    write(*,*)'There is an error in Lpsi_to_global: number of coarse points used',icheck
    write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
  end if

!##############################################################
! Now do fine region
!##############################################################

  icheck = 0
  start = Llr%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Llr%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  isegstart=Glr%wfd%nseg_c+1
  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)

     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax) cycle local_loop_f
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_f

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new global wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           do igrid=1,7
              do ispin = 1, nspin
                 Gindex = Gstart + (Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid + spinshift*(ispin-1)
                 Lindex = start+(icheck-1)*7+igrid + lincrement*norb*(ispin-1) 
                 psi(Gindex) = psi(Gindex) + lpsi(Lindex)
              end do
           end do
        end do
     end do global_loop_f
  end do local_loop_f

 ! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_f) then
    write(*,*)'There is an error in Lpsi_to_global: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Llr%wfd%nvctr_f
  end if

  i_all=-product(shape(keymask))*kind(keymask)
  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE Lpsi_to_global2


!> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
!! @warning       
!!    The quantity must not be stored in a compressed form.
subroutine global_to_local_parallel(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho,i3s,i3e)

 use module_base
 use module_types
 
 implicit none

! Arguments
 type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
 type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
 integer, intent(in) :: size_rho  ! size of rho array
 integer, intent(in) :: size_Lrho ! size of Lrho array
 integer, intent(in) :: nspin  !number of spins
 real(wp),dimension(size_rho),intent(in) :: rho  ! quantity in global region
 real(wp),dimension(size_Lrho),intent(out) :: Lrho ! piece of quantity in local region
 integer,intent(in):: i3s, i3e ! starting and ending indices on z direction (related to distribution of rho when parallel)

! Local variable
 integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
 integer :: indSmall, indSpin, indLarge ! indexes for the arrays
 
! Cut out a piece of the quantity (rho) from the global region (rho) and
! store it in a local region (Lrho).
 indSmall=0
 indSpin=0
 do ispin=1,nspin
     do ii3=i3s,i3e
         i3 = mod(ii3-1,Glr%d%n3i)+1
         do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
             i2 = mod(ii2-1,Glr%d%n2i)+1
             do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                 i1 = mod(ii1-1,Glr%d%n1i)+1
                 ! indSmall is the index in the local localization region
                 indSmall=indSmall+1
                 if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                    !DON'T NEED THIS IF ANYMORE. This initializes the buffers of locreg to zeros if outside the simulation box.
                     i3 <= Glr%d%n3i .and. i2 <= Glr%d%n2i .and. i1 <= Glr%d%n1i) then           !Should use periodic image instead... MUST FIX THIS.
                    ! indLarge is the index in the global localization region. 
                    indLarge=(i3-1)*Glr%d%n2i*Glr%d%n1i + (i2-1)*Glr%d%n1i + i1
                    Lrho(indSmall)=rho(indLarge+indSpin)
                 else
                    Lrho(indSmall)= 0.0_wp
                 end if
             end do
         end do
     end do
     indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
 end do

END SUBROUTINE global_to_local_parallel
