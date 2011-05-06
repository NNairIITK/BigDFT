!#############################################################################################################################################
!!****f* BigDFT/psi_to_locreg
!#############################################################################################################################################
!! FUNCTION: Tranform wavefunction between Global region and localisation region
!!
!! WARNING: 
!!         Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
!! SOURCE:
!!
subroutine psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,nlr,orbs,psi)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: nlr                  ! number of localization regions
  integer :: ilr           ! index of the localization region we are considering
  integer :: ldim          ! dimension of lpsi 
  type(orbitals_data),intent(in) :: orbs      ! orbital descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Olr  ! Localization grid descriptors 
  real(wp),dimension(orbs%npsidim),intent(in) :: psi       !Wavefunction (compressed format)
  real(wp),dimension(ldim),intent(inout) :: lpsi !Wavefunction in localization region
  !#############################################
  !local variables
  !############################################
  integer :: igrid,isegloc,isegG,ix,iorbs
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
  integer :: lfinc,Gfinc

! Define integers
  nseg = Olr(ilr)%wfd%nseg_c + Olr(ilr)%wfd%nseg_f
  lincrement = Olr(ilr)%wfd%nvctr_c + 7*Olr(ilr)%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0

! Initialize loc_psi
  call razero(lincrement*orbs%norbp*orbs%nspinor,lpsi)
 
! Get the keymask: shift for every segment of Llr (with respect to Glr)
  allocate(keymask(2,nseg),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  call shift_locreg_indexes(Glr,Olr(ilr),keymask,nseg)

!####################################################
! Do coarse region
!####################################################
  do isegloc = 1,Olr(ilr)%wfd%nseg_c
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
 
! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
     do isegG = 1,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keyg(1,isegG)
        Gmax = Glr%wfd%keyg(2,isegG)

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
 
        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           ! loop over the orbitals
           do iorbs=1,orbs%norbp*orbs%nspinor
              lpsi(icheck+lincrement*(iorbs-1))=psi(Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1))
           end do
        end do
     end do
  end do

! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Olr(ilr)%wfd%nvctr_c) then
    write(*,*)'There is an error in psi_to_locreg: number of coarse points used',icheck
    write(*,*)'is not equal to the number of coarse points in the region',Olr(ilr)%wfd%nvctr_c
  end if

!##############################################################
! Now do fine region
!##############################################################

  icheck = 0
  start = Olr(ilr)%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Olr(ilr)%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  do isegloc = Olr(ilr)%wfd%nseg_c+1,nseg
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
 
! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
     do isegG = Glr%wfd%nseg_c+1,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keyg(1,isegG)
        Gmax = Glr%wfd%keyg(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           icheck = icheck + 1
           do igrid=0,6
              do iorbs=1,orbs%norbp*orbs%nspinor
                 lpsi(start+icheck+lincrement*(iorbs-1)+igrid*lfinc)=&
&                psi(Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc)
              end do
           end do
        end do
     end do
  end do
  
 ! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Olr(ilr)%wfd%nvctr_f) then
    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Olr(ilr)%wfd%nvctr_f
  end if

  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE psi_to_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/shift_locreg_indexes
!#############################################################################################################################################
!! FUNCTION:
!!        Find the shift necessary for the indexes of every segment of Blr
!!        to make them compatible with the indexes of Alr. These shifts are
!!        returned in the array keymask(nseg), where nseg should be the number
!!        of segments in Blr.
!! WARNING: 
!!         This routine supposes that the region Blr is contained in the region Alr.
!!         This should always be the case, if we concentrate on the overlap between two regions.
!! SOURCE:
!!
subroutine shift_locreg_indexes(Alr,Blr,keymask,nseg)

  use module_base
  use module_types
 
 implicit none

!############################
! Arguments
!############################
 type(locreg_descriptors),intent(in) :: Alr,Blr   ! The two localization regions
 integer,intent(in) :: nseg
 integer,intent(out) :: keymask(2,nseg)
!#############################
! Local variable
!#############################
 integer :: iseg      !integer for the loop
 integer :: Bindex    !starting index of segments in Blr
 integer :: x,y,z     !coordinates of start of segments in Blr 
 integer :: shift(3)  !shift between the beginning of the segment in Blr and the origin of Alr
 integer ::  tmp

!Big loop on all segments
 do iseg=1,nseg

!##########################################
! For the Starting index
    Bindex = Blr%wfd%keyg(1,iseg)
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

    Bindex = Blr%wfd%keyg(2,iseg)
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
!%***

!#############################################################################################################################################
!!****f* BigDFT/global_to_local
!#############################################################################################################################################
!! FUNCTION: Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
!!        
!! WARNING: The quantity must not be stored in a compressed form.
!!
!! SOURCE:
!!
subroutine global_to_local(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho)

 use module_base
 use module_types
 
 implicit none

!############################
! Arguments
!############################
 type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
 type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
 integer, intent(in) :: size_rho  ! size of rho array
 integer, intent(in) :: size_Lrho ! size of Lrho array
 integer, intent(in) :: nspin  !number of spins
 real(wp),dimension(size_rho),intent(in) :: rho  ! quantity in global region
 real(wp),dimension(size_Lrho),intent(out) :: Lrho ! piece of quantity in local region
!#############################
! Local variable
!#############################
 integer :: ispin,i1,i2,i3  !integer for loops
 integer :: indSmall, indSpin, indLarge ! indexes for the arrays
 
! Cut out a piece of the quantity (rho) from the global region (rho) and
! store it in a local region (Lrho).
 indSmall=0
 indSpin=0
 do ispin=1,nspin
     do i3=Llr%ns3+1,Llr%d%n3i+Llr%ns3
         do i2=Llr%ns2+1,Llr%d%n2i+Llr%ns2
             do i1=Llr%ns1+1,Llr%d%n1i+Llr%ns1
                 ! indSmall is the index in the local localization region
                 indSmall=indSmall+1
                 ! indLarge is the index in the global localization region. 
                 indLarge=(i3-1)*Glr%d%n2i*Glr%d%n1i + (i2-1)*Glr%d%n1i + i1
                 Lrho(indSmall)=rho(indLarge+indSpin)
             end do
         end do
     end do
     indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
 end do

END SUBROUTINE global_to_local


!#############################################################################################################################################
!!****f* BigDFT/Lpsi_to_global
!#############################################################################################################################################
!! FUNCTION: Tranform wavefunction between localisation region and the global region
!!
!! WARNING: 
!!         Psi must be initialized to zero before entering this routine. Each Lpsi is added to the corresponding place in Global.
!!         Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
!! SOURCE:
!!
subroutine Lpsi_to_global(Glr,Gdim,Llr,lpsi,orbs,psi)

  use module_base
  use module_types

 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer :: Gdim          ! dimension of lpsi 
  type(orbitals_data),intent(in) :: orbs      ! orbital descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  real(wp),dimension(Gdim),intent(inout) :: psi       !Wavefunction (compressed format)
  real(wp),dimension(orbs%npsidim),intent(in) :: lpsi !Wavefunction in localization region
  !#############################################
  !local variables
  !############################################
  integer :: igrid,isegloc,isegG,ix,iorbs
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
  integer :: lfinc,Gfinc

! Define integers
  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0

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
        Gmin = Glr%wfd%keyg(1,isegG)
        Gmax = Glr%wfd%keyg(2,isegG)

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
           do iorbs=1,orbs%norbp*orbs%nspinor
              psi(Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1))=&
              psi(Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)) + lpsi(icheck+lincrement*(iorbs-1))
           end do
        end do
     end do
  end do

! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_c) then
    write(*,*)'There is an error in psi_to_locreg: number of coarse points used',icheck
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

        Gmin = Glr%wfd%keyg(1,isegG)
        Gmax = Glr%wfd%keyg(2,isegG)

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
              do iorbs=1,orbs%norbp*orbs%nspinor
                psi(Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc) = &
                psi(Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc) +&
                lpsi(start+icheck+lincrement*(iorbs-1)+igrid*lfinc)
              end do
           end do
        end do
     end do
  end do

 ! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_f) then
    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
    write(*,*)'is not equal to the number of fine points in the region',Llr%wfd%nvctr_f
  end if

  deallocate(keymask,stat=i_stat)
  call memocc(i_stat,i_all,'keymask',subname)

END SUBROUTINE Lpsi_to_global
!%***

