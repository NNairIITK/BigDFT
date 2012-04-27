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
 
! Cut out a piece of the quantity (rho) from the global region (rho) and
! store it in a local region (Lrho).
 indSmall=0
 indSpin=0
 do ispin=1,nspin
     ! WARNING: I added the factors 2.
     do ii3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
         i3 = ii3
         do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
             i2 = ii2
             do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                 i1 = ii1
                 ! indSmall is the index in the local localization region
                 indSmall=indSmall+1
                 if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                                       !This initializes the buffers of locreg to zeros if outside the simulation box.
                     i3 <= Glr%d%n3i+1 .and. i2 <= Glr%d%n2i+1 .and. i1 <= Glr%d%n1i+1) then       !Should use periodic image instead... MUST FIX THIS.
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
  integer :: igrid,isegloc,isegG,ix,iorbs
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
  integer :: lfinc,Gfinc,spinshift,ispin,Gindex

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

subroutine local_overlap_matrices(norbe,norb1,norb2,nvctrp,nspin,nspinor,ndim_hamovr,hamovr,psi,psi2,hpsi,&
                                  iorbst,iorbst2,imatrst)
  use module_base
  implicit none
  integer, intent(in) :: norbe          ! total number of orbitals for overlap region
  integer, intent(in) :: norb1          ! number of orbitals in first locreg
  integer, intent(in) :: norb2          ! number of orbitals in second locreg
  integer, intent(in) :: nvctrp,ndim_hamovr,nspin,nspinor
  integer, intent(in) :: iorbst,iorbst2
  integer, intent(inout) :: imatrst
  real(wp), dimension(nspin*ndim_hamovr,2), intent(out) :: hamovr
  real(wp), dimension(nvctrp*nspinor,norbe), intent(in) :: psi,psi2,hpsi
  !local variables
  integer :: ispin,ncomp,ncplx
  !WARNING: here nspin=1 for nspinor=4
  if(nspinor == 1) then
     ncplx=1
  elseif(nspinor == 2) then
     ncplx=2
     ncomp=1
  else if (nspinor == 4) then
     ncplx=2
     ncomp=2
  end if

     if (nspinor ==1) then
        call gemm('T','N',norb1,norb2,nvctrp,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
             hpsi(1,iorbst2),max(1,nvctrp),&
             0.0_wp,hamovr(imatrst,1),norb1)
        !here probably dsyrk can be used
        call gemm('T','N',norb1,norb2,nvctrp,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
             psi2(1,iorbst2),max(1,nvctrp),0.0_wp,hamovr(imatrst,2),norb1)
     else
        call c_gemm('C','N',norb1,norb2,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(1,iorbst),&
             max(1,ncomp*nvctrp),hpsi(1,iorbst2),max(1,ncomp*nvctrp),&
             (0.0_wp,0.0_wp),hamovr(imatrst,1),norb1)
        !here probably zherk can be used
        call c_gemm('C','N',norb2,norb1,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(1,iorbst),&
             max(1,ncomp*nvctrp),psi2(1,iorbst2),max(1,ncomp*nvctrp),&
             (0.0_wp,0.0_wp),hamovr(imatrst,2),norb1)
     end if

     imatrst =imatrst+ncplx*norb1*norb2

END SUBROUTINE local_overlap_matrices

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
 
     ! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)... DONE
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
           !!! loop over the orbitals
           !!do iorbs=1,orbs%norbp*orbs%nspinor
           !!   lpsi(icheck+lincrement*(iorbs-1))=psi(Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1))
           !!end do
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
 
     ! Could optimize the routine by looping only on Gsegs not looped on before (TO DO).. DONE
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
              !lpsi(start+(icheck-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyv(isegG)+ix-1)*7+offset+igrid)
              !!do iorbs=1,orbs%norbp*orbs%nspinor
              !!   lpsi(start+icheck+lincrement*(iorbs-1)+igrid*lfinc)=&
              !!   psi(Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc)
              !!end do
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
  integer :: igrid,isegloc,isegG,ix,iorbs
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

     ! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)... DONE
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
           !!do iorbs=1,norb*nspinor
           !!   do ispin=1,nspin
           !!      Gindex = Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+shift+spinshift*(ispin-1)
           !!      Lindex = icheck+lincrement*(iorbs-1)+lincrement*norb*(ispin-1)
           !!      psi(Gindex) = psi(Gindex) + lpsi(Lindex)
           !!   end do
           !!end do
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

! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
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
           !!do igrid=0,6
           !!   do iorbs=1,norb*nspinor
           !!     do ispin = 1, nspin
           !!        Gindex = Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc+&
           !!                 shift + spinshift*(ispin-1)
           !!        Lindex = start+icheck+lincrement*(iorbs-1)+igrid*lfinc + lincrement*norb*(ispin-1) 
           !!        psi(Gindex) = psi(Gindex) + lpsi(Lindex)
           !!     end do
           !!   end do
           !!end do
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
     ! WARNING: I added the factors 2.
     !do i3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
     do ii3=i3s,i3e
         i3 = ii3
         !if(Glr%geocode /='F' .and. i3<0) i3 = i3 + Glr%d%n3i
         !if(Glr%geocode /='F' .and. i3>Glr%d%n3i+1) i3 = modulo(i3, Glr%d%n3i+1)
         do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
             i2 = ii2
             !if(Glr%geocode =='P' .and. i2<0) i2 = i2 + Glr%d%n2i
             !if(Glr%geocode /='F' .and. i2>Glr%d%n2i+1) i2 = modulo(i2, Glr%d%n2i+1)
             do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                 i1 = ii1
                 !if(Glr%geocode /='F' .and. i1<0) i1 = i1 + Glr%d%n1i
                 !if(Glr%geocode /='F' .and. i1>Glr%d%n1i+1) i1 = modulo(i1, Glr%d%n1i+1)
                 ! indSmall is the index in the local localization region
                 indSmall=indSmall+1
                 if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                                  !This initializes the buffers of locreg to zeros if outside the simulation box.
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


!> "Inserts" a quantity which is stored in the localized region into the glocal region.
!!        
!! @warning
!!    The quantity must not be stored in a compressed form. The output (rho) must be initialized to zero
!!    before entering this subroutine.
subroutine local_to_global_parallel(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho,i3s,i3e)

 use module_base
 use module_types
 
 implicit none

! Arguments
 type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
 type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
 integer, intent(in) :: size_rho  ! size of rho array
 integer, intent(in) :: size_Lrho ! size of Lrho array
 integer, intent(in) :: nspin  !number of spins
 real(wp),dimension(size_rho),intent(out) :: rho  ! quantity in global region
 real(wp),dimension(size_Lrho),intent(in) :: Lrho ! piece of quantity in local region
 integer,intent(in):: i3s, i3e ! starting and ending indices on z direction (related to distribution of rho when parallel)

! Local variable
 integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
 integer :: indSmall, indSpin, indLarge ! indexes for the arrays
 
! Cut out a piece of the quantity (rho) from the global region (rho) and
! store it in a local region (Lrho).
 indSmall=0
 indSpin=0
 do ispin=1,nspin
     ! WARNING: I added the factors 2.
     !do i3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
     do ii3=i3s,i3e
         i3 = ii3
         !if(Glr%geocode /='F' .and. i3<0) i3 = i3 + Glr%d%n3i
         !if(Glr%geocode /='F' .and. i3>Glr%d%n3i+1) i3 = modulo(i3, Glr%d%n3i+1)
         do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
             i2=ii2
             !if(Glr%geocode =='P' .and. i2<0) i2 = i2 + Glr%d%n2i
             !if(Glr%geocode /='F' .and. i2>Glr%d%n2i+1) i2 = modulo(i2, Glr%d%n2i+1)
             do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                 i1=ii1
                ! if(Glr%geocode /='F' .and. i1<0) i1 = i1 + Glr%d%n1i
                ! if(Glr%geocode /='F' .and. i1>Glr%d%n1i+1) i1 = modulo(i1, Glr%d%n1i+1)
                 ! indSmall is the index in the local localization region
                 indSmall=indSmall+1
                 if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                                  !This initializes the buffers of locreg to zeros if outside the simulation box.
                     i3 <= Glr%d%n3i .and. i2 <= Glr%d%n2i .and. i1 <= Glr%d%n1i) then           !Should use periodic image instead... MUST FIX THIS.
                    ! indLarge is the index in the global localization region. 
                    indLarge=(i3-1)*Glr%d%n2i*Glr%d%n1i + (i2-1)*Glr%d%n1i + i1
                    !Lrho(indSmall)=rho(indLarge+indSpin)
                    rho(indLarge+indSpin)=Lrho(indSmall)
                 else
                    rho(indLarge+indSpin)= 0.0_wp
                 end if
             end do
         end do
     end do
     indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
 end do

END SUBROUTINE local_to_global_parallel



