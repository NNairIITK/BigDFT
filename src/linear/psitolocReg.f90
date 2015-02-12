!> @file
!! Wavefunction put into a localisation region
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
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
!$omp parallel do default(private) shared(Blr,nseg,Alr,keymask)
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
!$omp end parallel do

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




!> Tranform one wavefunction between Global region and localisation region
subroutine psi_to_locreg2(iproc, ldim, gdim, Llr, Glr, gpsi, lpsi)

  use module_base
  use module_types
 
 implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: iproc                  ! process ID
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
!  integer :: i_stat,i_all
  integer :: start,Gstart
  integer :: isegstart,istart

  call f_routine(id=subname)

! Define integers
  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0

! Initialize loc_psi
  call f_zero(lpsi)

! Get the keymask: shift for every segment of Llr (with respect to Glr)
! allocate(keymask(2,nseg),stat=i_stat)
  keymask = f_malloc((/ 2, nseg /),id='keymask')

  call shift_locreg_indexes(Glr,Llr,keymask,nseg)


!####################################################
! Do coarse region
!####################################################
  isegstart=1
  icheck = 0


!$omp parallel default(private) &
!$omp shared(icheck,lpsi,gpsi,Glr,Llr,keymask,lincrement,Gincrement,Gstart) &
!$omp firstprivate(isegstart,nseg)

  !$omp do reduction(+:icheck)
  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
     istart = llr%wfd%keyvloc(isegloc)-1
 
     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax) exit global_loop_c
        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c
        
        ! Define the offset between the two segments
        offset = lmin - Gmin
        if(offset < 0) then
           offset = 0
        end if
    
        ! Define the length of the two segments
        length = min(lmax,Gmax)-max(lmin,Gmin)

        icheck = icheck + (length + 1)
 
        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element

        do ix = 0,length
           istart = istart + 1
           lpsi(istart) = gpsi(Glr%wfd%keyvloc(isegG)+offset+ix)
        end do
     end do global_loop_c
  end do local_loop_c
  !$omp end do

! Check if the number of elements in loc_psi is valid
 ! if(icheck .ne. Llr%wfd%nvctr_c) then
   ! write(*,*)'There is an error in psi_to_locreg2: number of coarse points used',icheck
   ! write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
 ! end if

!##############################################################
! Now do fine region
!##############################################################

  !icheck = 0
  start = Llr%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c

  isegstart=Glr%wfd%nseg_c+1

  !$omp do reduction(+:icheck)
  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
     lmin = keymask(1,isegloc)
     lmax = keymask(2,isegloc)
     istart = llr%wfd%keyvloc(isegloc)-1
 
     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keygloc(1,isegG)
        Gmax = Glr%wfd%keygloc(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax)  exit global_loop_f
        if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        icheck = icheck + (length + 1)

        !Find the common elements and write them to the new localized wavefunction
        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
           istart = istart+1
           do igrid=1,7
              lpsi(start+(istart-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid)
           end do
        end do
     end do global_loop_f
  end do local_loop_f
  !$omp end do

  !$omp end parallel

 !! Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
    write(*,'(a,i0,a,i0)')'process ',iproc,': There is an error in psi_to_locreg: number of fine points used ',icheck
    write(*,'(a,i0)')'is not equal to the number of fine points in the region ',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
  end if



!  i_all=-product(shape(keymask))*kind(keymask)
! deallocate(keymask,stat=i_stat)
  call f_free(keymask)
  call f_release_routine()

END SUBROUTINE psi_to_locreg2


!> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
!! @warning       
!!    The quantity must not be stored in a compressed form.
subroutine global_to_local_parallel(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho,i1s,i1e,i2s,i2e,i3s,i3e,ni1,ni2, &
           i1shift, i2shift, i3shift, ise)

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
 integer,intent(in) :: i1s, i1e, i2s, i2e
 integer,intent(in) :: i3s, i3e ! starting and ending indices on z direction (related to distribution of rho when parallel)
 integer,intent(in) :: ni1, ni2 ! x and y extent of rho
 integer,intent(in) :: i1shift, i2shift, i3shift
 integer,dimension(6) :: ise

! Local variable
 integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
 integer :: indSmall, indSpin, indLarge ! indexes for the arrays
 integer :: ist2S,ist3S, ist2L, ist3L, istsa, ists, istl
 integer :: ii1shift, ii2shift, ii3shift, i1glob, i2glob, i3glob
 integer :: iii1, iii2, iii3

 write(*,'(a,8i8)') 'in global_to_local_parallel: i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2', i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2
 
 ! Cut out a piece of the quantity (rho) from the global region (rho) and
 ! store it in a local region (Lrho).
 indSmall=0
 indSpin=0
 ! Deactivate the spin for the moment
 do ispin=1,1!nspin
     !!$omp parallel default(none) &
     !!$omp shared(Glr,Llr,Lrho,rho,indSpin,i1s,i1e,i2s,i2e,i3s,i3e,ni1,ni2) &
     !!$omp private(ii1, ii2, ii3, i1, i2, i3, i1shift, i2shift, i3shift) &
     !!$omp private(ist3S, ist3L, istsa, ist2S, ist2L, ists, istl, indSmall, indLarge)
     !!$omp do
     do ii3=i3s,i3e
         !i3 = modulo(ii3-1,Glr%d%n3i)+1
         !i3glob = ii3+ise(5)-1
         i3glob = ii3+ise(5)-1
         i3 = modulo(i3glob-1,glr%d%n3i)+1
         !if (i3>modulo(i3e-1,glr%d%n3i)+1) then
         !if (i3glob>modulo(i3e-1,glr%d%n3i)+1) then
         if (modulo(ii3-1,glr%d%n3i)+1>modulo(i3e-1,glr%d%n3i)+1) then
             !This is a line before the wrap around, i.e. one needs a shift since 
             !the potential in the global region starts with the wrapped around part
             !i3shift = modulo(i3e-1,glr%d%n3i)
             !i3shift = min(i3shift,i3s-1)
             ii3shift = i3shift
         else
             ii3shift = 0
         end if
         if (i3glob<=glr%d%n3i) then
             iii3=ii3+i3shift
         else
             iii3=modulo(i3glob-1,glr%d%n3i)+1
         end if
         ist3S = (ii3-i3s)*Llr%d%n2i*Llr%d%n1i
         !ist3L = (i3+ii3shift-1)*ni2*ni1
         !ist3L = (i3+ii3shift-1)*ni2*ni1
         ist3L = (iii3-1)*ni2*ni1
         !write(*,'(a,10i8)') 'iproc, i3s, i3e, n3i, ii3, i3, i3glob, ii3shift, ist3l, nl', bigdft_mpi%iproc, i3s, i3e, glr%d%n3i, ii3, i3, i3glob, ii3shift, ist3l, size_rho
         write(*,'(a,11i8)') 'iproc, i3s, i3e, n3i, ii3, i3, i3glob, ii3shift, iii3, ist3l, nl', bigdft_mpi%iproc, i3s, i3e, glr%d%n3i, ii3, i3, i3glob, ii3shift, iii3, ist3l, glr%d%n1i*glr%d%n2i*glr%d%n3i
         istsa=ist3S-i1s+1
         do ii2=i2s,i2e
             !i2 = modulo(ii2-1,Glr%d%n2i)+1
             !i2glob = ii2+ise(3)-1
             i2glob = ii2+ise(3)-1
             i2 = modulo(i2glob-1,Glr%d%n2i)+1
             !if (i2>modulo(i2e-1,glr%d%n2i)+1) then
             !if (i2glob>modulo(i2e-1,glr%d%n2i)+1) then
             !!#if (i2>modulo(i2e-1,glr%d%n2i)+1) then
             !!#    !Same as above...
             !!#    !i2shift = modulo(i2e-1,glr%d%n2i)
             !!#    !i2shift = min(i2shift,i2s-1)
             !!#    ii2shift = i2shift
             !!#else
             !!#    ii2shift = 0
             !!#end if
             if (modulo(ii2-1,glr%d%n2i)+1>modulo(i2e-1,glr%d%n2i)+1) then
                 !This is a line before the wrap around, i.e. one needs a shift since 
                 !the potential in the global region starts with the wrapped around part
                 !i3shift = modulo(i3e-1,glr%d%n3i)
                 !i3shift = min(i3shift,i3s-1)
                 ii2shift = i2shift
             else
                 ii2shift = 0
             end if
             if (i2glob<=glr%d%n2i) then
                 iii2=ii2+i2shift
             else
                 iii2=modulo(i2glob-1,glr%d%n2i)+1
             end if
             ist2S = (ii2-i2s)*Llr%d%n1i 
             !ist2L = (i2+ii2shift-1)*ni1
             !ist2L = (i2+ii2shift-1)*ni1
             ist2L = (iii2-1)*ni1
             ists=istsa+ist2S
             istl=ist3L+ist2L
             do ii1=i1s,i1e
                 !i1 = modulo(ii1-1,Glr%d%n1i)+1
                 !i1glob = ii1+ise(1)-1
                 i1glob = ii1+ise(1)-1
                 i1 = modulo(i1glob-1,Glr%d%n1i)+1
                 !if (i1>modulo(i1e-1,glr%d%n1i)+1) then
                 !if (i1glob>modulo(i1e-1,glr%d%n1i)+1) then
                 !!#if (i1>modulo(i1e-1,glr%d%n1i)+1) then
                 !!#    !Same as above...
                 !!#    !i1shift = modulo(i1e-1,glr%d%n1i)
                 !!#    !i1shift = min(i1shift,i1s-1)
                 !!#    ii1shift = i1shift
                 !!#else
                 !!#    ii1shift = 0
                 !!#end if
                 if (modulo(ii1-1,glr%d%n1i)+1>modulo(i1e-1,glr%d%n1i)+1) then
                     !This is a line before the wrap around, i.e. one needs a shift since 
                     !the potential in the global region starts with the wrapped around part
                     !i3shift = modulo(i3e-1,glr%d%n3i)
                     !i3shift = min(i3shift,i3s-1)
                     ii1shift = i1shift
                 else
                     ii1shift = 0
                 end if
                 if (i1glob<=glr%d%n1i) then
                     iii1=ii1+i1shift
                 else
                     iii1=modulo(i1glob-1,glr%d%n1i)+1
                 end if
                 ! indSmall is the index in the local localization region
                 indSmall=ists+ii1
                 ! indLarge is the index in the global localization region. 
                 !indLarge= i1+ii1shift+istl
                 !indLarge= i1+ii1shift+istl
                 indLarge= iii1+istl
                 Lrho(indSmall)=rho(indLarge+indSpin)
                 write(666,'(a,8i8,2es16.8)') 'i1glob, i2glob, i3glob, iii1, iii2, iii3, indsmall, indlarge, val, testval', &
                     i1glob, i2glob, i3glob, iii1, iii2, iii3, indsmall, indlarge, Lrho(indSmall), real((i1glob+(i2glob-1)*glr%d%n1i+(i3glob-1)*glr%d%n1i*glr%d%n2i),kind=8)
             end do
         end do
     end do
     !!$omp end do
     !!$omp end parallel
     indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
 end do

END SUBROUTINE global_to_local_parallel

