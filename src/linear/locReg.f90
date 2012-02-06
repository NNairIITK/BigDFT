!determine a set of localisation regions from the centers and the radii.
!cut in cubes the global reference system
subroutine determine_locreg_periodic(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,calculateBounds)!,outofzone)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
  logical,dimension(nlr),intent(in):: calculateBounds
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg'
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz, calculate
  logical :: warningx,warningy,warningz
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3, iorb, jproc, jlr
  integer :: ierr 
  integer,dimension(3) :: outofzone
  real(gp) :: rx,ry,rz,cutoff  
  !!if (iproc == 0) then
  !!   write(*,*)'Inside determine_locreg_periodic:'
  !!end if


  !initialize out of zone and logicals
  outofzone (:) = 0     
  warningx = .false.
  warningy = .false.
  warningz = .false.  

  !determine the limits of the different localisation regions
  do ilr=1,nlr

     rx=cxyz(1,ilr)
     ry=cxyz(2,ilr)
     rz=cxyz(3,ilr)
     llr(ilr)%locregCenter(1)=rx
     llr(ilr)%locregCenter(2)=ry
     llr(ilr)%locregCenter(3)=rz

     cutoff=locrad(ilr)
     llr(ilr)%locrad=cutoff

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
!!!     if (iproc == 0 .and. verbose > 1) then
!!!        if ((iex - isx >= Glr%d%n1 - 14) .and. (warningx .eqv. .false.)) then
!!!           write(*,*)'Width of direction x :',(iex - isx)*hx,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n1*hx
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the '
!!!           write(*,*)'simulation box width. This is the only warning for x direction.'
!!!           warningx = .true.
!!!        end if
!!!        if ((iey - isy >= Glr%d%n2 - 14) .and. (warningy .eqv. .false.)) then
!!!           write(*,*)'Width of direction y :',(iey - isy)*hy,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n2*hy,'.'
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the width'
!!!           write(*,*)'of the simulation box. This is the only warning for y direction.'
!!!           warningy = .true.
!!!        end if
!!!        if ((iez - isz >= Glr%d%n3 - 14) .and. (warningz .eqv. .false.)) then
!!!           write(*,*)'Width of direction z :',(iez - isz)*hz,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n3*hz,'.'
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the width'
!!!           write(*,*)'of the simulation box. This is the only warning for z direction.'
!!!           warningz = .true.
!!!        end if 
!!!     end if

     ! Localization regions should always have free boundary conditions
     Llr(ilr)%geocode='F'

     !assign the starting/ending points and outofzone for the different
     ! geometries
     select case(Glr%geocode)
     case('F')
        isx=max(isx,Glr%ns1)
        isy=max(isy,Glr%ns2)
        isz=max(isz,Glr%ns3)

        iex=min(iex,Glr%ns1+Glr%d%n1)
        iey=min(iey,Glr%ns2+Glr%d%n2)
        iez=min(iez,Glr%ns3+Glr%d%n3)

     case('S')
        ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        isy=max(isy,Glr%ns2)
        iey=min(iey,Glr%ns2 + Glr%d%n2)
        outofzone(2) = 0

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if 
        end if

     case('P')
         ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        if (iey - isy >= Glr%d%n2) then       
           isy=Glr%ns2
           iey=Glr%ns2 + Glr%d%n2
         else
           isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
           iey= ln2 + isy
           if (iey > Glr%ns2+Glr%d%n2) then
              outofzone(2)=modulo(iey,Glr%d%n2+1)
           end if           
        end if

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if 
        end if
     end select

     !values for the starting point of the cube for wavelet grid
     Llr(ilr)%ns1=isx
     Llr(ilr)%ns2=isy
     Llr(ilr)%ns3=isz

     !dimensions of the localisation region
     Llr(ilr)%d%n1=iex-isx
     Llr(ilr)%d%n2=iey-isy
     Llr(ilr)%d%n3=iez-isz

     !assign outofzone
     Llr(ilr)%outofzone(:) = outofzone(:)

     ! Set the conditions for ext_buffers (conditions for buffer size)
     Gperx=(Glr%geocode /= 'F')
     Gpery=(Glr%geocode == 'P')
     Gperz=(Glr%geocode /= 'F')
     Lperx=(Llr(ilr)%geocode /= 'F')
     Lpery=(Llr(ilr)%geocode == 'P')
     Lperz=(Llr(ilr)%geocode /= 'F')

     !calculate the size of the buffers of interpolating function grid
     call ext_buffers(Gperx,Gnbl1,Gnbr1)
     call ext_buffers(Gpery,Gnbl2,Gnbr2)
     call ext_buffers(Gperz,Gnbl3,Gnbr3)
     call ext_buffers(Lperx,Lnbl1,Lnbr1)
     call ext_buffers(Lpery,Lnbl2,Lnbr2)
     call ext_buffers(Lperz,Lnbl3,Lnbr3)

     !starting point of the region for interpolating functions grid
     Llr(ilr)%nsi1= 2 * Llr(ilr)%ns1 - (Lnbl1 - Gnbl1)
     Llr(ilr)%nsi2= 2 * Llr(ilr)%ns2 - (Lnbl2 - Gnbl2)
     Llr(ilr)%nsi3= 2 * Llr(ilr)%ns3 - (Lnbl3 - Gnbl3)

     !dimensions of the fine grid inside the localisation region
     Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
     Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
     Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
     
     !NOTE: This will not work with symmetries (must change it)
     Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
     Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
     Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz

     !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
     Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
     Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
     Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31

!DEBUG
!     if (iproc == 0) then
!        write(*,*)'Description of zone:',ilr
!        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
!        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
!        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
!        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
!        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
!        write(*,*)'outofzone',ilr,':',outofzone(:)
!     end if
!DEBUG

    ! construct the wavefunction descriptors (wfd)
     call determine_wfd_periodicity(ilr,nlr,Glr,Llr)

     ! Sould check if nfu works properly... also relative to locreg!!
     !if the localisation region is isolated build also the bounds
     if (Llr(ilr)%geocode=='F') then
        ! Check whether the bounds shall be calculated. Do this only if the currect process handles
        ! orbitals in the current localization region.
        if(calculateBounds(ilr)) then
            call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                 Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                 Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
        end if
     end if
  end do !on ilr

END SUBROUTINE determine_locreg_periodic


!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg taking into account the pediodicity
!!          
!! WARNING: We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
!!         
!! SOURCE:
!!
subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  !#############################################
  !local variables
  !############################################
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  character(len=*), parameter :: subname='determine_wfd_periodicity'
  integer, allocatable :: keyvglob(:)  

   !starting point of locreg (always inside locreg)
   isdir(1) = Llr(ilr)%ns1
   isdir(2) = Llr(ilr)%ns2
   isdir(3) = Llr(ilr)%ns3
   !ending point of locreg (can be outside the simulation box)
   iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
   iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
   iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
   ! starting and ending point of fine grid in Global region
   Gifs(1) = Glr%d%nfl1 + Glr%ns1
   Gifs(2) = Glr%d%nfl2 + Glr%ns2
   Gifs(3) = Glr%d%nfl3 + Glr%ns3
   Gife(1) = Glr%d%nfu1 + Glr%ns1
   Gife(2) = Glr%d%nfu2 + Glr%ns2
   Gife(3) = Glr%d%nfu3 + Glr%ns3
   ! periodicity
   period(1) = Glr%d%n1
   period(2) = Glr%d%n2
   period(3) = Glr%d%n3

   ! Determine starting point of the fine grid in locreg
   do ii=1,3
      if (Llr(ilr)%outofzone(ii) > 0) then
         ! When periodicity, we must check for 2 different situations:
         ! (1) : starting of locreg before or in fine grid zone
         if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
         ! (2) : starting point after fine grid
         if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
      else
          Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
      end if 
   end do

   ! Determine ending point of the fine grid in locreg
   do ii=1,3
      if(Llr(ilr)%outofzone(ii) > 0) then
         !When periodicity, we must check for three different situations:
         ! (1) : ending of locreg before fine grid zone
         if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
         ! (2) : ending of locreg in fine grid zone
         if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
           Life(ii) = iedir(ii)-isdir(ii)
         end if
         ! (3) : ending of locreg after ending of fine grid zone
         if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
      else
         Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
      end if
   end do

   ! Assign values to Llr
   Llr(ilr)%d%nfl1 = Lifs(1)
   Llr(ilr)%d%nfl2 = Lifs(2)
   Llr(ilr)%d%nfl3 = Lifs(3)
   Llr(ilr)%d%nfu1 = Life(1)
   Llr(ilr)%d%nfu2 = Life(2)
   Llr(ilr)%d%nfu3 = Life(3)

   ! define the wavefunction descriptors inside the localisation region
   !coarse part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_c,Glr%wfd%nvctr_c,&
          Glr%wfd%keygloc(1,1),Glr%wfd%keyv(1),nseg_c,nvctr_c,Llr(ilr)%outofzone(:))
   !fine part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),nseg_f,nvctr_f,Llr(ilr)%outofzone(:))

   ! Assign the values to Llr
   Llr(ilr)%wfd%nseg_c = nseg_c
   Llr(ilr)%wfd%nseg_f = nseg_f
   Llr(ilr)%wfd%nvctr_c= nvctr_c
   Llr(ilr)%wfd%nvctr_f= nvctr_f

   !allocate the wavefunction descriptors following the needs
   call allocate_wfd(Llr(ilr)%wfd,subname)

   allocate(keyvglob(Llr(ilr)%wfd%nseg_c+Llr(ilr)%wfd%nseg_f))   !for tests should remove

   !Now, fill the descriptors:
   !coarse part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),&
        Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyv(1),&
        Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
        Llr(ilr)%wfd%keygloc(1,1),Llr(ilr)%wfd%keyglob(1,1),Llr(ilr)%wfd%keyv(1),&
        keyvglob(1),&
!        Llr(ilr)%wfd%keyvglob(1),&
        Llr(ilr)%outofzone(:))

   !fine part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
        Llr(ilr)%wfd%keygloc(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyglob(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        keyvglob(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
!        Llr(ilr)%wfd%keyvglob(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%outofzone(:))

   deallocate(keyvglob)

END SUBROUTINE determine_wfd_periodicity






subroutine determine_locregSphere(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,calculateBounds)!,outofzone)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
  logical,dimension(nlr),intent(in):: calculateBounds
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg'
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz, calculate
  logical :: warningx,warningy,warningz
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3, iorb, jproc, jlr
  integer :: ierr 
  integer,dimension(3) :: outofzone
  real(gp) :: rx,ry,rz,cutoff  


  !initialize out of zone and logicals
  outofzone (:) = 0     
  warningx = .false.
  warningy = .false.
  warningz = .false.  

  !determine the limits of the different localisation regions
  do ilr=1,nlr

     rx=cxyz(1,ilr)
     ry=cxyz(2,ilr)
     rz=cxyz(3,ilr)
     llr(ilr)%locregCenter(1)=rx
     llr(ilr)%locregCenter(2)=ry
     llr(ilr)%locregCenter(3)=rz

     cutoff=locrad(ilr)
     llr(ilr)%locrad=cutoff

     ! Determine the extrema of this localization regions (using only the coarse part, since this is always larger or equal than the fine part).
     call determine_boxbounds_sphere(glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, hx, hy, hz, &
          cutoff, llr(ilr)%locregCenter, &
           glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyv, isx, isy, isz, iex, iey, iez)

     ln1 = iex-isx
     ln2 = iey-isy
     ln3 = iez-isz


     ! Localization regions should always have free boundary conditions
     Llr(ilr)%geocode='F'

     !assign the starting/ending points and outofzone for the different
     ! geometries
     select case(Glr%geocode)
     case('F')
        isx=max(isx,Glr%ns1)
        isy=max(isy,Glr%ns2)
        isz=max(isz,Glr%ns3)

        iex=min(iex,Glr%ns1+Glr%d%n1)
        iey=min(iey,Glr%ns2+Glr%d%n2)
        iez=min(iez,Glr%ns3+Glr%d%n3)

     case('S')
        ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        isy=max(isy,Glr%ns2)
        iey=min(iey,Glr%ns2 + Glr%d%n2)
        outofzone(2) = 0

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if 
        end if

     case('P')
         ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        if (iey - isy >= Glr%d%n2) then       
           isy=Glr%ns2
           iey=Glr%ns2 + Glr%d%n2
         else
           isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
           iey= ln2 + isy
           if (iey > Glr%ns2+Glr%d%n2) then
              outofzone(2)=modulo(iey,Glr%d%n2+1)
           end if           
        end if

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if 
        end if
     end select

     !values for the starting point of the cube for wavelet grid
     Llr(ilr)%ns1=isx
     Llr(ilr)%ns2=isy
     Llr(ilr)%ns3=isz

     !dimensions of the localisation region
     Llr(ilr)%d%n1=iex-isx
     Llr(ilr)%d%n2=iey-isy
     Llr(ilr)%d%n3=iez-isz

     !assign outofzone
     Llr(ilr)%outofzone(:) = outofzone(:)

     ! Set the conditions for ext_buffers (conditions for buffer size)
     Gperx=(Glr%geocode /= 'F')
     Gpery=(Glr%geocode == 'P')
     Gperz=(Glr%geocode /= 'F')
     Lperx=(Llr(ilr)%geocode /= 'F')
     Lpery=(Llr(ilr)%geocode == 'P')
     Lperz=(Llr(ilr)%geocode /= 'F')

     !calculate the size of the buffers of interpolating function grid
     call ext_buffers(Gperx,Gnbl1,Gnbr1)
     call ext_buffers(Gpery,Gnbl2,Gnbr2)
     call ext_buffers(Gperz,Gnbl3,Gnbr3)
     call ext_buffers(Lperx,Lnbl1,Lnbr1)
     call ext_buffers(Lpery,Lnbl2,Lnbr2)
     call ext_buffers(Lperz,Lnbl3,Lnbr3)

     !starting point of the region for interpolating functions grid
     Llr(ilr)%nsi1= 2 * Llr(ilr)%ns1 - (Lnbl1 - Gnbl1)
     Llr(ilr)%nsi2= 2 * Llr(ilr)%ns2 - (Lnbl2 - Gnbl2)
     Llr(ilr)%nsi3= 2 * Llr(ilr)%ns3 - (Lnbl3 - Gnbl3)

     !dimensions of the fine grid inside the localisation region
     Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
     Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
     Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
     
     !NOTE: This will not work with symmetries (must change it)
     Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
     Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
     Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz

     !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
     Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
     Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
     Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31

!DEBUG
!     if (iproc == 0) then
!        write(*,*)'Description of zone:',ilr
!        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
!        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
!        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
!        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
!        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
!        write(*,*)'outofzone',ilr,':',outofzone(:)
!     end if
!DEBUG

    ! construct the wavefunction descriptors (wfd)
     call determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)

     ! Sould check if nfu works properly... also relative to locreg!!
     !if the localisation region is isolated build also the bounds
     if (Llr(ilr)%geocode=='F') then
        ! Check whether the bounds shall be calculated. Do this only if the currect process handles
        ! orbitals in the current localization region.
        if(calculateBounds(ilr)) then
            call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                 Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                 Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
        end if
     end if
  end do !on ilr

END SUBROUTINE determine_locregSphere



subroutine determine_locregSphere_parallel(iproc,nproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,calculateBounds)!,outofzone)
  use module_base
  use module_types
  use module_communicatetypes
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
  logical,dimension(nlr),intent(in):: calculateBounds
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg'
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz,communicate_bounds
  logical :: warningx,warningy,warningz
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3, iorb, jproc, jlr
  integer :: ierr, ii, istat, root
  integer,dimension(3) :: outofzone
  real(gp) :: rx,ry,rz,cutoff


  !initialize out of zone and logicals
  outofzone (:) = 0     
  warningx = .false.
  warningy = .false.
  warningz = .false.  

  ! Determine how many locregs one process handles at most
  ii=ceiling(dble(nlr)/dble(nproc))

  !determine the limits of the different localisation regions
  do ilr=1,nlr

     if(mod(ilr-1,nproc)==iproc) then
    
         rx=cxyz(1,ilr)
         ry=cxyz(2,ilr)
         rz=cxyz(3,ilr)
         llr(ilr)%locregCenter(1)=rx
         llr(ilr)%locregCenter(2)=ry
         llr(ilr)%locregCenter(3)=rz
    
         cutoff=locrad(ilr)
         llr(ilr)%locrad=cutoff
    
         ! Determine the extrema of this localization regions (using only the coarse part, since this is always larger or equal than the fine part).
         call determine_boxbounds_sphere(glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, hx, hy, hz, &
              cutoff, llr(ilr)%locregCenter, &
               glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyv, isx, isy, isz, iex, iey, iez)
    
         ln1 = iex-isx
         ln2 = iey-isy
         ln3 = iez-isz
    
    
         ! Localization regions should always have free boundary conditions
         Llr(ilr)%geocode='F'
    
         !assign the starting/ending points and outofzone for the different
         ! geometries
         select case(Glr%geocode)
         case('F')
            isx=max(isx,Glr%ns1)
            isy=max(isy,Glr%ns2)
            isz=max(isz,Glr%ns3)
    
            iex=min(iex,Glr%ns1+Glr%d%n1)
            iey=min(iey,Glr%ns2+Glr%d%n2)
            iez=min(iez,Glr%ns3+Glr%d%n3)
    
         case('S')
            ! Get starting and ending for x direction     
            if (iex - isx >= Glr%d%n1) then       
               isx=Glr%ns1
               iex=Glr%ns1 + Glr%d%n1
            else
               isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
               iex= ln1 + isx
               if (iex > Glr%ns1+Glr%d%n1) then
                  outofzone(1)=modulo(iex,Glr%d%n1+1)
               end if           
            end if
            
            ! Get starting and ending for y direction (perpendicular to surface)
            isy=max(isy,Glr%ns2)
            iey=min(iey,Glr%ns2 + Glr%d%n2)
            outofzone(2) = 0
    
            !Get starting and ending for z direction
            if (iez - isz >= Glr%d%n3) then
               isz=Glr%ns3 
               iez=Glr%ns3 + Glr%d%n3
            else
               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               iez= ln3 + isz
               if (iez > Glr%ns3+Glr%d%n3) then
                  outofzone(3)=modulo(iez,Glr%d%n3+1)
               end if 
            end if
    
         case('P')
             ! Get starting and ending for x direction     
            if (iex - isx >= Glr%d%n1) then       
               isx=Glr%ns1
               iex=Glr%ns1 + Glr%d%n1
            else
               isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
               iex= ln1 + isx
               if (iex > Glr%ns1+Glr%d%n1) then
                  outofzone(1)=modulo(iex,Glr%d%n1+1)
               end if           
            end if
            
            ! Get starting and ending for y direction (perpendicular to surface)
            if (iey - isy >= Glr%d%n2) then       
               isy=Glr%ns2
               iey=Glr%ns2 + Glr%d%n2
             else
               isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
               iey= ln2 + isy
               if (iey > Glr%ns2+Glr%d%n2) then
                  outofzone(2)=modulo(iey,Glr%d%n2+1)
               end if           
            end if
    
            !Get starting and ending for z direction
            if (iez - isz >= Glr%d%n3) then
               isz=Glr%ns3 
               iez=Glr%ns3 + Glr%d%n3
            else
               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               iez= ln3 + isz
               if (iez > Glr%ns3+Glr%d%n3) then
                  outofzone(3)=modulo(iez,Glr%d%n3+1)
               end if 
            end if
         end select
    
         !values for the starting point of the cube for wavelet grid
         Llr(ilr)%ns1=isx
         Llr(ilr)%ns2=isy
         Llr(ilr)%ns3=isz
    
         !dimensions of the localisation region
         Llr(ilr)%d%n1=iex-isx
         Llr(ilr)%d%n2=iey-isy
         Llr(ilr)%d%n3=iez-isz
    
         !assign outofzone
         Llr(ilr)%outofzone(:) = outofzone(:)
    
         ! Set the conditions for ext_buffers (conditions for buffer size)
         Gperx=(Glr%geocode /= 'F')
         Gpery=(Glr%geocode == 'P')
         Gperz=(Glr%geocode /= 'F')
         Lperx=(Llr(ilr)%geocode /= 'F')
         Lpery=(Llr(ilr)%geocode == 'P')
         Lperz=(Llr(ilr)%geocode /= 'F')
    
         !calculate the size of the buffers of interpolating function grid
         call ext_buffers(Gperx,Gnbl1,Gnbr1)
         call ext_buffers(Gpery,Gnbl2,Gnbr2)
         call ext_buffers(Gperz,Gnbl3,Gnbr3)
         call ext_buffers(Lperx,Lnbl1,Lnbr1)
         call ext_buffers(Lpery,Lnbl2,Lnbr2)
         call ext_buffers(Lperz,Lnbl3,Lnbr3)
    
         !starting point of the region for interpolating functions grid
         Llr(ilr)%nsi1= 2 * Llr(ilr)%ns1 - (Lnbl1 - Gnbl1)
         Llr(ilr)%nsi2= 2 * Llr(ilr)%ns2 - (Lnbl2 - Gnbl2)
         Llr(ilr)%nsi3= 2 * Llr(ilr)%ns3 - (Lnbl3 - Gnbl3)
    
         !dimensions of the fine grid inside the localisation region
         Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
         Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
         Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
         
         !NOTE: This will not work with symmetries (must change it)
         Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
         Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
         Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz
    
         !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
         Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
         Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
         Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31
    
    !DEBUG
    !     if (iproc == 0) then
    !        write(*,*)'Description of zone:',ilr
    !        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
    !        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
    !        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
    !        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
    !        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
    !        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
    !        write(*,*)'outofzone',ilr,':',outofzone(:)
    !     end if
    !DEBUG
    
        ! construct the wavefunction descriptors (wfd)
         call determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)
    
         ! Sould check if nfu works properly... also relative to locreg!!
         !if the localisation region is isolated build also the bounds
         if (Llr(ilr)%geocode=='F') then
            ! Check whether the bounds shall be calculated. Do this only if the currect process handles
            ! orbitals in the current localization region.
            !if(calculateBounds(ilr)) then
                call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                     Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                     Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
            !end if
         end if
     end if
  end do !on ilr


  ! Communicate the locregs
  do ilr=1,nlr
     root=mod(ilr-1,nproc)
     !!if(iproc==root) then
     !!    communicate_bounds=calculateBounds(ilr)
     !!end if
     !!call mpi_bcast(communicate_bounds, 1, mpi_logical, root, mpi_comm_world, ierr)
     !!if(iproc==root) write(*,'(a,3i5,2l3)') 'iproc, root, ilr, communicate_bounds, associated(llr(ilr)%bounds%sb%ibzzx_c)', iproc, root, ilr, communicate_bounds, associated(llr(ilr)%bounds%sb%ibzzx_c)
     !if(iproc==root) write(*,'(a,3i5,l3)') 'iproc, root, ilr, associated(llr(ilr)%bounds%sb%ibzzx_c)', iproc, root, ilr, associated(llr(ilr)%bounds%sb%ibzzx_c)
     !!call communicate_locreg_descriptors(iproc, root, llr(ilr), communicate_bounds)
     call communicate_locreg_descriptors(iproc, root, llr(ilr))
     if (Llr(ilr)%geocode=='F') then
        ! Check whether the bounds shall be calculated. Do this only if the currect process handles
        ! orbitals in the current localization region.
        if(.not.calculateBounds(ilr)) then
            !write(*,'(a,i0,a,i0)') 'process ',iproc,' deletes bounds for locreg ',ilr
            call deallocate_convolutions_bounds(llr(ilr)%bounds, subname)
        end if
    end if
  end do


END SUBROUTINE determine_locregSphere_parallel






!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg taking into account the pediodicity
!!          
!! WARNING: We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
!!         
!! SOURCE:
!!
subroutine determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  real(8),intent(in):: hx, hy, hz
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  !#############################################
  !local variables
  !############################################
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  character(len=*), parameter :: subname='determine_wfd_periodicity'

   !starting point of locreg (always inside locreg)
   isdir(1) = Llr(ilr)%ns1
   isdir(2) = Llr(ilr)%ns2
   isdir(3) = Llr(ilr)%ns3
   !ending point of locreg (can be outside the simulation box)
   iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
   iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
   iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
   ! starting and ending point of fine grid in Global region
   Gifs(1) = Glr%d%nfl1 + Glr%ns1
   Gifs(2) = Glr%d%nfl2 + Glr%ns2
   Gifs(3) = Glr%d%nfl3 + Glr%ns3
   Gife(1) = Glr%d%nfu1 + Glr%ns1
   Gife(2) = Glr%d%nfu2 + Glr%ns2
   Gife(3) = Glr%d%nfu3 + Glr%ns3
   ! periodicity
   period(1) = Glr%d%n1
   period(2) = Glr%d%n2
   period(3) = Glr%d%n3

   ! Determine starting point of the fine grid in locreg
   do ii=1,3
      if (Llr(ilr)%outofzone(ii) > 0) then
         ! When periodicity, we must check for 2 different situations:
         ! (1) : starting of locreg before or in fine grid zone
         if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
         ! (2) : starting point after fine grid
         if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
      else
          Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
      end if 
   end do

   ! Determine ending point of the fine grid in locreg
   do ii=1,3
      if(Llr(ilr)%outofzone(ii) > 0) then
         !When periodicity, we must check for three different situations:
         ! (1) : ending of locreg before fine grid zone
         if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
         ! (2) : ending of locreg in fine grid zone
         if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
           Life(ii) = iedir(ii)-isdir(ii)
         end if
         ! (3) : ending of locreg after ending of fine grid zone
         if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
      else
         Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
      end if
   end do

   ! Assign values to Llr
   Llr(ilr)%d%nfl1 = Lifs(1)
   Llr(ilr)%d%nfl2 = Lifs(2)
   Llr(ilr)%d%nfl3 = Lifs(3)
   Llr(ilr)%d%nfu1 = Life(1)
   Llr(ilr)%d%nfu2 = Life(2)
   Llr(ilr)%d%nfu3 = Life(3)

   ! define the wavefunction descriptors inside the localisation region
   !coarse part
   call num_segkeys_sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_c, Glr%wfd%keygloc(1,1), &
        Glr%wfd%keyv(1), &
        llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nvctr_c)

   !fine part
   call num_segkeys_sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        glr%wfd%nseg_f, Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        llr(ilr)%wfd%nseg_f, llr(ilr)%wfd%nvctr_f)


   !allocate the wavefunction descriptors following the needs
   call allocate_wfd(Llr(ilr)%wfd,subname)

   !Now, fill the descriptors:
   !coarse part
   call segkeys_Sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        llr(ilr)%wfd%nseg_c, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_c, Glr%wfd%keygloc(1,1), &
        Glr%wfd%keyv(1), &
        llr(ilr)%wfd%keygloc(1,1),llr(ilr)%wfd%keyglob(1,1), llr(ilr)%wfd%keyv(1))

   !fine part
   call segkeys_Sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        llr(ilr)%wfd%nseg_f, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_f, Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        llr(ilr)%wfd%keygloc(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyglob(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyv(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)))


END SUBROUTINE determine_wfdSphere


!#############################################################################################################################################
!!****f* BigDFT/num_segkeys_periodic
!#############################################################################################################################################
!! FUNCTION: Calculates the number of segments and elements in localisation region
!!          
!! WARNING:   
!!         
!! SOURCE:
!!
subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
  integer, dimension(3),intent(in) :: outofzone
  !local variables
  logical :: lseg,go1,go2,go3
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check

  nvctr_loc=0
  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0

  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     lseg=.false.
     ! overlap conditions if zone completely inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec)
     ! overlap conditions if zone as components in other periodic cells
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        nvctr_check=nvctr_check+1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)

        if (go1 .and. go2 .and. go3 ) then
           nvctr_loc=nvctr_loc+1
           if (.not. lseg) then
              nsrt=nsrt+1
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
     end if
  end do
  nseg_loc=nend

  !check
  if (nend /= nsrt) then
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif

  if (nvctr_check /= nvctr) then
     write(*,'(1x,a,2(i8))')&
          'ERROR: incorrect number of coarse points examined for reducing the localisation region',&
          nvctr_check,nvctr
     stop
  end if

END SUBROUTINE num_segkeys_periodic




subroutine num_segkeys_sphere(n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, nseg, nvctr)
  implicit none
  integer, intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nsegglob
  real(8),intent(in):: hx, hy, hz, locrad
  real(8),dimension(3),intent(in):: locregCenter
  integer,dimension(2,nsegglob),intent(in):: keygglob
  integer,dimension(nsegglob),intent(in):: keyvglob
  integer,intent(out):: nseg, nvctr
  !local variables
  logical :: segment
  integer :: i, i1, i2, i3, nstart, nend, i2old, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3
  real(8):: cut, dx,dy, dz


  nvctr=0
  nstart=0
  nend=0
  segment=.false.

  cut=locrad**2
  i2old=-1
  do iseg=1,nsegglob
      jj=keyvglob(iseg)
      j0=keygglob(1,iseg)
      j1=keygglob(2,iseg)
      ii=j0-1
      i3=ii/((n1+1)*(n2+1))
      ii=ii-i3*(n1+1)*(n2+1)
      i2=ii/(n1+1)
      i0=ii-i2*(n1+1)
      i1=i0+j1-j0

      ii2=i2+nl2glob
      ii3=i3+nl3glob

      do i=i0,i1
          ii1=i+nl1glob
          dz=((ii3*hz)-locregCenter(3))**2
          dy=((ii2*hy)-locregCenter(2))**2
          dx=((ii1*hx)-locregCenter(1))**2
          if(dx+dy+dz<=cut) then
              nvctr=nvctr+1
              if(.not.segment) then
                  nstart=nstart+1
                  segment=.true.
              end if
          else
              if(segment) then
                  nend=nend+1
                  segment=.false.
              end if
          end if
      end do
      !if(segment .and. i2/=i2old) then
      if(segment) then
          ! Always start a new segment if we come to a new line in y direction.
          nend=nend+1
          segment=.false.
      end if
      i2old=i2
  end do


  nseg=nstart

  !check
  if (nend /= nstart) then
     write(*,*) 'nend , nstart',nend,nstart
     stop 'nend <> nstart'
  endif


END SUBROUTINE num_segkeys_sphere





subroutine determine_boxbounds_sphere(n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, hx, hy, hz, locrad, locregCenter, &
           nsegglob, keygglob, keyvglob, ixmin, iymin, izmin, ixmax, iymax, izmax)
  implicit none
  integer, intent(in) :: n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, nsegglob
  real(8),intent(in):: hx, hy, hz, locrad
  real(8),dimension(3),intent(in):: locregCenter
  integer,dimension(2,nsegglob),intent(in):: keygglob
  integer,dimension(nsegglob),intent(in):: keyvglob
  integer,intent(out):: ixmin, iymin, izmin, ixmax, iymax, izmax
  !local variables
  logical :: segment
  integer :: i, i1, i2, i3, nstart, nend, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3
  real(8):: cut, dx,dy, dz
  !debug
  integer:: iiimin, isegmin
  iiimin=0
  isegmin=0

  ! Initialize the retun values
  ixmax=0
  iymax=0
  izmax=0
  ixmin=nl1glob+n1glob
  iymin=nl2glob+n2glob
  izmin=nl3glob+n3glob

  cut=locrad**2
  do iseg=1,nsegglob
      jj=keyvglob(iseg)
      j0=keygglob(1,iseg)
      j1=keygglob(2,iseg)
      ii=j0-1
      i3=ii/((n1glob+1)*(n2glob+1))
      ii=ii-i3*(n1glob+1)*(n2glob+1)
      i2=ii/(n1glob+1)
      i0=ii-i2*(n1glob+1)
      i1=i0+j1-j0

      ii2=i2+nl2glob
      ii3=i3+nl3glob

      do i=i0,i1
          ii1=i+nl1glob
          dz=((ii3*hz)-locregCenter(3))**2
          dy=((ii2*hy)-locregCenter(2))**2
          dx=((ii1*hx)-locregCenter(1))**2
          if(dx+dy+dz<=cut) then
              if(ii1>ixmax) ixmax=ii1
              if(ii2>iymax) iymax=ii2
              if(ii3>izmax) izmax=ii3
              if(ii1<ixmin) ixmin=ii1 ; iiimin=j0-1 ; isegmin=iseg
              if(ii2<iymin) iymin=ii2
              if(ii3<izmin) izmin=ii3
          end if
      end do
  end do


END SUBROUTINE determine_boxbounds_sphere





subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keygloc,keyglob,keyvloc,keyvglob,outofzone)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(3), intent(in) :: outofzone
  integer, dimension(nseg_loc), intent(out) :: keyvglob
  integer, dimension(nseg_loc), intent(out) :: keyvloc
  integer, dimension(2,nseg_loc), intent(out) :: keygloc
  integer, dimension(2,nseg_loc), intent(out) :: keyglob
  !local variables
  character(len=*),parameter :: subname = 'segkeys_periodic'
  logical :: go1,go2,go3,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l
  integer :: i_stat, i_all
  integer :: ngridp,ngridlob,loc
  integer, allocatable :: keyg_loc(:,:)

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc


  allocate(keyg_loc(2,nseg_loc),stat=i_stat)
  call memocc(i_stat,keyg_loc,'keyg_loc',subname)

  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     lseg=.false.

     ! intersection condition if zone inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec)
     ! intersection condition if zone has components outside simulation box (periodic)
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)
        if (go1 .and. go2 .and. go3) then
          !index of the compressed function
          i1l=i-i1sc
          if(outofzone(1) > 0 .and. i <= outofzone(1))i1l = i - i1sc + n1 + 1
          i2l=i2-i2sc
          if(outofzone(2) > 0 .and. i2 <= outofzone(2))i2l = i2 - i2sc + n2 + 1
          i3l=i3-i3sc
          if(outofzone(3) > 0 .and. i3 <= outofzone(3))i3l = i3 - i3sc + n3 + 1
          ngridp=i3l*((n1l+1)*(n2l+1)) + i2l*(n1l+1) + i1l+1
          ngridlob = i3 * ((n1+1)*(n2+1)) + i2 * (n1+1) + i + 1

          nvctr_check=nvctr_check+1
          if (.not. lseg) then
!             print *,'         check:',i,i2,i3,i1l,i2l,i3l,ngridp
             nsrt=nsrt+1
             keyg_loc(1,nsrt)=ngridp
             keyglob(1,nsrt)=ngridlob
             keyvglob(nsrt)=nvctr_check
          end if
          lseg=.true.
        else
           if (lseg) then
!              print *,'in        else:',i,i2,i3,i1l,i2l,i3l,ngridp
              nend=nend+1
              keyg_loc(2,nend)=ngridp
              keyglob(2,nend)=ngridlob
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
!        print *,'in second else:',i,i2,i3,i1l,i2l,i3l,ngridp
        nend=nend+1
        keyg_loc(2,nend)=ngridp
        keyglob(2,nend)=ngridlob
     end if
  end do

  !check
  if (nvctr_check /= nvctr_loc .or. nend /= nsrt .or. nend /= nseg_loc) then
     print *,'global region statistics:',nseg,nvctr
     write(*,*)&
          'ERROR: problem in segkeys_periodic  ',&
          'nvctr_check:',nvctr_check,'nvctr_loc:',nvctr_loc,&
          'nend:',nend,'nsrt:',nsrt,'nseg_loc:',nseg_loc
     stop
  end if

 ! Now build the keyvloc where we replace the segments in order for the loc
 do iseg = 1, nseg_loc
    !sorting the keyg_loc
    loc = minloc(keyg_loc(1,:),1)
    keygloc(1,iseg) = keyg_loc(1,loc)
    keygloc(2,iseg) = keyg_loc(2,loc)
    keyg_loc(1,loc) = maxval(keyg_loc) + 1
    keyvloc(iseg) = keyvglob(loc)
!    print *,'iseg,keyglob,keyvglob,keygloc,keyvloc',iseg,keyglob(1,iseg),keyvglob(iseg),keygloc(1,iseg),keyvloc(iseg)
 end do

 i_all = -product(shape(keyg_loc))*kind(keyg_loc)
 deallocate(keyg_loc, stat = i_stat)
 call memocc(i_stat,i_all,'keyg_loc',subname)

END SUBROUTINE segkeys_periodic
                                     
!> This routine generates the keyglob, keygloc and keyv for the localization regions using periodic boundary conditions
!! The keys are continous is the localization region (keygloc), while they are inverted global region (keyglob). 
!!subroutine segkeys_periodic_loc(Glr,Llr,logrid,keygloc,keyglob,keyv)
!!  implicit none
!!  type(locreg_descriptors), intent(in) :: Glr
!!  type(locreg_descriptors), intent(in) :: Llr
!!  integer, dimension(2,nseg), intent(out) :: keygloc
!!  integer, dimension(2,nseg), intent(out) :: keyglob
!!  integer, dimension(nseg), intent(out) :: keyv
!!  !local variables
!!  logical :: plogrid
!!  integer :: i1,i2,i3,i,j,k,nsrt,nend, mvctr
!!  integer :: ngridloc,ngridglob
!!
!!  mvctr = 0
!!  do k = 1, Llr%d%n3
!!     i3 = modulo(Llr%ns3 + k, Glr%d%n3)
!!     do j = 1, Llr%d%n2
!!        i2 = modulo(Llr%ns2 + k, Glr%d%n2)
!!        plogrid=.false.
!!        do i = 1, Llr%d%n1
!!           i1 = modulo(Llr%ns1 + k, Glr%d%n1)
!!           ngridloc = k*(Llr%d%n2+1)*(Llr%d%n1+1) + j*(Llr%d%n1+1) + i + 1
!!           ngridglob = i3*(Glr%d%n2+1)*(Glr%d%n1+1) + i2*(Glr%d%n1+1) + i1 +1
!!           if (logrid(i1,i2,i3)) then
!!              mvctr=mvctr+1
!!              if (.not. plogrid) then
!!                 nsrt=nsrt+1
!!                 keygloc(1,nsrt)=ngridloc
!!                 keyglob(1,nsrt)=ngridglob
!!                 keyv(nsrt)=mvctr
!!              endif
!!           else
!!              if (plogrid) then
!!                 nend=nend+1
!!                 keygloc(2,nend)=ngridloc-1
!!                 keyglob(2,nend)=ngridglob-1
!!              endif
!!           endif
!!           plogrid=logrid(i1,i2,i3)
!!        end do
!!        if (plogrid) then
!!           nend=nend+1
!!           keygloc(2,nend)=ngridloc
!!           keyglob(2,nend)=ngridglob
!!        endif
!!     end do
!!  end do
!!
!!  if (nend /= nsrt) then 
!!     write(*,*) 'nend , nsrt',nend,nsrt
!!     stop 'nend <> nsrt'
!!  endif
!!
!!END SUBROUTINE segkeys_periodic_loc



subroutine segkeys_Sphere(n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, keyg_loc, keyg_glob, keyv)
  implicit none
  integer,intent(in):: n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, nsegglob
  real(8):: hx, hy, hz, locrad
  real(8),dimension(3):: locregCenter
  integer,dimension(2,nsegglob),intent(in):: keygglob
  integer,dimension(nsegglob),intent(in):: keyvglob
  integer,dimension(2,nseg),intent(out):: keyg_loc, keyg_glob
  integer,dimension(nseg),intent(out):: keyv
  !local variables
  integer :: i, i1, i2, i3, nstart, nend, nvctr, igridpoint, igridglob, i2old, iseg, jj, j0, j1, ii, i0, n1l, n2l, n3l
  integer:: i1l, i2l, i3l, ii1, ii2, ii3
  real(8):: cut, dx, dy, dz
  logical:: segment


  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  !n1l=i1ec-i1sc
  !n2l=i2ec-i2sc
  !n3l=i3ec-i3sc
  n1l=nu1-nl1
  n2l=nu2-nl2
  n3l=nu3-nl3

  nvctr=0
  nstart=0
  nend=0
  segment=.false.

  cut=locrad**2
  i2old=-1
  do iseg=1,nsegglob
      jj=keyvglob(iseg)
      j0=keygglob(1,iseg)
      j1=keygglob(2,iseg)
      ii=j0-1
      i3=ii/((n1+1)*(n2+1))
      ii=ii-i3*(n1+1)*(n2+1)
      i2=ii/(n1+1)
      i0=ii-i2*(n1+1)
      i1=i0+j1-j0

      ii2=i2+nl2glob
      ii3=i3+nl3glob

      do i=i0,i1
          ii1=i+nl1glob
          dz=((ii3*hz)-locregCenter(3))**2
          dy=((ii2*hy)-locregCenter(2))**2
          dx=((ii1*hx)-locregCenter(1))**2
          i1l=ii1-nl1
          i2l=ii2-nl2
          i3l=ii3-nl3
          igridpoint=i3l*((n1l+1)*(n2l+1)) + i2l*(n1l+1) + i1l+1
          igridglob = ii3*(n1+1)*(n2+1) + ii2*(n1+1) + ii1 + 1 
          if(dx+dy+dz<=cut) then
              ! Check that we are not outside the global region
              if(ii1>nu1) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii1=',ii1,'>',nu1,'=nu1'
                  stop
              end if
              if(ii2>nu2) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii2=',ii2,'>',nu2,'=nu2'
                  stop
              end if
              if(ii3>nu3) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii3=',ii3,'>',nu3,'=nu3'
                  stop
              end if
              nvctr=nvctr+1
              if(.not.segment) then
                  nstart=nstart+1
                  keyg_loc(1,nstart)=igridpoint
                  keyg_glob(1,nstart)=igridglob
                  keyv(nstart)=nvctr
                  segment=.true.
              end if
          else
              if(segment) then
                  nend=nend+1
                  keyg_loc(2,nend)=igridpoint-1
                  keyg_glob(2,nend)=igridglob-1
                  segment=.false.
              end if
          end if
      end do
      if(segment) then
          ! Close the segment
          nend=nend+1
          keyg_loc(2,nend)=igridpoint
          keyg_glob(2,nend)=igridglob
          segment=.false.
      end if
      i2old=i2
  end do


  if (nend /= nstart) then
     write(*,*) 'nend , nstart',nend,nstart
     stop 'nend <> nstart'
  endif
  if (nseg /= nstart) then
     write(*,*) 'nseg , nstart',nseg,nstart
     stop 'nseg <> nstart'
  endif


END SUBROUTINE segkeys_Sphere







!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the number of intersection regions between locregs, taking into account the periodicity of the system.
!!          
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(out) :: isovrlp             ! Integer giving the number of overlaps (max 8 with periodicity)
  !########################################
  !Subroutine Array Arguments
  !########################################
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  !#############################################
  !local variables
  !############################################
  integer :: ii,azones,bzones,i_stat,i_all
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  character(len=*), parameter :: subname='get_number_of_overlap_region'
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
  end do

!write(*,*)'azones,bzones',azones,bzones
!write(*,*)'outofzone',alr,':',outofzone(:,alr)
!write(*,*)'outofzone',blr,':',outofzone(:,blr)

!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)

! Now check the number of overlapping zones
 isovrlp = 0
 do izones=1,azones
   do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
      if(go1 .and. go2 .and. go3) then
         isovrlp = isovrlp + 1
      end if
   end do
 end do

! Deallocation block
  i_all = -product(shape(astart))*kind(astart)
  deallocate(astart,stat=i_stat)
  call memocc(i_stat,i_all,'astart',subname)
  i_all = -product(shape(aend))*kind(aend)
  deallocate(aend,stat=i_stat)
  call memocc(i_stat,i_all,'aend',subname)
  i_all = -product(shape(bstart))*kind(bstart)
  deallocate(bstart,stat=i_stat)
  call memocc(i_stat,i_all,'bstart',subname)
  i_all = -product(shape(bend))*kind(bend)
  deallocate(bend,stat=i_stat)
  call memocc(i_stat,i_all,'bend',subname)
  
END SUBROUTINE get_number_of_overlap_region

!#############################################################################################################################################
!!****f* BigDFT/fracture_periodic_zone
!#############################################################################################################################################
!! FUNCTION: Divides the locreg into zones contained inside the simulation box, by applying the primitive vectors
!!           It returns: astart(3,nzones) which is the starting points of the different zones (max. 8)
!!                       aend(3,nzones) which is the ending points of the different zones (max. 8)
!!          
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine fracture_periodic_zone(nzones,Glr,Llr,outofzone,astart,aend)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: nzones
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(3),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  integer,dimension(3,nzones),intent(out) :: astart !
  integer,dimension(3,nzones),intent(out) :: aend !
  !#############################################
  !local variables
  !############################################
  integer :: ii,index,jj
  integer,dimension(3) :: alrs,alre,Gend,Gstart,period
  character(len=*), parameter :: subname='fracture_periodic_zone'
  
! Start and end of Global region
  Gstart(1) = Glr%ns1 
  Gstart(2) = Glr%ns2
  Gstart(3) = Glr%ns3  
  Gend(1) = Glr%ns1 + Glr%d%n1
  Gend(2) = Glr%ns2 + Glr%d%n2
  Gend(3) = Glr%ns3 + Glr%d%n3

! Periodicity of the system
  period(1) = Glr%d%n1 + 1
  period(2) = Glr%d%n2 + 1
  period(3) = Glr%d%n3 + 1

! Start and end of local region
  alrs(1) = Llr%ns1
  alrs(2) = Llr%ns2
  alrs(3) = Llr%ns3
  alre(1) = Llr%ns1 + Llr%d%n1
  alre(2) = Llr%ns2 + Llr%d%n2
  alre(3) = Llr%ns3 + Llr%d%n3

!assign the first zone (necessarily without shift) and initialize the rest
  do ii=1,3
     astart(ii,:) = alrs(ii)
     aend(ii,:) = min(Gend(ii),alre(ii))
  end do

!assign the other zones
  index = 2
  do ii=1,3
     if(outofzone(ii) > 0) then    !Translation: X,Y,Z
        astart(ii,index) =  Gstart(ii)
        aend(ii,index) = modulo(alre(ii),period(ii))
        index = index + 1
     end if 
     do jj=ii+1,3
        if(outofzone(ii) > 0 .and. outofzone(jj) > 0) then  !Translation: X+Y,X+Z,Y+Z
           astart(ii,index) = Gstart(ii)
           astart(jj,index) = Gstart(jj)
           aend(ii,index) = modulo(alre(ii),period(ii))
           aend(jj,index) = modulo(alre(jj),period(jj))
           index = index + 1
        end if
     end do
  end do

  if(outofzone(1) > 0 .and. outofzone(2) > 0 .and. outofzone(3) > 0 ) then ! Translation: X+Y+Z
     astart(1,index) = Gstart(1)
     astart(2,index) = Gstart(2)
     astart(3,index) = Gstart(3)
     aend(1,index) = modulo(alre(1),period(1))
     aend(2,index) = modulo(alre(2),period(2))
     aend(3,index) = modulo(alre(3),period(3))
  end if

END SUBROUTINE fracture_periodic_zone


!##############################################################################################################################################
!!****f* BigDFT/get_overlap_region
!##############################################################################################################################################
!! FUNCTION Given two localization regions, A and B, this routine returns a localization region corresponding to the intersection of A & B. 
!!
!! SOURCE
!!
subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(in) :: isovrlp              ! Number of overlap regions
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  !#############################################
  !local variables
  !############################################
  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
  character(len=*), parameter :: subname='get_overlap_region_periodic'
  !# NEW
  integer :: ii,azones,bzones,i_stat,i_all,index
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
  end do

!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)

! Now check the number of overlapping zones
  index = 0
  do izones=1,azones
    do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
      if(go1 .and. go2 .and. go3) then
        index = index + 1

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.
         
        ! Determine the limits of the overlap region
        isx = max(astart(1,izones),bstart(1,jzones))
        isy = max(astart(2,izones),bstart(2,jzones))
        isz = max(astart(3,izones),bstart(3,jzones))

        iex = min(aend(1,izones),bend(1,jzones))
        iey = min(aend(2,izones),bend(2,jzones))
        iez = min(aend(3,izones),bend(3,jzones))

!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!       This could change the values of the bounds, so do it here
!       for now, in sandbox,put free boundary to all zones
        Olr(index)%geocode = 'F'  

!       Values for the starting point of the cube
        Olr(index)%ns1 = isx
        Olr(index)%ns2 = isy
        Olr(index)%ns3 = isz

!       Dimensions of the overlap region
        Olr(index)%d%n1 = iex - isx 
        Olr(index)%d%n2 = iey - isy 
        Olr(index)%d%n3 = iez - isz 
    
!       Dimensions of the fine grid inside the overlap region
        if (isx < iex) then
           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
        else
           write(*,*)'1: Yet to be implemented (little effort?)'
           stop
        end if

        if (isy < iey) then
           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
        else
           write(*,*)'2: Yet to be implemented (little effort?)'
           stop
        end if

        if (isz < iez) then
           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
        else
           write(*,*)'3: Yet to be implemented (little effort?)'
           stop
        end if

!       Dimensions of the interpolating scaling function grid 
!       (geocode already taken into acount because it is simple)
        select case(Olr(index)%geocode)
        case('F')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
        case('S')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        case('P')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        end select
 
!       Now define the wavefunction descriptors inside the overlap region
!       First calculate the number of points and segments for the region
!       Coarse part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyv(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyv(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keygloc(1,1),Olr(index)%wfd%keyv(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,&
          Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keygloc(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
          Olr(index)%wfd%keyv(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))

!       If the localisation region is isolated build also the bounds
        if (Olr(index)%geocode=='F') then
           call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
             Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
             Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)
     
        end if
     end if ! go1 .and. go2 .and. go3
   end do !jzones
 end do !izones

! Deallocation block
  i_all = -product(shape(astart))*kind(astart)
  deallocate(astart,stat=i_stat)
  call memocc(i_stat,i_all,'astart',subname)
  i_all = -product(shape(aend))*kind(aend)
  deallocate(aend,stat=i_stat)
  call memocc(i_stat,i_all,'aend',subname)
  i_all = -product(shape(bstart))*kind(bstart)
  deallocate(bstart,stat=i_stat)
  call memocc(i_stat,i_all,'bstart',subname)
  i_all = -product(shape(bend))*kind(bend)
  deallocate(bend,stat=i_stat)
  call memocc(i_stat,i_all,'bend',subname)

! Check on the number of zones
  if (index /= isovrlp) then
      write(*,*)&
          'ERROR: problem in get_overlap_region_periodic ',&
          'index:',index,'not equal to isovrlp:',isovrlp,&
          'The number of overlap descriptors constructed does not',&
          'correspond to the number of overlap regions.'
     stop
  end if

END SUBROUTINE get_overlap_region_periodic
!%***

!determine a set of localisation regions from the centers and the radii.
!cut in cubes the global reference system
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
  character(len=*), parameter :: subname='check_linear_inputguess'
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

!##############################################################################################################################################
!!****f* BigDFT/get_overlap_region2
!##############################################################################################################################################
!! FUNCTION Given two localization regions, A and B, this routine returns a localization region corresponding to the intersection of A & B.
!!          This is the same as get_overlap_region_periodic, but does not allocate the bound arrays to save memory.
!!
!! SOURCE
!!
subroutine get_overlap_region_periodic2(alr,blr,Glr,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types

 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(in) :: isovrlp              ! Number of overlap regions
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors
  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  !#############################################
  !local variables
  !############################################
  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
  character(len=*), parameter :: subname='get_overlap_region_periodic2'
  !# NEW
  integer :: ii,azones,bzones,i_stat,i_all,index
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
  end do


!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)

! Now check the number of overlapping zones
  index = 0
  do izones=1,azones
    do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
      if(go1 .and. go2 .and. go3) then
        index = index + 1

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.

        ! Determine the limits of the overlap region
        isx = max(astart(1,izones),bstart(1,jzones))
        isy = max(astart(2,izones),bstart(2,jzones))
        isz = max(astart(3,izones),bstart(3,jzones))

        iex = min(aend(1,izones),bend(1,jzones))
        iey = min(aend(2,izones),bend(2,jzones))
        iez = min(aend(3,izones),bend(3,jzones))

!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!       This could change the values of the bounds, so do it here
!       for now, in sandbox,put free boundary to all zones
        Olr(index)%geocode = 'F'

!       Values for the starting point of the cube
        Olr(index)%ns1 = isx
        Olr(index)%ns2 = isy
        Olr(index)%ns3 = isz

!       Dimensions of the overlap region
        Olr(index)%d%n1 = iex - isx
        Olr(index)%d%n2 = iey - isy
        Olr(index)%d%n3 = iez - isz

!       Dimensions of the fine grid inside the overlap region
        if (isx < iex) then
           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
        else
           write(*,*)'4: Yet to be implemented (little effort?), isx, iex', isx, iex
           stop
        end if

        if (isy < iey) then
           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
        else
           write(*,*)'5: Yet to be implemented (little effort?)'
           stop
        end if

        if (isz < iez) then
           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
        else
           write(*,*)'6: Yet to be implemented (little effort?)'
           stop
        end if

!       Dimensions of the interpolating scaling function grid
!       (geocode already taken into acount because it is simple)
        select case(Olr(index)%geocode)
        case('F')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
        case('S')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        case('P')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        end select

!       Now define the wavefunction descriptors inside the overlap region
!       First calculate the number of points and segments for the region
!       Coarse part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyv(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyv(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keygloc(1,1),Olr(index)%wfd%keyv(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,&
          Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keygloc(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
          Olr(index)%wfd%keyv(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))

!       If the localisation region is isolated build also the bounds
        !!!if (Olr(index)%geocode=='F') then
        !!!   call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
        !!!     Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
        !!!     Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)

        !!!end if
     end if ! go1 .and. go2 .and. go3
   end do !jzones
 end do !izones

! Deallocation block
  i_all = -product(shape(astart))*kind(astart)
  deallocate(astart,stat=i_stat)
  call memocc(i_stat,i_all,'astart',subname)
  i_all = -product(shape(aend))*kind(aend)
  deallocate(aend,stat=i_stat)
  call memocc(i_stat,i_all,'aend',subname)
  i_all = -product(shape(bstart))*kind(bstart)
  deallocate(bstart,stat=i_stat)
  call memocc(i_stat,i_all,'bstart',subname)
  i_all = -product(shape(bend))*kind(bend)
  deallocate(bend,stat=i_stat)
  call memocc(i_stat,i_all,'bend',subname)

! Check on the number of zones
  if (index /= isovrlp) then
      write(*,'(a,i0,a,i0,a)')&
          'ERROR: problem in get_overlap_region_periodic; index: ',index,' not equal to isovrlp: ',isovrlp,'&
          &. The number of overlap descriptors constructed does not correspond to the number of overlap regions.'
     stop
  end if

END SUBROUTINE get_overlap_region_periodic2




subroutine get_overlap_region_periodic2Sphere(alr,blr,Glr,hx,hy,hz,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types

 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  real(8),intent(in):: hx, hy, hz
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(in) :: isovrlp              ! Number of overlap regions
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors
  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  !#############################################
  !local variables
  !############################################
  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
  character(len=*), parameter :: subname='get_overlap_region_periodic2'
  !# NEW
  integer :: ii,azones,bzones,i_stat,i_all,index
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
  end do


!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)

! Now check the number of overlapping zones
  index = 0
  do izones=1,azones
    do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
      if(go1 .and. go2 .and. go3) then
        index = index + 1

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.

        ! Determine the limits of the overlap region
        isx = max(astart(1,izones),bstart(1,jzones))
        isy = max(astart(2,izones),bstart(2,jzones))
        isz = max(astart(3,izones),bstart(3,jzones))

        iex = min(aend(1,izones),bend(1,jzones))
        iey = min(aend(2,izones),bend(2,jzones))
        iez = min(aend(3,izones),bend(3,jzones))

!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!       This could change the values of the bounds, so do it here
!       for now, in sandbox,put free boundary to all zones
        Olr(index)%geocode = 'F'

!       Values for the starting point of the cube
        Olr(index)%ns1 = isx
        Olr(index)%ns2 = isy
        Olr(index)%ns3 = isz

!       Dimensions of the overlap region
        Olr(index)%d%n1 = iex - isx
        Olr(index)%d%n2 = iey - isy
        Olr(index)%d%n3 = iez - isz

!       Dimensions of the fine grid inside the overlap region
        if (isx < iex) then
           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
        else
           write(*,*)'4: Yet to be implemented (little effort?), isx, iex', isx, iex
           stop
        end if

        if (isy < iey) then
           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
        else
           write(*,*)'5: Yet to be implemented (little effort?)'
           stop
        end if

        if (isz < iez) then
           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
        else
           write(*,*)'6: Yet to be implemented (little effort?)'
           stop
        end if

!       Dimensions of the interpolating scaling function grid
!       (geocode already taken into acount because it is simple)
        select case(Olr(index)%geocode)
        case('F')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
        case('S')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        case('P')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        end select

!       Now define the wavefunction descriptors inside the overlap region
!       First calculate the number of points and segments for the region
!       Coarse part:
        !call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
        ! Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
        ! Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
        !!call num_segkeys_overlapSphere(olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
        !!     olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
        !!     olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
        !!     hx, hy, hz, llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
        !!     olr(index)%wfd%nseg_c, olr(index)%wfd%nvctr_c)
        call num_segkeys_sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
             hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
             llr(alr)%wfd%nseg_c, llr(alr)%wfd%keygloc(1,1), &
             llr(alr)%wfd%keyv(1), &
             olr(index)%wfd%nseg_c, olr(index)%wfd%nvctr_c)
        
!       Fine part:
        !call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
        ! Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        ! Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        ! Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        ! Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)
        !!call num_segkeys_overlapSphere(olr(index)%d%nfl1, olr(index)%d%nfu1, &
        !!     olr(index)%d%nfl2, olr(index)%d%nfu2, &
        !!     olr(index)%d%nfl3, olr(index)%d%nfu3, &
        !!     hx, hy, hz, llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
        !!     olr(index)%wfd%nseg_f, olr(index)%wfd%nvctr_f)
        call num_segkeys_sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
             hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
             llr(alr)%wfd%nseg_f, llr(alr)%wfd%keygloc(1,llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
             llr(alr)%wfd%keyv(llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
             olr(index)%wfd%nseg_f, olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        !call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
        !  Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
        !  Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
        !  Olr(index)%wfd%keyg(1,1),Olr(index)%wfd%keyv(1))
        !!call segkeys_overlapSphere(olr(index)%d%n1, olr(index)%d%n2, olr(index)%d%n3, &
        !!     olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
        !!     olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
        !!     olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
        !!     olr(index)%wfd%nseg_c, hx, hy, hz, &
        !!     llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
        !!     olr(index)%wfd%keyg(1,1), olr(index)%wfd%keyv(1))
        call segkeys_Sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
             olr(index)%wfd%nseg_c, hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
             llr(alr)%wfd%nseg_c, llr(alr)%wfd%keygloc(1,1), &
             llr(alr)%wfd%keyv(1), &
             olr(index)%wfd%keygloc(1,1),olr(index)%wfd%keyglob(1,1), &
             olr(index)%wfd%keyv(1))

!       Fine part
        !call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
        !  Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        !  Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        !  Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        !  Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
        !  Olr(index)%wfd%keyg(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
        !  Olr(index)%wfd%keyv(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))
        !!call segkeys_overlapSphere(olr(index)%d%n1, olr(index)%d%n2, olr(index)%d%n3, &
        !!     olr(index)%d%nfl1, olr(index)%d%nfu1, &
        !!     olr(index)%d%nfl2, olr(index)%d%nfu2, &
        !!     olr(index)%d%nfl3, olr(index)%d%nfu3, &
        !!     olr(index)%wfd%nseg_f, hx, hy, hz, &
        !!     llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
        !!     olr(index)%wfd%keyg(1,olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
        !!     olr(index)%wfd%keyv(olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)))
        call segkeys_Sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
             olr(index)%wfd%nseg_f, hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
             llr(alr)%wfd%nseg_f, llr(alr)%wfd%keygloc(1,llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
             llr(alr)%wfd%keyv(llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
             olr(index)%wfd%keygloc(1,olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
             olr(index)%wfd%keyglob(1,olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
             olr(index)%wfd%keyv(olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)))

!       If the localisation region is isolated build also the bounds
        !!!if (Olr(index)%geocode=='F') then
        !!!   call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
        !!!     Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
        !!!     Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)

        !!!end if
     end if ! go1 .and. go2 .and. go3
   end do !jzones
 end do !izones

! Deallocation block
  i_all = -product(shape(astart))*kind(astart)
  deallocate(astart,stat=i_stat)
  call memocc(i_stat,i_all,'astart',subname)
  i_all = -product(shape(aend))*kind(aend)
  deallocate(aend,stat=i_stat)
  call memocc(i_stat,i_all,'aend',subname)
  i_all = -product(shape(bstart))*kind(bstart)
  deallocate(bstart,stat=i_stat)
  call memocc(i_stat,i_all,'bstart',subname)
  i_all = -product(shape(bend))*kind(bend)
  deallocate(bend,stat=i_stat)
  call memocc(i_stat,i_all,'bend',subname)

! Check on the number of zones
  if (index /= isovrlp) then
      write(*,'(a,i0,a,i0,a)')&
          'ERROR: problem in get_overlap_region_periodic; index: ',index,' not equal to isovrlp: ',isovrlp,'&
          &. The number of overlap descriptors constructed does not correspond to the number of overlap regions.'
     stop
  end if

END SUBROUTINE get_overlap_region_periodic2Sphere






subroutine parallel_repartition_locreg(iproc,nproc,nlr,nlr_par,islr_par)
  implicit none
  integer, intent(in) :: iproc,nproc,nlr
  integer, dimension(0:nproc-1), intent(out) :: nlr_par
  integer, dimension(0:nproc-1), intent(out) :: islr_par
  !local variables
  integer :: jproc,numlr,difflr,ind

  numlr = int(nlr/nproc)
  difflr = nlr - numlr*nproc

  nlr_par = numlr
  ind = 0
  do jproc=0,difflr-1
     nlr_par(jproc) = nlr_par(jproc) + 1   
     islr_par(jproc) = ind
     ind = ind + nlr_par(jproc)
  end do
  do jproc=difflr,nproc
     islr_par(jproc) = ind
     ind = ind + nlr_par(jproc)
  end do

END SUBROUTINE parallel_repartition_locreg


!!subroutine make_LLr_MpiType(Llr,nlr,mpiLlr)
!!use module_types
!!
!!implicit none
!!
!!integer,intent(in) :: nlr
!!type(locreg_descriptors), intent(in) :: Llr
!!integer, intent(out) :: mpiLlr
!!!Local variables
!!integer :: ierr
!!integer :: Llr_add
!!integer :: geo_add, hyb_add, ns1_add, ns2_add,ns3_add
!!integer :: n1_add,n2_add,n3_add,nfl1_add,nfl2_add,nfl3_add
!!integer :: nfu1_add,nfu2_add,nfu3_add,n1i_add,n2i_add,n3i_add
!!integer :: nvctrc_add,nvctrf_add,nsegc_add,nsegf_add
!!integer :: keyg_add,keyv_add
!!integer :: localnorb_add,outofzone_add,projflg_add
!!integer :: ibyzc_add,ibxzc_add,ibxyc_add,ibyzf_add,ibxzf_add,ibxyf_add
!!integer :: ibzzxc_add,ibyyzzc_add,ibxyff_add,ibzzxf_add,ibyyzzf_add
!!integer :: ibzxxc_add,ibxxyyc_add,ibyzff_add,ibzxxf_add,ibxxyyf_add
!!integer :: ibyyzzr_add
!!integer :: geo_off,hyb_off,ns1_off,ns2_off,ns3_off,n1_off,n2_off,n3_off,nfl1_off,nfl2_off,nfl3_off
!!integer :: nfu1_off,nfu2_off,nfu3_off,n1i_off,n2i_off,n3i_off,nvctrc_off,nvctrf_off,nsegc_off,nsegf_off 
!!integer :: keyg_off,keyv_off,ibyzc_off,ibxzc_off,ibxyc_off,ibyzf_off,ibxzf_off,ibxyf_off,ibzzxc_off
!!integer :: ibyyzzc_off,ibxyff_off,ibzzxf_off,ibyyzzf_off,ibzxxc_off,ibxxyyc_off,ibyzff_off,ibzxxf_off
!!integer :: ibxxyyf_off,ibyyzzr_off,localnorb_off,outofzone_off,projflg_off
!!
!! ! Get all the adresses of the components of Llr
!!  call MPI_Get_address(Llr, Llr_add,ierr)
!!  call MPI_Get_address(Llr%geocode, geo_add,ierr)
!!  call MPI_Get_address(Llr%hybrid_on, hyb_add,ierr)
!!  call MPI_Get_address(Llr%ns1, ns1_add,ierr)
!!  call MPI_Get_address(Llr%ns2, ns2_add,ierr)
!!  call MPI_Get_address(Llr%ns3, ns3_add,ierr)
!!  call MPI_Get_address(Llr%localnorb, localnorb_add,ierr)
!!  call MPI_Get_address(Llr%outofzone, outofzone_add,ierr)
!!  call MPI_Get_address(Llr%projflg, projflg_add,ierr)
!!  call MPI_Get_address(Llr%d%n1, n1_add,ierr)
!!  call MPI_Get_address(Llr%d%n2, n2_add,ierr)
!!  call MPI_Get_address(Llr%d%n3, n3_add,ierr)
!!  call MPI_Get_address(Llr%d%nfl1, nfl1_add,ierr)
!!  call MPI_Get_address(Llr%d%nfl2, nfl2_add,ierr)
!!  call MPI_Get_address(Llr%d%nfl3, nfl3_add,ierr)
!!  call MPI_Get_address(Llr%d%nfu1, nfu1_add,ierr)
!!  call MPI_Get_address(Llr%d%nfu2, nfu2_add,ierr)
!!  call MPI_Get_address(Llr%d%nfu3, nfu3_add,ierr)
!!  call MPI_Get_address(Llr%wfd%nvctr_c, nvctrc_add,ierr)
!!  call MPI_Get_address(Llr%wfd%nvctr_f, nvctrf_add,ierr)
!!  call MPI_Get_address(Llr%wfd%nseg_c, nsegc_add,ierr)
!!  call MPI_Get_address(Llr%wfd%nseg_f, nsegf_add,ierr)
!!  call MPI_Get_address(Llr%wfd%keyg, keyg_add,ierr)
!!  call MPI_Get_address(Llr%wfd%keyv, keyv_add,ierr)
!!  call MPI_Get_address(Llr%bounds%kb%ibyz_c, ibyzc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%kb%ibxz_c, ibxzc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%kb%ibxy_c, ibxyc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%kb%ibyz_f, ibyzf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%kb%ibxz_f, ibxzf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%kb%ibxy_f, ibxyf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%sb%ibzzx_c, ibzzxc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%sb%ibyyzz_c, ibyyzzc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%sb%ibxy_ff, ibxyff_add,ierr)
!!  call MPI_Get_address(Llr%bounds%sb%ibzzx_f, ibzzxf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%sb%ibyyzz_f, ibyyzzf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%gb%ibzxx_c, ibzxxc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%gb%ibxxyy_c, ibxxyyc_add,ierr)
!!  call MPI_Get_address(Llr%bounds%gb%ibyz_ff, ibyzff_add,ierr)
!!  call MPI_Get_address(Llr%bounds%gb%ibzxx_f, ibzxxf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%gb%ibxxyy_f, ibxxyyf_add,ierr)
!!  call MPI_Get_address(Llr%bounds%ibyyzz_r, ibyyzzr_add,ierr)
!!
!!!Calculate offsets
!!  geo_off = geo_add - Llr_add
!!  hyb_off = hyb_add - Llr_add
!!  ns1_off = ns1_add - Llr_add
!!  ns2_off = ns2_add - Llr_add
!!  ns3_off = ns3_add - Llr_add
!!  n1_off = n1_add - Llr_add
!!  n2_off = n2_add - Llr_add
!!  n3_off = n3_add - Llr_add
!!  localnorb_off = localnorb_add - Llr_add
!!  outofzone_off = outofzone_add - Llr_add
!!  projflg_off = projflg_add - Llr_add
!!  nfl1_off = nfl1_add - Llr_add
!!  nfl2_off = nfl2_add - Llr_add
!!  nfl3_off = nfl3_add - Llr_add
!!  nfu1_off = nfu1_add - Llr_add
!!  nfu2_off = nfu2_add - Llr_add
!!  nfu3_off = nfu3_add - Llr_add
!!  n1i_off = n1i_add - Llr_add
!!  n2i_off = n2i_add - Llr_add
!!  n3i_off = n3i_add - Llr_add
!!  nvctrc_off = nvctrc_add - Llr_add
!!  nvctrf_off = nvctrf_add - Llr_add
!!  nsegc_off = nsegc_add - Llr_add
!!  nsegf_off = nsegf_add - Llr_add
!!  keyg_off = keyg_add - Llr_add
!!  keyv_off = keyv_add - Llr_add
!!  ibyzc_off = ibyzc_add - Llr_add
!!  ibxzc_off = ibxzc_add - Llr_add
!!  ibxyc_off = ibxyc_add - Llr_add
!!  ibyzf_off = ibyzf_add - Llr_add
!!  ibxzf_off = ibxzf_add - Llr_add
!!  ibxyf_off = ibxyf_add - Llr_add
!!  ibzzxc_off = ibzzxc_add - Llr_add
!!  ibyyzzc_off = ibyyzzc_add - Llr_add
!!  ibxyff_off = ibxyff_add - Llr_add
!!  ibzzxf_off = ibzzxf_add - Llr_add
!!  ibyyzzf_off = ibyyzzf_add - Llr_add
!!  ibzxxc_off  = ibzxxc_add - Llr_add
!!  ibxxyyc_off = ibxxyyc_add - Llr_add
!!  ibyzff_off  = ibyzff_add - Llr_add
!!  ibzxxf_off  = ibzxxf_add - Llr_add
!!  ibxxyyf_off = ibxxyyf_add - Llr_add
!!  ibyyzzr_off = ibyyzzr_add - Llr_add
!!
!!! Give sizes of all the components
!!  geo_siz = kind(Llr%geocode)
!!  hyb_siz = kind(Llr%hybrid_on)
!!  ns1_siz = kind(Llr%ns1)
!!  ns2_siz = kind(Llr%ns2)
!!  ns3_siz = kind(Llr%ns3)
!!  n1_siz = kind(Llr%d%n1)
!!  n2_siz = kind(Llr%d%n2)
!!  n3_siz = kind(Llr%d%n3)
!!  localnorb_siz = product(shape(Llr%Localnorb))*kind(Llr%Localnorb)
!!  outofzone_siz = product(shape(Llr%outofzone))*kind(Llr%outofzone)
!!  projflg_siz = product(shape(Llr%projflg))*kind(Llr%projflg)
!!  nfl1_siz = kind(Llr%d%nfl1)
!!  nfl2_siz = kind(Llr%d%nfl2)
!!  nfl3_siz = kind(Llr%d%nfl3)
!!  nfu1_siz = kind(Llr%d%nfu1)
!!  nfu2_siz = kind(Llr%d%nfu2)
!!  nfu3_siz = kind(Llr%d%nfu3)
!!  n1i_siz = kind(Llr%d%n1i)
!!  n2i_siz = kind(Llr%d%n2i)
!!  n3i_siz = kind(Llr%d%n3i)
!!  nvctrc_siz = kind(Llr%wfd%nvctr_c)
!!  nvctrf_siz = kind(Llr%wfd%nvctr_f)
!!  nsegc_siz = kind(Llr%wfd%nseg_c)
!!  nsegf_siz = kind(Llr%wfd%nseg_f)
!!  keyg_siz = product(shape(Llr%keyg))*kind(Llr%wfd%keyg)
!!  keyv_siz = product(shape(Llr%keyv))*kind(Llr%wfd%keyv)
!!  ibyzc_siz = product(shape(Llr%bounds%kb%ibyz_c))*kind(Llr%bounds%kb%ibyz_c)
!!  ibxzc_siz = product(shape(Llr%bounds%kb%ibxz_c))*kind(Llr%bounds%kb%ibxz_c)
!!  ibxyc_siz = product(shape(Llr%bounds%kb%ibxy_c))*kind(Llr%bounds%kb%ibxy_c)
!!  ibyzf_siz = product(shape(Llr%bounds%kb%ibyz_f))*kind(Llr%bounds%kb%ibyz_f)
!!  ibxzf_siz = product(shape(Llr%bounds%kb%ibxz_f))*kind(Llr%bounds%kb%ibxz_f)
!!  ibxyf_siz = product(shape(Llr%bounds%kb%ibxyf))*kind(Llr%bounds%kb%ibxyf)
!!  ibzzxc_siz = product(shape(Llr%bounds%sb%ibzzxc))*kind(Llr%bounds%sb%ibzzxc)
!!  ibyyzzc_siz = product(shape(Llr%bounds%sb%ibyyzz_c))*kind(Llr%bounds%sb%ibyyzz_c)
!!  ibxyff_siz = product(shape(Llr%bounds%sb%ibxy_ff))*kind(Llr%bounds%sb%ibxy_ff)
!!  ibzzxf_siz = product(shape(Llr%bounds%sb%ibzzx_f))*kind(Llr%bounds%sb%ibzzx_f)
!!  ibyyzzf_siz = product(shape(Llr%bounds%sb%ibyyzz_f))*kind(Llr%bounds%sb%ibyyzz_f)
!!  ibzxxc_siz  = product(shape(Llr%bounds%gb%ibzxx_c))*kind(Llr%bounds%gb%ibzxx_c)
!!  ibxxyyc_siz = product(shape(Llr%bounds%gb%ibxxyy_c))*kind(Llr%bounds%gb%ibxxyy_c)
!!  ibyzff_siz  = product(shape(Llr%bounds%gb%ibyz_ff))*kind(Llr%bounds%gb%ibyz_ff)
!!  ibzxxf_siz  = product(shape(Llr%bounds%gb%ibzxx_f))*kind(Llr%bounds%gb%ibzxx_f)
!!  ibxxyyf_siz = product(shape(Llr%bounds%gb%ibxxyy_f))*kind(Llr%bounds%gb%ibxxyy_f)
!!  ibyyzzr_siz = product(shape(Llr%bounds%ibyyzz_r))*kind(Llr%bounds%ibyyzz_r)
!!
!!! Arrays
!!  array_of_block_lengths = (/ geo_siz,hyb_siz,ns1_siz,ns2_siz,ns3_siz,localnorb_siz,outofzone_siz,projflg_siz,&
!!                              n1_siz,n2_siz,n3_siz,nfl1_siz,nfl2_siz,nfl3_siz,nfu1_siz,nfu2_siz,nfu3_siz,&
!!                              n1i_siz,n2i_siz,n3i_siz,nvctrc_siz,nvctrf_siz,nsegc_siz,nsegf_siz,keyg_siz,keyv_siz,&
!!                              ibyzc_siz,ibxzc_siz,ibxyc_siz,ibyzf_siz,ibxzf_siz,ibxyf_siz,ibzzxc_siz,ibyyzzc_siz,&
!!                              ibxyff_siz,ibzzxf_siz,ibyyzzf_siz,ibzxxc_siz,ibxxyyc_siz,ibyzff_siz,ibzxxf_siz,ibxxyyf_siz,&
!!                              ibyyzzr_siz /)
!!  array_of_displacements = (/ geo_off,hyb_off,ns1_off,ns2_off,ns3_off,localnorb_off,outofzone_off,projflg_off,&
!!                              n1_off,n2_off,n3_off,nfl1_off,nfl2_off,nfl3_off,nfu1_off,nfu2_off,nfu3_off,&
!!                              n1i_off,n2i_off,n3i_off,nvctrc_off,nvctrf_off,nsegc_off,nsegf_off,keyg_off,keyv_off,& 
!!                              ibyzc_off,ibxzc_off,ibxyc_off,ibyzf_off,ibxzf_off,ibxyf_off,ibzzxc_off,ibyyzzc_off,& 
!!                              ibxyff_off,ibzzxf_off,ibyyzzf_off,ibzxxc_off,ibxxyyc_off,ibyzff_off,ibzxxf_off,ibxxyyf_off,& 
!!                              ibyyzzr_off /)
!!  array_of_types         = (//)
!!
!!
!!
!!
!!end subroutine


!determine a set of localisation regions from the centers and the radii.
!cut in cubes the global reference system
subroutine determine_locreg_parallel(iproc,nproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,orbs,calculateBounds)!,outofzone)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
  type(orbitals_data),intent(in) :: orbs
  logical,dimension(nlr),intent(in):: calculateBounds
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg_parallel'
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz
  logical :: warningx,warningy,warningz,calc
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez,iorb
  integer :: ln1,ln2,ln3
  integer :: iilr,ierr
  integer,dimension(3) :: outofzone
  integer,dimension(0:nproc-1) :: nlr_par,islr_par
  real(gp) :: rx,ry,rz,cutoff

  !!if (iproc == 0) then
  !!   write(*,*)'Inside determine_locreg_periodic:'
  !!end if

!  call parallel_repartition_locreg(iproc,nproc,nlr,nlr_par,islr_par)

  !initialize out of zone and logicals
  outofzone (:) = 0
  warningx = .false.
  warningy = .false.
  warningz = .false.

  !determine the limits of the different localisation regions
  do ilr=1,nlr
     !nullify all pointers
     nullify(Llr(ilr)%projflg)
     nullify(Llr(ilr)%wfd%keygloc)
     nullify(Llr(ilr)%wfd%keyglob)
     nullify(Llr(ilr)%wfd%keyv)
     nullify(Llr(ilr)%bounds%ibyyzz_r)
     nullify(Llr(ilr)%bounds%kb%ibyz_c)
     nullify(Llr(ilr)%bounds%kb%ibxz_c)
     nullify(Llr(ilr)%bounds%kb%ibxy_c)
     nullify(Llr(ilr)%bounds%kb%ibyz_f)
     nullify(Llr(ilr)%bounds%kb%ibxz_f)
     nullify(Llr(ilr)%bounds%kb%ibxy_f)
     nullify(Llr(ilr)%bounds%sb%ibzzx_c)
     nullify(Llr(ilr)%bounds%sb%ibyyzz_c)
     nullify(Llr(ilr)%bounds%sb%ibxy_ff)
     nullify(Llr(ilr)%bounds%sb%ibzzx_f)
     nullify(Llr(ilr)%bounds%sb%ibyyzz_f)
     nullify(Llr(ilr)%bounds%gb%ibzxx_c)
     nullify(Llr(ilr)%bounds%gb%ibxxyy_c)
     nullify(Llr(ilr)%bounds%gb%ibyz_ff)
     nullify(Llr(ilr)%bounds%gb%ibzxx_f)
     nullify(Llr(ilr)%bounds%gb%ibxxyy_f)

     calc=.false.
     do iorb=1,orbs%norbp
        if(ilr == orbs%inwhichLocreg(iorb+orbs%isorb)) calc=.true.
     end do
     if (.not. calc) cycle         !calculate only for the locreg on this processor, without repeating for same locreg

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
!!!     if (iproc == 0 .and. verbose > 1) then
!!!        if ((iex - isx >= Glr%d%n1 - 14) .and. (warningx .eqv. .false.)) then
!!!           write(*,*)'Width of direction x :',(iex - isx)*hx,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n1*hx
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the '
!!!           write(*,*)'simulation box width. This is the only warning for x direction.'
!!!           warningx = .true.
!!!        end if
!!!        if ((iey - isy >= Glr%d%n2 - 14) .and. (warningy .eqv. .false.)) then
!!!           write(*,*)'Width of direction y :',(iey - isy)*hy,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n2*hy,'.'
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the width'
!!!           write(*,*)'of the simulation box. This is the only warning for y direction.'
!!!           warningy = .true.
!!!        end if
!!!        if ((iez - isz >= Glr%d%n3 - 14) .and. (warningz .eqv. .false.)) then
!!!           write(*,*)'Width of direction z :',(iez - isz)*hz,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n3*hz,'.'
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the width'
!!!           write(*,*)'of the simulation box. This is the only warning for z direction.'
!!!           warningz = .true.
!!!        end if 
!!!     end if

     ! Localization regions should always have free boundary conditions
     Llr(ilr)%geocode='F'

     !assign the starting/ending points and outofzone for the different
     ! geometries
     select case(Glr%geocode)
     case('F')
        isx=max(isx,Glr%ns1)
        isy=max(isy,Glr%ns2)
        isz=max(isz,Glr%ns3)

        iex=min(iex,Glr%ns1+Glr%d%n1)
        iey=min(iey,Glr%ns2+Glr%d%n2)
        iez=min(iez,Glr%ns3+Glr%d%n3)

     case('S')
        ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if
        end if

        ! Get starting and ending for y direction (perpendicular to surface)
        isy=max(isy,Glr%ns2)
        iey=min(iey,Glr%ns2 + Glr%d%n2)
        outofzone(2) = 0

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if
        end if

     case('P')
         ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if
        end if

        ! Get starting and ending for y direction (perpendicular to surface)
        if (iey - isy >= Glr%d%n2) then
           isy=Glr%ns2
           iey=Glr%ns2 + Glr%d%n2
         else
           isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
           iey= ln2 + isy
           if (iey > Glr%ns2+Glr%d%n2) then
              outofzone(2)=modulo(iey,Glr%d%n2+1)
           end if
        end if

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if
        end if
     end select

     !values for the starting point of the cube for wavelet grid
     Llr(ilr)%ns1=isx
     Llr(ilr)%ns2=isy
     Llr(ilr)%ns3=isz

     !dimensions of the localisation region
     Llr(ilr)%d%n1=iex-isx
     Llr(ilr)%d%n2=iey-isy
     Llr(ilr)%d%n3=iez-isz

     !assign outofzone
     Llr(ilr)%outofzone(:) = outofzone(:)

     ! Set the conditions for ext_buffers (conditions for buffer size)
     Gperx=(Glr%geocode /= 'F')
     Gpery=(Glr%geocode == 'P')
     Gperz=(Glr%geocode /= 'F')
     Lperx=(Llr(ilr)%geocode /= 'F')
     Lpery=(Llr(ilr)%geocode == 'P')
     Lperz=(Llr(ilr)%geocode /= 'F')

     !calculate the size of the buffers of interpolating function grid
     call ext_buffers(Gperx,Gnbl1,Gnbr1)
     call ext_buffers(Gpery,Gnbl2,Gnbr2)
     call ext_buffers(Gperz,Gnbl3,Gnbr3)
     call ext_buffers(Lperx,Lnbl1,Lnbr1)
     call ext_buffers(Lpery,Lnbl2,Lnbr2)
     call ext_buffers(Lperz,Lnbl3,Lnbr3)

     !starting point of the region for interpolating functions grid
     Llr(ilr)%nsi1= 2 * Llr(ilr)%ns1 - (Lnbl1 - Gnbl1)
     Llr(ilr)%nsi2= 2 * Llr(ilr)%ns2 - (Lnbl2 - Gnbl2)
     Llr(ilr)%nsi3= 2 * Llr(ilr)%ns3 - (Lnbl3 - Gnbl3)

     !dimensions of the fine grid inside the localisation region
     Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
     Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
     Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz

     !NOTE: This will not work with symmetries (must change it)
     Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
     Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
     Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz

     !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
     Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
     Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
     Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31

!DEBUG
!     if (iproc == 0) then
!        write(*,*)'Description of zone:',ilr
!        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
!        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
!        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
!        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
!        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
!        write(*,*)'outofzone',ilr,':',outofzone(:)
!     end if
!DEBUG

    ! construct the wavefunction descriptors (wfd)
     call determine_wfd_periodicity(ilr,nlr,Glr,Llr)

     ! Sould check if nfu works properly... also relative to locreg!!
     !if the localisation region is isolated build also the bounds
     if (Llr(ilr)%geocode=='F') then
        ! Check whether the bounds shall be calculated. Do this only if the currect process handles
        ! orbitals in the current localization region.
        if(calculateBounds(ilr)) then
           call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
               Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
               Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
        end if
     end if
  end do !on iilr

!  call make_LLr_MpiType(Llr,nlr,mpiLlr)

!  call MPI_ALLREDUCE(Llr(1),Llr(1),nlr,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
  !after all localisation regions are determined draw them
  !call draw_locregs(nlr,hx,hy,hz,Llr)

END SUBROUTINE determine_locreg_parallel







! Count for each orbital and each process the number of overlapping orbitals.
subroutine determine_overlap_from_descriptors(iproc, nproc, orbs, lzd, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(out):: op
type(p2pComms),intent(out):: comon

! Local variables
integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
integer:: js1, je1, js2, je2, js3, je3, iiorb, istat, iall, noverlaps, ierr
logical:: ovrlpx, ovrlpy, ovrlpz
integer:: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp, i1, i2, jjorb, ii
logical,dimension(:,:,:),allocatable:: overlapMatrix
integer,dimension(:),allocatable:: noverlapsarr, displs, recvcnts, overlaps_comon
integer,dimension(:,:),allocatable:: overlaps_op
character(len=*),parameter:: subname='determine_overlap_from_descriptors'

allocate(overlapMatrix(orbs%norb,maxval(orbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, overlapMatrix, 'overlapMatrix', subname)
!!overlapMatrix=.false.
allocate(noverlapsarr(orbs%norbp), stat=istat)
call memocc(istat, noverlapsarr, 'noverlapsarr', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)
allocate(recvcnts(0:nproc-1), stat=istat)
call memocc(istat, recvcnts, 'recvcnts', subname)


overlapMatrix=.false.
    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
    ilrold=-1
    do iorb=1,orbs%norbp
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=orbs%isorb+iorb
        ilr=orbs%inWhichLocreg(iiorb)
        call get_start_and_end_indices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
            jlr=orbs%inWhichLocreg(jorb)
            call get_start_and_end_indices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
                ! Now explicitely check whether there is an overlap by using the descriptors.
                call overlapbox_from_descriptors(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
                     lzd%llr(jlr)%d%n1, lzd%llr(jlr)%d%n2, lzd%llr(jlr)%d%n3, &
                     lzd%glr%d%n1, lzd%glr%d%n2, lzd%glr%d%n3, &
                     lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
                     lzd%llr(jlr)%ns1, lzd%llr(jlr)%ns2, lzd%llr(jlr)%ns3, &
                     lzd%glr%ns1, lzd%glr%ns2, lzd%glr%ns3, &
                     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c, &
                     lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyv, lzd%llr(jlr)%wfd%keygloc, lzd%llr(jlr)%wfd%keyv, &
                     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)
                if(n1_ovrlp>0 .and. n2_ovrlp>0 .and. n3_ovrlp>0) then
                    ! There is really an overlap
                    overlapMatrix(jorb,iorb,iproc)=.true.
                    ioverlaporb=ioverlaporb+1
                    if(ilr/=ilrold) then
                        ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                        ! would get the same orbitals again. Therefore the counter is not increased
                        ! in that case.
                        ioverlapMPI=ioverlapMPI+1
                    end if
                else
                    overlapMatrix(jorb,iorb,iproc)=.false.
                end if
            else
                overlapMatrix(jorb,iorb,iproc)=.false.
            end if
        end do
        noverlapsarr(iorb)=ioverlaporb
        ilrold=ilr
    end do
    !comon%noverlaps(jproc)=ioverlapMPI
    noverlaps=ioverlapMPI

!call mpi_allreduce(overlapMatrix, orbs%norb*maxval(orbs%norb_par(:,0))*nproc, mpi_sum mpi_comm_world, ierr)

! Communicate op%noverlaps and comon%noverlaps

    if (nproc > 1) then
       call mpi_allgatherv(noverlapsarr, orbs%norbp, mpi_integer, op%noverlaps, orbs%norb_par, &
            orbs%isorb_par, mpi_integer, mpi_comm_world, ierr)
    else
       call vcopy(orbs%norb,noverlapsarr(1),1,op%noverlaps(1),1)
    end if
    do jproc=0,nproc-1
       recvcnts(jproc)=1
       displs(jproc)=jproc
    end do
    if (nproc > 1) then
       call mpi_allgatherv(noverlaps, 1, mpi_integer, comon%noverlaps, recvcnts, &
            displs, mpi_integer, mpi_comm_world, ierr)
    else
       comon%noverlaps=noverlaps
    end if



allocate(op%overlaps(maxval(op%noverlaps),orbs%norb), stat=istat)
call memocc(istat, op%overlaps, 'op%overlaps', subname)
!!allocate(comon%overlaps(maxval(comon%noverlaps),0:nproc-1), stat=istat)
!!call memocc(istat, comon%overlaps, 'comon%overlaps', subname)

allocate(overlaps_op(maxval(op%noverlaps),orbs%norbp), stat=istat)
call memocc(istat, overlaps_op, 'overlaps_op', subname)
allocate(overlaps_comon(comon%noverlaps(iproc)), stat=istat)
call memocc(istat, overlaps_comon, 'overlaps_comon', subname)


! Now we know how many overlaps have to be calculated, so determine which orbital overlaps
! with which one. This is essentially the same loop as above, but we use the array 'overlapMatrix'
! which indicates the overlaps.

! Initialize to some value which will never be used.
op%overlaps=-1
!!comon%overlaps=-1

iiorb=0
ioverlapMPI=0 ! counts the overlaps for the given MPI process.
ilrold=-1
do iorb=1,orbs%norbp
    ioverlaporb=0 ! counts the overlaps for the given orbital.
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    do jorb=1,orbs%norb
        jlr=orbs%inWhichLocreg(jorb)
        if(overlapMatrix(jorb,iorb,iproc)) then
            ioverlaporb=ioverlaporb+1
            !op%overlaps(ioverlaporb,iiorb)=jorb
            overlaps_op(ioverlaporb,iorb)=jorb
            if(ilr/=ilrold) then
                ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                ! would get the same orbitals again. Therefore the counter is not increased
                ! in that case.
                ioverlapMPI=ioverlapMPI+1
                !comon%overlaps(ioverlapMPI,iproc)=jorb
                overlaps_comon(ioverlapMPI)=jorb
            end if
        end if
    end do 
    ilrold=ilr
end do


displs(0)=0
recvcnts(0)=comon%noverlaps(0)
do jproc=1,nproc-1
    recvcnts(jproc)=comon%noverlaps(jproc)
    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
end do
!!if (nproc > 1) then
!!   call mpi_allgatherv(overlaps_comon, comon%noverlaps(iproc), mpi_integer, comon%overlaps, recvcnts, &
!!        displs, mpi_integer, mpi_comm_world, ierr)
!!else
!!   call vcopy(comon%noverlaps(iproc),overlaps_comon(1),1,comon%overlaps(1,0),1)
!!end if
ii=maxval(op%noverlaps)
displs(0)=0
recvcnts(0)=ii*orbs%norb_par(0,0)
do jproc=1,nproc-1
    recvcnts(jproc)=ii*orbs%norb_par(jproc,0)
    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
end do
if (nproc > 1) then
   call mpi_allgatherv(overlaps_op, ii*orbs%norbp, mpi_integer, op%overlaps, recvcnts, &
        displs, mpi_integer, mpi_comm_world, ierr)
else
   call vcopy(ii*orbs%norbp,overlaps_op(1,1),1,op%overlaps(1,1),1)
end if

do iorb=1,orbs%norb
end do


iall=-product(shape(overlapMatrix))*kind(overlapMatrix)
deallocate(overlapMatrix, stat=istat)
call memocc(istat, iall, 'overlapMatrix', subname)




op%noverlapsmaxp=maxval(op%noverlaps(orbs%isorb+1:orbs%isorb+orbs%norbp))

!allocate(op%olr(op%noverlapsmaxp,orbs%norbp), stat=istat)
!do i2=1,orbs%norbp
!    do i1=1,op%noverlapsmaxp
!        call nullify_locreg_descriptors(op%olr(i1,i2))
!    end do
!end do


!do iorb=1,orbs%norbp
!    iiorb=orbs%isorb_par(iproc)+iorb
!    ilr=orbs%inWhichLocreg(iiorb)
!    do jorb=1,op%noverlaps(iiorb)
!        jjorb=op%overlaps(jorb,iiorb)
!        jlr=orbs%inWhichLocreg(jjorb)
!    end do
!end do

iall=-product(shape(noverlapsarr))*kind(noverlapsarr)
deallocate(noverlapsarr, stat=istat)
call memocc(istat, iall, 'noverlapsarr', subname)

iall=-product(shape(overlaps_op))*kind(overlaps_op)
deallocate(overlaps_op, stat=istat)
call memocc(istat, iall, 'overlaps_op', subname)

iall=-product(shape(overlaps_comon))*kind(overlaps_comon)
deallocate(overlaps_comon, stat=istat)
call memocc(istat, iall, 'overlaps_comon', subname)

iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)

iall=-product(shape(recvcnts))*kind(recvcnts)
deallocate(recvcnts, stat=istat)
call memocc(istat, iall, 'recvcnts', subname)

end subroutine determine_overlap_from_descriptors



subroutine determine_overlapdescriptors_from_descriptors(llr_i, llr_j, glr, olr)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: llr_i, llr_j, glr
type(locreg_descriptors),intent(out):: olr

! Local variables
integer:: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp
character(len=*),parameter:: subname='determine_overlapdescriptors_from_descriptors'



! Determine the values describing the localization region.
! First coarse region.
call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, &
     llr_i%wfd%keygloc, llr_i%wfd%keyv, llr_j%wfd%keygloc, llr_j%wfd%keyv, &
     olr%d%n1, olr%d%n2, olr%d%n3, olr%ns1, olr%ns2, olr%ns3, olr%wfd%nseg_c)

! Now the fine region.
call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     llr_i%wfd%nseg_f, llr_j%wfd%nseg_f, &
     llr_i%wfd%keygloc(1,llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), llr_i%wfd%keyv(llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), &
     llr_j%wfd%keygloc(1,llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), llr_j%wfd%keyv(llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), &
     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, olr%wfd%nseg_f)

     ! Determine the boundary for the fine part.
     ! ns1_ovrlp etc is in global coordinates, but olr%d%nfl1 etc is in local coordinates, so correct this.
     olr%d%nfl1=ns1_ovrlp-olr%ns1
     olr%d%nfu1=olr%d%nfl1+n1_ovrlp
     olr%d%nfl2=ns2_ovrlp-olr%ns2
     olr%d%nfu2=olr%d%nfl2+n2_ovrlp
     olr%d%nfl3=ns2_ovrlp-olr%ns2
     olr%d%nfu3=olr%d%nfl3+n3_ovrlp

if(llr_i%geocode/=llr_j%geocode) then
    write(*,*) 'ERROR: llr_i%geocode/=llr_j%geocode'
    stop
end if
olr%geocode=llr_i%geocode

! Dimensions for interpolating scaling function grid
olr%d%n1i=2*olr%d%n1+31
olr%d%n2i=2*olr%d%n2+31
olr%d%n3i=2*olr%d%n3+31



! Allocate the descriptor structures
call allocate_wfd(olr%wfd,subname)


! some checks
call check_overlapregion(glr, llr_i, llr_j, olr)



! Fill the descriptors for the coarse part.
call overlapdescriptors_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     olr%d%n1, olr%d%n2, olr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     olr%ns1, olr%ns2, olr%ns3, &
     llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, olr%wfd%nseg_c, &
     llr_i%wfd%keygloc, llr_i%wfd%keyv, llr_j%wfd%keygloc, llr_j%wfd%keyv, &
     olr%wfd%keygloc, olr%wfd%keyv, &
     olr%wfd%nvctr_c)

! Fill the descriptors for the fine part.
call overlapdescriptors_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     olr%d%n1, olr%d%n2, olr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     olr%ns1, olr%ns2, olr%ns3, &
     llr_i%wfd%nseg_f, llr_j%wfd%nseg_f, olr%wfd%nseg_f, &
     llr_i%wfd%keygloc(1,llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), llr_i%wfd%keyv(llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), &
     llr_j%wfd%keygloc(1,llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), llr_j%wfd%keyv(llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), &
     olr%wfd%keygloc(1,olr%wfd%nseg_c+min(1,olr%wfd%nseg_f)), olr%wfd%keyv(olr%wfd%nseg_c+min(1,olr%wfd%nseg_f)), &
     olr%wfd%nvctr_f)



end subroutine determine_overlapdescriptors_from_descriptors





subroutine check_overlapregion(glr, llr_i, llr_j, olr)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: glr, llr_i, llr_j, olr

  if(olr%ns1<glr%ns1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1 = ', olr%ns1, ' < ', glr%ns1, '= glr%ns1'
      stop
  end if
  if(olr%ns2<glr%ns2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2 = ', olr%ns2, ' < ', glr%ns2, '= glr%ns2'
      stop
  end if
  if(olr%ns3<glr%ns3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3 = ', olr%ns3, ' < ', glr%ns3, '= glr%ns3'
      stop
  end if
  if(olr%ns1+olr%d%n1>glr%ns1+glr%d%n1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1+olr%d%n1 = ', olr%ns1+olr%d%n1, ' < ', glr%ns1+glr%d%n1, '= glr%ns1+glr%d%n1'
      stop
  end if
  if(olr%ns2+olr%d%n2>glr%ns2+glr%d%n2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2+olr%d%n2 = ', olr%ns2+olr%d%n2, ' < ', glr%ns2+glr%d%n2, '= glr%ns2+glr%d%n2'
      stop
  end if
  if(olr%ns3+olr%d%n3>glr%ns3+glr%d%n3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3+olr%d%n3 = ', olr%ns3+olr%d%n3, ' < ', glr%ns3+glr%d%n3, '= glr%ns3+glr%d%n3'
      stop
  end if
  
  if(olr%ns1<llr_i%ns1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1 = ', olr%ns1, ' < ', llr_i%ns1, '= llr_i%ns1'
      stop
  end if
  if(olr%ns2<llr_i%ns2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2 = ', olr%ns2, ' < ', llr_i%ns2, '= llr_i%ns2'
      stop
  end if
  if(olr%ns3<llr_i%ns3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3 = ', olr%ns3, ' < ', llr_i%ns3, '= llr_i%ns3'
      stop
  end if
  if(olr%ns1+olr%d%n1>llr_i%ns1+llr_i%d%n1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1+olr%d%n1 = ', olr%ns1+olr%d%n1, ' < ', llr_i%ns1+llr_i%d%n1, '= llr_i%ns1+llr_i%d%n1'
      stop
  end if
  if(olr%ns2+olr%d%n2>llr_i%ns2+llr_i%d%n2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2+olr%d%n2 = ', olr%ns2+olr%d%n2, ' < ', llr_i%ns2+llr_i%d%n2, '= llr_i%ns2+llr_i%d%n2'
      stop
  end if
  if(olr%ns3+olr%d%n3>llr_i%ns3+llr_i%d%n3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3+olr%d%n3 = ', olr%ns3+olr%d%n3, ' < ', llr_i%ns3+llr_i%d%n3, '= llr_i%ns3+llr_i%d%n3'
      stop
  end if
  
  if(olr%ns1<llr_j%ns1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1 = ', olr%ns1, ' < ', llr_j%ns1, '= llr_j%ns1'
      stop
  end if
  if(olr%ns2<llr_j%ns2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2 = ', olr%ns2, ' < ', llr_j%ns2, '= llr_j%ns2'
      stop
  end if
  if(olr%ns3<llr_j%ns3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3 = ', olr%ns3, ' < ', llr_j%ns3, '= llr_j%ns3'
      stop
  end if
  if(olr%ns1+olr%d%n1>llr_j%ns1+llr_j%d%n1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1+olr%d%n1 = ', olr%ns1+olr%d%n1, ' < ', llr_j%ns1+llr_j%d%n1, '= llr_j%ns1+llr_j%d%n1'
      stop
  end if
  if(olr%ns2+olr%d%n2>llr_j%ns2+llr_j%d%n2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2+olr%d%n2 = ', olr%ns2+olr%d%n2, ' < ', llr_j%ns2+llr_j%d%n2, '= llr_j%ns2+llr_j%d%n2'
      stop
  end if
  if(olr%ns3+olr%d%n3>llr_j%ns3+llr_j%d%n3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3+olr%d%n3 = ', olr%ns3+olr%d%n3, ' < ', llr_j%ns3+llr_j%d%n3, '= llr_j%ns3+llr_j%d%n3'
      stop
  end if


end subroutine check_overlapregion




!> Returns the starting and ending indices (on the coarse grid) of a given localization region.
!  Attention: There is a similar routine that uses "is1=lr%ns1+1" instead of "is1=lr%ns1".
!  Check why there is a difference of 1.
subroutine get_start_and_end_indices(lr, is1, ie1, is2, ie2, is3, ie3)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: lr
integer,intent(out):: is1, ie1, is2, ie2, is3, ie3

  is1=lr%ns1
  ie1=lr%ns1+lr%d%n1
  is2=lr%ns2
  ie2=lr%ns2+lr%d%n2
  is3=lr%ns3
  ie3=lr%ns3+lr%d%n3

end subroutine get_start_and_end_indices





!> Gives the dimensions of the overlap box resulting from the overlap of two wavefunction descriptors and
!> the number of segments of the resulting overlap descriptor.
!> Calling arguments: *_i refers to overlap region i (input)
!>                    *_j refers to overlap region j (input)
!>                    *_g refers to the global region (input)
!>                    *_k refers to the overlap region (output)
subroutine overlapbox_from_descriptors(n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, &
           ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, &
           nseg_i, nseg_j, &
           keyg_i, keyv_i, keyg_j, keyv_j, &
           n1_k, n2_k, n3_k, ns1_k, ns2_k, ns3_k, nseg_k)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g
integer,intent(in):: ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g
integer:: nseg_i, nseg_j
integer,dimension(2,nseg_i),intent(in):: keyg_i
integer,dimension(nseg_i),intent(in):: keyv_i
integer,dimension(2,nseg_j),intent(in):: keyg_j
integer,dimension(nseg_j),intent(in):: keyv_j
integer,intent(out):: n1_k, n2_k, n3_k, ns1_k, ns2_k, ns3_k, nseg_k

! Local variables
integer:: iseg, jseg, knvctr, istart, jstart, kstart, istartg, jstartg, kstartg
integer:: iend, jend, kend, iendg, jendg, kendg, transform_index
integer:: kxs, kys, kzs, kxe, kye, kze, kxemax, kyemax, kzemax
character(len=1):: increase


! Initialize the return values such that they represent a box with no volume
ns1_k=ns1_g+n1_g+1
ns2_k=ns2_g+n2_g+1
ns3_k=ns3_g+n3_g+1
n1_k=0
n2_k=0
n3_k=0
nseg_k=0

! Initialize some counters
iseg=1
jseg=1
kxemax=0
kyemax=0
kzemax=0


segment_loop: do

    ! Starting point in local coordinates
    istart=keyg_i(1,iseg)
    jstart=keyg_j(1,jseg)

    ! Get the global counterparts
    istartg=transform_index(istart, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
    jstartg=transform_index(jstart, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

    ! Ending point in local coordinates
    iend=keyg_i(2,iseg)
    jend=keyg_j(2,jseg)

    ! Get the global counterparts
    iendg=transform_index(iend, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
    jendg=transform_index(jend, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

    ! Determine starting and ending point of the common segment in global coordinates.
    kstartg=max(istartg,jstartg)
    kendg=min(iendg,jendg)

    ! Determine which segment counter should be increased.
    if((iendg<=jendg .and. iseg<nseg_i) .or. jseg==nseg_j) then
        increase='i'
    else if(jseg<nseg_j) then
        increase='j'
    end if

    ! Check whether this common segment has a non-zero length
    if(kendg-kstartg+1>0) then
        nseg_k=nseg_k+1

        ! Get the global coordinates of this segment
        call get_coordinates(kstartg, n1_g, n2_g, n3_g, kxs, kys, kzs)
        call get_coordinates(kendg, n1_g, n2_g, n3_g, kxe, kye, kze)

        ! Check whether this segment enlarges the overlap box
        if(kxs<ns1_k) ns1_k=kxs
        if(kys<ns2_k) ns2_k=kys
        if(kzs<ns3_k) ns3_k=kzs
        if(kxe>kxemax) kxemax=kxe
        if(kye>kyemax) kyemax=kye
        if(kze>kzemax) kzemax=kze

    end if

    ! Check whether all segments of both localization regions have been processed.
    if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

    ! Increase the segment index
    if(increase=='i') then
        iseg=iseg+1
    else if(increase=='j') then
        jseg=jseg+1
    end if

end do segment_loop

! n1_k etc is the length of the segment, but kxemax etc is the end position of the segment, 
! therefore subtract the starting position
n1_k=kxemax-ns1_k
n2_k=kyemax-ns2_k
n3_k=kzemax-ns3_k



end subroutine overlapbox_from_descriptors




!> Creates the overlap descriptor from two input overlap descriptors. The dimension of the overlap box (i.e.
!> ns1_k,ns2_k,ns3_k (starting indices) and n1_k,n2_k,n3_k (length) have to be determined earlier (using
!> the subroutine overlapbox_from_descriptors)).
!> Calling arguments: *_i refers to overlap region i (input)
!>                    *_j refers to overlap region j (input)
!>                    *_g refers to the global region (input)
!>                    *_k refers to the overlap region (input/output)
subroutine overlapdescriptors_from_descriptors(n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, &
           ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, ns1_k, ns2_k, ns3_k, &
           nseg_i, nseg_j, nseg_k, &
           keyg_i, keyv_i, keyg_j, keyv_j, keyg_k, keyv_k, nvctr_k)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k
integer,intent(in):: ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, ns1_k, ns2_k, ns3_k
integer,intent(in):: nseg_i, nseg_j, nseg_k
integer,dimension(2,nseg_i),intent(in):: keyg_i
integer,dimension(nseg_i),intent(in):: keyv_i
integer,dimension(2,nseg_j),intent(in):: keyg_j
integer,dimension(nseg_j),intent(in):: keyv_j
integer,dimension(2,nseg_k),intent(out):: keyg_k
integer,dimension(nseg_k),intent(out):: keyv_k
integer,intent(out):: nvctr_k

! Local variables
integer:: iseg, jseg, kseg, knvctr, istart, jstart, kstart, istartg, jstartg, kstartg
integer:: iend, jend, kend, iendg, jendg, kendg, transform_index
character(len=1):: increase

! Initialize some counters
iseg=1
jseg=1
kseg=0
knvctr=0
nvctr_k=0

segment_loop: do

    ! Starting point in local coordinates
    istart=keyg_i(1,iseg)
    jstart=keyg_j(1,jseg)

    ! Get the global counterparts
    istartg=transform_index(istart, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
    jstartg=transform_index(jstart, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

    ! Ending point in local coordinates
    iend=keyg_i(2,iseg)
    jend=keyg_j(2,jseg)

    ! Get the global counterparts
    iendg=transform_index(iend, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
    jendg=transform_index(jend, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

    ! Determine starting and ending point of the common segment in global coordinates.
    kstartg=max(istartg,jstartg)
    kendg=min(iendg,jendg)
    if((iendg<=jendg .and. iseg<nseg_i) .or. jseg==nseg_j) then
        increase='i'
    else if(jseg<nseg_j) then
        increase='j'
    end if

    ! Check whether this common segment has a non-zero length
    if(kendg-kstartg+1>0) then
        kseg=kseg+1

        ! Transform the starting and ending point to the overlap localization region.
        kstart=transform_index(kstartg, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, ns1_g-ns1_k, ns2_g-ns2_k, ns3_g-ns3_k)
        kend=transform_index(kendg, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, ns1_g-ns1_k, ns2_g-ns2_k, ns3_g-ns3_k)

        ! Assign the values to the descriptors
        keyg_k(1,kseg)=kstart
        keyg_k(2,kseg)=kend
        keyv_k(kseg)=knvctr+1

        knvctr=knvctr+kendg-kstartg+1
    end if

    ! Check whether all segments of both localization regions have been processed
    if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

    ! Increase the segment index
    if(increase=='i') then
        iseg=iseg+1
    else if(increase=='j') then
        jseg=jseg+1
    end if

end do segment_loop

nvctr_k=knvctr


end subroutine overlapdescriptors_from_descriptors




!> Transform an index from localization region A to localization region B.
! Calling arguments:
!   ist: index to transform (with respcet to coordinate system of locreg A)
!   n1a, n2a, n3a:             box size of locreg A
!   n1b, n2b, n3b:             box size of locreg B
!   nshift1, nshift2, nshift3: nsa-nsb, where nsa,nsb are the starting points of the boxes of A,B (for all 3 dimensions)
function transform_index(ist, n1a, n2a, n3a, n1b, n2b, n3b, nshift1, nshift2, nshift3)
implicit none

! Calling arguments
integer,intent(in):: ist, n1a, n2a, n3a, n1b, n2b, n3b, nshift1, nshift2, nshift3
integer:: transform_index

! Local variables
integer:: ii, ix, iy, iz, ixg, iyg, izg, istg

  ! Get the coordinates with respect to localization region A
  ii = ist - 1
  iz = ii / ((n1a+1) * (n2a+1))
  ii = ii - iz * ((n1a+1) * (n2a+1))
  iy = ii / (n1a+1)
  ix = ii - iy * (n1a+1)

  ! Transform ix, iy, iz to the coordinates with respect to localization region B
  izg = iz + nshift3
  iyg = iy + nshift2
  ixg = ix + nshift1

  ! Transform ist to its counterpart in the coordinate system of B
  istg = izg*(n1b+1)*(n2b+1) + iyg*(n1b+1) + ixg + 1
  
  ! Assign istg to the value that is passed back.
  transform_index=istg

end function transform_index





!> Get the coordinates of ist with respect to its localization region
! Calling arguments:
!   ist          index for which coordinates shall be calculated
!   n1, n2, n3   box sizes
!   ix, iy, iz   coordinates of ist
subroutine get_coordinates(ist, n1, n2, n3, ix, iy, iz)
implicit none

! Calling arguments
integer,intent(in):: ist, n1, n2, n3
integer,intent(out):: ix, iy, iz

! Local variable
integer:: ii

  ! Get the coordinates ix, iy, iz
  ii = ist - 1
  !ii = ist
  iz = ii / ((n1+1) * (n2+1))
  ii = ii - iz * ((n1+1) * (n2+1))
  iy = ii / (n1+1)
  ix = ii - iy * (n1+1)

end subroutine get_coordinates

function check_whether_locregs_overlap(llr_i, llr_j, glr)
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: llr_i, llr_j, glr
logical:: check_whether_locregs_overlap

! Local variables
integer:: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp

  ! Check whether there is an overlap by comparing the descriptors.
  call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
       llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
       glr%d%n1, glr%d%n2, glr%d%n3, &
       llr_i%ns1, llr_i%ns2, llr_i%ns3, &
       llr_j%ns1, llr_j%ns2, llr_j%ns3, &
       glr%ns1, glr%ns2, glr%ns3, &
       llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, &
       llr_i%wfd%keygloc, llr_i%wfd%keyv, llr_j%wfd%keygloc, llr_j%wfd%keyv, &
       n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)

  ! n1_ovrlp, n2_ovrlp, n3_ovrlp are the dimensions of the overlap localization regions.
  if(n1_ovrlp>0 .and. n2_ovrlp>0 .and. n3_ovrlp>0) then
      ! There is an overlap
      check_whether_locregs_overlap=.true.
  else
      ! There is no overlap
      check_whether_locregs_overlap=.false.
  end if
end function check_whether_locregs_overlap

subroutine check_overlap(Llr_i, Llr_j, Glr, overlap)
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: Llr_i, Llr_j, Glr
logical, intent(out) :: overlap

! Local variables
integer:: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp
logical:: go1, go2, go3
  
  ! Begin by checking if the boxes overlap
  overlap = .false.
  go1 = (Llr_i%ns1 < (Llr_j%ns1 + Llr_j%d%n1)) .and. ((Llr_i%ns1 + Llr_i%d%n1) > Llr_j%ns1 )
  go2 = (Llr_i%ns2 < (Llr_j%ns2 + Llr_j%d%n2)) .and. ((Llr_i%ns2 + Llr_i%d%n2) > Llr_j%ns2 )
  go3 = (Llr_i%ns3 < (Llr_j%ns3 + Llr_j%d%n3)) .and. ((Llr_i%ns3 + Llr_i%d%n3) > Llr_j%ns3 )
  if(go1 .and. go2 .and. go3) overlap = .true.

  ! If so, check if the descriptors overlap
  if (overlap) then
     ! Check whether there is an overlap by comparing the descriptors.
     call overlapbox_from_descriptors(Llr_i%d%n1, Llr_i%d%n2, Llr_i%d%n3, &
          Llr_j%d%n1, Llr_j%d%n2, Llr_j%d%n3, &
          Glr%d%n1, Glr%d%n2, Glr%d%n3, &
          Llr_i%ns1, Llr_i%ns2, Llr_i%ns3, &
          Llr_j%ns1, Llr_j%ns2, Llr_j%ns3, &
          Glr%ns1, Glr%ns2, Glr%ns3, &
          Llr_i%wfd%nseg_c, Llr_j%wfd%nseg_c, &
          Llr_i%wfd%keygloc, Llr_i%wfd%keyv, Llr_j%wfd%keygloc, Llr_j%wfd%keyv, &
          n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)

     ! n1_ovrlp, n2_ovrlp, n3_ovrlp are the dimensions of the overlap localization regions.
     if(n1_ovrlp>0 .and. n2_ovrlp>0 .and. n3_ovrlp>0) then
         ! There is an overlap
         overlap=.true.
     else
         ! There is no overlap
         overlap=.false.
     end if
  end if

end subroutine check_overlap

subroutine transform_keyglob_to_keygloc(Glr,Llr,nseg,keyglob,keygloc)
use module_base
use module_types
use module_interfaces
implicit none
type(locreg_descriptors),intent(in):: Glr, Llr
integer, intent(in) :: nseg
integer, dimension(2,nseg),intent(in) :: keyglob
integer, dimension(2,nseg),intent(out) :: keygloc
!local variables
integer :: i, j, j0, ii, iz, iy, ix
do i = 1 , 2
   do j = 1, nseg
      ! Writing keyglob in cartesian coordinates
      j0 = keyglob(i,j)
      ii = j0-1
      iz = ii/((Glr%d%n1+1)*(Glr%d%n2+1))
      ii = ii-iz*(Glr%d%n1+1)*(Glr%d%n2+1)
      iy = ii/(Glr%d%n1+1)
      ix = ii-iy*(Glr%d%n1+1)

      ! Checking consistency
      if(iz < Llr%ns3 .or. iy < Llr%ns2 .or. ix < Llr%ns1) stop 'transform_keyglob_to_keygloc : minimum overflow'
      if(iz > Llr%ns3+Llr%d%n3 .or. iy > Llr%ns2+Llr%d%n2 .or. ix > Llr%ns1+Llr%d%n1)&
         stop 'transform_keyglob_to_keygloc : maximum overflow'

      ! Using coordinates to write keygloc      
      keygloc(i,j) = (iz-Llr%ns3)*(Llr%d%n1+1)*(Llr%d%n2+1) + (iy-Llr%ns2)*(Llr%d%n1+1) + (ix-Llr%ns1) + 1
   end do
end do

end subroutine transform_keyglob_to_keygloc
