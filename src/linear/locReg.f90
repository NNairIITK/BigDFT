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
          Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),nseg_c,nvctr_c,Llr(ilr)%outofzone(:))
   !fine part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),nseg_f,nvctr_f,Llr(ilr)%outofzone(:))

   ! Assign the values to Llr
   Llr(ilr)%wfd%nseg_c = nseg_c
   Llr(ilr)%wfd%nseg_f = nseg_f
   Llr(ilr)%wfd%nvctr_c= nvctr_c
   Llr(ilr)%wfd%nvctr_f= nvctr_f

   !allocate the wavefunction descriptors following the needs
   call allocate_wfd(Llr(ilr)%wfd,subname)

   !Now, fill the descriptors:
   !coarse part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),&
        Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
        Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
        Llr(ilr)%wfd%keyg(1,1),Llr(ilr)%wfd%keyv(1),Llr(ilr)%outofzone(:))

   !fine part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
        Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),Llr(ilr)%outofzone(:))

END SUBROUTINE determine_wfd_periodicity


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
     write(*,'(1x,a,2(i6))')&
          'ERROR: incorrect number of coarse points examined for reducing the localisation region',&
          nvctr_check,nvctr
     stop
  end if

END SUBROUTINE num_segkeys_periodic

subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(3), intent(in) :: outofzone
  integer, dimension(nseg_loc), intent(out) :: keyv_loc
  integer, dimension(2,nseg_loc), intent(out) :: keyg_loc
  !local variables
  logical :: go1,go2,go3,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l
  integer :: ngridp

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

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

          nvctr_check=nvctr_check+1
          if (.not. lseg) then
!             print *,'         check:',i,i2,i3,i1l,i2l,i3l,ngridp
             nsrt=nsrt+1
             keyg_loc(1,nsrt)=ngridp
             keyv_loc(nsrt)=nvctr_check
          end if
          lseg=.true.
        else 
           if (lseg) then
!              print *,'in        else:',i,i2,i3,i1l,i2l,i3l,ngridp
              nend=nend+1
              keyg_loc(2,nend)=ngridp
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
!        print *,'in second else:',i,i2,i3,i1l,i2l,i3l,ngridp
        nend=nend+1
        keyg_loc(2,nend)=ngridp
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

END SUBROUTINE segkeys_periodic

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
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keyg(1,1),Olr(index)%wfd%keyv(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keyg(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
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
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keyg(1,1),Olr(index)%wfd%keyv(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keyg(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
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
      write(*,*)&
          'ERROR: problem in get_overlap_region_periodic ',&
          'index:',index,'not equal to isovrlp:',isovrlp,&
          'The number of overlap descriptors constructed does not',&
          'correspond to the number of overlap regions.'
     stop
  end if

END SUBROUTINE get_overlap_region_periodic2


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
subroutine determine_locreg_parallel(iproc,nproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,orbs)!,outofzone)
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
        call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
             Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
             Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
     end if
  end do !on iilr

!  call make_LLr_MpiType(Llr,nlr,mpiLlr)

!  call MPI_ALLREDUCE(Llr(1),Llr(1),nlr,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
  !after all localisation regions are determined draw them
  !call draw_locregs(nlr,hx,hy,hz,Llr)

END SUBROUTINE determine_locreg_parallel

