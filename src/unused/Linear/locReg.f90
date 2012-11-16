!!$subroutine determine_locregSphere(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,calculateBounds)!,outofzone)
!!$  use module_base
!!$  use module_types
!!$  implicit none
!!$  integer, intent(in) :: iproc
!!$  integer, intent(in) :: nlr
!!$  real(gp), intent(in) :: hx,hy,hz
!!$  type(locreg_descriptors), intent(in) :: Glr
!!$  real(gp), dimension(nlr), intent(in) :: locrad
!!$  real(gp), dimension(3,nlr), intent(in) :: cxyz
!!$  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
!!$  logical,dimension(nlr),intent(in):: calculateBounds
!!$!  integer, dimension(3,nlr),intent(out) :: outofzone
!!$  !local variables
!!$  character(len=*), parameter :: subname='determine_locreg'
!!$  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz, calculate
!!$  logical :: warningx,warningy,warningz
!!$  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
!!$  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
!!$  integer :: ilr,isx,isy,isz,iex,iey,iez
!!$  integer :: ln1,ln2,ln3, iorb, jproc, jlr
!!$  integer :: ierr 
!!$  integer,dimension(3) :: outofzone
!!$  real(gp) :: rx,ry,rz,cutoff  
!!$
!!$
!!$  !initialize out of zone and logicals
!!$  outofzone (:) = 0     
!!$  warningx = .false.
!!$  warningy = .false.
!!$  warningz = .false.  
!!$
!!$  !determine the limits of the different localisation regions
!!$  do ilr=1,nlr
!!$
!!$     rx=cxyz(1,ilr)
!!$     ry=cxyz(2,ilr)
!!$     rz=cxyz(3,ilr)
!!$     llr(ilr)%locregCenter(1)=rx
!!$     llr(ilr)%locregCenter(2)=ry
!!$     llr(ilr)%locregCenter(3)=rz
!!$
!!$     cutoff=locrad(ilr)
!!$     llr(ilr)%locrad=cutoff
!!$
!!$     ! Determine the extrema of this localization regions (using only the coarse part, since this is always larger or equal than the fine part).
!!$     call determine_boxbounds_sphere(glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, hx, hy, hz, &
!!$          cutoff, llr(ilr)%locregCenter, &
!!$           glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyvloc, isx, isy, isz, iex, iey, iez)
!!$
!!$     ln1 = iex-isx
!!$     ln2 = iey-isy
!!$     ln3 = iez-isz
!!$
!!$
!!$     ! Localization regions should always have free boundary conditions
!!$     Llr(ilr)%geocode='F'
!!$
!!$     !assign the starting/ending points and outofzone for the different
!!$     ! geometries
!!$     select case(Glr%geocode)
!!$     case('F')
!!$        isx=max(isx,Glr%ns1)
!!$        isy=max(isy,Glr%ns2)
!!$        isz=max(isz,Glr%ns3)
!!$
!!$        iex=min(iex,Glr%ns1+Glr%d%n1)
!!$        iey=min(iey,Glr%ns2+Glr%d%n2)
!!$        iez=min(iez,Glr%ns3+Glr%d%n3)
!!$
!!$     case('S')
!!$        ! Get starting and ending for x direction     
!!$        if (iex - isx >= Glr%d%n1) then       
!!$           isx=Glr%ns1
!!$           iex=Glr%ns1 + Glr%d%n1
!!$        else
!!$           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
!!$           iex= ln1 + isx
!!$           if (iex > Glr%ns1+Glr%d%n1) then
!!$              outofzone(1)=modulo(iex,Glr%d%n1+1)
!!$           end if           
!!$        end if
!!$        
!!$        ! Get starting and ending for y direction (perpendicular to surface)
!!$        isy=max(isy,Glr%ns2)
!!$        iey=min(iey,Glr%ns2 + Glr%d%n2)
!!$        outofzone(2) = 0
!!$
!!$        !Get starting and ending for z direction
!!$        if (iez - isz >= Glr%d%n3) then
!!$           isz=Glr%ns3 
!!$           iez=Glr%ns3 + Glr%d%n3
!!$        else
!!$           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
!!$           iez= ln3 + isz
!!$           if (iez > Glr%ns3+Glr%d%n3) then
!!$              outofzone(3)=modulo(iez,Glr%d%n3+1)
!!$           end if 
!!$        end if
!!$
!!$     case('P')
!!$         ! Get starting and ending for x direction     
!!$        if (iex - isx >= Glr%d%n1) then       
!!$           isx=Glr%ns1
!!$           iex=Glr%ns1 + Glr%d%n1
!!$        else
!!$           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
!!$           iex= ln1 + isx
!!$           if (iex > Glr%ns1+Glr%d%n1) then
!!$              outofzone(1)=modulo(iex,Glr%d%n1+1)
!!$           end if           
!!$        end if
!!$        
!!$        ! Get starting and ending for y direction (perpendicular to surface)
!!$        if (iey - isy >= Glr%d%n2) then       
!!$           isy=Glr%ns2
!!$           iey=Glr%ns2 + Glr%d%n2
!!$         else
!!$           isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
!!$           iey= ln2 + isy
!!$           if (iey > Glr%ns2+Glr%d%n2) then
!!$              outofzone(2)=modulo(iey,Glr%d%n2+1)
!!$           end if           
!!$        end if
!!$
!!$        !Get starting and ending for z direction
!!$        if (iez - isz >= Glr%d%n3) then
!!$           isz=Glr%ns3 
!!$           iez=Glr%ns3 + Glr%d%n3
!!$        else
!!$           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
!!$           iez= ln3 + isz
!!$           if (iez > Glr%ns3+Glr%d%n3) then
!!$              outofzone(3)=modulo(iez,Glr%d%n3+1)
!!$           end if 
!!$        end if
!!$     end select
!!$
!!$     !values for the starting point of the cube for wavelet grid
!!$     Llr(ilr)%ns1=isx
!!$     Llr(ilr)%ns2=isy
!!$     Llr(ilr)%ns3=isz
!!$
!!$     !dimensions of the localisation region
!!$     Llr(ilr)%d%n1=iex-isx
!!$     Llr(ilr)%d%n2=iey-isy
!!$     Llr(ilr)%d%n3=iez-isz
!!$
!!$     !assign outofzone
!!$     Llr(ilr)%outofzone(:) = outofzone(:)
!!$
!!$     ! Set the conditions for ext_buffers (conditions for buffer size)
!!$     Gperx=(Glr%geocode /= 'F')
!!$     Gpery=(Glr%geocode == 'P')
!!$     Gperz=(Glr%geocode /= 'F')
!!$     Lperx=(Llr(ilr)%geocode /= 'F')
!!$     Lpery=(Llr(ilr)%geocode == 'P')
!!$     Lperz=(Llr(ilr)%geocode /= 'F')
!!$
!!$     !calculate the size of the buffers of interpolating function grid
!!$     call ext_buffers(Gperx,Gnbl1,Gnbr1)
!!$     call ext_buffers(Gpery,Gnbl2,Gnbr2)
!!$     call ext_buffers(Gperz,Gnbl3,Gnbr3)
!!$     call ext_buffers(Lperx,Lnbl1,Lnbr1)
!!$     call ext_buffers(Lpery,Lnbl2,Lnbr2)
!!$     call ext_buffers(Lperz,Lnbl3,Lnbr3)
!!$
!!$     !starting point of the region for interpolating functions grid
!!$     Llr(ilr)%nsi1= 2 * Llr(ilr)%ns1 - (Lnbl1 - Gnbl1)
!!$     Llr(ilr)%nsi2= 2 * Llr(ilr)%ns2 - (Lnbl2 - Gnbl2)
!!$     Llr(ilr)%nsi3= 2 * Llr(ilr)%ns3 - (Lnbl3 - Gnbl3)
!!$
!!$     !dimensions of the fine grid inside the localisation region
!!$     Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
!!$     Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
!!$     Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
!!$     
!!$     !NOTE: This will not work with symmetries (must change it)
!!$     Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
!!$     Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
!!$     Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz
!!$
!!$     !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
!!$     Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
!!$     Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
!!$     Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31
!!$
!!$!DEBUG
!!$!     if (iproc == 0) then
!!$!        write(*,*)'Description of zone:',ilr
!!$!        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
!!$!        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
!!$!        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!!$!        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
!!$!        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
!!$!        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
!!$!        write(*,*)'outofzone',ilr,':',outofzone(:)
!!$!     end if
!!$!DEBUG
!!$
!!$    ! construct the wavefunction descriptors (wfd)
!!$     call determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)
!!$
!!$     ! Sould check if nfu works properly... also relative to locreg!!
!!$     !if the localisation region is isolated build also the bounds
!!$     if (Llr(ilr)%geocode=='F') then
!!$        ! Check whether the bounds shall be calculated. Do this only if the currect process handles
!!$        ! orbitals in the current localization region.
!!$        if(calculateBounds(ilr)) then
!!$            call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
!!$                 Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
!!$                 Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
!!$        end if
!!$     end if
!!$  end do !on ilr
!!$
!!$END SUBROUTINE determine_locregSphere

!!!> This routine generates the keyglob, keygloc and keyv for the localization regions using periodic boundary conditions
!!!! The keys are continous is the localization region (keygloc), while they are inverted global region (keyglob). 
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

!> Determines the number of intersection regions between locregs, 
!! taking into account the periodicity of the system.
!!subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr)!,outofzone)
!!
!!  use module_base
!!  use module_types
!! 
!!  implicit none
!!
!!  ! Subroutine Scalar Arguments
!!  integer, intent(in) :: alr,blr              ! index of the two localization regions
!!  integer, intent(in) :: nlr                  ! number of localization regions
!!  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!  integer, intent(out) :: isovrlp             ! Integer giving the number of overlaps (max 8 with periodicity)
  
!!  !Subroutine Array Arguments
!!!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
!!  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 

!!  !local variables
!!  integer :: ii,azones,bzones,i_stat,i_all
!!  integer :: izones,jzones
!!  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
!!  character(len=*), parameter :: subname='get_number_of_overlap_region'
!!  logical :: go1,go2,go3
!!
!!  azones = 1
!!  bzones = 1
!!! Calculate the number of regions to cut alr and blr
!!  do ii=1,3
!!     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
!!     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
!!  end do
!!
!!!write(*,*)'azones,bzones',azones,bzones
!!!write(*,*)'outofzone',alr,':',outofzone(:,alr)
!!!write(*,*)'outofzone',blr,':',outofzone(:,blr)
!!
!!!allocate astart and aend
!!  allocate(astart(3,azones),stat=i_stat)
!!  call memocc(i_stat,astart,'astart',subname)
!!  allocate(aend(3,azones),stat=i_stat)
!!  call memocc(i_stat,aend,'aend',subname)
!!
!!!FRACTURE THE FIRST LOCALIZATION REGION
!!  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)
!!
!!!allocate bstart and bend
!!  allocate(bstart(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bstart,'bstart',subname)
!!  allocate(bend(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bend,'bend',subname)
!!
!!!FRACTURE SECOND LOCREG
!!  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)
!!
!!! Now check the number of overlapping zones
!! isovrlp = 0
!! do izones=1,azones
!!   do jzones=1,bzones
!!      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
!!      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
!!      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
!!      if(go1 .and. go2 .and. go3) then
!!         isovrlp = isovrlp + 1
!!      end if
!!   end do
!! end do
!!
!!! Deallocation block
!!  i_all = -product(shape(astart))*kind(astart)
!!  deallocate(astart,stat=i_stat)
!!  call memocc(i_stat,i_all,'astart',subname)
!!  i_all = -product(shape(aend))*kind(aend)
!!  deallocate(aend,stat=i_stat)
!!  call memocc(i_stat,i_all,'aend',subname)
!!  i_all = -product(shape(bstart))*kind(bstart)
!!  deallocate(bstart,stat=i_stat)
!!  call memocc(i_stat,i_all,'bstart',subname)
!!  i_all = -product(shape(bend))*kind(bend)
!!  deallocate(bend,stat=i_stat)
!!  call memocc(i_stat,i_all,'bend',subname)
!!  
!!END SUBROUTINE get_number_of_overlap_region

!> Given two localization regions, A and B, this routine returns a localization region 
!! corresponding to the intersection of A & B. 
!!subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr)
!!
!!  use module_base
!!  use module_types
!! 
!! implicit none
!!
!!  ! Subroutine Scalar Arguments
!!  integer, intent(in) :: alr,blr              ! index of the two localization regions
!!  integer, intent(in) :: nlr                  ! number of localization regions
!!  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!  integer, intent(in) :: isovrlp              ! Number of overlap regions
  
!!  !Subroutine Array Arguments
!!  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
!!  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  
!!  !local variables
!!  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
!!  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
!!  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
!!  character(len=*), parameter :: subname='get_overlap_region_periodic'
!!  !# NEW
!!  integer :: ii,azones,bzones,i_stat,i_all,index
!!  integer :: izones,jzones
!!  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
!!  logical :: go1,go2,go3
!!
!!  azones = 1
!!  bzones = 1
!!! Calculate the number of regions to cut alr and blr
!!  do ii=1,3
!!     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
!!     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
!!  end do
!!
!!!allocate astart and aend
!!  allocate(astart(3,azones),stat=i_stat)
!!  call memocc(i_stat,astart,'astart',subname)
!!  allocate(aend(3,azones),stat=i_stat)
!!  call memocc(i_stat,aend,'aend',subname)
!!
!!!FRACTURE THE FIRST LOCALIZATION REGION
!!  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)
!!
!!!allocate bstart and bend
!!  allocate(bstart(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bstart,'bstart',subname)
!!  allocate(bend(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bend,'bend',subname)
!!
!!!FRACTURE SECOND LOCREG
!!  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)
!!
!!! Now check the number of overlapping zones
!!  index = 0
!!  do izones=1,azones
!!    do jzones=1,bzones
!!      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
!!      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
!!      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
!!      if(go1 .and. go2 .and. go3) then
!!        index = index + 1
!!
!!! Now construct the Overlap localization region descriptor
!!! only if there is an overlap. The following only works
!!! when the previous test is successful. Note also that
!!! isx, isy and isz are necessarily in the Glr by construction
!!! of the Llrs, so don't need to test them.
!!         
!!        ! Determine the limits of the overlap region
!!        isx = max(astart(1,izones),bstart(1,jzones))
!!        isy = max(astart(2,izones),bstart(2,jzones))
!!        isz = max(astart(3,izones),bstart(3,jzones))
!!
!!        iex = min(aend(1,izones),bend(1,jzones))
!!        iey = min(aend(2,izones),bend(2,jzones))
!!        iez = min(aend(3,izones),bend(3,jzones))
!!
!!!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!!!       This could change the values of the bounds, so do it here
!!!       for now, in sandbox,put free boundary to all zones
!!        Olr(index)%geocode = 'F'  
!!
!!!       Values for the starting point of the cube
!!        Olr(index)%ns1 = isx
!!        Olr(index)%ns2 = isy
!!        Olr(index)%ns3 = isz
!!
!!!       Dimensions of the overlap region
!!        Olr(index)%d%n1 = iex - isx 
!!        Olr(index)%d%n2 = iey - isy 
!!        Olr(index)%d%n3 = iez - isz 
!!    
!!!       Dimensions of the fine grid inside the overlap region
!!        if (isx < iex) then
!!           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
!!           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
!!        else
!!           write(*,*)'1: Yet to be implemented (little effort?)'
!!           stop
!!        end if
!!
!!        if (isy < iey) then
!!           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
!!           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
!!        else
!!           write(*,*)'2: Yet to be implemented (little effort?)'
!!           stop
!!        end if
!!
!!        if (isz < iez) then
!!           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
!!           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
!!        else
!!           write(*,*)'3: Yet to be implemented (little effort?)'
!!           stop
!!        end if
!!
!!!       Dimensions of the interpolating scaling function grid 
!!!       (geocode already taken into acount because it is simple)
!!        select case(Olr(index)%geocode)
!!        case('F')
!!          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
!!          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
!!          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
!!        case('S')
!!          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
!!          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
!!          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
!!        case('P')
!!          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
!!          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
!!          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
!!        end select
!! 
!!!       Now define the wavefunction descriptors inside the overlap region
!!!       First calculate the number of points and segments for the region
!!!       Coarse part:
!!        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),&
!!         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!!!       Fine part:
!!        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
!!         Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!         Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)
!!
!!!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
!!        call allocate_wfd(Olr(index)%wfd,subname)
!!
!!!       At last, fill the wavefunction descriptors
!!!       Coarse part
!!        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!          Glr%wfd%nseg_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),&
!!          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
!!          Olr(index)%wfd%keygloc(1,1),Olr(index)%wfd%keyvloc(1))
!!!       Fine part
!!        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!          Glr%wfd%nseg_f,&
!!          Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
!!          Olr(index)%wfd%keygloc(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
!!          Olr(index)%wfd%keyvloc(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))
!!
!!!       If the localisation region is isolated build also the bounds
!!        if (Olr(index)%geocode=='F') then
!!           call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
!!             Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
!!             Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)
!!     
!!        end if
!!     end if ! go1 .and. go2 .and. go3
!!   end do !jzones
!! end do !izones
!!
!!! Deallocation block
!!  i_all = -product(shape(astart))*kind(astart)
!!  deallocate(astart,stat=i_stat)
!!  call memocc(i_stat,i_all,'astart',subname)
!!  i_all = -product(shape(aend))*kind(aend)
!!  deallocate(aend,stat=i_stat)
!!  call memocc(i_stat,i_all,'aend',subname)
!!  i_all = -product(shape(bstart))*kind(bstart)
!!  deallocate(bstart,stat=i_stat)
!!  call memocc(i_stat,i_all,'bstart',subname)
!!  i_all = -product(shape(bend))*kind(bend)
!!  deallocate(bend,stat=i_stat)
!!  call memocc(i_stat,i_all,'bend',subname)
!!
!!! Check on the number of zones
!!  if (index /= isovrlp) then
!!      write(*,*)&
!!          'ERROR: problem in get_overlap_region_periodic ',&
!!          'index:',index,'not equal to isovrlp:',isovrlp,&
!!          'The number of overlap descriptors constructed does not',&
!!          'correspond to the number of overlap regions.'
!!     stop
!!  end if
!!
!!END SUBROUTINE get_overlap_region_periodic

!> Given two localization regions, A and B, this routine returns a localization region 
!! corresponding to the intersection of A & B.
!! This is the same as get_overlap_region_periodic, but does not allocate the bound arrays to save memory.
subroutine get_overlap_region_periodic2(alr,blr,Glr,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types

 implicit none

  ! Subroutine Scalar Arguments
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(in) :: isovrlp              ! Number of overlap regions
  
  !Subroutine Array Arguments
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors
  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  
  !local variables
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
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keygloc(1,1),Olr(index)%wfd%keyvloc(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,&
          Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keygloc(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
          Olr(index)%wfd%keyvloc(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))

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


!!subroutine get_overlap_region_periodic2Sphere(alr,blr,Glr,hx,hy,hz,isovrlp,Llr,nlr,Olr)
!!
!!  use module_base
!!  use module_types
!!
!! implicit none
!!
!!  ! Subroutine Scalar Arguments
!!  integer, intent(in) :: alr,blr              ! index of the two localization regions
!!  integer, intent(in) :: nlr                  ! number of localization regions
!!  real(8),intent(in):: hx, hy, hz
!!  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!  integer, intent(in) :: isovrlp              ! Number of overlap regions
  
!!  !Subroutine Array Arguments
!!  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors
!!  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  
!!  !local variables
!!  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
!!  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
!!  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
!!  character(len=*), parameter :: subname='get_overlap_region_periodic2'
!!  !# NEW
!!  integer :: ii,azones,bzones,i_stat,i_all,index
!!  integer :: izones,jzones
!!  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
!!  logical :: go1,go2,go3
!!
!!  azones = 1
!!  bzones = 1
!!! Calculate the number of regions to cut alr and blr
!!  do ii=1,3
!!     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
!!     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
!!  end do
!!
!!
!!!allocate astart and aend
!!  allocate(astart(3,azones),stat=i_stat)
!!  call memocc(i_stat,astart,'astart',subname)
!!  allocate(aend(3,azones),stat=i_stat)
!!  call memocc(i_stat,aend,'aend',subname)
!!
!!!FRACTURE THE FIRST LOCALIZATION REGION
!!  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)
!!
!!!allocate bstart and bend
!!  allocate(bstart(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bstart,'bstart',subname)
!!  allocate(bend(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bend,'bend',subname)
!!
!!!FRACTURE SECOND LOCREG
!!  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)
!!
!!! Now check the number of overlapping zones
!!  index = 0
!!  do izones=1,azones
!!    do jzones=1,bzones
!!      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
!!      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
!!      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
!!      if(go1 .and. go2 .and. go3) then
!!        index = index + 1
!!
!!! Now construct the Overlap localization region descriptor
!!! only if there is an overlap. The following only works
!!! when the previous test is successful. Note also that
!!! isx, isy and isz are necessarily in the Glr by construction
!!! of the Llrs, so don't need to test them.
!!
!!        ! Determine the limits of the overlap region
!!        isx = max(astart(1,izones),bstart(1,jzones))
!!        isy = max(astart(2,izones),bstart(2,jzones))
!!        isz = max(astart(3,izones),bstart(3,jzones))
!!
!!        iex = min(aend(1,izones),bend(1,jzones))
!!        iey = min(aend(2,izones),bend(2,jzones))
!!        iez = min(aend(3,izones),bend(3,jzones))
!!
!!!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!!!       This could change the values of the bounds, so do it here
!!!       for now, in sandbox,put free boundary to all zones
!!        Olr(index)%geocode = 'F'
!!
!!!       Values for the starting point of the cube
!!        Olr(index)%ns1 = isx
!!        Olr(index)%ns2 = isy
!!        Olr(index)%ns3 = isz
!!
!!!       Dimensions of the overlap region
!!        Olr(index)%d%n1 = iex - isx
!!        Olr(index)%d%n2 = iey - isy
!!        Olr(index)%d%n3 = iez - isz
!!
!!!       Dimensions of the fine grid inside the overlap region
!!        if (isx < iex) then
!!           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
!!           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
!!        else
!!           write(*,*)'4: Yet to be implemented (little effort?), isx, iex', isx, iex
!!           stop
!!        end if
!!
!!        if (isy < iey) then
!!           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
!!           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
!!        else
!!           write(*,*)'5: Yet to be implemented (little effort?)'
!!           stop
!!        end if
!!
!!        if (isz < iez) then
!!           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
!!           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
!!        else
!!           write(*,*)'6: Yet to be implemented (little effort?)'
!!           stop
!!        end if
!!
!!!       Dimensions of the interpolating scaling function grid
!!!       (geocode already taken into acount because it is simple)
!!        select case(Olr(index)%geocode)
!!        case('F')
!!          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
!!          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
!!          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
!!        case('S')
!!          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
!!          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
!!          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
!!        case('P')
!!          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
!!          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
!!          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
!!        end select
!!
!!!       Now define the wavefunction descriptors inside the overlap region
!!!       First calculate the number of points and segments for the region
!!!       Coarse part:
!!        !call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!        ! Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
!!        ! Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!!        !!call num_segkeys_overlapSphere(olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
!!        !!     olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
!!        !!     olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
!!        !!     hx, hy, hz, llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
!!        !!     olr(index)%wfd%nseg_c, olr(index)%wfd%nvctr_c)
!!        call num_segkeys_sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
!!             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
!!             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
!!             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
!!             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
!!             hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
!!             llr(alr)%wfd%nseg_c, llr(alr)%wfd%keygloc(1,1), &
!!             llr(alr)%wfd%keyvloc(1), &
!!             olr(index)%wfd%nseg_c, olr(index)%wfd%nvctr_c)
!!        
!!!       Fine part:
!!        !call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!        ! Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
!!        ! Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!        ! Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!        ! Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)
!!        !!call num_segkeys_overlapSphere(olr(index)%d%nfl1, olr(index)%d%nfu1, &
!!        !!     olr(index)%d%nfl2, olr(index)%d%nfu2, &
!!        !!     olr(index)%d%nfl3, olr(index)%d%nfu3, &
!!        !!     hx, hy, hz, llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
!!        !!     olr(index)%wfd%nseg_f, olr(index)%wfd%nvctr_f)
!!        call num_segkeys_sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
!!             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
!!             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
!!             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
!!             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
!!             hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
!!             llr(alr)%wfd%nseg_f, llr(alr)%wfd%keygloc(1,llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
!!             llr(alr)%wfd%keyvloc(llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
!!             olr(index)%wfd%nseg_f, olr(index)%wfd%nvctr_f)
!!
!!!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
!!        call allocate_wfd(Olr(index)%wfd,subname)
!!
!!!       At last, fill the wavefunction descriptors
!!!       Coarse part
!!        !call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!        !  Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
!!        !  Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
!!        !  Olr(index)%wfd%keyg(1,1),Olr(index)%wfd%keyv(1))
!!        !!call segkeys_overlapSphere(olr(index)%d%n1, olr(index)%d%n2, olr(index)%d%n3, &
!!        !!     olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
!!        !!     olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
!!        !!     olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
!!        !!     olr(index)%wfd%nseg_c, hx, hy, hz, &
!!        !!     llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
!!        !!     olr(index)%wfd%keyg(1,1), olr(index)%wfd%keyv(1))
!!        call segkeys_Sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
!!             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
!!             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
!!             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
!!             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
!!             olr(index)%wfd%nseg_c, hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
!!             llr(alr)%wfd%nseg_c, llr(alr)%wfd%keygloc(1,1), &
!!             llr(alr)%wfd%keyvloc(1), &
!!             olr(index)%wfd%keygloc(1,1),olr(index)%wfd%keyglob(1,1), &
!!             olr(index)%wfd%keyvloc(1),olr(index)%wfd%keyvglob(1))
!!
!!!       Fine part
!!        !call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
!!        !  Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
!!        !  Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!        !  Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!        !  Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
!!        !  Olr(index)%wfd%keyg(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
!!        !  Olr(index)%wfd%keyv(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))
!!        !!call segkeys_overlapSphere(olr(index)%d%n1, olr(index)%d%n2, olr(index)%d%n3, &
!!        !!     olr(index)%d%nfl1, olr(index)%d%nfu1, &
!!        !!     olr(index)%d%nfl2, olr(index)%d%nfu2, &
!!        !!     olr(index)%d%nfl3, olr(index)%d%nfu3, &
!!        !!     olr(index)%wfd%nseg_f, hx, hy, hz, &
!!        !!     llr(alr)%locrad, llr(blr)%locrad, llr(alr)%locregCenter, llr(blr)%locregCenter, &
!!        !!     olr(index)%wfd%keyg(1,olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
!!        !!     olr(index)%wfd%keyv(olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)))
!!        call segkeys_Sphere(llr(alr)%d%n1, llr(alr)%d%n2, llr(alr)%d%n3, &
!!             llr(alr)%ns1, llr(alr)%ns2, llr(alr)%ns3, &
!!             olr(index)%ns1, olr(index)%ns1+olr(index)%d%n1, &
!!             olr(index)%ns2, olr(index)%ns2+olr(index)%d%n2, &
!!             olr(index)%ns3, olr(index)%ns3+olr(index)%d%n3, &
!!             olr(index)%wfd%nseg_f, hx, hy, hz, llr(blr)%locrad, llr(blr)%locregCenter, &
!!             llr(alr)%wfd%nseg_f, llr(alr)%wfd%keygloc(1,llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
!!             llr(alr)%wfd%keyvloc(llr(alr)%wfd%nseg_c+min(1,llr(alr)%wfd%nseg_f)), &
!!             olr(index)%wfd%keygloc(1,olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
!!             olr(index)%wfd%keyglob(1,olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
!!             olr(index)%wfd%keyvloc(olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)), &
!!             olr(index)%wfd%keyvglob(olr(index)%wfd%nseg_c+min(1,olr(index)%wfd%nseg_f)))
!!
!!!       If the localisation region is isolated build also the bounds
!!        !!!if (Olr(index)%geocode=='F') then
!!        !!!   call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
!!        !!!     Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
!!        !!!     Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)
!!
!!        !!!end if
!!     end if ! go1 .and. go2 .and. go3
!!   end do !jzones
!! end do !izones
!!
!!! Deallocation block
!!  i_all = -product(shape(astart))*kind(astart)
!!  deallocate(astart,stat=i_stat)
!!  call memocc(i_stat,i_all,'astart',subname)
!!  i_all = -product(shape(aend))*kind(aend)
!!  deallocate(aend,stat=i_stat)
!!  call memocc(i_stat,i_all,'aend',subname)
!!  i_all = -product(shape(bstart))*kind(bstart)
!!  deallocate(bstart,stat=i_stat)
!!  call memocc(i_stat,i_all,'bstart',subname)
!!  i_all = -product(shape(bend))*kind(bend)
!!  deallocate(bend,stat=i_stat)
!!  call memocc(i_stat,i_all,'bend',subname)
!!
!!! Check on the number of zones
!!  if (index /= isovrlp) then
!!      write(*,'(a,i0,a,i0,a)')&
!!          'ERROR: problem in get_overlap_region_periodic; index: ',index,' not equal to isovrlp: ',isovrlp,'&
!!          &. The number of overlap descriptors constructed does not correspond to the number of overlap regions.'
!!     stop
!!  end if
!!
!!END SUBROUTINE get_overlap_region_periodic2Sphere



!!$subroutine parallel_repartition_locreg(iproc,nproc,nlr,nlr_par,islr_par)
!!$  implicit none
!!$  integer, intent(in) :: iproc,nproc,nlr
!!$  integer, dimension(0:nproc-1), intent(out) :: nlr_par
!!$  integer, dimension(0:nproc-1), intent(out) :: islr_par
!!$  !local variables
!!$  integer :: jproc,numlr,difflr,ind
!!$
!!$  numlr = int(nlr/nproc)
!!$  difflr = nlr - numlr*nproc
!!$
!!$  nlr_par = numlr
!!$  ind = 0
!!$  do jproc=0,difflr-1
!!$     nlr_par(jproc) = nlr_par(jproc) + 1   
!!$     islr_par(jproc) = ind
!!$     ind = ind + nlr_par(jproc)
!!$  end do
!!$  do jproc=difflr,nproc
!!$     islr_par(jproc) = ind
!!$     ind = ind + nlr_par(jproc)
!!$  end do
!!$
!!$END SUBROUTINE parallel_repartition_locreg


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

!!!> Count for each orbital and each process the number of overlapping orbitals.
!!subroutine determine_overlap_from_descriptors2(iproc, nproc, orbs_a, orbs_b, lzd_a, lzd_b, op, comon)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => determine_overlap_from_descriptors2
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs_a, orbs_b
!!type(local_zone_descriptors),intent(in):: lzd_a, lzd_b
!!type(overlapParameters),intent(out):: op
!!type(p2pComms),intent(out):: comon
!!! Local variables
!!integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold
!!!integer :: is1, ie1, is2, ie2, is3, ie3
!!!integer :: js1, je1, js2, je2, js3, je3
!!integer :: iiorb, istat, iall, noverlaps, ierr
!!!logical:: ovrlpx, ovrlpy, ovrlpz
!!logical :: isoverlap
!!integer:: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp, i1, i2, jjorb, ii
!!integer :: onseg
!!logical,dimension(:,:,:),allocatable:: overlapMatrix
!!integer,dimension(:),allocatable:: noverlapsarr, displs, recvcnts, overlaps_comon
!!integer,dimension(:,:),allocatable:: overlaps_op
!!integer,dimension(:,:,:),allocatable :: overlaps_nseg
!!!integer,dimension(:,:,:),allocatable :: iseglist, jseglist
!!character(len=*),parameter:: subname='determine_overlap_from_descriptors2'
!!
!!allocate(overlapMatrix(orbs_b%norb,maxval(orbs_a%norb_par(:,0)),0:nproc-1), stat=istat)
!!call memocc(istat, overlapMatrix, 'overlapMatrix', subname)
!!allocate(noverlapsarr(orbs_a%norbp), stat=istat)
!!call memocc(istat, noverlapsarr, 'noverlapsarr', subname)
!!allocate(displs(0:nproc-1), stat=istat)
!!call memocc(istat, displs, 'displs', subname)
!!allocate(recvcnts(0:nproc-1), stat=istat)
!!call memocc(istat, recvcnts, 'recvcnts', subname)
!!allocate(overlaps_nseg(orbs_b%norb,orbs_a%norbp,2), stat=istat)
!!call memocc(istat, overlaps_nseg, 'overlaps_nseg', subname)
!!
!!
!!    overlapMatrix=.false.
!!    overlaps_nseg = 0
!!    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!    ilrold=-1
!!    do iorb=1,orbs_a%norbp
!!        ioverlaporb=0 ! counts the overlaps for the given orbital.
!!        iiorb=orbs_a%isorb+iorb
!!        ilr=orbs_a%inWhichLocreg(iiorb)
!!!        call get_start_and_end_indices(lzd_a%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!        do jorb=1,orbs_b%norb
!!            jlr=orbs_b%inWhichLocreg(jorb)
!!            call check_overlap_cubic_periodic(lzd_a%Glr,lzd_a%llr(ilr),lzd_b%llr(jlr),isoverlap)
!!!            call get_start_and_end_indices(lzd_a%llr(jlr), js1, je1, js2, je2, js3, je3)
!!!            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!!            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!!            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!!            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!             if(isoverlap) then
!!                ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
!!                ! Now explicitely check whether there is an overlap by using the descriptors.
!!                !!call overlapbox_from_descriptors(lzd_a%llr(ilr)%d%n1, lzd_a%llr(ilr)%d%n2, lzd_a%llr(ilr)%d%n3, &
!!                !!     lzd_a%llr(jlr)%d%n1, lzd_a%llr(jlr)%d%n2, lzd_a%llr(jlr)%d%n3, &
!!                !!     lzd_a%glr%d%n1, lzd_a%glr%d%n2, lzd_a%glr%d%n3, &
!!                !!     lzd_a%llr(ilr)%ns1, lzd_a%llr(ilr)%ns2, lzd_a%llr(ilr)%ns3, &
!!                !!     lzd_a%llr(jlr)%ns1, lzd_a%llr(jlr)%ns2, lzd_a%llr(jlr)%ns3, &
!!                !!     lzd_a%glr%ns1, lzd_a%glr%ns2, lzd_a%glr%ns3, &
!!                !!     lzd_a%llr(ilr)%wfd%nseg_c, lzd_a%llr(jlr)%wfd%nseg_c, &
!!                !!     lzd_a%llr(ilr)%wfd%keygloc, lzd_a%llr(ilr)%wfd%keyvloc, lzd_a%llr(jlr)%wfd%keygloc, lzd_a%llr(jlr)%wfd%keyvloc, &
!!                !!     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)
!!                call check_overlap_from_descriptors_periodic(lzd_a%llr(ilr)%wfd%nseg_c, lzd_b%llr(jlr)%wfd%nseg_c,&
!!                     lzd_a%llr(ilr)%wfd%keyglob, lzd_b%llr(jlr)%wfd%keyglob, &
!!                     isoverlap, onseg)
!!                if(isoverlap) then
!!                    ! There is really an overlap
!!                    overlapMatrix(jorb,iorb,iproc)=.true.
!!                    ioverlaporb=ioverlaporb+1
!!                    overlaps_nseg(ioverlaporb,iorb,1)=onseg
!!                    if(ilr/=ilrold) then
!!                        ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
!!                        ! would get the same orbitals again. Therefore the counter is not increased
!!                        ! in that case.
!!                        ioverlapMPI=ioverlapMPI+1
!!                    end if
!!                else
!!                    overlapMatrix(jorb,iorb,iproc)=.false.
!!                end if
!!             else
!!                overlapMatrix(jorb,iorb,iproc)=.false.
!!             end if
!!        end do
!!        noverlapsarr(iorb)=ioverlaporb
!!        ilrold=ilr
!!    end do
!!    !comon%noverlaps(jproc)=ioverlapMPI
!!    noverlaps=ioverlapMPI
!!
!!!call mpi_allreduce(overlapMatrix, orbs_a%norb*maxval(orbs_a%norb_par(:,0))*nproc, mpi_sum mpi_comm_world, ierr)
!!
!!! Communicate op%noverlaps and comon%noverlaps
!!    if (nproc > 1) then
!!       call mpi_allgatherv(noverlapsarr, orbs_a%norbp, mpi_integer, op%noverlaps, orbs_a%norb_par, &
!!            orbs_a%isorb_par, mpi_integer, mpi_comm_world, ierr)
!!    else
!!       call vcopy(orbs_a%norb,noverlapsarr(1),1,op%noverlaps(1),1)
!!    end if
!!    do jproc=0,nproc-1
!!       recvcnts(jproc)=1
!!       displs(jproc)=jproc
!!    end do
!!    if (nproc > 1) then
!!       call mpi_allgatherv(noverlaps, 1, mpi_integer, comon%noverlaps, recvcnts, &
!!            displs, mpi_integer, mpi_comm_world, ierr)
!!    else
!!       comon%noverlaps=noverlaps
!!    end if
!!
!!
!!allocate(op%overlaps(maxval(op%noverlaps),orbs_a%norb), stat=istat)
!!call memocc(istat, op%overlaps, 'op%overlaps', subname)
!!!!allocate(comon%overlaps(maxval(comon%noverlaps),0:nproc-1), stat=istat)
!!!!call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
!!
!!allocate(overlaps_op(maxval(op%noverlaps),orbs_a%norbp), stat=istat)
!!call memocc(istat, overlaps_op, 'overlaps_op', subname)
!!allocate(overlaps_comon(comon%noverlaps(iproc)), stat=istat)
!!call memocc(istat, overlaps_comon, 'overlaps_comon', subname)
!!
!!! Now we know how many overlaps have to be calculated, so determine which orbital overlaps
!!! with which one. This is essentially the same loop as above, but we use the array 'overlapMatrix'
!!! which indicates the overlaps.
!!
!!! Initialize to some value which will never be used.
!!op%overlaps=-1
!!!!comon%overlaps=-1
!!
!!iiorb=0
!!ioverlapMPI=0 ! counts the overlaps for the given MPI process.
!!ilrold=-1
!!do iorb=1,orbs_a%norbp
!!    ioverlaporb=0 ! counts the overlaps for the given orbital.
!!    iiorb=orbs_a%isorb+iorb
!!    ilr=orbs_a%inWhichLocreg(iiorb)
!!    do jorb=1,orbs_b%norb
!!        jlr=orbs_b%inWhichLocreg(jorb)
!!        if(overlapMatrix(jorb,iorb,iproc)) then
!!            ioverlaporb=ioverlaporb+1
!!            ! Determine the number of segments of the fine grid in the overlap
!!            call check_overlap_from_descriptors_periodic(lzd_a%llr(ilr)%wfd%nseg_c, lzd_b%llr(jlr)%wfd%nseg_f,&
!!                 lzd_a%llr(ilr)%wfd%keyglob(1,1), lzd_b%llr(jlr)%wfd%keyglob(1,1+lzd_b%llr(jlr)%wfd%nseg_c), &
!!                 isoverlap, overlaps_nseg(ioverlaporb,iorb,2))
!!            !op%overlaps(ioverlaporb,iiorb)=jorb
!!            overlaps_op(ioverlaporb,iorb)=jorb
!!            if(ilr/=ilrold) then
!!                ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
!!                ! would get the same orbitals again. Therefore the counter is not increased
!!                ! in that case.
!!                ioverlapMPI=ioverlapMPI+1
!!                !comon%overlaps(ioverlapMPI,iproc)=jorb
!!                overlaps_comon(ioverlapMPI)=jorb
!!            end if
!!        end if
!!    end do 
!!    ilrold=ilr
!!end do
!!
!!! Allocate the overlap wavefunctions_descriptors
!!! and copy the nseg_c
!!op%noverlapsmaxp=maxval(op%noverlaps(orbs_a%isorb+1:orbs_a%isorb+orbs_a%norbp))
!!allocate(op%wfd_overlap(op%noverlapsmaxp,orbs_a%norbp), stat=istat)
!!do i2=1,orbs_a%norbp
!!    do i1=1,op%noverlapsmaxp
!!        call nullify_wavefunctions_descriptors(op%wfd_overlap(i1,i2))
!!        op%wfd_overlap(i1,i2)%nseg_c = overlaps_nseg(i1,i2,1)
!!        op%wfd_overlap(i1,i2)%nseg_f = overlaps_nseg(i1,i2,2)
!!        call allocate_wfd(op%wfd_overlap(i1,i2),subname)
!!    end do
!!end do
!!
!!!Now redo the loop for the keygs
!!iiorb=0
!!do iorb=1,orbs_a%norbp
!!    ioverlaporb=0 ! counts the overlaps for the given orbital.
!!    iiorb=orbs_a%isorb+iorb
!!    ilr=orbs_a%inWhichLocreg(iiorb)
!!    do jorb=1,orbs_b%norb
!!        jlr=orbs_b%inWhichLocreg(jorb)
!!        if(overlapMatrix(jorb,iorb,iproc)) then
!!            ioverlaporb=ioverlaporb+1
!!           ! Determine the keyglob, keyvglob, nvctr of the coarse grid
!!           call get_overlap_from_descriptors_periodic(lzd_a%llr(ilr)%wfd%nseg_c, lzd_b%llr(jlr)%wfd%nseg_c, &
!!                lzd_a%llr(ilr)%wfd%keyglob(1,1), lzd_b%llr(jlr)%wfd%keyglob(1,1),  &
!!                .true.,op%wfd_overlap(ioverlaporb,iorb)%nseg_c, op%wfd_overlap(ioverlaporb,iorb)%nvctr_c,&
!!                op%wfd_overlap(ioverlaporb,iorb)%keyglob(1,1), op%wfd_overlap(ioverlaporb,iorb)%keyvglob(1))
!!           ! Determine the keyglob, keyvglob, nvctr of the fine grid
!!           if(op%wfd_overlap(ioverlaporb,iorb)%nseg_f > 0) then
!!              call get_overlap_from_descriptors_periodic(lzd_a%llr(ilr)%wfd%nseg_c, lzd_b%llr(jlr)%wfd%nseg_f, &
!!                   lzd_a%llr(ilr)%wfd%keyglob(1,1), lzd_b%llr(jlr)%wfd%keyglob(1,1+lzd_b%llr(jlr)%wfd%nseg_c),  &
!!                   .true.,op%wfd_overlap(ioverlaporb,iorb)%nseg_f, op%wfd_overlap(ioverlaporb,iorb)%nvctr_f,&
!!                   op%wfd_overlap(ioverlaporb,iorb)%keyglob(1,op%wfd_overlap(ioverlaporb,iorb)%nseg_c+1), &
!!                   op%wfd_overlap(ioverlaporb,iorb)%keyvglob(op%wfd_overlap(ioverlaporb,iorb)%nseg_c+1))
!!           else
!!              op%wfd_overlap(ioverlaporb,iorb)%nvctr_f = 0
!!           end if
!!        end if
!!    end do 
!!end do
!!
!!
!!displs(0)=0
!!recvcnts(0)=comon%noverlaps(0)
!!do jproc=1,nproc-1
!!    recvcnts(jproc)=comon%noverlaps(jproc)
!!    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
!!end do
!!!!if (nproc > 1) then
!!!!   call mpi_allgatherv(overlaps_comon, comon%noverlaps(iproc), mpi_integer, comon%overlaps, recvcnts, &
!!!!        displs, mpi_integer, mpi_comm_world, ierr)
!!!!else
!!!!   call vcopy(comon%noverlaps(iproc),overlaps_comon(1),1,comon%overlaps(1,0),1)
!!!!end if
!!ii=maxval(op%noverlaps)
!!displs(0)=0
!!recvcnts(0)=ii*orbs_a%norb_par(0,0)
!!do jproc=1,nproc-1
!!    recvcnts(jproc)=ii*orbs_a%norb_par(jproc,0)
!!    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
!!end do
!!if (nproc > 1) then
!!   call mpi_allgatherv(overlaps_op, ii*orbs_a%norbp, mpi_integer, op%overlaps, recvcnts, &
!!        displs, mpi_integer, mpi_comm_world, ierr)
!!else
!!   call vcopy(ii*orbs_a%norbp,overlaps_op(1,1),1,op%overlaps(1,1),1)
!!end if
!!
!!!!do iorb=1,orbs_a%norb
!!!!end do
!!
!!
!!iall=-product(shape(overlapMatrix))*kind(overlapMatrix)
!!deallocate(overlapMatrix, stat=istat)
!!call memocc(istat, iall, 'overlapMatrix', subname)
!!
!!iall=-product(shape(noverlapsarr))*kind(noverlapsarr)
!!deallocate(noverlapsarr, stat=istat)
!!call memocc(istat, iall, 'noverlapsarr', subname)
!!
!!iall=-product(shape(overlaps_op))*kind(overlaps_op)
!!deallocate(overlaps_op, stat=istat)
!!call memocc(istat, iall, 'overlaps_op', subname)
!!
!!iall=-product(shape(overlaps_comon))*kind(overlaps_comon)
!!deallocate(overlaps_comon, stat=istat)
!!call memocc(istat, iall, 'overlaps_comon', subname)
!!
!!iall=-product(shape(displs))*kind(displs)
!!deallocate(displs, stat=istat)
!!call memocc(istat, iall, 'displs', subname)
!!
!!iall=-product(shape(recvcnts))*kind(recvcnts)
!!deallocate(recvcnts, stat=istat)
!!call memocc(istat, iall, 'recvcnts', subname)
!!
!!iall=-product(shape(overlaps_nseg))*kind(overlaps_nseg)
!!deallocate(overlaps_nseg, stat=istat)
!!call memocc(istat, iall, 'overlaps_nseg', subname)
!!
!!end subroutine determine_overlap_from_descriptors2

!> Returns the starting and ending indices (on the coarse grid) of a given localization region.
!! Attention: There is a similar routine that uses "is1=lr%ns1+1" instead of "is1=lr%ns1".
!! Check why there is a difference of 1.
!!subroutine get_start_and_end_indices(lr, is1, ie1, is2, ie2, is3, ie3)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(locreg_descriptors),intent(in):: lr
!!integer,intent(out):: is1, ie1, is2, ie2, is3, ie3
!!
!!  is1=lr%ns1
!!  ie1=lr%ns1+lr%d%n1
!!  is2=lr%ns2
!!  ie2=lr%ns2+lr%d%n2
!!  is3=lr%ns3
!!  ie3=lr%ns3+lr%d%n3
!!
!!end subroutine get_start_and_end_indices


!!$function check_whether_locregs_overlap(llr_i, llr_j, glr)
!!$use module_base
!!$use module_types
!!$use module_interfaces
!!$implicit none
!!$
!!$! Calling arguments
!!$type(locreg_descriptors),intent(in):: llr_i, llr_j, glr
!!$logical:: check_whether_locregs_overlap
!!$
!!$! Local variables
!!$integer:: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp
!!$
!!$  ! Check whether there is an overlap by comparing the descriptors.
!!$  call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
!!$       llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
!!$       glr%d%n1, glr%d%n2, glr%d%n3, &
!!$       llr_i%ns1, llr_i%ns2, llr_i%ns3, &
!!$       llr_j%ns1, llr_j%ns2, llr_j%ns3, &
!!$       glr%ns1, glr%ns2, glr%ns3, &
!!$       llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, &
!!$       llr_i%wfd%keygloc, llr_i%wfd%keyvloc, llr_j%wfd%keygloc, llr_j%wfd%keyvloc, &
!!$       n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)
!!$
!!$  ! n1_ovrlp, n2_ovrlp, n3_ovrlp are the dimensions of the overlap localization regions.
!!$  if(n1_ovrlp>0 .and. n2_ovrlp>0 .and. n3_ovrlp>0) then
!!$      ! There is an overlap
!!$      check_whether_locregs_overlap=.true.
!!$  else
!!$      ! There is no overlap
!!$      check_whether_locregs_overlap=.false.
!!$  end if
!!$end function check_whether_locregs_overlap


!> @file
!! Locallisation regions
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Determine a set of localisation regions from the centers and the radii.
!! cut in cubes the global reference system
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



  !determine the limits of the different localisation regions
  do ilr=1,nlr

     !initialize out of zone and logicals
     outofzone (:) = 0     
     warningx = .false.
     warningy = .false.
     warningz = .false.  

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
     if(Llr(ilr)%geocode == 'F') then
        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31
     else if(Llr(ilr)%geocode == 'S') then
        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+2
        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+2
     else
        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+2
        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+2
        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+2
     end if

!DEBUG
!!     if (iproc == 0) then
!!        write(*,*)'Description of zone:',ilr
!!        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
!!        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
!!        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!!        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
!!        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
!!        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
!!        write(*,*)'outofzone',ilr,':',outofzone(:)
!!     end if
!DEBUG

    ! construct the wavefunction descriptors (wfd)
     call determine_wfd_periodicity(ilr,nlr,Glr,Llr)

     ! Sould check if nfu works properly... also relative to locreg!!
     !if the localisation region is isolated build also the bounds
     if (Llr(ilr)%geocode=='F') then
        ! Check whether the bounds shall be calculated. Do this only if the currect process handles
        ! orbitals in the current localization region.
        if(calculateBounds(ilr)) then
!           print *,'===>ilr',ilr,Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,Llr(ilr)%outofzone
            call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                 Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                 Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
        end if
     end if
  end do !on ilr

END SUBROUTINE determine_locreg_periodic


