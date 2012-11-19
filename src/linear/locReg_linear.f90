!> Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg 
!! taking into account the pediodicity
!!          
!! @warning
!!    We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 
  
  !Subroutine Array Arguments
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr

  !local variables
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
          Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),nseg_c,nvctr_c,Llr(ilr)%outofzone(:))
   !fine part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),nseg_f,nvctr_f,Llr(ilr)%outofzone(:))

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
        Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),&
        Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
        Llr(ilr)%wfd%keygloc(1,1),Llr(ilr)%wfd%keyglob(1,1),Llr(ilr)%wfd%keyvloc(1),&
        Llr(ilr)%wfd%keyvglob(1),&
        Llr(ilr)%outofzone(:))

   !fine part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
        Llr(ilr)%wfd%keygloc(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyglob(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyvloc(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyvglob(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%outofzone(:))


END SUBROUTINE determine_wfd_periodicity


subroutine determine_locregSphere_parallel(iproc,nproc,nlr,cxyz,locrad,hx,hy,hz,at,orbs,Glr,Llr,calculateBounds)!,outofzone)
  use module_base
  use module_types
  use module_interfaces, except_this_one => determine_locregSphere_parallel

  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data),intent(in) :: at
  type(orbitals_data),intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
  logical,dimension(nlr),intent(in) :: calculateBounds
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg'
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz
  logical :: warningx,warningy,warningz,xperiodic,yperiodic,zperiodic
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3
  integer :: ii, root, ierr, iall, istat, iorb, iat, norb, norbu, norbd, nspin, iilr
  integer,dimension(3) :: outofzone
  integer,dimension(:),allocatable :: rootarr, norbsperatom, norbsperlocreg
  real(8),dimension(:,:),allocatable :: locregCenter
  type(orbitals_data) :: orbsder

  allocate(rootarr(nlr), stat=istat)
  call memocc(istat, rootarr, 'rootarr', subname)

  ! Determine how many locregs one process handles at most
  ii=ceiling(dble(nlr)/dble(nproc))
  !determine the limits of the different localisation regions
  rootarr=1000000000
  iilr=0

  call timing(iproc,'wfd_creation  ','ON')  
  do ilr=1,nlr
     !initialize out of zone and logicals
     outofzone (:) = 0     
     warningx = .false.
     warningy = .false.
     warningz = .false. 
     xperiodic = .false.
     yperiodic = .false.
     zperiodic = .false. 

     !if(mod(ilr-1,nproc)==iproc) then
     if(orbs%onwhichmpi(ilr)==iproc) then
     !if (ilr>orbs%isorb .and. iilr<orbs%norbp) then
     !    iilr=iilr+1
     !if(calculateBounds(ilr) .or. (mod(ilr-1,nproc)==iproc)) then 
         ! This makes sure that each locreg is only handled once by one specific processor.
    
         llr(ilr)%locregCenter(1)=cxyz(1,ilr)
         llr(ilr)%locregCenter(2)=cxyz(2,ilr)
         llr(ilr)%locregCenter(3)=cxyz(3,ilr)
    
         llr(ilr)%locrad=locrad(ilr)
  
         ! Determine the extrema of this localization regions (using only the coarse part, since this is always larger or equal than the fine part).
         call determine_boxbounds_sphere(glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, hx, hy, hz, &
              llr(ilr)%locrad, llr(ilr)%locregCenter, &
               glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyvloc, isx, isy, isz, iex, iey, iez)
    
         ln1 = iex-isx
         ln2 = iey-isy
         ln3 = iez-isz
    
    
         ! Localization regions should have free boundary conditions by default
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
               xperiodic = .true.
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
               zperiodic = .true.
            else
               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               iez= ln3 + isz
               if (iez > Glr%ns3+Glr%d%n3) then
                  outofzone(3)=modulo(iez,Glr%d%n3+1)
               end if 
            end if
            if(xperiodic .and. zperiodic) then
              Llr(ilr)%geocode = 'S'
            end if    

         case('P')
             ! Get starting and ending for x direction     
            if (iex - isx >= Glr%d%n1) then       
               isx=Glr%ns1
               iex=Glr%ns1 + Glr%d%n1
               xperiodic = .true.
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
               yperiodic = .true.
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
               zperiodic = .true.
            else
               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               iez= ln3 + isz
               if (iez > Glr%ns3+Glr%d%n3) then
                  outofzone(3)=modulo(iez,Glr%d%n3+1)
               end if 
            end if
            if(xperiodic .and. yperiodic .and. zperiodic ) then
              Llr(ilr)%geocode = 'P'
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
        rootarr(ilr)=iproc
        call determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)
     end if
  end do !on ilr
  call timing(iproc,'wfd_creation  ','OF') 

  ! Communicate the locregs
  ! This communication is uneffective. Instead of using bcast we should be using mpialltoallv.
  call timing(iproc,'comm_llr      ','ON')
  if (nproc > 1) then
     call mpiallred(rootarr(1), nlr, mpi_min, bigdft_mpi%mpi_comm, ierr)
     
     ! Communicate those parts of the locregs that all processes need.
     call communicate_locreg_descriptors_basics(iproc, nlr, rootarr, orbs, llr)


     ! Now communicate those parts of the locreg that only some processes need (the keys).
     ! For this we first need to create orbsder that describes the derivatives.
     call create_orbsder()

     ! Now communicate the keys
     call communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, orbsder, rootarr)

     call deallocate_orbitals_data(orbsder, subname)
  end if
  call timing(iproc,'comm_llr      ','OF')


!create the bound arrays for the locregs we need on the MPI tasks
  call timing(iproc,'calc_bounds   ','ON') 
  do ilr=1,nlr
         if (Llr(ilr)%geocode=='F' .and. calculateBounds(ilr) ) then
            call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                 Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                 Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
         end if
  end do
  call timing(iproc,'calc_bounds   ','OF') 



iall = -product(shape(rootarr))*kind(rootarr)
deallocate(rootarr,stat=istat)
call memocc(istat,iall,'rootarr',subname)



contains 
  subroutine create_orbsder()
    call nullify_orbitals_data(orbsder)
    allocate(norbsperatom(at%nat), stat=istat)
    call memocc(istat, norbsperatom, 'norbsperatom', subname)
    allocate(locregCenter(3,nlr), stat=istat)
    call memocc(istat, locregCenter, 'locregCenter', subname)
    allocate(norbsPerLocreg(nlr), stat=istat) 
    call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)
    norbsperatom=0
    do iorb=1,orbs%norb
        iat=orbs%onwhichatom(iorb)
        norbsperatom(iat)=norbsperatom(iat)+3
    end do
    norb=3*orbs%norb
    norbu=norb
    norbd=0
    nspin=1
    call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, nspin, orbs%nspinor,&
         orbs%nkpts, orbs%kpts, orbs%kwgts, orbsder,.true.) !simple repartition
    iall=-product(shape(orbsder%onwhichatom))*kind(orbsder%inWhichLocreg)
    deallocate(orbsder%onwhichatom, stat=istat)
    call memocc(istat, iall, 'orbsder%onwhichatom', subname)
                 
    call assignToLocreg2(iproc, nproc, orbsder%norb, orbsder%norb_par, at%nat, at%nat, &
         nspin, norbsPerAtom, cxyz, orbsder%onwhichatom)


    do ilr=1,nlr
        locregCenter(:,ilr)=llr(ilr)%locregCenter
    end do

    iall=-product(shape(orbsder%inWhichLocreg))*kind(orbsder%inWhichLocreg)
    deallocate(orbsder%inWhichLocreg, stat=istat)
    norbsPerLocreg=3

    call memocc(istat, iall, 'orbsder%inWhichLocreg', subname)
    call assignToLocreg2(iproc, nproc, orbsder%norb, orbsder%norb_par, at%nat, nlr, &
         nspin, norbsPerLocreg, locregCenter, orbsder%inwhichlocreg)

    iall=-product(shape(locregCenter))*kind(locregCenter)
    deallocate(locregCenter, stat=istat)
    call memocc(istat, iall, 'locregCenter', subname)

    iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
    deallocate(norbsPerLocreg, stat=istat)
    call memocc(istat, iall, 'norbsPerLocreg', subname)

    iall=-product(shape(norbsperatom))*kind(norbsperatom)
    deallocate(norbsperatom, stat=istat)
    call memocc(istat, iall, 'norbsperatom', subname)
  end subroutine create_orbsder

END SUBROUTINE determine_locregSphere_parallel


!> Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg 
!! taking into account the pediodicity
!!          
!! @warning
!!    We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
subroutine determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  real(kind=8),intent(in) :: hx, hy, hz
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 
  
  !Subroutine Array Arguments
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr

  !local variables
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period
  character(len=*), parameter :: subname='determine_wfdSphere'
!!  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements

   !starting point of locreg (always inside global locreg)
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
        Glr%wfd%keyvloc(1), &
        llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nvctr_c)

   !fine part
   call num_segkeys_sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        glr%wfd%nseg_f, Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
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
        Glr%wfd%keyvloc(1), &
        llr(ilr)%wfd%keygloc(1,1),llr(ilr)%wfd%keyglob(1,1), &
        llr(ilr)%wfd%keyvloc(1), llr(ilr)%wfd%keyvglob(1))

   !fine part
   call segkeys_Sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        llr(ilr)%wfd%nseg_f, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_f, Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        llr(ilr)%wfd%keygloc(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyglob(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyvloc(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyvglob(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)))


END SUBROUTINE determine_wfdSphere


!> Calculates the number of segments and elements in localisation region
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
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,intent(out) :: nseg, nvctr
  !local variables
  logical :: segment
  integer :: i, i1, i2, i3, nstart, nend, i2old, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3
  real(kind=8) :: cut, dx,dy, dz


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
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,intent(out) :: ixmin, iymin, izmin, ixmax, iymax, izmax
  !local variables
  integer :: i, i1, i2, i3, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3
  real(kind=8) :: cut, dx,dy, dz
  !debug
  integer :: iiimin, isegmin
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
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l
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
!print *,'iseg,keygloc,keyg_loc',iseg,keygloc(1,loc),keygloc(2,loc),keyg_loc(1,iseg),keyg_loc(2,iseg)
    keyg_loc(1,loc) = maxval(keyg_loc) + 1
    keyvloc(iseg) = keyvglob(loc)
!    print *,'iseg,keyglob,keyvglob,keygloc,keyvloc',iseg,keyglob(1,iseg),keyvglob(iseg),keygloc(1,iseg),keyvloc(iseg)
 end do

 i_all = -product(shape(keyg_loc))*kind(keyg_loc)
 deallocate(keyg_loc, stat = i_stat)
 call memocc(i_stat,i_all,'keyg_loc',subname)

END SUBROUTINE segkeys_periodic


subroutine segkeys_Sphere(n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, keyg_loc, keyg_glob, keyv_loc, keyv_glob)
  use module_base
  implicit none
  integer,intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, nsegglob
  real(kind=8) :: hx, hy, hz, locrad
  real(kind=8),dimension(3) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,dimension(2,nseg),intent(out) :: keyg_loc, keyg_glob
  integer,dimension(nseg),intent(out) :: keyv_loc, keyv_glob
  !local variables
  character(len=*),parameter :: subname = 'segkeys_Sphere'
  integer :: i, i1, i2, i3, nstart, nend, nvctr, igridpoint, igridglob, i2old, iseg, jj, j0, j1, ii, i0, n1l, n2l, n3l
  integer :: i1l, i2l, i3l, ii1, ii2, ii3, istat, iall, loc
  real(kind=8) :: cut, dx, dy, dz
  logical :: segment
  integer, allocatable :: keygloc(:,:)

  allocate(keygloc(2,nseg),stat=istat)
  call memocc(istat,keygloc,'keygloc',subname)

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
                  keygloc(1,nstart)=igridpoint
                  keyg_glob(1,nstart)=igridglob
                  keyv_glob(nstart)=nvctr
                  segment=.true.
              end if
          else
              if(segment) then
                  nend=nend+1
                  keygloc(2,nend)=igridpoint-1
                  keyg_glob(2,nend)=igridglob-1
                  segment=.false.
              end if
          end if
      end do
      if(segment) then
          ! Close the segment
          nend=nend+1
          keygloc(2,nend)=igridpoint
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

  ! Now build the keyvloc where we replace the segments in order for the loc
  do iseg = 1, nseg
     !sorting the keyg_loc
     loc = minloc(keygloc(1,:),1)
     keyg_loc(1,iseg) = keygloc(1,loc)
     keyg_loc(2,iseg) = keygloc(2,loc)
!    print *,'iseg,keygloc,keyg_loc',iseg,keygloc(1,loc),keygloc(2,loc),keyg_loc(1,iseg),keyg_loc(2,iseg)
     keygloc(1,loc) = maxval(keygloc) + 1
     keyv_loc(iseg) = keyv_glob(loc)
!    print *,'iseg,keyglob,keyvglob,keygloc,keyvloc',iseg,keyglob(1,iseg),keyvglob(iseg),keygloc(1,iseg),keyvloc(iseg)
  end do
  iall = -product(shape(keygloc))*kind(keygloc)
  deallocate(keygloc,stat=istat)
  call memocc(istat,iall,'keygloc',subname)

END SUBROUTINE segkeys_Sphere

!> Divides the locreg into zones contained inside the simulation box, by applying the primitive vectors
!! It returns: astart(3,nzones) which is the starting points of the different zones (max. 8)
!!             aend(3,nzones) which is the ending points of the different zones (max. 8)
subroutine fracture_periodic_zone(nzones,Glr,Llr,outofzone,astart,aend)

  use module_base
  use module_types
 
  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: nzones
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Localization grid descriptors 
  
  !Subroutine Array Arguments
  integer,dimension(3),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  integer,dimension(3,nzones),intent(out) :: astart !
  integer,dimension(3,nzones),intent(out) :: aend !
  
  !local variables
  integer :: ii,index,jj
  integer,dimension(3) :: alrs,alre,Gend,Gstart,period
  
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

!> Determine a set of localisation regions from the centers and the radii.
!! cut in cubes the global reference system
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
  logical,dimension(nlr),intent(in) :: calculateBounds
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz
  logical :: warningx,warningy,warningz,calc
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez,iorb
  integer :: ln1,ln2,ln3
  integer,dimension(3) :: outofzone
  real(gp) :: rx,ry,rz,cutoff
!!  integer :: iilr,ierr
!!  integer,dimension(0:nproc-1) :: nlr_par,islr_par

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
!     nullify(Llr(ilr)%projflg)
     nullify(Llr(ilr)%wfd%keygloc)
     nullify(Llr(ilr)%wfd%keyglob)
     nullify(Llr(ilr)%wfd%keyvloc)
     nullify(Llr(ilr)%wfd%keyvglob)
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

!  call MPI_ALLREDUCE(Llr(1),Llr(1),nlr,mpidtypg,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  !after all localisation regions are determined draw them
  !call draw_locregs(nlr,hx,hy,hz,Llr)

END SUBROUTINE determine_locreg_parallel



!> Count for each orbital and each process the number of overlapping orbitals.
subroutine determine_overlap_from_descriptors(iproc, nproc, orbs, orbsig, lzd, lzdig, op, comon)
use module_base
use module_types
use module_interfaces, except_this_one => determine_overlap_from_descriptors
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
type(orbitals_data),intent(in) :: orbs, orbsig
type(local_zone_descriptors),intent(in) :: lzd, lzdig
type(overlapParameters),intent(inout) :: op
type(p2pComms),intent(inout) :: comon
! Local variables
integer :: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold
!integer :: is1, ie1, is2, ie2, is3, ie3
!integer :: js1, je1, js2, je2, js3, je3
integer :: iiorb, istat, iall, noverlaps, ierr
!logical :: ovrlpx, ovrlpy, ovrlpz
logical :: isoverlap
integer :: i1, i2, ii
integer :: onseg
logical,dimension(:,:,:),allocatable :: overlapMatrix
integer,dimension(:),allocatable :: noverlapsarr, displs, recvcnts, overlaps_comon
integer,dimension(:,:),allocatable :: overlaps_op
integer,dimension(:,:,:),allocatable :: overlaps_nseg
!integer,dimension(:,:,:),allocatable :: iseglist, jseglist
character(len=*),parameter :: subname='determine_overlap_from_descriptors'

allocate(overlapMatrix(orbsig%norb,maxval(orbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, overlapMatrix, 'overlapMatrix', subname)
allocate(noverlapsarr(orbs%norbp), stat=istat)
call memocc(istat, noverlapsarr, 'noverlapsarr', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)
allocate(recvcnts(0:nproc-1), stat=istat)
call memocc(istat, recvcnts, 'recvcnts', subname)
allocate(overlaps_nseg(orbsig%norb,orbs%norbp,2), stat=istat)
call memocc(istat, overlaps_nseg, 'overlaps_nseg', subname)

    overlapMatrix=.false.
    overlaps_nseg = 0
    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
    ilrold=-1
    do iorb=1,orbs%norbp
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=orbs%isorb+iorb
        ilr=orbs%inWhichLocreg(iiorb)
!        call get_start_and_end_indices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbsig%norb
            jlr=orbsig%inWhichLocreg(jorb)
            call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzdig%llr(jlr),isoverlap)
!            call get_start_and_end_indices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
             !write(*,*) 'second: iproc, ilr, jlr, isoverlap', iproc, ilr, jlr, isoverlap
             if(isoverlap) then
                ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
                ! Now explicitely check whether there is an overlap by using the descriptors.
                !!call overlapbox_from_descriptors(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
                !!     lzd%llr(jlr)%d%n1, lzd%llr(jlr)%d%n2, lzd%llr(jlr)%d%n3, &
                !!     lzd%glr%d%n1, lzd%glr%d%n2, lzd%glr%d%n3, &
                !!     lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
                !!     lzd%llr(jlr)%ns1, lzd%llr(jlr)%ns2, lzd%llr(jlr)%ns3, &
                !!     lzd%glr%ns1, lzd%glr%ns2, lzd%glr%ns3, &
                !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c, &
                !!     lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, lzd%llr(jlr)%wfd%keygloc, lzd%llr(jlr)%wfd%keyvloc, &
                !!     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)
                call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_c,&
                     lzd%llr(ilr)%wfd%keyglob, lzdig%llr(jlr)%wfd%keyglob, &
                     isoverlap, onseg)
                if(isoverlap) then
                    ! There is really an overlap
                    overlapMatrix(jorb,iorb,iproc)=.true.
                    ioverlaporb=ioverlaporb+1
                    overlaps_nseg(ioverlaporb,iorb,1)=onseg
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

!call mpi_allreduce(overlapMatrix, orbs%norb*maxval(orbs%norb_par(:,0))*nproc, mpi_sum bigdft_mpi%mpi_comm, ierr)

! Communicate op%noverlaps and comon%noverlaps
    if (nproc > 1) then
       call mpi_allgatherv(noverlapsarr, orbs%norbp, mpi_integer, op%noverlaps, orbs%norb_par, &
            orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    else
       call vcopy(orbs%norb,noverlapsarr(1),1,op%noverlaps(1),1)
    end if
    do jproc=0,nproc-1
       recvcnts(jproc)=1
       displs(jproc)=jproc
    end do
    if (nproc > 1) then
       call mpi_allgatherv(noverlaps, 1, mpi_integer, comon%noverlaps, recvcnts, &
            displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    else
       comon%noverlaps=noverlaps
    end if


allocate(op%overlaps(maxval(op%noverlaps),orbs%norb), stat=istat)
call memocc(istat, op%overlaps, 'op%overlaps', subname)
!!allocate(comon%overlaps(maxval(comon%noverlaps),0:nproc-1), stat=istat)
!!call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
allocate(overlaps_op(maxval(op%noverlaps),orbs%norbp), stat=istat)
call memocc(istat, overlaps_op, 'overlaps_op', subname)
if(orbs%norbp>0) call to_zero(maxval(op%noverlaps)*orbs%norbp,overlaps_op(1,1))
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
    do jorb=1,orbsig%norb
        jlr=orbsig%inWhichLocreg(jorb)
        if(overlapMatrix(jorb,iorb,iproc)) then
            ioverlaporb=ioverlaporb+1
            ! Determine the number of segments of the fine grid in the overlap
            call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_f,&
                 lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1+lzdig%llr(jlr)%wfd%nseg_c), &
                 isoverlap, overlaps_nseg(ioverlaporb,iorb,2))
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

! Allocate the overlap wavefunctions_descriptors
! and copy the nseg_c
op%noverlapsmaxp=maxval(op%noverlaps(orbs%isorb+1:orbs%isorb+orbs%norbp))
allocate(op%wfd_overlap(op%noverlapsmaxp,orbs%norbp), stat=istat)
do i2=1,orbs%norbp
    do i1=1,op%noverlapsmaxp
        call nullify_wavefunctions_descriptors(op%wfd_overlap(i1,i2))
        op%wfd_overlap(i1,i2)%nseg_c = overlaps_nseg(i1,i2,1)
        op%wfd_overlap(i1,i2)%nseg_f = overlaps_nseg(i1,i2,2)
        call allocate_wfd(op%wfd_overlap(i1,i2),subname)
    end do
end do

!Now redo the loop for the keygs
iiorb=0
do iorb=1,orbs%norbp
    ioverlaporb=0 ! counts the overlaps for the given orbital.
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    do jorb=1,orbsig%norb
        jlr=orbsig%inWhichLocreg(jorb)
        if(overlapMatrix(jorb,iorb,iproc)) then
            ioverlaporb=ioverlaporb+1
           ! Determine the keyglob, keyvglob, nvctr of the coarse grid
           call get_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_c, &
                lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1),  &
                .true.,op%wfd_overlap(ioverlaporb,iorb)%nseg_c, op%wfd_overlap(ioverlaporb,iorb)%nvctr_c,&
                op%wfd_overlap(ioverlaporb,iorb)%keyglob(1,1), op%wfd_overlap(ioverlaporb,iorb)%keyvglob(1))
           ! Determine the keyglob, keyvglob, nvctr of the fine grid
           if(op%wfd_overlap(ioverlaporb,iorb)%nseg_f > 0) then
              call get_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_f, &
                   lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1+lzdig%llr(jlr)%wfd%nseg_c),  &
                   .true.,op%wfd_overlap(ioverlaporb,iorb)%nseg_f, op%wfd_overlap(ioverlaporb,iorb)%nvctr_f,&
                   op%wfd_overlap(ioverlaporb,iorb)%keyglob(1,op%wfd_overlap(ioverlaporb,iorb)%nseg_c+1), &
                   op%wfd_overlap(ioverlaporb,iorb)%keyvglob(op%wfd_overlap(ioverlaporb,iorb)%nseg_c+1))
           else
              op%wfd_overlap(ioverlaporb,iorb)%nvctr_f = 0
           end if
        end if
    end do 
end do


displs(0)=0
recvcnts(0)=comon%noverlaps(0)
do jproc=1,nproc-1
    recvcnts(jproc)=comon%noverlaps(jproc)
    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
end do
!!if (nproc > 1) then
!!   call mpi_allgatherv(overlaps_comon, comon%noverlaps(iproc), mpi_integer, comon%overlaps, recvcnts, &
!!        displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
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
        displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
else
   call vcopy(ii*orbs%norbp,overlaps_op(1,1),1,op%overlaps(1,1),1)
end if


iall=-product(shape(overlapMatrix))*kind(overlapMatrix)
deallocate(overlapMatrix, stat=istat)
call memocc(istat, iall, 'overlapMatrix', subname)

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

iall=-product(shape(overlaps_nseg))*kind(overlaps_nseg)
deallocate(overlaps_nseg, stat=istat)
call memocc(istat, iall, 'overlaps_nseg', subname)

end subroutine determine_overlap_from_descriptors

subroutine determine_overlapdescriptors_from_descriptors(llr_i, llr_j, glr, olr)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in) :: llr_i, llr_j, glr
type(locreg_descriptors),intent(out) :: olr

! Local variables
integer :: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp
character(len=*),parameter :: subname='determine_overlapdescriptors_from_descriptors'



! Determine the values describing the localization region.
! First coarse region.
call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, &
     llr_i%wfd%keygloc, llr_i%wfd%keyvloc, llr_j%wfd%keygloc, llr_j%wfd%keyvloc, &
     olr%d%n1, olr%d%n2, olr%d%n3, olr%ns1, olr%ns2, olr%ns3, olr%wfd%nseg_c)

! Now the fine region.
call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     llr_i%wfd%nseg_f, llr_j%wfd%nseg_f, &
     llr_i%wfd%keygloc(1,llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), llr_i%wfd%keyvloc(llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), &
     llr_j%wfd%keygloc(1,llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), llr_j%wfd%keyvloc(llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), &
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
     llr_i%wfd%keygloc, llr_i%wfd%keyvloc, llr_j%wfd%keygloc, llr_j%wfd%keyvloc, &
     olr%wfd%keygloc, olr%wfd%keyvloc, &
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
     llr_i%wfd%keygloc(1,llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), llr_i%wfd%keyvloc(llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), &
     llr_j%wfd%keygloc(1,llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), llr_j%wfd%keyvloc(llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), &
     olr%wfd%keygloc(1,olr%wfd%nseg_c+min(1,olr%wfd%nseg_f)), olr%wfd%keyvloc(olr%wfd%nseg_c+min(1,olr%wfd%nseg_f)), &
     olr%wfd%nvctr_f)



end subroutine determine_overlapdescriptors_from_descriptors


subroutine check_overlapregion(glr, llr_i, llr_j, olr)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in) :: glr, llr_i, llr_j, olr

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
integer,intent(in) :: n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g
integer,intent(in) :: ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g
integer :: nseg_i, nseg_j
integer,dimension(2,nseg_i),intent(in) :: keyg_i
integer,dimension(nseg_i),intent(in) :: keyv_i
integer,dimension(2,nseg_j),intent(in) :: keyg_j
integer,dimension(nseg_j),intent(in) :: keyv_j
integer,intent(out) :: n1_k, n2_k, n3_k, ns1_k, ns2_k, ns3_k, nseg_k

! Local variables
integer :: iseg, jseg, istart, jstart, istartg, jstartg, kstartg
integer :: iend, jend, iendg, jendg, kendg, transform_index
integer :: kxs, kys, kzs, kxe, kye, kze, kxemax, kyemax, kzemax
character(len=1) :: increase


! Initialize the return values such that they represent a box with no volume
ns1_k=ns1_g+n1_g+1
ns2_k=ns2_g+n2_g+1
ns3_k=ns3_g+n3_g+1
n1_k=0
n2_k=0
n3_k=0
nseg_k=0


! Quick return if possible
if(nseg_i==0 .or. nseg_j==0) return

! Initialize some counters
iseg=min(1,nseg_i)
jseg=min(1,nseg_j)
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


! check if Llrs overlap from there descriptors
! The periodicity is hidden in the fact that we are using the keyglobs
! which are correctly defined. 
subroutine check_overlap_from_descriptors_periodic(nseg_i, nseg_j, keyg_i, keyg_j,  &
           isoverlap, onseg)
use module_base
use module_types
implicit none
! Calling arguments
integer :: nseg_i, nseg_j
integer,dimension(2,nseg_i),intent(in) :: keyg_i
integer,dimension(2,nseg_j),intent(in) :: keyg_j
logical,intent(out) :: isoverlap
integer, intent(out) :: onseg
! Local variables
integer :: iseg, jseg, istart, jstart, kstartg
integer :: iend, jend, kendg, nseg_k


! Initialize some counters
iseg=1
jseg=1
nseg_k=0
isoverlap = .false.
onseg = 0  ! in case they don't overlap
segment_loop: do

    ! Starting point already in global coordinates
    istart=keyg_i(1,iseg)
    jstart=keyg_j(1,jseg)

    ! Ending point already in global coordinates
    iend=keyg_i(2,iseg)
    jend=keyg_j(2,jseg)
    ! Determine starting and ending point of the common segment in global coordinates.
    kstartg=max(istart,jstart)
    kendg=min(iend,jend)

    ! Check whether this common segment has a non-zero length
    if(kendg-kstartg+1>0) then
        isoverlap = .true.
        nseg_k=nseg_k+1
    end if

    ! Check whether all segments of both localization regions have been processed.
    if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

    ! Increase the segment index
    if((iend<=jend .and. iseg<nseg_i) .or. jseg==nseg_j) then
        iseg=iseg+1
    else if(jseg<nseg_j) then
        jseg=jseg+1
    end if

end do segment_loop

if(isoverlap) then
   onseg = nseg_k
end if

end subroutine check_overlap_from_descriptors_periodic

subroutine get_overlap_from_descriptors_periodic(nseg_i, nseg_j, keyg_i, keyg_j,  &                                                                                                                          
           isoverlap,onseg, onvctr, keyglob, keyvglob)
use module_base
use module_types
implicit none
! Calling arguments
integer :: nseg_i, nseg_j
integer,dimension(2,nseg_i),intent(in) :: keyg_i
integer,dimension(2,nseg_j),intent(in) :: keyg_j
logical, intent(in) :: isoverlap
integer, intent(in) :: onseg
integer, intent(out) :: onvctr
integer, dimension(2,max(onseg,1)),intent(out) :: keyglob
integer, dimension(max(onseg,1)), intent(out) :: keyvglob
! Local variables
integer :: iseg, jseg, knvctr, istart, jstart, kstartg
integer :: iend, jend, kendg, nseg_k

if(.not. isoverlap) then
!initialize the variable
onvctr = 0
return
end if

nseg_k = 0
knvctr = 1
iseg = 1
jseg= 1
segment_loop: do

    ! Starting point already in global coordinates
    istart=keyg_i(1,iseg)
    jstart=keyg_j(1,jseg)

    ! Ending point already in global coordinates
    iend=keyg_i(2,iseg)
    jend=keyg_j(2,jseg)

    ! Determine starting and ending point of the common segment in global coordinates.
    kstartg=max(istart,jstart)
    kendg=min(iend,jend)

    ! Check whether this common segment has a non-zero length
    ! and assign the keys
    if(kendg-kstartg+1>0) then
        nseg_k=nseg_k+1
        keyglob(1,nseg_k) = kstartg
        keyglob(2,nseg_k) = kendg
        keyvglob(nseg_k) =  knvctr
        knvctr = knvctr + kendg-kstartg+1
    end if

    ! Check whether all segments of both localization regions have been processed.
    if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

    ! Increase the segment index
    if((iend<=jend .and. iseg<nseg_i) .or. jseg==nseg_j) then
        iseg=iseg+1
    else if(jseg<nseg_j) then
        jseg=jseg+1
    end if

end do segment_loop

onvctr = knvctr-1
!check if everything matches
if(onseg .ne. nseg_k) then
  print *,'onseg ',onseg,' nseg_k ',nseg_k
  stop 'get_overlap_from_descriptors_periodic: number of segments not right'
end if

end subroutine get_overlap_from_descriptors_periodic
!> Creates the overlap descriptor from two input overlap descriptors. The dimension of the overlap box (i.e.
!! ns1_k,ns2_k,ns3_k (starting indices) and n1_k,n2_k,n3_k (length) have to be determined earlier (using
!! the subroutine overlapbox_from_descriptors)).
!! Calling arguments: *_i refers to overlap region i (input)
!!                    *_j refers to overlap region j (input)
!!                    *_g refers to the global region (input)
!!                    *_k refers to the overlap region (input/output)
!> SHOULD NOW CHANGE THIS BECAUSE GLOBAL COORDINATES ARE NOW ALWAYS AVAILABLE
subroutine overlapdescriptors_from_descriptors(n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, &
           ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, ns1_k, ns2_k, ns3_k, &
           nseg_i, nseg_j, nseg_k, &
           keyg_i, keyv_i, keyg_j, keyv_j, keyg_k, keyv_k, nvctr_k)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k
integer,intent(in) :: ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, ns1_k, ns2_k, ns3_k
integer,intent(in) :: nseg_i, nseg_j, nseg_k
integer,dimension(2,nseg_i),intent(in) :: keyg_i
integer,dimension(nseg_i),intent(in) :: keyv_i
integer,dimension(2,nseg_j),intent(in) :: keyg_j
integer,dimension(nseg_j),intent(in) :: keyv_j
integer,dimension(2,nseg_k),intent(out) :: keyg_k
integer,dimension(nseg_k),intent(out) :: keyv_k
integer,intent(out) :: nvctr_k

! Local variables
integer :: iseg, jseg, kseg, knvctr, istart, jstart, kstart, istartg, jstartg, kstartg
integer :: iend, jend, kend, iendg, jendg, kendg, transform_index
character(len=1) :: increase

! Initialize some counters
iseg=min(1,nseg_i)
jseg=min(1,nseg_j)
kseg=0
knvctr=0
nvctr_k=0

! Quick return if possible
if(nseg_i==0 .or. nseg_j==0) return

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
!! Calling arguments:
!!   ist: index to transform (with respcet to coordinate system of locreg A)
!!   n1a, n2a, n3a:             box size of locreg A
!!   n1b, n2b, n3b:             box size of locreg B
!!   nshift1, nshift2, nshift3: nsa-nsb, where nsa,nsb are the starting points of the boxes of A,B (for all 3 dimensions)
function transform_index(ist, n1a, n2a, n3a, n1b, n2b, n3b, nshift1, nshift2, nshift3)
implicit none

! Calling arguments
integer,intent(in) :: ist, n1a, n2a, n3a, n1b, n2b, n3b, nshift1, nshift2, nshift3
integer :: transform_index

! Local variables
integer :: ii, ix, iy, iz, ixg, iyg, izg, istg

  ! Get the coordinates with respect to localization region A
  ii = ist - 1
  iz = ii / ((n1a+1) * (n2a+1))
  ii = ii - iz * ((n1a+1) * (n2a+1))
  iy = ii / (n1a+1)
  ix = ii - iy * (n1a+1)

  ! Transform ix, iy, iz to the coordinates with respect to localization region B
  ! WARNING: MODULO ONLY CORRECT WHEN COMPARING TO GLR
  ! FOR NOW THIS FUNCTION ALWAYS DOES THIS... SHOULD DISAPEAR NOW THE KEYGLOB IS DEFINED
!!  izg = modulo(iz + nshift3, n3b)
!!  iyg = modulo(iy + nshift2, n2b)
!!  ixg = modulo(ix + nshift1, n1b)
  izg = iz + nshift3
  iyg = iy + nshift2
  ixg = ix + nshift1

  ! Transform ist to its counterpart in the coordinate system of B
  istg = izg*(n1b+1)*(n2b+1) + iyg*(n1b+1) + ixg + 1
  
  ! Assign istg to the value that is passed back.
  transform_index=istg

end function transform_index



!> Get the coordinates of ist with respect to its localization region
!! Calling arguments:
!!  ist          index for which coordinates shall be calculated
!!  n1, n2, n3   box sizes
!!  ix, iy, iz   coordinates of ist
subroutine get_coordinates(ist, n1, n2, n3, ix, iy, iz)
implicit none

! Calling arguments
integer,intent(in) :: ist, n1, n2, n3
integer,intent(out) :: ix, iy, iz

! Local variable
integer :: ii

  ! Get the coordinates ix, iy, iz
  ii = ist - 1
  !ii = ist
  iz = ii / ((n1+1) * (n2+1))
  ii = ii - iz * ((n1+1) * (n2+1))
  iy = ii / (n1+1)
  ix = ii - iy * (n1+1)

end subroutine get_coordinates

subroutine check_overlap(Llr_i, Llr_j, Glr, overlap)
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
type(locreg_descriptors),intent(in) :: Llr_i, Llr_j, Glr
logical, intent(out) :: overlap

! Local variables
integer :: onseg

  call check_overlap_cubic_periodic(Glr,Llr_i,Llr_j,overlap)
  if(overlap) then
    call check_overlap_from_descriptors_periodic(Llr_i%wfd%nseg_c, Llr_j%wfd%nseg_c,&
         Llr_i%wfd%keyglob, Llr_j%wfd%keyglob, overlap, onseg)
  end if

end subroutine check_overlap


subroutine transform_keyglob_to_keygloc(Glr,Llr,nseg,keyglob,keygloc)
use module_base
use module_types
use module_interfaces
implicit none
type(locreg_descriptors),intent(in) :: Glr, Llr
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

subroutine transform_ISFcoordinates(direction,iin1,iin2,iin3,Glr,Llr,iout1,iout2,iout3,ishift1, ishift2, ishift3)
use module_types
implicit none
integer, intent(in) :: direction                      !< integer specifying the direction of the transformation (0 from local to global, 1 from global to local)
integer, intent(in) :: iin1, iin2, iin3               !< the input coordinate
type(locreg_descriptors), intent(in) :: Glr           !< global region descriptors
type(locreg_descriptors), intent(in) :: Llr           !< Localization regiondescriptors
integer, intent(out) :: iout1, iout2, iout3           !< ouput coordinates
integer, intent(out) :: ishift1, ishift2, ishift3     !< shift between coordinates: in_coord + shift = out_coord

!This routines supposes that the specified coordinate is part of both regions.
if(direction==0) then
!from local to global
iout1 = iin1 -(Glr%nsi1 - Llr%nsi1)
ishift1 = Llr%nsi1 - Glr%nsi1
iout2 = iin2 -(Glr%nsi2 - Llr%nsi2)
ishift2 = Llr%nsi2 - Glr%nsi2
iout3 = iin3 -(Glr%nsi3 - Llr%nsi3)
ishift3 = Llr%nsi3 - Glr%nsi3
if(iout1 > Glr%nsi1 + Glr%d%n1i .and. Glr%geocode/='F') then
   iout1 = modulo(iout1,Glr%d%n1i+1)+Glr%nsi1
   ishift1 = ishift1 - Glr%d%n1i
end if
if(iout2 > Glr%nsi2 + Glr%d%n2i .and. Glr%geocode=='P') then
   iout2 = modulo(iout2,Glr%d%n2i+1)+Glr%nsi2
   ishift2 = ishift2 - Glr%d%n2i
end if
if(iout3 > Glr%nsi3 + Glr%d%n3i .and. Glr%geocode/='F') then
   iout3 = modulo(iout3,Glr%d%n3i+1)+Glr%nsi3
   ishift3 = ishift3 - Glr%d%n3i
end if

else if(direction==1) then
!from global to local
iout1 = iin1 -(Llr%nsi1 - Glr%nsi1)
ishift1 = Glr%nsi1 - Llr%nsi1
iout2 = iin2 -(Llr%nsi2 - Glr%nsi2)
ishift2 = Glr%nsi2 - Llr%nsi2
iout3 = iin3 -(Llr%nsi3 - Glr%nsi3)
ishift3 = Glr%nsi3 - Llr%nsi3
if(iout1 < 0 .and. Glr%geocode/='F') then
   iout1 = iin1 + Glr%d%n1i - Llr%nsi1
   ishift1 = ishift1 + Glr%d%n1i 
end if
if(iout2 < 0 .and. Glr%geocode=='P') then
   iout2 = iin2 + Glr%d%n2i - Llr%nsi2
   ishift2 = ishift2 + Glr%d%n2i 
end if
if(iout3 < 0 .and. Glr%geocode/='F') then
   iout3 = iin3 + Glr%d%n3i - Llr%nsi3
   ishift3 = ishift3 + Glr%d%n3i 
end if

else
STOP 'transform_ISFcoordinates: wrong descriptor, should be 0 or 1'
end if

end subroutine transform_ISFcoordinates

! Almost degenerate with get_number_of_overlap_region
! should merge the two... prefering this one since argument list is better 
subroutine check_overlap_cubic_periodic(Glr,Ilr,Jlr,isoverlap)
use module_types
use module_base
implicit none
type(locreg_descriptors), intent(in) :: Glr
type(locreg_descriptors), intent(in) :: Ilr
type(locreg_descriptors), intent(in) :: Jlr
logical, intent(out) :: isoverlap
!Local variables
integer :: azones,bzones,ii,izones,jzones !, i_stat, i_all
logical :: go1, go2, go3
integer,dimension(3,8) :: astart,bstart,aend,bend

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Ilr%outofzone(ii) > 0) azones = azones * 2
     if(Jlr%outofzone(ii) > 0) bzones = bzones * 2
  end do

!allocate astart and aend
!!  allocate(astart(3,azones),stat=i_stat)
!!  call memocc(i_stat,astart,'astart',subname)
!!  allocate(aend(3,azones),stat=i_stat)
!!  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Ilr,Ilr%outofzone,astart,aend)

!allocate bstart and bend
!!  allocate(bstart(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bstart,'bstart',subname)
!!  allocate(bend(3,bzones),stat=i_stat)
!!  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Jlr,Jlr%outofzone,bstart,bend)

! Now check if they overlap
  isoverlap = .false.
  loop_izones: do izones=1,azones
    do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
      if(go1 .and. go2 .and. go3) then
        isoverlap = .true.
        exit loop_izones
      end if
    end do
  end do loop_izones

! Deallocation block
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

end subroutine check_overlap_cubic_periodic


subroutine fracture_periodic_zone_ISF(nzones,Glr,Llr,outofzone,astart,aend)

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

! Start and end of Global region
  Gstart(1) = Glr%nsi1
  Gstart(2) = Glr%nsi2
  Gstart(3) = Glr%nsi3
  Gend(1) = Glr%nsi1 + Glr%d%n1i
  Gend(2) = Glr%nsi2 + Glr%d%n2i
  Gend(3) = Glr%nsi3 + Glr%d%n3i

! Periodicity of the system
  period(1) = Glr%d%n1i + 1
  period(2) = Glr%d%n2i + 1
  period(3) = Glr%d%n3i + 1

! Start and end of local region
  alrs(1) = Llr%nsi1
  alrs(2) = Llr%nsi2
  alrs(3) = Llr%nsi3
  alre(1) = Llr%nsi1 + Llr%d%n1i
  alre(2) = Llr%nsi2 + Llr%d%n2i
  alre(3) = Llr%nsi3 + Llr%d%n3i

  if(outofzone(1) <= 0 .and. outofzone(2) <= 0 .and. outofzone(3) <= 0 ) then !Nothing to do
     astart(1,1) = alrs(1)
     astart(2,1) = alrs(2)
     astart(3,1) = alrs(3)
     aend(1,1) = alre(1)
     aend(2,1) = alre(2)
     aend(3,1) = alre(3)
     return
  end if

!assign the first zone (necessarily without shift) and initiliaze the others
  do ii=1,3
     astart(ii,:) = alrs(ii)
     aend(ii,:) = min(Gend(ii),alre(ii))
  end do

!assign the other zones
  index = 2
  do ii=1,2
     if(outofzone(ii) > 0) then    !Translation: X,Y
        astart(ii,index) =  Gstart(ii)
        aend(ii,index) = modulo(alre(ii),period(ii))
        index = index + 1
     end if
     do jj=ii+1,2
        if(outofzone(ii) > 0 .and. outofzone(jj) > 0) then  !Translation: X+Y
           astart(ii,index) = Gstart(ii)
           astart(jj,index) = Gstart(jj)
           aend(ii,index) = modulo(alre(ii),period(ii))
           aend(jj,index) = modulo(alre(jj),period(jj))
           index = index + 1
        end if
     end do
  end do

END SUBROUTINE fracture_periodic_zone_ISF

