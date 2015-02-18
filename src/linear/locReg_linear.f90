!> @file
!!  Routines used by the linear scaling version
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg 
!! taking into account the pediodicity
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

  call f_routine(id='determine_wfd_periodicity')

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
   call allocate_wfd(Llr(ilr)%wfd)

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

   call f_release_routine()

END SUBROUTINE determine_wfd_periodicity


subroutine determine_locregSphere_parallel(iproc,nproc,nlr,hx,hy,hz,astruct,orbs,Glr,Llr,calculateBounds)!,outofzone)

  use module_base
  use module_types
  use module_interfaces, except_this_one => determine_locregSphere_parallel
  use communications, only: communicate_locreg_descriptors_basics, communicate_locreg_descriptors_keys

  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(atomic_structure),intent(in) :: astruct
  type(orbitals_data),intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  type(locreg_descriptors), dimension(nlr), intent(inout) :: Llr
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
  integer :: ii, iall, istat, iorb, iat, norb, norbu, norbd, nspin, jproc, iiorb
  integer,dimension(3) :: outofzone
  integer,dimension(:),allocatable :: rootarr, norbsperatom, norbsperlocreg, onwhichmpi
  real(8),dimension(:,:),allocatable :: locregCenter
  type(orbitals_data) :: orbsder

  call f_routine(id='determine_locregSphere_parallel')

  rootarr = f_malloc(nlr,id='rootarr')

  ! Determine how many locregs one process handles at most
  ii=ceiling(dble(nlr)/dble(nproc))
  !determine the limits of the different localisation regions
  rootarr=1000000000

  onwhichmpi=f_malloc(nlr,id='onwhichmpi')
  iiorb=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc,0)
          iiorb=iiorb+1
          onWhichMPI(iiorb)=jproc
      end do
  end do

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


     if(calculateBounds(ilr)) then 
         ! This makes sure that each locreg is only handled once by one specific processor.
    
         ! Determine the extrema of this localization regions (using only the coarse part, since this is always larger or equal than the fine part).
         call determine_boxbounds_sphere(glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, hx, hy, hz, &
              llr(ilr)%locrad, llr(ilr)%locregCenter, &
              glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyvloc, isx, isy, isz, iex, iey, iez)
         write(*,'(a,3i7)') 'ilr, isx, iex', ilr, isx, iex
    
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
               !isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
               !iex= ln1 + isx
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
               !isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               !iez= ln3 + isz
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
               !isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
               !iex= ln1 + isx
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
               !isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
               !iey= ln2 + isy
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
               !isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               !iez= ln3 + isz
               if (iez > Glr%ns3+Glr%d%n3) then
                  outofzone(3)=modulo(iez,Glr%d%n3+1)
               end if 
            end if
            if(xperiodic .and. yperiodic .and. zperiodic ) then
              Llr(ilr)%geocode = 'P'
            end if
         end select

         ! Make sure that the localization regions are not periodic
         if (xperiodic .or. yperiodic .or. zperiodic) then
             call f_err_throw('The size of the localization region '&
                 &//trim(yaml_toa(ilr,fmt='(i0)'))//&
                 &' is larger than that of the global region.&
                 & Reduce the localization radii or use the cubic version',&
                 & err_name='BIGDFT_RUNTIME_ERROR')
         end if
    
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
         !write(*,*) 'ilr, Llr(ilr)%nsi3',ilr, Llr(ilr)%nsi3

    
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


         ! Make sure that the extent of the interpolating functions grid for the
         ! locreg is not larger than the that of the global box.
         if (Llr(ilr)%d%n1i>Glr%d%n1i) then
             call f_err_throw('The interpolating functions grid in x dimension for locreg '&
                 &//trim(yaml_toa(ilr,fmt='(i0)'))//'('//trim(yaml_toa(Llr(ilr)%d%n1i,fmt='(i0)'))//')&
                 & is larger than that of the global region('//trim(yaml_toa(Glr%d%n1i,fmt='(i0)'))//').&
                 & Reduce the localization radii or use the cubic version',&
                 & err_name='BIGDFT_RUNTIME_ERROR')
         end if
         if (Llr(ilr)%d%n2i>Glr%d%n2i) then
             call f_err_throw('The interpolating functions grid in y dimension for locreg '&
                 &//trim(yaml_toa(ilr,fmt='(i0)'))//'('//trim(yaml_toa(Llr(ilr)%d%n2i,fmt='(i0)'))//')&
                 & is larger than that of the global region('//trim(yaml_toa(Glr%d%n2i,fmt='(i0)'))//').&
                 & Reduce the localization radii or use the cubic version',&
                 & err_name='BIGDFT_RUNTIME_ERROR')
         end if
         if (Llr(ilr)%d%n3i>Glr%d%n3i) then
             call f_err_throw('The interpolating functions grid in z dimension for locreg '&
                 &//trim(yaml_toa(ilr,fmt='(i0)'))//'('//trim(yaml_toa(Llr(ilr)%d%n3i,fmt='(i0)'))//')&
                 & is larger than that of the global region('//trim(yaml_toa(Glr%d%n3i,fmt='(i0)'))//').&
                 & Reduce the localization radii or use the cubic version',&
                 & err_name='BIGDFT_RUNTIME_ERROR')
         end if
    
    !DEBUG
         if (iproc == 0) then
            write(*,*)'Description of zone:',ilr
            write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
            write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
            write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
            write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
            write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
            write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
            write(*,*)'outofzone',ilr,':',outofzone(:)
         end if
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
     call mpiallred(rootarr(1), nlr, mpi_min, bigdft_mpi%mpi_comm)
     
     ! Communicate those parts of the locregs that all processes need.
     call communicate_locreg_descriptors_basics(iproc, nlr, rootarr, orbs, llr)

     !do ilr=1,nlr
     !    write(*,*) 'iproc, nseg_c', iproc, llr(ilr)%wfd%nseg_c
     !end do

     ! Now communicate those parts of the locreg that only some processes need (the keys).
     ! For this we first need to create orbsder that describes the derivatives.
     !call create_orbsder()

     !iiorb=0
     !onwhichmpider=f_malloc(orbsder%norb,id='onwhichmpider')
     !do jproc=0,nproc-1
     !   do iorb=1,orbsder%norb_par(jproc,0)
     !     iiorb=iiorb+1
     !     onWhichMPIder(iiorb)=jproc
     !   end do
     !end do

     ! Now communicate the keys
     call communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, rootarr, onwhichmpi)

     !call deallocate_orbitals_data(orbsder)
     !call f_free(onwhichmpider)
  end if
  call timing(iproc,'comm_llr      ','OF')

  !create the bound arrays for the locregs we need on the MPI tasks
  call timing(iproc,'calc_bounds   ','ON') 
  do ilr=1,nlr
         if (Llr(ilr)%geocode=='F' .and. calculateBounds(ilr) ) then
            write(*,*) 'calling locreg_bounds, ilr', ilr
            call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
                 Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
                 Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
         end if
  end do

  call timing(iproc,'calc_bounds   ','OF') 

  call f_free(rootarr)

  call f_free(onwhichmpi)
  call f_release_routine()


!!$  subroutine create_orbsder()
!!$    call nullify_orbitals_data(orbsder)
!!$    norbsperatom = f_malloc0(astruct%nat,id='norbsperatom')
!!$    locregCenter = f_malloc((/ 3, nlr /),id='locregCenter')
!!$    norbsPerLocreg = f_malloc(nlr,id='norbsPerLocreg')
!!$    do iorb=1,orbs%norb
!!$        iat=orbs%onwhichatom(iorb)
!!$        norbsperatom(iat)=norbsperatom(iat)+3
!!$    end do
!!$    norb=3*orbs%norb
!!$    norbu=norb
!!$    norbd=0
!!$    nspin=1
!!$    call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, nspin, orbs%nspinor,&
!!$         orbs%nkpts, orbs%kpts, orbs%kwgts, orbsder,.true.) !simple repartition
!!$    call f_free_ptr(orbsder%onwhichatom)
!!$
!!$    do ilr=1,nlr
!!$        locregCenter(:,ilr)=llr(ilr)%locregCenter
!!$    end do
!!$                 
!!$    call assignToLocreg2(iproc, nproc, orbsder%norb, orbsder%norb_par, astruct%nat, astruct%nat, &
!!$         nspin, norbsPerAtom, locregCenter, orbsder%onwhichatom)
!!$
!!$    call f_free_ptr(orbsder%inWhichLocreg)
!!$    norbsPerLocreg=3
!!$
!!$    call assignToLocreg2(iproc, nproc, orbsder%norb, orbsder%norb_par, astruct%nat, nlr, &
!!$         nspin, norbsPerLocreg, locregCenter, orbsder%inwhichlocreg)
!!$
!!$    call f_free(locregCenter)
!!$    call f_free(norbsPerLocreg)
!!$    call f_free(norbsperatom)
!!$  end subroutine create_orbsder

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
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period,Gics,Gice
  character(len=*), parameter :: subname='determine_wfdSphere'
!!  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  integer, allocatable :: keygloc_tmp(:,:)

   call f_routine(id=subname)

   !starting point of locreg (can be outside the simulation box)
   isdir(1) = Llr(ilr)%ns1
   isdir(2) = Llr(ilr)%ns2
   isdir(3) = Llr(ilr)%ns3
   !ending point of locreg (can be outside the simulation box)
   iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
   iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
   iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
   ! starting and ending point of coarse grid in Global region
   Gics(1) = Glr%ns1
   Gics(2) = Glr%ns2
   Gics(3) = Glr%ns3
   Gice(1) = Glr%ns1 + Glr%d%n1
   Gice(2) = Glr%ns2 + Glr%d%n2
   Gice(3) = Glr%ns3 + Glr%d%n3
   ! starting and ending point of fine grid in Global region
   Gifs(1) = Glr%d%nfl1 + Glr%ns1
   Gifs(2) = Glr%d%nfl2 + Glr%ns2
   Gifs(3) = Glr%d%nfl3 + Glr%ns3
   Gife(1) = Glr%d%nfu1 + Glr%ns1
   Gife(2) = Glr%d%nfu2 + Glr%ns2
   Gife(3) = Glr%d%nfu3 + Glr%ns3
   ! periodicity
   period(1) = Glr%d%n1+1
   period(2) = Glr%d%n2+1
   period(3) = Glr%d%n3+1

   !!! Determine starting point of the fine grid in locreg
   !!do ii=1,3
   !!   if (Llr(ilr)%outofzone(ii) > 0) then
   !!      ! When periodicity, we must check for 2 different situations:
   !!      ! (1) : starting of locreg before or in fine grid zone
   !!      if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
   !!      ! (2) : starting point after fine grid
   !!      if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
   !!   else
   !!       Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
   !!   end if 
   !!end do

   !!! Determine ending point of the fine grid in locreg
   !!do ii=1,3
   !!   if(Llr(ilr)%outofzone(ii) > 0) then
   !!      !When periodicity, we must check for three different situations:
   !!      ! (1) : ending of locreg before fine grid zone
   !!      if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
   !!      ! (2) : ending of locreg in fine grid zone
   !!      if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
   !!        Life(ii) = iedir(ii)-isdir(ii)
   !!      end if
   !!      ! (3) : ending of locreg after ending of fine grid zone
   !!      if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
   !!   else
   !!      Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
   !!   end if
   !!end do

   do ii=1,3
       ! Determine starting point of the fine grid in locreg. There are two possibilities:
       if (isdir(ii)<gics(ii)) then
           ! Start of the locreg locreg outside of the global box
           lifs(ii) = max(isdir(ii)+period(ii),gifs(ii)) - period(ii) - isdir(ii)
       else if(isdir(ii)>=gics(ii)) then
           ! Start of locreg inside of the global box
           lifs(ii) = max(isdir(ii),gifs(ii)) - isdir(ii)
       else
           stop 'cannot determine start of fine grid'
       end if

       ! Determine ending point of the fine grid in locreg. There are two possibilities:
       if (iedir(ii)>gice(ii)) then
           ! End of the locreg outside of the global box
           life(ii) = min(iedir(ii)-period(ii),gife(ii)) + period(ii) - isdir(ii)
       else if(iedir(ii)<=gice(ii)) then
           ! End of the locreg inside of the global box
           life(ii) = min(iedir(ii),gife(ii)) - isdir(ii)
       else
           stop 'cannot determine start of fine grid'
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
        hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_c, Glr%wfd%keygloc(1,1), &
        Glr%wfd%keyvloc(1), &
        llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nvctr_c)

   !fine part
   call num_segkeys_sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        glr%wfd%nseg_f, Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
        llr(ilr)%wfd%nseg_f, llr(ilr)%wfd%nvctr_f)

   write(*,'(a,2i8)') 'llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f', llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f

   !allocate the wavefunction descriptors following the needs
   call allocate_wfd(Llr(ilr)%wfd)

   !Now, fill the descriptors:

   keygloc_tmp = f_malloc((/ 2, (llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f) /),id='keygloc_tmp')

   !$omp parallel default(private) &
   !$omp shared(Glr,llr,hx,hy,hz,ilr,keygloc_tmp)  
   !$omp sections
   !$omp section

   !coarse part
   call segkeys_Sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        llr(ilr)%wfd%nseg_c, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_c, Glr%wfd%keygloc(1,1), &
        Glr%wfd%keyvloc(1), llr(ilr)%wfd%nvctr_c, &
        llr(ilr)%wfd%keygloc(1,1),llr(ilr)%wfd%keyglob(1,1), &
        llr(ilr)%wfd%keyvloc(1), llr(ilr)%wfd%keyvglob(1), &
        keygloc_tmp(1,1))

   !$omp section
   !fine part
   call segkeys_Sphere(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
        glr%ns1, glr%ns2, glr%ns3, &
        llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
        llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
        llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
        llr(ilr)%wfd%nseg_f, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
        Glr%wfd%nseg_f, Glr%wfd%keygloc(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), llr(ilr)%wfd%nvctr_f, &
        llr(ilr)%wfd%keygloc(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyglob(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyvloc(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        llr(ilr)%wfd%keyvglob(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)), &
        keygloc_tmp(1,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f)))  
   !$omp end sections
   !$omp end parallel

   call f_free(keygloc_tmp)

   call f_release_routine()


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
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check,n1p1,np

  nvctr_loc=0
  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0

  n1p1=n1+1
  np=n1p1*(n2+1)
  do iseg=1,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
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


subroutine num_segkeys_sphere(n1, n2, n3, nl1glob, nl2glob, nl3glob, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, nseg, nvctr)
  implicit none
  integer, intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nsegglob
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,intent(out) :: nseg, nvctr
  !local variables
  logical :: segment, inside
  integer :: i, i1, i2, i3, nstart, nend, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3, n1p1, np
  integer :: ij1, ij2, ij3, jj1, jj2, jj3
  real(kind=8) :: cut, dx,dy, dz


  nvctr=0
  nstart=0
  nend=0

  cut=locrad**2
  n1p1=n1+1
  np=n1p1*(n2+1)

  !$omp parallel default(none) &
  !$omp shared(nsegglob,keygglob,nl1glob,nl2glob,nl3glob,locregCenter) &
  !$omp shared(hx,hy,hz,cut,n1p1,np,nstart,nvctr,nend, n1, n2, n3) &
  !$omp private(iseg,jj,j0,j1,ii,i3,i2,i0,i1,ii2,ii3,ii1,i,dx,dy,dz,segment) &
  !$omp private(inside, ij1, ij2, ij3, jj1, jj2, jj3)
  segment=.false.
  !$omp do reduction(+:nstart,nvctr,nend)
  do iseg=1,nsegglob
      j0=keygglob(1,iseg)
      j1=keygglob(2,iseg)
      ii=j0-1
      i3=ii/np
      ii=ii-i3*np
      i2=ii/n1p1
      i0=ii-i2*n1p1
      i1=i0+j1-j0

      ii2=i2+nl2glob
      ii3=i3+nl3glob

      !dz=((ii3*hz)-locregCenter(3))**2
      !dy=((ii2*hy)-locregCenter(2))**2
      do i=i0,i1
          ii1=i+nl1glob
          !dx=((ii1*hx)-locregCenter(1))**2
          inside=.false.
          do ij3=-1,1
              jj3=ii3+ij3*(n3+1)
              dz=((jj3*hz)-locregCenter(3))**2
              do ij2=-1,1
                  jj2=ii2+ij2*(n2+1)
                  dy=((jj2*hy)-locregCenter(2))**2
                  do ij1=-1,1
                      jj1=ii1+ij1*(n1+1)
                      dx=((jj1*hx)-locregCenter(1))**2
                      if(dx+dy+dz<=cut) then
                          inside=.true.
                      end if
                  end do
              end do
          end do
          if(inside) then
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
      if(segment) then
          ! Always start a new segment if we come to a new line in y direction.
          nend=nend+1
          segment=.false.
      end if
  end do
  !$omp enddo
  !$omp end parallel

  nseg=nstart

  !check
  if (nend /= nstart) then
     write(*,*) 'nend , nstart',nend,nstart
     stop 'nend <> nstart'
  endif

END SUBROUTINE num_segkeys_sphere


subroutine determine_boxbounds_sphere(n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, hx, hy, hz, locrad, locregCenter, &
           nsegglob, keygglob, keyvglob, ixmin, iymin, izmin, ixmax, iymax, izmax)
  use dynamic_memory
  implicit none
  integer, intent(in) :: n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, nsegglob
  real(kind=8),intent(in) :: hx, hy, hz, locrad
  real(kind=8),dimension(3),intent(in) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,intent(out) :: ixmin, iymin, izmin, ixmax, iymax, izmax
  !local variables
  integer :: i, i1, i2, i3, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3, n1p1, np
  integer :: ij1, ij2 ,ij3, jj1, jj2, jj3
  real(kind=8) :: cut, dx,dy, dz
  !debug
  integer :: iiimin, isegmin

  call f_routine(id='determine_boxbounds_sphere')

  iiimin=0
  isegmin=0

  ! Initialize the return values
  ixmax=0
  iymax=0
  izmax=0
  ixmin=nl1glob+n1glob
  iymin=nl2glob+n2glob
  izmin=nl3glob+n3glob

  cut=locrad**2
  n1p1=n1glob+1
  np=n1p1*(n2glob+1)
  !$omp parallel default(none) &
  !$omp shared(nsegglob,keygglob,n1glob,n2glob,n3glob,nl1glob,nl2glob,nl3glob,locregCenter) &
  !$omp shared(ixmin,iymin,izmin,ixmax,iymax,izmax,hx,hy,hz,cut,n1p1,np) &
  !$omp private(iseg,jj,j0,j1,ii,i3,i2,i0,i1,ii2,ii3,ii1,i,dx,dy,dz,iiimin,isegmin) &
  !$omp private(ij1, ij2, ij3, jj1, jj2, jj3)
  !$omp do reduction(max:ixmax,iymax,izmax) reduction(min:ixmin,iymin,izmin)
  do iseg=1,nsegglob
      j0=keygglob(1,iseg)
      j1=keygglob(2,iseg)
      ii=j0-1
      i3=ii/np
      ii=ii-i3*np
      i2=ii/n1p1
      i0=ii-i2*n1p1
      i1=i0+j1-j0

      ii2=i2+nl2glob
      ii3=i3+nl3glob

      !dz=((ii3*hz)-locregCenter(3))**2
      !dy=((ii2*hy)-locregCenter(2))**2
      do i=i0,i1
          ii1=i+nl1glob
          do ij3=-1,1
              jj3=ii3+ij3*(n3glob+1)
              dz=((jj3*hz)-locregCenter(3))**2
              do ij2=-1,1
                  jj2=ii2+ij2*(n2glob+1)
                  dy=((jj2*hy)-locregCenter(2))**2
                  do ij1=-1,1
                      jj1=ii1+ij1*(n1glob+1)
                      dx=((jj1*hx)-locregCenter(1))**2
                      if(dx+dy+dz<=cut) then
                          ixmax=max(jj1,ixmax)
                          iymax=max(jj2,iymax)
                          izmax=max(jj3,izmax)
                          ixmin=min(jj1,ixmin)
                          iymin=min(jj2,iymin)
                          izmin=min(jj3,izmin)
                      end if
                  end do
              end do
          end do
          !dx=((ii1*hx)-locregCenter(1))**2
          !!dx=((ii1*hx)-locregCenter(1))**2
          !!if(dx+dy+dz<=cut) then
          !!    ixmax=max(ii1,ixmax)
          !!    iymax=max(ii2,iymax)
          !!    izmax=max(ii3,izmax)
          !!    ixmin=min(ii1,ixmin)
          !!    !if(ii1<ixmin) iiimin=j0-1 ; isegmin=iseg
          !!    iymin=min(ii2,iymin)
          !!    izmin=min(ii3,izmin)
          !!end if
      end do
  end do
  !$omp enddo
  !$omp end parallel

  call f_release_routine()

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
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l,n1p1,np,n1lp1,nlp
  integer :: i_stat, i_all
  integer :: ngridp,ngridlob,loc
  integer, allocatable :: keyg_loc(:,:)

  call f_routine('segkeys_periodic')

  !should be initialized
  ngridp=-1000
  ngridlob=-1000

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc


  keyg_loc = f_malloc((/ 2, nseg_loc /),id='keyg_loc')

  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  n1p1=n1+1
  np=n1p1*(n2+1)
  n1lp1=n1l+1
  nlp=n1lp1*(n2l+1)
  do iseg=1,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
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
          ngridp=i3l*nlp + i2l*n1lp1 + i1l+1
          ngridlob = i3 * np + i2 * n1p1 + i + 1

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

 call f_free(keyg_loc)

 call f_release_routine()

END SUBROUTINE segkeys_periodic


subroutine segkeys_Sphere(n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, hx, hy, hz, &
     locrad, locregCenter, &
     nsegglob, keygglob, keyvglob, nvctr_loc, keyg_loc, keyg_glob, keyv_loc, keyv_glob, keygloc)
  use module_base
  implicit none
  integer,intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nl1, nu1, nl2, nu2, nl3, nu3, nseg, nsegglob, nvctr_loc
  real(kind=8) :: hx, hy, hz, locrad
  real(kind=8),dimension(3) :: locregCenter
  integer,dimension(2,nsegglob),intent(in) :: keygglob
  integer,dimension(nsegglob),intent(in) :: keyvglob
  integer,dimension(2,nseg),intent(out) :: keyg_loc, keyg_glob
  integer,dimension(nseg),intent(out) :: keyv_loc, keyv_glob
  integer,dimension(2,nseg),intent(out) :: keygloc !tmp
  !local variables
  character(len=*),parameter :: subname = 'segkeys_Sphere'
  integer :: i, i1, i2, i3, nstart, nend, nvctr, igridpoint, igridglob, iseg, jj, j0, j1, ii, i0, n1l, n2l, n3l
  integer :: i1l, i2l, i3l, ii1, ii2, ii3, istat, iall, loc, n1p1, np, n1lp1, nlp, igridpointa, igridgloba
  integer :: ij1, ij2, ij3, jj1, jj2, jj3, ii1mod, ii2mod, ii3mod, ivctr, jvctr
  real(kind=8) :: cut, dx, dy, dz
  logical :: segment, inside
  !integer, allocatable :: keygloc(:,:)

  !call f_routine('segkeys_Sphere')

  !keygloc = f_malloc((/ 2, nseg /),id='keygloc')

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

  !can add openmp here too as segment always ends at end of y direction? 
  !problem is need nend value - can do a pre-scan to find seg value only as with init_collcom.
  !for now just do omp section
  cut=locrad**2
  n1p1=n1+1
  np=n1p1*(n2+1)
  n1lp1=n1l+1
  nlp=n1lp1*(n2l+1)
  ivctr=0
  jvctr=0
  do iseg=1,nsegglob
      j0=keygglob(1,iseg)
      j1=keygglob(2,iseg)
      ii=j0-1
      i3=ii/np
      ii=ii-i3*np
      i2=ii/n1p1
      i0=ii-i2*n1p1
      i1=i0+j1-j0
      ii2=i2+nl2glob
      ii3=i3+nl3glob
      dz=((ii3*hz)-locregCenter(3))**2
      dy=((ii2*hy)-locregCenter(2))**2
      !i2l=ii2-nl2
      !i3l=ii3-nl3
      !igridpointa=i3l*nlp+i2l*n1lp1+1
      igridgloba=ii3*np+ii2*n1p1+1 
      do i=i0,i1
          ii1=i+nl1glob
          dx=((ii1*hx)-locregCenter(1))**2
          i1l=ii1-nl1
          !igridpoint=igridpointa+i1l
          igridglob=igridgloba+ii1 
          inside=.false.
          do ij3=-1,1
              jj3=ii3+ij3*(n3+1)
              dz=((jj3*hz)-locregCenter(3))**2
              do ij2=-1,1
                  jj2=ii2+ij2*(n2+1)
                  dy=((jj2*hy)-locregCenter(2))**2
                  do ij1=-1,1
                      jj1=ii1+ij1*(n1+1)
                      dx=((jj1*hx)-locregCenter(1))**2
                      if(dx+dy+dz<=cut) then
                          if (inside) stop 'twice inside'
                          inside=.true.
                          ii1mod=jj1
                          ii2mod=jj2
                          ii3mod=jj3
                          i1l=jj1-nl1
                          i2l=jj2-nl2
                          i3l=jj3-nl3
                          igridpoint=i3l*nlp+i2l*n1lp1+i1l+1
                          !write(*,'(a,4i8)') 'i1l, i2l, i3l, igridpoint', i1l, i2l, i3l, igridpoint
                      end if
                  end do
              end do
          end do
          if(inside) then
              ! Check that we are not outside of the locreg region
              !ii1mod=modulo(ii1-1,n1)+1
              !ii2mod=modulo(ii2-1,n2)+1
              !ii3mod=modulo(ii3-1,n3)+1
              ivctr=ivctr+1
              if(ii1mod<nl1) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii1mod=',ii1mod,'<',nl1,'=nl1'
                  stop
              end if
              if(ii2mod<nl2) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii2mod=',ii2mod,'<',nl2,'=nl2'
                  stop
              end if
              if(ii3mod<nl3) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii3mod=',ii3mod,'<',nl3,'=nl3'
                  stop
              end if
              if(ii1mod>nu1) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii1mod=',ii1mod,'>',nu1,'=nu1'
                  stop
              end if
              if(ii2mod>nu2) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii2mod=',ii2mod,'>',nu2,'=nu2'
                  stop
              end if
              if(ii3mod>nu3) then
                  write(*,'(a,i0,a,i0,a)') 'ERROR: ii3mod=',ii3mod,'>',nu3,'=nu3'
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
                  jvctr=jvctr+keygloc(2,nend)-keygloc(1,nend)+1
              end if
          end if
      end do
      if(segment) then
          ! Close the segment
          nend=nend+1
          keygloc(2,nend)=igridpoint
          keyg_glob(2,nend)=igridglob
          segment=.false.
          jvctr=jvctr+keygloc(2,nend)-keygloc(1,nend)+1
      end if
  end do

  ! Some checks
  if (ivctr/=nvctr_loc) then
      write(*,*) 'ivctr, nvctr_loc', ivctr, nvctr_loc
      stop 'ivctr/=nvctr_loc'
  end if

  if (jvctr/=nvctr_loc) then
      write(*,*) 'jvctr, nvctr_loc', jvctr, nvctr_loc
      stop 'jvctr/=nvctr_loc'
  end if

  if (nend /= nstart) then
     write(*,*) 'nend , nstart',nend,nstart
     stop 'nend <> nstart'
  endif
  if (nseg /= nstart) then
     write(*,*) 'nseg , nstart',nseg,nstart
     stop 'nseg <> nstart'
  endif

  ! Now build the keyvloc where we replace the segments in order for the loc
  ivctr=0
  do iseg = 1, nseg
     !sorting the keyg_loc
     loc = minloc(keygloc(1,:),1)
     keyg_loc(1,iseg) = keygloc(1,loc)
     keyg_loc(2,iseg) = keygloc(2,loc)
!    print *,'iseg,keygloc,keyg_loc',iseg,keygloc(1,loc),keygloc(2,loc),keyg_loc(1,iseg),keyg_loc(2,iseg)
     keygloc(1,loc) = maxval(keygloc) + 1
     keyv_loc(iseg) = keyv_glob(loc)
     !write(*,'(a,7i8)') 'iseg,keyglob,keyvglob,keygloc,keyvloc',iseg,keyg_glob(1:2,iseg),keyv_glob(iseg),keyg_loc(1:2,iseg),keyv_loc(iseg)
     ivctr=ivctr+keyg_loc(2,iseg)-keyg_loc(1,iseg)+1
  end do
  !call f_free(keygloc)
  if (ivctr/=nvctr_loc) then
      write(*,*) 'ivctr, nvctr_loc', ivctr, nvctr_loc
      stop 'rearrangement check: ivctr/=nvctr_loc'
  end if

  ! Some checks
  ivctr=0
  write(*,*) 'nlp, n1lp1', nlp, n1lp1
  do iseg=1,nseg
     j0=keyg_loc(1,iseg)
     j1=keyg_loc(2,iseg)
     ii=j0-1
     i3=ii/nlp
     ii=ii-i3*nlp
     i2=ii/n1lp1
     i0=ii-i2*n1lp1
     i1=i0+j1-j0
     !if (i2<nl2) then
     !    write(*,'(a,2(i0,a))') 'ERROR: i2=',i2,'<',nl2,'=nl2' ; stop
     !end if
     if (i2>n2l) then
         write(*,'(a,2(i0,a))') 'ERROR: i2=',i2,'>',n2l,'=n2l' ; stop
     end if
     !if (i3<nl3) then
     !    write(*,'(a,2(i0,a))') 'ERROR: i3=',i3,'<',nl3,'=nl3' ; stop
     !end if
     if (i3>n3l) then
         write(*,'(a,2(i0,a))') 'ERROR: i3=',i3,'>',n3l,'=n3l' ; stop
     end if
     do i=i0,i1
        ivctr=ivctr+1
        !write(*,'(a,6i8)') 'j0, j1, ii, i, i2, i3', j0, j1, ii, i, i2, i3
        !if (i<nl1) then
        !    write(*,'(a,2(i0,a))') 'ERROR: i=',i,'<',nl1,'=nl1' ; stop
        !end if
        if (i>n1l) then
            write(*,'(a,2(i0,a))') 'ERROR: i=',i,'>',n1l,'=n1l' ; stop
        end if
     end do
  end do

  if (ivctr/=nvctr_loc) then
      write(*,*) 'ivctr, nvctr_loc', ivctr, nvctr_loc
      stop 'second check: ivctr/=nvctr_loc'
  end if

  !call f_release_routine()

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
  ! Check whether all segments of both localization regions have been processed.
  if(iseg>=nseg_i .and. jseg>=nseg_j) return

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


  subroutine check_overlap(Llr_i, Llr_j, Glr, overlap)
  use locregs, only: locreg_descriptors
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
  integer :: i, j, j0, ii, iz, iy, ix, n1p1, np

  call f_routine(id='transform_keyglob_to_keygloc')

  n1p1=Glr%d%n1+1
  np=n1p1*(Glr%d%n2+1)
  do i = 1 , 2
     do j = 1, nseg
        ! Writing keyglob in cartesian coordinates
        j0 = keyglob(i,j)
        ii = j0-1
        iz = ii/np
        ii = ii-iz*np
        iy = ii/n1p1
        ix = ii-iy*n1p1

        ! Checking consistency
        if(iz < Llr%ns3 .or. iy < Llr%ns2 .or. ix < Llr%ns1) stop 'transform_keyglob_to_keygloc : minimum overflow'
        if(iz > Llr%ns3+Llr%d%n3 .or. iy > Llr%ns2+Llr%d%n2 .or. ix > Llr%ns1+Llr%d%n1)&
           stop 'transform_keyglob_to_keygloc : maximum overflow'

        ! Using coordinates to write keygloc      
        keygloc(i,j) = (iz-Llr%ns3)*(Llr%d%n1+1)*(Llr%d%n2+1) + (iy-Llr%ns2)*(Llr%d%n1+1) + (ix-Llr%ns1) + 1
     end do
  end do

  call f_release_routine()

end subroutine transform_keyglob_to_keygloc


!> Almost degenerate with get_number_of_overlap_region
!! should merge the two... prefering this one since argument list is better 
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

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Ilr,Ilr%outofzone,astart,aend)

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

end subroutine check_overlap_cubic_periodic



!> Tranform wavefunction between localisation region and the global region
!!!!!#######!> This routine only works if both locregs have free boundary conditions.
!! @warning 
!! WARNING: Make sure psi is set to zero where Glr does not collide with Llr (or everywhere)
subroutine Lpsi_to_global2(iproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, lpsi, psi)

  use module_base
  use module_types

 implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in):: iproc
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
  !integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
  character(len=*), parameter :: subname='Lpsi_to_global'
  integer :: i_all
  integer :: start,Gstart,Lindex
  integer :: lfinc,Gfinc,spinshift,ispin,Gindex,isegstart
  integer :: istart
  !integer :: i_stat

  call f_routine(id=subname)

  !!! This routine is only intended for conversions between locregs with the same boundary conditions.
  !!if (glr%geocode/= 'F' .or. llr%geocode/='F') then
  !!    call f_err_throw('Lpsi_to_global2 can only be used for locregs with free boundary conditions', &
  !!         err_name='BIGDFT_RUNTIME_ERROR')
  !!end if

  if(nspin/=1) stop 'not fully implemented for nspin/=1!'

! Define integers
  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
  icheck = 0
  spinshift = Gdim / nspin
 
! Get the keymask: shift for every segment of Llr (with respect to Glr)
! allocate(keymask(2,nseg),stat=i_stat)
  !keymask = f_malloc((/2,nseg/),id='keymask')

  !call shift_locreg_indexes(Glr,Llr,keymask,nseg)
  !call shift_locreg_indexes_global(Glr,Llr,keymask,nseg)
  !!keymask = llr%wfd%keyglob

!####################################################
! Do coarse region
!####################################################
  isegstart=1

 
  !$omp parallel default(private) &
  !$omp shared(Glr,Llr, lpsi,icheck,psi,norb) &
  !$omp firstprivate(isegstart,nseg,lincrement,Gincrement,spinshift,nspin) 

  !$omp do reduction(+:icheck)
  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
     lmin = llr%wfd%keyglob(1,isegloc)
     lmax = llr%wfd%keyglob(2,isegloc)
     istart = llr%wfd%keyvglob(isegloc)-1

     
     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
        Gmin = Glr%wfd%keyglob(1,isegG)
        Gmax = Glr%wfd%keyglob(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        !if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax) exit global_loop_c

        !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_c
        if(lmin > Gmax)  cycle global_loop_c

        ! Define the offset between the two segments
        offset = lmin - Gmin
        if(offset < 0) then
           offset = 0
        end if

        ! Define the length of the two segments
        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new global wavefunction
        icheck = icheck + (length + 1)

        ! WARNING: index goes from 0 to length because it is the offset of the element

        do ix = 0,length     
           istart = istart + 1
           do ispin=1,nspin
              Gindex = Glr%wfd%keyvglob(isegG)+offset+ix+spinshift*(ispin-1)
              Lindex = istart+lincrement*norb*(ispin-1)
              psi(Gindex) = lpsi(Lindex) 
           end do
        end do
     end do global_loop_c
  end do local_loop_c
  !$omp end do


!##############################################################
! Now do fine region
!##############################################################

  start = Llr%wfd%nvctr_c
  Gstart = Glr%wfd%nvctr_c
  lfinc  = Llr%wfd%nvctr_f
  Gfinc = Glr%wfd%nvctr_f

  isegstart=Glr%wfd%nseg_c+1

  !$omp do reduction(+:icheck)
  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
     lmin = llr%wfd%keyglob(1,isegloc)
     lmax = llr%wfd%keyglob(2,isegloc)
     istart = llr%wfd%keyvglob(isegloc)-1

     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f

        Gmin = Glr%wfd%keyglob(1,isegG)
        Gmax = Glr%wfd%keyglob(2,isegG)

        ! For each segment in Llr check if there is a collision with the segment in Glr
        ! if not, cycle
        if(lmin > Gmax) then
            isegstart=isegG
        end if
        if(Gmin > lmax)  exit global_loop_f
        !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f
        if(lmin > Gmax)  cycle global_loop_f

        offset = lmin - Gmin
        if(offset < 0) offset = 0

        length = min(lmax,Gmax)-max(lmin,Gmin)

        !Find the common elements and write them to the new global wavefunction
        ! First set to zero those elements which are not copied. WARNING: will not work for npsin>1!!
 
        icheck = icheck + (length + 1)

        ! WARNING: index goes from 0 to length because it is the offset of the element
        do ix = 0,length
        istart = istart + 1
           do igrid=1,7
              do ispin = 1, nspin
                 Gindex = Gstart + (Glr%wfd%keyvglob(isegG)+offset+ix-1)*7+igrid + spinshift*(ispin-1)
                 Lindex = start+(istart-1)*7+igrid + lincrement*norb*(ispin-1) 
                 psi(Gindex) = lpsi(Lindex) 
              end do
           end do
        end do
     end do global_loop_f
  end do local_loop_f
  !$omp end do

  !$omp end parallel

  !Check if the number of elements in loc_psi is valid
  if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
    write(*,*)'There is an error in Lpsi_to_global2: sum of fine and coarse points used',icheck
    write(*,*)'is not equal to the sum of fine and coarse points in the region',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
    stop
  end if

  !!call f_free(keymask)

  call f_release_routine()

END SUBROUTINE Lpsi_to_global2
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
