module locregs_init
  implicit none

  private

  !> Public routines
  public :: determine_locregsphere_parallel

  contains


    subroutine determine_locregSphere_parallel(iproc,nproc,nlr,hx,hy,hz,astruct,orbs,Glr,Llr,calculateBounds)!,outofzone)
    
      use module_base
      use module_types
      !use module_interfaces, except_this_one => determine_locregSphere_parallel
      use communications, only: communicate_locreg_descriptors_basics, communicate_locreg_descriptors_keys
      use bounds, only: locreg_bounds 
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
      logical :: perx, pery, perz
    
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
    
      ! Periodicity in the three directions
      Gperx=(Glr%geocode /= 'F')
      Gpery=(Glr%geocode == 'P')
      Gperz=(Glr%geocode /= 'F')
    
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
             call determine_boxbounds_sphere(gperx, gpery, gperz, glr%d%n1, glr%d%n2, glr%d%n3, glr%ns1, glr%ns2, glr%ns3, &
                  hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
                  glr%wfd%nseg_c, glr%wfd%keygloc, glr%wfd%keyvloc, isx, isy, isz, iex, iey, iez)
             !write(*,'(a,3i7)') 'ilr, isx, iex', ilr, isx, iex
        
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
        
        !!DEBUG
        !     if (iproc == 0) then
        !        write(*,*)'Description of zone:',ilr
        !        write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
        !        write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
        !        write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
        !        write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
        !        write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
        !        write(*,*)'nsi:',Llr(ilr)%nsi1,Llr(ilr)%nsi2,Llr(ilr)%nsi3
        !        write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
        !        write(*,*)'outofzone',ilr,':',outofzone(:)
        !     end if
        !!DEBUG
        
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
         call mpiallred(rootarr(1), nlr, mpi_min, comm=bigdft_mpi%mpi_comm)
         
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
                !write(*,*) 'calling locreg_bounds, ilr', ilr
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



    subroutine determine_boxbounds_sphere(gperx, gpery, gperz, n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, &
               hx, hy, hz, locrad, locregCenter, &
               nsegglob, keygglob, keyvglob, ixmin, iymin, izmin, ixmax, iymax, izmax)
      use dynamic_memory
      implicit none
      logical,intent(in) :: gperx, gpery, gperz
      integer, intent(in) :: n1glob, n2glob, n3glob, nl1glob, nl2glob, nl3glob, nsegglob
      real(kind=8),intent(in) :: hx, hy, hz, locrad
      real(kind=8),dimension(3),intent(in) :: locregCenter
      integer,dimension(2,nsegglob),intent(in) :: keygglob
      integer,dimension(nsegglob),intent(in) :: keyvglob
      integer,intent(out) :: ixmin, iymin, izmin, ixmax, iymax, izmax
      !local variables
      integer :: i, i1, i2, i3, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3, n1p1, np
      integer :: ij1, ij2 ,ij3, jj1, jj2, jj3
      integer :: ijs1, ije1, ijs2, ije2, ijs3, ije3
      real(kind=8) :: cut, dx,dy, dz
      !debug
      integer :: iiimin, isegmin
    
      call f_routine(id='determine_boxbounds_sphere')
    
      ! For perdiodic boundary conditions, one has to check also in the neighboring
      ! cells (see in the loop below)
      if (gperx) then
          ijs1 = -1
          ije1 = 1
      else
          ijs1 = 0
          ije1 = 0
      end if
      if (gpery) then
          ijs2 = -1
          ije2 = 1
      else
          ijs2 = 0
          ije2 = 0
      end if
      if (gperz) then
          ijs3 = -1
          ije3 = 1
      else
          ijs3 = 0
          ije3 = 0
      end if
    
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
      !$omp shared(ixmin,iymin,izmin,ixmax,iymax,izmax,hx,hy,hz,cut,n1p1,np,ijs1,ije1,ijs2,ije2,ijs3,ije3) &
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
              do ij3=ijs3,ije3!-1,1
                  jj3=ii3+ij3*(n3glob+1)
                  dz=((jj3*hz)-locregCenter(3))**2
                  do ij2=ijs2,ije2!-1,1
                      jj2=ii2+ij2*(n2glob+1)
                      dy=((jj2*hy)-locregCenter(2))**2
                      do ij1=ijs1,ije1!-1,1
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


    subroutine num_segkeys_sphere(perx, pery, perz, n1, n2, n3, nl1glob, nl2glob, nl3glob, hx, hy, hz, &
         locrad, locregCenter, &
         nsegglob, keygglob, keyvglob, nseg, nvctr)
      implicit none
      logical,intent(in) :: perx, pery, perz 
      integer, intent(in) :: n1, n2, n3, nl1glob, nl2glob, nl3glob, nsegglob
      real(kind=8),intent(in) :: hx, hy, hz, locrad
      real(kind=8),dimension(3),intent(in) :: locregCenter
      integer,dimension(2,nsegglob),intent(in) :: keygglob
      integer,dimension(nsegglob),intent(in) :: keyvglob
      integer,intent(out) :: nseg, nvctr
      !local variables
      logical :: segment, inside
      integer :: i, i1, i2, i3, nstart, nend, iseg, jj, j0, j1, ii, i0, ii1, ii2, ii3, n1p1, np
      integer :: ij1, ij2, ij3, jj1, jj2, jj3, ijs1, ijs2, ijs3, ije1, ije2, ije3
      real(kind=8) :: cut, dx,dy, dz
    
    
      nvctr=0
      nstart=0
      nend=0
    
      cut=locrad**2
      n1p1=n1+1
      np=n1p1*(n2+1)
    
      ! For perdiodic boundary conditions, one has to check also in the neighboring
      ! cells (see in the loop below)
      if (perx) then
          ijs1 = -1
          ije1 = 1
      else
          ijs1 = 0
          ije1 = 0
      end if
      if (pery) then
          ijs2 = -1
          ije2 = 1
      else
          ijs2 = 0
          ije2 = 0
      end if
      if (perz) then
          ijs3 = -1
          ije3 = 1
      else
          ijs3 = 0
          ije3 = 0
      end if
    
      !$omp parallel default(none) &
      !$omp shared(nsegglob,keygglob,nl1glob,nl2glob,nl3glob,locregCenter) &
      !$omp shared(hx,hy,hz,cut,n1p1,np,nstart,nvctr,nend, n1, n2, n3, ijs1, ijs2, ijs3, ije1, ije2, ije3) &
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
              do ij3=ijs3,ije3!-1,1
                  jj3=ii3+ij3*(n3+1)
                  dz=((jj3*hz)-locregCenter(3))**2
                  do ij2=ijs2,ije2!-1,1
                      jj2=ii2+ij2*(n2+1)
                      dy=((jj2*hy)-locregCenter(2))**2
                      do ij1=ijs1,ije1!-1,1
                          jj1=ii1+ij1*(n1+1)
                          dx=((jj1*hx)-locregCenter(1))**2
                          !write(*,'(a,6i7,4es12.4)') 'ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut', ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut
                          if(dx+dy+dz<=cut) then
                              !write(*,'(a,6i7,4es12.4)') 'ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut', ii1, ii2, ii3, jj1, jj2, jj3, dx, dy, dz, cut
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


    subroutine segkeys_Sphere(perx, pery, perz, n1, n2, n3, nl1glob, nl2glob, nl3glob, &
         nl1, nu1, nl2, nu2, nl3, nu3, nseg, hx, hy, hz, &
         locrad, locregCenter, &
         nsegglob, keygglob, keyvglob, nvctr_loc, keyg_loc, keyg_glob, keyv_loc, keyv_glob, keygloc)
      use module_base
      use sparsematrix_init, only: distribute_on_threads
      implicit none
      logical,intent(in) :: perx, pery, perz
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
      integer :: ij1, ij2, ij3, jj1, jj2, jj3, ii1mod, ii2mod, ii3mod, ivctr, jvctr, kvctr, ijs1, ijs2, ijs3, ije1, ije2, ije3
      real(kind=8) :: cut, dx, dy, dz
      logical :: segment, inside
      integer,dimension(:,:),pointer :: ise
      integer :: ithread, jthread, nthread, ivctr_tot, jvctr_tot, nstart_tot, nend_tot, kthread, j, offset
      integer,dimension(:),allocatable :: nstartarr, keyv_last
      integer,dimension(:,:,:),allocatable :: keygloc_work, keyg_glob_work
      integer,dimension(:,:),allocatable :: keyv_glob_work
      !$ integer :: omp_get_thread_num
      !integer, allocatable :: keygloc(:,:)
    
      call f_routine('segkeys_Sphere')
    
      !keygloc = f_malloc((/ 2, nseg /),id='keygloc')
    
      !dimensions of the localisation region (O:nIl)
      ! must be smaller or equal to simulation box dimensions
      !n1l=i1ec-i1sc
      !n2l=i2ec-i2sc
      !n3l=i3ec-i3sc
      n1l=nu1-nl1
      n2l=nu2-nl2
      n3l=nu3-nl3
    
    
      ! For perdiodic boundary conditions, one has to check also in the neighboring
      ! cells (see in the loop below)
      if (perx) then
          ijs1 = -1
          ije1 = 1
      else
          ijs1 = 0
          ije1 = 0
      end if
      if (pery) then
          ijs2 = -1
          ije2 = 1
      else
          ijs2 = 0
          ije2 = 0
      end if
      if (perz) then
          ijs3 = -1
          ije3 = 1
      else
          ijs3 = 0
          ije3 = 0
      end if
    
      call distribute_on_threads(nsegglob, nthread, ise)
    
      keygloc_work = f_malloc((/1.to.2,1.to.nseg,0.to.nthread-1/),id='keygloc_work')
      keyg_glob_work = f_malloc((/1.to.2,1.to.nseg,0.to.nthread-1/),id='keyg_glob_work')
      keyv_glob_work = f_malloc((/1.to.nseg,0.to.nthread-1/),id='keyv_glob_work')
      keyv_last = f_malloc(0.to.nthread-1,id='keyv_last')
    
    
      nstartarr = f_malloc(0.to.nthread-1,id='nstartarr')
    
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
      kvctr=0
      nvctr=0
      nstart=0
      nend=0
      ivctr_tot = 0
      jvctr_tot = 0
      nstart_tot = 0
      nend_tot = 0
      segment=.false.
      ithread = 0
      !$omp parallel &
      !$omp default(none) &
      !$omp shared(ise, hx, hy, hz, keygglob, np, n1p1, nl1glob, nl2glob, nl3glob, locregCenter) &
      !$omp shared(keygloc_work, keyg_glob_work, keyv_glob_work, nstartarr, nl1, nl2, nl3, nu1, nu2, nu3) &
      !$omp shared(ijs3, ije3, ijs2, ije2, ijs1, ije1, n1, n2, n3, cut, n1lp1, nlp, nthread) &
      !$omp shared(keygloc, keyg_glob, keyv_glob, ivctr_tot, jvctr_tot, nstart_tot, nend_tot, keyv_last) &
      !$omp firstprivate(ithread, ivctr, jvctr, kvctr, nvctr, nstart, nend, segment) &
      !$omp private(iseg, j0, j1, ii, i3, i2, i1, i0, ii2, ii3, dz, dy, igridgloba, jj1) &
      !$omp private(i, ii1, dx, i1l, igridglob, inside, ij3, jj3, ij2, jj2, ij1, i2l, i3l) &
      !$omp private(ii1mod, ii2mod, ii3mod, igridpoint, offset, j, kthread)
      !jj1, )
      !$ ithread = omp_get_thread_num()
      do iseg=ise(1,ithread),ise(2,ithread)
      !do iseg=1,nsegglob
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
              do ij3=ijs3,ije3!-1,1
                  jj3=ii3+ij3*(n3+1)
                  dz=((jj3*hz)-locregCenter(3))**2
                  do ij2=ijs2,ije2!-1,1
                      jj2=ii2+ij2*(n2+1)
                      dy=((jj2*hy)-locregCenter(2))**2
                      do ij1=ijs1,ije1!-1,1
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
              !write(*,*) 'ii1, ii2, ii3, inside', ii1, ii2, ii3, inside
              if(inside) then
                  ! Check that we are not outside of the locreg region
                  ivctr=ivctr+1
                  kvctr=kvctr+1
                  !write(*,*) 'inside: kvctr, igridpoint', kvctr, igridpoint
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
                      keygloc_work(1,nstart,ithread)=igridpoint
                      keyg_glob_work(1,nstart,ithread)=igridglob
                      keyv_glob_work(nstart,ithread)=nvctr
                      segment=.true.
                  end if
              else
                  if(segment) then
                      nend=nend+1
                      keygloc_work(2,nend,ithread)=igridpoint!-1
                      keyg_glob_work(2,nend,ithread)=igridglob-1
                      !write(*,'(a,4i7)') 'outside: kvctr, igridpoint, keygloc(1:2,nend)', kvctr, igridpoint, keygloc(1:2,nend)
                      segment=.false.
                      jvctr=jvctr+keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
                      if (kvctr/=keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1) then
                          write(*,*) 'kvctr, keygloc(2,nend)-keygloc(1,nend)+1', &
                               kvctr, keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
                          stop 'kvctr/=keygloc(2,nend)-keygloc(1,nend)+1'
                      end if
                      kvctr=0
                  end if
              end if
          end do
          if(segment) then
              ! Close the segment
              nend=nend+1
              keygloc_work(2,nend,ithread)=igridpoint
              keyg_glob_work(2,nend,ithread)=igridglob
              segment=.false.
              jvctr=jvctr+keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
              if (kvctr/=keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1) then
                  write(*,*) 'kvctr, keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1', &
                      kvctr, keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1
                  stop 'kvctr/=keygloc_work(2,nend,ithread)-keygloc_work(1,nend,ithread)+1'
              end if
              kvctr=0
          end if
      end do
      ! Some checks
      if (nstart/=nend) call f_err_throw('nstart/=nend',err_name='BIGDFT_RUNTIME_ERROR')
      ! Number of segments calculated by ithread
      nstartarr(ithread) = nstart
      ! Number of elements calculated by ithread
      if (nstart>0) then
          keyv_last(ithread) = keyv_glob_work(nstart,ithread)+keyg_glob_work(2,nstart,ithread)-keyg_glob_work(1,nstart,ithread)
      else
          keyv_last(ithread) = 0
      end if
      !$omp barrier
      ii = 1
      do jthread=0,nthread-1
          if (ithread==jthread) then
              if (nstartarr(jthread)>0) then
                  call f_memcpy(n=2*nstartarr(jthread), src=keygloc_work(1,1,ithread), dest=keygloc(1,ii))
                  call f_memcpy(n=2*nstartarr(jthread), src=keyg_glob_work(1,1,ithread), dest=keyg_glob(1,ii))
                  offset = 0
                  do kthread=0,jthread-1
                      offset = offset + keyv_last(kthread)
                  end do
                  do j=1,nstartarr(jthread)
                      keyv_glob(ii+j-1) = keyv_glob_work(j,ithread) + offset
                  end do
                  !call f_memcpy(n=nstartarr(jthread), src=keyv_glob_work(1,ithread), dest=keyv_glob(ii))
              end if
          end if
          ii = ii + nstartarr(jthread)
      end do
    
      !$omp critical
          ivctr_tot = ivctr_tot + ivctr
          jvctr_tot = jvctr_tot + jvctr
          nstart_tot = nstart_tot + nstart
          nend_tot = nend_tot + nend
          !nseg_tot = nseg_tot + nseg
      !$omp end critical
      !$omp end parallel
    
      !write(*,*) 'nstartarr',nstartarr
      !do ii=1,nseg
      !    write(*,*) 'ii, keygloc(:,ii)', ii, keygloc(:,ii)
      !    write(*,*) 'ii, keyg_glob(:,ii)', ii, keyg_glob(:,ii)
      !    write(*,*) 'ii, keyv_glob(ii)', ii, keyv_glob(ii)
      !end do
    
    
      ! Some checks
      if (ivctr_tot/=nvctr_loc) then
          write(*,*) 'ivctr_tot, nvctr_loc', ivctr_tot, nvctr_loc
          stop 'ivctr_tot/=nvctr_loc'
      end if
    
      if (jvctr_tot/=nvctr_loc) then
          write(*,*) 'jvctr_tot, nvctr_loc', jvctr_tot, nvctr_loc
          stop 'jvctr_tot/=nvctr_loc'
      end if
    
      if (nend_tot /= nstart_tot) then
         write(*,*) 'nend_tot , nstart_tot',nend_tot,nstart_tot
         stop 'nend_tot <> nstart_tot'
      endif
      if (nseg /= nstart_tot) then
         write(*,*) 'nseg , nstart_tot',nseg,nstart_tot
         stop 'nseg <> nstart_tot'
      endif
    
      ! Now build the keyvloc where we replace the segments in order for the loc
      ivctr=0
      ii = maxval(keygloc)
      do iseg=1,nseg
         !sorting the keyg_loc
         loc = minloc(keygloc(1,:),1)
         keyg_loc(1,iseg) = keygloc(1,loc)
         keyg_loc(2,iseg) = keygloc(2,loc)
    !    print *,'iseg,keygloc,keyg_loc',iseg,keygloc(1,loc),keygloc(2,loc),keyg_loc(1,iseg),keyg_loc(2,iseg)
         keyv_loc(iseg) = keyv_glob(loc)
         !keygloc(1,loc) = maxval(keygloc) + 1
         keygloc(1,loc) = ii+iseg !just put to the maximal value
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
      !write(*,*) 'nlp, n1lp1', nlp, n1lp1
      !$omp parallel &
      !$omp default(none) &
      !$omp shared(nseg, keyg_loc, nlp, n1lp1, n1l, n2l, n3l, ivctr) &
      !$omp private(iseg, j0, j1, ii, i3, i2, i1, i0, i)
      !$omp do reduction(+:ivctr)
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
      !$omp end do
      !$omp end parallel
    
      if (ivctr/=nvctr_loc) then
          write(*,*) 'ivctr, nvctr_loc', ivctr, nvctr_loc
          stop 'second check: ivctr/=nvctr_loc'
      end if
    
      call f_free(keygloc_work)
      call f_free(keyg_glob_work)
      call f_free(keyv_glob_work)
      call f_free(keyv_last)
      call f_free(nstartarr)
      call f_free_ptr(ise)
    
      call f_release_routine()
    
    END SUBROUTINE segkeys_Sphere


    !> Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg 
    !! taking into account the pediodicity
    !!          
    !! @warning
    !!    We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
    subroutine determine_wfdSphere(ilr,nlr,Glr,hx,hy,hz,Llr)!,outofzone)
    
      use module_base
      use module_types
      use locregs, only: allocate_wfd
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
      logical :: perx, pery, perz
    
       call f_routine(id=subname)
    
       ! periodicity in the three directions
       perx=(glr%geocode /= 'F')
       pery=(glr%geocode == 'P')
       perz=(glr%geocode /= 'F')
    
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
       call num_segkeys_sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
            glr%ns1, glr%ns2, glr%ns3, &
            hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
            Glr%wfd%nseg_c, Glr%wfd%keygloc(1:,1:), &
            Glr%wfd%keyvloc(1:), &
            llr(ilr)%wfd%nseg_c, llr(ilr)%wfd%nvctr_c)
    
       !fine part
       call num_segkeys_sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
            glr%ns1, glr%ns2, glr%ns3, &
            hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
            glr%wfd%nseg_f, Glr%wfd%keygloc(1:,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):), &
            Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):), &
            llr(ilr)%wfd%nseg_f, llr(ilr)%wfd%nvctr_f)
    
       !write(*,'(a,2i8)') 'llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f', llr(ilr)%wfd%nvctr_c, llr(ilr)%wfd%nvctr_f
    
       !allocate the wavefunction descriptors following the needs
       call allocate_wfd(Llr(ilr)%wfd)
    
       !Now, fill the descriptors:
    
       keygloc_tmp = f_malloc((/ 2, (llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f) /),id='keygloc_tmp')
    
       !!$omp parallel default(private) &
       !!$omp shared(Glr,llr,hx,hy,hz,ilr,keygloc_tmp,perx,pery,perz)  
       !!$omp sections
       !!$omp section
    
       !coarse part
       call segkeys_Sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
            glr%ns1, glr%ns2, glr%ns3, &
            llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
            llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
            llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
            llr(ilr)%wfd%nseg_c, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
            Glr%wfd%nseg_c, Glr%wfd%keygloc(1:,1:), &
            Glr%wfd%keyvloc(1:), llr(ilr)%wfd%nvctr_c, &
            llr(ilr)%wfd%keygloc(1:,1:),llr(ilr)%wfd%keyglob(1:,1:), &
            llr(ilr)%wfd%keyvloc(1:), llr(ilr)%wfd%keyvglob(1:), &
            keygloc_tmp(1:,1:))
    
       !!$omp section
       !fine part
       call segkeys_Sphere(perx, pery, perz, Glr%d%n1, Glr%d%n2, Glr%d%n3, &
            glr%ns1, glr%ns2, glr%ns3, &
            llr(ilr)%ns1, llr(ilr)%ns1+llr(ilr)%d%n1, &
            llr(ilr)%ns2, llr(ilr)%ns2+llr(ilr)%d%n2, &
            llr(ilr)%ns3, llr(ilr)%ns3+llr(ilr)%d%n3, &
            llr(ilr)%wfd%nseg_f, hx, hy, hz, llr(ilr)%locrad, llr(ilr)%locregCenter, &
            Glr%wfd%nseg_f, Glr%wfd%keygloc(1:,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):),&
            Glr%wfd%keyvloc(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f):), llr(ilr)%wfd%nvctr_f, &
            llr(ilr)%wfd%keygloc(1:,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
            llr(ilr)%wfd%keyglob(1:,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
            llr(ilr)%wfd%keyvloc(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
            llr(ilr)%wfd%keyvglob(llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):), &
            keygloc_tmp(1:,llr(ilr)%wfd%nseg_c+min(1,llr(ilr)%wfd%nseg_f):))  
       !!$omp end sections
       !!$omp end parallel
    
       call f_free(keygloc_tmp)
    
       call f_release_routine()
    
    
    END SUBROUTINE determine_wfdSphere
    
    
    end module locregs_init
