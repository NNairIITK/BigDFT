!> @file
!!  Routines to create the localisation region
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Create the localisation region information for cubic code
subroutine create_Glr(geocode,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i,wfd,bounds,Glr)
  use locregs
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(in) :: bounds
  type(locreg_descriptors), intent(out) :: Glr

  Glr%geocode=geocode
  Glr%ns1=0
  Glr%ns2=0
  Glr%ns3=0
  Glr%d%n1=n1
  Glr%d%n2=n2
  Glr%d%n3=n3
  Glr%d%nfl1=nfl1
  Glr%d%nfl2=nfl2
  Glr%d%nfl3=nfl3
  Glr%d%nfu1=nfu1
  Glr%d%nfu2=nfu2
  Glr%d%nfu3=nfu3
  Glr%d%n1i=n1i
  Glr%d%n2i=n2i
  Glr%d%n3i=n3i
  Glr%wfd=wfd !it just associates the pointers
  if (geocode == 'F') then
     Glr%bounds=bounds
  end if
END SUBROUTINE create_Glr


!> Determine a set of localisation regions from the centers and the radii.
!! cut in cubes the global reference system
!!subroutine determine_locreg(nlr,cxyz,locrad,hx,hy,hz,Glr,Llr)
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: nlr
!!  real(gp), intent(in) :: hx,hy,hz
!!  type(locreg_descriptors), intent(in) :: Glr
!!  real(gp), dimension(nlr), intent(in) :: locrad
!!  real(gp), dimension(3,nlr), intent(in) :: cxyz
!!  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
!!  !local variables
!!  character(len=*), parameter :: subname='determine_locreg'
!!  logical :: perx,pery,perz
!!  integer :: ilr,isx,isy,isz,iex,iey,iez
!!  real(gp) :: rx,ry,rz,cutoff
!!
!!  !determine the limits of the different localisation regions
!!  do ilr=1,nlr
!!
!!     rx=cxyz(1,ilr)
!!     ry=cxyz(2,ilr)
!!     rz=cxyz(3,ilr)
!!
!!     cutoff=locrad(ilr)
!!
!!     isx=floor((rx-cutoff)/hx)
!!     isy=floor((ry-cutoff)/hy)
!!     isz=floor((rz-cutoff)/hz)
!!
!!     iex=ceiling((rx+cutoff)/hx)
!!     iey=ceiling((ry+cutoff)/hy)
!!     iez=ceiling((rz+cutoff)/hz)
!!
!!     !assign the geometric code to the localisation region
!!     select case(Glr%geocode)
!!     case('F')
!!        isx=max(isx,0)
!!        isy=max(isy,0)
!!        isz=max(isz,0)
!!        
!!        iex=min(iex,Glr%d%n1)
!!        iey=min(iey,Glr%d%n2)
!!        iez=min(iez,Glr%d%n3)
!!
!!        perx=.false.
!!        pery=.false.
!!        perz=.false.
!!        Llr(ilr)%geocode='F'
!!     case('S')
!!        if (iex - isx >= Glr%d%n1 - 14) then
!!           isx=0
!!           iex=Glr%d%n1
!!           perx=.true.
!!        else
!!           isx=modulo(isx,Glr%d%n1+1)
!!           iex=modulo(iex,Glr%d%n1+1)
!!           perx=.false.
!!        end if
!!
!!        isy=max(isy,0)
!!        iey=min(iey,Glr%d%n2)
!!        pery=.false.
!!
!!        !control the geometric code for the localisation region
!!        if (iez - isz >= Glr%d%n3 - 14) then
!!           isz=0
!!           iez=Glr%d%n3
!!           perz=.true.
!!        else
!!           isz=modulo(isz,Glr%d%n3+1)
!!           iez=modulo(iez,Glr%d%n3+1)
!!           perz=.false.
!!        end if
!!
!!        if (perx .and. perz) then
!!           Llr(ilr)%geocode='S'
!!        else if (.not.perx .and. .not.perz) then
!!           Llr(ilr)%geocode='F'
!!        else
!!           write(*,*)'ERROR: localisation region geometry not allowed'
!!           stop
!!        end if
!!
!!     case('P')
!!        if (iex - isx >= Glr%d%n1 - 14) then
!!           isx=0
!!           iex=Glr%d%n1
!!           perx=.true.
!!        else
!!           isx=modulo(isx,Glr%d%n1+1)
!!           iex=modulo(iex,Glr%d%n1+1)
!!           perx=.false.
!!        end if
!!        if (iey - isy >= Glr%d%n2 - 14) then
!!           isy=0
!!           iey=Glr%d%n2
!!           pery=.true.
!!        else
!!           isy=modulo(isy,Glr%d%n2+1)
!!           iey=modulo(iey,Glr%d%n2+1)
!!           pery=.false.
!!        end if
!!        if (iez - isz >= Glr%d%n3 - 14) then
!!           isz=0
!!           iez=Glr%d%n3
!!           perz=.true.
!!        else
!!           isz=modulo(isz,Glr%d%n3+1)
!!           iez=modulo(iez,Glr%d%n3+1)
!!           perz=.false.
!!        end if
!!        if (perx .and. perz .and. pery) then
!!           Llr(ilr)%geocode='P'
!!        else if (.not.perx .and. .not.perz .and. .not. pery) then
!!           Llr(ilr)%geocode='F'
!!        else if (perx .and. perz .and. .not. pery) then
!!           Llr(ilr)%geocode='S'
!!        else
!!           write(*,*)'ERROR: localisation region geometry not allowed'
!!           stop
!!        end if
!!
!!     end select
!!     
!!     !values for the starting point of the cube
!!     Llr(ilr)%ns1=isx
!!     Llr(ilr)%ns2=isy
!!     Llr(ilr)%ns3=isz
!!     !dimensions of the localisation region
!!     Llr(ilr)%d%n1=iex-isx
!!     Llr(ilr)%d%n2=iey-isy
!!     Llr(ilr)%d%n3=iez-isz
!!
!!     !dimensions of the fine grid inside the localisation region
!!     if (isx < iex) then
!!        Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx
!!        Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
!!     else
!!        write(*,*)'Yet to be implemented (little effort?)'
!!        stop
!!     end if
!!
!!     if (isy < iey) then
!!        Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
!!        Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
!!     else
!!        write(*,*)'Yet to be implemented (little effort?)'
!!        stop
!!     end if
!!
!!     if (isz < iez) then
!!        Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
!!        Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz
!!     else
!!        write(*,*)'Yet to be implemented (little effort?)'
!!        stop
!!     end if
!!     
!!     !dimensions of the interpolating scaling functions grid
!!     select case(Llr(ilr)%geocode)
!!     case('F')
!!        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
!!        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
!!        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31
!!     case('S')
!!        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+2
!!        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
!!        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+2
!!     case('P')
!!        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+2
!!        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+2
!!        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+2
!!     end select
!!
!!     
!!     !define the wavefunction descriptors inside the localisation region
!!     !calculate the number of point and segments for local localisation regions
!!     !coarse part
!!     call num_segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,&
!!          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
!!          Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c)
!!     !fine part
!!     call num_segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,&
!!          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
!!          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f)
!!
!!     !allocate the wavefunction descriptors following the needs
!!     call allocate_wfd(Llr(ilr)%wfd,subname)
!!
!!     !fill such descriptors
!!     !coarse part
!!     call segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,& !n(m)
!!          Glr%wfd%nseg_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
!!          Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
!!          Llr(ilr)%wfd%keyg(1,1),Llr(ilr)%wfd%keyv(1))
!!     !fine part
!!     call segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,& !n(m) 
!!          Glr%wfd%nseg_f,&
!!          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
!!          Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
!!          Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)))
!!
!!     !if the localisation region is isolated build also the bounds
!!     if (Llr(ilr)%geocode=='F') then
!!        call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
!!             Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
!!             Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
!!     end if
!!
!!  end do
!!
!!  !after all localisation regions are determined draw them
!!  !call draw_locregs(nlr,hx,hy,hz,Llr)
!!
!!END SUBROUTINE determine_locreg


subroutine draw_locregs(nlr,hx,hy,hz,Llr)
  use module_base
  use locregs
  implicit none
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
  !local variables
  character(len=*), parameter :: subname='draw_locregs'
  character(len=4) :: message
  integer :: i1,i2,i3,ilr,nvctr_tot
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  !calculate total number
  nvctr_tot=0
  do ilr=1,nlr
     nvctr_tot=nvctr_tot+Llr(ilr)%wfd%nvctr_c+Llr(ilr)%wfd%nvctr_f
  end do

  !open file for writing
  open(unit=22,file='locregs.xyz',status='unknown')
  write(22,*) nvctr_tot,' atomic'
  write(22,*)'coarse and fine points of all the different localisation regions'

  do ilr=1,nlr
     !define logrids
     logrid_c = f_malloc((/ 0.to.Llr(ilr)%d%n1, 0.to.Llr(ilr)%d%n2, 0.to.Llr(ilr)%d%n3 /),id='logrid_c')
     logrid_f = f_malloc((/ 0.to.Llr(ilr)%d%n1, 0.to.Llr(ilr)%d%n2, 0.to.Llr(ilr)%d%n3 /),id='logrid_f')

     call wfd_to_logrids(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,Llr(ilr)%wfd,&
          logrid_c,logrid_f)

     write(message,'(1a,i0)')'g',ilr
     do i3=0,Llr(ilr)%d%n3  
        do i2=0,Llr(ilr)%d%n2  
           do i1=0,Llr(ilr)%d%n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   message,real(i1+Llr(ilr)%ns1,gp)*hx,&
                   real(i2+Llr(ilr)%ns2,gp)*hy,real(i3+Llr(ilr)%ns3,gp)*hz
           enddo
        enddo
     end do
     write(message,'(1a,i0)')'G',ilr
     do i3=0,Llr(ilr)%d%n3 
        do i2=0,Llr(ilr)%d%n2 
           do i1=0,Llr(ilr)%d%n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   message,real(i1+Llr(ilr)%ns1,gp)*hx,&
                   real(i2+Llr(ilr)%ns2,gp)*hy,real(i3+Llr(ilr)%ns3,gp)*hz
           enddo
        enddo
     enddo


     call f_free(logrid_c)
     call f_free(logrid_f)
  end do

  !close file for writing
  close(unit=22)  
END SUBROUTINE draw_locregs


!> Calculates the bounds arrays needed for convolutions
subroutine locreg_bounds(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds)
  use module_base
  use locregs
  !use module_interfaces, except_this_one => locreg_bounds
  implicit none
  !Arguments
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(out) :: bounds
  !Local variables
  character(len=*), parameter :: subname='locreg_bounds'
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  call f_routine(id=subname)

  !define logrids
  logrid_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_c')
  logrid_f = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_f')
  
  call wfd_to_logrids(n1,n2,n3,wfd,logrid_c,logrid_f)

  !allocate and calculate kinetic bounds
  bounds%kb%ibyz_c = f_malloc_ptr((/ 1.to.2,0.to.n2,0.to.n3/),id='bounds%kb%ibyz_c')
  bounds%kb%ibxz_c = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n3/),id='bounds%kb%ibxz_c')
  bounds%kb%ibxy_c = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n2/),id='bounds%kb%ibxy_c')
  bounds%kb%ibyz_f = f_malloc_ptr((/ 1.to.2,0.to.n2,0.to.n3/),id='bounds%kb%ibyz_f')
  bounds%kb%ibxz_f = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n3/),id='bounds%kb%ibxz_f')
  bounds%kb%ibxy_f = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n2/),id='bounds%kb%ibxy_f')

  call make_bounds(n1,n2,n3,logrid_c,bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c)
  call make_bounds(n1,n2,n3,logrid_f,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)

  call f_free(logrid_c)
  call f_free(logrid_f)
  
  !allocate grow, shrink and real bounds
  bounds%gb%ibzxx_c = f_malloc_ptr((/ 1.to.2, 0.to.n3, -14.to.2*n1+16 /),id='bounds%gb%ibzxx_c')
  bounds%gb%ibxxyy_c = f_malloc_ptr((/ 1.to.2, -14.to.2*n1+16, -14.to.2*n2+16 /),id='bounds%gb%ibxxyy_c')
  bounds%gb%ibyz_ff = f_malloc_ptr((/ 1.to.2, nfl2.to.nfu2, nfl3.to.nfu3 /),id='bounds%gb%ibyz_ff')
  bounds%gb%ibzxx_f = f_malloc_ptr((/ 1.to.2, nfl3.to.nfu3, 2*nfl1-14.to.2*nfu1+16 /),id='bounds%gb%ibzxx_f')
  bounds%gb%ibxxyy_f = f_malloc_ptr((/ 1.to.2, 2*nfl1-14.to.2*nfu1+16, 2*nfl2-14.to.2*nfu2+16 /),id='bounds%gb%ibxxyy_f')

  bounds%sb%ibzzx_c = f_malloc_ptr((/ 1.to.2 , -14.to.2*n3+16 , 0.to.n1 /),id='bounds%sb%ibzzx_c')
  bounds%sb%ibyyzz_c = f_malloc_ptr((/ 1.to.2, -14.to.2*n2+16, -14.to.2*n3+16 /),id='bounds%sb%ibyyzz_c')
  bounds%sb%ibxy_ff = f_malloc_ptr((/ 1.to.2, nfl1.to.nfu1, nfl2.to.nfu2 /),id='bounds%sb%ibxy_ff')
  bounds%sb%ibzzx_f = f_malloc_ptr((/ 1.to.2, -14+2*nfl3.to.2*nfu3+16, nfl1.to.nfu1 /),id='bounds%sb%ibzzx_f')
  bounds%sb%ibyyzz_f = f_malloc_ptr((/ 1.to.2, -14+2*nfl2.to.2*nfu2+16, -14+2*nfl3.to.2*nfu3+16 /),id='bounds%sb%ibyyzz_f')

  bounds%ibyyzz_r = f_malloc_ptr((/ 1.to.2,-14.to.2*n2+16,-14.to.2*n3+16/),id='bounds%ibyyzz_r')

  call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       bounds%kb%ibxy_c,bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
       bounds%kb%ibxy_f,bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
       bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
       bounds%kb%ibyz_f,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,&
       bounds%ibyyzz_r)

  call f_release_routine()

END SUBROUTINE locreg_bounds








