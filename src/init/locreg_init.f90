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
  use module_base
  use module_types
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
  use module_types
  implicit none
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
  !local variables
  character(len=*), parameter :: subname='draw_locregs'
  character(len=4) :: message
  integer :: i1,i2,i3,ilr,nvctr_tot,i_stat,i_all
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
     allocate(logrid_c(0:Llr(ilr)%d%n1,0:Llr(ilr)%d%n2,0:Llr(ilr)%d%n3+ndebug),stat=i_stat)
     call memocc(i_stat,logrid_c,'logrid_c',subname)
     allocate(logrid_f(0:Llr(ilr)%d%n1,0:Llr(ilr)%d%n2,0:Llr(ilr)%d%n3+ndebug),stat=i_stat)
     call memocc(i_stat,logrid_f,'logrid_f',subname)

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


     i_all=-product(shape(logrid_c))*kind(logrid_c)
     deallocate(logrid_c,stat=i_stat)
     call memocc(i_stat,i_all,'logrid_c',subname)
     i_all=-product(shape(logrid_f))*kind(logrid_f)
     deallocate(logrid_f,stat=i_stat)
     call memocc(i_stat,i_all,'logrid_f',subname)
  end do

  !close file for writing
  close(unit=22)  
END SUBROUTINE draw_locregs


!> Calculates the bounds arrays needed for convolutions
subroutine locreg_bounds(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds)
  use module_base
  use module_types
  use module_interfaces, except_this_one => locreg_bounds
  implicit none
  !Arguments
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(out) :: bounds
  !Local variables
  character(len=*), parameter :: subname='locreg_bounds'
  integer :: i_stat,i_all
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f


  !define logrids
  allocate(logrid_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_c,'logrid_c',subname)
  allocate(logrid_f(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_f,'logrid_f',subname)
  
  call wfd_to_logrids(n1,n2,n3,wfd,logrid_c,logrid_f)

  !allocate and calculate kinetic bounds
  allocate(bounds%kb%ibyz_c(2,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%kb%ibyz_c,'bounds%kb%ibyz_c',subname)
  allocate(bounds%kb%ibxz_c(2,0:n1,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%kb%ibxz_c,'bounds%kb%ibxz_c',subname)
  allocate(bounds%kb%ibxy_c(2,0:n1,0:n2+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%kb%ibxy_c,'bounds%kb%ibxy_c',subname)
  allocate(bounds%kb%ibyz_f(2,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%kb%ibyz_f,'bounds%kb%ibyz_f',subname)
  allocate(bounds%kb%ibxz_f(2,0:n1,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%kb%ibxz_f,'bounds%kb%ibxz_f',subname)
  allocate(bounds%kb%ibxy_f(2,0:n1,0:n2+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%kb%ibxy_f,'bounds%kb%ibxy_f',subname)

  call make_bounds(n1,n2,n3,logrid_c,bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c)
  call make_bounds(n1,n2,n3,logrid_f,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c',subname)
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f',subname)
  
  !allocate grow, shrink and real bounds
  allocate(bounds%gb%ibzxx_c(2,0:n3,-14:2*n1+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%gb%ibzxx_c,'bounds%gb%ibzxx_c',subname)
  allocate(bounds%gb%ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%gb%ibxxyy_c,'bounds%gb%ibxxyy_c',subname)
  allocate(bounds%gb%ibyz_ff(2,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%gb%ibyz_ff,'bounds%gb%ibyz_ff',subname)
  allocate(bounds%gb%ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%gb%ibzxx_f,'bounds%gb%ibzxx_f',subname)
  allocate(bounds%gb%ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%gb%ibxxyy_f,'bounds%gb%ibxxyy_f',subname)

  allocate(bounds%sb%ibzzx_c(2,-14:2*n3+16,0:n1+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%sb%ibzzx_c,'bounds%sb%ibzzx_c',subname)
  allocate(bounds%sb%ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%sb%ibyyzz_c,'bounds%sb%ibyyzz_c',subname)
  allocate(bounds%sb%ibxy_ff(2,nfl1:nfu1,nfl2:nfu2+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%sb%ibxy_ff,'bounds%sb%ibxy_ff',subname)
  allocate(bounds%sb%ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%sb%ibzzx_f,'bounds%sb%ibzzx_f',subname)
  allocate(bounds%sb%ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%sb%ibyyzz_f,'bounds%sb%ibyyzz_f',subname)

  allocate(bounds%ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds%ibyyzz_r,'bounds%ibyyzz_r',subname)

  call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       bounds%kb%ibxy_c,bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
       bounds%kb%ibxy_f,bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
       bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
       bounds%kb%ibyz_f,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,&
       bounds%ibyyzz_r)

END SUBROUTINE locreg_bounds




!>   This subroutine define other wavefunctions descriptors starting from the original descriptors 
!!   and the limits of a given localisation region
!!   it also returns an array which is used to mask the compressed wavefunction into the new one
!! INPUTS
!!   @param ilocreg        localisation region to be considered
!!   @param nlocreg        total number of localisation regions
!!   @param n1,n2          original dimensions of the global box
!!   @param lrlims         array of limits of the localisation regions (global system coordinates)
!!   @param wfdg           global wavefunction descriptors structure
!! OUTPUT
!!   @param wfdl           local wavefunction descriptors structure in local system coordinates
!!   @param keymask        mask array for traducing the wavefunction in compressed form
!!                         to the wavefunction in compressed form for the local system
!!   @param ncountlocreg   array of elements for each localisation region
!!subroutine loc_wfd(ilocreg,nlocreg,n1,n2,lrlims,wfdg,wfdl,keymask,ncountlocreg)
!!  use module_base
!!  use module_types
!!  implicit none
!!  type(wavefunctions_descriptors), intent(in) :: wfdg
!!  integer, intent(in) :: ilocreg,nlocreg,n1,n2
!!  integer, dimension(2,3,nlocreg), intent(in) :: lrlims
!!  type(wavefunctions_descriptors), intent(out) :: wfdl
!!  integer, dimension(nlocreg), intent(out) :: ncountlocreg
!!  integer, dimension(:), pointer :: keymask
!!  !local variables
!!  character(len=*), parameter :: subname='loc_wfd'
!!  integer :: i_stat
!!  integer :: iloc,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nvctr_c,nseg_c,nseg_f,nvctr_f,ndimkey
!!
!!  !calculate the number of segments of the new descriptors for each localisation region
!!  !and the dimension of the array for the translation of the localisation regions
!!  ndimkey=0
!!  do iloc=1,nlocreg
!!     if (iloc /= ilocreg) then
!!        !coarse part
!!        call num_segkeys_loc(n1,n2,lrlims(1,1,iloc),lrlims(2,1,iloc),&
!!             lrlims(1,2,iloc),lrlims(2,2,iloc),lrlims(1,3,iloc),lrlims(2,3,iloc),&
!!             wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
!!             nseg_c,nvctr_c)
!!        !fine part
!!        call num_segkeys_loc(n1,n2,lrlims(1,1,iloc),lrlims(2,1,iloc),&
!!             lrlims(1,2,iloc),lrlims(2,2,iloc),lrlims(1,3,iloc),lrlims(2,3,iloc),&
!!             wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!             wfdg%keyv(wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!             nseg_f,nvctr_f)
!!        ncountlocreg(iloc)=nvctr_c+7*nvctr_f
!!        ndimkey=ndimkey+nseg_c+nseg_f
!!     end if
!!  end do
!!
!!  i1sc=lrlims(1,1,ilocreg)
!!  i1ec=lrlims(2,1,ilocreg)
!!  i2sc=lrlims(1,2,ilocreg)
!!  i2ec=lrlims(2,2,ilocreg)
!!  i3sc=lrlims(1,3,ilocreg)
!!  i3ec=lrlims(2,3,ilocreg)
!!
!!  !coarse part
!!  call num_segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
!!       wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
!!       wfdl%nseg_c,wfdl%nvctr_c)
!!
!!  !fine part
!!  call num_segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
!!       wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdg%keyv(wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdl%nseg_f,wfdl%nvctr_f)
!!
!!  ncountlocreg(ilocreg)=wfdl%nvctr_c+7*wfdl%nvctr_f
!!  ndimkey=ndimkey+wfdl%nseg_c+wfdl%nseg_f
!!
!!  call allocate_wfd(wfdl,subname)
!!
!!  allocate(keymask(ndimkey+ndebug),stat=i_stat)
!!  call memocc(i_stat,keymask,'keymask',subname)
!!
!!  !now fill the local wavefunction descriptors
!!  !and define the mask array for the wavefunction
!!  !coarse part
!!  call segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,& !n(m)
!!       wfdg%nseg_c,wfdg%keyg(1,1),wfdg%keyv(1),&
!!       wfdl%nseg_c,wfdl%nvctr_c,wfdl%keyg(1,1),wfdl%keyv(1))!,keymask(1))
!!
!!  !fine part
!!  call segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,& !n(m)
!!       wfdg%nseg_f,wfdg%keyg(1,wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdg%keyv(wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdl%nseg_f,wfdl%nvctr_f,wfdl%keyg(1,wfdl%nseg_c+min(1,wfdl%nseg_f)),&
!!       wfdl%keyv(wfdl%nseg_c+min(1,wfdl%nseg_f)))!,&
!!       !keymask(wfdg%nseg_c+1))
!!
!!  !a little check on the masking array
!!!!!  if (count(maskarr) /= wfdl%nvctr_c+7*wfdl%nvctr_f) then
!!!!!     write(*,'(1x,a)')'ERROR : Masking problem, check maskarr'
!!!!!     stop
!!!!!  end if
!!
!!END SUBROUTINE loc_wfd


!!subroutine build_keymask(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg_tot,keyg,keyv,&
!!     nseg_loc,keymask)
!!  implicit none
!!  integer, intent(in) :: n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg_tot,nseg_loc
!!  integer, dimension(nseg_tot), intent(in) :: keyv
!!  integer, dimension(2,nseg_tot), intent(in) :: keyg
!!  integer, dimension(2,nseg_loc), intent(out) :: keymask
!!  !local variables
!!  logical :: go,lseg
!!  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend
!!
!!  !start and end points
!!  nsrt=0
!!  nend=0
!!  do iseg=1,nseg_tot
!!     jj=keyv(iseg)
!!     j0=keyg(1,iseg)
!!     j1=keyg(2,iseg)
!!     ii=j0-1
!!     i3=ii/((n1+1)*(n2+1))
!!     ii=ii-i3*(n1+1)*(n2+1)
!!     i2=ii/(n1+1)
!!     i0=ii-i2*(n1+1)
!!     i1=i0+j1-j0
!!     go=(i3sc <= i3 .and. i3 <= i3ec) .and. (i2sc <= i2 .and. i2 <= i2ec)
!!     lseg=.false.
!!     do i=i0,i1
!!        !index of the compressed function
!!        ind=i-i0+jj
!!        if (go .and. (i1sc <= i .and. i <= i1ec)) then
!!           if (.not. lseg) then
!!              nsrt=nsrt+1
!!              keymask(1,nsrt)=ind
!!           end if
!!           lseg=.true.
!!        else
!!           if (lseg) then
!!              keymask(2,nend)=ind-1
!!              nend=nend+1
!!              lseg=.false. 
!!           end if
!!        end if
!!     end do
!!     if (lseg) then
!!        keymask(2,nend)=ind
!!        nend=nend+1
!!     end if
!!  end do
!!
!!  !check
!!  if (nend /= nsrt .or. nend /= nseg_loc) then
!!     write(*,'(1x,a,2(i6))')&
!!          'ERROR: problem in build_keymask',&
!!          nend,nsrt,nseg_loc
!!     stop
!!  end if
!!
!!END SUBROUTINE build_keymask









