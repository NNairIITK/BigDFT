!!****p* BigDFT/memguess
!! NAME
!!   memguess
!!
!! FUNCTION
!!  Test the input files and estimates the memory occupation versus the number
!!  of processors
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2007 CEA
!!
!! SOURCE
!!
program memguess

  use module_types
  !use libBigDFT

  implicit none
  integer, parameter :: ngx=31
  character(len=1) :: geocode
  character(len=20) :: tatonam,units
  logical :: calc_tail,output_grid,optimise
  integer :: ierror,nat,ntypes,nproc,n1,n2,n3,i_stat,i_all
  integer :: nelec,natsc,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  integer :: norb,norbu,norbd,norbe,norbp,norbsc
  integer :: iunit,ityp
  real(kind=8) :: alat1,alat2,alat3,peakmem,hx,hy,hz
  type(input_variables) :: in
  character(len=20), dimension(:), allocatable :: atomnames
  integer, dimension(:), allocatable :: iatype,nelpsp,nzatom,npspcode,iasctype
  integer, dimension(:,:), allocatable :: neleconf
  real(kind=8), dimension(:,:), allocatable :: rxyz,radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psppar, psiat
  real(kind=8), dimension(:,:), allocatable :: xp, occupat
  integer, dimension(:,:), allocatable :: nl
  logical, dimension(:), allocatable :: lfrztyp
  logical, dimension(:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng,norbsc_arr
  real(kind=8), dimension(:), allocatable :: occup,spinar

! Get arguments
  call getarg(1,tatonam)
  
  optimise=.false.
  if(trim(tatonam)=='') then
     write(*,'(1x,a)')&
          'Usage: ./memguess <nproc> [y]'
     write(*,'(1x,a)')&
          'Indicate the number of processes after the executable'
     write(*,'(1x,a)')&
          'You can put a "y" in the second argument (optional) if you want the '
     write(*,'(1x,a)')&
          '  grid to be plotted with V_Sim'
     write(*,'(1x,a)')&
          'You can also put an "o" if you want to rotate the molecule such that'
     write(*,'(1x,a)')&
          '  the volume of the simulation box is optimised'
     stop
  else
     read(unit=tatonam,fmt=*) nproc
     call getarg(2,tatonam)
     if(trim(tatonam)=='') then
        output_grid=.false.
     else if (trim(tatonam)=='y') then
        output_grid=.true.
        write(*,'(1x,a)')&
             'The system grid will be displayed in the "grid.xyz" file'
     else if (trim(tatonam)=='o') then
        optimise=.true.
        output_grid=.true.
        write(*,'(1x,a)')&
             'The optimised system grid will be displayed in the "grid.xyz" file'
     else
        write(*,'(1x,a)')&
             'Usage: ./memguess <nproc> [y]'
        write(*,'(1x,a)')&
             'Indicate the number of processes after the executable'
        write(*,'(1x,a)')&
             'ERROR: The only second argument which is accepted is "y"'
        stop
     end if
  end if

  !initialize memory counting
  call memocc(0,0,'count','start')

  !welcome screen
  call print_logo()

  !read number of atoms
  open(unit=99,file='posinp',status='old')
  read(99,*) nat,units
 
  allocate(rxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz))*kind(rxyz),'rxyz','memguess')
  allocate(iatype(nat),stat=i_stat)
  call memocc(i_stat,product(shape(iatype))*kind(iatype),'iatype','memguess')
  allocate(atomnames(100),stat=i_stat) 
  call memocc(i_stat,product(shape(atomnames))*kind(atomnames),'atomnames','memguess')
  allocate(lfrztyp(nat),stat=i_stat)
  call memocc(i_stat,product(shape(lfrztyp))*kind(lfrztyp),'lfrztyp','memguess')

  !read atomic positions
  call read_atomic_positions(0,99,units,nat,ntypes,iatype,atomnames,lfrztyp,rxyz)

  close(99)

  !new way of reading the input variables, use structures
  call read_input_variables(0,in)

  !hard-code the geocode variable for the moment
  geocode='F'
  call print_input_parameters(in)

  write(*,'(1x,a)')&
       '------------------------------------------------------------------ System Properties'
 
  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(psppar(0:4,0:6,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(psppar))*kind(psppar),'psppar','memguess')
  allocate(nelpsp(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nelpsp))*kind(nelpsp),'nelpsp','memguess')
  allocate(radii_cf(ntypes,2),stat=i_stat)
  call memocc(i_stat,product(shape(radii_cf))*kind(radii_cf),'radii_cf','memguess')
  allocate(npspcode(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(npspcode))*kind(npspcode),'npspcode','memguess')
  allocate(nzatom(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nzatom))*kind(nzatom),'nzatom','memguess')
  allocate(iasctype(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(iasctype))*kind(iasctype),'iasctype','memguess')

  call read_system_variables(0,nproc,nat,ntypes,in%nspin,in%ncharge,in%mpol,in%ixc,in%hgrid,&
       atomnames,iatype,psppar,radii_cf,npspcode,iasctype,nelpsp,nzatom,nelec,natsc,&
       norb,norbu,norbd,norbp,iunit)

! Allocations for the occupation numbers
  allocate(occup(norb),stat=i_stat)
  call memocc(i_stat,product(shape(occup))*kind(occup),'occup','memguess')
  allocate(spinar(norb),stat=i_stat)
  call memocc(i_stat,product(shape(spinar))*kind(spinar),'occup','memguess')
  
! Occupation numbers
  call input_occup(0,iunit,nelec,norb,norbu,norbd,in%nspin,in%mpol,occup,spinar)

! De-allocations
  i_all=-product(shape(occup))*kind(occup)
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup','memguess')
  i_all=-product(shape(spinar))*kind(spinar)
  deallocate(spinar,stat=i_stat)
  call memocc(i_stat,i_all,'spinar','memguess')
  
  !in the case in which the number of orbitals is not "trivial" check whether they are too many
  if ( max(norbu,norbd) /= ceiling(real(nelec,kind=4)/2.0) ) then
     ! Allocations for readAtomicOrbitals (check inguess.dat and psppar files + give norbe)
     allocate(xp(ngx,ntypes),stat=i_stat)
     call memocc(i_stat,product(shape(xp))*kind(xp),'xp','memguess')
     allocate(psiat(ngx,5,ntypes),stat=i_stat)
     call memocc(i_stat,product(shape(psiat))*kind(psiat),'psiat','memguess')
     allocate(occupat(5,ntypes),stat=i_stat)
     call memocc(i_stat,product(shape(occupat))*kind(occupat),'occupat','memguess')
     allocate(ng(ntypes),stat=i_stat)
     call memocc(i_stat,product(shape(ng))*kind(ng),'ng','memguess')
     allocate(nl(4,ntypes),stat=i_stat)
     call memocc(i_stat,product(shape(nl))*kind(nl),'nl','memguess')
     allocate(scorb(4,natsc),stat=i_stat)
     call memocc(i_stat,product(shape(scorb))*kind(scorb),'scorb','memguess')
     allocate(norbsc_arr(natsc+1),stat=i_stat)
     call memocc(i_stat,product(shape(norbsc_arr))*kind(norbsc_arr),'norbsc_arr','memguess')

     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(0,ngx,xp,psiat,occupat,ng,nl,nzatom,nelpsp,psppar,&
          & npspcode,norbe,norbsc,atomnames,ntypes,iatype,iasctype,nat,natsc,scorb,&
          & norbsc_arr)

     ! De-allocations
     i_all=-product(shape(xp))*kind(xp)
     deallocate(xp,stat=i_stat)
     call memocc(i_stat,i_all,'xp','memguess')
     i_all=-product(shape(psiat))*kind(psiat)
     deallocate(psiat,stat=i_stat)
     call memocc(i_stat,i_all,'psiat','memguess')
     i_all=-product(shape(occupat))*kind(occupat)
     deallocate(occupat,stat=i_stat)
     call memocc(i_stat,i_all,'occupat','memguess')
     i_all=-product(shape(ng))*kind(ng)
     deallocate(ng,stat=i_stat)
     call memocc(i_stat,i_all,'ng','memguess')
     i_all=-product(shape(nl))*kind(nl)
     deallocate(nl,stat=i_stat)
     call memocc(i_stat,i_all,'nl','memguess')
     i_all=-product(shape(scorb))*kind(scorb)
     deallocate(scorb,stat=i_stat)
     call memocc(i_stat,i_all,'scorb','memguess')
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr','memguess')

     ! Check the maximum number of orbitals
     if (in%nspin==1) then
        if (norb>norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals (',norb,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
     else
        if (norbu>norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals up (',norbu,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
        if (norbd>norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals down (',norbd,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
     end if

  end if

  i_all=-product(shape(lfrztyp))*kind(lfrztyp)
  deallocate(lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp','BigDFT')
  i_all=-product(shape(psppar))*kind(psppar)
  deallocate(psppar,stat=i_stat)
  call memocc(i_stat,i_all,'psppar','memguess')
  i_all=-product(shape(npspcode))*kind(npspcode)
  deallocate(npspcode,stat=i_stat)
  call memocc(i_stat,i_all,'npspcode','memguess')
  i_all=-product(shape(nelpsp))*kind(nelpsp)
  deallocate(nelpsp,stat=i_stat)
  call memocc(i_stat,i_all,'nelpsp','memguess')
  !no need of using nzatom array
  i_all=-product(shape(nzatom))*kind(nzatom)
  deallocate(nzatom,stat=i_stat)
  call memocc(i_stat,i_all,'nzatom','memguess')
  i_all=-product(shape(iasctype))*kind(iasctype)
  deallocate(iasctype,stat=i_stat)
  call memocc(i_stat,i_all,'iasctype','memguess')

  if (optimise) then
     call optimise_volume(nat,ntypes,iatype,atomnames,in%crmult,in%frmult,in%hgrid,rxyz,radii_cf)
  end if

! Determine size alat of overall simulation cell and shift atom positions
! then calculate the size in units of the grid space
  hx=In%hgrid
  hy=In%hgrid
  hz=In%hgrid
  call system_size(0,geocode,nat,ntypes,rxyz,radii_cf,in%crmult,in%frmult,hx,hy,hz,&
       iatype,atomnames,alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)

  call MemoryEstimator(geocode,nproc,in%idsx,n1,n2,n3,alat1,alat2,alat3,hx,hy,hz,nat,&
       ntypes,iatype,rxyz,radii_cf,in%crmult,in%frmult,norb,atomnames,output_grid,in%nspin,peakmem)

  i_all=-product(shape(atomnames))*kind(atomnames)
  deallocate(atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames','memguess')
  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf','memguess')
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz','memguess')
  i_all=-product(shape(iatype))*kind(iatype)
  deallocate(iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype','memguess')

  !finalize memory counting
  call memocc(0,0,'count','stop')
  
end program memguess
!!***

!rotate the molecule via an orthogonal matrix in order to minimise the
!volume of the cubic cell
subroutine optimise_volume(nat,ntypes,iatype,atomnames,crmult,frmult,hgrid,rxyz,radii_cf)
  implicit none
  integer, intent(in) :: nat,ntypes
  real(kind=8), intent(in) :: crmult,frmult,hgrid
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  integer, dimension(nat), intent(in) :: iatype
  real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
  real(kind=8), dimension(3,nat), intent(inout) :: rxyz
  !local variables
  integer :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,n1i,n2i,n3i,iat,i_all,i_stat,it,i
  real(kind=8) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax,alat1,alat2,alat3
  real(kind=8), dimension(3,3) :: urot
  real(kind=8), dimension(:,:), allocatable :: txyz

  allocate(txyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(txyz))*kind(txyz),'txyz','optimise_volume')

  call system_size(1,'F',nat,ntypes,rxyz,radii_cf,crmult,frmult,hgrid,hgrid,hgrid,iatype,&
       atomnames,alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)
  !call volume(nat,rxyz,vol)
  vol=alat1*alat2*alat3
  write(*,'(1x,a,1pe16.8)')'Initial volume (Bohr^3)',vol

  it=0
  diag=1.d-2 ! initial small diagonal element allows for search over all angles
  loop_rotations: do  ! loop over all trial rotations
     diag=diag*1.0001d0 ! increase diag to search over smaller angles
     it=it+1
     if (diag.gt.100.d0) exit loop_rotations ! smaller angle rotations do not make sense

     ! create a random orthogonal (rotation) matrix
     call random_number(urot)
     urot(:,:)=urot(:,:)-.5d0
     do i=1,3
        urot(i,i)=urot(i,i)+diag
     enddo

     s=urot(1,1)**2+urot(2,1)**2+urot(3,1)**2
     s=1.d0/sqrt(s)
     urot(:,1)=s*urot(:,1) 

     s=urot(1,1)*urot(1,2)+urot(2,1)*urot(2,2)+urot(3,1)*urot(3,2)
     urot(:,2)=urot(:,2)-s*urot(:,1)
     s=urot(1,2)**2+urot(2,2)**2+urot(3,2)**2
     s=1.d0/sqrt(s)
     urot(:,2)=s*urot(:,2) 

     s=urot(1,1)*urot(1,3)+urot(2,1)*urot(2,3)+urot(3,1)*urot(3,3)
     urot(:,3)=urot(:,3)-s*urot(:,1)
     s=urot(1,2)*urot(1,3)+urot(2,2)*urot(2,3)+urot(3,2)*urot(3,3)
     urot(:,3)=urot(:,3)-s*urot(:,2)
     s=urot(1,3)**2+urot(2,3)**2+urot(3,3)**2
     s=1.d0/sqrt(s)
     urot(:,3)=s*urot(:,3) 

     ! eliminate reflections
     if (urot(1,1).le.0.d0) urot(:,1)=-urot(:,1)
     if (urot(2,2).le.0.d0) urot(:,2)=-urot(:,2)
     if (urot(3,3).le.0.d0) urot(:,3)=-urot(:,3)

     ! apply the rotation to all atomic positions! 
     do iat=1,nat
        x=rxyz(1,iat) ; y=rxyz(2,iat) ; z=rxyz(3,iat)
        txyz(:,iat)=x*urot(:,1)+y*urot(:,2)+z*urot(:,3)
     enddo

     call system_size(1,'F',nat,ntypes,txyz,radii_cf,crmult,frmult,hgrid,hgrid,hgrid,iatype,&
          atomnames,alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)
     tvol=alat1*alat2*alat3
     !call volume(nat,txyz,tvol)
     if (tvol.lt.vol) then
        write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol,it,diag
        rxyz(:,:)=txyz(:,:)
        vol=tvol
        dmax=max(alat1,alat2,alat3)
        ! if box longest along x switch x and z
        if (alat1 == dmax)  then
           do  iat=1,nat
              tx=rxyz(1,iat)
              tz=rxyz(3,iat)

              rxyz(1,iat)=tz
              rxyz(3,iat)=tx
           enddo
           ! if box longest along y switch y and z
        else if (alat2 == dmax)  then
           do  iat=1,nat
              ty=rxyz(1,iat) ; tz=rxyz(3,iat)
              rxyz(1,iat)=tz ; rxyz(3,iat)=ty
           enddo
        endif
     endif
  end do loop_rotations

  i_all=-product(shape(txyz))*kind(txyz)
  deallocate(txyz,stat=i_stat)
  call memocc(i_stat,i_all,'txyz','optimise_volume')
end subroutine optimise_volume

