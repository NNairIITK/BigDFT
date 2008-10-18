!!****p* BigDFT/memguess
!! FUNCTION
!!  Test the input files and estimates the memory occupation versus the number
!!  of processors
!! AUTHOR
!!    Luigi Genovese
!! COPYRIGHT
!!    Copyright (C) 2007 CEA
!!
!! SOURCE
!!
program memguess

  use module_base
  use module_types
  use module_interfaces

  implicit none
  character(len=*), parameter :: subname='memguess'
  integer, parameter :: ngx=31
  character(len=20) :: tatonam,units
  logical :: calc_tail,output_grid,optimise
  integer :: ierror,nproc,n1,n2,n3,i_stat,i_all
  integer :: nelec,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  integer :: norb,norbu,norbd,norbe,norbp,norbsc
  integer :: iunit,ityp
  real(kind=8) :: peakmem,hx,hy,hz
  type(input_variables) :: in
  type(atoms_data) :: atoms
  type(nonlocal_psp_descriptors) :: nlpspd
  logical, dimension(:,:,:), allocatable :: logrid
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(kind=8), dimension(:,:), allocatable :: rxyz,radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psiat
  real(kind=8), dimension(:,:), allocatable :: xp, occupat
  integer, dimension(:,:), allocatable :: nl
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng
  real(kind=8), dimension(:), allocatable :: occup,spinsgn

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
  read(99,*) atoms%nat,units
 
  allocate(rxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',subname)

  !read atomic positions
  call read_atomic_positions(0,99,units,in,atoms,rxyz)

  close(99)

  !new way of reading the input variables, use structures
  call read_input_variables(0,'input.dat',in)

  call print_input_parameters(in)

  write(*,'(1x,a)')&
       '------------------------------------------------------------------ System Properties'
 
  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)
 
  call read_system_variables(0,nproc,in,atoms,radii_cf,nelec,norb,norbu,norbd,norbp,iunit)

! Allocations for the occupation numbers
  allocate(occup(norb+ndebug),stat=i_stat)
  call memocc(i_stat,occup,'occup',subname)
  allocate(spinsgn(norb+ndebug),stat=i_stat)
  call memocc(i_stat,spinsgn,'spinsgn',subname)
  
! Occupation numbers
  call input_occup(0,iunit,nelec,norb,norbu,norbd,in%nspin,in%mpol,occup,spinsgn)

! De-allocations
  i_all=-product(shape(occup))*kind(occup)
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup',subname)
  i_all=-product(shape(spinsgn))*kind(spinsgn)
  deallocate(spinsgn,stat=i_stat)
  call memocc(i_stat,i_all,'spinsgn',subname)
  
  !in the case in which the number of orbitals is not "trivial" check whether they are too many
  if ( max(norbu,norbd) /= ceiling(real(nelec,kind=4)/2.0) ) then
     ! Allocations for readAtomicOrbitals (check inguess.dat and psppar files + give norbe)
     allocate(xp(ngx,atoms%ntypes+ndebug),stat=i_stat)
     call memocc(i_stat,xp,'xp',subname)
     allocate(psiat(ngx,5,atoms%ntypes+ndebug),stat=i_stat)
     call memocc(i_stat,psiat,'psiat',subname)
     allocate(occupat(5,atoms%ntypes+ndebug),stat=i_stat)
     call memocc(i_stat,occupat,'occupat',subname)
     allocate(ng(atoms%ntypes+ndebug),stat=i_stat)
     call memocc(i_stat,ng,'ng',subname)
     allocate(nl(4,atoms%ntypes+ndebug),stat=i_stat)
     call memocc(i_stat,nl,'nl',subname)
     allocate(scorb(4,2,atoms%natsc+ndebug),stat=i_stat)
     call memocc(i_stat,scorb,'scorb',subname)
     allocate(norbsc_arr(atoms%natsc+1,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)

     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(0,ngx,xp,psiat,occupat,ng,nl,atoms,norbe,norbsc,in%nspin,&
          scorb,norbsc_arr)

     ! De-allocations
     i_all=-product(shape(xp))*kind(xp)
     deallocate(xp,stat=i_stat)
     call memocc(i_stat,i_all,'xp',subname)
     i_all=-product(shape(psiat))*kind(psiat)
     deallocate(psiat,stat=i_stat)
     call memocc(i_stat,i_all,'psiat',subname)
     i_all=-product(shape(occupat))*kind(occupat)
     deallocate(occupat,stat=i_stat)
     call memocc(i_stat,i_all,'occupat',subname)
     i_all=-product(shape(ng))*kind(ng)
     deallocate(ng,stat=i_stat)
     call memocc(i_stat,i_all,'ng',subname)
     i_all=-product(shape(nl))*kind(nl)
     deallocate(nl,stat=i_stat)
     call memocc(i_stat,i_all,'nl',subname)
     i_all=-product(shape(scorb))*kind(scorb)
     deallocate(scorb,stat=i_stat)
     call memocc(i_stat,i_all,'scorb',subname)
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr',subname)

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


  if (optimise) then
     call optimise_volume(atoms,in%crmult,in%frmult,in%hgrid,rxyz,radii_cf)
     write(*,'(1x,a)')'Writing optimised positions in file posout_000.xyz...'
     call wtposout(0,0.d0,rxyz,atoms)
  end if

! Determine size alat of overall simulation cell and shift atom positions
! then calculate the size in units of the grid space
  hx=in%hgrid
  hy=in%hgrid
  hz=in%hgrid

  call system_size(0,in%geocode,atoms,rxyz,radii_cf,in%crmult,in%frmult,hx,hy,hz,&
       in%alat1,in%alat2,in%alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)

  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  allocate(nlpspd%nseg_p(0:2*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nseg_p,'nseg_p',subname)
  allocate(nlpspd%nvctr_p(0:2*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nvctr_p,'nvctr_p',subname)
  allocate(nlpspd%nboxp_c(2,3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_c,'nboxp_c',subname)
  allocate(nlpspd%nboxp_f(2,3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_f,'nboxp_f',subname)
  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors(in%geocode,0,n1,n2,n3,hx,hy,hz,in%cpmult,in%fpmult,rxyz,radii_cf,&
       logrid,atoms,nlpspd)

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)
  i_all=-product(shape(nlpspd%nvctr_p))*kind(nlpspd%nvctr_p)
  deallocate(nlpspd%nvctr_p,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_p',subname)
  i_all=-product(shape(nlpspd%nseg_p))*kind(nlpspd%nseg_p)
  deallocate(nlpspd%nseg_p,stat=i_stat)
  call memocc(i_stat,i_all,'nseg_p',subname)
  i_all=-product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c)
  deallocate(nlpspd%nboxp_c,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_c',subname)
  i_all=-product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f)
  deallocate(nlpspd%nboxp_f,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_f',subname)
  i_all=-product(shape(atoms%lfrztyp))*kind(atoms%lfrztyp)
  deallocate(atoms%lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp',subname)
  i_all=-product(shape(atoms%natpol))*kind(atoms%natpol)
  deallocate(atoms%natpol,stat=i_stat)
  call memocc(i_stat,i_all,'natpol',subname)
  i_all=-product(shape(atoms%psppar))*kind(atoms%psppar)
  deallocate(atoms%psppar,stat=i_stat)
  call memocc(i_stat,i_all,'psppar',subname)
  i_all=-product(shape(atoms%npspcode))*kind(atoms%npspcode)
  deallocate(atoms%npspcode,stat=i_stat)
  call memocc(i_stat,i_all,'npspcode',subname)
  i_all=-product(shape(atoms%nelpsp))*kind(atoms%nelpsp)
  deallocate(atoms%nelpsp,stat=i_stat)
  call memocc(i_stat,i_all,'nelpsp',subname)
  !no need of using nzatom array
  i_all=-product(shape(atoms%nzatom))*kind(atoms%nzatom)
  deallocate(atoms%nzatom,stat=i_stat)
  call memocc(i_stat,i_all,'nzatom',subname)
  i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
  deallocate(atoms%iasctype,stat=i_stat)
  call memocc(i_stat,i_all,'iasctype',subname)


  call MemoryEstimator(in%geocode,nproc,in%idsx,n1,n2,n3,in%alat1,in%alat2,in%alat3,&
       hx,hy,hz,atoms%nat,atoms%ntypes,atoms%iatype,rxyz,radii_cf,in%crmult,in%frmult,&
       norb,nlpspd%nprojel,atoms%atomnames,output_grid,in%nspin,peakmem)

  i_all=-product(shape(atoms%atomnames))*kind(atoms%atomnames)
  deallocate(atoms%atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames',subname)
  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  i_all=-product(shape(atoms%iatype))*kind(atoms%iatype)
  deallocate(atoms%iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype',subname)

  !finalize memory counting
  call memocc(0,0,'count','stop')
  
end program memguess
!!***


!!****f* BigDFT/optimise_volume
!! FUNCTION
!!  Rotate the molecule via an orthogonal matrix in order to minimise the
!!  volume of the cubic cell
!! AUTHOR
!!    Luigi Genovese
!! SOURCE
!!
subroutine optimise_volume(atoms,crmult,frmult,hgrid,rxyz,radii_cf)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(kind=8), intent(in) :: crmult,frmult,hgrid
  real(kind=8), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
  !local variables
  character(len=*), parameter :: subname='optimise_volume'
  integer :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,n1i,n2i,n3i,iat,i_all,i_stat,it,i
  real(kind=8) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax,alat1,alat2,alat3
  real(kind=8), dimension(3,3) :: urot
  real(kind=8), dimension(:,:), allocatable :: txyz

  allocate(txyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)

  call system_size(1,'F',atoms,rxyz,radii_cf,crmult,frmult,hgrid,hgrid,hgrid,&
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)
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
     do iat=1,atoms%nat
        x=rxyz(1,iat) ; y=rxyz(2,iat) ; z=rxyz(3,iat)
        txyz(:,iat)=x*urot(:,1)+y*urot(:,2)+z*urot(:,3)
     enddo

     call system_size(1,'F',atoms,txyz,radii_cf,crmult,frmult,hgrid,hgrid,hgrid,&
          alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)
     tvol=alat1*alat2*alat3
     !call volume(nat,txyz,tvol)
     if (tvol.lt.vol) then
        write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol,it,diag
        rxyz(:,:)=txyz(:,:)
        vol=tvol
        dmax=max(alat1,alat2,alat3)
        ! if box longest along x switch x and z
        if (alat1 == dmax)  then
           do  iat=1,atoms%nat
              tx=rxyz(1,iat)
              tz=rxyz(3,iat)

              rxyz(1,iat)=tz
              rxyz(3,iat)=tx
           enddo
           ! if box longest along y switch y and z
        else if (alat2 == dmax .and. alat1 /= dmax)  then
           do  iat=1,atoms%nat
              ty=rxyz(2,iat) 
              tz=rxyz(3,iat)

              rxyz(2,iat)=tz 
              rxyz(3,iat)=ty
           enddo
        endif
     endif
  end do loop_rotations

  i_all=-product(shape(txyz))*kind(txyz)
  deallocate(txyz,stat=i_stat)
  call memocc(i_stat,i_all,'txyz',subname)
end subroutine optimise_volume
!!***
