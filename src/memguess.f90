!!****p* BigDFT/memguess
!! NAME
!!   memguess
!!
!! FUNCTION
!!  Test the input files and estimates the memory occupation versus the number
!!  of processeors
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
  logical :: calc_tail,output_grid
  character(len=20) :: tatonam,units
  integer :: ierror,nat,ntypes,nproc,n1,n2,n3,i_stat,i_all
  integer :: nelec,natsc,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  integer :: norb,norbu,norbd,norbe,norbp,norbsc
  integer :: iunit,ityp
  real(kind=8) :: alat1,alat2,alat3,peakmem
  type(input_variables) :: in
  character(len=20), dimension(:), allocatable :: atomnames
  integer, dimension(:), allocatable :: iatype,nelpsp,nzatom,npspcode,iasctype
  integer, dimension(:,:), allocatable :: neleconf
  real(kind=8), dimension(:,:), allocatable :: rxyz,radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psppar, psiat
  real(kind=8), dimension(:,:), allocatable :: xp, occupat
  integer, dimension(:,:), allocatable :: nl
  logical, dimension(:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng,norbsc_arr
  real(kind=8), dimension(:), allocatable :: occup,spinar

! Get arguments
  call getarg(1,tatonam)

  if(trim(tatonam)=='') then
     write(*,'(1x,a)')&
          'Usage: ./memguess <nproc> [y]'
     write(*,'(1x,a)')&
          'Indicate the number of processes after the executable'
     write(*,'(1x,a)')&
          'You can put a "y" in the second argument (optional) if you want the '
     write(*,'(1x,a)')&
          '  grid to be plotted with V_Sim'
     stop
  else
     read(unit=tatonam,fmt=*) nproc
     call getarg(2,tatonam)
     if(trim(tatonam)=='') then
        output_grid=.false.
     else if (trim(tatonam)=='y') then
        output_grid=.true.
        write(*,'(1x,a)')&
             'The system grid will be displayed in the "grid.ascii" file'
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
  write(*,'(1x,a,i0)') 'Number of atoms     = ',nat

  allocate(rxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz))*kind(rxyz),'rxyz','memguess')
  allocate(iatype(nat),stat=i_stat)
  call memocc(i_stat,product(shape(iatype))*kind(iatype),'iatype','memguess')
  allocate(atomnames(100),stat=i_stat) 
  call memocc(i_stat,product(shape(atomnames))*kind(atomnames),'atomnames','memguess')

  !read atomic positions
  call read_atomic_positions(0,99,units,nat,ntypes,iatype,atomnames,rxyz)

  close(99)

  write(*,'(1x,a,i0)') 'Number of atom types= ',ntypes

  do ityp=1,ntypes
     write(*,'(1x,a,i0,a,a)') 'Atoms of type ',ityp,' are ',trim(atomnames(ityp))
  enddo

  !new way of reading the input variables, use structures
  call read_input_variables(0,in)

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

  call read_system_variables(0,nproc,nat,ntypes,in%nspin,in%ncharge,in%mpol,atomnames,iatype,&
       psppar,radii_cf,npspcode,iasctype,nelpsp,nzatom,nelec,natsc,norb,norbu,norbd,norbp,iunit)

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


! Determine size alat of overall simulation cell and shift atom positions
! then calculate the size in units of the grid space
  call system_size(0,nat,ntypes,rxyz,radii_cf,in%crmult,in%frmult,in%hgrid,iatype,atomnames, &
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)


  call MemoryEstimator(nproc,in%idsx,n1,n2,n3,alat1,alat2,alat3,in%hgrid,nat,ntypes,iatype,&
          rxyz,radii_cf,in%crmult,in%frmult,norb,atomnames,output_grid,in%nspin,peakmem)

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
