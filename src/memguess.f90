program memguess

  use libBigDFT

  implicit none
  logical :: calc_tail,output_grid
  character(len=20) :: tatonam,units
  integer :: ierror,nat,ntypes,nproc,n1,n2,n3,i_stat,i_all
  integer :: nelec,norb,natsc,norbp,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  integer :: norbu,norbd,iunit,ityp
  real(kind=8) :: alat1,alat2,alat3
  type(input_variables) :: in
  character(len=20), dimension(:), allocatable :: atomnames
  integer, dimension(:), allocatable :: iatype,nelpsp,nzatom,npspcode,iasctype
  integer, dimension(:,:), allocatable :: neleconf
  real(kind=8), dimension(:,:), allocatable :: rxyz,radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psppar
  real(kind=8), dimension(:), allocatable :: occup,spinar

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

  ! read atomic positions
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

  allocate(occup(norb),stat=i_stat)
  call memocc(i_stat,product(shape(occup))*kind(occup),'occup','memguess')
  allocate(spinar(norb),stat=i_stat)
  call memocc(i_stat,product(shape(spinar))*kind(spinar),'occup','memguess')
  
! Occupation numbers
  call input_occup(0,iunit,nelec,norb,norbu,norbd,in%nspin,occup,spinar)

! Determine size alat of overall simulation cell and shift atom positions
! then calculate the size in units of the grid space
  call system_size(0,nat,ntypes,rxyz,radii_cf,in%crmult,in%frmult,in%hgrid,iatype,atomnames, &
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)


  i_all=-product(shape(psppar))*kind(psppar)
  deallocate(psppar,stat=i_stat)
  call memocc(i_stat,i_all,'psppar','memguess')
  i_all=-product(shape(nelpsp))*kind(nelpsp)
  deallocate(nelpsp,stat=i_stat)
  call memocc(i_stat,i_all,'nelpsp','memguess')
  i_all=-product(shape(occup))*kind(occup)
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup','memguess')
  i_all=-product(shape(spinar))*kind(spinar)
  deallocate(spinar,stat=i_stat)
  call memocc(i_stat,i_all,'spinar','memguess')
  !no need of using nzatom array
  i_all=-product(shape(nzatom))*kind(nzatom)
  deallocate(nzatom,stat=i_stat)
  call memocc(i_stat,i_all,'nzatom','memguess')
  i_all=-product(shape(iasctype))*kind(iasctype)
  deallocate(iasctype,stat=i_stat)
  call memocc(i_stat,i_all,'iasctype','memguess')

  call MemoryEstimator(nproc,in%idsx,n1,n2,n3,alat1,alat2,alat3,in%hgrid,nat,ntypes,iatype,&
          rxyz,radii_cf,in%crmult,in%frmult,norb,atomnames,output_grid,in%nspin)

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
