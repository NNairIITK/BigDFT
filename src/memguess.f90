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
  character(len=20) :: tatonam
  character(len=40) :: comment
  logical :: optimise,GPUtest,convert=.false.,exists
  integer :: nelec,ntimes,nproc,i_stat,ierror,i_all,output_grid
  integer :: norbe,norbsc,nspin,iorb,norbu,norbd,nspinor,norb
  integer :: norbgpu,nspin_ig
  real(kind=8) :: peakmem,hx,hy,hz
  type(input_variables) :: in
  type(atoms_data) :: atoms
  type(orbitals_data) :: orbs,orbstst
  type(locreg_descriptors) :: Glr
  type(nonlocal_psp_descriptors) :: nlpspd
  logical, dimension(:,:,:), allocatable :: logrid
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:,:), pointer :: rxyz
  real(kind=8), dimension(:,:), allocatable :: radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psiat
  real(kind=8), dimension(:,:), allocatable :: xp, occupat
  integer, dimension(:,:), allocatable :: nl
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng
  real(kind=8), dimension(:), allocatable :: locrad

! Get arguments
  call getarg(1,tatonam)
  
  optimise=.false.
  GPUtest=.false.
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
     write(*,'(1x,a)')&
          'In the case of a CUDAGPU calculation you can put "GPUtest" after the'
     write(*,'(1x,a)')&
          '  number of processors, followed by the number of repeats'
     write(*,'(1x,a)')&
          'You can also put "convert" keyword for converting input file into 1.3 format'
     stop
  else
     read(unit=tatonam,fmt=*) nproc
     call getarg(2,tatonam)
     if(trim(tatonam)=='') then
        output_grid=0
     else if (trim(tatonam)=='y') then
        output_grid=1
        write(*,'(1x,a)')&
             'The system grid will be displayed in the "grid.xyz" file'
     else if (trim(tatonam)=='o') then
        optimise=.true.
        output_grid=1
        write(*,'(1x,a)')&
             'The optimised system grid will be displayed in the "grid.xyz" file'
     else if (trim(tatonam)=='GPUtest') then
        GPUtest=.true.
        write(*,'(1x,a)')&
             'Perform the test with GPU, if present.'
        call getarg(3,tatonam)
        ntimes=1
        norbgpu=0
        read(tatonam,*,iostat=ierror)ntimes
        if (ierror==0) then
           write(*,'(1x,a,i0,a)')&
                'Repeat each calculation ',ntimes,' times.'
           call getarg(4,tatonam)
           read(tatonam,*,iostat=ierror)norbgpu
        end if
     else if (trim(tatonam)=='convert') then
        convert=.true.
        write(*,'(1x,a)')&
             'convert the input.dat file in the "input_convert.dft" (1.3 format)'
     else
        write(*,'(1x,a)')&
             'Usage: ./memguess <nproc> [y]'
        write(*,'(1x,a)')&
             'Indicate the number of processes after the executable'
        write(*,'(1x,a)')&
             'ERROR: The only second argument which is accepted are "y", "o", "convert" or "GPUtest"'
        stop
     end if
  end if

  !initialize memory counting
  call memocc(0,0,'count','start')

  !welcome screen
  call print_logo()

  !read number of atoms
  call read_atomic_file('posinp',0,atoms,rxyz)

  if (convert) then
     call read_input_variables(0,'input.dft',in)
     write(*,'(a)',advance='NO')' Conversion of the input file...'
     call dft_input_converter(in)
     write(*,*)' ...done'
  else
     call dft_input_variables(0,'input.dft',in)
     !read geometry optimsation input variables
     !inquire for the file needed for geometry optimisation
     !if not present, perform a simple geometry optimisation
     inquire(file="input.geopt",exist=exists)
     if (exists) then
        call geopt_input_variables(0,'input.geopt',in)
     else
        call geopt_input_variables_default(in)
     end if
  end if


  call print_input_parameters(in,atoms)


  write(*,'(1x,a)')&
       '------------------------------------------------------------------ System Properties'
 

  ! store PSP parameters
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call system_properties(0,nproc,in,atoms,orbs,radii_cf,nelec)

  i_all=-product(shape(orbs%norb_par))*kind(orbs%norb_par)
  deallocate(orbs%norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'orbs%norb_par',subname)

  if (optimise) then
     if (atoms%geocode =='F') then
        call optimise_volume(atoms,in%crmult,in%frmult,in%hx,in%hy,in%hz,rxyz,radii_cf)
     else
        call shift_periodic_directions(atoms,rxyz,radii_cf)
     end if
     write(*,'(1x,a)')'Writing optimised positions in file posout_000.[xyz,ascii]...'
     write(comment,'(a)')'POSITIONS IN OPTIMIZED CELL '
     call write_atomic_file('posopt',0.d0,rxyz,atoms,trim(comment))
     !call wtxyz('posopt',0.d0,rxyz,atoms,trim(comment))

  end if

  
  !in the case in which the number of orbitals is not "trivial" check whether they are too many
  if ( max(orbs%norbu,orbs%norbd) /= ceiling(real(nelec,kind=4)/2.0)) then
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
     allocate(locrad(atoms%nat+ndebug),stat=i_stat)
     call memocc(i_stat,locrad,'locrad',subname)

     if (in%nspin==4) then
        nspin_ig=1
     else
        nspin_ig=in%nspin
     end if


     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(0,ngx,xp,psiat,occupat,ng,nl,atoms,norbe,norbsc,nspin_ig,&
          scorb,norbsc_arr,locrad)

     if (in%nspin==4) then
        !in that case the number of orbitals doubles
        norbe=2*norbe
     end if

     ! De-allocations
     i_all=-product(shape(locrad))*kind(locrad)
     deallocate(locrad,stat=i_stat)
     call memocc(i_stat,i_all,'locrad',subname)
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
     if (in%nspin==1 .or. in%nspin==4) then
        if (orbs%norb>norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals (',orbs%norb,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
     else if (in%nspin == 2) then
        if (orbs%norbu > norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals up (',orbs%norbu,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
        if (orbs%norbd > norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals down (',orbs%norbd,&
                ') must not be greater than the number of orbitals (',norbe,&
                ') generated from the input guess.'
           stop
        end if
     end if

  end if

! Determine size alat of overall simulation cell and shift atom positions
! then calculate the size in units of the grid space
  hx=in%hx
  hy=in%hy
  hz=in%hz

  call system_size(0,atoms,rxyz,radii_cf,in%crmult,in%frmult,hx,hy,hz,Glr)

  ! De-allocations
  i_all=-product(shape(orbs%occup))*kind(orbs%occup)
  deallocate(orbs%occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup',subname)
  i_all=-product(shape(orbs%spinsgn))*kind(orbs%spinsgn)
  deallocate(orbs%spinsgn,stat=i_stat)
  call memocc(i_stat,i_all,'spinsgn',subname)
  i_all=-product(shape(orbs%kpts))*kind(orbs%kpts)
  deallocate(orbs%kpts,stat=i_stat)
  call memocc(i_stat,i_all,'orbs%kpts',subname)
  i_all=-product(shape(orbs%kwgts))*kind(orbs%kwgts)
  deallocate(orbs%kwgts,stat=i_stat)
  call memocc(i_stat,i_all,'orbs%kwgts',subname)
  i_all=-product(shape(orbs%iokpt))*kind(orbs%iokpt)
  deallocate(orbs%iokpt,stat=i_stat)
  call memocc(i_stat,i_all,'orbs%iokpt',subname)



  if (GPUtest .and. .not. GPUconv) then
     write(*,*)' ERROR: you can not put a GPUtest flag is there is no GPUrun.'
     stop
  end if
  if (GPUconv .and. atoms%geocode=='P' .and. GPUtest) then
     !test the hamiltonian in CPU or GPU
     !create the orbitals data structure for one orbital
     allocate(orbstst%norb_par(0:0+ndebug),stat=i_stat)
     call memocc(i_stat,orbstst%norb_par,'orbstst%norb_par',subname)
     !test orbitals
     nspin=1
     if (norbgpu == 0) then
        norb=orbs%norb
     else
        norb=norbgpu
     end if
     norbu=norb
     norbd=0
     nspinor=1

     call orbitals_descriptors(0,nproc,norb,norbu,norbd,nspinor,orbstst)
     allocate(orbstst%occup(orbstst%norb+ndebug),stat=i_stat)
     call memocc(i_stat,orbstst%occup,'orbstst%occup',subname)
     allocate(orbstst%spinsgn(orbstst%norb+ndebug),stat=i_stat)
     call memocc(i_stat,orbstst%spinsgn,'orbstst%spinsgn',subname)
     allocate(orbstst%eval(orbstst%norbp+ndebug),stat=i_stat)
     call memocc(i_stat,orbstst%eval,'orbstst%eval',subname)
     do iorb=1,orbstst%norbp
        orbstst%eval(iorb)=-0.5_gp
     end do

     do iorb=1,orbstst%norb
        orbstst%occup(iorb)=1.0_gp
        orbstst%spinsgn(iorb)=1.0_gp
     end do

     call createWavefunctionsDescriptors(0,hx,hy,hz,&
          atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr,orbstst)
     
     call compare_cpu_gpu_hamiltonian(0,1,atoms,orbstst,nspin,in%ncong,in%ixc,&
          Glr,hx,hy,hz,rxyz,ntimes)
     
     call deallocate_wfd(Glr%wfd,subname)

     i_all=-product(shape(orbstst%norb_par))*kind(orbstst%norb_par)
     deallocate(orbstst%norb_par,stat=i_stat)
     call memocc(i_stat,i_all,'orbstst%norb_par',subname)
     i_all=-product(shape(orbstst%spinsgn))*kind(orbstst%spinsgn)
     deallocate(orbstst%spinsgn,stat=i_stat)
     call memocc(i_stat,i_all,'orbstst%spinsgn',subname)
     i_all=-product(shape(orbstst%occup))*kind(orbstst%occup)
     deallocate(orbstst%occup,stat=i_stat)
     call memocc(i_stat,i_all,'orbstst%occup',subname)

     i_all=-product(shape(orbstst%eval))*kind(orbstst%eval)
     deallocate(orbstst%eval,stat=i_stat)
     call memocc(i_stat,i_all,'orbstst%eval',subname)
     i_all=-product(shape(orbstst%occup))*kind(orbstst%occup)

  end if


  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  allocate(nlpspd%nseg_p(0:2*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nseg_p,'nseg_p',subname)
  allocate(nlpspd%nvctr_p(0:2*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nvctr_p,'nvctr_p',subname)
  allocate(nlpspd%nboxp_c(2,3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_c,'nboxp_c',subname)
  allocate(nlpspd%nboxp_f(2,3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_f,'nboxp_f',subname)
  allocate(logrid(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors(0,Glr%d%n1,Glr%d%n2,Glr%d%n3,hx,hy,hz,&
       in%frmult,in%frmult,rxyz,radii_cf,logrid,atoms,nlpspd)

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
  i_all=-product(shape(atoms%ifrztyp))*kind(atoms%ifrztyp)
  deallocate(atoms%ifrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'ifrztyp',subname)
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


  call MemoryEstimator(atoms%geocode,nproc,in%idsx,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
       atoms%alat1,atoms%alat2,atoms%alat3,&
       hx,hy,hz,atoms%nat,atoms%ntypes,atoms%iatype,rxyz,radii_cf,in%crmult,in%frmult,&
       orbs%norb,nlpspd%nprojel,atoms%atomnames,output_grid,in%nspin,peakmem)

  !add the comparison between cuda hamiltonian and normal one if it is the case


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
!!    Stefan Goedecker, Luigi Genovese
!! SOURCE
!!
subroutine optimise_volume(atoms,crmult,frmult,hx,hy,hz,rxyz,radii_cf)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), intent(in) :: crmult,frmult,hx,hy,hz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  !local variables
  character(len=*), parameter :: subname='optimise_volume'
  integer :: iat,i_all,i_stat,it,i
  real(gp) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax
  type(locreg_descriptors) :: Glr
  real(gp), dimension(3,3) :: urot
  real(gp), dimension(:,:), allocatable :: txyz

  allocate(txyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)
  call system_size(1,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr)
  !call volume(nat,rxyz,vol)
  vol=atoms%alat1*atoms%alat2*atoms%alat3
  write(*,'(1x,a,1pe16.8)')'Initial volume (Bohr^3)',vol

  it=0
  diag=1.d-2 ! initial small diagonal element allows for search over all angles
  loop_rotations: do  ! loop over all trial rotations
     diag=diag*1.0001_gp ! increase diag to search over smaller angles
     it=it+1
     if (diag > 100._gp) exit loop_rotations ! smaller angle rotations do not make sense

     ! create a random orthogonal (rotation) matrix
     call random_number(urot)
     urot(:,:)=urot(:,:)-.5_gp
     do i=1,3
        urot(i,i)=urot(i,i)+diag
     enddo

     s=urot(1,1)**2+urot(2,1)**2+urot(3,1)**2
     s=1._gp/sqrt(s)
     urot(:,1)=s*urot(:,1) 

     s=urot(1,1)*urot(1,2)+urot(2,1)*urot(2,2)+urot(3,1)*urot(3,2)
     urot(:,2)=urot(:,2)-s*urot(:,1)
     s=urot(1,2)**2+urot(2,2)**2+urot(3,2)**2
     s=1._gp/sqrt(s)
     urot(:,2)=s*urot(:,2) 

     s=urot(1,1)*urot(1,3)+urot(2,1)*urot(2,3)+urot(3,1)*urot(3,3)
     urot(:,3)=urot(:,3)-s*urot(:,1)
     s=urot(1,2)*urot(1,3)+urot(2,2)*urot(2,3)+urot(3,2)*urot(3,3)
     urot(:,3)=urot(:,3)-s*urot(:,2)
     s=urot(1,3)**2+urot(2,3)**2+urot(3,3)**2
     s=1._gp/sqrt(s)
     urot(:,3)=s*urot(:,3) 

     ! eliminate reflections
     if (urot(1,1) <= 0._gp) urot(:,1)=-urot(:,1)
     if (urot(2,2) <= 0._gp) urot(:,2)=-urot(:,2)
     if (urot(3,3) <= 0._gp) urot(:,3)=-urot(:,3)

     ! apply the rotation to all atomic positions! 
     do iat=1,atoms%nat
        x=rxyz(1,iat) 
        y=rxyz(2,iat) 
        z=rxyz(3,iat)

        txyz(:,iat)=x*urot(:,1)+y*urot(:,2)+z*urot(:,3)
     enddo

     call system_size(1,atoms,txyz,radii_cf,crmult,frmult,hx,hy,hz,Glr)
     tvol=atoms%alat1*atoms%alat2*atoms%alat3
     !call volume(nat,txyz,tvol)
     if (tvol < vol) then
        write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol,it,diag
        rxyz(:,:)=txyz(:,:)
        vol=tvol
        dmax=max(atoms%alat1,atoms%alat2,atoms%alat3)
        ! if box longest along x switch x and z
        if (atoms%alat1 == dmax)  then
           do  iat=1,atoms%nat
              tx=rxyz(1,iat)
              tz=rxyz(3,iat)

              rxyz(1,iat)=tz
              rxyz(3,iat)=tx
           enddo
           ! if box longest along y switch y and z
        else if (atoms%alat2 == dmax .and. atoms%alat1 /= dmax)  then
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

!!****f* BigDFT/shift_periodic_directions
!! FUNCTION
!!  Add a shift in the periodic directions such that the system
!!  uses as less as possible the modulo operation
!! AUTHOR
!!    Luigi Genovese
!! SOURCE
!!
subroutine shift_periodic_directions(at,rxyz,radii_cf)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  !local variables
  character(len=*), parameter :: subname='shift_periodic_directions'
  integer :: iat,i_all,i_stat,i,ityp
  real(gp) :: vol,tvol,maxsh,shiftx,shifty,shiftz
  real(gp), dimension(:,:), allocatable :: txyz

  !calculate maximum shift between these values
  !this is taken as five times the coarse radius around atoms
  maxsh=0.0_gp
  do ityp=1,at%ntypes
     maxsh=max(maxsh,5_gp*radii_cf(ityp,1))
  end do
  
  allocate(txyz(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)

  call calc_vol(at%geocode,at%nat,rxyz,vol)

  if (at%geocode /= 'F') then
     loop_shiftx: do i=1,5000 ! loop over all trial rotations
        ! create a random orthogonal (rotation) matrix
        call random_number(shiftx)

        !apply the shift to all atomic positions taking into account the modulo operation
        do iat=1,at%nat
           txyz(1,iat)=modulo(rxyz(1,iat)+shiftx*maxsh,at%alat1)
           txyz(2,iat)=rxyz(2,iat)
           txyz(3,iat)=rxyz(3,iat)
        end do

        call calc_vol(at%geocode,at%nat,txyz,tvol)
        !print *,'vol',tvol

        if (tvol < vol) then
           write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
           rxyz(:,:)=txyz(:,:)
           vol=tvol
        endif
     end do loop_shiftx
  end if

  if (at%geocode == 'P') then
     loop_shifty: do i=1,5000 ! loop over all trial rotations
        ! create a random orthogonal (rotation) matrix
        call random_number(shifty)

        !apply the shift to all atomic positions taking into account the modulo operation
        do iat=1,at%nat
           txyz(1,iat)=rxyz(1,iat)
           txyz(2,iat)=modulo(rxyz(2,iat)+shifty*maxsh,at%alat2)
           txyz(3,iat)=rxyz(3,iat)
        end do

        call calc_vol(at%geocode,at%nat,txyz,tvol)

        if (tvol < vol) then
           write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
           rxyz(:,:)=txyz(:,:)
           vol=tvol
        endif
     end do loop_shifty
  end if

    if (at%geocode /= 'F') then
     loop_shiftz: do i=1,5000 ! loop over all trial rotations
        ! create a random orthogonal (rotation) matrix
        call random_number(shiftz)

        !apply the shift to all atomic positions taking into account the modulo operation
        do iat=1,at%nat
           txyz(1,iat)=rxyz(1,iat)
           txyz(2,iat)=rxyz(2,iat)
           txyz(3,iat)=modulo(rxyz(3,iat)+shiftz*maxsh,at%alat3)
        end do

        call calc_vol(at%geocode,at%nat,txyz,tvol)

        if (tvol < vol) then
           write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol
           rxyz(:,:)=txyz(:,:)
           vol=tvol
        endif
     end do loop_shiftz
  end if

  


  i_all=-product(shape(txyz))*kind(txyz)
  deallocate(txyz,stat=i_stat)
  call memocc(i_stat,i_all,'txyz',subname)

end subroutine shift_periodic_directions
!!***

subroutine calc_vol(geocode,nat,rxyz,vol)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: nat
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(gp), intent(out) :: vol
  !local variables
  integer :: iat
  real(gp) :: cxmin,cxmax,cymin,cymax,czmin,czmax

  !calculate the extremes of the boxes taking into account the spheres around the atoms
  cxmax=-1.e10_gp 
  cxmin=1.e10_gp

  cymax=-1.e10_gp 
  cymin=1.e10_gp

  czmax=-1.e10_gp 
  czmin=1.e10_gp

  do iat=1,nat
     cxmax=max(cxmax,rxyz(1,iat)) 
     cxmin=min(cxmin,rxyz(1,iat))

     cymax=max(cymax,rxyz(2,iat)) 
     cymin=min(cymin,rxyz(2,iat))
     
     czmax=max(czmax,rxyz(3,iat)) 
     czmin=min(czmin,rxyz(3,iat))
  enddo
  !print *,cxmax,cxmin,cymax,cymin,czmax,czmin
  !now calculate the volume for the periodic part
  if (geocode == 'P') then
     vol=(cxmax-cxmin)*(cymax-cymin)*(czmax-czmin)
  else if (geocode == 'S') then
     vol=(cxmax-cxmin)*(czmax-czmin)
  end if

end subroutine calc_vol


subroutine compare_cpu_gpu_hamiltonian(iproc,nproc,at,orbs,nspin,ixc,ncong,&
     lr,hx,hy,hz,rxyz,ntimes)
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  implicit none
  integer, intent(in) :: iproc,nproc,nspin,ncong,ixc,ntimes
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables
  character(len=*), parameter :: subname='compare_cpu_gpu_hamiltonian'
  logical :: rsflag
  integer :: icoeff,i_stat,i_all,i1,i2,i3,ispin,j
  integer :: iorb,n3d,n3p,n3pi,i3xcsh,i3s,jproc,nrhotot,nspinn,nvctrp
  real(kind=4) :: tt,t0,t1
  real(gp) :: ttd,x,y,z,r2,arg,sigma2,ekin_sum,epot_sum,ekinGPU,epotGPU,gnrm,gnrmGPU
  real(gp) :: Rden,Rham,Rgemm,Rsyrk,Rprec
  real(kind=8) :: CPUtime,GPUtime
  type(gaussian_basis) :: G
  type(GPU_pointers) :: GPU
  integer, dimension(:,:), allocatable :: nscatterarr
  real(wp), dimension(:,:,:,:), allocatable :: pot,rho
  real(wp), dimension(:,:), allocatable :: gaucoeffs,psi,hpsi
  real(wp), dimension(:,:,:), allocatable :: overlap

  !nullify the G%rxyz pointer
  nullify(G%rxyz)
  !extract the gaussian basis from the pseudowavefunctions
  call gaussian_pswf_basis(iproc,at,rxyz,G)
  
  allocate(gaucoeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

  !fill randomly the gaussian coefficients for the orbitals considered
  do iorb=1,orbs%norbp*orbs%nspinor
     do icoeff=1,G%ncoeff
        call random_number(tt)
        gaucoeffs(icoeff,iorb)=real(tt,wp)
     end do
  end do

  !allocate the wavefunctions
  allocate(psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  call razero(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,psi)


  !convert the gaussians in wavelets
  call gaussians_to_wavelets(iproc,nproc,at%geocode,orbs,lr%d,&
       hx,hy,hz,lr%wfd,G,gaucoeffs,psi)

  i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
  deallocate(gaucoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'gaucoeffs',subname)

  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)

  !allocate and initialise the potential and the density
  allocate(pot(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,pot,'pot',subname)
  allocate(rho(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,rho,'rho',subname)

  !here the potential can be used for building the density
  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  !normally nproc=1
  do jproc=0,nproc-1
     call PS_dim4allocation(at%geocode,'D',jproc,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,&
          n3d,n3p,n3pi,i3xcsh,i3s)
     nscatterarr(jproc,1)=n3d
     nscatterarr(jproc,2)=n3p
     nscatterarr(jproc,3)=i3s+i3xcsh-1
     nscatterarr(jproc,4)=i3xcsh
  end do

  !components of the charge density
  if (orbs%nspinor ==4) then
     nspinn=4
  else
     nspinn=nspin
  end if

  !flag for toggling the REDUCE_SCATTER stategy
  rsflag=.not. (ixc >= 11 .and. ixc <=16)

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rsflag) then
     nrhotot=0
     do jproc=0,nproc-1
        nrhotot=nrhotot+nscatterarr(jproc,1)
     end do
  else
     nrhotot=lr%d%n3i
  end if


  !allocate the necessary objects on the GPU
  !set initialisation of GPU part 
  call prepare_gpu_for_locham(lr%d%n1,lr%d%n2,lr%d%n3,nspin,hx,hy,hz,lr%wfd,orbs,GPU)

  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Density calculation'

  call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,pot,nproc)

  !for each of the orbitals treated by the processor build the partial densities
  call cpu_time(t0)
  do j=1,ntimes
     call local_partial_density(iproc,nproc,rsflag,nscatterarr,&
          nrhotot,lr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,orbs,&
          psi,pot)
  end do
  call cpu_time(t1)
  CPUtime=real(t1-t0,kind=8)


  !copy the wavefunctions on GPU
  do iorb=1,orbs%norbp
     call GPU_send((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,&
          psi(1,(iorb-1)*orbs%nspinor+1),GPU%psi(iorb),i_stat)
  end do

  !now the GPU part
  !for each of the orbitals treated by the processor build the partial densities
  call cpu_time(t0)
  do j=1,ntimes
     call gpu_locden(lr,nspin,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,orbs,GPU)
  end do
  call cpu_time(t1)
  GPUtime=real(t1-t0,kind=8)

  !receive the density on GPU
  call GPU_receive(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin,rho,GPU%rhopot,i_stat)

  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)

  !compare the results between the different actions of the hamiltonian
  !check the differences between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,pot,rho,&
       lr%d%n1i*lr%d%n2i*lr%d%n3i,ntimes,.false.,Rden)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)


  !here the grid spacings are the small ones
  sigma2=0.125_gp*((lr%d%n1i*hx)**2+(lr%d%n2i*hy)**2+(lr%d%n3i*hz)**2)
  do ispin=1,nspin
     do i3=1,lr%d%n3i
        z=hz*real(i3-lr%d%n3i/2-1,gp)
        do i2=1,lr%d%n2i
           y=hy*real(i2-lr%d%n2i/2-1,gp)
           do i1=1,lr%d%n1i
              x=hx*real(i1-lr%d%n1i/2-1,gp)
              !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
              r2=x**2+y**2+z**2
              arg=0.5d0*r2/sigma2
              ttd=dexp(-arg)

              pot(i1,i2,i3,ispin)=ttd
           end do
        end do
     end do
  end do

  !copy the potential on GPU
  call GPU_send(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin,pot,GPU%rhopot,i_stat)


  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Local Hamiltonian calculation'

  !warm-up
  !call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum) 

  !apply the CPU hamiltonian
  !take timings
  call cpu_time(t0)
  do j=1,ntimes
     call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum) 
  end do
  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)

  print *,'ekin,epot=',ekin_sum,epot_sum

  !WARNING: local hamiltonian overwrites the psis
  !warm-up
  !call gpu_locham(lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,orbs,GPU,ekinGPU,epotGPU)

  !apply the GPU hamiltonian
  !take timings
  call cpu_time(t0)
  do j=1,ntimes
     call gpu_locham(lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,orbs,GPU,ekinGPU,epotGPU)
  end do
  call cpu_time(t1)

  print *,'ekinGPU,epotGPU',ekinGPU,epotGPU

  GPUtime=real(t1-t0,kind=8)
  
  !receive the data of GPU
  do iorb=1,orbs%norbp
     call GPU_receive((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,&
          psi(1,(iorb-1)*orbs%nspinor+1),GPU%psi(iorb),i_stat)
  end do
  

  !compare the results between the different actions of the hamiltonian
  !check the differences between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,hpsi,psi,&
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),ntimes,.false.,Rham)


  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)

  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Linear Algebra (Blas)'
 
  !perform the scalar product between the hpsi wavefunctions
  !actually this is <hpsi|hpsi> it has no meaning.
  !this works only if nspinor==1
  allocate(overlap(orbs%norbp,orbs%norbp,2+ndebug),stat=i_stat)
  call memocc(i_stat,overlap,'overlap',subname)

  call cpu_time(t0)
  do j=1,ntimes
     call DGEMM('T','N',orbs%norbp,orbs%norbp,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1.0_wp,&
          psi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
          hpsi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),0.0_wp,&
          overlap(1,1,1),orbs%norbp)
  end do
  call cpu_time(t1)
  CPUtime=real(t1-t0,kind=8)


  call cpu_time(t0)
  do j=1,ntimes
     call GEMM('T','N',orbs%norbp,orbs%norbp,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1.0_wp,&
          psi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
          hpsi(1,1),(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),0.0_wp,&
          overlap(1,1,2),orbs%norbp)
  end do
  call cpu_time(t1)
  GPUtime=real(t1-t0,kind=8)

  !comparison between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,overlap(1,1,1),overlap(1,1,2),&
       orbs%norbp**2,ntimes,.false.,Rgemm)


  nvctrp=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f

  call cpu_time(t0)
  do j=1,ntimes
     call dsyrk('L','T',orbs%norbp,nvctrp,1.0_wp,psi(1,1),nvctrp,0.0_wp,&
          overlap(1,1,1),orbs%norbp)
  end do
  call cpu_time(t1)
  CPUtime=real(t1-t0,kind=8)


  call cpu_time(t0)
  do j=1,ntimes
     call syrk('L','T',orbs%norbp,nvctrp,1.0_wp,psi(1,1),nvctrp,0.0_wp,&
          overlap(1,1,2),orbs%norbp)
  end do
  call cpu_time(t1)
  GPUtime=real(t1-t0,kind=8)

  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,overlap(1,1,1),overlap(1,1,2),&
       orbs%norbp**2,ntimes,.false.,Rsyrk)

  i_all=-product(shape(overlap))*kind(overlap)
  deallocate(overlap,stat=i_stat)
  call memocc(i_stat,i_all,'overlap',subname)


  !-------------------now the same for preconditioning
  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Preconditioner'
  !the input function is psi
  call cpu_time(t0)
  do j=1,ntimes
     call preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm)
  end do
  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)
  print *,'gnrm',gnrm


  !GPU data are already on the card, must be only copied back
  !the input function is GPU%hpsi in that case
  call cpu_time(t0)
  do j=1,ntimes
     call gpu_precond(lr,hx,hy,hz,GPU,orbs%norbp,ncong,&
          orbs%eval(min(orbs%isorb+1,orbs%norb)),gnrmGPU)
  end do
  call cpu_time(t1)
  
  GPUtime=real(t1-t0,kind=8)
  print *,'gnrmGPU',gnrmGPU

  do iorb=1,orbs%norbp
     call GPU_receive((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,&
          psi(1,(iorb-1)*orbs%nspinor+1),GPU%psi(iorb),i_stat)
  end do


  !free the card at the end
  call free_gpu(GPU,orbs%norbp)

  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,hpsi,psi,&
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),ntimes,.false.,Rprec)


  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)

  write(*,'(1x,a,5(1x,f7.3))')'Ratios:',Rden,Rham,Rgemm,Rsyrk,Rprec
  
end subroutine compare_cpu_gpu_hamiltonian


subroutine compare_data_and_gflops(CPUtime,GPUtime,GFlopsfactor,&
     CPUdata,GPUdata,n,ntimes,dowrite,ratio)
  use module_base
  implicit none
  logical, intent(in) :: dowrite
  integer, intent(in) :: n,ntimes
  real(gp), intent(in) :: CPUtime,GPUtime,GFlopsfactor
  real(gp), intent(out) :: ratio
  real(wp), dimension(n), intent(in) :: CPUdata,GPUdata
  !local variables
  integer :: i
  real(gp) :: CPUGflops,GPUGflops,maxdiff,comp

  CPUGflops=GFlopsfactor*real(ntimes,gp)/(CPUtime*1.d9)
  GPUGflops=GFlopsfactor*real(ntimes,gp)/(GPUtime*1.d9)

  maxdiff=0.0_gp

  rewind(17)

  do i=1,n
     if (dowrite) write(17,'(i6,2(1pe24.17))')i,CPUdata(i),GPUdata(i)
     comp=abs(CPUdata(i)-GPUdata(i))
     maxdiff=max(maxdiff,comp)
  end do
  ratio=CPUtime/GPUtime
  if (maxdiff <= 1.d-12) then
     write(*,'(1x,a,1x,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
          'GPU/CPU ratio,Time,Gflops: CPU,GPU',&
          CPUtime/GPUtime,maxdiff,&
          CPUtime*1.d3/real(ntimes,kind=8),CPUGflops,&
          GPUtime*1.d3/real(ntimes,kind=8),GPUGflops
  else
     write(*,'(1x,a,1x,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
          'GPU/CPU ratio,Time,Gflops: CPU,GPU',&
          CPUtime/GPUtime,maxdiff,&
          CPUtime*1.d3/real(ntimes,kind=8),CPUGflops,&
          GPUtime*1.d3/real(ntimes,kind=8),GPUGflops,&
          '<<<< WARNING' 
  end if


end subroutine compare_data_and_gflops
