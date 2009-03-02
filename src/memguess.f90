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
  logical :: calc_tail,optimise
  integer :: ierror,nproc,i_stat,i_all,output_grid
  integer :: nelec
  integer :: norbe,norbsc,nvctrp,nspin,iorb,norbu,norbd,nspinor,norb
  integer :: iunit,ityp
  real(kind=8) :: peakmem,hx,hy,hz
  type(input_variables) :: in
  type(atoms_data) :: atoms
  type(orbitals_data) :: orbs,orbstst
  type(locreg_descriptors) :: Glr
  type(nonlocal_psp_descriptors) :: nlpspd
  logical, dimension(:,:,:), allocatable :: logrid
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(kind=8), dimension(:,:), allocatable, target :: rxyz
  real(kind=8), dimension(:,:), allocatable :: radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psiat
  real(kind=8), dimension(:,:), allocatable :: xp, occupat
  integer, dimension(:,:), allocatable :: nl
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng
  real(kind=8), dimension(:), allocatable :: occup,spinsgn,locrad

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
  read(99,*) atoms%nat,atoms%units
 
  allocate(rxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',subname)

  !read atomic positions
  call read_atomic_positions(0,99,atoms,rxyz)

  close(99)

  !new way of reading the input variables, use structures
  call read_input_variables(0,'input.dat',in)

  call print_input_parameters(in,atoms)

  write(*,'(1x,a)')&
       '------------------------------------------------------------------ System Properties'
 

  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call system_properties(0,nproc,in,atoms,orbs,radii_cf,nelec)

  i_all=-product(shape(orbs%norb_par))*kind(orbs%norb_par)
  deallocate(orbs%norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'orbs%norb_par',subname)

  if (optimise) then
     call optimise_volume(atoms,in%crmult,in%frmult,in%hgrid,rxyz,radii_cf)
     write(*,'(1x,a)')'Writing optimised positions in file posout_000.xyz...'
     call wtposout(0,0.d0,rxyz,atoms)
  end if

  
  !in the case in which the number of orbitals is not "trivial" check whether they are too many
  if ( max(orbs%norbu,orbs%norbd) /= ceiling(real(nelec,kind=4)/2.0) .and. in%nspin /=4 ) then
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


     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(0,ngx,xp,psiat,occupat,ng,nl,atoms,norbe,norbsc,in%nspin,&
          scorb,norbsc_arr,locrad)

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
     if (in%nspin==1) then
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
  hx=in%hgrid
  hy=in%hgrid
  hz=in%hgrid

  call system_size(0,atoms,rxyz,radii_cf,in%crmult,in%frmult,hx,hy,hz,Glr)

  ! De-allocations
  i_all=-product(shape(orbs%occup))*kind(orbs%occup)
  deallocate(orbs%occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup',subname)
  i_all=-product(shape(orbs%spinsgn))*kind(orbs%spinsgn)
  deallocate(orbs%spinsgn,stat=i_stat)
  call memocc(i_stat,i_all,'spinsgn',subname)


  if (GPUconv .and. atoms%geocode=='P') then
     !test the hamiltonian in CPU or GPU
     !create the orbitals data structure for one orbital
     allocate(orbstst%norb_par(0:0+ndebug),stat=i_stat)
     call memocc(i_stat,orbstst%norb_par,'orbstst%norb_par',subname)
     !test orbitals
     nspin=1
     norb=orbs%norb
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

     call createWavefunctionsDescriptors(0,nproc,hx,hy,hz,&
          atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr,orbstst,nvctrp)
     
     call compare_cpu_gpu_hamiltonian(0,1,atoms,orbstst,nspin,in%ncong,in%ixc,&
          Glr,hx,hy,hz,rxyz)
     
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
subroutine optimise_volume(atoms,crmult,frmult,hgrid,rxyz,radii_cf)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), intent(in) :: crmult,frmult,hgrid
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  !local variables
  character(len=*), parameter :: subname='optimise_volume'
  integer :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,n1i,n2i,n3i,iat,i_all,i_stat,it,i
  real(gp) :: x,y,z,vol,tx,ty,tz,tvol,s,diag,dmax
  type(locreg_descriptors) :: Glr
  real(gp), dimension(3,3) :: urot
  real(gp), dimension(:,:), allocatable :: txyz

  allocate(txyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)

  call system_size(1,atoms,rxyz,radii_cf,crmult,frmult,hgrid,hgrid,hgrid,Glr)
  !call volume(nat,rxyz,vol)
  vol=atoms%alat1*atoms%alat2*atoms%alat3
  write(*,'(1x,a,1pe16.8)')'Initial volume (Bohr^3)',vol

  it=0
  diag=1.d-2 ! initial small diagonal element allows for search over all angles
  loop_rotations: do  ! loop over all trial rotations
     diag=diag*1.0001_gp ! increase diag to search over smaller angles
     it=it+1
     if (diag.gt.100._gp) exit loop_rotations ! smaller angle rotations do not make sense

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

     call system_size(1,atoms,txyz,radii_cf,crmult,frmult,hgrid,hgrid,hgrid,Glr)
     tvol=atoms%alat1*atoms%alat2*atoms%alat3
     !call volume(nat,txyz,tvol)
     if (tvol.lt.vol) then
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


subroutine compare_cpu_gpu_hamiltonian(iproc,nproc,at,orbs,nspin,ixc,ncong,&
     lr,hx,hy,hz,rxyz)
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  implicit none
  integer, intent(in) :: iproc,nproc,nspin,ncong,ixc
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables
  character(len=*), parameter :: subname='compare_cpu_gpu_hamiltonian'
  logical :: rsflag
  integer :: icoeff,norb,norbu,norbd,nspinor,i_stat,i_all,i1,i2,i3,ispin,j,ntimes
  integer :: iorb,n3d,n3p,n3pi,i3xcsh,i3s,jproc,nrhotot,nspinn
  real(kind=4) :: tt,t0,t1
  real(wp) :: maxdiff,comp
  real(gp) :: ttd,x,y,z,r2,arg,sigma2,ekin_sum,epot_sum,ekinGPU,epotGPU,gnrm,gnrmGPU
  real(kind=8) :: CPUtime,CPUGflops,GPUtime,GPUGflops
  type(gaussian_basis) :: G
  type(GPU_pointers) :: GPU
  integer, dimension(:,:), allocatable :: nscatterarr
  real(wp), dimension(:,:,:,:), allocatable :: pot,psig,rho
  real(wp), dimension(:,:), allocatable :: gaucoeffs,psi,hpsi

  ntimes=1

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

  !gaucoeffs(1,1)=1.0_wp

  !allocate the wavefunctions
  allocate(psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  call razero(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f*orbs%nspinor*orbs%norbp,psi)

  !convert the gaussians in wavelets
  call gaussians_to_wavelets(iproc,nproc,at%geocode,orbs,lr%d,&
       hx,hy,hz,lr%wfd,G,gaucoeffs,psi)

  psi=1.d0/sqrt(real(lr%d%n1i*lr%d%n2i*lr%d%n3i,wp))

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

  call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,pot,nproc)

  !for each of the orbitals treated by the processor build the partial densities
  call cpu_time(t0)
  call local_partial_density(iproc,nproc,rsflag,nscatterarr,&
     nrhotot,lr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,orbs%nspinor,orbs%norbp,&
     orbs%occup(min(orbs%isorb+1,orbs%norb)),orbs%spinsgn(min(orbs%isorb+1,orbs%norb)),&
     psi,pot)
  call cpu_time(t1)
  CPUtime=real(t1-t0,kind=8)

  call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho,nproc)

  !now the GPU part
  !for each of the orbitals treated by the processor build the partial densities
  call cpu_time(t0)
  call local_partial_density(iproc,nproc,rsflag,nscatterarr,&
     nrhotot,lr,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,nspin,orbs%nspinor,orbs%norbp,&
     orbs%occup(min(orbs%isorb+1,orbs%norb)),orbs%spinsgn(min(orbs%isorb+1,orbs%norb)),&
     psi,rho)
  call cpu_time(t1)
  GPUtime=real(t1-t0,kind=8)


  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)

  !compare the results between the different actions of the hamiltonian
  !check the differences between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,pot,rho,&
       lr%d%n1i*lr%d%n2i*lr%d%n3i,ntimes)

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

  !set initialisation of GPU part 
  call prepare_gpu_for_locham(lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,lr%wfd,orbs,GPU)

  !copy the wavefunctions and the potential on GPU
  do iorb=1,orbs%norbp
     call GPU_send((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,&
          psi(1,(iorb-1)*orbs%nspinor+1),GPU%psi(iorb),i_stat)
  end do
  call GPU_send(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin,pot,GPU%pot,i_stat)


  write(*,'(1x,a)')repeat('-',34)//' CPU-GPU comparison: Local Hamiltonian calculation'

  !warm-up
  call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum) 

  !apply the CPU hamiltonian
  !take timings
  call cpu_time(t0)
  do j=1,ntimes
     call local_hamiltonian(iproc,orbs,lr,hx,hy,hz,nspin,pot,psi,hpsi,ekin_sum,epot_sum) 
  end do
  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)

  print *,'ekin,epot=',ekin_sum,epot_sum

  !warm-up
  call gpu_locham(lr%d%n1,lr%d%n2,lr%d%n3,hx,hy,hz,orbs,GPU,ekinGPU,epotGPU)

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
          psi(1,(iorb-1)*orbs%nspinor+1),GPU%hpsi(iorb),i_stat)
  end do
  

  !compare the results between the different actions of the hamiltonian
  !check the differences between the results
  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,hpsi,psi,&
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),ntimes)


  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)

  !-------------------now the same for preconditioning
  !the input function is psi
  call cpu_time(t0)
  do j=1,ntimes
     call preconditionall(iproc,nproc,orbs%norbp,lr,hx,hy,hz,ncong,orbs%nspinor,&
          orbs%eval(min(orbs%isorb+1,orbs%norb)),hpsi,gnrm)
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
          psi(1,(iorb-1)*orbs%nspinor+1),GPU%hpsi(iorb),i_stat)
  end do


  !free the card at the end
  call free_gpu(GPU,orbs%norbp)

  call compare_data_and_gflops(CPUtime,GPUtime,&
       8.d0*real(lr%d%n1*lr%d%n2*lr%d%n3,kind=8)*366.d0,hpsi,psi,&
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),ntimes)


  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)

  
end subroutine compare_cpu_gpu_hamiltonian


subroutine compare_data_and_gflops(CPUtime,GPUtime,GFlopsfactor,CPUdata,GPUdata,n,ntimes)
  use module_base
  implicit none
  integer, intent(in) :: n,ntimes
  real(gp), intent(in) :: CPUtime,GPUtime,GFlopsfactor
  real(wp), dimension(n), intent(in) :: CPUdata,GPUdata
  !local variables
  integer :: i
  real(gp) :: CPUGflops,GPUGflops,maxdiff,comp

  CPUGflops=GFlopsfactor*real(ntimes,gp)/(CPUtime*1.d9)
  GPUGflops=GFlopsfactor*real(ntimes,gp)/(GPUtime*1.d9)

  maxdiff=0.0_gp

  rewind(17)

  do i=1,n
     write(17,'(i6,2(1pe24.17))')i,CPUdata(i),GPUdata(i)
     comp=abs(CPUdata(i)-GPUdata(i))
     maxdiff=max(maxdiff,comp)
  end do
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
