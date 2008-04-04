subroutine MemoryEstimator(geocode,nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hx,hy,hz,nat,ntypes,&
     iatype,rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid,nspin,peakmem)

  use Poisson_Solver

  implicit none
  !Arguments
  character(len=1), intent(in) :: geocode
  logical, intent(in) :: output_grid
  integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin
  integer, dimension(nat), intent(in) :: iatype
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  real(kind=8), intent(in) :: hx,hy,hz,crmult,frmult,alat1,alat2,alat3
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
  real(kind=8), intent(out) :: peakmem
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f,norbp,nvctrp,i_all,i_stat
  integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,iat,i1,i2,i3
  real(kind=8) :: omemwf,omemker,omemden,omempot
  real(kind=8) :: tt,tmemker,tmemden,tmemps,tmemha,tminamount
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  ! determine localization region for all orbitals
  allocate(logrid_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_c))*kind(logrid_c),'logrid_c','memoryestimator')
  allocate(logrid_f(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_f))*kind(logrid_f),'logrid_f','memoryestimator')

  ! coarse grid quantities
  call fill_logrid(geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,1),crmult,hx,hy,hz,logrid_c)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,nvctr_c)

  ! fine grid quantities
  call fill_logrid(geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,2),frmult,hx,hy,hz,logrid_f)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,nvctr_f)

! Create the file grid.xyz to visualize the grid of functions
  if (output_grid) then
     open(unit=22,file='grid.xyz',status='unknown')
     write(22,*) nvctr_c+nvctr_f,' atomic'
     write(22,*)'complete simulation grid with low and high resolution points'
     do iat=1,nat
        write(22,'(a6,2x,3(1x,e12.5),3x)') &
             trim(atomnames(iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
     enddo
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  g ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     end do
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  G ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     enddo
     close(22)
  endif

  !here we must add the estimation for the projectors

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c','memoryestimator')
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f','memoryestimator')


  tt=dble(norb)/dble(nproc)
  norbp=int((1.d0-eps_mach*tt) + tt)
  tt=dble(nvctr_c+7*nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)

  !wavefunction memory per orbitals
  omemwf=real(nvctrp*nproc*8,kind=8)
  
  if (geocode == 'P') then
     call F_FFT_dimensions(2*n1+2,2*n2+2,2*n3+2,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc)
     n01=n1
     n02=n2
     n03=n3
     tt=1.d0
  else if (geocode == 'S') then
     call S_FFT_dimensions(2*n1+2,2*n2+31,2*n3+2,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc)
     n01=n1
     n02=2*n2+31
     n03=n3
     tt=real(2*n2,kind=8)/real(n02,kind=8)
  else if (geocode == 'F') then
     call F_FFT_dimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc)
     n01=2*n1+31
     n02=2*n2+31
     n03=2*n3+31
     tt=8.d0*real(n1*n2*n3,kind=8)/real(n01*n02*n03,kind=8)
  end if

  !density memory
  omemden=real(md3*md2/nproc,kind=8)*8.d0*real(md1*nspin,kind=8)
  !kernel memory
  omemker=real(nd2*nd3/nproc,kind=8)*8.d0*real(nd1,kind=8)
  !memory of full grid arrays
  omempot=real(n02*n03,kind=8)*8.d0*real(n01*nspin,kind=8)

  write(*,'(1x,a)')&
       '------------------------------------------------------------------ Memory Estimation'
  write(*,'(1x,a,i5,a,i6,a,3(i5))')&
       'Number of atoms=',nat,' Number of orbitals=',norb,' Sim. Box Dimensions= ', n1,n2,n3
  write(*,'(1x,a,i0,a)')&
       'Estimation performed for ',nproc,' processors.'
  write(*,'(1x,a)')&
       'Memory occupation for principal arrays:'
  write(*,'(1x,a,2(i6,a3))')&
       '              Poisson Solver Kernel (K):',int(omemker/1048576.d0),'MB',&
       ceiling((omemker-aint(omemker/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  write(*,'(1x,a,2(i6,a3))')&
       '             Poisson Solver Density (D):',int(omemden/1048576.d0),'MB',&
       ceiling((omemden-aint(omemden/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  write(*,'(1x,a,2(i6,a3))')&
       '    Single Wavefunction for one orbital:',int(omemwf/1048576.d0),'MB',&
       ceiling((omemwf-aint(omemwf/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  if(nproc > 1 ) omemwf=24.d0*real(norbp*nvctrp*nproc,kind=8)
  if(nproc == 1 ) omemwf=16.d0*real(norbp*nvctrp*nproc,kind=8)
  write(*,'(1x,a,2(i6,a3))')&
       '   All Wavefunctions for each processor:',int(omemwf/1048576.d0),'MB',&
       ceiling((omemwf-aint(omemwf/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  if(nproc > 1 ) omemwf=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
  if(nproc == 1 ) omemwf=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
  write(*,'(1x,a,2(i6,a3))')&
       '      Wavefunctions + DIIS per proc (W):',int(omemwf/1048576.d0),'MB',&
       ceiling((omemwf-aint(omemwf/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  write(*,'(1x,a,2(i6,a3))')&
       '   Arrays of full uncompressed grid (U):',int(omempot/1048576.d0),'MB',&
       ceiling((omempot-aint(omempot/1048576.d0)*1048576.d0)/1024.d0),'KB'  

  write(*,'(1x,a)')&
       'Estimation of Memory requirements for principal code sections:'
  write(*,'(1x,a)')&
       ' Kernel calculation | Density Construction | Poisson Solver | Hamiltonian application'
  if (nproc > 1) then 
     write(*,'(1x,a)')&
       '      ~19*K         |      W+(~4)*U+D+K    |    ~12*D+K+W   |      ~W+(~5)*U+D+K '
     tmemker=19.d0*omemker
     tmemden=omemwf+(3.d0+tt)*omempot+omemker
     tmemps=12.d0*omemden+omemwf+omemker
     tmemha=(3.d0+2.d0*tt)*omempot+omemwf+omemker
  else
     write(*,'(1x,a)')&
       '      ~11*K         |       ~W+(~3)*U      |    ~8*D+K+W    |         ~W+(~3)*U '
     tmemker=11.d0*omemker
     tmemden=omemwf+(2.d0+tt)*omempot+omemker
     tmemps=8.d0*omemden+omemwf+omemker
     tmemha=(2.d0+2.d0*tt)*omempot+omemwf+omemker
  end if
  write(*,'(1x,4(1x,i8,a))')&
       int(tmemker/1048576.d0),'MB         | ',int(tmemden/1048576.d0),'MB          |',&
       int(tmemps/1048576.d0),'MB     |     ',int(tmemha/1048576.d0),'MB'
  !estimation of the memory peak
  peakmem=max(tmemker,tmemden,tmemps,tmemha)
  write(*,'(1x,a,i0,a)')&
       'The overall memory requirement needed for this calculation is thus: ',&
       int(max(tmemker,tmemden,tmemps,tmemha)/1048576.d0),' MB'
  tminamount=real(3*(nvctr_c+7*nvctr_f)*8,kind=8)+3.d0*real(n01*n02,kind=8)+&
       (3.d0+2.d0*tt)*omempot
  write(*,'(1x,a)')&
       'By reducing the DIIS history and/or increasing the number of processors the amount of'
  write(*,'(1x,a,i0,a)')&
       ' memory can be reduced but for this system it will never be less than ',&
       int(tminamount/1048576.d0),' MB'

end subroutine MemoryEstimator
