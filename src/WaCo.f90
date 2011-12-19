program WaCo

   use module_base
   use module_types
   use module_interfaces, except_this_one => writeonewave
   use Poisson_Solver
   implicit none
   character :: filetype*4
   type(locreg_descriptors) :: Glr
   type(orbitals_data) :: orbs,orbsw,orbsv
   type(atoms_data) :: atoms
   type(input_variables) :: input
   type(workarr_sumrho) :: w
   type(communications_arrays), target :: comms, commsp,commsv,commsw
   real(gp), parameter :: b2a=0.5291772108
   real :: tcpu0,tcpu1,tel
   real(gp), dimension(3) :: shift
   real(gp) :: dist
   integer :: iproc, nproc, nproctiming, i_stat, nelec, ind, ierr, npsidim, npsidim2
   integer :: n_proj,nvctrp,npp,nvirtu,nvirtd,pshft,nbl1,nbl2,nbl3,nrpts
   integer :: ncount0,ncount1,ncount_rate,ncount_max,nbr1,nbr2,nbr3,iat,iformat
   integer :: int1, int2, int3, iwann, iiwann, iband, nwann, nband, plotwann
   integer :: ix, iy, iz, iw1, iw2, rem, ntot, isplotwann, plotwannp, jproc, ifile
   integer :: iter, nbuf,nsprd,ndiag
   character(len=*), parameter :: subname='WaCo'
   character(len=4) :: num, char1
   integer, allocatable :: wann_list(:), Zatoms(:), ncenters(:)
   real(gp), dimension(:,:), pointer :: rxyz, rxyz_old, cxyz,rxyz_wann
   real(gp), dimension(:,:), allocatable :: radii_cf
   real(gp), allocatable :: sprd(:)
   real(wp), allocatable :: psi(:,:),wann(:),wannr(:)
   real(wp), allocatable :: ham(:,:,:),hamr(:,:,:),Uham(:),evectors(:,:),evalues(:),work(:),iwork(:)
   real(wp), allocatable :: diag(:)
   character(len=60) :: radical, filename
   character(len=8) :: string
   logical :: perx, pery, perz, commented, same
   integer, dimension(:), allocatable :: plotwann_par,isplotwann_par
   !cube
   integer :: nx, ny, nz, nb, nb1, nb2, nk, inn
   integer, allocatable, dimension(:) :: Z
   real(kind=8) :: bx(3), by(3), bz(3), b1, b2, b3, r0x, r0y, r0z, vec(18)
   real(kind=8) :: xx, yy, zz
   real(kind=8), allocatable :: at_pos(:,:)
   real(kind=8), allocatable :: umn(:,:), rho(:,:), rhoprime(:,:)
   integer :: i, j, k, np,i_all
   character :: seedname*16
   logical :: calc_only_A
   real, dimension(3,3) :: real_latt, recip_latt
   integer :: n_kpts, n_poj, n_nnkpts, n_excb, n_at, n_bands, s
   integer :: n_occ, n_virt, n_virt_tot
   logical :: w_unk, w_sph, w_ang, w_rad, pre_check
   real, allocatable, dimension (:,:) :: kpts
   real(kind=8), allocatable, dimension (:,:) :: ctr_proj, x_proj, y_proj, z_proj
   integer, allocatable, dimension (:) :: l, mr, rvalue
   real, allocatable, dimension (:) :: zona
   integer, allocatable, dimension (:,:) :: k_plus_b
   integer, allocatable, dimension (:,:) :: G_vec
   integer, allocatable, dimension (:) :: excb
   integer, allocatable, dimension (:) :: virt_list, amnk_bands_sorted
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   logical :: file_exist,idemp
   interface
      subroutine read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)
        implicit none
        integer, intent(in) :: iproc
        character, intent(out) :: seedname*16, filetype*4
        integer, intent(out) :: n_occ, n_virt, n_virt_tot
        logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check
      end subroutine read_inter_header
      subroutine read_inter_list(iproc,n_virt, virt_list)
        implicit none
        integer, intent(in) :: n_virt,iproc
        integer, dimension(n_virt), intent(out) :: virt_list
      end subroutine read_inter_list
      subroutine scalar_kmeans_diffIG(nIG,crit,nel,vect,string,nbuf)
        implicit none
        integer, intent(in) :: nel,nIG
        real(kind=8),intent(in) :: crit
        real(kind=8), dimension(nel),intent(in) :: vect
        character(len=*),intent(in) :: string
        integer, intent(out) :: nbuf
      end subroutine scalar_kmeans_diffIG
   end interface

   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

   call memocc_set_memory_limit(memorylimit)

   ! Read a possible radical format argument.
   call get_command_argument(1, value = radical, status = i_stat)
   if (i_stat > 0) then
      write(radical, "(A)") "input"
   end if

   if (input%verbosity > 2) then
      nproctiming=-nproc !timing in debug mode                                                                                                                                                                 
   else
      nproctiming=nproc
   end if

   call timing(nproctiming,'WaCo_time.prc','IN')

   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max) 
   call timing(iproc,'Precondition  ','ON')   

   !###################################################################
   ! Initialise the variables for the wavefunctions
   !###################################################################
   call standard_inputfile_names(input,radical)
   call read_input_variables(iproc,'posinp',input, atoms, rxyz)

   if (iproc == 0) call print_general_parameters(nproc,input,atoms)

   allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)

   call system_properties(iproc,nproc,input,atoms,orbs,radii_cf,nelec)

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space
   call system_size(iproc,atoms,rxyz,radii_cf,input%crmult,input%frmult,input%hx,input%hy,input%hz,&
        Glr,shift)

   ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
   call createWavefunctionsDescriptors(iproc,input%hx,input%hy,input%hz,&
        atoms,rxyz,radii_cf,input%crmult,input%frmult,Glr)

   ! don't need radii_cf anymore
   i_all = -product(shape(radii_cf))*kind(radii_cf)
   deallocate(radii_cf,stat=i_stat)
   call memocc(i_stat,i_all,'radii_cf',subname)

   !#################################################################
   ! Read the input.inter file
   !#################################################################
   call read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)
   allocate(virt_list(n_virt),stat=i_stat)
   call memocc(i_stat,virt_list,'virt_list',subname)
   if (n_virt .ne. 0) then
      call read_inter_list(iproc, n_virt, virt_list)
   end if 

   !#####################################################################
   ! Read the Umn file
   !#####################################################################
   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .umn :               !' 
      write(*,*) '!==================================!'
   end if

   inquire(file=trim(seedname)//'.umn',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'.umn, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   open(11, file=trim(seedname)//'.umn', status='OLD')
   read(11,*) nwann, nband

   allocate(umn(nwann,nband),stat=i_stat)
   call memocc(i_stat,umn,'umn',subname)

   if(nband .ne. n_occ+n_virt) then
     if(iproc == 0) then
        write(*,'(A,1x,i4,1x,A)') 'ERROR : Number of bands in the input.inter',n_occ+n_virt,' not equal to'
        write(*,'(A,1x,i4)') 'the number of bands used in the Wannier construction:', nband 
     end if
     call mpi_finalize(i)
     stop 
   end if

   do iwann = 1,nwann
      do iband = 1,nband
         read(11,*) int1, int2, umn(iwann,iband)
      end do
   end do
   close(11)

   if(iproc==0) then
      write(*,*) 'Number of Wannier functions : ',nwann
      write(*,*) 'Number of orbitals used     : ',nband
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .umn : DONE          !' 
      write(*,*) '!==================================!'
   end if
   call timing(iproc,'Precondition  ','OF')

   !###################################################################
   ! Constructing the density matrix
   !###################################################################

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Constructing density matrix    !'
      write(*,*) '!==================================!'
      write(*,*)
   end if

   allocate(rho(nwann,nwann),stat=i_stat)
   call memocc(i_stat,rho,'rho',subname)
   call to_zero(nwann*nwann, rho(1,1))

   do i=1,nwann
      do j=1,nwann
         do iband=1,nband
            rho(i,j) = rho(i,j) + umn(i,iband)*umn(j,iband)
          end do
      end do
   end do

!!   if(iproc==0) then
!!      do i=1,nwann
!!         do j=1,nwann
!!            write(*,*) i , j ,rho(i,j)
!!         end do
!!      end do
!!   end if

   if(iproc==0) then
      write(*,*) '!===================================!'
      write(*,*) '!Constructing density matrix : DONE !'
      write(*,*) '!===================================!'
      write(*,*)
   end if

   !##################################################################
   ! Check idempotence of density matrix
   !##################################################################
   allocate(rhoprime(nwann,nwann),stat=i_stat)
   call memocc(i_stat,rhoprime,'rhoprime',subname)
   call to_zero(nwann*nwann, rhoprime(1,1))
   
   idemp = .true.
   do i = 1, nwann
      do j = 1, nwann
         do k = 1, nwann
            rhoprime(i,j) = rhoprime(i,j) + rho(i,k)*rho(k,j)
         end do
         if(abs(rhoprime(i,j) - rho(i,j))>1.0d-5) then
           write(*,*) 'Not indempotent',i,j,rhoprime(i,j), rho(i,j)
           idemp = .false.
         end if
      end do
   end do

   if(idemp .eqv. .false.) then
     stop 'Density matrix not idempotent'
   end if

   i_all = -product(shape(rhoprime))*kind(rhoprime)
   deallocate(rhoprime,stat=i_stat)
   call memocc(i_stat,i_all,'rhoprime',subname)
   i_all = -product(shape(rho))*kind(rho)
   deallocate(rho,stat=i_stat)
   call memocc(i_stat,i_all,'rho',subname)

   !#########################################################
   ! Read seedname.dat to make the Wannier plot list
   !#########################################################
   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .dat :               !'
      write(*,*) '!==================================!'
   end if

   inquire(file=trim(seedname)//'.dat',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'.dat, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   !wann_list will contain the list of occupied Wannier functions
   allocate(wann_list(nwann),stat=i_stat)
   call memocc(i_stat,wann_list,'wann_list',subname)
   wann_list = 0

   open(11, file=trim(seedname)//'.dat', status='OLD')

   ! Check if the line is commented
   plotwann = 0
   do i = 1 ,nwann
      read(11,*) string
      if(VERIFY(string, ' #', .false.) == 0) cycle
      plotwann = plotwann +1
      wann_list(plotwann) = i
   end do

   close(11)

   if(iproc == 0) then
      write(*,'(A)') 'Occupied/Plotted Wannier functions are:'
      do i=1, plotwann, 6
         write(*,'(6(2x,I4))') (wann_list(j), j=i,i+5)
      end do
   end if

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .dat : DONE          !'
      write(*,*) '!==================================!'
   end if

   !########################################################
   ! Bonding analysis
   !########################################################

  if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Bonding Analysis :           !'
      write(*,*) '!==================================!'
   end if

   !read the _centres.xyz file
   inquire(file=trim(seedname)//'_centres.xyz',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_centres.xyz, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   !wann_list will contain the list of occupied Wannier functions
   allocate(cxyz(3,plotwann),stat=i_stat)
   call memocc(i_stat,cxyz,'cxyz',subname)
   allocate(rxyz_wann(3,atoms%nat),stat=i_stat)
   call memocc(i_stat,rxyz_wann,'rxyz_wann',subname)
   allocate(sprd(plotwann),stat=i_stat)
   call memocc(i_stat,sprd,'sprd',subname)
   allocate(Zatoms(atoms%nat),stat=i_stat)
   call memocc(i_stat,Zatoms,'Zatoms',subname)  

   open(11, file=trim(seedname)//'_centres.xyz', status='OLD')

   !skip first two lines
   read(11,*)
   read(11,*)

   !now read the centers
   iiwann = 0
   do iwann = 1, nwann
      commented = .true.
      do i = 1, plotwann
         if(iwann == wann_list(i)) commented = .false.
      end do
      if (.not. commented) then
         iiwann = iiwann + 1
         read(11,*) char1, cxyz(1,iiwann), cxyz(2,iiwann), cxyz(3,iiwann)
      end if
      if (commented) read(11,*)  !just skip this line
   end do

   !now read atomic positions (already in rxyz but beware shifts...)
   do i = 1, atoms%nat
      read(11,*) char1, rxyz_wann(1,i), rxyz_wann(2,i), rxyz_wann(3,i)
   end do

   ! now read the spreads (the additionnal last value is total sprd)
   do iiwann = 1, plotwann+1
      read(11, *) char1, sprd(iiwann)
   end do
   close(11)

   allocate(ncenters(plotwann),stat=i_stat)
   call memocc(i_stat,ncenters,'ncenters',subname)

   ! Now calculate the bonding distances and ncenters
   if(iproc == 0) write(*,'(2x,A)') 'Number of atoms associated to the WFs'
   if(iproc == 0) write(*,'(3x,A,4x,A,3x,A,3x,A)') 'WF','Spr(ang^2)','Nc','Atom numbers:'
   do iwann = 1, plotwann
      ncenters(iwann) = 0.0d0
      iat = 0
      Zatoms = 0
      do i = 1, atoms%nat
         dist = (rxyz_wann(1,i)-cxyz(1,iwann))**2 + (rxyz_wann(2,i)-cxyz(2,iwann))**2 + (rxyz_wann(3,i)-cxyz(3,iwann))**2
         if (dist <= 1.64 * sprd(iwann)) then    !for normal distribution: 1=68%, 1.64=80%, 3=94%
            ncenters(iwann) = ncenters(iwann) +1
            iat = iat +1
            Zatoms(iat) = i
         end if
      end do
      if(iproc == 0) then
        write(*,'(I4,F14.6,2x,I4,6(2x,I4))') iwann, sprd(iwann), ncenters(iwann), (Zatoms(i_all),i_all=1,iat)
      end if
   end do

   call scalar_kmeans_diffIG(0,maxval(sprd)*1.0d-1,plotwann,sprd,'spread',nsprd)

   i_all = -product(shape(sprd))*kind(sprd)
   deallocate(sprd,stat=i_stat)
   call memocc(i_stat,i_all,'sprd',subname)
   i_all = -product(shape(Zatoms))*kind(Zatoms)
   deallocate(Zatoms,stat=i_stat)
   call memocc(i_stat,i_all,'Zatoms',subname)
   i_all = -product(shape(rxyz_wann))*kind(rxyz_wann)
   deallocate(rxyz_wann,stat=i_stat)
   call memocc(i_stat,i_all,'rxyz_wann',subname)


  if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Bonding Analysis : DONE      !' 
      write(*,*) '!==================================!'
   end if


   !##########################################################################
   ! Hamiltonian stuff
   !##########################################################################
   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Hamiltonian Analysis :       !' 
      write(*,*) '!==================================!'
   end if

   !read the _centres.xyz file
   inquire(file=trim(seedname)//'_hr.dat',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_hr.dat, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   !read the hamiltonian
   open(11, file=trim(seedname)//'_hr.dat', status='OLD')

   !skip first lines
   read(11,*)           !comment line containing the date 
   read(11,*)           !this line contains the number of wannier functions (already known)
   read(11,*) nrpts     ! this line contains the number of Wigner-Seitz grid points
   do i=1,nrpts,15
      read(11,*)        !skip the block of integers giving the degeneracy of the Wigner-Seitz grid points (15 integers per line)
   end do
  
!  if(iproc==0) write(*,*) 'This calculation needs:',(nrpts*nwann*nwann+nwann*(nwann+1)/2+1+&
!                6*nwann+nwann**2+3+5*nwann+nwann**2+nwann)*kind(work),'bytes'

   allocate(hamr(nrpts,nwann,nwann),stat=i_stat)
   call memocc(i_stat,hamr,'hamr',subname)

   do i = 1, nrpts 
      do iiwann = 1, nwann
         do iwann = 1, nwann
            read(11,*) int1, int2, int3, iw1, iw2, hamr(i,iw1,iw2)
         end do
      end do
   end do
   close(11)

   allocate(ham(nrpts,nwann,nwann),stat=i_stat)
   call memocc(i_stat,ham,'ham',subname)

   !Eliminate the unoccupied states
!   write(*,'(A)') 'Hamiltonian (without empty WFs)'
!   write(*,'(3x,A,4x,A,6x,A,6x,A,3x,A,4x,A)')'WF','WF','R[ham]','Nc(1)','Nc(2)','Dist[bohr]'
   do i = 1, nrpts
     do iiwann = 1, plotwann
           iw1 = wann_list(iiwann)
        do iwann = 1, plotwann
            iw2 = wann_list(iwann)
            dist = (cxyz(1,iiwann) - cxyz(1,iwann))**2 + (cxyz(2,iiwann) - cxyz(2,iwann))**2 + (cxyz(3,iiwann) - cxyz(3,iwann))**2
            ham(i,iiwann,iwann) = hamr(i,iw1,iw2)
!            if(abs(hamr(i,iw1,iw2)) > maxval(hamr)*1.0d-2) write(*,'(i4,2x,i4,2x,E14.6,2x,i4,2x,i4,2x,F14.6)'),iiwann,iwann,&
!                   hamr(i,iw1,iw2),ncenters(iiwann),ncenters(iwann),dist
        end do
     end do
   end do

  write(*,'(A)') 'Diagonal of the Hamiltonian (without empty WFs)'
  allocate(diag(plotwann),stat=i_stat)
  call memocc(i_stat,diag,'diag',subname)
  do i = 1, nrpts
     do iwann = 1, plotwann
        write(*,'(i4,2x,i4,2x,E14.6,2x,i4)')iwann,iwann,ham(i,iwann,iwann),ncenters(iwann)
        diag(iwann) = ham(i,iwann,iwann)
     end do
   end do

  call scalar_kmeans_diffIG(nsprd,2.0d-2,plotwann,diag,'diagonal',ndiag)

!! Diagonalizing the Hamiltonian yields the Kohn-Sham eigenvalues and orbitals (U^-1)
!!   allocate(Uham(nwann*(nwann+1)/2),stat=i_stat)
!!   call memocc(i_stat,Uham,'Uham',subname)
!!
!!   ! Copy the upper half of the hamiltonian
!!   !do i = 1, nrpts     ! Right now only one cell
!!   ind = 0
!!   do iiwann = 1, nwann
!!      do iwann = iiwann, nwann
!!           ind = ind +1
!!           Uham(ind) = hamr(1,iiwann,iwann)
!!           if(iproc==0)print *,iiwann,iwann,hamr(1,iiwann,iwann)
!!      end do
!!   end do
!!
!!   allocate(evalues(nwann),stat=i_stat)
!!   call memocc(i_stat,evalues,'evalues',subname)
!!   allocate(evectors(nwann,nwann),stat=i_stat)
!!   call memocc(i_stat,evectors,'evectors',subname)
!!   allocate(work(1+6*nwann+nwann**2),stat=i_stat)
!!   call memocc(i_stat,work,'work',subname)
!!   allocate(iwork(3+5*nwann),stat=i_stat)
!!   call memocc(i_stat,iwork,'iwork',subname)
!!
!!   !Now we can solve the matrix: if(ierr == 0), evectors and evalues contains the eigenvectors and eigenvalues 
!!   call dspevd('V', 'L', nwann, Uham, evalues, evectors, nwann, work, 1+6*nwann+nwann**2, iwork, 3+5*nwann, ierr)
!!
!!   i_all=-product(shape(work))*kind(work)
!!   deallocate(work,stat=i_stat)
!!   call memocc(i_stat,i_all,'work',subname)
!!   i_all=-product(shape(iwork))*kind(iwork)
!!   deallocate(iwork,stat=i_stat)
!!   call memocc(i_stat,i_all,'iwork',subname)
!!
!!
!!   !Print the eigenvalues and eigenvectors
!!   if(iproc == 0) then
!!     do iwann = 1, nwann
!!        write(*,*) 'Eigenvalue:',evalues(iwann)
!!        write(*,*) evectors(iwann,:)
!!     end do
!!   end if
!!
!!   i_all=-product(shape(evalues))*kind(evalues)
!!   deallocate(evalues,stat=i_stat)
!!   call memocc(i_stat,i_all,'evalues',subname)
!!   i_all=-product(shape(evectors))*kind(evectors)
!!   deallocate(evectors,stat=i_stat)
!!   call memocc(i_stat,i_all,'evectors',subname)
!!   i_all=-product(shape(Uham))*kind(Uham)
!!   deallocate(Uham,stat=i_stat)
!!   call memocc(i_stat,i_all,'Uham',subname)
   i_all=-product(shape(hamr))*kind(hamr)
   deallocate(hamr,stat=i_stat)
   call memocc(i_stat,i_all,'hamr',subname)
   i_all=-product(shape(ham))*kind(ham)
   deallocate(ham,stat=i_stat)
   call memocc(i_stat,i_all,'ham',subname)
   i_all=-product(shape(diag))*kind(diag)
   deallocate(diag,stat=i_stat)
   call memocc(i_stat,i_all,'diag',subname)
   i_all = -product(shape(cxyz))*kind(cxyz)
   deallocate(cxyz,stat=i_stat)
   call memocc(i_stat,i_all,'cxyz',subname)


   if(iproc==0) then 
      write(*,*) '!==================================!'
      write(*,*) '!     Hamiltonian Analysis : DONE  !' 
      write(*,*) '!==================================!'
   end if
call mpi_finalize(ierr)
stop
   !###########################################################################
   ! Set-up number of states used in the Wannier construction
   !###########################################################################
   nvirtu = nband
   nvirtd = 0
   if (input%nspin==2) nvirtd=nvirtu
   call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
       & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsw)
   call orbitals_communicators(iproc,nproc,Glr,orbsw,commsw)

   nvirtu = n_virt
   nvirtd = 0
   if (input%nspin==2) nvirtd=nvirtu
   call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
       & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv)

   !##########################################################################
   ! Read the Wavefunctions
   !##########################################################################

   ! Wavefunctions calculated by BigDFT already are normalized.
   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!  Reading the wavefunctions       !'
      write(*,*) '!==================================!'
   end if

   call timing(iproc,'CrtProjectors ','ON')

   ! assign the input_wf_format
   iformat = WF_FORMAT_NONE
   select case (filetype)
   case ("ETSF","etsf")
      iformat = WF_FORMAT_ETSF
   case ("BIN","bin")
      iformat = WF_FORMAT_BINARY
   case default
      if (iproc == 0) write(*,*)' WARNING: Missing specification of wavefunction files'
      stop
   end select

   npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsw%norbp*orbsw%nspinor,sum(commsw%ncntt(0:nproc-1)))
   allocate(psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbsw%norbp*orbsw%nspinor),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)
   allocate(rxyz_old(3,atoms%nat),stat=i_stat)  
   call memocc(i_stat,rxyz_old,'rxyz_old',subname)

   ! For the occupied orbitals, need to modifify norbp,isorb to match the total distributed scheme
   orbs%norbp = n_occ - orbsw%isorb
   if (orbsw%isorb + orbsw%norbp < n_occ ) orbs%norbp = orbsw%norbp
   if(orbsw%isorb > n_occ) orbs%norbp = 0
   orbs%isorb = orbsw%isorb
   if(associated(orbs%iokpt)) then
      i_all = -product(shape(orbs%iokpt))*kind(orbs%iokpt)
      deallocate(orbs%iokpt,stat=i_stat)
      call memocc(i_stat,i_all,'orbs%iokpt',subname)
   end if
   allocate(orbs%iokpt(orbs%norbp),stat=i_stat)
   call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)
   orbs%iokpt=1
   if(orbs%norbp > 0) then
      if(associated(orbs%eval)) nullify(orbs%eval)
      allocate(orbs%eval(orbs%norb*orbs%nkpts), stat=i_stat)
      call memocc(i_stat,orbs%eval,'orbs%eval',subname)
      filename=trim(input%dir_output) // 'wavefunction'
      call readmywaves(iproc,filename,iformat,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
         Glr%wfd,psi(1,1))
      i_all = -product(shape(orbs%eval))*kind(orbs%eval)
      deallocate(orbs%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbs%eval',subname)
   
   end if

   ! For the non-occupied orbitals, need to change norbp,isorb
   orbsv%norbp = orbsw%isorb + orbsw%norbp - n_occ
   if (orbsw%isorb + orbsw%norbp < n_occ ) orbsv%norbp = 0
   if (orbsw%isorb > n_occ) orbsv%norbp = orbsw%norbp
   orbsv%isorb = 0
   if(orbsw%isorb >= n_occ) orbsv%isorb = orbsw%isorb - n_occ
   if(associated(orbsv%iokpt)) then
      i_all = -product(shape(orbsv%iokpt))*kind(orbsv%iokpt)
      deallocate(orbsv%iokpt,stat=i_stat)
      call memocc(i_stat,i_all,'orbsv%iokpt',subname)
   end if
   allocate(orbsv%iokpt(orbsv%norbp),stat=i_stat)
   call memocc(i_stat,orbsv%iokpt,'orbsv%iokpt',subname)
   orbsv%iokpt=1

   ! read unoccupied wavefunctions
   if(orbsv%norbp > 0) then
      filename=trim(input%dir_output) // 'virtuals'
      if(associated(orbsv%eval)) nullify(orbsv%eval)
      allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
      call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)
      call readmywaves(iproc,filename,iformat,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
         Glr%wfd,psi(1,1+orbs%norbp),virt_list)
      i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
      deallocate(orbsv%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbsv%eval',subname)
   end if

   i_all = -product(shape(rxyz_old))*kind(rxyz_old)
   deallocate(rxyz_old,stat=i_stat)
   call memocc(i_stat,i_all,'rxyz_old',subname)
   i_all = -product(shape(virt_list))*kind(virt_list)
   deallocate(virt_list,stat=i_stat)
   call memocc(i_stat,i_all,'virt_list',subname)


   call timing(iproc,'CrtProjectors ','OF')

   if (iproc==0) then
      write(*,*) '!===================================!'
      write(*,*) '!  Reading the wavefunctions : DONE !'
      write(*,*) '!===================================!'
   end if

   !#########################################################
   ! Construct the Wannier functions
   !#########################################################

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Constructing WFs :           !'
      write(*,*) '!==================================!'
   end if


   allocate(wann(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),stat=i_stat)
   call memocc(i_stat,wann,'wann',subname) 
   allocate(wannr(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i),stat=i_stat)
   call memocc(i_stat,wannr,'wannr',subname) 
   call initialize_work_arrays_sumrho(Glr,w)

   perx=(Glr%geocode /= 'F')
   pery=(Glr%geocode == 'P')
   perz=(Glr%geocode /= 'F')
   call ext_buffers(perx,nbl1,nbr1)
   call ext_buffers(pery,nbl2,nbr2)
   call ext_buffers(perz,nbl3,nbr3)
   if(nbr1 > 0) nbr1 = nbr1 + 2
   if(nbr2 > 0) nbr2 = nbr2 + 2
   if(nbr3 > 0) nbr3 = nbr3 + 2
   ! Volumetric data in batches of 6 values per line, 'z'-direction first.
   rem=Glr%d%n3i-floor(Glr%d%n3i/6.d0)*6

   ! Separate plotwann
!!   allocate(plotwann_par(0:nproc-1),stat=i_stat)
!!   call memocc(i_stat,plotwann_par,'plotwann_par',subname)
!!   allocate(isplotwann_par(0:nproc-1),stat=i_stat)
!!   call memocc(i_stat,isplotwann_par,'isplotwann_par',subname)
!!   call parallel_repartition_with_kpoints(nproc,1,plotwann,plotwann_par)
!!   ntot=0
!!   do jproc=0,nproc-1
!!      isplotwann_par(jproc)=ntot
!!      ntot=ntot+plotwann_par(jproc)
!!   end do
!!   isplotwann = isplotwann_par(iproc) 
!!   plotwannp = plotwann_par(iproc)
!!   i_all = -product(shape(plotwann_par))*kind(plotwann_par)
!!   deallocate(plotwann_par,stat=i_stat)
!!   call memocc(i_stat,i_all,'plotwann_par',subname)
!!   i_all = -product(shape(isplotwann_par))*kind(isplotwann_par)
!!   deallocate(isplotwann_par,stat=i_stat)
!!   call memocc(i_stat,i_all,'isplotwann_par',subname)


   ! Now construct the WFs
   ifile = 12 + iproc
   do iiwann = 1, plotwann
      iwann = wann_list(iiwann)
      if(iwann == 0) stop 'this should not happen' 
      call to_zero(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, wann(1))
      do iband = 1, orbsw%norbp
         do i = 1, Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
            wann(i) = wann(i) + umn(iwann,iband+orbsw%isorb) * psi(i,iband)
         end do
      end do
      ! Construction of the Wannier function.
      call mpiallred(wann(1),Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(iproc == 0) then
         if(.true.) then
            !Put it in interpolating scaling functions
            call daub_to_isf(Glr,w,wann(1),wannr)
            !write it to cube file
            write(num,'(i4.4)') iiwann
            open(ifile, file=trim(seedname)//'_'//num//'.cube', status='unknown')
            write(ifile,*) ' CUBE file for ISF field'
            write(ifile,*) ' Case for'
            write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') atoms%nat, real(0.d0), real(0.d0), real(0.d0)
            write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n1i-(nbl1+nbr1), 0.5*input%hx, real(0.d0),  real(0.d0)
            write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n2i-(nbl2+nbr2), real(0.d0),  0.5*input%hy, real(0.d0)
            write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n3i-(nbl3+nbr3), real(0.d0),  real(0.d0),  0.5*input%hz
            do i=1, atoms%nat
               write(ifile,'(I4,1X,F12.6,3(1X,F12.6))') atoms%nzatom(atoms%iatype(i)), real(0.d0), (real(rxyz(j,i)), j=1,3)
            end do
!            do ix=Glr%d%n1i-nbr1,1+nbl1,-1
!               do iy=Glr%d%n2i-nbr2,1+nbl2,-1
!                  do iz=Glr%d%n3i-nbr3,1+nbl3,-1
            do ix=1+nbl1,Glr%d%n1i-nbr1
               do iy=1+nbl2,Glr%d%n2i-nbr2
                  do iz=1+nbl3,Glr%d%n3i-nbr3
                     ind = iz*(Glr%d%n2i*Glr%d%n1i) + iy*Glr%d%n1i + ix
                     write(ifile,'(E14.6)',advance='no') real(wannr(ind))
                     if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. Glr%d%n3i) ) .or. (iz .eq. 1) ) then
                        write(ifile,'(a)') ''
                     end if
                  end do
               end do
            end do
            close(ifile)
         else
           if(.false.) then
              open(ifile, file=trim(seedname)//'_'//num//'.bin', status='unknown')
              call writeonewave(ifile,.false.,iiwann,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms%nat,rxyz,  & 
                   Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),  & 
                   Glr%wfd%nseg_f,Glr%wfd%nvctr_f,Glr%wfd%keyg(1,Glr%wfd%nseg_c+1),Glr%wfd%keyv(Glr%wfd%nseg_c+1), & 
                   wann(1),wann(Glr%wfd%nvctr_c+1), 0.d0)
           else
              ! should be write_wave_etsf  (only one orbital)
              call write_waves_etsf(iproc,trim(seedname)//'_'//num//'.etsf',orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
                   input%hx,input%hy,input%hz,atoms,rxyz,Glr%wfd,wann)
           end if 
         end if
      end if
   end do  ! closing loop on iwann

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Constructing WFs : DONE      !'
      write(*,*) '!==================================!'
   end if


   call deallocate_work_arrays_sumrho(w)
   call deallocate_lr(Glr,subname)
   call deallocate_orbs(orbs,subname)
   call deallocate_orbs(orbsv,subname)
   call deallocate_orbs(orbsw,subname)
   call deallocate_comms(commsw,subname)
   call deallocate_atoms(atoms,subname)
   call free_input_variables(input)


   i_all = -product(shape(wann))*kind(wann)
   deallocate(wann,stat=i_stat)
   call memocc(i_stat,i_all,'wann',subname)
   i_all = -product(shape(wannr))*kind(wannr)
   deallocate(wannr,stat=i_stat)
   call memocc(i_stat,i_all,'wannr',subname)
   i_all = -product(shape(umn))*kind(umn)
   deallocate(umn,stat=i_stat)
   call memocc(i_stat,i_all,'umn',subname)
   i_all = -product(shape(psi))*kind(psi)
   deallocate(psi,stat=i_stat)
   call memocc(i_stat,i_all,'psi',subname)
   i_all = -product(shape(rxyz))*kind(rxyz)
   deallocate(rxyz,stat=i_stat)
   call memocc(i_stat,i_all,'rxyz',subname)
   i_all = -product(shape(wann_list))*kind(wann_list)
   deallocate(wann_list,stat=i_stat)
   call memocc(i_stat,i_all,'wann_list',subname)
   i_all = -product(shape(ncenters))*kind(ncenters)
   deallocate(ncenters,stat=i_stat)
   call memocc(i_stat,i_all,'ncenters',subname)


   !#########################################################
   ! Ending timing and MPI
   !#########################################################
   call timing(iproc,'             ','RE')

   call cpu_time(tcpu1)
   call system_clock(ncount1,ncount_rate,ncount_max)
   tel=dble(ncount1-ncount0)/dble(ncount_rate)
   if (iproc == 0) &
     write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0
    !finalize memory counting
    call memocc(0,0,'count','stop')
 
   ! Barrier suggested by support for titane.ccc.cea.fr, before finalise.
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
   call MPI_FINALIZE(ierr)

end program Waco

subroutine read_inter_list(iproc,n_virt, virt_list)
   
   ! This routine reads the list of virtual orbitals needed

   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt,iproc
   integer, dimension(n_virt), intent(out) :: virt_list

   ! Local variables
   integer :: i,j


   OPEN(11, FILE='input.inter', STATUS='OLD')

!   write(*,*) '!==================================!'
!   write(*,*) '!  Reading virtual orbitals list : !'
!   write(*,*) '!==================================!'

   do i=1,6
      read(11,*) ! Skip first lines
   end do
   read(11,*) (virt_list(j), j=1,n_virt)
   CLOSE(11)

   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!Reading virtual orbitals list done!'
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_inter_list

subroutine read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)

   ! This routine reads the first lines of a .inter file

   implicit none

   ! I/O variables
   integer, intent(in) :: iproc
   character, intent(out) :: seedname*16, filetype*4
   integer, intent(out) :: n_occ, n_virt, n_virt_tot
   logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check

   ! Local variables
   character :: char1*1, char2*1, char3*1, char4*1
   logical :: file_exist
   integer :: ierr

   ! Should check if it exists, if not, make a nice output message
   inquire(file="input.inter",exist=file_exist)
   if (.not. file_exist) then
      if(iproc == 0) then
         write(*,'(A)') 'ERROR : Input file, input.inter, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   OPEN(11, FILE='input.inter', STATUS='OLD')

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Reading input.inter header :   !'
      write(*,*) '!==================================!'
   end if

   w_unk=.false.
   w_sph=.false.
   w_ang=.false.
   w_rad=.false.

   ! First line
   read(11,*) seedname
   if(iproc==0)write(*,*) 'System studied : ', trim(seedname)

   ! Second line
   read(11,*) filetype
   if(iproc==0)write(*,*) 'file type : ', filetype

   ! Third line
   read(11,*) n_occ
   if(iproc==0)write(*,'(A30,I4)') 'Number of occupied orbitals :', n_occ

   ! Fourth line
   read(11,*) char1, n_virt_tot, n_virt
   if (char1=='T') then
      pre_check=.true.
      if(iproc==0)write(*,*) 'Pre-check before calculating Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A38,I4)') 'Total number of unnocupied orbitals :', n_virt_tot
   else
      pre_check=.false.
      if(iproc==0)write(*,*) 'Calculation of Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A39,I4)') 'Number of chosen unnocupied orbitals :', n_virt
   end if

   ! Fifth line
   read(11,*) char1, char2, char3, char4
   if (char1=='T') then
      w_unk=.true.
      if(iproc==0) write(*,*) 'You want to write a UNKp.s file'
   else if (char1 /= 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_unk'
      STOP
   end if
   if (char2=='T') then
      w_sph=.true.
      if(iproc==0) write(*,*) 'You want to write .cube files for spherical harmonics'
   else if (char2 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_sph'
      STOP
   end if
   if (char3=='T') then
      w_ang=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for angular parts of the spherical harmonics'
   else if (char3 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_ang'
      STOP
   end if
   if (char4=='T') then
      w_rad=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for radial parts of the spherical harmonics'
   else if (char4 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_rad'
      STOP
   end if

   CLOSE(11)

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '! Reading input.inter header done  !'
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_inter_header


subroutine scalar_kmeans_diffIG(nIG,crit,nel,vect,string,nbuf)
  use BigDFT_API
  use module_interfaces
  implicit none
  integer, intent(in) :: nel,nIG
  real(kind=8),intent(in) :: crit
  real(kind=8), dimension(nel),intent(in) :: vect
  character(len=*),intent(in) :: string
  integer, intent(out) :: nbuf
  ! local variables
  character(len=*), parameter :: subname='scalar_kmeans_diffIG'
  integer :: i_stat, i_all, i, iel, iiel, iter
  real :: r
  real(kind=8) :: minold, maxold
  logical :: same
  integer, allocatable :: degen(:), buf(:), oldbuf(:)
  real(kind=8), allocatable :: buffers(:), means(:), diff(:)
  

  ! Implement k-means for the diagonal 
  !1) Initial guess
  if(nIG == 0) then ! Maximum Difference between states (determines the number of clusters)
     allocate(buffers(nel),stat=i_stat)
     call memocc(i_stat,buffers,'buffers',subname)
     nbuf = 0
     buffers = 9999.99999
     do iel = 1, nel
        same = .false.
        do iiel = 1, nel
           if(abs(vect(iel) - buffers(iiel)) < crit) then
             same = .true.
           end if
        end do
        if(.not. same) then
           nbuf = nbuf +1
           buffers(nbuf) = vect(iel)
        end if
     end do
  else if (nIG > 0) then  ! Random IG
     nbuf = nIG
     allocate(buffers(nbuf),stat=i_stat)
     call memocc(i_stat,buffers,'buffers',subname)
     do i= 1, nbuf
        call init_random_seed(i)
        call random_number(r)          !returns a random number between [0,1]
        iel = int(r*nel)
        buffers(i) = vect(iel)
     end do
  end if

  ! Store the initial guess in means
  allocate(means(nbuf),stat=i_stat)
  call memocc(i_stat,means,'means',subname)
  do i = 1, nbuf
     means(i) = buffers(i)
  end do

  allocate(diff(nbuf),stat=i_stat)
  call memocc(i_stat,diff,'diff',subname)
  allocate(degen(nbuf),stat=i_stat)
  call memocc(i_stat,degen,'degen',subname)
  allocate(buf(nel),stat=i_stat)
  call memocc(i_stat,buf,'buf',subname)
  allocate(oldbuf(nel),stat=i_stat)
  call memocc(i_stat,oldbuf,'oldbuf',subname)
 
  buf = 0 
  iter = 0
  loop_iter: do
     iter = iter + 1
     ! Copy the buffers
     do iel = 1, nel
        oldbuf(iel) = buf(iel)
     end do

     !2) Associate the numbers
     do iel = 1, nel
        do i = 1, nbuf
           diff(i) = abs(vect(iel) - means(i))
        end do
        buf(iel) =  minloc(diff,dim=1)
     end do
   
     !4) Stop iterating if assignement does not change (ordering can change?)
     same = .true.
     do iel = 1, nel
        if (buf(iel) .ne. oldbuf(iel)) then
           same = .false.
           exit
        end if
     end do
     if(same) exit loop_iter
 
     !3) Calculate new means
     do i = 1, nbuf
        means(i) = 0.0
        degen(i) = 0
        do iel = 1, nel
           if(buf(iel) .ne. i) cycle
           means(i) = means(i) + vect(iel)
           degen(i) = degen(i) + 1
        end do
          means(i) = means(i)/max(degen(i),1)
     end do
  end do loop_iter

  write(*,'(A,x,i4.4,x,A)') 'Convergence reached in',iter,'iterations.'
  write(*,'(A,A,A,1x,i4.4,1x,A)') 'The ',trim(string),' can be clustered in',nbuf,'elements:'
  do i = 1, nbuf
     minold = huge(minold)
     maxold = -huge(maxold)
     do iel = 1, nel
        if(buf(iel) .ne. i) cycle
        if(vect(iel) < minold) minold = vect(iel)
        if(vect(iel) > maxold) maxold = vect(iel)
     end do
     write(*,'(E14.6,3x,A,i4,2x,A,2x,E14.6,2x,A,2x,E14.6)')   means(i), 'with',degen(i),'elements from:',minold,'to',maxold
  end do

  i_all = -product(shape(means))*kind(means)
  deallocate(means,stat=i_stat)
  call memocc(i_stat,i_all,'means',subname)
  i_all = -product(shape(buffers))*kind(buffers) 
  deallocate(buffers,stat=i_stat)
  call memocc(i_stat,i_all,'buffers',subname)
  i_all = -product(shape(diff))*kind(diff) 
  deallocate(diff,stat=i_stat)
  call memocc(i_stat,i_all,'diff',subname)
  i_all = -product(shape(degen))*kind(degen)
  deallocate(degen,stat=i_stat)
  call memocc(i_stat,i_all,'degen',subname)
  i_all = -product(shape(buf))*kind(buf)
  deallocate(buf,stat=i_stat)
  call memocc(i_stat,i_all,'buf',subname)
  i_all = -product(shape(oldbuf))*kind(oldbuf)
  deallocate(oldbuf,stat=i_stat)
  call memocc(i_stat,i_all,'oldbuf',subname)

end subroutine scalar_kmeans_diffIG

subroutine init_random_seed(shuffler)
  integer, intent(in) :: shuffler
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
          
  call random_seed(size = n)
  allocate(seed(n))
  
  call system_clock(count=clock)
  
  seed = clock*shuffler + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  
  deallocate(seed)
END SUBROUTINE init_random_seed

