!> @file
!!  Routines to do frequencies calculation by finite difference
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! @todo
!!  Add higher order for finite difference
!!  Maybe possibility to use Lanczos to determine lowest frequencies
!!  Indicate correct formulae for entropy


!> Calculate vibrational frequencies by frozen phonon approximation.
!! Use a file 'frequencies.res' to restart calculations.
program frequencies

   use module_base
   use module_types
   use module_interfaces
   use m_ab6_symmetry
   use yaml_output
   implicit none

   !Parameters
   character(len=*), parameter :: subname='frequencies'
   !if |freq1-freq2|<tol_freq, freq1 and freq2 are identical
   real(gp), parameter :: tol_freq=1.d-11
   !Temperature (300K)
   real(gp), parameter :: Temperature=300.0_gp

   character(len=4) :: cc
   !File unit
   integer, parameter :: u_hessian=20
   integer :: iproc,nproc,iat,jat,i,j,i_stat,i_all,ierr,infocode,ity,nconfig
   real(gp) :: etot,alat,dd,rmass,fnoise
   character(len=60) :: run_id
   !Input variables
   type(atoms_data) :: atoms
   type(input_variables) :: inputs
   type(restart_objects) :: rst
   !Atomic coordinates, forces
   real(gp), dimension(:), allocatable :: fxyz
   real(gp), dimension(:,:), allocatable :: rpos
   real(gp), dimension(:,:), allocatable :: fpos
   real(gp), dimension(:,:), pointer :: rxyz
   ! hessian, eigenvectors
   real(gp), dimension(:,:), allocatable :: hessian,vector_l,vector_r
   real(gp), dimension(:), allocatable :: eigen_r,eigen_i,sort_work
   !Array to sort eigenvalues
   integer, dimension(:), allocatable :: iperm
   !Array which indicates moves to calculate for a given direction
   integer, dimension(:), allocatable :: kmoves
   ! logical: .true. if already calculated
   logical, dimension(:,:), allocatable :: moves
   real(gp), dimension(:,:), allocatable :: energies
   real(gp), dimension(:,:,:), allocatable :: forces
   real(gp), dimension(3) :: freq_step
   real(gp), dimension(6) :: strten
   real(gp) :: zpenergy,freq_exp,freq2_exp,vibrational_entropy,vibrational_energy,total_energy
   real(gp) :: tel
   real :: tcpu0,tcpu1
   integer :: k,km,ii,jj,ik,imoves,order,n_order,ncount0,ncount1,ncount_rate,ncount_max
   logical :: exists
   integer, dimension(4) :: mpi_info

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(mpi_info,nconfig,run_id,ierr)

   !just for backward compatibility
   iproc=mpi_info(1)
   nproc=mpi_info(2)

   if (nconfig < 0) stop 'runs-file not supported for frequencies executable'

   !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
   ! Read all input files. This should be the sole routine which is called to initialize the run.
   call bigdft_set_input(trim(run_id),'posinp',rxyz,inputs,atoms)

   ! Read all input files.
   inquire(file="input.freq",exist=exists)
   if (.not. exists) then
      if (bigdft_mpi%iproc == 0) write(*,*)'ERROR: need file input.freq for vibrational frequencies calculations.'
      if(nproc/=0)   call MPI_FINALIZE(ierr)
      stop
   end if
   call frequencies_input_variables_new(bigdft_mpi%iproc,.true.,'input.freq',inputs)

   !Order of the finite difference scheme
   order = inputs%freq_order
   if (order == -1) then
      n_order = 1
      allocate(kmoves(n_order),stat=i_stat)
      kmoves = (/ -1 /)
   else if (order == 1) then
      n_order = 1
      allocate(kmoves(n_order),stat=i_stat)
      kmoves = (/ 1 /)
   else if (order == 2) then
      n_order = 2
      allocate(kmoves(n_order),stat=i_stat)
      kmoves = (/ -1, 1 /)
   else if (order == 3) then
      n_order = 4
      allocate(kmoves(n_order),stat=i_stat)
      kmoves = (/ -2, -1, 1, 2 /)
   else
      print *, "Frequencies: This order",order," is not implemented!"
      stop
   end if
   call memocc(i_stat,kmoves,'kmoves',subname)

   ! Allocations
   allocate(fxyz(3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,fxyz,'fxyz',subname)
   allocate(moves(n_order,0:3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,moves,'moves',subname)
   allocate(energies(n_order,0:3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,energies,'energies',subname)
   allocate(forces(3*atoms%nat,n_order,0:3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,forces,'forces',subname)
   allocate(rpos(3,atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,rpos,'rpos',subname)
   allocate(fpos(3*atoms%nat,n_order+ndebug),stat=i_stat)
   call memocc(i_stat,fpos,'fpos',subname)
   allocate(hessian(3*atoms%nat,3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,hessian,'hessian',subname)

   ! Initialize the Hessian
   hessian = 0.d0
   ! Initialize freq_step (step to move atoms)
   freq_step(1) = inputs%freq_alpha*inputs%hx
   freq_step(2) = inputs%freq_alpha*inputs%hy
   freq_step(3) = inputs%freq_alpha*inputs%hz

   call init_restart_objects(bigdft_mpi%iproc,inputs,atoms,rst,subname)

   !Initialize the moves using a restart file if present
   call frequencies_read_restart(atoms%nat,n_order,imoves,moves,energies,forces,freq_step,atoms%amu,etot)
   !Message
   if (bigdft_mpi%iproc == 0) then
      write(*,'(1x,a,i0,a,i0,a)') '=F=> There are ', imoves, ' moves already calculated over ', &
         &   n_order*3*atoms%nat,' frequencies.'
      write(*,*)
   end if

   !Reference state
   if (moves(1,0)) then
      fxyz = forces(:,1,0)
      infocode=0
   else
      call call_bigdft(nproc,bigdft_mpi%iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)
      call frequencies_write_restart(bigdft_mpi%iproc,0,0,0,rxyz,etot,fxyz,&
         &   n_order=n_order,freq_step=freq_step,amu=atoms%amu)
      moves(:,0) = .true.
      call restart_inputs(inputs)
   end if

   if (bigdft_mpi%iproc == 0) then
      write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
      !Print atomic forces
      !call write_forces(atoms,fxyz)
      !This file contains the hessian for post-processing: it is regenerated each time.
      open(unit=u_hessian,file='hessian.dat',status="unknown")
      write(u_hessian,'(a,3(1pe20.10))') '#step=',freq_step(:)
      write(u_hessian,'(a,100(1pe20.10))') '#--',etot,fxyz
   end if

   if (bigdft_mpi%iproc == 0) then
      write(*,*)
      write(*,'(1x,a,59("="))') '=Frequencies calculation '
   end if

   do iat=1,atoms%nat

      if (atoms%ifrztyp(iat) == 1) then
         if (bigdft_mpi%iproc == 0) write(*,"(1x,a,i0,a)") '=F:The atom ',iat,' is frozen.'
         cycle
      end if

      do i=1,3
         ii = i+3*(iat-1)
         if (i==1) then
            alat=atoms%alat1
            cc(3:4)='*x'
         else if (i==2) then
            alat=atoms%alat2
            cc(3:4)='*y'
         else
            alat=atoms%alat3
            cc(3:4)='*z'
         end if
         km = 0
         do ik=1,n_order
            k = kmoves(ik)
            !-1-> 1, 1 -> 2, y = ( x + 3 ) / 2
            km = km + 1
            if (moves(km,ii)) then
               !This move is already done. We use the values from the restart file.
               fpos(:,km) = forces(:,km,ii)
               cycle
            end if
            write(cc(1:2),"(i2)") k
            !Displacement
            dd=real(k,gp)*freq_step(i)
            !We copy atomic positions
            rpos=rxyz
            if (bigdft_mpi%iproc == 0) then
               write(*,"(1x,a,i0,a,a,a,1pe20.10,a)") &
                  &   '=F Move the atom ',iat,' in the direction ',cc,' by ',dd,' bohr'
            end if
            if (atoms%geocode == 'P') then
               rpos(i,iat)=modulo(rxyz(i,iat)+dd,alat)
            else if (atoms%geocode == 'S') then
               rpos(i,iat)=modulo(rxyz(i,iat)+dd,alat)
            else
               rpos(i,iat)=rxyz(i,iat)+dd
            end if
            call call_bigdft(nproc,bigdft_mpi%iproc,atoms,rpos,inputs,etot,fpos(:,km),strten,fnoise,&
                 rst,infocode)
            call frequencies_write_restart(bigdft_mpi%iproc,km,i,iat,rpos,etot,fpos(:,km))
            moves(km,ii) = .true.
            call restart_inputs(inputs)
         end do
         ! Build the Hessian
         do jat=1,atoms%nat
            rmass = amu_emass*sqrt(atoms%amu(atoms%iatype(iat))*atoms%amu(atoms%iatype(jat)))
            do j=1,3
               jj = j+3*(jat-1)
               !Force is -dE/dR
               if (order == -1) then
                  dd = - (fxyz(jj) - fpos(jj,1))/freq_step(i)
               else if (order == 1) then
                  dd = - (fpos(jj,1) - fxyz(jj))/freq_step(i)
               else if (order == 2) then
                  dd = - (fpos(jj,2) - fpos(jj,1))/(2.d0*freq_step(i))
               else if (order == 3) then
                  dd = - (fpos(jj,4) + fpos(jj,3) - fpos(jj,2) - fpos(jj,1))/(6.d0*freq_step(i))
               else
                  stop "BUG: frequencies this order is not defined"
               end if
               !if (abs(dd).gt.1.d-10) then
               hessian(jj,ii) = dd/rmass
               !end if
            end do
         end do
         if (bigdft_mpi%iproc == 0) write(u_hessian,'(i0,1x,i0,1x,100(1pe20.10))') i,iat,hessian(:,ii)
      end do
   end do

   close(unit=u_hessian)

   !Deallocations
   i_all=-product(shape(rpos))*kind(rpos)
   deallocate(rpos,stat=i_stat)
   call memocc(i_stat,i_all,'rpos',subname)
   i_all=-product(shape(fpos))*kind(fpos)
   deallocate(fpos,stat=i_stat)
   call memocc(i_stat,i_all,'fpos',subname)
   i_all=-product(shape(kmoves))*kind(kmoves)
   deallocate(kmoves,stat=i_stat)
   call memocc(i_stat,i_all,'kmoves',subname)

   !allocations
   allocate(eigen_r(3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,eigen_r,'eigen_r',subname)
   allocate(eigen_i(3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,eigen_i,'eigen_i',subname)
   allocate(vector_r(3*atoms%nat,3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,vector_r,'vector_r',subname)
   allocate(vector_l(3*atoms%nat,3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,vector_l,'vector_l',subname)
   allocate(sort_work(3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,sort_work,'sort_work',subname)
   allocate(iperm(3*atoms%nat+ndebug),stat=i_stat)
   call memocc(i_stat,iperm,'iperm',subname)

   !Start timing only for the last part
   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max)

   !Diagonalise the hessian matrix
   call solve(hessian,3*atoms%nat,eigen_r,eigen_i,vector_l,vector_r)
   !Sort eigenvalues in ascending order (use abinit routine sort_dp)
   sort_work=eigen_r
   do i=1,3*atoms%nat
      iperm(i)=i
   end do
   call sort_dp(3*atoms%nat,sort_work,iperm,tol_freq)

   if (bigdft_mpi%iproc == 0) then
      write(*,*)
      write(*,'(1x,a,81("="))') '=F '
      write(*,*)
      write(*,'(1x,a,1x,100(1pe20.10))') '=F: eigenvalues (real)      =',eigen_r(iperm(3*atoms%nat:1:-1))
      write(*,'(1x,a,1x,100(1pe20.10))') '=F: eigenvalues (imaginary) =',eigen_i(iperm(3*atoms%nat:1:-1))
      do i=1,3*atoms%nat
         if (eigen_r(i)<0.0_dp) then
            eigen_r(i)=-sqrt(-eigen_r(i))
         else
            eigen_r(i)= sqrt( eigen_r(i))
         end if
      end do
      write(*,'(1x,a,1x,100(1pe20.10))') '=F: frequencies (Hartree)   =',eigen_r(iperm(3*atoms%nat:1:-1))
      write(*,'(1x,a,1x,100(f13.2))')    '=F: frequencies (cm-1)      =',eigen_r(iperm(3*atoms%nat:1:-1))*Ha_cmm1
      !Build frequencies.xyz in descending order
      open(unit=15,file='frequencies.xyz',status="unknown")
      do i=3*atoms%nat,1,-1
         write(15,'(1x,i0,1x,1pe20.10,a)') atoms%nat,eigen_r(iperm(i))
         write(15,'(1x,a)') 'Frequency'
         do iat=1,atoms%nat
            ity=atoms%iatype(iat)
            do j=1,3
               write(15,'(1x,a,1x,100(1pe20.10))') &
                  &   atoms%atomnames(ity),vector_l(3*(iat-1)+j,iperm(i))
            end do
         end do
         !Blank line
         write(15,*)
      end do
      close(unit=15)
      !Vibrational entropy of the molecule
      ! See : http://www.codessa-pro.com/descriptors/thermodynamic/entropy.htm)
      !       http://www.ncsu.edu/chemistry/franzen/public_html/CH795N/lecture/XIV/XIV.html
      !Zero-point energy
      zpenergy = 0.0_gp
      vibrational_energy=0.0_gp
      vibrational_entropy=0.0_gp
      !iperm: ascending order
      !Remove almost zero frequencies
      do i=6,3*atoms%nat
         freq_exp=exp(eigen_r(iperm(i))*Ha_K/Temperature)
         freq2_exp=exp(-eigen_r(iperm(i))*Ha_K/(2.0_gp*Temperature))
         zpenergy=zpenergy+0.5_gp*eigen_r(iperm(i))
         vibrational_energy=vibrational_entropy+eigen_r(iperm(i))*(0.5_gp+1.0_gp/(freq_exp-1.0_gp))
         vibrational_entropy=vibrational_entropy + eigen_r(iperm(i))*freq2_exp/(1.0_gp-freq2_exp) - log(1.0_gp-freq2_exp)
      end do
      !Multiply by 1/kT
      vibrational_entropy=vibrational_entropy*Ha_K/Temperature
      total_energy=energies(1,0)+vibrational_energy
      write(*,'(1x,a,81("="))') '=F '
      write(*,'(1x,a,f13.2,1x,a,5x,1pe20.10,1x,a)') &
         &   '=F: Zero-point energy   =', zpenergy*Ha_cmm1, 'cm-1',zpenergy,'Hartree'
      write(*,'(1x,a,1pe22.10,a,0pf5.1,a)') &
         &   '=F: Vibrational entropy =', vibrational_entropy,' at ',Temperature,'K'
      write(*,'(1x,a,f13.2,1x,a,5x,1pe20.10,1x,a,0pf5.1,a)') &
         &   '=F: Vibrational  energy =', vibrational_energy*Ha_cmm1, 'cm-1',vibrational_energy,'Hartree at ',Temperature,'K'
      write(*,'(1x,a,1pe22.10,1x,a,0pf5.1,a)') &
         &   '=F: Total energy        =', total_energy,'Hartree at ',Temperature,'K'
   end if

!!$   !Deallocations
!!$   call deallocate_atoms(atoms,subname)
!!$
!!$
!!$!   call deallocate_local_zone_descriptors(rst%Lzd, subname)
!!$   call free_restart_objects(rst,subname)

   i_all=-product(shape(rxyz))*kind(rxyz)
   deallocate(rxyz,stat=i_stat)
   call memocc(i_stat,i_all,'rxyz',subname)
   i_all=-product(shape(fxyz))*kind(fxyz)
   deallocate(fxyz,stat=i_stat)
   call memocc(i_stat,i_all,'fxyz',subname)

   i_all=-product(shape(hessian))*kind(hessian)
   deallocate(hessian,stat=i_stat)
   call memocc(i_stat,i_all,'hessian',subname)
   i_all=-product(shape(eigen_r))*kind(eigen_r)
   deallocate(eigen_r,stat=i_stat)
   call memocc(i_stat,i_all,'eigen_r',subname)
   i_all=-product(shape(eigen_i))*kind(eigen_i)
   deallocate(eigen_i,stat=i_stat)
   call memocc(i_stat,i_all,'eigen_i',subname)
   i_all=-product(shape(vector_l))*kind(vector_l)
   deallocate(vector_l,stat=i_stat)
   call memocc(i_stat,i_all,'vector_l',subname)
   i_all=-product(shape(vector_r))*kind(vector_r)
   deallocate(vector_r,stat=i_stat)
   call memocc(i_stat,i_all,'vector_r',subname)

   i_all=-product(shape(iperm))*kind(iperm)
   deallocate(iperm,stat=i_stat)
   call memocc(i_stat,i_all,'iperm',subname)
   i_all=-product(shape(sort_work))*kind(sort_work)
   deallocate(sort_work,stat=i_stat)
   call memocc(i_stat,i_all,'sort_work',subname)

   i_all=-product(shape(moves))*kind(moves)
   deallocate(moves,stat=i_stat)
   call memocc(i_stat,i_all,'moves',subname)
   i_all=-product(shape(energies))*kind(energies)
   deallocate(energies,stat=i_stat)
   call memocc(i_stat,i_all,'energies',subname)
   i_all=-product(shape(forces))*kind(forces)
   deallocate(forces,stat=i_stat)
   call memocc(i_stat,i_all,'forces',subname)

   call free_restart_objects(rst,subname)

   call deallocate_atoms(atoms,subname) 

   call bigdft_free_input(inputs)


!!$   call free_input_variables(inputs)

   !Final timing
   call cpu_time(tcpu1)
   call system_clock(ncount1,ncount_rate,ncount_max)
   tel=real(ncount1-ncount0,kind=gp)/real(ncount_rate,kind=gp)
   if (bigdft_mpi%iproc == 0) &
      &   write( *,'(1x,a,1x,i4,2(1x,f12.2))') '=F: CPU time/ELAPSED time for root process ', bigdft_mpi%iproc,tel,tcpu1-tcpu0

!!$   !Finalize memory counting
!!$   call memocc(0,0,'count','stop')
!!$
!!$   call MPI_FINALIZE(ierr)

   call bigdft_finalize(ierr)

   contains

   subroutine solve(hessian,n,eigen_r,eigen_i,vector_l,vector_r)
      implicit none
      integer, intent(in) :: n
      real(gp), intent(inout) :: hessian(n,n)
      real(gp), intent(out) :: eigen_r(n),eigen_i(n),vector_l(n,n),vector_r(n,n)
      !Local variables
      character(len=*), parameter :: subname = "solve"
      integer :: info,lwork
      real(gp), dimension(:), allocatable :: work

      lwork=6*n
      allocate(work(lwork+ndebug),stat=i_stat)
      call memocc(i_stat,work,'work',subname)

      call dgeev('V','V',n,hessian,n,eigen_r,eigen_i,vector_l,n,vector_r,n,work,lwork,info)

      if (info /= 0) then
         write(*,'(1x,a,i0)') 'Error from the routine dgeev: info=',info
      end if

      !De-allocation
      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work',subname)

   END SUBROUTINE solve


   subroutine frequencies_read_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,etot)
      implicit none
      !Arguments
      integer, intent(in) :: nat,n_order
      logical, dimension(n_order,0:3*nat), intent(out) :: moves
      real(gp), dimension(n_order,0:3*nat), intent(out) :: energies
      real(gp), dimension(3*nat,n_order,0:1+3*nat), intent(out) :: forces
      real(gp), intent(in) :: freq_step(3)
      integer, intent(out) :: imoves
      real(gp), intent(out) :: etot
      real(gp), dimension(:), intent(out) :: amu
      !Local variables
      character(len=*), parameter :: subname = "frequencies_read_restart"
      logical :: exists
      integer, parameter :: iunit = 15
      character(len=*), parameter :: freq_form = '(1x,"=F=> ",a)'
      real(gp) :: steps(3)
      real(gp), dimension(:), allocatable :: rxyz,fxyz
      integer :: ierror,km,i,iat,ii,i_order
      !Initialize by default to false
      imoves=0
      moves = .false.
      !Test if the file does exist.
      if (bigdft_mpi%iproc == 0) then
         write(*,freq_form) 'Check if the file "frequencies.res" is present.'
      end if
      inquire(file='frequencies.res', exist=exists)
      if (.not.exists) then
         !There is no restart file.
         call razero(n_order*(3*nat+1),energies(1,0))
         if (bigdft_mpi%iproc == 0) write(*,freq_form) 'No "frequencies.res" file present.'
         return
      end if

      !Allocations
      allocate(rxyz(3*nat))
      call memocc(i_stat,rxyz,'rxyz',subname)
      allocate(fxyz(3*nat))
      call memocc(i_stat,fxyz,'fxyz',subname)

      !We read the file
      open(unit=iunit,file='frequencies.res',status='old',form='unformatted')
      !First line is data for coherency of the calculation
      read(unit=iunit,iostat=ierror) i_order,steps,amu
      if (ierror /= 0) then
         !Read error, we clean the file.
         if (bigdft_mpi%iproc == 0) then
            close(unit=iunit)
            write(*,freq_form) 'Erase the file "frequencies.res".'
            open(unit=iunit,file='frequencies.res',status='old',form='unformatted')
            write(iunit)
            close(unit=iunit)
         end if
         return
      else
         if (steps(1) /= freq_step(1) .or. &
            &   steps(2) /= freq_step(2) .or. &
         steps(3) /= freq_step(3)) then
         if (bigdft_mpi%iproc == 0) write(*,freq_form) 'The step to calculate frequencies is not the same: stop.'
         stop
      end if
      if (i_order > n_order) then
         if (bigdft_mpi%iproc == 0) then 
            write(*,freq_form) 'The number of points per direction is bigger in the "frequencies.res" file.'
            write(*,freq_form) 'Increase the order of the finite difference scheme'
         end if
         stop
      end if
   end if
   !Read the reference state
   read(unit=iunit,iostat=ierror) iat,etot,rxyz,fxyz
   if (ierror /= 0 .or. iat /= 0) then
      !Read error, we assume that it is not calculated
      if (bigdft_mpi%iproc == 0) write(*,freq_form) 'The reference state is not calculated in "frequencies.res" file.'
   else
      if (bigdft_mpi%iproc == 0) write(*,freq_form) 'The reference state is already calculated.'
      energies(:,0) = etot
      forces(:,1,0) = fxyz
      moves(:,0) = .true.
   end if
   do
      read(unit=iunit,iostat=ierror) km,i,iat,rxyz,etot,fxyz
      if (ierror /= 0) then
         !Read error, we exit
         exit
      end if
      ii = i + 3*(iat-1)
      imoves = imoves + 1
      energies(km,ii) = etot
      forces(:,km,ii) = fxyz
      moves(km,ii) = .true.
   end do
   close(unit=iunit)

   !Deallocations
   i_all=-product(shape(rxyz))*kind(rxyz)
   deallocate(rxyz)
   call memocc(i_stat,i_all,'rxyz',subname)
   i_all=-product(shape(fxyz))*kind(fxyz)
   deallocate(fxyz)
   call memocc(i_stat,i_all,'fxyz',subname)

END SUBROUTINE frequencies_read_restart


subroutine frequencies_write_restart(iproc,km,i,iat,rxyz,etot,fxyz,n_order,freq_step,amu)
   implicit none
   !Arguments
   integer, intent(in) :: iproc,km,i,iat
   real(gp), dimension(:,:), intent(in) :: rxyz
   real(gp), intent(in) :: etot
   real(gp), dimension(:), intent(in) :: fxyz
   integer, intent(in), optional :: n_order
   real(gp), intent(in), optional :: freq_step(3)
   real(gp), dimension(:), intent(in), optional :: amu
   !Local variables
   integer, parameter :: iunit = 15

   if (km == 0 .and. &
      &   .not.(present(n_order).and.present(freq_step).and.present(amu))) then
   if (iproc == 0) write(*,*) "Bug for use of frequencies_write_restart"
   stop
end if
if (iproc ==0 ) then
   !This file is used as a restart
   open(unit=iunit,file='frequencies.res',status="unknown",form="unformatted",position="append")
   if (km == 0) then
      write(unit=iunit) n_order,freq_step,amu
      write(unit=iunit) 0,etot,rxyz,fxyz
   else
      write(unit=iunit) km,i,iat,etot,rxyz,fxyz
   end if
   close(unit=iunit)
end if
  END SUBROUTINE frequencies_write_restart


  subroutine restart_inputs(inputs)
     implicit none
     !Argument
     type(input_variables), intent(inout) :: inputs
     inputs%inputPsiId=1
  END SUBROUTINE restart_inputs


END PROGRAM frequencies


!> Integrate forces (not used)
subroutine integrate_forces(iproc,n_moves) !n(c) energies,forces (arg:2,3)

   use module_base

   implicit none
   !Arguments
   integer, intent(in) :: iproc,n_moves
   !n(c) real(gp), intent(in) :: energies(n_moves)
   !n(c) real(gp), intent(in) :: forces(3*nat,n_moves)
   !Local variables
   character(len=*), parameter :: subname = "integrate_forces"
   real(gp), dimension(:), allocatable :: weight
   !n(c) real(gp) :: path
   integer :: i,i_stat,i_all

   !Allocation
   allocate(weight(n_moves+ndebug),stat=i_stat)
   call memocc(i_stat,weight,'weight',subname)

   !Prepare the array of the correct weights of the iteration steps
   if (mod(n_moves,2).ne.1) then
      if (iproc == 0) write(*,*) 'the number of iteration steps has to be odd'
      stop
   end if
   weight(1)=1.d0/3.d0
   weight(2)=4.d0/3.d0
   do i=3,n_moves-2,2
      weight(i)=2.d0/3.d0
      weight(i+1)=4.d0/3.d0
   enddo
   weight(n_moves)=1.d0/3.d0

   !Start integration
   !n(c) path = 0_gp
   do i=1,n_moves
   end do

   !De-allocation
   i_all=-product(shape(weight))*kind(weight)
   deallocate(weight,stat=i_stat)
   call memocc(i_stat,i_all,'weight',subname)

END SUBROUTINE integrate_forces
