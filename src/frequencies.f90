!> @file
!!  Routines to do frequencies calculation by finite difference
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
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
   real(gp), parameter :: Temperature=300.0_gp !< Temperature (300K)
   character(len=*), dimension(3), parameter :: cc = (/ 'x', 'y', 'z' /)
   !File unit
   integer, parameter :: u_hessian=20
   real(gp) :: alat,dd,rmass
   character(len=60) :: run_id
   !Input variables
   type(run_objects) :: runObj
   type(DFT_global_output) :: outs
   !Atomic coordinates, forces
   real(gp), dimension(:,:), allocatable :: rxyz0
   real(gp), dimension(:,:), allocatable :: fpos
   real(gp), dimension(:,:), allocatable :: hessian           !< Hessain matrix
   real(gp), dimension(:,:), allocatable :: vector_l,vector_r !< left and right eignevectors
   real(gp), dimension(:), allocatable :: eigen_r,eigen_i     !< Real and Imaginary part of the eigenvalues
   real(gp), dimension(:), allocatable :: sort_work           !< To sort the eigenvalues in ascending order
   integer, dimension(:), allocatable :: iperm   !< Array to sort eigenvalues
   integer, dimension(:), allocatable :: kmoves  !< Array which indicates moves to calculate for a given direction
   logical, dimension(:,:), allocatable :: moves !< logical: .true. if already calculated
   real(gp), dimension(:,:), allocatable :: energies
   real(gp), dimension(:,:,:), allocatable :: forces
   real(gp), dimension(3) :: freq_step
   real(gp) :: zpenergy,freq_exp,freq2_exp,vibrational_entropy,vibrational_energy,total_energy
   integer :: k,km,ii,jj,ik,imoves,order,n_order
   integer :: iproc,nproc,igroup,ngroups,icalc
   integer :: iat,jat,i,j,ierr,infocode,ity,nconfig
   logical :: exists
   integer, dimension(4) :: mpi_info

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(mpi_info,nconfig,run_id,ierr)

   if (nconfig < 0) stop 'runs-file not supported for frequencies executable'

   call f_routine(id=subname)

   !just for backward compatibility
   iproc=mpi_info(1)
   nproc=mpi_info(2)
   igroup=mpi_info(3)
   !number of groups
   ngroups=mpi_info(4)

   !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
   ! Read all input files. This should be the sole routine which is called to initialize the run.
   call run_objects_init_from_files(runObj, trim(run_id), 'posinp')

   ! Read all input files.
   inquire(file="input.freq",exist=exists)
   if (.not. exists) then
      if (bigdft_mpi%iproc == 0) write(*,*)'ERROR: need file input.freq for vibrational frequencies calculations.'
      if(nproc/=0)   call MPI_FINALIZE(ierr)
      stop
   end if
   call frequencies_input_variables_new(bigdft_mpi%iproc,.true.,'input.freq',runObj%inputs)

   !Order of the finite difference scheme
   order = runObj%inputs%freq_order
   if (order == -1) then
      n_order = 1
      kmoves = f_malloc(n_order,id='kmoves')
      kmoves = (/ -1 /)
   else if (order == 1) then
      n_order = 1
      kmoves = f_malloc(n_order,id='kmoves')
      kmoves = (/ 1 /)
   else if (order == 2) then
      n_order = 2
      kmoves = f_malloc(n_order,id='kmoves')
      kmoves = (/ -1, 1 /)
   else if (order == 3) then
      n_order = 4
      kmoves = f_malloc(n_order,id='kmoves')
      kmoves = (/ -2, -1, 1, 2 /)
   else
      print *, "Frequencies: This order",order," is not implemented!"
      stop
   end if

   ! Allocations
   call init_global_output(outs, runObj%atoms%astruct%nat)
   rxyz0 = f_malloc((/ 3, runObj%atoms%astruct%nat /), id = 'rxyz0')
   moves = f_malloc((/ 1.to.n_order, 0.to.3*runObj%atoms%astruct%nat /),id='moves')
   energies = f_malloc((/ 1.to.n_order, 0.to.3*runObj%atoms%astruct%nat /),id='energies')
   forces = f_malloc((/ 1.to.3*runObj%atoms%astruct%nat, 1.to.n_order, 0.to.3*runObj%atoms%astruct%nat /),id='forces')
   fpos = f_malloc((/ 3*runObj%atoms%astruct%nat, n_order /),id='fpos')
   hessian = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='hessian')

   ! Initialize the Hessian
   hessian = 0.d0
   ! Initialize freq_step (step to move atoms)
   freq_step(1) = runObj%inputs%freq_alpha*runObj%inputs%hx
   freq_step(2) = runObj%inputs%freq_alpha*runObj%inputs%hy
   freq_step(3) = runObj%inputs%freq_alpha*runObj%inputs%hz
   ! Reference positions.
   call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, rxyz0(1,1), 1)

   !Determination of the calculation id
   icalc=0

   !Initialize the moves using a restart file if present
   !Regenerate it if trouble and indicate if all calculations are done
   call frequencies_check_restart(runObj%atoms%astruct%nat,n_order,imoves,moves,energies,forces,freq_step,runObj%atoms%amu,infocode)
   !We ignore at this stage the infocode

   !Message
   if (bigdft_mpi%iproc == 0) then
      call yaml_map('(F) Frequency moves already calculated',imoves)
      call yaml_map('(F) Total Frequency moves',n_order*3*runObj%atoms%astruct%nat)
   end if

   !Reference state
   if (moves(1,0)) then
      call vcopy(3*outs%fdim, forces(1,1,0), 1, outs%fxyz(1,1), 1)
      outs%energy = energies(1,0)
      infocode=0
   else
      if (bigdft_mpi%iproc == 0) call yaml_comment('(F) Reference state calculation',hfill='=')
      call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
      call frequencies_write_restart(0,0,0,runObj%atoms%astruct%rxyz,outs%energy,outs%fxyz, &
           & n_order=n_order,freq_step=freq_step,amu=runObj%atoms%amu)
      moves(:,0) = .true.
      call restart_inputs(runObj%inputs)
   end if

   if (bigdft_mpi%iproc == 0) then
      call yaml_map('(F) Exit signal for Wavefunction Optimization Finished',infocode)
      !Print atomic forces
      !call write_forces(atoms,fxyz)
      !This file contains the hessian for post-processing: it is regenerated each time.
      open(unit=u_hessian,file='hessian.dat',status="unknown")
      write(u_hessian,'(a,3(1pe20.10))') '#step=',freq_step(:)
      write(u_hessian,'(a,100(1pe20.10))') '#--',outs%energy,outs%fxyz
   end if

   if (bigdft_mpi%iproc == 0) then
      call yaml_comment('(F) Start Frequencies calculation',hfill='=')
   end if

   do iat=1,runObj%atoms%astruct%nat

      if (runObj%atoms%astruct%ifrztyp(iat) == 1) then
         if (bigdft_mpi%iproc == 0) call yaml_comment('(F) The atom ' // trim(yaml_toa(iat)) // ' is frozen.')
         cycle
      end if

      do i=1,3
         ii = i+3*(iat-1)
         if (i==1) then
            !Along x axis
            alat=runObj%atoms%astruct%cell_dim(1)
         else if (i==2) then
            !Along y axis
            alat=runObj%atoms%astruct%cell_dim(2)
         else
            !Along z axis
            alat=runObj%atoms%astruct%cell_dim(3)
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
            !Displacement
            dd=real(k,gp)*freq_step(i)
            !We copy atomic positions
            call vcopy(3*runObj%atoms%astruct%nat, rxyz0(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
            if (bigdft_mpi%iproc == 0) then
               call yaml_open_map('(F) Move',flow=.true.)
                  call yaml_map('atom',      iat)
                  call yaml_map('direction', k)
                  call yaml_map('axis',      cc(i))
                  call yaml_map('displacement (Bohr)', dd,fmt='(1pe20.10)')
               call yaml_close_map()
            end if
            if (runObj%atoms%astruct%geocode == 'P') then
               runObj%atoms%astruct%rxyz(i,iat)=modulo(rxyz0(i,iat)+dd,alat)
            else if (runObj%atoms%astruct%geocode == 'S') then
               runObj%atoms%astruct%rxyz(i,iat)=modulo(rxyz0(i,iat)+dd,alat)
            else
               runObj%atoms%astruct%rxyz(i,iat)=rxyz0(i,iat)+dd
            end if
            call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
            call frequencies_write_restart(km,i,iat,runObj%atoms%astruct%rxyz,outs%energy,outs%fxyz)
            call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fpos(1,km), 1)
            moves(km,ii) = .true.
            call restart_inputs(runObj%inputs)
         end do
         ! Build the Hessian
         do jat=1,runObj%atoms%astruct%nat
            rmass = amu_emass*sqrt(runObj%atoms%amu(runObj%atoms%astruct%iatype(iat))* &
                 & runObj%atoms%amu(runObj%atoms%astruct%iatype(jat)))
            do j=1,3
               jj = j+3*(jat-1)
               !Force is -dE/dR
               if (order == -1) then
                  dd = - (outs%fxyz(j,jat) - fpos(jj,1))/freq_step(i)
               else if (order == 1) then
                  dd = - (fpos(jj,1) - outs%fxyz(j,jat))/freq_step(i)
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
   call f_free(fpos)
   call f_free(kmoves)

   !allocations
   eigen_r   = f_malloc(3*runObj%atoms%astruct%nat,id='eigen_r')
   eigen_i   = f_malloc(3*runObj%atoms%astruct%nat,id='eigen_i')
   vector_r  = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='vector_r')
   vector_l  = f_malloc((/ 3*runObj%atoms%astruct%nat, 3*runObj%atoms%astruct%nat /),id='vector_l')
   sort_work = f_malloc(3*runObj%atoms%astruct%nat,id='sort_work')
   iperm     = f_malloc(3*runObj%atoms%astruct%nat,id='iperm')

   !Diagonalise the hessian matrix
   call solve(hessian,3*runObj%atoms%astruct%nat,eigen_r,eigen_i,vector_l,vector_r)
   !Sort eigenvalues in ascending order (use abinit routine sort_dp)
   sort_work=eigen_r
   do i=1,3*runObj%atoms%astruct%nat
      iperm(i)=i
   end do
   call sort_dp(3*runObj%atoms%astruct%nat,sort_work,iperm,tol_freq)

   if (bigdft_mpi%iproc == 0) then
      call yaml_comment('(F) Frequencies results',hfill='=')
      call yaml_map('(F) Eigenvalues (real part)',eigen_r(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
      call yaml_map('(F) Eigenvalues (imag part)',eigen_i(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
      do i=1,3*runObj%atoms%astruct%nat
         if (eigen_r(i)<0.0_dp) then
            eigen_r(i)=-sqrt(-eigen_r(i))
         else
            eigen_r(i)= sqrt( eigen_r(i))
         end if
      end do
      call yaml_map('(F) Frequencies (Hartree)',eigen_r(iperm(3*runObj%atoms%astruct%nat:1:-1)),fmt='(1pe20.10)')
      call yaml_map('(F) Frequencies (cm-1)',   eigen_r(iperm(3*runObj%atoms%astruct%nat:1:-1))*Ha_cmm1,fmt='(f13.2)')
      !Build frequencies.xyz in descending order
      open(unit=15,file='frequencies.xyz',status="unknown")
      do i=3*runObj%atoms%astruct%nat,1,-1
         write(15,'(1x,i0,1x,1pe20.10,a)') runObj%atoms%astruct%nat,eigen_r(iperm(i))
         write(15,'(1x,a)') 'Frequency'
         do iat=1,runObj%atoms%astruct%nat
            ity=runObj%atoms%astruct%iatype(iat)
            do j=1,3
               write(15,'(1x,a,1x,100(1pe20.10))') &
                  &   runObj%atoms%astruct%atomnames(ity),vector_l(3*(iat-1)+j,iperm(i))
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
      do i=6,3*runObj%atoms%astruct%nat
         freq_exp=exp(eigen_r(iperm(i))*Ha_K/Temperature)
         freq2_exp=exp(-eigen_r(iperm(i))*Ha_K/(2.0_gp*Temperature))
         zpenergy=zpenergy+0.5_gp*eigen_r(iperm(i))
         vibrational_energy=vibrational_entropy+eigen_r(iperm(i))*(0.5_gp+1.0_gp/(freq_exp-1.0_gp))
         vibrational_entropy=vibrational_entropy + eigen_r(iperm(i))*freq2_exp/(1.0_gp-freq2_exp) - log(1.0_gp-freq2_exp)
      end do
      !Multiply by 1/kT
      vibrational_entropy=vibrational_entropy*Ha_K/Temperature
      total_energy=energies(1,0)+vibrational_energy
      call yaml_map('(F) Zero-point energy (cm-1 and Hartree)', (/ zpenergy*Ha_cmm1, zpenergy /),fmt='(1pe20.10)')
      call yaml_map('(F) Considered Temperature (Kelvin)',      Temperature,fmt='(f5.1)')
      call yaml_map('(F) Vibrational entropy =',                vibrational_entropy,fmt='(1pe22.10)')
      call yaml_map('(F) Vibrational Energy (cm-1 and Hartree)', &
           &                     (/ vibrational_energy*Ha_cmm1, vibrational_energy/),fmt='(1pe20.10)')
      call yaml_map('(F) Total Energy (Hartree)',               total_energy,fmt='(1pe22.10)')
   end if

   ! De-allocations
   call f_free(rxyz0)

   call deallocate_global_output(outs)

   call f_free(hessian)
   call f_free(eigen_r)
   call f_free(eigen_i)
   call f_free(vector_l)
   call f_free(vector_r)

   call f_free(iperm)
   call f_free(sort_work)

   call f_free(moves)
   call f_free(energies)
   call f_free(forces)

   call f_release_routine()

   call run_objects_free(runObj, subname)

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

      call f_routine(id=subname)
      lwork=6*n
      work=f_malloc(lwork+ndebug,id='work')

      call dgeev('V','V',n,hessian,n,eigen_r,eigen_i,vector_l,n,vector_r,n,work,lwork,info)

      if (info /= 0) then
         call yaml_warning('(F) Error from the routine dgeev: info=' // trim(yaml_toa(info)))
      end if

      !De-allocation
      call f_free(work)

   END SUBROUTINE solve


   !> Check the restart file
   subroutine frequencies_check_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      implicit none
      !Arguments
      integer, intent(in) :: nat     !< Number of atoms
      integer, intent(in) :: n_order !< Order of the finite difference
      logical, dimension(n_order,0:3*nat), intent(out) :: moves            !< Contains moves already done
      real(gp), dimension(n_order,0:3*nat), intent(out) :: energies        !< Energies of the already moves
      real(gp), dimension(3*nat,n_order,0:3*nat), intent(out) :: forces    !< Forces of the already moves
      real(gp), dimension(3), intent(in) :: freq_step     !< Frequency step in each direction
      integer, intent(out) :: imoves                      !< Number of frequency already calculated   
      real(gp), dimension(:), intent(out) :: amu          !< Atomic masses
      integer, intent(out) :: ierror                      !< 0 means all calculations are done
      !Local variables
      character(len=*), parameter :: subname = "frequencies_check_restart"
      !We read the file
      call frequencies_read_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      !if (ierror /= 0) then
      !   !If error, we write a new file
      !   call frequencies_write_new_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      !end if
      !Check also if all calculations are done.
   end subroutine frequencies_check_restart


   !> Read the restart file associated to the calculation of the frequencies
   subroutine frequencies_read_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      implicit none
      !Arguments
      integer, intent(in) :: nat     !< Number of atoms
      integer, intent(in) :: n_order !< Order of the finite difference
      logical, dimension(n_order,0:3*nat), intent(out) :: moves            !< Contains moves already done
      real(gp), dimension(n_order,0:3*nat), intent(out) :: energies        !< Energies of the already moves
      real(gp), dimension(3*nat,n_order,0:3*nat), intent(out) :: forces    !< Forces of the already moves
      real(gp), dimension(3), intent(in) :: freq_step     !< Frequency step in each direction
      integer, intent(out) :: imoves                      !< Number of frequency already calculated
      real(gp), dimension(:), intent(out) :: amu          !< Atomic masses
      integer, intent(out) :: ierror                      !< Error when reading the file
      !Local variables
      character(len=*), parameter :: subname = "frequencies_read_restart"
      logical :: exists
      integer, parameter :: iunit = 15
      real(gp), dimension(3) :: steps
      real(gp), dimension(:), allocatable :: rxyz,fxyz
      real(gp) :: etot
      integer :: km,i,iat,ii,i_order

      call f_routine(id=subname)
      !Initialize by default to false
      imoves=0
      moves = .false.

      !Test if the file does exist.
      inquire(file='frequencies.res', exist=exists)
      if (.not.exists) then
         !There is no restart file.
         call razero(n_order*(3*nat+1),energies(1,0))
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies.res" present',.false.)
         !Code error non zero
         ierror = -1
         return
      else
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) File "frequencies.res" present',.true.)
      end if

      !We read the file
      open(unit=iunit,file='frequencies.res',status='old',form='unformatted',iostat=ierror)
      !First line is data for coherency of the calculation
      read(unit=iunit,iostat=ierror) i_order,steps,amu
      if (ierror /= 0) then
         !Read error, we exit
         if (bigdft_mpi%iproc == 0) then
            close(unit=iunit)
            call yaml_warning('(F) Error when reading the first line of "frequencies.res"')
         end if
         call f_release_routine()
         return
      else
         if (steps(1) /= freq_step(1) .or. steps(2) /= freq_step(2) .or. steps(3) /= freq_step(3)) then
            if (bigdft_mpi%iproc == 0) call yaml_warning('(F) The step to calculate frequencies is not the same: stop.')
            stop
         end if

         if (i_order > n_order) then
            if (bigdft_mpi%iproc == 0) then 
               call yaml_warning('(F) The number of points per direction is bigger in the "frequencies.res" file.')
               call yaml_warning('(F) Increase the order of the finite difference scheme')
            end if
            stop
         end if
      end if
      
      !Allocations
      rxyz=f_malloc(3*nat,id='rxyz')
      fxyz=f_malloc(3*nat,id='fxyz')

      !Read the reference state
      read(unit=iunit,iostat=ierror) iat,etot,rxyz,fxyz
      if (ierror /= 0 .or. iat /= 0) then
         !Read error, we assume that it is not calculated
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) Reference state calculated in the "frequencies.res" file',.false.)
      else
         if (bigdft_mpi%iproc == 0) call yaml_map('(F) Reference state calculated in the "frequencies.res" file',.true.)
         energies(:,0) = etot
         forces(:,1,0) = fxyz
         moves(:,0) = .true.
      end if
      do
         read(unit=iunit,iostat=ierror) km,i,iat,rxyz,etot,fxyz
         if (ierror /= 0) then
            !Read error, we exit
            if (bigdft_mpi%iproc == 0) then
               close(unit=iunit)
               !Error if all moves are not read
               if (imoves < 3*nat+1) call yaml_warning('(F) The file "frequencies.res" is incomplete!')
            end if
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
      call f_free(rxyz)
      call f_free(fxyz)

      call f_release_routine()

   END SUBROUTINE frequencies_read_restart


   !> Write the full restart file
   subroutine frequencies_write_new_restart(nat,n_order,imoves,moves,energies,forces,freq_step,amu,ierror)
      implicit none
      !Arguments
      integer, intent(in) :: nat     !< Number of atoms
      integer, intent(in) :: n_order !< Order of the finite difference
      logical, dimension(n_order,0:3*nat), intent(in) :: moves         !< Contains moves already done
      real(gp), dimension(n_order,0:3*nat), intent(in) :: energies     !< Energies of the already moves
      real(gp), dimension(3*nat,n_order,0:3*nat), intent(in) :: forces !< Forces of the already moves
      real(gp), dimension(3), intent(in) :: freq_step    !< Frequency step in each direction
      integer, intent(out) :: imoves                     !< Number of frequency already calculated
      real(gp), dimension(:), intent(in) :: amu          !< Atomic masses
      integer, intent(out) :: ierror                     !< Error when reading the file
      !Local variables
      integer, parameter :: iunit = 15

      if (bigdft_mpi%iproc ==0 ) then
         !This file is used as a restart
         open(unit=iunit,file='frequencies.res',status="unknown",form="unformatted")
         write(unit=iunit) n_order,freq_step,amu

         write(unit=iunit) 0,outs%energy,rxyz0,outs%fxyz
         do iat=1,runObj%atoms%astruct%nat
            if (runObj%atoms%astruct%ifrztyp(iat) == 1) then
               if (bigdft_mpi%iproc == 0) call yaml_comment('(F) The atom ' // trim(yaml_toa(iat)) // ' is frozen.')
               cycle
            end if
            do i=1,3
               ii = i+3*(iat-1)
               km = 0
               do ik=1,n_order
                  km = km + 1
                  if (moves(km,ii)) then
                     write(unit=iunit) km,i,iat,outs%energy,rxyz0,outs%fxyz
                  end if
               end do
            end do
         end do
         close(unit=iunit)
      end if
   END SUBROUTINE frequencies_write_new_restart


   !> Write one move in the file restart (only moves==.true.)
   subroutine frequencies_write_restart(km,i,iat,rxyz,etot,fxyz,n_order,freq_step,amu)
      implicit none
      !Arguments
      integer, intent(in) :: km,i,iat
      real(gp), dimension(:,:), intent(in) :: rxyz
      real(gp), intent(in) :: etot
      real(gp), dimension(:,:), intent(in) :: fxyz
      integer, intent(in), optional :: n_order
      real(gp), intent(in), optional :: freq_step(3)
      real(gp), dimension(:), intent(in), optional :: amu
      !Local variables
      integer, parameter :: iunit = 15

      if (km == 0 .and. .not.(present(n_order).and.present(freq_step).and.present(amu))) then
         if (bigdft_mpi%iproc == 0) write(*,*) "Bug for use of frequencies_write_restart"
         if (bigdft_mpi%iproc == 0) call yaml_warning("(F) Bug for use of frequencies_write_restart")
         stop
      end if

      if (bigdft_mpi%iproc ==0 ) then
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
      integer :: i
   
      !Allocation
      weight = f_malloc(n_moves,id='weight')
   
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
      call f_free(weight)
   
   END SUBROUTINE integrate_forces

END PROGRAM frequencies
