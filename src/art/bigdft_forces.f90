!> @file
!!   Information to interface BigDFT with ART
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by:
!! -EM 2010, see ~/AUTHORS
!! -Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!


!> Module which contains information for BigDFT run inside ART
module bigdft_forces

   use module_base!, only : gp,wp,dp,Bohr_Ang
   use module_types
   use module_atoms
   use module_interfaces
   use defs, only : iproc
   implicit none

   private

   ! Storage of required variables for a SCF loop calculation.
   logical :: initialised = .false.
   logical :: first_time  = .True.
   integer :: nproc, me
   type(run_objects) :: runObj
   type(atoms_data) :: atoms_all !at are the quantum atoms
   real(gp), parameter :: ht2ev = 27.2113834_gp

   real(kind=8) :: gnrm_l
   real(kind=8) :: gnrm_h  ! For lanczos

   public :: bigdft_init_art
   public :: calcforce_bigdft
   public :: mingeo
   public :: bigdft_finalise
   public :: init_all_atoms
   public :: copy_atoms_object
   public :: prepare_quantum_atoms_Si

   ! development : At some key points we do not use the previously 
   ! calculated wave function
   logical, public :: new_wf
   public :: check_force_clean_wf
   integer, dimension(:), allocatable, public :: in_system   ! Constraint over atoms

   contains

   !> ART init_all_atoms
   !! Routine to initialize all the positions. uses BigDFT read files
   subroutine init_all_atoms( nat, typa, posa, const_, boxl, boxtype, nproc_, me_, file )

      implicit none

      !Arguments
      integer,      intent(out)               :: nat
      integer,      pointer                   :: typa(:)
      real(kind=8), pointer                   :: posa(:)
      integer,      pointer                   :: const_(:)
      real(kind=8), dimension(3), intent(out) :: boxl
      character(len=1), intent(out)           :: boxtype
      integer,      intent(in)                :: nproc_
      integer,      intent(in)                :: me_
      character(len=*), intent(in)            :: file

      !Local variables
      integer                                 :: i
      character(len=2)                        :: symbol
      character(len=10)                       :: name
      character(len=*), parameter             :: subname='init_all_atoms'
      !_______________________

      nproc = nproc_
      me = me_

      call set_astruct_from_file(file,me_,atoms_all%astruct)
      call allocate_atoms_data(atoms_all)
      !call allocate_atoms_nat(atoms_all)!, subname)
      !call allocate_atoms_ntypes(atoms_all)!, subname)
      nat = atoms_all%astruct%nat
      boxtype = atoms_all%astruct%geocode

      posa = f_malloc_ptr(3 * nat,id='posa')
      typa = f_malloc_ptr(nat,id='typa')
      const_ = f_malloc_ptr(nat,id='const_')

      do i = 1, nat, 1
         posa(i)           = atoms_all%astruct%rxyz(1, i) * Bohr_Ang
         posa(i + nat)     = atoms_all%astruct%rxyz(2, i) * Bohr_Ang
         posa(i + 2 * nat) = atoms_all%astruct%rxyz(3, i) * Bohr_Ang
         typa(i) = atoms_all%astruct%iatype(i)
      end do

      boxl(1) = atoms_all%astruct%cell_dim(1) * Bohr_Ang
      boxl(2) = atoms_all%astruct%cell_dim(2) * Bohr_Ang
      boxl(3) = atoms_all%astruct%cell_dim(3) * Bohr_Ang
      ! Blocked atoms 
      const_ = 0                          ! Initialization, everyone is free.
      const_(:) = atoms_all%astruct%ifrztyp(:)

      if ( iproc == 0 ) then
         do i=1,atoms_all%astruct%nat
            name=trim(atoms_all%astruct%atomnames(atoms_all%astruct%iatype(i)))
            if (name(3:3)=='_') then
               symbol=name(1:2)
            else if (name(2:2)=='_') then
               symbol=name(1:1)
            else
               symbol=name(1:2)
            end if
            !write(9,'(i3, a2,4x,3(1x,1pe24.17))')i, symbol,(atoms_all%astruct%rxyz(j,i),j=1,3)
         enddo
      end if

   END SUBROUTINE init_all_atoms


   !> ART bigdft_init_art
   !! Routine to initialize all BigDFT stuff
   subroutine bigdft_init_art( nat, me_, nproc_, my_gnrm,passivate,total_nb_atoms )
     use dictionaries
     use module_input_dicts
     use module_atoms, only: read_atomic_file=>set_astruct_from_file
      implicit none

      !Arguments
      integer,      intent(in) :: nat
      integer,      intent(in)  :: me_, nproc_
      real(kind=8), intent(in)  :: my_gnrm
      logical,      intent(in)  :: passivate
      integer,      intent(in)  :: total_nb_atoms

      !Local variables
      character(len=*), parameter :: subname='bigdft_init_art'
      real(gp),dimension(3*total_nb_atoms) :: posquant
      integer :: natoms_calcul, i_stat
      type(dictionary), pointer :: dict
      !_______________________

      call dict_init(dict)
      call read_input_dict_from_files("input", bigdft_mpi,dict)

      me = me_
      nproc = nproc_

      allocate(runObj%atoms)
      runObj%atoms = atoms_data_null()
      allocate(runObj%inputs)
      if (nat .eq. total_nb_atoms .and. .not. passivate) then 
         ! we just reread all atoms
         call read_atomic_file("posinp",me_,runObj%atoms%astruct)
      else 
         !uses the big object to prepare. everything should
         ! be alright in the object exept the length
         call prepare_quantum_atoms_Si(atoms_all,posquant,natoms_calcul)
         !we just copy it in a smaller vect
         call copy_atoms_object(atoms_all,runObj%atoms,runObj%atoms%astruct%rxyz,natoms_calcul,total_nb_atoms,posquant)
         call initialize_atomic_file(me_,runObj%atoms,runObj%atoms%astruct%rxyz)
      endif
      call astruct_merge_to_dict(dict // "posinp", runObj%atoms%astruct, runObj%atoms%astruct%rxyz)

      call atoms_file_merge_to_dict(dict)

      call standard_inputfile_names(runObj%inputs,'input')
      call inputs_from_dict(runObj%inputs, runObj%atoms, dict)

      call dict_free(dict)
      
      ! Transfer "at" data to ART variables.
      gnrm_l = runObj%inputs%gnrm_cv
      if ( my_gnrm == 1.0d0 ) then 
         gnrm_h = runObj%inputs%gnrm_cv 
      else
         gnrm_h = my_gnrm
      end if
      ! The BigDFT restart structure.
      allocate(runObj%rst)
      call init_restart_objects(me, runObj%inputs, runObj%atoms, runObj%rst, subname)

      runObj%radii_cf = f_malloc_ptr((/ runObj%atoms%astruct%ntypes, 3 /),id='runObj%radii_cf')
      call read_radii_variables(runObj%atoms, runObj%radii_cf, &
           & runObj%inputs%crmult, runObj%inputs%frmult, runObj%inputs%projrad)

   END SUBROUTINE bigdft_init_art


   !> ART calcforce_bigdft
   !! Calculation of forces
   subroutine calcforce_bigdft( posa, forca, boxl, energy, evalf_number, conv )

      implicit none

      !Arguments

      real(kind=8), intent(in),  dimension(3*runObj%atoms%astruct%nat), target :: posa
      real(kind=8), intent(out), dimension(3*runObj%atoms%astruct%nat), target :: forca
      real(kind=8), dimension(3), intent(inout)           :: boxl
      real(kind=8), intent(out)                           :: energy
      integer,      intent(inout)                         :: evalf_number
      logical,      intent(in)                            :: conv

      !Local variables
      integer  :: infocode, i, ierror 
      type(DFT_global_output) :: outs
      !_______________________

      if ( conv ) then                    ! Convergence criterion for the wavefunction optimization
         runObj%inputs%gnrm_cv = gnrm_h              ! in Lanczos procedure.              
      else 
         runObj%inputs%gnrm_cv = gnrm_l                                    
      end if
      ! We transfer acell into 'at'
      runObj%atoms%astruct%cell_dim = boxl/Bohr_Ang
      ! Need to transform posa into xcart
      ! 1D -> 2D array
      do i = 1, runObj%atoms%astruct%nat, 1
         runObj%atoms%astruct%rxyz(:, i) = (/ posa(i), posa(runObj%atoms%astruct%nat + i), &
              & posa(2 * runObj%atoms%astruct%nat + i) /) / Bohr_Ang
      end do

      call init_global_output(outs, runObj%atoms%astruct%nat)

      if ( first_time ) then              ! This is done by default at the beginning.


         runObj%inputs%inputPsiId = 0
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call call_bigdft(runObj, outs, nproc, me, infocode )
         evalf_number = evalf_number + 1

         runObj%inputs%inputPsiId = 1
         initialised   = .true.
         first_time    = .False.
         new_wf        = .False.

      else 

         if ( .not. initialised ) then
            write(0,*) "No previous call to bigdft_init_art(). On strike, refuse to work."
            write(*,*) "No previous call to bigdft_init_art(). On strike, refuse to work."
            stop
         end if

         if ( new_wf ) then               ! if true,  we do not use the previously 
            ! calculated wave function.
            runObj%inputs%inputPsiId = 0
         else 
            runObj%inputs%inputPsiId = 1
         end if

         ! Get into BigDFT
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call call_bigdft(runObj, outs, nproc, me, infocode )
         evalf_number = evalf_number + 1

      end if
      ! Energy in eV 
      energy = outs%energy * ht2ev
      ! box in ang
      boxl = runObj%atoms%astruct%cell_dim * Bohr_Ang

      ! zero forces for blocked atoms:
      ! This was already done in clean_forces (forces.f90).
      ! But, up to now, ART only works with totally frozen atoms
      ! ( i.e "f" ). Therefore, this is a safe action.
      do i = 1, runObj%atoms%astruct%nat, 1
         if ( runObj%atoms%astruct%ifrztyp(i) /= 0  .or. in_system(i) /= 0 ) outs%fxyz(:,i) = 0.0d0 
      end do 

      call center_f( outs%fxyz, runObj%atoms%astruct%nat )     ! We remove the net force over our free atomos.

      do i = 1, runObj%atoms%astruct%nat, 1                    ! Forces into ev/ang and in 1D array.
         forca( i )                        = outs%fxyz(1, i) * ht2ev / Bohr_Ang
         forca( runObj%atoms%astruct%nat + i )     = outs%fxyz(2, i) * ht2ev / Bohr_Ang
         forca( 2 * runObj%atoms%astruct%nat + i ) = outs%fxyz(3, i) * ht2ev / Bohr_Ang
      end do

      call deallocate_global_output(outs)

   END SUBROUTINE calcforce_bigdft


   !> Minimise geometry
   subroutine mingeo( posa, forca, boxl, evalf_number, total_energy, success )

      implicit none

      !Arguments
      real(kind=8), intent(inout), dimension(3*atoms_all%astruct%nat) :: posa
      real(kind=8), intent(in),    dimension(3*atoms_all%astruct%nat), target :: forca
      real(kind=8), intent(inout), dimension(3)     :: boxl
      integer,      intent(inout)                   :: evalf_number
      real(kind=8), intent(out)                     :: total_energy
      logical,      intent(out)                     :: success

      !Local variables
      integer :: i, ierror, ncount_bigdft
      type(DFT_global_output) :: outs

      if ( .not. initialised ) then
         write(0,*) "No previous call to bigdft_init_art(). On strike, refuse to work."
         write(*,*) "No previous call to bigdft_init_art(). On strike, refuse to work."
         stop
      end if

      success = .True.                    ! success will be .False. if:
      !ncount_bigdft > in%ncount_cluster_x-1

      runObj%inputs%gnrm_cv = gnrm_l                 ! For relaxation, we use always the default value in input.dft

      runObj%atoms%astruct%cell_dim = boxl/Bohr_Ang
      ! Need to transform posa into xcart
      ! 1D -> 2D array
      do i = 1, runObj%atoms%astruct%nat, 1
         runObj%atoms%astruct%rxyz(:, i) = (/ posa(i), posa(runObj%atoms%astruct%nat + i), &
              & posa(2 * runObj%atoms%astruct%nat + i) /) / Bohr_Ang
      end do

      call init_global_output(outs, runObj%atoms%astruct%nat)
      do i = 1, runObj%atoms%astruct%nat, 1
         outs%fxyz(:, i) = (/ forca(i), forca(runObj%atoms%astruct%nat + i), &
              & forca(2 * runObj%atoms%astruct%nat + i) /) * Bohr_Ang / ht2ev
      end do

      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call geopt(runObj, outs, nproc, me, ncount_bigdft )
      evalf_number = evalf_number + ncount_bigdft 
      if (ncount_bigdft > runObj%inputs%ncount_cluster_x-1) success = .False.

      total_energy = outs%energy * ht2ev
      ! box in ang
      boxl = runObj%atoms%astruct%cell_dim * Bohr_Ang
      ! Positions into ang.
      do i = 1, runObj%atoms%astruct%nat, 1
         posa(i)                        = runObj%atoms%astruct%rxyz(1, i) * Bohr_Ang
         posa(runObj%atoms%astruct%nat + i)     = runObj%atoms%astruct%rxyz(2, i) * Bohr_Ang
         posa(2 * runObj%atoms%astruct%nat + i) = runObj%atoms%astruct%rxyz(3, i) * Bohr_Ang
      end do

      call deallocate_global_output(outs)

   END SUBROUTINE mingeo


   !> Routine to finalise all BigDFT stuff
   subroutine bigdft_finalise ( )
      implicit none
      !Local variable
      character(len=*), parameter :: subname='bigdft_finalise'

      call run_objects_free(runObj, subname)
   END SUBROUTINE bigdft_finalise

   !> Removes the net force taking into account the blocked atoms
   subroutine center_f( vector, natoms )

      implicit none

      !Arguments
      integer, intent(in) :: natoms 
      real(kind=8), dimension(3,natoms), intent(inout), target :: vector

      !Local variables
      integer :: i
      integer :: natoms_f                              ! degrees of freedom
      real(kind=8) :: xtotal, ytotal, ztotal
      logical, dimension(natoms) :: mask

      ! degrees of freedom 
      mask = runObj%atoms%astruct%ifrztyp .eq. 0 .and. in_system .eq. 0
      natoms_f = count(mask)

      xtotal = 0.0d0
      ytotal = 0.0d0
      ztotal = 0.0d0

      ! Do over free atoms ( although, the frozen ones add zero ) 
      do i = 1, natoms
         if ( mask(i) ) then
            xtotal = xtotal + vector(1,i)
            ytotal = ytotal + vector(2,i)
            ztotal = ztotal + vector(3,i)
         end if 
      enddo 

      if (iproc==0) then 
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) net force over free region ', sqrt(xtotal**2 + ytotal**2 + ztotal**2)
      end if 

      ! The average is only over the degrees of freedom in each direction 
      xtotal = xtotal / natoms_f 
      ytotal = ytotal / natoms_f
      ztotal = ztotal / natoms_f

      if (iproc==0) then 
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) residual force along x(pa)=', xtotal  
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) residual force along y(pa)=', ytotal  
         write(*,'(a,1x,1p,e24.17,0p)') 'CENTER: 1) residual force along z(pa)=', ztotal  
      end if

      ! Do over free atoms 
      do i = 1, natoms, 1
         if ( mask(i) ) then
            vector(1,i) = vector(1,i) - xtotal
            vector(2,i) = vector(2,i) - ytotal
            vector(3,i) = vector(3,i) - ztotal
         end if 
      end do 

   END SUBROUTINE center_f


   subroutine copy_atoms_object(atoms1,atoms2,rxyz,nat,total_nb_atoms,posquant)
      use module_types
      use module_base
      implicit none


      type(atoms_data),intent(in) :: atoms1 
      type(atoms_data),intent(out) :: atoms2
      real(gp), dimension(:,:), pointer     :: rxyz
      integer, intent(in)         :: nat
      integer,intent(in) :: total_nb_atoms
      real(8), dimension(3*total_nb_atoms),intent(in) :: posquant
      integer :: i_stat

      character(len=*), parameter :: subname='copy_atoms_object'



      atoms2%astruct%units   = atoms1%astruct%units
      atoms2%astruct%nat  =    nat
      atoms2%astruct%cell_dim(1)  = atoms1%astruct%cell_dim(1)*Bohr_Ang
      atoms2%astruct%cell_dim(2)  = atoms1%astruct%cell_dim(2)*Bohr_Ang
      atoms2%astruct%cell_dim(3)  = atoms1%astruct%cell_dim(3)*Bohr_Ang

      atoms2%astruct%geocode = atoms1%astruct%geocode


      atoms2%astruct%input_polarization = f_malloc_ptr(atoms2%astruct%nat,id='atoms2%astruct%input_polarization')

      !also the spin polarisation and the charge are is fixed to zero by default
      !this corresponds to the value of 100
      !RULE natpol=charge*1000 + 100 + spinpol

      atoms2%astruct%input_polarization(:)=100
      atoms2%astruct%input_polarization(:) = atoms1%astruct%input_polarization(1:nat)

      atoms2%astruct%ifrztyp = f_malloc_ptr(atoms2%astruct%nat,id='atoms2%astruct%ifrztyp')


      !this array is useful for frozen atoms
      !no atom is frozen by default
      atoms2%astruct%ifrztyp(:)=0
      atoms2%astruct%ifrztyp(:)= atoms1%astruct%ifrztyp(1:nat)

      rxyz = f_malloc_ptr((/ 3, atoms2%astruct%nat /),id='rxyz')


      rxyz(1,:) = posquant(1:nat)
      rxyz(2,:) = posquant(1+total_nb_atoms:nat+total_nb_atoms)
      rxyz(3,:) = posquant(1+total_nb_atoms+total_nb_atoms:nat+total_nb_atoms+total_nb_atoms)

      atoms2%astruct%ntypes  = atoms1%astruct%ntypes




      atoms2%astruct%iatype = f_malloc_ptr(atoms2%astruct%nat,id='atoms2%astruct%iatype')
      atoms2%astruct%iatype(:) = atoms1%astruct%iatype(1:nat)

      allocate(atoms2%astruct%atomnames(atoms2%astruct%ntypes+ndebug),stat=i_stat)
      call memocc(i_stat,atoms2%astruct%atomnames,'atoms%astruct%atomnames',subname)
      atoms2%astruct%atomnames(1:atoms2%astruct%ntypes)=atoms1%astruct%atomnames(1:atoms2%astruct%ntypes)
      atoms2%astruct%inputfile_format = atoms1%astruct%inputfile_format

   END SUBROUTINE copy_atoms_object




   !Right now, this routine only prepares a region made mostly of Si
   !It only passivates the bottom
   !Modifications could be made to passivate the sides as well
   subroutine prepare_quantum_atoms_Si(atoms,posquant,nat)
      use module_types
      use module_base
      use defs
      implicit none

      type(atoms_data),intent(inout) :: atoms
      real(8), dimension(3*natoms),intent(out) :: posquant
      integer, intent(out)           :: nat

      real(8) :: x_min,x_max,y_min,y_max,z_min,z_max
      integer :: i,j,k
      logical, dimension(natoms) :: is_at_quantum
      real(8) :: xij,yij,zij,rij2
      logical :: have_hydro
      character(len=20), dimension(100) :: atomnames
      integer :: i_stat, i_all
      integer :: hydro_atom_type

      integer, dimension(natoms) :: numnei
      integer, dimension(natoms,maxnei) :: nei
      real(8), dimension(3) :: invbox

      character(len=*), parameter :: subname='prepare_quantum_atoms_Si'

      invbox = 1.0d0/box

      is_at_quantum = .false. !vectorial operation
      nat = 0
      posquant = 0.0d0

      !First : define a box which encompasses all of the atoms of the quantum region
      x_max = -10000000000000000.0d0
      y_max = -10000000000000000.0d0
      z_max = -10000000000000000.0d0

      x_min = 10000000000000000.0d0
      y_min = 10000000000000000.0d0
      z_min = 10000000000000000.0d0

      have_hydro = .false.
      if (passivate) then
         do i = 1,5   !!!type_name only goes to 5 and is hard coded!!!
            if (type_name(i) == "H") then
               have_hydro = .true.
               hydro_atom_type = i
            endif      
         end do
         if (.not. have_hydro) then
            atomnames(1:atoms%astruct%ntypes) = atoms%astruct%atomnames(1:atoms%astruct%ntypes)
            atoms%astruct%ntypes = atoms%astruct%ntypes +1
            i_all=-product(shape(atoms%astruct%atomnames))*kind(atoms%astruct%atomnames)
            deallocate(atoms%astruct%atomnames, stat=i_stat)
            call memocc(i_stat, i_all, 'atoms%astruct%atomnames', subname)
            atomnames(atoms%astruct%ntypes) = "H"
            atoms%astruct%atomnames(1:atoms%astruct%ntypes) = atomnames(1:atoms%astruct%ntypes)
            hydro_atom_type = atoms%astruct%ntypes
         endif
      endif


      call neighbours(natoms,pos,box,boundary,maxnei,numnei, nei)

      do i = 1,nbr_quantum
         is_at_quantum(i) = .true.
         nat = nat + 1
         if ( pos(i) < x_min) x_min = pos(i) -0.02d0
         if ( pos(i+natoms) < y_min) y_min = pos(i+natoms) -0.02d0
         if ( pos(i+natoms+natoms) < z_min) z_min = pos(i+natoms+natoms) -0.02d0

         if ( pos(i) > x_max) x_max = pos(i) +0.02d0
         if ( pos(i+natoms) > y_max) y_max = pos(i+natoms)  +0.02d0
         if ( pos(i+natoms+natoms) > z_max) z_max = pos(i+natoms+natoms) +0.02d0

         posquant(i) = pos(i)
         posquant(i+natoms) = pos(i+natoms)
         posquant(i+natoms+natoms) = pos(i+natoms+natoms)
      enddo

      do i = 1,nbr_quantum
         if (passivate) then 
            do j = 1,numnei(i)
               k = nei(i,j)
               if ( .not. is_at_quantum(k)) then 
                  xij = pos(k)-pos(i) - box(1) * nint((pos(k)-pos(i))*invbox(1))
                  yij = pos(k+natoms)-pos(i+natoms)
                  zij = pos(k+2*natoms)-pos(i+2*natoms) - box(3) * nint((pos(k+2*natoms)-pos(i+2*natoms))*invbox(3))
                  rij2 = xij*xij + yij*yij + zij*zij
                  if (rij2 .lt. 2.5d0*2.5d0) then
                     nat = nat + 1
                     posquant(nat) = pos(i) + 0.5d0*xij
                     posquant(nat+natoms) = pos(i+natoms) + 0.5d0*yij
                     posquant(nat+natoms+natoms) = pos(i+2*natoms) + 0.5d0*zij !we passivate with hydrogene at this distance
                     atoms%astruct%iatype(nat) = hydro_atom_type
                     atoms%astruct%ifrztyp(i) = 1  !this atom is frozen
                     atoms%astruct%ifrztyp(nat) = 1  !this one as well
                  endif
               endif   
            enddo
         endif
      enddo

   END SUBROUTINE prepare_quantum_atoms_Si


   subroutine check_force_clean_wf( posa, boxl, evalf_number, total_energy, success )

      implicit none

      real(kind=8), intent(inout), dimension(3*atoms_all%astruct%nat) :: posa
      real(kind=8), intent(inout), dimension(3)     :: boxl
      integer,      intent(inout)                   :: evalf_number
      real(kind=8), intent(out)                     :: total_energy
      logical,      intent(out)                     :: success

      !Local variables
      integer      :: infocode, i, ierror, ncount_bigdft 
      type(DFT_global_output) :: outs
      real(gp)     ::  fmax, fnrm
      !_______________________

      runObj%inputs%inputPsiId = 0 
      runObj%inputs%gnrm_cv = gnrm_l 
      ! We transfer acell into 'at'
      runObj%atoms%astruct%cell_dim = boxl/Bohr_Ang

      do i = 1, runObj%atoms%astruct%nat, 1
         runObj%atoms%astruct%rxyz(:, i) = (/ posa(i), posa(runObj%atoms%astruct%nat + i), &
              & posa(2 * runObj%atoms%astruct%nat + i) /) / Bohr_Ang
      end do

      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call call_bigdft(runObj, outs, nproc, me, infocode )
      evalf_number = evalf_number + 1
      runObj%inputs%inputPsiId = 1

      call fnrmandforcemax(outs%fxyz,fnrm,fmax, outs%fdim)

      if ( fmax > runObj%inputs%forcemax ) then

         if ( iproc == 0 ) then
            write(*,*) 'BART:check_force_clean_wf'
            write(*,*) 'BART: fmax =', fmax,'H/Bohr. We relax again!' 
         end if

         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call geopt(runObj, outs, nproc, me, ncount_bigdft )
         evalf_number = evalf_number + ncount_bigdft 
         if (ncount_bigdft > runObj%inputs%ncount_cluster_x-1) success = .False.

         ! and we clean again here
         runObj%inputs%inputPsiId = 0 
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call call_bigdft(runObj, outs, nproc, me, infocode )
         evalf_number = evalf_number + 1
         runObj%inputs%inputPsiId = 1

         total_energy = outs%energy * ht2ev
         ! box in ang
         boxl = runObj%atoms%astruct%cell_dim * Bohr_Ang
         ! Positions into ang.
         do i = 1, runObj%atoms%astruct%nat, 1
            posa(i)                        = runObj%atoms%astruct%rxyz(1, i) * Bohr_Ang
            posa(runObj%atoms%astruct%nat + i)     = runObj%atoms%astruct%rxyz(2, i) * Bohr_Ang
            posa(2 * runObj%atoms%astruct%nat + i) = runObj%atoms%astruct%rxyz(3, i) * Bohr_Ang
         end do

      end if 

      call deallocate_global_output(outs)

   END SUBROUTINE check_force_clean_wf


END MODULE bigdft_forces
