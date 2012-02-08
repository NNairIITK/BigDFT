!> @file
!! parameters from the interatomic potential
!! @author
!!   Fedwa El-Mellouhi, Normand Mousseau February 2006
!!   Copyright (C) 2010-2012 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> This module define the most important variables in the program 
module parameters 

   implicit none
   save

   ! Some fixed parameters  

   integer, parameter :: idim = 3                !< Dimension of the problem
   integer, parameter :: species = 1             !< Number of species
   integer, parameter :: max_spec = 5            !< Max number of atomic  species admissible
   integer, parameter :: max_neighbor_cells = 26 !< Max number of neighboring cells up to 3D
   integer, parameter :: max_cluster = 512       !< Max number of atoms in the cluster surrounding each atom
   integer, parameter :: NSPECIES = species

   real(kind=8), parameter :: PI = 3.14159265358979D0

END MODULE parameters


!> Module which define the radius cut of SW potential
module cutoff
   use parameters
   implicit none

   real(kind=8), dimension(2,2):: rcut2
END MODULE cutoff


!> This module defines the parameters for the Stillinger-Weber potential
module SWpotential
   use parameters 
   use cutoff 
   implicit none


   real(kind=8),dimension(:,:),allocatable :: SW_parameters
   integer,dimension(:), allocatable ::    angles_to_fit !angles to fit per atom
   integer,dimension(:,:,:),allocatable :: jk_of_angle
   real(kind=8),dimension(:,:), allocatable ::  optimum_angle
   integer      ::                         len_vect_to_fit
   integer      ::                         configs_to_fit

   ! 1   : real(kind=8), parameter :: SIGMA = 2.095037426D0
   ! 2   : real(kind=8), parameter :: A =  7.049556277D0
   ! 3   : real(kind=8), parameter :: ALPHA = 1.8d0 
   ! 4   : real(kind=8), parameter :: BETA = 0.60222455844D0
   !  real(kind=8), parameter :: EPSILON = 2.16823D0* 1.07d0 !  FE Balamane et al. PRB 46, 2250
   ! 5   : real(kind=8), parameter :: EPSILON = 2.16823D0 ! Original SW
   ! 6  :  integer, parameter :: P = 4
   ! 7  :  real(kind=8), parameter :: GAMMA = 1.2d0
   ! 8  :  real(kind=8), parameter :: LAMBDA = (21.0d0 * 1.0d0)
   ! 9 :    real(kind=8), parameter :: ONE_THIRD = 1.0d0/3.0d0

   real(kind=8)  :: SIGMA 
   real(kind=8)  :: A 
   real(kind=8)  :: ALPHA 
   real(kind=8)  :: BETA
   !  real(kind=8), parameter :: EPSILON = 2.16823D0* 1.07d0 !  FE Balamane et al. PRB 46, 2250
   real(kind=8)  :: EPSILON  ! Original SW
   real(kind=8)  :: P 
   real(kind=8)  :: GAMMA
   real(kind=8)  :: LAMBDA 
   real(kind=8) :: ONE_THIRD = 1.0d0/3.0d0
   real(kind=8) :: RCUT
   real(kind=8) :: A_EPS

   real(kind=8), dimension(2,2) :: A0, S0

   real(kind=8), allocatable,dimension(:,:) ::force_ref_fit  !!reference force for fit
   real(kind=8), allocatable,dimension(:,:) ::force_work_fit  !!working forces for fit
   real(kind=8), allocatable,dimension(:,:) :: pos_ref_fit  !!reference position for fit

END MODULE  SWpotential


!> This subroutine initializes the SW parameters 
subroutine init_potential_SW()
   use SWpotential
   use defs, only : natoms

   implicit none

   !the 10th bond is the reference
   if (.not. allocated(SW_parameters)) allocate(SW_parameters(natoms,8))
   if (.not. allocated( angles_to_fit)) allocate(angles_to_fit(natoms))!angles to fit per atom
   if (.not. allocated( jk_of_angle))  allocate(jk_of_angle(natoms,100,2)) !up to 100 angles to fit
   if (.not. allocated( optimum_angle)) allocate(optimum_angle(natoms,100))

   angles_to_fit = 0
   jk_of_angle = 0
   configs_to_fit = 4

   if (.not. allocated(force_ref_fit)) allocate(force_ref_fit(3*natoms,configs_to_fit))
   if (.not. allocated(pos_ref_fit)) allocate(pos_ref_fit(3*natoms,configs_to_fit))
   if (.not. allocated(force_work_fit)) allocate(force_work_fit(3*natoms,configs_to_fit))

   ! SW_parameters(:,1) =  2.095037426D0
   ! we change sigma so that cell size is 5.4 and not 5.43 to emulate LDA
   SW_parameters(:,1) = 2.095037426d0*5.465d0/5.43d0

   SW_parameters(:,2) =   7.049556277D0
   SW_parameters(:,3) = 1.8d0 
   SW_parameters(:,4) =  0.60222455844D0
   !  real(kind=8), parameter :: EPSILON = 2.16823D0* 1.07d0 !  FE Balamane et al. PRB 46, 2250
   SW_parameters(:,5) =  2.16823D0 ! Original SW
   SW_parameters(:,6) =  4.0d0
   SW_parameters(:,7) =  1.2d0
   SW_parameters(:,8) =  (21.0d0 * 1.0d0)
   optimum_angle = 1.0d0/3.0d0
   !SW_parameters(:,9:29) = (1.0d0/3.0d0)
   RCUT = SW_parameters(1,1) * SW_parameters(1,3)

   if (NSPECIES == 2) then
      S0(1,1)=3.6
      S0(1,2)=3.6
      S0(2,1)=3.6
      S0(2,2)=3.6
      ! for GaAs A0(1,1)=1.2*EPSILON
      ! or       A0(1,1)=1.5*EPSILON 
      A0(1,1)=0.0*EPSILON
      A0(2,2)=0.0*EPSILON
      A0(1,2)=0.0*EPSILON
      A0(2,1)=0.0*EPSILON

      rcut2(1,1)=RCUT*RCUT
      rcut2(1,2)=RCUT*RCUT
      rcut2(2,1)=RCUT*RCUT
      rcut2(2,2)=RCUT*RCUT
   else if (NSPECIES == 1 ) then
      S0(1,1) = 3.6d0
      A0(1,1) = 0.0d0
      rcut2      = 0.0d0
      rcut2(1,1) = RCUT*RCUT
   else
      stop
   endif
END SUBROUTINE init_potential_SW


subroutine reset_SW_potential()
   use SWpotential
   implicit none

   !  SW_parameters(:,1) =  2.095037426D0
   !  we change sigma so that cell size is 5.4 and not 5.43 to emulate LDA
   !  5.465 for GGA
   SW_parameters(:,1) = 2.095037426d0*5.465d0/5.43d0
   SW_parameters(:,2) =   7.049556277D0
   SW_parameters(:,3) = 1.8d0
   SW_parameters(:,4) =  0.60222455844D0
   !  real(kind=8), parameter :: EPSILON = 2.16823D0* 1.07d0 !  FE Balamane et al. PRB 46, 2250
   SW_parameters(:,5) =  2.16823D0 ! Original SW
   SW_parameters(:,6) =  4.0d0
   SW_parameters(:,7) =  1.2d0
   SW_parameters(:,8) =  (21.0d0 * 1.0d0)
   optimum_angle = 1.0d0/3.0d0
   !SW_parameters(:,9:29) = (1.0d0/3.0d0)
END SUBROUTINE


!> This subroutine fits the SW classical forces near the quantum region to a hybrid force
!! calculation
subroutine fit_SW_potential()
   use defs
   use SWpotential
   use random
   implicit none

   integer :: my_counter,trash_evalf
   integer :: i,j

   real(kind=8) , dimension(:),allocatable :: vector_to_fit

   real(kind=8) :: trash_energy
   real(kind=8),dimension(3) :: invbox
   real(kind=8), parameter :: tolerance = 0.40d0
   integer,parameter :: max_mc_steps = 20
   integer :: ierror
   real(kind=8) :: value,dummy
   real(kind=8) :: diff_square_force
   real(kind=8), dimension(0:nproc-1,0:nproc-1) :: sendbuff
   real(kind=8),dimension(0:nproc-1,0:nproc-1) :: recvbuff
   real(kind=8),dimension(3*natoms)     :: box_vec
   integer,dimension(natoms) :: numnei
   integer,dimension(natoms,maxnei) :: nei
   real(kind=8) :: ran3

   sendbuff = 9999999999.0d0
   recvbuff = 9999999999.0d0

   invbox = 1.0d0 / box

   if (iproc == 0) write(*,*) "we are at the beggining of the routine. atoms to fit :", nbr_to_fit
   trash_evalf = 0

   pos_ref_fit(:,1) = pos
   do i = 1,configs_to_fit

      energy_type = "BSW"

      call calcforce( NATOMS, pos_ref_fit(:,i), boxref, force_ref_fit(:,i),&
         &   trash_energy, trash_evalf, .false. )
      do j = 1,3*natoms
         dummy = (0.5d0-ran3())*0.1d0
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call MPI_Bcast(dummy,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

         pos_ref_fit(j,i+1) = pos(j) + dummy
      enddo
      energy_type = "OTF"
      call calcforce( NATOMS, pos_ref_fit(:,i), boxref, force_work_fit(:,i),&
         &   trash_energy, trash_evalf, .false. )

   enddo

   box_vec(1:natoms) = box(1)
   box_vec(1+natoms:2*natoms) = box(2)
   box_vec(1+natoms+natoms:3*natoms) = box(3)
   nei = 0 
   numnei  = 0
   if ( boundary .eq. 'P') pos = modulo(pos, box_vec) !apply PBC implement surface later on
   if ( boundary .eq. 'S') pos(1:natoms) = modulo(pos(1:natoms), box_vec(1:natoms))
   if ( boundary .eq. 'S') pos(2*natoms+1:3*natoms) = modulo(pos(2*natoms+1:3*natoms), box_vec(2*natoms+1:3*natoms))
   call neighbours(NATOMS,pos,box,boundary,maxnei,numnei,nei)

   call neighbour_angles(numnei,nei)

   len_vect_to_fit = 0

   do i = 1,natoms
      if(should_fit(i)) then
         len_vect_to_fit = len_vect_to_fit+ angles_to_fit(i)+8               
      endif
   enddo

   if (.not. allocated(vector_to_fit)) allocate(vector_to_fit(len_vect_to_fit))

   my_counter = 0

   do i = 1,natoms
      if(should_fit(i)) then
         vector_to_fit(1+my_counter:8+my_counter) = SW_parameters(i,1:8)
         my_counter = my_counter + 8
         do j = 1,angles_to_fit(i)
            vector_to_fit(1+my_counter) = optimum_angle(i,j)
            my_counter = my_counter + 1
         enddo
      endif
   enddo


   do i = 1,40

      call fit_sd(len_vect_to_fit,vector_to_fit,tolerance,value,numnei,nei)
      call mc_fit(len_vect_to_fit,vector_to_fit,tolerance,value,numnei,nei)
      if (iproc == 0) write(*,*) "current_value",value 
      sendbuff = value
      call MPI_Alltoall(sendbuff,nproc,MPI_REAL8, &
         &   recvbuff,nproc,MPI_REAL8, &
         &   MPI_COMM_WORLD,ierror)

      do j = 0,nproc-1
         if (recvbuff(iproc,j) < tolerance) then
            value = recvbuff(iproc,j)
            call MPI_Bcast(value,1,MPI_REAL8,j,MPI_COMM_WORLD,ierror)
            call MPI_Bcast(vector_to_fit,len_vect_to_fit,MPI_REAL8,j,MPI_COMM_WORLD,ierror)
            exit
         endif
      enddo
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      if (value < tolerance) then
         dummy = diff_square_force(vector_to_fit) !only there so SW_parameters = vector
         exit
      endif


   enddo

   if (value .ge. tolerance) then
      call MPI_Bcast(value,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(vector_to_fit,len_vect_to_fit,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      dummy = diff_square_force(vector_to_fit)
   endif

   do i = 1,natoms
      do j=1,8
         if (iproc ==0 ) write(1030+j,*) i,SW_parameters(i,j)
      enddo
   enddo

   !the SW_parameters are equaled to the vector_to_fit in the function

   if (allocated(vector_to_fit)) deallocate(vector_to_fit)

END SUBROUTINE fit_SW_potential


subroutine mc_fit(vec_size,vector,tolerance,current_value,numnei,nei)
   use defs,only : iproc,maxnei,natoms,should_fit
   use SWpotential,only : angles_to_fit
   use random
   implicit none

   integer,intent(in) :: vec_size
   real(kind=8), dimension(vec_size), intent(inout) :: vector
   real(kind=8), intent(in) :: tolerance
   real(kind=8), intent(out) :: current_value
   integer,dimension(natoms),intent(in) :: numnei
   integer,dimension(natoms,maxnei),intent(in) :: nei

   integer :: i,param,my_counter,this_atom
   real(kind=8) :: diff_square_force
   real(kind=8) :: diff_square_force_one

   real(kind=8) , dimension(vec_size) :: tmp_vector

   real(kind=8) :: ran3
   real(kind=8) :: ran_number
   real(kind=8) :: difference

   real(kind=8) :: step_size = 0.05d0
   integer,parameter :: max_mc_steps = 25


   difference = diff_square_force(vector)
   current_value = difference
   tmp_vector = vector
   do i = 1,max_mc_steps
      my_counter = 0
      this_atom = 1

      do param = 1,vec_size !we change all the parameters randomly
         ran_number = ran3() - 0.5d0
         vector(param) = vector(param) + ran_number*step_size*abs(vector(param))

         difference = diff_square_force_one(vector,numnei,nei,this_atom)

         if ( difference < current_value ) then
            tmp_vector(param) = vector(param)
            current_value = difference
            difference = diff_square_force(vector)
         else
            vector(param) = tmp_vector(param)
            difference = current_value
         endif
         my_counter = my_counter + 1
         if (my_counter == 8 + angles_to_fit(this_atom)) then
            my_counter = 0

            do
               this_atom = this_atom +1
               if (should_fit(this_atom)) exit
            enddo
         endif

      enddo

      if (current_value < tolerance) exit

   enddo

END SUBROUTINE mc_fit


subroutine fit_sd(vec_size,vector,tolerance,current_value,numnei,nei)
   use defs,only : iproc,maxnei,natoms,should_fit
   use SWpotential,only :angles_to_fit
   implicit none

   integer,intent(in) :: vec_size
   real(kind=8), dimension(vec_size), intent(inout) :: vector
   real(kind=8), intent(in) :: tolerance
   real(kind=8), intent(out) :: current_value

   integer,dimension(natoms),intent(in) :: numnei
   integer,dimension(natoms,maxnei),intent(in) :: nei

   real(kind=8) :: diff_square_force

   real(kind=8) :: difference
   real(kind=8), dimension(vec_size) :: tmp_vector,grad_vector
   real(kind=8) :: step
   integer :: iter,i
   integer, parameter :: MAX_ITER = 50
   real(kind=8), parameter :: STEPSIZE = 1.0D0  ! Size of displacement

   call deriv_diff_square_force(vector,grad_vector,numnei,nei)
   difference = diff_square_force(vector)
   current_value = difference 
   tmp_vector = vector
   step = STEPSIZE


   do iter = 1, MAX_ITER
      !tmp_vector = vector
      do i =1,vec_size
         if (abs(step*grad_vector(i)) < 0.1d0) then
            vector(i) = vector(i)-step*grad_vector(i)*abs(vector(i))
            elseif (grad_vector(i) > 0) then
            vector(i) = vector(i)-0.05d0*vector(i)
            elseif (grad_vector(i) < 0) then
            vector(i)= vector(i) +0.05d0*vector(i)
         endif
      enddo
      difference =  diff_square_force(vector)

      if(difference < current_value  ) then    
         current_value = difference
         step = 2.0d0 * step         
         tmp_vector = vector
         call deriv_diff_square_force(vector,grad_vector,numnei,nei)
      else 
         vector = tmp_vector
         step = 0.6d0*step
      endif

      if(difference < tolerance) exit
   end do
END SUBROUTINE fit_sd


function diff_square_force_one(P,numnei,nei,this_atom)

   use SWpotential,only : SW_parameters,force_ref_fit,len_vect_to_fit,optimum_angle,&
      &   angles_to_fit,pos_ref_fit,configs_to_fit,force_work_fit
   use defs,only : should_fit,natoms,pos,boxref,iproc,nbr_to_fit,maxnei,force
   implicit none

   real(kind=8) :: diff_square_force_one
   real(kind=8),dimension(len_vect_to_fit),intent(in) :: P
   integer,dimension(natoms),intent(in) :: numnei
   integer,dimension(natoms,maxnei),intent(in) :: nei
   integer,intent(in) ::this_atom

   integer i,j
   integer my_counter,that_atom
   integer ::  trash_evalf
   real(kind=8), dimension(3*natoms,configs_to_fit) :: force_tempo,tmp_force
   interface
      subroutine SWcalczone(nat,posa,boxl,tmp_force, this_atom,numnei,nei)
         use defs, only : maxnei
         integer, intent(in)                               :: nat
         real(kind=8), intent(in), dimension(3*nat) :: posa
         real(kind=8), dimension(3), intent(inout)          :: boxl
         integer, intent(in) :: this_atom
         real(kind=8), intent(out), dimension(3*nat), target:: tmp_force
         integer, dimension(nat),intent(in) :: numnei 
         integer, dimension(nat,maxnei),intent(in) :: nei 
      END SUBROUTINE SWcalczone
   end interface

   trash_evalf = 0
   my_counter = 0
   force_tempo = force_work_fit
   do i = 1, natoms
      if (should_fit(i)) then

         SW_parameters(i,1:8)=P(my_counter+1:my_counter+8)
         my_counter =my_counter + 8

         do j = 1,angles_to_fit(i)
            optimum_angle(i,j) = P(my_counter +1)
            my_counter = my_counter+1
         enddo
      endif
   enddo

   diff_square_force_one = 0.0d0
   do j = 1,configs_to_fit   
      call SWcalczone(natoms,pos_ref_fit(:,j),boxref,tmp_force(:,j), this_atom,numnei,nei)
      force_tempo(this_atom,j) = tmp_force(this_atom,j)
      force_tempo(this_atom+natoms,j) = tmp_force(this_atom+natoms,j)
      force_tempo(this_atom+natoms+natoms,j) = tmp_force(this_atom+natoms+natoms,j)

      do i = 1,numnei(this_atom)
         that_atom = nei(this_atom,i)
         force_tempo(that_atom,j) = tmp_force(that_atom,j)
         force_tempo(that_atom+natoms,j) = tmp_force(that_atom+natoms,j)
         force_tempo(that_atom+natoms+natoms,j) = tmp_force(that_atom+natoms+natoms,j)
      enddo

      do i = 1,natoms
         diff_square_force_one = diff_square_force_one + (force_ref_fit(i,j)-force_tempo(i,j))**2 + &
            &   (force_ref_fit(i+natoms,j)-force_tempo(i+natoms,j))**2 +(force_ref_fit(i+2*natoms,j)&
            &   -force_tempo(i+2*natoms,j))**2

      enddo
   enddo

END FUNCTION diff_square_force_one


function diff_square_force(P)

   use SWpotential,only : SW_parameters,force_ref_fit,len_vect_to_fit,optimum_angle &
      &   ,angles_to_fit,pos_ref_fit,force_work_fit,configs_to_fit
   use defs,only : should_fit,natoms,pos,boxref,iproc,nbr_to_fit,maxnei,force
   implicit none

   real(kind=8) :: diff_square_force
   real(kind=8),dimension(len_vect_to_fit),intent(in) :: P


   integer i,j
   integer my_counter
   real(kind=8) :: trash_energy
   integer ::  trash_evalf

   trash_evalf = 0
   my_counter = 0

   do i = 1, natoms
      if (should_fit(i)) then

         SW_parameters(i,1:8)=P(my_counter+1:my_counter+8)
         my_counter =my_counter + 8

         do j = 1,angles_to_fit(i)
            optimum_angle(i,j) = P(my_counter +1)
            my_counter = my_counter+1
         enddo
      endif
   enddo

   diff_square_force = 0.0d0 
   do j = 1,configs_to_fit

      call calcforce( NATOMS, pos_ref_fit(:,j), boxref, force_work_fit(:,j), trash_energy, trash_evalf, .false. )

      do i = 1,natoms
         diff_square_force = diff_square_force + (force_ref_fit(i,j)-force_work_fit(i,j))**2 + &
            &   (force_ref_fit(i+natoms,j)-force_work_fit(i+natoms,j))**2 &
            &   +(force_ref_fit(i+2*natoms,j)-force_work_fit(i+2*natoms,j))**2

      enddo

   enddo

END FUNCTION diff_square_force


!> Finds the gradient vector
subroutine deriv_diff_square_force(P,DF,numnei,nei)
   use SWpotential,only : len_vect_to_fit,angles_to_fit
   use defs,only : natoms,iproc,nbr_to_fit,maxnei,should_fit
   implicit none

   real(kind=8),dimension(len_vect_to_fit),intent(in) :: P
   real(kind=8), dimension(len_vect_to_fit),intent(out) :: DF

   integer,dimension(natoms),intent(in) :: numnei
   integer,dimension(natoms,maxnei),intent(in) :: nei

   real(kind=8),dimension(len_vect_to_fit) ::P_before,P_after
   real(kind=8) :: diff_square_force_one
   real(kind=8) :: valeur1,valeur2
   integer :: i,this_atom,my_counter

   my_counter = 0
   this_atom = 1
   do i = 1,len_vect_to_fit

      my_counter = my_counter + 1
      if (my_counter == 8 + angles_to_fit(this_atom)) then
         my_counter = 0

         do
            this_atom = this_atom +1
            if (should_fit(this_atom)) exit
         enddo
      endif

      P_before(i) = P(i) - 0.001d0
      P_after(i) = P(i) + 0.001d0

      valeur2 = diff_square_force_one(P_after,numnei,nei,this_atom)
      valeur1 = diff_square_force_one(P_before,numnei,nei,this_atom)

      DF(i) = (valeur2-valeur1)/0.002d0

   enddo

END SUBROUTINE deriv_diff_square_force


!> subroutine to compute forces of one atom
subroutine SWcalczone(nat,posa,boxl,tmp_force, this_atom,numnei,nei)

   use SWpotential
   use defs, only : boundary,maxnei,iproc,MPI_COMM_WORLD

   implicit none

   integer, intent(in)                               :: nat
   real(kind=8), intent(in), dimension(3*nat) :: posa
   real(kind=8), dimension(3), intent(inout)          :: boxl
   integer, intent(in) :: this_atom
   real(kind=8), intent(out), dimension(3*nat), target:: tmp_force
   integer, dimension(nat),intent(in) :: numnei 
   integer, dimension(nat,maxnei),intent(in) :: nei 

   real(kind=8), dimension(3*nat), target :: pos_normalised
   real(kind=8), dimension(3*nat) :: box_vec
   real(kind=8), dimension(3*nat),target :: pos
   real(kind=8), dimension(3)             :: box
   integer                           :: NATOMS

   real(kind=8), dimension(:), pointer :: xa, ya, za
   real(kind=8), dimension(:), pointer :: fxa, fya, fza

   integer :: i,j, i_id, j_id, ind_j, k, k_id, ind_k
   integer :: my_counter_i
   real(kind=8) :: invsig
   real(kind=8) :: xi, yi, zi, xij, yij, zij, rij, rij2
   real(kind=8) :: twobodyenergy, threebodyenergy
   real(kind=8) :: cos_x_ij, cos_y_ij, cos_z_ij, invrij, rhoij
   real(kind=8) :: one_o_a_ij, expo, gam_o_a_ij, exp_gam_ij, r_to_minusp
   real(kind=8) :: one_o_a2, term1, term2,fact3_ij
   real(kind=8) :: ffx,ffy,ffz, term_ij, term_ik
   real(kind=8) :: dcosdxj,dcosdyj,dcosdzj,dcosdxk,dcosdyk,dcosdzk
   real(kind=8) :: dhdcos_ij,dhdcos_ik,cos_p_1o3,cos_jik
   real(kind=8) :: cos_x_ik, cos_y_ik, cos_z_ik, one_o_a_ik, gam_o_a_ik, fact     
   real(kind=8) :: xik,yik,zik,rik, rik2, invrik,rhoik
   logical :: go_on

   NATOMS = nat
   pos = posa
   box = boxl

   box_vec(1:natoms) = boxl(1)
   box_vec(1+natoms:2*natoms) = boxl(2)
   box_vec(1+natoms+natoms:3*natoms) = boxl(3)


   if ( boundary .eq. 'P') pos = modulo(pos, box_vec) !apply PBC implement surface later on
   if ( boundary .eq. 'S') pos(1:natoms) = modulo(pos(1:natoms), box_vec(1:natoms))
   if ( boundary .eq. 'S') pos(2*natoms+1:3*natoms) = modulo(pos(2*natoms+1:3*natoms), box_vec(2*natoms+1:3*natoms))


   ! Generate a neighbour list

   ! We first set-up pointers for the x, y, z components in the position and forces
   pos_normalised = pos / box_vec

   xa => pos_normalised(1:NATOMS)
   ya => pos_normalised(NATOMS+1:2*NATOMS)
   za => pos_normalised(2*NATOMS+1:3*NATOMS)

   fxa => tmp_force(1:NATOMS)
   fya => tmp_force(NATOMS+1:2*NATOMS)
   fza => tmp_force(2*NATOMS+1:3*NATOMS)

   tmp_force = 0.0d0   ! Vectorial operation

   twobodyenergy = 0.0d0
   threebodyenergy = 0.0d0

   do i = 1,natoms 

      go_on = .false.
      if (i .eq. this_atom) go_on = .true.
      do ind_j = 1,numnei(this_atom)
         j = nei(this_atom,ind_j)
         if ( i .eq. j ) go_on = .true.
         do ind_k = 1,numnei(j)
            k = nei(j,ind_k)
            if ( i .eq. k ) go_on = .true.
         enddo

      enddo

      if (.not. go_on ) cycle

      !     i_id = types(i)
      i_id = 1
      xi = xa(i)
      yi = ya(i)
      zi = za(i)
      do ind_j=1, numnei(i)
         j = nei(i,ind_j)
         j_id = 1
         !j_id = types(j)

         ! Pair interactions of i and j
         ! Distance, with periodic boundary conditions
         if (boundary == "P") then
            xij = xa(j) - xi - 1.0d0 * nint( xa(j)-xi )
            yij = ya(j) - yi - 1.0d0 * nint( ya(j)-yi )
            zij = za(j) - zi - 1.0d0 * nint( za(j)-zi )
            elseif (boundary == "S") then
            xij = xa(j) - xi - 1.0d0 * nint( xa(j)-xi )
            yij = ya(j) - yi
            zij = za(j) - zi - 1.0d0 * nint( za(j)-zi )
         endif

         ! Rescale the lengths into Angstroems
         xij = xij * box(1)
         yij = yij * box(2)
         zij = zij * box(3)

         rij2 = xij*xij + yij*yij + zij*zij

         SIGMA = (SW_parameters(i,1)+SW_parameters(j,1))*0.5d0
         ALPHA = (SW_parameters(i,3)+SW_parameters(j,3))*0.5d0

         ! Check the cut-off before proceeding
         if( rij2 < ALPHA*SIGMA*ALPHA*SIGMA ) then

            !initialize the parameters

            A = (SW_parameters(i,2) + SW_parameters(j,2))*0.5d0
            EPSILON = (SW_parameters(i,5)+SW_Parameters(j,5))*0.5d0
            GAMMA = (SW_parameters(i,7)+SW_parameters(j,7))*0.5d0
            P = (SW_parameters(i,6)+SW_parameters(j,6))*0.5d0
            BETA = (SW_parameters(i,4)+SW_parameters(j,4))*0.5d0

            A_EPS = (A * EPSILON)   
            invsig = 1.0D0 / SIGMA

            rij = sqrt(rij2)
            invrij = 1.0d0 / rij
            rhoij = rij * invsig
            cos_x_ij = xij * invrij
            cos_y_ij = yij * invrij
            cos_z_ij = zij * invrij

            ! Some useful quantities
            one_o_a_ij =1.0/(rhoij-ALPHA)
            expo=exp(one_o_a_ij)
            gam_o_a_ij=GAMMA*one_o_a_ij
            exp_gam_ij=exp(gam_o_a_ij)
            r_to_minusp=rhoij ** (-1.0d0*P)

            ! Two body energy and force 
            term1=A_EPS*(BETA*r_to_minusp-1.0d0)*expo
            one_o_a2 = one_o_a_ij * one_o_a_ij
            term2=(one_o_a2*term1+A_EPS*P*BETA*r_to_minusp*expo/rhoij)*invsig

            ! Contribution to the binary repulsive term 
            if(rij <=  S0(i_id,j_id)) then
               term1=term1+A0(i_id,j_id)*(cos(PI*rij/S0(i_id,j_id))+1.0d0)
               term2=term2+A0(i_id,j_id)*PI/S0(i_id,j_id)* sin(PI*rij/S0(i_id,j_id))
            endif
            twobodyenergy=twobodyenergy + 0.5d0*term1;

            fxa(i)=fxa(i)-term2*cos_x_ij;
            fya(i)=fya(i)-term2*cos_y_ij;
            fza(i)=fza(i)-term2*cos_z_ij;

            ! Prepare for the three body term 
            fact3_ij=gam_o_a_ij*one_o_a_ij*invsig

            do ind_k = ind_j+1, numnei(i)
               ! Triplet interaction with i in the middle; all interactions
               k = nei(i,ind_k)
               !k_id = types(k)
               k_id = 1
               ! Distance, with periodic boundary conditions
               if (boundary == "P") then
                  xik = xa(k) - xi - 1.0d0 * nint( xa(k)-xi )
                  yik = ya(k) - yi - 1.0d0 * nint( ya(k)-yi )
                  zik = za(k) - zi - 1.0d0 * nint( za(k)-zi )
                  elseif (boundary == "S") then
                  xik = xa(k) - xi - 1.0d0 * nint( xa(k)-xi )
                  yik = ya(k) - yi
                  zik = za(k) - zi - 1.0d0 * nint( za(k)-zi )
               endif

               ! Rescale the lengths into Angstroems
               xik = xik * box(1)
               yik = yik * box(2)
               zik = zik * box(3)

               rik2 = xik*xik + yik*yik + zik*zik

               SIGMA = (SW_parameters(i,1) + SW_parameters(k,1))*0.5d0
               ALPHA = (SW_parameters(i,3) + SW_parameters(k,3))*0.5d0

               ! Check whether the distance is too large 
               if (rik2<ALPHA*SIGMA*ALPHA*SIGMA)  then

                  !initialize the parameters

                  EPSILON = (SW_parameters(i,5)+SW_Parameters(k,5))*0.5d0
                  GAMMA = (SW_parameters(i,7) + SW_parameters(k,7))*0.5d0
                  LAMBDA = (SW_parameters(i,8))


                  if(angles_to_fit(i) == 0) then                 
                     ONE_THIRD = 1.0d0/3.0d0
                  else 
                     do my_counter_i=1,angles_to_fit(i)
                        if ( jk_of_angle(i,my_counter_i,1) == j .and. jk_of_angle(i,my_counter_i,2) == k ) then
                           !   if (iproc ==0) write(*,*) "looking at angle", optimum_angle(i,my_counter_i) 
                           ONE_THIRD = optimum_angle(i,my_counter_i)    
                        endif
                     enddo

                  endif

                  invsig = 1.0D0 / SIGMA

                  rik=sqrt(rik2)
                  invrik=1.0/rik
                  rhoik=rik*invsig
                  cos_x_ik=xik*invrik
                  cos_y_ik=yik*invrik
                  cos_z_ik=zik*invrik

                  ! Some useful quantities 
                  one_o_a_ik=1.0/(rhoik-ALPHA)
                  gam_o_a_ik=GAMMA*one_o_a_ik
                  fact     =EPSILON*LAMBDA*exp(gam_o_a_ik)*exp_gam_ij

                  cos_jik  =cos_x_ij*cos_x_ik+cos_y_ij*cos_y_ik+ cos_z_ij*cos_z_ik
                  cos_p_1o3=cos_jik+ONE_THIRD

                  ! Energy (added only to central atom)  
                  threebodyenergy=threebodyenergy+fact*cos_p_1o3*cos_p_1o3

                  ! Force 
                  term_ij=fact*fact3_ij*cos_p_1o3*cos_p_1o3;
                  dhdcos_ij=2*fact*cos_p_1o3;
                  term_ik=fact*gam_o_a_ik*one_o_a_ik*cos_p_1o3*cos_p_1o3/SIGMA;
                  dhdcos_ik=2*fact*cos_p_1o3;

                  dcosdxj=(cos_x_ik-cos_jik*cos_x_ij)*invrij;
                  dcosdyj=(cos_y_ik-cos_jik*cos_y_ij)*invrij;
                  dcosdzj=(cos_z_ik-cos_jik*cos_z_ij)*invrij;

                  dcosdxk=(cos_x_ij-cos_jik*cos_x_ik)*invrik;
                  dcosdyk=(cos_y_ij-cos_jik*cos_y_ik)*invrik;
                  dcosdzk=(cos_z_ij-cos_jik*cos_z_ik)*invrik;

                  ffx=term_ij*cos_x_ij-dhdcos_ij*dcosdxj;
                  ffy=term_ij*cos_y_ij-dhdcos_ij*dcosdyj;
                  ffz=term_ij*cos_z_ij-dhdcos_ij*dcosdzj;
                  fxa(j)=fxa(j)+ffx;
                  fya(j)=fya(j)+ffy;
                  fza(j)=fza(j)+ffz;
                  fxa(i)=fxa(i)-ffx;
                  fya(i)=fya(i)-ffy;
                  fza(i)=fza(i)-ffz;
                  ffx=term_ik*cos_x_ik-dhdcos_ik*dcosdxk;
                  ffy=term_ik*cos_y_ik-dhdcos_ik*dcosdyk;
                  ffz=term_ik*cos_z_ik-dhdcos_ik*dcosdzk;
                  fxa(k)=fxa(k)+ffx;
                  fya(k)=fya(k)+ffy;
                  fza(k)=fza(k)+ffz;
                  fxa(i)=fxa(i)-ffx;
                  fya(i)=fya(i)-ffy;
                  fza(i)=fza(i)-ffz;
               endif
            end do
         endif
      end do

   end do

END SUBROUTINE SWcalczone


!> subroutine to compute forces and energies using Stillinger-Weber potential
subroutine SWcalcforce(nat,posa,boxl,tmp_force, pot_energy)

   use SWpotential
   use defs, only : boundary,maxnei,iproc

   implicit none

   integer, intent(in)                               :: nat
   real(kind=8), intent(in), dimension(3*nat) :: posa
   real(kind=8), dimension(3), intent(inout)          :: boxl
   real(kind=8), intent(out) :: pot_energy
   real(kind=8), intent(out), dimension(3*nat), target:: tmp_force

   integer, dimension(nat) :: numnei 
   integer, dimension(nat,maxnei) :: nei 
   real(kind=8), dimension(3*nat), target :: pos_normalised
   real(kind=8), dimension(3*nat) :: box_vec
   real(kind=8), dimension(3*nat),target :: pos
   real(kind=8), dimension(3)             :: box
   integer                           :: NATOMS

   real(kind=8), dimension(:), pointer :: xa, ya, za
   real(kind=8), dimension(:), pointer :: fxa, fya, fza

   integer :: i, j, i_id, j_id, ind_j, k, k_id, ind_k
   integer :: my_counter_i
   real(kind=8) :: invsig
   real(kind=8) :: xi, yi, zi, xij, yij, zij, rij, rij2
   real(kind=8) :: twobodyenergy, threebodyenergy
   real(kind=8) :: cos_x_ij, cos_y_ij, cos_z_ij, invrij, rhoij
   real(kind=8) :: one_o_a_ij, expo, gam_o_a_ij, exp_gam_ij, r_to_minusp
   real(kind=8) :: one_o_a2, term1, term2,fact3_ij
   real(kind=8) :: ffx,ffy,ffz, term_ij, term_ik
   real(kind=8) :: dcosdxj,dcosdyj,dcosdzj,dcosdxk,dcosdyk,dcosdzk
   real(kind=8) :: dhdcos_ij,dhdcos_ik,cos_p_1o3,cos_jik
   real(kind=8) :: cos_x_ik, cos_y_ik, cos_z_ik, one_o_a_ik, gam_o_a_ik, fact     
   real(kind=8) :: xik,yik,zik,rik, rik2, invrik,rhoik


   NATOMS = nat
   pos = posa
   box = boxl


   box_vec(1:natoms) = boxl(1)
   box_vec(1+natoms:2*natoms) = boxl(2)
   box_vec(1+natoms+natoms:3*natoms) = boxl(3)

   nei = 0 
   numnei  = 0


   if ( boundary .eq. 'P') pos = modulo(pos, box_vec) !apply PBC implement surface later on
   if ( boundary .eq. 'S') pos(1:natoms) = modulo(pos(1:natoms), box_vec(1:natoms))
   if ( boundary .eq. 'S') pos(2*natoms+1:3*natoms) = modulo(pos(2*natoms+1:3*natoms), box_vec(2*natoms+1:3*natoms))


   ! Generate a neighbour list

   call neighbours(NATOMS,pos,boxl,boundary,maxnei,numnei,nei)

   ! We first set-up pointers for the x, y, z components in the position and forces
   pos_normalised = pos / box_vec

   xa => pos_normalised(1:NATOMS)
   ya => pos_normalised(NATOMS+1:2*NATOMS)
   za => pos_normalised(2*NATOMS+1:3*NATOMS)

   fxa => tmp_force(1:NATOMS)
   fya => tmp_force(NATOMS+1:2*NATOMS)
   fza => tmp_force(2*NATOMS+1:3*NATOMS)

   tmp_force = 0.0d0   ! Vectorial operation

   twobodyenergy = 0.0d0
   threebodyenergy = 0.0d0


   do i=1, NATOMS
      !     i_id = types(i)
      i_id = 1
      xi = xa(i)
      yi = ya(i)
      zi = za(i)
      do ind_j=1, numnei(i)
         j = nei(i,ind_j)
         j_id = 1
         !j_id = types(j)

         ! Pair interactions of i and j
         ! Distance, with periodic boundary conditions
         if (boundary == "P") then
            xij = xa(j) - xi - 1.0d0 * nint( xa(j)-xi )
            yij = ya(j) - yi - 1.0d0 * nint( ya(j)-yi )
            zij = za(j) - zi - 1.0d0 * nint( za(j)-zi )
            elseif (boundary == "S") then
            xij = xa(j) - xi - 1.0d0 * nint( xa(j)-xi )
            yij = ya(j) - yi
            zij = za(j) - zi - 1.0d0 * nint( za(j)-zi )
         endif

         ! Rescale the lengths into Angstroems
         xij = xij * box(1)
         yij = yij * box(2)
         zij = zij * box(3)

         rij2 = xij*xij + yij*yij + zij*zij


         SIGMA = (SW_parameters(i,1)+SW_parameters(j,1))*0.5d0
         ALPHA = (SW_parameters(i,3)+SW_parameters(j,3))*0.5d0


         ! Check the cut-off before proceeding
         if( rij2 < ALPHA*SIGMA*ALPHA*SIGMA ) then

            !initialize the parameters

            A = (SW_parameters(i,2) + SW_parameters(j,2))*0.5d0
            EPSILON = (SW_parameters(i,5)+SW_Parameters(j,5))*0.5d0
            GAMMA = (SW_parameters(i,7)+SW_parameters(j,7))*0.5d0
            P = (SW_parameters(i,6)+SW_parameters(j,6))*0.5d0
            BETA = (SW_parameters(i,4)+SW_parameters(j,4))*0.5d0

            A_EPS = (A * EPSILON)   
            invsig = 1.0D0 / SIGMA

            rij = sqrt(rij2)
            invrij = 1.0d0 / rij
            rhoij = rij * invsig
            cos_x_ij = xij * invrij
            cos_y_ij = yij * invrij
            cos_z_ij = zij * invrij

            ! Some useful quantities
            one_o_a_ij =1.0/(rhoij-ALPHA)
            expo=exp(one_o_a_ij)
            gam_o_a_ij=GAMMA*one_o_a_ij
            exp_gam_ij=exp(gam_o_a_ij)
            r_to_minusp=rhoij ** (-1.0d0*P)

            ! Two body energy and force 
            term1=A_EPS*(BETA*r_to_minusp-1.0d0)*expo
            one_o_a2 = one_o_a_ij * one_o_a_ij
            term2=(one_o_a2*term1+A_EPS*P*BETA*r_to_minusp*expo/rhoij)*invsig

            ! Contribution to the binary repulsive term 
            if(rij <=  S0(i_id,j_id)) then
               term1=term1+A0(i_id,j_id)*(cos(PI*rij/S0(i_id,j_id))+1.0d0)
               term2=term2+A0(i_id,j_id)*PI/S0(i_id,j_id)* sin(PI*rij/S0(i_id,j_id))
            endif
            twobodyenergy=twobodyenergy + 0.5d0*term1;

            fxa(i)=fxa(i)-term2*cos_x_ij;
            fya(i)=fya(i)-term2*cos_y_ij;
            fza(i)=fza(i)-term2*cos_z_ij;

            ! Prepare for the three body term 
            fact3_ij=gam_o_a_ij*one_o_a_ij*invsig

            do ind_k = ind_j+1, numnei(i)
               ! Triplet interaction with i in the middle; all interactions
               k = nei(i,ind_k)
               !k_id = types(k)
               k_id = 1
               ! Distance, with periodic boundary conditions
               if (boundary == "P") then
                  xik = xa(k) - xi - 1.0d0 * nint( xa(k)-xi )
                  yik = ya(k) - yi - 1.0d0 * nint( ya(k)-yi )
                  zik = za(k) - zi - 1.0d0 * nint( za(k)-zi )
                  elseif (boundary == "S") then
                  xik = xa(k) - xi - 1.0d0 * nint( xa(k)-xi )
                  yik = ya(k) - yi
                  zik = za(k) - zi - 1.0d0 * nint( za(k)-zi )
               endif

               ! Rescale the lengths into Angstroems
               xik = xik * box(1)
               yik = yik * box(2)
               zik = zik * box(3)

               rik2 = xik*xik + yik*yik + zik*zik

               SIGMA = (SW_parameters(i,1) + SW_parameters(k,1))*0.5d0
               ALPHA = (SW_parameters(i,3) + SW_parameters(k,3))*0.5d0

               ! Check whether the distance is too large 
               if (rik2<ALPHA*SIGMA*ALPHA*SIGMA)  then

                  !initialize the parameters

                  EPSILON = (SW_parameters(i,5)+SW_Parameters(k,5))*0.5d0
                  GAMMA = (SW_parameters(i,7) + SW_parameters(k,7))*0.5d0
                  LAMBDA = (SW_parameters(i,8))


                  if(angles_to_fit(i) == 0) then                 
                     ONE_THIRD = 1.0d0/3.0d0
                  else 
                     do my_counter_i=1,angles_to_fit(i)
                        if ( jk_of_angle(i,my_counter_i,1) == j .and. jk_of_angle(i,my_counter_i,2) == k ) then
                           !   if (iproc ==0) write(*,*) "looking at angle", optimum_angle(i,my_counter_i) 
                           ONE_THIRD = optimum_angle(i,my_counter_i)    
                        endif
                     enddo

                  endif

                  invsig = 1.0D0 / SIGMA

                  rik=sqrt(rik2)
                  invrik=1.0/rik
                  rhoik=rik*invsig
                  cos_x_ik=xik*invrik
                  cos_y_ik=yik*invrik
                  cos_z_ik=zik*invrik

                  ! Some useful quantities 
                  one_o_a_ik=1.0/(rhoik-ALPHA)
                  gam_o_a_ik=GAMMA*one_o_a_ik
                  fact     =EPSILON*LAMBDA*exp(gam_o_a_ik)*exp_gam_ij

                  cos_jik  =cos_x_ij*cos_x_ik+cos_y_ij*cos_y_ik+ cos_z_ij*cos_z_ik
                  cos_p_1o3=cos_jik+ONE_THIRD

                  ! Energy (added only to central atom)  
                  threebodyenergy=threebodyenergy+fact*cos_p_1o3*cos_p_1o3

                  ! Force 
                  term_ij=fact*fact3_ij*cos_p_1o3*cos_p_1o3;
                  dhdcos_ij=2*fact*cos_p_1o3;
                  term_ik=fact*gam_o_a_ik*one_o_a_ik*cos_p_1o3*cos_p_1o3/SIGMA;
                  dhdcos_ik=2*fact*cos_p_1o3;

                  dcosdxj=(cos_x_ik-cos_jik*cos_x_ij)*invrij;
                  dcosdyj=(cos_y_ik-cos_jik*cos_y_ij)*invrij;
                  dcosdzj=(cos_z_ik-cos_jik*cos_z_ij)*invrij;

                  dcosdxk=(cos_x_ij-cos_jik*cos_x_ik)*invrik;
                  dcosdyk=(cos_y_ij-cos_jik*cos_y_ik)*invrik;
                  dcosdzk=(cos_z_ij-cos_jik*cos_z_ik)*invrik;

                  ffx=term_ij*cos_x_ij-dhdcos_ij*dcosdxj;
                  ffy=term_ij*cos_y_ij-dhdcos_ij*dcosdyj;
                  ffz=term_ij*cos_z_ij-dhdcos_ij*dcosdzj;
                  fxa(j)=fxa(j)+ffx;
                  fya(j)=fya(j)+ffy;
                  fza(j)=fza(j)+ffz;
                  fxa(i)=fxa(i)-ffx;
                  fya(i)=fya(i)-ffy;
                  fza(i)=fza(i)-ffz;
                  ffx=term_ik*cos_x_ik-dhdcos_ik*dcosdxk;
                  ffy=term_ik*cos_y_ik-dhdcos_ik*dcosdyk;
                  ffz=term_ik*cos_z_ik-dhdcos_ik*dcosdzk;
                  fxa(k)=fxa(k)+ffx;
                  fya(k)=fya(k)+ffy;
                  fza(k)=fza(k)+ffz;
                  fxa(i)=fxa(i)-ffx;
                  fya(i)=fya(i)-ffy;
                  fza(i)=fza(i)-ffz;
               endif
            end do
         endif
      end do
   end do
   pot_energy= twobodyenergy+threebodyenergy

END SUBROUTINE SWcalcforce
