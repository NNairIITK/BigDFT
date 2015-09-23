!> @file
!! Define some module for BigDFT+ART
!! @author
!!    Written by Laurent Karim Beland, UdeM 2011!!
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module BAYES used by ART (Activation Relaxation Technique)
module BAYES

   implicit none
   integer :: data_size
   integer, parameter :: configs_to_learn = 10
   integer :: configs_learned = 0

   real(kind=8), dimension(:), allocatable :: covariance_new_data
   real(kind=8), dimension(:,:), allocatable :: covariance_data_data,matrice_inv
   real(kind=8), dimension(:), allocatable :: work
   real(kind=8), dimension(:), allocatable :: data_func_x,data_func_y,data_func_z
   real(kind=8), dimension(:), allocatable :: data_coord_x,data_coord_y,data_coord_z,data_coord_depth

   integer, dimension(:), allocatable :: data_coord_nei
   integer, dimension(:), allocatable :: pivot

   real(kind=8), dimension(:), allocatable :: data_trans_x,data_trans_y,data_trans_z
   real(kind=8) :: delta=0.7d0
   real(kind=8), dimension(3,4) :: sigma=0.2d0
   real(kind=8) :: depth_toler = 1.0d0
   real(kind=8) :: prob_h!probability of hyperparameters being right

   real(kind=8), dimension(:), allocatable :: tmpdata_coord_x,tmpdata_coord_y,tmpdata_coord_z&
      &   ,tmpdata_coord_depth
   integer, dimension(:), allocatable :: tmpdata_coord_nei
   real(kind=8), dimension(:), allocatable :: tmpdata_func_x,tmpdata_func_y,tmpdata_func_z

END MODULE BAYES


!> ART add_bayes_data
subroutine add_bayes_data(pos_in)

   use random
   use bayes
   use defs
   implicit none
   real(kind=8), dimension(3*natoms),intent(in) :: pos_in
   integer :: config_size,j
   integer, dimension(natoms) :: numnei
   integer, dimension(natoms,maxnei) :: nei
   real(kind=8), dimension(3*natoms) :: force_tempo
   real(kind=8) :: trash_energy
   integer :: trash_evalf

   if ( allocated(tmpdata_coord_x)) deallocate(tmpdata_coord_x)
   if ( allocated(tmpdata_coord_y)) deallocate(tmpdata_coord_y)
   if ( allocated(tmpdata_coord_z)) deallocate(tmpdata_coord_z)
   if ( allocated(tmpdata_coord_nei)) deallocate(tmpdata_coord_nei)
   if ( allocated(tmpdata_coord_depth)) deallocate(tmpdata_coord_depth)

   if ( allocated(tmpdata_func_x)) deallocate(tmpdata_func_x)
   if ( allocated(tmpdata_func_y)) deallocate(tmpdata_func_y)
   if ( allocated(tmpdata_func_z)) deallocate(tmpdata_func_z)

   if (.not. allocated(tmpdata_coord_x)) allocate(tmpdata_coord_x(data_size))
   if (.not. allocated(tmpdata_coord_y)) allocate(tmpdata_coord_y(data_size))
   if (.not. allocated(tmpdata_coord_z)) allocate(tmpdata_coord_z(data_size))
   if (.not. allocated(tmpdata_coord_nei)) allocate(tmpdata_coord_nei(data_size))
   if (.not. allocated(tmpdata_coord_depth)) allocate(tmpdata_coord_depth(data_size))

   if (.not. allocated(tmpdata_func_x)) allocate(tmpdata_func_x(data_size))
   if (.not. allocated(tmpdata_func_y)) allocate(tmpdata_func_y(data_size))
   if (.not. allocated(tmpdata_func_z)) allocate(tmpdata_func_z(data_size))

   config_size = (nbr_quantum-nbr_quantum_trash)

   tmpdata_coord_x = data_coord_x
   tmpdata_coord_y = data_coord_y
   tmpdata_coord_z = data_coord_z
   tmpdata_coord_nei = data_coord_nei
   tmpdata_coord_depth = data_coord_depth
   tmpdata_func_x = data_func_x
   tmpdata_func_y = data_func_y
   tmpdata_func_z = data_func_z

   data_size = data_size + config_size

   call init_bayes()
   data_coord_x(1:data_size-config_size) = tmpdata_coord_x
   data_coord_y(1:data_size-config_size)  = tmpdata_coord_y
   data_coord_z(1:data_size-config_size)  = tmpdata_coord_z
   data_coord_nei(1:data_size-config_size)  = tmpdata_coord_nei
   data_coord_depth(1:data_size-config_size) = tmpdata_coord_depth   

   data_func_x(1:data_size-config_size)  = tmpdata_func_x
   data_func_y(1:data_size-config_size)  = tmpdata_func_y
   data_func_z(1:data_size-config_size)  = tmpdata_func_z

   call neighbour_bayes(natoms,pos_in,box,boundary,maxnei,numnei, nei)
   energy_type = "BSW" !we get the quantum data

   call calcforce( NATOMS, pos_in, boxref, force_tempo,&
      &   trash_energy, trash_evalf, .false. )
   do j = 1,config_size
      data_func_x(j+(data_size-config_size) ) = force_tempo(j)
      data_func_y(j+(data_size-config_size) ) = force_tempo(j+natoms)
      data_func_z(j+(data_size-config_size) ) = force_tempo(j+natoms+natoms)
   enddo

   energy_type = "SWP" !we could use any other function for the coordinates
   call calcforce( NATOMS, pos_in, boxref, force_tempo,&
      &   trash_energy, trash_evalf, .false. )

   do j = 1,config_size
      data_coord_x(j+(data_size-config_size)) = force_tempo(j)
      data_coord_y(j+(data_size-config_size)) = force_tempo(j+natoms)
      data_coord_z(j+(data_size-config_size)) = force_tempo(j+natoms+natoms)
      data_coord_nei(j+(data_size-config_size))=numnei(j)
      data_coord_depth(j+(data_size-config_size)) = pos_in(j+natoms+natoms)
   enddo

   energy_type = "BAY"

   !call bayes_learn_from_data(1) !x
   !call bayes_learn_from_data(2) !y
   !call bayes_learn_from_data(3) !z

   configs_learned = configs_learned + 1

END SUBROUTINE add_bayes_data


!> ART create_bayes_data
subroutine create_bayes_data()
   use random
   use bayes
   use defs
   use module_base
   implicit none
   integer, parameter :: configs_for_bayes = 1
   real(kind=8) :: dummy
   real(kind=8) :: ran3
   real(kind=8), dimension(3*natoms) :: pos_tempo,force_tempo
   integer :: i,j,config_size,trash_evalf,ierror
   real(kind=8) :: trash_energy
   integer, dimension(natoms,maxnei) :: nei
   integer, dimension(natoms) :: numnei

   data_size = (nbr_quantum-nbr_quantum_trash)*configs_for_bayes
   config_size = (nbr_quantum-nbr_quantum_trash)

   call init_bayes()

   pos_tempo = pos

   do i = 1,configs_for_bayes
      call neighbour_bayes(natoms,pos,box,boundary,maxnei,numnei, nei)
      energy_type = "BSW" !we get the quantum data

      call calcforce( NATOMS, pos, boxref, force_tempo,&
         &   trash_energy, trash_evalf, .false. )
      do j = 1,config_size
         data_func_x(j+(i-1)*config_size) = force_tempo(j)
         data_func_y(j+(i-1)*config_size) = force_tempo(j+natoms)
         data_func_z(j+(i-1)*config_size) = force_tempo(j+natoms+natoms)
      enddo

      energy_type = "SWP" !we could use any other function for the coordinates
      call calcforce( NATOMS, pos, boxref, force_tempo,&
         &   trash_energy, trash_evalf, .false. )

      do j = 1,config_size
         data_coord_x(j+(i-1)*config_size) = force_tempo(j)
         data_coord_y(j+(i-1)*config_size) = force_tempo(j+natoms)
         data_coord_z(j+(i-1)*config_size) = force_tempo(j+natoms+natoms)
         data_coord_nei(j+(i-1)*config_size)=numnei(j)
         data_coord_depth(j+(i-1)*config_size)=pos(j+natoms+natoms)
      enddo

      energy_type = "BAY"

      do j = 1,3*natoms
         dummy = (0.5d0-ran3())*0.04d0
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
         call MPI_Bcast(dummy,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

         pos(j) = pos_tempo(j) + dummy
      enddo
      energy_type = "BAY"

   enddo

   pos = pos_tempo

   !call bayes_learn_from_data(1) !x
   !call bayes_learn_from_data(2) !y
   !call bayes_learn_from_data(3) !z

   configs_learned = configs_learned + configs_for_bayes

END SUBROUTINE create_bayes_data


!> ART init_bayes
subroutine init_bayes()
   use bayes
   use defs
   implicit none

   !we deallocate the whole thing if needed
   if ( allocated(covariance_new_data)) deallocate(covariance_new_data)
   if ( allocated(data_coord_x)) deallocate(data_coord_x)
   if ( allocated(data_coord_y)) deallocate(data_coord_y)
   if ( allocated(data_coord_z)) deallocate(data_coord_z)
   if ( allocated(data_coord_nei)) deallocate(data_coord_nei)
   if ( allocated(data_coord_depth)) deallocate(data_coord_depth)

   if ( allocated(data_func_x)) deallocate(data_func_x)
   if ( allocated(data_func_y)) deallocate(data_func_y)
   if ( allocated(data_func_z)) deallocate(data_func_z)

   !if ( allocated(data_trans)) deallocate(data_trans)
   if ( allocated(work)) deallocate(work)
   if ( allocated(pivot)) deallocate(pivot)

   if ( allocated(covariance_data_data)) deallocate(covariance_data_data)
   if ( allocated(matrice_inv)) deallocate(matrice_inv)
   if ( allocated(data_trans_x)) deallocate(data_trans_x)
   if ( allocated(data_trans_y)) deallocate(data_trans_y)
   if ( allocated(data_trans_z)) deallocate(data_trans_z)

   ! we allocated all of the working vectors
   if (.not. allocated(covariance_new_data)) allocate(covariance_new_data(data_size))
   if (.not. allocated(data_coord_x)) allocate(data_coord_x(data_size))
   if (.not. allocated(data_coord_y)) allocate(data_coord_y(data_size)) 
   if (.not. allocated(data_coord_z)) allocate(data_coord_z(data_size))
   if (.not. allocated(data_coord_nei)) allocate(data_coord_nei(data_size))
   if (.not. allocated(data_coord_depth)) allocate(data_coord_depth(data_size))

   if (.not. allocated(data_func_x)) allocate(data_func_x(data_size))
   if (.not. allocated(data_func_y)) allocate(data_func_y(data_size))
   if (.not. allocated(data_func_z)) allocate(data_func_z(data_size))

   !if (.not. allocated(data_trans)) allocate(data_trans(data_size))
   if (.not. allocated(work)) allocate(work(data_size))
   if (.not. allocated(pivot)) allocate(pivot(data_size))

   if (.not. allocated(covariance_data_data)) allocate(covariance_data_data(data_size,data_size))
   if (.not. allocated(matrice_inv)) allocate(matrice_inv(data_size,data_size))
   if (.not. allocated(data_trans_x)) allocate(data_trans_x(data_size))
   if (.not. allocated(data_trans_y)) allocate(data_trans_y(data_size))
   if (.not. allocated(data_trans_z)) allocate(data_trans_z(data_size))

END SUBROUTINE init_bayes 


!> ART bayes_learn_from_data
subroutine bayes_learn_from_data(what_direction)
   use bayes
   implicit none
   !real(8), dimension(data_size),intent(in) :: datax_in !the abscysse coordinates of the data
   !real(8), dimension(data_size),intent(in) :: dataf_in !the value of the function of the data
   integer,intent(in) :: what_direction !1 = x, 2= y, 3= z
   integer i,j
   integer :: good
   real(kind=8) :: sqrd_exp,sqrd_total,det

   !we build the kernel
   do i = 1,data_size
      do j = i,data_size
         sqrd_total = 1.0d0
         call calc_kernel(data_coord_x(i),data_coord_x(j),sigma(what_direction,1),sqrd_exp)
         sqrd_total =sqrd_total * sqrd_exp
         call calc_kernel(data_coord_y(i),data_coord_y(j),sigma(what_direction,2),sqrd_exp)
         sqrd_total =sqrd_total * sqrd_exp
         call calc_kernel(data_coord_z(i),data_coord_z(j),sigma(what_direction,3),sqrd_exp)
         sqrd_total =sqrd_total * sqrd_exp
         call calc_kernel(data_coord_depth(i),data_coord_depth(j),sigma(what_direction,4),sqrd_exp)
         sqrd_total =sqrd_total * sqrd_exp*delta*delta
         if (data_coord_nei(i) .ne. data_coord_nei(j)) sqrd_total = 0.0d0
         if (abs(data_coord_depth(i) - data_coord_depth(j)) > depth_toler) sqrd_total = 0.0d0
         if (sqrd_total .lt. 0.01d0) sqrd_total = 0.0d0
         covariance_data_data(j,i) = sqrd_total
         covariance_data_data(i,j) = sqrd_total
      enddo
   enddo

   !we now inverse the matrix
   matrice_inv = covariance_data_data
   call DGETRF(data_size,data_size,matrice_inv,data_size,pivot,good)
   det = 0.0d0
   do i = 1,data_size
      det =det + log(abs(matrice_inv(i,i)))
   enddo

   call DGETRI(data_size,matrice_inv,data_size,pivot,work,data_size,good) 

   if (what_direction == 1) data_trans_x = matmul(data_func_x,matrice_inv)
   if (what_direction == 2) data_trans_y = matmul(data_func_y,matrice_inv)
   if (what_direction == 3) data_trans_z = matmul(data_func_z,matrice_inv)

   if (what_direction == 1) prob_h = exp(-0.5d0*dot_product(data_func_x,data_trans_x))&
      &   *exp(det)**(-0.5d0)
   if (what_direction == 2) prob_h = exp(-0.5d0*dot_product(data_func_y,data_trans_y))&
      &   *exp(det)**(-0.5d0)
   if (what_direction == 3) prob_h = exp(-0.5d0*dot_product(data_func_z,data_trans_z))&
      &   *exp(det)**(-0.5d0)

END SUBROUTINE bayes_learn_from_data


!> ART bayes_calc_one_force
subroutine bayes_calc_one_force(coord_in,coord_nei,what_direction,force_out,variance)
   use bayes
   implicit none

   real(kind=8), intent(in), dimension(4) :: coord_in
   integer,intent(in) :: coord_nei
   integer, intent(in) :: what_direction !1 = x, 2= y, 3= z
   real(kind=8), intent(out) :: force_out
   real(kind=8), intent(out) :: variance

   integer :: i
   real(8) :: sqrd_exp,sqrd_total

   do i = 1,data_size
      sqrd_total = 1.0d0
      call calc_kernel(coord_in(1),data_coord_x(i),sigma(what_direction,1),sqrd_exp)
      sqrd_total =sqrd_total * sqrd_exp
      call calc_kernel(coord_in(2),data_coord_y(i),sigma(what_direction,2),sqrd_exp)
      sqrd_total =sqrd_total * sqrd_exp
      call calc_kernel(coord_in(3),data_coord_z(i),sigma(what_direction,3),sqrd_exp)
      sqrd_total =sqrd_total * sqrd_exp
      call calc_kernel(coord_in(4),data_coord_depth(i),sigma(what_direction,4),sqrd_exp)
      sqrd_total =sqrd_total * sqrd_exp*delta*delta

      if (coord_nei .ne. data_coord_nei(i)) sqrd_total = 0.0d0
      if (abs(coord_in(4) - data_coord_depth(i)) > depth_toler) sqrd_total = 0.0d0
      ! if (iproc == 0) write(*,*) sqrd_total,i    
      covariance_new_data(i) = sqrd_total

   enddo 

   if (what_direction == 1) force_out = dot_product(covariance_new_data,data_trans_x)
   if (what_direction == 2) force_out = dot_product(covariance_new_data,data_trans_y)
   if (what_direction == 3) force_out = dot_product(covariance_new_data,data_trans_z)

   variance = dot_product(matmul(covariance_new_data,matrice_inv),covariance_new_data)
   variance = abs(delta*delta-variance)

END SUBROUTINE bayes_calc_one_force


!> ART calc_kernel
subroutine calc_kernel(x1,x2,sigma,kernel) !this subroutine calculates the kernel
   implicit none
   real(kind=8), intent(in) :: x1,x2,sigma
   real(kind=8), intent(out) :: kernel
   real(kind=8) :: unsursigma

   unsursigma = 1.0d0 / sigma
   kernel = 1.0d0 * exp(-0.5d0*(x1-x2)*(x1-x2)*unsursigma*unsursigma)

END SUBROUTINE calc_kernel


subroutine squared_exponential(x1,x2,SE)
   implicit none
   real(kind=8), intent(in) :: x1,x2
   real(kind=8), intent(out) :: SE

   SE = 4.0d0*exp(-2.5d0*(x1-x2)*(x1-x2))

END SUBROUTINE squared_exponential


!> ART calcforce_bayes
subroutine calcforce_bayes(nat,posa,boxl,forca, energy)

   use defs, only : nbr_quantum,nbr_quantum_trash,maxnei,boundary
   use bayes,only : sigma,delta,prob_h,configs_learned,configs_to_learn
   use random
   implicit none
   integer,      intent(in)                            :: nat
   real(kind=8), intent(in),  dimension(3*nat)         :: posa
   real(kind=8), dimension(3), intent(inout)           :: boxl
   real(kind=8), intent(out), dimension(3*nat)         :: forca
   real(kind=8), intent(out)                           :: energy

   real(kind=8), dimension(3*nat) :: force_tmp,force_sum,variance_tmp,variance_sum
   integer :: i,j,k
   real(kind=8), dimension(4) :: input_bayes
   real(kind=8) :: force_bayes,variancex
   integer, dimension(nat,maxnei) :: nei
   integer, dimension(nat) :: numnei

   real(kind=8) :: ran3,tmp_delta
   real(kind=8), dimension(3,4) :: tmp_sigma
   real(kind=8) :: former_prob_h 
   integer :: nb_integration = 400
   integer :: my_counter

   if (configs_learned < configs_to_learn) call add_bayes_data(posa)

   call SWcalcforce(nat,posa,boxl,forca, energy)
   call neighbour_bayes(nat,posa,boxl,boundary,maxnei,numnei, nei)

   !Replace i by 1 otherwise trouble for i uninitialized (TD, 11/2011)
   input_bayes(1) = forca(1)
   input_bayes(2) = forca(1+nat)
   input_bayes(3) = forca(1+nat+nat)
   input_bayes(4) = posa(1+nat+nat)

   !we do a SW calcforce to get the general coordinates of the atoms we want
   my_counter = 0
   former_prob_h = 0.0d0

   sigma = ran3() * 1.45d0+0.001d0
   force_tmp = 0.0d0
   force_sum = 0.0d0
   variance_tmp = 0.0d0
   variance_sum = 0.0d0

   do j = 1,nb_integration
      do k =1,2
         my_counter = my_counter + 1
         tmp_sigma = sigma
         tmp_delta = delta
         if (k == 1)sigma = ran3()*1.0d0 + 0.001d0
         if (k==2) delta = ran3()*2.0d0 + 0.001d0
         call bayes_learn_from_data(1)

         if (prob_h > former_prob_h) then
            former_prob_h = prob_h
            tmp_sigma = sigma
            tmp_delta = delta
         else
            if (ran3() < prob_h/former_prob_h ) then
               former_prob_h = prob_h
               tmp_sigma = sigma
               tmp_delta = delta
            else
               sigma = tmp_sigma
               delta = tmp_delta
               force_sum = force_sum+force_tmp
               variance_sum = variance_sum + variance_tmp
               cycle
            endif
         endif

         if (my_counter > 100) then
            do i = 1,(nbr_quantum-nbr_quantum_trash)
               call  bayes_calc_one_force(input_bayes,numnei(i),1,force_bayes,variancex)      
               force_sum(i) = force_sum(i)+force_bayes
               variance_sum(i) = variance_sum(i) + variancex
               force_tmp(i) = force_bayes
               variance_tmp(i) = variancex
            enddo 
         endif
      enddo
   enddo
   !force_sum = force_sum/real(nb_integration-10)
   !variance_sum = force_sum/real(nb_integration-10)

   force_tmp = 0.0d0
   variance_tmp = 0.0d0
   my_counter = 0
   former_prob_h = 0.0d0
   do j = 1,nb_integration
      do k =1,2
         my_counter = my_counter + 1
         tmp_sigma = sigma
         tmp_delta = delta
         sigma = ran3()*1.0d0 + 0.001d0
         delta = ran3()*2.0d0 + 0.001d0
         call bayes_learn_from_data(2)

         if (prob_h > former_prob_h) then
            former_prob_h = prob_h
            tmp_sigma = sigma
            tmp_delta = delta
         else
            if (ran3() < prob_h/former_prob_h ) then
               former_prob_h = prob_h
               tmp_sigma = sigma
               tmp_delta = delta
            else
               sigma = tmp_sigma
               delta = tmp_delta
               force_sum = force_sum+force_tmp
               variance_sum = variance_sum + variance_tmp
               cycle
            endif
         endif

         if (my_counter > 100) then
            do i = 1,(nbr_quantum-nbr_quantum_trash)
               call  bayes_calc_one_force(input_bayes,numnei(i),2,force_bayes,variancex)
               force_sum(i+nat) = force_sum(i+nat)+force_bayes
               variance_sum(i+nat) = variance_sum(i+nat) + variancex
               force_tmp(i+nat) = force_bayes
               variance_tmp(i+nat) = variancex
            enddo 
         endif
      enddo
   enddo

   force_tmp = 0.0d0
   variance_tmp = 0.0d0
   my_counter = 0
   former_prob_h = 0.0d0
   do j = 1,nb_integration
      do k =1,2
         my_counter = my_counter + 1
         tmp_sigma = sigma
         tmp_delta = delta
         sigma = ran3()*1.0d0 + 0.001d0
         delta = ran3()*2.0d0+0.001d0
         call bayes_learn_from_data(3)

         if (prob_h > former_prob_h) then
            former_prob_h = prob_h
            tmp_sigma = sigma
            tmp_delta = delta
         else
            if (ran3() < prob_h/former_prob_h ) then
               former_prob_h = prob_h
               tmp_sigma = sigma
               tmp_delta = delta
            else
               sigma = tmp_sigma
               delta = tmp_delta
               force_sum = force_sum+force_tmp
               variance_sum = variance_sum + variance_tmp
               cycle
            endif
         endif

         if (my_counter > 100) then
            do i = 1,(nbr_quantum-nbr_quantum_trash)
               call  bayes_calc_one_force(input_bayes,numnei(i),3,force_bayes,variancex)
               force_sum(i+nat+nat) = force_sum(i+nat+nat)+force_bayes
               variance_sum(i+nat+nat) = variance_sum(i+nat+nat) + variancex
               force_tmp(i+nat+nat) = force_bayes
               variance_tmp(i+nat+nat) = variancex
            enddo
         endif
      enddo
   enddo

   force_sum = force_sum/real(nb_integration-100,kind=8)
   variance_sum = force_sum/real(nb_integration-100,kind=8)

   do i = 1,(nbr_quantum-nbr_quantum_trash)
      forca(i) = force_sum(i)
      forca(i+nat) = force_sum(i+nat)
      forca(i+nat+nat) = force_sum(i+nat+nat)
   enddo

   !if (iproc == 0) write(*,*) "avg_variance",variance_totale

END SUBROUTINE calcforce_bayes
