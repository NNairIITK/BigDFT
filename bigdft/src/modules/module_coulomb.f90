module module_coulomb
  
	implicit none
	
	type, public :: coulomb_force_parameters
		integer :: force_count = 0
		real(8) :: gauss_width = 1.d0, q_amp = 4.d0
		real(8), allocatable, dimension(:) :: charges
	end type coulomb_force_parameters

	public :: set_coulomb_gauss_width
	public :: get_coulomb_gauss_width
	public :: alloc_dim_charges_arr
	public :: distribute_random_charges
	public :: set_q_amp
	public :: get_q_amp
	public :: reset_coulomb_force_count
	public :: get_coulomb_force_count
	public :: coulomb_energyandforces

	type(coulomb_force_parameters), save, private :: cfp

	contains

	subroutine set_coulomb_gauss_width(gauss_width)
		real(8) :: gauss_width
		cfp%gauss_width = gauss_width
	end subroutine set_coulomb_gauss_width
  
	function get_coulomb_gauss_width() result(output)
	  real(8) :: output
		output = cfp%gauss_width
	end function get_coulomb_gauss_width
	
	subroutine alloc_dim_charges_arr(nat)
	  integer :: nat
	  allocate(cfp%charges(nat))
	end subroutine alloc_dim_charges_arr
	
	subroutine distribute_charges(nat, charge)
	  integer :: nat
	  real(8), dimension(nat) :: charge
	  cfp%charges = charge
	end subroutine distribute_charges

	subroutine distribute_random_charges(nat)
		!Allocates charges between +q_amp and -q_amp on atoms randomly

		implicit none
		integer :: i, nat
		real(8) :: rnd, q_av
		
		q_av = 0.d0	
		do i = 1, nat
		  call random_number(rnd)
		  
		  cfp%charges(i) = -cfp%q_amp + 2 * cfp%q_amp * rnd
		  q_av = q_av + cfp%charges(i)
		end do

		q_av = q_av / nat      
		do i = 1, nat
		  cfp%charges(i) = cfp%charges(i) - q_av
		end do
		
	end subroutine distribute_random_charges

	subroutine set_q_amp(q_amp)
		real(8) :: q_amp
		cfp%q_amp = q_amp
	end subroutine set_q_amp
  
	function get_q_amp() result(output)
	  real(8) :: output
		output = cfp%q_amp
	end function get_q_amp

	subroutine reset_coulomb_force_count()
		cfp%force_count = 0
	end subroutine reset_coulomb_force_count

	function get_coulomb_force_count() result(output)
	  integer :: output
		output = cfp%force_count
	end function get_coulomb_force_count
  
	subroutine coulomb_energyandforces(nat, rxyz, fxyz, energy)
		implicit none
		integer :: i, j, nat
		real(8) :: rij, rij_x, rij_y, rij_z, ggamma, coulomb_force, energy
		real(8), dimension(:, :) :: rxyz, fxyz
		real(8), parameter :: pi = 4.d0*atan(1.d0)
	 
		
		!dimension cfp%charge(nat)
		
		!allocate (rxyz(3, nat), fxyz(3, nat))

		!gamma = gamma_ij = gamma_ii   in our case
		ggamma = 1.d0 / (dsqrt(2.d0) * cfp%gauss_width)

		cfp%force_count = cfp%force_count + 1

		do i = 1, nat-1
			do j = i+1, nat
				rij_x = rxyz(1, i) - rxyz(1, j)
				rij_y = rxyz(2, i) - rxyz(2, j)
				rij_z = rxyz(3, i) - rxyz(3, j)
				rij = dsqrt(rij_x**2 + rij_y**2 + rij_z**2)

				energy = energy + cfp%charges(i) * cfp%charges(j) * derf(ggamma * rij) / rij

				coulomb_force = cfp%charges(i) * cfp%charges(j) * (derf(ggamma * rij) * rij_x / (rij**3) - &
								(2.d0 / dsqrt (pi)) * (dexp(-(ggamma * rij)**2) * ggamma * rij_x) / (rij**2))
				fxyz(1, i) = fxyz(1, i) + coulomb_force
				fxyz(1, j) = fxyz(1, j) - coulomb_force

				coulomb_force = cfp%charges(i) * cfp%charges(j) * (derf(ggamma * rij) * rij_y / (rij**3) - &
								(2.d0 / dsqrt (pi)) * (dexp(-(ggamma * rij)**2) * ggamma * rij_y) / (rij**2))
				fxyz(2, i) = fxyz(2, i) + coulomb_force
				fxyz(2, j) = fxyz(2, j) - coulomb_force

				coulomb_force = cfp%charges(i) * cfp%charges(j) * (derf(ggamma * rij) * rij_z / (rij**3) - &
								(2.d0 / dsqrt (pi)) * (dexp(-(ggamma * rij)**2) * ggamma * rij_z) / (rij**2))
				fxyz(3, i) = fxyz(3, i) + coulomb_force
				fxyz(3, j) = fxyz(3, j) - coulomb_force
			end do
		end do
	end subroutine coulomb_energyandforces

end module module_coulomb
