!> Module defining the parameters for the Stillinger-Weber potential.
module SWpotential
  use module_defs, only: gp
  implicit none

  private

  real(gp) :: fcell !< Multiplying factor to adjust cell lattice.
  real(gp),dimension(:,:),allocatable :: SW_parameters
  integer,dimension(:), allocatable ::    angles_to_fit      !<angles to fit per atom
  integer,dimension(:,:,:),allocatable :: jk_of_angle
  real(gp),dimension(:,:), allocatable ::  optimum_angle
  integer      ::                         len_vect_to_fit
  integer      ::                         configs_to_fit

  ! 1   : real(gp), parameter :: SIGMA = 2.095037426D0
  ! 2   : real(gp), parameter :: A =  7.049556277D0
  ! 3   : real(gp), parameter :: ALPHA = 1.8d0 
  ! 4   : real(gp), parameter :: BETA = 0.60222455844D0
  !  real(gp), parameter :: EPSILON = 2.16823D0* 1.07d0 !  FE Balamane et al. PRB 46, 2250
  ! 5   : real(gp), parameter :: EPSILON = 2.16823D0 ! Original SW
  ! 6  :  integer, parameter :: P = 4
  ! 7  :  real(gp), parameter :: GAMMA = 1.2d0
  ! 8  :  real(gp), parameter :: LAMBDA = (21.0d0 * 1.0d0)
  ! 9 :    real(gp), parameter :: ONE_THIRD = 1.0d0/3.0d0

  real(gp), dimension(2,2) :: A0, S0, rcut2

  real(gp), allocatable,dimension(:,:) ::force_ref_fit  !!reference force for fit
  real(gp), allocatable,dimension(:,:) ::force_work_fit  !!working forces for fit
  real(gp), allocatable,dimension(:,:) :: pos_ref_fit  !!reference position for fit

  public :: SWcalcforce, init_potential_SW, free_potential_SW
  
contains

  !> This subroutine initializes the SW parameters.
  subroutine init_potential_SW(natoms, ntypes, acell, ixc)
    use module_xc
    implicit none
    integer, intent(in) :: natoms, ntypes
    real(gp), intent(in), optional :: acell
    integer, intent(in), optional :: ixc

    real(gp) :: RCUT, EPSILON
    character(len=500) :: name_ixc, name_xcpsp(2)

    !the 10th bond is the reference
    if (.not. allocated(SW_parameters)) allocate(SW_parameters(natoms,8))
    if (.not. allocated(angles_to_fit)) allocate(angles_to_fit(natoms))!angles to fit per atom
    if (.not. allocated(jk_of_angle)) allocate(jk_of_angle(natoms,100,2)) !up to 100 angles to fit
    if (.not. allocated(optimum_angle)) allocate(optimum_angle(natoms,100))

    angles_to_fit = 0
    jk_of_angle = 0
    configs_to_fit = 4

    if (.not. allocated(force_ref_fit)) allocate(force_ref_fit(3*natoms,configs_to_fit))
    if (.not. allocated(pos_ref_fit)) allocate(pos_ref_fit(3*natoms,configs_to_fit))
    if (.not. allocated(force_work_fit)) allocate(force_work_fit(3*natoms,configs_to_fit))

    ! SW_parameters(:,1) =  2.095037426_GP
    ! we change sigma so that cell size is 5.4 and not 5.43 to emulate LDA
    fcell = 1._gp
    if (present(acell)) then
       fcell = acell / 5.43_gp
    else if (present(ixc)) then
       if (ixc < 0) then
          call xc_get_name(name_ixc, ixc, XC_MIXED)
       else
          call xc_get_name(name_ixc, ixc, XC_ABINIT)
       end if
       call xc_get_name(name_xcpsp(1),  1, XC_ABINIT)
       call xc_get_name(name_xcpsp(2), 11, XC_ABINIT)
       if (name_ixc == name_xcpsp(1)) then
          fcell = 5.4_gp / 5.43_gp
       else if (name_ixc == name_xcpsp(2)) then
          fcell = 5.465_gp / 5.43_gp
       end if
    end if
    SW_parameters(:,1) = 2.095037426_gp * fcell

    SW_parameters(:,2) = 7.049556277_GP
    SW_parameters(:,3) = 1.8_gp 
    SW_parameters(:,4) =  0.60222455844_GP
    !  real(gp), parameter :: EPSILON = 2.16823_GP* 1.07_gp !  FE Balamane et al. PRB 46, 2250
    SW_parameters(:,5) = 2.16823_GP ! Original SW
    SW_parameters(:,6) = 4.0_gp
    SW_parameters(:,7) = 1.2_gp
    SW_parameters(:,8) = (21.0_gp * 1.0_gp)
    optimum_angle = 1.0_gp/3.0_gp
    !SW_parameters(:,9:29) = (1.0_gp/3.0_gp)
    RCUT = SW_parameters(1,1) * SW_parameters(1,3)
    EPSILON = SW_parameters(1, 5)

    if (ntypes == 2) then
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
    else if (ntypes == 1 ) then
       S0(1,1) = 3.6_gp
       A0(1,1) = 0.0_gp
       rcut2      = 0.0_gp
       rcut2(1,1) = RCUT*RCUT
    else
       stop
    endif
  END SUBROUTINE init_potential_SW


  subroutine reset_SW_potential()
    implicit none

    !  SW_parameters(:,1) =  2.095037426_GP
    !  we change sigma so that cell size is 5.4 and not 5.43 to emulate LDA
    !  5.465 for GGA
    SW_parameters(:,1) = 2.095037426_gp * fcell
    SW_parameters(:,2) =   7.049556277_GP
    SW_parameters(:,3) = 1.8_gp
    SW_parameters(:,4) =  0.60222455844_GP
    !  real(gp), parameter :: EPSILON = 2.16823_GP* 1.07_gp !  FE Balamane et al. PRB 46, 2250
    SW_parameters(:,5) =  2.16823_GP ! Original SW
    SW_parameters(:,6) =  4.0_gp
    SW_parameters(:,7) =  1.2_gp
    SW_parameters(:,8) =  (21.0_gp * 1.0_gp)
    optimum_angle = 1.0_gp/3.0_gp
    !SW_parameters(:,9:29) = (1.0_gp/3.0_gp)
  END SUBROUTINE reset_SW_potential

  subroutine free_potential_SW()
    implicit none
    
    if (allocated(SW_parameters)) deallocate(SW_parameters)
    if (allocated(angles_to_fit)) deallocate(angles_to_fit)
    if (allocated(jk_of_angle)) deallocate(jk_of_angle)
    if (allocated(optimum_angle)) deallocate(optimum_angle)

    if (allocated(force_ref_fit)) deallocate(force_ref_fit)
    if (allocated(pos_ref_fit)) deallocate(pos_ref_fit)
    if (allocated(force_work_fit)) deallocate(force_work_fit)
  end subroutine free_potential_SW

  !> subroutine to compute forces and energies using Stillinger-Weber potential
  subroutine SWcalcforce(astruct, rxyz, fxyz, pot_energy, infocode)
    use numerics, only: Bohr_Ang, pi, eV_Ha
    use module_defs, only: gp
    use module_atoms, only: atomic_structure, atomic_neighbours, &
         & astruct_neighbours_next, astruct_neighbours_iter, &
         & deallocate_atomic_neighbours
    use dynamic_memory
    use f_utils

    implicit none

    type(atomic_structure), intent(in) :: astruct
    real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
    real(gp), intent(out) :: pot_energy
    real(gp), intent(out), dimension(3, astruct%nat), target:: fxyz
    integer, intent(out) :: infocode

    type(atomic_neighbours) :: nei, nei2

    real(gp), dimension(3*astruct%nat), target :: pos_normalised
    real(gp), dimension(3*astruct%nat) :: box_vec
    real(gp), dimension(:), pointer :: pos
    real(gp), dimension(3) :: box
    integer :: NATOMS
    real(gp)  :: SIGMA 
    real(gp)  :: A 
    real(gp)  :: ALPHA 
    real(gp)  :: BETA
    !  real(gp), parameter :: EPSILON = 2.16823_GP* 1.07_gp !  FE Balamane et al. PRB 46, 2250
    real(gp)  :: EPSILON  ! Original SW
    real(gp)  :: P 
    real(gp)  :: GAMMA
    real(gp)  :: LAMBDA 
    real(gp) :: ONE_THIRD = 1.0_gp/3.0_gp
    real(gp) :: A_EPS

    real(gp), dimension(:), pointer :: xa, ya, za

    integer :: i, j, i_id, j_id, k, k_id
    integer :: my_counter_i
    real(gp) :: invsig
    real(gp) :: xi, yi, zi, xij, yij, zij, rij, rij2
    real(gp) :: twobodyenergy, threebodyenergy
    real(gp) :: cos_x_ij, cos_y_ij, cos_z_ij, invrij, rhoij
    real(gp) :: one_o_a_ij, expo, gam_o_a_ij, exp_gam_ij, r_to_minusp
    real(gp) :: one_o_a2, term1, term2,fact3_ij
    real(gp) :: ffx,ffy,ffz, term_ij, term_ik
    real(gp) :: dcosdxj,dcosdyj,dcosdzj,dcosdxk,dcosdyk,dcosdzk
    real(gp) :: dhdcos_ij,dhdcos_ik,cos_p_1o3,cos_jik
    real(gp) :: cos_x_ik, cos_y_ik, cos_z_ik, one_o_a_ik, gam_o_a_ik, fact     
    real(gp) :: xik,yik,zik,rik, rik2, invrik,rhoik

    infocode = 0

    NATOMS = astruct%nat
    pos = f_malloc_ptr(NATOMS * 3, id = "pos")
    pos(1:NATOMS) = rxyz(1, :) * Bohr_Ang
    pos(1 + NATOMS:2 * NATOMS) = rxyz(2, :) * Bohr_Ang
    pos(1 + 2 * NATOMS:3 * NATOMS) = rxyz(3, :) * Bohr_Ang
    box = astruct%cell_dim * Bohr_Ang
    if (astruct%geocode == "S") box(2) = 1._gp

    box_vec(1:natoms) = box(1)
    box_vec(1+natoms:2*natoms) = box(2)
    box_vec(1+natoms+natoms:3*natoms) = box(3)

    if ( astruct%geocode .eq. 'P') pos = modulo(pos, box_vec) !apply PBC implement surface later on
    if ( astruct%geocode .eq. 'S') pos(1:natoms) = modulo(pos(1:natoms), box_vec(1:natoms))
    if ( astruct%geocode .eq. 'S') pos(2*natoms+1:3*natoms) = modulo(pos(2*natoms+1:3*natoms), box_vec(2*natoms+1:3*natoms))

    ! Generate a neighbour list
    call astruct_neighbours(astruct, rxyz, nei)

    ! We first set-up pointers for the x, y, z components in the position and forces
    pos_normalised = pos / box_vec

    xa => pos_normalised(1:NATOMS)
    ya => pos_normalised(NATOMS+1:2*NATOMS)
    za => pos_normalised(2*NATOMS+1:3*NATOMS)

    call f_zero(fxyz)

    twobodyenergy = 0.0_gp
    threebodyenergy = 0.0_gp


    do i=1, NATOMS
       !     i_id = types(i)
       i_id = 1
       xi = xa(i)
       yi = ya(i)
       zi = za(i)

       call astruct_neighbours_iter(nei, i)
       do while(astruct_neighbours_next(nei, j))
          j_id = 1
          !j_id = types(j)

          ! Pair interactions of i and j
          ! Distance, with periodic boundary conditions
          if (astruct%geocode == "P") then
             xij = xa(j) - xi - 1.0_gp * nint( xa(j)-xi )
             yij = ya(j) - yi - 1.0_gp * nint( ya(j)-yi )
             zij = za(j) - zi - 1.0_gp * nint( za(j)-zi )
          elseif (astruct%geocode == "S") then
             xij = xa(j) - xi - 1.0_gp * nint( xa(j)-xi )
             yij = ya(j) - yi
             zij = za(j) - zi - 1.0_gp * nint( za(j)-zi )
          endif

          ! Rescale the lengths into Angstroems
          xij = xij * box(1)
          yij = yij * box(2)
          zij = zij * box(3)

          rij2 = xij*xij + yij*yij + zij*zij


          SIGMA = (SW_parameters(i,1)+SW_parameters(j,1))*0.5_gp
          ALPHA = (SW_parameters(i,3)+SW_parameters(j,3))*0.5_gp


          ! Check the cut-off before proceeding
          if( rij2 < ALPHA*SIGMA*ALPHA*SIGMA ) then

             !initialize the parameters

             A = (SW_parameters(i,2) + SW_parameters(j,2))*0.5_gp
             EPSILON = (SW_parameters(i,5)+SW_Parameters(j,5))*0.5_gp
             GAMMA = (SW_parameters(i,7)+SW_parameters(j,7))*0.5_gp
             P = (SW_parameters(i,6)+SW_parameters(j,6))*0.5_gp
             BETA = (SW_parameters(i,4)+SW_parameters(j,4))*0.5_gp

             A_EPS = (A * EPSILON)   
             invsig = 1.0_GP / SIGMA

             rij = sqrt(rij2)
             invrij = 1.0_gp / rij
             rhoij = rij * invsig
             cos_x_ij = xij * invrij
             cos_y_ij = yij * invrij
             cos_z_ij = zij * invrij

             ! Some useful quantities
             one_o_a_ij =1.0/(rhoij-ALPHA)
             expo=exp(one_o_a_ij)
             gam_o_a_ij=GAMMA*one_o_a_ij
             exp_gam_ij=exp(gam_o_a_ij)
             r_to_minusp=rhoij ** (-1.0_gp*P)

             ! Two body energy and force 
             term1=A_EPS*(BETA*r_to_minusp-1.0_gp)*expo
             one_o_a2 = one_o_a_ij * one_o_a_ij
             term2=(one_o_a2*term1+A_EPS*P*BETA*r_to_minusp*expo/rhoij)*invsig

             ! Contribution to the binary repulsive term 
             if(rij <=  S0(i_id,j_id)) then
                term1=term1+A0(i_id,j_id)*(cos(pi*rij/S0(i_id,j_id))+1.0_gp)
                term2=term2+A0(i_id,j_id)*pi/S0(i_id,j_id)* sin(pi*rij/S0(i_id,j_id))
             endif
             twobodyenergy=twobodyenergy + 0.5_gp*term1;

             fxyz(1, i) = fxyz(1, i)-term2*cos_x_ij;
             fxyz(2, i) = fxyz(2, i)-term2*cos_y_ij;
             fxyz(3, i) = fxyz(3, i)-term2*cos_z_ij;

             ! Prepare for the three body term 
             fact3_ij=gam_o_a_ij*one_o_a_ij*invsig

             nei2 = nei
             do while (astruct_neighbours_next(nei2, k))
                ! Triplet interaction with i in the middle; all interactions
                !k_id = types(k)
                k_id = 1
                ! Distance, with periodic boundary conditions
                if (astruct%geocode == "P") then
                   xik = xa(k) - xi - 1.0_gp * nint( xa(k)-xi )
                   yik = ya(k) - yi - 1.0_gp * nint( ya(k)-yi )
                   zik = za(k) - zi - 1.0_gp * nint( za(k)-zi )
                elseif (astruct%geocode == "S") then
                   xik = xa(k) - xi - 1.0_gp * nint( xa(k)-xi )
                   yik = ya(k) - yi
                   zik = za(k) - zi - 1.0_gp * nint( za(k)-zi )
                endif

                ! Rescale the lengths into Angstroems
                xik = xik * box(1)
                yik = yik * box(2)
                zik = zik * box(3)

                rik2 = xik*xik + yik*yik + zik*zik

                SIGMA = (SW_parameters(i,1) + SW_parameters(k,1))*0.5_gp
                ALPHA = (SW_parameters(i,3) + SW_parameters(k,3))*0.5_gp

                ! Check whether the distance is too large 
                if (rik2<ALPHA*SIGMA*ALPHA*SIGMA)  then

                   !initialize the parameters

                   EPSILON = (SW_parameters(i,5)+SW_Parameters(k,5))*0.5_gp
                   GAMMA = (SW_parameters(i,7) + SW_parameters(k,7))*0.5_gp
                   LAMBDA = (SW_parameters(i,8))


                   if(angles_to_fit(i) == 0) then                 
                      ONE_THIRD = 1.0_gp/3.0_gp
                   else 
                      do my_counter_i=1,angles_to_fit(i)
                         if ( jk_of_angle(i,my_counter_i,1) == j .and. jk_of_angle(i,my_counter_i,2) == k ) then
                            !   if (iproc ==0) write(*,*) "looking at angle", optimum_angle(i,my_counter_i) 
                            ONE_THIRD = optimum_angle(i,my_counter_i)    
                         endif
                      enddo

                   endif

                   invsig = 1.0_GP / SIGMA

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
                   fxyz(1, j) = fxyz(1, j)+ffx;
                   fxyz(2, j) = fxyz(2, j)+ffy;
                   fxyz(3, j) = fxyz(3, j)+ffz;
                   fxyz(1, i) = fxyz(1, i)-ffx;
                   fxyz(2, i) = fxyz(2, i)-ffy;
                   fxyz(3, i) = fxyz(3, i)-ffz;
                   ffx=term_ik*cos_x_ik-dhdcos_ik*dcosdxk;
                   ffy=term_ik*cos_y_ik-dhdcos_ik*dcosdyk;
                   ffz=term_ik*cos_z_ik-dhdcos_ik*dcosdzk;
                   fxyz(1, k) = fxyz(1, k)+ffx;
                   fxyz(2, k) = fxyz(2, k)+ffy;
                   fxyz(3, k) = fxyz(3, k)+ffz;
                   fxyz(1, i) = fxyz(1, i)-ffx;
                   fxyz(2, i) = fxyz(2, i)-ffy;
                   fxyz(3, i) = fxyz(3, i)-ffz;
                endif
             end do
          endif
       end do
    end do
    pot_energy = (twobodyenergy+threebodyenergy) * eV_Ha
    fxyz = fxyz / Bohr_Ang * eV_Ha

    call f_free_ptr(pos)
    call deallocate_atomic_neighbours(nei)
  END SUBROUTINE SWcalcforce

END MODULE  SWpotential
