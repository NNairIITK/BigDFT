module multipole
  use module_base
  use multipole_base, only: external_potential_descriptors, lmax
  implicit none

  private

  !> Public routines
  public :: interaction_multipoles_ions
  public :: potential_from_charge_multipoles
  public :: ionic_energy_of_external_charges
  public :: gaussian_density
  public :: support_function_gross_multipoles
  public :: calculate_dipole_moment
  public :: calculate_rpowerx_matrices
  public :: multipole_analysis_driver_new

  contains


    !> Calculate the interaction between the ions and the external multipoles.
    !! At the moment only the monopoles are taken into account.
    subroutine interaction_multipoles_ions(iproc, ep, at, eion, fion)
      use module_types, only: atoms_data
      use yaml_output, only: yaml_map
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc
      type(external_potential_descriptors),intent(in) :: ep
      type(atoms_data),intent(in) :: at
      real(gp),intent(inout) :: eion
      real(gp),dimension(3,at%astruct%nat),intent(inout) :: fion

      ! Local variables
      integer :: iat, ityp, impl
      real(gp) :: r, charge, emp

      !write(*,*) 'WARNING DEBUG HERE!!!!!!!!!!!!!!!!!!!!!!!!!'
      !return

      call f_routine(id='interaction_multipoles_ions')

      emp = 0.0_gp
      do iat=1,at%astruct%nat
          ityp=at%astruct%iatype(iat)
          do impl=1,ep%nmpl
              r = sqrt((at%astruct%rxyz(1,iat)-ep%mpl(impl)%rxyz(1))**2 + &
                       (at%astruct%rxyz(2,iat)-ep%mpl(impl)%rxyz(2))**2 + &
                       (at%astruct%rxyz(3,iat)-ep%mpl(impl)%rxyz(3))**2)
              if (associated(ep%mpl(impl)%qlm(0)%q)) then
                  ! For the multipoles, a positive value corresponds to a
                  ! negative charge! Therefore multiply by -1
                  charge = real(at%nelpsp(ityp),gp)*real(-1.0_gp*ep%mpl(impl)%qlm(0)%q(1),kind=gp)
                  emp = emp + charge/r
                  fion(1,iat) = fion(1,iat) + charge/(r**3)*(at%astruct%rxyz(1,iat)-ep%mpl(impl)%rxyz(1))
                  fion(2,iat) = fion(2,iat) + charge/(r**3)*(at%astruct%rxyz(2,iat)-ep%mpl(impl)%rxyz(2))
                  fion(3,iat) = fion(3,iat) + charge/(r**3)*(at%astruct%rxyz(3,iat)-ep%mpl(impl)%rxyz(3))
              end if
          end do
      end do


      if (iproc==0) then
          call yaml_map('Interaction energy ions multipoles',emp)
      end if
      eion = eion + emp


      call f_release_routine()

    end subroutine interaction_multipoles_ions


    !> Calculate the interaction between the external multipoles.
    !! At the moment only the monopoles are taken into account.
    subroutine ionic_energy_of_external_charges(iproc, ep, at, eion)
      use module_types, only: atoms_data
      use yaml_output, only: yaml_map
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc
      type(external_potential_descriptors),intent(in) :: ep
      type(atoms_data),intent(in) :: at
      real(gp),intent(inout) :: eion

      ! Local variables
      integer :: impl, jmpl
      real(gp) :: r, charge, emp

      !write(*,*) 'WARNING DEBUG HERE!!!!!!!!!!!!!!!!!!!!!!!!!'
      !return

      call f_routine(id='ionic_energy_of_external_charges')

      emp = 0.0_gp
      do impl=1,ep%nmpl
          do jmpl=impl+1,ep%nmpl
              r = sqrt((ep%mpl(impl)%rxyz(1)-ep%mpl(jmpl)%rxyz(1))**2 + &
                       (ep%mpl(impl)%rxyz(2)-ep%mpl(jmpl)%rxyz(2))**2 + &
                       (ep%mpl(impl)%rxyz(3)-ep%mpl(jmpl)%rxyz(3))**2)
              if (associated(ep%mpl(impl)%qlm(0)%q)) then
                  ! For the multipoles, a positive value corresponds to a
                  ! negative charge, therefore multiply by -1. Actually it doesn't matter
                  charge = real(-1.0_gp*ep%mpl(impl)%qlm(0)%q(1),kind=gp)*real(-1.0_gp*ep%mpl(jmpl)%qlm(0)%q(1),kind=gp)
                  emp = emp + charge/r
              end if
          end do
      end do

      if (iproc==0) then
          call yaml_map('Interaction energy multipoles multipoles',emp)
      end if
      eion = eion + emp

      call f_release_routine()

    end subroutine ionic_energy_of_external_charges


    !> Calculate the external potential arising from the multipoles of the charge density
    subroutine potential_from_charge_multipoles(iproc, nproc, at, denspot, ep, &
               is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, shift, &
               verbosity, ixc, lzd, pot, rxyz, ixyz0, write_directory, dipole_total, quadrupole_total, all_norms_ok, &
               rho_mp, pot_mp)
      use module_types, only: DFT_local_fields, local_zone_descriptors
      use Poisson_Solver, except_dp => dp, except_gp => gp
      use module_atoms, only: atoms_data
      use bounds, only: ext_buffers
      use yaml_output
      use io, only: plot_density
      use bounds, only: geocode_buffers
      use box, only: cell_periodic_dims
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, verbosity, ixc
      type(atoms_data),intent(in) :: at
      type(DFT_local_fields),intent(inout) :: denspot
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3
      real(gp),intent(in) :: hx, hy, hz
      real(gp),dimension(3),intent(in) :: shift !< global shift of the atomic positions
      type(local_zone_descriptors),intent(in) :: lzd
      real(gp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot
      real(kind=8),dimension(3,at%astruct%nat),intent(in),optional :: rxyz
      integer,dimension(3),intent(in),optional :: ixyz0
      character(len=*),intent(in),optional :: write_directory
      real(kind=8),dimension(3),intent(out),optional :: dipole_total
      real(kind=8),dimension(3,3),intent(out),optional :: quadrupole_total
      logical,intent(out),optional :: all_norms_ok
      real(kind=8),dimension(is1:ie1,is2:ie2,is3:ie3),intent(out),optional :: rho_mp, pot_mp

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l, m, ii, mm, nthread, ithread, ll
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp, ehart_ps, tt, ttt, gg, hhh, tt0, tt1, tt2
      real(dp),dimension(3) :: r
      logical, dimension(3) :: peri
      real(kind=8) :: dr
      real(dp),dimension(:,:,:),allocatable :: density, density_cores
      real(dp),dimension(:,:,:,:),allocatable :: density_loc, potential_loc
      real(kind=8),dimension(0:lmax) :: sigma
      real(8),dimension(:),allocatable :: monopole
      real(8),dimension(:,:),allocatable :: norm, dipole, quadrupole, norm_check
      real(kind=8),dimension(:,:,:),allocatable :: gaussians1, gaussians2, gaussians3
      logical,dimension(:),allocatable :: norm_ok
      real(kind=8),parameter :: norm_threshold = 1.d-2
      real(kind=8),dimension(0:lmax) :: max_error
      integer :: ixc_tmp, nzatom, npspcode, ilr, j1s, j1e, j2s, j2e, j3s, j3e, j1, j2, j3
      integer :: nbl1, nbl2, nbl3, nbr1, nbr2, nbr3, n3pi, i3s, lmax_avail, nl1, nl2, nl3
      integer,dimension(:),allocatable :: nelpsp, psp_source
      real(gp),dimension(0:4,0:6) :: psppar
      logical :: exists, found, found_non_associated, written
      logical :: perx, pery, perz
      logical,parameter :: use_iterator = .false.
      real(kind=8) :: cutoff, rholeaked, hxh, hyh, hzh, rx, ry, rz, qq, ttl, sig
      real(kind=8),dimension(3) :: center
      integer :: n1i, n2i, n3i, itype, ntype
      integer :: nmpx, nmpy, nmpz, ndensity, izion, ioffset, ishift, iat
      real(dp), dimension(:), allocatable  :: mpx,mpy,mpz
      real(kind=8),dimension(:),allocatable :: rmax
      !real(kind=8),parameter :: rmin=3.d-1
      real(kind=8) :: rmin
      real(kind=8),dimension(:),allocatable :: rloc
      integer,dimension(0:lmax) :: error_meaningful
      character(len=20),dimension(0:lmax) :: output_arr
      real(kind=8),dimension(:,:),allocatable :: rxyz_noshift
      integer,dimension(3) :: ixyz0_
      character(len=128) :: filename
      !$ integer  :: omp_get_thread_num,omp_get_max_threads

      call f_routine(id='potential_from_charge_multipoles')

      call f_zero(rholeaked)

      ! Conditions for periodicity
      perx=(at%astruct%geocode /= 'F')
      pery=(at%astruct%geocode == 'P')
      perz=(at%astruct%geocode /= 'F')
      if (perx) then
          j1s = -1
          j1e = 1
      else
          j1s = 0
          j1e = 0
      end if
      if (pery) then
          j2s = -1
          j2e = 1
      else
          j2s = 0
          j2e = 0
      end if
      if (perz) then
          j3s = -1
          j3e = 1
      else
          j3s = 0
          j3e = 0
      end if

      ! No need to do this if there are no multipoles given.
      multipoles_if: if (ep%nmpl>0) then

          hhh = hx*hy*hz

          ! Used for the calculations of the solid harmonics, see description there
          rmin = 2.0d0*hhh**(1.d0/3.d0)
    
          sigma(0) = 1.0d0
          sigma(1) = 0.8d0
          sigma(2) = 0.6d0
    
          density = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='density')
          density_cores = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='density_cores')
          
          nthread = 1
          !$ nthread = omp_get_max_threads()
          density_loc = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3,0.to.nthread-1/),id='density_loc')
          potential_loc = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3,0.to.nthread-1/),id='potential_loc')
    
          gaussians1 = f_malloc((/0.to.lmax,is1.to.ie1,1.to.ep%nmpl/),id='gaussians1')
          gaussians2 = f_malloc((/0.to.lmax,is2.to.ie2,1.to.ep%nmpl/),id='gaussians2')
          gaussians3 = f_malloc((/0.to.lmax,is3.to.ie3,1.to.ep%nmpl/),id='gaussians3')

          do ilr=1,lzd%nlr 
              if (lzd%Llr(ilr)%geocode/='F') then
                  call f_err_throw('support function locregs must always have free BC')
              end if
          end do
          call geocode_buffers('F', lzd%glr%geocode, nl1, nl2, nl3)
          call calculate_gaussian(is1, ie1, 1, nl1, lzd%glr%d%n1i, perx, hx, shift, ep, gaussians1)
          call calculate_gaussian(is2, ie2, 2, nl2, lzd%glr%d%n2i, pery, hy, shift, ep, gaussians2)
          call calculate_gaussian(is3, ie3, 3, nl3, lzd%glr%d%n3i, perz, hz, shift, ep, gaussians3)
    
    
          norm = f_malloc((/0.to.2,1.to.ep%nmpl/),id='norm')
          norm_check = f_malloc((/0.to.2,1.to.ep%nmpl/),id='norm_check')
          monopole = f_malloc(ep%nmpl,id='monopole')
          dipole = f_malloc((/3,ep%nmpl/),id='dipole')
          quadrupole = f_malloc((/5,ep%nmpl/),id='quadrupole')
          norm_ok = f_malloc0(ep%nmpl,id='norm_ok')


    
          ! First calculate the norm of the Gaussians for each multipole
          !norm = 0.d0
          call calculate_norm(nproc, is1, ie1, is2, ie2, is3, ie3, ep, &
               hhh, gaussians1, gaussians2, gaussians3, norm)
    
          ! Check whether they are ok.
          do impl=1,ep%nmpl
              norm_ok(impl) = .true.
              do l=0,lmax !
                  !write(*,*) 'impl, l, norm', impl, l, norm(l,impl)
                  if (abs(1.d0-norm(l,impl))>norm_threshold) then
                      norm_ok(impl) = .false.
                  end if
              end do
              !write(*,*) 'impl, norm_ok(impl)', impl, norm_ok(impl)
          end do

          if (present(all_norms_ok)) then
              all_norms_ok = all(norm_ok)
              if (.not.all_norms_ok) then
                  call f_free(density)
                  call f_free(density_cores)
                  call f_free(density_loc)
                  call f_free(potential_loc)
                  call f_free(gaussians1)
                  call f_free(gaussians2)
                  call f_free(gaussians3)
                  call f_free(norm)
                  call f_free(norm_check)
                  call f_free(monopole)
                  call f_free(dipole)
                  call f_free(quadrupole)
                  call f_free(norm_ok)
                  return
              end if
          end if
    
    
          ! Get the parameters for each multipole, required to compensate for the pseudopotential part
          !nzatom = f_malloc(ep%nmpl,id='nzatom')
          !nelpsp = f_malloc(ep%nmpl,id='nelpsp')
          rloc = f_malloc(ep%nmpl,id='rloc')
!!$         perx = (denspot%dpbox%geocode /= 'F')
!!$         pery = (denspot%dpbox%geocode == 'P')
!!$         perz = (denspot%dpbox%geocode /= 'F')
         peri=cell_periodic_dims(denspot%dpbox%mesh)
         perx =peri(1)
         pery =peri(2)
         perz =peri(3)

         n3pi = denspot%dpbox%n3pi
         i3s = denspot%dpbox%i3s + denspot%dpbox%i3xcsh
         hxh = denspot%dpbox%mesh%hgrids(1)
         hyh = denspot%dpbox%mesh%hgrids(2)
         hzh = denspot%dpbox%mesh%hgrids(3)
         n1i = denspot%dpbox%mesh%ndims(1)
         n2i = denspot%dpbox%mesh%ndims(2)
         n3i = denspot%dpbox%mesh%ndims(3)
         call ext_buffers(perx,nbl1,nbr1)
         call ext_buffers(pery,nbl2,nbr2)
         call ext_buffers(perz,nbl3,nbr3)
         !write(*,*) 'ep%nmpl, n1i, n2i, n3i', ep%nmpl, n1i, n2i, n3i
    
    
    
    
         ! Generate the density that comes from the pseudopotential atoms
         ndensity = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
         !psp_source = f_malloc(ep%nmpl,id='psp_source')
         !do impl=1,ep%nmpl
         !    ! Search the rloc and zion of the corresponding pseudopotential
         !    call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp(impl), psp_source(impl), rloc(impl))
         !end do

         !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
         !cutoff=10.0_gp*maxval(rloc(:))
         cutoff=10.0_gp*maxval(ep%mpl(:)%sigma(0))
         if (at%multipole_preserving) then
            !We want to have a good accuracy of the last point rloc*10
            cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
         end if
         !Separable function: do 1-D integrals before and store it.
         nmpx = (ceiling(cutoff/hxh) - floor(-cutoff/hxh)) + 1
         nmpy = (ceiling(cutoff/hyh) - floor(-cutoff/hyh)) + 1
         nmpz = (ceiling(cutoff/hzh) - floor(-cutoff/hzh)) + 1
         mpx = f_malloc( (/ 0 .to. nmpx /),id='mpx')
         mpy = f_malloc( (/ 0 .to. nmpy /),id='mpy')
         mpz = f_malloc( (/ 0 .to. nmpz /),id='mpz')

         do impl=1,ep%nmpl
             !! Search the rloc and zion of the corresponding pseudopotential
             !call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp(impl), psp_source(impl), rloc(impl))
             if(norm_ok(impl)) then
                 ! The following routine needs the shifted positions
                 rx = ep%mpl(impl)%rxyz(1) - shift(1)
                 ry = ep%mpl(impl)%rxyz(2) - shift(2)
                 rz = ep%mpl(impl)%rxyz(3) - shift(3)
                 call gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, &
                      rx, ry, rz, &
                      ep%mpl(impl)%sigma(0), ep%mpl(impl)%nzion, at%multipole_preserving, use_iterator, at%mp_isf, &
                      denspot%dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, ndensity, density_cores, rholeaked)
                 !!call gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, &
                 !!     rx, ry, rz, &
                 !!     ep%mpl(impl)%sigma(0), nelpsp(impl), at%multipole_preserving, use_iterator, at%mp_isf, &
                 !!     denspot%dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, ndensity, density_cores, rholeaked)
                 !!call gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, &
                 !!     rx, ry, rz, &
                 !!     rloc(impl), nelpsp(impl), at%multipole_preserving, use_iterator, at%mp_isf, &
                 !!     denspot%dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, ndensity, density_cores, rholeaked)
             end if
         end do
!! UNCOMMENT FOR TEST          do i3=is3,ie3
!! UNCOMMENT FOR TEST              do i2=is2,ie2
!! UNCOMMENT FOR TEST                  do i1=is1,ie1
!! UNCOMMENT FOR TEST                     !write(400+iproc,'(a,3i7,es18.6)') 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
!! UNCOMMENT FOR TEST                     write(500+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
!! UNCOMMENT FOR TEST                     tt = tt + density(i1,i2,i3)*hhh
!! UNCOMMENT FOR TEST                 end do
!! UNCOMMENT FOR TEST             end do
!! UNCOMMENT FOR TEST         end do
    
         call f_free(mpx)
         call f_free(mpy)
         call f_free(mpz)
    
    
          ! Calculate the density only within a sphere of radius rmax
          rmax = f_malloc(ep%nmpl,id='rmax')
          do impl=1,ep%nmpl
              !rmax(impl) = min(denspot%dpbox%mesh%ndims(1)*0.25d0*hx, &
              !                 denspot%dpbox%mesh%ndims(2)*0.25d0*hy, &
              !                 denspot%dpbox%mesh%ndims(3)*0.25d0*hz)
              rmax(impl) = min((denspot%dpbox%mesh%ndims(1)-31)*0.25d0*hx, &
                               (denspot%dpbox%mesh%ndims(2)-31)*0.25d0*hy, &
                               (denspot%dpbox%mesh%ndims(3)-31)*0.25d0*hz)
          end do
    
    
          norm_check = 0.d0
          monopole = 0.d0
          dipole = 0.d0
          quadrupole = 0.d0
          !$omp parallel &
          !$omp default(none) &
          !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, hhh, ep, shift, nthread, norm_ok) &
          !$omp shared(norm_check, monopole, dipole, quadrupole, density, density_loc, potential_loc) &
          !$omp shared (gaussians1, gaussians2, gaussians3, rmax, rmin) &
          !$omp shared (j1s, j1e, j2s, j2e, j3s, j3e, nl1, nl2, nl3, lzd) &
          !$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, l, gg, m, mm, tt, ttt, ttl, ithread, center, ll) &
          !$omp private(rnrm1, rnrm2, rnrm3, rnrm5, qq, ii, sig, lmax_avail, found_non_associated, j1, j2, j3, dr)
          ithread = 0
          !$ ithread = omp_get_thread_num()
          if (ithread<0 .or. ithread>nthread-1) then
              !SM: Is it possible to call f_err_throw within OpenMP? Anyway this condition should never be true...
              call f_err_throw('wrong value of ithread',err_name='BIGDFT_RUNTIME_ERROR')
              !LG: yes it is possible but not advised to, as running conditions might arise. we will not be sure of
              !! the actual status of the shared variable on exit as some threads might not have called the error
              !! or even the routine has been called more than once by different threads. 
              !! it is better to raise exceptions outside OMP parallel regions. BTW, by construction this error can never happen
              !! unless the OMP implementation is buggy.
          end if
          !$omp do
          do impl=1,ep%nmpl
              norm_if: if (norm_ok(impl)) then
                  ! Use the method based on the Gaussians
                  ! First determine the highest multipole coefficients which are available. It is required that
                  ! all "lower" multipoles are associated as well.
                  lmax_avail = 0
                  found_non_associated = .false.
                  do l=0,lmax
                      if (associated(ep%mpl(impl)%qlm(l)%q)) then
                          if (found_non_associated) then
                              call f_err_throw('The multipoles for l='//trim(yaml_toa(l))//&
                                   &' are associated, but there are lower multipoles which are &
                                   &not associated. This is not allowed',err_name='BIGDFT_RUNTIME_ERROR')
                          end if
                          lmax_avail = l
                      else
                          found_non_associated = .true.
                      end if
                  end do
                  i3loop: do i3=is3,ie3
                      if (maxval(gaussians3(:,i3,impl))<1.d-20) cycle i3loop
                      ii3 = i3 - nl3 -1
                      !z = real(ii3,kind=8)*hz + shift(3)
                      r(3) = huge(r(3))
                      do j3=j3s,j3e
                          dr = real(ii3+j3*lzd%glr%d%n3i,kind=8)*hz + shift(3) - ep%mpl(impl)%rxyz(3)
                          if (abs(dr)<abs(r(3))) r(3) = dr
                      end do
                      i2loop: do i2=is2,ie2
                          if (maxval(gaussians2(:,i2,impl))<1.d-20) cycle i2loop
                          ii2 = i2 - nl2 - 1
                          !y = real(ii2,kind=8)*hy + shift(2)
                          r(2) = huge(r(2))
                          do j2=j2s,j2e
                              dr = real(ii2+j2*lzd%glr%d%n2i,kind=8)*hy + shift(2) - ep%mpl(impl)%rxyz(2)
                              if (abs(dr)<abs(r(2))) r(2) = dr
                          end do
                          i1loop: do i1=is1,ie1
                              if (maxval(gaussians1(:,i1,impl))<1.d-20) cycle i1loop
                              ii1 = i1 - nl1 - 1
                              !x = real(ii1,kind=8)*hx + shift(1)
                              r(1) = huge(r(1))
                              do j1=j1s,j1e
                                  dr = real(ii1+j1*lzd%glr%d%n1i,kind=8)*hx + shift(1) - ep%mpl(impl)%rxyz(1)
                                  if (abs(dr)<abs(r(1))) r(1) = dr
                              end do
                              !r(1) = x - ep%mpl(impl)%rxyz(1)
                              !r(2) = y - ep%mpl(impl)%rxyz(2)
                              !r(3) = z - ep%mpl(impl)%rxyz(3)
                              rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                              tt = 0.d0
                              ttl = 0.d0
                              do l=0,lmax_avail
                                  ! Calculate the Gaussian as product of three 1D Gaussians
                                  gg = gaussians1(l,i1,impl)*gaussians2(l,i2,impl)*gaussians3(l,i3,impl)
                                  ! Additional modification to avoid divergence
                                  sig = ep%mpl(impl)%sigma(l)
                                  if (l==1) then
                                      gg = gg/(3.d0*sig**2)
                                  else if (l==2) then
                                      gg = gg/(15.d0*sig**4)
                                  end if
                                  norm_check(l,impl) = norm_check(l,impl) + gg*hhh*rnrm2**l
                                  !if (rnrm2<=rmax(impl)**2) then
                                      mm = 0
                                      do m=-l,l
                                          mm = mm + 1
                                          ! For the monopole term, the atomic core charge (which has been expressed using a Gaussian
                                          ! above) has to be added in order to compensate it. In addition the sign has to be
                                          ! switched since the charge density is a positive quantity.
                                          if (l==0) then
                                              !qq = -(ep%mpl(impl)%qlm(l)%q(mm) - real(nelpsp(impl),kind=8))
                                              qq = -(ep%mpl(impl)%qlm(l)%q(mm) - real(ep%mpl(impl)%nzion,kind=8))
                                              !qq = -ep%mpl(impl)%qlm(l)%q(mm)
                                          else
                                              qq = -ep%mpl(impl)%qlm(l)%q(mm)
                                          end if
                                          ttt = qq*&
                                                real(2*l+1,kind=8)*solid_harmonic(0, rmin, l, m, r(1), r(2), r(3))*&
                                                sqrt(4.d0*pi/real(2*l+1,kind=8))*gg!*sqrt(4.d0*pi_param)
                                          tt = tt + ttt
                                          ttl = ttl + ttt
                                      end do
                                  !end if
                              end do
                              density_loc(i1,i2,i3,ithread) = density_loc(i1,i2,i3,ithread) + tt
                              ! Again calculate the multipole values to verify whether they are represented exactly
                              ll = 0
                              m = 0
                              monopole(impl) = monopole(impl) + tt*hhh*&
                                               solid_harmonic(0,0.d0,ll,m,r(1),r(2),r(3))*&
                                               sqrt(4.d0*pi/real(2*ll+1,kind=8))
                              ll = 1
                              do m=-ll,ll
                                  ii = m + 2
                                  dipole(ii,impl) = dipole(ii,impl) + tt*hhh*&
                                                   solid_harmonic(0,0.d0,ll,m,r(1),r(2),r(3))*&
                                                   sqrt(4.d0*pi/real(2*ll+1,kind=8))
                              end do
                              ll = 2
                              do m=-ll,ll
                                  ii = m + 3
                                  quadrupole(ii,impl) = quadrupole(ii,impl) + tt*hhh*&
                                                   solid_harmonic(0,0.d0,ll,m,r(1),r(2),r(3))*&
                                                   sqrt(4.d0*pi/real(2*ll+1,kind=8))
                              end do
                          end do i1loop
                      end do i2loop
                  end do i3loop
              else norm_if
                  ! Use the method based on the analytic formula
                  do l=0,lmax
                      if (associated(ep%mpl(impl)%qlm(l)%q)) then
                          do i3=is3,ie3
                              ii3 = i3 - nl3 -1
                              !z = real(ii3,kind=8)*hz + shift(3)
                              r(3) = huge(r(3))
                              do j3=j3s,j3e
                                  dr = real(ii3+j3*lzd%glr%d%n3i,kind=8)*hz + shift(3) - ep%mpl(impl)%rxyz(3)
                                  if (abs(dr)<abs(r(3))) r(3) = dr
                              end do
                              do i2=is2,ie2
                                  ii2 = i2 - nl2 -1
                                  r(2) = huge(r(2))
                                  do j2=j2s,j2e
                                      dr = real(ii2+j2*lzd%glr%d%n2i,kind=8)*hy + shift(2) - ep%mpl(impl)%rxyz(2)
                                      if (abs(dr)<abs(r(2))) r(2) = dr
                                  end do
                                  do i1=is1,ie1
                                      ii1 = i1 - nl1 -1
                                      r(1) = huge(r(1))
                                      do j1=j1s,j1e
                                          dr = real(ii1+j1*lzd%glr%d%n1i,kind=8)*hx + shift(1) - ep%mpl(impl)%rxyz(1)
                                          if (abs(dr)<abs(r(1))) r(1) = dr
                                      end do
                                      !r(1) = x - ep%mpl(impl)%rxyz(1)
                                      !r(2) = y - ep%mpl(impl)%rxyz(2)
                                      !r(3) = z - ep%mpl(impl)%rxyz(3)
                                      rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                                      rnrm1 = sqrt(rnrm2)
                                      rnrm3 = rnrm1*rnrm2
                                      rnrm5 = rnrm3*rnrm2
                                      select case(l)
                                      case (0)
                                          tt = calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                      case (1)
                                          tt = calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                      case (2)
                                          tt = calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                      case (3)
                                          call f_err_throw('octupole not yet implemented', err_name='BIGDFT_RUNTIME_ERROR')
                                      case default
                                          call f_err_throw('Wrong value of l', err_name='BIGDFT_RUNTIME_ERROR')
                                      end select
                                      potential_loc(i1,i2,i3,ithread) = potential_loc(i1,i2,i3,ithread) + tt
                                  end do
                              end do
                          end do
                      end if
                  end do
              end if norm_if
          end do
          !$omp end do
          !$omp end parallel

          ! Write the PSP info
          !if (verbosity> 0 .and. iproc==0) call write_psp_source(ep, psp_source)
          !!ntype = 0
          !!do impl=1,ep%nmpl
          !!    ! Check whether the info for this type has already been written
          !!    written = .false.
          !!    do itype=1,ntype
          !!        if (trim(ep%mpl(impl)%sym)==trim(multipole_type_names(itype))) then
          !!            written = .true.
          !!            exit
          !!        end if
          !!    end do
          !!    if (.not. written) then
          !!        ntype = ntype + 1
          !!        if (ntype>nmax_multipole_types) call f_err_throw('More than 5000 different multipole types are not allowed')
          !!        multipole_type_names(ntype) = trim(ep%mpl(impl)%sym)
          !!        if (psp_source(impl)==0) then
          !!            call yaml_map(trim(ep%mpl(impl)%sym),'PSP of QM region')
          !!        else if (psp_source(impl)==1) then 
          !!            call yaml_map(trim(ep%mpl(impl)%sym),'built-in PSP')
          !!        end if
          !!    end if
          !!end do
          !call f_free(psp_source)
    
          if ((ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1) > 0) then
             do ithread=0,nthread-1
                ! Gather the total density
                call axpy((ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), 1.0_gp, density_loc(is1,is2,is3,ithread), 1, density(is1,is2,is3), 1)
                ! Gather the total potential, store it directly in pot
                call axpy((ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), 1.0_gp, potential_loc(is1,is2,is3,ithread), 1, pot(is1,is2,is3), 1)
             end do
          end if

             
          call f_free(density_loc)
          call f_free(potential_loc)
          call f_free(gaussians1)
          call f_free(gaussians2)
          call f_free(gaussians3)
    
!! UNCOMMENT FOR TEST          tt = 0.d0
!! UNCOMMENT FOR TEST          do i3=is3,ie3
!! UNCOMMENT FOR TEST              do i2=is2,ie2
!! UNCOMMENT FOR TEST                  do i1=is1,ie1
!! UNCOMMENT FOR TEST                      !write(400+iproc,'(a,3i7,es18.6)') 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
!! UNCOMMENT FOR TEST                      write(400+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
!! UNCOMMENT FOR TEST                      tt = tt + density(i1,i2,i3)*hhh
!! UNCOMMENT FOR TEST                  end do
!! UNCOMMENT FOR TEST              end do
!! UNCOMMENT FOR TEST          end do
!! UNCOMMENT FOR TEST          call plot_density(iproc,nproc,'data'//'multipoles'//'.cube',&
!! UNCOMMENT FOR TEST                        at,at%astruct%rxyz,denspot%pkernel,1,density)

          !write(*,*) 'DEBUG: tt',tt
          
          if (nproc>1) then
              call mpiallred(norm_check, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(monopole, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(dipole, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(quadrupole, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
    
    
          ! Check that the norm is the same as above
          do impl=1,ep%nmpl
              if (norm_ok(impl)) then
                  do l=0,lmax
                      tt = abs(1.d0-norm_check(l,impl))
                      !if (abs(norm(l,impl)-norm_check(l,impl))>1.d-10) then
                      !if (tt>1.d-2) then
                      !    write(*,*) 'ERROR', abs(norm(l,impl)-norm_check(l,impl)), norm_check(l,impl)
                      !    !call f_err_throw('The deviation from normalization of the radial function is too large: '//&
                      !    !    yaml_toa(tt,fmt='(es7.1)'), err_name='BIGDFT_RUNTIME_ERROR')
                      !end if
                      if (tt>1.d-2 .and. iproc==0 .and. verbosity>1) then
                          !write(*,*) 'ERROR', abs(norm(l,impl)-norm_check(l,impl)), norm_check(l,impl)
                          call yaml_warning('The deviation from normalization of the radial function is large: '//&
                              yaml_toa(tt,fmt='(es7.1)'))
                      end if
                  end do
              end if
          end do
    
          call f_free(rmax)
          call f_free(norm_check)
    
    
          if (iproc==0 .and. ep%nmpl > 0 .and. verbosity>0) then
              call yaml_mapping_open('Potential from multipoles')
              call yaml_map('Number of multipole centers',ep%nmpl)
              call yaml_map('Threshold for the norm of the Gaussians',norm_threshold)
              call yaml_map('Minimal radius for divion of the solid harmonics by r^{2l}',rmin)
              call yaml_sequence_open('Details for each multipole')
              do impl=1,ep%nmpl
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open(trim(yaml_toa(impl)))
                  if (norm_ok(impl)) then
                      call yaml_map('Method','Density based on Gaussians')
                      call yaml_map('Sigma of the Gaussians',ep%mpl(impl)%sigma(:),fmt='(f5.3)')
                      tt0 = 1.d0-norm(0,impl)
                      tt1 = 1.d0-norm(1,impl)
                      tt2 = 1.d0-norm(2,impl)
                      call yaml_map('Deviation from normalization for the radial function',&
                          (/tt0,tt1,tt2/),fmt='(1es8.1)')
                      max_error(:) = 0.d0
                      error_meaningful(:) = 0
                      do l=0,lmax
                          if (associated(ep%mpl(impl)%qlm(l)%q)) then
                              mm = 0
                              do m=-l,l
                                  mm = mm + 1
                                  if (l==0) then
                                      !qq = -(ep%mpl(impl)%qlm(l)%q(mm)-real(nelpsp(impl),kind=8))
                                      qq = -(ep%mpl(impl)%qlm(l)%q(mm)-real(ep%mpl(impl)%nzion,kind=8))
                                      !qq = -ep%mpl(impl)%qlm(l)%q(mm)
                                  else
                                      qq = -ep%mpl(impl)%qlm(l)%q(mm)
                                  end if
                                  if (abs(qq)>1.d-8) then
                                      select case (l)
                                      case (0)
                                          max_error(l) = max(max_error(l),monopole(impl)/qq)
                                      case (1)
                                          max_error(l) = max(max_error(l),dipole(mm,impl)/qq)
                                      case (2)
                                          max_error(l) = max(max_error(l),quadrupole(mm,impl)/qq)
                                      end select
                                  else
                                      error_meaningful(l) = 1
                                  end if
                              end do
                              ! Convert to percentaged deviation
                              if (error_meaningful(l)==0) then
                                  max_error(l) = 100.d0*(max_error(l)-1.d0)
                              end if
                          else
                              error_meaningful(l) = 2
                          end if
                      end do
                      do l=0,lmax
                          select case (error_meaningful(l))
                          case (0)
                              output_arr(l) = yaml_toa(max_error(l),fmt='(1es8.1)')
                          case (1)
                              output_arr(l) = 'not meaningful'
                          case (2)
                              output_arr(l) = 'no mltpole l='//yaml_toa(l)
                          end select
                      end do
                      call yaml_map('Maximal deviation from the original values in percent',output_arr(:))
                  else
                      call yaml_map('Method','Analytic expression')
                  end if
                  call yaml_mapping_close()
              end do
              call yaml_sequence_close()
              call yaml_mapping_close()
          end if
    
          if (ep%nmpl > 0) then


             if (present(rxyz) .and. present(dipole_total)) then
                 if (present(quadrupole_total)) then
                     call calculate_dipole_moment(denspot%dpbox, 1, at, rxyz, density, &
                          calculate_quadrupole=.true., &
                          dipole=dipole_total, quadrupole=quadrupole_total, quiet_=.true.)
                 else
                     call calculate_dipole_moment(denspot%dpbox, 1, at, rxyz, density, &
                          calculate_quadrupole=.true., &
                          dipole=dipole_total, quiet_=.true.)
                 end if
             end if

             if (present(rho_mp)) then
                 call f_memcpy(src=density, dest=rho_mp)
             end if

             ! Add the core contribution
             if ((ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)>0) then
                 call axpy((ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), 1.0_gp, density_cores(is1,is2,is3), 1, density(is1,is2,is3), 1)
             end if


             call H_potential('D',denspot%pkernel,density,denspot%V_ext,ehart_ps,0.0_dp,.false.,&
                  quiet='yes')!,rho_ion=denspot%rho_ion)

             if (present(pot_mp)) then
                 call f_memcpy(src=density, dest=pot_mp)
             end if

             !write(*,*) 'ehart_ps',ehart_ps
             !LG: attention to stack overflow here !
             !pot = pot + density
             call daxpy(size(density),1.0_gp,density,1,pot,1)
    !!$
    !!$         !what if this API for axpy? Maybe complicated to understand
    !!$         pot = f_axpy(1.d0,density)
    !!$         !otherwise this is more explicit, but more verbose
    !!$         call f_axpy(a=1.d0,x=density,y=pot)
    !!$         !or this one, for a coefficient of 1
    !!$         pot = .plus_equal. density
    !!$         !maybe this is the better solution?
    !!$         pot = pot .plus. density
    !!$         !as it might be generalized for multiplications and gemms
    !!$         pot= pot .plus. (1.d0 .times. density)
    !!$         !for two matrices in the gemm API
    !!$         C = alpha .times. (A .times. B) .plus. (beta .times. C)
    !!$         !which might be shortcut as
    !!$         C = alpha .t. (A .t. B) .p. (beta .t. C)
    !!$         ! and for transposition
    !!$         C = alpha .t. (A**'t' .t. B) .p. (beta .t. C)
    
          end if
    
          call f_free(norm)
          call f_free(monopole)
          call f_free(dipole)
          call f_free(quadrupole)
          call f_free(norm_ok)
    
    
          !!$$ UNCOMMENT FOR TEST ii = 0
          !!$$ UNCOMMENT FOR TEST do i3=is3,ie3
          !!$$ UNCOMMENT FOR TEST     ii3 = i3 - 15
          !!$$ UNCOMMENT FOR TEST     do i2=is2,ie2
          !!$$ UNCOMMENT FOR TEST         ii2 = i2 - 15
          !!$$ UNCOMMENT FOR TEST         do i1=is1,ie1
          !!$$ UNCOMMENT FOR TEST             ii1 = i1 - 15
          !!$$ UNCOMMENT FOR TEST             ii = ii + 1
          !!$$ UNCOMMENT FOR TEST             !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
          !!$$ UNCOMMENT FOR TEST             write(300+iproc,'(3(a,i6),a,es18.8)') 'i1= ',i1,' i2= ',i2,' i3= ',i3,' val= ',density(i1,i2,i3)
          !!$$ UNCOMMENT FOR TEST             !do impl=1,ep%nmpl
          !!$$ UNCOMMENT FOR TEST             !    r(1) = ep%mpl(impl)%rxyz(1) - x
          !!$$ UNCOMMENT FOR TEST             !    r(2) = ep%mpl(impl)%rxyz(2) - y
          !!$$ UNCOMMENT FOR TEST             !    r(3) = ep%mpl(impl)%rxyz(3) - z 
          !!$$ UNCOMMENT FOR TEST             !    rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
          !!$$ UNCOMMENT FOR TEST             !    rnrm1 = sqrt(rnrm2)
          !!$$ UNCOMMENT FOR TEST             !    tt = spherical_harmonic(l, m, x, y, z)*gaussian(sigma, rnrm1)
          !!$$ UNCOMMENT FOR TEST             !    density(i1,i2,i3) =+ tt
          !!$$ UNCOMMENT FOR TEST             !    !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
          !!$$ UNCOMMENT FOR TEST             !end do
          !!$$ UNCOMMENT FOR TEST         end do
          !!$$ UNCOMMENT FOR TEST     end do
          !!$$ UNCOMMENT FOR TEST end do

          if (present(ixyz0)) then
              ixyz0_(1:3) = ixyz0
          else
              ixyz0_(1:3) = -1
          end if

          if (any((/ixyz0_(1)>=0,ixyz0_(2)>=0,ixyz0_(3)>=0/))) then
              if (.not.all((/ixyz0_(1)>=0,ixyz0_(2)>=0,ixyz0_(3)>=0/))) then
                  call f_err_throw('The coordinates of the point through which &
                                   &the potential shall be plotted must all be non-zero')
              end if
              rxyz_noshift = f_malloc((/3,at%astruct%nat/),id='rxyz_noshift')
              do iat=1,at%astruct%nat
                  rxyz_noshift(1:3,iat) = at%astruct%rxyz(1:3,iat) - shift(1:3)
              end do
              ! Plot of the density, in particular along the axes through the point ixyz0_
              if (present(write_directory)) then
                  filename = trim(write_directory)//'mppot.cube'
              else
                  filename = 'mppot.cube'
              end if
              call plot_density(iproc,nproc,trim(filename),at,rxyz_noshift,denspot%pkernel,nspin=1,rho=density, &
                   ixyz0=ixyz0_)
              call f_free(rxyz_noshift)
          end if
    
          call f_free(density)
          call f_free(density_cores)
          !call f_free(nzatom)
          !call f_free(nelpsp)
          call f_free(rloc)
          !call f_free(npspcode)
          !call f_free(psppar)

      end if multipoles_if

      call f_release_routine()


      !!contains


      !!  function gaussian(sigma, r2) result(g)
      !!    use module_base, only: pi => pi_param
      !!    implicit none
      !!    ! Calling arguments
      !!    real(kind=8),intent(in) :: sigma, r2
      !!    real(kind=8) :: tt, g

      !!    ! Only calculate the Gaussian if the result will be larger than 10^-30
      !!    tt = r2/(2.d0*sigma**2)
      !!    if (tt<=69.07755279d0) then
      !!        g = safe_exp(-tt)
      !!        g = g/sqrt(2.d0*pi*sigma**2)**3
      !!    else
      !!        g = 0.d0
      !!    end if
      !!    !g = g/(sigma**3*sqrt(2.d0*pi)**3)

      !!  end function gaussian

    end subroutine potential_from_charge_multipoles




    function calc_monopole(q, rnrm1) result(mpm)
      implicit none
      ! Calling arguments
      real(dp),dimension(1),intent(in) :: q
      real(dp),intent(in) :: rnrm1
      real(dp) :: mpm

      mpm = -q(1)/rnrm1

    end function calc_monopole


    function calc_dipole(q, r, rnrm3) result(dpm)
      implicit none
      ! Calling arguments
      real(dp),dimension(3),intent(in) :: q
      real(dp),intent(in) :: rnrm3
      real(dp),dimension(3),intent(in) :: r
      real(dp) :: dpm
      real(kind=8),parameter :: factor = sqrt(4.d0*pi/3.d0)

      !dpm = q(1)*r(1) + q(2)*r(2) + q(3)*r(3)
      !dpm = factor*(q(3)*r(1) + q(1)*r(2) + q(2)*r(3))
      dpm = (q(3)*r(1) + q(1)*r(2) + q(2)*r(3))
      dpm = -dpm/rnrm3

    end function calc_dipole


    function calc_quadropole(q, r, rnrm5) result(qpm)
      implicit none
      ! Calling arguments
      real(dp),dimension(5),intent(in) :: q
      real(dp),intent(in) :: rnrm5
      real(dp),dimension(3),intent(in) :: r
      real(dp) :: qpm
      ! Local variables
      real(dp),dimension(3,3) :: qq
      real(kind=8),parameter :: sqrt3=sqrt(3.d0)
      !real(kind=8),parameter :: factor=sqrt(4.d0*pi/15.d0)
      real(kind=8),parameter :: factor=1.d0/sqrt3
      real(kind=8),dimension(3) :: Qr

      !!qq(1,1) = q(1)
      !!qq(2,1) = q(2)
      !!qq(3,1) = q(3)
      !!qq(1,2) = qq(2,1)
      !!qq(2,2) = q(4)
      !!qq(3,2) = q(5)
      !!qq(1,3) = qq(3,1)
      !!qq(2,3) = qq(3,2)
      !!qq(3,3) = 1.0_dp-qq(1,1)-qq(2,2)

      qq(1,1) = factor*(-sqrt3*q(3)+q(5))
      qq(2,1) = factor*q(1)
      qq(3,1) = factor*q(4)
      qq(1,2) = qq(2,1)
      qq(2,2) = factor*(-sqrt3*q(3)-q(5))
      qq(3,2) = factor*q(2)
      qq(1,3) = qq(3,1)
      qq(2,3) = qq(3,2)
      qq(3,3) = factor*2.d0*sqrt3*q(3)

      !qpm = qq(1,1)*r(1)*r(1) + &
      !      qq(2,1)*r(2)*r(1) + &
      !      qq(3,1)*r(3)*r(1) + &
      !      qq(1,2)*r(1)*r(2) + &
      !      qq(2,2)*r(2)*r(2) + &
      !      qq(3,2)*r(3)*r(2) + &
      !      qq(1,3)*r(1)*r(3) + &
      !      qq(2,3)*r(2)*r(3) + &
      !      qq(3,3)*r(3)*r(3)

      ! Calculate r^T*Q*r
      Qr(1) = qq(1,1)*r(1) + qq(2,1)*r(2) + qq(3,1)*r(3)
      Qr(2) = qq(1,2)*r(1) + qq(2,2)*r(2) + qq(3,2)*r(3)
      Qr(3) = qq(1,3)*r(1) + qq(2,3)*r(2) + qq(3,3)*r(3)
      qpm = r(1)*Qr(1) + r(2)*Qr(2) + r(3)*Qr(3)

      qpm = -0.5_dp*qpm/rnrm5

    end function calc_quadropole



    !> Calculate either:
    !! - S^1/2 * K * S^1/2, which is the kernel corresponding to a orthonormal set of support functions.
    !! - S^-1/2 * S * S^-1/2, which is the overlap corresponding to a orthonormal set of support functions.
    !! To keep it simple, always call the matrix in the middle matrix
    subroutine matrix_for_orthonormal_basis(iproc, nproc, meth_overlap, smats, smatl, &
               ovrlp, matrix, operation, weight_matrix_compr)
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, SPARSE_TASKGROUP, &
                                   matrices_null, assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   deallocate_matrices
      use matrix_operations, only: overlapPowerGeneral
      use sparsematrix, only: matrix_matrix_mult_wrapper, gather_matrix_from_taskgroups
      use yaml_output
      implicit none

      ! Calling arguments
      integer :: iproc, nproc,  meth_overlap
      type(sparse_matrix),intent(in) :: smats, smatl
      type(matrices),intent(in) :: matrix
      type(matrices),intent(in) :: ovrlp
      character(len=*),intent(in) :: operation
      real(kind=8),dimension(smatl%nvctrp_tg*smatl%nspin),intent(out) :: weight_matrix_compr

      ! Local variables
      type(matrices),dimension(1) :: inv_ovrlp
      real(kind=8),dimension(:),allocatable :: weight_matrix_compr_tg, proj_ovrlp_half_compr
      real(kind=8) :: max_error, mean_error
      integer :: ioperation
      integer, dimension(1) :: power

      call f_routine(id='matrix_for_orthonormal_basis')

      select case (trim(operation))
      case ('plus')
          ioperation = 2
      case ('minus')
          ioperation = -2
      case default
          call f_err_throw('wrong value of operation')
      end select

      !!if (iproc==0) then
      !!    call yaml_comment('Calculating matrix for orthonormal support functions',hfill='~')
      !!end if

      inv_ovrlp(1) = matrices_null()
      inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='inv_ovrlp(1)%matrix_compr')

      power(1)=ioperation
      call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
           meth_overlap, 1, power, -1, &
           imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
           ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.false., verbosity=0)
      !call f_free_ptr(ovrlp%matrix)

      proj_ovrlp_half_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='proj_mat_compr')
      !if (norbp>0) then
         call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
              matrix%matrix_compr, inv_ovrlp(1)%matrix_compr, proj_ovrlp_half_compr)
      !end if
      !weight_matrix_compr_tg = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='weight_matrix_compr_tg')
      !if (norbp>0) then
         call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
              inv_ovrlp(1)%matrix_compr, proj_ovrlp_half_compr, weight_matrix_compr)
      !end if
      call f_free(proj_ovrlp_half_compr)

      call deallocate_matrices(inv_ovrlp(1))

      !!! Maybe this can be improved... not really necessary to gather the entire matrix
      !!!weight_matrix_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='weight_matrix_compr')
      !!call gather_matrix_from_taskgroups(iproc, nproc, smatl, weight_matrix_compr_tg, weight_matrix_compr)

      !call f_free(weight_matrix_compr_tg)

      !!if (iproc==0) then
      !!    call yaml_comment('Kernel calculated',hfill='~')
      !!end if

      call f_release_routine()

    end subroutine matrix_for_orthonormal_basis




    subroutine write_multipoles_new(ep, ll, units, delta_rxyz, on_which_atom, scaled, monopoles_analytic)
      use yaml_output
      use numerics, only: Bohr_Ang
      use f_precisions, only: db => f_double
      implicit none
      
      ! Calling arguments
      !integer,dimension(nat),intent(in) :: iatype
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: ll
      character(len=*),intent(in) :: units
      real(kind=8),dimension(3,ep%nmpl),intent(in),optional :: delta_rxyz !< can be used to display the difference between the charge center 
                                                                      !! of a support function and its localization center
      integer,dimension(ep%nmpl),intent(in),optional :: on_which_atom !< can be used to display on which atom a given support function multipole is located
      real(kind=8),dimension(ep%nmpl),intent(in),optional :: scaled !< can be used to display by how muched the multipoles have been scaled
      real(kind=8),dimension(ep%nmpl),intent(in),optional :: monopoles_analytic !< can be used to sidplay also the "analytical"
                                                                                !! monopoles (i.e. the ones calculated directly with the overlap matrix, without numerical integration)
      
      ! Local variables
      character(len=20) :: atomname
      character(len=9) :: function_type
      integer :: i, impl, l, m, nit
      real(kind=8) :: factor, convert_units, tt!, get_normalization, get_test_factor
      real(kind=8),dimension(:,:,:),allocatable :: multipoles_tmp
      real(kind=8),dimension(-ll:ll,0:ll) :: multipoles
      logical :: present_delta_rxyz, present_on_which_atom, present_scaled, present_monopoles_analytic

      call f_routine(id='write_multipoles_new')

      present_delta_rxyz = present(delta_rxyz)
      present_on_which_atom = present(on_which_atom)
      present_scaled = present(scaled)
      present_monopoles_analytic = present(monopoles_analytic)


          ! See whether a conversion of the units is necessary
          select case (units)
          case ('angstroem','angstroemd0')
              convert_units = Bohr_Ang
          case ('atomic','atomicd0','bohr','bohrd0','reduced')
              convert_units = 1.0_db
          case default
              convert_units = 1.0_db
              call yaml_warning('units not recognized, no conversion done')
          end select


          call yaml_mapping_open('Multipole coefficients')
          call yaml_map('units',trim(units))
          tt = 0.d0
          do impl=1,ep%nmpl
              tt = tt + ep%mpl(impl)%qlm(0)%q(1)
          end do
          call yaml_map('global monopole',tt,fmt='(es13.6)')
          if (present_monopoles_analytic) then
              call yaml_map('global monopole analytic',sum(monopoles_analytic),fmt='(es13.6)')
          end if
          call yaml_sequence_open('values')
          do impl=1,ep%nmpl
              call yaml_sequence(advance='no')
              call yaml_map('sym',adjustl(trim(ep%mpl(impl)%sym))//' # '//adjustl(trim(yaml_toa(impl,fmt='(i4.4)'))))
              if (present_on_which_atom) then
                  call yaml_map('Atom number',on_which_atom(impl))
              end if
              call yaml_map('r',convert_units*ep%mpl(impl)%rxyz)
              call yaml_map('nzion',ep%mpl(impl)%nzion)
              if (present_delta_rxyz) then
                  call yaml_map('Delta r',convert_units*delta_rxyz(1:3,impl),fmt='(es13.6)')
              end if
              if (any(ep%mpl(impl)%sigma(0:ll)/=0.d0)) then
                  call yaml_map('sigma',ep%mpl(impl)%sigma(0:ll),fmt='(f5.3)')
              endif
              call f_zero(multipoles)
              if (present_monopoles_analytic) then
                  call yaml_map('q0 analytic',monopoles_analytic(impl),fmt='(1es13.6)')
              end if
              do l=0,ll
                  call yaml_map('q'//adjustl(trim(yaml_toa(l))),ep%mpl(impl)%qlm(l)%q(:),fmt='(1es13.6)')
                  multipoles(-l:l,l) = ep%mpl(impl)%qlm(l)%q(1:2*l+1)
                  call yaml_newline()
              end do
              if (present_scaled) then
                  call yaml_map('scaling factor',scaled(impl),fmt='(es9.2)')
              end if
              function_type = guess_type(ll, multipoles)
              call yaml_map('type',trim(function_type))
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()


      call f_release_routine()


          contains
            ! Try to guess the type (s, p, etc.) of a support function
            function guess_type(ll,mp) result(gt)
              implicit none
              ! Calling arguments
              integer,Intent(in) :: ll
              real(kind=8),dimension(-ll:ll,0:ll),intent(in) :: mp
              character(len=9) :: gt
              ! Local variables
              integer :: il, im, ilmax, immax
              real(kind=8) :: maxvalue1, maxvalue2

              
              ! A type is recognized if an element is at least four times as large as all other elements
              maxvalue1 = 0.d0 !the largest element
              maxvalue2 = 0.d0 !the second largest element
              do il=0,ll
                  do im=-il,il
                      if (abs(mp(im,il))>maxvalue2) then
                          maxvalue2 = abs(mp(im,il))
                      end if
                      if (abs(mp(im,il))>maxvalue1) then
                          maxvalue2 = maxvalue1
                          maxvalue1 = abs(mp(im,il))
                          ilmax = il
                          immax = im
                      end if
                  end do
              end do
              if (maxvalue1 > 4.d0*maxvalue2) then
                  ! type recognized
                  select case (ilmax)
                  case (0)
                      gt = 's'
                  case (1)
                      select case(immax)
                      case (-1)
                          gt = 'p_y'
                      case ( 0)
                          gt = 'p_z'
                      case ( 1)
                          gt = 'p_x'
                      end select
                  case (2)
                      select case(immax)
                      case (-2)
                          gt = 'd_xy'
                      case (-1)
                          gt = 'd_yz'
                      case ( 0)
                          gt = 'd_z^2'
                      case ( 1)
                          gt = 'd_xz'
                      case ( 2)
                          gt = 'd_x^2-y^2'
                      end select
                  end select
              else
                  ! type unknown
                  gt = 'unknown'
              end if
            end function guess_type


    end subroutine write_multipoles_new



    function get_test_factor(l, m) result(ff)
      implicit none
      integer,intent(in) :: l, m
      real(kind=8) :: ff
      select case (l)
      case (0)
          ff = 1.d0
      case (1)
          select case (m)
          case (-1)
              ff = 11.d0
          case ( 0)
              ff = 12.d0
          case ( 1)
              ff = 13.d0
          end select
      case (2)
          select case (m)
          case (-2)
              ff = 21.d0
          case (-1)
              ff = 22.d0
          case ( 0)
              ff = 23.d0
          case ( 1)
              ff = 24.d0
          case ( 2)
              ff = 25.d0
          end select
      end select
    end function get_test_factor


    !> Calculates the real spherical harmonic S_lm (multplied by a power or r) for given values of l, m, x, y, z.
    !! The functions are normalized to one within a sphere of radius rmax.
    !! The S_lm are a priori defined without any r, e.g. S_10 = sqrt(3/4pi)z
    !! r_exponent indicates how the function is multiplied by r: The final result is given by S_lm*r^(r_exponent*l), with the
    !! definition of the S_lm ar given above. r_exponent can also be negative, e.g. -1 to yield the "real" spherical harmonics.
    function spherical_harmonic(r_exponent, rmax, l, m, x, y, z) result(sh)
      use module_base, only: pi => pi_param
      implicit none
      ! Calling arguments
      integer,intent(in) :: r_exponent
      integer,intent(in) :: l, m
      real(kind=8),intent(in) :: rmax, x, y, z
      real(kind=8) :: sh

      ! Local variables
      integer,parameter :: l_max=2
      real(kind=8) :: r, r2, rnorm

      if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
      if (l>l_max) call f_err_throw('spherical harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
          err_name='BIGDFT_RUNTIME_ERROR')
      if (abs(m)>l) call f_err_throw('abs of m ('//trim(yaml_toa(m))//') must not be larger than l ('//trim(yaml_toa(l))//')', &
                    err_name='BIGDFT_RUNTIME_ERROR')


      ! Normalization for a sphere of radius rmax
      select case (l)
      case (0)
          ! No need for r, as l=0
          rnorm = sqrt(rmax**3/3.d0)
          !rnorm = 1.d0
          !sh = sqrt(4.d0*pi)*0.5d0*sqrt(1/pi)
          sh = 0.5d0*sqrt(1/pi)/rnorm
      case (1)
          rnorm = sqrt(rmax**5/5.d0)
          !rnorm = sqrt(rmax**2/2.d0)
          !rnorm = 1.d0
          r = sqrt(x**2+y**2+z**2)
          !r = 1.d0
          ! fix for small r (needs proper handling later...)
          if (r<1.d-3) r=1.d-3
          select case (m)
          case (-1)
              !sh = sqrt(4*pi/3.d0)*sqrt(3.d0/(4.d0*pi))*y/r
              !sh = sqrt(3.d0/(4.d0*pi))*y/r
              sh = sqrt(3.d0/(4.d0*pi))*y/rnorm
              !if (.not. with_r) sh = sh/r
          case (0)
              !sh = sqrt(4*pi/3.d0)*sqrt(3.d0/(4.d0*pi))*z/r
              !sh = sqrt(3.d0/(4.d0*pi))*z/r
              sh = sqrt(3.d0/(4.d0*pi))*z/rnorm
              !if (.not. with_r) sh = sh/r
          case (1)
              !sh = sqrt(4*pi/3.d0)*sqrt(3.d0/(4.d0*pi))*x/r
              !sh = sqrt(3.d0/(4.d0*pi))*x/r
              sh = sqrt(3.d0/(4.d0*pi))*x/rnorm
              !if (.not. with_r) sh = sh/r
          end select
          ! Multiply by r^{r_exp*l}
          if (r<1.d0) write(*,*) 'BEFORE: r, m, sh', r, m, sh
          sh = sh*r**r_exponent
          !sh = sh/r
          if (r<1.d0) write(*,*) 'AFTER: r, m, sh, x, y, z', r, m, sh, x, y, z
      case (2)
          rnorm = sqrt(rmax**7/7.d0)
          !rnorm = sqrt(rmax**3/3.d0)
          !rnorm = 1.d0
          r2 = x**2+y**2+z**2
          !r2=1.d0
          ! fix for small r2 (needs proper handling later...)
          if (r2==0.d0) r2=1.d-12
          select case (m)
          case (-2)
              !sh = sqrt(4.d0*pi/5.d0)*0.5d0*sqrt(15.d0/pi)*x*y/r2
              !sh = 0.5d0*sqrt(15.d0/pi)*x*y/r2
              sh = 0.5d0*sqrt(15.d0/pi)*x*y/rnorm
              !if (.not. with_r) sh = sh/r2
          case (-1)
              !sh = sqrt(4.d0*pi/5.d0)*0.5d0*sqrt(15.d0/pi)*y*z/r2
              !sh = 0.5d0*sqrt(15.d0/pi)*y*z/r2
              sh = 0.5d0*sqrt(15.d0/pi)*y*z/rnorm
              !if (.not. with_r) sh = sh/r2
          case (0)
              !sh = sqrt(4.d0*pi/5.d0)*0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
              !sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
              sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/rnorm
              !if (.not. with_r) sh = sh/r2
          case (1)
              !sh = sqrt(4.d0*pi/5.d0)*0.5d0*sqrt(15.d0/pi)*z*x/r2
              !sh = 0.5d0*sqrt(15.d0/pi)*z*x/r2
              sh = 0.5d0*sqrt(15.d0/pi)*z*x/rnorm
              !if (.not. with_r) sh = sh/r2
          case (2)
              !sh = sqrt(4.d0*pi/5.d0)*0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
              !sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
              sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/rnorm
              !if (.not. with_r) sh = sh/r2
          end select
          ! Multiply by r^{r_exp*l}
          sh = sh*r2**r_exponent
          !sh = sh/r2
      end select

      if (abs(sh)>10.d0) then
          write(*,*) 'LARGE VALUE', sh
      end if

    end function spherical_harmonic




    !> Creates the charge density of a Gaussian function, to be used for the local part
    !! of the pseudopotentials (gives the error function term when later processed by the Poisson solver).
    subroutine gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, rx, ry, rz, &
               rloc, zion, multipole_preserving, use_iterator, mp_isf, &
               dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, nrho, pot_ion, rholeaked)
      use module_base
      use module_dpbox, only: denspot_distribution, dpbox_iterator, DPB_POT_ION, dpbox_iter, dpbox_iter_next
      use gaussians, only: mp_exp
      implicit none
      ! Calling arguments
      logical,intent(in) :: perx, pery, perz
      integer,intent(in) :: n1i, n2i, n3i, nrho, i3s, n3pi
      real(kind=8),intent(in) :: rloc, rx, ry, rz, hxh, hyh, hzh
      integer,intent(in) :: nbl1, nbl2, nbl3
      integer,intent(in) :: zion !< ionic charge (integer!)
      logical,intent(in) :: multipole_preserving, use_iterator
      integer,intent(in) :: mp_isf !< interpolating scaling function order for the multipole preserving
      integer,intent(in) :: nmpx, nmpy, nmpz !< sizes of the temporary arrays; if too small the code stops
      real(kind=8),dimension(0:nmpx),intent(inout) :: mpx !< temporary array for the exponetials in x direction
      real(kind=8),dimension(0:nmpy),intent(inout) :: mpy !< temporary array for the exponetials in y direction
      real(kind=8),dimension(0:nmpz),intent(inout) :: mpz !< temporary array for the exponetials in z direction
      type(denspot_distribution),intent(in) :: dpbox
      real(dp),dimension(nrho),intent(inout) :: pot_ion
      real(kind=8),intent(inout) :: rholeaked
    
      ! Local variables
      real(kind=8) :: rlocinv2sq, charge, cutoff, xp, yp, zp
      integer,dimension(2,3) :: nbox
      integer :: i1, i2, i3, isx, iex, isy, iey, isz, iez, j1, j2, j3, ind
      type(dpbox_iterator) :: boxit
      real(gp),parameter :: mp_tiny = 1.e-30_gp
      logical :: gox, goy, goz
    
      call f_routine(id='gaussian_density')

    
      !rloc=at%psppar(0,0,atit%ityp)
      rlocinv2sq=0.5_gp/rloc**2
      charge=real(zion,gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
    
      !write(*,*) 'rloc, charge', rloc, charge
    
      !cutoff of the range
      cutoff=10.0_gp*rloc
      if (multipole_preserving) then
         !We want to have a good accuracy of the last point rloc*10
         !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
         cutoff=cutoff+max(hxh,hyh,hzh)*real(mp_isf,kind=gp)
      end if
      
      if (use_iterator) then
         nbox(1,1)=floor((rx-cutoff)/hxh)
         nbox(1,2)=floor((ry-cutoff)/hyh)
         nbox(1,3)=floor((rz-cutoff)/hzh)
         nbox(2,1)=ceiling((rx+cutoff)/hxh)
         nbox(2,2)=ceiling((ry+cutoff)/hyh)
         nbox(2,3)=ceiling((rz+cutoff)/hzh)
    
         ! Check whether the temporary arrays are large enough
         if (nbox(2,1)-nbox(1,1)>nmpx) then
             call f_err_throw('Temporary array in x direction too small',err_name='BIGDFT_RUNTIME_ERROR')
         end if
         if (nbox(2,2)-nbox(1,2)>nmpy) then
             call f_err_throw('Temporary array in y direction too small',err_name='BIGDFT_RUNTIME_ERROR')
         end if
         if (nbox(2,3)-nbox(1,3)>nmpz) then
             call f_err_throw('Temporary array in z direction too small',err_name='BIGDFT_RUNTIME_ERROR')
         end if
    
         !Separable function: do 1-D integrals before and store it.
         !mpx = f_malloc( (/ nbox(1,1).to.nbox(2,1) /),id='mpx')
         !mpy = f_malloc( (/ nbox(1,2).to.nbox(2,2) /),id='mpy')
         !mpz = f_malloc( (/ nbox(1,3).to.nbox(2,3) /),id='mpz')
         do i1=nbox(1,1),nbox(2,1)
            mpx(i1-nbox(1,1)) = mp_exp(hxh,rx,rlocinv2sq,i1,0,multipole_preserving)
         end do
         do i2=nbox(1,2),nbox(2,2)
            mpy(i2-nbox(1,2)) = mp_exp(hyh,ry,rlocinv2sq,i2,0,multipole_preserving)
         end do
         do i3=nbox(1,3),nbox(2,3)
            mpz(i3-nbox(1,3)) = mp_exp(hzh,rz,rlocinv2sq,i3,0,multipole_preserving)
         end do
         boxit = dpbox_iter(dpbox,DPB_POT_ION,nbox)
    
    
         do while(dpbox_iter_next(boxit))
            xp = mpx(boxit%ibox(1)-nbox(1,1)) * mpy(boxit%ibox(2)-nbox(1,2)) * mpz(boxit%ibox(3)-nbox(1,3))
            pot_ion(boxit%ind) = pot_ion(boxit%ind) - xp*charge
         end do
    
      else
         isx=floor((rx-cutoff)/hxh)
         isy=floor((ry-cutoff)/hyh)
         isz=floor((rz-cutoff)/hzh)
    
         iex=ceiling((rx+cutoff)/hxh)
         iey=ceiling((ry+cutoff)/hyh)
         iez=ceiling((rz+cutoff)/hzh)
    
         ! Check whether the temporary arrays are large enough
         if (iex-isx>nmpx) then
             call f_err_throw('Temporary array in x direction too small',err_name='BIGDFT_RUNTIME_ERROR')
         end if
         if (iey-isy>nmpy) then
             call f_err_throw('Temporary array in y direction too small',err_name='BIGDFT_RUNTIME_ERROR')
         end if
         if (iez-isz>nmpz) then
             call f_err_throw('Temporary array in z direction too small',err_name='BIGDFT_RUNTIME_ERROR')
         end if
    
         !Separable function: do 1-D integrals before and store it.
         !call mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,at%multipole_preserving,mpx,mpy,mpz)
         !mpx = f_malloc( (/ isx.to.iex /),id='mpx')
         !mpy = f_malloc( (/ isy.to.iey /),id='mpy')
         !mpz = f_malloc( (/ isz.to.iez /),id='mpz')
         do i1=isx,iex
            mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,multipole_preserving)
         end do
         do i2=isy,iey
            mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,multipole_preserving)
         end do
         do i3=isz,iez
            mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,multipole_preserving)
         end do
    
         do i3=isz,iez
            zp = mpz(i3-isz)
            if (abs(zp) < mp_tiny) cycle
            !call ind_positions(perz,i3,grid%n3,j3,goz) 
            call ind_positions_new(perz,i3,n3i,j3,goz) 
            j3=j3+nbl3+1
            do i2=isy,iey
               yp = zp*mpy(i2-isy)
               if (abs(yp) < mp_tiny) cycle
               !call ind_positions(pery,i2,grid%n2,j2,goy)
               call ind_positions_new(pery,i2,n2i,j2,goy)
               do i1=isx,iex
                  xp = yp*mpx(i1-isx)
                  if (abs(xp) < mp_tiny) cycle
                  !call ind_positions(perx,i1,grid%n1,j1,gox)
                  call ind_positions_new(perx,i1,n1i,j1,gox)
                  if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                     ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
                     pot_ion(ind)=pot_ion(ind)-xp*charge
                  else if (.not. goz ) then
                     rholeaked=rholeaked+xp*charge
                  endif
               enddo
            enddo
         enddo
    
    
      end if

      call f_release_routine()

    end subroutine gaussian_density


    !> Calculates the solid harmonic S_lm (possibly multplied by a power or r) for given values of l, m, x, y, z.
    !! They are normalized such that the integral over the angle gives r^2, i.e.
    !! \int d\Omega S_{lm}*S_{l'm'}/r^{2l} = r^2 \delta_{ll'}\delta_{mm'}
    !! r_exponent indicates how the function is multiplied by r: The final result is given by S_lm*r^(r_exponent*l), with the
    !! definition of the S_lm given above. r_exponent can also be negative, e.g. -1 to yield the "real" spherical harmonics.
    !! rmin gives the minimal radius that is used for the multiplication by r^(r_exponent*l) (can be used to avoid the divergence
    !! around r=0)
    function solid_harmonic(r_exponent, rmin, l, m, x, y, z) result(sh)
      use module_base, only: pi => pi_param
      implicit none
      ! Calling arguments
      integer,intent(in) :: r_exponent
      integer,intent(in) :: l, m
      real(kind=8),intent(in) :: rmin, x, y, z
      real(kind=8) :: sh

      ! Local variables
      integer,parameter :: l_max=2
      real(kind=8) :: r, r2, r2min

      if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
      if (l>l_max) call f_err_throw('solid harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
          err_name='BIGDFT_RUNTIME_ERROR')
      if (abs(m)>l) call f_err_throw('abs of m ('//trim(yaml_toa(m))//') must not be larger than l ('//trim(yaml_toa(l))//')', &
                    err_name='BIGDFT_RUNTIME_ERROR')


      select case (l)
      case (0)
          ! No need for r, as l=0
          sh = sqrt(1.d0/(4.d0*pi))
      case (1)
          r2 = x**2+y**2+z**2
          !r2min = rmin**2
          !r2 = max(r2,r2min)
          r = sqrt(r2)
          select case (m)
          case (-1)
              sh = sqrt(3.d0/(4.d0*pi))*y
          case (0)
              sh = sqrt(3.d0/(4.d0*pi))*z
          case (1)
              sh = sqrt(3.d0/(4.d0*pi))*x
          end select
          ! Multiply by r^{r_exp*l}
          sh = sh*r**r_exponent
      case (2)
          r2 = x**2+y**2+z**2
          !r2min = rmin**2
          !r2 = max(r2,r2min)
          select case (m)
          case (-2)
              sh = sqrt(15.d0/(4.d0*pi))*x*y
          case (-1)
              sh = sqrt(15.d0/(4.d0*pi))*y*z
          case (0)
              sh = sqrt(5.d0/(16.d0*pi))*(-x**2-y**2+2.d0*z**2)
          case (1)
              sh = sqrt(15.d0/(4.d0*pi))*z*x
          case (2)
              sh = sqrt(15.d0/(16.d0*pi))*(x**2-y**2)
          end select
          ! Multiply by r^{r_exp*l}
          sh = sh*r2**r_exponent
      end select

      !r2 = x**2+y**2+z**2
      !r = sqrt(r2)
      !if (r<1.d-1) then
      !    sh = 0.d0
      !end if

      !if (sh>10.d0) then
      !    r2 = x**2+y**2+z**2
      !    write(*,*) 'LARGE, value, r', sh, sqrt(r2)
      !end if

    end function solid_harmonic

    !!!!!!> Calculates the prefactor of the solid harmonics S_lm  for given values of l, m, x, y, z.
    !!!!!function prefactor_solid_harmonic(l, m, x, y, z) result(psh)
    !!!!!  use module_base, only: pi => pi_param
    !!!!!  implicit none
    !!!!!  ! Calling arguments
    !!!!!  integer,intent(in) :: l, m
    !!!!!  real(kind=8),intent(in) :: x, y, z
    !!!!!  real(kind=8) :: psh

    !!!!!  ! Local variables
    !!!!!  integer,parameter :: l_max=2

    !!!!!  if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
    !!!!!  if (l>l_max) call f_err_throw('solid harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
    !!!!!      err_name='BIGDFT_RUNTIME_ERROR')
    !!!!!  if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


    !!!!!  select case (l)
    !!!!!  case (0)
    !!!!!      psh = 1.d0
    !!!!!  case (1)
    !!!!!      psh = 1.d0
    !!!!!      select case (m)
    !!!!!      case (-1)
    !!!!!          sh = y
    !!!!!      case (0)
    !!!!!          sh = z
    !!!!!      case (1)
    !!!!!          sh = x
    !!!!!      end select
    !!!!!      ! Multiply by r^{r_exp*l}
    !!!!!      sh = sh*r**r_exponent
    !!!!!  case (2)
    !!!!!      r2 = x**2+y**2+z**2
    !!!!!      r2min = rmin**2
    !!!!!      r2 = max(r2,r2min)
    !!!!!      select case (m)
    !!!!!      case (-2)
    !!!!!          sh = sqrt(3.d0)*x*y
    !!!!!      case (-1)
    !!!!!          sh = sqrt(3.d0)*y*z
    !!!!!      case (0)
    !!!!!          sh = sqrt(0.25d0)*(-x**2-y**2+2*z**2)
    !!!!!      case (1)
    !!!!!          sh = sqrt(3.d0)*z*x
    !!!!!      case (2)
    !!!!!          sh = sqrt(0.75d0)*(x**2-y**2)
    !!!!!      end select
    !!!!!      ! Multiply by r^{r_exp*l}
    !!!!!      sh = sh*r2**r_exponent
    !!!!!  end select

    !!!!!  !r2 = x**2+y**2+z**2
    !!!!!  !r = sqrt(r2)
    !!!!!  !if (r<1.d-1) then
    !!!!!  !    sh = 0.d0
    !!!!!  !end if

    !!!!!  !if (sh>10.d0) then
    !!!!!  !    r2 = x**2+y**2+z**2
    !!!!!  !    write(*,*) 'LARGE, value, r', sh, sqrt(r2)
    !!!!!  !end if

    !!!!!end function solid_harmonic



    !> Determines the position of the gridpoint which is closest to a given point.
    function nearest_gridpoint(r, hh) result(ngp)
      implicit none
      ! Calling arguments
      real(kind=8),dimension(3),intent(in) :: r, hh
      real(kind=8),dimension(3) :: ngp
      ! Local variables
      integer :: i, ii
      real(kind=8) :: tt

      do i=1,3
          tt = r(i)/hh(i)
          ii = nint(tt)
          ngp(i) = real(ii,kind=8)*hh(i)
      end do

    end function nearest_gridpoint

    !>apply the Slm operator onto a set of support functions
    subroutine apply_Slm(l,m,geocode,hgrids,acell,psi_ob,nphi,Slmphi,integrate_in_sphere,centers)
      use module_base
      use locreg_operations
      use orbitalbasis
      use bounds, only: geocode_buffers
      implicit none
      integer, intent(in) :: l, m, nphi
      character(len=1), intent(in) :: geocode
      real(gp),dimension(3) :: hgrids,acell
      type(orbital_basis), intent(in) :: psi_ob
      real(wp),dimension(nphi),intent(out) :: Slmphi
      logical, intent(in), optional :: integrate_in_sphere
      real(gp), dimension(3,*), intent(in), optional :: centers
      !local variables
      logical :: perx,pery,perz,sphere
      integer :: npsir,ii1,ii2,ii3,nl1,nl2,nl3,i1,i2,i3,ind
      type(ket) :: psi_it
      type(workarr_sumrho) :: w
      real(wp) :: norm, rmax, tt, x, y, z
      real(wp), dimension(3) :: lrcntr
      real(wp),dimension(:),allocatable :: phi2r, sphi2r
      real(wp), dimension(:), pointer :: sphi_ptr
      
      call f_routine(id='apply_Slm')

      sphere=.false.
      if (present(integrate_in_sphere)) sphere=integrate_in_sphere
      ! Conditions for periodicity
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      !first search the maximum sizes of psir array
      npsir=1
      psi_it=orbital_basis_iterator(psi_ob)
      do while(ket_next_locreg(psi_it))
         npsir=max(npsir,psi_it%lr%d%n1i*psi_it%lr%d%n2i*psi_it%lr%d%n3i)
      end do

      phi2r = f_malloc(npsir,id='phi2r')
      sphi2r = f_malloc(npsir,id='sphi2r')

      call f_zero(Slmphi)
      !iterate over the orbital_basis
      psi_it=orbital_basis_iterator(psi_ob)
      do while(ket_next_locreg(psi_it))
         call initialize_work_arrays_sumrho(psi_it%lr,.true.,w)
         rmax = min(psi_it%lr%d%n1*0.5d0*hgrids(1),psi_it%lr%d%n2*0.5d0*hgrids(2),&
              psi_it%lr%d%n3*0.5d0*hgrids(3))+1.e-3_gp*maxval(hgrids)
         call geocode_buffers(psi_it%lr%geocode,geocode, nl1, nl2, nl3)
         if (present(centers)) then
            lrcntr=centers(:,psi_it%ilr)
         else
            lrcntr=psi_it%lr%locregcenter
         end if
         do while(ket_next(psi_it,ilr=psi_it%ilr))
            if (sphere) call f_zero(sphi2r)
            call daub_to_isf(psi_it%lr,w,psi_it%phi_wvl,phi2r)
            !$omp parallel default(none) &
            !$omp shared(psi_it, hgrids, lrcntr, acell, nl3, nl2, nl1) &
            !$omp shared(perz, pery, perx, sphere, rmax, sphi2r, phi2r, l, m) &
            !$omp private(i3, ii3, z, i2, ii2, y, i1, ii1, x, ind, tt)
            !$omp do
            do i3=1,psi_it%lr%d%n3i
               ii3 = psi_it%lr%nsi3 + i3 - nl3 - 1
               z=ii3*0.5d0*hgrids(3)-lrcntr(3)
               z=closest_image(z,acell(3),perz)
               do i2=1,psi_it%lr%d%n2i
                  ii2 = psi_it%lr%nsi2 + i2 - nl2 - 1
                  y=ii2*0.5d0*hgrids(2)-lrcntr(2)
                  y=closest_image(y,acell(2),pery)
                  do i1=1,psi_it%lr%d%n1i
                     ii1 = psi_it%lr%nsi1 + i1 - nl1 - 1
                     x=ii1*0.5d0*hgrids(1)-lrcntr(1)
                     x=closest_image(x,acell(1),perx)
                     ind = (i3-1)*psi_it%lr%d%n2i*psi_it%lr%d%n1i + (i2-1)*psi_it%lr%d%n1i + i1
                     if (sphere) then
                        if (x**2+y**2+z**2>rmax**2) cycle
                     end if
                     tt = solid_harmonic(0, 0.d0, l, m, x, y, z)
                     tt = tt*sqrt(4.d0*pi/real(2*l+1,gp))
                     sphi2r(ind) = tt*phi2r(ind)
                  end do
               end do
            end do
            !$omp end do
            !$omp end parallel
            sphi_ptr => ob_ket_map(Slmphi,psi_it)
            call isf_to_daub(psi_it%lr, w, sphi2r, sphi_ptr)
         end do
         !deallocations of work arrays
         call deallocate_work_arrays_sumrho(w)
      end do
      call f_free(phi2r)
      call f_free(sphi2r)

      call f_release_routine()

    end subroutine apply_Slm


    !>calculate the multipoles of phi
    subroutine Qlm_phi(lmax,geocode,hgrids,acell,psi_ob,Qlm,integrate_in_sphere,centers)
      use module_base
      use locreg_operations
      use orbitalbasis
      use bounds, only: geocode_buffers
      implicit none
      integer, intent(in) :: lmax
      character(len=1), intent(in) :: geocode
      real(gp),dimension(3) :: hgrids,acell
      type(orbital_basis), intent(in) :: psi_ob
      real(wp), dimension(-lmax:lmax,0:lmax,psi_ob%orbs%norbp), intent(out) :: Qlm
      logical, intent(in), optional :: integrate_in_sphere
      real(gp), dimension(3,*), intent(in), optional :: centers
      !local variables
      logical :: perx,pery,perz,sphere
      integer :: npsir,ii1,ii2,ii3,nl1,nl2,nl3,i1,i2,i3,ind,l,m
      type(ket) :: psi_it
      type(workarr_sumrho) :: w
      real(wp) :: norm, rmax, tt, x, y, z
      real(wp), dimension(3) :: lrcntr
      real(wp),dimension(:),allocatable :: phi2r
      real(wp), dimension(:), pointer :: sphi_ptr
      real(wp),dimension(-lmax:lmax,0:lmax) :: Qlm_work

      call f_routine(id='Qlm_phi')

      sphere=.false.
      if (present(integrate_in_sphere)) sphere=integrate_in_sphere
      ! Conditions for periodicity
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      !first search the maximum sizes of psir array
      npsir=1
      psi_it=orbital_basis_iterator(psi_ob)
      do while(ket_next_locreg(psi_it))
         npsir=max(npsir,psi_it%lr%d%n1i*psi_it%lr%d%n2i*psi_it%lr%d%n3i)
      end do

      call f_zero(Qlm)
      phi2r = f_malloc(npsir,id='phi2r')
      !iterate over the orbital_basis
      psi_it=orbital_basis_iterator(psi_ob)
      do while(ket_next_locreg(psi_it))
         call initialize_work_arrays_sumrho(psi_it%lr,.true.,w)
         rmax = min(psi_it%lr%d%n1*0.5d0*hgrids(1),psi_it%lr%d%n2*0.5d0*hgrids(2),&
              psi_it%lr%d%n3*0.5d0*hgrids(3))+1.e-3_gp*maxval(hgrids)
         call geocode_buffers(psi_it%lr%geocode,geocode, nl1, nl2, nl3)
         if (present(centers)) then
            lrcntr=centers(:,psi_it%ilr)
         else
            lrcntr=psi_it%lr%locregcenter
         end if
         do while(ket_next(psi_it,ilr=psi_it%ilr))
            call daub_to_isf(psi_it%lr,w,psi_it%phi_wvl,phi2r)
            call f_zero(Qlm_work)
            !$omp parallel default(none) &
            !$omp shared(psi_it, hgrids, lrcntr, acell, nl3, nl2, nl1) &
            !$omp shared(perz, pery, perx, sphere, rmax, Qlm_work, phi2r, lmax) &
            !$omp private(i3, ii3, z, i2, ii2, y, i1, ii1, x, ind, tt, l, m)
            !$omp do reduction(+: Qlm_work)
            do i3=1,psi_it%lr%d%n3i
               ii3 = psi_it%lr%nsi3 + i3 - nl3 - 1
               z=ii3*0.5d0*hgrids(3)-lrcntr(3)
               z=closest_image(z,acell(3),perz)
               do i2=1,psi_it%lr%d%n2i
                  ii2 = psi_it%lr%nsi2 + i2 - nl2 - 1
                  y=ii2*0.5d0*hgrids(2)-lrcntr(2)
                  y=closest_image(y,acell(2),pery)
                  do i1=1,psi_it%lr%d%n1i
                     ii1 = psi_it%lr%nsi1 + i1 - nl1 - 1
                     x=ii1*0.5d0*hgrids(1)-lrcntr(1)
                     x=closest_image(x,acell(1),perx)
                     ind = (i3-1)*psi_it%lr%d%n2i*psi_it%lr%d%n1i + (i2-1)*psi_it%lr%d%n1i + i1
                     if (sphere) then
                        if (x**2+y**2+z**2>rmax**2) cycle
                     end if
                     do l=0,lmax
                        do m=-l,l
                           tt = solid_harmonic(0, 0.d0, l, m, x, y, z)
                           tt = tt*sqrt(4.d0*pi/real(2*l+1,gp))
                           Qlm_work(m,l)=Qlm_work(m,l)+tt*phi2r(ind)
                        end do
                     end do
                  end do
               end do
            end do
            !$end do
            !$omp end parallel
            Qlm(-lmax:lmax,0:lmax,psi_it%iorbp) = Qlm_work(-lmax:lmax,0:lmax)
         end do
         !deallocations of work arrays
         call deallocate_work_arrays_sumrho(w)
      end do
      call f_free(phi2r)

      call f_release_routine()

    end subroutine Qlm_phi


    subroutine calculate_multipole_matrix(iproc, nproc, l, m, nphi, phi1, phi2, nphir, hgrids, &
               orbs, collcom, lzd, smmd, smat, locregcenter, ingegration_volume, multipole_matrix)
      use module_base
      use module_types, only: orbitals_data, comms_linear, local_zone_descriptors
      use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
      use sparsematrix_base, only: sparse_matrix, matrices, sparse_matrix_metadata
      use communications_base, only: TRANSPOSE_FULL
      use transposed_operations, only: calculate_overlap_transposed
      use communications, only: transpose_localized
      use bounds, only: geocode_buffers
      use orthonormalization, only: overlap_matrix
      use orbitalbasis
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, l, m, nphi, nphir
      real(kind=8),dimension(nphi),intent(in) :: phi1, phi2
      real(kind=8),dimension(3) :: hgrids
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(3,lzd%nlr),intent(in) :: locregcenter
      type(matrices),intent(inout) :: multipole_matrix
      character(len=*),intent(in) :: ingegration_volume

      ! Local variables
      integer :: ist, istr, iorb, i1, i2, i3, ii1, ii2, ii3, iiorb, ind, ilr, i, nl1, nl2, nl3,npsir
      integer :: i1mod, i2mod, i3mod, is1, ie1, is2, ie2, is3, ie3, ii, nd, nu, j1, j2, j3
      real(kind=8),dimension(:),allocatable :: phi2r, sphi2r, sphi2, phi1t_c, phi1t_f, sphi2t_c, sphi2t_f
      real(kind=8) :: norm, rmax, factor_normalization, tt, x, y, z
      type(workarr_sumrho) :: w
      real(kind=8) :: ddot, dr
      character(len=*),parameter :: sphere = 'sphere', box = 'box'
      logical :: integrate_in_sphere, perx, pery, perz
      integer :: j1s, j1e, j2s, j2e, j3s, j3e
      type(orbital_basis) :: psi_ob
      type(ket) :: psi_it
      real(gp), dimension(3) :: acell
      real(wp), dimension(:), pointer :: sphi_ptr

      call f_routine(id='calculate_multipole_matrix')

      ! Check the arguments
      if (trim(ingegration_volume)==sphere) then
          integrate_in_sphere = .true.
      else if (trim(ingegration_volume)==box) then
          integrate_in_sphere = .false.
      else
          call f_err_throw('wrong argument for ingegration_volume ('//trim(ingegration_volume)//')',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end if

      acell(1)=0.5_gp*hgrids(1)*Lzd%glr%d%n1i
      acell(2)=0.5_gp*hgrids(2)*Lzd%glr%d%n2i
      acell(3)=0.5_gp*hgrids(3)*Lzd%glr%d%n3i

      ! Transform back to wavelets
      sphi2 = f_malloc0(nphi,id='sphi2')

      call orbital_basis_associate(psi_ob,orbs=orbs,phis_wvl=phi2,Lzd=Lzd,id='calculate_multipole_matrix')

      call apply_Slm(l,m,smmd%geocode,hgrids,acell,psi_ob,nphi,sphi2,&
           integrate_in_sphere,centers=locregcenter)

      call orbital_basis_release(psi_ob)

      call overlap_matrix(phi1,nphi,lzd,orbs,collcom,smat,multipole_matrix,sphi2)      

      call f_free(sphi2)

      call f_release_routine()

    end subroutine calculate_multipole_matrix

    pure function closest_image(t,L,periodic) result(x)
      implicit none
      logical, intent(in) :: periodic
      real(gp), intent(in) :: t !< point
      real(gp), intent(in) :: L !< size of the simulation domain
      real(gp) :: x
      !local varaibles
      integer :: j,js,je
      real(gp) :: dx

      if (periodic) then
         js = -1
         je = 1
      else
         js = 0
         je = 0
      end if

      x = huge(x)
      do j=js,je
         dx = t + j*L
         if (abs(dx)<abs(x)) x = dx
      end do

    end function closest_image


     
    subroutine multipole_analysis_driver_new(iproc, nproc, lmax, ixc, smmd, smats, smatm, smatl, &
               ovrlp, ham, kernel, rxyz, method, do_ortho, projectormode, &
               calculate_multipole_matrices, do_check, write_multipole_matrices, &
               nphi, lphi, nphir, hgrids, orbs, collcom, collcom_sr, &
               lzd, at, denspot, orthpar, shift, multipole_matrix_in, ice_obj, filename)
      use module_base
      use module_types, only: orbitals_data, comms_linear, local_zone_descriptors, orthon_data, DFT_local_fields, comms_linear
      use sparsematrix_base, only: sparse_matrix, matrices, sparsematrix_malloc0, assignment(=), &
                                   sparsematrix_malloc, matrices_null, sparsematrix_malloc_ptr, deallocate_matrices, &
                                   SPARSE_TASKGROUP, sparse_matrix_metadata
      use sparsematrix_init, only: distribute_on_tasks
      use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix
      use sparsematrix_io, only: write_sparse_matrix
      use communications, only: transpose_localized
      use orthonormalization, only: orthonormalizelocalized,overlap_matrix
      use module_atoms, only: atoms_data
      use yaml_output
      use multipole_base, only: external_potential_descriptors_null, multipole_set_null, multipole_null, &
           deallocate_external_potential_descriptors
      use orbitalbasis
      use matrix_operations, only: overlapPowerGeneral
      !use Poisson_Solver, only: H_potential
      use Poisson_Solver, except_dp => dp, except_gp => gp
      use foe_base, only: foe_data
      use box
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, lmax, ixc
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smats
      type(sparse_matrix),intent(in) :: smatm
      type(sparse_matrix),intent(in) :: smatl
      type(matrices),intent(in) :: ovrlp
      type(matrices),intent(in) :: ham
      type(matrices),intent(in) :: kernel
      real(kind=8),dimension(3,smmd%nat),intent(in) :: rxyz
      character(len=*),intent(in) :: method
      character(len=*),intent(in) :: do_ortho
      character(len=*),intent(in) :: projectormode
      logical,intent(in) :: calculate_multipole_matrices, do_check, write_multipole_matrices
      integer,intent(in),optional :: nphi, nphir
      real(kind=8),dimension(:),intent(in),optional :: lphi
      real(kind=8),dimension(3),intent(in),optional :: hgrids
      type(orbitals_data),intent(in),optional :: orbs
      type(comms_linear),intent(in),optional :: collcom, collcom_sr
      type(local_zone_descriptors),intent(in),optional :: lzd
      type(atoms_data),intent(in),optional :: at
      type(DFT_local_fields),intent(inout),optional :: denspot
      type(orthon_data),intent(in),optional :: orthpar
      real(kind=8),dimension(3),intent(in),optional :: shift
      type(matrices),dimension(-lmax:lmax,0:lmax),intent(in),target,optional :: multipole_matrix_in
      type(foe_data),intent(inout),optional :: ice_obj
      character(len=*),intent(in),optional :: filename

      ! Local variables
      integer :: methTransformOverlap, iat, ind, ispin, ishift, iorb, jorb, iiorb, l, m, itype, natpx, isatx, kat, n, i, kkat
      integer :: ilr, impl, mm, lcheck, nelpsp, psp_source, j, lwork, ii
      integer, dimension(1) :: power
      logical :: can_use_transposed, all_norms_ok
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      real(kind=8),dimension(:),allocatable :: phi_ortho, Qmat, kernel_ortho, Qmat_tmp,Slmphi
      real(kind=8),dimension(:),allocatable :: eval, work, newoverlap, newmultipole_matrix_large, newoverlap_large
      real(kind=8),dimension(:,:),allocatable :: Qmat_tilde, kp, locregcenter, overlap_small, tmpmat, tempmat
      real(kind=8),dimension(:,:,:),pointer :: atomic_multipoles
      real(kind=8),dimension(:),pointer :: atomic_monopoles_analytic
      real(kind=8),dimension(:,:,:),allocatable :: test_pot
      real(kind=8),dimension(:,:,:,:),allocatable :: lmp_extracted
      real(kind=8),dimension(:,:),allocatable :: projx
      real(kind=8),dimension(:,:),allocatable :: kernel_extracted, multipole_extracted
      real(kind=8) :: q, tt, rloc, max_error, mean_error
      type(matrices) :: multipole_matrix
      !type(matrices),target :: multipole_matrix_
      type(matrices) :: newovrlp, ovrlp_large, multipole_matrix_large
      type(matrices),dimension(-1:1,0:1) :: lower_multipole_matrices
      type(matrices),dimension(1) :: inv_ovrlp
      logical :: perx, pery, perz
      logical,dimension(:,:),allocatable :: neighborx
      integer,dimension(:),allocatable :: nx
      character(len=20),dimension(:),allocatable :: names
      real(kind=8) :: rr1, rr2, rr3
      real(kind=8),dimension(3) :: dipole_check
      real(kind=8),dimension(3,3) :: quadrupole_check
      type(external_potential_descriptors) :: ep_check
      type(matrices),dimension(24) :: rpower_matrix
      type(orbital_basis) :: psi_ob
      real(gp), dimension(3) :: acell, center
      character(len=*),parameter :: no='no', yes='yes'
      type(external_potential_descriptors) :: ep
      !character(len=*),parameter :: projectormode='verynew'!'old'
      !character(len=*),parameter :: do_ortho = no!yes
      integer :: is1, ie1, is2, ie2, is3, ie3, ioffset, icheck, nmax
      real(kind=8),dimension(:,:,:),allocatable :: rho_exact, rho_mp, pot_exact, pot_mp
      integer,parameter :: ncheck = 5
      real(kind=8),dimension(ncheck),parameter :: check_threshold = [ 1.d-12 , &
                                                                      1.d-10 , &
                                                                      1.d-8 , &
                                                                      1.d-6 , &
                                                                      1.d-4]
      real(kind=8),dimension(ncheck) :: charge_error, charge_total, potential_error, potential_total
      type(cell) :: mesh
      character(len=2) :: lname, mname
      character(len=14) :: matname


      call f_routine(id='multipole_analysis_driver')

      perx=(smmd%geocode /= 'F')
      pery=(smmd%geocode == 'P')
      perz=(smmd%geocode /= 'F')

      ! Check that the proper optional arguments are present
      if (trim(do_ortho)==yes .and. calculate_multipole_matrices) then
          if (.not.present(nphi) .or. &
              .not.present(orbs) .or. &
              .not.present(lzd) .or. &
              .not.present(collcom) .or. &
              .not.present(lphi) .or. &
              .not.present(orthpar)) then
              call f_err_throw('do_ortho: not all required optional arguments are present')
          end if
      end if
      if (trim(projectormode)=='full') then
          if (.not.present(orbs) .or. &
              .not.present(lzd) .or. &
              .not.present(nphi) .or. &
              .not.present(lphi) .or. &
              .not.present(collcom) .or. &
              .not.present(collcom_sr)) then
              call f_err_throw('projectormode==full: not all required optional arguments are present')
          end if
      end if
      !!if (lmax>0) then
      !!    if (.not.present(orbs) .or. &
      !!        .not.present(lzd)) then
      !!        call f_err_throw('lmax>0: not all required optional arguments are present')
      !!    end if
      !!    if (.not.calculate_multipole_matrices) then
      !!        call f_err_throw('The multipole matrices must be calculated in-situ for lmax>0')
      !!    end if
      !!end if

      if (calculate_multipole_matrices) then
          if (.not.present(orbs) .or. &
              .not.present(lzd) .or. &
              .not.present(nphi) .or. &
              .not.present(nphir) .or. &
              .not.present(lphi) .or. &
              .not.present(hgrids) .or. &
              .not.present(collcom)) then
              call f_err_throw('calculate_multipole_matrices .true.: not all required optional arguments are present')
          end if
      else
          if (.not.present(multipole_matrix_in)) then
              call f_err_throw('multipole_matrix_in .false.: not all required optional arguments are present')
          end if
      end if

      if (do_check) then
          if (.not.present(denspot) .or. &
              .not.present(shift) .or. &
              .not.present(lzd) .or. &
              .not.present(at)) then
              call f_err_throw('calculate_multipole_matrices .true.: not all required optional arguments are present')
          end if
      end if

      if (present(lphi) .and. present(nphi)) then
          if (size(lphi)<nphi) then
              call f_err_throw('wrong size of lphi')
          end if
      end if


      if (iproc==0) then
          call yaml_comment('Atomic multipole analysis, new approach',hfill='=')
          call yaml_map('Method',trim(method))
          call yaml_map('Projector mode',trim(projectormode))
          call yaml_map('Orthogonalized support functions',trim(do_ortho))
      end if

      if (calculate_multipole_matrices) then
          call unitary_test_multipoles(iproc, nproc, nphi, nphir, orbs, lzd, smmd, smats, collcom, hgrids)
      end if

      ! Check the method
      if (trim(method)/='projector' .and. trim(method)/='loewdin') then
          call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (trim(method)=='projector') then
          call f_err_throw('projector method is deprecated',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (trim(do_ortho)/='no' .and. trim(do_ortho)/='yes') then
          call f_err_throw('wrong do_ortho',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      select case (trim(projectormode))
      case ('none','simple','full')
          ! everything ok
      case default
          call f_err_throw('wrong projectormode',err_name='BIGDFT_RUNTIME_ERROR')
      end select

      if (write_multipole_matrices) then
          if (.not.present(filename)) then
              call f_err_throw('filename not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      end if



      !!multipole_matrix_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='multipole_matrix_large')
      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='kernel_ortho')

      ! Calculate the "effective" kernel. For Mulliken, this is K*S, whereas for
      ! Loewdin it is S^-1/2*K*S^-1/2.
      if (do_ortho==yes) then
          methTransformOverlap = 1020
          call matrix_for_orthonormal_basis(iproc, nproc, methTransformOverlap, smats, smatl, &
               ovrlp, kernel, 'plus', kernel_ortho)
      else if (do_ortho==no) then
          ovrlp_large = matrices_null()
          ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, id='ovrlp_large%matrix_compr')
          call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
               smat_in=ovrlp%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
          ! Should use the highlevel wrapper...
          call matrix_matrix_mult_wrapper(iproc, nproc, smatl, kernel%matrix_compr, ovrlp_large%matrix_compr, kernel_ortho)
          call deallocate_matrices(ovrlp_large)
      end if

      ! Parallelization over the atoms
      call distribute_on_tasks(smmd%nat, bigdft_mpi%iproc, bigdft_mpi%nproc, natpx, isatx)

      neighborx = f_malloc((/smats%nfvctr,natpx/),id='neighborx')
      nx = f_malloc(natpx,id='nx')
      call determine_submatrix_sizes(natpx, isatx, smmd, smatl, neighborx, nx, nmax)
      projx = f_malloc((/nmax**2,natpx/),id='projx')

      ! Calculate the matrix for the projector matrix, which is S^-1 for Mulliken and Id for Loewdin.
      ! However, to be consistent (error cancellation of the inverse etc), we calculate for Loewdin the matrix as [S^-1/2*S*S^-1/2]^-1
      inv_ovrlp(1) = matrices_null()
      inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, id='inv_ovrlp%matrix_compr')
      newovrlp = matrices_null()
      newovrlp%matrix_compr = sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='newovrlp%matrix_compr')
      if (do_ortho==yes) then
          methTransformOverlap = 1020
          newoverlap_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='newoverlap_large')
          ovrlp_large = matrices_null()
          ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, id='ovrlp_large%matrix_compr')
          call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
               smat_in=ovrlp%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
          call matrix_for_orthonormal_basis(iproc, nproc, methTransformOverlap, smats, smatl, &
               ovrlp, ovrlp_large, 'minus', newoverlap_large)
          call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'large_to_small', &
               lmat_in=newoverlap_large, smat_out=newovrlp%matrix_compr)
          call deallocate_matrices(ovrlp_large)
          call f_free(newoverlap_large)
      else
          call f_memcpy(src=ovrlp%matrix_compr, dest=newovrlp%matrix_compr)
      end if
      power=1
      if (present(ice_obj)) then
          call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
               1020, 1, power, -1, &
               imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
               ovrlp_mat=newovrlp, inv_ovrlp_mat=inv_ovrlp, &
               check_accur=.false., ice_obj=ice_obj)
      else
          call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
               1020, 1, power, -1, &
               imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
               ovrlp_mat=newovrlp, inv_ovrlp_mat=inv_ovrlp, &
               check_accur=.false.)
      end if
      call deallocate_matrices(newovrlp)

      Qmat = sparsematrix_malloc(smatl,iaction=SPARSE_TASKGROUP,id='Qmat')
      atomic_multipoles = f_malloc0_ptr((/-lmax.to.lmax,0.to.lmax,1.to.smmd%nat/),id='atomic_multipoles')


      multipole_matrix = matrices_null()
      multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='multipole_matrix%matrix_compr')


      ! Choose as reference point the midpoint of the simulation cell, in order to avoid
      ! problems with periodic BC (in this way the reference point is always the same and never a periodic image)
      center(1:3) = 0.5d0*smmd%cell_dim(1:3)
      !!!!SM: This is not really well done...
      if (calculate_multipole_matrices) then
          locregcenter = f_malloc((/3,lzd%nlr/),id='locregcenter')
          do ilr=1,lzd%nlr
              locregcenter(1:3,ilr) = lzd%llr(ilr)%locregcenter(1:3) !+ (/1.d0,2.d0,3.d0/)
          end do
      !!!else
      !!!  locregcenter = f_malloc((/3,smmd%nfvctr/),id='locregcenter')
      !!!  do i=1,smmd%nfvctr
      !!!      iat = smmd%on_which_atom(i)
      !!!      locregcenter(1:3,i) = smmd%rxyz(1:3,iat)
      !!!  end do
      end if

      acell(1)=smmd%cell_dim(1)
      acell(2)=smmd%cell_dim(2)
      acell(3)=smmd%cell_dim(3)

      do l=0,lmax
          write(*,*) 'l, lmax',l, lmax
          do m=-l,l

              call f_zero(multipole_matrix%matrix_compr)

              ! Calculate the multipole matrix
              if (calculate_multipole_matrices) then
                  call calculate_multipole_matrix(iproc, nproc, l, m, nphi, lphi, lphi, nphir, hgrids, &
                       orbs, collcom, lzd, smmd, smats, locregcenter, 'box', multipole_matrix) 
                  if (write_multipole_matrices) then
                      write(lname,'(i0)') l
                      write(mname,'(i0)') m
                      matname = 'mpmat_'//trim(lname)//'_'//trim(mname)//'.bin'
                      call write_sparse_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
                           smats, multipole_matrix, &
                           filename=trim(filename//matname))
                  end if
              else
                  call f_memcpy(src=multipole_matrix_in(m,l)%matrix_compr, dest=multipole_matrix%matrix_compr)
              end if

              if (do_ortho==yes) then
                  ! Calculate S^-1/2*P*S^-1/2, where P is the multipole matrix
                  methTransformOverlap = 1020
                  newoverlap_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='newoverlap_large')
                  ovrlp_large = matrices_null()
                  ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, id='ovrlp_large%matrix_compr')
                  call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
                       smat_in=multipole_matrix%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
                  call matrix_for_orthonormal_basis(iproc, nproc, methTransformOverlap, smats, smatl, &
                       ovrlp, ovrlp_large, 'minus', newoverlap_large)
                  call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'large_to_small', &
                       lmat_in=newoverlap_large, smat_out=multipole_matrix%matrix_compr)
                  call deallocate_matrices(ovrlp_large)
                  call f_free(newoverlap_large)
              end if

              if (l<=1) then
                  lower_multipole_matrices(m,l) = matrices_null()
                  lower_multipole_matrices(m,l)%matrix_compr = &
                      sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='lower_multipole_matrix%matrix_compr')
                  call f_memcpy(src=multipole_matrix%matrix_compr,dest=lower_multipole_matrices(m,l)%matrix_compr)
              end if



              do kat=1,natpx
                  kkat = kat + isatx
                  n = nx(kat)
                  qmat_tilde = f_malloc((/n,n/),id='qmat_tilde')
                  kp = f_malloc((/n,n/),id='kp')
                  kernel_extracted = f_malloc((/n,n/),id='kernel_extracted')
                  multipole_extracted = f_malloc((/n,n/),id='multipole_extracted')
                  call extract_matrix(smats, multipole_matrix%matrix_compr, &
                      neighborx(1,kat), n, multipole_extracted)
                  ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
                  call dscal(n**2, -1.d0, multipole_extracted(1,1), 1)
                  call extract_matrix(smatl, kernel_ortho, neighborx(1,kat), n, kernel_extracted)
                  if (l>0) then
                      call correct_multipole_origin(smmd%nat, l, m, n, smmd%nfvctr, natpx, kat, kkat, &
                           smmd, smats, smmd%rxyz, neighborx, perx, pery, perz, acell, &
                           lower_multipole_matrices, multipole_extracted)
                  end if
                  !do i=1,n
                  !    do j=1,n
                  !        write(*,*) 'i, j, multipole_extracted(j,i)', i, j, multipole_extracted(j,i)
                  !    end do
                  !end do
                  if (trim(method)=='loewdin') then
                          call extract_matrix(smatl, inv_ovrlp(1)%matrix_compr, neighborx(1,kat), n, projx(1,kat))
                          iiorb = 0
                          do iorb=1,smats%nfvctr
                              if (neighborx(iorb,kat)) then
                                  iiorb = iiorb + 1
                                  if (smmd%on_which_atom(iorb)/=kkat) then
                                      do jorb=1,n
                                          projx((iiorb-1)*n+jorb,kat) = 0.d0
                                      end do
                                  end if
                              end if
                          end do
                  end if
                  !do i=1,n**2
                  !    write(*,*) 'i, j, projx(i,kat)', i, j, projx(i,kat)
                  !end do
                  call gemm('n', 'n', n, n, n, 1.d0, kernel_extracted(1,1), n, &
                       projx(1,kat), n, 0.d0, qmat_tilde(1,1), n)
                  call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, multipole_extracted(1,1), n, 0.d0, kp(1,1), n)
                  call f_free(kernel_extracted)
                  call f_free(multipole_extracted)
                  tt = 0.d0
                  do i=1,n
                      tt = tt + kp(i,i)
                  end do
                  atomic_multipoles(m,l,kkat) = tt
                  call f_free(qmat_tilde)
                  call f_free(kp)
              end do

          end do
      end do


      if (calculate_multipole_matrices) then
          call f_free(locregcenter)
      end if

      call mpiallred(atomic_multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)


      ! The monopole term should be the net charge, i.e. add the positive atomic charges
      do iat=1,smmd%nat
          itype = smmd%iatype(iat)
          q = real(smmd%nelpsp(itype),kind=8)
          atomic_multipoles(0,0,iat) = atomic_multipoles(0,0,iat) + q
      end do


      names = f_malloc_str(len(names),smmd%nat,id='names')
      do iat=1,smmd%nat
          itype = smmd%iatype(iat)
          names(iat) = smmd%atomnames(itype)
      end do


      ep = external_potential_descriptors_null()
      ep%nmpl = smmd%nat
      allocate(ep%mpl(ep%nmpl))
      do impl=1,ep%nmpl
          ep%mpl(impl) = multipole_set_null()
          allocate(ep%mpl(impl)%qlm(0:lmax))
          ep%mpl(impl)%rxyz = smmd%rxyz(1:3,impl)
          ep%mpl(impl)%sym = trim(names(impl))
          if (present(at)) then
              call get_psp_info(ep%mpl(impl)%sym, ixc, smmd, nelpsp, psp_source, rloc, at%psppar)
          else
              call get_psp_info(ep%mpl(impl)%sym, ixc, smmd, nelpsp, psp_source, rloc)
          end if
          if (psp_source/=0 .and. iproc==0) then
              call yaml_warning('Taking internal PSP information for multipole '//trim(yaml_toa(impl)))
          end if
          ep%mpl(impl)%nzion = nelpsp
          ep%mpl(impl)%sigma(0:lmax) = rloc
          do l=0,lmax
              ep%mpl(impl)%qlm(l) = multipole_null()
              !if (l>=3) cycle
              ep%mpl(impl)%qlm(l)%q = f_malloc0_ptr(2*l+1,id='q')
              mm = 0
              !if (l>0) cycle
              do m=-l,l
                  mm = mm + 1
                  ep%mpl(impl)%qlm(l)%q(mm) = atomic_multipoles(m,l,impl)
              end do
          end do
      end do

      if (iproc==0) then
          call yaml_comment('Final result of the multipole analysis',hfill='~')
          call write_multipoles_new(ep, lmax, smmd%units)
      end if


      if (do_check) then
          ! Calculate the total dipole moment resulting from the previously calculated multipoles.
          ! This is done by calling the following routine (which actually calculates the potential, but also
          ! has the option to calculate the dipole on the fly).
          test_pot = f_malloc((/size(denspot%V_ext,1),size(denspot%V_ext,2),size(denspot%V_ext,3)/),id='test_pot')
          if (iproc==0) call yaml_sequence_open('Checking the total multipoles based on the atomic multipoles')
          is1 = 1
          ie1 = denspot%dpbox%mesh%ndims(1)
          is2 = 1
          ie2 = denspot%dpbox%mesh%ndims(2)
          is3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1
          ie3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
                denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2)
          rho_exact = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='rho_exact')
          rho_mp = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='rho_mp')
          pot_exact = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='pot_exact')
          pot_mp = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='pot_mp')
          do lcheck=0,lmax
              ep_check = external_potential_descriptors_null()
              ep_check%nmpl = ep%nmpl
              allocate(ep_check%mpl(ep_check%nmpl))
              do impl=1,ep_check%nmpl
                  ep_check%mpl(impl) = multipole_set_null()
                  allocate(ep_check%mpl(impl)%qlm(0:lmax))
                  ep_check%mpl(impl)%rxyz = ep%mpl(impl)%rxyz
                  ep_check%mpl(impl)%sym = ep%mpl(impl)%sym
                  ep_check%mpl(impl)%nzion = ep%mpl(impl)%nzion
                  do l=0,lmax
                      ep_check%mpl(impl)%sigma(l) = ep%mpl(impl)%sigma(l)
                      ep_check%mpl(impl)%qlm(l) = multipole_null()
                      if (l>lcheck) cycle
                      ep_check%mpl(impl)%qlm(l)%q = f_malloc0_ptr(2*l+1,id='q')
                      mm = 0
                      do m=-l,l
                          mm = mm + 1
                          ep_check%mpl(impl)%qlm(l)%q(mm) = ep%mpl(impl)%qlm(l)%q(mm)
                      end do
                  end do
              end do
              call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
                   denspot%V_ext(1,1,1,1), 1, test_pot(1,1,1), 1)
              call potential_from_charge_multipoles(iproc, nproc, at, denspot, ep_check, 1, &
                   denspot%dpbox%mesh%ndims(1), 1, denspot%dpbox%mesh%ndims(2), &
                   denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
                   denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
                   denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
                   denspot%dpbox%mesh%hgrids(1),denspot%dpbox%mesh%hgrids(2),denspot%dpbox%mesh%hgrids(3), &
                   shift, verbosity=0, ixc=ixc, lzd=lzd, pot=test_pot, &
                   rxyz=rxyz, dipole_total=dipole_check, quadrupole_total=quadrupole_check, &
                   all_norms_ok=all_norms_ok, &
                   rho_mp=rho_mp, pot_mp=pot_mp)
              if (.not. all_norms_ok) then
                  call f_err_throw('When checking the previously calculated multipoles, all norms should be ok')
              end if
              dipole_check=dipole_check/Debye_AU  ! au2debye              

              !# NEW: compare the density and potential ##########################
              if (smatl%nspin/=1) then
                  call f_err_throw('Multipole analysis check not yet ready for nspin>1')
              end if
              ! Get the exact charge density
              ioffset = denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*&
                        denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,4)
              !write(*,*) 'MP: ioffset', ioffset
              call f_memcpy(n=(ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), &
                   src=denspot%rhov(ioffset+1), dest=rho_exact(is1,is2,is3))
              call f_memcpy(src=rho_exact, dest=pot_exact)
!!$              call H_potential('D',denspot%pkernel,pot_exact,denspot%V_ext,tt,0.0_dp,.true.,&
!!$                   quiet='yes')
              call Electrostatic_Solver(denspot%pkernel,pot_exact,pot_ion=denspot%V_ext)
              !mesh=cell_new(smmd%geocode,denspot%pkernel%ndims,denspot%pkernel%hgrids)
              call compare_charge_and_potential(denspot%dpbox%bitp,&!iproc, is1, ie1, is2, ie2, is3, ie3, &
                   smmd%nat, &
                   rho_exact, rho_mp, pot_exact, pot_mp, denspot%pkernel, rxyz, &
                   ncheck, check_threshold, charge_error, charge_total, potential_error, potential_total)
              !# NEW: compare the density and potential ##########################
              if (iproc==0) then
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open('Up to multipole l='//trim(yaml_toa(lcheck)))
                  call yaml_mapping_open('Electric Dipole Moment (Debye)')
                  call yaml_map('P vector',dipole_check(1:3),fmt='(1es13.4)')
                  call yaml_map('norm(P)',sqrt(sum(dipole_check**2)),fmt='(1es14.6)')
                  call yaml_mapping_close()
                  call yaml_mapping_open('Quadrupole Moment (AU)')
                  call yaml_map('Q matrix',quadrupole_check,fmt='(1es13.4)')
                  call yaml_map('trace',quadrupole_check(1,1)+quadrupole_check(2,2)+quadrupole_check(3,3),fmt='(es12.2)')
                  call yaml_mapping_close()
                  call yaml_sequence_open('Average relative error of resulting potential in the Exterior region')
                  !call yaml_sequence_open('density threshold for check')
                  do icheck=1,ncheck
                      !call yaml_mapping_open('density threshold for check',check_threshold(icheck))
                      call yaml_sequence(advance='no')
                      call yaml_mapping_open(flow=.true.)
                      call yaml_map('Thr',check_threshold(icheck),fmt='(es9.1)')
                      call yaml_map('Ext. Vol. %',charge_total(icheck)/&
                           (denspot%dpbox%mesh%volume_element*product(real(denspot%dpbox%mesh%ndims,gp))),&
                           fmt='(2pf5.1)')
                      call yaml_map('int(V)',potential_total(icheck),fmt='(es10.3)')
                      call yaml_map('Err %',potential_error(icheck)/potential_total(icheck),fmt='(2pf5.1)')
                      call yaml_map('int(rho)',charge_error(icheck),fmt='(es10.3)')
!!$                      !call yaml_mapping_open('density threshold for check is'//&
!!$                      !     &trim(yaml_toa(check_threshold(icheck),fmt='(es9.2)')))
!!$                      call yaml_mapping_open('rho',flow=.true.)
!!$                      call yaml_map('int(q-q_exact))',charge_error(icheck),fmt='(es10.3)')
!!$                      call yaml_map('int(q_exact)',charge_total(icheck),fmt='(es10.3)')
!!$                      call yaml_map('ratio',charge_error(icheck)/charge_total(icheck),fmt='(es10.3)')
!!$                      call yaml_mapping_close()
!!$                      call yaml_mapping_open('pot',flow=.true.)
!!$                      call yaml_map('int(V-V_exact))',potential_error(icheck),fmt='(es10.3)')
!!$                      call yaml_map('int(V_exact)',potential_total(icheck),fmt='(es10.3)')
!!$                      call yaml_map('ratio',potential_error(icheck)/potential_total(icheck),fmt='(es10.3)')
!!$                      call yaml_mapping_close()
                      call yaml_mapping_close()
                  end do
                  call yaml_sequence_close()
                  !call yaml_mapping_close()
                  call yaml_mapping_close()
               end if
               call deallocate_external_potential_descriptors(ep_check)
            end do
          if (iproc==0) call yaml_sequence_close()
          call f_free(test_pot)
          call f_free(rho_exact)
          call f_free(rho_mp)
          call f_free(pot_exact)
          call f_free(pot_mp)
      end if

      do l=0,min(1,lmax)
          do m=-l,l
              call deallocate_matrices(lower_multipole_matrices(m,l))
          end do
      end do

      if (trim(method)=='loewdin') then
          call deallocate_matrices(inv_ovrlp(1))
      end if

      call f_free_str(len(names),names)
      call deallocate_matrices(multipole_matrix)
      call f_free(kernel_ortho)
      call f_free(Qmat)
      !if (do_ortho==yes .and. calculate_multipole_matrices) then
      !    call f_free(phi_ortho)
      !end if
      call f_free(projx)
      call f_free(nx)
      call f_free(neighborx)
      call f_free_ptr(atomic_multipoles)
      !!call f_free(multipole_matrix_large)
      call deallocate_external_potential_descriptors(ep)

      if (iproc==0) then
          call yaml_comment('Atomic multipole analysis done',hfill='=')
      end if

      call f_release_routine()

  end subroutine multipole_analysis_driver_new



  subroutine extract_matrix(smat, matrix_compr, neighbor, n, matrix)
    use module_base
    use sparsematrix_base,only: sparse_matrix, matrices
    use sparsematrix_init, only: matrixindex_in_compressed
    implicit none

    ! Calling arguments
    type(sparse_matrix),intent(in) :: smat
    real(kind=8),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: matrix_compr
    logical,dimension(smat%nfvctr),intent(in) :: neighbor
    integer,intent(in) :: n
    real(kind=8),dimension(n,n),intent(out) :: matrix

    ! Local variables
    integer :: icheck, ii, jj, i, j, ind
    !logical :: optional_present
    integer,dimension(:),allocatable :: lookup

    call f_routine(id='extract_matrix')

    !optional_present = present(ilup)

    lookup = f_malloc(smat%nfvctr,id='lookup')
    ii = 0
    do i=1,smat%nfvctr
        if (neighbor(i)) then
            ii = ii + 1
            lookup(i) = ii
        end if
    end do


    icheck = 0
    ii = 0
    !SM: The function matrixindex_in_compressed is rather expensive, so probably worth to use OpenMP all the time
    !$omp parallel default(none) &
    !$omp shared(smat, neighbor, lookup, matrix, matrix_compr, icheck) &
    !$omp private(i, jj, ii, j, ind)
    !$omp do schedule(guided) reduction(+: icheck)
    do i=1,smat%nfvctr
        if (neighbor(i)) then
            jj = 0
            ii = lookup(i)
            do j=1,smat%nfvctr
                if (neighbor(j)) then
                    icheck = icheck + 1
                    jj = jj + 1
                    !if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                    ind =  matrixindex_in_compressed(smat, j, i)
                    if (ind>0) then
                        matrix(jj,ii) = matrix_compr(ind-smat%isvctrp_tg)
                    else
                        matrix(jj,ii) = 0.d0
                    end if
                    !if (optional_present) then
                    !    ilup(1,jj,ii) = j
                    !    ilup(2,jj,ii) = i
                    !end if
                end if
            end do
        end if
    end do
    !$omp end do
    !$omp end parallel
    if (icheck>n**2) then
        call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
            &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
    end if

    call f_free(lookup)

    call f_release_routine()

  end subroutine extract_matrix



    subroutine supportfunction_centers(nat, rxyz, nphidim, phi, nphirdim, &
               norb, norbp, isorb, in_which_locreg, lzd, com)
      use module_base
      use module_types, only: local_zone_descriptors
      use bounds, only: geocode_buffers
      use locreg_operations
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: nat, nphidim, nphirdim, norb, norbp, isorb
      integer,dimension(norb),intent(in) :: in_which_locreg
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      real(kind=8),dimension(nphidim),intent(in) :: phi
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(3,norbp),intent(out) :: com

      ! Local variables
      real(kind=8),dimension(:),allocatable :: psir
      type(workarr_sumrho) :: w
      integer :: ist, istr, iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, iat, iiat, l, m, nl1, nl2, nl3
      real(kind=8),dimension(-1:1) :: dipole
      real(kind=8) :: weight, tt, x, y, z, r2, hxh, hyh, hzh, q, qtot, monopole, r
      real(kind=8),parameter :: sigma2=0.1d0

      call f_routine(id='supportfunction_centers')

      ! Transform the support functions to real space
      psir = f_malloc(max(nphirdim,1),id='psir')
      ist=1
      istr=1
      do iorb=1,norbp
          iiorb=isorb+iorb
          ilr=in_which_locreg(iiorb)
          call initialize_work_arrays_sumrho(lzd%Llr(ilr),.true.,w)
          call daub_to_isf(lzd%Llr(ilr), w, phi(ist), psir(istr))
          call deallocate_work_arrays_sumrho(w)
          !write(*,'(a,4i8,es16.6)') 'INITIAL: iproc, iiorb, n, istr, ddot', &
          !    iproc, iiorb, lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, &
          !    istr, ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1)
          !testarr(1,iiorb) = ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1) 
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=nphirdim+1) then
          call f_err_throw('ERROR on process '//adjustl(trim(yaml_toa(bigdft_mpi%iproc)))//': istr/=nphirdim+1', &
               err_name='BIGDFT_RUNTIME_ERROR')
          stop
      end if

      hxh = 0.5d0*lzd%hgrids(1)
      hyh = 0.5d0*lzd%hgrids(2)
      hzh = 0.5d0*lzd%hgrids(3)

      istr = 1
      do iorb=1,norbp
          iiorb=isorb+iorb
          ilr=in_which_locreg(iiorb)
          call geocode_buffers(lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)
          !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
          com(1:3,iorb) = 0.d0
          weight = 0.d0
          do i3=1,lzd%llr(ilr)%d%n3i
              ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
              z = ii3*hzh
              do i2=1,lzd%llr(ilr)%d%n2i
                  ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
                  y = ii2*hyh
                  do i1=1,lzd%llr(ilr)%d%n1i
                      ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                      x = ii1*hxh
                      tt = psir(istr)**2
                      com(1,iorb) = com(1,iorb) + x*tt
                      com(2,iorb) = com(2,iorb) + y*tt
                      com(3,iorb) = com(3,iorb) + z*tt
                      weight = weight + tt
                      istr = istr + 1
                  end do
              end do
          end do
          !call yaml_map('weight',weight)
          com(1:3,iorb) = com(1:3,iorb)/weight

      end do

      call f_free(psir)

      call f_release_routine()

    end subroutine supportfunction_centers







  subroutine add_penalty_term(geocode, nfvctr, neighbor, rxyz, cell_dim, com, alpha, n, ovrlp, ham)
    use module_base
    implicit none
 
    ! Calling arguments
    character(len=1),intent(in) :: geocode
    integer,intent(in) :: nfvctr, n
    logical,dimension(nfvctr),intent(in) :: neighbor
    real(kind=8),dimension(3),intent(in) :: rxyz, cell_dim
    real(kind=8),intent(in) :: alpha
    real(kind=8),dimension(3,nfvctr),intent(in) :: com
    real(kind=8),dimension(n,n),intent(inout) :: ovrlp
    real(kind=8),dimension(n,n),intent(inout) :: ham

    ! Local variables
    logical :: perx, pery, perz
    integer :: is1, ie1, is2, ie2, is3, ie3, icheck, ii, i, jj, j, i1, i2, i3
    real(kind=8) :: rr2, x, y, z, ttx, tty, ttz, tt
 
    call f_routine(id='add_penalty_term')
 
    ! Determine the periodicity...
    !write(*,*) 'geocode',geocode
    perx=(geocode /= 'F')
    pery=(geocode == 'P')
    perz=(geocode /= 'F')
    if (perx) then
        is1 = -1
        ie1 = 1
    else
        is1 = 0
        ie1 = 0
    end if
    if (pery) then
        is2 = -1
        ie2 = 1
    else
        is2 = 0
        ie2 = 0
    end if
    if (perz) then
        is3 = -1
        ie3 = 1
    else
        is3 = 0
        ie3 = 0
    end if
 
 
    ! Add the penalty term
    icheck = 0
    ii = 0
    do i=1,nfvctr
        if (neighbor(i)) then
            jj = 0
            do j=1,nfvctr
                if (neighbor(j)) then
                    icheck = icheck + 1
                    jj = jj + 1
                    if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                    if (i==j) then
                        rr2 = huge(rr2)
                        do i3=is3,ie3
                            z = rxyz(3) + i3*cell_dim(3)
                            ttz = (com(3,i)-z)**2
                            do i2=is2,ie2
                                y = rxyz(2) + i2*cell_dim(2)
                                tty = (com(2,i)-y)**2
                                do i1=is1,ie1
                                    x = rxyz(1) + i1*cell_dim(1)
                                    ttx = (com(1,i)-x)**2
                                    tt = ttx + tty + ttz
                                    if (tt<rr2) then
                                        rr2 = tt
                                    end if
                                end do
                            end do
                        end do
                        !write(*,*) 'i, j, ii, jj, tt', ii, jj, alpha*rr2**3
                        !ham(jj,ii) = ham(jj,ii) + alpha*rr2**3*ovrlp(jj,ii)
                        if (jj==ii) then
                            ham(jj,ii) = ham(jj,ii) + alpha*rr2**3
                        end if
                    end if
                end if
            end do
        end if
    end do
    if (icheck>n**2) then
        call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
            &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
    end if

    call f_release_routine()

  end subroutine add_penalty_term




  subroutine add_penalty_term_new(geocode, nat, nfvctr, neighbor, rxyz, on_which_atom, &
             multipoles, cell_dim, com, alpha, n, ham, nmax, penalty_matrix)
    use module_base
    use multipole_base, only: lmax
    use module_base
    implicit none
 
    ! Calling arguments
    character(len=1),intent(in) :: geocode
    integer,intent(in) :: nat, nfvctr, n, nmax
    logical,dimension(nfvctr),intent(in) :: neighbor
    real(kind=8),dimension(3),intent(in) :: rxyz, cell_dim
    integer,dimension(nfvctr),intent(in) :: on_which_atom
    real(kind=8),dimension(-lmax:lmax,0:lmax,nfvctr) :: multipoles
    real(kind=8),dimension(3,nfvctr),intent(in) :: com
    real(kind=8),intent(in) :: alpha
    real(kind=8),dimension(n,n),intent(inout) :: ham
    real(kind=8),dimension(n,n),intent(out) :: penalty_matrix

    ! Local variables
    logical :: perx, pery, perz
    integer :: is1, ie1, is2, ie2, is3, ie3, icheck, ii, i, jj, j, i1, i2, i3
    integer :: il, jl, im, jm, ix, iy, iz, iat, jat
    real(kind=8) :: rr2, x, y, z, ttx, tty, ttz, tt, argi, argj, expi, expj, silim, sjljm, rr, xx, yy, zz
    real(kind=8),dimension(1:3) :: rip, rjp
    integer,parameter :: nn=25
    real(kind=8),parameter :: hh=0.35d0
    real(kind=8),parameter :: sigma2 = 1.d-1


    stop 'not working any more'
 
!!    call f_routine(id='add_penalty_term_new')
!!
!!    call f_zero(penalty_matrix)
!! 
!!    ! Determine the periodicity...
!!    !write(*,*) 'geocode',geocode
!!    perx=(geocode /= 'F')
!!    pery=(geocode == 'P')
!!    perz=(geocode /= 'F')
!!    if (perx) then
!!        is1 = -1
!!        ie1 = 1
!!    else
!!        is1 = 0
!!        ie1 = 0
!!    end if
!!    if (pery) then
!!        is2 = -1
!!        ie2 = 1
!!    else
!!        is2 = 0
!!        ie2 = 0
!!    end if
!!    if (perz) then
!!        is3 = -1
!!        ie3 = 1
!!    else
!!        is3 = 0
!!        ie3 = 0
!!    end if
!! 
!!    ! FOR THE MOMENT NOT WORKING FOR PERIODIC SYSTEMS! NEED TO TAKE THIS INTO ACCOUNT.
!!
!!
!!    ! Loop over all elements to be calculated
!!    icheck = 0
!!    ii = 0
!!    do i=1,nfvctr
!!        if (neighbor(i)) then
!!            iat = on_which_atom(i)
!!            jj = 0
!!            do j=1,nfvctr
!!                if (neighbor(j)) then
!!                    !!write(*,*) 'i, j', i, j
!!                    jat = on_which_atom(j)
!!                    icheck = icheck + 1
!!                    jj = jj + 1
!!                    if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
!!                    !if (i==j) then
!!                        ! distances from the atoms iat and jat (respectively the sup fun centered on them) to the one for which the projector should be calculated
!!                        rip(1:3) = com(1:3,i) - rxyz(1:3)
!!                        rjp(1:3) = com(1:3,j) - rxyz(1:3)
!!                        ! Do the summation over l,l' and m,m'
!!                        tt = 0.d0
!!                        do il=0,lmax
!!                            do jl=0,lmax
!!                                do im=-il,il
!!                                    do jm=-jl,jl
!!                                        ! Do the integration
!!                                        rr = 0.d0
!!                                        do ix=-nn,nn
!!                                            x = real(ix,kind=8)*hh
!!                                            xx = x + rxyz(1)
!!                                            do iy=-nn,nn
!!                                                y = real(iy,kind=8)*hh
!!                                                yy = y + rxyz(2)
!!                                                do iz=-nn,nn
!!                                                    z = real(iz,kind=8)*hh
!!                                                    zz = z + rxyz(3)
!!                                                    argi = ((x-rip(1))**2 + (y-rip(2))**2 + (z-rip(3))**2)/(2*sigma2)
!!                                                    argj = ((x-rjp(1))**2 + (y-rjp(2))**2 + (z-rjp(3))**2)/(2*sigma2)
!!                                                    !expi = safe_exp(-argi)/(2*pi*sigma2)**(3.d0/2.d0)
!!                                                    !expj = safe_exp(-argj)/(2*pi*sigma2)**(3.d0/2.d0)
!!                                                    expi = safe_exp(-argi)/(1.d0*pi*sigma2)**(3.d0/4.d0)
!!                                                    expj = safe_exp(-argj)/(1.d0*pi*sigma2)**(3.d0/4.d0)
!!                                                    silim = spherical_harmonic(il,im,xx,yy,zz)*sqrt(4.d0*pi)
!!                                                    sjljm = spherical_harmonic(jl,jm,xx,yy,zz)*sqrt(4.d0*pi)
!!                                                    !if (abs(argi)>1000.d0) write(*,*) 'WARNING argi'
!!                                                    !if (abs(argj)>1000.d0) write(*,*) 'WARNING argj'
!!                                                    !if (abs(expi)>1000.d0) write(*,*) 'WARNING expi'
!!                                                    !if (abs(expj)>1000.d0) write(*,*) 'WARNING expj'
!!                                                    !if (abs(silim)>1000.d0) write(*,*) 'WARNING silim'
!!                                                    !if (abs(sjljm)>1000.d0) write(*,*) 'WARNING sjljm'
!!                                                    !rr = rr + silim*expi*alpha*(x**2+y**2+z**2)*sjljm*expj*hh**3
!!                                                    !rr = rr + silim*expi*sjljm*expj*hh**3*sqrt(4.d0*pi)
!!                                                    !write(*,*) 'argi, expi, argj, expj', argi, argj, expi, expj
!!                                                    rr = rr + silim*expi*alpha*(x**2+y**2+z**2)**3*sjljm*expj*hh**3
!!                                                    !rr = rr + silim*expi*sjljm*expj*hh**3
!!                                                end do
!!                                            end do
!!                                        end do
!!                                        !write(*,*) 'i, j, il, im, jl, jm, rr', i, j, il, im, jl, jm, rr
!!                                        !tt = tt + multipoles(im,il,iat)*multipoles(jm,jl,jat)*rr
!!                                        !if (il==0 .and. jl==0) then
!!                                            tt = tt + multipoles(im,il,i)*multipoles(jm,jl,j)*rr
!!                                        !end if
!!                                        !if (abs(multipoles(im,il,iat))>1000.d0) write(*,*) 'WARNING multipoles(im,il,iat)'
!!                                        !if (abs(multipoles(jm,jl,jat))>1000.d0) write(*,*) 'WARNING multipoles(jm,jl,jat)'
!!                                    end do
!!                                end do
!!                            end do
!!                        end do
!!                        write(*,*) 'i, j, ii, jj, tt', ii, jj, tt
!!                        ham(jj,ii) = ham(jj,ii) + tt
!!                        penalty_matrix(jj,ii) = tt
!!                    !end if
!!                end if
!!            end do
!!        end if
!!    end do
!!    
!!    if (icheck>n**2) then
!!        call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
!!            &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
!!    end if
!!
!!    call f_release_routine()

  end subroutine add_penalty_term_new




  subroutine order_eigenvalues(n, eigenvalues, ids)
    use module_base
    use sort, only: QsortC
    implicit none

    ! Calling arguments
    integer,intent(in) :: n
    real(kind=8),dimension(n),intent(inout) :: eigenvalues
    integer,dimension(n),intent(inout) :: ids

    ! Local variables
    integer :: i, ind, ii
    real(kind=8) :: tt
    integer,dimension(:),allocatable :: lookup
    integer,dimension(:),allocatable :: ids_tmp

    call f_routine(id='order_eigenvalues')

    !! Order the eigenvalues and IDs
    !do i=1,n
    !    ! add i-1 since we are only searching in the subarray
    !    ind = minloc(eigenvalues(i:n),1) + (i-1)
    !    tt = eigenvalues(i)
    !    eigenvalues(i) = eigenvalues(ind)
    !    eigenvalues(ind) = tt
    !    ii = ids(i)
    !    ids(i) = ids(ind)
    !    ids(ind) = ii
    !end do

    !do i=1,n
    !    write(200+bigdft_mpi%iproc,*) eigenvalues(i), ids(i)
    !end do

    lookup = f_malloc(n,id='lookup')
    do i=1,n
        lookup(i) = i
    end do
    call QsortC(eigenvalues, lookup)
    ids_tmp = f_malloc(n,id='ids_tmp')
    call f_memcpy(src=ids, dest=ids_tmp)
    do i=1,n
        ind = lookup(i)
        ids(i) = ids_tmp(ind)
    end do

    !do i=1,n
    !    write(300+bigdft_mpi%iproc,*) eigenvalues(i), ids(i)
    !end do

    call f_free(lookup)
    call f_free(ids_tmp)

    call f_release_routine()

 end subroutine order_eigenvalues


 subroutine calculate_projector(n, ntot, nmax, kkat, ids, evals, coeff, occ_all, proj)
   use module_base
   implicit none

   ! Calling arguments
   integer :: n, ntot, nmax, kkat
   integer,dimension(ntot),intent(in) :: ids
   real(kind=8),dimension(ntot),intent(in) :: evals, occ_all
   real(kind=8),dimension(nmax,nmax),intent(in) :: coeff
   !real(kind=8),intent(in) :: kT, ef
   real(kind=8),dimension(n,n),intent(out) :: proj

   ! Local variables
   integer :: ij, ieval, i, j
   real(kind=8) :: occ

   call f_routine(id='calculate_projector')

   ij = 0
   do ieval=1,ntot
       if (ids(ieval)/=kkat) cycle
       ij = ij + 1
       !occ = 1.d0/(1.d0+safe_exp( (evals(ieval)-ef)*(1.d0/kT) ) )
       occ = occ_all(ieval)
       do i=1,n
           do j=1,n
               proj(j,i) = proj(j,i) + occ*coeff(j,ij)*coeff(i,ij)
           end do
      end do
   end do

   call f_release_routine()

 end subroutine calculate_projector


 subroutine unitary_test_multipoles(iproc, nproc, nphi, nphir, orbs, lzd, smmd, smat, collcom, hgrids)
   use module_base
   use module_types, only: orbitals_data, comms_linear, local_zone_descriptors, comms_linear
   use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_TASKGROUP, assignment(=), &
                                matrices_null, sparsematrix_malloc_ptr, deallocate_matrices, &
                                sparse_matrix_metadata
   use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
   use yaml_output
   use bounds, only: geocode_buffers
   use orbitalbasis
   implicit none
   ! Calling arguments
   integer,intent(in) :: iproc, nproc, nphi, nphir
   type(orbitals_data),intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: lzd
   type(sparse_matrix_metadata),intent(in) :: smmd
   type(sparse_matrix),intent(in) :: smat
   type(comms_linear),intent(in) :: collcom
   real(kind=8),dimension(3) :: hgrids
   ! Local variables
   integer :: iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, l, m, i, ind, ist, istr, ii, nl1, nl2, nl3
   real(kind=8) :: x, y, z, r2, r, factor, rmax, factor_normalization, val,sigma
   real(kind=8),dimension(:),allocatable :: phi2r, phi2, phi1r, phi1
   real(kind=8),dimension(:,:),allocatable :: locregcenter
   type(matrices) :: multipole_matrix
   type(workarr_sumrho) :: w
   real(kind=8),dimension(-lmax:lmax,0:lmax) :: errors
   real(kind=8),dimension(-lmax:lmax,0:lmax) :: values_orig
   real(kind=8),dimension(-lmax:lmax,0:lmax) :: values
   type(orbital_basis) :: psi_ob
   real(gp), dimension(3) :: acell
   real(wp), dimension(:,:,:), allocatable :: Qlm
   real(kind=8),dimension(:),allocatable :: gg1, gg2, gg3


   call f_routine(id='unitary_test_multipoles')

   if (iproc==0) then
       call yaml_comment('Unitary test of the multipole routines',hfill='~')
   end if

   multipole_matrix = matrices_null()
   multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(smat, SPARSE_TASKGROUP, id='multipole_matrix%matrix_compr')

   phi2r = f_malloc0(nphir,id='phi2r')
   phi1r = f_malloc(nphir,id='phi1r')
   phi1 = f_malloc0(nphi,id='phi1')

   locregcenter = f_malloc0((/3,lzd%nlr/),id='locregcenter')

   do ilr=1,lzd%nlr 
       if (lzd%Llr(ilr)%geocode/='F') then
           call f_err_throw('support function locregs must always have free BC')
       end if
   end do
   call geocode_buffers('F', lzd%glr%geocode, nl1, nl2, nl3)

  sigma=0.5d0
  r=1.d0 !not used...

   ist = 0
   do iorb=1,orbs%norbp
       iiorb = orbs%isorb + iorb
       ilr = orbs%inwhichlocreg(iiorb)
       !rmax = min(lzd%llr(ilr)%d%n1i*0.25d0*hgrids(1),lzd%llr(ilr)%d%n2i*0.25d0*hgrids(2),lzd%llr(ilr)%d%n3i*0.25d0*hgrids(3))
       rmax = min(lzd%llr(ilr)%d%n1*0.5d0*hgrids(1),lzd%llr(ilr)%d%n2*0.5d0*hgrids(2),lzd%llr(ilr)%d%n3*0.5d0*hgrids(3))
       factor_normalization = 0.5d0*lzd%hgrids(1)*0.5d0*lzd%hgrids(2)*0.5d0*lzd%hgrids(3) !*3.d0/(4.d0*pi*rmax**3)
       ! Since the radial function is constant and thus not decaying towards the boundaries of the integration sphere, the center
       ! of the integration volume must be on a gridpoint to avoid truncation artifacts.
       locregcenter(1:3,ilr) = get_closest_gridpoint(lzd%llr(ilr)%locregcenter,hgrids)

       gg1 = f_malloc(lzd%llr(ilr)%d%n1i,id='gg1')
       gg2 = f_malloc(lzd%llr(ilr)%d%n2i,id='gg2')
       gg3 = f_malloc(lzd%llr(ilr)%d%n3i,id='gg3')
       do i1=1,lzd%llr(ilr)%d%n1i
           ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
           x = ii1*0.5d0*lzd%hgrids(1) - locregcenter(1,ilr)
           gg1(i1) = safe_exp(-0.5d0*x**2/sigma**2)
       end do
       do i2=1,lzd%llr(ilr)%d%n2i
           ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
           y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2,ilr)
           gg2(i2) = safe_exp(-0.5d0*y**2/sigma**2)
       end do
       do i3=1,lzd%llr(ilr)%d%n3i
           ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
           z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3,ilr)
           gg3(i3) = safe_exp(-0.5d0*z**2/sigma**2)
       end do

       !$omp parallel default(none) &
       !$omp shared(lzd, nl1, nl2, nl3, ilr, locregcenter, sigma, phi2r, ist, factor_normalization) &
       !$omp shared(gg1, gg2, gg3, r) &
       !$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, ind, l, m, factor)
       !$omp do schedule(static)
       do i3=1,lzd%llr(ilr)%d%n3i
           ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
           z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3,ilr)
           do i2=1,lzd%llr(ilr)%d%n2i
               ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
               y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2,ilr)
               do i1=1,lzd%llr(ilr)%d%n1i
                   ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                   x = ii1*0.5d0*lzd%hgrids(1) - locregcenter(1,ilr)
                   ind = (i3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2-1)*lzd%llr(ilr)%d%n1i + i1
                   do l=0,lmax
                       do m=-l,l
                           factor = get_test_factor(l,m)*factor_normalization*sqrt(4.d0*pi*real(2*l+1,kind=8))/sigma**3
                           if (l==1) then
                              factor = factor/(3.d0*sigma**2)
                           else if (l==2) then
                              factor = factor/(15.d0*sigma**4)
                           end if
                           phi2r(ist+ind) = phi2r(ist+ind) + &
                                gg1(i1)*gg2(i2)*gg3(i3)*factor*solid_harmonic(0, r, l, m , x, y, z)/sqrt(twopi**3)
                       end do
                   end do
               end do
           end do
       end do
       !$omp end do
       !$omp end parallel
       call f_free(gg1)
       call f_free(gg2)
       call f_free(gg3)
       ist = ist + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
    end do

   if (nproc>1) then
       call mpiallred(locregcenter, mpi_sum, comm=bigdft_mpi%mpi_comm)
   end if

   ! Transform back to wavelets
   phi2 = f_malloc0(nphi,id='phi2')
   ist=1
   istr=1
   do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      call initialize_work_arrays_sumrho(lzd%llr(ilr),.true.,w)
      call isf_to_daub(lzd%llr(ilr), w, phi2r(istr), phi2(ist))
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
      istr = istr + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
   end do

   !alternative solution, less memory, less operations, less communications
   acell(1)=0.5_gp*hgrids(1)*Lzd%glr%d%n1i
   acell(2)=0.5_gp*hgrids(2)*Lzd%glr%d%n2i
   acell(3)=0.5_gp*hgrids(3)*Lzd%glr%d%n3i
   Qlm=f_malloc([-lmax .to. lmax ,0 .to. lmax,1 .to. orbs%norbp ],id='Qlm')
   call orbital_basis_associate(psi_ob,orbs=orbs,phis_wvl=phi2,Lzd=Lzd,id='unitary_test_multipoles')
   call Qlm_phi(lmax,smmd%geocode,hgrids,acell,psi_ob,Qlm,.false.,centers=locregcenter)
   call orbital_basis_release(psi_ob)
   call f_zero(values)
   do l=0,lmax
      do m=-l,l 
         val = 0.d0
         do iorb=1,orbs%norbp
            val = val + Qlm(m,l,iorb)
         end do
         values(m,l) = val/real(orbs%norb,kind=8)
      end do
   end do
   call f_free(Qlm)

  if (nproc > 1) call mpiallred(values,op=MPI_SUM,comm=bigdft_mpi%mpi_comm)
  do l=0,lmax
     do m=-l,l !to copy also zeros
        errors(m,l) = 100.d0*abs(values(m,l)/get_test_factor(l,m)-1.d0)
        values_orig(m,l) = get_test_factor(l,m)
     end do
  end do



!!$   ! Set phi1 to 1
!!$   phi1r(:) = 1.d0
!!$
!!$   ! Transform back to wavelets
!!$   phi2 = f_malloc0(nphi,id='phi2')
!!$   ist=1
!!$   istr=1
!!$   do iorb=1,orbs%norbp
!!$       iiorb=orbs%isorb+iorb
!!$       ilr=orbs%inwhichlocreg(iiorb)
!!$       call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.true.,w)
!!$       call isf_to_daub(lzd%llr(ilr), w, phi2r(istr), phi2(ist))
!!$       call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.false.,w)
!!$       call isf_to_daub(lzd%llr(ilr), w, phi1r(istr), phi1(ist))
!!$       call deallocate_work_arrays_sumrho(w)
!!$       ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!$       istr = istr + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
!!$   end do
!!$
!!$
!!$   !do ind=1,nphi
!!$   !    write(*,*) 'ind, val', ind, phi2(ind)
!!$   !end do
!!$
!!$
!!$
!!$   do l=0,lmax
!!$       do m=-l,l
!!$           call calculate_multipole_matrix(iproc, nproc, l, m, nphi, phi1, phi2, nphir, hgrids, &
!!$                orbs, collcom, lzd, smat, locregcenter, 'sphere', multipole_matrix) !==> values
!!$           val = 0.d0
!!$           do iorb=1,orbs%norb
!!$               iiorb = modulo(iorb-1,smat%nfvctr)+1
!!$               ind = matrixindex_in_compressed(smat, iorb, iorb)
!!$               val = val + multipole_matrix%matrix_compr(ind)
!!$               !write(*,*) 'l, m, iorb, ind, val', &
!!$               !    l, m, iorb, ind, multipole_matrix%matrix_compr(ind)
!!$           end do
!!$           values(m,l) = val/real(orbs%norb,kind=8)
!!$           errors(m,l) = 100.d0*abs(values(m,l)/get_test_factor(l,m)-1.d0)
!!$           values_orig(m,l) = get_test_factor(l,m)
!!$           !if (iproc==0) write(*,*) 'l, m, val, error', l, m, val, abs(val-get_test_factor(l,m))
!!$       end do
!!$   end do

   call f_free(locregcenter)

   if (iproc==0) then
       call yaml_mapping_open('Unitary check of the multipole calculations')
       call yaml_sequence_open('Original values')
       do l=0,lmax
           call yaml_sequence(advance='no')
           call yaml_map('q'//adjustl(trim(yaml_toa(l))),values_orig(-l:l,l),fmt='(1es16.8)')
       end do
       call yaml_mapping_close()
       call yaml_sequence_open('Calculated values')
       do l=0,lmax
           call yaml_sequence(advance='no')
           call yaml_map('q'//adjustl(trim(yaml_toa(l))),values(-l:l,l),fmt='(1es16.8)')
       end do
       call yaml_mapping_close()
       call yaml_sequence_open('Relative errors in percent')
       do l=0,lmax
           call yaml_sequence(advance='no')
           call yaml_map('q'//adjustl(trim(yaml_toa(l))),errors(-l:l,l),fmt='(1es10.2)')
       end do
       call yaml_sequence_close()
       call yaml_mapping_close()
   end if



   call f_free(phi1r)
   call f_free(phi2r)
   call f_free(phi1)
   call f_free(phi2)
   call deallocate_matrices(multipole_matrix)

   !call multipole_analysis_driver(iproc, nproc, 2, tmb%npsidim_orbs, tmb%psi, &
   !     max(tmb%collcom_sr%ndimpsi_c,1), at, tmb%lzd%hgrids, &
   !     tmb%orbs, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, tmb%collcom, tmb%lzd, &
   !     tmb%orthpar, tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%kernel_, rxyz, &
   !     method='projector')

   call f_release_routine()

 end subroutine unitary_test_multipoles


 !> Calculate the closes grid point to a given point r
 function get_closest_gridpoint(r,hgrids) result(cgp)
   implicit none
   ! Calling arguments
   real(kind=8),dimension(3),intent(in) :: r
   real(kind=8),dimension(3),intent(in) :: hgrids
   real(kind=8),dimension(3) :: cgp
   ! Local variables
   integer :: i, ii
   real(kind=8) :: tt

   do i=1,3
       tt = r(i)/hgrids(1)
       ii = nint(tt)
       cgp(i) = real(ii,kind=8)*hgrids(i)
   end do
 end function get_closest_gridpoint


 !> SM: similar to support_function_multipoles. This one calculates the "gross" multipoles (i.e. without taking into account the "core" contribution)
 subroutine support_function_gross_multipoles(iproc, nproc, tmb, atoms, shift, denspot)
   use module_base
   use module_types
   use locreg_operations
   use yaml_output
   use multipole_base, only: lmax
   use bounds, only: geocode_buffers
   use sparsematrix_base, only: matrices, matrices_null, sparsematrix_malloc_ptr, SPARSE_TASKGROUP, assignment(=), &
                                deallocate_matrices
   use orthonormalization, only: orthonormalizelocalized

   use communications_base, only: TRANSPOSE_FULL
   use transposed_operations, only: calculate_overlap_transposed
   use communications, only: transpose_localized
   use multipole_base, only: external_potential_descriptors, external_potential_descriptors_null, &
                             multipole_set_null, multipole_null, deallocate_external_potential_descriptors
   use orbitalbasis
   ! Calling arguments
   integer,intent(in) :: iproc, nproc
   type(DFT_wavefunction),intent(inout) :: tmb
   type(atoms_data),intent(in) :: atoms
   real(kind=8),dimension(3),intent(in) :: shift !< global shift of the atomic positions
   type(DFT_local_fields), intent(inout) :: denspot
 
   integer :: ist, istr, iorb, iiorb, ilr, i, iat, iter, itype, mm
   integer :: i1, i2, i3, ii1, ii2, ii3, nl1, nl2, nl3, ii, l, m, ind, iat_old, methTransformOverlap
   real(kind=8),dimension(:),allocatable :: rmax, phi1, phi1r, phi_ortho
   real(kind=8),dimension(:,:),allocatable :: delta_centers
   real(kind=8),dimension(:,:),allocatable :: center_locreg, center_orb
   real(kind=8),dimension(:),allocatable :: phir, phir_one
   real(kind=8) :: hxh, hyh, hzh, tt, x, y, z, weight, factor
   type(workarr_sumrho) :: w
   character(len=20) :: atomname
   character(len=20),dimension(:),allocatable :: names
   integer,dimension(:),allocatable :: iatype_tmp
   type(matrices) :: multipole_matrix
   real(kind=8),dimension(:,:,:),allocatable :: multipoles
   real(kind=8),dimension(:),allocatable :: scaled
   real(kind=8),dimension(:),pointer :: phit_c, phit_f
   logical :: can_use_transposed
   type(external_potential_descriptors) :: ep
   character(len=*),parameter :: no='none', onsite='on-site'
   character(len=*),parameter :: do_ortho = onsite
   type(orbital_basis) :: psi_ob
   real(gp), dimension(3) :: acell, center
   real(wp), dimension(:,:,:), allocatable :: Qlm

   call f_routine(id='support_function_gross_multipoles')

   phi_ortho = f_malloc(size(tmb%psi),id='phi_ortho')
   call f_memcpy(src=tmb%psi, dest=phi_ortho)
   if (do_ortho == no) then
       ! Do nothing
   else if (do_ortho == onsite) then
       phit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='phit_c')
       phit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='phit_f')
       methTransformOverlap = -2
       can_use_transposed = .false.
       call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
            1.d-8, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
            tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, &
            phi_ortho, phit_c, phit_f, &
            can_use_transposed)
       !!!@ TEST @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       !!call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
       !!     TRANSPOSE_FULL, phi_ortho, phit_c, phit_f, tmb%lzd)
       !!multipole_matrix = matrices_null()
       !!multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%s, SPARSE_FULL, id='multipole_matrix%matrix_compr')
       !!call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
       !!     phit_c, phit_c, phit_f, phit_f, tmb%linmat%s, multipole_matrix)
       !!if (iproc==0) then
       !!    do i=1,size(multipole_matrix%matrix_compr)
       !!        write(*,*) 'i, mat', i, multipole_matrix%matrix_compr(i)
       !!    end do
       !!end if
       !!call deallocate_matrices(multipole_matrix)
       !!!@ END TEST @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       call f_free_ptr(phit_c)
       call f_free_ptr(phit_f)
       ! END TEST ############################################################
   else
       call f_err_throw('wrong orthonormalisation method',err_name='BIGDFT_RUNTIME_ERROR')
   end if

 
   rmax = f_malloc0(tmb%orbs%norb,id='rmax')
   phir = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir')
   phir_one = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir_one')
   phir_one = 1.d0
 
   !call to_zero(3*tmb%orbs%norb, dipole_net(1,1))
   !call to_zero(9*tmb%orbs%norb, quadropole_net(1,1,1))


   center_locreg = f_malloc0((/3,tmb%lzd%nlr/),id='center_locreg')
   center_orb = f_malloc0((/3,tmb%lzd%nlr/),id='center_orb')
   multipole_matrix = matrices_null()
   multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%s, SPARSE_TASKGROUP, id='multipole_matrix%matrix_compr')

  ! Set phi1 to 1
  phi1r = f_malloc(max(tmb%collcom_sr%ndimpsi_c,1),id='phi1r')
  phi1 = f_malloc0(tmb%npsidim_orbs,id='phi1')

  multipoles = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.tmb%orbs%norb/),id='multipoles')

  scaled = f_malloc0(tmb%orbs%norb,id='scaled')

  phi1r(:) = 1.d0


  call f_zero(multipoles)
 
  ist=1
  istr=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iat=tmb%orbs%onwhichatom(iiorb)
      call initialize_work_arrays_sumrho(tmb%lzd%Llr(ilr),.true.,w)
      ! Transform the support function to real space
      call daub_to_isf(tmb%lzd%llr(ilr), w, phi_ortho(ist), phir(istr))
      call initialize_work_arrays_sumrho(tmb%lzd%llr(ilr),.false.,w)
      ! Transform the functions which is constantly one to wavelets
      call isf_to_daub(tmb%lzd%llr(ilr), w, phi1r(istr), phi1(ist))
      call deallocate_work_arrays_sumrho(w)

      call calculate_weight_center(tmb%lzd%llr(ilr), tmb%lzd%glr, tmb%lzd%hgrids, &
           phir(istr), center_locreg(1:3,ilr), center_orb(1:3,iiorb))
      ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i

!!$      mesh=cell_new(tmb%lzd%Llr(ilr)%geocode,&
!!$           [tmb%lzd%llr(ilr)%d%n1i,tmb%lzd%llr(ilr)%d%n2i,tmb%lzd%llr(ilr)%d%n3i],0.5_gp*tmb%lzd%hgrids)
!!$      boxit=box_iterator(mesh,origin=0.5d0*hgrids*[tmb%lzd%llr(ilr)%nsi1,tmb%lzd%llr(ilr)%nsi2,tmb%lzd%llr(ilr)%nsi3])
!!$      weight = 0.d0
!!$      do while(box_next_point(boxit))
!!$         tt = phir(boxit%ind+istr)**2
!!$         center_locreg(1,ilr) = center_locreg(1,ilr) + boxit%rxyz(1)*tt
!!$         center_locreg(2,ilr) = center_locreg(2,ilr) + boxit%rxyz(2)*tt
!!$         center_locreg(3,ilr) = center_locreg(3,ilr) + boxit%rxyz(3)*tt
!!$         weight = weight + tt
!!$      end do
  end do



  if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
      call f_err_throw('istr/=tmb%collcom_sr%ndimpsi_c+1')
      !write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
      !stop
  end if

  if (nproc>1) then
      call mpiallred(center_locreg, mpi_sum, comm=bigdft_mpi%mpi_comm)
      call mpiallred(center_orb, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

  hxh = 0.5d0*tmb%lzd%hgrids(1)
  hyh = 0.5d0*tmb%lzd%hgrids(2)
  hzh = 0.5d0*tmb%lzd%hgrids(3)
  factor = hxh*hyh*hzh

  !alternative solution, less memory, less operations, less communications
  acell(1)=0.5_gp*tmb%lzd%hgrids(1)*tmb%Lzd%glr%d%n1i
  acell(2)=0.5_gp*tmb%lzd%hgrids(2)*tmb%Lzd%glr%d%n2i
  acell(3)=0.5_gp*tmb%lzd%hgrids(3)*tmb%Lzd%glr%d%n3i
  Qlm=f_malloc([-lmax .to. lmax ,0 .to. lmax,1 .to. tmb%orbs%norbp ],id='Qlm')
  call orbital_basis_associate(psi_ob,orbs=tmb%orbs,&
       phis_wvl=phi_ortho,Lzd=tmb%Lzd,id='support_function_gross_multipoles')
  call Qlm_phi(lmax,tmb%linmat%smmd%geocode,tmb%lzd%hgrids,acell,psi_ob,Qlm,.false.,centers=center_locreg)
  call orbital_basis_release(psi_ob)
  do iorb=1,tmb%orbs%norbp
     iiorb = tmb%orbs%isorb + iorb
     do l=0,lmax
        do m=-lmax,lmax !to copy also zeros
           multipoles(m,l,iiorb) = Qlm(m,l,iorb)*factor
        end do
     end do
  end do
  call f_free(Qlm)

!!$  do l=0,lmax
!!$      do m=-l,l
!!$          call f_zero(multipole_matrix%matrix_compr)
!!$          ! Calculate the multipole matrix
!!$          call calculate_multipole_matrix(iproc, nproc, l, m, tmb%npsidim_orbs, phi1, phi_ortho, &
!!$               max(tmb%collcom_sr%ndimpsi_c,1), tmb%lzd%hgrids, &
!!$               tmb%orbs, tmb%collcom, tmb%lzd, tmb%linmat%s, center_locreg, 'box', multipole_matrix)! =>>multipoles
!!$          !write(*,*) 'multipole_matrix%matrix_compr(1)',multipole_matrix%matrix_compr(1)
!!$          ! Take the diagonal elements and scale by factor (anyway there is no really physical meaning in the actual numbers)
!!$          do iorb=1,tmb%orbs%norbp
!!$              iiorb = tmb%orbs%isorb + iorb
!!$              ind = matrixindex_in_compressed(tmb%linmat%s, iiorb, iiorb)
!!$              multipoles(m,l,iiorb) = multipole_matrix%matrix_compr(ind)*factor
!!$              !write(*,*) 'iorb, multipoles(:,:,iiorb)',iorb, multipoles(:,:,iiorb)
!!$          end do
!!$      end do
!!$  end do

  ! Normalize the multipoles such that the largest component has the magnitude 1
  do iorb=1,tmb%orbs%norbp
      iiorb = tmb%orbs%isorb + iorb
      tt = maxval(abs(multipoles(:,:,iiorb)))
      !write(*,*) 'iorb, tt', iorb, tt
      multipoles(:,:,iiorb) = multipoles(:,:,iiorb)/tt
      scaled(iiorb) = tt
  end do
 
 
  if (bigdft_mpi%nproc>1) then
      call mpiallred(multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
      call mpiallred(scaled, mpi_sum, comm=bigdft_mpi%mpi_comm)
      call mpiallred(rmax, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

 
  if (iproc==0) then
      call yaml_sequence_open('Gross support functions moments')
      call yaml_map('Orthonormalization',do_ortho)
      iatype_tmp = f_malloc(tmb%orbs%norb,id='iatype_tmp')
      delta_centers = f_malloc((/3,tmb%orbs%norb/),id='delta_centers')
      iat_old = -1
      names = f_malloc_str(len(names),tmb%orbs%norb,id='names')

      ep = external_potential_descriptors_null()
      ep%nmpl = tmb%orbs%norb
      allocate(ep%mpl(ep%nmpl))

      do iorb=1,tmb%orbs%norb
          iat = tmb%orbs%onwhichatom(iorb)
          if (iat/=iat_old) then
              ii = 1
          else
              ii = ii + 1
          end if
          iat_old = iat
          ilr = tmb%orbs%inwhichlocreg(iorb)
          itype = atoms%astruct%iatype(iat)
          iatype_tmp(iorb) = itype
          names(iorb) = trim(atoms%astruct%atomnames(itype))//'-'//adjustl(trim(yaml_toa(ii)))
          ! delta_centers gives the difference between the charge center and the localization center
          delta_centers(1:3,iorb) = center_locreg(1:3,ilr) - tmb%lzd%llr(ilr)%locregcenter(1:3)
          !write(*,*) 'iorb, ilr, center_locreg(1:3,ilr) - tmb%lzd%llr(ilr)%locregcenter(1:3)', &
          !           iorb, ilr, center_locreg(1:3,ilr) - tmb%lzd%llr(ilr)%locregcenter(1:3)
          !write(*,*) 'iorb, delta_centers(1:3,iorb)', iorb, delta_centers(1:3,iorb)
          ! Undo the global shift of the centers
          center_orb(1:3,iorb) = center_orb(1:3,iorb) + shift(1:3)

          ep%mpl(iorb) = multipole_set_null()
          allocate(ep%mpl(iorb)%qlm(0:lmax))
          ep%mpl(iorb)%rxyz = center_orb(1:3,iorb)
          ep%mpl(iorb)%sym = trim(names(iorb))
          do l=0,lmax
              ep%mpl(iorb)%qlm(l) = multipole_null()
              !if (l>=3) cycle
              ep%mpl(iorb)%qlm(l)%q = f_malloc_ptr(2*l+1,id='q')
              mm = 0
              do m=-l,l
                  mm = mm + 1
                  ep%mpl(iorb)%qlm(l)%q(mm) = multipoles(m,l,iorb)
              end do
          end do
      end do
      call write_multipoles_new(ep, lmax, atoms%astruct%units, &
           delta_centers, tmb%orbs%onwhichatom, scaled)
      call deallocate_external_potential_descriptors(ep)
      call f_free(delta_centers)
      call f_free(iatype_tmp)
      call f_free_str(len(names),names)
      call f_free(scaled)
      call yaml_sequence_close()
  end if

 
  call f_free(rmax)
  call f_free(phir)
  call f_free(phi1r)
  call f_free(phi1)
  call f_free(phir_one)
  call deallocate_matrices(multipole_matrix)
  call f_free(center_locreg)
  call f_free(center_orb)
  call f_free(multipoles)
  call f_free(phi_ortho)

  call f_release_routine()
 
 end subroutine support_function_gross_multipoles


 !!subroutine get_optimal_sigmas(iproc, nproc, nsigma, collcom_sr, smatl, kernel_, at, lzd, ep, shift, rxyz, ixc, denspot)
 !!  use module_base
 !!  use module_types, only: DFT_wavefunction, input_variables, DFT_local_fields, comms_linear, DFT_local_fields, &
 !!                          local_zone_descriptors
 !!  use sparsematrix_base, only: sparse_matrix, matrices
 !!  use module_atoms, only: atoms_data
 !!  use Poisson_Solver, only: H_potential
 !!  use rhopotential, only: sumrho_for_TMBs, corrections_for_negative_charge
 !!  use yaml_output
 !!  implicit none
 !!  ! Calling arguments
 !!  integer,intent(in) :: iproc, nproc, nsigma, ixc
 !!  type(comms_linear),intent(inout) :: collcom_sr
 !!  type(sparse_matrix),intent(in) :: smatl
 !!  type(matrices),intent(in) :: kernel_
 !!  type(atoms_data),intent(in) :: at
 !!  type(local_zone_descriptors),intent(in) :: lzd
 !!  type(external_potential_descriptors),intent(in) :: ep
 !!  real(kind=8),dimension(3),intent(in) :: shift
 !!  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
 !!  type(DFT_local_fields),intent(inout) :: denspot
 !!  ! Local variables
 !!  real(kind=8),dimension(:,:,:,:),allocatable :: test_pot
 !!  logical :: rho_negative, exists, found, all_norms_ok
 !!  real(kind=8) :: ehart_ps, diff, tt, diff_min, diff_dipole, diff_dipole_min, rdim
 !!  !integer,parameter :: nsigma=3
 !!  real(kind=8),parameter :: step=0.20d0
 !!  integer :: i1, i2, i3, isigma0, isigma1, isigma2, impl, l
 !!  integer :: nzatom, nelpsp, npspcode, itype, ioffset, ishift
 !!  real(gp),dimension(0:4,0:6) :: psppar
 !!  real(kind=8),dimension(:,:),allocatable :: sigmax
 !!  real(kind=8),dimension(0:lmax) :: factor, factorx, factor_min
 !!  real(kind=8),dimension(:),allocatable :: rhov_orig
 !!  real(kind=8),dimension(3) :: dipole_exact, dipole_trial
 !!  real(kind=8) :: rloc
 !!  integer,dimension(:),allocatable :: psp_source

 !!  call f_routine(id='get_optimal_sigmas')

 !!  if (iproc==0) call yaml_comment('Determine optimal sigmas for the radial Gaussians',hfill='~')

 !!  test_pot = f_malloc0((/size(denspot%V_ext,1),size(denspot%V_ext,2),size(denspot%V_ext,3),2/),id='test_pot')
 !!  rhov_orig = f_malloc(size(denspot%rhov),id='rhov_orig')

 !!  ! Keep the original value fo rhov, which contains the entire potential
 !!  call f_memcpy(src=denspot%rhov, dest=rhov_orig)

 !!  ! Calculate the correct electrostatic potential, i.e. electronic plus ionic part
 !!  call sumrho_for_TMBs(iproc, nproc, lzd%hgrids(1), lzd%hgrids(2), lzd%hgrids(3), &
 !!       collcom_sr, smatl, kernel_, &
 !!       denspot%dpbox%ndimrhopot, &
 !!       denspot%rhov, rho_negative)
 !!  if (rho_negative) then
 !!      call corrections_for_negative_charge(iproc, nproc, at, denspot)
 !!  end if

 !!  denspot%rho_work = f_malloc_ptr(denspot%dpbox%ndimrhopot,id='denspot%rho_work')
 !!  ioffset=lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%i3xcsh
 !!  if (denspot%dpbox%ndimrhopot>0) then
 !!      call vcopy(denspot%dpbox%ndimpot,denspot%rhov(ioffset+1),1,denspot%rho_work(1),1)
 !!      ! add the spin down part if present
 !!      if (denspot%dpbox%nrhodim==2) then
 !!          ishift=denspot%dpbox%ndimrhopot/denspot%dpbox%nrhodim !start of the spin down part
 !!          call axpy(denspot%dpbox%ndimpot, 1.d0, &
 !!                    denspot%rhov(ioffset+ishift+1), &
 !!                    1, denspot%rho_work(1),1)
 !!      end if
 !!  end if
 !!  !do l=1,size(denspot%rho_work)
 !!  !    write(100,*) denspot%rho_work(l)
 !!  !end do
 !!  !write(*,*) 'calculate dipole with rho_work'
 !!  call calculate_dipole_moment(denspot%dpbox, 1, at, rxyz, denspot%rho_work, &
 !!       calculate_quadrupole=.true., dipole=dipole_exact, quiet_=.true.)
 !!  call f_free_ptr(denspot%rho_work)

 !!  call H_potential('D',denspot%pkernel,denspot%rhov,denspot%V_ext,ehart_ps,0.0_dp,.true.,&
 !!       quiet=denspot%PSquiet)!,rho_ion=denspot%rho_ion)
 !!  call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
 !!       denspot%rhov(1), 1, test_pot(1,1,1,1), 1)

 !!  ! Get an initial guess for the sigmas (use rloc from the pseudopotential)
 !!  sigmax = f_malloc((/0.to.lmax,1.to.ep%nmpl/),id='sigmax')
 !!  psp_source = f_malloc(ep%nmpl,id='psp_source')
 !!  do impl=1,ep%nmpl
 !!      !ixc = 1
 !!      !if (iproc==0) then
 !!      !    call yaml_warning('WARNING: USE ixc = 1 IN GET_OPTIMAL_SIGMAS')
 !!      !end if
 !!      !call psp_from_data(ep%mpl(impl)%sym, nzatom, nelpsp, npspcode, ixc, psppar, exists)
 !!      !if (.not.exists) then
 !!      !    call f_err_throw('No PSP available for external multipole type '//trim(ep%mpl(impl)%sym), &
 !!      !         err_name='BIGDFT_INPUT_VARIABLES_ERROR')
 !!      !end if
 !!      !ep%mpl(impl)%sigma(0:lmax) = psppar(0,0)-min(0.9d0,step*real(nsigma/2,kind=8))*psppar(0,0)
 !!      !sigmax(0:lmax,impl) = psppar(0,0)
 !!      !!found = .false.
 !!      !!search_loop: do itype=1,at%astruct%ntypes
 !!      !!    if (trim(ep%mpl(impl)%sym)==trim(at%astruct%atomnames(itype))) then
 !!      !!        sigmax(0:lmax,impl) = 1.0d0*at%psppar(0,0,itype)
 !!      !!        found = .true.
 !!      !!        exit search_loop
 !!      !!    end if
 !!      !!end do search_loop
 !!      !!if (.not.found) then
 !!      !!    call f_err_throw('No PSP available for external multipole type '//trim(ep%mpl(impl)%sym), &
 !!      !!         err_name='BIGDFT_INPUT_VARIABLES_ERROR')
 !!      !!end if
 !!      call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp, psp_source(impl), rloc)
 !!      sigmax(0:lmax,impl) = 1.d0*rloc
 !!  end do
 !!  if (iproc==0) call write_psp_source(ep, psp_source)
 !!  call f_free(psp_source)

 !!  if (iproc==0) then
 !!      call yaml_sequence_open('Determine optimal sigmas')
 !!  end if
 !!  factorx(0:lmax) = max(0.1d0,1.d0-step*real(nsigma/2,kind=8))
 !!  ! The following loops are designed for lmax=2... stop otherwise
 !!  if (lmax>2) then
 !!      call f_err_throw('the maximal lmax possible is 2, but here we have '//trim(yaml_toa(lmax)),&
 !!           err_name='BIGDFT_RUNTIME_ERROR')
 !!  end if
 !!  diff_min = huge(diff_min)
 !!  diff_dipole_min = huge(diff_dipole_min)
 !!  factor_min(0:lmax) = 1.d0 !initialization
 !!  do isigma2=1,nsigma
 !!      do isigma1=1,nsigma
 !!          do isigma0=1,nsigma
 !!              factor(0) = factorx(0) + real(isigma0-1,kind=8)*step
 !!              factor(1) = factorx(1) + real(isigma1-1,kind=8)*step
 !!              factor(2) = factorx(2) + real(isigma2-1,kind=8)*step
 !!              do impl=1,ep%nmpl
 !!                  !ep%mpl(impl)%sigma(l) = ep%mpl(impl)%sigma(l) + step
 !!                  ep%mpl(impl)%sigma(0:lmax) = sigmax(0:lmax,impl)*factor(0:lmax)
 !!                  !if (iproc==0) write(*,*) 'impl, sigma', impl, ep%mpl(impl)%sigma(0:lmax)
 !!              end do
 !!              call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
 !!                   denspot%V_ext(1,1,1,1), 1, test_pot(1,1,1,2), 1)
 !!              call potential_from_charge_multipoles(iproc, nproc, at, denspot, ep, 1, &
 !!                   denspot%dpbox%ndims(1), 1, denspot%dpbox%ndims(2), &
 !!                   denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
 !!                   denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
 !!                   denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
 !!                   denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3), &
 !!                   shift, verbosity=0, ixc=ixc, lzd=lzd, pot=test_pot(:,:,:,2), &
 !!                   rxyz=rxyz, dipole_total=dipole_trial, all_norms_ok=all_norms_ok)

 !!              if (all_norms_ok) then
 !!                  diff_dipole = (dipole_exact(1)-dipole_trial(1))**2 + &
 !!                                (dipole_exact(2)-dipole_trial(2))**2 + &
 !!                                (dipole_exact(3)-dipole_trial(3))**2
 !!                  rdim = 1.d0/(real(size(denspot%V_ext,1),kind=8)*&
 !!                               real(size(denspot%V_ext,1),kind=8)*&
 !!                               real(size(denspot%V_ext,1),kind=8))
 !!                  diff = 0.d0
 !!                  do i3=1,size(denspot%V_ext,3)
 !!                      do i2=1,size(denspot%V_ext,2)
 !!                          do i1=1,size(denspot%V_ext,1)
 !!                              !write(800,*) 'i1, i2, i3, vals', i1, i2, i3, test_pot(i1,i2,i3,1), test_pot(i1,i2,i3,2)
 !!                              diff = diff + rdim*(test_pot(i1,i2,i3,1)-test_pot(i1,i2,i3,2))**2
 !!                          end do
 !!                      end do
 !!                  end do
 !!                  call mpiallred(diff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
 !!                  !!tt = diff/(real(size(denspot%V_ext,1),kind=8)*&
 !!                  !!           real(size(denspot%V_ext,1),kind=8)*&
 !!                  !!           real(size(denspot%V_ext,1),kind=8))
 !!              end if
 !!              if (iproc==0) then
 !!                  call yaml_sequence(advance='no')
 !!                  call yaml_mapping_open(flow=.true.)
 !!                  call yaml_map('rloc mult',factor,fmt='(f3.1)')
 !!                  call yaml_map('Gaussian norms ok',all_norms_ok)
 !!                  if (all_norms_ok) then
 !!                      call yaml_map('dipole norm diff (actual/min)',(/diff_dipole,diff_dipole_min/),fmt='(es9.3)')
 !!                      call yaml_map('avg pot diff (actual/min)',(/diff,diff_min/),fmt='(es9.3)')
 !!                  end if
 !!                  call yaml_mapping_close()
 !!              end if
 !!              if (all_norms_ok) then
 !!                  if (diff<diff_min) then
 !!                      !factor_min(0:lmax) = factor(0:lmax)
 !!                      diff_min = diff
 !!                  end if
 !!                  if (diff_dipole<diff_dipole_min) then
 !!                      factor_min(0:lmax) = factor(0:lmax)
 !!                      diff_dipole_min = diff_dipole
 !!                  end if
 !!              end if
 !!          end do
 !!      end do
 !!  end do
 !!  if (iproc==0) then
 !!      call yaml_sequence_close()
 !!  end if
 !!  if (iproc==0) call yaml_map('optimal sigma multiplication factors',factor_min,fmt='(f4.2)')
 !!  do impl=1,ep%nmpl
 !!      ep%mpl(impl)%sigma(0:lmax) = sigmax(0:lmax,impl)*factor_min(0:lmax)
 !!  end do

 !!  call f_memcpy(src=rhov_orig, dest=denspot%rhov)

 !!  call f_free(sigmax)
 !!  call f_free(test_pot)
 !!  call f_free(rhov_orig)

 !!  call f_release_routine()

 !!end subroutine get_optimal_sigmas

 subroutine calculate_gaussian(is, ie, idim, nl, nglob, periodic, hh, shift, ep, gaussian_array)
   use module_base
   use multipole_base, only: lmax, external_potential_descriptors
   implicit none

   ! Calling arguments
   integer,intent(in) :: is, ie, idim, nl, nglob
   logical,intent(in) :: periodic
   real(kind=8),intent(in) :: hh
   real(kind=8),dimension(3),intent(in) :: shift
   type(external_potential_descriptors),intent(in) :: ep
   real(kind=8),dimension(0:lmax,is:ie,ep%nmpl),intent(out) :: gaussian_array

   ! Local variables
   integer :: i, ii, impl, l, isx, iex, n, imod, nn, nu, nd, js, je, j
   real(kind=8) :: x, tt, sig, dr

   call f_routine(id='calculate_gaussian')

   !do impl=1,ep%nmpl
   !    write(*,*) 'idim, shift(idim), ep%mpl(impl)%rxyz(idim)', idim, shift(idim), ep%mpl(impl)%rxyz(idim)
   !end do

   call f_zero(gaussian_array)

   ! Calculate the boundaries of the Gaussian to be calculated. To make it simple, take always the maximum:
   ! - free BC: entire box
   ! - periodic BC: half of the box size, with periodic wrap around
   if (.not.periodic) then
       js = 0
       je = 0
   else
       js = -1
       je = 1
   end if

   !$omp parallel default(none) &
   !$omp shared(is, ie, hh, shift, idim, ep, gaussian_array, js, je, nl, nglob) &
   !$omp private(i, ii, x, impl, tt, l, sig, j, dr)
   !$omp do
   do impl=1,ep%nmpl
       do i=is,ie
           ii = i - nl - 1
           tt = huge(tt)
           do j=js,je
               dr = real(ii+j*nglob,kind=8)*hh + shift(idim) - ep%mpl(impl)%rxyz(idim)
               if (abs(dr)<abs(tt)) tt = dr
           end do
           tt = tt**2
           do l=0,lmax
               sig = ep%mpl(impl)%sigma(l)
               gaussian_array(l,i,impl) = gaussian(sig,tt)
           end do
       end do
   end do
   !$omp end do
   !$omp end parallel

   call f_release_routine()

   contains

     function gaussian(sigma, r2) result(g)
       use module_base, only: pi => pi_param
       implicit none
       ! Calling arguments
       real(kind=8),intent(in) :: sigma, r2
       real(kind=8) :: tt, g

       ! Only calculate the Gaussian if the result will be larger than 10^-30
       tt = r2/(2.d0*sigma**2)
       if (tt<=69.07755279d0) then
           g = safe_exp(-tt)
           g = g/sqrt(2.d0*pi*sigma**2)**1!3
       else
           g = 0.d0
       end if
       !g = g/(sigma**3*sqrt(2.d0*pi)**3)

     end function gaussian

 end subroutine calculate_gaussian

 subroutine calculate_norm(nproc, is1, ie1, is2, ie2, is3, ie3, ep, &
            hhh, gaussians1, gaussians2, gaussians3, norm)
   use module_base
   use multipole_base, only: external_potential_descriptors
   implicit none
   
   ! Calling arguments
   integer,intent(in) :: nproc, is1, ie1, is2, ie2, is3, ie3
   type(external_potential_descriptors),intent(in) :: ep
   real(kind=8),intent(in) :: hhh
   real(kind=8),dimension(0:lmax,is1:ie1,1:ep%nmpl),intent(in) :: gaussians1
   real(kind=8),dimension(0:lmax,is2:ie2,1:ep%nmpl),intent(in) :: gaussians2
   real(kind=8),dimension(0:lmax,is3:ie3,1:ep%nmpl),intent(in) :: gaussians3
   real(kind=8),dimension(0:2,ep%nmpl),intent(out) :: norm

   ! Local variables
   integer :: impl, i1, i2, i3, ii1, ii2, ii3, l
   real(kind=8) :: gg
   real(kind=8),dimension(0:lmax) :: gg23

   call f_routine(id='calculate_norm')

   call f_zero(norm)

   !$omp parallel default(none) &
   !$omp shared(ep, is1, ie1, is2, ie2, is3, ie3, norm, hhh) &
   !$omp shared(gaussians1, gaussians2, gaussians3) &
   !$omp private(impl, i1, i2, i3, ii1, ii2, ii3, l, gg23, gg) 
   !$omp do schedule(guided)
   do impl=1,ep%nmpl
       i3loop: do i3=is3,ie3
           if (maxval(gaussians3(:,i3,impl))<1.d-20) cycle i3loop
           ii3 = i3 - 15
           i2loop: do i2=is2,ie2
               if (maxval(gaussians2(:,i2,impl))<1.d-20) cycle i2loop
               ii2 = i2 - 15
               do l=0,lmax
                   gg23(l) = gaussians2(l,i2,impl)*gaussians3(l,i3,impl)
               end do
               i1loop: do i1=is1,ie1
                   if (maxval(gaussians1(:,i1,impl))<1.d-20) cycle i1loop
                   ii1 = i1 - 15
                   do l=0,lmax
                       ! Calculate the Gaussian as product of three 1D Gaussians
                       !gg = gaussians1(l,i1,impl)*gaussians2(l,i2,impl)*gaussians3(l,i3,impl)
                       gg = gaussians1(l,i1,impl)*gg23(l)
                       norm(l,impl) = norm(l,impl) + gg*hhh
                   end do
               end do i1loop
           end do i2loop
       end do i3loop
   end do
   !$omp end do
   !$omp end parallel

   ! Sum up the norms of the Gaussians.
   if (nproc>1) then
       call mpiallred(norm, mpi_sum, comm=bigdft_mpi%mpi_comm)
   end if

   call f_release_routine()

 end subroutine calculate_norm
 

!> Calculate the dipole of a Field given in the rho array.
!! The parallel distribution used is the one of the potential
subroutine calculate_dipole_moment(dpbox,nspin,at,rxyz,rho,calculate_quadrupole,dipole,quadrupole,quiet_)
  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  use yaml_output
  use box
  use numerics
  implicit none
  integer, intent(in) :: nspin
  type(denspot_distribution), intent(inout) :: dpbox
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(dpbox%mesh%ndims(1),dpbox%mesh%ndims(2),max(dpbox%n3p, 1),nspin), target, intent(in) :: rho
  !!!!logical :: is_net_charge !< true if the charge density is already the net charge (i.e. including the compensating core charge)
  logical,intent(in) :: calculate_quadrupole
  real(kind=8),dimension(3),intent(out),optional :: dipole
  real(kind=8),dimension(3,3),intent(out),optional :: quadrupole
  logical,intent(in),optional :: quiet_

!  integer :: ierr,n3p,nc1,nc2,nc3, nnc3, ii3, i3shift
  real(gp) :: q,qtot, delta_term,x,y,z,ri,rj,tt
  !integer  :: i1,i2,i3, nl1,nl2,nl3, n1i,n2i,n3i,  is, ie
  integer :: iat,ispin,i, j
  real(gp), dimension(3) :: dipole_el,dipole_cores,tmpdip,charge_center_cores
  real(gp),dimension(3,nspin) :: charge_center_elec
  real(gp), dimension(3,3) :: quadropole_el,quadropole_cores,tmpquadrop
  real(dp), dimension(:,:,:,:), pointer :: ele_rho
  logical :: quiet
!!$  real(dp), dimension(:,:,:,:), pointer :: rho_buf

  call f_routine(id='calculate_dipole_moment')

  if (present(quiet_)) then
      quiet = quiet_
  else
      quiet = .false.
  end if
  
!!$  n1i=dpbox%mesh%ndims(1)
!!$  n2i=dpbox%mesh%ndims(2)
!!$  n3i=dpbox%mesh%ndims(3)
!!$  n3p=dpbox%n3p
!!$
!!$
!!$  if (at%astruct%geocode /= 'F') then
!!$     nl1=1
!!$     !nl3=1
!!$     nc1=n1i
!!$     !nc3=n3i
!!$     nc3=n3p
!!$     nnc3=n3i
!!$     !is = 1
!!$     is = dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+1
!!$     ie = dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+dpbox%nscatterarr(dpbox%mpi_env%iproc,2)
!!$     i3shift = 1
!!$  else
!!$     nl1=15
!!$     !nl3=15
!!$     !nl3=max(1,15-dpbox%nscatterarr(dpbox%mpi_env%iproc,3))
!!$     nc1=n1i-31
!!$     !nc3=n3i-31
!!$     !nc3=n3p-31
!!$     is = max(dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+1,15)
!!$     ie = min(dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+dpbox%nscatterarr(dpbox%mpi_env%iproc,2),n3i-17)
!!$     nnc3=n3i-31
!!$     i3shift = 15
!!$     !write(*,*) 'iproc, is, ie, nl3, nc3, n3p', bigdft_mpi%iproc, is, ie, nl3, nc3, n3p
!!$  end if
!!$  nc3 = ie - is + 1 !number of z planes to be treated
!!$  nl3=max(1,i3shift-dpbox%nscatterarr(dpbox%mpi_env%iproc,3)) !offset within rho array
!!$  !value of the buffer in the y direction
!!$  if (at%astruct%geocode == 'P') then
!!$     nl2=1
!!$     nc2=n2i
!!$  else
!!$     nl2=15
!!$     nc2=n2i-31
!!$  end if

  qtot=0.d0
  call f_zero(dipole_cores)!(1:3)=0._gp
  call f_zero(charge_center_cores)
  do iat=1,at%astruct%nat
     !write(*,*) 'iat, rxyz(1:3,iat)',iat, rxyz(1:3,iat)
     q=at%nelpsp(at%astruct%iatype(iat))
     dipole_cores(1:3)=dipole_cores(1:3)+q * rxyz(1:3,iat)
     qtot=qtot+q
  end do
  !this defines the origin of the coordinate system
  if (qtot /=0.0_gp) charge_center_cores=dipole_cores/qtot

  !!write(*,*) 'dipole_cores',dipole_cores
  !!write(*,*) 'nc3',nc3

  !calculate electronic dipole and thus total dipole of the system
  call f_zero(dipole_el)!   (1:3)=0._gp
  do ispin=1,nspin
     !the iterator here is on the potential distribution
     do while(box_next_point(dpbox%bitp))
        q= - rho(dpbox%bitp%i,dpbox%bitp%j,dpbox%bitp%k-dpbox%bitp%i3s+1,ispin) *dpbox%mesh%volume_element
        !write(*,*) 'i1, i2, i3, nl1, nl2, nl3, q', i1, i2, i3, nl1, nl2, nl3, q
        qtot=qtot+q
        dipole_el=dipole_el+q*(dpbox%bitp%rxyz-charge_center_cores)
     end do
!!$     do i3=0,nc3 - 1
!!$        !ii3 = i3 + dpbox%nscatterarr(dpbox%mpi_env%iproc,3)
!!$        ii3 = i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3) - i3shift !real coordinate, without buffer
!!$        !write(*,*) 'iproc, i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3), ii3', &
!!$        !            bigdft_mpi%iproc, i3+nl3+dpbox%nscatterarr(dpbox%mpienv%iproc,3), ii3
!!$        do i2=0,nc2 - 1
!!$           do i1=0,nc1 - 1
!!$              !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
!!$              !q= ( ele_rho(ind,ispin) ) * hxh*hyh*hzh 
!!$              !q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
!!$              q= - rho(i1+nl1,i2+nl2,i3+nl3,ispin) *dpbox%mesh%volume_element
!!$              !write(*,*) 'i1, i2, i3, nl1, nl2, nl3, q', i1, i2, i3, nl1, nl2, nl3, q
!!$              qtot=qtot+q
!!$              dipole_el(1)=dipole_el(1)+ q* at%astruct%cell_dim(1)/real(nc1,dp)*i1 
!!$              dipole_el(2)=dipole_el(2)+ q* at%astruct%cell_dim(2)/real(nc2,dp)*i2
!!$              dipole_el(3)=dipole_el(3)+ q* at%astruct%cell_dim(3)/real(nnc3,dp)*ii3
!!$           end do
!!$        end do
!!$     end do
  !!write(*,*) 'iproc, dipole_el,sum(rho), qtot',bigdft_mpi%iproc,dipole_el,sum(rho), qtot
  end do

  !!call mpi_barrier(mpi_comm_world,ispin)
  call mpiallred(qtot, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
  call mpiallred(dipole_el, mpi_sum, comm=bigdft_mpi%mpi_comm)
  !!call mpi_barrier(mpi_comm_world,ispin)
  !!write(*,*) 'after allred: iproc, dipole_el,sum(rho), qtot',bigdft_mpi%iproc,dipole_el,sum(rho), qtot

  !!write(*,*) 'dipole_cores first', dipole_cores
  !!call mpi_barrier(mpi_comm_world,ispin)

  !quadrupole should be calculated with the shifted positions!
  quadrupole_if: if (calculate_quadrupole) then

      call f_zero(quadropole_cores)!(1:3,1:3)=0._gp
      do iat=1,at%astruct%nat
         q=at%nelpsp(at%astruct%iatype(iat))
         tmpdip=rxyz(:,iat)-charge_center_cores
         tt=square(dpbox%mesh,tmpdip)
          do i=1,3
             ri=rxyz(i,iat)-charge_center_cores(i)
             do j=1,3
                rj=rxyz(j,iat)-charge_center_cores(j)
                if (i==j) then
                   delta_term = tt
                else
                   delta_term=0.d0
                end if
                quadropole_cores(j,i) = quadropole_cores(j,i) + q*(3.d0*rj*ri-delta_term)
             end do
          end do
       end do

      ! charge center
!!$      charge_center_elec(1:3,1:nspin)=0.d0
!!$      do ispin=1,nspin
!!$         !LG: this is exactly the same calculation as before
!!$         qtot=0.d0
!!$         do while(box_next_point(dpbox%bitp))
!!$            q= - rho(dpbox%bitp%i,dpbox%bitp%j,dpbox%bitp%k-dpbox%bitp%i3s+1,ispin) *dpbox%mesh%volume_element
!!$            !write(*,*) 'i1, i2, i3, nl1, nl2, nl3, q', i1, i2, i3, nl1, nl2, nl3, q
!!$            qtot=qtot+q
!!$            charge_center_elec=charge_center_elec+q*dpbox%bitp%rxyz
!!$         end do

!!$          do i3=0,nc3 - 1
!!$             ii3 = i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3) - i3shift !real coordinate, without buffer
!!$              do i2=0,nc2 - 1
!!$                  do i1=0,nc1 - 1
!!$                      !q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
!!$                      q= - rho(i1+nl1,i2+nl2,i3+nl3,ispin) * dpbox%mesh%volume_element
!!$                      x=at%astruct%cell_dim(1)/real(nc1,dp)*i1
!!$                      y=at%astruct%cell_dim(2)/real(nc2,dp)*i2
!!$                      z=at%astruct%cell_dim(3)/real(nnc3,dp)*ii3
!!$                      charge_center_elec(1,ispin) = charge_center_elec(1,ispin) + q*x
!!$                      charge_center_elec(2,ispin) = charge_center_elec(2,ispin) + q*y
!!$                      charge_center_elec(3,ispin) = charge_center_elec(3,ispin) + q*z
!!$                      qtot=qtot+q
!!$                  end do
!!$              end do
!!$      end do
!!$          !!write(*,*) 'qtot',qtot
!!$         call mpiallred(qtot, 1, mpi_sum, comm=bigdft_mpi%mpi_comm) !LG: why two spins parallelized that way?
!!$         charge_center_elec(1:3,ispin)=charge_center_elec(1:3,ispin)/qtot
!!$      end do
!!$
!!$      call mpiallred(charge_center_elec, mpi_sum, comm=bigdft_mpi%mpi_comm)

       call f_zero(quadropole_el)
      do ispin=1,nspin
         do while(box_next_point(dpbox%bitp))
            q= - rho(dpbox%bitp%i,dpbox%bitp%j,dpbox%bitp%k-dpbox%bitp%i3s+1,ispin) *dpbox%mesh%volume_element
            tmpdip=dpbox%bitp%rxyz-charge_center_cores
            tt=square(dpbox%mesh,tmpdip)
            do i=1,3
               ri=dpbox%bitp%rxyz(i)-charge_center_cores(i)
               do j=1,3
                  rj=dpbox%bitp%rxyz(j)-charge_center_cores(j)
                  if (i==j) then
                     delta_term = tt
                  else
                     delta_term=0.d0
                  end if
                  quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
               end do
            end do
         end do
         
!!$          do i3=0,nc3 - 1
!!$             ii3 = i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3) - i3shift !real coordinate, without buffer
!!$              do i2=0,nc2 - 1
!!$                  do i1=0,nc1 - 1
!!$                      !q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
!!$                      q= - rho(i1+nl1,i2+nl2,i3+nl3,ispin) * dpbox%mesh%volume_element
!!$                      x=at%astruct%cell_dim(1)/real(nc1,dp)*i1
!!$                      y=at%astruct%cell_dim(2)/real(nc2,dp)*i2
!!$                      z=at%astruct%cell_dim(3)/real(nnc3,dp)*ii3
!!$                      do i=1,3
!!$                          select case (i)
!!$                          case (1)
!!$                              !ri=x-charge_center_cores(1)
!!$                              ri=x+(charge_center_cores(1)-charge_center_elec(1,ispin))
!!$                          case (2)
!!$                              !ri=y-charge_center_cores(2)
!!$                              ri=y+(charge_center_cores(2)-charge_center_elec(2,ispin))
!!$                          case (3)
!!$                              !ri=z-charge_center_cores(3)
!!$                              ri=z+(charge_center_cores(3)-charge_center_elec(3,ispin))
!!$                          case default
!!$                              stop 'wrong value of i'
!!$                          end select
!!$                          do j=1,3
!!$                              select case (j)
!!$                              case (1)
!!$                                  !rj=x-charge_center_cores(1)
!!$                                  rj=x+(charge_center_cores(1)-charge_center_elec(1,ispin))
!!$                              case (2)
!!$                                  !rj=y-charge_center_cores(2)
!!$                                  rj=y+(charge_center_cores(2)-charge_center_elec(2,ispin))
!!$                              case (3)
!!$                                  !rj=z-charge_center_cores(3)
!!$                                  rj=z+(charge_center_cores(3)-charge_center_elec(3,ispin))
!!$                              case default
!!$                                  stop 'wrong value of j'
!!$                              end select
!!$                              if (i==j) then
!!$                                  !delta_term = (x-charge_center_cores(1))**2 + &
!!$                                  !             (y-charge_center_cores(2))**2 + &
!!$                                  !             (z-charge_center_cores(3))**2
!!$                                  delta_term = (x+(charge_center_cores(1)-charge_center_elec(1,ispin)))**2 + &
!!$                                               (y+(charge_center_cores(2)-charge_center_elec(2,ispin)))**2 + &
!!$                                               (z+(charge_center_cores(3)-charge_center_elec(3,ispin)))**2
!!$                              else
!!$                                  delta_term=0.d0
!!$                              end if
!!$                              quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
!!$                          end do
!!$                      end do
!!$                  end do
!!$              end do
!!$          end do
      end do

      !!if (.not. is_net_charge) then
          call mpiallred(quadropole_el, mpi_sum, comm=bigdft_mpi%mpi_comm)
          tmpquadrop=quadropole_cores+quadropole_el
      !!else
      !!    tmpquadrop=quadropole_el
      !!end if

      if (present(quadrupole)) then
          quadrupole = tmpquadrop
      end if

      !!if (dpbox%mpi_env%nproc > 1) then
      !!   call f_free_ptr(ele_rho)
      !!else
      !!   nullify(ele_rho)
      !!end if

  end if quadrupole_if

  !!if (.not.is_net_charge) then
      tmpdip=dipole_el !dipole_cores+ !should not be needed as it is now included in the center of charge
  !!else
  !!    tmpdip=dipole_el
  !!end if
  !!write(*,*) 'tmpdip before',tmpdip
  !!call mpi_barrier(mpi_comm_world,ispin)
  !!write(*,*) 'tmpdip',tmpdip
  if (present(dipole)) dipole(1:3) = tmpdip(1:3)
  if(bigdft_mpi%iproc==0 .and. .not.quiet) then
     call yaml_map('Multipole analysis origin',charge_center_cores,fmt='(1pe14.6)')
     call yaml_mapping_open('Electric Dipole Moment (AU)')
       call yaml_map('P vector',tmpdip(1:3),fmt='(1pe13.4)')
       call yaml_map('norm(P)',sqrt(sum(tmpdip**2)),fmt='(1pe14.6)')
     call yaml_mapping_close()
     tmpdip=tmpdip/Debye_AU  ! au2debye
     call yaml_mapping_open('Electric Dipole Moment (Debye)')
       call yaml_map('P vector',tmpdip(1:3),fmt='(1pe13.4)')
       call yaml_map('norm(P)',sqrt(sum(tmpdip**2)),fmt='(1pe14.6)')
     call yaml_mapping_close()


      if (calculate_quadrupole) then
          call yaml_mapping_open('Quadrupole Moment (AU)')
            call yaml_map('Q matrix',tmpquadrop,fmt='(1pe13.4)')
           call yaml_map('trace',tmpquadrop(1,1)+tmpquadrop(2,2)+tmpquadrop(3,3),fmt='(es12.2)')
          call yaml_mapping_close()
      end if

  end if

  call f_release_routine()

END SUBROUTINE calculate_dipole_moment


subroutine calculate_rpowerx_matrices(iproc, nproc, nphi, nphir, lzd, orbs, collcom, phi, smat, rpower_matrix)
  use module_base
  use module_types, only: local_zone_descriptors, orbitals_data, comms_linear
  use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use transposed_operations, only: calculate_overlap_transposed
  use sparsematrix_base, only: sparse_matrix, matrices
  use bounds, only: geocode_buffers
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nphi, nphir
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  type(comms_linear),intent(in) :: collcom
  real(kind=8),dimension(nphi),intent(in) :: phi
  type(sparse_matrix),intent(in) :: smat
  type(matrices),dimension(24),intent(inout) :: rpower_matrix
  
  ! Local variables
  integer :: iorb, iiorb, ilr, iat, ii, i1, i2, i3, ii1, ii2, ii3, ist, istr, nl1, nl2, nl3, i
  type(workarr_sumrho) :: w
  real(kind=8),dimension(:),allocatable :: phir, phit_c, phit_f, xphit_c, xphit_f
  real(kind=8),dimension(:,:),allocatable :: xphi, xphir
  real(kind=8) :: hxh, hyh, hzh, x, y, z, r, r2

  call f_routine(id='calculate_rpowerx_matrices')

  xphi = f_malloc0((/nphi,24/),id='xphi')
  phir = f_malloc(nphir,id='phir')
  xphir = f_malloc0((/nphir,24/),id='xphir')

  ist=1
  istr=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      iat=orbs%onwhichatom(iiorb)
      call initialize_work_arrays_sumrho(lzd%Llr(ilr),.true.,w)
      ! Transform the support function to real space
      call daub_to_isf(lzd%llr(ilr), w, phi(ist), phir(istr))
      call initialize_work_arrays_sumrho(lzd%llr(ilr),.false.,w)

      ! NEW: CALCULATE THE WEIGHT CENTER OF THE SUPPORT FUNCTION ############################
      hxh = 0.5d0*lzd%hgrids(1)
      hyh = 0.5d0*lzd%hgrids(2)
      hzh = 0.5d0*lzd%hgrids(3)
      ii = istr
      call geocode_buffers(lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)
      do i3=1,lzd%llr(ilr)%d%n3i
          ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
          z = ii3*hzh
          do i2=1,lzd%llr(ilr)%d%n2i
              ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
              y = ii2*hyh
              do i1=1,lzd%llr(ilr)%d%n1i
                  ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                  x = ii1*hxh
                  r2 = x**2+y**2+z**2
                  xphir(ii,1) = x*phir(ii)
                  xphir(ii,2) = x**2*phir(ii)
                  xphir(ii,3) = x**3*phir(ii)
                  xphir(ii,4) = x**4*phir(ii)
                  xphir(ii,5) = y*phir(ii)
                  xphir(ii,6) = y**2*phir(ii)
                  xphir(ii,7) = y**3*phir(ii)
                  xphir(ii,8) = y**4*phir(ii)
                  xphir(ii,9) = z*phir(ii)
                  xphir(ii,10) = z**2*phir(ii)
                  xphir(ii,11) = z**3*phir(ii)
                  xphir(ii,12) = z**4*phir(ii)
                  xphir(ii,13) = x*y*phir(ii)
                  xphir(ii,14) = x**2*y*phir(ii)
                  xphir(ii,15) = x*y**2*phir(ii)
                  xphir(ii,16) = x**2*y**2*phir(ii)
                  xphir(ii,17) = x*z*phir(ii)
                  xphir(ii,18) = x**2*z*phir(ii)
                  xphir(ii,19) = x*z**2*phir(ii)
                  xphir(ii,20) = x**2*z**2*phir(ii)
                  xphir(ii,21) = y*z*phir(ii)
                  xphir(ii,22) = y**2*z*phir(ii)
                  xphir(ii,23) = y*z**2*phir(ii)
                  xphir(ii,24) = y**2*z**2*phir(ii)
                  ii = ii + 1
              end do
          end do
      end do
      ! Transform the functions back to wavelets
      do i=1,24
          call isf_to_daub(lzd%llr(ilr), w, xphir(istr,i), xphi(ist,i))
      end do
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
  end do

  ! Calculate the matrices
  phit_c = f_malloc(collcom%ndimind_c,id='phit_c')
  phit_f = f_malloc(7*collcom%ndimind_f,id='phit_f')
  xphit_c = f_malloc(collcom%ndimind_c,id='xphit_c')
  xphit_f = f_malloc(7*collcom%ndimind_f,id='xphit_f')
  call transpose_localized(iproc, nproc, nphi, orbs, collcom, &
       TRANSPOSE_FULL, phi, phit_c, phit_f, lzd)
  do i=1,24
      call transpose_localized(iproc, nproc, nphi, orbs, collcom, &
           TRANSPOSE_FULL, xphi(:,i), xphit_c, xphit_f, lzd)
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
           phit_c, xphit_c, phit_f, xphit_f, smat, rpower_matrix(i))
  end do
  !call transpose_localized(iproc, nproc, nphi, orbs, collcom, &
  !     TRANSPOSE_FULL, xphi(:,2), xphit_c, xphit_f, lzd)
  !call calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
  !     phit_c, xphit_c, xphit_f, xphit_f, smat, rpower_matrix(2))
  call f_free(phit_c)
  call f_free(phit_f)
  call f_free(xphit_c)
  call f_free(xphit_f)

  !!if (iproc==0) then
  !!    do iorb=1,smat%nvctr
  !!        write(*,*) 'i, val', iorb, rpower_matrix(1)%matrix_compr(iorb)
  !!    end do
  !!end if

  call f_free(xphi)
  call f_free(phir)
  call f_free(xphir)

  call f_release_routine()

end subroutine calculate_rpowerx_matrices


  function get_quartic_penalty(n, j, i, penaltymat, ovrlp, rxyz) result(gqp)
    implicit none

    ! Calling arguments
    integer,intent(in) :: n, j, i
    real(kind=8),dimension(n,n,24),intent(in) :: penaltymat
    real(kind=8),dimension(n,n),intent(in) :: ovrlp
    real(kind=8),dimension(3) :: rxyz
    real(kind=8) :: gqp

    gqp = penaltymat(j,i,4) - 4.d0*rxyz(1)*penaltymat(j,i,3) &
          + 6.d0*rxyz(1)**2*penaltymat(j,i,2) - 4.d0*rxyz(1)**3*penaltymat(j,i,1) &
          + rxyz(1)**4*ovrlp(j,i) &
          + penaltymat(j,i,8) - 4.d0*rxyz(2)*penaltymat(j,i,7) &
          + 6.d0*rxyz(2)**2*penaltymat(j,i,6) - 4.d0*rxyz(2)**3*penaltymat(j,i,5) &
          + rxyz(2)**4*ovrlp(j,i) &
          + penaltymat(j,i,12) - 4.d0*rxyz(3)*penaltymat(j,i,11) &
          + 6.d0*rxyz(3)**2*penaltymat(j,i,10) - 4.d0*rxyz(3)**3*penaltymat(j,i,9) &
          + rxyz(3)**4*ovrlp(j,i) &
          + 2.d0*(penaltymat(j,i,16) &
                  - 2.d0*rxyz(2)*penaltymat(j,i,14) &
                  + rxyz(2)**2*penaltymat(j,i,2) &
                  - 2.d0*rxyz(1)*penaltymat(j,i,15) &
                  + 4.d0*rxyz(1)*rxyz(2)*penaltymat(j,i,13) &
                  - 2.d0*rxyz(1)*rxyz(2)**2*penaltymat(j,i,1) &
                  + rxyz(1)**2*penaltymat(j,i,6) &
                  - 2.d0*rxyz(1)**2*rxyz(2)*penaltymat(j,i,5) &
                  + rxyz(1)**2*rxyz(2)**2*ovrlp(j,i) &
                  + penaltymat(j,i,20) &
                  - 2.d0*rxyz(3)*penaltymat(j,i,18) &
                  + rxyz(3)**2*penaltymat(j,i,2) &
                  - 2.d0*rxyz(1)*penaltymat(j,i,19) &
                  + 4.d0*rxyz(1)*rxyz(3)*penaltymat(j,i,17) &
                  - 2.d0*rxyz(1)*rxyz(3)**2*penaltymat(j,i,1) &
                  + rxyz(1)**2*penaltymat(j,i,10) &
                  - 2.d0*rxyz(1)**2*rxyz(3)*penaltymat(j,i,9) &
                  + rxyz(1)**2*rxyz(3)**2*ovrlp(j,i) &
                  + penaltymat(j,i,24) &
                  - 2.d0*rxyz(3)*penaltymat(j,i,22) &
                  + rxyz(3)**2*penaltymat(j,i,6) &
                  - 2.d0*rxyz(2)*penaltymat(j,i,23) &
                  + 4.d0*rxyz(2)*rxyz(3)*penaltymat(j,i,21) &
                  - 2.d0*rxyz(2)*rxyz(3)**2*penaltymat(j,i,5) &
                  + rxyz(2)**2*penaltymat(j,i,10) &
                  - 2.d0*rxyz(2)**2*rxyz(3)*penaltymat(j,i,9) &
                  + rxyz(2)**2*rxyz(3)**2*ovrlp(j,i) )

  end function get_quartic_penalty


  subroutine get_psp_info(sym, ixc, smmd, nelpsp, psp_source, rloc, psppar)
    use module_base
    use yaml_output
    use sparsematrix_base, only: sparse_matrix_metadata
    implicit none

    ! Calling arguments
    character(len=*),intent(in) :: sym
    integer,intent(in) :: ixc
    type(sparse_matrix_metadata),intent(in) :: smmd
    integer,intent(out) :: nelpsp, psp_source
    real(kind=8),intent(out) :: rloc
    real(kind=8),dimension(0:4,0:6,1:smmd%ntypes),intent(in),optional :: psppar

    ! Local variables
    integer :: itype, ixc_tmp, npspcode, nzatom
    logical :: found, exists
    real(gp),dimension(0:4,0:6) :: pspparx

    found = .false.
    if (present(psppar)) then
        search_loop: do itype=1,smmd%ntypes
            if (trim(sym)==trim(smmd%atomnames(itype))) then
                rloc = psppar(0,0,itype)
                nelpsp = smmd%nelpsp(itype)
                found = .true.
                psp_source = 0
                exit search_loop
            end if
        end do search_loop
    end if
    if (.not.found) then
        ixc_tmp = ixc
        call psp_from_data(trim(sym), nzatom, nelpsp, npspcode, ixc_tmp, pspparx, exists)
        if (exists) then
            rloc = pspparx(0,0)
            found = .true.
            psp_source = 1
        end if
    end if
    if (.not.found) then
        !call f_err_throw('No PSP available for external multipole type '//trim(sym), &
        !     err_name='BIGDFT_INPUT_VARIABLES_ERROR')
        call yaml_warning('No PSP available for external multipole type '//trim(sym))
        rloc = -1.d0
        nelpsp = -1
    end if
  end subroutine get_psp_info


  subroutine write_psp_source(ep, psp_source)
    use module_base
    use yaml_output
    implicit none

    ! Calling arguments
    type(external_potential_descriptors),intent(in) :: ep
    integer,dimension(ep%nmpl),intent(in) :: psp_source

    ! Local variables
    integer :: ntype, itype, impl
    logical :: written
    integer,parameter :: nmax_multipole_types = 5000
    character(len=20),dimension(nmax_multipole_types) :: multipole_type_names


    ntype = 0
    call yaml_sequence_open('Origin of the PSP data')
    do impl=1,ep%nmpl
        ! Check whether the info for this type has already been written
        written = .false.
        do itype=1,ntype
            if (trim(ep%mpl(impl)%sym)==trim(multipole_type_names(itype))) then
                written = .true.
                exit
            end if
        end do
        if (.not. written) then
            ntype = ntype + 1
            if (ntype>nmax_multipole_types) call f_err_throw('More than 5000 different multipole types are not allowed')
            multipole_type_names(ntype) = trim(ep%mpl(impl)%sym)
            call yaml_sequence(advance='no')
            if (psp_source(impl)==0) then
                call yaml_map(trim(ep%mpl(impl)%sym),'PSP of QM region')
            else if (psp_source(impl)==1) then 
                call yaml_map(trim(ep%mpl(impl)%sym),'built-in PSP')
            end if
        end if
    end do
    call yaml_sequence_close()

  end subroutine write_psp_source

    subroutine get_minmax_eigenvalues(iproc, smat, mat)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(smat%nvctr),intent(in) :: mat

      ! Local variables
      integer :: iseg, ii, i, lwork, info
      real(kind=8),dimension(:,:),allocatable :: tempmat
      real(kind=8),dimension(:),allocatable :: eval, work

      call f_routine(id='get_minmax_eigenvalues')

      tempmat = f_malloc0((/smat%nfvctr,smat%nfvctr/),id='tempmat')
      do iseg=1,smat%nseg
          ii=smat%keyv(iseg)
          do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
              tempmat(i,smat%keyg(1,2,iseg)) = mat(ii)
              ii = ii + 1
          end do
      end do
      !!if (iproc==0) then
      !!    do i=1,smat%nfvctr
      !!        do j=1,smat%nfvctr
      !!            write(*,'(a,2i6,es17.8)') 'i,j,val',i,j,tempmat(j,i)
      !!        end do
      !!    end do
      !!end if
      eval = f_malloc(smat%nfvctr,id='eval')
      lwork=100*smat%nfvctr
      work = f_malloc(lwork,id='work')
      call dsyev('n','l', smat%nfvctr, tempmat, smat%nfvctr, eval, work, lwork, info)
      if (iproc==0) write(*,*) 'eval',eval
      if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(smat%nfvctr)/),fmt='(es16.6)')

      call f_free(tempmat)
      call f_free(eval)
      call f_free(work)

      call f_release_routine()

    end subroutine get_minmax_eigenvalues


    subroutine correct_multipole_origin(nat, l, m, n, nlr, natpx, kat, kkat, &
               smmd, smats, rxyz, neighborx, perx, pery, perz, acell, &
               lower_multipole_matrices, multipole_extracted)
      use sparsematrix_base, only: sparse_matrix_metadata, sparse_matrix, matrices
      use module_types, only: orbitals_data
      implicit none

      ! Calling arguments
      integer,intent(in) :: nat, l, m, n, nlr, natpx, kat, kkat
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      type(sparse_matrix),intent(in) :: smats
      type(sparse_matrix_metadata),intent(in) :: smmd
      logical,dimension(smats%nfvctr,natpx),intent(in) :: neighborx
      logical,intent(in) :: perx, pery, perz
      real(kind=8),dimension(3),intent(in) :: acell
      type(matrices),dimension(-1:1,0:1),intent(in):: lower_multipole_matrices
      real(kind=8),dimension(n,n),intent(inout) :: multipole_extracted

      ! Local variables
      integer :: ii, ilr, i, j, iat
      real(kind=8) :: rr1, rr2, rr3
      real(kind=8),dimension(:,:,:,:),allocatable :: lmp_extracted
      real(kind=8),dimension(:,:),allocatable :: tmpmat

      call f_routine(id='correct_multipole_origin')

        if (l==1) then
            lmp_extracted = f_malloc((/1.to.n,1.to.n,0.to.0,0.to.0/),id='lmp_extracted')
            tmpmat = f_malloc((/n,n/),id='tmpmat')
            call extract_matrix(smats, lower_multipole_matrices(0,0)%matrix_compr, &
                 neighborx(1,kat), n, lmp_extracted(1,1,0,0))
            select case (m)
            case (-1)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr2 = closest_image(rxyz(2,kkat)-smmd%rxyz(2,iat),acell(2),pery)
                        do j=1,n
                            tmpmat(j,ii) = rr2*lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            case (0)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr3 = closest_image(rxyz(3,kkat)-smmd%rxyz(3,iat),acell(3),perz)
                        do j=1,n
                            tmpmat(j,ii) = rr3*lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            case (1)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr1 = closest_image(rxyz(1,kkat)-smmd%rxyz(1,iat),acell(1),perx)
                        do j=1,n
                            tmpmat(j,ii) = rr1*lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            end select
            call axpy(n**2, 1.d0, tmpmat(1,1), 1, multipole_extracted(1,1), 1)
            call f_free(lmp_extracted)
            call f_free(tmpmat)
        else if (l==2) then
            lmp_extracted = f_malloc((/1.to.n,1.to.n,-1.to.1,0.to.1/),id='lmp_extracted')
            tmpmat = f_malloc((/n,n/),id='tmpmat')
            call extract_matrix(smats, lower_multipole_matrices(0,0)%matrix_compr, &
                 neighborx(1,kat), n, lmp_extracted(1,1,0,0))
            do i=-1,1
                call extract_matrix(smats, lower_multipole_matrices(i,1)%matrix_compr, &
                     neighborx(1,kat), n, lmp_extracted(1,1,i,1))
            end do
            select case (m)
            case (-2)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr1 = closest_image(rxyz(1,kkat)-smmd%rxyz(1,iat),acell(1),perx)
                        rr2 = closest_image(rxyz(2,kkat)-smmd%rxyz(2,iat),acell(2),pery)
                        rr3 = closest_image(rxyz(3,kkat)-smmd%rxyz(3,iat),acell(3),perz)
                        do j=1,n
                            tmpmat(j,ii) = -sqrt(3.d0)*rr1*lmp_extracted(j,ii,-1,1) &
                                           -sqrt(3.d0)*rr2*lmp_extracted(j,ii,1,1) &
                                           +sqrt(3.d0)*rr1*rr2*lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            case (-1)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr1 = closest_image(rxyz(1,kkat)-smmd%rxyz(1,iat),acell(1),perx)
                        rr2 = closest_image(rxyz(2,kkat)-smmd%rxyz(2,iat),acell(2),pery)
                        rr3 = closest_image(rxyz(3,kkat)-smmd%rxyz(3,iat),acell(3),perz)
                        do j=1,n
                            tmpmat(j,ii) = -sqrt(3.d0)*rr2*lmp_extracted(j,ii,0,1) &
                                           -sqrt(3.d0)*rr3*lmp_extracted(j,ii,-1,1) &
                                           +sqrt(3.d0)*rr2*rr3*lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            case (0)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr1 = closest_image(rxyz(1,kkat)-smmd%rxyz(1,iat),acell(1),perx)
                        rr2 = closest_image(rxyz(2,kkat)-smmd%rxyz(2,iat),acell(2),pery)
                        rr3 = closest_image(rxyz(3,kkat)-smmd%rxyz(3,iat),acell(3),perz)
                        do j=1,n
                            tmpmat(j,ii) =  rr1*lmp_extracted(j,ii,1,1) &
                                           +rr2*lmp_extracted(j,ii,-1,1) &
                                           -2.d0*rr3*lmp_extracted(j,ii,0,1) &
                                           +0.5d0*(-rr1**2-rr2**2+&
                                             2.d0*rr3**2)&
                                             *lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            case (1)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr1 = closest_image(rxyz(1,kkat)-smmd%rxyz(1,iat),acell(1),perx)
                        rr2 = closest_image(rxyz(2,kkat)-smmd%rxyz(2,iat),acell(2),pery)
                        rr3 = closest_image(rxyz(3,kkat)-smmd%rxyz(3,iat),acell(3),perz)
                        do j=1,n
                            tmpmat(j,ii) = -sqrt(3.d0)*rr1*lmp_extracted(j,ii,0,1) &
                                           -sqrt(3.d0)*rr3*lmp_extracted(j,ii,1,1) &
                                           +sqrt(3.d0)*rr1*rr3&
                                             *lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            case (2)
                ii = 0
                do i=1,smats%nfvctr
                    if (neighborx(i,kat)) then
                        ii = ii + 1
                        iat = smmd%on_which_atom(i)
                        rr1 = closest_image(rxyz(1,kkat)-smmd%rxyz(1,iat),acell(1),perx)
                        rr2 = closest_image(rxyz(2,kkat)-smmd%rxyz(2,iat),acell(2),pery)
                        rr3 = closest_image(rxyz(3,kkat)-smmd%rxyz(3,iat),acell(3),perz)
                        do j=1,n
                            tmpmat(j,ii) = -sqrt(3.d0)*(rr1)*lmp_extracted(j,ii,1,1) &
                                           +sqrt(3.d0)*(rr2)*lmp_extracted(j,ii,-1,1) &
                                           +sqrt(0.75d0)*(rr1**2-rr2**2)&
                                             *lmp_extracted(j,ii,0,0)
                        end do
                    end if
                end do
            end select
            call axpy(n**2, -1.d0, tmpmat(1,1), 1, multipole_extracted(1,1), 1)
            call f_free(lmp_extracted)
            call f_free(tmpmat)
        end if

      call f_release_routine()

    end subroutine correct_multipole_origin


    ! SM: This is neither tested nor used...
    !!subroutine correct_multipole_origin_new(nat, l, m, n, nlr, natpx, nmaxx, kat, kkat, &
    !!           smats, orbs, rxyz, neighborx, perx, pery, perz, acell, &
    !!           lower_multipole_matrices, locregcenter, multipole_matrix)
    !!  use sparsematrix_base, only: sparse_matrix, matrices
    !!  use module_types, only: orbitals_data
    !!  implicit none

    !!  ! Calling arguments
    !!  integer,intent(in) :: nat, l, m, n, nlr, natpx, nmaxx, kat, kkat
    !!  real(kind=8),dimension(3,nat),intent(in) :: rxyz
    !!  type(sparse_matrix),intent(in) :: smats
    !!  type(orbitals_data),intent(in) :: orbs
    !!  logical,dimension(smats%nfvctr,natpx),intent(in) :: neighborx
    !!  logical,intent(in) :: perx, pery, perz
    !!  real(kind=8),dimension(3),intent(in) :: acell
    !!  type(matrices),dimension(-1:1,0:1),intent(in):: lower_multipole_matrices
    !!  real(kind=8),dimension(3,nlr),intent(in) :: locregcenter
    !!  type(matrices),intent(inout) :: multipole_matrix

    !!  ! Local variables
    !!  integer :: ii, ilr, i, j, iseg
    !!  real(kind=8) :: rr1, rr2, rr3
    !!  real(kind=8),dimension(:,:,:,:),allocatable :: lmp_extracted
    !!  real(kind=8),dimension(:,:),allocatable :: tmpmat

    !!  call f_routine(id='correct_multipole_origin')

    !!    if (l==1) then
    !!        !!lmp_extracted = f_malloc((/1.to.n,1.to.n,0.to.0,0.to.0/),id='lmp_extracted')
    !!        !!tmpmat = f_malloc((/n,n/),id='tmpmat')
    !!        !!call extract_matrix(smats, lower_multipole_matrices(0,0)%matrix_compr, &
    !!        !!     neighborx(1,kat), n, nmaxx, lmp_extracted(1,1,0,0))
    !!        select case (m)
    !!        case (-1)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = rr2*lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        + rr2*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        case (0)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = rr3*lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        + rr3*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        case (1)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = rr1*lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        + rr1*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        end select
    !!        !!call axpy(n**2, 1.d0, tmpmat(1,1), 1, multipole_extracted(1,1), 1)
    !!        !!call f_free(lmp_extracted)
    !!        !!call f_free(tmpmat)
    !!    else if (l==2) then
    !!        !!lmp_extracted = f_malloc((/1.to.n,1.to.n,-1.to.1,0.to.1/),id='lmp_extracted')
    !!        !!tmpmat = f_malloc((/n,n/),id='tmpmat')
    !!        !!call extract_matrix(smats, lower_multipole_matrices(0,0)%matrix_compr, &
    !!        !!     neighborx(1,kat), n, nmaxx, lmp_extracted(1,1,0,0))
    !!        !!do i=-1,1
    !!        !!    call extract_matrix(smats, lower_multipole_matrices(i,1)%matrix_compr, &
    !!        !!         neighborx(1,kat), n, nmaxx, lmp_extracted(1,1,i,1))
    !!        !!end do
    !!        select case (m)
    !!        case (-2)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!            !!        rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!            !!        rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = -sqrt(3.d0)*rr1*lmp_extracted(j,ii,-1,1) &
    !!            !!                           -sqrt(3.d0)*rr2*lmp_extracted(j,ii,1,1) &
    !!            !!                           +sqrt(3.d0)*rr1*rr2*lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!                rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!                rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr1*lower_multipole_matrices(-1,1)%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr2*lower_multipole_matrices(1,1)%matrix_compr(ii) &
    !!                        +sqrt(3.d0)*rr1*rr2*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        case (-1)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!            !!        rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!            !!        rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = -sqrt(3.d0)*rr2*lmp_extracted(j,ii,0,1) &
    !!            !!                           -sqrt(3.d0)*rr3*lmp_extracted(j,ii,-1,1) &
    !!            !!                           +sqrt(3.d0)*rr2*rr3*lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!                rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!                rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr2*lower_multipole_matrices(0,1)%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr3*lower_multipole_matrices(-1,1)%matrix_compr(ii) &
    !!                        +sqrt(3.d0)*rr2*rr3*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        case (0)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!            !!        rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!            !!        rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) =  rr1*lmp_extracted(j,ii,1,1) &
    !!            !!                           +rr2*lmp_extracted(j,ii,-1,1) &
    !!            !!                           -2.d0*rr3*lmp_extracted(j,ii,0,1) &
    !!            !!                           +0.5d0*(-rr1**2-rr2**2+&
    !!            !!                             2.d0*rr3**2)&
    !!            !!                             *lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!                rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!                rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                         +rr1*lower_multipole_matrices(1,1)%matrix_compr(ii) &
    !!                         +rr2*lower_multipole_matrices(-1,1)%matrix_compr(ii) &
    !!                         -2.d0*rr3*lower_multipole_matrices(0,1)%matrix_compr(ii) &
    !!                         +0.5d0*(-rr1**2-rr2**2+&
    !!                           2.d0*rr3**2)&
    !!                           *lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        case (1)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!            !!        rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!            !!        rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = -sqrt(3.d0)*rr1*lmp_extracted(j,ii,0,1) &
    !!            !!                           -sqrt(3.d0)*rr3*lmp_extracted(j,ii,1,1) &
    !!            !!                           +sqrt(3.d0)*rr1*rr3&
    !!            !!                             *lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!                rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!                rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr1*lower_multipole_matrices(0,1)%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr3*lower_multipole_matrices(1,1)%matrix_compr(ii) &
    !!                        +sqrt(3.d0)*rr1*rr3*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        case (2)
    !!            !!ii = 0
    !!            !!do i=1,smats%nfvctr
    !!            !!    if (neighborx(i,kat)) then
    !!            !!        ii = ii + 1
    !!            !!        ilr = orbs%inwhichlocreg(i)
    !!            !!        rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!            !!        rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!            !!        rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!            !!        do j=1,n
    !!            !!            tmpmat(j,ii) = -sqrt(3.d0)*(rr1)*lmp_extracted(j,ii,1,1) &
    !!            !!                           +sqrt(3.d0)*(rr2)*lmp_extracted(j,ii,-1,1) &
    !!            !!                           +sqrt(0.75d0)*(rr1**2-rr2**2)&
    !!            !!                             *lmp_extracted(j,ii,0,0)
    !!            !!        end do
    !!            !!    end if
    !!            !!end do
    !!            ii = 0
    !!            do iseg=1,smats%nseg
    !!                ilr = orbs%inwhichlocreg(smats%keyg(1,2,iseg))
    !!                rr1 = closest_image(rxyz(1,kkat)-locregcenter(1,ilr),acell(1),perx)
    !!                rr2 = closest_image(rxyz(2,kkat)-locregcenter(2,ilr),acell(2),pery)
    !!                rr3 = closest_image(rxyz(3,kkat)-locregcenter(3,ilr),acell(3),perz)
    !!                do i=smats%keyg(1,1,iseg),smats%keyg(2,1,iseg)
    !!                    ii = ii + 1
    !!                    multipole_matrix%matrix_compr(ii) = multipole_matrix%matrix_compr(ii) &
    !!                        -sqrt(3.d0)*rr1*lower_multipole_matrices(1,1)%matrix_compr(ii) &
    !!                        +sqrt(3.d0)*rr2*lower_multipole_matrices(-1,1)%matrix_compr(ii) &
    !!                        +sqrt(0.75d0)*(rr1**2-rr2**2)*lower_multipole_matrices(0,0)%matrix_compr(ii)
    !!                end do
    !!            end do
    !!        end select
    !!        !!call axpy(n**2, -1.d0, tmpmat(1,1), 1, multipole_extracted(1,1), 1)
    !!        !!call f_free(lmp_extracted)
    !!        !!call f_free(tmpmat)
    !!    end if

    !!  call f_release_routine()

    !!end subroutine correct_multipole_origin_new


    subroutine compare_charge_and_potential(boxit,nat,& !iproc, is1, ie1, is2, ie2, is3, ie3, nat, &
               rho_exact, rho_mp, pot_exact, pot_mp, kernel, rxyz, &
               ncheck, check_threshold, charge_error, external_volume, potential_error, potential_total)
    use PStypes, only: coulomb_operator
    use PSbox, only: PS_gather
    use box
    implicit none
    ! Calling arguments
    !integer,intent(in) :: iproc, is1, ie1, is2, ie2, is3, ie3
    type(box_iterator), intent(inout) :: boxit
    integer, intent(in) :: nat, ncheck
    type(coulomb_operator),intent(in) :: kernel
    real(kind=8),dimension(kernel%ndims(1)*kernel%ndims(2)*kernel%grid%n3p),intent(in) :: rho_exact, rho_mp, pot_exact, pot_mp
    real(kind=8),dimension(3,nat),intent(in) :: rxyz
    real(kind=8),dimension(ncheck),intent(in) :: check_threshold
    real(kind=8),dimension(ncheck),intent(out) :: charge_error, external_volume, potential_error, potential_total

    ! Local variables
    integer :: i1, i2, i3, iat, icheck,icnt,igood
    real(kind=8) :: qex, factor,vex
    real(kind=8),parameter :: min_distance = 2.0d0
!!$    logical,dimension(:,:,:),allocatable :: is_close

    call f_routine(id='compare_charge_and_potential')

!!$    ! Determine the grid points which are farther away from the atoms than the minimal distance
!!$    is_close = f_malloc((/0.to.kernel%ndims(1)-31-1,0.to.kernel%ndims(2)-31-1,is3.to.ie3/),id='is_close')
!!$    is_close(:,:,:) = .false.
!!$    do iat=1,nat
!!$        do i3=max(15,is3),min(ie3,kernel%ndims(3)-16-1)
!!$            z = (i3-15)*kernel%hgrids(3)
!!$            d = sqrt( (z-rxyz(3,iat))**2 )
!!$            if (d<=min_distance) then
!!$                ! From the viewpoint of the z coordinate, the grid points might be close to atom iat,
!!$                ! so also check the other directions.
!!$                do i2=0,kernel%ndims(2)-31-1
!!$                    y = i2*kernel%hgrids(2)
!!$                    d = sqrt( (y-rxyz(2,iat))**2 + (z-rxyz(3,iat))**2 )
!!$                    if (d<=min_distance) then
!!$                        ! From the viewpoint of the y and z coordinates, the grid points might be close to atom iat,
!!$                        ! so also check the other directions.
!!$                        do i1=0,kernel%ndims(1)-31-1
!!$                            x = i1*kernel%hgrids(1)
!!$                            d = sqrt( (x-rxyz(1,iat))**2 + (y-rxyz(2,iat))**2 + (z-rxyz(3,iat))**2 )
!!$                            if (d<=min_distance) then
!!$                                is_close(i1,i2,i3) = .true.
!!$                            end if
!!$                        end do
!!$                    end if
!!$                end do
!!$            end if
!!$        end do
!!$    end do

    call f_zero(charge_error)
    call f_zero(external_volume)
    call f_zero(potential_error)
    call f_zero(potential_total)

    !use the box iterator
    factor=boxit%mesh%volume_element
    icnt=0
    igood=0
    box_loop: do while(box_next_point(boxit))
       icnt=icnt+1
       do iat=1,nat
          if (distance(boxit%mesh,boxit%rxyz,rxyz(:,iat)) <= min_distance) cycle box_loop
       end do
       igood=igood+1
       ! Farther away from the atoms than the minimal distance
       qex = rho_exact(boxit%ind)
       vex = pot_exact(boxit%ind)
       do icheck=1,ncheck
          if (abs(qex)<check_threshold(icheck)) then
             ! Charge density smaller than the threshold
             !LG: it seems therefore normal than the density is not good as by hypothesis
             ! we only check that when it is small, thus the mp density will surely be smaller
             charge_error(icheck) = charge_error(icheck) + &
                  abs(rho_mp(boxit%ind))*factor !qex-
             external_volume(icheck) = external_volume(icheck)+ factor
             potential_error(icheck) = potential_error(icheck) + &
                  abs(vex-pot_mp(boxit%ind))*factor 
             potential_total(icheck) = potential_total(icheck) + &
                  abs(vex)*factor
          end if
       end do
    end do box_loop

!!$    do i3=0,kernel%ndims(3)-31-1
!!$        z = i3*kernel%hgrids(3)
!!$        do i2=0,kernel%ndims(2)-31-1
!!$            y = i2*kernel%hgrids(2)
!!$            do i1=0,kernel%ndims(1)-31-1
!!$                x = i1*kernel%hgrids(1)
!!$                dmin = huge(1.d0)
!!$                do iat=1,nat
!!$                    d = sqrt( (x-rxyz(1,iat))**2 + (y-rxyz(2,iat))**2 + (z-rxyz(3,iat))**2 )
!!$                    dmin = min(d,dmin)
!!$                end do
!!$                !write(200+iproc,*) 'x, y, z, q, v', x, y, z, &
!!$                !     rhog_exact(i1+15,i2+15,i3+15), potg_exact(i1+15,i2+15,i3+15)
!!$                !write(300+iproc,*) 'x, y, z, q, v', x, y, z, &
!!$                !     rhog_mp(i1+15,i2+15,i3+15), potg_mp(i1+15,i2+15,i3+15)
!!$                if (dmin>min_distance) then
!!$                    ! Farther away from the atoms than the minimal distance
!!$                    qex = rhog_exact(i1+15,i2+15,i3+15)
!!$                    do icheck=1,ncheck
!!$                        if (abs(qex)<check_threshold(icheck)) then
!!$                            ! Charge density smaller than the threshold
!!$                            charge_error(icheck) = charge_error(icheck) + &
!!$                                abs(qex-rhog_mp(i1+15,i2+15,i3+15))
!!$                            charge_total(icheck) = charge_total(icheck) + &
!!$                                abs(qex)
!!$                            potential_error(icheck) = potential_error(icheck) + &
!!$                                abs(potg_exact(i1+15,i2+15,i3+15)-potg_mp(i1+15,i2+15,i3+15))
!!$                            potential_total(icheck) = potential_total(icheck) + &
!!$                                abs(potg_exact(i1+15,i2+15,i3+15))
!!$                        end if
!!$                    end do
!!$                end if
!!$            end do
!!$        end do
!!$     end do

!!$    call f_free(is_close)

!!$    call dscal(ncheck, factor, charge_error(1), 1)
!!$    call dscal(ncheck, factor, charge_total(1), 1)
!!$    call dscal(ncheck, factor, potential_error(1), 1)
!!$    call dscal(ncheck, factor, potential_total(1), 1)
    if (bigdft_mpi%nproc > 1) then
       call mpiallred(charge_error, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(external_volume, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(potential_error, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(potential_total, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(icnt,1,op=mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(igood,1,op=mpi_sum, comm=bigdft_mpi%mpi_comm)
    end if

    call f_release_routine()

    end subroutine compare_charge_and_potential


    subroutine calculate_weight_center(llr, glr, hgrids, phir, center_locreg, center_orb)
      use module_base
      use module_types, only: locreg_descriptors
      use bounds, only: geocode_buffers
      implicit none

      ! Calling arguments
      type(locreg_descriptors),intent(in) :: llr, glr
      real(kind=8),dimension(3),intent(in) :: hgrids
      real(kind=8),dimension(llr%d%n1i*llr%d%n2i*llr%d%n3i),intent(in) :: phir
      real(kind=8),dimension(3),intent(out) :: center_locreg
      real(kind=8),dimension(3),intent(out) :: center_orb

      ! Local variables
      integer :: nl1, nl2, nl3, i1, i2, i3, ii1, ii2, ii3, ii
      real(kind=8) :: weight, hxh, hyh, hzh, x, y, z, tt
      real(kind=8),dimension(3) :: center

      call f_routine(id='calculate_weight_center')

      hxh = 0.5d0*hgrids(1)
      hyh = 0.5d0*hgrids(2)
      hzh = 0.5d0*hgrids(3)
      call geocode_buffers(llr%geocode, glr%geocode, nl1, nl2, nl3)
      weight = 0.d0
      center(1:3) = 0.0_gp
      !$omp parallel default(none) &
      !$omp shared(llr, nl1, nl2, nl3, hxh, hyh, hzh, phir, center, weight) &
      !$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, tt, ii)
      !$omp do reduction(+: center, weight)
      do i3=1,llr%d%n3i
          ii3 = llr%nsi3 + i3 - nl3 - 1
          z = ii3*hzh
          do i2=1,llr%d%n2i
              ii2 = llr%nsi2 + i2 - nl2 - 1
              y = ii2*hyh
              do i1=1,llr%d%n1i
                  ii1 = llr%nsi1 + i1 - nl1 - 1
                  x = ii1*hxh
                  ii = (i3-1)*llr%d%n2i*llr%d%n1i+(i2-1)*llr%d%n1i+i1
                  tt = phir(ii)**2
                  center(1) = center(1) + x*tt
                  center(2) = center(2) + y*tt
                  center(3) = center(3) + z*tt
                  weight = weight + tt
                  !!ii = ii + 1
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel
      center_locreg(1:3) = center(1:3)/weight
      center_orb(1:3) = center_locreg(1:3)

      call f_release_routine()

    end subroutine calculate_weight_center



    subroutine estimate_system_volume_from_density(boxit,nat,& !iproc, is1, ie1, is2, ie2, is3, ie3, nat, &
               rho_exact, rho_mp, pot_exact, pot_mp, kernel, rxyz, &
               ncheck, check_threshold, charge_error, external_volume, potential_error, potential_total)
    use PStypes, only: coulomb_operator
    use PSbox, only: PS_gather
    use box
    implicit none
    ! Calling arguments
    !integer,intent(in) :: iproc, is1, ie1, is2, ie2, is3, ie3
    type(box_iterator), intent(inout) :: boxit
    integer, intent(in) :: nat, ncheck
    type(coulomb_operator),intent(in) :: kernel
    real(kind=8),dimension(kernel%ndims(1)*kernel%ndims(2)*kernel%grid%n3p),intent(in) :: rho_exact, rho_mp, pot_exact, pot_mp
    real(kind=8),dimension(3,nat),intent(in) :: rxyz
    real(kind=8),dimension(ncheck),intent(in) :: check_threshold
    real(kind=8),dimension(ncheck),intent(out) :: charge_error, external_volume, potential_error, potential_total

    ! Local variables
    integer :: i1, i2, i3, iat, icheck,icnt,igood
    real(kind=8) :: qex, factor,vex
    real(kind=8),parameter :: min_distance = 2.0d0
!!$    logical,dimension(:,:,:),allocatable :: is_close

    call f_routine(id='compare_charge_and_potential')

!!$    ! Determine the grid points which are farther away from the atoms than the minimal distance
!!$    is_close = f_malloc((/0.to.kernel%ndims(1)-31-1,0.to.kernel%ndims(2)-31-1,is3.to.ie3/),id='is_close')
!!$    is_close(:,:,:) = .false.
!!$    do iat=1,nat
!!$        do i3=max(15,is3),min(ie3,kernel%ndims(3)-16-1)
!!$            z = (i3-15)*kernel%hgrids(3)
!!$            d = sqrt( (z-rxyz(3,iat))**2 )
!!$            if (d<=min_distance) then
!!$                ! From the viewpoint of the z coordinate, the grid points might be close to atom iat,
!!$                ! so also check the other directions.
!!$                do i2=0,kernel%ndims(2)-31-1
!!$                    y = i2*kernel%hgrids(2)
!!$                    d = sqrt( (y-rxyz(2,iat))**2 + (z-rxyz(3,iat))**2 )
!!$                    if (d<=min_distance) then
!!$                        ! From the viewpoint of the y and z coordinates, the grid points might be close to atom iat,
!!$                        ! so also check the other directions.
!!$                        do i1=0,kernel%ndims(1)-31-1
!!$                            x = i1*kernel%hgrids(1)
!!$                            d = sqrt( (x-rxyz(1,iat))**2 + (y-rxyz(2,iat))**2 + (z-rxyz(3,iat))**2 )
!!$                            if (d<=min_distance) then
!!$                                is_close(i1,i2,i3) = .true.
!!$                            end if
!!$                        end do
!!$                    end if
!!$                end do
!!$            end if
!!$        end do
!!$    end do

    call f_zero(charge_error)
    call f_zero(external_volume)
    call f_zero(potential_error)
    call f_zero(potential_total)

    !use the box iterator
    factor=boxit%mesh%volume_element
    icnt=0
    igood=0
    box_loop: do while(box_next_point(boxit))
       icnt=icnt+1
       do iat=1,nat
          if (distance(boxit%mesh,boxit%rxyz,rxyz(:,iat)) <= min_distance) cycle box_loop
       end do
       igood=igood+1
       ! Farther away from the atoms than the minimal distance
       qex = rho_exact(boxit%ind)
       vex = pot_exact(boxit%ind)
       do icheck=1,ncheck
          if (abs(qex)<check_threshold(icheck)) then
             ! Charge density smaller than the threshold
             !LG: it seems therefore normal than the density is not good as by hypothesis
             ! we only check that when it is small, thus the mp density will surely be smaller
             charge_error(icheck) = charge_error(icheck) + &
                  abs(rho_mp(boxit%ind))*factor !qex-
             external_volume(icheck) = external_volume(icheck)+ factor
             potential_error(icheck) = potential_error(icheck) + &
                  abs(vex-pot_mp(boxit%ind))*factor 
             potential_total(icheck) = potential_total(icheck) + &
                  abs(vex)*factor
          end if
       end do
    end do box_loop

!!$    do i3=0,kernel%ndims(3)-31-1
!!$        z = i3*kernel%hgrids(3)
!!$        do i2=0,kernel%ndims(2)-31-1
!!$            y = i2*kernel%hgrids(2)
!!$            do i1=0,kernel%ndims(1)-31-1
!!$                x = i1*kernel%hgrids(1)
!!$                dmin = huge(1.d0)
!!$                do iat=1,nat
!!$                    d = sqrt( (x-rxyz(1,iat))**2 + (y-rxyz(2,iat))**2 + (z-rxyz(3,iat))**2 )
!!$                    dmin = min(d,dmin)
!!$                end do
!!$                !write(200+iproc,*) 'x, y, z, q, v', x, y, z, &
!!$                !     rhog_exact(i1+15,i2+15,i3+15), potg_exact(i1+15,i2+15,i3+15)
!!$                !write(300+iproc,*) 'x, y, z, q, v', x, y, z, &
!!$                !     rhog_mp(i1+15,i2+15,i3+15), potg_mp(i1+15,i2+15,i3+15)
!!$                if (dmin>min_distance) then
!!$                    ! Farther away from the atoms than the minimal distance
!!$                    qex = rhog_exact(i1+15,i2+15,i3+15)
!!$                    do icheck=1,ncheck
!!$                        if (abs(qex)<check_threshold(icheck)) then
!!$                            ! Charge density smaller than the threshold
!!$                            charge_error(icheck) = charge_error(icheck) + &
!!$                                abs(qex-rhog_mp(i1+15,i2+15,i3+15))
!!$                            charge_total(icheck) = charge_total(icheck) + &
!!$                                abs(qex)
!!$                            potential_error(icheck) = potential_error(icheck) + &
!!$                                abs(potg_exact(i1+15,i2+15,i3+15)-potg_mp(i1+15,i2+15,i3+15))
!!$                            potential_total(icheck) = potential_total(icheck) + &
!!$                                abs(potg_exact(i1+15,i2+15,i3+15))
!!$                        end if
!!$                    end do
!!$                end if
!!$            end do
!!$        end do
!!$     end do

!!$    call f_free(is_close)

!!$    call dscal(ncheck, factor, charge_error(1), 1)
!!$    call dscal(ncheck, factor, charge_total(1), 1)
!!$    call dscal(ncheck, factor, potential_error(1), 1)
!!$    call dscal(ncheck, factor, potential_total(1), 1)
    if (bigdft_mpi%nproc > 1) then
       call mpiallred(charge_error, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(external_volume, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(potential_error, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(potential_total, mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(icnt,1,op=mpi_sum, comm=bigdft_mpi%mpi_comm)
       call mpiallred(igood,1,op=mpi_sum, comm=bigdft_mpi%mpi_comm)
    end if

    call f_release_routine()

    end subroutine estimate_system_volume_from_density


    subroutine determine_submatrix_sizes(natp, isat, smmd, smat, neighbor, nx, nmax)
      use sparsematrix_base, only: sparse_matrix_metadata, sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none

      ! Calling arguments
      integer,intent(in) :: natp, isat
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smat
      logical,dimension(smat%nfvctr,natp),intent(out) :: neighbor
      integer,dimension(:),intent(out) :: nx
      integer,intent(out) :: nmax

      ! Local variables
      integer :: kat, kkat, n, i, j, iat, iatold, inds

      neighbor(:,:) = .false.
      !ntot = 0
      nmax = 0
      do kat=1,natp
          iatold = 0
          kkat = kat + isat
          n = 0
          do i=1,smat%nfvctr
               iat = smmd%on_which_atom(i)
               ! Only do the following for the first TMB per atom
               if (iat==iatold) cycle
               iatold = iat
               if (iat==kkat) then
                   do j=1,smat%nfvctr
                       inds =  matrixindex_in_compressed(smat, j, i)
                       if (inds/=0) then
                          neighbor(j,kat) = .true.
                          n = n + 1
                       end if
                   end do
               end if
          end do
          nx(kat) = n
          !ntot = ntot + n
          nmax = max(nmax, n)
      end do

    end subroutine determine_submatrix_sizes

end module multipole
