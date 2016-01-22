module multipole
  use module_base
  use multipole_base, only: external_potential_descriptors,lmax
  implicit none

  private

  !> Public routines
  public :: interaction_multipoles_ions
  public :: potential_from_charge_multipoles
  public :: ionic_energy_of_external_charges
  public :: gaussian_density
  public :: multipole_analysis_driver
  public :: projector_for_charge_analysis
  public :: support_function_gross_multipoles
  public :: calculate_dipole_moment
  public :: calculate_rpowerx_matrices

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
    subroutine potential_from_charge_multipoles(iproc, nproc, at, denspot, ep, is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, shift, &
               verbosity, ixc, lzd, pot, rxyz, ixyz0, dipole_total, quadrupole_total, all_norms_ok)
      use module_types, only: DFT_local_fields, local_zone_descriptors
      use Poisson_Solver, except_dp => dp, except_gp => gp
      use module_atoms, only: atoms_data
      use bounds, only: ext_buffers
      use yaml_output
      use io, only: plot_density
      use bounds, only: geocode_buffers
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
      real(kind=8),dimension(3),intent(out),optional :: dipole_total
      real(kind=8),dimension(3,3),intent(out),optional :: quadrupole_total
      logical,intent(out),optional :: all_norms_ok

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l, m, ii, mm, nthread, ithread, ll
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp, ehart_ps, tt, ttt, gg, hhh, tt0, tt1, tt2
      real(dp),dimension(3) :: r
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
      !$ integer  :: omp_get_thread_num,omp_get_max_threads

      call f_routine(id='potential_from_charge_multipoles')

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
          nelpsp = f_malloc(ep%nmpl,id='nelpsp')
          rloc = f_malloc(ep%nmpl,id='rloc')
         perx = (denspot%dpbox%geocode /= 'F')
         pery = (denspot%dpbox%geocode == 'P')
         perz = (denspot%dpbox%geocode /= 'F')
         n3pi = denspot%dpbox%n3pi
         i3s = denspot%dpbox%i3s + denspot%dpbox%i3xcsh
         hxh = denspot%dpbox%hgrids(1)
         hyh = denspot%dpbox%hgrids(2)
         hzh = denspot%dpbox%hgrids(3)
         n1i = denspot%dpbox%ndims(1)
         n2i = denspot%dpbox%ndims(2)
         n3i = denspot%dpbox%ndims(3)
         call ext_buffers(perx,nbl1,nbr1)
         call ext_buffers(pery,nbl2,nbr2)
         call ext_buffers(perz,nbl3,nbr3)
         !write(*,*) 'ep%nmpl, n1i, n2i, n3i', ep%nmpl, n1i, n2i, n3i
    
    
    
    
         ! Generate the density that comes from the pseudopotential atoms
         ndensity = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
         psp_source = f_malloc(ep%nmpl,id='psp_source')
         do impl=1,ep%nmpl
             ! Search the rloc and zion of the corresponding pseudopotential
             call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp(impl), psp_source(impl), rloc(impl))
         end do

         !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
         cutoff=10.0_gp*maxval(rloc(:))
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
             ! Search the rloc and zion of the corresponding pseudopotential
             call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp(impl), psp_source(impl), rloc(impl))
             if(norm_ok(impl)) then
                 ! The following routine needs the shifted positions
                 rx = ep%mpl(impl)%rxyz(1) - shift(1)
                 ry = ep%mpl(impl)%rxyz(2) - shift(2)
                 rz = ep%mpl(impl)%rxyz(3) - shift(3)
                 call gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, &
                      rx, ry, rz, &
                      ep%mpl(impl)%sigma(0), nelpsp(impl), at%multipole_preserving, use_iterator, at%mp_isf, &
                      denspot%dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, ndensity, density_cores, rholeaked)
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
              !rmax(impl) = min(denspot%dpbox%ndims(1)*0.25d0*hx, &
              !                 denspot%dpbox%ndims(2)*0.25d0*hy, &
              !                 denspot%dpbox%ndims(3)*0.25d0*hz)
              rmax(impl) = min((denspot%dpbox%ndims(1)-31)*0.25d0*hx, &
                               (denspot%dpbox%ndims(2)-31)*0.25d0*hy, &
                               (denspot%dpbox%ndims(3)-31)*0.25d0*hz)
          end do
    
    
          norm_check = 0.d0
          monopole = 0.d0
          dipole = 0.d0
          quadrupole = 0.d0
          !$omp parallel &
          !$omp default(none) &
          !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, hhh, ep, shift, nthread, norm_ok) &
          !$omp shared(norm_check, monopole, dipole, quadrupole, density, density_loc, potential_loc) &
          !$omp shared (gaussians1, gaussians2, gaussians3, nelpsp, rmax, rmin) &
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
                  do i3=is3,ie3
                      ii3 = i3 - nl3 -1
                      !z = real(ii3,kind=8)*hz + shift(3)
                      r(3) = huge(r(3))
                      do j3=j3s,j3e
                          dr = real(ii3+j3*lzd%glr%d%n3i,kind=8)*hz + shift(3) - ep%mpl(impl)%rxyz(3)
                          if (abs(dr)<abs(r(3))) r(3) = dr
                      end do
                      do i2=is2,ie2
                          ii2 = i2 - nl2 - 1
                          !y = real(ii2,kind=8)*hy + shift(2)
                          r(2) = huge(r(2))
                          do j2=j2s,j2e
                              dr = real(ii2+j2*lzd%glr%d%n2i,kind=8)*hy + shift(2) - ep%mpl(impl)%rxyz(2)
                              if (abs(dr)<abs(r(2))) r(2) = dr
                          end do
                          do i1=is1,ie1
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
                                              qq = -(ep%mpl(impl)%qlm(l)%q(mm) - real(nelpsp(impl),kind=8))
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
                          end do
                      end do
                  end do
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
          if (verbosity> 0 .and. iproc==0) call write_psp_source(ep, psp_source)
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
          call f_free(psp_source)
    
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
                                      qq = -(ep%mpl(impl)%qlm(l)%q(mm)-real(nelpsp(impl),kind=8))
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

             ! Add the core contribution
             call axpy((ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), 1.0_gp, density_cores(is1,is2,is3), 1, density(is1,is2,is3), 1)


             call H_potential('D',denspot%pkernel,density,denspot%V_ext,ehart_ps,0.0_dp,.false.,&
                  quiet='yes')!,rho_ion=denspot%rho_ion)
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
              call plot_density(iproc,nproc,'mppot.cube',at,rxyz_noshift,denspot%pkernel,nspin=1,rho=density, &
                   ixyz0=ixyz0_)
              call f_free(rxyz_noshift)
          end if
    
          call f_free(density)
          call f_free(density_cores)
          !call f_free(nzatom)
          call f_free(nelpsp)
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

      dpm = q(1)*r(1) + q(2)*r(2) + q(3)*r(3)
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

      qq(1,1) = q(1)
      qq(2,1) = q(2)
      qq(3,1) = q(3)
      qq(1,2) = qq(2,1)
      qq(2,2) = q(4)
      qq(3,2) = q(5)
      qq(1,3) = qq(3,1)
      qq(2,3) = qq(3,2)
      qq(3,3) = 1.0_dp-qq(1,1)-qq(2,2)

      qpm = qq(1,1)*r(1)*r(1) + &
           qq(2,1)*r(2)*r(1) + &
           qq(3,1)*r(3)*r(1) + &
           qq(1,2)*r(1)*r(2) + &
           qq(2,2)*r(2)*r(2) + &
           qq(3,2)*r(3)*r(2) + &
           qq(1,3)*r(1)*r(3) + &
           qq(2,3)*r(2)*r(3) + &
           qq(3,3)*r(3)*r(3)
      qpm = -0.5_dp*qpm/rnrm5

    end function calc_quadropole



    !> Calculate S^1/2 * K * S^1/2, which is the kernel corresponding to a
    !! orthonormal set of support functions.
    subroutine kernel_for_orthonormal_basis(iproc, nproc, norbp, meth_overlap, smats, smatl, &
               ovrlp, kernel, weight_matrix_compr)
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, SPARSE_TASKGROUP, &
                                   matrices_null, assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   deallocate_matrices
      use matrix_operations, only: overlapPowerGeneral
      use sparsematrix, only: matrix_matrix_mult_wrapper, gather_matrix_from_taskgroups
      use yaml_output
      implicit none

      ! Calling arguments
      integer :: iproc, nproc, norbp, meth_overlap
      type(sparse_matrix),intent(in) :: smats, smatl
      type(matrices),intent(in) :: kernel
      type(matrices),intent(in) :: ovrlp
      real(kind=8),dimension(smatl%nvctrp_tg*smatl%nspin),intent(out) :: weight_matrix_compr

      ! Local variables
      type(matrices),dimension(1) :: inv_ovrlp
      real(kind=8),dimension(:),allocatable :: weight_matrix_compr_tg, proj_ovrlp_half_compr
      real(kind=8) :: max_error, mean_error

      call f_routine(id='kernel_for_orthonormal_basis')

      if (iproc==0) then
          call yaml_comment('Calculating kernel for orthonormal support functions',hfill='~')
      end if

      inv_ovrlp(1) = matrices_null()
      inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='inv_ovrlp(1)%matrix_compr')

      call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), -1, &
           imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
           ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
           max_error=max_error, mean_error=mean_error)
      !call f_free_ptr(ovrlp%matrix)

      proj_ovrlp_half_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='proj_mat_compr')
      !if (norbp>0) then
         call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
              kernel%matrix_compr, inv_ovrlp(1)%matrix_compr, proj_ovrlp_half_compr)
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

      if (iproc==0) then
          call yaml_comment('Kernel calculated',hfill='~')
      end if

      call f_release_routine()

    end subroutine kernel_for_orthonormal_basis




    subroutine write_multipoles_new(ep, units, delta_rxyz, on_which_atom, scaled, monopoles_analytic)
      use yaml_output
      use numerics, only: Bohr_Ang
      use f_precisions, only: db => f_double
      implicit none
      
      ! Calling arguments
      !integer,dimension(nat),intent(in) :: iatype
      type(external_potential_descriptors),intent(in) :: ep
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
      real(kind=8),dimension(-lmax:lmax,0:lmax) :: multipoles
      logical :: present_delta_rxyz, present_on_which_atom, present_scaled, present_monopoles_analytic

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
              if (present_delta_rxyz) then
                  call yaml_map('Delta r',convert_units*delta_rxyz(1:3,impl),fmt='(es13.6)')
              end if
              if (any(ep%mpl(impl)%sigma(0:lmax)/=0.d0)) then
                  call yaml_map('sigma',ep%mpl(impl)%sigma(0:lmax),fmt='(f5.3)')
              endif
              call f_zero(multipoles)
              if (present_monopoles_analytic) then
                  call yaml_map('q0 analytic',monopoles_analytic(impl),fmt='(1es13.6)')
              end if
              do l=0,lmax
                  call yaml_map('q'//adjustl(trim(yaml_toa(l))),ep%mpl(impl)%qlm(l)%q(:),fmt='(1es13.6)')
                  multipoles(-l:l,l) = ep%mpl(impl)%qlm(l)%q(1:2*l+1)
                  call yaml_newline()
              end do
              if (present_scaled) then
                  call yaml_map('scaling factor',scaled(impl),fmt='(es9.2)')
              end if
              function_type = guess_type(multipoles)
              call yaml_map('type',trim(function_type))
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()




          contains
            ! Try to guess the type (s, p, etc.) of a support function
            function guess_type(mp) result(gt)
              implicit none
              ! Calling arguments
              real(kind=8),dimension(-lmax:lmax,0:lmax),intent(in) :: mp
              character(len=9) :: gt
              ! Local variables
              integer :: il, im, ilmax, immax
              real(kind=8) :: maxvalue1, maxvalue2

              
              ! A type is recognized if an element is at least four times as large as all other elements
              maxvalue1 = 0.d0 !the largest element
              maxvalue2 = 0.d0 !the second largest element
              do il=0,lmax
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










    subroutine calculte_multipole_matrix(iproc, nproc, l, m, nphi, phi1, phi2, nphir, hgrids, &
               orbs, collcom, lzd, smat, locregcenter, ingegration_volume, multipole_matrix)
      use module_base
      use module_types, only: orbitals_data, comms_linear, local_zone_descriptors
      use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
      use sparsematrix_base, only: sparse_matrix, matrices
      use communications_base, only: TRANSPOSE_FULL
      use transposed_operations, only: calculate_overlap_transposed
      use communications, only: transpose_localized
      use bounds, only: geocode_buffers
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, l, m, nphi, nphir
      real(kind=8),dimension(nphi),intent(in) :: phi1, phi2
      real(kind=8),dimension(3) :: hgrids
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(3,lzd%nlr),intent(in) :: locregcenter
      type(matrices),intent(inout) :: multipole_matrix
      character(len=*),intent(in) :: ingegration_volume

      ! Local variables
      integer :: ist, istr, iorb, i1, i2, i3, ii1, ii2, ii3, iiorb, ind, ilr, i, nl1, nl2, nl3
      integer :: i1mod, i2mod, i3mod, is1, ie1, is2, ie2, is3, ie3, ii, nd, nu, j1, j2, j3
      real(kind=8),dimension(:),allocatable :: phi2r, sphi2r, sphi2, phi1t_c, phi1t_f, sphi2t_c, sphi2t_f
      real(kind=8) :: norm, rmax, factor_normalization, tt, x, y, z
      type(workarr_sumrho) :: w
      real(kind=8) :: ddot, dr
      character(len=*),parameter :: sphere = 'sphere', box = 'box'
      logical :: integrate_in_sphere, perx, pery, perz
      integer :: j1s, j1e, j2s, j2e, j3s, j3e


      call f_routine(id='calculte_multipole_matrix')

      ! Check the arguments
      if (trim(ingegration_volume)==sphere) then
          integrate_in_sphere = .true.
      else if (trim(ingegration_volume)==box) then
          integrate_in_sphere = .false.
      else
          call f_err_throw('wrong argument for ingegration_volume ('//trim(ingegration_volume)//')',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end if

      ! Transform the support functions to real space
      phi2r = f_malloc0(max(nphir,1),id='phi2r')
      sphi2r = f_malloc0(max(nphir,1),id='sphi2r')
      ist=1
      istr=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.true.,w)
          call daub_to_isf(lzd%llr(ilr), w, phi2(ist), phi2r(istr))
          !write(*,*) 'iorb, n, tt', iorb, &
          !     lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, &
          !     lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f, &
          !     ddot(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, phi2r(istr), 1, phi2r(istr), 1), &
          !     ddot(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f, phi2(ist), 1, phi2(ist), 1)
          call deallocate_work_arrays_sumrho(w)
          call deallocate_work_arrays_sumrho(w)
          ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
      end do
      !!do i=1,ist-1
      !!    write(556,*) i, phi_ortho(i)
      !!end do

      !do i=1,nphir
      !    write(700,*) i, phi2r(i)
      !end do
      !if(istr/=collcom_sr%ndimpsi_c+1) then
      !    write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
      !    stop
      !end if

      !write(*,*) 'after daub_to_isf'

      ! Conditions for periodicity
      perx=(smat%geocode /= 'F')
      pery=(smat%geocode == 'P')
      perz=(smat%geocode /= 'F')
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



      ! Apply the spherical harmonic
      ist = 0
      do iorb=1,orbs%norbp
          iiorb = orbs%isorb + iorb
          ilr = orbs%inwhichlocreg(iiorb)
          !rmax = min(lzd%llr(ilr)%d%n1i*0.25d0*hgrids(1),lzd%llr(ilr)%d%n2i*0.25d0*hgrids(2),lzd%llr(ilr)%d%n3i*0.25d0*hgrids(3))
          rmax = min(lzd%llr(ilr)%d%n1*0.5d0*hgrids(1),lzd%llr(ilr)%d%n2*0.5d0*hgrids(2),lzd%llr(ilr)%d%n3*0.5d0*hgrids(3))
          !write(*,*) 'iorb, ilr, rmax', iorb, ilr, rmax
          !write(*,*) 'zmin, zmax, locregcenter(3)',  (lzd%llr(ilr)%nsi3+1-14-1)*0.5d0*hgrids(3), &
          !            (lzd%llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-14-1)*0.5d0*hgrids(3), locregcenter(3,ilr)
          !write(*,*) 'ymin, ymax, locregcenter(2)',  (lzd%llr(ilr)%nsi2+1-14-1)*0.5d0*hgrids(2), &
          !            (lzd%llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i-14-1)*0.5d0*hgrids(2), locregcenter(2,ilr)
          !write(*,*) 'xmin, xmax, locregcenter(1)',  (lzd%llr(ilr)%nsi1+1-14-1)*0.5d0*hgrids(1), &
          !            (lzd%llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i-14-1)*0.5d0*hgrids(1), locregcenter(1,ilr)

          call geocode_buffers(lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)

          ! Calculate the boundaries:
          ! - free BC: entire box
          ! - periodic BC: half of the box size, with periodic wrap around
          if (.not.perx) then
              is1 = 1
              ie1 = lzd%llr(ilr)%d%n1i
          else
              nd = (lzd%llr(ilr)%d%n1i + 1)/2 - 1
              nu = lzd%llr(ilr)%d%n1i - nd -1
              !write(*,*) 'nu, nd, lzd%llr(ilr)%d%n1i',nu, nd, lzd%llr(ilr)%d%n1i
              if (nu+nd+1/=lzd%llr(ilr)%d%n1i) call f_err_throw('wrong values of nu and nd')
              ii1 = nint(locregcenter(1,ilr)/(0.5d0*hgrids(1)))
              is1 = ii1 - nd + nl1 + 1
              ie1 = ii1 + nu + nl1 + 1
              if (ie1-is1+1/=lzd%llr(ilr)%d%n1i) call f_err_throw('wrong values of is1 and ie1')
          end if
          if (.not.pery) then
              is2 = 1
              ie2 = lzd%llr(ilr)%d%n2i
          else
              nd = (lzd%llr(ilr)%d%n2i + 1)/2 - 1
              nu = lzd%llr(ilr)%d%n2i - nd -1
              if (nu+nd+1/=lzd%llr(ilr)%d%n2i) call f_err_throw('wrong values of nu and nd')
              ii2 = nint(locregcenter(2,ilr)/(0.5d0*hgrids(2)))
              is2 = ii2 - nd + nl2 + 1
              ie2 = ii2 + nu + nl2 + 1
              if (ie2-is2+1/=lzd%llr(ilr)%d%n2i) call f_err_throw('wrong values of is2 and ie2')
          end if
          if (.not.perz) then
              is3 = 1
              ie3 = lzd%llr(ilr)%d%n3i
          else
              nd = (lzd%llr(ilr)%d%n3i + 1)/2 - 1
              nu = lzd%llr(ilr)%d%n3i - nd -1
              if (nu+nd+1/=lzd%llr(ilr)%d%n3i) call f_err_throw('wrong values of nu and nd')
              ii3 = nint(locregcenter(3,ilr)/(0.5d0*hgrids(3)))
              is3 = ii3 - nd + nl3 + 1
              ie3 = ii3 + nu + nl3 + 1
              if (ie3-is3+1/=lzd%llr(ilr)%d%n3i) call f_err_throw('wrong values of is3 and ie3')
          end if

          !!write(*,*) 'perx', perx
          !!write(*,*) 'pery', pery
          !!write(*,*) 'perz', perz

          !!write(*,*) 'iorb, is1, ie1, ii1, n1, is2, ie2, ii2, n2, is3, ie3, ii3, n3', &
          !!    iorb, is1, ie1, ii1, lzd%llr(ilr)%d%n1i, is2, ie2, ii2, lzd%llr(ilr)%d%n2i, is3, ie3, ii3, lzd%llr(ilr)%d%n3i


          norm = 0.d0
          factor_normalization = sqrt(0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3))
          !do i3=is3,ie3
          do i3=1,lzd%llr(ilr)%d%n3i
              !j3 = i3 - is3 + 1
              !i3mod = modulo(i3-1,lzd%llr(ilr)%d%n3i) + 1
              ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
              !z = ii3*0.5d0*hgrids(3) - locregcenter(3,ilr)
              ! Search the closest locregcenter (might be in a periodically replicated cell)
              z = huge(z)
              do j3=j3s,j3e
                  dr = (ii3+j3*lzd%glr%d%n3i)*0.5d0*hgrids(3) - locregcenter(3,ilr)
                  if (abs(dr)<abs(z)) z = dr
              end do
              !write(*,*) 'is3, ie3, n3, i3, i3mod, ii3, z', is3, ie3, lzd%llr(ilr)%d%n3i, i3, i3mod, ii3, z
              !write(*,*) 'is3, ie3, n3, i3, ii3, z', is3, ie3, lzd%llr(ilr)%d%n3i, i3, ii3, z
              !do i2=is2,ie2
              do i2=1,lzd%llr(ilr)%d%n2i
                  !j2 = i2 - is2 + 1
                  !i2mod = modulo(i2-1,lzd%llr(ilr)%d%n2i) + 1
                  ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
                  !y = ii2*0.5d0*hgrids(2) - locregcenter(2,ilr)
                  ! Search the closest locregcenter (might be in a periodically replicated cell)
                  y = huge(y)
                  do j2=j2s,j2e
                      dr = (ii2+j2*lzd%glr%d%n2i)*0.5d0*hgrids(2) - locregcenter(2,ilr)
                      if (abs(dr)<abs(y)) y = dr
                  end do
                  !do i1=is1,ie1
                  do i1=1,lzd%llr(ilr)%d%n1i
                      !j1 = i1 - is1 + 1
                      !i1mod = modulo(i1-1,lzd%llr(ilr)%d%n1i) + 1
                      ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                      !x = ii1*0.5d0*hgrids(1) - locregcenter(1,ilr)
                      x = huge(x)
                      do j1=j1s,j1e
                          dr = (ii1+j1*lzd%glr%d%n1i)*0.5d0*hgrids(1) - locregcenter(1,ilr)
                          if (abs(dr)<abs(x)) x = dr
                      end do
                      !ind = (i3mod-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2mod-1)*lzd%llr(ilr)%d%n1i + i1mod
                      !ind = (j3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (j2-1)*lzd%llr(ilr)%d%n1i + j1
                      ind = (i3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2-1)*lzd%llr(ilr)%d%n1i + i1
                      !ind = (i3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2-1)*lzd%llr(ilr)%d%n1i + i1
                      !!if (i1/=i1mod .or. i2/=i2mod .or. i3/=i3mod) then
                      !!    write(*,*) 'iproc, is1, ie1, is2, ie2, is3, ie3, i1, i2, i3, i1mod, i2mod, i3mod', &
                      !!                iproc, is1, ie1, is2, ie2, is3, ie3, i1, i2, i3, i1mod, i2mod, i3mod
                      !!end if
                      !!if (iorb==orbs%norbp) then
                      !!    write(bigdft_mpi%iproc+20,*) 'is1, ie1, is2, ie2, is3, ie3, i1, i2, i3, i1mod, i2mod, i3mod', &
                      !!                                  is1, ie1, is2, ie2, is3, ie3, i1, i2, i3, i1mod, i2mod, i3mod
                      !!    write(bigdft_mpi%iproc+30,*) 'ind',ind
                      !!end if
                      if (integrate_in_sphere) then
                          if (x**2+y**2+z**2>rmax**2) cycle
                      end if
                      tt = solid_harmonic(0, 0.d0, l, m, x, y, z)
                      tt = tt*sqrt(4.d0*pi/real(2*l+1,kind=8))
                      sphi2r(ist+ind) = tt*phi2r(ist+ind)
                      !!if (iorb==orbs%norbp) then
                      !!    write(bigdft_mpi%iproc+40,*) 'is1, ie1, is2, ie2, is3, ie3, i1, i2, i3, i1mod, i2mod, i3mod', &
                      !!                                  is1, ie1, is2, ie2, is3, ie3, i1, i2, i3, i1mod, i2mod, i3mod
                      !!    write(bigdft_mpi%iproc+50,*) 'ind, phi2r(ist+ind), sphi2r(ist+ind)',ind, phi2r(ist+ind), sphi2r(ist+ind)
                      !!end if
                      !write(*,*) 'iorb, i1, i1, i2, tt, phi2r', iorb, i1, i2, i3, tt, phi2r(ist+ind)
                      ! For the calculation of the norm, do the integration always only in the sphere
                      if (x**2+y**2+z**2>rmax**2) cycle
                      norm = norm + (tt*factor_normalization)**2*&
                          real((2*l+3)*(2*l+1),kind=8)/(4.d0*pi*rmax**(2*l+3)) !normalization of a solid harmonic within a sphere of radius rmax... hopefully correct
                      !write(*,*) 'iorb, i1, i2, i3, tt, phi', iorb, i1, i2, i3, tt, phir(ist+ind)
                  end do
              end do
          end do
          ist = ist + ind
      end do



      ! Transform back to wavelets
      sphi2 = f_malloc0(nphi,id='sphi2')
      ist=1
      istr=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.true.,w)
          call isf_to_daub(lzd%llr(ilr), w, sphi2r(istr), sphi2(ist))
          !call isf_to_daub(lzd%llr(ilr), w, phi2r(istr), sphi2(ist))
          call deallocate_work_arrays_sumrho(w)
          !write(*,*) 'iorb, n, firsts, tt', iorb, &
          !     lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, &
          !     phi2r(istr), sphi2r(istr), &
          !     ddot(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, phi2r(istr), 1, sphi2r(istr), 1), &
          !     ddot(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, phi2r(istr), 1, phi2r(istr), 1), &
          !     ddot(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f, phi1(ist), 1, sphi2(ist), 1)
          !do i=1,lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
          !    write(*,*) i, phir(istr+i-1), sphir(istr+i-1)
          !end do
          ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
      end do
      !if(istr/=collcom_sr%ndimpsi_c+1) then
      !    write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
      !    stop
      !end if

      !write(*,*) 'after isf_to_daub'

      call f_free(phi2r)
      call f_free(sphi2r)



      ! Calculate the scalar products, i.e. the matrix <phi_ab|S_lm|phi_ab>
      phi1t_c = f_malloc(collcom%ndimind_c,id='phi1t_c')
      phi1t_f = f_malloc(7*collcom%ndimind_f,id='phi1t_f')
      sphi2t_c = f_malloc(collcom%ndimind_c,id='sphit2_c')
      sphi2t_f = f_malloc(7*collcom%ndimind_f,id='sphit2_f')
      call transpose_localized(iproc, nproc, nphi, orbs, collcom, &
           TRANSPOSE_FULL, phi1, phi1t_c, phi1t_f, lzd)
      call transpose_localized(iproc, nproc, nphi, orbs, collcom, &
           TRANSPOSE_FULL, sphi2, sphi2t_c, sphi2t_f, lzd)
      call calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
           phi1t_c, sphi2t_c, phi1t_f, sphi2t_f, smat, multipole_matrix)
      !call calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
      !     phi1t_c, phi1t_c, phi1t_f, phi1t_f, smat, multipole_matrix)

      !!do i=1,size(multipole_matrix%matrix_compr)
      !!    write(*,*) 'i, val', i, multipole_matrix%matrix_compr(i)
      !!end do

      !write(*,*) 'after overlap'

      call f_free(sphi2)
      call f_free(phi1t_c)
      call f_free(phi1t_f)
      call f_free(sphi2t_c)
      call f_free(sphi2t_f)

      call f_release_routine()

      !!call mpi_finalize(ii)
      !!stop

    end subroutine calculte_multipole_matrix



    subroutine multipole_analysis_driver(iproc, nproc, ll, nphi, lphi, nphir, at, hgrids, &
               orbs, smats, smatm, smatl, collcom, collcom_sr, lzd, denspot, orthpar, ovrlp, ham, kernel, rxyz, &
               method, do_ortho, shift, nsigma, ixc, ep)
      use module_base
      use module_types, only: orbitals_data, comms_linear, local_zone_descriptors, orthon_data, DFT_local_fields, comms_linear
      use sparsematrix_base, only: sparse_matrix, matrices, sparsematrix_malloc0, assignment(=), &
                                   sparsematrix_malloc, matrices_null, sparsematrix_malloc_ptr, deallocate_matrices, &
                                   SPARSE_TASKGROUP
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix_local
      use communications, only: transpose_localized
      use orthonormalization, only: orthonormalizelocalized
      use module_atoms, only: atoms_data
      use yaml_output
      use multipole_base, only: external_potential_descriptors_null, multipole_set_null, multipole_null, &
                                deallocate_external_potential_descriptors
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ll, nphi, nphir, nsigma, ixc
      real(kind=8),dimension(3) :: hgrids
      type(atoms_data),intent(in) :: at
      type(orbitals_data),intent(in) :: orbs
      type(sparse_matrix),intent(in) :: smats
      type(sparse_matrix),intent(in) :: smatm
      type(sparse_matrix),intent(in) :: smatl
      type(comms_linear),intent(in) :: collcom, collcom_sr
      type(local_zone_descriptors),intent(in) :: lzd
      type(DFT_local_fields),intent(inout) :: denspot
      type(orthon_data),intent(in) :: orthpar
      type(matrices),intent(in) :: ovrlp
      type(matrices),intent(in) :: ham
      type(matrices),intent(in) :: kernel
      real(kind=8),dimension(nphi),intent(in) :: lphi
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      character(len=*),intent(in) :: method
      character(len=*),intent(in) :: do_ortho
      real(kind=8),dimension(3),intent(in) :: shift
      type(external_potential_descriptors),intent(out) :: ep

      ! Local variables
      integer :: methTransformOverlap, iat, ind, ispin, ishift, iorb, iiorb, l, m, itype, natpx, isatx, nmaxx, kat, n, i, kkat
      integer :: ilr, impl, mm, lcheck, nelpsp, psp_source
      logical :: can_use_transposed, all_norms_ok
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      real(kind=8),dimension(:),allocatable :: phi_ortho, Qmat, kernel_ortho, multipole_matrix_large, Qmat_tmp
      real(kind=8),dimension(:,:),allocatable :: Qmat_tilde, kp, locregcenter, overlap_small
      real(kind=8),dimension(:,:,:),pointer :: atomic_multipoles
      real(kind=8),dimension(:),pointer :: atomic_monopoles_analytic
      real(kind=8),dimension(:,:,:),allocatable :: test_pot
      real(kind=8),dimension(:,:),pointer :: projx
      real(kind=8) :: q, tt, rloc
      type(matrices) :: multipole_matrix
      logical,dimension(:,:),pointer :: neighborx
      integer,dimension(:),pointer :: nx
      character(len=20),dimension(:),allocatable :: names
      real(kind=8),dimension(3) :: dipole_check
      real(kind=8),dimension(3,3) :: quadrupole_check
      type(external_potential_descriptors) :: ep_check
      character(len=*),parameter :: no='no', yes='yes'
      !character(len=*),parameter :: do_ortho = no!yes


      call f_routine(id='multipole_analysis_driver')


      if (iproc==0) then
          call yaml_comment('Atomic multipole analysis, new approach',hfill='=')
          call yaml_map('Method',trim(method))
          call yaml_map('Orthogonalized support functions',trim(do_ortho))
      end if

      call unitary_test_multipoles(iproc, nproc, nphi, nphir, orbs, lzd, smats, collcom, hgrids)

      ! Check the method
      if (trim(method)/='projector' .and. trim(method)/='loewdin') then
          call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (trim(do_ortho)/='no' .and. trim(do_ortho)/='yes') then
          call f_err_throw('wrong do_ortho',err_name='BIGDFT_RUNTIME_ERROR')
      end if


      if (trim(method)=='projector') then
          if (smatl%nspin/=1) then
              call f_err_throw('projector not tested for spin polarized calculations, better to stop here')
          end if
          ! Calculate the projector using the penalty term
          call projector_for_charge_analysis(at, smats, smatm, smatl, &
               ovrlp, ham, kernel, rxyz, calculate_centers=.false., write_output=.false., ortho=do_ortho, &
               natpx=natpx, isatx=isatx, nmaxx=nmaxx, nx=nx, projx=projx, neighborx=neighborx)
      end if

      multipole_matrix_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='multipole_matrix_large')
      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='kernel_ortho')

      if (do_ortho==yes) then
          ! Calculate the kernel for orthormal support functions
          methTransformOverlap = 20
          call kernel_for_orthonormal_basis(iproc, nproc, orbs%norbp, methTransformOverlap, smats, smatl, &
               ovrlp, kernel, kernel_ortho)
       !else if (do_ortho==no) then
       !    ! Calculate K*S, use multipole_matrix_large as workarray
       !    call transform_sparse_matrix(smats, smatl, 'small_to_large', &
       !         smat_in=ovrlp%matrix_compr, lmat_out=multipole_matrix_large)
       !    call matrix_matrix_mult_wrapper(iproc, nproc, smatl, kernel%matrix_compr, multipole_matrix_large, kernel_ortho)
       end if

      if (do_ortho==yes) then
          ! Orthogonalize the support functions
          can_use_transposed = .false.
          methTransformOverlap = 1020
          phit_c = f_malloc_ptr(collcom%ndimind_c,id='phit_c')
          phit_f = f_malloc_ptr(7*collcom%ndimind_f,id='phit_f')
          phi_ortho = f_malloc(nphi,id='phi_ortho')
          call vcopy(nphi, lphi(1), 1, phi_ortho(1), 1)
          if (iproc==0) then
              call yaml_comment('Orthonormalizing the support functions',hfill='~')
          end if
          call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
               1.d-8, nphi, orbs, lzd, &
               smats, smatl, collcom, orthpar, &
               phi_ortho, phit_c, phit_f, &
               can_use_transposed)
          call f_free_ptr(phit_c)
          call f_free_ptr(phit_f)
      end if


      Qmat = sparsematrix_malloc(smatl,iaction=SPARSE_TASKGROUP,id='Qmat')
      atomic_multipoles = f_malloc0_ptr((/-ll.to.ll,0.to.ll,1.to.at%astruct%nat/),id='atomic_multipoles')
      atomic_monopoles_analytic = f_malloc0_ptr(1.to.at%astruct%nat,id='atomic_monopoles_analytic')


      multipole_matrix = matrices_null()
      multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='multipole_matrix%matrix_compr')


      locregcenter = f_malloc((/3,lzd%nlr/),id='locregcenter')
      do ilr=1,lzd%nlr
          locregcenter(1:3,ilr) = lzd%llr(ilr)%locregcenter(1:3)
      end do

      do l=0,ll
          do m=-l,l

              call f_zero(multipole_matrix%matrix_compr)

              ! Calculate the multipole matrix
              if (do_ortho==yes) then
                  call calculte_multipole_matrix(iproc, nproc, l, m, nphi, phi_ortho, phi_ortho, nphir, hgrids, &
                       orbs, collcom, lzd, smats, locregcenter, 'box', multipole_matrix)
              else if (do_ortho==no) then
                  call calculte_multipole_matrix(iproc, nproc, l, m, nphi, lphi, lphi, nphir, hgrids, &
                       orbs, collcom, lzd, smats, locregcenter, 'box', multipole_matrix)
              end if

              call transform_sparse_matrix_local(smats, smatl, 'small_to_large', &
                   smatrix_compr_in=multipole_matrix%matrix_compr, lmatrix_compr_out=multipole_matrix_large)

              ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
              call dscal(smatl%nvctrp_tg*smatl%nspin, -1.d0, multipole_matrix_large(1), 1)

              ! Multiply the orthogonalized kernel with the multipole matrix
              call f_zero(Qmat)
              if (do_ortho==yes) then
                  call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                       kernel_ortho, multipole_matrix_large, Qmat)
               else if (do_ortho==no) then
                   !if (trim(method)=='projector') then
                   !    ! Calculate K*S, use Qmat as workarray
                   !    call transform_sparse_matrix(smats, smatl, 'small_to_large', &
                   !         smat_in=ovrlp%matrix_compr, lmat_out=Qmat)
                   !    call matrix_matrix_mult_wrapper(iproc, nproc, smatl, kernel%matrix_compr, Qmat, kernel_ortho)
                   !    call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                   !         kernel_ortho, multipole_matrix_large, Qmat)
                   !else if (trim(method)=='loewdin') then
                       call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                            kernel%matrix_compr, multipole_matrix_large, Qmat)
                   !end if
               end if

              if (trim(method)=='projector') then
                  do kat=1,natpx
                      kkat = kat + isatx
                      n = nx(kat)
                      qmat_tilde = f_malloc((/n,n/),id='qmat_tilde')
                      kp = f_malloc((/n,n/),id='kp')
                      call extract_matrix(smatl, qmat, neighborx(1:,kat), n, nmaxx, qmat_tilde)
                      call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, projx(1,kat), n, 0.d0, kp(1,1), n)
                      if (do_ortho==no) then
                          overlap_small = f_malloc((/n,n/),id='overlap_small')
                          call extract_matrix(smats, ovrlp%matrix_compr, neighborx(1:,kat), n, nmaxx, overlap_small)
                          call f_memcpy(src=kp,dest=qmat_tilde)
                          call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, overlap_small(1,1), n, 0.d0, kp(1,1), n)
                          call f_free(overlap_small)
                      end if
                      tt = 0.d0
                      do i=1,n
                          tt = tt + kp(i,i)
                      end do
                      atomic_multipoles(m,l,kkat) = tt
                      call f_free(qmat_tilde)
                      call f_free(kp)
                  end do
              else if (trim(method)=='loewdin') then
                  do ispin=1,smatl%nspin
                      ishift = (ispin-1)*smatl%nvctrp_tg
                      ! Need to do this in parallel (norbp), since the matrices might not be fully filled (matrix taskgroup etc.)
                      ! This should be carefull checked again.
                      do iorb=1,orbs%norbp
                          iiorb = modulo(orbs%isorb+iorb-1,smatl%nfvctr)+1
                          iat=smatl%on_which_atom(iiorb)
                          ind = matrixindex_in_compressed(smatl, iiorb, iiorb) - smatl%isvctrp_tg
                          ind = ind + ishift
                          atomic_multipoles(m,l,iat) = atomic_multipoles(m,l,iat) + qmat(ind)
                      end do
                  end do
              end if

              ! For the monopole, do the same with the exact overlap matrix
              if (l==0) then
                  call transform_sparse_matrix_local(smats, smatl, 'small_to_large', &
                       smatrix_compr_in=ovrlp%matrix_compr, lmatrix_compr_out=multipole_matrix_large)

                  ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
                  call dscal(smatl%nvctrp_tg*smatl%nspin, -1.d0, multipole_matrix_large(1), 1)

                  ! Multiply the kernel with the multipole matrix. Use the original kernel not the orthogonalized one
                  call f_zero(Qmat)
                  call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                       kernel%matrix_compr, multipole_matrix_large, Qmat)

                  !!if (trim(method)=='projector') then
                  !!    do kat=1,natpx
                  !!        kkat = kat + isatx
                  !!        n = nx(kat)
                  !!        Qmat_tilde = f_malloc((/n,n/),id='Qmat_tilde')
                  !!        kp = f_malloc((/n,n/),id='kp')
                  !!        call extract_matrix(smatl, Qmat, neighborx(1:,kat), n, nmaxx, Qmat_tilde)
                  !!        call gemm('n', 'n', n, n, n, 1.d0, Qmat_tilde(1,1), n, projx(1,kat), n, 0.d0, kp(1,1), n)
                  !!        tt = 0.d0
                  !!        do i=1,n
                  !!            tt = tt + kp(i,i)
                  !!        end do
                  !!        atomic_monopoles_analytic(kkat) = tt
                  !!        call f_free(Qmat_tilde)
                  !!        call f_free(kp)
                  !!    end do
                  !!else if (trim(method)=='loewdin') then
                  !!    do ispin=1,smatl%nspin
                  !!        ishift = (ispin-1)*smatl%nvctr
                  !!        do iorb=1,orbs%norb
                  !!            iiorb = modulo(iorb-1,smatl%nfvctr)+1
                  !!            iat=smatl%on_which_atom(iiorb)
                  !!            ind = matrixindex_in_compressed(smatl, iorb, iorb)
                  !!            ind = ind + ishift
                  !!            atomic_monopoles_analytic(iat) = atomic_monopoles_analytic(iat) + Qmat(ind)
                  !!        end do
                  !!    end do
                  !!end if
                  if (trim(method)=='projector') then
                      do kat=1,natpx
                          kkat = kat + isatx
                          n = nx(kat)
                          qmat_tilde = f_malloc((/n,n/),id='qmat_tilde')
                          kp = f_malloc((/n,n/),id='kp')
                          call extract_matrix(smatl, qmat, neighborx(1:,kat), n, nmaxx, qmat_tilde)
                          call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, projx(1,kat), n, 0.d0, kp(1,1), n)
                          if (do_ortho==no) then
                              overlap_small = f_malloc((/n,n/),id='overlap_small')
                              call extract_matrix(smats, ovrlp%matrix_compr, neighborx(1:,kat), n, nmaxx, overlap_small)
                              call f_memcpy(src=kp,dest=qmat_tilde)
                              call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, overlap_small(1,1), n, 0.d0, kp(1,1), n)
                              call f_free(overlap_small)
                          end if
                          tt = 0.d0
                          do i=1,n
                              tt = tt + kp(i,i)
                          end do
                          atomic_monopoles_analytic(kkat) = tt
                          call f_free(qmat_tilde)
                          call f_free(kp)
                      end do
                  else if (trim(method)=='loewdin') then
                      do ispin=1,smatl%nspin
                          ishift = (ispin-1)*smatl%nvctrp_tg
                          ! Need to do this in parallel (norbp), since the matrices might not be fully filled (matrix taskgroup etc.)
                          ! This should be carefull checked again.
                          do iorb=1,orbs%norbp
                              iiorb = modulo(orbs%isorb+iorb-1,smatl%nfvctr)+1
                              iat=smatl%on_which_atom(iiorb)
                              ! SM: verufy this taskgroup shift
                              ind = matrixindex_in_compressed(smatl, iiorb, iiorb) - smatl%isvctrp_tg
                              ind = ind + ishift
                              atomic_monopoles_analytic(iat) = atomic_monopoles_analytic(iat) + qmat(ind)
                          end do
                      end do
                  end if
              end if

          end do
      end do

      call f_free(locregcenter)

      !if (trim(method)=='projector') then
          call mpiallred(atomic_multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(atomic_monopoles_analytic, mpi_sum, comm=bigdft_mpi%mpi_comm)
      !end if


      ! The monopole term should be the net charge, i.e. add the positive atomic charges
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          q = real(at%nelpsp(itype),kind=8)
          atomic_multipoles(0,0,iat) = atomic_multipoles(0,0,iat) + q
          atomic_monopoles_analytic(iat) = atomic_monopoles_analytic(iat) + q
      end do

      names = f_malloc_str(len(names),at%astruct%nat,id='names')
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          names(iat) = at%astruct%atomnames(itype)
      end do


      ep = external_potential_descriptors_null()
      ep%nmpl = at%astruct%nat
      allocate(ep%mpl(ep%nmpl))
      do impl=1,ep%nmpl
          ep%mpl(impl) = multipole_set_null()
          allocate(ep%mpl(impl)%qlm(0:lmax))
          ep%mpl(impl)%rxyz = at%astruct%rxyz(1:3,impl)
          ep%mpl(impl)%sym = trim(names(impl))
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

      !!! Calculate the optimal sigmas
      !!call get_optimal_sigmas(iproc, nproc, nsigma, collcom_sr, smatl, kernel, at, lzd, ep, shift, rxyz, ixc, denspot)
      do impl=1,ep%nmpl
          call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp, psp_source, rloc)
          ep%mpl(impl)%sigma(0:lmax) = rloc
      end do

      if (iproc==0) then
          call yaml_comment('Final result of the multipole analysis',hfill='~')
          call write_multipoles_new(ep, at%astruct%units, monopoles_analytic=atomic_monopoles_analytic)
      end if

      call f_free_ptr(atomic_monopoles_analytic)


      ! Calculate the total dipole moment resulting from the previously calculated multipoles.
      ! This is done by calling the following routine (which actually calculates the potential, but also
      ! has the option to calculte the dipole on the fly).
      test_pot = f_malloc((/size(denspot%V_ext,1),size(denspot%V_ext,2),size(denspot%V_ext,3)/),id='test_pot')
      if (iproc==0) call yaml_sequence_open('Checking the total multipoles based on the atomic multipoles')
      do lcheck=0,lmax
          ep_check = external_potential_descriptors_null()
          ep_check%nmpl = ep%nmpl
          allocate(ep_check%mpl(ep_check%nmpl))
          do impl=1,ep_check%nmpl
              ep_check%mpl(impl) = multipole_set_null()
              allocate(ep_check%mpl(impl)%qlm(0:lmax))
              ep_check%mpl(impl)%rxyz = ep%mpl(impl)%rxyz
              ep_check%mpl(impl)%sym = ep%mpl(impl)%sym
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
               denspot%dpbox%ndims(1), 1, denspot%dpbox%ndims(2), &
               denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
               denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
               denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
               denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3), &
               shift, verbosity=0, ixc=ixc, lzd=lzd, pot=test_pot, &
               rxyz=rxyz, dipole_total=dipole_check, quadrupole_total=quadrupole_check, &
               all_norms_ok=all_norms_ok)
          if (.not. all_norms_ok) then
              call f_err_throw('When checking the previously calculated multipoles, all norms should be ok')
          end if
          dipole_check=dipole_check/0.393430307_gp  ! au2debye              
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
              call yaml_mapping_close()
              call deallocate_external_potential_descriptors(ep_check)
          end if
      end do
      if (iproc==0) call yaml_sequence_close()
      call f_free(test_pot)


      call f_free_str(len(names),names)
      call deallocate_matrices(multipole_matrix)
      call f_free(kernel_ortho)
      call f_free(Qmat)
      if (do_ortho==yes) then
          call f_free(phi_ortho)
      end if
      if (trim(method)=='projector') then
          call f_free_ptr(projx)
          call f_free_ptr(nx)
          call f_free_ptr(neighborx)
      end if
      call f_free_ptr(atomic_multipoles)
      call f_free_ptr(atomic_monopoles_analytic)
      call f_free(multipole_matrix_large)

      if (iproc==0) then
          call yaml_comment('Atomic multipole analysis done',hfill='=')
      end if

      call f_release_routine()

  end subroutine multipole_analysis_driver



    subroutine projector_for_charge_analysis(at, smats, smatm, smatl, &
               ovrlp_, ham_, kernel_, rxyz, calculate_centers, write_output, ortho, &
               lzd, nphirdim, psi, orbs, &
               multipoles, &
               natpx, isatx, nmaxx, nx, projx, neighborx, &
               rpower_matrix)
      use module_base
      use module_types, only: local_zone_descriptors, orbitals_data
      use module_atoms, only: atoms_data
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   sparsematrix_malloc, sparsematrix_malloc0, &
                                   sparsematrix_malloc_ptr, sparsematrix_malloc0_ptr, &
                                   SPARSE_TASKGROUP, assignment(=), &
                                   matrices_null, deallocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix
      use matrix_operations, only: overlapPowerGeneral, overlap_plus_minus_one_half_exact
      use yaml_output
      use multipole_base, only: lmax
      use io, only: write_partial_charges
      implicit none

      ! Calling arguments
      type(atoms_data),intent(in) :: at
      type(sparse_matrix),intent(in) :: smats, smatl
      type(sparse_matrix),intent(in) :: smatm
      type(matrices),intent(in) :: ovrlp_
      type(matrices),intent(in) :: ham_, kernel_
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      logical,intent(in) :: calculate_centers, write_output
      character(len=*),intent(in) :: ortho
      type(local_zone_descriptors),intent(in),optional :: lzd
      integer,intent(in),optional :: nphirdim
      real(kind=8),dimension(:),intent(in),optional :: psi
      type(orbitals_data),intent(in),optional :: orbs
      real(kind=8),dimension(-lmax:lmax,0:lmax,1:smats%nfvctr),intent(in),optional :: multipoles
      integer,intent(out),optional :: natpx, isatx, nmaxx
      integer,dimension(:),pointer,intent(out),optional :: nx
      real(kind=8),dimension(:,:),pointer,intent(out),optional :: projx
      logical,dimension(:,:),pointer,intent(out),optional :: neighborx
      type(matrices),dimension(24),intent(in),optional :: rpower_matrix

      ! Local variables
      integer :: kat, iat, jat, i, j, ii, jj, icheck, n, indm, inds, ntot, ist, ind, iq, itype, ieval, ij, nmax, indl, lwork
      integer :: k, l, iatold, isat, natp, kkat, istot, ntotp, i1, i2, i3, is1, ie1, is2, ie2, is3, ie3, j1, j2, j3, ikT, info
      integer :: ialpha
      real(kind=8) :: r2, cutoff2, rr2, tt, ef, q, occ, max_error, mean_error, rr2i, rr2j, ttxi, ttyi, ttzi, ttxj, ttyj, ttzj
      real(kind=8) :: tti, ttj, charge_net, charge_total
      real(kind=8) :: xi, xj, yi, yj, zi, zj, ttx, tty, ttz, xx, yy, zz, x, y, z
      real(kind=8),dimension(:),allocatable :: work, occ_all
      real(kind=8),dimension(:,:),allocatable :: com
      real(kind=8),dimension(:,:),allocatable :: ham, ovrlp, proj, ovrlp_tmp, ovrlp_minusonehalf, kp, ktilde
      real(kind=8),dimension(:,:,:),allocatable :: coeff_all, ovrlp_onehalf_all, penaltymat
      integer,dimension(:,:,:,:),allocatable :: ilup
      real(kind=8),dimension(:),allocatable :: eval, eval_all, ovrlp_large, tmpmat1, tmpmat2, kerneltilde, charge_per_atom
      real(kind=8),dimension(:,:,:),allocatable :: tmpmat2d
      integer,dimension(:),allocatable :: id_all, n_all, itmparr
      real(kind=8),dimension(3) :: rr
      logical,dimension(:,:),allocatable :: neighbor
      type(matrices),dimension(1) :: ovrlp_onehalf_
      logical :: perx, pery, perz, final, bound_low_ok, bound_up_ok
      !real(kind=8),parameter :: kT = 5.d-2
      real(kind=8) :: kT, ttt
      !real(kind=8),parameter :: alpha = 5.d-1
      real(kind=8) :: alpha, alpha_up, alpha_low, convergence_criterion
      real(kind=8),dimension(:,:,:),allocatable :: multipoles_fake, penalty_matrices
      real(kind=8),dimension(:),allocatable :: alpha_calc
      character(len=*),parameter :: mode='old'

      call f_routine(id='projector_for_charge_analysis')

      if (bigdft_mpi%iproc==0) then
          call yaml_comment('Calculate projector for multipole analysis',hfill='~')
      end if

      if (present(natpx) .or. present(isatx) .or. present(nmaxx) .or. &
          present(projx) .or. present(neighborx) .or. present(nx)) then
          if (.not. (present(natpx) .and. present(isatx) .and. present(nmaxx) .and. &
              present(projx) .and. present(neighborx) .and.  present(nx))) then
              call f_err_throw('not all optional arguments present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      end if

      kT = 1.d-2

      ! Convergence criterion: one million-th of the total charge
      tt = 0.d0
      do iat=1,at%astruct%nat
          tt = tt + real(at%nelpsp(at%astruct%iatype(iat)),kind=8)
      end do
      convergence_criterion = 1.d-6*abs(tt)

      ! Check the arguments
      if (calculate_centers) then
      !!    ! The centers of the support functions are already given
      !!    if (.not.present(com_)) then
      !!        call f_err_throw('com_ not present',err_name='BIGDFT_RUNTIME_ERROR')
      !!    end if
      !!    if (size(com_,1)/=3) then
      !!        call f_err_throw('wrong first dimension of com_',err_name='BIGDFT_RUNTIME_ERROR')
      !!    end if
      !!    if (size(com_,2)/=smats%nfvctr) then
      !!        call f_err_throw('wrong second dimension of com_',err_name='BIGDFT_RUNTIME_ERROR')
      !!    end if
      !!    com => com_
      !!else
          ! Must calculate the centers of the support functions
          if (.not.present(lzd)) then
              call f_err_throw('lzd not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(nphirdim)) then
              call f_err_throw('nphirdim not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(psi)) then
              call f_err_throw('psi not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(orbs)) then
              call f_err_throw('orbs not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      end if


      kerneltilde = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='kerneltilde')
      if (ortho=='yes') then
          ! Calculate S^1/2
          ovrlp_onehalf_(1) = matrices_null()
          ovrlp_onehalf_(1)%matrix_compr = &
              sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_onehalf_(1)%matrix_compr')
          call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, 1020, 1, (/2/), -1, &
                imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
                ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_(1), &
                check_accur=.true., max_error=max_error, mean_error=mean_error)

          ! Calculate S^1/2 * K * S^1/2 = Ktilde
          tmpmat1 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat1')
          !tmpmat2 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat2')
          call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
               kernel_%matrix_compr, ovrlp_onehalf_(1)%matrix_compr, tmpmat1)
          call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
               ovrlp_onehalf_(1)%matrix_compr, tmpmat1, kerneltilde)
      else if (ortho=='no') then
          ovrlp_large = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_large')
          call transform_sparse_matrix(smats, smatl, 'small_to_large', &
               smat_in=ovrlp_%matrix_compr, lmat_out=ovrlp_large)
          call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
               kernel_%matrix_compr, ovrlp_large, kerneltilde)
          call f_free(ovrlp_large)
      end if


      ! Parallelization over the number of atoms
      ii = at%astruct%nat/bigdft_mpi%nproc
      natp = ii
      jj = at%astruct%nat - bigdft_mpi%nproc*natp
      if (bigdft_mpi%iproc<jj) then
          natp = natp + 1
      end if
      isat = (bigdft_mpi%iproc)*ii + min(bigdft_mpi%iproc,jj)

      if (present(natpx)) natpx = natp
      if (present(isatx)) isatx = isat


      ! Determine the sum of the size of all submatrices (i.e. the total number of eigenvalues we will have)
      ! and the maximal value for one atom.
      neighbor = f_malloc((/smats%nfvctr,natp/),id='neighbor')
      neighbor(:,:) = .false.
      ntot = 0
      nmax = 0
      do kat=1,natp
          iatold = 0
          kkat = kat + isat
          n = 0
          do i=1,smats%nfvctr
               iat = smats%on_which_atom(i)
               ! Only do the following for the first TMB per atom
               if (iat==iatold) cycle
               iatold = iat
               if (iat==kkat) then
                   do j=1,smats%nfvctr
                       inds =  matrixindex_in_compressed(smats, j, i)
                       if (inds/=0) then
                          neighbor(j,kat) = .true.
                          n = n + 1
                       end if
                   end do
               end if
          end do
          ntot = ntot + n
          nmax = max(nmax, n)
      end do

      if (present(neighborx)) then
          neighborx = f_malloc_ptr((/smats%nfvctr,natpx/),id='neighborx')
          !call f_memcpy(src=neighbor,dest=neighborx)
          neighborx = neighbor
      end if

      if (present(nx)) then
          nx = f_malloc_ptr(natpx,id='nx')
      end if


      itmparr = f_malloc0(0.to.bigdft_mpi%nproc-1,id='itmparr')
      itmparr(bigdft_mpi%iproc) = ntot
      ntotp = ntot
      if (bigdft_mpi%nproc>1) then
          call mpiallred(itmparr, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(ntot, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(nmax, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
      end if
      istot = 0
      do i=0,bigdft_mpi%iproc-1
          istot = istot + itmparr(i)
      end do
      call f_free(itmparr)
      
      eval_all = f_malloc0(ntot,id='eval_all')
      occ_all = f_malloc0(ntot,id='occ_all')
      id_all = f_malloc0(ntot,id='id_all')
      coeff_all = f_malloc((/nmax,nmax,natp/),id='coeff_all')
      ovrlp_onehalf_all = f_malloc((/nmax,nmax,natp/),id='ovrlp_onehalf_all')
      ovrlp_minusonehalf = f_malloc((/nmax,nmax/),id='ovrlp_minusonehalf')
      ilup = f_malloc((/2,nmax,nmax,natp/),id='ilup')
      n_all = f_malloc(natp,id='n_all')

      penalty_matrices = f_malloc((/nmax,nmax,natp/),id='penalty_matrices')
      alpha_calc = f_malloc(natp,id='alpha_calc')


      ! Centers of the support functions
      com = f_malloc0((/3,smats%nfvctr/),id='com')
      if (calculate_centers) then
          if (orbs%norb>0) then
              !call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, tmb%collcom_sr%ndimpsi_c, &
              !     orbs%norb, orbs%norbp, orbs%isorb, orbs%in_which_locreg, lzd, com(1:,orbs%isorb+1:))
              call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, nphirdim, &
                   orbs%norb, orbs%norbp, orbs%isorb, orbs%inwhichlocreg, lzd, com(1:,orbs%isorb+1:))
              if (bigdft_mpi%nproc>1) then
                  call mpiallred(com, mpi_sum, comm=bigdft_mpi%mpi_comm)
              end if
          end if
      else
          do i=1,smats%nfvctr
              iat = smats%on_which_atom(i)
              com(1:3,i) = rxyz(1:3,iat)
          end do
      end if


      charge_per_atom = f_malloc0(at%astruct%nat,id='charge_per_atom')


      ! Calculate how many states should be included
      q = 0.d0
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          q = q + ceiling(0.5d0*real(at%nelpsp(itype),kind=8))
      end do
      iq = nint(q)
      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_open('Calculating projector for charge analysis')
          call yaml_map('convergence criterion',convergence_criterion)
          call yaml_map('maximal size of a submatrix',nmax)
          call yaml_sequence_open('Searching alpha for charge neutrality')
      end if

      ! Initial guess for the bisection bounds
      alpha_low = 1.d-3
      alpha_up = 1.d1
      bound_low_ok = .false.
      bound_up_ok = .false.

      eF = -1.d0

      alpha_loop: do ialpha=1,10000

          if (bigdft_mpi%iproc==0) then
              call yaml_sequence(advance='no')
          end if

          if (.not.bound_low_ok) then
              alpha = alpha_low
          else if (.not.bound_up_ok) then
              alpha = alpha_up
          else
              alpha = 0.5d0*(alpha_low+alpha_up)
          end if

          if (ialpha==0) alpha = 2.d-1

          charge_net = 0.d0
          call f_zero(eval_all)
          call f_zero(id_all)

          ist = 0
          do kat=1,natp
              kkat = kat + isat
    
              ! Determine the size of the submatrix
              n = 0
              do j=1,smats%nfvctr
                  if (neighbor(j,kat)) then
                      n = n + 1
                  end if
              end do
              n_all(kat) = n
              if (present(nx)) nx(kat) = n
    
    
              ! Extract the submatrices
              ham = f_malloc0((/n,n/),id='ham')
              ovrlp = f_malloc0((/n,n/),id='ovrlp')
              proj = f_malloc0((/n,n/),id='proj')
              penaltymat = f_malloc0((/n,n,24/),id='penaltymat')
              eval = f_malloc0((/n/),id='eval')
              call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1:,kat), n, nmax, ovrlp, ilup)
              call extract_matrix(smatm, ham_%matrix_compr, neighbor(1:,kat), n, nmax, ham)



              if (ortho=='yes') then
                  ! Calculate ovrlp^1/2 and ovrlp^-1/2. The last argument is wrong, clean this.
                  ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                  call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                  call overlap_plus_minus_one_half_exact(1, n, -1, .true., ovrlp_tmp, smats)
                  do i=1,n
                      call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
                  end do
                  call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                  call overlap_plus_minus_one_half_exact(1, n, -1, .false., ovrlp_tmp, smats)
                  do i=1,n
                      call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_minusonehalf(1,i), 1)
                  end do
                  call f_free(ovrlp_tmp)
    
                  ! Calculate S^-1/2 * H * S^-1/2
                  tmpmat2d = f_malloc((/n,n,1/),id='tmppmat2d')
                  call gemm('n', 'n', n, n, n, 1.d0, ham(1,1), n, ovrlp_minusonehalf(1,1), nmax, 0.d0, tmpmat2d(1,1,1), n)
                  call gemm('n', 'n', n, n, n, 1.d0, ovrlp_minusonehalf(1,1), nmax, tmpmat2d(1,1,1), n, 0.d0, ham(1,1), n)
                  call f_free(tmpmat2d)
              end if



              ! Add the penalty term
              if (mode=='verynew') then
                  ! @ NEW #####################################################################
                  if (.not.present(rpower_matrix)) then
                      call f_err_throw('rpower_matrix not present')
                  end if
                  if (.not.present(orbs)) then
                      call f_err_throw('orbs not present')
                  end if
                  write(*,*) 'call extract_matrix with penaltymat'
                  write(*,*) 'orbs%onwhichatom',orbs%onwhichatom
                  write(*,*) 'rxyz(:,kkat)',rxyz(:,kkat)
                  !write(*,*) 'HACK: SET ALPHA TO 0.5d0'
                  !alpha = 0.02d0
                  do i=1,24
                      call extract_matrix(smats, rpower_matrix(i)%matrix_compr, neighbor(1:,kat), n, nmax, penaltymat(:,:,i))
                  end do
                  !tt = sqrt(rxyz(1,kkat)**2+rxyz(2,kkat)**2+rxyz(3,kkat)**2)
                  tt = rxyz(1,kkat)**2 + rxyz(2,kkat)**2 + rxyz(3,kkat)**2
                  do i=1,n
                      do j=1,n
                          !if (i==j .and. orbs%onwhichatom(i)/=kkat) then
                          if (i==j) then
                              !!ttt = penaltymat(j,i,4) &
                              !!    - 2.d0*(rxyz(1,kkat)*penaltymat(j,i,1) &
                              !!          + rxyz(2,kkat)*penaltymat(j,i,2) &
                              !!          + rxyz(3,kkat)*penaltymat(j,i,3)) &
                              !!    + tt*ovrlp(j,i)
                              !!ttt = penaltymat(j,i,2) &
                              !!      + penaltymat(j,i,6) &
                              !!      + penaltymat(j,i,10) &
                              !!    - 2.d0*(rxyz(1,kkat)*penaltymat(j,i,1) &
                              !!          + rxyz(2,kkat)*penaltymat(j,i,5) &
                              !!          + rxyz(3,kkat)*penaltymat(j,i,9)) &
                              !!    + tt*ovrlp(j,i)
                              ttt = penaltymat(j,i,4) - 4.d0*rxyz(1,kkat)*penaltymat(j,i,3) &
                                    + 6.d0*rxyz(1,kkat)**2*penaltymat(j,i,2) - 4.d0*rxyz(1,kkat)**3*penaltymat(j,i,1) &
                                    + rxyz(1,kkat)**4*ovrlp(j,i) &
                                    + penaltymat(j,i,8) - 4.d0*rxyz(2,kkat)*penaltymat(j,i,7) &
                                    + 6.d0*rxyz(2,kkat)**2*penaltymat(j,i,6) - 4.d0*rxyz(2,kkat)**3*penaltymat(j,i,5) &
                                    + rxyz(2,kkat)**4*ovrlp(j,i) &
                                    + penaltymat(j,i,12) - 4.d0*rxyz(3,kkat)*penaltymat(j,i,11) &
                                    + 6.d0*rxyz(3,kkat)**2*penaltymat(j,i,10) - 4.d0*rxyz(3,kkat)**3*penaltymat(j,i,9) &
                                    + rxyz(3,kkat)**4*ovrlp(j,i) &
                                    + 2.d0*(penaltymat(j,i,16) &
                                            - 2.d0*rxyz(2,kkat)*penaltymat(j,i,14) &
                                            + rxyz(2,kkat)**2*penaltymat(j,i,2) &
                                            - 2.d0*rxyz(1,kkat)*penaltymat(j,i,15) &
                                            + 4.d0*rxyz(1,kkat)*rxyz(2,kkat)*penaltymat(j,i,13) &
                                            - 2.d0*rxyz(1,kkat)*rxyz(2,kkat)**2*penaltymat(j,i,1) &
                                            + rxyz(1,kkat)**2*penaltymat(j,i,6) &
                                            - 2.d0*rxyz(1,kkat)**2*rxyz(2,kkat)*penaltymat(j,i,5) &
                                            + rxyz(1,kkat)**2*rxyz(2,kkat)**2*ovrlp(j,i) &
                                            + penaltymat(j,i,20) &
                                            - 2.d0*rxyz(3,kkat)*penaltymat(j,i,18) &
                                            + rxyz(3,kkat)**2*penaltymat(j,i,2) &
                                            - 2.d0*rxyz(1,kkat)*penaltymat(j,i,19) &
                                            + 4.d0*rxyz(1,kkat)*rxyz(3,kkat)*penaltymat(j,i,17) &
                                            - 2.d0*rxyz(1,kkat)*rxyz(3,kkat)**2*penaltymat(j,i,1) &
                                            + rxyz(1,kkat)**2*penaltymat(j,i,10) &
                                            - 2.d0*rxyz(1,kkat)**2*rxyz(3,kkat)*penaltymat(j,i,9) &
                                            + rxyz(1,kkat)**2*rxyz(3,kkat)**2*ovrlp(j,i) &
                                            + penaltymat(j,i,24) &
                                            - 2.d0*rxyz(3,kkat)*penaltymat(j,i,22) &
                                            + rxyz(3,kkat)**2*penaltymat(j,i,6) &
                                            - 2.d0*rxyz(2,kkat)*penaltymat(j,i,23) &
                                            + 4.d0*rxyz(2,kkat)*rxyz(3,kkat)*penaltymat(j,i,21) &
                                            - 2.d0*rxyz(2,kkat)*rxyz(3,kkat)**2*penaltymat(j,i,5) &
                                            + rxyz(2,kkat)**2*penaltymat(j,i,10) &
                                            - 2.d0*rxyz(2,kkat)**2*rxyz(3,kkat)*penaltymat(j,i,9) &
                                            + rxyz(2,kkat)**2*rxyz(3,kkat)**2*ovrlp(j,i) )
                              if (kkat==1) then
                                  ttt = alpha*ttt
                              else
                                  !if (i==1 .or. i==6) ttt=3.d0*ttt
                                  ttt = 1.4641d0*alpha*ttt
                                  !ttt = 5.0d0*alpha*ttt
                              end if
                              !ttt = alpha*penaltymat(j,i,2) - alpha*2.d0*tt*penaltymat(j,i,1) + alpha*tt**2*ovrlp(j,i)
                              ham(j,i) = ham(j,i) + ttt
                              write(*,*) 'kkat, j, i, owa(j), owa(i), alpha, tt, pm1, pm2, pm3, pm4, ovrlp, ttt', &
                                          kkat, j, i, orbs%onwhichatom(j), orbs%onwhichatom(i), &
                                          alpha, tt, penaltymat(j,i,1), penaltymat(j,i,2), &
                                          penaltymat(j,i,3), penaltymat(j,i,4), &
                                          ovrlp(j,i), ttt
                          end if
                      end do
                  end do
                  ! @ END NEW #################################################################
              elseif (mode=='old') then
                  if (ortho=='yes') then
                      ! directly add the penalty terms to ham
                      call add_penalty_term(at%astruct%geocode, smats%nfvctr, neighbor(1:,kat), rxyz(1:,kkat), &
                           at%astruct%cell_dim, com, alpha, n, ovrlp, ham)
                   else if (ortho=='no') then
                          ! Calculate ovrlp^1/2. The last argument is wrong, clean this.
                          ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                          call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                          call overlap_plus_minus_one_half_exact(1, n, -1, .true., ovrlp_tmp, smats)
                          do i=1,n
                              call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
                          end do
                          call f_free(ovrlp_tmp)
                          ! Calculate the penaly term separately and then calculate S^1/2*penalty*S^1/2
                          tmpmat2d = f_malloc0((/n,n,2/),id='tmppmat2d')
                          call add_penalty_term(at%astruct%geocode, smats%nfvctr, neighbor(1:,kat), rxyz(1:,kkat), &
                               at%astruct%cell_dim, com, alpha, n, ovrlp, tmpmat2d(1,1,1))

                          ! Calculate S^1/2 * penalty * S^1/2
                          call gemm('n', 'n', n, n, n, 1.d0, tmpmat2d(1,1,1), n, &
                               ovrlp_onehalf_all(1,1,kat), nmax, 0.d0, tmpmat2d(1,1,2), n)
                          call gemm('n', 'n', n, n, n, 1.d0, ovrlp_onehalf_all(1,1,kat), nmax, &
                               tmpmat2d(1,1,2), n, 0.d0, tmpmat2d(1,1,1), n)
                          call axpy(n**2, 1.d0, tmpmat2d(1,1,1), 1, ham(1,1), 1)
                          call f_free(tmpmat2d)
                      end if
              else if (mode=='new') then
                  multipoles_fake = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.smats%nfvctr/),id='multipoles_fake')
                  multipoles_fake = 0.d0
                  multipoles_fake(0,0,:) = 1.d0
                  if (ialpha==1) then
                      if (present(multipoles)) then
                          write(*,*) 'call with multipoles'
                          call add_penalty_term_new(at%astruct%geocode, at%astruct%nat, smats%nfvctr, &
                               neighbor(1:,kat), rxyz(1:,kkat), smats%on_which_atom, &
                               multipoles, at%astruct%cell_dim, com, alpha, n, ham, &
                               nmax, penalty_matrices(1:n,1:n,kat))
                      else
                          write(*,*) 'call with multipoles_fake'
                          call add_penalty_term_new(at%astruct%geocode, at%astruct%nat, smats%nfvctr, &
                               neighbor(1:,kat), rxyz(1:,kkat), smats%on_which_atom, &
                               multipoles_fake, at%astruct%cell_dim, com, alpha, n, ham, &
                               nmax, penalty_matrices(1:n,1:n,kat))
                      end if
                      alpha_calc(kat) = alpha
                  else
                      tt = alpha/alpha_calc(kat)
                      !write(*,*) 'tt',tt
                      !do i=1,n
                      !    do j=1,n
                      !        write(*,*) 'i, j, penmat', i, j, penalty_matrices(j,i,kat)
                      !    end do
                      !end do
                      ham(1:n,1:n) = ham(1:n,1:n) + tt*penalty_matrices(1:n,1:n,kat)
                  end if
                  call f_free(multipoles_fake)
              end if

    
    
              !!call diagonalizeHamiltonian2(bigdft_mpi%iproc, n, ham, ovrlp, eval)
              lwork = 100*n
              work = f_malloc(lwork,id='work')
              if (ortho=='yes') then
                  call syev('v', 'l', n, ham(1,1), n, eval(1), work(1), lwork, info)
              else if (ortho=='no') then
                  ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                  call f_memcpy(src=ovrlp,dest=ovrlp_tmp)
                  call sygv(1, 'v', 'l', n, ham(1,1), n, ovrlp_tmp(1,1), n, eval(1), work(1), lwork, info)
                  call f_free(ovrlp_tmp)
              end if
              call f_free(work)
              do i=1,n
                  ii = ist + i
                  eval_all(istot+ii) = eval(i)
                  id_all(istot+ii) = kkat
                  call vcopy(n, ham(1,i), 1, coeff_all(1,i,kat), 1)
              end do
    
              ist = ist + n
    
              call f_free(ham)
              call f_free(ovrlp)
              call f_free(penaltymat)
              call f_free(proj)
              call f_free(eval)
    
          end do

    
          if (ist/=ntotp) call f_err_throw('ist/=ntotp',err_name='BIGDFT_RUNTIME_ERROR')
    
          if (bigdft_mpi%nproc>1) then
              call mpiallred(eval_all, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(id_all, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
    
    
          ! Order the eigenvalues and IDs
          call order_eigenvalues(ntot, eval_all, id_all)
        
        


          if (ialpha==1) then
              if (present(nmaxx)) nmaxx = maxval(n_all)
              if (present(projx)) projx = f_malloc_ptr((/nmaxx**2,natpx/),id='projx')
          end if
    
    
              !ikT = ikT + 1

              !kT = 1.d-1*abs(eval_all(1))
              if (mode=='verynew') then
                  kT = max(1.d-1*abs(eval_all(5)),1.d-1)
                  write(*,*) 'adjust kT to',kT
              end if
    
              call f_zero(charge_per_atom)
    
              ! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
              if (ialpha>=0) then
                  ef = eval_all(1)
                  do
                      ef = ef + 1.d-3
                      occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
                      if (abs(occ-1.d0)<1.d-8) exit
                  end do
                  !!write(*,*) 'HACK: INCREASE eF by 0.001d0'
                  !!eF = eF + 0.001d0
                  do ieval=1,ntot
                      ij = ij + 1
                      occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
                      occ_all(ieval) = occ
                  end do
                  !!!if (bigdft_mpi%iproc==0) then
                  !!!    call yaml_sequence_close()
                  !!!    call yaml_map('number of states to be occupied (without smearing)',iq)
                  !!!    call yaml_map('Pseudo Fermi level for occupations',ef)
                  !!!    call yaml_sequence_open('ordered eigenvalues and occupations')
                  !!!    ii = 0
                  !!!    do i=1,ntot
                  !!!        !occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
                  !!!        occ = occ_all(i)
                  !!!        if (.true. .or. occ>1.d-100) then
                  !!!            call yaml_sequence(advance='no')
                  !!!            call yaml_mapping_open(flow=.true.)
                  !!!            call yaml_map('eval',eval_all(i),fmt='(es13.4)')
                  !!!            call yaml_map('atom',id_all(i),fmt='(i5.5)')
                  !!!            call yaml_map('occ',occ,fmt='(1pg13.5e3)')
                  !!!            call yaml_mapping_close(advance='no')
                  !!!            call yaml_comment(trim(yaml_toa(i,fmt='(i5.5)')))
                  !!!        else
                  !!!            ii = ii + 1
                  !!!        end if
                  !!!    end do
                  !!!    if (ii>0) then
                  !!!        call yaml_sequence(advance='no')
                  !!!        call yaml_mapping_open(flow=.true.)
                  !!!        call yaml_map('remaining states',ii)
                  !!!        call yaml_map('occ','<1.d-100')
                  !!!        call yaml_mapping_close()
                  !!!    end if
                  !!!    call yaml_sequence_close()
                  !!!end if
              end if
        
              ! Calculate the projector. First for each single atom, then insert it into the big one.
              charge_total = 0.d0
              do kat=1,natp
                  kkat = kat + isat
                  n = n_all(kat)
                  proj = f_malloc0((/n,n/),id='proj')
                  call calculate_projector(n, ntot, nmax, kkat, id_all, eval_all, &
                       coeff_all(1:,1:,kat), occ_all, proj)
                  if (present(projx)) then
                      call vcopy(n**2, proj(1,1), 1, projx(1,kat), 1)
                  end if
    
    
                  !@ TEMPORARY ############################################
                  ! Extract ktilde
                  ktilde = f_malloc0((/n,n/),id='ktilde')
                  !if (ortho=='yes') then
                      call extract_matrix(smatl, kerneltilde, neighbor(1:,kat), n, nmax, ktilde)
                  !else if (ortho=='no') then
                  !    call extract_matrix(smatl, kernel_%matrix_compr, neighbor(1:,kat), n, nmax, ktilde)
                  !end if
                  kp = f_malloc((/n,n/),id='kp')
                  call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, proj(1,1), n, 0.d0, kp(1,1), n)
                  if (ortho=='no') then
                      call f_memcpy(src=kp,dest=ktilde)
                      ovrlp = f_malloc0((/n,n/),id='ovrlp')
                      call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1:,kat), n, nmax, ovrlp)
                      call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, ovrlp(1,1), n, 0.d0, kp(1,1), n)
                      call f_free(ovrlp)
                  end if
                  tt = 0
                  do i=1,n
                      tt = tt + kp(i,i)
                  end do
                  charge_per_atom(kkat) = tt
                  charge_total = charge_total + tt
                  call f_free(proj)
                  call f_free(ktilde)
                  call f_free(kp)
    
    
              end do
    
    

          if (bigdft_mpi%nproc>1) then
              call mpiallred(charge_per_atom, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(charge_total, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          if (bigdft_mpi%iproc==0) then
              !do iat=1,at%astruct%nat
              !    write(*,*) 'iat, cpa',iat,charge_per_atom(iat)
              !end do
              !write(*,*) 'charge_total',charge_total
          end if
          charge_net = 0.d0
          do iat=1,at%astruct%nat
              charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
          end do
          if (bigdft_mpi%iproc==0) then
              !write(*,*) 'net charge', charge_net
              call yaml_mapping_open(flow=.true.)
              call yaml_map('alpha',alpha,fmt='(es12.4)')
              call yaml_map('net charge',charge_net,fmt='(es12.4)')
              call yaml_map('bisection bounds ok',(/bound_low_ok,bound_up_ok/))
              call yaml_mapping_close()
          end if


          if (abs(charge_net)<convergence_criterion .or. ialpha==10000) then
          !if (.true.) then
              if (bigdft_mpi%iproc==0) then
                  call yaml_sequence_close()
                  call yaml_map('number of states to be occupied (without smearing)',iq)
                  call yaml_map('Pseudo Fermi level for occupations',ef)
                  call yaml_sequence_open('ordered eigenvalues and occupations')
                  ii = 0
                  do i=1,ntot
                      !occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
                      occ = occ_all(i)
                      if (occ>1.d-100) then
                          call yaml_sequence(advance='no')
                          call yaml_mapping_open(flow=.true.)
                          call yaml_map('eval',eval_all(i),fmt='(es13.4)')
                          call yaml_map('atom',id_all(i),fmt='(i5.5)')
                          call yaml_map('occ',occ,fmt='(1pg13.5e3)')
                          call yaml_mapping_close(advance='no')
                          call yaml_comment(trim(yaml_toa(i,fmt='(i5.5)')))
                      else
                          ii = ii + 1
                      end if
                  end do
                  if (ii>0) then
                      call yaml_sequence(advance='no')
                      call yaml_mapping_open(flow=.true.)
                      call yaml_map('remaining states',ii)
                      call yaml_map('occ','<1.d-100')
                      call yaml_mapping_close()
                  end if
                  call yaml_sequence_close()
              end if
              exit alpha_loop
          end if

          ! If we are still searching the boundaries for the bisection...
          if (.not.bound_low_ok) then
              if (charge_net<0.d0) then
                  ! this is a lower bound
                  alpha_low = alpha
                  bound_low_ok = .true.
              else
                  alpha_low = 0.5d0*alpha
              end if
              cycle alpha_loop
          else if (.not.bound_up_ok) then
              if (charge_net>0.d0) then
                  ! this is an upper bound
                  alpha_up = alpha
                  bound_up_ok = .true.
              else
                  alpha_up = 2.0d0*alpha
              end if
              cycle alpha_loop
          end if

          if (charge_net>0.d0) then
              ! Too few electrons, i.e. confinement should be smaller
              !alpha = alpha*0.80
              alpha_up = alpha
          else if (charge_net<0.d0) then
              ! Too many electrons, i.e. confinement should be larger
              !alpha = alpha*1.2
              alpha_low = alpha
          end if



      end do alpha_loop

      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_close()
          call yaml_comment('Projector calculation finished',hfill='~')
      end if

      !call f_free(tmpmat2)
      if (ortho=='yes') then
          call f_free(tmpmat1)
          call deallocate_matrices(ovrlp_onehalf_(1))
      end if
          call f_free(kerneltilde)
      call f_free(coeff_all)
      call f_free(ilup)
      call f_free(n_all)
      call f_free(ovrlp_minusonehalf)
      call f_free(penalty_matrices)
      call f_free(alpha_calc)
    
      !if (bigdft_mpi%iproc==0) then
      !    call yaml_mapping_close()
      !end if
    
    
    
      if (write_output .and. bigdft_mpi%iproc==0) then
          call write_partial_charges(at, charge_per_atom, write_gnuplot=.false.)
      end if

      call f_free(charge_per_atom)
      call f_free(neighbor)
      call f_free(eval_all)
      call f_free(occ_all)
      call f_free(id_all)
      call f_free(ovrlp_onehalf_all)
      call f_free(com)

      call f_release_routine()

  end subroutine projector_for_charge_analysis




  subroutine extract_matrix(smat, matrix_compr, neighbor, n, nmax, matrix, ilup)
    use module_base
    use sparsematrix_base,only: sparse_matrix, matrices
    use sparsematrix_init, only: matrixindex_in_compressed
    implicit none

    ! Calling arguments
    type(sparse_matrix),intent(in) :: smat
    real(kind=8),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: matrix_compr
    logical,dimension(smat%nfvctr),intent(in) :: neighbor
    integer,intent(in) :: n, nmax
    real(kind=8),dimension(n,n),intent(out) :: matrix
    integer,dimension(2,nmax,nmax),intent(out),optional :: ilup

    ! Local variables
    integer :: icheck, ii, jj, i, j, ind
    logical :: optional_present

    call f_routine(id='extract_matrix')

    optional_present = present(ilup)

    icheck = 0
    ii = 0
    do i=1,smat%nfvctr
        if (neighbor(i)) then
            jj = 0
            do j=1,smat%nfvctr
                if (neighbor(j)) then
                    icheck = icheck + 1
                    jj = jj + 1
                    if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                    ind =  matrixindex_in_compressed(smat, j, i)
                    if (ind>0) then
                        matrix(jj,ii) = matrix_compr(ind-smat%isvctrp_tg)
                    else
                        matrix(jj,ii) = 0.d0
                    end if
                    if (optional_present) then
                        ilup(1,jj,ii) = j
                        ilup(2,jj,ii) = i
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
          call initialize_work_arrays_sumrho(1,[lzd%Llr(ilr)],.true.,w)
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




 subroutine unitary_test_multipoles(iproc, nproc, nphi, nphir, orbs, lzd, smat, collcom, hgrids)
   use module_base
   use module_types, only: orbitals_data, comms_linear, local_zone_descriptors, comms_linear
   use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_TASKGROUP, assignment(=), &
                                matrices_null, sparsematrix_malloc_ptr, deallocate_matrices
   use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
   use sparsematrix_init, only: matrixindex_in_compressed
   use yaml_output
   use bounds, only: geocode_buffers
   implicit none
   ! Calling arguments
   integer,intent(in) :: iproc, nproc, nphi, nphir
   type(orbitals_data),intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: lzd
   type(sparse_matrix),intent(in) :: smat
   type(comms_linear),intent(in) :: collcom
   real(kind=8),dimension(3) :: hgrids
   ! Local variables
   integer :: iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, l, m, i, ind, ist, istr, ii, nl1, nl2, nl3
   real(kind=8) :: x, y, z, r2, r, factor, rmax, factor_normalization, val
   real(kind=8),dimension(:),allocatable :: phi2r, phi2, phi1r, phi1
   real(kind=8),dimension(:,:),allocatable :: locregcenter
   type(matrices) :: multipole_matrix
   type(workarr_sumrho) :: w
   real(kind=8),dimension(-lmax:lmax,0:lmax) :: errors
   real(kind=8),dimension(-lmax:lmax,0:lmax) :: values_orig
   real(kind=8),dimension(-lmax:lmax,0:lmax) :: values

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

   ist = 0
   do iorb=1,orbs%norbp
       iiorb = orbs%isorb + iorb
       ilr = orbs%inwhichlocreg(iiorb)
       !rmax = min(lzd%llr(ilr)%d%n1i*0.25d0*hgrids(1),lzd%llr(ilr)%d%n2i*0.25d0*hgrids(2),lzd%llr(ilr)%d%n3i*0.25d0*hgrids(3))
       rmax = min(lzd%llr(ilr)%d%n1*0.5d0*hgrids(1),lzd%llr(ilr)%d%n2*0.5d0*hgrids(2),lzd%llr(ilr)%d%n3*0.5d0*hgrids(3))
       factor_normalization = 3.d0/(4.d0*pi*rmax**3)*0.5d0*lzd%hgrids(1)*0.5d0*lzd%hgrids(2)*0.5d0*lzd%hgrids(3)
       ii = 0
       ! Since the radial function is constant and thus not decaying towards the boundaries of the integration sphere, the center
       ! of the integration volume must be on a gridpoint to avoid truncation artifacts.
       locregcenter(1:3,ilr) = get_closest_gridpoint((/lzd%llr(ilr)%locregcenter(1),&
                                                       lzd%llr(ilr)%locregcenter(2),&
                                                       lzd%llr(ilr)%locregcenter(3)/),&
                                                       hgrids)
       do i3=1,lzd%llr(1)%d%n3i
           ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
           z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3,ilr)
           do i2=1,lzd%llr(1)%d%n2i
               ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
               y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2,ilr)
               do i1=1,lzd%llr(1)%d%n1i
                   ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                   x = ii1*0.5d0*lzd%hgrids(1) - locregcenter(1,ilr)
                   r2 = x**2+y**2+z**2
                   if (r2>rmax**2) cycle
                   !r = sqrt(r2)
                   !r = max(0.5d0,r)
                   ind = (i3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2-1)*lzd%llr(ilr)%d%n1i + i1
                   ii = ii + 1
                   do l=0,lmax
                       do m=-l,l
                           factor = get_test_factor(l,m)*factor_normalization*sqrt(4.d0*pi*real(2*l+1,kind=8))
                           if (l==1) then
                               factor = factor*5.d0/(3.d0*rmax**2)
                           else if (l==2) then
                               factor = factor*7.d0/(3.d0*rmax**4)
                           end if
                           phi2r(ist+ind) = phi2r(ist+ind) + factor*solid_harmonic(0, r, l, m , x, y, z)
                       end do
                   end do
                   !write(*,*) 'i1, i2, i3, ist+ind, val', i1, i2, i3, ist+ind, phi2r(ist+ind)
               end do
           end do
       end do
       ist = ist + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
   end do

   if (nproc>1) then
       call mpiallred(locregcenter, mpi_sum, comm=bigdft_mpi%mpi_comm)
   end if


   ! Set phi1 to 1
   phi1r(:) = 1.d0

   ! Transform back to wavelets
   phi2 = f_malloc0(nphi,id='phi2')
   ist=1
   istr=1
   do iorb=1,orbs%norbp
       iiorb=orbs%isorb+iorb
       ilr=orbs%inwhichlocreg(iiorb)
       call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.true.,w)
       call isf_to_daub(lzd%llr(ilr), w, phi2r(istr), phi2(ist))
       call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.false.,w)
       call isf_to_daub(lzd%llr(ilr), w, phi1r(istr), phi1(ist))
       call deallocate_work_arrays_sumrho(w)
       ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
       istr = istr + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
   end do


   !do ind=1,nphi
   !    write(*,*) 'ind, val', ind, phi2(ind)
   !end do



   do l=0,lmax
       do m=-l,l
           call calculte_multipole_matrix(iproc, nproc, l, m, nphi, phi1, phi2, nphir, hgrids, &
                    orbs, collcom, lzd, smat, locregcenter, 'sphere', multipole_matrix)
           val = 0.d0
           do iorb=1,orbs%norb
               iiorb = modulo(iorb-1,smat%nfvctr)+1
               ind = matrixindex_in_compressed(smat, iorb, iorb)
               val = val + multipole_matrix%matrix_compr(ind)
               !write(*,*) 'l, m, iorb, ind, val', &
               !    l, m, iorb, ind, multipole_matrix%matrix_compr(ind)
           end do
           values(m,l) = val/real(orbs%norb,kind=8)
           errors(m,l) = 100.d0*abs(values(m,l)/get_test_factor(l,m)-1.d0)
           values_orig(m,l) = get_test_factor(l,m)
           !if (iproc==0) write(*,*) 'l, m, val, error', l, m, val, abs(val-get_test_factor(l,m))
       end do
   end do

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
   use sparsematrix_init, only: matrixindex_in_compressed
   use orthonormalization, only: orthonormalizelocalized

   use communications_base, only: TRANSPOSE_FULL
   use transposed_operations, only: calculate_overlap_transposed
   use communications, only: transpose_localized
   use multipole_base, only: external_potential_descriptors, external_potential_descriptors_null, &
                             multipole_set_null, multipole_null, deallocate_external_potential_descriptors
   
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
      call initialize_work_arrays_sumrho(1,[tmb%lzd%Llr(ilr)],.true.,w)
      ! Transform the support function to real space
      call daub_to_isf(tmb%lzd%llr(ilr), w, phi_ortho(ist), phir(istr))
      call initialize_work_arrays_sumrho(1,[tmb%lzd%llr(ilr)],.false.,w)
      ! Transform the functions which is constantly one to wavelets
      call isf_to_daub(tmb%lzd%llr(ilr), w, phi1r(istr), phi1(ist))
      call deallocate_work_arrays_sumrho(w)

      ! NEW: CALCULATE THE WEIGHT CENTER OF THE SUPPORT FUNCTION ############################
      hxh = 0.5d0*tmb%lzd%hgrids(1)
      hyh = 0.5d0*tmb%lzd%hgrids(2)
      hzh = 0.5d0*tmb%lzd%hgrids(3)
      ii = istr
      call geocode_buffers(tmb%lzd%Llr(ilr)%geocode, tmb%lzd%glr%geocode, nl1, nl2, nl3)
      weight = 0.d0
      do i3=1,tmb%lzd%llr(ilr)%d%n3i
          ii3 = tmb%lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
          z = ii3*hzh
          do i2=1,tmb%lzd%llr(ilr)%d%n2i
              ii2 = tmb%lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
              y = ii2*hyh
              do i1=1,tmb%lzd%llr(ilr)%d%n1i
                  ii1 = tmb%lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                  x = ii1*hxh
                  tt = phir(ii)**2
                  center_locreg(1,ilr) = center_locreg(1,ilr) + x*tt
                  center_locreg(2,ilr) = center_locreg(2,ilr) + y*tt
                  center_locreg(3,ilr) = center_locreg(3,ilr) + z*tt
                  weight = weight + tt
                  ii = ii + 1
              end do
          end do
      end do
      !write(*,*) 'iorb, weight, sum(phir)',iorb, weight, sum(phir)
      center_locreg(1:3,ilr) = center_locreg(1:3,ilr)/weight
      center_orb(1:3,iiorb) = center_locreg(1:3,ilr)
      !write(*,*) 'iorb, ilr, center_locreg(1:3,ilr), lzd%llr(ilr)%locregcenter(1:3)', &
      !            iorb, ilr, center_locreg(1:3,ilr), tmb%lzd%llr(ilr)%locregcenter(1:3)
      ! ######################################################################################
      ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i
  end do

  if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
      stop
  end if

  if (nproc>1) then
      call mpiallred(center_locreg, mpi_sum, comm=bigdft_mpi%mpi_comm)
      call mpiallred(center_orb, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

  factor = hxh*hyh*hzh

  do l=0,lmax
      do m=-l,l
          call f_zero(multipole_matrix%matrix_compr)
          ! Calculate the multipole matrix
          call calculte_multipole_matrix(iproc, nproc, l, m, tmb%npsidim_orbs, phi1, phi_ortho, &
               max(tmb%collcom_sr%ndimpsi_c,1), tmb%lzd%hgrids, &
               tmb%orbs, tmb%collcom, tmb%lzd, tmb%linmat%s, center_locreg, 'box', multipole_matrix)
          !write(*,*) 'multipole_matrix%matrix_compr(1)',multipole_matrix%matrix_compr(1)
          ! Take the diagonal elements and scale by factor (anyway there is no really physical meaning in the actual numbers)
          do iorb=1,tmb%orbs%norbp
              iiorb = tmb%orbs%isorb + iorb
              ind = matrixindex_in_compressed(tmb%linmat%s, iiorb, iiorb)
              multipoles(m,l,iiorb) = multipole_matrix%matrix_compr(ind)*factor
              !write(*,*) 'iorb, multipoles(:,:,iiorb)',iorb, multipoles(:,:,iiorb)
          end do
      end do
  end do

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
      call yamL_map('Orthonormalization',do_ortho)
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
      call write_multipoles_new(ep, atoms%astruct%units, &
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


 subroutine get_optimal_sigmas(iproc, nproc, nsigma, collcom_sr, smatl, kernel_, at, lzd, ep, shift, rxyz, ixc, denspot)
   use module_base
   use module_types, only: DFT_wavefunction, input_variables, DFT_local_fields, comms_linear, DFT_local_fields, &
                           local_zone_descriptors
   use sparsematrix_base, only: sparse_matrix, matrices
   use module_atoms, only: atoms_data
   use Poisson_Solver, only: H_potential
   use rhopotential, only: sumrho_for_TMBs, corrections_for_negative_charge
   use yaml_output
   implicit none
   ! Calling arguments
   integer,intent(in) :: iproc, nproc, nsigma, ixc
   type(comms_linear),intent(inout) :: collcom_sr
   type(sparse_matrix),intent(in) :: smatl
   type(matrices),intent(in) :: kernel_
   type(atoms_data),intent(in) :: at
   type(local_zone_descriptors),intent(in) :: lzd
   type(external_potential_descriptors),intent(in) :: ep
   real(kind=8),dimension(3),intent(in) :: shift
   real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
   type(DFT_local_fields),intent(inout) :: denspot
   ! Local variables
   real(kind=8),dimension(:,:,:,:),allocatable :: test_pot
   logical :: rho_negative, exists, found, all_norms_ok
   real(kind=8) :: ehart_ps, diff, tt, diff_min, diff_dipole, diff_dipole_min, rdim
   !integer,parameter :: nsigma=3
   real(kind=8),parameter :: step=0.20d0
   integer :: i1, i2, i3, isigma0, isigma1, isigma2, impl, l
   integer :: nzatom, nelpsp, npspcode, itype, ioffset, ishift
   real(gp),dimension(0:4,0:6) :: psppar
   real(kind=8),dimension(:,:),allocatable :: sigmax
   real(kind=8),dimension(0:lmax) :: factor, factorx, factor_min
   real(kind=8),dimension(:),allocatable :: rhov_orig
   real(kind=8),dimension(3) :: dipole_exact, dipole_trial
   real(kind=8) :: rloc
   integer,dimension(:),allocatable :: psp_source

   call f_routine(id='get_optimal_sigmas')

   if (iproc==0) call yaml_comment('Determine optimal sigmas for the radial Gaussians',hfill='~')

   test_pot = f_malloc0((/size(denspot%V_ext,1),size(denspot%V_ext,2),size(denspot%V_ext,3),2/),id='test_pot')
   rhov_orig = f_malloc(size(denspot%rhov),id='rhov_orig')

   ! Keep the original value fo rhov, which contains the entire potential
   call f_memcpy(src=denspot%rhov, dest=rhov_orig)

   ! Calculate the correct electrostatic potential, i.e. electronic plus ionic part
   call sumrho_for_TMBs(iproc, nproc, lzd%hgrids(1), lzd%hgrids(2), lzd%hgrids(3), &
        collcom_sr, smatl, kernel_, &
        denspot%dpbox%ndimrhopot, &
        denspot%rhov, rho_negative)
   if (rho_negative) then
       call corrections_for_negative_charge(iproc, nproc, at, denspot)
   end if

   denspot%rho_work = f_malloc_ptr(denspot%dpbox%ndimrhopot,id='denspot%rho_work')
   ioffset=lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%i3xcsh
   if (denspot%dpbox%ndimrhopot>0) then
       call vcopy(denspot%dpbox%ndimpot,denspot%rhov(ioffset+1),1,denspot%rho_work(1),1)
       ! add the spin down part if present
       if (denspot%dpbox%nrhodim==2) then
           ishift=denspot%dpbox%ndimrhopot/denspot%dpbox%nrhodim !start of the spin down part
           call axpy(denspot%dpbox%ndimpot, 1.d0, &
                     denspot%rhov(ioffset+ishift+1), &
                     1, denspot%rho_work(1),1)
       end if
   end if
   !do l=1,size(denspot%rho_work)
   !    write(100,*) denspot%rho_work(l)
   !end do
   !write(*,*) 'calculate dipole with rho_work'
   call calculate_dipole_moment(denspot%dpbox, 1, at, rxyz, denspot%rho_work, &
        calculate_quadrupole=.true., dipole=dipole_exact, quiet_=.true.)
   call f_free_ptr(denspot%rho_work)

   call H_potential('D',denspot%pkernel,denspot%rhov,denspot%V_ext,ehart_ps,0.0_dp,.true.,&
        quiet=denspot%PSquiet)!,rho_ion=denspot%rho_ion)
   call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
        denspot%rhov(1), 1, test_pot(1,1,1,1), 1)

   ! Get an initial guess for the sigmas (use rloc from the pseudopotential)
   sigmax = f_malloc((/0.to.lmax,1.to.ep%nmpl/),id='sigmax')
   psp_source = f_malloc(ep%nmpl,id='psp_source')
   do impl=1,ep%nmpl
       !ixc = 1
       !if (iproc==0) then
       !    call yaml_warning('WARNING: USE ixc = 1 IN GET_OPTIMAL_SIGMAS')
       !end if
       !call psp_from_data(ep%mpl(impl)%sym, nzatom, nelpsp, npspcode, ixc, psppar, exists)
       !if (.not.exists) then
       !    call f_err_throw('No PSP available for external multipole type '//trim(ep%mpl(impl)%sym), &
       !         err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       !end if
       !ep%mpl(impl)%sigma(0:lmax) = psppar(0,0)-min(0.9d0,step*real(nsigma/2,kind=8))*psppar(0,0)
       !sigmax(0:lmax,impl) = psppar(0,0)
       !!found = .false.
       !!search_loop: do itype=1,at%astruct%ntypes
       !!    if (trim(ep%mpl(impl)%sym)==trim(at%astruct%atomnames(itype))) then
       !!        sigmax(0:lmax,impl) = 1.0d0*at%psppar(0,0,itype)
       !!        found = .true.
       !!        exit search_loop
       !!    end if
       !!end do search_loop
       !!if (.not.found) then
       !!    call f_err_throw('No PSP available for external multipole type '//trim(ep%mpl(impl)%sym), &
       !!         err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       !!end if
       call get_psp_info(trim(ep%mpl(impl)%sym), ixc, at, nelpsp, psp_source(impl), rloc)
       sigmax(0:lmax,impl) = 1.d0*rloc
   end do
   if (iproc==0) call write_psp_source(ep, psp_source)
   call f_free(psp_source)

   if (iproc==0) then
       call yaml_sequence_open('Determine optimal sigmas')
   end if
   factorx(0:lmax) = max(0.1d0,1.d0-step*real(nsigma/2,kind=8))
   ! The following loops are designed for lmax=2... stop otherwise
   if (lmax>2) then
       call f_err_throw('the maximal lmax possible is 2, but here we have '//trim(yaml_toa(lmax)),&
            err_name='BIGDFT_RUNTIME_ERROR')
   end if
   diff_min = huge(diff_min)
   diff_dipole_min = huge(diff_dipole_min)
   factor_min(0:lmax) = 1.d0 !initialization
   do isigma2=1,nsigma
       do isigma1=1,nsigma
           do isigma0=1,nsigma
               factor(0) = factorx(0) + real(isigma0-1,kind=8)*step
               factor(1) = factorx(1) + real(isigma1-1,kind=8)*step
               factor(2) = factorx(2) + real(isigma2-1,kind=8)*step
               do impl=1,ep%nmpl
                   !ep%mpl(impl)%sigma(l) = ep%mpl(impl)%sigma(l) + step
                   ep%mpl(impl)%sigma(0:lmax) = sigmax(0:lmax,impl)*factor(0:lmax)
                   !if (iproc==0) write(*,*) 'impl, sigma', impl, ep%mpl(impl)%sigma(0:lmax)
               end do
               call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
                    denspot%V_ext(1,1,1,1), 1, test_pot(1,1,1,2), 1)
               call potential_from_charge_multipoles(iproc, nproc, at, denspot, ep, 1, &
                    denspot%dpbox%ndims(1), 1, denspot%dpbox%ndims(2), &
                    denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
                    denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
                    denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
                    denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3), &
                    shift, verbosity=0, ixc=ixc, lzd=lzd, pot=test_pot(:,:,:,2), &
                    rxyz=rxyz, dipole_total=dipole_trial, all_norms_ok=all_norms_ok)

               if (all_norms_ok) then
                   diff_dipole = (dipole_exact(1)-dipole_trial(1))**2 + &
                                 (dipole_exact(2)-dipole_trial(2))**2 + &
                                 (dipole_exact(3)-dipole_trial(3))**2
                   rdim = 1.d0/(real(size(denspot%V_ext,1),kind=8)*&
                                real(size(denspot%V_ext,1),kind=8)*&
                                real(size(denspot%V_ext,1),kind=8))
                   diff = 0.d0
                   do i3=1,size(denspot%V_ext,3)
                       do i2=1,size(denspot%V_ext,2)
                           do i1=1,size(denspot%V_ext,1)
                               !write(800,*) 'i1, i2, i3, vals', i1, i2, i3, test_pot(i1,i2,i3,1), test_pot(i1,i2,i3,2)
                               diff = diff + rdim*(test_pot(i1,i2,i3,1)-test_pot(i1,i2,i3,2))**2
                           end do
                       end do
                   end do
                   call mpiallred(diff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
                   !!tt = diff/(real(size(denspot%V_ext,1),kind=8)*&
                   !!           real(size(denspot%V_ext,1),kind=8)*&
                   !!           real(size(denspot%V_ext,1),kind=8))
               end if
               if (iproc==0) then
                   call yaml_sequence(advance='no')
                   call yaml_mapping_open(flow=.true.)
                   call yaml_map('rloc mult',factor,fmt='(f3.1)')
                   call yaml_map('Gaussian norms ok',all_norms_ok)
                   if (all_norms_ok) then
                       call yaml_map('dipole norm diff (actual/min)',(/diff_dipole,diff_dipole_min/),fmt='(es9.3)')
                       call yaml_map('avg pot diff (actual/min)',(/diff,diff_min/),fmt='(es9.3)')
                   end if
                   call yaml_mapping_close()
               end if
               if (all_norms_ok) then
                   if (diff<diff_min) then
                       !factor_min(0:lmax) = factor(0:lmax)
                       diff_min = diff
                   end if
                   if (diff_dipole<diff_dipole_min) then
                       factor_min(0:lmax) = factor(0:lmax)
                       diff_dipole_min = diff_dipole
                   end if
               end if
           end do
       end do
   end do
   if (iproc==0) then
       call yaml_sequence_close()
   end if
   if (iproc==0) call yaml_map('optimal sigma multiplication factors',factor_min,fmt='(f4.2)')
   do impl=1,ep%nmpl
       ep%mpl(impl)%sigma(0:lmax) = sigmax(0:lmax,impl)*factor_min(0:lmax)
   end do

   call f_memcpy(src=rhov_orig, dest=denspot%rhov)

   call f_free(sigmax)
   call f_free(test_pot)
   call f_free(rhov_orig)

   call f_release_routine()

 end subroutine get_optimal_sigmas

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

   call f_routine(id='calculate_norm')

   call f_zero(norm)

   !$omp parallel default(none) &
   !$omp shared(ep, is1, ie1, is2, ie2, is3, ie3, norm, hhh) &
   !$omp shared(gaussians1, gaussians2, gaussians3) &
   !$omp private(impl, i1, i2, i3, ii1, ii2, ii3, l, gg) 
   !$omp do
   do impl=1,ep%nmpl
       do i3=is3,ie3
           ii3 = i3 - 15
           do i2=is2,ie2
               ii2 = i2 - 15
               do i1=is1,ie1
                   ii1 = i1 - 15
                   do l=0,lmax
                       ! Calculate the Gaussian as product of three 1D Gaussians
                       gg = gaussians1(l,i1,impl)*gaussians2(l,i2,impl)*gaussians3(l,i3,impl)
                       norm(l,impl) = norm(l,impl) + gg*hhh
                   end do
               end do
           end do
       end do
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
  implicit none
  integer, intent(in) :: nspin
  type(denspot_distribution), intent(in) :: dpbox
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(dpbox%ndims(1),dpbox%ndims(2),max(dpbox%n3p, 1),nspin), target, intent(in) :: rho
  !!!!logical :: is_net_charge !< true if the charge density is already the net charge (i.e. including the compensating core charge)
  logical,intent(in) :: calculate_quadrupole
  real(kind=8),dimension(3),intent(out),optional :: dipole
  real(kind=8),dimension(3,3),intent(out),optional :: quadrupole
  logical,intent(in),optional :: quiet_

  integer :: ierr,n3p,nc1,nc2,nc3, nnc3, ii3, i3shift
  real(gp) :: q,qtot, delta_term,x,y,z,ri,rj
  integer  :: iat,i1,i2,i3, nl1,nl2,nl3, ispin,n1i,n2i,n3i, i, j, is, ie
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
  
  n1i=dpbox%ndims(1)
  n2i=dpbox%ndims(2)
  n3i=dpbox%ndims(3)
  n3p=dpbox%n3p


  if (at%astruct%geocode /= 'F') then
     nl1=1
     !nl3=1
     nc1=n1i
     !nc3=n3i
     nc3=n3p
     nnc3=n3i
     !is = 1
     is = dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+1
     ie = dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+dpbox%nscatterarr(dpbox%mpi_env%iproc,2)
     i3shift = 1
  else
     nl1=15
     !nl3=15
     !nl3=max(1,15-dpbox%nscatterarr(dpbox%mpi_env%iproc,3))
     nc1=n1i-31
     !nc3=n3i-31
     !nc3=n3p-31
     is = max(dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+1,15)
     ie = min(dpbox%nscatterarr(dpbox%mpi_env%iproc,3)+dpbox%nscatterarr(dpbox%mpi_env%iproc,2),n3i-17)
     nnc3=n3i-31
     i3shift = 15
     !write(*,*) 'iproc, is, ie, nl3, nc3, n3p', bigdft_mpi%iproc, is, ie, nl3, nc3, n3p
  end if
  nc3 = ie - is + 1 !number of z planes to be treated
  nl3=max(1,i3shift-dpbox%nscatterarr(dpbox%mpi_env%iproc,3)) !offset within rho array
  !value of the buffer in the y direction
  if (at%astruct%geocode == 'P') then
     nl2=1
     nc2=n2i
  else
     nl2=15
     nc2=n2i-31
  end if

  qtot=0.d0
  dipole_cores(1:3)=0._gp
  do iat=1,at%astruct%nat
     !write(*,*) 'iat, rxyz(1:3,iat)',iat, rxyz(1:3,iat)
     dipole_cores(1:3)=dipole_cores(1:3)+at%nelpsp(at%astruct%iatype(iat)) * rxyz(1:3,iat)
  end do
  !!write(*,*) 'dipole_cores',dipole_cores
  !!write(*,*) 'nc3',nc3

  dipole_el   (1:3)=0._gp
  do ispin=1,nspin
     do i3=0,nc3 - 1
        !ii3 = i3 + dpbox%nscatterarr(dpbox%mpi_env%iproc,3)
        ii3 = i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3) - i3shift !real coordinate, without buffer
        !write(*,*) 'iproc, i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3), ii3', &
        !            bigdft_mpi%iproc, i3+nl3+dpbox%nscatterarr(dpbox%mpi_env%iproc,3), ii3
        do i2=0,nc2 - 1
           do i1=0,nc1 - 1
              !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              !q= ( ele_rho(ind,ispin) ) * hxh*hyh*hzh 
              !q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
              q= - rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
              !write(*,*) 'i1, i2, i3, nl1, nl2, nl3, q', i1, i2, i3, nl1, nl2, nl3, q
              qtot=qtot+q
              dipole_el(1)=dipole_el(1)+ q* at%astruct%cell_dim(1)/real(nc1,dp)*i1 
              dipole_el(2)=dipole_el(2)+ q* at%astruct%cell_dim(2)/real(nc2,dp)*i2
              dipole_el(3)=dipole_el(3)+ q* at%astruct%cell_dim(3)/real(nnc3,dp)*ii3
           end do
        end do
     end do
  !!write(*,*) 'iproc, dipole_el,sum(rho), qtot',bigdft_mpi%iproc,dipole_el,sum(rho), qtot
  end do

  !!call mpi_barrier(mpi_comm_world,ispin)
  call mpiallred(qtot, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
  call mpiallred(dipole_el, mpi_sum, comm=bigdft_mpi%mpi_comm)
  !!call mpi_barrier(mpi_comm_world,ispin)
  !!write(*,*) 'after allred: iproc, dipole_el,sum(rho), qtot',bigdft_mpi%iproc,dipole_el,sum(rho), qtot

  !!write(*,*) 'dipole_cores first', dipole_cores
  !!call mpi_barrier(mpi_comm_world,ispin)


  if (calculate_quadrupole) then
      ! Quadrupole not yet parallelized

      if (at%astruct%geocode /= 'F') then
         nl3=1
         nc3=n3i
      else
         nl3=15
         nc3=n3i-31
      end if

      if (dpbox%mpi_env%nproc > 1) then
         !allocate full density in pot_ion array
         ele_rho = f_malloc_ptr((/ n1i, n2i, n3i, nspin /),id='ele_rho')
    
    !Commented out, it is enough to allocate the rho at 1
    !!$     ! rho_buf is used instead of rho for avoiding the case n3p=0 in 
    !!$     ! some procs which makes MPI_ALLGATHERV failed.
    !!$     if (n3p.eq.0) then
    !!$       allocate(rho_buf(n1i,n2i,n3p+1,nspin),stat=i_stat)
    !!$       call memocc(i_stat,rho_buf,'rho_buf',subname)
    !!$       rho_buf = 0.0_dp
    !!$     else
    !!$       allocate(rho_buf(n1i,n2i,n3p,nspin),stat=i_stat)
    !!$       call memocc(i_stat,rho_buf,'rho_buf',subname)
    !!$       rho_buf = rho
    !!$     endif  
    
    
         do ispin=1,nspin
            call MPI_ALLGATHERV(rho(1,1,1,ispin),n1i*n2i*n3p,&
                 mpidtypd,ele_rho(1,1,1,ispin),dpbox%ngatherarr(0,1),&
                 dpbox%ngatherarr(0,2),mpidtypd,dpbox%mpi_env%mpi_comm,ierr)
            !write(*,*) 'dpbox%ngatherarr(:,1)',dpbox%ngatherarr(:,1)
            !write(*,*) 'dpbox%ngatherarr(:,2)',dpbox%ngatherarr(:,2)
            !write(*,*) 'dpbox%nscatterarr(:,2)',dpbox%nscatterarr(:,2)
            !write(*,*) 'dpbox%nscatterarr(:,3)',dpbox%nscatterarr(:,3)
         end do
    
      else
         ele_rho => rho
      end if

      ! charge center
      charge_center_cores(1:3)=0.d0
      qtot=0.d0
      do iat=1,at%astruct%nat
          q=at%nelpsp(at%astruct%iatype(iat))
          charge_center_cores(1:3) = charge_center_cores(1:3) + q*at%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_cores=charge_center_cores/qtot


      quadropole_cores(1:3,1:3)=0._gp
      do iat=1,at%astruct%nat
          do i=1,3
              select case (i)
              case (1)
                  !ri=at%astruct%rxyz(1,iat)-charge_center_cores(1)
                  ri=at%astruct%rxyz(1,iat)
              case (2)
                  !ri=at%astruct%rxyz(2,iat)-charge_center_cores(2)
                  ri=at%astruct%rxyz(2,iat)
              case (3)
                  !ri=at%astruct%rxyz(3,iat)-charge_center_cores(3)
                  ri=at%astruct%rxyz(3,iat)
              case default
                  stop 'wrong value of i'
              end select
              do j=1,3
                  select case (j)
                  case (1)
                      !rj=at%astruct%rxyz(1,iat)-charge_center_cores(1)
                      rj=at%astruct%rxyz(1,iat)
                  case (2)
                      !rj=at%astruct%rxyz(2,iat)-charge_center_cores(2)
                      rj=at%astruct%rxyz(2,iat)
                  case (3)
                      !rj=at%astruct%rxyz(3,iat)-charge_center_cores(3)
                      rj=at%astruct%rxyz(3,iat)
                  case default
                      stop 'wrong value of j'
                  end select
                  if (i==j) then
                      !delta_term = (at%astruct%rxyz(1,iat)-charge_center_cores(1))**2 + &
                      !             (at%astruct%rxyz(2,iat)-charge_center_cores(2))**2 + &
                      !             (at%astruct%rxyz(3,iat)-charge_center_cores(3))**2
                      delta_term = at%astruct%rxyz(1,iat)**2 + &
                                   at%astruct%rxyz(2,iat)**2 + &
                                   at%astruct%rxyz(3,iat)**2
                  else
                      delta_term=0.d0
                  end if
                  q=at%nelpsp(at%astruct%iatype(iat))
                  quadropole_cores(j,i) = quadropole_cores(j,i) + q*(3.d0*rj*ri-delta_term)
              end do
          end do
      end do

      !!ele_rho=0.d0
      !!do iat=1,at%astruct%nat
      !!    i1=nint((at%astruct%rxyz(1,iat)/at%astruct%cell_dim(1))*real(nc1,dp))+nl1
      !!    i2=nint((at%astruct%rxyz(2,iat)/at%astruct%cell_dim(2))*real(nc2,dp))+nl2
      !!    i3=nint((at%astruct%rxyz(3,iat)/at%astruct%cell_dim(3))*real(nc3,dp))+nl3
      !!    if (bigdft_mpi%iproc==0) write(*,*) 'iat,i1,i2,i3',iat,i1,i2,i3
      !!    ele_rho(i1,i2,i3,1)=real(at%nelpsp(at%astruct%iatype(iat)))/product(dpbox%hgrids)
      !!end do


      ! charge center
      charge_center_elec(1:3,1:nspin)=0.d0
      do ispin=1,nspin
          qtot=0.d0
          do i3=0,nc3 - 1
              do i2=0,nc2 - 1
                  do i1=0,nc1 - 1
                      q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
                      x=at%astruct%cell_dim(1)/real(nc1,dp)*i1
                      y=at%astruct%cell_dim(2)/real(nc2,dp)*i2
                      z=at%astruct%cell_dim(3)/real(nc3,dp)*i3
                      charge_center_elec(1,ispin) = charge_center_elec(1,ispin) + q*x
                      charge_center_elec(2,ispin) = charge_center_elec(2,ispin) + q*y
                      charge_center_elec(3,ispin) = charge_center_elec(3,ispin) + q*z
                      qtot=qtot+q
                  end do
              end do
          end do
          !!write(*,*) 'qtot',qtot
          charge_center_elec(1:3,ispin)=charge_center_elec(1:3,ispin)/qtot
      end do

      quadropole_el(1:3,1:3)=0._gp
      do ispin=1,nspin
          do i3=0,nc3 - 1
              do i2=0,nc2 - 1
                  do i1=0,nc1 - 1
                      q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
                      x=at%astruct%cell_dim(1)/real(nc1,dp)*i1
                      y=at%astruct%cell_dim(2)/real(nc2,dp)*i2
                      z=at%astruct%cell_dim(3)/real(nc3,dp)*i3
                      do i=1,3
                          select case (i)
                          case (1)
                              !ri=x-charge_center_cores(1)
                              ri=x+(charge_center_cores(1)-charge_center_elec(1,ispin))
                          case (2)
                              !ri=y-charge_center_cores(2)
                              ri=y+(charge_center_cores(2)-charge_center_elec(2,ispin))
                          case (3)
                              !ri=z-charge_center_cores(3)
                              ri=z+(charge_center_cores(3)-charge_center_elec(3,ispin))
                          case default
                              stop 'wrong value of i'
                          end select
                          do j=1,3
                              select case (j)
                              case (1)
                                  !rj=x-charge_center_cores(1)
                                  rj=x+(charge_center_cores(1)-charge_center_elec(1,ispin))
                              case (2)
                                  !rj=y-charge_center_cores(2)
                                  rj=y+(charge_center_cores(2)-charge_center_elec(2,ispin))
                              case (3)
                                  !rj=z-charge_center_cores(3)
                                  rj=z+(charge_center_cores(3)-charge_center_elec(3,ispin))
                              case default
                                  stop 'wrong value of j'
                              end select
                              if (i==j) then
                                  !delta_term = (x-charge_center_cores(1))**2 + &
                                  !             (y-charge_center_cores(2))**2 + &
                                  !             (z-charge_center_cores(3))**2
                                  delta_term = (x+(charge_center_cores(1)-charge_center_elec(1,ispin)))**2 + &
                                               (y+(charge_center_cores(2)-charge_center_elec(2,ispin)))**2 + &
                                               (z+(charge_center_cores(3)-charge_center_elec(3,ispin)))**2
                              else
                                  delta_term=0.d0
                              end if
                              quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
                          end do
                      end do
                  end do
              end do
          end do
      end do

      !!if (.not. is_net_charge) then
          tmpquadrop=quadropole_cores+quadropole_el
      !!else
      !!    tmpquadrop=quadropole_el
      !!end if

      if (present(quadrupole)) then
          quadrupole = tmpquadrop
      end if

      if (dpbox%mpi_env%nproc > 1) then
         call f_free_ptr(ele_rho)
      else
         nullify(ele_rho)
      end if

  end if

  !!write(*,*) 'dipole_cores second', dipole_cores
  !!call mpi_barrier(mpi_comm_world,ispin)

  !!write(*,*) 'dipole_cores', dipole_cores
  !!call mpi_barrier(mpi_comm_world,ispin)
  !!write(*,*) 'after cores'
  !!write(*,*) 'dipole_el', dipole_el
  !!call mpi_barrier(mpi_comm_world,ispin)
  !!write(*,*) 'after el'

  !!if (.not.is_net_charge) then
      tmpdip=dipole_cores+dipole_el
  !!else
  !!    tmpdip=dipole_el
  !!end if
  !!write(*,*) 'tmpdip before',tmpdip
  !!call mpi_barrier(mpi_comm_world,ispin)
  !!write(*,*) 'tmpdip',tmpdip
  if (present(dipole)) dipole(1:3) = tmpdip(1:3)
  if(bigdft_mpi%iproc==0 .and. .not.quiet) then
     !dipole_el=dipole_el        !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     !dipole_cores=dipole_cores  !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     !write(*,*) 'dipole_cores', dipole_cores
     !write(*,*) 'dipole_el', dipole_el
     call yaml_mapping_open('Electric Dipole Moment (AU)')
       call yaml_map('P vector',tmpdip(1:3),fmt='(1pe13.4)')
       call yaml_map('norm(P)',sqrt(sum(tmpdip**2)),fmt='(1pe14.6)')
     call yaml_mapping_close()
     tmpdip=tmpdip/0.393430307_gp  ! au2debye              
     call yaml_mapping_open('Electric Dipole Moment (Debye)')
       call yaml_map('P vector',tmpdip(1:3),fmt='(1pe13.4)')
       call yaml_map('norm(P)',sqrt(sum(tmpdip**2)),fmt='(1pe14.6)')
     call yaml_mapping_close()


      if (calculate_quadrupole) then
          !call yaml_sequence_open('core quadropole')
          !do i=1,3
          !   call yaml_sequence(trim(yaml_toa(quadropole_cores(i,1:3),fmt='(es15.8)')))
          !end do
          !call yaml_sequence_close()

          !call yaml_sequence_open('electronic quadropole')
          !do i=1,3
          !   call yaml_sequence(trim(yaml_toa(quadropole_el(i,1:3),fmt='(es15.8)')))
          !end do
          !call yaml_sequence_close()

          !!call yaml_sequence_open('Quadrupole Moment (AU)')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(tmpquadrop(i,1:3),fmt='(es15.8)')))
          !!end do
          !!call yaml_map('trace',tmpquadrop(1,1)+tmpquadrop(2,2)+tmpquadrop(3,3),fmt='(es12.2)')
          !!call yaml_sequence_close()
          !call yaml_sequence_open('Quadrupole Moment (AU)')
          call yaml_mapping_open('Quadrupole Moment (AU)')
            call yaml_map('Q matrix',tmpquadrop,fmt='(1pe13.4)')
          !do i=1,3
          !   call yaml_sequence(trim(yaml_toa(tmpquadrop(i,1:3),fmt='(es15.8)')))
          !end do
           call yaml_map('trace',tmpquadrop(1,1)+tmpquadrop(2,2)+tmpquadrop(3,3),fmt='(es12.2)')
          !call yaml_sequence_close()
          call yaml_mapping_close()
      end if

  end if

  !call mpi_barrier(mpi_comm_world,ispin)
  !write(*,*) 'end calculate_dipole_moment'


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
      call initialize_work_arrays_sumrho(1,[lzd%Llr(ilr)],.true.,w)
      ! Transform the support function to real space
      call daub_to_isf(lzd%llr(ilr), w, phi(ist), phir(istr))
      call initialize_work_arrays_sumrho(1,[lzd%llr(ilr)],.false.,w)

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


  subroutine get_psp_info(sym, ixc, at, nelpsp, psp_source, rloc)
    use module_base
    use module_atoms, only: atoms_data
    implicit none

    ! Calling arguments
    character(len=*),intent(in) :: sym
    integer,intent(in) :: ixc
    type(atoms_data),intent(in) :: at
    integer,intent(out) :: nelpsp, psp_source
    real(kind=8),intent(out) :: rloc

    ! Local variables
    integer :: itype, ixc_tmp, npspcode, nzatom
    logical :: found, exists
    real(gp),dimension(0:4,0:6) :: psppar

    found = .false.
    search_loop: do itype=1,at%astruct%ntypes
        if (trim(sym)==trim(at%astruct%atomnames(itype))) then
            rloc = at%psppar(0,0,itype)
            nelpsp = at%nelpsp(itype)
            found = .true.
            psp_source = 0
            exit search_loop
        end if
    end do search_loop
    if (.not.found) then
        ixc_tmp = ixc
        call psp_from_data(trim(sym), nzatom, nelpsp, npspcode, ixc_tmp, psppar, exists)
        if (exists) then
            rloc = psppar(0,0)
            found = .true.
            psp_source = 1
        end if
    end if
    if (.not.found) then
        call f_err_throw('No PSP available for external multipole type '//trim(sym), &
             err_name='BIGDFT_INPUT_VARIABLES_ERROR')
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

end module multipole
