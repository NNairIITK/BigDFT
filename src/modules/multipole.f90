module multipole
  use module_base
  use multipole_base, only: external_potential_descriptors,lmax
  implicit none

  private

  !> Public routines
  public :: interaction_multipoles_ions
  public :: potential_from_charge_multipoles
  public :: potential_from_multipoles
!  public :: multipoles_from_density
  public :: ionic_energy_of_external_charges
  public :: multipole_analysis_core !should better be private...
 !public :: write_multipoles_new !should better be private
  public :: gaussian_density
  public :: multipole_analysis_driver
  public :: projector_for_charge_analysis
  public :: support_function_gross_multipoles

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
    subroutine potential_from_charge_multipoles(iproc, nproc, at, denspot, ep, is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, shift, pot)
      use module_types, only: DFT_local_fields
      use Poisson_Solver, except_dp => dp, except_gp => gp
      use module_atoms, only: atoms_data
      use bounds, only: ext_buffers
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(atoms_data),intent(in) :: at
      type(DFT_local_fields),intent(inout) :: denspot
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3
      real(gp),intent(in) :: hx, hy, hz
      real(gp),dimension(3),intent(in) :: shift !< global shift of the atomic positions
      real(gp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l, m, ii, mm, nthread, ithread
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp, ehart_ps, tt, ttt, gg, hhh, tt0, tt1, tt2
      real(dp),dimension(3) :: r
      real(dp),dimension(:,:,:),allocatable :: density
      real(dp),dimension(:,:,:,:),allocatable :: density_loc, potential_loc
      real(kind=8),dimension(0:lmax) :: sigma
      real(8),dimension(:),allocatable :: monopole
      real(8),dimension(:,:),allocatable :: norm, dipole, quadrupole, norm_check
      real(kind=8),dimension(:,:,:),allocatable :: gaussians1, gaussians2, gaussians3
      logical,dimension(:),allocatable :: norm_ok
      real(kind=8),parameter :: norm_threshold = 1.d-4
      real(kind=8),dimension(0:lmax) :: max_error
      integer :: ixc
      integer :: nbl1, nbl2, nbl3, nbr1, nbr2, nbr3, n3pi, i3s
      integer,dimension(:),allocatable :: nzatom, nelpsp, npspcode
      real(gp),dimension(:,:,:),allocatable :: psppar
      logical :: exists, perx, pery, perz
      logical,parameter :: use_iterator = .false.
      real(kind=8) :: cutoff, rholeaked, hxh, hyh, hzh, rx, ry, rz, qq, ttl
      real(kind=8),dimension(3) :: center
      integer :: n1i, n2i, n3i
      integer :: nmpx, nmpy, nmpz, ndensity
      real(dp), dimension(:), allocatable  :: mpx,mpy,mpz
      real(kind=8),dimension(:),allocatable :: rmax
      !real(kind=8),parameter :: rmin=3.d-1
      real(kind=8) :: rmin
      integer,dimension(0:lmax) :: error_meaningful
      character(len=20),dimension(0:lmax) :: output_arr
      !$ integer  :: omp_get_thread_num,omp_get_max_threads

      call f_routine(id='potential_from_charge_multipoles')

      ! No need to do this if there are no multipoles given.
      multipoles_if: if (ep%nmpl>0) then

          hhh = hx*hy*hz

          ! Used for the calculations of the solid harmonics, see description there
          rmin = 2.0d0*hhh**(1.d0/3.d0)
    
          !sigma(0) = 5.d0*hhh**(1.d0/3.d0) !5.d0*hhh**(1.d0/3.d0)
          !sigma(1) = 4.d0*hhh**(1.d0/3.d0)
          !sigma(2) = 2.d0*hhh**(1.d0/3.d0)
          !sigma(0) = 7.d0*hhh**(1.d0/3.d0) !5.d0*hhh**(1.d0/3.d0)
          !sigma(1) = 7.d0*hhh**(1.d0/3.d0)
          !sigma(2) = 7.d0*hhh**(1.d0/3.d0)
          sigma(0) = 1.0d0
          sigma(1) = 0.8d0
          sigma(2) = 0.6d0
    
          density = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='density')
          
          nthread = 1
          !$ nthread = omp_get_max_threads()
          density_loc = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3,0.to.nthread-1/),id='density')
          potential_loc = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3,0.to.nthread-1/),id='potential_loc')
    
          gaussians1 = f_malloc((/0.to.lmax,1.to.ep%nmpl,is1.to.ie1/),id='gaussians1')
          gaussians2 = f_malloc((/0.to.lmax,1.to.ep%nmpl,is2.to.ie2/),id='gaussians2')
          gaussians3 = f_malloc((/0.to.lmax,1.to.ep%nmpl,is3.to.ie3/),id='gaussians3')
    
          !$omp parallel default(none) &
          !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, shift, ep, sigma) &
          !$omp shared(gaussians1, gaussians2, gaussians3) &
          !$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, tt, l,impl)
          !$omp do
          do i3=is3,ie3
              ii3 = i3 - 15
              z = real(ii3,kind=8)*hz + shift(3)
              do impl=1,ep%nmpl
                  tt = z - ep%mpl(impl)%rxyz(3)
                  tt = tt**2
                  do l=0,lmax
                      gaussians3(l,impl,i3) = gaussian(sigma(l),tt)
                  end do
              end do
          end do
          !$omp end do
          !$omp do
          do i2=is2,ie2
              ii2 = i2 - 15
              y = real(ii2,kind=8)*hy + shift(2)
              do impl=1,ep%nmpl
                  tt = y - ep%mpl(impl)%rxyz(2)
                  tt = tt**2
                  do l=0,lmax
                      ! Undo the normalization for this Gaussian
                      gaussians2(l,impl,i2) = gaussian(sigma(l),tt)*sqrt(2.d0*pi_param*sigma(l)**2)**3
                  end do
              end do
          end do
          !$omp end do
          !$omp do
          do i1=is1,ie1
              ii1 = i1 - 15
              x = real(ii1,kind=8)*hx + shift(1)
              do impl=1,ep%nmpl
                  tt = x - ep%mpl(impl)%rxyz(1)
                  tt = tt**2
                  do l=0,lmax
                      ! Undo the normalization for this Gaussian
                      gaussians1(l,impl,i1) = gaussian(sigma(l),tt)*sqrt(2.d0*pi_param*sigma(l)**2)**3
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel
    
    
    
          norm = f_malloc((/0.to.2,1.to.ep%nmpl/),id='norm')
          norm_check = f_malloc((/0.to.2,1.to.ep%nmpl/),id='norm_check')
          monopole = f_malloc(ep%nmpl,id='monopole')
          dipole = f_malloc((/3,ep%nmpl/),id='dipole')
          quadrupole = f_malloc((/5,ep%nmpl/),id='quadrupole')
          norm_ok = f_malloc0(ep%nmpl,id='norm_ok')


    
          ! First calculate the norm of the Gaussians for each multipole
          norm = 0.d0
          monopole = 0.d0
          dipole = 0.d0
          quadrupole = 0.d0
          !$omp parallel &
          !$omp default(none) &
          !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, hhh, ep, shift, sigma, nthread, norm_ok) &
          !$omp shared(norm, monopole, dipole, quadrupole, density, density_loc, potential_loc) &
          !$omp shared (gaussians1, gaussians2, gaussians3) &
          !$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, l, gg, m, mm, tt, ttt, ithread) &
          !$omp private(rnrm1, rnrm2, rnrm3, rnrm5)
          ithread = 0
          !$ ithread = omp_get_thread_num()
          if (ithread<0 .or. ithread>nthread-1) then
              !SM: Is it possible to call f_err_throw within OpenMP? Anyway this condition should never be true...
              call f_err_throw('wrong value of ithread',err_name='BIGDFT_RUNTIME_ERROR')
          end if
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
                              gg = gaussians1(l,impl,i1)*gaussians2(l,impl,i2)*gaussians3(l,impl,i3)
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
    
          ! Check whether they are ok.
          do impl=1,ep%nmpl
              norm_ok(impl) = .true.
              do l=0,lmax !
                  if (abs(1.d0-norm(l,impl))>norm_threshold) then
                      norm_ok(impl) = .false.
                  end if
              end do
              !write(*,*) 'impl, norm_ok(impl)', impl, norm_ok(impl)
          end do
    
    
          ! Get the parameters for each multipole, required to compensate for the pseudopotential part
          nzatom = f_malloc(ep%nmpl,id='nzatom')
          nelpsp = f_malloc(ep%nmpl,id='nelpsp')
          npspcode = f_malloc(ep%nmpl,id='npspcode')
          psppar = f_malloc( (/0.to.4,0.to.6,1.to.ep%nmpl/),id='psppar')
          do impl=1,ep%nmpl
              ixc = 1
              if (iproc==0) then
                  call yaml_warning('WARNING: USE ixc = 1 IN POTENTIAL_FROM_CHARGE_MULTIPOLES')
              end if
              call psp_from_data(ep%mpl(impl)%sym, nzatom(impl), nelpsp(impl), npspcode(impl), ixc, psppar(:,:,impl), exists)
              if (.not.exists) then
                  call f_err_throw('No PSP available for external multipole type '//trim(ep%mpl(impl)%sym), &
                       err_name='BIGDFT_INPUT_VARIABLES_ERROR')
              end if
          end do
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
    
         !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
         cutoff=10.0_gp*maxval(psppar(0,0,:))
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
    
    
    
         ! Generate the density that comes from the pseudopotential atoms
         ndensity = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
         do impl=1,ep%nmpl
             ! Only do this if the norm of tha Gaussian is ok; otherwise the analytic
             ! expression using the net monopole is used.
             if(norm_ok(impl)) then
                 ! The following routine needs the shifted positions
                 rx = ep%mpl(impl)%rxyz(1) - shift(1)
                 ry = ep%mpl(impl)%rxyz(2) - shift(2)
                 rz = ep%mpl(impl)%rxyz(3) - shift(3)
                 !write(*,*) 'nelpsp(impl)',nelpsp(impl)
          !write(*,*) 'WARNING: GAUSSIAN_DENSITY COMMENTED!!!'
                 call gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, &
                      rx, ry, rz, &
                      psppar(0,0,impl), nelpsp(impl), at%multipole_preserving, use_iterator, at%mp_isf, &
                      denspot%dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, ndensity, density(is1:,is2:,is3:), rholeaked)
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
          !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, hhh, ep, shift, sigma, nthread, norm_ok) &
          !$omp shared(norm_check, monopole, dipole, quadrupole, density, density_loc, potential_loc) &
          !$omp shared (gaussians1, gaussians2, gaussians3, nelpsp, rmax, rmin) &
          !$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, l, gg, m, mm, tt, ttt, ttl, ithread, center) &
          !$omp private(rnrm1, rnrm2, rnrm3, rnrm5, qq, ii)
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
              do i3=is3,ie3
                  ii3 = i3 - 15
                  z = real(ii3,kind=8)*hz + shift(3)
                  !write(*,*) 'impl, z, ep%mpl(impl)%rxyz(3)', impl, z, ep%mpl(impl)%rxyz(3)
                  do i2=is2,ie2
                      ii2 = i2 - 15
                      y = real(ii2,kind=8)*hy + shift(2)
                      do i1=is1,ie1
                          ii1 = i1 - 15
                          x = real(ii1,kind=8)*hx + shift(1)
                          tt = 0.d0
                          center = nearest_gridpoint(ep%mpl(impl)%rxyz, (/hx,hy,hz/))
                          !write(*,*) 'rxyz, center', ep%mpl(impl)%rxyz, center
                          r(1) = x - ep%mpl(impl)%rxyz(1)
                          r(2) = y - ep%mpl(impl)%rxyz(2)
                          r(3) = z - ep%mpl(impl)%rxyz(3)
                          !r(1) = x - center(1)
                          !r(2) = y - center(2)
                          !r(3) = z - center(3)
                          rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                          if (norm_ok(impl)) then
                              ! Use the method based on the Gaussians
                              ttl = 0.d0
                              do l=0,lmax
                                  !gg = gaussian(sigma(l),rnrm2)
                                  ! Calculate the Gaussian as product of three 1D Gaussians
                                  gg = gaussians1(l,impl,i1)*gaussians2(l,impl,i2)*gaussians3(l,impl,i3)
                                  ! Additional modification to avoid divergence
                                  if (l==1) then
                                      !gg = gg*rnrm2/(3.d0*sigma(l)**2)
                                      gg = gg/(3.d0*sigma(l)**2)
                                  else if (l==2) then
                                      !gg = gg*rnrm2**2/(15.d0*sigma(l)**4)
                                      gg = gg/(15.d0*sigma(l)**4)
                                  end if
                                  !gg = gaussian(sigma(l),r(1)**2)*gaussian(sigma(l),r(2)**2)*gaussian(sigma(l),r(3)**2)
                                  !gg = gg*sqrt(2.d0*pi_param*sigma(l)**2)**6
                                  norm_check(l,impl) = norm_check(l,impl) + gg*hhh*rnrm2**l
                                  if (rnrm2<=rmax(impl)**2) then
                                  !if (.true.) then
                                      if (associated(ep%mpl(impl)%qlm(l)%q)) then
                                          mm = 0
                                          do m=-l,l
                                              mm = mm + 1
                                              !ttt = ep%mpl(impl)%qlm(l)%q(mm)*&
                                              !      spherical_harmonic(l, m, r(1), r(2), r(3))*gg*sqrt(4.d0*pi_param)
                                              ! For the monopole term, the atomic core charge (which has been expressed using a Gaussian
                                              ! above) has to be added in order to compensate it.
                                              if (l==0) then
                                                  qq = ep%mpl(impl)%qlm(l)%q(mm) + real(nelpsp(impl),kind=8)
                                              else
                                                  qq = ep%mpl(impl)%qlm(l)%q(mm)
                                              end if
                                              !ttt = qq*&
                                              !      spherical_harmonic(-2, rmax(impl), l, m, r(1), r(2), r(3))*gg!*sqrt(4.d0*pi_param)
                                              !ttt = qq*&
                                              !      solid_harmonic(-2, rmin, l, m, r(1), r(2), r(3))*gg!*sqrt(4.d0*pi_param)
                                              !ttt = qq*&
                                              !      real(2*l+1,kind=8)*solid_harmonic(-2, rmin, l, m, r(1), r(2), r(3))*&
                                              !      sqrt(4.d0*pi/real(2*l+1,kind=8))*gg!*sqrt(4.d0*pi_param)
                                              ttt = qq*&
                                                    real(2*l+1,kind=8)*solid_harmonic(0, rmin, l, m, r(1), r(2), r(3))*&
                                                    sqrt(4.d0*pi/real(2*l+1,kind=8))*gg!*sqrt(4.d0*pi_param)
                                              !ttt = qq*&
                                              !      spherical_harmonic(-2, rmax(impl), l, m, r(1), r(2), r(3))*gg*4.d0*pi_param
                                              !ttt = qq*&
                                              !      spherical_harmonic(-2, rmax(impl), l, m, r(1), r(2), r(3))**2
                                              !if (l==0) then
                                              !    ttt = ttt/sqrt(3.d0)
                                              !else if (l==1) then
                                              !    ttt = ttt/sqrt(5.d0)
                                              !else if (l==2) then
                                              !    ttt = ttt/sqrt(7.d0)
                                              !end if
                                              !ttt = ttt/get_normalization(rmax(impl), l, m)!*(0.5d0*sqrt(1/pi_param))
                                              ttt = ttt
                                              !ttt = ttt*get_normalization(rmax(impl), l, m)*4.d0*pi_param!*(0.5d0*sqrt(1/pi_param))
                                              !ttt = ttt*get_normalization(rmax(impl), l, m)!*sqrt(4.d0*pi_param)
                                              !ttt = ttt**2!*get_normalization(rmax(impl), l, m)!*sqrt(4.d0*pi_param)
                                              tt = tt + ttt
                                              ttl = ttl + ttt
                                          end do
                                      end if
                                  end if
                              end do
    
                              ! Again calculate the multipole values to verify whether they are represented exactly
                              do l=0,lmax
                                  do m=-l,l
                                      if (l==0) then
                                          !monopole(impl) = monopole(impl) + ttl*hhh!*sqrt(1.d0/(4.d0*pi))
                                          !monopole(impl) = monopole(impl) + ttl*hhh*solid_harmonic()
                                          monopole(impl) = monopole(impl) + ttl*hhh*&
                                                           solid_harmonic(0,0.d0,l,m,r(1),r(2),r(3))*&
                                                           sqrt(4.d0*pi/real(2*l+1,kind=8))
                                      else if (l==1) then
                                          if (m==-1) then
                                              !dipole(1,impl) = dipole(1,impl) + ttl*hhh*r(2)*sqrt(3.d0/(4.d0*pi))
                                              ii = 1
                                          else if (m==0) then
                                              !dipole(2,impl) = dipole(2,impl) + ttl*hhh*r(3)*sqrt(3.d0/(4.d0*pi))
                                              ii = 2
                                          else if (m==1) then
                                              !dipole(3,impl) = dipole(3,impl) + ttl*hhh*r(1)*sqrt(3.d0/(4.d0*pi))
                                              ii = 3
                                          end if
                                          dipole(ii,impl) = dipole(ii,impl) + ttl*hhh*&
                                                           solid_harmonic(0,0.d0,l,m,r(1),r(2),r(3))*&
                                                           sqrt(4.d0*pi/real(2*l+1,kind=8))
                                      else if (l==2) then
                                          if (m==-2) then
                                              !quadrupole(1,impl) = quadrupole(1,impl) + ttl*hhh*r(1)*r(2)*sqrt(15.d0/(4.d0*pi))
                                              ii = 1
                                          else if (m==-1) then
                                              !quadrupole(2,impl) = quadrupole(2,impl) + ttl*hhh*r(2)*r(3)*sqrt(15.d0/(4.d0*pi))
                                              ii = 2
                                          else if (m==0) then
                                              !quadrupole(3,impl) = quadrupole(3,impl) + ttl*hhh*&
                                              !                     (-r(1)**2-r(2)**2+2.d0*r(3)**2)*sqrt(5.d0/(16.d0*pi))
                                              ii = 3
                                          else if (m==1) then
                                              !quadrupole(4,impl) = quadrupole(4,impl) + ttl*hhh*r(1)*r(3)*sqrt(15.d0/(4.d0*pi))
                                              ii = 4
                                          else if (m==2) then
                                              !quadrupole(5,impl) = quadrupole(5,impl) + ttl*hhh*&
                                              !                     (r(1)**2-r(2)**2)*sqrt(15.d0/(16.d0*pi))
                                              ii = 5
                                          end if
                                          quadrupole(ii,impl) = quadrupole(ii,impl) + ttl*hhh*&
                                                           solid_harmonic(0,0.d0,l,m,r(1),r(2),r(3))*&
                                                           sqrt(4.d0*pi/real(2*l+1,kind=8))
                                      end if
                                  end do
                              end do
                              !density_try(i1,i2,i3,ithread) = density_try(i1,i2,i3,ithread) + tt
                              ! If the norm of the Gaussian is close to one, 
                              density_loc(i1,i2,i3,ithread) = density_loc(i1,i2,i3,ithread) + tt
                          else
                              ! Use the method based on the analytic formula
                              rnrm1 = sqrt(rnrm2)
                              rnrm3 = rnrm1*rnrm2
                              rnrm5 = rnrm3*rnrm2
                              tt = 0.0_dp
                              do l=0,lmax
                                  if (associated(ep%mpl(impl)%qlm(l)%q)) then
                                      select case(l)
                                      case (0)
                                          tt = tt + calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                          !write(*,'(a,3es12.4,es16.8)') 'x, y, z, monopole', x, y, z, calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                      case (1)
                                          tt = tt + calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                          !write(*,*) 'dipole', calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                      case (2)
                                          tt = tt + calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                          !write(*,*) 'quadrupole', calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                      case (3)
                                          call f_err_throw('octupole not yet implemented', err_name='BIGDFT_RUNTIME_ERROR')
                                          !multipole_terms(l) = calc_octopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                      case default
                                          call f_err_throw('Wrong value of l', err_name='BIGDFT_RUNTIME_ERROR')
                                      end select
                                  end if
                              end do
                              potential_loc(i1,i2,i3,ithread) = potential_loc(i1,i2,i3,ithread) + tt
                          end if
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel
    
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
                      if (tt>1.d-2) then
                          !write(*,*) 'ERROR', abs(norm(l,impl)-norm_check(l,impl)), norm_check(l,impl)
                          call f_err_throw('The deviation from normalization of the radial function is too large: '//&
                              yaml_toa(tt,fmt='(es7.1)'), err_name='BIGDFT_RUNTIME_ERROR')
                      end if
                      if (tt>1.d-4) then
                          !write(*,*) 'ERROR', abs(norm(l,impl)-norm_check(l,impl)), norm_check(l,impl)
                          call yaml_warning('The deviation from normalization of the radial function is large: '//&
                              yaml_toa(tt,fmt='(es7.1)'))
                      end if
                  end do
              end if
          end do
    
          call f_free(rmax)
          call f_free(norm_check)
    
    
          if (iproc==0 .and. ep%nmpl > 0) then
                  !do iat=1,nat
                  !    call yaml_sequence(advance='no')
                  !    atomname=atomnames(iatype(iat))
                  !    call yaml_sequence_open(trim(atomname))
                  !    do l=0,lmax
                  !        call yaml_sequence(advance='no')
                  !        !call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat),fmt='(1es16.8)')
                  !        !call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat)*sqrt(4.d0**(2*l+3)),fmt='(1es16.8)')
                  !        !do m=-l,l
                  !            !multipoles(m,l,iat) = multipoles(m,l,iat)*get_normalization(rmax, l, m)
                  !            !max_error = max(max_error,abs(multipoles(m,l,iat)-get_test_factor(l,m)))
                  !        !end do
                  !        call yaml_map('l='//yaml_toa(l),multipoles_tmp(-l:l,l,iat),fmt='(1es16.8)')
                  !        call yaml_newline()
                  !    end do
                  !    !call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
                  !    call yaml_sequence_close()
                  !end do
                  !call yaml_sequence_close()
              call yaml_mapping_open('Potential from multipoles')
              call yaml_map('Number of multipole centers',ep%nmpl)
              call yaml_map('Sigma of the Gaussians',sigma)
              call yaml_map('Threshold for the norm of the Gaussians',norm_threshold)
              call yaml_map('Minimal radius for divion of the solid harmonics by r^{2l}',rmin)
              call yaml_sequence_open('Details for each multipole')
              do impl=1,ep%nmpl
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open(trim(yaml_toa(impl)))
                  if (norm_ok(impl)) then
                      call yaml_map('Method','Density based on Gaussians')
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
                                  !if (l==0) then
                                  !    max_error(l) = max(max_error(l), &
                                  !                    abs(monopole(impl)-(ep%mpl(impl)%qlm(l)%q(mm)+real(nelpsp(impl),kind=8))))
                                  !else if (l==1) then
                                  !    max_error(l) = max(max_error(l),abs(dipole(mm,impl)-ep%mpl(impl)%qlm(l)%q(mm)))
                                  !else if (l==2) then
                                  !    max_error(l) = max(max_error(l),abs(quadrupole(mm,impl)-ep%mpl(impl)%qlm(l)%q(mm)))
                                  !end if
                                  if (l==0) then
                                      qq = ep%mpl(impl)%qlm(l)%q(mm)+real(nelpsp(impl),kind=8)
                                  else
                                      qq = ep%mpl(impl)%qlm(l)%q(mm)
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
                                  !!if (l==0) then
                                  !!    max_error(l) = max(max_error(l), &
                                  !!                    monopole(impl)/(ep%mpl(impl)%qlm(l)%q(mm)+real(nelpsp(impl),kind=8)))
                                  !!else if (l==1) then
                                  !!    max_error(l) = max(max_error(l),dipole(mm,impl)/ep%mpl(impl)%qlm(l)%q(mm))
                                  !!    !write(*,*) 'calc, orig', dipole(mm,impl), ep%mpl(impl)%qlm(l)%q(mm)
                                  !!else if (l==2) then
                                  !!    max_error(l) = max(max_error(l),quadrupole(mm,impl)/ep%mpl(impl)%qlm(l)%q(mm))
                                  !!    !write(*,*) 'calc, orig', quadrupole(mm,impl),ep%mpl(impl)%qlm(l)%q(mm)
                                  !!end if
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
                      !call yaml_map('Maximal deviation from the original values in percent',max_error(:),fmt='(1f6.1)')
                      call yaml_map('Maximal deviation from the original values in percent',output_arr(:))
                  else
                      call yaml_map('Method','Analytic expression')
                  end if
                  !call yaml_map('monopole',monopole(impl),fmt='(1es16.8)')
                  !call yaml_map('dipole',dipole(:,impl),fmt='(1es16.8)')
                  !call yaml_map('quadrupole',quadrupole(:,impl),fmt='(1es16.8)')
                  call yaml_mapping_close()
              end do
              call yaml_sequence_close()
              call yaml_mapping_close()
          end if
    
          if (ep%nmpl > 0) then
             call H_potential('D',denspot%pkernel,density,denspot%V_ext,ehart_ps,0.0_dp,.false.,&
                  quiet=denspot%PSquiet)!,rho_ion=denspot%rho_ion)
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
    
    
!!  UNCOMMENT FOR TEST          ii = 0
!!  UNCOMMENT FOR TEST          do i3=is3,ie3
!!  UNCOMMENT FOR TEST              ii3 = i3 - 15
!!  UNCOMMENT FOR TEST              do i2=is2,ie2
!!  UNCOMMENT FOR TEST                  ii2 = i2 - 15
!!  UNCOMMENT FOR TEST                  do i1=is1,ie1
!!  UNCOMMENT FOR TEST                      ii1 = i1 - 15
!!  UNCOMMENT FOR TEST                      ii = ii + 1
!!  UNCOMMENT FOR TEST                      write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
!!  UNCOMMENT FOR TEST                      !do impl=1,ep%nmpl
!!  UNCOMMENT FOR TEST                      !    r(1) = ep%mpl(impl)%rxyz(1) - x
!!  UNCOMMENT FOR TEST                      !    r(2) = ep%mpl(impl)%rxyz(2) - y
!!  UNCOMMENT FOR TEST                      !    r(3) = ep%mpl(impl)%rxyz(3) - z 
!!  UNCOMMENT FOR TEST                      !    rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
!!  UNCOMMENT FOR TEST                      !    rnrm1 = sqrt(rnrm2)
!!  UNCOMMENT FOR TEST                      !    tt = spherical_harmonic(l, m, x, y, z)*gaussian(sigma, rnrm1)
!!  UNCOMMENT FOR TEST                      !    density(i1,i2,i3) =+ tt
!!  UNCOMMENT FOR TEST                      !    !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
!!  UNCOMMENT FOR TEST                      !end do
!!  UNCOMMENT FOR TEST                  end do
!!  UNCOMMENT FOR TEST              end do
!!  UNCOMMENT FOR TEST          end do
    
          call f_free(density)
          call f_free(nzatom)
          call f_free(nelpsp)
          call f_free(npspcode)
          call f_free(psppar)

      end if multipoles_if

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
              g = g/sqrt(2.d0*pi*sigma**2)**3
          else
              g = 0.d0
          end if
          !g = g/(sigma**3*sqrt(2.d0*pi)**3)


        end function gaussian

        !!!> Calculates the real spherical harmonic for given values of l, m, x, y, z.
        !!function spherical_harmonic(l, m, x, y, z) result(sh)
        !!  use module_base, only: pi => pi_param
        !!  implicit none
        !!  ! Calling arguments
        !!  integer,intent(in) :: l, m
        !!  real(kind=8),intent(in) :: x, y, z
        !!  real(kind=8) :: sh

        !!  ! Local variables
        !!  integer,parameter :: l_max=2
        !!  real(kind=8) :: r, r2, rnorm
        !!  real(kind=8),parameter :: sqrt_1_over_pi = sqrt(1/pi)
        !!  real(kind=8),parameter :: sqrt_3_over_4pi = sqrt(3.d0/(4.d0*pi))
        !!  real(kind=8),parameter :: sqrt_15_over_pi = sqrt(15.d0/pi)
        !!  real(kind=8),parameter :: sqrt_5_over_pi = sqrt(5.d0/pi)


        !!  if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
        !!  if (l>l_max) call f_err_throw('spherical harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
        !!      err_name='BIGDFT_RUNTIME_ERROR')
        !!  if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


        !!  ! Normalization for a sphere of radius rmax
        !!  select case (l)
        !!  case (0)
        !!      sh = 0.5d0*sqrt_1_over_pi !0.5d0*sqrt(1/pi)
        !!  case (1)
        !!      r = sqrt(x**2+y**2+z**2)
        !!      ! fix for small r (needs proper handling later...)
        !!      if (r==0.d0) r=1.d-20
        !!      select case (m)
        !!      case (-1)
        !!          sh = sqrt_3_over_4pi*y/r !sqrt(3.d0/(4.d0*pi))*y/r
        !!      case (0)
        !!          sh = sqrt_3_over_4pi*z/r !sqrt(3.d0/(4.d0*pi))*z/r
        !!      case (1)
        !!          sh = sqrt_3_over_4pi*x/r !sqrt(3.d0/(4.d0*pi))*x/r
        !!      end select
        !!  case (2)
        !!      r2 = x**2+y**2+z**2
        !!      ! fix for small r2 (needs proper handling later...)
        !!      if (r2==0.d0) r2=1.d-20
        !!      select case (m)
        !!      case (-2)
        !!          sh = 0.5d0*sqrt_15_over_pi*x*y/r2 !0.5d0*sqrt(15.d0/pi)*x*y/r2
        !!      case (-1)
        !!          sh = 0.5d0*sqrt_15_over_pi*y*z/r2 !0.5d0*sqrt(15.d0/pi)*y*z/r2
        !!      case (0)
        !!          sh = 0.25d0*sqrt_5_over_pi*(-x**2-y**2+2*z**2)/r2 !0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
        !!      case (1)
        !!          sh = 0.5d0*sqrt_15_over_pi*z*x/r2 !0.5d0*sqrt(15.d0/pi)*z*x/r2
        !!      case (2)
        !!          sh = 0.25d0*sqrt_15_over_pi*(x**2-y**2)/r2 !0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
        !!      end select
        !!  end select


        !!end function spherical_harmonic

    end subroutine potential_from_charge_multipoles


    !> Calculate the external potential arising from the multipoles provided
    subroutine potential_from_multipoles(ep, is1, ie1, is2, ie2, is3, ie3, iis3, iie3, hx, hy, hz, shift, pot)
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3 !< parallelized box bounds
      integer,intent(in) :: iis3, iie3 !< non-parallelized box bounds (z direction)
      real(gp),intent(in) :: hx, hy, hz
      real(gp),dimension(3),intent(in) :: shift !< global shift of the atomic positions
      real(gp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp
      real(dp),dimension(3) :: r
      real(kind=8),parameter :: buffer = 1.d0
      real(kind=8) :: xxs, xxe, yys, yye, zzs, zze

      !stop 'deprecated'

      xxs = real(is1-15,kind=8)*hx + shift(1) - buffer
      xxe = real(ie1-15,kind=8)*hx + shift(1) + buffer
      yys = real(is2-15,kind=8)*hy + shift(2) - buffer
      yye = real(ie2-15,kind=8)*hy + shift(2) + buffer
      zzs = real(iis3-15,kind=8)*hz + shift(3) - buffer
      zze = real(iie3-15,kind=8)*hz + shift(3) + buffer

      !!$omp parallel &
      !!$omp default(none) &
      !!$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, ep, pot) &
      !!$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, rnrm1, rnrm2, rnrm3, rnrm5, l, mp)
      !!$omp do
      do impl=1,ep%nmpl
          ! Only take this atom into account if it lies outside of the simulation box (plus some buffer).
          ! Otherwise it has already been taken into account by the potential generated via the Poisson Solver
          if ( ep%mpl(impl)%rxyz(1) <= xxs .or. ep%mpl(impl)%rxyz(1) >= xxe .or. &
               ep%mpl(impl)%rxyz(2) <= yys .or. ep%mpl(impl)%rxyz(2) >= yye .or. &
               ep%mpl(impl)%rxyz(3) <= zzs .or. ep%mpl(impl)%rxyz(3) >= zze ) then
              !write(*,*) ep%mpl(impl)%rxyz(1) <= xxs
              !write(*,*) ep%mpl(impl)%rxyz(1) >= xxe
              !write(*,*) ep%mpl(impl)%rxyz(2) <= yys
              !write(*,*) ep%mpl(impl)%rxyz(2) >= yye
              !write(*,*) ep%mpl(impl)%rxyz(3) <= zzs
              !write(*,*) ep%mpl(impl)%rxyz(3) >= zze
              !write(*,*) 'ok for impl',impl, ep%mpl(impl)%rxyz(1:3), xxs, xxe, yys, yye, zzs, zze
              do i3=is3,ie3
                  ii3 = i3 - 15
                  z = real(ii3,kind=8)*hz + shift(3)
                  !write(*,'(a,i7,2es16.7)') 'i3, z, ep%mpl(1)%rxyz(3)', i3, z, ep%mpl(1)%rxyz(3)
                  do i2=is2,ie2
                      ii2 = i2 - 15
                      y = real(ii2,kind=8)*hy + shift(2)
                      do i1=is1,ie1
                          ii1 = i1 - 15
                          x = real(ii1,kind=8)*hx + shift(1)
                          r(1) = ep%mpl(impl)%rxyz(1) - x
                          r(2) = ep%mpl(impl)%rxyz(2) - y
                          r(3) = ep%mpl(impl)%rxyz(3) - z 
                          rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                          ! To avoid floating point exception
                          rnrm2=max(rnrm2,1.d-6)
                          rnrm1 = sqrt(rnrm2)
                          rnrm3 = rnrm1*rnrm2
                          rnrm5 = rnrm3*rnrm2
                          mp = 0.0_dp
                          do l=0,lmax
                              if (associated(ep%mpl(impl)%qlm(l)%q)) then
                                  select case(l)
                                  case (0)
                                      mp = mp + calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                      !write(*,'(a,3es12.4,es16.8)') 'x, y, z, monopole', x, y, z, calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                  case (1)
                                      mp = mp + calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                      !write(*,*) 'dipole', calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                  case (2)
                                      mp = mp + calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                      !write(*,*) 'quadrupole', calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                  case (3)
                                      call f_err_throw('octupole not yet implemented', err_name='BIGDFT_RUNTIME_ERROR')
                                      !multipole_terms(l) = calc_octopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                  case default
                                      call f_err_throw('Wrong value of l', err_name='BIGDFT_RUNTIME_ERROR')
                                  end select
                              end if
                          end do
                          pot(i1,i2,i3) = pot(i1,i2,i3) + mp
                          !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
                      end do
                  end do
              end do
              !!$omp end do
              !!$omp end parallel
          end if
      end do

    end subroutine potential_from_multipoles


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



!!!!    !> Calculates the multipole moments for each atom 
!!!!    subroutine multipoles_from_density(iproc, nproc, at, lzd, smats, smatl, orbs, &
!!!!               npsidim, lphi, norbsPerType, collcom, collcom_sr, orthpar, &
!!!!               ovrlp, kernel, meth_overlap, multipoles_out)
!!!!      use module_base
!!!!      use module_types
!!!!      use sparsematrix_base, only: sparsematrix_malloc0, SPARSE_FULL, assignment(=)
!!!!      use sparsematrix_init, only: matrixindex_in_compressed
!!!!      use orthonormalization, only: orthonormalizeLocalized
!!!!      use yaml_output
!!!!      use locreg_operations
!!!!      use bounds, only: geocode_buffers
!!!!      implicit none
!!!!
!!!!      ! Calling arguments
!!!!      integer,intent(in) :: iproc, nproc, npsidim, meth_overlap
!!!!      type(atoms_data),intent(in) :: at
!!!!      type(local_zone_descriptors),intent(in) :: lzd
!!!!      type(sparse_matrix),intent(inout) :: smats, smatl
!!!!      type(orbitals_data),intent(in) :: orbs
!!!!      real(kind=8),dimension(npsidim),intent(in) :: lphi
!!!!      integer,dimension(at%astruct%ntypes),intent(in) :: norbsPerType
!!!!      type(comms_linear),intent(in) :: collcom
!!!!      type(comms_linear),intent(in) :: collcom_sr
!!!!      type(orthon_data),intent(in) :: orthpar
!!!!      type(matrices),intent(in) :: kernel
!!!!      type(matrices),intent(inout) :: ovrlp !in principle also intent(in)
!!!!      real(kind=8),dimension(-lmax:lmax,0:lmax,1:at%astruct%nat),intent(out),optional :: multipoles_out
!!!!
!!!!      ! Local variables
!!!!      integer :: ist, istr, iorb, iiorb, ilr, ii, natp, isat, nr, jproc, iat, n, norb_get, istr_get
!!!!      integer :: window, ioffset, methTransformOverlap, l, m, iiat, ityp, norb_per_atom, i1, i2, i3, ind, jorb, jat
!!!!      integer :: ii1, ii2, ii3, jjorb, i, itype
!!!!      real(kind=8),dimension(:),allocatable :: psir
!!!!      real(kind=8),dimension(:),pointer :: phit_c, phit_f
!!!!      type(workarr_sumrho) :: w
!!!!      integer,dimension(:),allocatable :: nat_par, norb_list
!!!!      real(kind=8),dimension(:),allocatable :: psir_get, locrad, rmax
!!!!      real(kind=8),dimension(:,:),allocatable :: locregcenter, psir_get_fake
!!!!      real(kind=8),dimension(:,:,:,:),allocatable :: phi
!!!!      real(kind=8),dimension(:,:,:,:,:,:),allocatable :: sphi
!!!!      integer,dimension(:,:),allocatable :: comms
!!!!      logical :: can_use_transposed, arr_allocated
!!!!      real(kind=8) :: ddot, x, y, z, tt, rnorm, factor, max_error, q!, get_normalization, get_test_factor
!!!!      !real(kind=8) ,dimension(2,orbs%norb) :: testarr
!!!!      real(kind=8),dimension(:),allocatable :: kernel_ortho, phi_ortho
!!!!      real(kind=8),dimension(:,:),allocatable :: weight_centers
!!!!      integer,dimension(:),allocatable :: n1i, n2i, n3i, ns1i, ns2i, ns3i
!!!!      !real(kind=8),dimension(-lmax:lmax,0:lmax,at%astruct%nat) :: multipoles
!!!!      real(kind=8),dimension(:,:,:),allocatable :: multipoles
!!!!      real(kind=8) :: factor_normalization, hxh, hyh, hzh, weight
!!!!      character(len=20) :: atomname
!!!!      real(kind=8),dimension(-lmax:lmax,0:lmax) :: norm
!!!!      integer :: nl1, nl2, nl3
!!!!      !real(kind=8),parameter :: rmax=5.d0
!!!!      !testarr = 0.d0
!!!!
!!!!      call f_routine(id='multipoles_from_density')
!!!!
!!!!      if (iproc==0) call yaml_comment('Atomic multipole analysis, old approach',hfill='=')
!!!!
!!!!      multipoles = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.at%astruct%nat/),id='multipoles')
!!!!
!!!!
!!!!      call unitary_test()
!!!!
!!!!      if (iproc==0) then
!!!!          call yaml_mapping_open('Multipole analysis')
!!!!      end if
!!!!
!!!!      ! Orthogonalize the support functions
!!!!      can_use_transposed = .false.
!!!!      methTransformOverlap = 1020
!!!!      phit_c = f_malloc_ptr(collcom%ndimind_c,id='phit_c')
!!!!      phit_f = f_malloc_ptr(7*collcom%ndimind_f,id='phit_f')
!!!!      phi_ortho = f_malloc(npsidim,id='phi_ortho')
!!!!      call vcopy(npsidim, lphi(1), 1, phi_ortho(1), 1)
!!!!      if (iproc==0) then
!!!!          call yaml_map('Orthonormalizing support functions',.true.)
!!!!      end if
!!!!      call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
!!!!           1.d-8, npsidim, orbs, lzd, &
!!!!           smats, smatl, collcom, orthpar, &
!!!!           phi_ortho, phit_c, phit_f, &
!!!!           can_use_transposed)
!!!!      call f_free_ptr(phit_c)
!!!!      call f_free_ptr(phit_f)
!!!!
!!!!
!!!!      ! Transform the support functions to real space
!!!!      psir = f_malloc(max(collcom_sr%ndimpsi_c,1),id='psir')
!!!!      ist=1
!!!!      istr=1
!!!!      do iorb=1,orbs%norbp
!!!!          iiorb=orbs%isorb+iorb
!!!!          ilr=orbs%inwhichlocreg(iiorb)
!!!!          call initialize_work_arrays_sumrho(1,[lzd%Llr(ilr)],.true.,w)
!!!!          call daub_to_isf(lzd%Llr(ilr), w, phi_ortho(ist), psir(istr))
!!!!          call deallocate_work_arrays_sumrho(w)
!!!!          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
!!!!          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!!!      end do
!!!!      do i=1,ist-1
!!!!          write(445,*) i, phi_ortho(i)
!!!!      end do
!!!!      do i=1,istr-1
!!!!          write(444,*) i, psir(i)
!!!!      end do
!!!!      call f_free(phi_ortho)
!!!!      if(istr/=collcom_sr%ndimpsi_c+1) then
!!!!          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
!!!!          stop
!!!!      end if
!!!!
!!!!
!!!!      !! NEW: CALCULATE THE WEIGHT CENTER OF EACH SUPPORT FUNCTION ############################
!!!!      !weight_centers = f_malloc0((/3,orbs%norb/),id='weight_centers')
!!!!      !hxh = 0.5d0*lzd%hgrids(1)
!!!!      !hyh = 0.5d0*lzd%hgrids(2)
!!!!      !hzh = 0.5d0*lzd%hgrids(3)
!!!!
!!!!      !istr = 1
!!!!      !do iorb=1,orbs%norbp
!!!!      !    iiorb=orbs%isorb+iorb
!!!!      !    ilr=orbs%inwhichlocreg(iiorb)
!!!!      !    call geocode_buffers(lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)
!!!!      !    !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
!!!!      !    !com(1:3,iorb) = 0.d0
!!!!      !    weight = 0.d0
!!!!      !    do i3=1,lzd%llr(ilr)%d%n3i
!!!!      !        ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
!!!!      !        z = ii3*hzh
!!!!      !        do i2=1,lzd%llr(ilr)%d%n2i
!!!!      !            ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
!!!!      !            y = ii2*hyh
!!!!      !            do i1=1,lzd%llr(ilr)%d%n1i
!!!!      !                ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
!!!!      !                x = ii1*hxh
!!!!      !                tt = psir(istr)**2
!!!!      !                weight_centers(1,iiorb) = weight_centers(1,iiorb) + x*tt
!!!!      !                weight_centers(2,iiorb) = weight_centers(2,iiorb) + y*tt
!!!!      !                weight_centers(3,iiorb) = weight_centers(3,iiorb) + z*tt
!!!!      !                weight = weight + tt
!!!!      !                istr = istr + 1
!!!!      !            end do
!!!!      !        end do
!!!!      !    end do
!!!!      !    !call yaml_map('weight',weight)
!!!!      !    weight_centers(1:3,iorb) = weight_centers(1:3,iorb)/weight
!!!!      !end do
!!!!      !call mpiallred(weight_centers, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!!!      !if (iproc==0) then
!!!!      !    do iorb=1,orbs%norb
!!!!      !        write(*,'(a,i4,3es13.4)') 'iorb, weight_centers(1:3,iorb)', iorb, weight_centers(1:3,iorb)
!!!!      !    end do
!!!!      !end if
!!!!      !call f_free(weight_centers)
!!!!      !! ######################################################################################
!!!!
!!!!
!!!!      ! Switch to a partition over the atoms
!!!!      nat_par = f_malloc(0.to.nproc-1,id='nat_par')
!!!!      ii = at%astruct%nat/nproc
!!!!      nat_par(0:nproc-1) = ii
!!!!      ii = at%astruct%nat-nproc*ii
!!!!      nat_par(0:ii-1) = nat_par(0:ii-1) + 1
!!!!      if (sum(nat_par(:))/=at%astruct%nat) then
!!!!          call f_err_throw('wrong partition of the atoms',err_name='BIGDFT_RUNTIME_ERROR')
!!!!      end if
!!!!      natp = nat_par(iproc)
!!!!      if (iproc==0) then
!!!!          isat = 0
!!!!      else
!!!!          isat = sum(nat_par(0:iproc-1))
!!!!      end if
!!!!      call f_free(nat_par)
!!!!      comms = f_malloc((/4,orbs%norb/),id='comms')
!!!!      !write(*,'(a,i5,3x,13i4)') 'iproc, orbs%onwhichatom', iproc, orbs%onwhichatom
!!!!      nr = 0
!!!!      norb_get = 0
!!!!      istr = 1
!!!!      istr_get = 1
!!!!      do iat=isat+1,isat+natp
!!!!          do jproc=0,nproc-1
!!!!              istr = 0
!!!!              do iorb=1,orbs%norb_par(jproc,0)
!!!!                  iiorb = iorb + orbs%isorb_par(jproc)
!!!!                  ilr = orbs%inwhichlocreg(iiorb)
!!!!                  iiat = orbs%onwhichatom(iiorb)
!!!!                  n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!!!                  !write(*,'(a,5i8)') 'iproc, jproc, iiorb, iat', iproc, jproc, iiorb, iat
!!!!                  !if (iat>=isat+1 .and. iat<=isat+natp) then
!!!!                  if (iiat==iat) then
!!!!                      norb_get = norb_get + 1
!!!!                      comms(1,norb_get) = jproc
!!!!                      comms(2,norb_get) = n
!!!!                      comms(3,norb_get) = istr
!!!!                      comms(4,norb_get) = istr_get
!!!!                      nr = nr + n
!!!!                      !write(*,'(a,5i8)') 'iproc, jproc, n, istr, istr_get', iproc, jproc, n, istr, istr_get
!!!!                      istr = istr + n
!!!!                      istr_get = istr_get + n
!!!!                  else
!!!!                      istr = istr + n
!!!!                  end if
!!!!              end do
!!!!          end do
!!!!      end do
!!!!      !write(*,*) 'iproc, nr', iproc, nr
!!!!      !do iorb=1,norb_get
!!!!      !    write(*,'(a,2i7,4i9)') 'iproc, iorb, comm',iproc, iorb, comms(:,iorb)
!!!!      !end do
!!!!      psir_get = f_malloc(nr,id='psir_get')
!!!!      if (nproc>1) then
!!!!          window = mpiwindow(size(psir), psir(1), bigdft_mpi%mpi_comm)
!!!!          do iorb=1,norb_get
!!!!              jproc = comms(1,iorb)
!!!!              n = comms(2,iorb)
!!!!              ioffset = comms(3,iorb)
!!!!              istr = comms(4,iorb)
!!!!              !write(*,'(5(a,i0))') 'task ',iproc,' gets ',n,' elements at position ', &
!!!!              !                     istr,' from position ',ioffset+1,' on task ',jproc
!!!!              call mpiget(psir_get(istr), n, jproc, int(ioffset,kind=mpi_address_kind), window)
!!!!          end do
!!!!          call mpi_fenceandfree(window)
!!!!      else
!!!!          do iorb=1,norb_get
!!!!              n = comms(2,iorb)
!!!!              ioffset = comms(3,iorb)
!!!!              istr = comms(4,iorb)
!!!!              call vcopy(n, psir(ioffset+1), 1, psir_get(istr), 1)
!!!!          end do
!!!!      end if
!!!!      call f_free(psir)
!!!!      call f_free(comms)
!!!!      istr = 1
!!!!      !write(*,*) 'iproc, isat, natp', iproc, isat, natp
!!!!      do iorb=1,orbs%norb
!!!!          ilr = orbs%inwhichlocreg(iorb)
!!!!          iat = orbs%onwhichatom(iorb)
!!!!          if (iat>=isat+1 .and. iat<=isat+natp) then
!!!!              n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!!!              !write(*,'(a,4i8,es16.6)') 'iproc, iorb, n, istr, ddot', &
!!!!              !    iproc, iorb, n, istr, ddot(n, psir_get(istr), 1, psir_get(istr), 1)
!!!!              !testarr(2,iorb) = ddot(n, psir_get(istr), 1, psir_get(istr), 1)
!!!!              istr = istr + n
!!!!          end if
!!!!      end do
!!!!      !do iorb=1,size(psir_get)
!!!!      !    write(300+iproc,'(a,i7,es16.7)') 'i, psir_get(i)', iorb, psir_get(iorb)
!!!!      !end do
!!!!      !call mpiallred(testarr, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!!!      !if (iproc==0) then
!!!!      !    do iorb=1,orbs%norb
!!!!      !        write(*,*) 'DIFF, iorb, val', iorb, abs(testarr(1,iorb)-testarr(2,iorb))
!!!!      !    end do
!!!!      !end if
!!!!
!!!!
!!!!      ! Calculate the kernel for orthormal support functions
!!!!      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='kernel_ortho')
!!!!      call kernel_for_orthonormal_basis(iproc, nproc, orbs%norbp, meth_overlap, smats, smatl, &
!!!!           ovrlp, kernel, kernel_ortho)
!!!!
!!!!
!!!!      ! Apply the spherical harmonics to the suport functions.
!!!!      ! For the support functions on atom A we only need to apply 
!!!!      ! the spherical harmonics centered as well on atom A.
!!!!
!!!!
!!!!
!!!!      !lmax = 1
!!!!      n1i = f_malloc(orbs%norb,id='n1i')
!!!!      n2i = f_malloc(orbs%norb,id='n2i')
!!!!      n3i = f_malloc(orbs%norb,id='n3i')
!!!!      ns1i = f_malloc(orbs%norb,id='ns1i')
!!!!      ns2i = f_malloc(orbs%norb,id='ns2i')
!!!!      ns3i = f_malloc(orbs%norb,id='ns3i')
!!!!      locrad = f_malloc(orbs%norb,id='locrad')
!!!!      locregcenter = f_malloc((/3,orbs%norb/),id='locregcenter')
!!!!      do iorb=1,orbs%norb
!!!!          n1i(iorb) = lzd%llr(iorb)%d%n1i
!!!!          n2i(iorb) = lzd%llr(iorb)%d%n2i
!!!!          n3i(iorb) = lzd%llr(iorb)%d%n3i
!!!!          ns1i(iorb) = lzd%llr(iorb)%nsi1
!!!!          ns2i(iorb) = lzd%llr(iorb)%nsi2
!!!!          ns3i(iorb) = lzd%llr(iorb)%nsi3
!!!!          locrad(iorb) = lzd%llr(iorb)%locrad
!!!!          locregcenter(1:3,iorb) = lzd%llr(iorb)%locregcenter(1:3)
!!!!      end do
!!!!      rmax = f_malloc(at%astruct%nat,id='rmax')
!!!!      call multipole_analysis_core(iproc, nproc, natp, isat, at%astruct%nat, at%astruct%ntypes, orbs%norb, &
!!!!           0, at%astruct%iatype, norbsPerType, orbs%inwhichlocreg, orbs%onwhichatom, &
!!!!           n1i, n2i, n3i, ns1i, ns2i, ns3i, locrad, lzd%hgrids, locregcenter, &
!!!!           nr, psir_get, psir_get, &
!!!!           smatl%nvctr, kernel_ortho, 1, multipoles, rmax, 101, smatl)!, matrixindex)
!!!!      call f_free(psir_get_fake)
!!!!      call f_free(n1i)
!!!!      call f_free(n2i)
!!!!      call f_free(n3i)
!!!!      call f_free(ns1i)
!!!!      call f_free(ns2i)
!!!!      call f_free(ns3i)
!!!!      call f_free(locrad)
!!!!      call f_free(locregcenter)
!!!!
!!!!      ! The monopole term should be the net charge, i.e. add the positive atomic charges
!!!!      do iat=1,at%astruct%nat
!!!!          itype = at%astruct%iatype(iat)
!!!!          q = real(at%nelpsp(itype),kind=8)
!!!!          multipoles(0,0,iat) = multipoles(0,0,iat) + q
!!!!      end do
!!!!
!!!!      if (iproc==0) then
!!!!          !!call write_multipoles(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%atomnames, &
!!!!          !!     multipoles, rmax, lzd%hgrids, without_normalization=.false.)
!!!!          call write_multipoles_new(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%atomnames, &
!!!!               at%astruct%rxyz, at%astruct%units, multipoles)
!!!!      end if
!!!!      call f_free(rmax)
!!!!      call f_free(psir_get)
!!!!
!!!!      call f_free(kernel_ortho)
!!!!
!!!!      if (iproc==0) then
!!!!          call yaml_mapping_close()
!!!!      end if
!!!!
!!!!
!!!!      if (present(multipoles_out)) then
!!!!          call f_memcpy(src=multipoles,dest=multipoles_out)
!!!!      end if
!!!!      call f_free(multipoles)
!!!!
!!!!      if (iproc==0) then
!!!!          call yaml_comment('Atomic multipole analysis done',hfill='=')
!!!!      end if
!!!!
!!!!      call f_release_routine()
!!!!
!!!!
!!!!
!!!!      contains
!!!!
!!!!
!!!!        subroutine unitary_test()
!!!!          implicit none
!!!!          real(kind=8),dimension(1) :: rmax
!!!!          integer,parameter :: n1i=101, n2i=81, n3i=91
!!!!          integer,parameter :: nsi1=0, nsi2=10, nsi3=20
!!!!          !integer,parameter :: n1i=100, n2i=100, n3i=100
!!!!          !integer,parameter :: nsi1=0, nsi2=0, nsi3=0
!!!!          real(kind=8),dimension(3) :: locregcenter
!!!!          integer :: nr
!!!!          real(kind=8) :: factor_normalization, r2, r
!!!!
!!!!          if (iproc==0) then
!!!!              call yaml_mapping_open('Unitary test for multipoles')
!!!!          end if
!!!!
!!!!          locregcenter(1) = (ceiling(real(n1i,kind=8)/2.d0)+nsi1-14-1)*0.5d0*lzd%hgrids(1)
!!!!          locregcenter(2) = (ceiling(real(n2i,kind=8)/2.d0)+nsi2-14-1)*0.5d0*lzd%hgrids(2)
!!!!          locregcenter(3) = (ceiling(real(n3i,kind=8)/2.d0)+nsi3-14-1)*0.5d0*lzd%hgrids(3)
!!!!
!!!!          !psir_get_fake = f_malloc0((/lzd%llr(1)%d%n1i*lzd%llr(1)%d%n2i*lzd%llr(1)%d%n3i,2/),id='psir_get_fake')
!!!!          nr = n1i*n2i*n3i
!!!!          psir_get_fake = f_malloc0((/nr,2/),id='psir_get_fake')
!!!!          psir_get_fake(:,2) = 1.d0
!!!!          rmax(1) = min(n1i*0.25d0*lzd%hgrids(1), &
!!!!                        n2i*0.25d0*lzd%hgrids(2), &
!!!!                        n3i*0.25d0*lzd%hgrids(3))
!!!!          ! Normalized to within a sphere of radius rmax
!!!!          factor_normalization = 3.d0/(4.d0*pi*rmax(1)**3)*0.5d0*lzd%hgrids(1)*0.5d0*lzd%hgrids(2)*0.5d0*lzd%hgrids(3)
!!!!          do i3=1,n3i!lzd%llr(1)%d%n3i
!!!!              ii3 = nsi3 + i3 - 14 - 1
!!!!              z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3)
!!!!              do i2=1,n2i!lzd%llr(1)%d%n2i
!!!!                  ii2 = nsi2 + i2 - 14 - 1
!!!!                  y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2)
!!!!                  do i1=1,n1i!lzd%llr(1)%d%n1i
!!!!                      ii1 = nsi1 + i1 - 14 - 1
!!!!                      x = ii1*0.5d0*lzd%hgrids(1) - locregcenter(1)
!!!!                      r2 = x**2+y**2+z**2
!!!!                      if (r2>rmax(1)**2) cycle
!!!!                      r = sqrt(r2)
!!!!                      r = max(0.5d0,r)
!!!!                      ind = (i3-1)*n2i*n1i + (i2-1)*n1i + i1
!!!!                      do l=0,lmax
!!!!                          do m=-l,l
!!!!                              factor = get_test_factor(l,m)*factor_normalization*sqrt(4.d0*pi*real(2*l+1,kind=8))
!!!!                              ! The minus sign is necessary to compensate the fact that the multipole calculation uses itself a minus sign (since it
!!!!                              ! assumes that the charge density (which is a negative quantity) is positive)
!!!!                              psir_get_fake(ind,1) = psir_get_fake(ind,1) - factor*solid_harmonic(-2, 1.0d0, l, m , x, y, z)
!!!!                          end do
!!!!                      end do
!!!!                  end do
!!!!              end do
!!!!          end do
!!!!          ! Only do it for one MPI task
!!!!          if (iproc==0) then
!!!!              call multipole_analysis_core(0, 1, 1, 0, 1, 1, 1, &
!!!!                   0, (/1/), (/1/), (/1/), (/1/), &
!!!!                   (/n1i/), (/n2i/), (/n3i/), &
!!!!                   (/nsi1/), (/nsi2/), (/nsi3/), rmax, &
!!!!                   lzd%hgrids, locregcenter, &
!!!!                   nr, psir_get_fake(:,1), psir_get_fake(:,2), &
!!!!                   1, (/1.d0/), 1, multipoles, rmax, 102, matrixindex=(/1/))
!!!!          end if
!!!!
!!!!          if (iproc==0) then
!!!!              call write_multipoles_new(1, 1, (/1/), (/'testatom'/), (/0.d0,0.d0,0.d0/), 'fake', &
!!!!                   multipoles, check_values_=.true.)
!!!!              call yaml_mapping_close()
!!!!          end if
!!!!
!!!!        end subroutine unitary_test
!!!!
!!!!
!!!!    end subroutine multipoles_from_density



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
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(matrices),intent(in) :: kernel
      type(matrices),intent(inout) :: ovrlp
      real(kind=8),dimension(smatl%nvctr*smatl%nspin),intent(out) :: weight_matrix_compr

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
      weight_matrix_compr_tg = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='weight_matrix_compr_tg')
      !if (norbp>0) then
         call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
              inv_ovrlp(1)%matrix_compr, proj_ovrlp_half_compr, weight_matrix_compr_tg)
      !end if
      call f_free(proj_ovrlp_half_compr)

      call deallocate_matrices(inv_ovrlp(1))

      ! Maybe this can be improved... not really necessary to gather the entire matrix
      !weight_matrix_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='weight_matrix_compr')
      call gather_matrix_from_taskgroups(iproc, nproc, smatl, weight_matrix_compr_tg, weight_matrix_compr)

      call f_free(weight_matrix_compr_tg)

      if (iproc==0) then
          call yaml_comment('Kernel calculated',hfill='~')
      end if

      call f_release_routine()

    end subroutine kernel_for_orthonormal_basis


    !!!subroutine write_multipoles(nat, ntypes, iatype, atomnames, multipoles, rmax, hgrids, without_normalization, check_values_)
    !!!  use yaml_output
    !!!  implicit none
    !!!  
    !!!  ! Calling arguments
    !!!  integer,intent(in) :: nat, ntypes
    !!!  integer,dimension(nat),intent(in) :: iatype
    !!!  character(len=*),dimension(ntypes),intent(in) :: atomnames
    !!!  real(kind=8),dimension(-lmax:lmax,0:lmax,nat),intent(in) :: multipoles
    !!!  real(kind=8),dimension(nat),intent(in) :: rmax
    !!!  real(kind=8),dimension(3),intent(in) :: hgrids
    !!!  logical,intent(in) :: without_normalization
    !!!  logical,intent(in),optional :: check_values_
    !!!  
    !!!  ! Local variables
    !!!  character(len=20) :: atomname
    !!!  integer :: i, iat, l, m, nit
    !!!  real(kind=8) :: max_error, factor!, get_normalization, get_test_factor
    !!!  real(kind=8),dimension(:,:,:),allocatable :: multipoles_tmp
    !!!  logical :: check_values

    !!!  if (present(check_values_)) then
    !!!      check_values = check_values_
    !!!  else
    !!!      check_values = .false.
    !!!  end if

    !!!      multipoles_tmp = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.nat/),id='multipoles_tmp')

    !!!      if (without_normalization) then
    !!!          nit = 2
    !!!      else
    !!!          nit = 1
    !!!      end if

    !!!      factor = 0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3)

    !!!      max_error = 0.d0
    !!!      call yaml_mapping_open('Multipole coefficients')
    !!!      do i=1,nit
    !!!          if (i==1) then
    !!!              call yaml_map('normalized',.true.)
    !!!              call yaml_map('radius of normalization sphere',(/minval(rmax),maxval(rmax)/))
    !!!              call f_memcpy(src=multipoles, dest=multipoles_tmp)
    !!!          else if (i==2) then
    !!!              call yaml_map('normalized',.false.)
    !!!              do iat=1,nat
    !!!                  do l=0,lmax
    !!!                      do m=-l,l
    !!!                          !multipoles_tmp(m,l,iat) = multipoles(m,l,iat)/((get_normalization(rmax, l, m)*0.821583836)**2)
    !!!                          !multipoles_tmp(m,l,iat) = multipoles(m,l,iat)*((get_normalization(rmax, l, m)*0.106726871))
    !!!                          multipoles_tmp(m,l,iat) = multipoles(m,l,iat)*get_normalization(rmax(iat),l,m)**2*factor
    !!!                          max_error = max(max_error,abs(multipoles_tmp(m,l,iat)-get_test_factor(l,m)))
    !!!                          !write(*,'(a,3i5,2es14.5)') 'iat, l, m, multipoles(m,l,iat), ref', iat, l, m, multipoles(m,l,iat), get_test_factor(l,m)
    !!!                      end do
    !!!                  end do
    !!!              end do
    !!!          end if
    !!!          call yaml_sequence_open('Values')
    !!!          do iat=1,nat
    !!!              call yaml_sequence(advance='no')
    !!!              atomname=atomnames(iatype(iat))
    !!!              call yaml_sequence_open(trim(atomname))
    !!!              do l=0,lmax
    !!!                  call yaml_sequence(advance='no')
    !!!                  !call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat),fmt='(1es16.8)')
    !!!                  !call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat)*sqrt(4.d0**(2*l+3)),fmt='(1es16.8)')
    !!!                  if (check_values) then
    !!!                      max_error = max(max_error,abs(multipoles(m,l,iat)-get_test_factor(l,m)))
    !!!                  end if
    !!!                  !do m=-l,l
    !!!                      !multipoles(m,l,iat) = multipoles(m,l,iat)*get_normalization(rmax, l, m)
    !!!                      !max_error = max(max_error,abs(multipoles(m,l,iat)-get_test_factor(l,m)))
    !!!                  !end do
    !!!                  call yaml_map('l='//yaml_toa(l),multipoles_tmp(-l:l,l,iat),fmt='(1es16.8)')
    !!!                  call yaml_newline()
    !!!              end do
    !!!              !call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
    !!!              call yaml_sequence_close()
    !!!          end do
    !!!          call yaml_sequence_close()
    !!!          if (i==2) then
    !!!              call yaml_map('Maximal error from original values',max_error)
    !!!          end if
    !!!      end do
    !!!      call yaml_mapping_close()


    !!!      call f_free(multipoles_tmp)

    !!!end subroutine write_multipoles



    subroutine write_multipoles_new(nat, atomnames, rxyz, units, multipoles, &
               delta_rxyz, scaled)
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: nat
      !integer,dimension(nat),intent(in) :: iatype
      character(len=*),dimension(nat),intent(in) :: atomnames
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      character(len=*),intent(in) :: units
      real(kind=8),dimension(-lmax:lmax,0:lmax,nat),intent(in) :: multipoles
      real(kind=8),dimension(3,nat),intent(in),optional :: delta_rxyz !< can be used to display the difference between the charge center 
                                                                      !! of a support function and its localization center
      real(kind=8),dimension(nat),intent(in),optional :: scaled !< can be used to display by how muched the multipoles have been scaled
      
      ! Local variables
      character(len=20) :: atomname
      character(len=9) :: function_type
      integer :: i, iat, l, m, nit
      real(kind=8) :: factor, convert_units!, get_normalization, get_test_factor
      real(kind=8),dimension(:,:,:),allocatable :: multipoles_tmp
      logical :: present_delta_rxyz, present_scaled

      present_delta_rxyz = present(delta_rxyz)
      present_scaled = present(scaled)


          ! See whether a conversion of the units necessary
          select case (units)
          case ('angstroem','angstroemd0')
              convert_units = 0.52917721092_gp
          case ('atomic','atomicd0','bohr','bohrd0','reduced')
              convert_units = 1.d0
          case default
              convert_units = 1.d0
              call yaml_warning('units not recognized, no conversion done')
          end select


          call yaml_mapping_open('Multipole coefficients')
          call yaml_map('units for atomic positions',trim(units))
          call yaml_sequence_open('Values')
          do iat=1,nat
              call yaml_sequence(advance='no')
              atomname=trim(atomnames(iat))
              call yaml_map('sym',adjustl(trim(atomname))//' # '//adjustl(trim(yaml_toa(iat,fmt='(i4.4)'))))
              call yaml_map('r',convert_units*rxyz(1:3,iat))
              if (present_delta_rxyz) then
                  call yaml_map('Delta r',convert_units*delta_rxyz(1:3,iat),fmt='(es13.6)')
              end if
              do l=0,lmax
                  call yaml_map('q'//adjustl(trim(yaml_toa(l))),multipoles(-l:l,l,iat),fmt='(1es13.6)')
                  call yaml_newline()
              end do
              if (present_scaled) then
                  call yaml_map('scaling factor',scaled(iat),fmt='(es9.2)')
              end if
              function_type = guess_type(multipoles(:,:,iat))
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


    function get_normalization(rmax, l, m) result(fn)
      use module_base, only: pi => pi_param
      implicit none
      ! Calling arguments
      integer,intent(in) :: l, m
      real(kind=8),intent(in) :: rmax
      real(kind=8) :: fn

      ! Local variables
      integer,parameter :: l_max=2
      real(kind=8) :: r, r2, rnorm

      if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
      if (l>l_max) call f_err_throw('spherical harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
          err_name='BIGDFT_RUNTIME_ERROR')
      if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


      ! Normalization for a sphere of radius rmax
      select case (l)
      case (0)
          rnorm = sqrt(rmax**3/3.d0)
          fn = 0.5d0*sqrt(1/pi)/rnorm
      case (1)
          rnorm = sqrt(rmax**5/5.d0)
          !rnorm = sqrt(rmax**1/1.d0)
          select case (m)
          case (-1)
              fn = sqrt(3.d0/(4.d0*pi))/rnorm
          case (0)
              fn = sqrt(3.d0/(4.d0*pi))/rnorm
          case (1)
              fn = sqrt(3.d0/(4.d0*pi))/rnorm
          end select
      case (2)
          rnorm = sqrt(rmax**7/7.d0)
          select case (m)
          case (-2)
              fn = 0.5d0*sqrt(15.d0/pi)/rnorm
          case (-1)
              fn = 0.5d0*sqrt(15.d0/pi)/rnorm
          case (0)
              fn = 0.25d0*sqrt(5.d0/pi)/rnorm
          case (1)
              fn = 0.5d0*sqrt(15.d0/pi)/rnorm
          case (2)
              fn = 0.25d0*sqrt(15.d0/pi)/rnorm
          end select
      end select

    end function get_normalization


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
      if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


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


    subroutine multipole_analysis_core(iproc, nproc, natp, isat, nat, ntypes, norb, &
               r_exponent, iatype, norbsPerType, inwhichlocreg, onwhichatom, &
               n1i, n2i, n3i, nsi1, nsi2, nsi3, locrad, hgrids, locregcenter, &
               nr, psir1_get, psir2_get, &
               nvctr, matrix_compr, verbosity, &
               multipoles, rmax, get_index, smatl, matrixindex)
      use locreg_operations, only: workarr_sumrho      
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      !real(kind=8),intent(in) :: rmax
      integer,intent(in) :: natp, isat, nat, ntypes, norb, nr, nvctr, get_index, r_exponent
      integer,dimension(nat),intent(in) :: iatype
      integer,dimension(ntypes),intent(in) :: norbsPerType
      integer,dimension(norb),intent(in) :: inwhichlocreg, onwhichatom
      !logical,intent(in) :: with_rl !if true, the spherical harmonics are multiplied by r^l, otherwise not
      integer,dimension(norb),intent(in) :: n1i, n2i, n3i, nsi1, nsi2, nsi3
      real(kind=8),dimension(norb),intent(in) :: locrad
      real(kind=8),dimension(3),intent(in) :: hgrids
      real(kind=8),dimension(3,norb) :: locregcenter
      real(kind=8),dimension(nr),intent(in) :: psir1_get, psir2_get
      real(kind=8),dimension(nvctr),intent(in) :: matrix_compr
      integer,intent(in) ::verbosity
      real(kind=8),dimension(-lmax:lmax,0:lmax,nat),intent(out) :: multipoles
      real(kind=8),dimension(nat),intent(out) :: rmax
      type(sparse_matrix),intent(in),optional :: smatl
      integer,dimension(norb,norb),intent(in),optional :: matrixindex

      ! Local variables
      integer,parameter :: INDEX_AUTOMATIC=101
      integer,parameter :: INDEX_MANUALLY=102
      integer :: ist, istr, iorb, iiorb, ilr, ii, jproc, iat, n, norb_get, istr_get, i
      integer :: window, ioffset, methTransformOverlap, l, m, iiat, ityp, norb_per_atom, i1, i2, i3, ind, jorb, jat
      integer :: ii1, ii2, ii3, jjorb, iilr
      real(kind=8),dimension(:),allocatable :: psir
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      type(workarr_sumrho) :: w
      integer,dimension(:),allocatable :: nat_par, norb_list
      !real(kind=8),dimension(:),allocatable :: psir_get
      real(kind=8),dimension(:,:,:,:),allocatable :: phi1, phi2
      real(kind=8),dimension(:,:,:,:,:,:),allocatable :: sphi
      integer,dimension(:,:),allocatable :: comms
      logical :: can_use_transposed, arr_allocated
      real(kind=8) :: ddot, x, y, z, tt, rnorm, rnorm_maxdev, tt2
      !real(kind=8) ,dimension(2,orbs%norb) :: testarr
      real(kind=8),dimension(:),allocatable :: kernel_ortho, phi_ortho
      real(kind=8) :: factor_normalization
      character(len=20) :: atomname
      real(kind=8),dimension(-lmax:lmax,0:lmax) :: norm
      !real(kind=8),dimension(3) :: com
      real(kind=8) :: dnrm2
      !real(kind=8) :: rmax
      real(kind=8),dimension(7,7,4) :: khack

 khack(1:7,1,1) = (/ 2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,1) = (/-1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,1) = (/-3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,1) = (/-5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,1) = (/-5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,1) = (/ 7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,1) = (/-6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)

 khack(1:7,1,2) = (/  2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,2) = (/ -1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,2) = (/ -3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,2) = (/ -5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,2) = (/ -5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,2) = (/  7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,2) = (/ -6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)

 khack(1:7,1,3) = (/  2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,3) = (/ -1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,3) = (/ -3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,3) = (/ -5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,3) = (/ -5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,3) = (/  7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,3) = (/ -6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)

 khack(1:7,1,4) = (/  2.00E+00,-1.31E-04,-3.27E-05,-5.84E-08,-5.84E-08, 7.84E-05,-6.15E-09/)
 khack(1:7,2,4) = (/ -1.31E-04, 5.91E-01,-3.51E-01,-6.28E-04,-6.28E-04, 8.43E-01,-6.61E-05/)
 khack(1:7,3,4) = (/ -3.27E-05,-3.51E-01, 1.91E+00,-1.56E-04,-1.56E-04, 2.10E-01,-1.65E-05/)
 khack(1:7,4,4) = (/ -5.84E-08,-6.28E-04,-1.56E-04, 2.00E+00,-2.80E-07, 3.76E-04,-2.95E-08/)
 khack(1:7,5,4) = (/ -5.84E-08,-6.28E-04,-1.56E-04,-2.80E-07, 2.00E+00, 3.76E-04,-2.95E-08/)
 khack(1:7,6,4) = (/  7.84E-05, 8.43E-01, 2.10E-01, 3.76E-04, 3.76E-04, 1.50E+00, 3.95E-05/)
 khack(1:7,7,4) = (/ -6.15E-09,-6.61E-05,-1.65E-05,-2.95E-08,-2.95E-08, 3.95E-05, 2.00E+00/)


      !! Check that rmax remains within the box.
      !if (rmax>=0.5d0*(0.5d0*hgrids(1)*maxval(n1i)+0.5d0*hgrids(2)*maxval(n2i)+0.5d0*hgrids(3)*maxval(n3i))) then
      !    call f_err_throw('The radius for the multipole analysis is too small', err_name='BIGDFT_RUNTIME_ERROR')
      !end if



      rmax = 0.d0
      norb_list = f_malloc(maxval(norbsPerType(:)),id='norb_list')
      call f_zero(multipoles)
      rnorm_maxdev = 0.d0
      ist = 0
      !call mpi_barrier(mpi_comm_world,iat)
      do iat=1,natp
          iiat = iat + isat
          ityp=iatype(iiat)
          norb_per_atom = norbsPerType(ityp)
          arr_allocated = .false.
          iiorb = 0
          do iorb=1,norb
              ilr = inwhichlocreg(iorb)
              jat = onwhichatom(iorb)
              !rmax = locrad(iorb)
              if (jat==iiat) then
                  iilr = ilr
                  if (.not.arr_allocated) then
                      rmax(iiat) = min(n1i(ilr)*0.25d0*hgrids(1),n2i(ilr)*0.25d0*hgrids(2),n3i(ilr)*0.25d0*hgrids(3))
                      !rmax(iiat) = 7.d0
                      !write(*,*) 'iiat, rmax', iiat, rmax(iiat)
                      phi1 = f_malloc0((/1.to.n1i(ilr),1.to.n2i(ilr),1.to.n3i(ilr), &
                                       1.to.norb_per_atom/),id='ph1i')
                      phi2 = f_malloc0((/1.to.n1i(ilr),1.to.n2i(ilr),1.to.n3i(ilr), &
                                       1.to.norb_per_atom/),id='phi2')
                      sphi = f_malloc0((/1.to.n1i(ilr),1.to.n2i(ilr),1.to.n3i(ilr), &
                                        -lmax.to.lmax,0.to.lmax,1.to.norb_per_atom/),id='sphi')
                      !write(*,*) 'sphi allocated, natp, size',natp, size(sphi)
                      !write(*,'(a,4i9)') 'phi2 allocated, natp, n1, n2, n3',natp, n1i(ilr), n2i(ilr), n3i(ilr)
                      arr_allocated = .true.
                  else
                      ! Check the dimensions
                      if (n1i(ilr)/=size(phi1,1)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n2i(ilr)/=size(phi1,2)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n3i(ilr)/=size(phi1,3)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n1i(ilr)/=size(phi2,1)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n2i(ilr)/=size(phi2,2)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                      if (n3i(ilr)/=size(phi2,3)) call f_err_throw('wrong n1i',err_name='BIGDFT_RUNTIME_ERROR')
                  end if
                  iiorb = iiorb + 1
                  norb_list(iiorb) = iorb
                  ! Apply the spherical harmonic
                  rnorm = 0.d0
                  norm = 0.d0
                  factor_normalization = sqrt(0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3))
                  do i3=1,n3i(ilr)
                      ii3 = nsi3(ilr) + i3 - 14 - 1
                      z = ii3*0.5d0*hgrids(3) - locregcenter(3,ilr)
                      do i2=1,n2i(ilr)
                          ii2 = nsi2(ilr) + i2 - 14 - 1
                          y = ii2*0.5d0*hgrids(2) - locregcenter(2,ilr)
                          do i1=1,n1i(ilr)
                              ii1 = nsi1(ilr) + i1 - 14 - 1
                              x = ii1*0.5d0*hgrids(1) - locregcenter(1,ilr)
                              ind = (i3-1)*n2i(ilr)*n1i(ilr) + (i2-1)*n1i(ilr) + i1
                              phi1(i1,i2,i3,iiorb) = psir1_get(ist+ind)
                              phi2(i1,i2,i3,iiorb) = psir2_get(ist+ind)
                              if (x**2+y**2+z**2>rmax(iiat)**2) cycle
                              do l=0,lmax
                                  do m=-l,l
                                      tt = solid_harmonic(0, 0.d0, l, m, x, y, z)
                                      tt = tt*sqrt(4.d0*pi/real(2*l+1,kind=8))
                                      sphi(i1,i2,i3,m,l,iiorb) = tt*phi1(i1,i2,i3,iiorb)
                                      !write(*,*) 'i1, i1, i2, tt, phi1', i1, i2, i3, tt, phi1(i1,i2,i3,iiorb)
                                      norm(m,l) = norm(m,l) + (tt*factor_normalization)**2*&
                                          real((2*l+3)*(2*l+1),kind=8)/(4.d0*pi*rmax(iiat)**(2*l+3)) !normalization of a solid harmonic within a sphere of radius rmax... hopefully correct
                                      !write(*,*) 'iorb, i1, i2, i3, tt, phi', iorb, i1, i2, i3, tt, phi1(i1,i2,i3,iiorb)
                                  end do
                              end do
                          end do
                      end do
                  end do
                  do l=0,lmax
                      do m=-l,l
                          rnorm_maxdev = max(rnorm_maxdev,abs(1.d0-norm(m,l)))
                      end do
                  end do
                  ist = ist + ind
              end if
          end do
          ! Calculate the scalar product
          do l=0,lmax
              do m=-l,l
                  do iorb=1,norb_per_atom
                      iiorb = norb_list(iorb)
                      do jorb=1,norb_per_atom
                          jjorb = norb_list(jorb)
                          if (get_index==INDEX_AUTOMATIC) then
                              ind = matrixindex_in_compressed(smatl, iiorb, jjorb)
                          else if (get_index==INDEX_MANUALLY) then
                              ind = matrixindex(iiorb, jjorb)
                          end if
                          ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
                          tt = -ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), sphi(1,1,1,m,l,iorb), 1, phi2(1,1,1,jorb), 1)
                          !i = 0
                          !do i3=1,n3i(iilr)
                          !    do i2=1,n2i(iilr)
                          !        do i1=1,n1i(iilr)
                          !            i = i + 1
                          !            write(*,*) i, phi2(i1,i2,i3,jorb), sphi(i1,i2,i3,m,l,iorb)
                          !        end do
                          !    end do
                          !end do
                          ii = 0
                          do i3=1,n3i(iilr)
                              ii3 = nsi3(iilr) + i3 - 14 - 1
                              z = ii3*0.5d0*hgrids(3) - locregcenter(3,iilr)
                              do i2=1,n2i(iilr)
                                  ii2 = nsi2(iilr) + i2 - 14 - 1
                                  y = ii2*0.5d0*hgrids(2) - locregcenter(2,iilr)
                                  do i1=1,n1i(iilr)
                                      ii1 = nsi1(iilr) + i1 - 14 - 1
                                      x = ii1*0.5d0*hgrids(1) - locregcenter(1,iilr)
                                      if (x**2+y**2+z**2>rmax(iiat)**2) cycle
                                      !tt = tt + sphi(i1,i2,i3,m,l,iorb)*phi2(i1,i2,i3,jorb)
                                      ii = ii + 1
                                  end do
                              end do
                          end do
                          tt = tt
                          multipoles(m,l,iiat) = multipoles(m,l,iiat) + matrix_compr(ind)*tt
                          tt2 = ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), phi1(1,1,1,iorb), 1, phi2(1,1,1,jorb), 1)
                      end do
                  end do

              end do
          end do
          if (arr_allocated) then
              call f_free(phi1)
              call f_free(phi2)
              call f_free(sphi)
          end if
      end do
      if (nproc>1) then
          call mpiallred(multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(rnorm_maxdev, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
          call mpiallred(rmax, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (verbosity>0 .and. iproc==0) then
          call yaml_mapping_open('Calculation of solid harmonics')
          call yaml_map('Radius of integration sphere',(/minval(rmax),maxval(rmax)/),fmt='(es8.2)')
          call yaml_map('Maximal deviation from normalization',rnorm_maxdev,fmt='(es8.2)')
          call yaml_mapping_close()
      end if


      call f_free(norb_list)
    end subroutine multipole_analysis_core



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
      if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


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
               orbs, collcom, lzd, smat, locregcenter, multipole_matrix)
      use module_base
      use module_types, only: orbitals_data, comms_linear, local_zone_descriptors
      use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
      use sparsematrix_base, only: sparse_matrix, matrices
      use communications_base, only: TRANSPOSE_FULL
      use transposed_operations, only: calculate_overlap_transposed
      use communications, only: transpose_localized
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, l, m, nphi, nphir
      real(kind=8),dimension(nphi),intent(in) :: phi1, phi2
      real(kind=8),dimension(3) :: hgrids
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(inout) :: smat
      real(kind=8),dimension(3,lzd%nlr),intent(in) :: locregcenter
      type(matrices),intent(inout) :: multipole_matrix

      ! Local variables
      integer :: ist, istr, iorb, i1, i2, i3, ii1, ii2, ii3, iiorb, ind, ilr, i
      real(kind=8),dimension(:),allocatable :: phi2r, sphi2r, sphi2, phi1t_c, phi1t_f, sphi2t_c, sphi2t_f
      real(kind=8) :: norm, rmax, factor_normalization, tt, x, y, z
      type(workarr_sumrho) :: w
      real(kind=8) :: ddot


      call f_routine(id='calculte_multipole_matrix')

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


      ! Apply the spherical harmonic
      ist = 0
      do iorb=1,orbs%norbp
          iiorb = orbs%isorb + iorb
          ilr = orbs%inwhichlocreg(iiorb)
          !rmax = min(lzd%llr(ilr)%d%n1i*0.25d0*hgrids(1),lzd%llr(ilr)%d%n2i*0.25d0*hgrids(2),lzd%llr(ilr)%d%n3i*0.25d0*hgrids(3))
          rmax = min(lzd%llr(ilr)%d%n1*0.5d0*hgrids(1),lzd%llr(ilr)%d%n2*0.5d0*hgrids(2),lzd%llr(ilr)%d%n3*0.5d0*hgrids(3))
          norm = 0.d0
          factor_normalization = sqrt(0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3))
          do i3=1,lzd%llr(ilr)%d%n3i
              ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
              z = ii3*0.5d0*hgrids(3) - locregcenter(3,ilr)
              do i2=1,lzd%llr(ilr)%d%n2i
                  ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
                  y = ii2*0.5d0*hgrids(2) - locregcenter(2,ilr)
                  do i1=1,lzd%llr(ilr)%d%n1i
                      ii1 = lzd%llr(ilr)%nsi1 + i1 - 14 - 1
                      x = ii1*0.5d0*hgrids(1) - locregcenter(1,ilr)
                      ind = (i3-1)*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i + (i2-1)*lzd%llr(ilr)%d%n1i + i1
                      if (x**2+y**2+z**2>rmax**2) cycle
                      tt = solid_harmonic(0, 0.d0, l, m, x, y, z)
                      tt = tt*sqrt(4.d0*pi/real(2*l+1,kind=8))
                      sphi2r(ist+ind) = tt*phi2r(ist+ind)
                      !write(*,*) 'i1, i1, i2, tt, phi1', i1, i2, i3, tt, phi2r(ist+ind)
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
          call deallocate_work_arrays_sumrho(w)
          !write(*,*) 'iorb, n, firsts, tt', iorb, &
          !     lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, &
          !     phir(istr), sphir(istr), &
          !     ddot(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i, phir(istr), 1, sphir(istr), 1)
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

      !write(*,*) 'after overlap'

      call f_free(sphi2)
      call f_free(phi1t_c)
      call f_free(phi1t_f)
      call f_free(sphi2t_c)
      call f_free(sphi2t_f)

      call f_release_routine()

    end subroutine calculte_multipole_matrix



    subroutine multipole_analysis_driver(iproc, nproc, ll, nphi, lphi, nphir, at, hgrids, &
               orbs, smats, smatm, smatl, collcom, lzd, orthpar, ovrlp, ham, kernel, rxyz, &
               method)
      use module_base
      use module_types, only: orbitals_data, comms_linear, local_zone_descriptors, orthon_data
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, sparsematrix_malloc0, assignment(=), &
                                   sparsematrix_malloc, matrices_null, sparsematrix_malloc_ptr, deallocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix
      use communications, only: transpose_localized
      use orthonormalization, only: orthonormalizelocalized
      use module_atoms, only: atoms_data
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ll, nphi, nphir
      real(kind=8),dimension(3) :: hgrids
      type(atoms_data),intent(in) :: at
      type(orbitals_data),intent(in) :: orbs
      type(sparse_matrix),intent(inout) :: smats
      type(sparse_matrix),intent(in) :: smatm
      type(sparse_matrix),intent(inout) :: smatl
      type(comms_linear),intent(in) :: collcom
      type(local_zone_descriptors),intent(in) :: lzd
      type(orthon_data),intent(in) :: orthpar
      type(matrices),intent(inout) :: ovrlp
      type(matrices),intent(in) :: ham
      type(matrices),intent(in) :: kernel
      real(kind=8),dimension(nphi),intent(in) :: lphi
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      character(len=*),intent(in) :: method

      ! Local variables
      integer :: methTransformOverlap, iat, ind, ispin, ishift, iorb, iiorb, l, m, itype, natpx, isatx, nmaxx, kat, n, i, kkat
      integer :: ilr
      logical :: can_use_transposed
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      real(kind=8),dimension(:),allocatable :: phi_ortho, Qmat, kernel_ortho, multipole_matrix_large
      real(kind=8),dimension(:,:),allocatable :: Qmat_tilde, kp, locregcenter
      real(kind=8),dimension(:,:,:),allocatable :: atomic_multipoles
      real(kind=8),dimension(:,:),pointer :: projx
      real(kind=8) :: q, tt
      type(matrices) :: multipole_matrix
      logical,dimension(:,:),pointer :: neighborx
      integer,dimension(:),pointer :: nx
      character(len=20),dimension(:),allocatable :: names

      call f_routine(id='multipole_analysis_driver')

      if (iproc==0) then
          call yaml_comment('Atomic multipole analysis, new approach',hfill='=')
          call yaml_map('Method',trim(method))
      end if

      call unitary_test_multipoles(iproc, nproc, nphi, nphir, orbs, lzd, smats, collcom, hgrids)

      ! Check the method
      if (trim(method)/='projector' .and. trim(method)/='loewdin') then
          call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
      end if


      if (trim(method)=='projector') then
          ! Calculate the projector using the penalty term
          call projector_for_charge_analysis(at, smats, smatm, smatl, &
               ovrlp, ham, kernel, rxyz, calculate_centers=.false., write_output=.false., &
               natpx=natpx, isatx=isatx, nmaxx=nmaxx, nx=nx, projx=projx, neighborx=neighborx)
      end if

      ! Calculate the kernel for orthormal support functions
      methTransformOverlap = 20
      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='kernel_ortho')
      call kernel_for_orthonormal_basis(iproc, nproc, orbs%norbp, methTransformOverlap, smats, smatl, &
           ovrlp, kernel, kernel_ortho)

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

      !!i=1
      !!do iorb=1,orbs%norb
      !!   ilr = orbs%inwhichlocreg(iorb)
      !!   n = lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      !!   write(*,*) 'debug iorb tt', iorb, dot(n, phi_ortho(i), 1, phi_ortho(i), 1) 
      !!   i = i + n
      !!end do

      Qmat = sparsematrix_malloc(smatl,iaction=SPARSE_FULL,id='Qmat')
      atomic_multipoles = f_malloc0((/-ll.to.ll,0.to.ll,1.to.at%astruct%nat/),id='atomic_multipoles')

      multipole_matrix = matrices_null()
      multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(smats, SPARSE_FULL, id='multipole_matrix%matrix_compr')

      multipole_matrix_large = sparsematrix_malloc(smatl, SPARSE_FULL, id='multipole_matrix_large')

      locregcenter = f_malloc((/3,lzd%nlr/),id='locregcenter')
      do ilr=1,lzd%nlr
          locregcenter(1:3,ilr) = lzd%llr(ilr)%locregcenter(1:3)
      end do

      do l=0,ll
          do m=-l,l

              call f_zero(multipole_matrix%matrix_compr)

              ! Calculate the multipole matrix
              call calculte_multipole_matrix(iproc, nproc, l, m, nphi, phi_ortho, phi_ortho, nphir, hgrids, &
                   orbs, collcom, lzd, smats, locregcenter, multipole_matrix)

              call transform_sparse_matrix(smats, smatl, multipole_matrix%matrix_compr, multipole_matrix_large, 'small_to_large')


              ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
              call dscal(smatl%nvctr*smatl%nspin, -1.d0, multipole_matrix_large(1), 1)

              ! Multiply the orthogonalized kernel with the multipole matrix
              call f_zero(Qmat)
              call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                   kernel_ortho, multipole_matrix_large, Qmat)
              !write(*,*) 'after matrix_matrix_mult_wrapper'

              if (trim(method)=='projector') then
                  do kat=1,natpx
                      kkat = kat + isatx
                      n = nx(kat)
                      Qmat_tilde = f_malloc((/n,n/),id='Qmat_tilde')
                      kp = f_malloc((/n,n/),id='kp')
                      call extract_matrix(smatl, Qmat, neighborx(1:,kat), n, nmaxx, Qmat_tilde)
                      call gemm('n', 'n', n, n, n, 1.d0, Qmat_tilde(1,1), n, projx(1,kat), n, 0.d0, kp(1,1), n)
                      tt = 0.d0
                      do i=1,n
                          tt = tt + kp(i,i)
                      end do
                      atomic_multipoles(m,l,kkat) = tt
                      call f_free(Qmat_tilde)
                      call f_free(kp)
                  end do
              else if (trim(method)=='loewdin') then
                  do ispin=1,smatl%nspin
                      ishift = (ispin-1)*smatl%nvctr
                      do iorb=1,orbs%norb
                          iiorb = modulo(iorb-1,smatl%nfvctr)+1
                          iat=smatl%on_which_atom(iiorb)
                          ind = matrixindex_in_compressed(smatl, iorb, iorb)
                          ind = ind + ishift
                          !if (iproc==0) write(*,*) 'iorb, trace charge', iorb, weight_matrix_compr(ind)
                          !!write(*,*) 'l, m, iorb, iat, ind, val, val2', &
                          !!    l, m, iorb, iat, ind, Qmat(ind), multipole_matrix_large(ind)
                          atomic_multipoles(m,l,iat) = atomic_multipoles(m,l,iat) + Qmat(ind)
                      end do
                  end do
              end if

          end do
      end do

      call f_free(locregcenter)

      if (trim(method)=='projector') then
          call mpiallred(atomic_multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if


      ! The monopole term should be the net charge, i.e. add the positive atomic charges
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          q = real(at%nelpsp(itype),kind=8)
          atomic_multipoles(0,0,iat) = atomic_multipoles(0,0,iat) + q
      end do

      if (iproc==0) then
          call yaml_comment('Final result of the multipole analysis',hfill='~')
          names = f_malloc_str(len(names),at%astruct%nat,id='names')
          do iat=1,at%astruct%nat
              itype = at%astruct%iatype(iat)
              names(iat) = at%astruct%atomnames(itype)
          end do
          call write_multipoles_new(at%astruct%nat, names, &
               at%astruct%rxyz, at%astruct%units, &
               atomic_multipoles)
          call f_free_str(len(names),names)
      end if

      call deallocate_matrices(multipole_matrix)
      call f_free(kernel_ortho)
      call f_free(Qmat)
      call f_free(phi_ortho)
      if (trim(method)=='projector') then
          call f_free_ptr(projx)
          call f_free_ptr(nx)
          call f_free_ptr(neighborx)
      end if
      call f_free(atomic_multipoles)
      call f_free(multipole_matrix_large)

      if (iproc==0) then
          call yaml_comment('Atomic multipole analysis done',hfill='=')
      end if

      call f_release_routine()

  end subroutine multipole_analysis_driver



    subroutine projector_for_charge_analysis(at, smats, smatm, smatl, &
               ovrlp_, ham_, kernel_, rxyz, calculate_centers, write_output, &
               lzd, nphirdim, psi, orbs, &
               multipoles, &
               natpx, isatx, nmaxx, nx, projx, neighborx)
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
      type(sparse_matrix),intent(inout) :: smats, smatl !< should be intent(in)...
      type(sparse_matrix),intent(in) :: smatm
      type(matrices),intent(inout) :: ovrlp_ !< should be intent(in)...
      type(matrices),intent(in) :: ham_, kernel_
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      logical,intent(in) :: calculate_centers, write_output
      type(local_zone_descriptors),intent(in),optional :: lzd
      integer,intent(in),optional :: nphirdim
      real(kind=8),dimension(:),intent(in),optional :: psi
      type(orbitals_data),intent(in),optional :: orbs
      real(kind=8),dimension(-lmax:lmax,0:lmax,1:smats%nfvctr),intent(in),optional :: multipoles
      integer,intent(out),optional :: natpx, isatx, nmaxx
      integer,dimension(:),pointer,intent(out),optional :: nx
      real(kind=8),dimension(:,:),pointer,intent(out),optional :: projx
      logical,dimension(:,:),pointer,intent(out),optional :: neighborx

      ! Local variables
      integer :: kat, iat, jat, i, j, ii, jj, icheck, n, indm, inds, ntot, ist, ind, iq, itype, ieval, ij, nmax, indl, lwork
      integer :: k, l, iatold, isat, natp, kkat, istot, ntotp, i1, i2, i3, is1, ie1, is2, ie2, is3, ie3, j1, j2, j3, ikT, info
      integer :: ialpha
      real(kind=8) :: r2, cutoff2, rr2, tt, ef, q, occ, max_error, mean_error, rr2i, rr2j, ttxi, ttyi, ttzi, ttxj, ttyj, ttzj
      real(kind=8) :: tti, ttj, charge_net, charge_total
      real(kind=8) :: xi, xj, yi, yj, zi, zj, ttx, tty, ttz, xx, yy, zz, x, y, z
      real(kind=8),dimension(:),allocatable :: work
      real(kind=8),dimension(:,:),allocatable :: com
      real(kind=8),dimension(:,:),allocatable :: ham, ovrlp, proj, ovrlp_tmp, ovrlp_minusonehalf, kp, ktilde
      real(kind=8),dimension(:,:,:),allocatable :: coeff_all, ovrlp_onehalf_all
      integer,dimension(:,:,:,:),allocatable :: ilup
      real(kind=8),dimension(:),allocatable :: eval, eval_all, ovrlp_large, tmpmat1, tmpmat2, kerneltilde, charge_per_atom
      real(kind=8),dimension(:,:,:),allocatable :: tmpmat2d
      integer,dimension(:),allocatable :: id_all, n_all, itmparr
      real(kind=8),dimension(3) :: rr
      logical,dimension(:,:),allocatable :: neighbor
      type(matrices),dimension(1) :: ovrlp_onehalf_
      logical :: perx, pery, perz, final, bound_low_ok, bound_up_ok
      !real(kind=8),parameter :: kT = 5.d-2
      real(kind=8) :: kT
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


      ! Calculate S^1/2
      ovrlp_onehalf_(1) = matrices_null()
      ovrlp_onehalf_(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_onehalf_(1)%matrix_compr')
      call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, 1020, 1, (/2/), -1, &
            imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
            ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_(1), &
            check_accur=.true., max_error=max_error, mean_error=mean_error)

      ! Calculate S^1/2 * K * S^1/2 = Ktilde
      tmpmat1 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat1')
      !tmpmat2 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat2')
      kerneltilde = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='kerneltilde')
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           kernel_%matrix_compr, ovrlp_onehalf_(1)%matrix_compr, tmpmat1)
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           ovrlp_onehalf_(1)%matrix_compr, tmpmat1, kerneltilde)


      ! Determine the periodicity...
      perx=(smats%geocode /= 'F')
      pery=(smats%geocode == 'P')
      perz=(smats%geocode /= 'F')
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

      alpha_loop: do ialpha=1,100

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
              eval = f_malloc0((/n/),id='eval')
              call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1:,kat), n, nmax, ovrlp, ilup)
              call extract_matrix(smatm, ham_%matrix_compr, neighbor(1:,kat), n, nmax, ham)
              !!icheck = 0
              !!ii = 0
              !!do i=1,smats%nfvctr
              !!    if (neighbor(i,kat)) then
              !!        jj = 0
              !!        do j=1,smats%nfvctr
              !!            if (neighbor(j,kat)) then
              !!                icheck = icheck + 1
              !!                jj = jj + 1
              !!                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
              !!                inds =  matrixindex_in_compressed(smats, i, j)
              !!                if (inds>0) then
              !!                    ovrlp(jj,ii) = ovrlp_%matrix_compr(inds)
              !!                else
              !!                    ovrlp(jj,ii) = 0.d0
              !!                end if
              !!                indm =  matrixindex_in_compressed(smatm, i, j)
              !!                if (indm>0) then
              !!                    ham(jj,ii) = ham_%matrix_compr(indm)
              !!                else
              !!                    ham(jj,ii) = 0.d0
              !!                end if
              !!            end if
              !!        end do
              !!    end if
              !!end do
              !!if (icheck>n**2) then
              !!    call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
              !!        &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
              !!end if


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

              ! Add the penalty term
              if (mode=='old') then
                  call add_penalty_term(smats%geocode, smats%nfvctr, neighbor(1:,kat), rxyz(1:,kkat), &
                       at%astruct%cell_dim, com, alpha, n, ovrlp, ham)
              else if (mode=='new') then
                  multipoles_fake = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.smats%nfvctr/),id='multipoles_fake')
                  multipoles_fake = 0.d0
                  multipoles_fake(0,0,:) = 1.d0
                  if (ialpha==1) then
                      if (present(multipoles)) then
                          write(*,*) 'call with multipoles'
                          call add_penalty_term_new(smats%geocode, at%astruct%nat, smats%nfvctr, &
                               neighbor(1:,kat), rxyz(1:,kkat), smats%on_which_atom, &
                               multipoles, at%astruct%cell_dim, com, alpha, n, ham, &
                               nmax, penalty_matrices(1:n,1:n,kat))
                      else
                          write(*,*) 'call with multipoles_fake'
                          call add_penalty_term_new(smats%geocode, at%astruct%nat, smats%nfvctr, &
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

              !!icheck = 0
              !!ii = 0
              !!do i=1,smats%nfvctr
              !!    if (neighbor(i,kat)) then
              !!        jj = 0
              !!        do j=1,smats%nfvctr
              !!            if (neighbor(j,kat)) then
              !!                icheck = icheck + 1
              !!                jj = jj + 1
              !!                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
              !!                if (i==j) then
              !!                    rr2 = huge(rr2)
              !!                    do i3=is3,ie3
              !!                        z = rxyz(3,kkat) + i3*at%astruct%cell_dim(3)
              !!                        ttz = (com(3,i)-z)**2
              !!                        do i2=is2,ie2
              !!                            y = rxyz(2,kkat) + i2*at%astruct%cell_dim(2)
              !!                            tty = (com(2,i)-y)**2
              !!                            do i1=is1,ie1
              !!                                x = rxyz(1,kkat) + i1*at%astruct%cell_dim(1)
              !!                                ttx = (com(1,i)-x)**2
              !!                                tt = ttx + tty + ttz
              !!                                if (tt<rr2) then
              !!                                    rr2 = tt
              !!                                end if
              !!                            end do
              !!                        end do
              !!                    end do
              !!                    ham(jj,ii) = ham(jj,ii) + alpha*rr2**3*ovrlp(jj,ii)
              !!                end if
              !!                ilup(1,jj,ii,kat) = j
              !!                ilup(2,jj,ii,kat) = i
              !!            end if
              !!        end do
              !!    end if
              !!end do
              !!if (icheck>n**2) then
              !!    call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
              !!        &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
              !!end if
    
    
              !!call diagonalizeHamiltonian2(bigdft_mpi%iproc, n, ham, ovrlp, eval)
              lwork = 10*n
              work = f_malloc(lwork,id='work')
              call syev('v', 'l', n, ham(1,1), n, eval(1), work(1), lwork, info)
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
          !do i=1,ntot
          !    ! add i-1 since we are only searching in the subarray
          !    ind = minloc(eval_all(i:ntot),1) + (i-1)
          !    tt = eval_all(i)
          !    eval_all(i) = eval_all(ind)
          !    eval_all(ind) = tt
          !    ii = id_all(i)
          !    id_all(i) = id_all(ind)
          !    id_all(ind) = ii
          !end do
        
        


          if (ialpha==1) then
              if (present(nmaxx)) nmaxx = maxval(n_all)
              if (present(projx)) projx = f_malloc_ptr((/nmaxx**2,natpx/),id='projx')
          end if
    
    
          !!! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
          !!ef = eval_all(1)
          !!do
          !!    ef = ef + 1.d-3
          !!    occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
          !!    if (abs(occ-1.d0)<1.d-8) exit
          !!end do
          !!if (bigdft_mpi%iproc==0) then
          !!    call yaml_map('Pseudo Fermi level for occupations',ef)
          !!end if
        
          !!final = .true.
          !!ikT = 0
          !!kT_loop: do
    
              !ikT = ikT + 1
    
              call f_zero(charge_per_atom)
    
              ! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
              ef = eval_all(1)
              do
                  ef = ef + 1.d-3
                  occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
                  if (abs(occ-1.d0)<1.d-8) exit
              end do
        
              ! Calculate the projector. First for each single atom, then insert it into the big one.
              charge_total = 0.d0
              do kat=1,natp
                  kkat = kat + isat
                  n = n_all(kat)
                  proj = f_malloc0((/n,n/),id='proj')
                  call calculate_projector(n, ntot, nmax, kkat, id_all, eval_all, &
                       coeff_all(1:,1:,kat), ef, kT, proj)
                  if (present(projx)) then
                      call vcopy(n**2, proj(1,1), 1, projx(1,kat), 1)
                  end if
                  !ij = 0
                  !do ieval=1,ntot
                  !    if (id_all(ieval)/=kkat) cycle
                  !    ij = ij + 1
                  !    occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
                  !    do i=1,n
                  !        do j=1,n
                  !            proj(j,i) = proj(j,i) + occ*coeff_all(j,ij,kat)*coeff_all(i,ij,kat)
                  !        end do
                  !   end do
                  !end do
                  !tt = 0.d0
                  !do i=1,n
                  !    tt = tt + proj(i,i)
                  !end do
    
    
                  !@ TEMPORARY ############################################
                  ! Extract ktilde
                  ktilde = f_malloc0((/n,n/),id='ktilde')
                  call extract_matrix(smatl, kerneltilde, neighbor(1:,kat), n, nmax, ktilde)
                  kp = f_malloc((/n,n/),id='kp')
                  !ii = 0
                  !do i=1,smats%nfvctr
                  !    if (neighbor(i,kat)) then
                  !        jj = 0
                  !        do j=1,smats%nfvctr
                  !            if (neighbor(j,kat)) then
                  !                icheck = icheck + 1
                  !                jj = jj + 1
                  !                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                  !                indl =  matrixindex_in_compressed(smatl, i, j)
                  !                if (indl>0) then
                  !                    ktilde(jj,ii) = kerneltilde(indl)
                  !                else
                  !                    ktilde(jj,ii) = 0.d0
                  !                end if
                  !            end if
                  !        end do
                  !    end if
                  !end do
    
                  ! Calculate ktilde * proj
                  !write(*,*) 'ktilde'
                  !do i=1,n
                  !    write(*,'(10es9.2)') ktilde(i,1:n)
                  !end do
                  !!write(*,*) 'WARNING: COPY KHACK TO KTILDE'
                  !!ktilde = khack(1:7,1:7,kat)
                  call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, proj(1,1), n, 0.d0, kp(1,1), n)
                  tt = 0
                  do i=1,n
                      tt = tt + kp(i,i)
                  end do
                  !!if (bigdft_mpi%iproc==0) then
                  !!    do i=1,n
                  !!        do j=1,n
                  !!            write(*,'(a,2i5,3es13.3)') 'i, j, kt, proj, kp', i, j, ktilde(i,j), proj(j,i), kp(j,i)
                  !!        end do
                  !!    end do
                  !!    write(*,*) 'kkat, trace, sum(proj)', kkat, tt, sum(proj)
                  !!end if
                  charge_per_atom(kkat) = tt
                  !write(*,*) 'alpha, kkat, tt', alpha, kkat, tt
                  charge_total = charge_total + tt
                  call f_free(proj)
                  call f_free(ktilde)
                  call f_free(kp)
    
    
              end do
    
    
              !!if (final) exit kT_loop
    
              !!charge_net = 0.d0
              !!do iat=1,at%astruct%nat
              !!    charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
              !!end do
              !!!!if (bigdft_mpi%iproc==0) then
              !!!!    call yaml_map('kT, ef, net_charge',(/kT,ef,charge_net/))
              !!!!end if
              !!if (abs(charge_net)<1.d0 .or. ikT==100) then
              !!    final = .true.
              !!else if (charge_net<0.d0) then
              !!    kT = kT*0.95d0
              !!    !ef = ef + 1.d-3
              !!else if (charge_net>0.d0) then
              !!    kT = kT*1.05d0
              !!    !ef = ef - 1.d-3
              !!end if
    
          !!end do kT_loop

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


          if (abs(charge_net)<convergence_criterion) then
          !if (.true.) then
              if (bigdft_mpi%iproc==0) then
                  call yaml_sequence_close()
                  call yaml_map('number of states to be occupied (without smearing)',iq)
                  call yaml_map('Pseudo Fermi level for occupations',ef)
                  call yaml_sequence_open('ordered eigenvalues and occupations')
                  ii = 0
                  do i=1,ntot
                      occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
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

      call deallocate_matrices(ovrlp_onehalf_(1))
      call f_free(tmpmat1)
      !call f_free(tmpmat2)
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
      call f_free(id_all)
      call f_free(ovrlp_onehalf_all)
      call f_free(com)

      call f_release_routine()

  end subroutine projector_for_charge_analysis



  !!subroutine determine_atomic_charges(smat, nat, effective_kernel, charge_per_atom)
  !!  use module_base
  !!  use sparsematrix_base, only: sparse_matrix
  !!  use sparsematrix_init, only: matrixindex_in_compressed
  !!  implicit none
  !!
  !!  ! Calling arguments
  !!  type(sparse_matrix),intent(in) :: smat
  !!  integer,intent(in) :: nat
  !!  real(kind=8),dimension(smat%nvctr*smat%nspin),intent(in) :: effective_kernel
  !!  real(kind=8),dimension(nat),intent(out) :: charge_per_atom
  !!
  !!  ! Local variables
  !!  integer :: ispin, ishift, iorb, iiorb, ind, iat
  !!
  !!  call f_zero(charge_per_atom)
  !!  
  !!  do ispin=1,smat%nspin
  !!      ishift = (ispin-1)*smat%nvctr
  !!      do iorb=1,smat%nfvctr
  !!          iiorb = modulo(iorb-1,smat%nfvctr)+1
  !!          iat = smat%on_which_atom(iiorb)
  !!          ind = matrixindex_in_compressed(smat, iorb, iorb)
  !!          ind = ind + ishift
  !!          !write(*,*) 'iorb, ind, val', iorb, ind, effective_kernel(ind)
  !!          charge_per_atom(iat) = charge_per_atom(iat) + effective_kernel(ind)
  !!      end do
  !!  end do
  !!
  !!end subroutine determine_atomic_charges


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
                        ham(jj,ii) = ham(jj,ii) + alpha*rr2**3*ovrlp(jj,ii)
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


 subroutine calculate_projector(n, ntot, nmax, kkat, ids, evals, coeff, ef, kT, proj)
   use module_base
   implicit none

   ! Calling arguments
   integer :: n, ntot, nmax, kkat
   integer,dimension(ntot),intent(in) :: ids
   real(kind=8),dimension(ntot),intent(in) :: evals
   real(kind=8),dimension(nmax,nmax),intent(in) :: coeff
   real(kind=8),intent(in) :: kT, ef
   real(kind=8),dimension(n,n),intent(out) :: proj

   ! Local variables
   integer :: ij, ieval, i, j
   real(kind=8) :: occ

   call f_routine(id='calculate_projector')

   ij = 0
   do ieval=1,ntot
       if (ids(ieval)/=kkat) cycle
       ij = ij + 1
       occ = 1.d0/(1.d0+safe_exp( (evals(ieval)-ef)*(1.d0/kT) ) )
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
   use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, assignment(=), &
                                matrices_null, sparsematrix_malloc_ptr, deallocate_matrices
   use locreg_operations,only: workarr_sumrho, initialize_work_arrays_sumrho, deallocate_work_arrays_sumrho
   use sparsematrix_init, only: matrixindex_in_compressed
   use yaml_output
   implicit none
   ! Calling arguments
   integer,intent(in) :: iproc, nproc, nphi, nphir
   type(orbitals_data),intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: lzd
   type(sparse_matrix),intent(inout) :: smat
   type(comms_linear),intent(in) :: collcom
   real(kind=8),dimension(3) :: hgrids
   ! Local variables
   integer :: iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, l, m, i, ind, ist, istr, ii
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
   multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(smat, SPARSE_FULL, id='multipole_matrix%matrix_compr')

   phi2r = f_malloc0(nphir,id='phi2r')
   phi1r = f_malloc(nphir,id='phi1r')
   phi1 = f_malloc0(nphi,id='phi1')

   locregcenter = f_malloc0((/3,lzd%nlr/),id='locregcenter')

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
           ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
           z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3,ilr)
           do i2=1,lzd%llr(1)%d%n2i
               ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
               y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2,ilr)
               do i1=1,lzd%llr(1)%d%n1i
                   ii1 = lzd%llr(ilr)%nsi1 + i1 - 14 - 1
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
                    orbs, collcom, lzd, smat, locregcenter, multipole_matrix)
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
 subroutine support_function_gross_multipoles(iproc, nproc, tmb, atoms, denspot, multipoles)
   use module_base
   use module_types
   use locreg_operations
   use yaml_output
   use multipole_base, only: lmax
   use bounds, only: geocode_buffers
   use sparsematrix_base, only: matrices, matrices_null, sparsematrix_malloc_ptr, SPARSE_FULL, assignment(=), &
                                deallocate_matrices
   use sparsematrix_init, only: matrixindex_in_compressed
   
   ! Calling arguments
   integer,intent(in) :: iproc, nproc
   type(DFT_wavefunction),intent(inout) :: tmb
   type(atoms_data),intent(in) :: atoms
   type(DFT_local_fields), intent(inout) :: denspot
 
   integer :: ist, istr, iorb, iiorb, ilr, i, iat, iter, itype
   integer :: i1, i2, i3, ii1, ii2, ii3, nl1, nl2, nl3, ii, l, m, ind, iat_old
   real(kind=8),dimension(:),allocatable :: rmax, phi1, phi1r
   real(kind=8),dimension(:,:),allocatable :: delta_centers
   real(kind=8),dimension(:,:),allocatable :: charge_center
   real(kind=8),dimension(:),allocatable :: phir, phir_one
   real(kind=8) :: hxh, hyh, hzh, tt, x, y, z, weight
   type(workarr_sumrho) :: w
   character(len=20) :: atomname
   character(len=20),dimension(:),allocatable :: names
   integer,dimension(:),allocatable :: iatype_tmp
   type(matrices) :: multipole_matrix
   real(kind=8),dimension(:,:,:),allocatable :: multipoles
   real(kind=8),dimension(:),allocatable :: scaled
 
   call f_routine(id='support_function_multipoles')
 
   rmax = f_malloc0(tmb%orbs%norb,id='rmax')
   phir = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir')
   phir_one = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir_one')
   phir_one = 1.d0
 
   !call to_zero(3*tmb%orbs%norb, dipole_net(1,1))
   !call to_zero(9*tmb%orbs%norb, quadropole_net(1,1,1))


   charge_center = f_malloc0((/3,tmb%lzd%nlr/),id='charge_center')
   multipole_matrix = matrices_null()
   multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%s, SPARSE_FULL, id='multipole_matrix%matrix_compr')

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
      call daub_to_isf(tmb%lzd%llr(ilr), w, tmb%psi(ist), phir(istr))
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
                  charge_center(1,iiorb) = charge_center(1,iiorb) + x*tt
                  charge_center(2,iiorb) = charge_center(2,iiorb) + y*tt
                  charge_center(3,iiorb) = charge_center(3,iiorb) + z*tt
                  weight = weight + tt
                  ii = ii + 1
              end do
          end do
      end do
      charge_center(1:3,iiorb) = charge_center(1:3,iiorb)/weight
      !write(*,*) 'ilr, charge_center(1:3,iorb), lzd%llr(ilr)%locregcenter(1:3)', &
      !            ilr, charge_center(1:3,iorb), tmb%lzd%llr(ilr)%locregcenter(1:3)
      ! ######################################################################################
      ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i
  end do

  if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
      stop
  end if

  if (nproc>1) then
      call mpiallred(charge_center, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if


  do l=0,lmax
      do m=-l,l
          call f_zero(multipole_matrix%matrix_compr)
          ! Calculate the multipole matrix
          call calculte_multipole_matrix(iproc, nproc, l, m, tmb%npsidim_orbs, phi1, tmb%psi, &
               max(tmb%collcom_sr%ndimpsi_c,1), tmb%lzd%hgrids, &
               tmb%orbs, tmb%collcom, tmb%lzd, tmb%linmat%s, charge_center, multipole_matrix)
          ! Take the diagonal elements
          do iorb=1,tmb%orbs%norbp
              iiorb = tmb%orbs%isorb + iorb
              ind = matrixindex_in_compressed(tmb%linmat%s, iiorb, iiorb)
              multipoles(m,l,iiorb) = multipole_matrix%matrix_compr(ind)
          end do
      end do
  end do

  ! Normalize the multipoles such that the largest component has the magnitude 1
  do iorb=1,tmb%orbs%norbp
      iiorb = tmb%orbs%isorb + iorb
      tt = maxval(abs(multipoles(:,:,iiorb)))
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
      iatype_tmp = f_malloc(tmb%orbs%norb,id='iatype_tmp')
      delta_centers = f_malloc((/3,tmb%orbs%norb/),id='delta_centers')
      iat_old = -1
      names = f_malloc_str(len(names),tmb%orbs%norb,id='names')
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
          delta_centers(1:3,iorb) = charge_center(1:3,iorb) - tmb%lzd%llr(ilr)%locregcenter(1:3)
          !write(*,*) 'iorb, delta_centers(1:3,iorb)', iorb, delta_centers(1:3,iorb)
      end do
      call write_multipoles_new(tmb%orbs%norb, names, charge_center, &
           atoms%astruct%units, multipoles, &
           delta_centers, scaled)
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
  call f_free(charge_center)
  call f_free(multipoles)

  call f_release_routine()
 
 end subroutine support_function_gross_multipoles


end module multipole
