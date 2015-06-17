module multipole
  use module_base
  use multipole_base, only: external_potential_descriptors,lmax
  implicit none

  private

  !> Public routines
  public :: interaction_multipoles_ions
  public :: potential_from_charge_multipoles
  public :: potential_from_multipoles
  public :: multipoles_from_density

  contains


    subroutine interaction_multipoles_ions(ep, at, eion, fion)
      use module_types, only: atoms_data
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      type(atoms_data),intent(in) :: at
      real(gp),intent(inout) :: eion
      real(gp),dimension(3,at%astruct%nat),intent(inout) :: fion

      ! Local variables
      integer :: iat, ityp, impl
      real(gp) :: r, charge

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
                  eion = eion + charge/r
                  fion(1,iat) = fion(1,iat) + charge/(r**3)*(at%astruct%rxyz(1,iat)-ep%mpl(impl)%rxyz(1))
                  fion(2,iat) = fion(2,iat) + charge/(r**3)*(at%astruct%rxyz(2,iat)-ep%mpl(impl)%rxyz(2))
                  fion(3,iat) = fion(3,iat) + charge/(r**3)*(at%astruct%rxyz(3,iat)-ep%mpl(impl)%rxyz(3))
              end if
          end do
      end do

    end subroutine interaction_multipoles_ions


    !> Calculate the external potential arising from the multipoles of the charge density
    subroutine potential_from_charge_multipoles(iproc, nproc, denspot, ep, is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, shift, pot)
      use module_types, only: DFT_local_fields
      use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(DFT_local_fields),intent(inout) :: denspot
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3
      real(gp),intent(in) :: hx, hy, hz
      real(gp),dimension(3),intent(in) :: shift !< global shift of the atomic positions
      real(gp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l, m, ii, mm
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp, ehart_ps, tt, ttt
      real(dp),dimension(3) :: r
      real(dp),dimension(:,:,:),allocatable :: density
      real(kind=8),dimension(0:lmax) :: sigma
      real(8),dimension(ep%nmpl) :: monopole
      real(8),dimension(0:2,1:ep%nmpl) :: norm
      real(8),dimension(3,ep%nmpl) :: dipole
      real(8),dimension(5,ep%nmpl) :: quadrupole

      call f_routine(id='potential_from_charge_multipoles')

      sigma(0) = 5.d0*(hx*hy*hz)**(1.d0/3.d0)
      sigma(1) = 4.d0*(hx*hy*hz)**(1.d0/3.d0)
      sigma(2) = 2.d0*(hx*hy*hz)**(1.d0/3.d0)

      density = f_malloc0((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='density')

      !!$omp parallel &
      !!$omp default(none) &
      !!$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, ep, pot) &
      !!$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, rnrm1, rnrm2, rnrm3, rnrm5, l, mp)
      !!$omp do
      norm = 0.d0
      monopole = 0.d0
      dipole = 0.d0
      quadrupole = 0.d0
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
                  tt = 0.d0
                  do impl=1,ep%nmpl
                      r(1) = x - ep%mpl(impl)%rxyz(1)
                      r(2) = y - ep%mpl(impl)%rxyz(2)
                      r(3) = z - ep%mpl(impl)%rxyz(3)
                      rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                      rnrm1 = sqrt(rnrm2)
                      !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
                      !if (i1==is1 .and. i2==is2)  then
                      !    write(*,'(a,i5,3es16.7)') 'impl, rnrm1, g, norm', impl, rnrm1, gaussian(sigma, rnrm1), norm(impl)
                      !end if
                      do l=0,lmax
                          norm(l,impl) = norm(l,impl) + gaussian(sigma(l), rnrm1)*hx*hy*hz
                          if (associated(ep%mpl(impl)%qlm(l)%q)) then
                              mm = 0
                              do m=-l,l
                                  mm = mm + 1
                                  !if (l==0) then
                                  !    charge(impl) = charge(impl) + ep%mpl(impl)%qlm(l)%q(mm)*spherical_harmonic(l, m, r(1), r(2), r(3))*gaussian(sigma, rnrm1)*hx*hy*hz*sqrt(4.d0*pi_param)
                                  !end if
                                  !if (l==1) then
                                  !    dipole(mm,impl) = dipole(mm,impl) + ep%mpl(impl)%qlm(l)%q(mm)*spherical_harmonic(l, m, r(1), r(2), r(3))*gaussian(sigma, rnrm1)*hx*hy*hz*sqrt(4.d0*pi_param)
                                  !end if
                                  ttt = ep%mpl(impl)%qlm(l)%q(mm)*spherical_harmonic(l, m, r(1), r(2), r(3))*gaussian(sigma(l), rnrm1)*sqrt(4.d0*pi_param)
                                  tt = tt + ttt

                                  if (l==0) then
                                      monopole(impl) = monopole(impl) + ttt*hx*hy*hz
                                  else if (l==1) then
                                      if (m==-1) then
                                          dipole(1,impl) = dipole(1,impl) + ttt*hx*hy*hz*y
                                      else if (m==0) then
                                          dipole(2,impl) = dipole(2,impl) + ttt*hx*hy*hz*z
                                      else if (m==1) then
                                          dipole(3,impl) = dipole(3,impl) + ttt*hx*hy*hz*x
                                      end if
                                  else if (l==2) then
                                      if (m==-2) then
                                          quadrupole(1,impl) = quadrupole(1,impl) + ttt*hx*hy*hz*x*y
                                      else if (m==-1) then
                                          quadrupole(2,impl) = quadrupole(2,impl) + ttt*hx*hy*hz*y*z
                                      else if (m==0) then
                                          quadrupole(3,impl) = quadrupole(3,impl) + ttt*hx*hy*hz*(3*z**2-1.d0)
                                      else if (m==1) then
                                          quadrupole(4,impl) = quadrupole(4,impl) + ttt*hx*hy*hz*x*z
                                      else if (m==2) then
                                          quadrupole(5,impl) = quadrupole(5,impl) + ttt*hx*hy*hz*(x**2-y**2)
                                      end if
                                  end if
                              end do
                          end if
                      end do
                  end do
                  density(i1,i2,i3) = tt
                  !!do impl=1,ep%nmpl
                  !!    do l=0,lmax
                  !!        if (associated(ep%mpl(impl)%qlm(l)%q)) then
                  !!            mm = 0
                  !!            do m=-l,l
                  !!                mm = mm + 1
                  !!                if (l==0) then
                  !!                    charge(impl) = charge(impl) + tt*hx*hy*hz
                  !!                else if (l==1) then
                  !!                    if (m==-1) then
                  !!                        dipole(1,impl) = dipole(1,impl) + tt*hx*hy*hz*y
                  !!                    else if (m==0) then
                  !!                        dipole(2,impl) = dipole(1,impl) + tt*hx*hy*hz*z
                  !!                    else if (m==1) then
                  !!                        dipole(3,impl) = dipole(1,impl) + tt*hx*hy*hz*x
                  !!                    end if
                  !!                end if
                  !!            end do
                  !!        end if
                  !!    end do
                  !!end do
              end do
          end do
      end do
      if (nproc>0) then
          call mpiallred(norm, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(dipole, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(quadrupole, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (iproc==0) then
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
          call yaml_map('number of multipole centers',ep%nmpl)
          call yaml_map('sigma of the Gaussians',sigma)
          call yaml_sequence_open('Values')
          do impl=1,ep%nmpl
              call yaml_sequence(advance='no')
              call yaml_mapping_open(trim(yaml_toa(impl)))
              call yaml_map('norm of the Gaussians',norm,fmt='(1es16.8)')
              call yaml_map('monopole',monopole,fmt='(1es16.8)')
              call yaml_map('dipole',dipole,fmt='(1es16.8)')
              call yaml_map('quadrupole',quadrupole,fmt='(1es16.8)')
              call yaml_mapping_close()
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()
      end if
      !!$omp end do
      !!$omp end parallel


      call H_potential('D',denspot%pkernel,density,denspot%V_ext,ehart_ps,0.0_dp,.false.,&
           quiet=denspot%PSquiet,rho_ion=denspot%rho_ion)
      pot = pot + density



      ii = 0
      do i3=is3,ie3
          ii3 = i3 - 15
          do i2=is2,ie2
              ii2 = i2 - 15
              do i1=is1,ie1
                  ii1 = i1 - 15
                  ii = ii + 1
                  write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, density(i1,i2,i3)
                  !do impl=1,ep%nmpl
                  !    r(1) = ep%mpl(impl)%rxyz(1) - x
                  !    r(2) = ep%mpl(impl)%rxyz(2) - y
                  !    r(3) = ep%mpl(impl)%rxyz(3) - z 
                  !    rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                  !    rnrm1 = sqrt(rnrm2)
                  !    tt = spherical_harmonic(l, m, x, y, z)*gaussian(sigma, rnrm1)
                  !    density(i1,i2,i3) =+ tt
                  !    !write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
                  !end do
              end do
          end do
      end do

      call f_free(density)

      call f_release_routine()


      contains


        function gaussian(sigma, r) result(g)
          use module_base, only: pi => pi_param
          implicit none
          ! Calling arguments
          real(kind=8),intent(in) :: sigma, r
          real(kind=8) :: g
          
          g = exp(-r**2/(2.d0*sigma**2))
          g = g/sqrt(2.d0*pi*sigma**2)**3
          !g = g/(sigma**3*sqrt(2.d0*pi)**3)
        end function gaussian

        !> Calculates the real spherical harmonic for given values of l, m, x, y, z.
        function spherical_harmonic(l, m, x, y, z) result(sh)
          use module_base, only: pi => pi_param
          implicit none
          ! Calling arguments
          integer,intent(in) :: l, m
          real(kind=8),intent(in) :: x, y, z
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
              sh = 0.5d0*sqrt(1/pi)
          case (1)
              r = sqrt(x**2+y**2+z**2)
              ! fix for small r (needs proper handling later...)
              if (r==0.d0) r=1.d-20
              select case (m)
              case (-1)
                  sh = sqrt(3.d0/(4.d0*pi))*y/r
              case (0)
                  sh = sqrt(3.d0/(4.d0*pi))*z/r
              case (1)
                  sh = sqrt(3.d0/(4.d0*pi))*x/r
              end select
          case (2)
              r2 = x**2+y**2+z**2
              ! fix for small r2 (needs proper handling later...)
              if (r2==0.d0) r2=1.d-20
              select case (m)
              case (-2)
                  sh = 0.5d0*sqrt(15.d0/pi)*x*y/r2
              case (-1)
                  sh = 0.5d0*sqrt(15.d0/pi)*y*z/r2
              case (0)
                  sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
              case (1)
                  sh = 0.5d0*sqrt(15.d0/pi)*z*x/r2
              case (2)
                  sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
              end select
          end select

        end function spherical_harmonic

    end subroutine potential_from_charge_multipoles


    !> Calculate the external potential arising from the multipoles provided
    subroutine potential_from_multipoles(ep, is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, shift, pot)
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3
      real(gp),intent(in) :: hx, hy, hz
      real(gp),dimension(3),intent(in) :: shift !< global shift of the atomic positions
      real(gp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, ii1, ii2, ii3, impl, l
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp
      real(dp),dimension(3) :: r

      !!$omp parallel &
      !!$omp default(none) &
      !!$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, ep, pot) &
      !!$omp private(i1, i2, i3, ii1, ii2, ii3, x, y, z, impl, r, rnrm1, rnrm2, rnrm3, rnrm5, l, mp)
      !!$omp do
      write(*,*) 'shift',shift
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
                  do impl=1,ep%nmpl
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
                      write(300+bigdft_mpi%iproc,*) 'i1, i2, i3, val', i1, i2, i3, mp
                  end do
              end do
          end do
      end do
      !!$omp end do
      !!$omp end parallel


      contains


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


        !function calc_octopole()
        !end function calc_octopole

    end subroutine potential_from_multipoles



    subroutine multipoles_from_density(iproc, nproc, at, lzd, smats, smatl, orbs, &
               npsidim, lphi, norbsPerType, collcom, collcom_sr, orthpar, &
               ovrlp, kernel, meth_overlap)
      use module_base
      use module_types
      use module_interfaces
      use sparsematrix_base, only: sparsematrix_malloc0, SPARSE_FULL, assignment(=)
      use sparsematrix_init, only: matrixindex_in_compressed
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim, meth_overlap
      type(atoms_data),intent(in) :: at
      type(local_zone_descriptors),intent(in) :: lzd
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(orbitals_data),intent(in) :: orbs
      real(kind=8),dimension(npsidim),intent(in) :: lphi
      integer,dimension(at%astruct%ntypes),intent(in) :: norbsPerType
      type(comms_linear),intent(in) :: collcom
      type(comms_linear),intent(in) :: collcom_sr
      type(orthon_data),intent(in) :: orthpar
      type(matrices),intent(in) :: kernel
      type(matrices),intent(inout) :: ovrlp !in principle also intent(in)

      ! Local variables
      integer :: ist, istr, iorb, iiorb, ilr, ii, natp, isat, nr, jproc, iat, n, norb_get, istr_get
      integer :: window, ioffset, methTransformOverlap, l, m, iiat, ityp, norb_per_atom, i1, i2, i3, ind, jorb, jat
      integer :: ii1, ii2, ii3, jjorb, i
      real(kind=8),dimension(:),allocatable :: psir
      real(kind=8),dimension(:),pointer :: phit_c, phit_f
      type(workarr_sumrho) :: w
      integer,dimension(:),allocatable :: nat_par, norb_list
      real(kind=8),dimension(:),allocatable :: psir_get, locrad, rmax
      real(kind=8),dimension(:,:),allocatable :: locregcenter, psir_get_fake
      real(kind=8),dimension(:,:,:,:),allocatable :: phi
      real(kind=8),dimension(:,:,:,:,:,:),allocatable :: sphi
      integer,dimension(:,:),allocatable :: comms
      logical :: can_use_transposed, arr_allocated
      real(kind=8) :: ddot, x, y, z, tt, rnorm, factor, max_error!, get_normalization, get_test_factor
      !real(kind=8) ,dimension(2,orbs%norb) :: testarr
      real(kind=8),dimension(:),allocatable :: kernel_ortho, phi_ortho
      integer,dimension(:),allocatable :: n1i, n2i, n3i, ns1i, ns2i, ns3i
      real(kind=8),dimension(-lmax:lmax,0:lmax,at%astruct%nat) :: multipoles
      real(kind=8) :: factor_normalization
      character(len=20) :: atomname
      real(kind=8),dimension(-lmax:lmax,0:lmax) :: norm
      !real(kind=8),parameter :: rmax=5.d0
      !testarr = 0.d0

      call f_routine(id='multipoles_from_density')

      ! Orthogonalize the support functions
      can_use_transposed = .false.
      methTransformOverlap = 1020
      phit_c = f_malloc_ptr(collcom%ndimind_c,id='phit_c')
      phit_f = f_malloc_ptr(7*collcom%ndimind_f,id='phit_f')
      phi_ortho = f_malloc(npsidim,id='phi_ortho')
      call vcopy(npsidim, lphi(1), 1, phi_ortho(1), 1)
      call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
           1.d-8, npsidim, orbs, lzd, &
           smats, smatl, collcom, orthpar, &
           phi_ortho, phit_c, phit_f, &
           can_use_transposed)
      call f_free_ptr(phit_c)
      call f_free_ptr(phit_f)


      ! Transform the support functions to real space
      psir = f_malloc(max(collcom_sr%ndimpsi_c,1),id='psir')
      ist=1
      istr=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call initialize_work_arrays_sumrho(1,lzd%Llr(ilr),.true.,w)
          call daub_to_isf(lzd%Llr(ilr), w, phi_ortho(ist), psir(istr))
          call deallocate_work_arrays_sumrho(w)
          !write(*,'(a,4i8,es16.6)') 'INITIAL: iproc, iiorb, n, istr, ddot', &
          !    iproc, iiorb, lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, &
          !    istr, ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1)
          !testarr(1,iiorb) = ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1) 
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      call f_free(phi_ortho)
      if(istr/=collcom_sr%ndimpsi_c+1) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
          stop
      end if


      ! Switch to a partition over the atoms
      nat_par = f_malloc(0.to.nproc-1,id='nat_par')
      ii = at%astruct%nat/nproc
      nat_par(0:nproc-1) = ii
      ii = at%astruct%nat-nproc*ii
      nat_par(0:ii-1) = nat_par(0:ii-1) + 1
      if (sum(nat_par(:))/=at%astruct%nat) then
          call f_err_throw('wrong partition of the atoms',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      natp = nat_par(iproc)
      if (iproc==0) then
          isat = 0
      else
          isat = sum(nat_par(0:iproc-1))
      end if
      call f_free(nat_par)
      comms = f_malloc((/4,orbs%norb/),id='comms')
      !write(*,'(a,i5,3x,13i4)') 'iproc, orbs%onwhichatom', iproc, orbs%onwhichatom
      nr = 0
      norb_get = 0
      istr = 1
      istr_get = 1
      do iat=isat+1,isat+natp
          do jproc=0,nproc-1
              istr = 0
              do iorb=1,orbs%norb_par(jproc,0)
                  iiorb = iorb + orbs%isorb_par(jproc)
                  ilr = orbs%inwhichlocreg(iiorb)
                  iiat = orbs%onwhichatom(iiorb)
                  n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
                  !write(*,'(a,5i8)') 'iproc, jproc, iiorb, iat', iproc, jproc, iiorb, iat
                  !if (iat>=isat+1 .and. iat<=isat+natp) then
                  if (iiat==iat) then
                      norb_get = norb_get + 1
                      comms(1,norb_get) = jproc
                      comms(2,norb_get) = n
                      comms(3,norb_get) = istr
                      comms(4,norb_get) = istr_get
                      nr = nr + n
                      !write(*,'(a,5i8)') 'iproc, jproc, n, istr, istr_get', iproc, jproc, n, istr, istr_get
                      istr = istr + n
                      istr_get = istr_get + n
                  else
                      istr = istr + n
                  end if
              end do
          end do
      end do
      !write(*,*) 'iproc, nr', iproc, nr
      !do iorb=1,norb_get
      !    write(*,'(a,2i7,4i9)') 'iproc, iorb, comm',iproc, iorb, comms(:,iorb)
      !end do
      psir_get = f_malloc(nr,id='psir_get')
      if (nproc>1) then
          window = mpiwindow(size(psir), psir(1), bigdft_mpi%mpi_comm)
          do iorb=1,norb_get
              jproc = comms(1,iorb)
              n = comms(2,iorb)
              ioffset = comms(3,iorb)
              istr = comms(4,iorb)
              !write(*,'(5(a,i0))') 'task ',iproc,' gets ',n,' elements at position ', &
              !                     istr,' from position ',ioffset+1,' on task ',jproc
              call mpiget(psir_get(istr), n, jproc, int(ioffset,kind=mpi_address_kind), window)
          end do
          call mpi_fenceandfree(window)
      else
          do iorb=1,norb_get
              n = comms(2,iorb)
              ioffset = comms(3,iorb)
              istr = comms(4,iorb)
              call vcopy(n, psir(ioffset+1), 1, psir_get(istr), 1)
          end do
      end if
      call f_free(psir)
      call f_free(comms)
      istr = 1
      !write(*,*) 'iproc, isat, natp', iproc, isat, natp
      do iorb=1,orbs%norb
          ilr = orbs%inwhichlocreg(iorb)
          iat = orbs%onwhichatom(iorb)
          if (iat>=isat+1 .and. iat<=isat+natp) then
              n = lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
              !write(*,'(a,4i8,es16.6)') 'iproc, iorb, n, istr, ddot', &
              !    iproc, iorb, n, istr, ddot(n, psir_get(istr), 1, psir_get(istr), 1)
              !testarr(2,iorb) = ddot(n, psir_get(istr), 1, psir_get(istr), 1)
              istr = istr + n
          end if
      end do
      !do iorb=1,size(psir_get)
      !    write(300+iproc,'(a,i7,es16.7)') 'i, psir_get(i)', iorb, psir_get(iorb)
      !end do
      !call mpiallred(testarr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      !if (iproc==0) then
      !    do iorb=1,orbs%norb
      !        write(*,*) 'DIFF, iorb, val', iorb, abs(testarr(1,iorb)-testarr(2,iorb))
      !    end do
      !end if


      ! Calculate the kernel for orthormal support functions
      kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='kernel_ortho')
      call kernel_for_orthonormal_basis(iproc, nproc, orbs%norbp, meth_overlap, smats, smatl, &
           ovrlp, kernel, kernel_ortho)


      ! Apply the spherical harmonics to the suport functions.
      ! For the support functions on atom A we only need to apply 
      ! the spherical harmonics centered as well on atom A.

      call unitary_test()


      !lmax = 1
      n1i = f_malloc(orbs%norb,id='n1i')
      n2i = f_malloc(orbs%norb,id='n2i')
      n3i = f_malloc(orbs%norb,id='n3i')
      ns1i = f_malloc(orbs%norb,id='ns1i')
      ns2i = f_malloc(orbs%norb,id='ns2i')
      ns3i = f_malloc(orbs%norb,id='ns3i')
      locrad = f_malloc(orbs%norb,id='locrad')
      locregcenter = f_malloc((/3,orbs%norb/),id='locregcenter')
      do iorb=1,orbs%norb
          n1i(iorb) = lzd%llr(iorb)%d%n1i
          n2i(iorb) = lzd%llr(iorb)%d%n2i
          n3i(iorb) = lzd%llr(iorb)%d%n3i
          ns1i(iorb) = lzd%llr(iorb)%nsi1
          ns2i(iorb) = lzd%llr(iorb)%nsi2
          ns3i(iorb) = lzd%llr(iorb)%nsi3
          locrad(iorb) = lzd%llr(iorb)%locrad
          locregcenter(1:3,iorb) = lzd%llr(iorb)%locregcenter(1:3)
      end do
      rmax = f_malloc(at%astruct%nat,id='rmax')
      call multipole_analysis_core(iproc, nproc, natp, isat, at%astruct%nat, at%astruct%ntypes, orbs%norb, &
           at%astruct%iatype, norbsPerType, orbs%inwhichlocreg, orbs%onwhichatom, &
           n1i, n2i, n3i, ns1i, ns2i, ns3i, locrad, lzd%hgrids, locregcenter, &
           nr, psir_get, psir_get, &
           smatl%nvctr, kernel_ortho, multipoles, rmax, 101, smatl)!, matrixindex)
      call f_free(psir_get_fake)
      call f_free(n1i)
      call f_free(n2i)
      call f_free(n3i)
      call f_free(ns1i)
      call f_free(ns2i)
      call f_free(ns3i)
      call f_free(locrad)
      call f_free(locregcenter)

      if (iproc==0) then
          call write_multipoles(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%atomnames, &
               multipoles, rmax, lzd%hgrids, without_normalization=.false.)
      end if
      call f_free(rmax)
      call f_free(psir_get)

      call f_free(kernel_ortho)

      call f_release_routine()



      contains


        subroutine unitary_test()
          implicit none
          real(kind=8),dimension(1) :: rmax
          integer,parameter :: n1i=101, n2i=81, n3i=91
          integer,parameter :: nsi1=0, nsi2=10, nsi3=20
          real(kind=8),dimension(3) :: locregcenter
          integer :: nr

          locregcenter(1) = (ceiling(real(n1i,kind=8)/2.d0)+nsi1-14-1)*0.5d0*lzd%hgrids(1)
          locregcenter(2) = (ceiling(real(n2i,kind=8)/2.d0)+nsi2-14-1)*0.5d0*lzd%hgrids(2)
          locregcenter(3) = (ceiling(real(n3i,kind=8)/2.d0)+nsi3-14-1)*0.5d0*lzd%hgrids(3)

          !psir_get_fake = f_malloc0((/lzd%llr(1)%d%n1i*lzd%llr(1)%d%n2i*lzd%llr(1)%d%n3i,2/),id='psir_get_fake')
          nr = n1i*n2i*n3i
          psir_get_fake = f_malloc0((/nr,2/),id='psir_get_fake')
          psir_get_fake(:,2) = 1.d0
          rmax(1) = min(n1i*0.25d0*lzd%hgrids(1), &
                        n2i*0.25d0*lzd%hgrids(2), &
                        n3i*0.25d0*lzd%hgrids(3))
          !rmax(1) = 5.d0
          !call yaml_map('rmax for unitary test',rmax(1))
          !write(*,'(a,6i6,3es16.6)') 'n1, n2, n3, ns1, ns2, ns3, locreg',lzd%llr(1)%d%n1i,lzd%llr(1)%d%n2i,lzd%llr(1)%d%n3i,lzd%llr(1)%nsi1,lzd%llr(1)%nsi2,lzd%llr(1)%nsi3, lzd%llr(1)%locregcenter(1), lzd%llr(1)%locregcenter(2), lzd%llr(1)%locregcenter(3)
          !write(*,'(a,6i6,3es16.6)') 'n1, n2, n3, ns1, ns2, ns3, locreg',n1i,n2i,n3i,nsi1,nsi2,nsi3, locregcenter(1), locregcenter(2), locregcenter(3)
          do i3=1,n3i!lzd%llr(1)%d%n3i
              !ii3 = lzd%llr(1)%nsi3 + i3 - 14 - 1
              ii3 = nsi3 + i3 - 14 - 1
              !z = ii3*0.5d0*lzd%hgrids(3) - lzd%llr(1)%locregcenter(3)
              z = ii3*0.5d0*lzd%hgrids(3) - locregcenter(3)
              !write(*,'(a,2i9,es12.4)') 'Z COORD: i3, ii3, z', i3, ii3, z
              do i2=1,n2i!lzd%llr(1)%d%n2i
                  !ii2 = lzd%llr(1)%nsi2 + i2 - 14 - 1
                  ii2 = nsi2 + i2 - 14 - 1
                  !y = ii2*0.5d0*lzd%hgrids(2) - lzd%llr(1)%locregcenter(2)
                  y = ii2*0.5d0*lzd%hgrids(2) - locregcenter(2)
                  !write(*,'(a,2i9,es12.4)') 'Y COORD: i2, ii2, y', i2, ii2, y
                  do i1=1,n1i!lzd%llr(1)%d%n1i
                      !ii1 = lzd%llr(1)%nsi1 + i1 - 14 - 1
                      ii1 = nsi1 + i1 - 14 - 1
                      !x = ii1*0.5d0*lzd%hgrids(1) - lzd%llr(1)%locregcenter(1)
                      x = ii1*0.5d0*lzd%hgrids(1) - locregcenter(1)
                      !write(*,'(a,2i9,es12.4)') 'X COORD: i1, ii1, x', i1, ii1, x
                      if (x**2+y**2+z**2>rmax(1)**2) cycle
                      !ind = (i3-1)*lzd%llr(1)%d%n2i*lzd%llr(1)%d%n1i + (i2-1)*lzd%llr(1)%d%n1i + i1
                      ind = (i3-1)*n2i*n1i + (i2-1)*n1i + i1
                      !write(*,'(a,3es12.4,i9)') 'x, y, z, ind', x, y, z, ind
                      do l=0,lmax
                          do m=-l,l
                              factor = get_test_factor(l,m)
                              select case (l)
                              case (0)
                                  psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor
                              case (1)
                                  select case (m)
                                  case (-1)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*y
                                  case ( 0)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*z
                                  case ( 1)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*x
                                  end select
                              case (2)
                                  select case (m)
                                  case (-2)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*x*y
                                  case (-1)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*y*z
                                  case ( 0)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*(-x**2-y**2+2*z**2)
                                  case ( 1)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*z*x
                                  case ( 2)
                                      psir_get_fake(ind,1) = psir_get_fake(ind,1) + factor*(x**2-y**2)
                                  end select
                              end select
                          end do
                      end do
                      !write(200,*) 'ind, val', ind, psir_get_fake(ind,1)
                  end do
              end do
          end do
          ! Only do it for one MPI task
          if (iproc==0) then
              !call multipole_analysis_core(0, 1, 1, 0, 1, 1, 1, &
              !     (/1/), (/1/), (/1/), (/1/), &
              !     (/lzd%llr(1)%d%n1i/), (/lzd%llr(1)%d%n2i/), (/lzd%llr(1)%d%n3i/), &
              !     (/lzd%llr(1)%nsi1/), (/lzd%llr(1)%nsi2/), (/lzd%llr(1)%nsi3/), rmax, &
              !     lzd%hgrids, lzd%llr(1)%locregcenter, &
              !     nr, psir_get_fake(:,1), psir_get_fake(:,2), &
              !     1, (/1.d0/), multipoles, rmax, 102, matrixindex=(/1/))
              call multipole_analysis_core(0, 1, 1, 0, 1, 1, 1, &
                   (/1/), (/1/), (/1/), (/1/), &
                   (/n1i/), (/n2i/), (/n3i/), &
                   (/nsi1/), (/nsi2/), (/nsi3/), rmax, &
                   lzd%hgrids, locregcenter, &
                   nr, psir_get_fake(:,1), psir_get_fake(:,2), &
                   1, (/1.d0/), multipoles, rmax, 102, matrixindex=(/1/))
          end if
          !call mpi_barrier(mpi_comm_world,i3)

          if (iproc==0) then
              call write_multipoles(1, 1, (/1/), (/'testatom'/), multipoles, rmax, lzd%hgrids, without_normalization=.true.)
          end if

        end subroutine unitary_test


    end subroutine multipoles_from_density



    !> Calculate S^1/2 * K * S^1/2, which is the kernel corresponding to a
    !! orthonormal set of supoprt functions.
    subroutine kernel_for_orthonormal_basis(iproc, nproc, norbp, meth_overlap, smats, smatl, &
               ovrlp, kernel, weight_matrix_compr)
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, SPARSE_TASKGROUP, &
                                   matrices_null, assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   deallocate_matrices
      use matrix_operations, only: overlapPowerGeneral
      use sparsematrix, only: matrix_matrix_mult_wrapper, gather_matrix_from_taskgroups
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

    end subroutine kernel_for_orthonormal_basis


    subroutine write_multipoles(nat, ntypes, iatype, atomnames, multipoles, rmax, hgrids, without_normalization)
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: nat, ntypes
      integer,dimension(nat),intent(in) :: iatype
      character(len=*),dimension(ntypes),intent(in) :: atomnames
      real(kind=8),dimension(-lmax:lmax,0:lmax,nat),intent(in) :: multipoles
      real(kind=8),dimension(nat),intent(in) :: rmax
      real(kind=8),dimension(3),intent(in) :: hgrids
      logical,intent(in) :: without_normalization
      
      ! Local variables
      character(len=20) :: atomname
      integer :: i, iat, l, m, nit
      real(kind=8) :: max_error, factor!, get_normalization, get_test_factor
      real(kind=8),dimension(:,:,:),allocatable :: multipoles_tmp

          multipoles_tmp = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.nat/),id='multipoles_tmp')

          if (without_normalization) then
              nit = 2
          else
              nit = 1
          end if

          factor = 0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3)

          max_error = 0.d0
          call yaml_mapping_open('Multipole coefficients')
          do i=1,nit
              if (i==1) then
                  call yaml_map('normalized',.true.)
                  call yaml_map('radius of normalization sphere',(/minval(rmax),maxval(rmax)/))
                  call f_memcpy(src=multipoles, dest=multipoles_tmp)
              else if (i==2) then
                  call yaml_map('normalized',.false.)
                  do iat=1,nat
                      do l=0,lmax
                          do m=-l,l
                              !multipoles_tmp(m,l,iat) = multipoles(m,l,iat)/((get_normalization(rmax, l, m)*0.821583836)**2)
                              !multipoles_tmp(m,l,iat) = multipoles(m,l,iat)*((get_normalization(rmax, l, m)*0.106726871))
                              multipoles_tmp(m,l,iat) = multipoles(m,l,iat)*get_normalization(rmax(iat),l,m)**2*factor
                              max_error = max(max_error,abs(multipoles_tmp(m,l,iat)-get_test_factor(l,m)))
                              !write(*,'(a,3i5,2es14.5)') 'iat, l, m, multipoles(m,l,iat), ref', iat, l, m, multipoles(m,l,iat), get_test_factor(l,m)
                          end do
                      end do
                  end do
              end if
              call yaml_sequence_open('Values')
              do iat=1,nat
                  call yaml_sequence(advance='no')
                  atomname=atomnames(iatype(iat))
                  call yaml_sequence_open(trim(atomname))
                  do l=0,lmax
                      call yaml_sequence(advance='no')
                      !call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat),fmt='(1es16.8)')
                      !call yaml_map('l='//yaml_toa(l),multipoles(-l:l,l,iat)*sqrt(4.d0**(2*l+3)),fmt='(1es16.8)')
                      !do m=-l,l
                          !multipoles(m,l,iat) = multipoles(m,l,iat)*get_normalization(rmax, l, m)
                          !max_error = max(max_error,abs(multipoles(m,l,iat)-get_test_factor(l,m)))
                      !end do
                      call yaml_map('l='//yaml_toa(l),multipoles_tmp(-l:l,l,iat),fmt='(1es16.8)')
                      call yaml_newline()
                  end do
                  !call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
                  call yaml_sequence_close()
              end do
              call yaml_sequence_close()
              if (i==2) then
                  call yaml_map('max error from original values',max_error)
              end if
          end do
          call yaml_mapping_close()


          call f_free(multipoles_tmp)

    end subroutine write_multipoles


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


    !> Calculates the real spherical harmonic for given values of l, m, x, y, z.
    !! The functions are normalized to one within a sphere of radius rmax.
    function spherical_harmonic(rmax, l, m, x, y, z) result(sh)
      use module_base, only: pi => pi_param
      implicit none
      ! Calling arguments
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
          rnorm = sqrt(rmax**3/3.d0)
          !sh = sqrt(4.d0*pi)*0.5d0*sqrt(1/pi)
          sh = 0.5d0*sqrt(1/pi)/rnorm
      case (1)
          rnorm = sqrt(rmax**5/5.d0)
          r = sqrt(x**2+y**2+z**2)
          r = 1.d0
          ! fix for small r (needs proper handling later...)
          if (r==0.d0) r=1.d-20
          select case (m)
          case (-1)
              !sh = sqrt(4*pi/3.d0)*sqrt(3.d0/(4.d0*pi))*y/r
              !sh = sqrt(3.d0/(4.d0*pi))*y/r
              sh = sqrt(3.d0/(4.d0*pi))*y/rnorm
          case (0)
              !sh = sqrt(4*pi/3.d0)*sqrt(3.d0/(4.d0*pi))*z/r
              !sh = sqrt(3.d0/(4.d0*pi))*z/r
              sh = sqrt(3.d0/(4.d0*pi))*z/rnorm
          case (1)
              !sh = sqrt(4*pi/3.d0)*sqrt(3.d0/(4.d0*pi))*x/r
              !sh = sqrt(3.d0/(4.d0*pi))*x/r
              sh = sqrt(3.d0/(4.d0*pi))*x/rnorm
          end select
      case (2)
          rnorm = sqrt(rmax**7/7.d0)
          r2 = x**2+y**2+z**2
          r2=1.d0
          ! fix for small r2 (needs proper handling later...)
          if (r2==0.d0) r2=1.d-20
          select case (m)
          case (-2)
              !sh = sqrt(4.d0*pi/5.d0)*0.5d0*sqrt(15.d0/pi)*x*y/r2
              !sh = 0.5d0*sqrt(15.d0/pi)*x*y/r2
              sh = 0.5d0*sqrt(15.d0/pi)*x*y/rnorm
          case (-1)
              !sh = sqrt(4.d0*pi/5.d0)*0.5d0*sqrt(15.d0/pi)*y*z/r2
              !sh = 0.5d0*sqrt(15.d0/pi)*y*z/r2
              sh = 0.5d0*sqrt(15.d0/pi)*y*z/rnorm
          case (0)
              !sh = sqrt(4.d0*pi/5.d0)*0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
              !sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
              sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/rnorm
          case (1)
              !sh = sqrt(4.d0*pi/5.d0)*0.5d0*sqrt(15.d0/pi)*z*x/r2
              !sh = 0.5d0*sqrt(15.d0/pi)*z*x/r2
              sh = 0.5d0*sqrt(15.d0/pi)*z*x/rnorm
          case (2)
              !sh = sqrt(4.d0*pi/5.d0)*0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
              !sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
              sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/rnorm
          end select
      end select

    end function spherical_harmonic


    subroutine multipole_analysis_core(iproc, nproc, natp, isat, nat, ntypes, norb, &
               iatype, norbsPerType, inwhichlocreg, onwhichatom, &
               n1i, n2i, n3i, nsi1, nsi2, nsi3, locrad, hgrids, locregcenter, &
               nr, psir1_get, psir2_get, &
               nvctr, matrix_compr, multipoles, rmax, get_index, smatl, matrixindex)
      use module_types, only: workarr_sumrho
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      !real(kind=8),intent(in) :: rmax
      integer,intent(in) :: natp, isat, nat, ntypes, norb, nr, nvctr, get_index
      integer,dimension(nat),intent(in) :: iatype
      integer,dimension(ntypes),intent(in) :: norbsPerType
      integer,dimension(norb),intent(in) :: inwhichlocreg, onwhichatom
      integer,dimension(norb),intent(in) :: n1i, n2i, n3i, nsi1, nsi2, nsi3
      real(kind=8),dimension(norb),intent(in) :: locrad
      real(kind=8),dimension(3),intent(in) :: hgrids
      real(kind=8),dimension(3,norb) :: locregcenter
      real(kind=8),dimension(nr),intent(in) :: psir1_get, psir2_get
      real(kind=8),dimension(nvctr),intent(in) :: matrix_compr
      real(kind=8),dimension(-lmax:lmax,0:lmax,nat),intent(out) :: multipoles
      real(kind=8),dimension(nat),intent(out) :: rmax
      type(sparse_matrix),intent(in),optional :: smatl
      integer,dimension(norb,norb),intent(in),optional :: matrixindex

      ! Local variables
      integer,parameter :: INDEX_AUTOMATIC=101
      integer,parameter :: INDEX_MANUALLY=102
      integer :: ist, istr, iorb, iiorb, ilr, ii, jproc, iat, n, norb_get, istr_get
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
      !real(kind=8) :: rmax

      !! Check that rmax does remains within the box.
      !if (rmax>=0.5d0*(0.5d0*hgrids(1)*maxval(n1i)+0.5d0*hgrids(2)*maxval(n2i)+0.5d0*hgrids(3)*maxval(n3i))) then
      !    call f_err_throw('The radius for the multipole analysis is too small', err_name='BIGDFT_RUNTIME_ERROR')
      !end if


      if (iproc==0) call yaml_comment('Multipole analysis',hfill='-')

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
                  !write(*,'(a,6i6)') 'iat, iiat, ilr, n1i, n2i, n3i', iat, iiat, ilr, n1i(ilr), n2i(ilr), n3i(ilr)
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
                              !write(300,*) 'ind, val', ind, phi1(i1,i2,i3,iiorb)
                              !write(400,*) 'ind, val', ind, (1.d0 + 11.d0*y + 12.d0*z + 13.d0*x + &
                              !    21.d0*x*y + 22.d0*y*z + 23.d0*(-x**2-y**2+2*z**2) + 24.d0*z*x + 25.d0*(x**2-y**2))
                              !write(*,*) 'DIFF', &
                              !    phi1(i1,i2,i3,iiorb) - &
                              !    (1.d0 + 11.d0*y + 12.d0*z + 13.d0*x + &
                              !    21.d0*x*y + 22.d0*y*z + 23.d0*(-x**2-y**2+2*z**2) + 24.d0*z*x + 25.d0*(x**2-y**2))
                              !phi1(i1,i2,i3,iiorb) = 1.d0 + 11.d0*y + 12.d0*z + 13.d0*x + 21.d0*x*y + 22.d0*y*z + 23.d0*(-x**2-y**2+2*z**2) + 24.d0*z*x + 25.d0*(x**2-y**2)
                              !write(200+iproc,'(a,4i7,es16.7)') 'iiat, iorb, ist, ind, iorb, psir2_get(ist+ind)',&
                              !    iiat, iorb, ist, ind, psir2_get(ist+ind)
                              do l=0,lmax
                                  do m=-l,l
                                      tt = spherical_harmonic(rmax(iiat), l, m, x, y, z)
                                      !tt = factor_normalization*tt
                                      !sphi(i1,i2,i3,m,l,iiorb) = factor_normalization*tt*phi1(i1,i2,i3,iiorb)
                                      !norm(m,l) = norm(m,l) + (factor_normalization*tt)**2
                                      sphi(i1,i2,i3,m,l,iiorb) = tt*phi1(i1,i2,i3,iiorb)
                                      !sphi(i1,i2,i3,m,l,iiorb) = tt**2*phi1(i1,i2,i3,iiorb)
                                      !sphi(i1,i2,i3,m,l,iiorb) = phi1(i1,i2,i3,iiorb)
                                      !norm(m,l) = norm(m,l) + tt**2
                                      norm(m,l) = norm(m,l) + (factor_normalization*tt)**2
                                      !if (i1==1 .and. i2==1) write(*,'(a,2i5,3es13.3,3es19.8)') 'l, m, x, y, z, tt, phi, sphi', l, m, x, y, z, tt, phi(i1,i2,i3,iiorb), sphi(i1,i2,i3,iiorb,m,l)
                                  end do
                              end do
                          end do
                      end do
                  end do
                  do l=0,lmax
                      do m=-l,l
                          !norm(m,l) = norm(m,l) + (factor_normalization*tt)**2
                          rnorm_maxdev = max(rnorm_maxdev,abs(1.d0-norm(m,l)))
                      end do
                  end do
                  !call yaml_map('rnorm',rnorm)
                  !call yaml_map('norm',norm,fmt='(es14.6)')
                  !call yaml_map('4.d0/rnorm',4.d0/rnorm)
                  !sphi(:,:,:,:,:,iiorb) = sphi(:,:,:,:,:,iiorb)*sqrt(0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3))
                  !write(*,*) 'calling vscal,iproc, iat, iiorb', iproc, iat, iiorb
                  !call vscal(n1i(ilr)*n2i(ilr)*n3i(ilr)*(2*lmax+1)*(lmax+1), &
                  !     sqrt(0.5d0*hgrids(1)*0.5d0*hgrids(2)*0.5d0*hgrids(3)), sphi(1,1,1,-lmax,0,iiorb), 1)
                  !write(*,*) 'after vscal, iproc', iproc
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
                          !tt = ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), sphi(1,1,1,iorb,m,l), 1, phi2(1,1,1,jorb), 1)
                          !write(*,'(a,9i9)') 'call ddot, iproc, iat, iorb, jorb, n, size, n1, n2, n3', &
                          !    iproc, iat, iorb, jorb, n1i(iilr)*n2i(iilr)*n3i(iilr), size(phi2,1)*size(phi2,2)*size(phi2,3), &
                          !    n1i(iilr), n2i(iilr), n3i(iilr)
                          tt = ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), sphi(1,1,1,m,l,iorb), 1, phi2(1,1,1,jorb), 1)
                          !tt = 0.d0
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
                          !write(*,'(a,9i7,5es18.8)') 'iproc, iiat, m, l, iorb, iiorb, jorb, jjorb, ind, tt, K, sphi, phi1, phi2', &
                          !     iproc, iiat, m, l, iorb, iiorb, jorb, jjorb, ind, tt, matrix_compr(ind), &
                          !     ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), sphi(1,1,1,iorb,m,l), 1, sphi(1,1,1,jorb,m,l), 1), &
                          !     ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), phi1(1,1,1,iorb), 1, phi1(1,1,1,jorb), 1), &
                          !     ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), phi2(1,1,1,iorb), 1, phi2(1,1,1,jorb), 1)
                          !tt = tt/((get_normalization(rmax,0,0))**1)
                          tt = tt/((get_normalization(rmax(iiat),l,m))**1)
                          multipoles(m,l,iiat) = multipoles(m,l,iiat) + matrix_compr(ind)*tt
                          tt2 = ddot(n1i(iilr)*n2i(iilr)*n3i(iilr), phi1(1,1,1,iorb), 1, phi2(1,1,1,jorb), 1)
                          !tt = tt*real(ii,kind=8)
                          !write(*,'(a,5i8,2es16.8)') 'iproc, l, m, iorb, jorb, ddots', &
                          !    iproc, l, m, iorb, jorb, tt, tt2
                              !iproc, l, m, iorb, jorb, tt*real(ii**2,kind=8)/real(n1i(iilr)*n2i(iilr)*n3i(iilr),kind=8), tt2
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
      if (iproc==0) then
          call yaml_mapping_open('Calculation of solid harmonics')
          call yaml_map('radius for normalization',(/minval(rmax),maxval(rmax)/),fmt='(es8.2)')
          call yaml_map('max deviation from normalization',rnorm_maxdev,fmt='(es8.2)')
          call yaml_mapping_close()
      end if


      call f_free(norb_list)
    end subroutine multipole_analysis_core

end module multipole
