!> @file 
!! @brief Atomic program for pseudopotential calculations
!! @author
!!    Program for atomic calculations
!!    written by Sverre Froyen, February 1982
!!    while at UC Berkeley, California
!!
!!    then modified by
!!    Christian Hartwigsen, February 1998
!!    while at MPI Stuttgart, Germany
!!
!!    and further altered by
!!    Alex Willand, under the supervision of
!!    Stefan Goedecker, December 2010
!!    while at Universitaet Basel, Switzerland
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!
!! @warning
!!    Some parameters are set inside the program
!!    the most important ones are the tolerance for the self-
!!    consistency in the potential (set in the main program)
!!    and the accuracy in the eigenvalue for a given potential
!!    (set in difrel,difnrl)
!!
!! @todo
!!    Possibility to use hybrid functionals (need exchange interaction)
!!    atom is not parallelized: The tests are a problem...


!> Calculate the all-electron electronic structure for one atom
!! @ingroup PSEUDO
program ae_atom

   implicit none

   !Parameters
   integer, parameter :: nrmax=10000                !< Maximal radial grid points
   integer, parameter :: maxorb=60                  !< Maximal number of orbitals
   integer, parameter :: lmax=5                     !< Maximal orbital moment
   integer, parameter :: maxconf=19                 !< Maximal electronic configurations
   integer, parameter :: maxit_max=3000             !< Maximal number of iterations
   character(len=2), parameter :: stop_chain = 'st'
   real(kind=8), parameter :: tol=1.0d-11
   real(kind=8), parameter :: xmixo_min = 1.d-5     !< Minimal value of the mixing parameter
   logical, parameter :: debug=.false.              !< Debug flag
   !Local variables
   integer :: norb                                  !< Number of orbitals
   integer, dimension(maxorb) :: no                 !< For each orbital, n quantum number
   integer, dimension(maxorb) :: lo                 !< For each orbital, l quantum number
   real(kind=8), dimension(maxorb) :: so            !< spin (+/- 0.5, or 0 for unpolarized)
   real(kind=8), dimension(maxorb) :: zo            !< electrons
                                                   
   integer :: ncore                                 !< Number of core electrons
   real(kind=8) :: znuc                             !< Atomic number
                                                   
   integer :: nr                                    !< Mesh points
   real(kind=8) :: a,b                              !< Parameters to build the logarithmic mesh
   real(kind=8), dimension(nrmax) :: r              !< Radial mesh r(i) = a*(exp(b*(i-1))-1)
   real(kind=8), dimension(nrmax) :: rab            !< rab(i) = (r(i)+a)*b (integration grid)
                                                    !< c.hartwig: additional grids for modified integration
   real(kind=8), dimension(nrmax) :: rw             !< rw(i) = rab(i)*12.56637061435917d0*r(i)**2 (integration grid)
   real(kind=8), dimension(nrmax) :: rd             !< rd(i) = 1/rab(i) (integration grid)
   real(kind=8), dimension(nrmax) :: cdd            !< Charge density (spin down, 2 pi r**2 rho(r))
   real(kind=8), dimension(nrmax) :: cdu            !< Charge density (spin up,   2 pi r**2 rho(r))
   real(kind=8), dimension(nrmax) :: cdc            !< Core charge density (up to ncore orbitals)
   ! Potentials: always multiplied by r.
   real(kind=8), dimension(lmax,nrmax) :: viod,viou !< viod,u  ionic potential (down,up)
   real(kind=8), dimension(nrmax) :: vid,viu        !< vid,u   input screening potential (down,up)
   real(kind=8), dimension(nrmax) :: vod,vou        !< vod,u   output screening potential (down,up)

   !> etot(i) i=1,10 contains various contributions to the total energy.
   !!     (1)   Sum of eigenvalues ev
   !!     (2)   Sum of orbital kinetic energies ek
   !!     (3)   El-ion interaction from sum of orbital, potential energies ep
   !!     (4)   Electrostatic el-el interaction (from velect)
   !!     (5)   Vxc (exchange-correlation) correction to sum of eigenvalues (from velect)
   !!     (6)   3 * vc - 4 * ec, correction term for virial theorem, when correlation is included (from velect)
   !!     (7)   Exchange and correlation energy  (from velect)
   !!     (8)   Kinetic energy from eigenvalues  (1,3,4,5)
   !!     (9)   Potential energy
   !!     (10)  Total energy
   real(kind=8), dimension(10) :: etot

   real(kind=8), dimension(maxconf) :: econf         !< Total energy of the different electronic configurations
   real(kind=8), dimension(maxorb) :: ev             !< Eigenvalues
   real(kind=8), dimension(maxorb) :: ek             !< Kinetic energy for each orbital
   real(kind=8), dimension(maxorb) :: ep             !< Potential energy (Vionic*rho) for each orbital

   character(len=2) :: nameat,itype                  !< Name and type of atom
   character(len=1) :: ispp                          !< Spin, relativistic calculation ('n', 's', 'r')

   integer :: iXC                                    !< Exchange-Correlation parameter
   integer :: iter                                   !< Iteration number
   integer :: iconv                                  !< Convergence done or not
   integer :: nconf                                  !< Number of electronic configurations
   integer :: nspol                                  !< number of spin components (1 or 2 spin polarisation)
   integer :: ifcore                                 !< If core electrons (for NLCC ?)
   integer :: nvalo,ncoreo
   integer :: maxit                                  !< Maximal number of iterations
   real(kind=8) :: rcov                              !< Covalent radius
   real(kind=8) :: rprb                              !< Radius for the parabolic confinement potential
   real(kind=8) :: zcore,zel

   character(len=2) :: name_old,itype_old,cnum
   integer :: iXC_old
   logical :: abort
   real(kind=8) :: dvold,aa,a2,dcrc,ddcrc,dv,dvmax,weight,xmixo
   integer :: i,icon2,ii,iorb

   !Heuristic value for fluctuating GGAs
   name_old = '  '
   iXC_old = 0
   itype_old ='  '
   nconf = 0
   dvold = 1.0d10
   nr    = 1
   norb  = 1

   !Open the main input file: the routine input will read each configuration
   open(unit=35,file='atom.dat',status='unknown')

   !begin main loop over electronic configurations
   loop_configurations: do

      ! Read input data
      !r ...... radial mesh
      !nr ..... # mesh points
      !norb ... # orbitals
      !ncore .. # core orbitals (closed shells)
      !no ..... n quantum number
      !lo ..... l do.
      !so ..... spin (+/- 0.5, or 0 for unpolarized)
      !zo ..... # electrons
      !znuc ... atomic number
      call input(itype,iXC,ispp, &
      & nrmax,nr,a,b,r,rab,rw,rd,rprb,rcov, &
      & nameat,norb,ncore, &
      & maxorb,no,lo,so,zo, &
      & znuc,zel,zcore, &
      & nconf,  &
      & nvalo,ncoreo) 

      if (itype == stop_chain) exit loop_configurations

      if (nconf > maxconf) then
         write(6,'(a,1x,i0)') 'too many configurations, max. is:',maxconf
         stop
      end if

      !set up initial charge density.
      !cdd  charge density (spin down)
      !cdu  charge density (spin up)
      !cdc  core charge density (up to ncore orbitals)
      !cdd and cdu  =  2 pi r**2 rho(r)
      aa = sqrt(sqrt(znuc))/2.0d0+1.0d0
      a2 = zel/4.0d0*aa**3
      do i=1,nr
        cdd(i) = a2*exp(-aa*r(i))*r(i)**2
        cdu(i) = cdd(i)
      end do

      !set up ionic potentials

      call vionic(ifcore, &
         nr,r,rprb,lmax, &
         znuc, &
         viod,viou)

      !Potentials: always multiplied by r.
      !viod,u ..... ionic potential (down,up)
      !vid,u ...... input screening potential (down,up)
      !vod,u ...... output screening potential (down,up)

      !set up electronic potential

      !new variable nspol: spin channels for XC
      nspol=1
      if (ispp=='s') nspol=2

      call velect(0,0,iXC,nspol,ifcore,  &
         nr,r,rab,rw,rd, &
         zel,cdd,cdu,cdc,  &
         vod,vou,  &
         etot)

      do i=1,nr
         vid(i) = vod(i)
         viu(i) = vou(i)
         if (debug) write(*,'(a,i6,2(1pe22.15))') 'DEBUG: vid viu',i,vid(i),viu(i)
      end do

      !start iteration loop

      iconv = 0
      icon2 = 0
      maxit = maxit_max
      if (debug) then
         write(6,'(a)') 'DEBUG: enter max SCF iterations'
         read(*,*) maxit
      end if


      !empirical function
      xmixo = 1.0d0/log(znuc+7.0d0)

      !start of iteration loop
      do iter=1,maxit

         if (iter == maxit) iconv=1

         !compute orbitals (solve Schrodinger equation)
         if (icon2 == 0) then
            !finite difference solution (less accurate)
            call dsolv1( &
               nr,a,b,r,rab,lmax, &
               norb,no,lo,so,zo, &
               cdd,cdu, &
               viod,viou,vid,viu, &
               ev)

         else
            !predictor - corrector method (more accurate)
            call dsolv2(iter,iconv,iXC,ispp,ifcore, &
               nr,a,b,r,rab,lmax, &
               norb,ncore,no,lo,so,zo, &
               znuc,zcore,cdd,cdu,cdc,dcrc,ddcrc, &
               viod,viou,vid,viu, &
               ev,ek,ep,rcov,rprb,nconf)

         end if

         !etot ..... terms in Etotal
         !ev ....... eigenvalues
         !ek ....... kinetic energy for each orbital
         !ep ....... potential energy (Vionic*rho) for each orbital

         !set up output electronic potential from charge density
         call velect(iter,iconv,iXC,nspol,ifcore,  &
            nr,r,rab,rw,rd, &
            zel,cdd,cdu,cdc,  &
            vod,vou,  &
            etot)

         !check for convergence (Vout - Vin)
         if (iconv > 0) exit

         dvmax = 0.d0
         do i=2,nr
            dv = (vod(i)-vid(i))/(1.d0+vod(i)+vou(i))
            if (abs(dv) > dvmax) dvmax=abs(dv)
            dv = (vou(i)-viu(i))/(1.d0+vou(i)+vod(i))
            if (abs(dv) > dvmax) dvmax=abs(dv)
         end do
         iconv = 1
         icon2 = icon2+1
         if (dvmax >  tol) iconv=0
         if (dvmax >= dvold) xmixo=0.8d0*xmixo
         inquire(file='EXIT', exist=abort)
         if (abort) iconv=1
         !EXPERIMENTAL: why not in both directions?
         !if (dvmax <= dvold) xmixo=1.05d0*xmixo

         !For now, ignore convergence for at least the first 30 cycles
         !because we may want to switch to GGA thereafter
         if (iter < 40) iconv=0

         ! diverging - reduce mixing coefficient
         if (xmixo < xmixo_min) xmixo = xmixo_min
         dvold = dvmax
         write(6,'(1x,a,i5,1x,a,1pe19.12,1x,a,1pe19.12)') 'iter =', iter, 'dvmax =', dvmax, 'xmixo =', xmixo

         ! mix input and output electronic potentials
         call mixer(xmixo,nr,vid,viu,vod,vou)


      !end of iteration loop
      end do

      if (iconv == 0) then
         write(6,'(/," potential not converged - dvmax =",1pe11.4,"  xmixo = ",0pf6.3)') dvmax,xmixo
         call ext(1)
      end if

      !Find total energy
      call etotal(nameat,norb,no,lo,so,zo, etot,ev,ek,ep)

      if (name_old /= nameat .or. iXC_old /= iXC .or. itype_old /= itype ) call prdiff(nconf,econf)

      nconf = nconf + 1
      econf(nconf) = etot(10)
      if (nconf /= 1) write(6,'(//," excitation energy         =",f18.8,/,1x,45("-"))') etot(10)-econf(1)
      name_old = nameat
      iXC_old = iXC
      itype_old = itype


      if (nconf == 1) then
         ! it is better to append and not overwrite existing weights
         open(unit=60,file='input.weights',position='append')

         write(60,'(a,i3,a)') '----suggested weights for occup numbers of conf',nconf,'-----'
         write(60,'(/,a)') ' 1d0 1d5 1d0 1d0 1d0 1d3  psi(0), dEkin_wvlt, radii, hij, locality, E_exct'
        !     we put the overall weights for each configuration in atom.??.ae files
        !     such that the user can combine atomic data more easily.
        !     the old format was:
        !     if (nconf>1) then
        !        write(60,*) ('1.00 ',i=1,nconf),' weights for configurations'
        !        write(60,*) '0.00 ',('1.00 ',i=2,nconf),
        !    :        ' weights for excitation-energies '
        !     end if
         write(60,'(a)') '  n   l  so    eigval  chrg    dchrg  ddchrg  res    rnode dnode ddnode'
         do iorb=ncore+1,norb
            weight=0.0d0
            if (zo(iorb) > 1.0d-4) then
               write(60,'(2i4,1x,f5.2,tr3,a)') no(iorb),lo(iorb),so(iorb),  &
               '1.0e5   1.0e5   0.0e0  0.0e0   1.0e5  1.0e0 0.0e0 0.0e0'
            else
               write(60,'(2i4,1x,f5.2,tr3,a)') no(iorb),lo(iorb),so(iorb),  &
               '1.0e0   1.0e0   0.0e0  0.0e0   0.0e0  0.0e0 0.0e0 0.0e0'
            end if
         end do
         close(unit=60)
      end if

      !next electronic configuration of the atom
   end do loop_configurations !end loop over electronic configurations


   !All electronic configurations are calculated.

   !FITPAR, do not overwrite, append
   open(unit=60,file='input.fitpar',position='append')
   write(60,'(a)') ' fitting parameters appended by atom.f90: auto'
   close(unit=60)

   !and input.pseudo
   open(unit=60,file='input.pseudo',position='append')
   write(60,'(a)') 'input line written by atom: -plot -ng 20 -rij 2.0'
   close(unit=60)

   !append excitation energies (in hartree!) to file atom.ae
   !if more than one configuration
   if (nconf > 1) then
      do ii=0,nconf-1
         close(unit=40)
         write(cnum,'(i2.2)') ii
         open(unit=40, file='atom.'//cnum//'.ae', position='append')
         !write(40,*) 'Configuration ',nconf
         !write(40,*) 'Total reference energy',econf(1)
         write(40,'(f21.14,6x,a)') econf(ii+1),'total energy of this configuration'
         !write(40,*) 'Total energy of this configuration',econf(ii+1)
         write(40,'(a)') 'EXCITATION ENERGIES:'
         do i=1,nconf
            write(40,'(i2.2,1x,1pe25.17)') i,(econf(i)-econf(1))/2.d0
         end do
         close(unit=40)
      end do
   end if

   !call libxc_functionals_end()
   call prdiff(nconf,econf)
   call ext(0)

end program ae_atom


!> Print difference total energy between electronic configurations
subroutine prdiff(nconf,econf)
   implicit none
   !Arguments
   integer, intent(inout) :: nconf
   real(kind=8), dimension (nconf) :: econf
   !Local variables
   integer :: i,j
   if (nconf <= 1) then
      nconf = 0
      return
   end if
   write(6,'(/,a)') '---------------------------------------------'
   write(6,'(" Total energy differences",//,2x,19i9)') (i,i=0,nconf-1)
   do i=1,nconf
      write(6,'(1x,i2,1x,19f9.5)') i-1,(0.5d0*(econf(i)-econf(j)),j=1,i)
   end do
end subroutine prdiff


!> Subroutine computes the new exchange correlation potential
!! given the input and the output potential from the previous
!! iteration.
!! @todo
!! Improving the mixing algorithm
subroutine mixer(xmixo,nr,vid,viu,vod,vou)
   implicit none
   !Arguments
   integer, intent(in) :: nr
   real(kind=8), intent(in) :: xmixo
   real(kind=8), dimension(nr), intent(in) :: vod,vou
   real(kind=8), dimension(nr), intent(inout) :: vid,viu
   !Local variables
   real(kind=8) :: xmixi
   integer :: i

   xmixi = 1.d0 - xmixo
   do i=1,nr
      vid(i) = xmixo * vod(i) + xmixi * vid(i)
      viu(i) = xmixo * vou(i) + xmixi * viu(i)
   end do
end subroutine mixer


!> etotal computes the total energy from the electron charge density.
subroutine etotal(nameat,norb,no,lo,so,zo, etot,ev,ek,ep)
   implicit none
   !Arguments
   character(len=2), intent(in) :: nameat                    !< Name of the atom
   integer, intent(in) :: norb                               !< #orbitals
   integer, dimension(norb), intent(in) :: no, lo            !< quantum numbers
   real(kind=8), dimension(norb), intent(in) :: so, zo, ek   !< used only to display the values
   real(kind=8), dimension(norb), intent(out) :: ev, ep
   !> etot(i) i=1,10 contains various contributions to the total energy.
   !!     (1)   Sum of eigenvalues ev
   !!     (2)   Sum of orbital kinetic energies ek
   !!     (3)   El-ion interaction from sum of orbital, potential energies ep
   !!     (4)   Electrostatic el-el interaction (from velect)
   !!     (5)   Vxc (exchange-correlation) correction to sum of eigenvalues (from velect)
   !!     (6)   3 * vc - 4 * ec, correction term for virial theorem, when correlation is included (from velect)
   !!     (7)   Exchange and correlation energy  (from velect)
   !!     (8)   Kinetic energy from eigenvalues  (1,3,4,5)
   !!     (9)   Potential energy
   !!     (10)  Total energy
   real(kind=8), dimension(10), intent(out) :: etot
   !Local variables
   character(len=1), dimension(5) :: il
   real(kind=8) :: esh,vshift
   integer :: i

   ! sum up eigenvalues ev, kinetic energies ek, and
   ! el-ion interaction ep
   etot(1) = 0.d0
   etot(2) = 0.d0
   etot(3) = 0.d0
   !c.hartwig
   !subtract vshift
   vshift=-15.0d0
   do i=1,norb
      etot(1) = etot(1) + zo(i)*(ev(i)-vshift)
      etot(2) = etot(2) + zo(i)*ek(i)
      etot(3) = etot(3) + zo(i)*(ep(i)-vshift)
      !etot(1) = etot(1) + zo(i)*ev(i)
      !etot(2) = etot(2) + zo(i)*ek(i)
      !etot(3) = etot(3) + zo(i)*ep(i)
   end do

   !compute interaction shell - (nucleus-core)
   esh = 0.d0

   !kinetic energy
   etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)

   !potential energy
   etot(9) = etot(3) + etot(4) + etot(7) + esh

   !total energy
   etot(10) = etot(1) - etot(4) - etot(5) + etot(7) + esh

   !printout
   il(1) = 's'
   il(2) = 'p'
   il(3) = 'd'
   il(4) = 'f'
   il(5) = 'g'
   write(6,'(/,a3,1x,a,/,1x,27("-"),//,1x,a,9x,a,3x,a,6x,a,/)') &
        nameat,'output data for orbitals','nl    s      occ','eigenvalue','kinetic energy','pot energy'
   do i=1,norb
      !c.hartwig give energies in hartree
      ev(i) = ev(i) - vshift
      ep(i) = ep(i) - vshift
      write(6,'(1x,i1,a1,f6.1,f10.4,3f17.8)') no(i),il(lo(i)+1),so(i),zo(i),&
                          ev(i)/2.d0,ek(i)/2.d0,ep(i)/2.d0
   end do
   !c.hartwig give energies in hartree; no virial correction
   write(6,'(//," total energies",/,1x,14("-"),/,  &
   & /," sum of eigenvalues        =",f18.8,  &
   & /," kinetic energy from ek    =",f18.8,  &
   & /," el-ion interaction energy =",f18.8,  &
   & /," el-el  interaction energy =",f18.8,  &
   & /," vxc    correction         =",f18.8,  &
   & /," exchange + corr energy    =",f18.8,  &
   & /," kinetic energy from ev    =",f18.8,  &
   & /," potential energy          =",f18.8,/,1x,45("-"),  &
   & /," total energy              =",f18.8)') (etot(i)*.5d0,i=1,5), (etot(i)*.5d0,i=7,10)

end subroutine etotal


!> exit routine (i is a stop parameter)
!!    000-099 main     (0 is normal exit)
!!    100-199 input
!!    200-299 charge
!!    300-399 vionic
!!    400-499 velect
!!    500-599 dsolv1
!!    600-699 dsolv2   (including difnrl and difrel)
!!    700-799 etotal
!!    800-899 pseudo
subroutine ext(info)
   implicit none
   integer, intent(in) :: info  !< Stop parameter
   if (info /= 0) then
      write(6,'(17x,a,i3)') "stop parameter =", info
   end if
   stop
end subroutine ext


!> vionic sets up the ionic potential
!! note that viod,u is the ionic potential times r
!! Potentials: always multiplied by r.
!! viod,u ..... ionic potential (down,up)
subroutine vionic(ifcore,  &
      nr,r,rprb,lmax,  &
      znuc, &
      viod,viou)

   implicit none

   !Arguments
   integer, intent(out) :: ifcore                           !< if core
   integer, intent(in) :: nr                                !< \#radial mesh points
   real(kind=8), dimension(nr), intent(in) :: r             !< Radial mesh
   real(kind=8), intent(in) :: rprb                         !< Radius for the parabolic confinement potential
   real(kind=8), intent(in) :: znuc
   integer, intent(in) :: lmax                              !< l channel
   real(kind=8), dimension(lmax,nr), intent(out) :: viod    !< Ionic potential down
   real(kind=8), dimension(lmax,nr), intent(out) :: viou    !< ionic potential up
   !Local variables
   real(kind=8) :: vshift,rinvrprb2
   integer :: i,j

   !2*znuc part (Rydberg units)

   ifcore = 0
   rinvrprb2 = 1.d0/rprb**2
   do i=1,lmax
      do j=1,nr
         ! c.hartwig  add confining potential
         viod(i,j) = -2.0d0*(znuc -.5d0*(r(j)*rinvrprb2)**2*r(j))
         viou(i,j) = viod(i,j)
         !viou(i,j) = -2.0d0*(znuc -.5d0*(r(j)/rprb**2)**2*r(j))
         ! c.hartwig  shift potential to avoid positive eigenvalues
         ! and convergence problems
         vshift=-15.0d0*r(j)
         viod(i,j) = viod(i,j)+vshift
         viou(i,j) = viou(i,j)+vshift
      end do
   end do

end subroutine vionic


!> velect generates the electronic output potential from the electron charge density.
!! The ionic part is added in dsolve.
subroutine velect(iter,iconv,iXC,nspol,ifcore,  &
      nr,r,rab,rw,rd, &
      zel,cdd,cdu,cdc,  &
      vod,vou,  &
      etot)
   ! we need these modules to re-initialize libXC in case
   ! the first few iterations are done with LDA XC
   use libxcModule
   implicit none

   logical, parameter :: debug = .false.            !< Debug flag
   real(kind=8), parameter :: pi = 4.d0*atan(1.d0)
   !Arguments
   integer, intent(in) :: nr, ifcore, nspol, iconv, iter
   real(kind=8), dimension(10), intent(inout) :: etot
   real(kind=8), dimension(nr), intent(in) :: r,rab,rw,rd,cdc
   real(kind=8), dimension(nr), intent(out) :: cdd,cdu,vod,vou
   real(kind=8), intent(in) :: zel
   integer, intent(in) :: iXC
   !Local variables
   real(kind=8), dimension(nr) :: y,yp,ypp,s1,s2
   real(kind=8), dimension(3*nr) :: w
   !c.hartwig
   !convention for spol as in dsolv: spin down=1 and spin up=2
   real(kind=8), dimension(nr,nspol) :: rho,vxcgrd
   real(kind=8), dimension(nr) :: excgrd

   real(kind=8) :: a1,an,b1,bn,ehart,enexc,exc,exct,rhodw,rhoup
   real(kind=8) :: vxc,vxcd,vxcu,xlo,xnorm,rinvp4,rinv2
   integer :: i,ierr,ll,isx

   !fit cd/r by splines
   y(1) = 0.d0
   if (ifcore == 2) then
      do i=2,nr
         y(i) = (cdd(i)+cdu(i)+cdc(i))/r(i)
      end do
   else
      do i=2,nr
         y(i) = (cdd(i)+cdu(i))/r(i)
      end do
   end if

   isx = 0
   a1 = 0.d0
   an = 0.d0
   b1 = 0.d0
   bn = 0.d0
   call splift(r,y,yp,ypp,nr,w,ierr,isx,a1,b1,an,bn)

   !compute the integrals of cd/r and cd from r(1)=0 to r(i)
   xlo = 0.d0
   call spliq(r,y,yp,ypp,nr,xlo,r,nr,s2,ierr)
   !s2 ==    ans(i) = integral from xlo to xup(i)
   do i=1,nr
      ypp(i) = r(i)*ypp(i) + 2.d0*yp(i)
      yp(i)  = r(i)*yp(i)  + y(i)
      y(i)   = r(i)*y(i)
   end do
   call spliq(r,y,yp,ypp,nr,xlo,r,nr,s1,ierr)

   !check normalization

   xnorm = 0.d0
   if (zel /= 0.d0) xnorm = zel/s1(nr)


   !let us try this
   if (iter > 3 .and. abs(zel-s1(nr)) > 0.01d0) then
     if (zel < s1(nr)+1.d0 ) then
       write(6,'(/," warning *** charge density rescaled in",  &
         & " velect",/," iteration number",i4,3x,"scaling factor =",1pe10.3,/)') iter,xnorm
     else
       xnorm=.99d0*xnorm
       write(6,'(/," warning *** charge density partially rescaled in",  &
         & " velect",/," iteration number",i4,3x,"scaling factor =",1pe10.3,/)') iter,xnorm
     end if
   end if

   !rather than:
   !if (iter > 0 .and. abs(zel-s1(nr)) > 0.01d0) then
      !write(6,'(/,46h," warning *** charge density rescaled in velect", &
      !/,17h," iteration number",i4,3x,16h,"scaling factor =",g10.3,/)') iter,xnorm
   !end if


   !compute new hartree potential
   !renormalize the charge density

   if (debug) write(*,'(a,3(1pe22.15))') 'DEBUG: xnorm,s1(nr),s2(nr)',xnorm,s1(nr),s2(nr)
   do i=2,nr
      ! at this point, V is the same for spin up and down
      vod(i) = 2.d0 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
      !if (debug) write(*,'(a,i6,3(1pe22.15))') 'DEBUG: vod,s1,s2',i,vod(i),s1(i),s2(i)
      vou(i) = vod(i)
      cdd(i) = xnorm*cdd(i)
      !if (debug) write(*,'(a,i6,3(1pe22.15))') 'DEBUG: cdu cdd ',i,cdu(i),cdd(i)
      cdu(i) = xnorm*cdu(i)
   end do

   if (iconv == 1) then
      !compute hartree contribution to total energy
      !does not look spin polarized yet
      ehart = 0.d0
      ll = 4
      do i=2,nr
         ehart = ehart + ll * (cdd(i)+cdu(i)) * vod(i) * rab(i)
         !      ehart = ehart + ll * (cdd(i)*vod(i)+cdu(i)*vod(i))* rab(i)
         !      ^^ identical if vod=vou, which is the case for now
         ll = 6 - ll
      end do
      ehart = ehart / 6.d0
   end if

   !find derivatives of the charge density

   ! here was the functional specification section c
   rinvp4 = 1.d0/4.d0/pi
   if (nspol == 2) then
      do i=2,nr
         rinv2=rinvp4/r(i)**2
         rho(i,1)=rinv2*cdd(i)
         rho(i,2)=rinv2*cdu(i)
      end do
   else
      do i=2,nr
         rho(i,1)=rinvp4*(cdd(i)+cdu(i))/r(i)**2
      end do
   end if
    !do i=2,nr
       !so how do we define this line now:
       !rho(i)=(cdd(i)+cdu(i))/4.d0/pi/r(i)**2
       !CAREFUL: this looks clumsy: is this just some multiplying back and forth of rw(i)?!
       !if (nspol==2) then
          !rho(i,1)=(cdd(i))/4.d0/pi/r(i)**2
          !rho(i,2)=(cdu(i))/4.d0/pi/r(i)**2
       !else
          !rho(i,1)=(cdd(i)+cdu(i))/4.d0/pi/r(i)**2
       !end if
    !end do
    !some : added here Only for the first point!
    rho(1,:)=rho(2,:)-(rho(3,:)-rho(2,:))*r(2)/(r(3)-r(2))

    !CMK   this should avoid problems with XC-functionals
    !CMK   with kinks (BLYP,....)
    !if (iter<30) then
        !mfxcx=0
        !mfxcc=9
        !mgcc=0
        !mgcx=0
    !else if (iter==30) then

    if (iter==40.and.iXC/=-20) then
       write(6,*) 'Switching from LDA to the requested functional'
       write(6,*) ' iXC =',iXC,'nspol=',nspol
       call libxc_functionals_end()
       call libxc_functionals_init(iXC,nspol)
    end if

    !end if

    !hutter's routine
    !call evxc(nr,r,rho,vxcgrd,excgrd)
    !goedecker's routine
    !call ggaenergy_15(nspol,nr,r,rw,rd,rho,enexc,vxcgrd,excgrd)

    call driveXC(nspol,nr,r,rw,rd,rho,enexc,vxcgrd,excgrd)
    !rho and vxcgr are now of dimension  (ng,nspol)

    !c.hartwig modified integration
    exc=0.d0
    vxc=0.d0
    !need energy/potential in ryd


    ! This section was and is very inefficient.
    ! let us keep this style for now
    ! but not forget to clean it up later
    ! the factors of two that cancel each other
    ! are from the previous versions.
    if (nspol==1) then
       ! non-polarized case
       ! quite the same as in older versions
       do i=1,nr
          exct = 2.d0*excgrd(i)*rho(i,1)
          vxcd = 2.d0*vxcgrd(i,1)
          vxcu=vxcd
          rhodw=rho(i,1)/2.d0
          rhoup=rhodw
          vod(i) = vod(i) + vxcd
          vou(i) = vou(i) + vxcu
          vxc = vxc + (vxcd*rhodw + vxcu*rhoup) * rw(i)
          exc = exc + exct * rw(i)
       end do
    else
       ! spin polarized case
       ! same dirty style, but with two spin channels
       do i=1,nr
          exct =  2.d0*excgrd(i)*(rho(i,1)+rho(i,2))
          vxcd =  2.d0*vxcgrd(i,1)
          vxcu =  2.d0*vxcgrd(i,2)
          rhodw=rho(i,1)/2.d0
          rhoup=rho(i,2)/2.d0
          vod(i) = vod(i) + vxcd
          vou(i) = vou(i) + vxcu
          vxc = vxc + (vxcd*rhodw + vxcu*rhoup) * rw(i)
          exc = exc + exct * rw(i)
          ! write(18,*)vxc, vxcd
       end do
    end if

    !Finally update energy quantities
    etot(4) = ehart
    etot(5) = nspol*vxc
    etot(7) = exc

end subroutine velect


!> Subroutine to read input parameters  for one configuration and build the different integration grids
!!    ncore .. # core orbitals (closed shells)
!!    no ..... n quantum number
!!    lo ..... l do.
!!    so ..... spin (+/- 0.5, or 0 for unpolarized)
!!    zo ..... # electrons
!!    znuc ... atomic number
subroutine input(itype,iXC,ispp,  &
      nrmax,nr,a,b,r,rab,rw,rd,rprb,rcov, &
      nameat,norb,ncore, &
      maxorb,no,lo,so,zo,  &
      znuc,zel,zcore, &
      nconf, &
      nvalo,ncoreo)

   ! We need these modules to initialize libXC when reading iXC
   ! use defs_basis
   use libxcModule
   implicit none

   !Arguments
   integer, intent(in) :: nrmax                       !< Maximal number of radial mesh points
   integer, intent(out) :: nr                         !< # mesh points
   real(kind=8), dimension(nrmax), intent(out) :: r   !< Radial mesh r(i) = a*(exp(b*(i-1))-1)
   real(kind=8), dimension(nrmax), intent(out) :: rab !< rab(i) = (r(i)+a)*b (integration grid)
   real(kind=8), dimension(nrmax), intent(out) :: rw  !< rw(i) = rab(i)*12.56637061435917d0*r(i)**2
   real(kind=8), dimension(nrmax), intent(out) :: rd  !< rd(i) = 1/rab(i)
   real(kind=8), intent(out) :: a,b                   !< Parameters to build the logarithmic mesh
   real(kind=8), intent(out) :: rprb                  !< Radius for the parabolic confinement potential
   real(kind=8), intent(out) :: rcov                  !< Covalent radius
   integer, intent(in) :: maxorb                      !< Maximal number of orbitals
   integer, intent(out) :: norb                       !< Number of orbitals
   integer, dimension(maxorb), intent(out) :: no      !< First quantum number (n) for each orbital
   integer, dimension(maxorb), intent(out) :: lo      !< Second quantum number (l) for each orbital
   real(kind=8), dimension(maxorb), intent(out) :: so, zo
   character(len=2), intent(out) :: itype
   character(len=2), intent(out) :: nameat            !< Name of the atom
   character(len=1), intent(out) :: ispp              !< n for non relativistic calculations, r for relativistic, s for (relat) spin polarized calculations
   integer, intent(out) :: iXC                        !< integer code the eXchange-Correlation functional
   integer, intent(out) :: ncore
   real(kind=8), intent(out) :: znuc                  !< Nuclear charge (integer with the type real(kind=8))
   real(kind=8), intent(out) :: zel,zcore
   integer, intent(out) :: ncoreo,nvalo
   integer, intent(inout) :: nconf                    !< Electronic configuration id
                                                      !!    nconf==0 read all parameters for the mesh
                                                      !!    for any nconf read also the electronic configurations
   !Local variables
   integer, dimension(15), parameter :: nc = (/ 1,2,2,3,3,3,4,4,4,4,5,5,5,6,6 /)
   integer, dimension(15), parameter :: lc = (/ 0,0,1,0,1,2,0,1,2,3,0,1,2,0,1 /)
   character(len=1), parameter :: blank = ' '
   logical, parameter :: debug = .false.      !< Debug flag

   ! Spin polarization information
   integer, dimension(5) :: nomin = (/ 10, 10 ,10, 10, 10 /)
   ! For use in routine atomwr
   character(len=80) :: instrg
   character(len=3) :: irel
   character(len=3) :: name

   real(kind=8) :: charge
   real(kind=8) :: aa,bb,rmax,sc,si,xji,zd,zion,zu,zval
   integer :: i,ierr,j,j1,j2,jmax,li,ni,nspol,nval

   !Read all electron configuration
   itype='ae'

   !Loop for reading all electron configurations
   !i.e. the first line of the file atom.dat
   loop_next: do
      !Read instruction
      read(35,'(a)',iostat=ierr) instrg
      call error_iostat()
      if (ierr < 0) return
      !If blank lines, next line
      if (instrg == ' ') cycle
      !next configuration ?
      if (index(instrg,'NEXT CONFIGURATION') /= 0 ) exit
      !Ignore the lines if nconf > 1 otherwise read infomration about the mesh
      if (nconf == 0) exit
   end do loop_next

   !For the first electronic configuration only !!
   if (nconf == 0) then
      j1=1
      j2=2
      do i=len(instrg),1,-1
         if (instrg(i:i)/=' ') j1=i
      end do
      do i=len(instrg),j1,-1
         if (instrg(i:i)==' ') j2=i
      end do
      j2=j2-1
      nameat=instrg(j1:j2)
      if (j2==1) nameat(2:2)=' '
      read(35,'(a)',iostat=ierr) instrg
      call error_iostat()
      if (ierr < 0) return
      j1=1
      j2=2
      do i=len(instrg),1,-1
         if (instrg(i:i)/=' ') j1=i
      end do
      do i=len(instrg),j1,-1
         if (instrg(i:i)==' ') j2=i
      end do
      j2=j2-1

      ! READ iXC (Information about exchange-correlation functional)
      ! for now, only keep the two most commonly used functionals
      ! backwards compatible. Otherwise, require ABINITs iXC < 0
      if    (instrg(j1:j2)=='PADE') then
         iXC=-20
      elseif (instrg(j1:j2)=='PBE') then
         iXC=-101130
      else
         read(instrg(j1:j2),*,iostat=ierr) iXC
         if (ierr/=0) then
          write(6,'(1x,a)' ) 'Could not read the XC input in atom.dat'
          stop
         end if
      end if

      ! Read information about spin polarization
      read(35,'(a)',iostat=ierr) instrg
      call error_iostat()
      if (ierr < 0) return
      do i=len(instrg),1,-1
         if (instrg(i:i)/=' ') j1=i
      end do
      ispp=instrg(j1:j1)
      if (ispp=='R') ispp='r'
      if (ispp/='r'.and.ispp/='n'.and.ispp/='s') then
         write(6,'(1x,a)')       'The first non-blank character on line 3'
         write(6,'(1x,a)')       'of atom.dat must be one of'
         write(6,'(1x,a)')       'n: for non relativistic calculations'
         write(6,'(1x,a)')       'r: for relativistic calculations'
         write(6,'(1x,a)')       's: for (relat) spin polarized calculations'
         write(6,'(/,1x,a,1x,a)') 'Character found:',ispp
         write(6,'(1x,a)')        'Exiting.'
         stop
      end if

      ! if (ispp /= 's' .and. ispp  /= 'r')  ispp=blank
      ! spin-polarization needs relativistic calculation
      znuc=0.d0
      read(35,*,iostat=ierr) rmax,aa,bb
      call error_iostat()
      if (ierr < 0) return

      read(35,*,iostat=ierr) rcov,rprb
      call error_iostat()
      if (ierr < 0) return

      znuc=charge(nameat)

      ! Set up grid
      if (abs(rmax) < 0.00001) rmax=100.0d0 
      if (abs(aa)   < 0.00001) aa = 3.0d0
      if (abs(bb)   < 0.00001) bb = 40.0d0

      if (znuc == 0.0d0) then
         a = 10**(-aa)
      else
         a=exp(-aa)/znuc
         b = 1.d0/bb

         !modify grid-parameter, so that one grid-point matches
         !rcov exactly (r(i) is used in this part only for that)
         do i=1,nrmax
            ! For i = 0 r(i) = 0
            ! For i = 2 r(i) = a*(exp(b)-1)
            ! For i>2 r(i+1) + a =r(i)*exp(b) + a
            r(i) = a*(exp(b*(i-1))-1.d0)
            if (r(i) >= rcov) then
               a= rcov/(exp(b*(i-1))-1.d0)
               aa=-log(a*znuc)
               !We have finished!
               !write(*,'(1x,a,f21.15)') 'Adjusted value for aa',aa
               exit
            end if
         end do
         !Fail to match rocv
         if (i > nrmax) then
            write(6,'(/," Error in input - arraylimits", " for radial array exceeded",/)')
            stop 'input two'
         end if
      end if

      !Build all integration grids (r, rab, rw and rd)
      do i=1,nrmax
         if (i == nrmax) then
            write(6,'(/," error in input - arraylimits", " for radial array exceeded",/)')
            call ext(100)
         end if
        r(i) = a*(exp(b*(i-1))-1)
        rab(i) = (r(i)+a)*b
        ! c.hartwig: set up grids for modified integration
        rw(i) = b*(r(i)+a)
        rd(i) = 1.d0/rw(i)
        rw(i)=rw(i)*12.56637061435917d0*r(i)**2
        if (r(i) > rmax) exit
      end do

      !Set the number of grid points (< rmax)
      nr = i-1

      !Modify weights at end point for improved accuracy
      if (debug) then
         write(*,*) 'DEBUG OPTION: No modified weights at origin!'
      end if
      rw(1)=rw(1)*17.d0/48.d0
      rw(2)=rw(2)*59.d0/48.d0
      rw(3)=rw(3)*43.d0/48.d0
      rw(4)=rw(4)*49.d0/48.d0


      !read the number of core and valence orbitals
      read(35,*,iostat=ierr) ncore, nval
      call error_iostat()
      if (ierr < 0) return

      nvalo=nval
      ncoreo=ncore
      if (ncore > 15) then
         write(6,*) 'more than 15 core orbitals'
         call ext(101)
      end if

   end if !End of the reading for the first configuration


   ! For all configurations
   ncore = ncoreo
   nval = nvalo
   if (ispp == 'R') ispp='r'
   !if (ispp/='r') ispp=' '
   if (ispp /= 's' .and. ispp  /= 'r')  ispp=blank
   nspol=1
   if (ispp=='s') nspol=2

   !compute occupation numbers and orbital energies for the core

   !the following section is not quite clear.
   zcore = 0.d0
   if (ispp == blank) then
      jmax = 1
      sc = 0.0d0
   else
      jmax = 2
      sc = - 0.5d0
   end if

   norb = 0
   if (ncore > 0) then
      do i=1,ncore
         do j=1,jmax
            norb = norb + 1
            no(norb) = nc(i)
            lo(norb) = lc(i)
            so(norb) = sc
            zo(norb) = 2.d0*lo(norb)+1
            ! why not do the same in the 's' case?
            if (ispp == blank) zo(norb) = 2.d0*zo(norb)
            if (ispp == 'r')   zo(norb) = 2.d0*(lo(norb)+sc) + 1.d0

            ! there must be a reason that
            ! the convention for zo is
            ! 4l+2           'n' or ''
            ! 2l+(1 or 2)    'r'
            ! 2l+1           's'
            zcore = zcore + zo(norb)
            if (abs(zo(norb)) < 0.1d0) norb=norb-1
            if (ispp /= blank) sc=-sc
         end do
      end do
      ncore = norb
   end if

   norb = ncore
   !end of core orbital energies and occupations

   zval = 0.d0

   if (nval /= 0) then

      do i=1,nval
         read(35,*,iostat=ierr) ni,li,zd,zu
         call error_iostat()
         if (ierr < 0) return
         si = 0.d0
         if (ispp /= blank) si=0.5d0

         do j=1,jmax
            norb = norb + 1
            if (ispp /= blank) si=-si
            no(norb) = ni
            lo(norb) = li
            so(norb) = si
            zo(norb) = zd+zu

             ! c.hartwig
             if (zo(norb) == 0.0) zo(norb)=1.0d-20

             ! this is an experimental option:
             ! zd > 0 = zu  ---> use Hund s rule
             ! for auto assignment in polarized case:
             if (ispp == 's') then
                if (zu==0d0 .and. zd>0d0 .and. j==1 ) then
                    zd = min( dble(2*li+1), zo(norb) )
                    zu = zo(norb)-zd
                    ! write(*,*)"(Hunds rule)",ni,li,si,zd,zu
                end if
                if ( si < 0.1d0) then
                    zo(norb) = zd
                else
                    zo(norb) = zu
                end if
             end if

             ! assign occupation numbers for spin orbit coupling in the relativistic case 'r'
             if (ispp == 'r') zo(norb)=zo(norb)*(2*(li+si)+1)/(4*li+2)
             zval = zval + zo(norb)
             ! no s down orbital in the 'r' case
             if (ispp == 'r' .and. li+si < 0.d0) norb=norb-1
             if (norb == 0) cycle
             if (nomin(lo(norb)+1) > no(norb)) nomin(lo(norb)+1)=no(norb)
          end do
       end do

       nval = norb - ncore

       !abort if two orbitals are equal

       if (norb <= 0) call ext(110)
       do i = 1, (norb - 1)
          do j = (i + 1),norb
             if (no(i) == no(j) .and. lo(i) == lo(j)) then
                if (abs(so(i)-so(j)) < 1.0d-3) then
                   write(*,'(1x,a,3(i0,1x),f6.1)') 'i,no(i),lo(i),so(i):',i,no(i),lo(i),so(i)
                   write(*,'(1x,a,3(i0,1x),f6.1)') 'j,no(j),lo(j),so(j):',j,no(j),lo(j),so(j)
                   call ext(110+i)
                end if
              end if
          end do
          ! print*,'i,no(i),lo(i),so(i):',i,no(i),lo(i),so(i)
       end do

    end if !End of the initialization for the valence orbitals

    zion = znuc - zcore - zval
     !write(6,*)' zion = ',zion
     !write(6,*)' znuc = ',znuc
     !write(6,*)' zcore = ',zcore
     !write(6,*)' zval = ',zval
    zel = zval
    zel=zel+zcore

    write(6,'(1x,a2," all electron calculation  ",/,1x,27("-"),/)') nameat
    if (ispp == 'r') write(6,'(" r e l a t i v i s t i c ! !",/)')
    if (ispp /= 's') then
       name = 'non'
    else
       name = '   '
    end if
    write(6,'(" iXC = ",i7,3x,a3," spin-polarized",/)') iXC,name
    write(6, &
         '(" nuclear charge             =",f10.6,/,  &
         & " number of core orbitals    =",i3,/,  &
         & " number of valence orbitals =",i3,/,  &
         & " electronic charge          =",f10.6,/,  &
         & " ionic charge               =",f10.6,//)') znuc,ncore,nval,zel,zion

    write(6,'(" input data for orbitals",//,"  i    n    l    s     j     occ",/)')
    
    write(6,'(1x,a)') ' Using LDA for generating the input guess wfn'
    call libxc_functionals_init(-20,nspol)
    
    !write(6,*)' initializing libXC with iXC=',iXC,'spin',nspol
    !call libxc_functionals_init(iXC,nspol)
    
    xji = 0.d0
    do i=1,norb
       if (ispp == 'r') xji = lo(i) + so(i)
       write(6,'(1x,i2,2i5,2f6.1,f10.4)') i,no(i),lo(i),so(i),xji,zo(i)
    end do
    
    write(6,'(//," radial grid parameters",//," r(1) = .0 , r(2) =",1pe13.6," , ... , r(",&
       &      i4,") =",0pf12.8,/," a =",f12.8,"  b =",f12.8,/)') r(2),nr,r(nr),aa,bb
    
    irel   = 'nrl'
    if (ispp == 'r') irel = 'rel'
    if (ispp == 's') irel = 'spp'
     !do i = 1, norb
     !   write(text(i),'(1x,i1,a," s=",f4.1," (occ=",f6.3,") ",a)') no(i),spdf(lo(i)+1),so(i),zo(i),irel
     !end do

contains

   subroutine error_iostat()
      implicit none
      if (ierr > 0) then
         !Error of reading
          write(6,'(1x,a)') 'Error while reading atom.dat'
          stop
       else if (ierr < 0) then
          !Error: end of file
          write(6,'(1x,a)') 'Reached end of file atom.dat'
          !itype gives and error code
          itype='st'
          return
       end if
    end subroutine error_iostat

end subroutine input


 !> Function determines the nuclear charge of an element
double precision function charge(name)
   implicit none
   !Arguments
   character(len=2), intent(in) :: name
   !Local variables
   integer, parameter :: nelem = 103
   !> The periodic table
   character(len=2), dimension(nelem), parameter :: pertab = (/ &
    'H ','HE',  &
    'LI','BE','B ','C ','N ','O ','F ','NE',  &
    'NA','MG','AL','SI','P ','S ','CL','AR',  &
    'K ','CA',  &
         'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',  &
              'GA','GE','AS','SE','BR','KR',  &
    'RB','SR',  &
         'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',  &
              'IN','SN','SB','TE','I ','XE',  &
    'CS','BA',  &
         'LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',  &
                                  'HO','ER','TM','YB','LU',  &
              'HF','TA','W ','RE','OS','IR','PT','AU','HG',  &
              'TL','PB','BI','PO','AT','RN',  &
    'FR','RA',  &
         'AC','TH','PA','U ','NP','PU','AM','CM','BK','CF',  &
                                  'ES','FM','MD','NO','LR' /)
   character(len=2) :: elemnt
   integer, dimension(2) :: ic
   integer :: i

   ! convert the name to upper-case, and possibly left-justify
   ! code 97-122: lower case
   ! code 65-90:  upper case
   ! code 32:     blank

   !Read the two characters
   do i = 1,2
      ! get the ascii value
      ic(i) = ichar( name(i:i) )
      if (ic(i) >= 97 .and. ic(i) <= 122) then
         ! convert to upper case
         ic(i) = ic(i) - 32
      else if (ic(i) >= 65 .and. ic(i) <= 90) then
         ! upper-case - do nothing
      else if (ic(i) == 32) then
         ! 'space' - do nothing
      else if (ic(i) == 0) then
         ! 'nul' - replace by space
         ic(i) = 32
      else
         write (6,*) 'unrecognized element name:',name
         call ext(200)
      end if
   end do

   ! left justify
   if (ic(1) == 32) then
     ic(1) = ic(2)
     ic(2) = 32
   end if
   ! the standard name of the element:
   elemnt = char(ic(1))//char(ic(2))

   ! find the element in the periodic table
   do i = 1, nelem
      if (elemnt == pertab(i)) then
         charge = real(i,kind=8)
         return
      end if
   end do
   write(6,'(" could not locate name in list of elements:",/," name=",a," converted to=",a, &
      &       " ascii codes=",2i3)') name,elemnt,ic
   call ext (200)
end function charge


!> dsolv1 finds the non relativistic wave function
!! using finite differences and matrix diagonalization
!! initial guess for the eigenvalues need not be supplied
!! @todo
!!  Use Lapack instead of the routines from atom.splines.f90
subroutine dsolv1( &
      nr,a,b,r,rab,lmax, &
      norb,no,lo,so,zo, &
      cdd,cdu, &
      viod,viou,vid,viu, &
      ev)

   implicit none

   !Arguments
   integer, intent(in) :: nr,norb,lmax
   real(kind=8), dimension(nr), intent(in) :: r,rab,vid,viu
   real(kind=8), dimension(nr), intent(out) :: cdd,cdu
   real(kind=8), dimension(lmax,nr), intent(in) :: viod,viou
   integer, dimension(norb), intent(in) :: no,lo
   real(kind=8), dimension(norb) :: so,zo,ev
   real(kind=8), intent(in) :: a,b
   !Local variables!
   logical, parameter :: debug = .false.  !< Debug flag
   integer, parameter :: nomax = 10       !< Maxval(no)
   integer, dimension(2,5) :: nmax
   integer, dimension(10) :: ind
   real(kind=8), dimension(nr) :: dk,d,sd,sd2,rv1,rv2,rv3,rv4,rv5
   real(kind=8), dimension(10) :: e
   real(kind=8), dimension((nr-1)*nomax) :: z
   real(kind=8) :: bl,bu,c1,c2,denr,eps
   integer :: i,j,k,ki,l,llp,nrm,nvmax,ierr,kn

   !initialize charge density arrays
   d=0.0d0
   do i=1,nr
      cdd(i) = 0.d0
      cdu(i) = 0.d0
   end do
   nvmax = (nr-1)*nomax

   !find max n given l and s
   !zero spin is treated as down
   do i=1,2
      do j=1,lmax
         nmax(i,j) = 0
         do k=1,norb
            if (no(k) <= 0) cycle
            if (lo(k) /= j-1) cycle
            if ((so(k)-0.1d0)*(real(i,kind=8)-1.5d0) < 0.d0) cycle
            nmax(i,j)=no(k)
            if (no(k)*(nr-1) > nvmax) then
              write(*,'(2(1x,i0))') no(k),nr-1
              !Max dimension of z used by tinvit
              write(*,'(1x,i0,a,i0)') no(k)*(nr-1)," > ",nvmax
              call ext(500)
            end if
         end do
      end do
   end do

   !set up hamiltonian matrix for kinetic energy
   !only the diagonal depends on the potential
   c2 = -1.d0/b**2
   c1 = -2.d0*c2 + 0.25d0
   dk(1)  = c1 / (r(2)+a)**2
   sd(1)  = 0.d0
   sd2(1) = 0.d0
   do i=3,nr
      dk(i-1)  = c1 / (r(i)+a)**2
      sd(i-1)  = c2 / ((r(i)+a)*(r(i-1)+a))
      sd2(i-1) = sd(i-1)**2
   end do

   !start loop over spin down=1 and spin up=2
   nrm = nr - 1
   do i=1,2

      !start loop over s p d... states
      do j=1,lmax
         if (nmax(i,j) == 0) cycle
         llp = j*(j-1)
         do k=2,nr
            if (i == 1) &
               d(k-1) = dk(k-1) + (viod(j,k) + llp/r(k))/r(k) + vid(k)
            if (i == 2) &
               d(k-1) = dk(k-1) + (viou(j,k) + llp/r(k))/r(k) + viu(k)
            if (debug) then
               write(*,*) 'debug: vio u d (k)',k,viou(j,k),viod(j,k)
               write(*,*) 'debug: vi u d (k)',k,viu(k),vid(k)           !!! NaN
               write(*,*) 'debug: r (k)',k,r(k)
               write(*,*) 'debug: dk (k)',k-1,dk(k-1)
            end if
         end do

         ! diagonalize

         eps = -1.d0
         ! Find eigenvalues of a tridiagonal symmetric matrix between specific boundary indices using bisection
         call tridib(nrm,eps,d,sd,sd2,bl,bu,1,nmax(i,j),e,ind,ierr, rv4,rv5)
         if (ierr /= 0) write(6,'(/," ****** error  ierr =",i3,/)') ierr
         !Warning: z(nrm,nmax(i,j))
         ! Find eigenvectors of a tridiagonal symmetric matrix corresponding 
         ! to specified eigenvalues, using inverse iteration.
         call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr, rv1,rv2,rv3,rv4,rv5)
         if (ierr /= 0) write(6,'(/," ****** error  ierr =",i3,/)') ierr

         ! save energy levels and add to charge density
         ki = 1
         kn = 0
         do k=1,norb
            if (no(k) <= 0) cycle
            if (lo(k) /= j-1) cycle
            ! if spin(k) /= spin(i) cycle
            if ((so(k)-0.1d0)*(i-1.5d0) < 0.d0) cycle
            ev(k) = e(ki)
            ! write(6,*)'DSOLV1:',k,no(k),lo(k),so(k),ev(k)
            do l=2,nr
               denr = zo(k) * z(kn+l-1)**2 / rab(l)
               if (i == 1) cdd(l) = cdd(l) + denr
               if (i == 2) cdu(l) = cdu(l) + denr
            end do
            ki = ki + 1
            kn = kn + nrm
         end do

      end do ! end loop over s p and d states
   end do    ! end loop over spin down=1 and spin up=2

end subroutine dsolv1


!> dsolv2 finds the (non) relativistic wave function using
!! difnrl to integrate the Schroedinger equation or
!! difrel to integrate the Dirac equation
!! the energy level from the previous iteration is used
!! as initial guess, and it must therefore be reasonable
!! accurate.
subroutine dsolv2(iter,iconv,iXC,ispp,ifcore, &
      nr,a,b,r,rab,lmax,  &
      norb,ncore,no,lo,so,zo,  &
      znuc,zcore,cdd,cdu,cdc,dcrc,ddcrc,  &
      viod,viou,vid,viu, &
      ev,ek,ep,rcov,rprb,nconf)

   implicit none

   !Arguments
   integer, intent(in) :: iter,iconv,ifcore,ncore,lmax,nconf
   integer, intent(in) :: nr,norb
   real(kind=8), dimension(nr), intent(in) :: r,rab,vid,viu
   real(kind=8), dimension(nr), intent(out) :: cdd,cdu,cdc
   integer, dimension(norb), intent(in) :: no,lo
   real(kind=8), dimension(norb), intent(in) :: so, zo
   real(kind=8), dimension(lmax,nr), intent(in) :: viod,viou
   real(kind=8), dimension(norb), intent(out) :: ev,ek,ep
   real(kind=8), intent(out) :: dcrc,ddcrc
   real(kind=8), intent(in) :: a,b,znuc,zcore,rcov
   real(kind=8), intent(in) :: rprb                         !< Radius for the parabolic confinement potential
   character(len=1) :: ispp
   integer, intent(in) :: iXC

   !Local variables
   real(kind=8), dimension(nr) :: v,ar,br
   real(kind=8) :: denr
   integer :: i,j,lp,llp

   !initialize arrays for charge density
   do i=1,nr
      cdd(i) = 0.d0
      cdu(i) = 0.d0
      if (ifcore /= 1) cdc(i)=0.d0
   end do
   !and the moments of the core charge density
   dcrc =0d0
   ddcrc=0d0

   !start loop over orbitals
   !note that spin zero is treated as down

   do i=1,norb
      if (no(i) <= 0) cycle
      if (zo(i) == 0.d0 .and. iconv == 0) cycle
      if (ev(i) >= 0.d0) ev(i)=-1.d0

      !set up potential

      lp  = lo(i)+1
      llp = lo(i)*lp
      do j=2,nr
         if (so(i) < 0.1d0) v(j) = viod(lp,j)/r(j) + vid(j)
         if (so(i) > 0.1d0) v(j) = viou(lp,j)/r(j) + viu(j)
         if (ispp /= 'r') v(j) = v(j) + llp/r(j)**2
         !if (ispp == 'n') v(j) = v(j) + llp/r(j)**2
      end do

      !call integration routine
      if (ispp /= 'r' ) then
          call difnrl(iter,i,v,ar,br,  &
              lmax,nr,a,b,r,rab,  &
              norb,no,lo,so,  &
              znuc,  &
              viod,viou,vid,viu,ev)
      end if
      if (ispp == 'r' ) then
          call difrel(iter,i,v,ar,br,  &
              nr,r,rab,norb,  &
              no,lo,so,znuc,vid,viu,ev)
      end if

      !add to the charge density

      do j=1,nr
         denr = zo(i) * ar(j) * ar(j)
         !the relativistic case requires the minor component of the spinor to be added
         if (ispp == 'r') denr = denr + zo(i) * br(j) * br(j)
         if (so(i) < 0.1d0) cdd(j) = cdd(j) + denr
         if (so(i) > 0.1d0) cdu(j) = cdu(j) + denr
         if (ifcore /= 1 .and. i <= ncore) cdc(j)=cdc(j)+denr
      end do

      !compute various quantities if last iteration

      if (iconv == 1) then
          !orban is used to analyze and printout data about the orbital
          call orban(iXC,ispp,i,ar,br,  &
             nr,r,rab, &
             lmax,norb,ncore,no,lo,so,zo,  &
             znuc,zcore,cdd,cdu,cdc,dcrc,ddcrc,  &
             viod,viou,vid,viu, &
             v,ev,ek,ep,rcov,rprb,nconf)
      end if
   end do ! end loop over orbitals

end subroutine dsolv2


!> difnrl integrates the Schroedinger equation
!! it finds the eigenvalue ev, the wavefunction ar and the derivative br = d(ar)/dr
subroutine difnrl(iter,iorb,v,ar,br,lmax,  &
      nr,a,b,r,rab,norb,no,lo,so,znuc,viod,viou,  &
      vid,viu,ev)

   implicit none

   !Arguments
   integer, intent(in) :: iter
   integer, intent(in) :: iorb  !< Number of the considered orbital
   integer, intent(in) :: lmax
   integer, intent(in) :: norb,nr
   integer, dimension(norb), intent(in) :: no,lo
   real(kind=8), dimension(norb), intent(in) :: so
   real(kind=8), dimension(norb), intent(inout) :: ev     !< Eigenvalues
   real(kind=8), dimension(nr), intent(in) :: v,vid,viu   !< Potential quantities
   real(kind=8), dimension(nr), intent(in) :: r,rab       !< Integration grid
   real(kind=8), dimension(nr), intent(out) :: ar         !< Wavefunctions
   real(kind=8), dimension(nr), intent(out) :: br         !< derivatives br = d(ar)/dr
   real(kind=8), dimension(lmax,nr) :: viod,viou
   real(kind=8), intent(in) :: a,b,znuc
   !Local variables
   !> Maximum number of integration step
   integer, parameter :: max_itmax = 100
   !> Tolerance
   real(kind=8), parameter :: etol=-1.d-7
   real(kind=8), parameter :: tol=1.0d-14
   !> Arrays added to gain speed.
   real(kind=8), dimension(5) :: rabrlo,rlp
   real(kind=8), dimension(nr) :: rab2,fa,fb
   real(kind=8) :: brp,brctp,brc,bb,arp,arctp,arc,amc4,amc3,amc2,amc1,amc0,alf
   real(kind=8) :: aa,abc1,abc2,abc3,abc4,abc5,dev,emax,emin,evold,expzer,factor
   real(kind=8) :: fb0,fb1,temp,var0,vev,vzero,zeff
   integer :: icount,istop,itmax,j,j1,j2,j3,j4,J5,juflow,ll,nctp,ninf,ninf1,ninf2,ninf3,ninf4,lp
   integer :: nodes

   ! Machine dependent parameter-
   ! Require exp(-2*expzer) to be within the range of the machine
   ! IBM
   expzer = 3.7d2
   !Iris     expzer = 3.7E2
   !Apollo   expzer = 3.7D2
   !Sun      expzer = 3.7D2
   !Vax      expzer = 44.d0
   !ray      expzer =  2.8E3

   !for numerical stability:
   expzer = expzer/2.d0

   !integration coefficients
   abc1 = 1901.d0/720.d0
   abc2 = -1387.d0/360.d0
   abc3 = 109.d0/30.d0
   abc4 = -637.d0/360.d0
   abc5 = 251.d0/720.d0
   amc0 = 251.d0/720.d0
   amc1 = 323.d0/360.d0
   amc2 = -11.d0/30.d0
   amc3 = 53.d0/360.d0
   amc4 = -19.d0/720.d0

   lp = lo(iorb)+1
   ar(1) = 0.0d0
   if (lo(iorb) == 0) then
      br(1) = b*a
   else
      br(1) = 0.0d0
   end if
   do j=2,nr
      ar(j) = 0.0d0
   end do
   do j=2,nr
      br(j) =0.0d0
   end do
   do j=2,5
      rlp(j)=r(j)**lp
   end do
   do j=2,5
      rabrlo(j)=rab(j)*r(j)**lo(iorb)
   end do
   do j=1,nr
      rab2(j)=rab(j)*rab(j)
   end do

   !set underflow trap

   juflow=1
   do j=2,nr
      if (lp*abs(log(r(j))) >= expzer/2.d0) juflow = j
   end do

   !determine effective charge and vzero for startup of
   !outward integration
   !ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
   !aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)

   zeff = 0.0d0
   if (so(iorb) < 0.1d0 .and. viod(lp,2) < -0.1d0) zeff=znuc
   if (so(iorb) > 0.1d0 .and. viou(lp,2) < -0.1d0) zeff=znuc
   aa = -zeff/lp
   vzero = -2.d0*zeff*aa
   if (zeff == 0.0d0) then
      if (so(iorb) < 0.1d0 ) then
         vzero=vzero+viod(lp,2)/r(2)
      else
         vzero=vzero+viou(lp,2)/r(2)
      end if
   end if
   if (so(iorb) < 0.1d0) then
      vzero=vzero+vid(2)
   else
      vzero=vzero+viu(2)
   end if
   var0 = 0.0d0
   if (lo(iorb) == 0) var0=-2.d0*zeff
   if (lo(iorb) == 1) var0=2.0d0
   emax = 0.0d0
   emin = -200000.0d0
   if (ev(iorb) > emax) ev(iorb) = emax

   !Main loop for the integration
   main_loop: do itmax=1,max_itmax+1
      if (itmax >= max_itmax) write(6,'(" iorb =",i3," iter =",i3," ev =",1pe18.10," nodes =",i4)') iorb,iter,ev(iorb),nodes
      !Too many iterations !!
      if (itmax > max_itmax) return

      if (ev(iorb) > 0.d0) then
         write(6,'(//," error in difnrl - ev(",i2,") greater then v(infinty)")') iorb
         stop 'difnrl one'
      end if

      !Find practical infinity ninf and classical turning
      !point nctp for orbital
      do icount=0,100
         do j=nr,2,-1
            temp = v(j) - ev(iorb)
            if (temp < 0.0d0) temp = 0.0d0
            if (r(j)*sqrt(temp) < expzer) exit
         end do

         ninf=j
         nctp = ninf - 5
         do j=2,ninf-5
           if (v(j) < ev(iorb)) nctp = j
         end do
         if (ev(iorb) >= etol*10.d0) nctp=ninf-5
         if (ev(iorb) >= etol) ev(iorb)=0.0d0
         !Find it!
         if (nctp > 6) exit
         !Decrease ev(iorb) and start again
         ev(iorb) = 0.9d0*ev(iorb)
      end do
      if (icount > 100) then
         write(6,'(///,"error in difnrl - cannot find the classical ",/" turning point for orbital ",i2)') iorb
         stop 'difnrl two'
      end if

      !outward integration from 1 to nctp
      !startup

      bb = (vzero-ev(iorb))/(4*lp+2)
      do j=2,5
         ar(j) = rlp(j) * (1+(aa+bb*r(j))*r(j))
         br(j) = rabrlo(j) * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
      end do

      !Predictor-corrector array added.
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-ev(iorb))*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-ev(iorb))*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-ev(iorb))*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-ev(iorb))*ar(5)


      ! integration loop
      nodes = 0
      do j=6,nctp

         ! predictor (Adams-Bashforth)
         j1=j-1
         j2=j-2
         j3=j-3
         j4=j-4
         j5=j-5
         vev=v(j)-ev(iorb)
         arp = ar(j1) + abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+  &
               abc4*fa(j4)+abc5*fa(j5)
         brp = br(j1) + abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+  &
               abc4*fb(j4)+abc5*fb(j5)
         fb1 = b*brp + rab2(j)*vev*arp

         ! corrector (Adams-Moulton)
         arc = ar(j1) + amc0*brp+amc1*fa(j1)+amc2*fa(j2)+  &
               amc3*fa(j3)+amc4*fa(j4)
         brc = br(j1) + amc0*fb1+amc1*fb(j1)+amc2*fb(j2)+  &
               amc3*fb(j3)+amc4*fb(j4)
         fb0 = b*brc + rab2(j)*vev*arc

         ! error reduction step
         ar(j) = arc + amc0*(brc-brp)
         br(j) = brc + amc0*(fb0-fb1)
         fa(j) = br(j)
         fb(j) = b*br(j) + rab2(j)*vev*ar(j)
         
         ! count nodes - if no underflow
         if (j > juflow .and. ar(j)*ar(j-1) < 0.d0) nodes=nodes+1
      end do

      arctp = ar(nctp)
      brctp = br(nctp)

      !end outward integration

      ! if number of nodes correct, start inward integration
      ! else modify energy stepwise and try again
      if (nodes /= no(iorb)-lo(iorb)-1) then
         !     c.hartwig
         !         write(6,*) 'nodes,ev(iorb)',nodes,ev(iorb)
         if (nodes < no(iorb)-lo(iorb)-1) then

            !  too few nodes; increase ev
            if (ev(iorb) > emin) emin = ev(iorb)
            ev(iorb) = ev(iorb) - ev(iorb)/10.d0
         else

            !  too many nodes; decrease ev
            if (ev(iorb) < emax) emax = ev(iorb)
            ev(iorb) = ev(iorb) + ev(iorb)/10.d0
         end if
         !New integration step
         cycle main_loop
      end if

      !inward integration from ninf to nctp
      !startup

      do j=ninf,ninf-4,-1
         alf = v(j) - ev(iorb)
         if (alf < 0.d0) alf = 0.0d0
         alf = sqrt(alf)
         ar(j) = exp(-alf*r(j))
         br(j) = -rab(j)*alf*ar(j)
      end do

      !Array for predictor-corrector added.

      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*  &
       (v(ninf)-ev(iorb))*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*  &
             (v(ninf1)-ev(iorb))*ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*  &
             (v(ninf2)-ev(iorb))*ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*  &
             (v(ninf3)-ev(iorb))*ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*  &
             (v(ninf4)-ev(iorb))*ar(ninf4)

      !integration loop
      istop = ninf - nctp
      if (istop >= 5) then

         do j=ninf-5,nctp,-1

            !predictor (Adams-Bashforth)
            j1 = j + 1
            j2 = j + 2
            j3 = j + 3
            j4 = j + 4
            j5 = j + 5
            vev = v(j)-ev(iorb)
            arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+  &
                  abc4*fa(j4)+abc5*fa(j5))
            brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+  &
                  abc4*fb(j4)+abc5*fb(j5))
            fb0 = b*brp + rab2(j)*vev*arp

            !corrector (Adams-Moulton)
            arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+  &
                  amc3*fa(j3)+amc4*fa(j4))
            brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+  &
                  amc3*fb(j3)+amc4*fb(j4))

            fb1 = b*brc + rab2(j)*vev*arc

            !error reduction step
            ar(j) = arc - amc0*(brc-brp)
            br(j) = brc - amc0*(fb1-fb0)
            fa(j) = br(j)
            fb(j) = b*br(j) + rab2(j)*vev*ar(j)
         end do
         !end inward integration
      end if

      !rescale ar and br outside nctp to match ar(nctp) from
      !outward integration
      factor = arctp/ar(nctp)
      do j=nctp,ninf
         ar(j) = factor * ar(j)
         br(j) = factor * br(j)
      end do

      !find normalizing factor

      factor = 0.0d0
      ll = 4
      do j=2,ninf
         factor = factor + ll*ar(j)*ar(j)*rab(j)
         ll = 6 - ll
      end do
      factor = factor / 3.d0

      !modify eigenvalue ev
      dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
      if (5.d0*abs(dev) > -ev(iorb)) dev=sign(ev(iorb),dev)/5.d0
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) > emax) ev(iorb) = (evold + emax) / 2.d0
      if (ev(iorb) < emin) ev(iorb) = (evold + emin) / 2.d0
      !Check the convergence
      if (abs(dev) <= tol*(1.d0-ev(iorb))) exit main_loop

   end do main_loop ! end of the integration loop


   !normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
   factor = 1.d0 / sqrt(factor)
   do j=1,ninf
      ar(j) = factor*ar(j)
      br(j) = factor*br(j) / rab(j)
   end do

end subroutine difnrl


!> difrel integrates the relativistic Dirac equation
!! it finds the eigenvalue ev, the major and minor component
!! of the wavefunction, ar and br.  It uses an initial guess
!! for the eigenvalues from dsolv1
subroutine difrel(iter,iorb,v,ar,br,nr,r,rab,  &
      norb,no,lo,so,znuc,vid,viu,ev)

   implicit none

   !Arguments
   integer, intent(in) :: iter,iorb
   integer, intent(in) :: norb,nr
   integer, dimension(norb) :: no,lo
   real(kind=8), dimension(norb) :: so,ev
   real(kind=8), dimension(nr) :: v,ar,br,r,rab,vid,viu
   real(kind=8), intent(in) :: znuc

   !Local variables
   !> Maximum number of integration step
   integer, parameter :: max_itmax=100
   real(kind=8), parameter :: ai=2.d0*137.0360411d0
   real(kind=8), parameter :: ai2 = ai * ai
   !> Tolerances
   real(kind=8), parameter :: etol=-1.d-7, tol = 1.0d-14
   real(kind=8), dimension(nr) :: rabkar,rabai, fa,fb
   real(kind=8), dimension(5) :: rs
   real(kind=8) :: a1,a2,abc1,abc2,abc3,abc4,abc5,alf,amc0,amc1,amc2,amc3,amc4,arc,arin,arout,arp,arpin
   real(kind=8) :: az,b0,b1,b2,brc,brp,dev,emax,emin,evold,evv,arpout,evvai2,expzer,factor,faj,fbj,s,temp,vzero
   integer :: icount,istop,itmax,j,juflow,ka,ll,nctp,ninf,nodes

   !------Machine dependent parameter-
   !------Require exp(-2*expzer) to be within the range of the machine
   ! IBM
   expzer = 3.7D2
   !Iris     expzer =3.7E2
   !Apollo   expzer = 3.7E2
   !Sun      expzer = 3.7D2
   !Vax      expzer = 44.d0
   !ray      expzer = 2.8E3

   ! for numerical stability:
   expzer = expzer/2.d0


   ! integration coefficients
   abc1 = 1901.d0/720.d0
   abc2 = -1387.d0/360.d0
   abc3 = 109.d0/30.d0
   abc4 = -637.d0/360.d0
   abc5 = 251.d0/720.d0
   amc0 = 251.d0/720.d0
   amc1 = 323.d0/360.d0
   amc2 = -11.d0/30.d0
   amc3 = 53.d0/360.d0
   amc4 = -19.d0/720.d0

   az = znuc/(2.d0*ai)
   ka = lo(iorb)+1
   if (so(iorb) < 0.1d0 .and. lo(iorb) /= 0) ka=-lo(iorb)

   !determine effective charge and vzero for startup of
   !outward integration
   !ar = r**s * (1  + a1 r + a2 r**2 + ... )
   !br = r**s * (b0 + b1 r + b2 r**2 + ... )
   !s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
   !an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
         !/ (n ai (2 s + n))
   !bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))

   s = sqrt(real(ka*ka,kind=8)-az*az)
   if (ka > 0) then
      b0 = -az/(s+real(ka,kind=8))
   else
      b0 = (s-real(ka,kind=8))/az
   end if
   if (so(iorb) < 0.1d0) then
      vzero=vid(2)
   else
      vzero=viu(2)
   end if

   !Loop data calculated only once.
   !Set ar() and br() to zero.

   do j=1,nr
      ar(j) = 0.0d0
      br(j) = 0.0d0
   end do
   do j=2,nr
      rabkar(j)=rab(j)*real(ka,kind=8)/r(j)
   end do
   do j=2,nr
      rabai(j)=rab(j)/ai
   end do
   do j=2,5
      rs(j)=r(j)**s
   end do

   ! set the underflow trap
   juflow=1
   do j=2,nr
      if (s*abs(log(r(j))) >= expzer/2.d0) juflow = j
   end do


   emax = 0.0d0
   emin = -100000.0d0
   if (ev(iorb) > emax) ev(iorb) = emax

   !Main loop for the integration
   main_loop: do itmax=1,max_itmax+1
      if (itmax >= max_itmax) write(6,'(" iorb =",i3," iter =",i3," ev =",1pe18.10," nodes =",i4)') iorb,iter,ev(iorb),nodes
      !Too many iterations !!
      if (itmax > max_itmax) return

      if (ev(iorb) > 0.d0) then
        write(6,'(//," error in difrel - ev(",i2,") greater then v(infinty)")') iorb
        stop 'difrel one'
      end if

      !Find practical infinity ninf and classical turning
      !point nctp for orbital.
      do icount=0,100
         do j=nr,2,-1
            temp = v(j) - ev(iorb)
            if (temp < 0.0d0) temp = 0.0d0
            if (r(j)*sqrt(temp) < expzer) exit
         end do

         ninf=j
         nctp = ninf - 5
         do j=2,ninf-5
            if (v(j) < ev(iorb)) nctp = j
         end do
         if (ev(iorb) >= etol*100d0) nctp=ninf-5
         if (ev(iorb) >= etol) ev(iorb)=0.0d0
         !Find it!
         if (nctp > 6) exit
         !Decrease ev(iorb) and start again
         ev(iorb) = 0.9d0*ev(iorb)
      end do
      if (icount > 100) then
          write(6,'(//,a)') 'error in difrel - cannot find classical'
          write(6,'(a,i2)') 'turning point in orbital ',iorb
          stop 'difrel two'
      end if

      !Outward integration from 1 to nctp, startup.

      a1 = (az*(vzero-ev(iorb))-(s+1.d0+real(ka,kind=8))*(vzero-ev(iorb)-ai2)*b0)  &
         / (ai*(2.d0*s+1.d0))
      b1 = ((vzero-ev(iorb))-2*znuc*a1) / (ai*(s+1.d0+real(ka,kind=8)))
      a2 = (az*(vzero-ev(iorb))*a1-(s+2.d0+real(ka,kind=8))*(vzero-ev(iorb)-ai2)*b1)  &
         / (2.d0*ai*(2.d0*s+2))
      b2 = ((vzero-ev(iorb))*a1-2*znuc*a2) / (ai*(s+2.d0+real(ka,kind=8)))
      do j=2,5
         ar(j) = rs(j) * (1 +(a1+a2*r(j))*r(j))
         br(j) = rs(j) * (b0+(b1+b2*r(j))*r(j))
      end do
      fa(1) = 0.0d0
      fb(1) = 0.0d0
      fa(2) = rabkar(2)*ar(2)+(ev(iorb)-v(2)+ai2)*br(2)*rabai(2)
      fb(2) = -rabkar(2)*br(2)-(ev(iorb)-v(2))*ar(2)*rabai(2)
      fa(3) = rabkar(3)*ar(3)+(ev(iorb)-v(3)+ai2)*br(3)*rabai(3)
      fb(3) = -rabkar(3)*br(3)-(ev(iorb)-v(3))*ar(3)*rabai(3)
      fa(4) = rabkar(4)*ar(4)+(ev(iorb)-v(4)+ai2)*br(4)*rabai(4)
      fb(4) = -rabkar(4)*br(4)-(ev(iorb)-v(4))*ar(4)*rabai(4)
      fa(5) = rabkar(5)*ar(5)+(ev(iorb)-v(5)+ai2)*br(5)*rabai(5)
      fb(5) = -rabkar(5)*br(5)-(ev(iorb)-v(5))*ar(5)*rabai(5)

      !Integration loop.

      nodes = 0
      do j=6,nctp

         ! Predictor (Adams-Bashforth).
         evvai2=ev(iorb)-v(j)+ai2
         evv=ev(iorb)-v(j)
         arp = ar(j-1) + abc1*fa(j-1)+abc2*fa(j-2)+abc3*fa(j-3)  &
             + abc4*fa(j-4)+abc5*fa(j-5)
         brp = br(j-1) + abc1*fb(j-1)+abc2*fb(j-2)+abc3*fb(j-3)  &
             + abc4*fb(j-4)+abc5*fb(j-5)
         fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
         fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)

         ! Corrector (Adams-Moulton).
         arc = ar(j-1) + amc0*fa(j)+amc1*fa(j-1)+amc2*fa(j-2)  &
             + amc3*fa(j-3)+amc4*fa(j-4)
         brc = br(j-1) + amc0*fb(j)+amc1*fb(j-1)+amc2*fb(j-2)  &
             + amc3*fb(j-3)+amc4*fb(j-4)
         faj = rabkar(j)*arc+evvai2*brc*rabai(j)
         fbj = -rabkar(j)*brc-evv*arc*rabai(j)

         !  Error reduction step.
         ar(j) = arc + amc0*(faj-fa(j))
         br(j) = brc + amc0*(fbj-fb(j))
         fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
         fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)

         !  Count nodes - if no underflow.
        if (j>juflow.and.ar(j)*ar(j-1)<0.d0) nodes=nodes+1

      end do

       arout = ar(nctp)
       arpout = fa(nctp)

      ! End outward integration.
      ! If number of nodes correct, start inward integration
      ! else modify energy stepwise and try again.
      if (nodes /= no(iorb)-lo(iorb)-1) then
         !  too many nodes decrease ev
         if (nodes > no(iorb)-lo(iorb)-1) then
            if (ev(iorb) < emax) emax = ev(iorb)
            ev(iorb) = ev(iorb) + ev(iorb)/10.d0
         ! too few nodes increase ev
         else
            if (ev(iorb) > emin) emin = ev(iorb)
            ev(iorb) = ev(iorb) - ev(iorb)/10.d0
         end if
         !New integration step
         cycle main_loop
      end if

      !Inward integration from ninf to nctp startup.

      do j=ninf,ninf-4,-1
         alf = v(j) - ev(iorb)
         if (alf < 0.d0) alf = 0.d0
         alf = sqrt(alf)
         ar(j) = exp(-alf*r(j))
         br(j) = ai*(alf+real(ka,kind=8)/r(j))*ar(j)/(v(j)-ev(iorb)-ai2)
      end do
      fa(ninf) = rabkar(ninf)*ar(ninf)+  &
                (ev(iorb)-v(ninf)+ai2)*br(ninf)*rabai(ninf)
      fb(ninf) = -rabkar(ninf)*br(ninf)  &
                -(ev(iorb)-v(ninf))*ar(ninf)*rabai(ninf)
      fa(ninf-1) = rabkar(ninf-1)*ar(ninf-1)+  &
                (ev(iorb)-v(ninf-1)+ai2)*br(ninf-1)*rabai(ninf-1)
      fb(ninf-1) = -rabkar(ninf-1)*br(ninf-1)  &
                -(ev(iorb)-v(ninf-1))*ar(ninf-1)*rabai(ninf-1)
      fa(ninf-2) = rabkar(ninf-2)*ar(ninf-2)  &
                +(ev(iorb)-v(ninf-2)+ai2)*br(ninf-2)*rabai(ninf-2)
      fb(ninf-2) = -rabkar(ninf-2)*br(ninf-2)  &
                -(ev(iorb)-v(ninf-2))*ar(ninf-2)*rabai(ninf-2)
      fa(ninf-3) = rabkar(ninf-3)*ar(ninf-3)  &
                +(ev(iorb)-v(ninf-3)+ai2)*br(ninf-3)*rabai(ninf-3)
      fb(ninf-3) = -rabkar(ninf-3)*br(ninf-3)  &
                -(ev(iorb)-v(ninf-3))*ar(ninf-3)*rabai(ninf-3)
      fa(ninf-4) = rabkar(ninf-4)*ar(ninf-4)  &
                +(ev(iorb)-v(ninf-4)+ai2)*br(ninf-4)*rabai(ninf-4)
      fb(ninf-4) = -rabkar(ninf-4)*br(ninf-4)  &
                -(ev(iorb)-v(ninf-4))*ar(ninf-4)*rabai(ninf-4)

      !Integration loop.
      istop = ninf-nctp

      if (istop >= 5) then
         do j=ninf-5,nctp,-1

            !Predictor (Adams-Bashforth).
            evvai2=ev(iorb)-v(j)+ai2
            evv=ev(iorb)-v(j)
            arp = ar(j+1)-(abc1*fa(j+1)+abc2*fa(j+2)+abc3*fa(j+3)  &
             +abc4*fa(j+4)+abc5*fa(j+5))
            brp = br(j+1)-(abc1*fb(j+1)+abc2*fb(j+2)+abc3*fb(j+3)  &
             +abc4*fb(j+4)+abc5*fb(j+5))
            fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
            fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)

            !Corrector (Adams-Moulton).
            arc = ar(j+1)-(amc0*fa(j)+amc1*fa(j+1)+amc2*fa(j+2)  &
             +amc3*fa(j+3)+amc4*fa(j+4))
            brc = br(j+1)-(amc0*fb(j)+amc1*fb(j+1)+amc2*fb(j+2)  &
             +amc3*fb(j+3)+amc4*fb(j+4))
            faj = rabkar(j)*arc+evvai2*brc*rabai(j)
            fbj = -rabkar(j)*brc-evv*arc*rabai(j)

            !Error reduction step.
            ar(j) = arc + amc0*(faj-fa(j))
            br(j) = brc + amc0*(fbj-fb(j))
            fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
            fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
         end do
      end if

      arin = ar(nctp)
      arpin = fa(nctp)

      !End inward integration
      !Rescale ar and br outside nctp to match ar(nctp) from
      !outward integration.

      factor = arout/arin
      do j=nctp,ninf
         ar(j) = factor * ar(j)
         br(j) = factor * br(j)
      end do
      arpin = factor * arpin

      !Find the normalizing factor.

      factor = 0.0d0
      ll = 4
      do j=2,ninf
         factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
         ll = 6 - ll
      end do
      factor = factor / 3.d0

      !Modify the eigenvalue ev.

      dev = arout * (arpout-arpin) / (factor * rab(nctp))
      if (5.d0*abs(dev) > -ev(iorb)) dev=dsign(ev(iorb),dev)/5.d0
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) > emax) then
        ev(iorb) = (evold + emax) / 2.d0
      elseif (ev(iorb) < emin) then
        ev(iorb) = (evold + emin) / 2.d0
      end if
      if (abs(dev) <= tol*(1.d0-ev(iorb))) exit main_loop

   end do main_loop !end of the integration loop


   !Normalize the wavefunction.
   factor = 1.d0 / sqrt(factor)
   do j=1,ninf
      ar(j) = factor*ar(j)
      br(j) = factor*br(j)
   end do

end subroutine difrel


!> orban is used to analyze and printout data about the orbital
subroutine orban(iXC,ispp,iorb,ar,br, &
      nr,r,rab, &
      lmax,norb,ncore,no,lo,so,zo, &
      znuc,zcore,cdd,cdu,cdc,dcrc,ddcrc, &
      viod,viou,vid,viu, &
      v,ev,ek,ep,rcov,rprb,nconf)

   implicit none
   !Arguments
   character(len=1), intent(out) :: ispp
   integer, intent(in) :: ncore,lmax,nconf,iXC,iorb
   integer, intent(in) :: nr,norb
   real(kind=8), dimension(nr), intent(in) :: ar,br
   real(kind=8), dimension(nr), intent(in) :: r,rab
   integer, dimension(norb), intent(in) :: no,lo
   real(kind=8), dimension(norb), intent(in) :: so,zo,ev
   real(kind=8), dimension(norb), intent(out) :: ek,ep
   real(kind=8), dimension(nr), intent(in) :: cdd,cdu,cdc,vid,viu,v
   real(kind=8), dimension(lmax,nr), intent(in) :: viod,viou
   real(kind=8), intent(in) :: znuc,zcore, rcov
   real(kind=8), intent(in) :: rprb                         !< Radius for the parabolic confinement potential
   real(kind=8), intent(out) :: dcrc, ddcrc
   !Local variables
   real(kind=8), parameter :: ai=2.d0*137.0360411d0
   real(kind=8), parameter :: ai2 = ai * ai
   real(kind=8), dimension(10) :: rzero,rextr,aextr,bextr
   !character(len=10) :: name
   character(len=30) :: plotfile,orbname
   real(kind=8) :: a1,an,ar2,arp,arpm,b1,bn,br2,dena,denb,deni,expzer
   real(kind=8) :: ddd,dddd,sa2,syswght
   real(kind=8) :: temp,toplot,tt
   real(kind=8) :: zps,vshift
   integer :: i,ierr,i90,i99,ii,ircov,isx,j,jj,ka,ll,llp,lp,nextr,ninf,npoint,nzero
   !c.hartwig
   !work-arrays for integration, and xc-potential
   real(kind=8) :: ttxup,ttxlo,cmin,crcov,dcrcov,ddcrcov
   real(kind=8), dimension(nr) :: ttx,tty,ttyp,ttypp
   real(kind=8), dimension(3*nr) :: ttw
   
   character(len=1), dimension(5) :: il
   character(len=2) :: cnum
   character(len=8) :: dateymd


   ka = lo(iorb)+1
   lp = ka
   if (so(iorb) < 0.1d0 .and. lo(iorb) /= 0) ka=-lo(iorb)

   !compute zeroes and extrema

   nzero = 0
   nextr = 0
   rzero(1) = 0.d0
   arp = br(2)
   if (ispp == 'r' .and. so(iorb) < 0.1d0) &
         arp = ka*ar(2)/r(2) + (ev(iorb) - viod(lp,2)/r(2) - vid(2) + ai2) * br(2) / ai
   if (ispp == 'r' .and. so(iorb) > 0.1d0) &
         arp = ka*ar(2)/r(2) + (ev(iorb) - viou(lp,2)/r(2) - viu(2) + ai2) * br(2) / ai

   do i=3,nr
      if (nextr >= no(iorb)-lo(iorb)) exit
      if (ar(i)*ar(i-1) <= 0.d0) then
         ! zero
         nzero = nzero + 1
         rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i)) / (ar(i)-ar(i-1))
      end if
      arpm = arp
      arp = br(i)
      if (ispp == 'r' .and. so(iorb) < 0.1d0) &
            arp = ka*ar(i)/r(i) + (ev(iorb) - viod(lp,i)/r(i) - vid(i) + ai2) * br(i) / ai
      if (ispp == 'r' .and. so(iorb) > 0.1d0) &
            arp = ka*ar(i)/r(i) + (ev(iorb) - viou(lp,i)/r(i) - viu(i) + ai2) * br(i) / ai
      if (arp*arpm > 0.d0) cycle

      ! extremum
      nextr = nextr + 1
      if ((arp-arpm) /=0.0_8) then
         rextr(nextr) = (arp*r(i-1)-arpm*r(i)) / (arp-arpm)
         aextr(nextr) = (ar(i)+ar(i-1))/2 - (arp**2+arpm**2) * (r(i)-r(i-1)) / (4*(arp-arpm))
      end if
      bextr(nextr) = br(i)
   end do

   !find orbital kinetic and potential energy
   !the potential part includes only the interaction with
   !the nuclear part
   ek(iorb) = br(1)*br(1)*rab(1)
   ep(iorb) = 0.d0
   sa2 = 0.d0
   lp = lo(iorb)+1
   llp = lo(iorb)*lp
   ll = 2
   if (2*(nr/2) == nr) ll=4
   do ii=2,nr
      i = nr-ii+2
      ar2 = ar(i)*ar(i)
      br2 = br(i)*br(i)
      deni = ar2
      if (ispp == 'r') deni=deni+br2
      ek(iorb) = ek(iorb) + ll * (br2 + ar2*llp/r(i)**2)*rab(i)
      if (so(iorb) < 0.1d0) &
         ep(iorb) = ep(iorb) + ll * deni*viod(lp,i)*rab(i)/r(i)
      if (so(iorb) > 0.1d0) &
         ep(iorb) = ep(iorb) + ll * deni*viou(lp,i)*rab(i)/r(i)
      ll = 6 - ll
      if (sa2 > 0.10d0) cycle
      sa2 = sa2 + deni*rab(i)
      if (sa2 <= 0.01d0) i99 = i
      i90 = i
   end do

   ek(iorb) = ek(iorb) / 3
   ep(iorb) = ep(iorb) / 3
   if (ispp == 'r') ek(iorb) = 0.d0

   !fourier analyze orbital

   !if (iorb < ncore) return
   !kzero = 0
   !kextr = 0
   !iextr = 1
   !delg = 0.2d0*pi/r(i90)
   !do i=1,100
   !   g = delg * (i-1)
   !   cg(i) = 0.d0
   !   if (i == 1 .and. lp /= 1) cycle
   !   ll = 4
   !   do j=2,nr
   !      rg = r(j) * g
   !      bsl = 1.d0
   !      if (i  /= 1) bsl = sin(rg) / rg
   !      if (lp == 2) bsl = (bsl - cos(rg)) / rg
   !      if (lp == 3) bsl = 3.d0 * (bsl - cos(rg)) / rg**2 -  bsl
   !      cg(i) = cg(i) + ll * r(j) * ar(j) * bsl * rab(j)
   !      ll = 6 - ll
   !   end do
   !   cg(i) = cg(i) / (6.d0*pi**2)
   !   write(6,'(2i3,3f13.6)') lo(iorb),i,g,cg(i),cg(i)*g**2
   !   if (i == 1) cycle
   !   
   !   !find extremum
   !   
   !   if (abs(cg(i)) > abs(cg(iextr))) iextr = i
   !   if (i == 2) cycle
   !   
   !   !zero
   !   
   !   if (cg(i)*cg(i-1) > 0.d0) cycle
   !   
   !   !zero found - update arrays
   !   
   !   if (i-iextr < 4) exit
   !   kzero = kzero + 1
   !   gzero(kzero) = delg*(cg(i)*(i-2)-cg(i-1)*(i-1))/(cg(i)-cg(i-1))
   !   kextr = kextr + 1
   !   cextr(kextr) = Dlog10(abs(cg(iextr)))
   !   gextr(kextr) = delg * (iextr-1)
   !   if (kextr == 5) exit
   !   iextr = i
   !end do
   !if (i-iextr >= 4 .or. kextr /= 5) then
   !   kextr = kextr + 1
   !   cextr(kextr) = Dlog10(abs(cg(iextr)))
   !   gextr(kextr) = delg * iextr
   !   printout
   !   vshift=-15.d0
   !end if
   !if (iorb < ncore) return
   !write(6,'(/" n =",i2,"  l =",i2,"  s =",f4.1)') no(iorb),lo(iorb),so(iorb)
   !write(6,'(8x,"ev =",1pe15.8,"  ek =",1pe14.8,"  ep =",1pe15.8)') (ev(iorb)-vshift)/2.,ek(iorb)/2.,ep(iorb)/2.
   !name = 'a extr    '
   !write(6,'(8x,a10,8f8.3)') name,(aextr(i),i=1,nextr)
   !name = 'b extr    '
   !if (ispp == 'r') write(6,'(8x,a10,8f8.3)') name,(bextr(i),i=1,nextr)
   !name = 'r extr    '
   !write(6,'(8x,a10,8f8.3)') name,(rextr(i),i=1,nextr)
   !name = 'r zero    '
   !write(6,'(8x,a10,8f8.3)') name,(rzero(i),i=1,nzero)
   !name = 'r 90/99 % '
   !write(6,'(8x,a10,8f8.3)') name,r(i90),r(i99)
   !name = 'c extr lg '
   !write(6,'(8x,a10,8f8.3)') name,(cextr(i),i=1,kextr)
   !name = 'g extr    '
   !write(6,'(8x,a10,8f8.3)') name,(gextr(i),i=1,kextr)
   !name = 'g zero    '
   !write(6,'(8x,a10,8f8.3)') name,(gzero(i),i=1,kzero)

   !------Machine dependent parameter-
   !------Require exp(-2*expzer) to be within the range of the machine
   !IBM
   expzer = 3.7d2
   !c.hartwig for numerical stability:
   expzer = expzer/2.d0

   !Find practical infinity ninf and classical turning
   !point nctp for orbital.
   do j=nr,2,-1
      temp = v(j) - ev(iorb)
      if (temp < 0.d0) temp = 0.d0
      if (r(j)*sqrt(temp) < expzer) exit
   end do
   ninf=j

   !compute charge at rcov + higher moments
   !spline interpolation/integration


   !some additional points for the spline
   npoint=min(ninf+5,nr)
   !charge(rcov)= int_0^rcov g^2 r^2 dr + int_0^infinity f^2 r^2 dr

   a1=0
   an=0
   b1=0
   bn=0
   isx=0

   !ALEX: Why not also set
   ttw=0d0
   ttxlo=0d0

   !int_0^rcov g^2 r^2 dr
   do i=1,npoint
      ttx(i)=r(i)
      tty(i)=ar(i)*ar(i)
      if (r(i)<=rcov) ircov=i
   end do
   if (ircov > ninf) then
      ircov=ninf
      write(6,'(1x,a)')                 'warning: ircov > ninf ! (ircov set to ninf)'
      write(6,'(1x,a,i12,1x,a,f21.16)') '---> ninf=   ',ninf,' r(ninf)=',r(ninf)
      write(6,'(1x,a,i12,1x,a,f21.16)') '---> npoints=',npoint,' r(npoint)=',r(npoint)
   end if
   call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,a1,b1,an,bn)
   if (ierr/=1) write(6,*)'SPLIFT ERROR!',ttw !stop 'spliq'
   isx=1
   ttxup=ttx(ircov)
   call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,crcov,ierr)
   if (ierr/=1) write(6,*)'SPLIQ ERROR!' !stop 'spliq'
   if (ispp == 'r') then
   !int_0^infinity f^2 r^2 dr
      cmin=0.
      do i=1,npoint
         tty(i) = br(i)*br(i)
      end do
      call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,  &
           a1,b1,an,bn)
      if (ierr/=1) write(6,'(a)') 'SPLIFT ERROR!' !stop 'spliq'
       !ttxup=ttx(ircov)
       !call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,cmin,ierr)
       !if (ierr/=1) stop 'spliq'
       !write(*,*) 'crcov+cmin:',crcov+cmin
      ttxup=ttx(ninf)
      call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,cmin,ierr)
      if (ierr/=1) write(6,'(a)') 'SPLIQ ERROR!' !stop 'spliq'
       !write(*,*) 'crcov+cmin:',crcov+cmin
      crcov=crcov+cmin
   end if

   !dcharge      = int_0^infinity (f^2+g^2) r^4 dr

   ttxup=ttx(ninf)
   do i=1,npoint
      tty(i)=ar(i)*ar(i)
      if (ispp=='r')tty(i)=tty(i)+br(i)*br(i)
      tty(i)=tty(i)*r(i)*r(i)
   end do
   call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,a1,b1,an,bn)
   if (ierr/=1) write(6,*)'SPLIFT ERROR!' !stop 'spliq'
   !ddd =  = int_0^rcov (f^2+g^2) r^4 dr
   call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttx(ircov),  &
        1,ddd,ierr)
   if (ierr/=1) write(6,*)'SPLIQ ERROR!' !stop 'spliq'
   call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,dcrcov,ierr)
   if (ierr/=1) write(6,*)'SPLIQ ERROR!' !stop 'spliq'

   !int_0^infinity (f^2+g^2) r^6 dr

   do i=1,npoint
      tty(i)=tty(i)*r(i)*r(i)
   end do
   call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,a1,b1,an,bn)
   if (ierr/=1) write(6,*)'SPLIFT ERROR!' !stop 'spliq'
   call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,ddcrcov,ierr)
   if (ierr/=1) write(6,*)'SPLIQ ERROR!' !stop 'spliq'
   !dddd =  = int_0^rcov (f^2+g^2) r^6 dr
   call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttx(ircov),  &
        1,dddd,ierr)
   if (ierr/=1) write(6,*)'SPLIQ ERROR!' !stop 'spliq'

   nextr=01

   !printout

   !if (iorb==ncore+nval) then
   !nval is not uninitialized so presume that ncore+nval = norb
   if (iorb==norb) then
      write(plotfile, '(a,i0,a)') 'ae.pot.conf.',nconf ,'.plt'
      open(unit=37,file=trim(plotfile),status='unknown')
      write(37,'(20e20.10)') r(1), 0.0d0
      do j=2,nr
         if (ispp=='r') then
            toplot = 0.5d0*(vid(j)+viu(j))
         elseif (ispp=='s') then
            toplot = 0.5d0*(vid(j)+viu(j))
         else
            toplot = vid(j)
         end if
         toplot = toplot + viod(1,j)/r(j)
         write(37,'(2e20.10)') r(j), toplot
      end do
      close(unit=37)
   end if

   vshift=-15d0
   il(1) = 's'
   il(2) = 'p'
   il(3) = 'd'
   il(4) = 'f'
   il(5) = 'g'

   if (iorb==1) then
      write(6,'(/,1x,a,f20.16)') 'rcov         = ',rcov
      if (ispp /= 'r' ) then
         write(6,'(1x,a)') 'charge(rcov) = int_0^rcov psi^2 r^2 dr'
         write(6,'(1x,a)') 'dcharge      = int_0^infinity psi^2 r^4 dr'
         write(6,'(1x,a)') 'ddcharge     = int_0^infinity psi^2 r^6 dr'
      else
         write(6,'(1x,a)') 'charge(rcov) = int_0^rcov g^2 r^2 dr   +  int_0^infinity f^2 r^2 dr'
         write(6,'(1x,a)') 'dcharge      = int_0^infinity (f^2+g^2) r^4 dr '
         write(6,'(1x,a)') 'ddcharge     = int_0^infinity (f^2+g^2) r^6 dr '
      end if
      write(6,'(/,1x,a,5x,a,4x,a,4x,a,4x,a,5x,a)') &
      'nl   s    occ','eigenvalue','charge(rcov)','dcharge','ddcharge','r extr'
   end if

   !Collect 2nd and 4th moment of the core charge density for NCC
   dcrc = dcrc+zo(iorb)* dcrcov
   ddcrc=ddcrc+zo(iorb)*ddcrcov
   write(6,'(1x,i1,a1,f5.1,f8.3,2(1pe15.7),2(1pe12.5),9(f7.2))') &
        no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),(ev(iorb)-vshift)/2.d0,crcov,dcrcov,ddcrcov,(rextr(i),i=1,nextr)
   !write(6,*) 'drcov at rcov :',ddd
   !write(6,*) 'ddrcov at rcov:',dddd
   !name = 'r extr    '
   !write(6,'(5x,a10,9f7.2)') name,(rextr(i),i=1,nextr)

   !write data to files atom.ae for pseudopotential-fit
   !only valence electrons
   !if (ispp/='r') ispp='n'
   if (iorb>ncore) then
      if (iorb==ncore+1) then
         zps=0.d0
         do jj=iorb,norb
            zps=zps+zo(jj)
         end do

         if (nconf == 0 ) then
            !write the comment lines for psppar here, as we need the zion= zps for the first configuration
            !There will be no input guess for psppar, but only some clues
            !what input variables for the fit need to be added
            open(unit=50,file='psppar',form='formatted')

            if (ispp == '') ispp='n'
            write(50,'(a,2g10.3,1x,a)') ispp, rcov, rprb, 'method, rcov and rprb'
            !old format
            !write(50,'(a,a,2g10.3,1x,a)') ispp, ' 20 2.0',rcov, rprb, 'the first line contains some input data'
            write(50,'(a)') '^- suggested header for initial guess of psppar for fitting -^'
            write(50,'(a)') 'Then add the Goedecker psuedopotential (psppar) with the first following line:'
            call date_and_time(dateymd)
            write(50,'(1x,2i4,2x,a,23x,a)') int(znuc+.1),int(zps+.1),dateymd,' zatom, zion, date (yymmdd)'
            write(50,'(a)') '-- you can download pseudopotentials from http://www.abinit.org/downloads/psp-links --'
            close(unit=50)
            syswght=1.d0
         else
            syswght=1.d-2
         end if

         ! Write the configuration into atom.??.ae
         write(cnum,'(i2.2)') nconf
         open(unit=40,file='atom.'//cnum//'.ae',form='formatted')
         write(40,'(i12,1x,f21.16,1x,a)')     norb-ncore,syswght, 'orbitals, system weight'
         write(40,'(3(f21.16),1pg10.3,1x,a)') znuc,zps,rcov,rprb, 'znuc, zpseudo, rcov, rprb'
         select case(ispp)
         case('r')
            write(40,'(a)') 'relativistic calculation'
         case('s')
            write(40,'(a)') 'spin polarized calculation'
         case default
            write(40,'(a)') 'non relativistic calculation'
         end select
         write(40,'(i10,a)') iXC, '   iXC (ABINIT-libXC)'
         write(40,*) nr,        'number of gridpoints'
         write(40,'(3(4x,a),9x,a,23x,a,4(12x,a))') &
              '#','n','l','s','z', '    eval','    charge','  dcharge','ddcharge'
      end if
      !better use formatted output here!
      !in case many plot files are written,
      !append the file names using advance='no'  &
      write(40,'(5x,2i5,6e20.12,1x)')  &
          no(iorb),lo(iorb),  &
          so(iorb),zo(iorb),  &
         (ev(iorb)-vshift)/2.,crcov,dcrcov,ddcrcov
      !:        ' # n l s z eval, charge, dcharge, ddcharge'! residue'
                                                          !^?^
      !Here we used to write the same wfn plot to both, the atom.ae
      !file and to one addtional plot file per config and orbital.
      !Better keep atom.??.ae short and clean. Either include a
      !line with the plot filename or even put all plots in two
      !files, one for core, one for other orbitals.
      !do i=1,nr
         !if (ispp=='r') then
            !write(40,*) r(i),ar(i),br(i)
         !else
            !write(40,*) r(i),ar(i)
         !end if
      !end do

   end if

   !c.chartwig:
   !save data for plots
   write(cnum,'(i2.2)') nconf
   if (ispp=='r') then
      orbname='ae.'//  &
           char(ichar('0')+no(iorb))//  &
           il(lo(iorb)+1)//  &
           char(ichar('0')+int(2*(lo(iorb)+so(iorb))))//'by2'//  &
           '.conf'//cnum//'.dat'
   elseif (ispp=='n') then
      orbname='ae.'//  &
           char(ichar('0')+no(iorb))//il(lo(iorb)+1)//  &
           '.conf'//cnum//'.dat'
   elseif (so(iorb)>0) then
      orbname='ae.'//  &
          char(ichar('0')+no(iorb))//il(lo(iorb)+1)//'.up'//  &
           '.conf'//cnum//'.dat'
   else
      orbname='ae.'//  &
          char(ichar('0')+no(iorb))//il(lo(iorb)+1)//'.down'//  &
           '.conf'//cnum//'.dat'
   end if
   !Let us create only two plotfiles and put the old plot file
   !name on a comment line instead. Do not write the plot to atom.ae.

   !if you prefer one file per orbital and config,
   !append the file name to the orbitals line in atom.ae
   !if (iorb>ncore) write(40,'(1x,a)')trim(plotfile)

   !if this was the last orbital, then close the current atom file
   !Pb: nval not initialized: presume norb (TD)
   if (iorb==norb) then
      close(unit=40)
   end if

   dena=0
   denb=0
   i=iorb
   if (i>ncore)i=i-ncore

   !old convention: One plot per orbital, including all core states
   !new convention: dump all plots of one configuration in two files
   !ae.core.orbitals.plt and ae.orbitals.plt
   !Those two files  will be read by the pseudo fitting program
   !In case plots of other configurations are intersting, those will
   !be written into separate, optional files, e.g. ae.03.orbitals.plt
   if (iorb == 1 .and. iorb <= ncore) then
      if (nconf == 0) then
         !(nconf is incremented shortly after calling this routine)
         plotfile='ae.core.orbitals.plt'
      else
         write(plotfile,'(a,i2.2,a)')'ae.',nconf,'.core.orbs.plt'
      end if
      open(unit=33,file=plotfile,status='unknown')
   else if (iorb==ncore+1) then
      if (nconf==0) then
         plotfile='ae.orbitals.plt'
      else
         write(plotfile,'(a,i2.2,a)')'ae.',nconf,'.orbs.plt'
      end if
      open(unit=33,file=plotfile,status='unknown')
   end if
   write(33,'(3a,i3,a,i3)')'# ' ,trim(orbname),  &
             '; plot me every :::',i-1,'::',i-1
   if (ispp=='r') then
      write(33,'(a)') '# r , major , minor , den(major) , den(minor)'
   else
      write(33,'(a)') '# r , psi , den(major)'
   end if
   do i=1,npoint
      dena = dena +  ar(i)*ar(i)*rab(i)
      denb = denb +  br(i)*br(i)*rab(i)
      if (ispp=='r') then
         write(33,'(20e20.10)') r(i),ar(i),br(i),dena,denb
      else
        write(33,'(20e20.10)') r(i),ar(i),dena
      end if
   end do
   !add a blank line for readability and use of gnuplots "every" keyword
   write(33,*)
   !do not close this io unit unless we write to one file per orbital
   if (iorb==ncore) then
      !we close the file "ae.core.orbitals.plt"
      close(unit=33)
   else if (iorb==norb) then
      !we close the file "ae.orbitals.plt"
      close(unit=33)
   end if

   if (iorb==norb) then
      !Addition for Nonlinear Core Corrections:
      !write out the charge density of the core for plotting and fitting.
      open(unit=33,file='ae.core.dens.plt')
      write(33,'(a)') '# plot file for all electron charges'
      if ( zcore /= 0.0d0) then
         write(33,'(a,3e15.6,a)') '#',zcore,dcrc/zcore,ddcrc/zcore,' 0th, 2nd and 4th moment of core charge'
      end if
      write(33,'(40x,a)')      '# radial charge distributions rho(r)*4pi*r**2'
      write(33,'(4(a,14x),a)') '#',' r ','core','valence','total'
      do i=1,npoint
          tt=cdu(i)+cdd(i)
          write(33,'(4e20.12)') r(i),cdc(i),tt-cdc(i),tt
      end do
      close(unit=33)
   end if

end subroutine orban
