!> @file
!! atomic program for generating and optimizing HGH pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module (pseudo) giving the dimensions for static memory just as in the old f77 code
!! this can be cleaned up later, but allocatable arrays
!! will not really save much memory
module dims
   implicit none
   integer, parameter :: norbmx=40           !< Maximal number of orbitals
   integer, parameter :: nrmax=10000         !< Maximal number of grid points
   integer, parameter :: maxdim=30           !< Maximal number of parameters
   integer, parameter :: lmx=5               !< Maximum orbital quantum number
   integer, parameter :: lpmx= 4
   integer, parameter :: noccmx=7
   integer, parameter :: nsmx=2
   integer, parameter :: ngmx=32
   integer, parameter :: nintmx=5*(ngmx+14)
end module dims


!> Module (pseudo) holding the actual pseudopotential parameters (psppar)
!! as needed for the packing and optimization routines.
module pseudovars
   use dims
   implicit none
   !> The functional (libxc), spin densities, orbital spin components
   integer :: ixc, nspol, nspin 
   !!@note  nspin is just to check for the relativistic case. 
   !!       only if nspin==2 and nspol==1 then kij should be nonzero.
   !!       atomic number, ionic charge. technically, these could be non-integer.
   real(kind=8) :: znuc, zion 
   !> The local part
   real(kind=8) :: rloc, gpot(4) 
   !> The separable part: nr of angular momentum channels, radii, coeffs
   integer :: lpx
   real(kind=8) :: r_l(4),  hsep(6, 4, 2) !< Indices: ij, l, s
   !> An experimental extra length scale, not useful
   real(kind=8) :: r_l2(4)
   !> nlcc parameters: actually, gcore(2:4) are zero in practice.
   real(kind=8) ::  rcore, gcore(4)
   !> A backup is needed for packing and shifting the pseudopotential parameters
   real(kind=8) :: orloc, ogpot(4), orl(4), orl2(4), ohsep(6, 4, 2), orcore, ogcore(4)
end module pseudovars


!> Module (pseudo) defining the penalty variables for amoeaba
module penaltyvars
   use dims
   implicit none
   !> Simple time profiling for penalty evaluation
   real(kind=8) :: time(3)
   !> All-electron reference values 
   real(kind=8) :: ev(norbmx), crcov(norbmx), dcrcov(norbmx), ddcrcov(norbmx)
   real(kind=8) :: excitae
   integer :: no(norbmx),noae(norbmx),lo(norbmx),nomax(0:lmx),nomin(0:lmx)
   real(kind=8) :: so(norbmx),zo(norbmx)
   !> Some weights for the penalty
   real(kind=8) :: wghtconf, wghtexci, wghtsoft, wghtrad, wghthij, wghtloc, wghtp0
   !> Some variables for the softness estimate via wavelets
   integer :: nhpow, nhgrid
   real(kind=8) :: hgridmin, hgridmax, ampl, crmult, frmult
   !weights for orbital properties are actually part of gatomvars
end module penaltyvars


!> Module (pseudo) giving the parameters for amoeba/packing
module ppackvars
   use dims
   implicit none

   integer, parameter :: maxpar = 42    !< Max parameters: up to 5 local,  4*7 nonlocal, 5 nlcc, 4 r_l2

   logical, dimension(maxpar) :: lpack  !< .true. if the parameters is fitted

   !> The spack array gives the names of all possible input.fitpar
   !! the meaning of the index is the same as for lpack
   !> @note
   !!   hsep index convention:
   !!   hsep 1   -    h11
   !!   hsep 2   -    h12
   !!   hsep 3   -    h22
   !!   hsep 4   -    h13
   !!   hsep 5   -    h23
   !!   hsep 6   -    h33
   character(len=6), dimension(maxpar), parameter :: spack = (/ &
      'rloc  ', 'gpot1 ', 'gpot2 ', 'gpot3 ', 'gpot4 ', &
      'rcore ', 'gcore1', 'gcore2', 'gcore3', 'gcore4', &
      'rs    ', 'hs11  ', 'hs12  ', 'hs22  ', 'hs13  ', 'hs23  ', 'hs33  ', 'rs2   ', &
      'rp    ', 'hp11  ', 'hp12  ', 'hp22  ', 'hp13  ', 'hp23  ', 'hp33  ', 'rp2   ', &
      'rd    ', 'hd11  ', 'hd12  ', 'hd22  ', 'hd13  ', 'hd23  ', 'hd33  ', 'rd2   ', &
      'rf    ', 'hf11  ', 'hf12  ', 'hf22  ', 'hf13  ', 'hf23  ', 'hf33  ', 'rf2   ' /) 
   
   integer :: nfit    !< Number of free parameters for amoeba/packing
   !> These options are rarely in use, but still supported
   logical :: avgl1,avgl2,avgl3,ortprj,litprj
   !the psppar are packed into the current vertex p. 
   !pp is the simplex consisting of several vertices, each with a penalty value y
   !let us keep that as dummy arguments
   !real(kind=8) :: pp(maxdim, maxdim+1), p(maxdim+1), yp(maxdim+1)
end module ppackvars


!> Module (pseudo) holding user input and relevant output for gatom,
!! the routine that solves the ks equations for the pseudoatom. 
module gatomvars
   use dims
   implicit none
   integer :: itertot   !< The total of scf iterations during this run
   integer :: ntime     !< The number of calls to gatom
   real(kind=8) :: rprb !< Confining potential, charge integral radius
   real(kind=8) :: rcov !< Charge integral radius
   integer :: ng        !< Gaussian basis set:
   real(kind=8) :: rij
   logical :: denbas
   real(kind=8) :: xp(0:ngmx)
   integer :: lmax      !< max angular momentum number of considered
   integer :: lcx       !< max occupied orbitals
   !> radial grid, integer weights, finite difference
   integer :: nint
   real(kind=8) :: rr(nintmx), rw(nintmx), rd(nintmx)
   !> Hartree potential from gaussians on the real space grid
   real(kind=8) :: vh((lmx+1)*((ngmx+1)*(ngmx+2))/2, (lmx+1)*((ngmx+1)*(ngmx+2))/2)
   !> Transformation from gaussian matrix elements to realspace and vice versa
   real(kind=8) :: rmt(nintmx,((ngmx+1)*(ngmx+2))/2, lmx+1), &
                    ud(nintmx,((ngmx+1)*(ngmx+2))/2, lmx+1)
   !deprecated matrix for computing the gradient.
   real(kind=8) :: rmtg(nintmx,((ngmx+1)*(ngmx+2))/2, lmx+1)
   !highest number of distinct n quantum numbers per l channel, nr of orbitals
   integer :: noccmax, norb
   real(kind=8) :: occup(noccmx,lmx,nsmx) !< occupation numbers
   real(kind=8) :: aeval(noccmx,lmx,nsmx) !< eigenvalues
   real(kind=8) :: etotal                 !< DFT energy
   !orbitals, charge integrals, moments, residues, nodes
   real(kind=8) :: psi(0:ngmx,noccmx,lmx,nsmx), psir0, chrg(noccmx,lmx,nsmx),  &
        dhrg(noccmx,lmx,nsmx),ehrg(noccmx,lmx,nsmx), &
        res(noccmx,lmx,nsmx), wfnode(noccmx,lmx,nsmx,3)
   !orbital weights, passed to skip computation of unneeded quantities
   real(kind=8) :: wght(noccmx,lmx,nsmx,8)
end module gatomvars


!> Parallel fitting of hgh pseudopotentials using a simplex downhill method.
!!  Takes into account multiple atomic references, excitation energies and
!!  softness monitored by errors in 3d wavelet transformations.
!!  Uses mpi and libxc, supports collinear spin polarization as well as 
!!  nonlinear core corrections and has a gpu accelerated version of the
!!  wavelet part.
program pseudo
   use libxcModule
   use pseudovars
   use penaltyvars
   use ppackvars
   use gatomvars
   
   implicit none
  
   real(kind=8), parameter :: fourpi = 16.d0*datan(1.d0)
   real(kind=8), parameter :: sqrt2pi = dsqrt(fourpi*0.5d0) 

   real(kind=8) :: pp(maxdim*(maxdim+1)), yp(maxdim+1)
   logical :: plotwf, mixref, energ, verbose, info, exists, ldump
   !all electron orbital plots
   real(kind=8) :: rae(nrmax), gf(nrmax,norbmx,nsmx)
   !more plotting arrays
   real(kind=8) :: psiold(nrmax,noccmx,lmx+1,nsmx)
   !local variables for reading and reformatting of psppar
   real(kind=8) :: ofdcoef(3,4), psppar(6,4,2), hso(6), havg(6)
   ! just a temporary array for excitations
   real(kind=8) :: excit(200)
   character(len=1) :: il(5),methodps,methodae
   character(len=80) :: errmsg
   character(len=7), dimension(2) :: is
   !clean up later
   character(len=35) :: fname
   character(len=10) :: tname
   character(len=520) :: string
   character(len=80) :: filename,label
   integer :: namoeb,nsmplx
   ! integer*4:: datedmy(3)
   character(len=8) :: dateymd
   integer ::  nconfpaw, nchannelspaw, npawl, pawstn, pawstl, pawstp
   real(kind=8) :: pawrcovfact
   character(len=125) :: pawstatom
   logical :: fullac,igrad
   !For timings
   real(kind=8) :: t0,t
   real(kind=8) :: a0,a,a_grd,b_grd,dh,et,ftol,gt,penref,ppc1,ppc2,ppc3,ppr1,ppr2,ppr3
   real(kind=8) :: qcore,r,r2,ra,randnr,rcovp,rmax,rnrm1,rnrm2,rnrm3,rprbp,sign1,sign2,ss,sw
   real(kind=8) :: tt,ttdiff,ttmax,ttold,ttpsi,ttrmax,zcore,zionp,znucp
   real(kind=8) :: wave
   integer :: i,ierr,ierrpp,igrid,igf,ii,iiter,iline,iorb,ipspcod,j,j1,j2
   integer :: k,l,lpj,lq,ngpot,nloc,np,nprl,nsign,ntrymax,nw
   integer :: nproc,iproc,ispin,iter,ixcpp,lw,ngrid,nocc

   ! include 'func.inc'
   include 'mpif.h'

   ! initialize some counters
   ntime=0
   itertot=0
   
   ! and time profiling
   time=0d0
   call cpu_time(t0)
   
   
   ! this is a large number passed as reference when the penalty
   ! contributions shall be written out by the penalty subroutine.
   
   penref=1d100
   
   
   ! MPI initialization
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
   
   
   ! this might not work with all compilers:
   ! generate additional output files for each
   ! process iproc > 0 connected to unit 6,
   ! which serves as standard output for process 0.
   if (iproc > 0) then
      write(label,'(a,i2.2,a)')'proc.',iproc,'.out'
      open(6,file=trim(label))
   end if
   
   ! debug files: intended to split part of the output per process       
   ! write(label,'(a,i2.2,a)')'dbg.',iproc,'.out'
   ! open(66,file=trim(label))
   write(6,'(1x,a)') '*********************************************'
   write(6,'(1x,a)') '***              pseudo_2.5               ***'
   write(6,'(1x,a)') '***              fitting of               ***'
   write(6,'(1x,a)') '***   Goedecker type pseudopotentials     ***'
   write(6,'(1x,a)') '***   last changes:    september 2013     ***'
   
   if (iproc>0) then  
      write(6,'(1x,a,i2,a,i2,a)')            &
           '***   output file for process ',iproc, ' of ', nproc,'    ***'  
   elseif (nproc>1) then
      write(6,'(1x,a,i2,a)')   &
           '***   parallel run with ', nproc, ' processes      ***'  
   end if
   write(6,'(1x,a,/)') '*********************************************'
   
   ! optional info when calls to bash are available
   write(label,'(4a,i3.3)')'echo -n  ',  &
        ' started  on `hostname`',  &
        ' at `date` "... "',  &
        ' >> hostnames.',iproc
   call system(label)
   
   ! default values:
   namoeb=0
   nsmplx=2000
   ng=20
   rij=2.0
   fullac=.false.
   denbas=.false.
   avgl1=.false.
   avgl2=.false.
   avgl3=.false.
   plotwf=.false.
   ortprj=.false.
   litprj=.false.
   energ=.false.
   igrad=.false.
   verbose=.false.
   info=.false.
   mixref=.false.
   ldump=.false.
   
   
   nconfpaw=-1   !disable paw-patches unless requested via -paw
   npawl   = 3
   nchannelspaw = 4
   pawstatom=''
   pawrcovfact=1.0_8
   
   
   ! some defaults for new input related to the simplex
   dh=0.2d0
   ntrymax=200
   ftol=1d-7
   
   ! some defaults related to wavelet transformations:
   nhpow=14
   nhgrid=0
   hgridmin=0.3d0
   hgridmax=0.5d0
   ampl=0d0
   crmult=8d0
   frmult=8d0
   
   ! do not read command line options as in previous versions.
   ! this modification is needed for a parallel run.
   ! instead, read those options from a file input.pseudo
   ! the file will be parsed for keywords line by line.
   ! some keywords are to be followed with input values.
   ! for backwards compatibility, it is allowed to omit
   ! blank characters between keyword and value.
   
   
   write(6,'(1x,a)')   'Reading the file input.pseudo'
   write(6,'(1x,a,/)') '_____________________________'
   
   
   open(unit=11,file='input.pseudo')
   iline = 0
   loop_iline: do
      iline = iline + 1
      read(11,'(a)',iostat=ierr) string
      if (ierr < 0) then
         write(6,'(i3,a)') iline-1,' lines have been read.'
         
         if (iline == 1) then
            write(6,'(3x,a,i2)')  'could not read the file input.pseudo'
            write(6,'(/,1x,a,/)') 'possible options (with input n/float):'
            write(6,'(19x,a)')    '-ng      n number of gaussians - raise this for wider sigmamax'
            write(6,'(19x,a)')    '-rij     f divider for gaussians - raise this for denser sigmamin' 
            write(6,'(19x,a)')    '-denbas    dense basis with smaller steps simgan/simgan-1'
            write(6,'(19x,a)')    '-mixref    allow mixed ae reference data, e.g. various rprb'
            write(6,'(19x,a)')    '-plot      plot wfs after each iteration'
            write(6,'(19x,a)')    '-info      write gaussian coefficients of final wavefunctions to gauss.par'

            write(6,'(/,1x,a)') 'keywords for fitting:'
            write(6,'(19x,a)')    '-cy      n number of fitting cycles (-1< n< 10000)'
            write(6,'(19x,a)')    '-fit       equivalent to -cy 1'
            write(6,'(19x,a)')    '-maxiter n max number of iterations per simplex cycle'
            write(6,'(19x,a)')    '-inwidth f initial width of the random simplex'
            write(6,'(19x,a)')    '-stuck   n max steps to consider simplex stuck'
            write(6,'(19x,a)')    '-cvg     f convergence criteria for simplex spread'
            
            write(6,'(/,1x,a)') 'keywords for backwards compatibility:'
            write(6,'(19x,a)')    '-orth      "orthogonalisation" of the hij as in ref. Krack'
            write(6,'(19x,a)')    '-lith      offdiagonal hij depend on hii as in ref. HGH'
            write(6,'(19x,a)')    '-fullacc   use max. number of gaussians'
            write(6,'(19x,a)')    '-lnso      avg nonlocal potential =0 for l=n  '
            write(6,'(19x,a)')    '           (only for the highest projector)'
            
            write(6,'(/,1x,a)') 'keywords related to wavelets for softness:'
            write(6,'(19x,a)')    '-nh      n no of grids with different spacings'
            write(6,'(19x,a)')    '-hmin    f minimum value for the grid spacing '
            write(6,'(19x,a)')    '-hmax    f maximum value for the grid spacing'
            write(6,'(19x,a)')    '-nhpow   n power for weighting softness samples'
            write(6,'(19x,a)')    '-crmult  f localization radius for scaling functions'
            write(6,'(19x,a)')    '-frmult  f localization radius for wavelets'
            write(6,'(19x,a)')    '-offset  f offset between grid and atom'
            
            write(6,'(/,1x,a)') 'keywords related to paw patches:'
            write(6,'(19x,a)')    '-pawn          no opt.,calculate  pawpatch projectors for the nth configuration '
            write(6,'(19x,a)')    '-noflpawn      pawpatch patches for the first n ls (defaults to 3)'
            write(6,'(19x,a)')    '-nchannelspawn set number of paw projectors to n (defaults to 4)'
            write(6,'(19x,a)')    '-pawstatomname file named name will be read for initial wavefunction'
            write(6,'(19x,a)')    '-pawstnn       initial wave function has first quantum number n ' 
            write(6,'(19x,a)')    '-pawstll       initial wave function has angular momentum l' 
            write(6,'(19x,a)')    '-pawstpp       initial wave function is multiplied by r**p' 
            write(6,'(19x,a)')    '-pawrcovfactf  rbox for paw is equal to rcov*pawrcovfact. defaults to 1' 
         end if
         exit
      else if (ierr > 0) then
         write(6,'(a,i0,a)') 'line ', iline, ' skipped, reading error.'
         cycle
      end if

      ! Fix the number of cycles for the fit
      ii=index(string,'-cy')
      if (ii /= 0) then
         label=string(ii+3:min(ii+13,120))
         read(label,*,iostat=ierr) namoeb
         write(6,'(1x,i0,1x,a)') namoeb, 'fit cycles'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      ! Fix to one the number of fit cycle.
      ii=index(string,'-fit')
      if (ii /= 0) then
         namoeb=1
         write(6,'(1x,a)') 'do one fit cycle'
      endif
      ! Orthogonalize the projectors
      ii=index(string,'-orth')
      if (ii /= 0) then
         ortprj=.true.
         write(6,'(1x,a)') 'orthogonalize the projectors'
      endif
      ii=index(string,'-lith')
      if (ii /= 0) then
         litprj=.true.
         write(6,'(1x,a)') 'transform the projectors as in literature'
      endif
      ! Fix the maximum number of iterations for the simplex.
      ii=index(string,'-maxiter')
      if (ii /= 0) then
         label=string(ii+8:min(ii+18,120))
         read(label,*,iostat=ierr)nsmplx
         write(6,'(i4,a)') nsmplx, ' max. simplex iterations'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      ii=index(string,'-ng')
      if (ii /= 0) then
         label=string(ii+3:min(ii+13,120))
         read(label,*,iostat=ierr)ng
         write(6,'(a,i3,a)')' Basis set size:',ng,' Gaussians'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      ii=index(string,'-rij')
      if (ii /= 0) then
         label=string(ii+4:min(ii+14,120))
         read(label,*,iostat=ierr)rij
         write(6,'(a,f6.3)')' smallest Gaussian: sigma0 = rloc /',rij
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      ii=index(string,'-fullacc')
      if (ii /= 0) then
         ng=ngmx
         write(6,'(1x,a)') 'use max. number of gaussians'
      endif
      ii=index(string,'-plot')
      if (ii /= 0) then
         plotwf=.true.
         write(6,'(1x,a)') 'plot wfs after each iteration'
      endif
      ii=index(string,'-denbas')
      if (ii /= 0) then
         denbas=.true.
         write(6,'(1x,a)') 'use dense Gaussian basis'
      endif
      ii=index(string,'-mixref')
      if (ii /= 0) then
         mixref=.true.
         write(6,'(1x,a)') 'allow mixed reference data'
      endif
      ii=index(string,'-info')
      if (ii /= 0) then
         info=.true.
         write(6,'(1x,a)') 'write gaussian coefficients of final'
         write(6,'(1x,a)') 'wavefunctions to gauss.par'
      endif
      ii=index(string,'-l1so')
      if (ii /= 0) then
         avgl1=.true.
         write(6,'(1x,a)') 'average nonlocal potential zero for l=1'
         write(6,'(1x,a)') '(only for the highest projector )'
      endif
      ii=index(string,'-l2so')
      if (ii /= 0) then
         avgl2=.true.
         write(6,'(1x,a)') 'average nonlocal potential zero for l=2'
         write(6,'(1x,a)') '(only for the highest projector)'
      endif
         
      ii=index(string,'-l3so')
      if (ii /= 0) then
         avgl3=.true.
         write(6,'(1x,a)') 'average nonlocal potential zero for l=3'
         write(6,'(1x,a)') '(only for the highest projector)'
      endif

      ii=index(string,'-dump')
      if (ii /= 0) then
         ldump=.true.
         write(6,'(1x,a)') 'dumpfile for ekin.test.f requested'
      endif
      
      ii=index(string,'-inwidth')
      if (ii /= 0) then
         label=string(ii+8:min(ii+18,120))
         read(label,*,iostat=ierr)dh
         write(6,'(f12.4,a)')dh,' initial width of the simplex'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',  &
              'its default due to a reading error. check line',iline
      endif
      ii=index(string,'-stuck')
      if (ii /= 0) then
         label=string(ii+6:min(ii+16,120))
         read(label,*,iostat=ierr)ntrymax
         write(6,'(i4,a)')ntrymax,' max simplex steps when stuck'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',  &
              'its default due to a reading error. check line',iline
      endif
      ii=index(string,'-cvg')
      if (ii /= 0) then
         label=string(ii+4:min(ii+14,120))
         read(label,*,iostat=ierr)ftol
         write(6,'(1pe12.4,a)') ftol,' convergence criteria for the simplex'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',  &
              'its default due to a reading error. check line',iline
      endif
      
      
      ! more options related to and only needed with wavelet
      ! transformations, i.e. for positive weight on softness
      ii=index(string,'-nh')
      if (ii /= 0) then
         label=string(ii+3:min(ii+13,120))
         read(label,*,iostat=ierr)nhgrid
         write(6,'(i4,a)') nhgrid, ' samples of grid spacings for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      
      ii=index(string,'-hmin')
      if (ii /= 0) then
         label=string(ii+5:min(ii+15,120))
         read(label,*,iostat=ierr)hgridmin
         write(6,'(f12.4,a)') hgridmin, ' min grid spacing for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      
      ii=index(string,'-hmax')
      if (ii /= 0) then
         label=string(ii+5:min(ii+15,120))
         read(label,*,iostat=ierr)hgridmax
         write(6,'(f12.4,a)') hgridmax, ' max grid spacing for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      
      
      ii=index(string,'-nhpow')
      if (ii /= 0) then
         label=string(ii+6:min(ii+16,120))
         read(label,*,iostat=ierr)nhpow
         write(6,'(i4,a)') nhpow, ' power weighting for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      
      ii=index(string,'-offset')
      if (ii /= 0) then
         label=string(ii+6:min(ii+16,120))
         read(label,*,iostat=ierr)ampl
         write(6,'(f12.4,a)') ampl, ' offset for grid for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      
      ii=index(string,'-crmult')
      if (ii /= 0) then
         label=string(ii+7:min(ii+17,120))
         read(label,*,iostat=ierr)crmult
         write(6,'(f12.4,a)') crmult, ' coarse cutoff for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
      
      ii=index(string,'-frmult')
      if (ii /= 0) then
         label=string(ii+7:min(ii+17,120))
         read(label,*,iostat=ierr)frmult
         write(6,'(f12.4,a)') frmult, ' fine cutoff for wvlt'
         if (ierr /= 0) write(6,'(a,a,i3)')'above value was set to ',   &
              'its default due to a reading error. check line',iline
      endif
     
      !keywords for paw
      ii=index(string,'-paw')
      if (ii /= 0) then
         label=string(ii+4:min(ii+12,520))
         read(label,*) nconfpaw
         write(6,*)  'will calculate pawpatches for conf no ', nconfpaw
      endif
      ii=index(string,'-noflpaw')
      if (ii /= 0) then
         label=string(ii+8:min(ii+16,520))
         read(label,*) npawl
         write(6,*)  'will calculate paw patches for the', npawl, ' first ls '
         
      endif
      
      ii=index(string,'-nchannelspaw')
      if (ii /= 0) then
         label=string(ii+13:min(ii+21,520))
         read(label,*)  nchannelspaw
         write(6,*)  'will consider', nchannelspaw , ' paw channels  '
         
      endif
      
      ii=index(string,'-pawstatom')
      if (ii /= 0) then
         pawstatom=trim(string(ii+10:min(ii+130,520)))
         ii=index(pawstatom,' ')
         if (ii /= 0) then
            pawstatom=trim(pawstatom(:ii-1))
         endif
         write(6,*)  'will consider  ', trim(pawstatom) ,'file for reading the initial potential'
      endif
      
      ii=index(string,'-pawstn')
      if (ii /= 0) then
         label=string(ii+7:min(ii+15,520))
         read(label,*)  pawstn
         write(6,*)  ' n of st. wf. is ', pawstn
      endif
      
      ii=index(string,'-pawstl')
      if (ii /= 0) then
         label=string(ii+7:min(ii+15,520))
         read(label,*)  pawstl
         write(6,*)  ' l of st. wf. is ' , pawstl
      endif
      
      ii=index(string,'-pawstp')
      if (ii /= 0) then
         label=string(ii+7:min(ii+15,520))
         read(label,*)  pawstp
         write(6,*)  ' initial wf radial part is  multiplied by r**' , pawstp
      endif
      
      ii=index(string,'-pawrcovfact')
      if (ii /= 0) then
         label=string(ii+12:min(ii+20,520))
         print *, label
         read(label,*)  pawrcovfact
         write(6,*)  ' rbox is rcov   multiplied by ' , pawrcovfact
      endif
      
   end do loop_iline ! Loop over input lines from input.pseudo ends here
   close(unit=11)
   
   
   ! further output about input variables
   if (nconfpaw /= -1) then
      namoeb=0
      write(6,'(1x,a)') ' fitting disactivated because paw option is active'
   endif
   
   if (namoeb == 0) then
      write(6,'(1x,a)') 'do one pseudopotential calculation.'
      write(6,'(1x,a)') 'no fitting.'
   endif
   
   if (ortprj .and. litprj ) then
      write(6,'(1x,a)') 'use only one option -orth or -lith!'
      ! stop
   endif
   
   
   ! a little estimation on memory reqs with wvlt
   if (nhgrid>0) write(6,'(a,f5.2,a)')  &
        'the wavelet coeffs will use about ', 8.0/2**20*int(2*crmult/hgridmin)**3,' mb per orbital.'
   
   
   ! ----------------- read data from ae calculation -----------------
   write(filename,'(a,i2.2,a)') 'atom.',iproc,'.ae'
   write(6,'(/,1x,a)') 'Reading data from '//trim(filename)
   write(6,'(1x,a,/)') '____________________________'
   inquire(file=filename, exist=exists)
   ierr=0
   if (.not.exists) then
      ierr=2
      write(6,*)
      write(6,*) 'no such file! attempt to read atom.00.ae instead.' 
      write(6,*) 'this is a waste of resources, as you could treat'
      write(6,*) 'more excitations or different parameters (rprb,'
      write(6,*) 'rloc, or even ixc) with this number of processes'
      filename='atom.00.ae'
      inquire(file=filename, exist=exists)
      if (.not.exists) then
         write(6,*)'cannot proceed. file atom.00.ae not found. '
         ierr=3
      end if
   end if
   call errorhandler(ierr,iproc,nproc,'atomic reference file missing!')

   ! Read the file atom.??.ae
   open(unit=40,file=trim(filename),form='formatted',status='unknown')
   read(40,*,iostat=ierr) norb, wghtconf
   if (ierr /= 0) ierr=3
   call errorhandler(ierr,iproc,nproc,'error: 1st line of ae ref data')
   read(40,*,iostat=ierr) znucp,zionp,rcovp,rprbp
   if (ierr /= 0) ierr=3
   call errorhandler(ierr,iproc,nproc,'error: 2nd line of ae ref data')
   read(40,'(a)',iostat=ierr) label
   if (ierr /= 0) ierr=3
   call errorhandler(ierr,iproc,nproc,'error: 3rd line of ae ref data')
   j=1
   
   do i=len(label),1,-1
      if (label(i:i) /= ' ') j=i
   end do
   methodae=label(j:j)
   read(40,'(a)',iostat=ierr) label
   j1=1
   j2=2
   do i=len(label),1,-1
      if (label(i:i) /= ' ') j1=i
   end do
   do i=len(label),j1,-1
      if (label(i:i) == ' ') j2=i
   end do
   j2=j2-1
   
   ! reading of ixc
   ! for now, only keep the two most commonly used functionals
   ! backwards compatible. otherwise, require abinits ixc < 0
   
   ! icorrp=label(j1:j2)
   if    (label(j1:j2) == 'pade') then
      ixc=-20
   elseif (label(j1:j2) == 'pbe') then   
      ixc=-101130
   else
      read(label(j1:j2),*,iostat=ierr) ixc
      if (ierr /= 0) then
         write(6,'(1x,a)') 'could not read the xc input in '//trim(filename)
         stop
      end if
   end if
   
   write(6,*)
   ! read(40,'(t2,a)',iostat=ierr) methodae
   ! read(40,'(t2,a)',iostat=ierr) icorrp
   read(40,*,iostat=ierr) ngrid
   if (ierr /= 0 .or. ngrid .gt. nrmax ) ierr=3
   if (ngrid > nrmax )  write(6,'(1x,a)') 'input number value is to large.'
   call errorhandler(ierr,iproc,nproc,'error: nr gridpoints in ae ref')

   write(6,'(/,1x,a,i10)')     'pseudo states = ', norb
   write(6,'(1x,a,f10.3)')     'znuc          = ', znucp
   write(6,'(1x,a,f10.3)')     'zpseudo       = ', zionp
   write(6,'(1x,a,1pe10.3)')   'r_covalent    = ', rcovp
   write(6,'(1x,a,1pe10.3)')   'r_confining   = ', rprbp
   write(6,'(1x,a,a10)')       'method        = ', methodae
   write(6,'(1x,a,i10)')       'gridpoints    = ', ngrid

   il(1) = 's'
   il(2) = 'p'
   il(3) = 'd'
   il(4) = 'f'
   il(5) = 'g'
   nspin=1
   is(1) = '  so=0'
   if (methodae == 'r' .or. methodae == 's') then
      nspin=2
      is(1)= 'so=+0.5'
      is(2)= 'so=-0.5'
   endif
   ! for now, use this convention
   nspol=1
   if (methodae == 's') then
      nspol=2
   end if
   
   write(6,*)
   write(6,'(1x,a,i7,a,i5)') 'Initializing libXC with iXC =', ixc,'; nspol =',nspol
   ! the barriers are here because the init routine writes to stout
   ! with each process. need to find a way to fix this.
   if (nproc > 1) call mpi_barrier(mpi_comm_world,ierr)  
   call libxc_functionals_init(ixc,nspol)
   if (nproc > 1) call mpi_barrier(mpi_comm_world,ierr)  
   
   
   write(6,*)
   read(40,*) !some comment line
   write(6,*)' nl    s      occ        eigenvalue     charge(rcov)    '
   
   
   ! this file should hold all ae plots
   if (plotwf.and.iproc==0) open(41,file='ae.orbitals.plt')
   ! read the ae data
   do iorb=1,norb
      read(40,*,iostat=ierr) no(iorb),lo(iorb),so(iorb),zo(iorb),  &
           ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)!,plotfile
      ! write(6,*,iostat=ierr) no(iorb),lo(iorb),so(iorb),zo(iorb), &
      ! ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb),plotfile
      if (ierr /= 0) ierr=3
      call errorhandler(ierr,iproc,nproc,'Reading error in ae ref data')
      write(6,'(1x,i1,a1,f6.1,f10.4,2f16.10,2f16.7,a)') &
           no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),  &
           ev(iorb),crcov(iorb) !& 
      ! ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)
      if (plotwf.and.iproc==0) then
         ! use this statement if atom is compiled to write one file
         ! per configuration and orbital
         ! open(41,file=trim(plotfile))
         read(41,*)
         read(41,*)
         do igrid=1,ngrid
            read(41,*,iostat=ierr) rae(igrid),(gf(igrid,iorb,igf),  &
                 igf=1,nspin)  
            ! error handling in the loop is slow, but better give detailed feedback
            ! for now, we only plot the ground state, i.e. atom.00.ae
            if (ierr /= 0) ierr=2
            ! write(errmsg,'(a,a,a,i0,a,i0)')'error reading ae plots', trim(plotfile), 'orb',iorb,'pt',igrid
            ! call errorhandler(ierr,iproc,nproc,errmsg)
         end do
         read(41,*)
         ! don't close this unit when reading from one single file
         ! close(unit=41)
      end if
   end do
   if (plotwf) close(unit=41)
   lmax=0
   lcx=0
   do iorb=1,norb
      lmax=max(lo(iorb),lmax)
      if (zo(iorb).gt.1.d-10)lcx=max(lo(iorb),lcx)
   end do
   ! print*,'lmax=',lmax
   ! print*,'lcx=',lcx, '( charge > 1.0d-10)'
   if (lmax.gt.lmx+1)ierr=3
   call errorhandler(ierr,iproc,nproc,'array dimension problem:lmax')
   ! compute corresponding n-quantum numbers of the pseudostates
   ! no()   will contain n quantum numbers of the pseudostates afterwards
   ! noae() will contain n quantum numbers of the ae-states afterwards
   ! no() starts from n=1 for each(!) l
   noccmax=0
   do l=0,lmax
      nomin(l)=100
      nomax(l)=0
   end do
   do iorb=1,norb
      nomin(lo(iorb))=min(no(iorb),nomin(lo(iorb)))
      nomax(lo(iorb))=max(no(iorb),nomax(lo(iorb)))
   end do
   do iorb=1,norb
      noae(iorb)=no(iorb)
      no(iorb)=no(iorb)-nomin(lo(iorb))+ 1
   end do
   do l=0,lmax
      noccmax= max(noccmax,nomax(l)-nomin(l)+1)
   end do
   write(6,*) 'noccmax ', noccmax
   if (noccmax.gt.noccmx)ierr=3 
   call errorhandler(ierr,iproc,nproc,'array dimension problem: noccmax')
   ! print*,'noccmax=',noccmax
   do nocc=1,noccmax
      do l=0,lmax
         do ispin=1,nspin
            occup(nocc,l+1,ispin)=0.0d0
            aeval(nocc,l+1,ispin)=0.0d0
            chrg (nocc,l+1,ispin)=0.0d0
            dhrg (nocc,l+1,ispin)=0.0d0
            ehrg (nocc,l+1,ispin)=0.0d0
         end do
      end do
   end do
   !why not initialize all elements of occup with zeroes?
   occup=0d0
   do iorb=1,norb
      nocc=no(iorb)
      l=lo(iorb)
      ispin=1
      if (so(iorb).lt.0) ispin=2
      occup(nocc,l+1,ispin)=zo(iorb)
      do j=1,ngrid
         psiold(j,nocc,l+1,ispin) = 0.0d0
         ! use major comp. as reference
         if (rae(j) /= 0.0)  &
              psiold(j,nocc,l+1,ispin)=gf(j,iorb,1)/rae(j)
      end do
   end do
   
   write(6,'(/,1x,a)') 'All electron and pseudo-wfn quantum numbers'
   write(6,'(9x,a)')   'n(AE)   l   inl(PS)'
   do iorb=1,norb
      write(6,'(6x,3i6)')  noae(iorb),lo(iorb), no(iorb)
   end do
   
   ! read excitation energies from the last line of atom.??.ae
   
   if (nproc>1) then
      write(6,*)
      ! read etotal and exchange data with all processes
      ! it would be enough to broadcast etot of system 0,
      ! but this will not take much time and allow some 
      ! output and proofreading.
      excit=0d0
      read(40,*,iostat=ierr)excit(1)
      write(6,*)'all electron energy (ryd):',excit(1)
      if (ierr /= 0) then
         write(6,*)
         write(6,*)'               warning'
         write(6,*)'the last line of the atomic reference file'
         write(6,*)'must specify the total energy.'
      end if
      excit(nproc+1:2*nproc)=excit(1)
      call mpi_alltoall(  &
           excit(nproc+1:2*nproc),1,mpi_double_precision,   &
           excit(1:nproc),1,mpi_double_precision,  &
           mpi_comm_world,ierr)
      
      if (ierr /= 0) write(6,*)'           warning: mpi error'        
      ! total energies to excitation energies
      ! different conventions of etot in gatom and atom.f
      excit=0.5d0*(excit-excit(1))
      !write(6,*)'debug: no factor two in excitation energies'
      !excit=(excit-excit(1))
      excitae=excit(iproc+1)
      if (ierr /= 0)excitae=0d0
      write(6,*)
      write(6,*)'excitation energies (ae)'
      do i=0,nproc-1
         write(6,'(1x,i3,3(1x,e20.10))')  &
              i, excit(i+1)
      end do
      ierr=0
      if (excitae==0d0.and.iproc/=0)ierr=1
      call errorhandler(ierr,iproc,nproc,'excitation energy is zero and thus ignored')
   else
      excitae = 0d0
   end if
   ! done reading atomic ref data
   close(unit=40)
   
   ! weights will be read from input.weights
   !
   
   ! pseudo 2.4 was backwards compatible with this files
   ! reading conventions from pseudo2.2 and 2.3.
   ! for this version, this is not the case anymore!
   
   write(6,'(/,1x,a)') 'Reading penalty weights from file input.weights'
   write(6,'(1x,a,/)') '_______________________________________________'

   inquire(file='input.weights',exist=exists) 
   if (.not.exists) then
      write(6,'(1x,a)') 'The file input.weights is lacking.'
      write(6,'(1x,a)') 'This file is generated by the program atom.'
      if (nproc > 1) call mpi_finalize(ierr)
      stop
   end if
   open(unit=24,file='input.weights',form='formatted')
   read(24,*)
   wghtp0   = 0d0
   wghtsoft = 0d0
   wghtrad  = 0d0
   wghthij  = 0d0
   wghtexci = 0d0
   ! what about a different ordering for reading these?
   read(24,*,iostat=ierr) wghtp0,wghtsoft,wghtrad,wghthij,wghtloc,wghtexci
   ! no need to call error handler here, shared input file
   if (ierr /= 0) write(6,*) 'Reading error for weights of psi(r=0) and softness.'
   write(6,'(a,1pe10.3)') ' weight for psi(r=0)=0 is     ',wghtp0
   write(6,'(a,1pe10.3)') ' for ekin(gauss-wavelet)      ',wghtsoft
   write(6,'(a,1pe10.3)') ' for keeping radii wide       ',wghtrad
   write(6,'(a,1pe10.3)') ' for keeping hij small        ',wghthij
   write(6,'(a,1pe10.3)') ' for keeping vloc local       ',wghtloc
   write(6,'(a,1pe10.3)') ' and for excitation energies  ',wghtexci
   
   read(24,*) !comment line
   
   ! read the weights for eigenvalues and integrals into the array wght
   do iorb=1,norb
      nocc=no(iorb)
      l=lo(iorb)
      ispin=1
      ss = so(iorb)
      if (ss.lt.0) ispin=2
      read(24,*) nw,lw,sw,(wght(nocc,l+1,ispin,i),i=1,8)
      if (noae(iorb) /= nw .or. l /= lw .or. ss /= sw) then
         write(6,*) 'error in file input.weights'
         write(6,*) 'need weights for n,l,s:', noae(iorb),l,so(iorb)
         write(6,*) 'found            n,l,s:',nw,lw,sw
         if (nproc>1) call mpi_finalize(ierr)
         stop
      endif
   end do
   close(unit=24)
   
   !Open the file to dump the vertex points
   open(unit=99,file='vertex.dump')


   ! ---------------------------------------------------------------
   ! main loop begins here
   ! ---------------------------------------------------------------
   
   do iiter=1,max(namoeb,1)
      
      ! if the excitation energy has a nonzero weight
      ! we always need to compute the total energy. 
      ! the ground state energy is needed as a reference
      ! for the excitation  energies, thus iproc==zero
      energ=(wghtexci>0d0.or.(iproc==0.and.nproc/=1))
      
      
      ! read initial pseudopotential parameter from psppar
      ! test if parameter from atom.ae and psppar are consistent
      ! notice that fitting to inconsistent ae data might be desired
      ! in a parellel run, for instance several values for rprb.
      if (iiter == 1) then
         
         ! data for nlcc used to be read from a separate file "nlcc"... 
         ! the convention was changed to read from psppar only, using pspcod = 11
         ! i leave this block as comment in case one wants to try a
         ! more general analytic form of the core charge,
         ! e.g. a gaussian times a polynomial.
         
         !!! ! the conventions for the polynomials in bigdft and pseudo
         !!! ! are slightly different. values are checked for consistency.
         !!! 
         !!! ! initial values - negative radius means no nlcc is used 
         !!! 
         !!! rcore=-1d0
         !!! gcore(1:4)=0d0
         !!! rcorepp=-1d0
         !!! gcorepp(1:4)=0d0
         !!! 
         !!! ! read optional input file as used by bigdft
         !!! open(12,file='nlcc')
         !!! read(12,*,iostat=ierr)
         !!! read(12,*,iostat=ierr)
         !!! read(12,*,iostat=ierrnlcc)rcore,gcore
         !!! close(unit=12)
         !!! if (ierrnlcc==0) then
         !!! write(6,*)
         !!! write(6,*)'Reading core charge coefficients from nlcc'
         !!! write(6,*)'__________________________________________'
         !!! write(6,*)
         !!! ! pseudos convention is to scale the polynomial by rcore for fitting
         !!! ! convert coefficients to this convention here
         !!! write(6,'(a,2f16.8)')' rcore and c0      ',rcore, gcore(1)
         !!! write(6,'(a,3f16.8)')' c2 to c4          ', gcore(2:4)
         !!! do i=2,4
         !!! gcore(i)=gcore(i)*rcore**(2*i-2) 
         !!! end do
         !!! write(6,'(a,3f16.8)')' scaled by rcore:  ', gcore(2:4)
         !!! write(6,*)
         !!! end if
         
         write(6,'(/,1x,a)') 'Reading data from psppar'
         write(6,'(1x,a,/)') '________________________'
         
         
         ! the format has been changed from previous version
         ! to allow direct compatibility with abinit pseudos.
         
         ! output will always be in HGH-K format (pspcod=10)
         ! input formats are supported as gth (2) and hgh (3) and hgh-k (10)
         ! notice some additional data needs can be read from the first line.
         
         inquire(file='psppar',exist=exists) 
         if (.not.exists) then
            ! no need to call errorhandler, shared file
            write(6,'(1x,a)') 'The file psppar is lacking.'
            write(6,'(1x,a)') 'pseudo potentials are available from http://www.abinit.org/downloads/psp-links'
            if (nproc > 1) call mpi_finalize(ierr)
            stop
         end if
         
         ! First of all be sure to have nlcc disabled per default, i.e. rcore < 0d0
         gcore =  0d0
         rcore = -1d0

         open(unit=11,file='psppar',form='formatted',status='unknown')
         ! the first line is usually for comments only.
         ! pseudo uses it to write, read and compare some additional data
         
         ! read 1st line into a string
         read(11,'(a)',iostat=ierr) label
         ! then get optional data
         read(label,*,iostat=ierr) methodps , rcov, rprb
         ! the very first character (when reading in list format)
         ! defines the method:
         ! n   non relativistc
         ! r   relativistc
         ! s   spin polarized (and relativistic)
         
         ! ng and rij are the number and relative max width of the gaussians
         if (ierr /= 0) then
            write(6,'(/,1x,a)') 'Warning'
            write(6,'(1x,a)')   'The first line of psppar does not specify optional information '
            write(6,'(1x,a)')   'about the calculation type (polarized, non-/relativistic),'
            write(6,'(1x,a)')   'the confinement rprb and the integration radius rcov.'
            write(6,'(1x,a,/)') 'The values are taken from atom.ae files without proofreading.'
            ! methodps, rprb and rcov are taken from atom.00.ae 
            ! mixed reference could cause problems with the wavelet part
            ! when process 0 has a different no of spin channels
            ! than some of the other processes.
            methodps=methodae
            rprb=rprbp
            rcov=rcovp
            ! ng and rij are always taken from input.pseudo 
            ! if none were specified, the defaults are 20 2.0
            ! in earlier versions they were read from this line as well.
            ierr=0
         else
            ! read the calculation type from the label
            j=1
            do i=len(label),1,-1
               if (label(i:i) /= ' ') j=i
            end do
            methodps=label(j:j)
            ierr=0
            if (methodps /= 'r' .and. methodps /= 'n' .and. methodps /= 's') then
               write(6,'(1x,a)')   'The first non-blank character of psppar must be one of' 
               write(6,'(1x,a)')   'n: for nonrelativistic calculations'
               write(6,'(1x,a)')   'r: for relativisitc    calculations'
               write(6,'(1x,a,/)') 's: for spin polarized  calculations'
               write(6,'(1x,a,a)') 'character found:',methodps
               write(6,'(1x,a)')   'exiting.'
               call mpi_finalize(ierr)
               stop
            end if
         end if
         if (methodae /= methodps) ierr=3
         call errorhandler(ierr,iproc,nproc,'inconsistent spin treatment.')

         ! below option does not really make sense.
         ! it could actually be useful, but needs testing for
         ! allocation errors and other causes of trouble.
         
         ! if (mixref) then
         ! write(6,*) 'warning! continue program using methodps',
         ! :                 ' from psp.par'
         ! else
         ! write(6,*) 'option ''-mixref'' allows such settings'
         ! stop
         ! endif
        
         !Second line (date in third position useless)
         read(11,*,iostat=ierr) znuc, zion
         if (ierr /= 0) then
            ! no need to call error handler, shared input file
            ! thus some stop statements have not been eliminated here
            write(6,'(/,1x,a)') 'WARNING'
            write(6,'(1x,a)')   'Could not read nuclear and valence charges on the second line of psppar.'
            if (nproc>1) call mpi_finalize(ierr)
            stop
         end if
         
         
         !Third line 
         read(11,*,iostat=ierr) ipspcod, ixcpp
         if (ierr /= 0) then
            ! no need to call error handler, shared input file
            write(6,*)
            write(6,*)'             warning'
            write(6,*)'could not read psp format and ixc from'
            write(6,*)'the third line of psppar.'
            if (nproc>1) call mpi_finalize(ierr)
            stop
         end if
         
         ! for convenience: convert LDA and PBE ixc from abinit to libxc
         if (ixcpp==1) then
            write(6,*)'LDA PADE: ixc = 1 or -20 are equivalent.'
            ixcpp=-20
         end if
         if (ixcpp==11) then
            write(6,*)'PBE: ixc = 11 or -101130 are equivalent.'
            ixcpp=-101130
         end if
         
         if (ng .gt. ngmx ) then 
            write(6,*) 'gaussians: ',ng
            write(6,*) 'maximum is:',ngmx
            if (nproc>1) call mpi_finalize(ierr)
            stop
         endif
         if (noccmax.gt.ng+1) then
            write(6,*) 'noccmax>ng+1'
            write(6,*) 'ng+1,rij ',ng+1,rij
            if (nproc>1) call mpi_finalize(ierr)
            stop 
         end if
         
         
         ! already read from psppar
         ! read(23,*) rcov,rprb
         ierr=0
         if (rcov /= rcovp) then
            ierr=1
            write(6,*)'rcov from atom.ae and psppar not identical'
            write(6,*) 'atom.ae   rcov=',rcovp
            write(6,*) 'psppar    rcov=',rcov
            if (mixref) then
               write(6,*)'no excitation energies for mixed rcov'
               excitae=0d0
               rcov=rcovp
               ierr=1
            else
               ierr=3
               write(6,*) 'option ''-mixref'' allows mixed rcov'
            endif
         endif
         call errorhandler(ierr,iproc,nproc,'different integration radii found')
         
         ! this definitely needs more testing. gpu fails for some choices of rmult
         ! i.e. crmult = 2 frmult looks savest so far
         if (3*rcov > max(crmult,frmult) ) then
            ierr=1
            write(6,*)'                note'
            write(6,*)
            write(6,*)'the localization radii for the wavelet basis'
            write(6,*)'are smaller than 3 times rcov'
            write(6,*)'psppar       rcov=',rcov
            write(6,*)'input.pseudo  crmult=',crmult
            write(6,*)'input.pseudo  frmult=',frmult
            write(6,*)'if the curve dekin.orb.dat  gets bumpy,'
            write(6,*)'try raising these two radii'
         endif
         ! no error handler needed, same for all processes.
         
         
         ierr=0
         if (rprb /= rprbp) then
            write(6,*)'rprb in atomic reference differs',   &
                 ' from the value in psppar.'
            write(6,*) 'atom.ae   rprb=',rprbp
            write(6,*) 'psppar    rprb=',rprb
            
            if (mixref) then
               write(6,*)'no excitation energies for mixed rprb'
               excitae=0d0
               rprb=rprbp
               ierr=1
            else
               write(6,*) 'option ''-mixref'' allows mixed rprb'
               ierr=3
            endif
         endif
         write(errmsg,*)'different confining potential found' 
         call errorhandler(ierr,iproc,nproc,errmsg)
         
         is(1) = '  so=0'
         if (methodps == 'r'.or.methodps == 's') then
            nspin=2
            is(1)= 'so=+0.5'
            is(2)= 'so=-0.5'
         endif
         ! read(23,'(a)') label
         ! j1=1
         ! j2=2
         ! do i=len(label),1,-1
         ! if (label(i:i) /= ' ') j1=i
         ! end do
         ! do i=len(label),j1,-1
         ! if (label(i:i) == ' ') j2=i
         ! end do
         ! j2=j2-1
         ! icorr=label(j1:j2)
         ! if (icorr /= icorrp) then
         ierr=0
         if (ixc /= ixcpp) then
            write(6,*) 'contradiction in exchange correlation'
            write(6,*) 'atom.ae   ixc=',ixc
            write(6,*) 'psppar    ixc=',ixcpp
            if (mixref) then
               write(errmsg,*)'no excitation energy for this system'
               excitae=0d0
               ierr=1
            else
               write(errmsg,*) 'option ''-mixref'' allows mixed ixc'
               ierr=3
            endif
         endif
         call errorhandler(ierr,iproc,nproc,'mixed ixc in atomic reference')
         
         ! here we previously set the parameter/variables for the xc-functional(s)
         
         ! note: znuc must be consistent, mixed references make no sense here
         ! zion will be different as soon as we have ionic configurations. 
         
         ! read(11,*) znuc, zion, rloc, gpot(1),gpot(2),gpot(3),gpot(4)
         ierr=0
         if (znucp /= znuc) then
            write(6,*) 'znuc from atom.ae and psppar not identical'
            write(6,*) 'atom.ae   znuc=',znucp
            write(6,*) 'psppar    znuc=',znuc
            ierr=3
         endif
         call errorhandler(ierr,iproc,nproc,'nucleonic charge differs from ae data')
         ierr=0
         if (zionp /= zion) then
            write(6,*) 'zion from atom.ae and psppar not identical'
            write(6,*) 'atom.ae  zion=',zionp
            write(6,*) 'psppar   zion=',zion
            ierr=1
            ! zion=zionp
         endif
         call errorhandler(ierr,iproc,nproc,'valence charge differs from ae data')
         
         ! Read the pseudopotential the way bigdft does

         ! Be sure to have zeroes for undefined entries.
         psppar = 0d0
         
         if (ipspcod == 10 .or. ipspcod == 11) then
            write(6,*) 'HGH matrix format'
            ! HGH-K format: all projector elements given.
            ! dimensions explicitly defined for nonzero output.
            
            ! local part
            read(11,*) rloc,nloc,(gpot(j),j=1,nloc) 
            
            ! only read lpx from the line after the local part
            read(11,*) lpx
            
            ! lpx is here equivalent to nsep. previous versions used
            ! it as the max l quantum number, subracting one.
            ! here, 0 does not mean s only, but a purely local psppar.
            ! negative is no longer used for local, but for r_l2.
            if (lpx-1 .gt. lmx ) then
               write(6,*) 'array dimension problem: lpx,lpmx',lpx,lpmx
               if (nproc > 1) call mpi_finalize(ierr)
               stop
            end if
            ! separable part
            ! relativistic; hij are averages, kij separatons of hsep
            do l=1,lpx !l channels
               ! add some line to read r_l2 if nprl < 0
               read(11,'(a)') string
               read(string,*) r_l(l),nprl
               if (nprl > 0) then
                  read(string,*) r_l(l),nprl,  &
                       psppar(1,l,1),(psppar(j+2,l,1),  &
                       j=2,nprl)  !h_ij 1st line
                  do i=2,nprl
                     ! spin up
                     read(11,*) psppar(i,l,1),(psppar(i+j+1,l,1),  &
                          j=i+1,nprl)!remaining h_ij 
                  end do
                  ! there is an experimental feature to use two
                  ! different r_l for each l component in the spirit:
                  ! |pi(r/rl1)ylm >hij <pj(r/rl2)|  
                  ! disable r_l2, i.e set it equal to r_l
                  r_l2(l) = r_l(l)
               else
                  ! if nprl is negative, read r_l2 from the 2nd line of hij
                  nprl=-nprl
                  read(string,*) r_l(l),i,  &
                       psppar(1,l,1),(psppar(j+2,l,1),  &
                       j=2,nprl)  !h_ij 1st line
                  read(11,*)r_l2(l), psppar(2,l,1),(psppar(2+j+1,l,1),  &
                       j=2+1,nprl)!2nd line
                  if (nprl == 3) read(11,*) psppar(3,l,1)! thid line
               end if
               
               ! there are no kij the s-projector
               if (l==1) cycle
               do i=1,nprl
                  read(11,*) psppar(i,l,2),(psppar(i+j+1,l,2),  &
                       j=i+1,nprl)!all k_ij
               end do
            end do ! loop over l 

            if (ipspcod==11) then
               ! this psppar uses nlcc
               read(11,*,iostat=ierrpp) rcore, qcore 
               if (ierrpp/=0) then
                  write(6,*)' pspcod=11 implies nlcc data on the last line,'
                  write(6,*)' but rcore and qcore could not be read!' 
                  rcore= -1d0 
               else
                  !compute gcore(1) from qcore. gcore(2:4) are
                  !always zero, but we keep them for future testing.
                  gcore(1) = fourpi* qcore * (znuc-zion) / &
                       (sqrt2pi*rcore)**3
                  write(6,*)
                  write(6,*) 'nlcc data from the last line:'
                  write(6,*) '    width rcore', rcore
                  write(6,*) ' fraction qcore', qcore
                  write(6,*)
                  write(6,*) 'derived nlcc properties:'
                  write(6,*) '    core charge', (znuc-zion)*qcore
                  write(6,*) '   rhocore(r=0)', gcore(1)
                  write(6,*) '  analytic form           gaussian'
                  write(6,*) '    polynomials               none' 
                  write(6,*)
               end if
            end if
            
         elseif (ipspcod==3) then
            write(6,'(1x,a)') 'HGH diagonal format'
            ! HGH diagonal part case
            ! technically, lpx is fixed at the max value of
            lpx=4
            read(11,*) rloc,(gpot(j),j=1,4)
            read(11,*) r_l(1),psppar(1:3,1,1) 
            do l=2,4
               read(11,*) r_l(l),psppar(1:3,l,1) 
               read(11,*)        psppar(1:3,l,2) 
            end do
         elseif (ipspcod==2) then
            write(6,'(1x,a)') 'GTH format'
            ! GTh case
            ! technically, lpx is fixed at s and p
            lpx=2
            read(11,*) rloc,(gpot(j),j=1,4)
            read(11,*) r_l(1),psppar(1:2,1,1)
            read(11,*) r_l(2),psppar(1  ,2,1)
            ! for convenience, if we have no p projector:
            if (psppar(1,2,1)<1d-5) lpx=1
         else
            ! no need to call error handler, shared input file
            write(6,*)'               warning'
            write(6,*)'pspcod (1st number of 3rd line) read from' 
            write(6,*)'psppar is unknown or not supported.' 
            write(6,*)'supported are 2,3, or 10, not ',ipspcod
            if (nproc > 1) call mpi_finalize(ierr)
            stop
         end if
         
         ! done with reading psppar
         close(unit=11) 
         
         ! avoid radii equal zero, even for unused projectors. 
         do l=1,lpx
            if (r_l(l)==0d0) then
               write(6,*)'all r_l should be nonzero.'
               write(6,*)'the r_l of the ',il(l),'-projector has been'
               write(6,*)'adjusted from 0 to 1. check your psppar.'
               r_l(l)=1d0
               r_l2(l)=1d0
            end if
            if (r_l2(l)==0d0) r_l2(l) = r_l(l)
         end do
         if ( sum((r_l2-r_l)**2) > 1d-6) then
            write(6,*)
            write(6,*)'               note'
            write(6,*)'the separable part will use two length'
            write(6,*)'scales as read from psppar.'
            write(6,*)'  l             r_l            r_l2'
            do l=1,lpx
               write(6,'(i4,5x,2f15.5)')l-1,r_l(l), r_l2(l)
            end do
            write(6,*)
         end if
         
         ! and to be very sure
         ! in case lpx increses later
         r_l(lpx+1:lpmx)=1d0         
         r_l2(lpx+1:lpmx)=1d0         
         
         
         ! then rearrange pspcod into hsep accordingly and 
         ! convert hij,kij to hij(up,down)  
         
         ! pspcod  as read are packed upper diagonal row elements of
         
         ! h =  ((l+1) hup + l hdn)/(2l+1)
         ! k =      2( hup -   hdn)/(2l+1)
         
         ! we want hsep,  packed upper diagonal col elements of
         
         ! hup = h +   l/2   k
         ! hdn = h - (l+1)/2 k
         
         
         ! 1) get offdiagonal elements where needed
         ! 2) permute from row to col ordering
         ! 3) transform from h/k to up/down
         
         
         ! just to be sure no rubbish will be added 
         hsep=0d0
         
         if (ipspcod == 2) then !gth case
            ! offdiagonal elements are zero per definition.
            ! simply rearrange hij and fill zero elements
            do l=1,lpx
               hsep(1,l,1)=psppar(1,l,1)
               hsep(2,l,1)=0.0d0
               hsep(3,l,1)=psppar(2,l,1)
               hsep(4,l,1)=0.0d0
               hsep(5,l,1)=0.0d0
               hsep(6,l,1)=psppar(3,l,1)
               ! in the polarized or relativistic case,
               ! we assume all kij to be zero,
               ! i.e. same up and down projectors
               if (nspin==2) hsep(:,l,2)=hsep(:,l,1)
            end do
         elseif (ipspcod == 3) then !hgh diagonal case
            ! we need to compute the offdiagonal elements with the following coeffs
            ofdcoef(1,1)=-0.5d0*sqrt(3.d0/5.d0) !h2
            ofdcoef(2,1)=0.5d0*sqrt(5.d0/21.d0) !h4
            ofdcoef(3,1)=-0.5d0*sqrt(100.0d0/63.d0) !h5
            
            ofdcoef(1,2)=-0.5d0*sqrt(5.d0/7.d0) !h2
            ofdcoef(2,2)=1.d0/6.d0*sqrt(35.d0/11.d0) !h4
            ofdcoef(3,2)=-7.d0/3.d0*sqrt(1.d0/11.d0) !h5
            
            ofdcoef(1,3)=-0.5d0*sqrt(7.d0/9.d0) !h2
            ofdcoef(2,3)=0.5d0*sqrt(63.d0/143.d0) !h4
            ofdcoef(3,3)=-9.d0*sqrt(1.d0/143.d0) !h5
            
            ofdcoef(1,4)=0.0d0 !h2
            ofdcoef(2,4)=0.0d0 !h4
            ofdcoef(3,4)=0.0d0 !h5
            
            ! this could possibly be done in a simpler way ...
            
            do l=1,lpx
               do ispin=1,nspin
                  hsep(1,l,ispin)=psppar(1,l,ispin)
                  hsep(2,l,ispin)=psppar(2,l,ispin)*ofdcoef(1,l)
                  hsep(3,l,ispin)=psppar(2,l,ispin)
                  hsep(4,l,ispin)=psppar(3,l,ispin)*ofdcoef(2,l)
                  hsep(5,l,ispin)=psppar(3,l,ispin)*ofdcoef(3,l)
                  hsep(6,l,ispin)=psppar(3,l,ispin)
               end do
            end do
            
            ! in the nonrelativistic case, we are done.
            if (nspin==2) then
               ! in the polarized case, copy the missing s projector
               if (nspol==2)  hsep(:,1,2)=hsep(:,1,1)
               ! use psppar as a temporary array 
               psppar=hsep
               ! and then convert hij/kij to hij up/down 
               do l=2,lpx
                  ! l is index, angular momentum +one
                  ! up
                  hsep(1,l,1)=psppar(1,l,1) +.5d0*(l-1)* psppar(1,l,2)!h11
                  hsep(2,l,1)=psppar(2,l,1) +.5d0*(l-1)* psppar(2,l,2)!h12
                  hsep(3,l,1)=psppar(3,l,1) +.5d0*(l-1)* psppar(3,l,2)!h22
                  hsep(4,l,1)=psppar(4,l,1) +.5d0*(l-1)* psppar(4,l,2)!h13
                  hsep(5,l,1)=psppar(5,l,1) +.5d0*(l-1)* psppar(5,l,2)!h23
                  hsep(6,l,1)=psppar(6,l,1) +.5d0*(l-1)* psppar(6,l,2)!h33
                  ! down
                  hsep(1,l,2)=psppar(1,l,1) -.5d0* l   * psppar(1,l,2)!h11
                  hsep(2,l,2)=psppar(2,l,1) -.5d0* l   * psppar(2,l,2)!h12
                  hsep(3,l,2)=psppar(3,l,1) -.5d0* l   * psppar(3,l,2)!h22
                  hsep(4,l,2)=psppar(4,l,1) -.5d0* l   * psppar(4,l,2)!h13
                  hsep(5,l,2)=psppar(5,l,1) -.5d0* l   * psppar(5,l,2)!h23
                  hsep(6,l,2)=psppar(6,l,1) -.5d0* l   * psppar(6,l,2)!h33
               end do
            end if
         end if
         
         if (ipspcod>9) then !hgh-k or hgh case,
            ! psppar holds hij and kij in hghk convention
            ! fill hsep(up,dn) upper diagonal col by col, as needed for the fit
            
            ! for a nonrelativistic calculation, discard the kij elements
            if (nspin==1) then
               do l=1,lpx
                  hsep(1,l,1)=psppar(1,l,1) !h11
                  hsep(2,l,1)=psppar(4,l,1) !h12
                  hsep(3,l,1)=psppar(2,l,1) !h22
                  hsep(4,l,1)=psppar(5,l,1) !h13
                  hsep(5,l,1)=psppar(6,l,1) !h23
                  hsep(6,l,1)=psppar(3,l,1) !h33
               end do
            else
               ! relativistic or polarized calculation
               do l=1,lpx
                  ! l is the index, angular momentum +one
                  ! up
                  hsep(1,l,1)=psppar(1,l,1) +.5d0*(l-1)* psppar(1,l,2)!h11
                  hsep(2,l,1)=psppar(4,l,1) +.5d0*(l-1)* psppar(4,l,2)!h12
                  hsep(3,l,1)=psppar(2,l,1) +.5d0*(l-1)* psppar(2,l,2)!h22
                  hsep(4,l,1)=psppar(5,l,1) +.5d0*(l-1)* psppar(5,l,2)!h13
                  hsep(5,l,1)=psppar(6,l,1) +.5d0*(l-1)* psppar(6,l,2)!h23
                  hsep(6,l,1)=psppar(3,l,1) +.5d0*(l-1)* psppar(3,l,2)!h33
                  ! if nspol==1 and l==1
                  ! if (l==2-nspol)cycle
                  ! down
                  hsep(1,l,2)=psppar(1,l,1) -.5d0* l   * psppar(1,l,2)!h11
                  hsep(2,l,2)=psppar(4,l,1) -.5d0* l   * psppar(4,l,2)!h12
                  hsep(3,l,2)=psppar(2,l,1) -.5d0* l   * psppar(2,l,2)!h22
                  hsep(4,l,2)=psppar(5,l,1) -.5d0* l   * psppar(5,l,2)!h13
                  hsep(5,l,2)=psppar(6,l,1) -.5d0* l   * psppar(6,l,2)!h23
                  hsep(6,l,2)=psppar(3,l,1) -.5d0* l   * psppar(3,l,2)!h33
               end do
            end if
         end if
         
         
         ! output
         write(6,*)
         write(6,'(2f10.3,4x,a)') rcov,rprb,'rcov and rprb (charge integral and confinement)'
         if (methodps == 'r') then
            write(6,'(t30,a)')'relativistic calculation'
         elseif (methodps == 's') then
            write(6,'(t30,a)')'spin polarized calculation'
         else
            write(6,'(t30,a)')'non relativistic calculation'
         endif
         write(6,'(t2,i10,t30,a)')ixc ,'ixc for XC-functional'
         write(6,*) 
         write(6,*) 'local part'
         write(6,'(f5.0,f7.2,f7.3,4e11.3,t65,a)')  &
              znuc,zion,rloc,gpot(1),gpot(2),gpot(3),gpot(4),  &
              'znuc,zion, rloc, gpot() '
         if (lpx.ge.1) then
            write(6,*) 'nonlocal part in internal format'
            write(6,'(i4,t60,a)') lpx ,  &
                 'lpx, (Projectors for l=0..lpx)'
            do l=1,lpx
               write(6,*) il(l)//'-projector'
               write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l),  &
                    (hsep(i,l,1),i=1,6),'r_l(),hsep(), '//is(1)
               if (l.gt.1 .and. nspin == 2)  &
                    write(6,'(t8,6e11.3,t76,a)')  &
                    (hsep(i,l,2),i=1,6),'      hsep(), '//is(2)
            end do
         endif
      endif ! first iteration
      
      ! previous versions have read input.weights at each iteration step.
      
      ! let us sacrifice this feature for the sake of stability in
      ! parallel runs.
      
      
      ! calc. exponents of gaussians
      a0=rloc/rij
      ! take this for an crude initial fit:
      ! tt=2.d0**.4d0
      if (denbas) then
         ! fine fit:
         tt=sqrt(sqrt(2.d0))
      else
         ! normal fit:
         tt=2.d0**.3d0
      endif
      a=a0
      do i=0,ng
         ! 
         ! a=a0*tt**i
         xp(i)=.5d0/a**2
         a=a*tt
      end do
      write(6,*)
      write(6,*)'Gaussian basis'
      write(6,*)'______________'
      write(6,*)
      write(6,'(a,4e11.4)') ' amin,amax',a0,a
      write(6,'(a,t10,3(e11.4),a,2(e11.4))') ' gaussians ',  &
           xp(1),xp(2),xp(3),' .... ',xp(ng-1),xp(ng)  
      write(6,*)'gaussians:',ng
      ! set up radial grid
      nint=5*(ng+14)
      rmax=min(15.d0*rprb,120.d0)
      a_grd=a0/400.d0
      b_grd=log(rmax/a_grd)/(nint-2)
      call radgrid(nint,rr,rw,rd,a_grd,b_grd,rmax)
      write(6,'(a,t10,3(e11.4),a,2(e11.4))') ' r-grid: ',  &
           rr(1),rr(2),rr(3),' .... ',rr(nint-1),rr(nint)  
      write(6,*)'gridpoints:',nint
      write(6,*)
      call crtvh(ng,lcx,lmax,xp,vh,nint,rmt,rmtg,ud,rr)
      write(6,*)
      
     
      if (namoeb > 0) then
         
         ! some old, commmented out feature
         ! refine simplex only every 10.th step
         ! start of if-block
         ! if (mod(iter,10) == 0) then
         
         
         ! pack initial guess
         write(6,'(1x,a)')   'Reading fitting parameters from input.fitpar'
         write(6,'(1x,a,/)') '____________________________________________'
         
         verbose=.true.
         call ppack (verbose, pp(1), 'init')
         verbose=.false.
         !
         ! initial simplex
         !
         
         ! copy the packed initial vertex to the other vertices
         ! does not make much sense, though, as we copy zeroes only.
         ! keep it in case ppack.f is changed
         do i=1,nfit
            do j=1,nfit
               pp(j+i*nfit)=pp(j)
            end do
         end do
         
         ! this width is not hard coded anymore
         ! dh=0.2d0
         
         ! shift vertices number 2 to nfit+1 by random numbers in [-dh/2:dh/2]
         ! in this convention, only move the k-th component of vertex k+1
         
         ! the random numbers generated should be the same for all processes.
         ! though, let us enforce equality with mpi broadcast from process 0.
         if (iproc==0) then
            do i=1,nfit
               ! f90 intrinsic
               call random_number(randnr)
               pp(i+i*nfit)=pp(i)+dh*1.0d0*(randnr-.5d0)
               ! intel (ifc)
               ! pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand(0.0d0))-.5d0)
               ! ibm/dec/pgi
               ! pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand())-.5d0)
               ! cray
               ! pp(i+i*nfit)=pp(i)+dh*1.0d0*(ranf()-.5d0)
            end do
         end if
         call mpi_bcast(pp,nfit*(nfit+1),mpi_double_precision,0,  &
              mpi_comm_world,ierr)
         
         write(6,*)
         write(6,'(1x,a,i4)')'starting amoeba cycle',iiter
         write(6,*)          '_________________________'
         write(6,*)
         write(6,'(2a)')' penalty contributions of the initial ',  &
              'parameters and random simplex vertices'
         write(6,'(2a)')' _______________________________________',  &
              '____________________________________'
         write(6,*)
         if (nproc>1) then
            write(6,'(a)')'  amoeba   vertex    overall penalty'//  &
                 ! :               '      softness             narrow radii      '//  &
                 '      softness             psp empirical     '//  &
                 '   all configurations   all excitation       '//  &
                 'this configuration'
            write(6,'(a)')'    step   number    value          '//  &
                 '      gauss-wvlt           (hgridmax/r)**12  '//  &
                 '   e_ks, integrals      energies             '//  &
                 'and its excitation'
         else
            write(6,'(a)')'    step    vertex   penalty,       '//  &
                 ! :               '      gaus-wvlt            narrow radii      '//  &
                 '      gaus-wvlt            psp empirical     '//  &
                 '   eigenvalues etc'
         end if
         ! do i=2,nfit+1  would be the random vertices, 1 has no shift
         iter=0
         do i=1,nfit+1
            call penalty(energ,verbose,pp(1+(i-1)*nfit),yp(i),  &
                 iproc, nproc,iter,penref)
         end do
         write(99,*)'history of lowest vertices'
         
         
         write(6,*)
         write(6,'(2a)')' penalty contributions of the',  &
              ' currently lowest vertex'
         write(6,'(2a)')' _______________________________________',  &
              '_____________'
         write(6,*)
         if (nproc>1) then
            write(6,'(a)')'  amoeba   gatom     overall penalty'//  &
                 ! :               '      softness             narrow radii      '//  &
                 '      softness             psp empirical     '//  &
                 '   all configurations   all excitation       '//  &
                 'this configuration'
            write(6,'(a)')'    step   calls     value          '//  &
                 '      gauss-wvlt           (hgridmax/r)**12  '//  &
                 '   e_ks, integrals      energies             '//  &
                 'and its excitation'
         else
            write(6,'(a)')'    step    gatom    penalty,       '//  &
                 ! :               '      gaus-wvlt            narrow radii      '//  &
                 '      gaus-wvlt            psp empirical     '//  &
                 '   eigenvalues etc'
         end if
         
         
         
         ! the vertex number 1 holds the initial psp
         
         ! this call is not needed anymore, because
         ! we already wrote details of all vertices
         
         ! iter=0
         ! write(*,*)'test line, vertex 1 again'
         ! call penalty(energ,verbose,nfit,pp(1),yp(1),
         ! :          noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
         ! :          no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
         ! :          occup,aeval,chrg,dhrg,ehrg,res,wght,
         ! :          wfnode,psir0,wghtp0,
         ! :          rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,r_l2,hsep,
         ! :          vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
         ! :          avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
         ! the following lines differ from pseudo2.2
         ! :          iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,wghtloc,
         ! :          nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
         ! :          excitae,ntime,iter,itertot,penref,time)
         
         
         ! refine simplex only every 10.th step
         ! end of if-block
         ! endif
         !
         ! starting amoeba
         !
         ! ftol=1.d-7 is not hardcoded anymore
         call amoeba(pp(1), yp(1), nfit, ftol, iter, nsmplx, namoeb,&
              iproc, nproc, ntrymax, energ, verbose)
         ! write(6,*) 'finished amoeba with ',iter,'iterations'
      else
         
         ! no fitting was requested, evaluate penalty contributions once
         
         
         ! call penalty with the verbose flag to print the details
         verbose=.false.
         energ=.true.
         penref=-1d9
         
         call ppack (verbose, pp(1), 'init') 
         verbose=.true.
         
         if (nconfpaw/=-1) then
            call pawpatch(noccmax,noccmx,lmax,lmx,lpmx,nspin,nsmx,&
                 occup,aeval,&
                 rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,&
                 psi,&
                 rae,&
                 iproc,&
                 ngrid, &
                 nconfpaw, npawl, nchannelspaw, methodps, pawstatom,&
                 pawstn, pawstl, pawstp, pawrcovfact)
            
            stop " pawpatch normal stop"
         endif
         
         call penalty(energ,verbose,pp(1),yp(1),  &
              iproc,nproc, iter,penref)
         write(6,*)'resulting penalty value',yp(1)
         
         
         ! a test for debugging: calling gatom once more should not make any difference 
         ! call gatom(energ,verbose)
      endif
      ! fit or evaluate once
      
      
      ! print results
      write(6,'(/,1x,a)')                    'penalty contributions from this configuration'
      write(6,'(1x,a,/)')                    '_____________________________________________'
      write(6,'(2(tr10,a,e12.4),/)')         'psi(r=0) =',psir0,'; psi(0)*wght=',abs(psir0*wghtp0)
      write(6,'(a,t32,a,t42,a,t55,a,t64,a)') ' nl    s      occ','ae','pseudo','diff','diff*weight'
      
      do iorb=1,norb
         write(6,'(1x,i1,a1,f6.1,f10.4)') noae(iorb),il(lo(iorb)+1),so(iorb),zo(iorb)
         nocc=no(iorb)
         l=lo(iorb)
         ispin=1
         if (so(iorb).lt.0) ispin=2
         write(6,'(t10,a,t25,4(1pe12.4))') 'eigenvalue',  &
              ev(iorb),aeval(nocc,l+1,ispin), aeval(nocc,l+1,ispin)-ev(iorb),  &
              abs(wght(nocc,l+1,ispin,1)*(aeval(nocc,l+1,ispin)-ev(iorb)))
         write(6,'(t10,a,t25,4(1pe12.4))') 'charge',  &
              crcov(iorb),chrg(nocc,l+1,ispin), chrg(nocc,l+1,ispin)-crcov(iorb),  &
              abs(wght(nocc,l+1,ispin,2)*(chrg(nocc,l+1,ispin)-crcov(iorb)))
         if (wght(nocc,l+1,ispin,3) /= 0.0d0)  &
              write(6,'(t10,a,t25,4(1pe12.4))') 'dcharge',  &
              dcrcov(iorb),dhrg(nocc,l+1,ispin), 100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb)), &
              abs(wght(nocc,l+1,ispin,3))* 100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))
         if (wght(nocc,l+1,ispin,4) /= 0.0d0)  &
              write(6,'(t10,a,t25,4(1pe12.4))') 'echarge',  &
              ddcrcov(iorb),ehrg(nocc,l+1,ispin), 100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb)),  &
              abs(wght(nocc,l+1,ispin,4))* 100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))
         write(6,'(t10,a,t25,2(1pe24.4))') 'residue',  &
              res(nocc,l+1,ispin), abs(wght(nocc,l+1,ispin,5)*res(nocc,l+1,ispin))
         if (wght(nocc,l+1,ispin,6) /= 0.0d0)  &
              write(6,'(t10,a,t25,2(1pe24.4))') 'rnode',  &
              wfnode(nocc,l+1,ispin,1),  abs(wght(nocc,l+1,ispin,6)*wfnode(nocc,l+1,ispin,1))
         if (wght(nocc,l+1,ispin,7) /= 0.0d0)  &
              write(6,'(t10,a,t25,2(1pe24.4))') 'dnode',  &
              wfnode(nocc,l+1,ispin,2), abs(wght(nocc,l+1,ispin,7)*wfnode(nocc,l+1,ispin,2))
         if (wght(nocc,l+1,ispin,8) /= 0.0d0)  &
              write(6,'(t10,a,t25,2(1pe24.4),/)') 'ddnode',  &
              wfnode(nocc,l+1,ispin,3), abs(wght(nocc,l+1,ispin,8)*wfnode(nocc,l+1,ispin,3))
      end do
      write(6,'(1x,a)') 'diff for dcharg and echarge is given in (%)'
      
      ! always print out the psppar, also if no fitting was done.
      ! this is useful to convert to pspcod 10 without fitting,
      ! and to include some extra information.
      write(6,'(/,1x,a)') 'resulting pseudpotential parameters'
      write(6,'(1x,a,/)') '___________________________________'
      
      ! at this point, overwrite  the psppar in abinit pspcod=10 format
      ! do the output twice, to psppar and to the logfile(s)
      ! additional output beyond the standad format:
      ! some info about method, basis set, rcov and rprb on line 1
      ! if nlcc is used, coeffs are written after nsep on the same line
      ! if r_l2 are used, put them in front of the h12
      ! and give lpj a negative sign, usually -2.
      
      if (iproc==0) then
         open(unit=13,file='psppar')!,position='append')
         
         if (methodps=='r') then
            write(13,'(a)',advance='no') 'relativistic '
         elseif (methodps=='s') then
            write(13,'(a)',advance='no') 'spin-polarized '
         else
            write(13,'(a)',advance='no') 'nonrelativistic '
         end if
         write(13,'(2(2x,g10.3),4x,a)') rcov,rprb, 'method, rcov and rprb'
         
         call date_and_time(dateymd)
         write(13,'(1x,2i4,2x,a,23x,a)') int(znuc+.1),int(zion+.1),dateymd,' zatom, zion, date (yymmdd)'
         write( 6,'(1x,2i4,2x,a,23x,a)') int(znuc+.1),int(zion+.1),dateymd,' zatom, zion, date (yymmdd)'
         
         ! if nlcc was used, use pspcod=11 instead of 10
         if (rcore < 0d0) then
            ipspcod=10
         else
            ipspcod=11
         end if
         write(13,'(i5,i10,i2,a,17x,a)') ipspcod,ixc,lpx-1,  &
              ' 0 2002 0','pspcod, ixc, lmax, lloc, mmax, r2well'
         write( 6,'(i5,i10,i2,a,17x,a)') ipspcod,ixc,lpx-1,  &
              ' 0 2002 0','pspcod, ixc, lmax, lloc, mmax, r2well'
         
         ! determine the number of nonzero terms in the local potential
         ngpot=0
         do j=1,4
            if (gpot(j) /= 0.d0) ngpot=j
         end do
         if (ngpot == 0) then
            write(13,'(2x,f16.8,a)')rloc,' 0 rloc nloc ck (none)'
            write( 6,'(2x,f16.8,a)')rloc,' 0 rloc nloc ck (none)'
         else
            write(label,'(a,i1.1,a)')'(2x,f16.8,i3,',ngpot,'f16.8,5x,a)'
            write(13,label)rloc,ngpot,gpot(1:ngpot),  &
                 ' rloc nloc c1 .. cnloc'
            write( 6,label)rloc,ngpot,gpot(1:ngpot),  &
                 ' rloc nloc c1 .. cnloc'
         end if
         
         write(13,'(i5,38x,a)')lpx,'nsep' 
         write( 6,'(i5,38x,a)')lpx,'nsep' 
         
         ! write out the separable part
         do l=1,lpx
            ! l is a positive index, lq the l quantum number
            lq=l-1
            do i=1,6
               ! in the relativistic case and for l > s,
               ! we need to compute havg and hso from hup and hdown
               ! needed if: l > 1, nspin =2 and nspol = 1
               if ( l>1 .and. nspin > nspol ) then
                  
                  havg(i)=((lq+1)*hsep(i,l,1)+lq*hsep(i,l,2))  &
                       /(2*lq+1)
                  hso(i)=2*(hsep(i,l,1)-hsep(i,l,2))  &
                       /(2*lq+1)
               else
                  havg(i)=hsep(i,l,1)
                  hso(i)=0d0
               endif
            end do
            ! get the matrix dimension lpj = 1, 2 or 3
            lpj=1
            if (max(abs(havg(3)),abs(hso(3)))>1d-8)lpj=2
            if (max(abs(havg(6)),abs(hso(6)))>1d-8)lpj=3
            
            ! then see if the r_l2 feature was used, compare with r_l
            if ((r_l2(l)-r_l(l))**2>1d-8) then 
               ! this will be written in front of h12 
               write(label,'(f16.8)') r_l2(l)
               ! negative sign for the dimension integer 
               write(tname,'(i3)') -lpj
            else
               ! no r_l2 used -> 16 spaces 
               label='                ' 
               write(tname,'(i3)') lpj
            end if
            ! formatted output by case of lpj
            select case(lpj)
            case(1)
               ! note: r_l2 has no meaning when lpj is one.
               ! keep it, though, as we may want to add a convention later. 
               write(13,'(2x,f16.8,a,f16.8,6x,a)')  &
                    r_l(l),trim(tname), havg(1)  &
                    ,il(l)//'-projector'
               write( 6,'(2x,f16.8,i3,f16.8,6x,a)')  &
                    r_l(l),lpj, havg(1)  &
                    ,il(l)//'-projector'
               if (l.gt.1)write(13,'(21x,f16.8)')hso(1)
               if (l.gt.1)write( 6,'(21x,f16.8)')hso(1)
            case(2)
               write(13,'(2x,f16.8,a,2f16.8,6x,a)')  &
                    r_l(l),trim(tname), havg(1:2)  &
                    ,il(l)//'-projector'
               write( 6,'(2x,f16.8,a,2f16.8,6x,a)')  &
                    r_l(l),trim(tname), havg(1:2)  &
                    ,il(l)//'-projector'
               write(13,'(2x,a,35x,f16.8)')trim(label),havg(3)
               write( 6,'(2x,a,35x,f16.8)')trim(label),havg(3)
               if (l.gt.1) then
                  write(13,'(21x,2f16.8)') hso(1:2)
                  write( 6,'(21x,2f16.8)') hso(1:2)
                  write(13,'(37x, f16.8)') hso(3)
                  write( 6,'(37x, f16.8)') hso(3)
               end if
            case(3)
               write(13,'(2x,f16.8,a,3f16.8,6x,a)')  &
                    r_l(l),trim(tname),havg(1:2), havg(4)  &
                    ,il(l)//'-projector'
               write( 6,'(2x,f16.8,a,3f16.8,6x,a)')  &
                    r_l(l),trim(tname),havg(1:2), havg(4)  &
                    ,il(l)//'-projector'
               write(13,'(2x,a,35x,2f16.8)')  &
                    trim(label),havg(3),havg(5)
               write( 6,'(2x,a,35x,2f16.8)')  &
                    trim(label),havg(3),havg(5)
               write(13,'(53x,f16.8)') havg(6)
               write( 6,'(53x,f16.8)') havg(6)
               if (l.gt.1) then
                  write(13,'(21x,3f16.8)') hso(1:2),hso(4)
                  write( 6,'(21x,3f16.8)') hso(1:2),hso(4)
                  write(13,'(37x,2f16.8)') hso(3),hso(5)
                  write( 6,'(37x,2f16.8)') hso(3),hso(5)
                  write(13,'(53x, f16.8)') hso(6)
                  write( 6,'(53x, f16.8)') hso(6)
               end if
            end select
            ! dimension of hij
         end do
         ! loop over l
         if (ipspcod==11) then
            if (any(gcore(2:4)/=0d0)) write(*,*) &
                 'warning: core charge is not just a gaussian'
            !compute qcore from gcore(1)
            qcore = gcore(1) * (sqrt2pi*rcore)**3 / &
                 (fourpi  *   (znuc-zion) )
            !append nlcc line to psppar
            write(13,'(2x,f16.8,3x, f16.8,6x,a)') rcore, qcore, 'rcore, qcore (nlcc)'
            write( 6,'(2x,f16.8,3x, f16.8,6x,a)') rcore, qcore, 'rcore, qcore (nlcc)'
         end if
         close(unit=13)
         ! iproc is zero
      end if ! end of writing the psppar by process zero

      
      if (rcore>0d0) then
         write(6,*)
         write(6,*)'analytic core charge of the nlcc:', dble(znuc-zion)*qcore
         !much simpler than before. 
         !gcorepp(i)=gcore(i)*rcore**(2-2*i)
         !sqrt( 2.0d0*atan(1.0) )*(  &
              ! gcore(1)*rcore**3  &
              !+  3.0d0*gcorepp(2)*rcore**5  &
              !+ 15.0d0*gcorepp(3)*rcore**7  &
              !+105.0d0*gcorepp(4)*rcore**9)
         !fourpi = 16.d0*atan(1.d0)
         tt=0d0
         do k= 1,nint
            r2=(rr(k)/rcore)**2
            tt=tt+  &
                 exp(-.5d0*r2)/fourpi *rw(k)  *(  &
                 gcore(1) )! &
                 ! + gcore(2)*r2   &
                 ! + gcore(3)*r2**2   &
                 ! + gcore(4)*r2**3 )
         end do
         write(6,*)'value for the radial grid used:  ',tt
         write(6,*)
         
         
         
         ! and a plot to compare with the full ae core charge
         open(13,file='nlcc.gplt')
         write(13,*)'rcore=',rcore
         write(13,*)'qcore=',qcore
         write(13,*)'const=',(znuc-zion)/sqrt2pi**3
         write(13,*)'ref="ae.core.dens.plt"'
         write(13,*)'g(x)=exp(-0.5*(x/rcore)**2)'
         write(13,*)'a(x)=qcore*const/rcore**3'
         write(13,*)'rhocore(x)=4*pi*x**2 *a(x)*g(x)'
         write(13,*)"set xrange [0.1*rcore:5.0*rcore]"
         write(13,*)"p rhocore(x),\" 
         write(13,*)"  ref       t 'ae core',\"
         write(13,*)"  ref u 1:3 t 'ae val'"
         write(13,*)"show function"
         write(13,*)"pr 'to fit the core charge, select an area and'"
         write(13,*)"pr 'fit rhocore(x) ref via rcore, qcore'"
         close(unit=13)
      end if
      
      
      
      
      ! dumpfile for testing with another program
      if (ldump .and. iproc == 0) then
         write(6,*)'writing out a dumpfile of',  &
              8*4+size(xp)+size(psi)+size(occup),'byte'           
         open(13,file='dumpfile.bin',form='unformatted')
         ! open(13,file='dumpfile')
         write(13)ng,noccmax,lmax,lpx,  &
              lcx,nspin,xp,psi,occup
         close(unit=13)
      end if
      
      ! here we used to overwrite old values of 'psp.par' with the current ones
      ! there was an info flag used to append some gaussian coeffs to psp.par
      if (iproc == 0 .and. namoeb /= 0) then
         if (info) then
            open(unit=23,file='gauss.par',form='formatted')
            write(23,*) 'additional information (last calculation):'
            write(23,*) 'gaussian exponents:'
            write(23,*) 'xp():',(xp(i),i=1,ng)
            write(23,*) 'orbital coefficients ',  &
                 '(one orbital per column)'
            do l=0,lmax
               ! no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
               do ispin=1,max(nspol,min(2*l+1,nspin))
                  write(23,*) 'l,ispin:',l,ispin
                  write(23,*) 'psi n=1,noccmax(l):'
                  do i=0,ng
                     write(23,'(10e20.10)')  &
                          (psi(i,nocc,l+1,ispin),  &
                          nocc=1,nomax(l))
                  end do
               end do
            end do
         endif
         close(unit=23)
      endif
      
      ! plot wavefunctions (up to 5*rcov)   c.hartwigsen: version for gnuplot
      
      ! note: for now, only do this with config0 from atom.00.ae
      ! should be straight forward to write out all plots in parallel
      if (plotwf.and.iproc==0) then
         ! generate various files for plotting
         ! write out the local potential in real space
         open(17,file='local.pot')
         write(17,'(a)')'# r, vloc-vion, vion=-z/r, erf term, gpot term'
         do k=1,nint
            r=rr(k)
            gt=exp(-.5d0*(r/rloc)**2)*  &
                 (gpot(1) + gpot(2)*(r/rloc)**2+    &
                 gpot(3)*(r/rloc)**4 +  &
                 gpot(4)*(r/rloc)**6 )
            et= -zion*derf(r/(sqrt(2.d0)*rloc))/r
            
            write(17,'(5e18.8)') r,et+gt+zion/r,zion/r, et, gt
            !.5d0*(r/rprb**2)**2
         end do
         close(unit=17)
         
         ! write out some plotable kernel for the separable part.
         ! for now let us try the diagonal part in real space
         ! and the contour at r'=rcov for each l-component.
         ! ignore the case kij /= 0 
         do l=1,lpx
            if (iproc>1) exit
            lq=l-1
            rnrm1=1.d0/sqrt(.5d0*gamma(lq+1.5d0)*r_l(l)**(2*lq+3))
            rnrm2=1.d0/sqrt(.5d0*gamma(lq+3.5d0)*r_l2(l)**(2*lq+7))
            rnrm3=1.d0/sqrt(.5d0*gamma(lq+5.5d0)*r_l2(l)**(2*lq+11))
            open(17,file=trim(il(l))//'.kernel.pot')
            write(17,'(3(9x,a))')'#   r          ',&         
                 'v_'//il(l)//'(r,r)   ',&
                 'v_'//il(l)//'(r,rcov)'
            do k=1,nint
               r=rr(k)
               ppr1=rnrm1*r**lq    *exp(-.5d0*(r/r_l(l))**2)
               ppr2=rnrm2*r**(lq+2)*exp(-.5d0*(r/r_l2(l))**2)
               ppr3=rnrm3*r**(lq+4)*exp(-.5d0*(r/r_l2(l))**2)
               ppc1=rnrm1*rcov**lq    *exp(-.5d0*(rcov/r_l(l))**2)
               ppc2=rnrm2*rcov**(lq+2)*exp(-.5d0*(rcov/r_l2(l))**2)
               ppc3=rnrm3*rcov**(lq+4)*exp(-.5d0*(rcov/r_l2(l))**2)
               write(17,'(3f20.5)')r, &
                    ppr1*hsep(1,l,1)*ppr1   + &
                    ppr1*hsep(2,l,1)*ppr2*2 + &  
                    ppr2*hsep(3,l,1)*ppr2   + &
                    ppr1*hsep(4,l,1)*ppr3*2 + &
                    ppr2*hsep(5,l,1)*ppr3*2 + &  
                    ppr3*hsep(6,l,1)*ppr3 , &
                    ! ---------------------------
               ppc1*hsep(1,l,1)*ppr1   + &
                    ppc1*hsep(2,l,1)*ppr2*2 + &  
                    ppc2*hsep(3,l,1)*ppr2   + &
                    ppc1*hsep(4,l,1)*ppr3*2 + &
                    ppc2*hsep(5,l,1)*ppr3*2 + &  
                    ppc3*hsep(6,l,1)*ppr3 
            end do
            close(unit=17)
         end do
         
         ! plotting of the orbitals
         call detnp(ngrid,rae,5.d0*rcov,np)
         open(32,file='orbitals.gplt',form='formatted',status='unknown')
         write(32,'(a)') 'set style data lines'
         do iorb=1,norb
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            if (so(iorb).lt.0) ispin=2
            if (methodps == 'r') then
               fname= 'ps.'//char(ichar('0')+noae(iorb))  &
                    //il(lo(iorb)+1)  &
                    //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))  &
                    //'by2.dat'
               tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)  &
                    //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))  &
                    //'/2'
            else
               ! either nonrel or spin pol
               fname= 'ps.'//char(ichar('0')+noae(iorb))  &
                    //il(lo(iorb)+1)
               tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)
               if (so(iorb)<0) then
                  fname=trim(fname)//'.down'
                  tname=trim(tname)//'.down'
               elseif (so(iorb)>0) then
                  fname=trim(fname)//'.up'
                  tname=trim(tname)//'.up'
               end if
               fname=trim(fname)//'.dat'
            endif
            open(unit=2,file=trim(fname),  &
                 form='formatted',status='unknown')
            ! find outer max of psi (approx), search from 10 bohr down
            ra=10.d0
            ttrmax=ra
            ttmax= dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
            do i=100,0, -1
               ra= 0.1d0 * i
               ttpsi=dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
               ! print*,ra,ttpsi
               if ( ttpsi .gt. ttmax  &
                    .and. ttpsi .gt. 1.0d-4 ) then
                  ttmax=ttpsi
                  ttrmax=ra
               endif
               if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) exit
            end do
            ! ae/pseudo wfs should have the same sign for large r when plotted
            call detnp(ngrid,rae,ttrmax,nsign)
            ! never use first gridpoint (only relevant for h and he)
            if (nsign == 1) nsign=nsign+1
            tt=psiold(nsign,nocc,l+1,ispin)
            sign1=tt/abs(tt)
            tt= wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(nsign))
            sign2=tt/abs(tt)
            do i=2,np
               ttold=psiold(i,nocc,l+1,ispin)*sign1*rae(i)
               ttold=max(min(3.d0,ttold),-3.d0)
               tt=wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))  &
                    *sign2*rae(i)
               tt=max(min(3.d0,tt),-3.d0)
               ttdiff=psiold(i,nocc,l+1,ispin)*sign1-  &
                    wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))*sign2
               ttdiff= ttdiff*rae(i)
               ttdiff=log(max(abs(ttdiff),1.d-8))/log(10.d0)
               write(2,'(7e20.10)') rae(i),ttold,tt,ttdiff
               ! plot of the wavefunction and the higher derivatives
               ! write(2,'(7e20.10)') rae(i),psiold(i,nocc,l+1,ispin), &
               ! wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)), &
               ! dwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)), &
               ! ddwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))
            end do
            close(unit=2)
            write(32,'(a)')'  plot "'//fname(1:index(fname,' ')-1)  &
                 //'"     title "'//tname(1:index(tname,' ')-1)  &
                 //'"'
            write(32,'(a)')'replot "'//fname(1:index(fname,' ')-1)  &
                 //'" using 1:3 title "pseudo"'
            write(32,'(a)') 'pause -10 "hit return to continue"'
            write(32,'(a)')'replot "'//fname(1:index(fname,' ')-1)  &
                 //'" using 1:4 title "diff"'
            write(32,'(a)') 'pause -10 "hit return to continue"'
         end do
         write(32,'(a)') 'set nokey'
         close(unit=32)

         write(6,'(/,1x,a)') 'gnuplot script written for comparison of valence '
         write(6,'(1x,a,/)')'                and pseudovalence: orbitals.gplt'
      endif
      
      
      
      ! -----------------------------------------------
      if (iiter>namoeb)exit
   end do ! main loop end
   
   
   ! test orthogonality of the projectors
   if (lpx.ge.1.and.ortprj)  &
        call pj2test(hsep,lpx-1,lpmx,lmx,nspin,nsmx,r_l,is)
   !
   write(6,*)
   write(6,*) 'total scf-cycles:       ',itertot
   write(6,*) 'pseudoatom calculations:',ntime
   call mpi_finalize(ierr)
   call libxc_functionals_end()
   
   call cpu_time(t)
   write(6,'(/,1x,a)')      'simple time profiling'
   write(6,'(1x,a,/)')      '_____________________'
   write(6,'(1x,a,f9.3,a)') 'gatom  overall runtime ',time(1),' seconds'
   write(6,'(1x,a,f9.3,a)') 'mpi communication time ',time(2),' seconds'
   write(6,'(1x,a,f9.3,a)') 'wavelet libraries time ',time(3),' seconds'
   write(6,'(1x,a)')        '________________________________________'
   write(6,'(1x,a,f9.3,a)') '              cpu time ',t-t0,' seconds'
   write(6,'(/,1x,a)')      '______________________________________________________'
   write(6,'(1x,a,/)')      '                                              finished'
   ! end
   
   
   ! optional info when call to the shell is available
   write(label,'(a,i3.3)') 'echo finished at `date` >> hostnames.',iproc
   call system(label)
end program pseudo

