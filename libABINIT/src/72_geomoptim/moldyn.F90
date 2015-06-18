!{\src2tex{textfont=tt}}
!!****f* ABINIT/moldyn
!! NAME
!! moldyn
!!
!! FUNCTION
!! perform dynamics on ions according to ionmov (see help file)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2009 ABINIT group (DCA, XG, GMR, JCC, JYR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number ofidum k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis and boundary at this k point.
!!  nspinor=number of spinorial components of the wavefunctions
!!  mxfh=last dimension of the xfhist array
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data paw
!!                                                 tabulated data read at start
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!
!! SIDE EFFECTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  nxfh=actual number of (x,f) history pairs, see xfhist array.
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point.
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  rprim(3,3)=dimensionless real space primitive translations
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  vel(3,natom)=cartesian velocities at the initialisation; updated on output
!!  wffnew,wffnow=struct info for wf disk files.
!!  xfhist(3,natom+4,2,mxfh)=(x,f) history array,
!!                                 also includes acell, rprim and stress
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=work space for old xred
!!
!! NOTES
!! For ionmov=6 :
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), a velocity vector (in cartesian
!! coordinates), and unit cell parameters (acell and rprim - without
!! velocities in the present implementation), the
!! Verlet dynamics is performed, using the gradient
!! of the energy (atomic forces and
!! stress : fred or fcart and stress) as calculated by the routine scfcv.
!! Some atoms can be kept fixed, while the propagation of unit cell
!! parameters is only performed if optcell/=0.
!! No more than dtset%ntime steps are performed.
!! The time step is governed by dtion (contained in dtset)
!! Returned quantities are xred, and eventually acell and rprim (new ones!).
!!
!! For ionmov=7 :
!! Block every atom for which the scalar product of velocity and
!! forces is negative, in order to reach the minimum.
!! The convergence requirement on
!! the atomic forces, dtset%tolmxf,  allows an early exit.
!!
!! For  ionmov=8
!! See ionmov=6, but with a nose-hoover thermostat
!! Velocity verlet algorithm : Swope et al JCP 76 (1982) 637
!!
!! For  ionmov=9
!! Uses a Langevin dynamics algorithm :
!! see J. Chelikowsky, J. Phys. D : Appl Phys. 33(2000)R33
!!
!! For  ionmov=12
!!    Application of Gauss' principle of least constraint according
!!    to Fei Zhang's algorithm (J. Chem. Phys. 106, 1997, p.6102);
!!    see also Minary et al. (J. Chem. Phys. 118, 2003, p.2510)
!!
!! NOTE : there are many similarities between this routine
!! and brdmin.f, so that they might have to be maintained
!! together. Some common subroutines might be extracted from them.
!!
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      chkexi,fconv,initylmg,leave_new,metric,mkrdim,prtxvf,scfcv,status
!!      write_header_moldynnetcdf,write_moldynvaluenetcdf,wrtout,xfpack
!!      xredxcart
!!
!! SOURCE

subroutine moldyn(acell,amass,me,&
& mxfh,nxfh,natom,rprim,etotal,iexit,&
& optcell, ionmov, ntime, dtion, noseinert, mditemp, mdftemp, &
& friction, mdwall, nnos, qmass, bmass, vmass, iatfix, strtarget, &
& strprecon, strfact, tolmxf, &
& nsym, symrel, &
& vel,xfhist,fred,xred)

  use interfaces_14_hidewrite

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
! nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!scalars
 integer,intent(in) :: mxfh,natom, nsym, optcell, ionmov, ntime, nnos, me
 integer,intent(inout) :: nxfh
 integer,intent(out) :: iexit
 real(dp), intent(in) :: dtion, noseinert, mditemp, mdftemp, friction, mdwall
 real(dp), intent(in) :: strprecon, strfact, tolmxf, bmass, vmass
 real(dp), intent(inout) :: etotal
!arrays
 integer, intent(in) :: iatfix(3, natom)
 real(dp),intent(in) :: amass(natom), qmass(nnos), strtarget(6)
 real(dp),intent(inout) :: acell(3)
 real(dp),intent(in) :: symrel(3,3,nsym)
 real(dp),intent(inout) :: rprim(3,3),vel(3,natom)
 real(dp),intent(inout) :: xfhist(3,natom+4,2,mxfh),xred(3,natom),fred(3,natom)
!Local variables ------------------------------
!scalars
 integer,parameter :: level=5
 integer :: iatom,idim,idir,idir1,idir2,ierr
 integer :: ii,ipos,itime
 integer :: ndim,nstopped,option
 integer :: prtvel,prtvol
 real(dp) :: diag,ekin, ekin_corr,etotal_prev,favg,gnose, ktemp
 real(dp) :: massvol,snose,ucvol,ucvol0, ucvol_next,v2nose,xi_nose
 character(len=500) :: message
 type(mttk_type) :: mttk_vars
!arrays
 real(dp) :: acell0(3),acell_next(3),angle(3),gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3)
 real(dp) :: rprim_next(3,3),rprimd(3,3),rprimd0(3,3),rprimd_next(3,3)
 real(dp) :: strten(6)
 real(dp),allocatable :: fcart(:,:)
 real(dp),allocatable :: fcart_mold(:,:),fred_corrected(:,:)
 real(dp),allocatable :: hessin(:,:)
 real(dp),allocatable :: vel_nexthalf(:,:)
 real(dp),allocatable :: vel_prevhalf(:,:),vin(:),vin_next(:)
 real(dp),allocatable :: vin_prev(:),vout(:),xcart(:,:)
 real(dp),allocatable :: xcart_next(:,:),xred_prev(:,:),xred_next(:,:)

!************************************************************************
!Beginning of executable session
!***************************************************************************

!XG020711 : this line, tentatively added when setting up v3.4, is incorrect
!nxfh is initialized in the calling routine
!nxfh=0
 !!$ call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=0
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' moldyn : enter '
  call wrtout(std_out,message,'COLL')
 end if

 ipos=0

!dtion=time step for molecular dynamics in atomic time units
!(1 atomic time unit=2.418884e-17 seconds)

 ndim=3*natom
 if(optcell==1 .or. optcell==4 .or. optcell==5 .or. optcell==6)ndim=ndim+1
 if(optcell==2 .or. optcell==3)ndim=ndim+6
 if(optcell==7 .or. optcell==8 .or. optcell==9)ndim=ndim+3

 allocate(fcart(3,natom))
 allocate(fcart_mold(3,natom))
 allocate(fred_corrected(3,natom))
 allocate(vel_nexthalf(3,natom),vel_prevhalf(3,natom))
 allocate(vin(ndim),vin_next(ndim))
 allocate(vin_prev(ndim))
 allocate(vout(ndim))
 allocate(xcart(3,natom),xcart_next(3,natom))
 allocate(xred_next(3,natom),xred_prev(3,natom))

!Compute dimensional primitive translations rprimd, then metric tensor gmet
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Save initial values
 acell0(:)=acell(:)
 rprimd0(:,:)=rprimd(:,:)
 ucvol0=ucvol
 strten=0_dp

 if (ionmov==6 .and. maxval(vel) == zero) then
    call md_nose_init(amass, natom, mditemp, vel)
 end if

!Initialize input vectors : first vin, then vout
 option=1
 call xfpack(acell,acell0,fred,natom,ndim,&
& nsym,optcell,option,rprim,rprimd0,&
& strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)
 option=3
 call xfpack(acell,acell0,fred,natom,ndim,&
& nsym,optcell,option,rprim,rprimd0,&
& strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)

!Here, set up the matrix of transformation between forces and
!acceleration. Masses must be included here.
!Beside this feature, one could define
!a preconditioner, in which case it should
!be the inverse hessian, like in Broyden. This explains the
!name chosen for this transformation matrix. This would allow
!to find easily the optimal geometry with ionmov=7.
!The default, now implemented, corresponds to the identity matrix
!in cartesian coordinates, which makes use of metric tensor gmet
!in reduced coordinates.
 allocate(hessin(ndim,ndim))
 hessin(:,:)=0.0_dp
 do iatom=1,natom
  do idir1=1,3
   do idir2=1,3
!   Warning : implemented in reduced coordinates
    if ( iatfix(idir1,iatom) ==0 .and. iatfix(idir2,iatom) ==0 )then
     hessin(idir1+3*(iatom-1),idir2+3*(iatom-1))=&
&     gmet(idir1,idir2)/amass(iatom)
    end if
   end do
  end do
 end do
 if(optcell/=0)then
! These values might lead to too large changes in some cases ...
! No "mass" is included here
  diag=strprecon*30.0_dp/ucvol
  if(optcell==1)diag=diag/3.0_dp
  do idim=3*natom+1,ndim
   hessin(idim,idim)=diag
  end do
 end if

!-----------------------------------------------------------------------
!
!Iterative procedure (main loop)
!
 do itime=0,ntime

!!$  if(itime/=0)then
!!$   call status(itime,dtfil%filstat,iexit,level,'loop itime    ')
!!$  end if

! Check whether exiting was required by the user.
! If found then beat a hasty exit from loop on itime
!!$  openexit=1 ; if(dtset%chkexit==0) openexit=0
!!$  call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
!!$  if (iexit/=0) exit

  write(message, '(a,a,i4,a)' ) ch10,' MOLDYN STEP NUMBER ',itime,&
&  '  ------------------------------------------------------'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')

  if (ionmov==8.or.ionmov==9.or.ionmov==13) then
!  The temperature is linear between initial and final values
!  It is here converted from Kelvin to Hartree (kb_HaK)
   ktemp=(mditemp+((mdftemp-mditemp)/dble(ntime))*itime)*kb_HaK
  end if

! If not initialisation time step
  if(itime>0)then

!  Shift the data
   etotal_prev=etotal
   acell(:)=acell_next(:)
   rprim(:,:)=rprim_next(:,:)
   rprimd(:,:)=rprimd_next(:,:)
   ucvol=ucvol_next
   vel_prevhalf(:,:)=vel_nexthalf(:,:)
   vin_prev(:)=vin(:)
   vin(:)=vin_next(:)
   xcart(:,:)=xcart_next(:,:)
   xred_prev(:,:)=xred(:,:)
   xred(:,:)=xred_next(:,:)
  else

!  Get xred, and eventually acell, rprim and rprimd, from current vin
   option=2
   call xfpack(acell,acell0,fred,natom,ndim,&
&   nsym,optcell,option,rprim,rprimd0,&
&   strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)

!  End condition on itime
  end if

  !!$ call status(itime,dtfil%filstat,iexit,level,'call scfcv    ')

  if(optcell/=0)then

   call mkrdim(acell,rprim,rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  Write, but only to log file
   write(message, '(a,a,4(a,3es18.10,a),a,es18.10,a,a)' )&
&   ' Unit cell characteristics (before scfcv) :',ch10,&
&   '  acell=',acell(1:3),ch10,&
&   '  rprim=',rprim(1:3,1),ch10,&
&   '        ',rprim(1:3,2),ch10,&
&   '        ',rprim(1:3,3),ch10,&
&   '  ucvol=',ucvol,' Bohr^3',ch10
   call wrtout(std_out,message,'COLL')

  end if

! If metric has changed since the initialization, update the Ylm's
!!$  if (optcell/=0.and.psps%useylm==1.and.itime>0)then
!!$   !!$ call status(0,dtfil%filstat,iexit,level,'call initylmg ')
!!$   option=0;if (dtset%iscf>0) option=1
!!$   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
!!$&   npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
!!$  end if

  if((ionmov/=12.and.ionmov/=13.and.ionmov/=14).or.&
&  (ionmov==12.and.itime==0).or.(ionmov==13.and.itime==0)&
&  ) then
!  Compute LDA forces (big loop)
     if (itime /= 0 .or. maxval(fred) == zero) then
        call scfloop_main(acell, etotal, fcart, fred, itime, me, natom, rprimd, xred)
     else
        call xredxcart(natom,1,rprimd,fcart,fred)
     end if
     call xredxcart(natom,1,rprimd,xcart,xred)
!!$   iapp=-1
!!$   if(itime>0)iapp=itime
!!$   call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
!!$&   dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
!!$&   irrzon,kg,mpi_enreg,&
!!$&   nattyp,nfftf,npwarr,nspinor,occ,&
!!$&   pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&   pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&   scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
  end if
  !!$ call status(itime,dtfil%filstat,iexit,level,'call prtxvf   ')

! Output of acell and/or rprim ( and angles ! - should become a routine later)
  if(optcell/=0)then
   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0
   write(message, '(a,a,4(a,3es18.10,a))' )&
&   ' Unit cell characteristics :',ch10,&
&   '  acell=',acell(1:3),ch10,&
&   '  rprim=',rprim(1:3,1),ch10,&
&   '        ',rprim(1:3,2),ch10,&
&   '        ',rprim(1:3,3)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es18.10,a,a,a,3es18.10,a,a,a,3f13.8,a)' )&
&   '  ucvol=',ucvol,' Bohr^3',ch10,&
&   '  lengths=',sqrt(rmet(1,1)),sqrt(rmet(2,2)),sqrt(rmet(3,3)),' Bohr',&
&   ch10,'  angles (23,13,12)=',angle(1:3),' degrees'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
  end if

! Get rid off mean force on whole unit cell
  do idir=1,3
   favg=sum(fred(idir,:))/dble(natom)
   fred_corrected(idir,:)=fred(idir,:)-favg
!!$   if(dtset%jellslab/=0.and.idir==3) fred_corrected(idir,:)=fred(idir,:)
  end do

! Update xfhist
  nxfh=nxfh+1
  xfhist(:,1:natom,1,nxfh)=xred(:,:)
  xfhist(:,natom+1,1,nxfh)=acell(:)
  xfhist(:,natom+2:natom+4,1,nxfh)=rprim(1:3,1:3)
  xfhist(:,1:natom,2,nxfh)=fred_corrected(:,:)
  xfhist(:,natom+2,2,nxfh)=strten(1:3)
  xfhist(:,natom+3,2,nxfh)=strten(4:6)

! Store computed gradient in vout
  option=3
  call xfpack(acell,acell0,fred_corrected,&
&  natom,ndim,nsym,optcell,option,rprim,rprimd0,&
&  strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)

! ####### Test case ionmov=9  Langevin dynamics  ##########
  if (ionmov==9) then
     call md_langevin(amass, dtion, fcart, fcart_mold, friction, itime, ktemp, &
          & mditemp, mdwall, natom, rprimd, vel, xcart, xcart_next, xred_next)
     !  Impose no change of acell, ucvol, rprim, and rprimd
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     rprim_next(:,:)=rprim(:,:)
     rprimd_next(:,:)=rprimd(:,:)

!  #####   case ionmov==8  Nose dynamics ###########
  else if (ionmov==8) then
     call md_nose(amass, dtion, fcart, fcart_mold, gnose, itime, ktemp, mditemp, &
          & natom, noseinert, rprimd, snose, v2nose, vel, xcart, xcart_next, &
          & xi_nose, xred_next)
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     rprim_next(:,:)=rprim(:,:)
     rprimd_next(:,:)=rprimd(:,:)

!  #####   case ionmov==12 Isokinetic Ensemble ###########
  else if (ionmov==12) then
     call md_isokinetic(acell, amass, dtion, etotal, fcart, itime, natom, &
          & mditemp, me, rprimd, vel, vel_nexthalf, xcart, xcart_next, xred_next)
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     rprim_next(:,:)=rprim(:,:)
     rprimd_next(:,:)=rprimd(:,:)

!  #####   case ionmov==13 Reversible integrator of Martyna at al.###########
!  There are three sub cases according to the value of optcell
!  optcell=0 means isothermal, optcell==1:homogeneous cell fluctuations
!  optcell=2: full cell fluctuation in addition to temperature control.
  else if (ionmov==13) then
     call md_isothermal(acell, acell_next, amass, bmass, dtion, etotal, massvol, &
          & fcart, iatfix, itime, ktemp, mditemp, me, mttk_vars, natom, nnos, &
          & optcell, qmass, rprim, rprimd, rprim_next, rprimd_next, strten, &
          & strtarget, ucvol, ucvol_next, vel, vel_nexthalf, vmass, xcart, &
          & xcart_next, xred_next)

  else if(ionmov==14) then
     write(*,*) "TODO"
  
  !  ####### Begin case ionmov/=8 and ionmov/=9 and ionmov/=12 and ionmov/=13 and ionmov /= 14 ##########
  else
     call md_velocity_verlet(acell, acell_next, amass, dtion, fred_corrected, &
     & hessin, iatfix, itime, natom, optcell, rprim, &
     & rprim_next, rprimd, rprimd_next, &
     & ucvol, ucvol_next, vel, vel_nexthalf, vel_prevhalf, &
     & xcart, xcart_next, xred_next, xred_prev)
  end if

! Compute the ionic kinetic energy (no cell shape kinetic energy yet)
  ekin=0.0_dp
  do iatom=1,natom
   do idir=1,3
!   Warning : the fixing of atomis is implemented in reduced
!   coordinates, so that this expression is wrong
    if (iatfix(idir,iatom) == 0) then
     ekin=ekin+0.5_dp*amass(iatom)*vel(idir,iatom)**2
    end if
   end do
  end do

! Output coordinates, forces and velocities
  prtvel=1
  call prtxvf(fcart,iatfix,ab_out,natom,prtvel,vel,xcart)
  call prtxvf(fcart,iatfix, 06 ,natom,prtvel,vel,xcart)

! Here, stop the atoms for which the scalar product of velocity
! and force is negative, and recompute the kinetic energy.
  if(ionmov==7)then
     call md_quenched_stop_atoms(amass, dtion, ekin_corr, fcart, iatfix, itime, &
          & natom, nstopped, rprimd, vel, vel_prevhalf, vel_nexthalf, &
          & xcart, xcart_next, xred_next)
     !   Store xred_next, and eventual acell_next and rprim_next in vin
!!$     option=1
!!$     call xfpack(acell_next,acell0,fred_corrected,&
!!$          &    natom,ndim,nsym,optcell,option,rprim_next,rprimd0,&
!!$          &    strtarget,strten,symrel,ucvol_next,ucvol0,vin_next,vout,xred_next)
  end if

  if (ionmov==8) then
     call md_nose_finalise(etotal, gnose, itime, ktemp, noseinert, snose, &
          & v2nose, xi_nose)
  end if

! Output total energy in a format that can be captured easily
  write(message, '(a,i6,a,es22.14,a)' )&
&  ' At the end of Moldyn step',itime,', POT.En.=',etotal,' Ha.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  write(message, '(a,es22.14,a)' )&
&  '                                  KIN.En.=',&
&  ekin,' Ha.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  write(message, '(a,es22.14,a)' )&
&  '                              KIN+POT.En.=',&
&  etotal+ekin,' Ha.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  if(ionmov==7 .and. nstopped/=0)then
   write(message, '(a,es22.14,a)' )&
&   '                    corrected KIN+POT.En.=',&
&   etotal+ekin_corr,' Ha.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
  end if
  if(ionmov==9)then
   write(message, '(a,es22.14,2x,es22.14)' )&
&   '           TEMP         TKA =',&
&   ktemp/kb_HaK,ekin/(1.5_dp*natom*ktemp)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
  end if

  call scfloop_output(acell, etotal, ekin, fred, itime, me, natom, rprimd, vel, xred)

! Writing outpout netcdf
! the output is being done every nctime time
!!$  if ( iproc == 0) then
!!$   nb1 = size(strten)
!!$   nbdir = size(xcart,1)
!!$   if (dtset%nctime > 0)then
!!$    if (itime == 0)then
!!$     call  write_header_moldynnetcdf(dtfil, dtset, natom, nbdir, nb1 )
!!$    end if
!!$    if ( mod (itime, dtset%nctime ) == 0)then
!!$     ipos = ipos +1
!!$     call write_moldynvaluenetcdf(amass, ipos, dtfil, dtset, etotal, ekin, &
!!$&     natom, nbdir, nb1, xcart, vel, strten, rprimd, ucvol )
!!$     open(dtfil%unpos,file='POSABIN',status='replace',form='formatted')
!!$     do iatom=1,natom
!!$      if(iatom==1) then
!!$       write(dtfil%unpos,'(a7,3d18.5)') 'xred  ',(xred_next(idim,iatom),idim=1,3)
!!$      else
!!$       write(dtfil%unpos,'(3d18.5)') (xred_next(idim,iatom),idim=1,3)
!!$      end if
!!$     end do
!!$     do iatom=1,natom
!!$      if(iatom==1) then
!!$       write(dtfil%unpos,'(a7,3d18.5)') 'vel  ',(vel(idim,iatom),idim=1,3)
!!$      else
!!$       write(dtfil%unpos,'(3d18.5)') (vel(idim,iatom),idim=1,3)
!!$      end if
!!$     end do
!!$     close(dtfil%unpos)
!!$    end if
!!$   end if
!!$  end if

! Check whether forces and stresses are below tolerance; if so, exit
! from the itime loop
  if(ionmov==7)then
   iexit=0
   if(itime==ntime)iexit=1
   call fconv(fcart,iatfix,iexit,itime,natom,ntime,&
&   optcell,strfact,strtarget,strten,tolmxf)
   if (iexit/=0) exit
  end if

! Back to another iteration in the itime loop.
! Note that there are "exit" instructions inside the loop.
 end do

!-----------------------------------------------------------------------

 if(ionmov==8)then
  write(message, '(a,i6)')'Nb time steps',nxfh
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  do ierr=1,nxfh
   write(message, '(a,i6)')'step ',ierr
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   do iatom = 1, natom
    write(message, '(a,i6,3es22.14)')'atom ',iatom,xfhist(1:3,iatom,1,ierr)
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
   end do
  end do
 end if

 call xredxcart(natom,-1,rprimd,fcart,fred)

 deallocate(fcart)
 deallocate(fcart_mold)
 deallocate(fred_corrected,hessin)
 deallocate(vel_nexthalf,vel_prevhalf)
 deallocate(vin,vin_next,vin_prev,vout)
 deallocate(xcart,xcart_next,xred_prev,xred_next)
 if(ionmov==13) then
  deallocate(mttk_vars%glogs,mttk_vars%vlogs,mttk_vars%xlogs)
 end if
!Structured debugging : if prtvol=-level, stop here.
!!$ if(prtvol==-level)then
!!$  write(message,'(a1,a,a1,a,i1,a)') ch10,' moldyn : exit ',&
!!$&  ch10,'  prtvol=-',level,', debugging mode => stop '
!!$  call wrtout(std_out,message,'COLL')
!!$  call leave_new('COLL')
!!$ end if

 !!$ call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine moldyn
!!***
