!> @file
!!  Routines to do geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Geometry optimization, parametrisation routine.
subroutine geopt_init()
  use minpar
  implicit none

  parmin%approach  = 'unknown'
  parmin%iter      = 0
  parmin%iflag     = 0
  parmin%verbosity = 1
  parmin%MSAVE=7
  parmin%MP=16
  parmin%LP=16
  parmin%MAXFEV=10
  parmin%GTOL=9.d-1
  parmin%XTOL=1.d-15
  parmin%FTOL=1.d-6
  parmin%STPMIN=1.d-20
  parmin%STPMAX=20.d0
  parmin%DIAGCO=.FALSE.
  parmin%IWRITE=.FALSE.

END SUBROUTINE geopt_init


!> Geometry optimization, parametrisation routine.
subroutine geopt_set_verbosity(verb)
  use minpar
  implicit none
  !Arguments
  integer, intent(in) :: verb
  parmin%verbosity = verb
END SUBROUTINE geopt_set_verbosity


!> Geometry optimization
subroutine geopt(runObj,outs,nproc,iproc,ncount_bigdft)
  use module_base
  use module_interfaces, except_this_one => geopt
  use module_types
  use yaml_output
  use dictionaries
  use minpar
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  integer, intent(in) :: nproc,iproc
  integer, intent(out) :: ncount_bigdft
  !local variables
  logical :: fail
  integer :: ibfgs,ierr
  character(len=6) :: outfile, fmt
  character(len=5) :: fn4
  character(len=40) :: comment
  character(len=60) :: filename
  !-------------------------------------------

  !Geometry Initialization
  call geopt_init()

  filename=trim(runObj%inputs%dir_output)//'geopt.mon'
  !open(unit=16,file=filename,status='unknown',position='append')
  !without istat this opening should crash 
  !do nothing if the unit is already opened
  if (iproc == 0) call yaml_set_stream(unit=16,filename=trim(filename),tabbing=0,record_length=100,setdefault=.false.,istat=ierr)
  if (iproc ==0 ) call yaml_comment('Geopt file opened, name: '//trim(filename)//', timestamp: '//trim(yaml_date_and_time_toa()),&
       hfill='-',unit=16)
  !write(16,*) '----------------------------------------------------------------------------'

  if (iproc == 0 .and. parmin%verbosity > 0 .and. ierr /= 0)  &!write(16,'(a)')  & 
     call yaml_comment('Geometry optimization log file, grep for GEOPT for consistent output',unit=16)
  if (iproc == 0 .and. parmin%verbosity > 0) write(16,'(a)')  & 
      '# COUNT  IT  GEOPT_METHOD  ENERGY                 DIFF       FMAX       FNRM      FRAC*FLUC FLUC      ADD. INFO'

  !if (iproc ==0 .and. parmin%verbosity > 0) write(* ,'(a)') & 
  !    '# COUNT  IT  GEOPT_METHOD  ENERGY                 DIFF       FMAX       FNRM      FRAC*FLUC FLUC      ADD. INFO'

  !assign the geometry optimisation method
  parmin%approach=runObj%inputs%geopt_approach

  outs%strten=0 !not used for the moment

  outs%energy=0.d0
  ncount_bigdft=0
  if (iproc == 0) then
     outfile = 'posout'
     if (trim(parmin%approach)=='AB6MD') outfile = 'posmd '
     fmt = "(i4.4)"
     if (trim(parmin%approach)=='AB6MD') fmt = '(i5.5)'
     write(fn4,fmt) ncount_bigdft
     write(comment,'(a)')'INITIAL CONFIGURATION '
     call write_atomic_file(trim(runObj%inputs%dir_output)//trim(outfile)//'_'//trim(fn4),&
          & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),forces=outs%fxyz)
     call yaml_new_document()
     call yaml_comment('Geometry minimization using ' // trim(parmin%approach),hfill='-')
     call yaml_map('Begin of minimization using ',parmin%approach)
     !write(*,'(a,1x,a)') ' Begin of minimization using ',parmin%approach
  end if

  select case(trim(parmin%approach))

  case('LBFGS')
     ibfgs=0
     do
        ibfgs=ibfgs+1
        if (iproc ==0) call yaml_map('ENTERING LBFGS,ibfgs',ibfgs)
        call lbfgsdriver(runObj,outs,nproc,iproc,ncount_bigdft,fail)
        if (.not. fail .or. ibfgs .ge. 5) exit
     end do

     if (fail) then
        if (iproc ==0) call yaml_map('ENTERING CG after LBFGS failure,ibfgs',ibfgs)
        call conjgrad(runObj,outs,nproc,iproc,ncount_bigdft)
     end if

!  if(trim(parmin%approach)=='LBFGS') then
!
!     if (iproc ==0) write(*,*) '# ENTERING BFGS'
!
!     call bfgs(nproc,iproc,pos,fxyz,outs%energy,at,rst,in,ncount_bigdft,fail)
!
!     if (fail) then 
!        if (iproc ==0) write(*,*) '# ENTERING CG after BFGS failure'
!        call conjgrad(nproc,iproc,pos,at,outs%energy,fxyz,rst,in,ncount_bigdft)
!     end if
!
  case('BFGS','PBFGS')
     call bfgsdriver(runObj,outs,nproc,iproc,ncount_bigdft)
  case('SDCG')
     if (iproc ==0) call yaml_map('ENTERING CG',ncount_bigdft)
     call conjgrad(runObj,outs,nproc,iproc,ncount_bigdft)
  case('VSSD')
     if (iproc ==0) call yaml_map('ENTERING VSSD',ncount_bigdft)
     call vstepsd(runObj,outs,nproc,iproc,ncount_bigdft)
  case('FIRE')
     if (iproc ==0) call yaml_map('ENTERING FIRE',ncount_bigdft)
     call fire(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case('SBFGS')                                                  
   if (iproc ==0) call yaml_map('ENTERING SBFGS',ncount_bigdft)
   call sbfgs(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case('DIIS')   
     if (iproc ==0) call yaml_map('ENTERING DIIS',ncount_bigdft)
     call rundiis(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case('AB6MD')
     if (iproc ==0) call yaml_map('ENTERING Molecular Dynamics (ABINIT implementation)',ncount_bigdft)
     call ab6md(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case default
     call f_err_throw('Geometry optimization method undefined ('//trim(parmin%approach)//')',&
          err_name='BIGDFT_RUNTIME_ERROR')
     return
  end select

  if (iproc==0) call yaml_close_stream(unit=16)

  if (iproc==0) call yaml_map('End of minimization using ',parmin%approach)

  if (iproc==0) call finaliseCompress()

END SUBROUTINE geopt


!> Molecular Dynamics
subroutine ab6md(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  use module_base
  use module_types
  use scfloop_API
  use ab6_moldyn
!  use module_interfaces, only: memocc
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  logical, intent(out) :: fail
  !Local variables
  character(len=*), parameter :: subname='ab6md'
  ! 1 atomic mass unit, in electronic mass
  real(gp), parameter :: amu2emass=1.660538782d-27/9.10938215d-31
  integer, allocatable, dimension(:,:) :: iatfix
  real(gp), allocatable, dimension(:,:,:,:) :: xfhist
  real(gp), allocatable, dimension(:,:) :: vel, xred, fred
  real(gp), allocatable, dimension(:) :: amass
  real(gp), dimension(3) :: acell
  real(gp), dimension(3,3) :: rprim
  real(gp), dimension(3,3,1) :: symrel
  integer :: nxfh, iat, idim, iexit, i_stat,i_all

  ! We save pointers on data used to call bigdft() routine.
  if (ncount_bigdft == 0) then
     call scfloop_init(nproc, runObj)
  end if

  ! Prepare the objects used by ABINIT.
  allocate(amass(runObj%atoms%astruct%nat),stat=i_stat)
  call memocc(i_stat,amass,'amass',subname)
  allocate(xfhist(3, runObj%atoms%astruct%nat + 4, 2, runObj%inputs%ncount_cluster_x+1))
  call memocc(i_stat,xfhist,'xfhist',subname)
  allocate(vel(3, runObj%atoms%astruct%nat),stat=i_stat)
  call memocc(i_stat,vel,'vel',subname)
  allocate(xred(3, runObj%atoms%astruct%nat),stat=i_stat)
  call memocc(i_stat,xred,'xred',subname)
  allocate(iatfix(3, runObj%atoms%astruct%nat),stat=i_stat)
  call memocc(i_stat,iatfix,'iatfix',subname)
  allocate(fred(3, runObj%atoms%astruct%nat),stat=i_stat)
  call memocc(i_stat,fred,'fred',subname)

  nxfh = 0
  !acell = (/ runObj%atoms%astruct%cell_dim(1), runObj%atoms%astruct%cell_dim(2), runObj%atoms%astruct%cell_dim(3) /)
  acell(1)=runObj%atoms%astruct%cell_dim(1)
  acell(2)=runObj%atoms%astruct%cell_dim(2)
  acell(3)=runObj%atoms%astruct%cell_dim(3)
  rprim(:,:) = 0.0_gp
  rprim(1,1) = real(1, gp)
  rprim(2,2) = real(1, gp)
  rprim(3,3) = real(1, gp)
  symrel(:,:,1) = 0.0_gp
  symrel(1,1,1) = real(1, gp)
  symrel(2,2,1) = real(1, gp)
  symrel(3,3,1) = real(1, gp)
  do iat = 1, runObj%atoms%astruct%nat
     amass(iat) = amu2emass * runObj%atoms%amu(runObj%atoms%astruct%iatype(iat))
     do idim = 1, 3
        xred(idim, iat) = runObj%atoms%astruct%rxyz(idim, iat) / acell(idim)
        fred(idim, iat) = - outs%fxyz(idim, iat) / acell(idim)
     end do
     select case(runObj%atoms%astruct%ifrztyp(iat))
     case(0)
        iatfix(:, iat) = 0
     case(1)
        iatfix(:, iat) = 1
     case(2)
        iatfix(:, iat) = 0
        iatfix(2, iat) = 1
     case(3)
        iatfix(:, iat) = 1
        iatfix(2, iat) = 0
     end select
  end do

  !read the velocities from input file, if present
  call read_velocities(iproc,'velocities.xyz',runObj%atoms,vel)
  !vel(:,:) = zero

  ! Call the ABINIT routine.
  ! currently, we force optcell == 0
  call moldyn(acell, amass, iproc, runObj%inputs%ncount_cluster_x+1, nxfh, runObj%atoms%astruct%nat, &
       & rprim, outs%energy, iexit, &
       & 0, runObj%inputs%ionmov, runObj%inputs%ncount_cluster_x, runObj%inputs%dtion, runObj%inputs%noseinert, &
       & runObj%inputs%mditemp, runObj%inputs%mdftemp, runObj%inputs%friction, runObj%inputs%mdwall, runObj%inputs%nnos, &
       & runObj%inputs%qmass, runObj%inputs%bmass, runObj%inputs%vmass, iatfix, runObj%inputs%strtarget, &
       & runObj%inputs%strprecon, runObj%inputs%strfact, runObj%inputs%forcemax, &
       & 1, symrel, vel, xfhist, fred, xred)

  do iat = 1, runObj%atoms%astruct%nat, 1
     outs%fxyz(1, iat) = fred(1, iat) * acell(1)
     outs%fxyz(2, iat) = fred(2, iat) * acell(2)
     outs%fxyz(3, iat) = fred(3, iat) * acell(3)
  end do

  !De-Allocations
  i_all=-product(shape(fred))*kind(fred)
  deallocate(fred,stat=i_stat)
  call memocc(i_stat,i_all,'fred',subname)
  i_all=-product(shape(iatfix))*kind(iatfix)
  deallocate(iatfix,stat=i_stat)
  call memocc(i_stat,i_all,'iatfix',subname)
  i_all=-product(shape(xred))*kind(xred)
  deallocate(xred,stat=i_stat)
  call memocc(i_stat,i_all,'xred',subname)
  i_all=-product(shape(vel))*kind(vel)
  deallocate(vel,stat=i_stat)
  call memocc(i_stat,i_all,'vel',subname)
  i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
  call memocc(i_stat,i_all,'amass',subname)
  i_all=-product(shape(xfhist))*kind(xfhist)
  deallocate(xfhist,stat=i_stat)
  call memocc(i_stat,i_all,'xfhist',subname)

  fail = (iexit == 0)
END SUBROUTINE ab6md


!> MODIFIED version for refined time limit on restart of global.f90.
!! Only difference: Calls routine CPUtime(tt)
subroutine timeleft(tt)
  use module_base
  implicit none
  real(gp), intent(out) :: tt
  !local variables
  integer :: ierr
  real(kind=4) :: tcpu
  real(gp) :: timelimit

  open(unit=55,file='CPUlimit',status='unknown')
  read(55,*,iostat=ierr) timelimit ! in hours
  if(ierr/=0)timelimit=1d6
  close(55)
  call cpu_time(tcpu)
  tt=timelimit-real(tcpu,gp)/3600._gp ! in hours
END SUBROUTINE timeleft


subroutine convcheck(fmax,fluctfrac_fluct,forcemax,check) !n(c) fnrm (arg:1)
  use module_base
  implicit none
  real(gp), intent(in):: fmax, fluctfrac_fluct,forcemax !n(c) fnrm
!  logical, intent(out)::check
  integer, intent(inout)::check

!  check=.false.
!  if ( fmax < max(forcemax,fluctfrac_fluct)) check=.true.
  if ( fmax < max(forcemax,fluctfrac_fluct)) then 
    check=check+1
  else
    check=0
  endif

END SUBROUTINE convcheck


subroutine fnrmandforcemax(ff,fnrm,fmax,nat)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), intent(in):: ff(3,nat)
  real(gp), intent(out):: fnrm, fmax
  real(gp):: t1,t2,t3
  integer:: iat

  fmax=0._gp
  do iat=1,nat
     t1=ff(1,iat)**2
     t2=ff(2,iat)**2
     t3=ff(3,iat)**2
     fmax=max(fmax,sqrt(t1+t2+t3))
  enddo

  !this is the norm of the forces of non-blocked atoms
  !one has to discuss whether also the center mass shift should be added
  fnrm=dot(3*nat,ff(1,1),1,ff(1,1),1)
END SUBROUTINE fnrmandforcemax


subroutine fnrmandforcemax_old(ff,fnrm,fmax,at)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), intent(in):: ff(3,at%astruct%nat)
  real(gp), intent(out):: fnrm, fmax
  real(gp):: t1,t2,t3
  integer:: iat

  fmax=0._gp
  do iat=1,at%astruct%nat
     call frozen_alpha(at%astruct%ifrztyp(iat),1,ff(1,iat)**2,t1)
     call frozen_alpha(at%astruct%ifrztyp(iat),2,ff(2,iat)**2,t2)
     call frozen_alpha(at%astruct%ifrztyp(iat),3,ff(3,iat)**2,t3)
     fmax=max(fmax,sqrt(t1+t2+t3))
  enddo

  !This is the norm of the forces of non-blocked atoms
  call atomic_dot(at,ff,ff,fnrm)
END SUBROUTINE fnrmandforcemax_old


subroutine updatefluctsum(fnoise,fluct) !n(c) nat (arg:1)
  use module_base
  use module_types
  implicit none
  !n(c) integer, intent(in) :: nat
  real(gp),intent(in):: fnoise
  real(gp),intent(inout):: fluct

   if (fluct == 0.d0) then
     fluct=fnoise
   else
     fluct=.8d0*fluct+.2d0*fnoise
   endif

END SUBROUTINE updatefluctsum


!> Should we evaluate the translational force also with blocked atoms?
subroutine transforce(at,fxyz,sumx,sumy,sumz)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp),intent(in):: fxyz(3,at%astruct%nat)
  real(gp), intent(out) :: sumx,sumy,sumz
  integer :: iat

  !atomic_dot with one
  sumx=0._gp 
  sumy=0._gp 
  sumz=0._gp
  do iat=1,at%astruct%nat

     sumx=sumx+fxyz(1,iat) 
     sumy=sumy+fxyz(2,iat) 
     sumz=sumz+fxyz(3,iat)

  end do
END SUBROUTINE transforce


!> Should we evaluate the translational force also with blocked atoms?
subroutine transforce_forfluct(at,fxyz,sumx,sumy,sumz)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp),intent(in):: fxyz(3,at%astruct%nat)
  real(gp), intent(out) :: sumx,sumy,sumz
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  !atomic_dot with one
  sumx=0._gp 
  sumy=0._gp 
  sumz=0._gp
  do iat=1,at%astruct%nat

     call frozen_alpha(at%astruct%ifrztyp(iat),1,1.0_gp,alphax)
     call frozen_alpha(at%astruct%ifrztyp(iat),2,1.0_gp,alphay)
     call frozen_alpha(at%astruct%ifrztyp(iat),3,1.0_gp,alphaz)

     sumx=sumx+alphax*fxyz(1,iat) 
     sumy=sumy+alphay*fxyz(2,iat) 
     sumz=sumz+alphaz*fxyz(3,iat)

  end do
END SUBROUTINE transforce_forfluct


!> DIIS relax. Original source from ART from N. Mousseau.
!! Adaptations to BigDFT by D. Caliste.
!! WARNING: strten not minimized here
subroutine rundiis(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  logical, intent(out) :: fail
  !local variables
  character(len=*), parameter :: subname='rundiis'
  real(gp), dimension(:,:,:), allocatable  :: previous_forces
  real(gp), dimension(:,:,:), allocatable  :: previous_pos
  real(gp), dimension(:,:), allocatable :: product_matrix
  integer :: lter, maxter, i, i_err, n, nrhs, lwork, infocode, j, i_stat, i_all
  real(gp) :: fluct, fmax, fnrm,etotprev
  character(len = 4) :: fn4
  character(len = 40) :: comment
  ! Local variables for Lapack.
  integer, dimension(:), allocatable :: interchanges
  real(kind=8), dimension(:), allocatable :: work
  real(kind=8), dimension(:), allocatable :: solution
  real(kind=8), dimension(:,:), allocatable :: matrice
  !logical :: check
  integer :: check
  check=0
  fluct=0.0_gp

  ! We save pointers on data used to call bigdft() routine.
  allocate(previous_forces(3, outs%fdim, runObj%inputs%history+ndebug),stat=i_stat)
  call memocc(i_stat,previous_forces,'previous_forces',subname)
  allocate(previous_pos(3, RUNOBJ%ATOMS%astruct%NAT, runObj%inputs%history+ndebug),stat=i_stat)
  call memocc(i_stat,previous_pos,'previous_pos',subname)
  allocate(product_matrix(runObj%inputs%history, runObj%inputs%history+ndebug),stat=i_stat)
  call memocc(i_stat,product_matrix,'product_matrix',subname)


  ! Set to zero the arrays
  call to_zero(runObj%inputs%history**2,product_matrix)
  call to_zero(runObj%inputs%history*outs%fdim*3,previous_forces)
  call to_zero(runObj%inputs%history*runObj%atoms%astruct%nat*3,previous_pos)

  ! We set the first step and move to the second
  call vcopy(3 * outs%fdim, outs%fxyz (1,1), 1, previous_forces(1,1,1), 1)
  call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, previous_pos(1,1,1), 1)
  
  call axpy(3 * runObj%atoms%astruct%nat, runObj%inputs%betax, outs%fxyz(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
!!$  !always better to use the atomic_* routines to move atoms
!!$  !it performs modulo operation as well as constrained search
!!$  call atomic_axpy(at,x,runObj%inputs%betax,f,x)


  do lter = 1, runObj%inputs%ncount_cluster_x - 1

     maxter = min(lter, runObj%inputs%history)

     allocate(interchanges(maxter+1+ndebug),stat=i_stat)
     call memocc(i_stat,interchanges,'interchanges',subname)
     allocate(work((maxter+1) * (maxter+1)+ndebug),stat=i_stat)
     call memocc(i_stat,work,'work',subname)
     allocate(solution(maxter+1+ndebug),stat=i_stat)
     call memocc(i_stat,solution,'solution',subname)
     allocate(matrice(maxter+1, maxter+1))
     call memocc(i_stat,matrice,'matrice',subname)

     ! If lter is greater than maxvec, we move the previous solution up by one
     if (lter > maxter) then
        do i=2, maxter
           previous_forces(:,:,i-1) = previous_forces(:,:,i)
           previous_pos(:,:,i-1)    = previous_pos(:,:,i)
           do j=2, maxter
              product_matrix(i-1,j-1) = product_matrix(i,j)
           end do
        end do
     end if

     ! we first add the force to the previous_force vector and the
     ! position to the previous_pos vector
     call vcopy(3 * runObj%atoms%astruct%nat, outs%fxyz(1,1), 1, previous_forces(1,1,maxter), 1)
     call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, previous_pos(1,1,maxter), 1)

     ! And we add the scalar products to the matrix
     do i = 1, maxter
        product_matrix(i,maxter) = dot(3 * runObj%atoms%astruct%nat,previous_forces(1,1,i),1,outs%fxyz(1,1),1)
        product_matrix(maxter,i) = product_matrix(i,maxter)
     end do

     matrice(1:maxter,1:maxter) = product_matrix(1:maxter, 1:maxter)
     matrice(1:maxter,maxter+1) = 1.0d0
     matrice(maxter+1,1:maxter) = 1.0d0
     matrice(maxter+1,maxter+1) = 0.0d0

     solution(1:maxter) = 0.0d0
     solution(maxter+1) = 1.0d0

     ! We now need the routines from Lapack. We define a few values
     i_err = 0

     ! We call the routine for diagonalizing a tridiagonal  matrix
     n = maxter+1
     nrhs = 1    ! Number of right-hand arguments in B (Ax=B)
     lwork = n*n

     ! We prepare the upper triangular matrix for lapack
     !call dsysv('U',n, nrhs, matrice, n, interchanges, solution,n,work,lwork,i_err)
     call gesv(n, nrhs, matrice(1,1), n, interchanges(1), solution(1), n, i_err)

     i_all=-product(shape(interchanges))*kind(interchanges)
     deallocate(interchanges,stat=i_stat)
     call memocc(i_stat,i_all,'interchanges',subname)
     i_all=-product(shape(work))*kind(work)
     deallocate(work,stat=i_stat)
     call memocc(i_stat,i_all,'work',subname)
     i_all=-product(shape(matrice))*kind(matrice)
     deallocate(matrice,stat=i_stat)
     call memocc(i_stat,i_all,'matrice',subname)


     ! The solution that interests us is made of two parts
     call to_zero(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz)
     do i = 1, maxter
        call axpy(3 * runObj%atoms%astruct%nat, solution(i), previous_pos(1,1,i), 1, runObj%atoms%astruct%rxyz(1,1), 1)
     end do
!!$     !reput the modulo operation on the atoms
!!$     call atomic_axpy(at,x,0.0_gp,x,x)

     i_all=-product(shape(solution))*kind(solution)
     deallocate(solution,stat=i_stat)
     call memocc(i_stat,i_all,'solution',subname)

     runObj%inputs%inputPsiId=1
     etotprev=outs%energy

     call call_bigdft(runObj,outs,nproc,iproc,infocode)

!!$     if (iproc == 0) then
!!$        call transforce(at,f,sumx,sumy,sumz)
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
!!$     end if

     ncount_bigdft=ncount_bigdft+1

     call fnrmandforcemax(outs%fxyz,fnrm,fmax,outs%fdim)
     if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct) !n(m)
     call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,check) !n(m)

     if (iproc==0) then 
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,es11.3,3(1pe10.2),2x,i3)')  & 
          ncount_bigdft,lter,"GEOPT_DIIS",outs%energy,outs%energy-etotprev, &
          & fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,check

!        write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=', &
!             & fmax,'fnrm=',    fnrm    ,'fluct=', fluct
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'DIIS:fnrm= ',sqrt(fnrm)
        call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
             & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),forces=outs%fxyz)
     endif

     if(check > 5)then
        if (iproc==0) write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*runObj%inputs%frac_fluct,fluct
        if (iproc==0) write(16,*) 'DIIS converged'
        exit
     endif

     if(ncount_bigdft>runObj%inputs%ncount_cluster_x-1)  then 
      if (iproc==0)  &
           write(16,*) 'DIIS exited before the geometry optimization converged because more than ',& 
                            runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
      exit
     endif

     call axpy(3 * runObj%atoms%astruct%nat, runObj%inputs%betax, outs%fxyz(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
  end do

  i_all=-product(shape(previous_forces))*kind(previous_forces)
  deallocate(previous_forces,stat=i_stat)
  call memocc(i_stat,i_all,'previous_forces',subname)
  i_all=-product(shape(previous_pos))*kind(previous_pos)
  deallocate(previous_pos,stat=i_stat)
  call memocc(i_stat,i_all,'previous_pos',subname)
  i_all=-product(shape(product_matrix))*kind(product_matrix)
  deallocate(product_matrix,stat=i_stat)
  call memocc(i_stat,i_all,'product_matrix',subname)

  fail = (ncount_bigdft>runObj%inputs%ncount_cluster_x-1)
END SUBROUTINE rundiis


!> Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201 (2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University 
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine fire(runObj,outs,nproc,iproc,ncount_bigdft,fail) 
  use module_base
  use module_types
  use module_interfaces
  use minpar
  use yaml_output
  use communications_base

  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  logical, intent(inout) :: fail

  real(gp) :: fluct,fnrm
  real(gp) :: fmax,vmax,maxdiff
  integer :: check
  integer :: infocode,iat
  character(len=4) :: fn4
  character(len=40) :: comment

  !n(c) character(len=*), parameter :: subname='fire'

!Fire parameters:
  real(gp):: alpha,P,finc,fdec,falpha,alphastart,dt,dtmax,vnrm
  real(gp):: velcur(3*runObj%atoms%astruct%nat), velpred(3*runObj%atoms%astruct%nat),poscur(3*runObj%atoms%astruct%nat),&
       & pospred(3*runObj%atoms%astruct%nat),fcur(3*runObj%atoms%astruct%nat),fpred(3*runObj%atoms%astruct%nat),&
       & mass(3*runObj%atoms%astruct%nat)
  real(gp):: eprev,anoise !n(c) ecur
  integer:: Nmin,nstep,it

  fluct=0.0_gp
  check=0
  ! Set FIRE parameters
  Nmin=5
  finc=1.1_gp
  fdec=0.5_gp
  alphastart=0.25_gp
  anoise=1.e-8_gp

  alpha=alphastart
  falpha=0.99_gp
  nstep=1
  dt=runObj%inputs%dtinit

  dtmax=runObj%inputs%dtmax

  fail=.false.
  fnrm=1.e10_gp
  velcur=0.0_gp
  call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, poscur(1), 1)
  call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fcur(1), 1)
  mass=1.0_gp
  !n(c) ecur=etot
  eprev=0.0_gp

  Big_loop: do it=1,runObj%inputs%ncount_cluster_x-1
     do iat=1,3*runObj%atoms%astruct%nat
        pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5_gp*fcur(iat)/mass(iat)
     enddo

     runObj%inputs%inputPsiId=1
     call vcopy(3 * runObj%atoms%astruct%nat, pospred(1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
     call call_bigdft(runObj,outs,nproc,iproc,infocode)
     call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fpred(1), 1)
     ncount_bigdft=ncount_bigdft+1
     call fnrmandforcemax(fpred,fnrm,fmax,outs%fdim)
   !  call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,check) !n(m)

     do iat=1,3*runObj%atoms%astruct%nat
        velpred(iat)=velcur(iat)+0.5_gp*dt*(fpred(iat))/mass(iat)+0.5_gp*dt*fcur(iat)/mass(iat)
     enddo
     P=dot_product(fpred,velpred)
     call fnrmandforcemax(velpred,vnrm,vmax,runObj%atoms%astruct%nat)

     if (iproc == 0) then
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'FIRE:fnrm= ',sqrt(fnrm)
        call  write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4,&
             & outs%energy,pospred,runObj%atoms,trim(comment),forces=fpred)
     endif
     if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct) !n(m)

     if (iproc==0.and.parmin%verbosity > 0) then
         write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),  & 
         & 2x,a6,es8.2e1,2x,a3,es8.2e1,2x,a6,es9.2,2x,a6,I5,2x,a2,es9.2)') &
         & ncount_bigdft,it,"GEOPT_FIRE",outs%energy,outs%energy-eprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct, &
         & "alpha=",alpha, "dt=",dt, "vnrm=",sqrt(vnrm), "nstep=",nstep,"P=",P

         call yaml_open_map('Geometry')
            call yaml_map('Ncount_BigDFT',ncount_bigdft)
            call yaml_map('Geometry step',it)
            call yaml_map('Geometry Method','GEOPT_FIRE')
            call yaml_map('epred',(/ outs%energy,outs%energy-eprev /),fmt='(1pe21.14)')
            call geometry_output(fmax,fnrm,fluct)
            call yaml_map('Alpha', alpha, fmt='(es7.2e1)')
            call yaml_map('dt',dt, fmt='(es7.2e1)')
            call yaml_map('vnrm',sqrt(vnrm), fmt='(es8.2)')
            call yaml_map('nstep',nstep, fmt='(I5)')
            call yaml_map('P',P, fmt='(es9.2)')
            call geometry_output(fmax,fnrm,fluct)
         call yaml_close_map()
         !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2), & 
         !& 2x,a6,es7.2e1,2x,a3,es7.2e1,2x,a6,es8.2,2x,a6,I5,2x,a2,es9.2)') &
         !& ncount_bigdft,it,"GEOPT_FIRE",epred,epred-eprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct, &
         !& "alpha=",alpha, "dt=",dt, "vnrm=",sqrt(vnrm), "nstep=",nstep,"P=",P 
         !eprev=epred
     end if

     eprev=outs%energy

     call convcheck(fmax,fluct*runObj%inputs%frac_fluct, runObj%inputs%forcemax,check) !n(m)
     if (ncount_bigdft >= runObj%inputs%ncount_cluster_x-1) then
         !Too many iterations
         exit Big_loop
     end if

     if(check > 5) then
        if(iproc==0)  call yaml_map('Iterations when FIRE converged',it)
        !if(iproc==0)  write(16,'(a,i0,a)') "   FIRE converged in ",it," iterations"
        !Exit from the loop (the calculation is finished).
        exit Big_loop
     endif

!Update variables
     fcur=fpred
     poscur=pospred
!Normal verlet velocity update
!  velcur=velpred

!!FIRE Update
     call fnrmandforcemax(fpred,fnrm,fmax,outs%fdim)
     fnrm=sqrt(fnrm)
     call fnrmandforcemax(velpred,vnrm,vmax,runObj%atoms%astruct%nat)
     vnrm=sqrt(vnrm)
!Modified velocity update, suggested by Alireza
!  velcur(:)=(1.0_gp-alpha)*velpred(:)+fpred(:)*min(alpha*vnrm/fnrm,2.0_gp*runObj%inputs%betax)!alpha*fpred(:)/fnrm*vnrm
!Original FIRE velocitiy update
     velcur(:)=(1.0_gp-alpha)*velpred(:)+alpha*fpred(:)/fnrm*vnrm
     if(P > -anoise*vnrm .and. nstep > Nmin) then
        dt=min(dt*finc,dtmax)
!        alpha=max(alpha*falpha,0.1_gp) !Limit the decrease of alpha
        alpha=alpha*falpha
     elseif(P.le.-anoise*vnrm) then
        nstep=0
        dt=dt*fdec
        velcur=0.d0
        alpha=alphastart
     endif
     nstep=nstep+1

     ! Check velcur consistency.
     call check_array_consistency(maxdiff, nproc, velcur(1), &
          & 3 * runObj%atoms%astruct%nat, bigdft_mpi%mpi_comm)
     if (iproc==0 .and. maxdiff > epsilon(1.0_gp)) &
          call yaml_warning('Fire velocities not identical! '//&
          '(difference:'//trim(yaml_toa(maxdiff))//' )')

     !if (iproc==0) write(10,*) epred, vnrm*0.5d0
   end do Big_loop
        
! Output the final energy, atomic positions and forces
   call vcopy(3*runObj%atoms%astruct%nat, pospred(1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
   call vcopy(3*outs%fdim, fpred(1), 1, outs%fxyz(1,1), 1)

END SUBROUTINE fire


!> Display geometry quantities
subroutine geometry_output(fmax,fnrm,fluct)
   use module_base
   use yaml_output
   implicit none
   real(gp), intent(in) :: fmax    !< Maximal absolute value of atomic forces
   real(gp), intent(in) :: fnrm    !< Norm of atomic forces
   real(gp), intent(in) :: fluct   !< Fluctuation of atomic forces
   !write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=',fmax,'fnrm2=',fnrm,'fluct=', fluct
   call yaml_open_map('FORCES norm(Ha/Bohr)',flow=.true.)
      call yaml_map(' maxval',fmax,fmt='(1pe14.5)')
      call yaml_map('fnrm2',fnrm,fmt='(1pe14.5)')
      call yaml_map('fluct',fluct,fmt='(1pe14.5)')
   call yaml_close_map()
end subroutine geometry_output
