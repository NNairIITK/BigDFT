!> @file
!!  Routines to do geometry optmization
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
  use bigdft_run
  use yaml_output
  use minpar
  implicit none
  !Arguments
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  integer, intent(in) :: nproc,iproc
  integer, intent(out) :: ncount_bigdft
  !local variables
  integer, parameter :: ugeopt = 16
  logical :: fail
  integer :: ibfgs,ierr
  character(len=6) :: outfile, fmt
  character(len=5) :: fn4
  character(len=40) :: comment
  character(len=60) :: filename
  !-------------------------------------------
  call f_routine(id='geopt')
  !Geometry Initialization
  call geopt_init()

  filename=trim(runObj%inputs%dir_output)//'geopt.mon'
  !open(unit=ugeopt,file=filename,status='unknown',position='append')
  !without istat this opening should crash 
  !do nothing if the unit is already opened
  if (iproc == 0) &
     call yaml_set_stream(unit=ugeopt,filename=trim(filename),tabbing=0,record_length=100,setdefault=.false.,istat=ierr)
  if (iproc ==0 ) call yaml_comment('Geopt file opened, name: '//trim(filename)//', timestamp: '//trim(yaml_date_and_time_toa()),&
       hfill='-',unit=ugeopt)
  !write(ugeopt,*) '----------------------------------------------------------------------------'

  if (iproc == 0 .and. parmin%verbosity > 0 .and. ierr /= 0)  &!write(ugeopt,'(a)')  & 
     call yaml_comment('Geometry optimization log file, grep for GEOPT for consistent output',unit=ugeopt)
  if (iproc == 0 .and. parmin%verbosity > 0) write(ugeopt,'(a)')  & 
      '# COUNT  IT  GEOPT_METHOD  ENERGY                 DIFF       FMAX       FNRM      FRAC*FLUC FLUC      ADD. INFO'

  !if (iproc ==0 .and. parmin%verbosity > 0) write(* ,'(a)') & 
  !    '# COUNT  IT  GEOPT_METHOD  ENERGY                 DIFF       FMAX       FNRM      FRAC*FLUC FLUC      ADD. INFO'

  !assign the geometry optmization method
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
     write(comment,'(a)')'INITIAL_CONFIGURATION '
     call bigdft_write_atomic_file(runObj,outs,&
          trim(outfile)//'_'//trim(fn4),trim(comment))
!!$     call write_atomic_file(trim(runObj%inputs%dir_output)//trim(outfile)//'_'//trim(fn4),&
!!$          & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms%astruct%ixyz_int,runObj%atoms,trim(comment),forces=outs%fxyz)
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
  case('SQNM','SBFGS')
   if(trim(adjustl(parmin%approach))=='SBFGS')then
     if (iproc==0) &
          call yaml_warning('Keyword SBFGS is deprecated and will be removed in a future release. Use SQNM instead.')
   endif 
   if (iproc ==0) call yaml_map('ENTERING SQNM',ncount_bigdft)
   call sqnm(runObj,outs,nproc,iproc,parmin%verbosity,ncount_bigdft,fail)
  case('DIIS')   
     if (iproc ==0) call yaml_map('ENTERING DIIS',ncount_bigdft)
     call rundiis(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case('AB6MD')
     if (iproc ==0) call yaml_map('ENTERING Molecular Dynamics (ABINIT implementation)',ncount_bigdft)
     call ab6md(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case('SOCK')
     if (iproc ==0) call yaml_map('ENTERING Socket communication',ncount_bigdft)
     call f2fslave(runObj,outs,nproc,iproc,ncount_bigdft,fail) 
  case('LOOP')
     if (iproc ==0) call yaml_map('ENTERING LOOP mode',ncount_bigdft)
     call loop(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  case default
     call f_err_throw('Geometry optimization method undefined ('//trim(parmin%approach)//')',&
          err_name='BIGDFT_RUNTIME_ERROR')
  end select

  if (iproc==0) call yaml_close_stream(unit=ugeopt)

  if (iproc==0) call yaml_map('End of minimization using ',parmin%approach)

  if (iproc==0) call finaliseCompress()

  call f_release_routine()
END SUBROUTINE geopt

!> Loop mode
subroutine loop(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  use module_base
  use bigdft_run
  use yaml_output
  use module_input_keys, only: inputpsiid_set_policy
  use public_enums, only: ENUM_MEMORY
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  logical, intent(out) :: fail
  !Local variables
  real(gp) :: fnrm, fmax, fluct
  integer :: check, infocode, it

  fail=.false.
  fluct=0.0_gp
  check=0
  it = 0

  do
     !calculate the max of the forces
     call fnrmandforcemax(outs%fxyz,fnrm,fmax,outs%fdim)
     if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct)

     call convcheck(fmax,fluct*runObj%inputs%frac_fluct, &
          & runObj%inputs%forcemax,check)
     if (ncount_bigdft >= runObj%inputs%ncount_cluster_x-1) then
        !Too many iterations
        fail = .true.
        return
     end if

     if(check > 0) then
        if(iproc==0)  call yaml_map('Iterations when LOOP converged',it)
        return
     endif

     !runObj%inputs%inputPsiId=1
     call inputpsiid_set_policy(ENUM_MEMORY, runObj%inputs%inputPsiId)
     call bigdft_state(runObj, outs, infocode)
     ncount_bigdft = ncount_bigdft + 1

     it = it + 1
  end do
END SUBROUTINE loop

!> Molecular Dynamics
subroutine ab6md(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  use module_base
  use bigdft_run
  use scfloop_API
  use m_ab6_moldyn, only: abi_moldyn
!  use module_interfaces, only: memocc
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
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
  integer, dimension(3,3,1) :: symrel
  integer :: nxfh, iat, idim, iexit, i_stat,i_all

  call f_routine(id='ab6md')

  ! We save pointers on data used to call bigdft() routine.
  if (ncount_bigdft == 0) then
     call scfloop_init(nproc, runObj)
  end if
  ! Prepare the objects used by ABINIT.
  amass = f_malloc(runObj%atoms%astruct%nat,id='amass')
  xfhist = f_malloc((/ 3, runObj%atoms%astruct%nat + 4, 2, runObj%inputs%ncount_cluster_x+1 /),id='xfhist')
  vel = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='vel')
  xred = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='xred')
  iatfix = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='iatfix')
  fred = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='fred')

  nxfh = 0
  !acell = (/ runObj%atoms%astruct%cell_dim(1), runObj%atoms%astruct%cell_dim(2), runObj%atoms%astruct%cell_dim(3) /)
  acell=runObj%atoms%astruct%cell_dim
  rprim(:,:) = 0.0_gp
  rprim(1,1) = real(1, gp)
  rprim(2,2) = real(1, gp)
  rprim(3,3) = real(1, gp)
  symrel(:,:,1) = 0
  symrel(1,1,1) = 1
  symrel(2,2,1) = 1
  symrel(3,3,1) = 1
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
  call abi_moldyn(acell, amass, iproc, runObj%inputs%ncount_cluster_x+1, nxfh, runObj%atoms%astruct%nat, &
       rprim, outs%energy, iexit, &
       0, runObj%inputs%ionmov, runObj%inputs%ncount_cluster_x, runObj%inputs%dtion, runObj%inputs%noseinert, &
       runObj%inputs%mditemp, runObj%inputs%mdftemp, runObj%inputs%friction, runObj%inputs%mdwall, runObj%inputs%nnos, &
       runObj%inputs%qmass, runObj%inputs%bmass, runObj%inputs%vmass, iatfix, runObj%inputs%strtarget, &
       runObj%inputs%strprecon, runObj%inputs%strfact, runObj%inputs%forcemax, &
       1, symrel, vel, xfhist, fred, xred)
  do iat = 1, runObj%atoms%astruct%nat, 1
     outs%fxyz(1, iat) = fred(1, iat) * acell(1)
     outs%fxyz(2, iat) = fred(2, iat) * acell(2)
     outs%fxyz(3, iat) = fred(3, iat) * acell(3)
  end do

  !De-Allocations
  call f_free(fred)
  call f_free(iatfix)
  call f_free(xred)
  call f_free(vel)
  call f_free(amass)
  call f_free(xfhist)

  call scfloop_finalise()
  fail = (iexit == 0)

  call f_release_routine()
END SUBROUTINE ab6md


!> MODIFIED version for refined time limit on restart of global.f90.
!! Only difference: Calls routine CPUtime(tt)
subroutine timeleft(tt)
  use module_base
  implicit none
  real(gp), intent(out) :: tt
  !local variables
  integer :: ierr,unt
  real(kind=4) :: tcpu
  real(gp) :: timelimit

  unt=55
  call f_open_file(unit=unt,file='CPUlimit',status='unknown')
  read(unt,*,iostat=ierr) timelimit ! in hours
  if(ierr/=0)timelimit=1d6
  call f_close(unt)
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


!> DIIS relax. Original source from ART from N. Mousseau.
!! Adaptations to BigDFT by D. Caliste.
!! WARNING: strten not minimized here
subroutine rundiis(runObj,outs,nproc,iproc,ncount_bigdft,fail)
  use module_base
  use bigdft_run
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  logical, intent(out) :: fail
  !local variables
  integer, parameter :: ugeopt=16
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
  type(f_tree) :: f_info
  check=0
  fluct=0.0_gp

  f_info=f_tree_new()
  ! We save pointers on data used to call bigdft() routine.
  ! Set to zero the arrays
  previous_forces = f_malloc0((/ 3, outs%fdim, runObj%inputs%history /),id='previous_forces')
  previous_pos = f_malloc0((/ 3, RUNOBJ%ATOMS%astruct%NAT, runObj%inputs%history /),id='previous_pos')
  product_matrix = f_malloc0((/ runObj%inputs%history, runObj%inputs%history /),id='product_matrix')




  ! We set the first step and move to the second
  call vcopy(3 * outs%fdim, outs%fxyz (1,1), 1, previous_forces(1,1,1), 1)
  call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, previous_pos(1,1,1), 1)
  
  call axpy(3 * runObj%atoms%astruct%nat, runObj%inputs%betax, outs%fxyz(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
!!$  !always better to use the atomic_* routines to move atoms
!!$  !it performs modulo operation as well as constrained search
!!$  call atomic_axpy(at,x,runObj%inputs%betax,f,x)


  do lter = 1, runObj%inputs%ncount_cluster_x - 1

     maxter = min(lter, runObj%inputs%history)

     interchanges = f_malloc(maxter+1,id='interchanges')
     work = f_malloc((maxter+1) * (maxter+1),id='work')
     solution = f_malloc(maxter+1,id='solution')
     matrice = f_malloc((/ maxter+1, maxter+1 /),id='matrice')

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

     call f_free(interchanges)
     call f_free(work)
     call f_free(matrice)


     ! The solution that interests us is made of two parts
     call f_zero(runObj%atoms%astruct%rxyz)
     do i = 1, maxter
        call axpy(3 * runObj%atoms%astruct%nat, solution(i), previous_pos(1,1,i), 1, runObj%atoms%astruct%rxyz(1,1), 1)
     end do
!!$     !reput the modulo operation on the atoms
!!$     call atomic_axpy(at,x,0.0_gp,x,x)

     call f_free(solution)

     !runObj%inputs%inputPsiId=1
     call bigdft_set_input_policy(INPUT_POLICY_MEMORY, runObj)
     etotprev=outs%energy

     call bigdft_state(runObj,outs,infocode)

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
        write(ugeopt,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,es11.3,3(1pe10.2),2x,i3)')  & 
          ncount_bigdft,lter,"GEOPT_DIIS",outs%energy,outs%energy-etotprev, &
          & fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,check

        call f_tree_push(f_info//'etot',     yaml_toa((/ outs%energy,outs%energy-etotprev /),fmt='(1pe21.14)'))
        call f_tree_push(f_info//'Forces' ,yaml_toa( (/ fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct /), fmt='(1pe10.2)'))
        call geometry_output('GEOPT_DIIS',ncount_bigdft,lter,fmax,fnrm,fluct,f_info)

        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'DIIS:fnrm= ',sqrt(fnrm)
        call bigdft_write_atomic_file(runObj,outs,&
             'posout_'//fn4,trim(comment))
!!$        call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
!!$             & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms%astruct%ixyz_int,runObj%atoms,trim(comment),forces=outs%fxyz)
     endif

     if(check > 5)then
        if (iproc==0) then
           write(ugeopt,'(1x,a,3(1x,1pe14.5))') '# fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*runObj%inputs%frac_fluct,fluct
           write(ugeopt,*) '# DIIS converged'
        end if
        exit
     endif

     if(ncount_bigdft>runObj%inputs%ncount_cluster_x-1)  then 
      if (iproc==0)  &
           write(ugeopt,*) 'DIIS exited before the geometry optimization converged because more than ',& 
                            runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
      exit
     endif

     call axpy(3 * runObj%atoms%astruct%nat, runObj%inputs%betax, outs%fxyz(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
  end do

  call f_free(previous_forces)
  call f_free(previous_pos)
  call f_free(product_matrix)
  call f_tree_free(f_info)

  fail = (ncount_bigdft>runObj%inputs%ncount_cluster_x-1)
END SUBROUTINE rundiis


!> Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201 (2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University 
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine fire(runObj,outs,nproc,iproc,ncount_bigdft,fail) 
  use module_base
  use bigdft_run
  use minpar
  use yaml_output
  use communications_base

  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  logical, intent(inout) :: fail
  !Local variables
  integer, parameter :: ugeopt=16
  real(gp) :: fluct,fnrm
  real(gp) :: fmax,vmax,maxdiff
  integer :: check
  integer :: infocode,iat
  character(len=4) :: fn4
  character(len=40) :: comment
  type(f_tree) :: info

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

  info=f_tree_new()

  Big_loop: do it=1,runObj%inputs%ncount_cluster_x-1
     do iat=1,3*runObj%atoms%astruct%nat
        pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5_gp*fcur(iat)/mass(iat)
     enddo

     !runObj%inputs%inputPsiId=1
     call bigdft_set_input_policy(INPUT_POLICY_MEMORY, runObj)
     call bigdft_set_rxyz(runObj,rxyz_add=pospred(1))
!!$     call vcopy(3 * runObj%atoms%astruct%nat, pospred(1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
     call bigdft_state(runObj,outs,infocode)
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
        call bigdft_write_atomic_file(runObj,outs,&
             'posout_'//fn4,trim(comment))
!!$        call  write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4,&
!!$             & outs%energy,pospred,runObj%atoms%astruct%ixyz_int,runObj%atoms,trim(comment),forces=fpred)
     endif
     if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct) !n(m)

     if (iproc==0.and.parmin%verbosity > 0) then
         write(ugeopt,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),  & 
         & 2x,a6,es8.2e1,2x,a3,es8.2e1,2x,a6,es9.2,2x,a6,I5,2x,a2,es9.2)') &
         & ncount_bigdft,it,"GEOPT_FIRE",outs%energy,outs%energy-eprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct, &
         & "alpha=",alpha, "dt=",dt, "vnrm=",sqrt(vnrm), "nstep=",nstep,"P=",P

         call f_tree_push(info//'epred',yaml_toa([outs%energy,outs%energy-eprev],fmt='(1pe21.14)'))
         call f_tree_push(info//'Alpha',yaml_toa(alpha, fmt='(es7.2e1)'))
         call f_tree_push(info//'dt'   ,yaml_toa(dt, fmt='(es7.2e1)'))
         call f_tree_push(info//'vnrm' ,yaml_toa(sqrt(vnrm), fmt='(es8.2)'))
         call f_tree_push(info//'nstep',yaml_toa(nstep, fmt='(I5)'))
         call f_tree_push(info//'P'    ,yaml_toa(P, fmt='(es9.2)'))
         call geometry_output('GEOPT_FIRE',ncount_bigdft,it,fmax,fnrm,fluct,info)
         
!!$         call yaml_mapping_open('Geometry')
!!$            call yaml_map('Ncount_BigDFT',ncount_bigdft) !universal
!!$            call yaml_map('Geometry step',it)
!!$            call yaml_map('Geometry Method','GEOPT_FIRE')
!!$            call yaml_map('epred',(/ outs%energy,outs%energy-eprev /),fmt='(1pe21.14)')
!!$            call yaml_map('Alpha', alpha, fmt='(es7.2e1)')
!!$            call yaml_map('dt',dt, fmt='(es7.2e1)')
!!$            call yaml_map('vnrm',sqrt(vnrm), fmt='(es8.2)')
!!$            call yaml_map('nstep',nstep, fmt='(I5)')
!!$            call yaml_map('P',P, fmt='(es9.2)')
!!$            call geometry_output(fmax,fnrm,fluct)
!!$         call yaml_mapping_close()
     end if

     eprev=outs%energy

     call convcheck(fmax,fluct*runObj%inputs%frac_fluct, runObj%inputs%forcemax,check) !n(m)
     if (ncount_bigdft >= runObj%inputs%ncount_cluster_x-1) then
         !Too many iterations
         exit Big_loop
     end if

     if(check > 5) then
        if(iproc==0)  call yaml_map('Iterations when FIRE converged',it)
        !if(iproc==0)  write(ugeopt,'(a,i0,a)') "   FIRE converged in ",it," iterations"
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

!!$     ! Check velcur consistency (debug).
!!$     maxdiff=mpimaxdiff(velcur,root=0,comm=bigdft_mpi%mpi_comm)
!!$     !call check_array_consistency(maxdiff, nproc, velcur, bigdft_mpi%mpi_comm)
!!$     if (iproc==0 .and. maxdiff > epsilon(1.0_gp)) &
!!$          call yaml_warning('Fire velocities were not identical (broadcasting)! '//&
!!$          '(difference:'//trim(yaml_toa(maxdiff))//' )')

     !if (iproc==0) write(10,*) epred, vnrm*0.5d0
   end do Big_loop

   call f_tree_free(info)
! Output the final energy, atomic positions and forces
   call vcopy(3*runObj%atoms%astruct%nat, pospred(1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
   call vcopy(3*outs%fdim, fpred(1), 1, outs%fxyz(1,1), 1)

END SUBROUTINE fire


!> Display geometry quantities
subroutine geometry_output(method,ncount_bigdft,it,fmax,fnrm,fluct,extra_tree)
   use module_base
   use yaml_output
   implicit none
   character(len=*), intent(in) :: method !<description of the geometry
   integer, intent(in) :: ncount_bigdft !<number of external calls
   integer, intent(in) :: it !< number of the iteration
   real(gp), intent(in) :: fmax    !< Maximal absolute value of atomic forces
   real(gp), intent(in) :: fnrm    !< Norm of atomic forces
   real(gp), intent(in) :: fluct   !< Fluctuation of atomic forces
   type(f_tree), intent(in) :: extra_tree !<extra information on the geometry step

   call yaml_mapping_open('Geometry')
     call yaml_map('Geometry Method',method)
     call yaml_map('Ncount_BigDFT',ncount_bigdft) !universal
     call yaml_map('Geometry step',it)
     call yaml_mapping_open('FORCES norm(Ha/Bohr)',flow=.true.)
        call yaml_map(' maxval',fmax,fmt='(1pe14.5)')
        call yaml_map('fnrm2',fnrm,fmt='(1pe14.5)')
        call yaml_map('fluct',fluct,fmt='(1pe14.5)')
     call yaml_mapping_close()
     call f_tree_dump(extra_tree)
   call yaml_mapping_close()
end subroutine geometry_output

!************************************************************************************
!Socket communication
subroutine f2fslave(runObj,outs,nproc,iproc,ncount_bigdft,fail) 
  use module_base
  use bigdft_run
  use yaml_output
  use communications_base
  use f90sockets, only: open_socket, writebuffer, readbuffer
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  logical, intent(inout) :: fail
! Test fortran-to-fortran socket communication
! JMS, Mar.2015
  integer,parameter    :: MSGLEN=12
  ! this should be the IP address or fully-qualified host name on which the slave should connect
  ! in order to use UNIX sockets, this should be the name of the socket.
  ! these (and port in case of TCP sockets) should match the address and port on which the master is running
  character(len=1032)  :: host='127.0.0.1'
  integer              :: inet = 1   ! this has to be one to use a TCP/IP socket, zero to have a UNIX socket
  integer              :: port=21211 ! this is the port number (only used for TCP/IP sockets) 
  integer              :: socket
  character(len=MSGLEN):: header
  character,dimension(MSGLEN):: header_arr
  integer              :: nat,repid
  !Controlling variables
  CHARACTER*1024 :: msg
  INTEGER :: nmsg
  INTEGER :: nat_get
  LOGICAL :: isinit=.false., hasdata=.false., firststep=.true.
  !Lattice vectors and others  
  REAL*8 :: latvec(3,3), latvec_inv(3,3), strmat(3,3),etot,strten(6)
  REAL*8 :: fcart(3,runObj%atoms%astruct%nat),pos(3,runObj%atoms%astruct%nat),rxyz(3,runObj%atoms%astruct%nat)
  REAL*8,allocatable :: send_array(:), get_array(:)
  integer:: infocode,i,PsiId,str_index



  host=trim(runObj%inputs%sockhost)
  inet=runObj%inputs%sockinet
  port=runObj%inputs%sockport

  if(iproc==0) then
    host = TRIM(host)//achar(0)
    call open_socket( socket, inet, port, host )
  endif
  if (bigdft_mpi%nproc > 1)  call mpibcast(port,1,comm=bigdft_mpi%mpi_comm)
  
  do    ! receive-send iteration
    if (iproc==0) call readbuffer(socket, header, MSGLEN)
    if (iproc==0) call str2arr(header,header_arr,MSGLEN)
    if (bigdft_mpi%nproc > 1) call mpibcast(header_arr,comm=bigdft_mpi%mpi_comm)
    call arr2str(header,header_arr,MSGLEN)

    if (iproc==0) write(*,'(i6,a,a)') iproc, ' # SOCKET SLAVE: header received ',trim(header)
      if (trim(header) == "STATUS") then
         if(iproc==0) call send_status(header, MSGLEN, isinit)
      else if (trim(header) == "INIT") then
         if(iproc==0) call get_init   (header, MSGLEN, repid, isinit)
         isinit=.true.
         if (bigdft_mpi%nproc > 1)  call mpibcast(repid,1,comm=bigdft_mpi%mpi_comm)
!Although the number of atoms and the arrays assiciated with that side should be
!allocate already, we allocate fcart and pos here         
         if(iproc==0.and.runObj%atoms%astruct%nat.ne.repid) &
              call f_err_throw("The number of atoms on the master and slave side should be the same")
!Check if cell should be reset, in this case forget previous WF
         if(iproc==0) then
            str_index = index(msg(1:msglen),"CRESET",.true.)
            if(str_index.gt.0) then
               write(*,*) " MESSAGE: CRESET"
               PsiId=0
            else
               PsiId=1
            endif
          endif
          if (bigdft_mpi%nproc > 1)  call mpibcast(PsiId,1,comm=bigdft_mpi%mpi_comm)
          select case(PsiId)
          case(0)
             call bigdft_set_input_policy(INPUT_POLICY_SCRATCH,runObj)
          case(1)
             call bigdft_set_input_policy(INPUT_POLICY_MEMORY,runObj)
          end select
      else if (trim(header) == "POSDATA") then
         if(iproc==0 ) call get_data(pos,latvec,runObj%atoms%astruct%nat,nat_get);
         if (bigdft_mpi%nproc > 1) call mpibcast(pos,comm=bigdft_mpi%mpi_comm)
         if (bigdft_mpi%nproc > 1) call mpibcast(latvec,comm=bigdft_mpi%mpi_comm)
!Convert the units of positions and lattice vectors in a format the bigdft can understand
         call rxyz_int2cart(latvec,pos,rxyz,runObj%atoms%astruct%nat)
         call bigdft_set_rxyz(runObj,rxyz_add=rxyz(1,1))
!!!Compute the forces and stress here!!!
         call bigdft_state(runObj,outs,infocode)
!Get Stress
         call vcopy(6,outs%strten (1), 1, strten(1) ,1)
         strmat=0.d0
         strmat(1,1)=strten(1)
         strmat(2,2)=strten(2)
         strmat(3,3)=strten(3)
         strmat(2,1)=strten(6);strmat(1,2)=strmat(2,1)
         strmat(3,1)=strten(5);strmat(1,3)=strmat(3,1)
         strmat(3,2)=strten(4);strmat(2,3)=strmat(3,2)
!Get Force
         call vcopy(3*outs%fdim, outs%fxyz (1,1), 1, fcart(1,1), 1)
!Get Energy
         etot=outs%energy
!Now I have the data
         hasdata=.true.
      else if (trim(header)=="GETFORCE") then
         nmsg = 0
         if(iproc==0) call send_data(etot,fcart,strmat,runObj%atoms%astruct%nat,nmsg,msg)
         isinit = .false. ! resets init so that we will get replica index again at next step!
         hasdata= .false.
    elseif (trim(header)=="STOP") then
      exit
    elseif (trim(header)=="WAIT") then
      cycle
    endif
  enddo

contains
  subroutine send_status(header, MSGLEN, isinit)
  !Report the status to the master
  implicit none
  integer:: MSGLEN
  character(MSGLEN):: header
  logical:: isinit
            if (hasdata) then
               header="HAVEDATA    "
            else if (isinit) then
               header="READY       "
            else 
               header="NEEDINIT    "
            endif
            call writebuffer(socket,header,MSGLEN)               
            write(*,'(a,a)')   " # SOCKET SLAVE: header sent ", trim(header)
  end subroutine
  subroutine get_init(header, MSGLEN, repid, isinit)
  !Get initiallization string plus a repid
  implicit none
  integer:: MSGLEN,repid
  character(MSGLEN):: header
  logical:: isinit
         write(*,'(a)')   " # SOCKET SLAVE: initiallizing... "
         call readbuffer(socket, repid)        ! actually this reads the replica id         
         call readbuffer(socket, nmsg)         ! length of parameter string -- ignored at present!
         call readbuffer(socket, msg, nmsg)    ! the actual message
         write(*,'(a,a)') " # SOCKET SLAVE: initiallization string ", msg(1:nmsg)
         isinit=.true.
  end subroutine
  subroutine get_data(pos,latvec,nat,nat_get)
  ! Receives the positions & the cell data
  implicit none
  integer:: nat,nat_get,i
  real(8):: pos(3,nat),pos_cart(3,nat),latvec(3,3)
  real(8),allocatable:: get_array(:)
! first reads cell and the number of atoms
         write(*,'(a)')   " # SOCKET SLAVE: waiting for positions "
         allocate(get_array(9))
         call readbuffer(socket, get_array , 9)
         latvec = transpose(RESHAPE(get_array, (/3,3/)))       !cell vector      
         call readbuffer(socket, get_array, 9)
         latvec_inv = transpose(RESHAPE(get_array, (/3,3/)))   !inverse cell vector
         deallocate(get_array)
         call readbuffer(socket, nat_get)                      !number of atoms
         if (nat.ne.nat_get) &
              call f_err_throw("Received NAT not the same as the local NAT")
         allocate(get_array(3*nat_get))
         call readbuffer(socket, get_array, nat_get*3)
         pos_cart = RESHAPE(get_array, (/ 3 , nat /) ) 
         call rxyz_cart2int(latvec,pos,pos_cart,nat)
         deallocate(get_array)
         write(21,*) 'lat'
         write(21,*)  latvec(:,1)  
         write(21,*)  latvec(:,2)  
         write(21,*)  latvec(:,3)  
         write(21,*) 'pos'
         do i=1,nat
            write(21,*) pos(:,i)
         enddo

         write(*,'(a)')   " # SOCKET SLAVE: received positions "
  end subroutine
  subroutine send_data(etot,fcart,strmat,nat,nmsg,msg)
  ! Sends the energy, forces and stresses
  implicit none
  integer::nat,nmsg
  real(8):: fcart(3,nat),strmat(3,3),etot
  real(8),allocatable:: get_array(:)
  CHARACTER*1024 :: msg
 ! communicates energy info back to master
         write(*,'(a)')   " # SOCKET SLAVE: sending energy, forces and stress "
         call writebuffer(socket,"FORCEREADY  ",MSGLEN)         
         call writebuffer(socket,etot)
         allocate(send_array(3*nat))
         send_array=reshape(fcart,(/3*nat/))
         call writebuffer(socket,nat)            
         call writebuffer(socket,send_array,3*nat)
         deallocate(send_array)
         allocate(send_array(9))
         send_array=reshape(strmat,(/9/))
         call writebuffer(socket,send_array,9)
         deallocate(send_array)
         ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
         ! trajectory file. this is useful if you want to return additional information, e.g.
         ! atomic charges, wannier centres, etc. one must return the number of characters, then
         ! the string. here we just send back zero characters.
         call writebuffer(socket,nmsg)
         if(nmsg.gt.0) then
                 call writebuffer(socket,msg,nmsg)
         endif
  end subroutine
  subroutine str2arr(string,strarr,n)
  implicit none
  integer:: n,i
  character(n)::string
  character(1)::strarr(n)
  do i=1,n
    strarr(i)=string(i:i)
  enddo
  end subroutine
  subroutine arr2str(string,strarr,n)
  implicit none
  integer:: n,i
  character(n)::string
  character(1)::strarr(n)
  do i=1,n
    string(i:i)=strarr(i)
  enddo
  end subroutine
end subroutine f2fslave

!************************************************************************************

 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int

!************************************************************************************

 subroutine rxyz_int2cart(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3)
 integer:: nat,iat
 do iat=1,nat
  rxyzcart(:,iat)=matmul(latvec,rxyzint(:,iat))
 enddo
 end subroutine rxyz_int2cart

!************************************************************************************


 subroutine invertmat(mat,matinv,n)
 implicit none
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8),allocatable   :: WORK(:)
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
 integer               :: LDWORK
 !Here only for a 3*3 matrix
 if (n==3) then
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
 &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
 div=1.d0/div
      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 else
 !General n*n matrix 
 matinv=mat
 allocate(WORK(n))
 call  DGETRF( n, n, matinv, n, IPIV, INFO )
 if (info.ne.0) stop "Error in DGETRF"
 LDWORK=-1
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 LDWORK=int(WORK(1))
 deallocate(WORK)
 allocate(WORK(LDWORK))
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 if (info.ne.0) stop "Error in DGETRI"
 endif
 end subroutine

!************************************************************************************
