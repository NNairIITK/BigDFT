!> @file
!!  Routines to do geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>   Geometry optimization, parametrisation routine.
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


!>   Geometry optimization, parametrisation routine.
subroutine geopt_set_verbosity(verbosity_)
  use minpar
  implicit none
  !Arguments
  integer, intent(in) :: verbosity_
  parmin%verbosity = verbosity_
END SUBROUTINE geopt_set_verbosity


!>   Geometry optimization
subroutine geopt(nproc,iproc,pos,at,fxyz,strten,epot,rst,in,ncount_bigdft)
  use module_base
  use module_interfaces, except_this_one => geopt
  use module_types
  use yaml_output
  use minpar
  implicit none
  integer, intent(in) :: nproc,iproc
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: epot
  integer, intent(inout) :: ncount_bigdft
  real(gp), dimension(3*at%nat), intent(inout) :: pos
  real(gp), dimension(6), intent(inout) :: strten
  real(gp), dimension(3*at%nat), intent(inout) :: fxyz
  !local variables
  logical :: fail
  integer :: ibfgs
  character(len=6) :: outfile, fmt
  character*5 fn4
  character*40 comment
  !-------------------------------------------

  call geopt_init()
  if (iproc ==0 .and. parmin%verbosity > 0)  write(16,'(a)')  & 
     '# Geometry optimization log file, grep for GEOPT for consistent output'
  if (iproc ==0 .and. parmin%verbosity > 0) write(16,'(a)')  & 
      '# COUNT  IT  GEOPT_METHOD  ENERGY                 DIFF       FMAX       FNRM      FRAC*FLUC FLUC      ADD. INFO'

  if (iproc ==0 .and. parmin%verbosity > 0) write(* ,'(a)') & 
      '# COUNT  IT  GEOPT_METHOD  ENERGY                 DIFF       FMAX       FNRM      FRAC*FLUC FLUC      ADD. INFO'

  !assign the geometry optimisation method
  parmin%approach=in%geopt_approach

  strten=0 !not used for the moment

  epot=0.d0
  ncount_bigdft=0
  if (iproc == 0) then
     outfile = 'posout'
     if (trim(parmin%approach)=='AB6MD') outfile = 'posmd '
     fmt = "(i4.4)"
     if (trim(parmin%approach)=='AB6MD') fmt = '(i5.5)'
     write(fn4,fmt) ncount_bigdft
     write(comment,'(a)')'INITIAL CONFIGURATION '
     call write_atomic_file(trim(in%dir_output)//trim(outfile)//'_'//trim(fn4),epot,pos,at,trim(comment),forces=fxyz)
     write(*,'(a,1x,a)') ' Begin of minimization using ',parmin%approach
  end if

  if (trim(parmin%approach)=='LBFGS') then
  
     ibfgs=0
86   ibfgs=ibfgs+1
     if (iproc ==0) write(*,*) '# ENTERING LBFGS,ibfgs',ibfgs
     call lbfgsdriver(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)
     if (fail .and. ibfgs .lt. 5) goto 86

     if (fail) then
        if (iproc ==0) write(*,*) '# ENTERING CG after LBFGS failure'
        call conjgrad(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)
     end if


!  if(trim(parmin%approach)=='LBFGS') then
!
!     if (iproc ==0) write(*,*) '# ENTERING BFGS'
!
!     call bfgs(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)
!
!     if (fail) then 
!        if (iproc ==0) write(*,*) '# ENTERING CG after BFGS failure'
!        call conjgrad(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)
!     end if
!
  else if(trim(parmin%approach)=='BFGS' .or. trim(parmin%approach)=='PBFGS') then
     call bfgsdriver(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft)
  else if(trim(parmin%approach)=='SDCG') then

     if (iproc ==0) write(*,*) '# ENTERING CG'
!     call yaml_open_map('Geometry optimization')
     call conjgrad(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)
!     call yaml_close_map()

  else if(trim(parmin%approach)=='VSSD') then
 
     if (iproc ==0) write(*,*) '# ENTERING VSSD'
     call vstepsd(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)

  else if(trim(parmin%approach)=='FIRE') then

     if (iproc ==0) write(*,*) '# ENTERING FIRE'
     call fire(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft,fail)

  else if(trim(parmin%approach)=='DIIS') then
 
     if (iproc ==0) write(*,*) '# ENTERING DIIS'
     call rundiis(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)

  else if(trim(parmin%approach)=='AB6MD') then

     if (iproc ==0) write(*,*) '# ENTERING Molecular Dynamics (ABINIT implementation)'
     call ab6md(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)

  else
     write(*,*) 'ERROR: geometry optimization method undefined, exiting...',trim(parmin%approach)
     stop 
  endif
  if (iproc==0) write(*,'(a,1x,a)') 'End of minimization using ',parmin%approach

  if (iproc==0) call finaliseCompress()

END SUBROUTINE geopt

!>  Molecular Dynamics
subroutine ab6md(nproc,iproc,x,f,epot,at,rst,in,ncount_bigdft,fail)
  use module_base
  use module_types
  use scfloop_API
  use ab6_moldyn
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: epot
  real(gp), dimension(3*at%nat), intent(inout) :: x
  logical, intent(out) :: fail
  real(gp), dimension(3*at%nat), intent(inout) :: f
  !local variables
  !n(c) character(len=*), parameter :: subname='ab6md'
  ! 1 atomic mass unit, in electronic mass
  real(gp), parameter :: amu2emass=1.660538782d-27/9.10938215d-31
  integer :: nxfh, iat, idim, iexit
  integer, allocatable :: iatfix(:,:)
  real(gp) :: acell(3), rprim(3,3), symrel(3,3,1)
  real(gp), allocatable :: xfhist(:,:,:,:), amass(:), vel(:,:), xred(:,:), fred(:,:)

  ! We save pointers on data used to call bigdft() routine.
  if (ncount_bigdft == 0) then
     call scfloop_init(nproc, at, in, rst)
  end if

  ! Prepare the objects used by ABINIT.
  allocate(amass(at%nat))
  allocate(xfhist(3, at%nat + 4, 2, in%ncount_cluster_x+1))
  allocate(vel(3, at%nat))
  allocate(xred(3, at%nat))
  allocate(iatfix(3, at%nat))
  allocate(fred(3, at%nat))
  nxfh = 0
  !acell = (/ at%alat1, at%alat2, at%alat3 /)
  acell(1)=at%alat1
  acell(2)=at%alat2
  acell(3)=at%alat3
  rprim(:,:) = 0.0_gp
  rprim(1,1) = real(1, gp)
  rprim(2,2) = real(1, gp)
  rprim(3,3) = real(1, gp)
  symrel(:,:,1) = 0.0_gp
  symrel(1,1,1) = real(1, gp)
  symrel(2,2,1) = real(1, gp)
  symrel(3,3,1) = real(1, gp)
  do iat = 1, at%nat
     amass(iat) = amu2emass * at%amu(at%iatype(iat))
     do idim = 1, 3
        xred(idim, iat) = x((iat - 1) * 3 + idim) / acell(idim)
        fred(idim, iat) = - f((iat - 1) * 3 + idim) / acell(idim)
     end do
     if (at%ifrztyp(iat) == 0) then
        iatfix(:, iat) = 0
     else if (at%ifrztyp(iat) == 1) then
        iatfix(:, iat) = 1
     else if (at%ifrztyp(iat) == 2) then
        iatfix(:, iat) = 0
        iatfix(2, iat) = 1
     else if (at%ifrztyp(iat) == 3) then
        iatfix(:, iat) = 1
        iatfix(2, iat) = 0
     end if
  end do

  !read the velocities from input file, if present
  call read_velocities(iproc,'velocities.xyz',at,vel)
  !vel(:,:) = zero

  ! Call the ABINIT routine.
  ! currently, we force optcell == 0
  call moldyn(acell, amass, iproc, in%ncount_cluster_x+1, nxfh, at%nat, &
       & rprim, epot, iexit, &
       & 0, in%ionmov, in%ncount_cluster_x, in%dtion, in%noseinert, &
       & in%mditemp, in%mdftemp, in%friction, in%mdwall, in%nnos, &
       & in%qmass, in%bmass, in%vmass, iatfix, in%strtarget, &
       & in%strprecon, in%strfact, in%forcemax, &
       & 1, symrel, vel, xfhist, fred, xred)

  do iat = 1, at%nat, 1
     do idim = 1, 3, 1
        f(idim + 3 * (iat - 1)) = fred(idim, iat) * acell(idim)
     end do
  end do

  deallocate(fred)
  deallocate(iatfix)
  deallocate(xred)
  deallocate(vel)
  deallocate(amass)
  deallocate(xfhist)

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
  real(gp), intent(in):: ff(3,at%nat)
  real(gp), intent(out):: fnrm, fmax
  real(gp):: t1,t2,t3
  integer:: iat

!!!  t1=0._gp 
!!!  t2=0._gp 
!!!  t3=0._gp
  fmax=0._gp
  do iat=1,at%nat
     call frozen_alpha(at%ifrztyp(iat),1,ff(1,iat)**2,t1)
     call frozen_alpha(at%ifrztyp(iat),2,ff(2,iat)**2,t2)
     call frozen_alpha(at%ifrztyp(iat),3,ff(3,iat)**2,t3)
     fmax=max(fmax,sqrt(t1+t2+t3))
!!!     if (at%ifrztyp(iat) == 0) then
!!!        t1=t1+ff(1,iat)**2 
!!!        t2=t2+ff(2,iat)**2 
!!!        t3=t3+ff(3,iat)**2
        !in general fmax is measured with the inf norm
!!!     fmax=max(fmax,sqrt(ff(1,iat)**2+ff(2,iat)**2+ff(3,iat)**2))
        !fmax=max(fmax,abs(ff(1,iat)),abs(ff(2,iat)),abs(ff(3,iat)))
!!!     end if
  enddo

  !this is the norm of the forces of non-blocked atoms
  call atomic_dot(at,ff,ff,fnrm)
!!!  fnrm=t1+t2+t3
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
  real(gp),intent(in):: fxyz(3,at%nat)
  real(gp), intent(out) :: sumx,sumy,sumz
  integer :: iat

  !atomic_dot with one
  sumx=0._gp 
  sumy=0._gp 
  sumz=0._gp
  do iat=1,at%nat

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
  real(gp),intent(in):: fxyz(3,at%nat)
  real(gp), intent(out) :: sumx,sumy,sumz
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  !atomic_dot with one
  sumx=0._gp 
  sumy=0._gp 
  sumz=0._gp
  do iat=1,at%nat

     call frozen_alpha(at%ifrztyp(iat),1,1.0_gp,alphax)
     call frozen_alpha(at%ifrztyp(iat),2,1.0_gp,alphay)
     call frozen_alpha(at%ifrztyp(iat),3,1.0_gp,alphaz)

     sumx=sumx+alphax*fxyz(1,iat) 
     sumy=sumy+alphay*fxyz(2,iat) 
     sumz=sumz+alphaz*fxyz(3,iat)

  end do
END SUBROUTINE transforce_forfluct


!>  DIIS relax. Original source from ART from N. Mousseau.
!!  Adaptations to BigDFT by D. Caliste.
!!  WARNING: strten not minimized here
subroutine rundiis(nproc,iproc,x,f,epot,at,rst,in,ncount_bigdft,fail)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: epot
  real(gp), dimension(3*at%nat), intent(inout) :: x
  logical, intent(out) :: fail
  real(gp), dimension(3*at%nat), intent(inout) :: f
  !local variables
  character(len=*), parameter :: subname='rundiis'
  real(gp), dimension(6) :: strten
  real(gp), dimension(:,:), allocatable  :: previous_forces
  real(gp), dimension(:,:), allocatable  :: previous_pos
  real(gp), dimension(:,:), allocatable :: product_matrix
  integer :: lter, maxter, i, i_err, n, nrhs, lwork, infocode, j, i_stat, i_all
  real(gp) :: fluct, fmax, fnrm,fnoise,etotprev
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
  allocate(previous_forces(in%history, AT%NAT * 3+ndebug),stat=i_stat)
  call memocc(i_stat,previous_forces,'previous_forces',subname)
  allocate(previous_pos(in%history, AT%NAT * 3+ndebug),stat=i_stat)
  call memocc(i_stat,previous_pos,'previous_pos',subname)
  allocate(product_matrix(in%history, in%history+ndebug),stat=i_stat)
  call memocc(i_stat,product_matrix,'product_matrix',subname)


  !zet to zero the arrays
  call razero(in%history**2,product_matrix)
  call razero(in%history*at%nat*3,previous_forces)
  call razero(in%history*at%nat*3,previous_pos)

  ! We set the first step and move to the second
  previous_forces(1,:) = f(:)
  previous_pos(1,:) = x(:)

  x(:) = x(:) + in%betax * f(:)
!!$  !always better to use the atomic_* routines to move atoms
!!$  !it performs modulo operation as well as constrained search
!!$  call atomic_axpy(at,x,in%betax,f,x)


  do lter = 2, in%ncount_cluster_x

     maxter = min(lter, in%history)

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
           previous_forces(i-1,:) = previous_forces(i,:)
           previous_pos(i-1,:)    = previous_pos(i,:)
           do j=2, maxter
              product_matrix(i-1,j-1) = product_matrix(i,j)
           end do
        end do
     end if

     ! we first add the force to the previous_force vector and the
     ! position to the previous_pos vector
     previous_forces(maxter,:) = f(:)
     previous_pos(maxter,:) = x(:)

     ! And we add the scalar products to the matrix
     do i = 1, maxter
        product_matrix(i,maxter) = dot_product(previous_forces(i,:),f(:))
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
     call dsysv('U',n, nrhs, matrice, n, interchanges, solution,n,work,lwork,i_err)

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

     x(:) = 0.0d0
     do i = 1, maxter
        x(:) = x(:) + solution(i) * previous_pos(i,:)
     end do
!!$     !reput the modulo operation on the atoms
!!$     call atomic_axpy(at,x,0.0_gp,x,x)

     i_all=-product(shape(solution))*kind(solution)
     deallocate(solution,stat=i_stat)
     call memocc(i_stat,i_all,'solution',subname)

     in%inputPsiId=1
     etotprev=epot

     call call_bigdft(nproc,iproc,at,x,in,epot,f,strten,fnoise,rst,infocode)

!!$     if (iproc == 0) then
!!$        call transforce(at,f,sumx,sumy,sumz)
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
!!$     end if

     ncount_bigdft=ncount_bigdft+1

     call fnrmandforcemax(f,fnrm,fmax,at%nat)

     if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)

     call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,check) !n(m)

     if (iproc==0) then 
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,es11.3,3(1pe10.2),2x,i3)')  & 
          ncount_bigdft,lter,"GEOPT_DIIS",epot,epot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,check

!        write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=', &
!             & fmax,'fnrm=',    fnrm    ,'fluct=', fluct
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'DIIS:fnrm= ',sqrt(fnrm)
        call write_atomic_file(trim(in%dir_output)//'posout_'//fn4,epot,x,at,trim(comment),forces=f)
     endif

     if(check.gt.5)then
        if (iproc==0) write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*in%frac_fluct,fluct
        if (iproc==0) write(16,*) 'DIIS converged'
        exit
     endif

     if(ncount_bigdft>in%ncount_cluster_x-1)  then 
      if (iproc==0)  &
           write(16,*) 'DIIS exited before the geometry optimization converged because more than ',& 
                            in%ncount_cluster_x,' wavefunction optimizations were required'
      exit
     endif

     x(:) = x(:) + in%betax * f(:)
     !call atomic_axpy(at,x,in%betax,f,x)
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

  fail = (ncount_bigdft>in%ncount_cluster_x-1)
END SUBROUTINE rundiis

!! Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201 (2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University 
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine fire(nproc,iproc,rxyz,at,etot,fxyz,rst,in,ncount_bigdft,fail) 
  use module_base
  use module_types
  use module_interfaces
  use minpar

  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: etot
  real(gp), dimension(3*at%nat), intent(inout) :: rxyz
  logical, intent(inout) :: fail
  real(gp), dimension(3*at%nat), intent(inout) :: fxyz

  real(gp) :: fluct,fnrm,  fnoise
  real(gp) :: fmax,vmax
  integer :: check
  integer :: infocode,iat
  character(len=4) :: fn4
  character(len=40) :: comment
  real(gp), dimension(6) :: strten

  !n(c) character(len=*), parameter :: subname='fire'

!Fire parameters:
  real(gp):: alpha,P,finc,fdec,falpha,alphastart,dt,dtmax,vnrm
  real(gp):: velcur(3*at%nat), velpred(3*at%nat),poscur(3*at%nat),pospred(3*at%nat),fcur(3*at%nat),fpred(3*at%nat),mass(3*at%nat)
  real(gp):: epred,eprev,anoise !n(c) ecur
  integer:: Nmin,nstep,it

  fluct=0.0_gp
  check=0
!Set FIRE parameters
  Nmin=5
  finc=1.1_gp
  fdec=0.5_gp
  alphastart=0.25_gp
  anoise=1.e-8_gp


  alpha=alphastart
  falpha=0.99_gp
  nstep=1
  dt=in%dtinit

  dtmax=in%dtmax

  fail=.false.
  fnrm=1.e10_gp
  velcur=0.0_gp
  poscur=rxyz
  fcur=fxyz
  mass=1.0_gp
  !n(c) ecur=etot
  epred=etot


  Big_loop: do it=1,in%ncount_cluster_x-1
     do iat=1,3*at%nat
        pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5_gp*fcur(iat)/mass(iat)
     enddo

     in%inputPsiId=1
     call call_bigdft(nproc,iproc,at,pospred,in,epred,fpred,strten,fnoise,rst,infocode)
     ncount_bigdft=ncount_bigdft+1
     call fnrmandforcemax(fpred,fnrm,fmax,at%nat)
   !  call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,check) !n(m)

     do iat=1,3*at%nat
        velpred(iat)=velcur(iat)+0.5_gp*dt*(fpred(iat))/mass(iat)+0.5_gp*dt*fcur(iat)/mass(iat)
     enddo
     P=dot_product(fpred,velpred)
     call fnrmandforcemax(velpred,vnrm,vmax,at%nat)

     if (iproc == 0) then
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'FIRE:fnrm= ',sqrt(fnrm)
        call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,epred,pospred,at,trim(comment),forces=fpred)
     endif
     if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)
     if (iproc==0.and.parmin%verbosity > 0) & 
         write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),  & 
         &2x,a6,es7.2e1,2x,a3,es7.2e1,2x,a6,es8.2,2x,a6,I5,2x,a2,es9.2)') &
         &ncount_bigdft,it,"GEOPT_FIRE",epred,epred-eprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct, &
         &"alpha=",alpha, "dt=",dt, "vnrm=",sqrt(vnrm), "nstep=",nstep,"P=",P
     if (iproc==0.and.parmin%verbosity > 0) & 
         write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2), & 
         &2x,a6,es7.2e1,2x,a3,es7.2e1,2x,a6,es8.2,2x,a6,I5,2x,a2,es9.2)') &
         &ncount_bigdft,it,"GEOPT_FIRE",epred,epred-eprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct, &
         &"alpha=",alpha, "dt=",dt, "vnrm=",sqrt(vnrm), "nstep=",nstep,"P=",P 
         eprev=epred
     if (iproc==0.and.parmin%verbosity > 0) write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                             'FORCES norm(Ha/Bohr): maxval=',fmax,'fnrm2=',fnrm,'fluct=', fluct
     call convcheck(fmax,fluct*in%frac_fluct, in%forcemax,check) !n(m)
     if (ncount_bigdft >= in%ncount_cluster_x-1) then
         !Too many iterations
         exit Big_loop
     end if

     if(check.gt.5) then
        if(iproc==0)  write(16,'(a,i0,a)') "   FIRE converged in ",it," iterations"
        !Exit from the loop (the calculation is finished).
        exit Big_loop
     endif

!Update variables
     fcur=fpred
     poscur=pospred
!Normal verlet velocity update
!  velcur=velpred

!!FIRE Update
     call fnrmandforcemax(fpred,fnrm,fmax,at%nat)
     fnrm=sqrt(fnrm)
     call fnrmandforcemax(velpred,vnrm,vmax,at%nat)
     vnrm=sqrt(vnrm)
!Modified velocity update, suggested by Alireza
!  velcur(:)=(1.0_gp-alpha)*velpred(:)+fpred(:)*min(alpha*vnrm/fnrm,2.0_gp*in%betax)!alpha*fpred(:)/fnrm*vnrm
!Original FIRE velocitiy update
     velcur(:)=(1.0_gp-alpha)*velpred(:)+alpha*fpred(:)/fnrm*vnrm
     if(P.gt.-anoise*vnrm .and. nstep.gt.Nmin) then
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

     !if (iproc==0) write(10,*) epred, vnrm*0.5d0
   end do Big_loop


        
! Output the final energy, atomic positions and forces
   etot = epred
   rxyz = pospred
   fxyz = fpred

END SUBROUTINE fire



