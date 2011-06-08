!> @file
!!  Routines to do geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define the type parameterminimization
module minpar
  implicit none

  type parameterminimization
     !>general parameters for all methods
     character (len=10) :: approach
     integer :: iter
     integer :: iflag
     integer :: history
     !>parameters for print information
     integer :: verbosity
     integer :: MSAVE
     integer :: MP
     integer :: LP
     integer :: MAXFEV
     integer :: FINSTEP
     double precision :: ALPHA 
     double precision :: GTOL
     double precision :: XTOL
     double precision :: FTOL
     double precision :: STPMIN
     double precision :: STPMAX
     logical :: DIAGCO
     logical :: IWRITE
  end type parameterminimization

  type(parameterminimization) :: parmin

end module minpar


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
subroutine geopt(nproc,iproc,pos,at,fxyz,epot,rst,in,ncount_bigdft)
  use module_base
  use module_interfaces, except_this_one => geopt
  use module_types
  use minpar
  implicit none
  integer, intent(in) :: nproc,iproc
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: epot
  integer, intent(inout) :: ncount_bigdft
  real(gp), dimension(3*at%nat), intent(inout) :: pos
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

  epot=0.d0
  ncount_bigdft=0
  if (iproc == 0) then
     outfile = 'posout'
     if (trim(parmin%approach)=='AB6MD') outfile = 'posmd '
     fmt = "(i4.4)"
     if (trim(parmin%approach)=='AB6MD') fmt = '(i5.5)'
     write(fn4,fmt) ncount_bigdft
     write(comment,'(a)')'INITIAL CONFIGURATION '
     call write_atomic_file(trim(outfile)//'_'//trim(fn4),epot,pos,at,trim(comment))
     write(*,'(a,1x,a)') ' Begin of minimization using ',parmin%approach
  end if

  if (trim(parmin%approach)=='BFGS') then
  
     ibfgs=0
86   ibfgs=ibfgs+1
     if (iproc ==0) write(*,*) '# ENTERING BFGS,ibfgs',ibfgs
     call lbfgsdriver(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)
     if (fail .and. ibfgs .lt. 5) goto 86

     if (fail) then
        if (iproc ==0) write(*,*) '# ENTERING CG after BFGS failure'
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
  else if(trim(parmin%approach)=='SDCG') then

     if (iproc ==0) write(*,*) '# ENTERING CG'
     call conjgrad(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)

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

     if (iproc ==0) write(*,*) '# ENTERING Molecular Dynamic (ABINIT implementation)'
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
  real(gp), dimension(3*at%nat), intent(out) :: f
  !local variables
  character(len=*), parameter :: subname='ab6md'
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


!>  Conjugate gradient method
subroutine conjgrad(nproc,iproc,rxyz,at,etot,fxyz,rst,in,ncount_bigdft)
  use module_base
  use module_types
  use minpar
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  real(gp), intent(inout) :: etot
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: fxyz
  !local variables
  real(gp) ::fnoise
  character(len=*), parameter :: subname='conjgrad'  
  integer :: nfail,it,iat,i_all,i_stat,infocode, nitsd
  real(gp) :: anoise,fluct,avbeta,avnum,fnrm,etotprev,beta0,beta
  real(gp) :: y0,y1,tt,oben1,oben2,oben,unten,rlambda,tetot,fmax,tmp!,eprev
  real(gp), dimension(:,:), allocatable :: tpos,gpf,hh
!  logical::check
  integer::check
  character*4 fn4
  character*40 comment

  check=0
  allocate(tpos(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,tpos,'tpos',subname)
  allocate(gpf(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gpf,'gpf',subname)
  allocate(hh(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,hh,'hh',subname)

  anoise=1.e-4_gp
  fluct=0._gp 

  avbeta=0._gp
  avnum=0._gp
  nfail=0
  nitsd=500


  !start with a steepest descent algorithm
  call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,&
       fnrm,fnoise,in,in%forcemax,nitsd,fluct)
  if (ncount_bigdft >= in%ncount_cluster_x) then
      if (iproc==0 .and. parmin%verbosity > 0) &
           & write(16,*) 'SDCG exited before the geometry optimization converged because more than ',&
           in%ncount_cluster_x,' wavefunction optimizations were required'
     call close_and_deallocate
     return
  end if

  !calculate the max of the forces
  call fnrmandforcemax(fxyz,tmp,fmax,at%nat)

  !control whether the convergence criterion is reached after SD
  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
  if (check.gt.5) then
     if (iproc.eq.0) write(16,*) 'Converged before entering CG',iproc
     call close_and_deallocate
     return
  endif

  redo_cg: do
     do iat=1,at%nat
        hh(1,iat)=fxyz(1,iat)
        hh(2,iat)=fxyz(2,iat)
        hh(3,iat)=fxyz(3,iat)
     end do

     beta0=4._gp*in%betax
     it=0
     loop_cg: do
        it=it+1

        !C line minimize along hh ----
        !call atomic_axpy(at,rxyz,beta0,hh,tpos)
        tpos=rxyz+beta0*hh

        in%inputPsiId=1
        in%output_grid=0
        in%output_wf=.false.

        call call_bigdft(nproc,iproc,at,tpos,in,tetot,gpf,fnoise,rst,infocode)
!!$        if (iproc == 0) then
!!$           call transforce(at,gpf,sumx,sumy,sumz)
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
!!$        end if
        ncount_bigdft=ncount_bigdft+1


        !C projection of gradients at beta=0 and beta onto hh
        y0=0._gp
        y1=0._gp
        do iat=1,at%nat
           y0=y0+fxyz(1,iat)*hh(1,iat)+fxyz(2,iat)*hh(2,iat)+fxyz(3,iat)*hh(3,iat)
           y1=y1+gpf(1,iat)*hh(1,iat)+gpf(2,iat)*hh(2,iat)+gpf(3,iat)*hh(3,iat)
        end do
        tt=y0/(y0-y1)
!        if (iproc.eq.0) then
!           if (parmin%verbosity > 0) &
!                & write(16,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
!           write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
!        end if

        beta=beta0*max(min(tt,2._gp),-.25_gp)
        
        tpos=rxyz

        !call atomic_axpy(at,rxyz,beta,hh,rxyz)
        rxyz=rxyz+beta*hh

        avbeta=avbeta+beta/in%betax
        avnum=avnum+1._gp
        !if (iproc ==0)print *,'beta,avbeta,avnum',beta,avbeta,avnum 
        !C new gradient
        do iat=1,at%nat
           gpf(1,iat)=fxyz(1,iat)
           gpf(2,iat)=fxyz(2,iat)
           gpf(3,iat)=fxyz(3,iat)
        end do

        etotprev=etot
        call call_bigdft(nproc,iproc,at,rxyz,in,etot,fxyz,fnoise,rst,infocode)
!!$        if (iproc == 0) then
!!$           call transforce(at,fxyz,sumx,sumy,sumz)
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
!!$        end if
        ncount_bigdft=ncount_bigdft+1
        !if the energy goes up (a small tolerance of anoise is allowed)
        !switch back to SD
        if (etot > etotprev+anoise) then

           if (iproc.eq.0 .and. parmin%verbosity > 0) &
                & write(16,'(a,i5,2(1pe20.12))') 'switching back to SD:etot,etotprev',it,etot,etotprev
           if (iproc.eq.0) write(*,'(a,i5,2(1pe20.12))') ' switching back to SD:etot,etotprev',it,etot,etotprev
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,&
                fnrm,fnoise,in,in%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(fxyz,tmp,fmax,at%nat)
           call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
           if(check.gt.5) then
              if (iproc.eq.0) write(16,*) 'Converged in switch back SD',iproc
              call close_and_deallocate
              return
           endif
           cycle redo_cg
        endif

        if (iproc==0) then 
           write(fn4,'(i4.4)') ncount_bigdft
           write(comment,'(a,1pe10.3)')'CONJG:fnrm= ',sqrt(fnrm)
           call  write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))
        endif


        !if (iproc.eq.0) write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'CG ',ncount_bigdft,etot,sqrt(fnrm)

        if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)

!!$        call atomic_dot(at,gpf,gpf,unten)
!!$        call atomic_dot(at,gpf,fxyz,oben1)
!!$        call atomic_dot(at,fxyz,fxyz,oben2)
        unten=dot(3*at%nat,gpf(1,1),1,gpf(1,1),1)
        oben1=dot(3*at%nat,gpf(1,1),1,fxyz(1,1),1)
        oben2=dot(3*at%nat,fxyz(1,1),1,fxyz(1,1),1)

        oben=oben2-oben1

!!!        obenx=0._gp
!!!        obeny=0._gp
!!!        obenz=0._gp
!!!        unten=0._gp
!!!        do iat=1,at%nat
!!!           obenx=obenx+(fxyz(1,iat)-gpf(1,iat))*fxyz(1,iat)
!!!           obeny=obeny+(fxyz(2,iat)-gpf(2,iat))*fxyz(2,iat)
!!!           obenz=obenz+(fxyz(3,iat)-gpf(3,iat))*fxyz(3,iat)
!!!           unten=unten+gpf(1,iat)**2+gpf(2,iat)**2+gpf(3,iat)**2
!!!        end do

        call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
        if (iproc == 0) then
           if (parmin%verbosity > 0) then
!              write(16,'(i5,1x,e12.5,1x,e21.14,a,1x,e9.2)')it,sqrt(fnrm),etot,' GEOPT CG ',beta/in%betax
!              write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*in%frac_fluct,fluct
           write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
                ncount_bigdft,it,"GEOPT_CG  ",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"b/b0=",beta/in%betax
           write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
                ncount_bigdft,it,"GEOPT_CG  ",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"b/b0=",beta/in%betax
           end if
           write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                'FORCES norm(Ha/Bohr): maxval=',    fmax,'fnrm=',    fnrm   , 'fluct=',fluct
        end if

        call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

        if(check.gt.5) exit loop_cg

        if (ncount_bigdft.gt.in%ncount_cluster_x) then 
           if (iproc==0)  write(16,*) 'SDCG exited before the geometry optimization converged because more than ',&
                            in%ncount_cluster_x,' wavefunction optimizations were required'
           if (iproc.eq.0) write(*,*) 'ncount_bigdft in CG',ncount_bigdft
           exit loop_cg
        endif


        if(iproc==0)call timeleft(tt)
        call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_stat)
        if(tt<0)then
           call close_and_deallocate
           return
        endif

        !if no convergence is reached after 500 CG steps
        !switch back to SD
        if (it == 500) then
           if (iproc == 0) then
              if (parmin%verbosity > 0) write(16,*) &
                   'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
              write(*,*) 'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
           end if
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,&
                fnrm,fnoise,in,in%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(fxyz,tmp,fmax,at%nat)

           call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
           if(check.gt.5) then
              if (iproc == 0) write(16,*) 'Converged in back up SD',iproc
              call close_and_deallocate
              return
           endif

           nfail=nfail+1
           if (nfail >= 100) stop 'too many failures of CONJG'
           cycle redo_cg
        endif

        !rlambda=(obenx+obeny+obenz)/unten
        rlambda=oben/unten
        !call atomic_axpy_forces(at,fxyz,rlambda,hh,hh)
        do iat=1,at%nat
           hh(1,iat)=fxyz(1,iat)+rlambda*hh(1,iat)
           hh(2,iat)=fxyz(2,iat)+rlambda*hh(2,iat)
           hh(3,iat)=fxyz(3,iat)+rlambda*hh(3,iat)
        end do
     end do loop_cg
     exit redo_cg
  end do redo_cg

  !!        write(6,*) 'CG finished',it,fnrm,etot
  if (iproc == 0) then
     if (parmin%verbosity > 0) &
          & write(16,'(1x,a,f8.5,i5)') 'average CG stepsize in terms of betax',avbeta/avnum,iproc
     write(*,'(1x,a,f8.5,i5)') 'average CG stepsize in terms of betax',avbeta/avnum,iproc
  end if

 call close_and_deallocate

contains

  subroutine close_and_deallocate
    use module_base
    implicit none
    !    Close the file
    !close(unit=16)
    i_all=-product(shape(tpos))*kind(tpos)
    deallocate(tpos,stat=i_stat)
    call memocc(i_stat,i_all,'tpos',subname)
    i_all=-product(shape(gpf))*kind(gpf)
    deallocate(gpf,stat=i_stat)
    call memocc(i_stat,i_all,'gpf',subname)
    i_all=-product(shape(hh))*kind(hh)
    deallocate(hh,stat=i_stat)
    call memocc(i_stat,i_all,'hh',subname)
  END SUBROUTINE close_and_deallocate

END SUBROUTINE conjgrad


!>  Steepest descent method
subroutine steepdes(nproc,iproc,at,rxyz,etot,ff,rst,ncount_bigdft,&
     fnrm,fnoise,in,forcemax_sw,nitsd,fluct)
  use module_base
  use module_types
  use minpar
  !use module_interfaces
  implicit none
  integer, intent(in) :: nproc,iproc,nitsd
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  integer, intent(inout) :: ncount_bigdft
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  real(gp), intent(out) :: fnrm,fnoise
  real(gp), intent(inout) :: etot
  real(gp), intent(inout) :: fluct
  real(gp), dimension(3,at%nat), intent(out) ::ff
  real(gp), intent(in)::forcemax_sw
  !local variables
  character(len=*), parameter :: subname='steepdes'
  logical :: care,move_this_coordinate
  integer :: nsatur,iat,itot,itsd,i_stat,i_all,infocode,nbeqbx,i,ixyz,nr
  real(gp) :: etotitm2,fnrmitm2,etotitm1,fnrmitm1,anoise
  real(gp) :: fmax,de1,de2,df1,df2,beta,etotprev
  real(gp), allocatable, dimension(:,:) :: tpos
  character*4 fn4
  character*40 comment

  allocate(tpos(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,tpos,'tpos',subname)

  etotprev=etot
  anoise=0.e-4_gp
  fluct=0._gp

  beta=in%betax
  care=.true.
  nbeqbx=0
  nsatur=0
  etotitm2=1.e100_gp
  fnrmitm2=1.e100_gp
  etotitm1=1.e100_gp
  fnrmitm1=1.e100_gp

  do iat=1,at%nat
     tpos(1,iat)=rxyz(1,iat)
     tpos(2,iat)=rxyz(2,iat)
     tpos(3,iat)=rxyz(3,iat)
  end do

  itot=0

  nr=0
  do i=1,3*at%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(at%ifrztyp(iat),ixyz)) nr=nr+1
  enddo


  redo_sd: do
     if (ncount_bigdft.gt.in%ncount_cluster_x) then 
        if (iproc.eq.0) then
           write(*,'(a,i6,i6,i6)') 'SD FINISHED because ncount_bigdft > ncount_cluster_x',iproc,ncount_bigdft,in%ncount_cluster_x
           if (parmin%verbosity > 0) then
              write(16,'(a,i6,i6,i6)') &
                   'SD FINISHED because ncount_bigdft > ncount_cluster_x',iproc,ncount_bigdft,in%ncount_cluster_x
              write(16,'(a,i6,a)') 'SDCG exited before the geometry optimization converged because more than ', & 
                   in%ncount_cluster_x,' wavefunction optimizations were required'
           end if
        end if

        i_all=-product(shape(tpos))*kind(tpos)
        deallocate(tpos,stat=i_stat)
        call memocc(i_stat,i_all,'tpos',subname)
        return
     endif

     itsd=0
     loop_sd: do
        itsd=itsd+1
        itot=itot+1

        in%inputPsiId=1
        in%output_grid=0
        in%output_wf=.false.
        call call_bigdft(nproc,iproc,at,rxyz,in,etot,ff,fnoise,rst,infocode)
!!$        if (iproc == 0) then
!!$           call transforce(at,ff,sumx,sumy,sumz)
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
!!$        end if
        ncount_bigdft=ncount_bigdft+1

        !if the energy goes up (a small tolerance is allowed by anoise)
        !reduce the value of beta
        !this procedure stops in the case beta is much too small compared with the initial one
        if (care .and. etot > etotitm1+anoise) then
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do
           beta=.5_gp*beta
           if (iproc == 0) write(16,'(a,1x,e9.2,1x,i5,2(1x,e21.14))') &
                'SD reset, beta,itsd,etot,etotitm1= ',beta,itsd,etot,etotitm1
           if (beta <= 1.d-1*in%betax) then
              if (iproc == 0) write(16,*) &
                   'beta getting too small, do not care anymore if energy goes up'
              care=.false.
           endif
           cycle redo_sd
        endif

        call fnrmandforcemax(ff,fnrm,fmax,at%nat)
        if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)


        !first and second derivatives of the energy and of the norm of the forces
        !(in units of beta steps)
        de1=etot-etotitm1
        de2=etot-2._gp*etotitm1+etotitm2
        df1=fnrm-fnrmitm1
        df2=fnrm-2._gp*fnrmitm1+fnrmitm2

        if (iproc.eq.0) then
           if (parmin%verbosity > 0) &
                & write(16,'(a,6(1x,e10.3),1x,i2)') 'fmax, fnrm/fnrmitm1, de1<0 , de2>0 , df1<0 , df2>0 ,nsatur',  & 
                fmax, fnrm/fnrmitm1,de1,de2,df1,df2,nsatur
           write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=',    fmax,'fnrm=',    fnrm    ,'fluct=', fluct
        end if

        !control whether we are in a situation in which SD do not change things too much
        if (care .and. itsd >= 3 .and. beta == in%betax .and. &
             df1 < anoise .and. &                              !forces are decreasing
             de1 > -.1_gp .and. de1 < anoise .and. &            !energy slowly decreasing
!             fmax <= .03_gp .and. fnrm/fnrmitm1 > .50_gp .and. &   !norm of forces is saturating
             fnrm <= max(2.5d0*1.d-3*sqrt(real(nr,8)),0.008d0) .and. fnrm/fnrmitm1 > .50_gp .and. &   !norm of forces is saturating
             de2 > -2._gp*anoise .and. df2 > -2._gp*anoise) then !close to a local minimum (E&F)
           nsatur=nsatur+1
        else
           nsatur=0
        endif

        if (iproc == 0) then 
!           if (parmin%verbosity > 0) &
!                & write(16,'(i5,1x,e12.5,1x,e21.14,a)') itsd,sqrt(fnrm),etot,' GEOPT SD '
           write(fn4,'(i4.4)') ncount_bigdft 
           write(comment,'(a,1pe10.3)')'SD:fnrm= ',sqrt(fnrm)
           call write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))

           !write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'SD ',ncount_bigdft,etot,sqrt(fnrm)
        end if


        if (iproc==0 .and. parmin%verbosity > 0) then
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,I2)') &
             ncount_bigdft,itsd,"GEOPT_SD  ",etot, etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
             "b/b0=",beta/in%betax,"nsat=",nsatur
        write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,I2)') &
        &ncount_bigdft,itsd,"GEOPT_SD  ",etot, etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct, & 
        &"b/b0=",beta/in%betax,"nsat=",nsatur
        etotprev=etot 
!           write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',&
!                fnrm,fluct*in%frac_fluct,fluct
        end if

        !exit statements

        !the norm of the forces is of the same order of the fluctuations
        !or we are much too close to the fixed point from the SD viewpoint
        !or the forces are below the fixed tolerance

        if (iproc==0) then
            if (nsatur > 5  ) write(*,'(a,1x,i0)')&
                 'SD EXIT because nsatur > 5 ',nsatur 
            if (fmax < forcemax_sw ) write(*,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < forcemax_sw ', fmax , forcemax_sw
            if (fmax < fluct*in%frac_fluct ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < fluctuation ' ,fmax , fluct*in%frac_fluct
            if (fmax < fluct*in%forcemax ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < forcemax ' ,fmax , in%forcemax
        endif

        if ( nsatur > 5 .or. fmax < max(forcemax_sw,in%forcemax,fluct*in%frac_fluct))  exit loop_sd

        !maximum number of allowed iterations reached
        if (ncount_bigdft > in%ncount_cluster_x) then 
           if (iproc==0)  write(16,*) 'SDCG exited before the geometry optimization converged because more than ',& 
                            in%ncount_cluster_x,' wavefunction optimizations were required'
           if (iproc == 0) write(*,*) 'ncount_bigdft in SD2',ncount_bigdft,in%ncount_cluster_x
           exit loop_sd
        endif

        !maximum number of allowed SD steps reached, locally or globally
        if (itsd >= nitsd) then 
           if (iproc.eq.0) write(16,'(a,i5,1x,e10.3,1x,e21.14)') &
                'SD: NO CONVERGENCE:itsd,fnrm2,etot',itsd,fnrm,etot
           exit loop_sd
        endif
        if (itot >= nitsd) then
           if (iproc.eq.0) write(16,'(a,i5,i5,1x,e10.3,1x,e21.14)') &
                'SD: NO CONVERGENCE:itsd,itot,fnrm2,etot:',itsd,itot,fnrm,etot
           exit loop_sd
        endif

        etotitm2=etotitm1
        etotitm1=etot
        fnrmitm2=fnrmitm1
        fnrmitm1=fnrm
        
        beta=min(1.2_gp*beta,in%betax)
!!!        if (beta /= in%betax) nbeqbx=0
        if (beta == in%betax) then 
           !     if (iproc.eq.0) write(16,*) 'beta=betax'
           care=.true.
!!!           !if beta=betax since too many iterations (say 5),
!!!           !then betax can be increased
!!!           nbeqbx=nbeqbx+1
!!!           if (nbeqbx == 5) then
!!!              in%betax=1.2_gp*in%betax
!!!              nbeqbx=0
!!!           end if
        endif
!        if (iproc == 0 .and. parmin%verbosity > 0) write(16,*) 'beta=',beta

        tpos=rxyz
        !call atomic_axpy(at,rxyz,beta,ff,rxyz)
        call axpy(3*at%nat,beta,ff(1,1),1,rxyz(1,1),1)

!!!        do iat=1,at%nat
!!!           tpos(1,iat)=rxyz(1,iat)
!!!           tpos(2,iat)=rxyz(2,iat)
!!!           tpos(3,iat)=rxyz(3,iat)
!!!           if ( .not. at%lfrztyp(iat)) then
!!!              if (at%geocode == 'P') then
!!!                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*ff(1,iat),at%alat1)
!!!                 rxyz(2,iat)=modulo(rxyz(2,iat)+beta*ff(2,iat),at%alat2)
!!!                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*ff(3,iat),at%alat3)
!!!              else if (at%geocode == 'S') then
!!!                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*ff(1,iat),at%alat1)
!!!                 rxyz(2,iat)=rxyz(2,iat)+beta*ff(2,iat)
!!!                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*ff(3,iat),at%alat3)
!!!              else
!!!                 rxyz(1,iat)=rxyz(1,iat)+beta*ff(1,iat)
!!!                 rxyz(2,iat)=rxyz(2,iat)+beta*ff(2,iat)
!!!                 rxyz(3,iat)=rxyz(3,iat)+beta*ff(3,iat)
!!!              end if
!!!           end if
!!!        end do
     end do loop_sd
     exit redo_sd
  end do redo_sd

  if (iproc == 0 .and. parmin%verbosity > 0) write(16,*) 'SD FINISHED',iproc

  i_all=-product(shape(tpos))*kind(tpos)
  deallocate(tpos,stat=i_stat)
  call memocc(i_stat,i_all,'tpos',subname)

END SUBROUTINE steepdes


!> Variable step steepest descent
subroutine vstepsd(nproc,iproc,wpos,at,etot,ff,rst,in,ncount_bigdft)
  use module_base
  use module_types
  use minpar
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  real(gp), intent(inout) :: etot
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), dimension(3,at%nat), intent(inout) :: wpos
  real(gp), dimension(3,at%nat), intent(out) :: ff
  !local variables
  real(gp) ::fnoise
  character(len=*), parameter :: subname='vstepsd'  
  integer :: iat,i_all,i_stat,infocode, nitsd,itsd
  real(gp) :: anoise,fluct,fnrm,fnrmold,beta,betaxx,betalast,betalastold
  real(gp) :: etotold,fmax,scpr,curv,tt,etotprev
  real(gp), dimension(:,:), allocatable :: posold,ffold
  logical :: reset!,check
  integer :: check
  character(len=4) :: fn4
  character(len=40) :: comment

  check=0
  etotprev=etot
  allocate(posold(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,posold,'posold',subname)
  allocate(ffold(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,ffold,'ffold',subname)

  anoise=1.e-4_gp
  fluct=0._gp

  beta=in%betax
  itsd=0
  in%inputPsiId=1
  in%output_grid=0
  in%output_wf=.false.
  call call_bigdft(nproc,iproc,at,wpos,in,etotold,ffold,fnoise,rst,infocode)
  call fnrmandforcemax(ffold,fnrm,fmax,at%nat)   
  if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)
  if (iproc == 0) then
     if (parmin%verbosity > 0)   write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
     &ncount_bigdft,itsd,"GEOPT_VSSD",etotold,etotold-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"beta=",beta
     if (parmin%verbosity > 0)   write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
     &ncount_bigdft,itsd,"GEOPT_VSSD",etotold,etotold-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"beta=",beta
     etotprev=etotold
!!$     call transforce(at,ffold,sumx,sumy,sumz)                         
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  

     write(fn4,'(i4.4)') ncount_bigdft
     write(comment,'(a,1pe10.3)')'Initial VSSD:fnrm= ',sqrt(fnrm)
     call  write_atomic_file('posout_'//fn4,etotold,wpos,at,trim(comment))
!     if (parmin%verbosity > 0) &
!          & write(16,'(1x,e12.5,1x,e21.14,a,e10.3)')sqrt(fnrm),etotold,' GEOPT VSSD ',beta
  end if

  ncount_bigdft=ncount_bigdft+1


  fnrmold=0.0_gp
  do iat=1,at%nat
     fnrmold=fnrmold+ffold(1,iat)**2+ffold(2,iat)**2+ffold(3,iat)**2
     posold(1,iat)=wpos(1,iat)
     posold(2,iat)=wpos(2,iat)
     posold(3,iat)=wpos(3,iat)
     wpos(1,iat)=wpos(1,iat)+beta*ffold(1,iat)
     wpos(2,iat)=wpos(2,iat)+beta*ffold(2,iat)
     wpos(3,iat)=wpos(3,iat)+beta*ffold(3,iat)
  enddo
  !call atomic_dot(at,ffold,ffold,fnrmold)
  !call atomic_axpy(at,wpos,beta,ffold,wpos)
  betaxx=1.d100
  reset=.true.
  betalastold=in%betax

  nitsd=1000
  loop_ntsd: do itsd=1,nitsd
     in%inputPsiId=1
     in%output_grid=0
     in%output_wf=.false.
     call call_bigdft(nproc,iproc,at,wpos,in,etot,ff,fnoise,rst,infocode)
!!$     if (iproc == 0) then                                        
!!$        call transforce(at,ff,sumx,sumy,sumz)                         
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
!!$     end if
     ncount_bigdft=ncount_bigdft+1
     fnrm=0.d0
     scpr=0.d0
     do iat=1,at%nat
        fnrm=fnrm+ff(1,iat)**2+ff(2,iat)**2+ff(3,iat)**2
        scpr=scpr+ffold(1,iat)*ff(1,iat)+ffold(2,iat)*ff(2,iat)+ffold(3,iat)*ff(3,iat)
     enddo
!!$     call  atomic_dot(at,ff,ff,fnrm)
!!$     call  atomic_dot(at,ff,ffold,scpr)
     curv=(fnrmold-scpr)/(beta*fnrmold)
     betalast=.5d0/curv
     if (betalast.gt.0.d0) betaxx=min(betaxx,1.5d0*betalast)
     call fnrmandforcemax(ff,fnrm,fmax,at%nat)   
     if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)
!     if (iproc==0) write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',fnrm,fluct*in%frac_fluct,fluct
     call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
     if (check.gt.5) exit loop_ntsd
     if (ncount_bigdft >= in%ncount_cluster_x) then 
        if (iproc==0)  write(16,*) 'VSSD exited before the geometry optimization converged because more than ',& 
             in%ncount_cluster_x,' wavefunction optimizations were required'
        exit loop_ntsd
     endif


     if (etot > etotold) then
        reset=.true.
        beta=in%betax
        if (iproc == 0) write(16,*) 'new positions rejected, reduced beta',beta
        do iat=1,at%nat
           wpos(1,iat)=posold(1,iat)+beta*ffold(1,iat)
           wpos(2,iat)=posold(2,iat)+beta*ffold(2,iat)
           wpos(3,iat)=posold(3,iat)+beta*ffold(3,iat)
        enddo
        !call atomic_axpy(at,posold,beta,ffold,wpos)
     else
        reset=.false.
        if (betalast.gt.0.d0) then
           beta=max(min(beta*1.5d0,betalast),in%betax)
        else
           beta=1.25d0*beta
        endif

        if (iproc == 0) then
           write(fn4,'(i4.4)') ncount_bigdft-1
           write(comment,'(a,1pe10.3)')'VSSD:fnrm= ',sqrt(fnrm)
           call  write_atomic_file('posout_'//fn4,etot,wpos,at,trim(comment))
        endif

        do iat=1,at%nat
           posold(1,iat)=wpos(1,iat)
           posold(2,iat)=wpos(2,iat)
           posold(3,iat)=wpos(3,iat)
           wpos(1,iat)=wpos(1,iat)+beta*ff(1,iat)
           wpos(2,iat)=wpos(2,iat)+beta*ff(2,iat)
           wpos(3,iat)=wpos(3,iat)+beta*ff(3,iat)
           ffold(1,iat)=ff(1,iat)
           ffold(2,iat)=ff(2,iat)
           ffold(3,iat)=ff(3,iat)
        enddo
        !call atomic_axpy(at,wpos,beta,ff,wpos)
        etotold=etot
        fnrmold=fnrm
        betalastold=betalast
     endif
  
     if (iproc == 0.and.parmin%verbosity > 0) & 
          write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
          ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
          "beta=",beta,"last beta=",betalast
     if (iproc == 0.and.parmin%verbosity > 0) & 
          write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
          ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
     &"beta=",beta,"last beta=",betalast
     etotprev=etot
!     if (iproc == 0) write(16,'(i5,1x,e12.5,1x,e21.14,a,e10.3,1x,e10.3)') itsd,sqrt(fnrm),etot,' GEOPT VSSD ',beta,betalast
     if(iproc==0)call timeleft(tt)
     call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_stat)
     if(tt<0) exit loop_ntsd


  enddo loop_ntsd
  if (iproc == 0.and.parmin%verbosity > 0) & 
       write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
       ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
       "beta=",beta,"last beta=",betalast
  if (iproc == 0.and.parmin%verbosity > 0) & 
     &write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
     &ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
     &"beta=",beta,"last beta=",betalast
  if (iproc == 0 .and. itsd == nitsd+1) &
       write(16,'(a,i5,e9.2,e18.10,e9.2)') '---- SD FAILED  TO CONVERGE'
  if (iproc == 0) then
     if (parmin%verbosity > 0) then
        write(16,'(a,i5,e9.2,e18.10)') 'variable stepsize SD FINISHED,iter, force norm,energy',itsd,sqrt(fnrm),etot
        write(16,'(a,e9.2)') 'suggested value for stepsize:', betaxx
     end if
     write(fn4,'(i4.4)') ncount_bigdft-1
     write(comment,'(a,1pe10.3)')'VSSD:fnrm= ',sqrt(fnrm)
     call  write_atomic_file('posout_'//fn4,etot,wpos,at,trim(comment))
  endif


  i_all=-product(shape(posold))*kind(posold)
  deallocate(posold,stat=i_stat)
  call memocc(i_stat,i_all,'posold',subname)
  i_all=-product(shape(ffold))*kind(ffold)
  deallocate(ffold,stat=i_stat)
  call memocc(i_stat,i_all,'ffold',subname)


END SUBROUTINE vstepsd


subroutine convcheck(fnrm,fmax,fluctfrac_fluct,forcemax,check)
  use module_base
  implicit none
  real(gp), intent(in):: fnrm, fmax, fluctfrac_fluct,forcemax
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


subroutine updatefluctsum(nat,fnoise,fluct)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nat
  real(gp),intent(in):: fnoise
  real(gp),intent(inout):: fluct

   if (fluct.eq.0.d0) then
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
subroutine rundiis(nproc,iproc,x,f,epot,at,rst,in,ncount_bigdft,fail)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: epot
  real(gp), dimension(3*at%nat), intent(inout) :: x
  logical, intent(out) :: fail
  real(gp), dimension(3*at%nat), intent(out) :: f
  !local variables
  character(len=*), parameter :: subname='rundiis'
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

     call call_bigdft(nproc,iproc,at,x,in,epot,f,fnoise,rst,infocode)

!!$     if (iproc == 0) then
!!$        call transforce(at,f,sumx,sumy,sumz)
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
!!$     end if

     ncount_bigdft=ncount_bigdft+1

     call fnrmandforcemax(f,fnrm,fmax,at%nat)

     if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)

     if (iproc==0) then 
     write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,2(1pe11.3),3(1pe10.2))')  & 
          ncount_bigdft,lter,"GEOPT_DIIS",epot,epot-etotprev,fmax,sqrt(fnrm),fnrm,fluct*in%frac_fluct,fluct

!        write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=', &
!             & fmax,'fnrm=',    fnrm    ,'fluct=', fluct
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'DIIS:fnrm= ',sqrt(fnrm)
        call write_atomic_file('posout_'//fn4,epot,x,at,trim(comment))
     endif

     call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

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


!>   Driver for the LBFGS routine found on the Nocedal Homepage
!!   The subroutines have only been modified slightly, so a VIMDIFF will show all modifications!
!!   This is helpfull when we are looking for the source of problems during BFGS runs
subroutine lbfgsdriver(nproc,iproc,rxyz,fxyz,etot,at,rst,in,ncount_bigdft,fail) 
  use module_base
  use module_types
!  use par_driver
  use minpar
  implicit none
!  type(driverparameters)::par
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: etot
  real(gp), dimension(3*at%nat), intent(inout) :: rxyz
  logical, intent(out) :: fail
  real(gp), dimension(3*at%nat), intent(out) :: fxyz

  real(gp), dimension(3*at%nat):: txyz, sxyz
  real(gp) :: fluct,fnrm, fnoise
  real(gp) :: fmax
!  logical :: check
  integer :: check
  integer :: infocode,i,ixyz,iat,nitsd
  real(gp) :: fnormmax_sw,etotprev
  character(len=4) :: fn4
  character(len=40) :: comment
  logical :: move_this_coordinate

  integer:: n,nr,ndim
  integer:: NWORK
  real(gp),allocatable:: X(:),G(:),DIAG(:),W(:)
  real(gp):: F,EPS!,XTOL,GTOL,,STPMIN,STPMAX
  real(gp), dimension(3*at%nat) :: rxyz0,rxyzwrite
  INTEGER:: IPRINT(2),IFLAG,ICALL,M
  character(len=*), parameter :: subname='bfgs'
  integer :: i_stat,i_all

  check=0

!  call init_driver(par)     !Initialize the parameters
  parmin%finstep=0
  parmin%alpha=1.d0
  fail=.false.
  fnrm=1.d10
  nitsd=10!500                 !Maximum number of SD steps before entering BFGS
  fnormmax_sw=in%forcemax!1.e-2_gp      !SD till the max force comp is less than this value
  
  !Dummy variables
  txyz=0._gp
  sxyz=0._gp
  
  if (iproc==0)    write(*,*) 'Maximum number of SD steps used in the beginning: ',nitsd

  call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,fnrm,fnoise,in,&
       fnormmax_sw,nitsd,fluct)
  etotprev=etot
  rxyz0=rxyz     !Save initial positions, since the unconstrained degrees of freedom will be updated upon them
  rxyzwrite=rxyz
  call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
  !call fnrmandforcemax(fxyz,fnrm,fmax,at)
  !check if the convergence is reached after SD
  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

  if (check.gt.5) then
     if (iproc.eq.0) write(*,*) 'Converged before entering BFGS'
     return
  endif


  !Make a list of all degrees of freedom that should be passed to bfgs
  n=3*at%nat
  nr=0
  do i=1,3*at%nat
     iat=(i-1)/3+1
     ixyz=mod(i-1,3)+1
     if(move_this_coordinate(at%ifrztyp(iat),ixyz)) nr=nr+1
  enddo
  if(iproc==0) write(*,*) 'DOF: n,nr ',n,nr
     NDIM=nr
     NWORK=NDIM*(2*parmin%MSAVE +1)+2*parmin%MSAVE
      
     allocate(X(NDIM),stat=i_stat)
     call memocc(i_stat,X,'X',subname)
     allocate(G(NDIM),stat=i_stat)
     call memocc(i_stat,G,'G',subname)
     allocate(DIAG(NDIM),stat=i_stat)
     call memocc(i_stat,DIAG,'DIAG',subname)
     allocate(W(NWORK),stat=i_stat)
     call memocc(i_stat,W,'W',subname)

     call atomic_copymoving_forward(at,n,rxyz,nr,X)

     N=nr
     M=parmin%MSAVE
     IPRINT(1)= 1
     IPRINT(2)= 0
     F=etot
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.

     EPS=0.0_gp
     ICALL=0
     IFLAG=0

 20   CONTINUE
        if (parmin%IWRITE) then
           if (iproc == 0) then
              write(fn4,'(i4.4)') ncount_bigdft
              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
              call  write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))
           endif
           parmin%IWRITE=.false.
        endif
        rxyzwrite=rxyz

        if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)
   if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) & 
              &write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
              &ncount_bigdft,ICALL,"GEOPT_BFGS",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
              &,"BFGS-it=",parmin%finstep,"alpha=",parmin%alpha
   if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) & 
              & write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
              &ncount_bigdft,ICALL,"GEOPT_BFGS",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
              &,"BFGS-it=",parmin%finstep,"alpha=",parmin%alpha
              etotprev=etot
              if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                           'FORCES norm(Ha/Bohr): maxval=',fmax,'fnrm2=',fnrm,'fluct=', fluct
              call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
              if (ncount_bigdft >= in%ncount_cluster_x) goto 50
              close(16)
              open(unit=16,file='geopt.mon',status='unknown',position='APPEND')

      if(check.gt.5) then
         if(iproc==0)  write(16,'(a,i0,a)') "   BFGS converged in ",ICALL," iterations"
         if (iproc == 0) then
            write(fn4,'(i4.4)') ncount_bigdft
            write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
            call  write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))
         endif
         goto 100
      endif

      
      rxyz=rxyz0
      call atomic_copymoving_backward(at,nr,X,n,rxyz)
!      txyz=rxyz
!      alpha=0._gp
!      call atomic_axpy(at,txyz,alpha,sxyz,rxyz)
      in%inputPsiId=1
      in%output_grid=0
      in%output_wf=.false.
!      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,rst,infocode)
      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,fnoise,rst,infocode)
      if(ICALL.ne.0) ncount_bigdft=ncount_bigdft+1
      call atomic_copymoving_forward(at,n,fxyz,nr,G)
      etot=F
      G=-G
      call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
!      call fnrmandforcemax(fxyz,fnrm,fmax,at)

      CALL LBFGS(IPROC,IN,PARMIN,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
!     We allow at most the given number of evaluations of F and G
      if(ncount_bigdft>in%ncount_cluster_x-1)  then
        goto 100
      endif
      close(16)
      open(unit=16,file='geopt.mon',status='unknown',position='append')
      GO TO 20
  50  CONTINUE
        if (iproc==0) write(*,*) "# Error in BFGS, switching to SD and CG"
        if (iproc==0) write(16,*) "Error in BFGS, switching to SD and CG"
        rxyz(:)=rxyzwrite(:)
        fail=.true.
 100  CONTINUE
        
      i_all=-product(shape(X))*kind(X)
      deallocate(X,stat=i_stat)
      call memocc(i_stat,i_all,'X',subname)
      i_all=-product(shape(G))*kind(G)
      deallocate(G,stat=i_stat)
      call memocc(i_stat,i_all,'G',subname)
      i_all=-product(shape(DIAG))*kind(DIAG)
      deallocate(DIAG,stat=i_stat)
      call memocc(i_stat,i_all,'DIAG',subname)
      i_all=-product(shape(W))*kind(W)
      deallocate(W,stat=i_stat)
      call memocc(i_stat,i_all,'W',subname)

END SUBROUTINE lbfgsdriver


subroutine atomic_copymoving_forward(atoms,n,x,nr,xa)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,nr,i,iat,ixyz,ir
    real(8)::x(n),xa(nr)
    logical::move_this_coordinate
    ir=0
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            xa(ir)=x(i)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_forward


subroutine atomic_copymoving_backward(atoms,nr,xa,n,x)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,nr,i,iat,ixyz,ir
    real(8)::x(n),xa(nr)
    logical::move_this_coordinate
    ir=0
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            x(i)=xa(ir)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_backward


!>     This file contains the LBFGS algorithm and supporting routines
!!@todo     LUIGI: PLEASE CUT OUT THIS PART AND PUT IN A TABOO file

!     ----------------------------------------------------------------------
!
!     ****************
!     LBFGS SUBROUTINE
!     ****************
!
!      SUBROUTINE LBFGS(IPROC,IN,N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
!      SUBROUTINE LBFGS(IPROC,IN,PAR,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
      SUBROUTINE LBFGS(IPROC,IN,PARMIN,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
      use module_types
!
!      use par_driver , only:driverparameters 
      use minpar, only: parameterminimization
      IMPLICIT NONE
!      type(driverparameters)::par
      type(parameterminimization)::parmin
      INTEGER N,M,IPRINT(2),IFLAG,IPROC
      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      DOUBLE PRECISION F,EPS
!      DOUBLE PRECISION F,EPS,XTOL
!      LOGICAL DIAGCO
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
! 
!     This subroutine solves the unconstrained minimization problem
! 
!                      min F(x),    x= (x1,x2,...,xN),
!
!      using the limited memory BFGS method. The routine is especially
!      effective on problems involving a large number of variables. In
!      a typical iteration of this method an approximation Hk to the
!      inverse of the Hessian is obtained by applying M BFGS updates to
!      a diagonal matrix Hk0, using information from the previous M steps.
!      The user specifies the number M, which determines the amount of
!      storage required by the routine. The user may also provide the
!      diagonal matrices Hk0 if not satisfied with the default choice.
!      The algorithm is described in "On the limited memory BFGS method
!      for large scale optimization", by D. Liu and J. Nocedal,
!      Mathematical Programming B 45 (1989) 503-528.
! 
!      The user is required to calculate the function value F and its
!      gradient G. In order to allow the user complete control over
!      these computations, reverse  communication is used. The routine
!      must be called repeatedly under the control of the parameter
!      IFLAG. 
!
!      The steplength is determined at each iteration by means of the
!      line search routine MCVSRCH, which is a slight modification of
!      the routine CSRCH written by More' and Thuente.
! 
!      The calling statement is 
! 
!          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
! 
!      where
! 
!     N       is an INTEGER variable that must be set by the user to the
!             number of variables. It is not altered by the routine.
!             Restriction: N>0.
! 
!     M       is an INTEGER variable that must be set by the user to
!             the number of corrections used in the BFGS update. It
!             is not altered by the routine. Values of M less than 3 are
!             not recommended; large values of M will result in excessive
!             computing time. 3<= M <=7 is recommended. Restriction: M>0.
! 
!     X       is a DOUBLE PRECISION array of length N. On initial entry
!             it must be set by the user to the values of the initial
!             estimate of the solution vector. On exit with IFLAG=0, it
!             contains the values of the variables at the best point
!             found (usually a solution).
! 
!     F       is a DOUBLE PRECISION variable. Before initial entry and on
!             a re-entry with IFLAG=1, it must be set by the user to
!             contain the value of the function F at the point X.
! 
!     G       is a DOUBLE PRECISION array of length N. Before initial
!             entry and on a re-entry with IFLAG=1, it must be set by
!             the user to contain the components of the gradient G at
!             the point X.
! 
!     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
!             user  wishes to provide the diagonal matrix Hk0 at each
!             iteration. Otherwise it should be set to .FALSE., in which
!             case  LBFGS will use a default value described below. If
!             DIAGCO is set to .TRUE. the routine will return at each
!             iteration of the algorithm with IFLAG=2, and the diagonal
!              matrix Hk0  must be provided in the array DIAG.
! 
! 
!     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
!             then on initial entry or on re-entry with IFLAG=2, DIAG
!             it must be set by the user to contain the values of the 
!             diagonal matrix Hk0.  Restriction: all elements of DIAG
!             must be positive.
! 
!     IPRINT  is an INTEGER array of length two which must be set by the
!             user.
! 
!             IPRINT(1) specifies the frequency of the output:
!                IPRINT(1) < 0 : no output is generated,
!                IPRINT(1) = 0 : output only at first and last iteration,
!                IPRINT(1) > 0 : output every IPRINT(1) iterations.
! 
!             IPRINT(2) specifies the type of output generated:
!                IPRINT(2) = 0 : iteration count, number of function 
!                                evaluations, function value, norm of the
!                                gradient, and steplength,
!                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
!                                variables and  gradient vector at the
!                                initial point,
!                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
!                                variables,
!                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
! 
! 
!     EPS     is a positive DOUBLE PRECISION variable that must be set by
!             the user, and determines the accuracy with which the solution
!             is to be found. The subroutine terminates when
!
!                         ||G|| < EPS max(1,||X||),
!
!             where ||.|| denotes the Euclidean norm.
! 
!     XTOL    is a  positive DOUBLE PRECISION variable that must be set by
!             the user to an estimate of the machine precision (e.g.
!             10**(-16) on a SUN station 3/60). The line search routine will
!             terminate if the relative width of the interval of uncertainty
!             is less than XTOL.
! 
!     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
!             workspace for LBFGS. This array must not be altered by the
!             user.
! 
!     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
!             to the subroutine. A return with IFLAG<0 indicates an error,
!             and IFLAG=0 indicates that the routine has terminated without
!             detecting errors. On a return with IFLAG=1, the user must
!             evaluate the function F and gradient G. On a return with
!             IFLAG=2, the user must provide the diagonal matrix Hk0.
! 
!             The following negative values of IFLAG, detecting an error,
!             are possible:
! 
!              IFLAG=-1  The line search routine MCSRCH failed. The
!                        parameter INFO provides more detailed information
!                        (see also the documentation of MCSRCH):
!
!                       INFO = 0  IMPROPER INPUT PARAMETERS.
!
!                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
!                                 UNCERTAINTY IS AT MOST XTOL.
!
!                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
!                                 REQUIRED AT THE PRESENT ITERATION.
!
!                       INFO = 4  THE STEP IS TOO SMALL.
!
!                       INFO = 5  THE STEP IS TOO LARGE.
!
!                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
!                                 THERE MAY NOT BE A STEP WHICH SATISFIES
!                                 THE SUFFICIENT DECREASE AND CURVATURE
!                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
!
! 
!              IFLAG=-2  The i-th diagonal element of the diagonal inverse
!                        Hessian approximation, given in DIAG, is not
!                        positive.
!           
!              IFLAG=-3  Improper input parameters for LBFGS (N or M are
!                        not positive).
! 
!
!
!    ON THE DRIVER:
!
!    The program that calls LBFGS must contain the declaration:
!
!                       EXTERNAL LB2
!
!    LB2 is a BLOCK DATA that defines the default values of several
!    parameters described in the COMMON section. 
!
! 
! 
!    COMMON:
! 
!     The subroutine contains one common area, which the user may wish to
!    reference:
! 
!         COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
! 
!    MP  is an INTEGER variable with default value 6. It is used as the
!        unit number for the printing of the monitoring information
!        controlled by IPRINT.
! 
!    LP  is an INTEGER variable with default value 6. It is used as the
!        unit number for the printing of error messages. This printing
!        may be suppressed by setting LP to a non-positive value.
! 
!    GTOL is a DOUBLE PRECISION variable with default value 0.9, which
!        controls the accuracy of the line search routine MCSRCH. If the
!        function and gradient evaluations are inexpensive with respect
!        to the cost of the iteration (which is sometimes the case when
!        solving very large problems) it may be advantageous to set GTOL
!        to a small value. A typical small value is 0.1.  Restriction:
!        GTOL should be greater than 1.D-04.
! 
!    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
!        specify lower and uper bounds for the step in the line search.
!        Their default values are 1.D-20 and 1.D+20, respectively. These
!        values need not be modified unless the exponents are too large
!        for the machine being used, or unless the problem is extremely
!        badly scaled (in which case the exponents should be increased).
! 
!
!  MACHINE DEPENDENCIES
!
!        The only variables that are machine-dependent are XTOL,
!        STPMIN and STPMAX.
! 
!
!  GENERAL INFORMATION
! 
!    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH
! 
!    Input/Output  :  No input; diagnostic messages on unit MP and
!                     error messages on unit LP.
! 
! 
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!      DOUBLE PRECISION GTOL,ONE,ZERO,GNORM,DDOT,STP1,FTOL,STPMIN,STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM
      DOUBLE PRECISION ONE,ZERO,GNORM,DDOT,STP1,STP,YS,YY,SQ,YR,BETA,XNORM
!      INTEGER MP,LP,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
      INTEGER ITER,NFUN,POINT,ISPT,IYPT,INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
      LOGICAL FINISH
      type(input_variables), intent(inout) :: IN 

!
      SAVE
      DATA ONE,ZERO/1.0D+0,0.0D+0/
!
!     INITIALIZE
!     ----------
!
      IF(IFLAG.EQ.0) GO TO 10
      GO TO (172,100) IFLAG
  10  ITER= 0
      IF(N.LE.0.OR.M.LE.0) GO TO 196
      IF(parmin%GTOL.LE.1.D-04) THEN
        IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,245)
        parmin%GTOL=9.D-01
      ENDIF
      NFUN= 1
      POINT= 0
      FINISH= .FALSE.
      IF(parmin%DIAGCO) THEN
         DO 30 I=1,N
 30      IF (DIAG(I).LE.ZERO) GO TO 195
      ELSE
         DO 40 I=1,N
 40      DIAG(I)= 1.0D0
      ENDIF
!
!     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!     ---------------------------------------
!     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!         OTHER TEMPORARY INFORMATION.
!     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!         IN THE FORMULA THAT COMPUTES H*G.
!     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!         STEPS.
!     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!         GRADIENT DIFFERENCES.
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!
      ISPT= N+2*M
      IYPT= ISPT+N*M     
      DO 50 I=1,N
 50   W(ISPT+I)= -G(I)*DIAG(I)
      GNORM= DSQRT(DDOT(N,G,1,G,1))
!      STP1= ONE/GNORM
!      STP11=2.d-2/GNORM  !original convention
      STP1=in%betax
!     PARAMETERS FOR LINE SEARCH ROUTINE
!     They are now set at initialization, see module par_driver
!      FTOL= 1.0D-4
!      FTOL= 1.0D-6  
!      MAXFEV= 20
!      MAXFEV= 10
!
!      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PAR,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PARMIN,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
!
!    --------------------
!     MAIN ITERATION LOOP
!    --------------------
!
 80   ITER= ITER+1
      INFO=0
      BOUND=ITER-1
      IF(ITER.EQ.1) GO TO 165
      IF (ITER .GT. M)BOUND=M
!
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
      IF(.NOT.parmin%DIAGCO) THEN
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         DO 90 I=1,N
   90    DIAG(I)= YS/YY
      ELSE
         IFLAG=2
         RETURN
      ENDIF
 100  CONTINUE
      IF(parmin%DIAGCO) THEN
        DO 110 I=1,N
 110    IF (DIAG(I).LE.ZERO) GO TO 195
      ENDIF
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
      CP= POINT
      IF (POINT.EQ.0) CP=M
      W(N+CP)= ONE/YS
      DO 112 I=1,N
 112  W(I)= -G(I)
      CP= POINT
      DO 125 I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1)CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
 125  CONTINUE
!
      DO 130 I=1,N
 130  W(I)=DIAG(I)*W(I)
!
      DO 145 I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
         BETA= W(N+CP+1)*YR
         INMC=N+M+CP+1
         BETA= W(INMC)-BETA
         ISCN=ISPT+CP*N
         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
 145  CONTINUE
!
!     STORE THE NEW SEARCH DIRECTION
!     ------------------------------
!
       DO 160 I=1,N
 160   W(ISPT+POINT*N+I)= W(I)
!
!     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
!     BY USING THE LINE SEARCH ROUTINE MCSRCH
!     ----------------------------------------------------
 165  NFEV=0
      STP=ONE
      IF (ITER.EQ.1) STP=STP1
      DO 170 I=1,N
 170  W(I)=G(I)
 172  CONTINUE
!      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,XTOL,MAXFEV,INFO,NFEV,DIAG)
!      CALL MCSRCH(PAR,N,X,F,G,W(ISPT+POINT*N+1),STP,INFO,NFEV,DIAG)
      CALL MCSRCH(PARMIN,N,X,F,G,W(ISPT+POINT*N+1),STP,INFO,NFEV,DIAG)
      IF (INFO .EQ. -1) THEN
        IFLAG=1
        RETURN
      ENDIF
      IF (INFO .NE. 1) GO TO 190
      NFUN= NFUN + NFEV
!
!     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
!     -----------------------------------------
!
      NPT=POINT*N
      DO 175 I=1,N
      W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
 175  W(IYPT+NPT+I)= G(I)-W(I)
      POINT=POINT+1
      IF (POINT.EQ.M)POINT=0
!
!     TERMINATION TEST
!     ----------------
!
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      XNORM= DSQRT(DDOT(N,X,1,X,1))
      XNORM= DMAX1(1.0D0,XNORM)
!      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.
      IF (GNORM .LE. EPS) FINISH=.TRUE.
!
!      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PAR,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF(IPRINT(1).GE.0) CALL LB1(IPROC,PARMIN,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF (FINISH) THEN
         IFLAG=0
         RETURN
      ENDIF
      GO TO 80
!
!     ------------------------------------------------------------
!     END OF MAIN ITERATION LOOP. ERROR EXITS.
!     ------------------------------------------------------------
!
 190  IFLAG=-1
      IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,200) INFO
      RETURN
 195  IFLAG=-2
      IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,235) I
      RETURN
 196  IFLAG= -3
      IF(parmin%LP.GT.0.AND.IPROC==0) WRITE(parmin%LP,240)
!
!     FORMATS
!     -------
!
 200  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE',/&
               ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN',/&
               ' OF LINE SEARCH: INFO= ',I2,/&
               ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/&
               ' OR INCORRECT TOLERANCES')
 235  FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,&
            ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
 240  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M',&
            ' ARE NOT POSITIVE)')
 245  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',&
            / ' IT HAS BEEN RESET TO 9.D-01')
      RETURN
      END
!
!     LAST LINE OF SUBROUTINE LBFGS
!
!
!      SUBROUTINE LB1(IPROC,PAR,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      SUBROUTINE LB1(IPROC,PARMIN,IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
!      use par_driver , only:driverparameters
      use minpar, only: parameterminimization
!
!     -------------------------------------------------------------
!     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
!     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
!     -------------------------------------------------------------
      IMPLICIT NONE
!      type(driverparameters)::par
      type(parameterminimization)::parmin
!      INTEGER IPRINT(2),ITER,NFUN,LP,MP,N,M,IPROC
      INTEGER IPRINT(2),ITER,NFUN,N,M,IPROC,I
      DOUBLE PRECISION X(N),G(N),F,GNORM,STP!,GTOL,STPMIN,STPMAX
      LOGICAL FINISH
!      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
!
      IF (ITER.EQ.0)THEN
!      IF(IPROC==0)     WRITE(parmin%MP,10)
      IF(IPROC==0)     WRITE(parmin%MP,20) N,M
      IF(IPROC==0)     WRITE(parmin%MP,30)F,GNORM
                 IF (IPRINT(2).GE.1)THEN
                IF(IPROC==0)     WRITE(parmin%MP,40)
                IF(IPROC==0)     WRITE(parmin%MP,50) (X(I),I=1,N)
                IF(IPROC==0)     WRITE(parmin%MP,60)
                IF(IPROC==0)     WRITE(parmin%MP,50) (G(I),I=1,N)
                  ENDIF
!      IF(IPROC==0)     WRITE(parmin%MP,10)
!      IF(IPROC==0)     WRITE(parmin%MP,70)
      ELSE
          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN
              IF (IPRINT(1).NE.0)THEN
                   IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN
!                         IF(IPRINT(2).GT.1.AND.ITER.GT.1.AND.IPROC==0) WRITE(parmin%MP,70)
!                         IF(IPROC==0)   WRITE(parmin%MP,80)ITER,NFUN,F,GNORM,STP
!                         IF(IPROC==0)   write(parmin%MP,'(i5,1x,e12.5,1x,e21.14,a,i5,1x,e12.5)') ITER,gnorm,f,' GEOPT BFGS ', nfun,STP
                         IF(IPROC==0)   write(parmin%MP,'(a,2(i5,1x),1x,1pe21.14,1x,1pe12.5,1x,1E12.5)')  & 
                                        ' MIN ',ITER,nfun,f,gnorm,STP
                         parmin%FINSTEP=ITER
                         parmin%ALPHA=STP
                   parmin%IWRITE=.true.
                   ELSE
                         RETURN
                   ENDIF
              ELSE
                   IF( IPRINT(2).GT.1.AND.FINISH.AND.IPROC==0) WRITE(parmin%MP,70)
                   IF(IPROC==0)   WRITE(parmin%MP,80)ITER,NFUN,F,GNORM,STP

              ENDIF
              IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN
                    IF (FINISH)THEN
                    IF(IPROC==0)    WRITE(parmin%MP,90)
                    ELSE
                    IF(IPROC==0)    WRITE(parmin%MP,40)
                    ENDIF
                    IF(IPROC==0)  WRITE(parmin%MP,50)(X(I),I=1,N)
                  IF (IPRINT(2).EQ.3)THEN
                    IF(IPROC==0)  WRITE(parmin%MP,60)
                    IF(IPROC==0)  WRITE(parmin%MP,50)(G(I),I=1,N)
                  ENDIF
              ENDIF
            IF (FINISH) WRITE(parmin%MP,100)
      ENDIF
!
 20   FORMAT('  N=',I5,'   NUMBER OF CORRECTIONS=',I2,/,  '       INITIAL VALUES')
 30   FORMAT('  F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40   FORMAT(' VECTOR X= ')
 50   FORMAT(6(2X,1PD10.3))
 60   FORMAT(' GRADIENT VECTOR G= ')
 70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
 80   FORMAT(2(I4,1X),3X,3(1PD10.3,2X))
 90   FORMAT(' FINAL POINT X= ')
 100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.',/' IFLAG = 0')
!
      RETURN
      END
!     ******
!
!
!   ----------------------------------------------------------
!     DATA 
!   ----------------------------------------------------------
!
      BLOCK DATA LB2
      INTEGER LP,MP
      DOUBLE PRECISION GTOL,STPMIN,STPMAX
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
!      DATA MP,LP,GTOL,STPMIN,STPMAX/16,16,9.0D-01,1.0D-20,1.0D+20/
      DATA MP,LP,GTOL,STPMIN,STPMAX/16,16,9.0D-01,1.0D-20,20.d0/
      END
!
!
!!   ----------------------------------------------------------
!!
!      subroutine daxpy(n,da,dx,incx,dy,incy)
!!
!!     constant times a vector plus a vector.
!!     uses unrolled loops for increments equal to one.
!!     jack dongarra, linpack, 3/11/78.
!!
!      double precision dx(1),dy(1),da
!      integer i,incx,incy,ix,iy,m,mp1,n
!!
!      if(n.le.0)return
!      if (da .eq. 0.0d0) return
!      if(incx.eq.1.and.incy.eq.1)go to 20
!!
!!        code for unequal increments or equal increments
!!          not equal to 1
!!
!      ix = 1
!      iy = 1
!      if(incx.lt.0)ix = (-n+1)*incx + 1
!      if(incy.lt.0)iy = (-n+1)*incy + 1
!      do 10 i = 1,n
!        dy(iy) = dy(iy) + da*dx(ix)
!        ix = ix + incx
!        iy = iy + incy
!   10 continue
!      return
!!
!!        code for both increments equal to 1
!!        clean-up loop
!!
!   20 m = mod(n,4)
!      if( m .eq. 0 ) go to 40
!      do 30 i = 1,m
!        dy(i) = dy(i) + da*dx(i)
!   30 continue
!      if( n .lt. 4 ) return
!   40 mp1 = m + 1
!      do 50 i = mp1,n,4
!        dy(i) = dy(i) + da*dx(i)
!        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
!        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
!        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
!   50 continue
!      return
!      end
!
!
!!   ----------------------------------------------------------
!!
!      double precision function ddot(n,dx,incx,dy,incy)
!!
!!     forms the dot product of two vectors.
!!     uses unrolled loops for increments equal to one.
!!     jack dongarra, linpack, 3/11/78.
!!
!      double precision dx(1),dy(1),dtemp
!      integer i,incx,incy,ix,iy,m,mp1,n
!!
!      ddot = 0.0d0
!      dtemp = 0.0d0
!      if(n.le.0)return
!      if(incx.eq.1.and.incy.eq.1)go to 20
!!
!!        code for unequal increments or equal increments
!!          not equal to 1
!!
!      ix = 1
!      iy = 1
!      if(incx.lt.0)ix = (-n+1)*incx + 1
!      if(incy.lt.0)iy = (-n+1)*incy + 1
!      do 10 i = 1,n
!        dtemp = dtemp + dx(ix)*dy(iy)
!        ix = ix + incx
!        iy = iy + incy
!   10 continue
!      ddot = dtemp
!      return
!!
!!        code for both increments equal to 1
!!        clean-up loop
!!
!   20 m = mod(n,5)
!      if( m .eq. 0 ) go to 40
!      do 30 i = 1,m
!        dtemp = dtemp + dx(i)*dy(i)
!   30 continue
!      if( n .lt. 5 ) go to 60
!   40 mp1 = m + 1
!      do 50 i = mp1,n,5
!        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
!   50 continue
!   60 ddot = dtemp
!      return
!      end
!    ------------------------------------------------------------------
!
!     **************************
!     LINE SEARCH ROUTINE MCSRCH
!     **************************
!
!      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA)
!      SUBROUTINE MCSRCH(PAR,N,X,F,G,S,STP,INFO,NFEV,WA)
      SUBROUTINE MCSRCH(PARMIN,N,X,F,G,S,STP,INFO,NFEV,WA)
!      use par_driver , only:driverparameters
      use minpar, only: parameterminimization
      IMPLICIT NONE
!      type(driverparameters)::par
      type(parameterminimization)::parmin
      INTEGER N,INFO,NFEV!,MAXFEV
      DOUBLE PRECISION F,STP!,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
!      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      SAVE
!
!                     SUBROUTINE MCSRCH
!                
!     A slight modification of the subroutine CSRCH of More' and Thuente.
!     The changes are to allow reverse communication, and do not affect
!     the performance of the routine. 
!
!     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
!     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
!
!     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
!     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
!     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
!     MINIMIZER OF THE MODIFIED FUNCTION
!
!          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
!
!     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
!     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
!     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
!     CONTAINS A MINIMIZER OF F(X+STP*S).
!
!     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
!     THE SUFFICIENT DECREASE CONDITION
!
!           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
!
!     AND THE CURVATURE CONDITION
!
!           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
!
!     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
!     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
!     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
!     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
!     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
!     SATISFIES THE SUFFICIENT DECREASE CONDITION.
!
!     THE SUBROUTINE STATEMENT IS
!
!        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES.
!
!       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
!         X + STP*S.
!
!       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
!         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
!
!       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
!         OF F AT X + STP*S.
!
!       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
!         SEARCH DIRECTION.
!
!       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
!         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
!         STP CONTAINS THE FINAL ESTIMATE.
!
!       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
!         communication implementation GTOL is defined in a COMMON
!         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
!         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
!         SATISFIED.
!
!       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
!         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!         IS AT MOST XTOL.
!
!       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
!         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
!         communication implementatin they are defined in a COMMON
!         statement).
!
!       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
!         MAXFEV BY THE END OF AN ITERATION.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!
!         INFO = 0  IMPROPER INPUT PARAMETERS.
!
!         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
!
!         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
!                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
!
!         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!                   IS AT MOST XTOL.
!
!         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
!
!         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
!
!         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
!
!         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
!                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
!                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
!                   TOLERANCES MAY BE TOO SMALL.
!
!       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!         CALLS TO FCN.
!
!       WA IS A WORK ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       MCSTEP
!
!       FORTRAN-SUPPLIED...ABS,MAX,MIN
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!     JORGE J. MORE', DAVID J. THUENTE
!
!     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,&
            FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,&
            STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. parmin%FTOL .LT. ZERO .OR.&
         parmin%GTOL .LT. ZERO .OR. parmin%XTOL .LT. ZERO .OR. parmin%STPMIN .LT. ZERO &
         .OR. parmin%STPMAX .LT. parmin%STPMIN .OR. parmin%MAXFEV .LE. 0) RETURN
!
!     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
!     AND CHECK THAT S IS A DESCENT DIRECTION.
!
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) then
         write(parmin%LP,15)
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
!
!     INITIALIZE LOCAL VARIABLES.
!
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = parmin%FTOL*DGINIT
      WIDTH = parmin%STPMAX - parmin%STPMIN
      WIDTH1 = WIDTH/P5
      DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
!
!     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
!     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
!     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
!     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
!     THE INTERVAL OF UNCERTAINTY.
!     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
!     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
!
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
!
!     START OF ITERATION.
!
   30 CONTINUE
!
!        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
!        TO THE PRESENT INTERVAL OF UNCERTAINTY.
!
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
!
!        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
!
         STP = MAX(STP,parmin%STPMIN)
         STP = MIN(STP,parmin%STPMAX)
!
!        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
!        STP BE THE LOWEST POINT OBTAINED SO FAR.
!
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))&
           .OR. NFEV .GE. parmin%MAXFEV-1 .OR. INFOC .EQ. 0&
           .OR. (BRACKT .AND. STMAX-STMIN .LE. parmin%XTOL*STMAX)) STP = STX
!
!        EVALUATE THE FUNCTION AND GRADIENT AT STP
!        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
!        We return to main program to obtain F and G.
!
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
!
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
!
!        TEST FOR CONVERGENCE.
!
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))&
           .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. parmin%STPMAX .AND.&
            F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. parmin%STPMIN .AND.&
            (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. parmin%MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. parmin%XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. parmin%GTOL*(-DGINIT)) INFO = 1
!
!        CHECK FOR TERMINATION.
!
         IF (INFO .NE. 0) RETURN
!
!        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
!        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
!
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.&
            DG .GE. MIN(parmin%FTOL,parmin%GTOL)*DGINIT) STAGE1 = .FALSE.
!
!        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
!        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
!        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
!        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
!        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
!
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
!
!           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
!
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
!
!           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!           AND TO COMPUTE THE NEW STEP.
!
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,&
                      BRACKT,STMIN,STMAX,INFOC)
!
!           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
!
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
!
!           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!           AND TO COMPUTE THE NEW STEP.
!
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,&
                      BRACKT,STMIN,STMAX,INFOC)
            END IF
!
!        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
!        INTERVAL OF UNCERTAINTY.
!
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)&
              STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
!
!        END OF ITERATION.
!
         GO TO 30
!
!     LAST LINE OF SUBROUTINE MCSRCH.
!
      END
      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX,INFO)
      IMPLICIT NONE
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
!
!     SUBROUTINE MCSTEP
!
!     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
!     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
!     A MINIMIZER OF THE FUNCTION.
!
!     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
!     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
!     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
!     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
!     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
!     WITH ENDPOINTS STX AND STY.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
!                        STPMIN,STPMAX,INFO)
!
!     WHERE
!
!       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
!         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
!         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
!         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
!
!       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
!         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
!         UPDATED APPROPRIATELY.
!
!       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
!         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
!         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
!         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
!
!       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
!         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
!         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
!         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
!
!       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
!         AND UPPER BOUNDS FOR THE STEP.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
!         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
!         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!     JORGE J. MORE', DAVID J. THUENTE
!
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      INFO = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.&
         STP .GE. MAX(STX,STY))) .OR.&
         DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN
!
!     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
!
      SGND = DP*(DX/ABS(DX))
!
!     FIRST CASE. A HIGHER FUNCTION VALUE.
!     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
!     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
!     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
!
      IF (FP .GT. FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .LT. STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
           STPF = STPC + (STPQ - STPC)/2
           END IF
         BRACKT = .TRUE.
!
!     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
!     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
!     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
!     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
!
      ELSE IF (SGND .LT. 0.0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
            END IF
         BRACKT = .TRUE.
!
!     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
!     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
!     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
!     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
!     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
!     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
!     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
!
      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
!
!        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
!        TO INFINITY IN THE DIRECTION OF THE STEP.
!
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
            END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
         ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
            END IF
!
!     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
!     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
!     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
!
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
            END IF
         END IF
!
!     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
!     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
!
      IF (FP .GT. FX) THEN
         STY = STP
         FY = FP
         DY = DP
      ELSE
         IF (SGND .LT. 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
            END IF
         STX = STP
         FX = FP
         DX = DP
         END IF
!
!     COMPUTE THE NEW STEP AND SAFEGUARD IT.
!
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
            END IF
         END IF
      RETURN
!
!     LAST LINE OF SUBROUTINE MCSTEP.
!
      END


!! Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201 (2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University 
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine fire(nproc,iproc,rxyz,at,etot,fxyz,rst,in,ncount_bigdft,fail) 
  use module_base
  use module_types
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

  character(len=*), parameter :: subname='fire'

!Fire parameters:
  real(gp):: alpha,P,finc,fdec,falpha,alphastart,dt,dtmax,vnrm
  real(gp):: velcur(3*at%nat), velpred(3*at%nat),poscur(3*at%nat),pospred(3*at%nat),fcur(3*at%nat),fpred(3*at%nat),mass(3*at%nat)
  real(gp):: ecur,epred,eprev,anoise
  integer:: Nmin,nstep,it

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
  ecur=etot
  epred=etot


  Big_loop: do it=1,in%ncount_cluster_x-1
     do iat=1,3*at%nat
        pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5_gp*fcur(iat)/mass(iat)
     enddo

     in%inputPsiId=1
     in%output_grid=0
     in%output_wf=.false.
     call call_bigdft(nproc,iproc,at,pospred,in,epred,fpred,fnoise,rst,infocode)
     ncount_bigdft=ncount_bigdft+1
     call fnrmandforcemax(fpred,fnrm,fmax,at%nat)
   !  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

     do iat=1,3*at%nat
        velpred(iat)=velcur(iat)+0.5_gp*dt*(fpred(iat))/mass(iat)+0.5_gp*dt*fcur(iat)/mass(iat)
     enddo
     P=dot_product(fpred,velpred)
     call fnrmandforcemax(velpred,vnrm,vmax,at%nat)

     if (iproc == 0) then
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'FIRE:fnrm= ',sqrt(fnrm)
        call  write_atomic_file('posout_'//fn4,epred,pospred,at,trim(comment))
     endif
     if (fmax < 3.d-1) call updatefluctsum(at%nat,fnoise,fluct)
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
     call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
     if (ncount_bigdft >= in%ncount_cluster_x-1) then
         !Too many iterations
         exit Big_loop
     end if
     close(16)
     open(unit=16,file='geopt.mon',status='unknown',position='APPEND')

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

     if (iproc==0) write(10,*) epred, vnrm*0.5d0
   end do Big_loop


        
! Output the final energy, atomic positions and forces
   etot = epred
   rxyz = pospred
   fxyz = fpred

END SUBROUTINE fire



