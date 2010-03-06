!!****m* BigDFT/minimization
!! FUNCTION
!!   Define the type parameterminimization
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
module minpar
  implicit none

  type parameterminimization
     !general parameters for all methods
     character(10)::approach
     real(8)::fmaxtol
!     real(8)::eps=1.d-20
     integer::iter
     integer::iflag
     logical::converged
     !parameters for L-BFGS
     logical::diagco
     !parameters for DIIS
     integer::history
     !parameters for print information
     ! 0: no file output ; 1: normal
     integer :: verbosity
     integer::mp
     integer::lp
     integer::iprint(2)
     !parameters for line search routine
     integer::maxfev
     real(8)::ftol
     real(8)::gtol
     real(8)::stpmin
     real(8)::stpmax
     real(8)::xtol  !epsilon(1.d0)
     real(8)::betax 
  end type parameterminimization

  type(parameterminimization) :: parmin

  ! To be used in lbfgs()
  real(8) ::ys,a_t
  integer ::bound,info,nfun,nfev,point,iypt,ispt
  logical ::finish,new

  ! To be used in mcsrch()
  integer::infoc
  logical::brackt,stage1
  real(8)::dg,dginit,dgtest,dgx,dgy,finit,fx,fy
  real(8)::a_l,a_u,stmin,stmax,width,width1

end module minpar
!!***

!!****f* BigDFT/geopt_init
!! FUNCTION
!!   Geometry optimization, parametrisation routine.
!! SOURCE
!!
subroutine geopt_init()
  use minpar

  parmin%approach  = 'unknown'
  parmin%iter      = 0
  parmin%iflag     = 0
  parmin%verbosity = 1
  parmin%mp        = 6
  parmin%lp        = 6
  parmin%iprint    = (/1,0/)
  parmin%xtol      = 1.d-10
end subroutine geopt_init
!!***

!!****f* BigDFT/geopt_set_verbosity
!! FUNCTION
!!   Geometry optimization, parametrisation routine.
!! SOURCE
!!
subroutine geopt_set_verbosity(verbosity_)
  use minpar

  parmin%verbosity = verbosity_
end subroutine geopt_set_verbosity
!!***

!!****f* BigDFT/geopt
!! FUNCTION
!!   Geometry optimization
!! SOURCE
!!
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

  !-------------------------------------------
  character*4 fn4
  character*40 comment

  if (iproc ==0 .and. parmin%verbosity > 0) &
       & write(16,*) '---------------------------------------------------------------'

  !assign the geometry optimisation method
  parmin%approach=in%geopt_approach

  ncount_bigdft=0
  write(fn4,'(i4.4)') ncount_bigdft
  write(comment,'(a)')'INITIAL CONFIGURATION '
  if (iproc == 0) call write_atomic_file('posout_'//fn4,epot,pos,at,trim(comment))

  if (iproc ==0)  write(*,'(a,1x,a)') ' Begin of minimization using ',parmin%approach
  if(trim(parmin%approach)=='LBFGS') then

     if (iproc ==0) write(*,*) '# ENTERING BFGS'

     call bfgs(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)

     if (fail) then 
        if (iproc ==0) write(*,*) '# ENTERING CG after BFGS failure'
        call conjgrad(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)
     end if

  else if(trim(parmin%approach)=='SDCG') then

     if (iproc ==0) write(*,*) '# ENTERING CG'
     call conjgrad(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)

  else if(trim(parmin%approach)=='VSSD') then
 
     if (iproc ==0) write(*,*) '# ENTERING VSSD'
     call vstepsd(nproc,iproc,pos,at,epot,fxyz,rst,in,ncount_bigdft)

  else if(trim(parmin%approach)=='DIIS') then
 
     if (iproc ==0) write(*,*) '# ENTERING DIIS'
     call rundiis(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)

  else if(trim(parmin%approach)=='AB6MD') then
 
     if (iproc ==0) write(*,*) '# ENTERING Molecular Dynamic (ABINIT implementation)'
     call ab6md(nproc,iproc,pos,fxyz,epot,at,rst,in,ncount_bigdft,fail)

  else
     stop 'geometry optimization method undefined'
  endif
  if (iproc==0)   write(*,'(a,1x,a)') 'End of minimization using ',parmin%approach

END SUBROUTINE geopt
!!***


!!****f* BigDFT/ab6md
!! FUNCTION
!!  Molecular Dynamics
!! SOURCE
!!
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
  acell = (/ at%alat1, at%alat2, at%alat3 /)
  rprim(:,:) = real(0, gp)
  rprim(1,1) = real(1, gp)
  rprim(2,2) = real(1, gp)
  rprim(3,3) = real(1, gp)
  symrel(:,:,1) = real(0, gp)
  symrel(1,1,1) = real(1, gp)
  symrel(2,2,1) = real(1, gp)
  symrel(3,3,1) = real(1, gp)
  do iat = 1, at%nat, 1
     amass(iat) = amu2emass * at%amu(at%iatype(iat))
     do idim = 1, 3, 1
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
end subroutine ab6md
!!***


!!****f* BigDFT/bfgs
!! FUNCTION
!!  Broyden-Fletcher-Goldfarb-Shanno method
!! SOURCE
!!
subroutine bfgs(nproc,iproc,x,f,epot,at,rst,in,ncount_bigdft,fail)
  use module_base
  use module_types
  use minpar
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
  character(len=*), parameter :: subname='bfgs'
  real(gp) :: fluct,fnrm
  real(gp) ::sumx,sumy,sumz,fmax
  logical :: check
  integer :: n,i_all,i_stat,m,i,iter_old
  integer :: infocode,iwrite,nitsd,histlen,nfluct
  real(gp) :: fluctsum,shift,tmp,fnormmax_sw,limbfgs,erel
!  real(gp), dimension(3) :: alat
  real(gp), dimension(30) :: ehist
  real(gp), dimension(:), allocatable :: diag,work,xc,xdft,xwrite
  character*4 fn4
  character*40 comment

  parmin%betax=in%betax
!!!  !Read from input.geopt
!!!  open(84,file="input.geopt")
!!!  read(84,*) parmin%approach
!!!  close(84)
  parmin%approach=in%geopt_approach
  n=3*at%nat
  fail=.false.    
  parmin%converged=.false.
  iwrite=0

  !-------------------------------------------------------------------------------------
  ehist(:)=1.d10
  histlen=6                 !Maximum 30
  limbfgs=1.d-10            !Convergence in 5 consecutive steps in energy before bfgs is stopped
  m=3                       !BFGS history length
  nitsd=20                  !Maximum number of SD steps before entering BFGS
  shift=0.2d0               !Shift for diag of hess
  parmin%iflag=0            !Initialize bfgs
  iter_old=0                !Counter of iterations
  fnormmax_sw=1.d-2         !SD till the max force comp is less than this value
  !parmin%eps=1.d-20         !Original: 1.d-5
  parmin%fmaxtol=1.d-20     !Original: 1.d0
  parmin%ftol=1.d-6         !Original: 1.d-4
  parmin%gtol=9.d-1         !Original: 9.d-1
  parmin%stpmin=1.d-20      !Max size of alpha in line-min, orig value
  parmin%stpmax=1.d+20      !Min size of alpha in line-min, orig value
  parmin%maxfev=10          !Original: 20

     parmin%DIAGCO=.false.
     m=7
     parmin%stpmax=20.d0
  
  fluctsum=0._gp
  nfluct=0

  if (iproc==0)    write(*,*) 'Maximum number of SD steps used in the beginning: ',nitsd

  call steepdes(nproc,iproc,at,x,epot,f,rst,ncount_bigdft,fluctsum,nfluct,fnrm,in,&
       fnormmax_sw,nitsd,fluct)

  call fnrmandforcemax(f,tmp,fmax,at)
  !check if the convergence is reached after SD
  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

  if (check) then
     if (iproc.eq.0) write(*,*) 'Converged before entering BFGS'
     return 
  endif

  allocate(xc(n+ndebug),stat=i_stat)
  call memocc(i_stat,xc,'xc',subname)
  allocate(diag(n+ndebug),stat=i_stat)
  call memocc(i_stat,diag,'diag',subname)
  allocate(xdft(n+ndebug),stat=i_stat)
  call memocc(i_stat,xdft,'xdft',subname)
  allocate(xwrite(n+ndebug),stat=i_stat)
  call memocc(i_stat,xwrite,'xwrite',subname)
  allocate(work(n*(2*m+1)+2*m+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)
  xc(:)=x(:)
  xwrite(:)=x(:)

  fail=.false.
  bfgs_loop: do ! main BFGS loop

     if(ncount_bigdft==0 .or. parmin%iflag==1) then

        in%inputPsiId=1
        call call_bigdft(nproc,iproc,at,x,in,epot,f,rst,infocode)

        if (iproc == 0) then
           call transforce(at,f,sumx,sumy,sumz)
           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
        end if

        ncount_bigdft=ncount_bigdft+1
        if (ncount_bigdft==1) ehist(1)=epot
        xdft(:)=x(:)

        call fnrmandforcemax(f,fnrm,fmax,at)
        if (fmax < 3.d-1) call updatefluctsum(at,f,nfluct,fluctsum,fluct)

        call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

     endif

     call lbfgs(at,3*at%nat,m,x,xc,epot,f,diag,work,parmin,iproc,iwrite)

     if (iwrite==iter_old+1 ) then
        xwrite(:)=xdft(:)
        if (iproc==0) then 
           write(fn4,'(i4.4)') ncount_bigdft
           write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
           call write_atomic_file('posout_'//fn4,epot,xwrite,at,trim(comment))
        endif

        ehist(2:histlen)=ehist(1:histlen-1)
        ehist(1)=epot
        erel=0.d0
        do i=1,histlen-1
           erel=erel+(ehist(i+1)-ehist(i))/abs(ehist(i))
        enddo
        erel=erel/real(histlen-1,8)
        if (iproc==0 .and. parmin%verbosity > 0) &
             & write(16,*) "BFGS: erel energy convergence", erel 
        if (erel.lt.limbfgs) then
           if(iproc==0 .and. parmin%verbosity > 0) &
             & write(16,*) "BFGS: No progress in BFGS, switching to SD and CG", ncount_bigdft,sum(ehist)/5.d0.lt.1.d-7
           if(iproc==0) write(*,*) "# BFGS: No progress in BFGS, switching to SD and CG", ncount_bigdft,sum(ehist)/5.d0.lt.1.d-7
           x(:)=xwrite(:)
           fail=.true.
           exit
        endif


        ! If user realizes the BFGS does not work efficiently, a switch to CG can be obtained by putting SDCG in the input.geopt file
        open(84,file="input.geopt")
        read(84,*) parmin%approach
        close(84)
        if (trim(parmin%approach)=='SDCG') then
           if(iproc==0 .and. parmin%verbosity > 0) &
             & write(16,*) "BFGS: Manual switchback to SD and CG", ncount_bigdft 
           if(iproc==0) write(*,*) "# BFGS: Manual switchback to SD and CG", ncount_bigdft
           x(:)=xwrite(:)
           fail=.true.
           exit
        endif

        if (iproc==0 .and. parmin%verbosity > 0) &
             & write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*in%frac_fluct,fluct
        iter_old=iwrite
     endif

     if (iproc.eq.0 .AND. parmin%iflag /= 2)  write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=', &
          fmax,'fnrm=',    fnrm    ,'fluct=', fluct 

     if(check)then
        if (iproc==0 .and. parmin%verbosity > 0) then
           write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*in%frac_fluct,fluct
           write(16,*) 'BFGS converged'
        end if
        exit
     endif
     if(parmin%iflag<=0) then
        if (iproc==0) write(*,*) "# Error in BFGS, switching to SD and CG"
        if (iproc==0 .and. parmin%verbosity > 0) &
             & write(16,*) "Error in BFGS, switching to SD and CG"
        x(:)=xwrite(:)
        fail=.true. 
        !write(100+iproc,*) 'positions:',x
        exit
     endif

     if(ncount_bigdft>in%ncount_cluster_x-1)  then 
      if (iproc==0 .and. parmin%verbosity > 0) &
           & write(16,*) 'BFGS exited before the geometry optimization converged because more than ',& 
           in%ncount_cluster_x,' wavefunction optimizations were required'
      exit
     endif


  enddo bfgs_loop   ! main BFGS loop


  i_all=-product(shape(xc))*kind(xc)
  deallocate(xc,stat=i_stat)
  call memocc(i_stat,i_all,'xc',subname)

  i_all=-product(shape(diag))*kind(diag)
  deallocate(diag,stat=i_stat)
  call memocc(i_stat,i_all,'diag',subname)

  i_all=-product(shape(xdft))*kind(xdft)
  deallocate(xdft,stat=i_stat)
  call memocc(i_stat,i_all,'xdft',subname)

  i_all=-product(shape(xwrite))*kind(xwrite)
  deallocate(xwrite,stat=i_stat)
  call memocc(i_stat,i_all,'xwrite',subname)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

  !-------------------------------------------------------------------------------------
END SUBROUTINE bfgs
!!***


subroutine timeleft(tt)
  use module_base
  implicit none
  real(gp), intent(out) :: tt
  !local variables
  ! MODIFIED version for refined time limit on restart of global.f90.
  ! Only difference: Calls routine CPUtime(tt)
  integer :: ierr
  real(kind=4) :: tcpu
  real(gp) :: timelimit

  open(unit=55,file='CPUlimit',status='unknown')
  read(55,*,iostat=ierr) timelimit ! in hours
  if(ierr/=0)timelimit=1d6
  close(55)
  call cpu_time(tcpu)
  tt=timelimit-real(tcpu,gp)/3600._gp ! in hours
end subroutine timeleft


!!****f* BigDFT/conjgrad
!! FUNCTION
!!  Conjugate gradient method
!! SOURCE
!!
subroutine conjgrad(nproc,iproc,rxyz,at,etot,fxyz,rst,in,ncount_bigdft)
  use module_base
  use module_types
  use minpar
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  real(gp), intent(out) :: etot
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: fxyz
  !local variables
  real(gp) ::fluctsum 
  integer :: nfluct
  character(len=*), parameter :: subname='conjgrad'  
  integer :: nfail,it,iat,i_all,i_stat,infocode, nitsd
  real(gp) :: anoise,fluct,avbeta,avnum,fnrm,etotprec,beta0,beta
  real(gp) :: y0,y1,tt,sumx,sumy,sumz,oben1,oben2,oben,unten,rlambda,tetot,fmax,tmp
  real(gp), dimension(:,:), allocatable :: tpos,gpf,hh
  logical::check
  character*4 fn4
  character*40 comment

  allocate(tpos(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,tpos,'tpos',subname)
  allocate(gpf(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gpf,'gpf',subname)
  allocate(hh(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,hh,'hh',subname)

  anoise=1.e-4_gp
  fluctsum=0._gp 
  nfluct=0

  avbeta=0._gp
  avnum=0._gp
  nfail=0
  nitsd=500

  !start with a steepest descent algorithm
  call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,&
       fluctsum,nfluct,fnrm,in,in%forcemax,nitsd,fluct)
  if (ncount_bigdft >= in%ncount_cluster_x) then
      if (iproc==0 .and. parmin%verbosity > 0) &
           & write(16,*) 'SDCG exited before the geometry optimization converged because more than ',&
           in%ncount_cluster_x,' wavefunction optimizations were required'
     call close_and_deallocate
     return
  end if

  !calculate the max of the forces
  call fnrmandforcemax(fxyz,tmp,fmax,at)

  !control whether the convergence criterion is reached after SD
  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
  if (check) then
     if (iproc.eq.0) write(16,*) 'Converged before entering CG',iproc
     call close_and_deallocate
     return
  endif

  redo_cg: do
     etotprec=etot
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
        call atomic_axpy(at,rxyz,beta0,hh,tpos)

        in%inputPsiId=1
        in%output_grid=0
        in%output_wf=.false.
        call call_bigdft(nproc,iproc,at,tpos,in,tetot,gpf,rst,infocode)
        if (iproc == 0) then
           call transforce(at,gpf,sumx,sumy,sumz)
           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
        end if
        ncount_bigdft=ncount_bigdft+1

        !C projection of gradients at beta=0 and beta onto hh
        y0=0._gp
        y1=0._gp
        do iat=1,at%nat
           y0=y0+fxyz(1,iat)*hh(1,iat)+fxyz(2,iat)*hh(2,iat)+fxyz(3,iat)*hh(3,iat)
           y1=y1+gpf(1,iat)*hh(1,iat)+gpf(2,iat)*hh(2,iat)+gpf(3,iat)*hh(3,iat)
        end do
        tt=y0/(y0-y1)
        if (iproc.eq.0) then
           if (parmin%verbosity > 0) &
                & write(16,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
           write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
        end if

        beta=beta0*max(min(tt,2._gp),-.25_gp)
        
        tpos=rxyz
        call atomic_axpy(at,rxyz,beta,hh,rxyz)
!!!        do iat=1,at%nat
!!!           tpos(1,iat)=rxyz(1,iat)
!!!           tpos(2,iat)=rxyz(2,iat)
!!!           tpos(3,iat)=rxyz(3,iat)
!!!           if ( .not. at%lfrztyp(iat)) then
!!!              if (at%geocode == 'P') then
!!!                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*hh(1,iat),at%alat1)
!!!                 rxyz(2,iat)=modulo(rxyz(2,iat)+beta*hh(2,iat),at%alat2)
!!!                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*hh(3,iat),at%alat3)
!!!              else if (at%geocode == 'S') then
!!!                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*hh(1,iat),at%alat1)
!!!                 rxyz(2,iat)=rxyz(2,iat)+beta*hh(2,iat)
!!!                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*hh(3,iat),at%alat3)
!!!              else
!!!                 rxyz(1,iat)=rxyz(1,iat)+beta*hh(1,iat)
!!!                 rxyz(2,iat)=rxyz(2,iat)+beta*hh(2,iat)
!!!                 rxyz(3,iat)=rxyz(3,iat)+beta*hh(3,iat)
!!!              end if
!!!
!!!           end if
!!!        end do
        avbeta=avbeta+beta/in%betax
        avnum=avnum+1._gp

        !C new gradient
        do iat=1,at%nat
           gpf(1,iat)=fxyz(1,iat)
           gpf(2,iat)=fxyz(2,iat)
           gpf(3,iat)=fxyz(3,iat)
        end do

        call call_bigdft(nproc,iproc,at,rxyz,in,etot,fxyz,rst,infocode)
        if (iproc == 0) then
           call transforce(at,fxyz,sumx,sumy,sumz)
           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
        end if
        ncount_bigdft=ncount_bigdft+1
        !if the energy goes up (a small tolerance of anoise is allowed)
        !switch back to SD
        if (etot > etotprec+anoise) then

           if (iproc.eq.0 .and. parmin%verbosity > 0) &
                & write(16,*) 'switching back to SD:etot,etotprec',it,etot,etotprec
           if (iproc.eq.0) write(*,*) ' switching back to SD:etot,etotprec',it,etot,etotprec
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,&
                fluctsum,nfluct,fnrm,in,in%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(fxyz,tmp,fmax,at)
           call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
           if(check) then
              if (iproc.eq.0) write(16,*) 'Converged in switch back SD',iproc
              call close_and_deallocate
              return
           endif
           cycle redo_cg
        endif

        etotprec=etot
        if (iproc==0) then 
           write(fn4,'(i4.4)') ncount_bigdft
           write(comment,'(a,1pe10.3)')'CONJG:fnrm= ',sqrt(fnrm)
           call  write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))
        endif


        !if (iproc.eq.0) write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'CG ',ncount_bigdft,etot,sqrt(fnrm)

        if (fmax < 3.d-1) call updatefluctsum(at,fxyz,nfluct,fluctsum,fluct)

        call atomic_dot(at,gpf,gpf,unten)
        call atomic_dot(at,gpf,fxyz,oben1)
        call atomic_dot(at,fxyz,fxyz,oben2)
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

        call fnrmandforcemax(fxyz,fnrm,fmax,at)
        if (iproc.eq.0) then
           if (parmin%verbosity > 0) then
              write(16,'(i5,1x,e12.5,1x,e21.14,a,1x,e9.2)')it,sqrt(fnrm),etot,' GEOPT CG ',beta/in%betax
              write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*in%frac_fluct,fluct
           end if
           write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                'FORCES norm(Ha/Bohr): maxval=',    fmax,'fnrm=',    fnrm   , 'fluct=',fluct
        end if

        call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
        if(check) exit loop_cg

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
        if (it.eq.500) then
           if (iproc.eq.0) then
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
                fluctsum,nfluct,fnrm,in,in%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(fxyz,tmp,fmax,at)

           call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)
           if(check) then
              if (iproc.eq.0) write(16,*) 'Converged in back up SD',iproc
              call close_and_deallocate
              return
           endif

           nfail=nfail+1
           if (nfail.ge.100) stop 'too many failures of CONJG'
           cycle redo_cg
        endif

        !rlambda=(obenx+obeny+obenz)/unten
        rlambda=oben/unten
        call atomic_axpy_forces(at,fxyz,rlambda,hh,hh)
!!!        do iat=1,at%nat
!!!           hh(1,iat)=fxyz(1,iat)+rlambda*hh(1,iat)
!!!           hh(2,iat)=fxyz(2,iat)+rlambda*hh(2,iat)
!!!           hh(3,iat)=fxyz(3,iat)+rlambda*hh(3,iat)
!!!        end do
     end do loop_cg
     exit redo_cg
  end do redo_cg

  !!        write(6,*) 'CG finished',it,fnrm,etot
  if (iproc == 0) then
     if (parmin%verbosity > 0) &
          & write(16,*) 'average CG stepsize in terms of betax',avbeta/avnum,iproc
     write(*,*) 'average CG stepsize in terms of betax',avbeta/avnum,iproc
  end if

 call close_and_deallocate

contains

  subroutine close_and_deallocate
    use module_base
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
  end subroutine close_and_deallocate

end subroutine conjgrad
!!***


!!****f* BigDFT/steepdes
!! FUNCTION
!!  Steepedt descent method
!! SOURCE
!!
subroutine steepdes(nproc,iproc,at,rxyz,etot,ff,rst,ncount_bigdft,fluctsum,&
     nfluct,fnrm,in,forcemax_sw,nitsd,fluct)
  use module_base
  use module_types
  use minpar
  !use module_interfaces
  implicit none
  integer, intent(in) :: nproc,iproc,nitsd,nfluct
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  integer, intent(inout) :: ncount_bigdft
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  real(gp), intent(out) :: fnrm,etot
  real(gp), intent(inout) :: fluctsum,fluct
  real(gp), dimension(3,at%nat), intent(out) ::ff
  real(gp), intent(in)::forcemax_sw
  !local variables
  character(len=*), parameter :: subname='steepdes'
  logical :: care
  integer :: nsatur,iat,itot,itsd,i_stat,i_all,infocode,nbeqbx
  real(gp) :: etotitm2,fnrmitm2,etotitm1,fnrmitm1,anoise,sumx,sumy,sumz
  real(gp) :: fmax,de1,de2,df1,df2,beta
  real(gp), allocatable, dimension(:,:) :: tpos
  character*4 fn4
  character*40 comment

  allocate(tpos(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,tpos,'tpos',subname)

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

  redo_sd: do
     if (ncount_bigdft.gt.in%ncount_cluster_x) then 
        if (iproc.eq.0) then
           write(*,*) 'SD FINISHED because ncount_bigdft > ncount_cluster_x',iproc,ncount_bigdft,in%ncount_cluster_x
           if (parmin%verbosity > 0) then
              write(16,*) 'SD FINISHED because ncount_bigdft > ncount_cluster_x',iproc,ncount_bigdft,in%ncount_cluster_x
              write(16,*) 'SDCG exited before the geometry optimization converged because more than ', & 
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
        call call_bigdft(nproc,iproc,at,rxyz,in,etot,ff,rst,infocode)
        if (iproc == 0) then
           call transforce(at,ff,sumx,sumy,sumz)
           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
        end if
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

        call fnrmandforcemax(ff,fnrm,fmax,at)

        !first and second derivatives of the energy and of the norm of the forces
        !(in units of beta steps)
        de1=etot-etotitm1
        de2=etot-2._gp*etotitm1+etotitm2
        df1=fnrm-fnrmitm1
        df2=fnrm-2._gp*fnrmitm1+fnrmitm2

        if (fmax < 3.d-1) call updatefluctsum(at,ff,nfluct,fluctsum,fluct)
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
             fmax <= .03_gp .and. fnrm/fnrmitm1 > .50_gp .and. &   !norm of forces is saturating
             de2 > -2._gp*anoise .and. df2 > -2._gp*anoise) then !close to a local minimum (E&F)
           nsatur=nsatur+1
        else
           nsatur=0
        endif

        if (iproc.eq.0) then 
           if (parmin%verbosity > 0) &
                & write(16,'(i5,1x,e12.5,1x,e21.14,a)') itsd,sqrt(fnrm),etot,' GEOPT SD '
           write(fn4,'(i4.4)') ncount_bigdft 
           write(comment,'(a,1pe10.3)')'SD:fnrm= ',sqrt(fnrm)
           call write_atomic_file('posout_'//fn4,etot,rxyz,at,trim(comment))

           !write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'SD ',ncount_bigdft,etot,sqrt(fnrm)
        end if


        if (iproc==0 .and. parmin%verbosity > 0) then
           write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',&
                fnrm,fluct*in%frac_fluct,fluct
        end if

        !exit statements

        !the norm of the forces is of the same order of the fluctuations
        !or we are much too close to the fixed point from the SD viewpoint
        !or the forces are below the fixed tolerance

        if (iproc==0) then
            if (fnrm < fluct*in%frac_fluct ) write(*,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fnrm < fluct*frac_fluct' ,&
                 fnrm,fluct*in%frac_fluct
            if (fmax < forcemax_sw ) write(*,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < forcemax_sw ',&
                 fmax , forcemax_sw
            if (nsatur > 5  ) write(*,'(a,1x,i0)')&
                 'SD EXIT because nsatur > 5 ',nsatur 
            if (fnrm < fluct*in%frac_fluct ) write(16,'(a,2(1x,1pe24.17))') &
                 'SD EXIT because fnrm < fluct*frac_fluct',&
                 fnrm , fluct*in%frac_fluct
            if (fmax < forcemax_sw ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < forcemax_sw ' ,fmax , forcemax_sw
            if (nsatur > 5  ) write(16,'(a,1x,i0)') 'SD EXIT because nsatur > 5 ',nsatur 
        endif

        if (fnrm < fluct*in%frac_fluct .or. &
             nsatur > 5 .or. fmax < forcemax_sw) &
             exit loop_sd

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
                'SD: NO CONVERGENCE:itsd,fnrm,etot',itsd,fnrm,etot
           exit loop_sd
        endif
        if (itot >= nitsd) then
           if (iproc.eq.0) write(16,'(a,i5,i5,1x,e10.3,1x,e21.14)') &
                'SD: NO CONVERGENCE:itsd,itot,fnrm,etot:',itsd,itot,fnrm,etot
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
        if (iproc.eq.0 .and. parmin%verbosity > 0) write(16,*) 'beta=',beta

        tpos=rxyz
        call atomic_axpy(at,rxyz,beta,ff,rxyz)

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

  if (iproc.eq.0 .and. parmin%verbosity > 0) write(16,*) 'SD FINISHED',iproc

  i_all=-product(shape(tpos))*kind(tpos)
  deallocate(tpos,stat=i_stat)
  call memocc(i_stat,i_all,'tpos',subname)

END SUBROUTINE steepdes
!!***


subroutine vstepsd(nproc,iproc,wpos,at,etot,ff,rst,in,ncount_bigdft)
! variable step steepest descent
  use module_base
  use module_types
  use minpar
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  real(gp), intent(out) :: etot
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), dimension(3,at%nat), intent(inout) :: wpos
  real(gp), dimension(3,at%nat), intent(out) :: ff
  !local variables
  real(gp) ::fluctsum 
  integer :: nfluct
  character(len=*), parameter :: subname='vstepsd'  
  integer :: iat,i_all,i_stat,infocode, nitsd,itsd
  real(gp) :: anoise,fluct,fnrm,fnrmold,beta,betaxx,betalast,betalastold
  real(gp) :: etotold,fmax,scpr,curv,tt,sumx,sumy,sumz
  real(gp), dimension(:,:), allocatable :: posold,ffold
  logical reset,check
  character*4 fn4
  character*40 comment


  allocate(posold(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,posold,'posold',subname)
  allocate(ffold(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,ffold,'ffold',subname)

  anoise=1.e-4_gp
  fluctsum=0._gp 
  fluct=0._gp
  nfluct=0

  beta=in%betax

  in%inputPsiId=1
  in%output_grid=0
  in%output_wf=.false.
  call call_bigdft(nproc,iproc,at,wpos,in,etotold,ffold,rst,infocode)
  call fnrmandforcemax(ffold,fnrm,fmax,at)   
  if (fmax < 3.d-1) call updatefluctsum(at,ffold,nfluct,fluctsum,fluct)
  if (iproc == 0) then
     if (parmin%verbosity > 0) &
          & write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',fnrm,fluct*in%frac_fluct,fluct
     call transforce(at,ffold,sumx,sumy,sumz)                         
     write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
     write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
     write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  

     write(fn4,'(i4.4)') ncount_bigdft
     write(comment,'(a,1pe10.3)')'Initial VSSD:fnrm= ',sqrt(fnrm)
     call  write_atomic_file('posout_'//fn4,etotold,wpos,at,trim(comment))
     if (parmin%verbosity > 0) &
          & write(16,'(1x,e12.5,1x,e21.14,a,e10.3)')sqrt(fnrm),etotold,' GEOPT VSSD ',beta
  end if

  ncount_bigdft=ncount_bigdft+1


  fnrmold=0.d0
  do iat=1,at%nat
     !          fnrmold=fnrmold+ffold(1,iat)**2+ffold(2,iat)**2+ffold(3,iat)**2
     posold(1,iat)=wpos(1,iat)
     posold(2,iat)=wpos(2,iat)
     posold(3,iat)=wpos(3,iat)
     !          wpos(1,iat)=wpos(1,iat)+beta*ffold(1,iat)
     !          wpos(2,iat)=wpos(2,iat)+beta*ffold(2,iat)
     !          wpos(3,iat)=wpos(3,iat)+beta*ffold(3,iat)
  enddo
  call atomic_dot(at,ffold,ffold,fnrmold)
  call atomic_axpy(at,wpos,beta,ffold,wpos)
  betaxx=1.d100
  reset=.true.
  betalastold=in%betax


  nitsd=1000
  do itsd=1,nitsd
     in%inputPsiId=1
     in%output_grid=0
     in%output_wf=.false.
     call call_bigdft(nproc,iproc,at,wpos,in,etot,ff,rst,infocode)
     if (iproc == 0) then                                        
        call transforce(at,ff,sumx,sumy,sumz)                         
        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
     end if
     ncount_bigdft=ncount_bigdft+1
     fnrm=0.d0
     scpr=0.d0
     !          do iat=1,at%nat
     !          fnrm=fnrm+ff(1,iat)**2+ff(2,iat)**2+ff(3,iat)**2
     !          scpr=scpr+ffold(1,iat)*ff(1,iat)+ffold(2,iat)*ff(2,iat)+ffold(3,iat)*ff(3,iat)
     !          enddo
     call  atomic_dot(at,ff,ff,fnrm)
     call  atomic_dot(at,ff,ffold,scpr)
     curv=(fnrmold-scpr)/(beta*fnrmold)
     betalast=.5d0/curv
     if (reset) betaxx=min(betaxx,1.5d0*betalast)
     call fnrmandforcemax(ff,fnrm,fmax,at)   
     if (fmax < 3.d-1) call updatefluctsum(at,ff,nfluct,fluctsum,fluct)
     if (iproc==0) write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',fnrm,fluct*in%frac_fluct,fluct
     call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
     if (check) goto 100
     if (ncount_bigdft >= in%ncount_cluster_x) then 
        if (iproc==0)  write(16,*) 'VSSD exited before the geometry optimization converged because more than ',& 
             in%ncount_cluster_x,' wavefunction optimizations were required'
        goto 100
     endif


     if (etot.gt.etotold) then
        reset=.true.
        beta=in%betax
        if (iproc == 0) write(16,*) 'new positions rejected, reduced beta',beta
        !          do iat=1,at%nat
        !          wpos(1,iat)=posold(1,iat)+beta*ffold(1,iat)
        !          wpos(2,iat)=posold(2,iat)+beta*ffold(2,iat)
        !          wpos(3,iat)=posold(3,iat)+beta*ffold(3,iat)
        !          enddo
        call atomic_axpy(at,posold,beta,ffold,wpos)
     else
        reset=.false.
        if (betalast.gt.0) then
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
           !          wpos(1,iat)=wpos(1,iat)+beta*ff(1,iat)
           !          wpos(2,iat)=wpos(2,iat)+beta*ff(2,iat)
           !          wpos(3,iat)=wpos(3,iat)+beta*ff(3,iat)
           ffold(1,iat)=ff(1,iat)
           ffold(2,iat)=ff(2,iat)
           ffold(3,iat)=ff(3,iat)
        enddo
        call atomic_axpy(at,wpos,beta,ff,wpos)
        etotold=etot
        fnrmold=fnrm
        betalastold=betalast
     endif

     if (iproc == 0) write(16,'(i5,1x,e12.5,1x,e21.14,a,e10.3,1x,e10.3)') itsd,sqrt(fnrm),etot,' GEOPT VSSD ',beta,betalast
     if(iproc==0)call timeleft(tt)
     call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_stat)
     if(tt<0) goto 100


  enddo
  if (iproc == 0) write(16,'(a,i5,e9.2,e18.10,e9.2)') '---- SD FAILED  TO CONVERGE'
100 if (iproc == 0) then
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


subroutine convcheck(fnrm, fmax, fluctfrac_fluct, forcemax, check)
  use module_base
  implicit none
  real(gp), intent(in):: fnrm, fmax, fluctfrac_fluct,forcemax
  logical, intent(out)::check


  check=.false.
  if (fnrm < fluctfrac_fluct .or. &
       fmax < forcemax ) then
     check=.true.
  endif

END SUBROUTINE convcheck

subroutine fnrmandforcemax(ff,fnrm,fmax,at)
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
END SUBROUTINE fnrmandforcemax


subroutine updatefluctsum(at,fxyz,nfluct,fluctsum,fluct)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp),intent(inout):: fluctsum,fluct
  real(gp), dimension(3,at%nat), intent(in) :: fxyz
  integer, intent(inout):: nfluct
  !local variables
!!!  real(gp), save :: fluct_im1,fluct_im2,fluct_i
  real(gp) :: sumx,sumy,sumz


  call transforce_forfluct(at,fxyz,sumx,sumy,sumz)
  nfluct=nfluct+1

!!!  !limit the fluctuation value to n=3
!!!  !to be done only in the case of SDCG
!!!  !for BFGS it is best to consider everything
!!!  if (nfluct >= 3) then
!!!     fluct_im2=fluct_im1
!!!     fluct_im1=fluct_i
!!!  else if (nfluct == 2) then
!!!     fluct_im2=0.0_gp
!!!     fluct_im1=fluct_i
!!!  else if (nfluct == 1) then
!!!     fluct_im2=0.0_gp
!!!     fluct_im1=0.0_gp
!!!  end if
!!!  fluct_i=sumx**2+sumy**2+sumz**2
!!!  fluct=sqrt(real(nat,gp))*(fluct_i+fluct_im1+fluct_im2)/3.0_gp

  fluctsum=fluctsum+sumx**2+sumy**2+sumz**2
  !commented out, it increases the fluctuation artificially
  fluct=fluctsum*sqrt(real(at%nat,gp))/real(nfluct,gp)
END SUBROUTINE updatefluctsum

!should we evaluate the translational force also with blocked atoms?
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

!should we evaluate the translational force also with blocked atoms?
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



!*****************************************************************************************
subroutine lbfgs(at,n,m,x,xc,f,g,diag,w,parmin,iproc,iwrite)
  use module_base
  use module_types
  use minpar, only: parameterminimization,ys,a_t,bound, &
       & info,nfun,nfev,point,iypt,ispt,finish,new
  implicit none
  integer :: n,m,iproc,iwrite
  type(atoms_data), intent(in) :: at
  type(parameterminimization), intent(inout) :: parmin
  real(8)::x(n),xc(n),g(n),diag(n),w(n*(2*m+1)+2*m),f
  real(8)::one,zero,gnorm,stp1,xnorm,beta,yr,sq,yy
  integer::npt,cp,i,inmc,iycn,iscn
  data one,zero/1.0d+0,0.0d+0/

  if(parmin%iflag==0) then
     call init_lbfgs(at,n,m,g,diag,w,parmin,nfun,point,finish,stp1,ispt,iypt)
  endif
  new=.false.
  if(parmin%iflag==0) new=.true.
  !main iteration loop
  do 
     if(new) then
        parmin%iter=parmin%iter+1
        info=0
        bound=parmin%iter-1
        if(parmin%iter/=1) then
           if(parmin%iter>m) bound=m
           call atomic_dot(at,w(iypt+npt+1),w(ispt+npt+1),ys)
!!!           ys=ddot(n,w(iypt+npt+1),1,w(ispt+npt+1),1)
           if(.not.parmin%diagco) then
              call atomic_dot(at,w(iypt+npt+1),w(iypt+npt+1),yy)
!!!              yy=ddot(n,w(iypt+npt+1),1,w(iypt+npt+1),1)
              do i=1,n
                 diag(i)= ys/yy
              enddo
           else
              parmin%iflag=2
              return
           endif
        endif
     endif
     if((new .or. parmin%iflag==2) .and. parmin%iter/=1) then
        !write(*,*) 'ALI',iflag,iter,new
        if(parmin%diagco) then
           do i=1,n
              if(diag(i)<=zero) then
                 parmin%iflag=-2
                 write(*,'(a,i5,2a)') 'iflag=-2, the',i,'-th diagonal element ', &
                      'of the inverse hessian approximation is not positive'
                 return
              endif
           enddo
        endif
        !compute -h*g using the formula given in: nocedal, j. 1980,
        !"updating quasi-newton matrices with limited storage",
        !mathematics of computation, vol.24, no.151, pp. 773-782.
        cp= point
        if(point==0) cp=m
        w(n+cp)= one/ys
        do i=1,n
           w(i)=g(i)
        enddo
        cp=point
        do i=1,bound
           cp=cp-1
           if(cp==-1)cp=m-1
           call atomic_dot(at,w(ispt+cp*n+1),w,sq)
!!!           sq= ddot(n,w(ispt+cp*n+1),1,w,1)
           inmc=n+m+cp+1
           iycn=iypt+cp*n
           w(inmc)= w(n+cp+1)*sq
           call atomic_axpy(at,w,-w(inmc),w(iycn+1),w)
!!!           call daxpy(n,-w(inmc),w(iycn+1),1,w,1)
        enddo
        !to be used only in the case parmin%diagco==.false.
        !in that case diag=constant
        if (.not. parmin%diagco) then
           call atomic_axpy(at,w,diag(1)-1.0_gp,w,w)
        else
           !here the blocked atoms are treated as the others
           !also atomic_gemv can be used here
           do i=1,n
              w(i)=diag(i)*w(i)
           enddo
        end if


        do i=1,bound
           call atomic_dot(at,w(iypt+cp*n+1),w,yr)
!!!           yr=ddot(n,w(iypt+cp*n+1),1,w,1)
           beta= w(n+cp+1)*yr
           inmc=n+m+cp+1
           beta= w(inmc)-beta
           iscn=ispt+cp*n
           call atomic_axpy(at,w,beta,w(iscn+1),w)
!!!           call daxpy(n,beta,w(iscn+1),1,w,1)
           cp=cp+1
           if(cp==m) cp=0
        enddo
        do i=1,n
           w(ispt+point*n+i)=w(i) !store the new search direction
        enddo
        !obtain the one-dimensional minimizer of the function by using 
        !the line search routine mcsrch
     endif
     if(parmin%iflag/=1 .or. new) then
        nfev=0
        a_t=one
        if(parmin%iter==1) a_t=stp1
        do i=1,n
           w(i)=-g(i)
        enddo
     endif
     !write(*,*) 'a_t',a_t
     call mcsrch(at,n,x,f,g,w(ispt+point*n+1),a_t,info,nfev,diag,parmin)
     if(iproc==0) write(*,'(a,1pe12.5,i6)') "ALPHA LINESEARCH ", a_t, parmin%iter
     if (info==-1) then
        parmin%iflag=1
        return
     endif
     if(info/=1) then
        parmin%iflag=-1
        write(*,'(2a)') 'iflag=-1, line search failed. ', & 
             ' see documentation of routine mcsrch.'
        write(*,'(a,i2)') 'error return of line search: info= ',info
        write(*,'(2a)') 'possible causes: ', &
             ' function or gradient are incorrect or incorrect tolerances.'
        return
     endif
     nfun=nfun+nfev
     !compute the new step and gradient change
     npt=point*n
     call atomic_axpy(at,w(ispt+npt+1),a_t-1.0_gp,w(ispt+npt+1),w(ispt+npt+1))
     call atomic_axpy(at,-1.0_gp*w,-1.0_gp,g, w(iypt+npt+1))
!!!     do i=1,n
!!!        w(ispt+npt+i)=a_t*w(ispt+npt+i)
!!!        w(iypt+npt+i)=-g(i)-w(i)
!!!     enddo
     !point=point+1
     !if(point==m) point=0
     point=mod(point+1,m)
     !termination test
     call atomic_dot(at,g,g,gnorm)
     call atomic_dot(at,x,x,xnorm)
     gnorm=sqrt(gnorm)
     xnorm=sqrt(xnorm)
!!!     gnorm=sqrt(ddot(n,g,1,g,1))
!!!     xnorm=sqrt(ddot(n,x,1,x,1))
     xnorm=max(1.0d0,xnorm)
!     if(gnorm/xnorm<=parmin%eps) finish=.true.
     if(parmin%iprint(1)>=0) then
        if (iproc==0)  call lb1(nfun,gnorm,n,m,x,f,g,a_t,finish,parmin)
        !Keep correct transformed positions
        do i=1,n
           xc(i)=x(i)
        enddo
        iwrite=iwrite+1
     endif
     if(finish) then
        parmin%iflag=0
        return
     endif
     new=.true.
  enddo
END SUBROUTINE lbfgs
!*****************************************************************************************
subroutine init_lbfgs(at,n,m,g,diag,w,parmin,nfun,point,finish,stp1,ispt,iypt)
  use module_base
  use module_types
  use minpar, only:parameterminimization
  implicit none
  type(atoms_data), intent(in) :: at
  type(parameterminimization)::parmin
  integer::n,m,i
  real(8)::g(n),diag(n),w(n*(2*m+1)+2*m)
  integer::nfun,point,iypt,ispt
  real(8)::one,zero,gnorm,stp1
  logical::finish
  data one,zero/1.0d+0,0.0d+0/
  parmin%iter=0
  if(n<1 .or. m<1) then
     parmin%iflag= -3
     write(*,'(a)') 'iflag= -3, improper input parameters(n or m are less than 1)'
     return
  endif
  if(parmin%gtol<1.d-4) then
     if(parmin%lp>0) then
        write(*,'(a)') 'gtol is less than 1.d-4, it has been reset to 9.d-1'
     endif
     parmin%gtol=9.d-01
  endif
  nfun=1
  point=0
  finish=.false.
  if(parmin%diagco) then
     do i=1,n
        if(diag(i)<=zero) then
           parmin%iflag=-2
           write(*,'(a,i5,2a)') 'iflag=-2, the',i, '-th diagonal element of ', &
                'the inverse hessian approximation is not positive'
           return
        endif
     enddo
  else
     do i=1,n
        diag(i)= 1.d0
     enddo
  endif
  ispt=n+2*m
  iypt=ispt+n*m     

  !to be ussed only in the case parmin%diagco==.false.
  !in that case diag=1.
!!!  if (.not. parmin%diagco) then
!!!     !this line equals to w=g
!!!     call atomic_axpy(at,g,0.0_gp,g,w(ispt+1))
!!!  else
     !here the blocked atoms are treated as the others
  do i=1,n
     w(ispt+i)=g(i)*diag(i)
  enddo
!!!  end if
 
  call atomic_dot(at,g,g,gnorm)
  gnorm=dsqrt(gnorm)
!!!  gnorm=dsqrt(ddot(n,g,1,g,1))

  !    stp1=one/gnorm
  

  !stp1=2.d-2/gnorm  !original convention
  stp1=parmin%betax
end subroutine init_lbfgs


!!****f* BigDFT/lb1
!! FUNCTION
!!  This routine prints monitoring information. the frequency and
!!  amount of output are controlled by iprint.
!! SOURCE
!!
subroutine lb1(nfun,gnorm,n,m,x,f,g,a_t,finish,parmin_)
  use minpar, only: parameterminimization
  implicit none
  type(parameterminimization) :: parmin_
  integer::nfun,n,m,i
  real(8)::x(n),g(n),f,gnorm,a_t
  logical finish
  if(parmin_%iter==0)then
     write(parmin_%mp,'(a)') '*************************************************'
     write(parmin_%mp,'(a,i5,a,i2,a)') 'n=',n,' number of corrections=',m,' initial values'
     write(parmin_%mp,'(a,1pd10.3,a,1pd10.3)') ' f= ',f,'   gnorm= ',gnorm
     if (parmin_%iprint(2)>=1) then
        write(parmin_%mp,'(a)') ' vector x= '
        write(parmin_%mp,'(6(2x,1pd10.3))') (x(i),i=1,n)
        write(parmin_%mp,'(a)') ' gradient vector g= '
        write(parmin_%mp,'(6(2x,1pd10.3))') (g(i),i=1,n)
     endif
     write(parmin_%mp,'(a)') '*************************************************'
     write(parmin_%mp,'(a)') '   iter    nfn  func  gnorm steplength'
  else
     if((parmin_%iprint(1)==0) .and. (parmin_%iter/=1 .and. .not.finish)) return
     if(parmin_%iprint(1)/=0)then
        if(mod(parmin_%iter-1,parmin_%iprint(1)).eq.0.or.finish)then
           if(parmin_%iprint(2)>1 .and. parmin_%iter>1) write(parmin_%mp,'(a)') & 
                '   iter    nfn  func  gnorm steplength'
           write(16,'(i5,1x,e12.5,1x,e21.14,a,i5,1x,e12.5)') parmin_%iter,gnorm,f,' GEOPT BFGS ', nfun,a_t
           write(parmin_%mp,'(a,2(i4,1x),3x,3(1E24.15,2x))') 'MIN ',parmin_%iter,nfun,f,gnorm,a_t
        else
           return
        endif
     else
        if(parmin_%iprint(2)>1 .and. finish) write(parmin_%mp,'(a)') &
             '   iter    nfn  func  gnorm steplength'
        write(parmin_%mp,'(2(i4,1x),3x,3(1pd24.15,2x))') parmin_%iter,nfun,f,gnorm,a_t
     endif
     if(parmin_%iprint(2)==2 .or. parmin_%iprint(2)==3)then
        if(finish)then
           write(parmin_%mp,'(a)') ' final point x= '
        else
           write(parmin_%mp,'(a)') ' vector x= '
        endif
        write(parmin_%mp,'(6(2x,1pd10.3))')(x(i),i=1,n)
        if(parmin_%iprint(2)==3) then
           write(parmin_%mp,'(a)') ' gradient vector g= '
           write(parmin_%mp,'(6(2x,1pd10.3))')(g(i),i=1,n)
        endif
     endif
     if(finish) write(parmin_%mp,'(a)') &
          ' the minimization terminated without detecting errors. iflag = 0'
     if(finish) write(16,*) &
          ' BFGS terminated without detecting errors. iflag = 0'
  endif
  return
END SUBROUTINE lb1
!!***

!!****f* BigDFT/mcsrch
!! FUNCTION
!!  Line search routine
!! SOURCE
!!
subroutine mcsrch(at,n,x,f,g,s,a_t,info,nfev,wa,parmin)
  use module_base
  use module_types
  use minpar, only: parameterminimization,infoc,brackt,stage1,dg,&
       & dginit,dgtest,dgx,dgy,finit,fx,fy,a_l,a_u,stmin,stmax,width,width1
  implicit none
  type(atoms_data), intent(in) :: at
  type(parameterminimization), intent(inout) :: parmin
  integer, intent(inout) :: n,info,nfev
  real(8), intent(inout) :: f,a_t
  real(8), intent(inout) :: x(n),g(n),s(n),wa(n)
  logical :: yes
  integer :: j
  real(8) :: fm, fxm, fym, dgm, dgxm, dgym, ftest1
  real(8), parameter :: p5 = 0.5d0, p66 = 0.66d0, xtrapf = 4.0d0, zero = 0.d0

  yes=.true.
  if(info==-1) yes=.false.
  if(yes) then
     infoc = 1
     !check the input parameters for errors.
     if(n<1 .or. a_t<=zero .or. parmin%ftol<zero .or. parmin%gtol<zero .or. &
          parmin%xtol<zero .or. parmin%stpmin<zero .or. parmin%stpmax<parmin%stpmin .or. parmin%maxfev<1) return
     !compute the initial gradient in the search direction
     !and check that s is a descent direction.
     call atomic_dot(at,g,s,dginit)
     dginit=-dginit

!!!     dginit = zero
!!!     do j = 1, n
!!!        dginit = dginit - g(j)*s(j)
!!!     enddo
     if (dginit >=  zero) then
        write(parmin%lp,'(a,1pe24.17)') 'the search direction is not a descent direction',dginit
        return
     endif
     !initialize local variables.
     brackt=.false.
     stage1=.true.
     nfev=0
     finit=f
     dgtest=parmin%ftol*dginit
     width=parmin%stpmax - parmin%stpmin
     width1=width/p5
     do j = 1,n
        wa(j)=x(j)
     enddo
     a_l=zero
     fx=finit
     dgx=dginit
     a_u=zero
     fy=finit
     dgy=dginit
  endif
  !start of iteration.
  do
     if(yes) then
        !set the minimum and maximum steps to correspond
        !to the present interval of uncertainty.
        if(brackt) then
           stmin=min(a_l,a_u)
           stmax=max(a_l,a_u)
        else
           stmin=a_l
           stmax=a_t+xtrapf*(a_t-a_l)
        end if
        !force the step to be within the bounds stpmax and stpmin.
        a_t=max(a_t,parmin%stpmin)
        a_t=min(a_t,parmin%stpmax)
        !if an unusual termination is to occur then let
        !a_t be the lowest point obtained so far.
        if((brackt .and. (a_t .le. stmin .or. a_t .ge. stmax)) &
             .or. nfev .ge. parmin%maxfev-1 .or. infoc .eq. 0  &
             .or. (brackt .and. stmax-stmin .le. parmin%xtol*stmax)) a_t = a_l
        !evaluate the function and gradient at a_t
        !and compute the directional derivative.
        !we return to main program to obtain f and g.
        call atomic_axpy(at,wa,a_t,s,x)
!!!        do j = 1, n
!!!           x(j) = wa(j) + a_t*s(j)
!!!        enddo
        info=-1
        return
     endif
     info=0
     nfev = nfev + 1
     call atomic_dot(at,g,s,dg)
     dg=-dg
!!!     dg = zero
!!!     do j=1,n
!!!        dg=dg-g(j)*s(j)
!!!     enddo
     ftest1=finit+a_t*dgtest
     !test for convergence.
     if((brackt .and. (a_t<=stmin .or. a_t>=stmax))   .or. infoc==0) info=6
     if(a_t==parmin%stpmax .and. f<=ftest1 .and. dg<=dgtest) info=5
     if(a_t==parmin%stpmin .and. (f>ftest1 .or. dg>=dgtest)) info=4!; write(*,*) "After test for convergence"
     if(nfev>=parmin%maxfev) info=3
     if(brackt .and. stmax-stmin<=parmin%xtol*stmax) info=2
     if(f<=ftest1 .and. abs(dg)<=parmin%gtol*(-dginit)) info=1
     !check for termination.
     if(info/=0) return
     !in the first stage we seek a step for which the modified
     !function has a nonpositive value and nonnegative derivative.
     if(stage1 .and. f<=ftest1 .and. dg>=min(parmin%ftol,parmin%gtol)*dginit) stage1=.false.
     !a modified function is used to predict the step only if
     !we have not obtained a step for which the modified
     !function has a nonpositive function value and nonnegative
     !derivative, and if a lower function value has been
     !obtained but the decrease is not sufficient.
     if(stage1 .and. f<=fx .and. f>ftest1) then
        !define the modified function and derivative values.
        fm=f-a_t*dgtest
        fxm=fx-a_l*dgtest
        fym=fy-a_u*dgtest
        dgm=dg-dgtest
        dgxm=dgx-dgtest
        dgym=dgy-dgtest
        !call mcstep to update the interval of uncertainty and to compute the new step.
        call mcstep(a_l,fxm,dgxm,a_u,fym,dgym,a_t,fm,dgm, brackt,stmin,stmax,infoc) !,parmin)
        !reset the function and gradient values for f.
        fx=fxm+a_l*dgtest
        fy=fym+a_u*dgtest
        dgx=dgxm+dgtest
        dgy=dgym+dgtest
     else
        !call mcstep to update the interval of uncertainty and to compute the new step.
        call mcstep(a_l,fx,dgx,a_u,fy,dgy,a_t,f,dg,brackt,stmin,stmax,infoc) !,parmin)
     end if
     !force a sufficient decrease in the size of the
     !interval of uncertainty.
     if(brackt) then
        if(abs(a_u-a_l)>=p66*width1) a_t=a_l+p5*(a_u-a_l)
        width1=width
        width=abs(a_u-a_l)
     end if
     yes=.true.
  enddo
END SUBROUTINE mcsrch
!!***

subroutine mcstep(a_l,fx,dx,a_u,fy,dy,a_t,fp,dp,brackt,stpmin,stpmax,info) !,parmin)
  use minpar, only:parameterminimization
  implicit none
!  type(parameterminimization) :: parmin
  integer::info
  real(8)::a_l,fx,dx,a_u,fy,dy,a_t,fp,dp,stpmin,stpmax
  logical::brackt,bound
  real(8)::gamma,p,q,r,s,sgnd,a_c,stpf,a_q,theta
  info = 0
  !check the input parameters for errors.
  if ((brackt .and. (a_t<=min(a_l,a_u) .or.   a_t>=max(a_l,a_u))) .or. &
       dx*(a_t-a_l)>=0.d0 .or. stpmax<stpmin) return
  !determine if the derivatives have opposite sign.
  sgnd = dp*(dx/abs(dx))
  if(fp>fx) then
     !first case. a higher function value.
     !the minimum is bracketed. if the cubic step is closer
     !to a_l than the quadratic step, the cubic step is taken,
     !else the average of the cubic and quadratic steps is taken.
     info = 1
     bound = .true.
     call cal_a_c(a_l,fx,dx,a_t,fp,dp,a_c)
     a_q = a_l + ((dx/((fx-fp)/(a_t-a_l)+dx))/2)*(a_t - a_l)
     if (abs(a_c-a_l) .lt. abs(a_q-a_l)) then
        stpf = a_c
     else
        stpf = a_c + (a_q - a_c)/2
     end if
     brackt = .true.
  else if (sgnd .lt. 0.0) then
     !second case. a lower function value and derivatives of
     !opposite sign. the minimum is bracketed. if the cubic
     !step is closer to a_l than the quadratic (secant) step,
     !the cubic step is taken, else the quadratic step is taken.
     info=2
     bound=.false.
     call cal_a_c_2(a_l,fx,dx,a_t,fp,dp,a_c)
     a_q=a_t+(dp/(dp-dx))*(a_l-a_t)
     if (abs(a_c-a_t)>abs(a_q-a_t)) then
        stpf=a_c
     else
        stpf=a_q
     end if
     brackt=.true.
  elseif(abs(dp)<abs(dx)) then
     !third case. a lower function value, derivatives of the
     !same sign, and the magnitude of the derivative decreases.
     !the cubic step is only used if the cubic tends to infinity
     !in the direction of the step or if the minimum of the cubic
     !is beyond a_t. otherwise the cubic step is defined to be
     !either stpmin or stpmax. the quadratic (secant) step is also
     !computed and if the minimum is bracketed then the the step
     !closest to a_l is taken, else the step farthest away is taken.
     info=3
     bound=.true.
     call cal_a_c_3(a_l,fx,dx,a_t,fp,dp,stpmin,stpmax,a_c)
     a_q=a_t+(dp/(dp-dx))*(a_l-a_t)
     if (brackt) then
        if (abs(a_t-a_c)<abs(a_t-a_q)) then
           stpf=a_c
        else
           stpf=a_q
        end if
     else
        if (abs(a_t-a_c)>abs(a_t-a_q)) then
           stpf=a_c
        else
           stpf=a_q
        endif
     end if
  else
     !fourth case. a lower function value, derivatives of the
     !same sign, and the magnitude of the derivative does
     !not decrease. if the minimum is not bracketed, the step
     !is either stpmin or stpmax, else the cubic step is taken.
     info=4; write(*,*) "Fourth case scenario"
     bound=.false.
     if(brackt) then
        theta=3.d0*(fp-fy)/(a_u-a_t)+dy+dp
        s=max(abs(theta),abs(dy),abs(dp))
        gamma=s*sqrt((theta/s)**2-(dy/s)*(dp/s))
        if(a_t>a_u) gamma=-gamma
        p=(gamma-dp)+theta
        q=((gamma-dp)+gamma)+dy
        r=p/q
        a_c=a_t + r*(a_u-a_t)
        stpf=a_c
     elseif(a_t>a_l) then
        stpf=stpmax
     else
        stpf=stpmin
     end if
  end if
  !update the interval of uncertainty. this update does not
  !depend on the new step or the case analysis above.
  if(fp>fx) then
     a_u=a_t
     fy=fp
     dy=dp
  else
     if(sgnd<0.d0) then
        a_u=a_l
        fy=fx
        dy=dx
     end if
     a_l=a_t
     fx=fp
     dx=dp
  end if
  !compute the new step and safeguard it.
  stpf=min(stpmax,stpf)
  stpf=max(stpmin,stpf)
  a_t=stpf
  if(brackt .and. bound) then
     if(a_u>a_l) then
        a_t=min(a_l+0.66d0*(a_u-a_l),a_t)
     else
        a_t=max(a_l+0.66d0*(a_u-a_l),a_t)
     end if
  end if
  return
END SUBROUTINE mcstep
!*****************************************************************************************
subroutine cal_a_c(a_l,fx,dx,a_t,fp,dp,a_c)
  implicit none
  real(8)::a_l,fx,dx,a_t,fp,dp,a_c,theta,gamma,r,q,p
  theta=3.d0*(fx-fp)/(a_t-a_l)+dx+dp
  !s=max(abs(theta),abs(dx),abs(dp))
  gamma=sqrt(theta**2-dx*dp)
  !gamma=s*sqrt((theta/s)**2-(dx/s)*(dp/s))
  if(a_t<a_l) gamma=-gamma
  p=(gamma-dx)+theta
  q=((gamma-dx)+gamma)+dp
  r=p/q
  a_c=a_l+r*(a_t-a_l)
  !---------------------------------------------------------------
  !p=-((2.d0*fx-2.d0*fp-dx*a_l-dp*a_l+dx*a_t+dp*a_t)/(a_l-a_t)**3)
  !q=-((3.d0*fx-3.d0*fp-2.d0*dx*a_l-dp*a_l+2.d0*dx*a_t+dp*a_t)/(a_l-a_t)**2)
  !a_c=a_l+(-q+sqrt(q**2-3.d0*p*dx))/(3.d0*p)
END SUBROUTINE cal_a_c
!*****************************************************************************************
subroutine cal_a_c_2(a_l,fx,dx,a_t,fp,dp,a_c)
  implicit none
  real(8)::a_l,fx,dx,a_t,fp,dp,a_c,theta,s,gamma,r,q,p
  theta=3.d0*(fx-fp)/(a_t-a_l)+dx+dp
  s=max(abs(theta),abs(dx),abs(dp))
  gamma=s*sqrt((theta/s)**2-(dx/s)*(dp/s))
  if(a_t>a_l) gamma=-gamma
  p=(gamma-dp)+theta
  q=((gamma-dp)+gamma)+dx
  r=p/q
  a_c=a_t+r*(a_l-a_t)
END SUBROUTINE cal_a_c_2
!*****************************************************************************************
subroutine cal_a_c_3(a_l,fx,dx,a_t,fp,dp,stpmin,stpmax,a_c)
  implicit none
  real(8)::a_l,fx,dx,a_t,fp,dp,stpmin,stpmax,a_c,theta,s,gamma,r,q,p
  theta=3.d0*(fx-fp)/(a_t-a_l)+dx+dp
  s=max(abs(theta),abs(dx),abs(dp))
  gamma=s*sqrt(max(0.d0,(theta/s)**2-(dx/s)*(dp/s)))
  if(a_t>a_l) gamma=-gamma
  p=(gamma-dp)+theta
  q=(gamma+(dx-dp))+gamma
  r=p/q
  !the case gamma = 0 only arises if the cubic does not tend
  !to infinity in the direction of the step.
  if(r<0.d0 .and. gamma/=0.d0) then
     a_c=a_t+r*(a_l-a_t)
  elseif(a_t>a_l) then
     a_c=stpmax
  else
     a_c=stpmin
  endif
END SUBROUTINE cal_a_c_3
!*****************************************************************************************

!!****f* BigDFT/rundiis
!! FUNCTION
!!  DIIS relax. Original source from ART from N. Mousseau.
!!  Adaptations to BigDFT by D. Caliste.
!! SOURCE
!!
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
  integer :: lter, maxter, i, i_err, n, nrhs, lwork, infocode, j, i_stat, i_all,nfluct
  real(gp) :: sumx, sumy, sumz, fluctsum, fluct, fmax, fnrm
  character(len = 4) :: fn4
  character(len = 40) :: comment
  ! Local variables for Lapack.
  integer, dimension(:), allocatable :: interchanges
  real(8), dimension(:), allocatable :: work
  real(8), dimension(:), allocatable :: solution
  real(8), dimension(:,:), allocatable :: matrice
  logical :: check

  ! We save pointers on data used to call bigdft() routine.
  allocate(previous_forces(in%history, AT%NAT * 3+ndebug),stat=i_stat)
  call memocc(i_stat,previous_forces,'previous_forces',subname)
  allocate(previous_pos(in%history, AT%NAT * 3+ndebug),stat=i_stat)
  call memocc(i_stat,previous_pos,'previous_pos',subname)
  allocate(product_matrix(in%history, in%history+ndebug),stat=i_stat)
  call memocc(i_stat,product_matrix,'product_matrix',subname)

  fluctsum = 0._gp
  nfluct = 0

  ! We set the first step and move to the second
  previous_forces(1,:) = f(:)
  previous_pos(1,:) = x(:)

  !x(:) = x(:) + in%betax * f(:)
  !always better to use the atomic_* routines to move atoms
  !it performs modulo operation as well as constrained search
  call atomic_axpy(at,x,in%betax,f,x)

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
     if (lter .gt. maxter) then
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
     !reput the modulo operation on the atoms
     call atomic_axpy(at,x,0.0_gp,x,x)

     i_all=-product(shape(solution))*kind(solution)
     deallocate(solution,stat=i_stat)
     call memocc(i_stat,i_all,'solution',subname)

     in%inputPsiId=1
     call call_bigdft(nproc,iproc,at,x,in,epot,f,rst,infocode)

     if (iproc == 0) then
        call transforce(at,f,sumx,sumy,sumz)
        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
     end if

     ncount_bigdft=ncount_bigdft+1

     call fnrmandforcemax(f,fnrm,fmax,at)
     if (fmax < 3.d-1) call updatefluctsum(at,f,nfluct,fluctsum,fluct)

     if (iproc==0) then 
        write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=', &
             & fmax,'fnrm=',    fnrm    ,'fluct=', fluct
        write(fn4,'(i4.4)') ncount_bigdft
        write(comment,'(a,1pe10.3)')'DIIS:fnrm= ',sqrt(fnrm)
        call write_atomic_file('posout_'//fn4,epot,x,at,trim(comment))
     endif

     call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

     if(check)then
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

     !x(:) = x(:) + in%betax * f(:)
     call atomic_axpy(at,x,in%betax,f,x)
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
end subroutine rundiis
!!***
