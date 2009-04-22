!!****p* MINHOP
!!
!! DESCRIPTION
!!  Main program fro the minima hopping
!!
!! COPYRIGHT
!!    Copyright (C) 2008 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS
!!
!! SOURCE
!!
program MINHOP

  use module_base
  use module_types
  use module_interfaces

  implicit real(kind=8) (a-h,o-z)
  real(kind=4) :: tts
  logical :: newmin
  character(len=20) :: units,atmn
  character(len=80) :: line
  type(atoms_data) :: atoms
  type(input_variables) :: inputs_opt, inputs_md
  type(restart_objects) :: rst
  ! some variable definitions are not needed anymore due to the new restart variable type
  ! type(wavefunctions_descriptors) :: wfd

  logical, dimension(:), allocatable :: lfrztyp
  !  real(kind=8), dimension(:,:), allocatable :: rxyz_old

  !  real(kind=8), dimension(:), pointer :: eval
  !  real(kind=8), dimension(:), pointer :: psi

  integer, parameter :: npminx=200
  !C parameters for minima hopping
  integer, parameter :: mdmin=2
  real(kind=8), parameter :: beta1=1.10d0,beta2=1.10d0,beta3=1.d0/1.10d0
  real(kind=8), parameter :: alpha1=1.d0/1.10d0,alpha2=1.10d0
  integer, parameter :: nbuf=1000
  ! interval for writing intermediate results
  integer, parameter :: minter=1
  real(kind=8) :: elocmin(npminx),abuf(nbuf)
  real(kind=8), allocatable, dimension(:,:) :: pos,ff,wpos,vxyz,gg,earr
  real(kind=8), allocatable, dimension(:,:,:) :: poslocmin

  character(len=*), parameter :: subname='global'
  character(len=41) :: filename
  character(len=4) :: fn4
  character(len=5) :: fn5
  character(len=50) :: comment

  !include 'mpif.h'

  ! Start MPI version
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  !call system('echo $HOSTNAME')

  !initialize memory counting
  call memocc(0,iproc,'count','start')

  if (iproc.eq.0)then
     write(*,'(23x,a)')' '
     write(*,'(23x,a)')'      __  __ _ _  _ _   _  __  ___ '
     write(*,'(23x,a)')'     |  \/  |_| \| | |_| |/  \| _ \ '
     write(*,'(23x,a)')'     | |\/| |-|    |  _  | <> |  _/ '
     write(*,'(23x,a)')'     |_|  |_|_|_|\_|_| |_|\__/|_|     WITH'
     write(*,'(23x,a)')''
     write(*,'(23x,a)')''
     write(*,'(23x,a)')''
     call print_logo()
     write(*,'(23x,a)')'----> you can grep this file for # to compare with global.out'
     write(*,'(23x,a)')' #NOTE: this version reads nspin, mpol  from input.dat'
  end if

  open(unit=67,file='global.out')

  write(filename,'(A,i3.3,a)')'latest.pos.force.',iproc,'.dat'
  open(unit=111,file=filename)


  if (iproc.eq.0) write(67,'(a,3(1x,1pe11.4))') 'beta1,beta2,beta3',beta1,beta2,beta3
  if (iproc.eq.0) write(67,'(a,2(1x,1pe11.4))') 'alpha1,alpha2',alpha1,alpha2
  ratio=-log(alpha2)/log(alpha1)
  if (iproc.eq.0) write(67,'(a,2(1x,1pe10.3))') 'predicted fraction accepted, rejected', & 
       ratio/(1.d0+ratio), 1.d0/(1.d0+ratio)
  tt= beta2**ratio*beta3
  if (iproc.eq.0) write(67,'(a,2(1x,1pe11.4))') 'critical ratio (.ge. 1)',tt
  if (tt.lt.0.999999999999999d0) stop 'incompatible alpha s, beta s'
  if (iproc.eq.0) write(67,*) 'mdmin',mdmin


  call cpu_time(tcpu1)
  ! read  earr.dat
  open(unit=12,file='earr.dat',status='unknown')
  read(12,*) nlmin,nlminx
  read(12,*) eref
  if (iproc.eq.0) write(67,*) 'eref=',eref
  if (iproc.eq.0) write(7,*) 'eref=',eref
  read(12,*) accur
  if (iproc.eq.0) write(67,*) 'accuracy for rounding=',accur
  if (nlmin.gt.nlminx) stop 'nlmin>nlminx'
  allocate(earr(0:nlminx+nbuf,2+ndebug),stat=i_stat)
  call memocc(i_stat,earr,'earr',subname)
  earr(0,1)=-1.d100
  if (nlmin.eq.0) then 
     if (iproc.eq.0) write(67,*) 'New run with nlminx=',nlminx
     if (iproc.eq.0) write(*,*) '#New run with nlminx=',nlminx
  else
     if (iproc.eq.0) write(67,*) 'Restart run with nlmin, nlminx=',nlmin,nlminx
     if (iproc.eq.0) write(*,*) '#Restart run with nlmin, nlminx=',nlmin,nlminx
     do k=1,nlmin
        read(12,*) earr(k,1),earr(k,2)
        if (earr(k,1).lt.earr(k-1,1)) stop 'wrong ordering in earr.dat'
     enddo
     if (iproc.eq.0) write(67,*) ' read earr.dat'
     if (iproc.eq.0) write(*,*) '# read earr.dat'
  endif
  close(12)

  filename = 'poscur.xyz'
  open(unit=99,file=filename,status='old')
  if (nlmin.eq.0) then 
     read(99,*) atoms%nat,atoms%units
     tre_pos=1000.d0
  else
     read(99,*) atoms%nat,atoms%units,tre_pos
  endif
  if (iproc.eq.0) write(67,*) 'nat=',atoms%nat
  if (iproc.eq.0) write(*,*) '#nat=',atoms%nat

  allocate(pos(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,pos,'pos',subname)

  if (iproc.eq.0) write(67,'(a,a)') 'reading positions from ',filename
  if (iproc.eq.0) write(*,'(a,a)') ' # reading positions from ',filename

  call read_atomic_positions(iproc,99,atoms,pos)
  close(unit=99)
  !Read input parameters for geometry optimization 
  call read_input_variables(iproc,'input.dat',inputs_opt)
  call read_input_variables(iproc,'mdinput.dat',inputs_md)

  ! allocate other arrays
  allocate(ff(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,ff,'ff',subname)
  allocate(wpos(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,wpos,'wpos',subname)
  allocate(vxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,vxyz,'vxyz',subname)
  allocate(gg(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gg,'gg',subname)
  allocate(poslocmin(3,atoms%nat,npminx+ndebug),stat=i_stat)
  call memocc(i_stat,poslocmin,'poslocmin',subname)
  !        allocate(rxyz_old(3,atoms%nat+ndebug),stat=i_stat)
  !        call memocc(i_stat,rxyz_old,'rxyz_old','BigDFT')


  ! read input parameters
  open(unit=11,file='ioput',status='old')
  if (iproc.eq.0) write(67,'(a,a)') 'reading input from ',filename
  if (iproc.eq.0) write(*,'(a,a)') ' # reading input from ',filename
  read(11,*) ediff,ekinetic,dt
  read(11,*,iostat=ierr)irestart
  if(ierr/=0)irestart=0
  read(11,*,iostat=ierr)av_ekinetic
  if(ierr==0)then
     if (iproc.eq.0) write(*,'(a,a)') ' # reading full restart variables'
     read(11,*)av_ediff
     read(11,*)escape
     read(11,*)escape_sam
     read(11,*)escape_old
     read(11,*)escape_new
     read(11,*)hopp
     read(11,*)hopp_acc
     read(11,*)hopp_rej
     read(11,*)egap
     read(11,*)esep
     read(11,*)count_sd
     read(11,*)count_cg
     read(11,*)count_md
     read(11,*)nputback
  else
     if (iproc.eq.0) write(*,'(a,a)') ' # initialize counter variables'
     !C first local minimum
     av_ekinetic=0.d0
     av_ediff=0.d0
     escape=0.d0
     escape_sam=0.d0
     escape_old=0.d0
     escape_new=0.d0
     hopp=0.d0
     hopp_acc=0.d0
     hopp_rej=0.d0
     egap=1.d100
     esep=0.d0
     count_sd=0.d0
     count_cg=0.d0
     count_md=0.d0
     nputback=0
  end if
  close(11)

  open(unit=16,file='geopt.mon',status='unknown')

  open(unit=11,file='rand.inp',status='old')
  read(11,*) nrandoff
  close(11)
  if (iproc == 0) write(*,*) '# nrandoff ',nrandoff

  do  i=1,nrandoff
     call random_number(tts)
  enddo
  if (iproc == 0) write(67,'(a,1x,3(1x,1pe10.3))') 'In :ediff,ekinetic,dt',ediff,ekinetic,dt
  if (iproc == 0) write(*,'(a,1x,3(1x,1pe10.3))') ' # In :ediff,ekinetic,dt',ediff,ekinetic,dt


  ! If restart run read previously found energies
  nlmin_l=0
  if (nlmin > 0) then
     kk=1
     do
        write(fn5,'(i5.5)') kk
        filename = 'poslow'//fn5//'.xyz'
        if (nlmin_l.ge.npminx .or. kk.gt.min(nlmin,npminx)) then
           exit
        end if
        nlmin_l=nlmin_l+1
        open(unit=9,file=filename,status='old',iostat=ierror)
        if (ierror /= 0) then
           write(*,*) iproc,' COULD not read file ',filename
           exit
        end if
        read(9,*) natp,units,elocmin(nlmin_l)
        if (atoms%nat.ne.natp) stop   'nat <> natp'
        read(9,*) 
        do iat=1,atoms%nat
           read(9,*) atmn,poslocmin(1,iat,nlmin_l),poslocmin(2,iat,nlmin_l),poslocmin(3,iat,nlmin_l)
        enddo
        close(9)
        if (iproc.eq.0) write(67,*) 'read file',filename
        if (iproc.eq.0) write(*,*) '# read file',filename
        kk=kk+1
     end do
     if (iproc.eq.0) write(67,*) 'read ',nlmin_l,'poslow files'
     if (iproc.eq.0) write(*,*) '# read ',nlmin_l,'poslow files'
  endif

  ! open output files
  if (iproc.eq.0) then
     !write(fn,'(i3.3)') iproc
     filename = 'global.mon'
     write(67,*) 'iteration results written in:',filename
     write(*,*) '# iteration results written in:',filename
     open(unit=2,file=filename,status='unknown',position='append')
  endif

  ! (un)comment the correct 2 lines
  !        alphax=2.d-2  ! Si
  !        if (iproc.eq.0) write(67,*) alphax,'optimal stepsize for Si systems'
  !        alphax=1.d-3  ! LJ
  !       if (iproc.eq.0) write(67,*) alphax,'optimal stepsize for LJ systems'

  inputs_opt%inputPsiId=0


  call init_restart_objects(atoms,rst,subname)
  ! new way of initializing data

  if (atoms%geocode.eq.'P') & 
       call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,pos)
  call call_bigdft(nproc,iproc,atoms,pos,inputs_opt,e_pos,ff,rst,infocode)


  write(17,*) 'ENERGY ',e_pos
  energyold=1.d100
  ncount_bigdft=0


  !        if(irestart>0)then
  !          if (iproc.eq.0) write(*,*)&
  !              ' # WARNING: This version skips relaxation of poscur.xyz! '
  !           if (iproc.eq.0) write(*,*)&
  !              ' # Ready to continue an aborted optimization. Reading posout file ',irestart
  !           call rdposout(irestart,wpos,atoms%nat)
  !           ncout_cluster=irestart
  !           re_sm=earr(1,1)
  !           goto 5556
  !        end if

  if (iproc.eq.0) write(*,*)'# calling conjgrad for the first time here. energy ',e_pos

  if (atoms%geocode.eq.'P') & 
       call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,pos)

    call geopt(nproc,iproc,pos,atoms,ff,e_pos,rst,inputs_opt,ncount_bigdft)

  if (iproc.eq.0) write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for accurate initial conf'
  if (iproc.eq.0) write(*,*)'# ', ncount_bigdft,' Wvfnctn Opt. steps for accurate initial conf'
  nconjgr=0
  if (iproc == 0) then 
     tt=dnrm2(3*atoms%nat,ff,1)

     !call wtlmin(nconjgr,atoms%nat,e_pos-eref,tt,pos,atoms%iatype,atoms%atomnames,atoms%natpol)

     write(fn4,'(i4.4)') nconjgr
     write(comment,'(a,1pe10.3)')'fnrm= ',tt
     call wtxyz('poslocm_'//fn4,e_pos-eref,pos,atoms,trim(comment))

  endif
  nconjgr=nconjgr+1

  re_pos=round(e_pos-eref,accur)
  if (iproc.eq.0) write(67,'(a,1x,i3,3(1x,1pe17.10))') 'INPUT(relaxed): iproc, e_pos,re_pos,eref ',iproc,e_pos,re_pos,eref
  if (iproc.eq.0) write(*,'(a,1x,i3,3(1x,1pe17.10))') ' # INPUT(relaxed): iproc, e_pos,re_pos,eref ',iproc,e_pos,re_pos,eref


  if (nlmin.gt.0) then
     if (iproc.eq.0) write(67,'(a,2(1x,1pe24.17))') 'new/old energy for input file',re_pos,tre_pos
     if (iproc.eq.0) write(*,'(a,2(1x,1pe24.17))') '# new/old energy for input file',re_pos,tre_pos
     if (re_pos.ne.tre_pos .and. iproc.eq.0) write(67,*) 'WARNING: new/old energies for input file differ'
     if (re_pos.ne.tre_pos .and. iproc.eq.0) write(*,*) '# WARNING: new/old energies for input file differ'
  endif
  k_e_wpos=1
  if (nlmin.eq.0) then
     nlmin=1
     nlmin_l=1
     earr(1,1)=re_pos
     earr(1,2)=1.d0
     elocmin(1)=re_pos
     do iat=1,atoms%nat
        poslocmin(1,iat,1)=pos(1,iat) 
        poslocmin(2,iat,1)=pos(2,iat) 
        poslocmin(3,iat,1)=pos(3,iat) 
     enddo
  endif
  re_sm=min(re_pos,earr(1,1))
  if (iproc.eq.0) write(67,*) 'iproc,initial re_sm',iproc,re_sm
  if (iproc.eq.0) write(*,*) '# iproc,initial re_sm',iproc,re_sm
  if (iproc.eq.0) write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3))') hopp,escape,e_pos-eref,ediff,ekinetic
  !        write(200+iproc,'(2(1x,f10.0),1x,1pe24.17,2(1x,1pe10.3))') hopp,escape,e_pos-eref,ediff,ekinetic


  if(irestart>0)then
     if (iproc.eq.0) write(*,*)&
          ' # Ready to continue an aborted optimization. Reading posout file ',irestart
     !presumably better to use read_atomic_positions
     call rdposout(irestart,wpos,atoms%nat)
     ncount_bigdft=irestart
     inputs_opt%inputPsiId=0  !ALEX notices we need an input guess for the aborted geo
     goto 5556
  end if
  
  
  do iat=1,atoms%nat
     wpos(1,iat)=pos(1,iat) 
     wpos(2,iat)=pos(2,iat) 
     wpos(3,iat)=pos(3,iat)
  end do
  nlmin_old=nlmin

  !C outer (hopping) loop
  hopping_loop: do
     !1000 continue
     irestart=0! we can restart from a local minimum
     rewind(111)
     write(111,*)'*** process monitor file for nlmin',nlim
     ! check whether other processes send local minima data
     if (nlmin >= nlminx) then 
        write(67,*)iproc,'has  nlminx collected',nlmin
        write(*,*)'#:process', iproc,'has  nlminx collected',nlmin
        !             do i=1,nlmin ;  write(*,*) earr(i,1) ; enddo
        !goto 3000
        exit hopping_loop 
     endif
     !            Energy has reached target eref and global minimum is presumably found
     if (re_sm <= 1.d-3) then
        write(*,*)'# process', iproc,'success: relative energy < 0.001'
        !goto 3000
        exit hopping_loop 
     endif

     ! write intermediate results
     irestart=0
     if (nlmin-nlmin_old.ge.minter) then
        if (iproc == 0) then 
           write(67,*) ' Writing intermediate results at',nlmin_l,nlmin
           write(*,*) '#  Writing intermediate results at',nlmin_l,nlmin
        endif
        call winter(iproc,atoms,re_pos,pos,npminx,nlminx,nbuf,nlmin,nlmin_l,accur, & 
             earr,elocmin,poslocmin,eref,ediff,ekinetic,dt)
        nlmin_old=nlmin
     endif

     ! check whether global minimum was found
     if (re_sm.le.1.d-3) then
        write(*,*) iproc,'SUCCESS success'
        write(*,*)'# process', iproc,'success: relative energy < 0.001'
        !goto 3000
        exit hopping_loop 
     endif

     if(iproc==0)then
        call timeleft(tt)
        write(*,*) '# CPU time hours left',tt
     end if
     !call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     !if (tt<0) then
     !write(*,*) '# CPU time limit exceeded for process',iproc
     !goto 3000
     !endif
     ! Escape loop
5555 continue

     if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
     !         sometimes, syncronization errors occured, so let's try this ^^^
     do iat=1,atoms%nat
        wpos(1,iat)=pos(1,iat)
        wpos(2,iat)=pos(2,iat) 
        wpos(3,iat)=pos(3,iat)
     enddo

     escape=escape+1.d0
     call mdescape(mdmin,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,wpos, &
          nproc,iproc,atoms,rst,inputs_md)

     call fix_fragmentation(iproc,atoms,wpos,nputback)

     if(iproc==0)call  timeleft(tt)
     if(iproc==0)write(*,*)'# CPU time hours left',tt
     av_ekinetic=av_ekinetic+ekinetic
     ncount_bigdft=0

5556 continue ! entry point for restart of optimization at cluster step irestart+1
     if (atoms%geocode.eq.'P') & 
          call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,wpos)

      call geopt(nproc,iproc,wpos,atoms,ff,e_wpos,rst,inputs_md,ncount_bigdft)

     if(iproc==0)then
        call timeleft(tt)
        write(*,*) '# CPU time hours left',tt
     end if
     call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     if (tt<0) then
        write(*,*)&
             '# CPU time exceeded during approximate geometry relaxation, process',iproc
        ! allows to continue the next run by relaxing the aborted escape attempt
        irestart=ncount_bigdft  
        exit hopping_loop 
     end if

     if (iproc.eq.0) write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
     if (iproc.eq.0) write(*,*)'# ', ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
     !ncount_bigdft=0
     !ncount_cluster=0
     if (atoms%geocode.eq.'P') & 
          call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,wpos)
      call geopt(nproc,iproc,wpos,atoms,ff,e_wpos,rst,inputs_opt,ncount_bigdft)

     if(iproc==0)then
        call timeleft(tt)
        write(*,*) '# CPU time hours left',tt
     end if
     call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     if (tt<0) then
        !if (iproc.eq.0) write(*,*) '# CPU time exceeded during accurate geometry relaxation'
        write(*,*) '# CPU time exceeded during accurate geometry relaxation, process',iproc
        irestart=ncount_bigdft  ! allows to continue the next run by relaxing the aborted escape attempt
        !goto 3000
        exit hopping_loop 
     end if
     if (iproc.eq.0) write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for accurate geo. rel of MD conf.'
     if (iproc.eq.0) write(*,*)'# ', ncount_bigdft,' Wvfnctn Opt. steps for accurate geo. rel of MD conf.'
     if (iproc == 0) then 
        tt=dnrm2(3*atoms%nat,ff,1)
        !call wtlmin(nconjgr,atoms%nat,e_wpos-eref,tt,wpos,atoms%iatype,atoms%atomnames,atoms%natpol)
        write(fn4,'(i4.4)') nconjgr
        write(comment,'(a,1pe10.3)')'fnrm= ',tt
        call wtxyz('poslocm_'//fn4,e_wpos-eref,wpos,atoms,trim(comment))


     endif
  nconjgr=nconjgr+1
  re_wpos=round(e_wpos-eref,accur)
  if (iproc.eq.0) write(67,'(a,i3,i3,4(1x,1pe14.7))')  & 
       'nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos', & 
       nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos
  if (iproc.eq.0) write(*,'(a,i3,i3,4(1x,1pe14.7))')  & 
       '# nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos', & 
       nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos
  !          write(200+iproc,'(a,i3,i3,4(1x,1pe14.7))')  & 
  !                    'nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos', & 
  !                     nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos
  !C not escaped
  if (re_pos.eq.re_wpos) then
     escape_sam=escape_sam+1.d0
     esep=esep+(e_pos-e_wpos)**2
     ekinetic=ekinetic*beta1
     if (iproc.eq.0) call wtioput(ediff,ekinetic,dt)
     if (iproc.eq.0) write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),a)')  &
          hopp,escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     !           write(200+iproc,'(2(1x,f10.0),1x,1pe24.17,2(1x,1pe10.3),3(1x,0pf5.2),a)')  &
     !                       hopp,escape,e_wpos-eref,ediff,ekinetic, &
     !                       escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     !write(100+iproc,*) 'SAME:re_pos,re_wpos ',re_pos,re_wpos
     !write(100+iproc,*) 'SAME ',e_wpos-eref,ediff,ekinetic
     if (iproc.eq.0) write(*,'(a)')' # no escape from current minimum.'
     goto 5555
  endif



  !C continue since escaped
  !C  check whether new minimum
  call hunt(earr(1,1),min(nlmin,nlminx+nbuf),re_wpos,k_e_wpos)
  if (re_wpos.eq.earr(k_e_wpos,1)) then
     if (iproc.eq.0) write(67,'(a,i3,i3,i4,1x,1pe14.7)')  & 
          'nlmin_l,nlmin,k_e_wpos,re_wpos=earr',nlmin_l,nlmin,k_e_wpos,re_wpos
     if (iproc.eq.0) write(*,'(a,i3,i3,i4,1x,1pe14.7)')  & 
          ' # Revisited: nlmin_l,nlmin,k_e_wpos,re_wpos=earr',nlmin_l,nlmin,k_e_wpos,re_wpos
     newmin=.false.
     escape_old=escape_old+1.d0
     ! standard feedback
     ekinetic=ekinetic*beta2
     ! enhanced feedback mechanism for ekinetic                                       
     !            ekinetic=ekinetic*beta2*(1.d0+.33d0*log(earr(k_e_wpos,2)))
     ! enhanced feedback mechanism for ediff                              
     !            ediff=ediff*alpha2**(1.d0+.1d0*log(earr(k_e_wpos,2)))
     !            ekinetic=ekinetic*(1.d0+(beta2-1.d0)*(1.d0+1.d-1*earr(k_e_wpos,2)))
  else
     newmin=.true.
     escape_new=escape_new+1.d0
     ekinetic=ekinetic*beta3
     ! determine energy separation between distinct minima
     if (k_e_wpos+1.le.nlminx)  then
        tt=min(re_wpos-earr(k_e_wpos,1),earr(k_e_wpos+1,1)-re_wpos)
        if (tt.gt. accur*(1.1d0)) egap=min(egap,tt)
        !       write(*,*) 'tt,egap',tt,egap
     endif
     if (iproc.eq.0) call wtioput(ediff,ekinetic,dt)

     if (k_e_wpos.eq.0) then
        re_sm=min(re_sm,re_wpos)
        if (iproc == 0) then 
           write(67,'(a,i7,1x,1pe24.17,1x,i2)') 'new lowest ',nlmin,re_wpos,iproc
           write(*,'(a,i7,1x,1pe24.17,1x,i2)') '# new lowest ',nlmin,re_wpos,iproc
           !call wtbest(atoms%nat,re_wpos,wpos, &
           !     atoms%iatype,atoms%atomnames,atoms%natpol)
           call wtxyz('posbest',re_wpos,wpos,atoms,'')
        endif
     else
        if (iproc.eq.0) write(*,'(a,i3)')' # New local minimum, lower are ',k_e_wpos 
     endif
  endif

  hopp=hopp+1.d0

  !C master: Monte Carlo step for local minima hopping
  !write(67,*) 'e_wpos,e_pos',e_wpos,e_pos
  if(e_wpos-e_pos.lt.ediff) then
     !C          local minima accepted -------------------------------------------------------
     if (iproc.eq.0) write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),l3,a)')  &
          hopp,escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A '
     !           write(200+iproc,'(2(1x,f10.0),1x,1pe24.17,2(1x,1pe10.3),3(1x,0pf5.2),l,a)')  &
     !                       hopp,escape,e_wpos-eref,ediff,ekinetic, &
     !                       escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A '
     hopp_acc=hopp_acc+1.d0
     ediff=ediff*alpha1
     if (iproc.eq.0) call wtioput(ediff,ekinetic,dt)
     av_ediff=av_ediff+ediff
     e_pos=e_wpos
     re_pos=re_wpos
     do iat=1,atoms%nat
        pos(1,iat)=wpos(1,iat) 
        pos(2,iat)=wpos(2,iat) 
        pos(3,iat)=wpos(3,iat)
     enddo
     if (newmin) then
        !C            if new local minimum
        nlmin=nlmin+1
        nlmin_l=nlmin_l+1
        !C            add minimum to history list
        call insert(nlminx,nbuf,nlmin,k_e_wpos,re_wpos,earr(0,1))
        !                write(67,*) iproc,'EARR'
        !                do i=1,nlmin ;  write(*,*) earr(i,1),earr(i,2) ; enddo
        !C            save configuration if it is among the lowest ones in energy
        call save_low_conf(atoms%nat,nlmin_l,npminx,re_wpos,pos,elocmin,poslocmin)
     else
        !C            old minimum revisited
        earr(k_e_wpos,2)=earr(k_e_wpos,2)+1.d0
     endif
     !goto 1000
  else
     !C          local minima rejected -------------------------------------------------------
     inputs_opt%inputPsiId=0  !ALEX says: Better do an input guess for the next escape
     if (iproc.eq.0) write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),l3,a)')  &
          hopp,escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' R '
     if (iproc.eq.0) write(*,'(a,1pe21.14)')' # rejected: ew-e>ediff ',e_wpos-e_pos

     !           write(200+iproc,'(2(1x,f10.0),1x,1pe24.17,2(1x,1pe10.3),3(1x,0pf5.2),l,a)')  &
     !                       hopp,escape,e_wpos-eref,ediff,ekinetic, &
     !                       escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' R '
     hopp_rej=hopp_rej+1.d0
     ediff=ediff*alpha2
     if (iproc.eq.0) call wtioput(ediff,ekinetic,dt)
     !goto 1000
     !C                          ------------------------------------------------------------
  endif

end do hopping_loop

!3000 continue

  if (iproc == 0) then
     write(67,*) 'writing final results'
     write(*,*) '# writing final results'
  end if
  call winter(iproc,atoms,re_pos,pos,npminx,nlminx,nbuf,nlmin,nlmin_l,accur, & 
       earr,elocmin,poslocmin,eref,ediff,ekinetic,dt)
  if (iproc == 0) then
     write(67,*) ' found ',nlmin_l,nlmin,' minima'
     write(*,*) '# found ',nlmin_l,nlmin,' minima'
     !this section could be put in wtioput. However, here we only put
     !these variables when the progam exited nicely. Otherwise,
     !the counters here will be reinitialized on restart.
     open(11,file='ioput',position='append')
     write(11,*)irestart      ,'irestart for aborted geometry optimization'
     write(11,*)av_ekinetic   ,'av_ekinetic'
     write(11,*)av_ediff      ,'av_ediff   '
     write(11,*)escape        ,'escape     '
     write(11,*)escape_sam    ,'escape_sam '
     write(11,*)escape_old    ,'escape_old '
     write(11,*)escape_new    ,'escape_new '
     write(11,*)hopp          ,'hopp       '
     write(11,*)hopp_acc      ,'hopp_acc   '
     write(11,*)hopp_rej      ,'hopp_rej   '
     write(11,*)egap          ,'egap       '
     write(11,*)esep          ,'esep       '
     write(11,*)count_sd      ,'count_sd   '
     write(11,*)count_cg      ,'count_cg   '
     write(11,*)count_md      ,'count_md   '
     write(11,*)nputback      ,'nputback   '
     close(11)
  endif


  call cpu_time(tcpu2)
  if (iproc.eq.0) then
     !C ratios from all the global counters
     write(67,'(i2,1x,a,3(1x,1pe10.3))') iproc,'ratio stuck,same,old,new', &
          escape_sam/escape,escape_old/escape,escape_new/escape
     write(67,'(i2,1x,a,2(1x,1pe10.3))') iproc,'ratio acc,rej', hopp_acc/hopp,hopp_rej/hopp
     write(67,'(i2,1x,a,3(1x,f12.1))') iproc,'count_md,count_sd,count_cg',count_md,count_sd,count_cg
     write(67,'(i2,1x,a,1x,1pe10.3)') iproc,'cpu(hrs) ', (tcpu2-tcpu1)/3600.d0
     write(67,'(i2,1x,a,2(1x,1pe10.3))') &
          iproc,'average ediff, ekinetic',av_ediff/hopp,av_ekinetic/escape
     write(*,'(a,1x,i8)') 'number of configurations for which atoms escaped ',nputback

     tt=0.d0
     ss=0.d0
     do i=1,nlmin
        tt=max(tt,earr(i,2))
        ss=ss+earr(i,2)
     enddo
     write(67,'(i2,a,f8.0)') iproc,'  most frequent visits ',tt
     write(67,'(i2,a,1pe10.3)') iproc,'  av. numb. visits per minimum',ss/nlmin
     write(67,'(a,e9.2)') 'minimum energy separation between presumably different configurations',egap
     if (escape_sam.gt.0) then
        esep=sqrt(esep/escape_sam)
        write(67,'(a,e9.2)') 'average energy separation between presumably identical configurations',esep
     endif
  endif

  !deallocations as in BigDFT
  i_all=-product(shape(atoms%lfrztyp))*kind(atoms%lfrztyp)
  deallocate(atoms%lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp',subname)
  i_all=-product(shape(atoms%iatype))*kind(atoms%iatype)
  deallocate(atoms%iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype',subname)
  i_all=-product(shape(atoms%natpol))*kind(atoms%natpol)
  deallocate(atoms%natpol,stat=i_stat)
  call memocc(i_stat,i_all,'natpol',subname)
  i_all=-product(shape(atoms%atomnames))*kind(atoms%atomnames)
  deallocate(atoms%atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames',subname)

  call free_restart_objects(rst,subname)


  ! deallocation of global's variables


  i_all=-product(shape(pos))*kind(pos)
  deallocate(pos,stat=i_stat)
  call memocc(i_stat,i_all,'pos',subname)
  i_all=-product(shape(earr))*kind(earr)
  deallocate(earr,stat=i_stat)
  call memocc(i_stat,i_all,'earr',subname)
  i_all=-product(shape(ff))*kind(ff)
  deallocate(ff,stat=i_stat)
  call memocc(i_stat,i_all,'ff',subname)
  i_all=-product(shape(wpos))*kind(wpos)
  deallocate(wpos,stat=i_stat)
  call memocc(i_stat,i_all,'wpos',subname)
  i_all=-product(shape(vxyz))*kind(vxyz)
  deallocate(vxyz,stat=i_stat)
  call memocc(i_stat,i_all,'vxyz',subname)
  i_all=-product(shape(gg))*kind(gg)
  deallocate(gg,stat=i_stat)
  call memocc(i_stat,i_all,'gg',subname)
  i_all=-product(shape(poslocmin))*kind(poslocmin)
  deallocate(poslocmin,stat=i_stat)
  call memocc(i_stat,i_all,'poslocmin',subname)

  if (iproc.eq.0) write(67,'(a,1x,3(1x,1pe10.3))') 'Out:ediff,ekinetic,dt',ediff,ekinetic,dt
  close(2) 

  !  call deallocate_wfd(wfd,'BigDFT')
  !  i_all=-product(shape(psi))*kind(psi)
  !  deallocate(psi,stat=i_stat)
  !  call memocc(i_stat,i_all,'psi','BigDFT')
  !  i_all=-product(shape(eval))*kind(eval)
  !  deallocate(eval,stat=i_stat)
  !  call memocc(i_stat,i_all,'eval','BigDFT')
  !  i_all=-product(shape(rxyz_old))*kind(rxyz_old)
  !  deallocate(rxyz_old,stat=i_stat)
  !  call memocc(i_stat,i_all,'rxyz_old','BigDFT')

  !finalize memory counting
  call memocc(0,0,'count','stop')
  close(111)
  if (nproc > 1) call MPI_FINALIZE(ierr)


contains

  subroutine mdescape(mdmin,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,rxyz, &
       nproc,iproc,atoms,rst,inputs_md)!  &
    ! Does a MD run with the atomic positiosn rxyz
    use module_base
    use module_types
    use module_interfaces
    implicit real*8 (a-h,o-z)
    type(atoms_data) :: atoms
    type(restart_objects) :: rst
    dimension ff(3,atoms%nat),gg(3,atoms%nat),vxyz(3,atoms%nat),rxyz(3,atoms%nat),rxyz_old(3,atoms%nat)
    type(input_variables) :: inputs_md
    character(len=4) :: fn
    !type(wavefunctions_descriptors), intent(inout) :: wfd
    !real(kind=8), pointer :: psi(:), eval(:)

    if(iproc==0) write(*,*) '# MINHOP start soften '

    !C initialize positions,velocities, forces
    call gausdist(atoms%nat,rxyz,vxyz)
    inputs_md%inputPsiId=1
    !if(iproc==0)write(*,*)' #  no softening'
    call soften(ekinetic,e_pos,ff,gg,vxyz,dt,count_md,rxyz, &
         nproc,iproc,atoms,rst,inputs_md)
    call velopt(atoms,rxyz,ekinetic,vxyz)
    call zero(3*atoms%nat,gg)

    if(iproc==0) write(*,*) '# MINHOP start MD'
    !C inner (escape) loop
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
    md_loop: do istep=1,1000

       !C      Evolution of the system according to 'VELOCITY VERLET' algorithm
       rkin=0.d0
       do iat=1,atoms%nat
          if (.not. atoms%lfrztyp(iat)) then
             if (atoms%geocode == 'P') then
                rxyz(1,iat)=modulo(rxyz(1,iat) + dt*vxyz(1,iat) + (.5d0*dt*dt)*gg(1,iat),&
                     atoms%alat1)
                rxyz(2,iat)=modulo(rxyz(2,iat) + dt*vxyz(2,iat) + (.5d0*dt*dt)*gg(2,iat),&
                     atoms%alat2)
                rxyz(3,iat)=modulo(rxyz(3,iat) + dt*vxyz(3,iat) + (.5d0*dt*dt)*gg(3,iat),&
                     atoms%alat3)
             else if (atoms%geocode == 'S') then
                rxyz(1,iat)=modulo(rxyz(1,iat) + dt*vxyz(1,iat) + (.5d0*dt*dt)*gg(1,iat),&
                     atoms%alat1)
                rxyz(2,iat)=       rxyz(2,iat) + dt*vxyz(2,iat) + (.5d0*dt*dt)*gg(2,iat)
                rxyz(3,iat)=modulo(rxyz(3,iat) + dt*vxyz(3,iat) + (.5d0*dt*dt)*gg(3,iat),&
                     atoms%alat3)
             else if (atoms%geocode == 'F') then
                rxyz(1,iat)=rxyz(1,iat) + dt*vxyz(1,iat) + (.5d0*dt*dt)*gg(1,iat)
                rxyz(2,iat)=rxyz(2,iat) + dt*vxyz(2,iat) + (.5d0*dt*dt)*gg(2,iat)
                rxyz(3,iat)=rxyz(3,iat) + dt*vxyz(3,iat) + (.5d0*dt*dt)*gg(3,iat)
             end if
             rkin=rkin+vxyz(1,iat)**2+vxyz(2,iat)**2+vxyz(3,iat)**2
          end if
       enddo
       rkin=rkin*.5d0

       enmin2=enmin1
       enmin1=en0000
       !    if (iproc.eq.0) write(*,*) 'CLUSTER FOR  MD'
       inputs_md%inputPsiId=1
       if (istep > 2) inputs_md%itermax=50
       call call_bigdft(nproc,iproc,atoms,rxyz,inputs_md,e_rxyz,ff,rst,infocode)

       !call wtmd(istep,atoms%nat,e_rxyz,rxyz,atoms%iatype,atoms%atomnames,atoms%natpol)
       if (iproc == 0) then
!          write(fn,'(i4.4)') istep+int(count_md)
          write(fn,'(i4.4)') istep
          call wtxyz('posmd_'//fn,e_rxyz,rxyz,atoms,'')
       end if

       en0000=e_rxyz-e_pos
       if (istep >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (istep >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
       econs_max=max(econs_max,rkin+e_rxyz)
       econs_min=min(econs_min,rkin+e_rxyz)
       devcon=econs_max-econs_min
       if (iproc.eq.0) write(17,'(a,i5,1x,1pe17.10,2(1x,i2))') 'MD ',&
            istep,e_rxyz,nummax,nummin
       if (iproc.eq.0) write(*,'(a,i5,1x,1pe17.10,2(1x,i2))') ' #MD ',&
            istep,e_rxyz,nummax,nummin
       if (nummin.ge.mdmin .and. istep > 50) then
          if (nummax.ne.nummin .and. iproc == 0) &
               write(67,*) 'WARNING: nummin,nummax',nummin,nummax
          exit md_loop
       endif
       do iat=1,atoms%nat
          at1=ff(1,iat)
          at2=ff(2,iat)
          at3=ff(3,iat)
          !C Evolution of the velocities of the system
          if (.not. atoms%lfrztyp(iat)) then
             vxyz(1,iat)=vxyz(1,iat) + (.5d0*dt) * (at1 + gg(1,iat))
             vxyz(2,iat)=vxyz(2,iat) + (.5d0*dt) * (at2 + gg(2,iat))
             vxyz(3,iat)=vxyz(3,iat) + (.5d0*dt) * (at3 + gg(3,iat))
          end if
          !C Memorization of old forces
          gg(1,iat) = at1
          gg(2,iat) = at2
          gg(3,iat) = at3
       end do
    end do md_loop
    if (istep >=1000) then
       if (iproc == 0) write(67,*) 'TOO MANY MD STEPS'
       dt=2.d0*dt
    end if
    !save the value of count_md for the moment
    count_md=count_md+real(istep,gp)

    !C MD stopped, now do relaxation

    !  if (iproc.eq.0) write(67,*) 'EXIT MD',istep
    
    ! adjust time step to meet precision criterion
    devcon=devcon/(3*atoms%nat-3)
    if (iproc == 0) &
         write(66,'(a,2(1x,1pe11.4),1x,i5)')&
         'MD devcon ',devcon,devcon/ekinetic,istep
    if (devcon/ekinetic.lt.7.d-2) then
       write(66,*) 'MD:old,new dt',dt,dt*1.05d0
       dt=dt*1.05d0
    else
       write(66,*) 'MD:old,new dt',dt,dt/1.05d0
       dt=dt*(1.d0/1.05d0)
    endif
    
  end subroutine mdescape
  
  

  subroutine soften(ekinetic,e_pos,fxyz,gg,vxyz,dt,count_md,rxyz, &
       nproc,iproc,atoms,rst,inputs_md)! &
    use module_base
    use module_types
    use module_interfaces
    implicit real*8 (a-h,o-z)
    type(atoms_data) :: atoms
    dimension fxyz(3*atoms%nat),gg(3*atoms%nat),vxyz(3*atoms%nat),rxyz(3*atoms%nat),rxyz_old(3*atoms%nat)
    type(input_variables) :: inputs_md
    type(restart_objects) :: rst
    !  type(wavefunctions_descriptors), intent(inout) :: wfd
    !  real(kind=8), pointer :: psi(:), eval(:)
    !Local variables
    dimension wpos(3*atoms%nat)

    nit=20
    eps_vxyz=5.d-1
    alpha=inputs_md%betax

    !allocate(wpos(3,nat),fxyz(3,nat))

    inputs_md%inputPsiId=1
    if(iproc==0)write(*,*)'# soften initial step '
    call call_bigdft(nproc,iproc,atoms,rxyz,inputs_md,etot0,fxyz,rst,infocode)

    ! scale velocity to generate dimer 

    svxyz=0.d0
    do i=1,3*atoms%nat
       iat=(i-1)/3+1
       if (.not. atoms%lfrztyp(iat)) then
          svxyz=svxyz+vxyz(i)**2
       end if
    enddo
    svxyz=eps_vxyz/sqrt(svxyz)
    do i=1,3*atoms%nat
       vxyz(i)=svxyz*vxyz(i)
    enddo
    !equivalent to
    !      vxyz = vxyz*eps_vxyz/(dnrm2(3*atoms%nat,vxyz,1))


    do it=1,nit
       !       if(iproc==0)write(*,*)'w   =      r       +       v  '
       !update the positions, if the atom is not frozen
!!$       do i=1,3*atoms%nat
!!$          iat=(i-1)/3+1
!!$          if (atoms%lfrztyp(iat)) then
!!$             wpos(i)=rxyz(i)
!!$          else
!!$             wpos(i)=rxyz(i)+vxyz(i)
!!$          end if
!!$          !if(iproc==0)write(*,*)wpos(i),rxyz(i),vxyz(i)
!!$       end do
       
       do iat=1,atoms%nat
          if (atoms%lfrztyp(iat)) then
             wpos(3*(iat-1)+1)=rxyz(3*(iat-1)+1)
             wpos(3*(iat-1)+2)=rxyz(3*(iat-1)+2)
             wpos(3*(iat-1)+3)=rxyz(3*(iat-1)+3)
          else
             if (atoms%geocode == 'P') then
                wpos(3*(iat-1)+1)=modulo(rxyz(3*(iat-1)+1)+vxyz(3*(iat-1)+1),atoms%alat1)
                wpos(3*(iat-1)+2)=modulo(rxyz(3*(iat-1)+2)+vxyz(3*(iat-1)+2),atoms%alat2)
                wpos(3*(iat-1)+3)=modulo(rxyz(3*(iat-1)+3)+vxyz(3*(iat-1)+3),atoms%alat3)
             else if (atoms%geocode == 'S') then
                wpos(3*(iat-1)+1)=modulo(rxyz(3*(iat-1)+1)+vxyz(3*(iat-1)+1),atoms%alat1)
                wpos(3*(iat-1)+2)=       rxyz(3*(iat-1)+2)+vxyz(3*(iat-1)+2)
                wpos(3*(iat-1)+3)=modulo(rxyz(3*(iat-1)+3)+vxyz(3*(iat-1)+3),atoms%alat3)
             else if (atoms%geocode == 'F') then
                wpos(3*(iat-1)+1)=rxyz(3*(iat-1)+1)+vxyz(3*(iat-1)+1)
                wpos(3*(iat-1)+2)=rxyz(3*(iat-1)+2)+vxyz(3*(iat-1)+2)
                wpos(3*(iat-1)+3)=rxyz(3*(iat-1)+3)+vxyz(3*(iat-1)+3)
             end if
          end if
          !if(iproc==0)write(*,*)wpos(i),rxyz(i),vxyz(i)
       end do

       !        wpos=rxyz+vxyz
       call call_bigdft(nproc,iproc,atoms,wpos,inputs_md,etot,fxyz,rst,infocode)
       fd2=2.d0*(etot-etot0)/eps_vxyz**2

       sdf=0.d0
       svxyz=0.d0
       do i=1,3*atoms%nat
          iat=(i-1)/3+1
          if (.not. atoms%lfrztyp(iat)) then
             sdf=sdf+vxyz(i)*fxyz(i)
             svxyz=svxyz+vxyz(i)*vxyz(i)
          end if
       end do
       !equivalent to
       !        sdf=ddot(3*atoms%nat,vxyz(1,1),1,fxyz(1,1),1)
       !        svxyz=dnrm2(3*atoms%nat,vxyz(1,1),1)

       curv=-sdf/svxyz
       if (it.eq.1) curv0=curv

       res=0.d0
       tt=0.d0
       do i=1,3*atoms%nat
          iat=(i-1)/3+1
          if (.not. atoms%lfrztyp(iat)) then
             tt=tt+fxyz(i)**2
             fxyz(i)=fxyz(i)+curv*vxyz(i)
             res=res+fxyz(i)**2
          end if
       end do
       res=sqrt(res)
       tt=sqrt(tt)
       !equivalent to
       !        tt=dsqrt(dnrm2(3*atoms%nat,fxyz(1,1),1))
       !        call daxpy(3*atoms%nat,curv,vxyz(1,1),1,fxyz(1,1),1)
       !        res=dsqrt(dnrm2(3*atoms%nat,fxyz(1,1),1))

       if(iproc==0)write(*,'(a,i3,5(f12.5))')'# soften it, curv, fd2,dE and res:',&
            it, curv, fd2,etot-etot0,res
       !if(iproc==0)write(11,'(i5,4(1pe13.5),1pe18.10)')&
       !if(iproc==0)write(11,*)&
       !        it,tt,res,curv,fd2,etot-etot0

       do iat=1,atoms%nat
          if (.not. atoms%lfrztyp(iat)) then
             if (atoms%geocode == 'P') then
                wpos(3*(iat-1)+1)=modulo(wpos(3*(iat-1)+1)+alpha*fxyz(3*(iat-1)+1),atoms%alat1)
                wpos(3*(iat-1)+2)=modulo(wpos(3*(iat-1)+2)+alpha*fxyz(3*(iat-1)+2),atoms%alat2)
                wpos(3*(iat-1)+3)=modulo(wpos(3*(iat-1)+3)+alpha*fxyz(3*(iat-1)+3),atoms%alat3)
             else if (atoms%geocode == 'S') then
                wpos(3*(iat-1)+1)=modulo(wpos(3*(iat-1)+1)+alpha*fxyz(3*(iat-1)+1),atoms%alat1)
                wpos(3*(iat-1)+2)=       wpos(3*(iat-1)+2)+alpha*fxyz(3*(iat-1)+2)
                wpos(3*(iat-1)+3)=modulo(wpos(3*(iat-1)+3)+alpha*fxyz(3*(iat-1)+3),atoms%alat3)
             else if (atoms%geocode == 'F') then
                wpos(3*(iat-1)+1)=wpos(3*(iat-1)+1)+alpha*fxyz(3*(iat-1)+1)
                wpos(3*(iat-1)+2)=wpos(3*(iat-1)+2)+alpha*fxyz(3*(iat-1)+2)
                wpos(3*(iat-1)+3)=wpos(3*(iat-1)+3)+alpha*fxyz(3*(iat-1)+3)
             end if

          end if
       end do
!!$       do i=1,3*atoms%nat
!!$          iat=(i-1)/3+1
!!$          if (.not. atoms%lfrztyp(iat)) then
!!$             wpos(i)=wpos(i)+alpha*fxyz(i)
!!$          end if
!!$       end do

       do i=1,3*atoms%nat
          vxyz(i)=wpos(i)-rxyz(i)
       end do
       !equivalent to
       !        call daxpy(3*atoms%nat,alpha,fxyz(1),1,wpos(1),1)
       !        vxyz=wpos-rxyz

       call  elim_moment(atoms%nat,vxyz)
       call  elim_torque(atoms%nat,rxyz,vxyz)

       svxyz=0.d0
       do i=1,3*atoms%nat
          iat=(i-1)/3+1
          if (.not. atoms%lfrztyp(iat)) then
             svxyz=svxyz+vxyz(i)*vxyz(i)
          end if
       end do
       !equivalent to
       !        svxyz=dnrm2(3*atoms%nat,vxyz(1),1)**2
       if (res <= curv*eps_vxyz*5.d-1) exit
       svxyz=eps_vxyz/dsqrt(svxyz)

       do i=1,3*atoms%nat
          vxyz(i)=vxyz(i)*svxyz
       end do
       !equivalent to
       !        call dscal(3*atoms%nat,svxyz,vxyz(1),1)

    end do ! iter
    !if(iproc==0)close(11)
    !if(res>curv*eps_vxyz*5.d-1)write(*,*)'# process', iproc,' has no convergence in low_cur_dir',res

    !        call  moment(nat,vxyz)
    !        call  torque(nat,rxyz,vxyz)
    !        deallocate(wpos,fxyz)
  end subroutine soften


end program
!!***




subroutine insert(nlminx,nbuf,nlmin,k_e_wpos,re_wpos,earr)
  ! inserts the energy re_wpos at position k_e_wpos and shifts up all other energies
  implicit real*8 (a-h,o-z)
  dimension earr(0:nlminx+nbuf,2)
  do k=nlmin-1,k_e_wpos+1,-1
     earr(k+1,1)=earr(k,1)
     earr(k+1,2)=earr(k,2)
  enddo
  earr(k_e_wpos+1,1)=re_wpos
  earr(k_e_wpos+1,2)=1.d0
  return
end subroutine insert


subroutine save_low_conf(nat,nlmin_l,npminx,e_wpos,pos,elocmin,poslocmin)
  !C save configuration if it is among the lowest ones in energy
  implicit real*8 (a-h,o-z)
  dimension elocmin(npminx)
  dimension pos(3,nat),poslocmin(3,nat,npminx)

  if (nlmin_l.le.npminx) then
     kmax=nlmin_l
     elocmin(kmax)=e_wpos
     do iat=1,nat
        poslocmin(1,iat,kmax)=pos(1,iat)
        poslocmin(2,iat,kmax)=pos(2,iat)
        poslocmin(3,iat,kmax)=pos(3,iat)
     enddo
  else
     ! find configuration kmax that is highest in energy
     emax=-1.d100
     do k=1,npminx
        if (elocmin(k).gt.emax) then
           emax=elocmin(k)
           kmax=k
        endif
     enddo
     if (e_wpos.lt.elocmin(kmax)) then
        elocmin(kmax)=e_wpos
        do iat=1,nat
           poslocmin(1,iat,kmax)=pos(1,iat)
           poslocmin(2,iat,kmax)=pos(2,iat)
           poslocmin(3,iat,kmax)=pos(3,iat)
        enddo
     endif
  endif


  return
end subroutine save_low_conf


subroutine hunt(xx,n,x,jlo)
  implicit none
  !C x is in interval [xx(jlo),xx(jlow+1)[ ; xx(0)=-Infinity ; xx(n+1) = Infinity
  !Arguments
  integer :: jlo,n
  real(kind=8) :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  if (n.le.0) stop 'hunt'
  if (n.eq.1) then
     if (x.ge.xx(1)) then
        jlo=1
     else
        jlo=0
     endif
     return
  endif
  ascnd=xx(n).ge.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    continue
     jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    continue
     jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 continue
  if(jhi-jlo.eq.1)then
     if(x.eq.xx(n))jlo=n
     if(x.eq.xx(1))jlo=1
     return
  endif
  jm=(jhi+jlo)/2
  if(x.ge.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt


subroutine velopt(at,rxyz,ekinetic,vxyz)
  !     assigns initial velocities for the MD escape part
  use module_base
  use module_types
  implicit none
  !implicit real*8 (a-h,o-z)
  real(gp), intent(in) :: ekinetic
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(inout) :: vxyz
  !local variables
  integer :: iat
  real(gp) :: rkin,rkinsum,sclvel


  !! Either random velocity distribution 
  !        call randdist(nat,rxyz,vxyz)
  !! or Gauss velocity distribution
  !! or exponential  velocity distribution
  !        call expdist(nat,rxyz,vxyz)
  !! or localized velocities
  !        call localdist(nat,rxyz,vxyz)
  !! or global move velocities
  !        to be implemented

  ! Soften previous velocity distribution
  !call soften(iproc,nat,rxyz,vxyz,rayl0,rayl,res,it)
  !^^^^^^^^^^^^^^^^^^^^^^
  !IMPORTANT: This is done at the level of mdescape in this version! 


  !C      Kinetic energy of the random velocities
  rkinsum= 0.d0      
  do iat=1,at%nat
     if (.not. at%lfrztyp(iat)) then
        rkinsum= rkinsum+vxyz(1,iat)**2+vxyz(2,iat)**2+vxyz(3,iat)**2
     end if
  end do
  rkin=.5d0*rkinsum/(3*at%nat-3)
  !       write(*,*) 'rkin,ekinetic',rkin,ekinetic

  !C      Rescaling of velocities to get reference kinetic energy
  sclvel= dsqrt(ekinetic/rkin)
  do iat=1,at%nat
     if (.not. at%lfrztyp(iat)) then
        vxyz(1,iat)=vxyz(1,iat)*sclvel
        vxyz(2,iat)=vxyz(2,iat)*sclvel
        vxyz(3,iat)=vxyz(3,iat)*sclvel
     end if
  end do

end subroutine velopt

subroutine randdist(nat,rxyz,vxyz)
  implicit real*8 (a-h,o-z)
  real tt
  dimension vxyz(3*nat),rxyz(3*nat)
  ! create a random displacement vector without translational and angular moment
  do i=1,3*nat
     call random_number(tt)
     vxyz(i)=dble(tt-.5)
  enddo
  call  elim_moment(nat,vxyz)
  call  elim_torque(nat,rxyz,vxyz)
  return
end subroutine randdist



subroutine gausdist(nat,rxyz,vxyz)
  !C  generates 3*nat random numbers distributed according to  exp(-.5*vxyz**2)
  implicit real*8 (a-h,o-z)
  real s1,s2
  !C On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
  parameter(eps=1.d-8)
  dimension vxyz(3*nat),rxyz(3*nat)

  do i=1,3*nat-1,2
     call random_number(s1)
     t1=eps+(1.d0-2.d0*eps)*dble(s1)
     call random_number(s2)
     t2=dble(s2)
     tt=sqrt(-2.d0*log(t1))
     vxyz(i)=tt*cos(6.28318530717958648d0*t2)
     vxyz(i+1)=tt*sin(6.28318530717958648d0*t2)
  enddo
  call random_number(s1)
  t1=eps+(1.d0-2.d0*eps)*dble(s1)
  call random_number(s2)
  t2=dble(s2)
  tt=sqrt(-2.d0*log(t1))
  vxyz(3*nat)=tt*cos(6.28318530717958648d0*t2)

  call elim_moment(nat,vxyz)
  call  elim_torque(nat,rxyz,vxyz)
  return
end subroutine gausdist

subroutine expdist(nat,rxyz,vxyz)
  !C  generates n random numbers distributed according to  exp(-x)
  implicit real*8 (a-h,o-z)
  real ss
  !C On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
  parameter(eps=1.d-8)
  dimension rxyz(3*nat),vxyz(3*nat)

  do i=1,3*nat
     call random_number(ss)
     tt=eps+(1.d0-2.d0*eps)*dble(ss)
     vxyz(i)=log(tt)
  enddo

  call elim_moment(nat,vxyz)
  call  elim_torque(nat,rxyz,vxyz)

  return
end subroutine expdist




subroutine localdist(nat,rxyz,vxyz)
  implicit real*8 (a-h,o-z)
  real*4 ts
  parameter(nbix=20)
  dimension rxyz(3,nat), vxyz(3,nat),nbi(nbix)
  parameter( bondlength=2.7d0)

  nloop=0
100 continue
  nloop=nloop+1
  if (nloop.gt.2) write(*,*) 'nloop=',nloop
  if (nloop.gt.11) stop 'ERROR LOCALDIST'

  ! pick an atom iat randomly
  call random_number(ts)
  iat=min(nat,int(ts*nat+1.))
  !       write(*,*) 'iat=',iat

  ! find iat's neighbors
  inb=0
  do i=1,nat
     dd=(rxyz(1,iat)-rxyz(1,i))**2+(rxyz(2,iat)-rxyz(2,i))**2+(rxyz(3,iat)-rxyz(3,i))**2
     if (dd < bondlength**2 .and. i /= iat) then
        inb=inb+1; if (inb.gt.nbix) stop 'enlarge size of nbi' ; nbi(inb)=i
     endif
  enddo
  !        write(*,*) 'inb=',inb
  if (inb < 2 ) goto 100

  ! pick another atom jat that is a neighbor of iat
  call random_number(ts)
  j=min(inb,int(ts*inb+1.))
  jat=nbi(j)
  !       write(*,*) 'jat=',jat


  ! Choose velocities for remaining atoms (i.e. not iat and jat )
  ampl=2.d-1
  do i=1,nat
     call random_number(ts) ; ts=ts-.5
     vxyz(1,i)=dble(ts)*ampl
     call random_number(ts) ; ts=ts-.5
     vxyz(2,i)=dble(ts)*ampl
     call random_number(ts) ; ts=ts-.5
     vxyz(3,i)=dble(ts)*ampl
  enddo
  ! Finally choose velocities for iat and jat 
  ampl=2.d0
  i=iat
  call random_number(ts) ; ts=ts-.5
  vxyz(1,i)=dble(ts)*ampl
  call random_number(ts) ; ts=ts-.5
  vxyz(2,i)=dble(ts)*ampl
  call random_number(ts) ; ts=ts-.5
  vxyz(3,i)=dble(ts)*ampl
  i=jat
  call random_number(ts) ; ts=ts-.5
  vxyz(1,i)=dble(ts)*ampl
  call random_number(ts) ; ts=ts-.5
  vxyz(2,i)=dble(ts)*ampl
  call random_number(ts) ; ts=ts-.5
  vxyz(3,i)=dble(ts)*ampl

  !        write(*,'(i3,3(1pe12.4))') iat,(rxyz(i,iat)+.5d0*vxyz(i,iat),i=1,3)
  !        write(*,'(i3,3(1pe12.4))') jat,(rxyz(i,jat)+.5d0*vxyz(i,jat),i=1,3)

  return
end subroutine localdist




subroutine torque(nat,rxyz,vxyz)
  implicit real*8 (a-h,o-z)
  dimension rxyz(3,nat),vxyz(3,nat)

  ! center of mass
  cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
  do iat=1,nat
     cmx=cmx+rxyz(1,iat)
     cmy=cmy+rxyz(2,iat)
     cmz=cmz+rxyz(3,iat)
  enddo
  cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat

  ! torque
  tx=0.d0 ; ty=0.d0 ; tz=0.d0
  do iat=1,nat
     tx=tx+(rxyz(2,iat)-cmy)*vxyz(3,iat)-(rxyz(3,iat)-cmz)*vxyz(2,iat)
     ty=ty+(rxyz(3,iat)-cmz)*vxyz(1,iat)-(rxyz(1,iat)-cmx)*vxyz(3,iat)
     tz=tz+(rxyz(1,iat)-cmx)*vxyz(2,iat)-(rxyz(2,iat)-cmy)*vxyz(1,iat)
  enddo
  write(*,'(a,3(1pe11.3))') 'torque',tx,ty,tz

  return
end subroutine torque


subroutine elim_torque(nat,rxyz,vxyz)
  implicit real*8 (a-h,o-z)
  dimension rxyz(3,nat),vxyz(3,nat),t(3)

  ! center of mass
  cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
  do iat=1,nat
     cmx=cmx+rxyz(1,iat)
     cmy=cmy+rxyz(2,iat)
     cmz=cmz+rxyz(3,iat)
  enddo
  cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat

  do it=1,100

     ! torque and radii in planes
     t(1)=0.d0 ; t(2)=0.d0 ; t(3)=0.d0
     sx=0.d0 ; sy=0.d0 ; sz=0.d0
     do iat=1,nat
        t(1)=t(1)+(rxyz(2,iat)-cmy)*vxyz(3,iat)-(rxyz(3,iat)-cmz)*vxyz(2,iat)
        t(2)=t(2)+(rxyz(3,iat)-cmz)*vxyz(1,iat)-(rxyz(1,iat)-cmx)*vxyz(3,iat)
        t(3)=t(3)+(rxyz(1,iat)-cmx)*vxyz(2,iat)-(rxyz(2,iat)-cmy)*vxyz(1,iat)
        sx=sx+(rxyz(1,iat)-cmx)**2
        sy=sy+(rxyz(2,iat)-cmy)**2
        sz=sz+(rxyz(3,iat)-cmz)**2
     enddo

     if (t(1)**2+t(2)**2+t(3)**2.lt.1.d-22) return

     ii=0
     tmax=0.d0
     do i=1,3
        if (t(i)**2.gt.tmax**2) then 
           ii=i
           tmax=t(i)
        endif
     enddo

     !         write(*,'(i4,3(1pe11.3))') ii,t


     ! modify velocities
     if (ii.eq.1) then 
        cx=t(1)/(sz+sy)
        do iat=1,nat
           vxyz(2,iat)=vxyz(2,iat)+cx*(rxyz(3,iat)-cmz)
           vxyz(3,iat)=vxyz(3,iat)-cx*(rxyz(2,iat)-cmy)
        enddo
     else if(ii.eq.2) then 
        cy=t(2)/(sz+sx)
        do iat=1,nat
           vxyz(1,iat)=vxyz(1,iat)-cy*(rxyz(3,iat)-cmz)
           vxyz(3,iat)=vxyz(3,iat)+cy*(rxyz(1,iat)-cmx)
        enddo
     else if(ii.eq.3) then 
        cz=t(3)/(sy+sx)
        do iat=1,nat
           vxyz(1,iat)=vxyz(1,iat)+cz*(rxyz(2,iat)-cmy)
           vxyz(2,iat)=vxyz(2,iat)-cz*(rxyz(1,iat)-cmx)
        enddo
     else
        stop 'wrong ii'
     endif

  enddo
  write(*,'(a,3(1pe11.3))') 'WARNING REMAINING TORQUE',t

  return
end subroutine elim_torque



subroutine moment(nat,vxyz)
  implicit real*8 (a-h,o-z)
  dimension vxyz(3,nat)

  sx=0.d0 ; sy=0.d0 ; sz=0.d0
  do iat=1,nat
     sx=sx+vxyz(1,iat)
     sy=sy+vxyz(2,iat)
     sz=sz+vxyz(3,iat)
  enddo
  write(*,'(a,3(1pe11.3))') 'momentum',sx,sy,sz

  return
end subroutine moment


subroutine elim_moment(nat,vxyz)
  implicit real*8 (a-h,o-z)
  dimension vxyz(3,nat)

  sx=0.d0 ; sy=0.d0 ; sz=0.d0
  do iat=1,nat
     sx=sx+vxyz(1,iat)
     sy=sy+vxyz(2,iat)
     sz=sz+vxyz(3,iat)
  enddo
  sx=sx/nat ; sy=sy/nat ; sz=sz/nat
  do iat=1,nat
     vxyz(1,iat)=vxyz(1,iat)-sx
     vxyz(2,iat)=vxyz(2,iat)-sy
     vxyz(3,iat)=vxyz(3,iat)-sz
  enddo

  return
end subroutine elim_moment




subroutine fix_fragmentation(iproc,at,rxyz,nputback)
  use module_base
  use module_types
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: at
  integer, intent(inout) :: nputback
  real(gp), dimension(3,at%nat) :: rxyz
  !local variables
  real(gp), parameter :: bondlength=8.0_gp
  integer :: iat,nloop,ncluster,ii,jat,jj,kat,nadd
  real(gp) :: xi,yi,zi,xj,yj,zj,ddmin,dd,d1,d2,d3,tt
  ! automatic arrays
  logical, dimension(at%nat) :: belong

  nloop=1

  fragment_loop: do

     iat=1
     belong(iat)=.true.
     ncluster=1
     do iat=2,at%nat
        belong(iat)=.false.
     enddo

     !   ic=0
     form_cluster: do
        nadd=0
        do iat=1,at%nat
           xi=rxyz(1,iat) 
           yi=rxyz(2,iat) 
           zi=rxyz(3,iat)
           if (belong(iat)) then 
              do jat=1,at%nat
                 xj=rxyz(1,jat) ; yj=rxyz(2,jat) ; zj=rxyz(3,jat)
                 if ( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 <= (bondlength*1.25d0)**2) then 
                    if (.not. belong(jat)) nadd=nadd+1
                    belong(jat)=.true. 
                 endif
              end do
           endif
        end do
        ncluster=ncluster+nadd
        !     ic=ic+1 ; write(*,*) 'nadd,ncluster',ic,nadd,ncluster
        if (nadd == 0) exit form_cluster
     enddo form_cluster

     if (ncluster == at%nat) then 
        !   write(*,*) 'No fragmentation has occured',nloop
        return

     else
        nputback=nputback+1

        if (iproc.eq.0) then
           write(*,*) 'fragmentation occured',nloop,ncluster
           write(444,*) at%nat,iat
           write(444,*) ' fragmented configuration ', nputback
           do kat=1,at%nat
              write(444,*) ' LJ  ',rxyz(1,kat),rxyz(2,kat),rxyz(3,kat)
           enddo
        endif

        ! make sure the part that flew away is smaller than the cluster
        if (ncluster <= at%nat/2) then
           !     write(*,*) 'FLIP'
           do iat=1,at%nat
              belong(iat)=.not. belong(iat)
           enddo
        endif

        ! pull back the fragment of atoms that flew away
        ii=-99999
        do iat=1,at%nat
           if (.not. belong(iat)) then
              xi=rxyz(1,iat) 
              yi=rxyz(2,iat) 
              zi=rxyz(3,iat)
              ddmin=1.e100_gp
              jj=-99999
              do jat=1,at%nat
                 if (belong(jat)) then
                    xj=rxyz(1,jat) 
                    yj=rxyz(2,jat) 
                    zj=rxyz(3,jat)
                    dd= (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 
                    if (dd < ddmin) then 
                       jj=jat
                       ii=iat
                       ddmin=dd
                    endif
                 endif
              enddo
           endif
        enddo

        d1=rxyz(1,ii)-rxyz(1,jj)
        d2=rxyz(2,ii)-rxyz(2,jj)
        d3=rxyz(3,ii)-rxyz(3,jj)
        tt=bondlength/sqrt(d1**2+d2**2+d3**2)
        do iat=1,at%nat
           if (.not. belong(iat) .and. .not. at%lfrztyp(iat)) then
              if (at%geocode == 'P') then
                 rxyz(1,iat)=modulo(rxyz(1,iat)+d1*(tt-1.0d0),at%alat1)
                 rxyz(2,iat)=modulo(rxyz(2,iat)+d2*(tt-1.0d0),at%alat2)
                 rxyz(3,iat)=modulo(rxyz(3,iat)+d3*(tt-1.0d0),at%alat3)
              else if (at%geocode == 'S') then
                 rxyz(1,iat)=modulo(rxyz(1,iat)+d1*(tt-1.0d0),at%alat1)
                 rxyz(2,iat)=       rxyz(2,iat)+d2*(tt-1.0d0)
                 rxyz(3,iat)=modulo(rxyz(3,iat)+d3*(tt-1.0d0),at%alat3)
              else
                 rxyz(1,iat)=rxyz(1,iat)+d1*(tt-1.0d0)
                 rxyz(2,iat)=rxyz(2,iat)+d2*(tt-1.0d0)
                 rxyz(3,iat)=rxyz(3,iat)+d3*(tt-1.0d0)
              end if
           endif
        enddo

        if (iproc.eq.0) then
           write(444,*) at%nat, 'atomic ' 
           write(444,*) ' fixed configuration ', nputback
           do iat=1,at%nat
              write(444,*) ' LJ  ',rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
           enddo
        endif
        nloop=nloop+1
     endif
  end do fragment_loop

end subroutine fix_fragmentation


subroutine winter(iproc,at,re_pos,pos,npminx,nlminx,nbuf,nlmin,nlmin_l,accur, & 
     earr,elocmin,poslocmin,eref,ediff,ekinetic,dt)
  use module_base
  use module_types
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: npminx,nlminx,nbuf,nlmin,nlmin_l,iproc
  real(gp), intent(in) :: re_pos,eref,ediff,ekinetic,dt,accur
  type(atoms_data), intent(in) :: at
  real(gp), dimension(npminx), intent(in) :: elocmin
  real(gp), dimension(3,at%nat), intent(in) :: pos
  real(gp), dimension(0:nlminx+nbuf,2), intent(in) :: earr
  real(gp), dimension(3,at%nat,npminx), intent(in) :: poslocmin
  !local variables
  character(len=5) :: fn
  character(len=20) :: filename
  integer :: mm,k

  if (iproc == 0) then
     call wtxyz('poscur',re_pos,pos,at,'')
     write(*,*) ' wrote poscur.xyz for  RESTART'
  end if
     
  call wtpos(iproc,at,npminx,nlminx,nbuf,nlmin,nlmin_l,poslocmin,earr,elocmin)

  if (iproc == 0) then
     write(*,*) ' wrote poslow files'
     
     open(unit=12,file='earr.dat',status='unknown')
     mm=min(nlmin,nlminx+nbuf)
     write(12,'(2(i10),a)') mm,mm+5,&
          ' # of minima already found, # of minima to be found in consecutive run'
     write(12,'(e24.17,1x,a)') eref,'   eref'
     write(12,'(e24.17,1x,a)') accur,'   accur'
     do k=1,mm
        write(12,'(e24.17,1x,1pe17.10)') earr(k,1),earr(k,2)
     enddo
     write(*,*) ' wrote earr.dat for  RESTART'
     close(12)
     
     call  wtioput(ediff,ekinetic,dt)
     write(*,*) ' wrote ioput for  RESTART'
  end if

end subroutine winter


subroutine wtioput(ediff,ekinetic,dt)
  implicit real*8 (a-h,o-z)
  open(unit=11,file='ioput',status='unknown')
  write(11,'(3(1x,1pe24.17),a)') ediff,ekinetic,dt,' ediff, ekinetic and dt'
  close(11)
end subroutine wtioput


subroutine wtbest(nat,energy,pos,iatype,atomnames,natpol)
  implicit real*8 (a-h,o-z)
  character(41) filename
  character(20) atomnames
  character(3) fn
  dimension pos(3,nat),iatype(nat),atomnames(100)
  integer, dimension(nat):: natpol

  !C generate filename and open files
  filename = 'posbest.xyz'
  open(unit=9,file=filename,status='unknown')
  write(9,'(i4,2x,a,1pe24.17)') nat,'  atomic', energy
  write(9,'(e24.17)') energy
  do iat=1,nat
     write(9,'(a8,3x,3(1x,1pe24.17),3x,i5.5)') atomnames(iatype(iat)),&
          pos(1,iat),pos(2,iat),pos(3,iat),natpol(iat)-100
  enddo
  close(9)
  return
end subroutine wtbest

subroutine wtmd(istep,nat,energy,pos,iatype,atomnames,natpol)
  implicit real*8 (a-h,o-z)
  character(20) filename
  character(20) atomnames(100)
  character(4) fn
  dimension pos(3,nat),iatype(nat)
  integer, dimension(nat):: natpol

  !C generate filename and open files
  write(fn,'(i4.4)') istep
  filename = 'posmd_'//fn//'.xyz'
  open(unit=9,file=filename,status='unknown')
  write(9,'(i4,2x,a,1x,1pe24.17)') nat,'  atomic',energy
  write(9,'(e24.17)') energy
  do iat=1,nat
     write(9,'(a8,3x,3(1x,1pe24.17),3x,i5.5)') atomnames(iatype(iat)),&
          pos(1,iat),pos(2,iat),pos(3,iat),natpol(iat)-100
  enddo
  close(9)
  return
end subroutine wtmd



subroutine wtlmin(nconjgr,nat,energy,fnrm,pos,iatype,atomnames,natpol)
  implicit real*8 (a-h,o-z)
  character(20) filename
  character(20) atomnames(100)
  character(4) fn
  dimension pos(3,nat),iatype(nat)
  integer, dimension(nat):: natpol

  !C generate filename and open files
  write(fn,'(i4.4)') nconjgr
  filename = 'poslocm_'//fn//'.xyz'
  open(unit=9,file=filename,status='unknown')
  write(9,'(i4,2x,a,1x,1pe24.17)') nat,'  atomic',energy
  write(9,'(e24.17,1x,1pe10.3)') energy,fnrm
  do iat=1,nat
     write(9,'(a8,3x,3(1x,1pe24.17),3x,i5.5)') atomnames(iatype(iat)),&
          pos(1,iat),pos(2,iat),pos(3,iat),natpol(iat)-100
  enddo
  close(9)

end subroutine wtlmin



subroutine wtpos(iproc,at,npminx,nlminx,nbuf,nlmin,nlmin_l,pos,earr,elocmin)
  use module_base
  use module_types
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: npminx,nlminx,nbuf,nlmin,nlmin_l,iproc
  type(atoms_data), intent(in) :: at
  real(gp), dimension(npminx), intent(in) :: elocmin
  real(gp), dimension(0:nlminx+nbuf,2), intent(in) :: earr
  real(gp), dimension(3,at%nat,npminx), intent(in) :: pos
  !local variables
  character(len=5) :: fn
  character(len=17) :: filename
  character(len=20) :: atomnames
  integer :: k,kk

  !	do i=1,min(40,nlmin,nlminx+nbuf)
  !	write(*,*) i,earr(i,1)
  !        enddo

  do k=1,min(nlmin_l,npminx)
     !        write(*,*) 'k,elocmin(k)',k,elocmin(k)


     !C Classify the configuration in the global ranking
     kk=0
     find_kk : do
        kk=kk+1
        if (kk > min(nlmin,nlminx+nbuf)) then 
           if (iproc == 0) write(*,*) 'ranking error for',k
           stop 
        endif
        if (earr(kk,1) == elocmin(k)) exit find_kk
     end do find_kk

     if (kk <= npminx) then

        !        write(*,'(a,i2,i4,i4,1x,1pe21.14)') 'k,kk,elocmin(k)',k,kk,elocmin(k)

        !C generate filename and open files
        if (iproc == 0) then
           write(fn,'(i5.5)') kk
           call wtxyz('poslow'//fn,elocmin(k),pos(1,1,k),at,'')
        end if
     endif

  end do

end subroutine wtpos



subroutine zero(n,x)
  implicit real*8 (a-h,o-z)
  dimension x(n)
  do j=1,n
     x(j)=0.d0
  end do
  return
end subroutine zero


real*8 function round(enerd,accur)
  implicit none
  real*8 enerd,accur
  integer*8 ii
  ii=enerd/accur
  round=ii*accur
  !           write(*,'(a,1pe24.17,1x,i17,1x,1pe24.17)') 'enerd,ii,round',enerd,ii,round
  return
end function round



subroutine rdposout(igeostep,rxyz,nat)
  implicit none
  integer, intent(in) :: igeostep,nat
  real(kind=8), dimension(3,nat), intent(out) :: rxyz
  !local variables
  character(len=3) :: fn
  character(len=20) :: filename
  integer :: iat
  write(fn,'(i3.3)') igeostep
  !filename = 'posout_'//fn//'.ascii'
  filename = 'posout_'//fn//'.xyz'
  ! write(*,*)'# reading unrelaxed structure from ',filename
  open(unit=9,file=filename,status='old')
  read(9,*)fn!no need for header
  read(9,*)fn!same
  do iat=1,nat
     read(9,*)fn,rxyz(:,iat)!we know the atom types already
  enddo
  close(unit=9)
end subroutine rdposout

!routine for adjusting the dimensions with the center of mass in the middle
subroutine adjustrxyz(nat,alat1,alat2,alat3,rxyz)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp) ,intent(in) :: alat1,alat2,alat3
  real(gp), dimension(3,nat), intent(inout) :: rxyz
  !local variables
  integer :: iat,i 
  real(gp), dimension(3)  :: cent

  do i=1,3
     cent(i)=0.0_gp
  enddo
  do iat=1,nat
     do i=1,3
        cent(i)=cent(i)+rxyz(i,iat)
     enddo
  enddo
  do i=1,3
     cent(i)=cent(i)/real(nat,gp)
  enddo

  write(*,'(a,6(1x,e9.2))') 'old CM, shift',(cent(i),i=1,3),  & 
       -cent(1)+alat1*.5_gp,-cent(2)+alat2*.5_gp,-cent(3)+alat3*.5_gp

  do iat=1,nat
     rxyz(1,iat)=rxyz(1,iat)-cent(1)+alat1*.5_gp
     rxyz(2,iat)=rxyz(2,iat)-cent(2)+alat2*.5_gp
     rxyz(3,iat)=rxyz(3,iat)-cent(3)+alat3*.5_gp
  enddo
end subroutine adjustrxyz

