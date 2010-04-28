!!****p* MINHOP
!!
!! DESCRIPTION
!!  Main program fro the minima hopping
!!  New modified version 17th Nov 2009 Sandip De  
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
  use ab6_symmetry

  implicit real(kind=8) (a-h,o-z)
  real(kind=4) :: tts
  logical :: newmin,CPUcheck
  character(len=20) :: unitsp,units,atmn
  character(len=80) :: line
  type(atoms_data) :: atoms
  type(input_variables) :: inputs_opt, inputs_md
  type(restart_objects) :: rst
  character(len=20), dimension(:), allocatable :: atomnames
  integer, parameter :: npminx=200
  !C parameters for minima hopping
  integer, parameter :: mdmin=2
  real(kind=8), parameter :: beta1=1.10d0,beta2=1.10d0,beta3=1.d0/1.10d0
  real(kind=8), parameter :: alpha1=1.d0/1.10d0,alpha2=1.10d0
  real(kind=8) :: elocmin(npminx)
  real(kind=8), allocatable, dimension(:,:) ::ff,wpos,vxyz,gg,earr
  real(kind=8),allocatable, dimension(:,:,:):: poslocmin
  real(kind=8), dimension(:,:), pointer :: pos
  integer :: iproc,nproc,iat,ityp,j,i_stat,i_all,ierr,infocode
  character(len=*), parameter :: subname='global'
  character(len=41) :: filename
  character(len=4) :: fn4
  character(len=5) :: fn5
  character(len=50) :: comment
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem

  ! Start MPI version
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  !call system('echo $HOSTNAME')

  ! Read performance inputs if present
  call perf_input_variables(iproc,'input.perf',inputs_opt)

  ! Initialize memory counting
  !call memocc(0,iproc,'count','start')

  if (iproc == 0)then
     write(*,'(23x,a)')' NEW '
     write(*,'(23x,a)')'      __  __ _ _  _ _   _  __  ___ '
     write(*,'(23x,a)')'     |  \/  |_| \| | |_| |/  \| _ \ '
     write(*,'(23x,a)')'     | |\/| |-|    |  _  | <> |  _/ '
     write(*,'(23x,a)')'     |_|  |_|_|_|\_|_| |_|\__/|_|     WITH'
     write(*,'(23x,a)')''
     write(*,'(23x,a)')''
     write(*,'(23x,a)')''
     call print_logo()
     write(*,'(23x,a)')'----> you can grep this file for # to compare with global.out'
     write(*,'(23x,a)')' #NOTE: this version reads nspin, mpol from input.dat'
  end if

  open(unit=67,file='global.out')


  if (iproc == 0) write(67,'(a,3(1x,1pe11.4))') 'beta1,beta2,beta3',beta1,beta2,beta3
  if (iproc == 0) write(67,'(a,2(1x,1pe11.4))') 'alpha1,alpha2',alpha1,alpha2
  ratio=-log(alpha2)/log(alpha1)
  if (iproc == 0) write(67,'(a,2(1x,1pe10.3))') 'predicted fraction accepted, rejected', & 
       ratio/(1.d0+ratio), 1.d0/(1.d0+ratio)
  tt= beta2**ratio*beta3
  if (iproc == 0) write(67,'(a,2(1x,1pe11.4))') 'critical ratio (.ge. 1)',tt
  if (tt.lt.0.999999999999999d0) stop 'incompatible alpha s, beta s'
  if (iproc == 0) write(67,*) 'mdmin',mdmin


  call cpu_time(tcpu1)
  ! read  earr.dat
  open(unit=12,file='earr.dat',status='unknown')
  read(12,*) nlmin,nlminx
  read(12,*) eref
  if (iproc == 0) write(67,*) 'eref=',eref
  if (iproc == 0) write(7,*) 'eref=',eref
  read(12,*) accur
  if (iproc == 0) write(67,*) 'accuracy for rounding=',accur
  if (nlmin.gt.nlminx) stop 'nlmin>nlminx'
  allocate(earr(0:nlminx,2+ndebug),stat=i_stat)
  call memocc(i_stat,earr,'earr',subname)
  earr(0,1)=-1.d100
  if (nlmin == 0) then 
     if (iproc == 0) write(67,*) 'New run with nlminx=',nlminx
     if (iproc == 0) write(*,*) '#New run with nlminx=',nlminx
  else
     if (iproc == 0) write(67,*) 'Restart run with nlmin, nlminx=',nlmin,nlminx
     if (iproc == 0) write(*,*) '#Restart run with nlmin, nlminx=',nlmin,nlminx
     do k=1,nlmin
        read(12,*) earr(k,1),earr(k,2)
        if (earr(k,1).lt.earr(k-1,1)) stop 'wrong ordering in earr.dat'
     enddo
     if (iproc == 0) write(67,*) ' read earr.dat'
     if (iproc == 0) write(*,*) '# read earr.dat'
  endif
  close(12)

  call read_atomic_file('poscur',iproc,atoms,pos)
  !Read input parameters for geometry optimization 
  call default_input_variables(inputs_opt)
  call dft_input_variables(iproc,'input.dft',inputs_opt)
  call geopt_input_variables('input.geopt',inputs_opt)
  call kpt_input_variables(iproc,'input.kpt',inputs_opt,atoms)

  !read input parameters for molecular dynamics
  call default_input_variables(inputs_md)
  call dft_input_variables(iproc,'mdinput.dft',inputs_md)
  call geopt_input_variables('mdinput.geopt',inputs_md)
  call kpt_input_variables(iproc,'input.kpt',inputs_md,atoms)


  do iat=1,atoms%nat
     if (atoms%ifrztyp(iat) == 0) then
        call random_number(tt)
        pos(1,iat)=pos(1,iat)+inputs_opt%randdis*tt
        call random_number(tt)
        pos(2,iat)=pos(2,iat)+inputs_opt%randdis*tt
        call random_number(tt)
        pos(3,iat)=pos(3,iat)+inputs_opt%randdis*tt
     end if
  enddo
  
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
 
! read random offset
  open(unit=11,file='rand.inp')
  read(11,*) nrandoff
  !        write(*,*) 'nrandoff ',nrandoff
  close(11)
  do i=1,nrandoff
     call random_number(ts)
  enddo
  
  ! read input parameters
  write(filename,'(a6,i3.3)') 'ioput'   !,iproc
  open(unit=11,file='ioput',status='old')
  read(11,*) ediff,ekinetic,dt,nsoften
  close(11)
  !write(*,'(a,1x,i3,3(1x,e10.3),1x,i4)') 'In :iproc,ediff,ekinetic,dt,nsoften',iproc,ediff,ekinetic,dt,nsoften
  if (iproc == 0) write(*,'(a,1x,3(1x,e10.3),1x,i4)') 'In :ediff,ekinetic,dt,nsoften',ediff,ekinetic,dt,nsoften

  if (iproc == 0) open(unit=16,file='geopt.mon',status='unknown')

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
        read(9,*) natp,unitsp,elocmin(nlmin_l)
        if (atoms%nat.ne.natp) stop   'nat <> natp'
        if (trim(unitsp).ne.trim(atoms%units) .and. iproc.eq.0) write(*,*)  & 
                 '# different units in poslow and poscur file: ',trim(unitsp),' ',trim(atoms%units)
        read(9,*)
        do iat=1,atoms%nat
          read(9,*) atmn,t1,t2,t3
          if (atoms%units=='angstroem' .or. atoms%units=='angstroemd0') then ! if Angstroem convert to Bohr
              poslocmin(1,iat,nlmin_l)=t1/bohr 
              poslocmin(2,iat,nlmin_l)=t2/bohr 
              poslocmin(3,iat,nlmin_l)=t3/bohr
          else
              poslocmin(1,iat,nlmin_l)=t1
              poslocmin(2,iat,nlmin_l)=t2
              poslocmin(3,iat,nlmin_l)=t3
          endif
        enddo
        close(9)
        if (iproc == 0) write(67,*) 'read file',filename
        if (iproc == 0) write(*,*) '# read file',filename
        kk=kk+1
     end do
     if (iproc == 0) write(67,*) 'read ',nlmin_l,'poslow files'
     if (iproc == 0) write(*,*) '# read ',nlmin_l,'poslow files'
  endif


! open output files
         if (iproc==0) open(unit=2,file='global.mon',status='unknown',position='append')

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

!C first local minimum
  count_sd=0.d0
  count_cg=0.d0
  count_soft=0.d0
  count_md=0.d0
  nputback=0

  do iat=1,atoms%nat
     if (atoms%geocode == 'P') then
        pos(1,iat)=modulo(pos(1,iat),atoms%alat1)
        pos(2,iat)=modulo(pos(2,iat),atoms%alat2)
        pos(3,iat)=modulo(pos(3,iat),atoms%alat3)
     else if (atoms%geocode == 'S') then
        pos(1,iat)=modulo(pos(1,iat),atoms%alat1)
        pos(3,iat)=modulo(pos(3,iat),atoms%alat3)
     end if
  end do


  inputs_opt%inputPsiId=0
  call init_restart_objects(iproc,inputs_opt%iacceleration,atoms,rst,subname)
  call call_bigdft(nproc,iproc,atoms,pos,inputs_md,e_pos,ff,rst,infocode)

  write(17,*) 'ENERGY ',e_pos
  energyold=1.d100
  ncount_bigdft=0

  if (iproc == 0) write(*,*)'# calling conjgrad for the first time here. energy ',e_pos

!  if (atoms%geocode == 'P') & 
!       call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,pos)

  call geopt(nproc,iproc,pos,atoms,ff,e_pos,rst,inputs_md,ncount_bigdft)
  if (iproc == 0) then
     write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
     write(*,*) '# ', ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
  end if

  call geopt(nproc,iproc,pos,atoms,ff,e_pos,rst,inputs_opt,ncount_bigdft)
  if (iproc == 0) then
     write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for accurate initial conf'
     write(*,*) '# ', ncount_bigdft,' Wvfnctn Opt. steps for accurate initial conf'
  end if

  nconjgr=0
  if (iproc == 0) then 
     tt=dnrm2(3*atoms%nat,ff,1)
     write(fn4,'(i4.4)') nconjgr
     write(comment,'(a,1pe10.3)')'fnrm= ',tt
     call write_atomic_file('poslocm_'//fn4,e_pos-eref,pos,atoms,trim(comment))
  endif
  nconjgr=nconjgr+1

  re_pos=round(e_pos-eref,accur)
  if (iproc == 0) then
     write(67,'(a,1x,i3,3(1x,1pe17.10))') 'INPUT(relaxed): iproc, e_pos,re_pos,eref ',iproc,e_pos,re_pos,eref
     write(*,'(a,1x,i3,3(1x,1pe17.10))') ' # INPUT(relaxed): iproc, e_pos,re_pos,eref ',iproc,e_pos,re_pos,eref
  end if


  if (nlmin.gt.0) then
     if (iproc == 0) write(67,'(a,2(1x,1pe24.17))') 'new/old energy for input file',re_pos
     if (iproc == 0) write(*,'(a,2(1x,1pe24.17))') '# new/old energy for input file',re_pos
  endif
  k_e_wpos=1
  if (nlmin == 0) then
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
  if (iproc == 0) then
     write(67,*) 'iproc,initial re_sm',iproc,re_sm
     write(*,*) '# iproc,initial re_sm',iproc,re_sm
     write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3))') hopp,escape,e_pos-eref,ediff,ekinetic
  end if

  nlmin_old=nlmin
  CPUcheck=.false.

  !C outer (hopping) loop
!  hopping_loop: do
1000 continue
     if (nlmin >= nlminx) then 
        write(67,*)iproc,'has  nlminx collected',nlmin
        write(*,*)'#:process', iproc,'has  nlminx collected',nlmin
        !             do i=1,nlmin ;  write(*,*) earr(i,1) ; enddo
        goto 3000
     endif
     !            Energy has reached taregt eref and global minimum is presumably found
     if (re_sm <= 1.d-3) then
        write(*,*)'# process', iproc,'success: relative energy < 0.001'
        goto 3000
     endif

5555 continue

!C check whether CPU time exceeded
     tleft=1.d100
        if(iproc==0 .and. CPUcheck)then
        open(unit=55,file='CPUlimit_global',status='unknown')
        read(55,*,end=555) cpulimit ; cpulimit=cpulimit*3600
        write(*,'(a,i5,i3,2(1x,e9.2))') 'iproc,nlmin,tcpu2-tcpu1,cpulimit',iproc,nlmin,tcpu2-tcpu1,cpulimit
        close(55)
        call cpu_time(tcpu2)
        tleft=cpulimit-(tcpu2-tcpu1)
     end if
555    continue
       call MPI_BCAST(tleft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       if (tleft < 0.d0) then
       write(*,*) 'CPU time exceeded',tleft
       goto 3000
       endif
          CPUcheck=.true.

     do iat=1,atoms%nat
        wpos(1,iat)=pos(1,iat)
        wpos(2,iat)=pos(2,iat) 
        wpos(3,iat)=pos(3,iat)
     enddo

     escape=escape+1.d0
     call mdescape(nsoften,mdmin,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,wpos, &
          nproc,iproc,atoms,rst,inputs_md)

     call fix_fragmentation(iproc,atoms,wpos,nputback)

     av_ekinetic=av_ekinetic+ekinetic
     ncount_bigdft=0

!5556 continue ! entry point for restart of optimization at cluster step irestart+1
!     if (atoms%geocode == 'P') & 
!        call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,wpos)
     call geopt(nproc,iproc,wpos,atoms,ff,e_wpos,rst,inputs_md,ncount_bigdft)

     if (iproc == 0) write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
     if (iproc == 0) write(*,*)'# ', ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
     !ncount_bigdft=0
     !ncount_cluster=0
!     if (atoms%geocode == 'P') & 
!          call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,wpos)
      call geopt(nproc,iproc,wpos,atoms,ff,e_wpos,rst,inputs_opt,ncount_bigdft)

     if (iproc == 0) write(67,*) ncount_bigdft,' Wvfnctn Opt. steps for accurate geo. rel of MD conf.'
     if (iproc == 0) write(*,*)'# ', ncount_bigdft,' Wvfnctn Opt. steps for accurate geo. rel of MD conf.'
     if (iproc == 0) then 
        tt=dnrm2(3*atoms%nat,ff,1)
        !call wtlmin(nconjgr,atoms%nat,e_wpos-eref,tt,wpos,atoms%iatype,atoms%atomnames,atoms%natpol)
        write(fn4,'(i4.4)') nconjgr
        write(comment,'(a,1pe10.3)')'fnrm= ',tt
  call write_atomic_file('poslocm_'//fn4,e_wpos-eref,wpos,atoms,trim(comment))
     endif
  nconjgr=nconjgr+1
  re_wpos=round(e_wpos-eref,accur)
  if (iproc == 0) write(67,'(a,i3,i3,4(1x,1pe14.7))')  & 
       'nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos', nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos
  if (iproc == 0) write(*,'(a,i3,i3,4(1x,1pe14.7))')  & 
       '# nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos', nlmin_l,nlmin,e_wpos,e_pos,re_wpos,re_pos
  !C not escaped
  if (re_pos == re_wpos) then
     escape_sam=escape_sam+1.d0
     esep=esep+(e_pos-e_wpos)**2
     ekinetic=ekinetic*beta1
     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)
     if (iproc == 0) write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),a)')  &
          hopp,escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     if (iproc == 0) write(*,'(a)')' # no escape from current minimum.'
     goto 5555
  endif



  !C continue since escaped
  !C  check whether new minimum
  call hunt(earr(1,1),min(nlmin,nlminx),re_wpos,k_e_wpos)
  if (re_wpos == earr(k_e_wpos,1)) then
     if (iproc == 0) write(67,'(a,i3,i3,i4,1x,1pe14.7)')  & 
          'nlmin_l,nlmin,k_e_wpos,re_wpos=earr',nlmin_l,nlmin,k_e_wpos,re_wpos
     if (iproc == 0) write(*,'(a,i3,i3,i4,1x,1pe14.7)')  & 
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
     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)

     if (k_e_wpos == 0) then
        re_sm=min(re_sm,re_wpos)
        if (iproc == 0) then 
           write(67,'(a,i7,1x,1pe24.17,1x,i2)') 'new lowest ',nlmin,re_wpos,iproc
           write(*,'(a,i7,1x,1pe24.17,1x,i2)') '# new lowest ',nlmin,re_wpos,iproc
           call write_atomic_file('posbest',re_wpos,wpos,atoms,'')
        endif
     else
        if (iproc == 0) write(*,'(a,i3)')' # New local minimum, lower are ',k_e_wpos 
     endif
  endif

  hopp=hopp+1.d0

  !C master: Monte Carlo step for local minima hopping
  !write(67,*) 'e_wpos,e_pos',e_wpos,e_pos
  if(e_wpos-e_pos.lt.ediff) then
     !C          local minima accepted -------------------------------------------------------
     hopp_acc=hopp_acc+1.d0
     ediff=ediff*alpha1
     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)
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
        call insert(nlminx,nlmin,k_e_wpos,re_wpos,earr(0,1))
        !                write(67,*) iproc,'EARR'
        !                do i=1,nlmin ;  write(*,*) earr(i,1),earr(i,2) ; enddo
        !C            save configuration if it is among the lowest ones in energy
        call save_low_conf(atoms%nat,nlmin_l,npminx,re_wpos,pos,elocmin,poslocmin)
        k_e_wpos=k_e_wpos+1   ! update for the writing statement below
     else
        !C            old minimum revisited
        earr(k_e_wpos,2)=earr(k_e_wpos,2)+1.d0
     endif
      if (iproc == 0) then
     write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),l3,a,i5)')  &
          hopp,escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A ', int(earr(k_e_wpos,2))
! write intermediate results
      write(*,*) 'WINTER'
      call winter(atoms,re_pos,pos,npminx,nlminx,nlmin,nlmin_l,accur, & 
           earr,elocmin,poslocmin,eref,ediff,ekinetic,dt,nsoften)
      end if
      goto 1000
  else
     !C          local minima rejected -------------------------------------------------------
     inputs_opt%inputPsiId=0  !ALEX says: Better do an input guess for the next escape
     if (iproc == 0) write(2,'(2(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),l3,a,i5)')  &
          hopp,escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' R ', int(earr(k_e_wpos,2))
     if (iproc == 0) write(*,'(a,1pe21.14)')' # rejected: ew-e>ediff ',e_wpos-e_pos

     hopp_rej=hopp_rej+1.d0
     ediff=ediff*alpha2
     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)
     goto 1000
     !C                          ------------------------------------------------------------
  endif

!end do hopping_loop

3000 continue

  if (iproc == 0) then
     write(67,*) 'writing final results'
     write(*,*) '# writing final results'
     write(67,*) ' found ',nlmin_l,nlmin,' minima'
     write(*,*) '# found ',nlmin_l,nlmin,' minima'
  endif


  call cpu_time(tcpu2)
  if (iproc == 0) then
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
  call deallocate_atoms(atoms,subname)
  call free_restart_objects(rst,subname)
  call free_input_variables(inputs_opt)
  call free_input_variables(inputs_md)


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

  if (iproc == 0) write(67,'(a,1x,3(1x,1pe10.3))') 'Out:ediff,ekinetic,dt',ediff,ekinetic,dt
  if (iproc == 0) write(*,'(a,1x,3(1x,1pe10.3))') '# Out:ediff,ekinetic,dt',ediff,ekinetic,dt
  close(2) 

  !Finalize memory counting
  call memocc(0,0,'count','stop')
  if (nproc > 1) call MPI_FINALIZE(ierr)


contains


  subroutine mdescape(nsoften,mdmin,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,rxyz, &
       nproc,iproc,atoms,rst,inputs_md)!  &
    ! Does a MD run with the atomic positiosn rxyz
    use module_base
    use module_types
    use module_interfaces
    use ab6_symmetry
    implicit real*8 (a-h,o-z)
    type(atoms_data) :: atoms
    type(restart_objects) :: rst
    dimension ff(3,atoms%nat),gg(3,atoms%nat),vxyz(3,atoms%nat),rxyz(3,atoms%nat),rxyz_old(3,atoms%nat)
    type(input_variables) :: inputs_md
    character(len=4) :: fn,name
    !type(wavefunctions_descriptors), intent(inout) :: wfd
    !real(kind=8), pointer :: psi(:), eval(:)

    if(iproc==0) write(*,*) '# MINHOP start soften ',nsoften

    !C initialize positions,velocities, forces
    call randdist(atoms%nat,rxyz,vxyz)
    inputs_md%inputPsiId=1
    !if(iproc==0)write(*,*)' #  no softening'
    call soften(nsoften,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,rxyz, &
         nproc,iproc,atoms,rst,inputs_md)
    call velopt(atoms,rxyz,ekinetic,vxyz)
    call razero(3*atoms%nat,gg)

    if(iproc==0) write(*,*) '# MINHOP start MD'
    !C inner (escape) loop
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
    md_loop: do istep=1,200

       !C      Evolution of the system according to 'VELOCITY VERLET' algorithm
!!       rkin=0.d0
!!       do iat=1,atoms%nat
!!          if (.not. atoms%lfrztyp(iat)) then
!!             if (atoms%geocode == 'P') then
!!                rxyz(1,iat)=modulo(rxyz(1,iat) + dt*vxyz(1,iat) + (.5d0*dt*dt)*gg(1,iat),&
!!                     atoms%alat1)
!!                rxyz(2,iat)=modulo(rxyz(2,iat) + dt*vxyz(2,iat) + (.5d0*dt*dt)*gg(2,iat),&
!!                     atoms%alat2)
!!                rxyz(3,iat)=modulo(rxyz(3,iat) + dt*vxyz(3,iat) + (.5d0*dt*dt)*gg(3,iat),&
!!                     atoms%alat3)
!!             else if (atoms%geocode == 'S') then
!!                rxyz(1,iat)=modulo(rxyz(1,iat) + dt*vxyz(1,iat) + (.5d0*dt*dt)*gg(1,iat),&
!!                     atoms%alat1)
!!                rxyz(2,iat)=       rxyz(2,iat) + dt*vxyz(2,iat) + (.5d0*dt*dt)*gg(2,iat)
!!                rxyz(3,iat)=modulo(rxyz(3,iat) + dt*vxyz(3,iat) + (.5d0*dt*dt)*gg(3,iat),&
!!                     atoms%alat3)
!!             else if (atoms%geocode == 'F') then
!!                rxyz(1,iat)=rxyz(1,iat) + dt*vxyz(1,iat) + (.5d0*dt*dt)*gg(1,iat)
!!                rxyz(2,iat)=rxyz(2,iat) + dt*vxyz(2,iat) + (.5d0*dt*dt)*gg(2,iat)
!!                rxyz(3,iat)=rxyz(3,iat) + dt*vxyz(3,iat) + (.5d0*dt*dt)*gg(3,iat)
!!             end if
!!             rkin=rkin+vxyz(1,iat)**2+vxyz(2,iat)**2+vxyz(3,iat)**2
!!          end if
!!       enddo
      call atomic_axpy(atoms,rxyz,dt,vxyz,rxyz)
      call atomic_axpy(atoms,rxyz,.5d0*dt*dt,gg,rxyz)
      call atomic_dot(atoms,vxyz,vxyz,rkin)
       rkin=rkin*.5d0

       enmin2=enmin1
       enmin1=en0000
       !    if (iproc == 0) write(*,*) 'CLUSTER FOR  MD'
       inputs_md%inputPsiId=1
       if (istep > 2) inputs_md%itermax=50
       call call_bigdft(nproc,iproc,atoms,rxyz,inputs_md,e_rxyz,ff,rst,infocode)

       if (iproc == 0) then
          write(fn,'(i4.4)') istep
          call write_atomic_file('posmd_'//fn,e_rxyz,rxyz,atoms,'')
       end if

       en0000=e_rxyz-e_pos
       if (istep >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (istep >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
       econs_max=max(econs_max,rkin+e_rxyz)
       econs_min=min(econs_min,rkin+e_rxyz)
       devcon=econs_max-econs_min
       if (iproc == 0) write(17,'(a,i5,1x,1pe17.10,2(1x,i2))') 'MD ',&
            istep,e_rxyz,nummax,nummin
       if (iproc == 0) write(*,'(a,i5,1x,1pe17.10,2(1x,i2))') ' #MD ',&
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
!          if (.not. atoms%lfrztyp(iat)) then
             vxyz(1,iat)=vxyz(1,iat) + (.5d0*dt) * (at1 + gg(1,iat))
             vxyz(2,iat)=vxyz(2,iat) + (.5d0*dt) * (at2 + gg(2,iat))
             vxyz(3,iat)=vxyz(3,iat) + (.5d0*dt) * (at3 + gg(3,iat))
!          end if
          !C Memorization of old forces
          gg(1,iat) = at1
          gg(2,iat) = at2
          gg(3,iat) = at3
       end do
    end do md_loop
    if (istep >=200) then
       if (iproc == 0) write(67,*) 'TOO MANY MD STEPS'
       dt=2.d0*dt
    end if
    !save the value of count_md for the moment
    count_md=count_md+real(istep,gp)

    !C MD stopped, now do relaxation

    !  if (iproc == 0) write(67,*) 'EXIT MD',istep
    
    ! adjust time step to meet precision criterion
    devcon=devcon/(3*atoms%nat-3)
    if (iproc == 0) &
         write(66,'(a,2(1x,1pe11.4),1x,i5)')&
         'MD devcon ',devcon,devcon/ekinetic,istep
    if (devcon/ekinetic.lt.10.d-2) then
       if (iproc == 0) write(66,*) 'MD:old,new dt',dt,dt*1.05d0
       dt=dt*1.05d0
    else
       if (iproc == 0) write(66,*) 'MD:old,new dt',dt,dt/1.05d0
       dt=dt*(1.d0/1.05d0)
    endif
    
  END SUBROUTINE mdescape
  
  

  subroutine soften(nsoften,ekinetic,e_pos,fxyz,gg,vxyz,dt,count_md,rxyz, &
       nproc,iproc,atoms,rst,inputs_md)! &
    use module_base
    use module_types
    use module_interfaces
    use ab6_symmetry
    implicit real*8 (a-h,o-z)
    type(atoms_data) :: atoms
    dimension fxyz(3*atoms%nat),gg(3*atoms%nat),vxyz(3*atoms%nat),rxyz(3*atoms%nat),rxyz_old(3*atoms%nat)
    type(input_variables) :: inputs_md
    type(restart_objects) :: rst
    !Local variables
    dimension wpos(3*atoms%nat)

!    eps_vxyz=1.d-1*atoms%nat
    alpha=inputs_md%betax

    !allocate(wpos(3,nat),fxyz(3,nat))

    inputs_md%inputPsiId=1
    if(iproc==0)write(*,*)'# soften initial step '
    call call_bigdft(nproc,iproc,atoms,rxyz,inputs_md,etot0,fxyz,rst,infocode)

    ! scale velocity to generate dimer 

    call atomic_dot(atoms,vxyz,vxyz,svxyz)
!    svxyz=0.d0
!    do i=1,3*atoms%nat
!       iat=(i-1)/3+1
!       if (.not. atoms%lfrztyp(iat)) then
!          svxyz=svxyz+vxyz(i)**2
!       end if
!    enddo
    eps_vxyz=sqrt(svxyz)
    if(iproc == 0) write(*,*)'#  eps_vxyz=',eps_vxyz

    do it=1,nsoften
       
!       do iat=1,atoms%nat
!          if (atoms%lfrztyp(iat)) then
!             wpos(3*(iat-1)+1)=rxyz(3*(iat-1)+1)
!             wpos(3*(iat-1)+2)=rxyz(3*(iat-1)+2)
!             wpos(3*(iat-1)+3)=rxyz(3*(iat-1)+3)
!          else
!             if (atoms%geocode == 'P') then
!                wpos(3*(iat-1)+1)=modulo(rxyz(3*(iat-1)+1)+vxyz(3*(iat-1)+1),atoms%alat1)
!                wpos(3*(iat-1)+2)=modulo(rxyz(3*(iat-1)+2)+vxyz(3*(iat-1)+2),atoms%alat2)
!                wpos(3*(iat-1)+3)=modulo(rxyz(3*(iat-1)+3)+vxyz(3*(iat-1)+3),atoms%alat3)
!             else if (atoms%geocode == 'S') then
!                wpos(3*(iat-1)+1)=modulo(rxyz(3*(iat-1)+1)+vxyz(3*(iat-1)+1),atoms%alat1)
!                wpos(3*(iat-1)+2)=       rxyz(3*(iat-1)+2)+vxyz(3*(iat-1)+2)
!                wpos(3*(iat-1)+3)=modulo(rxyz(3*(iat-1)+3)+vxyz(3*(iat-1)+3),atoms%alat3)
!             else if (atoms%geocode == 'F') then
!                wpos(3*(iat-1)+1)=rxyz(3*(iat-1)+1)+vxyz(3*(iat-1)+1)
!                wpos(3*(iat-1)+2)=rxyz(3*(iat-1)+2)+vxyz(3*(iat-1)+2)
!                wpos(3*(iat-1)+3)=rxyz(3*(iat-1)+3)+vxyz(3*(iat-1)+3)
!             end if
!          end if
!       end do
       call atomic_axpy(atoms,rxyz,1.d0,vxyz,wpos)

       call call_bigdft(nproc,iproc,atoms,wpos,inputs_md,etot,fxyz,rst,infocode)
       fd2=2.d0*(etot-etot0)/eps_vxyz**2

!       sdf=0.d0
!       svxyz=0.d0
!       do i=1,3*atoms%nat
!          iat=(i-1)/3+1
!          if (.not. atoms%lfrztyp(iat)) then
!             sdf=sdf+vxyz(i)*fxyz(i)
!             svxyz=svxyz+vxyz(i)*vxyz(i)
!          end if
!       end do
       call atomic_dot(atoms,vxyz,vxyz,svxyz)
       call atomic_dot(atoms,vxyz,fxyz,sdf)

       curv=-sdf/svxyz
       if (it == 1) curv0=curv

!       res=0.d0
!       do i=1,3*atoms%nat
!          iat=(i-1)/3+1
!          if (.not. atoms%lfrztyp(iat)) then
!             fxyz(i)=fxyz(i)+curv*vxyz(i)
!             res=res+fxyz(i)**2
!          end if
!       end do
       call atomic_axpy_forces(atoms,fxyz,curv,vxyz,fxyz)
       call atomic_dot(atoms,fxyz,fxyz,res)
       res=sqrt(res)

       write(fn4,'(i4.4)') it
       write(comment,'(a,1pe10.3)')'res= ',res
       if (iproc == 0) call write_atomic_file('possoft_'//fn4,etot,wpos,atoms,trim(comment))

       if(iproc==0)write(*,'(a,i3,5(f12.5),f10.3)')'# soften it, curv, fd2,dE,res,eps_vxyz:',&
            it, curv, fd2,etot-etot0,res,eps_vxyz
       if (curv.lt.0.d0 .or. fd2.lt.0.d0) then
          if(iproc==0) write(*,*) '# NEGATIVE CURVATURE'
          exit
       end if
       if (etot-etot0.lt.1.d-2) eps_vxyz=eps_vxyz*1.2d0

!       do iat=1,atoms%nat
!          if (.not. atoms%lfrztyp(iat)) then
!             if (atoms%geocode == 'P') then
!                wpos(3*(iat-1)+1)=modulo(wpos(3*(iat-1)+1)+alpha*fxyz(3*(iat-1)+1),atoms%alat1)
!                wpos(3*(iat-1)+2)=modulo(wpos(3*(iat-1)+2)+alpha*fxyz(3*(iat-1)+2),atoms%alat2)
!                wpos(3*(iat-1)+3)=modulo(wpos(3*(iat-1)+3)+alpha*fxyz(3*(iat-1)+3),atoms%alat3)
!             else if (atoms%geocode == 'S') then
!                wpos(3*(iat-1)+1)=modulo(wpos(3*(iat-1)+1)+alpha*fxyz(3*(iat-1)+1),atoms%alat1)
!                wpos(3*(iat-1)+2)=       wpos(3*(iat-1)+2)+alpha*fxyz(3*(iat-1)+2)
!                wpos(3*(iat-1)+3)=modulo(wpos(3*(iat-1)+3)+alpha*fxyz(3*(iat-1)+3),atoms%alat3)
!             else if (atoms%geocode == 'F') then
!                wpos(3*(iat-1)+1)=wpos(3*(iat-1)+1)+alpha*fxyz(3*(iat-1)+1)
!                wpos(3*(iat-1)+2)=wpos(3*(iat-1)+2)+alpha*fxyz(3*(iat-1)+2)
!                wpos(3*(iat-1)+3)=wpos(3*(iat-1)+3)+alpha*fxyz(3*(iat-1)+3)
!             end if
!
!          end if
!       end do
       call atomic_axpy_forces(atoms,wpos,alpha,fxyz,wpos)

       do i=1,3*atoms%nat
          vxyz(i)=wpos(i)-rxyz(i)
       end do
       write(comment,'(a,1pe10.3)')'curv= ',curv
       if (iproc == 0) call write_atomic_file('posvxyz',0.d0,vxyz,atoms,trim(comment))
       call elim_moment(atoms%nat,vxyz)
       call elim_torque(atoms%nat,rxyz,vxyz)

!       svxyz=0.d0
!       do i=1,3*atoms%nat
!          iat=(i-1)/3+1
!          if (.not. atoms%lfrztyp(iat)) then
!             svxyz=svxyz+vxyz(i)*vxyz(i)
!          end if
!       end do
       call atomic_dot(atoms,vxyz,vxyz,svxyz)
       if (res <= curv*eps_vxyz*5.d-1) exit
       svxyz=eps_vxyz/dsqrt(svxyz)

       do i=1,3*atoms%nat
          vxyz(i)=vxyz(i)*svxyz
       end do

    end do ! iter
 

    !        deallocate(wpos,fxyz)
  END SUBROUTINE soften


end program MINHOP
!!***



subroutine insert(nlminx,nlmin,k_e_wpos,re_wpos,earr)
  ! inserts the energy re_wpos at position k_e_wpos and shifts up all other energies
  implicit real*8 (a-h,o-z)
  dimension earr(0:nlminx,2)
  do k=nlmin-1,k_e_wpos+1,-1
     earr(k+1,1)=earr(k,1)
     earr(k+1,2)=earr(k,2)
  enddo
  earr(k_e_wpos+1,1)=re_wpos
  earr(k_e_wpos+1,2)=1.d0
  return
END SUBROUTINE insert


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
END SUBROUTINE save_low_conf


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
  if (n == 1) then
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
  if(jhi-jlo == 1)then
     if(x == xx(n))jlo=n
     if(x == xx(1))jlo=1
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
  use ab6_symmetry
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
!     if (.not. at%lfrztyp(iat)) then
        rkinsum= rkinsum+vxyz(1,iat)**2+vxyz(2,iat)**2+vxyz(3,iat)**2
!     end if
  end do
  rkin=.5d0*rkinsum/(3*at%nat-3)
  !       write(*,*) 'rkin,ekinetic',rkin,ekinetic

  !C      Rescaling of velocities to get reference kinetic energy
  sclvel= dsqrt(ekinetic/rkin)
  do iat=1,at%nat
!     if (.not. at%lfrztyp(iat)) then
        vxyz(1,iat)=vxyz(1,iat)*sclvel
        vxyz(2,iat)=vxyz(2,iat)*sclvel
        vxyz(3,iat)=vxyz(3,iat)*sclvel
!     end if
  end do

END SUBROUTINE velopt


subroutine randdist(nat,rxyz,vxyz)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(3*nat), intent(in) :: rxyz
  real(gp), dimension(3*nat), intent(out) :: vxyz
  !local variables
  integer :: i,idum=0
  real(kind=4) :: tt,builtin_rand
  ! create a random displacement vector without translational and angular moment
  do i=1,3*nat
     !call random_number(tt)
     !add built-in random number generator
     tt=builtin_rand(idum)
     vxyz(i)=real(tt-.5,gp)*3.e-1_gp
  end do

  call elim_moment(nat,vxyz)
  call elim_torque(nat,rxyz,vxyz)
END SUBROUTINE randdist


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
END SUBROUTINE gausdist


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
END SUBROUTINE expdist


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
END SUBROUTINE localdist


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

END SUBROUTINE torque


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
     if (ii == 1) then 
        cx=t(1)/(sz+sy)
        do iat=1,nat
           vxyz(2,iat)=vxyz(2,iat)+cx*(rxyz(3,iat)-cmz)
           vxyz(3,iat)=vxyz(3,iat)-cx*(rxyz(2,iat)-cmy)
        enddo
     else if(ii == 2) then 
        cy=t(2)/(sz+sx)
        do iat=1,nat
           vxyz(1,iat)=vxyz(1,iat)-cy*(rxyz(3,iat)-cmz)
           vxyz(3,iat)=vxyz(3,iat)+cy*(rxyz(1,iat)-cmx)
        enddo
     else if(ii == 3) then 
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

END SUBROUTINE elim_torque


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

END SUBROUTINE moment


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

END SUBROUTINE elim_moment


subroutine fix_fragmentation(iproc,at,rxyz,nputback)
  use module_base
  use module_types
  use ab6_symmetry
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: at
  integer, intent(inout) :: nputback
  real(gp), dimension(3,at%nat) :: rxyz
  !local variables
  real(gp), parameter :: bondlength=8.0_gp
  integer :: iat,nloop,ncluster,ii,jat,jj,kat,nadd,ierr
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

        if (iproc == 0) then
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
!           if (.not. belong(iat) .and. .not. at%lfrztyp(iat)) then
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
!           endif
        enddo

        if (iproc == 0) then
           write(444,*) at%nat, 'atomic ' 
           write(444,*) ' fixed configuration ', nputback
           do iat=1,at%nat
              write(444,*) ' LJ  ',rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
           enddo
        endif
        nloop=nloop+1
     if (nloop.gt.4) call MPI_ABORT(MPI_COMM_WORLD,ierr)
     endif
  end do fragment_loop

END SUBROUTINE fix_fragmentation


subroutine winter(at,re_pos,pos,npminx,nlminx,nlmin,nlmin_l,accur, & 
     earr,elocmin,poslocmin,eref,ediff,ekinetic,dt,nsoften)
  use module_base
  use module_types
  use ab6_symmetry
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: npminx,nlminx,nlmin,nlmin_l,nsoften
  real(gp), intent(in) :: re_pos,eref,ediff,ekinetic,dt,accur
  type(atoms_data), intent(in) :: at
  real(gp), dimension(npminx), intent(in) :: elocmin
  real(gp), dimension(3,at%nat), intent(in) :: pos
  real(gp), dimension(0:nlminx,2), intent(in) :: earr
  real(gp), dimension(3,at%nat,npminx), intent(in) :: poslocmin
  !local variables
  character(len=5) :: fn
  character(len=20) :: filename
  integer :: mm,k

     call write_atomic_file('poscur',re_pos,pos,at,'')
     write(*,*) ' wrote poscur.xyz for  RESTART'
     
  call wtpos(at,npminx,nlminx,nlmin,nlmin_l,poslocmin,earr,elocmin)

     write(*,*) ' wrote poslow files'
     
     open(unit=12,file='earr.dat',status='unknown')
     mm=min(nlmin,nlminx)
     write(12,'(2(i10),a)') mm,mm+10,&
          ' # of minima already found, # of minima to be found in consecutive run'
     write(12,'(e24.17,1x,a)') eref,'   eref'
     write(12,'(e24.17,1x,a)') accur,'   accur'
     do k=1,mm
        write(12,'(e24.17,1x,1pe17.10)') earr(k,1),earr(k,2)
     enddo
     write(*,*) ' wrote earr.dat for  RESTART'
     close(12)
     
     call  wtioput(ediff,ekinetic,dt,nsoften)
     write(*,*) ' wrote ioput for  RESTART'

END SUBROUTINE winter


subroutine wtioput(ediff,ekinetic,dt,nsoften)
  implicit real*8 (a-h,o-z)
  open(unit=11,file='ioput',status='unknown')
  write(11,'(3(1x,1pe24.17)1x,i4,a)') ediff,ekinetic,dt,nsoften,' ediff, ekinetic dt and nsoften'
  close(11)
END SUBROUTINE wtioput


subroutine wtpos(at,npminx,nlminx,nlmin,nlmin_l,pos,earr,elocmin)
  use module_base
  use module_types
  use ab6_symmetry
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: npminx,nlminx,nlmin,nlmin_l
  type(atoms_data), intent(in) :: at
  real(gp), dimension(npminx), intent(in) :: elocmin
  real(gp), dimension(0:nlminx,2), intent(in) :: earr
  real(gp), dimension(3,at%nat,npminx), intent(in) :: pos
  !local variables
  character(len=5) :: fn
  character(len=17) :: filename
  character(len=20) :: atomnames
  integer :: k,kk,i

       write(*,*) 'nlmin,nlminx,nlmin_l,npminx',nlmin,nlminx,nlmin_l,npminx
       do i=1,min(40,nlmin,nlminx)
         write(*,'(i4,e24.17)') i,earr(i,1)
       enddo

  do k=1,min(nlmin_l,npminx)
             write(*,'(a,i4,e24.17)') 'k,elocmin(k)',k,elocmin(k)


     !C Classify the configuration in the global ranking
     kk=0
     find_kk : do
        kk=kk+1
        if (kk > min(nlmin,nlminx)) then 
           write(*,*) 'ranking error for',k
           stop 
        endif
!        if (earr(kk,1) == elocmin(k)) exit find_kk
        if (abs(earr(kk,1) - elocmin(k)) .lt. 1.d-12 ) then 
             write(*,*) 'match ',abs(earr(kk,1) - elocmin(k))
             exit find_kk
        endif
     end do find_kk

     if (kk <= npminx) then

        !        write(*,'(a,i2,i4,i4,1x,1pe21.14)') 'k,kk,elocmin(k)',k,kk,elocmin(k)

        !C generate filename and open files
           write(fn,'(i5.5)') kk
           call  write_atomic_file('poslow'//fn,elocmin(k),pos(1,1,k),at,'')
     endif

  end do

END SUBROUTINE wtpos


real*8 function round(enerd,accur)
  implicit none
  real*8 enerd,accur
  integer*8 ii
  ii=enerd/accur
  round=ii*accur
  !           write(*,'(a,1pe24.17,1x,i17,1x,1pe24.17)') 'enerd,ii,round',enerd,ii,round
  return
end function round


!subroutine rdposout(igeostep,rxyz,nat)
!  implicit none
!  integer, intent(in) :: igeostep,nat
!  real(kind=8), dimension(3,nat), intent(out) :: rxyz
!  !local variables
!  character(len=3) :: fn
!  character(len=20) :: filename
!  integer :: iat
!  write(fn,'(i3.3)') igeostep
!  !filename = 'posout_'//fn//'.ascii'
!  filename = 'posout_'//fn//'.xyz'
!  ! write(*,*)'# reading unrelaxed structure from ',filename
!  open(unit=9,file=filename,status='old')
!  read(9,*)fn!no need for header
!  read(9,*)fn!same
!  do iat=1,nat
!     read(9,*)fn,rxyz(:,iat)!we know the atom types already
!  enddo
!  close(unit=9)
!END SUBROUTINE rdposout


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
END SUBROUTINE adjustrxyz
