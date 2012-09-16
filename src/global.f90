!> @file
!!  Minima hopping program
!! @author
!!    New modified version 17th Nov 2009 Sandip De  
!!    Copyright (C) 2008-2011 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS


!> MINHOP
!!  Main program for the minima hopping
program MINHOP

  use module_base
  use module_types
  use module_interfaces
  use m_ab6_symmetry
  use yaml_output
  implicit real(kind=8) (a-h,o-z)
  real(kind=4) :: tts
  logical :: newmin,CPUcheck,occured,exist_poslocm
  character(len=20) :: unitsp,units,atmn
  character(len=80) :: line
  type(atoms_data) :: atoms,md_atoms
  type(input_variables) :: inputs_opt, inputs_md
  type(restart_objects) :: rst
  integer, parameter :: npminx=200
  !C parameters for minima hopping
  integer, parameter :: mdmin=2
  real(kind=8), parameter :: beta1=1.10d0,beta2=1.10d0,beta3=1.d0/1.10d0
  real(kind=8), parameter :: alpha1=1.d0/1.10d0,alpha2=1.10d0
  real(gp):: fnoise
  real(kind=8) :: elocmin(npminx)
  real(kind=8), allocatable, dimension(:,:) ::ff,wpos,vxyz,gg,earr,poshop
  real(kind=8), allocatable, dimension(:) :: rcov,evals
  real(kind=8),allocatable, dimension(:,:,:):: poslocmin
  real(kind=8), dimension(:,:), pointer :: pos,mdpos
  integer :: iproc,nproc,iat,ityp,j,i_stat,i_all,ierr,infocode,norbs_eval
  character(len=*), parameter :: subname='global'
  character(len=41) :: filename
  character(len=4) :: fn4
  character(len=5) :: fn5
  character(len=16) :: fn16
  character(len=50) :: comment
  real(gp), dimension(6) :: strten
!  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem

  ! Start MPI version
  call bigdft_mpi_init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  !call system('echo $HOSTNAME')

  !if (iproc ==0) call yaml_set_stream(unit=70,filename='global-logfile.yaml')
!  if (iproc ==0) call yaml_set_stream(record_length=92)!unit=70,filename='log.yaml')
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
     write(*,'(23x,a)')'----> you can grep this file for #MH to compare with global.out'
     write(*,'(23x,a)')' #MH NOTE: this version reads nspin, mpol from input.dat'
  end if

!  open(unit=67,file='global.out')


  if (iproc == 0) write(*,'(a,3(1x,1pe11.4))') '#MH beta1,beta2,beta3',beta1,beta2,beta3
  if (iproc == 0) write(*,'(a,2(1x,1pe11.4))') '#MH alpha1,alpha2',alpha1,alpha2
  !if (iproc == 0) write(*,'(a,2(1x,1pe10.3))') '#MH predicted fraction accepted, rejected', & 
  !     ratio/(1.d0+ratio), 1.d0/(1.d0+ratio)
  if (iproc == 0) write(*,*) '#MH mdmin',mdmin
  accepted=0.0d0


  call cpu_time(tcpu1)
    do i=1,npminx
      elocmin(i)=1.d100
    enddo                                                               
  ! read  earr.dat
  open(unit=12,file='earr.dat',status='unknown')
  read(12,*) nlmin,nlminx
  read(12,*) eref
  if (iproc == 0) write(*,*) '#MH eref=',eref
  read(12,*) accur
  if (iproc == 0) write(*,*) '#MH accuracy for rounding=',accur
  if (nlmin.gt.nlminx) stop 'nlmin>nlminx'
  allocate(earr(0:nlminx,2+ndebug),stat=i_stat)
  call memocc(i_stat,earr,'earr',subname)
  earr(0,1)=-1.d100
  if (nlmin == 0) then 
     if (iproc == 0) write(*,*) '#MH New run with nlminx=',nlminx
  else
     if (iproc == 0) write(*,*) '#MH Restart run with nlmin, nlminx=',nlmin,nlminx
     do k=1,nlmin
        read(12,*) earr(k,1),earr(k,2)
        if (earr(k,1).lt.earr(k-1,1)) stop 'wrong ordering in earr.dat'
     enddo
     if (iproc == 0) write(*,*) '#MH read earr.dat'
  endif
  close(12)

  call standard_inputfile_names(inputs_opt,'input',nproc)
  call standard_inputfile_names(inputs_md,'mdinput',nproc)

  call read_atomic_file('poscur',iproc,atoms,pos)

  !Read input parameters for geometry optimization 
  call read_input_parameters(iproc,inputs_opt,atoms,pos)

!!$  call default_input_variables(inputs_opt)
!!$  call dft_input_variables_new(iproc,'input.dft',inputs_opt)
!!$  call geopt_input_variables('input.geopt',inputs_opt)
!!$  call kpt_input_variables(iproc,'input.kpt',inputs_opt,atoms)

  !read input parameters for molecular dynamics
  call read_atomic_file('poscur',iproc,md_atoms,mdpos)
  call read_input_parameters(iproc,inputs_md,md_atoms,pos)
!!$  call default_input_variables(inputs_md)
!!$  call dft_input_variables_new(iproc,'mdinput.dft',inputs_md)
!!$  call geopt_input_variables('mdinput.geopt',inputs_md)
!!$  call kpt_input_variables(iproc,'input.kpt',inputs_md,atoms)


  !use only the atoms structure for the run
  call init_atomic_values((iproc == 0),md_atoms,inputs_md%ixc)
  call deallocate_atoms(md_atoms,subname) 
  i_all=-product(shape(mdpos))*kind(mdpos)
  deallocate(mdpos,stat=i_stat)
  call memocc(i_stat,i_all,'mdpos',subname)

  ! Read associated pseudo files. Based on the inputs_opt set
  call init_atomic_values((iproc == 0), atoms, inputs_opt%ixc)
  call read_atomic_variables(atoms, trim(inputs_opt%file_igpop),inputs_opt%nspin)

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
  allocate(poshop(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,poshop,'poshop',subname)
  allocate(rcov(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rcov,'rcov',subname)

  call give_rcov(iproc,atoms,atoms%nat,rcov)

! read random offset
  open(unit=11,file='rand.inp')
  read(11,*) nrandoff
  !        write(*,*) 'nrandoff ',nrandoff
  close(11)
  do i=1,nrandoff
     call random_number(ts)
  enddo

! open output files
       if (iproc==0) then 
          open(unit=2,file='global.mon',status='unknown',position='append')
          open(unit=16,file='geopt.mon',status='unknown')
       endif
  
  ! read input parameters
  write(filename,'(a6,i3.3)') 'ioput'   !,iproc
  open(unit=11,file='ioput',status='old')
  read(11,*) ediff,ekinetic,dt,nsoften
  close(11)
  !write(*,'(a,1x,i3,3(1x,e10.3),1x,i4)') 'In :iproc,ediff,ekinetic,dt,nsoften',iproc,ediff,ekinetic,dt,nsoften
  if (iproc == 0) write(*,'(a,1x,3(1x,e10.3),1x,i4)') 'In :ediff,ekinetic,dt,nsoften',ediff,ekinetic,dt,nsoften


  ! If restart run read previously found energies
  if (nlmin > 0) then
     npmin=0
     kk=1
     do
        write(fn5,'(i5.5)') kk
        filename = 'poslow'//fn5//'.xyz'
        if (npmin.ge.npminx .or. kk.gt.min(nlmin,npminx)) then
           exit
        end if
        open(unit=9,file=filename,status='old',iostat=ierror)
        if (ierror == 0) then
            npmin=npmin+1  
        else
!           write(*,*) iproc,' COULD not read file ',filename
           exit
        end if
        read(9,*) natp,unitsp,elocmin(npmin)
        if (atoms%nat.ne.natp) stop   'nat <> natp'
        if (trim(unitsp).ne.trim(atoms%units) .and. iproc.eq.0) write(*,*)  & 
                 '#MH different units in poslow and poscur file: ',trim(unitsp),' ',trim(atoms%units)
        read(9,*)
        do iat=1,atoms%nat
          read(9,*) atmn,t1,t2,t3
          if (atoms%units=='angstroem' .or. atoms%units=='angstroemd0') then ! if Angstroem convert to Bohr
              poslocmin(1,iat,npmin)=t1/bohr2ang 
              poslocmin(2,iat,npmin)=t2/bohr2ang 
              poslocmin(3,iat,npmin)=t3/bohr2ang
          else
              poslocmin(1,iat,npmin)=t1
              poslocmin(2,iat,npmin)=t2
              poslocmin(3,iat,npmin)=t3
          endif
        enddo
        close(9)
        if (iproc == 0) write(*,*) '#MH read file',filename
        kk=kk+1
     end do
     if (iproc == 0) then 
        write(*,*) 'read ',npmin,'#MH poslow files with energies'
        do ip=1,npmin
          write(*,*) elocmin(ip)
        enddo                                                              
     endif

     if (iproc == 0) write(*,*) '#MH read ',npmin,'poslow files'
  endif

  av_ekinetic=0.d0
  av_ediff=0.d0
  escape=0.d0
  escape_sam=0.d0
  escape_old=0.d0
  escape_new=0.d0
  rejected=0
  accepeted=0
 ! hopp=0.d0
 ! hopp_acc=0.d0
 ! hopp_rej=0.d0
  egap=1.d100
  esep=0.d0
  e_hop=1.d100
!C first local minimum
  count_sd=0.d0
  count_cg=0.d0
  count_soft=0.d0
  count_md=0.d0
  nputback=0

  inputs_opt%inputPsiId=0
  call init_restart_objects(iproc,inputs_opt%matacc,atoms,rst,subname)
  call call_bigdft(nproc,iproc,atoms,pos,inputs_md,e_pos,ff,strten,fnoise,rst,infocode)

  !example for retrieving the eigenvalues from this run
  norbs_eval=bigdft_get_number_of_orbitals(rst,i_stat)
  if (i_stat /= BIGDFT_SUCCESS) then
     write(*,*)'error (norbs), i_stat',i_stat
     stop
  end if
  allocate(evals(norbs_eval+ndebug),stat=i_stat)
  call memocc(i_stat,evals,'evals',subname)
  call bigdft_get_eigenvalues(rst,evals,i_stat)
  if (i_stat /= BIGDFT_SUCCESS) then
     write(*,*)'error(evals), i_stat',i_stat
     stop
  end if

  i_all=-product(shape(evals))*kind(evals)
  deallocate(evals,stat=i_stat)
  call memocc(i_stat,i_all,'evals',subname)


  if (iproc==0)write(17,*) 'ENERGY ',e_pos
  energyold=1.d100
  ncount_bigdft=0

  if (iproc == 0) write(*,*)'#MH calling conjgrad for the first time here. energy ',e_pos

!  if (atoms%geocode == 'P') & 
!       call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,pos)

  nconjgr=0
      do 
        write(fn16,'(a8,i4.4,a4)') "poslocm_",nconjgr,".xyz"
        inquire(file=fn16,exist=exist_poslocm)
        if (exist_poslocm) then
            nconjgr=nconjgr+1
        else
            exit
        endif
       enddo
       if (iproc == 0) write(*,*) '#MH number of poslocm files that exist already ',nconjgr


  call geopt(nproc,iproc,pos,atoms,ff,strten,e_pos,rst,inputs_md,ncount_bigdft)
  if (iproc == 0) then
     write(*,*) '#MH ', ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
  end if

  if (iproc == 0) then 
     tt=dnrm2(3*atoms%nat,ff,1)
     write(fn4,'(i4.4)') nconjgr
     write(comment,'(a,1pe10.3)')'fnrm= ',tt
     call write_atomic_file('posimed_'//fn4,e_pos-eref,pos,atoms,trim(comment),forces=ff)
  endif

  call geopt(nproc,iproc,pos,atoms,ff,strten,e_pos,rst,inputs_opt,ncount_bigdft)
  if (iproc == 0) then
     write(*,*) '#MH ', ncount_bigdft,' Wvfnctn Opt. steps for accurate initial conf'
  end if


  if (iproc == 0) then 
     tt=dnrm2(3*atoms%nat,ff,1)
     write(fn4,'(i4.4)') nconjgr
     write(comment,'(a,1pe10.3)')'fnrm= ',tt
     call write_atomic_file('poslocm_'//fn4,e_pos-eref,pos,atoms,trim(comment),forces=ff)
  endif
  nconjgr=nconjgr+1

  re_pos=round(e_pos-eref,accur)
  if (iproc == 0) then
     write(*,'(a,1x,i3,3(1x,1pe17.10))') ' #MH INPUT(relaxed): iproc, e_pos,re_pos,eref ',iproc,e_pos,re_pos,eref
  end if


  if (nlmin.gt.0) then
     if (iproc == 0) write(*,'(a,2(1x,1pe24.17))') '#MH new/old energy for input file',re_pos
  endif
  k_e_wpos=1
  if (nlmin == 0) then
     nlmin=1
     npmin=1
     earr(1,1)=re_pos
     earr(1,2)=1.d0
     elocmin(1)=re_pos
     do iat=1,atoms%nat
        poslocmin(1,iat,1)=pos(1,iat) 
        poslocmin(2,iat,1)=pos(2,iat) 
        poslocmin(3,iat,1)=pos(3,iat) 
     enddo
  else  ! the poscur file might have been modified by hand
!         check whether new minimum  
        call hunt_g(earr(1,1),min(nlmin,nlminx),re_pos,k_e_pos)
            if (re_pos.eq.earr(k_e_pos,1)) then  
              if (iproc == 0) write(*,*) '#MH initial minimum is old '
            else
              if (iproc == 0) write(*,*) '#MH initial minimum is new '
              nlmin=nlmin+1
!            add minimum to history list
              call insert(nlminx,nlmin,k_e_pos,re_pos,earr(0,1))
!            save configuration if it is among the lowest ones in energy 
              npmin=npmin+1
              call save_low_conf(atoms%nat,npmin,npminx,re_pos,pos,elocmin,poslocmin)
              k_e_wpos=k_e_wpos+1
            endif
  endif
  re_sm=min(re_pos,earr(1,1))
  if (iproc == 0) then
     write(*,*) '#MH iproc,initial re_sm',iproc,re_sm
     write(2,'((1x,f10.0),1x,1pe21.14,2(1x,1pe10.3))')escape,e_pos-eref,ediff,ekinetic
  end if

  nlmin_old=nlmin
  CPUcheck=.false.

  !C outer (hopping) loop
!  hopping_loop: do
1000 continue
     if (nlmin >= nlminx) then 
!        write(67,*)iproc,'has  nlminx collected',nlmin
        write(*,*)'#MH:process', iproc,'has  nlminx collected',nlmin
        !             do i=1,nlmin ;  write(*,*) earr(i,1) ; enddo
        goto 3000
     endif
     !            Energy has reached taregt eref and global minimum is presumably found
     if (re_sm <= 1.d-3) then
        write(*,*)'#MH process', iproc,'success: relative energy < 0.001'
        goto 3000
     endif

5555 continue

!C check whether CPU time exceeded
     tleft=1.d100
     call cpu_time(tcpu2)
     if(iproc==0 .and. CPUcheck)then
        open(unit=55,file='CPUlimit_global',status='unknown')
        read(55,*,end=555) cpulimit 
        cpulimit=cpulimit*3600
        write(*,'(a,i5,i3,2(1x,e9.2))') '#MH iproc,nlmin,tcpu2-tcpu1,cpulimit',iproc,nlmin,tcpu2-tcpu1,cpulimit
        tleft=cpulimit-(tcpu2-tcpu1)
       end if
555    continue
       close(55)
       call MPI_BCAST(tleft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       if (tleft < 0.d0) then
       if(iproc==0) write(*,*) '#MH CPU time exceeded',iproc,tleft
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

!     call fix_fragmentation(iproc,atoms,wpos,nputback)
     if (atoms%geocode == 'F') call fixfrag_posvel(iproc,atoms%nat,rcov,wpos,vxyz,1,occured)
     if (atoms%geocode == 'S') call fixfrag_posvel_slab(iproc,atoms%nat,rcov,wpos,vxyz,1)
     

     av_ekinetic=av_ekinetic+ekinetic
     ncount_bigdft=0

!5556 continue ! entry point for restart of optimization at cluster step irestart+1
!     if (atoms%geocode == 'P') & 
!        call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,wpos)
     call geopt(nproc,iproc,wpos,atoms,ff,strten,e_wpos,rst,inputs_md,ncount_bigdft)

     if (iproc == 0) write(*,*)'#MH ', ncount_bigdft,' Wvfnctn Opt. steps for approximate geo. rel of MD conf.'
     !ncount_bigdft=0
     !ncount_cluster=0
!     if (atoms%geocode == 'P') & 
!          call  adjustrxyz(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,wpos)

     if (iproc == 0) then 
        tt=dnrm2(3*atoms%nat,ff,1)
        !call wtlmin(nconjgr,atoms%nat,e_wpos-eref,tt,wpos,atoms%iatype,atoms%atomnames,atoms%natpol)
        write(fn4,'(i4.4)') nconjgr
        write(comment,'(a,1pe10.3)')'fnrm= ',tt
        call write_atomic_file('posimed_'//fn4,e_wpos-eref,wpos,atoms,trim(comment),forces=ff)
     endif

      call geopt(nproc,iproc,wpos,atoms,ff,strten,e_wpos,rst,inputs_opt,ncount_bigdft)

     if (iproc == 0) write(*,*)'#MH ', ncount_bigdft,' Wvfnctn Opt. steps for accurate geo. rel of MD conf.'
     if (iproc == 0) then 
        tt=dnrm2(3*atoms%nat,ff,1)
        !call wtlmin(nconjgr,atoms%nat,e_wpos-eref,tt,wpos,atoms%iatype,atoms%atomnames,atoms%natpol)
        write(fn4,'(i4.4)') nconjgr
        write(comment,'(a,1pe10.3)')'fnrm= ',tt
        call write_atomic_file('poslocm_'//fn4,e_wpos-eref,wpos,atoms,trim(comment),forces=ff)
     endif
  nconjgr=nconjgr+1
  re_wpos=round(e_wpos-eref,accur)
  if (iproc == 0) write(*,'(a,i3,i3,4(1x,1pe14.7))')  & 
       '#MH npmin,nlmin,e_wpos,e_pos,re_wpos,re_pos', npmin,nlmin,e_wpos,e_pos,re_wpos,re_pos
  !C not escaped
  if (re_pos == re_wpos) then
     escape_sam=escape_sam+1.d0
     esep=esep+(e_pos-e_wpos)**2
     ekinetic=ekinetic*beta1
     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)
     if (iproc == 0) write(2,'((1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),a)')  &
          escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     if (iproc == 0) write(*,'(a)')' #MH no escape from current minimum.'
     goto 5555
  endif



  !C continue since escaped
  !C  check whether new minimum
  call hunt_g(earr(1,1),min(nlmin,nlminx),re_wpos,k_e_wpos)
  if (re_wpos == earr(k_e_wpos,1)) then
     if (iproc == 0) write(*,'(a,i3,i3,i4,1x,1pe14.7)')  & 
          ' #MH Revisited: npmin,nlmin,k_e_wpos,re_wpos=earr',npmin,nlmin,k_e_wpos,re_wpos
     newmin=.false.
     escape_old=escape_old+1.d0
     earr(k_e_wpos,2)=earr(k_e_wpos,2)+1.d0
     ekinetic=ekinetic*beta2
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

     nlmin=nlmin+1
     call insert(nlminx,nlmin,k_e_wpos,re_wpos,earr(0,1))                                   
     k_e_wpos=k_e_wpos+1
     npmin=npmin+1
     call save_low_conf(atoms%nat,npmin,npminx,re_wpos,wpos,elocmin,poslocmin)
  endif

!  hopp=hopp+1.d0
  if (e_wpos.lt.e_hop) then                                                                                 
    e_hop=e_wpos                                                                                           
    re_hop=re_wpos                                                                                         
    nvisit=int(earr(k_e_wpos,2))                                                                           
    do iat=1,atoms%nat                                                                                           
      poshop(1,iat)=wpos(1,iat) ; poshop(2,iat)=wpos(2,iat) ; poshop(3,iat)=wpos(3,iat)                    
    enddo                                                                                                  
  endif                             

  !C master: Monte Carlo step for local minima hopping
  av_ediff=av_ediff+ediff
  if (e_hop-e_pos.lt.ediff) then 
     !C          local minima accepted -------------------------------------------------------
   accepted=accepted+1.d0                                                                                  
   e_pos=e_hop                                                                                             
   re_pos=re_hop     

     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)
     e_pos=e_hop
     re_pos=re_hop
     do iat=1,atoms%nat
        pos(1,iat)=poshop(1,iat) 
        pos(2,iat)=poshop(2,iat) 
        pos(3,iat)=poshop(3,iat)
     enddo
!     if (newmin) then
!         npmin=npmin+1
!        call save_low_conf(atoms%nat,npmin,npminx,re_pos,pos,elocmin,poslocmin)
!     endif
      if (iproc == 0) then
      if (re_wpos.eq.re_hop) then  
         write(2,'((1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),l3,a,i5)')  &
          escape,e_hop-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A ', nvisit
      else 
         write(2,'(1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),a,i5)')  &                                           
              escape,e_wpos-eref,ediff,ekinetic, &                                                               
              escape_sam/escape,escape_old/escape,escape_new/escape,' I ',int(earr(k_e_wpos,2))                  
         write(2,'(1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),a,i5)')  &                                           
              escape,e_hop-eref,ediff,ekinetic, &                                                                
              escape_sam/escape,escape_old/escape,escape_new/escape,' A ',nvisit                                 
      endif                                                                                                    
      endif                                                                                                    
      e_hop=1.d100                                                                                            
      ediff=ediff*alpha1                               
! write intermediate results
      if (iproc == 0) write(*,*) 'WINTER'
      if (iproc == 0) call winter(atoms,re_pos,pos,npminx,nlminx,nlmin,npmin,accur, & 
           earr,elocmin,poslocmin,eref,ediff,ekinetic,dt,nsoften)
      goto 1000
  else
     !C          local minima rejected -------------------------------------------------------
     inputs_opt%inputPsiId=0  !ALEX says: Better do an input guess for the next escape
     if (iproc == 0) write(2,'((1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),3(1x,0pf5.2),l3,a,i5)')  &
          escape,e_wpos-eref,ediff,ekinetic, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' R ', int(earr(k_e_wpos,2))
     if (iproc == 0) write(*,'(a,1pe21.14)')' #MH rejected: ew-e>ediff ',e_wpos-e_pos

     rejected=rejected+1.d0
     ediff=ediff*alpha2
     if (iproc == 0) call wtioput(ediff,ekinetic,dt,nsoften)
     goto 1000
     !C                          ------------------------------------------------------------
  endif

!end do hopping_loop

3000 continue

  if (iproc == 0) then
     write(*,*) '#MH writing final results'
     write(*,*) '#MH found in total ',nlmin,' minima'
     write(*,*) '#MH Accepted ',accepted,' minima'
     call winter(atoms,re_pos,pos,npminx,nlminx,nlmin,npmin,accur, & 
           earr,elocmin,poslocmin,eref,ediff,ekinetic,dt,nsoften)
  endif


  call cpu_time(tcpu2)
  if (iproc == 0) then
     !C ratios from all the global counters
     write(*,'(i2,1x,a,3(1x,1pe10.3))') iproc,'#MH ratio stuck,same,old,new', &
          escape_sam/escape,escape_old/escape,escape_new/escape
     write(*,'(i2,1x,a,2(1x,1pe10.3))') iproc,'#MH ratio acc,rej',accepted/(accepted+rejected),rejected/(accepted+rejected)
     write(*,'(i2,1x,a,3(1x,f12.1))') iproc,'#MH count_md,count_sd,count_cg',count_md,count_sd,count_cg
     write(*,'(i2,1x,a,1x,1pe10.3)') iproc,'cpu(hrs) ', (tcpu2-tcpu1)/3600.d0
     write(*,'(i2,1x,a,2(1x,1pe10.3))') &
          iproc,'#MH average ediff, ekinetic',av_ediff/(accepted+rejected),av_ekinetic/escape
     write(*,'(a,1x,i8)') '#MH number of configurations for which atoms escaped ',nputback

     tt=0.d0
     ss=0.d0
     do i=1,nlmin
        tt=max(tt,earr(i,2))
        ss=ss+earr(i,2)
     enddo
     write(*,'(i2,a,f8.0)') iproc,' #MH  most frequent visits ',tt
     write(*,'(i2,a,1pe10.3)') iproc,'#MH   av. numb. visits per minimum',ss/nlmin
     write(*,'(a,e9.2)') '#MH minimum energy separation between presumably different configurations',egap
     if (escape_sam.gt.0) then
        esep=sqrt(esep/escape_sam)
        write(*,'(a,e9.2)') '#MH average energy separation between presumably identical configurations',esep
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

  i_all=-product(shape(poshop))*kind(poshop)
  deallocate(poshop,stat=i_stat)
  call memocc(i_stat,i_all,'poshop',subname)

  i_all=-product(shape(rcov))*kind(rcov)
  deallocate(rcov,stat=i_stat)
  call memocc(i_stat,i_all,'rcov',subname)
  if (iproc == 0) write(*,'(a,1x,3(1x,1pe10.3))') '#MH Out:ediff,ekinetic,dt',ediff,ekinetic,dt
  close(2) 

  !Finalize memory counting
  call memocc(0,0,'count','stop')
  call MPI_FINALIZE(ierr)


contains


  !> Does a MD run with the atomic positiosn rxyz
  subroutine mdescape(nsoften,mdmin,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,rxyz, &
       nproc,iproc,atoms,rst,inputs_md)!  &
    use module_base
    use module_types
    use module_interfaces
    use m_ab6_symmetry
    implicit real*8 (a-h,o-z)
    type(atoms_data) :: atoms
    type(restart_objects) :: rst
    dimension ff(3,atoms%nat),gg(3,atoms%nat),vxyz(3,atoms%nat),rxyz(3,atoms%nat),rxyz_old(3,atoms%nat),strten(6)
    type(input_variables) :: inputs_md
    character(len=4) :: fn,name
    logical :: move_this_coordinate
    !type(wavefunctions_descriptors), intent(inout) :: wfd
    !real(kind=8), pointer :: psi(:), eval(:)

    if(iproc==0) write(*,*) '#MH MINHOP start soften ',nsoften

    !C initialize positions,velocities, forces

  !! Either random velocity distribution 
  !        call randdist(nat,rxyz,vxyz)
  !! or Gauss velocity distribution
  !! or exponential  velocity distribution
  !        call expdist(nat,rxyz,vxyz)
  !! or localized velocities
  !        call localdist(nat,rxyz,vxyz)
    call randdist(atoms%nat,rxyz,vxyz)
    inputs_md%inputPsiId=1
    !if(iproc==0)write(*,*)' #MH  no softening'
  ! Soften previous velocity distribution
    call soften(nsoften,ekinetic,e_pos,ff,gg,vxyz,dt,count_md,rxyz, &
         nproc,iproc,atoms,rst,inputs_md)
  ! put velocities for frozen degrees of freedom to zero
       ndfree=0.d0
       ndfroz=0.d0
  do iat=1,atoms%nat
  do ixyz=1,3
  if ( move_this_coordinate(atoms%ifrztyp(iat),ixyz) ) then
       ndfree=ndfree+1
  else
       ndfroz=ndfroz+1
       vxyz(ixyz,iat)=0.d0
  endif
  enddo
  enddo
  ! normalize velocities to target ekinetic
    call velnorm(atoms,rxyz,(ekinetic*ndfree)/(ndfree+ndfroz),vxyz)
    call razero(3*atoms%nat,gg)

    if(iproc==0) call torque(atoms%nat,rxyz,vxyz)

    if(iproc==0) write(*,*) '#MH MINHOP start MD',ndfree,ndfroz
    !C inner (escape) loop
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
    istepnext=5
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
!      call atomic_axpy(atoms,rxyz,dt,vxyz,rxyz)
!      call atomic_axpy(atoms,rxyz,.5d0*dt*dt,gg,rxyz)
!      call atomic_dot(atoms,vxyz,vxyz,rkin)
call daxpy(3*atoms%nat,dt,vxyz(1,1),1,rxyz(1,1),1)                                                           
call daxpy(3*atoms%nat,0.5_gp*dt*dt,gg(1,1),1,rxyz(1,1),1)                                                   

rkin=dot(3*atoms%nat,vxyz(1,1),1,vxyz(1,1),1)     
       rkin=rkin*.5d0

       enmin2=enmin1
       enmin1=en0000
       !    if (iproc == 0) write(*,*) 'CLUSTER FOR  MD'
       inputs_md%inputPsiId=1
       call call_bigdft(nproc,iproc,atoms,rxyz,inputs_md,e_rxyz,ff,strten,fnoise,rst,infocode)

       if (iproc == 0) then
          write(fn,'(i4.4)') istep
          call write_atomic_file(trim(inputs_md%dir_output)//'posmd_'//fn,e_rxyz,rxyz,atoms,'',forces=ff)
       end if

       en0000=e_rxyz-e_pos
       if (istep >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (istep >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
       econs_max=max(econs_max,rkin+e_rxyz)
       econs_min=min(econs_min,rkin+e_rxyz)
       devcon=econs_max-econs_min
       if (iproc == 0) write(17,'(a,i5,1x,1pe17.10,2(1x,i2))') 'MD ',&
            istep,e_rxyz,nummax,nummin
       if (iproc == 0) write(*,'(a,i5,1x,1pe17.10,2(1x,i2))') ' #MH MD ',&
            istep,e_rxyz,nummax,nummin
       if (nummin.ge.mdmin) then
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
      do iat=1,atoms%nat
      write(1000+iproc,*) "IN  ",wpos(:,iat)
      enddo
   if (atoms%geocode == 'S') then 
      call fixfrag_posvel_slab(iproc,atoms%nat,rcov,wpos,vxyz,2)
      do iat=1,atoms%nat
      write(1000+iproc,*) "OUT ",wpos(:,iat)
      enddo
      flush(1000+iproc)
   else if (atoms%geocode == 'F') then
     if (istep == istepnext) then 
           call fixfrag_posvel(iproc,atoms%nat,rcov,rxyz,vxyz,2,occured)
        if (occured) then 
          istepnext=istep+4
        else
          istepnext=istep+1
        endif
     endif
   endif
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
    use m_ab6_symmetry
    implicit real*8 (a-h,o-z)
    type(atoms_data) :: atoms
    dimension fxyz(3*atoms%nat),gg(3*atoms%nat),vxyz(3*atoms%nat),rxyz(3*atoms%nat),rxyz_old(3*atoms%nat)
    type(input_variables) :: inputs_md
    type(restart_objects) :: rst
    !Local variables
    dimension wpos(3*atoms%nat),strten(6)

!    eps_vxyz=1.d-1*atoms%nat
    alpha=inputs_md%betax

    !allocate(wpos(3,nat),fxyz(3,nat))

    inputs_md%inputPsiId=1
    if(iproc==0)write(*,*)'#MH soften initial step '
    call call_bigdft(nproc,iproc,atoms,rxyz,inputs_md,etot0,fxyz,strten,fnoise,rst,infocode)

    ! scale velocity to generate dimer 

!    call atomic_dot(atoms,vxyz,vxyz,svxyz)
    svxyz=0.d0
    do i=1,3*atoms%nat
       iat=(i-1)/3+1
       if (atoms%ifrztyp(iat) == 0) then
          svxyz=svxyz+vxyz(i)**2
       end if
    enddo
    eps_vxyz=sqrt(svxyz)
    if(iproc == 0) write(*,*)'#MH  eps_vxyz=',eps_vxyz

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
!      call atomic_axpy(atoms,rxyz,1.d0,vxyz,wpos)
       wpos=rxyz+vxyz
       call call_bigdft(nproc,iproc,atoms,wpos,inputs_md,etot,fxyz,strten,fnoise,rst,infocode)
       fd2=2.d0*(etot-etot0)/eps_vxyz**2

       sdf=0.d0
       svxyz=0.d0
       do i=1,3*atoms%nat
          iat=(i-1)/3+1
          if (atoms%ifrztyp(iat) == 0) then
             sdf=sdf+vxyz(i)*fxyz(i)
             svxyz=svxyz+vxyz(i)*vxyz(i)
          end if
       end do
!       call atomic_dot(atoms,vxyz,vxyz,svxyz)
!       call atomic_dot(atoms,vxyz,fxyz,sdf)

       curv=-sdf/svxyz
       if (it == 1) curv0=curv

       res=0.d0
       do i=1,3*atoms%nat
          iat=(i-1)/3+1
          if (atoms%ifrztyp(iat) == 0) then
             fxyz(i)=fxyz(i)+curv*vxyz(i)
             res=res+fxyz(i)**2
          end if
       end do
!       call atomic_axpy_forces(atoms,fxyz,curv,vxyz,fxyz)
!       call atomic_dot(atoms,fxyz,fxyz,res)
       res=sqrt(res)

       write(fn4,'(i4.4)') it
       write(comment,'(a,1pe10.3)')'res= ',res
       if (iproc == 0) &
            call write_atomic_file(trim(inputs_md%dir_output)//'possoft_'//fn4,&
            etot,wpos,atoms,trim(comment),forces=fxyz)

       if(iproc==0)write(*,'(a,i3,5(f12.5),f10.3)')'#MH soften it, curv, fd2,dE,res,eps_vxyz:',&
            it, curv, fd2,etot-etot0,res,eps_vxyz
       if (curv.lt.0.d0 .or. fd2.lt.0.d0) then
          if(iproc==0) write(*,*) '#MH NEGATIVE CURVATURE'
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
!       call atomic_axpy_forces(atoms,wpos,alpha,fxyz,wpos)
        call daxpy(3*atoms%nat,alpha,fxyz(1),1,wpos(1),1)
       do i=1,3*atoms%nat
          vxyz(i)=wpos(i)-rxyz(i)
       end do
       write(comment,'(a,1pe10.3)')'curv= ',curv
       if (iproc == 0) &
            call write_atomic_file(trim(inputs_md%dir_output)//'posvxyz',0.d0,vxyz,atoms,trim(comment),forces=fxyz)
       call elim_moment(atoms%nat,vxyz)
       call elim_torque_reza(atoms%nat,rxyz,vxyz)

       svxyz=0.d0
       do i=1,3*atoms%nat
          iat=(i-1)/3+1
          if (atoms%ifrztyp(iat) == 0) then
             svxyz=svxyz+vxyz(i)*vxyz(i)
          end if
       end do
!      call atomic_dot(atoms,vxyz,vxyz,svxyz)
       if (res <= curv*eps_vxyz*5.d-1) exit
       svxyz=eps_vxyz/dsqrt(svxyz)

       do i=1,3*atoms%nat
          vxyz(i)=vxyz(i)*svxyz
       end do

    end do ! iter
 

    !        deallocate(wpos,fxyz)
  END SUBROUTINE soften

end program MINHOP


!> Inserts the energy re_wpos at position k_e_wpos and shifts up all other energies
subroutine insert(nlminx,nlmin,k_e_wpos,re_wpos,earr)
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


!> C save configuration if it is among the lowest ones in energy
subroutine save_low_conf(nat,npmin,npminx,e_wpos,pos,elocmin,poslocmin)
  implicit real*8 (a-h,o-z)
  dimension elocmin(npminx)
  dimension pos(3,nat),poslocmin(3,nat,npminx)

  if (npmin.le.npminx) then
     kmax=npmin
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


!> C x is in interval [xx(jlo),xx(jlow+1)[ ; xx(0)=-Infinity ; xx(n+1) = Infinity
subroutine hunt_g(xx,n,x,jlo)
  implicit none
  !Arguments
  integer :: jlo,n
  real(kind=8) :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  if (n.le.0) stop 'hunt_g'
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
END SUBROUTINE hunt_g


!>  assigns initial velocities for the MD escape part
subroutine velnorm(at,rxyz,ekinetic,vxyz)
  use module_base
  use module_types
  use m_ab6_symmetry
  implicit none
  !implicit real*8 (a-h,o-z)
  real(gp), intent(in) :: ekinetic
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(inout) :: vxyz
  !local variables
  integer :: iat
  real(gp) :: rkin,rkinsum,sclvel

  !C      Kinetic energy of the initial velocities
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

END SUBROUTINE velnorm


!> create a random displacement vector without translational and angular moment
subroutine randdist(nat,rxyz,vxyz)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(3*nat), intent(in) :: rxyz
  real(gp), dimension(3*nat), intent(out) :: vxyz
  !local variables
  integer :: i,idum=0
  real(kind=4) :: tt,builtin_rand
  do i=1,3*nat
     !call random_number(tt)
     !add built-in random number generator
     tt=builtin_rand(idum)
     vxyz(i)=real(tt-.5,gp)*3.e-1_gp
  end do

  call elim_moment(nat,vxyz)
  call elim_torque_reza(nat,rxyz,vxyz)
END SUBROUTINE randdist


!>  generates 3*nat random numbers distributed according to  exp(-.5*vxyz**2)
subroutine gausdist(nat,rxyz,vxyz)
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
  call  elim_torque_reza(nat,rxyz,vxyz)
  return
END SUBROUTINE gausdist


!>  generates n random numbers distributed according to  exp(-x)
subroutine expdist(nat,rxyz,vxyz)
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
  call  elim_torque_reza(nat,rxyz,vxyz)

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
  use module_base, only: gp
  implicit real*8 (a-h,o-z)
  dimension rxyz(3,nat),vxyz(3,nat)

  ! center of mass
  cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
  do iat=1,nat
     cmx=cmx+rxyz(1,iat)
     cmy=cmy+rxyz(2,iat)
     cmz=cmz+rxyz(3,iat)
  enddo
  cmx=cmx/real(nat,gp) 
  cmy=cmy/real(nat,gp) 
  cmz=cmz/real(nat,gp)

  ! torque
  tx=0.d0 ; ty=0.d0 ; tz=0.d0
  do iat=1,nat
     tx=tx+(rxyz(2,iat)-cmy)*vxyz(3,iat)-(rxyz(3,iat)-cmz)*vxyz(2,iat)
     ty=ty+(rxyz(3,iat)-cmz)*vxyz(1,iat)-(rxyz(1,iat)-cmx)*vxyz(3,iat)
     tz=tz+(rxyz(1,iat)-cmx)*vxyz(2,iat)-(rxyz(2,iat)-cmy)*vxyz(1,iat)
  enddo
  write(*,'(a,3(1pe11.3))') '#MH torque',tx,ty,tz

END SUBROUTINE torque


!subroutine elim_torque(nat,rxyz,vxyz)
!  implicit real*8 (a-h,o-z)
!  dimension rxyz(3,nat),vxyz(3,nat),t(3)
!
!  ! center of mass
!  cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
!  do iat=1,nat
!     cmx=cmx+rxyz(1,iat)
!     cmy=cmy+rxyz(2,iat)
!     cmz=cmz+rxyz(3,iat)
!  enddo
!  cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat
!
!  do it=1,100
!
!     ! torque and radii in planes
!     t(1)=0.d0 ; t(2)=0.d0 ; t(3)=0.d0
!     sx=0.d0 ; sy=0.d0 ; sz=0.d0
!     do iat=1,nat
!        t(1)=t(1)+(rxyz(2,iat)-cmy)*vxyz(3,iat)-(rxyz(3,iat)-cmz)*vxyz(2,iat)
!        t(2)=t(2)+(rxyz(3,iat)-cmz)*vxyz(1,iat)-(rxyz(1,iat)-cmx)*vxyz(3,iat)
!        t(3)=t(3)+(rxyz(1,iat)-cmx)*vxyz(2,iat)-(rxyz(2,iat)-cmy)*vxyz(1,iat)
!        sx=sx+(rxyz(1,iat)-cmx)**2
!        sy=sy+(rxyz(2,iat)-cmy)**2
!        sz=sz+(rxyz(3,iat)-cmz)**2
!     enddo
!
!     if (t(1)**2+t(2)**2+t(3)**2.lt.1.d-22) return
!
!     ii=0
!     tmax=0.d0
!     do i=1,3
!        if (t(i)**2.gt.tmax**2) then 
!           ii=i
!           tmax=t(i)
!        endif
!     enddo
!
!     !         write(*,'(i4,3(1pe11.3))') ii,t
!
!     ! modify velocities
!     if (ii == 1) then 
!        cx=t(1)/(sz+sy)
!        do iat=1,nat
!           vxyz(2,iat)=vxyz(2,iat)+cx*(rxyz(3,iat)-cmz)
!           vxyz(3,iat)=vxyz(3,iat)-cx*(rxyz(2,iat)-cmy)
!        enddo
!     else if(ii == 2) then 
!        cy=t(2)/(sz+sx)
!        do iat=1,nat
!           vxyz(1,iat)=vxyz(1,iat)-cy*(rxyz(3,iat)-cmz)
!           vxyz(3,iat)=vxyz(3,iat)+cy*(rxyz(1,iat)-cmx)
!        enddo
!     else if(ii == 3) then 
!        cz=t(3)/(sy+sx)
!        do iat=1,nat
!           vxyz(1,iat)=vxyz(1,iat)+cz*(rxyz(2,iat)-cmy)
!           vxyz(2,iat)=vxyz(2,iat)-cz*(rxyz(1,iat)-cmx)
!        enddo
!     else
!        stop 'wrong ii'
!     endif
!
!  enddo
!  write(*,'(a,3(1pe11.3))') 'WARNING REMAINING TORQUE',t
!
!END SUBROUTINE elim_torque


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


subroutine winter(at,re_pos,pos,npminx,nlminx,nlmin,npmin,accur, & 
     earr,elocmin,poslocmin,eref,ediff,ekinetic,dt,nsoften)
  use module_base
  use module_types
  use module_interfaces
  use m_ab6_symmetry
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: npminx,nlminx,nlmin,npmin,nsoften
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

     call wtpos(at,npminx,nlminx,nlmin,npmin,poslocmin,earr,elocmin)

END SUBROUTINE winter


subroutine wtioput(ediff,ekinetic,dt,nsoften)
  implicit real*8 (a-h,o-z)
  open(unit=11,file='ioput',status='unknown')
  write(11,'(3(1x,1pe24.17)1x,i4,a)') ediff,ekinetic,dt,nsoften,' ediff, ekinetic dt and nsoften'
  close(11)
END SUBROUTINE wtioput


subroutine wtpos(at,npminx,nlminx,nlmin,npmin,pos,earr,elocmin)
  use module_base
  use module_types
  use module_interfaces
  use m_ab6_symmetry
  implicit none
  !implicit real*8 (a-h,o-z)
  integer, intent(in) :: npminx,nlminx,nlmin,npmin
  type(atoms_data), intent(in) :: at
  real(gp), dimension(npminx), intent(in) :: elocmin
  real(gp), dimension(0:nlminx,2), intent(in) :: earr
  real(gp), dimension(3,at%nat,npminx), intent(in) :: pos
  !local variables
  character(len=5) :: fn
  character(len=17) :: filename
  integer :: k,kk,i

       write(*,*) 'nlmin,nlminx,npmin,npminx',nlmin,nlminx,npmin,npminx
       do i=1,min(40,nlmin,nlminx)
         write(*,'(i4,e24.17)') i,earr(i,1)
       enddo

  do k=1,min(npmin,npminx)
             write(*,'(a,i4,e24.17)') 'k,elocmin(k)',k,elocmin(k)


     !C Classify the configuration in the global ranking
     kk=0
     find_kk : do
        kk=kk+1
        if (kk > min(nlmin,nlminx)) then 
           write(*,*) 'ranking error for',k
           stop 
        endif
        if (earr(kk,1) == elocmin(k)) exit find_kk
!        if (abs(earr(kk,1) - elocmin(k)) .lt. 1.d-12 ) then 
!             write(*,*) 'match ',abs(earr(kk,1) - elocmin(k))
!             exit find_kk
!        endif
     end do find_kk

     if (kk <= npminx) then

                write(*,'(a,i4,i4,1x,1pe21.14)') 'k,kk,elocmin(k)',k,kk,elocmin(k)

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
!  ! write(*,*)'#MH reading unrelaxed structure from ',filename
!  open(unit=9,file=filename,status='old')
!  read(9,*)fn!no need for header
!  read(9,*)fn!same
!  do iat=1,nat
!     read(9,*)fn,rxyz(:,iat)!we know the atom types already
!  enddo
!  close(unit=9)
!END SUBROUTINE rdposout


!> routine for adjusting the dimensions with the center of mass in the middle
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


! subroutine fix_fragmentation(iproc,at,rxyz,nputback)
!   use module_base
!   use module_types
!   use m_ab6_symmetry
!   implicit none
!   !implicit real*8 (a-h,o-z)
!   integer, intent(in) :: iproc
!   type(atoms_data), intent(in) :: at
!   integer, intent(inout) :: nputback
!   real(gp), dimension(3,at%nat) :: rxyz
!   !local variables
!   real(gp), parameter :: bondlength=8.0_gp
!   integer :: iat,nloop,ncluster,ii,jat,jj,kat,nadd,ierr
!   real(gp) :: xi,yi,zi,xj,yj,zj,ddmin,dd,d1,d2,d3,tt
!   ! automatic arrays
!   logical, dimension(at%nat) :: belong
! 
!   nloop=1
! 
!   fragment_loop: do
! 
!      iat=1
!      belong(iat)=.true.
!      ncluster=1
!      do iat=2,at%nat
!         belong(iat)=.false.
!      enddo
! 
!      !   ic=0
!      form_cluster: do
!         nadd=0
!         do iat=1,at%nat
!            xi=rxyz(1,iat) 
!            yi=rxyz(2,iat) 
!            zi=rxyz(3,iat)
!            if (belong(iat)) then 
!               do jat=1,at%nat
!                  xj=rxyz(1,jat) ; yj=rxyz(2,jat) ; zj=rxyz(3,jat)
!                  if ( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 <= (bondlength*1.25d0)**2) then 
!                     if (.not. belong(jat)) nadd=nadd+1
!                     belong(jat)=.true. 
!                  endif
!               end do
!            endif
!         end do
!         ncluster=ncluster+nadd
!         !     ic=ic+1 ; write(*,*) 'nadd,ncluster',ic,nadd,ncluster
!         if (nadd == 0) exit form_cluster
!      enddo form_cluster
! 
!      if (ncluster == at%nat) then 
!         !   write(*,*) 'No fragmentation has occured',nloop
!         return
! 
!      else
!         nputback=nputback+1
! 
!         if (iproc == 0) then
!            write(*,*) '#MH fragmentation occured',nloop,ncluster
!            write(444,*) at%nat,ncluster
!            write(444,*) ' fragmented configuration ', nputback
!            do kat=1,at%nat
!               write(444,*) ' LJ  ',rxyz(1,kat),rxyz(2,kat),rxyz(3,kat)
!            enddo
!         endif
! 
! 
!         ! make sure the part that flew away is smaller than the cluster
!         if (ncluster <= at%nat/2) then
!            !     write(*,*) 'FLIP'
!            do iat=1,at%nat
!               belong(iat)=.not. belong(iat)
!            enddo
!         endif
! 
!         ! pull back the fragment of atoms that flew away
!         ii=-99999
!         do iat=1,at%nat
!            if (.not. belong(iat)) then
!               xi=rxyz(1,iat) 
!               yi=rxyz(2,iat) 
!               zi=rxyz(3,iat)
!               ddmin=1.e100_gp
!               jj=-99999
!               do jat=1,at%nat
!                  if (belong(jat)) then
!                     xj=rxyz(1,jat) 
!                     yj=rxyz(2,jat) 
!                     zj=rxyz(3,jat)
!                     dd= (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 
!                     if (dd < ddmin) then 
!                        jj=jat
!                        ii=iat
!                        ddmin=dd
!                     endif
!                  endif
!               enddo
!            endif
!         enddo
! 
!         d1=rxyz(1,ii)-rxyz(1,jj)
!         d2=rxyz(2,ii)-rxyz(2,jj)
!         d3=rxyz(3,ii)-rxyz(3,jj)
!         tt=bondlength/sqrt(d1**2+d2**2+d3**2)
!         do iat=1,at%nat
!            if (.not. belong(iat) ) then  !.and. .not. at%lfrztyp(iat)) then
!               if (at%geocode == 'P') then
! stop  '------ P ----------'
!                  rxyz(1,iat)=modulo(rxyz(1,iat)-d1*(tt),at%alat1)
!                  rxyz(2,iat)=modulo(rxyz(2,iat)-d2*(tt),at%alat2)
!                  rxyz(3,iat)=modulo(rxyz(3,iat)-d3*(tt),at%alat3)
!               else if (at%geocode == 'S') then
! stop  '------ S ----------'
!                  rxyz(1,iat)=modulo(rxyz(1,iat)-d1*(tt),at%alat1)
!                  rxyz(2,iat)=       rxyz(2,iat)-d2*(tt)
!                  rxyz(3,iat)=modulo(rxyz(3,iat)-d3*(tt),at%alat3)
!               else
!                  rxyz(1,iat)=rxyz(1,iat)-d1*(tt)
!                  rxyz(2,iat)=rxyz(2,iat)-d2*(tt)
!                  rxyz(3,iat)=rxyz(3,iat)-d3*(tt)
!               end if
!            endif
!         enddo
! 
!         if (iproc == 0) then
!            write(444,*) at%nat, 'atomic ' 
!            write(444,*) ' fixed configuration ', nputback,sqrt(d1**2+d2**2+d3**2),ii,jj
!            do iat=1,at%nat
!               write(444,'(a5,3(e15.7),l1)') ' LJ  ',rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),belong(iat)
!            enddo
!         endif
!         nloop=nloop+1
!      if (nloop.gt.4) then 
!           write(*,*)"#MH fragmentation could not be fixed",nloop
!           call MPI_ABORT(MPI_COMM_WORLD,ierr)
!      endif
!      endif
!   end do fragment_loop
! 
! END SUBROUTINE fix_fragmentation



!    implicit real*8 (a-h,o-z)
!    integer option
!    logical occured
!    character(5) atomname
!    parameter(nat=12,iproc=0,option=1)
!    dimension pos(3,nat),vel(3,nat)
!
!
!    open(unit=1,file='T.xyz')
!    read(1,*) natp
!    if (nat .ne. natp) stop' natp'
!    read(1,*) 
!    do iat=1,nat
!    read(1,*) atomname,(pos(i,iat),i=1,3)
!    enddo
!    close(1)
! 
!
!    open(unit=444,file='t.xyz')
!    call fixfrag_posvel(iproc,nat,pos,vel,1,occured)
!    write(*,*) 'occured',occured
!    close(444)
!   
!    end
!


subroutine fixfrag_posvel(iproc,nat,rcov,pos,vel,option,occured)
!This subroutine can perform two tasks.
!ATTENTION: it will only work on free BC!!!
!
!option=1
!The atoms in pos are analyzed and, if there is a fragmentation occuring, the
!main fragment will be identified and all neighboring fragments will be moved towards the nearest
!atom of the main fragment. The array pos will then be updated and returned.
!
!option=2
!The fragments are identified and the center of mass of all fragments are computed separately.
!The center of mass of all cluster is also computed.
!Then, the velocities are modified in such a way that the projection of the velocities 
!along the vector pointing towards the center of mass of all fragments are inverted 
!!use module_base
!!use module_types
!!use m_ab6_symmetry
implicit none
integer, intent(in) :: iproc,nat
!type(atoms_data), intent(in) :: at
real(8),dimension(3,nat), INTENT(INOUT) :: pos
real(8),dimension(3,nat), INTENT(INOUT) :: vel
real(8),dimension(nat), INTENT(IN) :: rcov
integer, INTENT(IN):: option
integer :: nfrag, nfragold
logical :: occured,niter
real(8)::  dist, mindist, angle, vec(3), cmass(3), velcm(3), bondlength, bfactor,rnrmi,scpr
real(8):: ekin,vcm1,vcm2,vcm3,ekin0,scale
real(8), allocatable:: cm_frags(:,:), vel_frags(:,:)
integer::iat, jat, nmax(1), imin(2),ifrag
integer, allocatable:: fragcount(:)
integer, allocatable:: nat_frags(:)
integer, dimension(nat):: fragarr
logical, allocatable:: invert(:)

!The bondlength (in atomic units) is read from file input.bondcut
! OPTION 1: System is considered to be fragmented if the minimal distance between two atoms in the fragment is more than 2.0*bondlength
!           The two fragment are then brought together such that the minimal distance equals 1.5*bondlength
! OPTION  : System is considered to be fragmented if the minimal distance between two atoms in the fragment is more than 2.0*bondlength
!           the velocities are then inverted
!open(unit=43,file="input.bondcut")
!read(43,*) bondlength
!close(43)

if (option == 1) then 
   bfactor=1.5d0
else if (option == 2) then 
   bfactor=2.d0
else
   stop 'wrong option'
endif



fragarr(:)=0                     !Array, which atom belongs to which fragment
nfrag=0                       !Number of fragments

!Check for correct input
if (option.ne.1 .and. option.ne.2) stop "Wrong option in fixfrag_refvels"

!Calculate number of fragments and fragmentlist of the atoms
do
   nfragold=nfrag
   do iat=1,nat                !Check the first atom that isn't part of a cluster yet
      if(fragarr(iat)==0) then
         nfrag=nfrag+1
         fragarr(iat)=nfrag
         exit
      endif
   enddo
   if (nfragold==nfrag) exit
7000 niter=.false.
   do iat=1,nat                !Check if all the other atoms are part of the current cluster
      do jat=1,nat
         bondlength=rcov(iat)+rcov(jat)
         if(nfrag==fragarr(iat) .AND. jat.ne.iat .AND. fragarr(jat)==0) then
            dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
            if(dist<(bfactor*bondlength)**2) then
               fragarr(jat)=nfrag
               niter=.true.
            endif
         endif
      enddo
   enddo
   if(niter) then
      goto 7000
   endif
enddo


!   if(iproc==0) write(*,*) '#MH nfrag=',nfrag
occured=.false.
if(nfrag.ne.1) then          !"if there is fragmentation..."
   occured=.true.
   if(iproc==0) write(*,*) '#MH FIX: Number of Fragments counted with option', nfrag,option

   if (option==1) then !OPTION=1, FIX FRAGMENTATION
      !   if(nfrag.ne.1) then          !"if there is fragmentation..."

      !Find out which fragment is the main cluster
      allocate(fragcount(nfrag))
      fragcount=0
      do ifrag=1,nfrag
         do iat=1,nat
            if(fragarr(iat)==ifrag) then
               fragcount(ifrag)=fragcount(ifrag)+1
            endif
         enddo
      enddo
      nmax=maxloc(fragcount(:))
      if(iproc==0) write(*,*) '#MH FIX: The main Fragment index is', nmax(1)

      !Find the minimum distance between the clusters
      do ifrag=1,nfrag
         mindist=1.d100
         if(ifrag.ne.nmax(1)) then
            do iat=1,nat
               if(fragarr(iat)==ifrag) then
                  do jat=1,nat
                     if(fragarr(jat)==nmax(1)) then
                        dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
                        if(dist<mindist**2) then
                           mindist=sqrt(dist)
                           imin(1)=jat  !Atom with minimal distance in main fragment
                           imin(2)=iat   !Atom with minimal distance in fragment ifrag
                        endif
                     endif
                  enddo
               endif
            enddo

            if (iproc == 0) then
               write(444,*) nat, 'atomic '
               write(444,*) 'A fragmented configuration ',imin(1),imin(2)
               do iat=1,nat
                  write(444,'(a5,3(e15.7),l1)') ' Mg  ',pos(1,iat),pos(2,iat),pos(3,iat)
               enddo
            endif


            vec(:)=pos(:,imin(1))-pos(:,imin(2))
            bondlength=rcov(imin(1))+rcov(imin(2))
            do iat=1,nat        !Move fragments back towards the main fragment 
               if(fragarr(iat)==ifrag) then
                  pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.5d0*bondlength)/mindist)
                  fragarr(iat)=nmax(1)
               endif
            enddo


         endif
      enddo
      deallocate(fragcount)
      if(iproc==0) write(*,*) '#MH FIX: Fragmentation fixed! Keep on hopping...'
      if (iproc == 0) then
         write(444,*) nat, 'atomic '
         write(444,*) ' fixed configuration '
         do iat=1,nat
            write(444,'(a5,3(e15.7),l1)') ' Mg  ',pos(1,iat),pos(2,iat),pos(3,iat)
         enddo
      endif

      !   endif
   elseif(option==2) then !OPTION=2, INVERT VELOCITIES
      !   if(nfrag.ne.1) then          !"if there is fragmentation..."
      if(iproc==0) write(*,*) "#MH FIX: Preparing to invert velocities, option:",option
      !Compute center of mass of all fragments and the collectiove velocity of each fragment
      allocate(cm_frags(3,nfrag),vel_frags(3,nfrag),nat_frags(nfrag))
      allocate(invert(nfrag))
      cm_frags(:,:)=0.d0
      vel_frags(:,:)=0.d0
      nat_frags(:)=0         !number of atoms per fragment
      cmass(:)=0.d0
      velcm(:)=0.d0
      do iat=1,nat
         ifrag=fragarr(iat)
         nat_frags(ifrag)=nat_frags(ifrag)+1
         cm_frags(:,ifrag)=cm_frags(:,ifrag)+pos(:,iat)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)+vel(:,iat)
      enddo

      do ifrag=1,nfrag
         cm_frags(:,ifrag)=cm_frags(:,ifrag)/real(nat_frags(ifrag),8)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)/real(nat_frags(ifrag),8)
         cmass(:)=cmass(:)+cm_frags(:,ifrag)*nat_frags(ifrag)/real(nat,8)
         velcm(:)=velcm(:)+vel_frags(:,ifrag)*nat_frags(ifrag)/real(nat,8)
      enddo
      if (iproc==0) write(*,*) '#MH CM VELOCITY',sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
      if (velcm(1)**2+velcm(2)**2+velcm(3)**2.gt.1.d-24) then
         if (iproc==0) write(*,*) '#MH NONZERO CM VELOCITY'
      endif


      ! now cm_frags contains the unit vector pointing from the center of mass of the entire system to the center of mass of the fragment
      do ifrag=1,nfrag
         cm_frags(:,ifrag)=cm_frags(:,ifrag)-cmass(:)
         rnrmi=1.d0/sqrt(cm_frags(1,ifrag)**2+cm_frags(2,ifrag)**2+cm_frags(3,ifrag)**2)
         cm_frags(1,ifrag)=cm_frags(1,ifrag)*rnrmi
         cm_frags(2,ifrag)=cm_frags(2,ifrag)*rnrmi
         cm_frags(3,ifrag)=cm_frags(3,ifrag)*rnrmi
         angle=cm_frags(1,ifrag)*vel_frags(1,ifrag)+cm_frags(2,ifrag)*vel_frags(2,ifrag)+cm_frags(3,ifrag)*vel_frags(3,ifrag)
         rnrmi=1.d0/sqrt(vel_frags(1,ifrag)**2+vel_frags(2,ifrag)**2+vel_frags(3,ifrag)**2)
         angle=angle*rnrmi
         if (angle.gt.0.d0) then
            invert(ifrag)=.true.
         else
            invert(ifrag)=.false.
         endif
         if (iproc==0) write(*,*) '#MH ifrag, angle ',ifrag, angle,invert(ifrag)
      enddo
      !Decompose each atomic velocity into an component parallel and perpendicular to the cm_frags  vector and inter the 
      !paralle part if it point away from the CM

      !Check kinetic energy before inversion
      ekin0=0.d0
      vcm1=0.d0
      vcm2=0.d0
      vcm3=0.d0
      do iat=1,nat
         ekin0=ekin0+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
         vcm1=vcm1+vel(1,iat)
         vcm2=vcm2+vel(2,iat)
         vcm3=vcm3+vel(3,iat)
      enddo
      if (iproc==0) write(*,'(a,e14.7,3(e10.3))') '#MH EKIN CM before invert',ekin0,vcm1,vcm2,vcm3
      if (iproc==0) call torque(nat,pos,vel)
      !Checkend kinetic energy before inversion

      do iat=1,nat
         ! inversions  by fragment group
         ifrag=fragarr(iat)
         if (invert(ifrag)) then
            scpr=cm_frags(1,ifrag)*vel(1,iat)+cm_frags(2,ifrag)*vel(2,iat)+cm_frags(3,ifrag)*vel(3,iat)
            vel(:,iat)=vel(:,iat)-scpr*cm_frags(:,ifrag)*2.d0
         endif
      enddo

      call elim_moment(nat,vel)
      call elim_torque_reza(nat,pos,vel)

      ! scale velocities to regain initial ekin0
      ekin=0.d0
      do iat=1,nat
         ekin=ekin+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
      enddo
      scale=sqrt(ekin0/ekin)
      do iat=1,nat
         vel(1,iat)=vel(1,iat)*scale
         vel(2,iat)=vel(2,iat)*scale
         vel(3,iat)=vel(3,iat)*scale
      enddo

      !Check kinetic energy after inversion
      ekin=0.d0
      vcm1=0.d0
      vcm2=0.d0
      vcm3=0.d0
      do iat=1,nat
         ekin=ekin+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
         vcm1=vcm1+vel(1,iat)
         vcm2=vcm2+vel(2,iat)
         vcm3=vcm3+vel(3,iat)
      enddo
      if (iproc==0) write(*,'(a,e14.7,3(e10.3))') '#MH EKIN CM after  invert',ekin,vcm1,vcm2,vcm3
      if (iproc==0) call torque(nat,pos,vel)
      !Checkend kinetic energy after inversion

      !Check angle  after inversion
      vel_frags(:,:)=0.d0
      do iat=1,nat
         ifrag=fragarr(iat)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)+vel(:,iat)
      enddo
      do ifrag=1,nfrag
         angle=cm_frags(1,ifrag)*vel_frags(1,ifrag)+cm_frags(2,ifrag)*vel_frags(2,ifrag)+cm_frags(3,ifrag)*vel_frags(3,ifrag)
         rnrmi=1.d0/sqrt(vel_frags(1,ifrag)**2+vel_frags(2,ifrag)**2+vel_frags(3,ifrag)**2)
         angle=angle*rnrmi
         if (iproc==0) write(*,*) '#MH ifrag, angle a invert',ifrag, angle
      enddo
      !Checkend kinetic energy after inversion


      deallocate(cm_frags,vel_frags,nat_frags)
      deallocate(invert)
      !   endif
      !else
      !   stop "Wrong option within ff-rv"
      if(iproc==0) write(*,*) "#MH FIX: Velocity component towards the center of mass inverted! Keep on hopping..."
   endif
endif
end subroutine fixfrag_posvel



!implicit real(8) (a-h,o-z)
!parameter(nat=36)
!dimension rcov(nat), pos(3,nat),vel(3,nat)
!character*2 atn(nat)
!
!    iproc=0
!    open(unit=1,file='p.xyz')
!    read(1,*) 
!    read(1,*) 
!    do iat=1,nat
!    read(1,*) atn(iat),(pos(i,iat),i=1,3)
!    if (atn(iat).eq."Ca") then 
!    rcov(iat)=3.d0
!    else
!    rcov(iat)=1.35d0
!    endif
!    enddo
!    write(*,*) rcov
!
!  call fixfrag_posvel_slab(iproc,nat,rcov,pos,vel,2)
!  call fixfrag_posvel_slab(iproc,nat,rcov,pos,vel,1)
!
!    open(unit=2,file='t.xyz')
!    write(2,*) nat
!    write(2,*) 
!    do iat=1,nat
!    write(2,*) atn(iat),(pos(i,iat),i=1,3)
!    enddo
!
!  end
!

subroutine fixfrag_posvel_slab(iproc,nat,rcov,pos,vel,option)
!This subroutine points the velocities towards the surface if an atom is too far away from the surface with surface boundary conditions
!
implicit none
integer, intent(in) :: iproc,nat,option
!type(atoms_data), intent(in) :: at
real(8),dimension(3,nat), INTENT(INOUT) :: pos
real(8),dimension(3,nat), INTENT(INOUT) :: vel
real(8),dimension(nat), INTENT(IN) :: rcov
integer :: iat,i,ic,ib,ilow,ihigh,icen,mm,mj,jat
real(8) :: ymin, ylow,yhigh,dx,dy,dz,dl,dist,distmin,d

integer, dimension(-100:1000):: ygrid
logical ,dimension(nat) :: onsurface


! empty space = 0
    do i=-100,1000 
    ygrid(i)=0
    enddo

    ymin=1.d100 
    do iat=1,nat
        ymin=min(ymin,pos(2,iat)) 
    enddo

! occupied space= nonzero
    do iat=1,nat
        ic=nint((pos(2,iat)-ymin)*4.d0)  ! ygrid spacing=.25
         ib=2.0d0*rcov(iat)*4.d0
         if (ic-ib.lt.-100) stop "#MH error fixfrag_slab -100"
         if (ic+ib.gt.1000) stop "#MH error fixfrag_slab 1000"
         do i=ic-ib,ic+ib
         ygrid(i)=ygrid(i)+1
         enddo
    enddo

! find center of slab
    mm=0
    do i=-100,1000
    if (ygrid(i) .gt. mm) then
        icen=i
        mm=ygrid(i)
    endif
    enddo

! find border between empty and occupied space
    do i=icen,-100,-1
    if (ygrid(i).eq.0) then
        ilow=i
        exit
    endif
    enddo

    do i=icen,1000
    if (ygrid(i).eq.0) then
        ihigh=i
        exit
    endif
    enddo


    ylow=ymin+ilow*.25d0
    yhigh=ymin+ihigh*.25d0
    if (iproc.eq.0) write(*,*) "#MH ylow,ycen,yhigh",ylow,ymin+icen*.25d0,yhigh
             write(1000+iproc,*) "#MH ylow,ycen,yhigh",ylow,ymin+icen*.25d0,yhigh

if (option.eq.2) then

    do iat=1,nat
         if (pos(2,iat).lt.ylow-rcov(iat)) then 
             vel(2,iat)=abs(vel(2,iat))
             if (iproc.eq.0) write(*,*) "#MH velocity made positive for atom",iat
             write(1000+iproc,*) "#MH velocity made positive for atom",iat,pos(:,iat)
         endif
         if (pos(2,iat).gt.yhigh+rcov(iat)) then 
             vel(2,iat)=-abs(vel(2,iat))
             if (iproc.eq.0) write(*,*) "#MH velocity made negative for atom",iat
             write(1000+iproc,*) "#MH velocity made negative for atom",iat,pos(:,iat)
         endif
    enddo
             flush(1000+iproc) 

else if (option.eq.1) then
1000 continue
    do iat=1,nat
         if (pos(2,iat).lt.ylow-rcov(iat) .or. pos(2,iat).gt.yhigh+rcov(iat)) then 
         onsurface(iat)=.false.
         else
         onsurface(iat)=.true.
         endif
    enddo
    do iat=1,nat
         if (onsurface(iat) .eqv. .false.) then 
             distmin=1.d100
            do jat=1,nat
            if (jat.ne.iat .and. onsurface(jat)) then
              dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
              dist=sqrt(dist)-1.25d0*rcov(iat)-1.25d0*rcov(jat)
              if (dist.lt.distmin) then 
                distmin=dist
                mj=jat
              endif
            endif
            enddo
            if (iproc.eq.0) write(*,*) iat,mj,distmin
            if (distmin.gt.0.d0) then
                dx=pos(1,iat)-pos(1,mj)
                dy=pos(2,iat)-pos(2,mj)
                dz=pos(3,iat)-pos(3,mj)
                dl=sqrt(dx**2+dy**2+dz**2)
                d=distmin+0.1d0*(rcov(iat)+rcov(mj))
                dx=dx*(d/dl)
                dy=dy*(d/dl)
                dz=dz*(d/dl)
                if (iproc.eq.0) write(*,*) "#MH moving atom",iat,pos(:,iat)
                pos(1,iat)=pos(1,iat)-dx
                pos(2,iat)=pos(2,iat)-dy
                pos(3,iat)=pos(3,iat)-dz
                if (iproc.eq.0) write(*,*) "#MH moved atom",iat,pos(:,iat)
                onsurface(iat)=.true.
                goto 1000
            endif
         endif
    enddo
else 
    stop "invalid option for fixfrag_slab"
endif

end subroutine fixfrag_posvel_slab




subroutine give_rcov(iproc,atoms,nat,rcov)
  !    use module_base
  use module_types
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nat
  type(atoms_data), intent(in) :: atoms
  real(kind=8), intent(out) :: rcov(nat)
  !Local variables
  integer :: iat

  do iat=1,nat
     if (trim(atoms%atomnames(atoms%iatype(iat)))=='H') then
        rcov(iat)=0.75d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='He') then
        rcov(iat)=0.75d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Li') then
        rcov(iat)=3.40d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Be') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='B' ) then
        rcov(iat)=1.55d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='C' ) then
        rcov(iat)=1.45d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='N' ) then
        rcov(iat)=1.42d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='O' ) then
        rcov(iat)=1.38d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='F' ) then
        rcov(iat)=1.35d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ne') then
        rcov(iat)=1.35d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Na') then
        rcov(iat)=3.40d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Mg') then
        rcov(iat)=2.65d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Al') then
        rcov(iat)=2.23d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Si') then
        rcov(iat)=2.09d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='P' ) then
        rcov(iat)=2.00d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='S' ) then
        rcov(iat)=1.92d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Cl') then
        rcov(iat)=1.87d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ar') then
        rcov(iat)=1.80d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='K' ) then
        rcov(iat)=4.00d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ca') then
        rcov(iat)=3.00d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Sc') then
        rcov(iat)=2.70d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ti') then
        rcov(iat)=2.70d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='V' ) then
        rcov(iat)=2.60d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Cr') then
        rcov(iat)=2.60d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Mn') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Fe') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Co') then
        rcov(iat)=2.40d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ni') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Cu') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Zn') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ga') then
        rcov(iat)=2.10d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ge') then
        rcov(iat)=2.40d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='As') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Se') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Br') then
        rcov(iat)=2.20d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Kr') then
        rcov(iat)=2.20d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Rb') then
        rcov(iat)=4.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Sr') then
        rcov(iat)=3.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Y' ) then
        rcov(iat)=3.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Zr') then
        rcov(iat)=3.00d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Nb') then
        rcov(iat)=2.92d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Mo') then
        rcov(iat)=2.83d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Tc') then
        rcov(iat)=2.75d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ru') then
        rcov(iat)=2.67d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Rh') then
        rcov(iat)=2.58d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Pd') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ag') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Cd') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='In') then
        rcov(iat)=2.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Sn') then
        rcov(iat)=2.66d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Sb') then
        rcov(iat)=2.66d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Te') then
        rcov(iat)=2.53d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='I' ) then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Xe') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Cs') then
        rcov(iat)=4.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Pa') then
        rcov(iat)=4.00d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='La') then
        rcov(iat)=3.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ce') then
        rcov(iat)=3.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Pr') then
        rcov(iat)=3.44d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Nd') then
        rcov(iat)=3.38d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Pm') then
        rcov(iat)=3.33d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Sm') then
        rcov(iat)=3.27d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Eu') then
        rcov(iat)=3.21d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Gd') then
        rcov(iat)=3.15d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Td') then
        rcov(iat)=3.09d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Dy') then
        rcov(iat)=3.03d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ho') then
        rcov(iat)=2.97d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Er') then
        rcov(iat)=2.92d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Tm') then
        rcov(iat)=2.92d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Yb') then
        rcov(iat)=2.80d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Lu') then
        rcov(iat)=2.80d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Hf') then
        rcov(iat)=2.90d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ta') then
        rcov(iat)=2.70d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='W' ) then
        rcov(iat)=2.60d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Re') then
        rcov(iat)=2.60d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Os') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Ir') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Pt') then
        rcov(iat)=2.60d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Au') then
        rcov(iat)=2.70d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Hg') then
        rcov(iat)=2.80d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Tl') then
        rcov(iat)=2.50d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Pb') then
        rcov(iat)=3.30d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Bi') then
        rcov(iat)=2.90d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Po') then
        rcov(iat)=2.80d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='At') then
        rcov(iat)=2.60d0
     else if (trim(atoms%atomnames(atoms%iatype(iat)))=='Rn') then
        rcov(iat)=2.60d0
     else
        write(*,*) 'no covalent radius stored for this atomtype ',  & 
             trim(atoms%atomnames(atoms%iatype(iat)))
     endif
     if (iproc.eq.0) write(*,*) 'RCOV:',trim(atoms%atomnames(atoms%iatype(iat))),rcov(iat)
  enddo
end subroutine give_rcov
