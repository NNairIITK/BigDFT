program BigDFT

  use module_types
!  use libBigDFT


  !implicit real(kind=8) (a-h,o-z)
  !as a general policy, I will put "implicit none" by assuming the same
  !name convention as "implicit real(kind=8) (a-h,o-z)"
  !such that the implicit statement can be commented at will

  implicit none
  include 'mpif.h'
  ! For parallel MPI execution set parallel=.true., for serial parallel=.false.
  ! this statement wil be changed by using the MPIfake.f90 file
  include 'parameters.h'
  logical :: output_wf,output_grid,calc_tail
  character(len=20) :: units
  character(len=80) :: line
  integer :: iproc,nproc,nat,ntypes,n1,n2,n3,iat,ityp,j,i_stat,i_all,ierr,infocode
  integer :: ncount_cluster
  integer :: norb,norbp
  real(kind=8) :: energy,etot,energyold,beta,sumx,sumy,sumz,tt
  !input variables
  type(input_variables) :: inputs
  type(wavefunctions_descriptors) :: wfd

  logical, dimension(:), allocatable :: lfrztyp
  character(len=6), dimension(:), allocatable :: frzsymb
  character(len=20), dimension(:), allocatable :: atomnames
  ! atomic types
  integer, dimension(:), allocatable :: iatype
  ! atomic coordinates, forces
  real(kind=8), dimension(:,:), allocatable :: rxyz,fxyz,rxyz_old
 
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:,:), pointer :: psi

  !$      interface
  !$        integer ( kind=4 ) function omp_get_num_threads ( )
  !$        end function omp_get_num_threads
  !$      end interface
  !$      interface
  !$        integer ( kind=4 ) function omp_get_thread_num ( )
  !$        end function omp_get_thread_num
  !$      end interface

  ! Start MPI in parallel version
  if (parallel) then
     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  else
     nproc=1
     iproc=0
  endif

  !initialize memory counting
  call memocc(0,iproc,'count','start')

  !$omp parallel private(iam)  shared (npr)
  !$       iam=omp_get_thread_num()
  !$       if (iam.eq.0) npr=omp_get_num_threads()
  !$       write(*,*) 'iproc,iam,npr',iproc,iam,npr
  !$omp end parallel

  !welcome screen
  if (iproc==0) call print_logo()

  !read number of atoms
  open(unit=99,file='posinp',status='old')
  read(99,*) nat,units
  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atoms     = ',nat

  allocate(rxyz_old(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz_old))*kind(rxyz_old),'rxyz_old','BigDFT')
  allocate(rxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz))*kind(rxyz),'rxyz','BigDFT')
  allocate(iatype(nat),stat=i_stat)
  call memocc(i_stat,product(shape(iatype))*kind(iatype),'iatype','BigDFT')
  allocate(fxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(fxyz))*kind(fxyz),'fxyz','BigDFT')
  allocate(atomnames(100),stat=i_stat) 
  call memocc(i_stat,product(shape(atomnames))*kind(atomnames),'atomnames','BigDFT')

  ! read atomic positions
  call read_atomic_positions(iproc,99,units,nat,ntypes,iatype,atomnames,rxyz)

  close(99)

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atom types= ',ntypes

  do ityp=1,ntypes
     if (iproc.eq.0) &
          write(*,'(1x,a,i0,a,a)') 'Atoms of type ',ityp,' are ',trim(atomnames(ityp))
  enddo

  !new way of reading the input variables, use structures
  call read_input_variables(iproc,inputs)
 
  !this array is useful for frozen atoms
  !we should modify the way in which they must be entrered
  allocate(lfrztyp(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(lfrztyp))*kind(lfrztyp),'lfrztyp','BigDFT')
  lfrztyp(:)=.false.

!!$  !  Read the first line of "input.dat"
!!$  open(unit=1,file='input.dat',status='old')
!!$  !Only the first line for the main routine (the program)
!!$  ! ngeostep Number of steps of geometry optimisation (default 500)
!!$  ! ampl     Amplitude for random displacement away from input file geometry 
!!$  ! (usually equilibrium geom.)
!!$  ! betax    Geometry optimisation
!!$  read(1,'(a80)')line
!!$  close(unit=1)
!!$  
!!$  allocate(lfrztyp(ntypes),stat=i_stat)
!!$  call memocc(i_stat,product(shape(lfrztyp))*kind(lfrztyp),'lfrztyp','BigDFT')
!!$
!!$  
!!$  read(line,*,iostat=ierror)ngeostep,ampl,betax,nfrztyp
!!$  if (ierror /= 0) then
!!$     read(line,*)ngeostep,ampl,betax
!!$     lfrztyp(:)=.false.
!!$  else
!!$     allocate(frzsymb(nfrztyp),stat=i_stat)
!!$     call memocc(i_stat,product(shape(frzsymb))*kind(frzsymb),'frzsymb','BigDFT')
!!$     read(line,*,iostat=ierror)ngeostep,ampl,betax,nfrztyp,(frzsymb(i),i=1,nfrztyp)
!!$     lfrztyp(:)=.false.
!!$     do ityp=1,nfrztyp
!!$        seek_frozen: do jtyp=1,ntypes
!!$           if (trim(frzsymb(ityp))==trim(atomnames(jtyp))) then
!!$              lfrztyp(jtyp)=.true.
!!$              exit seek_frozen
!!$           end if
!!$        end do seek_frozen
!!$     end do
!!$     i_all=-product(shape(frzsymb))*kind(frzsymb)
!!$     deallocate(frzsymb,stat=i_stat)
!!$     call memocc(i_stat,i_all,'frzsymb','BigDFT')
!!$  end if

  do iat=1,nat
     if (.not. lfrztyp(iatype(iat))) then
        call random_number(tt)
        rxyz(1,iat)=rxyz(1,iat)+inputs%randdis*tt
        call random_number(tt)
        rxyz(2,iat)=rxyz(2,iat)+inputs%randdis*tt
        call random_number(tt)
        rxyz(3,iat)=rxyz(3,iat)+inputs%randdis*tt
     end if
  enddo

  call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz,&
       psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,inputs,infocode)

  if (inputs%ncount_cluster_x > 1) then
     if (iproc ==0 ) write(*,"(a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
     ! geometry optimization
     !    betax=2.d0   ! Cincodinine
     !    betax=4.d0  ! Si H_4
     !   betax=7.5d0  ! silicon systems
     !    betax=10.d0  !  Na_Cl clusters
     ncount_cluster=1
     beta=inputs%betax
     energyold=1.d100
     call conjgrad(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,rxyz,etot,fxyz,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,inputs)
  end if

  if (iproc.eq.0) then
     sumx=0.d0
     sumy=0.d0
     sumz=0.d0
     write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom'
     do iat=1,nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
             iat,trim(atomnames(iatype(iat))),(fxyz(j,iat),j=1,3)
        sumx=sumx+fxyz(1,iat)
        sumy=sumy+fxyz(2,iat)
        sumz=sumz+fxyz(3,iat)
     enddo
     write(*,'(1x,a)')'the sum of the forces is'
     write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
     write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
     write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
  endif

  !deallocations
  !the lfrztyp must be restored in the same way as before
  i_all=-product(shape(lfrztyp))*kind(lfrztyp)
  deallocate(lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp','BigDFT')

  call deallocate_wfd(wfd,'BigDFT')

  i_all=-product(shape(atomnames))*kind(atomnames)
  deallocate(atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames','BigDFT')
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi','BigDFT')
  i_all=-product(shape(eval))*kind(eval)
  deallocate(eval,stat=i_stat)
  call memocc(i_stat,i_all,'eval','BigDFT')
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz','BigDFT')
  i_all=-product(shape(rxyz_old))*kind(rxyz_old)
  deallocate(rxyz_old,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz_old','BigDFT')
  i_all=-product(shape(iatype))*kind(iatype)
  deallocate(iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype','BigDFT')
  i_all=-product(shape(fxyz))*kind(fxyz)
  deallocate(fxyz,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz','BigDFT')


  !finalize memory counting
  call memocc(0,0,'count','stop')

  if (parallel) call MPI_FINALIZE(ierr)

 contains

   subroutine conjgrad(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,etot,gg, &
         psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,in)
     use module_types
     implicit real(kind=8) (a-h,o-z)
     implicit integer (i-n) !this line is added in view of the changement to implicit none
     type(wavefunctions_descriptors) :: wfd
     logical :: parallel
     integer :: iatype(nat)
     logical :: lfrztyp(ntypes)
     type(input_variables) :: in
     character(len=20) :: atomnames(100)
     real(kind=8), pointer :: psi(:,:), eval(:)
     integer, pointer :: keyv(:), keyg(:,:)
     dimension wpos(3,nat),gg(3,nat),rxyz_old(3,nat)
     real(kind=8), allocatable, dimension(:,:) :: tpos,gp,hh

     allocate(tpos(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(tpos))*kind(tpos),'tpos','conjgrad')
     allocate(gp(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(gp))*kind(gp),'gp','conjgrad')
     allocate(hh(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(hh))*kind(hh),'hh','conjgrad')



     anoise=1.d-4
     fluct=-1.d100
     flucto=-1.d100

!    Open a log file for conjgrad
     open(unit=16,file='conjgrad.prc',status='unknown')

     if (in%betax <= 0.d0) then
        call detbetax(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,&
             psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in)
     endif

     avbeta=0.d0
     avnum=0.d0
     nfail=0
     !        call steepdes(nat,fnrmtol,betax,alat,wpos,gg,etot,count_sd)
     call steepdes(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,etot,gg,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,fluct,flucto,fluctoo,fnrm,in)

     if (fnrm.lt.sqrt(1.d0*nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3.d0) then
        if (iproc.eq.0) write(16,*) 'Converged before entering CG',iproc
        return
     endif

12345   continue

     etotprec=etot
     do iat=1,nat
        hh(1,iat)=gg(1,iat)
        hh(2,iat)=gg(2,iat)
        hh(3,iat)=gg(3,iat)
     end do

     beta0=4.d0*in%betax
     it=0
1000 it=it+1

     !C line minimize along hh ----

     do iat=1,nat
        if (lfrztyp(iatype(iat))) then
           tpos(1,iat)=wpos(1,iat)
           tpos(2,iat)=wpos(2,iat)
           tpos(3,iat)=wpos(3,iat)
        else
           tpos(1,iat)=wpos(1,iat)+beta0*hh(1,iat)
           tpos(2,iat)=wpos(2,iat)+beta0*hh(2,iat)
           tpos(3,iat)=wpos(3,iat)+beta0*hh(3,iat)
        end if
     end do
!        call energyandforces(nat,alat,tpos,gp,tetot,count_cg)
     in%inputPsiId=1
     in%output_grid=.false.
     in%output_wf=.false.
     call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,tpos,tetot,gp,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)

     do iat=1,nat
        rxyz_old(1,iat)=tpos(1,iat) 
        rxyz_old(2,iat)=tpos(2,iat) 
        rxyz_old(3,iat)=tpos(3,iat)
     enddo
     ncount_cluster=ncount_cluster+1

     !C projection of gradients at beta=0 and beta onto hh
     y0=0.d0
     y1=0.d0
     do iat=1,nat
        y0=y0+gg(1,iat)*hh(1,iat)+gg(2,iat)*hh(2,iat)+gg(3,iat)*hh(3,iat)
        y1=y1+gp(1,iat)*hh(1,iat)+gp(2,iat)*hh(2,iat)+gp(3,iat)*hh(3,iat)
     end do
     tt=y0/(y0-y1)
     if (iproc.eq.0) then
        write(16,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
        write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
     end if

     beta=beta0*max(min(tt,2.d0),-.25d0)
     !        beta=beta0*max(min(tt,1.d0),-.25d0)
     !        beta=beta0*max(min(tt,10.d0),-1.d0)
     !        beta=beta0*tt
     do iat=1,nat
        tpos(1,iat)=wpos(1,iat)
        tpos(2,iat)=wpos(2,iat)
        tpos(3,iat)=wpos(3,iat)
        if ( .not. lfrztyp(iatype(iat))) then
           wpos(1,iat)=wpos(1,iat)+beta*hh(1,iat)
           wpos(2,iat)=wpos(2,iat)+beta*hh(2,iat)
           wpos(3,iat)=wpos(3,iat)+beta*hh(3,iat)
        end if
     end do
     avbeta=avbeta+beta/in%betax
     avnum=avnum+1.d0

     !C new gradient
     do iat=1,nat
        gp(1,iat)=gg(1,iat)
        gp(2,iat)=gg(2,iat)
        gp(3,iat)=gg(3,iat)
     end do

!        call energyandforces(nat,alat,wpos,gg,etot,count_cg)
     call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,wpos,etot,gg,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)

     do iat=1,nat
        rxyz_old(1,iat) = wpos(1,iat) 
        rxyz_old(2,iat) = wpos(2,iat) 
        rxyz_old(3,iat) = wpos(3,iat)
     enddo

     ncount_cluster=ncount_cluster+1
     sumx=0.d0 
     sumy=0.d0 
     sumz=0.d0

     do iat=1,nat
        sumx=sumx+gg(1,iat) 
        sumy=sumy+gg(2,iat) 
        sumz=sumz+gg(3,iat)
     end do

     if (etot.gt.etotprec+anoise) then

        if (iproc.eq.0) write(16,*) 'switching back to SD:etot,etotprec',it,etot,etotprec
        if (iproc.eq.0) write(*,*) 'switching back to SD:etot,etotprec',it,etot,etotprec
        do iat=1,nat
           wpos(1,iat)=tpos(1,iat)
           wpos(2,iat)=tpos(2,iat)
           wpos(3,iat)=tpos(3,iat)
        end do

     call steepdes(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,etot,gg,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,fluct,flucto,fluctoo,fnrm,in)

        goto 12345

     endif
     etotprec=etot
     if (iproc.eq.0) call wtposout(ncount_cluster,etot,nat,wpos,atomnames,iatype)
     fluctoo=flucto
     flucto=fluct
     fluct=sumx**2+sumy**2+sumz**2

     obenx=0.d0
     obeny=0.d0
     obenz=0.d0
     unten=0.d0
     fnrm=0.d0
     do iat=1,nat
        obenx=obenx+(gg(1,iat)-gp(1,iat))*gg(1,iat)
        obeny=obeny+(gg(2,iat)-gp(2,iat))*gg(2,iat)
        obenz=obenz+(gg(3,iat)-gp(3,iat))*gg(3,iat)
        unten=unten+gp(1,iat)**2+gp(2,iat)**2+gp(3,iat)**2
        if (.not. lfrztyp(iatype(iat))) then
           fnrm=fnrm+gg(1,iat)**2+gg(2,iat)**2+gg(3,iat)**2
        end if
     end do
     if (iproc.eq.0) then
        write(16,'(i5,1x,e12.5,1x,e21.14,a,1x,e9.2)')it,sqrt(fnrm),etot,' CG ',beta/in%betax
        write(16,*) 'fnrm2,flucts',&
             fnrm,sqrt(real(nat,kind=8))*(fluct+flucto+fluctoo)*in%frac_fluct/3.d0
     end if

     if (fnrm.lt.sqrt(1.d0*nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3.d0) goto 2000

     if (ncount_cluster.gt.in%ncount_cluster_x) then 
        if (iproc.eq.0) write(*,*) 'ncount_cluster in CG',ncount_cluster
        goto 2000
     endif

     if (it.eq.500) then
        if (iproc.eq.0) then
           write(16,*) &
                'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
           write(*,*) 'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
        end if
        do iat=1,nat
           wpos(1,iat)=tpos(1,iat)
           wpos(2,iat)=tpos(2,iat)
           wpos(3,iat)=tpos(3,iat)
        end do

        call steepdes(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,etot,gg,&
             psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,&
             fluct,flucto,fluctoo,fnrm,in)
        
        nfail=nfail+1
        if (nfail.ge.100) stop 'too many failures of CONJG'
        goto 12345
     endif
     
     rlambda=(obenx+obeny+obenz)/unten
     do iat=1,nat
        hh(1,iat)=gg(1,iat)+rlambda*hh(1,iat)
        hh(2,iat)=gg(2,iat)+rlambda*hh(2,iat)
        hh(3,iat)=gg(3,iat)+rlambda*hh(3,iat)
     end do
     goto 1000    
2000 continue

!!        write(6,*) 'CG finished',it,fnrm,etot
        if (iproc.eq.0) then
           write(16,*) 'average CG stepsize in terms of betax',avbeta/avnum,iproc
           write(*,*) 'average CG stepsize in terms of betax',avbeta/avnum,iproc
        end if


!    Close the file
     close(unit=16)

     i_all=-product(shape(tpos))*kind(tpos)
     deallocate(tpos,stat=i_stat)
     call memocc(i_stat,i_all,'tpos','conjgrad')
     i_all=-product(shape(gp))*kind(gp)
     deallocate(gp,stat=i_stat)
     call memocc(i_stat,i_all,'gp','conjgrad')
     i_all=-product(shape(hh))*kind(hh)
     deallocate(hh,stat=i_stat)
     call memocc(i_stat,i_all,'hh','conjgrad')

   end subroutine conjgrad


   subroutine steepdes(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,etot,ff,&
             psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,&
             fluct,flucto,fluctoo,fnrm,in)
     use module_types
     implicit real(kind=8) (a-h,o-z)
     implicit integer (i-n) !this line is added in view of the changement to implicit none
     logical :: parallel
     integer :: iatype(nat)
     logical :: lfrztyp(ntypes)
     type(input_variables) :: in
     type(wavefunctions_descriptors) :: wfd
     character(len=20) :: atomnames(100)
     real(kind=8), pointer :: psi(:,:), eval(:)
     integer, pointer :: keyv(:), keyg(:,:)
     dimension wpos(3,nat),ff(3,nat),rxyz_old(3,nat)
     real(kind=8), allocatable, dimension(:,:) :: tpos
     logical care

     allocate(tpos(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(tpos))*kind(tpos),'tpos','steepdes')

     anoise=0.d-4

     beta=in%betax
     care=.true.
     nsatur=0
     etotitm2=1.d100
     fnrmitm2=1.d100
     etotitm1=1.d100
     fnrmitm1=1.d100

     do iat=1,nat
        tpos(1,iat)=wpos(1,iat)
        tpos(2,iat)=wpos(2,iat)
        tpos(3,iat)=wpos(3,iat)
     end do

     itot=0
12345 continue

     if (ncount_cluster.gt.in%ncount_cluster_x) then 
        if (iproc.eq.0) write(*,*) 'ncount_cluster in SD1',ncount_cluster
        goto 2000
     endif

     nitsd=500
     itsd=0
1000 itsd=itsd+1
     itot=itot+1

!        call energyandforces(nat,alat,wpos,ff,etot,count_sd)

     in%inputPsiId=1
     in%output_grid=.false.
     in%output_wf=.false.
     call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,wpos,etot,ff,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)

     ncount_cluster=ncount_cluster+1
     do iat=1,nat
        rxyz_old(1,iat)=wpos(1,iat)
        rxyz_old(2,iat)=wpos(2,iat)
        rxyz_old(3,iat)=wpos(3,iat)
     enddo
     if (care .and. etot.gt.etotitm1+anoise) then
        do iat=1,nat
           wpos(1,iat)=tpos(1,iat)
           wpos(2,iat)=tpos(2,iat)
           wpos(3,iat)=tpos(3,iat)
        end do
        beta=.5d0*beta
        if (iproc.eq.0) write(16,'(a,1x,e9.2,1x,i5,2(1x,e21.14))') &
             'SD reset, beta,itsd,etot,etotitm1= ',beta,itsd,etot,etotitm1
        if (beta.le.1.d-1*in%betax) then
           if (iproc.eq.0) write(16,*) &
           'beta getting too small, do not care anymore if energy goes up'
           care=.false.
        endif
        goto 12345
     endif
        
     t1=0.d0 
     t2=0.d0 
     t3=0.d0
     do iat=1,nat
        if (.not. lfrztyp(iatype(iat))) then
           t1=t1+ff(1,iat)**2 
           t2=t2+ff(2,iat)**2 
           t3=t3+ff(3,iat)**2
        end if
     enddo

     fnrm=t1+t2+t3
     de1=etot-etotitm1
     de2=etot-2.d0*etotitm1+etotitm2
     df1=fnrm-fnrmitm1
     df2=fnrm-2.d0*fnrmitm1+fnrmitm2
     if (iproc.eq.0) write(16,'(5(1x,e11.4),1x,i3)') fnrm/fnrmitm1, de1,de2,df1,df2,nsatur

     if (care .and. itsd.ge.3 .and. beta.eq.in%betax .and. fnrm/fnrmitm1.gt..8d0 .and. &
          de1.gt.-.1d0 .and. fnrm.le..1d0 .and. & 
          de1.lt.anoise .and. df1.lt.anoise .and. &
          de2.gt.-2.d0*anoise .and. df2.gt.-2.d0*anoise) then 
        nsatur=nsatur+1
     else
        nsatur=0
     endif
     if (iproc.eq.0) then 
        write(16,'(i5,1x,e12.5,1x,e21.14,a)') itsd,sqrt(fnrm),etot,' SD '
        call wtposout(ncount_cluster,etot,nat,wpos,atomnames,iatype)
     end if

     fluctoo=flucto
     flucto=fluct
     sumx=0.d0 
     sumy=0.d0 
     sumz=0.d0
     do iat=1,nat
        sumx=sumx+ff(1,iat) 
        sumy=sumy+ff(2,iat) 
        sumz=sumz+ff(3,iat)
     end do
     fluct=sumx**2+sumy**2+sumz**2
     if (iproc.eq.0) write(16,*) &
          'fnrm2,flucts',fnrm,sqrt(1.d0*nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3.d0

     if (fnrm < sqrt(1.d0*nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3.d0 .or. nsatur > 5 ) &
          goto 2000

     if (ncount_cluster > in%ncount_cluster_x) then 
        if (iproc.eq.0) write(*,*) 'ncount_cluster in SD2',ncount_cluster
        goto 2000
     endif

     if (itsd >= nitsd) then 
        if (iproc.eq.0) write(16,'(a,i5,1x,e10.3,1x,e21.14)') &
        'SD: NO CONVERGENCE:itsd,fnrm,etot',itsd,fnrm,etot
        goto 2000
     endif
     if (itot >= nitsd) then
        if (iproc.eq.0) write(16,'(a,i5,i5,1x,e10.3,1x,e21.14)') &
        'SD: NO CONVERGENCE:itsd,itot,fnrm,etot:',itsd,itot,fnrm,etot
        goto 2000
     endif
     
     etotitm2=etotitm1
     etotitm1=etot
     fnrmitm2=fnrmitm1
     fnrmitm1=fnrm

     beta=min(1.2d0*beta,in%betax)
     
     if (beta == in%betax) then 
        if (iproc.eq.0) write(16,*) 'beta=betax'
        care=.true.
     endif
     if (iproc.eq.0) write(16,*) 'beta=',beta
     do iat=1,nat
        tpos(1,iat)=wpos(1,iat)
        tpos(2,iat)=wpos(2,iat)
        tpos(3,iat)=wpos(3,iat)
        if ( .not. lfrztyp(iatype(iat))) then
           wpos(1,iat)=wpos(1,iat)+beta*ff(1,iat)
           wpos(2,iat)=wpos(2,iat)+beta*ff(2,iat)
           wpos(3,iat)=wpos(3,iat)+beta*ff(3,iat)
        end if
     end do
     goto 1000
2000        continue
     if (iproc.eq.0) write(16,*) 'SD FINISHED',iproc
     
     i_all=-product(shape(tpos))*kind(tpos)
     deallocate(tpos,stat=i_stat)
     call memocc(i_stat,i_all,'tpos','steepdes')

   end subroutine steepdes

   subroutine detbetax(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,pos,&
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in)
     ! determines stepsize betax
     use module_types
     implicit real(kind=8) (a-h,o-z)
     implicit integer (i-n) !this line is added in view of the changement to implicit none
     logical :: parallel
     integer :: iatype(nat)
     logical :: lfrztyp(ntypes)
     type(input_variables) :: in
     type(wavefunctions_descriptors) :: wfd
     character(len=20) :: atomnames(100)
     real(kind=8), pointer :: psi(:,:), eval(:)
     integer, pointer :: keyv(:), keyg(:,:)
     dimension pos(3,nat),alat(3),rxyz_old(3,nat)
     real(kind=8), allocatable, dimension(:,:) :: tpos,ff,gg

     allocate(tpos(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(tpos))*kind(tpos),'tpos','detbetax')
     allocate(ff(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(ff))*kind(ff),'ff','detbetax')
     allocate(gg(3,nat),stat=i_stat)
     call memocc(i_stat,product(shape(gg))*kind(gg),'gg','detbetax')
    
     beta0=abs(in%betax)
     beta=1.d100

     nsuc=0
100  continue
     !        call energyandforces(nat,alat,pos,ff,etotm1,count)

     in%inputPsiId=1
     in%output_grid=.false.
     in%output_wf=.false.
     call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,pos,etotm1,ff,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
     ncount_cluster=ncount_cluster+1
     do iat=1,nat
        rxyz_old(1,iat) = pos(1,iat)
        rxyz_old(2,iat) = pos(2,iat) 
        rxyz_old(3,iat) = pos(3,iat)
     enddo
     sum=0.d0
     
     do iat=1,nat
        tpos(1,iat)=pos(1,iat)
        tpos(2,iat)=pos(2,iat)
        tpos(3,iat)=pos(3,iat)
        t1=beta0*ff(1,iat)
        t2=beta0*ff(2,iat)
        t3=beta0*ff(3,iat)
        sum=sum+t1**2+t2**2+t3**2
        if (.not. lfrztyp(iatype(iat))) then
           pos(1,iat)=pos(1,iat)+t1
           pos(2,iat)=pos(2,iat)+t2
           pos(3,iat)=pos(3,iat)+t3
        end if
     enddo
        !        call energyandforces(nat,alat,pos,gg,etot0,count)

     call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,pos,etot0,gg,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)

     ncount_cluster=ncount_cluster+1
     do iat=1,nat
        rxyz_old(1,iat)=pos(1,iat) 
        rxyz_old(2,iat)=pos(2,iat) 
        rxyz_old(3,iat)=pos(3,iat)
     enddo

     if (etot0.gt.etotm1) then
        do iat=1,nat
           pos(1,iat)=tpos(1,iat)
           pos(2,iat)=tpos(2,iat)
           pos(3,iat)=tpos(3,iat)
        enddo
        beta0=.5d0*beta0
        if (iproc.eq.0) write(16,*) 'beta0 reset',beta0
        goto 100
     endif
     do iat=1,nat
        if (lfrztyp(iatype(iat))) then
           tpos(1,iat)=pos(1,iat)
           tpos(2,iat)=pos(2,iat)
           tpos(3,iat)=pos(3,iat)
        else
           tpos(1,iat)=pos(1,iat)+beta0*ff(1,iat)
           tpos(2,iat)=pos(2,iat)+beta0*ff(2,iat)
           tpos(3,iat)=pos(3,iat)+beta0*ff(3,iat)
        end if
     enddo
        !        call energyandforces(nat,alat,tpos,gg,etotp1,count)
     call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,tpos,etotp1,gg,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)

     ncount_cluster=ncount_cluster+1
     do iat=1,nat
        rxyz_old(1,iat)=tpos(1,iat) 
        rxyz_old(2,iat)=tpos(2,iat) 
        rxyz_old(3,iat)=tpos(3,iat)
     enddo

     if (iproc.eq.0) write(16,'(a,3(1x,e21.14))') 'etotm1,etot0,etotp1',etotm1,etot0,etotp1
     der2=(etotp1+etotm1-2.d0*etot0)/sum
     tt=.25d0/der2
     beta0=.125d0/der2
     if (iproc.eq.0) write(16,*) 'der2,tt=',der2,tt
     if (der2.gt.0.d0) then
        nsuc=nsuc+1
        beta=min(beta,.5d0/der2)
     endif
     do iat=1,nat
        if ( .not. lfrztyp(iatype(iat))) then
           pos(1,iat)=pos(1,iat)+tt*ff(1,iat)
           pos(2,iat)=pos(2,iat)+tt*ff(2,iat)
           pos(3,iat)=pos(3,iat)+tt*ff(3,iat)
        end if
     enddo

!!$     if (count > 100.d0) then
!!$        if (iproc.eq.0) write(16,*) 'CANNOT DETERMINE betax '
!!$        stop
!!$     endif
     if (nsuc < 3) goto 100
        
     in%betax=beta
        
     if (iproc.eq.0) write(16,*) 'betax=',in%betax
        
     i_all=-product(shape(tpos))*kind(tpos)
     deallocate(tpos,stat=i_stat)
     call memocc(i_stat,i_all,'tpos','detbetax')
     i_all=-product(shape(ff))*kind(ff)
     deallocate(ff,stat=i_stat)
     call memocc(i_stat,i_all,'ff','detbetax')
     i_all=-product(shape(gg))*kind(gg)
     deallocate(gg,stat=i_stat)
     call memocc(i_stat,i_all,'gg','detbetax')

   end subroutine detbetax

   
   subroutine call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz,&
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
     use module_types
     implicit none
     type(input_variables) :: in
     type(wavefunctions_descriptors) :: wfd
     logical, intent(in) :: parallel
     integer, intent(in) :: iproc,nproc,nat,ntypes,norbp,norb
     integer, intent(inout) :: infocode,n1,n2,n3
     integer :: i_stat,i_all,ierr,inputPsiId_orig
     real(kind=8), intent(out) :: energy
     character(len=20), dimension(100), intent(in) :: atomnames
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(3,nat), intent(in) :: rxyz_old
     real(kind=8), dimension(3,nat), intent(inout) :: rxyz
     real(kind=8), dimension(3,nat), intent(out) :: fxyz
     real(kind=8), dimension(:), pointer :: eval
     real(kind=8), dimension(:,:), pointer :: psi

     inputPsiId_orig=in%inputPsiId

     loop_cluster: do

        if (inputPsiId_orig == 0 .and. associated(psi)) then
           i_all=-product(shape(psi))*kind(psi)
           deallocate(psi,stat=i_stat)
           call memocc(i_stat,i_all,'psi','call_cluster')
           i_all=-product(shape(eval))*kind(eval)
           deallocate(eval,stat=i_stat)
           call memocc(i_stat,i_all,'eval','call_cluster')

           call deallocate_wfd(wfd,'call_cluster')
        end if

        call cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz,&
             psi,wfd,norbp,norb,eval,&
             n1,n2,n3,rxyz_old,in,infocode)

        if (in%inputPsiId==1 .and. infocode==2) then
           in%inputPsiId=0
        else if (in%inputPsiId == 0 .and. infocode==3) then
           if (iproc.eq.0) then
              write(*,'(1x,a)')'Convergence error, cannot proceed.'
              write(*,'(1x,a)')' writing positions in file posout_999.ascii then exiting'
              call wtposout(999,energy,nat,rxyz,atomnames,iatype)
           end if

           i_all=-product(shape(psi))*kind(psi)
           deallocate(psi,stat=i_stat)
           call memocc(i_stat,i_all,'psi','call_cluster')
           i_all=-product(shape(eval))*kind(eval)
           deallocate(eval,stat=i_stat)
           call memocc(i_stat,i_all,'eval','call_cluster')

           call deallocate_wfd(wfd,'call_cluster')
           !finalize memory counting (there are still the positions and the forces allocated)
           call memocc(0,0,'count','stop')

           if (parallel) call MPI_FINALIZE(ierr)

           stop
        else
           exit loop_cluster
        end if

     end do loop_cluster

     !preserve the previous value
     inputPsiId_orig=in%inputPsiId

   end subroutine call_cluster

 end program BigDFT


 subroutine wtposout(igeostep,energy,nat,rxyz,atomnames,iatype)

   implicit none
   integer, intent(in) :: igeostep,nat
   real(kind=8), intent(in) :: energy
   character(len=20), dimension(100), intent(in) :: atomnames(100)
   integer, dimension(nat), intent(in) :: iatype
   real(kind=8), dimension(3,nat), intent(in) :: rxyz
   !local variables
   character(len=3) :: fn
   character(len=20) :: filename
   integer :: iat,j
   real(kind=8) :: xmax,ymax,zmax

   write(fn,'(i3.3)') igeostep
   filename = 'posout_'//fn//'.ascii'
   open(unit=9,file=filename)
   xmax=0.d0 ; ymax=0.d0 ; zmax=0.d0
   do iat=1,nat
      xmax=max(rxyz(1,iat),xmax)
      ymax=max(rxyz(2,iat),ymax)
      zmax=max(rxyz(3,iat),zmax)
   enddo
   write(9,*) nat,' atomic ', energy,igeostep
   write(9,*) xmax+5.d0, 0.d0, ymax+5.d0
   write(9,*) 0.d0, 0.d0, zmax+5.d0
   do iat=1,nat
      write(9,'(3(1x,e21.14),2x,a10)') (rxyz(j,iat),j=1,3),atomnames(iatype(iat))
   enddo
   close(unit=9)

 end subroutine wtposout
