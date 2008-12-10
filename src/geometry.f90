subroutine timeleft(tt)
! MODIFIED version for refined time limit on restart of global.f90.
! Only difference: Calls routine CPUtime(tt)
    implicit real*8 (a-h,o-z)
          open(unit=55,file='CPUlimit',status='unknown')
          read(55,*,iostat=ierr) timelimit ! in hours
          if(ierr/=0)timelimit=1d6
          close(55)
          call cpu_time(tcpu)
          tt=timelimit-tcpu/3600d0 ! in hours
end subroutine timeleft

subroutine conjgrad(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_cluster,in)
  use module_base
  use module_types
  use module_interfaces, except_this_one => conjgrad
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_cluster
  real(gp), intent(out) :: etot
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: fxyz
  !local variables
  character(len=*), parameter :: subname='conjgrad'  
  integer :: nfail,it,iat,i_all,i_stat,infocode
  real(gp) :: anoise,fluct,flucto,fluctoo,avbeta,avnum,fnrm,etotprec,beta0,beta
  real(gp) :: y0,y1,tt,sumx,sumy,sumz,obenx,obeny,obenz,unten,rlambda,tetot,forcemax
  real(gp), dimension(:,:), allocatable :: tpos,gpf,hh

  allocate(tpos(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,tpos,'tpos',subname)
  allocate(gpf(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gpf,'gpf',subname)
  allocate(hh(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,hh,'hh',subname)

  anoise=1.e-4_gp
  fluct=0._gp
  flucto=0._gp
  fluctoo=0._gp

  !write the first position
  if (iproc.eq.0) call  wtposout(ncount_cluster,etot,rxyz,at)
  !    Open a log file for conjgrad
  open(unit=16,file='conjgrad.prc',access='append')

  if (in%betax <= 0._gp) then
     call detbetax(nproc,iproc,at,rxyz,rst,in)
  endif

  avbeta=0._gp
  avnum=0._gp
  nfail=0
 
  !start with a steepest descent algorithm
  call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_cluster,fluct,flucto,fluctoo,fnrm,in)

  !calculate the max of the forces
  forcemax=maxval(abs(fxyz))

  !control whether the convergence criterion is reached after SD
  if (fnrm < sqrt(1._gp*at%nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp .or. &
       forcemax < in%forcemax ) then
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

        do iat=1,at%nat
           if (at%lfrztyp(iat)) then
              tpos(1,iat)=rxyz(1,iat)
              tpos(2,iat)=rxyz(2,iat)
              tpos(3,iat)=rxyz(3,iat)
           else
              if (at%geocode == 'P') then
                 tpos(1,iat)=modulo(rxyz(1,iat)+beta0*hh(1,iat),at%alat1)
                 tpos(2,iat)=modulo(rxyz(2,iat)+beta0*hh(2,iat),at%alat2)
                 tpos(3,iat)=modulo(rxyz(3,iat)+beta0*hh(3,iat),at%alat3)
              else if (at%geocode == 'S') then
                 tpos(1,iat)=modulo(rxyz(1,iat)+beta0*hh(1,iat),at%alat1)
                 tpos(2,iat)=rxyz(2,iat)+beta0*hh(2,iat)
                 tpos(3,iat)=modulo(rxyz(3,iat)+beta0*hh(3,iat),at%alat3)
              else
                 tpos(1,iat)=rxyz(1,iat)+beta0*hh(1,iat)
                 tpos(2,iat)=rxyz(2,iat)+beta0*hh(2,iat)
                 tpos(3,iat)=rxyz(3,iat)+beta0*hh(3,iat)
              end if
           end if
        end do

        in%inputPsiId=1
        in%output_grid=0
        in%output_wf=.false.
        call call_bigdft(nproc,iproc,at,tpos,in,tetot,gpf,rst,infocode)

!!$        !useless, already done in call_bigdft routine
!!$        do iat=1,at%nat
!!$           rxyz_old(1,iat)=tpos(1,iat) 
!!$           rxyz_old(2,iat)=tpos(2,iat) 
!!$           rxyz_old(3,iat)=tpos(3,iat)
!!$        enddo

        ncount_cluster=ncount_cluster+1

        !C projection of gradients at beta=0 and beta onto hh
        y0=0._gp
        y1=0._gp
        do iat=1,at%nat
           y0=y0+fxyz(1,iat)*hh(1,iat)+fxyz(2,iat)*hh(2,iat)+fxyz(3,iat)*hh(3,iat)
           y1=y1+gpf(1,iat)*hh(1,iat)+gpf(2,iat)*hh(2,iat)+gpf(3,iat)*hh(3,iat)
        end do
        tt=y0/(y0-y1)
        if (iproc.eq.0) then
           write(16,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
           write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
        end if

        beta=beta0*max(min(tt,2._gp),-.25_gp)
        !        beta=beta0*max(min(tt,1._gp),-.25_gp)
        !        beta=beta0*max(min(tt,10._gp),-1._gp)
        !        beta=beta0*tt
        do iat=1,at%nat
           tpos(1,iat)=rxyz(1,iat)
           tpos(2,iat)=rxyz(2,iat)
           tpos(3,iat)=rxyz(3,iat)
           if ( .not. at%lfrztyp(iat)) then
              if (at%geocode == 'P') then
                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*hh(1,iat),at%alat1)
                 rxyz(2,iat)=modulo(rxyz(2,iat)+beta*hh(2,iat),at%alat2)
                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*hh(3,iat),at%alat3)
              else if (at%geocode == 'S') then
                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*hh(1,iat),at%alat1)
                 rxyz(2,iat)=rxyz(2,iat)+beta*hh(2,iat)
                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*hh(3,iat),at%alat3)
              else
                 rxyz(1,iat)=rxyz(1,iat)+beta*hh(1,iat)
                 rxyz(2,iat)=rxyz(2,iat)+beta*hh(2,iat)
                 rxyz(3,iat)=rxyz(3,iat)+beta*hh(3,iat)
              end if

           end if
        end do
        avbeta=avbeta+beta/in%betax
        avnum=avnum+1._gp

        !C new gradient
        do iat=1,at%nat
           gpf(1,iat)=fxyz(1,iat)
           gpf(2,iat)=fxyz(2,iat)
           gpf(3,iat)=fxyz(3,iat)
        end do

        call call_bigdft(nproc,iproc,at,rxyz,in,etot,fxyz,rst,infocode)

!!$        !useless, already done in call_bigdft routine
!!$        do iat=1,at%nat
!!$           rxyz_old(1,iat) = rxyz(1,iat) 
!!$           rxyz_old(2,iat) = rxyz(2,iat) 
!!$           rxyz_old(3,iat) = rxyz(3,iat)
!!$        enddo

        ncount_cluster=ncount_cluster+1
        sumx=0._gp 
        sumy=0._gp 
        sumz=0._gp

        do iat=1,at%nat
           sumx=sumx+fxyz(1,iat) 
           sumy=sumy+fxyz(2,iat) 
           sumz=sumz+fxyz(3,iat)
        end do

        !if the energy goes up (a small tolerance of anoise is allowed)
        !switch back to SD
        if (etot > etotprec+anoise) then

           if (iproc.eq.0) write(16,*) 'switching back to SD:etot,etotprec',it,etot,etotprec
           if (iproc.eq.0) write(*,*) ' switching back to SD:etot,etotprec',it,etot,etotprec
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_cluster,&
                fluct,flucto,fluctoo,fnrm,in)

           !calculate the max of the forces
           forcemax=maxval(abs(fxyz))

           if (fnrm < sqrt(1._gp*at%nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp .or. &
               forcemax < in%forcemax ) then
              if (iproc.eq.0) write(16,*) 'Converged in switch back SD',iproc
              call close_and_deallocate
              return
           endif
           cycle redo_cg
        endif

        etotprec=etot
        if (iproc.eq.0) call wtposout(ncount_cluster,etot,rxyz,at)
        !if (iproc.eq.0) write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'CG ',ncount_cluster,etot,sqrt(fnrm)
        fluctoo=flucto
        flucto=fluct
        fluct=sumx**2+sumy**2+sumz**2

        obenx=0._gp
        obeny=0._gp
        obenz=0._gp
        unten=0._gp
        fnrm=0._gp
        forcemax=0._gp
        do iat=1,at%nat
           obenx=obenx+(fxyz(1,iat)-gpf(1,iat))*fxyz(1,iat)
           obeny=obeny+(fxyz(2,iat)-gpf(2,iat))*fxyz(2,iat)
           obenz=obenz+(fxyz(3,iat)-gpf(3,iat))*fxyz(3,iat)
           unten=unten+gpf(1,iat)**2+gpf(2,iat)**2+gpf(3,iat)**2
           if (.not. at%lfrztyp(iat)) then
              fnrm=fnrm+fxyz(1,iat)**2+fxyz(2,iat)**2+fxyz(3,iat)**2
              forcemax=max(forcemax,abs(fxyz(1,iat)),abs(fxyz(2,iat)),abs(fxyz(3,iat)))
           end if
        end do
        if (iproc.eq.0) then
           write(16,'(i5,1x,e12.5,1x,e21.14,a,1x,e9.2)')it,sqrt(fnrm),etot,' CG ',beta/in%betax
           write(16,*) 'fnrm2,flucts',&
                fnrm,sqrt(real(at%nat,kind=8))*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp
           write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                'FORCES norm(Ha/Bohr): maxval=',    forcemax,'fnrm=',    fnrm    ,&
                'fluct=',sqrt(real(at%nat,kind=8))*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp
        end if

        if (fnrm < sqrt(1._gp*at%nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp .or. &
             forcemax < in%forcemax) exit loop_cg

        if (ncount_cluster.gt.in%ncount_cluster_x) then 
           if (iproc.eq.0) write(*,*) 'ncount_cluster in CG',ncount_cluster
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
              write(16,*) &
                   'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
              write(*,*) 'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
           end if
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_cluster,&
                fluct,flucto,fluctoo,fnrm,in)

           !calculate the max of the forces
           forcemax=maxval(abs(fxyz))

           if (fnrm < sqrt(1._gp*at%nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp .or. &
                forcemax < in%forcemax) then
              if (iproc.eq.0) write(16,*) 'Converged in back up SD',iproc
              call close_and_deallocate
              return
           endif

           nfail=nfail+1
           if (nfail.ge.100) stop 'too many failures of CONJG'
           cycle redo_cg
        endif

        rlambda=(obenx+obeny+obenz)/unten
        do iat=1,at%nat
           hh(1,iat)=fxyz(1,iat)+rlambda*hh(1,iat)
           hh(2,iat)=fxyz(2,iat)+rlambda*hh(2,iat)
           hh(3,iat)=fxyz(3,iat)+rlambda*hh(3,iat)
        end do
     end do loop_cg
     exit redo_cg
  end do redo_cg

  !!        write(6,*) 'CG finished',it,fnrm,etot
  if (iproc.eq.0) then
     write(16,*) 'average CG stepsize in terms of betax',avbeta/avnum,iproc
     write(*,*) 'average CG stepsize in terms of betax',avbeta/avnum,iproc
  end if

  call close_and_deallocate

contains

  subroutine close_and_deallocate
    use module_base
    !    Close the file
    close(unit=16)
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
  
  subroutine steepdes(nproc,iproc,at,rxyz,etot,ff,rst,ncount_cluster,&
       fluct,flucto,fluctoo,fnrm,in)
    use module_base
    use module_types
    use module_interfaces
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: at
    type(input_variables), intent(inout) :: in
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_cluster
    real(gp), dimension(3,at%nat), intent(inout) :: rxyz
    real(gp), intent(out) :: fluct,flucto,fluctoo,fnrm,etot
    real(gp), dimension(3,at%nat), intent(out) ::ff
    !local variables
    character(len=*), parameter :: subname='steepdes'
    logical :: care
    integer :: nsatur,iat,itot,nitsd,itsd
    real(gp) :: etotitm2,fnrmitm2,etotitm1,fnrmitm1,anoise,sumx,sumy,sumz
    real(gp) :: forcemax,t1,t2,t3,de1,de2,df1,df2
    real(gp), allocatable, dimension(:,:) :: tpos

    allocate(tpos(3,at%nat+ndebug),stat=i_stat)
    call memocc(i_stat,tpos,'tpos',subname)

    anoise=0.d-4

    beta=in%betax
    care=.true.
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
       if (ncount_cluster.gt.in%ncount_cluster_x) then 
          if (iproc.eq.0) then
             write(*,*) 'ncount_cluster in SD1',ncount_cluster
             write(16,*) 'SD FINISHED',iproc
          end if

          i_all=-product(shape(tpos))*kind(tpos)
          deallocate(tpos,stat=i_stat)
          call memocc(i_stat,i_all,'tpos',subname)
          return
       endif

       nitsd=500
       itsd=0
       loop_sd: do
          itsd=itsd+1
          itot=itot+1

          in%inputPsiId=1
          in%output_grid=0
          in%output_wf=.false.
          call call_bigdft(nproc,iproc,at,rxyz,in,etot,ff,rst,infocode)

          ncount_cluster=ncount_cluster+1

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

          t1=0._gp 
          t2=0._gp 
          t3=0._gp
          forcemax=0._gp
          do iat=1,at%nat
             if (.not. at%lfrztyp(iat)) then
                t1=t1+ff(1,iat)**2 
                t2=t2+ff(2,iat)**2 
                t3=t3+ff(3,iat)**2
                forcemax=max(forcemax,abs(ff(1,iat)),abs(ff(2,iat)),abs(ff(3,iat)))
             end if
          enddo

          !this is the norm of the forces of non-blocked atoms
          fnrm=t1+t2+t3

          !first and second derivatives of the energy and of the norm of the forces
          !(in units of beta steps)
          de1=etot-etotitm1
          de2=etot-2._gp*etotitm1+etotitm2
          df1=fnrm-fnrmitm1
          df2=fnrm-2._gp*fnrmitm1+fnrmitm2

          if (iproc.eq.0) then
             write(16,'(5(1x,e11.4),1x,i3)') fnrm/fnrmitm1, de1,de2,df1,df2,nsatur
             write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                  'FORCES norm(Ha/Bohr): maxval=',    forcemax,'fnrm=',    fnrm    ,'fluct=', &
                  sqrt(real(at%nat,kind=8))*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp
          end if

          !control whether we are in a situation in which SD do not change things too much
          if (care .and. itsd >= 3 .and. beta == in%betax .and. &
               df1 < anoise .and. &                              !forces are decreasing
               de1 > -.1_gp .and. de1 < anoise .and. &            !energy slowly decreasing
               fnrm <= .1_gp .and. fnrm/fnrmitm1 > .8_gp .and. &   !norm of forces is saturating
               de2 > -2._gp*anoise .and. df2 > -2._gp*anoise) then !close to a local minimum (E&F)
             nsatur=nsatur+1
          else
             nsatur=0
          endif

          if (iproc.eq.0) then 
             write(16,'(i5,1x,e12.5,1x,e21.14,a)') itsd,sqrt(fnrm),etot,' SD '
             call wtposout(ncount_cluster,etot,rxyz,at)
             !write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'SD ',ncount_cluster,etot,sqrt(fnrm)
          end if

          fluctoo=flucto
          flucto=fluct
          sumx=0.0_gp 
          sumy=0.0_gp 
          sumz=0.0_gp
          do iat=1,at%nat
             sumx=sumx+ff(1,iat) 
             sumy=sumy+ff(2,iat) 
             sumz=sumz+ff(3,iat)
          end do
          fluct=sumx**2+sumy**2+sumz**2
          if (iproc.eq.0) write(16,*) &
               'fnrm2,flucts',fnrm,sqrt(1.0_gp*at%nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3.0_gp

          !exit statements

          !the norm of the forces is of the same order of the fluctuations
          !or we are much too close to the fixed point from the SD viewpoint
          !or the forces are below the fixed tolerance
          if (fnrm < sqrt(1._gp*at%nat)*(fluct+flucto+fluctoo)*in%frac_fluct/3._gp .or. &
               nsatur > 5 .or. forcemax < in%forcemax) &
               exit loop_sd

          !maximum number of allowed iterations reached
          if (ncount_cluster > in%ncount_cluster_x) then 
             if (iproc.eq.0) write(*,*) 'ncount_cluster in SD2',ncount_cluster
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

          if (beta == in%betax) then 
             if (iproc.eq.0) write(16,*) 'beta=betax'
             care=.true.
          endif
          if (iproc.eq.0) write(16,*) 'beta=',beta
          do iat=1,at%nat
             tpos(1,iat)=rxyz(1,iat)
             tpos(2,iat)=rxyz(2,iat)
             tpos(3,iat)=rxyz(3,iat)
             if ( .not. at%lfrztyp(iat)) then
                if (at%geocode == 'P') then
                   rxyz(1,iat)=modulo(rxyz(1,iat)+beta*ff(1,iat),at%alat1)
                   rxyz(2,iat)=modulo(rxyz(2,iat)+beta*ff(2,iat),at%alat2)
                   rxyz(3,iat)=modulo(rxyz(3,iat)+beta*ff(3,iat),at%alat3)
                else if (at%geocode == 'S') then
                   rxyz(1,iat)=modulo(rxyz(1,iat)+beta*ff(1,iat),at%alat1)
                   rxyz(2,iat)=rxyz(2,iat)+beta*ff(2,iat)
                   rxyz(3,iat)=modulo(rxyz(3,iat)+beta*ff(3,iat),at%alat3)
                else
                   rxyz(1,iat)=rxyz(1,iat)+beta*ff(1,iat)
                   rxyz(2,iat)=rxyz(2,iat)+beta*ff(2,iat)
                   rxyz(3,iat)=rxyz(3,iat)+beta*ff(3,iat)
                end if
             end if
          end do
       end do loop_sd
       exit redo_sd
    end do redo_sd

    if (iproc.eq.0) write(16,*) 'SD FINISHED',iproc

    i_all=-product(shape(tpos))*kind(tpos)
    deallocate(tpos,stat=i_stat)
    call memocc(i_stat,i_all,'tpos',subname)

  end subroutine steepdes

  subroutine detbetax(nproc,iproc,at,pos,rst,in)
    ! determines stepsize betax
    use module_types
    use module_interfaces
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: at
    type(input_variables), intent(inout) :: in
    type(restart_objects), intent(inout) :: rst
    real(gp), dimension(3,at%nat), intent(inout) :: pos
    !local variables
    integer :: nsuc
    real(gp) :: beta0,beta,sum,t1,t2,t3,etotm1,etot0,etotp1,der2,tt
    real(gp), dimension(:,:), allocatable :: tpos,ff,gg

    allocate(tpos(3,at%nat+ndebug),stat=i_stat)
    call memocc(i_stat,tpos,'tpos',subname)
    allocate(ff(3,at%nat+ndebug),stat=i_stat)
    call memocc(i_stat,ff,'ff',subname)
    allocate(gg(3,at%nat+ndebug),stat=i_stat)
    call memocc(i_stat,gg,'gg',subname)

    beta0=abs(in%betax)
    beta=1.e100_gp

    nsuc=0
    loop_detbeta: do
  
       in%inputPsiId=1
       in%output_grid=0
       in%output_wf=.false.
       call call_bigdft(nproc,iproc,at,pos,in,etotm1,ff,rst,infocode)

       ncount_cluster=ncount_cluster+1
       sum=0.0_gp

       do iat=1,at%nat
          tpos(1,iat)=pos(1,iat)
          tpos(2,iat)=pos(2,iat)
          tpos(3,iat)=pos(3,iat)
          t1=beta0*ff(1,iat)
          t2=beta0*ff(2,iat)
          t3=beta0*ff(3,iat)
          sum=sum+t1**2+t2**2+t3**2
          if (.not. at%lfrztyp(iat)) then
             if (at%geocode == 'P') then
                pos(1,iat)=modulo(pos(1,iat)+t1,at%alat1)
                pos(2,iat)=modulo(pos(2,iat)+t2,at%alat2)
                pos(3,iat)=modulo(pos(3,iat)+t3,at%alat3)
             else if (at%geocode == 'S') then
                pos(1,iat)=modulo(pos(1,iat)+t1,at%alat1)
                pos(2,iat)=pos(2,iat)+t2
                pos(3,iat)=modulo(pos(3,iat)+t3,at%alat3)
             else
                pos(1,iat)=pos(1,iat)+t1
                pos(2,iat)=pos(2,iat)+t2
                pos(3,iat)=pos(3,iat)+t3
             end if
          end if
       enddo

       call call_bigdft(nproc,iproc,at,pos,in,etot0,gg,rst,infocode)

       ncount_cluster=ncount_cluster+1

       if (etot0.gt.etotm1) then
          do iat=1,at%nat
             pos(1,iat)=tpos(1,iat)
             pos(2,iat)=tpos(2,iat)
             pos(3,iat)=tpos(3,iat)
          enddo
          beta0=.5_gp*beta0
          if (iproc.eq.0) write(16,*) 'beta0 reset',beta0
          cycle loop_detbeta
       endif
       do iat=1,at%nat
          if (at%lfrztyp(iat)) then
             tpos(1,iat)=pos(1,iat)
             tpos(2,iat)=pos(2,iat)
             tpos(3,iat)=pos(3,iat)
          else
             if (at%geocode == 'P') then
                tpos(1,iat)=modulo(pos(1,iat)+beta0*ff(1,iat),at%alat1)
                tpos(2,iat)=modulo(pos(2,iat)+beta0*ff(2,iat),at%alat2)
                tpos(3,iat)=modulo(pos(3,iat)+beta0*ff(3,iat),at%alat3)
             else if (at%geocode == 'S') then
                tpos(1,iat)=modulo(pos(1,iat)+beta0*ff(1,iat),at%alat1)
                tpos(2,iat)=pos(2,iat)+beta0*ff(2,iat)
                tpos(3,iat)=modulo(pos(3,iat)+beta0*ff(3,iat),at%alat3)
             else
                tpos(1,iat)=pos(1,iat)+beta0*ff(1,iat)
                tpos(2,iat)=pos(2,iat)+beta0*ff(2,iat)
                tpos(3,iat)=pos(3,iat)+beta0*ff(3,iat)
             end if
          end if
       enddo

       call call_bigdft(nproc,iproc,at,tpos,in,etotp1,gg,rst,infocode)

       ncount_cluster=ncount_cluster+1

       if (iproc.eq.0) write(16,'(a,3(1x,e21.14))') 'etotm1,etot0,etotp1',etotm1,etot0,etotp1
       der2=(etotp1+etotm1-2.0_gp*etot0)/sum
       tt=.25_gp/der2
       beta0=.125_gp/der2
       if (iproc.eq.0) write(16,*) 'der2,tt=',der2,tt
       if (der2.gt.0.0_gp) then
          nsuc=nsuc+1
          beta=min(beta,.5_gp/der2)
       endif
       do iat=1,at%nat
          if ( .not. at%lfrztyp(iat)) then
             if (at%geocode == 'P') then
                pos(1,iat)=modulo(pos(1,iat)+tt*ff(1,iat),at%alat1)
                pos(2,iat)=modulo(pos(2,iat)+tt*ff(2,iat),at%alat2)
                pos(3,iat)=modulo(pos(3,iat)+tt*ff(3,iat),at%alat3)
             else if (at%geocode == 'S') then
                pos(1,iat)=modulo(pos(1,iat)+tt*ff(1,iat),at%alat1)
                pos(2,iat)=pos(2,iat)+tt*ff(2,iat)
                pos(3,iat)=modulo(pos(3,iat)+tt*ff(3,iat),at%alat3)
             else
                pos(1,iat)=pos(1,iat)+tt*ff(1,iat)
                pos(2,iat)=pos(2,iat)+tt*ff(2,iat)
                pos(3,iat)=pos(3,iat)+tt*ff(3,iat)
             end if
          end if
       enddo

       if (nsuc < 3) cycle loop_detbeta

       exit loop_detbeta
    end do loop_detbeta
    in%betax=beta

    if (iproc.eq.0) write(16,*) 'betax=',in%betax

    i_all=-product(shape(tpos))*kind(tpos)
    deallocate(tpos,stat=i_stat)
    call memocc(i_stat,i_all,'tpos',subname)
    i_all=-product(shape(ff))*kind(ff)
    deallocate(ff,stat=i_stat)
    call memocc(i_stat,i_all,'ff',subname)
    i_all=-product(shape(gg))*kind(gg)
    deallocate(gg,stat=i_stat)
    call memocc(i_stat,i_all,'gg',subname)

  end subroutine detbetax

end subroutine conjgrad

subroutine wtposout(igeostep,energy,rxyz,atoms)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: igeostep
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !local variables
  character(len=2) :: symbol
  character(len=3) :: suffix
  character(len=3) :: fn
  character(len=10) :: name
  character(len=20) :: filename
  integer :: iat,j,ichg,ispol
  real(gp) :: xmax,ymax,zmax

  write(fn,'(i3.3)') igeostep
  !filename = 'posout_'//fn//'.ascii'
  filename = 'posout_'//fn//'.xyz'
  open(unit=9,file=filename)
  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,atoms%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  write(9,*) atoms%nat,' atomicd0 '!, energy,igeostep
  !write(9,*) xmax+5.d0, 0.d0, ymax+5.d0
  !write(9,*) 0.d0, 0.d0, zmax+5.d0
  if (atoms%geocode == 'P') then
     write(9,'(a,3(1x,1pe21.14))')'periodic',atoms%alat1,atoms%alat2,atoms%alat3
  else if (atoms%geocode == 'S') then
     write(9,'(a,3(1x,1pe21.14))')'surface',atoms%alat1,atoms%alat2,atoms%alat3
  else
     write(9,*)' energy,igeostep ', energy,igeostep
  end if
  do iat=1,atoms%nat
     name=trim(atoms%atomnames(atoms%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
        suffix=name(4:6)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
        suffix=name(3:5)
     else
        symbol=name(1:2)
        suffix=' '
     end if

     call charge_and_spol(atoms%natpol(iat),ichg,ispol)

     !takes into account the blocked atoms and the input polarisation
     if (atoms%lfrztyp(iat) .and. ispol == 0 .and. ichg == 0 ) then
        write(9,'(a2,4x,3(1x,1pe21.14),2x,a4)')symbol,(rxyz(j,iat),j=1,3),'   f'
     else if (atoms%lfrztyp(iat) .and. ispol /= 0 .and. ichg == 0) then
        write(9,'(a2,4x,3(1x,1pe21.14),i7,2x,a4)')symbol,(rxyz(j,iat),j=1,3),&
             ispol,'   f'
     else if (atoms%lfrztyp(iat) .and. ichg /= 0) then
        write(9,'(a2,4x,3(1x,1pe21.14),2(i7),2x,a4)')symbol,(rxyz(j,iat),j=1,3),&
             ispol,ichg,'   f'
     else if (ispol /= 0 .and. ichg == 0) then
        write(9,'(a2,4x,3(1x,1pe21.14),i7)')symbol,(rxyz(j,iat),j=1,3),ispol
     else if (ichg /= 0) then
        write(9,'(a2,4x,3(1x,1pe21.14),2(i7))')symbol,(rxyz(j,iat),j=1,3),ispol,ichg
     else
        write(9,'(a2,4x,3(1x,1pe21.14),2x,a4)')symbol,(rxyz(j,iat),j=1,3)
     end if
  enddo
  close(unit=9)

end subroutine wtposout
