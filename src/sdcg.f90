!> @file
!!  Routines to do Conjugate gradient geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>  Conjugate gradient method
subroutine conjgrad(nproc,iproc,rxyz,at,etot,fxyz,rst,in,ncount_bigdft)
  use module_base
  use module_types
  use module_interfaces
  use minpar
  use yaml_output
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
  real(gp), dimension(6) :: strten
  real(gp), dimension(:,:), allocatable :: tpos,gpf,hh
!  logical::check
  integer :: check
  character(len=4) :: fn4
  character(len=40) :: comment

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
     call close_and_deallocate()
     return
  end if

  !calculate the max of the forces
  call fnrmandforcemax(fxyz,tmp,fmax,at%nat)

  !control whether the convergence criterion is reached after SD
  call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,check) !n(m)
  if (check.gt.5) then
     if (iproc == 0) then
        call yaml_map('Converged before entering CG',ncount_bigdft)
        write(16,*) 'Converged before entering CG',iproc
     end if
     call close_and_deallocate()
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
        call call_bigdft(nproc,iproc,at,tpos,in,tetot,gpf,strten,fnoise,rst,infocode)
!!$        if (iproc == 0) then
!!$           call transforce(at,gpf,sumx,sumy,sumz)
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
!!$        end if
        ncount_bigdft=ncount_bigdft+1

        ! Projection of gradients at beta=0 and beta onto hh
        y0=0._gp
        y1=0._gp
        do iat=1,at%nat
           y0=y0+fxyz(1,iat)*hh(1,iat)+fxyz(2,iat)*hh(2,iat)+fxyz(3,iat)*hh(3,iat)
           y1=y1+gpf(1,iat)*hh(1,iat)+gpf(2,iat)*hh(2,iat)+gpf(3,iat)*hh(3,iat)
        end do
        tt=y0/(y0-y1)
!        if (iproc == 0) then
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
        call call_bigdft(nproc,iproc,at,rxyz,in,etot,fxyz,strten,fnoise,rst,infocode)
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

           if (iproc == 0 .and. parmin%verbosity > 0) then
              write(16,'(a,i5,2(1pe20.12))') 'switching back to SD:etot,etotprev',it,etot,etotprev
              call yaml_open_map('Switching back to SD',flow=.true.)
                 call yaml_map('It',it)
                 call yaml_map('Etot',etot)
                 call yaml_map('Etotprev',etotprev)
              call yaml_close_map()
           end if
           !if (iproc == 0) write(*,'(a,i5,2(1pe20.12))') ' switching back to SD:etot,etotprev',it,etot,etotprev
           do iat=1,at%nat
              rxyz(1,iat)=tpos(1,iat)
              rxyz(2,iat)=tpos(2,iat)
              rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,&
                fnrm,fnoise,in,in%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(fxyz,tmp,fmax,at%nat)
           call convcheck(fmax,fluct*in%frac_fluct, in%forcemax,check) !n(m) 
           if(check.gt.5) then
              if (iproc == 0) write(16,*) 'Converged in switch back SD',iproc
              call close_and_deallocate()
              return
           endif
           cycle redo_cg
        endif

        if (iproc==0) then 
           write(fn4,'(i4.4)') ncount_bigdft
           write(comment,'(a,1pe10.3)')'CONJG:fnrm= ',sqrt(fnrm)
           call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,rxyz,at,trim(comment),forces=fxyz)
        endif

        !if (iproc == 0) write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'CG ',ncount_bigdft,etot,sqrt(fnrm)

        if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)

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
           call yaml_open_map('Geometry')
              if (parmin%verbosity > 0) then
                 ! write(16,'(i5,1x,e12.5,1x,e21.14,a,1x,e9.2)')it,sqrt(fnrm),etot,' GEOPT CG ',beta/in%betax
                 ! write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*in%frac_fluct,fluct
                 write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
                      & ncount_bigdft,it,"GEOPT_CG  ",etot,etot-etotprev, &
                      & fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"b/b0=",beta/in%betax

                 call yaml_map('Ncount_BigDFT',ncount_bigdft)
                 call yaml_map('Iteration',it)
                 call yaml_map('Geometry Method','GEOPT_CG')
                 call yaml_map('etot',(/ etot,etot-etotprev /),fmt='(1pe21.14)')
                 call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct /), fmt='(1pe10.2)')
                 call yaml_map('b/b0', beta/in%betax, fmt='(1pe8.2e1)')
                 !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
                 !     ncount_bigdft,it,"GEOPT_CG  ",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"b/b0=",beta/in%betax
              end if
              !write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
              !     'FORCES norm(Ha/Bohr): maxval=',    fmax,'fnrm=',    fnrm   , 'fluct=',fluct
              call geometry_output(fmax,fnrm,fluct)
           call yaml_close_map()
        end if

        call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,check) !n(m)

        if(check.gt.5) exit loop_cg

        if (ncount_bigdft.gt.in%ncount_cluster_x) then 
           if (iproc==0)  write(16,*) 'SDCG exited before the geometry optimization converged because more than ',&
                            in%ncount_cluster_x,' wavefunction optimizations were required'
           if (iproc == 0) call yaml_map('ncount_bigdft in CG',ncount_bigdft)
           !if (iproc == 0) write(*,*) 'ncount_bigdft in CG',ncount_bigdft
           exit loop_cg
        endif

        if(iproc==0)call timeleft(tt)
        call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,bigdft_mpi%mpi_comm,i_stat)
        if(tt<0)then
           call close_and_deallocate()
           return
        endif

        !if no convergence is reached after 500 CG steps
        !switch back to SD
        if (it == 500) then
           if (iproc == 0) then
              if (parmin%verbosity > 0) write(16,*) &
                   'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
              call yaml_open_map('NO conv in CG after 500 its: switching back to SD')
                 call yaml_map('Iteration',it)
                 call yaml_map('fnrm',fnrm)
                 call yaml_map('Total energy',etot)
              call yaml_close_map()
              !write(*,*) 'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
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

           call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,check) !n(m)
           if(check.gt.5) then
              if (iproc == 0) write(16,*) 'Converged in back up SD',iproc
              call close_and_deallocate()
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

  !! write(6,*) 'CG finished',it,fnrm,etot
  if (iproc == 0) then
     if (parmin%verbosity > 0) then
        write(16,'(1x,a,f8.5,i5)') 'average CG stepsize in terms of betax',avbeta/avnum,iproc
     end if
     call yaml_map('Average CG stepsize in terms of betax',avbeta/avnum,fmt='(f8.5)')
     !write(*,'(1x,a,f8.5,i5)') 'average CG stepsize in terms of betax',avbeta/avnum,iproc
  end if

 call close_and_deallocate()

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


!> Steepest descent method
subroutine steepdes(nproc,iproc,at,rxyz,etot,ff,rst,ncount_bigdft,&
     fnrm,fnoise,in,forcemax_sw,nitsd,fluct)
  use module_base
  use module_types
  use module_interfaces
  use minpar
  use yaml_output
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
  real(gp), dimension(6) :: strten
  real(gp), allocatable, dimension(:,:) :: tpos
  character(len=4) :: fn4
  character(len=40) :: comment

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
        if (iproc == 0) then
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
        call call_bigdft(nproc,iproc,at,rxyz,in,etot,ff,strten,fnoise,rst,infocode)
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
        if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)

        !first and second derivatives of the energy and of the norm of the forces
        !(in units of beta steps)
        de1=etot-etotitm1
        de2=etot-2._gp*etotitm1+etotitm2
        df1=fnrm-fnrmitm1
        df2=fnrm-2._gp*fnrmitm1+fnrmitm2

        if (iproc == 0) then
           if (parmin%verbosity > 0) then
              write(16,'(a,6(1x,e10.3),1x,i2)') 'fmax, fnrm/fnrmitm1, de1<0 , de2>0 , df1<0 , df2>0 ,nsatur',  & 
                & fmax, fnrm/fnrmitm1,de1,de2,df1,df2,nsatur
           end if
           call geometry_output(fmax,fnrm,fluct)
           !write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))') 'FORCES norm(Ha/Bohr): maxval=',    fmax,'fnrm=',    fnrm    ,'fluct=', fluct
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
           call write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,rxyz,at,trim(comment),forces=ff)

           !write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'SD ',ncount_bigdft,etot,sqrt(fnrm)
        end if

        if (iproc == 0) then
           if (parmin%verbosity > 0) then
              write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,I2)') &
                   ncount_bigdft,itsd,"GEOPT_SD  ",etot, etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
                   "b/b0=",beta/in%betax,"nsat=",nsatur
              call yaml_open_map('Geometry')
                 call yaml_map('Ncount_BigDFT',ncount_bigdft)
                 call yaml_map('Iteration',itsd)
                 call yaml_map('Geometry Method','GEOPT_SD')
                 call yaml_map('etot',(/ etot,etot-etotprev /),fmt='(1pe21.14)')
                 call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct /), fmt='(1pe10.2)')
                 call yaml_map('b/b0', beta/in%betax, fmt='(1pe8.2e1)')
                 call yaml_map('nsat',nsatur)
                 call geometry_output(fmax,fnrm,fluct)
              call yaml_close_map()
              !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,I2)') &
              !& ncount_bigdft,itsd,"GEOPT_SD  ",etot, etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct, & 
              !& "b/b0=",beta/in%betax,"nsat=",nsatur
              !write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',&
              !         fnrm,fluct*in%frac_fluct,fluct
           end if
        end if
        etotprev=etot

        !exit statements

        !the norm of the forces is of the same order of the fluctuations
        !or we are much too close to the fixed point from the SD viewpoint
        !or the forces are below the fixed tolerance

        if (iproc==0) then
            !if (nsatur > 5) write(*,'(a,1x,i0)') 'SD EXIT because nsatur > 5 ',nsatur 
            !if (fmax < forcemax_sw ) write(*,'(a,2(1x,1pe24.17))')&
            !     'SD EXIT because fmax < forcemax_sw ', fmax , forcemax_sw
            if (nsatur > 5) call yaml_map('SD EXIT because nsatur > 5 ',nsatur)
            if (fmax < forcemax_sw ) &
               & call yaml_map('SD EXIT because fmax < forcemax_sw ', (/ fmax , forcemax_sw /), fmt='(1pe24.17)')
            if (fmax < fluct*in%frac_fluct ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < fluctuation ' ,fmax , fluct*in%frac_fluct
!            if (fmax < fluct*in%forcemax ) write(16,'(a,2(1x,1pe24.17))')&
            if (fmax < in%forcemax ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < forcemax ' ,fmax , in%forcemax
        endif

        if ( nsatur > 5 .or. fmax < max(forcemax_sw,in%forcemax,fluct*in%frac_fluct))  exit loop_sd

        !maximum number of allowed iterations reached
        if (ncount_bigdft > in%ncount_cluster_x) then 
           if (iproc==0)  write(16,*) 'SDCG exited before the geometry optimization converged because more than ',& 
                            in%ncount_cluster_x,' wavefunction optimizations were required'
           if (iproc == 0) call yaml_map('ncount_bigdft in SD2', (/ ncount_bigdft,in%ncount_cluster_x /))
           !if (iproc == 0) write(*,*) 'ncount_bigdft in SD2',ncount_bigdft,in%ncount_cluster_x
           exit loop_sd
        endif

        !maximum number of allowed SD steps reached, locally or globally
        if (itsd >= nitsd) then 
           if (iproc == 0) write(16,'(a,i5,1x,e10.3,1x,e21.14)') &
                'SD: NO CONVERGENCE:itsd,fnrm2,etot',itsd,fnrm,etot
           exit loop_sd
        endif
        if (itot >= nitsd) then
           if (iproc == 0) write(16,'(a,i5,i5,1x,e10.3,1x,e21.14)') &
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
           !     if (iproc == 0) write(16,*) 'beta=betax'
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
  use module_interfaces
  use minpar
  use yaml_output
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
  real(gp) :: fluct,fnrm,fnrmold,beta,betaxx,betalast !n(c) anoise,betalastold 
  real(gp) :: etotold,fmax,scpr,curv,tt,etotprev
  real(gp), dimension(6) :: strten
  real(gp), dimension(:,:), allocatable :: posold,ffold
  !n(c) logical :: reset!,check
  integer :: check
  character(len=4) :: fn4
  character(len=40) :: comment

  check=0
  etotprev=etot
  allocate(posold(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,posold,'posold',subname)
  allocate(ffold(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,ffold,'ffold',subname)

  !n(c) anoise=1.e-4_gp
  fluct=0._gp

  beta=in%betax
  itsd=0
  in%inputPsiId=1
  call call_bigdft(nproc,iproc,at,wpos,in,etotold,ffold,strten,fnoise,rst,infocode)
  call fnrmandforcemax(ffold,fnrm,fmax,at%nat)   
  if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)
  if (iproc == 0) then
     if (parmin%verbosity > 0) then
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
          & ncount_bigdft,itsd,"GEOPT_VSSD",etotold,etotold-etotprev,fmax,sqrt(fnrm), &
          & fluct*in%frac_fluct,fluct,"beta=",beta
        call yaml_open_map('Geometry')
           call yaml_map('Ncount_BigDFT',ncount_bigdft)
           call yaml_map('Iteration',itsd)
           call yaml_map('Geometry Method','GEOPT_VSSD')
           call yaml_map('etotold',(/ etotold,etotold-etotprev /),fmt='(1pe21.14)')
           call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct /), fmt='(1pe10.2)')
           call yaml_map('beta', beta, fmt='(1pe8.2e1)')
        call yaml_close_map()
        !if (parmin%verbosity > 0)   write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
        !&ncount_bigdft,itsd,"GEOPT_VSSD",etotold,etotold-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,"beta=",beta
     end if
     etotprev=etotold
!!$     call transforce(at,ffold,sumx,sumy,sumz)                         
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  

     write(fn4,'(i4.4)') ncount_bigdft
     write(comment,'(a,1pe10.3)')'Initial VSSD:fnrm= ',sqrt(fnrm)
     call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etotold,wpos,at,trim(comment),forces=ffold)
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
  !n(c) reset=.true.
  !n(c) betalastold=in%betax

  nitsd=1000
  loop_ntsd: do itsd=1,nitsd
     in%inputPsiId=1
     call call_bigdft(nproc,iproc,at,wpos,in,etot,ff,strten,fnoise,rst,infocode)
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
     if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)
!     if (iproc==0) write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',fnrm,fluct*in%frac_fluct,fluct
     call convcheck(fmax,fluct*in%frac_fluct, in%forcemax,check) !n(m)
     if (check > 5) exit loop_ntsd
     if (ncount_bigdft >= in%ncount_cluster_x) then 
        if (iproc==0)  write(16,*) 'VSSD exited before the geometry optimization converged because more than ',& 
             in%ncount_cluster_x,' wavefunction optimizations were required'
        exit loop_ntsd
     endif


     if (etot > etotold) then
        !n(c) reset=.true.
        beta=in%betax
        if (iproc == 0) write(16,*) 'new positions rejected, reduced beta',beta
        do iat=1,at%nat
           wpos(1,iat)=posold(1,iat)+beta*ffold(1,iat)
           wpos(2,iat)=posold(2,iat)+beta*ffold(2,iat)
           wpos(3,iat)=posold(3,iat)+beta*ffold(3,iat)
        enddo
        !call atomic_axpy(at,posold,beta,ffold,wpos)
     else
        !n(c) reset=.false.
        if (betalast.gt.0.d0) then
           beta=max(min(beta*1.5d0,betalast),in%betax)
        else
           beta=1.25d0*beta
        endif

        if (iproc == 0) then
           write(fn4,'(i4.4)') ncount_bigdft-1
           write(comment,'(a,1pe10.3)')'VSSD:fnrm= ',sqrt(fnrm)
           call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,wpos,at,trim(comment),forces=ff)
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
        !n(c) betalastold=betalast
     endif
  
     if (iproc == 0.and.parmin%verbosity > 0) then
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
        ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
        "beta=",beta,"last beta=",betalast
     end if

     if (iproc == 0.and.parmin%verbosity > 0) then
        call yaml_open_map('Geometry')
           call yaml_map('Ncount_BigDFT',ncount_bigdft)
           call yaml_map('Iteration',itsd)
           call yaml_map('Geometry Method','GEOPT_VSSD')
           call yaml_map('etot',(/ etot,etot-etotprev /),fmt='(1pe21.14)')
           call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct /), fmt='(1pe10.2)')
           call yaml_map('beta', beta, fmt='(1pe8.2e1)')
           call yaml_map('last beta', betalast, fmt='(1pe8.2e1)')
           call geometry_output(fmax,fnrm,fluct)
        call yaml_close_map()
        !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
        !ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
        !&"beta=",beta,"last beta=",betalast
    end if

    etotprev=etot
!     if (iproc == 0) write(16,'(i5,1x,e12.5,1x,e21.14,a,e10.3,1x,e10.3)') itsd,sqrt(fnrm),etot,' GEOPT VSSD ',beta,betalast
     if(iproc==0)call timeleft(tt)
     call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,bigdft_mpi%mpi_comm,i_stat)
     if(tt<0) exit loop_ntsd


  enddo loop_ntsd

  if (iproc == 0.and.parmin%verbosity > 0) & 
       write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
       ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
       "beta=",beta,"last beta=",betalast

  if (iproc == 0.and.parmin%verbosity > 0) then
     call yaml_open_map('Geometry')
        call yaml_map('Ncount_BigDFT',ncount_bigdft)
        call yaml_map('Iteration',itsd)
        call yaml_map('Geometry Method','GEOPT_VSSD')
        call yaml_map('etot',(/ etot,etot-etotprev /),fmt='(1pe21.14)')
        call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct /), fmt='(1pe10.2)')
        call yaml_map('beta', beta, fmt='(1pe8.2e1)')
        call yaml_map('last beta', betalast, fmt='(1pe8.2e1)')
        call geometry_output(fmax,fnrm,fluct)
     call yaml_close_map()
     !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
     !&ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct,& 
     !&"beta=",beta,"last beta=",betalast
  end if

  if (iproc == 0 .and. itsd == nitsd+1) &
       write(16,'(a,i5,e9.2,e18.10,e9.2)') '---- SD FAILED  TO CONVERGE'
  if (iproc == 0) then
     if (parmin%verbosity > 0) then
        write(16,'(a,i5,e9.2,e18.10)') 'variable stepsize SD FINISHED,iter, force norm,energy',itsd,sqrt(fnrm),etot
        write(16,'(a,e9.2)') 'suggested value for stepsize:', betaxx
     end if
     write(fn4,'(i4.4)') ncount_bigdft-1
     write(comment,'(a,1pe10.3)')'VSSD:fnrm= ',sqrt(fnrm)
     call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,wpos,at,trim(comment),forces=ff)
  endif


  i_all=-product(shape(posold))*kind(posold)
  deallocate(posold,stat=i_stat)
  call memocc(i_stat,i_all,'posold',subname)
  i_all=-product(shape(ffold))*kind(ffold)
  deallocate(ffold,stat=i_stat)
  call memocc(i_stat,i_all,'ffold',subname)


END SUBROUTINE vstepsd
