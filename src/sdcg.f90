!> @file
!!  Routines to do Conjugate gradient geometry optmization
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>  Conjugate gradient method
subroutine conjgrad(runObj,outs,nproc,iproc,ncount_bigdft)
  use module_base
  use module_types
  use module_interfaces
  use minpar
  use yaml_output
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  !local variables
  type(DFT_global_output) :: l_outs
  character(len=*), parameter :: subname='conjgrad'  
  integer :: nfail,it,iat,i_all,i_stat,infocode, nitsd
  real(gp) :: anoise,fluct,avbeta,avnum,fnrm,etotprev,beta0,beta
  real(gp) :: y0,y1,tt,oben1,oben2,oben,unten,rlambda,fmax,tmp!,eprev
  real(gp), dimension(:,:), allocatable :: tpos,hh
!  logical::check
  integer :: check
  character(len=4) :: fn4
  character(len=40) :: comment

  check=0
  tpos = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='tpos')
  hh = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='hh')
  call init_global_output(l_outs, runObj%atoms%astruct%nat)

  anoise=1.e-4_gp
  fluct=0._gp

  avbeta=0._gp
  avnum=0._gp
  nfail=0
  nitsd=500

  !start with a steepest descent algorithm
  call steepdes(runObj,outs,nproc,iproc,ncount_bigdft,fnrm,runObj%inputs%forcemax,nitsd,fluct)
  if (ncount_bigdft >= runObj%inputs%ncount_cluster_x) then
     if (iproc==0 .and. parmin%verbosity > 0) &
        & write(16,*) 'SDCG exited before the geometry optimization converged because more than ',&
          runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
     call close_and_deallocate()
     return
  end if

  !calculate the max of the forces
  call fnrmandforcemax(outs%fxyz,tmp,fmax,outs%fdim)

  !control whether the convergence criterion is reached after SD
  call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,check) !n(m)
  if (check.gt.5) then
     if (iproc == 0) then
        call yaml_map('Converged before entering CG',ncount_bigdft)
        write(16,*) 'Converged before entering CG',iproc
     end if
     call close_and_deallocate()
     return
  endif

  redo_cg: do
     call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, hh(1,1), 1)

     beta0=4._gp*runObj%inputs%betax
     it=0
     loop_cg: do
        it=it+1
        ! We save rxyz in tpos
        call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, tpos(1,1), 1)

        !C line minimize along hh ----
        !call atomic_axpy(at,rxyz,beta0,hh,tpos)
        call axpy(3 * runObj%atoms%astruct%nat, beta0, hh(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)

        runObj%inputs%inputPsiId=1
        call call_bigdft(runObj,l_outs,nproc,iproc,infocode)
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
        do iat=1,outs%fdim
           y0=y0+outs%fxyz(1,iat)*hh(1,iat)+outs%fxyz(2,iat)*hh(2,iat)+outs%fxyz(3,iat)*hh(3,iat)
           y1=y1+l_outs%fxyz(1,iat)*hh(1,iat)+l_outs%fxyz(2,iat)*hh(2,iat)+l_outs%fxyz(3,iat)*hh(3,iat)
        end do
        tt=y0/(y0-y1)
!        if (iproc == 0) then
!           if (parmin%verbosity > 0) &
!                & write(16,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
!           write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
!        end if

        beta=beta0*max(min(tt,2._gp),-.25_gp)
        
        ! We put back rxyz to original
        call vcopy(3 * runObj%atoms%astruct%nat, tpos(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)

        !call atomic_axpy(at,rxyz,beta,hh,rxyz)
        call axpy(3 * runObj%atoms%astruct%nat, beta, hh(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)

        avbeta=avbeta+beta/runObj%inputs%betax
        avnum=avnum+1._gp
        !if (iproc ==0)print *,'beta,avbeta,avnum',beta,avbeta,avnum 
        !C new gradient
        do iat=1,outs%fdim
           l_outs%fxyz(1,iat)=outs%fxyz(1,iat)
           l_outs%fxyz(2,iat)=outs%fxyz(2,iat)
           l_outs%fxyz(3,iat)=outs%fxyz(3,iat)
        end do

        etotprev=outs%energy
        call call_bigdft(runObj,outs,nproc,iproc,infocode)
!!$        if (iproc == 0) then
!!$           call transforce(at,fxyz,sumx,sumy,sumz)
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy
!!$           write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz
!!$        end if
        ncount_bigdft=ncount_bigdft+1
        !if the energy goes up (a small tolerance of anoise is allowed)
        !switch back to SD
        if (outs%energy > etotprev+anoise) then

           if (iproc == 0 .and. parmin%verbosity > 0) then
              write(16,'(a,i5,2(1pe20.12))') 'switching back to SD:etot,etotprev',it,outs%energy,etotprev
              call yaml_open_map('Switching back to SD',flow=.true.)
                 call yaml_map('It',it)
                 call yaml_map('Etot',outs%energy)
                 call yaml_map('Etotprev',etotprev)
              call yaml_close_map()
           end if
           !if (iproc == 0) write(*,'(a,i5,2(1pe20.12))') ' switching back to SD:etot,etotprev',it,etot,etotprev
           do iat=1,runObj%atoms%astruct%nat
              runObj%atoms%astruct%rxyz(1,iat)=tpos(1,iat)
              runObj%atoms%astruct%rxyz(2,iat)=tpos(2,iat)
              runObj%atoms%astruct%rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(runObj,outs,nproc,iproc,ncount_bigdft,fnrm,runObj%inputs%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(outs%fxyz,tmp,fmax,outs%fdim)
           call convcheck(fmax,fluct*runObj%inputs%frac_fluct, runObj%inputs%forcemax,check) !n(m) 
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
           call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
                & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),&
                forces=outs%fxyz)
        endif

        !if (iproc == 0) write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'CG ',ncount_bigdft,etot,sqrt(fnrm)

        if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct) !n(m)

!!$        call atomic_dot(at,gpf,gpf,unten)
!!$        call atomic_dot(at,gpf,fxyz,oben1)
!!$        call atomic_dot(at,fxyz,fxyz,oben2)
        unten=dot(3*outs%fdim,l_outs%fxyz(1,1),1,l_outs%fxyz(1,1),1)
        oben1=dot(3*outs%fdim,l_outs%fxyz(1,1),1,outs%fxyz(1,1),1)
        oben2=dot(3*outs%fdim,outs%fxyz(1,1),1,outs%fxyz(1,1),1)

        oben=oben2-oben1

!!!        obenx=0._gp
!!!        obeny=0._gp
!!!        obenz=0._gp
!!!        unten=0._gp
!!!        do iat=1,outs%fdim
!!!           obenx=obenx+(fxyz(1,iat)-gpf(1,iat))*fxyz(1,iat)
!!!           obeny=obeny+(fxyz(2,iat)-gpf(2,iat))*fxyz(2,iat)
!!!           obenz=obenz+(fxyz(3,iat)-gpf(3,iat))*fxyz(3,iat)
!!!           unten=unten+gpf(1,iat)**2+gpf(2,iat)**2+gpf(3,iat)**2
!!!        end do

        call fnrmandforcemax(outs%fxyz,fnrm,fmax,outs%fdim)
        if (iproc == 0) then
           call yaml_open_map('Geometry')
              if (parmin%verbosity > 0) then
                 ! write(16,'(i5,1x,e12.5,1x,e21.14,a,1x,e9.2)')it,sqrt(fnrm),etot,' GEOPT CG ',beta/runObj%inputs%betax
                 ! write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct', fnrm,fluct*runObj%inputs%frac_fluct,fluct
                 write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
                      ncount_bigdft,it,"GEOPT_CG  ",outs%energy,outs%energy-etotprev, &
                      & fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,&
                      fluct,"b/b0=",beta/runObj%inputs%betax
                 
                 call yaml_map('Ncount_BigDFT',ncount_bigdft)
                 call yaml_map('Iteration',it)
                 call yaml_map('Geometry Method','GEOPT_CG')
                 call yaml_map('etot',(/ outs%energy,outs%energy-etotprev /),fmt='(1pe21.14)')
                 call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct /), fmt='(1pe10.2)')
                 call yaml_map('b/b0', beta/runObj%inputs%betax, fmt='(1pe8.2e1)')
                 !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
                 !     ncount_bigdft,it,"GEOPT_CG  ",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,&
                 !     fluct,"b/b0=",beta/runObj%inputs%betax
              end if
              !write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
              !     'FORCES norm(Ha/Bohr): maxval=',    fmax,'fnrm=',    fnrm   , 'fluct=',fluct
              call geometry_output(fmax,fnrm,fluct)
           call yaml_close_map()
        end if

        call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,check) !n(m)

        if(check.gt.5) exit loop_cg

        if (ncount_bigdft.gt.runObj%inputs%ncount_cluster_x) then 
           if (iproc==0)  write(16,*) 'SDCG exited before the geometry optimization converged because more than ',&
                            runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
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
                   'NO conv in CG after 500 its: switching back to SD',it,fnrm,outs%energy
              call yaml_open_map('NO conv in CG after 500 its: switching back to SD')
                 call yaml_map('Iteration',it)
                 call yaml_map('fnrm',fnrm)
                 call yaml_map('Total energy',outs%energy)
              call yaml_close_map()
              !write(*,*) 'NO conv in CG after 500 its: switching back to SD',it,fnrm,etot
           end if
           do iat=1,runObj%atoms%astruct%nat
              runObj%atoms%astruct%rxyz(1,iat)=tpos(1,iat)
              runObj%atoms%astruct%rxyz(2,iat)=tpos(2,iat)
              runObj%atoms%astruct%rxyz(3,iat)=tpos(3,iat)
           end do

           call steepdes(runObj,outs,nproc,iproc,ncount_bigdft,fnrm,runObj%inputs%forcemax,nitsd,fluct)

           !calculate the max of the forces
           call fnrmandforcemax(outs%fxyz,tmp,fmax,outs%fdim)

           call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,check) !n(m)
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
        do iat=1,outs%fdim
           hh(1,iat)=outs%fxyz(1,iat)+rlambda*hh(1,iat)
           hh(2,iat)=outs%fxyz(2,iat)+rlambda*hh(2,iat)
           hh(3,iat)=outs%fxyz(3,iat)+rlambda*hh(3,iat)
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
    call f_free(tpos)
    call f_free(hh)
    call deallocate_global_output(l_outs)
  END SUBROUTINE close_and_deallocate

END SUBROUTINE conjgrad


!> Steepest descent method
subroutine steepdes(runObj,outs,nproc,iproc,ncount_bigdft,fnrm,forcemax_sw,nitsd,fluct)
  use module_base
  use module_types
  use module_interfaces
  use minpar
  use yaml_output
  !use module_interfaces
  implicit none
  integer, intent(in) :: nproc,iproc,nitsd
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  integer, intent(inout) :: ncount_bigdft
  real(gp), intent(out) :: fnrm
  real(gp), intent(inout) :: fluct
  real(gp), intent(in)::forcemax_sw
  !local variables
  character(len=*), parameter :: subname='steepdes'
  logical :: care,move_this_coordinate
  integer :: nsatur,iat,itot,itsd,i_stat,i_all,infocode,nbeqbx,i,ixyz,nr
  real(gp) :: etotitm2,fnrmitm2,etotitm1,fnrmitm1,anoise
  real(gp) :: fmax,de1,de2,df1,df2,beta,etotprev
  real(gp), allocatable, dimension(:,:) :: tpos
  character(len=4) :: fn4
  character(len=40) :: comment

  tpos = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='tpos')

  etotprev=outs%energy
  anoise=0.e-4_gp
  fluct=0._gp

  beta=runObj%inputs%betax
  care=.true.
  nbeqbx=0
  nsatur=0
  etotitm2=1.e100_gp
  fnrmitm2=1.e100_gp
  etotitm1=1.e100_gp
  fnrmitm1=1.e100_gp

  call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, tpos(1,1), 1)

  itot=0

  nr=0
  do i=1,3*runObj%atoms%astruct%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(runObj%atoms%astruct%ifrztyp(iat),ixyz)) nr=nr+1
  enddo


  redo_sd: do
     if (ncount_bigdft.gt.runObj%inputs%ncount_cluster_x) then 
        if (iproc == 0) then
           write(*,'(a,i6,i6,i6)') 'SD FINISHED because ncount_bigdft > ncount_cluster_x', &
                & iproc,ncount_bigdft,runObj%inputs%ncount_cluster_x
           if (parmin%verbosity > 0) then
              write(16,'(a,i6,i6,i6)') &
                   'SD FINISHED because ncount_bigdft > ncount_cluster_x',iproc,ncount_bigdft,runObj%inputs%ncount_cluster_x
              write(16,'(a,i6,a)') 'SDCG exited before the geometry optimization converged because more than ', & 
                   runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
           end if
        end if

        call f_free(tpos)
        return
     endif

     itsd=0
     loop_sd: do
        itsd=itsd+1
        itot=itot+1

        runObj%inputs%inputPsiId=1
        call call_bigdft(runObj,outs,nproc,iproc,infocode)
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
        if (care .and. outs%energy > etotitm1+anoise) then
           call vcopy(3 * runObj%atoms%astruct%nat, tpos(1,1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
           beta=.5_gp*beta
           if (iproc == 0) write(16,'(a,1x,e9.2,1x,i5,2(1x,e21.14))') &
                'SD reset, beta,itsd,etot,etotitm1= ',beta,itsd,outs%energy,etotitm1
           if (beta <= 1.d-1*runObj%inputs%betax) then
              if (iproc == 0) write(16,*) &
                   'beta getting too small, do not care anymore if energy goes up'
              care=.false.
           endif
           cycle redo_sd
        endif

        call fnrmandforcemax(outs%fxyz,fnrm,fmax,outs%fdim)
        if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct) !n(m)

        !first and second derivatives of the energy and of the norm of the forces
        !(in units of beta steps)
        de1=outs%energy-etotitm1
        de2=outs%energy-2._gp*etotitm1+etotitm2
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
        if (care .and. itsd >= 3 .and. beta == runObj%inputs%betax .and. &
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
           call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
                & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),&
                forces=outs%fxyz)

           !write(17,'(a,i5,1x,e17.10,1x,e9.2)') 'SD ',ncount_bigdft,etot,sqrt(fnrm)
        end if

        if (iproc == 0) then
           if (parmin%verbosity > 0) then
              write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,I2)') &
                   ncount_bigdft,itsd,"GEOPT_SD  ",outs%energy, outs%energy-etotprev, &
                   & fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,& 
                   "b/b0=",beta/runObj%inputs%betax,"nsat=",nsatur
              call yaml_open_map('Geometry')
                 call yaml_map('Ncount_BigDFT',ncount_bigdft)
                 call yaml_map('Iteration',itsd)
                 call yaml_map('Geometry Method','GEOPT_SD')
                 call yaml_map('etot',(/ outs%energy,outs%energy-etotprev /),fmt='(1pe21.14)')
                 call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct /), fmt='(1pe10.2)')
                 call yaml_map('b/b0', beta/runObj%inputs%betax, fmt='(1pe8.2e1)')
                 call yaml_map('nsat',nsatur)
                 call geometry_output(fmax,fnrm,fluct)
              call yaml_close_map()
              !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,I2)') &
              !& ncount_bigdft,itsd,"GEOPT_SD  ",etot, etot-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct, & 
              !& "b/b0=",beta/runObj%inputs%betax,"nsat=",nsatur
              !write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',&
              !         fnrm,fluct*runObj%inputs%frac_fluct,fluct
           end if
        end if
        etotprev=outs%energy

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
            if (fmax < fluct*runObj%inputs%frac_fluct ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < fluctuation ' ,fmax , fluct*runObj%inputs%frac_fluct
!            if (fmax < fluct*runObj%inputs%forcemax ) write(16,'(a,2(1x,1pe24.17))')&
            if (fmax < runObj%inputs%forcemax ) write(16,'(a,2(1x,1pe24.17))')&
                 'SD EXIT because fmax < forcemax ' ,fmax , runObj%inputs%forcemax
        endif

        if ( nsatur > 5 .or. fmax < max(forcemax_sw,runObj%inputs%forcemax,fluct*runObj%inputs%frac_fluct))  exit loop_sd

        !maximum number of allowed iterations reached
        if (ncount_bigdft > runObj%inputs%ncount_cluster_x) then 
           if (iproc==0)  write(16,*) 'SDCG exited before the geometry optimization converged because more than ',& 
                            runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
           if (iproc == 0) call yaml_map('ncount_bigdft in SD2', (/ ncount_bigdft,runObj%inputs%ncount_cluster_x /))
           !if (iproc == 0) write(*,*) 'ncount_bigdft in SD2',ncount_bigdft,runObj%inputs%ncount_cluster_x
           exit loop_sd
        endif

        !maximum number of allowed SD steps reached, locally or globally
        if (itsd >= nitsd) then 
           if (iproc == 0) write(16,'(a,i5,1x,e10.3,1x,e21.14)') &
                'SD: NO CONVERGENCE:itsd,fnrm2,etot',itsd,fnrm,outs%energy
           exit loop_sd
        endif
        if (itot >= nitsd) then
           if (iproc == 0) write(16,'(a,i5,i5,1x,e10.3,1x,e21.14)') &
                'SD: NO CONVERGENCE:itsd,itot,fnrm2,etot:',itsd,itot,fnrm,outs%energy
           exit loop_sd
        endif

        etotitm2=etotitm1
        etotitm1=outs%energy
        fnrmitm2=fnrmitm1
        fnrmitm1=fnrm
        
        beta=min(1.2_gp*beta,runObj%inputs%betax)
!!!        if (beta /= runObj%inputs%betax) nbeqbx=0
        if (beta == runObj%inputs%betax) then 
           !     if (iproc == 0) write(16,*) 'beta=betax'
           care=.true.
!!!           !if beta=betax since too many iterations (say 5),
!!!           !then betax can be increased
!!!           nbeqbx=nbeqbx+1
!!!           if (nbeqbx == 5) then
!!!              runObj%inputs%betax=1.2_gp*runObj%inputs%betax
!!!              nbeqbx=0
!!!           end if
        endif
!        if (iproc == 0 .and. parmrunObj%inputs%verbosity > 0) write(16,*) 'beta=',beta

        call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, tpos(1,1), 1)
        !call atomic_axpy(at,rxyz,beta,ff,rxyz)
        call axpy(3 * runObj%atoms%astruct%nat,beta,outs%fxyz(1,1),1,runObj%atoms%astruct%rxyz(1,1),1)

!!!        do iat=1,runObj%atoms%nat
!!!           tpos(1,iat)=rxyz(1,iat)
!!!           tpos(2,iat)=rxyz(2,iat)
!!!           tpos(3,iat)=rxyz(3,iat)
!!!           if ( .not. runObj%atoms%lfrztyp(iat)) then
!!!              if (runObj%atoms%geocode == 'P') then
!!!                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*ff(1,iat),runObj%atoms%alat1)
!!!                 rxyz(2,iat)=modulo(rxyz(2,iat)+beta*ff(2,iat),runObj%atoms%alat2)
!!!                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*ff(3,iat),runObj%atoms%alat3)
!!!              else if (runObj%atoms%geocode == 'S') then
!!!                 rxyz(1,iat)=modulo(rxyz(1,iat)+beta*ff(1,iat),runObj%atoms%alat1)
!!!                 rxyz(2,iat)=rxyz(2,iat)+beta*ff(2,iat)
!!!                 rxyz(3,iat)=modulo(rxyz(3,iat)+beta*ff(3,iat),runObj%atoms%alat3)
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

  call f_free(tpos)

END SUBROUTINE steepdes


!> Variable step steepest descent
subroutine vstepsd(runObj,outs,nproc,iproc,ncount_bigdft)
  use module_base
  use module_types
  use module_interfaces
  use minpar
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  !local variables
  character(len=*), parameter :: subname='vstepsd'  
  type(DFT_global_output) :: outsold
  integer :: iat,i_all,i_stat,infocode, nitsd,itsd
  real(gp) :: fluct,fnrm,fnrmold,beta,betaxx,betalast !n(c) anoise,betalastold 
  real(gp) :: fmax,scpr,curv,tt,etotprev
  real(gp), dimension(:,:), allocatable :: posold
  !n(c) logical :: reset!,check
  integer :: check
  character(len=4) :: fn4
  character(len=40) :: comment

  check=0
  etotprev=outs%energy
  posold = f_malloc((/ 3, runObj%atoms%astruct%nat /),id='posold')
  call init_global_output(outsold, runObj%atoms%astruct%nat)

  !n(c) anoise=1.e-4_gp
  fluct=0._gp

  beta=runObj%inputs%betax
  itsd=0
  runObj%inputs%inputPsiId=1
  call call_bigdft(runObj,outsold,nproc,iproc,infocode)
  call fnrmandforcemax(outsold%fxyz,fnrm,fmax,outsold%fdim)
  if (fmax < 3.d-1) call updatefluctsum(outsold%fnoise,fluct) !n(m)
  if (iproc == 0) then
     if (parmin%verbosity > 0) then
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
          & ncount_bigdft,itsd,"GEOPT_VSSD",outsold%energy,outsold%energy-etotprev,fmax,sqrt(fnrm), &
          & fluct*runObj%inputs%frac_fluct,fluct,"beta=",beta
        call yaml_open_map('Geometry')
           call yaml_map('Ncount_BigDFT',ncount_bigdft)
           call yaml_map('Iteration',itsd)
           call yaml_map('Geometry Method','GEOPT_VSSD')
           call yaml_map('etotold',(/ outsold%energy,outsold%energy-etotprev /),fmt='(1pe21.14)')
           call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct /), fmt='(1pe10.2)')
           call yaml_map('beta', beta, fmt='(1pe8.2e1)')
        call yaml_close_map()
        !if (parmrunObj%inputs%verbosity > 0)   write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1)') &
        !&ncount_bigdft,itsd,"GEOPT_VSSD",etotold,etotold-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,"beta=",beta
     end if
     etotprev=outsold%energy
!!$     call transforce(at,ffold,sumx,sumy,sumz)                         
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  

     write(fn4,'(i4.4)') ncount_bigdft
     write(comment,'(a,1pe10.3)')'Initial VSSD:fnrm= ',sqrt(fnrm)
     call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
          & outsold%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),&
          forces=outsold%fxyz)
!     if (parmin%verbosity > 0) &
!          & write(16,'(1x,e12.5,1x,e21.14,a,e10.3)')sqrt(fnrm),etotold,' GEOPT VSSD ',beta
  end if

  ncount_bigdft=ncount_bigdft+1


  fnrmold=0.0_gp
  do iat=1,outsold%fdim
     fnrmold=fnrmold+outsold%fxyz(1,iat)**2+outsold%fxyz(2,iat)**2+outsold%fxyz(3,iat)**2
     posold(1,iat)=runObj%atoms%astruct%rxyz(1,iat)
     posold(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
     posold(3,iat)=runObj%atoms%astruct%rxyz(3,iat)
     runObj%atoms%astruct%rxyz(1,iat)=runObj%atoms%astruct%rxyz(1,iat)+beta*outsold%fxyz(1,iat)
     runObj%atoms%astruct%rxyz(2,iat)=runObj%atoms%astruct%rxyz(2,iat)+beta*outsold%fxyz(2,iat)
     runObj%atoms%astruct%rxyz(3,iat)=runObj%atoms%astruct%rxyz(3,iat)+beta*outsold%fxyz(3,iat)
  enddo
  !call atomic_dot(at,ffold,ffold,fnrmold)
  !call atomic_axpy(at,wpos,beta,ffold,wpos)
  betaxx=1.d100
  !n(c) reset=.true.
  !n(c) betalastold=runObj%inputs%betax

  nitsd=1000
  loop_ntsd: do itsd=1,nitsd
     runObj%inputs%inputPsiId=1
     call call_bigdft(runObj,outs,nproc,iproc,infocode)
!!$     if (iproc == 0) then                                        
!!$        call transforce(at,outs%fxyz,sumx,sumy,sumz)                         
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$        write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
!!$     end if
     ncount_bigdft=ncount_bigdft+1
     fnrm=0.d0
     scpr=0.d0
     do iat=1,outs%fdim
        fnrm=fnrm+outs%fxyz(1,iat)**2+outs%fxyz(2,iat)**2+outs%fxyz(3,iat)**2
        scpr=scpr+outsold%fxyz(1,iat)*outs%fxyz(1,iat)+outsold%fxyz(2,iat)*outs%fxyz(2,iat)+outsold%fxyz(3,iat)*outs%fxyz(3,iat)
     enddo
!!$     call  atomic_dot(at,ff,ff,fnrm)
!!$     call  atomic_dot(at,ff,ffold,scpr)
     curv=(fnrmold-scpr)/(beta*fnrmold)
     betalast=.5d0/curv
     if (betalast.gt.0.d0) betaxx=min(betaxx,1.5d0*betalast)
     call fnrmandforcemax(outs%fxyz,fnrm,fmax,outs%fdim)
     if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct) !n(m)
!     if (iproc==0) write(16,'(1x,a,3(1x,1pe14.5))') 'fnrm2,fluct*frac_fluct,fluct',fnrm,fluct*runObj%inputs%frac_fluct,fluct
     call convcheck(fmax,fluct*runObj%inputs%frac_fluct, runObj%inputs%forcemax,check) !n(m)
     if (check > 5) exit loop_ntsd
     if (ncount_bigdft >= runObj%inputs%ncount_cluster_x) then 
        if (iproc==0)  write(16,*) 'VSSD exited before the geometry optimization converged because more than ',& 
             runObj%inputs%ncount_cluster_x,' wavefunction optimizations were required'
        exit loop_ntsd
     endif


     if (outs%energy > outsold%energy) then
        !n(c) reset=.true.
        beta=runObj%inputs%betax
        if (iproc == 0) write(16,*) 'new positions rejected, reduced beta',beta
        do iat=1,runObj%atoms%astruct%nat
           runObj%atoms%astruct%rxyz(1,iat)=posold(1,iat)+beta*outsold%fxyz(1,iat)
           runObj%atoms%astruct%rxyz(2,iat)=posold(2,iat)+beta*outsold%fxyz(2,iat)
           runObj%atoms%astruct%rxyz(3,iat)=posold(3,iat)+beta*outsold%fxyz(3,iat)
        enddo
        !call atomic_axpy(at,posold,beta,ffold,wpos)
     else
        !n(c) reset=.false.
        if (betalast.gt.0.d0) then
           beta=max(min(beta*1.5d0,betalast),runObj%inputs%betax)
        else
           beta=1.25d0*beta
        endif

        if (iproc == 0) then
           write(fn4,'(i4.4)') ncount_bigdft-1
           write(comment,'(a,1pe10.3)')'VSSD:fnrm= ',sqrt(fnrm)
           call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
                & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),&
                forces=outs%fxyz)
        endif

        do iat=1,runObj%atoms%astruct%nat
           posold(1,iat)=runObj%atoms%astruct%rxyz(1,iat)
           posold(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
           posold(3,iat)=runObj%atoms%astruct%rxyz(3,iat)
           runObj%atoms%astruct%rxyz(1,iat)=runObj%atoms%astruct%rxyz(1,iat)+beta*outs%fxyz(1,iat)
           runObj%atoms%astruct%rxyz(2,iat)=runObj%atoms%astruct%rxyz(2,iat)+beta*outs%fxyz(2,iat)
           runObj%atoms%astruct%rxyz(3,iat)=runObj%atoms%astruct%rxyz(3,iat)+beta*outs%fxyz(3,iat)
           outsold%fxyz(1,iat)=outs%fxyz(1,iat)
           outsold%fxyz(2,iat)=outs%fxyz(2,iat)
           outsold%fxyz(3,iat)=outs%fxyz(3,iat)
        enddo
        !call atomic_axpy(at,wpos,beta,ff,wpos)
        outsold%energy=outs%energy
        fnrmold=fnrm
        !n(c) betalastold=betalast
     endif
  
     if (iproc == 0.and.parmin%verbosity > 0) then
        write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
        ncount_bigdft,itsd,"GEOPT_VSSD",outs%energy,outs%energy-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,& 
        "beta=",beta,"last beta=",betalast
     end if

     if (iproc == 0.and.parmin%verbosity > 0) then
        call yaml_open_map('Geometry')
           call yaml_map('Ncount_BigDFT',ncount_bigdft)
           call yaml_map('Iteration',itsd)
           call yaml_map('Geometry Method','GEOPT_VSSD')
           call yaml_map('etot',(/ outs%energy,outs%energy-etotprev /),fmt='(1pe21.14)')
           call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct /), fmt='(1pe10.2)')
           call yaml_map('beta', beta, fmt='(1pe8.2e1)')
           call yaml_map('last beta', betalast, fmt='(1pe8.2e1)')
           call geometry_output(fmax,fnrm,fluct)
        call yaml_close_map()
        !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
        !ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,& 
        !&"beta=",beta,"last beta=",betalast
    end if

    etotprev=outs%energy
!     if (iproc == 0) write(16,'(i5,1x,e12.5,1x,e21.14,a,e10.3,1x,e10.3)') itsd,sqrt(fnrm),etot,' GEOPT VSSD ',beta,betalast
     if(iproc==0)call timeleft(tt)
     call MPI_BCAST(tt,1,MPI_DOUBLE_PRECISION,0,bigdft_mpi%mpi_comm,i_stat)
     if(tt<0) exit loop_ntsd


  enddo loop_ntsd

  if (iproc == 0.and.parmin%verbosity > 0) & 
       write(16,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
       ncount_bigdft,itsd,"GEOPT_VSSD",outs%energy,outs%energy-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,& 
       "beta=",beta,"last beta=",betalast

  if (iproc == 0.and.parmin%verbosity > 0) then
     call yaml_open_map('Geometry')
        call yaml_map('Ncount_BigDFT',ncount_bigdft)
        call yaml_map('Iteration',itsd)
        call yaml_map('Geometry Method','GEOPT_VSSD')
        call yaml_map('etot',(/ outs%energy,outs%energy-etotprev /),fmt='(1pe21.14)')
        call yaml_map('Forces', (/ fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct /), fmt='(1pe10.2)')
        call yaml_map('beta', beta, fmt='(1pe8.2e1)')
        call yaml_map('last beta', betalast, fmt='(1pe8.2e1)')
        call geometry_output(fmax,fnrm,fluct)
     call yaml_close_map()
     !write(* ,'(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,1pe8.2E1,2x,a,1pe8.2E1)') &
     !&ncount_bigdft,itsd,"GEOPT_VSSD",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*runObj%inputs%frac_fluct,fluct,& 
     !&"beta=",beta,"last beta=",betalast
  end if

  if (iproc == 0 .and. itsd == nitsd+1) &
       write(16,'(a,i5,e9.2,e18.10,e9.2)') '---- SD FAILED  TO CONVERGE'
  if (iproc == 0) then
     if (parmin%verbosity > 0) then
        write(16,'(a,i5,e9.2,e18.10)') 'variable stepsize SD FINISHED,iter, force norm,energy',itsd,sqrt(fnrm),outs%energy
        write(16,'(a,e9.2)') 'suggested value for stepsize:', betaxx
     end if
     write(fn4,'(i4.4)') ncount_bigdft-1
     write(comment,'(a,1pe10.3)')'VSSD:fnrm= ',sqrt(fnrm)
     call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
          & outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms,trim(comment),&
          forces=outs%fxyz)
  endif


  call f_free(posold)
  call deallocate_global_output(outsold)

END SUBROUTINE vstepsd
