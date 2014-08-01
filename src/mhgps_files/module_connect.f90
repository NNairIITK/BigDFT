!> @file
!! module implementing the connection algorithm(s)
!!     
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_connect

contains
!=====================================================================
subroutine connect_recurs(nat,alat,rcov,,nbond,iconnect,rxyz1,rxyz2)
    use module_base
    use module_ls_rmsd
    implicit none
    !parameters
    integer, intent(in) :: nat
    integer, intent(in) :: nbond
    real(gp), intent(in) :: rcov(nat)
    real(gp), intent(in) :: alat(3)
    integer, intent(in)       :: iconnect(2,nbond)
    real(gp), intent(inout) :: rxyz1(3,nat), rxyz2(3,nat)
    !internal
    real(gp) :: saddle(3,nat)
    real(gp) :: rxyzleft(3,nat)
    real(gp) :: rxyzright(3,nat)
    real(gp) :: fxyz(3,nat)
    real(gp) :: rotforce(3,nat)
    real(gp) :: minmode(3,nat)
    real(gp) :: etot,displ=0._gp,ener_count=0._gp
    logical :: converged =.false.

    !rmsd alignment
    call superimpose(nat,rxyz1,rxyz2)

    !get input guess
    call get_ts_guess(nat,alat,rxyz1,rxyz2,saddle,minmode)

    !compute saddle
    call findsad(nat,alat,rcov,nbond,iconnect,&
                saddle,etot,fxyz,minmode,displ,ener_count,&
                  rotforce,converged)

    !pushoff and minimize left and right
    call pushoff(nat,saddle,minmode,leftmin,rightmin)

    !check if relaxed structures are identical with saddle itself

    !is minimum, obtained by relaxation from left bar end identical to left
    !input minimum?
!    lnl=equal(nid,en_delta,fp_delta,epot1,leftminener(nsad_local),fp1,leftminfp(1,nsad_local),'min')

    !is minimum obtained by relaxation from right bar end identical to right
    !input minimum?
!    rnr=equal(nid,en_delta,fp_delta,epot2,rightminener(nsad_local),fp2,rightminfp(1,nsad_local),'min')

    !is minimum obtained by relaxation from left bar end identical to right
    !input minimum?
!    lnr=equal(nid,en_delta,fp_delta,epot2,leftminener(nsad_local),fp2,leftminfp(1,nsad_local),'min')

    !is minimum obtained by relaxation from right bar end identical to left
    !input minimum?
!    rnl=equal(nid,en_delta,fp_delta,epot1,rightminener(nsad_local),fp1,rightminfp(1,nsad_local),'min')

!    if((lnl .and. rnr) .or. (lnr .and. rnl))then!connection done
!        return
!    endif

!    if(lnl .and. (.not. rnr))then!connect right input min with right relaxed barend
!        call
!connect(nat,alat,nid,bpPreset,rcov,xat,en_delta,fp_delta,en_delta_sp,fp_delta_sp,st,fmaxtol,fnrmtol,alphax,fire_dt_max,nsadmax&
!            &,saddle,leftmin,rightmin,saddleener,leftminener,rightminener,saddlefp,leftminfp,rightminfp&
!            &,nsad,rightmin(1,nsad_local),rxyz2,rightminener(nsad_local),epot2,rightminfp(1,nsad_local),fp2,connected)
!        return
!    endif
!
!    if(rnr .and. (.not. lnl))then!connect left relaxed bar end with left input min
!        call
!connect(nat,alat,nid,bpPreset,rcov,xat,en_delta,fp_delta,en_delta_sp,fp_delta_sp,st,fmaxtol,fnrmtol,alphax,fire_dt_max,nsadmax&
!            &,saddle,leftmin,rightmin,saddleener,leftminener,rightminener,saddlefp,leftminfp,rightminfp&
!            &,nsad,rxyz1,leftmin(1,nsad_local),epot1,leftminener(nsad_local),fp1,leftminfp(1,nsad_local),connected)
!        return
!    endif
!
!    if(lnr .and. .not. rnl)then!connect right relaxed bar end with left input min
!        call
!connect(nat,alat,nid,bpPreset,rcov,xat,en_delta,fp_delta,en_delta_sp,fp_delta_sp,st,fmaxtol,fnrmtol,alphax,fire_dt_max,nsadmax&
!            &,saddle,leftmin,rightmin,saddleener,leftminener,rightminener,saddlefp,leftminfp,rightminfp&
!            &,nsad,rxyz1,rightmin(1,nsad_local),epot1,rightminener(nsad_local),fp1,rightminfp(1,nsad_local),connected)
!        return
!    endif
!
!    if(.not. lnr .and. rnl)then!connect left relaxed bar end with right input min
!        call
!connect(nat,alat,nid,bpPreset,rcov,xat,en_delta,fp_delta,en_delta_sp,fp_delta_sp,st,fmaxtol,fnrmtol,alphax,fire_dt_max,nsadmax&
!            &,saddle,leftmin,rightmin,saddleener,leftminener,rightminener,saddlefp,leftminfp,rightminfp&
!            &,nsad,rxyz2,leftmin(1,nsad_local),epot2,leftminener(nsad_local),fp2,leftminfp(1,nsad_local),connected)
!        return
!    endif
!
!
!    if((.not. lnl) .and. (.not. rnr))then!connect left input min with left relaxed bar end  and right input min with right relaxed bar end
!        call
!connect(nat,alat,nid,bpPreset,rcov,xat,en_delta,fp_delta,en_delta_sp,fp_delta_sp,st,fmaxtol,fnrmtol,alphax,fire_dt_max,nsadmax&
!            &,saddle,leftmin,rightmin,saddleener,leftminener,rightminener,saddlefp,leftminfp,rightminfp&
!            &,nsad,rxyz1,leftmin(1,nsad_local),epot1,leftminener(nsad_local),fp1,leftminfp(1,nsad_local),connected)
!        call connect(nat,alat,nid,bpPreset,rcov,xat,en_delta,fp_delta,en_delta_sp,fp_delta_sp,st,fmaxtol,fnrmtol,alphax,fire_dt_max,nsadmax&
!            &,saddle,leftmin,rightmin,saddleener,leftminener,rightminener,saddlefp,leftminfp,rightminfp&
!            &,nsad,rightmin(1,nsad_local),rxyz2,rightminener(nsad_local),epot2,rightminfp(1,nsad_local),fp2,connected)
!        return
!    endif

    !should not happen:
!    write(100,*)'ERROR: none of the checks in connect subroutine matched! STOP'
!    write(100,*)'lnl, lnr, rnr, rnl:'
!    write(100,*)lnl,lnr,rnr,rnl
!    stop


end subroutine
!=====================================================================
subroutine pushoff(nat,saddle,minmode,left,right)
    use module_base
    use module_misc
    use module_global_variables, only: saddle_stepoff
    implicit none
    !parameters 
    integer, intent(in) :: nat
    real(gp), intent(in) :: sad
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(out) :: left(3,nat)
    real(gp), intent(out) :: right(3,nat)
    !internal
    real(gp), intent(out) :: step(3,nat)

    !functions
    real(gp) :: dnrm2
!debug check
if(.not. almostequal(1._gp,dnrm2(3*nat,minmode(1,1),1),4))&
stop'minmode not normalized'

    step = saddle_stepoff*minmode
    left = saddle - step
    right = saddle + step

end subroutine
!=====================================================================


end module
