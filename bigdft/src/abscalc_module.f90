!> @file
!!  Type and routines of interest in XANE calculation
!! @author
!!    Copyright (C) 2009-2014 BigDFT group (AM, LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module module_abscalc
  use module_base
  use module_types
  use psp_projectors_base, only: free_DFT_PSP_projectors
  use gaussians, only: deallocate_gwf
  implicit none

  !> Contains all array necessary to apply preconditioning projectors 
  type, public :: pcproj_data_type
     type(DFT_PSP_projectors) :: pc_nl
!     real(gp), pointer :: pc_proj(:)
     integer, dimension(:), pointer :: ilr_to_mproj, iproj_to_l
     real(gp), dimension(:), pointer ::  iproj_to_ene
     real(gp), dimension(:), pointer ::  iproj_to_factor
     integer, dimension(:), pointer :: iorbtolr
     integer :: mprojtot
     type(gaussian_basis)  :: G          
     real(gp), dimension(:), pointer :: gaenes
     real(gp) :: ecut_pc
     logical :: DistProjApply
  end type pcproj_data_type

  !> Contains all array necessary to apply preconditioning projectors 
  type, public :: PAWproj_data_type
     type(DFT_PSP_projectors) :: paw_nl
     integer , pointer , dimension(:) :: iproj_to_paw_nchannels
     integer , pointer , dimension(:) :: ilr_to_mproj, iproj_to_l
     integer , pointer , dimension(:) :: iprojto_imatrixbeg

     integer :: mprojtot
!     real(gp), pointer :: paw_proj(:)
     type(gaussian_basis_c)  :: G          
     integer, pointer :: iorbtolr(:)
     logical :: DistProjApply
  end type PAWproj_data_type

  interface
!!$     subroutine createPcProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs, &
!!$          &   radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
!!$          pcproj_data , Glr)
!!$
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       integer, intent(in) :: iproc,n1,n2,n3
!!$       real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
!!$       type(atoms_data), intent(in) :: at
!!$       type(orbitals_data), intent(in) :: orbs
!!$
!!$       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!!$       real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
!!$       real(gp), intent(in):: ecut_pc
!!$
!!$       type(pcproj_data_type) ::pcproj_data
!!$
!!$       type(locreg_descriptors),  intent(in):: Glr
!!$
!!$     END SUBROUTINE createPcProjectorsArrays


!!$     subroutine applyPCprojectors(orbs,at,&
!!$          &   hx,hy,hz,Glr,PPD,psi,hpsi, dotest)
!!$
!!$       use module_base
!!$       use module_types
!!$       type(orbitals_data), intent(inout) :: orbs
!!$       type(atoms_data) :: at
!!$       real(gp), intent(in) :: hx,hy,hz
!!$       type(locreg_descriptors), intent(in) :: Glr
!!$       type(pcproj_data_type) ::PPD
!!$       real(wp), dimension(:), pointer :: psi, hpsi
!!$       logical, optional :: dotest
!!$     END SUBROUTINE applyPCprojectors
!!$
!!$
!!$     subroutine applyPAWprojectors(orbs,at,&
!!$          &   hx,hy,hz,Glr,PAWD,psi,hpsi,  paw_matrix, dosuperposition , &
!!$          sup_iatom, sup_l, sup_arraym) !, sup_arraychannel)
!!$
!!$       use module_base
!!$       use module_types
!!$
!!$       type(orbitals_data), intent(inout) :: orbs
!!$       type(atoms_data) :: at
!!$       real(gp), intent(in) :: hx,hy,hz
!!$       type(locreg_descriptors), intent(in) :: Glr
!!$       type(pawproj_data_type) ::PAWD
!!$       real(wp), dimension(:), pointer :: psi, hpsi, paw_matrix
!!$       logical dosuperposition
!!$       integer, optional :: sup_iatom, sup_l
!!$       real(wp) , dimension(:), pointer, optional :: sup_arraym !, sup_arraychannel
!!$
!!$     END SUBROUTINE applyPAWprojectors
     subroutine gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
          &   zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
          aeval,ng,psi,res,chrg,&
          &   Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw, PAWpatch, &
          psipsigrid)
       use module_base, only: gp

       implicit real(gp) (a-h,o-z)
       logical :: noproj, readytoexit
       integer, parameter :: n_int=1000
       dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
            &   hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
            gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
            &   pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
            potgrd(n_int),&
            &   rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
            occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
            &   vh(0:ng,0:ng,4,0:ng,0:ng,4),&
            res(noccmax,lmax+1),xp(0:ng),& 
            psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol),&
            &   psigrid_naked_2(Ngrid,Nsol), projgrid(Ngrid,3), &
            rhogrid(Ngrid), potgrid(Ngrid), psigrid_not_fitted(Ngrid,Nsol),&
            &   psigrid_not_fitted_2(Ngrid,Nsol),&
            vxcgrid(Ngrid), &
            &   Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
            H(Nsol, Nsol), &
            &   H_2(Nsol, Nsol), &
            Hcorrected(Nsol, Nsol), &
            &   Hadd(Nsol, Nsol), Egrid_tmp(Nsol),Egrid_tmp_2(Nsol), Etofit(Nsol), &
            Soverlap(Nsol,Nsol), Tpsigrid(Nsol,Ngrid ),Tpsigrid_dum(Nsol, Ngrid),valuesatp(Nsol), &
            &   PAWpatch(Npaw, Npaw ), Spsitildes(Npaw, Npaw), genS(Nsol,Nsol), genH(Nsol,Nsol) , dumH(Nsol,Nsol)

       real(gp) , optional :: psipsigrid(Ngrid, Nsol)


       real(gp) :: rgrid(Ngrid), ene_m, ene_p, factadd, rcond, fixfact
       real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid), dumgrid3(Ngrid)
       logical dofit
       integer real_start, iocc, iwork(Nsol), INFO, volta, ngrid_box_2
       character(1) EQUED
       integer ipiv(Nsol), Npaw
     END SUBROUTINE gatom_modified

     subroutine abs_generator_modified(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,&
          &   psi, aeval, occup, psp_modifier, &
          Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw,  PAWpatch , psipsigrid )

       use module_base, only: gp
       implicit none
       integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,noccmax, lmax, Nsol, labs, Ngrid,  Ngrid_box
       real(gp), dimension(0:4,0:6), intent(in) :: psppar
       !! real(gp), dimension(:,:), intent(in) :: psppar
       integer, intent(in) :: psp_modifier, Npaw

       real(gp), dimension(ng+1), intent(out) :: expo

       integer, parameter :: n_int=1000

       real(gp), dimension(0:ng,noccmax,lmax+1), intent(out) :: psi, Egrid(Nsol),&
            &   rgrid(Ngrid), psigrid(Ngrid,Nsol  )
       real(gp),   intent(out), optional  :: psipsigrid(Ngrid,Nsol  )
       real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
       real(gp):: PAWpatch(Npaw,Npaw)

       !local variables
     END SUBROUTINE abs_generator_modified

!!$     subroutine cg_spectra(iproc,nproc,at,hx,hy,hz,rxyz,&
!!$          radii_cf,nlpsp,lr,ngatherarr,ndimpot,potential,&
!!$          energs,nspin,GPU,in_iat_absorber,in , PAWD  )! aggiunger a interface
!!$       !n(c) use module_base
!!$       use module_types
!!$       implicit none
!!$       integer  :: iproc,nproc,ndimpot,nspin
!!$       real(gp)  :: hx,hy,hz
!!$       type(atoms_data), target :: at
!!$       type(DFT_PSP_projectors), target :: nlpsp
!!$       type(locreg_descriptors), target :: lr
!!$       integer, dimension(0:nproc-1,2), target :: ngatherarr 
!!$       real(gp), dimension(3,at%astruct%nat), target :: rxyz
!!$       real(gp), dimension(at%astruct%ntypes,3), intent(in), target ::  radii_cf
!!$       real(wp), dimension(max(ndimpot,1),nspin), target :: potential
!!$       type(energy_terms), intent(inout) :: energs
!!$       type(GPU_pointers), intent(inout) , target :: GPU
!!$       integer, intent(in) :: in_iat_absorber
!!$       type(pawproj_data_type), target ::PAWD
!!$
!!$       type(input_variables),intent(in) :: in
!!$
!!$     END SUBROUTINE cg_spectra

  end interface


contains
  

  subroutine deallocate_pawproj_data(pawproj_data)
    use module_base
    implicit none
    type(pawproj_data_type), intent(inout) :: pawproj_data

    if(associated(pawproj_data%paw_nl%proj)) then

       call f_free_ptr(pawproj_data% ilr_to_mproj)
       call f_free_ptr(pawproj_data%  iproj_to_l)
       call f_free_ptr(pawproj_data%  iproj_to_paw_nchannels)
       call f_free_ptr(pawproj_data%  iprojto_imatrixbeg)
       call f_free_ptr(pawproj_data%  iorbtolr)

       call free_DFT_PSP_projectors(pawproj_data%paw_nl)

       if(pawproj_data%DistProjApply) then
          call deallocate_gwf_c(pawproj_data%G)
       endif
    end if

  END SUBROUTINE deallocate_pawproj_data


  !> deallocate_pcproj_data
  subroutine deallocate_pcproj_data(pcproj_data)
    use module_base
    implicit none
    type(pcproj_data_type), intent(inout) :: pcproj_data
    
    if(associated(pcproj_data%pc_nl%proj)) then
       call f_free_ptr( pcproj_data% ilr_to_mproj)
       call f_free_ptr(  pcproj_data% iproj_to_ene)
       call f_free_ptr(  pcproj_data% iproj_to_factor)
       call f_free_ptr(pcproj_data%  iproj_to_l)
       call f_free_ptr(pcproj_data%  iorbtolr)
       call f_free_ptr(pcproj_data%  gaenes)
       call free_DFT_PSP_projectors(pcproj_data%pc_nl)

       if(pcproj_data%DistProjApply) then
          call deallocate_gwf(pcproj_data%G)
       endif


    end if

  END SUBROUTINE deallocate_pcproj_data


  subroutine applyPAWprojectors(orbs,at,&
       &   hx,hy,hz,Glr,PAWD,psi,hpsi,  paw_matrix, dosuperposition , &
       &   sup_iatom, sup_l, sup_arraym)
    use module_base
    use module_types
    implicit none

    type(orbitals_data), intent(inout) :: orbs
    type(atoms_data) :: at
    real(gp), intent(in) :: hx,hy,hz
    type(locreg_descriptors), intent(in) :: Glr
    type(pawproj_data_type) ::PAWD
    real(wp), dimension(:), pointer :: psi, hpsi, paw_matrix
    logical ::  dosuperposition
    integer , optional :: sup_iatom, sup_l
    real(wp) , dimension(:), pointer, optional :: sup_arraym 
    ! local variables
    character(len=*), parameter :: subname='applyPAWprojectors'
    integer :: ikpt, istart_ck, ispsi_k, isorb,ieorb, ispsi, iproj, istart_c,&
         &   mproj, mdone, ispinor, istart_c_i, mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, &
         &   jseg_c, iproj_old, iorb, ncplx, l, jorb, lsign, ncplx_global
    real(gp) :: eproj_spinor
    real(gp) :: kx,ky,kz

    integer , parameter :: dotbuffersize = 1000
    real(dp)  :: dotbuffer(dotbuffersize), dotbufferbis(dotbuffersize)
    integer :: ibuffer, ichannel, nchannels, imatrix
    logical :: lfound_sup
    integer :: iat, old_istart_c, iatat , m, nspinor

    if (orbs%norbp.gt.0) then


       !apply the projectors  k-point of the processor
       !starting k-point
       ikpt=orbs%iokpt(1)
       istart_ck=1
       ispsi_k=1
       imatrix=1

!!$ check that the coarse wavelets cover the whole box
       if(Glr%wfd%nvctr_c .ne. &
            ((Glr%d%n1+1)*(Glr%d%n2+1)*(Glr%d%n3+1))) then
          print *," WARNING : coarse wavelets dont cover the whole box "
       endif

       if(dosuperposition) then 
          lfound_sup=.false.
       endif
       ncplx_global=min(orbs%nspinor,2)

       loop_kpt: do

          call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

          ! loop over all my orbitals
          do iorb=isorb,ieorb
             istart_c=istart_ck
             iproj=1
             iat=0

             do iatat=1, at%astruct%nat
                if (  at%paw_NofL(at%astruct%iatype(iatat)).gt.0  ) then
                   iat=iat+1
                   istart_c_i=istart_c
                   iproj_old=iproj
                   ispsi=ispsi_k
!!!! do iorb=isorb,ieorb


                   mproj= PAWD%ilr_to_mproj(iat)

                   if( ikpt .ne. orbs%iokpt(iorb) ) then
                      STOP " ikpt .ne. orbs%iokpt(iorb) in applypawprojectors " 
                   end if
                   kx=orbs%kpts(1,ikpt)
                   ky=orbs%kpts(2,ikpt)
                   kz=orbs%kpts(3,ikpt)
                   call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)

                   do ispinor=1,orbs%nspinor,ncplx_global
                      eproj_spinor=0.0_gp
                      if (ispinor >= 2) istart_c=istart_c_i
                      call plr_segs_and_vctrs(PAWD%paw_nl%pspd(iat)%plr,&
                           mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
                      jseg_c=1

                      mdone=0
                      iproj=iproj_old
                      if(mproj>0) then
                         if(  PAWD%DistProjApply) then
                            jorb=1
                            do while(jorb<=PAWD%G%ncoeff .and. PAWD%iorbtolr(jorb)/= iat) 
                               jorb=jorb+1
                            end do
                            if(jorb<PAWD%G%ncoeff) then
                               call fillPawProjOnTheFly(PAWD,Glr,iat,hx,hy,hz,&
                                    &   kx,ky,kz,&
                                    &   jorb,istart_c,at%astruct%geocode,at,iatat ) 
                            endif
                         end if
                      endif

                      do while(mdone< mproj)
                         lsign = PAWD%iproj_to_l(iproj)
                         l=abs(lsign)
                         ibuffer = 0
                         nchannels =  PAWD% iproj_to_paw_nchannels(iproj)
                         imatrix=PAWD%iprojto_imatrixbeg(iproj)
!!$
!!$                    if(.not. dosuperposition) then
!!$                       print *, "applying paw for l= ", l,&
!!$                            "  primo elemento ", paw_matrix(PAWD%iprojto_imatrixbeg(iproj))
!!$                    end if
                         old_istart_c=istart_c
                         do ichannel=1, nchannels
                            do m=1,2*l-1
                               ibuffer=ibuffer+1
                               if(ibuffer.gt.dotbuffersize ) then
                                  STOP 'ibuffer.gt.dotbuffersize'
                               end if

                               if( .not. dosuperposition .and. lsign>0 ) then
                                  call wpdot_wrap(ncplx,  &
                                       Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,&
                                       Glr%wfd%nseg_c,Glr%wfd%nseg_f,&
                                       Glr%wfd%keyvglob(1),Glr%wfd%keyglob(1,1),&
                                       psi(ispsi+&
                                       (ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)),&
                                       mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                                       PAWD%paw_nl%pspd(iat)%plr%wfd%keyvglob(jseg_c),&
                                       PAWD%paw_nl%pspd(iat)%plr%wfd%keyglob(1,jseg_c),&
                                       PAWD%paw_nl%proj(istart_c),&
                                       dotbuffer( ibuffer ) )
                               end if
                               ibuffer=ibuffer + (ncplx-1)

!!$                          !! TTTTTTTTTTTTTTTTTTTt TEST TTTTTTTTTTTTTTTTTTT
!!$                          call wpdot_wrap(ncplx,  &
!!$                               mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PAWD%paw_nlpspd%keyv_p(jseg_c),&
!!$                               PAWD%paw_nlpspd%keyg_p(1,jseg_c),PAWD%paw_proj(istart_c),& 
!!$                               mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PAWD%paw_nlpspd%keyv_p(jseg_c),&
!!$                               PAWD%paw_nlpspd%keyg_p(1,jseg_c),&
!!$                               PAWD%paw_proj(istart_c),&
!!$                               eproj_spinor)
!!$                          print *, "TEST:  THE PROJECTOR ichannel = ", ichannel, " m=",m, " HAS SQUARED MODULUS  " ,eproj_spinor 

!!$                          !! plot -------------------------------------------------------------
!!$                          Plr%d%n1 = Glr%d%n1
!!$                          Plr%d%n2 = Glr%d%n2
!!$                          Plr%d%n3 = Glr%d%n3
!!$                          Plr%geocode = at%astruct%geocode
!!$                          Plr%wfd%nvctr_c  =PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
!!$                          Plr%wfd%nvctr_f  =PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
!!$                          Plr%wfd%nseg_c   =PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)
!!$                          Plr%wfd%nseg_f   =PAWD%paw_nlpspd%nseg_p(2*iat  )-PAWD%paw_nlpspd%nseg_p(2*iat-1)
!!$                          call allocate_wfd(Plr%wfd,subname)
!!$                          Plr%wfd%keyv(:)  = PAWD%paw_nlpspd%keyv_p( PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:&
!!$                               PAWD%paw_nlpspd%nseg_p(2*iat)   )
!!$                          Plr%wfd%keyg(1:2, :)  = PAWD%paw_nlpspd%keyg_p( 1:2,  PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:&
!!$                               PAWD%paw_nlpspd%nseg_p(2*iat)   )
!!$
!!$                          !! ---------------  use this to plot projectors
!!$                          write(orbname,'(A,i4.4)')'paw_',iproj
!!$                          Plr%bounds = Glr%bounds
!!$                          Plr%d          = Glr%d
!!$                          call plot_wf_cube(orbname,at,Plr,hx,hy,hz,rxyz, PAWD%paw_proj(istart_c) ,"1234567890" ) 
!!$                          !! END plot ----------------------------------------------------------


                               !! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

                               istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                               iproj=iproj+1
                               mdone=mdone+1
                            end do
                         end do

!!$                    call DGEMM('N','N', nchannels ,(2*l-1)  , nchannels  ,&
!!$                         1.0d0 , paw_matrix(imatrix) , nchannels ,&
!!$                         dotbuffer  ,(2*l-1), 0.0D0 , dotbufferbis  ,(2*l-1))


                         if( .not. dosuperposition) then
                            if(lsign>0) then
                               call DGEMM('N','N',(2*l-1)*ncplx  , nchannels , nchannels  ,&
                                    &   1.0d0 ,dotbuffer , (2*l-1)*ncplx ,&
                                    &   paw_matrix(imatrix)  ,nchannels , 0.0D0 , dotbufferbis  ,(2*l-1)*ncplx )
                            else
                               dotbufferbis=0.0_wp
                            endif
                         else
                            !print *,'here',nchannels
                            if( sup_iatom .eq. iatat .and. (-sup_l) .eq. lsign ) then
                               do ichannel=1, nchannels
                                  do m=1,2*l-1
                                     dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx           ) = 0.0_gp ! keep this before
                                     dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx -(ncplx-1)) = sup_arraym(m)
                                  end do
                               enddo
                               lfound_sup=.true.
                            else
                               do ichannel=1, nchannels
                                  do m=1,2*l-1
                                     dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx           ) = 0.0_gp 
                                     dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx -(ncplx-1)) = 0.0_gp
                                  end do
                               enddo
                            endif
                         endif


                         ibuffer=0
                         iproj    =  iproj  - nchannels * ( 2*l-1 )
                         istart_c = old_istart_c

                         do ichannel=1, nchannels
                            do m=1,2*l-1
                               ibuffer=ibuffer+1

                               call waxpy_wrap(ncplx,dotbufferbis( ibuffer ) ,&
                                    mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                                 &   PAWD%paw_nlpspd%keyv_p(jseg_c),PAWD%paw_nlpspd%keyg_p(1,jseg_c),&
                                    PAWD%paw_nl%pspd(iat)%plr%wfd%keyvglob(jseg_c),&
                                    PAWD%paw_nl%pspd(iat)%plr%wfd%keyglob(1,jseg_c),&
                                    PAWD%paw_nl%proj(istart_c),&
                                    Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,Glr%wfd%nseg_c,Glr%wfd%nseg_f,&
                                    Glr%wfd%keyvglob(1),Glr%wfd%keyglob(1,1),&
                                    hpsi(ispsi+(ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)  )&
                                    )


                               istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                               iproj=iproj+1
                               ibuffer=ibuffer + (ncplx-1)
                            end do
                         end do
                      end do

                      mdone=0
!!$ iproj=iproj_old
                      istart_c=istart_c_i
                      ispsi=ispsi+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*nspinor

                   end do

                   if( PAWD%DistProjApply ) then
                      istart_c=1
                   else
                      istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*mproj*ncplx
                   endif

                end if
             end do

             ispsi_k=ispsi
          end do
          istart_ck=istart_c

          if(  dosuperposition ) then
             if(.not. lfound_sup) then
                print *, " initial state not found in routine ",subname
                STOP 
             endif
          endif


          if (ieorb == orbs%norbp) exit loop_kpt
          ikpt=ikpt+1


       end do loop_kpt
    end if
  END SUBROUTINE applyPAWprojectors


  subroutine applyPCprojectors(orbs,at,&
       &   hx,hy,hz,Glr,PPD,psi,hpsi, dotest)

    use module_base
    use module_types

    implicit none

    type(orbitals_data), intent(inout) :: orbs
    type(atoms_data) :: at
    real(gp), intent(in) :: hx,hy,hz
    type(locreg_descriptors), intent(in) :: Glr
    type(pcproj_data_type) ::PPD
    real(wp), dimension(:), pointer :: psi, hpsi
    logical, optional :: dotest


    ! local variables
    character(len=*), parameter :: subname='applyPCprojectors'
    character(len=11) :: orbname
    type(locreg_descriptors) :: Plr
    integer :: ikpt, istart_ck, ispsi_k, isorb,ieorb, ispsi, iproj, istart_c,&
         &   mproj, mdone, ispinor, istart_c_i, mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, &
         &   jseg_c, iproj_old, iorb, ncplx, l, i, jorb
    real(gp) eproj_spinor,psppar_aux(0:4, 0:6)
    integer :: iat , nspinor

    !apply the projectors  k-point of the processor
    !starting k-point
    ikpt=orbs%iokpt(1)
    istart_ck=1
    ispsi_k=1
    loop_kpt: do

       call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

       ! loop over all my orbitals
       istart_c=1
       iproj=1

       do iat=1,at%astruct%nat
          istart_c_i=istart_c
          iproj_old=iproj
          ispsi=ispsi_k
          do iorb=isorb,ieorb

             mproj= PPD%ilr_to_mproj(iat)

             call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)


             do ispinor=1,orbs%nspinor,ncplx
                eproj_spinor=0.0_gp

                if (ispinor >= 2) istart_c=istart_c_i

                call plr_segs_and_vctrs(PPD%pc_nl%pspd(iat)%plr,&
                     mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
                jseg_c=1

                mdone=0
                iproj=iproj_old

                if(mproj>0) then
                   if(  PPD%DistProjApply) then
                      jorb=1
                      do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)/= iat) 
                         jorb=jorb+1
                      end do
                      if(jorb<PPD%G%ncoeff) then
                         call fillPcProjOnTheFly(PPD,Glr,iat,at,hx,hy,hz,&
                              jorb,PPD%ecut_pc,istart_c) 
                      endif
                   end if
                endif


                do while(mdone< mproj)

                   l = PPD%iproj_to_l(iproj)

                   i=1
                   psppar_aux=0.0_gp
                   !! psppar_aux(l,i)=1.0_gp/PPD%iproj_to_ene(iproj)
                   !! psppar_aux(l,i)=1.0_gp  ! *iorb
                   psppar_aux(l,i)=PPD%iproj_to_factor(iproj)  


                   call applyprojector(ncplx,l,i, psppar_aux(0,0), 2 ,&
                        Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c,&
                        Glr%wfd%nseg_f,&
                        Glr%wfd%keyvglob(1),Glr%wfd%keyglob(1,1),&
                        mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                        PPD%pc_nl%pspd(iat)%plr%wfd%keyvglob(jseg_c),&
                        PPD%pc_nl%pspd(iat)%plr%wfd%keyglob(1,jseg_c),&
!!$                       PPD%pc_nlpspd%keyv_p(jseg_c),PPD%pc_nlpspd%keyg_p(1,jseg_c),&
                        PPD%pc_nl%proj(istart_c),&
                        psi(ispsi+ (ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)  ),&
                        hpsi(ispsi+(ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)  ),&
                        eproj_spinor)


                   if(iorb==1) then         
                      if( present(dotest) ) then
                         eproj_spinor=0.0_gp
!!$                        call wpdot_wrap(ncplx,  &
!!$                             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                             PPD%pc_nlpspd%keyv_p(jseg_c),&
!!$                             PPD%pc_nlpspd%keyg_p(1,jseg_c),&
!!$                             PPD%pc_proj(istart_c),& 
!!$                             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                             PPD%pc_nlpspd%keyv_p(jseg_c),&
!!$                             PPD%pc_nlpspd%keyg_p(1,jseg_c),&
!!$                             PPD%pc_proj(istart_c),&
!!$                             eproj_spinor)
                         print *, " IL PROIETTORE HA MODULO QUADRO  " ,eproj_spinor 
                         if(dotest) then
                            !! ---------------  use this to plot projectors
                            write(orbname,'(A,i4.4)')'pc_',iproj                     
                            Plr%d%n1=Glr%d%n1
                            Plr%d%n2=Glr%d%n2
                            Plr%d%n3=Glr%d%n3
                            Plr%geocode=at%astruct%geocode

                            call plr_segs_and_vctrs(PPD%pc_nl%pspd(iat)%plr,&
                                 Plr%wfd%nseg_c,Plr%wfd%nseg_f,&
                                 Plr%wfd%nvctr_c,Plr%wfd%nvctr_f)                  
!!$                           Plr%wfd%nvctr_c  =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!!$                           Plr%wfd%nvctr_f  =PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!!$                           Plr%wfd%nseg_c   =PPD%pc_nlpspd%nseg_p(2*iat-1 )-PPD%pc_nlpspd%nseg_p(2*iat-2)
!!$                           Plr%wfd%nseg_f   =PPD%pc_nlpspd%nseg_p(2*iat  ) -PPD%pc_nlpspd%nseg_p(2*iat-1)
                           ! call allocate_wfd(Plr%wfd)
!!$                           Plr%wfd%keyv(:)=PPD%pc_nlpspd%keyv_p(PPD%pc_nlpspd%nseg_p(2*iat-2)+1:&
!!$                              &   PPD%pc_nlpspd%nseg_p(2*iat)   )
!!$                           Plr%wfd%keyg(1:2, :)  = PPD%pc_nlpspd%keyg_p( 1:2,  PPD%pc_nlpspd%nseg_p(2*iat-2)+1:&
!!$                              &   PPD%pc_nlpspd%nseg_p(2*iat)   )
                            Plr%bounds = Glr%bounds
                            Plr%d          = Glr%d                    
                            !! call plot_wf_cube(orbname,at,Plr,hx,hy,hz,rxyz, PPD%pc_proj(istart_c) ,"1234567890" ) 
                            !call deallocate_wfd(Plr%wfd)
                         endif
                      endif

                   endif
                   istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx
                   iproj=iproj+(2*l-1)
                   mdone=mdone+(2*l-1)
                end do
             end  do
             istart_c=istart_c_i

             if( present(dotest) ) then
                if(dotest) then
                   eproj_spinor=0.0_gp
!!$                  call wpdot_wrap(ncplx,  &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), hpsi(ispsi + 0 ), &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), hpsi(ispsi + 0 ),  &
!!$                     &   eproj_spinor)
                   print *, "hpsi  HA MODULO QUADRO  " ,eproj_spinor 
                   eproj_spinor=0.0_gp
!!$                  call wpdot_wrap(ncplx,  &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), psi(ispsi + 0 ), &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), psi(ispsi + 0 ),  &
!!$                     &   eproj_spinor)
                   print *, "psi  HA MODULO QUADRO  " ,eproj_spinor         
                   !! CECCARE IPROJ = mproj tot, istart_c=nelproj 
                   write(orbname,'(A,i4.4)')'pcorb_',iorb
                   !! call plot_wf_cube(orbname,at,Glr,hx,hy,hz,rxyz,hpsi(ispsi + 0 ),"dopoprec.." ) ! solo spinore 1
                end if
             endif

             ispsi=ispsi+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*nspinor

          end do


          if( PPD%DistProjApply ) then
             istart_c=1
          else
             istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*mproj
          endif

       end do

       !! istart_ck=istart_c  non si incrementa

       if (ieorb == orbs%norbp) exit loop_kpt
       ikpt=ikpt+1
       ispsi_k=ispsi
    end do loop_kpt

  END SUBROUTINE applyPCprojectors

end module module_abscalc
