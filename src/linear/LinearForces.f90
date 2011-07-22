!>  Calculates the nonlocal forces on all atoms arising from the wavefunctions 
!!  belonging to iproc and adds them to the force array
!!   recalculate the projectors at the end if refill flag is .true.
subroutine Linearnonlocal_forces(iproc,Lzd,hx,hy,hz,at,rxyz,&
     orbs,proj,psi,fsep,refill)
  use module_base
  use module_types
  implicit none
  !Arguments-------------
  type(atoms_data), intent(in) :: at
  logical, intent(in) :: refill
  integer, intent(in) :: iproc
  real(gp), intent(in) :: hx,hy,hz
  type(linear_zone_descriptors) :: Lzd
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
  real(wp), dimension(Lzd%Gnlpspd%nprojel), intent(inout) :: proj
  real(gp), dimension(3,at%nat), intent(inout) :: fsep
  !local variables--------------
  character(len=*), parameter :: subname='Linearnonlocal_forces'
  integer :: istart_c,iproj,iat,ityp,i,j,l,m,ilr,ilr2
  integer :: mbseg_c,mbseg_f,jseg_c,jseg_f,iatom
  integer :: mbvctr_c,mbvctr_f,iorb,nwarnings,nspinor,ispinor,jorbd
  real(gp) :: offdiagcoeff,hij,sp0,spi,sp0i,sp0j,spj,orbfac
  integer :: idir,i_all,i_stat,ncplx,icplx,isorb,ikpt,ieorb,istart_ck,ispsi_k,ispsi,jorb
  real(gp), dimension(2,2,3) :: offdiagarr
  real(gp), dimension(:,:), allocatable :: fxyz_orb
  real(dp), dimension(:,:,:,:,:,:,:), allocatable :: scalprod
  
  !quick return if no orbitals on this processor
  if (orbs%norbp == 0) return


  !always put complex scalprod
  !also nspinor for the moment is the biggest as possible
  allocate(scalprod(2,0:3,7,3,4,at%nat,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,scalprod,'scalprod',subname)
  call razero(2*4*7*3*4*at%nat*orbs%norbp*orbs%nspinor,scalprod)

  !calculate the coefficients for the off-diagonal terms
  do l=1,3
     do i=1,2
        do j=i+1,3
           offdiagcoeff=0.0_gp
           if (l==1) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(3._gp/5._gp)
                 if (j==3) offdiagcoeff=0.5_gp*sqrt(5._gp/21._gp)
              else
                 offdiagcoeff=-0.5_gp*sqrt(100._gp/63._gp)
              end if
           else if (l==2) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(5._gp/7._gp)
                 if (j==3) offdiagcoeff=1._gp/6._gp*sqrt(35._gp/11._gp)
              else
                 offdiagcoeff=-7._gp/3._gp*sqrt(1._gp/11._gp)
              end if
           else if (l==3) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(7._gp/9._gp)
                 if (j==3) offdiagcoeff=0.5_gp*sqrt(63._gp/143._gp)
              else
                 offdiagcoeff=-9._gp*sqrt(1._gp/143._gp)
              end if
           end if
           offdiagarr(i,j-i,l)=offdiagcoeff
        end do
     end do
  end do

  if (DistProjApply) then
     !apply the projectors on the fly for each k-point of the processor
     !starting k-point
     ikpt=orbs%iokpt(1)
     ispsi_k=1
     jorb=0
     loop_kptD: do

        call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

        call ncplx_kpt(ikpt,orbs,ncplx)

        !calculate the contribution for each orbital
        !here the nspinor contribution should be adjusted
        ! loop over all my orbitals
        nwarnings=0 !not used, simply initialised 
        ispsi=ispsi_k
        
        do iorb=isorb,ieorb
           ilr = orbs%inWhichLocreg(iorb+orbs%isorb)   !check for the locreg
           do ispinor=1,nspinor,ncplx
              jorb=jorb+1
              iproj = 0 !should be equal to four times nproj at the end
              jseg_c = 1
              iatom = 0
              do iat=1,at%nat
                 if (Lzd%Llr(ilr)%projflg(iat) == 0) cycle
                 iatom = iatom + 1
                 mbseg_c=Lzd%Lnlpspd(ilr)%nseg_p(2*iatom-1)!-Lzd%Lnlpspd(ilr)%nseg_p(2*iat-2)
                 mbseg_f=Lzd%Lnlpspd(ilr)%nseg_p(2*iatom  )!-Lzd%Lnlpspd(ilr)%nseg_p(2*iat-1)
                 jseg_f = jseg_c + Lzd%Lnlpspd(ilr)%nseg_p(2*iatom - 1)
                 mbvctr_c=Lzd%Lnlpspd(ilr)%nvctr_p(2*iatom-1)!-Lzd%Lnlpspd(ilr)%nvctr_p(2*iat-2)
                 mbvctr_f=Lzd%Lnlpspd(ilr)%nvctr_p(2*iatom  )!-Lzd%Lnlpspd(ilr)%nvctr_p(2*iat-1)
                 ityp=at%iatype(iat) 
                 do idir=0,3
                    istart_c = 1
                    !calculate projectors
                    call Linearatom_projector(ikpt,iatom,idir,istart_c,iproj,&
                         Lzd%Llr(ilr),hx,hy,hz,rxyz,at,orbs,Lzd%Lnlpspd(ilr),proj,nwarnings)
                    istart_c=1        !because it is changed in Linearatom_projector
                    do l=1,4
                       do i=1,3
                          if (at%psppar(l,i,ityp) /= 0.0_gp) then
                             do m=1,2*l-1
                                call wpdot_wrap(ncplx,&
                                     Lzd%Llr(ilr)%wfd%nvctr_c,Lzd%Llr(ilr)%wfd%nvctr_f,Lzd%Llr(ilr)%wfd%nseg_c,&
                                     Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%keyv,Lzd%Llr(ilr)%wfd%keyg,psi(ispsi),&
                                     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                                     Lzd%Lnlpspd(ilr)%keyv_p(jseg_c),Lzd%Lnlpspd(ilr)%keyg_p(1,jseg_c),&
                                     proj(istart_c),&
                                     scalprod(1,idir,m,i,l,iat,jorb))
                                istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                             end do
                          end if
                       end do
                    end do
                 end do
                 jseg_c = jseg_c + Lzd%Lnlpspd(ilr)%nseg_p(2*iatom - 1)+ Lzd%Lnlpspd(ilr)%nseg_p(2*iatom)
              end do
              ispsi=ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*ncplx
           end do
        end do

        if (ieorb == orbs%norbp) exit loop_kptD
        ikpt=ikpt+1
        ispsi_k=ispsi
     end do loop_kptD

  else
     !calculate all the scalar products for each direction and each orbitals
     do idir=0,3

! Displaced to inside the loop on iorb because this loop decides the ilr
! Should not be a problem, since for linear scaling, each wavefunction can be in a separate locreg
!        if (idir /= 0) then !for the first run the projectors are already allocated
!           call fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,idir)
!        end if

        !apply the projectors  k-point of the processor
        !starting k-point
        ikpt=orbs%iokpt(1)
        istart_ck=1
        ispsi_k=1
        jorb=0
        ilr2 = 0
        loop_kpt: do

           call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

           call ncplx_kpt(ikpt,orbs,ncplx)
           ! calculate the scalar product for all the orbitals
           ispsi=ispsi_k
           do iorb=isorb,ieorb

              ilr = orbs%inWhichLocreg(iorb+orbs%isorb)   !check for the locreg

              !must recalculate the projectors for each different locreg and each idir
              if (ilr .ne. ilr2) then
                 ilr2 = ilr
                 call Linearfill_projectors(iproc,Lzd%Llr(ilr),hx,hy,hz,at,orbs,rxyz,Lzd%Lnlpspd(ilr),proj,idir)
              end if

              do ispinor=1,nspinor,ncplx
                 jorb=jorb+1
                 ! loop over all projectors of this k-point
                 iproj=0
                 istart_c=istart_ck
                 jseg_c = 1
                 iatom = 0
                 do iat=1,at%nat
                    if (Lzd%Llr(ilr)%projflg(iat) == 0) cycle
                    iatom = iatom + 1
                    mbseg_c=Lzd%Lnlpspd(ilr)%nseg_p(2*iatom-1)!-Lzd%Lnlpspd(ilr)%nseg_p(2*iat-2)
                    mbseg_f=Lzd%Lnlpspd(ilr)%nseg_p(2*iatom  )!-Lzd%Lnlpspd(ilr)%nseg_p(2*iat-1)
                    !jseg_c=Lzd%Lnlpspd(ilr)%nseg_p(2*iat-2)+1
                    !jseg_f=Lzd%Lnlpspd(ilr)%nseg_p(2*iat-1)+1
                    jseg_f = jseg_c + Lzd%Lnlpspd(ilr)%nseg_p(2*iatom - 1)
                    mbvctr_c=Lzd%Lnlpspd(ilr)%nvctr_p(2*iatom-1)!-Lzd%Lnlpspd(ilr)%nvctr_p(2*iat-2)
                    mbvctr_f=Lzd%Lnlpspd(ilr)%nvctr_p(2*iatom  )!-Lzd%Lnlpspd(ilr)%nvctr_p(2*iat-1)
                    ityp=at%iatype(iat)
                    do l=1,4
                       do i=1,3
                          if (at%psppar(l,i,ityp) /= 0.0_gp) then
                             do m=1,2*l-1
                                iproj=iproj+1
                                call wpdot_wrap(ncplx,&
                                     Lzd%Llr(ilr)%wfd%nvctr_c,Lzd%Llr(ilr)%wfd%nvctr_f,Lzd%Llr(ilr)%wfd%nseg_c,&
                                     Lzd%Llr(ilr)%wfd%nseg_f,&
                                     Lzd%Llr(ilr)%wfd%keyv,Lzd%Llr(ilr)%wfd%keyg,psi(ispsi),  &
                                     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                                     Lzd%Lnlpspd(ilr)%keyv_p(jseg_c),Lzd%Lnlpspd(ilr)%keyg_p(1,jseg_c),&
                                     proj(istart_c),scalprod(1,idir,m,i,l,iat,jorb))
                                istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                             end do
                          end if
                       end do
                    end do
                    jseg_c = jseg_c + Lzd%Lnlpspd(ilr)%nseg_p(2*iatom - 1)+ Lzd%Lnlpspd(ilr)%nseg_p(2*iatom)
                 end do
                 ispsi=ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*ncplx
              end do
              if (iproj /= Lzd%Lnlpspd(ilr)%nproj) stop '1:applyprojectors'
           end do
           istart_ck=istart_c
           if (ieorb == orbs%norbp) exit loop_kpt
           ikpt=ikpt+1
           ispsi_k=ispsi
        end do loop_kpt
! Test disabled for now
!        if (istart_ck-1  /= nlpspd%nprojel) stop '2:applyprojectors'

     end do

     !restore the projectors in the proj array (for on the run forces calc., tails or so)
     ! Chosen to be in the global box.
     if (refill) then
        call Linearfill_projectors(iproc,Lzd%Glr,hx,hy,hz,at,orbs,rxyz,Lzd%Gnlpspd,proj,0)
     end if

  end if

  allocate(fxyz_orb(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fxyz_orb,'fxyz_orb',subname)

  !apply the projectors  k-point of the processor
  !starting k-point
  ikpt=orbs%iokpt(1)
  jorb=0
  loop_kptF: do

     call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

     call ncplx_kpt(ikpt,orbs,ncplx)

     ! loop over all my orbitals for calculating forces
     do iorb=isorb,ieorb
        ilr = orbs%inWhichLocreg(iorb+orbs%isorb)   !check for the locreg

        ! loop over all projectors
        call razero(3*at%nat,fxyz_orb)
        do ispinor=1,nspinor,ncplx
           jorb=jorb+1
           do iat=1,at%nat
              if (Lzd%Llr(ilr)%projflg(iat) == 0) cycle
              ityp=at%iatype(iat)
              do l=1,4
                 do i=1,3
                    if (at%psppar(l,i,ityp) /= 0.0_gp) then
                       do m=1,2*l-1
                          do icplx=1,ncplx
                             ! scalar product with the derivatives in all the directions
                             sp0=real(scalprod(icplx,0,m,i,l,iat,jorb),gp)
                             do idir=1,3
                                spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                fxyz_orb(idir,iat)=fxyz_orb(idir,iat)+&
                                     at%psppar(l,i,ityp)*sp0*spi
                             end do
                          end do
                       end do
                    end if
                 end do
              end do
              !HGH case, offdiagonal terms
              if (at%npspcode(ityp) == 3 .or. at%npspcode(ityp) == 10) then
                 do l=1,3 !no offdiagoanl terms for l=4 in HGH-K case
                    do i=1,2
                       if (at%psppar(l,i,ityp) /= 0.0_gp) then
                          loop_j: do j=i+1,3
                             if (at%psppar(l,j,ityp) == 0.0_gp) exit loop_j
                             !offdiagonal HGH term
                             if (at%npspcode(ityp) == 3) then !traditional HGH convention
                                hij=offdiagarr(i,j-i,l)*at%psppar(l,j,ityp)
                             else !HGH-K convention
                                hij=at%psppar(l,i+j+1,ityp)
                             end if
                             do m=1,2*l-1
                                !F_t= 2.0*h_ij (<D_tp_i|psi><psi|p_j>+<p_i|psi><psi|D_tp_j>)
                                !(the two factor is below)
                                do icplx=1,ncplx
                                   sp0i=real(scalprod(icplx,0,m,i,l,iat,jorb),gp)
                                   sp0j=real(scalprod(icplx,0,m,j,l,iat,jorb),gp)
                                   do idir=1,3
                                      spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                      spj=real(scalprod(icplx,idir,m,j,l,iat,jorb),gp)
                                      fxyz_orb(idir,iat)=fxyz_orb(idir,iat)+&
                                           hij*(sp0j*spi+spj*sp0i)
                                   end do
                                end do
                             end do
                          end do loop_j
                       end if
                    end do
                 end do
              end if
           end do
        end do

        !orbital-dependent factor for the forces
        orbfac=orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*2.0_gp

        ! For the basis functions, use coeff to build the linear combinations
        do iat=1,at%nat
           fsep(1,iat)=fsep(1,iat)+orbfac*fxyz_orb(1,iat)
           fsep(2,iat)=fsep(2,iat)+orbfac*fxyz_orb(2,iat)
           fsep(3,iat)=fsep(3,iat)+orbfac*fxyz_orb(3,iat)
        end do

     end do
     if (ieorb == orbs%norbp) exit loop_kptF
     ikpt=ikpt+1
     ispsi_k=ispsi
  end do loop_kptF


!!!  do iat=1,at%nat
!!!     write(20+iat,'(1x,i5,1x,3(1x,1pe12.5))') &
!!!          iat,fsep(1,iat),fsep(2,iat),fsep(3,iat)
!!!  end do

  i_all=-product(shape(fxyz_orb))*kind(fxyz_orb)
  deallocate(fxyz_orb,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz_orb',subname)

  i_all=-product(shape(scalprod))*kind(scalprod)
  deallocate(scalprod,stat=i_stat)
  call memocc(i_stat,i_all,'scalprod',subname)

END SUBROUTINE Linearnonlocal_forces

!>   Fill the proj array with the PSP projectors or their derivatives, following idir value
subroutine Linearfill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,idir)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,idir
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors),intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(out) :: proj
  !Local variables
  integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: istart_c,iat,iproj,nwarnings,ikpt,iskpt,iekpt

  if (iproc.eq.0 .and. nlpspd%nproj /=0 .and. idir==0) write(*,'(1x,a)',advance='no') &
       'Calculating wavelets expansion of projectors...'
  !warnings related to the projectors norm
  nwarnings=0
  !allocate these vectors up to the maximum size we can get
  istart_c=1

  !create projectors for any of the k-point hosted by the processor
  !starting kpoint
  if (orbs%norbp > 0) then
     iskpt=orbs%iokpt(1)
     iekpt=orbs%iokpt(orbs%norbp)
  else
     iskpt=1
     iekpt=1
  end if

  do ikpt=iskpt,iekpt
     iproj=0
     do iat=1,at%nat
        if (lr%projflg(iat) == 0) cycle
        !this routine is defined to uniformise the call for on-the-fly application
        call Linearatom_projector(ikpt,iat,idir,istart_c,iproj,&
             lr,hx,hy,hz,rxyz,at,orbs,nlpspd,proj,nwarnings)
     enddo
     if (iproj /= nlpspd%nproj) stop 'incorrect number of projectors created'
     ! projector part finished
  end do
  if (istart_c-1 /= nlpspd%nprojel) stop 'incorrect once-and-for-all psp generation'

  if (iproc == 0 .and. nlpspd%nproj /=0 .and. idir == 0) then
     if (nwarnings == 0) then
        write(*,'(1x,a)')'done.'
     else
        write(*,'(1x,a,i0,a)')'found ',nwarnings,' warnings.'
        write(*,'(1x,a)')'Some projectors may be too rough.'
        write(*,'(1x,a,f6.3)')&
             'Consider the possibility of modifying hgrid and/or the localisation radii.'
     end if
  end if

END SUBROUTINE Linearfill_projectors

subroutine Linearatom_projector(ikpt,iat,idir,istart_c,iproj,&
     lr,hx,hy,hz,rxyz,at,orbs,nlpspd,proj,nwarnings)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iat,idir,ikpt
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors),intent(in) :: lr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, intent(inout) :: istart_c,iproj,nwarnings
  real(wp), dimension(nlpspd%nprojel), intent(out) :: proj
  !Local variables
  integer :: ityp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,jseg_c,l,i,ncplx
  real(gp) :: kx,ky,kz

  !features of the k-point ikpt
  kx=orbs%kpts(1,ikpt)
  ky=orbs%kpts(2,ikpt)
  kz=orbs%kpts(3,ikpt)

  !evaluate the complexity of the k-point
  if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
     ncplx=1
  else
     ncplx=2
  end if

  ityp=at%iatype(iat)
  mbvctr_c=nlpspd%nvctr_p(2*iat-1)!-nlpspd%nvctr_p(2*iat-2)
  mbvctr_f=nlpspd%nvctr_p(2*iat  )!-nlpspd%nvctr_p(2*iat-1)

  mbseg_c=nlpspd%nseg_p(2*iat-1)!-nlpspd%nseg_p(2*iat-2)
  mbseg_f=nlpspd%nseg_p(2*iat  )!-nlpspd%nseg_p(2*iat-1)

  jseg_c = 1
  do l=1,iat-1
     jseg_c = jseg_c +  nlpspd%nseg_p(2*l - 1) + nlpspd%nseg_p(2*l) 
  end do

  !decide the loop bounds
  do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
     do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
        if (at%psppar(l,i,ityp) /= 0.0_gp) then
           call local_projector(at%geocode,at%atomnames(ityp),iat,idir,l,i,&
                at%psppar(l,0,ityp),rxyz(1,iat),lr,&
                hx,hy,hz,kx,ky,kz,ncplx,&
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                nlpspd%keyv_p(jseg_c),nlpspd%keyg_p(1,jseg_c),proj(istart_c),nwarnings)
           iproj=iproj+2*l-1
           istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx
           !print *,'iproc,istart_c,nlpspd%nprojel',istart_c,nlpspd%nprojel,ncplx
           if (istart_c > nlpspd%nprojel+1) stop 'istart_c > nprojel+1'
        endif
     enddo
  enddo
END SUBROUTINE Linearatom_projector

