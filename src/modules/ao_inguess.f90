!> @file
!! Medium-level routines associated to the generation of Atomic Orbitals inputguess
!! wavefunctions
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!> Handling of input guess creation from basis of atomic orbitals
module ao_inguess
  use module_base, only: gp,memocc,f_err_raise,ndebug,to_zero

  implicit none

  
  integer, parameter :: nmax_ao=6 !<maximum allowed value of principal quantum number for the electron configuration
  integer, parameter :: lmax_ao=3 !<maximum value of the angular momentum for the electron configuration
  integer, parameter :: nelecmax_ao=32 !<size of the interesting values of the compressed atomic input polarization

  private:: nmax_ao,lmax_ao,nelecmax_ao

contains

  subroutine iguess_generator(izatom,ielpsp,zion,psppar,npspcode,ngv,ngc,nlccpar,ng,nl,&
       &   nmax_occ,noccmax,lmax,occup,expo,psiat,enlargerprb,quartic_prefactor,gaenes_aux)
    implicit none
    logical, intent(in) :: enlargerprb
    integer, intent(in) :: ng,npspcode,nmax_occ,lmax,noccmax,ielpsp,izatom,ngv,ngc
    real(gp), intent(in) :: zion
    integer, dimension(lmax+1), intent(in) :: nl
    real(gp), dimension(0:4,0:6), intent(in) :: psppar
    real(gp), dimension(0:4,max((ngv*(ngv+1)/2)+(ngc*(ngc+1)/2),1)), intent(in) :: nlccpar
    real(gp), dimension(noccmax,lmax+1), intent(in) :: occup
    real(gp), dimension(ng+1), intent(out) :: expo
    real(gp), dimension(ng+1,nmax_occ), intent(out) :: psiat
    real(gp),intent(in),optional:: quartic_prefactor
    real(gp), dimension(nmax_occ),intent(out), optional :: gaenes_aux
    !local variables
    character(len=*), parameter :: subname='iguess_generator'
    integer, parameter :: n_int=100
    real(gp), parameter :: fact=4.0_gp
    !character(len=2) :: symbol
    integer :: lpx!,nsccode,mxpl,mxchg
    integer :: l,i,j,iocc,i_all,i_stat,iorder
    real(gp) :: alpz,alpl,rprb,rcov,rij,a,a0,a0in,tt!,ehomo,amu
    !integer, dimension(6,4) :: neleconf
    !real(kind=8), dimension(6,4) :: neleconf
    real(gp), dimension(4) :: gpot
    real(gp), dimension(noccmax,lmax+1) :: aeval,chrg,res
    real(gp), dimension(:), allocatable :: xp,alps
    real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef
    real(gp), dimension(:,:,:), allocatable :: psi
    real(gp), dimension(:,:,:,:), allocatable :: rmt

    !filename = 'psppar.'//trim(atomname)
    if (present(gaenes_aux)) call to_zero(nmax_occ,gaenes_aux(1))


    lpx=0
    lpx_determination: do i=1,4
       if (psppar(i,0) == 0.0_gp) then
          exit lpx_determination
       else
          lpx=i-1
       end if
    end do lpx_determination

    allocate(alps(lpx+1+ndebug),stat=i_stat)
    call memocc(i_stat,alps,'alps',subname)
    allocate(hsep(6,lpx+1+ndebug),stat=i_stat)
    call memocc(i_stat,hsep,'hsep',subname)

    !assignation of radii and coefficients of the local part
    alpz=psppar(0,0)
    alpl=psppar(0,0)
    alps(1:lpx+1)=psppar(1:lpx+1,0)
    gpot(1:4)=psppar(0,1:4)

    !assignation of the coefficents for the nondiagonal terms
    if (npspcode == 2) then !GTH case
       do l=1,lpx+1
          hsep(1,l)=psppar(l,1)
          hsep(2,l)=0.0_gp
          hsep(3,l)=psppar(l,2)
          hsep(4,l)=0.0_gp
          hsep(5,l)=0.0_gp
          hsep(6,l)=psppar(l,3)
       end do
    else if (npspcode == 3) then !HGH case
       allocate(ofdcoef(3,4+ndebug),stat=i_stat)
       call memocc(i_stat,ofdcoef,'ofdcoef',subname)

       ofdcoef(1,1)=-0.5_gp*sqrt(3._gp/5._gp) !h2
       ofdcoef(2,1)=0.5_gp*sqrt(5._gp/21._gp) !h4
       ofdcoef(3,1)=-0.5_gp*sqrt(100.0_gp/63._gp) !h5

       ofdcoef(1,2)=-0.5_gp*sqrt(5._gp/7._gp) !h2
       ofdcoef(2,2)=1._gp/6._gp*sqrt(35._gp/11._gp) !h4
       ofdcoef(3,2)=-7._gp/3._gp*sqrt(1._gp/11._gp) !h5

       ofdcoef(1,3)=-0.5_gp*sqrt(7._gp/9._gp) !h2
       ofdcoef(2,3)=0.5_gp*sqrt(63._gp/143._gp) !h4
       ofdcoef(3,3)=-9._gp*sqrt(1._gp/143._gp) !h5

       ofdcoef(1,4)=0.0_gp !h2
       ofdcoef(2,4)=0.0_gp !h4
       ofdcoef(3,4)=0.0_gp !h5

       !define the values of hsep starting from the pseudopotential file
       do l=1,lpx+1
          hsep(1,l)=psppar(l,1)
          hsep(2,l)=psppar(l,2)*ofdcoef(1,l)
          hsep(3,l)=psppar(l,2)
          hsep(4,l)=psppar(l,3)*ofdcoef(2,l)
          hsep(5,l)=psppar(l,3)*ofdcoef(3,l)
          hsep(6,l)=psppar(l,3)
       end do
       i_all=-product(shape(ofdcoef))*kind(ofdcoef)
       deallocate(ofdcoef,stat=i_stat)
       call memocc(i_stat,i_all,'ofdcoef',subname)
    else if (npspcode == 10 .or. npspcode == 7 .or. npspcode == 12) then !HGH-K case
       ! For PAW this is just the initial guess
       do l=1,lpx+1
          hsep(1,l)=psppar(l,1) !h11
          hsep(2,l)=psppar(l,4) !h12
          hsep(3,l)=psppar(l,2) !h22
          hsep(4,l)=psppar(l,5) !h13
          hsep(5,l)=psppar(l,6) !h23
          hsep(6,l)=psppar(l,3) !h33
       end do
    end if

    !!Just for extracting the covalent radius and rprb
    call atomic_info(izatom,ielpsp,rcov=rcov,rprb=rprb)
!    call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
    !!write(*,*) 'WARNING: multiply rprb with 5!!'
    !!rprb=rprb*5.d0


    if(present(quartic_prefactor)) then
       tt=rprb
       if(quartic_prefactor>0.d0) then
          ! There is a non-zero confinement
          rprb=(1.d0/(2.d0*quartic_prefactor))**.25d0
       else
          ! No confinement is used. Adjust rprb such that the quartic potential has at r=12 the same
          ! value as the parabolic potential
          rprb=144.d0**.25d0*tt
       end if
       !if(iproc==0) write(*,'(2(a,es12.3))') 'quartic potential for AO: modify rprb from ',tt,' to ',rprb
       !write(*,'(2(a,es12.3))') 'quartic potential for AO: modify rprb from ',tt,' to ',rprb
    end if

    if (enlargerprb) then
       !experimental
       rprb=100.0_gp
    end if

    !  occup(:,:)=0.0_gp
    !   do l=0,lmax-1
    !     iocc=0
    !     do i=1,6
    !        if (elecorbs(i,l+1) > 0.0_gp) then
    !           iocc=iocc+1
    !           !print *,'elecorbs',i,l,elecorbs(i,l+1),noccmax
    !            if (iocc > noccmax) stop 'iguess_generator: noccmax too small'
    !           occup(iocc,l+1)=elecorbs(i,l+1)
    !        endif
    !     end do
    !     nl(l+1)=iocc
    !  end do

    !allocate arrays for the gatom routine
    allocate(vh(4*(ng+1)**2,4*(ng+1)**2+ndebug),stat=i_stat)
    call memocc(i_stat,vh,'vh',subname)
    allocate(psi(0:ng,noccmax,lmax+ndebug),stat=i_stat)
    call memocc(i_stat,psi,'psi',subname)
    allocate(xp(0:ng+ndebug),stat=i_stat)
    call memocc(i_stat,xp,'xp',subname)
    allocate(rmt(n_int,0:ng,0:ng,lmax+ndebug),stat=i_stat)
    call memocc(i_stat,rmt,'rmt',subname)

    !can be switched on for debugging
    !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
    !     'Input Guess Generation for atom',trim(atomname),&
    !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb

    rij=3._gp
    ! exponents of gaussians
    a0in=alpz
    a0=a0in/rij
    !       tt=sqrt(sqrt(2._gp))
    tt=2._gp**.3_gp
    do i=0,ng
       a=a0*tt**i
       xp(i)=.5_gp/a**2
    end do

    ! initial guess
    do l=0,lmax-1
       do iocc=1,noccmax
          do i=0,ng
             psi(i,iocc,l+1)=0.0_gp
          end do
       end do
    end do

    call crtvh(ng,lmax-1,xp,vh,rprb,fact,n_int,rmt)
    if(present(quartic_prefactor)) then
       iorder=4
    else
       iorder=2
    end if
    call gatom(rcov,rprb,lmax-1,lpx,noccmax,occup,&
         zion,alpz,gpot,alpl,hsep,alps,ngv,ngc,nlccpar,vh,xp,rmt,fact,n_int,&
         aeval,ng,psi,res,chrg,iorder)
    
    !post-treatment of the inguess data
    do i=1,ng+1
       expo(i)=sqrt(0.5_gp/xp(i-1))
    end do

    i=0
    do l=1,4
       do iocc=1,nl(l)
          i=i+1
          !occupat(i)=occup(iocc,l)
          do j=1,ng+1
             psiat(j,i)=psi(j-1,iocc,l)
             if (present(gaenes_aux)) gaenes_aux(i) = aeval(iocc,l)
          end do
       end do
    end do

    i_all=-product(shape(vh))*kind(vh)
    deallocate(vh,stat=i_stat)
    call memocc(i_stat,i_all,'vh',subname)
    i_all=-product(shape(psi))*kind(psi)
    deallocate(psi,stat=i_stat)
    call memocc(i_stat,i_all,'psi',subname)
    i_all=-product(shape(xp))*kind(xp)
    deallocate(xp,stat=i_stat)
    call memocc(i_stat,i_all,'xp',subname)
    i_all=-product(shape(rmt))*kind(rmt)
    deallocate(rmt,stat=i_stat)
    call memocc(i_stat,i_all,'rmt',subname)
    i_all=-product(shape(hsep))*kind(hsep)
    deallocate(hsep,stat=i_stat)
    call memocc(i_stat,i_all,'hsep',subname)
    i_all=-product(shape(alps))*kind(alps)
    deallocate(alps,stat=i_stat)
    call memocc(i_stat,i_all,'alps',subname)

  END SUBROUTINE iguess_generator

  !> retrieve the information from the atom.
  !! different information can be obtained according to the usage which is needed
  subroutine atomic_info(zatom,zion,symbol,elconf,amu,rcov,rprb,ehomo,nsccode,maxpol,maxchg)
    use yaml_output, only: yaml_toa
    implicit none
    ! Arguments
    integer, intent(in) :: zatom            !< Z number of atom
    integer, intent(in) :: zion             !< Number of valence electrons of the ion (PSP should be in agreement)
    character(len=2), intent(out), optional :: symbol  !< Atomic symbol of Z, from the periodic table of elements
    double precision, intent(out), optional :: rcov        !< Covalent radius, atomic units
    double precision, intent(out), optional :: rprb        !< Parabolic radius for the input guess, of interest in the subroutine "gatom"
    double precision, intent(out), optional :: ehomo       !< Highest occupied molecular orbital energy, atomic units,
    !! See <a>http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html</a>
    double precision, intent(out), optional :: amu         !< Atomic mass unit (use values coming from ABINIT/11util/atmdata.F90)
    double precision, dimension(:,0:), optional :: elconf            !< Occupation number (electron configuration of the PSP atom)
                                                           !! assumed-shape, intent(out) not specified as did not found in the norm if it is legal
    integer, intent(out), optional :: nsccode !< Semicore orbitals, indicated as an integer.
    !! The integer is the n_s + 4*n_p + 16* n_d + 64* n_f
    !! where n_l are the number of semicore orbitals for a given angular momentum
    !! starting from the lower level of course
    integer, intent(out), optional :: maxpol    !< Maximum spin polarisation to be placed on the atom
    integer, intent(out), optional :: maxchg   !< Maximum charge to be placed on the atom

    !local variables
    character(len=2) :: symbol_
    integer :: nsccode_,mxpl_,mxchg_
    double precision :: rprb_,ehomo_,rcov_,amu_
    double precision, dimension(nmax_ao,0:lmax_ao) :: releconf !<these dimensions have to be modified in the following

    !extract all the information from the tabulated values of eleconf-inc.f90 file
    call eleconf(zatom,zion,symbol_,rcov_,rprb_,ehomo_,releconf,nsccode_,mxpl_,mxchg_,amu_)

    !then assign the requested values
    if (present(elconf)) then
       if (f_err_raise(any(shape(elconf) /= shape(releconf)),&
            'Electron Configuration array has wrong shape, found '//&
            trim(yaml_toa(shape(elconf)))//', needed '//&
            trim(yaml_toa(shape(releconf)))//'.',&
            err_name='BIGDFT_RUNTIME_ERROR')) return
       elconf=releconf
    end if

    if (present(symbol)) symbol=symbol_
    if (present(amu))       amu=amu_
    if (present(rcov))     rcov=rcov_
    if (present(rprb))     rprb=rprb_
    if (present(ehomo))   ehomo=ehomo_
    if (present(nsccode)) nsccode=nsccode_
    if (present(maxpol))   maxpol=mxpl_
    if (present(maxchg))   maxpol=mxchg_
    
  end subroutine atomic_info

  !> Correct the electronic configuration for a given atomic charge
  subroutine correct_semicore(nmax,lmax,ichg,neleconf,eleconf,nsccode)
    use module_base
    implicit none
    integer, intent(in) :: nmax,lmax,ichg
    real(kind=8) , dimension(nmax,0:lmax), intent(in) :: neleconf
    !integer, dimension(nmax,0:lmax), intent(in) :: neleconf
    real(gp), dimension(nmax,0:lmax), intent(out) :: eleconf
    integer, intent(inout) :: nsccode
    !local variables
    logical :: inocc
    integer :: i,l,nchgres,ichgp,nlsc
    real(gp) :: atchg

    !convert the array in real numbers
    do i=1,nmax
       do l=0,lmax
          eleconf(i,l)=real(neleconf(i,l),gp)
       end do
    end do

    nchgres=ichg !residual charge
    if (ichg >0) then
       !place the charge on the atom starting from the non closed shells
       do i=nmax,1,-1
          do l=lmax,0,-1
             if (neleconf(i,l) /= 2*(2*l+1) .and. neleconf(i,l) /= 0) then
                ichgp=min(nint(neleconf(i,l)),nchgres)
                nchgres=nchgres-ichgp
                eleconf(i,l)=eleconf(i,l)-real(ichgp,gp)
             end if
          end do
       end do
       if (nchgres /= 0) then
          !localise the highest occupied shell and charge it
          do i=nmax,1,-1
             do l=lmax,0,-1
                if (nint(eleconf(i,l)) == 2*(2*l+1)) then
                   ichgp=min(nint(eleconf(i,l)),nchgres)
                   nchgres=nchgres-ichgp
                   eleconf(i,l)=eleconf(i,l)-real(ichgp,gp)
                end if
             end do
          end do
          !!        !charge only unoccupied shells 
          !!        print *,'Atom ',symbol,': cannot charge occupied shells for the moment'
          !!        stop
       end if
    else if (ichg < 0) then
       !place the charge on the atom starting from the non closed shells
       do i=nmax,1,-1
          do l=lmax,0,-1
             if (neleconf(i,l) /= 0) then
                ichgp=min(2*(2*l+1)-nint(neleconf(i,l)),-nchgres)
                nchgres=nchgres+ichgp
                eleconf(i,l)=eleconf(i,l)+real(ichgp,gp)
             end if
          end do
       end do
       if (nchgres /= 0) then
          !localise the highest unoccupied shell and charge it
          inocc=.false.
          do i=1,nmax
             do l=0,lmax
                !once found the first occupied shell search for the first unoccpied
                if (inocc .and. nint(eleconf(i,l)) == 0) then
                   ichgp=min(2*(2*l+1),-nchgres)
                   nchgres=nchgres+ichgp
                   eleconf(i,l)=eleconf(i,l)+real(ichgp,gp)
                end if
                inocc=eleconf(i,l) /= 0.0_gp
             end do
          end do
          !!        !charge only occupied shells 
          !!        print *,'Atom ',symbol,': cannot charge unoccupied shells for the moment'
          !!        stop
       end if

    end if

    atchg=0.0_gp
    if (ichg /= 0) then
       !correct the semicore informations for a charged atom
       nsccode=0
       do l=0,lmax
          nlsc=0
          do i=1,nmax
             atchg=atchg+eleconf(i,l)
             if (eleconf(i,l) == real(2*(2*l+1),gp)) then
                nlsc=nlsc+1
                !if (nlsc <= 2) nsccode=nsccode+4**l
             end if
          end do
       end do
       if (atchg==0.0_gp) then
          write(*,*)'ERROR: an Atom must have input charge'
          stop
       end if
    end if



!!!  !if the atom has only closed shells we can treat it as semicore atom (commented)
!!!  isccode=nsccode
!!!  do l=lmax,0,-1
!!!     !control whether it is already semicore
!!!     itmp=isccode/((lmax+1)**l)
!!!     isccode=isccode-itmp*((lmax+1)**l)
!!!     !print *,'symbol',symbol,l,itmp,isccode,itmp*(lmax**l)
!!!     do i=1,nmax
!!!        if (neleconf(i,l) == 2*(2*l+1)) then
!!!           if (itmp==1) then
!!!              itmp=0
!!!              cycle
!!!           else
!!!               nsccode=nsccode+4**l !the maximum occupied is noccmax=2
!!!           end if
!!!        end if
!!!     end do
!!!  end do
  END SUBROUTINE correct_semicore

  !>  Calculate the occupation number for any of the orbitals
  subroutine at_occnums(ipolres,nspin,nspinor,nmax,lmax,nelecmax,eleconf,occupIG)
    use module_base
    implicit none
    integer, intent(in) :: nspinor,nspin,nmax,lmax,nelecmax
    real(gp), dimension(nmax,lmax), intent(in) :: eleconf
    integer, intent(inout) :: ipolres
    real(gp), dimension(nelecmax), intent(out) :: occupIG
    !local variables
    logical :: polarised
    integer :: iocc,ipolorb,norbpol_nc,i,l,m,noncoll,icoll,ispin, ipolsign
    real(gp) :: shelloccup,occshell,occres,rnl

    !in the non-collinear case the number of orbitals doubles
    if (nspinor == 4) then
       noncoll=2
    else
       noncoll=1
    end if

    call razero(nelecmax,occupIG)
    !call to_zero(nelecmax, occupIG(1))

    !here we should define the array of the occupation numbers
    !such array can then be redefined on the parent routines and then used as input
    iocc=0
    polarised=.false.
    !the sign is always the same
    if (ipolres >= 0) then
       ipolsign=1
    else
       ipolsign=-1
    end if
    do l=1,lmax
       iocc=iocc+1
       rnl=0.0_gp !real since it goes in occupIG
       do i=1,nmax
          if (eleconf(i,l) > 0.0_gp) then
             rnl=rnl+1.0_gp
          endif
       end do
       occupIG(iocc)=rnl
       !print *,'rnl,l',l,rnl,eleconf(:,l)
       do i=1,nmax
          if (eleconf(i,l) > 0.0_gp) then  
             shelloccup=eleconf(i,l)
             !decide the polarisation of the orbital by changing the population
             if (nint(shelloccup) /=  2*(2*l-1) ) then
                !this is a polarisable orbital
                polarised=.true.
                !assuming that the control of the allowed polarisation is already done

                ipolorb=ipolsign*min(abs(ipolres),  ((2*l-1) - abs( (2*l-1)- int(shelloccup) ) )  )
                ipolres=ipolres-ipolorb
             else
                !check for odd values of the occupation number
                if (mod(nint(shelloccup),2) /= 0) then
                   write(*,'(1x,a)')&
                        &   'The occupation number in the case of closed shells must be even'
                   stop
                end if
             end if

             if( polarised .AND. nspinor==4 .and. ipolorb /=0) then
                stop " in non-collinear case at_moments must be used for polarising, not natpol input"  
             endif

             do ispin=1,nspin
                occshell=shelloccup                 
                if (nspin==2 .or. nspinor==4) then
                   if (polarised) then
                      occshell=0.5_gp*(occshell+real(1-2*(ispin-1),gp)*ipolorb)
                   else
                      occshell=0.5_gp*occshell
                   end if
                end if

                !residue for the occupation number, to be used for
                !non-collinear case 
                occres=occshell
                !number of orbitals which will be polarised in this shell
                norbpol_nc=2*l-1
                do m=1,2*l-1
                   !each orbital has two electrons in the case of the 
                   !non-collinear case
                   do icoll=1,noncoll !non-trivial only for nspinor=4
                      iocc=iocc+1
                      !the occupation number rule changes for non-collinear
                      if (nspinor == 4) then
                         !for each orbital of the shell, use the Hund rule
                         !for determining the occupation
                         !if the occupation is one the orbital is not polarised
                         !otherwise it can be polarised via the polarisation
                         !indicated by atmoments
                         if (ceiling(occres) >= real(2*l-1,gp)) then
                            occupIG(iocc)=1.0_gp
                            if (icoll==2) then
                               occres=occres-1.0_gp
                               norbpol_nc=norbpol_nc-1
                            end if
                         else
                            if (icoll ==1) then
                               occupIG(iocc)=2.0_gp*occres/real(norbpol_nc,gp)
                            else
                               occupIG(iocc)=0.0_gp
                            end if
                         end if
                      else
                         occupIG(iocc)=occshell/real(2*l-1,gp)
                      end if
                   end do
                end do
             end do
          end if
       end do
    end do
  END SUBROUTINE at_occnums


  include 'eleconf-inc.f90'

end module ao_inguess
