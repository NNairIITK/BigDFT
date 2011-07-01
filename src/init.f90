!> @file
!!  Routines to initialize the information about localisation regions
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>   Calculates the descriptor arrays and nvctrp
!!   Calculates also the bounds arrays needed for convolutions
!!   Refers this information to the global localisation region descriptor
subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
     crmult,frmult,Glr,output_grid)
  use module_base
  use module_types
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: iproc
  real(gp), intent(in) :: hx,hy,hz,crmult,frmult
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  type(locreg_descriptors), intent(inout) :: Glr
  logical, intent(in), optional :: output_grid
  !local variables
  character(len=*), parameter :: subname='createWavefunctionsDescriptors'
  integer :: i_all,i_stat,i1,i2,i3,iat
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  logical :: my_output_grid
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  !assign the dimensions to improve (a little) readability
  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3
  nfl1=Glr%d%nfl1
  nfl2=Glr%d%nfl2
  nfl3=Glr%d%nfl3
  nfu1=Glr%d%nfu1
  nfu2=Glr%d%nfu2
  nfu3=Glr%d%nfu3

  !allocate kinetic bounds, only for free BC
  if (atoms%geocode == 'F') then
     allocate(Glr%bounds%kb%ibyz_c(2,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%kb%ibyz_c,'Glr%bounds%kb%ibyz_c',subname)
     allocate(Glr%bounds%kb%ibxz_c(2,0:n1,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%kb%ibxz_c,'Glr%bounds%kb%ibxz_c',subname)
     allocate(Glr%bounds%kb%ibxy_c(2,0:n1,0:n2+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%kb%ibxy_c,'Glr%bounds%kb%ibxy_c',subname)
     allocate(Glr%bounds%kb%ibyz_f(2,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%kb%ibyz_f,'Glr%bounds%kb%ibyz_f',subname)
     allocate(Glr%bounds%kb%ibxz_f(2,0:n1,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%kb%ibxz_f,'Glr%bounds%kb%ibxz_f',subname)
     allocate(Glr%bounds%kb%ibxy_f(2,0:n1,0:n2+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%kb%ibxy_f,'Glr%bounds%kb%ibxy_f',subname)
  end if

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------- Wavefunctions Descriptors Creation'
  end if

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_c,'logrid_c',subname)
  allocate(logrid_f(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_f,'logrid_f',subname)

  ! coarse grid quantities
  call fill_logrid(atoms%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%nat,&
       atoms%ntypes,atoms%iatype,rxyz,radii_cf(1,1),crmult,hx,hy,hz,logrid_c)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c,Glr%wfd%nvctr_c)
  if (iproc == 0) write(*,'(2(1x,a,i10))') &
       'Coarse resolution grid: Number of segments= ',Glr%wfd%nseg_c,'points=',Glr%wfd%nvctr_c

  if (atoms%geocode == 'F') then
     call make_bounds(n1,n2,n3,logrid_c,Glr%bounds%kb%ibyz_c,Glr%bounds%kb%ibxz_c,Glr%bounds%kb%ibxy_c)
  end if

  if (atoms%geocode == 'P' .and. .not. Glr%hybrid_on .and. Glr%wfd%nvctr_c /= (n1+1)*(n2+1)*(n3+1) ) then
     if (iproc ==0)then
        write(*,*)&
          ' ERROR: the coarse grid does not fill the entire periodic box'
        write(*,*)&
          '          errors due to translational invariance breaking may occur'
        !stop
     end if
     if (GPUconv) then
!        if (iproc ==0)then
           write(*,*)&
                '          The code should be stopped for a GPU calculation     '
           write(*,*)&
                '          since density is not initialised to 10^-20               '
!        end if
        stop
     end if
  end if

  call fill_logrid(atoms%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%nat,&
       atoms%ntypes,atoms%iatype,rxyz,radii_cf(1,2),frmult,hx,hy,hz,logrid_f)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f,Glr%wfd%nvctr_f)
  if (iproc == 0) write(*,'(2(1x,a,i10))') & 
       '  Fine resolution grid: Number of segments= ',Glr%wfd%nseg_f,'points=',Glr%wfd%nvctr_f
  if (atoms%geocode == 'F') then
     call make_bounds(n1,n2,n3,logrid_f,Glr%bounds%kb%ibyz_f,Glr%bounds%kb%ibxz_f,Glr%bounds%kb%ibxy_f)
  end if

  ! allocations for arrays holding the wavefunctions and their data descriptors
  call allocate_wfd(Glr%wfd,subname)

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1))

  ! fine grid quantities
  if (Glr%wfd%nseg_f > 0) then
     call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f,Glr%wfd%keyg(1,Glr%wfd%nseg_c+1), &
          & Glr%wfd%keyv(Glr%wfd%nseg_c+1))
  end if

! Create the file grid.xyz to visualize the grid of functions
  my_output_grid = .false.
  if (present(output_grid)) my_output_grid = output_grid
  if (my_output_grid) then
     open(unit=22,file='grid.xyz',status='unknown')
     write(22,*) Glr%wfd%nvctr_c+Glr%wfd%nvctr_f+atoms%nat,' atomic'
     if (atoms%geocode=='F') then
        write(22,*)'complete simulation grid with low and high resolution points'
     else if (atoms%geocode =='S') then
        write(22,'(a,2x,3(1x,1pe24.17))')'surface',atoms%alat1,atoms%alat2,atoms%alat3
     else if (atoms%geocode =='P') then
        write(22,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%alat1,atoms%alat2,atoms%alat3
     end if
     do iat=1,atoms%nat
        write(22,'(a6,2x,3(1x,e12.5),3x)') &
             trim(atoms%atomnames(atoms%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
     enddo
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  g ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     end do
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  G ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     enddo
     close(22)
  endif

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c',subname)
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f',subname)

  !for free BC admits the bounds arrays
  if (atoms%geocode == 'F') then

     !allocate grow, shrink and real bounds
     allocate(Glr%bounds%gb%ibzxx_c(2,0:n3,-14:2*n1+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%gb%ibzxx_c,'Glr%bounds%gb%ibzxx_c',subname)
     allocate(Glr%bounds%gb%ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%gb%ibxxyy_c,'Glr%bounds%gb%ibxxyy_c',subname)
     allocate(Glr%bounds%gb%ibyz_ff(2,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%gb%ibyz_ff,'Glr%bounds%gb%ibyz_ff',subname)
     allocate(Glr%bounds%gb%ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%gb%ibzxx_f,'Glr%bounds%gb%ibzxx_f',subname)
     allocate(Glr%bounds%gb%ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%gb%ibxxyy_f,'Glr%bounds%gb%ibxxyy_f',subname)

     allocate(Glr%bounds%sb%ibzzx_c(2,-14:2*n3+16,0:n1+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%sb%ibzzx_c,'Glr%bounds%sb%ibzzx_c',subname)
     allocate(Glr%bounds%sb%ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%sb%ibyyzz_c,'Glr%bounds%sb%ibyyzz_c',subname)
     allocate(Glr%bounds%sb%ibxy_ff(2,nfl1:nfu1,nfl2:nfu2+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%sb%ibxy_ff,'Glr%bounds%sb%ibxy_ff',subname)
     allocate(Glr%bounds%sb%ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%sb%ibzzx_f,'Glr%bounds%sb%ibzzx_f',subname)
     allocate(Glr%bounds%sb%ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%sb%ibyyzz_f,'Glr%bounds%sb%ibyyzz_f',subname)

     allocate(Glr%bounds%ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
     call memocc(i_stat,Glr%bounds%ibyyzz_r,'Glr%bounds%ibyyzz_r',subname)

     call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          Glr%bounds%kb%ibxy_c,Glr%bounds%sb%ibzzx_c,Glr%bounds%sb%ibyyzz_c,&
          Glr%bounds%kb%ibxy_f,Glr%bounds%sb%ibxy_ff,Glr%bounds%sb%ibzzx_f,Glr%bounds%sb%ibyyzz_f,&
          Glr%bounds%kb%ibyz_c,Glr%bounds%gb%ibzxx_c,Glr%bounds%gb%ibxxyy_c,&
          Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f,&
          Glr%bounds%ibyyzz_r)

  end if

  if ( atoms%geocode == 'P' .and. Glr%hybrid_on) then
     call make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,Glr%bounds,Glr%wfd)
     call make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          Glr%bounds%kb%ibxy_f,Glr%bounds%sb%ibxy_ff,Glr%bounds%sb%ibzzx_f,Glr%bounds%sb%ibyyzz_f,&
          Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f)
  endif

  !assign geocode and the starting points
  Glr%geocode=atoms%geocode

END SUBROUTINE createWavefunctionsDescriptors


!>   Determine localization region for all projectors, but do not yet fill the descriptor arrays
subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
     radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  type(nonlocal_psp_descriptors), intent(out) :: nlpspd
  real(wp), dimension(:), pointer :: proj
  !local variables
  character(len=*), parameter :: subname='createProjectorsArrays'
  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj
  integer :: iat,i_stat,i_all,iseg
  logical, dimension(:,:,:), allocatable :: logrid
  
  allocate(nlpspd%nseg_p(0:2*at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nseg_p,'nlpspd%nseg_p',subname)
  allocate(nlpspd%nvctr_p(0:2*at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nvctr_p,'nlpspd%nvctr_p',subname)
  allocate(nlpspd%nboxp_c(2,3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_c,'nlpspd%nboxp_c',subname)
  allocate(nlpspd%nboxp_f(2,3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_f,'nlpspd%nboxp_f',subname)

  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,radii_cf,&
       logrid,at,orbs,nlpspd)

  ! allocations for arrays holding the projectors and their data descriptors
  allocate(nlpspd%keyg_p(2,nlpspd%nseg_p(2*at%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%keyg_p,'nlpspd%keyg_p',subname)
  allocate(nlpspd%keyv_p(nlpspd%nseg_p(2*at%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%keyv_p,'nlpspd%keyv_p',subname)
  allocate(proj(nlpspd%nprojel+ndebug),stat=i_stat)
  call memocc(i_stat,proj,'proj',subname)

  ! After having determined the size of the projector descriptor arrays fill them
  do iat=1,at%nat
     call numb_proj(at%iatype(iat),at%ntypes,at%psppar,at%npspcode,mproj)
     if (mproj.ne.0) then 

        ! coarse grid quantities
        nl1=nlpspd%nboxp_c(1,1,iat) 
        nl2=nlpspd%nboxp_c(1,2,iat) 
        nl3=nlpspd%nboxp_c(1,3,iat) 

        nu1=nlpspd%nboxp_c(2,1,iat)
        nu2=nlpspd%nboxp_c(2,2,iat)
        nu3=nlpspd%nboxp_c(2,3,iat)
        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,3),cpmult,hx,hy,hz,logrid)

        iseg=nlpspd%nseg_p(2*iat-2)+1
        mseg=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)

        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
             logrid,mseg,nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg))

        ! fine grid quantities
        nl1=nlpspd%nboxp_f(1,1,iat)
        nl2=nlpspd%nboxp_f(1,2,iat)
        nl3=nlpspd%nboxp_f(1,3,iat)

        nu1=nlpspd%nboxp_f(2,1,iat)
        nu2=nlpspd%nboxp_f(2,2,iat)
        nu3=nlpspd%nboxp_f(2,3,iat)
        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
        iseg=nlpspd%nseg_p(2*iat-1)+1
        mseg=nlpspd%nseg_p(2*iat)-nlpspd%nseg_p(2*iat-1)
        if (mseg > 0) then
           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                logrid,mseg,nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg))
        end if
     endif
  enddo

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)

  !fill the projectors if the strategy is a distributed calculation
  if (.not. DistProjApply) then
     !calculate the wavelet expansion of projectors
     call fill_projectors(iproc,n1,n2,n3,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,0)
  end if

END SUBROUTINE createProjectorsArrays





!!****f* BigDFT/fillPcProjOnTheFly
!! FUNCTION
!>   fill the preconditioning projectors for a given atom 
!! SOURCE
!!
subroutine fillPcProjOnTheFly(PPD, Glr, iat, at, hx,hy,hz,startjorb,ecut_pc,   initial_istart_c ) 
  use module_interfaces
  use module_base
  use module_types
  implicit none
  character(len=11) :: orbname
  type(pcproj_data_type),  intent(in) ::PPD
  type(locreg_descriptors),  intent(in):: Glr
  integer, intent(in)  ::iat, startjorb
  real(gp), intent(in) ::  ecut_pc, hx,hy,hz
  !! real(gp), pointer :: gaenes(:)
  integer, intent(in) :: initial_istart_c
  type(atoms_data), intent(in) :: at

  ! local variables  
  type(locreg_descriptors) :: Plr
  real(gp) kx, ky, kz
  integer :: jorb, ncplx, istart_c
  real(wp), dimension(PPD%G%ncoeff ) :: Gocc
  character(len=*), parameter :: subname='fillPcProjOnTheFly'

  istart_c=initial_istart_c
  
  Plr%d%n1 = Glr%d%n1
  Plr%d%n2 = Glr%d%n2
  Plr%d%n3 = Glr%d%n3
  Plr%geocode = at%geocode
  

  Plr%wfd%nvctr_c  =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
  Plr%wfd%nvctr_f  =PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
  Plr%wfd%nseg_c   =PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
  Plr%wfd%nseg_f   =PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)
  
  call allocate_wfd(Plr%wfd,subname)
  
  
  Plr%wfd%keyv(:)  = PPD%pc_nlpspd%keyv_p(  PPD%pc_nlpspd%nseg_p(2*iat-2)+1:  PPD%pc_nlpspd%nseg_p(2*iat)   )
  Plr%wfd%keyg(1:2, :)  = PPD%pc_nlpspd%keyg_p( 1:2,  PPD%pc_nlpspd%nseg_p(2*iat-2)+1:  PPD%pc_nlpspd%nseg_p(2*iat)   )
  
  kx=0.0_gp
  ky=0.0_gp
  kz=0.0_gp
  
  Gocc=0.0_wp

  jorb=startjorb
  
  do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat) 
     if( PPD%gaenes(jorb)<ecut_pc) then
   
        Gocc(jorb)=1.0_wp
        
        ncplx=1
        call gaussians_to_wavelets_orb(ncplx,Plr,hx,hy,hz,kx,ky,kz,PPD%G,&
             Gocc(1),  PPD%pc_proj(istart_c)  )
        Gocc(jorb)=0.0_wp
        
        
        !! ---------------  use this to plot projectors
!!$              write(orbname,'(A,i4.4)')'pc_',iproj
!!$              Plr%bounds = Glr%bounds
!!$              Plr%d          = Glr%d
!!$              call plot_wf_cube(orbname,at,Plr,hx,hy,hz,PPD%G%rxyz, PPD%pc_proj(istart_c) ,"1234567890" ) 
        
        istart_c=istart_c + (   Plr%wfd%nvctr_c    +   7*Plr%wfd%nvctr_f   )


     endif
     jorb=jorb+1
     
     if(jorb> PPD%G%ncoeff) exit
     
  end do
  
  call deallocate_wfd(Plr%wfd,subname)
  
end subroutine fillPcProjOnTheFly
!!***


!!****f* BigDFT/fillPawProjOnTheFly
!! FUNCTION
!>   fill the preconditioning projectors for a given atom 
!! SOURCE
!!
subroutine fillPawProjOnTheFly(PAWD, Glr, iat,  hx,hy,hz,kx,ky,kz,startjorb,   initial_istart_c, geocode, at, iatat) 
  use module_interfaces
  use module_base
  use module_types
  implicit none
  character(len=11) :: orbname
  type(pawproj_data_type),  intent(in) ::PAWD
  type(locreg_descriptors),  intent(in):: Glr
  integer, intent(in)  ::iat, startjorb
  real(gp), intent(in) ::   hx,hy,hz,kx,ky,kz
  integer, intent(in) :: initial_istart_c
  character(len=1), intent(in) :: geocode
  type(atoms_data) :: at
  integer :: iatat

  ! local variables  
  type(locreg_descriptors) :: Plr
  integer :: jorb, ncplx, istart_c
  real(wp), dimension(PAWD%G%ncoeff ) :: Gocc
  character(len=*), parameter :: subname='fillPawProjOnTheFly'


  !!Just for extracting the covalent radius and rprb
  integer :: nsccode,mxpl,mxchg
  real(gp) ::amu,rprb,ehomo,rcov
  character(len=2) :: symbol
  real(kind=8), dimension(6,4) :: neleconf
  
  istart_c=initial_istart_c
  
  Plr%d%n1 = Glr%d%n1
  Plr%d%n2 = Glr%d%n2
  Plr%d%n3 = Glr%d%n3
  Plr%geocode = geocode
  
  Plr%wfd%nvctr_c  =PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
  Plr%wfd%nvctr_f  =PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
  Plr%wfd%nseg_c   =PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)
  Plr%wfd%nseg_f   =PAWD%paw_nlpspd%nseg_p(2*iat  )-PAWD%paw_nlpspd%nseg_p(2*iat-1)
  
  call allocate_wfd(Plr%wfd,subname)
  
  Plr%wfd%keyv(:)  = PAWD%paw_nlpspd%keyv_p(  PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:  PAWD%paw_nlpspd%nseg_p(2*iat)   )
  Plr%wfd%keyg(1:2, :)  = PAWD%paw_nlpspd%keyg_p( 1:2,  PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:  PAWD%paw_nlpspd%nseg_p(2*iat)   )
  
  if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
     ncplx=1
  else
     ncplx=2
  end if
  
  Gocc=0.0_wp

  jorb=startjorb
  
  !!Just for extracting the covalent radius 
  call eleconf(at%nzatom( at%iatype(iatat)), at%nelpsp(at%iatype(iatat)) ,  &
       symbol, rcov, rprb, ehomo,neleconf, nsccode, mxpl, mxchg, amu)
  !!  this 1.5 factor is the same as in file xabsorber.f90, routine GetExcitedOrbitalAsG
  rcov=rcov*1.5_gp

  do while( jorb<=PAWD%G%ncoeff         .and. PAWD%iorbtolr(jorb)== iat)      
     Gocc(jorb)=1.0_wp

     call gaussians_c_to_wavelets_orb(ncplx,Plr,hx,hy,hz,kx,ky,kz,PAWD%G,&
          Gocc(1),  PAWD%paw_proj(istart_c), rcov  )

     Gocc(jorb)=0.0_wp
!!$     !! ---------------  use this to plot projectors
!!$              write(orbname,'(A,i4.4)')'paw2_',jorb
!!$              Plr%bounds = Glr%bounds
!!$              Plr%d          = Glr%d
!!$              call plot_wf_cube(orbname,PAWD%at,Plr,hx,hy,hz,PAWD%G%rxyz, PAWD%paw_proj(istart_c) ,"1234567890" ) 

     istart_c=istart_c + (   Plr%wfd%nvctr_c    +   7*Plr%wfd%nvctr_f   ) * ncplx


     jorb=jorb+1

     if(jorb> PAWD%G%ncoeff) exit
     
  end do
  
  call deallocate_wfd(Plr%wfd,subname)
  
end subroutine fillPawProjOnTheFly
!!***






!!****f* BigDFT/createPcProjectorsArrays
!! FUNCTION
!>   Determine localization region for all preconditioning projectors, but do not yet fill the descriptor arrays
!! SOURCE
!!
subroutine createPcProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
     radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
     PPD, Glr)
  use module_interfaces
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs

  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  real(gp), intent(in):: ecut_pc

  type(pcproj_data_type) ::PPD

  type(locreg_descriptors),  intent(in):: Glr

  
  !local variables
  character(len=*), parameter :: subname='createPcProjectorsArrays'
  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj, mvctr
  integer :: iat,i_stat,i_all,iseg, istart_c
  logical, dimension(:,:,:), allocatable :: logrid


  integer :: ng
  logical :: enlargerprb
  real(wp), dimension(:), pointer :: Gocc

  integer, pointer :: iorbto_l(:)
  integer, pointer :: iorbto_m(:)
  integer, pointer :: iorbto_ishell(:)
  integer, pointer :: iorbto_iexpobeg(:)
  
  integer :: ncplx, nspin
  real(gp) :: kx,ky,kz
  integer ::  jorb
  integer :: iproj, startjorb
  type(locreg_descriptors) :: Plr
  character(len=11) :: orbname
  real(gp) Pcpmult
  integer :: mprojtot, nvctr_c, nvctr_f
  real(gp) rdum
  integer nprojel_tmp
  
  
  integer mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, &
       jseg_c, idum

  Pcpmult=1.5*cpmult

  ng=21
  enlargerprb = .false.
  nspin=1


  nullify(PPD%G%rxyz)
  call gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,PPD%G,Gocc, PPD%gaenes, &
      PPD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg  )  


  ! allocated  : gaenes, Gocc , PPD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg


!!$ ========================================================================================


  !---------

  allocate(PPD%pc_nlpspd%nseg_p(0:2*at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_nlpspd%nseg_p,'pc_nlpspd%nseg_p',subname)
  allocate(PPD%pc_nlpspd%nvctr_p(0:2*at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_nlpspd%nvctr_p,'pc_nlpspd%nvctr_p',subname)
  allocate(PPD%pc_nlpspd%nboxp_c(2,3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_nlpspd%nboxp_c,'pc_nlpspd%nboxp_c',subname)
  allocate(PPD%pc_nlpspd%nboxp_f(2,3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_nlpspd%nboxp_f,'pc_nlpspd%nboxp_f',subname)

  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)


  call localize_projectors(iproc,n1,n2,n3,hx,hy,hz,Pcpmult,fpmult,rxyz,radii_cf,&
       logrid,at,orbs,PPD%pc_nlpspd)

  ! the above routine counts atomic projector and the number of their element for psp
  ! We must therefore correct , later, nlpspd%nprojel  and nlpspd%nproj
  !-------------------

  ! allocations for arrays holding the projectors and their data descriptors


  allocate(PPD%pc_nlpspd%keyg_p(2,PPD%pc_nlpspd%nseg_p(2*at%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_nlpspd%keyg_p,'pc_nlpspd%keyg_p',subname)


  allocate(PPD%pc_nlpspd%keyv_p(PPD%pc_nlpspd%nseg_p(2*at%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_nlpspd%keyv_p,'pc_nlpspd%keyv_p',subname)



!!$  -- this one delayed, waiting for the correct pc_nlpspd%nprojel, pc_nlpspd%nproj
!!$  --
!!$  allocate(pc_proj(pc_nlpspd%nprojel+ndebug),stat=i_stat)
!!$  call memocc(i_stat,pc_proj,'pc_proj',subname)
  PPD%ecut_pc=ecut_pc

  PPD%pc_nlpspd%nprojel=0
  PPD%pc_nlpspd%nproj  =0

!!$ =========================================================================================  

  mprojtot=0
  jorb=1  
  ! After having determined the size of the projector descriptor arrays fill them
  do iat=1,at%nat

     mproj=0

     do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat)

        if( PPD%gaenes(jorb)<ecut_pc) then
           mproj=mproj+1
        endif
        if(jorb==PPD%G%ncoeff) exit
        jorb=jorb+1
     end do

     mprojtot=mprojtot+mproj

     PPD%pc_nlpspd%nproj=PPD%pc_nlpspd%nproj+mproj


     if (mproj.ne.0) then 

        nprojel_tmp=0

        ! coarse grid quantities
        nl1=PPD%pc_nlpspd%nboxp_c(1,1,iat) 
        nl2=PPD%pc_nlpspd%nboxp_c(1,2,iat) 
        nl3=PPD%pc_nlpspd%nboxp_c(1,3,iat) 
            
        nu1=PPD%pc_nlpspd%nboxp_c(2,1,iat)
        nu2=PPD%pc_nlpspd%nboxp_c(2,2,iat)
        nu3=PPD%pc_nlpspd%nboxp_c(2,3,iat)

        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,3),Pcpmult,hx,hy,hz,logrid)

        iseg=PPD%pc_nlpspd%nseg_p(2*iat-2)+1
        mseg=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)

        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
             logrid,mseg,PPD%pc_nlpspd%keyg_p(1,iseg),PPD%pc_nlpspd%keyv_p(iseg))

        mvctr =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)


        nprojel_tmp =nprojel_tmp +mproj*mvctr
        
        ! fine grid quantities
        nl1=PPD%pc_nlpspd%nboxp_f(1,1,iat)
        nl2=PPD%pc_nlpspd%nboxp_f(1,2,iat)
        nl3=PPD%pc_nlpspd%nboxp_f(1,3,iat)
            
        nu1=PPD%pc_nlpspd%nboxp_f(2,1,iat)
        nu2=PPD%pc_nlpspd%nboxp_f(2,2,iat)
        nu3=PPD%pc_nlpspd%nboxp_f(2,3,iat)
        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
        iseg=PPD%pc_nlpspd%nseg_p(2*iat-1)+1
        mseg=PPD%pc_nlpspd%nseg_p(2*iat)-PPD%pc_nlpspd%nseg_p(2*iat-1)
        if (mseg > 0) then
           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                logrid,mseg,PPD%pc_nlpspd%keyg_p(1,iseg),PPD%pc_nlpspd%keyv_p(iseg))

           mvctr =PPD%pc_nlpspd%nvctr_p(2*iat)-PPD%pc_nlpspd%nvctr_p(2*iat-1)
           
           nprojel_tmp=nprojel_tmp+mproj*mvctr*7

        end if
        
        if( PPD%DistProjApply)  then
           PPD%pc_nlpspd%nprojel=max(PPD%pc_nlpspd%nprojel,nprojel_tmp   )
        else
           PPD%pc_nlpspd%nprojel= PPD%pc_nlpspd%nprojel+nprojel_tmp 
        endif


     endif

  enddo

     
  allocate(PPD%pc_proj(PPD%pc_nlpspd%nprojel+ndebug),stat=i_stat)
  call memocc(i_stat,PPD%pc_proj,'pc_proj',subname)

  allocate(PPD%ilr_to_mproj(at%nat  +ndebug ) , stat=i_stat)
  call memocc(i_stat,PPD%ilr_to_mproj,'ilr_to_mproj',subname)

  allocate(PPD%iproj_to_ene(mprojtot +ndebug ) , stat=i_stat)
  call memocc(i_stat ,PPD%iproj_to_ene,'iproj_to_ene',subname)

  allocate(PPD%iproj_to_factor(mprojtot +ndebug ) , stat=i_stat)
  call memocc(i_stat ,PPD%iproj_to_factor,'iproj_to_factor',subname)

  allocate(PPD%iproj_to_l(mprojtot +ndebug ) , stat=i_stat)
  call memocc(i_stat ,PPD%iproj_to_l,'iproj_to_l',subname)
  
  PPD%mprojtot=mprojtot

  
  startjorb=1
  jorb=1
  istart_c=1
  Gocc(:)=0.0_wp


  iproj=0
  do iat=1,at%nat

     mproj=0
     do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat)
        if( PPD%gaenes(jorb)<ecut_pc) then
           mproj=mproj+1
        endif
        if(jorb==PPD%G%ncoeff) exit
        jorb=jorb+1
     end do

     PPD%ilr_to_mproj(iat)=mproj
     if( mproj>0) then
        nvctr_c  =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
        nvctr_f  =PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)

        jorb=startjorb
        do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat) 
           if( PPD%gaenes(jorb)<ecut_pc) then
              iproj=iproj+1
              PPD%iproj_to_ene(iproj) = PPD%gaenes(jorb)
              PPD%iproj_to_l(iproj)   = iorbto_l(jorb)

              istart_c=istart_c + (   nvctr_c    +   7*nvctr_f   )
           endif
           jorb=jorb+1

           if(jorb> PPD%G%ncoeff) exit
        end do
        

        if( .not. PPD%DistProjApply) then
           istart_c= istart_c-mproj*(nvctr_c+7*nvctr_f)

           call fillPcProjOnTheFly(PPD, Glr, iat, at, hx,hy,hz, startjorb,ecut_pc ,  istart_c ) 
           istart_c= istart_c+mproj*(nvctr_c+7*nvctr_f)

!!$
!!$           ncplx=1
!!$           rdum=0.0_gp
!!$
!!$           mbvctr_c=PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!!$           mbvctr_f=PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!!$           
!!$           mbseg_c=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
!!$           mbseg_f=PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)
!!$           jseg_c=PPD%pc_nlpspd%nseg_p(2*iat-2)+1
!!$              
!!$           do idum=1, 9
!!$              call wpdot_wrap(ncplx,  &
!!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PPD%pc_nlpspd%keyv_p(jseg_c),&
!!$                   PPD%pc_nlpspd%keyg_p(1,jseg_c),PPD%pc_proj(istart_c-idum*(nvctr_c+7*nvctr_f)),& 
!!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PPD%pc_nlpspd%keyv_p(jseg_c),&
!!$                   PPD%pc_nlpspd%keyg_p(1,jseg_c),&
!!$                   PPD%pc_proj(istart_c-idum*(nvctr_c+7*nvctr_f)),&
!!$                   rdum)
!!$           end do
        endif

     end if

     !! aggiunger condizione su istartc_c per vedere se e nelprj

     startjorb=jorb

  enddo

  if( .not. PPD%DistProjApply) then
     call deallocate_gwf(PPD%G,subname)
  endif

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)

  i_all=-product(shape(Gocc))*kind(Gocc)
  deallocate(Gocc,stat=i_stat)
  call memocc(i_stat,i_all,'Gocc',subname)

!!$  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
!!$  deallocate(iorbtolr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'iorbtolr',subname)

  i_all=-product(shape(iorbto_l))*kind(iorbto_l)
  deallocate(iorbto_l,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_l',subname)

  i_all=-product(shape(iorbto_m))*kind(iorbto_m)
  deallocate(iorbto_m,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_m',subname)

  i_all=-product(shape(iorbto_ishell))*kind(iorbto_ishell)
  deallocate(iorbto_ishell,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_ishell',subname)


  i_all=-product(shape(iorbto_iexpobeg))*kind(iorbto_iexpobeg)
  deallocate(iorbto_iexpobeg,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_iexpobeg',subname)


END SUBROUTINE createPcProjectorsArrays
!!***






!!****f* BigDFT/createPawProjectorsArrays
!! FUNCTION
!>   Determine localization region for all preconditioning projectors, but do not yet fill the descriptor arrays
!! SOURCE
!!
subroutine createPawProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
     radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
     PAWD, Glr)
  use module_interfaces
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs

  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  real(gp), intent(in):: ecut_pc

  type(PAWproj_data_type) ::PAWD

  type(locreg_descriptors),  intent(in):: Glr
  
  !local variables
  character(len=*), parameter :: subname='createPawProjectorsArrays'

  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj, mvctr
  integer :: iat,i_stat,i_all,iseg, istart_c
  logical, dimension(:,:,:), allocatable :: logrid


  integer :: ng
  logical :: enlargerprb
  real(wp), dimension(:), pointer :: Gocc

  integer, pointer :: iorbto_l(:)
  integer, pointer :: iorbto_paw_nchannels(:)
  integer, pointer :: iorbto_m(:)
  integer, pointer :: iorbto_ishell(:)
  integer, pointer :: iorbto_iexpobeg(:)
  
  integer :: ncplx, nspin
  real(gp) :: kx,ky,kz
  integer ::  jorb
  integer :: iproj, startjorb
  type(locreg_descriptors) :: Plr
  character(len=11) :: orbname
  real(gp) Pcpmult
  integer :: mprojtot, nvctr_c, nvctr_f
  real(gp) rdum
  integer iatat
  
  integer mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, &
       jseg_c, idum
  integer :: ikpt,iskpt,iekpt

  Pcpmult=1.0*cpmult


  nullify(PAWD%G%rxyz)

  call gaussian_pswf_basis_for_paw(iproc,nspin,at,rxyz,PAWD%G, &
      PAWD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg  ,&
      iorbto_paw_nchannels, PAWD%iprojto_imatrixbeg )  


  allocate(Gocc(PAWD%G%ncoeff+ndebug),stat=i_stat)
  call memocc(i_stat,Gocc,'Gocc',subname)
  call razero(PAWD%G%ncoeff,Gocc)

  ! allocated  : gaenes, Gocc , PAWD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg, iorbto_paw_nchannels

!!$ ========================================================================================
  !---------

  allocate(PAWD%paw_nlpspd%nseg_p(0:2*PAWD%G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_nlpspd%nseg_p,'pc_nlpspd%nseg_p',subname)
  allocate(PAWD%paw_nlpspd%nvctr_p(0:2*PAWD%G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_nlpspd%nvctr_p,'pc_nlpspd%nvctr_p',subname)
  allocate(PAWD%paw_nlpspd%nboxp_c(2,3,PAWD%G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_nlpspd%nboxp_c,'pc_nlpspd%nboxp_c',subname)
  allocate(PAWD%paw_nlpspd%nboxp_f(2,3,PAWD%G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_nlpspd%nboxp_f,'pc_nlpspd%nboxp_f',subname)

  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors_paw(iproc,n1,n2,n3,hx,hy,hz,Pcpmult,1*fpmult,rxyz,radii_cf,&
       logrid,at,orbs,PAWD)

  ! the above routine counts atomic projector and the number of their element for psp
  ! We must therefore correct , later, nlpspd%nprojel  and nlpspd%nproj
  !-------------------

  ! allocations for arrays holding the projectors and their data descriptors

  allocate(PAWD%paw_nlpspd%keyg_p(2,PAWD%paw_nlpspd%nseg_p(2*PAWD%G%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_nlpspd%keyg_p,'pc_nlpspd%keyg_p',subname)

  allocate(PAWD%paw_nlpspd%keyv_p(PAWD%paw_nlpspd%nseg_p(2*PAWD%G%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_nlpspd%keyv_p,'pc_nlpspd%keyv_p',subname)

!!$  -- this one delayed, waiting for the correct pc_nlpspd%nprojel, pc_nlpspd%nproj
!!$  --
!!$  allocate(pc_proj(pc_nlpspd%nprojel+ndebug),stat=i_stat)
!!$  call memocc(i_stat,pc_proj,'pc_proj',subname)
  allocate(PAWD%paw_proj(PAWD%paw_nlpspd%nprojel+ndebug),stat=i_stat)
  call memocc(i_stat,PAWD%paw_proj,'paw_proj',subname)

  allocate(PAWD%ilr_to_mproj(PAWD%G%nat  +ndebug ) , stat=i_stat)
  call memocc(i_stat,PAWD%ilr_to_mproj,'ilr_to_mproj',subname)

  allocate(PAWD%iproj_to_l(PAWD%paw_nlpspd%nproj +ndebug ) , stat=i_stat)
  call memocc(i_stat ,PAWD%iproj_to_l,'iproj_to_l',subname)

  allocate(PAWD%iproj_to_paw_nchannels( PAWD%paw_nlpspd%nproj+ndebug ) , stat=i_stat)
  call memocc(i_stat ,PAWD%iproj_to_paw_nchannels,'iproj_to_paw_nchannels',subname)

!!$ =========================================================================================  

  jorb=1  
  ! After having determined the size of the projector descriptor arrays fill them
  iat=0
  do iatat=1, at%nat
     if (  at%paw_NofL(at%iatype(iatat)).gt.0  ) then
        iat=iat+1
        mproj=0
        do while( jorb<=PAWD%G%ncoeff         .and. PAWD%iorbtolr(jorb)== iat)
           mproj=mproj+1
           if(jorb==PAWD%G%ncoeff) exit
           jorb=jorb+1
        end do

        PAWD%paw_nlpspd%nproj=PAWD%paw_nlpspd%nproj+mproj
        if (mproj.ne.0) then 

           ! coarse grid quantities
           nl1=PAWD%paw_nlpspd%nboxp_c(1,1,iat) 
           nl2=PAWD%paw_nlpspd%nboxp_c(1,2,iat) 
           nl3=PAWD%paw_nlpspd%nboxp_c(1,3,iat) 

           nu1=PAWD%paw_nlpspd%nboxp_c(2,1,iat)
           nu2=PAWD%paw_nlpspd%nboxp_c(2,2,iat)
           nu3=PAWD%paw_nlpspd%nboxp_c(2,3,iat)

           call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                at%ntypes,at%iatype(iatat),rxyz(1,iatat),radii_cf(1,3),Pcpmult,hx,hy,hz,logrid)

           iseg=PAWD%paw_nlpspd%nseg_p(2*iat-2)+1
           mseg=PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)

           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                logrid,mseg,PAWD%paw_nlpspd%keyg_p(1,iseg),PAWD%paw_nlpspd%keyv_p(iseg))
           mvctr =PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)

           ! fine grid quantities
           nl1=PAWD%paw_nlpspd%nboxp_f(1,1,iat)
           nl2=PAWD%paw_nlpspd%nboxp_f(1,2,iat)
           nl3=PAWD%paw_nlpspd%nboxp_f(1,3,iat)

           nu1=PAWD%paw_nlpspd%nboxp_f(2,1,iat)
           nu2=PAWD%paw_nlpspd%nboxp_f(2,2,iat)
           nu3=PAWD%paw_nlpspd%nboxp_f(2,3,iat)
           call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                at%ntypes,at%iatype(iatat),rxyz(1,iatat),radii_cf(1,2),1*fpmult,hx,hy,hz,logrid)
           iseg=PAWD%paw_nlpspd%nseg_p(2*iat-1)+1
           mseg=PAWD%paw_nlpspd%nseg_p(2*iat)-PAWD%paw_nlpspd%nseg_p(2*iat-1)
           if (mseg > 0) then
              call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                   logrid,mseg,PAWD%paw_nlpspd%keyg_p(1,iseg),PAWD%paw_nlpspd%keyv_p(iseg))
              mvctr =PAWD%paw_nlpspd%nvctr_p(2*iat)-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
           end if


        endif
     endif
  enddo

  if (orbs%norbp > 0) then
     iskpt=orbs%iokpt(1)
     iekpt=orbs%iokpt(orbs%norbp)
  else
     iskpt=1
     iekpt=1
  end if

  istart_c=1
  do ikpt=iskpt,iekpt     

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
     
     startjorb=1
     jorb=1
     Gocc(:)=0.0_wp
     iproj=0

     iat=0
     do iatat=1, at%nat
        if (  at%paw_NofL(at%iatype(iatat)).gt.0  ) then
           iat=iat+1
           mproj=0
           do while( jorb<=PAWD%G%ncoeff         .and. PAWD%iorbtolr(jorb)== iat)
              mproj=mproj+1
              if(jorb==PAWD%G%ncoeff) exit
              jorb=jorb+1
           end do

           PAWD%ilr_to_mproj(iat)=mproj
           if( mproj>0) then
              nvctr_c  =PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
              nvctr_f  =PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)

              jorb=startjorb
              do while( jorb<=PAWD%G%ncoeff  .and. PAWD%iorbtolr(jorb)== iat) 
                 iproj=iproj+1
                 PAWD%iproj_to_l(iproj)   = iorbto_l(jorb)
                 PAWD%iproj_to_paw_nchannels(iproj)   = iorbto_paw_nchannels(jorb)
                 istart_c=istart_c + (   nvctr_c    +   7*nvctr_f   )*ncplx
                 jorb=jorb+1
                 if(jorb> PAWD%G%ncoeff) exit
              end do
              if( .not. PAWD%DistProjApply) then
                 istart_c= istart_c-mproj*(nvctr_c+7*nvctr_f)*ncplx
                 call fillPawProjOnTheFly(PAWD, Glr, iat,  hx,hy,hz, kx,ky,kz, startjorb,&
                      istart_c, at%geocode , at, iatat) 
                 istart_c= istart_c+mproj*(nvctr_c+7*nvctr_f)*ncplx
              endif
           end if
           startjorb=jorb
        end if
     enddo
  enddo
  if (istart_c-1 /= PAWD%paw_nlpspd%nprojel) stop 'incorrect once-and-for-all psp generation'


  if( .not. PAWD%DistProjApply) then
     call deallocate_gwf_c(PAWD%G,subname)
  endif

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)

  i_all=-product(shape(Gocc))*kind(Gocc)
  deallocate(Gocc,stat=i_stat)
  call memocc(i_stat,i_all,'Gocc',subname)

!!$  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
!!$  deallocate(iorbtolr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'iorbtolr',subname)

  i_all=-product(shape(iorbto_l))*kind(iorbto_l)
  deallocate(iorbto_l,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_l',subname)

  i_all=-product(shape(iorbto_paw_nchannels))*kind(iorbto_paw_nchannels)
  deallocate(iorbto_paw_nchannels,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_paw_nchannels',subname)





  i_all=-product(shape(iorbto_m))*kind(iorbto_m)
  deallocate(iorbto_m,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_m',subname)

  i_all=-product(shape(iorbto_ishell))*kind(iorbto_ishell)
  deallocate(iorbto_ishell,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_ishell',subname)


  i_all=-product(shape(iorbto_iexpobeg))*kind(iorbto_iexpobeg)
  deallocate(iorbto_iexpobeg,stat=i_stat)
  call memocc(i_stat,i_all,'iorbto_iexpobeg',subname)


END SUBROUTINE createPawProjectorsArrays
!!***




!>   input guess wavefunction diagonalization
subroutine input_wf_diag(iproc,nproc,at,&
     orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
     nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
     nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU,input)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, except_this_one => input_wf_diag
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,ixc,symObj
  integer, intent(inout) :: nspin,nvirt
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(inout) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  type(GPU_pointers), intent(inout) :: GPU
  type(input_variables):: input
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
  type(gaussian_basis), intent(out) :: G !basis for davidson IG
  real(wp), dimension(:), pointer :: psi,hpsi,psit,rhocore
  real(dp), dimension(:), pointer :: pkernel,pkernelseq
  integer, intent(in) ::potshortcut
  integer, dimension(*), intent(in) :: irrzon
  real(dp), dimension(*), intent(in) :: phnons
  !local variables
  character(len=*), parameter :: subname='input_wf_diag'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0,ncplx
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eexctX,eproj_sum,etol,accurex
  type(orbitals_data) :: orbse
  type(communications_arrays) :: commse
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc,passmat
  real(gp), dimension(:), allocatable :: locrad
  type(locreg_descriptors), dimension(:), allocatable :: Llr
  real(wp), dimension(:), pointer :: pot
  real(wp), dimension(:,:,:), pointer :: psigau

  allocate(norbsc_arr(at%natsc+1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  !spin for inputguess orbitals
  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       orbs,orbse,norbsc_arr,locrad,G,psigau,eks)

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbse,commse,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbse,commse)  

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  !check the communication distribution
  !call check_communications(iproc,nproc,orbse,Glr,commse)

  !once the wavefunction coefficients are known perform a set 
  !of nonblocking send-receive operations to calculate overlap matrices

!!!  !create mpirequests array for controlling the success of the send-receive operation
!!!  allocate(mpirequests(nproc-1+ndebug),stat=i_stat)
!!!  call memocc(i_stat,mpirequests,'mpirequests',subname)
!!!
!!!  call nonblocking_transposition(iproc,nproc,G%ncoeff,orbse%isorb+orbse%norbp,&
!!!       orbse%nspinor,psigau,orbse%norb_par,mpirequests)

  !experimental part for building the localisation regions
  if (at%geocode == 'F') then
     !allocate the array of localisation regions
     allocate(Llr(at%nat+ndebug),stat=i_stat)
     !call memocc(i_stat,Llr,'Llr',subname)

     !print *,'locrad',locrad

     call determine_locreg(at%nat,rxyz,locrad,hx,hy,hz,Glr,Llr)

     do iat=1,at%nat
        call deallocate_lr(Llr(iat),subname)
!!$        call deallocate_wfd(Llr(iat)%wfd,subname)
!!$        if (Llr(iat)%geocode=='F') then
!!$           call deallocate_bounds(Llr(iat)%bounds,subname)
!!$        end if
     end do

     !i_all=-product(shape(Llr))*kind(Llr)
     deallocate(Llr,stat=i_stat) !these allocation are special
     !call memocc(i_stat,i_all,'Llr',subname)
  end if

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(psi(orbse%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  !allocate arrays for the GPU if a card is present
  switchGPUconv=.false.
  switchOCLconv=.false.
  if (GPUconv .and. potshortcut ==0 ) then
     call prepare_gpu_for_locham(Glr%d%n1,Glr%d%n2,Glr%d%n3,nspin_ig,&
          hx,hy,hz,Glr%wfd,orbse,GPU)
  else if (OCLconv .and. potshortcut ==0) then
     call allocate_data_OCL(Glr%d%n1,Glr%d%n2,Glr%d%n3,at%geocode,&
          nspin_ig,hx,hy,hz,Glr%wfd,orbse,GPU)
     if (iproc == 0) write(*,*)&
          'GPU data allocated'
  else if (GPUconv .and. potshortcut >0 ) then
     switchGPUconv=.true.
     GPUconv=.false.
  else if (OCLconv .and. potshortcut >0 ) then
     switchOCLconv=.true.
     OCLconv=.false.
  end if


  !use only the part of the arrays for building the hamiltonian matrix
  call gaussians_to_wavelets_new(iproc,nproc,Glr,orbse,hx,hy,hz,G,&
       psigau(1,1,min(orbse%isorb+1,orbse%norb)),psi)


  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  !application of the hamiltonian for gaussian based treatment
  call sumrho(iproc,nproc,orbse,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
       & Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,nspin,GPU, &
       & symObj, irrzon, phnons)
     
  !-- if spectra calculation uses a energy dependent potential
  !    input_wf_diag will write (to be used in abscalc)
  !    the density to the file electronic_density.cube
  !  The writing is activated if  5th bit of  in%potshortcut is on.
  if( iand( potshortcut,16)==0 .and. potshortcut /= 0) then
     call plot_density_cube_old(at%geocode,'electronic_density',&
          iproc,nproc,Glr%d%n1,Glr%d%n2,Glr%d%n3,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nscatterarr(iproc,2),  & 
          nspin,hxh,hyh,hzh,at,rxyz,ngatherarr,rhopot(1+nscatterarr(iproc,4)*Glr%d%n1i*Glr%d%n2i))
  endif
  !---
  
  if(orbs%nspinor==4) then
     !this wrapper can be inserted inside the poisson solver 
     call PSolverNC(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
          nscatterarr(iproc,1),& !this is n3d
          ixc,hxh,hyh,hzh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
  else
     !Allocate XC potential
     if (nscatterarr(iproc,2) >0) then
        allocate(potxc(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*nspin+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     else
        allocate(potxc(1+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     end if

     call XC_potential(at%geocode,'D',iproc,nproc,&
          Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
          rhopot,eexcu,vexcu,nspin,rhocore,potxc)


     if( iand(potshortcut,4)==0) then
        call H_potential(at%geocode,'D',iproc,nproc,&
             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
     endif


     !sum the two potentials in rhopot array
     !fill the other part, for spin, polarised
     if (nspin == 2) then
        call dcopy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),rhopot(1),1,&
             rhopot(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)+1),1)
     end if
     !spin up and down together with the XC part
     call axpy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*nspin,1.0_dp,potxc(1),1,&
          rhopot(1),1)


     i_all=-product(shape(potxc))*kind(potxc)
     deallocate(potxc,stat=i_stat)
     call memocc(i_stat,i_all,'potxc',subname)

  end if

!!!  if (nproc == 1) then
!!!     !calculate the overlap matrix as well as the kinetic overlap
!!!     !in view of complete gaussian calculation
!!!     allocate(ovrlp(G%ncoeff*G%ncoeff),stat=i_stat)
!!!     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!!     allocate(tmp(G%ncoeff,orbse%norb),stat=i_stat)
!!!     call memocc(i_stat,tmp,'tmp',subname)
!!!     allocate(smat(orbse%norb,orbse%norb),stat=i_stat)
!!!     call memocc(i_stat,smat,'smat',subname)
!!!
!!!     !overlap calculation of the gaussian matrix
!!!     call gaussian_overlap(G,G,ovrlp)
!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!
!!!     !print overlap matrices
!!!     do i=1,orbse%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!     end do
!!!
!!!     !overlap calculation of the kinetic operator
!!!     call kinetic_overlap(G,G,ovrlp)
!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!
!!!     !print overlap matrices
!!!     tt=0.0_wp
!!!     do i=1,orbse%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!        tt=tt+smat(i,i)
!!!     end do
!!!     print *,'trace',tt
!!!
!!!     !overlap calculation of the kinetic operator
!!!     call cpu_time(t0)
!!!     call potential_overlap(G,G,rhopot,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!          ovrlp)
!!!     call cpu_time(t1)
!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!
!!!     !print overlap matrices
!!!     tt=0.0_wp
!!!     do i=1,orbse%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!        tt=tt+smat(i,i)
!!!     end do
!!!     print *,'trace',tt
!!!     print *, 'time',t1-t0
!!!
!!!     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!!     deallocate(ovrlp,stat=i_stat)
!!!     call memocc(i_stat,i_all,'ovrlp',subname)
!!!     i_all=-product(shape(tmp))*kind(tmp)
!!!     deallocate(tmp,stat=i_stat)
!!!     call memocc(i_stat,i_all,'tmp',subname)
!!!     i_all=-product(shape(smat))*kind(smat)
!!!     deallocate(smat,stat=i_stat)
!!!     call memocc(i_stat,i_all,'smat',subname)
!!!  end if

  if(potshortcut>0) then
!!$    if (GPUconv) then
!!$       call free_gpu(GPU,orbs%norbp)
!!$    end if
     if (switchGPUconv) then
        GPUconv=.true.
     end if
     if (switchOCLconv) then
        OCLconv=.true.
     end if

     call deallocate_orbs(orbse,subname)
     
     !deallocate the gaussian basis descriptors
     call deallocate_gwf(G,subname)
    
     i_all=-product(shape(psigau))*kind(psigau)
     deallocate(psigau,stat=i_stat)
     call memocc(i_stat,i_all,'psigau',subname)
     call deallocate_comms(commse,subname)
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr',subname)
    return 
  end if

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(hpsi(orbse%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  !call dcopy(orbse%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp

  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,nspin,&
       orbse%norb,orbse%norbp,ngatherarr,rhopot,pot)

  call HamiltonianApplication(iproc,nproc,at,orbse,hx,hy,hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,pot,&
       psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)

  !deallocate potential
  call free_full_potential(nproc,pot,subname)

!!!  !calculate the overlap matrix knowing that the original functions are gaussian-based
!!!  allocate(thetaphi(2,G%nat+ndebug),stat=i_stat)
!!!  call memocc(i_stat,thetaphi,'thetaphi',subname)
!!!  thetaphi=0.0_gp
!!!
!!!  !calculate the scalar product between the hamiltonian and the gaussian basis
!!!  allocate(hpsigau(G%ncoeff,orbse%norbp+ndebug),stat=i_stat)
!!!  call memocc(i_stat,hpsigau,'hpsigau',subname)
!!!
!!!
!!!  call wavelets_to_gaussians(at%geocode,orbse%norbp,Glr%d%n1,Glr%d%n2,Glr%d%n3,G,&
!!!       thetaphi,hx,hy,hz,Glr%wfd,hpsi,hpsigau)
!!!
!!!  i_all=-product(shape(thetaphi))*kind(thetaphi)
!!!  deallocate(thetaphi,stat=i_stat)
!!!  call memocc(i_stat,i_all,'thetaphi',subname)

  accurex=abs(eks-ekin_sum)
  !tolerance for comparing the eigenvalues in the case of degeneracies
  etol=accurex/real(orbse%norbu,gp)
  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
  if (iproc == 0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  endif

!!!  call Gaussian_DiagHam(iproc,nproc,at%natsc,nspin,orbs,G,mpirequests,&
!!!       psigau,hpsigau,orbse,etol,norbsc_arr)


!!!  i_all=-product(shape(mpirequests))*kind(mpirequests)
!!!  deallocate(mpirequests,stat=i_stat)
!!!  call memocc(i_stat,i_all,'mpirequests',subname)

!!!  i_all=-product(shape(hpsigau))*kind(hpsigau)
!!!  deallocate(hpsigau,stat=i_stat)
!!!  call memocc(i_stat,i_all,'hpsigau',subname)

  !free GPU if it is the case
  if (GPUconv) then
     call free_gpu(GPU,orbse%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbse,nspin_ig)
  end if

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')&
       'Input Wavefunctions Orthogonalization:'

  !nullify psit (will be created in DiagHam)
  nullify(psit)

  !psivirt can be eliminated here, since it will be allocated before davidson
  !with a gaussian basis
!!$  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
!!$       psi,hpsi,psit,orbse,commse,etol,norbsc_arr,orbsv,psivirt)

  !allocate the passage matrix for transforming the LCAO wavefunctions in the IG wavefucntions
  ncplx=1
  if (orbs%nspinor > 1) ncplx=2
  allocate(passmat(ncplx*orbs%nkptsp*(orbse%norbu*orbs%norbu+orbse%norbd*orbs%norbd)+ndebug),stat=i_stat)
  call memocc(i_stat,passmat,'passmat',subname)
  !print '(a,10i5)','iproc,passmat',iproc,ncplx*orbs%nkptsp*(orbse%norbu*orbs%norbu+orbse%norbd*orbs%norbd),&
  !     orbs%nspinor,orbs%nkptsp,orbse%norbu,orbse%norbd,orbs%norbu,orbs%norbd

  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
       psi,hpsi,psit,input%orthpar,passmat,orbse,commse,etol,norbsc_arr)

  i_all=-product(shape(passmat))*kind(passmat)
  deallocate(passmat,stat=i_stat)
  call memocc(i_stat,i_all,'passmat',subname)

  if (input%itrpmax > 1 .or. input%Tel > 0.0_gp) then
     !use the eval array of orbse structure to save the original values
     allocate(orbse%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,orbse%eval,'orbse%eval',subname)
     
     call dcopy(orbs%norb*orbs%nkpts,orbs%eval(1),1,orbse%eval(1),1)

     !add a small displacement in the eigenvalues
     do iorb=1,orbs%norb*orbs%nkpts
        tt=builtin_rand(idum)
        orbs%eval(iorb)=orbs%eval(iorb)*(1.0_gp+max(input%Tel,1.0e-3_gp)*real(tt,gp))
     end do

     !correct the occupation numbers wrt fermi level
     call evaltoocc(iproc,nproc,.false.,input%Tel,orbs)

     !restore the occupation numbers
     call dcopy(orbs%norb*orbs%nkpts,orbse%eval(1),1,orbs%eval(1),1)

     i_all=-product(shape(orbse%eval))*kind(orbse%eval)
     deallocate(orbse%eval,stat=i_stat)
     call memocc(i_stat,i_all,'orbse%eval',subname)
  end if

  call deallocate_comms(commse,subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  if (iproc == 0) then
     !gaussian estimation valid only for Free BC
     if (at%geocode == 'F') then
        write(*,'(1x,a,1pe9.2)') 'expected accuracy in energy ',accurex
        write(*,'(1x,a,1pe9.2)') &
          'expected accuracy in energy per orbital ',accurex/real(orbs%norb,kind=8)
        !write(*,'(1x,a,1pe9.2)') &
        !     'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
     end if
  endif

  !here we can define the subroutine which generates the coefficients for the virtual orbitals
  call deallocate_gwf(G,subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)

  call deallocate_orbs(orbse,subname)
     
END SUBROUTINE input_wf_diag
