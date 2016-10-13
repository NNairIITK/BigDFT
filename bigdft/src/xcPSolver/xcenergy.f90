!> @file
!!  Exact-exchange routines
!! @author
!!    Copyright (C) 2002-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate the array of the core density for the atom iat
subroutine calc_rhocore_iat(iproc,atoms,ityp,rx,ry,rz,cutoff,hxh,hyh,hzh,&
     n1i,n2i,n3i,i3s,n3d,charge,rhocore) 
  use module_defs, only: gp,dp,wp
  use module_types
  use yaml_output
  use bounds, only: ext_buffers
  use numerics, only: oneo4pi => oneofourpi
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,i3s,n3d,iproc,ityp 
  real(gp), intent(in) :: rx,ry,rz,cutoff,hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  real(dp), dimension(n1i*n2i*n3d,0:9), intent(inout) :: rhocore
  real(gp),intent(out) :: charge
  !local variables
  !n(c) character(len=*), parameter :: subname='calc_rhocore'
!  real(gp), parameter :: oneo4pi=.079577471545947_wp
  logical :: gox,goy,perx,pery,perz
  integer :: ig,ngv,ngc,isx,isy,isz,iex,iey,iez
  integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,ilcc,islcc
  integer :: i1,i2,i3,j1,j2,j3,ind
  real(gp) :: x,y,z,r2,rhov,rhoc,chv,chc
  real(gp) :: charge_from_gaussians,spherical_gaussian_value
  real(gp) :: drhoc,drhov,drhodr2

  !find the correct position of the nlcc parameters
  call nlcc_start_position(ityp,atoms,ngv,ngc,islcc)

  ilcc=islcc
  chv=0.0_gp
  do ig=1,(ngv*(ngv+1))/2
     ilcc=ilcc+1
     !read(79,*)rhovxp(ig),(rhovc(ig,j),j=1,4)
     chv=chv+charge_from_gaussians(atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc))
  end do
  chv=sqrt(2.0_gp*atan(1.0_gp))*chv

  chc=0.0_gp
  do ig=1,(ngc*(ngc+1))/2
     ilcc=ilcc+1
     !read(79,*)rhocxp(ig),(rhocc(ig,j),j=1,4)
     chc=chc+charge_from_gaussians(atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc))
     !rhocc(ig,1)*rhocxp(ig)**3+3.0_gp*rhocc(ig,2)*rhocxp(ig)**5+&
     !     15.0_gp*rhocc(ig,3)*rhocxp(ig)**7+105.0_gp*rhocc(ig,4)*rhocxp(ig)**9
  end do
  chc=sqrt(2.0_gp*atan(1.0_gp))*chc

  !close(unit=79)

  if (iproc == 0) &
       & call yaml_map('NLCC, Calculate core density for atom',trim(atoms%astruct%atomnames(ityp)))

!  if (iproc == 0) call yaml_map('Analytic core charge',chc-chv,fmt='(f12.6)')
  charge=chc-chv
  !if (iproc == 0) write(*,'(1x,a,f12.6)',advance='no')' analytic core charge: ',chc-chv

  !conditions for periodicity in the three directions
  perx=(atoms%astruct%geocode /= 'F')
  pery=(atoms%astruct%geocode == 'P')
  perz=(atoms%astruct%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (n3d >0) then

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     do i3=isz,iez
        z=real(i3,kind=8)*hzh-rz
        !call ind_positions(perz,i3,n3,j3,goz)
        if (perz) then
           j3=modulo(i3,n3i)
           if (j3 -n3i >= i3s-1 ) j3=j3-n3i !to wrap with negative values
           if (j3 + n3i < i3s+n3d-1) j3=j3+n3i
        else
           j3=i3
        end if
        !in periodic case nbl3=0
        j3=j3+nbl3+1

        !if (atoms%astruct%geocode /= 'F' .and. j3 >= modulo(i3s,n3i) .and. i3s < 0) j3=j3-n3i
        if (j3 >= i3s .and. j3 <= i3s+n3d-1) then
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              !call ind_positions(pery,i2,n2,j2,goy)
              call ind_positions_new(pery,i2,n2i,j2,goy)
              if (goy) then
                 do i1=isx,iex
                    x=real(i1,kind=8)*hxh-rx
                    !call ind_positions(perx,i1,n1,j1,gox)
                    call ind_positions_new(perx,i1,n1i,j1,gox)
                    if (gox) then
                       r2=x**2+y**2+z**2
                       !here we can sum up the gaussians for the
                       !valence density and the core density
                       !restart again from the previously calculated index
                       ilcc=islcc
                       rhov=0.0_dp
                       drhov=0.0_dp
                       do ig=1,(ngv*(ngv+1))/2
                          ilcc=ilcc+1
                          rhov=rhov+&
                               spherical_gaussian_value(r2,atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc),0)
                          !derivative wrt r2
                          drhov=drhov+&
                               spherical_gaussian_value(r2,atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc),1)

                          !arg=r2/rhovxp(ig)**2
                          !(rhovc(ig,1)+r2*rhovc(ig,2)+r2**2*rhovc(ig,3)+r2**3*rhovc(ig,4))*&
                          !     exp(-0.5_gp*arg)
                       end do
                       rhoc=0.0_dp
                       drhoc=0.0_dp
                       do ig=1,(ngc*(ngc+1))/2
                          ilcc=ilcc+1
                          !arg=r2/rhocxp(ig)**2
                          rhoc=rhoc+&
                               spherical_gaussian_value(r2,atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc),0)
                          drhoc=drhoc+&
                               spherical_gaussian_value(r2,atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc),1)
                          !(rhocc(ig,1)+r2*rhocc(ig,2)+r2**2*rhocc(ig,3)+r2**3*rhocc(ig,4))*&
                          !     exp(-0.5_gp*arg)
                       end do

                       !if (j3 >= i3s .and. j3 <= i3s+n3d-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       rhocore(ind,0)=rhocore(ind,0)+oneo4pi*(rhoc-rhov)
                       drhodr2=drhoc-drhov
                       rhocore(ind,1)=rhocore(ind,1)+x*drhodr2*oneo4pi
                       rhocore(ind,2)=rhocore(ind,2)+y*drhodr2*oneo4pi
                       rhocore(ind,3)=rhocore(ind,3)+z*drhodr2*oneo4pi
!stress components
                       rhocore(ind,4)=rhocore(ind,4)+x*x*drhodr2*oneo4pi
                       rhocore(ind,5)=rhocore(ind,5)+y*y*drhodr2*oneo4pi
                       rhocore(ind,6)=rhocore(ind,6)+z*z*drhodr2*oneo4pi
                       rhocore(ind,7)=rhocore(ind,7)+y*z*drhodr2*oneo4pi
                       rhocore(ind,8)=rhocore(ind,8)+x*z*drhodr2*oneo4pi
                       rhocore(ind,9)=rhocore(ind,9)+x*y*drhodr2*oneo4pi

!!$                 !print out the result, to see what happens
!!$                 if (z==0.0_gp .and. y==0.0_gp) then
!!$                    write(16,'(3(1x,i0),10(1pe25.17))')j1+1+nbl1,j2+1+nbl2,j3,rhocore(ind),rhoc,rhov,r2,x
!!$                 end if
                    endif
                 enddo
              end if
           enddo
        end if
     enddo
  end if
  
END SUBROUTINE calc_rhocore_iat

!> Accumulate the contribution of atom iat to core density for PAW atoms.
subroutine mkcore_paw_iat(iproc,atoms,ityp,rx,ry,rz,cutoff,hxh,hyh,hzh,&
     n1i,n2i,n3i,i3s,n3d,core_mesh,rhocore, ncmax, ifftsph, rr, rcart, raux)
  use module_defs, only: dp, gp, pi_param
  use module_base, only: bigdft_mpi
  use module_dpbox, only: denspot_distribution
  use module_atoms
  use dynamic_memory
  use dictionaries, only: f_err_raise
  use f_utils

  use m_pawrad,  only : pawrad_type, pawrad_init, pawrad_free
  use m_paw_numeric, only: paw_sort_dp, paw_splint
  use bounds, only: ext_buffers

  implicit none

  !Arguments ---------------------------------------------
  !scalars
  integer, intent(in) :: n1i,n2i,n3i,i3s,n3d,iproc,ityp, ncmax
  real(gp), intent(in) :: rx,ry,rz,cutoff,hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  type(pawrad_type), intent(in) ::core_mesh
  real(dp), dimension(n1i*n2i*n3d,0:9), intent(inout) :: rhocore
  integer, dimension(ncmax), intent(inout) :: ifftsph
  real(gp), dimension(ncmax), intent(inout) :: rr, raux
  real(gp), dimension(3,ncmax), intent(inout) :: rcart

  !Local variables ------------------------------
  !scalars
  !buffer to be added at the end of the last dimension of an array to control bounds_check
  integer :: i1,i2,i3,ii,jj
  integer :: iex,iey,iez,isx,isy,isz,ind,j1,j2,j3
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  integer :: nfgd,nfgd_r0
  real(gp) :: factor, ucvol,rr2,r2shp
  real(gp) :: xx,yy,zz
  logical :: perx,pery,perz,gox,goy,goz
  integer, dimension(:), allocatable :: iperm

  ! *************************************************************************

  !mpi
  ucvol = n1i * n2i * n3i / hxh / hyh / hzh

  !Conditions for periodicity in the three directions
  perx=(atoms%astruct%geocode /= 'F')
  pery=(atoms%astruct%geocode == 'P')
  perz=(atoms%astruct%geocode /= 'F')

  !Compute values of external buffers
  call ext_buffers(perx,nbl1,nbr1) 
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  iperm = f_malloc(ncmax, id = "iperm")

  if (n3d >0) then
     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)
     
     r2shp = cutoff * cutoff
     !    
     do i3=isz,iez
        zz=real(i3,gp)*hzh-rz
        call ind_positions_new(perz,i3,n3i,j3,goz)
        j3=j3+nbl3+1
        if (j3 < i3s .or. j3 > i3s+n3d-1) cycle

        !      Initialize counters
        nfgd=0
        nfgd_r0=0
        !      
        do i2=isy,iey
           yy=real(i2,gp)*hyh-ry
           call ind_positions_new(pery,i2,n2i,j2,goy)
           !        
           do i1=isx,iex
              xx=real(i1,gp)*hxh-rx
              call ind_positions_new(perx,i1,n1i,j1,gox)
              rr2=xx**2+yy**2+zz**2
              if (goy .and. gox .and. rr2<=r2shp) then
                 ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i

                 nfgd=nfgd+1
                 rr(nfgd)=(rr2)**0.5
                 ifftsph(nfgd)=ind
                 rcart(1, nfgd) = xx
                 rcart(2, nfgd) = yy
                 rcart(3, nfgd) = zz
                 iperm(nfgd) = nfgd
              end if !j3..
           end do !i1
        end do !i2

        !      All of the following  could be done inside or outside the loops (i2,i1,i3)
        !      Outside the loops: the memory consuption increases.
        !      Inside the inner loop: the time of calculation increases.
        !      Here, I chose to do it here, somewhere in the middle.
        if(nfgd==0)      cycle

        call paw_sort_dp(nfgd, rr(1), iperm(1), 1.d-16)
        !  Evaluate spline fit of core charge density
        !  from tcoredens(:,1) and tcoredens(:,3)
        call paw_splint(atoms%pawtab(ityp)%core_mesh_size, core_mesh%rad,&
             & atoms%pawtab(ityp)%tcoredens(1,1), atoms%pawtab(ityp)%tcoredens(1,3),&
             & nfgd,rr(1),raux(1))
        !  Accumulate contributions to core density on the entire cell
        do ii=1,nfgd
           jj=ifftsph(iperm(ii))
           rhocore(jj, 0) = rhocore(jj, 0) + raux(ii)
        end do
        !  Evaluate spline fit of first derivative core charge density
        !  from tcoredens(:,2) and tcoredens(:,4)
        call paw_splint(atoms%pawtab(ityp)%core_mesh_size, core_mesh%rad,&
             & atoms%pawtab(ityp)%tcoredens(1,2), atoms%pawtab(ityp)%tcoredens(1,4),&
             & nfgd,rr(1),raux(1))
        !  Accumulate contributions to derivative of core density on the entire cell
        do ii=1,nfgd
           jj=ifftsph(iperm(ii))
           if (rr(ii) > 0._gp) then
              factor = 0.5_gp / rr(ii)
              xx = rcart(1, iperm(ii)) * factor
              yy = rcart(2, iperm(ii)) * factor
              zz = rcart(3, iperm(ii)) * factor
           else
              xx = 1._gp; yy = 1._gp; zz = 1._gp
           end if
           rhocore(jj,1) = rhocore(jj,1) + xx * raux(ii)
           rhocore(jj,2) = rhocore(jj,2) + yy * raux(ii)
           rhocore(jj,3) = rhocore(jj,3) + zz * raux(ii)
           !stress components
           rhocore(jj,4) = rhocore(jj,4) + rcart(1, iperm(ii)) * xx * raux(ii)
           rhocore(jj,5) = rhocore(jj,5) + rcart(2, iperm(ii)) * yy * raux(ii)
           rhocore(jj,6) = rhocore(jj,6) + rcart(3, iperm(ii)) * zz * raux(ii)
           rhocore(jj,7) = rhocore(jj,7) + rcart(2, iperm(ii)) * zz * raux(ii)
           rhocore(jj,8) = rhocore(jj,8) + rcart(1, iperm(ii)) * zz * raux(ii)
           rhocore(jj,9) = rhocore(jj,9) + rcart(1, iperm(ii)) * yy * raux(ii)
        end do
     end do !i3
  end if

  call f_free(iperm)

end subroutine mkcore_paw_iat

!> Given a charge density, calculates the exchange-correlation potential
!! SYNOPSIS
!!    @param nproc       number of processors
!!    @param iproc       label of the process,from 0 to nproc-1
!!    @param n01,n02,n03 global dimension in the three directions. They are the same no matter if the 
!!                datacode is in 'G' or in 'D' position.
!!    @param ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. 
!!                The XC functional codes follow the ABINIT convention.
!!    @param hgrids    grid spacings. For the isolated BC case for the moment they are supposed to 
!!                be equal in the three directions
!!    @param rho         Main input array. it represents the density values on the grid points
!!    @param potxc       Main output array, the values on the grid points of the XC potential
!!    @param exc,vxc     XC energy and integral of @f$\rho V_{xc}@f$ respectively
!! @warning
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
subroutine XC_potential(geocode,datacode,iproc,nproc,mpi_comm,n01,n02,n03,xcObj,hgrids,&
     rho,exc,vxc,nspin,rhocore,rhohat,potxc,xcstr,dvxcdrho)
  use module_base
  !Rename dp into except_dp in order to avoid conflict with the definitions provided by module_base
  use Poisson_Solver, except_dp => dp, except_gp => gp
  !Idem
  use module_interfaces, only: calc_gradient
  use module_xc
  use module_types, only: TCAT_EXCHANGECORR
  use abi_interfaces_xc_lowlevel, only: abi_mkdenpos
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
  integer, intent(in) :: iproc,nproc,n01,n02,n03,mpi_comm
  integer, intent(in) :: nspin !< Value of the spin-polarisation
  real(gp), dimension(3), intent(in) :: hgrids
  type(xc_info), intent(in) :: xcObj
  real(gp), intent(out) :: exc,vxc
  real(dp), dimension(*), intent(inout) :: rho
  real(wp), dimension(:,:,:,:), pointer :: rhocore !associated if useful
  real(wp), dimension(:,:,:,:), pointer :: rhohat !associated if useful
  real(wp), dimension(*), intent(out) :: potxc
  real(dp), dimension(6), intent(out) :: xcstr
  real(dp), dimension(:,:,:,:), target, intent(out), optional :: dvxcdrho
  !local variables
  character(len=*), parameter :: subname='XC_potential'
  logical :: wrtmsg
  !n(c) integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: ierr,i,j,i3s_fake,i3xcsh_fake
  integer :: i1,i2,i3,istart,iend,iwarn,i3start,jend,jproc
  integer :: nxc,nwbl,nwbr,nxt,nwb,nxcl,nxcr,ispin,istden,istglo
  integer :: ndvxc,order
  !real(dp),parameter:: tol14=0.00000000000001_dp
  real(dp),parameter:: tol20=0.0000000000000000000001_dp
  real(dp) :: eexcuLOC,vexcuLOC,vexcuRC
  integer, dimension(:,:), allocatable :: gather_arr
  real(dp), dimension(:), allocatable :: rho_G
  real(dp), dimension(:,:,:,:,:), allocatable :: gradient
  real(dp), dimension(:,:,:,:), allocatable :: vxci
  real(gp), dimension(:), allocatable :: energies_mpi
  real(dp), dimension(:,:,:,:), pointer :: dvxci
  real(dp), dimension(6) :: wbstr, rhocstr

  !call timing(iproc,'Exchangecorr  ','ON')
  call f_routine(id='XC_potential')
  call f_timing(TCAT_EXCHANGECORR,'ON')

  call f_zero(xcstr)
  call f_zero(wbstr)
  call f_zero(rhocstr)

  wrtmsg=.false.
  !calculate the dimensions wrt the geocode
  select case(geocode)
  case('P')
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, periodic BC, dimensions: ',n01,n02,n03,'   proc',nproc,'   ixc:',xcObj%ixc,' ... '
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
  case('S')
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, surfaces BC, dimensions: ',n01,n02,n03,'   proc',nproc,'   ixc:',xcObj%ixc,' ... '
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  case('F')
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, free  BC, dimensions: ',n01,n02,n03,'   proc',nproc,'   ixc:',xcObj%ixc,' ... '
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  case('W')
     if (iproc==0 .and. wrtmsg) &
          write(*,'(1x,a,3(i5),a,i5,a,i7,a)',advance='no')&
          'PSolver, wires  BC, dimensions: ',n01,n02,n03,'   proc',nproc,'   ixc:',xcObj%ixc,' ... '
     call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  case default
     stop 'XC potential: geometry code not admitted'
  end select

  !dimension for exchange-correlation (different in the global or distributed case)
  !let us calculate the dimension of the portion of the rho array to be passed 
  !to the xc routine
  !this portion will depend on the need of calculating the gradient or not, 
  !and whether the White-Bird correction must be inserted or not 
  !(absent only in the LB ixc=13 case)
  
  !nxc is the effective part of the third dimension that is being processed
  !nxt is the dimension of the part of rho that must be passed to the gradient routine
  !nwb is the dimension of the part of rho in the wb-postprocessing routine
  !note: nxc <= nwb <= nxt
  !the dimension are related by the values of nwbl and nwbr
  !      nxc+nxcl+nxcr-2 = nwb
  !      nwb+nwbl+nwbr = nxt
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)

  call xc_dimensions(geocode,xc_isgga(xcObj),(xcObj%ixc/=13),istart,iend,&
       & m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s_fake,i3xcsh_fake)
  nwb=nxcl+nxc+nxcr-2
  nxt=nwbr+nwb+nwbl

  !quick return if no Semilocal XC potential is required (Hartree or Hartree-Fock)
  if (any(xcObj%ixc == [XC_HARTREE,XC_HARTREE_FOCK,XC_NO_HARTREE])) then
     if (datacode == 'G') then
        call f_zero(n01*n02*n03*nspin,potxc(1))
     else
        call f_zero(n01*n02*nxc*nspin,potxc(1))
     end if
     if (nspin == 2) call axpy(n01*n02*nxc,1.d0,rho(n01*n02*nxc+1),1,rho(1),1)
     exc=0.0_gp
     vxc=0.0_gp
     call f_timing(TCAT_EXCHANGECORR,'OF')
     !call timing(iproc,'Exchangecorr  ','OF')
     call f_release_routine()
     return
  end if
  
!!$ ! !if rhocore is associated we should add it on the charge density
!!$ ! if (associated(rhocore)) then
!!$ !    if (nspin == 1) then
!!$ !       !sum the complete core density for non-spin polarised calculations
!!$ !       call axpy(m1*m3*nxt,1.0_wp,rhocore(1,1,1,1),1,rho(1),1)
!!$ !    else if (nspin==2) then
!!$ !       !for spin-polarised calculation consider half per spin index
!!$ !       call axpy(m1*m3*nxt,0.5_wp,rhocore(1,1,1,1),1,rho(1),1)
!!$ !       call axpy(m1*m3*nxt,0.5_wp,rhocore(1,1,1,1),1,rho(1+m1*m3*nxt),1)
!!$ !    end if
!!$ ! end if
 
  if (datacode=='G') then
     !starting address of rho in the case of global i/o
     i3start=istart+2-nxcl-nwbl
     if((nspin==2 .and. nproc>1) .or. i3start <=0 .or. i3start+nxt-1 > n03 ) then
        !allocation of an auxiliary array for avoiding the shift of the density
        rho_G = f_malloc(m1*m3*nxt*2,id='rho_G')
        !here we put the modulo of the results for the non-isolated GGA
        do ispin=1,nspin
           do i3=1,nxt
              do i2=1,m3
                 do i1=1,m1
                    i=i1+(i2-1)*m1+(i3-1)*m1*m3+(ispin-1)*m1*m3*nxt
                    j=i1+(i2-1)*n01+(modulo(i3start+i3-2,n03))*n01*n02+(ispin-1)*n01*n02*n03
                    rho_G(i)=rho(j)
                 end do
              end do
           end do
        end do
     end if
  else if (datacode == 'D') then
     !distributed i/o
     i3start=1
  else
     stop 'PSolver: datacode not admitted'
  end if

  !if rhohat is present, substract it from charge density
  if (associated(rhohat) .and. m1*m3*nxc > 0) then
     if (datacode=='G' .and. (i3start <=0 .or. i3start+nxc-1 > n03 )) then
        call axpy(m1*m3*nxc,-1._wp,rhohat(1,1,1,1),1,rho_G(1),1)
        if (nspin==2) call axpy(m1*m3*nxc,-1._wp,rhohat(1,1,1,2),1,&
             rho_G(1+m1*m3*nxt),1)
     else
        call axpy(m1*m3*nxc,-1._wp,rhohat(1,1,1,1),1,rho(1+n01*n02*i3xcsh_fake),1)
        if (nspin==2) call axpy(m1*m3*nxc,-1._wp,rhohat(1,1,1,2),1,&
             rho(1+n01*n02*i3xcsh_fake+m1*m3*nxt),1)
     end if
  end if

  !print *,'density must go from',min(istart+1,m2),'to',iend,'with n2/2=',n2/2
  !print *,'        it goes from',i3start+nwbl+nxcl-1,'to',i3start+nxc-1
  !print *,'istart',i3start,nxcl,nwbl,istart,iproc,geocode

  !rescale the density to apply that to ABINIT routines
  if (nspin==1) then !rho_g does not enter here
     if (datacode=='G' .and. (i3start <=0 .or. i3start+nxt-1 > n03 )) then
        call vscal(m1*m3*nxt,0.5_dp,rho_G(1),1)
     else
        call vscal(m1*m3*nxt,0.5_dp,rho(1+n01*n02*(i3start-1)),1)
     end if
  end if
  !allocate array for XC potential enlarged for the WB procedure
  vxci = f_malloc((/ m1, m3, max(1, nwb), nspin /),id='vxci')

  !allocate the array of the second derivative of the XC energy if it is needed
  !Allocations of the exchange-correlation terms, depending on the ixc value
  if (present(dvxcdrho)) then
     if (nspin==1) then 
        order=-2
     else
        order=2
     end if
     ndvxc = size(dvxcdrho, 4)
     dvxci => dvxcdrho
  else
     order=1
     ndvxc=0
     !here ndebug is not put since it is taken form dvxcdrho (not very good)
     dvxci = f_malloc_ptr((/ m1, m3, max(1, nwb), ndvxc /),id='dvxci')
  end if

  !calculate gradient
  if (xc_isgga(xcObj) .and. nxc > 0) then
     !computation of the gradient
     gradient = f_malloc((/ 1.to.m1, 1.to.m3, 1.to.nwb, 1.to.2*nspin-1, 0.to.3 /),id='gradient')

     !!the calculation of the gradient will depend on the geometry code
     !this operation will also modify the density arrangment for a GGA calculation
     !in parallel and spin-polarised, since ABINIT routines need to calculate
     !the XC terms for spin up and then spin down
     if(datacode=='G' .and. &
          ((nspin==2 .and. nproc > 1) .or. i3start <=0 .or. i3start+nxt-1 > n03 )) then
        call calc_gradient(geocode,m1,m3,nxt,nwb,nwbl,nwbr,rho_G,nspin,&
             hgrids(1),hgrids(2),hgrids(3),gradient,rhocore)
     else
        call calc_gradient(geocode,m1,m3,nxt,nwb,nwbl,nwbr,rho(1+n01*n02*(i3start-1)),nspin,&
             hgrids(1),hgrids(2),hgrids(3),gradient,rhocore)
     end if
  else
     gradient = f_malloc((/ 1, 1, 1, 1, 1 /),id='gradient')
     !add rhocore to the density
     if (associated(rhocore) .and. m1*m3*nxt > 0) then
        if (datacode=='G' .and. (i3start <=0 .or. i3start+nxt-1 > n03 )) then
           call axpy(m1*m3*nxt,0.5_wp,rhocore(1,1,1,1),1,rho_G(1),1)
           if (nspin==2) call axpy(m1*m3*nxt,0.5_wp,rhocore(1,1,1,1),1,&
                rho_G(1+m1*m3*nxt),1)
        else
           call axpy(m1*m3*nxt,0.5_wp,rhocore(1,1,1,1),1,rho(1+n01*n02*(i3start-1)),1)
           if (nspin==2) call axpy(m1*m3*nxt,0.5_wp,rhocore(1,1,1,1),1,&
                rho(1+n01*n02*(i3start-1)+m1*m3*nxt),1)
        end if
     end if
  end if


  !if (present(dvxcdrho)) then
  !   write(*,*)'Array of second derivatives of Exc allocated, dimension',ndvxc,m1,m3,nwb
  !end if
  if (istart+1 <= m2) then 
     if(datacode=='G' .and. &
          ((nspin==2 .and. nproc > 1) .or. i3start <=0 .or. i3start+nxt-1 > n03 )) then
        !allocation of an auxiliary array for avoiding the shift
        call xc_energy_new(geocode,m1,m3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
             xcObj,hgrids(1),hgrids(2),hgrids(3),rho_G,gradient,vxci,&
             eexcuLOC,vexcuLOC,order,ndvxc,dvxci,nspin,wbstr)
        !restoring the density on the original form
        do ispin=1,nspin
           do i3=1,nxt
              do i2=1,m3
                 do i1=1,m1
                    i=i1+(i2-1)*m1+(i3-1)*m1*m3+(ispin-1)*m1*m3*nxt
                    j=i1+(i2-1)*n01+(modulo(i3start+i3-2,n03))*n01*n02+(ispin-1)*n01*n02*n03
                    rho(j)=rho_G(i)
                 end do
              end do
           end do
        end do
        call f_free(rho_G)
     else
        call xc_energy_new(geocode,m1,m3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
             xcObj,hgrids(1),hgrids(2),hgrids(3),rho(1+n01*n02*(i3start-1)),gradient,vxci,&
             eexcuLOC,vexcuLOC,order,ndvxc,dvxci,nspin,wbstr)
     end if
  else
     !presumably the vxc should be initialised
     call vscal(m1*m3*nwb*nspin,0.0_dp,vxci(1,1,1,1),1)
     eexcuLOC=0.0_dp
     vexcuLOC=0.0_dp
  end if

  !deallocate gradient here
  call f_free(gradient)


  !the value of the shift depends of the distributed i/o or not
  if ((datacode=='G' .and. nproc == 1) .or. datacode == 'D') then
     !copy the relevant part of vxci on the output potxc
     call vcopy(m1*m3*nxc,vxci(1,1,nxcl,1),1,potxc(1),1)
     if (nspin == 2) then
        call vcopy(m1*m3*nxc,vxci(1,1,nxcl,2),1,potxc(1+m1*m3*nxc),1)
     end if
  end if
 
  !if (iproc == 0) print *,'n03,nxt,nxc,geocode,datacode',n03,nxt,nxc,geocode,datacode

  !recollect the final data, and build the total charge density
  !no spin index anymore
  if (datacode == 'G') then
     do i3=nxc,1,-1
        do i2=1,n02
           do i1=1,n01
              rho(i1+(i2-1)*n01+(i3+istart-1)*n01*n02)=&
                   rho(i1+(i2-1)*n01+modulo(i3-1-1+i3start,n03)*n01*n02)
           end do
        end do
     end do
  end if

  !if rhohat is present, add it to charge density
  if(associated(rhohat) .and. m1*m3*nxc .gt.0) then
     call axpy(m1*m3*nxc,1._wp,rhohat(1,1,1,1),1,rho(1),1)
     if (nspin==2) call axpy(m1*m3*nxc,1._wp,rhohat(1,1,1,2),1,&
          rho(1+m1*m3*nxt),1)
     !This will add V_xc(n) nhat to V_xc(n) n, to get: V_xc(n) (n+nhat)
     ! Only if usexcnhat /= 0, not the default case.
     !call add_to_vexcu(rhohat)
  end if

  !if rhocore is associated we then remove it from the charge density
  !and subtract its contribution from the evaluation of the XC potential integral vexcu
  if (associated(rhocore) .and. m1*m3*nxc .gt.0) then
     !at this stage the density is not anymore spin-polarised
     !sum the complete core density for non-spin polarised calculations
     call axpy(m1*m3*nxc,-1.0_wp,rhocore(1,1,i3xcsh_fake+1,1),1,rho(1),1)
     call substract_from_vexcu(rhocore)
!print *,' aaaa', vexcuRC,vexcuLOC,eexcuLOC
  end if

  !gathering the data to obtain the distribution array
  !evaluating the total ehartree,eexcu,vexcu
  if (nproc > 1) then
     !allocate(energies_mpi(2+ndebug),stat=i_stat)
     !call memocc(i_stat,energies_mpi,'energies_mpi',subname)
     energies_mpi = f_malloc(4,id='energies_mpi')

     energies_mpi(1)=eexcuLOC
     energies_mpi(2)=vexcuLOC
     call mpiallred(energies_mpi(1), 2,MPI_SUM,comm=mpi_comm,recvbuf=energies_mpi(3))
     exc=energies_mpi(3)
     vxc=energies_mpi(4)

!XC-stress term
  if (geocode == 'P') then

        if (associated(rhocore)) then
        call calc_rhocstr(rhocstr,nxc,nxt,m1,m3,i3xcsh_fake,nspin,potxc,rhocore)
        call mpiallred(rhocstr,MPI_SUM,comm=mpi_comm)
        rhocstr=rhocstr/real(n01*n02*n03,dp)
        end if

     xcstr(1:3)=(exc-vxc)/real(n01*n02*n03,dp)/product(hgrids)!hx/hy/hz
     call mpiallred(wbstr,MPI_SUM,comm=mpi_comm)
     wbstr=wbstr/real(n01*n02*n03,dp)
     xcstr(:)=xcstr(:)+wbstr(:)+rhocstr(:)
  end if
     !i_all=-product(shape(energies_mpi))*kind(energies_mpi)
     !deallocate(energies_mpi,stat=i_stat)
     !call memocc(i_stat,i_all,'energies_mpi',subname)
     call f_free(energies_mpi)  

     if (datacode == 'G') then
        !building the array of the data to be sent from each process
        !and the array of the displacement
        if (present(dvxcdrho)) then
           write(*,*)'ERROR: gathering of 2nd derivative of Exc not implemented yet!'
           stop
        end if

        gather_arr = f_malloc((/ 0.to.nproc-1, 1.to.2 /),id='gather_arr')
        do jproc=0,nproc-1
           istart=min(jproc*(md2/nproc),m2-1)
           jend=max(min(md2/nproc,m2-md2/nproc*jproc),0)
           gather_arr(jproc,1)=m1*m3*jend
           gather_arr(jproc,2)=m1*m3*istart
           !print *,'TOINSPECT',iproc,jproc,istart,jend,istart+jend,m2
        end do

        !gather all the results in the same rho array
        istart=min(iproc*(md2/nproc),m2-1)

        istden=1+n01*n02*istart
        istglo=1
        do ispin=1,nspin
           if (ispin==2) then
              istden=istden+n01*n02*n03
              istglo=istglo+n01*n02*n03
           end if
           call MPI_ALLGATHERV(vxci(1,1,nxcl,ispin),gather_arr(iproc,1),mpidtypw,&
                potxc(istglo),gather_arr(0,1),gather_arr(0,2),mpidtypw,&
                mpi_comm,ierr)
        end do

        call f_free(gather_arr)

     end if

  else
     exc=real(eexcuLOC,gp)
     vxc=real(vexcuLOC,gp)

!XC-stress term
   if (geocode == 'P') then
        if (associated(rhocore)) then
        call calc_rhocstr(rhocstr,nxc,nxt,m1,m3,i3xcsh_fake,nspin,potxc,rhocore)
        rhocstr=rhocstr/real(n01*n02*n03,dp)
        end if
    xcstr(1:3)=(exc-vxc)/real(n01*n02*n03,dp)/product(hgrids)!hx/hy/hz
    wbstr=wbstr/real(n01*n02*n03,dp)
    xcstr(:)=xcstr(:)+wbstr(:)+rhocstr(:)
   end if

  end if

  call f_free(vxci)

  if (.not.present(dvxcdrho)) then
     call f_free_ptr(dvxci)
  end if

  call f_timing(TCAT_EXCHANGECORR,'OF')
  call f_release_routine()
  !call timing(iproc,'Exchangecorr  ','OF')
  !if (iproc==0 .and. wrtmsg) write(*,'(a)')'done.'

contains
subroutine substract_from_vexcu(rhoin)
 implicit none
 real(wp),dimension(:,:,:,:),intent(in)::rhoin
 
 vexcuRC=0.0_gp
 do i3=1,nxc
    do i2=1,m3
       do i1=1,m1
          !do i=1,nxc*m3*m1
          i=i1+(i2-1)*m1+(i3-1)*m1*m3
          vexcuRC=vexcuRC+rhoin(i1,i2,i3+i3xcsh_fake,1)*potxc(i)
       end do
    end do
 end do
 if (nspin==2) then
    do i3=1,nxc
       do i2=1,m3
          do i1=1,m1
             !do i=1,nxc*m3*m1
             !vexcuRC=vexcuRC+rhocore(i+m1*m3*i3xcsh_fake)*potxc(i+m1*m3*nxc)
             i=i1+(i2-1)*m1+(i3-1)*m1*m3
             vexcuRC=vexcuRC+rhoin(i1,i2,i3+i3xcsh_fake,1)*potxc(i+m1*m3*nxc)
          end do
       end do
    end do
    !divide the results per two because of the spin multiplicity
    vexcuRC=0.5*vexcuRC
 end if
 vexcuRC=vexcuRC*product(hgrids)!real(hx*hy*hz,gp)
 !subtract this value from the vexcu
 vexcuLOC=vexcuLOC-vexcuRC

end subroutine substract_from_vexcu

!This routine is slightly different than 'substract_from_vexcu',
! since rhoin has spin up and spin down components (there is no 0.5 factor).
subroutine add_to_vexcu(rhoin)
 implicit none
 real(wp),dimension(:,:,:,:),intent(in)::rhoin

 vexcuRC=0.0_gp
 do i3=1,nxc
    do i2=1,m3
       do i1=1,m1
          !do i=1,nxc*m3*m1
          i=i1+(i2-1)*m1+(i3-1)*m1*m3
          vexcuRC=vexcuRC+rhoin(i1,i2,i3+i3xcsh_fake,1)*potxc(i)
       end do
    end do
 end do
 if (nspin==2) then
    do i3=1,nxc
       do i2=1,m3
          do i1=1,m1
             !do i=1,nxc*m3*m1
             !vexcuRC=vexcuRC+rhocore(i+m1*m3*i3xcsh_fake)*potxc(i+m1*m3*nxc)
             i=i1+(i2-1)*m1+(i3-1)*m1*m3
             vexcuRC=vexcuRC+rhoin(i1,i2,i3+i3xcsh_fake,2)*potxc(i+m1*m3*nxc)
          end do
       end do
    end do
 end if
 vexcuRC=vexcuRC*product(hgrids)!real(hx*hy*hz,gp)
 !add this value to the vexcu
 vexcuLOC=vexcuLOC+vexcuRC

end subroutine add_to_vexcu


END SUBROUTINE XC_potential


!> Calculate the XC terms from the given density in a distributed way.
!! it assign also the proper part of the density to the zf array 
!! which will be used for the core of the FFT procedure.
!! Following the values of ixc and of sumpion, the array pot_ion is either summed or assigned
!! to the XC potential, or even ignored.
!!
!! @warning
!!    The dimensions of pot_ion must be compatible with geocode, datacode and ixc.
!!    Since the arguments of these routines are indicated with the *,
!!    it is IMPERATIVE to refer to PSolver routine for the correct allocation sizes.
subroutine xc_energy_new(geocode,m1,m3,nxc,nwb,nxt,nwbl,nwbr,&
     nxcl,nxcr,xc,hx,hy,hz,rho,gradient,vxci,exc,vxc,order,ndvxc,dvxci,nspden,wbstr)

  use module_base
  use module_xc

  implicit none

  !Arguments
  !> Indicates the boundary conditions (BC) of the problem:
  !!   'F' free BC, isolated systems.
  !!       The program calculates the solution as if the given density is
  !!       "alone" in R^3 space.
  !!   'S' surface BC, isolated in y direction, periodic in xz plane                
  !!       The given density is supposed to be periodic in the xz plane,
  !!       so the dimensions in these direction mus be compatible with the FFT
  !!       Beware of the fact that the isolated direction is y!
  !!   'P' periodic BC.
  !!       The density is supposed to be periodic in all the three directions,
  !!       then all the dimensions must be compatible with the FFT.
  !!       No need for setting up the kernel.
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: m1,m3     !< Global dimensions in the three directions.
  integer, intent(in) :: nxc       !< Value of the effective distributed dimension in the third direction
  integer, intent(in) :: nwb       !< Enlarged dimension for calculating the WB correction
  integer, intent(in) :: nxt       !< Enlarged dimension for calculating the GGA case 
                                   !! (further enlarged for compatibility with WB correction if it is the case)
  integer, intent(in) :: nwbl,nwbr !< nwb=nxc+nxcl+nxcr-2, nwb+nwbl+nwbr=nxt.
  integer, intent(in) :: nxcl,nxcr !< Shifts in the three directions to be compatible with the relation
  !> eXchange-Correlation code. Indicates the XC functional to be used 
  !!   for calculating XC energies and potential. 
  !!   ixc=0 indicates that no XC terms are computed. 
  !!   The XC functional codes follow the ABINIT convention or if negative the libXC one.
  type(xc_info), intent(in) :: xc
  integer, intent(in) :: order,ndvxc,nspden
  real(gp), intent(in) :: hx,hy,hz                            !< Grid spacings. 
  real(dp), dimension(*), intent(in) :: gradient              !< of size 1 if not needed
  real(dp), dimension(m1,m3,nxt,nspden), intent(inout) :: rho !< Density in the distributed format, also in spin-polarised
  real(dp), dimension(m1,m3,nwb,nspden), intent(out) :: vxci
  real(dp), dimension(m1,m3,nwb,ndvxc), intent(out) :: dvxci
  real(dp), intent(out) :: exc,vxc                            !< XC energy and integral of @f$\rho V_{xc}@f$ respectively
  real(dp), dimension(6), intent(inout) :: wbstr

  !Local variables----------------
  character(len=*), parameter :: subname='xc_energy'
  real(dp), dimension(:,:,:), allocatable :: exci
  real(dp), dimension(:,:,:,:), allocatable :: dvxcdgr
  !real(dp), dimension(:,:,:,:,:), allocatable :: gradient
  real(dp) :: elocal,vlocal,rhov,sfactor
  integer :: npts,offset,ispden
  integer :: i1,i2,i3,j1,j2,j3,jp2,jppp2
  logical :: use_gradient

  call f_routine(id='xc_energy_new')

  !check for the dimensions
  if (nwb/=nxcl+nxc+nxcr-2 .or. nxt/=nwbr+nwb+nwbl) then
     print *,'the XC dimensions are not correct'
     print *,'nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr',nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
     stop
  end if

  !starting point of the density array for the GGA cases in parallel
  offset=nwbl+1
  !divide by two the density to applicate it in the ABINIT xc routines
  use_gradient = xc_isgga(xc)

  if (use_gradient) then
!!$     !computation of the gradient
!!$     allocate(gradient(m1,m3,nwb,2*nspden-1,0:3+ndebug),stat=i_stat)
!!$     call memocc(i_stat,gradient,'gradient',subname)
!!$
!!$     !!the calculation of the gradient will depend on the geometry code
!!$     !this operation will also modify the density arrangment for a GGA calculation
!!$     !in parallel and spin-polarised, since ABINIT routines need to calculate
!!$     !the XC terms for spin up and then spin down
!!$     call calc_gradient(geocode,m1,m3,nxt,nwb,nwbl,nwbr,rho,nspden,&
!!$          real(hx,dp),real(hy,dp),real(hz,dp),gradient)

     dvxcdgr = f_malloc((/ m1, m3, nwb, 3 /),id='dvxcdgr')
  else
!!$     allocate(gradient(1,1,1,1,1+ndebug),stat=i_stat)
!!$     call memocc(i_stat,gradient,'gradient',subname)
     dvxcdgr = f_malloc((/ 1, 1, 1, 1 /),id='dvxcdgr')
  end if
  
  !Allocations
  exci = f_malloc((/ m1, m3, nwb /),id='exci')

  !this part can be commented out if you don't want to use ABINIT modules
  !of course it must be substituted with an alternative XC calculation
  npts=m1*m3*nwb

  !do a separate calculation of the grid to allow for OMP parallelisation
  ! Do the calculation.
  if (abs(order) == 1) then
     call xc_getvxc(xc, npts,exci,nspden,rho(1,1,offset,1),vxci,gradient,dvxcdgr)
!!MM
!     do i3=1,nxt
!        write(*,*) nspden, rho(m1/2,m3/2,i3,nspden)
!     end do
!!MM
  else if (abs(order) == 2) then
     call xc_getvxc(xc, npts,exci,nspden,rho(1,1,offset,1),vxci,gradient,dvxcdgr,dvxci)
  end if
  wbstr(:)=0._dp
  if (use_gradient) then
     ! Do not calculate the White-Bird term in the Leeuwen Baerends XC case
     if (xc%ixc /= 13 .and. xc%ixc /= -160) then
        call vxcpostprocessing(geocode,m1,m3,nwb,nxc,nxcl,nxcr,nspden,3,gradient,&
             real(hx,dp),real(hy,dp),real(hz,dp),dvxcdgr,vxci,wbstr)
     end if
!print *,wbstr
     !restore the density array in the good position if it was shifted for the parallel GGA
     !operation not necessarily needed, but related to the fact that the array has three
     !indices which make it difficult to treat
     !one should convert the operations with one indices arrays
     if (nspden==2 .and. nxt /= nwb) then
        j3=nwb+1
        do i3=nwb-nwbr,1,-1
           j3=j3-1
           do i2=1,m3
              do i1=1,m1
                 rho(i1,i2,nwbl+j3,2)=rho(i1,i2,i3,2)
              end do
           end do
        end do
        do i3=nxt,nwb+nwbl+1,-1 !we have nwbr points
           j3=j3-1
           do i2=1,m3
              do i1=1,m1
                 rho(i1,i2,nwbl+j3,2)=rho(i1,i2,i3,1)
              end do
           end do
        end do
     end if
  end if
  !end of the part that can be commented out

  if (allocated(dvxcdgr)) then
     call f_free(dvxcdgr)
  end if
!!$  if (allocated(gradient)) then
!!$     i_all=-product(shape(gradient))*kind(gradient)
!!$     deallocate(gradient,stat=i_stat)
!!$     call memocc(i_stat,i_all,'gradient',subname)
!!$  end if
  !     rewind(300)
  !     do ispden=1,nspden
  !        do i3=1,nxt
  !           do i2=1,m3
  !              do i1=1,m1
  !                 write(300,'(f18.12)') rho(i1,i2,i3,ispden)
  !              end do
  !           end do
  !        end do
  !     end do

  !this part should be put out from this routine due to the Global distribution code
  exc=0.0_dp
  vxc=0.0_dp
  sfactor=1.0_dp
  if(nspden==1) sfactor=2.0_dp

  !compact the rho array into the total charge density
  !try to use dot and vcopy routines, more general
  ! e.g. exc=dot(m1*m3*nxc,exci(1,1,nxcl),1,rho(1,1,offset+nxcl-1,ispden),1)

  ispden=1
  do jp2=1,nxc
     j2=offset+jp2+nxcl-2
     jppp2=jp2+nxcl-1
     do j3=1,m3
        do j1=1,m1
           rhov=rho(j1,j3,j2,ispden)
           elocal=exci(j1,j3,jppp2)
           vlocal=vxci(j1,j3,jppp2,ispden)
           exc=exc+elocal*rhov
           vxc=vxc+vlocal*rhov
           rho(j1,j3,jp2,1)=sfactor*rhov!restore the original normalization
           !potxc(j1,j3,jp2,ispden)=real(vlocal,wp)
        end do
     end do
  end do
  !spin-polarised case
  if (nspden==2) then
     ispden=2
     do jp2=1,nxc
        j2=offset+jp2+nxcl-2
        jppp2=jp2+nxcl-1
        do j3=1,m3
           do j1=1,m1
              rhov=rho(j1,j3,j2,ispden)
              elocal=exci(j1,j3,jppp2)
              vlocal=vxci(j1,j3,jppp2,ispden)
              exc=exc+elocal*rhov
              vxc=vxc+vlocal*rhov
              rho(j1,j3,jp2,1)=rho(j1,j3,jp2,1)+sfactor*rhov
              !potxc(j1,j3,jp2,ispden)=real(vlocal,dp)
           end do
        end do
     end do
  end if

  !the two factor is due to the 
  !need of using the density of states in abinit routines
  exc=sfactor*real(hx*hy*hz,dp)*exc
  vxc=sfactor*real(hx*hy*hz,dp)*vxc

  !De-allocations
  call f_free(exci)
!  call MPI_BARRIER(bigdft_mpi%mpi_comm,i_stat)
!stop

  call f_release_routine()

END SUBROUTINE xc_energy_new


!>    Calculate the XC terms from the given density in a distributed way.
!!    it assign also the proper part of the density to the zf array 
!!    which will be used for the core of the FFT procedure.
!!    Following the values of ixc and of sumpion, the array pot_ion is either summed or assigned
!!    to the XC potential, or even ignored.
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    m1,m2,m3    global dimensions in the three directions.
!!    md1,md2,md3 dimensions of the arrays compatible with the FFT in the three directions.
!!    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
!!                for calculating XC energies and potential. 
!!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
!!                the ABINIT convention.
!!    hx,hy,hz    grid spacings. 
!!    rhopot      density in the distributed format.
!!    karray      kernel of the poisson equation. It is provided in distributed case, with
!!                dimensions that are related to the output of the PS_dim4allocation routine
!!                it MUST be created by following the same geocode as the Poisson Solver.
!!    pot_ion     additional external potential that is added to the output, 
!!                when the XC parameter ixc/=0. It is always provided in the distributed form,
!!                clearly without the overlapping terms which are needed only for the XC part
!!    exc,vxc     XC energy and integral of @f$\rho V_{xc}@f$ respectively
!!    nxc         value of the effective distributed dimension in the third direction
!!    nwb         enlarged dimension for calculating the WB correction
!!    nxt         enlarged dimension for calculating the GGA case 
!!                (further enlarged for compatibility with WB correction if it is the case)
!!    nwbl,nwbr
!!    nxcl,nxcr   shifts in the three directions to be compatible with the relation
!!                nxc+nxcl+nxcr-2=nwb, nwb+nwbl+nwbr=nxt.
!!    sumpion     logical value which states whether to sum pot_ion to the final result or not
!!                if sumpion==.true. zfionxc will be pot_ion+vxci
!!                if sumpion==.false. zfionxc will be vxci
!!                this value is ignored when ixc=0. In that case zfionxc is untouched
!!    zf          output array corresponding to the density which can be passed to FFT part
!!    zfionxc     output array which will contain pot_ion+vxci or vxci, following sumpion
!!
!! @warning
!!    The dimensions of pot_ion must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to refer to PSolver routine for the correct allocation sizes.
subroutine xc_energy(geocode,m1,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
     nxcl,nxcr,xc,hx,hy,hz,rhopot,pot_ion,sumpion,zf,zfionxc,exc,vxc,nproc,nspden)

  use module_base
  use module_xc
  use abi_interfaces_xc_lowlevel, only: abi_drivexc,abi_size_dvxc
  use module_interfaces, only: calc_gradient

  implicit none

  !Arguments----------------------
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: sumpion
  integer, intent(in) :: m1,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,nproc,nspden
  integer, intent(in) :: nwbl,nwbr
  real(gp), intent(in) :: hx,hy,hz
  type(xc_info), intent(in) :: xc
  real(dp), dimension(m1,m3,nxt,nspden), intent(inout) :: rhopot
  real(wp), dimension(*), intent(in) :: pot_ion
  real(dp), dimension(md1,md3,md2/nproc), intent(out) :: zf
  real(wp), dimension(md1,md3,md2/nproc,nspden), intent(out) :: zfionxc
  real(dp), intent(out) :: exc,vxc

  !Local variables----------------
  character(len=*), parameter :: subname='xc_energy'
  real(dp), dimension(:,:,:), allocatable :: exci,d2vxci
  real(dp), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(dp), dimension(:,:,:,:,:), allocatable :: gradient
  real(dp) :: elocal,vlocal,rho,potion,sfactor
  integer :: npts,order,offset,ispden
  integer :: i1,i2,i3,j1,j2,j3,jp2,jpp2,jppp2
  integer :: ndvxc,nvxcdgr,ngr2,nd2vxc
  logical :: use_gradient
  real(dp),dimension(6) :: wbstr
  real(dp), dimension(:,:,:,:), pointer :: rhocore_fake
  !check for the dimensions
  if (nwb/=nxcl+nxc+nxcr-2 .or. nxt/=nwbr+nwb+nwbl) then
     print *,'the XC dimensions are not correct'
     print *,'nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr',nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
     stop
  end if

  nullify(rhocore_fake)

  !these are always the same
  order=1
  
  !starting point of the density array for the GGA cases in parallel
  offset=nwbl+1
  if (xc%ixc/=0) then
     !divide by two the density to applicate it in the ABINIT xc routines
     if(nspden==1) then
        do i3=1,nxt
           do i2=1,m3
              do i1=1,m1
                 rhopot(i1,i2,i3,nspden)=.5_dp*rhopot(i1,i2,i3,nspden)
              end do
           end do
        end do
     end if
!     rewind(301)
!     do ispden=1,nspden
!        do i3=1,nxt
!           do i2=1,m3
!              do i1=1,m1
!                 write(301,'(f18.12)') rhopot(i1,i2,i3,ispden)
!              end do
!           end do
!        end do
!     end do
     use_gradient = xc_isgga(xc)

     !Allocations of the exchange-correlation terms, depending on the ixc value
     nd2vxc=1
     call abi_size_dvxc(xc%ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)

     if (use_gradient) then
        !computation of the gradient
        gradient = f_malloc((/ 1.to.m1, 1.to.m3, 1.to.nwb, 1.to.2*nspden-1, 0.to.3 /),id='gradient')

        !!the calculation of the gradient will depend on the geometry code
        !if (geocode=='F') then
           nullify(rhocore_fake)
           !!call calc_gradient(geocode,m1,m3,nxt,nwb,nwbl,nwbr,rhopot,nspden,&
           !!     real(hx,dp),real(hy,dp),real(hz,dp),gradient)
           call calc_gradient(geocode,m1,m3,nxt,nwb,nwbl,nwbr,rhopot,nspden,&
                real(hx,dp),real(hy,dp),real(hz,dp),gradient,rhocore_fake)
        !else
        !print *,'geocode=',geocode,&
        !     ':the calculation of the gradient is still to be performed in this case'
        !stop
        !end if

     end if

     !Allocations
     exci = f_malloc((/ m1, m3, nwb /),id='exci')
     vxci = f_malloc((/ m1, m3, nwb, nspden /),id='vxci')

     if (ndvxc/=0) then
        dvxci = f_malloc((/ m1, m3, nwb, ndvxc /),id='dvxci')
     end if
     if (nvxcdgr/=0) then
        dvxcdgr = f_malloc((/ m1, m3, nwb, nvxcdgr /),id='dvxcdgr')
     end if
     if ((xc%ixc==3 .or. (xc%ixc>=7 .and. xc%ixc<=15)) .and. order==3) then
        d2vxci = f_malloc((/ m1, m3, nwb /),id='d2vxci')
     end if

     if (.not.allocated(gradient) .and. nxc/=nxt ) then
        print *,'xc_energy: if nxt/=nxc the gradient must be allocated'
        stop
     end if

     !this part can be commented out if you don't want to use ABINIT modules
     !of course it must be substituted with an alternative XC calculation
     npts=m1*m3*nwb
     !let us apply ABINIT routines
     !case with gradient
     if (xc%ixc >= 11 .and. xc%ixc <= 16) then
        if (order**2 <= 1 .or. xc%ixc == 16) then
           if (xc%ixc /= 13) then             
              call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                   &grho2_updn=gradient,vxcgr=dvxcdgr) 
           else
              call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                   &grho2_updn=gradient) 
           end if
        else if (order /= 3) then
           if (xc%ixc /= 13) then             
              call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                   &dvxc=dvxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
           else
              call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                   &dvxc=dvxci,grho2_updn=gradient) 
           end if
        else if (order == 3) then
           if (xc%ixc /= 13) then             
              call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                   &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
           else
              call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                   &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient) 
           end if
        end if

        !cases without gradient
     else if (xc%ixc >= 0) then
        if (order**2 <=1 .or. xc%ixc >= 31 .and. xc%ixc<=34) then
           call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr)
        else if (order==3 .and. (xc%ixc==3 .or. xc%ixc>=7 .and. xc%ixc<=10)) then
           call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                &dvxc=dvxci,d2vxc=d2vxci)
        else
           call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                &dvxc=dvxci)
        end if
        !case with libXC, with and without gradient
     else if (xc%ixc < 0) then
        call abi_drivexc(exci,xc%ixc,npts,nspden,order,rhopot(1,1,offset,1),vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,     &
             &      grho2_updn=gradient,vxcgr=dvxcdgr)
     end if

     if (use_gradient) then
        !do not calculate the White-Bird term in the Leeuwen Baerends XC case
        if (xc%ixc/=13) then
           call vxcpostprocessing(geocode,m1,m3,nwb,nxc,nxcl,nxcr,nspden,nvxcdgr,gradient,&
                real(hx,dp),real(hy,dp),real(hz,dp),dvxcdgr,vxci,wbstr)
        end if

        !restore the density array in the good position if it was shifted for the parallel GGA
        !operation not necessarily needed
        if (nspden==2 .and. nxt /= nwb) then
           j3=nwb+1
           do i3=nwb-nwbr,1,-1
              j3=j3-1
              do i2=1,m3
                 do i1=1,m1
                    rhopot(i1,i2,nwbl+j3,2)=rhopot(i1,i2,i3,2)
                 end do
              end do
           end do
           do i3=nxt,nwb+nwbl+1,-1 !we have nwbr points
              j3=j3-1
              do i2=1,m3
                 do i1=1,m1
                    rhopot(i1,i2,nwbl+j3,2)=rhopot(i1,i2,i3,1)
                 end do
              end do
           end do
        end if
     end if
     !end of the part that can be commented out

     if (allocated(dvxci)) then
        call f_free(dvxci)
     end if
     if (allocated(dvxcdgr)) then
        call f_free(dvxcdgr)
     end if
     if (allocated(d2vxci)) then
        call f_free(d2vxci)
     end if
     if (allocated(gradient)) then
        call f_free(gradient)
     end if

!     rewind(300)
!     do ispden=1,nspden
!        do i3=1,nxt
!           do i2=1,m3
!              do i1=1,m1
!                 write(300,'(f18.12)') rhopot(i1,i2,i3,ispden)
!              end do
!           end do
!        end do
!     end do


     exc=0.0_dp
     vxc=0.0_dp
     sfactor=1.0_dp
     if(nspden==1) sfactor=2.0_dp
     if (sumpion) then
        !summing the xc potential into the zfionxc array with pot_ion
        ispden=1
        do jp2=1,nxc
           j2=offset+jp2+nxcl-2
           jppp2=jp2+nxcl-1
           do j3=1,m3
              do j1=1,m1
                 jpp2=j1+(j3-1)*m1+(jp2-1)*m1*m3
                 rho=rhopot(j1,j3,j2,ispden)
                 potion=pot_ion(jpp2)
                 elocal=exci(j1,j3,jppp2)
                 vlocal=vxci(j1,j3,jppp2,ispden)
                 exc=exc+elocal*rho
                 vxc=vxc+vlocal*rho
                 zf(j1,j3,jp2)=sfactor*rho !restore the original normalization
                 zfionxc(j1,j3,jp2,ispden)=potion+real(vlocal,wp)
              end do
              do j1=m1+1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.0_wp
              end do
           end do
           do j3=m3+1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.0_wp
              end do
           end do
        end do
        do jp2=nxc+1,md2/nproc
           do j3=1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.d0
                 zfionxc(j1,j3,jp2,ispden)=0.0_wp
              end do
           end do
        end do
        !spin-polarised case
        if (nspden==2) then
           ispden=2
           do jp2=1,nxc
              j2=offset+jp2+nxcl-2
              jppp2=jp2+nxcl-1
              do j3=1,m3
                 do j1=1,m1
                    jpp2=j1+(j3-1)*m1+(jp2-1)*m1*m3
                    rho=rhopot(j1,j3,j2,ispden)
                    potion=pot_ion(jpp2)
                    elocal=exci(j1,j3,jppp2)
                    vlocal=vxci(j1,j3,jppp2,ispden)
                    exc=exc+elocal*rho
                    vxc=vxc+vlocal*rho
                    zf(j1,j3,jp2)=zf(j1,j3,jp2)+sfactor*rho !restore original normalization
                    zfionxc(j1,j3,jp2,ispden)=potion+real(vlocal,wp)
                 end do
                 do j1=m1+1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.0_wp
                 end do
              end do
              do j3=m3+1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.0_wp
                 end do
              end do
           end do
           do jp2=nxc+1,md2/nproc
              do j3=1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.0_wp
                 end do
              end do
           end do
        end if
     else
        !the zfionxc aray contain only the XC potential. pot_ion is ignored
        ispden=1
        do jp2=1,nxc
           j2=offset+jp2+nxcl-2
           jppp2=jp2+nxcl-1
           do j3=1,m3
              do j1=1,m1
                 rho=rhopot(j1,j3,j2,ispden)
                 elocal=exci(j1,j3,jppp2)
                 vlocal=vxci(j1,j3,jppp2,ispden)
                 exc=exc+elocal*rho
                 vxc=vxc+vlocal*rho
                 zf(j1,j3,jp2)=sfactor*rho!restore the original normalization
                 zfionxc(j1,j3,jp2,ispden)=real(vlocal,wp)
              end do
              do j1=m1+1,md1
                 zf(j1,j3,jp2)=0.0_dp
                 zfionxc(j1,j3,jp2,ispden)=0.0_wp
              end do
           end do
           do j3=m3+1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.0_dp
                 zfionxc(j1,j3,jp2,ispden)=0.0_wp
              end do
           end do
        end do
        do jp2=nxc+1,md2/nproc
           do j3=1,md3
              do j1=1,md1
                 zf(j1,j3,jp2)=0.0_dp
                 zfionxc(j1,j3,jp2,ispden)=0.0_wp
              end do
           end do
        end do
        !spin-polarised case
        if (nspden==2) then
           ispden=2
           do jp2=1,nxc
              j2=offset+jp2+nxcl-2
              jppp2=jp2+nxcl-1
              do j3=1,m3
                 do j1=1,m1
                    rho=rhopot(j1,j3,j2,ispden)
                    elocal=exci(j1,j3,jppp2)
                    vlocal=vxci(j1,j3,jppp2,ispden)
                    exc=exc+elocal*rho
                    vxc=vxc+vlocal*rho
                    zf(j1,j3,jp2)=zf(j1,j3,jp2)+sfactor*rho!restore the original normalization
                    zfionxc(j1,j3,jp2,ispden)=real(vlocal,dp)
                 end do
                 do j1=m1+1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.0_wp
                 end do
              end do
              do j3=m3+1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.0_wp
                 end do
              end do
           end do
           do jp2=nxc+1,md2/nproc
              do j3=1,md3
                 do j1=1,md1
                    zfionxc(j1,j3,jp2,ispden)=0.0_wp
                 end do
              end do
           end do
        end if
     end if
     !the two factor is due to the 
     !need of using the density of states in abinit routines
     exc=sfactor*real(hx*hy*hz,dp)*exc
     vxc=sfactor*real(hx*hy*hz,dp)*vxc

     !De-allocations
     call f_free(exci)
     call f_free(vxci)


  else

     !case without XC terms
     !distributing the density in the zf array
     exc=0.0_dp
     vxc=0.0_dp
     do jp2=1,nxc
        j2=offset+jp2+nxcl-2
        do j3=1,m3
           do j1=1,m1
              zf(j1,j3,jp2)=rhopot(j1,j3,j2,1)
           end do
           do j1=m1+1,md1
              zf(j1,j3,jp2)=0.0_dp
           end do
        end do
        do j3=m3+1,md3
           do j1=1,md1
              zf(j1,j3,jp2)=0.0_dp
           end do
        end do
     end do
     do jp2=nxc+1,md2/nproc
        do j3=1,md3
           do j1=1,md1
              zf(j1,j3,jp2)=0.0_dp
           end do
        end do
     end do

  end if

END SUBROUTINE xc_energy


!> Correct the XC potential with the White-Bird formula, to be used for the 
!! GGA case. Works either in parallel of in serial, by proper change of the 
!! arguments.
subroutine vxcpostprocessing(geocode,n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr,&
gradient,hx,hy,hz,dvxcdgr,wb_vxc,wbstr)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr
  real(dp), intent(in) :: hx,hy,hz
  real(dp), dimension(n01,n02,n03,2*nspden-1,0:3), intent(in) :: gradient
  real(dp), dimension(n01,n02,n03,nvxcdgr), intent(in) :: dvxcdgr
  real(dp), dimension(n01,n02,n03,nspden), intent(inout) :: wb_vxc
  real(dp), dimension(6),intent(inout) :: wbstr
  !Local variables
  character(len=*), parameter :: subname='vxcpostprocessing'
  integer :: i1,i2,i3,dir_i
  real(dp) :: dnexcdgog,grad_i,rho_up,rho_down,rho_tot
  real(dp), dimension(:,:,:,:,:), allocatable :: f_i

  call f_routine(id='vxcpostprocessing')

  !Body

  f_i = f_malloc((/ n01, n02, n03, 3, nspden /),id='f_i')
  
  !let us first treat the case nspden=1
  if (nspden == 1) then
     !Let us construct the object we have to manipulate with another gradient
     if (nvxcdgr == 3) then
        do dir_i=1,3
           !Let us construct the object we have to manipulate with another gradient
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    dnexcdgog=0.5_dp*dvxcdgr(i1,i2,i3,1) + dvxcdgr(i1,i2,i3,3)
                    grad_i=2.0_dp*gradient(i1,i2,i3,1,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=dnexcdgog*grad_i
                 end do
              end do
           end do
        end do
     else
        do dir_i=1,3
           !Let us construct the object we have to manipulate with another gradient
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    dnexcdgog=0.5_dp*dvxcdgr(i1,i2,i3,1)
                    grad_i=2.0_dp*gradient(i1,i2,i3,1,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=dnexcdgog*grad_i
                 end do
              end do
           end do
        end do
     end if

  !then the spin-polarized case
  else

     if (nvxcdgr == 3) then
        do dir_i=1,3
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    rho_up=gradient(i1,i2,i3,1,dir_i)  !rho_ instead of grad_ for ABINIT comp.
                    rho_down=gradient(i1,i2,i3,2,dir_i)
                    rho_tot=gradient(i1,i2,i3,3,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=rho_up*dvxcdgr(i1,i2,i3,1)+&
                         rho_tot*dvxcdgr(i1,i2,i3,3)
                    f_i(i1,i2,i3,dir_i,2)=rho_down*dvxcdgr(i1,i2,i3,2)+&
                         rho_tot*dvxcdgr(i1,i2,i3,3)
                 end do
              end do
           end do
        end do
     else
        do dir_i=1,3
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    rho_up=gradient(i1,i2,i3,1,dir_i)
                    rho_down=gradient(i1,i2,i3,2,dir_i)
                    f_i(i1,i2,i3,dir_i,1)=rho_up*dvxcdgr(i1,i2,i3,1)
                    f_i(i1,i2,i3,dir_i,2)=rho_down*dvxcdgr(i1,i2,i3,2)
                 end do
              end do
           end do
        end do
     end if

  !end of spin-polarized if statement
  end if

! wb-stress
  if (geocode == 'P') call wb_stress(n01,n02,n03,n3eff,wbl,nspden,f_i,gradient,wbstr)

  !let us now calculate the gradient and correct the result
  call wb_correction(geocode,n01,n02,n03,n3eff,wbl,wbr,f_i,hx,hy,hz,nspden,wb_vxc)


  call f_free(f_i)

  call f_release_routine()

END SUBROUTINE vxcpostprocessing

subroutine wb_stress(n1,n2,n3,n3eff,wbl,nsp,f_i,gradient,wbstr)
 use module_base
 implicit none
  integer, intent(in) :: n1,n2,n3,nsp,n3eff,wbl
  real(dp), dimension(n1,n2,n3,3,nsp) :: f_i
  real(dp), dimension(n1,n2,n3,2*nsp-1,0:3), intent(in) :: gradient
 real(dp),dimension(6),intent(out) :: wbstr
 integer :: i1,i2,i3,isp

wbstr=0._dp
!seq: 11 22 33 23 13 12
do isp=1,nsp
           do i3=wbl,n3eff+wbl-1
              do i2=1,n2
                 do i1=1,n1
        wbstr(1) = wbstr(1)-gradient(i1,i2,i3,isp,1)*f_i(i1,i2,i3,1,isp)
        wbstr(2) = wbstr(2)-gradient(i1,i2,i3,isp,2)*f_i(i1,i2,i3,2,isp)
        wbstr(3) = wbstr(3)-gradient(i1,i2,i3,isp,3)*f_i(i1,i2,i3,3,isp)
        wbstr(6) = wbstr(6)-gradient(i1,i2,i3,isp,1)*f_i(i1,i2,i3,2,isp)
        wbstr(5) = wbstr(5)-gradient(i1,i2,i3,isp,1)*f_i(i1,i2,i3,3,isp)
        wbstr(4) = wbstr(4)-gradient(i1,i2,i3,isp,2)*f_i(i1,i2,i3,3,isp)
                 end do
              end do
           end do
end do
! unpol case: up+dn=2*up
if (nsp==1) wbstr=wbstr*2._gp

!wbstress=wbstress/real(n1*n2*n3,gp)

!write(*,*) 'WB-correction to stress'
!write(*,*) wbstress(1:3)
!write(*,*) wbstress(4:6)

end subroutine wb_stress

subroutine calc_rhocstr(rhocstr,nxc,nxt,m1,m3,i3xcsh_fake,nspin,potxc,rhocore)
 use module_base
 implicit none
 real(dp),dimension(6),intent(out) :: rhocstr
 integer,intent(in) :: nxc,nxt,m1,m3,i3xcsh_fake,nspin
 real(wp), dimension(m1,m3,nxt,10), intent(in) :: rhocore
 real(wp), dimension(m1*m3*nxc*nspin), intent(in) :: potxc
 integer :: i,i1,i2,i3,isp

rhocstr=0._dp
!seq: 11 22 33 23 13 12

do isp=0,nspin-1
     do i3=1,nxc
        do i2=1,m3
           do i1=1,m1
              i=i1+(i2-1)*m1+(i3-1)*m1*m3
rhocstr(1)=rhocstr(1)+rhocore(i1,i2,i3+i3xcsh_fake,5)*potxc(i+m1*m3*nxc*isp)
rhocstr(2)=rhocstr(2)+rhocore(i1,i2,i3+i3xcsh_fake,6)*potxc(i+m1*m3*nxc*isp)
rhocstr(3)=rhocstr(3)+rhocore(i1,i2,i3+i3xcsh_fake,7)*potxc(i+m1*m3*nxc*isp)
rhocstr(4)=rhocstr(4)+rhocore(i1,i2,i3+i3xcsh_fake,8)*potxc(i+m1*m3*nxc*isp)
rhocstr(5)=rhocstr(5)+rhocore(i1,i2,i3+i3xcsh_fake,9)*potxc(i+m1*m3*nxc*isp)
rhocstr(6)=rhocstr(6)+rhocore(i1,i2,i3+i3xcsh_fake,10)*potxc(i+m1*m3*nxc*isp)
           end do
        end do
     end do
end do

if (nspin==1) rhocstr(:)=rhocstr(:)*2._gp

end subroutine calc_rhocstr
