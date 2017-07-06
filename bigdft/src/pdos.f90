!> @file
!!  Partial DOS analysis routines
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
subroutine spatially_resolved_dos(ob,output_dir)
  use module_base
  use orbitalbasis
  use locreg_operations
  use box, only: box_iterator
  implicit none
  character(len=*), intent(in) :: output_dir
  type(orbital_basis), intent(inout) :: ob
  !real(wp), dimension(ob%orbs%npsidim_orbs), intent(in) :: hpsi
  !local variables
  logical :: assert_cubic
  integer :: ispinor
  type(workarr_sumrho) :: w
  type(ket) :: psi_it
  type(box_iterator) :: box_it
  real(wp), dimension(:,:), allocatable :: psir,epsx,epsy,epsz!,hpsir
  real(wp), dimension(:), pointer :: psi_ptr!,hpsi_ptr

  if (ob%orbs%nspinor > 2) &
       call f_err_throw('Spatially resolved DoS not available for spinors',&
       err_name='BIGDFT_INPUT_VARIABLES_ERROR')

  call f_routine(id='spatially_resolved_dos')

  !this is only menaningful for the cubic code
  assert_cubic=.true.
  psi_it=orbital_basis_iterator(ob)
  loop_glr: do while(ket_next_locreg(psi_it))
     call f_assert(assert_cubic,&
          'Spatially resolved density of states only valid with nlr=1')
     call initialize_work_arrays_sumrho(psi_it%lr,.true.,w)
     !real space basis function, per orbital components
     psir = f_malloc(1.to.psi_it%ket_dims,id='psir')
     !hpsir = f_malloc(1.to.psi_it%ket_dims,id='hpsir')
     epsx=f_malloc([psi_it%lr%d%n1i,ob%orbs%norbp],id='epsx')
     epsy=f_malloc([psi_it%lr%d%n2i,ob%orbs%norbp],id='epsy')
     epsz=f_malloc([psi_it%lr%d%n3i,ob%orbs%norbp],id='epsz')
     do while(ket_next(psi_it,ilr=psi_it%ilr))
        do ispinor=1,ob%orbs%nspinor
           psi_ptr=>ob_subket_ptr(psi_it,ispinor)
           !hpsi_ptr=>ob_ket_map(hpsi,psi_it,ispinor)
           call daub_to_isf(psi_it%lr,w,psi_ptr,psir(1,ispinor))
           !call daub_to_isf(psi_it%lr,w,hpsi_ptr,hpsir(1,ispinor))
        end do
        call calculate_sdos(psi_it%nspinor,psi_it%lr%bit,psir,&!hpsir,&
             epsx(1,psi_it%iorbp),epsy(1,psi_it%iorbp),epsz(1,psi_it%iorbp))
     end do
     !deallocations of work arrays
     call f_free(psir)
     !call f_free(hpsir)
     call deallocate_work_arrays_sumrho(w)
     assert_cubic=.false. !we should exit now as we need the iterator
     box_it=psi_it%lr%bit
     exit loop_glr
  end do loop_glr
  !this should guarantee that the parallel arrays have been allocated
  if (assert_cubic) then
     epsx=f_malloc([1,1],id='epsx')
     epsy=f_malloc([1,1],id='epsy')
     epsz=f_malloc([1,1],id='epsz')
     box_it=ob%td%Glr%bit
  end if

  call write_sdos(box_it,ob%orbs%norbp,ob%orbs%norb*ob%orbs%nkpts,&
       epsx,epsy,epsz,output_dir)

  call orbital_basis_release(ob)

  call f_free(epsx)
  call f_free(epsy)
  call f_free(epsz)

  call f_release_routine()

end subroutine spatially_resolved_dos

!> write the real part of the spatially resolved density of states in the arrays epsx,y,z
subroutine calculate_sdos(nspinor,boxit,psi,&!hpsi,&
     epsx,epsy,epsz)
  use box
  use module_defs
  use f_utils, only: f_zero
  implicit none
  type(box_iterator) :: boxit
  integer, intent(in) :: nspinor
  real(wp), dimension(boxit%mesh%ndim,nspinor), intent(in) :: psi!,hpsi
  real(wp), dimension(boxit%mesh%ndims(1)), intent(out) :: epsx
  real(wp), dimension(boxit%mesh%ndims(2)), intent(out) :: epsy
  real(wp), dimension(boxit%mesh%ndims(3)), intent(out) :: epsz
  !local variables
  real(wp) :: tt
  
  !here we have to initialize the iterator over the total box
  call f_zero(epsx)
  call f_zero(epsy)
  call f_zero(epsz)

  do while(box_next_z(boxit))
     do while(box_next_y(boxit))
        do while(box_next_x(boxit))
           tt=psi(boxit%ind,1)**2!hpsi(boxit%ind,1)
           if (nspinor==2) tt=tt+psi(boxit%ind,2)**2!hpsi(boxit%ind,2)
           epsz(boxit%k)=epsz(boxit%k)+tt
           epsy(boxit%j)=epsy(boxit%j)+tt
           epsx(boxit%i)=epsx(boxit%i)+tt
        end do
     end do
     !print *,boxit%i,boxit%j,boxit%k,boxit%rxyz
  end do

  !print *,'sums',sum(epsx),sum(epsy),sum(epsz)

end subroutine calculate_sdos

subroutine write_sdos(bit,norbp,norb,epsx,epsy,epsz,output_dir)
  use module_base
  use box
  use yaml_output
  implicit none
  type(box_iterator) :: bit
  character(len=*), intent(in) :: output_dir
  integer, intent(in) :: norbp,norb
  real(wp), dimension(bit%mesh%ndims(1),norbp), intent(in) :: epsx
  real(wp), dimension(bit%mesh%ndims(2),norbp), intent(in) :: epsy
  real(wp), dimension(bit%mesh%ndims(3),norbp), intent(in) :: epsz
  !local variables
  integer :: unt
  real(wp), dimension(:), pointer :: epsx_tot,epsy_tot,epsz_tot

  !here we may gather the results and write the file
  epsx_tot=>mpigathered(epsx,comm=bigdft_mpi%mpi_comm)
  epsy_tot=>mpigathered(epsy,comm=bigdft_mpi%mpi_comm)
  epsz_tot=>mpigathered(epsz,comm=bigdft_mpi%mpi_comm)

  unt=82
  if (bigdft_mpi%iproc==0) then
     call yaml_sequence_open('SDos files',advance='no')
     call yaml_comment('Domain directions resolutions (a,b,c)')
     call f_open_file(unt,output_dir//'sdos_x.dat')
     do while(box_next_x(bit))
        write(unt,'(1pg26.16e3)',advance='no')bit%rxyz(1)
        call dump_sdos_line(unt,bit%i,bit%mesh%ndims(1),norb,epsx_tot)
     end do
     call f_close(unt)
     call yaml_sequence(output_dir//'sdos_x.dat')

     call f_open_file(unt,output_dir//'sdos_y.dat')
     do while(box_next_y(bit))
        write(unt,'(1pg26.16e3)',advance='no')bit%rxyz(2)
        call dump_sdos_line(unt,bit%j,bit%mesh%ndims(2),norb,epsy_tot)
     end do
     call f_close(unt)
     call yaml_sequence(output_dir//'sdos_y.dat')

     call f_open_file(unt,output_dir//'sdos_z.dat')
     do while(box_next_z(bit))
        write(unt,'(1pg26.16e3)',advance='no')bit%rxyz(3)
        call dump_sdos_line(unt,bit%k,bit%mesh%ndims(3),norb,epsz_tot)
     end do
     call f_close(unt)
     call yaml_sequence(output_dir//'sdos_z.dat')
     
     call yaml_sequence_close()

  end if

  call f_free_ptr(epsx_tot)
  call f_free_ptr(epsy_tot)
  call f_free_ptr(epsz_tot)

end subroutine write_sdos

subroutine dump_sdos_line(unt,i,ni,norb,epsi)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: unt,i,ni,norb
  real(wp), dimension(ni,norb), intent(in) :: epsi
  !local variables
  integer :: iorb
  do iorb=1,norb-1
     write(unt,'(1pg26.16e3)',advance='no')epsi(i,iorb)
  end do
  write(unt,'(1pg26.16e3)')epsi(i,norb)

end subroutine dump_sdos_line


!> Perform all the projection associated to local variables
subroutine local_analysis(iproc,nproc,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt)
   use module_base
   use module_types
   use module_interfaces, only: gaussian_pswf_basis
   use gaussians, only: deallocate_gwf
   use locregs
   implicit none
   integer, intent(in) :: iproc,nproc
   real(gp), intent(in) :: hx,hy,hz
   type(locreg_descriptors), intent(in) :: lr
   type(orbitals_data), intent(in) :: orbs,orbsv
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), dimension(:), pointer :: psi,psivirt
   !local variables
   character(len=*), parameter :: subname='local_analysis'
   integer :: norbpv
   !type(input_variables) :: inc
   !type(atoms_data) :: atc
   type(gaussian_basis) :: G
   real(wp), dimension(:,:), allocatable :: allpsigau,dualcoeffs
   real(gp), dimension(:,:), allocatable :: thetaphi
   real(gp), dimension(:), pointer :: Gocc
   !real(gp), dimension(:,:), pointer :: cxyz

   !the number of virtual orbitals in parallel is known only if orbsv%norb>0
   if (orbsv%norb >0) then
      norbpv=orbsv%norbp
   else
      norbpv=0
   end if

   !define the local basis starting from the input files
   !this is done to allow the calculations of charges also in points which
   !are different from the atoms.
   !NOTE: this means that the MCPA can be done only on SP calculations
   !call read_input_variables(iproc,'posinp','input.dft','','','',inc,atc,cxyz)
   !call read_system_variables('input.occup',iproc,inc,atc,radii_cf_fake,nelec,&
   !     norb,norbu,norbd,iunit)

   !shift the positions with the same value of the original positions
   !  do iat=1,atc%nat
   !     cxyz(1,iat)=cxyz(1,iat)-shift(1)
   !     cxyz(2,iat)=cxyz(2,iat)-shift(2)
   !     cxyz(3,iat)=cxyz(3,iat)-shift(3)
   !  end do

   nullify(G%rxyz)

   !extract the gaussian basis from the pseudowavefunctions
   !call gaussian_pswf_basis(31,.false.,iproc,inc%nspin,atc,cxyz,G,Gocc)
   call gaussian_pswf_basis(31,.false.,iproc,orbs%nspin,at,rxyz,G,Gocc)

   thetaphi = f_malloc0((/ 2, G%nat /),id='thetaphi')
   !call to_zero(2*G%nat,thetaphi)
   allpsigau = f_malloc((/ G%ncoeff*orbs%nspinor, orbs%norbp+norbpv /),id='allpsigau')
!print *,'there'
   !this routine should be simplified like gaussians_to_wavelets
   call wavelets_to_gaussians(lr%geocode,orbs%norbp,orbs%nspinor,&
        lr%d%n1,lr%d%n2,lr%d%n3,G,thetaphi,hx,hy,hz,lr%wfd,psi,allpsigau)
!print *,'here'
   !the same can be done for virtual orbitals if orbsv%norb > 0
   if (orbsv%norb > 0) then
      call wavelets_to_gaussians(lr%geocode,norbpv,orbsv%nspinor,&
           lr%d%n1,lr%d%n2,lr%d%n3,G,thetaphi,hx,hy,hz,lr%wfd,psivirt,&
           allpsigau(1,orbs%norbp+min(1,norbpv)))
   end if
   !calculate dual coefficients
   dualcoeffs = f_malloc(src=allpsigau,id='dualcoeffs')
   !if (G%ncoeff*orbs%nspinor*(orbs%norbp+norbpv)>0) then
   !    call vcopy(G%ncoeff*orbs%nspinor*(orbs%norbp+norbpv),allpsigau(1,1),1,dualcoeffs(1,1),1)
   !end if
   !build dual coefficients
   call dual_gaussian_coefficients(orbs%nspinor*(orbs%norbp+norbpv),G,dualcoeffs)


   !here we can calculate the Mulliken charge population
   !for any of the elements of the basis, ordered by angular momentum
   !do that only for the occupied orbitals
   call mulliken_charge_population(iproc,nproc,orbs,Gocc,G,allpsigau,dualcoeffs)

   !also partial density of states can be analysed here
   call gaussian_pdos(iproc,nproc,orbs,G,allpsigau,dualcoeffs) !n(m)

   call deallocate_gwf(G)
   nullify(G%rxyz)

   call f_free(allpsigau)
   call f_free(dualcoeffs)


   call f_free(thetaphi)

   !deallocate the auxiliary structures for the calculations
   !call deallocate_atoms(atc,subname) 
   !call free_input_variables(inc)
   call f_free_ptr(Gocc)

END SUBROUTINE local_analysis


!> Calculate Mulliken charge population
subroutine mulliken_charge_population(iproc,nproc,orbs,Gocc,G,coeff,duals)
  use module_base
  use module_types
  use gaussians
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(gaussian_basis), intent(in) :: G
  real(gp), dimension(G%ncoeff), intent(in) :: Gocc
  real(wp), dimension(G%ncoeff,orbs%nspinor,orbs%norbp), intent(in) :: coeff,duals
  !local variables
  character(len=*), parameter :: subname='mulliken_charge_population'
  character(len=11) :: shname
  integer :: icoeff,ishell,iexpo,iat,l,ng,iorb,isat,m,ispin,ig,nchannels
  integer :: ispinor,i
  real(wp) :: msum,rad,radnorm,r,sumch,mnrm
  real(wp), dimension(2) :: msumiat
  real(wp), dimension(3) :: mi,mtot
  real(wp), dimension(:,:), allocatable :: mchg,magn
 
  !allocate both for spins up and down
  mchg = f_malloc((/ G%ncoeff, 2 /),id='mchg')

  magn = f_malloc((/ G%ncoeff, 3 /),id='magn')


  !for any of the orbitals calculate the Mulliken charge
  do icoeff=1,G%ncoeff
     mchg(icoeff,1)=0.0_wp
     mchg(icoeff,2)=0.0_wp
 
     magn(icoeff,1)=0.0_wp
     magn(icoeff,2)=0.0_wp
     magn(icoeff,3)=0.0_wp

     !print '(a,100(1pe12.5))','icoeff,iorb',coeff(icoeff,:)
     !print '(a,100(1pe12.5))','idualc,iorb',duals(icoeff,:)
     do iorb=1,orbs%norbp
        if (orbs%spinsgn(orbs%isorb+iorb) == 1.0_gp .and. orbs%nspinor /= 4) then
           ispin=1
        else if (orbs%nspinor /= 4) then
           ispin=2
        end if
        !reduce the charge on site
        sumch=0.0_gp
        do ispinor=1,orbs%nspinor
           sumch=sumch+coeff(icoeff,ispinor,iorb)*duals(icoeff,ispinor,iorb)
        end do
        !reduce the magnetisation
        if (orbs%nspinor == 4) then
           mi(1)=coeff(icoeff,1,iorb)*duals(icoeff,3,iorb)+coeff(icoeff,3,iorb)*duals(icoeff,1,iorb)+&
                (coeff(icoeff,4,iorb)*duals(icoeff,2,iorb)+coeff(icoeff,2,iorb)*duals(icoeff,4,iorb))
           mi(2)=coeff(icoeff,2,iorb)*duals(icoeff,3,iorb)+coeff(icoeff,3,iorb)*duals(icoeff,2,iorb)-&
                (coeff(icoeff,1,iorb)*duals(icoeff,4,iorb)+coeff(icoeff,4,iorb)*duals(icoeff,1,iorb))
           mi(3)=coeff(icoeff,1,iorb)*duals(icoeff,1,iorb)+coeff(icoeff,2,iorb)*duals(icoeff,2,iorb)-&
                (coeff(icoeff,3,iorb)*duals(icoeff,3,iorb)+coeff(icoeff,4,iorb)*duals(icoeff,4,iorb))
        else
           mi(1)=0.0_wp
           mi(2)=0.0_wp
           mi(3)=orbs%spinsgn(orbs%isorb+iorb)*sumch          
        end if
        do i=1,3
           magn(icoeff,i)=magn(icoeff,i)+&
             orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*mi(i)
        end do
        if (orbs%nspinor /= 4) then
           mchg(icoeff,ispin)=mchg(icoeff,ispin)+&
                orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*sumch
        else
           !here the mchg represent the majority and minority spins respectively
           !modulus of m
           mnrm=nrm2(3,mi(1),1)
           !majority
           mchg(icoeff,1)=mchg(icoeff,1)+&
                orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*0.5_gp*(sumch+mnrm)
           !minority
           mchg(icoeff,2)=mchg(icoeff,2)+&
                orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*0.5_gp*(sumch-mnrm)
        end if

        !if no spin polarisation equals up and down spin quantities
     end do
     if (orbs%nspin ==1) then
        mchg(icoeff,1)=0.5_wp*mchg(icoeff,1)
        mchg(icoeff,2)=mchg(icoeff,1)
     end if
  end do

  !reduce the results
  if (nproc > 1) then
     call mpiallred(mchg,MPI_SUM,comm=bigdft_mpi%mpi_comm)
     call mpiallred(magn,MPI_SUM,comm=bigdft_mpi%mpi_comm)
  end if

  if (iproc == 0) then
     call yaml_comment('Mulliken Charge Population Analysis',hfill='-')
     call yaml_sequence_open('Mulliken Charge Population Analysis')
     call yaml_newline()
  end if

!  do iorb=1,orbs%norbp  
!     msum=0.0_wp
!     do icoeff=1,G%ncoeff
!        msum=msum+coeff(icoeff,iorb)*duals(icoeff,iorb)
!     end do
!     print *,'total sum,iorb',iorb,msum,&
!          orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)
!  end do

  !print the results as a function of the shell
  ishell=0
  iexpo=1
  icoeff=1
  msum=0.0_wp
  mtot(1)=0.0_wp
  mtot(2)=0.0_wp
  mtot(3)=0.0_wp
  do iat=1,G%nat
     if (iproc ==0) then
        call yaml_newline()
        call yaml_comment('Atom No.'//adjustl(trim(yaml_toa(iat,fmt='(i4.4)'))))
        call yaml_sequence(advance='no')
        !     call yaml_mapping_open()!flow=.true.)
     end if
     msumiat(1)=0.0_wp
     msumiat(2)=0.0_wp
     mi(1)=0.0_wp
     mi(2)=0.0_wp
     mi(3)=0.0_wp
     nchannels=0
     sumch=0.0_gp
     do isat=1,G%nshell(iat)
        ishell=ishell+1
        ng=G%ndoc(ishell)
        l=G%nam(ishell)
        !calculate mean radius (a.u.)
        rad=0.0_wp
        radnorm=0.0_wp
        do ig=1,ng
           r=G%xp(1,iexpo)
           rad=rad+(G%psiat(1,iexpo))**2*r
           radnorm=radnorm+(G%psiat(1,iexpo))**2
           iexpo=iexpo+1
        end do
        rad=rad/radnorm
        do m=1,2*l-1
           call shell_name(l,m,shname)
           msumiat(1)=msumiat(1)+mchg(icoeff,1)
           msumiat(2)=msumiat(2)+mchg(icoeff,2)
           if (iproc == 0) then
              call yaml_mapping_open(trim(shname))!, flow=.true.)
              if (orbs%nspinor /= 4) then
                    call yaml_map('Rad',rad,fmt='(f8.5)')
                    call yaml_map('Chg (up,down)', mchg(icoeff,1:2), fmt='(f8.5)')
                    call yaml_map('Partial Chg',sum(mchg(icoeff,1:2)),fmt='(f8.5)')
                    call yaml_map('Mag Pol',mchg(icoeff,1)-mchg(icoeff,2),fmt='(f8.5)')
                    call yaml_map('Net Chg',Gocc(icoeff)-(mchg(icoeff,1)+mchg(icoeff,2)),fmt='(f8.5)')
              else
                 mi(1)=mi(1)+magn(icoeff,1)
                 mi(2)=mi(2)+magn(icoeff,2)
                 mi(3)=mi(3)+magn(icoeff,3)
                    call yaml_map('Rad',rad,fmt='(f8.5)')
                    call yaml_map('Chg (Maj,Min)', mchg(icoeff,1:2), fmt='(f8.5)')
                    call yaml_map('Partial Chg',sum(mchg(icoeff,1:2)),fmt='(f8.5)')
                    call yaml_map('Mag Comp',magn(icoeff,1:3),fmt='(f8.5)')
                    call yaml_map('Net Chg',Gocc(icoeff)-(mchg(icoeff,1)+mchg(icoeff,2)),fmt='(f8.5)')
              end if
              call yaml_mapping_close()
           end if
           sumch=sumch+Gocc(icoeff)
           icoeff=icoeff+1
           nchannels=nchannels+1
        end do
     end do
     if (iproc == 0) then
        if (orbs%nspinor == 4) then
           call yaml_comment('Chg (Maj)| Chg (Min)  |Partial Chg| Mag Comp |  Net Chg')
        end if
        if (orbs%nspinor /= 4) then
           call yaml_map('Center Quantities', &
                & (/ msumiat(1),msumiat(2), msumiat(1)+msumiat(2), msumiat(1)-msumiat(2), &
                & sumch-(msumiat(1)+msumiat(2)) /), fmt='(f7.4)')
        else
           mtot(1)=mtot(1)+mi(1)
           mtot(2)=mtot(2)+mi(2)
           mtot(3)=mtot(3)+mi(3)
           call yaml_map('Center Quantities', &
                & (/ msumiat(1),msumiat(2), msumiat(1)+msumiat(2), mi(1), &
                & sumch-(msumiat(1)+msumiat(2)), &
                & mi(2), mi(3) /), fmt='(f7.4)')
        end if
     end if
     msum=msum+msumiat(1)+msumiat(2)
  end do

  if (iproc==0)call yaml_sequence_close()

  if (iproc == 0) then
     call yaml_map('Total Charge considered on the centers',msum,fmt='(f21.12)')
     if (orbs%nspinor==4) call yaml_map('Projected Magnetic density orientation',mtot,fmt='(f9.5)')
     !write(*,'(8x,a,f21.12)')'    Total Charge considered on the centers: ',msum
     !if (orbs%nspinor==4) write(*,'(7x,a,3(f9.5))')'    Projected Magnetic density orientation: ',mtot
     
  end if
  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  call f_free(mchg)

  call f_free(magn)

END SUBROUTINE mulliken_charge_population


subroutine gaussian_pdos(iproc,nproc,orbs,G,coeff,duals) !n(c) Gocc (arg:4)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(gaussian_basis), intent(in) :: G
   !n(c) real(gp), dimension(G%ncoeff), intent(in) :: Gocc
   real(wp), dimension(G%ncoeff,orbs%norbp), intent(in) :: coeff,duals
   !local variables
   character(len=*), parameter :: subname='gaussian_pdos'
   integer :: icoeff,ierr,iorb,ikpt !n(c) ispin
   integer :: jproc!,nspin
   real(wp) :: rsum,tnorm
   integer, dimension(:), allocatable :: norb_displ
   real(wp), dimension(:), allocatable :: work
   real(wp), dimension(:,:), allocatable :: pdos


   !allocate both for spins up and down
   pdos = f_malloc((/ G%ncoeff+1, orbs%norb*orbs%nkpts /),id='pdos')

   !for any of the orbitals calculate the Mulliken charge
!   nspin=1
   do icoeff=1,G%ncoeff
      do iorb=1,orbs%norbp
         !useful only for finding the spins
!         if (orbs%spinsgn(orbs%isorb+iorb) == 1.0_gp) then
            !n(c) ispin=1
!         else
!            nspin=2
            !n(c) ispin=2
!         end if
         pdos(icoeff,orbs%isorb+iorb)=coeff(icoeff,iorb)*duals(icoeff,iorb)
      end do
   end do
!   if (iproc==0) write(*,*) 'ey :', orbs%spinsgn(orbs%isorb+iorb), nspin

   !gather the results to the root process
   if (nproc > 1) then
      norb_displ = f_malloc(0.to.nproc-1,id='norb_displ')

      work = f_malloc(max((G%ncoeff+1)*orbs%norb_par(iproc, 0), 1),id='work')

      call vcopy((G%ncoeff+1)*orbs%norb_par(iproc,0),pdos(1,min(orbs%isorb+1,orbs%norb*orbs%nkpts)),1,&
           work(1),1)

      norb_displ(0)=0
      do jproc=1,nproc-1
         norb_displ(jproc)=norb_displ(jproc-1)+orbs%norb_par(jproc-1,0)
      end do

      call MPI_GATHERV(work(1),(G%ncoeff+1)*orbs%norb_par(iproc,0),mpidtypw,&
         &   pdos(1,1),(G%ncoeff+1)*orbs%norb_par(:,0),(G%ncoeff+1)*norb_displ,mpidtypw,&
         &   0,bigdft_mpi%mpi_comm,ierr)

      call f_free(work)
      call f_free(norb_displ)
   end if

   !now the results have to be written
   if (iproc == 0) then
      !renormalize the density of states to 10 (such as to gain a digit)
      tnorm=5.0_wp*real(orbs%nspin,wp)
      do iorb=1,orbs%norb*orbs%nkpts
         rsum=0.0_wp
         do icoeff=1,G%ncoeff
            rsum=rsum+pdos(icoeff,iorb)
            pdos(icoeff,iorb)=pdos(icoeff,iorb)*tnorm/real(G%ncoeff,wp)
         end do
         pdos(G%ncoeff+1,iorb)=tnorm-rsum*tnorm/real(G%ncoeff,wp)
      end do

      !first spin up, then spin down
      if (orbs%nspin == 2) then
         open(unit=12,file='pdos-up.dat',status='unknown')
      else
         open(unit=12,file='pdos.dat',status='unknown')
      end if
     write(12,'(a,a13,5x,i6,a)')  & 
          '# band', ' energy (eV),  ',G%ncoeff,' partial densities of states ' 
     do ikpt=1,orbs%nkpts
        do iorb=1,orbs%norbu
           write(12,'(i5,es14.5,5x,1000es14.5)')iorb,orbs%eval(iorb+(ikpt-1)*orbs%norb)*Ha_eV,&
                pdos(1:G%ncoeff,iorb+(ikpt-1)*orbs%norb)
        end do
        close(unit=12)
        if (orbs%norbd /= 0) then
           open(unit=12,file='pdos-down.dat',status='unknown')
           write(12,'(a,a13,5x,i6,a)')  & 
                '# band', ' energy (eV),  ',G%ncoeff,' partial densities of states ' 
           do iorb=orbs%norbu+1,orbs%norbu+orbs%norbd
              write(12,'(i5,es14.5,5x,1000es14.5)')iorb-orbs%norbu,&
                   orbs%eval(iorb+(ikpt-1)*orbs%norb)*Ha_eV,pdos(1:G%ncoeff+1,iorb+(ikpt-1)*orbs%norb)
           end do
        end if
     end do
   end if

   call f_free(pdos)

END SUBROUTINE gaussian_pdos


subroutine shell_name(l,m,name)
   implicit none
   integer, intent(in) :: l,m
   character(len=11), intent(out) :: name

   select case(l)
   case(1)
      name(1:1)='s'
      select case(m)
      case(1)
         name(2:11)='          '
      case default
         stop 'wrong m'
      end select
   case(2)
      name(1:1)='p'
      select case(m)
      case(1)
         name(2:11)='x         '
      case(2)
         name(2:11)='y         '
      case(3)
         name(2:11)='z         '
      case default
         stop 'wrong m'
      end select
   case(3)
      name(1:1)='d'        
      select case(m)
      case(1)
         name(2:11)='yz        '
      case(2)
         name(2:11)='xz        '
      case(3)
         name(2:11)='xy        '
      case(4)
         name(2:11)='x2-y2     '
      case(5)
         name(2:11)='2z2-r2    '
      case default
         stop 'wrong m'
      end select
   case(4)
      name(1:1)='f'        
      select case(m)
      case(1)
         name(2:11)='x(r2-5z2) '
      case(2)
         name(2:11)='y(r2-5z2) '
      case(3)
         name(2:11)='z(3r2-5z2)'
      case(4)
         name(2:11)='x(x2-3y2) '
      case(5)
         name(2:11)='y(y2-3x2) '
      case(6)
         name(2:11)='z(x2-y2)  '
      case(7)
         name(2:11)='xyz       '
      case default
         stop 'wrong m'
      end select
   case default
      stop 'l not recognized'
   end select

END SUBROUTINE shell_name


!> Perform a total DOS output.
subroutine global_analysis(orbs,wf,occopt,filename)
   use module_base
   use module_types
   !use public_enums
   use fermi_level, only: SMEARING_DIST_FERMI
   implicit none
   type(orbitals_data), intent(in) :: orbs
   real(gp), intent(in) :: wf
   integer , intent(in) :: occopt
   character(len = *), intent(in) :: filename

   integer, parameter :: DOS = 123456
   integer :: ikpt, iorb, index, i
   real(wp) :: minE, maxE, e


   ! We define a Gnuplot file.
   open(unit = DOS, file = trim(filename), action = "write")

   minE = 999_dp
   maxE = -999_dp

   write(DOS, "(A)") "#!/usr/bin/gnuplot"
   write(DOS, "(A)") "# DOS generated by BigDFT."
   write(DOS, "(A)")
   write(DOS, "(A)") "# Comment out to generate a EPS file"
   write(DOS, "(A)") "# set term postscript enhanced"
   write(DOS, "(A)") '# set output "dos.eps"'
   write(DOS, "(A)")
   write(DOS, "(A)") "# Comment out to generate a PNG file"
   write(DOS, "(A)") "# set term png font DejaVuSerif 10 size 600,480"
   write(DOS, "(A)") '# set output "dos.png"'
   write(DOS, "(A)")
   write(DOS, "(A)") "# This is the smearing value used in the calculation."
   write(DOS, "(A,F12.8,A)") "w = ", wf*Ha_eV,"  # eV"
  !write(DOS, "(A,F12.8,A)") "T = ", wf*Ha_K," K"
   write(DOS, "(A)")
   write(DOS, "(A)") "# This is the smearing function used in the calculation."
   if (occopt == SMEARING_DIST_FERMI) then
     write(DOS, "(A,F7.4,A,F12.6,A)") 'set title "Density of states, Fermi-Dirac smearing w = ', &
          & wf*Ha_eV, 'eV, E_f = ', orbs%efermi*Ha_eV , 'eV"'
     write(DOS, "(A)") "f(eb,E)  = 1 / (1 + exp((eb-E)/w))"
     write(DOS, "(A)") "df(eb,E) = 1 / (2 + exp((eb-E)/w) + exp((E-eb)/w)) / w"
   !elseif (occopt == SMEARING_DIST_ERF) then  
   else  ! to be changed for cold smearing and ... 
     write(DOS, "(A,F7.4,A,F12.6,A)") 'set title "Density of states, erf smearing w = ', &
          & wf*Ha_eV, 'eV,  E_f = ', orbs%efermi*Ha_eV , 'eV"'
     write(DOS, "(A)") "f(eb,E)  = 0.5 * (1 - erf((E - eb) / w))"
     write(DOS, "(A)") "df(eb,E) = exp(-((E - eb) / w) ** 2) / w / sqrt(pi)"
   end if
   write(DOS, "(A)")
   write(DOS, "(A)") "U(E) = " // char(92)
   do ikpt = 1, orbs%nkpts
      write(DOS, "(A,F12.8,A)") "  ", orbs%kwgts(ikpt), " * (" // char(92)
      index = 1
      do iorb = 0, (orbs%norbu - 1) / 6
         write(DOS, "(A)", advance = "NO") "   "
         do i = 1, 6
            e = orbs%eval(index+(ikpt-1)*orbs%norb)
            e = e*Ha_eV
            minE = min(e, minE)
            maxE = max(e, maxE)
            write(DOS, "(A,F12.8,A)", advance = "NO") "df(", e, ",E)"
            index = index + 1
            if (index > orbs%norbu) exit
            write(DOS, "(A)", advance = "NO") " + "
         end do
         write(DOS, "(A)") char(92)
      end do
      if (ikpt < orbs%nkpts) then
         write(DOS, "(A)") "  ) + " // char(92)
      else
         write(DOS, "(A)") "  )"
      end if
   end do
   if (orbs%norbd > 0) then
      write(DOS, "(A)")
      write(DOS, "(A)") "D(E) = " // char(92)
      do ikpt = 1, orbs%nkpts
         write(DOS, "(A,F12.8,A)") "  ", orbs%kwgts(ikpt), " * (" // char(92)
         index = orbs%norbu + 1
         do iorb = 0, (orbs%norbd - 1) / 6
            write(DOS, "(A)", advance = "NO") "   "
            do i = 1, 6
               e = orbs%eval(index+(ikpt-1)*orbs%norb)
               e = e*Ha_eV
               minE = min(e, minE)
               maxE = max(e, maxE)
               write(DOS, "(A,F12.8,A)", advance = "NO") "df(", e, ",E)"
               index = index + 1
               if (index > orbs%norb) exit
               write(DOS, "(A)", advance = "NO") " + "
            end do
            write(DOS, "(A)") char(92)
         end do
         if (ikpt < orbs%nkpts) then
            write(DOS, "(A)") "  ) + " // char(92)
         else
            write(DOS, "(A)") "  )"
         end if
      end do
   end if
   write(DOS, "(A)") "set samples 2500"
   write(DOS, "(A)") "set key bottom left"
   write(DOS, "(A)") 'set xlabel "Energy (eV)"'
   write(DOS, "(A)") 'set ylabel "States per unit cell per eV"'
   !write(DOS, "(A)") 'set ylabel "Electrons per unit cell per eV"'
   write(DOS, "(A,F12.6,A,F12.6,A)") "set arrow from ", orbs%efermi*Ha_eV , &
       & ",graph 0.95 to ", orbs%efermi*Ha_eV , ",graph 0.05 lt 0"
   write(DOS, "(A,F12.6,A)") "set label at  ", orbs%efermi*Ha_eV , &
       & ",graph 0.96  center 'E_f'"
   write(DOS, "(A,F12.8,A,F12.8,A)")  "plot [", minE-0.1*(maxE-minE) , &
      &   ":", maxE+0.1*(maxE-minE) , "] " // char(92)
   if (orbs%norbd > 0) then
      write(DOS, "(A)")  '  U(x) t "Spin up", -D(x) t "Spin down"'
   else
      write(DOS, "(A)")  "  2 * U(x) notitle"
   end if
   write(DOS, "(A)")  "pause -1"

   close(DOS)
END SUBROUTINE global_analysis
