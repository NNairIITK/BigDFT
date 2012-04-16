!> @file
!!  Program RISM
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


program rism
  use BigDFT_API
  implicit none
  character(len=*), parameter :: subname='rism'
  integer :: n1i,n2i,n3i,iproc,nproc,i_stat,i_all,nelec
  integer :: n3d,n3p,n3pi,i3xcsh,i3s,nlr,iat,nspin,istat
  real(gp) :: hxh,hyh,hzh
  type(atoms_data) :: atoms
  type(input_variables) :: in
  type(orbitals_data) :: orbs
  type(locreg_descriptors) :: Glr
  type(rho_descriptors)  :: rhodsc
  real(gp), dimension(3) :: shift
  integer, dimension(:), allocatable :: iatlr
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(gp), dimension(:), allocatable :: atchgs,radii
  real(gp), dimension(:,:), allocatable :: radii_cf
  real(gp), dimension(:,:), pointer :: rxyz
  real(dp), dimension(:,:,:,:), pointer :: rho,pot,pot_ion
  character(len=5) :: gridformat
  character(len=60) :: radical

  !for the moment no need to have parallelism
  iproc=0
  nproc=1

  call memocc_set_memory_limit(memorylimit)

  ! Read a possible radical format argument.
  call get_command_argument(1, value = radical, status = istat)
  if (istat > 0) then
     write(radical, "(A)") "input"
  end if


  !initalise the varaibles for the calculation
  !standard names
  call standard_inputfile_names(in,radical,nproc)
  call read_input_variables(iproc,'posinp',in,atoms,rxyz)
  write(gridformat, "(A)") ""
  select case (in%output_denspot_format)
     case (output_denspot_FORMAT_ETSF)
        write(gridformat, "(A)") ".etsf"
     case (output_denspot_FORMAT_CUBE)
        write(gridformat, "(A)") ".bin"
  end select

  if (iproc == 0) then
     call print_general_parameters(nproc,in,atoms)
  end if
       
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)

  !subtract the charge of the dummy atoms from the total number of electrons
  do iat=1,atoms%nat
     if (atoms%atomnames(atoms%iatype(iat)) == 'DD') then
        nelec=nelec-atoms%nelpsp(atoms%iatype(iat))
     end if
  end do

  print *,'nelec',nelec
  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,atoms,rxyz,radii_cf,in%crmult,in%frmult,in%hx,in%hy,in%hz,&
       Glr,shift)

  call createWavefunctionsDescriptors(iproc,in%hx,in%hy,in%hz,&
       atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr)

  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)

  call createDensPotDescriptors(iproc,nproc,atoms,Glr%d,hxh,hyh,hzh,&  !n(?) hxh, hyh, hzh (no value assigned)
       rxyz,in%crmult,in%frmult,radii_cf,in%nspin,'D',0,in%rho_commun,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)

  !read the files from the .cubes on the disk
  call read_density('electronic_density' // gridformat,atoms%geocode,&
       n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho)

  call read_density('hartree_potential' // gridformat,atoms%geocode,&
       n1i,n2i,n3i,nspin,hxh,hyh,hzh,pot)

  !not needed for the moment
  call read_density('ionic_potential' // gridformat,atoms%geocode,&
       n1i,n2i,n3i,nspin,hxh,hyh,hzh,pot_ion)

  !also extract the number of atoms and the positions
  !for the atomic densities, if they are present

  !start the atomic charges calculation
  allocate(atchgs(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atchgs,'atchgs',subname)

  allocate(iatlr(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,iatlr,'iatlr',subname)

  !decide three radii per atom
  nlr=0
  do iat=1,atoms%nat
     iatlr(iat)=3
     nlr=nlr+3
  end do
  !number of long range basis functions

  allocate(radii(nlr+ndebug),stat=i_stat)
  call memocc(i_stat,radii,'radii',subname)

  call assign_atomic_radii(atoms%nat,iatlr,nlr,radii)

  call atomic_charges(iproc,nproc,rxyz,iatlr,radii,atoms,nlr,nelec,Glr,ngatherarr,&
       hxh,hyh,hzh,n3p,i3s+i3xcsh,rho,atchgs) !n(m)

  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)
  i_all=-product(shape(ngatherarr))*kind(ngatherarr)
  deallocate(ngatherarr,stat=i_stat)
  call memocc(i_stat,i_all,'ngatherarr',subname)

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)


  call deallocate_lr(Glr,subname)
  
  call deallocate_orbs(orbs,subname)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)
  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
  i_all=-product(shape(radii))*kind(radii)
  deallocate(radii,stat=i_stat)
  call memocc(i_stat,i_all,'radii',subname)
  i_all=-product(shape(atchgs))*kind(atchgs)
  deallocate(atchgs,stat=i_stat)
  call memocc(i_stat,i_all,'atchgs',subname)

  i_all=-product(shape(iatlr))*kind(iatlr)
  deallocate(iatlr,stat=i_stat)
  call memocc(i_stat,i_all,'iatlr',subname)
  
  call deallocate_atoms(atoms,subname) 
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  call free_input_variables(in)

  call deallocate_rho_descriptors(rhodsc,subname)

  !finalize memory counting
  call memocc(0,0,'count','stop')


end program rism


subroutine assign_atomic_radii(nat,iatlr,nlr,radii)
  use module_base
  implicit none
  integer, intent(in) :: nat,nlr
  integer, dimension(nat), intent(in) :: iatlr
  real(gp), dimension(nlr), intent(out) :: radii
  !local variables
  integer :: iat,i,ilr
  real(gp) :: baseradius,width

  baseradius=1.0_gp
  width=1.5_gp

  ilr=0
  do iat=1,nat
     do i=1,iatlr(iat)
        ilr=ilr+1
        if (i==1) then
           radii(ilr)=baseradius
        else
           radii(ilr)=radii(ilr-1)*width
        end if
     end do
  end do

  if (ilr /= nlr) stop 'ERROR ilr'

END SUBROUTINE assign_atomic_radii


!>   Calculate atomic charges using Lee, York and Yang method 
!!   But with a basis similar to Blochl one 
!!   Refs: J.Chem.Phys. 102(19),7549 (1995)
!!         J.Chem.Phys. 103(17),7422 (1995) 
!!   use a basis of error functions centered on the atoms, with atom-defined radii
!!   and also short-range functions are allowed, as well as dummy atoms
subroutine atomic_charges(iproc,nproc,rxyz,iatlr,radii,atoms,nlr,nelec,lr,ngatherarr,& !n(c) pot (arg:l-1)
     hxh,hyh,hzh,n3p,i3s,rho,C)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,n3p,i3s,nelec,nlr
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(atoms_data), intent(in) :: atoms
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  integer, dimension(atoms%nat), intent(in) :: iatlr
  !n(c) real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p), intent(in) :: pot
  real(dp), dimension(lr%d%n1i*lr%d%n2i*n3p), intent(inout) :: rho
  real(gp), dimension(nlr) :: radii
  real(gp), dimension(3,atoms%nat) :: rxyz
  real(gp), dimension(atoms%nat) :: C
  !local variables
  character(len=*), parameter :: subname='atomic_charges'
  logical, parameter :: shortrange=.false.
  integer :: iat,info,nwork,i_all,i_stat,nbasis,i,j,ilr
  real(gp) :: ddotu,ddotv,gammafac
  type(gaussian_basis) :: Gpswf,Glongr
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: D,rhoarr,u,v,work,Caux
  real(gp), dimension(:,:), allocatable :: Hlrlr,Hwork,H,ovrlp
  real(gp), dimension(:), pointer :: Gocc

  if (atoms%geocode /= 'F') then
     write(*,*)'ERROR: the atomic charges can be calculated only in isolated BC!'
     stop
  end if
  !calculate the number of basis sets
  if (shortrange) then
     !consider also short-range behaviour
     nullify(Gpswf%rxyz)
     call gaussian_pswf_basis(31,.false.,iproc,1,atoms,rxyz,Gpswf,Gocc)
     if (associated(Gocc)) then
        i_all=-product(shape(Gocc))*kind(Gocc)
        deallocate(Gocc,stat=i_stat)
        call memocc(i_stat,i_all,'Gocc',subname)
        nullify(Gocc)
     end if
     nbasis=nlr+Gpswf%ncoeff
     !call gaussian_rism_basis_new(nlr,radii,rxyz,Glongr)
  else
     nbasis=nlr
  end if

  !allocate all the arrays
  allocate(H(nbasis,nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,H,'H',subname)
  allocate(rhoarr(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,rhoarr,'rhoarr',subname)
  allocate(D(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,D,'D',subname)
  allocate(u(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,u,'u',subname)
  allocate(v(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,v,'v',subname)
  allocate(Caux(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,Caux,'Caux',subname)

  if (shortrange) then
     !calculate the overlap matrix as well as the kinetic overlap
     !in view of complete gaussian calculation
     allocate(ovrlp(Gpswf%ncoeff,Gpswf%ncoeff+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     
     !overlap calculation of the kinetic operator
     call kinetic_overlap(Gpswf,Gpswf,ovrlp)
  
     !fill the last part of the H matrix
     !for the same structure the kinetic overlap is symmetric
     !H should be 1/4pi Delta but we have just calculated -1/2 Delta
     !so we have to divide by 2pi.
     do j=1,Gpswf%ncoeff
        do i=1,j-1
           H(nlr+i,nlr+j)=-1.0_gp/(8.0_gp*atan(1.0_gp))*ovrlp(i,j)
           H(nlr+j,nlr+i)=-1.0_gp/(8.0_gp*atan(1.0_gp))*ovrlp(i,j)
        end do
        !diagonal elements
        H(nlr+j,nlr+j)=-1.0_gp/(8.0_gp*atan(1.0_gp))*ovrlp(j,j)
     end do

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)

     allocate(ovrlp(Gpswf%ncoeff,Glongr%ncoeff+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     
     !overlap between longrange basis and short-range basis
     call gaussian_overlap(Gpswf,Glongr,ovrlp)
     
     !fill the block off-diagonal part of the H matrix
     print *,'coeffs',Gpswf%ncoeff,Glongr%ncoeff
     do j=1,Glongr%ncoeff
        do i=1,Gpswf%ncoeff
           if (Gpswf%ncoeff == Glongr%ncoeff .and. i > j) then
              H(nlr+i,j)=ovrlp(j,i)
              H(j,nlr+i)=ovrlp(j,i)
           else
              H(nlr+i,j)=ovrlp(i,j)
              H(j,nlr+i)=ovrlp(i,j)
           end if
        end do
     end do
     
     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
  end if

  !here nlr the number of long range basis functions
  !calculate H matrix
  allocate(Hlrlr(nlr,nlr+ndebug),stat=i_stat)
  call memocc(i_stat,Hlrlr,'Hlrlr',subname)

  call two_center_two_electrons_analytic(nlr,atoms%nat,iatlr,radii,rxyz,Hlrlr)
 
  !fill the first part of the H matrix
  do j=1,nlr
     do i=1,nlr
        H(i,j)=Hlrlr(i,j)
     end do
  end do

  i_all=-product(shape(Hlrlr))*kind(Hlrlr)
  deallocate(Hlrlr,stat=i_stat)
  call memocc(i_stat,i_all,'Hlrlr',subname)

  if (iproc == 0) then
     do iat=1,nbasis
        write(*,'(a,i0,10(1pe15.7))')'H',iat,H(1:iat,iat)
     end do
  end if
 
  !calculate the long range part of the density
  call calculate_rho_longrange(iproc,nproc,atoms,nlr,iatlr,radii,rxyz,hxh,hyh,hzh,&
       lr%d%n1,lr%d%n2,lr%d%n3,n3p,i3s,lr%d%n1i,lr%d%n2i,rho,rhoarr)

  if (iproc == 0) then
     do iat=1,nlr
        write(*,'(a,i0,10(1pe15.7))')'rhoarrrho',iat,rhoarr(iat)
     end do
  end if

  if (shortrange) then
     !calculate the shortrange part of the rho array
     call calculate_rho_shortrange(iproc,nproc,atoms,lr,Gpswf,hxh,hyh,hzh,rxyz,ngatherarr,&
          rho,rhoarr(nlr+1))
     if (iproc == 0) then
        do iat=1,nbasis
           write(*,'(a,i0,10(1pe15.7))')'rhoarr',iat,rhoarr(iat)
        end do
     end if

     nullify(Gpswf%rxyz)
     call deallocate_gwf(Gpswf,subname)
     nullify(Glongr%rxyz)
     call deallocate_gwf(Glongr,subname)
  end if


  !initalise D array
  D=0.0_gp
  do iat=1,nlr
     D(iat)=1.0_gp
  end do


  !calculate workspace array for dsysv
  !temporary allocation of the work array, workspace query in dsysv
  allocate(Hwork(nbasis,nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,Hwork,'Hwork',subname)
  allocate(iwork(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,iwork,'iwork',subname)
  allocate(work(100+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  call dcopy(nbasis*nbasis,H,1,Hwork,1)
  call dcopy(nbasis,rhoarr,1,u,1)
  !workspace query
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),u(1),nbasis,work(1),-1,info)
  nwork=work(1)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  !here we start with the linear algebra
  call dcopy(nbasis,rhoarr,1,u,1)
  call dcopy(nbasis*nbasis,H,1,Hwork,1)
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),u(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of u',info
  end if

  !here we start with the linear algebra
  call dcopy(nbasis,D,1,v,1)
  call dcopy(nbasis*nbasis,H,1,Hwork,1)
!!$  if (iproc == 0) print *,'here',v(:)
!!$  if (iproc == 0) print '(a,4(1pe15.7))','000',Hwork(1,1),Hwork(1,2),Hwork(2,1),Hwork(2,2)
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),v(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of v',info
  end if
!!$  if (iproc == 0) print '(a,4(1pe15.7))','there',v(:)
!!$  if (iproc == 0) print '(a,4(1pe15.7))','AAA',H(1,1),H(1,2),H(2,1),H(2,2)
!!$  !determinant of the matrix, temporary
!!$  gammafac=H(1,1)*H(2,2)-H(1,2)*H(2,1)
!!$  if (iproc == 0) print '(a,4(1pe15.7))','res',(H(2,2)-H(1,2))/gammafac,&
!!$       (H(1,1)-H(2,1))/gammafac,real(nelec,gp)
  ddotu=dot(nbasis,D(1),1,u(1),1)
  ddotv=dot(nbasis,D(1),1,v(1),1)
  
  gammafac=(real(nelec,gp)+ddotu)/ddotv

  !!zero has to be put when the potential is the deformation potential
  !gammafac=ddotu/ddotv

  !calculate the eigenvalues of H
  call dcopy(nbasis*nbasis,H,1,Hwork,1)
  !print *,'nwork',nwork,3*nbasis-1
  call dsyev('N','U',nbasis,Hwork(1,1),nbasis,v(1),work(1),nwork,info)

  if (iproc == 0) then
     do iat=1,nbasis
        write(*,'(a,i0,10(1pe15.7))')'eigenvalues',iat,v(iat)
     end do
     write(*,*)'Condition Number (valid only if all eigenvalues are negative):',v(1)/v(nbasis)
  end if



!!$  if (iproc == 0 )print *,'gamma',gammafac,nelec,v(:)

  !rescale D
  call vscal(nbasis,gammafac,D(1),1)

  !initalise C
  call dcopy(nbasis,D,1,Caux,1)
  
  !fill it 
  call axpy(nbasis,-1.0_gp,rhoarr(1),1,Caux(1),1)

  !solve the system
  call dsysv('U',nbasis,1,H(1,1),nbasis,iwork(1),Caux(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of charges',info
  end if

  !for any atom calculates the sum of the long range term to have the charge
  ilr=0
  do iat=1,atoms%nat
     C(iat)=0.0_gp
     do i=1,iatlr(iat)
        ilr=ilr+1
        C(iat)=C(iat)+Caux(ilr)
        print *,'Caux,iat,i',Caux(ilr),iat,i
     end do
  end do
  !print the charges
  if (iproc == 0) then
     print *,'nelec',nelec
     do iat=1,atoms%nat
        write(*,'(1x,a,i4,2(1x,f12.7))')'atom, charge',iat,C(iat),C(iat)-real(atoms%nelpsp(atoms%iatype(iat)),gp)
     end do
  end if

  i_all=-product(shape(iwork))*kind(iwork)
  deallocate(iwork,stat=i_stat)
  call memocc(i_stat,i_all,'iwork',subname)
  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  i_all=-product(shape(Hwork))*kind(Hwork)
  deallocate(Hwork,stat=i_stat)
  call memocc(i_stat,i_all,'Hwork',subname)

  i_all=-product(shape(H))*kind(H)
  deallocate(H,stat=i_stat)
  call memocc(i_stat,i_all,'H',subname)

  i_all=-product(shape(rhoarr))*kind(rhoarr)
  deallocate(rhoarr,stat=i_stat)
  call memocc(i_stat,i_all,'rhoarr',subname)

  i_all=-product(shape(D))*kind(D)
  deallocate(D,stat=i_stat)
  call memocc(i_stat,i_all,'D',subname)

  i_all=-product(shape(u))*kind(u)
  deallocate(u,stat=i_stat)
  call memocc(i_stat,i_all,'u',subname)

  i_all=-product(shape(v))*kind(v)
  deallocate(v,stat=i_stat)
  call memocc(i_stat,i_all,'v',subname)

  i_all=-product(shape(Caux))*kind(Caux)
  deallocate(Caux,stat=i_stat)
  call memocc(i_stat,i_all,'Caux',subname)


END SUBROUTINE atomic_charges


subroutine two_center_two_electrons_analytic(nlr,nat,iatlr,radii,rxyz,H)
  use module_base
  implicit none
  integer, intent(in) :: nat,nlr
  integer, dimension(nat), intent(in) :: iatlr
  real(gp), dimension(nlr), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(gp), dimension(nlr,nlr), intent(out) :: H
  !local variables
  integer :: iat,jat,i,j,ilr,jlr
  real(gp) :: ra2,ra2pb2,rab2,erfor
  real(gp), dimension(3) :: A

  !the loop is so fast that we just do not need to process only the upeer triangular part
  ilr=0
  do iat=1,nat
     A(1)=rxyz(1,iat)
     A(2)=rxyz(2,iat)
     A(3)=rxyz(3,iat)
     do i=1,iatlr(iat)
        ilr=ilr+1
        ra2=radii(ilr)**2
        jlr=0
        do jat=1,nat
           rab2=(rxyz(1,jat)-A(1))**2
           rab2=rab2+(rxyz(2,jat)-A(2))**2
           rab2=rab2+(rxyz(3,jat)-A(3))**2
           do j=1,iatlr(jat)
              jlr=jlr+1
              ra2pb2=ra2+radii(jlr)**2
              H(ilr,jlr)=-erfor(sqrt(rab2),sqrt(ra2pb2))
           end do
        end do
     end do
  end do

  if (ilr /= nlr .or. jlr /= nlr) then
     write(*,*)'ilr or jlr error',ilr,jlr,nlr
     stop
  end if

  !!copy the same values on the lower triangular part
  !do ilr=1,nlr
  !   do jlr=ilr+1,nlr
  !      H(jlr,ilr)=H(ilr,jlr)
  !   end do
  !end do

END SUBROUTINE two_center_two_electrons_analytic


subroutine calculate_rho_longrange(iproc,nproc,at,nlr,iatlr,radii,rxyz,hxh,hyh,hzh,&
     n1,n2,n3,n3pi,i3s,n1i,n2i,rho,rhoarr)
  use module_base
  use module_types
  implicit none
  !Arguments---------
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,nlr
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  integer, dimension(at%nat), intent(in) :: iatlr
  real(gp), dimension(nlr), intent(in) :: radii
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(n1i*n2i*n3pi), intent(inout) :: rho
  real(gp), dimension(nlr), intent(out) :: rhoarr
  !Local variables---------
  logical :: perx,pery,perz,gox,goy,goz
  integer :: i1,i2,i3,ind,iat,ierr,i,ilr
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,j1,j2,j3
  real(gp) :: cutoff,rloc,Rel !n(c) pi
  real(gp) :: rx,ry,rz,x,y,z,r2,charge,erfor

  !n(c) pi=4.d0*atan(1.d0)

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  charge=0.0_gp
  ilr=0
  do iat=1,at%nat
     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)

     !cutoff=10.d0*rloc
     cutoff=0.0_gp!cutofrac!*rloc

     do i=1,iatlr(iat)
        rhoarr(ilr+i)=0.0_gp
     end do
     if (n3pi >0 ) then
        do i3=-nbl3,2*n3+1+nbr3
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           if (j3 >= i3s .and. j3 <= i3s+n3pi-1 ) then
              do i2=-nbl2,2*n2+1+nbr2
                 y=real(i2,kind=8)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 if (goy) then
                    do i1=-nbl1,2*n1+1+nbr1
                       x=real(i1,kind=8)*hxh-rx
                       call ind_positions(perx,i1,n1,j1,gox)
                       r2=x**2+y**2+z**2
                       if (gox .and. r2 > cutoff**2 ) then
                          ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                          Rel=rho(ind)
                          do i=1,iatlr(iat)
                             rloc=radii(ilr+i)
                             rhoarr(ilr+i)=rhoarr(ilr+i)+Rel*erfor(sqrt(r2),rloc)
                             charge=charge+erfor(sqrt(r2),rloc)
                          end do
                       end if
                    end do
                 end if
              end do
           end if
        end do
     end if
     !final results
     do i=1,iatlr(iat)
        rhoarr(ilr+i)=(hxh*hyh*hzh)*rhoarr(ilr+i)
     end do
     ilr=ilr+iatlr(iat)

  end do
  charge=hxh*hyh*hzh*charge
  !reduce the results
  if (nproc > 1) then
     call mpiallred(rhoarr(1),nlr,MPI_SUM,MPI_COMM_WORLD,ierr)
     !the same should be done for the charge
     call mpiallred(charge,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  if (iproc == 0) write(*,*)'Charge:',charge

END SUBROUTINE calculate_rho_longrange


!> erf(r/(sqrt(2)rl)/r
function erfor(r,rl)
  use module_base
  implicit none
  real(gp), intent(in) :: r,rl
  real(gp) :: erfor
  !local variables
  real(gp) :: sq2rl,pi,derf_val

  pi=4.0_gp*atan(1.0_gp)

  sq2rl=sqrt(2.0_gp)*rl
  
  if (r == 0.0_gp) then
     erfor=2.0_gp/(sqrt(pi)*sq2rl)
  else
     call derf_ab(derf_val,r/sq2rl)
     erfor=derf_val/r
  end if

end function erfor
  

!>   Gaussian basis associated to the long-range term of rism calculation
subroutine gaussian_rism_basis(nat,radii,rxyz,G)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(nat), intent(in) :: radii
  real(gp), dimension(3,nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: G  
  !local variables
  character(len=*), parameter :: subname='gaussian_psp_basis'
  real(gp), parameter :: oneo2pi3halves=0.0634936359342409697857633_gp
  integer :: iat,nshell,iexpo,l,ishell,i_stat

  G%nat=nat
  G%rxyz => rxyz
  allocate(G%nshell(G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)
 
  G%nshltot=0
  do iat=1,G%nat
     nshell=1
     G%nshell(iat)=nshell
     G%nshltot=G%nshltot+nshell
  end do

  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)

  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  do iat=1,G%nat
     do l=1,4 
        if (l==1) then
           ishell=ishell+1
           G%ndoc(ishell)=1
           G%nam(ishell)=l
           G%nexpo=G%nexpo+1
           G%ncoeff=G%ncoeff+2*l-1
        end if
     enddo
  end do

  !allocate and assign the exponents and the coefficients
  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  ishell=0
  iexpo=0
  do iat=1,G%nat
     do l=1,4 
        if (l==1) then
           ishell=ishell+1
           iexpo=iexpo+1
           G%psiat(iexpo)=-oneo2pi3halves/radii(iat)**3
           G%xp(iexpo)=radii(iat)
        end if
     end do
  end do

END SUBROUTINE gaussian_rism_basis


!> calculate the second part, by expressing the atomic wavefunctions on a real grid
subroutine calculate_rho_shortrange(iproc,nproc,at,lr,Gpswf,hxh,hyh,hzh,rxyz,ngatherarr,&
     rho,rhoarr)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  type(gaussian_basis), intent(in) :: Gpswf
  type(locreg_descriptors), intent(in) :: lr
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(ngatherarr(iproc,1)), intent(inout) :: rho
  real(dp), dimension(Gpswf%ncoeff), intent(out) :: rhoarr
  !local variables
  character(len=*), parameter :: subname='calculate_rho_shortrange'
  logical, parameter :: deltarho=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,iat,isx,isy,isz,iex,iey,iez,n3pi,i3s
  integer :: i1,i2,i3,j1,j2,j3,ityp,ind,i_stat,i_all,ierr,jorb,isorb,jproc,kproc
  real(gp) :: pi,rx,ry,rz,charge,rloc,cutoff,x,y,z,r2,arg,xp,hfac
  type(orbitals_data) :: orbspswf
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncoeff_par
  real(wp), dimension(:,:), allocatable :: psigaup,psi
  real(wp), dimension(:,:,:), allocatable :: psir
  real(dp), dimension(:,:,:), allocatable :: rhotot

  !here we have to add the density the compensation value and gather the complete result
  !----extracted from createIonicPotential
  pi=4.0_gp*atan(1.0_gp)
  ! Ionic energy (can be calculated for all the processors)

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  n3pi=ngatherarr(iproc,1)/(lr%d%n1i*lr%d%n2i)
  i3s=ngatherarr(iproc,2)/(lr%d%n1i*lr%d%n2i)

  if (n3pi >0 .and. deltarho) then

     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        cutoff=10.0_gp*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        do i3=isz,iez
           z=real(i3,gp)*hzh-rz
           call ind_positions(perz,i3,lr%d%n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,gp)*hyh-ry
              call ind_positions(pery,i2,lr%d%n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,gp)*hxh-rx
                 call ind_positions(perx,i1,lr%d%n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5_gp*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*lr%d%n1i+(j3-i3s+1-1)*lr%d%n1i*lr%d%n2i
                    rho(ind)=rho(ind)-xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if
  !-----end of extraction

  !the gathering should be done here
  allocate(rhotot(lr%d%n1i,lr%d%n2i,lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,rhotot,'rhotot',subname)

  if (nproc > 1) then
     call MPI_ALLGATHERV(rho,ngatherarr(iproc,1),&
          mpidtypw,rhotot,ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
  else
     !not efficient, the copy can be saved
     call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i,rho,1,rhotot,1)
  end if

  allocate(ncoeff_par(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ncoeff_par,'ncoeff_par',subname)


  !determine the number of basis functions treated by each processor
  ncoeff_par(0:nproc-1,1)=0
  do jorb=1,Gpswf%ncoeff
     jproc=mod(jorb-1,nproc)
     ncoeff_par(jproc,1)=ncoeff_par(jproc,1)+1
  end do
  !starting point
  do jproc=0,nproc-1
     isorb=0
     do kproc=0,jproc-1
        isorb=isorb+ncoeff_par(kproc,1)
     end do
     ncoeff_par(jproc,2)=isorb
  end do
  isorb=ncoeff_par(iproc,2)

  !fake orbitals descriptor to calculate the wavelet expansion
  call orbitals_descriptors(0,1,ncoeff_par(iproc,1),ncoeff_par(iproc,1),0,1,1,1,&
       (/0.0_gp,0.0_gp,0.0_gp/),(/1.0_gp /),orbspswf,.false.)

  !allocate the wavefunctions in wavelets and in gaussians
  allocate(psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncoeff_par(iproc,1)+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  allocate(psigaup(Gpswf%ncoeff,ncoeff_par(iproc,1)+ndebug),stat=i_stat)
  call memocc(i_stat,psigaup,'psigaup',subname)

  !the coefficients are elements of the identity matrix
  call razero(Gpswf%ncoeff*ncoeff_par(iproc,1),psigaup)
  do jorb=1,ncoeff_par(iproc,1)
     psigaup(isorb+jorb,jorb)=1.0_gp
  end do

  !calculate the atomic wavefunctions in Daubechies basis
  !put nproc=1 since the orbitals descriptors is not parallelised by choice
  call gaussians_to_wavelets_new_h(0,1,lr,orbspswf,2.0_gp*hxh,2.0_gp*hyh,2.0_gp*hzh,Gpswf,&
       psigaup,psi)

  call deallocate_orbs(orbspswf,subname)

  !initalize work arrays for decompressing the wavefunction
  call initialize_work_arrays_sumrho(lr,w)

  allocate(psir(lr%d%n1i,lr%d%n2i,lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if

  do jorb=1,ncoeff_par(iproc,1)

     hfac=(hxh*hyh*hzh)
     call daub_to_isf(lr,w,psi(1,jorb),psir)

     !now the integral with the density
     !the square root of the volume units should be put because the psi is normalised to 
     !have unit norm without volume unit
     rhoarr(jorb+isorb)=sqrt(hfac)*dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,rhotot(1,1,1),1,psir(1,1,1),1)
     !check the norm of the wavefunctions
!!$     print *,'norm',jorb+isorb,&
!!$          dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir(1,1,1),1,psir(1,1,1),1),&
!!$          rhoarr(jorb+isorb),&
!!$          sum(rhotot),hfac

     
  end do

  !gather the results in the same array
  if (nproc > 1) then
     call MPI_ALLGATHERV(rhoarr(isorb+min(1,ncoeff_par(iproc,1))),ncoeff_par(iproc,1),&
          mpidtypw,rhoarr(1),ncoeff_par(0,1),&
          ncoeff_par(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
  end if


  !deallocate the deallocatable
  call deallocate_work_arrays_sumrho(w)

  i_all=-product(shape(rhotot))*kind(rhotot)
  deallocate(rhotot,stat=i_stat)
  call memocc(i_stat,i_all,'rhotot',subname)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

  i_all=-product(shape(psigaup))*kind(psigaup)
  deallocate(psigaup,stat=i_stat)
  call memocc(i_stat,i_all,'psigaup',subname)

  i_all=-product(shape(ncoeff_par))*kind(ncoeff_par)
  deallocate(ncoeff_par,stat=i_stat)
  call memocc(i_stat,i_all,'ncoeff_par',subname)
  
END SUBROUTINE calculate_rho_shortrange
