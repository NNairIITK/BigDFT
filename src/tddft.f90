subroutine tddft_casida(iproc,nproc,atoms,rxyz,hxh,hyh,hzh,n3p,n3parr,Glr,orbs,orbsv,i3s,fxc,pkernelseq,psi,psiv)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,n3p,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data), intent(in) :: orbs,orbsv
  type(locreg_descriptors), intent(in) :: Glr
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: fxc
  real(wp), dimension(*), intent(in) :: pkernelseq
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(wp), dimension(orbsv%npsidim_orbs), intent(in) :: psiv
  !local variables
  character(len=*), parameter :: subname='tddft_casida'
  integer :: i_all,i_stat
  real(gp), dimension(3) :: chargec
  real(wp), dimension(:), allocatable :: psirocc,psirvirt

  !temporary call to the coupling matrix calculation
  allocate(psirocc(max(max(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*orbs%norbp,&
       n3parr(0)*orbs%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psirocc,'psirocc',subname)

  allocate(psirvirt(max(max(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*orbsv%norbp,&
       n3parr(0)*orbsv%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psirvirt,'psirvirt',subname)

  call prepare_psirocc(iproc,nproc,Glr,orbs,n3p,n3parr,psi,psirocc)

  call prepare_psirocc(iproc,nproc,Glr,orbsv,n3p,n3parr,psiv,psirvirt)

  call center_of_charge(atoms,rxyz,chargec)

  call coupling_matrix_prelim(iproc,nproc,atoms%geocode,orbs%nspin,Glr,orbs,orbsv,&
       i3s,n3p,hxh,hyh,hzh,chargec,pkernelseq,fxc,psirocc,psirvirt)

  i_all=-product(shape(psirocc))*kind(psirocc)
  deallocate(psirocc,stat=i_stat)
  call memocc(i_stat,i_all,'psirocc',subname)

  i_all=-product(shape(psirvirt))*kind(psirvirt)
  deallocate(psirvirt,stat=i_stat)
  call memocc(i_stat,i_all,'psirvirt',subname)

end subroutine tddft_casida
