!!!> Choose among the wavefunctions a subset of them
!!!! Rebuild orbital descriptors for the new space and allocate the psi_as wavefunction
!!!! By hypothesis the work array is big enough to contain both wavefunctions
!!!! This routine has to be tested
!!subroutine select_active_space(iproc,nproc,orbs,comms,mask_array,Glr,orbs_as,comms_as,psi,psi_as)
!!   use module_base
!!   use module_types
!!   use module_interfaces, except_this_one => select_active_space
!!   use communications_base, only: comms_cubic
!!   use communications_init, only: orbitals_communicators
!!   implicit none
!!   integer, intent(in) :: iproc,nproc
!!   type(orbitals_data), intent(in) :: orbs
!!   type(locreg_descriptors), intent(in) :: Glr
!!   type(comms_cubic), intent(in) :: comms
!!   logical, dimension(orbs%norb*orbs%nkpts), intent(in) :: mask_array
!!   real(wp), dimension(orbs%npsidim_comp), intent(in) :: psi
!!   type(orbitals_data), intent(out) :: orbs_as
!!   type(comms_cubic), intent(out) :: comms_as
!!   real(wp), dimension(:), pointer :: psi_as
!!   !local variables
!!   character(len=*), parameter :: subname='select_active_space'
!!   integer :: iorb,ikpt,norbu_as,norbd_as,icnt,ikptp,ispsi,ispsi_as
!!   integer :: i_stat,nvctrp
!!
!!   !count the number of orbitals of the active space
!!   norbu_as=-1
!!   norbd_as=-1
!!   do ikpt=1,orbs%nkpts
!!      icnt=0
!!      do iorb=1,orbs%norbu
!!         if (mask_array(iorb+(ikpt-1)*orbs%norb)) icnt=icnt+1
!!      end do
!!      if (norbu_as /= icnt .and. norbu_as /= -1) then
!!         write(*,*)'ERROR(select_active_space): the mask array should define always the same norbu'
!!         stop
!!      end if
!!      norbu_as=icnt
!!      icnt=0
!!      do iorb=orbs%norbu+1,orbs%norbu+orbs%norbd
!!         if (mask_array(iorb+(ikpt-1)*orbs%norb)) icnt=icnt+1
!!      end do
!!      if (norbd_as /= icnt .and. norbd_as /= -1) then
!!         write(*,*)'ERROR(select_active_space): the mask array should define always the same norbd'
!!         stop
!!      end if
!!      norbd_as=icnt
!!   end do
!!
!!   !allocate the descriptors of the active space
!!   call orbitals_descriptors(iproc,nproc,norbu_as+norbd_as,norbu_as,norbd_as, &
!!        orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbs_as,&
!!        .false.,basedist=orbs%norb_par(0:,1))
!!   !allocate communications arrays for virtual orbitals
!!   call orbitals_communicators(iproc,nproc,Glr,orbs_as,comms_as,basedist=comms_as%nvctr_par(0:,1))  
!!   !allocate array of the eigenvalues
!!   allocate(orbs_as%eval(orbs_as%norb*orbs_as%nkpts+ndebug),stat=i_stat)
!!   call memocc(i_stat,orbs_as%eval,'orbs_as%eval',subname)
!!
!!   !fill the orbitals array with the values and the wavefunction in transposed form
!!   icnt=0
!!   do iorb=1,orbs%nkpts*orbs%norb
!!      if (mask_array(iorb)) then
!!         icnt=icnt+1
!!         orbs_as%eval(icnt)=orbs%eval(iorb)
!!      end if
!!   end do
!!   if (icnt/=orbs_as%norb*orbs_as%nkpts) stop 'ERROR(select_active_space): icnt/=orbs_as%norb*orbs_as%nkpts'
!!
!!   allocate(psi_as(orbs_as%npsidim_comp+ndebug),stat=i_stat)
!!   call memocc(i_stat,psi_as,'psi_as',subname)
!!
!!   ispsi=1
!!   do ikptp=1,orbs%nkptsp
!!      ikpt=orbs%iskpts+ikptp
!!      nvctrp=comms%nvctr_par(iproc,ikpt) 
!!      !this should be identical in both the distributions
!!      if (nvctrp /= comms_as%nvctr_par(iproc,ikpt)) then
!!         write(*,*)'ERROR(select_active_space): the component distrbution is not identical'
!!         stop
!!      end if
!!
!!      !put all the orbitals which match the active space
!!      ispsi=1
!!      ispsi_as=1
!!      do iorb=1,orbs%norb
!!         if (mask_array(iorb+(ikpt-1)*orbs%norb)) then
!!            call vcopy(nvctrp,psi(ispsi),1,psi_as(ispsi_as),1)
!!            ispsi_as=ispsi_as+nvctrp*orbs_as%nspinor
!!         end if
!!         ispsi=ispsi+nvctrp*orbs%nspinor
!!      end do
!!   end do
!!
!!END SUBROUTINE select_active_space

