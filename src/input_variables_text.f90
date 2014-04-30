!> @file
!!  Routines to read and print input variables from old fixed-format files
!!  To be declared deprecated in version 1.8 of the code
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Define the name of the input files
subroutine standard_inputfile_names(in, radical)
  use module_types
  use module_base
  use yaml_output
  implicit none
  type(input_variables), intent(inout) :: in
  character(len = *), intent(in) :: radical

  call set_inputfile(in%file_lin, radical,    "lin")
  call set_inputfile(in%file_frag, radical,   "frag")

  in%files = INPUTS_NONE
END SUBROUTINE standard_inputfile_names


!> Set and check the input file
!! if radical is empty verify if the file input.ext exists. 
!! otherwise search for radical.ext
!! if the so defined file is not existing, then filename becomes default.ext
subroutine set_inputfile(filename, radical, ext)
  implicit none
  character(len = *), intent(in) :: radical, ext
  character(len = 100), intent(out) :: filename
  
  logical :: exists

  write(filename, "(A)") ""
  if (trim(radical) == "") then
     write(filename, "(A,A,A)") "input", ".", trim(ext)
  else
     write(filename, "(A,A,A)") trim(radical), ".", trim(ext)
  end if

  inquire(file=trim(filename),exist=exists)
  if (.not. exists .and. (trim(radical) /= "input" .and. trim(radical) /= "")) &
       & write(filename, "(A,A,A)") "default", ".", trim(ext)
end subroutine set_inputfile

subroutine input_variables_from_old_text_format(in, atoms, run_name)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  character(len = *), intent(in) :: run_name

  integer :: ierr

  call set_inputfile(in%file_lin, run_name,    "lin")
  call set_inputfile(in%file_frag, run_name,   "frag")

  ! To avoid race conditions where procs create the default file and other test its
  ! presence, we put a barrier here.
  if (bigdft_mpi%nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)

  ! Linear scaling (if given)
  !in%lin%fragment_calculation=.false. ! to make sure that if we're not doing a linear calculation we don't read fragment information
  call lin_input_variables_new(bigdft_mpi%iproc,in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR, trim(in%file_lin),in,atoms)

  ! Fragment information (if given)
  call fragment_input_variables(bigdft_mpi%iproc,(in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR).and.in%lin%fragment_calculation,trim(in%file_frag),in,atoms)
END SUBROUTINE input_variables_from_old_text_format

!> Read linear input parameters
subroutine lin_input_variables_new(iproc,dump,filename,in,atoms)
  use module_base, dictionary_value_lgt => max_field_length
  use module_types
  use module_input
  use yaml_output, only: yaml_map
  use module_input_keys
  implicit none
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: filename
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  logical, intent(in) :: dump
  !local variables
  logical :: exists
  character(len=*), parameter :: subname='lin_input_variables'
  character(len=256) :: comments
  !logical,dimension(atoms%astruct%ntypes) :: parametersSpecified
  logical :: dummy_bool!found
  character(len=20):: atomname
  integer :: itype, jtype, ios, ierr, iat, npt, iiorb, iorb, nlr, istat
  integer :: dummy_int
  integer, dimension(2) :: dummy_iarr
  real(gp), dimension(2) :: dummy_darr
  real(gp):: ppao, ppl, pph, lrl, lrh, kco_FOE, kco, dummy_real
!  real(gp),dimension(atoms%astruct%ntypes) :: locradType, locradType_lowaccur, locradType_highaccur
  type(dictionary), pointer :: dict,dict_basis
  character(len=dictionary_value_lgt) :: dummy_char,dummy_char2


  !this variable is temporay for the moment, it should integrate the input dictionary
  call dict_init(dict)

  !Linear input parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Linear Parameters')  
  if (exists) in%files = in%files + INPUTS_LIN

  ! number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)
  comments='number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)'
  call input_var(dummy_int,'2',ranges=(/1,2/),comment=comments)
  call dict_set(dict//LIN_GENERAL//HYBRID,dummy_int==1)
  ! number of iterations
  comments = 'outer loop iterations (low, high)'
  call input_var(dummy_int,'15',dict//LIN_GENERAL//NIT//0,ranges=(/0,100000/))
  call input_var(dummy_int,'1',dict//LIN_GENERAL//NIT//1,ranges=(/0,100000/),comment=comments)

  comments = 'basis iterations (low, high)'
  call input_var(dummy_int,'12',dict//LIN_BASIS//NIT//0,ranges=(/0,100000/))
  call input_var(dummy_int,'50',dict//LIN_BASIS//NIT//1,ranges=(/0,100000/),comment=comments)

  comments = 'kernel iterations (low, high) - directmin only'
  call input_var(dummy_int,'1',dict//LIN_KERNEL//NSTEP//0,ranges=(/0,1000/))
  call input_var(dummy_int,'1',dict//LIN_KERNEL//NSTEP//1,ranges=(/0,1000/),comment=comments)

  comments = 'density iterations (low, high)'
  call input_var(dummy_int,'15',dict//LIN_KERNEL//NIT//0,ranges=(/0,1000/))
  call input_var(dummy_int,'15',dict//LIN_KERNEL//NIT//1,ranges=(/0,1000/),comment=comments)

  ! DIIS history lengths
  comments = 'DIIS history for basis (low, high)'
  call input_var(dummy_int,'5',dict//LIN_BASIS//IDSX//0,ranges=(/0,100/))
  call input_var(dummy_int,'0',dict//LIN_BASIS//IDSX//1,ranges=(/0,100/),comment=comments)

  comments = 'DIIS history for kernel (low, high) - directmin only'
  call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX_COEFF//0,ranges=(/0,100/))
  call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX_COEFF//1,ranges=(/0,100/),comment=comments)

  comments = 'DIIS history for density mixing (low, high)'
  call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX//0,ranges=(/0,100/))
  call input_var(dummy_int,'0',dict//LIN_KERNEL//IDSX//1,ranges=(/0,100/),comment=comments)

  ! mixing parameters
  comments = 'density mixing parameter (low, high)'
  call input_var(dummy_real,'.5d0',dict//LIN_KERNEL//ALPHAMIX//0,ranges=(/0.d0,1.d0/))
  call input_var(dummy_real,'.5d0',dict//LIN_KERNEL//ALPHAMIX//1,ranges=(/0.d0,1.d0/),comment=comments)

  ! Convergence criteria
  comments = 'outer loop convergence (low, high)'
  call input_var(dummy_real,'1.d-8' ,dict//LIN_GENERAL//RPNRM_CV//0,ranges=(/0.d0,1.d0/))
  call input_var(dummy_real,'1.d-12',dict//LIN_GENERAL//RPNRM_CV//1,ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'basis convergence (low, high) ; early stop TMB optimization, dynamic gnrm (experimental mode only)'
  call input_var(dummy_real,'1.d-3',dict//LIN_BASIS//GNRM_CV//0,ranges=(/0.0_gp,1.0_gp/))
  call input_var(dummy_real,'1.d-5',dict//LIN_BASIS//GNRM_CV//1,ranges=(/0.0_gp,1.0_gp/))

  call input_var(dummy_real,'1.d-4',dict//LIN_BASIS//DELTAE_CV,ranges=(/0.0_gp,1.0_gp/))
  call input_var(dummy_real,'1.d-4',dict//LIN_BASIS//GNRM_DYN,ranges=(/0.0_gp,1.0_gp/),comment=comments)

  !these variables seems deprecated, put them to their default value

  comments = 'factor to reduce the confinement. Only used for hybrid mode.'
  call input_var(dummy_real,'0.5d0',dict//LIN_GENERAL//CONF_DAMPING,ranges=(/-1.d100,1.d0/),comment=comments)

  comments = 'kernel convergence (low, high) - directmin only'
  call input_var(dummy_real,'1.d-5',dict//LIN_KERNEL//GNRM_CV_COEFF//0,ranges=(/0.d0,1.d0/))
  call input_var(dummy_real,'1.d-5',dict//LIN_KERNEL//GNRM_CV_COEFF//1,ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'density convergence (low, high)'
  call input_var(dummy_real,'1.d-13',dict//LIN_KERNEL//RPNRM_CV//0,ranges=(/0.d0,1.d0/))
  call input_var(dummy_real,'1.d-13',dict//LIN_KERNEL//RPNRM_CV//1,ranges=(/0.d0,1.d0/),comment=comments)

  !these variables seems deprecated, put them to their default value
  comments = 'convergence criterion on density to fix TMBS'
  call input_var(dummy_real,'1.d-10',ranges=(/0.d0,1.d0/),comment=comments)

  ! Miscellaneous
  comments='mixing method: 100 (direct minimization), 101 (simple dens mixing), 102 (simple pot mixing), 103 (FOE)'
  call input_var(dummy_int,'100',ranges=(/100,103/),comment=comments)
  select case(dummy_int)
  case(100)
     call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'DIRMIN')
  case(101) 
     call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'DIAG')
     call dict_set(dict//LIN_KERNEL//MIXING_METHOD,'DEN')
  case(102)      
     call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'DIAG')
     call dict_set(dict//LIN_KERNEL//MIXING_METHOD,'POT')
  case(103)
     call dict_set(dict//LIN_KERNEL//LINEAR_METHOD,'FOE')
  end select

  comments = 'initial step size for basis optimization (DIIS, SD)' ! DELETE ONE
  call input_var(dummy_real,'1.d0',dict//LIN_BASIS//ALPHA_DIIS,ranges=(/0.0_gp,10.0_gp/))
  call input_var(dummy_real,'1.d0',dict//LIN_BASIS//ALPHA_SD,ranges=(/0.0_gp,10.0_gp/),comment=comments)

  comments = 'initial step size for kernel update (SD), curve fitting for alpha update - directmin only'
  call input_var(dummy_real,'1.d0',dict//LIN_KERNEL//ALPHA_SD_COEFF,ranges=(/0.0_gp,10.0_gp/))
  call input_var(dummy_bool,'F',dict//LIN_KERNEL//ALPHA_FIT_COEFF,comment=comments)

  comments = 'lower and upper bound for the eigenvalue spectrum (FOE). Will be adjusted automatically if chosen too small'
  call input_var(dummy_real,'-.5d0',dict//LIN_KERNEL//EVAL_RANGE_FOE//0,ranges=(/-10.d0,-1.d-10/))
  call input_var(dummy_real,'-.5d0',dict//LIN_KERNEL//EVAL_RANGE_FOE//1,ranges=(/1.d-10,10.d0/),comment=comments)

  comments='number of iterations in the preconditioner, order of Taylor approximations'
  call input_var(dummy_int,'5',dict//LIN_BASIS//NSTEP_PREC,ranges=(/1,100/))
  call input_var(in%lin%order_taylor,'1',ranges=(/1,100/),comment=comments)
  
  comments = '0-> exact Loewdin, 1-> taylor expansion; &
             &in orthoconstraint: correction for non-orthogonality (0) or no correction (1)'
  call input_var(dummy_int,'1',dict//LIN_GENERAL//TAYLOR_ORDER,ranges=(/-100,10000/))
  call input_var(dummy_int,'1',dict//LIN_BASIS//CORRECTION_ORTHOCONSTRAINT,comment=comments)

  !!this variable seems deprecated
  !call input_var(dummy_int,'1',ranges=(/0,1/),comment=comments)

  comments='fscale: length scale over which complementary error function decays from 1 to 0'
  call input_var(dummy_real,'1.d-2',dict//LIN_KERNEL//FSCALE_FOE,ranges=(/0.d0,1.d0/),comment=comments)

  !plot basis functions: true or false
  comments='Output basis functions: 0 no output, 1 formatted output, 2 Fortran bin, 3 ETSF ;'//&
           'calculate dipole ; pulay correction (old and new); diagonalization at the end (dmin, FOE)'
  call input_var(dummy_int,'0',dict//LIN_GENERAL//OUTPUT_WF,ranges=(/0,3/))
  call input_var(dummy_bool,'F',dict//LIN_GENERAL//CALC_DIPOLE)
  call input_var(dummy_bool,'T',dict//LIN_GENERAL//CALC_PULAY//0)
  call input_var(dummy_bool,'F',dict//LIN_GENERAL//CALC_PULAY//1)

!  in%lin%pulay_correction=dummy_bool
!  call input_var(in%lin%new_pulay_correction,'F')
  call input_var(dummy_bool,'F',dict//LIN_GENERAL//SUBSPACE_DIAG,comment=comments)



  !!!!Came to here for the definition of the input variables. Fill the input variable structure with what found so far
  !filling of the input variable (to be moved)

  !list of variables declared as deprecated, which most likely have to be removed of filled with dummy values
  !in%lin%correctionOrthoconstraint
  !in%lin%support_functions_converged
  !in%lin%deltaenergy_multiplier_TMBexit
  !in%lin%deltaenergy_multiplier_TMBfix


  !fragment calculation and transfer integrals: true or false
  comments='fragment calculation; calculate transfer_integrals; constrained DFT calculation; extra states to optimize (dmin only)'
  !these should becode dummy variables to build dictionary
  call input_var(in%lin%fragment_calculation,'F')
  call input_var(in%lin%calc_transfer_integrals,'F')
  call input_var(in%lin%constrained_dft,'F')
  call input_var(in%lin%extra_states,'0',ranges=(/0,10000/),comment=comments)

  ! Now read in the parameters specific for each atom type.
  comments = 'Atom name, number of basis functions per atom, prefactor for confinement potential,'//&
       'localization radius, kernel cutoff, kernel cutoff FOE'
  !itype = 1
  read_basis: do !while(itype <= atoms%astruct%ntypes) 
     if (exists) then
        call input_var(atomname,'C',input_iostat=ios)
        if (ios /= 0) exit read_basis
        dict_basis=>dict//LIN_BASIS_PARAMS//trim(atomname)
     else
        call input_var(atomname,'C')!trim(atoms%astruct%atomnames(1)))
        dict_basis=>dict//LIN_BASIS_PARAMS! default values
        !itype = itype + 1
     end if

     !number of basis functions for this atom type
     call input_var(npt,'1',dict_basis//NBASIS,ranges=(/1,100/),input_iostat=ios)
     call input_var(ppao,'1.2d-2',dict_basis//AO_CONFINEMENT,&
          ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(ppl,'1.2d-2',dict_basis//CONFINEMENT//0,&
          ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(pph,'5.d-5',dict_basis//CONFINEMENT//1,&
          ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(lrl,'10.d0',dict_basis//RLOC//0,&
          ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(lrh,'10.d0',dict_basis//RLOC//1,&
          ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(kco,'12.d0',dict_basis//RLOC_KERNEL,&
          ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(kco_FOE,'20.d0',dict_basis//RLOC_KERNEL_FOE,&
          ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios,comment=comments)

     if (.not. exists) exit read_basis !default has been filled
  end do read_basis

  call input_free((iproc == 0) .and. dump)

  !input dictionary is ready
  !if (iproc==0) call yaml_map('Dictionary for lin',dict)

!!  dummy_bool=dict//LIN_GENERAL//HYBRID
!!  if (dummy_bool) then
!!     in%lin%nlevel_accuracy=1
!!  else
!!     in%lin%nlevel_accuracy=2
!!  end if
!!  dummy_iarr=dict//LIN_GENERAL//NIT
!!  in%lin%nit_lowaccuracy=dummy_iarr(1)
!!  in%lin%nit_highaccuracy=dummy_iarr(2)
!!  dummy_darr=dict//LIN_GENERAL//RPNRM_CV
!!  in%lin%lowaccuracy_conv_crit =dummy_darr(1)
!!  in%lin%highaccuracy_conv_crit=dummy_darr(2)
!!  in%lin%reduce_confinement_factor=dict//LIN_GENERAL//CONF_DAMPING
!!  in%lin%methTransformOverlap=dict//LIN_GENERAL//TAYLOR_ORDER
!!  in%lin%plotBasisFunctions=dict//LIN_GENERAL//OUTPUT_WF
!!  in%lin%calc_dipole=dict//LIN_GENERAL//CALC_DIPOLE
!!  in%lin%pulay_correction=dict//LIN_GENERAL//CALC_PULAY//0
!!  in%lin%new_pulay_correction=dict//LIN_GENERAL//CALC_PULAY//1
!!  in%lin%diag_end=dict//LIN_GENERAL//SUBSPACE_DIAG
!!
!!
!!  dummy_iarr=dict//LIN_BASIS//NIT
!!  in%lin%nItBasis_lowaccuracy =dummy_iarr(1) 
!!  in%lin%nItBasis_highaccuracy=dummy_iarr(2)
!!  dummy_iarr=dict//LIN_BASIS//IDSX
!!  in%lin%DIIS_hist_lowaccur =dummy_iarr(1)
!!  in%lin%DIIS_hist_highaccur=dummy_iarr(2)
!!  dummy_darr=dict//LIN_BASIS//GNRM_CV
!!  in%lin%convCrit_lowaccuracy =dummy_darr(1)
!!  in%lin%convCrit_highaccuracy=dummy_darr(2)
!!  in%lin%early_stop=dict//LIN_BASIS//DELTAE_CV
!!  in%lin%gnrm_dynamic=dict//LIN_BASIS//GNRM_DYN
!!  in%lin%alphaDIIS=dict//LIN_BASIS//ALPHA_DIIS
!!  in%lin%alphaSD=dict//LIN_BASIS//ALPHA_SD
!!  in%lin%nItPrecond=dict//LIN_BASIS//NSTEP_PREC
!!  in%lin%support_functions_converged=dict//LIN_BASIS//FIX_BASIS
!!
!!  !filling of input variable
!!  dummy_char=dict//LIN_KERNEL//LINEAR_METHOD
!!  select case(trim(dummy_char))
!!  case('DIRMIN')
!!     in%lin%scf_mode=100
!!  case('DIAG')
!!     dummy_char2=dict//LIN_KERNEL//MIXING_METHOD
!!     select case(trim(dummy_char2))
!!     case('DEN')
!!        in%lin%scf_mode=101
!!     case('POT')
!!        in%lin%scf_mode=102
!!     end select
!!  case('FOE')
!!     in%lin%scf_mode=103
!!  end select
!!  dummy_iarr=dict//LIN_KERNEL//NIT
!!  in%lin%nItSCCWhenFixed_lowaccuracy =dummy_iarr(1)
!!  in%lin%nItSCCWhenFixed_highaccuracy=dummy_iarr(2)
!!  dummy_iarr=dict//LIN_KERNEL//NSTEP
!!  in%lin%nItdmin_lowaccuracy =dummy_iarr(1)
!!  in%lin%nItdmin_highaccuracy=dummy_iarr(2)
!!  dummy_iarr=dict//LIN_KERNEL//IDSX_COEFF
!!  in%lin%dmin_hist_lowaccuracy =dummy_iarr(1)
!!  in%lin%dmin_hist_highaccuracy=dummy_iarr(2)
!!  dummy_iarr=dict//LIN_KERNEL//IDSX
!!  in%lin%mixHist_lowaccuracy =dummy_iarr(1)
!!  in%lin%mixHist_highaccuracy=dummy_iarr(2)
!!  dummy_darr=dict//LIN_KERNEL//ALPHAMIX
!!  in%lin%alpha_mix_lowaccuracy =dummy_darr(1)
!!  in%lin%alpha_mix_highaccuracy=dummy_darr(2)
!!  dummy_darr=dict//LIN_KERNEL//GNRM_CV_COEFF
!!  in%lin%convCritdmin_lowaccuracy =dummy_darr(1)
!!  in%lin%convCritdmin_highaccuracy=dummy_darr(2)
!!  dummy_darr=dict//LIN_KERNEL//RPNRM_CV
!!  in%lin%convCritMix_lowaccuracy =dummy_darr(1)
!!  in%lin%convCritMix_highaccuracy=dummy_darr(2)
!!  in%lin%alphaSD_coeff=dict//LIN_KERNEL//ALPHA_SD_COEFF
!!  in%lin%curvefit_dmin=dict//LIN_KERNEL//ALPHA_FIT_COEFF
!!  dummy_darr=dict//LIN_KERNEL//EVAL_RANGE_FOE
!!  in%lin%evlow =dummy_darr(1)
!!  in%lin%evhigh=dummy_darr(2)
!!  in%lin%fscale=dict//LIN_KERNEL//FSCALE_FOE

  !in%lin%correctionOrthoconstraint=1 !to be checked later

  ! not sure whether to actually make this an input variable or not so just set to false for now
  in%lin%diag_start=.false.

  ! It is not possible to use both the old and the new Pulay correction at the same time
  if (in%lin%pulay_correction .and. in%lin%new_pulay_correction) then
     stop 'It is not possible to use both the old and the new Pulay correction at the same time!'
  end if

  ! Allocate lin pointers and atoms%rloc

!  call nullifyInputLinparameters(in%lin)
!  call allocateBasicArraysInputLin(in%lin, atoms%astruct%ntypes)

  !!!first fill all the types by the default, then override by per-type values
  !!do jtype=1,atoms%astruct%ntypes
  !!   dict_basis => dict_iter(dict//BASIS_PARAMS)
  !!   do while(associated(dict_basis))
  !!      call basis_params_set_dict(dict_basis,in%lin,jtype)
  !!      dict_basis => dict_next(dict_basis)
  !!   end do
  !!   !then check if the objects exists in separate specifications
  !!   if (trim(atoms%astruct%atomnames(jtype)) .in. dict//BASIS_PARAMS) then
  !!      dict_basis => &
  !!           dict_iter(dict//BASIS_PARAMS//trim(atoms%astruct%atomnames(jtype)))
  !!   end if
  !!   do while(associated(dict_basis))
  !!      call basis_params_set_dict(dict_basis,in%lin,jtype)
  !!      dict_basis => dict_next(dict_basis)
  !!   end do
  !!end do
  !then perform extra allocations
  nlr=0
  do iat=1,atoms%astruct%nat
      itype=atoms%astruct%iatype(iat)
      nlr=nlr+in%lin%norbsPerType(itype)
  end do

  allocate(in%lin%locrad(nlr),stat=istat)
  call memocc(istat,in%lin%locrad,'in%lin%locrad',subname)

  allocate(in%lin%locrad_kernel(nlr),stat=istat)
  call memocc(istat,in%lin%locrad_kernel,'in%lin%locrad_kernel',subname)

  allocate(in%lin%locrad_mult(nlr),stat=istat)
  call memocc(istat,in%lin%locrad_mult,'in%lin%locrad_mult',subname)

  allocate(in%lin%locrad_lowaccuracy(nlr),stat=istat)
  call memocc(istat,in%lin%locrad_lowaccuracy,'in%lin%locrad_lowaccuracy',subname)

  allocate(in%lin%locrad_highaccuracy(nlr),stat=istat)
  call memocc(istat,in%lin%locrad_highaccuracy,'in%lin%locrad_highaccuracy',subname)

  ! Assign the localization radius to each atom.
  iiorb=0
  do iat=1,atoms%astruct%nat
      itype=atoms%astruct%iatype(iat)
      do iorb=1,in%lin%norbsPerType(itype)
          iiorb=iiorb+1
          in%lin%locrad(iiorb)=in%lin%locrad_type(itype,1)
          in%lin%locrad_kernel(iiorb)=in%lin%kernel_cutoff(itype)
          in%lin%locrad_mult(iiorb)=in%lin%kernel_cutoff_FOE(itype)
          in%lin%locrad_lowaccuracy(iiorb)=in%lin%locrad_type(itype,1) 
          !locradType_lowaccur(itype)
          in%lin%locrad_highaccuracy(iiorb)=in%lin%locrad_type(itype,2)
          !locradType_highaccur(itype)
      end do
  end do
!!$  if (.not. exists) then
!!$     call dict_free(dict)
!!$     call input_free(.false.)
!!$     return
!!$  end if

  call dict_free(dict)

END SUBROUTINE lin_input_variables_new

!> Read fragment input parameters
subroutine fragment_input_variables(iproc,dump,filename,in,atoms)
  use module_base
  use module_types
  use module_input
  implicit none
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: filename
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  logical, intent(in) :: dump
  !local variables
  !character(len=*), parameter :: subname='fragment_input_variables'
  logical :: exists
  character(len=256) :: comments
  integer :: ifrag, frag_num, ierr
  real(gp) :: charge
  type(dictionary), pointer :: dict

  !Linear input parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Fragment Parameters') 
  if (exists .and. dump) in%files = in%files + INPUTS_FRAG

  if (.not. exists .and. in%lin%fragment_calculation) then ! we should be doing a fragment calculation, so this is a problem
     write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' is missing and fragment calculation was specified"
     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
     stop
  end if

  !to be continued after having fixed linear variables
!!$  !iteration over the dictionary
!!$  !count the number of reference fragments
!!$  dict_tmp=>dict_iter(dict)
!!$  do while(associated(dict_tmp))
!!$     select case(trim(dict_key(dict_tmp)))
!!$     case(TRANSFER_INTEGRALS)
!!$        frag%calc_transfer_integrals=dict_tmp
!!$     case(CONSTRAINED_DFT)
!!$        frag%constrained_dft=dict_tmp
!!$        ncharged=dict_size(dict_tmp)
!!$     case default
!!$        frag%nfrag_ref=frag%nfrag_ref+1
!!$        !count the number of fragments for this reference
!!$        call count_local_fragments(dict_tmp,ncount)
!!$        frag%nfrag=frag%nfrag+ncount
!!$     end select
!!$     dict_tmp=>dict_next(dict_tmp)
!!$  end do
!!$
!!$  call allocateInputFragArrays(frag)
!!$
!!$  dict_tmp=>dict_iter(dict)
!!$  frag_num=0
!!$  do while(associated(dict_tmp))
!!$     select case(trim(dict_key(dict_tmp)))
!!$     case(TRANSFER_INTEGRALS)
!!$     case(CONSTRAINED_DFT)
!!$        !iterate over the charge
!!$     case default
!!$        frag_num=frag_num+1
!!$        frag%label(frag_num)=repeat(' ',len(frag%label(frag_num)))
!!$        frag%label(frag_num)=trim(dict_key(dict_tmp))
!!$        !update directory name
!!$        in%frag%dirname(frag_num)='data-'//trim(in%frag%label(frag_num))//'/'
!!$        call count_local_fragments(dict_tmp,icount,frag_index=frag%frag_index,frag_id=frag_num)
!!$     end select
!!$     dict_tmp=>dict_next(dict_tmp)
!!$  end do
!!$
!!$  subroutine count_local_fragments(dict_tmp,icount,frag_index,frag_id)
!!$    implicit none
!!$    integer, intent(out) :: icount
!!$    type(dictionary), pointer :: dict_tmp
!!$    integer, dimension(:), intent(inout), optional :: frag_index
!!$    integer, intent(in), optional :: frag_id
!!$    
!!$    !local variables
!!$    integer :: idum,istart
!!$    type(dictionary), pointer :: d_tmp
!!$    character(len=max_field_length) :: val
!!$
!!$    !iteration over the whole list
!!$    icount=0
!!$    istart=0
!!$    d_tmp=>dict_iter(dict_tmp)
!!$    do while(associated(d_tmp))
!!$       val=d_tmp
!!$       !if string is a integer consider it
!!$       if (is_atoi(val)) then
!!$          read(val,*) idum
!!$          if (f_err_raise(istart==idum,'error in entering fragment ids',&
!!$               err_name='BIGDFT_INPUT_VARIABLES_ERROR')) return
!!$          if (istart /=0) then
!!$             icount=idum-istart
!!$             do i=istart,idum
!!$                frag_index(i)=frag_id
!!$             end do
!!$             istart=0
!!$          else
!!$             icount=icount+1
!!$             frag_index(idum)=frag_id
!!$          end if
!!$       else if (f_err_raise(adjustl(trim(val))/='...',&
!!$            'the only allowed values in the fragment list are integers or "..." string',&
!!$            err_name='BIGDFT_INPUT_VARIABLES_ERROR')) then
!!$          return
!!$       else
!!$          istart=idum
!!$       end if
!!$       d_tmp=>dict_next(d_tmp)
!!$    end do
!!$  end subroutine count_local_fragments
!!$
!!$
!!$  call dict_init(dict)

  
  ! number of reference fragments
  comments='# number of fragments in reference system, number of fragments in current system'
  call input_var(in%frag%nfrag_ref,'1',ranges=(/1,100000/))
  call input_var(in%frag%nfrag,'1',ranges=(/1,100000/),comment=comments)
  
  ! Allocate fragment pointers
  call nullifyInputFragParameters(in%frag)
  call allocateInputFragArrays(in%frag)

  ! ADD A SENSIBLE DEFAULT AND ALLOW FOR USER NOT TO SPECIFY FRAGMENT NAMES
  comments = '#  reference fragment number i, fragment label'
  do ifrag=1,in%frag%nfrag_ref
    call input_var(frag_num,'1',ranges=(/1,in%frag%nfrag_ref/))
    if (frag_num/=ifrag) then
        write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' has an error when specifying&
             & the reference fragments"
       call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
       stop
    end if
    call input_var(in%frag%label(frag_num),' ',comment=comments)
    in%frag%label(frag_num)=trim(in%frag%label(frag_num))
    ! keep dirname blank if this isn't a fragment calculation
    if (len(trim(in%frag%label(frag_num)))>=1) then
       in%frag%dirname(frag_num)='data-'//trim(in%frag%label(frag_num))//'/'
    else
       in%frag%dirname(frag_num)=''
    end if
  end do

  comments = '# fragment number j, reference fragment i this corresponds to, charge on this fragment'
  do ifrag=1,in%frag%nfrag
    call input_var(frag_num,'1',ranges=(/1,in%frag%nfrag/))
    if (frag_num/=ifrag) then
        write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' has an error when specifying&
             & the system fragments"
       call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
       stop
    end if
    call input_var(in%frag%frag_index(frag_num),'1',ranges=(/0,100000/))
    call input_var(charge,'0.d0',ranges=(/-500.d0,500.d0/),comment=comments)
    in%frag%charge(frag_num)=charge
    !call input_var(in%frag%charge(frag_num),'1',ranges=(/-500,500/),comment=comments)
  end do

  call input_free((iproc == 0) .and. dump)

END SUBROUTINE fragment_input_variables
