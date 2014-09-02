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
  implicit none
  type(input_variables), intent(inout) :: in
  character(len = *), intent(in) :: radical

  call set_inputfile(in%file_lin, radical,    "lin")
  call set_inputfile(in%file_frag, radical,   "frag")

  !in%files = INPUTS_NONE
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

subroutine fragment_variables_from_old_text_format(in, run_name)
  use dictionaries
  use module_base, only: bigdft_mpi
  use module_types
  use input_old_text_format
  implicit none
  type(input_variables), intent(inout) :: in
  character(len = *), intent(in) :: run_name
  type(dictionary), pointer :: dict_frag
  integer :: ierr

  call set_inputfile(in%file_lin, run_name,    "lin")
  call set_inputfile(in%file_frag, run_name,   "frag")

  ! To avoid race conditions where procs create the default file and other test its
  ! presence, we put a barrier here.
  if (bigdft_mpi%nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)

  ! Linear scaling (if given)
  !in%lin%fragment_calculation=.false. ! to make sure that if we're not doing a linear calculation we don't read fragment information
  call fragment_input_variables_check(bigdft_mpi%iproc,in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR, trim(in%file_lin),in%lin)

  ! Fragment information (if given)
  nullify(dict_frag)
  call fragment_input_variables_from_text_format(bigdft_mpi%iproc,(in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR).and.in%lin%fragment_calculation,&
       trim(in%file_frag),in%lin%fragment_calculation,dict_frag)

  call frag_from_dict(dict_frag,in%frag)  
  call dict_free(dict_frag)
END SUBROUTINE fragment_variables_from_old_text_format

!> Read fragment input parameters from old format
subroutine fragment_input_variables_check(iproc,dump,filename,lin)
  use module_base, dictionary_value_lgt => max_field_length
  use module_types
  use module_input
  use yaml_output, only: yaml_map
  use module_input_keys
  implicit none
  integer, intent(in) :: iproc
  !integer, intent(inout) :: files
  character(len=*), intent(in) :: filename
  type(linearInputParameters), intent(inout) :: lin
  !type(input_variables), intent(inout) :: in
  !type(atoms_data), intent(in) :: atoms
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
  !if (exists) files = files + INPUTS_LIN

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

  !comments='number of iterations in the preconditioner, order of Taylor approximations'
  comments='number of iterations in the preconditioner'
  call input_var(dummy_int,'5',dict//LIN_BASIS//NSTEP_PREC,ranges=(/1,100/),comment=comments)
  !call input_var(in%lin%order_taylor,'1',ranges=(/1,100/),comment=comments)
  
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

  !fragment calculation and transfer integrals: true or false
  comments='fragment calculation; calculate transfer_integrals; constrained DFT calculation; extra states to optimize (dmin only)'
  !these should becode dummy variables to build dictionary
  call input_var(lin%fragment_calculation,'F')
  call input_var(lin%calc_transfer_integrals,'F')
  call input_var(lin%constrained_dft,'F')
  call input_var(lin%extra_states,'0',ranges=(/0,10000/),comment=comments)

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

  call dict_free(dict)

  ! not sure whether to actually make this an input variable or not so just set to false for now
  lin%diag_start=.false.

  ! It is not possible to use both the old and the new Pulay correction at the same time
  if (lin%pulay_correction .and. lin%new_pulay_correction) then
     stop 'It is not possible to use both the old and the new Pulay correction at the same time!'
  end if

END SUBROUTINE fragment_input_variables_check

