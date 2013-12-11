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

  !set prefix name of the run (as input by defaut for input.dft)
  in%run_name=repeat(' ',len(in%run_name))
  if (trim(radical) /= 'input') in%run_name=trim(radical)

  call set_inputfile(in%file_occnum, radical, "occ")
  call set_inputfile(in%file_igpop, radical,  "occup")
  call set_inputfile(in%file_lin, radical,    "lin")
  call set_inputfile(in%file_frag, radical,   "frag")

  if (trim(radical) == "input") then
     in%dir_output="data" // trim(bigdft_run_id_toa())
  else
     in%dir_output="data-"//trim(radical)!//trim(bigdft_run_id_toa())
  end if

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

!> Read linear input parameters
subroutine lin_input_variables_new(iproc,dump,filename,in,atoms)
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
  logical :: exists
  character(len=*), parameter :: subname='lin_input_variables'
  character(len=256) :: comments
  logical,dimension(atoms%astruct%ntypes) :: parametersSpecified
  logical :: found
  character(len=20):: atomname
  integer :: itype, jtype, ios, ierr, iat, npt, iiorb, iorb, nlr, istat
  real(gp):: ppao, ppl, pph, lrl, lrh, kco
  real(gp),dimension(atoms%astruct%ntypes) :: locradType, locradType_lowaccur, locradType_highaccur
  type(dictionary), pointer :: dict

  !this variable is temporay for the moment, it should integrate the input dictionary
  call dict_init(dict)

  !Linear input parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Linear Parameters')  
  if (exists) in%files = in%files + INPUTS_LIN

  ! number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)
  comments='number of accuracy levels: either 2 (for low/high accuracy) or 1 (for hybrid mode)'
  call input_var(in%lin%nlevel_accuracy,'2',dict//'nlevel_accuracy',ranges=(/1,2/),comment=comments)

  ! number of iterations
  comments = 'outer loop iterations (low, high)'
  call input_var(in%lin%nit_lowaccuracy,'15',dict//'nit_lowaccuracy',ranges=(/0,100000/))
  call input_var(in%lin%nit_highaccuracy,'1',dict//'nit_highaccuracy',ranges=(/0,100000/),comment=comments)

  comments = 'basis iterations (low, high)'
  call input_var(in%lin%nItBasis_lowaccuracy,'12',ranges=(/0,100000/))
  call input_var(in%lin%nItBasis_highaccuracy,'50',ranges=(/0,100000/),comment=comments)

  comments = 'kernel iterations (low, high) - directmin only'
  call input_var(in%lin%nItdmin_lowaccuracy,'1',ranges=(/0,1000/))
  call input_var(in%lin%nItdmin_highaccuracy,'1',ranges=(/0,1000/),comment=comments)

  comments = 'density iterations (low, high)'
  call input_var(in%lin%nItSCCWhenFixed_lowaccuracy,'15',ranges=(/0,1000/))
  call input_var(in%lin%nItSCCWhenFixed_highaccuracy,'15',ranges=(/0,1000/),comment=comments)

  ! DIIS history lengths
  comments = 'DIIS history for basis (low, high)'
  call input_var(in%lin%DIIS_hist_lowaccur,'5',ranges=(/0,100/))
  call input_var(in%lin%DIIS_hist_highaccur,'0',ranges=(/0,100/),comment=comments)

  comments = 'DIIS history for kernel (low, high) - directmin only'
  call input_var(in%lin%dmin_hist_lowaccuracy,'0',ranges=(/0,100/))
  call input_var(in%lin%dmin_hist_highaccuracy,'0',ranges=(/0,100/),comment=comments)

  comments = 'DIIS history for density mixing (low, high)'
  call input_var(in%lin%mixHist_lowaccuracy,'0',ranges=(/0,100/))
  call input_var(in%lin%mixHist_highaccuracy,'0',ranges=(/0,100/),comment=comments)

  ! mixing parameters
  comments = 'density mixing parameter (low, high)'
  call input_var(in%lin%alpha_mix_lowaccuracy,'.5d0',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%alpha_mix_highaccuracy,'.5d0',ranges=(/0.d0,1.d0/),comment=comments)

  ! Convergence criteria
  comments = 'outer loop convergence (low, high)'
  call input_var(in%lin%lowaccuracy_conv_crit,'1.d-8',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%highaccuracy_conv_crit,'1.d-12',ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'basis convergence (low, high) ; early stop TMB optimization (experimental mode only)'
  call input_var(in%lin%convCrit_lowaccuracy,'1.d-3',ranges=(/0.0_gp,1.0_gp/))
  call input_var(in%lin%convCrit_highaccuracy,'1.d-5',ranges=(/0.0_gp,1.0_gp/))
  call input_var(in%lin%early_stop,'1.d-4',ranges=(/0.0_gp,1.0_gp/),comment=comments)
  
  comments = 'multiplier to (exit one TMB optimization, fix TMB completely). Only used for hybrid mode'
  call input_var(in%lin%deltaenergy_multiplier_TMBexit,'1.d0',ranges=(/1.d-5,1.d1/))
  call input_var(in%lin%deltaenergy_multiplier_TMBfix,'1.d0',ranges=(/1.d-5,1.d1/),comment=comments)

  comments = 'factor to reduce the confinement. Only used for hybrid mode.'
  call input_var(in%lin%reduce_confinement_factor,'0.5d0',ranges=(/-1.d100,1.d0/),comment=comments)

  comments = 'kernel convergence (low, high) - directmin only'
  call input_var(in%lin%convCritdmin_lowaccuracy,'1.d-5',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%convCritdmin_highaccuracy,'1.d-5',ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'density convergence (low, high)'
  call input_var(in%lin%convCritMix_lowaccuracy,'1.d-13',ranges=(/0.d0,1.d0/))
  call input_var(in%lin%convCritMix_highaccuracy,'1.d-13',ranges=(/0.d0,1.d0/),comment=comments)

  comments = 'convergence criterion on density to fix TMBS'
  call input_var(in%lin%support_functions_converged,'1.d-10',ranges=(/0.d0,1.d0/),comment=comments)

  ! Miscellaneous
  comments='mixing method: 100 (direct minimization), 101 (simple dens mixing), 102 (simple pot mixing), 103 (FOE)'
  call input_var(in%lin%scf_mode,'100',ranges=(/100,103/),comment=comments)

  comments = 'initial step size for basis optimization (DIIS, SD)' ! DELETE ONE
  call input_var(in%lin%alphaDIIS,'1.d0',ranges=(/0.0_gp,10.0_gp/))
  call input_var(in%lin%alphaSD,'1.d0',ranges=(/0.0_gp,10.0_gp/),comment=comments)

  comments = 'initial step size for kernel update (SD), curve fitting for alpha update - directmin only'
  call input_var(in%lin%alphaSD_coeff,'1.d0',ranges=(/0.0_gp,10.0_gp/))
  call input_var(in%lin%curvefit_dmin,'F',comment=comments)

  comments = 'lower and upper bound for the eigenvalue spectrum (FOE). Will be adjusted automatically if chosen too small'
  call input_var(in%lin%evlow,'-.5d0',ranges=(/-10.d0,-1.d-10/))
  call input_var(in%lin%evhigh,'-.5d0',ranges=(/1.d-10,10.d0/),comment=comments)

  comments='number of iterations in the preconditioner, order of Taylor approximations'
  call input_var(in%lin%nItPrecond,'5',ranges=(/1,100/))
  call input_var(in%lin%order_taylor,'1',ranges=(/1,3/),comment=comments)
  
  comments = '0-> exact Loewdin, 1-> taylor expansion; &
             &in orthoconstraint: correction for non-orthogonality (0) or no correction (1)'
  call input_var(in%lin%methTransformOverlap,'1',ranges=(/-1,1/))
  call input_var(in%lin%correctionOrthoconstraint,'1',ranges=(/0,1/),comment=comments)

  comments='fscale: length scale over which complementary error function decays from 1 to 0'
  call input_var(in%lin%fscale,'1.d-2',ranges=(/0.d0,1.d0/),comment=comments)

  !plot basis functions: true or false
  comments='Output basis functions: 0 no output, 1 formatted output, 2 Fortran bin, 3 ETSF ;'//&
           'calculate dipole ; pulay correction; diagonalization at the end (dmin, FOE)'
  call input_var(in%lin%plotBasisFunctions,'0',ranges=(/0,3/))
  call input_var(in%lin%calc_dipole,'F')
  call input_var(in%lin%pulay_correction,'T')
  call input_var(in%lin%diag_end,'F',comment=comments)
  ! not sure whether to actually make this an input variable or not so just set to false for now
  in%lin%diag_start=.false.


  !fragment calculation and transfer integrals: true or false
  comments='fragment calculation; calculate transfer_integrals; constrained DFT calculation; extra states to optimize (dmin only)'
  call input_var(in%lin%fragment_calculation,'F')
  call input_var(in%lin%calc_transfer_integrals,'F')
  call input_var(in%lin%constrained_dft,'F')
  call input_var(in%lin%extra_states,'0',ranges=(/0,10000/),comment=comments)

  ! Allocate lin pointers and atoms%rloc
  call nullifyInputLinparameters(in%lin)
  call allocateBasicArraysInputLin(in%lin, atoms%astruct%ntypes)
  
  ! Now read in the parameters specific for each atom type.
  comments = 'Atom name, number of basis functions per atom, prefactor for confinement potential,'//&
             'localization radius, kernel cutoff'
  parametersSpecified=.false.
  itype = 1
  do
     !Check at the beginning to permit natom=0
     if (itype > atoms%astruct%ntypes) exit
     if (exists) then
        call input_var(atomname,'C',input_iostat=ios)
        if (ios /= 0) exit
     else
        call input_var(atomname,trim(atoms%astruct%atomnames(itype)))
        itype = itype + 1
     end if
     call input_var(npt,'1',ranges=(/1,100/),input_iostat=ios)
     call input_var(ppao,'1.2d-2',ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(ppl,'1.2d-2',ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(pph,'5.d-5',ranges=(/0.0_gp,1.0_gp/),input_iostat=ios)
     call input_var(lrl,'10.d0',ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(lrh,'10.d0',ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios)
     call input_var(kco,'20.d0',ranges=(/1.0_gp,10000.0_gp/),input_iostat=ios,comment=comments)
     ! The reading was successful. Check whether this atom type is actually present.
     found=.false.
     do jtype=1,atoms%astruct%ntypes
        if(trim(atomname)==trim(atoms%astruct%atomnames(jtype))) then
           found=.true.
           parametersSpecified(jtype)=.true.
           in%lin%norbsPerType(jtype)=npt
           in%lin%potentialPrefac_ao(jtype)=ppao
           in%lin%potentialPrefac_lowaccuracy(jtype)=ppl
           in%lin%potentialPrefac_highaccuracy(jtype)=pph
           locradType(jtype)=lrl
           in%lin%locrad_type(jtype)=lrl
           locradType_lowaccur(jtype)=lrl
           locradType_highaccur(jtype)=lrh
           atoms%rloc(jtype,:)=locradType(jtype)
           in%lin%kernel_cutoff(jtype)=kco
        end if
     end do
     if(.not.found) then
        if(iproc==0 .and. dump) write(*,'(1x,3a)') "WARNING: you specified informations about the atomtype '",trim(atomname), &
             "', which is not present in the file containing the atomic coordinates."
     end if
  end do
  found  = .true.
  do jtype=1,atoms%astruct%ntypes
     found = found .and. parametersSpecified(jtype)
  end do
  if (.not. found) then
     ! The parameters were not specified for all atom types.
     if(iproc==0) then
        write(*,'(1x,a)',advance='no') "ERROR: the file 'input.lin' does not contain the parameters&
             & for the following atom types:"
        do jtype=1,atoms%astruct%ntypes
           if(.not.parametersSpecified(jtype)) write(*,'(1x,a)',advance='no') trim(atoms%astruct%atomnames(jtype))
        end do
     end if
     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
     stop
  end if

  nlr=0
  do iat=1,atoms%astruct%nat
      itype=atoms%astruct%iatype(iat)
      nlr=nlr+in%lin%norbsPerType(itype)
  end do
  allocate(in%lin%locrad(nlr),stat=istat)
  call memocc(istat,in%lin%locrad,'in%lin%locrad',subname)
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
          in%lin%locrad(iiorb)=locradType(itype)
          in%lin%locrad_lowaccuracy(iiorb)=locradType_lowaccur(itype)
          in%lin%locrad_highaccuracy(iiorb)=locradType_highaccur(itype)
      end do
  end do
  

  call dict_free(dict)

  call input_free((iproc == 0) .and. dump)
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

  !Linear input parameters
  call input_set_file(iproc,dump,trim(filename),exists,'Fragment Parameters') 
  if (exists .and. dump) in%files = in%files + INPUTS_FRAG

  if (.not. exists .and. in%lin%fragment_calculation) then ! we should be doing a fragment calculation, so this is a problem
     write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' is missing and fragment calculation was specified"
     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
     stop
  end if

  ! number of reference fragments
  comments='# number of fragments in reference system, number of fragments in current system'
  call input_var(in%frag%nfrag_ref,'1',ranges=(/1,100000/))
  call input_var(in%frag%nfrag,'1',ranges=(/1,100000/),comment=comments)
  
  ! Allocate fragment pointers
  call nullifyInputFragParameters(in%frag)
  call allocateInputFragArrays(in%frag)

  !comments = '# reference fragment number i, number of atoms in reference fragment i, '//&
  !           'number of atoms in corresponding environment'
  !do ifrag=1,in%frag%nfrag_ref
  !  call input_var(frag_num,'1',ranges=(/1,in%frag%nfrag_ref/))
  !  if (frag_num/=ifrag) then
  !      write(*,'(1x,a)',advance='no') "ERROR: the file 'input.frag' has an error when specifying&
  !           & the reference fragments"
  !     call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  !     stop
  !  end if
  !  call input_var(in%frag%frag_info(frag_num,1),'1',ranges=(/1,100000/))
  !  call input_var(in%frag%frag_info(frag_num,2),'0',ranges=(/0,100000/),comment=comments)
  !end do

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
