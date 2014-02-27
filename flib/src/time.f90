!> @file
!!  Define routines for timing
!! @author
!!    Copyright (C) 2010-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module containing variables used for the timing for BigDFT
module time_profiling
  use dictionaries
  use dynamic_memory, only: f_time
  implicit none

  !private to be put

  integer, parameter :: ncat_max=144,ncls_max=7   ! define timimg categories and classes
  character(len=14), dimension(ncls_max), parameter :: clss = (/ &
       'Communications'    ,  &
       'Convolutions  '    ,  &
       'Linear Algebra'    ,  &
       'Other         '    ,  &
       'Potential     '    ,  &
       'Initialization'    ,  &
       'Finalization  '    /)
  character(len=14), dimension(3,ncat_max), parameter :: cats = reshape((/ &
       !       Name           Class       Operation Kind
       'ReformatWaves ','Initialization' ,'Small Convol  ' ,  &  !< Reformatting of input waves
       'CrtDescriptors','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of descriptor arrays
       'CrtLocPot     ','Initialization' ,'Miscellaneous ' ,  &  !< Calculation of local potential
       'CrtProjectors ','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of projectors
       'CrtPcProjects ','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of preconditioning projectors
       'CrtPawProjects','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of abscalc-pawprojectors
       'ApplyLocPotKin','Convolutions  ' ,'OpenCL ported ' ,  &  !< Application of PSP, kinetic energy
       'ApplyProj     ','Other         ' ,'RMA pattern   ' ,  &  !< Application of nonlocal PSP
       'Precondition  ','Convolutions  ' ,'OpenCL ported ' ,  &  !< Precondtioning
       'Rho_comput    ','Convolutions  ' ,'OpenCL ported ' ,  &  !< Calculation of charge density (sumrho) computation
       'Rho_commun    ','Communications' ,'AllReduce grid' ,  &  !< Calculation of charge density (sumrho) communication
       'Pot_commun    ','Communications' ,'AllGathrv grid' ,  &  !< Communication of potential
       'Pot_comm start','Communications' ,'MPI_types/_get' ,  &  !< Communication of potential
       'Un-TransSwitch','Other         ' ,'RMA pattern   ' ,  &  !< Transposition of wavefunction, computation
       'Un-TransComm  ','Communications' ,'ALLtoALLV     ' ,  &  !< Transposition of wavefunction, communication
       'GramS_comput  ','Linear Algebra' ,'DPOTRF        ' ,  &  !< Gram Schmidt computation        
       'GramS_commun  ','Communications' ,'ALLReduce orbs' ,  &  !< Gram Schmidt communication
       'LagrM_comput  ','Linear Algebra' ,'DGEMM         ' ,  &  !< Lagrange Multipliers computation
       'LagrM_commun  ','Communications' ,'ALLReduce orbs' ,  &  !< Lagrange Multipliers communication
       'Diis          ','Other         ' ,'Other         ' ,  &  
       'PSolv_comput  ','Potential     ' ,'3D FFT        ' ,  &  
       'PSolv_commun  ','Communications' ,'ALLtoALL      ' ,  &  
       'PSolvKernel   ','Initialization' ,'Miscellaneous ' ,  &  
       'Exchangecorr  ','Potential     ' ,'Miscellaneous ' ,  &  
       'Forces        ','Finalization  ' ,'Miscellaneous ' ,  &  
       'Tail          ','Finalization  ' ,'Miscellaneous ' ,  &
       'Loewdin_comput','Linear Algebra' ,'              ' ,  &
       'Loewdin_commun','Communications' ,'ALLReduce orbs' ,  &
       'Chol_commun   ','Communications' ,'              ' ,  &
       'Chol_comput   ','Linear Algebra' ,'ALLReduce orbs' ,  &
       'GS/Chol_comput','Linear Algebra' ,'              ' ,  &
       'GS/Chol_commun','Communications' ,'ALLReduce orbs' ,  &
       'Input_comput  ','Initialization' ,'Miscellaneous ' ,  &
       'Input_commun  ','Communications' ,'ALLtoALL+Reduc' ,  &
       'Davidson      ','Finalization  ' ,'Complete SCF  ' ,  &
       'check_IG      ','Initialization' ,'Linear Scaling' ,  &
       'constrc_locreg','Initialization' ,'Miscellaneous ' ,  &
       'wavefunction  ','Initialization' ,'Miscellaneous ' ,  &
       'create_nlpspd ','Initialization' ,'RMA pattern   ' ,  &
       'p2pOrtho_post ','Communications' ,'irecv / irsend' ,  &
       'p2pOrtho_wait ','Communications' ,'mpi_waitany   ' ,  &
       'lovrlp_comm   ','Communications' ,'mpi_allgatherv' ,  &
       'lovrlp_comp   ','Linear Algebra' ,'many ddots    ' ,  &
       'lovrlp_compr  ','Other         ' ,'cut out zeros ' ,  &
       'lovrlp_uncompr','Other         ' ,'insert zeros  ' ,  &
       'extract_orbs  ','Other         ' ,'copy to sendb ' ,  &
       'lovrlp^-1/2   ','Linear Algebra' ,'exact or appr ' ,  &
       'lovrlp^-1/2old','Linear Algebra' ,'exact or appr ' ,  &
       'lovrlp^-1/2com','Linear Algebra' ,'exact or appr ' ,  &
       'lovrlp^-1/2par','Linear Algebra' ,'exact or appr ' ,  &
       'build_lincomb ','Linear Algebra' ,'many daxpy    ' ,  &
       'convolQuartic ','Convolutions  ' ,'No OpenCL     ' ,  &
       'p2pSumrho_wait','Communications' ,'mpi_test/wait ' ,  &
       'sumrho_TMB    ','Other         ' ,'port to GPU?  ' ,  &
       'TMB_kernel    ','Linear Algebra' ,'dgemm         ' ,  &
       'diagonal_seq  ','Linear Algebra' ,'dsygv         ' ,  &
       'diagonal_par  ','Linear Algebra' ,'pdsygvx       ' ,  &
       'lovrlp^-1     ','Linear Algebra' ,'exact or appr ' ,  &
       'lagmat_orthoco','Linear Algebra' ,'dgemm seq/par ' ,  &
       'optimize_DIIS ','Other         ' ,'Other         ' ,  &
       'optimize_SD   ','Other         ' ,'Other         ' ,  &
       'mix_linear    ','Other         ' ,'Other         ' ,  &
       'mix_DIIS      ','Other         ' ,'Other         ' ,  &
       'ig_matric_comm','Communications' ,'mpi p2p       ' ,  &
       'wf_signals    ','Communications' ,'Socket transf.' ,  &
       'energs_signals','Communications' ,'Socket transf.' ,  &
       'rhov_signals  ','Communications' ,'Socket transf.' ,  &
       'init_locregs  ','Initialization' ,'Miscellaneous ' ,  &
       'init_commSumro','Initialization' ,'Miscellaneous ' ,  &
       'init_commPot  ','Initialization' ,'Miscellaneous ' ,  &
       'init_commOrtho','Initialization' ,'Miscellaneous ' ,  &
       'init_inguess  ','Initialization' ,'Miscellaneous ' ,  &
       'init_matrCompr','Initialization' ,'Miscellaneous ' ,  &
       'init_collcomm ','Initialization' ,'Miscellaneous ' ,  &
       'init_collco_sr','Initialization' ,'Miscellaneous ' ,  &
       'init_orbs_lin ','Initialization' ,'Miscellaneous ' ,  &
       'init_repart   ','Initialization' ,'Miscellaneous ' ,  &
       'initMatmulComp','Initialization' ,'Miscellaneous ' ,  &
       'Pot_after_comm','Other         ' ,'global_to_loca' ,  & 
       'Init to Zero  ','Other         ' ,'Memset        ' ,  &
       'calc_kernel   ','Other         ' ,'Miscellaneous ' ,  &
       'commun_kernel ','Communications' ,'mpi_allgatherv' ,  &
       'getlocbasinit ','Other         ' ,'Miscellaneous ' ,  &
       'updatelocreg1 ','Other         ' ,'Miscellaneous ' ,  &
       'linscalinit   ','Other         ' ,'Miscellaneous ' ,  &
       'commbasis4dens','Communications' ,'Miscellaneous ' ,  &
       'eglincomms    ','Communications' ,'Miscellaneous ' ,  &
       'allocommsumrho','Communications' ,'Miscellaneous ' ,  &
       'ovrlptransComp','Other         ' ,'Miscellaneous ' ,  &
       'ovrlptransComm','Communications' ,'mpi_allreduce ' ,  &
       'lincombtrans  ','Other         ' ,'Miscellaneous ' ,  &
       'glsynchham1   ','Other         ' ,'Miscellaneous ' ,  &
       'glsynchham2   ','Other         ' ,'Miscellaneous ' ,  &
       'gauss_proj    ','Other         ' ,'Miscellaneous ' ,  &
       'sumrho_allred ','Communications' ,'mpiallred     ' ,  &
       'deallocprec   ','Other         ' ,'Miscellaneous ' ,  &
       'large2small   ','Other         ' ,'Miscellaneous ' ,  &
       'small2large   ','Other         ' ,'Miscellaneous ' ,  &
       'renormCoefCom1','Linear Algebra' ,'Miscellaneous ' ,  &
       'renormCoefCom2','Linear Algebra' ,'Miscellaneous ' ,  &
       'renormCoefComm','Communications' ,'Miscellaneous ' ,  &
       'waitAllgatKern','Other         ' ,'Miscellaneous ' ,  &
       'UnBlockPot    ','Other         ' ,'Overlap comms ' ,  &
       'UnBlockDen    ','Other         ' ,'Overlap comms ' ,  &
       'global_local  ','Initialization' ,'Unknown       ' ,  &
       'wfd_creation  ','Other         ' ,'Miscellaneous ' ,  & 
       'comm_llr      ','Communications' ,'Miscellaneous ' ,  &
       'AllocationProf','Other         ' ,'Allocate arrs ' ,  &
       'dirmin_lagmat1','Linear Algebra' ,'grad calc     ' ,  &
       'dirmin_lagmat2','Linear Algebra' ,'allgatherv    ' ,  &
       'dirmin_dgesv  ','Linear Algebra' ,'dgesv/pdgesv  ' ,  &
       'dirmin_sddiis ','Linear Algebra' ,'Miscellaneous ' ,  &
       'dirmin_allgat ','Linear Algebra' ,'allgatherv    ' ,  &
       'dirmin_sdfit  ','Linear Algebra' ,'allgatherv etc' ,  &
       'chebyshev_comp','Linear Algebra' ,'matmul/matadd ' ,  &
       'chebyshev_comm','Communications' ,'allreduce     ' ,  &
       'chebyshev_coef','Other         ' ,'Miscellaneous ' ,  &
       'FOE_auxiliary ','Other         ' ,'Miscellaneous ' ,  &
       'FOE_init      ','Other         ' ,'Miscellaneous ' ,  &
       'compress_uncom','Other         ' ,'Miscellaneous ' ,  &
       'norm_trans    ','Other         ' ,'Miscellaneous ' ,  &
       'misc          ','Other         ' ,'Miscellaneous ' ,  &
       'sparse_copy   ','Other         ' ,'Miscellaneous ' ,  &
       'constraineddft','Other         ' ,'Miscellaneous ' ,  &
       'transfer_int  ','Other         ' ,'Miscellaneous ' ,  &
       'Reformatting  ','Initialization' ,'Interpolation ' ,  &
       'restart_wvl   ','Initialization' ,'inguess    rst' ,  &
       'restart_rsp   ','Initialization' ,'inguess    rst' ,  &
       'check_sumrho  ','Initialization' ,'unitary check ' ,  &
       'check_pot     ','Initialization' ,'unitary check ' ,  &
       'ApplyLocPot   ','Convolutions  ' ,'OpenCL ported ' ,  &
       'ApplyLocKin   ','Convolutions  ' ,'OpenCL ported ' ,  &
       'kernel_init   ','Other         ' ,'Fragment calc ' ,  &
       'calc_energy   ','Linear Algebra' ,'allred etc    ' ,  &
       'new_pulay_corr','Other         ' ,'Pulay forces  ' ,  &
       'dev_from_unity','Other         ' ,'Miscellaneous ' ,  &
       'ks_residue    ','Linear Algebra' ,'Miscellaneous ' ,  &
       'weightanalysis','Linear Algebra' ,'Fragment calc ' ,  &
       'tmbrestart    ','Initialization' ,'Miscellaneous ' ,  &
       'readtmbfiles  ','Initialization' ,'Miscellaneous ' ,  &
       'readisffiles  ','Initialization' ,'Miscellaneous ' ,  &
       'purify_kernel ','Linear Algebra' ,'dgemm         ' ,  &
       'potential_dims','Other         ' ,'auxiliary     ' ,  &
       'calc_bounds   ','Other         ' ,'Miscellaneous ' /),(/3,ncat_max/))
  logical :: parallel,init,debugmode
  integer :: ncaton,nproc = 0,nextra,ncat_stopped
  real(kind=8) :: time0,t0
  real(kind=8), dimension(ncat_max+1) :: timesum=0.0d0
  real(kind=8), dimension(ncat_max) :: pctimes=0.0d0 !total times of the partial counters
  character(len=10), dimension(ncat_max) :: pcnames=repeat(' ',10) !names of the partial counters, to be assigned
  character(len=128) :: filename_time
  !integer indicating the file unit
  integer, parameter :: timing_unit=60
  !integer indicating the tabbing of the file unit
  integer, parameter :: tabfile=25

  !categories definition
  type(dictionary), pointer :: dict_timing_categories=>null()  !< categories
  type(dictionary), pointer :: dict_timing_groups=>null()  !< groups

  !entries of the dictionary 
  character(len=*), parameter :: catname='Category'
  character(len=*), parameter :: grpname='Group'
  character(len=*), parameter :: catinfo='Info'

  !>number of groups 
  integer :: timing_ncls=0!ncls_max
  !>number of categories
  integer :: timing_ncat=0!ncat_max
  !>number of partial counters
  integer :: timing_nctr=0

  !unspecified timings
  integer, save :: TIMING_CAT_UNSPEC
  !Error codes
  integer, save :: TIMING_INVALID


  contains

    !> check if the module has been initialized
    subroutine check_initialization()
      implicit none
      if (.not. associated(dict_timing_categories)) then !call f_err_initialize()
         write(0,*)'Timing library not initialized, f_lib_initialized should be called'
         call f_err_severe()
      end if
    end subroutine check_initialization


    !>for the moment the timing callback is a severe error.
    !! we should decide what to do to override this
    subroutine f_timing_callback()
      use yaml_output, only: yaml_warning,yaml_map
      use dictionaries, only: f_err_severe
      implicit none
      call yaml_warning('An error occured in timing module')
      call yaml_map('Dictionary of category groups',dict_timing_groups)
      call yaml_map('Dictionary of active categories',dict_timing_categories)
      call f_err_severe()
    end subroutine f_timing_callback

    !> initialize error codes of timing module
    subroutine timing_errors()
      use dictionaries, only: f_err_define
      implicit none
       call f_err_define(err_name='TIMING_INVALID',err_msg='Error in timing routines',&
            err_id=TIMING_INVALID,&
            err_action='Control the order of the timing routines called',&
            callback=f_timing_callback)
    end subroutine timing_errors

    !> define a class, which is a group of categories
    subroutine f_timing_category_group(grp_name,grp_info)
      implicit none
      character(len=*), intent(in) :: grp_name !< name of the class
      character(len=*), intent(in) :: grp_info !< description of it
      !local variables

      call check_initialization()

      if (grp_name .in. dict_timing_groups) then
         call f_err_throw('The timing category group '//grp_name//' has already been defined',&
              err_id=TIMING_INVALID)
      else
         timing_ncls=timing_ncls+1
      end if
      !in case of dry run override the commentary nonetheless  
      call set(dict_timing_groups//grp_name,grp_info)

    end subroutine f_timing_category_group

    !> define a new timing category with its description
    subroutine f_timing_category(cat_name,grp_name,cat_info,cat_id)
      implicit none
      character(len=*), intent(in) :: cat_name !< name of the category
      character(len=*), intent(in) :: grp_name !<class to which category belongs (see f_timing_class)
      character(len=*), intent(in) :: cat_info !< description of it
      integer, intent(out) :: cat_id !< id of the defined class, to be used for reference
      !local variables
      type(dictionary), pointer :: dict_cat
     
      call check_initialization()

      if (.not. (grp_name .in. dict_timing_groups)) then
         call f_err_throw('The timing category group '//grp_name//' has not been defined',&
              err_id=TIMING_INVALID)
         return
      end if

      !then proceed to the definition of the category
      cat_id=dict_len(dict_timing_categories)

      call dict_init(dict_cat)
      call set(dict_cat//catname,cat_name)
      call set(dict_cat//grpname,grp_name)
      call set(dict_cat//catinfo,cat_info)

      call add(dict_timing_categories,dict_cat)
      timing_ncat=timing_ncat+1
    end subroutine f_timing_category

    !initialize the timing by putting to zero all the chronometers
    subroutine f_timing_initialize()
      implicit none
      !initialize errors
      call timing_errors()
      !create the general category for unspecified timings
      call dict_init(dict_timing_groups)
      call dict_init(dict_timing_categories)

      !define the main groups and categories
      call f_timing_category_group('NULL','Nullified group to contain unspecifed category')
      call f_timing_category('UNSPECIFIED','NULL',&
           'Unspecified category collecting garbage timings',TIMING_CAT_UNSPEC)

      !to be moved somewhere else
      call timing_initialize_categories()
    end subroutine f_timing_initialize

    !finalizee the timing by putting to zero all the chronometers
    subroutine f_timing_finalize()
      implicit none
      !create the general category for unspecified timings
      call dict_free(dict_timing_categories)
      call dict_free(dict_timing_groups)
      !put to zero the number of classes and the number of categories
      timing_ncls=0
      timing_ncat=0
    end subroutine f_timing_finalize

    !> this routine should go in the bigdft_init routine as the categories are specific to BigDFT
    subroutine timing_initialize_categories()
      implicit none
      !local variables
      integer :: icls,icat
      integer, dimension(ncat_max) :: cat_ids
      
      !initialize groups
      do icls=1,ncls_max
         call f_timing_category_group(trim(clss(icls)),'Empty description for the moment')
      end do
      !initialize categories
      do icat=1,ncat_max
         call f_timing_category(trim(cats(1,icat)),trim(cats(2,icat)),trim(cats(3,icat)),&
              cat_ids(icat))
      end do
      !then fill the cat ids into parameters
    end subroutine timing_initialize_categories

    !re-initialize the timing by putting to zero all the chronometers (old action IN)
    subroutine f_timing_reset(filename,master,verbose_mode)
      use yaml_output, only: yaml_new_document
      implicit none
      logical, intent(in) :: master !<true if the task is the one responsible for file writing
      character(len=*), intent(in) :: filename !<name of the file where summary have to be written
      !>Toggle verbose mode in the module. In case of parallel execution all the processors will
      !! print out their information on the counters.
      logical, intent(in), optional :: verbose_mode 
      !local variables
      integer :: icat,ictr,i,iunit_def
      integer(kind=8) :: itns

      !global timer
      itns=f_time()

      call check_initialization()

      !check if some categories have been initialized
      if (f_err_raise(timing_ncat==0,'No timing categories have been initialized, no counters to reset .'//&
           'Use f_timing_category(_group) routine(s)',err_id=TIMING_INVALID)) return

      time0=real(itns,kind=8)*1.d-9
      !reset partial counters and categories
      do i=1,timing_ncat
         timesum(i)=0.d0
      end do
      do ictr=1,timing_nctr
         pctimes(ictr)=0.d0
      enddo
      !store filename where report have to be written
      filename_time(1:len(filename_time))=filename
      !store debug mode
      if (present(verbose_mode)) then
         debugmode=verbose_mode
      else
         debugmode=.false.
      end if

      !no category has been initialized so far
      init=.false.
      ncat_stopped=0 !no stopped category
      timing_nctr=0 !no partial counters activated
      !initialize the document
      if (master) then
         call timing_open_stream(iunit_def)
         call yaml_new_document() !in principle is active only when the document is released
         call timing_close_stream(iunit_def)
      end if
    end subroutine f_timing_reset


    !>perform a checkpoint of the chronometer with a partial counter
    !! the last active category is halted and a summary of the timing 
    !! is printed out
    subroutine f_timing_checkpoint(ctr_name,mpi_comm)
      !use yaml_output, only: yaml_map
      implicit none
      !> name of the partial counter for checkpoint identification
      character(len=*), intent(in) :: ctr_name 
      !> handle of the mpi_communicator associated with the checkpoint 
      integer, intent(in), optional :: mpi_comm 
      !local variables
      integer :: i
      integer(kind=8) :: itns 

      !global timer
      itns=f_time()

      call check_initialization()

      !stop partial counters and restart from the beginning
      if (init) then
         print *, 'ERROR: TIMING IS INITIALIZED BEFORE PARTIAL RESULTS'
         stop 
      endif
      timing_nctr=timing_nctr+1
      if (f_err_raise(timing_nctr > ncat_max,&
           'Max No. of partial counters reached',err_id=TIMING_INVALID)) return
      !name of the category
      pcnames(timing_nctr)=trim(ctr_name)
      !total time elapsed in the category
      timesum(ncat_max+1)=real(itns,kind=8)*1.d-9-time0
      pctimes(timing_nctr)=timesum(ncat_max+1)
      !here iproc is the communicator
      call sum_results(ncat_max,ncls_max,nextra,&
           debugmode,mpi_comm,pcnames(timing_nctr),timesum,clss,cats)
      !call sum_results(iproc_true,iproc,pcnames(ncounters))
      !reset all timings
      time0=real(itns,kind=8)*1.d-9
      do i=1,timing_ncat
         timesum(i)=0.d0
      enddo
      
    end subroutine f_timing_checkpoint

    !>stop the timing and dump information of the partial counters
    subroutine f_timing_stop(mpi_comm)
      implicit none
      integer, intent(in) :: mpi_comm !< communicator for the results
      !local variables
      integer(kind=8) :: itns

      !global timer
      itns=f_time()

      if (init) then
         print *, 'TIMING IS INITIALIZED BEFORE RESULTS'
         stop 
      endif

      if (timing_nctr == 0) then !no partial counters selected
         timesum(ncat_max+1)=real(itns,kind=8)*1.d-9-time0
         !here iproc is the communicator
         call sum_results(ncat_max,ncls_max,&
              nextra,&
              debugmode,mpi_comm,'ALL',timesum,clss,cats)
      else !consider only the results of the partial counters
         call sum_counters(pctimes,pcnames,timing_nctr,mpi_comm,&
              debugmode)
      end if

    end subroutine f_timing_stop

    !> The same timing routine but with system_clock (in case of a supported specs)
    subroutine f_timing(category,action)
      use dictionaries, only: f_err_raise,f_err_throw
      use dynamic_memory, only: f_time
      use yaml_output, only: yaml_toa
      implicit none
      !Variables
      character(len=*), intent(in) :: category
      character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
      !Local variables
      logical :: catfound
      integer :: i,ierr,cat_id
      integer :: nthreads,jproc,namelen
      integer(kind=8) :: itns
      real(kind=8) :: t1
      character(len=128) :: cattmp

      !first of all, read the time
      itns=f_time()

      !find category in the old scheme
      call find_category(category,cat_id)

      cattmp=dict_timing_categories//cat_id//catname
      if (f_err_raise(trim(cattmp)/=trim(category),'Error in category '//&
           trim(yaml_toa(cat_id))//' (name='//trim(category)//' ), found '//&
           trim(cattmp)//' instead',err_id=TIMING_INVALID)) return

      select case(action)
      case('ON')
         !if some other category was initalized before, return (no action)
         if (init) return
         t0=real(itns,kind=8)*1.d-9
         init=.true.
         ncaton=cat_id !category which has been activated
      case('OF')
         if (cat_id==ncaton) then !switching off the good category
            if (f_err_raise(.not. init,'Timing category '//&
                 trim(cattmp)//' not initialized',err_id=TIMING_INVALID)) &
                 return
            t1=real(itns,kind=8)*1.d-9
            timesum(cat_id)=timesum(cat_id)+t1-t0
            init=.false.
         else !no action except for misuse of interrupts
            if (ncat_stopped /=0) &
                 call f_err_throw('INTERRUPTS SHOULD NOT BE HALTED BY OF',&
                 err_id=TIMING_INVALID)
            return
            !interrupt the active category and replace it by the proposed one
         end if
      case('IR') !interrupt category
         if (ncat_stopped /=0) then
            call f_err_throw('Category No. '//trim(yaml_toa(ncat_stopped))//&
                 ' already exclusively initialized',err_id=TIMING_INVALID)
            return
         end if
         !time
         t1=real(itns,kind=8)*1.d-9
         if (init) then !there is already something active
            !stop the active counter
            timesum(ncaton)=timesum(ncaton)+t1-t0
            ncat_stopped=ncaton
         else
            init=.true.
            ncat_stopped=-1
         end if
         ncaton=cat_id
         t0=t1
      case('RS') !resume the interrupted category
         if (f_err_raise(ncat_stopped ==0,'It appears no category has to be resumed',&
              err_id=TIMING_INVALID)) return
         if (f_err_raise(cat_id /= ncaton,'The category id '//trim(yaml_toa(cat_id))//' is not active',&
              err_id=TIMING_INVALID)) return
         !time
         t1=real(itns,kind=8)*1.d-9
         timesum(cat_id)=timesum(cat_id)+t1-t0
         if (ncat_stopped == -1) then
            init =.false. !restore normal counter
         else
            ncaton=ncat_stopped
            t0=t1
         end if
         ncat_stopped=0
      case default
         call f_err_throw('TIMING ACTION UNDEFINED',err_id=TIMING_INVALID)
         !print *,action,cat_id,ncaton,trim(category)
         !stop 
      end select

    contains
      
      subroutine find_category(category,ii)
        implicit none
        character(len=*), intent(in) :: category
        integer, intent(out) :: ii !< id of the found category
        !local variables
        integer :: i
        !controls if the category exists
        catfound=.false.
        ii=0
        do i=1,timing_ncat
           if (trim(category) == trim(cats(1,i))) then
              ii=i
              exit
           endif
        enddo
        if (ii==0) then
           call f_err_throw('Timing routine error,'//&
             ' the category '//trim(category)//' requested for action '//&
             trim(action)//' has not been found',err_id=TIMING_INVALID)
        end if
      end subroutine find_category

    END SUBROUTINE f_timing

    !>opens the file of the timing unit
    subroutine timing_open_stream(iunit_def)
      use yaml_output, only: yaml_get_default_stream,yaml_set_stream
      implicit none
      integer, intent(out) :: iunit_def !< previous default unit
      !first get the default stream
      call yaml_get_default_stream(iunit_def)
      if (iunit_def /= timing_unit) then
         call yaml_set_stream(unit=timing_unit,filename=trim(filename_time),&
              record_length=120,tabbing=tabfile)
      end if

    end subroutine timing_open_stream


    !> close the stream and restore old default unit
    subroutine timing_close_stream(iunit_def)
      use yaml_output, only: yaml_set_default_stream,yaml_close_stream
      implicit none
      integer, intent(in) :: iunit_def !< previous default unit
      !local variables
      integer :: ierr

      call yaml_set_default_stream(iunit_def,ierr)
      !close the previous one
      if (iunit_def /= timing_unit) call yaml_close_stream(unit=timing_unit)
    end subroutine timing_close_stream

    !> dump the line of the timings for time.yaml form
    subroutine timing_dump_line(name,tabbing,pc,secs,unit,loads)
      use yaml_output
      implicit none
      integer, intent(in) :: tabbing !<vlue of the tabbing for pretty printing
      double precision, intent(in) :: pc !< percent of the time for line id
      double precision, intent(in) :: secs !< seconds spent for line id
      character(len=*), intent(in) :: name !< id of the line printed
      integer, intent(in), optional :: unit !< @copydoc yaml_output::doc::unit
      !> extra info containing the load for each task for the line id, 
      !! calculated with respect to the average value given by secs
      double precision, dimension(0:), intent(in), optional :: loads 
      !local variables
      character(len=*), parameter :: fmt_pc='(f5.1)'
      character(len=*), parameter :: fmt_secs='(1pg9.2)'
      character(len=*), parameter :: fmt_extra='(f5.2)'
      integer :: iextra,nextra,unt

      unt=0
      if (present(unit)) unt=unit
      !determine the presence of extra information
      nextra=0
      if (present(loads)) nextra=size(loads)

      call yaml_open_sequence(name,flow=.true.,tabbing=tabbing,unit=unt)
      call yaml_sequence(yaml_toa(pc,fmt=fmt_pc),unit=unt)
      call yaml_sequence(yaml_toa(secs,fmt=fmt_secs),unit=unt)
      do iextra=0,nextra-1
         call yaml_sequence(yaml_toa(loads(iextra),fmt=fmt_extra),unit=unt)
      end do
      call yaml_close_sequence(unit=unt)

    end subroutine timing_dump_line

    !>put the average value of timeall in the timesum array
    !then rewrite each element with the deviation from it (in debug mode)
    !in normal mode write only the max and min deviations (only in parallel)
    subroutine timing_data_synthesis(nproc,ncats,timeall,timesum)
      implicit none
      integer, intent(in) :: nproc,ncats
      real(kind=8), dimension(ncats,0:nproc-1), intent(inout) :: timeall
      real(kind=8), dimension(ncats), intent(out) :: timesum
      !local variables
      integer :: icat,jproc
      real(kind=8) :: tmin,tmax

      do icat=1,ncats
         timesum(icat)=0.d0
         do jproc=0,nproc-1
            timesum(icat)=timesum(icat)+timeall(icat,jproc)
         end do
         timesum(icat)=timesum(icat)/real(nproc,kind=8)
         if (timesum(icat)>0.d0) then
            if (debugmode) then
               do jproc=0,nproc-1
                  timeall(icat,jproc)=timeall(icat,jproc)/timesum(icat)
               end do
            else if (nproc >1) then
               tmax=0.0d0
               tmin=1.0d300
               do jproc=0,nproc-1
                  tmax=max(timeall(icat,jproc),tmax)
                  tmin=min(timeall(icat,jproc),tmin)
               end do
               timeall(icat,0)=tmax/timesum(icat)
               timeall(icat,1)=tmin/timesum(icat)
            end if
         end if
      end do
    end subroutine timing_data_synthesis

  end module time_profiling

  subroutine sum_counters(pctimes,pcnames,ncounters,mpi_comm,debugmode)
    use yaml_output
    use dynamic_memory
    use time_profiling, only: timing_unit,timing_dump_line,timing_data_synthesis,&
         timing_open_stream,timing_close_stream

  implicit none
  include 'mpif.h'
  logical, intent(in) :: debugmode
  integer, intent(in) :: mpi_comm,ncounters
  real(kind=8), dimension(ncounters), intent(in) :: pctimes
  character(len=10), dimension(ncounters), intent(in) :: pcnames
  !local variables
  integer, parameter :: tabfile=25
  logical :: parallel
  integer :: i,ierr,iproc,jproc,icat,nthreads,namelen,iunit_def,nproc
  real(kind=8) :: pc
  
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  double precision, dimension(:,:), allocatable :: timecnt 
  character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
  !$ integer :: omp_get_max_threads

  ! Not initialised case.
  if (mpi_comm==MPI_COMM_NULL) return
  call f_routine(id='sum_counters')
  
  call MPI_COMM_SIZE(mpi_comm,nproc,ierr)
  parallel=nproc>1

  nodename=f_malloc_str(MPI_MAX_PROCESSOR_NAME,0.to.nproc-1,id='nodename')
  timecnt=f_malloc((/1.to.ncounters,0.to.nproc/),id='timecnt')

  if (parallel) then 
     call MPI_COMM_RANK(mpi_comm,iproc,ierr)
     call MPI_GATHER(pctimes,ncounters,MPI_DOUBLE_PRECISION,&
          timecnt,ncounters,MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
     if (debugmode) then
        !initalise nodenames
        do jproc=0,nproc-1
           nodename(jproc)=repeat(' ',MPI_MAX_PROCESSOR_NAME)
        end do

        call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)

        !gather the result between all the process
        call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
             nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
             mpi_comm,ierr)
     end if

  else
     do i=1,ncounters
        timecnt(i,0)=pctimes(i)
     end do
     iproc=0
  endif

  if (iproc == 0) then
     !synthesis of the counters
     call timing_data_synthesis(nproc,ncounters,timecnt,timecnt(1,nproc))

     call timing_open_stream(iunit_def)
!!$     !first get the default stream
!!$     call yaml_get_default_stream(iunit_def)
!!$     if (iunit_def /= timing_unit) then
!!$        call yaml_set_stream(unit=timing_unit,filename=trim(filename_time),record_length=120,&
!!$             istat=ierr,tabbing=tabfile)
!!$     end if
     call yaml_open_map('SUMMARY',advance='no')
     call yaml_comment('     % ,  Time (s)',tabbing=tabfile)
     
     !sum all the information by counters
     do i=1,ncounters
        pc=100.d0*timecnt(i,nproc)/sum(timecnt(1:ncounters,nproc))
        call timing_dump_line(trim(pcnames(i)),tabfile,pc,timecnt(i,nproc))
     end do
     call timing_dump_line('Total',tabfile,100.d0,sum(timecnt(1:ncounters,nproc)))
     call yaml_close_map() !summary
     call yaml_open_map('CPU parallelism')
     call yaml_map('MPI_tasks',nproc)
     nthreads = 0
     !$  nthreads=omp_get_max_threads()
     if (nthreads /= 0) call yaml_map('OMP threads',nthreads)
     call yaml_close_map()
     if (debugmode) then
        call yaml_open_sequence('Hostnames')
        do jproc=0,nproc-1
           call yaml_sequence(trim(nodename(jproc)))
        end do
        call yaml_close_sequence()
     end if
     call yaml_map('Report timestamp',trim(yaml_date_and_time_toa()))
     !restore the default stream
     call timing_close_stream(iunit_def)
!!$     call yaml_set_default_stream(iunit_def,ierr)
!!$     !close the previous one
!!$     if (iunit_def /= timing_unit) call yaml_close_stream(unit=timing_unit)
  end if
  call f_free(timecnt)
  call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)
  call f_release_routine()
end subroutine sum_counters


subroutine sum_results(ncat,ncls,nextra,&
     debugmode,mpi_comm,message,timesum,clss,cats)
  use dynamic_memory
  use yaml_output
  use time_profiling, only: timing_unit,timing_dump_line,timing_data_synthesis,&
         timing_open_stream,timing_close_stream
  implicit none
  include 'mpif.h'
  logical, intent(in) :: debugmode
  integer, intent(in) :: mpi_comm,ncat,ncls,nextra
  character(len=*), intent(in) :: message
 character(len=14), dimension(ncls), intent(in) :: clss
  character(len=14), dimension(3,ncat), intent(in) :: cats
!  logical, intent(inout) :: newfile
  real(kind=8), dimension(ncat+1), intent(inout) :: timesum
   !local variables
  integer, parameter :: tabfile=25
  logical :: parallel
  integer :: i,ierr,j,icls,icat,jproc,iextra,iproc,iunit_def,nproc
  real(kind=8) :: total_pc,pc
  integer, dimension(ncat) :: isort
  real(kind=8), dimension(:,:), allocatable :: timecls
  real(kind=8), dimension(:,:), allocatable :: timeall

  ! Not initialised case.
  if (mpi_comm==MPI_COMM_NULL) return
  call f_routine(id='sum_results')

  call MPI_COMM_SIZE(mpi_comm,nproc,ierr)
  parallel=nproc>1
  !allocate total timings
  timecls=f_malloc((/1.to.ncls,0.to.nproc/),id='timecls')
  timeall=f_malloc((/1.to.ncat+1,0.to.nproc-1/),id='timeall')

  if (parallel) then
     call MPI_COMM_RANK(mpi_comm,iproc,ierr)
     call MPI_GATHER(timesum,ncat+1,MPI_DOUBLE_PRECISION,&
          timeall,ncat+1,MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
  else
     do i=1,ncat+1
        timeall(i,0)=timesum(i)
     end do
     iproc=0
  endif
  if (iproc == 0) then
     !regroup the data for each category in any processor
     !this has to be done via the dictionary of the categories
     do icls=1,ncls
        timecls(icls,0:nproc)=0.d0 
        do icat=1,ncat
           if(trim(cats(2,icat))==clss(icls)) then
              do jproc=0,nproc-1
                 timecls(icls,jproc)=timecls(icls,jproc)+timeall(icat,jproc)
              end do
           end if
        end do
     end do

     !synthesis of the categories
     call timing_data_synthesis(nproc,ncat+1,timeall,timesum)
     !synthesis of the classes
     call timing_data_synthesis(nproc,ncls,timecls,timecls(1,nproc))

     !calculate the summary of the category
     call sort_positions(ncat,timesum,isort)

!!!! !use yaml to write time.yaml
     call timing_open_stream(iunit_def)
!!$     !first get the default stream
!!$     call yaml_get_default_stream(iunit_def)
!!$     if (iunit_def /= timing_unit) then
!!$        call yaml_set_stream(unit=timing_unit,filename=trim(filename_time),record_length=120,&
!!$             istat=ierr,tabbing=tabfile)
!!$     end if
     !start the writing of the file

!!$     if (newfile) then
!!$        call yaml_new_document() !in principle is active only when the document is released
!!$        newfile=.false.
!!$     end if

         call yaml_open_map(trim(message),advance='no')
         if (.not. parallel) then
            call yaml_comment('     % ,  Time (s)',tabbing=tabfile)
         else if (debugmode) then
            call yaml_comment('     % ,  Time (s), Load per MPI proc (relative) ',tabbing=tabfile)
         else
            call yaml_comment('     % ,  Time (s), Max, Min Load (relative) ',tabbing=tabfile)
         end if
         call yaml_open_map('Classes')
         total_pc=0.d0
         do icls=1,ncls
            pc=0.0d0
            if (timesum(ncat+1)/=0.d0) &
                 pc=100.d0*timecls(icls,nproc)/timesum(ncat+1)
            total_pc=total_pc+pc
            !only nonzero classes are printed out
            if (timecls(icls,nproc) /= 0.d0) then
               call timing_dump_line(trim(clss(icls)),tabfile,pc,timecls(icls,nproc),&
                    loads=timecls(icls,0:nextra-1))
            end if
         end do
         call timing_dump_line('Total',tabfile,total_pc,timesum(ncat+1),&
                    loads=timeall(ncat+1,0:nextra-1))
         call yaml_open_map('Categories',advance='no')
         call yaml_comment('In order of time consumption')
         do j=1,ncat
            i=isort(j)
            pc=0.d0
            !only nonzero categories are printed out
            if (timesum(i) /= 0.d0) then
               if (timesum(ncat+1)/=0.d0)&
                    pc=100.d0*timesum(i)/timesum(ncat+1)
               call timing_dump_line(trim(cats(1,i)),tabfile,pc,timesum(i),&
                    loads=timeall(i,0:nextra-1))
               call yaml_map('Class',trim(cats(2,i)))
               call yaml_map('Info',trim(cats(3,i)))
            end if
         enddo
         call yaml_close_map() !categories
         call yaml_close_map() !counter
         !restore the default stream
         call timing_close_stream(iunit_def)
!!$         call yaml_set_default_stream(iunit_def,ierr)
!!$         !close the previous one
!!$         if (iunit_def /= timing_unit) call yaml_close_stream(unit=timing_unit)
      endif

      call f_free(timecls)
      call f_free(timeall)
      call f_release_routine()

END SUBROUTINE sum_results


!> The same timing routine but with system_clock (in case of a supported specs)
subroutine timing(iproc,category,action)
  use time_profiling, ncat => ncat_max, ncls => ncls_max

  implicit none

  include 'mpif.h'
  !Variables
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: category
  character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
  !Local variables
  logical :: catfound
  integer :: i,ierr,ii,iproc_true
  integer :: nthreads,jproc,namelen
  integer(kind=8) :: itns
  !cputime routine gives a real
  !real :: total,total0,time,time0
  real(kind=8) :: pc,t1
!!$  real(kind=8), dimension(ncounters,0:nproc) :: timecnt !< useful only at the very end
!!$  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
!!$  character(len=MPI_MAX_PROCESSOR_NAME), dimension(0:nproc-1) :: nodename

!$ integer :: omp_get_max_threads

  !modification of the timing to see if it works
  select case(action)
  case('IN') ! INIT
     !to be changed IMMEDIATELY
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc_true,ierr)
     call f_timing_reset(filename=category,master=iproc_true==0,&
          verbose_mode=(abs(iproc)==2 .or. iproc<-1))
  case('PR')
     call f_timing_checkpoint(ctr_name=category,mpi_comm=iproc)
  case('RE') ! RESULT
     call f_timing_stop(mpi_comm=iproc)
  case default
     call f_timing(category,action)   
  end select
  
!!$
!!$  !first of all, read the time
!!$  !call system_clock(itime,count_rate,count_max)
!!$  call nanosec(itns)
!!$
!!$  ! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
!!$  if (action.eq.'IN') then  ! INIT
!!$     !!no need of using system clock for the total time (presumably more than a millisecond)
!!$     !call cpu_time(total0)
!!$     filename_time=repeat(' ',128)
!!$     time0=real(itns,kind=8)*1.d-9
!!$     do i=1,ncat
!!$        timesum(i)=0.d0
!!$        pctimes(i)=0.d0
!!$     enddo
!!$     !in this case iproc stands for nproc
!!$     parallel=abs(iproc) > 1!trim(category).eq.'parallel'
!!$     nproc=abs(iproc)
!!$     filename_time=trim(category)
!!$     newfile=.true.
!!$     init=.false.
!!$     debugmode=(nproc == 2) .or. iproc < -1
!!$     if (nproc >=2) then
!!$        nextra=nproc
!!$        if (.not. debugmode) nextra=2
!!$        write(strextra,'(i5)')nextra
!!$        formatstring='1x,f5.1,a,1x,1pe9.2,a,'//trim(strextra)//'(1x,0pf5.2,a)'
!!$     else
!!$        nextra=0
!!$        formatstring='1x,f5.1,a,1x,1pe9.2,a'
!!$     end if
!!$     ncat_stopped=0 !no stopped category
!!$     ncounters=0
!!$
!!$  else if (action.eq.'PR') then !stop partial counters and restart from the beginning
!!$     if (init) then
!!$        print *, 'ERROR: TIMING IS INITIALIZED BEFORE PARTIAL RESULTS'
!!$        stop 
!!$     endif
!!$     ncounters=ncounters+1
!!$     if (ncounters > ncat) then
!!$        print *, 'It is not allowed to have more partial counters that categories; ncat=',ncat
!!$        stop
!!$     end if
!!$     !name of the category
!!$     pcnames(ncounters)=trim(category)
!!$     !total time elapsed in the category
!!$     timesum(ncat+1)=real(itns,kind=8)*1.d-9-time0
!!$     pctimes(ncounters)=timesum(ncat+1)
!!$     !here iproc is the communicator
!!$     call sum_results(newfile,ncat,ncls,nextra,filename_time,&
!!$          debugmode,iproc,pcnames(ncounters),timesum,clss,cats)
!!$     !call sum_results(iproc_true,iproc,pcnames(ncounters))
!!$     !reset all timings
!!$     time0=real(itns,kind=8)*1.d-9
!!$     do i=1,ncat
!!$        timesum(i)=0.d0
!!$     enddo
!!$
!!$  else if (action.eq.'RE') then ! RESULT
!!$     if (init) then
!!$        print *, 'TIMING IS INITIALIZED BEFORE RESULTS'
!!$        stop 
!!$     endif
!!$     if (ncounters == 0) then !no partial counters selected
!!$        timesum(ncat+1)=real(itns,kind=8)*1.d-9-time0
!!$        !here iproc is the communicator
!!$        call sum_results(newfile,ncat,ncls,nextra,filename_time,&
!!$             debugmode,iproc,'ALL',timesum,clss,cats)
!!$        !call sum_results(iproc_true,iproc,'ALL')
!!$     else !consider only the results of the partial counters
!!$        call sum_counters(pctimes,pcnames,ncounters,iproc,&
!!$             debugmode,filename_time)
!!$     end if
!!$  else
!!$     !controls if the category exists
!!$     catfound=.false.
!!$     do i=1,ncat
!!$        if (trim(category) == trim(cats(1,i))) then
!!$           ii=i
!!$           catfound=.true.
!!$           exit
!!$        endif
!!$     enddo
!!$     if (.not. catfound) then
!!$        print *, 'ACTION  ',action
!!$        write(*,*) 'category, action',category, action
!!$        call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$        stop 'TIMING CATEGORY NOT DEFINED'
!!$     end if
!!$
!!$     if (action == 'ON') then  ! ON
!!$        !some other category was initalized before, overriding
!!$!if (iproc==0) print*,'timing on: ',trim(category)
!!$        if (init) return
!!$        t0=real(itns,kind=8)*1.d-9
!!$        init=.true.
!!$        ncaton=ii !category which has been activated
!!$     else if (action == 'OF' .and. ii==ncaton) then  ! OFF
!!$        if (.not. init) then
!!$           print *, cats(1,ii), 'not initialized'
!!$           stop 
!!$        endif
!!$!if (iproc==0) print*,'timing OFF: ',trim(category)
!!$        t1=real(itns,kind=8)*1.d-9
!!$        timesum(ii)=timesum(ii)+t1-t0
!!$        init=.false.
!!$     else if (action == 'OF' .and. ii/=ncaton) then
!!$        if (ncat_stopped /=0) stop 'INTERRUPTS SHOULD NOT BE HALTED BY OF'
!!$!if (iproc==0) print*,'timing2 OFF: ',trim(category)
!!$        !some other category was initalized before, taking that one
!!$        return
!!$    !interrupt the active category and replace it by the proposed one
!!$     else if (action == 'IR') then
!!$        if (ncat_stopped /=0) then
!!$           print *, cats(1,ncat_stopped), 'already exclusively initialized'
!!$           stop
!!$        end if
!!$        !time
!!$        t1=real(itns,kind=8)*1.d-9
!!$        if (init) then !there is already something active
!!$           !stop the active counter
!!$           timesum(ncaton)=timesum(ncaton)+t1-t0
!!$           ncat_stopped=ncaton
!!$        else
!!$           init=.true.
!!$           ncat_stopped=-1
!!$        end if
!!$        ncaton=ii
!!$        t0=t1
!!$
!!$     else if (action == 'RS') then !resume the interrupted category
!!$        if (ncat_stopped ==0) then
!!$           stop 'NOTHING TO RESUME'
!!$        end if
!!$        if (ii /= ncaton) stop 'WRONG RESUMED CATEGORY'
!!$        !time
!!$        t1=real(itns,kind=8)*1.d-9
!!$        timesum(ii)=timesum(ii)+t1-t0
!!$        if (ncat_stopped == -1) then
!!$           init =.false. !restore normal counter
!!$        else
!!$           ncaton=ncat_stopped
!!$           t0=t1
!!$        end if       
!!$        ncat_stopped=0
!!$     else
!!$        print *,action,ii,ncaton,trim(category)
!!$        stop 'TIMING ACTION UNDEFINED'
!!$     endif
!!$
!!$  endif

END SUBROUTINE timing


subroutine sort_positions(n,a,ipiv)
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: a
  integer, dimension(n), intent(out) :: ipiv
  !local variables
  integer :: i,j,jmax,imax
  real(kind=8) :: locmax

  !neutral permutation
  do i=1,n
     ipiv(i)=i
  end do
  !find the order for all the arrays
  do j=1,n
  !search the maximum
     locmax=-1.d300
     do i=j,n
        if (locmax < a(ipiv(i))) then
           locmax=a(ipiv(i))
           jmax=ipiv(i)
           imax=i
        end if
     end do
     !swap the position with j
     ipiv(imax)=ipiv(j) !throw in the present element
     ipiv(j)=jmax       !take out the present maximum
  end do
  !do i=1,n
  !   print *,'a',i,a(i),ipiv(i),a(ipiv(i))
  !end do
  !stop
end subroutine sort_positions

