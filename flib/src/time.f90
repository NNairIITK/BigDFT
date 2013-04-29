!> @file
!!  Define routines for timing
!! @author
!!    Copyright (C) 2010-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>    Contains variables used a timing for BigDFT
module timeData
!  use module_defs, only: mpi_environment, bigdft_mpi
  implicit none
  integer, parameter :: ncat=121,ncls=7   ! define timimg categories and classes
  character(len=14), dimension(ncls), parameter :: clss = (/ &
       'Communications'    ,  &
       'Convolutions  '    ,  &
       'Linear Algebra'    ,  &
       'Other         '    ,  &
       'Potential     '    ,  &
       'Initialization'    ,  &
       'Finalization  '    /)
  character(len=14), dimension(3,ncat), parameter :: cats = reshape((/ &
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
       'dirmin_lagmat1','Linear Algebra' ,'allgatherv etc' ,  &
       'dirmin_lagmat2','Linear Algebra' ,'allreduce etc ' ,  &
       'dirmin_dgesv  ','Linear Algebra' ,'dgesv/pdgesv  ' ,  &
       'dirmin_sddiis ','Linear Algebra' ,'allreduce etc ' ,  &
       'chebyshev_comp','Linear Algebra' ,'matmul/matadd ' ,  &
       'chebyshev_comm','Communications' ,'allreduce     ' ,  &
       'chebyshev_coef','Other         ' ,'Miscellaneous ' ,  &
       'FOE_auxiliary ','Other         ' ,'Miscellaneous ' ,  &
       'FOE_init      ','Other         ' ,'Miscellaneous ' ,  &
       'compress_uncom','Other         ' ,'Miscellaneous ' ,  &
       'norm_trans    ','Other         ' ,'Miscellaneous ' ,  &
       'misc          ','Other         ' ,'Miscellaneous ' ,  &
       'sparse_copy   ','Other         ' ,'Miscellaneous ' ,  &
       'calc_bounds   ','Other         ' ,'Miscellaneous ' /),(/3,ncat/))
  logical :: parallel,init,newfile,debugmode
  integer :: ncounters, ncaton,nproc = 0,nextra,ncat_stopped
  real(kind=8) :: time0,t0
  real(kind=8), dimension(ncat+1) :: timesum
  real(kind=8), dimension(ncat) :: pctimes !total times of the partial counters
  character(len=10), dimension(ncat) :: pcnames !names of the partial counters, to be assigned
  character(len=50) :: formatstring,strextra
  character(len=128) :: filename_time

  contains

    subroutine sum_results(iproc,mpi_comm,message)
      implicit none
      include 'mpif.h'
      character(len=*), intent(in) :: message
      integer, intent(in) :: iproc,mpi_comm
      !local variables
      integer :: i,ierr,j,icls,icat,jproc,iextra

      real(kind=8) :: total_pc,pc
      integer, dimension(ncat) :: isort
      real(kind=8), dimension(ncls,0:nproc) :: timecls
      real(kind=8), dimension(ncat+1,0:nproc-1) :: timeall

      ! Not initialised case.
      if (nproc == 0) return

      if (parallel) then 
         call MPI_GATHER(timesum,ncat+1,MPI_DOUBLE_PRECISION,&
              timeall,ncat+1,MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
      else
         do i=1,ncat+1
            timeall(i,0)=timesum(i)
         end do
      endif
      if (iproc == 0) then
        

         !regroup the data for each category in any processor
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
         call data_synthesis(parallel,debugmode,nproc,ncat+1,timeall,timesum)
         !synthesis of the classes
         call data_synthesis(parallel,debugmode,nproc,ncls,timecls,timecls(1,nproc))

         !calculate the summary of the category
         call sort_positions(ncat,timesum,isort)
!!$         iunit=60
         open(unit=60,file=trim(filename_time),status='unknown',position='append')
!!$                  
!!$         !first get the default stream
!!$         call yaml_get_default_stream(iunit_def)
!!$         if (iunit_def /= iunit) then
!!$            call yaml_set_stream(unit=iunit,tabbing=0,record_length=100,istat=iostat)
!!$            if (iostat /=0) then
!!$               call yaml_set_default_stream(iunit,ierr)
!!$            end if
!!$            !if the stream was not already present just set back the default to iunit_def
!!$         end if
!!$         if (newfile) then
!!$            !start the writing of the file
!!$            call yaml_new_document()
!!$            newfile=.false.
!!$         end if
!!$         call yaml_open_map(trim(message),advance='no')
!!$         if (.not. parallel) then
!!$            call yaml_comment('     % ,  Time (s)')
!!$         else if (debugmode) then
!!$            call yaml_comment('     % ,  Time (s), Load per MPI proc (relative) ')
!!$         else
!!$            call yaml_comment('     % ,  Time (s), Max, Min Load (relative) ')
!!$         end if
!!$         call yaml_open_map('Classes')
!!$         total_pc=0.d0
!!$         do icls=1,ncls
!!$            pc=0.0d0
!!$            if (timesum(ncat+1)/=0.d0) pc=100.d0*timecls(icls,nproc)/timesum(ncat+1)
!!$            total_pc=total_pc+pc
!!$            call yaml_open_sequence(trim(clss(icls)),flow=.true.)
!!$              call yaml_sequence(yaml_toa(pc,fmt='(f5.1)'))
!!$              call yaml_sequence(yaml_toa(timecls(icls,nproc),fmt='(1pg9.2)'))
!!$              do iextra=0,nextra-1
!!$                 call yaml_sequence(yaml_toa(timecls(icls,iextra),fmt='(f5.2)'))
!!$              end do
!!$            call yaml_close_sequence()
!!$         end do
!!$         total_pc=0.d0
!!$         do icls=1,ncls
!!$            pc=0.0d0
!!$            if (timesum(ncat+1)/=0.d0) pc=100.d0*timecls(icls,nproc)/timesum(ncat+1)
!!$            total_pc=total_pc+pc
!!$            write(60,'(4x,a,t21,a,'//trim(formatstring)//')') trim(clss(icls))//':','[',&
!!$                 pc,',',timecls(icls,nproc),&
!!$                 (',',timecls(icls,iextra),iextra=0,nextra-1),']'
!!$         end do
!!$         write(60,'(4x,a,t21,a,'//trim(formatstring)//')') 'Total:','[',&
!!$              total_pc,',',timesum(ncat+1),&
!!$              (',',timeall(ncat+1,iextra),iextra=0,nextra-1),']'
!!$         call yaml_close_map() !classes
!!$
!!$
!!$         call yaml_close_map() !counter
!!$         !restore the default stream
!!$         if (iostat==0) then
!!$            call yaml_set_default_stream(iunit_def,ierr)
!!$         end if

         if (newfile) then
            write(60,'(a)')'---'
            newfile=.false.
         end if
         if (.not. parallel) then
            write(60,'(a,t16,a)')trim(message)//':','   #     % ,  Time (s)' 
         else if (debugmode) then
            write(60,'(a,t16,a)')trim(message)//':','   #     % ,  Time (s), Load per MPI proc (relative) ' 
         else
            write(60,'(a,t16,a)')trim(message)//':','   #     % ,  Time (s), Max, Min Load (relative) ' 
         end if
         !sum all the information by class
         write(60,'(2x,a)')'Classes:'
         total_pc=0.d0
         do icls=1,ncls
            pc=0.0d0
            if (timesum(ncat+1)/=0.d0) pc=100.d0*timecls(icls,nproc)/timesum(ncat+1)
            total_pc=total_pc+pc
            write(60,'(4x,a,t21,a,'//trim(formatstring)//')') trim(clss(icls))//':','[',&
                 pc,',',timecls(icls,nproc),&
                 (',',timecls(icls,iextra),iextra=0,nextra-1),']'
         end do
         write(60,'(4x,a,t21,a,'//trim(formatstring)//')') 'Total:','[',&
                 total_pc,',',timesum(ncat+1),&
                    (',',timeall(ncat+1,iextra),iextra=0,nextra-1),']'
         !Write all relevant categories
         write(60,'(2x,a)')'Categories:'
         do j=1,ncat
            i=isort(j)
            pc=0.d0
            if (timesum(i) /= 0.d0) then
               if (timesum(ncat+1)/=0.d0) pc=100.d0*timesum(i)/timesum(ncat+1)
               write(60,'(4x,a)') trim(cats(1,i))//':'
               write(60,'(t12,a,1x,a,'//trim(formatstring)//')')&
                    ' Data:  ','[',pc,',',timesum(i),&
                    (',',timeall(i,iextra),iextra=0,nextra-1),']'
               write(60,'(t12,a,1x,a)')' Class: ',trim(cats(2,i))
               write(60,'(t12,a,1x,a)')' Info:  ',trim(cats(3,i))
            end if

         enddo
         close(unit=60)
      endif

    END SUBROUTINE sum_results

end module timeData


!> The same timing routine but with system_clock (in case of a supported specs)
subroutine timing(iproc,category,action)
  use timeData

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
  real(kind=8), dimension(ncounters,0:nproc) :: timecnt !< useful only at the very end
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  character(len=MPI_MAX_PROCESSOR_NAME), dimension(0:nproc-1) :: nodename

!$ integer :: omp_get_max_threads

  !first of all, read the time
  !call system_clock(itime,count_rate,count_max)
  call nanosec(itns)

  ! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
  if (action.eq.'IN') then  ! INIT
     !!no need of using system clock for the total time (presumably more than a millisecond)
     !call cpu_time(total0)
     filename_time=repeat(' ',128)
     time0=real(itns,kind=8)*1.d-9
     do i=1,ncat
        timesum(i)=0.d0
        pctimes(i)=0.d0
     enddo
     !in this case iproc stands for nproc
     parallel=abs(iproc) > 1!trim(category).eq.'parallel'
     nproc=abs(iproc)
     filename_time=trim(category)
     newfile=.true.
     init=.false.
     debugmode=(nproc == 2) .or. iproc < -1
     if (nproc >=2) then
        nextra=nproc
        if (.not. debugmode) nextra=2
        write(strextra,'(i5)')nextra
        formatstring='1x,f5.1,a,1x,1pe9.2,a,'//trim(strextra)//'(1x,0pf5.2,a)'
     else
        nextra=0
        formatstring='1x,f5.1,a,1x,1pe9.2,a'
     end if
     ncat_stopped=0 !no stopped category
     ncounters=0

  else if (action.eq.'PR') then !stop partial counters and restart from the beginning
     if (init) then
        print *, 'ERROR: TIMING IS INITIALIZED BEFORE PARTIAL RESULTS'
        stop 
     endif
     !here iproc is the communicator
     if (parallel) then
        call MPI_COMM_RANK(iproc,iproc_true,ierr)
     else
        iproc_true = 0
     end if
     ncounters=ncounters+1
     if (ncounters > ncat) then
        print *, 'It is not allowed to have more partial counters that categories; ncat=',ncat
        stop
     end if
     !name of the category
     pcnames(ncounters)=trim(category)
     !total time elapsed in the category
     timesum(ncat+1)=real(itns,kind=8)*1.d-9-time0
     pctimes(ncounters)=timesum(ncat+1)
     call sum_results(iproc_true,iproc,pcnames(ncounters))
     !reset all timings
     time0=real(itns,kind=8)*1.d-9
     do i=1,ncat
        timesum(i)=0.d0
     enddo

  else if (action.eq.'RE') then ! RESULT
     if (init) then
        print *, 'TIMING IS INITIALIZED BEFORE RESULTS'
        stop 
     endif
     !here iproc is the communicator
     if (parallel) then
        call MPI_COMM_RANK(iproc,iproc_true,ierr)
     else
        iproc_true = 0
     end if

     if (ncounters == 0) then !no partial counters selected
        timesum(ncat+1)=real(itns,kind=8)*1.d-9-time0
        call sum_results(iproc_true,iproc,'ALL')
     else !consider only the results of the partial counters
        if (parallel) then 
           call MPI_GATHER(pctimes,ncounters,MPI_DOUBLE_PRECISION,&
                timecnt,ncounters,MPI_DOUBLE_PRECISION,0,iproc,ierr)
           if (debugmode) then
              !initalise nodenames
              do jproc=0,nproc-1
                 nodename(jproc)=repeat(' ',MPI_MAX_PROCESSOR_NAME)
              end do
              
              call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
              
              !gather the result between all the process
              call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
                   nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
                   iproc,ierr)
           end if

        else
           do i=1,ncounters
              timecnt(i,0)=pctimes(i)
           end do
        endif
          
        if (iproc_true == 0) then
           open(unit=60,file=trim(filename_time),status='unknown',position='append')
           write(60,'(a,t14,a)')'SUMMARY:','   #     % ,  Time (s)'

           !synthesis of the counters
           call data_synthesis(parallel,debugmode,nproc,ncounters,timecnt,timecnt(1,nproc))

           !sum all the information by class
           do i=1,ncounters
              pc=100.d0*timecnt(i,nproc)/sum(timecnt(1:ncounters,nproc))
              write(60,'(2x,a,t19,a,'//trim(formatstring)//')') trim(pcnames(i))//':','[',&
                   pc,',',timecnt(i,nproc),']'
           end do
           write(60,'(2x,a,t19,a,1x,f5.1,a,1x,1pe9.2,a)') 'Total:','[',&
                100.d0,',',sum(timecnt(1:ncounters,nproc)),']'
           !write the number of processors and the number of OpenMP threads
           nthreads = 0
           !$  nthreads=omp_get_max_threads()
           write(60,'(2x,a)')'CPU Parallelism:'
           write(60,'(t10,a,1x,i6)')'MPI procs: ',nproc
           write(60,'(t10,a,1x,i6)')'OMP thrds: ',nthreads
           if (debugmode) then
              write(60,'(t10,a)')'Hostnames:'
              do jproc=0,nproc-1
                 write(60,'(t10,a)')'  - '//trim(nodename(jproc))
              end do
           end if
           close(unit=60)
        end if
     end if
  else
     !controls if the category exists
     catfound=.false.
     do i=1,ncat
        if (trim(category) == trim(cats(1,i))) then
           ii=i
           catfound=.true.
           exit
        endif
     enddo
     if (.not. catfound) then
        print *, 'ACTION  ',action
        write(*,*) 'category, action',category, action
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        stop 'TIMING CATEGORY NOT DEFINED'
     end if

     if (action == 'ON') then  ! ON
        !some other category was initalized before, overriding
!if (iproc==0) print*,'timing on: ',trim(category)
        if (init) return
        t0=real(itns,kind=8)*1.d-9
        init=.true.
        ncaton=ii !category which has been activated
     else if (action == 'OF' .and. ii==ncaton) then  ! OFF
        if (.not. init) then
           print *, cats(1,ii), 'not initialized'
           stop 
        endif
!if (iproc==0) print*,'timing OFF: ',trim(category)
        t1=real(itns,kind=8)*1.d-9
        timesum(ii)=timesum(ii)+t1-t0
        init=.false.
     else if (action == 'OF' .and. ii/=ncaton) then
        if (ncat_stopped /=0) stop 'INTERRUPTS SHOULD NOT BE HALTED BY OF'
!if (iproc==0) print*,'timing2 OFF: ',trim(category)
        !some other category was initalized before, taking that one
        return
    !interrupt the active category and replace it by the proposed one
     else if (action == 'IR') then
        if (ncat_stopped /=0) then
           print *, cats(1,ncat_stopped), 'already exclusively initialized'
           stop
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
        ncaton=ii
        t0=t1

     else if (action == 'RS') then !resume the interrupted category
        if (ncat_stopped ==0) then
           stop 'NOTHING TO RESUME'
        end if
        if (ii /= ncaton) stop 'WRONG RESUMED CATEGORY'
        !time
        t1=real(itns,kind=8)*1.d-9
        timesum(ii)=timesum(ii)+t1-t0
        if (ncat_stopped == -1) then
           init =.false. !restore normal counter
        else
           ncaton=ncat_stopped
           t0=t1
        end if       
        ncat_stopped=0
     else
        print *,action,ii,ncaton,trim(category)
        stop 'TIMING ACTION UNDEFINED'
     endif

  endif

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

!>put the average value of timeall in the timesum array
!then rewrite each element with the deviation from it (in debug mode)
!in normal mode write only the max and min deviations (only in parallel)
subroutine data_synthesis(parallel,debugmode,nproc,ncats,timeall,timesum)
  implicit none
  logical, intent(in) :: parallel,debugmode
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
        else if (parallel) then
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
end subroutine data_synthesis
