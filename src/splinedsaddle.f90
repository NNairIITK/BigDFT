!> @file
!! Contains routines for the executable splined_saddle (search of saddle point)
!! @author
!!    Copyright (C) 2010-2013 BigDFT group (AG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Program splined_saddle: Methot determine the saddle point
program splined_saddle

  use module_base
  use module_types
  use module_interfaces
  use m_ab6_symmetry
  use yaml_output
  implicit none
  character(len=*), parameter :: subname='BigDFT'
  integer :: iproc,nproc,ierr,infocode
  integer :: ncount_bigdft
!!$ logical :: exist_list
  !input/output variables
  type(run_objects) :: runObj
  type(DFT_global_output) :: outs
  character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
  character(len=60) :: run_id
  !character(len=60) :: filename
  ! atomic coordinates, forces
  integer :: iconfig,nconfig!,istat
  real(gp), dimension(:,:), allocatable :: ratsp,fatsp 
  integer, dimension(4) :: mpi_info
  !include 'mpif.h' !non-BigDFT

  call f_lib_initialize()
  !-finds the number of taskgroup size
  !-initializes the mpi_environment for each group
  !-decides the radical name for each run
  call bigdft_init(mpi_info,nconfig,run_id,ierr)

  !just for backward compatibility
  iproc=mpi_info(1)
  nproc=mpi_info(2)


   !allocate arrays of run ids
  allocate(arr_radical(abs(nconfig)))
  allocate(arr_posinp(abs(nconfig)))
   
  !here we call  a routine which
  ! Read a possible radical format argument.
  call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierr)
  
  do iconfig=1,abs(nconfig)
      if (modulo(iconfig-1,mpi_info(4))==mpi_info(3)) then

         ! Read all input files. This should be the sole routine which is called to initialize the run.
         call run_objects_init_from_files(runObj, arr_radical(iconfig),arr_posinp(iconfig))
         call init_global_output(outs, runObj%atoms%astruct%nat)

!!$     !welcome screen
!!$!     if (iproc==0) call print_logo()
!!$
!!$     ! Read all input files.
!!$     !standard names
!!$     call standard_inputfile_names(inputs,radical,nproc)
!!$     call read_input_variables(iproc,nproc,arr_posinp(iconfig),inputs, atoms, rxyz,nconfig,radical,istat)
!!$     !-----------------------------------------------------------
!!$     !-----------------------------------------------------------
!!$     if (iproc == 0) then
!!$        call print_general_parameters(inputs,atoms)
!!$     end if

     open(unit=16,file=trim(runObj%inputs%dir_output)//'geopt.mon',status='unknown',position='append')
     if (iproc ==0 ) write(16,*) '----------------------------------------------------------------------------'


     !initialize memory counting
     !call memocc(0,iproc,'count','start')

     !if other steps are supposed to be done leave the last_run to minus one
     !otherwise put it to one
     !if (inputs%last_run == -1 .and. inputs%ncount_cluster_x <=1) then
     if (runObj%inputs%last_run == -1 .and. runObj%inputs%ncount_cluster_x <=1 .or.  runObj%inputs%ncount_cluster_x <= 1) then
        runObj%inputs%last_run = 1
     end if
 
     call call_bigdft(runObj, outs, nproc,iproc,infocode)

     if (runObj%inputs%ncount_cluster_x > -1) then
        if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
        
        allocate(ratsp(3,runObj%atoms%astruct%nat),fatsp(3,outs%fdim))
        ratsp(1:3,1:runObj%atoms%astruct%nat)=runObj%atoms%astruct%rxyz(1:3,1:runObj%atoms%astruct%nat)
        fatsp(1:3,1:outs%fdim)=outs%fxyz(1:3,1:outs%fdim)
        outs%energy=outs%energy
        call givemesaddle(outs%energy,ratsp,fatsp,16,nproc,iproc,runObj%atoms,runObj%rst,runObj%inputs,ncount_bigdft)
        close(16)
        deallocate(ratsp,fatsp)

        ! geometry optimization
        !call geopt(nproc,iproc,rxyz,atoms,fxyz,etot,rst,inputs,ncount_bigdft)
        !filename=trim('final_'//trim(arr_posinp(iconfig)))
        !if (iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,' ')
     end if

     !if there is a last run to be performed do it now before stopping
     if (runObj%inputs%last_run == -1) then
        runObj%inputs%last_run = 1
        call call_bigdft(runObj, outs, nproc,iproc,infocode)
     end if

     if (iproc == 0) call write_forces(runObj%atoms,outs%fxyz)

     call deallocate_global_output(outs)
     call run_objects_free(runObj, subname)

!!$
!!$     call deallocate_atoms(atoms,subname) 
!!$
!!$!     call deallocate_local_zone_descriptors(rst%Lzd, subname) 
!!$     if(inputs%linear /= INPUT_IG_OFF .and. inputs%linear /= INPUT_IG_LIG) &
!!$          & call deallocateBasicArraysInput(inputs%lin)
!!$
!!$     call free_restart_objects(rst,subname)
!!$
!!$
!!$     call free_input_variables(inputs)
!!$     !-----------------------------------------------------------
!!$
!!$     !finalize memory counting
!!$     call memocc(0,0,'count','stop')
  end if
enddo !loop over iconfig

  deallocate(arr_posinp,arr_radical)

  call bigdft_finalize(ierr)
  call f_lib_finalize()

end program splined_saddle


!> Module for splined saddle
module minimization_sp
    implicit none
    type parameterminimization_sp
        !general parameters for all methods
        integer::ifile=6
        character(10)::approach='unknown'
        real(kind=8)::fmaxtol=-1.d0
        real(kind=8)::eps=1.d-5
        integer::iter=0
        integer::iflag=0
        logical::converged
        real(kind=8)::epotitm1
        real(kind=8)::avgalpha
        real(kind=8)::avgnum
        integer::iprint=1

        !parameters for SD and SDDIIS
        real(kind=8)::alpha=-1.d0
        real(kind=8)::alphax=-1.d0
        integer::maxforcecall=10000
        real(kind=8)::anoise=-1.d0

        !parameters for SD
        integer::itsd
        integer::nitsd=-1
        integer::nsatur=-1
        integer::isatur
        real(kind=8)::alphamin=-1.d0
        real(kind=8)::alphamax=-1.d0
        real(kind=8)::fnrmtolsatur=-1.d0
        real(kind=8)::epotitm2
        real(kind=8)::fnrmitm1
        real(kind=8)::fnrmitm2
        logical::sdsaturated
        logical::sdminimum
        logical::care

        !parameters for FIRE
        real(kind=8)::dt=-1.d0
        real(kind=8)::dtmax=-1.d0
        real(kind=8)::finc=-1.d0
        real(kind=8)::fdec=-1.d0
        real(kind=8)::falpha=-1.d0
        real(kind=8)::alphastart=-1.d0
        integer::ndowntol=-1
        integer::itfire=0

        !parameters for DIIS
        integer::idsx
        integer, allocatable::ipiv(:)
        real(kind=8), allocatable::a(:,:,:),b(:)
        real(kind=8)::emin
        real(kind=8)::fnrmlowest
        logical::diisminimum
        logical::diisdivergence
        integer::itdiis
        integer::ld
        integer::nd

        !parameters for print information
        integer::mp=6
        integer::lp=6

        !parameters for line search routine
        integer::maxfev=20
        real(kind=8)::ftol=1.d-4
        real(kind=8)::gtol=9.d-1
        real(kind=8)::stpmin=1.d-20
        real(kind=8)::stpmax=1.d+20
        !real(kind=8)::xtol=epsilon(xtol)
        integer::info
    end type parameterminimization_sp
end module minimization_sp


!> Module used by the program splined_saddle
module modulesplinedsaddle
    implicit none
    type parametersplinedsaddle 
        !integer, parameter::npmax=20
        !integer::napmax=50
        real(kind=8)::s(0:200)
        real(kind=8)::h(200)
        real(kind=8)::c(0:200)
        real(kind=8)::y(0:200)
        real(kind=8)::e1(200-1)
        real(kind=8)::e2(200-2)
        real(kind=8)::cv(0:50)
        real(kind=8)::ex(0:50)
        real(kind=8)::exd(0:50)
        real(kind=8)::sv(0:50)
        real(kind=8)::e1v(49)
        real(kind=8)::e2v(48)
        real(kind=8)::hv(50)
        real(kind=8)::a(50)
        real(kind=8)::b(50)
        real(kind=8)::tmax
        real(kind=8)::htol=1.d-6  !1.d-6
        real(kind=8)::vdtol=1.d-5  !1.d-10
        real(kind=8)::exends(2)
        real(kind=8)::exends_b(2)
        integer::ifile=6
        integer::ns=2
        integer::ns2=0  !45
        integer::npv
        integer::ncount_ll
        integer::ncount
        real(kind=8)::time_ll
        real(kind=8)::time
        real(kind=8)::epotci
        logical::granot
        character(20)::hybrid
        character(20)::doneb
        character(20)::docineb
        character(20)::pickbestanchorpoints
        logical::do_fill_ex_exd
        character(20)::runstat
        character(10)::typintpol
    end type parametersplinedsaddle
end module modulesplinedsaddle


subroutine givemesaddle(epot_sp,ratsp,fatsp,ifile,nproc,iproc,atoms,rst,inputs,ncount_bigdft)
  use dictionaries
    use module_base
    use module_interfaces
    use module_types
    use module_input_keys, only:DFT_VARIABLES, KPT_VARIABLES
    use minimization_sp, only:parameterminimization_sp  !Reza
    use modulesplinedsaddle, only:parametersplinedsaddle
    use module_input_dicts
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    real(gp), dimension(:,:), allocatable :: rxyz_2
    real(gp), dimension(:,:), allocatable :: x,f,xneb,rxyz_tmp,x_t
    !type(DFT_global_output), dimension(:), allocatable :: outs
    integer :: np,np_neb,np_t,iat,ifile
    type(parameterminimization_sp)::parmin_neb,parmin
    real(gp) ::epot_sp,ratsp(3,atoms%astruct%nat),fatsp(3,atoms%astruct%nat)
    character(len=20) :: tatonam
    integer::n,nr,istat,infocode,ixyz,i,mm1,mm2,mm3
    real(kind=8)::fnrm,fnrm1,fnrm2,tt1,tt2,tt3,time1,time2
    type(parametersplinedsaddle)::pnow
    type(dictionary), pointer :: dict
    type(input_variables), target :: ll_inputs
    type(atoms_data) :: ll_atoms
    type(run_objects) :: ll_runObj, runObj
    type(DFT_global_output), dimension(2) :: outends
    !character(50)::ssm
    character(len=20)::filename
    logical::move_this_coordinate
    character(40)::comment
    integer, parameter::ndeb1=0,ndeb2=0
    !---------------------------------------------------------------------------
    !pnow%ncount=1
    !pnow%ncount_ll=0
    nullify(dict)
    ncount_bigdft=0
    pnow%ifile=ifile
    parmin%ifile=ifile
    parmin_neb%ifile=ifile
    n=3*atoms%astruct%nat
    nr=0
    pnow%time_ll=0.0d0
    do i=1,3*atoms%astruct%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) nr=nr+1
    enddo
    if(iproc==0) write(*,'(a,i0,1x,i0)') 'DOF: n,nr ',n,nr
    !---------------------------------------------------------------------------
    call readinputsplsad(iproc,np,np_neb,parmin,parmin_neb,pnow)
    if(iproc==0)    then
        write(*,*) 'number of anchor points ',np
        write(*,*) 'degree of freedom: n,nr ',n,nr
    endif
    !-----------------------------------------------------------
    ! We read the input variable files.
    if(trim(pnow%hybrid)=='yes') then
       call read_input_dict_from_files("ll_input", bigdft_mpi,dict)
    else
       call read_input_dict_from_files("input", bigdft_mpi,dict)
    endif
    ! We add the atomic data.
    call astruct_merge_to_dict(dict // "posinp", atoms%astruct, atoms%astruct%rxyz)
    call atoms_file_merge_to_dict(dict)
    call atomic_data_file_merge_to_dict(dict, "Atomic occupation", "input.occup")
    ! We parse the dictionary.
    call inputs_from_dict(ll_inputs, ll_atoms, dict)
    call dict_free(dict)
    call deallocate_atoms(ll_atoms, "givemesaddle")
    
    !-----------------------------------------------------------
    allocate(rxyz_2(3,atoms%astruct%nat+ndeb1))
    call dmemocc(3*(atoms%astruct%nat),3*(atoms%astruct%nat+ndeb1),rxyz_2,'rxyz_2')
    allocate(rxyz_tmp(3,atoms%astruct%nat+ndeb1))
    call dmemocc(3*(atoms%astruct%nat),3*(atoms%astruct%nat+ndeb1),rxyz_tmp,'rxyz_tmp')
    allocate(f(n,0:np+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating f'
    call dmemocc(n*(np+1),n*(np+1+ndeb2),f,'f')
    allocate(x(n,0:np+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating x'
    call dmemocc(n*(np+1),n*(np+1+ndeb2),x,'x')
    allocate(xneb(n,0:np_neb+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xneb'
    call dmemocc(n*(np_neb+1),n*(np_neb+1+ndeb2),xneb,'xneb')
    call init_global_output(outends(1), atoms%astruct%nat)
    call init_global_output(outends(2), atoms%astruct%nat)
    !if(iproc==0) write(*,*) 'ALIREZA-01'
    !---------------------------------------------------------------------------
    if(trim(pnow%runstat)=='restart') then
        allocate(x_t(n,0:100+ndeb2),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating x_t'
        call dmemocc(n*(100+1),n*(100+1+ndeb2),x_t,'x_t')
        filename='anchorposinp.xyz' 
        call readanchorpoints(n,np_t,x_t,filename,atoms)
        if(np_t==np) then
            x(1:n,0:np)=x_t(1:n,0:np_t)
        else
            call change_np(n,np_t,x_t,atoms,np,x)
        endif
        deallocate(x_t,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x_t'
    else
        open(unit=1336,file='posinp.xyz',status='old') !read atomic positions
        read(1336,*) 
        read(1336,*) 
        do iat=1,atoms%astruct%nat
          read(1336,*) tatonam,rxyz_tmp(1,iat),rxyz_tmp(2,iat),rxyz_tmp(3,iat)
          if(atoms%astruct%units == 'angstroemd0' .or. atoms%astruct%units == 'angstroem') then
            rxyz_tmp(1:3,iat)=rxyz_tmp(1:3,iat)/0.5291772108_gp !non-BigDFT
          endif
        enddo
        close(1336)
        open(unit=1336,file='posinp2.xyz',status='old') !read atomic positions
        read(1336,*) 
        read(1336,*) 
        do iat=1,atoms%astruct%nat
          read(1336,*) tatonam,rxyz_2(1,iat),rxyz_2(2,iat),rxyz_2(3,iat)
          if(atoms%astruct%units == 'angstroemd0' .or. atoms%astruct%units == 'angstroem') then
            rxyz_2(1:3,iat)=rxyz_2(1:3,iat)/0.5291772108_gp !non-BigDFT
          endif
        enddo
        close(1336)
        rxyz_2(1:3,1:atoms%astruct%nat)=rxyz_2(1:3,1:atoms%astruct%nat)&
             +ratsp(1:3,1:atoms%astruct%nat)-rxyz_tmp(1:3,1:atoms%astruct%nat)
        if(trim(pnow%doneb)=='yes') then
            call initializepoints(atoms,n,ratsp,rxyz_2,np_neb,xneb)
            call change_np(n,np_neb,xneb,atoms,np,x)
        else
            call initializepoints(atoms,n,ratsp,rxyz_2,np,x)
        endif
    endif
    if(iproc==0) call writepathway(n,np,x,'pathinp.xyz',atoms)
    if(iproc==0) write(*,*) 'ALIREZA-02'
    !---------------------------------------------------------------------------
    call vcopy(n, fatsp(1,1), 1, outends(1)%fxyz(1,1), 1)
    pnow%exends_b(1)=epot_sp
    call atomic_dot(atoms,outends(1)%fxyz(1,1),outends(1)%fxyz(1,1),fnrm1)
    fnrm1=sqrt(fnrm1)
    call run_objects_nullify(runObj)
    call run_objects_associate(runObj, inputs, atoms, rst, x(1,np))
    !if(iproc==0) write(*,*) 'ALIREZA-03'
    call cpu_time(time1)
    call call_bigdft(runObj,outends(2),nproc,iproc,infocode)
    call cpu_time(time2)
    call run_objects_free_container(runObj)
    pnow%exends_b(2)=outends(2)%energy
    ncount_bigdft=ncount_bigdft+1
    pnow%ncount=2
    pnow%time=2.d0*(time2-time1)
    call atomic_dot(atoms,outends(2)%fxyz(1,1),outends(2)%fxyz(1,1),fnrm2)
    fnrm2=sqrt(fnrm2)
    !if(iproc==0) write(*,*) 'ALIREZA-04'
    if(iproc==0) then
        write(pnow%ifile,'(a,4e24.15)') 'ENDs: epot1,fnrm1,epot2,fnrm2 ', &
            pnow%exends_b(1),fnrm1,pnow%exends_b(2),fnrm2
        write(*,'(a,4e24.15)') 'ENDs: epot1,fnrm1,epot2,fnrm2 ', &
            pnow%exends_b(1),fnrm1,pnow%exends_b(2),fnrm2
    endif
    !---------------------------------------------------------------------------
    if(trim(pnow%hybrid)=='yes') then
       call run_objects_nullify(ll_runObj)
       call run_objects_associate(ll_runObj, ll_inputs, atoms, rst, x(1,0))
       call cpu_time(time1)
       call call_bigdft(ll_runObj,outends(1),nproc,iproc,infocode)
       pnow%exends(1) = outends(1)%energy
       call cpu_time(time2)
       ncount_bigdft=ncount_bigdft+1
       pnow%ncount_ll=1
       pnow%time_ll=pnow%time_ll+(time2-time1)
       call vcopy(n, x(1,np), 1, ll_runObj%atoms%astruct%rxyz(1,1), 1)
       call cpu_time(time1)
       call call_bigdft(ll_runObj,outends(2),nproc,iproc,infocode)
       call cpu_time(time2)
       call run_objects_free_container(ll_runObj)
       pnow%exends(2) = outends(2)%energy
       ncount_bigdft=ncount_bigdft+1
       pnow%ncount_ll=pnow%ncount_ll+1
       pnow%time_ll=pnow%time_ll+(time2-time1)
    else
        pnow%exends(1)=pnow%exends_b(1)
        pnow%exends(2)=pnow%exends_b(2)
    endif
    !---------------------------------------------------------------------------
    !if(trim(pnow%runstat)=='new') then
    if(trim(pnow%doneb)=='yes') then
        if(trim(pnow%runstat)=='restart') call change_np(n,np,x,atoms,np_neb,xneb)
        parmin_neb%alphax=1.d0*parmin%alphax !non-BigDFT
        parmin_neb%alphamin=5.d-2*parmin_neb%alphax
        parmin_neb%alphamax=2.0d0*parmin_neb%alphax
        parmin_neb%approach='BFGS' !SD or SDDIIS
        parmin_neb%alpha=1.d0*parmin%alphax
        call initminimize(parmin_neb)
        parmin_neb%maxforcecall=200 !30  !10
        parmin_neb%fnrmtolsatur=1.d-4 !5.d-2
        pnow%ex(0)=pnow%exends(1)
        pnow%ex(np_neb)=pnow%exends(2)
        call neb(n,nr,np_neb,xneb,parmin_neb,pnow, &
            nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)  
        call finalminimize(parmin_neb)
        call change_np(n,np_neb,xneb,atoms,np,x)
    endif
    !------------------------------------------------------
    if(trim(pnow%pickbestanchorpoints)=='yes') then
        !call pickbestanchors(n,np,x,fends,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
        call pickbestanchors2(n,np,x,outends,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    endif
    !call improvepeak(n,nr,np,x,fends,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    !------------------------------------------------------
    !------------------------------------------------------
    if(trim(pnow%docineb)=='no') then
    parmin_neb%stpmax=10.d0
    parmin_neb%eps=1.d-8
    parmin_neb%ftol=1.d-8
    parmin_neb%gtol=9.9d-1
    parmin%fnrmtolsatur=2.d-2 !5.d-2
    parmin%alphamin=1.d-1*parmin%alphax
    parmin%alphamax=3.d0*parmin%alphax
    call initminimize(parmin)
    call splinedsaddle(n,nr,np,x,epot_sp,f,ratsp,parmin,outends,pnow, & 
        nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft,fatsp)  
    endif
    if(iproc==0) call writepathway(n,np,x,'pathout.xyz',atoms)
    call finalminimize(parmin)
    if (iproc==0) then
        call atomic_dot(atoms,fatsp,fatsp,fnrm);fnrm=sqrt(fnrm)
       write(comment,'(a,1pe10.3)')'CONJG:fnrm= ',fnrm
       call  write_atomic_file('saddle',epot_sp,ratsp,atoms,trim(comment))
    endif
    if(iproc==0) then
        mm1=pnow%ncount
        mm2=pnow%ncount_ll
        tt1=pnow%time/pnow%ncount
        tt2=pnow%time_ll/pnow%ncount_ll
        tt3=tt1/tt2
        mm3=mm1+int(real(mm2,8)/tt3)
        write(*,'(a,3i5,3es15.5)') 'SP-TIMINGS: ',mm1,mm2,mm3,tt1,tt2,tt3
    endif
    !if(iproc==0) then
    !    write(55,*) atoms%astruct%nat
    !    write(55,*)
    !    do iat=1,atoms%astruct%nat
    !        write(55,'(a,3e24.15)') ' Si ',ratsp(1,iat),ratsp(2,iat),ratsp(3,iat)
    !    enddo
    !endif
     !-----------------------------------------------------------
     call free_input_variables(ll_inputs)
     !-----------------------------------------------------------
    deallocate(f,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating f'
    deallocate(x,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x'
    deallocate(xneb,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xneb'
    deallocate(rxyz_2,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating rxyz_2'
    deallocate(rxyz_tmp,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating rxyz_tmp'
    call deallocate_global_output(outends(1))
    call deallocate_global_output(outends(2))
end subroutine givemesaddle


subroutine change_np(n,np1,x1,atoms,np2,x2)
    use module_types
    !use modulesplinedsaddle, only:parametersplinedsaddle
    !use energyandforces, only:calenergyforces
    implicit none
    integer::n,np1,np2,i,ip,mp,iat,ixyz
    real(kind=8)::x1(n,0:np1),x2(n,0:np2)
    real(kind=8)::s(0:100),h(100),y(0:100),e1(200-1),e2(200-2),c(0:200),tt,dt,ed_tt,edd_tt
    type(atoms_data), intent(inout) :: atoms
    logical::move_this_coordinate
    !x_t(1:n,0:np)=x(1:n,0:np)
    !deallocate(x,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x'
    !allocate(x(n,0:50+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating x'
    !call dmemocc(n*(50+1),n*(50+1+ndeb2),x,'x')
    if(np1>100 .or. np2>100) stop 'ERROR: np1>100 .or. np2>100'
    if(np1==np2) then
        x2(1:n,0:np2)=x1(1:n,0:np1)
    else
        call equalarclengthparametrization(atoms,n,np1,x1,s,h)
        call factor_cubic(np1,h,e1,e2)
        dt=s(np1)/real(np2,8)
        x2(1:n,0)=x1(1:n,0)
        x2(1:n,np2)=x1(1:n,np1)
        do i=1,n
            !if(i<=nr) then
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                y(0:np1)=x1(i,0:np1)
                call inter_cubic(np1,y,h,e1,e2,c)
                do ip=1,np2-1
                    tt=dt*ip
                    call calindex(np1,s,tt,mp,'change_np')
                    call ffdfdd_cubic(np1,y,s,mp,h(mp),tt,c,x2(i,ip),ed_tt,edd_tt)
                enddo
            else
                x2(i,1:np2-1)=x1(i,0)
            endif
        enddo
    endif
    !np=np2
    !deallocate(f,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating f'
    !allocate(f(n,0:np+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating f'
    !call dmemocc(n*(np+1),n*(np+1+ndeb2),f,'f')
end subroutine change_np


subroutine improvepeak(n,nr,np,x,outends,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    !use energyandforces, only:calenergyforces
    implicit none
    integer::n,nr,np,i,ip,istat,npv,nproc,iproc,mp,lp,iat,ixyz,iter,ncount_bigdft,infocode
    type(DFT_global_output), dimension(2), intent(in) :: outends
    real(kind=8)::x(n,0:np),time1,time2 !,f(n,0:np),calnorm
    real(kind=8)::ed_tt,edd_tt,tarr(100),diff,proj,fnrm !n(c) dt
    real(kind=8), allocatable::xt(:)
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    type(parametersplinedsaddle)::pnow,pold
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    logical::move_this_coordinate
    type(run_objects) :: runObj
    type(DFT_global_output) :: outs

    if(mod(np+pnow%ns2,2)==0) then
        npv=np+pnow%ns2+4
    else
        npv=np+pnow%ns2+3
    endif
    allocate(xt(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xt.'
    call dmemocc(n,n+ndeb1,xt,'xt')
    call init_global_output(outs, n / 3)
    call equalarclengthparametrization(atoms,n,np,x,pnow%s,pnow%h)
    call factor_cubic(np,pnow%h,pnow%e1,pnow%e2)
    call fill_ex_exd(0,n,np,x,outends,npv,pnow,pold,xt,outs%fxyz,nproc,iproc,atoms,rst,ll_inputs,&
         ncount_bigdft)
    !call guessinitialtmax_hermite(npv,pnow)
    call guessinitialtmax_cubic(npv,pnow)
    !call calindex(np,pnow%s,8.8165d-01,ip) !CAUTIOUS
    !n(c) dt=pnow%s(np)/np
    diff=1.d10
    do ip=1,np-1
        tarr(ip)=pnow%s(ip)
        if(abs(tarr(ip)-pnow%tmax)<diff) then
            diff=abs(tarr(ip)-pnow%tmax)
            mp=ip
        endif
    enddo
    if(iproc==0) then
        write(*,*) 'MP ',mp,pnow%tmax/pnow%s(np)
    endif
    call calindex(np,pnow%s,pnow%tmax,lp,'improvepeak')
    do i=1,n
        !if(i<=nr) then
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            pnow%y(0:np)=x(i,0:np)
            call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
            call ffdfdd_cubic(np,pnow%y,pnow%s,lp,pnow%h(lp),pnow%tmax,pnow%c,x(i,mp),ed_tt,edd_tt)
        else
            x(i,1:np-1)=x(i,0)
        endif
    enddo
    call run_objects_nullify(runObj)
    call run_objects_associate(runObj, ll_inputs, atoms, rst)
    lp=mp+1
    do iter=1,10
        !call calenergyforces(n,x(1,lp),epot,ft)
        call cpu_time(time1)
        call vcopy(n, x(1,lp), 1, runObj%atoms%astruct%rxyz(1,1), 1)
        call call_bigdft(runObj,outs,nproc,iproc,infocode)
        call cpu_time(time2)
        ncount_bigdft=ncount_bigdft+1
        pnow%ncount_ll=pnow%ncount_ll+1
        pnow%time_ll=pnow%time_ll+(time2-time1)
        !fnrm=DNRM2(2,ft,1)
        call atomic_dot(atoms,outs%fxyz(1,1),outs%fxyz(1,1),fnrm);fnrm=sqrt(fnrm)
        xt(1:n)=x(1:n,lp)-x(1:n,mp)
        call normalizevector2(nr,xt)
        !proj=DDOT(nr,ft,1,xt,1)
        call atomic_dot(atoms,outs%fxyz(1,1),xt,proj)
        if(iproc==0) then
            write(*,*) 'REZA ',x(1,3),x(2,3)
            write(*,*) 'proj ',proj,sqrt(fnrm**2-proj**2)
        endif
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
               !ft(i)=ft(i)-proj*xt(i) !-1.d0*xt(i)
               outs%fxyz(ixyz, iat) = outs%fxyz(ixyz, iat) - proj * xt(i)
               !x(i,lp)=x(i,lp)+1.d-1*ft(i)
               x(i,lp)=x(i,lp)+1.d-1*outs%fxyz(ixyz, iat)
            endif
        enddo
    enddo
    call run_objects_free_container(runObj)
    deallocate(xt,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xt.'
    call deallocate_global_output(outs)
end subroutine improvepeak


subroutine pickbestanchors2(n,np,x,outends,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    !use energyandforces, only:calenergyforces
    implicit none
    integer::n,np,i,ip,istat,npv,nproc,iproc,mp,ncount_bigdft,ixyz,iat,icycle,ncycle
    type(DFT_global_output), dimension(2), intent(in) :: outends
    real(kind=8)::x(n,0:np) !,f(n,0:np),calnorm
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    real(kind=8)::tt,t1,t2,ed_tt,edd_tt,dt
    real(kind=8), allocatable::xt(:),ft(:)
    type(parametersplinedsaddle)::pnow,pold
    real(kind=8)::ttmin,ttmax,emin,e1,e2,exo(0:100),exn(0:100),area,areatot
    real(kind=8)::s_t(0:100)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    logical::move_this_coordinate
    if(mod(np+pnow%ns2,2)==0) then
        npv=np+pnow%ns2+4
    else
        npv=np+pnow%ns2+3
    endif
    allocate(xt(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xt.'
    call dmemocc(n,n+ndeb1,xt,'xt')
    allocate(ft(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ft.'
    call dmemocc(n,n+ndeb1,ft,'ft')
    call equalarclengthparametrization(atoms,n,np,x,pnow%s,pnow%h)
    call factor_cubic(np,pnow%h,pnow%e1,pnow%e2)
    call fill_ex_exd(0,n,np,x,outends,npv,pnow,pold,xt,ft,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    !call guessinitialtmax_hermite(npv,pnow)
    !-------------------------------------------------------------
    pold=pnow
    pold%npv=npv
    npv=np
    pnow%sv(0)=0.d0
    pnow%sv(npv)=pnow%s(np)
    emin=min(pold%ex(0),pold%ex(pold%npv))
    exo(0:pold%npv)=pold%ex(0:pold%npv)-emin
    areatot=0.d0
    do ip=1,pold%npv
        areatot=areatot+0.5d0*(exo(ip-1)+exo(ip))*(pold%sv(ip)-pold%sv(ip-1))
    enddo
    dt=pold%sv(pold%npv)/npv
    pnow%sv(0)=0.d0
    do ip=1,npv-1
        pnow%sv(ip)=dt*ip
    enddo

    ncycle=100
    exn(0)=exo(0)
    exn(npv)=exo(pold%npv)
    do icycle=1,ncycle
        pnow%sv(npv)=pold%sv(pold%npv) !-1.d-10
        do ip=1,npv-1
            call calindex(pold%npv,pold%sv,pnow%sv(ip),mp,'estimate_sv')
            t1=pold%sv(mp-1) ; e1=exo(mp-1)
            t2=pold%sv(mp  ) ; e2=exo(mp  )
            exn(ip)=(e2-e1)/(t2-t1)*(pnow%sv(ip)-t1)+e1
        enddo
        do ip=1,npv
            area=0.5d0*(exn(ip-1)+exn(ip))*(pnow%sv(ip)-pnow%sv(ip-1))
            e1=(exn(npv)-exn(0))/pnow%sv(npv)*pnow%sv(ip-1)+exn(0)
            e2=(exn(npv)-exn(0))/pnow%sv(npv)*pnow%sv(ip  )+exn(0)
            !area=area-0.5d0*(e1+e2)*(pnow%sv(ip)-pnow%sv(ip-1))
            if(area<areatot/npv) then
                pnow%hv(ip)=(pnow%sv(ip)-pnow%sv(ip-1))*1.02d0
            else
                pnow%hv(ip)=(pnow%sv(ip)-pnow%sv(ip-1))*0.98d0
            endif
            !pnow%hv(ip)=1.d0/area
        enddo
        !pnow%sv(0)=0.d0
        do ip=1,npv
            pnow%sv(ip)=pnow%sv(ip-1)+pnow%hv(ip)
        enddo
        tt=pnow%sv(npv)
        pnow%sv(0:npv)=pnow%sv(0:npv)*pold%sv(pold%npv)/tt
        do ip=1,npv
            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
        enddo
        ttmax=1.50d0*pold%sv(pold%npv)/npv
        ttmin=0.70d0*pold%sv(pold%npv)/npv
        do ip=1,npv
            tt=max(min(pnow%hv(ip),ttmax),ttmin)
            pnow%sv(ip)=pnow%sv(ip-1)+tt
        enddo
        tt=pnow%sv(npv)
        pnow%sv(0:npv)=pnow%sv(0:npv)*pold%sv(pold%npv)/tt
        do ip=1,npv
            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
        enddo
    enddo

    do ip=1,npv-1
        call calindex(pold%npv,pold%sv,pnow%sv(ip),mp,'estimate_sv')
        t1=pold%sv(mp-1) ; e1=exo(mp-1)
        t2=pold%sv(mp  ) ; e2=exo(mp  )
        exn(ip)=(e2-e1)/(t2-t1)*(pnow%sv(ip)-t1)+e1
    enddo
    exn(0:npv)=exn(0:npv)+emin

    tt=pnow%sv(npv)
    s_t(0:npv)=pnow%sv(0:npv)*pnow%s(np)/tt
    !do ip=1,npv
    !    pnow%hv(ip)=s_t(ip)-s_t(ip-1)
    !enddo
    !-------------------------------------------------------------
    !dt=pnow%s(np)/np
    do i=1,n
        !if(i<=nr) then
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            pnow%y(0:np)=x(i,0:np)
            call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
            !diff=1.d10
            !do ip=1,np-1
            !    tarr(ip)=dt*ip
            !    if(abs(tarr(ip)-pnow%tmax)<diff) then
            !        diff=abs(tarr(ip)-pnow%tmax)
            !        mp=ip
            !    endif
            !enddo
            !tarr(mp)=pnow%tmax
            do ip=1,np-1
                call calindex(np,pnow%s,s_t(ip),mp,'pickbestanchors2')
                call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),s_t(ip),pnow%c,x(i,ip),ed_tt,edd_tt)
            enddo
        else
            x(i,1:np-1)=x(i,0)
        endif
    enddo
    deallocate(xt,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xt.'
    deallocate(ft,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ft.'
end subroutine pickbestanchors2


subroutine pickbestanchors(n,np,x,outends,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    !use energyandforces, only:calenergyforces
    implicit none
    integer::n,np,i,ip,istat,npv,nproc,iproc,mp,ncount_bigdft,ixyz,iat
    type(DFT_global_output), dimension(2), intent(in) :: outends
    real(kind=8)::x(n,0:np) !,f(n,0:np),calnorm
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    real(kind=8)::ed_tt,edd_tt,tarr(100),diff,dt
    real(kind=8), allocatable::xt(:),ft(:)
    type(parametersplinedsaddle)::pnow,pold
    integer, parameter :: ndeb1=0
    !integer, parameter :: ndeb2=0
    logical::move_this_coordinate
    if(mod(np+pnow%ns2,2)==0) then
        npv=np+pnow%ns2+4
    else
        npv=np+pnow%ns2+3
    endif
    allocate(xt(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xt.'
    call dmemocc(n,n+ndeb1,xt,'xt')
    allocate(ft(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ft.'
    call dmemocc(n,n+ndeb1,ft,'ft')
    call equalarclengthparametrization(atoms,n,np,x,pnow%s,pnow%h)
    call factor_cubic(np,pnow%h,pnow%e1,pnow%e2)
    call fill_ex_exd(0,n,np,x,outends,npv,pnow,pold,xt,ft,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    !call fill_ex_exd(istep,n,np,x,fends,npv,pnow,pold,xt,ft,nproc,iproc,atoms,rst,inputs,ncount_bigdft)
    call guessinitialtmax_hermite(npv,pnow)
    !call calindex(np,pnow%s,8.8165d-01,ip)
    dt=pnow%s(np)/np
    do i=1,n
        !if(i<=nr) then
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            pnow%y(0:np)=x(i,0:np)
            call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
            diff=1.d10
            do ip=1,np-1
                tarr(ip)=dt*ip
                if(abs(tarr(ip)-pnow%tmax)<diff) then
                    diff=abs(tarr(ip)-pnow%tmax)
                    mp=ip
                endif
            enddo
            tarr(mp)=pnow%tmax
            do ip=1,np-1
                call calindex(np,pnow%s,tarr(ip),mp,'pickbestanchors')
                call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),tarr(ip),pnow%c,x(i,ip),ed_tt,edd_tt)
            enddo
        else
            x(i,1:np-1)=x(i,0)
        endif
    enddo
    deallocate(xt,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xt.'
    deallocate(ft,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ft.'
end subroutine pickbestanchors


subroutine readinputsplsad(iproc,np,np_neb,parmin,parmin_neb,pnow)
    use minimization_sp, only:parameterminimization_sp
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    type(parameterminimization_sp)::parmin_neb,parmin
    type(parametersplinedsaddle)::pnow
    integer::iproc,np,np_neb,ios,iline,ich,ios_t,ios_open
    character(256)::strline
    character(20)::command
    character(1)::ch
    character(10)::typintpol
    character(256)::str1
    pnow%runstat='new'
    parmin_neb%fmaxtol=2.d-2
    parmin%approach='SD'
    parmin%alphax=0.5d0
    parmin%fmaxtol=2.d-4
    parmin%maxforcecall=100
    parmin%dt=0.03d0
    np=3
    np_neb=-1
    pnow%ns2=0
    pnow%vdtol=1.d-1
    pnow%htol=2.d-2
    pnow%typintpol='cubic'
    pnow%pickbestanchorpoints='unknown'
    pnow%doneb='unknown'
    pnow%docineb='no'
    pnow%hybrid='no'
    open(unit=1376,file='input.splsad',status='old',iostat=ios_open)
    if(ios_open/=0) then
        if(iproc==0) write(*,*) &
            'WARNING: input.splsad not prived so all parameters are set to default'
    else
        do iline=1,1000
            strline=''
            do ich=1,256
                read(1376,'(a1)',advance='no',iostat=ios) ch
                !write(21,*) iline,ich,ios,ch
                if(ios<0) exit
                strline(ich:ich)=ch
            enddo
            if(ios==-1) exit
            strline(ich:ich)=achar(32)
            command=''
            read(strline,*,iostat=ios_t) command
            ch=command(1:1)
            if(ch=='#') cycle
            if(command=='runstat') then
                read(strline,*,iostat=ios_t) command,pnow%runstat
                if(trim(pnow%runstat)=='new') then
                    if(trim(pnow%pickbestanchorpoints)=='unknown') then
                        pnow%pickbestanchorpoints='yes'
                    endif
                    if(trim(pnow%doneb)=='unknown') then
                        pnow%doneb='yes'
                    endif
                endif
            elseif(command=='fmaxtolneb') then
                read(strline,*,iostat=ios_t) command,parmin_neb%fmaxtol
            elseif(command=='approach') then
                read(strline,*,iostat=ios_t) command,parmin%approach
            elseif(command=='alphax') then
                read(strline,*,iostat=ios_t) command,parmin%alphax
            elseif(command=='fmaxtol') then
                read(strline,*,iostat=ios_t) command,parmin%fmaxtol
            elseif(command=='maxforcecall') then
                read(strline,*,iostat=ios_t) command,parmin%maxforcecall
            elseif(command=='np') then
                read(strline,*,iostat=ios_t) command,np
            elseif(command=='npneb') then
                read(strline,*,iostat=ios_t) command,np_neb
            elseif(command=='ns2') then
                read(strline,*,iostat=ios_t) command,pnow%ns2
            elseif(command=='vdtol') then
                read(strline,*,iostat=ios_t) command,pnow%vdtol
            elseif(command=='dt') then
                read(strline,*,iostat=ios_t) command,parmin%dt
            elseif(command=='htol') then
                read(strline,*,iostat=ios_t) command,pnow%htol
            elseif(command=='hybrid') then
                read(strline,*,iostat=ios_t) command,pnow%hybrid
            elseif(command=='doneb') then
                read(strline,*,iostat=ios_t) command,pnow%doneb
            elseif(command=='docineb') then
                read(strline,*,iostat=ios_t) command,pnow%docineb
            elseif(command=='pickbestanchorpoints') then
                read(strline,*,iostat=ios_t) command,pnow%pickbestanchorpoints
            elseif(command=='typintpol') then
                read(strline,*,iostat=ios_t) command,typintpol
                if(trim(typintpol)=='quintic') then
                    pnow%typintpol='quintic'
                elseif(trim(typintpol)=='cubic') then
                    pnow%typintpol='cubic'
                else
                    pnow%typintpol='cubic'
                    str1='WARNING: not a correct keyword for typintpol, set to default'
                    if(iproc==0) then
                        write(pnow%ifile,*) trim(str1)
                        write(*         ,*) trim(str1)
                    endif
                endif
            endif
            !write(*,'(iline,1x,a,3f)') ios,command,a,b,c
        enddo
        close(1376)
    endif
    if(np_neb==-1) np_neb=np
    if(trim(pnow%runstat)=='new') then
        if(trim(pnow%pickbestanchorpoints)=='unknown') then
            pnow%pickbestanchorpoints='yes'
        endif
        if(trim(pnow%doneb)=='unknown') then
            pnow%doneb='yes'
        endif
    else
        if(trim(pnow%pickbestanchorpoints)=='unknown') then
            pnow%pickbestanchorpoints='no'
        endif
        if(trim(pnow%doneb)=='unknown') then
            pnow%doneb='no'
        endif
    endif
    if(trim(pnow%docineb)=='yes') then
        pnow%doneb='yes'
        pnow%pickbestanchorpoints='no'
    endif
    if(iproc==0) then
        write(*,*) '------------ parameters of splined saddle method ----------'
        write(*,*) 'SPINFO: runstat ',trim(pnow%runstat)
        write(*,*) 'SPINFO: hybrid ',trim(pnow%hybrid)
        write(*,*) 'SPINFO: fmaxtolneb ',parmin_neb%fmaxtol
        write(*,*) 'SPINFO: approach ',parmin%approach
        write(*,*) 'SPINFO: alphax ',parmin%alphax
        write(*,*) 'SPINFO: fmaxtol ',parmin%fmaxtol
        write(*,*) 'SPINFO: maxforcecall ',parmin%maxforcecall
        write(*,*) 'SPINFO: dt ',parmin%dt
        write(*,*) 'SPINFO: np ',np
        write(*,*) 'SPINFO: npneb ',np_neb
        write(*,*) 'SPINFO: ns2 ',pnow%ns2
        write(*,*) 'SPINFO: vdtol ',pnow%vdtol
        write(*,*) 'SPINFO: htol ',pnow%htol
        write(*,*) 'SPINFO: typintpol ',trim(pnow%typintpol)
        write(*,*) 'SPINFO: doneb ',trim(pnow%doneb)
        write(*,*) 'SPINFO: docineb ',trim(pnow%docineb)
        write(*,*) 'SPINFO: pickbestanchorpoints ',trim(pnow%pickbestanchorpoints)
        write(*,*) '-----------------------------------------------------------'
    endif
end subroutine readinputsplsad


subroutine neb(n,nr,np,x,parmin,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use minimization_sp, only:parameterminimization_sp
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    integer::n,nr,np,ip,icall,istat,it,nwork,nra
    real(kind=8)::x(n,0:np)
    real(kind=8)::fnrm,fspmax,fnrmtot
    real(kind=8), allocatable::work(:)
    real(kind=8), allocatable::xa(:,:),fa(:,:)
    type(parameterminimization_sp)::parmin
    type(parametersplinedsaddle)::pnow !n(c) pold
    type(DFT_global_output), dimension(1:np - 1) :: outs
    
    integer, parameter::ndeb1=0,ndeb2=0
    !logical::move_this_coordinate
    parmin%converged=.false.
    if(iproc==0) then
    write(pnow%ifile,'(a,1x,a)') 'begin of minimization_sp using ',parmin%approach
    write(*,'(a,1x,a)') 'begin of minimization_sp using ',parmin%approach
    endif
    if(parmin%approach=='unknown') then
        if(iproc==0) then
        write(pnow%ifile,*) 'The minimize routine returns becuase method is not specified.'
        write(*,*) 'The minimize routine returns becuase method is not specified.'
        endif
        return
    endif
    !-------------------------------------------------------------------------------------
    allocate(xa(nr,np-1+ndeb2),fa(nr,np-1+ndeb2),stat=istat)
    call dmemocc(nr*(np-1),nr*(np-1+ndeb2),xa,'xa')
    call dmemocc(nr*(np-1),nr*(np-1+ndeb2),fa,'fa')
    if(istat/=0) stop 'ERROR: failure allocating xa or fa.'
    do ip=1,np-1
       call init_global_output(outs(ip), atoms%astruct%nat)
        call atomic_copymoving_forward(atoms,n,x(1,ip),nr,xa(1,ip))
    enddo
    !xa(1:nr,1:np-1)=x(1:nr,1:np-1)
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='SD') then
        nwork=2*n*(np-1)
        allocate(work(nwork+ndeb1),stat=istat)
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        if(istat/=0) stop 'ERROR: failure allocating work.'
        parmin%sdsaturated=.false.
        parmin%converged=.false.
        parmin%sdminimum=.true.
        icall=0
        do it=1,parmin%maxforcecall
            call nebforce(n,np,x,outs,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
            call calmaxforcecomponentanchors(atoms,np,outs,fnrm,fspmax)
            do ip = 1, np - 1, 1
               call vcopy(n, outs(ip)%fxyz(1,1), 1, fa(1,ip), 1)
            end do
            call checkconvergence(parmin,fspmax)
            !call sdminimum(atoms,iproc,n,np,nr*(np-1),xa,fa,fnrmtot,parmin,nwork,work)
            call sdminimum(iproc,nr*(np-1),nr*(np-1),xa,fa,fnrmtot,parmin,nwork,work)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            !x(1:nr,1:np-1)=xa(1:nr,1:np-1)
            !n(c) pold=pnow
            if(parmin%iflag<0 .or. parmin%converged) exit
            icall=icall+1
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='SDDIIS') then
        nwork=(3*parmin%idsx+3)*nr*(np-1) !2*n+nr
        allocate(work(nwork+ndeb1),stat=istat)
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        if(istat/=0) stop 'ERROR: failure allocating work.'
        parmin%sdsaturated=.false.
        parmin%converged=.false.
        parmin%sdminimum=.true.
        parmin%diisminimum=.false.
        icall=0
        do it=1,parmin%maxforcecall
            call nebforce(n,np,x,outs,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
            call calmaxforcecomponentanchors(atoms,np,outs,fnrm,fspmax)
            do ip = 1, np - 1, 1
               call vcopy(n, outs(ip)%fxyz(1,1), 1, fa(1,ip), 1)
            end do
            call checkconvergence(parmin,fspmax)
            if(it==4) then
                parmin%sdsaturated=.true.
                parmin%sdminimum=.false.
                parmin%diisminimum=.true.
                parmin%iflag=0
            endif
            if(iproc==0) write(*,*) 'REZA ',parmin%sdminimum,it,parmin%sdsaturated
            if(parmin%sdminimum) then
                !call sdminimum(atoms,iproc,n,np,nr*(np-1),xa,fa,fnrmtot,parmin,nwork,work)
            call sdminimum(iproc,nr*(np-1),nr*(np-1),xa,fa,fnrmtot,parmin,nwork,work)
                if(iproc==0 .and. parmin%converged) then
                    write(pnow%ifile,'(a)') 'converged before starting DIIS'
                    write(*,'(a)') 'converged before starting DIIS'
                endif
                if(parmin%itsd>parmin%nitsd .and. .not. parmin%sdsaturated) then
                if(iproc==0) then
                    write(pnow%ifile,'(a)') 'SD did not saturate, so diisminimum can not continue.'
                    write(*,'(a)') 'SD did not saturate, so diisminimum can not continue.'
                endif
                endif 
            !elseif(.not. (parmin%iflag==0 .and.  parmin%converged)) then
            elseif(parmin%diisminimum) then
                !call diisminimum(iproc,nr*(np-1),xa,fnrmtot,fa,parmin,nwork,work)
                call diisminimum(iproc,nr*(np-1),nr*(np-1),xa,fnrmtot,fa,parmin,nwork,work)
            endif
            if(parmin%sdsaturated .or. (parmin%iflag==0 .and. parmin%converged)) then
                parmin%sdminimum=.false.
                parmin%diisminimum=.true.
            endif
            if(parmin%diisdivergence) then
                parmin%diisdivergence=.false.
                parmin%sdsaturated=.false.
                parmin%sdminimum=.true.
                parmin%diisminimum=.false.
            endif 
            if(parmin%iflag==0 .and. parmin%converged) then
                parmin%sdminimum=.false.
                parmin%diisminimum=.false.
            endif
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            !x(1:nr,1:np-1)=xa(1:nr,1:np-1)
            !n(c) pold=pnow
            icall=icall+1
            if(parmin%iflag<0 .or. (parmin%iflag==0 .and.  parmin%converged)) then
                parmin%alpha=-1.d0
                exit
            endif
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='BFGS') then
        !stop 'FIX xold in call  calvmaxanchorforces'
        parmin%iflag=0
        !allocate(xold(n,0:np+ndeb2),stat=istat)
        !if(istat/=0) stop 'ERROR: failure allocating xold.'
        !call dmemocc(n*(np+1),n*(np+1+ndeb2),xold,'xold')
        nra=nr*(np-1)
        nwork=nra*nra+3*nra+3*nra*nra+3*nra
        allocate(work(nwork+ndeb1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        icall=0
        do it=1,parmin%maxforcecall
            !call calvmaxanchorforces(icall,n,np,x,xold,fends,etmax,f,xtmax,pnow,pold,fatsp, &
            !    nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            call nebforce(n,np,x,outs,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
            !call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            !call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            !call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,outs(ip)%fxyz,nr,fa(1,ip))
            enddo
            call calmaxforcecomponentanchors(atoms,np,outs,fnrm,fspmax)
            call checkconvergence(parmin,fspmax)
            call bfgs_splsad(iproc,nr*(np-1),xa,fnrmtot,fa,nwork,work,parmin)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            if(parmin%converged) exit
            if(parmin%iflag<=0) exit
            !n(c) pold=pnow
            icall=icall+1
            if(icall>1000) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        !deallocate(xold,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xold.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='FIRE') then
        parmin%iflag=0
        allocate(work(3*nr*(np-1)+ndeb1),stat=istat)
        call dmemocc(3*nr*(np-1),3*nr*(np-1)+ndeb1,work,'work')
        !allocate(xold(n,0:np+ndeb2),stat=istat)
        !call dmemocc(n*(np+1),n*(np+1+ndeb2),xold,'xold')
        parmin%dt=0.02d0
        !parmin%dt=0.01d0 !non-BigDFT
        if(istat/=0) stop 'ERROR: failure allocating work.'
        icall=0
        do it=1,parmin%maxforcecall
            !call calenergyforces(iproc,n,x,f,epot)
            !call calvmaxanchorforces(icall,n,np,x,xold,fends,etmax,f,xtmax,pnow,pold,fatsp, &
            !    parmin,nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            call nebforce(n,np,x,outs,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
            !fa(1:nr,1:np-1)=f(1:nr,1:np-1)
            !!if(iproc==0 .and. it==1) then
            !if(iproc==0) then
            !do ip=1,np-1
            !do i=1,n
            !    write(8000+it-1,'(2i4,1es20.10)') ip,i,f(i,ip)
            !enddo
            !enddo
            !close(8000+it-1)
            !endif
            !if(it<=10) then
            !call perpendicularforce(n,np,x,f,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
            !endif
            !!if(iproc==0 .and. it==1) then
            !if(iproc==0) then
            !do ip=1,np-1
            !do i=1,n
            !    write(9000+it-1,'(2i4,1es20.10)') ip,i,f(i,ip)
            !enddo
            !enddo
            !close(9000+it-1)
            !endif
            !xold(1:n,0:np)=x(1:n,0:np)
            !call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            !call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            !call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,outs(ip)%fxyz,nr,fa(1,ip))
            enddo
            call calmaxforcecomponentanchors(atoms,np,outs,fnrm,fspmax)
            call checkconvergence(parmin,fspmax)
            call fire_splsad(iproc,nr*(np-1),xa,fnrmtot,fa,work,parmin)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            if(parmin%converged) exit
            if(parmin%iflag<=0) exit
            !n(c) pold=pnow
            icall=icall+1
            if(icall>1000) exit
            !if(parmin%iflag<=0) exit
            !if(parmin%iflag<0 .or. parmin%converged) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        !deallocate(xold,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xold.'
    endif
    !-------------------------------------------------------------------------------------
    deallocate(xa,fa,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xa or fa.'
    do ip=1,np-1
       call deallocate_global_output(outs(ip))
    enddo
    if(iproc==0) then
    write(pnow%ifile,'(a,1x,a)') 'end of minimization_sp using ',parmin%approach
    write(*         ,'(a,1x,a)') 'end of minimization_sp using ',parmin%approach
    endif
end subroutine neb


subroutine atomic_copymoving_forward(atoms,n,x,nr,xa)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,nr,i,iat,ixyz,ir
    real(kind=8)::x(n),xa(nr)
    logical::move_this_coordinate
    ir=0
    do i=1,3*atoms%astruct%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            ir=ir+1
            xa(ir)=x(i)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
end subroutine atomic_copymoving_forward


subroutine atomic_copymoving_backward(atoms,nr,xa,n,x)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,nr,i,iat,ixyz,ir
    real(kind=8)::x(n),xa(nr)
    logical::move_this_coordinate
    ir=0
    do i=1,3*atoms%astruct%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            ir=ir+1
            x(i)=xa(ir)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
end subroutine atomic_copymoving_backward


subroutine calmaxforcecomponentsub(atoms,f,fnrm,fspmax)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::i,iat,ixyz
    real(kind=8)::f(3*atoms%astruct%nat),fnrm,fspmax
    logical::move_this_coordinate
    fspmax=0.d0
    fnrm=0.d0
    do i=1,3*atoms%astruct%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            fnrm=fnrm+f(i)**2
            fspmax=max(fspmax,abs(f(i)))
        endif
    enddo
    fnrm=sqrt(fnrm)
end subroutine calmaxforcecomponentsub


subroutine calmaxforcecomponentanchors(atoms,np,outs,fnrm,fspmax)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::np,i,ip,iat,ixyz
    real(kind=8)::fnrm,fspmax
    type(DFT_global_output), dimension(1:np-1), intent(in) :: outs
    logical::move_this_coordinate
    fspmax=0.d0
    fnrm=0.d0
    do i=1,3*atoms%astruct%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            do ip=1,np-1
                fnrm=fnrm+outs(ip)%fxyz(ixyz,iat)**2
                fspmax=max(fspmax,abs(outs(ip)%fxyz(ixyz,iat)))
            enddo
        endif
    enddo
    fnrm=sqrt(fnrm)
end subroutine calmaxforcecomponentanchors


!subroutine testwrite(nat,np,ratall)
!    implicit none
!    integer::nat,iat,np,ip
!    real(kind=8)::ratall(3,nat,0:np)
!    do ip=0,np
!    write(99,*) nat
!    write(99,*) 
!    do iat=1,nat
!        write(99,*) ratall(1:3,iat,ip)
!    enddo
!    enddo
!    close(99)
!end subroutine testwrite


subroutine nebforce(n,np,x,outs,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    integer::n,np,i,ip,istat,infocode
    type(DFT_global_output), dimension(1:np-1), intent(inout) :: outs
    real(kind=8)::x(n,0:np)
    real(kind=8)::tt,t1,t2,springcons,fnrmtot,time1,time2,fnrmarr(99),fspmaxarr(99)!,DNRM2
    real(kind=8), allocatable::tang(:,:)
    type(parametersplinedsaddle)::pnow
    logical::move_this_coordinate
    integer::iat,ixyz,mp
    integer, parameter::ndeb1=0,ndeb2=0
    type(run_objects) :: runObj

    allocate(tang(n,0:np+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating tang.'
    call dmemocc(n*(np+1),n*(np+1+ndeb2),tang,'tang')
    call run_objects_nullify(runObj)
    call run_objects_associate(runObj, ll_inputs, atoms, rst)
    call dmemocc(3*atoms%astruct%nat,3*(atoms%astruct%nat+ndeb1),runObj%atoms%astruct%rxyz,'runObj%atoms%astruct%rxyz')
    do ip=1,np-1
        call vcopy(n, x(1,ip), 1, runObj%atoms%astruct%rxyz(1,1), 1)
        call cpu_time(time1)
        call call_bigdft(runObj,outs(ip),nproc,iproc,infocode)
        pnow%ex(ip) = outs(ip)%energy
        call cpu_time(time2)
        ncount_bigdft=ncount_bigdft+1
        pnow%ncount_ll=pnow%ncount_ll+1
        pnow%time_ll=pnow%time_ll+(time2-time1)
        call calmaxforcecomponentanchors(atoms,2,outs(ip:ip),fnrmarr(ip),fspmaxarr(ip))
        !fnrmarr(ip)=DNRM2(n,f(1,ip),1) !HERE
    enddo
    call run_objects_free_container(runObj)
    call caltangentupwind(n,np,x,pnow%ex,tang)
    springcons=5.d-2
    !springcons=1.d-1 !non-BigDFT
    if(trim(pnow%docineb)=='yes') then
        mp=0
        pnow%epotci=pnow%ex(0)
        do ip=1,np
            if(pnow%ex(ip)>pnow%epotci) then
                mp=ip
                pnow%epotci=pnow%ex(ip)
            endif
        enddo
        if(iproc==0) write(*,'(a,i5,i3,es24.15,2es15.5)') 'mp,epotci ', &
            ncount_bigdft,mp,pnow%epotci,fnrmarr(mp),fspmaxarr(mp)
        if(mp==0 .or. mp==np) then
            !if(iproc==0) write(*,*) 'mp,exmax ',mp,pnow%epotci
            stop 'ERROR: highest energy image in cineb is one of the two ends.'
        endif
    endif
    fnrmtot=0.d0
    do ip=1,np-1
        call atomic_dot(atoms,outs(ip)%fxyz(1,1),tang(1,ip),tt)
        tt=-tt
        t1=0.d0;t2=0.d0
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                t1=t1+(x(i,ip+1)-x(i,ip))**2
                t2=t2+(x(i,ip)-x(i,ip-1))**2
            endif
        enddo
        t1=sqrt(t1);t2=sqrt(t2)
        if(trim(pnow%docineb)=='yes'.and. ip==mp) then
            tt=2.d0*tt
        else
            tt=tt+springcons*(t1-t2)
        endif
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
               outs(ip)%fxyz(ixyz, iat) = outs(ip)%fxyz(ixyz, iat) + tt * tang(i,ip)
            endif
        enddo
        call atomic_dot(atoms,outs(ip)%fxyz(1,1),outs(ip)%fxyz(1,1),tt);tt=sqrt(tt)
        fnrmtot=fnrmtot+tt
    enddo
    deallocate(tang,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating tang.'
end subroutine nebforce


subroutine splinedsaddle(n,nr,np,x,etmax,f,xtmax,parmin,outends,pnow,nproc, &
    iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft,fatsp)
    use module_base
    use module_interfaces
    use module_types
    use minimization_sp, only:parameterminimization_sp
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: inputs,ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    type(DFT_global_output), dimension(2), intent(in) :: outends
    integer::n,nr,np,i,ip,icall,istat,it,nwork,nra
    real(kind=8)::x(n,0:np),f(n,0:np),etmax,xtmax(n),fatsp(n)
    real(kind=8)::fspmax,fspnrm,barrier1,barrier2
    real(kind=8), allocatable::work(:)
    real(kind=8), allocatable::xa(:,:),fa(:,:),xold(:,:) !,fsp(:)
    type(parameterminimization_sp)::parmin
    type(parametersplinedsaddle)::pnow,pold
    integer, parameter::ndeb1=0,ndeb2=0
    parmin%converged=.false.
    if(iproc==0) then
        write(pnow%ifile,'(a,1x,a)') 'begin of minimization_sp using ',parmin%approach
        write(*,'(a,1x,a)') 'begin of minimization_sp using ',parmin%approach
    endif
    if(parmin%approach=='unknown') then
        if(iproc==0) then
            write(pnow%ifile,*) 'ERROR: minimization method not specified.'
            write(*,*) 'ERROR: minimization method not specified.'
        endif
        return
    endif
    !-------------------------------------------------------------------------------------
    allocate(xa(nr,np-1+ndeb2),fa(nr,np-1+ndeb2),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating xa or fa.'
    call dmemocc(nr*(np-1),nr*(np-1+ndeb2),xa,'xa')
    call dmemocc(nr*(np-1),nr*(np-1+ndeb2),fa,'fa')
    !-------------------------------------------------------------------------------------
    do ip=1,np-1
        call atomic_copymoving_forward(atoms,n,x(1,ip),nr,xa(1,ip))
    enddo
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='SD') then
        !stop 'FIX xold in call  calvmaxanchorforces'
        nwork=2*n*(np-1)
        allocate(work(nwork+ndeb1),stat=istat)
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        allocate(xold(n,0:np+ndeb2),stat=istat)
        call dmemocc(n*(np+1),n*(np+1+ndeb2),xold,'xold')
        if(istat/=0) stop 'ERROR: failure allocating xold.'
        parmin%sdsaturated=.false.
        parmin%sdminimum=.true.
        parmin%iflag=0
        icall=0
        do it=1,parmin%maxforcecall
            call calvmaxanchorforces(icall,n,np,x,xold,outends,etmax,f,xtmax,pnow,pold,fatsp, &
                nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
            enddo
            !fa(1:nr,1:np-1)=f(1:nr,1:np-1)
            !call sdminimum(atoms,iproc,n,np,nr*(np-1),xa,fa,etmax,parmin,nwork,work)
            call sdminimum(iproc,nr*(np-1),nr*(np-1),xa,fa,etmax,parmin,nwork,work)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            !x(1:nr,1:np-1)=xa(1:nr,1:np-1)
            pold=pnow
            if(parmin%iflag<0 .or. parmin%converged) exit
            icall=icall+1
        enddo
        deallocate(xold,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xold.'
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif !end of if statement for approach=='SD'
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='SDDIIS') then
        stop 'FIX xold in call  calvmaxanchorforces'
        nwork=(3*parmin%idsx+3)*nr*(np-1) !2*n+nr
        allocate(work(nwork+ndeb1),stat=istat)
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        if(istat/=0) stop 'ERROR: failure allocating work.'
        parmin%sdsaturated=.false.
        parmin%sdminimum=.true.
        parmin%diisminimum=.false.
        parmin%diisdivergence=.false.
        parmin%iflag=0
        icall=0
        do it=1,parmin%maxforcecall
            call calvmaxanchorforces(icall,n,np,x,xold,outends,etmax,f,xtmax,pnow,pold, &
                fatsp,nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
            enddo
            call testparmin(iproc,it,parmin,'before-sdminimum')
            if(parmin%sdminimum) then
                call sdminimum(iproc,nr*(np-1),nr*(np-1),xa,fa,etmax,parmin,nwork,work)
                if(iproc==0 .and. parmin%converged) then
                    write(pnow%ifile,*) 'converged before starting DIIS'
                    write(*,*) 'converged before starting DIIS'
                endif
                if(parmin%itsd>parmin%nitsd .and. .not. parmin%sdsaturated) then
                if(iproc==0) then
                    write(pnow%ifile,'(a)') 'SD did not saturate, so diisminimum can not continue.'
                    write(*,'(a)') 'SD did not saturate, so diisminimum can not continue.'
                endif
                endif 
            endif
            call testparmin(iproc,it,parmin,'after-sdminimum')
            if(.not. parmin%sdminimum .and. .not. parmin%converged) parmin%diisminimum=.true.
            call testparmin(iproc,it,parmin,'before-diisminimum')
            if(parmin%diisminimum) then
                call diisminimum(iproc,nr*(np-1),nr*(np-1),xa,etmax,fa,parmin,nwork,work)
            endif
            call testparmin(iproc,it,parmin,'after-diisminimum-1')
            !if(parmin%sdsaturated .or. (parmin%iflag==0 .and. parmin%converged))then
            !    parmin%sdminimum=.false.
            !    parmin%diisminimum=.true.
            !endif
            if(parmin%diisdivergence)then
                parmin%diisdivergence=.false.
                parmin%sdsaturated=.false.
                parmin%diisminimum=.false.
                parmin%sdminimum=.true.
            endif 
            call testparmin(iproc,it,parmin,'after-diisminimum-2')
            if(parmin%iflag==0 .or. parmin%converged) then
                parmin%sdminimum=.false.
                parmin%diisminimum=.false.
            endif
            call testparmin(iproc,it,parmin,'after-diisminimum-3')
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            pold=pnow
            icall=icall+1
            if(parmin%iflag<0 .or. (parmin%iflag==0 .and.  parmin%converged)) then
                parmin%alpha=-1.d0
                exit
            endif
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='BFGS') then
        !stop 'FIX xold in call  calvmaxanchorforces'
        parmin%iflag=0
        allocate(xold(n,0:np+ndeb2),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating xold.'
        call dmemocc(n*(np+1),n*(np+1+ndeb2),xold,'xold')
        nra=nr*(np-1)
        nwork=nra*nra+3*nra+3*nra*nra+3*nra
        allocate(work(nwork+ndeb1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        icall=0
        do it=1,parmin%maxforcecall
            call calvmaxanchorforces(icall,n,np,x,xold,outends,etmax,f,xtmax,pnow,pold,fatsp, &
                nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
            enddo
            call bfgs_splsad(iproc,nr*(np-1),xa,etmax,fa,nwork,work,parmin)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            if(parmin%converged) exit
            if(parmin%iflag<=0) exit
            pold=pnow
            icall=icall+1
            if(icall>1000) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        deallocate(xold,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xold.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='DFP') then
        !stop 'FIX xold in call  calvmaxanchorforces'
        parmin%iflag=0
        allocate(xold(n,0:np+ndeb2),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating xold.'
        call dmemocc(n*(np+1),n*(np+1+ndeb2),xold,'xold')
        nra=nr*(np-1)
        nwork=nra*nra+3*nra+3*nra*nra+2*nra
        allocate(work(nwork+ndeb1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        call dmemocc(nwork,nwork+ndeb1,work,'work')
        icall=0
        do it=1,parmin%maxforcecall
            call calvmaxanchorforces(icall,n,np,x,xold,outends,etmax,f,xtmax,pnow,pold,fatsp, &
                nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
            enddo
            call dfp_splsad(iproc,nr*(np-1),xa,etmax,fa,nwork,work,parmin)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            if(parmin%converged) exit
            if(parmin%iflag<=0) exit
            pold=pnow
            icall=icall+1
            if(icall>1000) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        deallocate(xold,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xold.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(parmin%approach)=='FIRE') then
        parmin%iflag=0
        allocate(work(3*nr*(np-1)+ndeb1),stat=istat)
        call dmemocc(3*nr*(np-1),3*nr*(np-1)+ndeb1,work,'work')
        allocate(xold(n,0:np+ndeb2),stat=istat)
        call dmemocc(n*(np+1),n*(np+1+ndeb2),xold,'xold')
        if(istat/=0) stop 'ERROR: failure allocating work.'
        icall=0
        do it=1,parmin%maxforcecall
            !call calenergyforces(iproc,n,x,f,epot)
            call calvmaxanchorforces(icall,n,np,x,xold,outends,etmax,f,xtmax,pnow,pold,fatsp, &
                nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
            !if(iproc==0 .and. it==1) then
            if(iproc==0) then
            do ip=1,np-1
            do i=1,n
                write(8000+it-1,'(2i4,1es20.10)') ip,i,f(i,ip)
            enddo
            enddo
            close(8000+it-1)
            endif
            if(it<=-1) then
            call perpendicularforce(n,np,x,f,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
            endif
            !if(iproc==0 .and. it==1) then
            if(iproc==0) then
            do ip=1,np-1
            do i=1,n
                write(9000+it-1,'(2i4,1es20.10)') ip,i,f(i,ip)
            enddo
            enddo
            close(9000+it-1)
            endif
            !xold(1:n,0:np)=x(1:n,0:np)
            call calmaxforcecomponentsub(atoms,fatsp,fspnrm,fspmax)
            call reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
            call checkconvergence(parmin,fspmax)
            do ip=1,np-1
                call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
            enddo
            call fire_splsad(iproc,nr*(np-1),xa,etmax,fa,work,parmin)
            do ip=1,np-1
                call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
            enddo
            if(parmin%converged) exit
            if(parmin%iflag<=0) exit
            pold=pnow
            icall=icall+1
            if(icall>1000) exit
            !if(parmin%iflag<=0) exit
            !if(parmin%iflag<0 .or. parmin%converged) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        deallocate(xold,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xold.'
    endif
    !-------------------------------------------------------------------------------------
    if(iproc==0) then
        barrier1=27.2113845d0*(etmax-pnow%exends_b(1))
        barrier2=27.2113845d0*(etmax-pnow%exends_b(2))
        write(pnow%ifile,'(a,2f15.5)') 'barrier heights in eV',barrier1,barrier2
        write(*         ,'(a,2f15.5)') 'barrier heights in eV',barrier1,barrier2
    endif
    deallocate(xa,fa,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xa or fa.'
    if(iproc==0) then
        write(pnow%ifile,'(a,1x,a)') 'end of minimization_sp using ',parmin%approach
        write(*         ,'(a,1x,a)') 'end of minimization_sp using ',parmin%approach
    endif
end subroutine splinedsaddle


subroutine bfgs_splsad(iproc,nr,x,epot,f,nwork,work,parmin)
    !use minimization, only:parameterminimization
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::iproc,nr,nwork,mf,my,ms,nrsqtwo,iw1,iw2,iw3,iw4,info,i,j,l,mx
    real(kind=8)::x(nr),f(nr),epot,work(nwork)
    !real(kind=8), allocatable::eval(:),umat(:)
    type(parameterminimization_sp)::parmin
    real(kind=8)::DDOT,tt1,tt2,de,fnrm,calnorm,fmax,calmaxforcecomponent,beta
    real(kind=8)::tt3,tt4,tt5,tt6
    real(kind=8), save::epotold,alpha,alphamax,zeta
    logical, save::reset
    integer, save::isatur
    integer:: ierr
    include 'mpif.h'
    if(nwork/=nr*nr+3*nr+3*nr*nr+3*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    nrsqtwo=nr*nr*2
    mf=nr*nr+1       !for force of previous iteration in wiki notation
    my=mf+nr         !for y_k in wiki notation
    ms=my+nr         !for s_k in wiki notation
    iw1=ms+nr        !work array to keep the hessian untouched
    iw2=iw1+nr*nr    !for work array of DSYTRF
    iw3=iw2+nrsqtwo  !for p_k in wiki notation
    mx =iw3+nr       !for position of previous iteration
    iw4=mx+nr        !for eigenvalues of inverse og hessian
    if(parmin%iflag==0) then
        parmin%iflag=1
        parmin%iter=0
        epotold=epot
        alpha=1.d-1
        reset=.false.
        alphamax=0.9d0
        zeta=1.d0
        isatur=0
    else
        parmin%iter=parmin%iter+1
    endif
    de=epot-epotold
    fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
    if(iproc==0) then
    write(*,'(a10,i4,es23.15,es11.3,2es12.5,1es12.4)') &
        'BFGSMIN   ',parmin%iter,epot,de,fnrm,fmax,zeta
    endif
    !if(parmin%iter==602) then
    !    do i=1,nr/3
    !        write(31,*) x(i*3-2),x(i*3-1),x(i*3-0)
    !    enddo
    !    stop
    !endif
    if(fmax<parmin%fmaxtol) then
        parmin%converged=.true.
        parmin%iflag=0
        if(iproc==0) then
        write(*,'(a,i4,es23.15,2es12.5)') &
            'BFGS FINISHED: itfire,epot,fnrm,fmax ',parmin%iter,epot,fnrm,fmax
        endif
        return
    endif

    !if(de>0.d0 .and. zeta>1.d-1) then
    if(de>5.d-2) then
        epot=epotold
        x(1:nr)=work(mx:mx-1+nr)
        f(1:nr)=work(mf:mf-1+nr)
        reset=.true.
        !alpha=max(alpha*0.5d0/1.1d0,1.d-2)
        zeta=max(zeta*1.d-1,1.d-5)
        isatur=0
    else
        !zeta=1.d0
        !if(zeta>1.d-1) zeta=min(zeta*1.1d0,1.d0)
        zeta=min(zeta*1.1d0,1.d0)
        isatur=isatur+1
    endif
    if(parmin%iter==0 .or. reset) then
        !reset=.false.
        if(isatur>=10) then
            reset=.false.
            !alpha=5.d-1
        endif
        work(1:nr*nr)=0.d0
        do i=1,nr
            work(i+(i-1)*nr)=zeta*parmin%alphax
        enddo
        work(iw3:iw3-1+nr)=zeta*parmin%alphax*f(1:nr)
    else
        work(ms:ms-1+nr)=x(1:nr)-work(mx:mx-1+nr)
        work(my:my-1+nr)=work(mf:mf-1+nr)-f(1:nr)
        tt1=DDOT(nr,work(my),1,work(ms),1)
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(my-1+j)
            enddo
            work(iw2-1+i)=tt2
        enddo
        tt2=DDOT(nr,work(my),1,work(iw2),1)
        !write(21,*) parmin%iter,tt1,tt2
        !tt1=max(tt1,1.d-2)
        do i=1,nr
            do j=i,nr
                l=i+(j-1)*nr
                work(l)=work(l)+(tt1+tt2)*work(ms-1+i)*work(ms-1+j)/tt1**2- &
                    (work(iw2-1+i)*work(ms-1+j)+work(iw2-1+j)*work(ms-1+i))/tt1
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        !do i=1,nr
        !    tt2=0.d0
        !    do j=1,nr
        !        tt2=tt2+work(j+(i-1)*nr)*f(j)
        !    enddo
        !    work(iw3-1+i)=tt2
        !enddo
        !write(31,*) zeta
        work(iw1:iw1-1+nr*nr)=work(1:nr*nr)
        if(iproc==0) then
        call DSYEV('V','L',nr,work(iw1),nr,work(iw4),work(iw2),nrsqtwo,info)
        if(info/=0) stop 'DSYEV'
        endif
        call MPI_BCAST(work(iw1),nr*nr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(work(iw4),nr   ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        tt1=work(iw4+0)    ; tt2=work(iw4+1)    ; tt3=work(iw4+2)
        tt4=work(iw4+nr-3) ; tt5=work(iw4+nr-2) ; tt6=work(iw4+nr-1)
        if(iproc==0) then
        write(41,'(i5,6es15.5)') parmin%iter,tt1,tt2,tt3,tt4,tt5,tt6
        endif
        work(iw3:iw3-1+nr)=0.d0
        if(parmin%iter<3) then
            beta=5.d0/parmin%alphax
        elseif(parmin%iter<6) then
            beta=1.d0/parmin%alphax
        else
            beta=1.d-2/parmin%alphax
        endif
        do j=1,nr
            tt1=DDOT(nr,work(iw1+nr*(j-1)),1,f,1)
            tt2=1.d0/sqrt(1.d0/work(iw4-1+j)**2+beta**2)
            do i=1,nr
                work(iw3-1+i)=work(iw3-1+i)+tt1*work(iw1-1+i+nr*(j-1))*tt2
            enddo
        enddo
    endif
    epotold=epot
    work(mf:mf-1+nr)=f(1:nr)
    work(mx:mx-1+nr)=x(1:nr)
    !if(isatur>5) alpha=min(alphamax,alpha*1.1d0)
    alpha=min(alphamax,alpha*1.1d0)
    x(1:nr)=x(1:nr)+alpha*work(iw3:iw3-1+nr)
end subroutine bfgs_splsad


subroutine dfp_splsad(iproc,nr,x,epot,f,nwork,work,parmin)
    !use minimization, only:parameterminimization
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::iproc,nr,nwork,mf,my,ms,nrsqtwo,iw1,iw2,iw3,info,i,j,l,mx
    real(kind=8)::x(nr),f(nr),epot,work(nwork)
    integer, allocatable::ipiv(:)
    !type(parameterminimization)::parmin
    type(parameterminimization_sp)::parmin
    real(kind=8)::DDOT,tt1,tt2,de,fnrm,calnorm,fmax,calmaxforcecomponent,dx
    real(kind=8), save::epotold,alpha,alphamax,zeta,zetaold
    logical, save::reset
    if(nwork/=nr*nr+3*nr+3*nr*nr+2*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    mf=nr*nr+1 !for force of previous iteration in wiki notation
    my=nr*nr+nr+1 !for y_k in wiki notation
    ms=nr*nr+2*nr+1 !for s_k in wiki notation
    iw1=nr*nr+3*nr+1 !work array to keep the hessian untouched
    iw2=nr*nr+3*nr+nr*nr+1 !for work array of DSYTRF
    iw3=nr*nr+3*nr+3*nr*nr+1 !for p_k in wiki notation
    mx= nr*nr+3*nr+3*nr*nr+nr+1 !for position of previous iteration
    nrsqtwo=nr*nr*2
    if(parmin%iflag==0) then
        parmin%iflag=1
        parmin%iter=0
        epotold=epot
        alpha=7.d-1
        reset=.false.
        alphamax=1.d0
        zeta=1.d0/parmin%alphax
        !zeta=1000.d0
        zetaold=0.d0
    else
        parmin%iter=parmin%iter+1
    endif
    de=epot-epotold
    fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
    if(iproc==0) then
    write(*,'(a10,i4,es23.15,es11.3,2es12.5,1es12.4)') &
        'DFPMIN    ',parmin%iter,epot,de,fnrm,fmax,alpha
    endif
    !if(parmin%iter==1714) then
    !    do i=1,nr/3
    !        write(31,*) x(i*3-2),x(i*3-1),x(i*3-0)
    !    enddo
    !    stop
    !endif
    if(fmax<parmin%fmaxtol) then
        parmin%converged=.true.
        parmin%iflag=0
        if(iproc==0) then
        write(*,'(a,i4,es23.15,2es12.5)') &
            'DFP FINISHED: itfire,epot,fnrm,fmax ',parmin%iter,epot,fnrm,fmax
        return
        endif
    endif

    if(de>5.d-3) then
    !if(de>5.d1) then !non-BigDFT
    !if(de>1.d0) then !non-BigDFT
        epot=epotold
        x(1:nr)=work(mx:mx-1+nr)
        f(1:nr)=work(mf:mf-1+nr)
        reset=.true.
        !alpha=max(alpha*0.5d0/1.1d0,1.d-2)
        !alpha=alpha*0.5d0/1.1d0
        alpha=alpha*0.99d0/1.1d0
        if(alpha<1.d-2) then
            if(iproc==0) then
            write(*,'(a)') 'ERROR: it is unreasonable to continue minimization '
            write(*,'(a)') '       since the stepsize is too small. This can   '
            write(*,'(a)') '       occur if the gradient is wrong or if the    '
            write(*,'(a)') '       input standard stepsize is too large.       '
            endif
            stop
        endif
    endif
    if(de>1.d-2) then !non-BigDFT
        zeta=min(zeta*2.d0,1.d2/parmin%alphax)
    else
        zeta=max(zeta*0.9d0,1.d-2/parmin%alphax)
    endif
    if(parmin%iter<10 .or. reset) then
        reset=.false.
        work(1:nr*nr)=0.d0
        do i=1,nr
            work(i+(i-1)*nr)=2.d0/parmin%alphax
        enddo
        work(iw3:iw3-1+nr)=parmin%alphax*f(1:nr)/2.d0
    else
        work(ms:ms-1+nr)=x(1:nr)-work(mx:mx-1+nr)
        work(my:my-1+nr)=work(mf:mf-1+nr)-f(1:nr)
        !write(21,*) parmin%iter,DNRM2(nr,work(my),1)
        tt1=DDOT(nr,work(my),1,work(ms),1)
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(ms-1+j)
            enddo
            work(iw2-1+i)=tt2
        enddo
        tt2=DDOT(nr,work(ms),1,work(iw2),1)
        do i=1,nr
            do j=i,nr
                l=i+(j-1)*nr
                work(l)=work(l)+work(my-1+i)*work(my-1+j)/tt1-work(iw2-1+i)*work(iw2-1+j)/tt2
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        allocate(ipiv(nr))
        work(iw1:iw1-1+nr*nr)=work(1:nr*nr)
        do i=1,nr
            work(iw1+i+nr*(i-1))=work(iw1+i+nr*(i-1))+5.d2 !zeta
            !work(i+nr*(i-1))=work(i+nr*(i-1))-zetaold+zeta
        enddo
        zetaold=zeta
        write(31,'(i4,2f15.5)') parmin%iter,zeta,zetaold
        work(iw3:iw3-1+nr)=f(1:nr)
        !http://alcinoe.net/fortran/optim/optim.f90.html
        !http://www.netlib.no/netlib/lapack/double/dsytrf.f
        !http://www.netlib.no/netlib/lapack/double/dsytrs.f
        call DSYTRF('L',nr,work(iw1),nr,ipiv,work(iw2),nrsqtwo,info)
        if(info/=0) then;write(*,*) 'ERROR: DSYTRF failed: info',info;stop;endif
        call DSYTRS('L',nr,1,work(iw1),nr,ipiv,work(iw3),nr,info)
        if(info/=0) then;write(*,*) 'ERROR: DSYTRS failed: info',info;stop;endif
        deallocate(ipiv)
    endif
    epotold=epot
    work(mf:mf-1+nr)=f(1:nr)
    work(mx:mx-1+nr)=x(1:nr)
    alpha=min(alphamax,alpha*1.1d0)
    !x(1:nr)=x(1:nr)+alpha*work(iw3:iw3-1+nr)
    do i=1,nr
        dx=alpha*work(iw3-1+i)
        x(i)=x(i)+sign(min(abs(dx),5.d0),dx)
    enddo
end subroutine dfp_splsad


subroutine reportcalvmaxanchorforces(iproc,icall,n,np,x,etmax,fspnrm,fspmax,pnow,atoms,ncount_bigdft)
    use modulesplinedsaddle, only:parametersplinedsaddle
    use module_types
    implicit none
    integer::iproc,icall,n,np,ncount_bigdft
    real(kind=8)::x(n,0:np),etmax,fspnrm,fspmax,cbh1,cbh2
    type(parametersplinedsaddle)::pnow
    type(atoms_data), intent(inout) :: atoms
    character(25), parameter::frt1='(a,2i5,2es24.15,a,2f15.5)'
    character(len=20)::filename
    character(len=3)::fn
    if(iproc==0) then
        cbh1=27.2113845d0*(etmax-pnow%exends_b(1))
        cbh2=27.2113845d0*(etmax-pnow%exends_b(2))
        write(* ,frt1) 'SP ',ncount_bigdft,icall,fspnrm,fspmax,' fort52 ',cbh1,cbh2
        !write(52,frt1) 'SP ',ncount_bigdft,fspnrm,fspmax,' fort52 ',cbh1,cbh2
        write(fn,'(i3.3)') icall
        filename='anchorpoints'//fn//'.xyz' 
        call writeanchorpoints(n,np,x,filename,atoms)
    endif
end subroutine reportcalvmaxanchorforces


subroutine testparmin(iproc,it,parmin,str)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::iproc,it,ii
    type(parameterminimization_sp)::parmin
    character(*)::str
    logical::l2,l3,l4,l5,l6
    if(iproc==0) then
        ii=parmin%iflag
        l2=parmin%converged
        l3=parmin%sdminimum
        l4=parmin%sdsaturated
        l5=parmin%diisminimum
        l6=parmin%diisdivergence
        write(*,'(a,2i4,5l4,1x,a)') 'ALIREZA ',it,ii,l2,l3,l4,l5,l6,str
    endif
end subroutine testparmin


subroutine perpendicularforce(n,np,x,f,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    integer::n,np,i,ip,istat,infocode,mp
    real(kind=8)::x(n,0:np),f(n,0:np),epotarr(0:100)
    type(parametersplinedsaddle)::pnow
    real(kind=8)::tt,fnrm,fnrmmax,time1,time2
    real(kind=8), allocatable::tang(:,:)
    logical::move_this_coordinate
    integer::iat,ixyz
    integer, parameter::ndeb1=0,ndeb2=0
    type(run_objects) :: runObj
    type(DFT_global_output), dimension(0:np) :: outs

    allocate(tang(n,0:np+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating tang.'
    call dmemocc(n*(np+1),n*(np+1+ndeb2),tang,'tang')
    call run_objects_nullify(runObj)
    call run_objects_associate(runObj, ll_inputs, atoms, rst)
    mp=-1
    fnrmmax=0.d0
    do ip=1,np-1
        call atomic_dot(atoms,f(1,ip),f(1,ip),fnrm) ; fnrm=sqrt(fnrm)
        if(fnrm>fnrmmax) then
            fnrmmax=fnrm
            mp=ip
        endif
    enddo
    if(mp==-1) stop 'ERROR: in perpendicularforce mp==-1'
    !do ip=1,np-1
    do ip=mp,mp
       call init_global_output(outs(ip), atoms%astruct%nat)
        call vcopy(n, x(1,ip), 1, runObj%atoms%astruct%rxyz(1,1), 1)
        call cpu_time(time1)
        call call_bigdft(runObj,outs(ip),nproc,iproc,infocode)
        epotarr(ip) = outs(ip)%energy
        call cpu_time(time2)
        ncount_bigdft=ncount_bigdft+1
        pnow%ncount_ll=pnow%ncount_ll+1
        pnow%time_ll=pnow%time_ll+(time2-time1)
    enddo
    call run_objects_free_container(runObj)
    epotarr(0)=pnow%ex(0)
    epotarr(np)=pnow%ex(pnow%npv)
    call caltangentupwind(n,np,x,epotarr,tang)
    !do ip=1,np-1
    do ip=mp,mp
        call atomic_dot(atoms,outs(ip)%fxyz(1,1),tang(1,ip),tt)
        tt=-tt
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                f(i,ip)=f(i,ip)+1.d-1*(outs(ip)%fxyz(ixyz,iat)+tt*tang(i,ip))
            endif
        enddo
        call deallocate_global_output(outs(ip))
    enddo
    deallocate(tang,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating tang.'
end subroutine perpendicularforce


subroutine calvmaxanchorforces(istep,n,np,x,xold,outends,etmax,f,xtmax,pnow,pold,ftmax, &
    nproc,iproc,atoms,rst,inputs,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use minimization_sp, only:parameterminimization_sp
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: inputs,ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    type(DFT_global_output), dimension(2), intent(in) :: outends
    integer::n,np,mp,i,ip,j,infocode
    real(kind=8)::x(n,0:np),xold(n,0:np),f(n,0:np),xtmax(n),ftmax(n)
    integer::istat,istep
    type(parametersplinedsaddle)::pnow,pold
    !type(parameterminimization_sp)::parmin
    real(kind=8)::etmax,tt,time1,time2
    real(kind=8), allocatable::dd(:,:,:)
    integer, parameter::ndeb2=0 !n(c) ndeb1
    type(run_objects) :: runObj
    type(DFT_global_output) :: outs
    !----------------------------------------
    allocate(dd(n,n,np-1+ndeb2),stat=istat);if(istat/=0) stop 'ERROR: failure allocating dd.'
    if(istep==0) xold(1:n,0:np)=x(1:n,0:np)
    call dmemocc(n*n*(np-1),n*n*(np-1+ndeb2),dd,'dd')
    call equalarclengthparametrization(atoms,n,np,x,pnow%s,pnow%h)
    call factor_cubic(np,pnow%h,pnow%e1,pnow%e2)
    !call caltmax(n,np,x,etmax,xtmax,ftmax,pnow,pold)
    call caltmax2(istep,n,np,x,xold,outends,etmax,xtmax,ftmax,pnow,pold,nproc,iproc,&
        atoms,rst,ll_inputs,ncount_bigdft)
    if(trim(pnow%hybrid)=='yes') then
       call run_objects_nullify(runObj)
       call run_objects_associate(runObj, inputs, atoms, rst, xtmax(1))
       call init_global_output(outs, atoms%astruct%nat)
       inputs%inputPsiId=0
       call cpu_time(time1)
       call call_bigdft(runObj,outs,nproc,iproc,infocode)
       call cpu_time(time2)
       call run_objects_free_container(runObj)
       etmax = outs%energy
       call deallocate_global_output(outs, ftmax(1))
       ncount_bigdft=ncount_bigdft+1
       pnow%ncount=pnow%ncount+1
       pnow%time=pnow%time+(time2-time1)
    endif
    call calindex(np,pnow%s,pnow%tmax,mp,'calvmaxanchorforces')
    !write(*,'(a,i5,2e24.15)') 'SP ',istep,xtmax(1),xtmax(2)-1.d0
    !write(31,*) xtmax(1),xtmax(2)
    !call checkpathway(n,nr,np,x,npv,pnow,pold)
    !if(trim(parmin%approach)=='SDDIIS') then
    !    pnow%granot=.false.
    !else
    !    pnow%granot=.true.
    !endif
    !write(*,*) 'REZA ',ftmax(1),ftmax(2)
    if(pnow%granot) then
        call prepdd(atoms,n,np,x,pnow%e1,pnow%e2,pnow%h,pnow%s,mp,pnow%tmax,dd)
        do ip=1,np-1
            do i=1,n
                !dd(i,ip)=1.d0
                tt=0.d0
                do j=1,n
                    tt=tt+ftmax(j)*dd(j,i,ip)
                    !write(61,*) dd(j,i,ip)
                enddo
                f(i,ip)=tt
            enddo
        enddo
                !write(61,*) 
    else
        !write(*,*) 'granot',pnow%granot
        !call projectoutperpendicularforce(n,nr,np,x,f,pnow)
    endif
    !xold(1:n,0:np)=x(1:n,0:np)
    !stop
    deallocate(dd,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating dd.'
end subroutine calvmaxanchorforces


!subroutine projectoutperpendicularforce(n,nr,np,x,f,pnow)
!    use modulesplinedsaddle, only:parametersplinedsaddle
!    implicit none
!    integer::n,nr,np,npv,nflat,ip,i,istat,mp
!    real(kind=8)::x(n,0:np),f(n,0:np),t1,t2,dt,t,mydot,tt,epot
!    type(parametersplinedsaddle)::pnow
!    real(kind=8), allocatable::tang(:)
!    real(kind=8), allocatable::ft(:)
!    allocate(ft(n),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ft.'
!    allocate(tang(n),stat=istat);if(istat/=0) stop 'ERROR: failure allocating tang.'
!    do ip=1,np-1
!        t=pnow%s(ip)
!        mp=ip
!        do i=1,nr
!            pnow%y(0:np)=x(i,0:np)
!            call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
!            call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),t,pnow%c,t1,tang(i),t2)
!        enddo
!        do i=nr+1,n
!            tang(i)=0.d0
!        enddo
!        call calenergyforces(n,x(1,ip),epot,ft)
!        call normalizevector2(nr,tang)
!        tt=mydot(nr,tang,ft)
!        f(1:nr,ip)=ft(1:nr) -1.d0*tt*tang(1:nr)
!        f(nr+1:n,ip)=0.d0
!    enddo
!    deallocate(tang,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating tang.'
!    deallocate(ft,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ft.'
!end subroutine projectoutperpendicularforce


subroutine checkpathway(iproc,istep,n,np,x,xold,pnow)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,n,np,ip,i
    real(kind=8)::x(n,0:np),xold(n,0:np),dmax
    type(parametersplinedsaddle)::pnow
    dmax=0.d0
    do ip=1,np-1
        do i=1,n
            dmax=max(dmax,abs(x(i,ip)-xold(i,ip)))
        enddo
    enddo
    if(iproc==0) write(*,'(a,i4,1es14.5)') 'fort53 ',istep,dmax
    !if(dmax>5.d-3) then
    !if(dmax>2.d-2) then
    if(dmax>5.d-3) then !non-BigDFT
        pnow%do_fill_ex_exd=.true.
    else
        pnow%do_fill_ex_exd=.false.
    endif
end subroutine checkpathway


!subroutine caltmax(n,np,x,etmax,xt,ft,pnow,pold)
!    use modulesplinedsaddle, only:parametersplinedsaddle
!    implicit none
!    integer::n,np,mparr(2),iepotmax
!    real(kind=8)::x(n,0:np),xt(n),ft(n)
!    type(parametersplinedsaddle)::pnow,pold
!    real(kind=8)::dt,epotmax,etmax,calnorm,zbrent,tl,tr,t1,t2
!    integer::istat,i,ip,mp,npv
!    real(kind=8), allocatable::xall(:,:)
!    npv=pnow%ns*np
!    allocate(xall(n,npv-1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xall.'
!    pnow%ex(0)=pnow%exends(1)
!    pnow%ex(npv)=pnow%exends(2)
!    pnow%sv(0)=0.d0
!    pnow%sv(npv)=pnow%s(np)
!    do i=1,n
!        pnow%y(0:np)=x(i,0:np)
!        call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
!        dt=pnow%s(np)/npv
!        do ip=1,npv-1
!            pnow%sv(ip)=dt*ip 
!            call calindex(np,pnow%s,pnow%sv(ip),mp)
!            call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),pnow%sv(ip),pnow%c,xall(i,ip),t1,t2)
!        enddo
!    enddo
!    epotmax=max(pnow%exends(1),pnow%exends(2)) !-1.d50
!    iepotmax=-1
!    do ip=1,npv-1
!        call calenergyforces(n,xall(1,ip),pnow%ex(ip),ft)
!        !write(*,'(a5,i6,3f24.15)') 'along',ip,pnow%ex(ip),xall(1,ip),xall(2,ip)
!        if(pnow%ex(ip)>epotmax) then
!            epotmax=pnow%ex(ip)
!            iepotmax=ip
!        endif
!    enddo
!    if(iepotmax==-1) stop 'ERROR: iepotmax=-1'
!    write(*,*) 'epotmax,iepotmax',epotmax,iepotmax,npv
!    tl=pnow%sv(iepotmax-1)
!    tr=pnow%sv(iepotmax+1)
!    !write(*,*) 'tts',tl,tr
!    call calindex(np,pnow%s,tl,mparr(1))
!    call calindex(np,pnow%s,tr,mparr(2))
!    !write(*,*) 'mparr(1),mparr(2)',mparr(1),mparr(2)
!    pnow%tmax=zbrent(tl,tr,1.d-5,n,np,x,pnow,mparr,xt,ft,etmax)
!    !write(*,*) 'brent tmax',pnow%tmax
!    deallocate(xall,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xall.'
!end subroutine caltmax


subroutine caltmax2(istep,n,np,x,xold,outends,epot,xt,ft,pnow,pold,nproc,iproc,atoms, &
        rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    type(DFT_global_output), dimension(2), intent(in) :: outends
    integer::istep,n,np,ip,mp,mpv,iter,npv,ipv,ibad,ibadold
    character(len=20)::filename
    character(len=3)::fn
    real(kind=8)::x(n,0:np),xold(n,0:np),xt(n),ft(n)
    type(parametersplinedsaddle)::pnow,pold
    real(kind=8)::epot,alpha,oneisideal,vdold,vdtol
    !n(c) character(32), parameter::frt =  '(3i5,e24.15,e15.6,e13.5,4e12.4)'
    character(37), parameter::frt2='(a,3i5,es24.15,es15.6,es13.5,4es12.4)'
    integer, save::iii=-1
    real(kind=8)::fnrm,vd,vdd,vdc,vddc,vq,vdq,vddq
    npv=np+pnow%ns2
    if(istep==0) npv=min(int(1.5d0*real((np+pnow%ns2),8)),30)
    pnow%granot=.true.
    !--------------------------------------------------------
    if(iproc==0) then
        write(*,'(a,i4)',advance='no') 'fort57 ',istep
        do ip=1,np
            write(*,'(1es10.2)',advance='no') pnow%h(ip)
        enddo
        write(*,'(a,1es13.5)') ' snp ',pnow%s(np)
    endif
    !--------------------------------------------------------
    if(iproc==0) then
        write(fn,'(i3.3)') istep
        filename='path'//fn//'.xyz' 
        call writepathway(n,np,x,filename,atoms)
    endif
    !--------------------------------------------------------
    !call epot_along_traj(istep,n,nr,np,x,npv,pnow,nproc,iproc,atoms,rst,inputs,ncount_bigdft)
    !--------------------------------------------------------
    call checkpathway(iproc,istep,n,np,x,xold,pnow)
    !pnow%do_fill_ex_exd=.true.
!1358 continue
    !if(istep<3 .or. pnow%do_fill_ex_exd .or. pold%npv>np+pnow%ns2) then
    if(istep<1 .or. pnow%do_fill_ex_exd) then
        call fill_ex_exd(istep,n,np,x,outends,npv,pnow,pold,xt,ft,nproc,iproc,atoms, &
            rst,ll_inputs,ncount_bigdft)
        !call guessinitialtmax_quintic(npv,pnow,iproc) !non-BigDFT
        !call guessinitialtmax_hermite(npv,pnow)
        call guessinitialtmax_cubic(npv,pnow) !non-BigDFT
        xold(1:n,0:np)=x(1:n,0:np)
    else
        npv=pold%npv
        pnow%ex(0:npv)=pold%ex(0:npv)
        pnow%exd(0:npv)=pold%exd(0:npv)
        pnow%sv(0)=pold%sv(0)/pold%s(np)*pnow%s(np)
        do ipv=1,npv
            pnow%sv(ipv)=pold%sv(ipv)/pold%s(np)*pnow%s(np)
            pnow%hv(ipv)=pnow%sv(ipv)-pnow%sv(ipv-1)
        enddo
        pnow%tmax=pold%tmax/pold%s(np)*pnow%s(np)
        !call guessinitialtmax_hermite(npv,pnow)
    endif
    iii=iii+1
    call write_v_of_t(iproc,istep,npv,pnow,iii,'vogt')
    if(iproc==0) then
        write(1000+iii,'(a,i5)') '#istep ',istep
        do ipv=0,npv
            write(1000+iii,'(3es24.15)') pnow%sv(ipv),pnow%ex(ipv),pnow%exd(ipv)
        enddo
        close(1000+iii)
    endif
    call calindex(np,pnow%s,pnow%tmax,mp,'caltmax2_1')
    call func(pnow%tmax,epot,vd,n,np,x,pnow,mp,xt,ft,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    call calindex(npv,pnow%sv,pnow%tmax,mpv,'caltmax2_2')
    call insertpoint(npv,epot,vd,mpv,np,pnow)
    !--------------------------------------------------------
    call write_v_of_t(iproc,istep,npv,pnow,iii,'nogt')
    if(iproc==0) then
        write(*,'(a,i4)',advance='no') 'fort54 ',istep
        do ipv=0,npv
            write(*,'(1es14.5)',advance='no') pnow%sv(ipv)
        enddo
        write(*,'(1es14.5)') pnow%tmax
        !-----------------------------------------
        write(*,'(a,i4)',advance='no') 'fort55 ',istep
        do ipv=0,npv
            write(*,'(1es14.5)',advance='no') pnow%ex(ipv)
        enddo
        write(*,'(1es14.5)') epot
        !-----------------------------------------
        write(2000+iii,'(a,i5)') '#istep ',istep
        do ipv=0,npv
            write(2000+iii,'(3es24.15)') pnow%sv(ipv),pnow%ex(ipv),pnow%exd(ipv)
        enddo
        close(2000+iii)
    endif
    !--------------------------------------------------------
    !call polish_sv(npv,pnow)
    !call calindex(npv,pnow%sv,pnow%tmax,mpv)
    !--------------------------------------------------------
    if(pnow%typintpol=='cubic') then
        !call calvcubic(iproc,istep,npv,pnow,mpv,vc,vdc,vddc)
        !call calv_hermite(iproc,istep,npv,pnow,mpv,vc,vdc,vddc)
        call calv_quadratic(iproc,istep,npv,pnow,mpv,vdc,vddc)

        !call factor_cubic(npv,pnow%hv,pnow%e1v,pnow%e2v)
        !call inter_cubic(npv,pnow%ex,pnow%hv,pnow%e1v,pnow%e2v,pnow%cv)
        !call ffdfdd_cubic(npv,pnow%ex,pnow%sv,mpv,pnow%hv(mpv),pnow%tmax,pnow%cv,vc,vdc,vddc)
        !if(iproc==0 .and. vddc>0.d0) then
        !    write(pnow%ifile,'(a)') 'Not enough number of points for maximization, vddc>0'
        !    write(*         ,'(a)') 'Not enough number of points for maximization, vddc>0'
        !endif
        !alpha=min(5.d-1/abs(vddc),100.d0)
        alpha=min(5.d-1/abs(vddc),100.d0)
        vdd=vddc
    else
        call calvquintic(iproc,istep,npv,pnow,mpv,vq,vdq,vddq)
        !call factor_inter_quintic(npv,pnow%hv,pnow%ex,pnow%exd,pnow%a,pnow%b)
        !call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,pnow%tmax,vq,vdq,vddq)
        !if(iproc==0 .and. vddq>0.d0) then
        !    write(pnow%ifile,'(a)') 'Not enough number of points for maximization, vddq>0'
        !    write(*         ,'(a)') 'Not enough number of points for maximization, vddq>0'
        !endif
        alpha=min(5.d-1/abs(vddq),100.d0)
        vdd=vddq
    endif
    !if(.not. vdd<0.d0) then
    !    if(pnow%do_fill_ex_exd) then
    !        npv=2*npv
    !    else
    !        npv=np+pnow%ns2
    !    endif
    !    pnow%do_fill_ex_exd=.true.
    !    goto 1358
    !endif
    !----------------------------
    call atomic_dot(atoms,ft,ft,fnrm);fnrm=sqrt(fnrm)
    if(iproc==0) then
        !write(51,frt)           istep,0,npv,epot,vd,fnrm,vdd,alpha,pnow%tmax,pnow%tmax/pnow%s(np)
        write(*,frt2) 'fort51 ',istep,0,npv,epot,vd,fnrm,vdd,alpha,pnow%tmax,pnow%tmax/pnow%s(np)
    endif
    !--------------------------------------------------------
    !call epot_along_traj(istep,n,nr,np,x,npv,pnow,nproc,iproc,atoms,rst,inputs,ncount_bigdft)
    vdtol=max(pnow%vdtol*fnrm,5.d-5)
    !vdtol=max(pnow%vdtol*fnrm,1.d-8) !non-BigDFT
    if((abs(vd)<vdtol .and. vdd<0.d0) .or. (epot-max(pnow%exends(1),pnow%exends(2)))>1.d0) then 
    !if((abs(vd)<vdtol .and. vdd<0.d0) .or. (epot-max(pnow%exends(1),pnow%exends(2)))>1.d2) then !non-BigDFT
        pnow%npv=npv
        return
    endif
    vdold=vd
    ibad=0
    ibadold=0
    oneisideal=1.d0
    ll_inputs%inputPsiId=1
    do iter=1,50-(np+pnow%ns2+1)
        pnow%tmax=pnow%tmax-sign(min(abs(alpha*vd),0.05d0),alpha*vd)
        pnow%tmax=min(max(pnow%tmax,pnow%sv(0)),pnow%sv(npv))
        call calindex(np,pnow%s,pnow%tmax,mp,'caltmax2_3')
        !-------------------------
        call func(pnow%tmax,epot,vd,n,np,x,pnow,mp,xt,ft,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
        
        call calindex(npv,pnow%sv,pnow%tmax,mpv,'caltmax2_4')
        call insertpoint(npv,epot,vd,mpv,np,pnow)
        !-------------------------
        !call polish_sv(npv,pnow)
        !call calindex(npv,pnow%sv,pnow%tmax,mpv)
        !-------------------------
        !do ipv=0,npv
        !    write(1000*(iter+1)+istep,'(3es24.15)') pnow%sv(ipv),pnow%ex(ipv),pnow%exd(ipv)
        !enddo
        !close(1000*(iter+1)+istep)
        !-------------------------
        if(vd*vdold<0.d0) ibad=ibad+1
        if(ibad>ibadold) then
            oneisideal=0.7d0*oneisideal
            ibadold=ibad
        endif
        if(pnow%typintpol=='cubic') then
            !call calvcubic(iproc,istep,npv,pnow,mpv,vc,vdc,vddc)
            !call calv_hermite(iproc,istep,npv,pnow,mpv,vc,vdc,vddc)
            call calv_quadratic(iproc,istep,npv,pnow,mpv,vdc,vddc)

            !alpha=min(10.d-1/abs(vddc),200.d0)
            alpha=min(oneisideal/abs(vddc),200.d0)
            !alpha=min(8.d-1/abs(vddc),200.d0)
            vdd=vddc
        else
            call calvquintic(iproc,istep,npv,pnow,mpv,vq,vdq,vddq)
            !alpha=min(10.d-1/abs(vddq),200.d0)
            alpha=min(oneisideal/abs(vddq),200.d0)
            vdd=vddq
        endif
        vdold=vd
        !-------------------------
        call atomic_dot(atoms,ft,ft,fnrm);fnrm=sqrt(fnrm)
        if(iproc==0) then
            !write(51,frt)           istep,iter,npv,epot,vd,fnrm,vdd,alpha,pnow%tmax,pnow%tmax/pnow%s(np)
            write(*,frt2) 'fort51 ',istep,iter,npv,epot,vd,fnrm,vdd,alpha,pnow%tmax,pnow%tmax/pnow%s(np)
        endif
        vdtol=max(pnow%vdtol*fnrm,5.d-5)
        !vdtol=max(pnow%vdtol*fnrm,1.d-8) !non-BigDFT
        if((abs(vd)<vdtol .and. vdd<0.d0) .or. (epot-max(pnow%exends(1),pnow%exends(2)))>1.d0) then
        !if((abs(vd)<vdtol .and. vdd<0.d0) .or. (epot-max(pnow%exends(1),pnow%exends(2)))>1.d2) then !non-BigDFT
            exit
        endif
    enddo
    pnow%npv=npv
    ll_inputs%inputPsiId=0
    !if(iproc==0) &
    !write(54,'(3i5,e24.15,7e15.6)') istep,iter,npv,epot,vd,ft(1:2),vddq,alpha,pnow%tmax,pnow%tmax/pnow%s(np)
end subroutine caltmax2


subroutine polish_sv(npv,pnow)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::npv,ipv,lpvarr(100),lpvn,mpv
    type(parametersplinedsaddle)::pnow
    real(kind=8)::diff
    lpvn=0
    do ipv=1,npv-1
        !ipv=ipv+1
        !if(ipv==npv) 
        diff=abs(pnow%ex(ipv)-pnow%tmax)
        if(pnow%ex(ipv)<pnow%ex(ipv-1) .and. pnow%ex(ipv)<pnow%ex(ipv+1) .and. diff>1.d-5) then
            lpvn=lpvn+1
            lpvarr(lpvn)=ipv
        endif
    enddo
    if(lpvn==0) return
    do mpv=1,lpvn
        do ipv=lpvarr(mpv),npv-1
            pnow%sv(ipv)=pnow%sv(ipv+1)
            pnow%ex(ipv)=pnow%ex(ipv+1)
            pnow%exd(ipv)=pnow%exd(ipv+1)
        enddo
        npv=npv-1
        lpvarr(1:lpvn)=lpvarr(1:lpvn)-1
    enddo
    do ipv=1,npv
        pnow%hv(ipv)=pnow%sv(ipv)-pnow%sv(ipv-1)
    enddo
end subroutine polish_sv


subroutine write_v_of_t(iproc,istep,npv,pnow,icall,fn)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,npv,mpv,ipv
    type(parametersplinedsaddle)::pnow
    real(kind=8)::vc,vdc,vddc,dt,t
    !character(100)::text
    character(8)::filename
    character(*)::fn
    integer::icall
    !integer, save::icall=0
    !icall=icall+1
    call factor_cubic(npv,pnow%hv,pnow%e1v,pnow%e2v)
    call inter_cubic(npv,pnow%ex,pnow%hv,pnow%e1v,pnow%e2v,pnow%cv)
    write(filename,'(a4,a1,i3.3)') fn,'.',icall
    if(iproc==0) open(unit=2010,file=filename,status='replace')
    if(iproc==0) write(2010,'(a,i5)') '#istep = ',istep
    dt=pnow%sv(npv)/500.d0
    do ipv=0,500
        t=ipv*dt
        if(ipv==0) mpv=1
        if(ipv==500) mpv=npv
        if(ipv==500) t=pnow%sv(npv)
        if(.not. (ipv==0 .or. ipv==500)) call calindex(npv,pnow%sv,t,mpv,'write_v_of_t')
        call ffdfdd_cubic(npv,pnow%ex,pnow%sv,mpv,pnow%hv(mpv),t,pnow%cv,vc,vdc,vddc)
        !call ffdfdd_hermite(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,t,mpv,vc,vdc,vddc)
        if(iproc==0) write(2010,'(i5,4es19.10,i4)') ipv,t,vc,vdc,vddc,npv
    enddo
    if(iproc==0) close(2010)
end subroutine write_v_of_t


subroutine ffdfdd_hermite(n,s,h,e,ed,t,m,v,vd,vdd)
    implicit none
    integer::n,m
    real(kind=8)::s(0:n),h(n),e(0:n),ed(0:n),t,v,vd,vdd
    real(kind=8)::tt,a0,a1,a2,a3,hminv
    if(m<1 .or. m>n) stop 'ERROR: inconsistency in cubic Hermite spline m<1 .or. m>n'
    a0=e(m-1)
    a1=ed(m-1)*h(m)
    a2=3.d0*(-e(m-1)+e(m))-(2.d0*ed(m-1)+ed(m))*h(m)
    a3=2.d0*(e(m-1)-e(m))+(ed(m-1)+ed(m))*h(m)
    hminv=1.d0/h(m)
    tt=(t-s(m-1))*hminv
    v=a0+tt*(a1+tt*(a2+tt*a3))
    vd=(a1+tt*(2.d0*a2+tt*3.d0*a3))*hminv
    vdd=((2.d0*a2+tt*6.d0*a3))*hminv**2
end subroutine ffdfdd_hermite


subroutine fdd_quadratic(n,s,h,e,t,m,vd,vdd)
    implicit none
    integer::n,m
    real(kind=8)::s(0:n),h(n),e(0:n),t,vd,vdd
    !real(kind=8)::avg,w1,w2
    real(kind=8)::f0,f1,f2,h1,h2,a0,a1,a2,t1
    h1=h(m)
    h2=h(m+1)
    t1=s(m)
    f0=e(m-1)
    f1=e(m+0)
    f2=e(m+1)
    a0=(f1*(h1+h2)*(h1-t1)*(h2+t1)+t1*(f2*h1*(-h1+t1)+f0*h2*(h2+t1)))/(h1*h2*(h1+h2))
    a1=(f2*h1*(h1-2.d0*t1)-f1*(h1+h2)*(h1-h2-2*t1)-f0*h2*(h2+2.d0*t1))/(h1*h2*(h1+h2))
    a2=(f2*h1+f0*h2-f1*(h1+h2))/(h1*h2*(h1+h2))
    vd=a1+t*a2
    !vdd=a1+t*a2*2.d0
    vdd=a2*2.d0
end subroutine fdd_quadratic


subroutine calv_quadratic(iproc,istep,npv,pnow,mpv,vdc,vddc)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,npv,mpv
    type(parametersplinedsaddle)::pnow
    real(kind=8)::vdc,vddc
    character(100)::text
    if(mpv==npv) mpv=npv-1
    call fdd_quadratic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%tmax,mpv,vdc,vddc)
    text='Not enough number of points for maximization, vddc_quadratic>0'
    if(iproc==0 .and. vddc>0.d0) then
        write(pnow%ifile,'(i4,1x,a,i3)') istep,trim(text),npv
        write(*         ,'(i4,1x,a,i3)') istep,trim(text),npv
    endif
end subroutine calv_quadratic


subroutine calv_hermite(iproc,istep,npv,pnow,mpv,vc,vdc,vddc)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,npv,mpv
    type(parametersplinedsaddle)::pnow
    real(kind=8)::vc,vdc,vddc
    character(100)::text
    call ffdfdd_hermite(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%tmax,mpv,vc,vdc,vddc)
    text='Not enough number of points for maximization, vddc_hermite>0'
    if(iproc==0 .and. vddc>0.d0) then
        write(pnow%ifile,'(i4,1x,a,i3)') istep,trim(text),npv
        write(*         ,'(i4,1x,a,i3)') istep,trim(text),npv
    endif
end subroutine calv_hermite


subroutine guessinitialtmax_hermite(npv,pnow)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::npv,ipvt,istat,iroot,nroot,npvt
    type(parametersplinedsaddle)::pnow
    real(kind=8)::p1,p2,p3,t1,t2,hi,discriminant,v,vcmax,roots(50)
    real(kind=8), allocatable::svt(:),hvt(:),ext(:),exdt(:),e1vt(:),e2vt(:),cvt(:)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    npvt=npv  !10*npv
    allocate(svt(0:npvt+ndeb1),hvt(npvt+ndeb1),ext(0:npvt+ndeb1),exdt(0:npvt+ndeb1), &
        e1vt(npvt-1+ndeb1),e2vt(npvt-2+ndeb1),cvt(0:npvt+ndeb1),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating one of svt,hvt,ext,exdt,e1vt,e2vt,cvt.'
    call dmemocc(npvt+1,npvt+1+ndeb1,svt,'svt')
    call dmemocc(npvt,npvt+ndeb1,hvt,'hvt')
    call dmemocc(npvt+1,npvt+1+ndeb1,ext,'ext')
    call dmemocc(npvt+1,npvt+1+ndeb1,exdt,'exdt')
    call dmemocc(npvt-1,npvt-1+ndeb1,e1vt,'e1vt')
    call dmemocc(npvt-2,npvt-2+ndeb1,e2vt,'e2vt')
    call dmemocc(npvt+1,npvt+1+ndeb1,cvt,'cvt')
    !call factor_inter_quintic(npv,pnow%hv,pnow%ex,pnow%exd,pnow%a,pnow%b)
    ext(0)=pnow%ex(0)
    ext(npvt)=pnow%ex(npv)
    exdt(0)=pnow%exd(0)
    exdt(npvt)=pnow%exd(npv)
    svt(0)=pnow%sv(0)
    svt(npvt)=pnow%sv(npv)
    !dt=pnow%sv(npv)/(npvt)
    do ipvt=1,npvt-1
        !t=dt*ipvt
        svt(ipvt)=pnow%sv(ipvt)  !t
        hvt(ipvt)=pnow%hv(ipvt)  !svt(ipvt)-svt(ipvt-1)
        ext(ipvt)=pnow%ex(ipvt)
        exdt(ipvt)=pnow%exd(ipvt)
        !call calindex(npv,pnow%sv,t,mpv)
        !call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,t,ext(ipvt),t1,t2)
    enddo
    hvt(npvt)=pnow%hv(npv)  !svt(npvt)-svt(npvt-1)
    !ffdfdd_hermite(n,s,h,e,ed,t,m,v,vd,vdd)
    !call factor_cubic(npvt,hvt,e1vt,e2vt)
    !call inter_cubic(npvt,ext,hvt,e1vt,e2vt,cvt)
    nroot=0
    do ipvt=1,npvt
        hi=hvt(ipvt)
        !p3=(cvt(ipvt)-cvt(ipvt-1))/hi
        !p2=3.d0*(cvt(ipvt-1)*svt(ipvt)-cvt(ipvt)*svt(ipvt-1))/hi
        !p1=(3.d0*(cvt(ipvt)*svt(ipvt-1)**2-cvt(ipvt-1)*svt(ipvt)**2)+ext(ipvt)-ext(ipvt-1))/hi+hi*(cvt(ipvt-1)-cvt(ipvt))
        !p0=e(ipvt-1)
        p1=exdt(ipvt-1)*hi
        p2=3.d0*(-ext(ipvt-1)+ext(ipvt))-(2.d0*exdt(ipvt-1)+exdt(ipvt))*hi
        p3=2.d0*(ext(ipvt-1)-ext(ipvt))+(exdt(ipvt-1)+exdt(ipvt))*hi
        discriminant=p2**2-3.d0*p1*p3
        !write(*,*) 'discriminant',istep,ipvt,discriminant
        if(.not. discriminant<0.d0) then
            t1=(-p2+sqrt(discriminant))/(3.d0*p3)
            t2=(-p2-sqrt(discriminant))/(3.d0*p3)
            !t1=t1*(svt(ipvt)-svt(ipvt-1))+svt(ipvt-1)
            !t2=t2*(svt(ipvt)-svt(ipvt-1))+svt(ipvt-1)
            t1=t1*hvt(ipvt)+svt(ipvt-1)
            t2=t2*hvt(ipvt)+svt(ipvt-1)
            if(.not. t1<svt(ipvt-1) .and. .not. t1>svt(ipvt)) then
                !write(*,'(a,2i5,4f15.6)') 'roots: t1',istep,ipvt,t1,svt(ipvt-1),svt(ipvt),svt(npv)
                nroot=nroot+1
                roots(nroot)=t1
            endif
            if(.not. t2<svt(ipvt-1) .and. .not. t2>svt(ipvt)) then
                !write(*,'(a,2i5,4f15.6)') 'roots: t2',istep,ipvt,t2,svt(ipvt-1),svt(ipvt),svt(npv)
                nroot=nroot+1
                roots(nroot)=t2
            endif
        endif
    enddo
    if(nroot==0) stop 'ERROR: not enough points to interpolate v(t)'
    !--------------------------------------------------------
    vcmax=-1.d50
    do iroot=1,nroot
        !write(*,*) 'roots(iroot)',roots(iroot),svt(npv)
        call calindex(npvt,svt,roots(iroot),ipvt,'guessinitialtmax_hermite')
        !call ffdfdd_cubic(npvt,ext,svt,ipvt,hvt(ipvt),roots(iroot),cvt,vc,t1,t2)
        call ffdfdd_hermite(npvt,svt,hvt,ext,exdt,roots(iroot),ipvt,v,t1,t2)
        if(v>vcmax) then
            vcmax=v
            pnow%tmax=roots(iroot)
        endif
    enddo
    !call calindex(npv,pnow%sv,pnow%tmax,mpv)
    !call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,pnow%tmax,t1,t2,vddq)
    !if(vddq>0.d0) write(*,*) 'ERROR: vddq<0, use more points to find the maximum point.'
    deallocate(svt,hvt,ext,exdt,e1vt,e2vt,cvt,stat=istat)
    if(istat/=0) stop 'ERROR: failure deallocating one of svt,hvt,ext,e1vt,e2vt,cvt.'
end subroutine guessinitialtmax_hermite


subroutine calvcubic(iproc,istep,npv,pnow,mpv,vc,vdc,vddc)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,npv,mpv
    type(parametersplinedsaddle)::pnow
    real(kind=8)::vc,vdc,vddc
    character(100)::text
    call factor_cubic(npv,pnow%hv,pnow%e1v,pnow%e2v)
    call inter_cubic(npv,pnow%ex,pnow%hv,pnow%e1v,pnow%e2v,pnow%cv)
    call ffdfdd_cubic(npv,pnow%ex,pnow%sv,mpv,pnow%hv(mpv),pnow%tmax,pnow%cv,vc,vdc,vddc)
    text='Not enough number of points for maximization, vddc>0'
    if(iproc==0 .and. vddc>0.d0) then
        write(pnow%ifile,'(i4,1x,a,i3)') istep,trim(text),npv
        write(*         ,'(i4,1x,a,i3)') istep,trim(text),npv
    endif
end subroutine calvcubic


subroutine calvquintic(iproc,istep,npv,pnow,mpv,vq,vdq,vddq)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,npv,mpv
    type(parametersplinedsaddle)::pnow
    real(kind=8)::vq,vdq,vddq
    character(100)::text
    call factor_inter_quintic(npv,pnow%hv,pnow%ex,pnow%exd,pnow%a,pnow%b)
    call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv), &
        mpv,pnow%tmax,vq,vdq,vddq)
    text='Not enough number of points for maximization, vddq>0'
    if(iproc==0 .and. vddq>0.d0) then
        write(pnow%ifile,'(i4,1x,a,i3)') istep,trim(text),npv
        write(*         ,'(i4,1x,a,i3)') istep,trim(text),npv
    endif
end subroutine calvquintic


subroutine fill_ex_exd(istep,n,np,x,outends,npv,pnow,pold,xt,ft,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    type(DFT_global_output), dimension(2), intent(in) :: outends
    integer::istep,n,np,ip,mp,istat,i,npv,infocode
    real(kind=8)::x(n,0:np),xt(n),ft(n) !n(c) dt
    type(parametersplinedsaddle)::pnow,pold
    real(kind=8)::t1,tt,time1,time2
    real(kind=8), allocatable::tang(:)
    logical::move_this_coordinate
    integer::iat,ixyz
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    type(run_objects) :: runObj
    type(DFT_global_output) :: outs

    allocate(tang(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating tang.'
    call dmemocc(n,n+ndeb1,tang,'tang')
    call run_objects_nullify(runObj)
    call run_objects_associate(runObj, ll_inputs, atoms, rst)
    call init_global_output(outs, atoms%astruct%nat)
    pnow%ex(0)=pnow%exends(1)
    pnow%ex(npv)=pnow%exends(2)
    !test points along path will be distributed uniformly except one which is 
    !close to the pold%tmax will be replaced by pold%tmax
    !n(c) dt=pnow%s(np)/npv
    call estimate_sv(iproc,istep,np,npv,pnow,pold)
!    pnow%sv(0)=0.d0
!    pnow%sv(npv)=pnow%s(np)
!    if(istep>0) then
!        !-------------------------------------------------------------
!        !mpt=-1
!        !tnearest=1.d50
!        !do ip=1,npv-1
!        !    if(abs(pold%tmax/pold%s(np)*pnow%s(np)-dt*ip)<tnearest) then
!        !        mpt=ip
!        !        tnearest=abs(pold%tmax/pold%s(np)*pnow%s(np)-dt*ip)
!        !    endif
!        !enddo
!        !if(mpt==-1) stop 'ERROR: strange tmax of previous outer loop iteration'
!        !-------------------------------------------------------------
!        if(istep==6) then
!        do ip=0,pold%npv
!            write(301,'(3es24.15)') pold%sv(ip),pold%ex(ip),pold%exd(ip)
!        enddo
!        close(301)
!        endif
!        !pnow%ex(0)=
!        dtt=pold%sv(pold%npv)/npv
!        pnow%ex(0)=pold%ex(0)
!        pnow%ex(npv)=pold%ex(pold%npv)
!        pnow%exd(0)=pold%exd(0)
!        pnow%exd(npv)=pold%exd(pold%npv)
!        do ip=1,npv-1
!            tt=ip*dtt
!            call calindex(pold%npv,pold%sv,tt,mp)
!            tt1=(tt-pold%sv(mp-1))*(pold%ex(mp)-pold%ex(mp-1))/(pold%sv(mp)-pold%sv(mp-1))
!            pnow%ex(ip)=tt1+pold%ex(mp-1)
!            tt2=(tt-pold%sv(mp-1))*(pold%exd(mp)-pold%exd(mp-1))/(pold%sv(mp)-pold%sv(mp-1))
!            pnow%exd(ip)=tt2+pold%exd(mp-1)
!        enddo
!        exmax=-1.d50
!        exmin=+1.d50
!        do ip=0,npv
!            exmax=max(exmax,pnow%ex(ip))
!            exmin=min(exmin,pnow%ex(ip))
!        enddo
!        !exmax=exmax+1.d-2*(exmax-(pnow%ex(0)+pnow%ex(npv))*0.5d0)
!        !exmin=min(pnow%ex(0),pnow%ex(npv))
!        pnow%sv(0)=0.d0
!        do ip=1,npv
!            !pnow%hv(ip)=exp(2.d0*(pold%ex(ip)+pold%ex(ip-1)-2.d0*pold%ex(0)))
!            !tt=exp(exmax-pnow%ex(ip))
!            tt=exp(2.d0*(exmin-pnow%ex(ip))/(exmax-exmin))
!            pnow%sv(ip)=pnow%sv(ip-1)+tt/max(abs(pnow%exd(ip)),1.d-3)
!        enddo
!        tt=pnow%sv(npv)
!        pnow%sv(0:npv)=pnow%sv(0:npv)*pnow%s(np)/tt
!        if(istep==6) then
!        do ip=0,npv
!            write(302,'(3es24.15)') pnow%sv(ip),pnow%ex(ip),pnow%exd(ip)
!        enddo
!        close(302)
!        endif
!        !-------------------------------------------------------------
!    else
!        do ip=0,npv
!            pnow%sv(ip)=dt*ip
!        enddo
!    endif
    do ip=0,npv
        !pnow%sv(ip)=dt*ip
        !if(istep>0 .and. mpt==ip) pnow%sv(ip)=pold%tmax/pold%s(np)*pnow%s(np)
        if(ip==0) mp=1
        if(ip==npv) then
            mp=np
            pnow%sv(npv)=pnow%s(np)
        endif
        if(.not. (ip==0 .or. ip==npv)) call calindex(np,pnow%s,pnow%sv(ip),mp,'fill_ex_exd')
        if(ip>0) pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)


        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                pnow%y(0:np)=x(i,0:np)
                call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
                call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),pnow%sv(ip),pnow%c,xt(i),tang(i),t1)
            else
                xt(i)=x(i,0)
                tang(i)=0.d0
            endif
        enddo
        !do i=nr+1,n
        !    xt(i)=x(i,0)
        !    tang(i)=0.d0
        !enddo
        if(ip==0) then
            call atomic_dot(atoms,outends(1)%fxyz(1,1),tang,tt)
            pnow%exd(ip)=-tt
            !pnow%exd(ip)=-mydot(n,fends(1,1),tang) 
        elseif(ip==npv) then
            call atomic_dot(atoms,outends(2)%fxyz(1,1),tang,tt)
            pnow%exd(ip)=-tt
            !pnow%exd(ip)=-mydot(n,fends(1,2),tang) 
        else
            !call calenergyforces(n,xt,pnow%ex(ip),ft)
           call vcopy(n, xt(1), 1, runObj%atoms%astruct%rxyz(1,1), 1)
            call cpu_time(time1)
            call call_bigdft(runObj,outs,nproc,iproc,infocode)
            pnow%ex(ip) = outs%energy
            call cpu_time(time2)
            call vcopy(n, outs%fxyz(1,1), 1, ft(1), 1)
            ncount_bigdft=ncount_bigdft+1
            pnow%ncount_ll=pnow%ncount_ll+1
            pnow%time_ll=pnow%time_ll+(time2-time1)
            call atomic_dot(atoms,ft,tang,tt)
            pnow%exd(ip)=-tt
            !pnow%exd(ip)=-mydot(n,ft,tang) 
        endif
    enddo
    call run_objects_free_container(runObj)
    deallocate(tang,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating tang.'
    call deallocate_global_output(outs)
end subroutine fill_ex_exd

subroutine estimate_sv(iproc,istep,np,npv,pnow,pold)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::iproc,istep,np,npv,ip,mp,icycle,ncycle
    type(parametersplinedsaddle)::pnow,pold
    real(kind=8)::dt,tt,t1,t2
    !real(kind=8)::el,er,edl,edr
    real(kind=8)::ttmin,ttmax,emin,e1,e2,exo(0:100),exn(0:100),area,areatot
    pnow%sv(0)=0.d0
    pnow%sv(npv)=pnow%s(np)
    if(istep>0) then
        !-------------------------------------------------------------
        !mpt=-1
        !tnearest=1.d50
        !do ip=1,npv-1
        !    if(abs(pold%tmax/pold%s(np)*pnow%s(np)-dt*ip)<tnearest) then
        !        mpt=ip
        !        tnearest=abs(pold%tmax/pold%s(np)*pnow%s(np)-dt*ip)
        !    endif
        !enddo
        !if(mpt==-1) stop 'ERROR: strange tmax of previous outer loop iteration'
        !-------------------------------------------------------------
        !if(istep==6) then
        if(iproc==0) then
        do ip=0,pold%npv
            write(300+istep,'(3es24.15)') pold%sv(ip),pold%ex(ip),pold%exd(ip)
        enddo
        close(300+istep)
        endif
        !endif
        !-------------------------------------------------------------
        emin=min(pold%ex(0),pold%ex(pold%npv))
        exo(0:pold%npv)=pold%ex(0:pold%npv)-emin
        areatot=0.d0
        do ip=1,pold%npv
            areatot=areatot+0.5d0*(exo(ip-1)+exo(ip))*(pold%sv(ip)-pold%sv(ip-1))
        enddo
        dt=pold%sv(pold%npv)/npv
        pnow%sv(0)=0.d0
        do ip=1,npv-1
            pnow%sv(ip)=dt*ip
        enddo

        ncycle=100
        exn(0)=exo(0)
        exn(npv)=exo(pold%npv)
        do icycle=1,ncycle
            pnow%sv(npv)=pold%sv(pold%npv) !-1.d-10
            do ip=1,npv-1
                call calindex(pold%npv,pold%sv,pnow%sv(ip),mp,'estimate_sv')
                t1=pold%sv(mp-1) ; e1=exo(mp-1)
                t2=pold%sv(mp  ) ; e2=exo(mp  )
                exn(ip)=(e2-e1)/(t2-t1)*(pnow%sv(ip)-t1)+e1
            enddo
            do ip=1,npv
                area=0.5d0*(exn(ip-1)+exn(ip))*(pnow%sv(ip)-pnow%sv(ip-1))
                e1=(exn(npv)-exn(0))/pnow%sv(npv)*pnow%sv(ip-1)+exn(0)
                e2=(exn(npv)-exn(0))/pnow%sv(npv)*pnow%sv(ip  )+exn(0)
                !area=area-0.5d0*(e1+e2)*(pnow%sv(ip)-pnow%sv(ip-1))
                if(area<areatot/npv) then
                    pnow%hv(ip)=(pnow%sv(ip)-pnow%sv(ip-1))*1.02d0
                else
                    pnow%hv(ip)=(pnow%sv(ip)-pnow%sv(ip-1))*0.98d0
                endif
                !pnow%hv(ip)=1.d0/area
            enddo
            !pnow%sv(0)=0.d0
            do ip=1,npv
                pnow%sv(ip)=pnow%sv(ip-1)+pnow%hv(ip)
            enddo
            tt=pnow%sv(npv)
            pnow%sv(0:npv)=pnow%sv(0:npv)*pold%sv(pold%npv)/tt
            do ip=1,npv
                pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
            enddo
            ttmax=4.00d0*pold%sv(pold%npv)/npv
            ttmin=0.25d0*pold%sv(pold%npv)/npv
            do ip=1,npv
                tt=max(min(pnow%hv(ip),ttmax),ttmin)
                pnow%sv(ip)=pnow%sv(ip-1)+tt
            enddo
            tt=pnow%sv(npv)
            pnow%sv(0:npv)=pnow%sv(0:npv)*pold%sv(pold%npv)/tt
            do ip=1,npv
                pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
            enddo
        enddo

        do ip=1,npv-1
            call calindex(pold%npv,pold%sv,pnow%sv(ip),mp,'estimate_sv')
            t1=pold%sv(mp-1) ; e1=exo(mp-1)
            t2=pold%sv(mp  ) ; e2=exo(mp  )
            exn(ip)=(e2-e1)/(t2-t1)*(pnow%sv(ip)-t1)+e1
        enddo
        exn(0:npv)=exn(0:npv)+emin

        !-------------------------------------------------------------
!        dt=pold%sv(pold%npv)/npv
!        do ip=0,npv-1
!            pnow%sv(ip)=dt*ip
!        enddo
!        pnow%sv(npv)=pold%sv(pold%npv)-1.d-10
!        pnow%hv(1:npv)=0.d0
!        mpl=1
!        e_end_min=min(pold%ex(0),pold%ex(pold%npv))
!        do ip=0,pold%npv
!            pold%ex(ip)=pold%ex(ip)-e_end_min !+1.d-2
!        enddo
!        ttl=pold%ex(0) !-e_end_min
!        do ip=1,npv
!            call calindex(pold%npv,pold%sv,pnow%sv(ip),mpr)
!            do jp=mpl,mpr
!                el=pold%ex(jp-1) ; er=pold%ex(jp)
!                edl=pold%exd(jp-1) ; edr=pold%exd(jp)
!                h=pold%sv(jp)-pold%sv(jp-1) ; hinv=1.d0/h
!                a0=el
!                a1=edl*h
!                a2=3.d0*(-el+er)-(2.d0*edl+edr)*h
!                a3=2.d0*(el-er)+(edl+edr)*h
!                !t1=pold%sv(jp)
!                !if(jp==mpl) t1=pnow%sv(ip-1)
!                !t2=pold%sv(jp)
!                !if(jp==mpl) t1=pnow%sv(ip-1)
!                t1=max(pold%sv(jp-1),pnow%sv(ip-1))
!                t2=min(pold%sv(jp)  ,pnow%sv(ip)  )
!                tt1=(t1-pold%sv(jp-1))*hinv
!                tt2=(t2-pold%sv(jp-1))*hinv
!                vi1=tt1*(a0+tt1*(a1/2.d0+tt1*(a2/3.d0+tt1*a3/4.d0)))*h
!                vi2=tt2*(a0+tt2*(a1/2.d0+tt2*(a2/3.d0+tt2*a3/4.d0)))*h
!                !tt=
!                pnow%hv(ip)=pnow%hv(ip)+vi2-vi1
!            enddo
!            ttr=pnow%sv(ip)*(pold%ex(pold%npv)-pold%ex(0))/pold%sv(pold%npv)+pold%ex(0)
!            !ttr=ttr-e_end_min
!
!            !tt1=pold%ex(ip-1)-tt
!            !tt2=pold%ex(ip  )-tt
!            !pnow%hv(ip)=pnow%hv(ip)-0.5d0*(pold%ex(ip-1)+pold%ex(ip))*(pold%ex)
!            pnow%hv(ip)=pnow%hv(ip)-0.5d0*(ttl+ttr)*(pnow%sv(ip)-pnow%sv(ip-1))
!            mpl=mpr
!            ttl=ttr
!        enddo
!        !tt=4.d0*abs(pold%ex(0)-pold%ex(pold%npv))*pold%sv(pold%npv)/pold%npv
!        do ip=1,npv
!            if(pnow%hv(ip)<0.d0) write(*,*) 'WARNING: probably not sufficient points.'
!            !pnow%sv(ip)=pnow%sv(ip-1)+1.d0/max(abs(pnow%hv(ip))**0.8d0,tt)
!            pnow%sv(ip)=pnow%sv(ip-1)+1.d0/abs(pnow%hv(ip))**0.8d0
!            !pnow%sv(ip)=pnow%sv(ip-1)+1.d0/abs(pnow%hv(ip))**2
!            !pnow%sv(ip)=pnow%sv(ip-1)+exp(-abs(pnow%hv(ip)))
!        enddo
!        tt=pnow%sv(npv)
!        pnow%sv(0:npv)=pnow%sv(0:npv)*pnow%s(np)/tt
!        do ip=1,npv
!            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
!        enddo
!        tt=1.2d0*pnow%s(np)/npv
!        do ip=1,npv
!            pnow%hv(ip)=min(tt,pnow%hv(ip))
!            pnow%sv(ip)=pnow%sv(ip-1)+pnow%hv(ip)
!        enddo
!        tt=pnow%sv(npv)
!        pnow%sv(0:npv)=pnow%sv(0:npv)*pnow%s(np)/tt
!        do ip=1,npv
!            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
!        enddo
        !-------------------------------------------------------------
!        emin=min(pold%ex(0),pold%ex(pold%npv))+0.5d0*abs(pold%ex(pold%npv)-pold%ex(0))
!        !emin=max(pold%ex(0),pold%ex(pold%npv)) !-1.d0*abs(pold%ex(pold%npv)-pold%ex(0))
!        emax=-1.d10
!        do ip=0,pold%npv
!            emax=max(emax,pold%ex(ip))
!        enddo
!        dt=pold%sv(pold%npv)/npv
!        pnow%sv(0)=0.d0
!        do ip=1,npv-1
!            pnow%sv(ip)=dt*ip
!        enddo
!        do icycle=1,10
!        pnow%sv(npv)=pold%sv(pold%npv)-1.d-10
!        do ip=1,npv
!            !t=dt*ip
!            !if(ip==npv) then
!            !    mp=pold%npv
!            !else
!            !    call calindex(pold%npv,pold%sv,t,mp)
!            !endif
!            call calindex(pold%npv,pold%sv,pnow%sv(ip),mp,'estimate_sv')
!            t1=pold%sv(mp-1) ; e1=pold%ex(mp-1)
!            t2=pold%sv(mp  ) ; e2=pold%ex(mp  )
!            elin=(e2-e1)/(t2-t1)*(pnow%sv(ip)-t1)+e1
!            tt=(elin-emin)/(emax-emin)
!            !pnow%hv(ip)=exp(-5.d0*tt)
!            !pnow%hv(ip)=1.d0/(10.d-1+5.d-1*tt*8)
!            pnow%hv(ip)=1.d0/(1.d-6+5.d0*(tt**2)**0.4d0)
!        enddo
!        pnow%sv(0)=0.d0
!        do ip=1,npv
!            pnow%sv(ip)=pnow%sv(ip-1)+pnow%hv(ip)
!        enddo
!        tt=pnow%sv(npv)
!        pnow%sv(0:npv)=pnow%sv(0:npv)*pold%sv(pold%npv)/tt
!        do ip=1,npv
!            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
!        enddo
!        ttmax=3.0d0*pold%sv(pold%npv)/npv
!        ttmin=0.6d0*pold%sv(pold%npv)/npv
!        do ip=1,npv
!            tt=max(min(pnow%hv(ip),ttmax),ttmin)
!            pnow%sv(ip)=pnow%sv(ip-1)+tt
!        enddo
!        tt=pnow%sv(npv)
!        pnow%sv(0:npv)=pnow%sv(0:npv)*pold%sv(pold%npv)/tt
!        do ip=1,npv
!            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
!        enddo
!        enddo
!
!        tt=pnow%sv(npv)
!        pnow%sv(0:npv)=pnow%sv(0:npv)*pnow%s(np)/tt
!        do ip=1,npv
!            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
!        enddo
        !-------------------------------------------------------------
        !if(istep==6) then
        if(iproc==0) then
        do ip=0,npv
            write(400+istep,'(2es24.15)') pnow%sv(ip),exn(ip)
        enddo
        close(400+istep)
        endif
        !endif
        !-------------------------------------------------------------
        tt=pnow%sv(npv)
        pnow%sv(0:npv)=pnow%sv(0:npv)*pnow%s(np)/tt
        do ip=1,npv
            pnow%hv(ip)=pnow%sv(ip)-pnow%sv(ip-1)
        enddo
    else
        dt=pnow%s(np)/npv
        do ip=0,npv
            pnow%sv(ip)=dt*ip
        enddo
    endif
end subroutine estimate_sv


!subroutine epot_along_traj(istep,n,nr,np,x,npv,pnow,nproc,iproc,atoms,rst,inputs,ncount_bigdft)
!    use module_base
!    use module_interfaces
!    use module_types
!    use modulesplinedsaddle, only:parametersplinedsaddle
!    implicit none
!    integer, intent(in) :: nproc,iproc
!    type(atoms_data), intent(inout) :: atoms
!    type(input_variables), intent(inout) :: inputs
!    type(restart_objects), intent(inout) :: rst
!    integer, intent(inout) :: ncount_bigdft
!    integer::istep,n,nr,np,mpv,istat,i,npv,mp,ip,infocode
!    real(kind=8)::x(n,0:np),epot,dt,t,vc,vq,t1,t2,fnoise
!    type(parametersplinedsaddle)::pnow
!    character(3)::fn
!    character(20)::filename
!    real(kind=8), allocatable::xt(:)
!    real(kind=8), allocatable::ft(:),x_bigdft(:)
!    logical::move_this_coordinate
!    integer::iat,ixyz
!    integer, parameter::ndeb1=0,ndeb2=0
!    allocate(x_bigdft(n+ndeb1))
!    call dmemocc(n,n+ndeb1,x_bigdft,'x_bigdft')
!    allocate(xt(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xt.'
!    call dmemocc(n,n+ndeb1,xt,'xt')
!    allocate(ft(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ft.'
!    call dmemocc(n,n+ndeb1,ft,'ft')
!    if(iproc==0) then
!        write(fn,'(i3.3)') istep
!        filename='energy'//fn
!        open(unit=1358,file=filename,status='replace')
!    endif
!    mpv=1;dt=pnow%sv(npv)/20.d0  !*(1.d0-epsilon(dt))
!    call factor_cubic(npv,pnow%hv,pnow%e1v,pnow%e2v)
!    call inter_cubic(npv,pnow%ex,pnow%hv,pnow%e1v,pnow%e2v,pnow%cv)
!    call factor_inter_quintic(npv,pnow%hv,pnow%ex,pnow%exd,pnow%a,pnow%b)
!    do ip=0,19
!        t=dt*ip 
!        if(ip==20) then
!            mpv=npv
!            mp=np
!        else
!            call calindex(npv,pnow%sv,t,mpv)
!            call calindex(np,pnow%s,t,mp)
!        endif
!        call ffdfdd_cubic(npv,pnow%ex,pnow%sv,mpv,pnow%hv(mpv),t,pnow%cv,vc,t1,t2)
!        call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,t,vq,t1,t2)
!        do i=1,n
!            iat=(i-1)/3+1
!            ixyz=mod(i-1,3)+1
!            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
!                pnow%y(0:np)=x(i,0:np)
!                call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
!                call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),t,pnow%c,xt(i),t1,t2)
!            else
!                xt(i)=x(i,0)
!            endif
!        enddo
!        !do i=nr+1,n
!        !    xt(i)=x(i,0)
!        !enddo
!        !call calenergyforces(n,xt,epot,ft)
!        x_bigdft(1:n)=xt(1:n)
!        call call_bigdft(nproc,iproc,atoms,x_bigdft,inputs,epot,ft,fnoise,rst,infocode)
!        ncount_bigdft=ncount_bigdft+1
!        if(iproc==0) then
!            write(1358,'(4e20.10)') t/pnow%sv(npv),vc,vq,epot
!        endif
!    enddo
!    if(iproc==0) then
!        close(1358)
!    endif
!    deallocate(x_bigdft)
!    deallocate(xt,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xt.'
!    deallocate(ft,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ft.'
!end subroutine epot_along_traj


subroutine insertpoint(npv,epot,vd,mpv,np,pnow)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::npv,mpv,np,ip,lpv
    real(kind=8)::epot,vd
    type(parametersplinedsaddle)::pnow
    if(pnow%tmax-pnow%sv(mpv-1)<(pnow%sv(mpv)-pnow%tmax)) then
        lpv=mpv-1
    else
        lpv=mpv
    endif
    if(lpv==0) lpv=mpv
    if(lpv==npv) lpv=mpv-1
    !if((pnow%tmax-pnow%sv(mpv-1))<pnow%htol*pnow%s(np)) then
    !    pnow%sv(mpv-1)=pnow%tmax
    !    pnow%ex(mpv-1)=epot
    !    pnow%exd(mpv-1)=-vd
    !    if(mpv>1) pnow%hv(mpv-1)=pnow%sv(mpv-1)-pnow%sv(mpv-2)
    !    pnow%hv(mpv)=pnow%sv(mpv)-pnow%sv(mpv-1)
    !elseif((pnow%sv(mpv)-pnow%tmax)<pnow%htol*pnow%s(np)) then
    !    pnow%sv(mpv)=pnow%tmax
    !    pnow%ex(mpv)=epot
    !    pnow%exd(mpv)=-vd
    !    pnow%hv(mpv)=pnow%sv(mpv)-pnow%sv(mpv-1)
    !    if(mpv<npv) pnow%hv(mpv+1)=pnow%sv(mpv+1)-pnow%sv(mpv)
    if(abs(pnow%sv(lpv)-pnow%tmax)<pnow%htol*pnow%s(np)) then
        pnow%sv(lpv)=pnow%tmax
        pnow%ex(lpv)=epot
        pnow%exd(lpv)=-vd
        pnow%hv(lpv)=pnow%sv(lpv)-pnow%sv(lpv-1)
        pnow%hv(lpv+1)=pnow%sv(lpv+1)-pnow%sv(lpv)
        mpv=lpv
    else
        do ip=npv,mpv,-1
            pnow%hv(ip+1)=pnow%hv(ip)
            pnow%sv(ip+1)=pnow%sv(ip)
            pnow%ex(ip+1)=pnow%ex(ip)
            pnow%exd(ip+1)=pnow%exd(ip)
        enddo
        pnow%sv(mpv)=pnow%tmax
        pnow%ex(mpv)=epot
        pnow%exd(mpv)=-vd
        pnow%hv(mpv)=pnow%sv(mpv)-pnow%sv(mpv-1)
        pnow%hv(mpv+1)=pnow%sv(mpv+1)-pnow%sv(mpv)
        npv=npv+1
    endif
end subroutine insertpoint


subroutine guessinitialtmax_cubic(npv,pnow)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::npv,ipvt,istat,iroot,nroot,npvt
    type(parametersplinedsaddle)::pnow
    real(kind=8)::p1,p2,p3,t1,t2,hi,discriminant,vc,vcmax,roots(50)
    real(kind=8), allocatable::svt(:),hvt(:),ext(:),e1vt(:),e2vt(:),cvt(:)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    npvt=npv  !10*npv
    allocate(svt(0:npvt+ndeb1),hvt(npvt+ndeb1),ext(0:npvt+ndeb1),e1vt(npvt-1+ndeb1),e2vt(npvt-2+ndeb1),cvt(0:npvt+ndeb1),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating one of svt,hvt,ext,e1vt,e2vt,cvt.'
    call dmemocc(npvt+1,npvt+1+ndeb1,svt,'svt')
    call dmemocc(npvt,npvt+ndeb1,hvt,'hvt')
    call dmemocc(npvt+1,npvt+1+ndeb1,ext,'ext')
    call dmemocc(npvt-1,npvt-1+ndeb1,e1vt,'e1vt')
    call dmemocc(npvt-2,npvt-2+ndeb1,e2vt,'e2vt')
    call dmemocc(npvt+1,npvt+1+ndeb1,cvt,'cvt')
    !call factor_inter_quintic(npv,pnow%hv,pnow%ex,pnow%exd,pnow%a,pnow%b)
    ext(0)=pnow%ex(0)
    ext(npvt)=pnow%ex(npv)
    svt(0)=pnow%sv(0)
    svt(npvt)=pnow%sv(npv)
    !dt=pnow%sv(npv)/(npvt)
    do ipvt=1,npvt-1
        !t=dt*ipvt
        svt(ipvt)=pnow%sv(ipvt)  !t
        hvt(ipvt)=pnow%hv(ipvt)  !svt(ipvt)-svt(ipvt-1)
        ext(ipvt)=pnow%ex(ipvt)
        !call calindex(npv,pnow%sv,t,mpv)
        !call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,t,ext(ipvt),t1,t2)
    enddo
    hvt(npvt)=pnow%hv(npv)  !svt(npvt)-svt(npvt-1)
    call factor_cubic(npvt,hvt,e1vt,e2vt)
    call inter_cubic(npvt,ext,hvt,e1vt,e2vt,cvt)
    nroot=0
    do ipvt=1,npvt
        hi=hvt(ipvt)
        p3=(cvt(ipvt)-cvt(ipvt-1))/hi
        p2=3.d0*(cvt(ipvt-1)*svt(ipvt)-cvt(ipvt)*svt(ipvt-1))/hi
        p1=(3.d0*(cvt(ipvt)*svt(ipvt-1)**2-cvt(ipvt-1)*svt(ipvt)**2)+ext(ipvt)-ext(ipvt-1))/hi+hi*(cvt(ipvt-1)-cvt(ipvt))
        discriminant=p2**2-3.d0*p1*p3
        !write(*,*) 'discriminant',istep,ipvt,discriminant
        if(.not. discriminant<0.d0) then
            t1=(-p2+sqrt(discriminant))/(3.d0*p3)
            t2=(-p2-sqrt(discriminant))/(3.d0*p3)
            if(.not. t1<svt(ipvt-1) .and. .not. t1>svt(ipvt)) then
                !write(*,'(a,2i5,4f15.6)') 'roots: t1',istep,ipvt,t1,svt(ipvt-1),svt(ipvt),svt(npv)
                nroot=nroot+1
                roots(nroot)=t1
            endif
            if(.not. t2<svt(ipvt-1) .and. .not. t2>svt(ipvt)) then
                !write(*,'(a,2i5,4f15.6)') 'roots: t2',istep,ipvt,t2,svt(ipvt-1),svt(ipvt),svt(npv)
                nroot=nroot+1
                roots(nroot)=t2
            endif
        endif
    enddo
    if(nroot==0) stop 'ERROR: not enough points to interpolate v(t)'
    !--------------------------------------------------------
    vcmax=-1.d50
    do iroot=1,nroot
        !write(*,*) 'roots(iroot)',roots(iroot),svt(npv)
        call calindex(npvt,svt,roots(iroot),ipvt,'guessinitialtmax_cubic')
        call ffdfdd_cubic(npvt,ext,svt,ipvt,hvt(ipvt),roots(iroot),cvt,vc,t1,t2)
        if(vc>vcmax) then
            vcmax=vc
            pnow%tmax=roots(iroot)
        endif
    enddo
    !call calindex(npv,pnow%sv,pnow%tmax,mpv)
    !call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,pnow%tmax,t1,t2,vddq)
    !if(vddq>0.d0) write(*,*) 'ERROR: vddq<0, use more points to find the maximum point.'
    deallocate(svt,hvt,ext,e1vt,e2vt,cvt,stat=istat)
    if(istat/=0) stop 'ERROR: failure deallocating one of svt,hvt,ext,e1vt,e2vt,cvt.'
end subroutine guessinitialtmax_cubic


subroutine guessinitialtmax_quintic(npv,pnow,iproc)
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer::npv,mpv,ipvt,istat,iroot,nroot,npvt,iproc
    type(parametersplinedsaddle)::pnow
    real(kind=8)::p1,p2,p3,t1,t2,hi,discriminant,t,dt,vc,vcmax,vddq,roots(50)
    real(kind=8), allocatable::svt(:),hvt(:),ext(:),e1vt(:),e2vt(:),cvt(:)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    npvt=10*npv
    allocate(svt(0:npvt+ndeb1),hvt(npvt+ndeb1),ext(0:npvt+ndeb1),e1vt(npvt-1+ndeb1),e2vt(npvt-2+ndeb1),cvt(0:npvt+ndeb1),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating one of svt,hvt,ext,e1vt,e2vt,cvt.'
    call dmemocc(npvt+1,npvt+1+ndeb1,svt,'svt')
    call dmemocc(npvt,npvt+ndeb1,hvt,'hvt')
    call dmemocc(npvt+1,npvt+1+ndeb1,ext,'ext')
    call dmemocc(npvt-1,npvt-1+ndeb1,e1vt,'e1vt')
    call dmemocc(npvt-2,npvt-2+ndeb1,e2vt,'e2vt')
    call dmemocc(npvt+1,npvt+1+ndeb1,cvt,'cvt')
    call factor_inter_quintic(npv,pnow%hv,pnow%ex,pnow%exd,pnow%a,pnow%b)
    ext(0)=pnow%ex(0)
    ext(npvt)=pnow%ex(npv)
    svt(0)=pnow%sv(0)
    svt(npvt)=pnow%sv(npv)
    dt=pnow%sv(npv)/(npvt)
    do ipvt=1,npvt-1
        t=dt*ipvt
        svt(ipvt)=t
        hvt(ipvt)=svt(ipvt)-svt(ipvt-1)
        call calindex(npv,pnow%sv,t,mpv,'guessinitialtmax_quintic_1')
        call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,t,ext(ipvt),t1,t2)
    enddo
    hvt(npvt)=svt(npvt)-svt(npvt-1)
    call factor_cubic(npvt,hvt,e1vt,e2vt)
    call inter_cubic(npvt,ext,hvt,e1vt,e2vt,cvt)
    nroot=0
    do ipvt=1,npvt
        hi=hvt(ipvt)
        p3=(cvt(ipvt)-cvt(ipvt-1))/hi
        p2=3.d0*(cvt(ipvt-1)*svt(ipvt)-cvt(ipvt)*svt(ipvt-1))/hi
        p1=(3.d0*(cvt(ipvt)*svt(ipvt-1)**2-cvt(ipvt-1)*svt(ipvt)**2)+ext(ipvt)-ext(ipvt-1))/hi+hi*(cvt(ipvt-1)-cvt(ipvt))
        discriminant=p2**2-3.d0*p1*p3
        !write(*,*) 'discriminant',istep,ipvt,discriminant
        if(.not. discriminant<0.d0) then
            t1=(-p2+sqrt(discriminant))/(3.d0*p3)
            t2=(-p2-sqrt(discriminant))/(3.d0*p3)
            if(.not. t1<svt(ipvt-1) .and. .not. t1>svt(ipvt)) then
                !write(*,'(a,2i5,4f15.6)') 'roots: t1',istep,ipvt,t1,svt(ipvt-1),svt(ipvt),svt(npv)
                nroot=nroot+1
                roots(nroot)=t1
            endif
            if(.not. t2<svt(ipvt-1) .and. .not. t2>svt(ipvt)) then
                !write(*,'(a,2i5,4f15.6)') 'roots: t2',istep,ipvt,t2,svt(ipvt-1),svt(ipvt),svt(npv)
                nroot=nroot+1
                roots(nroot)=t2
            endif
        endif
    enddo
    if(nroot==0) stop 'ERROR: not enough points to interpolate v(t)'
    !--------------------------------------------------------
    vcmax=-1.d50
    do iroot=1,nroot
        !write(*,*) 'roots(iroot)',roots(iroot),svt(npv)
        call calindex(npvt,svt,roots(iroot),ipvt,'guessinitialtmax_quintic_2')
        call ffdfdd_cubic(npvt,ext,svt,ipvt,hvt(ipvt),roots(iroot),cvt,vc,t1,t2)
        if(vc>vcmax) then
            vcmax=vc
            pnow%tmax=roots(iroot)
        endif
    enddo
    call calindex(npv,pnow%sv,pnow%tmax,mpv,'guessinitialtmax_quintic_3')
    call ffdfdd_quintic(npv,pnow%sv,pnow%hv,pnow%ex,pnow%exd,pnow%a(mpv),pnow%b(mpv),mpv,pnow%tmax,t1,t2,vddq)
    if(vddq>0.d0 .and. iproc==0) then
            write(pnow%ifile,*) 'ERROR: vddq<0, use more points to find the maximum point.'
            write(*,*) 'ERROR: vddq<0, use more points to find the maximum point.'
    endif
    deallocate(svt,hvt,ext,e1vt,e2vt,cvt,stat=istat)
    if(istat/=0) stop 'ERROR: failure deallocating one of svt,hvt,ext,e1vt,e2vt,cvt.'
end subroutine guessinitialtmax_quintic


subroutine factor_inter_quintic(n,h,y,d,a,b)
    implicit none
    integer::n,i,j,k,info
    integer, parameter::kl=2,ku=2
    real(kind=8)::h(n),y(0:n),d(0:n),a(n),b(n),t1,t2,t3,t4
    real(kind=8), allocatable::c(:),v(:,:)
    integer, allocatable::ipiv(:)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    allocate(c(2*n+ndeb1),v(2*kl+ku+1,2*n+ndeb1),ipiv(2*n+ndeb1))
    call dmemocc(2*n,2*n+ndeb1,c,'c')
    call dmemocc((2*kl+ku+1)*(2*n),(2*kl+ku+1)*(2*n+ndeb1),v,'v')
    call imemocc(2*n,2*n+ndeb1,ipiv,'ipiv')
    c(1)=-6.d0*d(0)*h(1)-6.d0*d(1)*h(1)-12.d0*y(0)+12.d0*y(1)
    do i=1,n-1
        t1=-2.d0*(2.d0*d(i)+d(i+1))*h(i)-2.d0*(2.d0*d(i)+d(i-1))*h(i+1)
        c(2*i)=t1*h(i)*h(i+1)-6.d0*h(i+1)**2*(-y(i)+y(i-1))+6.d0*h(i)**2*(-y(i)+y(i+1))
        t2=2.d0*(2.d0*d(i)+d(i+1))*h(i)**2  +2.d0*(d(i)+2.d0*d(i+1))*h(i)**2
        t3=2.d0*(2.d0*d(i)+d(i-1))*h(i+1)**2+2.d0*(d(i)+2.d0*d(i-1))*h(i+1)**2
        t4=12.d0*h(i+1)**3*(y(i)-y(i-1))+12.d0*h(i)**3*(y(i)-y(i+1))
        c(2*i+1)=(t2-t3)*h(i)*h(i+1)+t4

        !c(2*i)=-2.d0*(2.d0*d(i) + d(i+1))*h(i)**2*h(i+1) - 2.d0*(2.d0*d(i) + d(i-1))*h(i)*h(i+1)**2 - &
        !    6.d0*h(i+1)**2*(-y(i) + y(i-1)) + 6.d0*h(i)**2*(-y(i) + y(i+1))
        !c(2*i+1)=2.d0*(2*d(i) + d(i+1))*h(i)**3*h(i+1) + 2.d0*(d(i) + 2.d0*d(i+1))*h(i)**3*h(i+1) - &
        !    2.d0*(2.d0*d(i) + d(i-1))*h(i)*h(i+1)**3 - 2.d0*(d(i) + 2.d0*d(i-1))*h(i)*h(i+1)**3 + &
        !    6.d0*h(i+1)**3*(y(i) - y(i-1)) - 6.d0*h(i+1)**3*(-y(i) + y(i-1)) + 6.d0*h(i)**3*(y(i) - &
        !    y(i+1)) - 6.d0*h(i)**3*(-y(i) + y(i+1))
    enddo
    c(2*n)=-6.d0*d(n)*h(n)-6.d0*d(n-1)*h(n)+12.d0*y(n)-12.d0*y(n-1)
    j=1
    v(kl+1,j)=0.d0
    v(kl+2,j)=0.d0
    v(kl+3,j)=-18.d0*h(1)**2
    v(kl+4,j)=6.d0*h(1)**2*h(2)**2
    v(kl+5,j)=42.d0*h(1)**2*h(2)**3
    j=2
    v(kl+1,j)=0.d0
    v(kl+2,j)=-42.d0*h(1)**2
    v(kl+3,j)=4.d0*h(1)**2*h(2)**2
    v(kl+4,j)=18.d0*h(1)**2*h(2)**3
    v(kl+5,j)=0.d0
    do k=2,n-1
        j=2*k-1
        v(kl+1,j)=0.d0
        v(kl+2,j)=-4.d0*h(k-1)**2*h(k)**2
        v(kl+3,j)=18.d0*h(k-1)**3*h(k)**2
        v(kl+4,j)=6.d0*h(k+1)**2*h(k)**2
        v(kl+5,j)=42.d0*h(k+1)**3*h(k)**2
        j=2*k
        v(kl+1,j)=-6.d0*h(k-1)**2*h(k)**2
        v(kl+2,j)=42.d0*h(k-1)**3*h(k)**2
        v(kl+3,j)=4.d0*h(k+1)**2*h(k)**2
        v(kl+4,j)=18.d0*h(k+1)**3*h(k)**2
        v(kl+5,j)=0.d0
    enddo
    j=2*n-1
    v(kl+1,j)=0.d0
    v(kl+2,j)=-4.d0*h(n-1)**2*h(n)**2
    v(kl+3,j)=18.d0*h(n-1)**3*h(n)**2
    v(kl+4,j)=42.d0*h(n)**2
    v(kl+5,j)=0.d0
    j=2*n
    v(kl+1,j)=-6.d0*h(n-1)**2*h(n)**2
    v(kl+2,j)=42.d0*h(n-1)**3*h(n)**2
    v(kl+3,j)=18.d0*h(n)**2
    v(kl+4,j)=0.d0
    v(kl+5,j)=0.d0
    call dgbsv(2*n,kl,ku,1,v,2*kl+ku+1,ipiv,c,2*n,info)
    if(info/=0) write(*,*) 'ERROR: solution of dgbsv failed: info',info
    do i=1,n
        a(i)=c(2*i-1)
        b(i)=c(2*i)
    enddo
    deallocate(c,v,ipiv)
end subroutine factor_inter_quintic


subroutine ffdfdd_quintic(n,s,h,y,d,a,b,i,t,v,vd,vdd)
    implicit none
    integer::n,i
    real(kind=8)::y(0:n),d(0:n),s(0:n),h(n),a,b,t,v,vd,vdd
    real(kind=8)::p5,p4,p3,p2,p1,p0,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
    real(kind=8)::t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26
    real(kind=8)::tt !n(c) f1,fd1,f2,fd2
    t0=-y(i) + y(i-1)
    t1=2.d0*s(i) + s(i-1)
    t2=s(i) + 2.d0*s(i-1)
    t3=2.d0*d(i) + d(i-1)
    t4=d(i) + 2.d0*d(i-1)
    t5=7.d0*s(i)**3 + 2.d0*s(i-1)**3
    t6=2.d0*s(i)**3 + 7.d0*s(i-1)**3
    t7=3.d0*s(i)**5   + 2.d0*h(i)**4*t1 - h(i)**2*t5
    t8=-3.d0*s(i-1)**5 - 2.d0*h(i)**4*t2 + h(i)**2*t6
    t9=t4*s(i)**3 + t3*s(i-1)**3
    t10=d(i-1)*t1 + d(i)*t2
    t11=3.d0*(s(i)**3 + s(i-1)**3)
    t12=-t2*y(i) + t1*y(i-1)
    p0=(-h(i)*t9+h(i)**3*t10+b*t7+a*t8-t11*t0+3.d0*h(i)**2*t12)/(3.d0*h(i)**3) 
    t13=d(i) + d(i-1)
    t14=t4*s(i)**2 + t3*s(i-1)**2
    t15=7.d0*s(i)**2 + 2.d0*s(i-1)**2
    t16=2.d0*s(i)**2 + 7.d0*s(i-1)**2
    t17=-2.d0*h(i)**4 - 5.d0*s(i)**4 + h(i)**2*t15
    t18=2.d0*h(i)**4 + 5.d0*s(i-1)**4 - h(i)**2*t16
    t19=3.d0*(s(i)**2 + s(i-1)**2)*t0
    t20=y(i) - y(i-1)
    p1=(-t13*h(i)**3+h(i)*t14+b*t17+a*t18+3.d0*h(i)**2*t20+t19)/h(i)**3
    t21=7.d0*s(i) + 2.d0*s(i-1)
    t22=2.d0*s(i) + 7.d0*s(i-1)
    t23=10.d0*s(i)**3 - h(i)**2*t21
    t24=-10.d0*s(i-1)**3 + h(i)**2*t22
    t25=3.d0*(s(i) + s(i-1))*t0
    p2=(-h(i)*t10+b*t23+a*t24-t25)/h(i)**3
    t26=-3.d0*h(i)**2 + 10.d0*s(i-1)**2
    p3=(h(i)*(t13+3.d0*b*h(i))-10.d0*b*s(i)**2+a*t26+2.d0*t0)/h(i)**3
    p4=(5.d0*b*s(i) - 5.d0*a*s(i-1))/h(i)**3
    p5=(a - b)/h(i)**3

    v=((((p5*t+p4)*t+p3)*t+p2)*t+p1)*t+p0
    vd=(((5.d0*p5*t+4.d0*p4)*t+3.d0*p3)*t+2.d0*p2)*t+p1
    vdd=((20.d0*p5*t+12.d0*p4)*t+6.d0*p3)*t+2.d0*p2

    tt=s(i-1)
    !n(c) f1=((((p5*tt+p4)*tt+p3)*tt+p2)*tt+p1)*tt+p0
    !n(c) fd1=(((5.d0*p5*tt+4.d0*p4)*tt+3.d0*p3)*tt+2.d0*p2)*tt+p1
    tt=s(i)
    !n(c) f2=((((p5*tt+p4)*tt+p3)*tt+p2)*tt+p1)*tt+p0
    !n(c) fd2=(((5.d0*p5*tt+4.d0*p4)*tt+3.d0*p3)*tt+2.d0*p2)*tt+p1
    !write(91,'(7e20.10)') abs(y(i-1)-f1),abs(d(i-1)-fd1),abs(y(i)-f2),abs(d(i)-fd2),s(i-1),t,s(i)

end subroutine ffdfdd_quintic


subroutine calindex(np,s,t,ip,strcall)
    implicit none
    integer::np,ip
    real(kind=8)::s(0:np),t
    character(*)::strcall
    character(100)::strprt
    do ip=1,np
        if(t<s(ip)) exit
    enddo
    strprt='ERROR: '//'in '//trim(strcall)//' inconsistent string length: ip,np,t,snp '
    if(ip==np+1) then
        !write(*,'(a,2i4,2es20.10)') 'ERROR: inconsistent string length: ip,np ',ip,np,t,s(ip)
        write(*,'(a,2i4,2es20.10)') trim(strprt),ip,np,t,s(ip)
        stop
    endif
end subroutine calindex


subroutine prepdd(atoms,n,np,x,e1,e2,h,s,mp,tmax,dd)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,np,mp,istat,i,info,ip,j,jp
    real(kind=8)::x(n,0:np),e1(np-1),e2(np-2),h(np),s(0:np),dd(n,n,np-1),tmax
    real(kind=8), allocatable::ainv(:,:)
    real(kind=8), allocatable::cd1(:)
    real(kind=8), allocatable::cd2(:)
    !real(kind=8), allocatable::cd3(:)
    !real(kind=8), allocatable::cd4(:)
    real(kind=8), allocatable::c(:)
    real(kind=8), allocatable::ddt(:)
    real(kind=8), allocatable::yi(:)
    real(kind=8), allocatable::yj(:)
    logical::move_this_coordinate
    integer::ixyz,iat,jxyz,jat
    integer, parameter::ndeb1=0
    !integer, parameter::ndeb2=0
    allocate(cd1(np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating cd1.'
    call dmemocc(np-1,np-1+ndeb1,cd1,'cd1')
    allocate(cd2(np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating cd2.'
    call dmemocc(np-1,np-1+ndeb1,cd2,'cd2')
    !allocate(cd3(np-1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating cd3.'
    !allocate(cd4(np-1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating cd4.'
    allocate(c(0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating c.'
    call dmemocc(np+1,np+1+ndeb1,c,'c')
    allocate(ddt(np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ddt.'
    call dmemocc(np-1,np-1+ndeb1,ddt,'ddt')
    allocate(yi(0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating yi.'
    call dmemocc(np+1,np+1+ndeb1,yi,'yi')
    allocate(yj(0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating yj.'
    call dmemocc(np+1,np+1+ndeb1,yj,'yj')
    allocate(ainv(np-1,np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ainv.'
    call dmemocc((np-1)*(np-1),(np-1)*(np-1+ndeb1),ainv,'ainv')
    ainv(1:np-1,1:np-1)=0.d0
    do ip=1,np-2
        !ainv(ip,ip)=e1(ip)
        !ainv(ip,ip+1)=e2(ip)
        ainv(ip,ip)=2.d0*(h(ip+1)+h(ip))
        ainv(ip,ip+1)=h(ip+1)
        !if(ip>1) ainv(ip,ip-1)=e2(ip)
    enddo
    ainv(np-1,np-1)=2.d0*(h(np)+h(np-1)) 
    call dpotrf('U',np-1,ainv,np-1,info)
    if(info/=0) write(*,*) 'ERROR, dpotrf: factorization failed: info',info
    call dpotri('U',np-1,ainv,np-1,info)
    if(info/=0) write(*,*) 'ERROR, dpotri: inversion failed: info',info
    do ip=1,np-1
        do jp=1,ip-1
            ainv(ip,jp)=ainv(jp,ip)
        enddo
    enddo
    dd(1:n,1:n,1:np-1)=0.d0
    do j=1,n
        jat=(j-1)/3+1
        jxyz=mod(j-1,3)+1
        if(.not. move_this_coordinate(atoms%astruct%ifrztyp(jat),jxyz)) cycle
        yj(0:np)=x(j,0:np)
        call inter_cubic(np,yj,h,e1,e2,c)
        !t1=0.d0
        !do ip=1,np-1
        !    t1=t1+ainv(1,ip)*((yj(ip+1)-yj(ip))/h(ip+1)-(yj(ip)-yj(ip-1))/h(ip))
        !enddo
        !t2=0.d0
        !do ip=1,np-1
        !    t2=t2+ainv(2,ip)*((yj(ip+1)-yj(ip))/h(ip+1)-(yj(ip)-yj(ip-1))/h(ip))
        !enddo

        !t1=ainv(1,1)*((yj(2)-yj(1))/h(2)-(yj(1)-yj(0))/h(1))+ainv(1,2)*((yj(3)-yj(2))/h(3)-(yj(2)-yj(1))/h(2))
        !t2=ainv(2,1)*((yj(2)-yj(1))/h(2)-(yj(1)-yj(0))/h(1))+ainv(2,2)*((yj(3)-yj(2))/h(3)-(yj(2)-yj(1))/h(2))
        !write(64,'(4e24.15)') c(1:2),t1,t2
        !write(64,'(2e24.15)') e1(1),e1(2)
        !write(64,'(e24.15)') e2(1)
        !write(64,'(2e24.15)') ainv(1,1),ainv(1,2)
        !write(64,'(2e24.15)') ainv(2,1),ainv(2,2)
        !write(64,*) 
        !write(*,*) 'ALIREZA ',ainv(1,1)*((yj(2)-yj(1))/h(2)-(yj(1)-yj(0))/h(1)),c(1)
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(.not. move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) cycle
            yi(0:np)=x(i,0:np)
            !call prepcd1cd2(np,h,mp,yi,yj,cd1,cd2,ainv)
            call prepcd3cd4(np,h,mp,ainv,i,j,yi,yj,cd1,cd2)
            !write(*,*) 'ALI ',mp,cd1(1),cd2(1)
            !if(i==1 .and. j==1) then
            !    cd1(1)=0.d0 !-0.125d0
            !    cd2(1)=0.d0
            !else
            !    cd1(1)=0.d0
            !    cd2(1)=0.d0
            !endif
            !cd3(1:np-1)=cd3(1:np-1)+cd1(1:np-1)
            !cd4(1:np-1)=cd4(1:np-1)+cd2(1:np-1)
            call qdq(np,s,mp,tmax,c,h,i,j,yi,yj,cd1,cd2,ddt)
            dd(j,i,1:np-1)=ddt(1:np-1)
        enddo
    enddo
    deallocate(cd1,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating cd1.'
    deallocate(cd2,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating cd2.'
    !deallocate(cd3,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating cd3.'
    !deallocate(cd4,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating cd4.'
    deallocate(c,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating c.'
    deallocate(ddt,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ddt.'
    deallocate(yi,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating yi.'
    deallocate(yj,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating yj.'
    deallocate(ainv,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ainv.'
end subroutine prepdd


subroutine prepcd3cd4(np,h,mp,ainv,i,j,yi,yj,cd1,cd2)
    implicit none
    integer::np,mp,istat,i,j,ip,jp
    real(kind=8)::h(np),yi(0:np),yj(0:np),cd1(np-1),cd2(np-1),ainv(np-1,np-1)
    real(kind=8)::t1,t2,t3,t4,t5,t6,t7,t8,t9,tt1,tt2,delta
    real(kind=8), allocatable::ainvd(:,:)
    real(kind=8)::hip,hipp1,yip,yipp1,yipm1
    real(kind=8)::ainvdmpip,ainvdmpipp1,ainvdmpipm1,ainvdmpm1ip,ainvdmpm1ipp1,ainvdmpm1ipm1
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    allocate(ainvd(0:np,0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ainvd.'
    call dmemocc((np+1)*(np+1),(np+1)*(np+1+ndeb1),ainvd,'ainvd')
    ainvd(0:np,0:np)=0.d0
    do jp=1,np-1
        do ip=1,np-1
            ainvd(ip,jp)=ainv(ip,jp)
        enddo
    enddo
    do ip=1,np-1
        t1=delta(i,j)/h(ip)
        t2=delta(i,j)/h(ip+1)
        t3=(yj(ip)-yj(ip-1))*(yi(ip)-yi(ip-1))/h(ip)**3
        t4=(yj(ip+1)-yj(ip))*(yi(ip+1)-yi(ip))/h(ip+1)**3
        !write(99,'(a,6e20.10)') 'RR ',t1,t2,t3,t4,ainvd(mp,ip),-ainvd(mp,ip)*(t1+t2-t3-t4)
        cd1(ip)=ainvd(mp,ip-1)*(t1-t3)*(1.d0-delta(ip,1))-ainvd(mp,ip)*(t1+t2-t3-t4)+ &
            ainvd(mp,ip+1)*(t2-t4)*(1.d0-delta(ip,np-1))
        cd2(ip)=ainvd(mp-1,ip-1)*(t1-t3)*(1.d0-delta(ip,1))-ainvd(mp-1,ip)*(t1+t2-t3-t4)+ &
            ainvd(mp-1,ip+1)*(t2-t4)*(1.d0-delta(ip,np-1))
    enddo
    !return
    do ip=1,np-1
        hip=h(ip)
        hipp1=h(ip+1)
        yip=yi(ip)
        yipp1=yi(ip+1)
        yipm1=yi(ip-1)
        ainvdmpip=ainvd(mp,ip)
        ainvdmpipp1=ainvd(mp,ip+1)
        ainvdmpipm1=ainvd(mp,ip-1)
        ainvdmpm1ip=ainvd(mp-1,ip)
        ainvdmpm1ipp1=ainvd(mp-1,ip+1)
        ainvdmpm1ipm1=ainvd(mp-1,ip-1)
        tt1=0.d0
        tt2=0.d0
        do jp=1,np-1
            t1=ainvdmpip  *1.d0*(yip-yipm1)/hip  *ainvd(ip-1,jp)  
            t2=ainvdmpipp1*1.d0*(yipp1-yip)/hipp1*ainvd(ip,jp)    
            t3=ainvdmpip  *2.d0*(yip-yipm1)/hip  *ainvd(ip,jp)   
            t4=ainvdmpipp1*2.d0*(yipp1-yip)/hipp1*ainvd(ip+1,jp)  
            t5=ainvdmpipm1*2.d0*(yip-yipm1)/hip  *ainvd(ip-1,jp)  
            t6=ainvdmpip  *2.d0*(yipp1-yip)/hipp1*ainvd(ip,jp)    
            t7=ainvdmpipm1*1.d0*(yip-yipm1)/hip  *ainvd(ip,jp)    
            t8=ainvdmpip  *1.d0*(yipp1-yip)/hipp1*ainvd(ip+1,jp) 
            t9=((yj(jp+1)-yj(jp))/h(jp+1)-(yj(jp)-yj(jp-1))/h(jp))
            tt1=tt1+(-t1+t2-t3+t4-t5+t6-t7+t8)*t9
            t1=ainvdmpm1ip  *1.d0*(yip-yipm1)/hip  *ainvd(ip-1,jp)  
            t2=ainvdmpm1ipp1*1.d0*(yipp1-yip)/hipp1*ainvd(ip,jp)    
            t3=ainvdmpm1ip  *2.d0*(yip-yipm1)/hip  *ainvd(ip,jp)    
            t4=ainvdmpm1ipp1*2.d0*(yipp1-yip)/hipp1*ainvd(ip+1,jp)  
            t5=ainvdmpm1ipm1*2.d0*(yip-yipm1)/hip  *ainvd(ip-1,jp)  
            t6=ainvdmpm1ip  *2.d0*(yipp1-yip)/hipp1*ainvd(ip,jp)    
            t7=ainvdmpm1ipm1*1.d0*(yip-yipm1)/hip  *ainvd(ip,jp)    
            t8=ainvdmpm1ip  *1.d0*(yipp1-yip)/hipp1*ainvd(ip+1,jp) 
            t9=((yj(jp+1)-yj(jp))/h(jp+1)-(yj(jp)-yj(jp-1))/h(jp))
            tt2=tt2+(-t1+t2-t3+t4-t5+t6-t7+t8)*t9
        enddo
        cd1(ip)=cd1(ip)+tt1
        cd2(ip)=cd2(ip)+tt2
        !write(55,*) mp,cd1(ip),cd2(ip)
    enddo
    deallocate(ainvd,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ainvd.'
end subroutine prepcd3cd4


subroutine prepcd1cd2(np,h,mp,yi,yj,cd1,cd2,ainv)
    implicit none
    integer::np,mp,istat,ip,jp
    real(kind=8)::h(np),yi(0:np),yj(0:np),cd1(np-1),cd2(np-1),ainv(np-1,np-1),t1,t2,t3,t4
    real(kind=8), allocatable::ainvd(:,:)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    allocate(ainvd(0:np,0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating ainvd.'
    call dmemocc((np+1)*(np+1),(np+1)*(np+1+ndeb1),ainvd,'ainvd')
    ainvd(0:np,0:np)=0.d0
    do jp=1,np-1
        do ip=1,np-1
            ainvd(ip,jp)=ainv(ip,jp)
        enddo
    enddo
    do ip=1,np-1
        t1=1.d0/h(ip)
        t2=1.d0/h(ip+1)
        t3=2.d0*(yj(ip)-yj(ip-1))*(yi(ip)-yi(ip-1))/h(ip)**3
        t4=2.d0*(yj(ip+1)-yj(ip))*(yi(ip+1)-yi(ip))/h(ip+1)**3
        cd1(ip)=ainvd(mp,ip-1)*(t1-t3)-ainvd(mp,ip)*(t1+t2-t3-t4)+ainvd(mp,ip+1)*(t2-t4) 
        cd2(ip)=ainvd(mp-1,ip-1)*(t1-t3)-ainvd(mp-1,ip)*(t1+t2-t3-t4)+ainvd(mp-1,ip+1)*(t2-t4) 
    enddo
    deallocate(ainvd,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating ainvd.'

    !if(np==2) then
    !    if(mp==1) then
    !        cd2(1)=0.d0
    !    else
    !        cd2(1)=-ainv(mp-1,1)*(1.d0/h(2)+1.d0/h(1))
    !    endif
    !    if(mp==np) then
    !        cd1(1)=0.d0
    !    else
    !        cd1(1)=-ainv(mp,1)*(1.d0/h(2)+1.d0/h(1))
    !    endif
    !else
    !    if(mp==1) then
    !        cd2(1:np-1)=0.d0
    !    else
    !        cd2(1)=-ainv(mp-1,1)*(1.d0/h(2)+1.d0/h(1))+ainv(mp-1,2)/h(2)
    !        do ip=2,np-2
    !            cd2(ip)=ainv(mp-1,ip-1)/h(ip)-ainv(mp-1,ip)*(1.d0/h(ip+1)+ &
    !                1.d0/h(ip))+ainv(mp-1,ip+1)/h(ip+1)
    !        enddo
    !        cd2(np-1)=ainv(np-1,np-2)/h(np-1)-ainv(np-1,np-1)*(1.d0/h(np)+1.d0/h(np-1))
    !    endif
    !    if(mp==np) then
    !        cd1(1:np-1)=0.d0
    !    else
    !        cd1(1)=-ainv(mp,1)*(1.d0/h(2)+1.d0/h(1))+ainv(mp,2)/h(2)
    !        do ip=2,np-2
    !            cd1(ip)=ainv(mp,ip-1)/h(ip)-ainv(mp,ip)*(1.d0/h(ip+1)+1.d0/h(ip))+ &
    !                ainv(mp,ip+1)/h(ip+1)
    !        enddo
    !        cd1(np-1)=ainv(mp,np-2)/h(np-1)-ainv(mp,np-1)*(1.d0/h(np)+1.d0/h(np-1))
    !    endif
    !endif
end subroutine prepcd1cd2


subroutine func(tt,epot,ett,n,np,x,pnow,mp,xt,ft,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
    use module_base
    use module_interfaces
    use module_types
    use modulesplinedsaddle, only:parametersplinedsaddle
    implicit none
    integer, intent(in) :: nproc,iproc
    type(atoms_data), intent(inout) :: atoms
    type(input_variables), intent(inout) :: ll_inputs
    type(restart_objects), intent(inout) :: rst
    integer, intent(inout) :: ncount_bigdft
    integer::n,np,mp,i,istat,infocode
    type(parametersplinedsaddle)::pnow
    real(kind=8)::tt,ett,x(n,0:np),epot,xt(n),ft(n),t1,time1,time2
    real(kind=8), allocatable::tang(:)
    logical::move_this_coordinate
    integer::ixyz,iat
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    type(run_objects) :: runObj
    type(DFT_global_output) :: outs
    allocate(tang(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating tang.'
    call dmemocc(n,n+ndeb1,tang,'tang')
    do i=1,n
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
            pnow%y(0:np)=x(i,0:np)
            call inter_cubic(np,pnow%y,pnow%h,pnow%e1,pnow%e2,pnow%c)
            call ffdfdd_cubic(np,pnow%y,pnow%s,mp,pnow%h(mp),tt,pnow%c,xt(i),tang(i),t1)
        else
            xt(i)=x(i,0)
            tang(i)=0.d0
        endif
    enddo
    !do i=nr+1,n
    !    xt(i)=x(i,0)
    !    tang(i)=0.d0
    !enddo
    !call calenergyforces(n,xt,epot,ft)
    call run_objects_nullify(runObj)
    call run_objects_associate(runObj, ll_inputs, atoms, rst, xt(1))
    call init_global_output(outs, atoms%astruct%nat)
    call cpu_time(time1)
    call call_bigdft(runObj,outs,nproc,iproc,infocode)
    call cpu_time(time2)
    call run_objects_free_container(runObj)
    epot = outs%energy
    call deallocate_global_output(outs, ft(1))
    ncount_bigdft=ncount_bigdft+1
    pnow%ncount_ll=pnow%ncount_ll+1
    pnow%time_ll=pnow%time_ll+(time2-time1)
    call atomic_dot(atoms,ft,tang,ett)
    !ett=mydot(n,ft,tang)
    !write(*,'(a20,2f24.15,e24.15)') 'inside: tt,epot,ett',tt,epot,ett
    deallocate(tang,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating tang.'
end subroutine func


subroutine equalarclengthparametrization(atoms,n,np,x,s,h)
    use module_types
    implicit none
    type(atoms_data), intent(in) :: atoms
    integer::n,np
    real(kind=8)::x(n,0:np),s(0:np),h(np),tt
    integer::i,ip,ixyz,iat
    logical::move_this_coordinate
    s(0)=0.d0
    do ip=1,np
        tt=0.d0
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                tt=tt+(x(i,ip)-x(i,ip-1))**2
            endif
        enddo
        tt=sqrt(tt)
        h(ip)=tt;s(ip)=s(ip-1)+tt
    enddo
    !snpinv=1.d0/s(np);h(1:np)=h(1:np)*snpinv;s(1:np)=s(1:np)*snpinv
end subroutine equalarclengthparametrization


subroutine factor_cubic(np,h,e1,e2)
    implicit none
    integer::np
    integer::ip,info
    real(kind=8)::h(np),e1(np-1),e2(np-2) 
    do ip=1,np-2;e1(ip)=2.d0*(h(ip+1)+h(ip));e2(ip)=h(ip+1);enddo
    e1(np-1)=2.d0*(h(np)+h(np-1)) !+h(np)
    !e1(1)=e1(1)+h(1)
    call dpttrf(np-1,e1,e2,info)
    if(info/=0) write(*,*) 'ERROR: factorization failed: info,np',info,np
end subroutine factor_cubic


subroutine inter_cubic(np,y,h,e1,e2,c)
    implicit none
    integer::np
    integer::ip,info
    real(kind=8)::y(0:np),h(np),e1(np-1),e2(np-2),c(0:np)
    !real(kind=8)::tt,dt,b(4),ipiv(4),hi,bt0,btn,p0,p1,p2,p3
    do ip=1,np-1;c(ip)=(y(ip+1)-y(ip))/h(ip+1)-(y(ip)-y(ip-1))/h(ip);enddo
    !write(*,*) c(1:np-1)
    call dpttrs(np-1,1,e1,e2,c(1),np-1,info)
    !write(*,*) e1(1:np-1)
    !stop
    if(info/=0) write(*,*) 'ERROR: solution of dpttrs failed: info',info
    c(0)=0.d0;c(np)=0.d0
end subroutine inter_cubic


subroutine qdq(np,s,mp,tmax,c,h,i,j,yi,yj,cd1,cd2,dd)
    implicit none
    integer::np,mp,i,j,ip,istat
    real(kind=8)::s(0:np),tmax,c(0:np),h(np),yi(0:np),yj(0:np),cd1(np-1),cd2(np-1),dd(np-1)
    real(kind=8)::p0,p1,p2,p3,delta,t1,t2,t3,t4,t5,t6,t7
    real(kind=8), allocatable::sd1(:)
    real(kind=8), allocatable::sd2(:)
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    allocate(sd1(np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating sd1.'
    call dmemocc(np-1,np-1+ndeb1,sd1,'sd1')
    allocate(sd2(np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating sd2.'
    call dmemocc(np-1,np-1+ndeb1,sd2,'sd2')
    call calsd1sd2(np,mp,yi,h,sd1,sd2)
    do ip=1,np-1
        t3=(yi(mp)-yi(mp-1))*(delta(mp,ip)-delta(mp-1,ip))/h(mp)
        p3=(cd1(ip)-cd2(ip))/h(mp)-(c(mp)-c(mp-1))*t3/h(mp)**2
        p2=3.d0*((cd2(ip)*s(mp)+c(mp-1)*sd1(ip)-cd1(ip)*s(mp-1)-c(mp)*sd2(ip))/h(mp)-&
            (c(mp-1)*s(mp)-c(mp)*s(mp-1))*t3/h(mp)**2)

        t1=3.d0*((cd1(ip)*s(mp-1)**2+2.d0*c(mp)*s(mp-1)*sd2(ip)-cd2(ip)*s(mp)**2-&
            2.d0*c(mp-1)*s(mp)*sd1(ip))/h(mp)-(c(mp)*s(mp-1)**2-c(mp-1)*s(mp)**2)*t3/h(mp)**2)
        t4=delta(i,j)*delta(mp,ip)/h(mp)-yj(mp)*t3/h(mp)**2-c(mp)*t3-h(mp)*cd1(ip)
        t5=delta(i,j)*delta(mp-1,ip)/h(mp)-yj(mp-1)*t3/h(mp)**2-c(mp-1)*t3-h(mp)*cd2(ip)
        p1=t1+t4-t5
        !p1=(t1+delta(mp,ip)-delta(mp-1,ip))/h(mp)+h(mp)*(cd2(ip)-cd1(ip))
        t2=(cd2(ip)*s(mp)**3+3.d0*c(mp-1)*s(mp)**2*sd1(ip)-cd1(ip)*s(mp-1)**3-3.d0*c(mp)*s(mp-1)**2*sd2(ip))/h(mp)-&
            (c(mp-1)*s(mp)**3-c(mp)*s(mp-1)**3)*t3/h(mp)**2
        t6=-t4*s(mp-1)-(yj(mp)/h(mp)-h(mp)*c(mp))*sd2(ip)
        t7=t5*s(mp)+(yj(mp-1)/h(mp)-h(mp)*c(mp-1))*sd1(ip)
        p0=t2+t6+t7
        !p0=t2/h(mp)+h(mp)*(cd1(ip)*s(mp-1)-cd2(ip)*s(mp))
        dd(ip)=((p3*tmax+p2)*tmax+p1)*tmax+p0
    enddo
    deallocate(sd1,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating sd1.'
    deallocate(sd2,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating sd2.'
end subroutine qdq


subroutine calsd1sd2(np,mp,yi,h,sd1,sd2)
    implicit none
    integer::np,mp,ip
    real(kind=8)::yi(0:np),h(np),sd1(np-1),sd2(np-1)
    do ip=1,mp-1
        sd1(ip)=(yi(ip)-yi(ip-1))/h(ip)-(yi(ip+1)-yi(ip))/h(ip+1)
    enddo
    if(mp<np) sd1(mp)=(yi(mp)-yi(mp-1))/h(mp)
    if(mp<np-1) sd1(mp+1:np-1)=0.d0
    do ip=1,mp-2
        sd2(ip)=(yi(ip)-yi(ip-1))/h(ip)-(yi(ip+1)-yi(ip))/h(ip+1)
    enddo
    if(mp>1) sd2(mp-1)=(yi(mp-1)-yi(mp-2))/h(mp-1)
    if(mp<np) sd2(mp:np-1)=0.d0
end subroutine calsd1sd2


function delta(i,j)
    implicit none
    integer::i,j
    real(kind=8)::delta
    if(i==j) then
        delta=1.d0
    else
        delta=0.d0
    endif
end function delta


subroutine ffdfdd_cubic(np,y,s,mp,hmp,t,c,f,fd,fdd)
    implicit none
    integer::np,mp
    real(kind=8)::y(0:np),s(0:np),hmp,c(0:np),t,p0,p1,p2,p3,f,fd,fdd
    if(mp<1 .or. mp>np) stop 'ERROR: invalid mp in cubic evaluation'
    !hmp=s(mp)-s(mp-1)
    p3=(c(mp)-c(mp-1))/hmp
    p2=3.d0*(c(mp-1)*s(mp)-c(mp)*s(mp-1))/hmp
    p1=(3.d0*(c(mp)*s(mp-1)**2-c(mp-1)*s(mp)**2)+y(mp)-y(mp-1))/hmp+hmp*(c(mp-1)-c(mp))
    p0=(c(mp-1)*s(mp)**3-c(mp)*s(mp-1)**3-y(mp)*s(mp-1)+y(mp-1)*s(mp))/hmp+ &
        hmp*(c(mp)*s(mp-1)-c(mp-1)*s(mp))
    f=((p3*t+p2)*t+p1)*t+p0
    fd=(3.d0*p3*t+2.d0*p2)*t+p1
    fdd=6.d0*p3*t+2.d0*p2
end subroutine ffdfdd_cubic


subroutine caltangentupwind(n,np,x,ex,tang)
    implicit none
    integer::n,np,ip
    real(kind=8)::x(n,0:np),ex(0:np),tang(n,0:np)
    real(kind=8)::tmp1,tmp2,e1,e2
    do ip=1,np-1
        if(ex(ip+1)>ex(ip) .and. ex(ip)>ex(ip-1)) then
            tang(1:n,ip)=x(1:n,ip+1)-x(1:n,ip)
        elseif(ex(ip+1)<ex(ip) .and. ex(ip)<ex(ip-1)) then
            tang(1:n,ip)=x(1:n,ip)-x(1:n,ip-1)
        else
            tmp1=abs(ex(ip-1)-ex(ip))
            tmp2=abs(ex(ip)-ex(ip+1))
            e1=max(tmp1,tmp2)
            e2=min(tmp1,tmp2)
            if(ex(ip+1)>ex(ip-1)) then
                tang(1:n,ip)=e1*(x(1:n,ip+1)-x(1:n,ip))+e2*(x(1:n,ip)-x(1:n,ip-1))
            else
                tang(1:n,ip)=e2*(x(1:n,ip+1)-x(1:n,ip))+e1*(x(1:n,ip)-x(1:n,ip-1))
            endif
        endif
    enddo
    tang(1:n,0)=tang(1:n,1)
    tang(1:n,np)=tang(1:n,np-1)
    do ip=0,np
        call normalizevector2(n,tang(1,ip))
    enddo
end subroutine caltangentupwind


subroutine normalizevector2(n,v)
    implicit none
    integer::n,i
    real(kind=8)::v(n),vnrm
    vnrm=0.d0
    do i=1,n
        vnrm=vnrm+v(i)**2
    enddo
    vnrm=sqrt(vnrm);v(1:n)=v(1:n)/vnrm
end subroutine normalizevector2


subroutine initminimize(parmin)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    type(parameterminimization_sp)::parmin
    integer::istat
    character(2)::tapp1,tapp2
    character(4)::tapp3
    integer, parameter::ndeb1=0
    !integer, parameter::ndeb2=0
    tapp1(1:2)=parmin%approach(1:2)
    if(len(trim(parmin%approach))==4) tapp2(1:2)=parmin%approach(3:4)
    if(len(trim(parmin%approach))==6) tapp3(1:4)=parmin%approach(3:6)
    !write(*,*) tapp1
    !write(*,*) tapp2
    !write(*,*) tapp3
    !write(*,*) len(trim(parmin%approach))
    if(parmin%fmaxtol<0.d0) stop 'ERROR: fmaxtol<0, maybe it is not set.'
    if(tapp1=='SD') then
        !if(parmin%alphax<0.d0) stop 'ERROR: alphax<0, maybe it is not set.'
        !parmin%alphamin=1.d-2*parmin%alphax
        !parmin%alphamax=2.d0*parmin%alphax
        !parmin%fnrmtolsatur=1.d0 !parmin%fmaxtol**0.1d0
        !parmin%nitsd=10000
        !parmin%nsatur=2
    endif
    if(tapp1=='SD' .or. tapp1=='CG' .or. tapp2=='CG') then
        if(parmin%anoise<0.d0) parmin%anoise=1.d-12 !epsilon(parmin%anoise)
    endif
    if(tapp3=='DIIS') then
        parmin%idsx=20
        allocate(parmin%a(parmin%idsx+1,parmin%idsx+1,3+ndeb1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating parmin%a.'
        call dmemocc((parmin%idsx+1)*(parmin%idsx+1)*(3),(parmin%idsx+1)*(parmin%idsx+1)*(3+ndeb1),parmin%a,'parmin%a')
        allocate(parmin%b(parmin%idsx+1+ndeb1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating parmin%b.'
        call dmemocc(parmin%idsx+1,parmin%idsx+1+ndeb1,parmin%b,'parmin%b')
        allocate(parmin%ipiv(parmin%idsx+1+ndeb1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating parmin%ipiv.'
        call imemocc(parmin%idsx+1,parmin%idsx+1+ndeb1,parmin%ipiv,'parmin%ipiv')
    endif
end subroutine initminimize


subroutine finalminimize(parmin)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    type(parameterminimization_sp)::parmin
    integer::istat
    !n(c) character(2)::tapp1,tapp2
    character(4)::tapp3
    !n(c) tapp1(1:2)=parmin%approach(1:2)
    !n(c) if(len(trim(parmin%approach))==4) tapp2(1:2)=parmin%approach(3:4)
    if(len(trim(parmin%approach))==6) tapp3(1:4)=parmin%approach(3:6)
    if(tapp3=='DIIS') then
        parmin%idsx=10
        deallocate(parmin%a,stat=istat)
        if(istat/=0) stop 'ERROR: failure deallocating parmin%a.'
        deallocate(parmin%b,stat=istat)
        if(istat/=0) stop 'ERROR: failure deallocating parmin%b.'
        deallocate(parmin%ipiv,stat=istat)
        if(istat/=0) stop 'ERROR: failure deallocating parmin%ipiv.'
    endif
end subroutine finalminimize


subroutine checkconvergence(parmin,fspmax)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    !integer::iproc
    real(kind=8)::fspmax
    type(parameterminimization_sp)::parmin
    if(fspmax<parmin%fmaxtol) then
        parmin%converged=.true.
    endif
end subroutine checkconvergence


!routine for moving atomic positions, takes into account the 
!frozen atoms and the size of the cell
!synopsis: xyzo=xyzi
!all the shift are inserted into the box if there are periodic directions
!if the atom are frozen they are not moved
subroutine atomic_copycoord(atoms,xyzi,xyzo)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: xyzi
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: xyzo
  !local variables
  integer :: iat,ix,iy,iz
  logical::move_this_coordinate

  do iat=1,atoms%astruct%nat
     ix=0;iy=0;iz=0
     if(move_this_coordinate(atoms%astruct%ifrztyp(iat),1)) ix=1
     if(move_this_coordinate(atoms%astruct%ifrztyp(iat),2)) iy=1
     if(move_this_coordinate(atoms%astruct%ifrztyp(iat),3)) iz=1
     if (atoms%astruct%geocode == 'P') then
        if(ix==1) xyzo(1,iat)=modulo(xyzi(1,iat),atoms%astruct%cell_dim(1))
        if(iy==1) xyzo(2,iat)=modulo(xyzi(2,iat),atoms%astruct%cell_dim(2))
        if(iz==1) xyzo(3,iat)=modulo(xyzi(3,iat),atoms%astruct%cell_dim(3))
     else if (atoms%astruct%geocode == 'S') then
        if(ix==1) xyzo(1,iat)=modulo(xyzi(1,iat),atoms%astruct%cell_dim(1))
        if(iy==1) xyzo(2,iat)=xyzi(2,iat)
        if(iz==1) xyzo(3,iat)=modulo(xyzi(3,iat),atoms%astruct%cell_dim(3))
     else
        if(ix==1) xyzo(1,iat)=xyzi(1,iat)
        if(iy==1) xyzo(2,iat)=xyzi(2,iat)
        if(iz==1) xyzo(3,iat)=xyzi(3,iat)
     end if
  end do
end subroutine atomic_copycoord


subroutine initsdminimum(n,nr,x,parmin,nwork,work)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::n,nr,nwork
    !integer::ip,ii
    !real(kind=8)::x(n*(np-1)),work(nwork)
    real(kind=8)::x(n),work(nwork)
    type(parameterminimization_sp)::parmin
    if(parmin%alphax<0.d0) stop 'ERROR: alphax<0, maybe it is not set.'
    if(parmin%alphamin<0.d0) parmin%alphamin=1.d-1*parmin%alphax
    if(parmin%alphamax<0.d0) parmin%alphamax=2.d0*parmin%alphax
    if(parmin%fnrmtolsatur<0.d0) parmin%fnrmtolsatur=1.d0 !parmin%fmaxtol**0.1d0
    if(parmin%nsatur<0) parmin%nsatur=2
    if(parmin%nitsd<0) parmin%nitsd=10000
    if(parmin%anoise<0.d0) parmin%anoise=1.d-12 !epsilon(parmin%anoise)
    parmin%iflag=1
    parmin%itsd=0
    if(parmin%alpha<0.d0) parmin%alpha=1.d0*parmin%alphax
    parmin%sdsaturated=.false.
    parmin%epotitm1=1.d50;parmin%epotitm2=1.d50
    parmin%fnrmitm1=1.d50;parmin%fnrmitm2=1.d50
    parmin%care=.true.
    parmin%isatur=0
    !do ip=1,np-1
    !    ii=1+(ip-1)*n
    !    call atomic_copycoord(atoms,work(ii),x(ii))
    !enddo
    !work(1:n*(np-1))=x(1:n*(np-1))
    !work(n*(np-1)+1:2*n*(np-1))=0.d0
    work(1:nr)=x(1:nr)
    work(nr+1:2*nr)=0.d0
    !write(*,'(a,i,e)') 'nsatur,fnrmtolsatur',parmin%nsatur,parmin%fnrmtolsatur
end subroutine initsdminimum


subroutine fire_splsad(iproc,nr,x,epot,f,work,parmin)
    !use minimization, only:parameterminimization
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::iproc,nr
    real(kind=8)::x(nr),epot,f(nr),de,DDOT,fnrm,fmax,vnrm,dt,p,calnorm
    real(kind=8)::tt,calmaxforcecomponent,vnrmmax
    real(kind=8)::work(3*nr) !1:nr velocities, nr+1:2*nr previous force
    real(kind=8), save::epotold,alpha
    integer, save::ndown
    type(parameterminimization_sp)::parmin
    if(parmin%iflag==0) then
        if(parmin%dt<-0.d0) then
            if(iproc==0) &
            write(*,*) 'ERROR: time step in FIRE method must be set by user'
            return
        endif
        if(parmin%dtmax<0.d0) parmin%dtmax=30.d0*parmin%dt
        if(parmin%finc<0.d0) parmin%finc=1.20d0
        if(parmin%fdec<0.d0) parmin%fdec=0.5d0
        if(parmin%falpha<0.d0) parmin%falpha=0.97d0  !0.99d0
        if(parmin%alphastart<0.d0) parmin%alphastart=0.5d0  !0.1d0  !0.05d0
        if(parmin%alphax<0.d0) then
            if(iproc==0) &
            write(*,*) 'ERROR: alphax in FIRE method must be set by user'
            return
        endif
        if(parmin%ndowntol<0) parmin%ndowntol=1
        !work(1:nr)=0.d0
        work(1:nr)=2.d-1*parmin%alphax*f(1:nr)
        !work(1:nr)=0.5d0*parmin%alphax*f(1:nr)
        !work(1:nr)=50.d0*parmin%alphax*f(1:nr)
        !do i=1,nr
        !    call random_number(tt)
        !    work(i)=2.d0*parmin%alphax*(f(i)+5.d-1*tt)
        !enddo
        work(nr+1:2*nr)=-f(1:nr)
        epotold=epot
        alpha=parmin%alphastart
        ndown=0
        parmin%iflag=1
    endif
    dt=parmin%dt
    if(ndown/=0) then
        work(1:nr)=work(1:nr)+0.5d0*dt*(f(1:nr)+work(nr+1:2*nr))
    endif
    p=DDOT(nr,f,1,work,1)
    vnrm=calnorm(nr,work);fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
    de=epot-epotold
    if(iproc==0) &
    write(*,'(a10,i4,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4)') &
        'FIREMIN   ',parmin%itfire,epot,de,fnrm,fmax,vnrm,dt,alpha,ndown,p
    !if(iproc==0) &
    !write(parmin%ifile,'(a10,i4,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4)') &
    !    'FIREMIN   ',parmin%itfire,epot,de,fnrm,fmax,vnrm,dt,alpha,ndown,p
    if(fmax<parmin%fmaxtol) then
        parmin%converged=.true.
        parmin%iflag=0
        if(iproc==0) &
        write(*,'(a,i4,es23.15,2es12.5)') &
            'FIRE FINISHED: itfire,epot,fnrm,fmax ',parmin%itfire,epot,fnrm,fmax
        return
    endif
    parmin%itfire=parmin%itfire+1
    epotold=epot
    work(2*nr+1:3*nr)=x(1:nr)
    x(1:nr)=x(1:nr)+dt*work(1:nr)+0.5d0*dt**2*f(1:nr)
    !f(1:nr)=f(1:nr)/fnrm
    !tt=min(alpha*vnrm/fnrm,3.d0*parmin%alphax)
    !tt=min(alpha*vnrm/fnrm,1.d-2*parmin%alphax)
    tt=min(alpha*vnrm/fnrm,5.d-1*parmin%alphax)
    work(1:nr)=(1.d0-alpha)*work(1:nr)+tt*f(1:nr) !min(alpha*vnrm,3.d0*parmin%alphax*fnrm)
    !-------------------------------------------------------
    tt=calnorm(nr,work)
    if(iproc==0) write(*,'(a,2es19.10)') 'fort56 ',tt,fnrm
    vnrmmax=20.d0*fnrm
    if(tt>vnrmmax) then
        work(1:nr)=work(1:nr)*vnrmmax/tt
    endif
    !-------------------------------------------------------
    if(.not. p<0.d0) then
    !if(p>-1.d-6) then
    !if(p>-1.d-6 .or. parmin%itfire<=40) then
        if(ndown>parmin%ndowntol) then
            parmin%dt=min(parmin%finc*parmin%dt,parmin%dtmax)
            alpha=max(parmin%falpha*alpha,2.d-1)
            !alpha=parmin%falpha*alpha
        endif
        ndown=ndown+1
    else
        parmin%dt=parmin%fdec*parmin%dt
        x(1:nr)=work(2*nr+1:3*nr)
        f(1:nr)=work(nr+1:2*nr)
        work(1:nr)=0.d0
        !work(1:nr)=0.05d0*parmin%alphax*f(1:nr)
        alpha=parmin%alphastart
        ndown=0
    endif
    work(nr+1:2*nr)=f(1:nr)
  end subroutine fire_splsad


subroutine sdminimum(iproc,n,nr,x,f,epot,parmin,nwork,work)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::iproc,n,nr,nwork
    real(kind=8)::x(n),f(n),epot,work(nwork),calmaxforcecomponent,calnorm,fmax,fnrm
    type(parameterminimization_sp)::parmin
    real(kind=8)::de1,de2,df1,df2
    logical::xmoved
    if(parmin%iflag==0) call initsdminimum(n,nr,x,parmin,nwork,work)
    fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
    de1=epot-parmin%epotitm1;de2=epot-2.d0*parmin%epotitm1+parmin%epotitm2
    df1=fnrm-parmin%fnrmitm1;df2=fnrm-2.d0*parmin%fnrmitm1+parmin%fnrmitm2
    if(parmin%itsd==0) de1=0.d0
    xmoved=.true.
    if(parmin%care .and. de1>parmin%anoise) then
        if(parmin%alpha<parmin%alphamin .and. parmin%itsd/=0) then
        if(iproc==0) then
        write(parmin%ifile,'(a)') 'alpha getting too small, do not care anymore if energy goes up'
        write(*,'(a)') 'alpha getting too small, do not care anymore if energy goes up'
        endif
            parmin%care=.false.
         else
            x(1:nr)=work(1:nr)
            f(1:nr)=work(nr+1:nr+nr)
            xmoved=.false.
        endif
    endif
    !write(*,'(5(1x,e11.4),1x,i3)') fnrm/parmin%fnrmitm1, de1,de2,df1,df2,isatur
    !if(parmin%care .and. parmin%itsd>5 .and. parmin%alpha==alphax .and. fnrm/parmin%fnrmitm1>0.8d0 &
    if(parmin%care .and. parmin%itsd>1 .and. fnrm/parmin%fnrmitm1>0.5d0 .and. de1>-0.1d0 &
        .and. fnrm<parmin%fnrmtolsatur .and. de1<parmin%anoise .and. df1<parmin%anoise &
        .and. de2>-2.d0*parmin%anoise .and. df2>-2.d0*parmin%anoise) then 
        parmin%isatur=parmin%isatur+1
    else
        parmin%isatur=0
    endif
    if(iproc==0) then
    write(parmin%ifile,'(a10,i4,e23.15,e11.3,2e12.5,e12.4,i5,l2)') 'SDMIN     ',parmin%itsd,epot,de1, &
        fnrm,fmax,0.5d0*(sign(1.d0,real(parmin%itsd-1,8))+1.d0)*parmin%alpha/parmin%alphax,parmin%isatur,xmoved
    write(*,'(a10,i4,e23.15,e11.3,2e12.5,e12.4,i5,l2)') 'SDMIN     ',parmin%itsd,epot,de1, &
        fnrm,fmax,0.5d0*(sign(1.d0,real(parmin%itsd-1,8))+1.d0)*parmin%alpha/parmin%alphax,parmin%isatur,xmoved
    endif
    if(parmin%converged) then
        parmin%sdminimum=.false.
        parmin%iflag=0
        if(iproc==0) then
            write(parmin%ifile,'(a,i4,e23.15,2e12.5)') &
                'SD FINISHED: itsd,epot,fnrm,fmax',parmin%itsd,epot,fnrm,fmax
            write(*,'(a,i4,e23.15,2e12.5)') &
                'SD FINISHED: itsd,epot,fnrm,fmax',parmin%itsd,epot,fnrm,fmax
        endif
        return
    endif
    if(parmin%isatur>=parmin%nsatur .and. .not. parmin%sdsaturated) then
        parmin%sdsaturated=.true.
        parmin%sdminimum=.false.
        if(iproc==0) then
            write(parmin%ifile,'(a,i4,e23.15,2e12.5)') &
                'SD SATURATED: itsd,epot,fnrm,fmax',parmin%itsd,epot,fnrm,fmax
            write(*,'(a,i4,e23.15,2e12.5)') &
                'SD SATURATED: itsd,epot,fnrm,fmax',parmin%itsd,epot,fnrm,fmax
        endif
        if(trim(parmin%approach)/='SD') then
            parmin%iflag=0
        endif
            return
    endif
    !if(parmin%isatur>parmin%nsatur) then
    !    parmin%sdsaturated=.true.
    !    parmin%iflag=0
    !    write(*,'(a,i4,e23.15,2e12.5)') &
    !        'SD SATURATED: itsd,epot,fnrm,fmax',parmin%itsd,epot,fnrm,fmax
    !    return
    !endif
    if(parmin%care .and. de1>parmin%anoise) then
        parmin%alpha=5.d-1*parmin%alpha
    endif
    if(xmoved) then
        parmin%epotitm2=parmin%epotitm1;parmin%epotitm1=epot
        parmin%fnrmitm2=parmin%fnrmitm1;parmin%fnrmitm1=fnrm
        !parmin%alpha=min(1.05d0*parmin%alpha,parmin%alphamax)
        parmin%alpha=min(1.20d0*parmin%alpha,parmin%alphamax)
        work(1:nr)=x(1:nr)
        work(nr+1:nr+nr)=f(1:nr)
    endif
    x(1:nr)=x(1:nr)+parmin%alpha*f(1:nr)
    if(.not. parmin%care .and. parmin%alpha>2.d0*parmin%alphamin) then
        parmin%care=.true.
    if(iproc==0) then
            write(parmin%ifile,'(a)') 'sdminimum starts to care whether energy goes up' 
            write(*,'(a)') 'sdminimum starts to care whether energy goes up' 
    endif
    endif
    parmin%itsd=parmin%itsd+1
    if(parmin%itsd>parmin%nitsd) then 
        parmin%iflag=-1
        x(1:nr)=work(1:nr)
        f(1:nr)=work(nr+1:nr+nr)
    if(iproc==0) then
            write(parmin%ifile,'(a,e23.15,e12.5)') 'SD: NO CONVERGENCE: fnrm,epot ',epot,fmax
            write(*,'(a,e23.15,e12.5)') 'SD: NO CONVERGENCE: fnrm,epot ',epot,fmax
    endif
    endif
end subroutine sdminimum


subroutine diisminimum(iproc,n,nr,x,epot,f,parmin,nwork,work)
    use minimization_sp, only:parameterminimization_sp
    implicit none
    integer::n,nr,nwork,i,info,id,jd,iproc
    real(kind=8)::x(n),f(n),epot,work(nwork),fnrm,dnrm2,ddot,fmax,calmaxforcecomponent
    type(parameterminimization_sp)::parmin
    character(28), parameter::frt1='(a10,i4,e23.15,e11.3,2e12.5)'
    if(parmin%iflag==0) then
        parmin%iflag=1;parmin%itdiis=0;parmin%epotitm1=epot
        parmin%emin=1.d100;parmin%fnrmlowest=1.d100;parmin%ld=0;parmin%nd=0
        !if(iproc==0) then
        !    write(parmin%ifile,*) 'iflag=0 as it must be.'
        !    write(*           ,*) 'iflag=0 as it must be.'
        !endif
    endif
    fnrm=dnrm2(nr,f,1)
    if(epot>parmin%emin+1.d-2*abs(parmin%emin) .or. fnrm>2.d0*parmin%fnrmlowest) then 
        if(iproc==0) then
            write(parmin%ifile,*) 'DIVERGENCE in DIIS, switch back to SD',parmin%itdiis
            write(*           ,*) 'DIVERGENCE in DIIS, switch back to SD',parmin%itdiis
        endif
        parmin%diisdivergence=.true.
        call dcopy(nr,work((3*parmin%idsx+2)*nr+1),1,x,1)
        !call sdminimum(0,n,n,x,fnrmtol,f,epot,sdconverged)
        !parmin%emin=1.d100;parmin%fnrmlowest=1.d100;parmin%ld=0;parmin%nd=0;parmin%epotitm1=epot
        parmin%sdsaturated=.false.
        parmin%iflag=0
        return
    endif
    parmin%nd=min(parmin%nd+1,parmin%idsx)
    parmin%ld=mod(parmin%ld,parmin%idsx)+1
    !call dcopy(n,x,1,xh(1,parmin%ld),1)
    call dcopy(nr,x,1,work(parmin%ld*nr+1),1)
    if(epot<parmin%emin) then 
        parmin%emin=epot;parmin%fnrmlowest=fnrm
        call dcopy(nr,x,1,work((3*parmin%idsx+2)*nr+1),1)
    endif
    call dcopy(nr,f,1,work((parmin%idsx+1)*nr+parmin%ld*nr+1),1)
    fmax=calmaxforcecomponent(nr,f)
    if(iproc==0) then
        write(parmin%ifile,frt1) 'DIISMIN   ',parmin%itdiis,epot,epot-parmin%epotitm1,fnrm,fmax
        write(*           ,frt1) 'DIISMIN   ',parmin%itdiis,epot,epot-parmin%epotitm1,fnrm,fmax
    endif
    if(parmin%converged) then
        parmin%iflag=0
        return
    endif
    !if(fmax<parmin%fmaxtol) then
    !    if(iproc==0) then
    !            write(parmin%ifile,'(a,i4,e23.15,2e12.5)') &
    !        'DIIS finished: ',parmin%itdiis,epot,fnrm,fmax
    !            write(*,'(a,i4,e23.15,2e12.5)') &
    !        'DIIS finished: ',parmin%itdiis,epot,fnrm,fmax
    !    endif
    !    parmin%iflag=0
    !    parmin%converged=.true.
    !    return
    !endif
    call dcopy(nr,f,1,work((2*parmin%idsx+2)*nr+(parmin%ld-1)*nr+1),1)
    !set up DIIS matrix (upper triangle)
    if(parmin%itdiis>parmin%idsx-1) then !shift left up matrix
        do i=1,parmin%idsx-1;parmin%a(1:i,i,1)=parmin%a(2:i+1,i+1,1);enddo
    endif
    !calculate new line, use b as work array for summation
    do id=1,parmin%nd
        jd=mod(parmin%ld+id-1,parmin%nd)+1
        parmin%a(id,parmin%nd,1)=ddot(nr,work((2*parmin%idsx+2)*nr+(parmin%ld-1)*nr+1),1, &
        work((2*parmin%idsx+2)*nr+(jd-1)*nr+1),1)
    enddo
    do i=1,parmin%nd;parmin%a(i,i:parmin%nd,2)=parmin%a(i,i:parmin%nd,1);enddo !copy to work array
    parmin%a(1:parmin%nd,parmin%nd+1,2)=1.d0;parmin%a(parmin%nd+1,parmin%nd+1,2)=0.d0 !prepare boundary elements
    parmin%b(1:parmin%nd)=0.d0;parmin%b(parmin%nd+1)=1.d0 !prepare right hand side
    if(parmin%itdiis>0) then !solve linear system:(LAPACK)
        call dsysv('U',parmin%nd+1,1,parmin%a(1,1,2),parmin%idsx+1,parmin%ipiv,parmin%b, &
            parmin%idsx+1,parmin%a(1,1,3),(parmin%idsx+1)**2,info)
        if(info/=0) then;write(parmin%ifile,*) 'ERROR: dsysv failed: info',info;stop;endif
    else
        parmin%b(1)=1.d0
    endif
    if(iproc==0) write(parmin%ifile,'(a,11(1pe11.2))')'DIIS weights',parmin%b(1:parmin%nd+1)
    !xh(1:nr,0)=0.d0;fh(1:nr,0)=0.d0 !new guess
    work(1:nr)=0.d0;work((parmin%idsx+1)*nr+1:(parmin%idsx+1)*nr+nr)=0.d0 !new guess
    do id=1,parmin%nd
        jd=mod(parmin%ld+id-1,parmin%nd)+1
        !xh(1:nr,0)=xh(1:nr,0)+b(id)*xh(1:nr,jd)
        work(1:nr)=work(1:nr)+parmin%b(id)*work(jd*nr+1:jd*nr+nr)  
        !fh(1:nr,0)=fh(1:nr,0)+b(id)*fh(1:nr,jd)
        work((parmin%idsx+1)*nr+1:(parmin%idsx+1)*nr+nr)=work((parmin%idsx+1)*nr+1:(parmin%idsx+1)*nr+nr)+&
            parmin%b(id)*work((parmin%idsx+1)*nr+jd*nr+1:(parmin%idsx+1)*nr+jd*nr+nr)
    enddo
    x(1:nr)=work(1:nr)+work((parmin%idsx+1)*nr+1:(parmin%idsx+1)*nr+nr)*parmin%alphax*2.d0
    parmin%epotitm1=epot
    parmin%itdiis=parmin%itdiis+1
end subroutine diisminimum


function mydot(n,v1,v2)
    implicit none
    integer::n,i
    real(kind=8)::v1(n),v2(n),mydot
    mydot=0.d0;do i=1,n;mydot=mydot+v1(i)*v2(i);enddo
end function mydot


function calmaxforcecomponent(n,v)
    implicit none
    integer::n,i
    real(kind=8)::v(n),calmaxforcecomponent
    calmaxforcecomponent=0.d0
    do i=1,n;calmaxforcecomponent=max(calmaxforcecomponent,abs(v(i)));enddo
end function calmaxforcecomponent


function calnorm(n,v)
    implicit none
    integer::n,i
    real(kind=8)::v(n),calnorm
    calnorm=0.d0
    do i=1,n;calnorm=calnorm+v(i)**2;enddo;calnorm=sqrt(calnorm)
end function calnorm


subroutine writepathway(n,np,x,filename,atoms)
    use module_types
    implicit none
    integer::n,np,jp,iat
    real(kind=8)::x(n,0:np),ed_tt,edd_tt,dtt,tt,xyz(3)
    integer::istat,ip,i
    type(atoms_data), intent(in) :: atoms
    character(len=10) :: name
    character(len=2) :: symbol
    !integer::istep
    !character(1)::fn
    character(*)::filename
    real(kind=8), allocatable::xt(:)
    real(kind=8), allocatable::s(:)
    real(kind=8), allocatable::h(:)
    real(kind=8), allocatable::e1(:)
    real(kind=8), allocatable::e2(:)
    real(kind=8), allocatable::y(:)
    real(kind=8), allocatable::c(:)
    logical::move_this_coordinate
    integer::ixyz
    integer, parameter::ndeb1=0 !n(c) ndeb2=0
    allocate(s(0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating s.'
    call dmemocc(np+1,np+1+ndeb1,s,'s')
    allocate(h(np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating h.'
    call dmemocc(np,np+ndeb1,h,'h')
    allocate(e1(np-1+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating e1.'
    call dmemocc(np-1,np-1+ndeb1,e1,'e1')
    allocate(e2(np-2+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating e2.'
    call dmemocc(np-2,np-2+ndeb1,e2,'e2')
    allocate(y(0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating y.'
    call dmemocc(np+1,np+1+ndeb1,y,'y')
    allocate(c(0:np+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating c.'
    call dmemocc(np+1,np+1+ndeb1,c,'c')
    allocate(xt(n+ndeb1),stat=istat);if(istat/=0) stop 'ERROR: failure allocating xt.'
    call dmemocc(n,n+ndeb1,xt,'xt')
    !-------------------------------------------------------------------------------------------
    call equalarclengthparametrization(atoms,n,np,x,s,h)
    call factor_cubic(np,h,e1,e2)
    !write(fn,'(i1.1)') iposout
    !filename='posout'//fn//'.xyz'
    open(unit=1388,file=filename,status='replace')
    ip=1;dtt=s(np)/100.d0*(1.d0-1.d-14)
    do jp=0,100
        tt=dtt*jp 
        call calindex(np,s,tt,ip,'writepathway')
        !write(*,*) 'jp,ip',jp,ip
        write(1388,*) n/3
        write(1388,*) jp,'   '
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                y(0:np)=x(i,0:np)
                call inter_cubic(np,y,h,e1,e2,c)
                call ffdfdd_cubic(np,y,s,ip,h(ip),tt,c,xt(i),ed_tt,edd_tt)
            else
                xt(i)=x(i,0)
            endif
            if(mod(i,3)==0) then
                iat=i/3
                name=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
                if (name(3:3)=='_') then
                   symbol=name(1:2)
                   !suffix=name(4:6)
                else if (name(2:2)=='_') then
                   symbol=name(1:1)
                   !suffix=name(3:5)
                else
                   symbol=name(1:2)
                   !suffix=' '
                end if
                if(atoms%astruct%units=='angstroemd0' .or. atoms%astruct%units=='angstroem') then
                    xyz(1:3)=xt(i-2:i)*0.5291772108_gp !non-BigDFT
                else
                    xyz(1:3)=xt(i-2:i)
                endif
                write(1388,'(a,1x,3e24.15)') symbol,xyz(1),xyz(2),xyz(3) !,xt(i-2:i)
            endif
        enddo
    enddo
    close(1388)
    !-------------------------------------------------------------------------------------------
    deallocate(xt,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating xt.'
    deallocate(y,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating y.'
    deallocate(s,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating s.'
    deallocate(h,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating h.'
    deallocate(e1,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating e1.'
    deallocate(e2,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating e2.'
    deallocate(c,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating c.'
end subroutine writepathway


subroutine writeanchorpoints(n,np,x,filename,atoms)
    use module_types
    implicit none
    integer::n,np,iat
    real(kind=8)::x(n,0:np),xyz(3)
    integer::ip,i
    type(atoms_data), intent(in) :: atoms
    character(len=10) :: name
    character(len=2) :: symbol
    character(*)::filename
    !-------------------------------------------------------------------------------------------
    open(unit=1389,file=filename,status='replace')
    write(1389,*) np
    do ip=0,np
        write(1389,*) n/3
        write(1389,*) ip,'   '
        do i=1,n
            if(mod(i,3)==0) then
                iat=i/3
                name=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
                if (name(3:3)=='_') then
                   symbol=name(1:2)
                   !suffix=name(4:6)
                else if (name(2:2)=='_') then
                   symbol=name(1:1)
                   !suffix=name(3:5)
                else
                   symbol=name(1:2)
                   !suffix=' '
                end if
                if(atoms%astruct%units=='angstroemd0' .or. atoms%astruct%units=='angstroem') then
                    xyz(1:3)=x(i-2:i,ip)*0.5291772108_gp !non-BigDFT
                else
                    xyz(1:3)=x(i-2:i,ip)
                endif
                write(1389,'(a,1x,3e24.15)') symbol,xyz(1),xyz(2),xyz(3)
            endif
        enddo
    enddo
    close(1389)
    !-------------------------------------------------------------------------------------------
end subroutine writeanchorpoints


subroutine readanchorpoints(n,np,x,filename,atoms)
    use module_types
    implicit none
    integer::n,np,ip,i !n(c) iat
    real(kind=8)::x(n,0:100),xyz(3)
    type(atoms_data), intent(in) :: atoms
    !integer::i
    character(len=10) :: tname
    character(*)::filename
    !-------------------------------------------------------------------------------------------
    open(unit=1390,file=filename,status='old')
    read(1390,*) np
    if(np>100) stop 'ERROR: in readanchorpoints np>100'
    do ip=0,np
        read(1390,*)
        read(1390,*)
        do i=1,n
            if(mod(i,3)==0) then
                !n(c) iat=i/3
                read(1390,*) tname,xyz(1),xyz(2),xyz(3) !x(i-2:i,ip)
                if(atoms%astruct%units=='angstroemd0' .or. atoms%astruct%units=='angstroem') then
                    x(i-2:i,ip)=xyz(1:3)/0.5291772108_gp !non-BigDFT
                else
                    x(i-2:i,ip)=xyz(1:3)
                endif
            endif
        enddo
    enddo
    close(1390)
    !-------------------------------------------------------------------------------------------
end subroutine readanchorpoints


subroutine initializepoints(atoms,n,x1,x2,np,x)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer::n,np,ip,i,iat,ixyz
    real(kind=8)::x1(n),x2(n),x(n,0:np),dt,t,tt
    logical::move_this_coordinate
    dt=1.d0/np
    !linear interpolation between the two ends.
    x(1:n,0)=x1(1:n)
    do ip=1,np-1
        t=dt*ip
        !if(np==2) t=0.3d0
        !x(1:n,np)=(1.d0-t)*x1(1:n)+t*x2(1:n)
        !NOTICE: atomic_copycoord does not copy frozen atoms
        !call atomic_copycoord(atoms,x(1,np),x(1,ip))

        !tt=-tan(2.d0*t-1.d0)**3/(2.d0*tan(1.d0)**3)+0.5d0
        !tt=1.d0/(1.d0+exp((t-0.5d0)*12.d0));
        tt=1.d0-t
        !if(ip<=np/2) then
        !    tt=1.d0-0.6d0*ip/real(np,8)
        !else
        !    tt=0.6d0*(np-ip)/real(np,8)
        !endif
        do i=1,n
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(move_this_coordinate(atoms%astruct%ifrztyp(iat),ixyz)) then
                x(i,ip)=tt*x1(i)+(1.d0-tt)*x2(i)
            else
                x(i,ip)=x1(i)
            endif
        enddo
    enddo
    x(1:n,np)=x2(1:n)
    !---------------------------------------------
end subroutine initializepoints

!debug mode desactivated, should be inserted in f_malloc functionalities
subroutine dmemocc(n1,n2,v,chv)
    use module_base
    implicit none
    integer::n1,n2,i
    real(kind=8)::v(n2)
    character(*)::chv
    if(n1>n2) then
        write(777,'(2a)') 'ERROR: n1>n2 while allocating ',chv
    elseif(n1==n2) then
        return
    else
        write(777,'(2a)') 'padding the double-precision array ',chv
        do i=n1+1,n2
            v(i)=0.0d0!d_nan()
        enddo
    endif
end subroutine dmemocc


subroutine imemocc(n1,n2,v,chv)
    use module_base
    implicit none
    !Arguments
    integer :: n1,n2,i,v(n2)
    character(len=*) :: chv

    if(n1>n2) then
        write(888,'(2a)') 'ERROR: n1>n2 while allocating ',chv
    elseif(n1==n2) then
        return
    else
        write(888,'(2a)') 'padding the integer array ',chv
        do i=n1+1,n2
            v(i)=0!int(r_nan())
        enddo
    endif
end subroutine imemocc
