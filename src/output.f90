!>  @file
!!  File where most relevant screen output are collected
!!  Routines which are present in this file should have *all* arguments as intent(in)
!!  Also, the master process only should acces these routines
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Display the logo of BigDFT subroutine print_logo()
subroutine print_logo()
  use module_base
  use yaml_output
  implicit none
  integer :: namelen,ierr
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  integer :: nthreads
  !integer, parameter :: ln = 1024
!$ integer :: omp_get_max_threads

!  call yaml_comment('Daubechies Wavelets for DFT Pseudopotential Calculations',hfill='=')


  call yaml_mapping_open('Code logo')
!logo of BigDFT, new version
  call yaml_scalar('"__________________________________ A fast and precise DFT wavelet code')
  call yaml_scalar('|     |     |     |     |     |                                        ')
  call yaml_scalar('|     |     |     |     |     |      BBBB         i       gggggg       ')
  call yaml_scalar('|_____|_____|_____|_____|_____|     B    B               g             ')
  call yaml_scalar('|     |  :  |  :  |     |     |    B     B        i     g              ')
  call yaml_scalar('|     |-0+--|-0+--|     |     |    B    B         i     g        g     ')
  call yaml_scalar('|_____|__:__|__:__|_____|_____|___ BBBBB          i     g         g    ')
  call yaml_scalar('|  :  |     |     |  :  |     |    B    B         i     g         g    ')
  call yaml_scalar('|--+0-|     |     |-0+--|     |    B     B     iiii     g         g    ')                                  
  call yaml_scalar('|__:__|_____|_____|__:__|_____|    B     B        i      g        g    ')
  call yaml_scalar('|     |  :  |  :  |     |     |    B BBBB        i        g      g     ')
  call yaml_scalar('|     |-0+--|-0+--|     |     |    B        iiiii          gggggg      ')
  call yaml_scalar('|_____|__:__|__:__|_____|_____|__BBBBB                                 ')
  call yaml_scalar('|     |     |     |  :  |     |                           TTTTTTTTT    ')
  call yaml_scalar('|     |     |     |--+0-|     |  DDDDDD          FFFFF        T        ')
  call yaml_scalar('|_____|_____|_____|__:__|_____| D      D        F        TTTT T        ')
  call yaml_scalar('|     |     |     |  :  |     |D        D      F        T     T        ')
  call yaml_scalar('|     |     |     |--+0-|     |D         D     FFFF     T     T        ')
  call yaml_scalar('|_____|_____|_____|__:__|_____|D___      D     F         T    T        ')
  call yaml_scalar('|     |     |  :  |     |     |D         D     F          TTTTT        ')
  call yaml_scalar('|     |     |--+0-|     |     | D        D     F         T    T        ')
  call yaml_scalar('|_____|_____|__:__|_____|_____|          D     F        T     T        ')                                       
  call yaml_scalar('|     |     |     |     |     |         D               T    T         ')
  call yaml_scalar('|     |     |     |     |     |   DDDDDD       F         TTTT          ')
  call yaml_scalar('|_____|_____|_____|_____|_____|______                    www.bigdft.org   "') 

  !old version 
!!$  call yaml_scalar('      TTTT         F       DDDDD    ')     
!!$  call yaml_scalar('     T    T               D         ')     
!!$  call yaml_scalar('    T     T        F     D          ')     
!!$  call yaml_scalar('    T    T         F     D        D ')     
!!$  call yaml_scalar('    TTTTT          F     D         D')     
!!$  call yaml_scalar('    T    T         F     D         D')     
!!$  call yaml_scalar('    T     T        F     D         D')     
!!$  call yaml_scalar('    T      T       F     D         D')     
!!$  call yaml_scalar('    T     T     FFFF     D         D')     
!!$  call yaml_scalar('    T TTTT         F      D        D')     
!!$  call yaml_scalar('    T             F        D      D ')     
!!$  call yaml_scalar('TTTTTTTTT    FFFFF          DDDDDD  ')     
!!$  call yaml_scalar('  gggggg          iiiii    BBBBBBBBB')     
!!$  call yaml_scalar(' g      g        i             B    ')     
!!$  call yaml_scalar('g        g      i         BBBB B    ')     
!!$  call yaml_scalar('g         g     iiii     B     B    ')     
!!$  call yaml_scalar('g         g     i       B      B    ')     
!!$  call yaml_scalar('g         g     i        B     B    ')     
!!$  call yaml_scalar('g         g     i         B    B    ')     
!!$  call yaml_scalar('g         g     i          BBBBB    ')     
!!$  call yaml_scalar(' g        g     i         B    B    ')     
!!$  call yaml_scalar('          g     i        B     B    ')     
!!$  call yaml_scalar('         g               B    B     ')     
!!$  call yaml_scalar('    ggggg       i         BBBB      ') 
  call yaml_mapping_close()

  call yaml_map('Reference Paper','The Journal of Chemical Physics 129, 014109 (2008)')
  call yaml_map('Version Number',package_version)
  call yaml_map('Timestamp of this run',yaml_date_and_time_toa())

  call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
  if (ierr ==0) call yaml_map('Root process Hostname',trim(nodename_local))
  call yaml_map('Number of MPI tasks',bigdft_mpi%nproc)
  
  nthreads = 0
!$  nthreads=omp_get_max_threads()
  call yaml_map('OpenMP parallelization',nthreads>0)
  if (nthreads > 0) then
     call yaml_map('Maximal OpenMP threads per MPI task',nthreads)
  endif

END SUBROUTINE print_logo


subroutine print_configure_options()
  use yaml_output
  implicit none
  integer, parameter :: ln = 1024
  character(len = ln), dimension(4) :: buf

  call yaml_comment('Code compiling options',hfill='-')
  call yaml_mapping_open("Compilation options")
  call bigdft_config_get_user_args(buf(1), ln)
  call yaml_map("Configure arguments", '"'//trim(buf(1))//'"')
  call bigdft_config_get_compilers(buf(1), buf(2), buf(3), ln)
  call yaml_map("Compilers (CC, FC, CXX)", buf(1:3))
  call bigdft_config_get_compiler_flags(buf(1), buf(2), buf(3), buf(4), ln)
  call yaml_mapping_open("Compiler flags")
  if (len_trim(buf(1))>0) call yaml_map("CFLAGS",   trim(buf(1)))
  if (len_trim(buf(2))>0) call yaml_map("FCFLAGS",  trim(buf(2)))
  if (len_trim(buf(3))>0) call yaml_map("CXXFLAGS", trim(buf(3)))
  if (len_trim(buf(4))>0) call yaml_map("CPPFLAGS", trim(buf(4)))
  call yaml_mapping_close()
!!$  call bigdft_config_get_linker(buf(1), buf(2), buf(3), buf(4), ln)
!!$  call yaml_mapping_open("Linker")
!!$   call yaml_map("LD",      trim(buf(1)))
!!$   call yaml_map("LDFLAGS", trim(buf(2)))
!!$   call yaml_map("LIBS",    trim(buf(3)))
!!$   call yaml_map("Full linking options", trim(buf(4)))
!!$  call yaml_mapping_close()

 call yaml_mapping_close()

end subroutine print_configure_options


!> Print all general parameters
subroutine print_general_parameters(in,atoms)
  use module_base
  use module_types
  use defs_basis
  use yaml_output
  use module_input_keys, only: input_keys_equal
  implicit none
  !Arguments
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms

  integer :: iat, i
  character(len=len(in%run_name)) :: prefix
  character(len = 11) :: potden
  character(len = 12) :: dos

  ! Output for atoms
  if (trim(in%run_name) == '') then
     call yaml_comment('Input Atomic System (file: posinp.'//trim(atoms%astruct%inputfile_format)//')',hfill='-')
     prefix = 'input'
  else
     prefix = in%run_name
     call yaml_comment('Input Atomic System (file: '//trim(prefix)//'.'//trim(atoms%astruct%inputfile_format)//')',hfill='-')
  end if

  ! Atomic systems
  call yaml_mapping_open('Atomic System Properties')
     call yaml_map('Number of atomic types', atoms%astruct%ntypes, fmt='(i0)')
     call yaml_map('Number of atoms', atoms%astruct%nat, fmt='(i0)')
     if (atoms%astruct%nat > 0) then
        call yaml_map('Types of atoms',atoms%astruct%atomnames)
        ! Fixed positions
        if (maxval(atoms%astruct%ifrztyp) /= 0) then
           call yaml_sequence_open('Fixed atoms',flow=.true.)
           ! The fixed atom column
           do iat=1,atoms%astruct%nat
              if (atoms%astruct%ifrztyp(iat) /= 0) then
                 call yaml_sequence('at.' // trim(yaml_toa(iat,fmt='(i4.4)')) // &
                      & '(' // trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))) // ')' &
                      & // trim(yaml_toa(atoms%astruct%ifrztyp(iat),fmt='(i0)')))
              end if
           end do
           call yaml_sequence_close()
        end if
     end if
     !Boundary Conditions
     select case(atoms%astruct%geocode)
     case('P')
        call yaml_map('Boundary Conditions','Periodic',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
        call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
             atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
     case('S')
        call yaml_map('Boundary Conditions','Surface',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
        call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
             atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
     case('W')
        call yaml_map('Boundary Conditions','Wire',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
        call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
             atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
     case('F')
        call yaml_map('Boundary Conditions','Free',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
     end select

     !Symmetries
     call yaml_map('Number of Symmetries',atoms%astruct%sym%nSym)
     call yaml_map('Space group',trim(atoms%astruct%sym%spaceGroup))
  call yaml_mapping_close()

  !Geometry imput Parameters
  if (in%ncount_cluster_x > 0) then
     call yaml_comment('Geometry optimization Input Parameters (file: '//trim(prefix)//'.geopt)',hfill='-')
     call yaml_mapping_open('Geometry Optimization Parameters')
        call yaml_map('Maximum steps',in%ncount_cluster_x)
        call yaml_map('Algorithm', in%geopt_approach)
        call yaml_map('Random atomic displacement', in%randdis, fmt='(1pe7.1)')
        call yaml_map('Fluctuation in forces',in%frac_fluct,fmt='(1pe7.1)')
        call yaml_map('Maximum in forces',in%forcemax,fmt='(1pe7.1)')
        call yaml_map('Steepest descent step',in%betax,fmt='(1pe7.1)')
        if (trim(in%geopt_approach) == "DIIS") then
           call yaml_map('DIIS history', in%history)
        end if
     call yaml_mapping_close()
     if (input_keys_equal(trim(in%geopt_approach),"AB6MD")) then
        call yaml_mapping_open('Molecular Dynamics Parameters')
           call yaml_map('ionmov',in%ionmov)
           call yaml_map('dtion', in%dtion,fmt='(0pf7.0)')
           if (in%ionmov > 7) then
              call yaml_map('Start Temperature', in%mditemp, fmt='(f5.0)')
              call yaml_map('Stop Temperature',  in%mdftemp, fmt='(f5.0)')
           end if
           if (in%ionmov == 8) then
              call yaml_map('Nose inertia', in%noseinert,fmt='(f15.5)')
           else if (in%ionmov == 9) then
              call yaml_map('Friction', in%friction,fmt='(f15.5)')
              call yaml_map('MD wall',in%mdwall,fmt='(f15.5)')
           else if (in%ionmov == 13) then
              call yaml_map('nnos', in%nnos,fmt='(f15.5)')
              call yaml_map('qmass',in%qmass,fmt='(f15.5)')
              call yaml_map('bmass',in%bmass,fmt='(f15.5)')
              call yaml_map('vmass',in%vmass,fmt='(f15.5)')
           end if
        call yaml_mapping_close()
     end if
     !write(*,'(1x,a)') '--- (file: input.geopt) ------------------------------------- Geopt Input Parameters'
     !write(*, "(A)")   "       Generic param.              Geo. optim.                MD param."
     !write(*, "(1x,a,i7,1x,a,1x,a,1pe7.1,1x,a,1x,a,i7)") &
     !     & "      Max. steps=", in%ncount_cluster_x, "|", &
     !     & "Fluct. in forces=", in%frac_fluct,       "|", &
     !     & "          ionmov=", in%ionmov
     !write(*, "(1x,a,a7,1x,a,1x,a,1pe7.1,1x,a,1x,a,0pf7.0)") &
     !     & "       algorithm=", in%geopt_approach, "|", &
     !     & "  Max. in forces=", in%forcemax,       "|", &
     !     & "           dtion=", in%dtion
     !if (trim(in%geopt_approach) /= "DIIS") then
     !   write(*, "(1x,a,1pe7.1,1x,a,1x,a,1pe7.1,1x,a)", advance="no") &
     !        & "random at.displ.=", in%randdis, "|", &
     !        & "  steep. descent=", in%betax,   "|"
     !else
     !   write(*, "(1x,a,1pe7.1,1x,a,1x,a,1pe7.1,2x,a,1I2,1x,a)", advance="no") &
     !        & "random at.displ.=", in%randdis,           "|", &
     !        & "step=", in%betax, "history=", in%history, "|"
     !end if
     !if (in%ionmov > 7) then
     !   write(*, "(1x,a,1f5.0,1x,a,1f5.0)") &
     !        & "start T=", in%mditemp, "stop T=", in%mdftemp
     !else
     !   write(*,*)
     !end if
     !
     !if (in%ionmov == 8) then
     !   write(*,'(1x,a,f15.5)') "TODO: pretty printing!", in%noseinert
     !else if (in%ionmov == 9) then
     !   write(*,*) "TODO: pretty printing!", in%friction
     !   write(*,*) "TODO: pretty printing!", in%mdwall
     !else if (in%ionmov == 13) then
     !   write(*,*) "TODO: pretty printing!", in%nnos
     !   write(*,*) "TODO: pretty printing!", in%qmass
     !   write(*,*) "TODO: pretty printing!", in%bmass, in%vmass
     !end if
  end if

  !Output for K points
  if (atoms%astruct%geocode /= 'F') then
     call yaml_comment('K points description (Reduced and Brillouin zone coordinates, Weight)',hfill='-')
     !write(*,'(1x,a)') '--- (file: input.kpt) ----------------------------------------------------- k-points'
     if (in%disableSym .and. in%gen_nkpt > 1) then
        call yaml_warning('symmetries have been disabled, k points are not irreductible.')
        !write(*, "(1x,A)") "WARNING: symmetries have been disabled, k points are not irreductible."
     end if
     call yaml_sequence_open('K points')!,advance='no')
     !call yaml_comment('Reduced coordinates  BZ coordinates  weight',hfill=' ')
     !write(*, "(1x,a)")    "       red. coordinates         weight       id        BZ coordinates"
     do i = 1, in%gen_nkpt, 1
        call yaml_sequence(advance='no')
        call yaml_mapping_open(flow=.true.)
          call yaml_map( 'Rc', &
             & in%gen_kpt(:, i) * atoms%astruct%cell_dim / two_pi,&
             & fmt='(f7.4)')
          call yaml_map( 'Bz', &
             & in%gen_kpt(:, i), &
             & fmt='(f7.4)')
          call yaml_map('Wgt',in%gen_wkpt(i),fmt='(f6.4)')
        call yaml_mapping_close(advance='no')
        call yaml_comment(trim(yaml_toa(i,fmt='(i4.4)')))
        !write(*, "(1x,3f9.5,2x,f9.5,5x,I4,1x,3f9.5)") &
        !     & in%kpt(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), atoms%astruct%cell_dim(3) /) / two_pi, &
        !     & in%wkpt(i), i, in%kpt(:, i)
     end do
     call yaml_sequence_close()

     if (in%nkptv > 0) then
        call yaml_sequence_open('K points for band structure calculation')
        !write(*, "(1x,a)")    " K points for band structure calculation"
        !write(*, "(1x,a)")    "       red. coordinates         weight       id        BZ coordinates"
        do i = 1, in%nkptv, 1
          call yaml_sequence(advance='no')
          call yaml_mapping_open(trim(yaml_toa(i,fmt='(i0)')),flow=.true.)
          call yaml_map( 'Red C.', &
             & in%kptv(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), &
             & atoms%astruct%cell_dim(3) /) / two_pi,&
             & fmt='(f9.5)')
          call yaml_map( 'Bz C.', &
             & in%kptv(:, i), &
             & fmt='(f9.5)')
          call yaml_map('Weight',1.0d0 / real(size(in%kptv, 2), gp),fmt='(f9.5)')
          call yaml_mapping_close()
        !   write(*, "(1x,3f9.5,2x,f9.5,5x,I4,1x,3f9.5)") &
        !        & in%kptv(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), atoms%astruct%cell_dim(3) /) / two_pi, &
        !        & 1.0d0 / real(size(in%kptv, 2), gp), i, in%kptv(:, i)
        end do
        call yaml_sequence_close()
     end if
  end if

  ! Printing for mixing parameters.
  if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) then
     if (in%iscf < 10) then
        write(potden, "(A)") "potential"
     else
        write(potden, "(A)") "density"
     end if
     call yaml_comment('Mixing (file: '//trim(prefix)//'.mix)',hfill='-')
     call yaml_mapping_open('Mixing parameters')
        call yaml_map('Target',trim(potden))
        call yaml_map('Additional bands', in%norbsempty)
        call yaml_map('Mixing Coefficients', in%alphamix,fmt='(0pe10.2)')
        call yaml_map('Scheme',modulo(in%iscf, 10))
        call yaml_map('Electronic temperature',in%Tel,fmt='(1pe12.2)')
        call yaml_map('DIIS',in%alphadiis,fmt='(0pe12.2)')
        call yaml_map('Maximum iterations',in%itrpmax)
        call yaml_map('Occupied scheme',trim(smearing_names(in%occopt)))
        call yaml_map('Rp norm',in%rpnrm_cv,fmt='(1pe12.2)')
        if (in%verbosity > 2) then
           write(dos, "(A)") "dos.gnuplot"
        else
           write(dos, "(A)") "no verb. < 3"
        end if
        call yaml_map('output DOS',trim(dos))
     call yaml_mapping_close()
     !write(*,'(1x,a)') '--- (file: input.mix) ------------------------------------------------------- Mixing'
     !write(*,"(1x,A12,A12,1x,A1,1x,A12,I12,1x,A1,1x,A11,F10.2)") &
     !     & "     Target=", potden,        "|", &
     !     & " Add. bands=", in%norbsempty, "|", &
     !     & "    Coeff.=", in%alphamix
     !write(*,"(1x,A12,I12,1x,A1,1x,A12,1pe12.2,1x,A1,1x,A11,0pe10.2)") &
     !     & "     Scheme=", modulo(in%iscf, 10), "|", &
     !     & "Elec. temp.=", in%Tel,              "|", &
     !     & "      DIIS=", in%alphadiis
     !write(*,"(1x,A12,I12,1x,A1,1x,A12,A12,1x,A1)") &
     !     & "  Max iter.=", in%itrpmax,    "|", &
     !     & "Occ. scheme=", smearing_names(in%occopt), "|"
     !if (in%verbosity > 2) then
     !   write(dos, "(A)") "dos.gnuplot"
     !else
     !   write(dos, "(A)") "no verb. < 3"
     !end if
     !write(*,"(1x,A12,1pe12.2,1x,A1,1x,2A12,1x,A1)") &
     !    & "   Rp norm.=", in%rpnrm_cv,    "|", " output DOS=", dos, "|"
  end if

END SUBROUTINE print_general_parameters


!> Print all dft input parameters
subroutine print_dft_parameters(in,atoms)
  use module_base
  use module_types
  use yaml_output
  use module_xc
  implicit none
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms

  call yaml_comment('Input parameters',hfill='-')

  call yaml_mapping_open('DFT parameters')
    call yaml_mapping_open('eXchange Correlation')
      call yaml_map('XC ID',in%ixc,fmt='(i8)',label='ixc')
      if (in%ixc < 0) then
         call xc_dump(in%ixc, XC_MIXED, in%nspin)
      else
         call xc_dump(in%ixc, XC_ABINIT, in%nspin)
      end if
      if (in%nspin>=2) then
         call yaml_map('Polarisation',in%mpol)
         !write(*,'(1x,a,i7,1x,a)')&
         !     'Polarisation=',in%mpol, '|'
      end if
      if (in%nspin==4) then
         call yaml_map('Spin polarization','non collinear')
      else if (in%nspin==2) then
         call yaml_map('Spin polarization','collinear')
      else if (in%nspin==1) then
         call yaml_map('Spin polarization',.false.)
      end if
    call yaml_mapping_close()

    if (in%ncharge > 0) call yaml_map('Net Charge (Ions-Electrons)',in%ncharge,fmt='(i8)')
    if (sqrt(sum(in%elecfield(:)**2)) > 0.0_gp) &
      call yaml_map('External Electric Field (Ha/a0)',in%elecfield(:),fmt='(1pe8.1)')
  call yaml_mapping_close()

  call yaml_mapping_open('Basis set definition')
    call yaml_map('Suggested Grid Spacings (a0)', (/in%hx,in%hy,in%hz/),fmt='(f5.2)')
    call yaml_map('Coarse and Fine Radii Multipliers', (/in%crmult,in%frmult/),fmt='(f4.1)')
  call yaml_mapping_close()


  call yaml_mapping_open('Self-Consistent Cycle Parameters')
    call yaml_mapping_open('Wavefunction')
      call yaml_map('Gradient Norm Threshold',in%gnrm_cv,fmt='(1pe8.1)',label='gnrm_cv')
      call yaml_map('CG Steps for Preconditioner',in%ncong,fmt='(i5)')
      call yaml_map('DIIS History length',in%idsx)
      call yaml_map('Max. Wfn Iterations',in%itermax,label='itermax')
      call yaml_map('Max. Subspace Diagonalizations',in%nrepmax)
      call yaml_map('Input wavefunction policy',  trim(input_psi_names(in%inputPsiId)), advance="no")
      call yaml_comment(trim(yaml_toa(in%inputPsiId)))
      call yaml_map('Output wavefunction policy', trim(wf_format_names(in%output_wf_format)), advance="no")
      call yaml_comment(trim(yaml_toa(in%output_wf_format)))
      call yaml_map('Output grid policy',trim(output_denspot_names(in%output_denspot)),advance='no')
      call yaml_comment(trim(yaml_toa(in%output_denspot)))
      call yaml_map('Output grid format',trim(output_denspot_format_names(in%output_denspot_format)),advance='no')
      call yaml_comment(trim(yaml_toa(in%output_denspot_format)))
      call yaml_map('Virtual orbitals',in%nvirt,fmt='(i0)')
      call yaml_map('Number of plotted density orbitals',abs(in%nplot),fmt='(i0)')
  
    call yaml_mapping_close()

    call yaml_mapping_open('Density/Potential')
       call yaml_map('Max. Iterations',in%itrpmax)
    call yaml_mapping_close()
  call yaml_mapping_close()

  if (atoms%astruct%geocode == 'F') then
     call yaml_mapping_open('Post Optimization Parameters')

     call yaml_mapping_open('Finite-Size Effect estimation')
     call yaml_map('Scheduled',(in%rbuf > 0.0_gp))
     if (in%rbuf > 0.0_gp) then
        call yaml_map('Extension',in%rbuf,fmt='(f4.1)')
        call yaml_map('No. of CG steps',in%ncongt)
     end if
     call yaml_mapping_close()
     call yaml_mapping_close()
  end if

!!$  write(*,'(1x,a)')&
!!$       '--- (file: input.dft) --------------------------------------------- Input Parameters'
!!$  write(*,'(1x,a)')&
!!$       '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
!!$  write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
!!$       '  Max. hgrid=',in%hx,   '|  Coarse Wfs.=',in%crmult,'| Wavefns Conv.=',in%gnrm_cv,&
!!$       '| Calculate=',(in%rbuf > 0.0_gp)
!!$  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i5,a,i2,1x,a,f4.1)')&
!!$       '       XC id=',in%ixc,     '|    Fine Wfs.=',in%frmult,'| Max. N. Iter.=',in%itermax,&
!!$       'x',in%nrepmax,'| Extension=',in%rbuf
!!$  write(*,'(1x,a,i7,1x,a,1x,a,i8,1x,a,i4)')&
!!$       'total charge=',in%ncharge, '|                   ','| CG Prec.Steps=',in%ncong,&
!!$       '|  CG Steps=',in%ncongt
!!$  write(*,'(1x,a,1pe7.1,1x,a,1x,a,i8)')&
!!$       ' elec. field=',sqrt(sum(in%elecfield(:)**2)),'|                   ','| DIIS Hist. N.=',in%idsx


  !write(*, "(1x,A19,I5,A,1x,A1,1x,A19,I6,A)") &
  !     & "Input wf. policy=", in%inputPsiId, " (" // input_psi_names(in%inputPsiId) // ")", "|", &
  !     & "Output wf. policy=", in%output_wf_format, " (" // wf_format_names(in%output_wf_format) // ")"

  !write(*, "(1x,A19,I5,A,1x,A1,1x,A19,I6,A)") &
  !     & "Output grid policy=", in%output_denspot, "   (" // output_denspot_names(in%output_denspot) // ")", "|", &
  !     & "Output grid format=", in%output_denspot_format, &
  !     "         (" // output_denspot_format_names(in%output_denspot_format) // ")"

END SUBROUTINE print_dft_parameters


!> Write input parameters
subroutine write_input_parameters(in)!,atoms)
  use module_base
  use module_types
  use yaml_output
  implicit none
  type(input_variables), intent(in) :: in
!  type(atoms_data), intent(in) :: atoms
  !local variables
  character(len = 11) :: potden
  !start yaml output
  !call yaml_indent_map('Physical System Parameters')
  !write(70,'(a)')repeat(' ',yaml_indent)//'Physical System Parameters:'
  !yaml_indent=yaml_indent+3
!  write(70,'(a,t55,a)')repeat(' ',yaml_indent)//'Boundary Conditions:',atoms%astruct%geocode
!  if (atoms%astruct%geocode /= 'F')write(70,'(a,t55,a,3(1x,f5.3,a))')&
!       repeat(' ',yaml_indent)//'Box Sizes (a0):','[',atoms%astruct%cell_dim(1),',',atoms%astruct%cell_dim(2),',',atoms%astruct%cell_dim(3),' ]'

!  yaml_indent=yaml_indent-3


  if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) then
     !write(70,'(a)')repeat(' ',yaml_indent)//'Mixing Parameters:'
     !yaml_indent=yaml_indent+3
       if (in%iscf < 10) then
        write(potden, "(A)") "potential"
     else
        write(potden, "(A)") "density"
     end if
     write(70,'(a,t55,a)')'Target:',potden
!     write(70,'(a,t55,I12)')'Scheme:',modulo(in%iscf, 10)
!!$     write(*,"(1x,A12,A12,1x,A1,1x,A12,I12,1x,A1,1x,A11,F10.2)") &
!!$          & "     Target=", potden,        "|", &
!!$          & " Add. bands=", in%norbsempty, "|", &
!!$          & "    Coeff.=", in%alphamix
!!$     write(*,"(1x,A12,I12,1x,A1,1x,A12,1pe12.2,1x,A1,1x,A11,0pe10.2)") &
!!$          & "     Scheme=", modulo(in%iscf, 10), "|", &
!!$          & "Elec. temp.=", in%Tel,              "|", &
!!$          & "      DIIS=", in%alphadiis
!!$     write(*,"(1x,A12,I12,1x,A1,1x,A12,A12,1x,A1)") &
!!$          & "  Max iter.=", in%itrpmax,    "|", &
!!$          & "Occ. scheme=", smearing_names(in%occopt), "|"
!!$     if (in%verbosity > 2) then
!!$        write(dos, "(A)") "dos.gnuplot"
!!$     else
!!$        write(dos, "(A)") "no verb. < 3"
!!$     end if
!!$     write(*,"(1x,A12,1pe12.2,1x,A1,1x,2A12,1x,A1)") &
!!$          & "   Rp norm.=", in%rpnrm_cv,    "|", " output DOS=", dos, "|"
  end if
  !write(70,'(a)')repeat(' ',yaml_indent)//'Post Optimization Treatments:'
  if (in%rbuf > 0.0_gp) then
     !write(70,'(a)')repeat(' ',yaml_indent)//'Finite-Size Correction Estimation:'
     !write(70,'(a,t55,f4.1)')repeat(' ',yaml_indent)//'Radius (a0):',in%rbuf
     !write(70,'(a,t55,i4)')repeat(' ',yaml_indent)//'CG Steps for the FS Correction:',in%ncongt
  end if
  stop
end subroutine write_input_parameters


!> Write the energies for a given iteration
subroutine write_energies(iter,iscf,energs,gnrm,gnrm_zero,comment,only_energies)
  use module_base
  use module_types
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: iter !< Iteration Id
  integer, intent(in) :: iscf
  type(energy_terms), intent(in) :: energs
  real(gp), intent(in) :: gnrm,gnrm_zero
  character(len=*), intent(in) :: comment
  logical,intent(in),optional :: only_energies
  !local variables
  logical :: write_only_energies

  if (present(only_energies)) then
      write_only_energies=only_energies
  else
      write_only_energies=.false.
  end if

  if (len(trim(comment)) > 0 .and. .not.write_only_energies) then
     if (verbose >0) call yaml_newline()
     call write_iter()
     if (verbose >0) call yaml_comment(trim(comment))    
  end if

  if (iscf < 1 .and. verbose > 0) then
     call yaml_newline()
     call yaml_mapping_open('Energies',flow=.true.)
  !call yaml_flow_map()
  !call yaml_indent_map('Energies')
     if (energs%ekin /= 0.0_gp)&
          call yaml_map('Ekin',energs%ekin,fmt='(1pe18.11)')
     if (energs%epot /= 0.0_gp)&
          call yaml_map('Epot',energs%epot,fmt='(1pe18.11)')
     if (energs%eproj /= 0.0_gp)&
          call yaml_map('Enl',energs%eproj,fmt='(1pe18.11)')
     if (energs%eh /= 0.0_gp)&
          call yaml_map('EH',energs%eh,fmt='(1pe18.11)')
     if (energs%exc /= 0.0_gp)&
          call yaml_map('EXC',energs%exc,fmt='(1pe18.11)')
     if (energs%evxc /= 0.0_gp)&
          call yaml_map('EvXC',energs%evxc,fmt='(1pe18.11)')
     if (energs%eexctX /= 0.0_gp)&
          call yaml_map('EexctX',energs%eexctX,fmt='(1pe18.11)')
     if (energs%evsic /= 0.0_gp)&
          call yaml_map('EvSIC',energs%evsic,fmt='(1pe18.11)')
     if (len(trim(comment)) > 0) then
        if (energs%eion /= 0.0_gp)&
             call yaml_map('Eion',energs%eion,fmt='(1pe18.11)')
        if (energs%edisp /= 0.0_gp)&
             call yaml_map('Edisp',energs%edisp,fmt='(1pe18.11)')
        if (energs%excrhoc /= 0.0_gp)&
             call yaml_map('Exc(rhoc)',energs%excrhoc,fmt='(1pe18.11)')
        if (energs%eTS /= 0.0_gp)&
             call yaml_map('TS',energs%eTS,fmt='(1pe18.11)')

     end if
     call yaml_mapping_close()
  end if

  if (.not.write_only_energies) then
     call yaml_newline()
     if (len(trim(comment)) == 0) then
        call write_iter()
        if (verbose >0) call yaml_newline()
     else if (verbose > 1) then
        call yaml_map('SCF criterion',iscf,fmt='(i6)')
     end if
  end if


  contains

    subroutine write_iter()
      implicit none
      if (iter > 0) call yaml_map('iter',iter,fmt='(i6)')
      if (iscf > 1) then
         call yaml_map('tr(H)',energs%trH,fmt='(1pe24.17)')
      else
         if (energs%eTS==0.0_gp) then
            call yaml_map('EKS',energs%energy,fmt='(1pe24.17)')
         else
            call yaml_map('FKS',energs%energy,fmt='(1pe24.17)')
         end if
      end if
      if (gnrm > 0.0_gp) call yaml_map('gnrm',gnrm,fmt='(1pe9.2)')
      if (gnrm_zero > 0.0_gp) &
           call yaml_map('gnrm0',gnrm_zero,fmt='(1pe8.1)')
      if (iscf > 1) then
         if (energs%trH_prev /=0.0_gp) &
              call yaml_map('D',energs%trH-energs%trH_prev,fmt='(1pe9.2)')
      else
         if (energs%e_prev /=0.0_gp) &
              call yaml_map('D',energs%energy-energs%e_prev,fmt='(1pe9.2)')
      end if

    end subroutine write_iter
end subroutine write_energies


!> Write the eigenvalues-related information
subroutine write_eigenvalues_data(etol,orbs,mom_vec)
  use module_base
  use module_types
  use yaml_output
  implicit none
  real(gp), intent(in) :: etol
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(:,:,:), pointer :: mom_vec
  !local variables
  logical :: degup,degdw
  integer :: ikptw,iorb,ikpt,jorb,isorb,nwrtmsg,ndegen
  real(gp) :: spinsignw,mx,my,mz,mpol,tolerance
  character(len=64) :: message
  character(len=150) :: commentline
  real(wp), dimension(2) :: preval

  commentline=repeat(' ',len(commentline))

  if (etol > 1.0_gp) then
     tolerance=0.0_gp
  else
     tolerance=etol
  end if
  
  ! Calculate and print the magnetisation, no matter the verbosity
  if (orbs%nspin == 2) then
     mpol = 0._gp
     do ikpt=1,orbs%nkpts
        isorb = (ikpt - 1) * orbs%norb
        do iorb = 1, orbs%norbu
           mpol = mpol + orbs%occup(isorb + iorb) * orbs%kwgts(ikpt)
        end do
        do iorb = orbs%norbu + 1, orbs%norb, 1
           mpol = mpol - orbs%occup(isorb + iorb) * orbs%kwgts(ikpt)
        end do
     end do
     call yaml_map("Total magnetization",mpol,fmt='(f9.6)')
  end if

  if (verbose > 1) then
     call yaml_comment('Eigenvalues and New Occupation Numbers')

     call yaml_sequence_open('Orbitals',flow=.true.)
     call yaml_newline()

     do ikpt=1,orbs%nkpts
        if (orbs%nkpts > 1 .and. orbs%nspinor >= 2) then
           write(commentline,"(1x,A,I4.4,A,3F12.6)") &
                &   "Kpt #", ikpt, " BZ coord. = ", orbs%kpts(:, ikpt)
           !write(*,'(a)')trim(commentline)
           call yaml_comment(trim(commentline))
           ikptw=ikpt
        else
           ikptw=UNINITIALIZED(1)
        end if
        preval=0.0_wp
        nwrtmsg=0
        ndegen=0
        isorb = (ikpt - 1) * orbs%norb
        if (orbs%nspin==1.or.orbs%nspinor==4) then
           spinsignw=UNINITIALIZED(1.0_gp)
           do iorb=1,orbs%norb
              if (orbs%nspinor ==4 .and. associated(mom_vec)) then
                 mx=(mom_vec(2,iorb,1)/mom_vec(1,iorb,1))
                 my=(mom_vec(3,iorb,1)/mom_vec(1,iorb,1))
                 mz=(mom_vec(4,iorb,1)/mom_vec(1,iorb,1))
              else
                 mx=UNINITIALIZED(1.0_gp)
                 my=UNINITIALIZED(1.0_gp)
                 mz=UNINITIALIZED(1.0_gp)
              end if
              degup = find_degeneracy_up(iorb+isorb)
              call yaml_sequence()
              call write_orbital_data(orbs%eval(isorb + iorb),orbs%occup(isorb+iorb),&
                   spinsignw,ikptw,mx,my,mz)
              !yaml output (carriage return)
              if (iorb == orbs%norb .and. ikpt == orbs%nkpts) then
                 call yaml_sequence_close(advance='no')
                 !print *,'there',nwrtmsg,message
              end if
              call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
              if (nwrtmsg==1) then
                 call yaml_comment(adjustl(message))
              else
                 call yaml_newline()
                 !call yaml_stream_attributes()
              end if
           end do
        else
           mx=UNINITIALIZED(1.0_gp)
           my=UNINITIALIZED(1.0_gp)
           mz=UNINITIALIZED(1.0_gp)

           do iorb=1,min(orbs%norbu,orbs%norbd)
              jorb=orbs%norbu+iorb
              call yaml_sequence()
              call write_orbital_data(orbs%eval(isorb + iorb),orbs%occup(isorb+iorb),&
                   1.0_gp,ikptw,mx,my,mz)
              call yaml_sequence()
              call write_orbital_data(orbs%eval(isorb + jorb),orbs%occup(isorb+jorb),&
                   -1.0_gp,ikptw,mx,my,mz)
              !yaml output (carriage return)
              degup=find_degeneracy_up(iorb+isorb)
              degdw=find_degeneracy_up(jorb+isorb)
              nwrtmsg=0
              if (degup .or. degdw) nwrtmsg=1
              if (degup .and. degdw) message='  <-deg->  '
              if (iorb == orbs%norbu .and. orbs%norbu==orbs%norbd .and. ikpt == orbs%nkpts) then
                 call yaml_sequence_close(advance='no')
              end if
              call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
              if (nwrtmsg==1) then
                 call yaml_comment(adjustl(message))
              else

                 call yaml_newline()
              end if

           end do
           if (orbs%norbu > orbs%norbd) then
              do iorb=orbs%norbd+1,orbs%norbu
                 call yaml_sequence()
                 call write_orbital_data(orbs%eval(isorb+iorb),orbs%occup(isorb+iorb),&
                      1.0_gp,ikptw,mx,my,mz)
                 !yaml output (carriage return)
                 degup = find_degeneracy_up(iorb+isorb)
                 if (iorb == orbs%norbu .and. ikpt == orbs%nkpts) then
                    call yaml_sequence_close(advance='no')
                 end if
                 call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
                 if (nwrtmsg==1) then
                    call yaml_comment(adjustl(message))
                 else
                    call yaml_newline()
                 end if
              end do
           else if (orbs%norbd > orbs%norbu) then
              do iorb=2*orbs%norbu+1,orbs%norbu+orbs%norbd
                 call yaml_sequence()
                 call write_orbital_data(orbs%eval(isorb+iorb),orbs%occup(isorb+iorb),&
                      -1.0_gp,ikptw,mx,my,mz)
                 !yaml output (carriage return)
                 degdw = find_degeneracy_down(iorb+isorb)
                 if (iorb == orbs%norbu+orbs%norbd .and. ikpt == orbs%nkpts) then
                    call yaml_sequence_close(advance='no')
                 end if
                 call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
                 if (nwrtmsg==1) then
                    call yaml_comment(adjustl(message))
                 else
                    call yaml_newline()
                 end if
              end do
           end if
        end if
     end do
     ! Close the map of Eigenvalues and New Occupations Numbers
     !call yaml_mapping_close()
  end if
  !find fermi level
  if (orbs%efermi /= uninitialized(orbs%efermi)) then
     call yaml_map('Fermi Energy',orbs%efermi,fmt='(1pe21.14)')
  end if


contains

  function find_degeneracy_up(iorb) result(wrt)
    implicit none
    integer, intent(in) :: iorb
!    logical find_degeneracy_up
    !local variables
    logical :: wrt

    if (nwrtmsg==0 .and. iorb+1 < orbs%norbu) then
       wrt = abs(orbs%eval(iorb)-orbs%eval(iorb+1)) <= tolerance
       if (wrt) preval(1)=orbs%eval(iorb)
    else
       wrt = abs(orbs%eval(iorb)-preval(1)) <= tolerance
    end if
    !print *,'etol',etol,orbs%eval(iorb),preval(1),wrt,iorb
    nwrtmsg=0

    !disable degeneracy finder, problems to be fixed
    wrt=.false.

    if (wrt) then
       nwrtmsg=1
       message='  <-deg    '
    end if

    if (.not. wrt) preval(1)=orbs%eval(iorb)
  end function find_degeneracy_up

  function find_degeneracy_down(iorb) result(wrt)
    implicit none
    integer, intent(in) :: iorb
!    logical find_degeneracy_down
    !local variables
    logical :: wrt

    if (nwrtmsg==0 .and. iorb+1 < orbs%norb) then
       wrt = abs(orbs%eval(iorb)-orbs%eval(iorb+1)) <= tolerance
       if (wrt) preval(2)=orbs%eval(iorb)
    else
       wrt = abs(orbs%eval(iorb)-preval(2)) <= tolerance
    end if

    !disable degeneracy finder, problems to be fixed
    wrt=.false.


    nwrtmsg=0
    if (wrt) then
       nwrtmsg=1
       message='    deg->  '
    end if

    if (.not. wrt) preval(2)=orbs%eval(iorb)
   
  end function find_degeneracy_down

END SUBROUTINE write_eigenvalues_data


!>Writing rules, control if the last eigenvector is degenerate
!!do this for each spin
!!for each spin it is supposed that only the last group is not completely passed
!!and also that the components of each of the group but the last are the same for up and 
!!down polarisation. Do not work properly in the other cases
subroutine write_ig_eigenvectors(etol,orbse,nspin,norb,norbu,norbd)
   use module_base
   use module_types
   use yaml_output
   implicit none
   integer, intent(in) :: nspin,norb,norbu,norbd
   real(gp), intent(in) :: etol
   type(orbitals_data), intent(in) :: orbse
   !local variables
   character(len=64) :: message
   character(len=25) :: gapstring
   integer :: iorb,ndegen,nwrtmsg,ikpt,iorbst,ikptw
   real(gp) :: HLIGgap,mx,my,mz,spinsignw
   real(wp), dimension(2) :: preval
  character(len=150) :: commentline

  commentline=repeat(' ',len(commentline))


   !loop over all the k-points of the IG 
   iorbst=0 !starting orbital in the k-points distribution
   !check if norbu and norbd are equal
   if (nspin==2 .and. orbse%norbu /= orbse%norbd) then
      write(*,*)'ERROR (write_ig_eigenvectors): the IG orbs structure should have norbu=norbd',orbse%norbu,orbse%norbd
      stop
   end if

  call yaml_sequence_open('Input Guess Orbitals',flow=.true.)!,advance='no')
  call yaml_newline()

  !always without spinors in the IG
  mx=UNINITIALIZED(1.0_gp)
  my=UNINITIALIZED(1.0_gp)
  mz=UNINITIALIZED(1.0_gp)


   do ikpt=1,orbse%nkpts
      if (orbse%nkpts > 1 .and. orbse%nspinor >= 2) then
         write(commentline,"(1x,A,I4.4,A,3F12.6)") &
              &   "Kpt #", ikpt, " BZ coord. = ", orbse%kpts(:, ikpt)
         !write(*,'(a)')trim(commentline)
         call yaml_comment(trim(commentline))
         ikptw=ikpt
      else
         ikptw=UNINITIALIZED(1)
      end if

      preval=0.0_wp
      nwrtmsg=0
      ndegen=0
      do iorb=1,orbse%norbu
         if (nspin==1) then
            spinsignw=UNINITIALIZED(1.0_gp)
            if (nwrtmsg==1) then
               if (abs(orbse%eval(iorb+iorbst)-preval(1)) <= etol) then
                  !degeneracy found
                  message='  <- found degeneracy'
                  ndegen=ndegen+1
               else
                  nwrtmsg=0
               end if
            end if
            if (abs(iorb - norb) <= 5) then
               nwrtmsg=1
               message=' <- '
            end if
            if (iorb == norb) then
               !calculate the IG HOMO-LUMO gap
               if(norb<orbse%norbu) then 
                  HLIGgap=orbse%eval(iorb+1+iorbst)-orbse%eval(iorb+iorbst)
                  write(gapstring,'(a,f8.4,a)') ', H-L IG gap: ',HLIGgap*Ha_eV,' eV'
               else
                  gapstring=''
               end if
               nwrtmsg=1
               message=' <- Last InputGuess eval'//gapstring
               preval(1)=orbse%eval(iorb+iorbst)
            end if
            if (iorb-1 == norb) then
               nwrtmsg=1
               message=' <- First virtual eval '
            end if
            call yaml_sequence()
            call write_orbital_data(orbse%eval(iorbst+iorb),&
                 orbse%occup(iorbst+iorb),spinsignw,ikptw,mx,my,mz)

            if (nwrtmsg == 1) then
               !write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
               !   &   'evale(',iorb,')=',orbse%eval(iorb+iorbst),trim(message)
            else
               !if ((iorb <= 5 .or. iorb >= orbse%norbu-5) .or. verbose > 0) & 
               !write(*,'(1x,a,i0,a,1x,1pe21.14)') &
               !   &   'evale(',iorb,')=',orbse%eval(iorb+iorbst)
            end if
         else
            if (nwrtmsg==1) then
               if (abs(orbse%eval(iorb+iorbst)-preval(1)) <= etol .and. &
                  &   abs(orbse%eval(iorb+orbse%norbu+iorbst)-preval(2)) <= etol) then
               !degeneracy found
               message='  <-deg->  '
               !ndegen=ndegen+1 removed, only for non magnetized cases
            else if (abs(orbse%eval(iorb+iorbst)-preval(1)) <= etol) then
               !degeneracy found
               message='  <-deg    '
            else if (abs(orbse%eval(iorb+orbse%norbu+iorbst)-preval(2)) <= etol) then
               !degeneracy found
               message='    deg->  '
            else
               nwrtmsg=0
            end if
         end if
         if (iorb == norbu .and. iorb == norbd) then
            nwrtmsg=1
            message='  <-Last-> ' 
            preval(1)=orbse%eval(iorb+iorbst)
            preval(2)=orbse%eval(iorb+orbse%norbu+iorbst)
         else if (iorb == norbu) then
            nwrtmsg=1
            message='  <-Last   '
            preval(1)=orbse%eval(iorb+iorbst)
         else if (iorb == norbd) then
            nwrtmsg=1
            message='    Last-> '
            preval(2)=orbse%eval(iorb+orbse%norbu+iorbst)
         end if
         if ((iorb <= 5 .or. iorb >= orbse%norbu-5) .or. verbose > 0) then
            call yaml_sequence()
            call write_orbital_data(orbse%eval(iorb+iorbst),&
                 orbse%occup(iorb+iorbst),&
                 1.0_gp,ikptw,mx,my,mz)
            call yaml_sequence()
            call write_orbital_data(orbse%eval(iorb+iorbst+orbse%norbu),&
                 orbse%occup(iorb+iorbst+orbse%norbu),&
                 -1.0_gp,ikptw,mx,my,mz)
         end if
         
         if (nwrtmsg==1) then
            !write(*,'(1x,a,i4,a,1x,1pe21.14,a12,a,i4,a,1x,1pe21.14)') &
            !   &   'evale(',iorb,',u)=',orbse%eval(iorb+iorbst),message,&
            !   &   'evale(',iorb,',d)=',orbse%eval(iorb+orbse%norbu+iorbst)
         else
            !if ((iorb <= 5 .or. iorb >= orbse%norbu-5) .or. verbose > 0) & 
            !write(*,'(1x,a,i4,a,1x,1pe21.14,12x,a,i4,a,1x,1pe21.14)') &
            !   &   'evale(',iorb,',u)=',orbse%eval(iorb+iorbst),&
            !   &   'evale(',iorb,',d)=',orbse%eval(iorb+orbse%norbu+iorbst)
         end if
      end if
      if (iorb == orbse%norbu .and. ikpt == orbse%nkpts) then
         call yaml_sequence_close(advance='no')
      end if
      if (nwrtmsg==1) then
         call yaml_comment(adjustl(message))
      else
         call yaml_newline()
      end if
   end do
   !increment k-points shift
   iorbst=iorbst+orbse%norb
end do

if (orbse%efermi /= uninitialized(orbse%efermi)) then
   call yaml_map('Fermi Energy',orbse%efermi,fmt='(1pe21.11)')
end if

!call yaml_stream_attributes()

END SUBROUTINE write_ig_eigenvectors


!> Write orbital information with NO advance
subroutine write_orbital_data(eval,occup,spinsign,ikpt,mx,my,mz)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: ikpt !< k-point id 
  real(gp), intent(in) :: eval !< orbital energy
  real(gp), intent(in) :: occup !< orbital occupation number
  real(gp), intent(in) :: spinsign !< orbital spin (collinear and averaged)
  real(gp), intent(in) :: mx,my,mz !< spin magnetisation directions
  !local variables
  logical :: smallfmt


  call yaml_mapping_open(flow=.true.)
  !change the format if the spin and the k-point are initialized at the same way
  smallfmt=ikpt /= UNINITIALIZED(ikpt) .and. spinsign /= UNINITIALIZED(spinsign)

  if (smallfmt) then
     call yaml_map('e',eval,fmt='(1pe13.6)')
  else
     call yaml_map('e',eval,fmt='(1pe19.12)')
  end if

  if (occup /= UNINITIALIZED(occup)) then
     if (smallfmt) then
        call yaml_map('f',occup,fmt='(f5.3)')
     else
        call yaml_map('f',occup,fmt='(f6.4)')
     end if
  end if
  if (spinsign /= UNINITIALIZED(spinsign)) call yaml_map('s',int(spinsign),fmt='(i2)')
  if (ikpt /= UNINITIALIZED(ikpt)) then
     !if (int(spinsign)==-1) then
     !   call yaml_stream_attributes()
     !end if
     call yaml_map('k',ikpt,fmt='(i5)')
  end if
  if (mx /= UNINITIALIZED(mx) .and. my /= UNINITIALIZED(my) .and. mz /= UNINITIALIZED(mz)) &
     call yaml_map('M',(/mx,my,mz/),fmt='(f8.5)')
  call yaml_mapping_close(advance='no')
 
END SUBROUTINE write_orbital_data


!> Write DIIS weights
subroutine write_diis_weights(ncplx,idsx,ngroup,nkpts,itdiis,rds)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: ncplx,idsx,ngroup,nkpts,itdiis
  real(tp), dimension(ncplx,idsx+1,ngroup,nkpts), intent(in) :: rds
  !local variables
  integer :: j,igroup,ikpt
  character(len=2) :: mesupdw
  if (verbose < 10) then  
     !we restrict the printing to the first k point only.
     if (ngroup==1) then
        if (verbose >0) then
!!$           write(*,'(1x,a,2x,18(1x,1pe9.2))')&
!!$             'DIIS wgts:',reshape(rds(1:ncplx,1:itdiis+1,1,1),&
!!$             (/ncplx*(itdiis+1)/))!,&
           !yaml output
           call yaml_newline()
           call yaml_sequence_open('DIIS weights',flow=.true.)
           do j=1,itdiis+1
              call yaml_sequence(yaml_toa(rds(1,j,1,1),fmt='(1pe9.2)'))
           end do
           call yaml_sequence_close()
           !call yaml_map('DIIS weights',&
           !     (/(rds(1:ncplx,j,1,1),j=1,itdiis+1)/),fmt='(1pe9.2)')
        end if
!        write(70,'(1x,a,1pe9.2)',advance='no')'DIIS wgts: [ ',rds(1:ncplx,1,1,1)
        do j=2,itdiis+1
!           write(70,'(a,1pe9.2)',advance='no')', ',rds(1:ncplx,j,1,1)
        end do
!        write(70,'(a)')']'
        !'(',ttr,tti,')'
     else if (verbose >0) then
        do igroup=1,ngroup
           if (igroup==1) mesupdw='up'
           if (igroup==2) mesupdw='dw'
           write(*,'(1x,a,2x,18(1x,1pe9.2))')'DIIS wgts'//mesupdw//':',&
                (rds(1:ncplx,j,igroup,1),j=1,itdiis+1)
        end do
     end if
  else if (verbose >0) then
     do ikpt = 1, nkpts
        if (ngroup==1) then
           write(*,'(1x,a,I3.3,a,2x,9(1x,(1pe9.2)))')'DIIS wgts (kpt #', ikpt, &
                & ')',(rds(1:ncplx,j,1,ikpt),j=1,itdiis+1)
        else
           do igroup=1,ngroup
              if (igroup==1) mesupdw='up'
              if (igroup==2) mesupdw='dw'
              write(*,'(1x,a,I3.3,a,2x,9(1x,a,2(1pe9.2),a))')'DIIS wgts (kpt #', ikpt, &
                   & ')'//mesupdw//':',('(',rds(1:ncplx,j,igroup,ikpt),')',j=1,itdiis+1)
           end do
        end if
     end do
  end if
END SUBROUTINE write_diis_weights

!> Print gnrms (residue per orbital)
subroutine write_gnrms(nkpts,norb,gnrms)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: norb,nkpts
  real(wp), dimension(norb,nkpts), intent(in) :: gnrms
  !local variables
  integer :: ikpt,iorb

  call yaml_newline()
  call yaml_sequence_open('Residues per orbital',flow=.true.)
  call yaml_newline()

  do ikpt=1,nkpts
     if (nkpts > 1) call yaml_comment('Kpt #'//adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))))
     do iorb=1,norb
        call yaml_sequence(trim(yaml_toa(gnrms(iorb,ikpt),fmt='(1pe19.12)')),advance='no')
        if (ikpt == nkpts .and. iorb == norb)   call yaml_sequence_close(advance='no')
        call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')))
     end do
  end do
  
END SUBROUTINE write_gnrms


!> Print the atomic forces
subroutine write_forces(atoms,fxyz)
   use module_base
   use module_types
   use yaml_output
   implicit none
   !Arguments
   type(atoms_data), intent(in) :: atoms                !< Atoms data
   real(gp), dimension(3,atoms%astruct%nat), intent(in) :: fxyz !< Atomic forces
   !Local variables
   real(gp) :: sumx,sumy,sumz
   integer :: iat

   sumx=0.d0
   sumy=0.d0
   sumz=0.d0
   call yaml_comment('Atomic Forces',hfill='-')
   call yaml_sequence_open('Atomic Forces (Ha/Bohr)')
   do iat=1,atoms%astruct%nat
      call yaml_sequence(advance='no')
      call yaml_mapping_open(flow=.true.)
      call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),fxyz(1:3,iat),fmt='(1pe20.12)')
      !call yaml_map('AU',fxyz(1:3,iat),fmt='(1pe20.12)')
      !call yaml_map('eV/A',fxyz(1:3,iat)*Ha_eV/Bohr_Ang,fmt='(1pe9.2)')
      call yaml_mapping_close(advance='no')
      call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
!      write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
!      iat,trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),(fxyz(j,iat),j=1,3)
      sumx=sumx+fxyz(1,iat)
      sumy=sumy+fxyz(2,iat)
      sumz=sumz+fxyz(3,iat)
   end do
   call yaml_sequence_close()
   !$$        if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
   !$$           write(*,'(1x,a)')'the sum of the forces is'
   !$$           write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
   !$$           write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
   !$$           write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
   !$$        end if
END SUBROUTINE write_forces


!> Write stress tensor matrix
subroutine write_strten_info(fullinfo,strten,volume,pressure,message)
  use module_base
  use yaml_output
  implicit none
  logical, intent(in) :: fullinfo
  real(gp), intent(in) :: volume,pressure
  character(len=*), intent(in) :: message
  real(gp), dimension(6), intent(in) :: strten
  !local variables
  
  call yaml_sequence_open(trim(message)//' stress tensor matrix (Ha/Bohr^3)')
  call yaml_sequence(yaml_toa((/strten(1),strten(6),strten(5)/),fmt='(1pg20.12)'))
  call yaml_sequence(yaml_toa((/strten(6),strten(2),strten(4)/),fmt='(1pg20.12)'))
  call yaml_sequence(yaml_toa((/strten(5),strten(4),strten(3)/),fmt='(1pg20.12)'))
  call yaml_sequence_close()
  !write(*,'(1x,a)')'Stress Tensor, '//trim(message)//' contribution (Ha/Bohr^3):'
  !write(*,'(1x,t10,10x,a,t30,10x,a,t50,10x,a)')'x','y','z'
  !write(*,'(1x,a,t10,1pe20.12,t30,1pe20.12,t50,1pe20.12)')'x',strten(1),strten(6),strten(5)
  !write(*,'(1x,a,t30,1pe20.12,t50,1pe20.12)')'y',strten(2),strten(4)
  !write(*,'(1x,a,t50,1pe20.12)')'z',strten(3)

  if (fullinfo) then
     call yaml_mapping_open('Pressure')
     call yaml_map('Ha/Bohr^3',pressure,fmt='(1pg22.14)')
     call yaml_map('GPa',pressure*AU_GPa,fmt='(1pg14.6)')
     call yaml_map('PV (Ha)',pressure*volume,fmt='(1pg22.14)')
     call yaml_mapping_close()
     !write(*,'(1x,a,1pe22.14,a,1pe14.6,a,1pe22.14)')'Pressure:',pressure,&
     !     ' (',pressure*AU_GPa,' GPa), P V:',pressure*volume
  end if

END SUBROUTINE write_strten_info


!> Assign some of the physical system variables
!! Performs also some cross-checks with other variables
!! The pointers in atoms structure have to be associated or nullified.
subroutine print_atomic_variables(atoms, hmax, ixc, dispersion)
  use module_base
  use module_types
  use module_atoms, only: RADII_SOURCE
  use module_xc
  use vdwcorrection
  use yaml_output
  use psp_projectors, only: PSPCODE_HGH,PSPCODE_HGH_K,PSPCODE_HGH_K_NLCC,&
       PSPCODE_PAW,PSPCODE_GTH
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: hmax
  integer, intent(in) :: ixc, dispersion
  !Local variables
  logical :: nonloc
  integer :: i,j,l,ityp,iat,natyp,mproj,inlcc
  real(gp) :: minrad
  real(gp), dimension(3,3) :: hij
  real(gp), dimension(2,2,3) :: offdiagarr
  character(len=500) :: name_xc1, name_xc2

  !If no atoms...
  if (atoms%astruct%ntypes == 0) return

  !print the pseudopotential matrices
  do l=1,3
     do i=1,2
        do j=i+1,3
           offdiagarr(i,j-i,l)=0._gp
           if (l==1) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(3._gp/5._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(5._gp/21._gp)
              else
                 offdiagarr(i,j-i,l)=-0.5_gp*sqrt(100._gp/63._gp)
              end if
           else if (l==2) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(5._gp/7._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=1._gp/6._gp*sqrt(35._gp/11._gp)
              else
                 offdiagarr(i,j-i,l)=-7._gp/3._gp*sqrt(1._gp/11._gp)
              end if
           else if (l==3) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(7._gp/9._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(63._gp/143._gp)
              else
                 offdiagarr(i,j-i,l)=-9._gp*sqrt(1._gp/143._gp)
              end if
           end if
        end do
     end do
  end do

!  write(*,'(1x,a)')&
  !       '------------------------------------ Pseudopotential coefficients (Upper Triangular)'
  inlcc=0
  call yaml_comment('System Properties',hfill='-')
  call yaml_sequence_open('Properties of atoms in the system')
  do ityp=1,atoms%astruct%ntypes
     call yaml_sequence(advance='no')
     call yaml_map('Symbol',trim(atoms%astruct%atomnames(ityp)),advance='no')
     call yaml_comment('Type No. '//trim(yaml_toa(ityp,fmt='(i2.2)')))
     call yaml_map('No. of Electrons',atoms%nelpsp(ityp))
     natyp=0
     do iat=1,atoms%astruct%nat
        if (atoms%astruct%iatype(iat) == ityp) natyp=natyp+1
     end do
     call yaml_map('No. of Atoms',natyp)

     call yaml_mapping_open('Radii of active regions (AU)')!,flow=.true.)
       call yaml_map('Coarse',atoms%radii_cf(ityp,1),fmt='(f8.5)')
       call yaml_map('Fine',atoms%radii_cf(ityp,2),fmt='(f8.5)')
       call yaml_map('Coarse PSP',atoms%radii_cf(ityp,3),fmt='(f8.5)')
       call yaml_map('Source',RADII_SOURCE(atoms%iradii_source(ityp)))
       !if (atoms%radii_cf(ityp, 1) == UNINITIALIZED(1.0_gp)) then
       !   call yaml_map('Source','Hard-Coded')
       !else
       !   call yaml_map('Source','PSP File')
       !end if
     call yaml_mapping_close()

     minrad=1.e10_gp
     do i=0,4
        if (atoms%psppar(i,0,ityp)/=0._gp) then
           minrad=min(minrad,atoms%psppar(i,0,ityp))
        end if
     end do
     if (atoms%radii_cf(ityp,2) /=0.0_gp) then
        call yaml_map('Grid Spacing threshold (AU)',2.5_gp*minrad,fmt='(f5.2)')
     else
        call yaml_map('Grid Spacing threshold (AU)',1.25_gp*minrad,fmt='(f5.2)')
     end if
     !control whether the grid spacing is too high
     if (hmax > 2.5_gp*minrad) then
        call yaml_warning('Chosen Grid spacings seem too high for the '// &
           & trim(atoms%astruct%atomnames(ityp))//' atom type. At you own risk!')
     end if

     select case(atoms%npspcode(ityp))
     case(PSPCODE_GTH)
        call yaml_map('Pseudopotential type','GTH')
     case(PSPCODE_HGH)
        call yaml_map('Pseudopotential type','HGH')
     case(PSPCODE_HGH_K)
        call yaml_map('Pseudopotential type','HGH-K')
     case(PSPCODE_HGH_K_NLCC)
        call yaml_map('Pseudopotential type','HGH-K + NLCC')
     case(PSPCODE_PAW)
        call yaml_map('Pseudopotential type','PAW + HGH')
     end select
     if (atoms%psppar(0,0,ityp)/=0) then
        call yaml_mapping_open('Local Pseudo Potential (HGH convention)')
          call yaml_map('Rloc',atoms%psppar(0,0,ityp),fmt='(f9.5)')
          call yaml_map('Coefficients (c1 .. c4)',atoms%psppar(0,1:4,ityp),fmt='(f9.5)')
        call yaml_mapping_close()
     end if
     !nlcc term
     if (atoms%npspcode(ityp) == PSPCODE_HGH_K_NLCC) then
        inlcc=inlcc+1
        call yaml_mapping_open('Non Linear Core Correction term')
            call yaml_map('Rcore',atoms%nlccpar(0,inlcc),fmt='(f9.5)')
            call yaml_map('Core charge',atoms%nlccpar(1,inlcc),fmt='(f9.5)')
        call yaml_mapping_close()
     end if
     !see if nonlocal terms are present
     nonloc=.false.
     verify_nl: do l=1,3
        do i=3,0,-1
           j=i
           if (atoms%psppar(l,i,ityp) /= 0._gp) exit
        end do
        if (j /=0) then
           nonloc=.true.
           exit verify_nl
        end if
     end do verify_nl
     if (nonloc) then
        call yaml_sequence_open('NonLocal PSP Parameters')
        do l=1,3
           do i=3,0,-1
              j=i
              if (atoms%psppar(l,i,ityp) /= 0._gp) exit
           end do
           if (j /=0) then
              call yaml_sequence(advance='no')
              call yaml_map('Channel (l)',l-1)
              call yaml_map('Rloc',atoms%psppar(l,0,ityp),fmt='(f9.5)')
              hij=0._gp
              do i=1,j
                 hij(i,i)=atoms%psppar(l,i,ityp)
              end do
              if (atoms%npspcode(ityp) == PSPCODE_HGH) then !traditional HGH convention
                 hij(1,2)=offdiagarr(1,1,l)*atoms%psppar(l,2,ityp)
                 hij(1,3)=offdiagarr(1,2,l)*atoms%psppar(l,3,ityp)
                 hij(2,3)=offdiagarr(2,1,l)*atoms%psppar(l,3,ityp)
              else if (atoms%npspcode(ityp) == PSPCODE_HGH_K &
                  .or. atoms%npspcode(ityp) == PSPCODE_HGH_K_NLCC) then !HGH-K convention
                 hij(1,2)=atoms%psppar(l,4,ityp)
                 hij(1,3)=atoms%psppar(l,5,ityp)
                 hij(2,3)=atoms%psppar(l,6,ityp)
              end if
              call yaml_sequence_open('h_ij matrix')
                call yaml_sequence(trim(yaml_toa(hij(1,1:3),fmt='(f9.5)')))
                call yaml_sequence(trim(yaml_toa((/hij(1,2),hij(2,2),hij(2,3)/),fmt='(f9.5)')))
                call yaml_sequence(trim(yaml_toa((/hij(1,3),hij(2,3),hij(3,3)/),fmt='(f9.5)')))
              call yaml_sequence_close()
           end if
        end do
        call yaml_sequence_close()
     end if
     ! PAW case.
     if (atoms%npspcode(ityp) == PSPCODE_PAW) then
        call yaml_map('No. of gaussians', atoms%pawtab(ityp)%wvl%pngau)
        call yaml_map('complex coefficients (1..5)', atoms%pawtab(ityp)%wvl%parg(:,1:5))
        call yaml_map('complex factors (1..5)', atoms%pawtab(ityp)%wvl%pfac(:,1:5))
     end if
     mproj = 0
     do l=1,4 
        do i=1,3 
           if (atoms%psppar(l,i,ityp) /= 0.0_gp) mproj=mproj+2*l-1
        enddo
     enddo
     !call numb_proj(ityp,atoms%astruct%ntypes,atoms%psppar,atoms%npspcode,mproj)
     call yaml_map('No. of projectors',mproj)

     !control if the PSP is calculated with the same XC value
     if (atoms%ixcpsp(ityp) < 0) then
        call xc_get_name(name_xc1, atoms%ixcpsp(ityp), XC_MIXED)
     else
        call xc_get_name(name_xc1, atoms%ixcpsp(ityp), XC_ABINIT)
     end if
     if (ixc < 0) then
        call xc_get_name(name_xc2, ixc, XC_MIXED)
     else
        call xc_get_name(name_xc2, ixc, XC_ABINIT)
     end if
     call yaml_map('PSP XC','"'//trim(name_xc1)//'"')
     if (trim(name_xc1) /= trim(name_xc2)) then
        call yaml_warning('PSP generated with a different XC. Input XC is "'//trim(name_xc2) // '"')
     end if
  end do
  call yaml_sequence_close()
!!!  tt=dble(norb)/dble(nproc)
!!!  norbp=int((1.d0-eps_mach*tt) + tt)
!!!  !if (verb.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp

  call vdwcorrection_warnings(atoms, dispersion, ixc)

  ! if linear scaling applied with more then InputGuess, then go read input.lin for radii
  !  if (in%linear /= 'OFF' .and. in%linear /= 'LIG') then
  !     lin%nlr=atoms%astruct%nat
  !     call allocateBasicArrays(atoms, lin)
  !     call readLinearParameters(verb, nproc, lin, atoms, atomNames)
  !  end if
END SUBROUTINE print_atomic_variables


!> Display an estimation of the occupied memory
subroutine print_memory_estimation(mem)
  use module_types
  use yaml_output
  implicit none
  type(memory_estimation), intent(in) :: mem

  call yaml_comment('Estimation of Memory Consumption',hfill='-')
  call yaml_mapping_open('Memory requirements for principal quantities (MiB.KiB)')
    call yaml_map('Subspace Matrix',trim(MibdotKib(mem%submat)),advance='no')
      call yaml_comment('(Number of Orbitals:'//trim(yaml_toa(mem%norb))//')',tabbing=50)
    call yaml_map('Single orbital',trim(MibdotKib(mem%oneorb)),advance='no')
      call yaml_comment('(Number of Components:'//trim(yaml_toa(mem%ncomponents))//')',tabbing=50)
    call yaml_map('All (distributed) orbitals',trim(MibdotKib(mem%allpsi_mpi)),advance='no')
      call yaml_comment('(Number of Orbitals per MPI task:'//trim(yaml_toa(mem%norbp))//')',tabbing=50)
    call yaml_map('Wavefunction storage size',trim(MibdotKib(mem%psistorage)),advance='no')
      call yaml_comment('(DIIS/SD workspaces included)',tabbing=50)
    call yaml_map('Nonlocal Pseudopotential Arrays',trim(MibdotKib(mem%projarr)))
    call yaml_map('Full Uncompressed (ISF) grid',trim(MibdotKib(mem%grid)))
    call yaml_map('Workspaces storage size',trim(MibdotKib(mem%workarr)))
  call yaml_mapping_close()

  call yaml_mapping_open('Accumulated memory requirements during principal run stages (MiB.KiB)')
     call yaml_map('Kernel calculation',trim(MibdotKib(mem%kernel)))
     call yaml_map('Density Construction',trim(MibdotKib(mem%density)))
     call yaml_map('Poisson Solver',trim(MibdotKib(mem%psolver)))
     call yaml_map('Hamiltonian application',trim(MibdotKib(mem%ham)))
     call yaml_map('Orbitals Orthonormalization',trim(MibdotKib(mem%ham+mem%submat)))
!           call yaml_comment('Wfn, Work, Den, Ker ',tabbing=50)
  call yaml_mapping_close()
  call yaml_map('Estimated Memory Peak (MB)',yaml_toa(mega(mem%peak)))

contains

  function mega(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    integer(kind=8) :: mega
    mega=int(omemory/1048576.d0,kind=8)
  end function mega

  function kappa(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    integer :: kappa
    kappa=ceiling((omemory-aint(omemory/1048576.d0)*1048576.d0)/1024.d0)
  end function kappa

  function MiBdotKiB(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    character(len=50) MiBdotKiB

    MiBdotKiB=repeat(' ',len(MiBdotKiB))

    MiBdotKiB=trim(adjustl(yaml_toa(int(mega(omemory)))))//'.'//&
         trim(adjustl(yaml_toa(int(kappa(omemory)))))
    
  end function MiBdotKiB

END SUBROUTINE print_memory_estimation


!> Display information about the box and the grid
subroutine print_atoms_and_grid(Glr, atoms, rxyz, shift, hx, hy, hz)
  use module_defs
  use module_types
  use yaml_output
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3), intent(in) :: shift
  real(gp), intent(in) :: hx, hy, hz
  !Local variables
  integer :: iat, iunit

  if (atoms%astruct%ntypes > 0) then
     call yaml_comment('Atom Positions',hfill='-')
     call yaml_sequence_open('Atomic positions within the cell (Atomic and Grid Units)')
     do iat=1,atoms%astruct%nat
        call yaml_sequence(advance='no')
        call yaml_mapping_open(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),flow=.true.)
        call yaml_map('AU',rxyz(1:3,iat),fmt='(1pg12.5)')
        call yaml_map('GU',(/rxyz(1,iat)/hx,rxyz(2,iat)/hy,rxyz(3,iat)/hz/),fmt='(1pg12.5)')
        call yaml_mapping_close(advance='no')
        call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
     enddo
     call yaml_sequence_close()
     call yaml_map('Rigid Shift Applied (AU)',(/-shift(1),-shift(2),-shift(3)/),fmt='(1pg12.5)')
     ! New version
     call yaml_mapping_open('Atomic structure')
     call yaml_get_default_stream(unit = iunit)
     call wtyaml(iunit, UNINITIALIZED(1.d0), rxyz, atoms%astruct, .false., rxyz, &
          .true., shift, (/ hx, hy, hz /))
     call yaml_mapping_close()
  end if
  call yaml_comment('Grid properties',hfill='-')
  call yaml_map('Box Grid spacings',(/hx,hy,hz/),fmt='(f7.4)')
  call yaml_mapping_open('Sizes of the simulation domain')
  call yaml_map('AU',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3)/),fmt='(1pg12.5)')
  call yaml_map('Angstroem',(/atoms%astruct%cell_dim(1)*Bohr_Ang,&
       atoms%astruct%cell_dim(2)*Bohr_Ang,atoms%astruct%cell_dim(3)*Bohr_Ang/),fmt='(1pg12.5)')
  call yaml_map('Grid Spacing Units',(/Glr%d%n1,Glr%d%n2,Glr%d%n3/),fmt='(i4)')
  call yaml_mapping_open('High resolution region boundaries (GU)',flow=.false.)
  call yaml_map('From',(/Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3/),fmt='(i4)')
  call yaml_map('To',(/Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3/),fmt='(i4)')
  call yaml_mapping_close()
  call yaml_mapping_close()
  call yaml_map('High Res. box is treated separately',Glr%hybrid_on)
END SUBROUTINE print_atoms_and_grid

!> Write atomic file in yaml format
subroutine wtyaml(iunit,energy,rxyz,astruct,wrtforces,forces, &
     & wrtlog, shift, hgrids)
  use module_base, only: f_err_throw
  use module_defs, only: Bohr_Ang, gp, UNINITIALIZED
  use yaml_output
  use module_atoms, only: atomic_structure,frozen_itof
  use ao_inguess, only: charge_and_spol
  implicit none
  !Arguments
  logical, intent(in) :: wrtforces !< True if write the atomic forces
  logical, intent(in) :: wrtlog
  integer, intent(in) :: iunit
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(in) :: energy
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz,forces
  real(gp), dimension(3), intent(in) :: shift, hgrids
  !local variables
  logical :: reduced, perx, pery, perz
  character(len=4) :: frzchain
  real(gp), dimension(3) :: xred
  real(gp) :: factor
  integer :: iat,ichg,ispol

  reduced=.false.
  Units: select case(trim(astruct%units))
  case('angstroem','angstroemd0')
     call yaml_map('Units','angstroem', unit = iunit)
     factor=Bohr_Ang
  case('reduced')
     if (.not. wrtlog) then
        call yaml_map('Units','reduced', unit = iunit)
        reduced=.true.
     end if
     factor = 1.0_gp
  case('atomic','atomicd0','bohr','bohrd0')
     ! Default
     factor=1.0_gp
     !call yaml_map('Units','bohr')
  case default
     call f_err_throw('Writing the atomic file. Error, unknown units ("'// trim(astruct%units)//'")', & 
          & err_name='BIGDFT_RUNTIME_ERROR')
  end select Units

  !cell information
  perx = .false.
  pery = .false.
  perz = .false.
  BC :select case(astruct%geocode)
  case('S')
     call yaml_sequence_open('Cell', flow=.true., unit = iunit)
     call yaml_sequence(yaml_toa(astruct%cell_dim(1)*factor), unit = iunit) !x
     call yaml_sequence('.inf', unit = iunit)             !y
     call yaml_sequence(yaml_toa(astruct%cell_dim(3)*factor), unit = iunit) !z
     call yaml_sequence_close(unit = iunit)
     !angdeg to be added
     perx = .true.
     pery = .false.
     perz = .true.
  case('W')
     call yaml_sequence_open('Cell', flow=.true., unit = iunit)
     call yaml_sequence('.inf', unit = iunit)             !x
     call yaml_sequence('.inf', unit = iunit)             !y
     call yaml_sequence(yaml_toa(astruct%cell_dim(3)*factor), unit = iunit) !z
     call yaml_sequence_close(unit = iunit)
     perx = .false.
     pery = .false.
     perz = .true.
  case('P')
     call yaml_map('Cell',(/astruct%cell_dim(1)*factor, &
          & astruct%cell_dim(2)*factor, astruct%cell_dim(3)*factor/), unit = iunit)
     !angdeg to be added
     perx = .true.
     pery = .true.
     perz = .true.
  case('F')
     ! Default
     !call yaml_map('BC','free')
  end select BC

  !Write atomic positions
  call yaml_sequence_open('Positions', unit = iunit)
  do iat=1,astruct%nat
     call yaml_sequence(advance='no', unit = iunit)
     if (extra_info(iat)) then
        call yaml_mapping_open(flow=.true., unit = iunit)
     end if
     xred=rxyz(:,iat)
     if (reduced) then
        if (perx) xred(1)=rxyz(1,iat)/astruct%cell_dim(1)
        if (pery) xred(2)=rxyz(2,iat)/astruct%cell_dim(2)
        if (perz) xred(3)=rxyz(3,iat)/astruct%cell_dim(3)
     else
        !Multiply by the factor to have the right units
        xred = xred*factor
     end if
     if (wrtlog) then
        call print_one_atom(trim(astruct%atomnames(astruct%iatype(iat))),&
             xred,hgrids,iat)
!!$        call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),&
!!$             & xred,fmt="(g18.10)", unit = iunit, advance = "no")
!!$        xred(1:3) = rxyz(1:3,iat) / hgrids
!!$        write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, iat
!!$        call yaml_comment(gu, unit = iunit)
     else
        call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),&
             & xred,fmt="(1pg25.17)", unit = iunit)
     end if
     if (extra_info(iat)) then
        call charge_and_spol(astruct%input_polarization(iat),ichg,ispol)
        if (ispol /=0) call yaml_map('IGSpin',ispol, unit = iunit)
        if (ichg /=0) call yaml_map('IGChg',ichg, unit = iunit)
        if (astruct%ifrztyp(iat) /= 0) then
           call frozen_itof(astruct%ifrztyp(iat),frzchain)
           call yaml_map('Frozen',frzchain, unit = iunit)
        end if
        call yaml_mapping_close(unit = iunit)
     end if
  end do
  call yaml_sequence_close(unit = iunit) !positions

  !Write atomic forces
  if (wrtforces) then
     call yaml_sequence_open('Forces (Ha/Bohr)', unit = iunit)
     do iat=1,astruct%nat
        call yaml_sequence(advance='no', unit = iunit)
        call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),forces(:,iat),fmt='(1pg25.17)', unit = iunit)
     end do
     call yaml_sequence_close(unit = iunit) !values
  end if
  if (wrtlog) then
     call yaml_map('Rigid Shift Applied (AU)',(/-shift(1),-shift(2),-shift(3)/),fmt='(1pg12.5)')
  else
     call yaml_mapping_open('Properties', unit = iunit)
     call yaml_map('Timestamp',yaml_date_and_time_toa(), unit = iunit)
     if (energy /= 0.0_gp .and. energy /= UNINITIALIZED(energy)) then
        call yaml_map("Energy (Ha)", energy, unit = iunit)
     end if
     call yaml_mapping_close(unit = iunit) !properties
  end if

contains

  function extra_info(iat)
    implicit none
    integer, intent(in) :: iat
    logical extra_info
    extra_info=astruct%input_polarization(iat) /=100 .or. astruct%ifrztyp(iat)/=0
  end function extra_info

  subroutine print_one_atom(atomname,rxyz,hgrids,id)
    implicit none
    integer, intent(in) :: id
    character(len=*), intent(in) :: atomname
    double precision, dimension(3), intent(in) :: rxyz,hgrids
    !local variables
    character(len=*), parameter :: fmtat='(1pg18.10)',fmtg='(F7.2)',fmti='(i4.4)'
    integer :: i

    call yaml_sequence_open(atomname,flow=.true.)
    do i=1,3
       call yaml_sequence(yaml_toa(rxyz(i),fmt=fmtat))
    end do
    call yaml_sequence_close(advance='no')
    call yaml_comment(trim(yaml_toa(rxyz/hgrids/factor,fmt=fmtg))//trim(yaml_toa(id,fmt=fmti))) !we can also put tabbing=

  end subroutine print_one_atom

END SUBROUTINE wtyaml

subroutine print_wfd(wfd)
  use module_types, only: wavefunctions_descriptors
  use yaml_output
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd

  call yaml_mapping_open('Wavefunctions Descriptors, full simulation domain')
  !write(*,'(1x,a)')&
  !   &   '------------------------------------------------- Wavefunctions Descriptors Creation'

  !write(*,'(2(1x,a,i10))') &
  !     &   'Coarse resolution grid: Number of segments= ',Glr%wfd%nseg_c,'points=',Glr%wfd%nvctr_c
  call yaml_mapping_open('Coarse resolution grid')!,flow=.true.)
  call yaml_map('No. of segments',wfd%nseg_c)
  call yaml_map('No. of points',wfd%nvctr_c)
  call yaml_mapping_close()

  !write(*,'(2(1x,a,i10))')
  !'  Fine resolution grid: Number of segments= ',Glr%wfd%nseg_f,'points=',Glr%wfd%nvctr_f
  call yaml_mapping_open('Fine resolution grid')!,flow=.true.)
  call yaml_map('No. of segments',wfd%nseg_f)
  call yaml_map('No. of points',wfd%nvctr_f)
  call yaml_mapping_close()
  call yaml_mapping_close()
  
END SUBROUTINE print_wfd

subroutine print_nlpsp(nlpsp)
  use module_defs, only: gp
  use module_types, only: DFT_PSP_projectors
  use yaml_output
  implicit none
  type(DFT_PSP_projectors), intent(in) :: nlpsp
  !local variables
  integer :: iat,ilr,sizemask,maxmask,totmask,totpack

  call yaml_mapping_open('NonLocal PSP Projectors Descriptors')
  if (nlpsp%on_the_fly) then
     call yaml_map('Creation strategy','On-the-fly')
  else
     call yaml_map('Creation strategy','Once-and-for-all')
  end if
  call yaml_map('Total number of projectors',nlpsp%nproj)
  call yaml_map('Total number of components',nlpsp%nprojel)
  call yaml_map('Percent of zero components',nint(100.0_gp*nlpsp%zerovol))
  !calculate the amount of memory spent in the descriptor for the wavefunction
  maxmask=0
  totmask=0
  totpack=0
  do iat=1,nlpsp%natoms
     if (nlpsp%pspd(iat)%mproj>0) then
        totpack=max(totpack,nlpsp%pspd(iat)%plr%wfd%nvctr_c+&
             7*nlpsp%pspd(iat)%plr%wfd%nvctr_f)
     end if
     sizemask=0
     if (associated(nlpsp%pspd(iat)%tolr)) then
        !do ilr=1,nlpsp%pspd(iat)%nlr
        do ilr=1,size(nlpsp%pspd(iat)%tolr)
           sizemask=sizemask+&
                nlpsp%pspd(iat)%tolr(ilr)%nmseg_c+nlpsp%pspd(iat)%tolr(ilr)%nmseg_f
        end do
     end if
     maxmask=max(maxmask,sizemask)
     totmask=totmask+sizemask
  end do
  totpack=totpack*4
  if (associated(nlpsp%scpr)) totpack=totpack+size(nlpsp%scpr)
  if (associated(nlpsp%cproj)) totpack=totpack+size(nlpsp%cproj)*2
  if (totpack /=0) &
       call yaml_map('Size of workspaces',totpack)
  if (maxmask /=0) &
       call yaml_map('Maximum size of masking arrays for a projector',3*maxmask)
  if (totmask /=0) &
       call yaml_map('Cumulative size of masking arrays',3*totmask)

  call yaml_mapping_close()
END SUBROUTINE print_nlpsp


!> Display information about the electronic orbitals
subroutine print_orbitals(orbs, geocode)
  use module_types, only: orbitals_data
  use module_defs, only: gp
  use yaml_output
  implicit none
  type(orbitals_data), intent(in) :: orbs
  character(len = 1), intent(in) :: geocode
  
  integer :: jproc, nproc, jpst, norbme, norbyou, nelec
  integer :: ikpts, iorb1, iorb
  real(gp) :: rocc

  nelec = int(sum(orbs%occup) + 1d-12) / orbs%nkpts

  call yaml_comment('Occupation Numbers',hfill='-')
  call yaml_map('Total Number of Electrons',nelec,fmt='(i8)')

  ! Number of orbitals
  if (orbs%nspin==1) then
     call yaml_map('Spin treatment','Averaged')
     if (mod(nelec,2).ne.0) then
        call yaml_warning('Odd number of electrons, no closed shell system')
        !write(*,'(1x,a)') 'WARNING: odd number of electrons, no closed shell system'
     end if
  else if(orbs%nspin==4) then
     call yaml_map('Spin treatment','Spinorial (non-collinearity possible)')
  else 
     call yaml_map('Spin treatment','Collinear')
  end if

  !distribution of wavefunction arrays between processors
  !tuned for the moment only on the cubic distribution
  call yaml_mapping_open('Orbitals Repartition')
  jpst=0
  nproc = size(orbs%norb_par, 1)
  do jproc=0,nproc-1
     norbme=orbs%norb_par(jproc,0)
     norbyou=orbs%norb_par(min(jproc+1,nproc-1),0)
     if (norbme /= norbyou .or. jproc == nproc-1) then
        call yaml_map('MPI tasks '//trim(yaml_toa(jpst,fmt='(i0)'))//'-'//trim(yaml_toa(jproc,fmt='(i0)')),norbme,fmt='(i0)')
        !write(*,'(3(a,i0),a)')&
        !     ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' orbitals '
        jpst=jproc+1
     end if
  end do
  !write(*,'(3(a,i0),a)')&
  !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' orbitals '
  call yaml_mapping_close()
  
  call yaml_map('Total Number of Orbitals',orbs%norb,fmt='(i8)')

  !No orbs finished
  if (orbs%norb == 0) return

  call yaml_sequence_open('Input Occupation Numbers')
  do ikpts=1,orbs%nkpts
     if (geocode /= 'F') then
        call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpts,fmt='(i4.4)'))) // ' BZ coord. = ' // &
        & trim(yaml_toa(orbs%kpts(:, ikpts),fmt='(f12.6)')))
     end if
     call yaml_sequence(advance='no')
     call yaml_mapping_open('Occupation Numbers',flow=.true.)
     !write(*,'(1x,a,t28,i8)') 'Total Number of Orbitals',norb
     iorb1=1
     rocc=orbs%occup(1+(ikpts-1)*orbs%norb)
     do iorb=1,orbs%norb
        if (orbs%occup(iorb+(ikpts-1)*orbs%norb) /= rocc) then
           if (iorb1 == iorb-1) then
              call yaml_map('Orbital No.'//trim(yaml_toa(iorb1)),rocc,fmt='(f6.4)')
              !write(*,'(1x,a,i0,a,f6.4)') 'occup(',iorb1,')= ',rocc
           else
           call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
                adjustl(trim(yaml_toa(iorb-1))),rocc,fmt='(f6.4)')
           !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',iorb-1,')= ',rocc
           end if
           rocc=orbs%occup(iorb+(ikpts-1)*orbs%norb)
           iorb1=iorb
        end if
     enddo
     if (iorb1 == orbs%norb) then
        call yaml_map('Orbital No.'//trim(yaml_toa(orbs%norb)),orbs%occup(ikpts*orbs%norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,f6.4)') 'occup(',norb,')= ',occup(norb)
     else
        call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
             adjustl(trim(yaml_toa(orbs%norb))),orbs%occup(ikpts*orbs%norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',norb,')= ',occup(norb)
     end if
     call yaml_mapping_close()
  end do
  call yaml_sequence_close()
END SUBROUTINE print_orbitals
