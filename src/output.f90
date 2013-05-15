!>  @file
!!  File where most relavant screen output are collected
!!  Routines which are present in this file should have *all* arguments as intent(in)
!!  Also, the master process only should acces these routines
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Display the logo of BigDFT 
subroutine print_logo()
  use module_base
  use yaml_output
  implicit none
  integer :: namelen,ierr
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  integer :: nthreads
  integer, parameter :: ln = 1024
  character(len = ln), dimension(4) :: buf
!$ integer :: omp_get_max_threads

!  call yaml_comment('Daubechies Wavelets for DFT Pseudopotential Calculations',hfill='=')


  call yaml_open_map('Code logo')
  call yaml_scalar('      TTTT         F       DDDDD    ')     
  call yaml_scalar('     T    T               D         ')     
  call yaml_scalar('    T     T        F     D          ')     
  call yaml_scalar('    T    T         F     D        D ')     
  call yaml_scalar('    TTTTT          F     D         D')     
  call yaml_scalar('    T    T         F     D         D')     
  call yaml_scalar('    T     T        F     D         D')     
  call yaml_scalar('    T      T       F     D         D')     
  call yaml_scalar('    T     T     FFFF     D         D')     
  call yaml_scalar('    T TTTT         F      D        D')     
  call yaml_scalar('    T             F        D      D ')     
  call yaml_scalar('TTTTTTTTT    FFFFF          DDDDDD  ')     
  !call yaml_scalar()'-----------------------------------)----'
  call yaml_scalar('  gggggg          iiiii    BBBBBBBBB')     
  call yaml_scalar(' g      g        i             B    ')     
  call yaml_scalar('g        g      i         BBBB B    ')     
  call yaml_scalar('g         g     iiii     B     B    ')     
  call yaml_scalar('g         g     i       B      B    ')     
  call yaml_scalar('g         g     i        B     B    ')     
  call yaml_scalar('g         g     i         B    B    ')     
  call yaml_scalar('g         g     i          BBBBB    ')     
  call yaml_scalar(' g        g     i         B    B    ')     
  call yaml_scalar('          g     i        B     B    ')     
  call yaml_scalar('         g               B    B     ')     
  call yaml_scalar('    ggggg       i         BBBB      ') 
  call yaml_close_map()

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
  call yaml_open_map("Compilation options")
  call bigdft_config_get_user_args(buf(1), ln)
  call yaml_map("Configure arguments", '"'//trim(buf(1))//'"')
  call bigdft_config_get_compilers(buf(1), buf(2), buf(3), ln)
  call yaml_map("Compilers (CC, FC, CXX)", buf(1:3))
  call bigdft_config_get_compiler_flags(buf(1), buf(2), buf(3), buf(4), ln)
  call yaml_open_map("Compiler flags")
   call yaml_map("CFLAGS",   trim(buf(1)))
   call yaml_map("FCFLAGS",  trim(buf(2)))
   call yaml_map("CXXFLAGS", trim(buf(3)))
   call yaml_map("CPPFLAGS", trim(buf(4)))
  call yaml_close_map()
!!$  call bigdft_config_get_linker(buf(1), buf(2), buf(3), buf(4), ln)
!!$  call yaml_open_map("Linker")
!!$   call yaml_map("LD",      trim(buf(1)))
!!$   call yaml_map("LDFLAGS", trim(buf(2)))
!!$   call yaml_map("LIBS",    trim(buf(3)))
!!$   call yaml_map("Full linking options", trim(buf(4)))
!!$  call yaml_close_map()
 
 call yaml_close_map()

end subroutine print_configure_options


!> Print all general parameters
subroutine print_general_parameters(nproc,in,atoms)
  use module_base
  use module_types
  use defs_basis
  use m_ab6_symmetry
  use yaml_output
  use module_input, only: case_insensitive_equiv
  implicit none
  !Arguments
  integer, intent(in) :: nproc
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms

  integer :: nSym, ierr, iat, i
  integer :: sym(3, 3, AB6_MAX_SYMMETRIES)
  integer :: symAfm(AB6_MAX_SYMMETRIES)
  real(gp) :: transNon(3, AB6_MAX_SYMMETRIES)
  real(gp) :: genAfm(3)
  character(len=15) :: spaceGroup
  character(len=len(in%run_name)) :: prefix
  integer :: spaceGroupId, pointGroupMagn
  !integer, parameter :: maxLen = 50, width = 24
  !character(len = width), dimension(maxLen) :: at, fixed, add
  !integer :: ityp,lg
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
  call yaml_open_map('Atomic System Properties')
     call yaml_map('Number of atomic types', atoms%astruct%ntypes, fmt='(i0)')
     call yaml_map('Number of atoms', atoms%astruct%nat, fmt='(i0)')
     if (atoms%astruct%nat > 0) then
        call yaml_map('Types of atoms',atoms%astruct%atomnames)
        !call yaml_map('Types of atoms',flow=.true.)
        !do ityp=1,atoms%astruct%ntypes
        !   call yaml_sequence(trim(atoms%astruct%atomnames(ityp)))
        !end do
        ! Fixed positions
        if (maxval(atoms%astruct%ifrztyp) /= 0) then
           call yaml_open_sequence('Fixed atoms',flow=.true.)
           ! The fixed atom column
           do iat=1,atoms%astruct%nat
              if (atoms%astruct%ifrztyp(iat) /= 0) then
                 call yaml_sequence('at.' // trim(yaml_toa(iat,fmt='(i4.4)')) // &
                      & '(' // trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))) // ')' &
                      & // trim(yaml_toa(atoms%astruct%ifrztyp(iat),fmt='(i0)')))
              end if
           end do
           call yaml_close_sequence()
        end if
     end if
     !Boundary Conditions
     !call yaml_map('Geometry Code',trim(atoms%astruct%geocode))
     if (atoms%astruct%geocode == 'P') then
        call yaml_map('Boundary Conditions','Periodic',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
        call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
             atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
     else if (atoms%astruct%geocode ==  'S') then
        call yaml_map('Boundary Conditions','Surface',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
        call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
             atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
     else if (atoms%astruct%geocode == 'F') then
        call yaml_map('Boundary Conditions','Free',advance='no')
        call yaml_comment('Code: '//atoms%astruct%geocode)
     end if
     !Symmetries
     if (atoms%astruct%geocode /= 'F' .and. .not. in%disableSym) then
        call symmetry_get_matrices(atoms%astruct%sym%symObj, nSym, sym, transNon, symAfm, ierr)
        call symmetry_get_group(atoms%astruct%sym%symObj, spaceGroup, &
             & spaceGroupId, pointGroupMagn, genAfm, ierr)
        if (ierr == AB6_ERROR_SYM_NOT_PRIMITIVE) write(spaceGroup, "(A)") "not prim."
     else 
        nSym = 0
        spaceGroup = 'disabled'
     end if
     call yaml_map('Number of Symmetries',nSym)
     call yaml_map('Space group',trim(spaceGroup))
  call yaml_close_map()

  !write(*,'(1x,a,a,a)') '--- (file: posinp.', &
  !     & atoms%astruct%inputfile_format, ') --------------------------------------- in atomic system'
  !write(*, "(A)")   "   Atomic system                  Fixed positions           Additional data"
  !do i = 1, maxLen
  !   write(at(i), "(a)") " "
  !   write(fixed(i), "(a)") " "
  !   write(add(i), "(a)") " "
  !end do
  !write(fixed(1), '(a)') "No fixed atom"
  !write(add(1), '(a)') "No symmetry for open BC"
  
  ! The atoms column
  !write(at(1), '(a,a)')  "Bound. C.= ", atoms%astruct%geocode
  !write(at(2), '(a,i5)') "N. types = ", atoms%astruct%ntypes
  !write(at(3), '(a,i5)') "N. atoms = ", atoms%astruct%nat
  !lg = 12
  !i = 4
  !write(at(i),'(a)' )    "Types    = "
  !do ityp=1,atoms%astruct%ntypes - 1
  !   if (lg + 4 + len(trim(atoms%astruct%atomnames(ityp))) >= width) then
  !      i = i + 1
  !      lg = 12
  !      write(at(i),'(a)') "           "
  !   end if
  !   write(at(i)(lg:),'(3a)') "'", trim(atoms%astruct%atomnames(ityp)), "', "
  !   lg = lg + 4 + len(trim(atoms%astruct%atomnames(ityp)))
  !end do
  !!Case no atom
  !if (atoms%astruct%ntypes > 0) then
  !   if (lg + 2 + len(trim(atoms%astruct%atomnames(ityp))) >= width) then
  !      i = i + 1
  !      lg = 12
  !      write(at(i),'(a)') "           "
  !   end if
  !   write(at(i)(lg:),'(3a)') "'", trim(atoms%astruct%atomnames(ityp)), "'"
  !end if

  ! The fixed atom column
  !i = 1
  !do iat=1,atoms%astruct%nat
  !   if (atoms%astruct%ifrztyp(iat)/=0) then
  !      if (i > maxLen) exit
  !      write(fixed(i),'(a,i4,a,a,a,i3)') &
  !           "at.", iat,' (', &
  !           & trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),&
  !           ') ',atoms%astruct%ifrztyp(iat)
  !      i = i + 1
  !   end if
  !end do
  !if (i > maxLen) write(fixed(maxLen), '(a)') " (...)"

  ! The additional data column
  !if (atoms%astruct%geocode /= 'F' .and. .not. in%disableSym) then
  !   call symmetry_get_matrices(atoms%astruct%sym%symObj, nSym, sym, transNon, symAfm, ierr)
  !   call symmetry_get_group(atoms%astruct%sym%symObj, spaceGroup, &
  !        & spaceGroupId, pointGroupMagn, genAfm, ierr)
  !   if (ierr == AB6_ERROR_SYM_NOT_PRIMITIVE) write(spaceGroup, "(A)") "not prim."
  !   write(add(1), '(a,i0)')       "N. sym.   = ", nSym
  !   write(add(2), '(a,a,a)')      "Sp. group = ", trim(spaceGroup)
  !else if (atoms%astruct%geocode /= 'F' .and. in%disableSym) then
  !   write(add(1), '(a)')          "N. sym.   = disabled"
  !   write(add(2), '(a)')          "Sp. group = disabled"
  !else
  !   write(add(1), '(a)')          "N. sym.   = free BC"
  !   write(add(2), '(a)')          "Sp. group = free BC"
  !end if
  !i = 3
  !if (in%nvirt > 0) then
  !   write(add(i), '(a,i5,a)')     "Virt. orb.= ", in%nvirt, " orb."
  !   write(add(i + 1), '(a,i5,a)') "Plot dens.= ", abs(in%nplot), " orb."
  !else
  !   write(add(i), '(a)')          "Virt. orb.= none"
  !   write(add(i + 1), '(a)')      "Plot dens.= none"
  !end if
  !i = i + 2
  !if (in%nspin==4) then
  !   write(add(i),'(a)')           "Spin pol. = non-coll."
  !else if (in%nspin==2) then
  !   write(add(i),'(a)')           "Spin pol. = collinear"
  !else if (in%nspin==1) then
  !   write(add(i),'(a)')           "Spin pol. = no"
  !end if

  ! Printing
  !do i = 1, maxLen
  !   if (len(trim(at(i))) > 0 .or. len(trim(fixed(i))) > 0 .or. len(trim(add(i))) > 0) then
  !      write(*,"(1x,a,1x,a,1x,a,1x,a,1x,a)") at(i), "|", fixed(i), "|", add(i)
  !   end if
  !end do

  !Geometry imput Parameters
  if (in%ncount_cluster_x > 0) then
     call yaml_comment('Geometry optimization Input Parameters (file: '//trim(prefix)//'.geopt)',hfill='-')
     call yaml_open_map('Geometry Optimization Parameters')
        call yaml_map('Maximum steps',in%ncount_cluster_x)
        call yaml_map('Algorithm', in%geopt_approach)
        call yaml_map('Random atomic displacement', in%randdis, fmt='(1pe7.1)')
        call yaml_map('Fluctuation in forces',in%frac_fluct,fmt='(1pe7.1)')
        call yaml_map('Maximum in forces',in%forcemax,fmt='(1pe7.1)')
        call yaml_map('Steepest descent step',in%betax,fmt='(1pe7.1)')
        if (trim(in%geopt_approach) == "DIIS") then
           call yaml_map('DIIS history', in%history)
        end if
     call yaml_close_map()
     if (case_insensitive_equiv(trim(in%geopt_approach),"AB6MD")) then
        call yaml_open_map('Molecular Dynamics Parameters')
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
        call yaml_close_map()
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
     if (in%disableSym .and. in%nkpt > 1) then
        call yaml_warning('symmetries have been disabled, k points are not irreductible.')
        !write(*, "(1x,A)") "WARNING: symmetries have been disabled, k points are not irreductible."
     end if
     call yaml_open_sequence('K points')!,advance='no')
     !call yaml_comment('Reduced coordinates  BZ coordinates  weight',hfill=' ')
     !write(*, "(1x,a)")    "       red. coordinates         weight       id        BZ coordinates"
     do i = 1, in%nkpt, 1
        call yaml_sequence(advance='no')
        call yaml_open_map(flow=.true.)
          call yaml_map( 'Rc', &
             & in%kpt(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), &
             & atoms%astruct%cell_dim(3) /) / two_pi,&
             & fmt='(f7.4)')
          call yaml_map( 'Bz', &
             & in%kpt(:, i), &
             & fmt='(f7.4)')
          call yaml_map('Wgt',in%wkpt(i),fmt='(f6.4)')
        call yaml_close_map(advance='no')
        call yaml_comment(trim(yaml_toa(i,fmt='(i4.4)')))
        !write(*, "(1x,3f9.5,2x,f9.5,5x,I4,1x,3f9.5)") &
        !     & in%kpt(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), atoms%astruct%cell_dim(3) /) / two_pi, &
        !     & in%wkpt(i), i, in%kpt(:, i)
     end do
     call yaml_close_sequence()

     if (in%nkptv > 0) then
        call yaml_open_sequence('K points for band structure calculation')
        !write(*, "(1x,a)")    " K points for band structure calculation"
        !write(*, "(1x,a)")    "       red. coordinates         weight       id        BZ coordinates"
        do i = 1, in%nkptv, 1
          call yaml_sequence(advance='no')
          call yaml_open_map(trim(yaml_toa(i,fmt='(i0)')),flow=.true.)
          call yaml_map( 'Red C.', &
             & in%kptv(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), &
             & atoms%astruct%cell_dim(3) /) / two_pi,&
             & fmt='(f9.5)')
          call yaml_map( 'Bz C.', &
             & in%kptv(:, i), &
             & fmt='(f9.5)')
          call yaml_map('Weight',1.0d0 / real(size(in%kptv, 2), gp),fmt='(f9.5)')
          call yaml_close_map()
        !   write(*, "(1x,3f9.5,2x,f9.5,5x,I4,1x,3f9.5)") &
        !        & in%kptv(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), atoms%astruct%cell_dim(3) /) / two_pi, &
        !        & 1.0d0 / real(size(in%kptv, 2), gp), i, in%kptv(:, i)
        end do
        call yaml_close_sequence()
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
     call yaml_open_map('Mixing parameters')
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
     call yaml_close_map()
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

  call yaml_open_map('DFT parameters')
    call yaml_open_map('eXchange Correlation')
      call yaml_map('XC ID',in%ixc,fmt='(i8)',label='ixc')
      call xc_dump()
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
    call yaml_close_map()

    if (in%ncharge > 0) call yaml_map('Net Charge (Ions-Electrons)',in%ncharge,fmt='(i8)')
    if (sqrt(sum(in%elecfield(:)**2)) > 0.0_gp) &
      call yaml_map('External Electric Field (Ha/a0)',in%elecfield(:),fmt='(1pe8.1)')
  call yaml_close_map()

  call yaml_open_map('Basis set definition')
    call yaml_map('Suggested Grid Spacings (a0)', (/in%hx,in%hy,in%hz/),fmt='(f5.2)')
    call yaml_map('Coarse and Fine Radii Multipliers', (/in%crmult,in%frmult/),fmt='(f4.1)')
  call yaml_close_map()


  call yaml_open_map('Self-Consistent Cycle Parameters')
    call yaml_open_map('Wavefunction')
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
  
    call yaml_close_map()

    call yaml_open_map('Density/Potential')
       call yaml_map('Max. Iterations',in%itrpmax)
    call yaml_close_map()
  call yaml_close_map()

  if (atoms%astruct%geocode == 'F') then
     call yaml_open_map('Post Optimization Parameters')

     call yaml_open_map('Finite-Size Effect estimation')
     call yaml_map('Scheduled',(in%rbuf > 0.0_gp))
     if (in%rbuf > 0.0_gp) then
        call yaml_map('Extension',in%rbuf,fmt='(f4.1)')
        call yaml_map('No. of CG steps',in%ncongt)
     end if
     call yaml_close_map()
     call yaml_close_map()
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


subroutine write_energies(iter,iscf,energs,gnrm,gnrm_zero,comment)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: iter,iscf
  type(energy_terms), intent(in) :: energs
  real(gp), intent(in) :: gnrm,gnrm_zero
  character(len=*), intent(in) :: comment
  !local variables

  if (len(trim(comment)) > 0) then
     if (verbose >0) call yaml_newline()
     call write_iter()
     if (verbose >0) call yaml_comment(trim(comment))    
  end if

  if (iscf < 1 .and. verbose > 0) then
     call yaml_newline()
     call yaml_open_map('Energies',flow=.true.)
  !call yaml_flow_map()
  !call yaml_indent_map('Energies')
     if (energs%ekin /= 0.0_gp)&
          call yaml_map('Ekin',energs%ekin,fmt='(1pe18.11)')
     if (energs%epot /= 0.0_gp)&
          call yaml_map('Epot',energs%epot,fmt='(1pe18.11)')
     if (energs%eproj /= 0.0_gp)&
          call yaml_map('Enl',energs%eproj,fmt='(1pe18.11)')
     if (energs%eh  /= 0.0_gp)&
          call yaml_map('EH',energs%eh,fmt='(1pe18.11)')
     if (energs%exc  /= 0.0_gp)&
     call yaml_map('EXC',energs%exc,fmt='(1pe18.11)')
     if (energs%evxc  /= 0.0_gp)&
          call yaml_map('EvXC',energs%evxc,fmt='(1pe18.11)')
     if (energs%eexctX  /= 0.0_gp)&
          call yaml_map('EexctX',energs%eexctX,fmt='(1pe18.11)')
     if (energs%evsic  /= 0.0_gp)&
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
     call yaml_close_map()
     call yaml_newline()
  end if

  if (len(trim(comment)) == 0) then
     call write_iter()
     if (verbose >0) call yaml_newline()
  else if (verbose > 1) then
     call yaml_map('SCF criterion',iscf,fmt='(i6)')
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
subroutine write_eigenvalues_data(nproc,etol,orbs,mom_vec)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: nproc
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

  if (verbose > 1) then
     call yaml_comment('Eigenvalues and New Occupation Numbers')
     !call yaml_open_map('Eigenvalues and New Occupation Numbers')
     !write(*,'(1x,a)')&
     !     &   '--------------------------------------- Kohn-Sham Eigenvalues and Occupation Numbers'
     ! Calculate and print the magnetisation
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
        !write(*,"(1x,A,f9.6)")"Total magnetisation: ", mpol
        call yaml_map("Total magnetization",mpol,fmt='(f9.6)')
     end if
     !if (orbs%nspinor ==4) then
     !   write(*,'(1x,a)')&
     !        &   '           Eigenvalue                                      m_x       m_y       m_z'
     !end if

     call yaml_open_sequence('Orbitals',flow=.true.)
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
                 call yaml_close_sequence(advance='no')
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
                 call yaml_close_sequence(advance='no')
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
                    call yaml_close_sequence(advance='no')
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
                 if (iorb == orbs%norbd .and. ikpt == orbs%nkpts) then
                    call yaml_close_sequence(advance='no')
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
     !call yaml_close_map()
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

  call yaml_open_sequence('Input Guess Orbitals',flow=.true.)!,advance='no')
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
         call yaml_close_sequence(advance='no')
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


  call yaml_open_map(flow=.true.)
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
  call yaml_close_map(advance='no')
 
END SUBROUTINE write_orbital_data


!> Write diis weights
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
           call yaml_open_sequence('DIIS weights',flow=.true.)
           do j=1,itdiis+1
              call yaml_sequence(yaml_toa(rds(1,j,1,1),fmt='(1pe9.2)'))
           end do
           call yaml_close_sequence()
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
  call yaml_open_sequence('Residues per orbital',flow=.true.)
  call yaml_newline()

  do ikpt=1,nkpts
     if (nkpts > 1) call yaml_comment('Kpt #'//adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))))
     do iorb=1,norb
        call yaml_sequence(trim(yaml_toa(gnrms(iorb,ikpt),fmt='(1pe19.12)')),advance='no')
        if (ikpt == nkpts .and. iorb == norb)   call yaml_close_sequence(advance='no')
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
   call yaml_open_sequence('Atomic Forces (Ha/Bohr)')
   do iat=1,atoms%astruct%nat
      call yaml_sequence(advance='no')
      call yaml_open_map(flow=.true.)
      call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),fxyz(1:3,iat),fmt='(1pe20.12)')
      !call yaml_map('AU',fxyz(1:3,iat),fmt='(1pe20.12)')
      !call yaml_map('eV/A',fxyz(1:3,iat)*Ha_eV/Bohr_Ang,fmt='(1pe9.2)')
      call yaml_close_map(advance='no')
      call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
!      write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
!      iat,trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),(fxyz(j,iat),j=1,3)
      sumx=sumx+fxyz(1,iat)
      sumy=sumy+fxyz(2,iat)
      sumz=sumz+fxyz(3,iat)
   end do
   call yaml_close_sequence()
   !$$        if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
   !$$           write(*,'(1x,a)')'the sum of the forces is'
   !$$           write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
   !$$           write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
   !$$           write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
   !$$        end if
END SUBROUTINE write_forces


!> Print the electronic configuration, with the semicore orbitals
subroutine print_eleconf(nspin,nspinor,noccmax,nelecmax,lmax,aocc,nsccode)
   use module_base
   use yaml_output
   implicit none
   integer, intent(in) :: nelecmax,nsccode,noccmax,lmax,nspinor,nspin
   real(gp), dimension(nelecmax), intent(in) :: aocc
   !local variables
   character(len=10) :: tmp
   character(len=500) :: string
   integer :: i,m,iocc,icoll,inl,noncoll,l,ispin,is,nl,niasc,lsc,nlsc,ntmp,iss
   logical, dimension(4,2) :: scorb

   !if non-collinear it is like nspin=1 but with the double of orbitals
   if (nspinor == 4) then
      noncoll=2
   else
      noncoll=1
   end if
   scorb=.false.
   if (nsccode/=0) then !the atom has some semicore orbitals
      niasc=nsccode
      do lsc=4,1,-1
         nlsc=niasc/4**(lsc-1)
         do i=1,nlsc
            scorb(lsc,i)=.true.
         end do
         niasc=niasc-nlsc*4**(lsc-1)
      end do
   end if

   call yaml_open_map('Electronic configuration',flow=.true.)

   !initalise string
   string=repeat(' ',len(string))

   is=1
   do i=1,noccmax
      iocc=0
      do l=1,lmax
         iocc=iocc+1
         nl=nint(aocc(iocc))
         do inl=1,nl
            !write to the string the angular momentum
            if (inl == i) then
               iss=is
               if (scorb(l,inl)) then
                  string(is:is)='('
                  is=is+1
               end if
               select case(l)
               case(1)
                  string(is:is)='s'
               case(2)
                  string(is:is)='p'
               case(3)
                  string(is:is)='d'
               case(4)
                  string(is:is)='f'
               case default
                  stop 'l not admitted'
               end select
               is=is+1
               if (scorb(l,inl)) then
                  string(is:is)=')'
                  is=is+1
               end if
               call yaml_open_sequence(string(iss:is))
            end if
            do ispin=1,nspin
               do m=1,2*l-1
                  do icoll=1,noncoll !non-trivial only for nspinor=4
                     iocc=iocc+1
                     !write to the string the value of the occupation numbers
                     if (inl == i) then
                        call write_fraction_string(l,aocc(iocc),tmp,ntmp)
                        string(is:is+ntmp-1)=tmp(1:ntmp)
                        call yaml_sequence(tmp(1:ntmp))
                        is=is+ntmp
                     end if
                  end do
               end do
            end do
            if (inl == i) then
               string(is:is+2)=' , '
               is=is+3
               call yaml_close_sequence()
            end if
         end do
      end do
   end do

   !write(*,'(2x,a,1x,a,1x,a)',advance='no')' Elec. Configuration:',trim(string),'...'

   call yaml_close_map()

END SUBROUTINE print_eleconf


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
  
  call yaml_open_sequence(trim(message)//' stress tensor matrix (Ha/Bohr^3)')
  call yaml_sequence(yaml_toa((/strten(1),strten(6),strten(5)/),fmt='(1pg20.12)'))
  call yaml_sequence(yaml_toa((/strten(6),strten(2),strten(4)/),fmt='(1pg20.12)'))
  call yaml_sequence(yaml_toa((/strten(5),strten(4),strten(3)/),fmt='(1pg20.12)'))
  call yaml_close_sequence()
  !write(*,'(1x,a)')'Stress Tensor, '//trim(message)//' contribution (Ha/Bohr^3):'
  !write(*,'(1x,t10,10x,a,t30,10x,a,t50,10x,a)')'x','y','z'
  !write(*,'(1x,a,t10,1pe20.12,t30,1pe20.12,t50,1pe20.12)')'x',strten(1),strten(6),strten(5)
  !write(*,'(1x,a,t30,1pe20.12,t50,1pe20.12)')'y',strten(2),strten(4)
  !write(*,'(1x,a,t50,1pe20.12)')'z',strten(3)

  if (fullinfo) then
     call yaml_open_map('Pressure')
     call yaml_map('Ha/Bohr^3',pressure,fmt='(1pg22.14)')
     call yaml_map('GPa',pressure*AU_GPa,fmt='(1pg14.6)')
     call yaml_map('PV (Ha)',pressure*volume,fmt='(1pg22.14)')
     call yaml_close_map()
     !write(*,'(1x,a,1pe22.14,a,1pe14.6,a,1pe22.14)')'Pressure:',pressure,&
     !     ' (',pressure*AU_GPa,' GPa), P V:',pressure*volume
  end if

END SUBROUTINE write_strten_info
