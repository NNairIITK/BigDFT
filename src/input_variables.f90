!!****f* BigDFT/print_logo
!! FUNCTION
!!    Display the logo of BigDFT 
!! COPYRIGHT
!!    Copyright (C) 2007-2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine print_logo()
  implicit none
  write(*,'(23x,a)')'      TTTT         F       DDDDD    '
  write(*,'(23x,a)')'     T    T               D         '
  write(*,'(23x,a)')'    T     T        F     D          '
  write(*,'(23x,a)')'    T    T         F     D        D '
  write(*,'(23x,a)')'    TTTTT          F     D         D'
  write(*,'(23x,a)')'    T    T         F     D         D'
  write(*,'(23x,a)')'    T     T        F     D         D'
  write(*,'(23x,a)')'    T      T       F     D         D'
  write(*,'(23x,a)')'    T     T     FFFF     D         D'
  write(*,'(23x,a)')'    T TTTT         F      D        D'
  write(*,'(23x,a)')'    T             F        D      D '
  write(*,'(23x,a)')'TTTTTTTTT    FFFFF          DDDDDD  ' 
  !write(*,'(23x,a)')'---------------------------------------'
  write(*,'(23x,a)')'  gggggg          iiiii    BBBBBBBBB'
  write(*,'(23x,a)')' g      g        i             B    '
  write(*,'(23x,a)')'g        g      i         BBBB B    '
  write(*,'(23x,a)')'g         g     iiii     B     B    '
  write(*,'(23x,a)')'g         g     i       B      B    '
  write(*,'(23x,a)')'g         g     i        B     B    '
  write(*,'(23x,a)')'g         g     i         B    B    '
  write(*,'(23x,a)')'g         g     i          BBBBB    '
  write(*,'(23x,a)')' g        g     i         B    B    '  
  write(*,'(23x,a)')'          g     i        B     B    ' 
  write(*,'(23x,a)')'         g               B    B     '
  write(*,'(23x,a)')'    ggggg       i         BBBB                     (Ver 1.3.0)'
  write(*,'(1x,a)')&
       '------------------------------------------------------------------------------------'
  write(*,'(1x,a)')&
       '|              Daubechies Wavelets for DFT Pseudopotential Calculations            |'
  write(*,'(1x,a)')&
       '------------------------------------------------------------------------------------'
  write(*,'(1x,a)')&
       '                                  The Journal of Chemical Physics 129, 014109 (2008)'
end subroutine print_logo
!!***

!!****f* BigDFT/dft_input_variables
!! FUNCTION
!!    Read the input variables needed for the DFT calculation
!!    The variables are divided in two groups:
!!    "cruising" variables -- general DFT run
!!    "brakeing" variables -- for the last run, once relaxation is achieved
!!                            of for a single-point calculation
!!    Every argument should be considered as mandatory
!! SOURCE
!!
subroutine dft_input_variables(iproc,filename,in)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  character(len=7) :: cudagpu
  character(len=100) :: line
  integer :: ierror,ierrfrc,iconv,iblas,iline,initerror

  ! Read the input variables.
  open(unit=1,file=filename,status='old')

  !line number, to control the input values
  iline=0
  !grid spacings
  read(1,*,iostat=ierror) in%hx,in%hy,in%hz
  call check()
  !coarse and fine radii around atoms
  read(1,*,iostat=ierror) in%crmult,in%frmult
  call check()
  !XC functional (ABINIT XC codes)
  read(1,*,iostat=ierror) in%ixc
  call check()
  !charged system, electric field (intensity and start-end points)
  read(1,'(a100)')line
  read(line,*,iostat=ierror) in%ncharge,in%ef(1)
  if (ierror == 0 .and. in%ef(1) /= 0.0_gp) then
     read(line,*,iostat=ierror) in%ncharge,in%ef(1),in%ef(2),in%ef(3)
  else
     in%ef(2)=0.0_gp
     in%ef(3)=0.0_gp
  end if
  call check()
  read(1,*,iostat=ierror) in%nspin,in%mpol
  call check()
  read(1,*,iostat=ierror) in%gnrm_cv
  call check()
  read(1,*,iostat=ierror) in%itermax,in%nrepmax
  call check()
  read(1,*,iostat=ierror) in%ncong,in%idsx
  call check()
  read(1,*,iostat=ierror) in%dispersion
  call check()
  !read the line for force the CUDA GPU calculation for all processors
  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) cudagpu
  iline=iline+1
  if (ierrfrc == 0 .and. cudagpu=='CUDAGPU') then
     call init_lib(iproc,initerror,iconv,iblas,GPUshare)
   !  iconv = 0
   !  iblas = 0
     if (initerror == 1) then
        stop
     end if
    ! GPUshare=.true.
     if (iconv == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUconv=.true.
     end if
     if (iblas == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUblas=.true.
     end if
  end if

  !now the varaibles which are to be used only for the last run
  read(1,*,iostat=ierror) in%inputPsiId,in%output_wf,in%output_grid
  call check()
  !project however the wavefunction on gaussians if asking to write them on disk
  in%gaussian_help=(in%inputPsiId >= 10)! commented .or. in%output_wf 
  !switch on the gaussian auxiliary treatment 
  !and the zero of the forces
  if (in%inputPsiId == 10) then
     in%inputPsiId=0
  end if
  read(1,*,iostat=ierror) in%rbuf,in%ncongt
  call check()
  in%calc_tail=(in%rbuf > 0.0_gp)

  !davidson treatment
  read(1,*,iostat=ierror) in%nvirt,in%nplot
  call check()

  !x-adsorber treatment (in progress)
  read(1,*,iostat=ierror)  in%iat_absorber
  call check()

  !electrostatic treatment of the vacancy (experimental)
  read(1,*,iostat=ierror)  in%nvacancy,in%read_ref_den,in%correct_offset,in%gnrm_sw
  call check()

  !verbosity of the output
  read(1,*,iostat=ierror)  in%verbosity
  call check()


  !performs some check: for the moment Davidson treatment is allowed only for spin-unpolarised
  !systems, while in principle it should work immediately
  if (in%nspin/=1 .and. in%nvirt/=0) then
     if (iproc==0) then
        write(*,'(1x,a)')'ERROR: Davidson treatment allowed only for non spin-polarised systems'
     end if
     stop
  end if
 
  close(unit=1,iostat=ierror)

  if (in%nvirt > 0 .and. iproc ==0) then
     !read virtual orbital and plotting request
     write(*,'(1x,a,i0)')'Virtual orbitals ',in%nvirt
     write(*,'(1x,a,i0,a)')'Output for density plots is requested for ',in%nplot,' orbitals'
  end if
  if (in%nspin==4) then
     if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation: YES (Non-collinear)'
  else if (in%nspin==2) then
     if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation: YES (Collinear)'
  else if (in%nspin==1) then
     if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation:  NO '
  else
     if (iproc == 0) write(*,'(1x,a,i0)')'Wrong spin polarisation id: ',in%nspin
     stop
  end if

contains

  subroutine check()
    iline=iline+1
    if (ierror/=0) then
       if (iproc == 0) write(*,'(1x,a,a,a,i3)') &
            'Error while reading the file "',trim(filename),'", line=',iline
       stop
    end if
  end subroutine check

end subroutine dft_input_variables
!!***

!!****f* BigDFT/geopt_input_variables_default
!! FUNCTION
!!    Assign default values for GEOPT variables
!! SOURCE
!!
subroutine geopt_input_variables_default(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(out) :: in

  !put some fake values for the geometry optimsation case
  in%geopt_approach='SDCG'
  in%ncount_cluster_x=0
  in%frac_fluct=1.0_gp
  in%forcemax=0.0_gp
  in%randdis=0.0_gp
  in%betax=2.0_gp

end subroutine geopt_input_variables_default


!!****f* BigDFT/geopt_input_variables
!! FUNCTION
!!    Read the input variables needed for the geometry optimisation
!!    Every argument should be considered as mandatory
!! SOURCE
!!
subroutine geopt_input_variables(iproc,filename,in)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  character(len=7) :: cudagpu
  character(len=100) :: line
  integer :: ierror,ierrfrc,iconv,iblas,iline,initerror

  ! Read the input variables.
  open(unit=1,file=filename,status='old')

  !line number, to control the input values
  iline=0

  read(1,*,iostat=ierror) in%geopt_approach
  call check()
  read(1,*,iostat=ierror) in%ncount_cluster_x
  call check()
  read(1,*,iostat=ierrfrc) in%frac_fluct,in%forcemax
  call check()
  read(1,*,iostat=ierror) in%randdis
  call check()
  read(1,*,iostat=ierror) in%betax
  call check()

 
  close(unit=1,iostat=ierror)

  if (iproc == 0) then
     write(*,'(1x,a,i0)') 'Max. number of wavefnctn optim ',in%ncount_cluster_x
     write(*,'(1x,a,1pe10.2)') 'Convergence criterion for forces: fraction of noise ',&
          in%frac_fluct
     write(*,'(1x,a,1pe10.2)') '                                : maximal component ',&
          in%forcemax
     write(*,'(1x,a,1pe10.2)') 'Random displacement amplitude ',in%randdis
     write(*,'(1x,a,1pe10.2)') 'Steepest descent step ',in%betax
  end if

contains

  subroutine check()
    iline=iline+1
    if (ierror/=0) then
       if (iproc == 0) write(*,'(1x,a,a,a,i3)') &
            'Error while reading the file "',trim(filename),'", line=',iline
       stop
    end if
  end subroutine check

end subroutine geopt_input_variables
!!***

!!****f* BigDFT/dft_input_converter
!! FUNCTION
!!  Convert the format of input variables
!! SOURCE
!!
subroutine dft_input_converter(in)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  !local variables
  character(len=7) :: cudagpu
  character(len=100) :: line
  integer :: ierror,ierrfrc,iconv,iblas,iline,initerror

  ! Read the input variables.
  open(unit=1,file='input_convert.dft',status='new')

  !line number, to control the input values
  iline=0
  !grid spacings
  line=''
  line=' hx,hy,hz: grid spacing in the three directions'
  write(1,'(3(f6.3),a)') in%hx,in%hy,in%hz,trim(line)
  !coarse and fine radii around atoms
  line=''
  line=' crmult, frmult: c(f)rmult*radii_cf(*,1(2)) gives the coarse (fine)radius around each atom'
  write(1,'(2(1x,f4.1),a)') in%crmult,in%frmult,trim(line)
  line=''
  line=' ixc: exchange-correlation parameter (LDA=1,PBE=11)'
  !XC functional (ABINIT XC codes)
  write(1,'(i3,a)') in%ixc,trim(line)

  line=''
  line=' ncharge: charge of the system, Electric field'
  write(1,'(i3,3(f6.3),a)') in%ncharge,in%ef(1),in%ef(2),in%ef(3),trim(line)

  line=''
  line=' nspin=1 non-spin polarization, mpol=total magnetic moment'
  write(1,'(2(i3),a)') in%nspin,in%mpol,trim(line)

  line=''
  line=' gnrm_cv: convergence criterion gradient'
  write(1,'(1pe7.0,a)') in%gnrm_cv,trim(line)
  
  line=''
  line=' itermax,nrepmax: maximum number of wavefunction optimizations and of re-diagonalised runs'
  write(1,'(2(i3),a)') in%itermax,in%nrepmax,trim(line)
  
  line=''
  line=' ncong, idsx: # CG iterations for the preconditioning equation, length of the diis history'
  write(1,'(2(i3),a)') in%ncong,in%idsx,trim(line)
  
  line=''
  line=' dispersion correction functional (values 1,2,3), 0=no correction'
  write(1,'(i3,a)') in%dispersion,trim(line)
  
  line=''
  line=' write "CUDAGPU" on this line to use GPU acceleration (GPU.config file is needed)'
  write(1,'(a)') trim(line)

  !now the varaibles which are to be used only for the last run
  line=''
  line=' InputPsiId, output_wf, output_grid'
  write(1,*) in%inputPsiId,in%output_wf,in%output_grid,trim(line)
  

  line=''
  line=' calc_tail, rbuf, ncongt: calculate tails,length of the tail (AU),# tail CG iterations'
  if (in%calc_tail) then
     write(1,'(f4.1,i4,a)') in%rbuf,in%ncongt,trim(line)
  else
     write(1,'(f4.1,i4,a)') 0.0_gp,in%ncongt,trim(line)
  end if


  !davidson treatment
  line=''
  line=' davidson treatment, no. of virtual orbitals, no of plotted orbitals'
  write(1,'(2(i3),a)') in%nvirt,in%nplot,trim(line)
  
  line=''
  line=' x-ray adsorber treatment'
  !x-adsorber treatment (in progress)
  write(1,'(i3,a)')  in%iat_absorber,trim(line)
  

  line=''
  line='0 .false. .false. 0.d0 vacancy: atom no., read_ref_den, correct_offset, gnrm_sw'
  !electrostatic treatment of the vacancy (experimental)
  write(1,*) trim(line)


  line=''
  line=' 2   verbosity of the output 0=low, 2=high'
  !electrostatic treatment of the vacancy (experimental)
  write(1,*) trim(line)
   
  close(unit=1)
end subroutine dft_input_converter
!!***





!!****f* BigDFT/read_input_variables
!! FUNCTION
!!    Read the input variables in the file 'input.dft'
!! SOURCE
!!
subroutine read_input_variables(iproc,filename,in)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  character(len=7) :: cudagpu
  character(len=100) :: line
  integer :: ierror,ierrfrc,iconv,iblas,iline,initerror

  ! Read the input variables.
  open(unit=1,file=filename,status='old')

  iline=0
  !read the line for force the CUDA GPU calculation for all processors
  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) cudagpu
  if (ierrfrc == 0 .and. cudagpu=='CUDAGPU') then
     call init_lib(iproc,initerror,iconv,iblas,GPUshare)
   !  iconv = 0
   !  iblas = 0
     

   !  call set_cpu_gpu_aff(iproc,iconv,iblas)
   ! GPUshare=.false.
  !   call init_gpu_sharing(initerror) !to fix the number of gpu and mpi tasks per node, we have to fil the inter_node.config file
     if (initerror == 1) then
        stop
     end if
    ! GPUshare=.true.
     if (iconv == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUconv=.true.
     end if
     if (iblas == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUblas=.true.
     end if
     read(1,*,iostat=ierror) in%ncount_cluster_x
     call check()
  else
     read(line,*,iostat=ierror) in%ncount_cluster_x
     call check()
  end if

  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) in%frac_fluct,in%forcemax
  if (ierrfrc /= 0) then
     read(line,*,iostat=ierror) in%frac_fluct
     in%forcemax=0.0_gp
  end if
  call check()
  read(1,*,iostat=ierror) in%randdis
  call check()
  read(1,*,iostat=ierror) in%betax
  call check()
  read(1,*,iostat=ierror) in%hx,in%hy,in%hz
  call check()
  read(1,*,iostat=ierror) in%crmult
  call check()
  read(1,*,iostat=ierror) in%frmult
  call check()

  read(1,*,iostat=ierror) in%ixc
  call check()
  read(1,'(a100)')line
  read(line,*,iostat=ierror) in%ncharge,in%ef(1)
  if (ierror == 0 .and. in%ef(1) /= 0.0_gp) then
     read(line,*,iostat=ierror) in%ncharge,in%ef(1),in%ef(2),in%ef(3)
  else
     in%ef(2)=0.0_gp
     in%ef(3)=0.0_gp
  end if
  call check()
  read(1,*,iostat=ierror) in%gnrm_cv
  call check()
  read(1,'(a100)')line
  read(line,*,iostat=ierror) in%itermax,in%nrepmax
  if (ierror == 0) then
     !read(line,*,iostat=ierror) in%ncharge,in%ef(1),in%ef(2),in%ef(3)
  else
     read(line,*,iostat=ierror)in%itermax
     in%nrepmax=10
  end if
  call check()
  read(1,*,iostat=ierror) in%ncong
  call check()
  read(1,*,iostat=ierror) in%idsx
  call check()
  read(1,*,iostat=ierror) in%calc_tail
  call check()
  read(1,*,iostat=ierror) in%rbuf
  call check()
  read(1,*,iostat=ierror) in%ncongt
  call check()
  read(1,*,iostat=ierror) in%nspin,in%mpol
  call check()
  read(1,*,iostat=ierror) in%inputPsiId,in%output_wf,in%output_grid
  call check()

  !project however the wavefunction on gaussians if asking to write them on disk
  in%gaussian_help=(in%inputPsiId >= 10)! commented .or. in%output_wf 
  !switch on the gaussian auxiliary treatment 
  !and the zero of the forces
  if (in%inputPsiId == 10) then
     in%inputPsiId=0
  end if

  ! qoh: Try to read dispersion input variable
  read(1,'(a100)',iostat=ierror)line
  if (ierror == 0) then
     if (index(line,"dispersion") /= 0) then 
        read(line,*,iostat=ierror) in%dispersion
        !add reading lines for Davidson treatment 
        !(optional for backward compatibility)
        read(1,*,iostat=ierror) in%nvirt, in%nplot
     else 
        in%dispersion = 0
        read(line,*,iostat=ierror) in%nvirt, in%nplot
     endif   
     ! add reading for absorbing atom. iat_absorber=0 ( default ) means no absorption calculation
     if(ierror==0) then
        read(1,*,iostat=ierror)  in%iat_absorber
        if(ierror/=0) then
           in%iat_absorber=0
        endif
     else
        in%iat_absorber=0
     endif

  ! AMmodif end

  else
     in%dispersion = 0
     in%nvirt=0
     in%nplot=0
     in%iat_absorber=0
 end if

  !performs some check: for the moment Davidson treatment is allowed only for spin-unpolarised
  !systems
  if (in%nspin/=1 .and. in%nvirt/=0) then
     if (iproc==0) then
        write(*,'(1x,a)')'ERROR: Davidson treatment allowed only for non spin-polarised systems'
     end if
     stop
  end if



 
  close(unit=1,iostat=ierror)

  if (iproc == 0) then
     write(*,'(1x,a,i0)') 'Max. number of wavefnctn optim ',in%ncount_cluster_x
     write(*,'(1x,a,1pe10.2)') 'Convergence criterion for forces: fraction of noise ',&
          in%frac_fluct
     write(*,'(1x,a,1pe10.2)') '                                : maximal component ',&
          in%forcemax
     write(*,'(1x,a,1pe10.2)') 'Random displacement amplitude ',in%randdis
     write(*,'(1x,a,1pe10.2)') 'Steepest descent step ',in%betax
     if (in%nvirt > 0) then
        !read virtual orbital and plotting request
        write(*,'(1x,a,i0)')'Virtual orbitals ',in%nvirt
        write(*,'(1x,a,i0,a)')'Output for density plots is requested for ',in%nplot,' orbitals'
     end if
  end if

     if (in%nspin==4) then
        if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation: YES (Non-collinear)'
     else if (in%nspin==2) then
        if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation: YES (Collinear)'
     else if (in%nspin==1) then
        if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation:  NO '
     else
        if (iproc == 0) write(*,'(1x,a,i0)')'Wrong spin polarisation id: ',in%nspin
        stop
     end if

contains

  subroutine check()
    iline=iline+1
    if (ierror/=0) then
       if (iproc == 0) write(*,'(1x,a,a,a,i3)') &
            'Error while reading the file "',trim(filename),'", line=',iline
       stop
    end if
  end subroutine check

end subroutine read_input_variables
!!***


!!****f* BigDFT/print_input_parameters
!! FUNCTION
!!    Print all input parameters
!! SOURCE
!!
subroutine print_input_parameters(in,atoms)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms

  write(*,'(1x,a)')&
       '------------------------------------------------------------------- Input Parameters'
  write(*,'(1x,a)')&
       '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
  write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
       'Max. hgrid  =',in%hx,   '|  Coarse Wfs.=',in%crmult,'| Wavefns Conv.=',in%gnrm_cv,&
       '| Calculate=',in%calc_tail
  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i5,a,i2,1x,a,f4.1)')&
       '       XC id=',in%ixc,     '|    Fine Wfs.=',in%frmult,'| Max. N. Iter.=',in%itermax,&
       'x',in%nrepmax,'| Extension=',in%rbuf
  write(*,'(1x,a,i7,1x,a,1x,a,i8,1x,a,i4)')&
       'total charge=',in%ncharge, '|                   ','| CG Prec.Steps=',in%ncong,&
       '|  CG Steps=',in%ncongt
  write(*,'(1x,a,1pe7.1,1x,a,1x,a,i8)')&
       ' elec. field=',in%ef(1),'|                   ','| DIIS Hist. N.=',in%idsx
  if (in%nspin>=2) then
     write(*,'(1x,a,i7,1x,a)')&
          'Polarisation=',2*in%mpol, '|'
  end if
  if (atoms%geocode /= 'F') then
     write(*,'(1x,a,1x,a,3(1x,1pe12.5))')&
          '  Geom. Code=    '//atoms%geocode//'   |',&
          '  Box Sizes (Bohr) =',atoms%alat1,atoms%alat2,atoms%alat3

  end if
end subroutine print_input_parameters
!!***

!!****f* BigDFT/read_atomic_file
!! FUNCTION
!!    Read atomic file
!! SOURCE
!!
subroutine read_atomic_file(file,iproc,at,rxyz)
  use module_base
  use module_types
  use module_interfaces, except_this_one => read_atomic_file
  implicit none
  character(len=*), intent(in) :: file
  integer, intent(in) :: iproc
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(:,:), pointer :: rxyz
  !local variables
  character(len=*), parameter :: subname='read_atomic_file'
  integer :: i_stat, l
  logical :: file_exists
  character(len = 128) :: filename

  file_exists = .false.

  ! Test posinp.xyz
  if (.not. file_exists) then
     inquire(FILE = file//'.xyz', EXIST = file_exists)
     if (file_exists) write(filename, "(A)") file//'.xyz'!"posinp.xyz"
     write(at%format, "(A)") "xyz"
  end if
  ! Test posinp.ascii
  if (.not. file_exists) then
     inquire(FILE = file//'.ascii', EXIST = file_exists)
     if (file_exists) write(filename, "(A)") file//'.ascii'!"posinp.ascii"
     write(at%format, "(A)") "ascii"
  end if
  ! Test the name directly
  if (.not. file_exists) then
     inquire(FILE = file, EXIST = file_exists)
     if (file_exists) write(filename, "(A)") file
     l = len(file)
     if (file(l-3:l) == ".xyz") then
        write(at%format, "(A)") "xyz"
     else if (file(l-5:l) == ".ascii") then
        write(at%format, "(A)") "ascii"
     else
        write(*,*) "Atomic input file format not recognised."
        write(*,*) " File should be *.ascii or *.xyz."
        stop
     end if
  end if

  if (.not. file_exists) then
     write(*,*) "Atomic input file not found."
     write(*,*) " Files looked for were '"//file//"'.ascii, '"//file//".xyz' and '"//file//"'."
     stop 
  end if

  open(unit=99,file=trim(filename),status='old')

  if (at%format == "xyz") then
     read(99,*) at%nat,at%units
 
     allocate(rxyz(3,at%nat+ndebug),stat=i_stat)
     call memocc(i_stat,rxyz,'rxyz',subname)

     !read atomic positions
     call read_atomic_positions(iproc,99,at,rxyz)
  else if (at%format == "ascii") then
     !read atomic positions
     call read_atomic_ascii(iproc,99,at,rxyz)
  end if

  close(99)
end subroutine read_atomic_file
!!***

!!****f* BigDFT/read_atomic_positions
!! FUNCTION
!!    Read atomic positions
!! SOURCE
!!
subroutine read_atomic_positions(iproc,ifile,at,rxyz)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(3,at%nat), intent(out) :: rxyz
  !local variables
  character(len=*), parameter :: subname='read_atomic_positions'
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem
  character(len=2) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl,dowrite
  integer :: nateq,iat,jat,ityp,i,ierror,ierrsfx,i_stat,j
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3
! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0
  character(len=20), dimension(100) :: atomnames

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atoms     = ',at%nat

  allocate(at%iatype(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%iatype,'at%iatype',subname)
  allocate(at%ifrztyp(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%ifrztyp,'at%ifrztyp',subname)
  allocate(at%natpol(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%natpol,'at%natpol',subname)

  !controls if the positions are provided with machine precision
  if (at%units == 'angstroemd0' .or. at%units== 'atomicd0' .or. at%units== 'bohrd0') then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  !this array is useful for frozen atoms
  !no atom is frozen by default
  at%ifrztyp(:)=0
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  at%natpol(:)=100

  !read from positions of .xyz format, but accepts also the old .ascii format
  read(ifile,'(a150)')line

!!$  !old format, still here for backward compatibility
!!$  !admits only simple precision calculation
!!$  read(line,*,iostat=ierror) rx,ry,rz,tatonam

!!$  !in case of old format, put geocode to F and alat to 0.
!!$  if (ierror == 0) then
!!$     at%geocode='F'
!!$     alat1d0=0.0_gp
!!$     alat2d0=0.0_gp
!!$     alat3d0=0.0_gp
!!$  else
  if (lpsdbl) then
     read(line,*,iostat=ierrsfx) tatonam,alat1d0,alat2d0,alat3d0
  else
     read(line,*,iostat=ierrsfx) tatonam,alat1,alat2,alat3
  end if
  if (ierrsfx == 0) then
     if (trim(tatonam)=='periodic') then
        at%geocode='P'
     else if (trim(tatonam)=='surface') then 
        at%geocode='S'
        at%alat2=0.0_gp
     else !otherwise free bc
        at%geocode='F'
        at%alat1=0.0_gp
        at%alat2=0.0_gp
        at%alat3=0.0_gp
     end if
     if (.not. lpsdbl) then
        alat1d0=real(alat1,gp)
        alat2d0=real(alat2,gp)
        alat3d0=real(alat3,gp)
     end if
  else
     at%geocode='F'
     alat1d0=0.0_gp
     alat2d0=0.0_gp
     alat3d0=0.0_gp
  end if
!!$  end if

  !reduced coordinates are possible only with periodic units
  if (at%units == 'reduced' .and. at%geocode == 'F') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are not allowed with isolated BC'
  end if

  !convert the values of the cell sizes in bohr
  if (at%units=='angstroem' .or. at%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     at%alat1=alat1d0/bohr
     at%alat2=alat2d0/bohr
     at%alat3=alat3d0/bohr
  else if  (at%units=='atomic' .or. at%units=='bohr'  .or.&
       at%units== 'atomicd0' .or. at%units== 'bohrd0') then
     at%alat1=alat1d0
     at%alat2=alat2d0
     at%alat3=alat3d0
  else if (at%units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
     at%alat1=real(alat1,gp)
     at%alat2=real(alat2,gp)
     at%alat3=real(alat3,gp)
  else
     write(*,*) 'length units in input file unrecognized'
     write(*,*) 'recognized units are angstroem or atomic = bohr'
     stop 
  endif

  at%ntypes=0
  do iat=1,at%nat
!!$     if (ierror == 0) then
!!$        !old case of ascii file, added for backward compatibility
!!$        if (iat /= 1) read(ifile,*) rx,ry,rz,tatonam
!!$     else
     !xyz input file, allow extra information
     read(ifile,'(a150)')line 
     if (lpsdbl) then
        read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,extra
     else
        read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,extra
     end if
     !print *,line
     call find_extra_info(line,extra)
     call parse_extra_info(iproc,iat,extra,at)

     tatonam=trim(symbol)
!!$     end if
     if (lpsdbl) then
        rxyz(1,iat)=rxd0
        rxyz(2,iat)=ryd0
        rxyz(3,iat)=rzd0
     else
        rxyz(1,iat)=real(rx,gp)
        rxyz(2,iat)=real(ry,gp)
        rxyz(3,iat)=real(rz,gp)
     end if
     if (at%units == 'reduced') then !add treatment for reduced coordinates
        rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
        if (at%geocode == 'P') rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
        rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
     else if (at%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
        rxyz(2,iat)=modulo(rxyz(2,iat),alat2d0)
        rxyz(3,iat)=modulo(rxyz(3,iat),alat3d0)
     else if (at%geocode == 'S') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
        rxyz(3,iat)=modulo(rxyz(3,iat),alat3d0)
     end if
 
     do ityp=1,at%ntypes
        if (tatonam == atomnames(ityp)) then
           at%iatype(iat)=ityp
           goto 200
        endif
     enddo
     at%ntypes=at%ntypes+1
     if (at%ntypes > 100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     at%iatype(iat)=at%ntypes
200  continue

     if (at%units=='angstroem' .or. at%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/bohr
        enddo
     else if (at%units == 'reduced') then 
        rxyz(1,iat)=rxyz(1,iat)*at%alat1
        if (at%geocode == 'P') rxyz(2,iat)=rxyz(2,iat)*at%alat2
        rxyz(3,iat)=rxyz(3,iat)*at%alat3
     endif
  enddo

  !now that ntypes is determined allocate at%atomnames and copy the values
  allocate(at%atomnames(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%atomnames,'at%atomnames',subname)
  at%atomnames(1:at%ntypes)=atomnames(1:at%ntypes)

  !control atom positions
  call check_atoms_positions(iproc,at,rxyz)

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atom types= ',at%ntypes

  do ityp=1,at%ntypes
     if (iproc.eq.0) &
          write(*,'(1x,a,i0,a,a)') 'Atoms of type ',ityp,' are ',trim(at%atomnames(ityp))
  enddo

  do iat=1,at%nat
     if (iproc.eq.0 .and. at%ifrztyp(iat)/=0) &
          write(*,'(1x,a,i0,a,a,a,i3)') &
          'FIXED Atom N.:',iat,', Name: ',trim(at%atomnames(at%iatype(iat))),&
          ', ifrztyp= ',at%ifrztyp(iat)
  enddo

end subroutine read_atomic_positions
!!***


subroutine check_atoms_positions(iproc,at,rxyz)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables
  logical :: dowrite
  integer :: iat,nateq,jat,j

  nateq=0
  do iat=1,at%nat
     do jat=iat+1,at%nat
        if ((rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+&
             (rxyz(3,iat)-rxyz(3,jat))**2 ==0.0_gp) then
           nateq=nateq+1
           write(*,'(1x,a,2(i0,a,a6,a))')'ERROR: atoms ',iat,&
                ' (',trim(at%atomnames(at%iatype(iat))),') and ',&
                jat,' (',trim(at%atomnames(at%iatype(jat))),&
                ') have the same positions'
        end if
     end do
  end do
  if (nateq /= 0) then
     if (iproc == 0) then
        write(*,'(1x,a)')'Control your posinp file, cannot proceed'
        write(*,'(1x,a)',advance='no')&
             'Writing tentative alternative positions in the file posinp_alt...'
        open(unit=9,file='posinp_alt')
        write(9,'(1x,a)')' ??? atomicd0'
        write(9,*)
        do iat=1,at%nat
           dowrite=.true.
           do jat=iat+1,at%nat
              if ((rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+&
                   (rxyz(3,iat)-rxyz(3,jat))**2 ==0.0_gp) then
                 dowrite=.false.
              end if
           end do
           if (dowrite) & 
                write(9,'(a2,4x,3(1x,1pe21.14))')trim(at%atomnames(at%iatype(iat))),&
                (rxyz(j,iat),j=1,3)
        end do
        close(9)
        write(*,'(1x,a)')' done.'
        write(*,'(1x,a)')' Replace ??? in the file heading with the actual atoms number'               
     end if
     stop
  end if
end subroutine check_atoms_positions


!!****f* BigDFT/find_extra_info
!! FUNCTION
!!    Find extra information
!!
!! SOURCE
!!
subroutine find_extra_info(line,extra)
  implicit none
  character(len=150), intent(in) :: line
  character(len=50), intent(out) :: extra
  !local variables
  logical :: space
  integer :: i,nspace
  i=1
  space=.true.
  nspace=-1
  find_space : do
     !toggle the space value for each time
     if (line(i:i) == ' ' .neqv. space) then
        nspace=nspace+1
        space=.not. space
     end if
     !print *,line(i:i),nspace
     if (nspace==8) then
        extra=line(i:min(150,i+49))
        exit find_space
     end if
     if (i==150) then
        extra='nothing'
        exit find_space
     end if
     i=i+1
  end do find_space

end subroutine find_extra_info
!!***


!!****f* BigDFT/parse_extra_info
!! FUNCTION
!!    Parse extra information
!!
!! SOURCE
!!
subroutine parse_extra_info(iproc,iat,extra,at)
  use module_types
  implicit none
  character(len=50), intent(in) :: extra
  integer, intent(in) :: iat,iproc
  type(atoms_data), intent(out) :: at
  !local variables
  character(len=4) :: suffix
  logical :: go
  integer :: ierr,ierr1,ierr2,nspol,nchrg,nsgn
  !case with all the information
  !print *,iat,'ex'//trim(extra)//'ex'
  read(extra,*,iostat=ierr)nspol,nchrg,suffix
  if (extra == 'nothing') then !case with empty information
     nspol=0
     nchrg=0
     suffix='    '
  else if (ierr /= 0) then !case with partial information
     read(extra,*,iostat=ierr1)nspol,suffix
     if (ierr1 /=0) then
        call valid_frzchain(trim(extra),go)
        if (go) then
           suffix=trim(extra)
           nspol=0
           nchrg=0
        else
           read(extra,*,iostat=ierr2)nspol
           if (ierr2 /=0) then
              call error
           end if
           suffix='    '
           nchrg=0
        end if
     else
        nchrg=0
        call valid_frzchain(trim(extra),go)
        if (.not. go) then
           read(suffix,*,iostat=ierr2)nchrg
           if (ierr2 /= 0) then
              call error
           else
              suffix='    '
           end if
        end if
     end if
  end if

  !now assign the array, following the rule
  if(nchrg>=0) then
     nsgn=1
  else
     nsgn=-1
  end if
  at%natpol(iat)=1000*nchrg+nsgn*100+nspol

  !print *,'natpol atomic',iat,at%natpol(iat)

  !convert the suffix into ifrztyp
  call frozen_ftoi(suffix,at%ifrztyp(iat))

!!$  if (trim(suffix) == 'f') then
!!$     !the atom is considered as blocked
!!$     at%ifrztyp(iat)=1
!!$  end if

contains

 subroutine error
   if (iproc == 0) then
      print *,extra
      write(*,'(1x,a,i0,a)')&
           'ERROR in input file for atom number ',iat,&
           ': after 4th column you can put the input polarisation(s) or the frzchain: f,fxz,fy'
   end if
   stop
 end subroutine error
  
end subroutine parse_extra_info
!!***

!!****f* BigDFT/read_atomic_ascii
!! FUNCTION
!!    Read atomic positions of ascii files.
!! SOURCE
!!
subroutine read_atomic_ascii(iproc,ifile,at,rxyz)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(:,:), pointer :: rxyz
  !local variables
  character(len=*), parameter :: subname='read_atomic_ascii'
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem
  character(len=2) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl,dowrite
  integer :: nateq,iat,jat,ityp,i,i_stat,j,nlines
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3,alat4,alat5,alat6
! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0,alat4d0,alat5d0,alat6d0
  character(len=20), dimension(100) :: atomnames
  ! Store the file.
  character(len = 150), dimension(5000) :: lines

  ! First pass to store the file in a string buffer.
  nlines = 1
  do
     read(ifile,'(a150)', iostat = i_stat) lines(nlines)
     if (i_stat /= 0) then
        exit
     end if
     nlines = nlines + 1
     if (nlines > 5000) then
        if (iproc==0) write(*,*) 'Atomic input file too long (> 5000 lines).'
        stop 
     end if
  end do
  nlines = nlines - 1

  if (nlines < 4) then
     if (iproc==0) write(*,*) 'Error in ASCII file format, file has less than 4 lines.'
     stop 
  end if

  ! Try to determine the number atoms and the keywords.
  write(at%units, "(A)") "bohr"
  at%geocode = 'P'
  at%nat     = 0
  do i = 4, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (line(1:1) /= '#' .and. line(1:1) /= '!' .and. len(trim(line)) /= 0) then
        at%nat = at%nat + 1
     else if (line(1:9) == "#keyword:" .or. line(1:9) == "!keyword:") then
        if (index(line, 'bohr') > 0)        write(at%units, "(A)") "bohr"
        if (index(line, 'bohrd0') > 0)      write(at%units, "(A)") "bohrd0"
        if (index(line, 'atomic') > 0)      write(at%units, "(A)") "atomicd0"
        if (index(line, 'angstroem') > 0)   write(at%units, "(A)") "angstroem"
        if (index(line, 'angstroemd0') > 0) write(at%units, "(A)") "angstroemd0"
        if (index(line, 'reduced') > 0)     write(at%units, "(A)") "reduced"
        if (index(line, 'periodic') > 0) at%geocode = 'P'
        if (index(line, 'surface') > 0)  at%geocode = 'S'
        if (index(line, 'freeBC') > 0)   at%geocode = 'F'
     end if
  end do
  
  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atoms     = ',at%nat

  allocate(at%iatype(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%iatype,'at%iatype',subname)
  allocate(at%ifrztyp(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%ifrztyp,'at%ifrztyp',subname)
  allocate(at%natpol(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%natpol,'at%natpol',subname)
  allocate(rxyz(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%natpol,'rxyz',subname)

  !controls if the positions are provided with machine precision
  if (index(at%units, 'd0') > 0) then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  !this array is useful for frozen atoms
  !no atom is frozen by default
  at%ifrztyp(:)=0
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  at%natpol(:)=100

  ! Read the box definition
  at%alat1 = 0.0_gp
  at%alat2 = 0.0_gp
  at%alat3 = 0.0_gp
  if (lpsdbl) then
     read(lines(2),*) alat1d0,alat2d0,alat3d0
     read(lines(3),*) alat4d0,alat5d0,alat6d0
     if (alat2d0 /= 0.d0 .or. alat4d0 /= 0.d0 .or. alat5d0 /= 0.d0) then
        if (iproc==0) write(*,*) 'Only orthorombic boxes are possible.'
        stop 
     end if
     at%alat1 = real(alat1d0,gp)
     at%alat2 = real(alat3d0,gp)
     at%alat3 = real(alat6d0,gp)
  else
     read(lines(2),*) alat1,alat2,alat3
     read(lines(3),*) alat4,alat5,alat6
     if (alat2 /= 0. .or. alat4 /= 0. .or. alat5 /= 0.) then
        if (iproc==0) write(*,*) 'Only orthorombic boxes are possible.'
        if (iproc==0) write(*,*) ' but alat2, alat4 and alat5 = ', alat2, alat4, alat5
        stop 
     end if
     at%alat1 = real(alat1,gp)
     at%alat2 = real(alat3,gp)
     at%alat3 = real(alat6,gp)
  end if
  if (at%geocode == 'S') then
     at%alat2 = 0.0_gp
  else if (at%geocode == 'F') then
     at%alat1 = 0.0_gp
     at%alat2 = 0.0_gp
     at%alat3 = 0.0_gp
  end if
  
  !reduced coordinates are possible only with periodic units
  if (at%units == 'reduced' .and. at%geocode /= 'P') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are only allowed with fully periodic BC'
  end if

  !convert the values of the cell sizes in bohr
  if (at%units=='angstroem' .or. at%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     at%alat1 = at%alat1 / bohr
     at%alat2 = at%alat2 / bohr
     at%alat3 = at%alat3 / bohr
  endif

  at%ntypes=0
  iat = 1
  do i = 4, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (line(1:1) /= '#' .and. line(1:1) /= '!' .and. len(trim(line)) /= 0) then
        write(extra, "(A)") "nothing"
        if (lpsdbl) then
           read(line,*, iostat = i_stat) rxd0,ryd0,rzd0,symbol,extra
           if (i_stat /= 0) read(line,*) rxd0,ryd0,rzd0,symbol
        else
           read(line,*, iostat = i_stat) rx,ry,rz,symbol,extra
           if (i_stat /= 0) read(line,*) rx,ry,rz,symbol
        end if
        call find_extra_info(line,extra)
        call parse_extra_info(iproc,iat,extra,at)

        tatonam=trim(symbol)

        if (lpsdbl) then
           rxyz(1,iat)=rxd0
           rxyz(2,iat)=ryd0
           rxyz(3,iat)=rzd0
        else
           rxyz(1,iat)=real(rx,gp)
           rxyz(2,iat)=real(ry,gp)
           rxyz(3,iat)=real(rz,gp)
        end if

        if (at%units == 'reduced') then !add treatment for reduced coordinates
           rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
           rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
           rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
        else if (at%geocode == 'P') then
           rxyz(1,iat)=modulo(rxyz(1,iat),at%alat1)
           rxyz(2,iat)=modulo(rxyz(2,iat),at%alat2)
           rxyz(3,iat)=modulo(rxyz(3,iat),at%alat3)
        else if (at%geocode == 'S') then
           rxyz(1,iat)=modulo(rxyz(1,iat),at%alat1)
           rxyz(3,iat)=modulo(rxyz(3,iat),at%alat3)
        end if
 
        do ityp=1,at%ntypes
           if (tatonam == atomnames(ityp)) then
              at%iatype(iat)=ityp
              goto 200
           endif
        enddo
        at%ntypes=at%ntypes+1
        if (at%ntypes > 100) stop 'more than 100 atomnames not permitted'
        atomnames(ityp)=tatonam
        at%iatype(iat)=at%ntypes
200     continue

        if (at%units=='angstroem' .or. at%units=='angstroemd0') then
           ! if Angstroem convert to Bohr
           do j=1,3 
              rxyz(j,iat)=rxyz(j,iat)/bohr
           enddo
        else if (at%units == 'reduced') then 
           rxyz(1,iat)=rxyz(1,iat)*at%alat1
           rxyz(2,iat)=rxyz(2,iat)*at%alat2
           rxyz(3,iat)=rxyz(3,iat)*at%alat3
        endif
        iat = iat + 1
     end if
  enddo

  !now that ntypes is determined allocate at%atomnames and copy the values
  allocate(at%atomnames(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%atomnames,'at%atomnames',subname)
  at%atomnames(1:at%ntypes)=atomnames(1:at%ntypes)

  !control atom positions
  call check_atoms_positions(iproc,at,rxyz)

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atom types= ',at%ntypes

  do ityp=1,at%ntypes
     if (iproc.eq.0) &
          write(*,'(1x,a,i0,a,a)') 'Atoms of type ',ityp,' are ',trim(at%atomnames(ityp))
  enddo

  do iat=1,at%nat
     if (iproc.eq.0 .and. at%ifrztyp(iat) /=0) &
          write(*,'(1x,a,i0,a,a)') 'FIXED Atom N.:',iat,', Name: ',trim(at%atomnames(at%iatype(iat)))
  enddo

end subroutine read_atomic_ascii
!!***

!!****f* BigDFT/charge_and_spol
!! FUNCTION
!!   Calculate the charge and the spin polarisation to be placed on a given atom
!!   RULE: natpol = c*1000 + sgn(c)*100 + s: charged and polarised atom (charge c, polarisation s)
!! SOURCE
!!
subroutine charge_and_spol(natpol,nchrg,nspol)
  implicit none
  integer, intent(in) :: natpol
  integer, intent(out) :: nchrg,nspol
  !local variables
  integer :: nsgn

  nchrg=natpol/1000
  if (nchrg>=0) then
     nsgn=1
  else
     nsgn=-1
  end if

  nspol=natpol-1000*nchrg-nsgn*100

end subroutine charge_and_spol
!!***

subroutine write_atomic_file(filename,energy,rxyz,atoms,comment)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename,comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz

  if (atoms%format == "xyz") then
     call wtxyz(filename,energy,rxyz,atoms,comment)
  else if (atoms%format == "ascii") then
     call wtascii(filename,energy,rxyz,atoms,comment)
  else
     write(*,*) "Error, unknown file format."
     stop
  end if
end subroutine write_atomic_file

subroutine wtxyz(filename,energy,rxyz,atoms,comment)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename,comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !local variables
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem
  character(len=2) :: symbol
  character(len=10) :: name
  character(len=11) :: units
  character(len=50) :: extra
  integer :: iat,j,ichg,ispol
  real(gp) :: xmax,ymax,zmax,factor

  open(unit=9,file=filename//'.xyz')
  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,atoms%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(atoms%units) == 'angstroem' .or. trim(atoms%units) == 'angstroemd0') then
     factor=bohr
     units='angstroemd0'
  else
     factor=1.0_gp
     units='atomicd0'
  end if


  write(9,'(i6,2x,a,2x,1pe24.17,2x,a)') atoms%nat,trim(units),energy,comment

  if (atoms%geocode == 'P') then
     write(9,'(a,3(1x,1pe24.17))')'periodic',&
          atoms%alat1*factor,atoms%alat2*factor,atoms%alat3*factor
  else if (atoms%geocode == 'S') then
     write(9,'(a,3(1x,1pe24.17))')'surface',&
          atoms%alat1*factor,atoms%alat2*factor,atoms%alat3*factor
  else
     write(9,*)'free'
  end if
  do iat=1,atoms%nat
     name=trim(atoms%atomnames(atoms%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:2)
     end if

     call write_extra_info(extra,atoms%natpol(iat),atoms%ifrztyp(iat))

     write(9,'(a2,4x,3(1x,1pe24.17),2x,a50)')symbol,(rxyz(j,iat)*factor,j=1,3),extra

  enddo
  close(unit=9)

end subroutine wtxyz

subroutine wtascii(filename,energy,rxyz,atoms,comment)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename,comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !local variables
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem
  character(len=2) :: symbol
  character(len=50) :: extra
  character(len=10) :: name
  integer :: iat,j,ichg,ispol
  real(gp) :: xmax,ymax,zmax,factor

  open(unit=9,file=filename//'.ascii')
  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,atoms%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(atoms%units) == 'angstroem' .or. trim(atoms%units) == 'angstroemd0') then
     factor=bohr
  else
     factor=1.0_gp
  end if

  write(9, "(A,A)") "# BigDFT file - ", trim(comment)
  write(9, "(3e24.17)") atoms%alat1, 0.d0, atoms%alat2
  write(9, "(3e24.17)") 0.d0,        0.d0, atoms%alat3

  write(9, "(A,A)") "#keyword: ", trim(atoms%units)
  if (atoms%geocode == 'P') write(9, "(A)") "#keyword: periodic"
  if (atoms%geocode == 'S') write(9, "(A)") "#keyword: surface"
  if (atoms%geocode == 'F') write(9, "(A)") "#keyword: freeBC"
  if (energy /= 0.d0) then
     write(9, "(A,e24.17)") "# Total energy (Ht): ", energy
  end if

  do iat=1,atoms%nat
     name=trim(atoms%atomnames(atoms%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:2)
     end if

     call write_extra_info(extra,atoms%natpol(iat),atoms%ifrztyp(iat))     

     write(9,'(3(1x,1pe24.17),2x,a2,2x,a50)') (rxyz(j,iat)*factor,j=1,3),symbol,extra
  end do
  close(unit=9)

end subroutine wtascii

!write the extra info necessary for the output file
subroutine write_extra_info(extra,natpol,ifrztyp)
  use module_base
  implicit none 
  integer, intent(in) :: natpol,ifrztyp
  character(len=50), intent(out) :: extra
  !local variables
  character(len=4) :: frzchain
  integer :: ispol,ichg

  call charge_and_spol(natpol,ichg,ispol)

  call frozen_itof(ifrztyp,frzchain)
  
  !takes into account the blocked atoms and the input polarisation
  if (ispol == 0 .and. ichg == 0 ) then
     write(extra,'(2x,a4)')frzchain
  else if (ispol /= 0 .and. ichg == 0) then
     write(extra,'(i7,2x,a4)')ispol,frzchain
  else if (ichg /= 0) then
     write(extra,'(2(i7),2x,a4)')ispol,ichg,frzchain
  else
     write(extra,'(2x,a4)') ''
  end if
  
end subroutine write_extra_info

subroutine frozen_itof(ifrztyp,frzchain)
  implicit none
  integer, intent(in) :: ifrztyp
  character(len=4), intent(out) :: frzchain

  if (ifrztyp == 0) then
     frzchain='    '
  else if (ifrztyp == 1) then
     frzchain='   f'
  else if (ifrztyp == 2) then
     frzchain='  fy'
  else if (ifrztyp == 3) then
     frzchain=' fxz'
  end if
        
end subroutine frozen_itof

subroutine valid_frzchain(frzchain,go)
  character(len=*), intent(in) :: frzchain
  logical, intent(out) :: go

  go= trim(frzchain) == 'f' .or. &
       trim(frzchain) == 'fy' .or. &
       trim(frzchain) == 'fxz'
  
end subroutine valid_frzchain

subroutine frozen_ftoi(frzchain,ifrztyp)
  implicit none
  character(len=4), intent(in) :: frzchain
  integer, intent(out) :: ifrztyp

  if (trim(frzchain)=='') then
     ifrztyp = 0
  else if (trim(frzchain)=='f') then
     ifrztyp = 1
  else if (trim(frzchain)=='fy') then
     ifrztyp = 2
  else if (trim(frzchain)=='fxz') then
     ifrztyp = 3
  end if
        
end subroutine frozen_ftoi

!calculate the coefficient for moving atoms following the ifrztyp
subroutine frozen_alpha(ifrztyp,ixyz,alpha,alphai)
  use module_base
  implicit none
  integer, intent(in) :: ifrztyp,ixyz
  real(gp), intent(in) :: alpha
  real(gp), intent(out) :: alphai
  !local variables
  logical :: move_this_coordinate

  if (move_this_coordinate(ifrztyp,ixyz)) then
     alphai=alpha
  else
     alphai=0.0_gp
  end if
 
end subroutine frozen_alpha

!routine for moving atomic positions, takes into account the 
!frozen atoms and the size of the cell
!synopsis: rxyz=txyz+alpha*sxyz
!all the shift are inserted into the box if there are periodic directions
!if the atom are frozen they are not moved
subroutine atomic_axpy(at,txyz,alpha,sxyz,rxyz)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: alpha
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,at%nat), intent(inout) :: rxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  do iat=1,at%nat
     !adjust the moving of the atoms following the frozen direction
     call frozen_alpha(at%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(at%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(at%ifrztyp(iat),3,alpha,alphaz)

     if (at%geocode == 'P') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),at%alat1)
        rxyz(2,iat)=modulo(txyz(2,iat)+alphay*sxyz(2,iat),at%alat2)
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),at%alat3)
     else if (at%geocode == 'S') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),at%alat1)
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),at%alat3)
     else
        rxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
     end if
  end do

end subroutine atomic_axpy

!routine for moving atomic positions, takes into account the 
!frozen atoms and the size of the cell
!synopsis: fxyz=txyz+alpha*sxyz
!update the forces taking into account the frozen atoms
!do not apply the modulo operation on forces
subroutine atomic_axpy_forces(at,txyz,alpha,sxyz,fxyz)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: alpha
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,at%nat), intent(inout) :: fxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz
  
  do iat=1,at%nat
     !adjust the moving of the forces following the frozen direction
     call frozen_alpha(at%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(at%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(at%ifrztyp(iat),3,alpha,alphaz)

     fxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
     fxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
     fxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
  end do
  
end subroutine atomic_axpy_forces


!calculate the scalar product between atomic positions by considering
!only non-blocked atoms
subroutine atomic_dot(at,x,y,scpr)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: x,y
  real(gp), intent(out) :: scpr
  !local variables
  integer :: iat
  real(gp) :: scpr1,scpr2,scpr3
  real(gp) :: alphax,alphay,alphaz

  scpr=0.0_gp

  do iat=1,at%nat
     call frozen_alpha(at%ifrztyp(iat),1,1.0_gp,alphax)
     call frozen_alpha(at%ifrztyp(iat),2,1.0_gp,alphay)
     call frozen_alpha(at%ifrztyp(iat),3,1.0_gp,alphaz)
     scpr1=alphax*x(1,iat)*y(1,iat)
     scpr2=alphay*x(2,iat)*y(2,iat)
     scpr3=alphaz*x(3,iat)*y(3,iat)
     scpr=scpr+scpr1+scpr2+scpr3
  end do
  
end subroutine atomic_dot

!z=alpha*A*x + beta* y
subroutine atomic_gemv(at,m,alpha,A,x,beta,y,z)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: m
  real(gp), intent(in) :: alpha,beta
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: x
  real(gp), dimension(m), intent(in) :: y
  real(gp), dimension(m,3,at%nat), intent(in) :: A
  real(gp), dimension(m), intent(out) :: z
  !local variables
  integer :: iat,i,j
  real(gp) :: mv,alphai
  
  do i=1,m
     mv=0.0_gp
     do iat=1,at%nat
        do j=1,3
           call frozen_alpha(at%ifrztyp(iat),j,A(i,j,iat),alphai)
           mv=mv+alphai*x(j,iat)
        end do
     end do
     z(i)=alpha*mv+beta*y(i)
  end do

end subroutine atomic_gemv

!the function which controls all the moving positions
function move_this_coordinate(ifrztyp,ixyz)
  use module_base
  implicit none
  integer, intent(in) :: ixyz,ifrztyp
  logical :: move_this_coordinate
  
  move_this_coordinate= &
       ifrztyp == 0 .or. &
       (ifrztyp == 2 .and. ixyz /=2) .or. &
       (ifrztyp == 3 .and. ixyz ==2)
       
end function move_this_coordinate

!synopsis: rxyz=txyz+alpha*sxyz
subroutine atomic_coordinate_axpy(at,ixyz,iat,t,alphas,r)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ixyz,iat
  real(gp), intent(in) :: t,alphas
  type(atoms_data), intent(in) :: at
  real(gp), intent(out) :: r
  !local variables
  logical :: move_this_coordinate,periodize
  real(gp) :: alat,alphai

  if (ixyz == 1) then
     alat=at%alat1
  else if (ixyz == 2) then
     alat=at%alat2
  else if (ixyz == 3) then
     alat=at%alat3
  end if
  
  periodize= at%geocode == 'P' .or. &
       (at%geocode == 'S' .and. ixyz /= 2)

  call frozen_alpha(at%ifrztyp(iat),ixyz,alphas,alphai)

  if (periodize) then
     r=modulo(t+alphai,alat)
  else
     r=t+alphai
  end if

end subroutine atomic_coordinate_axpy
