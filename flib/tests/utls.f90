!> @file
!! Routine to tests f_utils module
!! @example utls.f90
!! Examples using the @ref f_utils objects (units and timers)
!! @author
!!    Copyright (C) 2013-2015 BigDFT group
!!    This file is distributed oneder the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine f_utils_test()
  use f_precisions
  use f_utils
  use yaml_output
  use f_enums
  use dictionaries, only: f_loc
  implicit none
  !local variables
  type(f_enumerator) :: greetings=f_enumerator('Greetings',10) 
  type(f_enumerator) :: f1=f_enumerator('Ciao',1)              
  type(f_enumerator) :: f2=f_enumerator('Hello',2)             
  type(f_enumerator) :: f3=f_enumerator('Pizza',3)             
  integer :: unt,unt2,u
!  double precision :: t0
  integer(kind=8) :: i0
  real(f_simple), dimension(3) :: r1
  real(f_double), dimension(3) :: r2
  real(f_quadruple), dimension(3) :: r4
  complex(f_simple), dimension(3) :: c1
  complex(f_double), dimension(3) :: c2
  complex(f_quadruple), dimension(3) :: c4
  integer(f_short), dimension(3) :: is
  integer(f_integer), dimension(3) :: i4
  integer(f_long), dimension(3) :: il
  logical(f_byte), dimension(3) :: lb
  character(len=256) :: path
  logical, dimension(3) :: l

  r4=real(10.0,f_quadruple)

  call yaml_map('Long Integer kind',f_long)
  call yaml_map('Normal Integer kind',f_integer)
  call yaml_map('Short Integer kind',f_short)

  call yaml_map('Quadruple precision Real kind',f_quadruple)
  call yaml_map('Double precision Real kind',f_integer)
  call yaml_map('Single precision Real kind',f_short)


!  call expq(r4(1),r4(2))

  r2=10.d0
!  call expq(r2(1),r2(2))

  !check the precisions
  call yaml_map('Kinds of reals',[f_loc(r1(2))-f_loc(r1(1)),f_loc(r2(2))-f_loc(r2(1)),f_loc(r4(2))-f_loc(r4(1))])
  call yaml_map('Kinds of complex',[f_loc(c1(2))-f_loc(c1(1)),f_loc(c2(2))-f_loc(c2(1)),f_loc(c4(2))-f_loc(c4(1))])
  call yaml_map('Kinds of integer',[f_loc(is(2))-f_loc(is(1)),f_loc(i4(2))-f_loc(i4(1)),f_loc(il(2))-f_loc(il(1))])
  call yaml_map('Kinds of logicals',[f_loc(lb(2))-f_loc(lb(1)),f_loc(l(2))-f_loc(l(1))])

!  print *,'test',r4(1),r4(2),'scale',scale(100.0_quadruple,2),scale(100.0_quadruple,4)!,exp(r4(1))
!  print *,'test',r2(1),r2(2),exp(r2(1)),exponent(r2(2)),scale(r2(2),exponent(r2(2)))!,scale(r2(1))

  !see if the attribute can live outside a given scope
  call f_enum_attr(f1,attr=greetings)
  call f_enum_attr(f2,attr=greetings)

  call yaml_map('Enum1 char',str(f1))
  call yaml_map('Enum1 int',int(f1))
  call yaml_map('Enum1 check',f1=='Ciao')

  call yaml_map('Greetings 1a',f1 .hasattr. 'Greetings') 
  call yaml_map('Greetings 1b',f1 .hasattr. 10)
  call yaml_map('Greetings 1c',f1 .hasattr. greetings) 

  call yaml_map('Greetings 2a',f2 .hasattr. 'Greetings') 
  call yaml_map('Greetings 2b',f2 .hasattr. 10)
  call yaml_map('Greetings 2c',f2 .hasattr. greetings) 

  call yaml_map('Greetings 3a',f3 .hasattr. 'Greetings') 
  call yaml_map('Greetings 3b',f3 .hasattr. 10)
  call yaml_map('Greetings 3c',f3 .hasattr. greetings) 

  !wait one second
  !t0=dble(f_time())*1.d-9
  i0=f_time()
  call yaml_map('Absolute time before pause (since epoch)',yaml_walltime_toa(i0))
  call f_pause(1)
  call yaml_map('Time spent after pause (s)',yaml_walltime_toa(f_time()-i0))!dble(f_time())*1.d-9-t0)
  !call yaml_map('Conversion of walltimes in standard form',yaml_walltime_toa(f_time()))
  !open files and get free_units
  unt=f_get_free_unit()
  call yaml_map('First unit which can be opened',unt)
  call f_open_file(unt,'Testopen')
  call yaml_map('Opened file in unit',unt)
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Second unit which can be opened',unt2)
  call f_open_file(unt2,'Testopen2')
  call yaml_map('Opened file in unit',unt2)
  call f_close(unt)
  call f_delete_file('Testopen')
  call f_delete_file('Testopen2')
  !again in another unit
  unt2=f_get_free_unit()
  call yaml_map('Third unit which can be opened',unt2)

  call yaml_mapping_open('Bastian test')
  u=f_get_free_unit()                                                   
  call yaml_map('First file',u)
  open(unit=u,file='test') 
  !call f_open_file(u,'test')
  call yaml_map('First file opened',u)
  unt=f_get_free_unit()                                                 
  call yaml_map('Second file',unt)
  open(unit=unt,file='test2')
  call yaml_map('Second file opened',unt)
  unt2=f_get_free_unit()                                      
  call yaml_map('Third file',unt2)
  open(unit=unt2,file='test3')     
  call yaml_map('Third file opened',unt2)
  call f_delete_file('test')
  call f_delete_file('test2')
  call f_delete_file('test3')
  call yaml_mapping_close()
  call yaml_map('If this value is 7 then all files have been correctly closed',f_get_free_unit())

  !create a directory (we should add the test to remove it)
  !call f_mkdir('testdir',path)

end subroutine f_utils_test

!!$subroutine expq(x,res)
!!$  use f_precisions
!!$  implicit none
!!$  real(quadruple), intent(in) :: x
!!$  real(quadruple), intent(out) :: res
!!$  !local variables
!!$  real(quadruple), parameter :: log2=log(2.0_quadruple)
!!$  integer(long) :: i
!!$  real(quadruple) :: r,a
!!$
!!$  !easy case
!!$  if (x==real(0.0,quadruple)) then
!!$     res=real(1.0,quadruple)
!!$  else
!!$     !first determine the rest
!!$     r=x/log2
!!$     a=aint(r,quadruple)
!!$     r=r-a
!!$     r=2**r
!!$     i=int(a,long)
!!$     res=scale(r,i)
!!$  end if
!!$
!!$!  print *,'res',res,exp(x)
!!$  
!!$end subroutine expq

subroutine f_inputfile_test()
  use dictionaries
  use yaml_output
  use f_input_file
  use f_precisions, only: f_cr
  implicit none
  
  !> input definitions as they can be defined by the developers
  character(len=*), parameter :: inputdef='                                                         '//f_cr//&
       '  dft:                                                                                      '//f_cr//&
       '    DESCRIPTION: Density Functional Theory parameters                                       '//f_cr//&
       '    hgrids: #a variable with range                                                          '//f_cr//&
       '      COMMENT: Grid spacing in the three directions (bohr)                                  '//f_cr//&
       '      DESCRIPTION: |                                                                        '//f_cr//&
       '       Grid spacing in three directions (Bohr units) of the coarse mesh.  '//f_cr//&
       '       A scalar can also be given as 0.45.                                                  '//f_cr//&
       '      RANGE: [0., 2.]                                                                       '//f_cr//&
       '      default: [0.45, 0.45, 0.45]                                                           '//f_cr//&
       '      fast: [0.55, 0.55, 0.55]                                                              '//f_cr//&
       '      accurate: [0.30, 0.30, 0.30]                                                          '//f_cr//&
       '    ixc: #a variable with several profiles                                                  '//f_cr//&
       '      COMMENT: Exchange-correlation parameter (LDA=1,PBE=11)                                '//f_cr//&
       '      DESCRIPTION: Determine the exchange-correlation functional.                           '//f_cr//&
       '      default: 1                                                                            '//f_cr//&
       '      #Here follow a number of possibilities for the different XC functionals               '//f_cr//&
       '      LDA (ABINIT): 1                                                                       '//f_cr//&
       '      PBE (ABINIT): 11                                                                      '//f_cr//&
       '      LDA: -20                                                                              '//f_cr//&
       '      PBE: -101130                                                                          '//f_cr//&
       '      PBE0: -406                                                                            '//f_cr//&
       '      B97-D: -170 #to be verified                                                           '//f_cr//&
       '      B3LYP: -402                                                                           '//f_cr//&
       '      HF: 100 #Hartree-Fock                                                                 '//f_cr//&
       '  geopt:                                                                                    '//f_cr//&
       '    DESCRIPTION: Parameters for the geometry relaxation and molecular dynamics              '//f_cr//&
       '    method:  # a variable with exclusive                                                    '//f_cr//&
       '      COMMENT: Geometry optimisation method                                                 '//f_cr//&
       '      EXCLUSIVE:                                                                            '//f_cr//&
       '        none:   No geometry optimization                                                    '//f_cr//&
       '        SDCG:   A combination of Steepest Descent and Conjugate Gradient                    '//f_cr//&
       '        VSSD:   Variable Stepsize Steepest Descent method                                   '//f_cr//&
       '        LBFGS:  Limited-memory BFGS                                                         '//f_cr//&
       '        BFGS:   Broyden-Fletcher-Goldfarb-Shanno                                            '//f_cr//&
       '        PBFGS:  Same as BFGS with an initial Hessian obtained from a force field            '//f_cr//&
       '        AB6MD:  Molecular dynamics from ABINIT                                              '//f_cr//&
       '        DIIS:   Direct inversion of iterative subspace                                      '//f_cr//&
       '        FIRE:   Fast Inertial Relaxation Engine as described by Bitzek et al.               '//f_cr//&
       '        NEB:    Nudged Elastic Band                                                         '//f_cr//&
       '        SBFGS:  SQNM minimizer, keyword deprecated, will be replaced by SQNM in future relea'//f_cr//&
       '        SQNM:   Stabilized quasi-Newton minimzer                                            '//f_cr//&
       '      default: none                                                                         '//f_cr//&
       '    ncount_cluster_x:                                                                       '//f_cr//&
       '      COMMENT: Maximum number of force evaluations                                          '//f_cr//&
       '      RANGE: [0, 2000]                                                                      '//f_cr//&
       '      PROFILE_FROM: method                                                                  '//f_cr//&
       '      default: 50                                                                           '//f_cr//&
       '      none: 1                                                                               '//f_cr//&
       '    betax:                                                                                  '//f_cr//&
       '      COMMENT: Stepsize for the geometry optimization                                       '//f_cr//&
       '      CONDITION:                                                                            '//f_cr//&
       '        MASTER_KEY: method                                                                  '//f_cr//&
       '        WHEN:                                                                               '//f_cr//&
       '        - SDCG                                                                              '//f_cr//&
       '        - VSSD                                                                              '//f_cr//&
       '        - LBFGS                                                                             '//f_cr//&
       '        - BFGS                                                                              '//f_cr//&
       '        - PBFGS                                                                             '//f_cr//&
       '        - DIIS                                                                              '//f_cr//&
       '        - FIRE                                                                              '//f_cr//&
       '        - NEB                                                                               '//f_cr//&
       '        - SBFGS                                                                             '//f_cr//&
       '        - SQNM                                                                              '//f_cr//&
       '        - none                                                                              '//f_cr//&
       '      PROFILE_FROM: method                                                                  '//f_cr//&
       '      RANGE: [0., 100.]                                                                     '//f_cr//&
       '      default: 4.                                                                           '//f_cr//&
       '      DIIS: 2.                                                                              '//f_cr//&
       '      NEB: 0.5                                                                              '//f_cr

  !> dictionary of the input definitions
  type(dictionary), pointer :: inputdefinitions
  !> dictionary of the imports
  type(dictionary), pointer :: dict_profiles
  !> input file
  type(dictionary), pointer :: input
  !> minuimal input file
  type(dictionary), pointer :: input_minimal
  type(dictionary), pointer :: as_is,nested
  character(len=*), parameter :: example1='                                                         '//f_cr//&
       'dft:             '//f_cr//&
       ' hgrids: 0.45    '//f_cr//&
       ' ixc: B3LYP     '//f_cr//&
       ' bidon: 2       '//f_cr//&
       'geopt:           '//f_cr//&
       ' method: DIIS   '//f_cr
  character(len=*), parameter :: profiles='                                                         '//f_cr//&
       ' simple:  '//f_cr//&
       '   dft:             '//f_cr//&
       '     hgrids: 0.45    '//f_cr//&
       '     ixc: LDA     '//f_cr//&
       ' geopt:'//f_cr//&
       '    geopt: {method: DIIS} '//f_cr
  character(len=*), parameter :: example2='import: [geopt, simple]'
  

  !first, insert the definitiions of the input file in the dictionary
  call parse_dict(inputdefinitions,inputdef)
  call yaml_map('Initial input variables',inputdefinitions)
  !then parse the user's input file
  call parse_dict(input,example1)

  call yaml_map('User input file',input)

  !complete input file, should crash now
  call f_err_open_try()
  call input_file_complete(inputdefinitions,input)
  if (f_err_check()) then
     call f_dump_all_errors(-1)
     call dict_remove(input//'dft','bidon')
  end if
  !now try again
  call input_file_complete(inputdefinitions,input)
  call f_err_close_try()

  call yaml_map('Completed input file',input)
  call input_file_dump(input)

  !only user defined
  call input_file_dump(input,.true.)

  !retrieve input minimal
  nullify(as_is,nested)
  call input_file_minimal(inputdefinitions,input,input_minimal,nested,as_is)
  call yaml_map('Input minimal',input_minimal)

  !and redo the same stuff
  call input_file_complete(inputdefinitions,input_minimal)
  call yaml_map('Completed minimal input file',input_minimal)

  call dict_free(input,input_minimal)

  !now redo same thing for example2
  !then parse the user's input file
  call parse_dict(dict_profiles,profiles)

  !now redo same thing for example2
  !then parse the user's input file
  call parse_dict(input,example2)

  call yaml_map('User input file 2',input)
  call input_file_complete(inputdefinitions,input,imports=dict_profiles)
  
  call yaml_map('Completed input file 2',input)
  call input_file_dump(input)
  !retrieve input minimal
  nullify(as_is,nested)
  call input_file_minimal(inputdefinitions,input,input_minimal,nested,as_is)
  call yaml_map('Input minimal 2',input_minimal)

  call dict_free(inputdefinitions,input,input_minimal,dict_profiles)

contains
  subroutine parse_dict(dict,str)
    use yaml_parse
    implicit none
    type(dictionary), pointer :: dict
    character(len=*), intent(in) :: str
    !local variales
    type(dictionary), pointer :: tmp


    call yaml_parse_from_string(tmp,str)
    dict => tmp .pop. 0
    call dict_free(tmp)

  end subroutine parse_dict
end subroutine f_inputfile_test
