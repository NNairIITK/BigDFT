!>  Modules which contains information about fragment input variables
!! @author
!!    Copyright (C) 2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module fragment_base
  use public_keys
  use public_enums
  implicit none
  
  private

  !> Contains all parameters for the calculation of the fragments
  type, public :: fragmentInputParameters
     integer :: nfrag_ref
     integer :: nfrag
     integer, dimension(:), pointer :: frag_index         !< Array matching system fragments to reference fragments
     real(kind=8), dimension(:), pointer :: charge        !< Array giving the charge on each fragment for constrained DFT calculations
     !integer, dimension(:,:), pointer :: frag_info       !< Array giving number of atoms in fragment and environment for reference fragments
     character(len=100), dimension(:), pointer :: label   !< Array of fragment names
     character(len=100), dimension(:), pointer :: dirname !< Array of fragment directories, blank if not a fragment calculation
  end type fragmentInputParameters

  public :: allocateInputFragArrays,nullifyInputFragParameters,frag_from_dict
  public :: deallocateInputFragArrays,default_fragmentInputParameters,dict_from_frag

contains

  !> Nullify the parameters related to the fragments
  pure subroutine nullifyInputFragParameters(input_frag)
    implicit none

    ! Calling arguments
    type(fragmentInputParameters),intent(out) :: input_frag

    input_frag%nfrag_ref=-1
    input_frag%nfrag=-1
    !default scalar variables
    nullify(input_frag%frag_index)
    nullify(input_frag%charge)
    !nullify(input_frag%frag_info)
    nullify(input_frag%label)
    nullify(input_frag%dirname)

  end subroutine nullifyInputFragParameters

  subroutine frag_from_dict(dict,frag)
    use yaml_output, only: yaml_map,is_atoi
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    type(fragmentInputParameters), intent(out) :: frag
    !local variables
    integer :: frag_num,ncount,ncharged,ifrag
    character(len=max_field_length) :: frg_key
    type(dictionary), pointer :: dict_tmp,tmp2

    !now simulate the reading of the fragment structure from the dictionary
!!$  !to be continued after having fixed linear variables
    !iteration over the dictionary
    !count the number of reference fragments
    call nullifyInputFragParameters(frag)
    frag%nfrag_ref=0
    frag%nfrag=0
    !some sanity checks have to be added to this section
    dict_tmp=>dict_iter(dict)
    do while(associated(dict_tmp))
       select case(trim(dict_key(dict_tmp)))
       case(TRANSFER_INTEGRALS)
          !frag%calc_transfer_integrals=dict_tmp
       case(CONSTRAINED_DFT)
          ncharged=dict_size(dict_tmp)
          !constrained_dft=ncharged > 0
       case default
          frag%nfrag_ref=frag%nfrag_ref+1
          !count the number of fragments for this reference
          call count_local_fragments(dict_tmp,ncount)
          frag%nfrag=frag%nfrag+ncount
       end select
       dict_tmp=>dict_next(dict_tmp)
    end do

    if (frag%nfrag*frag%nfrag_ref == 0) then
       call f_err_throw('Fragment dictionary for the input is invalid',&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
       return
    end if
    call allocateInputFragArrays(frag)
    frag%charge=0.d0 !f_memset to be used

    dict_tmp=>dict_iter(dict)
    frag_num=0
    do while(associated(dict_tmp))
       select case(trim(dict_key(dict_tmp)))
       case(TRANSFER_INTEGRALS)
       case(CONSTRAINED_DFT)
          tmp2=>dict_iter(dict_tmp)
          !iterate over the charge TODO
          do while(associated(tmp2))
             frg_key=dict_key(tmp2)
             if (index(frg_key,FRAGMENT_NO) == 0 .or. &
                  .not. is_atoi(frg_key(len(FRAGMENT_NO)+1:))) then
                call f_err_throw('Invalid key in '//CONSTRAINED_DFT//&
                     ' section, key= '//trim(frg_key),&
                     err_name='BIGDFT_INPUT_VARIABLES_ERROR')
             end if
             read(frg_key(len(FRAGMENT_NO)+1:),*)ifrag
             frag%charge(ifrag)=tmp2
             tmp2=>dict_next(tmp2)
          end do
       case default
          frag_num=frag_num+1
          frag%label(frag_num)=repeat(' ',len(frag%label(frag_num)))
          frag%label(frag_num)=trim(dict_key(dict_tmp))
          !update directory name
          frag%dirname(frag_num)='data-'//trim(frag%label(frag_num))//'/'
          call count_local_fragments(dict_tmp,ncount,frag_index=frag%frag_index,frag_id=frag_num)
       end select
       dict_tmp=>dict_next(dict_tmp)
    end do

  contains

    subroutine count_local_fragments(dict_tmp,icount,frag_index,frag_id)
      use yaml_strings, only: is_atoi
      implicit none
      integer, intent(out) :: icount
      type(dictionary), pointer :: dict_tmp
      integer, dimension(:), intent(inout), optional :: frag_index
      integer, intent(in), optional :: frag_id

      !local variables
      integer :: idum,istart,i
      type(dictionary), pointer :: d_tmp
      character(len=max_field_length) :: val

      !iteration over the whole list
      icount=0
      istart=0
      d_tmp=>dict_iter(dict_tmp)
      do while(associated(d_tmp))
         val=d_tmp
         !if string is a integer consider it
         if (is_atoi(val)) then
            idum=d_tmp
            if (f_err_raise(istart>=idum,'error in entering fragment ids',&
                 err_name='BIGDFT_INPUT_VARIABLES_ERROR')) return
            if (istart /=0) then
               icount=icount+idum-istart
               if (present(frag_index) .and. present(frag_id)) then
                  do i=istart,idum
                     frag_index(i)=frag_id
                  end do
               end if
               istart=0
            else
               icount=icount+1
               if (present(frag_index) .and. present(frag_id)) frag_index(idum)=frag_id
            end if
         else if (f_err_raise(adjustl(trim(val))/='...',&
              'the only allowed values in the fragment list are integers or "..." string',&
              err_name='BIGDFT_INPUT_VARIABLES_ERROR')) then
            return
         else
            istart=idum
         end if
         d_tmp=>dict_next(d_tmp)
      end do
    end subroutine count_local_fragments

  end subroutine frag_from_dict

  !> Allocate the arrays for the input related to the fragment
  subroutine allocateInputFragArrays(input_frag)
    use dynamic_memory
    implicit none

    ! Calling arguments
    type(fragmentInputParameters),intent(inout) :: input_frag

    ! Local variables
    character(len=*),parameter :: subname='allocateInputFragArrays'

    input_frag%frag_index = f_malloc_ptr(input_frag%nfrag,id='input_frag%frag_index')
    input_frag%charge = f_malloc_ptr(input_frag%nfrag,id='input_frag%charge')

    input_frag%label=f_malloc_str_ptr(int(len(input_frag%label),kind=4),&
         input_frag%nfrag_ref,id='input_frag%label')
    !f_malloc0_str_ptr should be used here
    input_frag%dirname=f_malloc_str_ptr(int(len(input_frag%dirname),kind=4),&
         input_frag%nfrag_ref,id='input_frag%label')

    !set the variables to their default value

  end subroutine allocateInputFragArrays


  !> Deallocate the arrays related to the input for the fragments
  subroutine deallocateInputFragArrays(input_frag)
    use dynamic_memory
    implicit none

    ! Calling arguments
    type(fragmentInputParameters),intent(inout) :: input_frag

    ! Local variables
    character(len=*),parameter :: subname='deallocateInputFragArrays'

    call f_free_ptr(input_frag%frag_index)
    call f_free_ptr(input_frag%charge)

    call f_free_str_ptr(len(input_frag%label),input_frag%label)
    call f_free_str_ptr(len(input_frag%dirname),input_frag%dirname)

  end subroutine deallocateInputFragArrays

  !> initialize fragment input parameters to their default value
  subroutine default_fragmentInputParameters(frag)
    implicit none
    type(fragmentInputParameters),intent(out) :: frag

    !first nullify
    call nullifyInputFragParameters(frag)
    !then set defaults
    frag%nfrag_ref=1
    frag%nfrag=1
    !then allocate
    call allocateInputFragArrays(frag)
    !and fill to neutral values
    frag%label(1)=repeat(' ',len(frag%label))
    frag%dirname(1)=repeat(' ',len(frag%dirname))
    frag%frag_index(1)=1
    frag%charge(1)=0.0d0

  end subroutine default_fragmentInputParameters

  !> routine to build dictionary of fragment for purposes of backward compatibility with the old format
  subroutine dict_from_frag(frag,dict_frag)
    use yaml_strings, only: yaml_toa
    use dictionaries, dict_set => set
    implicit none
    type(fragmentInputParameters), intent(in) :: frag
    type(dictionary), pointer :: dict_frag
    !local variables
    integer :: ifrag
    type(dictionary), pointer :: frag_list

    !create a dictionary with the given information
    call dict_init(dict_frag)
    do ifrag=1,frag%nfrag_ref
       !build the list of fragments associated to the reference
       call build_frag_list(ifrag,frag%nfrag,frag%frag_index,frag_list)
       call dict_set(dict_frag//trim(frag%label(ifrag)),frag_list)
    end do
    !then set the constrained DFT elements in case there is a charge
    do ifrag=1,frag%nfrag
       if (frag%charge(ifrag) /= 0.d0) &
            call dict_set(dict_frag//CONSTRAINED_DFT//(FRAGMENT_NO//&
            trim(adjustl(yaml_toa(ifrag)))),frag%charge(ifrag),fmt='(f7.2)')
    end do

  contains

    subroutine build_frag_list(ifrag_ref,nfrag,frag_index,frag_list)
      implicit none
      integer, intent(in) :: nfrag,ifrag_ref
      integer, dimension(nfrag), intent(in) :: frag_index
      type(dictionary), pointer :: frag_list
      !local variables
      logical :: agree,willagree
      integer :: ifrag,istart,iend

      call dict_init(frag_list)

      !for each segment find the starting and ending point 
      ifrag=1
      istart=0
      iend=0
      find_segment: do while(ifrag <= nfrag)
         agree=frag_index(ifrag)== ifrag_ref
         willagree= ifrag < nfrag
         if (willagree) willagree=frag_index(ifrag+1) == ifrag_ref
         if (agree) then
            !starting point if it has not been found yet
            if (istart==0) then
               istart=ifrag
            end if
            if (willagree) then
               iend=0
            else
               iend=ifrag
            end if
         end if
         if (istart*iend /= 0) then
            if (istart==iend) then
               call add(frag_list,istart)
            else if (iend == istart+1) then
               call add(frag_list,istart)
               call add(frag_list,iend)
            else
               call add(frag_list,istart)
               call add(frag_list,'...')
               call add(frag_list,iend)
            end if
            istart=0
            iend=0
         end if
         ifrag=ifrag+1
      end do find_segment
    end subroutine build_frag_list
  end subroutine dict_from_frag

end module fragment_base
