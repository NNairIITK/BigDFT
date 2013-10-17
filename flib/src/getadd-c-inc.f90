!> @file
!! Include fortran file for metadata adresses of characters
!! file included in module metadata_interfaces of getadd.f90
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!we should verify if allocatable array of dimension * is valid in most compilers
subroutine getc1(array,iadd)
  implicit none
  character(len=*), dimension(:), allocatable, intent(in) :: array
  integer(kind=8), intent(out) :: iadd
end subroutine getc1

!!$subroutine getc1_2(array,iadd)
!!$  implicit none
!!$  character(len=2), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_2
!!$
!!$subroutine getc1_3(array,iadd)
!!$  implicit none
!!$  character(len=3), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_3
!!$
!!$subroutine getc1_4(array,iadd)
!!$  implicit none
!!$  character(len=4), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_4
!!$
!!$subroutine getc1_5(array,iadd)
!!$  implicit none
!!$  character(len=5), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_5
!!$
!!$subroutine getc1_6(array,iadd)
!!$  implicit none
!!$  character(len=6), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_6
!!$
!!$subroutine getc1_7(array,iadd)
!!$  implicit none
!!$  character(len=7), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_7
!!$
!!$subroutine getc1_8(array,iadd)
!!$  implicit none
!!$  character(len=8), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_8
!!$
!!$subroutine getc1_9(array,iadd)
!!$  implicit none
!!$  character(len=9), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_9
!!$
!!$subroutine getc1_10(array,iadd)
!!$  implicit none
!!$  character(len=10), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_10
!!$
!!$subroutine getc1_11(array,iadd)
!!$  implicit none
!!$  character(len=11), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_11
!!$
!!$subroutine getc1_12(array,iadd)
!!$  implicit none
!!$  character(len=12), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_12
!!$
!!$subroutine getc1_13(array,iadd)
!!$  implicit none
!!$  character(len=13), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_13
!!$
!!$subroutine getc1_14(array,iadd)
!!$  implicit none
!!$  character(len=14), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_14
!!$
!!$subroutine getc1_15(array,iadd)
!!$  implicit none
!!$  character(len=15), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_15
!!$
!!$subroutine getc1_16(array,iadd)
!!$  implicit none
!!$  character(len=16), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_16
!!$
!!$subroutine getc1_17(array,iadd)
!!$  implicit none
!!$  character(len=17), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_17
!!$
!!$subroutine getc1_18(array,iadd)
!!$  implicit none
!!$  character(len=18), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_18
!!$
!!$subroutine getc1_19(array,iadd)
!!$  implicit none
!!$  character(len=19), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_19
!!$
!!$subroutine getc1_20(array,iadd)
!!$  implicit none
!!$  character(len=20), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_20
!!$
!!$subroutine getc1_21(array,iadd)
!!$  implicit none
!!$  character(len=21), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_21
!!$
!!$subroutine getc1_22(array,iadd)
!!$  implicit none
!!$  character(len=22), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_22
!!$
!!$subroutine getc1_23(array,iadd)
!!$  implicit none
!!$  character(len=23), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_23
!!$
!!$subroutine getc1_24(array,iadd)
!!$  implicit none
!!$  character(len=24), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_24
!!$
!!$subroutine getc1_25(array,iadd)
!!$  implicit none
!!$  character(len=25), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_25
!!$
!!$subroutine getc1_26(array,iadd)
!!$  implicit none
!!$  character(len=26), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_26
!!$
!!$subroutine getc1_27(array,iadd)
!!$  implicit none
!!$  character(len=27), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_27
!!$
!!$subroutine getc1_28(array,iadd)
!!$  implicit none
!!$  character(len=28), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_28
!!$
!!$subroutine getc1_29(array,iadd)
!!$  implicit none
!!$  character(len=29), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_29
!!$
!!$subroutine getc1_30(array,iadd)
!!$  implicit none
!!$  character(len=30), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_30
!!$
!!$subroutine getc1_31(array,iadd)
!!$  implicit none
!!$  character(len=31), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_31
!!$
!!$subroutine getc1_32(array,iadd)
!!$  implicit none
!!$  character(len=32), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_32
!!$
!!$subroutine getc1_33(array,iadd)
!!$  implicit none
!!$  character(len=33), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_33
!!$
!!$subroutine getc1_34(array,iadd)
!!$  implicit none
!!$  character(len=34), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_34
!!$
!!$subroutine getc1_35(array,iadd)
!!$  implicit none
!!$  character(len=35), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_35
!!$
!!$subroutine getc1_36(array,iadd)
!!$  implicit none
!!$  character(len=36), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_36
!!$
!!$subroutine getc1_37(array,iadd)
!!$  implicit none
!!$  character(len=37), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_37
!!$
!!$subroutine getc1_38(array,iadd)
!!$  implicit none
!!$  character(len=38), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_38
!!$
!!$subroutine getc1_39(array,iadd)
!!$  implicit none
!!$  character(len=39), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_39
!!$
!!$subroutine getc1_40(array,iadd)
!!$  implicit none
!!$  character(len=40), dimension(:), allocatable, intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$end subroutine getc1_40
!!$
