!!****m* ABINIT/interfaces_linalg
!! NAME
!! interfaces_linalg
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory prereqs/linalg
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_linalg

 implicit none

interface
 subroutine caxpy(n,ca,cx,incx,cy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  complex :: ca
  complex :: cx(*)
  complex :: cy(*)
 end subroutine caxpy
end interface

interface
 subroutine  ccopy(n,cx,incx,cy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  complex :: cx(*)
  complex :: cy(*)
 end subroutine ccopy
end interface

interface
 complex function cdotc(n,cx,incx,cy,incy)
 implicit none
 integer :: incx
 integer :: incy
 integer :: n
 complex :: cx(*)
 complex :: cy(*)
end function cdotc
end interface

interface
 complex function cdotu(n,cx,incx,cy,incy)
 implicit none
 integer :: incx
 integer :: incy
 integer :: n
 complex :: cx(*)
 complex :: cy(*)
end function cdotu
end interface

interface
 SUBROUTINE CGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: KL
  integer :: KU
  integer :: LDA
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: TRANS
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CGBMV
end interface

interface
 SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: TRANSA
  character*1 :: TRANSB
  complex :: A( LDA, * )
  complex :: B( LDB, * )
  complex :: C( LDC, * )
 end subroutine CGEMM
end interface

interface
 SUBROUTINE CGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: TRANS
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CGEMV
end interface

interface
 SUBROUTINE CGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CGERC
end interface

interface
 SUBROUTINE CGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CGERU
end interface

interface
 SUBROUTINE CHBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: K
  integer :: LDA
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: UPLO
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CHBMV
end interface

interface
 SUBROUTINE CHEMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: SIDE
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: B( LDB, * )
  complex :: C( LDC, * )
 end subroutine CHEMM
end interface

interface
 SUBROUTINE CHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: UPLO
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CHEMV
end interface

interface
 SUBROUTINE CHER  ( UPLO, N, ALPHA, X, INCX, A, LDA )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  real :: ALPHA
  character*1 :: UPLO
  complex :: X( * )
  complex :: A( LDA, * )
 end subroutine CHER
end interface

interface
 SUBROUTINE CHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  complex :: ALPHA
  character*1 :: UPLO
  complex :: X( * )
  complex :: Y( * )
  complex :: A( LDA, * )
 end subroutine CHER2
end interface

interface
 SUBROUTINE CHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: N
  complex :: ALPHA
  real :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: B( LDB, * )
  complex :: C( LDC, * )
 end subroutine CHER2K
end interface

interface
 SUBROUTINE CHERK ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: C( LDC, * )
 end subroutine CHERK
end interface

interface
 SUBROUTINE CHPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: UPLO
  complex :: AP( * )
  complex :: X( * )
  complex :: Y( * )
 end subroutine CHPMV
end interface

interface
 SUBROUTINE CHPR  ( UPLO, N, ALPHA, X, INCX, AP )
  implicit none
  integer :: INCX
  integer :: N
  real :: ALPHA
  character*1 :: UPLO
  complex :: AP( * )
  complex :: X( * )
 end subroutine CHPR
end interface

interface
 SUBROUTINE CHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  complex :: ALPHA
  character*1 :: UPLO
  complex :: AP( * )
  complex :: X( * )
  complex :: Y( * )
 end subroutine CHPR2
end interface

interface
 subroutine crotg(ca,cb,c,s)
  implicit none
  real :: c
  complex :: ca
  complex :: cb
  complex :: s
 end subroutine crotg
end interface

interface
 subroutine  cscal(n,ca,cx,incx)
  implicit none
  integer :: incx
  integer :: n
  complex :: ca
  complex :: cx(*)
 end subroutine cscal
end interface

interface
 subroutine  csrot (n,cx,incx,cy,incy,c,s)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  real :: c
  real :: s
  complex :: cx(1)
  complex :: cy(1)
 end subroutine csrot
end interface

interface
 subroutine  csscal(n,sa,cx,incx)
  implicit none
  integer :: incx
  integer :: n
  real :: sa
  complex :: cx(*)
 end subroutine csscal
end interface

interface
 subroutine  cswap (n,cx,incx,cy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  complex :: cx(*)
  complex :: cy(*)
 end subroutine cswap
end interface

interface
 SUBROUTINE CSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: SIDE
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: B( LDB, * )
  complex :: C( LDC, * )
 end subroutine CSYMM
end interface

interface
 SUBROUTINE CSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: B( LDB, * )
  complex :: C( LDC, * )
 end subroutine CSYR2K
end interface

interface
 SUBROUTINE CSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: N
  complex :: ALPHA
  complex :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: C( LDC, * )
 end subroutine CSYRK
end interface

interface
 SUBROUTINE CTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: K
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: X( * )
  complex :: A( LDA, * )
 end subroutine CTBMV
end interface

interface
 SUBROUTINE CTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: K
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: X( * )
  complex :: A( LDA, * )
 end subroutine CTBSV
end interface

interface
 SUBROUTINE CTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: AP( * )
  complex :: X( * )
 end subroutine CTPMV
end interface

interface
 SUBROUTINE CTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: AP( * )
  complex :: X( * )
 end subroutine CTPSV
end interface

interface
 SUBROUTINE CTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  complex :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: B( LDB, * )
 end subroutine CTRMM
end interface

interface
 SUBROUTINE CTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: X( * )
  complex :: A( LDA, * )
 end subroutine CTRMV
end interface

interface
 SUBROUTINE CTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  complex :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  complex :: A( LDA, * )
  complex :: B( LDB, * )
 end subroutine CTRSM
end interface

interface
 SUBROUTINE CTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex :: X( * )
  complex :: A( LDA, * )
 end subroutine CTRSV
end interface

interface
 double precision function dasum(n,dx,incx)
 implicit none
 integer :: incx
 integer :: n
 double precision :: dx(*)
end function dasum
end interface

interface
 subroutine daxpy(n,da,dx,incx,dy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  double precision :: da
  double precision :: dx(*)
  double precision :: dy(*)
 end subroutine daxpy
end interface

interface
 double precision function dcabs1(z)
 implicit none
 complex(kind=8) :: z
end function dcabs1
end interface

interface
 subroutine  dcopy(n,dx,incx,dy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  double precision :: dx(*)
  double precision :: dy(*)
 end subroutine dcopy
end interface

interface
 double precision function ddot(n,dx,incx,dy,incy)
 implicit none
 integer :: incx
 integer :: incy
 integer :: n
 double precision :: dx(*)
 double precision :: dy(*)
end function ddot
end interface

interface
 SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: KL
  integer :: KU
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: ALPHA
  double precision :: BETA
  character*1 :: TRANS
  double precision :: X( * )
  double precision :: Y( * )
  double precision :: A( LDA, * )
 end subroutine DGBMV
end interface

interface
 SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  double precision :: ALPHA
  double precision :: BETA
  character*1 :: TRANSA
  character*1 :: TRANSB
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
  double precision :: C( LDC, * )
 end subroutine DGEMM
end interface

interface
 SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: ALPHA
  double precision :: BETA
  character*1 :: TRANS
  double precision :: X( * )
  double precision :: Y( * )
  double precision :: A( LDA, * )
 end subroutine DGEMV
end interface

interface
 SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: ALPHA
  double precision :: X( * )
  double precision :: Y( * )
  double precision :: A( LDA, * )
 end subroutine DGER
end interface

interface
 DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
 implicit none
 integer :: INCX
 integer :: N
 double precision :: X( * )
end function DNRM2
end interface

interface
 subroutine  drot (n,dx,incx,dy,incy,c,s)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  double precision :: c
  double precision :: s
  double precision :: dx(*)
  double precision :: dy(*)
 end subroutine drot
end interface

interface
 subroutine drotg(da,db,c,s)
  implicit none
  double precision :: c
  double precision :: da
  double precision :: db
  double precision :: s
 end subroutine drotg
end interface

interface
 SUBROUTINE DROTM (N,DX,INCX,DY,INCY,DPARAM)
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  double precision :: DPARAM(5)
  double precision :: DX(1)
  double precision :: DY(1)
 end subroutine DROTM
end interface

interface
 SUBROUTINE DROTMG (DD1,DD2,DX1,DY1,DPARAM)
  implicit none
  double precision :: DD1
  double precision :: DD2
  double precision :: DX1
  double precision :: DY1
  double precision :: DPARAM(5)
 end subroutine DROTMG
end interface

interface
 SUBROUTINE DSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: K
  integer :: LDA
  integer :: N
  double precision :: ALPHA
  double precision :: BETA
  character*1 :: UPLO
  double precision :: X( * )
  double precision :: Y( * )
  double precision :: A( LDA, * )
 end subroutine DSBMV
end interface

interface
 subroutine  dscal(n,da,dx,incx)
  implicit none
  integer :: incx
  integer :: n
  double precision :: da
  double precision :: dx(*)
 end subroutine dscal
end interface

interface
 DOUBLE PRECISION FUNCTION DSDOT (N, SX, INCX, SY, INCY)
 implicit none
 integer :: INCX
 integer :: INCY
 integer :: N
 real :: SX(*)
 real :: SY(*)
end function DSDOT
end interface

interface
 subroutine  dswap (n,dx,incx,dy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  double precision :: dx(*)
  double precision :: dy(*)
 end subroutine dswap
end interface

interface
 SUBROUTINE DSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  double precision :: ALPHA
  double precision :: BETA
  character*1 :: SIDE
  character*1 :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
  double precision :: C( LDC, * )
 end subroutine DSYMM
end interface

interface
 SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  double precision :: ALPHA
  double precision :: BETA
  character*1 :: UPLO
  double precision :: X( * )
  double precision :: Y( * )
  double precision :: A( LDA, * )
 end subroutine DSYMV
end interface

interface
 SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  double precision :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DTRMM
end interface

interface
 SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  double precision :: X( * )
  double precision :: A( LDA, * )
 end subroutine DTRMV
end interface

interface
 SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  double precision :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DTRSM
end interface

interface
 SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  double precision :: X( * )
  double precision :: A( LDA, * )
 end subroutine DTRSV
end interface

interface
 double precision function dzasum(n,zx,incx)
 implicit none
 integer :: incx
 integer :: n
 complex(kind=8) :: zx(*)
end function dzasum
end interface

interface
 DOUBLE PRECISION FUNCTION DZNRM2( N, X, INCX )
 implicit none
 integer :: INCX
 integer :: N
 complex(kind=8) :: X( * )
end function DZNRM2
end interface

interface
 integer function icamax(n,cx,incx)
 implicit none
 integer :: incx
 integer :: n
 complex :: cx(*)
end function icamax
end interface

interface
 integer function idamax(n,dx,incx)
 implicit none
 integer :: incx
 integer :: n
 double precision :: dx(*)
end function idamax
end interface

interface
 integer function isamax(n,sx,incx)
 implicit none
 integer :: incx
 integer :: n
 real :: sx(*)
end function isamax
end interface

interface
 integer function izamax(n,zx,incx)
 implicit none
 integer :: incx
 integer :: n
 complex(kind=8) :: zx(*)
end function izamax
end interface

interface
 real function sasum(n,sx,incx)
 implicit none
 integer :: incx
 integer :: n
 real :: sx(*)
end function sasum
end interface

interface
 subroutine saxpy(n,sa,sx,incx,sy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  real :: sa
  real :: sx(*)
  real :: sy(*)
 end subroutine saxpy
end interface

interface
 real function scasum(n,cx,incx)
 implicit none
 integer :: incx
 integer :: n
 complex :: cx(*)
end function scasum
end interface

interface
 REAL             FUNCTION SCNRM2( N, X, INCX )
 implicit none
 integer :: INCX
 integer :: N
 complex :: X( * )
end function SCNRM2
end interface

interface
 subroutine scopy(n,sx,incx,sy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  real :: sx(*)
  real :: sy(*)
 end subroutine scopy
end interface

interface
 real function sdot(n,sx,incx,sy,incy)
 implicit none
 integer :: incx
 integer :: incy
 integer :: n
 real :: sx(*)
 real :: sy(*)
end function sdot
end interface

interface
 REAL FUNCTION SDSDOT (N, SB, SX, INCX, SY, INCY)
 implicit none
 integer :: INCX
 integer :: INCY
 integer :: N
 real :: SB
 real :: SX(*)
 real :: SY(*)
end function SDSDOT
end interface

interface
 SUBROUTINE SGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: KL
  integer :: KU
  integer :: LDA
  integer :: M
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: TRANS
  real :: X( * )
  real :: Y( * )
  real :: A( LDA, * )
 end subroutine SGBMV
end interface

interface
 SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: TRANSA
  character*1 :: TRANSB
  real :: A( LDA, * )
  real :: B( LDB, * )
  real :: C( LDC, * )
 end subroutine SGEMM
end interface

interface
 SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: TRANS
  real :: X( * )
  real :: Y( * )
  real :: A( LDA, * )
 end subroutine SGEMV
end interface

interface
 SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  real :: ALPHA
  real :: X( * )
  real :: Y( * )
  real :: A( LDA, * )
 end subroutine SGER
end interface

interface
 REAL             FUNCTION SNRM2 ( N, X, INCX )
 implicit none
 integer :: INCX
 integer :: N
 real :: X( * )
end function SNRM2
end interface

interface
 subroutine srot (n,sx,incx,sy,incy,c,s)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  real :: c
  real :: s
  real :: sx(*)
  real :: sy(*)
 end subroutine srot
end interface

interface
 subroutine srotg(sa,sb,c,s)
  implicit none
  real :: c
  real :: s
  real :: sa
  real :: sb
 end subroutine srotg
end interface

interface
 SUBROUTINE SROTM (N,SX,INCX,SY,INCY,SPARAM)
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  real :: SPARAM(5)
  real :: SX(1)
  real :: SY(1)
 end subroutine SROTM
end interface

interface
 SUBROUTINE SROTMG (SD1,SD2,SX1,SY1,SPARAM)
  implicit none
  real :: SD1
  real :: SD2
  real :: SX1
  real :: SY1
  real :: SPARAM(5)
 end subroutine SROTMG
end interface

interface
 SUBROUTINE SSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: K
  integer :: LDA
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: UPLO
  real :: X( * )
  real :: Y( * )
  real :: A( LDA, * )
 end subroutine SSBMV
end interface

interface
 subroutine sscal(n,sa,sx,incx)
  implicit none
  integer :: incx
  integer :: n
  real :: sa
  real :: sx(*)
 end subroutine sscal
end interface

interface
 SUBROUTINE SSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: UPLO
  real :: AP( * )
  real :: X( * )
  real :: Y( * )
 end subroutine SSPMV
end interface

interface
 SUBROUTINE SSPR  ( UPLO, N, ALPHA, X, INCX, AP )
  implicit none
  integer :: INCX
  integer :: N
  real :: ALPHA
  character*1 :: UPLO
  real :: AP( * )
  real :: X( * )
 end subroutine SSPR
end interface

interface
 SUBROUTINE SSPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  real :: ALPHA
  character*1 :: UPLO
  real :: AP( * )
  real :: X( * )
  real :: Y( * )
 end subroutine SSPR2
end interface

interface
 subroutine sswap (n,sx,incx,sy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  real :: sx(*)
  real :: sy(*)
 end subroutine sswap
end interface

interface
 SUBROUTINE SSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: SIDE
  character*1 :: UPLO
  real :: A( LDA, * )
  real :: B( LDB, * )
  real :: C( LDC, * )
 end subroutine SSYMM
end interface

interface
 SUBROUTINE SSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: UPLO
  real :: X( * )
  real :: Y( * )
  real :: A( LDA, * )
 end subroutine SSYMV
end interface

interface
 SUBROUTINE SSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  real :: ALPHA
  character*1 :: UPLO
  real :: X( * )
  real :: A( LDA, * )
 end subroutine SSYR
end interface

interface
 SUBROUTINE SSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  real :: ALPHA
  character*1 :: UPLO
  real :: X( * )
  real :: Y( * )
  real :: A( LDA, * )
 end subroutine SSYR2
end interface

interface
 SUBROUTINE SSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  real :: A( LDA, * )
  real :: B( LDB, * )
  real :: C( LDC, * )
 end subroutine SSYR2K
end interface

interface
 SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: N
  real :: ALPHA
  real :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  real :: A( LDA, * )
  real :: C( LDC, * )
 end subroutine SSYRK
end interface

interface
 SUBROUTINE STBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: K
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  real :: X( * )
  real :: A( LDA, * )
 end subroutine STBMV
end interface

interface
 SUBROUTINE STBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: K
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  real :: X( * )
  real :: A( LDA, * )
 end subroutine STBSV
end interface

interface
 SUBROUTINE STPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  real :: AP( * )
  real :: X( * )
 end subroutine STPMV
end interface

interface
 SUBROUTINE STPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  real :: AP( * )
  real :: X( * )
 end subroutine STPSV
end interface

interface
 SUBROUTINE STRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  real :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  real :: A( LDA, * )
  real :: B( LDB, * )
 end subroutine STRMM
end interface

interface
 SUBROUTINE STRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  real :: X( * )
  real :: A( LDA, * )
 end subroutine STRMV
end interface

interface
 SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  real :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  real :: A( LDA, * )
  real :: B( LDB, * )
 end subroutine STRSM
end interface

interface
 SUBROUTINE STRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  real :: X( * )
  real :: A( LDA, * )
 end subroutine STRSV
end interface

interface
 complex(kind=8) function zdotu(n,zx,incx,zy,incy)
 implicit none
 integer :: incx
 integer :: incy
 integer :: n
 complex(kind=8) :: zx(*)
 complex(kind=8) :: zy(*)
end function zdotu
end interface

interface
 SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  double precision :: C
  double precision :: S
  complex(kind=8) :: CX( * )
  complex(kind=8) :: CY( * )
 end subroutine ZDROT
end interface

interface
 SUBROUTINE ZGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: KL
  integer :: KU
  integer :: LDA
  integer :: M
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: TRANS
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGBMV
end interface

interface
 SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: M
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGERU
end interface

interface
 SUBROUTINE ZHBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: K
  integer :: LDA
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZHBMV
end interface

interface
 SUBROUTINE ZHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
  BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZHEMV
end interface

interface
 SUBROUTINE ZHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: LDA
  integer :: N
  complex(kind=8) :: ALPHA
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZHER2
end interface

interface
 SUBROUTINE ZHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,&  
  C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: N
  complex(kind=8) :: ALPHA
  double precision :: BETA
  character :: TRANS
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZHER2K
end interface

interface
 SUBROUTINE ZHPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: UPLO
  complex(kind=8) :: AP( * )
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
 end subroutine ZHPMV
end interface

interface
 SUBROUTINE ZHPR  ( UPLO, N, ALPHA, X, INCX, AP )
  implicit none
  integer :: INCX
  integer :: N
  double precision :: ALPHA
  character*1 :: UPLO
  complex(kind=8) :: AP( * )
  complex(kind=8) :: X( * )
 end subroutine ZHPR
end interface

interface
 SUBROUTINE ZHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  complex(kind=8) :: ALPHA
  character*1 :: UPLO
  complex(kind=8) :: AP( * )
  complex(kind=8) :: X( * )
  complex(kind=8) :: Y( * )
 end subroutine ZHPR2
end interface

interface
 subroutine zrotg(ca,cb,c,s)
  implicit none
  double precision :: c
  complex(kind=8) :: ca
  complex(kind=8) :: cb
  complex(kind=8) :: s
 end subroutine zrotg
end interface

interface
 subroutine  zswap (n,zx,incx,zy,incy)
  implicit none
  integer :: incx
  integer :: incy
  integer :: n
  complex(kind=8) :: zx(*)
  complex(kind=8) :: zy(*)
 end subroutine zswap
end interface

interface
 SUBROUTINE ZSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: SIDE
  character*1 :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZSYMM
end interface

interface
 SUBROUTINE ZSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZSYR2K
end interface

interface
 SUBROUTINE ZSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
  BETA, C, LDC )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZSYRK
end interface

interface
 SUBROUTINE ZTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: K
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZTBMV
end interface

interface
 SUBROUTINE ZTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: K
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZTBSV
end interface

interface
 SUBROUTINE ZTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: AP( * )
  complex(kind=8) :: X( * )
 end subroutine ZTPMV
end interface

interface
 SUBROUTINE ZTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: AP( * )
  complex(kind=8) :: X( * )
 end subroutine ZTPSV
end interface

interface
 SUBROUTINE ZTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
  B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  complex(kind=8) :: ALPHA
  character*1 :: DIAG
  character*1 :: SIDE
  character*1 :: TRANSA
  character*1 :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
 end subroutine ZTRMM
end interface

interface
 SUBROUTINE ZTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZTRMV
end interface

interface
 SUBROUTINE ZTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
  implicit none
  integer :: INCX
  integer :: LDA
  integer :: N
  character*1 :: DIAG
  character*1 :: TRANS
  character*1 :: UPLO
  complex(kind=8) :: X( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZTRSV
end interface

interface
 SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  integer :: IPIV( * )
  complex :: A( LDA, * )
 end subroutine CGETF2
end interface

interface
 SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  integer :: IPIV( * )
  complex :: A( LDA, * )
 end subroutine CGETRF
end interface

interface
 SUBROUTINE CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  integer :: IPIV( * )
  complex :: WORK( * )
  complex :: A( LDA, * )
 end subroutine CGETRI
end interface

interface
 SUBROUTINE CHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,&  
  INFO )
  implicit none
  integer :: INFO
  integer :: LDZ
  integer :: N
  character :: JOBZ
  character :: UPLO
  complex :: AP( * )
  real :: RWORK( * )
  real :: W( * )
  complex :: WORK( * )
  complex :: Z( LDZ, * )
 end subroutine CHPEV
end interface

interface
 SUBROUTINE CHPTRD( UPLO, N, AP, D, E, TAU, INFO )
  implicit none
  integer :: INFO
  integer :: N
  character :: UPLO
  complex :: AP( * )
  real :: D( * )
  real :: E( * )
  complex :: TAU( * )
 end subroutine CHPTRD
end interface

interface
 COMPLEX FUNCTION CLADIV( X, Y )
 implicit none
 complex :: X
 complex :: Y
end function CLADIV
end interface

interface
 REAL             FUNCTION CLANHP( NORM, UPLO, N, AP, WORK )
 implicit none
 integer :: N
 character :: NORM
 character :: UPLO
 complex :: AP( * )
 real :: WORK( * )
end function CLANHP
end interface

interface
 SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
  implicit none
  integer :: INCV
  integer :: LDC
  integer :: M
  integer :: N
  character :: SIDE
  complex :: TAU
  complex :: V( * )
  complex :: WORK( * )
  complex :: C( LDC, * )
 end subroutine CLARF
end interface

interface
 SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )
  implicit none
  integer :: INCX
  integer :: N
  complex :: ALPHA
  complex :: TAU
  complex :: X( * )
 end subroutine CLARFG
end interface

interface
 SUBROUTINE CLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
  implicit none
  integer :: LDA
  integer :: M
  integer :: N
  character :: DIRECT
  character :: PIVOT
  character :: SIDE
  real :: C( * )
  real :: S( * )
  complex :: A( LDA, * )
 end subroutine CLASR
end interface

interface
 SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )
  implicit none
  integer :: INCX
  integer :: N
  real :: SCALE
  real :: SUMSQ
  complex :: X( * )
 end subroutine CLASSQ
end interface

interface
 SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX )
  implicit none
  integer :: INCX
  integer :: K1
  integer :: K2
  integer :: LDA
  integer :: N
  integer :: IPIV( * )
  complex :: A( LDA, * )
 end subroutine CLASWP
end interface

interface
 SUBROUTINE CLAZRO( M, N, ALPHA, BETA, A, LDA )
  implicit none
  integer :: LDA
  integer :: M
  integer :: N
  complex :: ALPHA
  complex :: BETA
  complex :: A( LDA, * )
 end subroutine CLAZRO
end interface

interface
 SUBROUTINE CSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDZ
  integer :: N
  character :: COMPZ
  real :: D( * )
  real :: E( * )
  real :: WORK( * )
  complex :: Z( LDZ, * )
 end subroutine CSTEQR
end interface

interface
 SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: DIAG
  character :: UPLO
  complex :: A( LDA, * )
 end subroutine CTRTRI
end interface

interface
 SUBROUTINE CUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: M
  integer :: N
  complex :: TAU( * )
  complex :: WORK( * )
  complex :: A( LDA, * )
 end subroutine CUNG2L
end interface

interface
 SUBROUTINE CUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: M
  integer :: N
  complex :: TAU( * )
  complex :: WORK( * )
  complex :: A( LDA, * )
 end subroutine CUNG2R
end interface

interface
 SUBROUTINE CUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDQ
  integer :: N
  character :: UPLO
  complex :: AP( * )
  complex :: TAU( * )
  complex :: WORK( * )
  complex :: Q( LDQ, * )
 end subroutine CUPGTR
end interface

interface
 SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,&  
  LDU, C, LDC, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDC
  integer :: LDU
  integer :: LDVT
  integer :: N
  integer :: NCC
  integer :: NCVT
  integer :: NRU
  character :: UPLO
  double precision :: D( * )
  double precision :: E( * )
  double precision :: WORK( * )
  double precision :: C( LDC, * )
  double precision :: U( LDU, * )
  double precision :: VT( LDVT, * )
 end subroutine DBDSQR
end interface

interface
 SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: D( * )
  double precision :: E( * )
  double precision :: TAUP( * )
  double precision :: TAUQ( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGEBD2
end interface

interface
 SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,&  
  INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: D( * )
  double precision :: E( * )
  double precision :: TAUP( * )
  double precision :: TAUQ( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGEBRD
end interface

interface
 SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGELQ2
end interface

interface
 SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGELQF
end interface

interface
 SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,&  
  WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDB
  integer :: LWORK
  integer :: M
  integer :: N
  integer :: NRHS
  integer :: RANK
  double precision :: RCOND
  double precision :: S( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DGELSS
end interface

interface
 SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGEQR2
end interface

interface
 SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGEQRF
end interface

interface
 SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&  
  WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDU
  integer :: LDVT
  integer :: LWORK
  integer :: M
  integer :: N
  character :: JOBU
  character :: JOBVT
  double precision :: S( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: U( LDU, * )
  double precision :: VT( LDVT, * )
 end subroutine DGESVD
end interface

interface
 SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  integer :: IPIV( * )
  double precision :: A( LDA, * )
 end subroutine DGETF2
end interface

interface
 SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  integer :: IPIV( * )
  double precision :: A( LDA, * )
 end subroutine DGETRF
end interface

interface
 SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  integer :: IPIV( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DGETRI
end interface

interface
 SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDQ
  integer :: N
  character :: UPLO
  double precision :: AP( * )
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: Q( LDQ, * )
 end subroutine DOPGTR
end interface

interface
 SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORG2L
end interface

interface
 SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORG2R
end interface

interface
 SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  character :: VECT
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORGBR
end interface

interface
 SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORGL2
end interface

interface
 SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORGLQ
end interface

interface
 SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORGQL
end interface

interface
 SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DORGQR
end interface

interface
 SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  character :: UPLO
  double precision :: TAU( * )
  double precision :: A( LDA, * )
  double precision :: WORK( LWORK )
 end subroutine DORGTR
end interface

interface
 SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
  WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: M
  integer :: N
  character :: SIDE
  character :: TRANS
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: C( LDC, * )
 end subroutine DORM2R
end interface

interface
 SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,&  
  LDC, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: LWORK
  integer :: M
  integer :: N
  character :: SIDE
  character :: TRANS
  character :: VECT
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: C( LDC, * )
 end subroutine DORMBR
end interface

interface
 SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
  WORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: M
  integer :: N
  character :: SIDE
  character :: TRANS
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: C( LDC, * )
 end subroutine DORML2
end interface

interface
 SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
  WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: LWORK
  integer :: M
  integer :: N
  character :: SIDE
  character :: TRANS
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: C( LDC, * )
 end subroutine DORMLQ
end interface

interface
 SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
  WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LDC
  integer :: LWORK
  integer :: M
  integer :: N
  character :: SIDE
  character :: TRANS
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: C( LDC, * )
 end subroutine DORMQR
end interface

interface
 SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDB
  integer :: N
  integer :: NRHS
  character :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DPOSV
end interface

interface
 SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: UPLO
  double precision :: A( LDA, * )
 end subroutine DPOTF2
end interface

interface
 SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: UPLO
  double precision :: A( LDA, * )
 end subroutine DPOTRF
end interface

interface
 SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDB
  integer :: N
  integer :: NRHS
  character :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DPOTRS
end interface

interface
 SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
  implicit none
  integer :: INFO
  integer :: N
  character :: UPLO
  double precision :: AP( * )
 end subroutine DPPTRF
end interface

interface
 SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: N
  character :: UPLO
  double precision :: AP( * )
  double precision :: BP( * )
 end subroutine DSPGST
end interface

interface
 SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,&  
  INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: LDZ
  integer :: N
  character :: JOBZ
  character :: UPLO
  double precision :: AP( * )
  double precision :: BP( * )
  double precision :: W( * )
  double precision :: WORK( * )
  double precision :: Z( LDZ, * )
 end subroutine DSPGV
end interface

interface
 SUBROUTINE DRSCL( N, SA, SX, INCX )
  implicit none
  integer :: INCX
  integer :: N
  double precision :: SA
  double precision :: SX( * )
 end subroutine DRSCL
end interface

interface
 SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDZ
  integer :: N
  character :: JOBZ
  character :: UPLO
  double precision :: AP( * )
  double precision :: W( * )
  double precision :: WORK( * )
  double precision :: Z( LDZ, * )
 end subroutine DSPEV
end interface

interface
 SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )
  implicit none
  integer :: INFO
  integer :: N
  character :: UPLO
  double precision :: AP( * )
  double precision :: D( * )
  double precision :: E( * )
  double precision :: TAU( * )
 end subroutine DSPTRD
end interface

interface
 SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E,&  
  M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,&  
  INFO )
  implicit none
  integer :: IL
  integer :: INFO
  integer :: IU
  integer :: M
  integer :: N
  integer :: NSPLIT
  double precision :: ABSTOL
  character :: ORDER
  character :: RANGE
  double precision :: VL
  double precision :: VU
  double precision :: D( * )
  double precision :: E( * )
  integer :: IBLOCK( * )
  integer :: ISPLIT( * )
  integer :: IWORK( * )
  double precision :: W( * )
  double precision :: WORK( * )
 end subroutine DSTEBZ
end interface

interface
 SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDZ
  integer :: N
  character :: COMPZ
  double precision :: D( * )
  double precision :: E( * )
  double precision :: WORK( * )
  double precision :: Z( LDZ, * )
 end subroutine DSTEQR
end interface

interface
 SUBROUTINE DSTERF( N, D, E, INFO )
  implicit none
  integer :: INFO
  integer :: N
  double precision :: D( * )
  double precision :: E( * )
 end subroutine DSTERF
end interface

interface
 SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  character :: JOBZ
  character :: UPLO
  double precision :: W( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DSYEV
end interface

interface
 SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: LDA
  integer :: LDB
  integer :: N
  character :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DSYGS2
end interface

interface
 SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: LDA
  integer :: LDB
  integer :: N
  character :: UPLO
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DSYGST
end interface

interface
 SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,&  
  LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: LDA
  integer :: LDB
  integer :: LWORK
  integer :: N
  character :: JOBZ
  character :: UPLO
  double precision :: W( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DSYGV
end interface

interface
 SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,&  
  LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDB
  integer :: LWORK
  integer :: N
  integer :: NRHS
  character :: UPLO
  integer :: IPIV( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DSYSV
end interface

interface
 SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: UPLO
  double precision :: D( * )
  double precision :: E( * )
  double precision :: TAU( * )
  double precision :: A( LDA, * )
 end subroutine DSYTD2
end interface

interface
 SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: UPLO
  integer :: IPIV( * )
  double precision :: A( LDA, * )
 end subroutine DSYTF2
end interface

interface
 SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  character :: UPLO
  double precision :: D( * )
  double precision :: E( * )
  double precision :: TAU( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DSYTRD
end interface

interface
 SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  character :: UPLO
  integer :: IPIV( * )
  double precision :: WORK( * )
  double precision :: A( LDA, * )
 end subroutine DSYTRF
end interface

interface
 SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDB
  integer :: N
  integer :: NRHS
  character :: UPLO
  integer :: IPIV( * )
  double precision :: A( LDA, * )
  double precision :: B( LDB, * )
 end subroutine DSYTRS
end interface

interface
 SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: DIAG
  character :: UPLO
  double precision :: A( LDA, * )
 end subroutine DTRTI2
end interface

interface
 SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: DIAG
  character :: UPLO
  double precision :: A( LDA, * )
 end subroutine DTRTRI
end interface

interface
 DOUBLE PRECISION FUNCTION DZSUM1( N, CX, INCX )
 implicit none
 integer :: INCX
 integer :: N
 complex(kind=8) :: CX( * )
end function DZSUM1
end interface

interface
 INTEGER FUNCTION IZMAX1( N, CX, INCX )
 implicit none
 integer :: INCX
 integer :: N
 complex(kind=8) :: CX( * )
end function IZMAX1
end interface

interface
 SUBROUTINE SSTERF( N, D, E, INFO )
  implicit none
  integer :: INFO
  integer :: N
  real :: D( * )
  real :: E( * )
 end subroutine SSTERF
end interface

interface
 SUBROUTINE ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,&  
  LDU, C, LDC, RWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDC
  integer :: LDU
  integer :: LDVT
  integer :: N
  integer :: NCC
  integer :: NCVT
  integer :: NRU
  character :: UPLO
  double precision :: D( * )
  double precision :: E( * )
  double precision :: RWORK( * )
  complex(kind=8) :: C( LDC, * )
  complex(kind=8) :: U( LDU, * )
  complex(kind=8) :: VT( LDVT, * )
 end subroutine ZBDSQR
end interface

interface
 SUBROUTINE ZGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,&  
  INFO )
  implicit none
  integer :: IHI
  integer :: ILO
  integer :: INFO
  integer :: LDV
  integer :: M
  integer :: N
  character :: JOB
  character :: SIDE
  double precision :: SCALE( * )
  complex(kind=8) :: V( LDV, * )
 end subroutine ZGEBAK
end interface

interface
 SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
  implicit none
  integer :: IHI
  integer :: ILO
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: JOB
  double precision :: SCALE( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGEBAL
end interface

interface
 SUBROUTINE ZGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: D( * )
  double precision :: E( * )
  complex(kind=8) :: TAUP( * )
  complex(kind=8) :: TAUQ( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGEBD2
end interface

interface
 SUBROUTINE ZGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,&  
  INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  double precision :: D( * )
  double precision :: E( * )
  complex(kind=8) :: TAUP( * )
  complex(kind=8) :: TAUQ( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGEBRD
end interface

interface
 SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,&  
  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDVS
  integer :: LWORK
  integer :: N
  integer :: SDIM
  character :: JOBVS
  logical :: SELECT
  character :: SORT
  logical :: BWORK( * )
  double precision :: RWORK( * )
  complex(kind=8) :: W( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: VS( LDVS, * )
 end subroutine ZGEES
end interface

interface
 SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,&  
  WORK, LWORK, RWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDVL
  integer :: LDVR
  integer :: LWORK
  integer :: N
  character :: JOBVL
  character :: JOBVR
  double precision :: RWORK( * )
  complex(kind=8) :: W( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: VL( LDVL, * )
  complex(kind=8) :: VR( LDVR, * )
 end subroutine ZGEEV
end interface

interface
 SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: IHI
  integer :: ILO
  integer :: INFO
  integer :: LDA
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGEHD2
end interface

interface
 SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: IHI
  integer :: ILO
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: WORK( LWORK )
 end subroutine ZGEHRD
end interface

interface
 SUBROUTINE ZGELQ2( M, N, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGELQ2
end interface

interface
 SUBROUTINE ZGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGELQF
end interface

interface
 SUBROUTINE ZGEQR2( M, N, A, LDA, TAU, WORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGEQR2
end interface

interface
 SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGEQRF
end interface

interface
 SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&  
  WORK, LWORK, RWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDU
  integer :: LDVT
  integer :: LWORK
  integer :: M
  integer :: N
  character :: JOBU
  character :: JOBVT
  double precision :: RWORK( * )
  double precision :: S( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: U( LDU, * )
  complex(kind=8) :: VT( LDVT, * )
 end subroutine ZGESVD
end interface

interface
 SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: M
  integer :: N
  integer :: IPIV( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGETF2
end interface

interface
 SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  integer :: IPIV( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZGETRI
end interface

interface
 SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LDB
  integer :: N
  integer :: NRHS
  character :: TRANS
  integer :: IPIV( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
 end subroutine ZGETRS
end interface

interface
 SUBROUTINE ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: LDA
  integer :: LDB
  integer :: N
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
 end subroutine ZHEGS2
end interface

interface
 SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: UPLO
  double precision :: D( * )
  double precision :: E( * )
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZHETD2
end interface

interface
 SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: LWORK
  integer :: N
  character :: UPLO
  double precision :: D( * )
  double precision :: E( * )
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZHETRD
end interface

interface
 SUBROUTINE ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
  implicit none
  integer :: INFO
  integer :: ITYPE
  integer :: N
  character :: UPLO
  complex(kind=8) :: AP( * )
  complex(kind=8) :: BP( * )
 end subroutine ZHPGST
end interface

interface
 SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )
  implicit none
  integer :: INFO
  integer :: N
  character :: UPLO
  complex(kind=8) :: AP( * )
  double precision :: D( * )
  double precision :: E( * )
  complex(kind=8) :: TAU( * )
 end subroutine ZHPTRD
end interface

interface
 SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,&  
  WORK, LWORK, INFO )
  implicit none
  integer :: IHI
  integer :: ILO
  integer :: INFO
  integer :: LDH
  integer :: LDZ
  integer :: LWORK
  integer :: N
  character :: COMPZ
  character :: JOB
  complex(kind=8) :: W( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZHSEQR
end interface

interface
 SUBROUTINE ZLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,&  
  LDY )
  implicit none
  integer :: LDA
  integer :: LDX
  integer :: LDY
  integer :: M
  integer :: N
  integer :: NB
  double precision :: D( * )
  double precision :: E( * )
  complex(kind=8) :: TAUP( * )
  complex(kind=8) :: TAUQ( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: X( LDX, * )
  complex(kind=8) :: Y( LDY, * )
 end subroutine ZLABRD
end interface

interface
 SUBROUTINE ZLACGV( N, X, INCX )
  implicit none
  integer :: INCX
  integer :: N
  complex(kind=8) :: X( * )
 end subroutine ZLACGV
end interface

interface
 SUBROUTINE ZLACON( N, V, X, EST, KASE )
  implicit none
  integer :: KASE
  integer :: N
  double precision :: EST
  complex(kind=8) :: V( N )
  complex(kind=8) :: X( N )
 end subroutine ZLACON
end interface

interface
 SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
  implicit none
  integer :: LDA
  integer :: LDB
  integer :: M
  integer :: N
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
 end subroutine ZLACPY
end interface

interface
 complex(kind=8)   FUNCTION ZLADIV( X, Y )
 implicit none
 complex(kind=8) :: X
 complex(kind=8) :: Y
end function ZLADIV
end interface

interface
 SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,&  
  IHIZ, Z, LDZ, INFO )
  implicit none
  integer :: IHI
  integer :: IHIZ
  integer :: ILO
  integer :: ILOZ
  integer :: INFO
  integer :: LDH
  integer :: LDZ
  integer :: N
  logical :: WANTT
  logical :: WANTZ
  complex(kind=8) :: W( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZLAHQR
end interface

interface
 SUBROUTINE ZLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDT
  integer :: LDY
  integer :: N
  integer :: NB
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: T( LDT, NB )
  complex(kind=8) :: TAU( NB )
  complex(kind=8) :: Y( LDY, NB )
 end subroutine ZLAHRD
end interface

interface
 SUBROUTINE ZLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
  implicit none
  integer :: K
  integer :: LDA
  integer :: LDT
  integer :: LDY
  integer :: N
  integer :: NB
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: T( LDT, NB )
  complex(kind=8) :: TAU( NB )
  complex(kind=8) :: Y( LDY, NB )
 end subroutine ZLAHR2
end interface

interface
 DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
 implicit none
 integer :: LDA
 integer :: M
 integer :: N
 character :: NORM
 double precision :: WORK( * )
 complex(kind=8) :: A( LDA, * )
end function ZLANGE
end interface

interface
 DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
 implicit none
 integer :: LDA
 integer :: N
 character :: NORM
 character :: UPLO
 double precision :: WORK( * )
 complex(kind=8) :: A( LDA, * )
end function ZLANHE
end interface

interface
 DOUBLE PRECISION FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK )
 implicit none
 integer :: N
 character :: NORM
 character :: UPLO
 complex(kind=8) :: AP( * )
 double precision :: WORK( * )
end function ZLANHP
end interface

interface
 DOUBLE PRECISION FUNCTION ZLANHS( NORM, N, A, LDA, WORK )
 implicit none
 integer :: LDA
 integer :: N
 character :: NORM
 double precision :: WORK( * )
 complex(kind=8) :: A( LDA, * )
end function ZLANHS
end interface

interface
 SUBROUTINE ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,&  
  IHIZ, Z, LDZ, WORK, LWORK, INFO )
  implicit none
  integer :: IHI
  integer :: IHIZ
  integer :: ILO
  integer :: ILOZ
  integer :: INFO
  integer :: LDH
  integer :: LDZ
  integer :: LWORK
  integer :: N
  logical :: WANTT
  logical :: WANTZ
  complex(kind=8) :: W( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZLAQR0
end interface

interface
 SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )
  implicit none
  integer :: LDH
  integer :: N
  complex(kind=8) :: S1
  complex(kind=8) :: S2
  complex(kind=8) :: V( * )
  complex(kind=8) :: H( LDH, * )
 end subroutine ZLAQR1
end interface

interface
 SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,&  
  IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,&  
  NV, WV, LDWV, WORK, LWORK )
  implicit none
  integer :: IHIZ
  integer :: ILOZ
  integer :: KBOT
  integer :: KTOP
  integer :: LDH
  integer :: LDT
  integer :: LDV
  integer :: LDWV
  integer :: LDZ
  integer :: LWORK
  integer :: N
  integer :: ND
  integer :: NH
  integer :: NS
  integer :: NV
  integer :: NW
  logical :: WANTT
  logical :: WANTZ
  complex(kind=8) :: SH( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: T( LDT, * )
  complex(kind=8) :: V( LDV, * )
  complex(kind=8) :: WV( LDWV, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZLAQR2
end interface

interface
 SUBROUTINE ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,&  
  IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,&  
  NV, WV, LDWV, WORK, LWORK )
  implicit none
  integer :: IHIZ
  integer :: ILOZ
  integer :: KBOT
  integer :: KTOP
  integer :: LDH
  integer :: LDT
  integer :: LDV
  integer :: LDWV
  integer :: LDZ
  integer :: LWORK
  integer :: N
  integer :: ND
  integer :: NH
  integer :: NS
  integer :: NV
  integer :: NW
  logical :: WANTT
  logical :: WANTZ
  complex(kind=8) :: SH( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: T( LDT, * )
  complex(kind=8) :: V( LDV, * )
  complex(kind=8) :: WV( LDWV, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZLAQR3
end interface

interface
 SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,&  
  IHIZ, Z, LDZ, WORK, LWORK, INFO )
  implicit none
  integer :: IHI
  integer :: IHIZ
  integer :: ILO
  integer :: ILOZ
  integer :: INFO
  integer :: LDH
  integer :: LDZ
  integer :: LWORK
  integer :: N
  logical :: WANTT
  logical :: WANTZ
  complex(kind=8) :: W( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZLAQR4
end interface

interface
 SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,&  
  H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,&  
  WV, LDWV, NH, WH, LDWH )
  implicit none
  integer :: IHIZ
  integer :: ILOZ
  integer :: KACC22
  integer :: KBOT
  integer :: KTOP
  integer :: LDH
  integer :: LDU
  integer :: LDV
  integer :: LDWH
  integer :: LDWV
  integer :: LDZ
  integer :: N
  integer :: NH
  integer :: NSHFTS
  integer :: NV
  logical :: WANTT
  logical :: WANTZ
  complex(kind=8) :: S( * )
  complex(kind=8) :: H( LDH, * )
  complex(kind=8) :: U( LDU, * )
  complex(kind=8) :: V( LDV, * )
  complex(kind=8) :: WH( LDWH, * )
  complex(kind=8) :: WV( LDWV, * )
  complex(kind=8) :: Z( LDZ, * )
 end subroutine ZLAQR5
end interface

interface
 SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
  implicit none
  integer :: INCV
  integer :: LDC
  integer :: M
  integer :: N
  character :: SIDE
  complex(kind=8) :: TAU
  complex(kind=8) :: V( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZLARF
end interface

interface
 SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,&  
  T, LDT, C, LDC, WORK, LDWORK )
  implicit none
  integer :: K
  integer :: LDC
  integer :: LDT
  integer :: LDV
  integer :: LDWORK
  integer :: M
  integer :: N
  character :: DIRECT
  character :: SIDE
  character :: STOREV
  character :: TRANS
  complex(kind=8) :: C( LDC, * )
  complex(kind=8) :: T( LDT, * )
  complex(kind=8) :: V( LDV, * )
  complex(kind=8) :: WORK( LDWORK, * )
 end subroutine ZLARFB
end interface

interface
 SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
  implicit none
  integer :: INCX
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: TAU
  complex(kind=8) :: X( * )
 end subroutine ZLARFG
end interface

interface
 SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
  implicit none
  integer :: K
  integer :: LDT
  integer :: LDV
  integer :: N
  character :: DIRECT
  character :: STOREV
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: T( LDT, * )
  complex(kind=8) :: V( LDV, * )
 end subroutine ZLARFT
end interface

interface
 SUBROUTINE ZLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
  implicit none
  integer :: LDC
  integer :: M
  integer :: N
  character :: SIDE
  complex(kind=8) :: TAU
  complex(kind=8) :: V( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZLARFX
end interface

interface
 SUBROUTINE ZLARTG( F, G, CS, SN, R )
  implicit none
  double precision :: CS
  complex(kind=8) :: F
  complex(kind=8) :: G
  complex(kind=8) :: R
  complex(kind=8) :: SN
 end subroutine ZLARTG
end interface

interface
 SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: KL
  integer :: KU
  integer :: LDA
  integer :: M
  integer :: N
  double precision :: CFROM
  double precision :: CTO
  character :: TYPE
  complex(kind=8) :: A( LDA, * )
 end subroutine ZLASCL
end interface

interface
 SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
  implicit none
  integer :: LDA
  integer :: M
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
 end subroutine ZLASET
end interface

interface
 SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
  implicit none
  integer :: LDA
  integer :: M
  integer :: N
  character :: DIRECT
  character :: PIVOT
  character :: SIDE
  double precision :: C( * )
  double precision :: S( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZLASR
end interface

interface
 SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
  implicit none
  integer :: INCX
  integer :: N
  double precision :: SCALE
  double precision :: SUMSQ
  complex(kind=8) :: X( * )
 end subroutine ZLASSQ
end interface

interface
 SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
  implicit none
  integer :: INCX
  integer :: K1
  integer :: K2
  integer :: LDA
  integer :: N
  integer :: IPIV( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZLASWP
end interface

interface
 SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
  implicit none
  integer :: LDA
  integer :: LDW
  integer :: N
  integer :: NB
  character :: UPLO
  double precision :: E( * )
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: W( LDW, * )
 end subroutine ZLATRD
end interface

interface
 SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,&  
  CNORM, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: DIAG
  character :: NORMIN
  double precision :: SCALE
  character :: TRANS
  character :: UPLO
  double precision :: CNORM( * )
  complex(kind=8) :: X( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZLATRS
end interface

interface
 SUBROUTINE ZLAZRO( M, N, ALPHA, BETA, A, LDA )
  implicit none
  integer :: LDA
  integer :: M
  integer :: N
  complex(kind=8) :: ALPHA
  complex(kind=8) :: BETA
  complex(kind=8) :: A( LDA, * )
 end subroutine ZLAZRO
end interface

interface
 SUBROUTINE ZPOTF2( UPLO, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
 end subroutine ZPOTF2
end interface

interface
 SUBROUTINE ZPPTRF( UPLO, N, AP, INFO )
  implicit none
  integer :: INFO
  integer :: N
  character :: UPLO
  complex(kind=8) :: AP( * )
 end subroutine ZPPTRF
end interface

interface
 SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
  implicit none
  integer :: INCX
  integer :: INCY
  integer :: N
  double precision :: C
  complex(kind=8) :: S
  complex(kind=8) :: CX( * )
  complex(kind=8) :: CY( * )
 end subroutine ZROT
end interface

interface
 SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,&  
  LDVR, MM, M, WORK, RWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDT
  integer :: LDVL
  integer :: LDVR
  integer :: M
  integer :: MM
  integer :: N
  character :: HOWMNY
  character :: SIDE
  double precision :: RWORK( * )
  logical :: SELECT( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: T( LDT, * )
  complex(kind=8) :: VL( LDVL, * )
  complex(kind=8) :: VR( LDVR, * )
 end subroutine ZTREVC
end interface

interface
 SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
  implicit none
  integer :: IFST
  integer :: ILST
  integer :: INFO
  integer :: LDQ
  integer :: LDT
  integer :: N
  character :: COMPQ
  complex(kind=8) :: Q( LDQ, * )
  complex(kind=8) :: T( LDT, * )
 end subroutine ZTREXC
end interface

interface
 SUBROUTINE ZTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,&  
  SEP, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: LDQ
  integer :: LDT
  integer :: LWORK
  integer :: M
  integer :: N
  character :: COMPQ
  character :: JOB
  double precision :: S
  double precision :: SEP
  logical :: SELECT( * )
  complex(kind=8) :: W( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: Q( LDQ, * )
  complex(kind=8) :: T( LDT, * )
 end subroutine ZTRSEN
end interface

interface
 SUBROUTINE ZTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,&  
  LDC, SCALE, INFO )
  implicit none
  integer :: INFO
  integer :: ISGN
  integer :: LDA
  integer :: LDB
  integer :: LDC
  integer :: M
  integer :: N
  double precision :: SCALE
  character :: TRANA
  character :: TRANB
  complex(kind=8) :: A( LDA, * )
  complex(kind=8) :: B( LDB, * )
  complex(kind=8) :: C( LDC, * )
 end subroutine ZTRSYL
end interface

interface
 SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: DIAG
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
 end subroutine ZTRTI2
end interface

interface
 SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
  implicit none
  integer :: INFO
  integer :: LDA
  integer :: N
  character :: DIAG
  character :: UPLO
  complex(kind=8) :: A( LDA, * )
 end subroutine ZTRTRI
end interface

interface
 SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
  implicit none
  integer :: INFO
  integer :: K
  integer :: LDA
  integer :: LWORK
  integer :: M
  integer :: N
  complex(kind=8) :: TAU( * )
  complex(kind=8) :: WORK( * )
  complex(kind=8) :: A( LDA, * )
 end subroutine ZUNGQR
end interface

end module interfaces_linalg
!!***
