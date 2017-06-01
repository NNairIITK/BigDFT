!	 Copyright (c) 2012 The Regents of the University of California,
!	 through Lawrence Berkeley National Laboratory.  
!
!  Author: Lin Lin
!	 
!  This file is part of PEXSI. All rights reserved.
!
!	 Redistribution and use in source and binary forms, with or without
!	 modification, are permitted provided that the following conditions are met:
!
!	 (1) Redistributions of source code must retain the above copyright notice, this
!	 list of conditions and the following disclaimer.
!	 (2) Redistributions in binary form must reproduce the above copyright notice,
!	 this list of conditions and the following disclaimer in the documentation
!	 and/or other materials provided with the distribution.
!	 (3) Neither the name of the University of California, Lawrence Berkeley
!	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
!	 be used to endorse or promote products derived from this software without
!	 specific prior written permission.
!
!	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
!	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!	 You are under no obligation whatsoever to provide any bug fixes, patches, or
!	 upgrades to the features, functionality or performance of the source code
!	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
!	 available either publicly, or directly to Lawrence Berkeley National
!	 Laboratory, without imposing a separate written license agreement for such
!	 Enhancements, then you hereby grant the following license: a non-exclusive,
!	 royalty-free perpetual license to install, use, modify, prepare derivative
!	 works, incorporate into other computer software, distribute, and sublicense
!	 such enhancements or derivative works thereof, in binary and source code form.
!

module pexsi_base
  use, intrinsic :: iso_c_binding
  implicit none

  private

  public :: f_ppexsi_options


  ! Struct for PPEXSIOptions
  ! NOTE: The order and the type of the parameters must be strictly the same as in PPEXSIOptions in 
  ! c_pexsi_interface.h
  type, bind(C) :: f_ppexsi_options
    real(c_double)         :: temperature
    real(c_double)         :: gap
    real(c_double)         :: deltaE
    integer(c_int)         :: numPole
    integer(c_int)         :: isInertiaCount
    integer(c_int)         :: maxPEXSIIter
    real(c_double)         :: muMin0
    real(c_double)         :: muMax0
    real(c_double)         :: mu0
    real(c_double)         :: muInertiaTolerance
    real(c_double)         :: muInertiaExpansion
    real(c_double)         :: muPEXSISafeGuard
    real(c_double)         :: numElectronPEXSITolerance
    integer(c_int)         :: matrixType
    integer(c_int)         :: isSymbolicFactorize
    integer(c_int)         :: isConstructCommPattern
    integer(c_int)         :: ordering
    integer(c_int)         :: npSymbFact
    integer(c_int)         :: symmetric
    integer(c_int)         :: verbosity
  end type f_ppexsi_options

end module pexsi_base
