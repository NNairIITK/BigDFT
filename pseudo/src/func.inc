C hello emacs, this is -*- fortran -*-
C     ==================================================================
      INTEGER       MFXCX,MFXCC,MGCX,MGCC,MHFX,MTAU,MSRX,MSELF
      COMMON/FUNC1/ MFXCX,MFXCC,MGCX,MGCC,MHFX,MTAU,MSRX,MSELF
C     ==================================================================
      REAL*8        SALPHA,BBETA,BETAPP,SRXA
      COMMON/FUNC2/ SALPHA,BBETA,BETAPP,SRXA
C     ==================================================================
      REAL*8        PXLDA,PXGC,PCLDA,PCGC,PHFX
      COMMON/FUNC3/ PXLDA,PXGC,PCLDA,PCGC,PHFX
C     ==================================================================
      REAL*8        F0HF,F0HFM,ALHFX
      COMMON/FUNC4/ F0HF,F0HFM,ALHFX
C     ==================================================================
      INTEGER       NLLI,NLLL,NLLR,MXCSET
      PARAMETER     (MXCSET=5)
      PARAMETER     (NLLI=6)
      PARAMETER     (NLLL=6)
      PARAMETER     (NLLR=5)
      INTEGER       IXCPAR(NLLI,MXCSET)
      LOGICAL       LXCPAR(NLLL,MXCSET)
      REAL*8        RXCPAR(NLLR,MXCSET)
      COMMON/FUXCI/ IXCPAR
      COMMON/FUXCL/ LXCPAR
      COMMON/FUXCR/ RXCPAR
      INTEGER       IFXC0,IFXCK,IFXCR
      PARAMETER     (IFXC0=1)
      PARAMETER     (IFXCK=2)
      PARAMETER     (IFXCR=3)
C     ------------------------------------------------------------------
C     IXCPAR(1,*) = MFXCX
C     IXCPAR(2,*) = MFXCC
C     IXCPAR(3,*) = MGCX
C     IXCPAR(4,*) = MGCC
C     IXCPAR(5,*) = MHFX
C     IXCPAR(6,*) = MTAU
C     ------------------------------------------------------------------
C     LXCPAR(1,*) = TGC
C     LXCPAR(2,*) = TGCX
C     LXCPAR(3,*) = TGCC
C     LXCPAR(4,*) = TOLDCODE
C     LXCPAR(5,*) = THYBRID
C     LXCPAR(6,*) = TTAU
C     ------------------------------------------------------------------
C     RXCPAR(1,*) = PXLDA
C     RXCPAR(2,*) = PCLDA
C     RXCPAR(3,*) = PXGC
C     RXCPAR(4,*) = PCGC
C     RXCPAR(5,*) = PHFX
C     ==================================================================
