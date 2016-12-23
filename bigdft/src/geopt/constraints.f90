!........ SHAKE AND RATTLE IN bigDFT  -----------

module constraints
      implicit none
      type, public :: cons_data
          integer                    ::  nconstraints              
          integer                    ::  algorithm                
          integer                    ::  diis_size     
          REAL(KIND=8), DIMENSION(:,:), POINTER :: CONSTRAINTS_TYPE(:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: COORDINATE_ATOMLIST(:,:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: TARG(:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: cc(:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: ffcv(:,:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: sigma(:,:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: sigma_gd(:,:)
          REAL(KIND=8), DIMENSION(:,:), POINTER :: dt2bym(:)
          real(kind=8)               ::  const_conv                
          integer                    ::  fflage
          
      end type 
      contains


      subroutine read_constraints(ndim,natoms,cvdata)
        use module_base    !to use constants
        implicit none
        integer , intent(in) :: ndim
        integer , intent(in) :: natoms
        type(cons_data), intent(inout) :: cvdata 
        !
        integer :: i,ii
        integer :: I0 
        integer, parameter  :: lmax = 128
        character(len=lmax) :: line           
        !character(len=lmax) :: mthd

        !
      cvdata%nconstraints=0
      open(133,file='constraints.in',status='old',iostat=i0) ! we will read this file from md_run.f90 
      if (i0/=0) return 
 
      read (133,*) cvdata%nconstraints, cvdata%const_conv, cvdata%algorithm, cvdata%diis_size
            !print*,cvdata%nconstraints, cvdata%const_conv, cvdata%algorithm, cvdata%diis_size
      if (cvdata%algorithm==1) then 
             print*, "Algorithm 1 selected: Lapack LU Decomposition for solving Linear System (for shake and rattle_) "

      else if (cvdata%algorithm==2) then
             print *, "Algorithm 2 selected: Conjugate Gradient for solving Linear System (for shake and rattle) "
      else if (cvdata%algorithm==3) then
             print *, "Algorithm 3 selected: DIIS for solving Linear System(for shake and rattle)  "
      end if 

        cvdata%CONSTRAINTS_TYPE  =  f_malloc_ptr([cvdata%nconstraints],id='cvdata%CONSTRAINTS_TYPE')
        cvdata%COORDINATE_ATOMLIST = f_malloc_ptr([cvdata%nconstraints,4],id='cvdata%COORDINATE_ATOMLIST')
        cvdata%TARG = f_malloc_ptr([cvdata%nconstraints],id='cvdata%TARG')
       cvdata%cc = f_malloc_ptr([cvdata%nconstraints],id='cvdata%cc')
        cvdata%ffcv = f_malloc_ptr([ndim*natoms,cvdata%nconstraints],id='cvdata%ffcv')
        cvdata%sigma = f_malloc_ptr([cvdata%nconstraints,1],id='cvdata%sigma')
        cvdata%sigma_gd = f_malloc_ptr([ndim*natoms,cvdata%nconstraints],id='cvdata%sigma_gd')
        cvdata%dt2bym =f_malloc_ptr([ndim*natoms],id='cvdata%dt2bym')

      do i=1,cvdata%nconstraints
       read (133,'(128A)') line(1:lmax)
        if (index(line,"DIST")/=0) THEN
          cvdata%CONSTRAINTS_TYPE(i)=1
          I0=index(line,"DIST")
          I0=I0+4
              READ(line(I0:lmax),*) cvdata%COORDINATE_ATOMLIST(i,1:2),cvdata%TARG(i)  
             cvdata%TARG(i)= cvdata%TARG(i)/Bohr_ang
        ELSE if (index(line,"ANGLE")/=0)THEN
              cvdata%CONSTRAINTS_TYPE(i)=2
              I0=index(line,"ANGLE")
              READ(LINE(I0+5:128),*) cvdata%COORDINATE_ATOMLIST(i,1:3),cvdata%TARG(i)
              cvdata%TARG(i) = (cvdata%TARG(i)*pi)/180.d0
        ELSE if (index(line,"DIHEDRAL")/=0)THEN
              cvdata%CONSTRAINTS_TYPE(i)=3
              I0=index(line,"DIHEDRAL")
              READ(LINE(I0+8:128),*) cvdata%COORDINATE_ATOMLIST(i,1:4),cvdata%TARG(i)
              cvdata%TARG(i) = (cvdata%TARG(i)*pi)/180.d0
        ELSE
              STOP "UNKONWN CONSTR."
        end if 
      end do
      close (133)
      end subroutine read_constraints 

! ---------

      subroutine  deallocate_cvdata(cvdata)
      use module_base
      implicit none
      type(cons_data), intent(inout) :: cvdata

      call f_free_ptr(cvdata%CONSTRAINTS_TYPE)
      call f_free_ptr(cvdata%COORDINATE_ATOMLIST)
      call f_free_ptr(cvdata%TARG)
      call f_free_ptr(cvdata%cc)
      call f_free_ptr(cvdata%ffcv)
      call f_free_ptr(cvdata%sigma)
      call f_free_ptr(cvdata%sigma_gd)
      call f_free_ptr(cvdata%dt2bym)
      end subroutine 


!     -----------------------------

      subroutine collective_coordinate(ndim,natoms,rxyz,cvdata)  !cc and fcv out  
         implicit none
         integer, intent(in)           :: ndim
         integer, intent(in)           :: natoms
         real(kind=8),intent(in)       :: rxyz(:,:)  
       type(cons_data), intent(inout)  :: cvdata   

         !........local varaibles.....
           integer :: i,j
           integer :: icv   
           integer :: l     
           real(kind=8),allocatable :: gcv(:,:) 
                                                       
                                                      
           real(kind=8) :: ATA(3)  
           real(kind=8) :: ATB(3) 
           real(kind=8) :: ATC(3) 
           real(kind=8) :: ATD(3) 
           real(kind=8) :: DAB
           real(kind=8) :: EBA(3)
           real(kind=8) :: EBC(3)
           real(kind=8) :: COST
           real(kind=8) :: THETA
           real(kind=8) :: ECD(3)
           real(kind=8) :: cost1 
           REAL(KIND=8) :: cost2

           REAL(KIND=8) :: THETA1
           REAL(KIND=8) :: THETA2
           REAL(KIND=8) :: SINT1
           REAL(KIND=8) :: SINT2
           REAL(KIND=8) :: BABC(3)

           REAL(KIND=8) :: COSP

           REAL(KIND=8) :: PHI
           REAL(KIND=8) :: ONE

           REAL(KIND=8) :: CDCB(3)
        !..........................
        do i=1, cvdata%nconstraints
            cvdata%cc(i)=0.d0
        end do 
        do i=1, natoms*ndim 
            do j=1, cvdata%nconstraints
              cvdata%ffcv(i,j)=0.d0
            end do 
        end do 

        !.... ALLOCATION OF LOCAL VARIABLE.....
      allocate( gcv(ndim,natoms) ) 
        !.....................................
      DO ICV = 1, cvdata%nconstraints 
        !   ** bond length **
           IF(cvdata%CONSTRAINTS_TYPE(icv)==1)THEN
            ATA(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,1))
            ATB(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,2))
              CALL DRAB(ATA,ATB,DAB)
                   cvdata%CC(ICV) = DAB
              CALL DRAB_DR(ATA,ATB,GCV)  !DRAB_DR evalutes the grad of bond assocated with atom ata and atb
                 do l=1,2
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+1,icv)=GCV(1,l) 
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+2,icv)=GCV(2,l)
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+3,icv)=GCV(3,l)
                 end do 
     !   **     angle        **
         ELSE IF(cvdata%CONSTRAINTS_TYPE(icv)==2)THEN
                ATA(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,1))
                ATB(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,2))
                ATC(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,3))
                  CALL ERAB(ATB,ATA,EBA)
                  CALL ERAB(ATB,ATC,EBC)
                  COST  = DOT_PRODUCT(EBA(1:3),EBC(1:3))
                  THETA = ACOS(COST)
                               !    IF(THETA<1.D-6)THEN
                               !         WRITE(*,'(A)')' !!! WARNING !!! '
                               !         WRITE(*,*)' THETA =',THETA
                               !         WRITE(*,'(A)')' THETA<1.D-6; THETA SET TO ',EPSI
                               !          THETA = EPSI
                               !      end if 
                  cvdata%CC(ICV)=THETA
                  CALL DTHETA_DR(ATA,ATB,ATC,EBA,EBC,THETA,GCV)
                 do l=1,3
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+1,icv)=GCV(1,l)
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+2,icv)=GCV(2,l)
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+3,icv)=GCV(3,l)
                 end do
         !   **   For dihedral angle  **
          ELSE IF(cvdata%CONSTRAINTS_TYPE(icv)==3)THEN
                 ATA(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,1)) 
                 ATB(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,2))
                 ATC(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,3))
                 ATD(1:ndim)=rxyz(1:ndim,cvdata%COORDINATE_ATOMLIST(icv,4))
                CALL ERAB(ATB,ATA,EBA)
                CALL ERAB(ATB,ATC,EBC)
                CALL ERAB(ATC,ATD,ECD)
              COST1  =      DOT_PRODUCT(EBA(1:3),EBC(1:3))
              COST2  = -1.0*DOT_PRODUCT(ECD(1:3),EBC(1:3))
             THETA1 = ACOS(COST1)
             THETA2 = ACOS(COST2)
             SINT1  = SIN(THETA1)
             SINT2  = SIN(THETA2)
              CALL ACRB(EBA,EBC,BABC) !evaluates the cross product of two vectors
              EBC(1:3) = -1.0*EBC(1:3)
              CALL ACRB(EBC,ECD,CDCB)
              EBC(1:3) = -1.0*EBC(1:3)
              COSP   = DOT_PRODUCT(BABC(1:3),CDCB(1:3))
              COSP   = COSP/SINT1/SINT2
              PHI    = ACOS(COSP)
              ONE    = 1.0
              ONE    = ONE*SIGN(ONE,DOT_PRODUCT(EBA(1:3),CDCB(1:3)))
                PHI    = PHI*ONE
                          !   IF(ABS(PHI)<1.0D-6)THEN
                          !   WRITE(*,*)' !!! WARNING !!! '
                          !   WRITE(*,*)' PHI =', PHI,'; PHI<1e-6; PHI SET TO ',EPSI
                          !   PHI  = EPSI
                          !   END IF
               cvdata%CC(ICV)=PHI
               CALL DPHI_DR(ATA,ATB,ATC,ATD,EBA,EBC,ECD, &
                     PHI,THETA1,THETA2,GCV)
                 do l=1, 4
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+1,icv)=GCV(1,l)
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+2,icv)=GCV(2,l)
                   cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+3,icv)=GCV(3,l)
                 end do
           END IF 
      END DO 


      end subroutine collective_coordinate  

       subroutine calculate_sigma(ndim,natoms,rxyz,cvdata)  !sigma , sigma_dr  out
        implicit none
           integer, intent(in)              ::  ndim 
           integer, intent(in)              ::  natoms
           real(kind=8), intent(in)         ::  rxyz(ndim,natoms)
           type(cons_data), intent(inout)  ::  cvdata 
        !-----     local varaiable ----
              integer :: icv,i,j
              integer :: l


      call  collective_coordinate(ndim,natoms,rxyz,cvdata)  
        do i =1, cvdata%nconstraints
           cvdata%sigma(i,1)=0.d0
        end do 

           do icv=1 , cvdata%nconstraints
                 if (cvdata%CONSTRAINTS_TYPE(icv)==1) then
                    cvdata%sigma(icv,1) = (cvdata%cc(icv) )**2.d0 - (cvdata%targ(icv) )**2.d0
                 else if (cvdata%CONSTRAINTS_TYPE(icv)==2) then
                    cvdata%sigma(icv,1) = (cvdata%cc(icv) )**2.d0 - (cvdata%targ(icv) )**2.d0
                 else if  (cvdata%CONSTRAINTS_TYPE(icv)==3) then
                    cvdata%sigma(icv,1) = (cvdata%cc(icv) )**2.d0 - (cvdata%targ(icv) )**2.d0
                 end  if 
           end do 
      

       do i=1, natoms*ndim
          do j=1, cvdata%nconstraints  
             cvdata%sigma_gd(i,j) = 0.d0
          end do 
       end do  

           do icv=1 , cvdata%nconstraints
                 if (cvdata%CONSTRAINTS_TYPE(icv)==1) then

                   do l=1,2
                     cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+1,icv)=&
                     2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+1,icv)

                     cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+2,icv)=&
                     2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+2,icv)

                     cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+3,icv)=&
                     2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+3,icv)
                   end do


                 else if (cvdata%CONSTRAINTS_TYPE(icv)==2) then

                   do l=1,3
                     cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+1,icv)=&
                     2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+1,icv)

                     cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+2,icv)=&
                     2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+2,icv)

                     cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+3,icv)=&
                     2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+3,icv)
                   end do

                 else if  (cvdata%CONSTRAINTS_TYPE(icv)==3) then


                   do l=1,4
                    cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+1,icv)=&
                    2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+1,icv)

                    cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+2,icv)=&
                    2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+2,icv)

                    cvdata%sigma_gd(3*( cvdata%COORDINATE_ATOMLIST(icv,l) -1)+3,icv)=&
                    2.d0*cvdata%cc(icv)*cvdata%ffcv(3*(cvdata%COORDINATE_ATOMLIST(icv,l)-1)+3,icv)
                   end do

                 end  if
           end do



       end subroutine



!          ________________________________

        subroutine  build_Amatrix(ndim,natoms,dt2bym,dsdr,sdr,AAA,cvdata)  ! amat out
        implicit none
        integer,intent(in)            :: ndim
        integer,intent(in)            :: natoms
        type(cons_data),intent(inout)  ::  cvdata

        real(kind=8),intent(in)       :: dt2bym(:)
        real(kind=8),intent(inout)    :: AAA(cvdata%nconstraints,cvdata%nconstraints)
        real(kind=8),intent(in)       :: dsdr(natoms*ndim,cvdata%nconstraints)
        real(kind=8),intent(in)       :: sdr(natoms*ndim,cvdata%nconstraints)

        integer :: i,j,k,l,p,q !local varible
           DO i=1,cvdata%nconstraints
              DO j=1,cvdata%nconstraints
                 AAA(i,j)=0.D0
                  DO k=1,natoms*ndim
                    AAA(i,j) = AAA(i,j) + sdr(k,i)*dsdr(k,j)*dt2bym(k)
                  ENDDO
              ENDDO
           ENDDO

!           print *,"A matrix "
!           do i=1,cvdata%nconstraints
!               print*,"A(",i,"1:",cvdata%nconstraints,")", AAA(i,1:cvdata%nconstraints)
!           end do

!stop
        end subroutine
!_____________________________________________________________
     subroutine   update_coodinates(ndim,dt2bym,natoms,lll,sdr,rxyz,cvdata)   !  rxyz inout
     implicit none
     type(cons_data)  ::  cvdata
     integer , intent(in)          :: ndim
     real(kind=8), intent(in)      ::dt2bym(:)
     integer , intent(in)          :: natoms
     real(kind=8), intent(in)      :: lll(cvdata%nconstraints)
     real(kind=8), intent(in)      :: sdr(natoms*ndim,cvdata%nconstraints)
     real(kind=8), intent(inout)   :: rxyz(:,:)

     !... local variable........
     integer :: i
     integer :: icv
     integer :: l
    do i=1, natoms
        do icv=1,cvdata%nconstraints
           rxyz(1,i) = rxyz(1,i) - dt2bym(3*(i-1)+1)*lll(icv)*sdr( 3*(i-1)+1 ,icv)      
           rxyz(2,i) = rxyz(2,i) - dt2bym(3*(i-1)+2)*lll(icv)*sdr( 3*(i-1)+2 ,icv)
           rxyz(3,i) = rxyz(3,i) - dt2bym(3*(i-1)+3)*lll(icv)*sdr( 3*(i-1)+3 ,icv)
        end do
    end do
    end subroutine


        
          SUBROUTINE DRAB(VECTA,VECTB,DAB)   
     !     -----------------------------------------------------------------
     !     Evaluates the distance |R_AB| when 
     !     vector_a and vector_b are given
     !     Creation: NNN (12.01.05)
     !     -----------------------------------------------------------------
           IMPLICIT NONE
           REAL(kind=8), INTENT(IN)   :: VECTA(3), VECTB(3)
           REAL(kind=8), INTENT(OUT)  :: DAB
           REAL(kind=8)  :: RAB(3)
           RAB(1)  = VECTB(1)-VECTA(1)
           RAB(2)  = VECTB(2)-VECTA(2)
           RAB(3)  = VECTB(3)-VECTA(3)
           DAB     = RAB(1) * RAB(1) + RAB(2) * RAB(2) +  RAB(3) * RAB(3)
           DAB     = DSQRT(DAB)
          END  SUBROUTINE DRAB
     !___________________________________
         SUBROUTINE DRAB_DR(ATA,ATB,GCV)
     !     -----------------------------------------------------------------
     !     Evaluates the gradient of distance |R_AB| when 
     !     vector_a and vector_b are given
     !     -----------------------------------------------------------------
           IMPLICIT NONE
           REAL*8 :: ATA(3),ATB(3),GCV(3,2)
     !
           REAL*8 :: EAB(3)  
     !
           CALL  ERAB(ATA,ATB,EAB)
           GCV(1:3,1)=-EAB(1:3)
           GCV(1:3,2)= EAB(1:3)
        END SUBROUTINE DRAB_DR
     !___________________________________

      SUBROUTINE ERAB(VECTA,VECTB,EAB)
      !     -----------------------------------------------------------------
      !     Evaluates the unit_distance_vector e_AB when 
      !     vector_a and vector_b are given
      !     Creation: NNN (12.01.05)
      !     -----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, INTENT(IN)   :: VECTA(3), VECTB(3)
      REAL*8, INTENT(OUT)  :: EAB(3)
      REAL*8   :: RAB(3),DAB
                       
      RAB(1)  = VECTB(1)-VECTA(1)   
      RAB(2)  = VECTB(2)-VECTA(2)
      RAB(3)  = VECTB(3)-VECTA(3)
      DAB     = RAB(1) * RAB(1) + RAB(2) * RAB(2) + RAB(3) * RAB(3)
      DAB     = DSQRT(DAB)
      EAB(1)  = RAB(1)/DAB
      EAB(2)  = RAB(2)/DAB
      EAB(3)  = RAB(3)/DAB
      END subroutine ERAB


!
      SUBROUTINE DTHETA_DR(ATA,ATB,ATC,EBA,EBC,THETA,GCV)
!     -----------------------------------------------------------------
!     Evaluates the gradient of angle theta  
!     -----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 :: ATA(3),ATB(3),ATC(3),EBA(3),EBC(3),GCV(3,4),THETA
!
      REAL*8 :: RBA,RBC,COST,SINT
      COST=COS(THETA)
      SINT=SIN(THETA)
      CALL DRAB(ATA,ATB,RBA)
      CALL DRAB(ATB,ATC,RBC)
      GCV(1:3,1)=-(EBC(1:3)-EBA(1:3)*COST)/RBA/SINT
      GCV(1:3,2)= (EBC(1:3)-EBA(1:3)*COST)/RBA/SINT + &
                  (EBA(1:3)-EBC(1:3)*COST)/RBC/SINT
      GCV(1:3,3)=-(EBA(1:3)-EBC(1:3)*COST)/RBC/SINT
      END SUBROUTINE
!

!
      SUBROUTINE ACRB(VECTA,VECTB,CPAB)
!     -----------------------------------------------------------------
!     Evaluates the the cross product of vector_a and vector_b
!     when vector_a and vector_b are given
!     Creation: NNN (12.01.05)
!-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: VECTA(3), VECTB(3)
      REAL*8, INTENT(OUT) :: CPAB(3)
      CPAB(1) = VECTA(2)*VECTB(3)-VECTA(3)*VECTB(2)
      CPAB(2) = VECTA(3)*VECTB(1)-VECTA(1)*VECTB(3)
      CPAB(3) = VECTA(1)*VECTB(2)-VECTA(2)*VECTB(1)
      END subroutine
!

        !....................................


      SUBROUTINE DPHI_DR(ATA,ATB,ATC,ATD,EBA,EBC,ECD, &
                         PHI,THETA1,THETA2,GCV)
!     -----------------------------------------------------------------
!     Evaluates the gradient of dihedral phi  
!     -----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 :: ATA(3),ATB(3),ATC(3),ATD(3),EBA(3),EBC(3),ECD(3), &
              GCV(3,4),PHI,THETA1,THETA2
!
      REAL*8 :: FACT,FAC(3),FAC2(3),COST1,COST2,SINT1,SINT2,SINP, &
              RBA,RBC,RCD
!
      COST1       = COS(THETA1)
      SINT1       = SIN(THETA1)
      COST2       = COS(THETA2)
      SINT2       = SIN(THETA2)
      SINP        = SIN(PHI)
!
      CALL DRAB(ATA,ATB,RBA)
      CALL DRAB(ATC,ATB,RBC)
      CALL DRAB(ATC,ATD,RCD)
!
      FACT      = DOT_PRODUCT(EBA(1:3),ECD(1:3))
!
      GCV(1:3,1)=(ECD(1:3)-FACT*EBA(1:3))/RBA
      FAC(1:3)  =-(EBC(1:3)-EBA(1:3)*COST1)/RBA/SINT1
      GCV(1:3,1)=GCV(1:3,1)-FAC(1:3)*SINT1*COST2- &
                            FAC(1:3)*(COST1/SINT1)*(COST1*COST2+FACT)
!
      GCV(1:3,2)=-(ECD(1:3)-FACT*EBA(1:3))/RBA
      FAC(1:3)  =(EBC(1:3)-EBA(1:3)*COST1)/RBA/SINT1 + &
                 (EBA(1:3)-EBC(1:3)*COST1)/RBC/SINT1
      FAC2(1:3) =-(ECD(1:3)+EBC(1:3)*COST2)/RBC/SINT2
      GCV(1:3,2)=GCV(1:3,2)-SINT1*COST2*FAC(1:3)- &
                  SINT2*COST1*FAC2(1:3)-   &
                  ((COST1/SINT1)*FAC(1:3)+ &
                   (COST2/SINT2)*FAC2(1:3))* &
                   (COST1*COST2+FACT)
      GCV(1:3,3)=(-EBA(1:3)+FACT*ECD(1:3))/RCD
      FAC(1:3)  =-(EBA(1:3)-EBC(1:3)*COST1)/RBC/SINT1
      FAC2(1:3) =(ECD(1:3)+EBC(1:3)*COST2)/RBC/SINT2 - &
                 (EBC(1:3)+ECD(1:3)*COST2)/RCD/SINT2
      GCV(1:3,3)=GCV(1:3,3)-SINT1*COST2*FAC(1:3)- &
                 SINT2*COST1*FAC2(1:3)- &
                 ((COST1/SINT1)*FAC(1:3)+ &
                 (COST2/SINT2)*FAC2(1:3))* &
                 (COST1*COST2+FACT)
      GCV(1:3,4)=(EBA(1:3)-FACT*ECD(1:3))/RCD
      FAC(1:3)  =-(-EBC(1:3)-ECD(1:3)*COST2)/RCD/SINT2
      GCV(1:3,4)=GCV(1:3,4)-FAC(1:3)*SINT2*COST1- &
                            FAC(1:3)*(COST2/SINT2)*(COST1*COST2+FACT)
      FACT      = -SINT1*SINT2*SINP
      GCV(1:3,1)= GCV(1:3,1)/FACT
      GCV(1:3,2)= GCV(1:3,2)/FACT
      GCV(1:3,3)= GCV(1:3,3)/FACT
      GCV(1:3,4)= GCV(1:3,4)/FACT
      END subroutine


!-------------------------------------------------------------------
      subroutine evalute_dt2bym(ndim,natoms,dt,amass,ddtm)
      implicit none
         integer,intent(inout)        ::  ndim               
         integer,intent(in)           ::  natoms             
         real(kind=8),intent(in)      ::  dt          
         real(kind=8),intent(in)      ::  amass(natoms)
         real(kind=8),intent(inout)   ::  ddtm(natoms*ndim)
         integer :: i
      do i=1,natoms
        ddtm(3*(i-1)+1)=dt*dt/2.d0/amass(i)
        ddtm(3*(i-1)+2)=dt*dt/2.d0/amass(i)
        ddtm(3*(i-1)+3)=dt*dt/2.d0/amass(i)
      end do
      
      end subroutine 

subroutine conjugate_gradient(ldd,AA,kk,kkd,noptmax,bbb,cvdata)  
implicit none 
  type(cons_data)  ::  cvdata
  integer :: kk,kkd,noptmax
  real(kind=8) ::  ldd(cvdata%nconstraints,noptmax), AA(cvdata%nconstraints, cvdata%nconstraints),bbb(cvdata%nconstraints)
  !local variable
  real (kind=8) :: dotpr,rsold1,rsold,alpha,rsnew,beta
  real(kind=8),allocatable   :: r(:),test1(:,:),ap(:),sol(:),ld1(:,:),AAAA1(:,:),p1(:,:),ap1(:,:)
  integer :: i,j, ii,jj

allocate(r(cvdata%nconstraints))
allocate(sol(cvdata%nconstraints))
allocate(AAAA1(cvdata%nconstraints, cvdata%nconstraints ))
allocate(ld1(cvdata%nconstraints,1))
allocate(test1(cvdata%nconstraints,1))
allocate(p1(cvdata%nconstraints,1))
allocate(ap1(cvdata%nconstraints,1))
 cvdata%fflage=0 
 do i=1, cvdata%nconstraints
  p1(i,1) = 0.d0
  ap1(i,1)= 0.d0
  test1(i,1)=0.d0
  ld1(i,1)=0.d0
 end do  
          call dcopy(cvdata%nconstraints*cvdata%nconstraints,AA,1,AAAA1,1)
          call dgemm('N','N' ,cvdata%nconstraints,1,cvdata%nconstraints,1.d0 &
                      ,AAAA1,cvdata%nconstraints ,ld1,cvdata%nconstraints &
                         ,0.d0 ,test1,cvdata%nconstraints)
 do i=1, cvdata%nconstraints
  r(i)=0.d0
  r(i) = bbb(i)  - test1(i,1)
 end do 

 call dcopy(cvdata%nconstraints,r,1,p1(:,1),1)
 rsold=0.d0
 rsold1=0.d0
 do ii=1,cvdata%nconstraints
   rsold= rsold + r(ii)*r(ii)
 end do 
 kk=0
CG_loop:do
   kk=kk+1
   if (kk.gt.noptmax) stop
            call dcopy(cvdata%nconstraints*cvdata%nconstraints,AA,1,AAAA1,1)
            call dgemm('N','N' ,cvdata%nconstraints,1,cvdata%nconstraints,1.d0 &
                      ,AAAA1,cvdata%nconstraints ,p1,cvdata%nconstraints &
                         ,0.d0 ,ap1,cvdata%nconstraints)
   dotpr=0.d0
   do ii=1,cvdata%nconstraints
    dotpr= dotpr + p1(ii,1)*ap1(ii,1)
   end do
   alpha= rsold/dotpr
   do i=1,cvdata%nconstraints   
    ld1(i,1)=ld1(i,1) + alpha*p1(i,1)
   end do 
   do i=1, cvdata%nconstraints   
    ldd(i,kk)=ld1(i,1)
   end do
   do i=1, cvdata%nconstraints 
    r(i)=r(i)- alpha*ap1(i,1)
   end do 
   rsnew=0.d0

   do ii=1,cvdata%nconstraints
     rsnew= rsnew + r(ii)*r(ii)
   end do

    if (cvdata%algorithm==2) then
       if (sqrt(rsnew)<1.0e-10) exit CG_loop
    end if 
   if (cvdata%algorithm==3) then
       if (kk==kkd)  exit CG_loop   
   end if
   if (sqrt(rsnew)<1.0e-10) then
      cvdata%fflage=1                 
      exit CG_loop
   end if 
   beta=rsnew/rsold

   do i=1, cvdata%nconstraints   
     p1(i,1)=r(i) + beta*p1(i,1)
   end do 

   rsold=rsnew

   if (kk>noptmax-2) then
    print*, "optimization cannot be completed by CG"
    print*, "CG has exceeded max. number of iterations" 
    stop
   endif 
end do CG_loop
 deallocate(r)
 deallocate(sol)
 deallocate(AAAA1)
 deallocate(ld1)
 deallocate(test1)
 deallocate(p1)
 deallocate(ap1)
end subroutine 

!===================================
subroutine DIIS(kk,ldd,noptmax,cvdata)
implicit none 
 type(cons_data),intent (inout)  ::  cvdata
 integer :: kk,noptmax 
 real(kind=8) :: ldd(cvdata%nconstraints,noptmax)
 !local variable 
 real(kind=8),allocatable   ::  pp(:,:),bib(:,:),ccc(:),ccc1(:)  
 real(kind=8),allocatable   ::  sol(:),dpp(:)    
 integer :: l,k,i,j,flag,kkd,ii
 integer :: INFO
 integer,allocatable  :: ipiv3(:)
 integer :: icount
 kkd=kk
 !DIIS Algorithm
allocate(bib(kkd,kkd))
allocate(ccc(kkd))
allocate(ccc1(kkd))
allocate(ipiv3(kkd)) 
allocate(dpp(cvdata%nconstraints))
allocate(pp(cvdata%nconstraints,noptmax-1)) 
allocate(sol(cvdata%nconstraints)) 
icount=0
diis_loop:do
   if (cvdata%fflage==1) then
    exit diis_loop
   endif 
     do l=1, kk-1
        do i=1, cvdata%nconstraints
         pp(i,l) = ldd(i,l+1) - ldd(i,l)
        end do
     end do
     bib=0.d0
          do i=1,kk-1
              do j=1,kk-1
                     bib(i,j)= 0.d0
                     do ii=1,cvdata%nconstraints
                        bib(i,j)= bib(i,j) + pp(ii,i)*pp(ii,j)
                     end do
              end do
          end do

           do l=1,kk
              bib(kk,l)=-1
           end do
           do l=1,kk
              bib(l,kk)=-1
           end do
           bib(kk,kk)=0.d0
           ccc=0.d0
           ccc(kk)=-1
           ccc1=ccc
        call  DGESV(kk,1,bib,kk,ipiv3,ccc1,kk,INFO) 
        sol=0.d0
          do i=1, kk-1
            do j=1, cvdata%nconstraints
             sol(j)= sol(j) + ccc1(j)*ldd(j,i)
            end do 
          end do
        dpp=0.d0
        do i=1,kk-1
           do j=1, cvdata%nconstraints
           dpp(j)= dpp(j) + ccc1(i)*pp(j,i)
           end do 
        end do
          flag = 0
          do i=1,cvdata%nconstraints
            if (abs(dpp(i) ) <= 1.0e-10 ) then
            flag = flag+1
            end if
          end do
          if (flag == cvdata%nconstraints) exit diis_loop
          do i=1, kk-1
             do j=1, cvdata%nconstraints
               ldd(j,i)  =  ldd(j,i+1)
             end do 
           end do 

           do i=1, cvdata%nconstraints            
             ldd(i,kk) =  sol(i)
           end do 
   icount = icount + 1
   if (icount>noptmax-2) then
    print*, "Optimization cannot be completed by DISS"
    print*, "DIIS has exceeded the maxi. number of iteration"
    stop
   endif


end do diis_loop

 deallocate(bib)
 deallocate(ccc)
 deallocate(ccc1)
 deallocate(ipiv3)
 deallocate(dpp)
 deallocate(pp) 
 deallocate(sol) 
end subroutine DIIS

!  pure subroutine nullify_cvdata_data(cvdata)
!    !use m_pawang, only: pawang_nullify
!    implicit none
!    type(cons_data), intent(out)  ::  cvdata
!    cvdata%nconstraints =-1
!    cvdata%algorithm =-1
!    cvdata%diis_size = -1
!    cvdata%const_conv = -1.d0
!    cvdata%fflage = -1
!
!  end subroutine

!......................................
end module  !_____________________________________

!!      subroutine cvdata_new(cvdata)
!!      implicit none
!!      type(cons_data), pointer  ::  cvdata
!!      type(cons_data), pointer, save  ::  intern
!!      end subroutine

!----------------------------------------------------------------------------------------------------------------------------
   subroutine shake(ndim,natoms,dt,amass,rxyz,drxyz0,cvdata)
      use module_base
      use constraints


      implicit none 
         type(cons_data),intent(inout)   ::  cvdata
         integer,intent(inout)           ::  ndim             
         integer,intent(inout)           ::  natoms          
         real(kind=8),intent(inout)      ::  dt          
         real(kind=8),intent(inout)      ::  amass(natoms)
         real(kind=8),intent(inout)      ::  rxyz(ndim,natoms),drxyz0(ndim,natoms)
      

         !....... local variables ..........
         integer :: i,flag,kkd 
         integer :: k,kk,l,j   
         integer :: iopt=0
         integer, parameter :: noptmax=999
 REAL(KIND=8), DIMENSION(:,:), POINTER ::  AA(:,:),bbb(:) , dsigma_gd(:,:) ,ld(:),ldd(:,:)

         !..................................
         integer :: N,NRHS 
         integer :: LDA,LDB
         integer :: INFO
 REAL(KIND=8), DIMENSION(:,:), POINTER :: IPIV(:)

         cvdata%fflage=0 
         kkd=cvdata%diis_size
         N=cvdata%nconstraints
         NRHS = 1
         LDA = cvdata%nconstraints
         LDB = cvdata%nconstraints


        IPIV  = f_malloc_ptr([N],id='IPIV')
        bbb  = f_malloc_ptr([cvdata%nconstraints],id='bbb')
        AA  = f_malloc_ptr([cvdata%nconstraints,cvdata%nconstraints],id='AA')
        ld  = f_malloc_ptr([cvdata%nconstraints],id='ld')
        ldd  = f_malloc_ptr([cvdata%nconstraints,noptmax],id='ldd')
        dsigma_gd  = f_malloc_ptr([natoms*ndim,cvdata%nconstraints],id='dsigma_gd')



     call evalute_dt2bym(ndim,natoms,dt,amass,cvdata%dt2bym)   
     call calculate_sigma(ndim,natoms,drxyz0,cvdata)   
     call dcopy(natoms*ndim*cvdata%nconstraints,cvdata%sigma_gd,1,dsigma_gd,1) 
     call calculate_sigma(ndim,natoms,drxyz0,cvdata)  

       iopt=0
       opt_loop: do
         iopt=iopt+1
          if(iopt.gt.noptmax) &
                  stop 'constraint optimization is not completed'
          call calculate_sigma(ndim,natoms,rxyz,cvdata)         
          flag = 0
          do i=1,cvdata%nconstraints
          if (abs(cvdata%sigma(i,1))<= cvdata%const_conv) then
          flag = flag+1
          end if
          end do 
          if (flag == cvdata%nconstraints) exit opt_loop
          call dcopy(cvdata%nconstraints,cvdata%sigma(:,1),1,bbb,1)
          if (cvdata%algorithm==1) call  build_Amatrix(ndim,natoms,cvdata%dt2bym,dsigma_gd,cvdata%sigma_gd,AA,cvdata)  
          !if (cvdata%algorithm/=1) call  build_Amatrix(ndim,natoms,cvdata%dt2bym,cvdata%sigma_gd,cvdata%sigma_gd,AA,cvdata)  
          if (cvdata%algorithm/=1) call  build_Amatrix(ndim,natoms,cvdata%dt2bym,dsigma_gd,dsigma_gd,AA,cvdata)  
          do i=1,cvdata%nconstraints
           ld(i)=0.d0
             do j=1, noptmax
             ldd(i,j) =  0.d0
             end do 
          end do 
          !linear sys. of eq. AX=B
         !================
        
          
       if (cvdata%algorithm==1) then
         
          call  DGESV( N,NRHS, AA, LDA, IPIV, cvdata%sigma, LDB, INFO  ) 
          IF( INFO.GT.0 ) THEN
            WRITE(*,*)'The diagonal element of the triangular factor of A,'
            WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
            WRITE(*,*)'A is singular; the solution could not be computed.'
            STOP
          END IF
          call dcopy(cvdata%nconstraints,cvdata%sigma,1,ld,1)
        end if 
        if (cvdata%algorithm==2) then
          call conjugate_gradient(ldd,AA,kk,kkd,noptmax,bbb,cvdata) 
          call dcopy(cvdata%nconstraints,ldd(:,kk),1,ld,1)
        end if 
        if (cvdata%algorithm==3) then
          call conjugate_gradient(ldd,AA,kk,kkd,noptmax,bbb,cvdata) 
          call DIIS(kk,ldd,noptmax,cvdata) 
          call dcopy(cvdata%nconstraints,ldd(:,kk),1,ld,1) 
        end if 
        !=================
          call  update_coodinates(ndim,cvdata%dt2bym,natoms,ld,dsigma_gd,rxyz,cvdata)

    
     end do opt_loop 
 
     call f_free_ptr(IPIV)
      call f_free_ptr(bbb)
      call f_free_ptr(AA)
      call f_free_ptr(ld)
      call f_free_ptr(ldd)
      call f_free_ptr(dsigma_gd)

      print*,"Number of Shake Iteration per MD:", iopt


    end subroutine shake


!------------------------------------------------------------

      subroutine  rattle(ndim,natoms,dt,amass,fxyz,drxyz,rxyz,vxyz,cvdata) 
      use module_base
      use constraints
      implicit none
          type(cons_data),intent(inout)  ::  cvdata
!         integer,intent(inout)           ::  istep
         integer,intent(inout)           ::  ndim        
         integer,intent(inout)           ::  natoms             
         real(kind=8),intent(inout)      ::  dt          
         real(kind=8),intent(inout)      ::  amass(natoms)
         real(kind=8),intent(inout)      ::  fxyz(ndim,natoms)
         real(kind=8),intent(inout)      ::  drxyz(ndim,natoms) 
         real(kind=8),intent(inout)      ::  rxyz(ndim,natoms)
         real(kind=8),intent(inout)      ::  vxyz(ndim,natoms)
         integer, parameter :: noptmax=999


     !-------local variables-------
     integer :: i,k,icv,p,kk,j

 REAL(KIND=8), DIMENSION(:,:), POINTER :: vxyz0(:)
 REAL(KIND=8), DIMENSION(:,:), POINTER ::ffxyz(:)
 REAL(KIND=8), DIMENSION(:,:), POINTER :: rhs1(:)
 REAL(KIND=8), DIMENSION(:,:), POINTER ::  rhs2(:)
 REAL(KIND=8), DIMENSION(:,:), POINTER :: rhs(:,:)
 REAL(KIND=8), DIMENSION(:,:), POINTER ::AA(:,:)
!
 REAL(KIND=8), DIMENSION(:,:), POINTER :: ld(:)
 REAL(KIND=8), DIMENSION(:,:), POINTER ::  ldd(:,:)

    !----------------------------
         integer :: N,NRHS
         integer :: LDA,LDB
         integer :: INFO
!         integer,allocatable  :: IPIV(:)
 REAL(KIND=8), DIMENSION(:,:), POINTER ::  IPIV(:)
         N=cvdata%nconstraints
         NRHS = 1
         LDA = cvdata%nconstraints
         LDB = cvdata%nconstraints

        IPIV  = f_malloc_ptr([N],id='IPIV')
        ldd  = f_malloc_ptr([cvdata%nconstraints,noptmax],id='ldd')
        AA  = f_malloc_ptr([cvdata%nconstraints, cvdata%nconstraints],id='AA')
        ld  = f_malloc_ptr([cvdata%nconstraints],id='ld')
        vxyz0  = f_malloc_ptr([natoms*ndim],id='vxyz0')
        ffxyz  = f_malloc_ptr([natoms*ndim],id='ffxyz')
        rhs1  = f_malloc_ptr([cvdata%nconstraints],id='rhs1')
        rhs2  = f_malloc_ptr([cvdata%nconstraints],id='rhs2')
        rhs  = f_malloc_ptr([cvdata%nconstraints,1],id='rhs')



     !---------------------------
   
     
      do i=1 , natoms
         do j=1, 3
            vxyz0(3*(i-1)+j) = 0.d0
            vxyz0(3*(i-1)+j) = vxyz(j,i)
           !vxyz0(3*(i-1)+2) = vxyz(2,i)
           !vxyz0(3*(i-1)+3) = vxyz(3,i)
         end do  
      end do
           
           do icv=1 , cvdata%nconstraints
                 if (cvdata%CONSTRAINTS_TYPE(icv)==1) then
                        do i=1 , 2 
                       p=cvdata%COORDINATE_ATOMLIST(icv,i)
                        vxyz0(3*(p-1)+1) = (rxyz(1,p) - drxyz(1,p) ) /dt
                        vxyz0(3*(p-1)+2) = (rxyz(2,p) - drxyz(2,p) ) /dt
                        vxyz0(3*(p-1)+3) = (rxyz(3,p) - drxyz(3,p) ) /dt
                        end do
                 else if (cvdata%CONSTRAINTS_TYPE(icv)==2) then
                        do i=1 , 3
                       p=cvdata%COORDINATE_ATOMLIST(icv,i)
                        vxyz0(3*(p-1)+1) = (rxyz(1,p) - drxyz(1,p) ) /dt
                        vxyz0(3*(p-1)+2) = (rxyz(2,p) - drxyz(2,p) ) /dt
                        vxyz0(3*(p-1)+3) = (rxyz(3,p) - drxyz(3,p) ) /dt

                        end do
                 else if  (cvdata%CONSTRAINTS_TYPE(icv)==3) then
                        do i=1 , 4
                       p=cvdata%COORDINATE_ATOMLIST(icv,i)
                        vxyz0(3*(p-1)+1) =1* (rxyz(1,p) - drxyz(1,p) ) /dt
                        vxyz0(3*(p-1)+2) =1* (rxyz(2,p) - drxyz(2,p) ) /dt
                        vxyz0(3*(p-1)+3) =1* (rxyz(3,p) - drxyz(3,p) ) /dt
                        end do
                 end  if
           end do

      do i=1 , natoms
      ffxyz(3*(i-1)+1) = fxyz(1,i)*amass(i) 
      ffxyz(3*(i-1)+2) = fxyz(2,i)*amass(i) 
      ffxyz(3*(i-1)+3) = fxyz(3,i)*amass(i) 
      end do

      call calculate_sigma(ndim,natoms,rxyz,cvdata) 
      do i=1, cvdata%nconstraints
       rhs1(i)=0.d0
        do k=1,ndim*natoms
          rhs1(i) = rhs1(i) + dt*vxyz0(k)*cvdata%sigma_gd(k,i) 
        end do 
      end do 

      do i=1, cvdata%nconstraints
       rhs2(i)=0.d0
        do k=1,ndim*natoms
            rhs2(i) = rhs2(i) + cvdata%dt2bym(k)*ffxyz(k)*cvdata%sigma_gd(k,i)
        end do 
      end do
      do i=1, cvdata%nconstraints
      rhs(i,1)=0.d0
      rhs(i,1) = rhs1(i)  +  rhs2(i)
      end do


      call  build_Amatrix(ndim,natoms,cvdata%dt2bym,cvdata%sigma_gd,cvdata%sigma_gd,AA,cvdata)  
      do i=1 , cvdata%nconstraints
      ld(i)=0.d0
      end do 
      ! linear sys. of eq. AX=B
   if (cvdata%algorithm==1) then
      call  DGESV( N,NRHS, AA, LDA, IPIV,rhs(:,1) , LDB, INFO  ) 
      call dcopy(cvdata%nconstraints*1,rhs(:,1),1,ld,1)  
   end if


   if (cvdata%algorithm==2) then
      call conjugate_gradient(ldd,AA,kk,cvdata%diis_size,noptmax,rhs(:,1),cvdata) 
      call dcopy(cvdata%nconstraints*1,ldd(:,kk),1,ld,1) 
   end if

   if (cvdata%algorithm==3) then
      call conjugate_gradient(ldd,AA,kk,cvdata%diis_size,noptmax,rhs(:,1),cvdata) 
      call DIIS(kk,ldd,noptmax,cvdata)
      call dcopy(cvdata%nconstraints*1,ldd(:,kk),1,ld,1)  
   end if

       do p=1,natoms
           vxyz0(3*(p-1)+1) = vxyz0(3*(p-1)+1) + cvdata%dt2bym(3*(p-1)+1)*ffxyz(3*(p-1)+1)/dt
           vxyz0(3*(p-1)+2) = vxyz0(3*(p-1)+2) + cvdata%dt2bym(3*(p-1)+2)*ffxyz(3*(p-1)+2)/dt
           vxyz0(3*(p-1)+3) = vxyz0(3*(p-1)+3) + cvdata%dt2bym(3*(p-1)+3)*ffxyz(3*(p-1)+3)/dt
       end do 

      do i=1, natoms
        do icv=1,cvdata%nconstraints
           vxyz0(3*(i-1)+1) = vxyz0(3*(i-1)+1) - ld(icv)*cvdata%sigma_gd( 3*(i-1)+1 ,icv)*cvdata%dt2bym(3*(i-1)+1)/dt
           vxyz0(3*(i-1)+2) = vxyz0(3*(i-1)+2) - ld(icv)*cvdata%sigma_gd( 3*(i-1)+2 ,icv)*cvdata%dt2bym(3*(i-1)+2)/dt
           vxyz0(3*(i-1)+3) = vxyz0(3*(i-1)+3) - ld(icv)*cvdata%sigma_gd( 3*(i-1)+3 ,icv)*cvdata%dt2bym(3*(i-1)+3)/dt
        end do
      end do

      do i=1 , natoms
       do j=1 , 3
       vxyz(j,i)=0.d0
       vxyz(j,i)= vxyz0(3*(i-1)+j) 
       end do 
      end do

      call f_free_ptr(vxyz0)
      call f_free_ptr(ffxyz)
      call f_free_ptr(rhs1)
      call f_free_ptr(rhs2)
      call f_free_ptr(rhs)
      call f_free_ptr(AA)
      call f_free_ptr(ld)
      call f_free_ptr(IPIV)
      call f_free_ptr(ldd)


     end subroutine rattle
!----------------------------------------------------- End ---------------------------------------------------------------
