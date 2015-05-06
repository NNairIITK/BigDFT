!parameters are from
!T Deutsch et al 1995 J. Phys.: Condens. Matter 7 6407. doi:10.1088/0953-8984/7/32/007

MODULE module_atomes_rgl
!~~~~~~~~~~~~~~~~~~~
!-----------------------------------------------------------------------
!-----    Modifications:
!-----    ~~~~~~~~~~~~~~
!-- 23/ 4/2001 : cochonné le module initial en enlevant « nat », etc. [F.L.]
!---------------
!-----  Auteur : Frédéric LANÇON --- 17 septembre 1997
!-----------------------------------------------------------------------
!.......................................................................
   IMPLICIT NONE
!.......................................................................

   INTEGER, PARAMETER :: ntype_max = 1  ! nombre maximum de types d'atomes.

   INTEGER, SAVE :: ntype      ! Nombre de types d'atomes.
   CHARACTER(LEN=8), SAVE, DIMENSION(ntype_max) :: ATOME_NOM
   ! Nom des types d'atomes.
!........................................................................

END MODULE module_atomes_rgl

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

MODULE Interactions_de_paire_rgl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~
!.......................................................................

  USE module_atomes_rgl
!.......................................................................
  IMPLICIT NONE
!.......................................................................
!----- Paramètres pour interactions de paire.

!--  Pour chaque type de paire :
  ! E0_paire est le parametre multiplicatif du terme de paire [énergie] ;
  ! R0_paire contientles distances de premier voisins [longueur] ;
  ! P_paire est l'exposant pour une loi puissance ou exponentielle.
  DOUBLE PRECISION E0_paire, R0_paire, P_paire
  COMMON/SUB_PAIRE/ E0_paire(ntype_max, ntype_max)         &
       &          , R0_paire(ntype_max, ntype_max)   &
       &          ,  P_paire(ntype_max, ntype_max)
  SAVE  /SUB_PAIRE/

  ! mP_E0surR0_paire     = -P E0 / R0          (pour dérivée 1ère) ;
  ! P_Pp1_E0surR02_paire = P (P+1) E0 / R0**2  (pour dérivée 2ème).
  DOUBLE PRECISION mP_E0surR0_paire, P_Pp1_E0surR02_paire
  COMMON/SUB_PAIRE_AUX/ mP_E0surR0_paire(ntype_max, ntype_max)   &
       &              , P_Pp1_E0surR02_paire(ntype_max, ntype_max)
  SAVE  /SUB_PAIRE_AUX/

!.......................................................................
!----- Coupure et raccordement :
!-- "RACCORD_paire" et "COUPE_paire" sont les bornes de l'intervalle
  ! de coupure. "COUPE_2_paire" = COUPE_paire^2   [longueur] ;
  ! A_paire, BLJB_paire et C_paire sont les coefficients du
  ! polynome de raccordement.
  DOUBLE PRECISION RACCORD_paire, COUPE_paire, COUPE_2_paire   &
       &               , A_paire, B_paire, C_paire
  COMMON/SUB_PAIRE_CUT/ RACCORD_paire(ntype_max, ntype_max)   &
       &                , COUPE_paire(ntype_max, ntype_max)   &
       &              , COUPE_2_paire(ntype_max, ntype_max)   &
       &                   ,  A_paire(ntype_max, ntype_max)   &
       &                   ,  B_paire(ntype_max, ntype_max)   &
       &                   ,  C_paire(ntype_max, ntype_max)
  SAVE  /SUB_PAIRE_CUT/

 CONTAINS
!***********************************************************************
!----- Fonctions correspondant aux
!----- termes de base définissant les interactions, ainsi que leurs
!----- dérivées premières et secondes.
!-----
!----- 1) potentiel de paire sans coupure :
!-----       Phi_infini ; Phi_P_infini ; Phi_S_infini ;
!----- 2) potentiel de paire avec coupure douce ;
!-----       Phi, Phi_P, Phi_S ;
!-----------------------------------------------------------------------
!-----    Modifications:
!-----    ~~~~~~~~~~~~~~
!-----
!---------------
!-----  Auteur : Frédéric LANÇON --- 17 septembre 1997
!***********************************************************************
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Phi_infini( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- interaction de paire sans coupure (potentiel de paire nominal).
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Phi_infini
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION e0, r0, p   ! pour alléger les notations.
!=======================================================================
      e0 = E0_paire(type_k,type_l)
      r0 = R0_paire(type_k,type_l)
      p  =  P_paire(type_k,type_l)

      Phi_infini =  e0 * EXP( p * ( 1 - rkl / r0 ) )

!=======================================================================
   END FUNCTION Phi_infini

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Phi_P_infini( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- Dérivée première de Phi_infini (sans coupure).
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Phi_P_infini
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION e0_P, r0, p   ! pour alléger les notations.
!=======================================================================
      e0_P = mP_E0surR0_paire(type_k,type_l)
           ! mP_E0surR0_paire = -p e0 / r0
      r0 = R0_paire(type_k,type_l)
      p  =  P_paire(type_k,type_l)

      Phi_P_infini =  e0_P * EXP( p * ( 1 - rkl / r0 ) )

!=======================================================================
   END FUNCTION Phi_P_infini

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Phi_S_infini( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- Dérivée seconde de Phi_infini (sans coupure).
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Phi_S_infini
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION e0_P, r0, p   ! pour alléger les notations.
!=======================================================================
      e0_P = mP_E0surR0_paire(type_k,type_l)
           ! mP_E0surR0_paire = -p e0 / r0
      r0 = R0_paire(type_k,type_l)
      p  =  P_paire(type_k,type_l)

      Phi_S_infini =  (-p*  e0_P / r0) * EXP( p * ( 1 - rkl / r0 ) )
!=======================================================================
   END FUNCTION Phi_S_infini

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Phi( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- interaction de paire avec coupure "douce".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Phi
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!----- Fonction :
!     DOUBLE PRECISION Phi_infini    ! potentiel de paire nominal.
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION a, b, c, x
!=======================================================================
        IF( rkl .LT. RACCORD_paire(type_k,type_l) )   THEN
           ! potentiel de paire nominal :
           Phi = Phi_infini( rkl, type_k, type_l )

        ELSE IF( rkl .LT. COUPE_paire(type_k,type_l) )   THEN
           ! Raccordement par un polynome :
           a = A_paire(type_k,type_l)
           b = B_paire(type_k,type_l)
           c = C_paire(type_k,type_l)
           x = ( rkl - COUPE_paire(type_k,type_l) )
           Phi = x**3 *(a * rkl**2 + b * rkl + c)
        ELSE
           Phi = 0.D0
        END IF

!=======================================================================
   END FUNCTION Phi

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Phi_P( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!-----  Dérivée première de Phi avec coupure "douce".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Phi_P
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!----- Fonction :
!     DOUBLE PRECISION Phi_P_infini    ! potentiel de paire nominal.
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION a, b, c, x
!=======================================================================
        IF( rkl .LT. RACCORD_paire(type_k,type_l) )   THEN
           ! potentiel de paire nominal :
           Phi_P = Phi_P_infini( rkl, type_k, type_l )

        ELSE IF( rkl .LT. COUPE_paire(type_k,type_l) )   THEN
           ! Raccordement par un polynome :
           a = A_paire(type_k,type_l)
           b = B_paire(type_k,type_l)
           c = C_paire(type_k,type_l)
           x = ( rkl - COUPE_paire(type_k,type_l) )
           Phi_P =  x**2 *(                                    &
     &                      3.D0 *(a * rkl**2 + b * rkl + c)   &
     &                    +   x  *(2.D0 * a * rkl + b)         &
     &                     )
        ELSE
           Phi_P = 0.D0
        END IF

!=======================================================================
   END FUNCTION Phi_P
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Phi_S( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!-----  Dérivée seconde de Phi avec coupure "douce".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Phi_S
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!----- Fonction :
!     DOUBLE PRECISION Phi_S_infini    ! potentiel de paire nominal.
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION a, b, c, x
!=======================================================================
        IF( rkl .LT. RACCORD_paire(type_k,type_l) )   THEN

           Phi_S = Phi_S_infini( rkl, type_k, type_l )

        ELSE IF( rkl .LT. COUPE_paire(type_k,type_l) )   THEN
           ! Raccordement par un polynome :
           a = A_paire(type_k,type_l)
           b = B_paire(type_k,type_l)
           c = C_paire(type_k,type_l)
           x = ( rkl - COUPE_paire(type_k,type_l) )
           Phi_S =  x *(  6.D0 *(a * rkl**2 + b * rkl + c)   &
     &                + x * (6.D0  *(2.D0 * a * rkl + b)     &
     &                    + x * 2.D0 * a )                   &
     &                     )
        ELSE
           Phi_S = 0.D0
        END IF

!=======================================================================
   END FUNCTION Phi_S


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   SUBROUTINE DEBUT_Phi
!***********************************************************************
!***** Lecture des constantes et initialisations .
!***********************************************************************
      USE module_atomes_rgl
!...............................................................................
      IMPLICIT NONE
!.......................................................................
!---- Variables locales
   INTEGER :: i, j
   REAL(KIND=8) alpha, alpha3, alpha4, alpha5
   REAL(KIND=8) TRU
   REAL(KIND=8) fc0, fc1, fc2
!=======================================================================
!-----------------------------------------------------------------------
!--- Potentiel liaison forte au second moment pour Au (14 Février 1994) :
   P_paire(1, 1) = 10.3651922828969D0
   E0_paire(1, 1) = 0.22269500104212D0
   RACCORD_paire(1, 1) = 1.6585D0
   COUPE_paire(1, 1) = 1.980D0
!-----------------------------------------------------------------------
!----- Parametres du potentiel de paire:
!-----               ~~~~~~~~~~~~~~~~~~
      WRITE(*,*)
      WRITE(*,*)   &
           'Potentiel  de paire = E0* (r/r0)^(-p),'
      WRITE(*,*)   &
           '~~~~~~~~~~~~~~~~~~~ '
      WRITE(*,*)   &
           'raccorde a 0 par un polynome P entre T/r0 (P=f ) et C/r0 (P=0).'

!----- Exposant:
      WRITE(*,*)'Valeur de p pour chaque type de paire :'
!!$      READ(*,*)   ( P_paire( i, i) , i = 1 , ntype)   &
!!$     &         , (( P_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)  ( P_paire( i, i) , i = 1 , ntype)   &
     &         , (( P_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            P_paire(i,j) = P_paire(j,i)
         END DO
      END DO

!----- Coefficients multiplicatifs

      WRITE(*,*)'Entrez E0 pour chaque type de paire :'
!!$      READ(*,*)   ( E0_paire( i, i) , i = 1 , ntype)   &
!!$     &         , (( E0_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)  ( E0_paire( i, i) , i = 1 , ntype)   &
     &        , (( E0_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            E0_paire(i,j) = E0_paire(j,i)
         END DO
      END DO

      DO  j = 1, ntype
         DO  i = 1, ntype
            E0_paire(i,j) = 2.D0 * E0_paire(i,j)
            ! Le 2 vient du fait que le programme ne somme que sur la moitie
            ! des paires, alors que la formule utilisee somme sur toutes les
            ! paires atomiques.
         END DO
      END DO

      ! Grandeurs auxiliaires pour les dérivées :
      DO  j = 1, ntype
         DO  i = 1, ntype
           ! pour dérivée 1ère :
            mP_E0surR0_paire(i, j) = - P_paire(i, j)   &
     &                             *  E0_paire(i, j) / R0_paire(i, j)
           ! pour dérivée 2ème :
            P_Pp1_E0surR02_paire(i, j) =   P_paire(i, j)         &
     &                                 * ( P_paire(i, j) + 1 )   &
     &                                 *  E0_paire(i, j)         &
     &                                 /  R0_paire(i, j)**2
         END DO
      END DO

!----- Valeurs des raccordements et des coupures
      WRITE(*,*)'   Donnez T/r0 et C/r0 pour chaque type de paire.'
      WRITE(*,*)'Raccords en T_ij/r0 :'
!!$      READ(*,*)   &
!!$     &        ( RACCORD_paire( i, i) , i = 1 , ntype)   &
!!$     &     , (( RACCORD_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)   &
     &       ( RACCORD_paire( i, i) , i = 1 , ntype)   &
     &    , (( RACCORD_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            RACCORD_paire(i,j) = RACCORD_paire(j,i)
         END DO
      END DO

      DO i = 1, ntype
         DO j = 1, ntype
            RACCORD_paire(j,i) = RACCORD_paire(j,i) * R0_paire(j,i)
         END DO
      END DO

      WRITE(*,*)'Coupures en C_ij/r0 :'
!!$      READ(*,*)   &
!!$     &        ( COUPE_paire( i, i) , i = 1 , ntype)   &
!!$     &     , (( COUPE_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)   &
     &       ( COUPE_paire( i, i) , i = 1 , ntype)   &
     &    , (( COUPE_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            COUPE_paire(i,j) = COUPE_paire(j,i)
         END DO
      END DO

      DO i = 1, ntype
         DO j = 1, ntype
            COUPE_paire(j,i) = COUPE_paire(j,i) * R0_paire(j,i)
            COUPE_2_paire(j,i) = COUPE_paire(j,i)**2
         END DO
      END DO
!-----------------------------------------------------------------------
!----- Raccordement du potentiel de paire entre raccord et coupure :
!----- Calcul des coefficients du polynome:
      PRINT *
      WRITE(*,*)'Polynomes de raccordement des termes de paire:'
      DO 8 j = 1, ntype
      DO 9 i = 1, ntype
         alpha = RACCORD_paire(i,j) - COUPE_paire(i,j)
         !  alpha en Angstroem.
         alpha3 = alpha * alpha * alpha
         alpha4 = alpha * alpha3
         alpha5 = alpha * alpha4
         TRU = RACCORD_paire(i,j)
!         TRU en Angstroem.

         fc0 =    Phi_infini( RACCORD_paire(i,j) , i, j )
         fc1 =  Phi_P_infini( RACCORD_paire(i,j) , i, j )   ! dérivée 1ère.
         fc2 =  Phi_S_infini( RACCORD_paire(i,j) , i, j )   ! dérivée 2ème.
         A_paire(i,j) = 0.5D0 * fc2 / alpha3      &
                      -  3.D0 * fc1 / alpha4      &
                      +  6.D0 * fc0 / alpha5
         B_paire(i,j) =         fc1 / alpha3      &
                      -  3.D0 * fc0 / alpha4      &
                      -  2.D0 * A_paire(i,j) * TRU
         C_paire(i,j) =         fc0 / alpha3               &
                      -         A_paire(i,j) * TRU * TRU   &
                      -         B_paire(i,j) * TRU
         WRITE(*,400) i, j, COUPE_paire(i,j)   &
                    , A_paire(i,j), B_paire(i,j), C_paire(i,j)
9     END DO
8     END DO
400   FORMAT(   &
       'Polynome de raccordement (paire ',I1,',',I1,' ) :',/   &
     , 'P(X) = (X -',1PG13.6, ')**3 (', G13.6                  &
     , '* X**2 + ', G13.6  , '* X + ', G13.6  , ' )'  )

!-----------------------------------------------------------------------
!=======================================================================
    END SUBROUTINE DEBUT_Phi
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 END MODULE Interactions_de_paire_rgl



 MODULE Interactions_a_N_corps_rgl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
!.......................................................................

   USE module_atomes_rgl
!.......................................................................
   IMPLICIT NONE
!.......................................................................
!----- Paramètres pour interactions à N corps.

!--  Pour chaque type de paire :
   ! E0_Ncorps est le parametre multiplicatif du terme de paire [énergie] ;
   ! R0_Ncorps contientles distances de premier voisins [longueur] ;
   ! P_Ncorps est l'exposant pour une loi puissance ou exponentielle.
   DOUBLE PRECISION E0_Ncorps, R0_Ncorps, P_Ncorps
   COMMON/SUB_NCORPS/ E0_Ncorps(ntype_max, ntype_max)   &
       &           , R0_Ncorps(ntype_max, ntype_max)    &
       &           ,  P_Ncorps(ntype_max, ntype_max)
   SAVE  /SUB_NCORPS/

   ! mP_E0surR0_Ncorps     = -P E0 / R0         (pour dérivée 1ère) ;
   ! P_Pp1_E0surR02_Ncorps = P (P+1) E0 / R0**2 (pour dérivée 2ème).
   DOUBLE PRECISION mP_E0surR0_Ncorps, P_Pp1_E0surR02_Ncorps
   COMMON/SUB_NCORPS_AUX/ mP_E0surR0_Ncorps(ntype_max, ntype_max)   &
       &           , P_Pp1_E0surR02_Ncorps(ntype_max, ntype_max)
   SAVE  /SUB_NCORPS_AUX/

!.......................................................................
!----- Coupure et raccordement :
!-- "RACCORD_Ncorps" et "COUPE_Ncorps" sont les bornes de l'intervalle
   ! de coupure. "COUPE_2_Ncorps" = COUPE_Ncorps^2   [longueur] ;
   ! A_Ncorps, BLJB_Ncorps et C_Ncorps sont les coefficients du
   ! polynome de raccordement.
   DOUBLE PRECISION RACCORD_Ncorps, COUPE_Ncorps, COUPE_2_Ncorps   &
        &               , A_Ncorps, B_Ncorps, C_Ncorps
   COMMON/SUB_NCORPS_CUT/ RACCORD_Ncorps(ntype_max, ntype_max)   &
        &                 , COUPE_Ncorps(ntype_max, ntype_max)   &
        &               , COUPE_2_Ncorps(ntype_max, ntype_max)   &
        &                    ,  A_Ncorps(ntype_max, ntype_max)   &
        &                    ,  B_Ncorps(ntype_max, ntype_max)   &
        &                    ,  C_Ncorps(ntype_max, ntype_max)
   SAVE  /SUB_NCORPS_CUT/


 CONTAINS
!***********************************************************************
!----- Fonctions correspondant aux
!----- termes de base définissant les interactions, ainsi que leurs
!----- dérivées premières et secondes.
!-----
!----- 1) terme de paire sans coupure des interactions à N corps :
!-----       Upsilon_infini ; Upsilon_P_infini ; Upsilon_S_infini ;
!----- 2) "embedded function" F des interactions à N corps :
!-----       F ; F_P ; F_S ;
!----- 3) terme de paire avec coupure douce des interactions à N corps :
!-----       Upsilon ; Upsilon_P ; Upsilon_S.
!-----------------------------------------------------------------------
!----- Peut être utilisé avec des interactions à
!----- N corps nulles (modif. dans F_P et F_S)  (F. Lançon)
!-----------------------------------------------------------------------
!-----    Modifications:
!-----    ~~~~~~~~~~~~~~
!-----
!---------------
!-----  Auteur : Frédéric LANÇON --- 17 septembre 1997
!***********************************************************************
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Upsilon_infini( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- Potentiel à N corps : terme de paire sans coupure.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Upsilon_infini
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION e0, r0, p   ! pour alléger les notations.
!=======================================================================
      e0 = E0_Ncorps(type_k,type_l)
      r0 = R0_Ncorps(type_k,type_l)
      p  =  P_Ncorps(type_k,type_l)

      Upsilon_infini =  e0 * EXP( 2.D0 * p * ( 1.D0 - rkl / r0 ) )

!=======================================================================
   END FUNCTION Upsilon_infini

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Upsilon_P_infini( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- Potentiel à N corps : dérivée première de Upsilon_infini (sans coupure).
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Upsilon_P_infini
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION e0_P, r0, p   ! pour alléger les notations.
!=======================================================================
      e0_P = mP_E0surR0_Ncorps(type_k,type_l)
           ! mP_E0surR0_Ncorps = -p e0 / r0
      r0 = R0_Ncorps(type_k,type_l)
      p  =  P_Ncorps(type_k,type_l)

      Upsilon_P_infini =  2.D0 * e0_P   &
     &                 * EXP( 2.D0 * p * ( 1.D0 - rkl / r0 ) )

!=======================================================================
   END FUNCTION Upsilon_P_infini

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Upsilon_S_infini( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- Potentiel à N corps : dérivée seconde de Upsilon_infini (sans coupure).
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Upsilon_S_infini
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION e0_P, r0, p   ! pour alléger les notations.
!=======================================================================
      e0_P = mP_E0surR0_Ncorps(type_k,type_l)
           ! mP_E0surR0_Ncorps = -p e0 / r0
      r0 = R0_Ncorps(type_k,type_l)
      p  =  P_Ncorps(type_k,type_l)

      Upsilon_S_infini =  (- 4.D0 * e0_P * p / r0)   &
     &                 * EXP( 2.D0 * p * ( 1.D0 - rkl / r0 ) )
!=======================================================================
   END FUNCTION Upsilon_S_infini

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION F( rho_k )
!-----------------------------------------------------------------------
!----- Potentiel à N corps : "embedded function".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION F
!----- En entree :
      DOUBLE PRECISION rho_k    ! RHO(k).
!=======================================================================

      F =  - SQRT( rho_k )
!=======================================================================
   END FUNCTION F

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION F_P( rho_k )
!-----------------------------------------------------------------------
!----- Potentiel à N corps : dérivée première de l'"embedded function".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION F_P
!----- En entree :
      DOUBLE PRECISION rho_k    ! RHO(k).
!=======================================================================

      IF ( rho_k .NE. 0.D0 )   THEN
         F_P =  - 0.5D0 / SQRT( rho_k )
      ELSE
         ! L'atome k est isolé, ce facteur est infini et multiplié par la
         ! suite par zéro (car en facteur des interactions), alors... :
         F_P =  1.D0
      ENDIF

!=======================================================================
   END FUNCTION F_P

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION F_S( rho_k )
!-----------------------------------------------------------------------
!----- Potentiel à N corps : dérivée deuxième de l'"embedded function".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION F_S
!----- En entree :
      DOUBLE PRECISION rho_k    ! RHO(k).
!=======================================================================

      IF ( rho_k .NE. 0.D0 )   THEN
         F_S =   0.25D0 / SQRT( rho_k )**3
      ELSE
         ! L'atome k est isolé, ce facteur est infini et multiplié par la
         ! suite par zéro (car en facteur des interactions), alors... :
         F_S = - 1.D0
      ENDIF
!=======================================================================
   END FUNCTION F_S

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Upsilon( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!----- interaction à N corps : terme de paire avec coupure "douce".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Upsilon
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!----- Fonction :
!     DOUBLE PRECISION Upsilon_infini    ! fonction sans coupure.
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION a, b, c, x
!=======================================================================
        IF( rkl .LT. RACCORD_Ncorps(type_k,type_l) )   THEN
           ! terme de paire nominal :
           Upsilon = Upsilon_infini( rkl, type_k, type_l )

        ELSE IF( rkl .LT. COUPE_Ncorps(type_k,type_l) )   THEN
           ! Raccordement par un polynome :
           a = A_Ncorps(type_k,type_l)
           b = B_Ncorps(type_k,type_l)
           c = C_Ncorps(type_k,type_l)
           x = ( rkl - COUPE_Ncorps(type_k,type_l) )
           Upsilon =  x**3 *(a * rkl**2 + b * rkl + c)
        ELSE
           Upsilon = 0.D0
        END IF

!=======================================================================
   END FUNCTION Upsilon

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Upsilon_P( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!-----  Dérivée première de Upsilon avec coupure "douce".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Upsilon_P
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!----- Fonction :
!     DOUBLE PRECISION Upsilon_P_infini    ! potentiel de paire nominal.
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION a, b, c, x
!=======================================================================
        IF( rkl .LT. RACCORD_Ncorps(type_k,type_l) )   THEN
           ! terme de paire nominal :
           Upsilon_P = Upsilon_P_infini( rkl, type_k, type_l )

        ELSE IF( rkl .LT. COUPE_Ncorps(type_k,type_l) )   THEN
           ! Raccordement par un polynome :
           a = A_Ncorps(type_k,type_l)
           b = B_Ncorps(type_k,type_l)
           c = C_Ncorps(type_k,type_l)
           x = ( rkl - COUPE_Ncorps(type_k,type_l) )
           Upsilon_P =  x**2 *(                                    &
     &                          3.D0 *(a * rkl**2 + b * rkl + c)   &
     &                        +   x  *(2.D0 * a * rkl + b)         &
     &                         )
        ELSE
           Upsilon_P = 0.D0
        END IF

!=======================================================================
   END FUNCTION Upsilon_P
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   FUNCTION Upsilon_S( rkl, type_k, type_l )
!-----------------------------------------------------------------------
!-----  Dérivée seconde de Upsilon avec coupure "douce".
!-----------------------------------------------------------------------
      IMPLICIT NONE
!.......................................................................
!----- En sortie :
      DOUBLE PRECISION Upsilon_S
!----- En entree :
      DOUBLE PRECISION rkl      ! distance en Å entre 2 atomes (k et l).
      INTEGER type_k, type_l    ! types d'atome
!----- Fonction :
!     DOUBLE PRECISION Upsilon_S_infini    ! potentiel de paire nominal.
!.......................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION a, b, c, x
!=======================================================================
        IF( rkl .LT. RACCORD_Ncorps(type_k,type_l) )   THEN

           Upsilon_S = Upsilon_S_infini( rkl, type_k, type_l )

        ELSE IF( rkl .LT. COUPE_Ncorps(type_k,type_l) )   THEN
           ! Raccordement par un polynome :
           a = A_Ncorps(type_k,type_l)
           b = B_Ncorps(type_k,type_l)
           c = C_Ncorps(type_k,type_l)
           x = ( rkl - COUPE_Ncorps(type_k,type_l) )
           Upsilon_S =  x *(  6.D0 *(a * rkl**2 + b * rkl + c)   &
     &                    + x * (6.D0  *(2.D0 * a * rkl + b)     &
     &                        + x * 2.D0 * a )                   &
     &                         )
        ELSE
           Upsilon_S = 0.D0
        END IF

!=======================================================================
   END FUNCTION Upsilon_S

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   SUBROUTINE DEBUT_Upsilon
!***********************************************************************
!***** Lecture des constantes et initialisations .
!***********************************************************************
      USE module_atomes_rgl
!...............................................................................
      IMPLICIT NONE
!---- Variables locales
   INTEGER :: i, j
   REAL(KIND=8) TRU
   REAL(KIND=8) alp, alp3, alp4, alp5
   REAL(KIND=8) ro0, ro1, ro2
!=======================================================================
!-----------------------------------------------------------------------
!--- Potentiel liaison forte au second moment pour Au (14 Février 1994) :
   P_Ncorps(1, 1) = 4.2873410848603D0
   E0_Ncorps(1, 1) = 1.8521959188095D0
   RACCORD_Ncorps(1, 1) = 1.59073509D0
   COUPE_Ncorps(1, 1) = 1.980D0
!-----------------------------------------------------------------------
!----- Parametres du potentiel a N corps:
!-----               ~~~~~~~~~~~~~~~~~~~
      WRITE(*,*)
      WRITE(*,*) &
           'Potentiel a N corps = Beta**2 * (r/r0)^(-q)'
      WRITE(*,*) &
           '~~~~~~~~~~~~~~~~~~~ '

!----- Exposant:
      WRITE(*,*)'Valeur de q pour chaque type de paire i-j :'
!!$      READ(*,*)   ( P_Ncorps( i, i) , i = 1 , ntype)   &
!!$     &         , (( P_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)  ( P_Ncorps( i, i) , i = 1 , ntype)   &
     &         , (( P_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            P_Ncorps(i,j) = P_Ncorps(j,i)
         END DO
      END DO

!----- Coefficients multiplicatifs

      WRITE(*,*)'Entrez beta pour chaque type de paire i-j :'
!!$      READ(*,*)   ( E0_Ncorps( i, i) , i = 1 , ntype)   &
!!$     &         , (( E0_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)  ( E0_Ncorps( i, i) , i = 1 , ntype)   &
     &        , (( E0_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            E0_Ncorps(i,j) = E0_Ncorps(j,i)
         END DO
      END DO

      DO  j = 1, ntype
         DO  i = 1, ntype
            E0_Ncorps(i,j) = E0_Ncorps(i,j)**2
            ! Dans la formule, le beta est eleve au carre...
         END DO
      END DO

      ! Grandeurs auxiliaires pour les dérivées :
      DO  j = 1, ntype
         DO  i = 1, ntype
           ! pour dérivée 1ère :
            mP_E0surR0_Ncorps(i, j) = - P_Ncorps(i, j)   &
     &                              *  E0_Ncorps(i, j) / R0_Ncorps(i, j)
           ! pour dérivée 2ème :
            P_Pp1_E0surR02_Ncorps(i, j) =   P_Ncorps(i, j)       &
     &                                  * ( P_Ncorps(i, j) + 1 ) &
     &                                  *  E0_Ncorps(i, j)       &
     &                                  /  R0_Ncorps(i, j)**2
         END DO
      END DO

!----- Valeurs des raccordements et des coupures
      WRITE(*,*)' Donnez T/r0 et C/r0 pour chaque type de paire.'
      WRITE(*,*)'Raccords en T_ij/r0 :'
!!$      READ(*,*)   &
!!$     &      ( RACCORD_Ncorps( i, i) , i = 1 , ntype)   &
!!$     &   , (( RACCORD_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)   &
     &      ( RACCORD_Ncorps( i, i) , i = 1 , ntype)   &
     &   , (( RACCORD_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            RACCORD_Ncorps(i,j) = RACCORD_Ncorps(j,i)
         END DO
      END DO
      DO i = 1, ntype
         DO j = 1, ntype
            ! Passage d'unite "rayon atomique" en unité "[longueur]" (Å) :
            RACCORD_Ncorps(j,i) = RACCORD_Ncorps(j,i) * R0_Ncorps(j,i)
         END DO
      END DO

      WRITE(*,*)'Coupures en C_ij/r0 :'
!!$      READ(*,*)   &
!!$     &        ( COUPE_Ncorps( i, i) , i = 1 , ntype)   &
!!$     &     , (( COUPE_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)   &
     &       ( COUPE_Ncorps( i, i) , i = 1 , ntype)   &
     &    , (( COUPE_Ncorps( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      DO i = 2, ntype
         DO j = 1, i-1
            COUPE_Ncorps(i,j) = COUPE_Ncorps(j,i)
         END DO
      END DO
      DO i = 1, ntype
         DO j = 1, ntype
            ! Passage d'unite "rayon atomique" en unité "[longueur]" (Å) :
            COUPE_Ncorps(j,i) = COUPE_Ncorps(j,i) * R0_Ncorps(j,i)
            COUPE_2_Ncorps(j,i) = COUPE_Ncorps(j,i)**2
         END DO
      END DO
!-----------------------------------------------------------------------
!----- Polynome de raccordement de la fonction rho
!----- Calcul des coefficients A_Ncorps, B_Ncorps et C_Ncorps
!----- P(X)=(X-RCUT)**3 x (A_Ncorps x X**2 + B_Ncorps x X + C_Ncorps)

      WRITE(*,*)'Polynomes de raccordement des termes a N corps :'
      DO 10 j = 1, ntype
         DO 11 i = 1, ntype
            alp = RACCORD_Ncorps(i,j) - COUPE_Ncorps(i,j)
            ! alp en Angstroem.
            alp3 = alp*alp*alp
            alp4 = alp3*alp
            alp5 = alp4*alp
            TRU = RACCORD_Ncorps(i,j)
            ! TRU en Angstroem.

            ! ro0, ro1 et ro2 sont les valeurs en RACCORD_Ncorps
            ! de Upsilon et de ses derivees :
            ro0 =   Upsilon_infini( RACCORD_Ncorps(i,j) , i, j )   ! fc
            ro1 = Upsilon_P_infini( RACCORD_Ncorps(i,j) , i, j )   ! fc'
            ro2 = Upsilon_S_infini( RACCORD_Ncorps(i,j) , i, j )   ! fc"
            A_Ncorps(i,j) = 0.5D0*ro2/alp3  &
                          - 3.D0*ro1/alp4 + 6.D0*ro0/alp5
            B_Ncorps(i,j) = ro1/alp3 - 3.D0*ro0/alp4   &
                          - 2.D0*TRU*A_Ncorps(i,j)
            C_Ncorps(i,j) = ro0/alp3   &
                          - A_Ncorps(i,j)*TRU*TRU - B_Ncorps(i,j)*TRU
           WRITE(*,400) i, j, COUPE_Ncorps(i,j)   &
                      , A_Ncorps(i,j), B_Ncorps(i,j), C_Ncorps(i,j)
11       END DO
10    END DO
      PRINT *

400   FORMAT(   &
       'Polynome de raccordement (paire ',I1,',',I1,' ) :',/   &
     , 'P(X) = (X -',1PG13.6, ')**3 (', G13.6                  &
     , '* X**2 + ', G13.6  , '* X + ', G13.6  , ' )'  )
!=======================================================================
   END SUBROUTINE DEBUT_Upsilon

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 END MODULE Interactions_a_N_corps_rgl



 MODULE module_interactions_rgl
!~~~~~~~~~~~~~~~~~~~~~~~~~~
!-----------------------------------------------------------------------
!.......................................................................

   USE Interactions_de_paire_rgl
   USE Interactions_a_N_corps_rgl
!.......................................................................
   IMPLICIT NONE
!.......................................................................

   REAL(KIND=8), SAVE :: cutmax, cutmax2
   ! cutmax est la porté maximum des interactions (en angstrom).
!.......................................................................
!=======================================================================
 CONTAINS


   SUBROUTINE DEBUT_Interactions
!***********************************************************************
!***** Lecture des constantes et initialisations .
!***********************************************************************
!.......................................................................

   USE module_atomes_rgl
!...............................................................................
   IMPLICIT NONE
!...............................................................................
!---- Variables locales
   INTEGER :: i, j
!=======================================================================
!-----------------------------------------------------------------------
!--- Potentiel liaison forte au second moment pour Au (14 Février 1994) :
   ntype = 1                 ! nombre de types d'atomes.
   ATOME_NOM(1) = "Au"
   R0_paire(1, 1) = 2.874D0
!-----------------------------------------------------------------------
!----- Entrée de la distance entre premier voisin :

      WRITE(*,*)
      WRITE(*,*) 'Taille "atomique" r0 pour chaque type de paire :'
!!$      WRITE(*,*) ' Donnez : '   &
!!$     &         ,  ('('//   &
!!$     &             TRIM( ATOME_NOM(i) ) //'-'// TRIM( ATOME_NOM(i) )   &
!!$     &             //') '   &
!!$     &             , i = 1, ntype)   &
!!$     &         , (('('//   &
!!$     &             TRIM( ATOME_NOM(i) ) //'-'// TRIM( ATOME_NOM(j) )   &
!!$     &             //') '   &
!!$     &             , j = i+1, ntype)   &
!!$     &             , i = 1, ntype)
!!$      READ(*,*)   ( R0_paire( i, i) , i = 1 , ntype)   &
!!$     &         , (( R0_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
      WRITE(*,*)  ( R0_paire( i, i) , i = 1 , ntype)   &
     &        , (( R0_paire( i, j) , j = i+1 , ntype) , i = 1 , ntype)
         DO i = 2, ntype
            DO j = 1, i-1
               R0_paire(i,j) = R0_paire(j,i)
            END DO
         END DO
         DO j = 1, ntype
            DO i = 1, ntype
               R0_Ncorps(i,j) = R0_paire(i,j)
            END DO
         END DO
!-----------------------------------------------------------------------

      CALL DEBUT_Phi
      CALL DEBUT_Upsilon

!-----------------------------------------------------------------------
!-----  Recapitulation des coupures:
      PRINT *
      WRITE(*,*)'Recapitulation des raccords en Angstroem et relatif :'
      WRITE(*,*)'   --termes  de paire--'
      DO i = 1 , ntype
         DO j = i , ntype
           WRITE(*,*) '   '//ATOME_NOM(i)//'- '//ATOME_NOM(j)   &
     &              , RACCORD_paire(i,j), ' Å ; '               &
     &              , RACCORD_paire(i,j)/R0_paire(i,j), ' r0'
         END DO
      END DO

      PRINT *
      WRITE(*,*)'   --termes a N corps--'
      DO i = 1 , ntype
         DO j = i , ntype
           WRITE(*,*) '   '//ATOME_NOM(i)//'- '//ATOME_NOM(j)   &
     &              , RACCORD_Ncorps(i,j), ' Å ; '              &
     &              , RACCORD_Ncorps(i,j)/R0_Ncorps(i,j), ' r0'
         END DO
      END DO

      PRINT *
      WRITE(*,*)'Recapitulation des coupures en Angstroem et relatif :'
      WRITE(*,*)'   --termes  de paire--'
      DO i = 1 , ntype
         DO j = i , ntype
           WRITE(*,*) '   '//ATOME_NOM(i)//'- '//ATOME_NOM(j)   &
     &              , COUPE_paire(i,j), ' Å ; '                 &
     &              , COUPE_paire(i,j)/R0_paire(i,j), ' r0'
         END DO
      END DO

      PRINT *
      WRITE(*,*)'   --termes a N corps--'
      DO i = 1 , ntype
         DO j = i , ntype
           WRITE(*,*) '   '//ATOME_NOM(i)//'- '//ATOME_NOM(j)   &
     &              , COUPE_Ncorps(i,j), ' Å ; '                &
     &              , COUPE_Ncorps(i,j)/R0_Ncorps(i,j), ' r0'
         END DO
      END DO
      PRINT *

!-----------------------------------------------------------------------
!  Determination du rayon de coupure maximum.

      cutmax = MAX( MAXVAL( COUPE_paire(1:ntype,1:ntype))   &
                  , MAXVAL(COUPE_Ncorps(1:ntype,1:ntype)) )
      cutmax2 = cutmax**2
!=======================================================================
   END SUBROUTINE DEBUT_Interactions

!=======================================================================
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
END MODULE module_interactions_rgl

module module_RosatoGuillopeLegrand_Au

logical, save :: initialized

contains

!=======================================================================!
!-----------------                                        --------------!
!=======================================================================!
	 SUBROUTINE ENERGETIQUE( nat, n_max_de_voisins, voisinage, dxx, dyy, dzz, XYZ, FORCE, energie )

!------------------------------------------------------------------------------
!----- Calcule:  Calcule l'ï¿½nergie et diffï¿½rentes dï¿½rivï¿½es.
!------------------------------------------------------------------------------
!-----     Modifications:
!-----    ~~~~~~~~~~~~~~
!--- 23/ 4/2001 : sous-programme mis sous une forme indï¿½pendante pour Stefan
!---       Goedecker. [Frï¿½dï¿½ric LANï¿½ON]
!--- 29/ 9/1997 : passage en Fortran 90 (Frï¿½dï¿½ric LANï¿½ON).
!---------------
!-----  Auteur : Frï¿½dï¿½ric LANï¿½ON --- 25 juin 1997
!------------------------------------------------------------------------------
!----- Modules :

   USE module_atomes_rgl

   USE module_interactions_rgl
   ! interactions de paire,
   ! et fonction ï¿½ N corps.

!...............................................................................
   IMPLICIT NONE
!...............................................................................
   INTEGER, INTENT(IN)                                     :: nat             ! nombre total d'atomes.
   INTEGER, INTENT(IN)                                     :: n_max_de_voisins! nombre maximum de voisins par atome.
   INTEGER, DIMENSION(0:n_max_de_voisins, nat), INTENT(IN) :: voisinage
   REAL(KIND=8), INTENT(IN)                     :: dxx, dyy, dzz     ! boï¿½te pï¿½riodique orthogonale.
   REAL(KIND=8), DIMENSION(3, nat), INTENT(IN)  :: XYZ               ! positions atomiques.
   REAL(KIND=8), INTENT(OUT)                    :: energie
   REAL(KIND=8), DIMENSION(3, nat), INTENT(OUT) :: FORCE
!...............................................................................
!--- Inspirï¿½ de MODULE module_atomes_rgl
!               ~~~~~~~~~~~~~~~~~~~

!   INTEGER, DIMENSION(nat) :: ITYPE ! alliage du pauvre -> 1 sorte d'atomes
   INTEGER, DIMENSION(140000) :: ITYPE ! alliage du pauvre -> 1 sorte d'atomes
   ! types des NAT atomes.
!--- Fin de (Inspirï¿½ de END MODULE module_atomes_rgl)
!..............................................................................
!...............................................................................
!--- Inspirï¿½ de MODULE module_energetique
!               ~~~~~~~~~~~~~~~~~~~~~~~~~
   REAL(KIND=8) :: energie_paire, energie_Ncorps
!!$   REAL(KIND=8), DIMENSION(3, 3) :: CONTRAINTES
!!$   REAL(KIND=8), SAVE :: force_max, force_moy

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: RHO
   ! RHO(i) : "densitï¿½ ï¿½lectronique effective"

   INTEGER  :: n_paires         ! nombre de paires en interaction...
   INTEGER, ALLOCATABLE, DIMENSION(:) :: PAIRE_1, PAIRE_2
   ! ï¿½PAIRE_1(i)--PAIRE_2(i)ï¿½ : i-ï¿½me paire d'atomes en interaction.

!--  R_PAIRE(i)   : distance en ï¿½ entre les atomes PAIRE_1(i) et PAIRE_2(i)
   ! U_PAIRE(:,i) : vecteur unitaire (Xlk, Ylk, Zlk) avec
   !            xlk = [X(k)-X(l)]/Rkl . Attention : Xlk = - Xkl !
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: R_PAIRE
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: U_PAIRE

!--- Fin de (Inspirï¿½ de END MODULE module_energetique)
!..............................................................................
!----- Variables locales :
!      ~~~~~~~~~~~~~~~~~~~
   INTEGER  :: n_paires_max  ! nombre maximum de paires en interaction...
   INTEGER k, l, i, first_i, last_i, type_k, type_l, alpha, beta
   REAL(KIND=8) :: Upsilon_kl
   REAL(KIND=8) :: dE_sur_drkl
!!$   REAL(KIND=8) :: volume

   REAL(KIND=8) :: xk, yk, zk
   ! Coordonnees de l'atome k.
   REAL(KIND=8) :: xkl, ykl, zkl, rkl, rkl2
   ! Grandeurs entre les atomes k et l.

	allocate(U_PAIRE(3,nat*n_max_de_voisins),  R_PAIRE(nat*n_max_de_voisins),  RHO(nat*n_max_de_voisins),  &
                 PAIRE_1(nat*n_max_de_voisins),  PAIRE_2(nat*n_max_de_voisins))
!	write(66,*) 'Mbytes of alloction: ',8*6*nat*n_max_de_voisins*1.d-6
!..............................................................................
!==============================================================================
!------------------------------------------------------------------------------
   ITYPE = 1 ! alliage du pauvre -> 1 sorte d'atomes
!------------------------------------------------------------------------------
!----- Calcul des paires rï¿½ellement en interaction :

!   write(*,*) 'cutmax2=',cutmax2

   n_paires_max = nat * n_max_de_voisins
   n_paires = 0

   DO  k = 1, nat
   ! parcours de la liste des atomes
      type_k = ITYPE (k)

      xk = XYZ(1, k)
      yk = XYZ(2, k)
      zk = XYZ(3, k)

      DO  i = 1, voisinage(0, k)
      ! Parcours de la liste des voisins l de l'atome k
         l = voisinage( i, k )
         IF( l < k )   CYCLE

         xkl = PLUS_PROCHE( XYZ(1, l) - xk, dxx )
         rkl2 = xkl**2
         IF(rkl2 > cutmax2)   CYCLE
         ykl = PLUS_PROCHE( XYZ(2, l) - yk, dyy )
         rkl2 = rkl2 + ykl**2
         IF(rkl2 > cutmax2)   CYCLE
         zkl = PLUS_PROCHE( XYZ(3, l) - zk, dzz )
         rkl2 = rkl2 + zkl**2
         IF(rkl2 > cutmax2)   CYCLE
         ! il ne reste plus que ce qui est dans la sphï¿½re de voisinage :

         n_paires = n_paires + 1
         IF( n_paires > n_paires_max )  &
              STOP "ERREUR sur les nombres de voisins :-("
         PAIRE_1(n_paires) = k
         PAIRE_2(n_paires) = l
         R_PAIRE(n_paires) = SQRT( rkl2 )     ! en ï¿½
         U_PAIRE(1, n_paires) = xkl / R_PAIRE(n_paires)
         U_PAIRE(2, n_paires) = ykl / R_PAIRE(n_paires)
         U_PAIRE(3, n_paires) = zkl / R_PAIRE(n_paires)
      END DO

      ! on va passer ï¿½ l'atome numï¿½ro k+1.
   END DO


!      WRITE(*,*)        n_paires, " paires d'atomes rï¿½ellement en interaction."
!!$   WRITE(unit_com,*) n_paires, " paires d'atomes rï¿½ellement en interaction."
!------------------------------------------------------------------------------
!----- Calcule :  "the local charge density" RHO(i)
!----- RHO(i)  :   densitï¿½ ï¿½lectronique effective.

   RHO = 0.D0

   DO  i = 1, n_paires
   ! parcours de la liste des paires d'atomes en interaction.
      k = PAIRE_1(i)
      l = PAIRE_2(i)

      Upsilon_kl = Upsilon (R_PAIRE (i), ITYPE(k), ITYPE(l) )
      RHO (k) = RHO (k) + Upsilon_kl
      RHO (l) = RHO (l) + Upsilon_kl

   END DO

!-----------------------------------------------------------------------
!----- Calcul de l'ï¿½nergie :

!-- Terme ï¿½ N corps :

   energie_Ncorps = 0.D0
   DO  k = 1, nat
      energie_Ncorps = energie_Ncorps + F( RHO(k) )
   END DO

!-- Terme de paire :

   energie_paire = 0.D0
   DO  i = 1, n_paires
      ! parcours de la liste des paires d'atomes en interaction.
      k = PAIRE_1(i)
      l = PAIRE_2(i)

      energie_paire = energie_paire + Phi( R_PAIRE(i), ITYPE(k), ITYPE(l) )
!   write(*,*)  Phi( R_PAIRE(i), ITYPE(k), ITYPE(l) )

   END DO

   energie_paire = energie_paire
   energie_Ncorps = energie_Ncorps
   energie = energie_paire + energie_Ncorps

!!$   WRITE ( unit_com,* )
!!$   WRITE ( unit_com,* ) 'ï¿½nergie par atome : '
!!$   WRITE ( unit_com,* ) 'contribution de paire  : ', energie_paire
!!$   WRITE ( unit_com,* ) 'contribution ï¿½ N corps : ', energie_Ncorps
!!$   WRITE ( unit_com,* ) '    ï¿½nergie totale     : ', energie

!-----------------------------------------------------------------------
!----- Calcul des forces ( = -grad(E) ) et des contraintes :

   FORCE = 0.D0
!!$   CONTRAINTES = 0.D0

   DO  i = 1, n_paires
   ! parcours de la liste des paires d'atomes en interaction.
      rkl = R_PAIRE(i)
      k = PAIRE_1(i)
      l = PAIRE_2(i)
      type_k = ITYPE (k)
      type_l = ITYPE (l)

      dE_sur_drkl =     Phi_P( rkl, type_k, type_l )+ Upsilon_P( rkl, type_k, type_l ) * ( F_P(RHO(k)) + F_P(RHO(l)) )
                            ! >-pair                         !\                                    ! >- N corps

      ! les forces, c'est l'opposï¿½ du gradiant :
      FORCE (:, k) = FORCE (:, k) + dE_sur_drkl * U_PAIRE(:, i)
      FORCE (:, l) = FORCE (:, l) - dE_sur_drkl * U_PAIRE(:, i)
                                ! ^
                                ! xlk = - xkl = - U_LISTE(1, i)

!!$      DO beta = 1, 3 ;  DO alpha = 1, 3
!!$            ! calcul des contraintes :
!!$            CONTRAINTES (alpha, beta) = CONTRAINTES (alpha, beta)    &
!!$                      + dE_sur_drkl * rkl * U_PAIRE (alpha, i)       &
!!$                                          * U_PAIRE (beta , i)
!!$      ENDDO ;  ENDDO

   END DO


!!$   force_max = MAXVAL(ABS(FORCE))
!!$   force_moy = SQRT ( SUM(FORCE**2) / nat )

!!$   WRITE ( unit_com,* )
!!$   WRITE ( unit_com,* ) 'Force sur les atomes : '
!!$   WRITE ( unit_com,* ) '  composante maximale  : ', force_max
!!$   WRITE ( unit_com,* ) 'force (module) moyenne : ', force_moy

!!$   volume = dxx * dyy * dzz
!!$   CONTRAINTES = CONTRAINTES / volume

!!$   WRITE ( unit_com,* )
!!$   WRITE ( unit_com,* ) 'Contrainte sigma : '
!!$   WRITE ( unit_com &
!!$        &, " (' /', 3 (1PG15.6) , '  \',  / , '| ', 3 (1PG15.6)       &
!!$        &, '  |',  / , '\ ', 3 (1PG15.6) , ' /') ")  ( (CONTRAINTES (alpha,&
!!$        & beta) , alpha = 1, 3) , beta = 1, 3)
!!$   WRITE ( unit_com,* ) 'Pression : ', - (CONTRAINTES (1, 1) +            &
!!$      CONTRAINTES (2, 2) + CONTRAINTES (3, 3) )

!------------------------------------------------------------------------------
!==============================================================================
      deallocate(U_PAIRE,R_PAIRE,RHO,PAIRE_1,PAIRE_2)
 CONTAINS

   !===========================================================================
   FUNCTION PLUS_PROCHE( x, divise )
   ! divise doit ï¿½tre positif (n'est pas testï¿½).
   !---------------------------------------------------------------------------
   !----- Pour calculer les diffï¿½rences de coordonnï¿½es entre atomes en
   !----- prenant leurs images pï¿½riodiques les plus proches.
   !---------------------------------------------------------------------------
   !----- Dï¿½but de gï¿½nï¿½rique le 25/ 9/97 (Frï¿½dï¿½ric Lanï¿½on).
   !----- Auteur : Frï¿½dï¿½ric Lanï¿½on, le 27 juin 1997
   !---------------------------------------------------------------------------
     IMPLICIT NONE
   !...........................................................................
     REAL(KIND=8), INTENT(IN) :: x
     REAL(KIND=8), INTENT(IN), OPTIONAL :: divise
     REAL(KIND=8) :: PLUS_PROCHE
   !===========================================================================
     IF( PRESENT(divise) )   THEN
        PLUS_PROCHE = x - NINT( x / divise ) * divise
     ELSE   ! on est en coordonnï¿½e rï¿½duite (divise = 1 )
        PLUS_PROCHE = x - NINT( x )
     END IF
   !===========================================================================
   END FUNCTION PLUS_PROCHE
   !===========================================================================
!==============================================================================


 END SUBROUTINE ENERGETIQUE

        subroutine verlet_q(cut,nat,nnbrx,alat,pos,lst)
        implicit real*8 (a-h,o-z)
        dimension pos(3,nat),lst(0:nnbrx,nat),alat(3)

!	write(*,*) 'nat ',nat
!	write(*,*) 'cut ',cut
!	write(*,*) 'nnbrx ',nnbrx
!	write(*,*) 'alat ',alat
!	write(*,*) 'pos ',pos
!	write(*,*) 'INPUT LST'
!	do iat=1,nat
!	  write(*,'(i4,35(i3))') (lst(ibrx,iat), ibrx=0,nnbrx)
!	enddo

        cut2=cut**2

        do 63 i=1,nat
             lst0=0
           do 66 j=1,nat
             if (j.ne.i) then
               xrel1= pos(1,j)-pos(1,i)
               if (abs(xrel1-alat(1)).lt.abs(xrel1)) xrel1=xrel1-alat(1)
               if (abs(xrel1+alat(1)).lt.abs(xrel1)) xrel1=xrel1+alat(1)
                 xrel2= pos(2,j)-pos(2,i)
               if (abs(xrel2-alat(2)).lt.abs(xrel2)) xrel2=xrel2-alat(2)
               if (abs(xrel2+alat(2)).lt.abs(xrel2)) xrel2=xrel2+alat(2)
                 xrel3= pos(3,j)-pos(3,i)
               if (abs(xrel3-alat(3)).lt.abs(xrel3)) xrel3=xrel3-alat(3)
               if (abs(xrel3+alat(3)).lt.abs(xrel3)) xrel3=xrel3+alat(3)
                 rr2=xrel1**2 + xrel2**2 + xrel3**2
               if ( rr2 .le. cut2 ) then
                 lst0=lst0+1
               if (lst0.gt.nnbrx) stop 'enlarge nnbrx'
                 lst(lst0,i)=j
               endif
            endif
66      continue
        lst(0,i)=lst0
!         print*, 'lst0',lst0
63	continue

!do iat=1,nat
!write(*,'(i4,35(i3))') (lst(ibrx,iat), ibrx=0,nnbrx)
!enddo


        return
        end
!----------------------------------------------------------------------------!
	subroutine energyandforces_rgl_au(nat,alat, rxyz,fxyz,etot)
	use module_interactions_rgl
	implicit none
	integer :: nat
	real*8  :: rxyz(3,nat), alat(3), etot,fxyz(3,nat)
	integer :: n_max_de_voisins
	integer, allocatable :: voisinage( : , : )
	real*8  :: cut
        if(.not.initialized)then                                                                                                        
            call f_err_throw('Potential "lensic" not initialized',&                                                                     
                 err_name='BIGDFT_RUNTIME_ERROR')                                                                                       
        endif 
	n_max_de_voisins=50
	ALLOCATE(voisinage(0:n_max_de_voisins, nat))
        voisinage=0
	cut=1.7d0*2.874d0
!	write(*,*) 'rxyz ',rxyz
	call verlet_q(cut,nat,n_max_de_voisins,alat,rxyz,voisinage)
!	write(*,*) 'voisinage ',voisinage
	CALL ENERGETIQUE( nat, n_max_de_voisins, voisinage,alat(1), alat(2), alat(3), rxyz, fxyz, etot )
	deallocate(voisinage)
	return
	end
subroutine init_RosatoGuillopeLegrand_Au(nat,astruct,paramset,paramfile,geocode)
    use module_base
    use yaml_output
    use module_atoms
    use module_interactions_rgl
    implicit none
    !parameter
    integer, intent(in) :: nat
    type(atomic_structure), intent(in) :: astruct
    character(len=*), intent(in) :: paramset
    character(len=*), intent(in) :: paramfile
    character(len=*), intent(in) :: geocode
    initialized=.false. 
    call yaml_comment('Initializing RGLAU',hfill='-')
    call yaml_mapping_open('Using Au parameters from'//&
         ' T Deutsch et al 1995 J. Phys.: Condens. Matter 7 6407. doi:10.1088/0953-8984/7/32/007')
    !list parameters
    !!call yaml_map('D(Cl<->Cl) [Bohr^8 Hartree]', parameters(4,3),  fmt='(1pe10.4)')
    call yaml_mapping_close()
    if(trim(paramfile)/='none')then
        inquire(file=trim(adjustl(paramfile)),exist=exists)                                                                         
        if(.not.exists)then                                                                                                         
            call f_err_throw('Parameter file '//trim(adjustl(paramfile))//&                                                         
                 ' not found.')                                                                                                     
        endif                                                                                                                       
        u=f_get_free_unit()                                                                                                         
        open(unit=u,file=trim(adjustl(paramfile)))                                                                                  
            read(u,*) alat_int(1), alat_int(2), alat_int(3)                                                                         
        close(u)                                                                                                                    
        if (units=='angstroem' .or. units=='angstroemd0') then                                                                      
           ! if Angstroem convert to Bohr                                                                                           
           alat_int(1)=alat_int(1)/Bohr_Ang                                                                                         
           alat_int(2)=alat_int(2)/Bohr_Ang                                                                                         
           alat_int(3)=alat_int(3)/Bohr_Ang                                                                                         
        else if  (units=='atomic' .or. units=='bohr'  .or.&                                                                         
             units== 'atomicd0' .or. units== 'bohrd0') then                                                                         
             !do nothing                                                                                                            
        else if (units == 'reduced') then                                                                                           
           !assume that for reduced coordinates cell size is in bohr                                                                
        else                                                                                                                        
           call f_err_throw('Length units in input file unrecognized.' // &                                                         
                'Recognized units are angstroem or atomic = bohr',err_id=BIGDFT_INPUT_VARIABLES_ERROR)                              
           return                                                                                                                   
           !write(*,*) 'length units in input file unrecognized'                                                                    
           !write(*,*) 'recognized units are angstroem or atomic = bohr'                                                            
           !stop                                                                                                                    
        endif                                                                                                                       
        call yaml_map('Read cell parameters from file',trim(adjustl(paramfile)))                                                    
        initialized=.true. 
    else
        select case(trim(paramset))
        case('default')
            call f_err_throw('No "default" parameter set for RosatoGuillopeLegrand_Au defined.')
        case default
            call f_err_throw('Following parameter set for RosatoGuillopeLegrand_Au force field '//&
                'is unknown: '//trim(adjustl(paramset)))
        end select
    endif
    call DEBUT_Phi
    call DEBUT_Upsilon
    call DEBUT_Interactions
end subroutine
end module module_RosatoGuillopeLegrand_Au
