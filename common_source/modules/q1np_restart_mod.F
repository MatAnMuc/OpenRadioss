Chd|====================================================================
Chd|  Q1NP_RESTART_MOD              common_source/modules/q1np_restart_mod.F
Chd|====================================================================
      MODULE Q1NP_RESTART_MOD
      IMPLICIT NONE
C-----------------------------------------------
C   m y _ r e a l
C-----------------------------------------------
#include "my_real.inc"
C-----------------------------------------------
      INTEGER :: NUMELQ1NP_G      = 0
      INTEGER :: SKQ1NP_G         = 0
      INTEGER :: SIQ1NP_G         = 0
      INTEGER :: SQ1NPBULK_G      = 0
      INTEGER :: SQ1NPCTRL_SHARED_G = 0
      INTEGER :: SQ1NPCTRL_L_G    = 0
      INTEGER :: SQ1NPWEIGHT_L_G  = 0
      INTEGER :: SQ1NPKNOT_L_G    = 0
      INTEGER :: TABVINT_LEN_G    = 0
      INTEGER :: NKQ1NP_G         = 0
      LOGICAL, ALLOCATABLE :: IS_Q1NP_REPLACED_G(:)
C
C     Q1NP Gauss quadrature scheme variables
C
      INTEGER :: Q1NP_NP_U_G = 0  ! Number of Gauss points in U direction
      INTEGER :: Q1NP_NP_V_G = 0  ! Number of Gauss points in V direction
      INTEGER :: Q1NP_NP_T_G = 0  ! Number of Gauss points in T direction
      my_real, DIMENSION(:), ALLOCATABLE :: Q1NP_GP_U_G  ! Gauss point coordinates in U direction
      my_real, DIMENSION(:), ALLOCATABLE :: Q1NP_GP_V_G  ! Gauss point coordinates in V direction
      my_real, DIMENSION(:), ALLOCATABLE :: Q1NP_GP_T_G  ! Gauss point coordinates in T direction
      my_real, DIMENSION(:), ALLOCATABLE :: Q1NP_GW_U_G  ! Gauss weights in U direction
      my_real, DIMENSION(:), ALLOCATABLE :: Q1NP_GW_V_G  ! Gauss weights in V direction
      my_real, DIMENSION(:), ALLOCATABLE :: Q1NP_GW_T_G  ! Gauss weights in T direction
C
C     Q1NP grid dimensions and knot vectors
C
      INTEGER :: Q1NP_NX_G = 0  ! Number of elements in U direction
      INTEGER :: Q1NP_NY_G = 0  ! Number of elements in V direction
      my_real, DIMENSION(:), POINTER :: Q1NP_KTAB_G => NULL()  ! Knot vectors (U then V)

      CONTAINS

      SUBROUTINE RESET_Q1NP_COUNTS()
        NUMELQ1NP_G      = 0
        SKQ1NP_G         = 0
        SIQ1NP_G         = 0
        SQ1NPBULK_G      = 0
        SQ1NPCTRL_SHARED_G = 0
        SQ1NPCTRL_L_G    = 0
        SQ1NPWEIGHT_L_G  = 0
        SQ1NPKNOT_L_G    = 0
        TABVINT_LEN_G    = 0
        NKQ1NP_G         = 0
        Q1NP_NP_U_G      = 0
        Q1NP_NP_V_G      = 0
        Q1NP_NP_T_G      = 0
        Q1NP_NX_G        = 0
        Q1NP_NY_G        = 0
        Q1NP_KTAB_G      => NULL()
        IF (ALLOCATED(Q1NP_GP_U_G)) DEALLOCATE(Q1NP_GP_U_G)
        IF (ALLOCATED(Q1NP_GP_V_G)) DEALLOCATE(Q1NP_GP_V_G)
        IF (ALLOCATED(Q1NP_GP_T_G)) DEALLOCATE(Q1NP_GP_T_G)
        IF (ALLOCATED(Q1NP_GW_U_G)) DEALLOCATE(Q1NP_GW_U_G)
        IF (ALLOCATED(Q1NP_GW_V_G)) DEALLOCATE(Q1NP_GW_V_G)
        IF (ALLOCATED(Q1NP_GW_T_G)) DEALLOCATE(Q1NP_GW_T_G)
      END SUBROUTINE RESET_Q1NP_COUNTS

      SUBROUTINE RESET_Q1NP_STATE()
        CALL RESET_Q1NP_COUNTS()
      END SUBROUTINE RESET_Q1NP_STATE

      SUBROUTINE SET_Q1NP_COUNTS(NUMELQ1NP_IN,
     .                           SKQ1NP,SIQ1NP,SQ1NPBULK,
     .                           SQ1NPCTRL_SHARED,SQ1NPCTRL_L,
     .                           SQ1NPWEIGHT_L,SQ1NPKNOT_L)
        INTEGER, INTENT(IN) :: NUMELQ1NP_IN
        INTEGER, INTENT(IN) :: SKQ1NP,SIQ1NP,SQ1NPBULK
        INTEGER, INTENT(IN) :: SQ1NPCTRL_SHARED,SQ1NPCTRL_L
        INTEGER, INTENT(IN) :: SQ1NPWEIGHT_L,SQ1NPKNOT_L

        NUMELQ1NP_G      = NUMELQ1NP_IN
        SKQ1NP_G         = SKQ1NP
        SIQ1NP_G         = SIQ1NP
        SQ1NPBULK_G      = SQ1NPBULK
        SQ1NPCTRL_SHARED_G = SQ1NPCTRL_SHARED
        SQ1NPCTRL_L_G    = SQ1NPCTRL_L
        SQ1NPWEIGHT_L_G  = SQ1NPWEIGHT_L
        SQ1NPKNOT_L_G    = SQ1NPKNOT_L
      END SUBROUTINE SET_Q1NP_COUNTS

      SUBROUTINE SET_Q1NP_TABVINT_LEN(LEN)
        INTEGER, INTENT(IN) :: LEN
        TABVINT_LEN_G = LEN
      END SUBROUTINE SET_Q1NP_TABVINT_LEN

      SUBROUTINE Q1NP_INIT_GAUSS_SCHEME_STARTER(NP_U, NP_V, NP_T)
        INTEGER, INTENT(IN) :: NP_U, NP_V, NP_T
        INTEGER :: I
        
C       Set number of Gauss points
        Q1NP_NP_U_G = NP_U
        Q1NP_NP_V_G = NP_V
        Q1NP_NP_T_G = NP_T
        
C       Deallocate if already allocated
        IF (ALLOCATED(Q1NP_GP_U_G)) DEALLOCATE(Q1NP_GP_U_G)
        IF (ALLOCATED(Q1NP_GP_V_G)) DEALLOCATE(Q1NP_GP_V_G)
        IF (ALLOCATED(Q1NP_GP_T_G)) DEALLOCATE(Q1NP_GP_T_G)
        IF (ALLOCATED(Q1NP_GW_U_G)) DEALLOCATE(Q1NP_GW_U_G)
        IF (ALLOCATED(Q1NP_GW_V_G)) DEALLOCATE(Q1NP_GW_V_G)
        IF (ALLOCATED(Q1NP_GW_T_G)) DEALLOCATE(Q1NP_GW_T_G)
        
C       Allocate arrays
        IF (NP_U > 0) THEN
          ALLOCATE(Q1NP_GP_U_G(NP_U))
          ALLOCATE(Q1NP_GW_U_G(NP_U))
          CALL Q1NP_GAUSS_1D(NP_U, Q1NP_GP_U_G, Q1NP_GW_U_G)
        ENDIF
        
        IF (NP_V > 0) THEN
          ALLOCATE(Q1NP_GP_V_G(NP_V))
          ALLOCATE(Q1NP_GW_V_G(NP_V))
          CALL Q1NP_GAUSS_1D(NP_V, Q1NP_GP_V_G, Q1NP_GW_V_G)
        ENDIF
        
        IF (NP_T > 0) THEN
          ALLOCATE(Q1NP_GP_T_G(NP_T))
          ALLOCATE(Q1NP_GW_T_G(NP_T))
          CALL Q1NP_GAUSS_1D(NP_T, Q1NP_GP_T_G, Q1NP_GW_T_G)
        ENDIF
        
      END SUBROUTINE Q1NP_INIT_GAUSS_SCHEME_STARTER
      
      SUBROUTINE Q1NP_GAUSS_1D(N, GP, GW)
        INTEGER, INTENT(IN) :: N
        my_real, INTENT(OUT) :: GP(N), GW(N)
        INTEGER :: I
        
C       Simple Gauss-Legendre quadrature for [-1,1]
C       For small orders, use hardcoded values
        IF (N == 1) THEN
          GP(1) = 0.0
          GW(1) = 2.0
        ELSEIF (N == 2) THEN
          GP(1) = -0.577350269189626
          GP(2) =  0.577350269189626
          GW(1) = 1.0
          GW(2) = 1.0
        ELSEIF (N == 3) THEN
          GP(1) = -0.774596669241483
          GP(2) =  0.0
          GP(3) =  0.774596669241483
          GW(1) = 0.555555555555556
          GW(2) = 0.888888888888889
          GW(3) = 0.555555555555556
        ELSEIF (N == 4) THEN
          GP(1) = -0.861136311594053
          GP(2) = -0.339981043584856
          GP(3) =  0.339981043584856
          GP(4) =  0.861136311594053
          GW(1) = 0.347854845137454
          GW(2) = 0.652145154862546
          GW(3) = 0.652145154862546
          GW(4) = 0.347854845137454
        ELSEIF (N == 5) THEN
          GP(1) = -0.906179845938664
          GP(2) = -0.538469310105683
          GP(3) =  0.0
          GP(4) =  0.538469310105683
          GP(5) =  0.906179845938664
          GW(1) = 0.236926885056189
          GW(2) = 0.478628670499366
          GW(3) = 0.568888888888889
          GW(4) = 0.478628670499366
          GW(5) = 0.236926885056189
        ELSE
C         For higher orders, use uniform distribution as fallback
          DO I = 1, N
            GP(I) = -1.0 + 2.0 * REAL(I-1) / REAL(N-1)
            GW(I) = 2.0 / REAL(N)
          ENDDO
        ENDIF
        
      END SUBROUTINE Q1NP_GAUSS_1D

      END MODULE Q1NP_RESTART_MOD

