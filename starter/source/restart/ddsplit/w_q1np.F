Chd|====================================================================
Chd|  W_Q1NP_INT / W_Q1NP_REAL   source/restart/ddsplit/w_q1np.F
Chd|-- called by -----------
Chd|        DDSPLIT                       source/restart/ddsplit/ddsplit.F
Chd|-- calls ---------------
Chd|        WRITE_I_C                     ../common_source/tools/input_output/write_routtines.c
Chd|        WRITE_DB                      source/restart/ddsplit/wrrest.F
Chd|====================================================================
      SUBROUTINE W_Q1NP_INT(NUMELQ1NP_IN,NKQ1NP,KQ1NP_TAB,
     .                      IQ1NP_TAB,IQ1NP_BULK_TAB,
     .                      LEN_IA)
      USE Q1NP_RESTART_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "units_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NUMELQ1NP_IN,NKQ1NP
      INTEGER, INTENT(IN) :: KQ1NP_TAB(NKQ1NP,*)
      INTEGER, INTENT(IN) :: IQ1NP_TAB(*)
      INTEGER, INTENT(IN) :: IQ1NP_BULK_TAB(*)
      INTEGER, INTENT(INOUT) :: LEN_IA
C-----------------------------------------------
C
      INTEGER :: BASE_INT_COUNT, CTRL_INT_COUNT, BULK_INT_COUNT

      IF (NUMELQ1NP_IN <= 0) RETURN

      BASE_INT_COUNT = NKQ1NP*NUMELQ1NP_IN
      CTRL_INT_COUNT = MAX(0,SIQ1NP_G)
      BULK_INT_COUNT = MAX(0,SQ1NPBULK_G)

      WRITE(IOUT,'(A)') ' '
      WRITE(IOUT,'(A)') ' ** W_Q1NP_INT: writing restart blocks'
      WRITE(IOUT,'(A,I10,A)') '    NUMELQ1NP = ', NUMELQ1NP_IN, '  ! Number of Q1Np elements'
      WRITE(IOUT,'(A,I10,A)') '    NKQ1NP    = ', NKQ1NP,    '  ! Q1Np metadata per element'
      WRITE(IOUT,'(A,I10,A)') '    SKQ1NP_G  = ', SKQ1NP_G,  '  ! Size of KQ1NP_TAB (expected NKQ1NP*NUMELQ1NP)'
      WRITE(IOUT,'(A,I10,A)') '    SIQ1NP_G  = ', SIQ1NP_G,  '  ! Number of control point connectivity ints to write'
      WRITE(IOUT,'(A,I10,A)') '    SQ1NPBULK_G = ', SQ1NPBULK_G, '  ! Number of legacy HEX8 bottom nodes to write'

      IF (SKQ1NP_G > 0 .AND. SKQ1NP_G /= BASE_INT_COUNT) THEN
        WRITE(IOUT,'(A,I10,A,I10)') 
     .    '    WARNING: Expected base ints=', BASE_INT_COUNT,
     .    ' but SKQ1NP_G=', SKQ1NP_G
      ENDIF

      CALL WRITE_I_C(KQ1NP_TAB,BASE_INT_COUNT)
      LEN_IA = LEN_IA + BASE_INT_COUNT

      IF (CTRL_INT_COUNT > 0) THEN
        CALL WRITE_I_C(IQ1NP_TAB,CTRL_INT_COUNT)
        LEN_IA = LEN_IA + CTRL_INT_COUNT
      ENDIF

      IF (BULK_INT_COUNT > 0) THEN
        CALL WRITE_I_C(IQ1NP_BULK_TAB,BULK_INT_COUNT)
        LEN_IA = LEN_IA + BULK_INT_COUNT
      ENDIF

      RETURN
      END

      SUBROUTINE W_Q1NP_REAL(Q1NP_WTAB,Q1NP_KTAB,
     .                       Q1NP_CPTAB,LEN_AM)
      USE Q1NP_RESTART_MOD
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "units_c.inc"
C-----------------------------------------------
      my_real, INTENT(IN) :: Q1NP_WTAB(*)
      my_real, INTENT(IN) :: Q1NP_KTAB(*)
      my_real, INTENT(IN) :: Q1NP_CPTAB(3,*)
      INTEGER, INTENT(INOUT) :: LEN_AM
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: WEIGHT_COUNT, KNOT_COUNT, CP_COUNT
      INTEGER :: NGAUSS_META
      INTEGER :: MAX_GAUSS_META
      my_real, ALLOCATABLE :: GAUSS_META(:)

      WEIGHT_COUNT = MAX(0,SQ1NPWEIGHT_L_G)
      KNOT_COUNT   = MAX(0,SQ1NPKNOT_L_G)
      CP_COUNT     = MAX(0,SQ1NPCTRL_SHARED_G)

      CALL Q1NP_REAL_DATA_CHECK(WEIGHT_COUNT,KNOT_COUNT,CP_COUNT,
     .                          Q1NP_WTAB,Q1NP_KTAB,Q1NP_CPTAB)

      IF (WEIGHT_COUNT+KNOT_COUNT+CP_COUNT <= 0) RETURN

      WRITE(IOUT,'(A)') ' '
      WRITE(IOUT,'(A)') ' ** W_Q1NP_REAL: writing restart blocks'
      WRITE(IOUT,'(A,I10,A)') '    SQ1NPWEIGHT_L_G      = ', SQ1NPWEIGHT_L_G, '  ! Weights count'
      WRITE(IOUT,'(A,I10,A)') '    SQ1NPKNOT_L_G        = ', SQ1NPKNOT_L_G, '  ! Knot count'      
      WRITE(IOUT,'(A,I10,A)') '    SQ1NPCTRL_SHARED_G   = ', SQ1NPCTRL_SHARED_G, '  ! Shared CPs'   
      WRITE(IOUT,'(A,I10,A)') '    3*SQ1NPCTRL_SHARED_G = ', 3*SQ1NPCTRL_SHARED_G, '  ! xyz reals'   

C     Pack Gauss scheme (orders + 1D points/weights) into real buffer
      NGAUSS_META = 0
      IF (Q1NP_NP_U_G>0 .OR. Q1NP_NP_V_G>0 .OR. Q1NP_NP_T_G>0) THEN
        NGAUSS_META = 3 + 2*Q1NP_NP_U_G + 2*Q1NP_NP_V_G + 2*Q1NP_NP_T_G
        MAX_GAUSS_META = NGAUSS_META
        ALLOCATE(GAUSS_META(MAX_GAUSS_META))
        GAUSS_META(1) = REAL(Q1NP_NP_U_G)
        GAUSS_META(2) = REAL(Q1NP_NP_V_G)
        GAUSS_META(3) = REAL(Q1NP_NP_T_G)
        IF (Q1NP_NP_U_G > 0) THEN
          GAUSS_META(4:3+Q1NP_NP_U_G) = Q1NP_GP_U_G(1:Q1NP_NP_U_G)
          GAUSS_META(4+Q1NP_NP_U_G:3+2*Q1NP_NP_U_G) =
     .      Q1NP_GW_U_G(1:Q1NP_NP_U_G)
        ENDIF    
        IF (Q1NP_NP_V_G > 0) THEN
          GAUSS_META(4+2*Q1NP_NP_U_G:
     .               3+2*Q1NP_NP_U_G+Q1NP_NP_V_G) =
     .      Q1NP_GP_V_G(1:Q1NP_NP_V_G)
          GAUSS_META(4+2*Q1NP_NP_U_G+Q1NP_NP_V_G:
     .               3+2*Q1NP_NP_U_G+2*Q1NP_NP_V_G) =
     .      Q1NP_GW_V_G(1:Q1NP_NP_V_G)
        ENDIF
        IF (Q1NP_NP_T_G > 0) THEN
          GAUSS_META(4+2*Q1NP_NP_U_G+2*Q1NP_NP_V_G:
     .               3+2*Q1NP_NP_U_G+2*Q1NP_NP_V_G+Q1NP_NP_T_G) =
     .      Q1NP_GP_T_G(1:Q1NP_NP_T_G)
          GAUSS_META(4+2*Q1NP_NP_U_G+2*Q1NP_NP_V_G+Q1NP_NP_T_G:
     .               3+2*Q1NP_NP_U_G+2*Q1NP_NP_V_G+2*Q1NP_NP_T_G) =
     .      Q1NP_GW_T_G(1:Q1NP_NP_T_G)
        ENDIF
      ENDIF

      IF (NGAUSS_META > 0) THEN
        CALL WRITE_DB(GAUSS_META,NGAUSS_META)
        LEN_AM = LEN_AM + NGAUSS_META
        DEALLOCATE(GAUSS_META)
      ENDIF

      IF (WEIGHT_COUNT > 0) THEN
        CALL WRITE_DB(Q1NP_WTAB,WEIGHT_COUNT)
        LEN_AM = LEN_AM + WEIGHT_COUNT
      ENDIF

      IF (KNOT_COUNT > 0) THEN
        CALL WRITE_DB(Q1NP_KTAB,KNOT_COUNT)
        LEN_AM = LEN_AM + KNOT_COUNT
      ENDIF

      IF (CP_COUNT > 0) THEN
        CALL WRITE_DB(Q1NP_CPTAB,3*CP_COUNT)
        LEN_AM = LEN_AM + 3*CP_COUNT
      ENDIF

      RETURN
      END

      SUBROUTINE Q1NP_REAL_DATA_CHECK(WEIGHT_COUNT,KNOT_COUNT,CP_COUNT,
     .                                Q1NP_WTAB,Q1NP_KTAB,Q1NP_CPTAB)
      USE Q1NP_RESTART_MOD
#include      "implicit_f.inc"
#include      "units_c.inc"
      INTEGER, INTENT(IN) :: WEIGHT_COUNT,KNOT_COUNT,CP_COUNT
      my_real, INTENT(IN) :: Q1NP_WTAB(*)
      my_real, INTENT(IN) :: Q1NP_KTAB(*)
      my_real, INTENT(IN) :: Q1NP_CPTAB(3,*)
      LOGICAL :: HAS_WEIGHT_DATA,HAS_KNOT_DATA,HAS_CP_DATA
      INTEGER :: I
      my_real :: DATA_TOL

      DATA_TOL = 1.0D-12
      HAS_WEIGHT_DATA = .FALSE.
      HAS_KNOT_DATA   = .FALSE.
      HAS_CP_DATA     = .FALSE.

      IF (WEIGHT_COUNT > 0) THEN
        DO I=1,WEIGHT_COUNT
          IF (ABS(Q1NP_WTAB(I)) > DATA_TOL) THEN
            HAS_WEIGHT_DATA = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT. HAS_WEIGHT_DATA) THEN
          WRITE(IOUT,'(A)') 
     .      ' ** WARNING: Q1Np weight buffer allocated but empty.'
        ENDIF
      ELSE
        WRITE(IOUT,'(A)') 
     .    ' ** WARNING: Q1Np weights count is zero.'
      ENDIF

      IF (KNOT_COUNT > 0) THEN
        DO I=1,KNOT_COUNT
          IF (ABS(Q1NP_KTAB(I)) > DATA_TOL) THEN
            HAS_KNOT_DATA = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT. HAS_KNOT_DATA) THEN
          WRITE(IOUT,'(A)') 
     .      ' ** WARNING: Q1Np knot buffer allocated but empty.'
        ENDIF
      ELSE
        WRITE(IOUT,'(A)') 
     .    ' ** WARNING: Q1Np knot count is zero.'
      ENDIF

      IF (CP_COUNT > 0) THEN
        DO I=1,CP_COUNT
          IF (ABS(Q1NP_CPTAB(1,I)) > DATA_TOL .OR.
     .        ABS(Q1NP_CPTAB(2,I)) > DATA_TOL .OR.
     .        ABS(Q1NP_CPTAB(3,I)) > DATA_TOL) THEN
            HAS_CP_DATA = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT. HAS_CP_DATA) THEN
          WRITE(IOUT,'(A)') 
     .      ' ** WARNING: Q1Np control-point buffer allocated but empty.'
        ENDIF
      ELSE
        WRITE(IOUT,'(A)') 
     .    ' ** WARNING: Q1Np control-point count is zero.'
      ENDIF

      IF (HAS_WEIGHT_DATA .AND. HAS_KNOT_DATA .AND. HAS_CP_DATA) THEN
        WRITE(IOUT,'(A)') ' ** Q1Np real data buffers contain values.'
      ENDIF

      RETURN
      END

