!====================================================================
!  W_Q1NP_INT / W_Q1NP_REAL   starter/source/restart/ddsplit/w_q1np.F
!====================================================================
      MODULE W_Q1NP_MOD
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
        USE PRECISION_MOD, ONLY : WP
        USE Q1NP_RESTART_MOD
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
        IMPLICIT NONE
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Public interface
! ----------------------------------------------------------------------------------------------------------------------
        PUBLIC :: W_Q1NP_INT
        PUBLIC :: W_Q1NP_REAL
        PUBLIC :: Q1NP_REAL_DATA_CHECK

      CONTAINS

        SUBROUTINE W_Q1NP_INT(NUMELQ1NP_IN, NKQ1NP, KQ1NP_TAB, &
          IQ1NP_TAB, IQ1NP_BULK_TAB,                            &
          LEN_IA)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Dummy arguments
! ----------------------------------------------------------------------------------------------------------------------
          INTEGER, INTENT(IN)    :: NUMELQ1NP_IN, NKQ1NP
          INTEGER, INTENT(IN)    :: KQ1NP_TAB(NKQ1NP, *)
          INTEGER, INTENT(IN)    :: IQ1NP_TAB(*)
          INTEGER, INTENT(IN)    :: IQ1NP_BULK_TAB(*)
          INTEGER, INTENT(INOUT) :: LEN_IA
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          INTEGER :: BASE_INT_COUNT, CTRL_INT_COUNT, BULK_INT_COUNT
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          IF (NUMELQ1NP_IN <= 0) RETURN

          BASE_INT_COUNT = NKQ1NP * NUMELQ1NP_IN
          CTRL_INT_COUNT = MAX(0, SIQ1NP_G)
          BULK_INT_COUNT = MAX(0, SQ1NPBULK_G)

          WRITE(*, '(A)') ' '
          WRITE(*, '(A)') ' ** W_Q1NP_INT: writing restart blocks'
          WRITE(*, '(A,I10,A)') '    NUMELQ1NP = ', NUMELQ1NP_IN, &
            '  ! Number of Q1Np elements'
          WRITE(*, '(A,I10,A)') '    NKQ1NP    = ', NKQ1NP,       &
            '  ! Q1Np metadata per element'
          WRITE(*, '(A,I10,A)') '    SKQ1NP_G  = ', SKQ1NP_G,     &
            '  ! Size of KQ1NP_TAB (expected NKQ1NP*NUMELQ1NP)'
          WRITE(*, '(A,I10,A)') '    SIQ1NP_G  = ', SIQ1NP_G,     &
            '  ! Number of control point connectivity ints to write'
          WRITE(*, '(A,I10,A)') '    SQ1NPBULK_G = ', SQ1NPBULK_G,&
            '  ! Number of legacy HEX8 bottom nodes to write'

          IF (SKQ1NP_G > 0 .AND. SKQ1NP_G /= BASE_INT_COUNT) THEN
            WRITE(*, '(A,I10,A,I10)')                           &
              '    WARNING: Expected base ints=', BASE_INT_COUNT,  &
              ' but SKQ1NP_G=', SKQ1NP_G
          END IF

          CALL WRITE_I_C(KQ1NP_TAB, BASE_INT_COUNT)
          LEN_IA = LEN_IA + BASE_INT_COUNT

          IF (CTRL_INT_COUNT > 0) THEN
            CALL WRITE_I_C(IQ1NP_TAB, CTRL_INT_COUNT)
            LEN_IA = LEN_IA + CTRL_INT_COUNT
          END IF

          IF (BULK_INT_COUNT > 0) THEN
            CALL WRITE_I_C(IQ1NP_BULK_TAB, BULK_INT_COUNT)
            LEN_IA = LEN_IA + BULK_INT_COUNT
          END IF

        END SUBROUTINE W_Q1NP_INT

        SUBROUTINE W_Q1NP_REAL(Q1NP_WTAB, Q1NP_KTAB, Q1NP_CPTAB, LEN_AM)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Dummy arguments
! ----------------------------------------------------------------------------------------------------------------------
          REAL(KIND=WP), INTENT(IN)    :: Q1NP_WTAB(*)
          REAL(KIND=WP), INTENT(IN)    :: Q1NP_KTAB(*)
          REAL(KIND=WP), INTENT(IN)    :: Q1NP_CPTAB(3, *)
          INTEGER,       INTENT(INOUT) :: LEN_AM
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          INTEGER        :: WEIGHT_COUNT, KNOT_COUNT, CP_COUNT
          INTEGER        :: NGAUSS_META
          INTEGER        :: MAX_GAUSS_META
          REAL(KIND=WP), ALLOCATABLE :: GAUSS_META(:)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          WEIGHT_COUNT = MAX(0, SQ1NPWEIGHT_L_G)
          KNOT_COUNT   = MAX(0, SQ1NPKNOT_L_G)
          CP_COUNT     = MAX(0, SQ1NPCTRL_SHARED_G)

          CALL Q1NP_REAL_DATA_CHECK(WEIGHT_COUNT, KNOT_COUNT, CP_COUNT, &
            Q1NP_WTAB, Q1NP_KTAB, Q1NP_CPTAB)

          IF (WEIGHT_COUNT + KNOT_COUNT + CP_COUNT <= 0) RETURN

          WRITE(*, '(A)') ' '
          WRITE(*, '(A)') ' ** W_Q1NP_REAL: writing restart blocks'
          WRITE(*, '(A,I10,A)') '    SQ1NPWEIGHT_L_G      = ',     &
            SQ1NPWEIGHT_L_G, '  ! Weights count'
          WRITE(*, '(A,I10,A)') '    SQ1NPKNOT_L_G        = ',     &
            SQ1NPKNOT_L_G, '  ! Knot count'
          WRITE(*, '(A,I10,A)') '    SQ1NPCTRL_SHARED_G   = ',     &
            SQ1NPCTRL_SHARED_G, '  ! Shared CPs'
          WRITE(*, '(A,I10,A)') '    3*SQ1NPCTRL_SHARED_G = ',     &
            3 * SQ1NPCTRL_SHARED_G, '  ! xyz reals'

          IF (WEIGHT_COUNT > 0) THEN
            CALL WRITE_DB(Q1NP_WTAB, WEIGHT_COUNT)
            LEN_AM = LEN_AM + WEIGHT_COUNT
          END IF

          IF (KNOT_COUNT > 0) THEN
            CALL WRITE_DB(Q1NP_KTAB, KNOT_COUNT)
            LEN_AM = LEN_AM + KNOT_COUNT
          END IF

          IF (CP_COUNT > 0) THEN
            CALL WRITE_DB(Q1NP_CPTAB, 3 * CP_COUNT)
            LEN_AM = LEN_AM + 3 * CP_COUNT
          END IF

        END SUBROUTINE W_Q1NP_REAL

        SUBROUTINE Q1NP_REAL_DATA_CHECK(WEIGHT_COUNT, KNOT_COUNT, CP_COUNT, &
          Q1NP_WTAB, Q1NP_KTAB, Q1NP_CPTAB)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Dummy arguments
! ----------------------------------------------------------------------------------------------------------------------
          INTEGER,       INTENT(IN) :: WEIGHT_COUNT, KNOT_COUNT, CP_COUNT
          REAL(KIND=WP), INTENT(IN) :: Q1NP_WTAB(*)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_KTAB(*)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_CPTAB(3, *)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          LOGICAL        :: HAS_WEIGHT_DATA, HAS_KNOT_DATA, HAS_CP_DATA
          INTEGER        :: I
          REAL(KIND=WP)  :: DATA_TOL
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          DATA_TOL        = 1.0_wp * 1.0e-12_wp
          HAS_WEIGHT_DATA = .FALSE.
          HAS_KNOT_DATA   = .FALSE.
          HAS_CP_DATA     = .FALSE.

          IF (WEIGHT_COUNT > 0) THEN
            DO I = 1, WEIGHT_COUNT
              IF (ABS(Q1NP_WTAB(I)) > DATA_TOL) THEN
                HAS_WEIGHT_DATA = .TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT. HAS_WEIGHT_DATA) THEN
              WRITE(*, '(A)')                                           &
                ' ** WARNING: Q1Np weight buffer allocated but empty.'
            END IF
          ELSE
            WRITE(*, '(A)')                                            &
              ' ** WARNING: Q1Np weights count is zero.'
          END IF

          IF (KNOT_COUNT > 0) THEN
            DO I = 1, KNOT_COUNT
              IF (ABS(Q1NP_KTAB(I)) > DATA_TOL) THEN
                HAS_KNOT_DATA = .TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT. HAS_KNOT_DATA) THEN
              WRITE(*, '(A)')                                          &
                ' ** WARNING: Q1Np knot buffer allocated but empty.'
            END IF
          ELSE
            WRITE(*, '(A)')                                            &
              ' ** WARNING: Q1Np knot count is zero.'
          END IF

          IF (CP_COUNT > 0) THEN
            DO I = 1, CP_COUNT
              IF (ABS(Q1NP_CPTAB(1, I)) > DATA_TOL .OR.                   &
                  ABS(Q1NP_CPTAB(2, I)) > DATA_TOL .OR.                   &
                  ABS(Q1NP_CPTAB(3, I)) > DATA_TOL) THEN
                HAS_CP_DATA = .TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT. HAS_CP_DATA) THEN
              WRITE(*, '(A)')                                          &
                ' ** WARNING: Q1Np control-point buffer allocated but empty.'
            END IF
          ELSE
            WRITE(*, '(A)')                                            &
              ' ** WARNING: Q1Np control-point count is zero.'
          END IF

          IF (HAS_WEIGHT_DATA .AND. HAS_KNOT_DATA .AND. HAS_CP_DATA) THEN
            WRITE(*, '(A)') ' ** Q1Np real data buffers contain values.'
          END IF

        END SUBROUTINE Q1NP_REAL_DATA_CHECK

      END MODULE W_Q1NP_MOD

