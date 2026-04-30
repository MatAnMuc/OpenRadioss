!||====================================================================
!||    q1np_contact_driver               ../engine/source/interfaces/int26/q1np_contact_driver.F90
!||--- called by ------------------------------------------------------
!||    i7mainf                           ../engine/source/interfaces/int07/i7mainf.F
!||--- calls      -----------------------------------------------------
!||    q1np_contact_broad_phase_check_proximity
!||         ../engine/source/interfaces/int26/q1np_contact_algorithms.F90
!||--- uses       -----------------------------------------------------
!||    q1np_restart_mod                  ../common_source/modules/q1np_restart_mod.F90
!||    restmod                           ../engine/share/modules/restart_mod.F
!||====================================================================
      MODULE Q1NP_CONTACT_DRIVER_MOD
        USE PRECISION_MOD, ONLY : WP
        USE Q1NP_RESTART_MOD
        USE RESTMOD, ONLY : KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB
        USE Q1NP_CONTACT_ALGORITHMS_MOD, ONLY : Q1NP_CONTACT_BROAD_PHASE_CHECK_PROXIMITY
        IMPLICIT NONE
        PRIVATE

        INTEGER, SAVE :: Q1NP_CONTACT_INT7_LAST_NCYCLE = -1
        LOGICAL, SAVE :: Q1NP_CONTACT_INT7_ALREADY = .FALSE.

        PUBLIC :: Q1NP_CONTACT_DATA_READY
        PUBLIC :: Q1NP_CONTACT_DRIVER_INT7

      CONTAINS

!=======================================================================
!   Q1NP_CONTACT_DATA_READY
!   True when Q1NP elements and restart tables exist.
!=======================================================================
        LOGICAL FUNCTION Q1NP_CONTACT_DATA_READY()
          Q1NP_CONTACT_DATA_READY = (NUMELQ1NP_G > 0 .AND. &
     &        Q1NP_NKNOT_SETS_G >= 2 .AND. &
     &        ALLOCATED(KQ1NP_TAB) .AND. &
     &        ALLOCATED(IQ1NP_TAB) .AND. &
     &        ALLOCATED(Q1NP_KTAB))
        END FUNCTION Q1NP_CONTACT_DATA_READY

!=======================================================================
!   Q1NP_CONTACT_DRIVER_INT7
!   INT7 entry: runs broad+narrow+penalty at most once per NCYCLE.
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_DRIVER_INT7(NCYCLE, NUMNOD, X, A, STIFN)
          INTEGER, INTENT(IN) :: NCYCLE, NUMNOD
          REAL(KIND=WP), INTENT(IN) :: X(3,NUMNOD)
          REAL(KIND=WP), INTENT(INOUT) :: A(3,NUMNOD), STIFN(NUMNOD)
          LOGICAL :: HIT

          IF (.NOT. Q1NP_CONTACT_DATA_READY()) RETURN

          IF (NCYCLE /= Q1NP_CONTACT_INT7_LAST_NCYCLE) THEN
            Q1NP_CONTACT_INT7_LAST_NCYCLE = NCYCLE
            Q1NP_CONTACT_INT7_ALREADY = .FALSE.
          END IF
          IF (Q1NP_CONTACT_INT7_ALREADY) RETURN

          CALL Q1NP_CONTACT_BROAD_PHASE_CHECK_PROXIMITY( &
     &      KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &      X, NUMNOD, NUMELQ1NP_G, A, STIFN, HIT)
          Q1NP_CONTACT_INT7_ALREADY = .TRUE.
        END SUBROUTINE Q1NP_CONTACT_DRIVER_INT7

      END MODULE Q1NP_CONTACT_DRIVER_MOD
