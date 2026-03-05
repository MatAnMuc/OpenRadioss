!||====================================================================
!||    sts_gp_state_mod  ../engine/source/interfaces/ists/ists_gp_state_mod.F90
!||--------------------------------------------------------------------
!||  Global Gauss point state for new STS contact
!||  - 1D arrays indexed by a global Gauss point index:
!||      gp_index = (pair_index-1) * (ip*ip) + local_gp_index + 1
!||    where local_gp_index = (z-1)*ip + (q-1)
!||====================================================================
      MODULE sts_gp_state_mod

      IMPLICIT NONE

!----------------------------------------------------------------------
! Public state
!----------------------------------------------------------------------
      INTEGER :: MAX_GLOBAL_GP = 0

! 1D global Gauss point coordinates (convective space)
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI1_GLOBAL
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI2_GLOBAL
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI1_GLOBAL_PREV
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI2_GLOBAL_PREV

! Period counters for convective coordinates
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI1_PERIOD
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI2_PERIOD

! Tangential traction history (1D, per global Gauss point)
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_TTRIAL1_HIST
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_TTRIAL2_HIST

! Sticking/sliding and initialization flags
      LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: GP_IS_STICKING
      LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: GP_INITIALIZED

!----------------------------------------------------------------------
! Interface
!----------------------------------------------------------------------
      CONTAINS

!----------------------------------------------------------------------
! Initialize or resize global Gauss point state arrays
!
! Arguments:
!   MAX_STS_SIZE_ACTUAL : maximum number of STS contact pairs
!   IP_MAX              : quadrature order per direction
!----------------------------------------------------------------------
      SUBROUTINE sts_gp_state_init(MAX_STS_SIZE_ACTUAL, IP_MAX)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAX_STS_SIZE_ACTUAL
      INTEGER, INTENT(IN) :: IP_MAX

      INTEGER :: NEW_MAX_GLOBAL_GP

!     Compute required global GP count
      NEW_MAX_GLOBAL_GP = MAX_STS_SIZE_ACTUAL * IP_MAX * IP_MAX

!     If no pairs, release any allocated storage
      IF (NEW_MAX_GLOBAL_GP <= 0) THEN
        IF (ALLOCATED(GP_XI1_GLOBAL)) THEN
          DEALLOCATE(GP_XI1_GLOBAL)
          DEALLOCATE(GP_XI2_GLOBAL)
          DEALLOCATE(GP_XI1_GLOBAL_PREV)
          DEALLOCATE(GP_XI2_GLOBAL_PREV)

          DEALLOCATE(GP_XI1_PERIOD)
          DEALLOCATE(GP_XI2_PERIOD)

          DEALLOCATE(GP_TTRIAL1_HIST)
          DEALLOCATE(GP_TTRIAL2_HIST)

          DEALLOCATE(GP_IS_STICKING)
          DEALLOCATE(GP_INITIALIZED)
        END IF
        MAX_GLOBAL_GP = 0
        RETURN
      END IF

!     Reallocate if size changed or not yet allocated
      IF (.NOT. ALLOCATED(GP_XI1_GLOBAL) .OR. NEW_MAX_GLOBAL_GP /= MAX_GLOBAL_GP) THEN

!       Deallocate previous storage if present
        IF (ALLOCATED(GP_XI1_GLOBAL)) THEN
          DEALLOCATE(GP_XI1_GLOBAL)
          DEALLOCATE(GP_XI2_GLOBAL)
          DEALLOCATE(GP_XI1_GLOBAL_PREV)
          DEALLOCATE(GP_XI2_GLOBAL_PREV)

          DEALLOCATE(GP_XI1_PERIOD)
          DEALLOCATE(GP_XI2_PERIOD)

          DEALLOCATE(GP_TTRIAL1_HIST)
          DEALLOCATE(GP_TTRIAL2_HIST)

          DEALLOCATE(GP_IS_STICKING)
          DEALLOCATE(GP_INITIALIZED)
        END IF

!       Allocate with new size
        MAX_GLOBAL_GP = NEW_MAX_GLOBAL_GP

        ALLOCATE(GP_XI1_GLOBAL      (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI2_GLOBAL      (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI1_GLOBAL_PREV (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI2_GLOBAL_PREV (MAX_GLOBAL_GP))

        ALLOCATE(GP_XI1_PERIOD      (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI2_PERIOD      (MAX_GLOBAL_GP))

        ALLOCATE(GP_TTRIAL1_HIST    (MAX_GLOBAL_GP))
        ALLOCATE(GP_TTRIAL2_HIST    (MAX_GLOBAL_GP))

        ALLOCATE(GP_IS_STICKING     (MAX_GLOBAL_GP))
        ALLOCATE(GP_INITIALIZED     (MAX_GLOBAL_GP))

!       Initialize state
        GP_XI1_GLOBAL       = 0.0D0
        GP_XI2_GLOBAL       = 0.0D0
        GP_XI1_GLOBAL_PREV  = 0.0D0
        GP_XI2_GLOBAL_PREV  = 0.0D0

        GP_XI1_PERIOD       = 0
        GP_XI2_PERIOD       = 0

        GP_TTRIAL1_HIST     = 0.0D0
        GP_TTRIAL2_HIST     = 0.0D0

        GP_IS_STICKING      = .FALSE.
        GP_INITIALIZED      = .FALSE.

      END IF

      END SUBROUTINE sts_gp_state_init

      END MODULE sts_gp_state_mod

