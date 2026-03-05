!||====================================================================
!||    REMAP_SEGMENTS_STS  ../engine/source/interfaces/ists/ists_remap_segments.F90
!||--- called by ------------------------------------------------------
!||    i7mainf              ../engine/source/interfaces/int07/i7mainf.F
!||--- calls ---------------------------------------------------------
!||    (none – local mapping only)
!||====================================================================
      SUBROUTINE STS_REMAP_SEGMENTS(INTBUF_TAB, ITAB, X, CAND_SEC_SEG, &
     &  IRECT, CONT_ELEMENT, COUNT, IGRSURF, CAND_SEC_SEG_ID, &
     &  CAND_MST_SEG_ID, MAX_STS_SIZE_ACTUAL)
!-----------------------------------------------
!   M o d u l e s
!----------------------------------------------- 
      USE INTBUFDEF_MOD
      USE GROUPDEF_MOD
!-----------------------------------------------
!   M o d u l e s   /   I m p l i c i t   T y p e s
!-----------------------------------------------
      use constant_mod
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      TYPE (SURF_)   , DIMENSION(*)   :: IGRSURF
      INTEGER ITAB(*)
      INTEGER IRECT(4,*)
      my_real X(3,*)
      INTEGER CAND_SEC_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_MST_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      my_real CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      INTEGER COUNT, MAX_STS_SIZE_ACTUAL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, N, NI
      INTEGER candidate, candidateM
      LOGICAL :: duplicate
      INTEGER, ALLOCATABLE :: IGRSURF_S_TEMP(:,:), IGRSURF_M_TEMP(:,:)
      INTEGER, ALLOCATABLE :: found_segments(:,:)
      INTEGER, ALLOCATABLE :: ADJA(:,:)
!-----------------------------------------------
!   S o u r c e   L i n e s
!-----------------------------------------------
      ! Safety checks
      IF (IGRSURF(1)%NSEG <= 0 .OR. IGRSURF(2)%NSEG <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      IF (.NOT. ALLOCATED(IGRSURF(1)%NODES) .OR. &
     &    .NOT. ALLOCATED(IGRSURF(2)%NODES)) THEN
        COUNT = 0
        RETURN
      END IF
      IF (INTBUF_TAB(1)%S_CAND_N <= 0) THEN
        COUNT = 0
        RETURN
      END IF

      ! Copy surface nodes to temporary arrays
      ALLOCATE(IGRSURF_S_TEMP(IGRSURF(1)%NSEG, 4))
      ALLOCATE(IGRSURF_M_TEMP(IGRSURF(2)%NSEG, 4))
      IGRSURF_S_TEMP = IGRSURF(1)%NODES
      IGRSURF_M_TEMP = IGRSURF(2)%NODES

      ! Initialize
      ALLOCATE(found_segments(2, MAX_STS_SIZE_ACTUAL))
      found_segments = -1  ! Initialize with -1 to indicate not found   
      N = 1
      CAND_SEC_SEG = -1
      CAND_MST_SEG = -1
      ! Compare candidate array with IGRSURF_TEMP table to derive the respective index
      COUNT = 0

      ! Map candidate nodes to segment pairs
      DO I = 1, INTBUF_TAB(1)%S_CAND_N
        candidate = INTBUF_TAB(1)%NSV(INTBUF_TAB(1)%CAND_N(I))
        candidateM = INTBUF_TAB(1)%CAND_E(I)
        
        ! Find which secondary segment contains this candidate node
        DO J = 1, IGRSURF(1)%NSEG
          ! Check whether the current candidate node (candidate) is among the four corner nodes
          ! of the current secondary segment (J). If yes, this segment is considered relevant.
          IF (ANY(candidate == IGRSURF_S_TEMP(J, 1:4))) THEN
            
            ! Check if this segment pair already exists
            duplicate = .FALSE.
            DO K = 1, COUNT
              IF (CAND_SEC_SEG(K) == J .AND. CAND_MST_SEG(K) == candidateM) THEN
                duplicate = .TRUE.
                EXIT
              END IF
            END DO
            
            ! Add new unique segment pair
            IF (.NOT. duplicate) THEN
              COUNT = COUNT + 1
              IF (COUNT > MAX_STS_SIZE_ACTUAL) THEN
                COUNT = COUNT - 1
                EXIT
              END IF
              
              found_segments(1, N) = J
              found_segments(2, N) = candidateM

              CAND_SEC_SEG(COUNT) = found_segments(1, N)
              CAND_MST_SEG(COUNT) = found_segments(2, N)

              N = N + 1
            END IF
          END IF
        END DO
      END DO
      
      ! Adjacency matrix - filter out invalid entries
      ALLOCATE(ADJA(COUNT, 2))
      J = 0
      DO I = 1, COUNT
        IF (CAND_SEC_SEG(I) /= 0 .AND. CAND_MST_SEG(I) /= 0) THEN
          J = J + 1
          ADJA(J, 1) = CAND_SEC_SEG(I)
          ADJA(J, 2) = CAND_MST_SEG(I)
        END IF
      END DO
      COUNT = K

      COUNT = J
      J = 1
      
      ! Build segment ID arrays
      DO I = 1, COUNT
        ! Secondary segment nodes
        CAND_SEC_SEG_ID(I, 1) = CAND_SEC_SEG(I)
        CAND_SEC_SEG_ID(I, 2:5) = IGRSURF(1)%NODES(CAND_SEC_SEG(I), 1:4)

        ! Primary segment nodes
        CAND_MST_SEG_ID(I, 1) = CAND_MST_SEG(I)
        CAND_MST_SEG_ID(I, 2:5) = IRECT(1:4, CAND_MST_SEG(I))
      END DO

      ! Store coordinates: Primary (1-4), Secondary (5-8)
      ! PRIMARY -> FIRST (1-4)
      DO I = 1, COUNT
        J = 1
        DO K = 2, 5
          NI = CAND_MST_SEG_ID(I, K)
          CONT_ELEMENT(I, 1, J) = X(1, NI)  ! X
          CONT_ELEMENT(I, 2, J) = X(2, NI)  ! Y
          CONT_ELEMENT(I, 3, J) = X(3, NI)  ! Z
          J = J + 1
        END DO
      END DO
          
      ! SECONDARY -> SECOND (5-8)
      DO I = 1, COUNT
        J = 5
        DO K = 2, 5
          NI = CAND_SEC_SEG_ID(I, K)
          CONT_ELEMENT(I, 1, J) = X(1, NI)  ! X
          CONT_ELEMENT(I, 2, J) = X(2, NI)  ! Y
          CONT_ELEMENT(I, 3, J) = X(3, NI)  ! Z
          J = J + 1
        END DO
      END DO
      
      ! Cleanup
      DEALLOCATE(IGRSURF_S_TEMP)
      DEALLOCATE(IGRSURF_M_TEMP)
      DEALLOCATE(found_segments)
      DEALLOCATE(ADJA)
      
      RETURN
      END