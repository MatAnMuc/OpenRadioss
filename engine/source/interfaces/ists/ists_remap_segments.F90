!||====================================================================
!||    REMAP_SEGMENTS_STS  ../engine/source/interfaces/ists/ists_remap_segments.F90
!||--- called by ------------------------------------------------------
!||    i7mainf              ../engine/source/interfaces/int07/i7mainf.F
!||--- calls ---------------------------------------------------------
!||    (none – local mapping only)
!||====================================================================
      SUBROUTINE STS_REMAP_SEGMENTS(INTBUF_TAB, X, NUMNOD, NRTM, CAND_SEC_SEG, &
     &  JLT, CAND_N_CUR, CAND_E_CUR, IRECT, CONT_ELEMENT, COUNT, &
     &  IGRSURF, CAND_SEC_SEG_ID, CAND_MST_SEG_ID, &
     &  MAX_STS_SIZE_ACTUAL, NSURF_LOCAL, SEC_SURF_ID, MST_SURF_ID)
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
      TYPE(INTBUF_STRUCT_) INTBUF_TAB
      TYPE (SURF_)   , DIMENSION(NSURF_LOCAL)   :: IGRSURF
      INTEGER JLT, NUMNOD, NRTM, CAND_N_CUR(JLT), CAND_E_CUR(JLT)
      INTEGER IRECT(4,NRTM)
      my_real X(3,NUMNOD)
      INTEGER CAND_SEC_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_MST_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      my_real CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      INTEGER COUNT, MAX_STS_SIZE_ACTUAL
      INTEGER NSURF_LOCAL, SEC_SURF_ID, MST_SURF_ID
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, N, NI
      INTEGER candidate, candidateM
      INTEGER BEST_OVERLAP, OVERLAP, SEC_SURF_IDX
      LOGICAL :: duplicate
      INTEGER, ALLOCATABLE :: IGRSURF_S_TEMP(:,:)
      INTEGER, ALLOCATABLE :: found_segments(:,:)
      INTEGER, ALLOCATABLE :: ADJA(:,:)
!-----------------------------------------------
!   S o u r c e   L i n e s
!-----------------------------------------------
      SEC_SURF_IDX = SEC_SURF_ID
      IF (INTBUF_TAB%S_NSV <= 0 .OR. JLT <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      IF (SEC_SURF_IDX <= 0 .OR. SEC_SURF_IDX > NSURF_LOCAL) THEN
        BEST_OVERLAP = 0
        SEC_SURF_IDX = 0
        DO I = 1, NSURF_LOCAL
          IF (I == MST_SURF_ID) CYCLE
          IF (IGRSURF(I)%NSEG <= 0) CYCLE
          IF (.NOT. ALLOCATED(IGRSURF(I)%NODES)) CYCLE
          OVERLAP = 0
          DO J = 1, IGRSURF(I)%NSEG
            DO K = 1, 4
              candidate = IGRSURF(I)%NODES(J, K)
              IF (candidate <= 0) CYCLE
              IF (ANY(INTBUF_TAB%NSV(1:INTBUF_TAB%S_NSV) == candidate)) THEN
                OVERLAP = OVERLAP + 1
              END IF
            END DO
          END DO
          IF (OVERLAP > BEST_OVERLAP) THEN
            BEST_OVERLAP = OVERLAP
            SEC_SURF_IDX = I
          END IF
        END DO
      END IF
      IF (SEC_SURF_IDX <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      ! Safety checks
      IF (IGRSURF(SEC_SURF_IDX)%NSEG <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      IF (.NOT. ALLOCATED(IGRSURF(SEC_SURF_IDX)%NODES)) THEN
        COUNT = 0
        RETURN
      END IF
      IF (JLT <= 0) THEN
        COUNT = 0
        RETURN
      END IF

      ! Copy surface nodes to temporary arrays
      ALLOCATE(IGRSURF_S_TEMP(IGRSURF(SEC_SURF_IDX)%NSEG, 4))
      IGRSURF_S_TEMP = IGRSURF(SEC_SURF_IDX)%NODES

      ! Initialize
      ALLOCATE(found_segments(2, MAX_STS_SIZE_ACTUAL))
      found_segments = -1  ! Initialize with -1 to indicate not found   
      N = 1
      CAND_SEC_SEG = -1
      CAND_MST_SEG = -1
      ! Compare candidate array with IGRSURF_TEMP table to derive the respective index
      COUNT = 0

      ! Map candidate nodes to segment pairs
      DO I = 1, JLT
        IF (CAND_N_CUR(I) <= 0 .OR. CAND_N_CUR(I) > INTBUF_TAB%S_NSV) CYCLE
        candidate = INTBUF_TAB%NSV(CAND_N_CUR(I))
        candidateM = CAND_E_CUR(I)
        IF (candidateM <= 0) CYCLE
        
        ! Find which secondary segment contains this candidate node
        DO J = 1, IGRSURF(SEC_SURF_IDX)%NSEG
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

      IF (COUNT <= 0) THEN
        DEALLOCATE(IGRSURF_S_TEMP)
        DEALLOCATE(found_segments)
        RETURN
      END IF

      ! Adjacency matrix - filter out invalid entries
      ALLOCATE(ADJA(COUNT, 2))
      J = 0
      DO I = 1, COUNT
        IF (CAND_SEC_SEG(I) > 0 .AND. CAND_MST_SEG(I) > 0) THEN
          J = J + 1
          ADJA(J, 1) = CAND_SEC_SEG(I)
          ADJA(J, 2) = CAND_MST_SEG(I)
        END IF
      END DO
      COUNT = J
      J = 1
      
      ! Build segment ID arrays
      DO I = 1, COUNT
        ! Secondary segment nodes
        CAND_SEC_SEG_ID(I, 1) = CAND_SEC_SEG(I)
        CAND_SEC_SEG_ID(I, 2:5) = IGRSURF(SEC_SURF_IDX)%NODES(CAND_SEC_SEG(I), 1:4)

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
      DEALLOCATE(found_segments)
      DEALLOCATE(ADJA)
      
      RETURN
      END