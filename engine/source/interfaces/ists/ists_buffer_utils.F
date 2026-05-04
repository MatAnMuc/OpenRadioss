C ##############################################################################
C
C WRITE_IGRSURF: Debug output for surface and interface data
C   Prints secondary surface nodes from IGRSURF and INTBUF_TAB%IRECTS
C   OPTIONAL FOR DEBUGGING

      SUBROUTINE WRITE_IGRSURF(INTBUF_TAB,IGRSURF,ITAB) 
C-----------------------------------------------
C   M o d u l e s
C----------------------------------------------- 
      USE INTBUFDEF_MOD
      USE GROUPDEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      TYPE (SURF_)   , DIMENSION(NSURF)   :: IGRSURF
      INTEGER ITAB(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER num_rows(2)
      INTEGER num_irects_size
      INTEGER I, J
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
      write(*,*) "=========================================="
      write(*,*) "WRITE_IGRSURF: Surface and Interface Data"
      write(*,*) "=========================================="
      
C-----------------------------------------------
C Output IGRSURF(1)%NODES: Secondary Surface Segments
C-----------------------------------------------
      IF (NSURF >= 1 .AND. ALLOCATED(IGRSURF(1)%NODES)) THEN
          num_rows = SHAPE(IGRSURF(1)%NODES)
          
          write(*,*)
          write(*,*) "SECONDARY SURFACE (IGRSURF(1)):"
          write(*,*) "  Title: ", IGRSURF(1)%TITLE
          write(*,*) "  Number of Segments: ", num_rows(1)
          write(*,*) "  Nodes per Segment: ", num_rows(2)
          write(*,*)
          write(*,*) "  Segment Node IDs (Local):"
          write(*,*) "  Segment | Node1 | Node2 | Node3 | Node4"
          write(*,*) "  --------|-------|-------|-------|-------"
          
          do I = 1, num_rows(1)
              write(*,'(T3, I7, T13, I6, T20, I6, T27, I6, T34, I6)') 
     &              I, (IGRSURF(1)%NODES(I,J), J=1,4)
          enddo
          
          write(*,*)
          write(*,*) "  Segment Node IDs (Global):"
          write(*,*) "  Segment | Node1 | Node2 | Node3 | Node4"
          write(*,*) "  --------|-------|-------|-------|-------"
          
          do I = 1, num_rows(1)
              write(*,'(T3, I7, T13, I6, T20, I6, T27, I6, T34, I6)') 
     &              I, (ITAB(IGRSURF(1)%NODES(I,J)), J=1,4)
          enddo
      ELSE
          write(*,*) "  WARNING: IGRSURF(1)%NODES not allocated or invalid"
      ENDIF
      
C-----------------------------------------------
C Output INTBUF_TAB(1)%IRECTS: Flattened Node Array
C-----------------------------------------------
      write(*,*)
      write(*,*) "INTBUF_TAB(1)%IRECTS (Flattened Array):"
      
      IF (ALLOCATED(INTBUF_TAB(1)%IRECTS)) THEN
          num_irects_size = SIZE(INTBUF_TAB(1)%IRECTS)
          write(*,*) "  Array Size: ", num_irects_size
          write(*,*) "  S_IRECTS: ", INTBUF_TAB(1)%S_IRECTS
          
          IF (num_irects_size > 0) THEN
              write(*,*)
              write(*,*) "  Node IDs (Local):"
              write(*,'(T3, 10I6)') (INTBUF_TAB(1)%IRECTS(I), I=1, 
     &              MIN(num_irects_size, 50))  ! Limit to first 50 for readability
              IF (num_irects_size > 50) THEN
                  write(*,*) "  ... (showing first 50 of ", 
     &                  num_irects_size, " elements)"
              ENDIF
              
              write(*,*)
              write(*,*) "  Node IDs (Global):"
              write(*,'(T3, 10I6)') (ITAB(INTBUF_TAB(1)%IRECTS(I)), 
     &              I=1, MIN(num_irects_size, 50))
              IF (num_irects_size > 50) THEN
                  write(*,*) "  ... (showing first 50 of ", 
     &                  num_irects_size, " elements)"
              ENDIF
          ELSE
              write(*,*) "  WARNING: IRECTS array is empty"
          ENDIF
      ELSE
          write(*,*) "  WARNING: INTBUF_TAB(1)%IRECTS not allocated"
      ENDIF
      
      write(*,*)
      write(*,*) "=========================================="
      
C-----------------------------------------------
      RETURN
      END

C ##############################################################################
C
C CREATE_IRECTS_TEMP: Copy secondary surface nodes to INTBUF_TAB%IRECTS
C   This creates a temporary array of secondary segment node IDs
C   that will be used for contact detection. The 2D array IGRSURF(1)%NODES
C   (segments × nodes) is flattened into a 1D array INTBUF_TAB(1)%IRECTS.
C
C
      SUBROUTINE CREATE_IRECTS_TEMP(INTBUF_TAB,IGRSURF,NIN,ISURF)
C-----------------------------------------------
C   M o d u l e s
C----------------------------------------------- 
      USE INTBUFDEF_MOD
      USE GROUPDEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      TYPE (SURF_)   , DIMENSION(NSURF)   :: IGRSURF
      INTEGER NIN, ISURF
      ! ISURF is the internal id of the group of secondary nodes.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER num_rows(2)
      INTEGER I, J, K
      INTEGER NSEG, NNOD_PER_SEG
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
C     Safety checks: Validate surface data
      IF (NSURF < 1 .OR. ISURF < 1 .OR. ISURF > NSURF) THEN
          write(*,*) "CREATE_IRECTS_TEMP: ERROR - NSURF < 1"
          RETURN
      ENDIF
      
      IF (IGRSURF(ISURF)%NSEG <= 0) THEN
          write(*,*) "CREATE_IRECTS_TEMP: WARNING - No segments in surface"
          RETURN
      ENDIF
      
      IF (.NOT. ALLOCATED(IGRSURF(ISURF)%NODES)) THEN
          write(*,*) "CREATE_IRECTS_TEMP: ERROR - IGRSURF(1)%NODES not allocated"
          RETURN
      ENDIF

C     Get dimensions of NODES array
      num_rows = SHAPE(IGRSURF(ISURF)%NODES)
      NSEG = num_rows(1)
      NNOD_PER_SEG = num_rows(2)
      
C     Validate array dimensions
      IF (NSEG /= IGRSURF(ISURF)%NSEG) THEN
          write(*,*) "CREATE_IRECTS_TEMP: WARNING - Dimension mismatch: ",
     &              "NSEG = ", IGRSURF(ISURF)%NSEG, ", but array has ", NSEG, " rows"
      ENDIF
      
      IF (NNOD_PER_SEG /= 4) THEN
          write(*,*) "CREATE_IRECTS_TEMP: WARNING - Expected 4 nodes per segment, ",
     &              "but found ", NNOD_PER_SEG
      ENDIF

C     Deallocate existing IRECTS if present
      IF (ALLOCATED(INTBUF_TAB(NIN)%IRECTS)) THEN
          DEALLOCATE(INTBUF_TAB(NIN)%IRECTS)
      END IF

C     Calculate total number of nodes = number of segments × nodes per segment
      INTBUF_TAB(NIN)%S_IRECTS = NSEG * NNOD_PER_SEG
      
C     Allocate IRECTS array
      ALLOCATE(INTBUF_TAB(NIN)%IRECTS(INTBUF_TAB(NIN)%S_IRECTS))

C     Copy nodes from 2D array (segments × nodes) to 1D array
C     Order: Segment 1 (nodes 1-4), Segment 2 (nodes 1-4), ...
      K = 1
      do I = 1, NSEG
          do J = 1, NNOD_PER_SEG
C             Skip zero nodes (if any) - they indicate missing nodes in segment
              IF (IGRSURF(ISURF)%NODES(I,J) > 0) THEN
                  INTBUF_TAB(NIN)%IRECTS(K) = IGRSURF(ISURF)%NODES(I,J)
                  K = K + 1
              ENDIF
          enddo
      enddo
      
C     Update S_IRECTS to actual number of nodes copied (if zeros were skipped)
      IF (K - 1 < INTBUF_TAB(NIN)%S_IRECTS) THEN
          write(*,*) "CREATE_IRECTS_TEMP: WARNING - Some zero nodes skipped. ",
     &              "Expected ", INTBUF_TAB(NIN)%S_IRECTS, " nodes, got ", K-1
          INTBUF_TAB(NIN)%S_IRECTS = K - 1
C         Reallocate to actual size (optional - can keep larger size)
      ENDIF

      RETURN
      END

C ##############################################################################
C
C REMOVE_DUPLI: Remove duplicate node IDs from INTBUF_TAB%IRECTS
C   Creates a unique list of secondary nodes for contact search algorithm.
C   The algorithm uses a mask-based approach to identify duplicates and
C   extracts only unique node IDs into INTBUF_TAB%NSV.
C
      SUBROUTINE REMOVE_DUPLI(INTBUF_TAB,NIN)
C-----------------------------------------------
C   M o d u l e s
C----------------------------------------------- 
      USE INTBUFDEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      INTEGER NIN
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: unique_arr
      LOGICAL, DIMENSION(:), ALLOCATABLE :: mask
      INTEGER I, J, N, SIZ_IRECTS
      INTEGER num_duplicates
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
C     Safety check: IRECTS must be allocated
      IF (.NOT. ALLOCATED(INTBUF_TAB(NIN)%IRECTS)) THEN
          write(*,*) "REMOVE_DUPLI: WARNING - INTBUF_TAB(1)%IRECTS not allocated"
          RETURN
      ENDIF
      
      SIZ_IRECTS = SIZE(INTBUF_TAB(NIN)%IRECTS)
      IF (SIZ_IRECTS <= 0) THEN
          write(*,*) "REMOVE_DUPLI: WARNING - IRECTS array is empty"
          RETURN
      ENDIF
      
C     Validate S_IRECTS matches actual array size
      IF (INTBUF_TAB(NIN)%S_IRECTS /= SIZ_IRECTS) THEN
          write(*,*) "REMOVE_DUPLI: WARNING - Size mismatch: S_IRECTS = ",
     &              INTBUF_TAB(NIN)%S_IRECTS, ", but array size = ", SIZ_IRECTS
      ENDIF
      
C     Allocate working arrays
      ALLOCATE(mask(SIZ_IRECTS))
      
C     Build mask: mark duplicates as .FALSE.
C     First element is always unique (no previous elements to compare)
      mask = .TRUE.
      num_duplicates = 0
      
      do i = 2, SIZ_IRECTS
C         Skip zero or invalid node IDs
          IF (INTBUF_TAB(NIN)%IRECTS(i) <= 0) THEN
              mask(i) = .FALSE.
              num_duplicates = num_duplicates + 1
              CYCLE
          ENDIF
          
C         Check if current element appears in previous elements
          do j = 1, i-1
C             Skip comparison with invalid nodes
              IF (INTBUF_TAB(NIN)%IRECTS(j) <= 0) CYCLE
              
              if (INTBUF_TAB(NIN)%IRECTS(j) == INTBUF_TAB(NIN)%IRECTS(i)) THEN
                  mask(i) = .FALSE.
                  num_duplicates = num_duplicates + 1
                  EXIT  ! Found duplicate, no need to check further
              endif
          enddo
      enddo
      
C     Count unique elements
      n = COUNT(mask)
      
C     Validate that we found at least one unique element
      IF (n <= 0) THEN
          write(*,*) "REMOVE_DUPLI: ERROR - No unique nodes found!"
          IF (ALLOCATED(INTBUF_TAB(NIN)%NSV)) THEN
              DEALLOCATE(INTBUF_TAB(NIN)%NSV)
          ENDIF
          INTBUF_TAB(NIN)%S_NSV = 0
          DEALLOCATE(mask)
          RETURN
      ENDIF
      
C     Allocate unique_arr with exact size needed
      ALLOCATE(unique_arr(n))
      
C     Extract unique elements manually
      J = 1
      do I = 1, SIZ_IRECTS
          if (mask(I)) THEN
C             Additional safety check: skip zero nodes
              IF (INTBUF_TAB(NIN)%IRECTS(I) > 0) THEN
                  unique_arr(J) = INTBUF_TAB(NIN)%IRECTS(I)
                  J = J + 1
              ENDIF
          endif
      enddo
      
C     Update actual count (in case zeros were skipped)
      n = J - 1
      
C     Update NSV array with unique nodes
      IF (n > 0) THEN
C         Deallocate existing NSV if present
          IF (ALLOCATED(INTBUF_TAB(NIN)%NSV)) THEN
              DEALLOCATE(INTBUF_TAB(NIN)%NSV)
          ENDIF
          
          INTBUF_TAB(NIN)%S_NSV = n
          ALLOCATE(INTBUF_TAB(NIN)%NSV(INTBUF_TAB(NIN)%S_NSV))
          INTBUF_TAB(NIN)%NSV(1:n) = unique_arr(1:n)
          
C         Optional: Debug output
          !write(*,*) "REMOVE_DUPLI: Found ", n, " unique nodes out of ", 
          !&          SIZ_IRECTS, " total (", num_duplicates, " duplicates removed)"
      ELSE
C         No valid unique elements found - set empty array
          write(*,*) "REMOVE_DUPLI: WARNING - No valid unique nodes found"
          IF (ALLOCATED(INTBUF_TAB(NIN)%NSV)) THEN
              DEALLOCATE(INTBUF_TAB(NIN)%NSV)
          ENDIF
          INTBUF_TAB(NIN)%S_NSV = 0
      ENDIF
      
C     Clean up working arrays
      IF (ALLOCATED(unique_arr)) DEALLOCATE(unique_arr)
      DEALLOCATE(mask)

      RETURN
      END
