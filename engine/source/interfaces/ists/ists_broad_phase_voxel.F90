!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!||====================================================================
!||    sts_broad_phase_voxel_mod   ../engine/source/interfaces/ists/ists_broad_phase_voxel.F90
!||====================================================================
!
!   STS voxel broad-phase: produce (master_seg, secondary_seg) candidate
!   pairs from two surfaces (IGRSURF entries) without using the legacy
!   INT7 candidate arrays.
!
!   Algorithm mirrors Q1NP's voxel kernel:
!     1) Sample NSAMPLES_PER_SEG points per segment of each surface
!        (4 corners + centroid by default).
!     2) Build a uniform voxel grid over the (padded) master point cloud.
!     3) For every secondary point that lies inside the padded box, query
!        the 3x3x3 cell neighborhood and emit unique (mst_seg, sec_seg)
!        pairs into the existing STS_CONTACTS_ASSEMBLE arrays.
!
!   Output formats match what STS_REMAP_SEGMENTS produces today, so
!   downstream STS_CONTACTS_ASSEMBLE consumes the data unchanged:
!     CAND_SEC_SEG_ID(I, 1)   = secondary segment index in IGRSURF(SEC)
!     CAND_SEC_SEG_ID(I, 2:5) = secondary segment node IDs
!     CAND_MST_SEG_ID(I, 1)   = master   segment index in IGRSURF(MST)
!     CAND_MST_SEG_ID(I, 2:5) = master   segment node IDs
!     CONT_ELEMENT(I, 1:3, 1:4) = master   (primary)   coordinates
!     CONT_ELEMENT(I, 1:3, 5:8) = secondary            coordinates
!
      MODULE STS_BROAD_PHASE_VOXEL_MOD
        USE PRECISION_MOD, ONLY : WP
        USE CONSTANT_MOD,  ONLY : ZERO, ONE, TWO, THREE, FOUR, HALF
        USE GROUPDEF_MOD,  ONLY : SURF_
        IMPLICIT NONE
        PRIVATE
!-----------------------------------------------------------------------
!       Tuning constants (mirroring Q1NP defaults)
!-----------------------------------------------------------------------
!       Trigger tolerance derived from the user GAP value.
        REAL(KIND=WP), PARAMETER, PUBLIC :: STS_VOXEL_TRIGGER_FACTOR   = ONE
!       Voxel cell size = factor * trigger tolerance.
        REAL(KIND=WP), PARAMETER, PUBLIC :: STS_VOXEL_CELL_SIZE_FACTOR = FOUR
!       Search padding around master AABB.
        REAL(KIND=WP), PARAMETER, PUBLIC :: STS_VOXEL_PADDING_FACTOR   = TWO
!       Lower bound used when GAP is zero/too small.
        REAL(KIND=WP), PARAMETER, PUBLIC :: STS_VOXEL_GAP_FALLBACK     = 1.0E-6_WP
!       Default number of samples per segment (4 corners + centroid).
        INTEGER, PARAMETER, PUBLIC :: STS_VOXEL_NSAMPLES_PER_SEG = 5
!
        PUBLIC :: STS_VOXEL_BROAD_PHASE
!
      CONTAINS
!=======================================================================
!   STS_VOXEL_BROAD_PHASE
!
!   Top-level entry: builds segment sample clouds, runs voxel search,
!   emits unique segment pairs and corresponding 8-node coordinate
!   block. Returns COUNT = number of valid pairs (clamped to
!   MAX_STS_SIZE_ACTUAL). OVERFLOW is set to .TRUE. if the storage was
!   saturated and not all candidate pairs could be stored.
!=======================================================================
        SUBROUTINE STS_VOXEL_BROAD_PHASE( &
     &      IGRSURF, NSURF, SEC_SURF_IDX, MST_SURF_IDX, &
     &      X, NUMNOD, GAP, MAX_STS_SIZE_ACTUAL, &
     &      CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &      COUNT, OVERFLOW)
!-----------------------------------------------
!         Dummy arguments
!-----------------------------------------------
          INTEGER, INTENT(IN)  :: NSURF
          TYPE(SURF_), DIMENSION(NSURF), INTENT(IN) :: IGRSURF
          INTEGER, INTENT(IN)  :: SEC_SURF_IDX, MST_SURF_IDX
          INTEGER, INTENT(IN)  :: NUMNOD
          REAL(KIND=WP), INTENT(IN) :: X(3, NUMNOD)
          REAL(KIND=WP), INTENT(IN) :: GAP
          INTEGER, INTENT(IN)  :: MAX_STS_SIZE_ACTUAL
          INTEGER, INTENT(OUT) :: CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          INTEGER, INTENT(OUT) :: CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          REAL(KIND=WP), INTENT(OUT) :: CONT_ELEMENT(MAX_STS_SIZE_ACTUAL, 3, 8)
          INTEGER, INTENT(OUT) :: COUNT
          LOGICAL, INTENT(OUT) :: OVERFLOW
!-----------------------------------------------
!         Local variables
!-----------------------------------------------
          INTEGER :: NSEG_SEC, NSEG_MST
          INTEGER :: NPTS_S, NPTS_M
          REAL(KIND=WP), ALLOCATABLE :: PTS_S(:,:), PTS_M(:,:)
          INTEGER, ALLOCATABLE :: SEG_OF_PT_S(:), SEG_OF_PT_M(:)
          REAL(KIND=WP) :: GAP_EFF, TRIGGER_TOL
!-----------------------------------------------
!         Initialization
!-----------------------------------------------
          COUNT = 0
          OVERFLOW = .FALSE.
          CAND_SEC_SEG_ID(:,:) = 0
          CAND_MST_SEG_ID(:,:) = 0
          CONT_ELEMENT(:,:,:)  = ZERO
!
          IF (SEC_SURF_IDX <= 0 .OR. SEC_SURF_IDX > NSURF) RETURN
          IF (MST_SURF_IDX <= 0 .OR. MST_SURF_IDX > NSURF) RETURN
          IF (MAX_STS_SIZE_ACTUAL <= 0) RETURN
!
          NSEG_SEC = IGRSURF(SEC_SURF_IDX)%NSEG
          NSEG_MST = IGRSURF(MST_SURF_IDX)%NSEG
          IF (NSEG_SEC <= 0 .OR. NSEG_MST <= 0) RETURN
          IF (.NOT. ALLOCATED(IGRSURF(SEC_SURF_IDX)%NODES)) RETURN
          IF (.NOT. ALLOCATED(IGRSURF(MST_SURF_IDX)%NODES)) RETURN
!
          GAP_EFF = MAX(STS_VOXEL_GAP_FALLBACK, ABS(GAP))
          TRIGGER_TOL = MAX(STS_VOXEL_GAP_FALLBACK, &
     &                      STS_VOXEL_TRIGGER_FACTOR * GAP_EFF)
!-----------------------------------------------
!         Build sample point clouds
!-----------------------------------------------
          ALLOCATE(PTS_S(3, NSEG_SEC * STS_VOXEL_NSAMPLES_PER_SEG))
          ALLOCATE(SEG_OF_PT_S(NSEG_SEC * STS_VOXEL_NSAMPLES_PER_SEG))
          CALL STS_VOXEL_BUILD_SEG_POINTS( &
     &        IGRSURF, NSURF, SEC_SURF_IDX, X, NUMNOD, &
     &        STS_VOXEL_NSAMPLES_PER_SEG, &
     &        PTS_S, SEG_OF_PT_S, NPTS_S)
!
          ALLOCATE(PTS_M(3, NSEG_MST * STS_VOXEL_NSAMPLES_PER_SEG))
          ALLOCATE(SEG_OF_PT_M(NSEG_MST * STS_VOXEL_NSAMPLES_PER_SEG))
          CALL STS_VOXEL_BUILD_SEG_POINTS( &
     &        IGRSURF, NSURF, MST_SURF_IDX, X, NUMNOD, &
     &        STS_VOXEL_NSAMPLES_PER_SEG, &
     &        PTS_M, SEG_OF_PT_M, NPTS_M)
!
          IF (NPTS_S <= 0 .OR. NPTS_M <= 0) THEN
            DEALLOCATE(PTS_S, PTS_M, SEG_OF_PT_S, SEG_OF_PT_M)
            RETURN
          END IF
!-----------------------------------------------
!         Voxel pair search and pair emission
!-----------------------------------------------
          CALL STS_VOXEL_PAIR_SEARCH( &
     &        IGRSURF, NSURF, SEC_SURF_IDX, MST_SURF_IDX, X, NUMNOD, &
     &        PTS_S, SEG_OF_PT_S, NPTS_S, &
     &        PTS_M, SEG_OF_PT_M, NPTS_M, &
     &        TRIGGER_TOL, MAX_STS_SIZE_ACTUAL, &
     &        CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &        COUNT, OVERFLOW)
!
          DEALLOCATE(PTS_S, PTS_M, SEG_OF_PT_S, SEG_OF_PT_M)
!
        END SUBROUTINE STS_VOXEL_BROAD_PHASE
!=======================================================================
!   STS_VOXEL_BUILD_SEG_POINTS
!
!   Sample NSAMPLES_PER_SEG points per segment of IGRSURF(SURF_IDX).
!   Layout per segment: 4 corner nodes followed by the centroid, so
!   NSAMPLES_PER_SEG must be either 4 or 5. Any other value falls back
!   to corners-only and the centroid sample is skipped.
!=======================================================================
        SUBROUTINE STS_VOXEL_BUILD_SEG_POINTS( &
     &      IGRSURF, NSURF, SURF_IDX, X, NUMNOD, &
     &      NSAMPLES_PER_SEG, PTS, SEG_OF_PT, NPTS_OUT)
          INTEGER, INTENT(IN)  :: NSURF, SURF_IDX
          TYPE(SURF_), DIMENSION(NSURF), INTENT(IN) :: IGRSURF
          INTEGER, INTENT(IN)  :: NUMNOD, NSAMPLES_PER_SEG
          REAL(KIND=WP), INTENT(IN) :: X(3, NUMNOD)
          REAL(KIND=WP), INTENT(OUT) :: PTS(:, :)
          INTEGER, INTENT(OUT) :: SEG_OF_PT(:)
          INTEGER, INTENT(OUT) :: NPTS_OUT
!
          INTEGER :: ISEG, NSEG, K, NID, NVALID
          REAL(KIND=WP) :: XC, YC, ZC, INV_NV
!
          NPTS_OUT = 0
          NSEG = IGRSURF(SURF_IDX)%NSEG
          IF (NSEG <= 0) RETURN
!
          DO ISEG = 1, NSEG
            XC = ZERO
            YC = ZERO
            ZC = ZERO
            NVALID = 0
!
!           4 corner samples
            DO K = 1, 4
              NID = IGRSURF(SURF_IDX)%NODES(ISEG, K)
              IF (NID <= 0 .OR. NID > NUMNOD) CYCLE
              NPTS_OUT = NPTS_OUT + 1
              PTS(1, NPTS_OUT) = X(1, NID)
              PTS(2, NPTS_OUT) = X(2, NID)
              PTS(3, NPTS_OUT) = X(3, NID)
              SEG_OF_PT(NPTS_OUT) = ISEG
              XC = XC + X(1, NID)
              YC = YC + X(2, NID)
              ZC = ZC + X(3, NID)
              NVALID = NVALID + 1
            END DO
!
!           Centroid sample (only if requested and at least one valid corner).
            IF (NSAMPLES_PER_SEG >= 5 .AND. NVALID > 0) THEN
              INV_NV = ONE / REAL(NVALID, WP)
              NPTS_OUT = NPTS_OUT + 1
              PTS(1, NPTS_OUT) = XC * INV_NV
              PTS(2, NPTS_OUT) = YC * INV_NV
              PTS(3, NPTS_OUT) = ZC * INV_NV
              SEG_OF_PT(NPTS_OUT) = ISEG
            END IF
          END DO
!
        END SUBROUTINE STS_VOXEL_BUILD_SEG_POINTS
!=======================================================================
!   STS_VOXEL_PAIR_SEARCH
!
!   Voxel grid neighborhood search:
!     1) Compute master AABB and pad it by SEARCH_PADDING.
!     2) Coarse AABB-vs-AABB rejection: if minimum distance > padding,
!        return COUNT = 0.
!     3) Build a linked-list voxel grid containing master sample points.
!     4) For each secondary point inside the padded box, look up the
!        3x3x3 cell neighborhood and emit unique (mst_seg, sec_seg)
!        pairs into the output arrays.
!=======================================================================
        SUBROUTINE STS_VOXEL_PAIR_SEARCH( &
     &      IGRSURF, NSURF, SEC_SURF_IDX, MST_SURF_IDX, X, NUMNOD, &
     &      PTS_S, SEG_OF_PT_S, NPTS_S, &
     &      PTS_M, SEG_OF_PT_M, NPTS_M, &
     &      TRIGGER_TOL, MAX_STS_SIZE_ACTUAL, &
     &      CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &      COUNT, OVERFLOW)
          INTEGER, INTENT(IN)  :: NSURF, SEC_SURF_IDX, MST_SURF_IDX
          TYPE(SURF_), DIMENSION(NSURF), INTENT(IN) :: IGRSURF
          INTEGER, INTENT(IN)  :: NUMNOD
          REAL(KIND=WP), INTENT(IN) :: X(3, NUMNOD)
          INTEGER, INTENT(IN)  :: NPTS_S, NPTS_M
          REAL(KIND=WP), INTENT(IN) :: PTS_S(3, NPTS_S), PTS_M(3, NPTS_M)
          INTEGER, INTENT(IN)  :: SEG_OF_PT_S(NPTS_S), SEG_OF_PT_M(NPTS_M)
          REAL(KIND=WP), INTENT(IN) :: TRIGGER_TOL
          INTEGER, INTENT(IN)  :: MAX_STS_SIZE_ACTUAL
          INTEGER, INTENT(INOUT) :: CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          INTEGER, INTENT(INOUT) :: CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          REAL(KIND=WP), INTENT(INOUT) :: CONT_ELEMENT(MAX_STS_SIZE_ACTUAL, 3, 8)
          INTEGER, INTENT(INOUT) :: COUNT
          LOGICAL, INTENT(INOUT) :: OVERFLOW
!
          REAL(KIND=WP) :: XMIN_M, XMAX_M, YMIN_M, YMAX_M, ZMIN_M, ZMAX_M
          REAL(KIND=WP) :: XMIN_S, XMAX_S, YMIN_S, YMAX_S, ZMIN_S, ZMAX_S
          REAL(KIND=WP) :: TOL, CELL_SIZE, SEARCH_PADDING
          REAL(KIND=WP) :: SPAN_X, SPAN_Y, SPAN_Z
          REAL(KIND=WP) :: SPAN_X_SAFE, SPAN_Y_SAFE, SPAN_Z_SAFE
          REAL(KIND=WP) :: DX, DY, DZ, AABB_DIST
          REAL(KIND=WP) :: RX, RY, RZ
          REAL(KIND=WP), PARAMETER :: EPS_SPAN = 1.0E-12_WP
          INTEGER :: NBX, NBY, NBZ, NVOXELS
          INTEGER :: IM, IS, IX, IY, IZ, JX, JY, JZ, CELLID, JJ
          INTEGER :: IX_LO, IX_HI, IY_LO, IY_HI, IZ_LO, IZ_HI
          INTEGER, ALLOCATABLE :: VOXEL(:), NEXT_PT(:)
          INTEGER :: SEG_S, SEG_M
!
!         Build master AABB
          XMIN_M = PTS_M(1, 1); XMAX_M = PTS_M(1, 1)
          YMIN_M = PTS_M(2, 1); YMAX_M = PTS_M(2, 1)
          ZMIN_M = PTS_M(3, 1); ZMAX_M = PTS_M(3, 1)
          DO IM = 2, NPTS_M
            IF (PTS_M(1,IM) < XMIN_M) XMIN_M = PTS_M(1,IM)
            IF (PTS_M(1,IM) > XMAX_M) XMAX_M = PTS_M(1,IM)
            IF (PTS_M(2,IM) < YMIN_M) YMIN_M = PTS_M(2,IM)
            IF (PTS_M(2,IM) > YMAX_M) YMAX_M = PTS_M(2,IM)
            IF (PTS_M(3,IM) < ZMIN_M) ZMIN_M = PTS_M(3,IM)
            IF (PTS_M(3,IM) > ZMAX_M) ZMAX_M = PTS_M(3,IM)
          END DO
!
!         Build secondary AABB
          XMIN_S = PTS_S(1, 1); XMAX_S = PTS_S(1, 1)
          YMIN_S = PTS_S(2, 1); YMAX_S = PTS_S(2, 1)
          ZMIN_S = PTS_S(3, 1); ZMAX_S = PTS_S(3, 1)
          DO IS = 2, NPTS_S
            IF (PTS_S(1,IS) < XMIN_S) XMIN_S = PTS_S(1,IS)
            IF (PTS_S(1,IS) > XMAX_S) XMAX_S = PTS_S(1,IS)
            IF (PTS_S(2,IS) < YMIN_S) YMIN_S = PTS_S(2,IS)
            IF (PTS_S(2,IS) > YMAX_S) YMAX_S = PTS_S(2,IS)
            IF (PTS_S(3,IS) < ZMIN_S) ZMIN_S = PTS_S(3,IS)
            IF (PTS_S(3,IS) > ZMAX_S) ZMAX_S = PTS_S(3,IS)
          END DO
!
          TOL = MAX(EPS_SPAN, TRIGGER_TOL)
          CELL_SIZE      = STS_VOXEL_CELL_SIZE_FACTOR * TOL
          SEARCH_PADDING = STS_VOXEL_PADDING_FACTOR  * TOL
!
!         Coarse AABB rejection (use unpadded master AABB vs secondary AABB)
          DX = MAX(ZERO, MAX(XMIN_S - XMAX_M, XMIN_M - XMAX_S))
          DY = MAX(ZERO, MAX(YMIN_S - YMAX_M, YMIN_M - YMAX_S))
          DZ = MAX(ZERO, MAX(ZMIN_S - ZMAX_M, ZMIN_M - ZMAX_S))
          AABB_DIST = SQRT(DX*DX + DY*DY + DZ*DZ)
          IF (AABB_DIST > SEARCH_PADDING) RETURN
!
!         Pad master AABB
          XMIN_M = XMIN_M - SEARCH_PADDING
          XMAX_M = XMAX_M + SEARCH_PADDING
          YMIN_M = YMIN_M - SEARCH_PADDING
          YMAX_M = YMAX_M + SEARCH_PADDING
          ZMIN_M = ZMIN_M - SEARCH_PADDING
          ZMAX_M = ZMAX_M + SEARCH_PADDING
!
          SPAN_X = XMAX_M - XMIN_M
          SPAN_Y = YMAX_M - YMIN_M
          SPAN_Z = ZMAX_M - ZMIN_M
          SPAN_X_SAFE = MAX(SPAN_X, EPS_SPAN)
          SPAN_Y_SAFE = MAX(SPAN_Y, EPS_SPAN)
          SPAN_Z_SAFE = MAX(SPAN_Z, EPS_SPAN)
!
          NBX = MAX(1, CEILING(SPAN_X / CELL_SIZE))
          NBY = MAX(1, CEILING(SPAN_Y / CELL_SIZE))
          NBZ = MAX(1, CEILING(SPAN_Z / CELL_SIZE))
          NVOXELS = NBX * NBY * NBZ
!
          ALLOCATE(VOXEL(NVOXELS))
          ALLOCATE(NEXT_PT(NPTS_M))
          VOXEL(1:NVOXELS) = 0
          NEXT_PT(1:NPTS_M) = 0
!
!         Fill voxel grid with master points (linked-list per cell)
          DO IM = 1, NPTS_M
            RX = REAL(NBX, WP) * (PTS_M(1,IM) - XMIN_M) / SPAN_X_SAFE
            RY = REAL(NBY, WP) * (PTS_M(2,IM) - YMIN_M) / SPAN_Y_SAFE
            RZ = REAL(NBZ, WP) * (PTS_M(3,IM) - ZMIN_M) / SPAN_Z_SAFE
            IF (RX /= RX) RX = ZERO
            IF (RY /= RY) RY = ZERO
            IF (RZ /= RZ) RZ = ZERO
            IX = INT(RX)
            IY = INT(RY)
            IZ = INT(RZ)
            IX = MAX(0, MIN(NBX - 1, IX))
            IY = MAX(0, MIN(NBY - 1, IY))
            IZ = MAX(0, MIN(NBZ - 1, IZ))
            CELLID = IZ * NBX * NBY + IY * NBX + IX + 1
            CELLID = MAX(1, MIN(NVOXELS, CELLID))
            NEXT_PT(IM) = VOXEL(CELLID)
            VOXEL(CELLID) = IM
          END DO
!
!         Query: for each secondary point inside the padded box, scan the
!         3x3x3 neighborhood and emit unique pairs.
          DO IS = 1, NPTS_S
            IF (PTS_S(1,IS) < XMIN_M) CYCLE
            IF (PTS_S(1,IS) > XMAX_M) CYCLE
            IF (PTS_S(2,IS) < YMIN_M) CYCLE
            IF (PTS_S(2,IS) > YMAX_M) CYCLE
            IF (PTS_S(3,IS) < ZMIN_M) CYCLE
            IF (PTS_S(3,IS) > ZMAX_M) CYCLE
!
            RX = REAL(NBX, WP) * (PTS_S(1,IS) - XMIN_M) / SPAN_X_SAFE
            RY = REAL(NBY, WP) * (PTS_S(2,IS) - YMIN_M) / SPAN_Y_SAFE
            RZ = REAL(NBZ, WP) * (PTS_S(3,IS) - ZMIN_M) / SPAN_Z_SAFE
            IF (RX /= RX) RX = ZERO
            IF (RY /= RY) RY = ZERO
            IF (RZ /= RZ) RZ = ZERO
            IX = INT(RX)
            IY = INT(RY)
            IZ = INT(RZ)
            IX = MAX(0, MIN(NBX - 1, IX))
            IY = MAX(0, MIN(NBY - 1, IY))
            IZ = MAX(0, MIN(NBZ - 1, IZ))
!
            IX_LO = MAX(0, IX - 1)
            IX_HI = MIN(NBX - 1, IX + 1)
            IY_LO = MAX(0, IY - 1)
            IY_HI = MIN(NBY - 1, IY + 1)
            IZ_LO = MAX(0, IZ - 1)
            IZ_HI = MIN(NBZ - 1, IZ + 1)
!
            SEG_S = SEG_OF_PT_S(IS)
!
            DO JZ = IZ_LO, IZ_HI
              DO JY = IY_LO, IY_HI
                DO JX = IX_LO, IX_HI
                  CELLID = JZ * NBX * NBY + JY * NBX + JX + 1
                  CELLID = MAX(1, MIN(NVOXELS, CELLID))
                  JJ = VOXEL(CELLID)
                  DO WHILE (JJ > 0)
                    IF (JJ > NPTS_M) EXIT
                    SEG_M = SEG_OF_PT_M(JJ)
                    IF (SEG_S > 0 .AND. SEG_M > 0) THEN
                      CALL STS_VOXEL_EMIT_PAIR( &
     &                  IGRSURF, NSURF, SEC_SURF_IDX, MST_SURF_IDX, &
     &                  X, NUMNOD, SEG_S, SEG_M, &
     &                  MAX_STS_SIZE_ACTUAL, &
     &                  CAND_SEC_SEG_ID, CAND_MST_SEG_ID, &
     &                  CONT_ELEMENT, COUNT, OVERFLOW)
                      IF (OVERFLOW) EXIT
                    END IF
                    JJ = NEXT_PT(JJ)
                    IF (JJ < 0) EXIT
                  END DO
                  IF (OVERFLOW) EXIT
                END DO
                IF (OVERFLOW) EXIT
              END DO
              IF (OVERFLOW) EXIT
            END DO
            IF (OVERFLOW) EXIT
          END DO
!
          DEALLOCATE(VOXEL, NEXT_PT)
!
        END SUBROUTINE STS_VOXEL_PAIR_SEARCH
!=======================================================================
!   STS_VOXEL_EMIT_PAIR
!
!   Append a single (sec_seg, mst_seg) pair to the output arrays unless
!   it already exists. On overflow (COUNT == MAX_STS_SIZE_ACTUAL) the
!   OVERFLOW flag is set and no further pairs are stored.
!=======================================================================
        SUBROUTINE STS_VOXEL_EMIT_PAIR( &
     &      IGRSURF, NSURF, SEC_SURF_IDX, MST_SURF_IDX, X, NUMNOD, &
     &      SEC_SEG, MST_SEG, MAX_STS_SIZE_ACTUAL, &
     &      CAND_SEC_SEG_ID, CAND_MST_SEG_ID, CONT_ELEMENT, &
     &      COUNT, OVERFLOW)
          INTEGER, INTENT(IN) :: NSURF, SEC_SURF_IDX, MST_SURF_IDX
          TYPE(SURF_), DIMENSION(NSURF), INTENT(IN) :: IGRSURF
          INTEGER, INTENT(IN) :: NUMNOD
          REAL(KIND=WP), INTENT(IN) :: X(3, NUMNOD)
          INTEGER, INTENT(IN) :: SEC_SEG, MST_SEG
          INTEGER, INTENT(IN) :: MAX_STS_SIZE_ACTUAL
          INTEGER, INTENT(INOUT) :: CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          INTEGER, INTENT(INOUT) :: CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL, 5)
          REAL(KIND=WP), INTENT(INOUT) :: CONT_ELEMENT(MAX_STS_SIZE_ACTUAL, 3, 8)
          INTEGER, INTENT(INOUT) :: COUNT
          LOGICAL, INTENT(INOUT) :: OVERFLOW
!
          INTEGER :: K, J, NID
!
          IF (OVERFLOW) RETURN
!
!         Quick duplicate check (linear scan; expected COUNT is small for
!         realistic STS interfaces).
          DO K = 1, COUNT
            IF (CAND_SEC_SEG_ID(K, 1) == SEC_SEG .AND. &
     &          CAND_MST_SEG_ID(K, 1) == MST_SEG) RETURN
          END DO
!
          IF (COUNT >= MAX_STS_SIZE_ACTUAL) THEN
            OVERFLOW = .TRUE.
            RETURN
          END IF
!
          COUNT = COUNT + 1
!
!         Secondary segment metadata
          CAND_SEC_SEG_ID(COUNT, 1) = SEC_SEG
          DO K = 1, 4
            CAND_SEC_SEG_ID(COUNT, K + 1) = &
     &          IGRSURF(SEC_SURF_IDX)%NODES(SEC_SEG, K)
          END DO
!
!         Master segment metadata
          CAND_MST_SEG_ID(COUNT, 1) = MST_SEG
          DO K = 1, 4
            CAND_MST_SEG_ID(COUNT, K + 1) = &
     &          IGRSURF(MST_SURF_IDX)%NODES(MST_SEG, K)
          END DO
!
!         Master (primary) coordinates -> CONT_ELEMENT(I, 1:3, 1:4)
          J = 1
          DO K = 2, 5
            NID = CAND_MST_SEG_ID(COUNT, K)
            IF (NID > 0 .AND. NID <= NUMNOD) THEN
              CONT_ELEMENT(COUNT, 1, J) = X(1, NID)
              CONT_ELEMENT(COUNT, 2, J) = X(2, NID)
              CONT_ELEMENT(COUNT, 3, J) = X(3, NID)
            ELSE
              CONT_ELEMENT(COUNT, 1:3, J) = ZERO
            END IF
            J = J + 1
          END DO
!
!         Secondary coordinates -> CONT_ELEMENT(I, 1:3, 5:8)
          J = 5
          DO K = 2, 5
            NID = CAND_SEC_SEG_ID(COUNT, K)
            IF (NID > 0 .AND. NID <= NUMNOD) THEN
              CONT_ELEMENT(COUNT, 1, J) = X(1, NID)
              CONT_ELEMENT(COUNT, 2, J) = X(2, NID)
              CONT_ELEMENT(COUNT, 3, J) = X(3, NID)
            ELSE
              CONT_ELEMENT(COUNT, 1:3, J) = ZERO
            END IF
            J = J + 1
          END DO
!
        END SUBROUTINE STS_VOXEL_EMIT_PAIR
!
      END MODULE STS_BROAD_PHASE_VOXEL_MOD
