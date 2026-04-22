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
!Copyright>        Commercial Alternative: Altair Radioss
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!||====================================================================
!||    q1np_contact_algorithms           ../engine/source/interfaces/int26/q1np_contact_algorithms.F90
!||--- called by ------------------------------------------------------
!||    resol                             ../engine/source/engine/resol.F
!||--- calls      -----------------------------------------------------
!||    q1np_evaluate_nurbs_top_surface_point  ../engine/source/elements/solid/solid_q1np/q1np_nurbs_surface_eval_mod.F90
!||--- uses       -----------------------------------------------------
!||    precision_mod                     ../common_source/modules/precision_mod.F
!||    constant_mod                      ../common_source/modules/constant_mod.F
!||    q1np_restart_mod                  ../common_source/modules/q1np_restart_mod.F90
!||    q1np_nurbs_surface_eval_mod       ../engine/source/elements/solid/solid_q1np/q1np_nurbs_surface_eval_mod.F90
!||====================================================================
!
!   Q1NP NURBS contact: broad phase, narrow phase, and penalty forces.
!
!   Builds point clouds on the NURBS top surfaces of two Q1NP element
!   sets (identified by knot_set_id = 1 and 2), performs voxel-based
!   proximity detection, Newton projection to find penetrating pairs,
!   and scatters normal penalty forces onto NURBS control-point nodes.
!
      MODULE Q1NP_CONTACT_ALGORITHMS_MOD
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
        USE PRECISION_MOD, ONLY : WP
        USE CONSTANT_MOD , ONLY : ZERO, ONE, TWO, THREE, TEN, HALF
        USE Q1NP_RESTART_MOD
        USE RESTMOD, ONLY : KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB
        USE Q1NP_NURBS_SURFACE_EVALUATION_MOD, ONLY : &
     &      Q1NP_EVALUATE_NURBS_TOP_SURFACE_POINT, &
     &      Q1NP_EVALUATE_NURBS_TOP_SURFACE_POINT_AND_DERIVS, &
     &      Q1NP_EVALUATE_NURBS_SHAPE_VALUES
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
        IMPLICIT NONE
        PRIVATE
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Tuning parameters
! ----------------------------------------------------------------------------------------------------------------------
!       Broad-phase trigger tolerance (voxel-search).
!       Keep this somewhat larger than the penalty reference gap so
!       contact candidates are still tracked once noticeable overlap starts.
        REAL(KIND=WP), PARAMETER, PUBLIC :: Q1NP_CONTACT_VOXEL_TRIGGER_TOLERANCE = 0.03_WP
!       Keep voxel cells wider than the minimum 2*tolerance rule so the
!       broad phase remains robust when the sampled B-cloud is sparse.
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_VOXEL_CELL_SIZE_FACTOR = 4.0_WP
!       Conservative cap for voxel-derived A-candidates stored per B point.
!       If the local neighborhood exceeds this cap, narrow phase falls back
!       to the original global nearest search for that B point.
        INTEGER, PARAMETER :: Q1NP_CONTACT_MAX_CANDIDATES_PER_B = 32
        LOGICAL, PARAMETER :: Q1NP_CONTACT_DEBUG_CANDIDATES = .FALSE.

!       Number of sample points per parametric direction on each element surface.
!       These samples are used as interior quadrature-style points, not
!       nodal/end-point samples, to avoid duplicate contact pairs on
!       shared element edges and corners.
        INTEGER, PARAMETER :: Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_U = 3
        INTEGER, PARAMETER :: Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_V = 3

!       STS-style penalty tuning for NURBS contact:
!       pair stiffness is derived from both sides' nodal STIFN, then
!       blended, clamped, and scaled by a gap factor.
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_PENALTY_SCALE = 1.0_WP
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_PENALTY_KMIN  = 1.0E+2_WP
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_PENALTY_KMAX  = 1.0E+4_WP
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_GAP_REFERENCE = 0.01_WP
!       Additional stabilizers for contact law.
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_FAC_MAX = 3.0_WP
!       Effective penetration clipping ratio relative to GAP_REFERENCE.
!       Set <= 0 to disable clipping.
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_PENETRATION_RATIO_MAX = -ONE
!       Optional force cap per contact pair (set <= 0 to disable).
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_FORCE_CAP = -ONE
!       Broad-phase scheduling: keep disabled by default while tuning
!       contact response to avoid delayed re-activation.
        LOGICAL, PARAMETER :: Q1NP_CONTACT_ENABLE_ADAPTIVE_SKIP = .TRUE.
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_SKIP_SCALE = 8.0_WP
        REAL(KIND=WP), PARAMETER :: Q1NP_CONTACT_SKIP_EXPONENT = 1.5_WP
        INTEGER, PARAMETER :: Q1NP_CONTACT_SKIP_MAX = 200

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Contact pair type
! ----------------------------------------------------------------------------------------------------------------------
        TYPE Q1NP_CONTACT_PAIR
          REAL(KIND=WP) :: PENETRATION
          REAL(KIND=WP) :: NORMAL(3)
          REAL(KIND=WP) :: XI_PROJ, ETA_PROJ
          INTEGER :: ELEM_A
          REAL(KIND=WP) :: XI_SRC, ETA_SRC
          INTEGER :: ELEM_B
        END TYPE Q1NP_CONTACT_PAIR

        PUBLIC :: Q1NP_CONTACT_BROAD_PHASE_CHECK_PROXIMITY

!       Adaptive call-skipping state:
!       number of future cycles to skip before running the voxel search again.
        INTEGER, SAVE :: Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING = 0

      CONTAINS

!=======================================================================
!   Q1NP_CONTACT_BROAD_PHASE_CHECK_PROXIMITY
!
!   Top-level entry point from Q1NP_CONTACT_DRIVER_INT7 (Type 7 / future INT26).
!   Builds both surface point clouds, computes minimum distance via
!   voxel search, runs narrow-phase projection to collect penetrating
!   pairs, computes penalty forces, and scatters them to A / STIFN.
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_BROAD_PHASE_CHECK_PROXIMITY( &
     &      KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, X_COORDS, NUMNOD, &
     &      NUMELQ1NP, TT, A, STIFN, PROXIMITY_DETECTED)
!C----------------------------------------------------------------------
!C   D u m m y   A r g u m e n t s
!C----------------------------------------------------------------------
          INTEGER, INTENT(IN) :: KQ1NP_TAB(:,:)
          INTEGER, INTENT(IN) :: IQ1NP_TAB(:)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_KTAB(:)
          INTEGER, INTENT(IN) :: NUMNOD, NUMELQ1NP
          REAL(KIND=WP), INTENT(IN) :: TT
          REAL(KIND=WP), INTENT(IN) :: X_COORDS(3,NUMNOD)
          REAL(KIND=WP), INTENT(INOUT) :: A(3,NUMNOD)
          REAL(KIND=WP), INTENT(INOUT) :: STIFN(NUMNOD)
          LOGICAL, INTENT(OUT) :: PROXIMITY_DETECTED
!C----------------------------------------------------------------------
!C   L o c a l   V a r i a b l e s
!C----------------------------------------------------------------------
          INTEGER :: NPTS_A, NPTS_B, IDX_A, IDX_B
          INTEGER :: MAX_PTS, N_PAIRS
          REAL(KIND=WP) :: D_MIN, DIST_RATIO
          REAL(KIND=WP), ALLOCATABLE :: SURF_POINTS_A(:,:)
          REAL(KIND=WP), ALLOCATABLE :: SURF_POINTS_B(:,:)
          INTEGER, ALLOCATABLE :: ELEM_IDS_A(:), ELEM_IDS_B(:)
          REAL(KIND=WP), ALLOCATABLE :: XI_A(:), ETA_A(:)
          REAL(KIND=WP), ALLOCATABLE :: XI_B(:), ETA_B(:)
          TYPE(Q1NP_CONTACT_PAIR), ALLOCATABLE :: PAIRS(:)
          INTEGER, ALLOCATABLE :: CANDIDATE_IA(:,:)
          INTEGER, ALLOCATABLE :: CANDIDATE_COUNT(:)
          LOGICAL, ALLOCATABLE :: CANDIDATE_OVERFLOW(:)
          INTEGER :: ICAND, IB_DBG

          PROXIMITY_DETECTED = .FALSE.

          IF (NUMELQ1NP <= 0) RETURN
          IF (Q1NP_NKNOT_SETS_G < 2) RETURN

          IF (Q1NP_CONTACT_ENABLE_ADAPTIVE_SKIP .AND. &
     &        Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING > 0) THEN
            Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING = &
     &          Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING - 1
            RETURN
          END IF

!  ----------------------------------------------------------------------------------------------------------------------
!                                                BUILD POINT CLOUDS
!  ----------------------------------------------------------------------------------------------------------------------
          MAX_PTS = NUMELQ1NP * Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_U &
     &                        * Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_V
          ALLOCATE(SURF_POINTS_A(3, MAX_PTS))
          ALLOCATE(SURF_POINTS_B(3, MAX_PTS))
          ALLOCATE(ELEM_IDS_A(MAX_PTS), XI_A(MAX_PTS), ETA_A(MAX_PTS))
          ALLOCATE(ELEM_IDS_B(MAX_PTS), XI_B(MAX_PTS), ETA_B(MAX_PTS))

          ! Build point cloud for surface A
          CALL Q1NP_CONTACT_BROAD_PHASE_BUILD_SURFACE_POINTS( &
     &        KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, NUMELQ1NP,    &
     &        X_COORDS, NUMNOD, 1,                             &
     &        Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_U,            &
     &        Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_V,            &
     &        SURF_POINTS_A, NPTS_A,                           &
     &        ELEM_IDS_A, XI_A, ETA_A)
          ! Build point cloud for surface B
          CALL Q1NP_CONTACT_BROAD_PHASE_BUILD_SURFACE_POINTS( &
     &        KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, NUMELQ1NP,    &
     &        X_COORDS, NUMNOD, 2,                             &
     &        Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_U,            &
     &        Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_V,            &
     &        SURF_POINTS_B, NPTS_B,                           &
     &        ELEM_IDS_B, XI_B, ETA_B)

          ALLOCATE(PAIRS(MAX(1, NPTS_B)))
          ALLOCATE(CANDIDATE_IA(Q1NP_CONTACT_MAX_CANDIDATES_PER_B, MAX(1, NPTS_B)))
          ALLOCATE(CANDIDATE_COUNT(MAX(1, NPTS_B)))
          ALLOCATE(CANDIDATE_OVERFLOW(MAX(1, NPTS_B)))
          CANDIDATE_IA(:,:) = 0
          CANDIDATE_COUNT(:) = 0
          CANDIDATE_OVERFLOW(:) = .FALSE.

!  ----------------------------------------------------------------------------------------------------------------------
!                                                VOXEL SEARCH
!  ----------------------------------------------------------------------------------------------------------------------
          CALL Q1NP_CONTACT_BROAD_PHASE_VOXEL_MIN_DISTANCE( &
     &        SURF_POINTS_A, NPTS_A, SURF_POINTS_B, NPTS_B,  &
     &        D_MIN, IDX_A, IDX_B,                           &
     &        CANDIDATE_IA, CANDIDATE_COUNT, CANDIDATE_OVERFLOW)
          ! Returns D_MIN, IDX_A, IDX_B
          ! D_MIN: minimum distance between the two point clouds
          ! CANDIDATE_IA: voxel-derived A-candidates for B
          ! CANDIDATE_COUNT: number of stored candidates for B
          ! CANDIDATE_OVERFLOW: candidate list exceeded storage
!   ----------------------------------------------------------------------------------------------------------------------
!                                       NARROW PHASE: NURBS->NURBS PROJECTION
!   ----------------------------------------------------------------------------------------------------------------------
          N_PAIRS = 0
          IF (IDX_A > 0 .AND. IDX_B > 0) THEN
            CALL Q1NP_CONTACT_NARROW_PHASE_PROJECT( &
     &          KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, X_COORDS, &
     &          NUMNOD, &
     &          SURF_POINTS_A, NPTS_A,                       &
     &          ELEM_IDS_A, XI_A, ETA_A,                     &
     &          SURF_POINTS_B, NPTS_B,                       &
     &          ELEM_IDS_B, XI_B, ETA_B,                     &
     &          CANDIDATE_IA, CANDIDATE_COUNT,               &
     &          CANDIDATE_OVERFLOW,                          &
     &          PAIRS, N_PAIRS, PROXIMITY_DETECTED)
          END IF
          
!   ----------------------------------------------------------------------------------------------------------------------
!                                       PENALTY FORCE CALCULATION AND ASSEMBLY
!   ----------------------------------------------------------------------------------------------------------------------
          IF (N_PAIRS > 0) THEN
            CALL Q1NP_CONTACT_COMPUTE_PENALTY_FORCES( &
     &          PAIRS, N_PAIRS, &
     &          KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &          X_COORDS, NUMNOD, A, STIFN)
          END IF
!   ----------------------------------------------------------------------------------------------------------------------
!                                                ADAPTIVE SKIP SCHEDULING
!   ----------------------------------------------------------------------------------------------------------------------
          IF (Q1NP_CONTACT_ENABLE_ADAPTIVE_SKIP) THEN
            IF (.NOT. PROXIMITY_DETECTED) THEN
              DIST_RATIO = D_MIN / Q1NP_CONTACT_VOXEL_TRIGGER_TOLERANCE
              IF (DIST_RATIO <= ONE) THEN
                Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING = 0
              ELSE
                Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING = MIN(Q1NP_CONTACT_SKIP_MAX, &
     &            INT(Q1NP_CONTACT_SKIP_SCALE * (DIST_RATIO - ONE)**Q1NP_CONTACT_SKIP_EXPONENT))
              END IF
            ELSE
              Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING = 0
            END IF
          ELSE
            Q1NP_CONTACT_BROAD_PHASE_SKIP_REMAINING = 0
          END IF

          DEALLOCATE(SURF_POINTS_A, SURF_POINTS_B)
          DEALLOCATE(ELEM_IDS_A, ELEM_IDS_B)
          DEALLOCATE(XI_A, ETA_A, XI_B, ETA_B)
          DEALLOCATE(CANDIDATE_IA, CANDIDATE_COUNT, CANDIDATE_OVERFLOW)
          DEALLOCATE(PAIRS)

        END SUBROUTINE Q1NP_CONTACT_BROAD_PHASE_CHECK_PROXIMITY

!=======================================================================
!   Q1NP_CONTACT_BROAD_PHASE_BUILD_SURFACE_POINTS
!
!   Build a point cloud on the NURBS top surface of all Q1NP elements
!   belonging to a given knot set (i.e. one of the two contact surfaces).
!
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_BROAD_PHASE_BUILD_SURFACE_POINTS( &
     &      KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, NUMELQ1NP,         &
     &      X_COORDS, NUMNOD, KNOT_SET_ID_FILTER,                 &
     &      NGP_SURF_U, NGP_SURF_V, SURF_POINTS, NPTS_TOTAL,     &
     &      PT_ELEM_IDS, PT_XI, PT_ETA)
!C----------------------------------------------------------------------
!C   D u m m y   A r g u m e n t s
!C----------------------------------------------------------------------
!C     SURF_POINTS(3,:)        - output point cloud coordinates                 (OUT)
!C     NPTS_TOTAL              - number of points actually written              (OUT)
!C     PT_ELEM_IDS(optional)   - element index for each point (optional)        (OUT)
!C     PT_XI(optional)         - xi  parametric coord for each point (optional) (OUT)
!C     PT_ETA(optional)        - eta parametric coord for each point (optional) (OUT)
!C----------------------------------------------------------------------
          INTEGER, INTENT(IN) :: KQ1NP_TAB(:,:)
          INTEGER, INTENT(IN) :: IQ1NP_TAB(:)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_KTAB(:)
          INTEGER, INTENT(IN) :: NUMELQ1NP, NUMNOD
          INTEGER, INTENT(IN) :: KNOT_SET_ID_FILTER
          INTEGER, INTENT(IN) :: NGP_SURF_U, NGP_SURF_V
          REAL(KIND=WP), INTENT(IN)  :: X_COORDS(3,NUMNOD)
          REAL(KIND=WP), INTENT(OUT) :: SURF_POINTS(3,NPTS_TOTAL)
          INTEGER, INTENT(OUT) :: NPTS_TOTAL
          INTEGER, INTENT(OUT), OPTIONAL :: PT_ELEM_IDS(:)
          REAL(KIND=WP), INTENT(OUT), OPTIONAL :: PT_XI(:), PT_ETA(:)
!C----------------------------------------------------------------------
!C   L o c a l   V a r i a b l e s
!C----------------------------------------------------------------------
          INTEGER :: IEL, IU, IV, KNOT_SET_ID
          INTEGER :: P_CUR, Q_CUR, NCTRL, CP_OFFSET
          INTEGER :: ELEM_U_IDX, ELEM_V_IDX
          INTEGER :: NX_CUR, NY_CUR, U_LEN, V_LEN
          INTEGER :: IPT
          INTEGER, PARAMETER :: MAX_CTRL = 50
          INTEGER :: CTRL_IDS(MAX_CTRL)
          REAL(KIND=WP) :: XI_SAMPLE, ETA_SAMPLE
          REAL(KIND=WP) :: XYZ(3)
          REAL(KIND=WP), ALLOCATABLE :: U_KNOT_LOCAL(:), V_KNOT_LOCAL(:)

          NPTS_TOTAL = 0

          DO IEL = 1, NUMELQ1NP
!           KQ1NP_TAB(15,:) stores the knot_set_id assigned by the
!           starter when GENQ1NP was called for surface A (=1) or B (=2).
            KNOT_SET_ID = KQ1NP_TAB(15, IEL)
            IF (KNOT_SET_ID /= KNOT_SET_ID_FILTER) CYCLE

            P_CUR      = KQ1NP_TAB(8, IEL)
            Q_CUR      = KQ1NP_TAB(9, IEL)
            NCTRL      = KQ1NP_TAB(3, IEL)
            CP_OFFSET  = KQ1NP_TAB(4, IEL)
            ELEM_U_IDX = KQ1NP_TAB(6, IEL)
            ELEM_V_IDX = KQ1NP_TAB(7, IEL)

!           Derive grid dimensions for this element's knot set
            NX_CUR = KQ1NP_TAB(12, IEL)
            NY_CUR = KQ1NP_TAB(13, IEL)
            IF (NX_CUR <= 0 .OR. NY_CUR <= 0) THEN
              IF (Q1NP_NKNOT_SETS_G > 0 .AND. &
     &            KNOT_SET_ID > 0 .AND. &
     &            KNOT_SET_ID <= Q1NP_NKNOT_SETS_G) THEN
                NX_CUR = Q1NP_NX_SET_G(KNOT_SET_ID)
                NY_CUR = Q1NP_NY_SET_G(KNOT_SET_ID)
              ELSE
                NX_CUR = Q1NP_NX_G
                NY_CUR = Q1NP_NY_G
              END IF
            END IF

!           Knot vector lengths: n_knots = n_elements + 2*degree + 1
            U_LEN = NX_CUR + 2*P_CUR + 1
            V_LEN = NY_CUR + 2*Q_CUR + 1
            ALLOCATE(U_KNOT_LOCAL(U_LEN), V_KNOT_LOCAL(V_LEN))

!           Extract knot vectors from the RESTMOD pool (Q1NP_KTAB).
!           Same logic as q1np_forc3 lines 478-491.
            IF (Q1NP_NKNOT_SETS_G > 0 .AND. &
     &          KNOT_SET_ID > 0 .AND. &
     &          KNOT_SET_ID <= Q1NP_NKNOT_SETS_G .AND. &
     &          ALLOCATED(Q1NP_KTAB_OFF_G)) THEN
              U_KNOT_LOCAL(1:U_LEN) = &
     &            Q1NP_KTAB(Q1NP_KTAB_OFF_G(KNOT_SET_ID) : &
     &                      Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN - 1)
              V_KNOT_LOCAL(1:V_LEN) = &
     &            Q1NP_KTAB(Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN : &
     &                      Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN + V_LEN - 1)
            ELSE
              U_KNOT_LOCAL(1:U_LEN) = Q1NP_KTAB(1:U_LEN)
              V_KNOT_LOCAL(1:V_LEN) = Q1NP_KTAB(U_LEN+1:U_LEN+V_LEN)
            END IF

!           Gather control-point IDs for this element
            DO IPT = 1, MIN(NCTRL, MAX_CTRL)
              CTRL_IDS(IPT) = IQ1NP_TAB(CP_OFFSET + IPT - 1)
            END DO

            DO IV = 1, NGP_SURF_V
!             Gauss-Legendre abscissae in [-1, +1] for eta direction.
              SELECT CASE (NGP_SURF_V)
              CASE (1)
                ETA_SAMPLE = ZERO
              CASE (2)
                IF (IV == 1) THEN
                  ETA_SAMPLE = -ONE / SQRT(THREE)
                ELSE
                  ETA_SAMPLE =  ONE / SQRT(THREE)
                END IF
              CASE (3)
                IF (IV == 1) THEN
                  ETA_SAMPLE = -SQRT(3.0_WP/5.0_WP)
                ELSEIF (IV == 2) THEN
                  ETA_SAMPLE = ZERO
                ELSE
                  ETA_SAMPLE =  SQRT(3.0_WP/5.0_WP)
                END IF
              CASE DEFAULT
!               Fallback for unsupported orders: uniform midpoint sampling.
                ETA_SAMPLE = -ONE + TWO * (REAL(IV, WP) - 0.5_WP) / &
     &                       REAL(NGP_SURF_V, WP)
              END SELECT

              DO IU = 1, NGP_SURF_U
!               Gauss-Legendre abscissae in [-1, +1] for xi direction.
                SELECT CASE (NGP_SURF_U)
                CASE (1)
                  XI_SAMPLE = ZERO
                CASE (2)
                  IF (IU == 1) THEN
                    XI_SAMPLE = -ONE / SQRT(THREE)
                  ELSE
                    XI_SAMPLE =  ONE / SQRT(THREE)
                  END IF
                CASE (3)
                  IF (IU == 1) THEN
                    XI_SAMPLE = -SQRT(3.0_WP/5.0_WP)
                  ELSEIF (IU == 2) THEN
                    XI_SAMPLE = ZERO
                  ELSE
                    XI_SAMPLE =  SQRT(3.0_WP/5.0_WP)
                  END IF
                CASE DEFAULT
!                 Fallback for unsupported orders: uniform midpoint sampling.
                  XI_SAMPLE = -ONE + TWO * (REAL(IU, WP) - 0.5_WP) / &
     &                        REAL(NGP_SURF_U, WP)
                END SELECT

                CALL Q1NP_EVALUATE_NURBS_TOP_SURFACE_POINT( &
     &              XI_SAMPLE, ETA_SAMPLE,                   &
     &              P_CUR, Q_CUR,                            &
     &              ELEM_U_IDX, ELEM_V_IDX,                  &
     &              U_KNOT_LOCAL, V_KNOT_LOCAL,              &
     &              MIN(NCTRL, MAX_CTRL), CTRL_IDS,          &
     &              X_COORDS, NUMNOD, XYZ)

                NPTS_TOTAL = NPTS_TOTAL + 1
                SURF_POINTS(1:3, NPTS_TOTAL) = XYZ(1:3)
                IF (PRESENT(PT_ELEM_IDS)) PT_ELEM_IDS(NPTS_TOTAL) = IEL
                IF (PRESENT(PT_XI))       PT_XI(NPTS_TOTAL)  = XI_SAMPLE
                IF (PRESENT(PT_ETA))      PT_ETA(NPTS_TOTAL) = ETA_SAMPLE
              END DO
            END DO

            DEALLOCATE(U_KNOT_LOCAL, V_KNOT_LOCAL)
          END DO

        END SUBROUTINE Q1NP_CONTACT_BROAD_PHASE_BUILD_SURFACE_POINTS

!=======================================================================
!   Q1NP_CONTACT_BROAD_PHASE_VOXEL_MIN_DISTANCE
!
!   Voxel minimum Euclidean distance between two 3D point clouds.
!
!   Algorithm:
!     1. Compute axis-aligned bounding box (AABB) of point cloud A and B.
!     2. Choose a conservative B-box padding and voxel size based on the trigger tolerance.
!     3. Compute a coarse AABB-to-AABB distance and return early when
!        the point clouds are clearly farther apart than the padded
!        search window.
!     4. Build a uniform B-point voxel grid over the padded box.
!     5. For every A-point that lies inside that grid,
!        look up the 3x3x3 cell neighborhood and find the closest B-point.
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_BROAD_PHASE_VOXEL_MIN_DISTANCE( &
     &      SURF_POINTS_A, NPTS_A, SURF_POINTS_B, NPTS_B,       &
     &      D_MIN_OUT, IDX_A_OUT, IDX_B_OUT,                    &
     &      CANDIDATE_IA, CANDIDATE_COUNT, CANDIDATE_OVERFLOW)
!C----------------------------------------------------------------------
!C   D u m m y   A r g u m e n t s
!C----------------------------------------------------------------------
!C     SURF_POINTS_A(3,NPTS_A)  - point cloud of surface A              (IN)
!C     NPTS_A                   - number of points in cloud A           (IN)
!C     SURF_POINTS_B(3,NPTS_B)  - point cloud of surface B              (IN)
!C     NPTS_B                   - number of points in cloud B           (IN)
!C     D_MIN_OUT                - minimum Euclidean distance            (OUT)
!C     IDX_A_OUT                - index in A of the closest pair        (OUT)
!C     IDX_B_OUT                - index in B of the closest pair        (OUT)
!C     CANDIDATE_IA(:,IB)       - voxel-derived A-candidates for B      (OUT)
!C     CANDIDATE_COUNT(IB)      - number of stored candidates for B     (OUT)
!C     CANDIDATE_OVERFLOW(IB)   - candidate list exceeded storage       (OUT)
!C----------------------------------------------------------------------
          INTEGER, INTENT(IN)  :: NPTS_A, NPTS_B
          REAL(KIND=WP), INTENT(IN)  :: SURF_POINTS_A(3, NPTS_A)
          REAL(KIND=WP), INTENT(IN)  :: SURF_POINTS_B(3, NPTS_B)
          REAL(KIND=WP), INTENT(OUT) :: D_MIN_OUT
          INTEGER, INTENT(OUT) :: IDX_A_OUT, IDX_B_OUT
          INTEGER, INTENT(OUT) :: CANDIDATE_IA(:,:)
          INTEGER, INTENT(OUT) :: CANDIDATE_COUNT(:)
          LOGICAL, INTENT(OUT) :: CANDIDATE_OVERFLOW(:)
!C----------------------------------------------------------------------
!C   L o c a l   V a r i a b l e s
!C----------------------------------------------------------------------
          REAL(KIND=WP) :: XMIN_A, XMAX_A, YMIN_A, YMAX_A
          REAL(KIND=WP) :: ZMIN_A, ZMAX_A
          REAL(KIND=WP) :: XMIN_B, XMAX_B, YMIN_B, YMAX_B
          REAL(KIND=WP) :: ZMIN_B, ZMAX_B
          REAL(KIND=WP) :: TOL, CELL_SIZE, SEARCH_PADDING
          REAL(KIND=WP) :: SPAN_X, SPAN_Y, SPAN_Z
          INTEGER :: NBX, NBY, NBZ, NVOXELS
          INTEGER, ALLOCATABLE :: VOXEL(:), NEXT_PT(:)
          INTEGER :: IB, IA, IX, IY, IZ, CELLID, JJ
          INTEGER :: IX_LO, IX_HI, IY_LO, IY_HI, IZ_LO, IZ_HI
          INTEGER :: JX, JY, JZ
          REAL(KIND=WP) :: DX, DY, DZ, DIST_SQ, D_MIN_SQ
          REAL(KIND=WP), PARAMETER :: EPS_SPAN = 1.0E-12_WP
          REAL(KIND=WP) :: SPAN_X_SAFE, SPAN_Y_SAFE, SPAN_Z_SAFE
          REAL(KIND=WP) :: RX, RY, RZ
          INTEGER :: NCAND

          D_MIN_SQ  = HUGE(ONE)
          IDX_A_OUT = 0
          IDX_B_OUT = 0
          TOL       = Q1NP_CONTACT_VOXEL_TRIGGER_TOLERANCE
          CANDIDATE_IA(:,:) = 0
          CANDIDATE_COUNT(:) = 0
          CANDIDATE_OVERFLOW(:) = .FALSE.

!   ----- Step 1a: bounding box of A points -----
          XMIN_A = SURF_POINTS_A(1, 1)
          XMAX_A = SURF_POINTS_A(1, 1)
          YMIN_A = SURF_POINTS_A(2, 1)
          YMAX_A = SURF_POINTS_A(2, 1)
          ZMIN_A = SURF_POINTS_A(3, 1)
          ZMAX_A = SURF_POINTS_A(3, 1)
          DO IA = 2, NPTS_A
            IF (SURF_POINTS_A(1,IA) < XMIN_A) XMIN_A = SURF_POINTS_A(1,IA)
            IF (SURF_POINTS_A(1,IA) > XMAX_A) XMAX_A = SURF_POINTS_A(1,IA)
            IF (SURF_POINTS_A(2,IA) < YMIN_A) YMIN_A = SURF_POINTS_A(2,IA)
            IF (SURF_POINTS_A(2,IA) > YMAX_A) YMAX_A = SURF_POINTS_A(2,IA)
            IF (SURF_POINTS_A(3,IA) < ZMIN_A) ZMIN_A = SURF_POINTS_A(3,IA)
            IF (SURF_POINTS_A(3,IA) > ZMAX_A) ZMAX_A = SURF_POINTS_A(3,IA)
          END DO

!   ----- Step 1b: bounding box of B points -----
          XMIN_B = SURF_POINTS_B(1, 1)
          XMAX_B = SURF_POINTS_B(1, 1)
          YMIN_B = SURF_POINTS_B(2, 1)
          YMAX_B = SURF_POINTS_B(2, 1)
          ZMIN_B = SURF_POINTS_B(3, 1)
          ZMAX_B = SURF_POINTS_B(3, 1)
          DO IB = 2, NPTS_B
            IF (SURF_POINTS_B(1,IB) < XMIN_B) XMIN_B = SURF_POINTS_B(1,IB)
            IF (SURF_POINTS_B(1,IB) > XMAX_B) XMAX_B = SURF_POINTS_B(1,IB)
            IF (SURF_POINTS_B(2,IB) < YMIN_B) YMIN_B = SURF_POINTS_B(2,IB)
            IF (SURF_POINTS_B(2,IB) > YMAX_B) YMAX_B = SURF_POINTS_B(2,IB)
            IF (SURF_POINTS_B(3,IB) < ZMIN_B) ZMIN_B = SURF_POINTS_B(3,IB)
            IF (SURF_POINTS_B(3,IB) > ZMAX_B) ZMAX_B = SURF_POINTS_B(3,IB)
          END DO
!   ----- Step 2: choose padding / cell size from the trigger tolerance -----
          CELL_SIZE = 4 * TOL
          SEARCH_PADDING = 2 * TOL
          ! Pad the bounding box of point cloud B by the search padding
          XMIN_B = XMIN_B - SEARCH_PADDING
          XMAX_B = XMAX_B + SEARCH_PADDING
          YMIN_B = YMIN_B - SEARCH_PADDING
          YMAX_B = YMAX_B + SEARCH_PADDING
          ZMIN_B = ZMIN_B - SEARCH_PADDING
          ZMAX_B = ZMAX_B + SEARCH_PADDING

!   ----- Step 3: coarse AABB rejection before voxel allocation -----
          DX = MAX(ZERO, MAX(XMIN_A - (XMAX_B - SEARCH_PADDING), &
     &                       (XMIN_B + SEARCH_PADDING) - XMAX_A))
          DY = MAX(ZERO, MAX(YMIN_A - (YMAX_B - SEARCH_PADDING), &
     &                       (YMIN_B + SEARCH_PADDING) - YMAX_A))
          DZ = MAX(ZERO, MAX(ZMIN_A - (ZMAX_B - SEARCH_PADDING), &
     &                       (ZMIN_B + SEARCH_PADDING) - ZMAX_A))
          D_MIN_OUT = SQRT(DX*DX + DY*DY + DZ*DZ)

!         If the sampled point-cloud AABBs are already farther apart
!         than the padded search window, skip the voxel grid build.
          IF (D_MIN_OUT > SEARCH_PADDING) THEN
            RETURN
          END IF

!   ----- Step 4: voxel grid dimensions -----
          SPAN_X = XMAX_B - XMIN_B
          SPAN_Y = YMAX_B - YMIN_B
          SPAN_Z = ZMAX_B - ZMIN_B
          SPAN_X_SAFE = MAX(SPAN_X, EPS_SPAN)
          SPAN_Y_SAFE = MAX(SPAN_Y, EPS_SPAN)
          SPAN_Z_SAFE = MAX(SPAN_Z, EPS_SPAN)
          ! Compute the number of voxels in each direction
          NBX = MAX(1, CEILING(SPAN_X / CELL_SIZE))
          NBY = MAX(1, CEILING(SPAN_Y / CELL_SIZE))
          NBZ = MAX(1, CEILING(SPAN_Z / CELL_SIZE))
          NVOXELS = NBX * NBY * NBZ

!   ----- Step 5: allocate and fill voxel grid with B points -----
          ALLOCATE(VOXEL(NVOXELS))
          ALLOCATE(NEXT_PT(NPTS_B))
          VOXEL(1:NVOXELS) = 0
          NEXT_PT(1:NPTS_B) = 0

          ! Insert every B-point into its voxel cell
          DO IB = 1, NPTS_B
            RX = REAL(NBX, WP) * (SURF_POINTS_B(1,IB) - XMIN_B) / SPAN_X_SAFE
            RY = REAL(NBY, WP) * (SURF_POINTS_B(2,IB) - YMIN_B) / SPAN_Y_SAFE
            RZ = REAL(NBZ, WP) * (SURF_POINTS_B(3,IB) - ZMIN_B) / SPAN_Z_SAFE
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
            NEXT_PT(IB) = VOXEL(CELLID)
            VOXEL(CELLID) = IB
          END DO

!   ----- Step 6: query -- for each A-point, check 3x3x3 neighborhood -----
          DO IA = 1, NPTS_A
!           Skip A-points outside the padded B bounding box
            IF (SURF_POINTS_A(1,IA) < XMIN_B) CYCLE
            IF (SURF_POINTS_A(1,IA) > XMAX_B) CYCLE
            IF (SURF_POINTS_A(2,IA) < YMIN_B) CYCLE
            IF (SURF_POINTS_A(2,IA) > YMAX_B) CYCLE
            IF (SURF_POINTS_A(3,IA) < ZMIN_B) CYCLE
            IF (SURF_POINTS_A(3,IA) > ZMAX_B) CYCLE

            RX = REAL(NBX, WP) * (SURF_POINTS_A(1,IA) - XMIN_B) / SPAN_X_SAFE
            RY = REAL(NBY, WP) * (SURF_POINTS_A(2,IA) - YMIN_B) / SPAN_Y_SAFE
            RZ = REAL(NBZ, WP) * (SURF_POINTS_A(3,IA) - ZMIN_B) / SPAN_Z_SAFE
            IF (RX /= RX) RX = ZERO
            IF (RY /= RY) RY = ZERO
            IF (RZ /= RZ) RZ = ZERO
            IX = INT(RX)
            IY = INT(RY)
            IZ = INT(RZ)
            IX = MAX(0, MIN(NBX - 1, IX))
            IY = MAX(0, MIN(NBY - 1, IY))
            IZ = MAX(0, MIN(NBZ - 1, IZ))

            IX_LO = MAX(0, IX - 1)
            IX_HI = MIN(NBX - 1, IX + 1)
            IY_LO = MAX(0, IY - 1)
            IY_HI = MIN(NBY - 1, IY + 1)
            IZ_LO = MAX(0, IZ - 1)
            IZ_HI = MIN(NBZ - 1, IZ + 1)

            DO JZ = IZ_LO, IZ_HI
              DO JY = IY_LO, IY_HI
                DO JX = IX_LO, IX_HI
                  CELLID = JZ * NBX * NBY + JY * NBX + JX + 1
                  CELLID = MAX(1, MIN(NVOXELS, CELLID))
                  JJ = VOXEL(CELLID)
                  DO WHILE (JJ > 0)
                    IF (JJ > NPTS_B) EXIT
                    NCAND = CANDIDATE_COUNT(JJ)
                    IF (NCAND < SIZE(CANDIDATE_IA, 1)) THEN
                      CANDIDATE_COUNT(JJ) = NCAND + 1
                      CANDIDATE_IA(NCAND + 1, JJ) = IA
                    ELSE
                      CANDIDATE_OVERFLOW(JJ) = .TRUE.
                    END IF
                    DX = SURF_POINTS_A(1,IA) - SURF_POINTS_B(1,JJ)
                    DY = SURF_POINTS_A(2,IA) - SURF_POINTS_B(2,JJ)
                    DZ = SURF_POINTS_A(3,IA) - SURF_POINTS_B(3,JJ)
                    DIST_SQ = DX*DX + DY*DY + DZ*DZ
                    IF (DIST_SQ < D_MIN_SQ) THEN
                      D_MIN_SQ  = DIST_SQ 
                      ! IDX_A_OUT: index of the closest point in surface A
                      ! IDX_B_OUT: index of the closest point in surface B
                      IDX_A_OUT = IA
                      IDX_B_OUT = JJ
                    END IF
                    JJ = NEXT_PT(JJ)
                    IF (JJ < 0) EXIT
                  END DO
                END DO
              END DO
            END DO
          END DO

          DEALLOCATE(VOXEL, NEXT_PT)

        END SUBROUTINE Q1NP_CONTACT_BROAD_PHASE_VOXEL_MIN_DISTANCE

!=======================================================================
!   Q1NP_CONTACT_NARROW_PHASE_PROJECT
!
!   Narrow phase projection:
!     For each set-2 sample point, test all voxel-derived candidates
!     (or fall back to the global nearest search), keep the best valid
!     penetrating projection, and store one contact pair per B point.
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_NARROW_PHASE_PROJECT( &
     &      KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, X_COORDS, &
     &      NUMNOD, &
     &      SURF_POINTS_A, NPTS_A,                       &
     &      ELEM_IDS_A, XI_A, ETA_A,                     &
     &      SURF_POINTS_B, NPTS_B,                       &
     &      ELEM_IDS_B, XI_B, ETA_B,                     &
     &      CANDIDATE_IA, CANDIDATE_COUNT,               &
     &      CANDIDATE_OVERFLOW,                          &
     &      CONTACT_PAIRS, N_PAIRS, PENETRATION_DETECTED)
          INTEGER, INTENT(IN) :: KQ1NP_TAB(:,:)
          INTEGER, INTENT(IN) :: IQ1NP_TAB(:)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_KTAB(:)
          REAL(KIND=WP), INTENT(IN) :: X_COORDS(3,NUMNOD)
          INTEGER, INTENT(IN) :: NUMNOD
          REAL(KIND=WP), INTENT(IN) :: SURF_POINTS_A(3,NPTS_A)
          INTEGER, INTENT(IN) :: NPTS_A
          INTEGER, INTENT(IN) :: ELEM_IDS_A(:)
          REAL(KIND=WP), INTENT(IN) :: XI_A(:), ETA_A(:)
          REAL(KIND=WP), INTENT(IN) :: SURF_POINTS_B(3,NPTS_B)
          INTEGER, INTENT(IN) :: NPTS_B
          INTEGER, INTENT(IN) :: ELEM_IDS_B(:)
          REAL(KIND=WP), INTENT(IN) :: XI_B(:), ETA_B(:)
          INTEGER, INTENT(IN) :: CANDIDATE_IA(:,:)
          INTEGER, INTENT(IN) :: CANDIDATE_COUNT(:)
          LOGICAL, INTENT(IN) :: CANDIDATE_OVERFLOW(:)
          TYPE(Q1NP_CONTACT_PAIR), INTENT(OUT) :: CONTACT_PAIRS(:)
          INTEGER, INTENT(OUT) :: N_PAIRS
          LOGICAL, INTENT(OUT) :: PENETRATION_DETECTED

          INTEGER :: IA, IB, ICAND
          INTEGER :: BEST_IA
          REAL(KIND=WP) :: DX, DY, DZ, DIST_SQ, BEST_DSQ
          REAL(KIND=WP) :: X_SRC(3)

          REAL(KIND=WP) :: XI_PROJ, ETA_PROJ, XYZ_PROJ(3), PROJ_DIST
          REAL(KIND=WP) :: SIGNED_PENETRATION
          REAL(KIND=WP) :: NORMAL_VEC(3)
          REAL(KIND=WP) :: RESIDUAL
          LOGICAL :: PROJ_VALID
          INTEGER :: N_ITER
          REAL(KIND=WP) :: BEST_XI_PROJ, BEST_ETA_PROJ
          REAL(KIND=WP) :: BEST_SIGNED_PENETRATION, BEST_PROJ_DIST
          REAL(KIND=WP) :: BEST_NORMAL_VEC(3)
          LOGICAL :: BEST_PROJ_VALID

          N_PAIRS = 0
          PENETRATION_DETECTED = .FALSE.
          IF (NPTS_A < 1 .OR. NPTS_B < 1) RETURN

          DO IB = 1, NPTS_B
            X_SRC(1:3) = SURF_POINTS_B(1:3, IB)

            BEST_IA  = 0
            BEST_DSQ = HUGE(ONE)
            BEST_PROJ_VALID = .FALSE.
            BEST_SIGNED_PENETRATION = HUGE(ONE)
            BEST_PROJ_DIST = HUGE(ONE)
            IF (CANDIDATE_COUNT(IB) > 0 .AND. .NOT. CANDIDATE_OVERFLOW(IB)) THEN
              DO ICAND = 1, CANDIDATE_COUNT(IB)
                IA = CANDIDATE_IA(ICAND, IB)
                IF (IA < 1 .OR. IA > NPTS_A) CYCLE
                CALL Q1NP_CONTACT_PROJECT_POINT_NEWTON( &
     &            X_SRC, KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &            X_COORDS, NUMNOD, ELEM_IDS_A(IA),       &
     &            XI_A(IA), ETA_A(IA),                    &
     &            XI_PROJ, ETA_PROJ, XYZ_PROJ, PROJ_DIST, &
     &            SIGNED_PENETRATION, NORMAL_VEC,         &
     &            RESIDUAL, N_ITER, PROJ_VALID)
                IF (.NOT. PROJ_VALID) CYCLE
                IF (SIGNED_PENETRATION >= ZERO) CYCLE
                IF ((.NOT. BEST_PROJ_VALID) .OR. &
     &              (SIGNED_PENETRATION < BEST_SIGNED_PENETRATION) .OR. &
     &              (SIGNED_PENETRATION == BEST_SIGNED_PENETRATION .AND. &
     &               PROJ_DIST < BEST_PROJ_DIST)) THEN
                  BEST_PROJ_VALID = .TRUE.
                  BEST_IA = IA
                  BEST_XI_PROJ = XI_PROJ
                  BEST_ETA_PROJ = ETA_PROJ
                  BEST_SIGNED_PENETRATION = SIGNED_PENETRATION
                  BEST_PROJ_DIST = PROJ_DIST
                  BEST_NORMAL_VEC(1:3) = NORMAL_VEC(1:3)
                END IF
              END DO
            ELSE
              DO IA = 1, NPTS_A
                DX = X_SRC(1) - SURF_POINTS_A(1,IA)
                DY = X_SRC(2) - SURF_POINTS_A(2,IA)
                DZ = X_SRC(3) - SURF_POINTS_A(3,IA)
                DIST_SQ = DX*DX + DY*DY + DZ*DZ
                IF (DIST_SQ < BEST_DSQ) THEN
                  BEST_DSQ = DIST_SQ
                  BEST_IA  = IA
                END IF
              END DO

              IF (BEST_IA == 0) CYCLE
              
              ! Run Newton projection on the nearest point
              CALL Q1NP_CONTACT_PROJECT_POINT_NEWTON( &
     &            X_SRC, KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &            X_COORDS, NUMNOD, ELEM_IDS_A(BEST_IA),   &
     &            XI_A(BEST_IA), ETA_A(BEST_IA),            &
     &            XI_PROJ, ETA_PROJ, XYZ_PROJ, PROJ_DIST,   &
     &            SIGNED_PENETRATION, NORMAL_VEC,            &
     &            RESIDUAL, N_ITER, PROJ_VALID)
              IF (.NOT. PROJ_VALID) CYCLE
              IF (SIGNED_PENETRATION >= ZERO) CYCLE

              ! Store contact pair if projection is valid and signed penetration is negative
              BEST_PROJ_VALID = .TRUE.
              BEST_XI_PROJ = XI_PROJ
              BEST_ETA_PROJ = ETA_PROJ
              BEST_SIGNED_PENETRATION = SIGNED_PENETRATION
              BEST_PROJ_DIST = PROJ_DIST
              BEST_NORMAL_VEC(1:3) = NORMAL_VEC(1:3)
            END IF
            IF (.NOT. BEST_PROJ_VALID) CYCLE
            IF (N_PAIRS >= SIZE(CONTACT_PAIRS)) EXIT

            ! Store contact all detected contact pairs
            N_PAIRS = N_PAIRS + 1
            CONTACT_PAIRS(N_PAIRS)%PENETRATION = BEST_SIGNED_PENETRATION
            CONTACT_PAIRS(N_PAIRS)%NORMAL(1:3) = BEST_NORMAL_VEC(1:3)
            CONTACT_PAIRS(N_PAIRS)%XI_PROJ     = BEST_XI_PROJ
            CONTACT_PAIRS(N_PAIRS)%ETA_PROJ    = BEST_ETA_PROJ
            CONTACT_PAIRS(N_PAIRS)%ELEM_A      = ELEM_IDS_A(BEST_IA)
            CONTACT_PAIRS(N_PAIRS)%XI_SRC      = XI_B(IB)
            CONTACT_PAIRS(N_PAIRS)%ETA_SRC     = ETA_B(IB)
            CONTACT_PAIRS(N_PAIRS)%ELEM_B      = ELEM_IDS_B(IB)
          END DO

          PENETRATION_DETECTED = (N_PAIRS > 0)

        END SUBROUTINE Q1NP_CONTACT_NARROW_PHASE_PROJECT

!=======================================================================
!   Q1NP_CONTACT_COMPUTE_PENALTY_FORCES
!
!   For each penetrating contact pair, evaluate NURBS shape functions
!   on both surfaces and scatter penalty forces to control-point nodes.
!
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_COMPUTE_PENALTY_FORCES( &
     &      CONTACT_PAIRS, N_PAIRS, &
     &      KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &      X_COORDS, NUMNOD, A, STIFN)
          TYPE(Q1NP_CONTACT_PAIR), INTENT(IN) :: CONTACT_PAIRS(:)
          INTEGER, INTENT(IN) :: N_PAIRS
          INTEGER, INTENT(IN) :: KQ1NP_TAB(:,:)
          INTEGER, INTENT(IN) :: IQ1NP_TAB(:)
          REAL(KIND=WP), INTENT(IN)    :: Q1NP_KTAB(:)
          REAL(KIND=WP), INTENT(IN)    :: X_COORDS(3,NUMNOD)
          INTEGER, INTENT(IN)          :: NUMNOD
          REAL(KIND=WP), INTENT(INOUT) :: A(3,NUMNOD)
          REAL(KIND=WP), INTENT(INOUT) :: STIFN(NUMNOD)

          INTEGER, PARAMETER :: MAX_CTRL = 50
          REAL(KIND=WP), PARAMETER :: EPS = 1.0E-10_WP
          REAL(KIND=WP), PARAMETER :: PREC = 1.0E-6_WP
          INTEGER :: IP, K, GID, ELEM_IDX
          INTEGER :: P_CUR, Q_CUR, NCTRL
          INTEGER :: ELEM_U_IDX, ELEM_V_IDX
          INTEGER :: CTRL_IDS_A(MAX_CTRL), CTRL_IDS_B(MAX_CTRL)
          REAL(KIND=WP), ALLOCATABLE :: U_KNOT_LOCAL(:), V_KNOT_LOCAL(:)
          REAL(KIND=WP), ALLOCATABLE :: STIFN_BASE(:)
          REAL(KIND=WP) :: NVAL_A(MAX_CTRL), NVAL_B(MAX_CTRL)
          REAL(KIND=WP) :: F_PEN(3), F_MAG
          REAL(KIND=WP) :: F_MAX, PENETR_MAX
          REAL(KIND=WP) :: PEN_ABS, PEN_EFF, FAC, D1
          REAL(KIND=WP) :: K_A, K_B, K_PAIR, GAP_REF
          REAL(KIND=WP) :: PAIR_WEIGHT, D1_WEIGHTED

          IF (N_PAIRS < 1) RETURN

          F_MAX     = ZERO
          PENETR_MAX = ZERO

          GAP_REF = Q1NP_CONTACT_GAP_REFERENCE
          PAIR_WEIGHT = ONE / REAL(Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_U * &
     &                             Q1NP_CONTACT_BROAD_PHASE_NGP_SURF_V, WP)
          ALLOCATE(STIFN_BASE(NUMNOD))
          STIFN_BASE(1:NUMNOD) = STIFN(1:NUMNOD)

          DO IP = 1, N_PAIRS
            PEN_ABS = ABS(CONTACT_PAIRS(IP)%PENETRATION)
            IF (PEN_ABS <= ZERO) CYCLE
            IF (Q1NP_CONTACT_PENETRATION_RATIO_MAX > ZERO) THEN
              PEN_EFF = MIN(PEN_ABS, &
     &                      Q1NP_CONTACT_PENETRATION_RATIO_MAX * GAP_REF)
            ELSE
              PEN_EFF = PEN_ABS
            END IF
            IF (PEN_EFF <= ZERO) CYCLE

!           --- Target surface A: evaluate shape values and side stiffness ---
            ELEM_IDX = CONTACT_PAIRS(IP)%ELEM_A
            CALL Q1NP_CONTACT_EXTRACT_ELEM_DATA( &
     &          ELEM_IDX, KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &          P_CUR, Q_CUR, NCTRL, ELEM_U_IDX, ELEM_V_IDX, &
     &          CTRL_IDS_A, U_KNOT_LOCAL, V_KNOT_LOCAL)

            CALL Q1NP_EVALUATE_NURBS_SHAPE_VALUES( &
     &          CONTACT_PAIRS(IP)%XI_PROJ, &
     &          CONTACT_PAIRS(IP)%ETA_PROJ, &
     &          P_CUR, Q_CUR, ELEM_U_IDX, ELEM_V_IDX, &
     &          U_KNOT_LOCAL, V_KNOT_LOCAL, &
     &          MIN(NCTRL, MAX_CTRL), NVAL_A)

            K_A = ZERO
            DO K = 1, MIN(NCTRL, MAX_CTRL)
              GID = CTRL_IDS_A(K)
              IF (GID <= 0 .OR. GID > NUMNOD) CYCLE
              K_A = K_A + ABS(NVAL_A(K)) * &
     &              MIN(ABS(STIFN_BASE(GID)), Q1NP_CONTACT_PENALTY_KMAX)
            END DO
            DEALLOCATE(U_KNOT_LOCAL, V_KNOT_LOCAL)

!           --- Source surface B: evaluate shape values and side stiffness ---
            ELEM_IDX = CONTACT_PAIRS(IP)%ELEM_B
            CALL Q1NP_CONTACT_EXTRACT_ELEM_DATA( &
     &          ELEM_IDX, KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &          P_CUR, Q_CUR, NCTRL, ELEM_U_IDX, ELEM_V_IDX, &
     &          CTRL_IDS_B, U_KNOT_LOCAL, V_KNOT_LOCAL)

            CALL Q1NP_EVALUATE_NURBS_SHAPE_VALUES( &
     &          CONTACT_PAIRS(IP)%XI_SRC, &
     &          CONTACT_PAIRS(IP)%ETA_SRC, &
     &          P_CUR, Q_CUR, ELEM_U_IDX, ELEM_V_IDX, &
     &          U_KNOT_LOCAL, V_KNOT_LOCAL, &
     &          MIN(NCTRL, MAX_CTRL), NVAL_B)

            K_B = ZERO
            DO K = 1, MIN(NCTRL, MAX_CTRL)
              GID = CTRL_IDS_B(K)
              IF (GID <= 0 .OR. GID > NUMNOD) CYCLE
              K_B = K_B + ABS(NVAL_B(K)) * &
     &              MIN(ABS(STIFN_BASE(GID)), Q1NP_CONTACT_PENALTY_KMAX)
            END DO
            DEALLOCATE(U_KNOT_LOCAL, V_KNOT_LOCAL)

!           --- STS-style pair stiffness blending + clamping ---
            IF (K_A > ZERO .AND. K_B > ZERO) THEN
              K_PAIR = (K_A * K_B) / MAX(EPS, (K_A + K_B))
            ELSE
              K_PAIR = MAX(K_A, K_B)
            END IF
            K_PAIR = Q1NP_CONTACT_PENALTY_SCALE * K_PAIR
            K_PAIR = MAX(Q1NP_CONTACT_PENALTY_KMIN, &
     &                   MIN(K_PAIR, Q1NP_CONTACT_PENALTY_KMAX))

!           --- STS-like gap scaling ---
            IF (GAP_REF > EPS .AND. PEN_EFF < GAP_REF) THEN
              FAC = GAP_REF / MAX(EPS, (GAP_REF - PEN_EFF))
              FAC = MIN(FAC, Q1NP_CONTACT_FAC_MAX)
            ELSE
              FAC = Q1NP_CONTACT_FAC_MAX
            END IF
            D1 = 0.5 * K_PAIR * FAC
            IF (D1 <= ZERO) CYCLE
            D1_WEIGHTED = PAIR_WEIGHT * D1

            F_MAG = D1_WEIGHTED * PEN_EFF
            IF (Q1NP_CONTACT_FORCE_CAP > ZERO) THEN
              F_MAG = MIN(F_MAG, Q1NP_CONTACT_FORCE_CAP)
            END IF
            F_PEN(1:3) = F_MAG * CONTACT_PAIRS(IP)%NORMAL(1:3)

            IF (F_MAG > F_MAX) F_MAX = F_MAG
            IF (PEN_ABS > PENETR_MAX) PENETR_MAX = PEN_ABS

!           --- Scatter to A: -N_k * F ---
            DO K = 1, MIN(KQ1NP_TAB(3, CONTACT_PAIRS(IP)%ELEM_A), MAX_CTRL)
              GID = CTRL_IDS_A(K)
              IF (GID <= 0 .OR. GID > NUMNOD) CYCLE
              A(1, GID) = A(1, GID) - NVAL_A(K) * F_PEN(1)
              A(2, GID) = A(2, GID) - NVAL_A(K) * F_PEN(2)
              A(3, GID) = A(3, GID) - NVAL_A(K) * F_PEN(3)
              STIFN(GID) = STIFN(GID) + D1_WEIGHTED * ABS(NVAL_A(K))
            END DO

!           --- Scatter to B: +N_j * F ---
            DO K = 1, MIN(KQ1NP_TAB(3, CONTACT_PAIRS(IP)%ELEM_B), MAX_CTRL)
              GID = CTRL_IDS_B(K)
              IF (GID <= 0 .OR. GID > NUMNOD) CYCLE
              A(1, GID) = A(1, GID) + NVAL_B(K) * F_PEN(1)
              A(2, GID) = A(2, GID) + NVAL_B(K) * F_PEN(2)
              A(3, GID) = A(3, GID) + NVAL_B(K) * F_PEN(3)
              STIFN(GID) = STIFN(GID) + D1_WEIGHTED * ABS(NVAL_B(K))
            END DO
          END DO

          DEALLOCATE(STIFN_BASE)

        END SUBROUTINE Q1NP_CONTACT_COMPUTE_PENALTY_FORCES

!=======================================================================
!   Q1NP_CONTACT_EXTRACT_ELEM_DATA
!
!   Helper: extract element metadata, control-point IDs, and knot
!   vectors for a given Q1NP element index.  The caller must
!   DEALLOCATE U_KNOT_LOCAL and V_KNOT_LOCAL after use.
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_EXTRACT_ELEM_DATA( &
     &      ELEM_IDX, KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &      P_OUT, Q_OUT, NCTRL_OUT, &
     &      ELEM_U_IDX_OUT, ELEM_V_IDX_OUT, &
     &      CTRL_IDS_OUT, U_KNOT_OUT, V_KNOT_OUT)
          INTEGER, INTENT(IN) :: ELEM_IDX
          INTEGER, INTENT(IN) :: KQ1NP_TAB(:,:)
          INTEGER, INTENT(IN) :: IQ1NP_TAB(:)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_KTAB(:)
          INTEGER, INTENT(OUT) :: P_OUT, Q_OUT, NCTRL_OUT
          INTEGER, INTENT(OUT) :: ELEM_U_IDX_OUT, ELEM_V_IDX_OUT
          INTEGER, PARAMETER :: MAX_CTRL = 50
          INTEGER, INTENT(OUT) :: CTRL_IDS_OUT(MAX_CTRL)
          REAL(KIND=WP), ALLOCATABLE, INTENT(OUT) :: U_KNOT_OUT(:)
          REAL(KIND=WP), ALLOCATABLE, INTENT(OUT) :: V_KNOT_OUT(:)

          INTEGER :: CP_OFFSET, KNOT_SET_ID
          INTEGER :: NX_CUR, NY_CUR, U_LEN, V_LEN, IPT

          P_OUT          = KQ1NP_TAB(8, ELEM_IDX)
          Q_OUT          = KQ1NP_TAB(9, ELEM_IDX)
          NCTRL_OUT      = KQ1NP_TAB(3, ELEM_IDX)
          CP_OFFSET      = KQ1NP_TAB(4, ELEM_IDX)
          ELEM_U_IDX_OUT = KQ1NP_TAB(6, ELEM_IDX)
          ELEM_V_IDX_OUT = KQ1NP_TAB(7, ELEM_IDX)
          KNOT_SET_ID    = KQ1NP_TAB(15, ELEM_IDX)

          NX_CUR = KQ1NP_TAB(12, ELEM_IDX)
          NY_CUR = KQ1NP_TAB(13, ELEM_IDX)
          IF (NX_CUR <= 0 .OR. NY_CUR <= 0) THEN
            IF (Q1NP_NKNOT_SETS_G > 0 .AND. &
     &          KNOT_SET_ID > 0 .AND. &
     &          KNOT_SET_ID <= Q1NP_NKNOT_SETS_G) THEN
              NX_CUR = Q1NP_NX_SET_G(KNOT_SET_ID)
              NY_CUR = Q1NP_NY_SET_G(KNOT_SET_ID)
            ELSE
              NX_CUR = Q1NP_NX_G
              NY_CUR = Q1NP_NY_G
            END IF
          END IF

          U_LEN = NX_CUR + 2*P_OUT + 1
          V_LEN = NY_CUR + 2*Q_OUT + 1
          ALLOCATE(U_KNOT_OUT(U_LEN), V_KNOT_OUT(V_LEN))

          IF (Q1NP_NKNOT_SETS_G > 0 .AND. &
     &        KNOT_SET_ID > 0 .AND. &
     &        KNOT_SET_ID <= Q1NP_NKNOT_SETS_G .AND. &
     &        ALLOCATED(Q1NP_KTAB_OFF_G)) THEN
            U_KNOT_OUT(1:U_LEN) = &
     &          Q1NP_KTAB(Q1NP_KTAB_OFF_G(KNOT_SET_ID) : &
     &                    Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN - 1)
            V_KNOT_OUT(1:V_LEN) = &
     &          Q1NP_KTAB(Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN : &
     &                    Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN + V_LEN - 1)
          ELSE
            U_KNOT_OUT(1:U_LEN) = Q1NP_KTAB(1:U_LEN)
            V_KNOT_OUT(1:V_LEN) = Q1NP_KTAB(U_LEN+1:U_LEN+V_LEN)
          END IF

          DO IPT = 1, MIN(NCTRL_OUT, MAX_CTRL)
            CTRL_IDS_OUT(IPT) = IQ1NP_TAB(CP_OFFSET + IPT - 1)
          END DO

        END SUBROUTINE Q1NP_CONTACT_EXTRACT_ELEM_DATA

!=======================================================================
!   Q1NP_CONTACT_PROJECT_POINT_NEWTON
!
!   Project a source point X_SRC onto the NURBS top surface of a
!   specific Q1NP element using Newton iteration on the orthogonality
!   conditions:  f1 = (S - X_SRC) . dS/dXI  = 0
!                f2 = (S - X_SRC) . dS/dETA = 0
!
!   Returns parametric coordinates, projected point, distance,
!   signed penetration, unit normal, residual, and validity flag.
!=======================================================================
        SUBROUTINE Q1NP_CONTACT_PROJECT_POINT_NEWTON( &
     &      X_SRC, KQ1NP_TAB, IQ1NP_TAB, Q1NP_KTAB, &
     &      X_COORDS, NUMNOD, ELEM_IDX, &
     &      XI_SEED, ETA_SEED, &
     &      XI_OUT, ETA_OUT, XYZ_PROJ, PROJ_DIST, &
     &      SIGNED_PENETRATION_OUT, NORMAL_OUT, &
     &      RESIDUAL_OUT, N_ITER_OUT, VALID)
          REAL(KIND=WP), INTENT(IN)  :: X_SRC(3)
          INTEGER, INTENT(IN) :: KQ1NP_TAB(:,:)
          INTEGER, INTENT(IN) :: IQ1NP_TAB(:)
          REAL(KIND=WP), INTENT(IN) :: Q1NP_KTAB(:)
          REAL(KIND=WP), INTENT(IN) :: X_COORDS(3,NUMNOD)
          INTEGER, INTENT(IN) :: NUMNOD, ELEM_IDX
          REAL(KIND=WP), INTENT(IN)  :: XI_SEED, ETA_SEED
          REAL(KIND=WP), INTENT(OUT) :: XI_OUT, ETA_OUT
          REAL(KIND=WP), INTENT(OUT) :: XYZ_PROJ(3), PROJ_DIST
          REAL(KIND=WP), INTENT(OUT) :: SIGNED_PENETRATION_OUT
          REAL(KIND=WP), INTENT(OUT) :: NORMAL_OUT(3)
          REAL(KIND=WP), INTENT(OUT) :: RESIDUAL_OUT
          INTEGER, INTENT(OUT) :: N_ITER_OUT
          LOGICAL, INTENT(OUT) :: VALID

          INTEGER, PARAMETER :: MAX_ITER = 5
          REAL(KIND=WP), PARAMETER :: NEWTON_TOL = 1.0E-10_WP
          INTEGER, PARAMETER :: MAX_CTRL = 50

          INTEGER :: ITER, KNOT_SET_ID
          INTEGER :: P_CUR, Q_CUR, NCTRL, CP_OFFSET
          INTEGER :: ELEM_U_IDX, ELEM_V_IDX
          INTEGER :: NX_CUR, NY_CUR, U_LEN, V_LEN, IPT
          INTEGER :: CTRL_IDS(MAX_CTRL)
          REAL(KIND=WP), ALLOCATABLE :: U_KNOT_LOCAL(:), V_KNOT_LOCAL(:)
          REAL(KIND=WP) :: XI, ETA
          REAL(KIND=WP) :: S(3), SU(3), SV(3), DIFF(3)
          REAL(KIND=WP) :: NORMAL(3), NORM_N
          REAL(KIND=WP) :: F1, F2, RES_NORM
          REAL(KIND=WP) :: A11, A12, A21, A22, DET_J
          REAL(KIND=WP) :: D_XI, D_ETA

          VALID = .FALSE.

!         Extract element metadata
          P_CUR      = KQ1NP_TAB(8, ELEM_IDX)
          Q_CUR      = KQ1NP_TAB(9, ELEM_IDX)
          NCTRL      = KQ1NP_TAB(3, ELEM_IDX)
          CP_OFFSET  = KQ1NP_TAB(4, ELEM_IDX)
          ELEM_U_IDX = KQ1NP_TAB(6, ELEM_IDX)
          ELEM_V_IDX = KQ1NP_TAB(7, ELEM_IDX)
          KNOT_SET_ID = KQ1NP_TAB(15, ELEM_IDX)

          NX_CUR = KQ1NP_TAB(12, ELEM_IDX)
          NY_CUR = KQ1NP_TAB(13, ELEM_IDX)
          IF (NX_CUR <= 0 .OR. NY_CUR <= 0) THEN
            IF (Q1NP_NKNOT_SETS_G > 0 .AND. &
     &          KNOT_SET_ID > 0 .AND. &
     &          KNOT_SET_ID <= Q1NP_NKNOT_SETS_G) THEN
              NX_CUR = Q1NP_NX_SET_G(KNOT_SET_ID)
              NY_CUR = Q1NP_NY_SET_G(KNOT_SET_ID)
            ELSE
              NX_CUR = Q1NP_NX_G
              NY_CUR = Q1NP_NY_G
            END IF
          END IF

          U_LEN = NX_CUR + 2*P_CUR + 1
          V_LEN = NY_CUR + 2*Q_CUR + 1
          ALLOCATE(U_KNOT_LOCAL(U_LEN), V_KNOT_LOCAL(V_LEN))

          IF (Q1NP_NKNOT_SETS_G > 0 .AND. &
     &        KNOT_SET_ID > 0 .AND. &
     &        KNOT_SET_ID <= Q1NP_NKNOT_SETS_G .AND. &
     &        ALLOCATED(Q1NP_KTAB_OFF_G)) THEN
            U_KNOT_LOCAL(1:U_LEN) = &
     &          Q1NP_KTAB(Q1NP_KTAB_OFF_G(KNOT_SET_ID) : &
     &                    Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN - 1)
            V_KNOT_LOCAL(1:V_LEN) = &
     &          Q1NP_KTAB(Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN : &
     &                    Q1NP_KTAB_OFF_G(KNOT_SET_ID) + U_LEN + V_LEN - 1)
          ELSE
            U_KNOT_LOCAL(1:U_LEN) = Q1NP_KTAB(1:U_LEN)
            V_KNOT_LOCAL(1:V_LEN) = Q1NP_KTAB(U_LEN+1:U_LEN+V_LEN)
          END IF

          DO IPT = 1, MIN(NCTRL, MAX_CTRL)
            CTRL_IDS(IPT) = IQ1NP_TAB(CP_OFFSET + IPT - 1)
          END DO

          XI  = XI_SEED
          ETA = ETA_SEED

          DO ITER = 1, MAX_ITER
            CALL Q1NP_EVALUATE_NURBS_TOP_SURFACE_POINT_AND_DERIVS( &
     &          XI, ETA, P_CUR, Q_CUR, &
     &          ELEM_U_IDX, ELEM_V_IDX, &
     &          U_KNOT_LOCAL, V_KNOT_LOCAL, &
     &          MIN(NCTRL, MAX_CTRL), CTRL_IDS, &
     &          X_COORDS, NUMNOD, S, SU, SV)

            DIFF(1:3) = S(1:3) - X_SRC(1:3)

!           Orthogonality residuals
            F1 = DIFF(1)*SU(1) + DIFF(2)*SU(2) + DIFF(3)*SU(3)
            F2 = DIFF(1)*SV(1) + DIFF(2)*SV(2) + DIFF(3)*SV(3)

            RES_NORM = SQRT(F1*F1 + F2*F2)
            IF (RES_NORM < NEWTON_TOL) EXIT

!           Jacobian of the residual system (first-order approximation):
!             A11 = Su . Su,  A12 = Su . Sv
!             A21 = Sv . Su,  A22 = Sv . Sv
            A11 = SU(1)*SU(1) + SU(2)*SU(2) + SU(3)*SU(3)
            A12 = SU(1)*SV(1) + SU(2)*SV(2) + SU(3)*SV(3)
            A21 = A12
            A22 = SV(1)*SV(1) + SV(2)*SV(2) + SV(3)*SV(3)

            DET_J = A11*A22 - A12*A21
            IF (ABS(DET_J) < 1.0E-30_WP) EXIT

!           Solve 2x2 system: [A] * [d_xi, d_eta]^T = -[f1, f2]^T
            D_XI  = -(A22*F1 - A12*F2) / DET_J
            D_ETA = -(A11*F2 - A21*F1) / DET_J

            XI  = XI  + D_XI
            ETA = ETA + D_ETA

!           Clamp to parent domain [-1, +1]
            XI  = MAX(-ONE, MIN(ONE, XI))
            ETA = MAX(-ONE, MIN(ONE, ETA))
          END DO

!         Final evaluation at converged point
          CALL Q1NP_EVALUATE_NURBS_TOP_SURFACE_POINT_AND_DERIVS( &
     &        XI, ETA, P_CUR, Q_CUR, &
     &        ELEM_U_IDX, ELEM_V_IDX, &
     &        U_KNOT_LOCAL, V_KNOT_LOCAL, &
     &        MIN(NCTRL, MAX_CTRL), CTRL_IDS, &
     &        X_COORDS, NUMNOD, S, SU, SV)

          DIFF(1:3) = S(1:3) - X_SRC(1:3)
          F1 = DIFF(1)*SU(1) + DIFF(2)*SU(2) + DIFF(3)*SU(3)
          F2 = DIFF(1)*SV(1) + DIFF(2)*SV(2) + DIFF(3)*SV(3)
          RES_NORM = SQRT(F1*F1 + F2*F2)

          XI_OUT   = XI
          ETA_OUT  = ETA
          XYZ_PROJ = S
          PROJ_DIST = SQRT(DIFF(1)**2 + DIFF(2)**2 + DIFF(3)**2)
          NORMAL(1) = SU(2)*SV(3) - SU(3)*SV(2)
          NORMAL(2) = SU(3)*SV(1) - SU(1)*SV(3)
          NORMAL(3) = SU(1)*SV(2) - SU(2)*SV(1)
          NORM_N = SQRT(NORMAL(1)**2 + NORMAL(2)**2 + NORMAL(3)**2)
          IF (NORM_N > 1.0E-30_WP) THEN
            NORMAL(1:3) = NORMAL(1:3) / NORM_N
            SIGNED_PENETRATION_OUT = (X_SRC(1)-S(1))*NORMAL(1) + &
     &                               (X_SRC(2)-S(2))*NORMAL(2) + &
     &                               (X_SRC(3)-S(3))*NORMAL(3)
          ELSE
            SIGNED_PENETRATION_OUT = ZERO
          END IF
          NORMAL_OUT(1:3) = NORMAL(1:3)
          RESIDUAL_OUT = RES_NORM
          N_ITER_OUT   = ITER

!         Converged = residual small AND parametric coords inside domain
          VALID = (RES_NORM < NEWTON_TOL * 1.0E3_WP) .AND. &
     &            (ABS(XI) <= ONE) .AND. (ABS(ETA) <= ONE)

          DEALLOCATE(U_KNOT_LOCAL, V_KNOT_LOCAL)

        END SUBROUTINE Q1NP_CONTACT_PROJECT_POINT_NEWTON

      END MODULE Q1NP_CONTACT_ALGORITHMS_MOD
