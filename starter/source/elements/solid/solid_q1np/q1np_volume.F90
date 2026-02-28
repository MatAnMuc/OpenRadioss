!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        The Free Software Foundation, either version 3 of the License, or
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
!Chd|====================================================================
!Chd|  Q1NP_VOLUME                   source/elements/solid/solid_q1np/q1np_volume.F
!Chd|====================================================================
!C=======================================================================
!C   Volume computation for Q1NP enriched elements
!C
!C   This module implements numerical volume integration for Q1NP elements
!C   using Gauss quadrature over enriched shape functions (NURBS top +
      module q1np_volume_mod
        use message_mod
        use q1np_restart_mod
        use precision_mod, only : WP
        use constant_mod, only : ZERO, ONE, HALF, FOURTH
        implicit none
        integer, parameter :: IDEBUG_Q1NP_VOL = 0
      contains
!C
!C   Routines (volume):
!C   1. Q1NP_GET_KNOT_VECTORS       - Extract U,V knot vectors from Q1NP_KTAB
!C   2. Q1NP_DERS_BASIS_FUNS       - B-spline basis and derivatives (Cox-de Boor)
!C   3. Q1NP_SHAPE_FUNCTIONS       - Q1NP shape functions (NURBS top + bilinear bottom)
!C   4. Q1NP_JACOBIAN              - Jacobian matrix and determinant
!C   5. Q1NP_COMPUTE_VOLUME_ELEMENT - Volume integration via Gauss quadrature
!C   6. Q1NP_FIND_SPAN             - Knot span index for parameter value
!C   7. Q1NP_DEBUG_BASIS_AT_UV     - Debug print of basis at (u,v)
!C   8. Q1NP_BERNSTEIN_BASIS       - Bernstein basis (single-span case)
!C   9. Q1NP_KNOT_SINGLE_SPAN      - True if knot vector has no interior knots
!C  10. Q1NP_BASIS_ROW_AT_UV       - One row of design matrix at (u,v)
!C=======================================================================

!C=======================================================================
!C   Main volume integration over element routine
!C=======================================================================
        subroutine q1np_compute_volume_element(iel_q1np, &
     &                                         kq1np_tab, iq1np_tab, &
     &                                         iq1np_bulk_tab, q1np_ktab, &
     &                                         x, nx, ny, vol_el)
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: IEL_Q1NP, NX, NY
INTEGER, INTENT(IN) :: KQ1NP_TAB(15,*), IQ1NP_TAB(*)
INTEGER, INTENT(IN) :: IQ1NP_BULK_TAB(*)
real(kind=WP), INTENT(IN) :: Q1NP_KTAB(*), X(3,*)
real(kind=WP), INTENT(OUT) :: VOL_EL
! IEL_Q1NP: element index
! NX, NY: number of control points in U, V direction
! KQ1NP_TAB: element metadata
! IQ1NP_TAB: control point indices
! IQ1NP_BULK_TAB: bulk node indices
! Q1NP_KTAB: knot vector
! X: coordinates of nodes
! VOL_EL: volume of element
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: P, Q, NCTRL, OFFSET_CTRL, OFFSET_BULK
INTEGER :: ELEM_U, ELEM_V
INTEGER :: NKNOT_U, NKNOT_V
INTEGER :: I, J, K, IU, IV
INTEGER :: N_TOP, N_TOTAL
INTEGER, ALLOCATABLE :: NODE_IDS(:)
real(kind=WP) :: XI, ETA, ZETA, WG
real(kind=WP), ALLOCATABLE :: NVAL(:)
real(kind=WP), ALLOCATABLE :: DN_LOCAL(:,:)
real(kind=WP), ALLOCATABLE :: XNODE(:,:)
real(kind=WP), ALLOCATABLE :: U(:), V(:)
real(kind=WP) :: JMAT(3,3), DETJ
real(kind=WP) :: DETJ_MIN, DETJ_MAX
real(kind=WP) :: SUM_W
real(kind=WP) :: SUM_N, SUM_N_MIN, SUM_N_MAX
real(kind=WP) :: DETJ2D, DET2D_MIN, DET2D_MAX, DETJ2D_GP
real(kind=WP) :: AREA2D, SUM_W2D
INTEGER :: IC, IB
real(kind=WP) :: XI_DBG(4), ETA_DBG(4), ZETA_DBG
!C----------------------------------------------------------------------
!C   Get order (P,Q), control-point count, offsets, and element spans
!C----------------------------------------------------------------------
P           = KQ1NP_TAB(8, IEL_Q1NP)
Q           = KQ1NP_TAB(9, IEL_Q1NP)
NCTRL       = KQ1NP_TAB(3, IEL_Q1NP)
OFFSET_CTRL = KQ1NP_TAB(4, IEL_Q1NP)
OFFSET_BULK = KQ1NP_TAB(14,IEL_Q1NP)
ELEM_U      = KQ1NP_TAB(6, IEL_Q1NP)
ELEM_V      = KQ1NP_TAB(7, IEL_Q1NP)

!C----------------------------------------------------------------------
!C   Allocate shape-function and node arrays (N_TOP control pts + 4 bulk)
!C----------------------------------------------------------------------
N_TOP   = NCTRL
N_TOTAL = N_TOP + 4
ALLOCATE(NODE_IDS(N_TOTAL))
ALLOCATE(NVAL(N_TOTAL))
ALLOCATE(DN_LOCAL(N_TOTAL,3))
ALLOCATE(XNODE(3,N_TOTAL))

!C----------------------------------------------------------------------
!C   Allocate and extract U,V knot vectors from global Q1NP_KTAB
!C----------------------------------------------------------------------
NKNOT_U = NX + 2*P + 1
NKNOT_V = NY + 2*Q + 1
ALLOCATE(U(NKNOT_U))
ALLOCATE(V(NKNOT_V))

CALL Q1NP_GET_KNOT_VECTORS(NX, NY, P, Q, Q1NP_KTAB, U, V)

!C----------------------------------------------------------------------
!C   Build node list: first NCTRL control points, then 4 bulk nodes
!C----------------------------------------------------------------------
!C     Control points (top surface)
DO I = 1, NCTRL
  NODE_IDS(I) = IQ1NP_TAB(OFFSET_CTRL + I - 1)
end do
      
!C     Bulk nodes (bottom face)
DO I = 1, 4
  NODE_IDS(N_TOP + I) = IQ1NP_BULK_TAB(OFFSET_BULK + I - 1)
end do

!C     Extract coordinates
DO K = 1, N_TOTAL
  XNODE(1,K) = X(1, NODE_IDS(K))
  XNODE(2,K) = X(2, NODE_IDS(K))
  XNODE(3,K) = X(3, NODE_IDS(K))
end do

!C----------------------------------------------------------------------
!C   Optional: 2D surface area of TOP NURBS surface (ZETA=+1) for debug
!C   A2D = ∫_{-1}^{1}∫_{-1}^{1} det(J_2D) dxi deta
!C----------------------------------------------------------------------
IF (IDEBUG_Q1NP_VOL > 0) THEN
  AREA2D    = ZERO
  SUM_W2D   = ZERO
  DET2D_MIN = HUGE(ZERO)
  DET2D_MAX = -HUGE(ZERO)

  DO IU = 1, Q1NP_NP_U_G
    XI = Q1NP_GP_U_G(IU)
    DO IV = 1, Q1NP_NP_V_G
      ETA  = Q1NP_GP_V_G(IV)
      ZETA = ONE

      WG = Q1NP_GW_U_G(IU) * Q1NP_GW_V_G(IV)

      CALL Q1NP_SHAPE_FUNCTIONS(XI, ETA, ZETA, P, Q, U, V, &
     &                          ELEM_U, ELEM_V, NVAL, DN_LOCAL)
      CALL Q1NP_JACOBIAN(DN_LOCAL, XNODE, N_TOTAL, JMAT, DETJ)

      DETJ2D = JMAT(1,1) * JMAT(2,2) - JMAT(1,2) * JMAT(2,1)

      AREA2D    = AREA2D    + WG * DETJ2D
      SUM_W2D   = SUM_W2D   + WG
      DET2D_MIN = MIN(DET2D_MIN, DETJ2D)
      DET2D_MAX = MAX(DET2D_MAX, DETJ2D)
    END DO
  END DO

  IF (IEL_Q1NP == 1) THEN
    WRITE(*,'(A,I6,4(A,1P,E12.5))') &
   &  'Q1NP 2D AREA DBG: IEL_Q1NP=', IEL_Q1NP, &
   &  ' A2D=', AREA2D, &
   &  ' SUM_W2D=', SUM_W2D, &
   &  ' DET2D_MIN=', DET2D_MIN, &
   &  ' DET2D_MAX=', DET2D_MAX
  END IF
END IF

!C----------------------------------------------------------------------
!C   Volume via 2D surface integral at ZETA=+1 (NURBS top surface).
!C
!C   V = 2 * ∫∫ det(J(ξ,η,ζ=+1)) dξ dη
!C
!C   At ζ=+1 the in-plane Jacobian rows come purely from the NURBS
!C   surface (NZ_bot=0), while ∂r/∂ζ = ½(r_top−r_bot) gives the
!C   element height vector.  The scalar triple product
!C   (r_top−r_bot)·(t₁×t₂) = 2·det(J)|_{ζ=+1} is the exact volume
!C   element of the ruled body between the NURBS top and bilinear
!C   bottom surfaces.  This avoids the parametric mixing that would
!C   otherwise reduce det(J) at intermediate ζ values.
!C----------------------------------------------------------------------
VOL_EL     = ZERO
DETJ_MIN   = HUGE(DETJ)
DETJ_MAX   = -HUGE(DETJ)
SUM_W      = ZERO
SUM_N_MIN  = HUGE(SUM_N)
SUM_N_MAX  = -HUGE(SUM_N)

ZETA = ONE

DO IU = 1, Q1NP_NP_U_G
  XI = Q1NP_GP_U_G(IU)
  DO IV = 1, Q1NP_NP_V_G
    ETA = Q1NP_GP_V_G(IV)
    WG  = Q1NP_GW_U_G(IU) * Q1NP_GW_V_G(IV)

    CALL Q1NP_SHAPE_FUNCTIONS(XI, ETA, ZETA, P, Q, U, V, &
     &                        ELEM_U, ELEM_V, NVAL, DN_LOCAL)

    SUM_N = ZERO
    DO K = 1, N_TOTAL
      SUM_N = SUM_N + NVAL(K)
    END DO
    SUM_N_MIN = MIN(SUM_N_MIN, SUM_N)
    SUM_N_MAX = MAX(SUM_N_MAX, SUM_N)

    CALL Q1NP_JACOBIAN(DN_LOCAL, XNODE, N_TOTAL, JMAT, DETJ)

    VOL_EL   = VOL_EL + WG * DETJ
    SUM_W    = SUM_W  + WG
    DETJ_MIN = MIN(DETJ_MIN, DETJ)
    DETJ_MAX = MAX(DETJ_MAX, DETJ)
  END DO
END DO

VOL_EL = 2 * VOL_EL

!C----------------------------------------------------------------------
!C   Optional debug output: volume diagnostics
!C----------------------------------------------------------------------
IF (IDEBUG_Q1NP_VOL >= 2) THEN
  WRITE(*,'(A,I6,6(A,1P,E12.5))') &
     & 'Q1NP VOL DBG: IEL_Q1NP=', IEL_Q1NP, &
     & ' VOL_EL=', VOL_EL, &
     & ' SUM_W=', SUM_W, &
     & ' DETJ_MIN=', DETJ_MIN, &
     & ' DETJ_MAX=', DETJ_MAX, &
     & ' SUM_N_MIN=', SUM_N_MIN, &
     & ' SUM_N_MAX=', SUM_N_MAX
END IF

!C----------------------------------------------------------------------
!C   Deallocate temporary arrays
!C----------------------------------------------------------------------
DEALLOCATE(U, V)
DEALLOCATE(NODE_IDS, NVAL, DN_LOCAL, XNODE)

          return ! Volume of element computed
        end subroutine q1np_compute_volume_element




!C=======================================================================
!C   Extract U,V knot vectors from Q1NP_KTAB
!C=======================================================================
        subroutine q1np_get_knot_vectors(nx, ny, p, q, q1np_ktab, u, v)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
          integer, intent(in) :: nx, ny, p, q
          real(kind=WP), intent(in) :: q1np_ktab(*)
          real(kind=WP), intent(out) :: u(*), v(*)
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
          integer :: nknot_u, nknot_v, iv_offset, i
!C----------------------------------------------------------------------
!C   Knot vector length: n_knots = n_control + 2*order + 1 (open B-spline)
!C----------------------------------------------------------------------
          nknot_u = nx + 2*p + 1
          nknot_v = ny + 2*q + 1
          iv_offset = nknot_u   ! V starts after U in Q1NP_KTAB

!C=======================================================================
!C   2. Extract U knot vector (first NKNOT_U entries)
!C=======================================================================
          do i = 1, nknot_u
            u(i) = q1np_ktab(i)
          end do

!C=======================================================================
!C   3. Extract V knot vector (starts at NKNOT_U + 1)
!C=======================================================================
          do i = 1, nknot_v
            v(i) = q1np_ktab(iv_offset + i)
          end do

          return ! U and V knot vectors extracted
        end subroutine q1np_get_knot_vectors

!C=======================================================================
!C   B-spline basis functions and derivatives computation
!C=======================================================================
SUBROUTINE Q1NP_DERS_BASIS_FUNS(SPAN, UVAL, P, U, NDERS, DERS)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: SPAN, P, NDERS
real(kind=WP), INTENT(IN) :: UVAL
real(kind=WP), INTENT(IN) :: U(:)
real(kind=WP), INTENT(OUT) :: DERS(0:NDERS,0:P)
! SPAN: element span index
! P: polynomial order
! NDERS: number of derivatives to compute
! UVAL: value of the knot vector at the current span
! U: knot vector
! DERS: output array for basis functions and derivatives
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: J, R, K, RK, PK, J1, J2
real(kind=WP) :: LEFT(0:P), RIGHT(0:P)
real(kind=WP) :: NDU(0:P,0:P)
real(kind=WP) :: A(0:1,0:P)
real(kind=WP) :: TEMP, DENOM, D, RFACT
INTEGER :: S1, S2
real(kind=WP), PARAMETER :: TOL_NDU = 1.0E-15

!C----------------------------------------------------------------------
!C   Build NDU: Cox-de Boor recursion. NDU(j,r) = knot diff; NDU(r,j) = basis.
!C----------------------------------------------------------------------
NDU(0,0) = ONE

DO J = 1, P
  LEFT(J)  = UVAL - U(SPAN + 1 - J)
  RIGHT(J) = U(SPAN + J) - UVAL
  TEMP = ZERO

  DO R = 0, J-1
    DENOM = RIGHT(R+1) + LEFT(J-R)
    IF (ABS(DENOM) .LE. TOL_NDU) THEN
      NDU(J,R) = ZERO
    ELSE
      NDU(J,R) = DENOM
    end if

    IF (ABS(NDU(J,R)) .GT. TOL_NDU) THEN
      DENOM = NDU(R,J-1) / NDU(J,R)
    ELSE
      DENOM = ZERO
    end if

    NDU(R,J) = TEMP + RIGHT(R+1) * DENOM
    TEMP     = LEFT(J-R) * DENOM
  end do

  NDU(J,J) = TEMP
end do

!C=======================================================================
!C   2. Extract basis functions DERS(0,J) = NDU(J,P)
!C=======================================================================
DO J = 0, P
  DERS(0,J) = NDU(J,P)
end do

IF (NDERS <= 0) THEN
  RETURN
end if

!C=======================================================================
!C   3. Compute derivatives
!C=======================================================================
DO R = 0, P
  S1 = 0
  S2 = 1
  A(0,0) = ONE
  DO J = 1, P
    A(0,J) = ZERO
    A(1,J) = ZERO
  end do

  DO K = 1, NDERS
    D  = ZERO
    RK = R - K
    PK = P - K

    IF (R >= K) THEN
      DENOM = NDU(PK+1,RK)
      IF (ABS(DENOM) .GT. TOL_NDU) THEN
        A(S2,0) = A(S1,0) / DENOM
      ELSE
        A(S2,0) = ZERO
      end if
      D = A(S2,0) * NDU(RK,PK)
    end if

    IF (RK >= -1) THEN
      J1 = 1
    ELSE
      J1 = -RK
    end if

    IF (R-1 <= PK) THEN
      J2 = K-1
    ELSE
      J2 = P-R
    end if

    DO J = J1, J2
      DENOM = NDU(PK+1,RK+J)
      IF (ABS(DENOM) .GT. TOL_NDU) THEN
        A(S2,J) = (A(S1,J) - A(S1,J-1)) / DENOM
      ELSE
        A(S2,J) = ZERO
      end if
      D = D + A(S2,J) * NDU(RK+J,PK)
    end do

    IF (R <= PK) THEN
      DENOM = NDU(PK+1,R)
      IF (ABS(DENOM) .GT. TOL_NDU) THEN
        A(S2,K) = -A(S1,K-1) / DENOM
      ELSE
        A(S2,K) = ZERO
      end if
      D = D + A(S2,K) * NDU(R,PK)
    end if

    DERS(K,R) = D

    S1 = S2
    S2 = 1 - S1
  end do
end do

!C=======================================================================
!C   4. Multiply by factorial factors
!C=======================================================================
RFACT = REAL(P)
DO K = 1, NDERS
  DO J = 0, P
    DERS(K,J) = DERS(K,J) * RFACT
  end do
  RFACT = RFACT * REAL(P-K)
end do

          return ! Basis functions and derivatives computed
      ! DERS(0,I): basis function
      ! DERS(1,I): derivative of basis function with respect to U
        end subroutine q1np_ders_basis_funs



!C=======================================================================
!C   Q1NP shape functions and derivatives computation (NURBS top + bilinear bottom)
!C=======================================================================
        subroutine q1np_shape_functions(xi, eta, zeta, p, q, u, v, &
     &                                  elem_u, elem_v, n, dn_local)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: P, Q, ELEM_U, ELEM_V
real(kind=WP), INTENT(IN) :: XI, ETA, ZETA
real(kind=WP), INTENT(IN) :: U(:), V(:)
real(kind=WP), INTENT(OUT) :: N(:), DN_LOCAL(:,:)
! XI, ETA, ZETA: local coordinates in parent element
! P, Q: polynomial order
! U, V: knot vectors
! ELEM_U, ELEM_V: element span indices in U, V direction
! N: output array for shape functions
! DN_LOCAL: output array for derivatives of shape functions
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: N_TOP, N_TOTAL, SU, SV
INTEGER :: I, J, IDX
real(kind=WP) :: XI_LOC, ETA_LOC
real(kind=WP) :: AU, BU, AV, BV
real(kind=WP) :: UVAL, VVAL
real(kind=WP) :: KSPAN_U, KSPAN_V
real(kind=WP) :: DU_DXI, DV_DETA
real(kind=WP) :: NU_DERS(0:1,0:P), NV_DERS(0:1,0:Q)
real(kind=WP) :: NU(0:P), DNU_DU(0:P)
real(kind=WP) :: NV(0:Q), DNV_DV(0:Q)
real(kind=WP) :: N_TOP_ARY((P+1)*(Q+1))
real(kind=WP) :: DN_TOP_DU((P+1)*(Q+1))
real(kind=WP) :: DN_TOP_DV((P+1)*(Q+1))
real(kind=WP) :: N_BOT(4), DN_BOT_DU(4), DN_BOT_DV(4)
real(kind=WP) :: NZ_TOP, NZ_BOT, DNZ_TOP, DNZ_BOT
real(kind=WP) :: INV_KSPAN_UV
!C=======================================================================
!C   1. Initialize N_TOP and N_TOTAL (number of shape functions on top and bottom surfaces)
!C=======================================================================
N_TOP   = (P+1)*(Q+1)   ! 9 CP on top surface
N_TOTAL = N_TOP + 4     ! 9 CP + 4 BULK nodes

!C=======================================================================
!C   2. Map from parent [-1,1] to local [0,1]
!C=======================================================================
XI_LOC  = HALF * (XI   + ONE)
ETA_LOC = HALF * (ETA  + ONE)

!C----------------------------------------------------------------------
!C   Element span in U,V and current knot span [AU,BU] x [AV,BV]
!C   ELEM_U,ELEM_V are 0-based element indices; convert to span index:
!C     For open uniform B-splines, spans run from P+1 to NE+P.
!C     With 1-based element index I=ELEM_U+1, span = P + I.
!C----------------------------------------------------------------------
SU = P + ELEM_U + 1
SV = Q + ELEM_V + 1
AU = U(SU)
BU = U(SU+1)
AV = V(SV)
BV = V(SV+1)

!C----------------------------------------------------------------------
!C   Map parent (XI,ETA) in [-1,1] to (UVAL,VVAL) in current knot span
!C----------------------------------------------------------------------
UVAL = AU + (BU - AU) * XI_LOC
VVAL = AV + (BV - AV) * ETA_LOC

!C=======================================================================
!C   5. Compute scaling factors for derivatives
!C   Mapping: XI in [-1,1] -> XI_LOC in [0,1] -> [AU,BU]
!C   dU/dXI = (BU-AU) * dXI_LOC/dXI = (BU-AU)/2
!C=======================================================================
KSPAN_U = BU - AU
KSPAN_V = BV - AV
DU_DXI  = HALF * KSPAN_U
DV_DETA = HALF * KSPAN_V

!C=======================================================================
!C   6. Evaluate B-spline basis functions and derivatives
!C=======================================================================
CALL Q1NP_DERS_BASIS_FUNS(SU, UVAL, P, U, 1, NU_DERS)
CALL Q1NP_DERS_BASIS_FUNS(SV, VVAL, Q, V, 1, NV_DERS)

DO I = 0, P
  NU(I)     = NU_DERS(0,I)
  DNU_DU(I) = NU_DERS(1,I)
end do
DO J = 0, Q
  NV(J)     = NV_DERS(0,J)
  DNV_DV(J) = NV_DERS(1,J)
end do

!C=======================================================================
!C   7. Tensor product for top surface
!C=======================================================================
IDX = 1
DO J = 0, Q
  DO I = 0, P
    N_TOP_ARY(IDX)   = NU(I)     * NV(J)
    DN_TOP_DU(IDX)   = DNU_DU(I) * NV(J)
    DN_TOP_DV(IDX)   = NU(I)     * DNV_DV(J)
    IDX = IDX + 1
  end do
end do

!C=======================================================================
!C   8. Linear blending in through-thickness direction
!C=======================================================================
NZ_TOP  = HALF * (ONE + ZETA)
NZ_BOT  = HALF * (ONE - ZETA)
DNZ_TOP = HALF
DNZ_BOT = -HALF

!C=======================================================================
!C   9. Bottom corner shape functions in (u,v) on patch [AU,BU] x [AV,BV]
!C     Same parametrization as top so in-plane scale (DU_DXI,DV_DETA) matches.
!C     Corner order: 1=(AU,AV), 2=(BU,AV), 3=(BU,BV), 4=(AU,BV).
!C=======================================================================
INV_KSPAN_UV = ONE / (KSPAN_U * KSPAN_V)
N_BOT(1) = (BU - UVAL) * (BV - VVAL) * INV_KSPAN_UV
N_BOT(2) = (UVAL - AU) * (BV - VVAL) * INV_KSPAN_UV
N_BOT(3) = (UVAL - AU) * (VVAL - AV) * INV_KSPAN_UV
N_BOT(4) = (BU - UVAL) * (VVAL - AV) * INV_KSPAN_UV

DN_BOT_DU(1) = -(BV - VVAL) * INV_KSPAN_UV
DN_BOT_DU(2) =  (BV - VVAL) * INV_KSPAN_UV
DN_BOT_DU(3) =  (VVAL - AV) * INV_KSPAN_UV
DN_BOT_DU(4) = -(VVAL - AV) * INV_KSPAN_UV

DN_BOT_DV(1) = -(BU - UVAL) * INV_KSPAN_UV
DN_BOT_DV(2) = -(UVAL - AU) * INV_KSPAN_UV
DN_BOT_DV(3) =  (UVAL - AU) * INV_KSPAN_UV
DN_BOT_DV(4) =  (BU - UVAL) * INV_KSPAN_UV

!C=======================================================================
!C   10. Combine top and bottom
!C=======================================================================
!C   Top surface nodes (NURBS). DN_LOCAL(I,J) = dN_I / d xi_J
DO I = 1, N_TOP
  N(I) = N_TOP_ARY(I) * NZ_TOP
  DN_LOCAL(I,1) = DN_TOP_DU(I) * DU_DXI  * NZ_TOP
  DN_LOCAL(I,2) = DN_TOP_DV(I) * DV_DETA * NZ_TOP
  DN_LOCAL(I,3) = N_TOP_ARY(I) * DNZ_TOP
end do

!C     Bottom: same (u,v) chain as top -> use DN_BOT_DU*DU_DXI, DN_BOT_DV*DV_DETA
DO I = 1, 4
  N(N_TOP + I) = N_BOT(I) * NZ_BOT
  DN_LOCAL(N_TOP + I,1) = DN_BOT_DU(I) * DU_DXI  * NZ_BOT
  DN_LOCAL(N_TOP + I,2) = DN_BOT_DV(I) * DV_DETA * NZ_BOT
  DN_LOCAL(N_TOP + I,3) = N_BOT(I)       * DNZ_BOT
end do

          return ! Q1NP shape functions N and derivatives DN_LOCAL computation
        end subroutine q1np_shape_functions

!C=======================================================================
!C   Jacobian matrix and determinant computation
!C=======================================================================
SUBROUTINE Q1NP_JACOBIAN(DN_LOCAL, XNODE, NNODE, J, DETJ)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: NNODE
real(kind=WP), INTENT(IN) :: DN_LOCAL(NNODE,3), XNODE(3,NNODE)
real(kind=WP), INTENT(OUT) :: J(3,3), DETJ
! NNODE: number of nodes
! DN_LOCAL: derivatives of shape functions
! XNODE: coordinates of nodes
! J: Jacobian matrix
! DETJ: determinant of Jacobian matrix
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: K
!C=======================================================================
!C   1. Assemble Jacobian matrix: J(i,j) = sum_k XNODE(i,k) * DN_LOCAL(k,j)
!C=======================================================================
J = ZERO
DO K = 1, NNODE
  J(1,1) = J(1,1) + XNODE(1,K) * DN_LOCAL(K,1)
  J(1,2) = J(1,2) + XNODE(1,K) * DN_LOCAL(K,2)
  J(1,3) = J(1,3) + XNODE(1,K) * DN_LOCAL(K,3)

  J(2,1) = J(2,1) + XNODE(2,K) * DN_LOCAL(K,1)
  J(2,2) = J(2,2) + XNODE(2,K) * DN_LOCAL(K,2)
  J(2,3) = J(2,3) + XNODE(2,K) * DN_LOCAL(K,3)

  J(3,1) = J(3,1) + XNODE(3,K) * DN_LOCAL(K,1)
  J(3,2) = J(3,2) + XNODE(3,K) * DN_LOCAL(K,2)
  J(3,3) = J(3,3) + XNODE(3,K) * DN_LOCAL(K,3)

end do

!C=======================================================================
!C   2. Compute determinant (3x3 explicit formula)
!C=======================================================================
          detj = j(1,1)*(j(2,2)*j(3,3) - j(2,3)*j(3,2)) - &
     &           j(1,2)*(j(2,1)*j(3,3) - j(2,3)*j(3,1)) + &
     &           j(1,3)*(j(2,1)*j(3,2) - j(2,2)*j(3,1))

          return ! Jacobian matrix and determinant computed
        end subroutine q1np_jacobian


!C=======================================================================
!C   Q1NP_FIND_SPAN: find knot span index i such that U(i) <= UVAL < U(i+1).
!C   Used by least-squares fitting and basis evaluation. DERS_BASIS_FUNS
!C   requires SPAN with SPAN+1+P <= NK, hence upper bound SPAN = NK-P-1.
!C=======================================================================
SUBROUTINE Q1NP_FIND_SPAN(U, NK, P, UVAL, SPAN)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: NK, P
INTEGER, INTENT(OUT) :: SPAN
real(kind=WP), INTENT(IN) :: U(*), UVAL
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: I
!C----------------------------------------------------------------------
IF (UVAL >= U(NK)) THEN
  SPAN = NK - P - 1
  RETURN
end if
SPAN = 1
DO I = 2, NK - 1
  IF (UVAL >= U(I) .AND. UVAL < U(I+1)) THEN
    SPAN = I
    RETURN
  end if
end do
          return
        end subroutine q1np_find_span

!C=======================================================================
!C   Q1NP_DEBUG_BASIS_AT_UV: debug output for basis at (UU,VV).
!C   Prints spans and sum of U/V basis (should be 1). Clamps (UU,VV) off
!C   0/1 to avoid degenerate Cox-de Boor at boundaries.
!C=======================================================================
SUBROUTINE Q1NP_DEBUG_BASIS_AT_UV(UU, VV, U, V, P, Q, NCP_U, NCP_V)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: P, Q, NCP_U, NCP_V
real(kind=WP), INTENT(IN) :: UU, VV
real(kind=WP), INTENT(IN) :: U(:), V(:)
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: SPAN_U, SPAN_V, NKU, NKV
real(kind=WP) :: NU_DERS(0:0,0:10), NV_DERS(0:0,0:10)
real(kind=WP) :: UU_EVAL, VV_EVAL
real(kind=WP), PARAMETER :: EPS_BOUND = 1.0E-5
!C----------------------------------------------------------------------
NKU = SIZE(U)
NKV = SIZE(V)
UU_EVAL = UU
VV_EVAL = VV
IF (UU_EVAL .LE. ZERO)  UU_EVAL = EPS_BOUND
IF (UU_EVAL .GE. ONE)   UU_EVAL = ONE - EPS_BOUND
IF (VV_EVAL .LE. ZERO)  VV_EVAL = EPS_BOUND
IF (VV_EVAL .GE. ONE)   VV_EVAL = ONE - EPS_BOUND
CALL Q1NP_FIND_SPAN(U, NCP_U + P + 1, P, UU_EVAL, SPAN_U)
CALL Q1NP_FIND_SPAN(V, NCP_V + Q + 1, Q, VV_EVAL, SPAN_V)
IF (Q1NP_KNOT_SINGLE_SPAN(U, NKU)) THEN
  CALL Q1NP_BERNSTEIN_BASIS(P, UU_EVAL, NU_DERS(0,0:P))
ELSE
  CALL Q1NP_DERS_BASIS_FUNS(SPAN_U, UU_EVAL, P, U, 0, NU_DERS)
end if
IF (Q1NP_KNOT_SINGLE_SPAN(V, NKV)) THEN
  CALL Q1NP_BERNSTEIN_BASIS(Q, VV_EVAL, NV_DERS(0,0:Q))
ELSE
  CALL Q1NP_DERS_BASIS_FUNS(SPAN_V, VV_EVAL, Q, V, 0, NV_DERS)
end if
          print *, 'Q1NP DEBUG: u,v=', uu, vv, ' span=', span_u, span_v, &
     &        ' sum(NU)=', sum(nu_ders(0,0:p)), ' sum(NV)=', sum(nv_ders(0,0:q))
          return
        end subroutine q1np_debug_basis_at_uv

!C=======================================================================
!C   Q1NP_BERNSTEIN_BASIS: Bernstein basis of degree P at T in [0,1].
!C   B(i) = C(P,i)*T^i*(1-T)^(P-i). Used when knot vector has no interior
!C   knots (single span [0,1]); Cox-de Boor would give 0/0 there.
!C=======================================================================
SUBROUTINE Q1NP_BERNSTEIN_BASIS(P, T, B)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: P
real(kind=WP), INTENT(IN) :: T
real(kind=WP), INTENT(OUT) :: B(0:P)
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: I, K
real(kind=WP) :: COEF, TI, T1I
!C----------------------------------------------------------------------
IF (P .LE. 0) THEN
  B(0) = ONE
  RETURN
end if
DO I = 0, P
  COEF = ONE
  DO K = 1, I
    COEF = COEF * REAL(P - K + 1) / REAL(K)
  end do
  TI  = MERGE(ONE, T**I, I .EQ. 0)
  T1I = MERGE(ONE, (ONE - T)**(P - I), I .EQ. P)
  b(i) = COEF * TI * T1I
end do
          return
        end subroutine q1np_bernstein_basis

!C=======================================================================
!C   Q1NP_KNOT_SINGLE_SPAN: .TRUE. if knot vector has no interior knots
!C   (only 0 and 1). Then Bernstein basis is used instead of Cox-de Boor.
!C=======================================================================
LOGICAL FUNCTION Q1NP_KNOT_SINGLE_SPAN(V, NK)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: NK
real(kind=WP), INTENT(IN) :: V(NK)
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: I
real(kind=WP), PARAMETER :: TOL = 1.0E-9
!C----------------------------------------------------------------------
Q1NP_KNOT_SINGLE_SPAN = .TRUE.
DO I = 1, NK
  IF (V(I) .GT. TOL .AND. V(I) .LT. (ONE - TOL)) THEN
    Q1NP_KNOT_SINGLE_SPAN = .FALSE.
    RETURN
  end if
end do
          return
        end function q1np_knot_single_span

!C=======================================================================
!C   Q1NP_BASIS_ROW_AT_UV: fill one row of design matrix for least-squares.
!C   A_ROW(1:NCP_U*NCP_V) = tensor-product B-spline basis at (UU,VV).
!C   Column index: (JJ-1)*NCP_U + II for control point (II,JJ).
!C=======================================================================
        subroutine q1np_basis_row_at_uv(uu, vv, u, v, p, q, ncp_u, ncp_v, &
     &                                  a_row)
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
INTEGER, INTENT(IN) :: P, Q, NCP_U, NCP_V
real(kind=WP), INTENT(IN) :: UU, VV
real(kind=WP), INTENT(IN) :: U(:), V(:)
real(kind=WP), INTENT(OUT) :: A_ROW(*)
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
INTEGER :: SPAN_U, SPAN_V, COL, II, JJ, I, J, NKU, NKV
real(kind=WP) :: NU_DERS(0:0,0:10), NV_DERS(0:0,0:10)
real(kind=WP) :: UU_EVAL, VV_EVAL
LOGICAL :: U_SINGLE_SPAN, V_SINGLE_SPAN
real(kind=WP), PARAMETER :: EPS_PARAM = 1.0E-5
!C----------------------------------------------------------------------
!C   Clamp (UU,VV) off 0/1 to avoid 0/0 in Cox-de Boor at boundaries
!C----------------------------------------------------------------------
UU_EVAL = UU
VV_EVAL = VV
IF (UU_EVAL .LE. ZERO) UU_EVAL = EPS_PARAM
IF (UU_EVAL .GE. ONE)  UU_EVAL = ONE - EPS_PARAM
IF (VV_EVAL .LE. ZERO) VV_EVAL = EPS_PARAM
IF (VV_EVAL .GE. ONE)  VV_EVAL = ONE - EPS_PARAM

DO COL = 1, NCP_U * NCP_V
  A_ROW(COL) = ZERO
end do
CALL Q1NP_FIND_SPAN(U, NCP_U + P + 1, P, UU_EVAL, SPAN_U)
CALL Q1NP_FIND_SPAN(V, NCP_V + Q + 1, Q, VV_EVAL, SPAN_V)
NKU = SIZE(U)
NKV = SIZE(V)
U_SINGLE_SPAN = Q1NP_KNOT_SINGLE_SPAN(U, NKU)
V_SINGLE_SPAN = Q1NP_KNOT_SINGLE_SPAN(V, NKV)

IF (U_SINGLE_SPAN) THEN
  CALL Q1NP_BERNSTEIN_BASIS(P, UU_EVAL, NU_DERS(0,0:P))
ELSE
  CALL Q1NP_DERS_BASIS_FUNS(SPAN_U, UU_EVAL, P, U, 0, NU_DERS)
end if
IF (V_SINGLE_SPAN) THEN
  CALL Q1NP_BERNSTEIN_BASIS(Q, VV_EVAL, NV_DERS(0,0:Q))
ELSE
  CALL Q1NP_DERS_BASIS_FUNS(SPAN_V, VV_EVAL, Q, V, 0, NV_DERS)
end if

DO J = 0, Q
  JJ = SPAN_V - Q + J
  IF (JJ < 1 .OR. JJ > NCP_V) CYCLE
  DO I = 0, P
    II = SPAN_U - P + I
    IF (II < 1 .OR. II > NCP_U) CYCLE
    COL = (JJ - 1) * NCP_U + II
    A_ROW(COL) = NU_DERS(0,I) * NV_DERS(0,J)
  end do
end do
          return
        end subroutine q1np_basis_row_at_uv

!C=======================================================================
      end module q1np_volume_mod
