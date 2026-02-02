Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        The Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
Chd|====================================================================
Chd|  Q1NP_VOLUME                   source/elements/solid/solid_q1np/q1np_volume.F
Chd|====================================================================
C=======================================================================
C   Volume computation for Q1NP enriched elements
C
C   This module implements numerical volume integration for Q1NP elements
C   using Gauss quadrature over enriched shape functions (NURBS top +
C   bilinear bottom nodes). All routines are in MODULE Q1NP_VOLUME_MOD.
C
C   Routines:
C   1. Q1NP_GET_KNOT_VECTORS - Extract U,V knot vectors from Q1NP_KTAB
C   2. Q1NP_DERS_BASIS_FUNS  - B-spline basis functions and derivatives computation
C   3. Q1NP_SHAPE_FUNCTIONS  - Q1NP shape functions (NURBS + bilinear)
C   4. Q1NP_JACOBIAN         - Jacobian matrix and determinant
C   5. Q1NP_COMPUTE_VOLUME_ELEMENT - Main volume integration routine
C=======================================================================
      MODULE Q1NP_VOLUME_MOD
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
C=======================================================================
      CONTAINS


C=======================================================================
C   Extract U,V knot vectors from Q1NP_KTAB
C=======================================================================
      SUBROUTINE Q1NP_GET_KNOT_VECTORS(NX, NY, P, Q, Q1NP_KTAB, U, V)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NX, NY, P, Q
      my_real, INTENT(IN) :: Q1NP_KTAB(*)
      my_real, INTENT(OUT) :: U(*), V(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NKNOT_U, NKNOT_V, IV_OFFSET, I
C=======================================================================
C   1. Calculate knot vector lengths
C=======================================================================
      NKNOT_U = NX + 2*P + 1
      NKNOT_V = NY + 2*Q + 1
      IV_OFFSET = NKNOT_U

C=======================================================================
C   2. Extract U knot vector (first NKNOT_U entries)
C=======================================================================
      DO I = 1, NKNOT_U
        U(I) = Q1NP_KTAB(I)
      ENDDO

C=======================================================================
C   3. Extract V knot vector (starts at NKNOT_U + 1)
C=======================================================================
      DO I = 1, NKNOT_V
        V(I) = Q1NP_KTAB(IV_OFFSET + I)
      ENDDO

      RETURN ! U and V knot vectors extracted
      END

C=======================================================================
C   B-spline basis functions and derivatives computation
C=======================================================================
      SUBROUTINE Q1NP_DERS_BASIS_FUNS(SPAN, UVAL, P, U, NDERS, DERS)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: SPAN, P, NDERS
      my_real, INTENT(IN) :: UVAL
      my_real, INTENT(IN) :: U(:)
      my_real, INTENT(OUT) :: DERS(0:NDERS,0:P)
      ! SPAN: element span index
      ! P: polynomial order
      ! NDERS: number of derivatives to compute
      ! UVAL: value of the knot vector at the current span
      ! U: knot vector
      ! DERS: output array for basis functions and derivatives
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, R, K, RK, PK, J1, J2
      my_real :: LEFT(0:P), RIGHT(0:P)
      my_real :: NDU(0:P,0:P)
      my_real :: A(0:1,0:P)
      my_real :: TEMP, DEN, D, RFACT
      INTEGER :: S1, S2

C=======================================================================
C   1. Build NDU table (upper triangle: basis functions, lower: knot diffs)
C=======================================================================
      NDU(0,0) = ONE

      DO J = 1, P
        LEFT(J)  = UVAL - U(SPAN + 2 - J)
        RIGHT(J) = U(SPAN + 1 + J) - UVAL
        TEMP = ZERO

        DO R = 0, J - 1
          NDU(J,R) = RIGHT(R+1) + LEFT(J-R)
          IF (ABS(NDU(J,R)) .LT. 1.0E-15) THEN
            DEN = ZERO
          ELSE
            DEN = NDU(R,J-1) / NDU(J,R)
          END IF
          NDU(R,J) = TEMP + RIGHT(R+1) * DEN
          TEMP     = LEFT(J-R) * DEN
        END DO
        NDU(J,J) = TEMP
      END DO

C=======================================================================
C   2. Extract basis functions DERS(0,J) = NDU(J,P)
C=======================================================================
      DO J = 0, P
        DERS(0,J) = NDU(J,P)
      ENDDO

C=======================================================================
C   Compute derivatives if requested
C=======================================================================
      IF (NDERS > 0) THEN
        DO R = 0, P
          S1 = 0
          S2 = 1
          A(0,0) = ONE
          A(1,0:P) = ZERO

          DO K = 1, NDERS
            D  = ZERO
            RK = R - K
            PK = P - K

            IF (R >= K) THEN
              IF (ABS(NDU(PK+1,RK)) < 1.0E-15) THEN
                A(S2,0) = ZERO
              ELSE
                A(S2,0) = A(S1,0) / NDU(PK+1,RK)
              ENDIF
              D = A(S2,0) * NDU(RK,PK)
            ENDIF

            IF (RK >= -1) THEN
              J1 = 1
            ELSE
              J1 = -RK
            ENDIF
            IF (R - 1 <= PK) THEN
              J2 = K - 1
            ELSE
              J2 = P - R
            ENDIF

            DO J = J1, J2
              DEN = NDU(PK+1,RK+J)
              IF (ABS(DEN) < 1.0E-15) THEN
                A(S2,J) = ZERO
              ELSE
                A(S2,J) = (A(S1,J) - A(S1,J-1)) / DEN
              ENDIF
              D = D + A(S2,J) * NDU(RK+J,PK)
            ENDDO

            IF (R <= PK) THEN
              DEN = NDU(PK+1,R)
              IF (ABS(DEN) < 1.0E-15) THEN
                A(S2,K) = ZERO
              ELSE
                A(S2,K) = -A(S1,K-1) / DEN
              ENDIF
              D = D + A(S2,K) * NDU(R,PK)
            ENDIF

            DERS(K,R) = D

            ! Toggle rows
            S1 = S2
            S2 = 1 - S1
          ENDDO
        ENDDO

C=======================================================================
C   3. Multiply by factorial factors
C=======================================================================
        RFACT = P
        DO K = 1, NDERS
          DO J = 0, P
            DERS(K,J) = DERS(K,J) * REAL(RFACT)
          ENDDO
          RFACT = RFACT * (P - K)
        ENDDO
      ENDIF

      RETURN ! Basis functions and derivatives computed
      ! DERS(0,I): basis function
      ! DERS(1,I): derivative of basis function with respect to U
      END



C=======================================================================
C   Q1NP shape functions and derivatives computation (NURBS top + bilinear bottom)
C=======================================================================
      SUBROUTINE Q1NP_SHAPE_FUNCTIONS(XI, ETA, ZETA, P, Q, U, V,
     .                                ELEM_U, ELEM_V, N, DN_LOCAL)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: P, Q, ELEM_U, ELEM_V
      my_real, INTENT(IN) :: XI, ETA, ZETA
      my_real, INTENT(IN) :: U(:), V(:)
      my_real, INTENT(OUT) :: N(:), DN_LOCAL(:,:)
      ! XI, ETA, ZETA: local coordinates in parent element
      ! P, Q: polynomial order
      ! U, V: knot vectors
      ! ELEM_U, ELEM_V: element span indices in U, V direction
      ! N: output array for shape functions
      ! DN_LOCAL: output array for derivatives of shape functions
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: N_TOP, N_TOTAL, SU, SV
      INTEGER :: I, J, IDX
      my_real :: XI_LOC, ETA_LOC
      my_real :: AU, BU, AV, BV
      my_real :: UVAL, VVAL
      my_real :: KSPAN_U, KSPAN_V
      my_real :: DU_DXI, DV_DETA
      my_real :: NU_DERS(0:1,0:P), NV_DERS(0:1,0:Q)
      my_real :: NU(0:P), DNU_DU(0:P)
      my_real :: NV(0:Q), DNV_DV(0:Q)
      my_real :: N_TOP_ARY((P+1)*(Q+1))
      my_real :: DN_TOP_DU((P+1)*(Q+1))
      my_real :: DN_TOP_DV((P+1)*(Q+1))
      my_real :: N_BOT(4), DN_BOT_DXI(4), DN_BOT_DETA(4)
      my_real :: NZ_TOP, NZ_BOT, DNZ_TOP, DNZ_BOT
C=======================================================================
C   1. Initialize N_TOP and N_TOTAL (number of shape functions on top and bottom surfaces)
C=======================================================================
      N_TOP   = (P+1)*(Q+1)   ! 9 CP on top surface
      N_TOTAL = N_TOP + 4     ! 9 CP + 4 BULK nodes

C=======================================================================
C   2. Map from parent [-1,1] to local [0,1]
C=======================================================================
      XI_LOC  = HALF * (XI   + ONE)
      ETA_LOC = HALF * (ETA  + ONE)

C=======================================================================
C   3. Determine element spans and knot span boundaries
C=======================================================================
      SU = P + ELEM_U
      SV = Q + ELEM_V

      AU = U(SU+1)
      BU = U(SU+1+1)
      AV = V(SV+1)
      BV = V(SV+1+1)

C=======================================================================
C   4. Map to actual knot span
C=======================================================================
      UVAL = AU + (BU - AU) * XI_LOC
      VVAL = AV + (BV - AV) * ETA_LOC

C=======================================================================
C   5. Compute scaling factors for derivatives
C=======================================================================
      KSPAN_U = BU - AU
      KSPAN_V = BV - AV
      DU_DXI  = HALF * KSPAN_U
      DV_DETA = HALF * KSPAN_V

C=======================================================================
C   6. Evaluate B-spline basis functions and derivatives
C=======================================================================
      CALL Q1NP_DERS_BASIS_FUNS(SU, UVAL, P, U, 1, NU_DERS)
      CALL Q1NP_DERS_BASIS_FUNS(SV, VVAL, Q, V, 1, NV_DERS)

      DO I = 0, P
        NU(I)     = NU_DERS(0,I)
        DNU_DU(I) = NU_DERS(1,I)
      ENDDO
      DO J = 0, Q
        NV(J)     = NV_DERS(0,J)
        DNV_DV(J) = NV_DERS(1,J)
      ENDDO

C=======================================================================
C   7. Tensor product for top surface
C=======================================================================
      IDX = 1
      DO J = 0, Q
        DO I = 0, P
          N_TOP_ARY(IDX)   = NU(I)     * NV(J)
          DN_TOP_DU(IDX)   = DNU_DU(I) * NV(J)
          DN_TOP_DV(IDX)   = NU(I)     * DNV_DV(J)
          IDX = IDX + 1
        ENDDO
      ENDDO

C=======================================================================
C   8. Linear blending in through-thickness direction
C=======================================================================
      NZ_TOP  = HALF * (ONE + ZETA)
      NZ_BOT  = HALF * (ONE - ZETA)
      DNZ_TOP = HALF
      DNZ_BOT = -HALF

C=======================================================================
C   9. Bottom corner shape functions (bilinear Q1)
C=======================================================================
      N_BOT(1) = FOURTH * (ONE - XI) * (ONE - ETA)
      N_BOT(2) = FOURTH * (ONE + XI) * (ONE - ETA)
      N_BOT(3) = FOURTH * (ONE + XI) * (ONE + ETA)
      N_BOT(4) = FOURTH * (ONE - XI) * (ONE + ETA)

      DN_BOT_DXI(1)  = -FOURTH * (ONE - ETA)
      DN_BOT_DXI(2)  =  FOURTH * (ONE - ETA)
      DN_BOT_DXI(3)  =  FOURTH * (ONE + ETA)
      DN_BOT_DXI(4)  = -FOURTH * (ONE + ETA)

      DN_BOT_DETA(1) = -FOURTH * (ONE - XI)
      DN_BOT_DETA(2) = -FOURTH * (ONE + XI)
      DN_BOT_DETA(3) =  FOURTH * (ONE + XI)
      DN_BOT_DETA(4) =  FOURTH * (ONE - XI)

C=======================================================================
C   10. Combine top and bottom
C=======================================================================
C   Top surface nodes (NURBS). DN_LOCAL(I,J) = dN_I / d xi_J
      DO I = 1, N_TOP
        N(I) = N_TOP_ARY(I) * NZ_TOP
        DN_LOCAL(I,1) = DN_TOP_DU(I) * DU_DXI  * NZ_TOP
        DN_LOCAL(I,2) = DN_TOP_DV(I) * DV_DETA * NZ_TOP
        DN_LOCAL(I,3) = N_TOP_ARY(I) * DNZ_TOP
      ENDDO

C     Bottom corner nodes (order in IQ1NP_BULK_TAB: (-1,-1), (+1,-1),
C     (+1,+1), (-1,+1) to match N_BOT(1)..N_BOT(4), see genq1np Step 7c)
      DO I = 1, 4
        N(N_TOP + I) = N_BOT(I) * NZ_BOT
        DN_LOCAL(N_TOP + I,1) = DN_BOT_DXI(I)  * NZ_BOT
        DN_LOCAL(N_TOP + I,2) = DN_BOT_DETA(I) * NZ_BOT
        DN_LOCAL(N_TOP + I,3) = N_BOT(I)       * DNZ_BOT
      ENDDO

      RETURN ! Q1NP shape functions N and derivatives DN_LOCAL computation 
      END

C=======================================================================
C   Jacobian matrix and determinant computation
C=======================================================================
      SUBROUTINE Q1NP_JACOBIAN(DN_LOCAL, XNODE, NNODE, J, DETJ)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NNODE
      my_real, INTENT(IN) :: DN_LOCAL(NNODE,3), XNODE(3,NNODE)
      my_real, INTENT(OUT) :: J(3,3), DETJ
      ! NNODE: number of nodes
      ! DN_LOCAL: derivatives of shape functions
      ! XNODE: coordinates of nodes
      ! J: Jacobian matrix
      ! DETJ: determinant of Jacobian matrix
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K
C=======================================================================
C   1. Assemble Jacobian matrix: J(i,j) = sum_k XNODE(i,k) * DN_LOCAL(k,j)
C=======================================================================
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

      ENDDO

C=======================================================================
C   2. Compute determinant (3x3 explicit formula)
C=======================================================================
      DETJ = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) -
     .       J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) +
     .       J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))

      RETURN ! Jacobian matrix and determinant computed
      END



C=======================================================================
C   Main volume integration over element routine
C=======================================================================
      SUBROUTINE Q1NP_COMPUTE_VOLUME_ELEMENT(IEL_Q1NP,
     .                                       KQ1NP_TAB, IQ1NP_TAB,
     .                                       IQ1NP_BULK_TAB, Q1NP_KTAB,
     .                                       X, NX, NY, VOL_EL)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Q1NP_RESTART_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include "param_c.inc"
#include "com04_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: IEL_Q1NP, NX, NY
      INTEGER, INTENT(IN) :: KQ1NP_TAB(15,*), IQ1NP_TAB(*)
      INTEGER, INTENT(IN) :: IQ1NP_BULK_TAB(*)
      my_real, INTENT(IN) :: Q1NP_KTAB(*), X(3,*)
      my_real, INTENT(OUT) :: VOL_EL
      ! IEL_Q1NP: element index
      ! NX, NY: number of control points in U, V direction
      ! KQ1NP_TAB: element metadata
      ! IQ1NP_TAB: control point indices
      ! IQ1NP_BULK_TAB: bulk node indices
      ! Q1NP_KTAB: knot vector
      ! X: coordinates of nodes
      ! VOL_EL: volume of element
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: P, Q, NCTRL, OFFSET_CTRL, OFFSET_BULK
      INTEGER :: ELEM_U, ELEM_V
      INTEGER :: NKNOT_U, NKNOT_V
      INTEGER :: I, J, K, IU, IV, IT
      INTEGER :: N_TOP, N_TOTAL
      INTEGER, ALLOCATABLE :: NODE_IDS(:)
      my_real :: XI, ETA, ZETA, WG
      my_real, ALLOCATABLE :: NVAL(:)
      my_real, ALLOCATABLE :: DN_LOCAL(:,:)
      my_real, ALLOCATABLE :: XNODE(:,:)
      my_real, ALLOCATABLE :: U(:), V(:)
      my_real :: JMAT(3,3), DETJ
C=======================================================================
C   1. Extract element metadata from KQ1NP_TAB
C=======================================================================
      P           = KQ1NP_TAB(8, IEL_Q1NP)
      Q           = KQ1NP_TAB(9, IEL_Q1NP)
      NCTRL       = KQ1NP_TAB(3, IEL_Q1NP)
      OFFSET_CTRL = KQ1NP_TAB(4, IEL_Q1NP)
      OFFSET_BULK = KQ1NP_TAB(14,IEL_Q1NP)
      ELEM_U      = KQ1NP_TAB(6, IEL_Q1NP)
      ELEM_V      = KQ1NP_TAB(7, IEL_Q1NP)

C=======================================================================
C   2. Allocate arrays with correct sizes
C=======================================================================
      N_TOP   = NCTRL
      N_TOTAL = N_TOP + 4
      ALLOCATE(NODE_IDS(N_TOTAL))
      ALLOCATE(NVAL(N_TOTAL))
      ALLOCATE(DN_LOCAL(N_TOTAL,3))
      ALLOCATE(XNODE(3,N_TOTAL))

C=======================================================================
C   3. Allocate knot vectors
C=======================================================================
      NKNOT_U = NX + 2*P + 1
      NKNOT_V = NY + 2*Q + 1
      ALLOCATE(U(NKNOT_U))
      ALLOCATE(V(NKNOT_V))

C=======================================================================
C   4. Extract knot vectors from Q1NP_KTAB
C=======================================================================
      CALL Q1NP_GET_KNOT_VECTORS(NX, NY, P, Q, Q1NP_KTAB, U, V)
  
C=======================================================================
C   5. Collect node IDs and coordinates
C=======================================================================
C     Control points (top surface)
      DO I = 1, NCTRL
        NODE_IDS(I) = IQ1NP_TAB(OFFSET_CTRL + I - 1)
      ENDDO
      
C     Bulk nodes (bottom face)
      DO I = 1, 4
        NODE_IDS(N_TOP + I) = IQ1NP_BULK_TAB(OFFSET_BULK + I - 1)
      ENDDO

C     Extract coordinates
      DO K = 1, N_TOTAL
        XNODE(1,K) = X(1, NODE_IDS(K))
        XNODE(2,K) = X(2, NODE_IDS(K))
        XNODE(3,K) = X(3, NODE_IDS(K))
      ENDDO

C=======================================================================
C   6. Calculate volume of element using Gauss quadrature
C=======================================================================
      VOL_EL = ZERO

      DO IU = 1, Q1NP_NP_U_G ! 3 Gauss points in u direction
        XI = Q1NP_GP_U_G(IU)
        DO IV = 1, Q1NP_NP_V_G ! 3 Gauss points in v direction
          ETA = Q1NP_GP_V_G(IV)
          DO IT = 1, Q1NP_NP_T_G ! 2 Gauss points in t direction
            ZETA = Q1NP_GP_T_G(IT)
            WG   = Q1NP_GW_U_G(IU) * Q1NP_GW_V_G(IV) * Q1NP_GW_T_G(IT)

C           Calculate Q1NP shape functions and derivatives
            CALL Q1NP_SHAPE_FUNCTIONS(XI, ETA, ZETA, P, Q, U, V,
     .                                ELEM_U, ELEM_V, NVAL, DN_LOCAL)

C           Calculate Jacobian matrix and determinant
            CALL Q1NP_JACOBIAN(DN_LOCAL, XNODE, N_TOTAL, JMAT, DETJ)

C           Accumulate volume of element
            VOL_EL = VOL_EL + WG * DETJ
          ENDDO
        ENDDO
      ENDDO

C=======================================================================
C   7. Cleanup
C=======================================================================
      DEALLOCATE(U, V)
      DEALLOCATE(NODE_IDS, NVAL, DN_LOCAL, XNODE)

      RETURN ! Volume of element computed
      END
C=======================================================================
      END MODULE Q1NP_VOLUME_MOD
