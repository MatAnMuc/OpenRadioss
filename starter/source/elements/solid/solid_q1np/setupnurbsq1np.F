Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
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
Chd|  SETUPNURBSQ1NP                 source/elements/solid/solid_q1np/setupnurbsq1np.F
Chd|====================================================================
C=======================================================================
C   Set up NURBS knot vectors and weights for Q1Np elements
C
C   This routine generates open uniform knot vectors following the
C   algorithm from element_lib.py:open_uniform_knot_vector.
C
C   Algorithm for open uniform knot vectors:
C   - For ne elements and degree p, the knot vector has length: ne + 2*p + 1
C   - Structure: (p+1) zeros, interior knots, (p+1) ones
C   - Interior knots: evenly spaced from 0 to 1
C   - Number of control points: ncp = ne + p
C
C   Example for ne=2, p=2:
C     U = [0, 0, 0, 0.5, 1, 1, 1]
C     (3 zeros, 1 interior knot at 0.5, 3 ones)
C
C   The routine generates:
C   1. U knot vector (for u direction, length = NX + 2*P + 1)
C   2. V knot vector (for v direction, length = NY + 2*Q + 1)
C   3. Weights array (all 1.0 for non-rational B-splines)
C
C   Knot vectors are stored concatenated: U knots first, then V knots
C=======================================================================
      SUBROUTINE SETUPNURBSQ1NP(NX,NY,P,Q,
     .                          Q1NP_KTAB,NCP_U,NCP_V,
     .                          Q1NP_WTAB,NWEIGHT_MAX)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     NX,NY     - grid dimensions (number of elements in u,v directions)
C     P,Q       - NURBS degrees (hardcoded: p=2, q=2 for testing)
C     Q1NP_KTAB - knot vectors array (output: U knots then V knots)
C     NCP_U     - number of control points in u direction (output)
C     NCP_V     - number of control points in v direction (output)
C     Q1NP_WTAB - weights array (output: all 1.0 for non-rational)
C     NWEIGHT_MAX- maximum number of weights available in Q1NP_WTAB
C-----------------------------------------------
      INTEGER NX,NY,P,Q
      INTEGER NCP_U,NCP_V
      INTEGER NWEIGHT_MAX
      my_real Q1NP_KTAB(*),Q1NP_WTAB(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I
      INTEGER NE_U,NE_V          ! Number of elements in u,v directions
      INTEGER NKNOT_U,NKNOT_V    ! Knot vector lengths
      INTEGER IV_OFFSET          ! Offset for V knot vector in Q1NP_KTAB
      my_real DELTA_U,DELTA_V    ! Spacing for interior knots
C=======================================================================
C   Calculate NURBS parameters
C=======================================================================
      NE_U = NX
      NE_V = NY
C     Number of control points: ncp = ne + p
      NCP_U = NE_U + P
      NCP_V = NE_V + Q
C     Knot vector length: len(U) = ne + 2*p + 1
      NKNOT_U = NE_U + 2*P + 1
      NKNOT_V = NE_V + 2*Q + 1

C=======================================================================
C   Generate U knot vector: (p+1) zeros, interior knots, (p+1) ones
C=======================================================================
C     Leading zeros: (p+1) zeros at start
      DO I=1,P+1
        Q1NP_KTAB(I) = ZERO
      ENDDO

C     Interior knots: evenly spaced from 0 to 1
C     For ne elements, interior knots are at: 1/ne, 2/ne, ..., (ne-1)/ne
      IF (NE_U > 1) THEN
        DELTA_U = ONE / REAL(NE_U)
        DO I=1,NE_U-1
          Q1NP_KTAB(P+1+I) = REAL(I) * DELTA_U
        ENDDO
      ENDIF

C     Trailing ones: (p+1) ones at end
      DO I=1,P+1
        Q1NP_KTAB(P+1+NE_U+I-1) = ONE
      ENDDO

C=======================================================================
C   Generate V knot vector: (q+1) zeros, interior knots, (q+1) ones
C=======================================================================
C     V knot vector starts after U knot vector
      IV_OFFSET = NKNOT_U

C     Leading zeros: (q+1) zeros at start
      DO I=1,Q+1
        Q1NP_KTAB(IV_OFFSET+I) = ZERO
      ENDDO

C     Interior knots: evenly spaced from 0 to 1
      IF (NE_V > 1) THEN
        DELTA_V = ONE / REAL(NE_V)
        DO I=1,NE_V-1
          Q1NP_KTAB(IV_OFFSET+Q+1+I) = REAL(I) * DELTA_V
        ENDDO
      ENDIF

C     Trailing ones: (q+1) ones at end
      DO I=1,Q+1
        Q1NP_KTAB(IV_OFFSET+Q+1+NE_V+I-1) = ONE
      ENDDO

C=======================================================================
C   Initialize weights to 1.0 (non-rational B-splines)
C=======================================================================
C     All weights are set to 1.0, making these non-rational B-splines
C     rather than full NURBS. For rational surfaces, weights would vary.
      DO I=1,MIN(NCP_U*NCP_V,NWEIGHT_MAX)
        Q1NP_WTAB(I) = ONE
      ENDDO

      RETURN
      END

