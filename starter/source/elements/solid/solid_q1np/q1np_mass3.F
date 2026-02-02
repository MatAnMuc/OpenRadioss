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
Chd|  Q1NP_MASS3                    source/elements/solid/solid_q1np/q1np_mass3.F
Chd|====================================================================
C=======================================================================
C   Calculate lumped mass for Q1Np enriched elements
C
C   This routine distributes mass for Q1Np elements:
C   - Total mass equals the underlying HEX8 element mass
C   - 50% of mass distributed to 4 bulk nodes (12.5% each)
C   - 50% of mass distributed equally to all control point nodes
C   - HEX8 top nodes (nodes 5-8) explicitly set to zero mass
C
C   NOTE: This is a temporary implementation. Later, proper shape
C   function integration will be used for accurate mass distribution.
C=======================================================================
      SUBROUTINE Q1NP_MASS3(
     .     RHO, MS, PARTSAV, X, V,
     .     IPART, MSS, VOLU,
     .     MSNF, MSSF, IN,
     .     VR, INS, WMA, RHOCP, MCP,
     .     MCPS, MSSA, RHOF, FRAC, FILL,
     .     KQ1NP_TAB, IQ1NP_TAB, IQ1NP_BULK_TAB,
     .     IXS, NUMELQ1NP_IN, PM)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
      USE Q1NP_RESTART_MOD
      USE Q1NP_VOLUME_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "elements.inc" 
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "param_c.inc"
#include      "units_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NUMELQ1NP_IN
      INTEGER, INTENT(IN) :: KQ1NP_TAB(15,*)
      INTEGER, INTENT(IN) :: IQ1NP_TAB(*)
      INTEGER, INTENT(IN) :: IQ1NP_BULK_TAB(*)
      INTEGER, INTENT(IN) :: IXS(NIXS,*)
      INTEGER, INTENT(IN) :: IPART(*)
      my_real, INTENT(INOUT) :: MS(*), PARTSAV(20,*), X(3,*), V(3,*)
      my_real, INTENT(INOUT) :: MSS(8,*), MSSF(8,*), MSNF(*)
      my_real, INTENT(INOUT) :: IN(*), VR(3,*), INS(8,*), WMA(*)
      my_real, INTENT(INOUT) :: RHOCP(*), MCP(*), MCPS(8,*)
      my_real, INTENT(INOUT) :: MSSA(*)
      my_real, INTENT(IN) :: RHO(*), VOLU(*), RHOF(*), FRAC(*), FILL(*)
      my_real, INTENT(IN) :: PM(NPROPM,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER IEL, IEL_HEX8, IP
      INTEGER I, J, K, N_CTRL, OFFSET_CTRL, OFFSET_BULK
      INTEGER NODE_BULK(4), NODE_CP, NODE_HEX8_TOP(4)
      INTEGER MID
      my_real MASS_TOTAL, MASS_BULK, MASS_CP, MASS_PER_CP
      my_real MASS_PER_BULK
      my_real XX, YY, ZZ, XY, YZ, ZX
      my_real VOLU_Q1NP
      INTEGER, PARAMETER :: IDEBUG_Q1NP = 1
C=======================================================================
C   Early return if no Q1Np elements
C=======================================================================
      IF (NUMELQ1NP_IN == 0) RETURN

C=======================================================================
C   Loop over all Q1Np elements
C=======================================================================
      DO IEL = 1, NUMELQ1NP_IN
C       Get original HEX8 element index
        IEL_HEX8 = KQ1NP_TAB(10, IEL)
        IF (IEL_HEX8 <= 0 .OR. IEL_HEX8 > NUMELS) CYCLE

C       Get material ID and part ID
        MID = KQ1NP_TAB(1, IEL)
        IP = KQ1NP_TAB(11, IEL)
        IF (IP <= 0) CYCLE

C       Compute Q1NP volume exactly here (NURBS + bilinear integration).

        CALL Q1NP_COMPUTE_VOLUME_ELEMENT(IEL,
     .      KQ1NP_TAB, IQ1NP_TAB, IQ1NP_BULK_TAB,
     .      Q1NP_KTAB_G, X, Q1NP_NX_G, Q1NP_NY_G,
     .      VOLU_Q1NP)
C       write(*,*) 'IEL: ', IEL, 'VOLU_Q1NP: ', VOLU_Q1NP
        
C       Get density from PM array (not from GBUF%RHO which may not be initialized yet)
C       PM(1, MID) contains the material density
        IF (MID > 0) THEN
          MASS_TOTAL = PM(1, MID) * VOLU_Q1NP
C         Apply fill factor if available
          IF (FILL(IEL_HEX8) > ZERO) THEN
            MASS_TOTAL = FILL(IEL_HEX8) * MASS_TOTAL
          ENDIF
        ELSE
C         Fallback: use RHO array if MID is invalid (should not happen)
          IF (FILL(IEL_HEX8) > ZERO) THEN
            MASS_TOTAL = FILL(IEL_HEX8) * RHO(IEL_HEX8) * VOLU_Q1NP
          ELSE
            MASS_TOTAL = RHO(IEL_HEX8) * VOLU_Q1NP
          ENDIF
        ENDIF

C=======================================================================
C       Step 1: Distribute mass to bulk nodes and control nodes
C=======================================================================
C       Distribute 50% to bulk nodes, 50% to control nodes
        MASS_BULK = HALF * MASS_TOTAL
        MASS_CP = HALF * MASS_TOTAL

C=======================================================================
C       Step 2: Get bulk nodes (4 bottom nodes from HEX8)
C=======================================================================
        OFFSET_BULK = KQ1NP_TAB(14, IEL)
        DO I = 1, 4
          NODE_BULK(I) = IQ1NP_BULK_TAB(OFFSET_BULK + I - 1)
        ENDDO

C=======================================================================
C       Step 3: Distribute bulk mass equally (12.5% per node)
C=======================================================================
        MASS_PER_BULK = MASS_BULK * FOURTH
        
        DO I = 1, 4
          IF (NODE_BULK(I) > 0 .AND. NODE_BULK(I) <= NUMNOD) THEN
            MS(NODE_BULK(I)) = MS(NODE_BULK(I)) + MASS_PER_BULK
          ENDIF
        ENDDO
        
C=======================================================================
C       Step 4: Get control point nodes and distribute mass to control nodes
C=======================================================================
        N_CTRL = KQ1NP_TAB(3, IEL)
        OFFSET_CTRL = KQ1NP_TAB(4, IEL)
        IF (N_CTRL > 0) THEN
          MASS_PER_CP = MASS_CP / REAL(N_CTRL)
          DO I = 1, N_CTRL
            NODE_CP = IQ1NP_TAB(OFFSET_CTRL + I - 1)
            IF (NODE_CP > 0 .AND. NODE_CP <= NUMNOD) THEN
              MS(NODE_CP) = MS(NODE_CP) + MASS_PER_CP
            ENDIF
          ENDDO
        ENDIF

C=======================================================================
C       Step 5: Update PARTSAV (part-level mass and inertia)
C=======================================================================
C       Note: SMASS3 is skipped for replaced elements, so we just add the Q1Np mass
        PARTSAV(1,IP) = PARTSAV(1,IP) + MASS_TOTAL

        DO I = 1, 4
          IF (NODE_BULK(I) > 0 .AND. NODE_BULK(I) <= NUMNOD) THEN
            PARTSAV(2,IP) = PARTSAV(2,IP) + MASS_PER_BULK * X(1, NODE_BULK(I))
            PARTSAV(3,IP) = PARTSAV(3,IP) + MASS_PER_BULK * X(2, NODE_BULK(I))
            PARTSAV(4,IP) = PARTSAV(4,IP) + MASS_PER_BULK * X(3, NODE_BULK(I))
          ENDIF
        ENDDO
        DO I = 1, N_CTRL
          NODE_CP = IQ1NP_TAB(OFFSET_CTRL + I - 1)
          IF (NODE_CP > 0 .AND. NODE_CP <= NUMNOD) THEN
            PARTSAV(2,IP) = PARTSAV(2,IP) + MASS_PER_CP * X(1, NODE_CP)
            PARTSAV(3,IP) = PARTSAV(3,IP) + MASS_PER_CP * X(2, NODE_CP)
            PARTSAV(4,IP) = PARTSAV(4,IP) + MASS_PER_CP * X(3, NODE_CP)
          ENDIF
        ENDDO

C=======================================================================
C       Step 6: Calculate inertia contributions (simplified - using nodal positions)
C=======================================================================
        XX = ZERO
        YY = ZERO
        ZZ = ZERO
        XY = ZERO
        YZ = ZERO
        ZX = ZERO
        DO I = 1, 4
          IF (NODE_BULK(I) > 0 .AND. NODE_BULK(I) <= NUMNOD) THEN
            XX = XX + MASS_PER_BULK * X(1, NODE_BULK(I)) * X(1, NODE_BULK(I))
            YY = YY + MASS_PER_BULK * X(2, NODE_BULK(I)) * X(2, NODE_BULK(I))
            ZZ = ZZ + MASS_PER_BULK * X(3, NODE_BULK(I)) * X(3, NODE_BULK(I))
            XY = XY + MASS_PER_BULK * X(1, NODE_BULK(I)) * X(2, NODE_BULK(I))
            YZ = YZ + MASS_PER_BULK * X(2, NODE_BULK(I)) * X(3, NODE_BULK(I))
            ZX = ZX + MASS_PER_BULK * X(3, NODE_BULK(I)) * X(1, NODE_BULK(I))
          ENDIF
        ENDDO
        DO I = 1, N_CTRL
          NODE_CP = IQ1NP_TAB(OFFSET_CTRL + I - 1)
          IF (NODE_CP > 0 .AND. NODE_CP <= NUMNOD) THEN
            XX = XX + MASS_PER_CP * X(1, NODE_CP) * X(1, NODE_CP)
            YY = YY + MASS_PER_CP * X(2, NODE_CP) * X(2, NODE_CP)
            ZZ = ZZ + MASS_PER_CP * X(3, NODE_CP) * X(3, NODE_CP)
            XY = XY + MASS_PER_CP * X(1, NODE_CP) * X(2, NODE_CP)
            YZ = YZ + MASS_PER_CP * X(2, NODE_CP) * X(3, NODE_CP)
            ZX = ZX + MASS_PER_CP * X(3, NODE_CP) * X(1, NODE_CP)
          ENDIF
        ENDDO
        PARTSAV(5,IP) = PARTSAV(5,IP) + (YY + ZZ)
        PARTSAV(6,IP) = PARTSAV(6,IP) + (ZZ + XX)
        PARTSAV(7,IP) = PARTSAV(7,IP) + (XX + YY)
        PARTSAV(8,IP) = PARTSAV(8,IP) - XY
        PARTSAV(9,IP) = PARTSAV(9,IP) - YZ
        PARTSAV(10,IP) = PARTSAV(10,IP) - ZX

C=======================================================================
C       Step 7: Calculate momentum contributions
C=======================================================================
        DO I = 1, 4
          IF (NODE_BULK(I) > 0 .AND. NODE_BULK(I) <= NUMNOD) THEN
            PARTSAV(11,IP) = PARTSAV(11,IP) + MASS_PER_BULK * V(1, NODE_BULK(I))
            PARTSAV(12,IP) = PARTSAV(12,IP) + MASS_PER_BULK * V(2, NODE_BULK(I))
            PARTSAV(13,IP) = PARTSAV(13,IP) + MASS_PER_BULK * V(3, NODE_BULK(I))
          ENDIF
        ENDDO
        DO I = 1, N_CTRL
          NODE_CP = IQ1NP_TAB(OFFSET_CTRL + I - 1)
          IF (NODE_CP > 0 .AND. NODE_CP <= NUMNOD) THEN
            PARTSAV(11,IP) = PARTSAV(11,IP) + MASS_PER_CP * V(1, NODE_CP)
            PARTSAV(12,IP) = PARTSAV(12,IP) + MASS_PER_CP * V(2, NODE_CP)
            PARTSAV(13,IP) = PARTSAV(13,IP) + MASS_PER_CP * V(3, NODE_CP)
          ENDIF
        ENDDO

C=======================================================================
C       Step 8: Calculate kinetic energy contribution
C=======================================================================
        DO I = 1, 4
          IF (NODE_BULK(I) > 0 .AND. NODE_BULK(I) <= NUMNOD) THEN
            PARTSAV(14,IP) = PARTSAV(14,IP) + HALF * MASS_PER_BULK *
     .        (V(1, NODE_BULK(I))*V(1, NODE_BULK(I)) +
     .         V(2, NODE_BULK(I))*V(2, NODE_BULK(I)) +
     .         V(3, NODE_BULK(I))*V(3, NODE_BULK(I)))
          ENDIF
        ENDDO
        DO I = 1, N_CTRL
          NODE_CP = IQ1NP_TAB(OFFSET_CTRL + I - 1)
          IF (NODE_CP > 0 .AND. NODE_CP <= NUMNOD) THEN
            PARTSAV(14,IP) = PARTSAV(14,IP) + HALF * MASS_PER_CP *
     .        (V(1, NODE_CP)*V(1, NODE_CP) +
     .         V(2, NODE_CP)*V(2, NODE_CP) +
     .         V(3, NODE_CP)*V(3, NODE_CP))
          ENDIF
        ENDDO

C=======================================================================
C       Step 9: Print optional debug information
C=======================================================================
        If (IDEBUG_Q1NP > 0) THEN
          PRINT*, ' '
          PRINT*, 'Q1NP DEBUG: Element ', IEL, ' (Global ID: ', IXS(NIXS,IEL_HEX8), ')'
          PRINT*, '  MASS_TOTAL = ', MASS_TOTAL
          PRINT*, '  MASS_BULK  = ', MASS_BULK
          PRINT*, '  MASS_CP    = ', MASS_CP
          PRINT*, '  Bulk nodes:'
          DO I = 1, 4
            IF (NODE_BULK(I) > 0 .AND. NODE_BULK(I) <= NUMNOD) THEN
              PRINT*, '    Node ', NODE_BULK(I), ' MS = ', MS(NODE_BULK(I))
            ENDIF
          ENDDO
          PRINT*, '  Control nodes (N_CTRL = ', N_CTRL, ')'
          DO I = 1, N_CTRL
            NODE_CP = IQ1NP_TAB(OFFSET_CTRL + I - 1)
            IF (NODE_CP > 0 .AND. NODE_CP <= NUMNOD) THEN
              PRINT*, '    Node ', NODE_CP, ' MS = ', MS(NODE_CP)
            ENDIF
          ENDDO
        ENDIF

      ENDDO

      RETURN
      END

C=======================================================================
C   Helper routine to calculate HEX8 volume from coordinates (unused)
C=======================================================================
      SUBROUTINE Q1NP_HEX8_VOLUME(IXS_ELEM, X, VOLU_OUT)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "elements.inc" 
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: IXS_ELEM(NIXS)
      my_real, INTENT(IN) :: X(3,*)  ! 2D array: X(1,N), X(2,N), X(3,N) for node N
      my_real, INTENT(OUT) :: VOLU_OUT
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8
      my_real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4
      my_real X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7, X8, Y8, Z8
      my_real V1(3), V2(3), V3(3), V4(3), V5(3), V6(3)
      my_real DET1, DET2, DET3, DET4, DET5, DET6
C=======================================================================
C   Get node IDs (IXS indices: 2-9 are nodes 1-8)
      I1 = IXS_ELEM(2)
      I2 = IXS_ELEM(3)
      I3 = IXS_ELEM(4)
      I4 = IXS_ELEM(5)
      I5 = IXS_ELEM(6)
      I6 = IXS_ELEM(7)
      I7 = IXS_ELEM(8)
      I8 = IXS_ELEM(9)

C   Get coordinates from 2D array
C   Node N coordinates: X(1,N), X(2,N), X(3,N)
      X1 = X(1, I1); Y1 = X(2, I1); Z1 = X(3, I1)
      X2 = X(1, I2); Y2 = X(2, I2); Z2 = X(3, I2)
      X3 = X(1, I3); Y3 = X(2, I3); Z3 = X(3, I3)
      X4 = X(1, I4); Y4 = X(2, I4); Z4 = X(3, I4)
      X5 = X(1, I5); Y5 = X(2, I5); Z5 = X(3, I5)
      X6 = X(1, I6); Y6 = X(2, I6); Z6 = X(3, I6)
      X7 = X(1, I7); Y7 = X(2, I7); Z7 = X(3, I7)
      X8 = X(1, I8); Y8 = X(2, I8); Z8 = X(3, I8)

C   Calculate volume using 5-tetrahedron decomposition
C   Split hex into 5 tetrahedra and sum volumes
C   Tetrahedron volume = |det(v1, v2, v3)| / 6
C   where v1, v2, v3 are edge vectors from one vertex

C   Tetrahedron 1: (1,2,3,5)
      V1(1) = X2 - X1; V1(2) = Y2 - Y1; V1(3) = Z2 - Z1
      V2(1) = X3 - X1; V2(2) = Y3 - Y1; V2(3) = Z3 - Z1
      V3(1) = X5 - X1; V3(2) = Y5 - Y1; V3(3) = Z5 - Z1
      DET1 = V1(1)*(V2(2)*V3(3) - V2(3)*V3(2)) -
     .       V1(2)*(V2(1)*V3(3) - V2(3)*V3(1)) +
     .       V1(3)*(V2(1)*V3(2) - V2(2)*V3(1))
      DET1 = ABS(DET1) / SIX

C   Tetrahedron 2: (1,3,4,5)
      V1(1) = X3 - X1; V1(2) = Y3 - Y1; V1(3) = Z3 - Z1
      V2(1) = X4 - X1; V2(2) = Y4 - Y1; V2(3) = Z4 - Z1
      V3(1) = X5 - X1; V3(2) = Y5 - Y1; V3(3) = Z5 - Z1
      DET2 = V1(1)*(V2(2)*V3(3) - V2(3)*V3(2)) -
     .       V1(2)*(V2(1)*V3(3) - V2(3)*V3(1)) +
     .       V1(3)*(V2(1)*V3(2) - V2(2)*V3(1))
      DET2 = ABS(DET2) / SIX

C   Tetrahedron 3: (1,4,2,6)
      V1(1) = X4 - X1; V1(2) = Y4 - Y1; V1(3) = Z4 - Z1
      V2(1) = X2 - X1; V2(2) = Y2 - Y1; V2(3) = Z2 - Z1
      V3(1) = X6 - X1; V3(2) = Y6 - Y1; V3(3) = Z6 - Z1
      DET3 = V1(1)*(V2(2)*V3(3) - V2(3)*V3(2)) -
     .       V1(2)*(V2(1)*V3(3) - V2(3)*V3(1)) +
     .       V1(3)*(V2(1)*V3(2) - V2(2)*V3(1))
      DET3 = ABS(DET3) / SIX

C   Tetrahedron 4: (1,5,6,7)
      V1(1) = X5 - X1; V1(2) = Y5 - Y1; V1(3) = Z5 - Z1
      V2(1) = X6 - X1; V2(2) = Y6 - Y1; V2(3) = Z6 - Z1
      V3(1) = X7 - X1; V3(2) = Y7 - Y1; V3(3) = Z7 - Z1
      DET4 = V1(1)*(V2(2)*V3(3) - V2(3)*V3(2)) -
     .       V1(2)*(V2(1)*V3(3) - V2(3)*V3(1)) +
     .       V1(3)*(V2(1)*V3(2) - V2(2)*V3(1))
      DET4 = ABS(DET4) / SIX

C   Tetrahedron 5: (1,7,8,4)
      V1(1) = X7 - X1; V1(2) = Y7 - Y1; V1(3) = Z7 - Z1
      V2(1) = X8 - X1; V2(2) = Y8 - Y1; V2(3) = Z8 - Z1
      V3(1) = X4 - X1; V3(2) = Y4 - Y1; V3(3) = Z4 - Z1
      DET5 = V1(1)*(V2(2)*V3(3) - V2(3)*V3(2)) -
     .       V1(2)*(V2(1)*V3(3) - V2(3)*V3(1)) +
     .       V1(3)*(V2(1)*V3(2) - V2(2)*V3(1))
      DET5 = ABS(DET5) / SIX

C   Tetrahedron 6: (1,5,7,8)
      V1(1) = X5 - X1; V1(2) = Y5 - Y1; V1(3) = Z5 - Z1
      V2(1) = X7 - X1; V2(2) = Y7 - Y1; V2(3) = Z7 - Z1
      V3(1) = X8 - X1; V3(2) = Y8 - Y1; V3(3) = Z8 - Z1
      DET6 = V1(1)*(V2(2)*V3(3) - V2(3)*V3(2)) -
     .       V1(2)*(V2(1)*V3(3) - V2(3)*V3(1)) +
     .       V1(3)*(V2(1)*V3(2) - V2(2)*V3(1))
      DET6 = ABS(DET6) / SIX

C   Total volume = sum of all tetrahedra
      VOLU_OUT = DET1 + DET2 + DET3 + DET4 + DET5 + DET6

      RETURN
      END

