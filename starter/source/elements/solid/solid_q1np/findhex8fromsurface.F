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
Chd|  FINDHEX8FROMSURF               source/elements/solid/solid_q1np/findhex8fromsurf.F
Chd|====================================================================
C=======================================================================
C   Map surface segments to HEX8 elements by matching top face nodes
C
C   This routine finds the HEX8 element that corresponds to a given
C   surface segment by matching the surface segment's 4 nodes to the
C   top face nodes of HEX8 elements.
C
C   Algorithm:
C   1. Loop through all HEX8 elements
C   2. Extract top face nodes (nodes 5-8) from each HEX8
C   3. Match surface segment nodes to HEX8 top face nodes
C   4. Use both node ID matching and coordinate-based matching
C      (with tolerance) to handle potential ID mismatches
C   5. Return the first HEX8 element where all 4 nodes match
C
C   HEX8 node ordering (standard):
C   - Bottom face (zeta=-1): nodes 1-4 (IXS indices 2-5)
C   - Top face (zeta=+1):    nodes 5-8 (IXS indices 6-9)
C=======================================================================
      SUBROUTINE FINDHEX8FROMSURF(ISEG,NODES_SURF,IXS,
     .                            X,IEL_HEX8,TOL)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "param_c.inc"
#include      "com04_c.inc"
#include      "elements.inc" 
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     ISEG        - surface segment index
C     NODES_SURF  - surface segment node IDs (4 nodes)
C     IXS         - HEX8 element connectivity array
C     X           - node coordinates array
C     IEL_HEX8    - HEX8 element index (output, 0 if not found)
C     TOL         - tolerance for coordinate matching
C     NUMELS and NUMNOD come from com04_c.inc
C-----------------------------------------------
      INTEGER ISEG,IEL_HEX8
      INTEGER NODES_SURF(4)
      INTEGER IXS(NIXS,*)
      my_real X(3,*)
      my_real TOL
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER IEL,I,J
      INTEGER NODES_HEX8(4)       ! Top face node IDs from HEX8
      INTEGER MATCH_COUNT         ! Number of matching nodes found
      my_real DX,DY,DZ            ! Coordinate differences
      my_real DSQ                 ! Squared distance
      my_real TOL_SQ              ! Squared tolerance for coordinate matching
C=======================================================================
C   Initialize: no element found yet
C=======================================================================
      IEL_HEX8 = 0
      TOL_SQ = TOL * TOL

C=======================================================================
C   Loop through all HEX8 elements to find matching top face
C=======================================================================
      DO IEL=1,NUMELS
C       Extract top face nodes from HEX8 (nodes 5-8, indices 6-9 in IXS)
        DO I=1,4
          NODES_HEX8(I) = IXS(6+I-1,IEL)
        ENDDO

C       Match surface segment nodes to HEX8 top face nodes
C       Try both node ID matching and coordinate-based matching
        MATCH_COUNT = 0
        DO I=1,4
          DO J=1,4
C           First try: exact node ID match
            IF (NODES_SURF(I) == NODES_HEX8(J)) THEN
              MATCH_COUNT = MATCH_COUNT + 1
              GOTO 100
            ENDIF
C           Second try: coordinate-based matching (for ID mismatches)
            IF (NODES_SURF(I) > 0 .AND. NODES_SURF(I) <= NUMNOD .AND.
     .          NODES_HEX8(J) > 0 .AND. NODES_HEX8(J) <= NUMNOD) THEN
              DX = X(1,NODES_SURF(I)) - X(1,NODES_HEX8(J))
              DY = X(2,NODES_SURF(I)) - X(2,NODES_HEX8(J))
              DZ = X(3,NODES_SURF(I)) - X(3,NODES_HEX8(J))
              DSQ = DX*DX + DY*DY + DZ*DZ
              IF (DSQ <= TOL_SQ) THEN
                MATCH_COUNT = MATCH_COUNT + 1
                GOTO 100
              ENDIF
            ENDIF
          ENDDO
  100     CONTINUE
        ENDDO

C       If all 4 nodes match, this is the corresponding HEX8 element
        IF (MATCH_COUNT >= 4) THEN
          IEL_HEX8 = IEL
          RETURN
        ENDIF
      ENDDO

C=======================================================================
C   No matching element found
C=======================================================================
      IEL_HEX8 = 0
      RETURN
      END

