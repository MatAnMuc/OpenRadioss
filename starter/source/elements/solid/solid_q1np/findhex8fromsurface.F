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
C   Map surface segments to HEX8 elements by matching any face
C
C   This routine finds the HEX8 element that corresponds to a given
C   surface segment by matching the segment 4 nodes to any of the
C   six faces of HEX8 elements. Returns element index and the four
C   node IDs of the opposite face (bulk side). 
C
C   Algorithm Overview:
C     - Iterate through all HEX8 elements; for each, check all six faces.
C     - For each face (4 nodes), match the input surface segment’s four nodes by comparing node IDs
C     - When all four match a face, return element index and fill NODES_BULK
C       with the 4 nodes of the opposite face. If no match, IEL_HEX8 = 0.
C   HEX8 face convention (IXS(2:9)=nodes 1-8): Face 1: 1,2,3,4; Face 2: 5,6,7,8;
C   Face 3: 1,2,6,5; Face 4: 2,3,7,6; Face 5: 3,4,8,7; Face 6: 4,1,5,8.
C
C=======================================================================
      SUBROUTINE FINDHEX8FROMSURF(ISEG,NODES_SURF,IXS,
     .                            X,IEL_HEX8,NODES_BULK,TOL)
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
C     NODES_BULK  - output: 4 node IDs of the face opposite to matched face
C     TOL         - tolerance for coordinate matching
C-----------------------------------------------
      INTEGER ISEG,IEL_HEX8
      INTEGER NODES_SURF(4)
      INTEGER NODES_BULK(4)
      INTEGER IXS(NIXS,*)
      my_real X(3,*)
      my_real TOL
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER IEL,I,J,IFACE,IFOPP
      INTEGER NODES_HEX8(4)       ! Current face node IDs from HEX8
      INTEGER MATCH_COUNT         ! Number of matching nodes found
      INTEGER FACE_IXS(4,6)       ! IXS row index for each face's 4 nodes
      INTEGER OPPOSITE(6)         ! Opposite face index for each face
      my_real DX,DY,DZ
      my_real DSQ
      my_real TOL_SQ
C     HEX8 faces: IXS(2:9)=nodes 1-8. Face 1: 1,2,3,4; Face 2: 5,6,7,8;
C     Face 3: 1,2,6,5; Face 4: 2,3,7,6; Face 5: 3,4,8,7; Face 6: 4,1,5,8
      DATA FACE_IXS /
     .  2,3,4,5, 6,7,8,9, 2,3,7,6, 3,4,8,7, 4,5,9,8, 5,2,6,9 /
      DATA OPPOSITE / 2,1,5,6,3,4 /
C=======================================================================
C   1. Initialize: no element found yet
C=======================================================================
      IEL_HEX8 = 0
      TOL_SQ = TOL * TOL

C=======================================================================
C   2. Loop through all HEX8 elements and all six faces
C=======================================================================
      DO IEL=1,NUMELS
        DO IFACE=1,6
C         Extract the 4 nodes of this face
          DO I=1,4
            NODES_HEX8(I) = IXS(FACE_IXS(I,IFACE),IEL)
          ENDDO

C         Match surface segment nodes to this face (ID)
          MATCH_COUNT = 0
          DO I=1,4
            DO J=1,4
              IF (NODES_SURF(I) == NODES_HEX8(J)) THEN
                MATCH_COUNT = MATCH_COUNT + 1
                GOTO 100
              ENDIF
            ENDDO
C           No match for this surface node: skip to next face
            GOTO 101
  100       CONTINUE
          ENDDO

          IF (MATCH_COUNT >= 4) THEN
            IEL_HEX8 = IEL
            IFOPP = OPPOSITE(IFACE)
            DO I=1,4
              NODES_BULK(I) = IXS(FACE_IXS(I,IFOPP),IEL)
            ENDDO
            RETURN
          ENDIF
  101     CONTINUE
        ENDDO
      ENDDO

C=======================================================================
C   No matching element found
C=======================================================================
      IEL_HEX8 = 0
      RETURN
      END

