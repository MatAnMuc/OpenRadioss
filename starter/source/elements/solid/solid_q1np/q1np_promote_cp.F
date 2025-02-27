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
Copyright>        Commercial Alternative: Altair Engineering Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
Chd|====================================================================
Chd|  Q1NP_PROMOTE_CP_TO_NODES        source/elements/solid/solid_q1np/q1np_promote_cp.F
Chd|====================================================================
C=======================================================================
C   Promote Q1Np control points to real mesh nodes
C
C   This routine converts control points stored in Q1NP_CPTAB into
C   real mesh nodes by:
C   1. Assigning unique global node IDs starting from NUMNOD+1
C   2. Copying coordinates to X array
C   3. Initializing velocities and other nodal arrays
C   4. Updating IQ1NP_TAB to use the new node IDs
C   5. Tracking which HEX8 top nodes should have zero mass
C
C
C=======================================================================
      SUBROUTINE Q1NP_PROMOTE_CP_TO_NODES(
     .                   NUMNOD_IN, NUMNOD_OUT,
     .                   X, V, D, MS, IN, VR, DR,
     .                   Q1NP_CPTAB, IQ1NP_TAB, IQ1NP_BULK_TAB, KQ1NP_TAB,
     .                   IXS, NUMELQ1NP_IN, NUMNOD_CP_ADDED,
     .                   SX, SV, SD, SMS, SIN, SVR, SDR)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
      USE Q1NP_RESTART_MOD, ONLY: SQ1NPCTRL_SHARED_G
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "param_c.inc"
#include      "units_c.inc"
#include      "elements.inc" 
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NUMNOD_IN
      INTEGER, INTENT(OUT) :: NUMNOD_OUT
      INTEGER, INTENT(IN) :: NUMELQ1NP_IN
      INTEGER, INTENT(OUT) :: NUMNOD_CP_ADDED
      INTEGER, INTENT(INOUT) :: IQ1NP_TAB(*)
      INTEGER, INTENT(IN) :: IQ1NP_BULK_TAB(*)
      INTEGER, INTENT(INOUT) :: KQ1NP_TAB(15,*)
      INTEGER, INTENT(IN) :: IXS(NIXS,*)
      my_real, INTENT(INOUT) :: X(*), V(*), D(*), MS(*), IN(*)
      my_real, INTENT(INOUT) :: VR(*), DR(*)
      my_real, INTENT(IN) :: Q1NP_CPTAB(3,*)
C
C     Explicit sizes of nodal arrays (for debug / bounds checking)
      INTEGER, INTENT(IN) :: SX, SV, SD, SMS, SIN, SVR, SDR
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: IDEBUG_Q1NP = 0
      INTEGER I, IEL, CP_IDX, CP_IDX_OLD, NEW_NODE_ID
      INTEGER OFFSET_CTRL, N_CTRL
      INTEGER, DIMENSION(:), ALLOCATABLE :: CP_TO_NODE_MAP
      INTEGER MAX_CP_ID, MAX_NEW_NODE
      INTEGER REQ_SX, REQ_SV, REQ_SD, REQ_SMS, REQ_SIN, REQ_SVR, REQ_SDR

C=======================================================================
C   Step 1: Determine number of unique control points
C=======================================================================
      
      IF (NUMELQ1NP_IN == 0 .OR. SQ1NPCTRL_SHARED_G == 0) THEN
        NUMNOD_OUT = NUMNOD_IN
        NUMNOD_CP_ADDED = 0
        RETURN
      ENDIF

      MAX_CP_ID = SQ1NPCTRL_SHARED_G
      MAX_NEW_NODE = NUMNOD_IN + MAX_CP_ID

C=======================================================================
C   Debug bounds checks for nodal arrays
C   Ensure that the maximum indices used later in this routine are
C   within the declared sizes passed from the caller.
C=======================================================================
      REQ_SX   = 3*(MAX_NEW_NODE-1) + 3
      REQ_SV   = REQ_SX
      REQ_SD   = 5*MAX_NEW_NODE
      REQ_SMS  = MAX_NEW_NODE
      REQ_SIN  = MAX_NEW_NODE
      REQ_SVR  = 3*MAX_NEW_NODE
      REQ_SDR  = 3*MAX_NEW_NODE

      IF (REQ_SX > SX .OR. REQ_SV > SV .OR.
     .    REQ_SD > SD .OR.
     .    (REQ_SMS > SMS .AND. SMS > 0) .OR.
     .    (REQ_SIN  > SIN  .AND. SIN  > 0) .OR.
     .    (REQ_SVR > SVR .AND. SVR > 0) .OR.
     .    (REQ_SDR > SDR .AND. SDR > 0)) THEN
        WRITE(IOUT,'(A)') 
     .  ' Q1Np ERROR: Q1NP_PROMOTE_CP_TO_NODES array bounds check failed'
        WRITE(IOUT,'(A,3I12)') 
     .  '  NUMNOD_IN, MAX_CP_ID, MAX_NEW_NODE =',
     .   NUMNOD_IN, MAX_CP_ID, MAX_NEW_NODE
        WRITE(IOUT,'(A,7I12)')
     .  '  Required sizes SX,SV,SD,SMS,SIN,SVR,SDR =',
     .   REQ_SX,REQ_SV,REQ_SD,REQ_SMS,REQ_SIN,REQ_SVR,REQ_SDR
        WRITE(IOUT,'(A,7I12)')
     .  '  Actual   sizes SX,SV,SD,SMS,SIN,SVR,SDR =',
     .   SX,SV,SD,SMS,SIN,SVR,SDR
C       Also print to standard output for interactive debugging (guarded)
        IF (IDEBUG_Q1NP > 0) THEN
          WRITE(*,'(A)') 
     .    ' Q1Np ERROR: Q1NP_PROMOTE_CP_TO_NODES array bounds check failed'
          WRITE(*,'(A,3I12)') 
     .    '  NUMNOD_IN, MAX_CP_ID, MAX_NEW_NODE =',
     .     NUMNOD_IN, MAX_CP_ID, MAX_NEW_NODE
          WRITE(*,'(A,7I12)')
     .    '  Required sizes SX,SV,SD,SMS,SIN,SVR,SDR =',
     .     REQ_SX,REQ_SV,REQ_SD,REQ_SMS,REQ_SIN,REQ_SVR,REQ_SDR
          WRITE(*,'(A,7I12)')
     .    '  Actual   sizes SX,SV,SD,SMS,SIN,SVR,SDR =',
     .     SX,SV,SD,SMS,SIN,SVR,SDR
        ENDIF
        CALL ANCMSG(MSGID=402,
     .              MSGTYPE=MSGERROR,
     .              ANMODE=ANINFO,
     .              C1='Q1NP_PROMOTE_CP_TO_NODES: array bounds')
        NUMNOD_OUT = NUMNOD_IN
        NUMNOD_CP_ADDED = 0
        RETURN
      ENDIF

      ALLOCATE(CP_TO_NODE_MAP(MAX_CP_ID), STAT=I)
      IF (I /= 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .              MSGTYPE=MSGERROR,
     .              C1='CP_TO_NODE_MAP')
        RETURN
      ENDIF
      CP_TO_NODE_MAP = 0

C=======================================================================
C   Step 2: Assign new node IDs to each control point
C=======================================================================
      NUMNOD_CP_ADDED = 0
      DO CP_IDX = 1, MAX_CP_ID
        NEW_NODE_ID = NUMNOD_IN + CP_IDX
        CP_TO_NODE_MAP(CP_IDX) = NEW_NODE_ID
        NUMNOD_CP_ADDED = NUMNOD_CP_ADDED + 1
      ENDDO

C=======================================================================
C   Step 3: Copy control point coordinates to Node array X
C=======================================================================
      DO CP_IDX = 1, MAX_CP_ID
        NEW_NODE_ID = CP_TO_NODE_MAP(CP_IDX)
        X(3*(NEW_NODE_ID-1)+1) = Q1NP_CPTAB(1, CP_IDX)
        X(3*(NEW_NODE_ID-1)+2) = Q1NP_CPTAB(2, CP_IDX)
        X(3*(NEW_NODE_ID-1)+3) = Q1NP_CPTAB(3, CP_IDX)

C       Initialize kinematic arrays to zero (only if present)
        V(3*(NEW_NODE_ID-1)+1) = ZERO
        V(3*(NEW_NODE_ID-1)+2) = ZERO
        V(3*(NEW_NODE_ID-1)+3) = ZERO
        D(3*(NEW_NODE_ID-1)+1) = ZERO
        D(3*(NEW_NODE_ID-1)+2) = ZERO
        D(3*(NEW_NODE_ID-1)+3) = ZERO
        IF (SMS > 0) MS(NEW_NODE_ID) = ZERO
        IF (SIN > 0) IN(NEW_NODE_ID) = ZERO
        IF (SVR > 0) THEN
          VR(3*(NEW_NODE_ID-1)+1) = ZERO
          VR(3*(NEW_NODE_ID-1)+2) = ZERO
          VR(3*(NEW_NODE_ID-1)+3) = ZERO
        ENDIF
        IF (SDR > 0) THEN
          DR(3*(NEW_NODE_ID-1)+1) = ZERO
          DR(3*(NEW_NODE_ID-1)+2) = ZERO
          DR(3*(NEW_NODE_ID-1)+3) = ZERO
        ENDIF
      ENDDO

C=======================================================================
C   Debug: Print control point node coordinates in X array
C=======================================================================
C      IF (IDEBUG_Q1NP > 0) THEN
C        WRITE(*,*)
C        WRITE(*,'(A)') 'Q1NP DEBUG: Control point node coordinates in X array'
C        DO CP_IDX = 1, MAX_CP_ID
C          NEW_NODE_ID = CP_TO_NODE_MAP(CP_IDX)
C          WRITE(*,'(A,I6,A,I10,A,3(1X,1PE15.8))') 
C     .      '  CP(',CP_IDX,') -> node ',NEW_NODE_ID,' :',
C     .      X(3*(NEW_NODE_ID-1)+1),
C     .      X(3*(NEW_NODE_ID-1)+2),
C     .      X(3*(NEW_NODE_ID-1)+3)
C        ENDDO
C      ENDIF

C=======================================================================
C   Step 4: Update IQ1NP_TAB to use new node IDs
C=======================================================================
      DO IEL = 1, NUMELQ1NP_IN
        N_CTRL = KQ1NP_TAB(3, IEL)
        OFFSET_CTRL = KQ1NP_TAB(4, IEL)
        DO I = 1, N_CTRL
          CP_IDX_OLD = IQ1NP_TAB(OFFSET_CTRL + I - 1)
          IF (CP_IDX_OLD > 0 .AND. CP_IDX_OLD <= MAX_CP_ID) THEN
            IQ1NP_TAB(OFFSET_CTRL + I - 1) = CP_TO_NODE_MAP(CP_IDX_OLD)
          ENDIF
        ENDDO
      ENDDO

C=======================================================================
C   Debug: Print control point to node map
C=======================================================================
      IF (IDEBUG_Q1NP > 0) THEN
        WRITE(*,*)
        WRITE(*,'(A)') '=============================================='
        WRITE(*,'(A)') 'Control point to FEM node map'
        WRITE(*,'(A)') '=============================================='
        DO CP_IDX = 1, MAX_CP_ID
          WRITE(*,'(A,I6,A,I10)') '  CP ',CP_IDX,' -> FEM node ',
     .      CP_TO_NODE_MAP(CP_IDX)
          ENDDO
        WRITE(*,*)
        WRITE(*,'(A)') '=============================================='
        WRITE(*,*)
      ENDIF

C=======================================================================
C   Step 5: Note: HEX8 top nodes (IXS indices 6-9) will be looked up
C   from IXS in the mass routine using KQ1NP_TAB(10,*) which stores
C   the original HEX8 index. No need to store them here.
C=======================================================================

C=======================================================================
C   Step 5 Update NUMNOD
C=======================================================================
      NUMNOD_OUT = NUMNOD_IN + NUMNOD_CP_ADDED

C=======================================================================
C   Clean:up and output
C=======================================================================
      DEALLOCATE(CP_TO_NODE_MAP)
      
      IF (NUMNOD_CP_ADDED > 0) THEN
        WRITE(IOUT,'(A,I10,A)') 
     .    ' Q1Np: Promoted ', NUMNOD_CP_ADDED, 
     .    ' control points to mesh nodes'
        WRITE(IOUT,'(A,I10)') 
     .    ' New total node count: ', NUMNOD_OUT
      ENDIF

      RETURN
      END

