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
Chd|  Q1NP_REALLOCATE_NODE_ARRAYS    source/elements/solid/solid_q1np/q1np_reallocate_arrays.F
Chd|====================================================================
C=======================================================================
C   Reallocate node arrays to accommodate new control point nodes
C
C   This routine reallocates node arrays (X, V, D, MS, etc.) to a larger
C   size, preserving existing data. The reallocation is handled automatically
C   using MOVE_ALLOC, which transfers the allocation from temporary arrays
C   to the actual arrays.
C
C   Process:
C   1. Copy existing data from original arrays to temporary arrays (new size)
C   2. Initialize new indices (SX_OLD+1 to SX_NEW, etc.) to zero
C   3. Use MOVE_ALLOC to transfer temporary arrays to actual arrays
C      (this automatically reallocates the actual arrays to the new size)
C
C=======================================================================
      SUBROUTINE Q1NP_REALLOCATE_NODE_ARRAYS(
     .                   X, V, D, MS, IN, VR, DR, XYZREF,
     .                   SX_OLD, SV_OLD, SD_OLD, SMS_OLD, SIN_OLD, SVR_OLD, SDR_OLD,
     .                   SX_NEW, SV_NEW, SD_NEW, SMS_NEW, SIN_NEW, SVR_NEW, SDR_NEW,
     .                   NUMNOD_OLD, NUMNOD_NEW)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: SX_OLD, SV_OLD, SD_OLD, SMS_OLD, SIN_OLD, SVR_OLD, SDR_OLD
      INTEGER, INTENT(IN) :: SX_NEW, SV_NEW, SD_NEW, SMS_NEW, SIN_NEW, SVR_NEW, SDR_NEW
      INTEGER, INTENT(IN) :: NUMNOD_OLD, NUMNOD_NEW
      my_real, INTENT(INOUT), ALLOCATABLE :: X(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: V(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: D(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: MS(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: IN(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: VR(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: DR(:)
      my_real, INTENT(INOUT), ALLOCATABLE :: XYZREF(:)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      my_real, ALLOCATABLE :: X_TMP(:)
      my_real, ALLOCATABLE :: V_TMP(:)
      my_real, ALLOCATABLE :: D_TMP(:)
      my_real, ALLOCATABLE :: MS_TMP(:)
      my_real, ALLOCATABLE :: IN_TMP(:)
      my_real, ALLOCATABLE :: VR_TMP(:)
      my_real, ALLOCATABLE :: DR_TMP(:)
      my_real, ALLOCATABLE :: XYZREF_TMP(:)
      INTEGER I, STAT
      INTEGER, PARAMETER :: IDEBUG_Q1NP = 2
C=======================================================================

C   Allocate temporary arrays with new sizes
      ALLOCATE(X_TMP(SX_NEW), STAT=STAT)
      IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                          MSGTYPE=MSGERROR, C1='X_TMP')
      ALLOCATE(V_TMP(SV_NEW), STAT=STAT)
      IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                          MSGTYPE=MSGERROR, C1='V_TMP')
      IF(SD_NEW > 0) THEN
        ALLOCATE(D_TMP(SD_NEW), STAT=STAT)
        IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                            MSGTYPE=MSGERROR, C1='D_TMP')
      ENDIF
      ALLOCATE(MS_TMP(SMS_NEW), STAT=STAT)
      IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                          MSGTYPE=MSGERROR, C1='MS_TMP')
      IF(SIN_NEW > 0) THEN
        ALLOCATE(IN_TMP(SIN_NEW), STAT=STAT)
        IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                            MSGTYPE=MSGERROR, C1='IN_TMP')
      ENDIF
      IF(SVR_NEW > 0) THEN
        ALLOCATE(VR_TMP(SVR_NEW), STAT=STAT)
        IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                            MSGTYPE=MSGERROR, C1='VR_TMP')
      ENDIF
      IF(SDR_NEW > 0) THEN
        ALLOCATE(DR_TMP(SDR_NEW), STAT=STAT)
        IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                            MSGTYPE=MSGERROR, C1='DR_TMP')
      ENDIF
      ALLOCATE(XYZREF_TMP(SX_NEW), STAT=STAT)
      IF(STAT /= 0) CALL ANCMSG(MSGID=268,ANMODE=ANINFO,
     .                          MSGTYPE=MSGERROR, C1='XYZREF_TMP')

C=======================================================================
C   Step 1: Initialize new arrays to zero
C=======================================================================
      X_TMP = ZERO
      V_TMP = ZERO
      IF(SD_NEW > 0) D_TMP = ZERO
      MS_TMP = ZERO
      IF(SIN_NEW > 0) IN_TMP = ZERO
      IF(SVR_NEW > 0) VR_TMP = ZERO
      IF(SDR_NEW > 0) DR_TMP = ZERO
      XYZREF_TMP = ZERO

C=======================================================================
C   Step 2: Copy existing data from original arrays to temporary arrays
C=======================================================================
      IF(SX_OLD > 0) THEN
C       Copy data using array section assignment
        X_TMP(1:MIN(SX_OLD, SX_NEW)) = X(1:MIN(SX_OLD, SX_NEW))
      ENDIF
      IF(SV_OLD > 0) THEN
        V_TMP(1:MIN(SV_OLD, SV_NEW)) = V(1:MIN(SV_OLD, SV_NEW))
      ENDIF
      IF(SMS_OLD > 0) THEN
        MS_TMP(1:MIN(SMS_OLD, SMS_NEW)) = MS(1:MIN(SMS_OLD, SMS_NEW))
      ENDIF
      IF(SD_OLD > 0 .AND. SD_NEW > 0) THEN
        D_TMP(1:MIN(SD_OLD, SD_NEW)) = D(1:MIN(SD_OLD, SD_NEW))
      ENDIF
      IF(SIN_OLD > 0 .AND. SIN_NEW > 0) THEN
        IN_TMP(1:MIN(SIN_OLD, SIN_NEW)) = IN(1:MIN(SIN_OLD, SIN_NEW))
      ENDIF
      IF(SVR_OLD > 0 .AND. SVR_NEW > 0) THEN
        VR_TMP(1:MIN(SVR_OLD, SVR_NEW)) = VR(1:MIN(SVR_OLD, SVR_NEW))
      ENDIF
      IF(SDR_OLD > 0 .AND. SDR_NEW > 0) THEN
        DR_TMP(1:MIN(SDR_OLD, SDR_NEW)) = DR(1:MIN(SDR_OLD, SDR_NEW))
      ENDIF
      IF(SX_OLD > 0) THEN
        XYZREF_TMP(1:MIN(SX_OLD, SX_NEW)) = XYZREF(1:MIN(SX_OLD, SX_NEW))
      ENDIF
      
C=======================================================================
C   Step 3: Deallocate original arrays and move temporary arrays to actual arrays
C=======================================================================
      IF(ALLOCATED(X)) DEALLOCATE(X)
      CALL MOVE_ALLOC(X_TMP, X)
      
      IF(ALLOCATED(V)) DEALLOCATE(V)
      CALL MOVE_ALLOC(V_TMP, V)
      
      IF(SD_NEW > 0) THEN
        IF(ALLOCATED(D)) DEALLOCATE(D)
        CALL MOVE_ALLOC(D_TMP, D)
      ENDIF
      
      IF(ALLOCATED(MS)) DEALLOCATE(MS)
      CALL MOVE_ALLOC(MS_TMP, MS)
      
      IF(SIN_NEW > 0) THEN
        IF(ALLOCATED(IN)) DEALLOCATE(IN)
        CALL MOVE_ALLOC(IN_TMP, IN)
      ENDIF
      
      IF(SVR_NEW > 0) THEN
        IF(ALLOCATED(VR)) DEALLOCATE(VR)
        CALL MOVE_ALLOC(VR_TMP, VR)
      ENDIF
      
      IF(SDR_NEW > 0) THEN
        IF(ALLOCATED(DR)) DEALLOCATE(DR)
        CALL MOVE_ALLOC(DR_TMP, DR)
      ENDIF
      
      IF(ALLOCATED(XYZREF)) DEALLOCATE(XYZREF)
      CALL MOVE_ALLOC(XYZREF_TMP, XYZREF)
      
C=======================================================================
C   Note: Q1NP_REALLOCATE_NODE_ARRAYS is no longer used. Reallocation
C   logic has been moved directly into lectur.F. This subroutine is
C   kept only as documentation of the previous approach.
C=======================================================================  
      RETURN
      END SUBROUTINE Q1NP_REALLOCATE_NODE_ARRAYS

