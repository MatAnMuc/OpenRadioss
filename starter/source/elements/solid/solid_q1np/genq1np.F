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
Chd|  GENQ1NP                       source/elements/solid/solid_q1np/genq1np.F
Chd|====================================================================
C=======================================================================
C   Generate Q1Np elements from surface definition
C   
C   This routine generates Q1Np (enriched hexahedral) elements internally
C   from an existing surface definition
C
C   Algorithm:
C   1. Extract surface segments from hardcoded surface ID (temporary solution)
C   2. Determine grid layout (NX x NY) from number of segments
C   3. Map surface segments to underlying HEX8 elements
C   4. Extract bulk nodes from HEX8 elements
C   5. Generate shared NURBS control points via procedural grid
C   6. Set up NURBS knot vectors and weights (open uniform)
C   7. Create Q1Np element connectivity arrays (CP + FEM bulk nodes separated)
C
C   Note: Currently uses hardcoded surface ID=1 and NURBS degrees p=q=2
C         for testing purposes. Control points are mapped directly from
C         surface node coordinates.
C=======================================================================
      SUBROUTINE GENQ1NP(IGRSURF,IXS,X,ITAB,IPART,
     .                   IPARTS,
     .                   KQ1NP_TAB,IQ1NP_TAB,IQ1NP_BULK_TAB,
     .                   Q1NP_WTAB,Q1NP_KTAB,Q1NP_CPTAB,
     .                   NWEIGHT_MAX,NUMELQ1NP_OUT)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
      USE GROUPDEF_MOD
      USE Q1NP_RESTART_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "param_c.inc"
#include      "scr17_c.inc"
#include      "units_c.inc"
#include      "elements.inc" 
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     IGRSURF    - Surface definitions array (input)
C     IXS        - HEX8 element connectivity array (input)
C     X          - Node coordinates array (input)
C     ITAB       - User to system node ID mapping (input)
C     IPART      - Part definitions array (input)
C     KQ1NP_TAB  - Q1Np element properties array (output)
C                  KQ1NP_TAB(1,*) = Material ID
C                  KQ1NP_TAB(2,*) = Property ID
C                  KQ1NP_TAB(3,*) = Number of control points
C                  KQ1NP_TAB(4,*) = Stores the starting index in IQ1NP_TAB for element i�s control-point connectivity
C                  KQ1NP_TAB(5,*) = Element ID
C                  KQ1NP_TAB(6,*) = Element index in u direction (0-based)
C                  KQ1NP_TAB(7,*) = Element index in v direction (0-based)
C                  KQ1NP_TAB(8,*) = NURBS degree p
C                  KQ1NP_TAB(9,*) = NURBS degree q
C                  KQ1NP_TAB(14,*)= Offset to bulk nodes in IQ1NP_BULK_TAB
C     IQ1NP_TAB  - Control point connectivity array (output)
C                  Stores (p+1)*(q+1) control-point IDs per element
C     IQ1NP_BULK_TAB - Bulk node connectivity array (output)
C                  Stores the 4 HEX8 bulk nodes per element
C     Q1NP_WTAB  - NURBS weights array (output, all 1.0 for non-rational)
C     Q1NP_KTAB  - NURBS knot vectors array (output)
C                  First NKNOT_U entries: U knot vector
C                  Next NKNOT_V entries: V knot vector
C     NUMELQ1NP_OUT  - Number of generated Q1Np elements (output)
C-----------------------------------------------
      INTEGER, PARAMETER :: IDEBUG_Q1NP = 1
C-----------------------------------------------
      TYPE (SURF_)   , DIMENSION(NSURF)   :: IGRSURF
      INTEGER IXS(NIXS,*)
      my_real X(3,*)
      INTEGER ITAB(*)
      INTEGER IPART(LIPART1,*)
      INTEGER KQ1NP_TAB(15,*),IQ1NP_TAB(*),IQ1NP_BULK_TAB(*)
      my_real Q1NP_WTAB(*),Q1NP_KTAB(*),Q1NP_CPTAB(3,*)
      INTEGER NWEIGHT_MAX
      INTEGER NUMELQ1NP_OUT
      INTEGER, INTENT(IN) :: IPARTS(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C     Hardcoded parameters for testing
      INTEGER ISURF_ID
      PARAMETER (ISURF_ID=1)      ! Surface ID to use for generation
      INTEGER P,Q
      PARAMETER (P=2,Q=2)         ! NURBS degrees (quadratic in u and v)
      INTEGER NBULKQ1NP
      PARAMETER (NBULKQ1NP=4)     ! Number of bulk nodes per Q1Np element
      INTEGER NKQ1NP
      PARAMETER (NKQ1NP=15)       ! Number of fields in KQ1NP_TAB per element

      INTEGER ISEG,IEL,IEL_HEX8
      INTEGER NODES_SURF(4)       ! Surface segment node IDs
      INTEGER NODES_BULK(4)       ! Bulk face node IDs from HEX8
      INTEGER NSEG                ! Number of surface segments
      INTEGER NX,NY               ! Grid dimensions (elements in u,v directions)
      INTEGER NCP_U,NCP_V         ! Number of control points in u,v directions
      INTEGER NODE_ID             ! Node ID for debug output
      INTEGER CP_NODE             ! Control point node ID
      INTEGER CP_IDX,CP_IDY       ! Control point grid indices
      INTEGER NCTRL               ! Control points per element: (p+1)*(q+1)
      INTEGER I,J,II,JJ           ! Loop indices
      INTEGER IEL_Q1NP            ! Q1Np element counter = 0
      INTEGER OFFSET_CTRL         ! Offset in IQ1NP_TAB for control points = 1
      INTEGER OFFSET_BULK         ! Offset in IQ1NP_BULK_TAB for bulk nodes = 1
      INTEGER OFFSET_KTAB         ! Offset in Q1NP_KTAB (unused, kept for future)
      INTEGER ELEM_U,ELEM_V       ! Element indices in u,v (unused, kept for future)
      INTEGER MID,PID             ! Material and property IDs = passed from underlying HEX8
      INTEGER ELEM_ID             ! Re-used element ID from HEX8
      INTEGER IEL_ORIG             ! Original HEX8 element local index

C     Control point map: maps (i,j) grid position to node ID
C     TODO: Fit NURBS surface to actual geometry for proper control points
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: CP_MAP
      INTEGER MAX_CP_U,MAX_CP_V   ! Dimensions of control point map
      my_real XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      my_real DX,DY,ZTOP,BIGVAL
      my_real, DIMENSION(:), ALLOCATABLE :: GRID_U,GRID_V
      INTEGER CP_COUNTER
      LOGICAL HAS_SURF_NODES
C=======================================================================
C   Step 1: Validate surface ID and get number of segments
C=======================================================================
C     Debug: print surface information
C=======================================================================
      IF (IDEBUG_Q1NP > 0) THEN
        PRINT*, 'Q1NP DEBUG: NSURF = ', NSURF
        PRINT*, 'Q1NP DEBUG: ISURF_ID = ', ISURF_ID
      ENDIF
      
      IF (ISURF_ID < 1 .OR. ISURF_ID > NSURF) THEN
        IF (IDEBUG_Q1NP > 0) THEN
          PRINT*, 'Q1NP DEBUG: Invalid surface ID'
        ENDIF
        CALL ANCMSG(MSGID=402,
     .             MSGTYPE=MSGERROR,
     .             ANMODE=ANINFO,
     .             C1="Q1NP",
     .             I1=ISURF_ID,
     .             I2=NSURF)
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF

      NSEG = IGRSURF(ISURF_ID)%NSEG
      IF (IDEBUG_Q1NP > 0) THEN
        PRINT*, 'Q1NP DEBUG: NSEG = ', NSEG
        PRINT*, 'Q1NP DEBUG: NSEG_IGE = ', IGRSURF(ISURF_ID)%NSEG_IGE
        PRINT*, 'Q1NP DEBUG: Surface ID = ', IGRSURF(ISURF_ID)%ID
      ENDIF
      
C     IMPORTANT: Surface segments might not be read yet when GENQ1NP is called
C     The surface exists (NSURF=1) but NSEG=0 means segments aren't populated
C     This suggests GENQ1NP is called too early in the starter pipeline
C     
C     For now, check both NSEG and NSEG_IGE fields
      IF (NSEG <= 0 .AND. IGRSURF(ISURF_ID)%NSEG_IGE <= 0) THEN
        IF (IDEBUG_Q1NP > 0) THEN
          PRINT*, 'Q1NP DEBUG: ERROR - No surface segments found!'
          PRINT*, 'Q1NP DEBUG: NSEG = ', NSEG, ', NSEG_IGE = ', IGRSURF(ISURF_ID)%NSEG_IGE
          PRINT*, 'Q1NP DEBUG: Surface segments may not be read yet.'
          PRINT*, 'Q1NP DEBUG: GENQ1NP might be called too early in the pipeline.'
          PRINT*, 'Q1NP DEBUG: Check if surfaces are read before GENQ1NP is called.'
        ENDIF
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF
      
C     Use NSEG if available, otherwise try NSEG_IGE
      IF (NSEG <= 0) THEN
        NSEG = IGRSURF(ISURF_ID)%NSEG_IGE
        IF (IDEBUG_Q1NP > 0) THEN
          PRINT*, 'Q1NP DEBUG: Using NSEG_IGE = ', NSEG
        ENDIF
      ENDIF

C=======================================================================
C   Step 2: Determine grid dimensions from surface segments
C=======================================================================
C     For now, assume square grid: nx = ny = sqrt(nseg)
C     This works for regular rectangular grids
C     TODO: Extract actual grid layout from segment arrangement
C           by analyzing node connectivity patterns
      NX = INT(SQRT(REAL(NSEG)))
      NY = NX
      IF (NX*NY < NSEG) NX = NX + 1
      IF (NX*NY < NSEG) NY = NY + 1

C=======================================================================
C   Step 3: Calculate NURBS parameters
C=======================================================================
C     Number of control points: ncp = ne + p
C     For ne elements and degree p, we need ne+p control points
      NCP_U = NX + P
      NCP_V = NY + Q
C     Number of control points per element (tensor product)
      NCTRL = (P+1)*(Q+1)

C=======================================================================
C   Step 4: Allocate and initialize control point map
C=======================================================================
C     Control point map stores node IDs for shared control points
C     Dimensions: (NCP_U x NCP_V) for the entire grid
      MAX_CP_U = NCP_U
      MAX_CP_V = NCP_V
      ALLOCATE(CP_MAP(MAX_CP_U,MAX_CP_V),STAT=IEL)
      IF (IEL /= 0) THEN
        CALL ANCMSG(MSGID=268,
     .             ANMODE=ANINFO,
     .             MSGTYPE=MSGERROR,
     .             C1='CP_MAP')
        RETURN
      ENDIF

C     Initialize control point map to zero
      DO J=1,MAX_CP_V
        DO I=1,MAX_CP_U
          CP_MAP(I,J) = 0
        ENDDO
      ENDDO
C=======================================================================
C   Step 5: Generate procedural control points
C=======================================================================
C     Control points are shared across elements in a tensor product layout.
C     Logic:
C       - Determine the planar bounding box from the surface nodes
C       - Build evenly spaced parametric coordinates via linspace
C       - Assign a unique CP ID per (u,v) location and store coordinates
C     Q1NP_CPTAB are the coordinates
C     CP_MAP is the control point map
C     TODO: Later this will be replace by best-fist
      BIGVAL = HUGE(ONE)
      XMIN =  BIGVAL
      YMIN =  BIGVAL
      ZMIN =  BIGVAL
      XMAX = -BIGVAL
      YMAX = -BIGVAL
      ZMAX = -BIGVAL
      HAS_SURF_NODES = .FALSE.

C     Bounding Box
      DO ISEG=1,NSEG
        DO II=1,4
          NODE_ID = IGRSURF(ISURF_ID)%NODES(ISEG,II)
          IF (NODE_ID > 0) THEN
            HAS_SURF_NODES = .TRUE.
            XMIN = MIN(XMIN, X(1,NODE_ID))
            XMAX = MAX(XMAX, X(1,NODE_ID))
            YMIN = MIN(YMIN, X(2,NODE_ID))
            YMAX = MAX(YMAX, X(2,NODE_ID))
            ZMIN = MIN(ZMIN, X(3,NODE_ID))
            ZMAX = MAX(ZMAX, X(3,NODE_ID))
          ENDIF
        ENDDO
      ENDDO
      IF (.NOT. HAS_SURF_NODES) THEN
C       Fallback: use unit square if no surface nodes were found
        XMIN = ZERO
        XMAX = ONE
        YMIN = ZERO
        YMAX = ONE
        ZMIN = ZERO
        ZMAX = ONE
      ENDIF

      ZTOP = ZMAX
      IF (ALLOCATED(GRID_U)) DEALLOCATE(GRID_U)
      IF (ALLOCATED(GRID_V)) DEALLOCATE(GRID_V)
      ALLOCATE(GRID_U(NCP_U))
      ALLOCATE(GRID_V(NCP_V))

C     Build evenly spaced parametric coordinates via linspace
      IF (NCP_U > 1) THEN
        DX = (XMAX - XMIN) / REAL(NCP_U-1)
      ELSE
        DX = ZERO
      ENDIF
      DO I=1,NCP_U
        GRID_U(I) = XMIN + REAL(I-1)*DX
      ENDDO

      IF (NCP_V > 1) THEN
        DY = (YMAX - YMIN) / REAL(NCP_V-1)
      ELSE
        DY = ZERO
      ENDIF
      DO J=1,NCP_V
        GRID_V(J) = YMIN + REAL(J-1)*DY
      ENDDO

      CP_COUNTER = 0
      DO J=1,NCP_V
        DO I=1,NCP_U
          CP_COUNTER = CP_COUNTER + 1
          CP_MAP(I,J) = CP_COUNTER
          Q1NP_CPTAB(1,CP_COUNTER) = GRID_U(I)
          Q1NP_CPTAB(2,CP_COUNTER) = GRID_V(J)
          Q1NP_CPTAB(3,CP_COUNTER) = ZTOP
        ENDDO
      ENDDO

C=======================================================================
C   Step 6: Debug: Print control point map
C=======================================================================
      IF (IDEBUG_Q1NP > 0) THEN
        PRINT*, ' '
        PRINT*, 'Q1NP DEBUG: Control Point Map (', MAX_CP_U, 'x', MAX_CP_V, '):'
        DO J=1,MAX_CP_V
          DO I=1,MAX_CP_U
            IF (CP_MAP(I,J) > 0) THEN
              CP_NODE = CP_MAP(I,J)
              WRITE(*,102) I,J,CP_NODE,
     .                     Q1NP_CPTAB(1,CP_NODE),
     .                     Q1NP_CPTAB(2,CP_NODE),
     .                     Q1NP_CPTAB(3,CP_NODE)
            ELSE
              PRINT*, '  CP_MAP(', I, ',', J, ') = 0'
            ENDIF
          ENDDO
        ENDDO
      ENDIF

C=======================================================================
C   Step 6: Set up NURBS knot vectors and weights
C=======================================================================
C     Generate open uniform knot vectors and initialize weights to 1.0
C     Knot vectors are stored in Q1NP_KTAB: U knots first, then V knots
      CALL SETUPNURBSQ1NP(NX,NY,P,Q,
     .                    Q1NP_KTAB,NCP_U,NCP_V,
     .                    Q1NP_WTAB,NWEIGHT_MAX)

C=======================================================================
C   Step 7: Generate Q1Np elements
C=======================================================================
C     Initialize offsets and counters
C     Control points connectivity and bulk nodes stored in separate arrays
      OFFSET_CTRL = 1
      OFFSET_BULK = 1
      IEL_Q1NP = 0

C     Loop over grid positions to create Q1Np elements
      DO J=1,NY
        DO I=1,NX
          ISEG = (J-1)*NX + I
          IF (ISEG > NSEG) GOTO 200

C         Get surface segment nodes (4 nodes forming a quad)
          DO II=1,4
            NODES_SURF(II) = IGRSURF(ISURF_ID)%NODES(ISEG,II)
          ENDDO

C         Find corresponding HEX8 element by matching any face; get bulk (opposite) face nodes
          CALL FINDHEX8FROMSURF(ISEG,NODES_SURF,IXS,
     .                          X,IEL_HEX8,NODES_BULK)

          IF (IEL_HEX8 <= 0) THEN
C           Error: no matching HEX8 element found
            CALL ANCMSG(MSGID=402,
     .                 MSGTYPE=MSGERROR,
     .                 ANMODE=ANINFO,
     .                 C1="Q1NP",
     .                 I1=ISEG,
     .                 I2=IEL_HEX8)
            GOTO 200
          ENDIF

C         Get material and property IDs from HEX8 element
C         Q1Np elements inherit these from the underlying HEX8
          MID = IXS(1,IEL_HEX8)
          PID = IXS(NIXS,IEL_HEX8)

C         Increment element counter and assign element ID
          IEL_Q1NP = IEL_Q1NP + 1
          IEL_ORIG = IEL_HEX8
          ELEM_ID = IXS(NIXS,IEL_HEX8)

C=======================================================================
C         Step 7a:   Store element properties in KQ1NP_TAB
C=======================================================================
          KQ1NP_TAB(1,IEL_Q1NP) = MID              ! Material ID
          KQ1NP_TAB(2,IEL_Q1NP) = PID              ! Property ID
          KQ1NP_TAB(3,IEL_Q1NP) = NCTRL            ! Number of control points
          KQ1NP_TAB(4,IEL_Q1NP) = OFFSET_CTRL      ! Offset to control points
          KQ1NP_TAB(5,IEL_Q1NP) = ELEM_ID          ! Element ID (original HEX)
          KQ1NP_TAB(6,IEL_Q1NP) = I-1              ! Element index u (0-based)
          KQ1NP_TAB(7,IEL_Q1NP) = J-1              ! Element index v (0-based)
          KQ1NP_TAB(8,IEL_Q1NP) = P                ! NURBS degree p
          KQ1NP_TAB(9,IEL_Q1NP) = Q                ! NURBS degree q
          KQ1NP_TAB(14,IEL_Q1NP) = OFFSET_BULK     ! Offset into IQ1NP_BULK_TAB
          KQ1NP_TAB(10,IEL_Q1NP) = IEL_ORIG        ! Local HEX8 index
          KQ1NP_TAB(11,IEL_Q1NP) = IPARTS(IEL_ORIG)! Owning part ID

C=======================================================================
C         Step 7b: Store control point connectivity (tensor product order)
C=======================================================================
C         Element at (i,j) uses control points: cp_map(i+ii, j+jj)
C         for ii in [0,p], jj in [0,q] (Fortran: ii=0..p, jj=0..q)
          DO JJ=0,Q
            DO II=0,P
              IF (I+II <= MAX_CP_U .AND. J+JJ <= MAX_CP_V) THEN
                IQ1NP_TAB(OFFSET_CTRL) = CP_MAP(I+II,J+JJ)
              ELSE
                IQ1NP_TAB(OFFSET_CTRL) = 0
              ENDIF
              OFFSET_CTRL = OFFSET_CTRL + 1
            ENDDO
          ENDDO

C=======================================================================
C         Step 7c: Store bulk nodes in IQ1NP_BULK_TAB in parametric order
C         (-1,-1), (+1,-1), (+1,+1), (-1,+1)
C=======================================================================
          IQ1NP_BULK_TAB(OFFSET_BULK)   = NODES_BULK(1)
          IQ1NP_BULK_TAB(OFFSET_BULK+1) = NODES_BULK(4)
          IQ1NP_BULK_TAB(OFFSET_BULK+2) = NODES_BULK(3)
          IQ1NP_BULK_TAB(OFFSET_BULK+3) = NODES_BULK(2)
          OFFSET_BULK = OFFSET_BULK + 4

C=======================================================================
C         Step 7d: Debug: Print element information (disabled – see Doku)
C=======================================================================
      IF (IDEBUG_Q1NP > 0) THEN
        PRINT*, ' '
        PRINT*, 'Q1NP DEBUG: '
        PRINT*, '  Element: ', IEL_Q1NP, ' (Global ID: ', ELEM_ID, ')'
        PRINT*, '  Grid position: (', I-1, ',', J-1, ')'
        PRINT*, '  Material ID:  ', MID
        PRINT*, '  Property ID:  ', PID
        !PRINT*, '  HEX8 element: ', IEL_HEX8, ' (Global ID: ', IXS(NIXS,IEL_HEX8), ')'
        PRINT*, '  Bulk nodes:'
        DO II=1,4
          NODE_ID = NODES_BULK(II)
          PRINT '(A,I7,A,F6.2,A,F6.2,A,F6.2,A)', '    Node ', NODE_ID, ': (', 
     .          X(1,NODE_ID), ', ', X(2,NODE_ID), ', ', X(3,NODE_ID), ')'
        ENDDO
        PRINT*, '  Control points (NURBS-defined surface):'
        DO JJ=0,Q
          DO II=0,P
            CP_IDX = I + II
            CP_IDY = J + JJ
            IF (CP_IDX <= MAX_CP_U .AND. CP_IDY <= MAX_CP_V) THEN
              CP_NODE = CP_MAP(CP_IDX,CP_IDY)
              PRINT 101, '    CP[', II, ',', JJ, '] = CP ', 
     .                  CP_NODE, ': (', Q1NP_CPTAB(1,CP_NODE), ', ', 
     .                  Q1NP_CPTAB(2,CP_NODE), ', ', 
     .                  Q1NP_CPTAB(3,CP_NODE), ')'
            ELSE
              PRINT*, '    CP[', II, ',', JJ, '] = 0 (not set)'
            ENDIF
          ENDDO
        ENDDO
      ENDIF

        ENDDO   ! end DO I=1,NX
      ENDDO     ! end DO J=1,NY

  200 CONTINUE

C     Set output: number of generated elements
      NUMELQ1NP_OUT = IEL_Q1NP

C     Inform Q1NP_RESTART_MOD about total control-point buffer length
      CALL SET_Q1NP_TABVINT_LEN(OFFSET_CTRL-1)

C=======================================================================
C   Print Q1Np element information
C=======================================================================
C     Always print Q1Np section header for visibility
      WRITE(IOUT,300)
      
      IF (NUMELQ1NP_OUT > 0) THEN
        WRITE(IOUT,301) NUMELQ1NP_OUT, NX, NY, NCP_U, NCP_V
        WRITE(IOUT,302)
        DO I=1,MIN(50,NUMELQ1NP_OUT)
          WRITE(IOUT,303) I, KQ1NP_TAB(5,I), KQ1NP_TAB(1,I),
     .                    KQ1NP_TAB(2,I), KQ1NP_TAB(6,I), KQ1NP_TAB(7,I)
        ENDDO
        IF (NUMELQ1NP_OUT > 50) THEN
          WRITE(IOUT,304) NUMELQ1NP_OUT - 50
        ENDIF
      ELSE
        WRITE(IOUT,305)
      ENDIF

 100  FORMAT(A,I0,A,F10.2,A,F10.2,A,F10.2,A)
 101  FORMAT(A,I0,A,I0,A,I0,A,F10.2,A,F10.2,A,F10.2,A)
 102  FORMAT('  CP_MAP(',I2,',',I2,') = CP ',I6,
     .       ': (',F8.2,', ',F8.2,', ',F8.2,')')
 300  FORMAT(/' Q1NP ENRICHED ELEMENTS'/
     +        ' ----------------------'/
     +        ' Generated from surface definition')
 301  FORMAT(' Number of Q1Np elements: ',I10/
     +       ' Grid dimensions: ',I3,' x ',I3/
     +       ' Control points: ',I3,' x ',I3)
 302  FORMAT('    LOC-EL     GLO-EL      MATER       GEOM    ELEM-U    ELEM-V'/
     +       '    NODES LIST')
 303  FORMAT(6(I10,1X))
 304  FORMAT(' ... and ',I10,' more elements')
 305  FORMAT(/' Q1NP ENRICHED ELEMENTS'/
     +        ' ----------------------'/
     +        ' No Q1Np elements generated')

C=======================================================================
C   Cleanup
C=======================================================================
      IF (ALLOCATED(CP_MAP)) DEALLOCATE(CP_MAP)
      IF (ALLOCATED(GRID_U)) DEALLOCATE(GRID_U)
      IF (ALLOCATED(GRID_V)) DEALLOCATE(GRID_V)

      RETURN
      END

