!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
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
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!Chd|====================================================================
!Chd|  GENQ1NP                       source/elements/solid/solid_q1np/genq1np.F90
!Chd|====================================================================
!C=======================================================================
!C
!C   This routine generates Q1Np (enriched hexahedral) elements
!C   from a surface-based definition.
!C
!C   Algorithm:
!C   1. Validate the requested surface ID and count its segments.
!C   2. Build a structured surface grid (NX x NY) from segment connectivity.
!C   3. Map surface segments to underlying HEX8 elements.
!C   4. Extract bulk face nodes from the associated HEX8 elements.
!C   5. Set up NURBS knot vectors and weights and fit NURBS control points
!C      to the surface via least-squares / Tikhonov (depending on method).
!C   6. Generate Q1Np element connectivity (NURBS control points + bulk nodes).
!C
!C   Note: Currently uses a fixed surface ID (ISURF_ID=1) and NURBS degrees
!C         P=Q=2. These can be generalized in future extensions.
!C=======================================================================
      module genq1np_mod
        use message_mod
        use groupdef_mod
        use q1np_restart_mod
        use q1np_volume_mod
        use q1np_cholesky_mod, only : cholesky_solve_q1np
        use q1np_export_csv_mod
        use q1np_surf_grid_mod
        use setupnurbsq1np_mod
        use findhex8fromsurface_mod, only : findhex8fromsurf
        use precision_mod, only : WP
        use constant_mod, only : ZERO, ONE, HALF, TWO
        implicit none
!C     Control of optional CSV debug export (kept ON by default)
        logical :: q1np_export_csv = .TRUE.
      contains
        subroutine genq1np(igrsurf, ixs, x, iparts,&
     &      nsurf, nixs, lipart1, numnod, numels, iout,&
     &      kq1np_tab, iq1np_tab, iq1np_bulk_tab,&
     &      q1np_wtab, q1np_ktab, q1np_cptab, nweight_max, numelq1np_out)
!C-----------------------------------------------
!C   M o d u l e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   I m p l i c i t   T y p e s
!C-----------------------------------------------
!C-----------------------------------------------
!C   C o m m o n   B l o c k s
!C-----------------------------------------------
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
!C     IGRSURF    - Surface definitions array (input)
!C     IXS        - HEX8 element connectivity array (input)
!C     X          - Node coordinates array (input)
!C     KQ1NP_TAB  - Q1Np element properties array (output)
!C                  KQ1NP_TAB(1,*) = Material ID
!C                  KQ1NP_TAB(2,*) = Property ID
!C                  KQ1NP_TAB(3,*) = Number of control points
!C                  KQ1NP_TAB(4,*) = Starting index in IQ1NP_TAB for element control-point connectivity
!C                  KQ1NP_TAB(5,*) = Element ID
!C                  KQ1NP_TAB(6,*) = Element index in u direction (0-based)
!C                  KQ1NP_TAB(7,*) = Element index in v direction (0-based)
!C                  KQ1NP_TAB(8,*) = NURBS degree p
!C                  KQ1NP_TAB(9,*) = NURBS degree q
!C                  KQ1NP_TAB(14,*)= Offset to bulk nodes in IQ1NP_BULK_TAB
!C     IQ1NP_TAB  - Control point connectivity array (output)
!C                  Stores (p+1)*(q+1) control-point IDs per element
!C     IQ1NP_BULK_TAB - Bulk node connectivity array (output)
!C                  Stores the 4 HEX8 bulk nodes per element
!C     Q1NP_WTAB  - NURBS weights array (output, all 1.0 for non-rational)
!C     Q1NP_KTAB  - NURBS knot vectors array (output)
!C                  First NKNOT_U entries: U knot vector
!C                  Next NKNOT_V entries: V knot vector
!C     NUMELQ1NP_OUT  - Number of generated Q1Np elements (output)
!C-----------------------------------------------
      INTEGER, PARAMETER :: IDEBUG_Q1NP = 0
!C     CP fitting method:
!C       0 = 3D least-squares fitting with subdivision samples
!C       1 = Tikhonov-regularized least-squares using grid nodes only
      INTEGER, PARAMETER :: IQ1NP_CP_METHOD = 0
      INTEGER, PARAMETER :: DIV = 2 ! subdivision per quad for method 0 (fitting samples)
!C-----------------------------------------------
          type(surf_), dimension(nsurf), intent(in) :: igrsurf
          integer, intent(in) :: ixs(nixs, *)
          real(kind=WP), intent(in) :: x(3, *)
          integer, intent(in) :: nsurf, nixs, lipart1, numnod, numels, iout
          integer, intent(in) :: nweight_max
          integer, intent(out) :: numelq1np_out
          integer, intent(in) :: iparts(*)
          integer, intent(inout) :: kq1np_tab(:,:), iq1np_tab(:), iq1np_bulk_tab(:)
          real(kind=WP), intent(inout) :: q1np_wtab(:), q1np_ktab(:), q1np_cptab(:,:)
!-------------------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!C-----------------------------------------------
!C     Hardcoded parameters for testing
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
      INTEGER I,J,K,II,JJ         ! Loop indices
      INTEGER IEL_Q1NP            ! Q1Np element counter = 0
      INTEGER OFFSET_CTRL         ! Offset in IQ1NP_TAB for control points = 1
      INTEGER OFFSET_BULK         ! Offset in IQ1NP_BULK_TAB for bulk nodes = 1
      INTEGER MID,PID             ! Material and property IDs = passed from underlying HEX8
      INTEGER ELEM_ID             ! Re-used element ID from HEX8
      INTEGER IEL_ORIG            ! Original HEX8 element local index

!C     Control point map: maps (i,j) control-point grid position to linear CP index
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: CP_MAP
      INTEGER MAX_CP_U,MAX_CP_V   ! Dimensions of control point map
      INTEGER CP_COUNTER

!C     Method least-squares fit: surface grid, data points, normal equations
      INTEGER NDATA,NCP,NKNOT_U,NKNOT_V

!C     Surface grid from connectivity (segment -> (I,J), grid corner -> node)
      INTEGER, ALLOCATABLE :: SEG_I(:), SEG_J(:)
      INTEGER, ALLOCATABLE :: GRID_NODE(:,:), GRID_TO_SEG(:,:)
      INTEGER IERR_GRID

!C=======================================================================
!C   Step 1: Validate surface ID and get number of segments
!C=======================================================================
!C     Debug: print surface information
!C=======================================================================
      IF (IDEBUG_Q1NP > 0) THEN
        PRINT*, 'Q1NP DEBUG: Number of surfaces NSURF = ', NSURF
        PRINT*, 'Q1NP DEBUG: Surface ID ISURF_ID = ', ISURF_ID
      ENDIF

      IF (ISURF_ID < 1 .OR. ISURF_ID > NSURF) THEN
        CALL ANCMSG(MSGID=402, &
     &      MSGTYPE=MSGERROR, &
     &      ANMODE=ANINFO, &
     &      C1="Q1NP", &
     &      I1=ISURF_ID, &
     &      I2=NSURF)
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF

      NSEG = IGRSURF(ISURF_ID)%NSEG

!C
!C     For now, check both NSEG and NSEG_IGE fields
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

!C=======================================================================
!C   Step 2: Determine grid dimensions and segment/grid mapping from connectivity
!C=======================================================================
!C     Allocate work arrays for surface grid: segment indices (SEG_I, SEG_J),
!C     grid-to-segment map (GRID_TO_SEG), and grid node connectivity (GRID_NODE).
      ALLOCATE(SEG_I(NSEG), SEG_J(NSEG), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='SEG_I')
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF
      ALLOCATE(GRID_TO_SEG(NSEG, NSEG), STAT=IEL)
      IF (IEL .NE. 0) THEN
        DEALLOCATE(SEG_I, SEG_J)
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='GRID_TO_SEG')
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF
!C     GRID_NODE: leading dimension NSEG*4 so (NX+1,NY+1) always fits (NX*NY=NSEG).
      ALLOCATE(GRID_NODE(NSEG*4, NSEG*4), STAT=IEL)
      IF (IEL .NE. 0) THEN
        DEALLOCATE(SEG_I, SEG_J, GRID_TO_SEG)
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='GRID_NODE')
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF

!C     Build surface grid from segment connectivity; fills GRID_NODE(1:NX+1,1:NY+1)
!C     and GRID_TO_SEG(1:NX,1:NY). NX, NY are the element counts in u and v.
      CALL Q1NP_BUILD_SURF_GRID(IGRSURF(ISURF_ID), NSEG, &
     &      NX, NY, SEG_I, SEG_J, &
     &      GRID_NODE, GRID_TO_SEG, IERR_GRID)

      IF (IERR_GRID .NE. 0) THEN
        IF (ALLOCATED(SEG_I)) DEALLOCATE(SEG_I)
        IF (ALLOCATED(SEG_J)) DEALLOCATE(SEG_J)
        IF (ALLOCATED(GRID_NODE)) DEALLOCATE(GRID_NODE)
        IF (ALLOCATED(GRID_TO_SEG)) DEALLOCATE(GRID_TO_SEG)
        CALL ANCMSG(MSGID=402, &
     &      MSGTYPE=MSGERROR, &
     &      ANMODE=ANINFO, &
     &      C1="Q1NP", &
     &      I1=ISURF_ID, &
     &      I2=IERR_GRID)
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF

!C=======================================================================
!C   Step 3: Calculate NURBS parameters
!C=======================================================================
!C     Number of control points: ncp = ne + p
!C     For ne elements and degree p, we need ne+p control points
      NCP_U = NX + P
      NCP_V = NY + Q
!C     Number of control points per element (tensor product)
      NCTRL = (P+1)*(Q+1)

!C=======================================================================
!C   Set up NURBS knot vectors and weights (baseline for fitting)
!C=======================================================================
      CALL SETUPNURBSQ1NP(NX,NY,P,Q, &
     &      Q1NP_KTAB,NCP_U,NCP_V, &
     &      Q1NP_WTAB,NWEIGHT_MAX)


      NKNOT_U = NX + 2*P + 1
      NKNOT_V = NY + 2*Q + 1

      ! Both U and V knot vectors are stored in Q1NP_KTAB
      ! NCP_U and NCP_V are the number of control points in u and v directions
      ! NCP is the total number of control points
      NCP = NCP_U * NCP_V

!C     Subdivision for sampling (method 0 only): DIV*DIV interior sample points
!C     per quad (ensures NDATA >= NCP for fitting). NDATA is the number of data
!C     points for the least-squares fitting. Note: DIV=1 already yields one
!C     interior sample at element center (no duplicates with grid nodes)

      IF (IQ1NP_CP_METHOD .EQ. 0) THEN
        NDATA = (NX+1)*(NY+1) + NX*NY*DIV*DIV
      ELSE
!C       Method 1 (Tikhonov) uses only surface grid nodes as data points
        NDATA = (NX+1)*(NY+1)
      ENDIF

!C=======================================================================
!C   Step 4: Allocate and initialize control point map
!C=======================================================================
!C     CP_MAP(i,j) = linear control-point index 1..NCP for grid (i,j).
!C     Used for element connectivity (which CPs belong to each Q1Np element).
      MAX_CP_U = NCP_U
      MAX_CP_V = NCP_V
      ALLOCATE(CP_MAP(MAX_CP_U,MAX_CP_V), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='CP_MAP')
        NUMELQ1NP_OUT = 0
        RETURN
      ENDIF
      CP_COUNTER = 0
      DO J=1,NCP_V
        DO I=1,NCP_U
          CP_COUNTER = CP_COUNTER + 1
          CP_MAP(I,J) = CP_COUNTER
        ENDDO
      ENDDO

!C=======================================================================
!C   Step 5: Least-squares fitting to derive the position of the NURBS control points
!C=======================================================================
      CALL Q1NP_FIT_CONTROL_POINTS( NX, NY, P, Q,           &
     &     NCP_U, NCP_V, NCP, NDATA, NKNOT_U, NKNOT_V,      &
     &     GRID_NODE, X, Q1NP_KTAB, CP_MAP, Q1NP_CPTAB,     &
     &     IDEBUG_Q1NP, IQ1NP_CP_METHOD, DIV )
!C=======================================================================
!C   Step 6: Generate Q1Np elements from fitted NURBS surface and HEX8 bulk mesh
!C=======================================================================
!C     Offsets into IQ1NP_TAB (control points) and IQ1NP_BULK_TAB (bulk nodes)
      OFFSET_CTRL = 1
      OFFSET_BULK = 1
      IEL_Q1NP = 0

!C     Loop over grid positions to create Q1Np elements
      elem_loop: do j = 1, ny
        do i = 1, nx
          ISEG = GRID_TO_SEG(I,J)
          IF (ISEG .LE. 0) exit elem_loop

!C         Get surface segment nodes (4 nodes forming a quad)
          DO II=1,4
            NODES_SURF(II) = IGRSURF(ISURF_ID)%NODES(ISEG,II)
          ENDDO

!C         Find corresponding HEX8 element by matching any face; get bulk (opposite) face nodes
          CALL findhex8fromsurf(NODES_SURF,IXS,IEL_HEX8,NODES_BULK,NIXS,NUMELS)

          IF (IEL_HEX8 <= 0) THEN
!C           Error: no matching HEX8 element found
            CALL ANCMSG(MSGID=402, &
     &      MSGTYPE=MSGERROR, &
     &      ANMODE=ANINFO, &
     &      C1="Q1NP", &
     &      I1=ISEG, &
     &                  I2=IEL_HEX8)
            exit elem_loop
          ENDIF

!C         Get material and property IDs from HEX8 element
!C         Q1Np elements inherit these from the underlying HEX8
          MID = IXS(1,IEL_HEX8)
          PID = IXS(NIXS,IEL_HEX8)

!C         Increment element counter and assign element ID
          IEL_Q1NP = IEL_Q1NP + 1
          IEL_ORIG = IEL_HEX8
          ELEM_ID = IXS(NIXS,IEL_HEX8)

!C=======================================================================
!C         Step 6a:   Store element properties in KQ1NP_TAB
!C=======================================================================
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

!C=======================================================================
!C         Step 6b: Store control point connectivity (tensor product order)
!C=======================================================================
!C         Element at (i,j) uses control points: cp_map(i+ii, j+jj)
!C         for ii in [0,p], jj in [0,q] (Fortran: ii=0..p, jj=0..q)
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

!C=======================================================================
!C         Step 6c: Store bulk nodes in IQ1NP_BULK_TAB in parametric order
!C         (-1,-1), (+1,-1), (+1,+1), (-1,+1)
!C=======================================================================
          IQ1NP_BULK_TAB(OFFSET_BULK)   = NODES_BULK(1)
          IQ1NP_BULK_TAB(OFFSET_BULK+1) = NODES_BULK(4)
          IQ1NP_BULK_TAB(OFFSET_BULK+2) = NODES_BULK(3)
          IQ1NP_BULK_TAB(OFFSET_BULK+3) = NODES_BULK(2)
          OFFSET_BULK = OFFSET_BULK + 4

!C=======================================================================
!C         Step 6d: Debug: Print element information
!C=======================================================================
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
          PRINT '(A,I7,A,F6.2,A,F6.2,A,F6.2,A)', '    Node ', NODE_ID, ': (', &
     &      X(1,NODE_ID), ', ', X(2,NODE_ID), ', ', X(3,NODE_ID), ')'
        ENDDO
        PRINT*, '  Control points (NURBS-defined surface):'
        DO JJ=0,Q
          DO II=0,P
            CP_IDX = I + II
            CP_IDY = J + JJ
            IF (CP_IDX <= MAX_CP_U .AND. CP_IDY <= MAX_CP_V) THEN
              CP_NODE = CP_MAP(CP_IDX,CP_IDY)
              PRINT 101, '    CP[', II, ',', JJ, '] = CP ', &
     &      CP_NODE, ': (', Q1NP_CPTAB(1,CP_NODE), ', ', &
     &      Q1NP_CPTAB(2,CP_NODE), ', ', &
     &      Q1NP_CPTAB(3,CP_NODE), ')'
            ELSE
              PRINT*, '    CP[', II, ',', JJ, '] = 0 (not set)'
            ENDIF
          ENDDO
        ENDDO
      ENDIF

        end do   ! i
      end do elem_loop   ! j

      NUMELQ1NP_OUT = IEL_Q1NP
!C=======================================================================
!C   Optional: export NURBS control points and bulk nodes to CSV (visualization)
!C=======================================================================
      IF (q1np_export_csv) THEN
        CALL Q1NP_EXPORT_NURBS_CSV(NCP_U,NCP_V,MAX_CP_U,MAX_CP_V, &
     &      Q1NP_CPTAB,CP_MAP, &
     &      NUMELQ1NP_OUT,KQ1NP_TAB, &
     &      IQ1NP_TAB,IQ1NP_BULK_TAB,NUMNOD)
        CALL Q1NP_EXPORT_BULK_NODES_CSV(NUMELQ1NP_OUT, &
     &      KQ1NP_TAB,IQ1NP_BULK_TAB,X)
      ENDIF

!C     Inform Q1NP_RESTART_MOD of total control-point buffer length (for restart I/O)
      CALL SET_Q1NP_TABVINT_LEN(OFFSET_CTRL-1)

!C=======================================================================
!C   Print Q1Np element information
!C=======================================================================
!C     Always print Q1Np section header for visibility
      WRITE(IOUT,300)

      IF (NUMELQ1NP_OUT > 0) THEN
        WRITE(IOUT,301) NUMELQ1NP_OUT, NX, NY, NCP_U, NCP_V
        WRITE(IOUT,302)
        DO I=1,MIN(50,NUMELQ1NP_OUT)
          WRITE(IOUT,303) I, KQ1NP_TAB(5,I), KQ1NP_TAB(1,I), &
     &      KQ1NP_TAB(2,I), KQ1NP_TAB(6,I), KQ1NP_TAB(7,I)
        ENDDO
        IF (NUMELQ1NP_OUT > 50) THEN
          WRITE(IOUT,304) NUMELQ1NP_OUT - 50
        ENDIF
      ELSE
        WRITE(IOUT,305)
      ENDIF

 100  FORMAT(A,I0,A,F10.2,A,F10.2,A,F10.2,A)
 101  FORMAT(A,I0,A,I0,A,I0,A,F10.2,A,F10.2,A,F10.2,A)
 300  FORMAT(/' Q1NP ENRICHED ELEMENTS'/ &
     &        ' ----------------------'/ &
     &        ' Generated from surface definition')
 301  FORMAT(' Number of Q1Np elements: ',I10/ &
     &       ' Grid dimensions: ',I3,' x ',I3/ &
     &       ' Control points: ',I3,' x ',I3)
 302  FORMAT('    LOC-EL     GLO-EL      MATER       GEOM    ELEM-U    ELEM-V'/ &
     &       '    NODES LIST')
 303  FORMAT(6(I10,1X))
 304  FORMAT(' ... and ',I10,' more elements')
 305  FORMAT(/' Q1NP ENRICHED ELEMENTS'/ &
     &        ' ----------------------'/ &
     &        ' No Q1Np elements generated')

!C=======================================================================
!C   Cleanup: deallocate grid and CP map work arrays
!C=======================================================================
      IF (ALLOCATED(CP_MAP)) DEALLOCATE(CP_MAP)
      IF (ALLOCATED(SEG_I)) DEALLOCATE(SEG_I)
      IF (ALLOCATED(SEG_J)) DEALLOCATE(SEG_J)
      IF (ALLOCATED(GRID_NODE)) DEALLOCATE(GRID_NODE)
      IF (ALLOCATED(GRID_TO_SEG)) DEALLOCATE(GRID_TO_SEG)
      RETURN
        end subroutine genq1np





!C=======================================================================
!C   Fit NURBS control points to the surface (least-squares / Tikhonov)
!C=======================================================================
           SUBROUTINE Q1NP_FIT_CONTROL_POINTS( NX, NY, P, Q,           &
     &     NCP_U, NCP_V, NCP, NDATA, NKNOT_U, NKNOT_V,            &
     &     GRID_NODE, X, Q1NP_KTAB, CP_MAP, Q1NP_CPTAB,           &
     &     IDEBUG_IN, ICP_METHOD_IN, DIV_IN )

        USE precision_mod,   ONLY : WP
        USE constant_mod,    ONLY : ZERO, ONE, HALF, TWO
        USE q1np_cholesky_mod, ONLY : cholesky_solve_q1np
        USE q1np_export_csv_mod
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: NX, NY, P, Q
        INTEGER, INTENT(IN) :: NCP_U, NCP_V, NCP, NDATA
        INTEGER, INTENT(IN) :: NKNOT_U, NKNOT_V
        INTEGER, INTENT(IN) :: GRID_NODE(:,:)
        REAL(WP), INTENT(IN) :: X(3,*)
        REAL(WP), INTENT(IN) :: Q1NP_KTAB(*)
        INTEGER, INTENT(IN) :: CP_MAP(:,:)
        INTEGER, INTENT(IN) :: IDEBUG_IN, ICP_METHOD_IN, DIV_IN
        REAL(WP), INTENT(INOUT) :: Q1NP_CPTAB(3,*)

        ! Local arrays
        REAL(WP), DIMENSION(:,:,:), ALLOCATABLE :: X_GRID
        INTEGER,  DIMENSION(:),      ALLOCATABLE :: GRID_NODE_RES
        REAL(WP), DIMENSION(:,:),    ALLOCATABLE :: DATA_PT
        REAL(WP), DIMENSION(:),      ALLOCATABLE :: U_PARAM, V_PARAM
        REAL(WP), DIMENSION(:,:),    ALLOCATABLE :: ATA, ATD, C_CP
        REAL(WP), DIMENSION(:),      ALLOCATABLE :: A_ROW
        REAL(WP), DIMENSION(:),      ALLOCATABLE :: U_KNOT, V_KNOT
        REAL(WP), DIMENSION(:,:),    ALLOCATABLE :: L_REG, LTL, ATA_REG
        REAL(WP), DIMENSION(:,:),    ALLOCATABLE :: C_CP_BEST
        REAL(WP), DIMENSION(:),      ALLOCATABLE :: TMPVEC

        ! Local scalars
        INTEGER :: IEL
        INTEGER :: I, J, K
        INTEGER :: IDATA, COL, COL1, COL2
        INTEGER :: II_SUB, JJ_SUB
        INTEGER :: NLAM, ILAM
        INTEGER :: NODE_ID, CP_COUNTER
        REAL(WP) :: UU, VV, S, T
        REAL(WP) :: DTD, RES2, RMS, BEST_RMS, LAMBDA_CUR, LAMBDA_BEST


!C     Allocate arrays for least-squares fit
      ! X_GRID: surface grid
      ! DATA_PT: data points
      ! ATA: left-hand side of normal equations
      ! ATD: right-hand side of normal equations
      ! A_ROW: basis functions at a given data point
      ! C_CP: control points
      ! U_KNOT: knot vector for u
      ! V_KNOT: knot vector for v

      ALLOCATE(X_GRID(3,NX+1,NY+1), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='X_GRID')
        RETURN
      ENDIF

      ALLOCATE(GRID_NODE_RES((NX+1)*(NY+1)), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='GRID_NODE_RES')
        DEALLOCATE(X_GRID)
        RETURN
      ENDIF

      ALLOCATE(DATA_PT(3,NDATA),U_PARAM(NDATA),V_PARAM(NDATA), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='DATA_PT/U_PARAM/V_PARAM')
        DEALLOCATE(X_GRID,GRID_NODE_RES)
        RETURN
      ENDIF

      ALLOCATE(ATA(NCP,NCP),ATD(NCP,3),A_ROW(NCP),C_CP(NCP,3), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='ATA/ATD/A_ROW/C_CP')
        DEALLOCATE(X_GRID,GRID_NODE_RES,DATA_PT,U_PARAM,V_PARAM)
        RETURN
      ENDIF

      ALLOCATE(U_KNOT(NKNOT_U),V_KNOT(NKNOT_V), STAT=IEL)
      IF (IEL .NE. 0) THEN
        CALL ANCMSG(MSGID=268,ANMODE=ANINFO,MSGTYPE=MSGERROR,C1='U_KNOT/V_KNOT')
        DEALLOCATE(X_GRID,GRID_NODE_RES,DATA_PT,U_PARAM,V_PARAM)
        DEALLOCATE(ATA,ATD,A_ROW,C_CP)
        RETURN
      ENDIF

!C     Fill X_GRID and GRID_NODE_RES from connectivity-derived grid (GRID_NODE)
      K=1
      DO J=1,NY+1
        DO I=1,NX+1
          NODE_ID = GRID_NODE(I,J)
          IF (NODE_ID .GT. 0) THEN
            X_GRID(1,I,J) = X(1,NODE_ID)
            X_GRID(2,I,J) = X(2,NODE_ID)
            X_GRID(3,I,J) = X(3,NODE_ID)
            GRID_NODE_RES(K) = NODE_ID
            K = K + 1
          ENDIF
        ENDDO
      ENDDO

!C=======================================================================
!C   Optional: export original surface nodes to CSV (for visualization/debug)
!C=======================================================================
      IF (q1np_export_csv) THEN
        CALL Q1NP_EXPORT_SURFACE_NODES_CSV(NX,NY,X_GRID,GRID_NODE_RES)
      ENDIF

!C     Fill DATA_PT with grid node coordinates; U_PARAM, V_PARAM = parametric (u,v) in [0,1]
      IDATA = 0
      DO J=1,NY+1
        DO I=1,NX+1
          IDATA = IDATA + 1
          DATA_PT(1,IDATA) = X_GRID(1,I,J)
          DATA_PT(2,IDATA) = X_GRID(2,I,J)
          DATA_PT(3,IDATA) = X_GRID(3,I,J)
          IF (NX > 0) U_PARAM(IDATA) = REAL(I-1)/REAL(NX)
          IF (NY > 0) V_PARAM(IDATA) = REAL(J-1)/REAL(NY)
        ENDDO
      ENDDO

!C     Extract knot vectors for basis function evaluation
      CALL Q1NP_GET_KNOT_VECTORS(NX, NY, P, Q, Q1NP_KTAB, U_KNOT, V_KNOT)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C   METHOD 0: 3D plain least-squares (subdivision points)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (ICP_METHOD_IN  .EQ. 0) THEN

!C     ------------------------------------------------------------------
!C     For method 0, add subdivision samples (DIV x DIV per quad) to
!C     DATA_PT, U_PARAM, V_PARAM for a richer least-squares fit. Method 1
!C     (Tikhonov) uses only the grid-node data filled above.
!C
!C     For subdivision points inside each quad (cell):
!C       - Parametric domain in u and v: [0,1]
!C       - S,T in [0,1] parameterize position in each quad (cell)
!C         S = local u within quad, T = local v within quad
!C         S = (II_SUB + 0.5) / DIV,    T = (JJ_SUB + 0.5) / DIV
!C       - The global surface parameters U_PARAM, V_PARAM are:
!C           U_PARAM = (I-1 + S) / NX,   V_PARAM = (J-1 + T) / NY
!C         where (I,J) are the quad indices (1-based), NX, NY = surface grid dims
!C
!C     The surface position at each (S,T) is bilinearly interpolated from the four corner nodes of the quad:
!C         X(S,T) = (1-S) * (1-T) * node (I,   J  )
!C                + S     * (1-T) * node (I+1, J  )
!C                + (1-S) * T     * node (I,   J+1)
!C                + S     * T     * node (I+1, J+1)

!C     The goal here is to assemble the normal equations for a least-squares NURBS surface fit.
!C     - ATA  (matrix): Accumulates the outer products of the NURBS basis
!C       functions evaluated at every data point.
!C     - ATD  (matrix): Accumulates the products of the NURBS basis functions and the coordinates (X,Y,Z) of each data point.
!C
!C     In least-squares surface fitting, we are solving for the control points C such that:
!C         minimize || A * C - D ||^2,
!C     where:
!C       - A is the matrix of basis functions evaluated at all sampled (u,v) parametric points,
!C       - D is the corresponding data point coordinates in space,
!C       - C is the set of unknown NURBS control points (to be found).
!C
!C     The normal equations come from setting the derivative to zero and yield:
!C         (A^T * A) * C = (A^T * D)
!C     or expressed in code variables:
!C         ATA * C_CP = ATD
!C
!C     The loops below form ATA and ATD by summing over all sampled points, evaluating how much each basis function is
!C     present at each point (A_ROW), and weighting appropriately by the data's coordinates.
!C     ------------------------------------------------------------------
        IF (IDEBUG_IN > 0) THEN
          WRITE(*,*) ''
          WRITE(*,*) 'Q1NP DEBUG: FITTING METHOD 0 - plain least-squares (subdivision points)'
        ENDIF
!C     Fill DATA_PT with subdivision points (bilinear interpolation in each quad)
        do j = 1, ny
          do i = 1, nx
            DO JJ_SUB=0,DIV_IN-1
              DO II_SUB=0,DIV_IN-1
                S = (REAL(II_SUB,WP) + HALF)/REAL(DIV_IN,WP)
                T = (REAL(JJ_SUB,WP) + HALF)/REAL(DIV_IN,WP)
                IDATA = IDATA + 1

                DATA_PT(1,IDATA) = &
     &      (ONE-S)*(ONE-T)*X_GRID(1,I,J) + &
     &      S*(ONE-T)*X_GRID(1,I+1,J) + &
     &      (ONE-S)*T*X_GRID(1,I,J+1) + &
     &      S*T*X_GRID(1,I+1,J+1)

                DATA_PT(2,IDATA) = &
     &      (ONE-S)*(ONE-T)*X_GRID(2,I,J) + &
     &      S*(ONE-T)*X_GRID(2,I+1,J) + &
     &      (ONE-S)*T*X_GRID(2,I,J+1) + &
     &      S*T*X_GRID(2,I+1,J+1)

                DATA_PT(3,IDATA) = &
     &      (ONE-S)*(ONE-T)*X_GRID(3,I,J) + &
     &      S*(ONE-T)*X_GRID(3,I+1,J) + &
     &      (ONE-S)*T*X_GRID(3,I,J+1) + &
     &      S*T*X_GRID(3,I+1,J+1)

                IF (NX > 0) U_PARAM(IDATA) = (REAL(I-1)+S)/REAL(NX)
                IF (NY > 0) V_PARAM(IDATA) = (REAL(J-1)+T)/REAL(NY)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        IF (IDEBUG_IN > 1) THEN
          WRITE(*,*) ''
          WRITE(*,*) 'Q1NP DEBUG: after filling data points, DIV = ', DIV_IN
          WRITE(*,'(A)') ' DATAPOINT            X        Y        Z        U       V'
          DO I = 1, NDATA
            WRITE(*,'(A,I4,A,2X,F7.2,2X,F7.2,2X,F7.2,2X,F7.2,2X,F7.2)') &
     &      'Q1NP pt(', I, '): ', &
     &      DATA_PT(1,I), DATA_PT(2,I), DATA_PT(3,I), &
     &      U_PARAM(I), V_PARAM(I)
          ENDDO
        ENDIF

!C     Initialize ATA and ATD
      DO COL1=1,NCP
        DO COL2=1,NCP
          ATA(COL1,COL2) = ZERO
        ENDDO
        ATD(COL1,1) = ZERO
        ATD(COL1,2) = ZERO
        ATD(COL1,3) = ZERO
      ENDDO

!C     Update ATA and ATD
      DO IDATA=1,NDATA
        UU = U_PARAM(IDATA)
        VV = V_PARAM(IDATA)

        ! A_ROW is the basis functions at the data point (UU,VV)
        CALL Q1NP_BASIS_ROW_AT_UV(UU, VV, U_KNOT, V_KNOT, P, Q, &
     &      NCP_U, NCP_V, A_ROW)

        ! Update ATA and ATD
        DO COL1=1,NCP
          DO COL2=1,NCP
            ATA(COL1,COL2) = ATA(COL1,COL2) + A_ROW(COL1)*A_ROW(COL2)
          ENDDO
          ATD(COL1,1) = ATD(COL1,1) + A_ROW(COL1)*DATA_PT(1,IDATA)
          ATD(COL1,2) = ATD(COL1,2) + A_ROW(COL1)*DATA_PT(2,IDATA)
          ATD(COL1,3) = ATD(COL1,3) + A_ROW(COL1)*DATA_PT(3,IDATA)
        ENDDO
      ENDDO

!C     Solve normal equations ATA * C_CP = ATD via Cholesky (3 RHS = X,Y,Z)
      CALL CHOLESKY_SOLVE_Q1NP(NCP, ATA, NCP, ATD, 3, C_CP)
      !return value is C_CP, (C_CP = inverse of ATA * ATD)

!C     Write solved control points into Q1NP_CPTAB (CP_MAP already set in Step 4)
      DO CP_COUNTER=1,NCP
        Q1NP_CPTAB(1,CP_COUNTER) = C_CP(CP_COUNTER,1)
        Q1NP_CPTAB(2,CP_COUNTER) = C_CP(CP_COUNTER,2)
        Q1NP_CPTAB(3,CP_COUNTER) = C_CP(CP_COUNTER,3)
      ENDDO

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C   METHOD 1: Tikhonov-regularized least-squares (grid nodes only)
!C   Minimizes ||A*C - D||^2 + lambda*||L*C||^2; L = 1D second-difference on CPs.
!C   Lambda is chosen from a log-spaced set to minimize RMS fit error.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ELSE IF (ICP_METHOD_IN  .EQ. 1) THEN
        IF (IDEBUG_IN > 0) THEN
          WRITE(*,*) ''
          WRITE(*,*) 'Q1NP DEBUG: FITTING METHOD 1 - Tikhonov-regularized least-squares'
        ENDIF

!C       Allocate regularization arrays: L_REG (second-diff), LTL = L^T*L, ATA_REG, C_CP_BEST, TMPVEC
        ALLOCATE(L_REG(NCP,NCP),LTL(NCP,NCP),ATA_REG(NCP,NCP), &
     &      C_CP_BEST(NCP,3),TMPVEC(NCP))

!C       Build normal equations ATA, ATD and data norm DTD (NDATA = (NX+1)*(NY+1) for method 1)
        DO COL1=1,NCP
          DO COL2=1,NCP
            ATA(COL1,COL2) = ZERO
          ENDDO
          ATD(COL1,1) = ZERO
          ATD(COL1,2) = ZERO
          ATD(COL1,3) = ZERO
        ENDDO
        DTD = ZERO

        DO IDATA=1,NDATA
          UU = U_PARAM(IDATA)
          VV = V_PARAM(IDATA)

          CALL Q1NP_BASIS_ROW_AT_UV(UU, VV, U_KNOT, V_KNOT, P, Q, &
     &      NCP_U, NCP_V, A_ROW)

          DO COL1=1,NCP
            DO COL2=1,NCP
              ATA(COL1,COL2) = ATA(COL1,COL2) + &
     &      A_ROW(COL1)*A_ROW(COL2)
            ENDDO
            ATD(COL1,1) = ATD(COL1,1) + A_ROW(COL1)*DATA_PT(1,IDATA)
            ATD(COL1,2) = ATD(COL1,2) + A_ROW(COL1)*DATA_PT(2,IDATA)
            ATD(COL1,3) = ATD(COL1,3) + A_ROW(COL1)*DATA_PT(3,IDATA)
          ENDDO

          DTD = DTD &
     &      + DATA_PT(1,IDATA)*DATA_PT(1,IDATA) &
     &      + DATA_PT(2,IDATA)*DATA_PT(2,IDATA) &
     &      + DATA_PT(3,IDATA)*DATA_PT(3,IDATA)
        ENDDO

!C       L_REG = 1D second-difference matrix (tridiagonal); LTL = L_REG^T * L_REG
        DO COL1=1,NCP
          DO COL2=1,NCP
            L_REG(COL1,COL2) = ZERO
            LTL(COL1,COL2)   = ZERO
          ENDDO
        ENDDO
        DO COL1=1,NCP
          L_REG(COL1,COL1) = -TWO
          IF (COL1 .GT. 1)   L_REG(COL1,COL1-1) = ONE
          IF (COL1 .LT. NCP) L_REG(COL1,COL1+1) = ONE
        ENDDO
        DO COL1=1,NCP
          DO COL2=1,NCP
            DO COL=1,NCP
              LTL(COL1,COL2) = LTL(COL1,COL2) &
     &      + L_REG(COL,COL1)*L_REG(COL,COL2)
            ENDDO
          ENDDO
        ENDDO

!C       Tikhonov solve (ATA + lambda*LTL)*C = ATD; try NLAM log-spaced lambdas in [1e-6, 1e-1]
        NLAM        = 15
        LAMBDA_BEST = -ONE
        BEST_RMS    = -ONE
        DO ILAM = 1, NLAM
          LAMBDA_CUR = 1.0D-6 * &
     &      ( (1.0D-1/1.0D-6) ** &
     &      ( REAL(ILAM-1,kind(ONE)) / REAL(NLAM-1,kind(ONE)) ) )
          DO COL1=1,NCP
            DO COL2=1,NCP
              ATA_REG(COL1,COL2) = &
     &      ATA(COL1,COL2) + LAMBDA_CUR * LTL(COL1,COL2)
            ENDDO
          ENDDO
          CALL CHOLESKY_SOLVE_Q1NP(NCP, ATA_REG, NCP, ATD, 3, C_CP)

!C         Residual ||A*C - D||^2 = C^T*ATA*C - 2*C^T*ATD + DTD (per dimension, then sum)
          RES2 = ZERO
          DO COL=1,3
            DO COL1=1,NCP
              TMPVEC(COL1) = ZERO
            ENDDO
            DO COL1=1,NCP
              DO COL2=1,NCP
                TMPVEC(COL1) = TMPVEC(COL1) + &
     &      ATA(COL1,COL2)*C_CP(COL2,COL)
              ENDDO
            ENDDO
            RES2 = RES2 &
     &      + DOT_PRODUCT(C_CP(1:NCP,COL),TMPVEC(1:NCP)) &
     &      - TWO*DOT_PRODUCT(C_CP(1:NCP,COL),ATD(1:NCP,COL))
          ENDDO
          RES2 = RES2 + DTD
          IF (RES2 .GT. ZERO) THEN
            RMS = SQRT(RES2 / REAL(3*NDATA,kind(ONE)))
          ELSE
            RMS = ZERO
          ENDIF

          IF (BEST_RMS .LT. ZERO .OR. RMS .LT. BEST_RMS) THEN
            BEST_RMS    = RMS
            LAMBDA_BEST = LAMBDA_CUR
            DO COL1=1,NCP
              C_CP_BEST(COL1,1) = C_CP(COL1,1)
              C_CP_BEST(COL1,2) = C_CP(COL1,2)
              C_CP_BEST(COL1,3) = C_CP(COL1,3)
            ENDDO
          ENDIF
        ENDDO

!C       Store best control points in Q1NP_CPTAB (same layout as Method 0)
        DO CP_COUNTER=1,NCP
          Q1NP_CPTAB(1,CP_COUNTER) = C_CP_BEST(CP_COUNTER,1)
          Q1NP_CPTAB(2,CP_COUNTER) = C_CP_BEST(CP_COUNTER,2)
          Q1NP_CPTAB(3,CP_COUNTER) = C_CP_BEST(CP_COUNTER,3)
        ENDDO
        DEALLOCATE(L_REG,LTL,ATA_REG,C_CP_BEST,TMPVEC)
      ENDIF

!C     Deallocate fit and working arrays
      DEALLOCATE(X_GRID, GRID_NODE_RES, DATA_PT, U_PARAM, V_PARAM)
      DEALLOCATE(ATA, ATD, A_ROW, C_CP, U_KNOT, V_KNOT)
      END SUBROUTINE Q1NP_FIT_CONTROL_POINTS
      
      end module genq1np_mod