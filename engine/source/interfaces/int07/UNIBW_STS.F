      SUBROUTINE REMAP_SEGMENTS_STS(INTBUF_TAB, ITAB, X, CAND_SEC_SEG, IRECT, 
     &  CONT_ELEMENT, COUNT, IGRSURF, CAND_SEC_SEG_ID, CAND_MST_SEG_ID,
     &  MAX_STS_SIZE_ACTUAL)
C-----------------------------------------------
C   M o d u l e s
C----------------------------------------------- 
      USE INTBUFDEF_MOD
      USE GROUPDEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
#include      "mvsiz_p.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      TYPE (SURF_)   , DIMENSION(NSURF)   :: IGRSURF
      INTEGER ITAB(*)
      INTEGER IRECT(4,*)
      my_real X(3,*)
      INTEGER CAND_SEC_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_MST_SEG(MAX_STS_SIZE_ACTUAL)
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      my_real CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      INTEGER COUNT, MAX_STS_SIZE_ACTUAL
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I, J, K, N, NI
      INTEGER candidate, candidateM
      LOGICAL :: duplicate
      INTEGER, ALLOCATABLE :: IGRSURF_S_TEMP(:,:), IGRSURF_M_TEMP(:,:)
      INTEGER, ALLOCATABLE :: found_segments(:,:)
      INTEGER, ALLOCATABLE :: ADJA(:,:)
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
      ! Safety checks
      IF (IGRSURF(1)%NSEG <= 0 .OR. IGRSURF(2)%NSEG <= 0) THEN
        COUNT = 0
        RETURN
      END IF
      IF (.NOT. ALLOCATED(IGRSURF(1)%NODES) .OR. 
     &    .NOT. ALLOCATED(IGRSURF(2)%NODES)) THEN
        COUNT = 0
        RETURN
      END IF
      IF (INTBUF_TAB(1)%S_CAND_N <= 0) THEN
        COUNT = 0
        RETURN
      END IF

      ! Copy surface nodes to temporary arrays
      ALLOCATE(IGRSURF_S_TEMP(IGRSURF(1)%NSEG, 4))
      ALLOCATE(IGRSURF_M_TEMP(IGRSURF(2)%NSEG, 4))
      IGRSURF_S_TEMP = IGRSURF(1)%NODES
      IGRSURF_M_TEMP = IGRSURF(2)%NODES

      ! Initialize
      ALLOCATE(found_segments(2, MAX_STS_SIZE_ACTUAL))
      found_segments = -1  ! Initialize with -1 to indicate not found   
      N = 1
      CAND_SEC_SEG = -1
      CAND_MST_SEG = -1  
      ! Compare candidate array with IGRSURF_TEMP table to derive the respective index
      COUNT = 0
      
      ! Map candidate nodes to segment pairs
      DO I = 1, INTBUF_TAB(1)%S_CAND_N
        candidate = INTBUF_TAB(1)%NSV(INTBUF_TAB(1)%CAND_N(I))
        candidateM = INTBUF_TAB(1)%CAND_E(I)
        
        ! Find which secondary segment contains this candidate node
        DO J = 1, IGRSURF(1)%NSEG
          ! Check whether the current candidate node (candidate) is among the four corner nodes
          ! of the current secondary segment (J). If yes, this segment is considered relevant.
          IF (ANY(candidate == IGRSURF_S_TEMP(J, 1:4))) THEN
            
            ! Check if this segment pair already exists
            duplicate = .FALSE.
            DO K = 1, COUNT
              IF (CAND_SEC_SEG(K) == J .AND. CAND_MST_SEG(K) == candidateM) THEN
                duplicate = .TRUE.
                EXIT
              END IF
            END DO

            ! Add new unique segment pair            
            IF (.NOT. duplicate) THEN
              COUNT = COUNT + 1
              IF (COUNT > MAX_STS_SIZE_ACTUAL) THEN
                COUNT = COUNT - 1
                EXIT
              END IF
              
              found_segments(1, N) = J
              found_segments(2, N) = candidateM

              CAND_SEC_SEG(COUNT) = found_segments(1, N)
              CAND_MST_SEG(COUNT) = found_segments(2, N)

              N = N + 1
            END IF
          END IF
        END DO
      END DO
      
      ! Adjacency matrix - filter out invalid entries
      ALLOCATE(ADJA(COUNT, 2))
      J = 0
      DO I = 1, COUNT
        IF (CAND_SEC_SEG(I) /= 0 .AND. CAND_MST_SEG(I) /= 0) THEN
          J = J + 1
          ADJA(J, 1) = CAND_SEC_SEG(I)
          ADJA(J, 2) = CAND_MST_SEG(I)
        END IF
      END DO
      
      COUNT = J
      J = 1
      
      ! Build segment ID arrays
      DO I = 1, COUNT
        ! Secondary segment nodes
        CAND_SEC_SEG_ID(I, 1) = CAND_SEC_SEG(I)
        CAND_SEC_SEG_ID(I, 2:5) = IGRSURF(1)%NODES(CAND_SEC_SEG(I), 1:4)
        
        ! Primary segment nodes
        CAND_MST_SEG_ID(I, 1) = CAND_MST_SEG(I)
        CAND_MST_SEG_ID(I, 2:5) = IRECT(1:4, CAND_MST_SEG(I))
      END DO

      ! Store coordinates: Primary (1-4), Secondary (5-8)
      ! PRIMARY -> FIRST (1-4)
      DO I = 1, COUNT
        J = 1
        DO K = 2, 5
          NI = CAND_MST_SEG_ID(I, K)
          CONT_ELEMENT(I, 1, J) = X(1, NI)  ! X
          CONT_ELEMENT(I, 2, J) = X(2, NI)  ! Y
          CONT_ELEMENT(I, 3, J) = X(3, NI)  ! Z
          J = J + 1
        END DO
      END DO
      
      ! SECONDARY -> SECOND (5-8)
      DO I = 1, COUNT
        J = 5
        DO K = 2, 5
          NI = CAND_SEC_SEG_ID(I, K)
          CONT_ELEMENT(I, 1, J) = X(1, NI)  ! X
          CONT_ELEMENT(I, 2, J) = X(2, NI)  ! Y
          CONT_ELEMENT(I, 3, J) = X(3, NI)  ! Z
          J = J + 1
        END DO
      END DO
      
      ! Cleanup
      DEALLOCATE(IGRSURF_S_TEMP)
      DEALLOCATE(IGRSURF_M_TEMP)
      DEALLOCATE(found_segments)
      DEALLOCATE(ADJA)
      
      RETURN
      END

      
       SUBROUTINE PROCESS_STS_CONTACTS(CONT_ELEMENT, COUNT, OPTION,
     & CAND_MST_SEG_ID, CAND_SEC_SEG_ID,
     & load_arr, node_id_load, L_out, IMPACT_glob, STIF, MAX_STS_SIZE_ACTUAL,
     & FRICC, FRIC_COEFS, VISCFFRIC, XMU, MFROT, IFQ, DT1, DT12,
     & V, MS, CAND_F, ALPHA0, IFPEN, INTTH, QFRICT, INTBUF_TAB, GAP,
     & XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST,
     & GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     & GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     & GP_XI1_PERIOD, GP_XI2_PERIOD,
     & GP_INITIALIZED,
     & GP_TTRIAL1_HIST, GP_TTRIAL2_HIST,
     & GP_IS_STICKING)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE INTBUFDEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
#include      "mvsiz_p.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      REAL*8 CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      my_real STIF(MVSIZ)
      INTEGER COUNT, OPTION
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      REAL*8 load_arr(MAX_STS_SIZE_ACTUAL,8,4)
      INTEGER node_id_load(MAX_STS_SIZE_ACTUAL*8)
      INTEGER L_out, IMPACT_glob, MAX_STS_SIZE_ACTUAL
      my_real FRICC(MVSIZ), FRIC_COEFS(MVSIZ,10), VISCFFRIC(MVSIZ)
      my_real XMU(MVSIZ), DT1, DT12, ALPHA0
      my_real V(3,*), MS(*)
      my_real CAND_F(8,*)
      INTEGER IFPEN(*)
      INTEGER MFROT, IFQ, INTTH
      my_real QFRICT
      my_real GAP  ! Gap value from user input
      REAL*8 XI1_HIST(MAX_STS_SIZE_ACTUAL,2,2)  ! History of xi1 per Gauss point
      REAL*8 XI2_HIST(MAX_STS_SIZE_ACTUAL,2,2)  ! History of xi2 per Gauss point
      REAL*8 TTRIAL1_HIST(MAX_STS_SIZE_ACTUAL,2,2)  ! History of T_trial(1) per Gauss point
      REAL*8 TTRIAL2_HIST(MAX_STS_SIZE_ACTUAL,2,2)  ! History of T_trial(2) per Gauss point
      REAL*8 GP_XI1_GLOBAL(*)
      REAL*8 GP_XI2_GLOBAL(*)
      REAL*8 GP_XI1_GLOBAL_PREV(*)
      REAL*8 GP_XI2_GLOBAL_PREV(*)
      INTEGER GP_XI1_PERIOD(*)
      INTEGER GP_XI2_PERIOD(*)
      LOGICAL GP_INITIALIZED(*)
      REAL*8 GP_TTRIAL1_HIST(*)
      REAL*8 GP_TTRIAL2_HIST(*)
      LOGICAL GP_IS_STICKING(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I, J, K, L, IMPACT
      REAL*8 XUPD(3,8)
      REAL*8 p_load_new(24)
      REAL*8 node_stiff(8)
      REAL*8 p_friction(24)  ! Friction forces (separate output)
      my_real EFRICT_LOC
      INTEGER node_ids(8)  ! Node IDs for velocity interpolation
C-----------------------------------------------
C   I n i t i a l i z a t i o n
C-----------------------------------------------
      IMPACT_glob = 0
      
      ! Initialize counters
      K = 1
      L = 1
C-----------------------------------------------
C   M a i n   L o o p
C-----------------------------------------------
      DO I = 1, COUNT
        IMPACT = 0
        XUPD = CONT_ELEMENT(I, 1:3, 1:8)
      
        ! Get node IDs for velocity interpolation
        DO J = 1, 4
          node_ids(J)   = CAND_MST_SEG_ID(I, J+1)   ! Primary nodes
          node_ids(J+4) = CAND_SEC_SEG_ID(I, J+1)   ! Secondary nodes
        ENDDO
        XMU(1) = FRICC(1) ! Friction coefficient mu
        
        ! Call GEN_STS with friction calculation integrated
        ! Note: STIF is used for both normal and tangential stiffness
        ! (STIF0 = STIF for now, until separate tangential stiffness is available)
        CALL GEN_STS(XUPD, STIF, p_load_new, IMPACT, I, 
     &                    node_stiff, OPTION,
     &                    V, MS, FRICC, FRIC_COEFS, VISCFFRIC, XMU, 
     &                    MFROT, IFQ, ALPHA0, CAND_F, IFPEN, STIF,
     &                    p_friction, EFRICT_LOC, QFRICT, 
     &                    INTTH, node_ids, .TRUE.,
     &                    XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST, 
     &                    MAX_STS_SIZE_ACTUAL, GAP,
     &                    CAND_SEC_SEG_ID,
     &                    GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     &                    GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     &                    GP_XI1_PERIOD, GP_XI2_PERIOD,
     &                    GP_INITIALIZED,
     &                    GP_TTRIAL1_HIST, GP_TTRIAL2_HIST,
     &                    GP_IS_STICKING)
      
        IF (IMPACT == 1) THEN
          IMPACT_glob = 1
          
          ! Save node IDs: Primary (1-4), Secondary (5-8)
          node_id_load(K:K+3) = CAND_MST_SEG_ID(I, 2:5)
          node_id_load(K+4:K+7) = CAND_SEC_SEG_ID(I, 2:5)
          K = K + 8
      
          ! Store forces: Primary (1-4), Secondary (5-8)
          DO J = 1, 4
            ! Primary forces
            load_arr(L, J, 1:3) = p_load_new(3*(J-1) + 1 : 3*J)
            ! Secondary forces
            load_arr(L, J + 4, 1:3) = p_load_new(12 + 3*(J-1) + 1 : 12 + 3*J)
          ENDDO
      
          ! Store stiffness info for all 8 nodes
          load_arr(L, 1:8, 4) = node_stiff(1:8)
      
          L = L + 1
          
        ENDIF
      ENDDO
      
      L_out = L
      END SUBROUTINE PROCESS_STS_CONTACTS
      
************************************************************************
c----------------------------------------------------------------------c
c.... general subroutines
c----------------------------------------------------------------------c
************************************************************************

      INTEGER FUNCTION GET_GLOBAL_GP_INDEX(PAIR_INDEX, Z, Q, IP_MAX)
C-----------------------------------------------------------------------
C Compute global Gauss point index from STS pair index and
C local Gauss point indices (z, q) for a given quadrature order IP_MAX.
C
C gp_index = (pair_index-1) * (IP_MAX*IP_MAX) + local_gp_index + 1
C local_gp_index = (z-1)*IP_MAX + (q-1)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER PAIR_INDEX, Z, Q, IP_MAX
      INTEGER LOCAL_GP_INDEX

      LOCAL_GP_INDEX = (Z - 1) * IP_MAX + (Q - 1)
      GET_GLOBAL_GP_INDEX = (PAIR_INDEX - 1) * (IP_MAX*IP_MAX)
     &                    + LOCAL_GP_INDEX + 1

      RETURN
      END

      subroutine GEN_STS(XUPD, STIF, p, IMPACT, EL_NR, node_stiff, OPTION,
     &                        V, MS, FRICC, FRIC_COEFS, VISCFFRIC, XMU, MFROT, IFQ,
     &                        ALPHA0, CAND_F, IFPEN, STIF0,
     &                        p_friction, EFRICT, QFRICT, INTTH, node_ids,
     &                        CALC_FRICTION, XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST, 
     &                        MAX_STS_SIZE, GAP,
     &                        CAND_SEC_SEG_ID,
     &                        GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     &                        GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     &                        GP_XI1_PERIOD, GP_XI2_PERIOD,
     &                        GP_INITIALIZED,
     &                        GP_TTRIAL1_HIST, GP_TTRIAL2_HIST,
     &                        GP_IS_STICKING)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
#include      "mvsiz_p.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
#include      "com08_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     XUPD     : Coordinates of contact element (3,8)
C                Nodes 1-4 are Primary (master) nodes
C                Nodes 5-8 are Secondary (slave) nodes
C     STIF     : Contact stiffness parameter (normal)
C     p        : Output contact forces (24 components) - normal + friction combined
C     IMPACT   : Output flag indicating if penetration was detected
C     EL_NR    : Temporary Element number (for debugging)
C     node_stiff: Output nodal stiffness values (8 components)
C     OPTION   : 0 = Gauss quadrature, 1 = Lobatto quadrature
C     V        : Nodal velocities (3,*)
C     MS       : Nodal masses (*)
C     FRICC, FRIC_COEFS, VISCFFRIC, XMU: Friction parameters
C     MFROT    : Friction model type
C     IFQ      : Friction formulation flag
C     ALPHA0   : Friction parameter
C     CAND_F   : Stored friction forces (for incremental formulation)
C     IFPEN    : Penetration flag array
C     STIF0    : Tangential contact stiffness
C     p_friction: Output friction forces (24 components) - separate output
C     EFRICT   : Friction energy (output)
C     QFRICT   : Total friction energy (INTENT(INOUT))
C     INTTH    : Thermal interface flag
C     node_ids : Node IDs for velocity interpolation (8 components)
C     CALC_FRICTION: Flag to enable/disable friction calculation
C     XI1_HIST : History of xi1 parametric coordinates per Gauss point
C     XI2_HIST : History of xi2 parametric coordinates per Gauss point
C     MAX_STS_SIZE: Maximum size for history arrays
C     GAPV: Gap value from user input
C-----------------------------------------------
      INTEGER IMPACT, OPTION, EL_NR, MFROT, IFQ, INTTH
      INTEGER node_ids(8)
      LOGICAL CALC_FRICTION
      my_real STIF(MVSIZ), STIF0(MVSIZ)
      real*8  p(24), p_friction(24)
      real*8  XUPD(3,8)
      real*8  node_stiff(8)
      my_real V(3,*), MS(*)
      my_real FRICC(MVSIZ), FRIC_COEFS(MVSIZ,10), VISCFFRIC(MVSIZ)
      my_real XMU(MVSIZ), ALPHA0
      my_real CAND_F(8,*)
      INTEGER IFPEN(*)
      my_real EFRICT, QFRICT
      INTEGER MAX_STS_SIZE  ! Maximum size for history arrays
      REAL*8 XI1_HIST(MAX_STS_SIZE,2,2)  ! History of xi1 per Gauss point
      REAL*8 XI2_HIST(MAX_STS_SIZE,2,2)  ! History of xi2 per Gauss point
      REAL*8 TTRIAL1_HIST(MAX_STS_SIZE,2,2)  ! History of T_trial(1) per Gauss point
      REAL*8 TTRIAL2_HIST(MAX_STS_SIZE,2,2)  ! History of T_trial(2) per Gauss poin
      my_real GAP  ! Gap value from user input
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE,5)
      REAL*8  GP_XI1_GLOBAL(*)
      REAL*8  GP_XI2_GLOBAL(*)
      REAL*8  GP_XI1_GLOBAL_PREV(*)
      REAL*8  GP_XI2_GLOBAL_PREV(*)
      INTEGER GP_XI1_PERIOD(*)
      INTEGER GP_XI2_PERIOD(*)
      LOGICAL GP_INITIALIZED(*)
      REAL*8  GP_TTRIAL1_HIST(*)
      REAL*8  GP_TTRIAL2_HIST(*)
      LOGICAL GP_IS_STICKING(*)
C     interface to global gp index function
      INTEGER GET_GLOBAL_GP_INDEX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER i, j, z, q, ip
      real*8  xi1, xi2
      real*8  penetr, PENE, GAPV, FAC
      my_real d1
      real*8  a(3,24), daxi1(3,24), daxi2(3,24)
      real*8  daeta1(3,24), daeta2(3,24)
      real*8  rhoxi1(3), rhoxi2(3)
      real*8  m_ij(2,2), detm, mij(2,2), detmPrimary
      real*8  norm(3)
      real*8  pm(24), pm_friction(24)
      real*8  eta1(10), eta2(10), wi1(10), wi2(10)
      real*8  energy
      real*8  N_xi(3,4), N_eta(3,4)
      real*8  VX, VY, VZ, FTN
      real*8  FXT, FYT, FZT, PHI, FN
      INTEGER INDEX_CAND
      INTEGER gp_index, secondary_el_id
      INTEGER MAX_GLOBAL_GP_LOCAL

      INTEGER, PARAMETER :: DEBUG_NODE = 1
      CHARACTER*7 FRICTION_QUAD_TYPE
      CHARACTER*7 STICK_STATE

      real*8  dxi1, dxi2  ! Convective coordinate increments (velocity-based)
      real*8  T_trial(2), T_real(2), T_trialabs
      
      ! Gauss quadrature for friction calculation
      real*8  eta1_gauss(10), eta2_gauss(10), wi1_gauss(10), wi2_gauss(10)
      real*8  xi1_gauss, xi2_gauss
      real*8  m_ij_gauss(2,2), detm_gauss
      real*8  rhoxi1_gauss(3), rhoxi2_gauss(3)
      real*8  a_gauss(3,24), daxi1_gauss(3,24), daxi2_gauss(3,24)
      real*8  daeta1_gauss(3,24), daeta2_gauss(3,24)
      logical gauss_valid
      
      ! Friction calculation variables (selected projection)
      real*8  xi1_fric, xi2_fric
      real*8  m_ij_fric(2,2), detm_fric
      real*8  rhoxi1_fric(3), rhoxi2_fric(3)
      real*8  wi1_fric, wi2_fric, eta1_fric, eta2_fric
      logical use_gauss_for_friction
C-----------------------------------------------
C   I n i t i a l i z a t i o n
C-----------------------------------------------
      IMPACT = 0
      ip = 2 ! Quadrature order
      
      ! Calculate maximum global GP index for bounds checking
      MAX_GLOBAL_GP_LOCAL = MAX_STS_SIZE * ip * ip
      
      ! Get quadrature points and weights
      
      ! Always initialize Gauss quadrature for friction calculation
      call GaussPt(ip, eta1_gauss, wi1_gauss)
      call GaussPt(ip, eta2_gauss, wi2_gauss)    
      
      ! Initialize quadrature points and weights for normal contact calculation
      IF (OPTION == 0) THEN
        ! Gauss: reuse Gauss quadrature points and weights
        eta1 = eta1_gauss
        eta2 = eta2_gauss
        wi1 = wi1_gauss
        wi2 = wi2_gauss
      ELSE
        ! Lobatto: use Lobatto quadrature points and weights
        call LobattoPt(ip, eta1, wi1)
        call LobattoPt(ip, eta2, wi2)
      ENDIF

      ! Initialize force arrays
      DO i=1,24
        pm(i) = 0.d0
        pm_friction(i) = 0.d0
        p_friction(i) = 0.d0
      ENDDO
      GAPV = GAP ! Use user-defined Gapmin value
      energy = 0.0d0
      EFRICT = 0.d0
      XMU(1) = FRICC(1) ! Friction coefficient mu
C-----------------------------------------------
C   M a i n   C o m p u t a t i o n
C-----------------------------------------------
C     Loop over integration points
      DO z=1,ip
        DO q=1,ip
          
          ! Project Secondary surface to Primary surface at current Gauss/Lobatto point
          call project112(XUPD, xi1, xi2, eta1(z), eta2(q))
          
          ! Check if projection is valid
          IF ((dabs(xi1) .GT. 1.05d0) .OR. (dabs(xi2) .GT. 1.05d0)) THEN
            CYCLE
          ENDIF
          
          ! Build position and derivative matrices
          call pos112(a, daxi1, daxi2, daeta1, daeta2, xi1, xi2, 
     &               eta1(z), eta2(q))

          ! Calculate surface geometry and metrics
          call surfgeom112(XUPD, daxi1, daxi2, daeta1, daeta2, norm,
     &                    rhoxi1, rhoxi2, m_ij, detm, mij, detmPrimary)
          
          ! ==== GAUSS PROJECTION FOR FRICTION ====
          IF (OPTION == 0) THEN
            ! Gauss algorithm: current projection IS Gauss projection
            ! Reuse coordinates and geometry
            xi1_gauss = xi1
            xi2_gauss = xi2
            gauss_valid = .TRUE.
            
            ! Reuse geometry arrays (already calculated above)
            DO i = 1, 3
              DO j = 1, 24
                a_gauss(i,j) = a(i,j)
                daxi1_gauss(i,j) = daxi1(i,j)
                daxi2_gauss(i,j) = daxi2(i,j)
                daeta1_gauss(i,j) = daeta1(i,j)
                daeta2_gauss(i,j) = daeta2(i,j)
              ENDDO
            ENDDO
            DO i = 1, 3
              rhoxi1_gauss(i) = rhoxi1(i)
              rhoxi2_gauss(i) = rhoxi2(i)
            ENDDO
            m_ij_gauss(1,1) = m_ij(1,1)
            m_ij_gauss(1,2) = m_ij(1,2)
            m_ij_gauss(2,1) = m_ij(2,1)
            m_ij_gauss(2,2) = m_ij(2,2)
            detm_gauss = detm
          ELSE
            ! Lobatto algorithm: calculate separate Gauss projection
            call project112(XUPD, xi1_gauss, xi2_gauss, 
     &                     eta1_gauss(z), eta2_gauss(q))
            
            ! Check if Gauss projection is valid for friction calculation
            gauss_valid = (dabs(xi1_gauss) .LE. 1.05d0 .AND. 
     &                     dabs(xi2_gauss) .LE. 1.05d0)
            
            IF (gauss_valid) THEN
              ! Calculate surface geometry for Gauss projection
              call pos112(a_gauss, daxi1_gauss, daxi2_gauss, 
     &                   daeta1_gauss, daeta2_gauss, 
     &                   xi1_gauss, xi2_gauss, 
     &                   eta1_gauss(z), eta2_gauss(q))
              
              call surfgeom112(XUPD, daxi1_gauss, daxi2_gauss, 
     &                        daeta1_gauss, daeta2_gauss, 
     &                        norm, rhoxi1_gauss, rhoxi2_gauss, 
     &                        m_ij_gauss, detm_gauss, mij, detmPrimary)
            ENDIF
          ENDIF
          
          ! ==== Normal impact =====
          ! Compute penetration at Gauss point
          call penetr112(XUPD, penetr, norm, a)
          
          ! Check for penetration
          PENE = penetr - GAPV
          
          ! No penetration - skip to next Gauss point
          IF (PENE .GT. 0.d0) THEN
            CYCLE
          ENDIF
          
          ! Penetration detected
          IMPACT = 1
          penetr = PENE

          ! Calculate penalty parameter
          d1 = STIF(1) * 1.0d0
          FAC = GAPV / MAX(EM10, (GAPV+PENE))
          d1 = 0.5d0 * d1 * FAC

          ! Extract nodal stiffness and set penetration flag if penetration detected
          IF (IMPACT == 1) THEN
            node_stiff = d1 ! TODO: Review
            INDEX_CAND = EL_NR ! INDEX_CAND
            IFPEN(INDEX_CAND) = 1  ! Penetration flag (1 = penetrated, 0 = not penetrated)
          ENDIF      
          
          ! Accumulate energy
          energy = energy + 0.5d0 * d1 * penetr**2 * wi1(z) * wi2(q) * 
     &            dsqrt(detm)
          
          ! Compute residual forces (normal component)
          DO i=1,24
            DO j=1,3
              pm(i) = pm(i) + d1 * penetr * 
     &                a(j,i) * norm(j) * wi1(z) * wi2(q) * dsqrt(detm)
            ENDDO
          ENDDO
          !XMU(1) = 0.6
          ! ------------------------------------------------------------------          
          ! ===== FRICTION CALCULATION =====
          IF (XMU(1) .GT. 0.0d0) THEN
            
            ! Map to global Gauss point index for history-based friction
            gp_index = GET_GLOBAL_GP_INDEX(EL_NR, z, q, ip)
            if (gp_index .LE. 0 .OR. gp_index .GT. MAX_GLOBAL_GP_LOCAL) then
              ! Invalid index - skip friction calculation for this GP
              CYCLE
            endif

            ! Determine which projection to use for friction
            ! Set values to appropriate projection (Gauss preferred, fallback to current)
            IF (gauss_valid) THEN
              ! ==== USE GAUSS PROJECTION (preferred) ====
              use_gauss_for_friction = .TRUE.
              xi1_fric = xi1_gauss
              xi2_fric = xi2_gauss
              m_ij_fric(1,1) = m_ij_gauss(1,1)
              m_ij_fric(1,2) = m_ij_gauss(1,2)
              m_ij_fric(2,1) = m_ij_gauss(2,1)
              m_ij_fric(2,2) = m_ij_gauss(2,2)
              detm_fric = detm_gauss
              rhoxi1_fric(1) = rhoxi1_gauss(1)
              rhoxi1_fric(2) = rhoxi1_gauss(2)
              rhoxi1_fric(3) = rhoxi1_gauss(3)
              rhoxi2_fric(1) = rhoxi2_gauss(1)
              rhoxi2_fric(2) = rhoxi2_gauss(2)
              rhoxi2_fric(3) = rhoxi2_gauss(3)
              wi1_fric = wi1_gauss(z)
              wi2_fric = wi2_gauss(q)
              eta1_fric = eta1_gauss(z)
              eta2_fric = eta2_gauss(q)
            ELSE IF (dabs(xi1) .LE. 1.05d0 .AND. dabs(xi2) .LE. 1.05d0) THEN
              ! ==== FALLBACK: USE CURRENT PROJECTION ====
              use_gauss_for_friction = .FALSE.
              xi1_fric = xi1
              xi2_fric = xi2
              m_ij_fric(1,1) = m_ij(1,1)
              m_ij_fric(1,2) = m_ij(1,2)
              m_ij_fric(2,1) = m_ij(2,1)
              m_ij_fric(2,2) = m_ij(2,2)
              detm_fric = detm
              rhoxi1_fric(1) = rhoxi1(1)
              rhoxi1_fric(2) = rhoxi1(2)
              rhoxi1_fric(3) = rhoxi1(3)
              rhoxi2_fric(1) = rhoxi2(1)
              rhoxi2_fric(2) = rhoxi2(2)
              rhoxi2_fric(3) = rhoxi2(3)
              wi1_fric = wi1(z)
              wi2_fric = wi2(q)
              eta1_fric = eta1(z)
              eta2_fric = eta2(q)
            ELSE
              ! Both projections invalid - skip friction for this GP
              CYCLE
            ENDIF

            ! ==== COMMON FRICTION CALCULATION (using selected projection) ====
            ! Calculate convective coordinate increments
            call tangentvel118_global(xi1_fric, xi2_fric, dxi1, dxi2,
     &          EL_NR, z, q, ip, MAX_STS_SIZE, MAX_GLOBAL_GP_LOCAL,
     &          GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     &          GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     &          GP_XI1_PERIOD, GP_XI2_PERIOD,
     &          GP_INITIALIZED,
     &          GP_TTRIAL1_HIST, GP_TTRIAL2_HIST)

            ! Calculate tangential traction using selected metric tensor
            T_trial(1) = GP_TTRIAL1_HIST(gp_index) - 
     &                   d1*(m_ij_fric(1,1)*dxi1 + m_ij_fric(1,2)*dxi2)
            T_trial(2) = GP_TTRIAL2_HIST(gp_index) - 
     &                   d1*(m_ij_fric(2,1)*dxi1 + m_ij_fric(2,2)*dxi2)

            ! Magnitude of trial tangential traction using selected metric
            T_trialabs = 0.0d0
            DO i = 1, 2
              DO j = 1, 2
                T_trialabs = T_trialabs + T_trial(i) * T_trial(j) * 
     &                       m_ij_fric(i, j)
              END DO
            END DO
            T_trialabs = DSQRT(MAX(EM30, T_trialabs))

            ! ===== FRICTION YIELD FUNCTION =====
            FN = d1 * DABS(penetr)  ! Normal force at the Gauss point from penalty method
            PHI = T_trialabs - XMU(1) * FN

            ! ===== Sticking/Sliding logic =====
            IF (PHI .LE. 0.0d0) THEN
              ! ===== STICKING =====
              T_real(1) = T_trial(1)
              T_real(2) = T_trial(2)
              GP_IS_STICKING(gp_index) = .TRUE.
            ELSE
              ! ===== SLIDING =====
              T_real(1) = XMU(1) * FN * T_trial(1) / MAX(1.d-30, T_trialabs)
              T_real(2) = XMU(1) * FN * T_trial(2) / MAX(1.d-30, T_trialabs)
              GP_IS_STICKING(gp_index) = .FALSE.
            ENDIF

            ! Update global history of tangential traction with final value
            GP_TTRIAL1_HIST(gp_index) = T_real(1)
            GP_TTRIAL2_HIST(gp_index) = T_real(2)

            ! ===== CONVERT FROM 2D TANGENT PLANE TO 3D =====
            ! Convert using selected tangent vectors
            FXT = T_real(1)*rhoxi1_fric(1) + T_real(2)*rhoxi2_fric(1)
            FYT = T_real(1)*rhoxi1_fric(2) + T_real(2)*rhoxi2_fric(2)
            FZT = T_real(1)*rhoxi1_fric(3) + T_real(2)*rhoxi2_fric(3)

            ! ===== DEBUG OUTPUT FOR SELECTED NODE =====
            IF (DEBUG_NODE .GT. 0) THEN
              STICK_STATE = 'SLIDE  '
              IF (GP_IS_STICKING(gp_index)) STICK_STATE = 'STICK  '
              IF (use_gauss_for_friction) THEN
                FRICTION_QUAD_TYPE = 'GAUSS  '
              ELSE
                FRICTION_QUAD_TYPE = 'LOBATTO'
              ENDIF
              DO j = 1, 8
                IF (node_ids(j) .EQ. DEBUG_NODE) THEN
                  WRITE(6,1000) EL_NR, node_ids(j), z, q,
     &                 FRICTION_QUAD_TYPE,
     &                 xi1_fric, xi2_fric, dxi1, dxi2,
     &                 GP_XI1_PERIOD(gp_index),
     &                 GP_XI2_PERIOD(gp_index),
     &                 FXT, FYT, FZT,
     &                 STICK_STATE
                ENDIF
              ENDDO
            ENDIF

            ! ===== ACCUMULATE FRICTION FORCES TO NODES =====
            call shape112(xi1_fric, xi2_fric, N_xi)
            call shape112(eta1_fric, eta2_fric, N_eta)
            ! Primary nodes (1-4): subtract friction
            DO j=1,4
              pm_friction((j-1)*3+1) = pm_friction((j-1)*3+1) - 
     &                               N_xi(1,j) * FXT * wi1_fric * 
     &                               wi2_fric * dsqrt(detm_fric)
              pm_friction((j-1)*3+2) = pm_friction((j-1)*3+2) - 
     &                               N_xi(1,j) * FYT * wi1_fric * 
     &                               wi2_fric * dsqrt(detm_fric)
              pm_friction((j-1)*3+3) = pm_friction((j-1)*3+3) - 
     &                               N_xi(1,j) * FZT * wi1_fric * 
     &                               wi2_fric * dsqrt(detm_fric)
            ENDDO
            
            ! Secondary nodes (5-8): add friction
            DO j=1,4
              pm_friction(12+(j-1)*3+1) = pm_friction(12+(j-1)*3+1) + 
     &                                  N_eta(1,j) * FXT * wi1_fric * 
     &                                  wi2_fric * dsqrt(detm_fric)
              pm_friction(12+(j-1)*3+2) = pm_friction(12+(j-1)*3+2) + 
     &                                  N_eta(1,j) * FYT * wi1_fric * 
     &                                  wi2_fric * dsqrt(detm_fric)
              pm_friction(12+(j-1)*3+3) = pm_friction(12+(j-1)*3+3) + 
     &                                  N_eta(1,j) * FZT * wi1_fric * 
     &                                  wi2_fric * dsqrt(detm_fric)
            ENDDO
            
            ! Energy calculation
            ! Calculate tangential velocity from relative velocity
            ! Project V_rel onto tangent plane
            ! TODO: Review this energy calculation
            !FTN = V_rel(1)*norm(1) + V_rel(2)*norm(2) + V_rel(3)*norm(3)
            !VX = V_rel(1) - FTN*norm(1)
            !VY = V_rel(2) - FTN*norm(2)
            !VZ = V_rel(3) - FTN*norm(3)
            !EFRICT = EFRICT + DT1 * (VX*FXT + VY*FYT + VZ*FZT) * 
            !&              wi1(z) * wi2(q) * dsqrt(detm)
          ENDIF
          ! ===== END FRICTION CALCULATION =====
        ENDDO
      ENDDO
                      
C-----------------------------------------------
C   F i n a l   R e s u l t s
C-----------------------------------------------
      ! Set output forces
      ! Combine normal and friction forces
      DO i=1,24
        p(i) = -pm(i) + pm_friction(i)
        p_friction(i) = pm_friction(i)
      ENDDO
      
      ! Update total friction energy
      IF (CALC_FRICTION) THEN
        QFRICT = QFRICT + EFRICT
      ENDIF

 1000 FORMAT('STS-DBG EL=',I6,' NODE=',I10,' z=',I2,' q=',I2,
     &       ' quad=',A7,
     &       ' xi=',2F12.6,
     &       ' dxi=',2F12.6,
     &       ' per=',2I6,
     &       ' FT=',3E12.4,
     &       ' state=',A7)

      RETURN
      END

C-----------------------------------------------
C Project point to surface parameterization
C-----------------------------------------------         
      subroutine project112(xupd,xi1,xi2,eta1,eta2)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     xupd   : Coordinates of contact element (3,8)
C     xi1,xi2: Output parametric coordinates on Primary surface
C     eta1,eta2: Input parametric coordinates on Secondary surface
C-----------------------------------------------          
      real*8  xupd(3,8)     
      real*8  xi1,xi2,eta1,eta2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER i, j, iter
      real*8  shape(3,4), N_eta(3,4)
      real*8  rho(3), rhoxi1(3), rhoxi2(3)
      real*8  xsl(3), dxi1, dxi2
      real*8  e, f(2)
      real*8  a(3,24), daxi1(3,24), daxi2(3,24)
      real*8  daeta1(3,24), daeta2(3,24)
      real*8  norm(3)
      real*8  m_ij(2,2), detm, detmPrimary, mij(2,2)
C-----------------------------------------------
C   Initialization - set initial guess
C-----------------------------------------------
      xi1 = 0.d0
      xi2 = 0.d0
C-----------------------------------------------
C   Get shape functions for Secondary surface
C-----------------------------------------------
      call shape112(eta1, eta2, N_eta)
C-----------------------------------------------
C   Compute Secondary point position in global coords
C-----------------------------------------------
      xsl(1) = 0.d0
      xsl(2) = 0.d0
      xsl(3) = 0.d0
      DO j=1,4
        xsl(1) = xsl(1) + N_eta(1,j)*xupd(1,j+4)
        xsl(2) = xsl(2) + N_eta(1,j)*xupd(2,j+4)
        xsl(3) = xsl(3) + N_eta(1,j)*xupd(3,j+4)
      ENDDO

C-----------------------------------------------
C   Newton's iteration for projection
C-----------------------------------------------
      DO iter=1,1
        ! Get shape functions and derivatives at current xi
        call shape112(xi1, xi2, shape)
         
        ! Setup position and derivatives matrices
        call pos112(a, daxi1, daxi2, daeta1, daeta2, xi1, xi2, eta1, eta2)
         
        ! Calculate surface geometry and metrics
        call surfgeom112(xupd, daxi1, daxi2, daeta1, daeta2, norm, 
     &                   rhoxi1, rhoxi2, m_ij, detm, mij, detmPrimary)

        ! Compute current point on Primary surface
        DO i=1,3
          rho(i) = 0.d0
          DO j=1,4
            rho(i) = rho(i) + xupd(i,j)*shape(1,j)
          ENDDO
        ENDDO
         
        ! Compute projection residual and sensitivities
        e = 0.d0
        f(1) = 0.d0
        f(2) = 0.d0
        DO i=1,3
          ! Surface curvature term
          e = e + (xupd(i,1)-xupd(i,2)+xupd(i,3)-xupd(i,4))*
     &        (xsl(i)-rho(i))*0.25d0
          ! Tangential distance terms
          f(1) = f(1) + (xsl(i)-rho(i))*rhoxi1(i)
          f(2) = f(2) + (xsl(i)-rho(i))*rhoxi2(i)
        ENDDO
      
        ! Update determinant with curvature correction
        detmPrimary = detmPrimary - e**2 + 2.d0*m_ij(1,2)*e
        
        ! Calculate Newton update steps
        dxi1 = (m_ij(2,2)*f(1) + (e-m_ij(1,2))*f(2))/detmPrimary
        dxi2 = (m_ij(1,1)*f(2) + (e-m_ij(2,1))*f(1))/detmPrimary
        
        ! Update solution
        xi1 = xi1 + dxi1
        xi2 = xi2 + dxi2
      ENDDO
      
      RETURN
      END

C-----------------------------------------------
C Shape functions for quad element
C-----------------------------------------------
C-----------------------------------------------
C Shape functions for quad element
C-----------------------------------------------
      subroutine shape112(xi1, xi2, shape)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     xi1,xi2: Parametric coordinates (-1 to 1)
C     shape: Output shape functions and derivatives
C            shape(1,i): Shape function values
C            shape(2,i): Derivatives w.r.t. xi1
C            shape(3,i): Derivatives w.r.t. xi2
C-----------------------------------------------
      real*8  shape(3,4)
      real*8  xi1, xi2
C-----------------------------------------------
C   Bilinear shape functions for quadrilateral element
C-----------------------------------------------
      ! Shape function values - standard bilinear quad
      shape(1,1) = 0.25d0*(1.d0-xi1)*(1.d0-xi2)  ! N1 at (-1,-1)
      shape(1,2) = 0.25d0*(1.d0+xi1)*(1.d0-xi2)  ! N2 at (1,-1)
      shape(1,3) = 0.25d0*(1.d0+xi1)*(1.d0+xi2)  ! N3 at (1,1)
      shape(1,4) = 0.25d0*(1.d0-xi1)*(1.d0+xi2)  ! N4 at (-1,1)

      ! Derivatives w.r.t. xi1
      shape(2,1) = -0.25d0*(1.d0-xi2)  ! dN1/dxi1
      shape(2,2) =  0.25d0*(1.d0-xi2)  ! dN2/dxi1
      shape(2,3) =  0.25d0*(1.d0+xi2)  ! dN3/dxi1
      shape(2,4) = -0.25d0*(1.d0+xi2)  ! dN4/dxi1

      ! Derivatives w.r.t. xi2
      shape(3,1) = -0.25d0*(1.d0-xi1)  ! dN1/dxi2
      shape(3,2) = -0.25d0*(1.d0+xi1)  ! dN2/dxi2
      shape(3,3) =  0.25d0*(1.d0+xi1)  ! dN3/dxi2
      shape(3,4) =  0.25d0*(1.d0-xi1)  ! dN4/dxi2

      RETURN
      END
  
C-----------------------------------------------
C Position matrix assembly
C-----------------------------------------------
      subroutine pos112(a, daxi1, daxi2, daeta1, daeta2, xi1, xi2, eta1, eta2)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     a, daxi1, daxi2, daeta1, daeta2: Output matrices for contact calculation
C     xi1, xi2: Parametric coordinates on Primary surface 
C     eta1, eta2: Parametric coordinates on Secondary surface
C-----------------------------------------------
      real*8  a(3,24), daxi1(3,24), daxi2(3,24), daeta1(3,24), daeta2(3,24)
      real*8  xi1, xi2, eta1, eta2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real*8  N_xi(3,4), N_eta(3,4)
      INTEGER i, j
C-----------------------------------------------
C   Get shape functions at parametric coordinates
C-----------------------------------------------
      call shape112(xi1, xi2, N_xi)
      call shape112(eta1, eta2, N_eta)
      
C-----------------------------------------------
C   Initialize all arrays to zero
C-----------------------------------------------
      DO i=1,3
        DO j=1,24
          a(i,j) = 0.d0
          daxi1(i,j) = 0.d0
          daxi2(i,j) = 0.d0
          daeta1(i,j) = 0.d0
          daeta2(i,j) = 0.d0
        ENDDO
      ENDDO
      
C-----------------------------------------------
C   Fill position and derivative matrices
C-----------------------------------------------
      ! Primary nodes (1-4): position matrix and derivatives
      DO i=1,4
        ! Position matrix (x, y, z components)
        a(1,(i-1)*3+1) = -N_xi(1,i)
        a(2,(i-1)*3+2) = -N_xi(1,i)
        a(3,(i-1)*3+3) = -N_xi(1,i)
        
        ! Derivatives w.r.t. xi1 and xi2
        daxi1(1,(i-1)*3+1) = -N_xi(2,i)
        daxi2(1,(i-1)*3+1) = -N_xi(3,i)
        daxi1(2,(i-1)*3+2) = -N_xi(2,i)
        daxi2(2,(i-1)*3+2) = -N_xi(3,i)
        daxi1(3,(i-1)*3+3) = -N_xi(2,i)
        daxi2(3,(i-1)*3+3) = -N_xi(3,i)
      ENDDO
      
      ! Secondary nodes (5-8): position matrix and derivatives
      DO i=1,4
        ! Position matrix (x, y, z components)
        a(1,12+(i-1)*3+1) = N_eta(1,i)
        a(2,12+(i-1)*3+2) = N_eta(1,i)
        a(3,12+(i-1)*3+3) = N_eta(1,i)
        
        ! Derivatives w.r.t. eta1 and eta2
        daeta1(1,12+(i-1)*3+1) = N_eta(2,i)
        daeta2(1,12+(i-1)*3+1) = N_eta(3,i)
        daeta1(2,12+(i-1)*3+2) = N_eta(2,i)
        daeta2(2,12+(i-1)*3+2) = N_eta(3,i)
        daeta1(3,12+(i-1)*3+3) = N_eta(2,i)
        daeta2(3,12+(i-1)*3+3) = N_eta(3,i)
      ENDDO

      RETURN
      END
      
C-----------------------------------------------
C Surface geometry calculation
C-----------------------------------------------
C-----------------------------------------------
C Surface geometry calculation
C-----------------------------------------------
      subroutine surfgeom112(xupd, daxi1, daxi2, daeta1, daeta2, 
     &                      norm, rhoxi1, rhoxi2, m_ij, detm, mij, detmPrimary)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     xupd: Element coordinates
C     daxi*, daeta*: Derivative matrices
C     norm: Output normal vector
C     rhoxi*: Output surface tangent vectors
C     m_ij, detm, mij, detmPrimary: Output metric parameters
C-----------------------------------------------
      real*8  xupd(3,8)
      real*8  daxi1(3,24), daxi2(3,24), daeta1(3,24), daeta2(3,24)
      real*8  rhoxi1(3), rhoxi2(3), reta1(3), reta2(3)
      real*8  m_ij(2,2), detm, mij(2,2), mSecondary_ij(2,2), detmPrimary
      real*8  norm(3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER i, j, k
      real*8  inv_sqrt_detm
      real*8  xupd1d(24)
      
C-----------------------------------------------
C   Convert 2D coordinates to 1D for matrix operations
C-----------------------------------------------
      k = 0
      DO i=1,8
        DO j=1,3
          k = k + 1
          xupd1d(k) = xupd(j,i)
        ENDDO
      ENDDO

C-----------------------------------------------
C   Calculate Secondary surface tangent vectors
C-----------------------------------------------
      DO i=1,3
        reta1(i) = 0.d0
        reta2(i) = 0.d0
        DO j=1,24
          reta1(i) = reta1(i) - daeta1(i,j)*xupd1d(j)
          reta2(i) = reta2(i) - daeta2(i,j)*xupd1d(j)
        ENDDO
      ENDDO

C-----------------------------------------------
C   Calculate Secondary surface metric tensor
C-----------------------------------------------
      mSecondary_ij(1,1) = 0.d0
      mSecondary_ij(1,2) = 0.d0
      mSecondary_ij(2,1) = 0.d0
      mSecondary_ij(2,2) = 0.d0

      DO i=1,3
        mSecondary_ij(1,1) = mSecondary_ij(1,1) + reta1(i)**2
        mSecondary_ij(1,2) = mSecondary_ij(1,2) + reta1(i)*reta2(i)
        mSecondary_ij(2,1) = mSecondary_ij(2,1) + reta1(i)*reta2(i)
        mSecondary_ij(2,2) = mSecondary_ij(2,2) + reta2(i)**2
      ENDDO

      detm = mSecondary_ij(1,1)*mSecondary_ij(2,2) - 
     &       mSecondary_ij(2,1)*mSecondary_ij(1,2)

C-----------------------------------------------
C   Calculate Primary surface tangent vectors
C-----------------------------------------------
      DO i=1,3
        rhoxi1(i) = 0.d0
        rhoxi2(i) = 0.d0
        DO j=1,24
          rhoxi1(i) = rhoxi1(i) - daxi1(i,j)*xupd1d(j)
          rhoxi2(i) = rhoxi2(i) - daxi2(i,j)*xupd1d(j)
        ENDDO
      ENDDO

C-----------------------------------------------
C   Calculate Primary surface metric tensor
C-----------------------------------------------
      m_ij(1,1) = 0.d0
      m_ij(1,2) = 0.d0
      m_ij(2,1) = 0.d0
      m_ij(2,2) = 0.d0

      DO i=1,3
        m_ij(1,1) = m_ij(1,1) + rhoxi1(i)**2
        m_ij(1,2) = m_ij(1,2) + rhoxi1(i)*rhoxi2(i)
        m_ij(2,1) = m_ij(2,1) + rhoxi1(i)*rhoxi2(i)
        m_ij(2,2) = m_ij(2,2) + rhoxi2(i)**2
      ENDDO

      detmPrimary = m_ij(1,1)*m_ij(2,2) - m_ij(2,1)*m_ij(1,2)

C-----------------------------------------------
C   Calculate inverse metric tensor
C-----------------------------------------------
      ! Safety check for division by zero
      IF (dabs(detmPrimary) .LT. 1.0d-12) THEN
        ! Use identity matrix if determinant is too small
        mij(1,1) = 1.0d0
        mij(1,2) = 0.d0
        mij(2,1) = 0.d0
        mij(2,2) = 1.0d0
      ELSE
        mij(1,1) = m_ij(2,2) / detmPrimary
        mij(1,2) = -m_ij(1,2) / detmPrimary
        mij(2,1) = -m_ij(2,1) / detmPrimary
        mij(2,2) = m_ij(1,1) / detmPrimary
      ENDIF

C-----------------------------------------------
C   Calculate normal vector (cross product of tangents)
C-----------------------------------------------
      ! Safety check for zero determinant
      IF (dabs(detmPrimary) .LT. 1.0d-12) THEN
        ! Default normal if surface is degenerate
        norm(1) = 0.d0
        norm(2) = 0.d0
        norm(3) = 1.0d0
      ELSE
        inv_sqrt_detm = 1.0d0 / dsqrt(detmPrimary)
        
        ! Normal = (rhoxi1 × rhoxi2) / |rhoxi1 × rhoxi2|
        norm(1) = (rhoxi1(2)*rhoxi2(3) - rhoxi2(2)*rhoxi1(3)) * inv_sqrt_detm
        norm(2) = (rhoxi1(3)*rhoxi2(1) - rhoxi2(3)*rhoxi1(1)) * inv_sqrt_detm
        norm(3) = (rhoxi1(1)*rhoxi2(2) - rhoxi2(1)*rhoxi1(2)) * inv_sqrt_detm
      ENDIF

      RETURN
      END

      
C-----------------------------------------------
C Penetration computation
C-----------------------------------------------
      subroutine penetr112(xupd, penetr, norm, a)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C     xupd: Element coordinates (3,8)
C     penetr: Output penetration distance
C     norm: Normal vector
C     a: Position matrix (3,24)
C-----------------------------------------------
      real*8  xupd(3,8)
      real*8  penetr
      real*8  norm(3)
      real*8  a(3,24)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER i, j, k
      real*8  xupd1d(24)
      
C-----------------------------------------------
C   Convert 2D coordinates to 1D for matrix operations
C-----------------------------------------------
      k = 0
      DO i=1,8
        DO j=1,3
          k = k + 1
          xupd1d(k) = xupd(j,i)
        ENDDO
      ENDDO
      
C-----------------------------------------------
C   Calculate penetration as dot product: penetr = x · a · norm
C-----------------------------------------------
      penetr = 0.d0
      DO i=1,3
        DO j=1,24
          penetr = penetr + xupd1d(j)*a(i,j)*norm(i)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------
C Gaussian quadrature points and weights
C-----------------------------------------------
      subroutine GaussPt(norder, xtab, weight)
C-----------------------------------------------
C   Provides Gauss-Legendre quadrature points and weights
C   for integration on [-1,1]
C
C   Parameters:
C     norder: Number of integration points (1-10)
C     xtab: Output array of quadrature points
C     weight: Output array of quadrature weights
C-----------------------------------------------
      implicit none
      integer norder, i
      real*8 xtab(norder), weight(norder)
      
      real*8, parameter :: ONE = 1.0d0
      real*8, parameter :: SQRT3 = 1.73205080756887729352d0
      
C-----------------------------------------------
C   Initialize arrays
C-----------------------------------------------
      DO i=1,norder
        xtab(i) = 0.0d0
        weight(i) = 0.0d0
      ENDDO
      
C-----------------------------------------------
C   Set quadrature points and weights based on order
C-----------------------------------------------
      IF (norder == 1) THEN
        xtab(1) = 0.0d0
        weight(1) = 2.0d0
        
      ELSE IF (norder == 2) THEN
        xtab(1) = -ONE/SQRT3
        xtab(2) = ONE/SQRT3
        weight(1) = ONE
        weight(2) = ONE
        
      ELSE IF (norder == 3) THEN
        xtab(1) = -0.774596669241483377d0
        xtab(2) = 0.0d0
        xtab(3) = 0.774596669241483377d0
        weight(1) = 5.0d0/9.0d0
        weight(2) = 8.0d0/9.0d0
        weight(3) = 5.0d0/9.0d0
        
      ELSE IF (norder == 4) THEN
        xtab(1) = -0.861136311594052575d0
        xtab(2) = -0.339981043584856264d0
        xtab(3) = 0.339981043584856264d0
        xtab(4) = 0.861136311594052575d0
        weight(1) = 0.347854845137453857d0
        weight(2) = 0.652145154862546142d0
        weight(3) = 0.652145154862546142d0
        weight(4) = 0.347854845137453857d0
        
      ELSE IF (norder == 5) THEN
        xtab(1) = -0.906179845938663992d0
        xtab(2) = -0.538469310105683091d0
        xtab(3) = 0.0d0
        xtab(4) = 0.538469310105683091d0
        xtab(5) = 0.906179845938663992d0
        weight(1) = 0.236926885056189087d0
        weight(2) = 0.478628670499366468d0
        weight(3) = 0.568888888888888888d0
        weight(4) = 0.478628670499366468d0
        weight(5) = 0.236926885056189087d0
        
      ELSE
        ! Unsupported order - fallback to 3-point rule
        IF (norder > 0 .AND. norder <= 10) THEN
          ! Use 3-point rule for unsupported orders
          xtab(1) = -0.774596669241483377d0
          xtab(2) = 0.0d0
          xtab(3) = 0.774596669241483377d0
          weight(1) = 5.0d0/9.0d0
          weight(2) = 8.0d0/9.0d0
          weight(3) = 5.0d0/9.0d0
        ELSE
          ! Invalid order - set to default
          xtab(1) = 0.0d0
          weight(1) = 2.0d0
        END IF
      ENDIF
      
      RETURN
      END
      
C-----------------------------------------------
C Lobatto quadrature points and weights
C-----------------------------------------------
      subroutine LOBATTOPt(norder, xtab, weight)
C-----------------------------------------------
C   Provides Lobatto quadrature points and weights
C   for integration on [-1,1]
C
C   Parameters:
C     norder: Number of integration points (2-10)
C     xtab: Output array of quadrature points
C     weight: Output array of quadrature weights
C-----------------------------------------------
      implicit none
      integer norder, i
      real*8 xtab(norder), weight(norder)

      real*8, parameter :: ONE = 1.0d0

C-----------------------------------------------
C   Initialize arrays
C-----------------------------------------------
      DO i=1,norder
        xtab(i) = 0.0d0
        weight(i) = 0.0d0
      ENDDO
      
C-----------------------------------------------
C   Set quadrature points and weights based on order
C-----------------------------------------------
      IF (norder == 2) THEN
        xtab(1) = -ONE
        xtab(2) = ONE
        weight(1) = ONE
        weight(2) = ONE
        
      ELSE IF (norder == 3) THEN
        xtab(1) = -ONE
        xtab(2) = 0.0d0
        xtab(3) = ONE
        weight(1) = ONE / 3.0d0
        weight(2) = 4.0d0 / 3.0d0
        weight(3) = ONE / 3.0d0
        
      ELSE IF (norder == 4) THEN
        xtab(1) = -ONE
        xtab(2) = -0.447213595499957939d0
        xtab(3) = 0.447213595499957939d0
        xtab(4) = ONE
        weight(1) = ONE / 6.0d0
        weight(2) = 5.0d0 / 6.0d0
        weight(3) = 5.0d0 / 6.0d0
        weight(4) = ONE / 6.0d0
        
      ELSE IF (norder == 5) THEN
        xtab(1) = -ONE
        xtab(2) = -0.654653670707977143d0
        xtab(3) = 0.0d0
        xtab(4) = 0.654653670707977143d0
        xtab(5) = ONE
        weight(1) = ONE / 10.0d0
        weight(2) = 49.0d0 / 90.0d0
        weight(3) = 32.0d0 / 45.0d0
        weight(4) = 49.0d0 / 90.0d0
        weight(5) = ONE / 10.0d0
        
      ELSE
        ! Unsupported order - fallback to 3-point rule
        IF (norder > 0 .AND. norder <= 10) THEN
          xtab(1) = -ONE
          xtab(2) = 0.0d0
          xtab(3) = ONE
          weight(1) = ONE / 3.0d0
          weight(2) = 4.0d0 / 3.0d0
          weight(3) = ONE / 3.0d0
        ELSE
          ! Invalid order - set to default 2-point rule
          xtab(1) = -ONE
          xtab(2) = ONE
          weight(1) = ONE
          weight(2) = ONE
        END IF
      ENDIF
      
      RETURN
      END
      

C-----------------------------------------------
C   Compute convective increments in local coordinates for friction
C-----------------------------------------------
      subroutine tangentvel118(xi1, xi2, dxi1, dxi2,
     &                         EL_NR, z, q, ip,
     &                         XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST, MAX_STS_SIZE)

      implicit none

C-- Input arguments
      real*8  xi1, xi2
      integer EL_NR, z, q, ip, MAX_STS_SIZE

C-- In/out history arrays
      real*8  XI1_HIST(MAX_STS_SIZE,2,2)
      real*8  XI2_HIST(MAX_STS_SIZE,2,2)
      real*8  TTRIAL1_HIST(MAX_STS_SIZE,2,2)
      real*8  TTRIAL2_HIST(MAX_STS_SIZE,2,2)

C-- Output arguments
      real*8  dxi1, dxi2

C-- Local variables
      logical first_contact

      ! Check if EL_NR is out of range
      if (EL_NR .LE. 0 .OR. EL_NR .GT. MAX_STS_SIZE) then
        dxi1 = 0.d0
        dxi2 = 0.d0
        return
      endif
      ! Check if z or q is out of range
      if (z .LE. 0 .OR. z .GT. ip .OR. q .LE. 0 .OR. q .GT. ip) then
        dxi1 = 0.d0
        dxi2 = 0.d0
        return
      endif

C-----------------------------------------------
C  First contact? -> no relative displacement
C-----------------------------------------------
      first_contact = (XI1_HIST(EL_NR,z,q) == 0.d0 .and.
     &                 XI2_HIST(EL_NR,z,q) == 0.d0)

      if (first_contact) then
C       At the very first contact: no relative displacement
        dxi1 = 0.d0
        dxi2 = 0.d0
        TTRIAL1_HIST(EL_NR,z,q) = 0.d0
        TTRIAL2_HIST(EL_NR,z,q) = 0.d0
      else
C       Incremental displacement between current and previous step
        dxi1 = xi1 - XI1_HIST(EL_NR,z,q)
        dxi2 = xi2 - XI2_HIST(EL_NR,z,q)
      endif

C-----------------------------------------------
C  Store current xi values as "previous" for next step
C-----------------------------------------------
      XI1_HIST(EL_NR,z,q) = xi1
      XI2_HIST(EL_NR,z,q) = xi2

      return
      end


C-----------------------------------------------
C   Global convective increments with period tracking
C-----------------------------------------------
      subroutine tangentvel118_global(xi1, xi2, dxi1, dxi2,
     &     EL_NR, z, q, ip, MAX_STS_SIZE_ACTUAL, MAX_GLOBAL_GP,
     &     GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     &     GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     &     GP_XI1_PERIOD, GP_XI2_PERIOD,
     &     GP_INITIALIZED,
     &     GP_TTRIAL1_HIST, GP_TTRIAL2_HIST)

      implicit none

C-- Input arguments
      real*8  xi1, xi2
      integer EL_NR, z, q, ip, MAX_STS_SIZE_ACTUAL, MAX_GLOBAL_GP

C-- Global GP state arrays (1D)
      real*8  GP_XI1_GLOBAL(*)
      real*8  GP_XI2_GLOBAL(*)
      real*8  GP_XI1_GLOBAL_PREV(*)
      real*8  GP_XI2_GLOBAL_PREV(*)
      integer GP_XI1_PERIOD(*)
      integer GP_XI2_PERIOD(*)
      logical GP_INITIALIZED(*)

C-- Global traction history arrays (1D)
      real*8  GP_TTRIAL1_HIST(*)
      real*8  GP_TTRIAL2_HIST(*)

      real*8  dxi1, dxi2

C-- Local variables
      integer gp_index
      real*8  xi1_prev_local, xi2_prev_local
      real*8  tol
      INTEGER GET_GLOBAL_GP_INDEX

      tol = 1.0d-6

      dxi1 = 0.d0
      dxi2 = 0.d0

      if (EL_NR .LE. 0 .OR. EL_NR .GT. MAX_STS_SIZE_ACTUAL) return
      if (z .LE. 0 .OR. z .GT. ip .OR. q .LE. 0 .OR. q .GT. ip) return

      gp_index = GET_GLOBAL_GP_INDEX(EL_NR, z, q, ip)

      if (gp_index .LE. 0 .OR. gp_index .GT. MAX_GLOBAL_GP) then
        ! Invalid index - return zero increments
        return
      endif

      if (.NOT. GP_INITIALIZED(gp_index)) then
C       First contact: initialize global state and reset traction history
        GP_XI1_GLOBAL(gp_index)      = xi1
        GP_XI2_GLOBAL(gp_index)      = xi2
        GP_XI1_GLOBAL_PREV(gp_index) = xi1
        GP_XI2_GLOBAL_PREV(gp_index) = xi2
        GP_XI1_PERIOD(gp_index)      = 0
        GP_XI2_PERIOD(gp_index)      = 0

        GP_TTRIAL1_HIST(gp_index)    = 0.d0
        GP_TTRIAL2_HIST(gp_index)    = 0.d0

        GP_INITIALIZED(gp_index)     = .TRUE.
        dxi1 = 0.d0
        dxi2 = 0.d0
      else
C       Reconstruct previous local coordinates from global + period
        xi1_prev_local = GP_XI1_GLOBAL(gp_index)
     &                  - 2.0d0*GP_XI1_PERIOD(gp_index)
        xi2_prev_local = GP_XI2_GLOBAL(gp_index)
     &                  - 2.0d0*GP_XI2_PERIOD(gp_index)

C       Update period counters based on boundary crossing
        call handle_element_transition(xi1_prev_local, xi1,
     &       xi2_prev_local, xi2,
     &       GP_XI1_PERIOD(gp_index), GP_XI2_PERIOD(gp_index), tol)

C       Update global coordinates
        GP_XI1_GLOBAL_PREV(gp_index) = GP_XI1_GLOBAL(gp_index)
        GP_XI2_GLOBAL_PREV(gp_index) = GP_XI2_GLOBAL(gp_index)

        GP_XI1_GLOBAL(gp_index) = xi1 + 2.0d0*GP_XI1_PERIOD(gp_index)
        GP_XI2_GLOBAL(gp_index) = xi2 + 2.0d0*GP_XI2_PERIOD(gp_index)

        dxi1 = GP_XI1_GLOBAL(gp_index) - GP_XI1_GLOBAL_PREV(gp_index)
        dxi2 = GP_XI2_GLOBAL(gp_index) - GP_XI2_GLOBAL_PREV(gp_index)

C       Check if dxi_global is suspiciously large (indicates lost contact)
C       If dxi > 0.2, treat as new contact and reinitialize
C       Note: Normal dxi should be very small (e.g., 1e-7 to 1e-6 per timestep)
        if (ABS(dxi1) > 0.2d0 .OR. ABS(dxi2) > 0.2d0) then
C         Reinitialize as new contact
          GP_XI1_GLOBAL(gp_index)      = xi1
          GP_XI2_GLOBAL(gp_index)      = xi2
          GP_XI1_GLOBAL_PREV(gp_index) = xi1
          GP_XI2_GLOBAL_PREV(gp_index) = xi2
          GP_XI1_PERIOD(gp_index)      = 0
          GP_XI2_PERIOD(gp_index)      = 0
          GP_TTRIAL1_HIST(gp_index)    = 0.d0
          GP_TTRIAL2_HIST(gp_index)    = 0.d0
          dxi1 = 0.d0
          dxi2 = 0.d0
        endif
      endif  ! End of if (.NOT. GP_INITIALIZED)

      return
      end

C-----------------------------------------------
C   Handle element transition for global xi tracking
C-----------------------------------------------
      subroutine handle_element_transition(xi1_local_prev, xi1_local_new,
     &     xi2_local_prev, xi2_local_new,
     &     gp_xi1_period, gp_xi2_period, tol)

      implicit none

C-- Input arguments
      real*8  xi1_local_prev, xi1_local_new
      real*8  xi2_local_prev, xi2_local_new
      real*8  tol

C-- In/out arguments
      integer gp_xi1_period
      integer gp_xi2_period

C-----------------------------------------------
C  Detect crossings of +/-1 in xi1 and xi2 and
C  update period counters accordingly.
C-----------------------------------------------

C-- Check transition in xi1 direction
      if (dabs(xi1_local_prev - 1.0d0) < tol) then
C       Was at +1 boundary, likely crossed to next element
        if (xi1_local_new < 0.0d0) then
C         Crossed from +1 to -1 (forward)
          gp_xi1_period = gp_xi1_period + 1
        endif
      else if (dabs(xi1_local_prev + 1.0d0) < tol) then
C       Was at -1 boundary, likely crossed to previous element
        if (xi1_local_new > 0.0d0) then
C         Crossed from -1 to +1 (backward)
          gp_xi1_period = gp_xi1_period - 1
        endif
      endif

C-- Check transition in xi2 direction (same logic)
      if (dabs(xi2_local_prev - 1.0d0) < tol) then
        if (xi2_local_new < 0.0d0) then
          gp_xi2_period = gp_xi2_period + 1
        endif
      else if (dabs(xi2_local_prev + 1.0d0) < tol) then
        if (xi2_local_new > 0.0d0) then
          gp_xi2_period = gp_xi2_period - 1
        endif
      endif

      return
      end