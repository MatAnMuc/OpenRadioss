!||====================================================================
!||    STS_CONTACT_EVAL_PAIR   ../engine/source/interfaces/ists/ists_contact_eval_pair.F90
!||--- called by ------------------------------------------------------
!||    STS_CONTACTS_ASSEMBLE   ../engine/source/interfaces/ists/ists_contacts_assemble.F90
!||--- calls ---------------------------------------------------------
!||    sts_gausspt             ../engine/source/interfaces/ists/ists_sts_gausspt.F90
!||    sts_lobattopt           ../engine/source/interfaces/ists/ists_sts_lobattopt.F90
!||    sts_project             ../engine/source/interfaces/ists/ists_sts_project.F90
!||    sts_pos                 ../engine/source/interfaces/ists/ists_pos.F90
!||    sts_surfgeom            ../engine/source/interfaces/ists/ists_sts_surfgeom.F90
!||    sts_penetr              ../engine/source/interfaces/ists/ists_sts_penetr.F90
!||    sts_tangentvel_global   ../engine/source/interfaces/ists/ists_tangentvel.F90
!||    sts_shape               ../engine/source/interfaces/ists/ists_shape_fct.F90
!||====================================================================
      subroutine STS_CONTACT_EVAL_PAIR(XUPD, STIF, p, IMPACT, EL_NR, node_stiff, OPTION, &
     &                   V, MS, FRICC, FRIC_COEFS, VISCFFRIC, XMU, MFROT, &
     &                   IFQ, ALPHA0, CAND_F, IFPEN, STIF0, &
     &                   p_friction, EFRICT, QFRICT, INTTH, node_ids, &
     &                   CALC_FRICTION, XI1_HIST, XI2_HIST, &
     &                   TTRIAL1_HIST, TTRIAL2_HIST, MAX_STS_SIZE, NUMNOD, GAP, &
     &                   CAND_SEC_SEG_ID)
!-----------------------------------------------
!   M o d u l e s   /   I m p l i c i t   T y p e s
!-----------------------------------------------
      use constant_mod
      use sts_gp_state_mod
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!     XUPD     : Coordinates of contact element (3,8)
!                Nodes 1-4 are Primary (master) nodes
!                Nodes 5-8 are Secondary (slave) nodes
!     STIF     : Contact stiffness parameter (normal)
!     p        : Output contact forces (24 components) - normal + friction combined
!     IMPACT   : Output flag indicating if penetration was detected
!     EL_NR    : Temporary Element number (for debugging)
!     node_stiff: Output nodal stiffness values (8 components)
!     OPTION   : 0 = Gauss quadrature, 1 = Lobatto quadrature
!     V        : Nodal velocities (3,:)
!     MS       : Nodal masses (:)
!     FRICC, FRIC_COEFS, VISCFFRIC, XMU: Friction parameters
!     MFROT    : Friction model type
!     IFQ      : Friction formulation flag
!     ALPHA0   : Friction parameter
!     CAND_F   : Stored friction forces (for incremental formulation)
!     IFPEN    : Penetration flag array
!     STIF0    : Tangential contact stiffness
!     p_friction: Output friction forces (24 components) - separate output
!     EFRICT   : Friction energy (output)
!     QFRICT   : Total friction energy (INTENT(INOUT))
!     INTTH    : Thermal interface flag
!     node_ids : Node IDs for velocity interpolation (8 components)
!     CALC_FRICTION: Flag to enable/disable friction calculation
!     XI1_HIST : History of xi1 parametric coordinates per Gauss point
!     XI2_HIST : History of xi2 parametric coordinates per Gauss point
!     MAX_STS_SIZE: Maximum size for history arrays
!     GAPV: Gap value from user input
!-----------------------------------------------
      INTEGER IMPACT, OPTION, EL_NR, MFROT, IFQ, INTTH
      INTEGER MAX_STS_SIZE, NUMNOD  ! Maximum size for history arrays / number of nodes
      INTEGER node_ids(8)
      LOGICAL CALC_FRICTION
      my_real STIF(MVSIZ), STIF0(MVSIZ)
      real*8  p(24), p_friction(24)
      real*8  XUPD(3,8)
      real*8  node_stiff(8)
      my_real V(3,NUMNOD), MS(NUMNOD)
      my_real CAND_F(8,MAX_STS_SIZE)
      INTEGER IFPEN(MAX_STS_SIZE)
      my_real FRICC(MVSIZ), FRIC_COEFS(MVSIZ,10), VISCFFRIC(MVSIZ)
      my_real XMU(MVSIZ), ALPHA0
      my_real EFRICT, QFRICT
      REAL*8 XI1_HIST(MAX_STS_SIZE,2,2)  ! History of xi1 per Gauss point
      REAL*8 XI2_HIST(MAX_STS_SIZE,2,2)  ! History of xi2 per Gauss point
      REAL*8 TTRIAL1_HIST(MAX_STS_SIZE,2,2)  ! History of T_trial(1) per Gauss point
      REAL*8 TTRIAL2_HIST(MAX_STS_SIZE,2,2)  ! History of T_trial(2) per Gauss poin
      my_real GAP  ! Gap value from user input
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE,5)
!     interface to global gp index function
      INTEGER GET_GLOBAL_GP_INDEX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
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
!-----------------------------------------------
!   I n i t i a l i z a t i o n
!-----------------------------------------------
      IMPACT = 0
      ip = 2 ! Quadrature order
      
      ! Calculate maximum global GP index for bounds checking
      MAX_GLOBAL_GP_LOCAL = MAX_STS_SIZE * ip * ip
      
      ! Get quadrature points and weights
      
      ! Always initialize Gauss quadrature for friction calculation
      call sts_gausspt(ip, eta1_gauss, wi1_gauss)
      call sts_gausspt(ip, eta2_gauss, wi2_gauss)    
      
      ! Initialize quadrature points and weights for normal contact calculation
      IF (OPTION == 0) THEN
        ! Gauss: reuse Gauss quadrature points and weights
        eta1 = eta1_gauss
        eta2 = eta2_gauss
        wi1 = wi1_gauss
        wi2 = wi2_gauss
      ELSE
        ! Lobatto: use Lobatto quadrature points and weights
        call sts_lobattopt(ip, eta1, wi1)
        call sts_lobattopt(ip, eta2, wi2)
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
!-----------------------------------------------
!   M a i n   C o m p u t a t i o n
!-----------------------------------------------
!     Loop over integration points
      DO z=1,ip
        DO q=1,ip
          
          ! Project Secondary surface to Primary surface at current Gauss/Lobatto point
          call sts_project(XUPD, xi1, xi2, eta1(z), eta2(q))
          
          ! Check if projection is valid
          IF ((dabs(xi1) .GT. 1.05d0) .OR. (dabs(xi2) .GT. 1.05d0)) THEN
            CYCLE
          ENDIF
          
          ! Build position and derivative matrices
          call sts_pos(a, daxi1, daxi2, daeta1, daeta2, xi1, xi2, &
     &               eta1(z), eta2(q))
          
          ! Calculate surface geometry and metrics
          call sts_surfgeom(XUPD, daxi1, daxi2, daeta1, daeta2, norm, &
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
            call sts_project(XUPD, xi1_gauss, xi2_gauss, &
     &                     eta1_gauss(z), eta2_gauss(q))
            
            ! Check if Gauss projection is valid for friction calculation
            gauss_valid = (dabs(xi1_gauss) .LE. 1.05d0 .AND. &
     &                     dabs(xi2_gauss) .LE. 1.05d0)
            
            IF (gauss_valid) THEN
              ! Calculate surface geometry for Gauss projection
              call sts_pos(a_gauss, daxi1_gauss, daxi2_gauss, &
     &                   daeta1_gauss, daeta2_gauss, &
     &                   xi1_gauss, xi2_gauss, &
     &                   eta1_gauss(z), eta2_gauss(q))
              
              call sts_surfgeom(XUPD, daxi1_gauss, daxi2_gauss, &
     &                        daeta1_gauss, daeta2_gauss, &
     &                        norm, rhoxi1_gauss, rhoxi2_gauss, &
     &                        m_ij_gauss, detm_gauss, mij, detmPrimary)
            ENDIF
          ENDIF
          
          ! ==== Normal impact =====
          ! Compute penetration at Gauss point
          call sts_penetr(XUPD, penetr, norm, a)
          
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
          energy = energy + 0.5d0 * d1 * penetr**2 * wi1(z) * wi2(q) * &
     &            dsqrt(detm)
          
          ! Compute residual forces (normal component)
          DO i=1,24
            DO j=1,3
              pm(i) = pm(i) + d1 * penetr * &
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
            call sts_tangentvel_global(xi1_fric, xi2_fric, dxi1, dxi2, &
     &          EL_NR, z, q, ip, MAX_STS_SIZE)

            ! Calculate tangential traction using selected metric tensor
            T_trial(1) = GP_TTRIAL1_HIST(gp_index) - &
     &                   d1*(m_ij_fric(1,1)*dxi1 + m_ij_fric(1,2)*dxi2)
            T_trial(2) = GP_TTRIAL2_HIST(gp_index) - &
     &                   d1*(m_ij_fric(2,1)*dxi1 + m_ij_fric(2,2)*dxi2)

            ! Magnitude of trial tangential traction using selected metric
            T_trialabs = 0.0d0
            DO i = 1, 2
              DO j = 1, 2
                T_trialabs = T_trialabs + T_trial(i) * T_trial(j) * &
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
                  WRITE(6,1000) EL_NR, node_ids(j), z, q, &
     &                 FRICTION_QUAD_TYPE, &
     &                 xi1_fric, xi2_fric, dxi1, dxi2, &
     &                 GP_XI1_PERIOD(gp_index), &
     &                 GP_XI2_PERIOD(gp_index), &
     &                 FXT, FYT, FZT, &
     &                 STICK_STATE
                ENDIF
              ENDDO
            ENDIF

            ! ===== ACCUMULATE FRICTION FORCES TO NODES =====
            call sts_shape(xi1_fric, xi2_fric, N_xi)
            call sts_shape(eta1_fric, eta2_fric, N_eta)
            ! Primary nodes (1-4): subtract friction
            DO j=1,4
              pm_friction((j-1)*3+1) = pm_friction((j-1)*3+1) - &
     &                               N_xi(1,j) * FXT * wi1_fric * &
     &                               wi2_fric * dsqrt(detm_fric)
              pm_friction((j-1)*3+2) = pm_friction((j-1)*3+2) - &
     &                               N_xi(1,j) * FYT * wi1_fric * &
     &                               wi2_fric * dsqrt(detm_fric)
              pm_friction((j-1)*3+3) = pm_friction((j-1)*3+3) - &
     &                               N_xi(1,j) * FZT * wi1_fric * &
     &                               wi2_fric * dsqrt(detm_fric)
            ENDDO
            
            ! Secondary nodes (5-8): add friction
            DO j=1,4
              pm_friction(12+(j-1)*3+1) = pm_friction(12+(j-1)*3+1) + &
     &                                  N_eta(1,j) * FXT * wi1_fric * &
     &                                  wi2_fric * dsqrt(detm_fric)
              pm_friction(12+(j-1)*3+2) = pm_friction(12+(j-1)*3+2) + &
     &                                  N_eta(1,j) * FYT * wi1_fric * &
     &                                  wi2_fric * dsqrt(detm_fric)
              pm_friction(12+(j-1)*3+3) = pm_friction(12+(j-1)*3+3) + &
     &                                  N_eta(1,j) * FZT * wi1_fric * &
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
         
!-----------------------------------------------
!   F i n a l   R e s u l t s
!-----------------------------------------------
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

 1000 FORMAT('STS-DBG EL=',I6,' NODE=',I10,' z=',I2,' q=',I2, &
     &       ' quad=',A7, &
     &       ' xi=',2F12.6, &
     &       ' dxi=',2F12.6, &
     &       ' per=',2I6, &
     &       ' FT=',3E12.4, &
     &       ' state=',A7)

      RETURN
      END