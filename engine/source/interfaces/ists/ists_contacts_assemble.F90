!||====================================================================
!||    STS_CONTACTS_ASSEMBLE  ../engine/source/interfaces/ists/ists_contacts_assemble.F90
!||--- called by ------------------------------------------------------
!||    i7mainf              ../engine/source/interfaces/int07/i7mainf.F
!||--- calls ---------------------------------------------------------
!||    STS_CONTACT_EVAL_PAIR    ../engine/source/interfaces/ists/ists_CONTACT_EVAL_PAIR.F90
!||====================================================================
       SUBROUTINE STS_CONTACTS_ASSEMBLE(CONT_ELEMENT, COUNT, OPTION, &
     & CAND_MST_SEG_ID, CAND_SEC_SEG_ID, &
     & load_arr, node_id_load, L_out, IMPACT_glob, STIF, &
     & MAX_STS_SIZE_ACTUAL, FRICC, FRIC_COEFS, VISCFFRIC, XMU, MFROT, &
     & IFQ, DT1, DT12, V, MS, CAND_F, ALPHA0, IFPEN, INTTH, QFRICT, &
     & INTBUF_TAB, GAP, XI1_HIST, XI2_HIST, TTRIAL1_HIST, &
     & TTRIAL2_HIST)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE INTBUFDEF_MOD
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
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, L, IMPACT
      REAL*8 XUPD(3,8)
      REAL*8 p_load_new(24)
      REAL*8 node_stiff(8)
      REAL*8 p_friction(24)  ! Friction forces (separate output)
      my_real EFRICT_LOC
      INTEGER node_ids(8)  ! Node IDs for velocity interpolation
!-----------------------------------------------
!   I n i t i a l i z a t i o n
!-----------------------------------------------
      IMPACT_glob = 0
      
      ! Safety check
      IF (COUNT <= 0) THEN
        L_out = 1
        RETURN
      END IF
      
      ! Initialize counters
      K = 1
      L = 1
!-----------------------------------------------
!   M a i n   L o o p
!-----------------------------------------------
      DO I = 1, COUNT
        IMPACT = 0
        XUPD = CONT_ELEMENT(I, 1:3, 1:8)
      
        ! Get node IDs for velocity interpolation
        DO J = 1, 4
          node_ids(J)   = CAND_MST_SEG_ID(I, J+1)   ! Primary nodes
          node_ids(J+4) = CAND_SEC_SEG_ID(I, J+1)   ! Secondary nodes
        ENDDO
        XMU(1) = FRICC(1) ! Friction coefficient mu
      
        ! Call STS_CONTACT_EVAL_PAIRwith friction calculation integrated
        ! Note: STIF is used for both normal and tangential stiffness
        ! (STIF0 = STIF for now, until separate tangential stiffness is available)
        CALL STS_CONTACT_EVAL_PAIR(XUPD, STIF, p_load_new, IMPACT, I, &
     &                    node_stiff, OPTION, &
     &                    V, MS, FRICC, FRIC_COEFS, VISCFFRIC, XMU, &
     &                    MFROT, IFQ, ALPHA0, CAND_F, IFPEN, STIF, &
     &                    p_friction, EFRICT_LOC, QFRICT, &
     &                    INTTH, node_ids, .TRUE., &
     &                    XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST, &
     &                    MAX_STS_SIZE_ACTUAL, GAP, &
     &                    CAND_SEC_SEG_ID)
      
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
          
          ! Safety check - prevent array overflow
          IF (L > MAX_STS_SIZE_ACTUAL .OR. K > MAX_STS_SIZE_ACTUAL*8) THEN
            EXIT
          END IF
        ENDIF
      ENDDO
      
      L_out = L
      END SUBROUTINE STS_CONTACTS_ASSEMBLE
      
!************************************************************************
!----------------------------------------------------------------------!
!.... general subroutines
!----------------------------------------------------------------------!
!************************************************************************
      
      INTEGER FUNCTION GET_GLOBAL_GP_INDEX(PAIR_INDEX, Z, Q, IP_MAX)
!-----------------------------------------------------------------------
! Compute global Gauss point index from STS pair index and
! local Gauss point indices (z, q) for a given quadrature order IP_MAX.
!
! gp_index = (pair_index-1) * (IP_MAX*IP_MAX) + local_gp_index + 1
! local_gp_index = (z-1)*IP_MAX + (q-1)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER PAIR_INDEX, Z, Q, IP_MAX
      INTEGER LOCAL_GP_INDEX

      LOCAL_GP_INDEX = (Z - 1) * IP_MAX + (Q - 1)
      GET_GLOBAL_GP_INDEX = (PAIR_INDEX - 1) * (IP_MAX*IP_MAX) + LOCAL_GP_INDEX + 1

      RETURN
      END