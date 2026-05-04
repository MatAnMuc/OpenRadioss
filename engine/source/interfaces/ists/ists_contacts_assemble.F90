!||====================================================================
!||    STS_CONTACTS_ASSEMBLE  ../engine/source/interfaces/ists/ists_contacts_assemble.F90
!||--- called by ------------------------------------------------------
!||    i7mainf              ../engine/source/interfaces/int07/i7mainf.F
!||--- calls ---------------------------------------------------------
!||    STS_CONTACT_EVAL_PAIR    ../engine/source/interfaces/ists/ists_CONTACT_EVAL_PAIR.F90
!||====================================================================
      SUBROUTINE STS_CONTACTS_ASSEMBLE(CONT_ELEMENT, COUNT, OPTION, STS_INTERFACE_ID, NCYCLE_IN, TIME_CUR, &
     & CAND_MST_SEG_ID, CAND_SEC_SEG_ID, &
     & load_arr, node_id_load, L_out, IMPACT_glob, STIF, &
     & MAX_STS_SIZE_ACTUAL, FRICC, XMU, IFPEN, QFRICT, GAP)
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
!      TYPE(INTBUF_STRUCT_) INTBUF_TAB(*)
      REAL*8 CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8)
      my_real STIF(MVSIZ)
      INTEGER COUNT, OPTION, STS_INTERFACE_ID, NCYCLE_IN
      REAL*8 TIME_CUR
      INTEGER CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      INTEGER CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5)
      REAL*8 load_arr(MAX_STS_SIZE_ACTUAL,8,4)
      INTEGER node_id_load(MAX_STS_SIZE_ACTUAL*8)
      INTEGER L_out, IMPACT_glob, MAX_STS_SIZE_ACTUAL
      my_real FRICC(MVSIZ)
      my_real XMU(MVSIZ)
      INTEGER IFPEN(MAX_STS_SIZE_ACTUAL)     
      my_real QFRICT
      my_real GAP  ! Gap value from user input
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, L, IMPACT
      INTEGER LUX_STS
      REAL*8 XUPD(3,8)
      REAL*8 p_load_new(24)
      REAL*8 node_stiff(8)
      REAL*8 p_friction(24)  ! Friction forces (separate output)
      REAL*8 pair_max_penetration
      my_real EFRICT_LOC
      INTEGER node_ids(8)  ! Node IDs for velocity interpolation
      REAL*8 fx_prim, fy_prim, fz_prim, fx_sec, fy_sec, fz_sec
      REAL*8 fxf_prim, fyf_prim, fzf_prim, fxf_sec, fyf_sec, fzf_sec
      REAL*8 k_pair_avg, fn_pair_est
      LOGICAL FILE_EXISTS_STS, STS_CSV_INITIALIZED
      LOGICAL, PARAMETER :: CSV_OUTPUT_ENABLED = .FALSE.
      LOGICAL, PARAMETER :: STS_DEBUG_PRINT_FORCES = .FALSE.
      SAVE STS_CSV_INITIALIZED
      DATA STS_CSV_INITIALIZED /.FALSE./
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

      IF (CSV_OUTPUT_ENABLED) THEN
        IF (.NOT. STS_CSV_INITIALIZED) THEN
          INQUIRE(FILE='sts_contact_forces.csv', EXIST=FILE_EXISTS_STS)
          IF (FILE_EXISTS_STS) THEN
            OPEN(NEWUNIT=LUX_STS, FILE='sts_contact_forces.csv', &
     &           STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
          ELSE
            OPEN(NEWUNIT=LUX_STS, FILE='sts_contact_forces.csv', &
     &           STATUS='NEW', ACTION='WRITE')
            WRITE(LUX_STS,'(A)') &
     &        'cycle,time,interface_id,entity_id,fx,fy,fz,force_norm,fn,ft,n_pairs,max_penetration'
          ENDIF
          STS_CSV_INITIALIZED = .TRUE.
        ELSE
          OPEN(NEWUNIT=LUX_STS, FILE='sts_contact_forces.csv', &
     &         STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        ENDIF
      END IF
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
        XMU(1) = FRICC(MIN(I,MVSIZ)) ! Friction coefficient mu
      
        ! Call STS_CONTACT_EVAL_PAIRwith friction calculation integrated
        ! Note: STIF is used for both normal and tangential stiffness
        ! (STIF0 = STIF for now, until separate tangential stiffness is available)
        CALL STS_CONTACT_EVAL_PAIR(XUPD, STIF, p_load_new, IMPACT, I, &
     &                    node_stiff, OPTION, &
     &                    FRICC, XMU, IFPEN, &
     &                    p_friction, EFRICT_LOC, QFRICT, node_ids, &
     &                    .TRUE., MAX_STS_SIZE_ACTUAL, GAP, &
     &                    pair_max_penetration)
      
        IF (IMPACT == 1) THEN
          IMPACT_glob = 1

!         Export two rows per contact pair:
!         - surface 1 (primary nodes 1..4)
!         - surface 2 (secondary nodes 5..8)
          fx_prim = 0.0D0
          fy_prim = 0.0D0
          fz_prim = 0.0D0
          fxf_prim = 0.0D0
          fyf_prim = 0.0D0
          fzf_prim = 0.0D0
          DO J = 1, 4
            fx_prim = fx_prim + p_load_new(3*(J-1)+1)
            fy_prim = fy_prim + p_load_new(3*(J-1)+2)
            fz_prim = fz_prim + p_load_new(3*(J-1)+3)
            fxf_prim = fxf_prim + p_friction(3*(J-1)+1)
            fyf_prim = fyf_prim + p_friction(3*(J-1)+2)
            fzf_prim = fzf_prim + p_friction(3*(J-1)+3)
          ENDDO

          fx_sec = 0.0D0
          fy_sec = 0.0D0
          fz_sec = 0.0D0
          fxf_sec = 0.0D0
          fyf_sec = 0.0D0
          fzf_sec = 0.0D0
          DO J = 5, 8
            fx_sec = fx_sec + p_load_new(3*(J-1)+1)
            fy_sec = fy_sec + p_load_new(3*(J-1)+2)
            fz_sec = fz_sec + p_load_new(3*(J-1)+3)
            fxf_sec = fxf_sec + p_friction(3*(J-1)+1)
            fyf_sec = fyf_sec + p_friction(3*(J-1)+2)
            fzf_sec = fzf_sec + p_friction(3*(J-1)+3)
          ENDDO

          IF (CSV_OUTPUT_ENABLED) THEN
            CALL STS_CONTACT_EXPORT_CSV_PAIR(LUX_STS, NCYCLE_IN, TIME_CUR, STS_INTERFACE_ID, &
     &          CAND_MST_SEG_ID(I,1), CAND_SEC_SEG_ID(I,1), pair_max_penetration, &
     &          fx_prim, fy_prim, fz_prim, fxf_prim, fyf_prim, fzf_prim, &
     &          fx_sec, fy_sec, fz_sec, fxf_sec, fyf_sec, fzf_sec)
          END IF
      
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
      IF (CSV_OUTPUT_ENABLED) THEN
        CLOSE(LUX_STS)
      END IF
      
      L_out = L
      END SUBROUTINE STS_CONTACTS_ASSEMBLE
      
      SUBROUTINE STS_CONTACT_EXPORT_CSV_PAIR(LUX_STS, NCYCLE_IN, TIME_CUR, STS_INTERFACE_ID, &
     & MST_ENTITY_ID, SEC_ENTITY_ID, pair_max_penetration, &
     & FX_PRIM, FY_PRIM, FZ_PRIM, FXF_PRIM, FYF_PRIM, FZF_PRIM, &
     & FX_SEC, FY_SEC, FZ_SEC, FXF_SEC, FYF_SEC, FZF_SEC)
      IMPLICIT NONE
      INTEGER LUX_STS, NCYCLE_IN, STS_INTERFACE_ID, MST_ENTITY_ID, SEC_ENTITY_ID
      REAL*8 TIME_CUR, pair_max_penetration
      REAL*8 FX_PRIM, FY_PRIM, FZ_PRIM, FXF_PRIM, FYF_PRIM, FZF_PRIM
      REAL*8 FX_SEC, FY_SEC, FZ_SEC, FXF_SEC, FYF_SEC, FZF_SEC
      REAL*8 force_norm, fn_mag, ft_mag

      force_norm = SQRT(FX_PRIM*FX_PRIM + FY_PRIM*FY_PRIM + FZ_PRIM*FZ_PRIM)
      ft_mag = SQRT(FXF_PRIM*FXF_PRIM + FYF_PRIM*FYF_PRIM + FZF_PRIM*FZF_PRIM)
      fn_mag = SQRT(MAX(0.0D0, (FX_PRIM-FXF_PRIM)*(FX_PRIM-FXF_PRIM) + &
     &    (FY_PRIM-FYF_PRIM)*(FY_PRIM-FYF_PRIM) + (FZ_PRIM-FZF_PRIM)*(FZ_PRIM-FZF_PRIM)))
      WRITE(LUX_STS,'(I0,'','',ES23.15,'','',I0,'','',I0)',ADVANCE='NO') &
     &  NCYCLE_IN, TIME_CUR, ABS(STS_INTERFACE_ID), MST_ENTITY_ID
      WRITE(LUX_STS,901) FX_PRIM, FY_PRIM, FZ_PRIM, force_norm, fn_mag, ft_mag, 1, &
     &    pair_max_penetration

      force_norm = SQRT(FX_SEC*FX_SEC + FY_SEC*FY_SEC + FZ_SEC*FZ_SEC)
      ft_mag = SQRT(FXF_SEC*FXF_SEC + FYF_SEC*FYF_SEC + FZF_SEC*FZF_SEC)
      fn_mag = SQRT(MAX(0.0D0, (FX_SEC-FXF_SEC)*(FX_SEC-FXF_SEC) + &
     &    (FY_SEC-FYF_SEC)*(FY_SEC-FYF_SEC) + (FZ_SEC-FZF_SEC)*(FZ_SEC-FZF_SEC)))
      WRITE(LUX_STS,'(I0,'','',ES23.15,'','',I0,'','',I0)',ADVANCE='NO') &
     &  NCYCLE_IN, TIME_CUR, -ABS(STS_INTERFACE_ID), SEC_ENTITY_ID
      WRITE(LUX_STS,901) FX_SEC, FY_SEC, FZ_SEC, force_norm, fn_mag, ft_mag, 1, &
     &    pair_max_penetration

 901  FORMAT ( ',', ES23.15, ',', ES23.15, ',', ES23.15, ',', ES23.15, ',', &
     &    ES23.15, ',', ES23.15, ',', I0, ',', ES23.15)
      END SUBROUTINE STS_CONTACT_EXPORT_CSV_PAIR
      
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
