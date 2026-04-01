!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
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
!Copyright>        Commercial Alternative: Altair Radioss
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!||====================================================================
!||    q1np_forc3                       ../engine/source/elements/solid/solid_q1np/q1np_forc3.F90
!||--- called by ------------------------------------------------------
!||    forint                           ../engine/source/elements/forint.F
!||--- calls      -----------------------------------------------------
!||    ig3daverage                      ../engine/source/elements/ige3d/ig3daverage.F
!||    ige3dbilan                       ../engine/source/elements/ige3d/ige3dbilan.F
!||    mmain                            ../engine/source/materials/mat_share/mmain.F90
!||    q1np_get_knot_vectors            ../common_source/modules/q1np_geom_mod.F90
!||    q1np_init_gauss_scheme_starter   ../common_source/modules/q1np_restart_mod.F90
!||    q1np_jacobian                    ../common_source/modules/q1np_geom_mod.F90
!||    q1np_shape_functions             ../common_source/modules/q1np_geom_mod.F90
!||    smalla3                          ../engine/source/elements/solid/solide/smalla3.F
!||    smallb3                          ../engine/source/elements/solid/solide/smallb3.F
!||    srho3                            ../engine/source/elements/solid/solide/srho3.F
!||    srota3                           ../engine/source/elements/solid/solide/srota3.F
!||    sstra3                           ../engine/source/elements/solid/solide/sstra3.F
!||--- uses       -----------------------------------------------------
!||    ale_connectivity_mod             ../common_source/modules/ale/ale_connectivity_mod.F
!||    com08_mod                        ../engine/share/modules/com08_mod.F
!||    constant_mod                     ../common_source/modules/constant_mod.F
!||    debug_mod                        ../engine/share/modules/debug_mod.F
!||    dt_mod                           ../engine/source/modules/dt_mod.F
!||    elbufdef_mod                     ../common_source/modules/mat_elem/elbufdef_mod.F90
!||    element_mod                      ../common_source/modules/elements/element_mod.F90
!||    glob_therm_mod                   ../common_source/modules/mat_elem/glob_therm_mod.F90
!||    mat_elem_mod                     ../common_source/modules/mat_elem/mat_elem_mod.F90
!||    mmain_mod                        ../engine/source/materials/mat_share/mmain.F90
!||    nlocal_reg_mod                   ../common_source/modules/nlocal_reg_mod.F
!||    output_mod                       ../common_source/modules/output/output_mod.F90
!||    param_c_mod                      ../engine/share/modules/param_c_mod.F
!||    q1np_geom_mod                    ../common_source/modules/q1np_geom_mod.F90
!||    q1np_restart_mod                 ../common_source/modules/q1np_restart_mod.F90
!||    restmod                          ../engine/share/modules/restart_mod.F
!||    sensor_mod                       ../common_source/modules/sensor_mod.F90
!||    table_mod                        ../engine/share/modules/table_mod.F
!||    timer_mod                        ../engine/source/system/timer_mod.F90
!||====================================================================
! Calculation of the internal forces for Q1NP elements
!=======================================================================
      SUBROUTINE Q1NP_FORC3(TIMERS, OUTPUT, &
     &                      ELBUF_STR, PM, GEO, IXS, IGEO, X, A, V, W, &
     &                      FV, ALE_CONNECT, IPARG, TF, NPF, BUFMAT, PARTSAV, NPART, &
     &                      NLOC_DMG, STIFN, OFFSET, IPARTS, NEL, DT2T, NELTST, &
     &                      ITYPTST, IPM, ITASK, GRESAV, GRTH, IGRTH, MSSA, DMELS, &
     &                      TABLE, IPRI, MAT_ELEM, NG, H3D_STRAIN, SVIS, GLOB_THERM, &
     &                      SNPC, NUMGEO, NUMNOD, NUMELS, NUMELQ, NGROUP, SBUFMAT, STF, &
     &                      NUMMAT, NTABLE, NSVOIS, IRESP, IDEL7NOK, MAXFUNC, USERL_AVAIL, &
     &                      IMON_MAT, IMPL_S, IDYNA, IDTMIN, DT, SENSORS)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
      USE TIMER_MOD
      USE OUTPUT_MOD, ONLY : OUTPUT_
      USE ALE_CONNECTIVITY_MOD
      USE ELBUFDEF_MOD
      USE MAT_ELEM_MOD
      USE MMAIN_MOD
      USE NLOCAL_REG_MOD
      USE TABLE_MOD
      USE DT_MOD
      USE GLOB_THERM_MOD
      USE CONSTANT_MOD, ONLY : ZERO, ONE, HALF
      USE DEBUG_MOD,    ONLY : ITAB_DEBUG
      USE RESTMOD,      ONLY : IQ1NP_TAB, IQ1NP_BULK_TAB, KQ1NP_TAB, Q1NP_KTAB
      USE Q1NP_RESTART_MOD
      USE Q1NP_GEOM_MOD
      USE SENSOR_MOD
      USE PARAM_C_MOD
      USE COM08_MOD
      USE ELEMENT_MOD, ONLY : NIXS

      IMPLICIT NONE
  !=======================================================================
  !  Q1NP_FORC3: internal forces calculation for Q1NP solids
  !
  !  Sequential hierarchy 
  !    1) Identify Q1NP element IDs for the group and derive degrees (P,Q)
  !    2) Reconstruct Q1NP parametric grid + initialize Gauss scheme
  !    3) Allocate and initialize vectorized work arrays
  !    4) (TT=0) snapshot Gauss volumes (VOL -> VOL0DP) for restart/diagnostics
  !    5) Build node/group bookkeeping and gather X/V into element-local arrays
  !    6) Initialize element fields (OFF, RHO0, DELTAX, buffers)
  !    7) SMALLA3: small-rotation coordinate transform for small strain storage
  !    8) Gauss integration loop (IU,IV,IT):
  !         - reset per-GP fields
  !         - geometry (DN, detJ, VOLN/VOLG, etc.)
  !         - material/strain update through engine pipeline
  !         - accumulate nodal internal forces in element-local buffers
  !    9) SMALLB3: update OFF/GBUF%OFF scaling after integration
  !   10) Assemble element-local forces into global A + update dt + stress averages
  !=======================================================================

  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "my_real.inc"
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
      TYPE(TIMER_),  INTENT(INOUT) :: TIMERS
      TYPE(OUTPUT_), INTENT(INOUT) :: OUTPUT
      TYPE(ELBUF_STRUCT_), TARGET, INTENT(INOUT) :: ELBUF_STR
      INTEGER,            INTENT(IN)    :: NG
      INTEGER,            INTENT(IN)    :: H3D_STRAIN
      INTEGER,            INTENT(IN)    :: SNPC
      INTEGER,            INTENT(IN)    :: NUMGEO
      INTEGER,            INTENT(IN)    :: NUMNOD
      INTEGER,            INTENT(IN)    :: NUMELS
      INTEGER,            INTENT(IN)    :: NUMELQ
      INTEGER,            INTENT(IN)    :: NGROUP
      INTEGER,            INTENT(IN)    :: SBUFMAT
      INTEGER,            INTENT(IN)    :: STF
      INTEGER,            INTENT(IN)    :: NUMMAT
      INTEGER,            INTENT(IN)    :: NTABLE
      INTEGER,            INTENT(IN)    :: NSVOIS
      INTEGER,            INTENT(IN)    :: IRESP
      INTEGER,            INTENT(IN)    :: MAXFUNC
      INTEGER,            INTENT(IN)    :: USERL_AVAIL
      INTEGER,            INTENT(IN)    :: IMON_MAT
      INTEGER,            INTENT(IN)    :: IMPL_S
      INTEGER,            INTENT(IN)    :: IDYNA
      INTEGER,            INTENT(IN)    :: IDTMIN(102)
      INTEGER,            INTENT(IN)    :: NEL, ITASK, NPART
      INTEGER,            INTENT(INOUT) :: NELTST, ITYPTST
      INTEGER,            INTENT(IN)    :: OFFSET, IPRI
      INTEGER,            INTENT(INOUT) :: IDEL7NOK
      INTEGER,            INTENT(IN)    :: NPF(SNPC)
      INTEGER,            INTENT(IN)    :: IPARTS(*) !SIZE NEL?
      INTEGER,            INTENT(IN)    :: IPM(NPROPMI,NUMMAT)
      INTEGER,            INTENT(IN)    :: IXS(NIXS,NUMELS)
      INTEGER,            INTENT(IN)    :: IPARG(NPARG,NGROUP)
      INTEGER,            INTENT(IN)    :: IGEO(NPROPGI,NUMGEO)
      INTEGER,            INTENT(INOUT) :: IGRTH(*)
      my_real,            INTENT(INOUT) :: PM(NPROPM,NUMMAT)
      my_real,            INTENT(INOUT) :: GEO(NPROPG,NUMGEO)
      my_real,            INTENT(IN)    :: X(3,NUMNOD)
      my_real,            INTENT(INOUT) :: A(3,NUMNOD)
      my_real,            INTENT(INOUT) :: V(3,NUMNOD)
      my_real,            INTENT(IN)    :: W(3,NUMNOD)
      my_real,            INTENT(INOUT) :: STIFN(NUMNOD)
      my_real,            INTENT(INOUT) :: TF(STF)
      my_real,            INTENT(INOUT) :: BUFMAT(SBUFMAT)
      my_real,            INTENT(INOUT) :: FV(*) !SIZE MAXFUNC?
      my_real,            INTENT(INOUT) :: PARTSAV(NPSAV,NPART)
      my_real,            INTENT(INOUT) :: GRESAV(*)
      my_real,            INTENT(INOUT) :: GRTH(*)
      my_real,            INTENT(INOUT) :: DT2T
      my_real,            INTENT(INOUT) :: MSSA(*)
      my_real,            INTENT(INOUT) :: DMELS(*)
      TYPE(TTABLE),       INTENT(INOUT) :: TABLE(*)
      TYPE(MAT_ELEM_),    INTENT(INOUT) :: MAT_ELEM
      TYPE(NLOCAL_STR_),  INTENT(INOUT) :: NLOC_DMG
      TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
      my_real, DIMENSION(MVSIZ,6), INTENT(INOUT) :: SVIS
      TYPE(GLOB_THERM_),  INTENT(INOUT) :: GLOB_THERM
      TYPE(DT_),          INTENT(IN)    :: DT
      TYPE(SENSORS_),     INTENT(INOUT) :: SENSORS
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
      INTEGER :: I, K, IU, IV, IT, IEL, IQ1NP
      INTEGER :: NFT_G, MAX_NNODE, TOTAL_NODE_REF, P, Q, U_LEN, V_LEN
      INTEGER :: NX_FOUND, NY_FOUND, NX_CAND, NY_CAND
      INTEGER :: NKNOT_U, NKNOT_V, NGP_Q1NP, IPT_Q1NP, IERR
      INTEGER :: NUM_GROUP_NODE, GPOS, LOCAL_ID
      INTEGER :: IEXPAN, ISTRAIN, ILAY, SZ_IX, TH_STRAIN
      INTEGER :: Q1NP_IDS(MVSIZ), II(6)
      INTEGER, ALLOCATABLE :: NCTRL_ELEM(:)
      INTEGER, ALLOCATABLE :: ELEM_U(:), ELEM_V(:), PID_ELEM(:)
      INTEGER, ALLOCATABLE :: MAT_ID_ELEM(:), NGL_ELEM(:)
      INTEGER, ALLOCATABLE :: NODE_GID(:,:), NODE_LID(:,:), NODE_POS(:,:)
      INTEGER, ALLOCATABLE :: GROUP_GID(:), GROUP_LID(:)
      TYPE(G_BUFEL_) ,POINTER :: GBUF
      TYPE(L_BUFEL_) ,POINTER :: LBUF
      TYPE(ELBUF_STRUCT_) :: ELBUF_TAB_LOCAL(1)
      my_real :: XI, ETA, ZETA, GPW
      my_real, ALLOCATABLE :: X_ELEM(:,:,:), V_ELEM(:,:,:)
      my_real, ALLOCATABLE :: F_INT_ELEM(:,:,:)
      my_real, ALLOCATABLE :: MASS_ELEM(:,:), STIG_ELEM(:,:)
      my_real, ALLOCATABLE :: U_KNOT(:,:), V_KNOT(:,:)
      my_real, ALLOCATABLE :: NVAL(:), DN_LOCAL(:,:), DN_GLOBAL(:,:)
      my_real, ALLOCATABLE :: MATB_GP(:,:)
      my_real, ALLOCATABLE :: VX_BAL(:,:), VY_BAL(:,:), VZ_BAL(:,:)
      my_real, ALLOCATABLE :: XX_BAL(:,:), YY_BAL(:,:), ZZ_BAL(:,:)
      my_real, ALLOCATABLE :: VGAUSS(:,:)
      my_real :: STI(MVSIZ), OFF(MVSIZ), RHO0(MVSIZ)
      my_real :: VOLN(MVSIZ), VOLG(MVSIZ), DVOL(MVSIZ), VD2(MVSIZ)
      my_real :: DELTAX(MVSIZ), DIVDE(MVSIZ)
      my_real :: VIS(MVSIZ), QVIS(MVSIZ), CXX(MVSIZ)
      my_real :: S1(MVSIZ), S2(MVSIZ), S3(MVSIZ), S4(MVSIZ), S5(MVSIZ), S6(MVSIZ)
      my_real :: DXX(MVSIZ), DYY(MVSIZ), DZZ(MVSIZ)
      my_real :: DXY(MVSIZ), DYX(MVSIZ), DYZ(MVSIZ), DZY(MVSIZ), DZX(MVSIZ), DXZ(MVSIZ)
      my_real :: D4(MVSIZ), D5(MVSIZ), D6(MVSIZ)
      my_real :: WXX(MVSIZ), WYY(MVSIZ), WZZ(MVSIZ)
      my_real :: AJ1(MVSIZ), AJ2(MVSIZ), AJ3(MVSIZ), AJ4(MVSIZ), AJ5(MVSIZ), AJ6(MVSIZ)
      my_real :: VDX(MVSIZ), VDY(MVSIZ), VDZ(MVSIZ), MUVOID(MVSIZ)
      my_real :: SSP_EQ(MVSIZ), AIRE(MVSIZ), SIGY(MVSIZ), ET(MVSIZ)
      my_real :: BUFVOIS(MVSIZ), R3_DAM(MVSIZ), AMU(MVSIZ)
      my_real :: MFXX(MVSIZ), MFXY(MVSIZ), MFXZ(MVSIZ), MFYX(MVSIZ), MFYY(MVSIZ)
      my_real :: MFYZ(MVSIZ), MFZX(MVSIZ), MFZY(MVSIZ), MFZZ(MVSIZ)
      my_real :: GAMA(MVSIZ,6), FR_WAV(MVSIZ), TEMPEL(MVSIZ), DIE(MVSIZ)
      my_real :: VARNL(MVSIZ), CONDE(MVSIZ)
      my_real :: FVD2(MVSIZ), FDELTAX(MVSIZ), FSSP(MVSIZ), FQVIS(MVSIZ)
      my_real :: FHEAT(MVSIZ)
      DOUBLE PRECISION :: VOLDP(MVSIZ)
      my_real :: DUMMY_FLUX(1,1)
      my_real :: DTFAC1(102), DTMIN1(102), PERCENT_ADDMASS
      my_real :: DT_STOP_PERCENT_ADDMASS, MASS0_START, PERCENT_ADDMASS_OLD
  !-----------------------------------------------
  !   VECT01 common block re-declaration
  !-----------------------------------------------
      INTEGER :: V01_LFT, V01_LLT, V01_NFT, V01_MTN, V01_IAD,        &
     &           V01_ITY, V01_NPT, V01_JALE, V01_ISMSTR, V01_JEUL,   &
     &           V01_JTUR, V01_JTHE, V01_JLAG, V01_JMULT, V01_JHBE,  &
     &           V01_JIVF, V01_NVAUX, V01_JPOR, V01_JCVT, V01_JSPH,  &
     &           V01_JCLOSE, V01_JPLASOL, V01_IREP, V01_IINT,        &
     &           V01_IHET, V01_IGTYP, V01_ISORTH, V01_ISORTHG,       &
     &           V01_ISRAT, V01_ISROT, V01_ICSEN, V01_IFAILURE,      &
     &           V01_JSMS, V01_ISPH2SOL, V01_IPARTSPH, V01_IGRE,     &
     &           V01_IFORMDT
      COMMON /VECT11/ V01_LFT, V01_LLT, V01_NFT, V01_MTN, V01_IAD,  &
     &                V01_ITY, V01_NPT, V01_JALE, V01_ISMSTR,        &
     &                V01_JEUL, V01_JTUR, V01_JTHE, V01_JLAG,        &
     &                V01_JMULT, V01_JHBE, V01_JIVF, V01_NVAUX,      &
     &                V01_JPOR, V01_JCVT, V01_JSPH, V01_JCLOSE,      &
     &                V01_JPLASOL, V01_IREP, V01_IINT, V01_IHET,     &
     &                V01_IGTYP, V01_ISORTH, V01_ISORTHG, V01_ISRAT, &
     &                V01_ISROT, V01_ICSEN, V01_IFAILURE, V01_JSMS,  &
     &                V01_ISPH2SOL, V01_IPARTSPH, V01_IGRE,          &
     &                V01_IFORMDT

  !=======================================================================
  !   S o u r c e  L i n e s
  !=======================================================================
  !-----------------------------------------------------------------------
  !  (1) Identify Q1NP element IDs for the group and derive degrees (P,Q)
  !-----------------------------------------------------------------------
      NFT_G = IPARG(3,NG)
      MAX_NNODE = 0
      TOTAL_NODE_REF = 0

      DO K = 1, NEL
        DO IQ1NP = 1, NUMELQ1NP_G
          IF (KQ1NP_TAB(5,IQ1NP) == IXS(NIXS,NFT_G + K)) EXIT
        END DO

        Q1NP_IDS(K) = IQ1NP
        MAX_NNODE = MAX(MAX_NNODE, KQ1NP_TAB(3,IQ1NP) + 4)
        TOTAL_NODE_REF = TOTAL_NODE_REF + KQ1NP_TAB(3,IQ1NP) + 4
        P = KQ1NP_TAB(8,IQ1NP)
        Q = KQ1NP_TAB(9,IQ1NP)
      END DO

  !-----------------------------------------------------------------------
  !  (2) Reconstruct Q1NP parametric grid + initialize Gauss scheme
  !-----------------------------------------------------------------------
      CALL Q1NP_REBUILD_GRID(Q1NP_NX_G, Q1NP_NY_G, P, Q)

      U_LEN = Q1NP_NX_G + 2*P + 1
      V_LEN = Q1NP_NY_G + 2*Q + 1

      CALL Q1NP_INIT_GAUSS_SCHEME_STARTER(P + 1, Q + 1, 2)

      ! Total Number of Gauss points in the Q1NP element
      NGP_Q1NP = Q1NP_NP_U_G * Q1NP_NP_V_G * Q1NP_NP_T_G

  !-----------------------------------------------------------------------
  !  (3) Allocate and initialize vectorized work arrays
  !-----------------------------------------------------------------------
      GBUF => ELBUF_STR%GBUF
      ELBUF_TAB_LOCAL(1) = ELBUF_STR

      ALLOCATE(NCTRL_ELEM(NEL))
      ALLOCATE(ELEM_U(NEL), ELEM_V(NEL), PID_ELEM(NEL))
      ALLOCATE(MAT_ID_ELEM(NEL), NGL_ELEM(NEL))
      ALLOCATE(NODE_GID(MAX_NNODE,NEL), NODE_LID(MAX_NNODE,NEL), NODE_POS(MAX_NNODE,NEL))
      ALLOCATE(GROUP_GID(TOTAL_NODE_REF), GROUP_LID(TOTAL_NODE_REF))
      ALLOCATE(X_ELEM(3,MAX_NNODE,NEL), V_ELEM(3,MAX_NNODE,NEL))
      ALLOCATE(F_INT_ELEM(3,MAX_NNODE,NEL))
      ALLOCATE(MASS_ELEM(MAX_NNODE,NEL), STIG_ELEM(MAX_NNODE,NEL))
      ALLOCATE(U_KNOT(Q1NP_NX_G + 2*P + 1, NEL))
      ALLOCATE(V_KNOT(Q1NP_NY_G + 2*Q + 1, NEL))
      ALLOCATE(NVAL(MAX_NNODE), DN_LOCAL(MAX_NNODE,3), DN_GLOBAL(MAX_NNODE,3))
      ALLOCATE(MATB_GP(3*MAX_NNODE,NEL))
      ALLOCATE(VGAUSS(NGP_Q1NP,NEL))

  !-----------------------------------------------------------------------
  !  (4) Snapshot the Gauss-point volume into VOL0DP
  !-----------------------------------------------------------------------
      IF (TT == ZERO) THEN
        DO IT = 1, Q1NP_NP_T_G
          DO IV = 1, Q1NP_NP_V_G
            DO IU = 1, Q1NP_NP_U_G
              LBUF => ELBUF_STR%BUFLY(1)%LBUF(IU,IV,IT)
              LBUF%VOL0DP(1:NEL) = LBUF%VOL(1:NEL)
            END DO
          END DO
        END DO
      END IF

  !-----------------------------------------------------------------------
  !  (5) Build node/group bookkeeping and gather X/V into element-local arrays
  !-----------------------------------------------------------------------
      CALL Q1NP_INIT_NODE_MAP()

  !-----------------------------------------------------------------------
  !  (6) Initialize element fields (OFF, RHO0, DELTAX, buffers, etc.)
  !-----------------------------------------------------------------------
      CALL Q1NP_INIT_ELEM_FIELDS()

  !-----------------------------------------------------------------------
  !  (7) Small-rotation coordinate transform for small strain storage
  !-----------------------------------------------------------------------
      CALL SMALLA3(GBUF%SMSTR, GBUF%OFF, OFF, WXX, WYY, WZZ, NEL, V01_ISMSTR, V01_JLAG)

  !-----------------------------------------------------------------------
  !  (8) Full Gauss integration (IU,IV,IT) over all lanes (IEL=1..NEL)
  !     Each GP contributes: geometry + material law + internal force accumulation
  !-----------------------------------------------------------------------

      IPT_Q1NP = 0
      DO IT = 1, Q1NP_NP_T_G
        ZETA = Q1NP_GP_T_G(IT)
        DO IV = 1, Q1NP_NP_V_G
          ETA = Q1NP_GP_V_G(IV)
          DO IU = 1, Q1NP_NP_U_G
            XI = Q1NP_GP_U_G(IU)

            ! Gauss point weight
            GPW = Q1NP_GW_U_G(IU) * Q1NP_GW_V_G(IV) * Q1NP_GW_T_G(IT)
            IPT_Q1NP = IPT_Q1NP + 1

            ! Reset the Gauss point fields
            CALL Q1NP_RESET_GP()

            ! Compute the geometry terms (derivatives + gauss point volume)
            CALL Q1NP_GP_GEOM(XI, ETA, ZETA, GPW, IU, IV, IT, IERR)
            IF (IERR /= 0) RETURN

            ! Compute the material and strain evaluation
            CALL Q1NP_GP_MAT(IU, IV, IT)

            ! Integrate the Gauss-point stress contributions into the internal force buffers
            CALL Q1NP_ACCUM_NFORCE(IU, IV, IT)
          END DO
        END DO
      END DO

  !-----------------------------------------------------------------------
  !  (9) Update OFF/GBUF%OFF scaling after Gauss integration
  !-----------------------------------------------------------------------
      CALL SMALLB3(GBUF%OFF, OFF, NEL, V01_ISMSTR)

  !-----------------------------------------------------------------------
  !  (10) Assemble element-local forces into global A + update dt + stress averages
  !-----------------------------------------------------------------------
      CALL Q1NP_ASSEMBLE_FINT() ! Assemble internal forces into the global array A
      CALL Q1NP_AVG_SIG_BILAN() ! Update the element-local stress tensor and compute the stress average and bilans

      CONTAINS

!=======================================================================
! Node/group bookkeeping and gather X/V into element-local arrays
!=======================================================================
        SUBROUTINE Q1NP_INIT_NODE_MAP()
          NODE_GID = 0 ! Global node ID
          NODE_LID = 0 ! Local node ID
          NODE_POS = 0 ! Position of the node in the group
          GROUP_GID = 0 ! Global group ID
          GROUP_LID = 0 ! Local group ID
          X_ELEM = ZERO ! Element-local coordinate
          V_ELEM = ZERO ! Element-local velocity
          F_INT_ELEM = ZERO ! Element-local internal force
          MASS_ELEM = ZERO ! Element-local mass
          STIG_ELEM = ZERO ! Element-local stiffness
          U_KNOT = ZERO ! Element-local knot vector
          V_KNOT = ZERO ! Element-local knot vector
          MATB_GP = ZERO ! Material matrix for the Gauss point
          VGAUSS = ZERO ! Gauss point values
          NUM_GROUP_NODE = 0 ! Number of group nodes

          ! Extract common U,V knot vectors once per group; they are identical for all Q1NP elements in this group
          CALL Q1NP_GET_KNOT_VECTORS(Q1NP_NX_G, Q1NP_NY_G, P, Q, &
     &                               Q1NP_KTAB, U_KNOT(1:U_LEN,1), &
     &                               V_KNOT(1:V_LEN,1))

          ! Build the element-local node lists (control points + 4 bulk nodes)
          DO IEL = 1, NEL
            IQ1NP = Q1NP_IDS(IEL)
            NCTRL_ELEM(IEL) = KQ1NP_TAB(3,IQ1NP)

            ELEM_U(IEL) = KQ1NP_TAB(6,IQ1NP)
            ELEM_V(IEL) = KQ1NP_TAB(7,IQ1NP)

            PID_ELEM(IEL) = KQ1NP_TAB(2,IQ1NP)
            MAT_ID_ELEM(IEL) = KQ1NP_TAB(1,IQ1NP)

            NGL_ELEM(IEL) = KQ1NP_TAB(5,IQ1NP)

            ! Reuse the group-level U,V knot vectors for each element lane
            U_KNOT(1:U_LEN,IEL) = U_KNOT(1:U_LEN,1)
            V_KNOT(1:V_LEN,IEL) = V_KNOT(1:V_LEN,1)

            ! Get the control point IDs
            DO K = 1, NCTRL_ELEM(IEL)
              NODE_GID(K,IEL) = IQ1NP_TAB(KQ1NP_TAB(4,IQ1NP) + K - 1)
            END DO
            ! Get the bulk node IDs
            DO K = 1, 4
              NODE_GID(NCTRL_ELEM(IEL) + K, IEL) = IQ1NP_BULK_TAB(KQ1NP_TAB(14,IQ1NP) + K - 1)
            END DO

            ! Get the group node IDs
            DO K = 1, NCTRL_ELEM(IEL) + 4
              GPOS = FIND_GROUP_NODE(NODE_GID(K,IEL), GROUP_GID, NUM_GROUP_NODE)
              IF (GPOS <= 0) THEN
                NUM_GROUP_NODE = NUM_GROUP_NODE + 1
                GROUP_GID(NUM_GROUP_NODE) = NODE_GID(K,IEL)
                GPOS = NUM_GROUP_NODE
              END IF
              NODE_POS(K,IEL) = GPOS
            END DO
          END DO

          ! Build local-id lookup table
          DO I = 1, NUM_GROUP_NODE
            GROUP_LID(I) = FIND_GROUP_NODE(GROUP_GID(I), ITAB_DEBUG, SIZE(ITAB_DEBUG))
          END DO

          ! Element coordinate and velocity arrays
          DO IEL = 1, NEL
            DO K = 1, NCTRL_ELEM(IEL) + 4
              LOCAL_ID = GROUP_LID(NODE_POS(K,IEL))
              NODE_LID(K,IEL) = LOCAL_ID
              X_ELEM(1,K,IEL) = X(1,LOCAL_ID)
              X_ELEM(2,K,IEL) = X(2,LOCAL_ID)
              X_ELEM(3,K,IEL) = X(3,LOCAL_ID)
              V_ELEM(1,K,IEL) = V(1,LOCAL_ID)
              V_ELEM(2,K,IEL) = V(2,LOCAL_ID)
              V_ELEM(3,K,IEL) = V(3,LOCAL_ID)
            END DO
          END DO
        END SUBROUTINE Q1NP_INIT_NODE_MAP

!=======================================================================
! Initialize element fields (OFF, RHO0, DELTAX, buffers, etc.)
! Compute characteristic length per element that is later reused in material law and time-step logic.
!=======================================================================
        SUBROUTINE Q1NP_INIT_ELEM_FIELDS()
          INTEGER :: I_INIT, IEL_LOCAL

          ! Initialize the arrays for the Q1NP element
          DO I_INIT = 1, 6
            II(I_INIT) = NEL * (I_INIT - 1)
          END DO

          ISTRAIN = 1
          ILAY = 1
          IEXPAN = IPARG(49,NG)
          TH_STRAIN = 0
          SZ_IX = NUMELQ + NUMELS + NSVOIS

          STI = ZERO
          OFF = ONE
          RHO0 = ZERO
          VOLN = ZERO
          VOLG = ZERO
          DVOL = ZERO
          VD2 = ZERO
          DELTAX = ZERO
          DIVDE = ZERO
          VIS = ZERO
          QVIS = ZERO
          CXX = ZERO
          S1 = ZERO
          S2 = ZERO
          S3 = ZERO
          S4 = ZERO
          S5 = ZERO
          S6 = ZERO
          DXX = ZERO
          DYY = ZERO
          DZZ = ZERO
          DXY = ZERO
          DYX = ZERO
          DYZ = ZERO
          DZY = ZERO
          DZX = ZERO
          DXZ = ZERO
          D4 = ZERO
          D5 = ZERO
          D6 = ZERO
          WXX = ZERO
          WYY = ZERO
          WZZ = ZERO
          AJ1 = ZERO
          AJ2 = ZERO
          AJ3 = ZERO
          AJ4 = ZERO
          AJ5 = ZERO
          AJ6 = ZERO
          VDX = ZERO
          VDY = ZERO
          VDZ = ZERO
          MUVOID = ZERO
          SSP_EQ = ZERO
          AIRE = ZERO
          SIGY = ZERO
          ET = ZERO
          BUFVOIS = ZERO
          R3_DAM = ZERO
          AMU = ZERO
          MFXX = ZERO
          MFXY = ZERO
          MFXZ = ZERO
          MFYX = ZERO
          MFYY = ZERO
          MFYZ = ZERO
          MFZX = ZERO
          MFZY = ZERO
          MFZZ = ZERO
          GAMA = ZERO
          FR_WAV = ZERO
          TEMPEL = ZERO
          DIE = ZERO
          VARNL = ZERO
          CONDE = ZERO
          FVD2 = ZERO
          FDELTAX = ZERO
          FSSP = ZERO
          FQVIS = ZERO
          DUMMY_FLUX = ZERO

          IF (ASSOCIATED(GBUF%SIG)) GBUF%SIG(1:6*NEL) = ZERO
          IF (ASSOCIATED(GBUF%EINT)) GBUF%EINT(1:NEL) = ZERO
          IF (ASSOCIATED(GBUF%RHO)) GBUF%RHO(1:NEL) = ZERO
          IF (ASSOCIATED(GBUF%QVIS)) GBUF%QVIS(1:NEL) = ZERO
          IF (GBUF%G_PLA > 0 .AND. ASSOCIATED(GBUF%PLA)) GBUF%PLA(1:NEL) = ZERO
          IF (IPARG(40,NG) > 0 .AND. ASSOCIATED(GBUF%EPSD)) GBUF%EPSD(1:NEL) = ZERO

          DO IEL_LOCAL = 1, NEL
            OFF(IEL_LOCAL) = GBUF%OFF(IEL_LOCAL)
            RHO0(IEL_LOCAL) = PM(1,MAT_ID_ELEM(IEL_LOCAL))
            CALL Q1NP_CHAR_LEN(IEL_LOCAL, DELTAX(IEL_LOCAL))
          END DO
        END SUBROUTINE Q1NP_INIT_ELEM_FIELDS

!=======================================================================
! Reset Gauss point fields to zero
!=======================================================================
        SUBROUTINE Q1NP_RESET_GP()
          VIS(1:NEL) = ZERO
          QVIS(1:NEL) = ZERO
          CXX(1:NEL) = ZERO
          MATB_GP(1:3*MAX_NNODE,1:NEL) = ZERO
          S1(1:NEL) = ZERO
          S2(1:NEL) = ZERO
          S3(1:NEL) = ZERO
          S4(1:NEL) = ZERO
          S5(1:NEL) = ZERO
          S6(1:NEL) = ZERO
          DXX(1:NEL) = ZERO
          DYY(1:NEL) = ZERO
          DZZ(1:NEL) = ZERO
          DXY(1:NEL) = ZERO
          DYX(1:NEL) = ZERO
          DYZ(1:NEL) = ZERO
          DZY(1:NEL) = ZERO
          DZX(1:NEL) = ZERO
          DXZ(1:NEL) = ZERO
          D4(1:NEL) = ZERO
          D5(1:NEL) = ZERO
          D6(1:NEL) = ZERO
          WXX(1:NEL) = ZERO
          WYY(1:NEL) = ZERO
          WZZ(1:NEL) = ZERO
          AJ1(1:NEL) = ZERO
          AJ2(1:NEL) = ZERO
          AJ3(1:NEL) = ZERO
          AJ4(1:NEL) = ZERO
          AJ5(1:NEL) = ZERO
          AJ6(1:NEL) = ZERO
          VDX(1:NEL) = ZERO
          VDY(1:NEL) = ZERO
          VDZ(1:NEL) = ZERO
          MUVOID(1:NEL) = ZERO
          SSP_EQ(1:NEL) = ZERO
          AIRE(1:NEL) = ZERO
          SIGY(1:NEL) = ZERO
          ET(1:NEL) = ZERO
          BUFVOIS(1:NEL) = ZERO
          R3_DAM(1:NEL) = ZERO
          AMU(1:NEL) = ZERO
          MFXX(1:NEL) = ZERO
          MFXY(1:NEL) = ZERO
          MFXZ(1:NEL) = ZERO
          MFYX(1:NEL) = ZERO
          MFYY(1:NEL) = ZERO
          MFYZ(1:NEL) = ZERO
          MFZX(1:NEL) = ZERO
          MFZY(1:NEL) = ZERO
          MFZZ(1:NEL) = ZERO
          GAMA(1:NEL,1:6) = ZERO
          FR_WAV(1:NEL) = ZERO
          TEMPEL(1:NEL) = ZERO
          DIE(1:NEL) = ZERO
          VARNL(1:NEL) = ZERO
          CONDE(1:NEL) = ZERO
          FVD2(1:NEL) = ZERO
          FDELTAX(1:NEL) = ZERO
          FSSP(1:NEL) = ZERO
          FQVIS(1:NEL) = ZERO
          STI(1:NEL) = ZERO
        END SUBROUTINE Q1NP_RESET_GP

!=======================================================================
! Geometry terms (derivatives + Gauss point volume)
!=======================================================================
        SUBROUTINE Q1NP_GP_GEOM(XI, ETA, ZETA, GPW, IU, IV, IT, IERR_OUT)
          my_real, INTENT(IN)  :: XI, ETA, ZETA, GPW
          INTEGER, INTENT(IN)  :: IU, IV, IT
          INTEGER, INTENT(OUT) :: IERR_OUT
          INTEGER :: IEL_LOCAL, NNODE_LOCAL, K_LOCAL
          my_real :: JMAT_LOCAL(3,3), JINV_LOCAL(3,3), DETJ_LOCAL
          INTEGER :: IERR_LOCAL

          IERR_OUT = 0

          DO IEL_LOCAL = 1, NEL
            NNODE_LOCAL = NCTRL_ELEM(IEL_LOCAL) + 4

            CALL Q1NP_SHAPE_FUNCTIONS(XI, ETA, ZETA, P, Q, &
     &                                  U_KNOT(1:U_LEN,IEL_LOCAL), V_KNOT(1:V_LEN,IEL_LOCAL), &
     &                                  ELEM_U(IEL_LOCAL), ELEM_V(IEL_LOCAL), &
     &                                  NVAL(1:NNODE_LOCAL), DN_LOCAL(1:NNODE_LOCAL,1:3))

            CALL Q1NP_JACOBIAN(DN_LOCAL(1:NNODE_LOCAL,1:3), X_ELEM(1:3,1:NNODE_LOCAL,IEL_LOCAL), &
     &                           NNODE_LOCAL, JMAT_LOCAL, DETJ_LOCAL, JINV_LOCAL, DN_GLOBAL(1:NNODE_LOCAL,1:3), IERR_LOCAL)

            IF (IERR_LOCAL /= 0) THEN
              write (*,*) 'TT = ', TT
              WRITE(*,'(A,I10,A,I10,A,I4,A,I4,A,I4)') &
     &            'Q1NP ERROR: singular/non-positive Jacobian in element ', NGL_ELEM(IEL_LOCAL), &
     &            ' group ', NG, ' at GP ', IU, ',', IV, ',', IT
              IERR_OUT = IERR_LOCAL
              RETURN
            END IF

            VOLN(IEL_LOCAL) = GPW * DETJ_LOCAL ! Gauss point volume
            VOLDP(IEL_LOCAL) = DBLE(VOLN(IEL_LOCAL)) ! Double precision Gauss point volume
            VGAUSS(IPT_Q1NP,IEL_LOCAL) = VOLN(IEL_LOCAL) ! Gauss point volume
            VOLG(IEL_LOCAL) = VOLG(IEL_LOCAL) + VOLN(IEL_LOCAL) ! Total Gauss point volume

            IF (IDTMIN(101) == 1 .AND. Q1NP_IS_ACTIVE(IEL_LOCAL)) THEN
              DO K_LOCAL = 1, NNODE_LOCAL
                MASS_ELEM(K_LOCAL,IEL_LOCAL) = MASS_ELEM(K_LOCAL,IEL_LOCAL) + &
     &                                        PM(89,MAT_ID_ELEM(IEL_LOCAL)) * NVAL(K_LOCAL) * VOLN(IEL_LOCAL)
              END DO
            END IF

            CALL Q1NP_FILL_MATB(IEL_LOCAL, NNODE_LOCAL, DN_GLOBAL(1:NNODE_LOCAL,1:3))
          END DO

          CALL Q1NP_EVAL_DEF()
        END SUBROUTINE Q1NP_GP_GEOM

!=======================================================================
! Material law + strain/stress evaluation at the Gauss point
!=======================================================================
        SUBROUTINE Q1NP_GP_MAT(IU, IV, IT)
          INTEGER, INTENT(IN) :: IU, IV, IT

          LBUF => ELBUF_STR%BUFLY(1)%LBUF(IU,IV,IT)

          ! rotates/updates stress with spin terms (WXX/WYY/WZZ)
          CALL SROTA3(LBUF%SIG, S1, S2, S3, S4, S5, S6, WXX, WYY, WZZ, &
     &                  NEL, V01_MTN, IPARG(9,NG))

          ! volumetric strain increment
          DIVDE(1:NEL) = DT1 * (DXX(1:NEL) + DYY(1:NEL) + DZZ(1:NEL))

          ! call the density update
          CALL SRHO3(PM, LBUF%VOL, LBUF%RHO, LBUF%EINT, DIVDE, DUMMY_FLUX, FV, &
     &                 VOLN, DVOL, NGL_ELEM, MAT_ID_ELEM, OFF, IPARG(64,NG), &
     &                 GBUF%TAG22, VOLDP, LBUF%VOL0DP, AMU, GBUF%OFF, NEL, &
     &                 V01_MTN, V01_JALE, V01_ISMSTR, V01_JEUL, V01_JLAG, 1, 1, 0)

          ! main material-law driver
          CALL MMAIN(TIMERS, OUTPUT, ELBUF_TAB_LOCAL, 1, PM, GEO, ALE_CONNECT, IXS, IPARG, &
     &                 V, TF, NPF, BUFMAT, STI, X, DT2T, NELTST, ITYPTST, OFFSET, NEL, W, &
     &                 OFF, PID_ELEM, MAT_ID_ELEM, NGL_ELEM, VOLN, VD2, DVOL, DELTAX, VIS, &
     &                 QVIS, CXX, S1, S2, S3, S4, S5, S6, DXX, DYY, DZZ, D4, D5, D6, WXX, &
     &                 WYY, WZZ, AJ1, AJ2, AJ3, AJ4, AJ5, AJ6, VDX, VDY, VDZ, MUVOID, SSP_EQ, &
     &                 AIRE, SIGY, ET, BUFVOIS, LBUF%PLA, R3_DAM, AMU, MFXX, MFXY, MFXZ, MFYX, &
     &                 MFYY, MFYZ, MFZX, MFZY, MFZZ, IPM, GAMA, FR_WAV, DXY, DYX, DYZ, DZY, &
     &                 DZX, DXZ, ISTRAIN, TEMPEL, DIE, IEXPAN, ILAY, MSSA, DMELS, IU, IV, IT, &
     &                 TABLE, FVD2, FDELTAX, FSSP, FQVIS, IPARG(:,NG), IGEO, CONDE, ITASK, &
     &                 NLOC_DMG, VARNL, MAT_ELEM, H3D_STRAIN, V01_JPLASOL, V01_JSPH, MVSIZ, &
     &                 SNPC, STF, SBUFMAT, GLOB_THERM, SVIS, SZ_IX, IRESP, 0, TH_STRAIN, &
     &                 NGROUP, TT, DT1, NTABLE, NUMELQ, NUMMAT, NUMGEO, NUMNOD, NUMELS, &
     &                 IDEL7NOK, IDTMIN, MAXFUNC, IMON_MAT, USERL_AVAIL, IMPL_S, IDYNA, DT, &
     &                 FHEAT, SENSORS)

          ! strain history update routine. Used to store the small-strain tensor
          ! Update LBUF%STRA (stored strain components) from current deformation-rate terms
          CALL SSTRA3(DXX, DYY, DZZ, D4, D5, D6, LBUF%STRA, WXX, WYY, WZZ, OFF, NEL, V01_JCVT)
        END SUBROUTINE Q1NP_GP_MAT

!=======================================================================
! Fill the element-local material basis matrix
!=======================================================================
        SUBROUTINE Q1NP_FILL_MATB(IEL_LOCAL, NNODE_LOCAL, DRDX_LOCAL) 
          INTEGER, INTENT(IN) :: IEL_LOCAL, NNODE_LOCAL
          my_real, INTENT(IN) :: DRDX_LOCAL(NNODE_LOCAL,3)
          INTEGER :: K_LOCAL, IAD_LOCAL

          DO K_LOCAL = 1, NNODE_LOCAL
            IAD_LOCAL = (K_LOCAL - 1) * 3
            MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) = DRDX_LOCAL(K_LOCAL,1)
            MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) = DRDX_LOCAL(K_LOCAL,2)
            MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) = DRDX_LOCAL(K_LOCAL,3)
          END DO
        END SUBROUTINE Q1NP_FILL_MATB

!=======================================================================  
! Evaluate deformation-rate/kinematic terms used for stress update
! computes deformation-rates (DXX, DYY, DZZ, DXY, DYX, DYZ, DZY, DZX, DXZ)
! applies time-step corrections (DT1D2_LOCAL)
!=======================================================================
        SUBROUTINE Q1NP_EVAL_DEF()
          INTEGER :: IEL_LOCAL, K_LOCAL, IAD_LOCAL
          my_real :: DT1D2_LOCAL, AAA_LOCAL

          DO IEL_LOCAL = 1, NEL
            IF (.NOT. Q1NP_IS_ACTIVE(IEL_LOCAL)) CYCLE

            DXX(IEL_LOCAL) = ZERO
            DYY(IEL_LOCAL) = ZERO
            DZZ(IEL_LOCAL) = ZERO
            DXY(IEL_LOCAL) = ZERO
            DYX(IEL_LOCAL) = ZERO
            DYZ(IEL_LOCAL) = ZERO
            DZY(IEL_LOCAL) = ZERO
            DZX(IEL_LOCAL) = ZERO
            DXZ(IEL_LOCAL) = ZERO
            WXX(IEL_LOCAL) = ZERO
            WYY(IEL_LOCAL) = ZERO
            WZZ(IEL_LOCAL) = ZERO
            D4(IEL_LOCAL) = ZERO
            D5(IEL_LOCAL) = ZERO
            D6(IEL_LOCAL) = ZERO

            DO K_LOCAL = 1, NCTRL_ELEM(IEL_LOCAL) + 4
              IAD_LOCAL = (K_LOCAL - 1) * 3
              DXX(IEL_LOCAL) = DXX(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) * V_ELEM(1,K_LOCAL,IEL_LOCAL)
              DYY(IEL_LOCAL) = DYY(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) * V_ELEM(2,K_LOCAL,IEL_LOCAL)
              DZZ(IEL_LOCAL) = DZZ(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) * V_ELEM(3,K_LOCAL,IEL_LOCAL)
              DXY(IEL_LOCAL) = DXY(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) * V_ELEM(1,K_LOCAL,IEL_LOCAL)
              DYX(IEL_LOCAL) = DYX(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) * V_ELEM(2,K_LOCAL,IEL_LOCAL)
              DYZ(IEL_LOCAL) = DYZ(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) * V_ELEM(2,K_LOCAL,IEL_LOCAL)
              DZY(IEL_LOCAL) = DZY(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) * V_ELEM(3,K_LOCAL,IEL_LOCAL)
              DZX(IEL_LOCAL) = DZX(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) * V_ELEM(3,K_LOCAL,IEL_LOCAL)
              DXZ(IEL_LOCAL) = DXZ(IEL_LOCAL) + MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) * V_ELEM(1,K_LOCAL,IEL_LOCAL)
            END DO
          END DO

          DT1D2_LOCAL = HALF * DT1
          DO IEL_LOCAL = 1, NEL
            IF (.NOT. Q1NP_IS_ACTIVE(IEL_LOCAL)) CYCLE

            DXX(IEL_LOCAL) = DXX(IEL_LOCAL) - DT1D2_LOCAL * &
     &                       (DXX(IEL_LOCAL)*DXX(IEL_LOCAL) + DYX(IEL_LOCAL)*DYX(IEL_LOCAL) + DZX(IEL_LOCAL)*DZX(IEL_LOCAL))
            DYY(IEL_LOCAL) = DYY(IEL_LOCAL) - DT1D2_LOCAL * &
     &                       (DYY(IEL_LOCAL)*DYY(IEL_LOCAL) + DZY(IEL_LOCAL)*DZY(IEL_LOCAL) + DXY(IEL_LOCAL)*DXY(IEL_LOCAL))
            DZZ(IEL_LOCAL) = DZZ(IEL_LOCAL) - DT1D2_LOCAL * &
     &                       (DZZ(IEL_LOCAL)*DZZ(IEL_LOCAL) + DXZ(IEL_LOCAL)*DXZ(IEL_LOCAL) + DYZ(IEL_LOCAL)*DYZ(IEL_LOCAL))

            AAA_LOCAL = DT1D2_LOCAL * &
     &                  (DXX(IEL_LOCAL)*DXY(IEL_LOCAL) + DYX(IEL_LOCAL)*DYY(IEL_LOCAL) + DZX(IEL_LOCAL)*DZY(IEL_LOCAL))
            DXY(IEL_LOCAL) = DXY(IEL_LOCAL) - AAA_LOCAL
            DYX(IEL_LOCAL) = DYX(IEL_LOCAL) - AAA_LOCAL
            D4(IEL_LOCAL) = DXY(IEL_LOCAL) + DYX(IEL_LOCAL)

            AAA_LOCAL = DT1D2_LOCAL * &
     &                  (DYY(IEL_LOCAL)*DYZ(IEL_LOCAL) + DZY(IEL_LOCAL)*DZZ(IEL_LOCAL) + DXY(IEL_LOCAL)*DXZ(IEL_LOCAL))
            DYZ(IEL_LOCAL) = DYZ(IEL_LOCAL) - AAA_LOCAL
            DZY(IEL_LOCAL) = DZY(IEL_LOCAL) - AAA_LOCAL
            D5(IEL_LOCAL) = DYZ(IEL_LOCAL) + DZY(IEL_LOCAL)

            AAA_LOCAL = DT1D2_LOCAL * &
     &                  (DZZ(IEL_LOCAL)*DZX(IEL_LOCAL) + DXZ(IEL_LOCAL)*DXX(IEL_LOCAL) + DYZ(IEL_LOCAL)*DYX(IEL_LOCAL))
            DXZ(IEL_LOCAL) = DXZ(IEL_LOCAL) - AAA_LOCAL
            DZX(IEL_LOCAL) = DZX(IEL_LOCAL) - AAA_LOCAL
            D6(IEL_LOCAL) = DXZ(IEL_LOCAL) + DZX(IEL_LOCAL)

            WXX(IEL_LOCAL) = DT1D2_LOCAL * (DZY(IEL_LOCAL) - DYZ(IEL_LOCAL))
            WYY(IEL_LOCAL) = DT1D2_LOCAL * (DXZ(IEL_LOCAL) - DZX(IEL_LOCAL))
            WZZ(IEL_LOCAL) = DT1D2_LOCAL * (DYX(IEL_LOCAL) - DXY(IEL_LOCAL))
          END DO
        END SUBROUTINE Q1NP_EVAL_DEF

!=======================================================================
! Integrate nodal internal forces from stresses into element-local buffers
!=======================================================================
        SUBROUTINE Q1NP_ACCUM_FINT(IEL_LOCAL, SIG1_IN, SIG2_IN, SIG3_IN, SIG4_IN, SIG5_IN, SIG6_IN, &
     &                                 FX_SUM_OUT, FY_SUM_OUT, FZ_SUM_OUT)
          INTEGER, INTENT(IN) :: IEL_LOCAL
          my_real, INTENT(IN) :: SIG1_IN, SIG2_IN, SIG3_IN, SIG4_IN, SIG5_IN, SIG6_IN
          my_real, INTENT(OUT) :: FX_SUM_OUT, FY_SUM_OUT, FZ_SUM_OUT
          INTEGER :: K_LOCAL, IAD_LOCAL
          my_real :: SUMX_LOCAL, SUMY_LOCAL, SUMZ_LOCAL
          my_real :: FX_LOCAL, FY_LOCAL, FZ_LOCAL
          my_real :: STIN_LOCAL, AA_LOCAL
          my_real :: FCOMP_LOCAL(3,6)

          FX_SUM_OUT = ZERO
          FY_SUM_OUT = ZERO
          FZ_SUM_OUT = ZERO
          SUMX_LOCAL = ZERO
          SUMY_LOCAL = ZERO
          SUMZ_LOCAL = ZERO

          DO K_LOCAL = 1, NCTRL_ELEM(IEL_LOCAL) + 4
            IAD_LOCAL = (K_LOCAL - 1) * 3
            SUMX_LOCAL = SUMX_LOCAL + ABS(MATB_GP(IAD_LOCAL + 1,IEL_LOCAL))
            SUMY_LOCAL = SUMY_LOCAL + ABS(MATB_GP(IAD_LOCAL + 2,IEL_LOCAL))
            SUMZ_LOCAL = SUMZ_LOCAL + ABS(MATB_GP(IAD_LOCAL + 3,IEL_LOCAL))
          END DO

          AA_LOCAL = PM(89,MAT_ID_ELEM(IEL_LOCAL)) * SSP_EQ(IEL_LOCAL) * SSP_EQ(IEL_LOCAL)

          DO K_LOCAL = 1, NCTRL_ELEM(IEL_LOCAL) + 4
            IAD_LOCAL = (K_LOCAL - 1) * 3
            FCOMP_LOCAL = ZERO
            FCOMP_LOCAL(1,1) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) * SIG1_IN
            FCOMP_LOCAL(2,2) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) * SIG2_IN
            FCOMP_LOCAL(3,3) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) * SIG3_IN
            FCOMP_LOCAL(1,4) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) * SIG4_IN
            FCOMP_LOCAL(2,4) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) * SIG4_IN
            FCOMP_LOCAL(2,5) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) * SIG5_IN
            FCOMP_LOCAL(3,5) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 2,IEL_LOCAL) * SIG5_IN
            FCOMP_LOCAL(1,6) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 3,IEL_LOCAL) * SIG6_IN
            FCOMP_LOCAL(3,6) = -VOLN(IEL_LOCAL) * MATB_GP(IAD_LOCAL + 1,IEL_LOCAL) * SIG6_IN

            FX_LOCAL = SUM(FCOMP_LOCAL(1,1:6))
            FY_LOCAL = SUM(FCOMP_LOCAL(2,1:6))
            FZ_LOCAL = SUM(FCOMP_LOCAL(3,1:6))

            F_INT_ELEM(1,K_LOCAL,IEL_LOCAL) = F_INT_ELEM(1,K_LOCAL,IEL_LOCAL) + FX_LOCAL
            F_INT_ELEM(2,K_LOCAL,IEL_LOCAL) = F_INT_ELEM(2,K_LOCAL,IEL_LOCAL) + FY_LOCAL
            F_INT_ELEM(3,K_LOCAL,IEL_LOCAL) = F_INT_ELEM(3,K_LOCAL,IEL_LOCAL) + FZ_LOCAL
            FX_SUM_OUT = FX_SUM_OUT + FX_LOCAL
            FY_SUM_OUT = FY_SUM_OUT + FY_LOCAL
            FZ_SUM_OUT = FZ_SUM_OUT + FZ_LOCAL

            STIN_LOCAL = HALF * VOLN(IEL_LOCAL) * &
     &                   (ABS(MATB_GP(IAD_LOCAL + 1,IEL_LOCAL)) * SUMX_LOCAL + &
     &                    ABS(MATB_GP(IAD_LOCAL + 2,IEL_LOCAL)) * SUMY_LOCAL + &
     &                    ABS(MATB_GP(IAD_LOCAL + 3,IEL_LOCAL)) * SUMZ_LOCAL)
            STIG_ELEM(K_LOCAL,IEL_LOCAL) = STIG_ELEM(K_LOCAL,IEL_LOCAL) + STIN_LOCAL * AA_LOCAL

          END DO
        END SUBROUTINE Q1NP_ACCUM_FINT

!=======================================================================
! Uses stress (SIG1..SIG6) to compute and accumulate internal nodal forces via Q1NP_ACCUM_FIN
!=======================================================================
        SUBROUTINE Q1NP_ACCUM_NFORCE(IU, IV, IT)
          INTEGER, INTENT(IN) :: IU, IV, IT
          INTEGER :: IEL_LOCAL
          my_real :: SIG1_LOCAL, SIG2_LOCAL, SIG3_LOCAL
          my_real :: SIG4_LOCAL, SIG5_LOCAL, SIG6_LOCAL
          my_real :: FX_SUM, FY_SUM, FZ_SUM

          LBUF => ELBUF_STR%BUFLY(1)%LBUF(IU,IV,IT)

          DO IEL_LOCAL = 1, NEL
            IF (.NOT. Q1NP_IS_ACTIVE(IEL_LOCAL)) CYCLE

            CALL Q1NP_BUILD_SIG(IEL_LOCAL, SIG1_LOCAL, SIG2_LOCAL, SIG3_LOCAL, &
     &                                      SIG4_LOCAL, SIG5_LOCAL, SIG6_LOCAL)

            CALL Q1NP_ACCUM_FINT(IEL_LOCAL, SIG1_LOCAL, SIG2_LOCAL, SIG3_LOCAL, SIG4_LOCAL, SIG5_LOCAL, &
     &                               SIG6_LOCAL, FX_SUM, FY_SUM, FZ_SUM)
          END DO
        END SUBROUTINE Q1NP_ACCUM_NFORCE

!=======================================================================
! Assemble element-local internal forces into the global array A
!=======================================================================
        SUBROUTINE Q1NP_ASSEMBLE_FINT()
          INTEGER :: IEL_LOCAL, K_LOCAL, LOCAL_ID_LOCAL

          DO IEL_LOCAL = 1, NEL
            IF (.NOT. Q1NP_IS_ACTIVE(IEL_LOCAL)) CYCLE

            DO K_LOCAL = 1, NCTRL_ELEM(IEL_LOCAL) + 4
              LOCAL_ID_LOCAL = NODE_LID(K_LOCAL,IEL_LOCAL)
              A(1,LOCAL_ID_LOCAL) = A(1,LOCAL_ID_LOCAL) + F_INT_ELEM(1,K_LOCAL,IEL_LOCAL)
              A(2,LOCAL_ID_LOCAL) = A(2,LOCAL_ID_LOCAL) + F_INT_ELEM(2,K_LOCAL,IEL_LOCAL)
              A(3,LOCAL_ID_LOCAL) = A(3,LOCAL_ID_LOCAL) + F_INT_ELEM(3,K_LOCAL,IEL_LOCAL)
              STIFN(LOCAL_ID_LOCAL) = STIFN(LOCAL_ID_LOCAL) + STIG_ELEM(K_LOCAL,IEL_LOCAL)
            END DO
          END DO
        END SUBROUTINE Q1NP_ASSEMBLE_FINT
       
!=======================================================================
! Gauss stress average + optional bilans + small-strain housekeeping
! TODO: #REVIEW: Q1NP_AVG_SIG_BILAN
!=======================================================================
        SUBROUTINE Q1NP_AVG_SIG_BILAN()
          INTEGER :: IU_LOCAL, IV_LOCAL, IT_LOCAL
          INTEGER :: IEL_LOCAL, K_LOCAL
          INTEGER :: IPT_Q1NP_LOCAL
          LOGICAL :: CAN_BILAN_LOCAL

          IPT_Q1NP_LOCAL = 0
          DO IT_LOCAL = 1, Q1NP_NP_T_G
            DO IV_LOCAL = 1, Q1NP_NP_V_G
              DO IU_LOCAL = 1, Q1NP_NP_U_G
                IPT_Q1NP_LOCAL = IPT_Q1NP_LOCAL + 1
                LBUF => ELBUF_STR%BUFLY(1)%LBUF(IU_LOCAL,IV_LOCAL,IT_LOCAL)
                CALL IG3DAVERAGE(LBUF%SIG, GBUF%SIG, LBUF%VOL, GBUF%VOL, LBUF%RHO, LBUF%EINT, &
     &                             GBUF%EINT, GBUF%RHO, VGAUSS(IPT_Q1NP_LOCAL,:), VOLG, LBUF%PLA, GBUF%PLA, &
     &                             GBUF%G_PLA, LBUF%EPSD, GBUF%EPSD, NEL, IPARG(40,NG))
              END DO
            END DO
          END DO

          CAN_BILAN_LOCAL = .TRUE.
          DO IEL_LOCAL = 2, NEL
            IF (NCTRL_ELEM(IEL_LOCAL) /= NCTRL_ELEM(1)) THEN
              CAN_BILAN_LOCAL = .FALSE.
              EXIT
            END IF
          END DO

          IF (CAN_BILAN_LOCAL .AND. IPRI > 0) THEN
            ALLOCATE(VX_BAL(NCTRL_ELEM(1) + 4,MVSIZ), VY_BAL(NCTRL_ELEM(1) + 4,MVSIZ), VZ_BAL(NCTRL_ELEM(1) + 4,MVSIZ))
            ALLOCATE(XX_BAL(NCTRL_ELEM(1) + 4,MVSIZ), YY_BAL(NCTRL_ELEM(1) + 4,MVSIZ), ZZ_BAL(NCTRL_ELEM(1) + 4,MVSIZ))

            DO IEL_LOCAL = 1, NEL
              DO K_LOCAL = 1, NCTRL_ELEM(IEL_LOCAL) + 4
                VX_BAL(K_LOCAL,IEL_LOCAL) = V_ELEM(1,K_LOCAL,IEL_LOCAL)
                VY_BAL(K_LOCAL,IEL_LOCAL) = V_ELEM(2,K_LOCAL,IEL_LOCAL)
                VZ_BAL(K_LOCAL,IEL_LOCAL) = V_ELEM(3,K_LOCAL,IEL_LOCAL)
                XX_BAL(K_LOCAL,IEL_LOCAL) = X_ELEM(1,K_LOCAL,IEL_LOCAL)
                YY_BAL(K_LOCAL,IEL_LOCAL) = X_ELEM(2,K_LOCAL,IEL_LOCAL)
                ZZ_BAL(K_LOCAL,IEL_LOCAL) = X_ELEM(3,K_LOCAL,IEL_LOCAL)
              END DO
            END DO

            CALL IGE3DBILAN(PARTSAV, GBUF%EINT, GBUF%RHO, VOLG, VX_BAL, VY_BAL, VZ_BAL, IPARTS, &
     &                      GBUF%VOL, GRESAV, GRTH, IGRTH, XX_BAL, YY_BAL, ZZ_BAL, NCTRL_ELEM(1) + 4, &
     &                      ITASK, IPARG(1,NG), SENSORS)

            DEALLOCATE(VX_BAL, VY_BAL, VZ_BAL, XX_BAL, YY_BAL, ZZ_BAL)
          END IF
        END SUBROUTINE Q1NP_AVG_SIG_BILAN

!=======================================================================
! Build the element-local stress tensor
!=======================================================================
        SUBROUTINE Q1NP_BUILD_SIG(IEL_LOCAL, SIG1_OUT, SIG2_OUT, SIG3_OUT, SIG4_OUT, SIG5_OUT, SIG6_OUT)
          INTEGER, INTENT(IN) :: IEL_LOCAL
          my_real, INTENT(OUT) :: SIG1_OUT, SIG2_OUT, SIG3_OUT, SIG4_OUT, SIG5_OUT, SIG6_OUT
          SIG1_OUT = LBUF%SIG(II(1)+IEL_LOCAL) + SVIS(IEL_LOCAL,1) - QVIS(IEL_LOCAL)
          SIG2_OUT = LBUF%SIG(II(2)+IEL_LOCAL) + SVIS(IEL_LOCAL,2) - QVIS(IEL_LOCAL)
          SIG3_OUT = LBUF%SIG(II(3)+IEL_LOCAL) + SVIS(IEL_LOCAL,3) - QVIS(IEL_LOCAL)
          SIG4_OUT = LBUF%SIG(II(4)+IEL_LOCAL) + SVIS(IEL_LOCAL,4)
          SIG5_OUT = LBUF%SIG(II(5)+IEL_LOCAL) + SVIS(IEL_LOCAL,5)
          SIG6_OUT = LBUF%SIG(II(6)+IEL_LOCAL) + SVIS(IEL_LOCAL,6)
        END SUBROUTINE Q1NP_BUILD_SIG

!=======================================================================
! Check if the element is active
!=======================================================================
        LOGICAL FUNCTION Q1NP_IS_ACTIVE(IEL_LOCAL)
          INTEGER, INTENT(IN) :: IEL_LOCAL
          Q1NP_IS_ACTIVE = .TRUE.
          IF (GBUF%OFF(IEL_LOCAL) <= ZERO) Q1NP_IS_ACTIVE = .FALSE.
          IF (OFF(IEL_LOCAL) <= ZERO) Q1NP_IS_ACTIVE = .FALSE.
        END FUNCTION Q1NP_IS_ACTIVE

!=======================================================================
! Calculate the characteristic length of the element
! Used in MMAIN_Q1NP to update the element time-step
!=======================================================================
        SUBROUTINE Q1NP_CHAR_LEN(IEL_LOCAL, DELTAX_OUT)
          INTEGER, INTENT(IN) :: IEL_LOCAL
          my_real, INTENT(OUT) :: DELTAX_OUT
          INTEGER :: N_TOP_LOCAL
          INTEGER :: TOP_1_LOCAL, TOP_2_LOCAL, TOP_3_LOCAL, TOP_4_LOCAL
          INTEGER :: BOT_1_LOCAL, BOT_2_LOCAL, BOT_3_LOCAL, BOT_4_LOCAL
          my_real, PARAMETER :: SPAN_SCALE = 0.25
          my_real :: VOL_LENGTH_LOCAL
          my_real :: LU_TOP_1, LU_TOP_2, LV_TOP_1, LV_TOP_2
          my_real :: LU_BOT_1, LU_BOT_2, LV_BOT_1, LV_BOT_2
          my_real :: LT_1, LT_2, LT_3, LT_4
          my_real :: SPAN_U, SPAN_V, SPAN_T
          my_real :: INV_P_LOCAL, INV_Q_LOCAL

          N_TOP_LOCAL = NCTRL_ELEM(IEL_LOCAL)
          INV_P_LOCAL = ONE / REAL(P, KIND(ONE))
          INV_Q_LOCAL = ONE / REAL(Q, KIND(ONE))
          TOP_1_LOCAL = 1
          TOP_2_LOCAL = P + 1
          TOP_3_LOCAL = N_TOP_LOCAL
          TOP_4_LOCAL = N_TOP_LOCAL - P

          BOT_1_LOCAL = N_TOP_LOCAL + 1
          BOT_2_LOCAL = N_TOP_LOCAL + 2
          BOT_3_LOCAL = N_TOP_LOCAL + 3
          BOT_4_LOCAL = N_TOP_LOCAL + 4

          LU_TOP_1 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_1_LOCAL, TOP_2_LOCAL)
          LU_TOP_2 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_4_LOCAL, TOP_3_LOCAL)
          LV_TOP_1 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_1_LOCAL, TOP_4_LOCAL)
          LV_TOP_2 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_2_LOCAL, TOP_3_LOCAL)

          LU_BOT_1 = Q1NP_NODE_DIST(IEL_LOCAL, BOT_1_LOCAL, BOT_2_LOCAL)
          LU_BOT_2 = Q1NP_NODE_DIST(IEL_LOCAL, BOT_4_LOCAL, BOT_3_LOCAL)
          LV_BOT_1 = Q1NP_NODE_DIST(IEL_LOCAL, BOT_1_LOCAL, BOT_4_LOCAL)
          LV_BOT_2 = Q1NP_NODE_DIST(IEL_LOCAL, BOT_2_LOCAL, BOT_3_LOCAL)

          LT_1 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_1_LOCAL, BOT_1_LOCAL)
          LT_2 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_2_LOCAL, BOT_2_LOCAL)
          LT_3 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_3_LOCAL, BOT_3_LOCAL)
          LT_4 = Q1NP_NODE_DIST(IEL_LOCAL, TOP_4_LOCAL, BOT_4_LOCAL)

          SPAN_U = SPAN_SCALE * (LU_TOP_1 + LU_TOP_2 + LU_BOT_1 + LU_BOT_2) * INV_P_LOCAL
          SPAN_V = SPAN_SCALE * (LV_TOP_1 + LV_TOP_2 + LV_BOT_1 + LV_BOT_2) * INV_Q_LOCAL
          SPAN_T = SPAN_SCALE * (LT_1 + LT_2 + LT_3 + LT_4)

          DELTAX_OUT =  MIN(SPAN_U, MIN(SPAN_V, SPAN_T))

          IF (DELTAX_OUT <= EPSILON(ONE)) THEN
            VOL_LENGTH_LOCAL = MAX(ABS(GBUF%VOL(IEL_LOCAL)), EPSILON(ONE)) ** &
     &                         (ONE / REAL(3, KIND(ONE)))
            DELTAX_OUT = VOL_LENGTH_LOCAL
          END IF

        END SUBROUTINE Q1NP_CHAR_LEN

!=======================================================================
! Calculate the distance between two nodes in the element
!=======================================================================
        my_real FUNCTION Q1NP_NODE_DIST(IEL_LOCAL, NODE_A_LOCAL, NODE_B_LOCAL)
          INTEGER, INTENT(IN) :: IEL_LOCAL, NODE_A_LOCAL, NODE_B_LOCAL
          my_real :: DX_LOCAL, DY_LOCAL, DZ_LOCAL

          DX_LOCAL = X_ELEM(1,NODE_B_LOCAL,IEL_LOCAL) - X_ELEM(1,NODE_A_LOCAL,IEL_LOCAL)
          DY_LOCAL = X_ELEM(2,NODE_B_LOCAL,IEL_LOCAL) - X_ELEM(2,NODE_A_LOCAL,IEL_LOCAL)
          DZ_LOCAL = X_ELEM(3,NODE_B_LOCAL,IEL_LOCAL) - X_ELEM(3,NODE_A_LOCAL,IEL_LOCAL)

          Q1NP_NODE_DIST = SQRT(DX_LOCAL*DX_LOCAL + DY_LOCAL*DY_LOCAL + DZ_LOCAL*DZ_LOCAL)
        END FUNCTION Q1NP_NODE_DIST

        ! ----------------------------------------------------------------------
        ! Find the group node ID for a given node ID
        ! ----------------------------------------------------------------------
        INTEGER FUNCTION FIND_GROUP_NODE(GID, GIDS, NGID)
          INTEGER, INTENT(IN) :: GID, NGID
          INTEGER, INTENT(IN) :: GIDS(:)
          INTEGER :: ITMP
          FIND_GROUP_NODE = 0
          DO ITMP = 1, NGID
            IF (GIDS(ITMP) == GID) THEN
              FIND_GROUP_NODE = ITMP
              RETURN
            END IF
          END DO
        END FUNCTION FIND_GROUP_NODE

!=======================================================================
! Rebuild the grid for a given number of control points
!=======================================================================
        SUBROUTINE Q1NP_REBUILD_GRID(NX_OUT, NY_OUT, P_IN, Q_IN)
          INTEGER, INTENT(OUT) :: NX_OUT, NY_OUT
          INTEGER, INTENT(IN)  :: P_IN, Q_IN
          NX_OUT = 0
          NY_OUT = 0
          IF (SQ1NPCTRL_SHARED_G <= 0 .OR. SQ1NPKNOT_L_G <= 0) RETURN
          DO NX_CAND = 1, SQ1NPCTRL_SHARED_G
            IF (MOD(SQ1NPCTRL_SHARED_G, NX_CAND) /= 0) CYCLE
            NY_CAND = SQ1NPCTRL_SHARED_G / NX_CAND
            NX_FOUND = NX_CAND - P_IN
            NY_FOUND = NY_CAND - Q_IN
            IF (NX_FOUND <= 0 .OR. NY_FOUND <= 0) CYCLE
            NKNOT_U = NX_FOUND + 2*P_IN + 1
            NKNOT_V = NY_FOUND + 2*Q_IN + 1
            IF (NKNOT_U + NKNOT_V == SQ1NPKNOT_L_G) THEN
              NX_OUT = NX_FOUND
              NY_OUT = NY_FOUND
              RETURN
            END IF
          END DO
        END SUBROUTINE Q1NP_REBUILD_GRID

      END SUBROUTINE Q1NP_FORC3
