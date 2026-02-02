Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
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
!||====================================================================
!||    i7mainf                     ../engine/source/interfaces/int07/i7mainf.F
!||--- called by ------------------------------------------------------
!||    intfop2                     ../engine/source/interfaces/interf/intfop2.F
!||--- calls      -----------------------------------------------------
!||    ancmsg                      ../engine/source/output/message/message.F
!||    arret                       ../engine/source/system/arret.F
!||    frictionparts_model_isot    ../engine/source/interfaces/int07/frictionparts_model.F
!||    frictionparts_model_ortho   ../engine/source/interfaces/int07/frictionparts_model.F
!||    i17main_pena                ../engine/source/interfaces/int17/i17main_pena.F
!||    i18dst3                     ../engine/source/interfaces/int18/i18dst3.F
!||    i18for3                     ../engine/source/interfaces/int18/i18for3.F
!||    i7ass33                     ../engine/source/interfaces/int07/i7ass3.F
!||    i7cdcor3                    ../engine/source/interfaces/int07/i7cdcor3.F
!||    i7cor3                      ../engine/source/interfaces/int07/i7cor3.F
!||    i7curvsz                    ../engine/source/interfaces/int07/i7curv.F
!||    i7dst3                      ../engine/source/interfaces/int07/i7dst3.F
!||    i7for3                      ../engine/source/interfaces/int07/i7for3.F
!||    i7norm                      ../engine/source/interfaces/int07/i7curv.F
!||    i7norme                     ../engine/source/interfaces/int07/i7rcurv.F
!||    i7normn                     ../engine/source/interfaces/int07/i7rcurv.F
!||    i7normnp                    ../engine/source/interfaces/int07/i7rcurv.F
!||    i7normp                     ../engine/source/interfaces/int07/i7curv.F
!||    i7rcurv                     ../engine/source/interfaces/int07/i7rcurv.F
!||    i7therm                     ../engine/source/interfaces/int07/i7therm.F
!||    my_barrier                  ../engine/source/system/machine.F
!||    spmd_exch_n                 ../engine/source/mpi/generic/spmd_exch_n.F
!||    spmd_i7curvsz               ../engine/source/mpi/interfaces/spmd_i7curvsz.F
!||    startime                    ../engine/source/system/timer_mod.F90
!||    stoptime                    ../engine/source/system/timer_mod.F90
!||    sum_6_float_sens            ../engine/source/system/parit.F
!||--- uses       -----------------------------------------------------
!||    ale_connectivity_mod        ../common_source/modules/ale/ale_connectivity_mod.F
!||    elbufdef_mod                ../common_source/modules/mat_elem/elbufdef_mod.F90
!||    groupdef_mod                ../common_source/modules/groupdef_mod.F
!||    h3d_mod                     ../engine/share/modules/h3d_mod.F
!||    i18dst3_mod                 ../engine/source/interfaces/int18/i18dst3.F
!||    i18for3_mod                 ../engine/source/interfaces/int18/i18for3.F
!||    intbuf_fric_mod             ../common_source/modules/interfaces/intbuf_fric_mod.F90
!||    intbufdef_mod               ../common_source/modules/interfaces/intbufdef_mod.F90
!||    message_mod                 ../engine/share/message_module/message_mod.F
!||    multi_fvm_mod               ../common_source/modules/ale/multi_fvm_mod.F90
!||    output_mod                  ../common_source/modules/output/output_mod.F90
!||    timer_mod                   ../engine/source/system/timer_mod.F90
!||====================================================================
      SUBROUTINE I7MAINF(OUTPUT, TIMERS,
     1                   IPARI          ,X           ,A              ,ALE_CONNECTIVITY ,XCELL      ,
     2                   ICODT          ,FSAV        ,V              ,MS               ,DT2T       ,
     3                   NELTST         ,ITYPTST     ,ITAB           ,STIFN            ,FSKYI      ,
     4                   ISKY           ,FCONT       ,NIN            ,LINDMAX          ,KINET      ,
     5                   JTASK          ,NB_JLT      ,NB_JLT_NEW     ,NB_STOK_N        ,ELBUF_TAB  ,
     6                   NISKYFI        ,NEWFRONT    ,NSTRF          ,SECFCUM          ,IGROUPS    ,
     7                   ICONTACT       ,VISCN       ,NUM_IMP        ,
     9                   NS_IMP         ,NE_IMP      ,IND_IMP        ,FSAVSUB          ,NRTMDIM    ,
     A                   IGRBRIC        ,
     B                   EMINX          ,IXS         ,IXS16          ,IXS20            ,FNCONT     ,
     C                   FTCONT         ,IAD_ELEM    ,FR_ELEM        ,RCONTACT         ,ACONTACT   ,
     D                   PCONTACT       ,TEMP        ,FTHE           ,FTHESKYI         ,
     E                   PM             ,IPARG       ,IAD17          ,MSKYI_SMS        ,ISKYI_SMS  ,
     F                   NODNX_SMS      ,MS0         ,QFRICINT       ,NPC              ,TF         ,
     G                   CONDN          ,CONDNSKYI   ,INTBUF_TAB     ,NODADT_THERM     ,THEACCFACT , 
     H                   FBSAV6         ,ISENSINT    ,DIMFB          ,IXIG3D           ,KXIG3D     ,
     I                   WIGE           ,KNOT        ,IGEO           ,MULTI_FVM        ,H3D_DATA   ,
     J                   INTBUF_FRIC_TAB,KNOTLOCPC   ,KNOTLOCEL      ,ITASK            ,TAGNCONT   ,
     K                   KLOADPINTER    ,LOADPINTER  ,LOADP_HYD_INTER,DGAPLOADINT      ,S_LOADPINTER,
     L                   INTEREFRIC     ,S_XCELL_REMOTE,XCELL_REMOTE )
C======================================================================
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE TIMER_MOD
      USE ELBUFDEF_MOD 
      USE MESSAGE_MOD
      USE INTBUFDEF_MOD
      USE MULTI_FVM_MOD
      USE H3D_MOD
      USE INTBUF_FRIC_MOD
      USE GROUPDEF_MOD
      USE ALE_CONNECTIVITY_MOD
      USE OUTPUT_MOD
      USE I18DST3_MOD , ONLY : I18DST3
      USE I18FOR3_MOD , ONLY : I18FOR3
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
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "com08_c.inc"
#include      "param_c.inc"
#include      "warn_c.inc"
#include      "task_c.inc"
#include      "parit_c.inc"
#include      "timeri_c.inc"
#include      "impl1_c.inc"
#include      "macro.inc"
#include      "comlock.inc"
#include      "tabsiz_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      TYPE(TIMER_), INTENT(INOUT) :: TIMERS
      TYPE(OUTPUT_), INTENT(INOUT) :: OUTPUT
      TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECTIVITY
      TYPE (ELBUF_STRUCT_), DIMENSION(NGROUP) :: ELBUF_TAB
      INTEGER NELTST,ITYPTST,NIN,NEWFRONT,NSTRF(*),NRTMDIM, IAD17, ISENSINT(*),DIMFB,IGROUPS(NUMELS)
      INTEGER IPARI(NPARI,NINTER), ICODT(*),ICONTACT(*),ITAB(*), ISKY(*), KINET(*),IPARG(NPARG,*),TAGNCONT(NLOADP_HYD_INTER,NUMNOD)
      INTEGER NB_JLT,NB_JLT_NEW,NB_STOK_N,JTASK,NISKYFI, LINDMAX,ITASK
      INTEGER NUM_IMP,NS_IMP(*),NE_IMP(*),IND_IMP(*)
      INTEGER IXS(*)  ,IXS16(*) ,IXS20(*)
      INTEGER IAD_ELEM(2,*),FR_ELEM(*), ISKYI_SMS(*), NODNX_SMS(*),NPC(*),KXIG3D(NIXIG3D,*),IXIG3D(*),IGEO(NPROPGI,*)
      INTEGER, INTENT(IN) :: S_LOADPINTER
      INTEGER, INTENT(IN) :: KLOADPINTER(NINTER+1),LOADPINTER(S_LOADPINTER),LOADP_HYD_INTER(NLOADP_HYD)
      INTEGER, INTENT(IN) :: INTEREFRIC
      INTEGER, INTENT(IN) :: NODADT_THERM
      my_real, intent(in) :: theaccfact
      my_real, INTENT(IN) :: DGAPLOADINT(S_LOADPINTER)
      my_real  EMINX(*), XCELL(3,SXCELL)
      my_real DT2T,X(3,*), A(3,*), FSAV(*), V(3,*),MS(*),STIFN(*),FSKYI(LSKYI,4),FCONT(3,*),MS0(*),
     .        SECFCUM(7,NUMNOD,NSECT),VISCN(*), FSAVSUB(*),
     .        FNCONT(3,*), FTCONT(3,*), RCONTACT(*), ACONTACT(*),
     .        PCONTACT(*),TEMP(*),FTHE(*),FTHESKYI(LSKYI),PM(NPROPM,*),
     .        MSKYI_SMS(*), QFRICINT(*),TF(*),CONDN(*),CONDNSKYI(LSKYI),WIGE(*),KNOT(*)
      my_real KNOTLOCPC(*),KNOTLOCEL(*)
      DOUBLE PRECISION FBSAV6(12,6,DIMFB)
      INTEGER, INTENT(in) :: S_XCELL_REMOTE
      my_real, DIMENSION(S_XCELL_REMOTE), INTENT(in) :: XCELL_REMOTE
      TYPE(INTBUF_STRUCT_) INTBUF_TAB
      TYPE(MULTI_FVM_STRUCT), INTENT(INOUT) :: MULTI_FVM
      TYPE(H3D_DATABASE) :: H3D_DATA
      TYPE(INTBUF_FRIC_STRUCT_), TARGET, DIMENSION(NINTERFRIC) ::  INTBUF_FRIC_TAB
      TYPE (GROUP_)  , DIMENSION(NGRBRIC) :: IGRBRIC
C-----------------------------------------------
C     N E W   V a r i a b l e s
C-----------------------------------------------
      my_real, DIMENSION(:,:,:), ALLOCATABLE :: CONT_ELEMENT  !contact element array
      my_real, DIMENSION(:,:,:), ALLOCATABLE :: load_arr     !load array
      INTEGER, DIMENSION(:,:),     ALLOCATABLE :: CAND_SEC_SEG_ID   !secondary segment nodes
      INTEGER, DIMENSION(:,:),     ALLOCATABLE :: CAND_MST_SEG_ID   !primary segment nodes
      INTEGER, DIMENSION(:),       ALLOCATABLE :: CAND_SEC_SEG      !secondary segment nodes
      INTEGER, DIMENSION(:),       ALLOCATABLE :: node_id_load      !node id array (8 nodes per element)
      INTEGER COUNT !count for the number of contact elements found
      INTEGER ASS_FLAG !Assembly flag, 0 if old Assembly, 1 if new assembly algorithm called
      INTEGER :: MAX_STS_SIZE_ACTUAL !actual calculated sts array size needed      
      INTEGER OPTION !option for the contact algorithm
      INTEGER IMPACT_glob !global impact flag, 0 if no impact, 1 if impact
      LOGICAL :: STS_ACTIVE !true if STS is active, false if NTS is active
      INTEGER :: IP_MAX
C-----------------------------------------------
C   H i s t o r y   V a r i a b l e s   f o r   F r i c t i o n
C-----------------------------------------------
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: XI1_HIST  ! History of xi1 per Gauss point
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: XI2_HIST  ! History of xi2 per Gauss point
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: TTRIAL1_HIST  ! History of T_trial(1) per Gauss point
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: TTRIAL2_HIST  ! History of T_trial(2) per Gauss point
      INTEGER, SAVE :: MAX_STS_SIZE_HIST = 0  ! Stored size for history arrays
      INTEGER, SAVE :: IP_MAX_HIST = 0        ! Stored quadrature order
C-----------------------------------------------
C   G l o b a l   G a u s s   P o i n t   S t a t e   (N e w   S T S)
C-----------------------------------------------
C  1D arrays indexed by a global Gauss point index:
C     gp_index = secondary_element_id * (IP_MAX*IP_MAX) + local_gp_index
C  where local_gp_index = (z-1)*IP_MAX + (q-1)
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_GLOBAL_ID
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_SECONDARY_EL_ID
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_MASTER_EL_ID
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_MASTER_EL_ID_PREV

      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI1_GLOBAL
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI2_GLOBAL
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI1_GLOBAL_PREV
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI2_GLOBAL_PREV

      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI1_PERIOD
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: GP_XI2_PERIOD

      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_TTRIAL1_HIST
      REAL*8 , DIMENSION(:), ALLOCATABLE, SAVE :: GP_TTRIAL2_HIST

      LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: GP_IS_STICKING
      LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: GP_INITIALIZED
      INTEGER :: MAX_GLOBAL_GP
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I, J, I_STOK, JLT_NEW, JLT , NFT, IVIS2,
     .        IBC,  ISECIN, IBAG, IADM,
     .        IGAP, INACTI, IFQ, MFROT, IGSTI, NISUB, ITIED,
     .        JLT_TIED, NB_LOC, I_STOK_LOC,DEBUT,
     .        ILAGM, LENR, LENT, MAXCC,INTTH,IFORM, H, IERROR,
     .        SYM_FLAG_TYPE19, INTFRIC ,NSETPRTS ,NPARTFRIC ,IORTHFRIC,
     .        NFORTH ,NFISOT
      INTEGER IX1(MVSIZ), IX2(MVSIZ), IX3(MVSIZ), IX4(MVSIZ),
     .        NSVG(MVSIZ), CN_LOC(MVSIZ),CE_LOC(MVSIZ),
     .        CAND_N_N(MVSIZ),CAND_E_N(MVSIZ),KINI(MVSIZ),
     .        INDEX2(LINDMAX),
     .        ISDSIZ(NSPMD+1),IRCSIZ(NSPMD+1),ITAG(NUMNOD),
     .        IELECI(MVSIZ), NSMS(MVSIZ), IPARTFRICSI(MVSIZ),
     .        IPARTFRICMI(MVSIZ),INDEXISOT(MVSIZ),INDEXORTH(MVSIZ),
     .        IREP_FRICMI(MVSIZ)
      my_real STARTT, FRIC, GAP, STOPT,VISC,VISCF,STIGLO,GAPMIN, KMIN, KMAX, GAPMAX,RSTIF,FHEATS,FHEATM,TINT,DTMINI
C-----------------------------------------------
      my_real NX1(MVSIZ), NX2(MVSIZ), NX3(MVSIZ), NX4(MVSIZ),
     .        NY1(MVSIZ), NY2(MVSIZ), NY3(MVSIZ), NY4(MVSIZ),
     .        NZ1(MVSIZ), NZ2(MVSIZ), NZ3(MVSIZ), NZ4(MVSIZ),
     .        LB1(MVSIZ), LB2(MVSIZ), LB3(MVSIZ), LB4(MVSIZ),
     .        LC1(MVSIZ), LC2(MVSIZ), LC3(MVSIZ), LC4(MVSIZ),
     .        P1(MVSIZ), P2(MVSIZ), P3(MVSIZ), P4(MVSIZ),
     .        X1(MVSIZ), X2(MVSIZ), X3(MVSIZ), X4(MVSIZ),
     .        Y1(MVSIZ), Y2(MVSIZ), Y3(MVSIZ), Y4(MVSIZ),
     .        Z1(MVSIZ), Z2(MVSIZ), Z3(MVSIZ), Z4(MVSIZ),
     .        XI(MVSIZ), YI(MVSIZ), ZI(MVSIZ), STIF(MVSIZ),
     .        H1(MVSIZ), H2(MVSIZ), H3(MVSIZ), H4(MVSIZ),
     .        GAPV(MVSIZ),VXI(MVSIZ),VYI(MVSIZ),VZI(MVSIZ),MSI(MVSIZ),
     .        TEMPI(MVSIZ),PHI(MVSIZ),AREASI(MVSIZ),
     .        EFRICT(MVSIZ),FNI(MVSIZ),CONDINT(MVSIZ),MSI_APPEND(2*MVSIZ)
      my_real RCURVI(MVSIZ), ANGLMI(MVSIZ), ANGLT, PADM
      INTEGER NRTMFT, NRTMLT, NMNFT, NMNLT, NRADM,IRSTH,IFRIC
      my_real NNX1(MVSIZ), NNX2(MVSIZ), NNX3(MVSIZ), NNX4(MVSIZ),
     .        NNY1(MVSIZ), NNY2(MVSIZ), NNY3(MVSIZ), NNY4(MVSIZ),
     .        NNZ1(MVSIZ), NNZ2(MVSIZ), NNZ3(MVSIZ), NNZ4(MVSIZ),
     .        CMAJ(MVSIZ)
      INTEGER ICURV,SFSAVPARIT,NINLOADP
      my_real DRAD, DRAD2, FRAD,XTHE,XFRIC,DGAPLOADP
      my_real, DIMENSION(:,:,:), ALLOCATABLE :: FSAVPARIT
      my_real KS(MVSIZ),K1(MVSIZ),K2(MVSIZ),K3(MVSIZ),K4(MVSIZ),
     .        CS(MVSIZ),C1(MVSIZ),C2(MVSIZ),C3(MVSIZ),C4(MVSIZ),
     .        KT(MVSIZ),C(MVSIZ),CF(MVSIZ),
     .        FXI(MVSIZ), FYI(MVSIZ), FZI(MVSIZ),
     .        FX1(MVSIZ), FX2(MVSIZ), FX3(MVSIZ), FX4(MVSIZ),
     .        FY1(MVSIZ), FY2(MVSIZ), FY3(MVSIZ), FY4(MVSIZ),
     .        FZ1(MVSIZ), FZ2(MVSIZ), FZ3(MVSIZ), FZ4(MVSIZ),
     .        PHI1(MVSIZ),PHI2(MVSIZ),PHI3(MVSIZ),PHI4(MVSIZ)
      my_real XFILTR_FRIC,FRIC_COEFS(MVSIZ,10),VISCFFRIC(MVSIZ),FRICC(MVSIZ),
     .        FRIC_COEFS2(MVSIZ,10),VISCFFRIC2(MVSIZ),FRICC2(MVSIZ),
     .        DIR1(MVSIZ,3),DIR2(MVSIZ,3),DIR_FRICMI(MVSIZ,2),
     .        XMU(MVSIZ), ! Effective friction coefficients for STS
     .        ALPHA0, QFRICT
      INTEGER, DIMENSION(:) ,POINTER  :: TABCOUPLEPARTS_FRIC
      INTEGER, DIMENSION(:) ,POINTER  :: TABPARTS_FRIC 
      INTEGER, DIMENSION(:) ,POINTER  :: ADPARTS_FRIC 
      INTEGER, DIMENSION(:) ,POINTER  :: IFRICORTH 
      my_real, DIMENSION(:) ,POINTER  :: TABCOEF_FRIC
      INTEGER,TARGET, DIMENSION(1):: TABCOUPLEPARTS_FRIC_BID
      INTEGER,TARGET, DIMENSION(1):: TABPARTS_FRIC_BID
      INTEGER,TARGET, DIMENSION(1):: ADPARTS_FRIC_BID
      INTEGER,TARGET, DIMENSION(1):: IFRICORTH_BID
      my_real,TARGET, DIMENSION(1):: TABCOEF_FRIC_BID
      INTEGER :: NSN, NTY, NOINT
      INTEGER :: NRTM
      LOGICAL TYPE18
      INTEGER :: IDIR
C-----------------------------------------------
C   S o u r c e   L i n e s
C-----------------------------------------------
      TYPE18 = .FALSE.
      NRTM = IPARI(4,NIN)
      NSN = IPARI(5,NIN)
      NTY = IPARI(7,NIN)
      IBC = IPARI(11,NIN)
      IVIS2 = IPARI(14,NIN)
      IF(IPARI(33,NIN) == 1) RETURN  !precondition
      NOINT = IPARI(15,NIN)
      IGAP = IPARI(21,NIN)
      INACTI = IPARI(22,NIN) !inact=7 <=> /INTER/TYPE18
      ISECIN = IPARI(28,NIN)
      MFROT = IPARI(30,NIN)
      IFQ = IPARI(31,NIN)
      IBAG = IPARI(32,NIN)
      IGSTI = IPARI(34,NIN)
      NISUB = IPARI(36,NIN)
      ICURV = IPARI(39,NIN)
      ITIED = IPARI(85,NIN)
      SYM_FLAG_TYPE19 = IPARI(71,NIN)
C heat interface
      INTTH = IPARI(47,NIN)
      IFORM = IPARI(48,NIN)
      IRSTH = IPARI(42,NIN)
C adaptive meshing
      IADM = IPARI(44,NIN)
      NRADM = IPARI(49,NIN)
      PADM = INTBUF_TAB%VARIABLES(24)
      ANGLT = INTBUF_TAB%VARIABLES(25)
      STIGLO = -INTBUF_TAB%STFAC(1)
      STARTT = INTBUF_TAB%VARIABLES(3)
      STOPT = INTBUF_TAB%VARIABLES(11)
      IF(STARTT > TT) RETURN !precondition
      IF(TT > STOPT)  RETURN !precondition
      FRIC = INTBUF_TAB%VARIABLES(1)
      GAP = INTBUF_TAB%VARIABLES(2)
      GAPMIN = INTBUF_TAB%VARIABLES(13)
      VISC = INTBUF_TAB%VARIABLES(14)
      VISCF= INTBUF_TAB%VARIABLES(15)
      IFRIC = 0
      IF(INTTH > 0) IFRIC =IPARI(50,NIN)
      XFRIC  = INTBUF_TAB%VARIABLES(34)
      GAPMAX = INTBUF_TAB%VARIABLES(16)
      KMIN   = INTBUF_TAB%VARIABLES(17)
      KMAX   = INTBUF_TAB%VARIABLES(18)
      RSTIF  = INTBUF_TAB%VARIABLES(20)
      FHEATS = INTBUF_TAB%VARIABLES(21)
      TINT   = INTBUF_TAB%VARIABLES(22)
      DTMINI = INTBUF_TAB%VARIABLES(41)
      FHEATM = INTBUF_TAB%VARIABLES(23)
      XTHE   = INTBUF_TAB%VARIABLES(33)
      FRAD   = ZERO
      DRAD   = ZERO
      DRAD2  = ZERO
      IF(NTY == 7)THEN
        FRAD  = INTBUF_TAB%VARIABLES(31)
        DRAD  = INTBUF_TAB%VARIABLES(32)
        DRAD2 = DRAD*DRAD
      ENDIF

C--- Corresponding Friction model  
      INTFRIC=IPARI(72,NIN)
      IORTHFRIC = 0
      NSETPRTS = 0
      NPARTFRIC = 0 
      XFILTR_FRIC = ZERO
      
      IF(INTFRIC /= 0) THEN 
         TABCOUPLEPARTS_FRIC => INTBUF_FRIC_TAB(INTFRIC)%TABCOUPLEPARTS_FRIC
         TABCOEF_FRIC => INTBUF_FRIC_TAB(INTFRIC)%TABCOEF_FRIC
         TABPARTS_FRIC => INTBUF_FRIC_TAB(INTFRIC)%TABPARTS_FRIC
         ADPARTS_FRIC => INTBUF_FRIC_TAB(INTFRIC)%ADPARTS_FRIC
         XFILTR_FRIC = INTBUF_FRIC_TAB(INTFRIC)%XFILTR_FRIC
         NSETPRTS = INTBUF_FRIC_TAB(INTFRIC)%NSETPRTS
         NPARTFRIC = INTBUF_FRIC_TAB(INTFRIC)%S_TABPARTS_FRIC
         IORTHFRIC = INTBUF_FRIC_TAB(INTFRIC)%IORTHFRIC
         IFRICORTH => INTBUF_FRIC_TAB(INTFRIC)%IFRICORTH
      ELSE
         TABCOUPLEPARTS_FRIC => TABCOUPLEPARTS_FRIC_BID
         TABPARTS_FRIC => TABPARTS_FRIC_BID
         TABCOEF_FRIC => TABCOEF_FRIC_BID
         ADPARTS_FRIC => ADPARTS_FRIC_BID
         IFRICORTH => IFRICORTH_BID
         IF (IFQ/=0) XFILTR_FRIC = INTBUF_TAB%XFILTR(1)
      ENDIF

      NINLOADP = IPARI(95,NIN) ! load pressure related to inter
      DGAPLOADP = INTBUF_TAB%VARIABLES(46)

      IF(NTY == 7)THEN
        IF(INACTI == 7) TYPE18 = .TRUE.
        IF(TYPE18 .AND. IPARI(34,NIN) == -2) RETURN
!----------------------------------------------------
! Compute nodal normals using quadratic curvature
!----------------------------------------------------
        IF(ICURV == 3)THEN
         CALL MY_BARRIER()
         IF(JTASK == 1)THEN
           ALLOCATE(INTBUF_TAB%NODNORM_NORMAL (3,NUMNOD))
           IF(IPARIT == 0)THEN
             CALL I7NORM(IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD,X,INTBUF_TAB%NODNORM_NORMAL,IPARI(6,NIN),INTBUF_TAB%MSR)
             IF(NSPMD > 1)THEN
               LENR = IAD_ELEM(1,NSPMD+1)-IAD_ELEM(1,1)
               CALL SPMD_EXCH_N(INTBUF_TAB%NODNORM_NORMAL,IAD_ELEM,FR_ELEM,LENR)
             END IF
           ELSE
! Parallel assembly (parith/on using SPMD): optimize if necessary
             LENR = IAD_ELEM(1,NSPMD+1)-IAD_ELEM(1,1)
             IF(NSPMD > 1)THEN
! Memory space estimation for parith/on (spmd communication)
               CALL SPMD_I7CURVSZ(
     1           IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD,IAD_ELEM,FR_ELEM,
     2           ISDSIZ      ,IRCSIZ           ,ITAG  ,LENR    ,LENT   ,
     3           MAXCC       )
             ELSE
! Memory space estimation for parith/on (without spmd)
               CALL I7CURVSZ(
     1           IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD,ITAG,LENT,MAXCC)
             END IF
             CALL I7NORMP(
     1         IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD ,X    ,INTBUF_TAB%NODNORM_NORMAL,
     2         IPARI(6,NIN),INTBUF_TAB%MSR  ,LENT   ,MAXCC,ISDSIZ    ,
     3         IRCSIZ      ,IAD_ELEM        ,FR_ELEM,ITAG )
           END IF
         END IF
         CALL MY_BARRIER()
        ENDIF
!----------------------------------------------------
! Curvature radius: compute normalized nodal normals
! Warning / Precondition : IADM /= 0 and ICURV /= 0 cannot be used together (Starter error).
!----------------------------------------------------
        IF(IADM /= 0)THEN
         CALL MY_BARRIER()
         IF(JTASK==1)THEN
           ALLOCATE(INTBUF_TAB%MODRCURV(NRTMDIM),INTBUF_TAB%MODANGLM(NRTMDIM))
           ALLOCATE(INTBUF_TAB%NODNORM_NORMAL (3,NUMNOD))
           IF(IPARIT == 0)THEN
             CALL I7NORMN(
     .         IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD,X ,INTBUF_TAB%NODNORM_NORMAL,
     .                   IPARI(6,NIN),INTBUF_TAB%MSR)
             IF(NSPMD > 1)THEN
               LENR = IAD_ELEM(1,NSPMD+1)-IAD_ELEM(1,1)
               CALL SPMD_EXCH_N(INTBUF_TAB%NODNORM_NORMAL,IAD_ELEM,FR_ELEM,LENR)
             END IF
           ELSE
C Traitement d'assemblage parith/on spmd a optimiser si besoin
             LENR = IAD_ELEM(1,NSPMD+1)-IAD_ELEM(1,1)
             IF(NSPMD > 1)THEN
C Sizing memory space for parith/on and spmd
               CALL SPMD_I7CURVSZ(
     1           IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD,IAD_ELEM,FR_ELEM,
     2           ISDSIZ      ,IRCSIZ           ,ITAG  ,LENR    ,LENT   ,
     3           MAXCC       )
             ELSE
C Sizing memory space for parith/on without spmd
               CALL I7CURVSZ(
     1           IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD,ITAG,LENT,MAXCC)
             END IF
             CALL I7NORMNP(
     1         IPARI(4,NIN),INTBUF_TAB%IRECTM,NUMNOD ,X    ,INTBUF_TAB%NODNORM_NORMAL,
     2         IPARI(6,NIN),INTBUF_TAB%MSR,LENT   ,MAXCC,ISDSIZ    ,
     3         IRCSIZ      ,IAD_ELEM     ,FR_ELEM,ITAG )
           END IF
         END IF
         CALL MY_BARRIER()
         NMNFT=1+(JTASK-1)*IPARI(6,NIN)/NTHREAD
         NMNLT=JTASK*IPARI(6,NIN)/NTHREAD
         CALL I7NORME(NMNFT,NMNLT,INTBUF_TAB%NODNORM_NORMAL,INTBUF_TAB%MSR)
         CALL MY_BARRIER()
         NRTMFT=1+(JTASK-1)*IPARI(4,NIN)/NTHREAD
         NRTMLT=JTASK*IPARI(4,NIN)/NTHREAD
         CALL I7RCURV(NRTMFT, NRTMLT ,X ,INTBUF_TAB%NODNORM_NORMAL ,INTBUF_TAB%IRECTM ,
     .                INTBUF_TAB%MODRCURV , NRADM  ,INTBUF_TAB%MODANGLM ,ANGLT )
         CALL MY_BARRIER()
        END IF
C----------------------------------------------------
C
        I_STOK = INTBUF_TAB%I_STOK(1)
        DEBUT = 0
        I_STOK_LOC = 0
        
C static decoupage
        IF(INACTI /= 7 .OR. IBAG == 0)THEN
          NB_LOC = I_STOK / NTHREAD
          IF (JTASK == NTHREAD) THEN
            I_STOK_LOC = I_STOK-NB_LOC*(NTHREAD-1)
          ELSE
            I_STOK_LOC = NB_LOC
          ENDIF
          DEBUT = (JTASK-1)*NB_LOC
        END IF
        I_STOK = 0

        IF (IMPL_S == 1) THEN
          NUM_IMP = 0
          VISC  =ZERO
          VISCF =ZERO
        ENDIF
C
C updating ISTOK
C
        IF(ITIED == 0)THEN
        IF (INACTI == 5 .OR. INACTI == 6 .OR. TYPE18)THEN
          DO I = DEBUT+1, DEBUT+I_STOK_LOC
            IF(INTBUF_TAB%CAND_N(I)<0) THEN
              I_STOK = I_STOK + 1
              INDEX2(I_STOK) = I
C inbuf == cand_n
              INTBUF_TAB%CAND_N(I) = -INTBUF_TAB%CAND_N(I)
            ELSE
C Discount at 0 of Cand_P
              INTBUF_TAB%CAND_P(I) = ZERO
            ENDIF
          ENDDO
        ELSE
          DO I = DEBUT+1, DEBUT+I_STOK_LOC
            IF(INTBUF_TAB%CAND_N(I)<0) THEN
              I_STOK = I_STOK + 1
              INDEX2(I_STOK) = I
C inbuf == cand_n
              INTBUF_TAB%CAND_N(I) = -INTBUF_TAB%CAND_N(I)
            ENDIF
          ENDDO
        ENDIF
        ELSE ! IF(ITIED == 0)THEN
          IF (INACTI == 5 .OR. INACTI == 6)THEN
            DO I = DEBUT+1, DEBUT+I_STOK_LOC
              IF(INTBUF_TAB%CAND_F(8*(I-1)+1)/=ZERO) THEN
                I_STOK = I_STOK + 1
                INDEX2(I_STOK) = I
                INTBUF_TAB%CAND_N(I) = IABS(INTBUF_TAB%CAND_N(I))
              ELSEIF(INTBUF_TAB%CAND_N(I)<0) THEN
                I_STOK = I_STOK + 1
                INDEX2(I_STOK) = I
C inbuf == cand_n
                INTBUF_TAB%CAND_N(I) = -INTBUF_TAB%CAND_N(I)
              ELSE
C Discount at 0 of Cand_P
                INTBUF_TAB%CAND_P(I) = ZERO
                INTBUF_TAB%CAND_F(8*(I-1)+1:8*(I-1)+8) = ZERO
              ENDIF
            ENDDO
          ELSE
            DO I = DEBUT+1, DEBUT+I_STOK_LOC
              IF(INTBUF_TAB%CAND_F(8*(I-1)+1)/=ZERO) THEN
                I_STOK = I_STOK + 1
                INDEX2(I_STOK) = I
                INTBUF_TAB%CAND_N(I) = IABS(INTBUF_TAB%CAND_N(I))
              ELSEIF(INTBUF_TAB%CAND_N(I)<0) THEN
                I_STOK = I_STOK + 1
                INDEX2(I_STOK) = I
C inbuf == cand_n
                INTBUF_TAB%CAND_N(I) = -INTBUF_TAB%CAND_N(I)
              ELSE
                INTBUF_TAB%CAND_F(8*(I-1)+1:8*(I-1)+8) = ZERO
              ENDIF
            ENDDO
          END IF
        ENDIF

c------------------------------------------------
        IF (DEBUG(3) >= 1) THEN
          NB_JLT = NB_JLT + I_STOK_LOC
          NB_STOK_N = NB_STOK_N + I_STOK
        ENDIF

        SFSAVPARIT = 0
        DO I=1,NISUB+1
          IF(ISENSINT(I)/=0) THEN
            SFSAVPARIT = SFSAVPARIT + 1
          ENDIF
        ENDDO
        IF (SFSAVPARIT /= 0) THEN
          ALLOCATE(FSAVPARIT(NISUB+1,11,I_STOK),STAT=IERROR)
          IF(IERROR/=0) THEN
            CALL ANCMSG(MSGID=19,ANMODE=ANINFO,C1='(/INTER/TYPE7)')
           CALL ARRET(2)
          ENDIF
          DO J=1,I_STOK
            DO I=1,11
              DO H=1,NISUB+1
                FSAVPARIT(H,I,J) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ALLOCATE(FSAVPARIT(0,0,0),STAT=IERROR)
          IF(IERROR/=0) THEN
            CALL ANCMSG(MSGID=19,ANMODE=ANINFO,C1='(/INTER/TYPE7)')
           CALL ARRET(2)
          ENDIF
        ENDIF
        
        !--------------------------------------------------------
        ! FORCE COMPUTATION : 
        !    inter 7 & inter 18
        !--------------------------------------------------------
C!$OMP ORDERED
        DO NFT = 0 , I_STOK - 1 , NVSIZ
          JLT = MIN( NVSIZ, I_STOK - NFT )
          ! preparing retained candidates
          CALL I7CDCOR3(JLT,INDEX2(NFT+1),INTBUF_TAB%CAND_E,INTBUF_TAB%CAND_N,CAND_E_N,CAND_N_N)
          ! cand_n & cand_e replaced with cand_n_n & cand_e_n
          CALL I7COR3(
     1  JLT       ,X       ,INTBUF_TAB%IRECTM,INTBUF_TAB%NSV,CAND_E_N,
     2  CAND_N_N  ,INTBUF_TAB%STFM,INTBUF_TAB%STFNS,X1 ,X2           ,
     3  X3        ,X4         ,Y1        ,Y2           ,Y3           ,
     4  Y4        ,Z1         ,Z2        ,Z3           ,Z4           ,
     5  XI        ,YI         ,ZI        ,STIF         ,IX1          ,
     6  IX2       ,IX3        ,IX4       ,NSVG         ,IGAP         ,
     7      GAP       ,INTBUF_TAB%GAP_S,INTBUF_TAB%GAP_M,GAPV ,
     9  MS        ,VXI        ,VYI       ,
     A  VZI       ,MSI        ,NSN       ,V            ,KINET        ,
     B  KINI      ,NTY        ,NIN       ,IGSTI        ,KMIN         ,  
     C  KMAX      ,GAPMAX     ,GAPMIN    ,IADM         ,INTBUF_TAB%MODRCURV        ,
     D  RCURVI    ,INTBUF_TAB%MODANGLM      ,ANGLMI    ,INTTH        ,TEMP         ,
     E  TEMPI     ,PHI        ,INTBUF_TAB%AREAS,INTBUF_TAB%IELEC,AREASI ,
     F  IELECI    ,NODNX_SMS  ,NSMS      ,INTBUF_TAB%GAP_SL,INTBUF_TAB%GAP_ML,
     G  INTFRIC  ,INTBUF_TAB%IPARTFRICS,IPARTFRICSI,INTBUF_TAB%IPARTFRICM,IPARTFRICMI,
     H  IORTHFRIC,INTBUF_TAB%IREP_FRICM,INTBUF_TAB%DIR_FRICM,IREP_FRICMI,DIR_FRICMI  )
          ! -- PIT

          JLT_NEW = 0
          JLT_TIED= 0

          IF(TYPE18)THEN
            CALL I18DST3(
     1        JLT           ,CAND_N_N  ,CAND_E_N          , CN_LOC                    ,CE_LOC                     ,
     2        X1            ,X2        ,X3                , X4                        ,Y1                         ,
     3        Y2            ,Y3        ,Y4                , Z1                        ,Z2                         ,
     4        Z3            ,Z4        ,XI                , YI                        ,ZI                         ,
     5        NX1           ,NX2       ,NX3               , NX4                       ,NY1                        ,
     6        NY2           ,NY3       ,NY4               , NZ1                       ,NZ2                        ,
     7        NZ3           ,NZ4       ,LB1               , LB2                       ,LB3                        ,
     8        LB4           ,LC1       ,LC2               , LC3                       ,LC4                        ,
     9        P1            ,P2        ,P3                , P4                        ,IX1                        ,
     A        IX2           ,IX3       ,IX4               , NSVG                      ,STIF                       ,
     B        JLT_NEW       ,GAPV      ,INTBUF_TAB%CAND_P ,ALE_CONNECTIVITY%NE_CONNECT,
     C        INDEX2(NFT+1) ,VXI       ,VYI               ,ITAB                       ,XCELL                      ,
     D        VZI           ,MSI       ,KINI              ,
     E        IGAP          ,MULTI_FVM ,S_XCELL_REMOTE ,XCELL_REMOTE)

            IF( MULTI_FVM%IS_INT18_LAW151 ) THEN
              MSI_APPEND(1:MVSIZ) = ZERO
              MSI_APPEND(MVSIZ+1:2*MVSIZ) = MSI(1:MVSIZ)
            ELSE
              MSI_APPEND(1:MVSIZ) = MSI(1:MVSIZ)
              MSI_APPEND(MVSIZ+1:2*MVSIZ) = ZERO
            ENDIF

            JLT = JLT_NEW
            IF (IMONM > 0 .AND. JTASK == 1) CALL STARTIME(TIMERS,20)
            IF(JLT_NEW /= 0) THEN
              IPARI(29,NIN) = 1
              IDIR = IPARI(34,NIN)
              IF (DEBUG(3) >= 1)NB_JLT_NEW = NB_JLT_NEW + JLT_NEW
              CALL I18FOR3(OUTPUT,
     1          JLT               ,A                 ,V                 ,IBC                ,ICODT             ,
     2          FSAV              ,GAP               ,FRIC              ,MS                 ,VISC              ,
     3          VISCF             ,NOINT             ,INTBUF_TAB%STFNS  ,ITAB               ,CN_LOC            ,
     4          STIGLO            ,STIFN             ,STIF              ,FSKYI              ,ISKY              ,
     5          NX1               ,NX2               ,NX3               ,NX4                ,NY1               ,
     6          NY2               ,NY3               ,NY4               ,NZ1                ,NZ2               ,
     7          NZ3               ,NZ4               ,LB1               ,LB2                ,LB3               ,
     8          LB4               ,LC1               ,LC2               ,LC3                ,LC4               ,
     9          P1                ,P2                ,P3                ,P4                 ,FCONT             ,
     B          IX1               ,IX2               ,IX3               ,IX4                ,NSVG              ,
     C          IVIS2             ,NELTST            ,ITYPTST           ,DT2T               ,IXS               ,
     D          GAPV              ,INTBUF_TAB%CAND_P ,INDEX2(NFT+1)     ,NISKYFI           ,
     E          KINET             ,NEWFRONT          ,ISECIN            ,NSTRF              ,SECFCUM           ,
     F          X                 ,INTBUF_TAB%IRECTM ,CE_LOC            ,MFROT              ,IFQ               ,
     G          INTBUF_TAB%FRIC_P ,INTBUF_TAB%FTSAVX ,INTBUF_TAB%FTSAVY ,INTBUF_TAB%FTSAVZ  ,
     H          INTBUF_TAB%IFPEN  , ICONTACT         ,IGROUPS           ,IPARG              ,
     J          VISCN             ,VXI               ,VYI               ,VZI                ,MSI_APPEND        ,
     K          KINI              ,NIN               ,NISUB             ,INTBUF_TAB%LISUB   ,INTBUF_TAB%ADDSUBS,
     L          INTBUF_TAB%ADDSUBM,INTBUF_TAB%LISUBS ,INTBUF_TAB%LISUBM ,FSAVSUB            ,INTBUF_TAB%CAND_N ,
     M          IPARI(33,NIN)     ,IPARI(39,NIN)     ,FNCONT            ,MS0                ,
     N          JTASK             ,ISENSINT          ,FSAVPARIT         ,NFT                ,MULTI_FVM         ,
     O          H3D_DATA          ,ELBUF_TAB         ,IDIR             )
            ENDIF
            
       ELSE ! NOT TYPE18


C     The maximum number of contact pairs is limited by the number of candidates
C
      MAX_STS_SIZE_ACTUAL = 1  ! Minimum size
      IF (INTBUF_TAB%S_CAND_N > 0) THEN
        ! Maximum contact pairs cannot exceed number of candidates
        MAX_STS_SIZE_ACTUAL = MAX(MAX_STS_SIZE_ACTUAL, INTBUF_TAB%S_CAND_N)
      ENDIF
      ! Add safety factor (some candidates might map to same segment pair, but we need buffer)
      MAX_STS_SIZE_ACTUAL = MAX_STS_SIZE_ACTUAL * 2  ! Safety factor      

      ! Adjust working array sizes to match MAX_STS_SIZE_ACTUAL
      IF (ALLOCATED(CONT_ELEMENT)) DEALLOCATE(CONT_ELEMENT)
      IF (ALLOCATED(load_arr))     DEALLOCATE(load_arr)
      IF (ALLOCATED(CAND_SEC_SEG_ID)) DEALLOCATE(CAND_SEC_SEG_ID)
      IF (ALLOCATED(CAND_MST_SEG_ID)) DEALLOCATE(CAND_MST_SEG_ID)
      IF (ALLOCATED(CAND_SEC_SEG))    DEALLOCATE(CAND_SEC_SEG)
      IF (ALLOCATED(node_id_load))    DEALLOCATE(node_id_load)

      ALLOCATE(CONT_ELEMENT(MAX_STS_SIZE_ACTUAL,3,8))
      ALLOCATE(load_arr(MAX_STS_SIZE_ACTUAL,8,4))
      ALLOCATE(CAND_SEC_SEG_ID(MAX_STS_SIZE_ACTUAL,5))
      ALLOCATE(CAND_MST_SEG_ID(MAX_STS_SIZE_ACTUAL,5))
      ALLOCATE(CAND_SEC_SEG(MAX_STS_SIZE_ACTUAL))
      ALLOCATE(node_id_load(MAX_STS_SIZE_ACTUAL*8))

C-----------------------------------------------
C   Allocate/Reallocate History Arrays for Friction
C-----------------------------------------------
      ! Quadrature order (currently fixed to 2, but make it flexible)
      IP_MAX = 2  ! Current quadrature order
      
      ! Maximum number of Gauss points for STS pairs (all secondary elements)
      ! MAX_GLOBAL_GP = MAX_STS_SIZE_ACTUAL * IP_MAX * IP_MAX
      ! Used for 1D global Gauss point state arrays
      MAX_GLOBAL_GP = MAX_STS_SIZE_ACTUAL * IP_MAX * IP_MAX
      
      ! Check if reallocation is needed for element-based history arrays
      IF (.NOT. ALLOCATED(XI1_HIST) .OR. 
     &    MAX_STS_SIZE_ACTUAL > MAX_STS_SIZE_HIST .OR.
     &    IP_MAX /= IP_MAX_HIST) THEN
        
        ! Deallocate if already allocated
        IF (ALLOCATED(XI1_HIST)) THEN
          DEALLOCATE(XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST)
        ENDIF
        
        ! Allocate with new size
        ALLOCATE(XI1_HIST(MAX_STS_SIZE_ACTUAL, IP_MAX, IP_MAX))
        ALLOCATE(XI2_HIST(MAX_STS_SIZE_ACTUAL, IP_MAX, IP_MAX))
        ALLOCATE(TTRIAL1_HIST(MAX_STS_SIZE_ACTUAL, IP_MAX, IP_MAX))
        ALLOCATE(TTRIAL2_HIST(MAX_STS_SIZE_ACTUAL, IP_MAX, IP_MAX))

        ! Initialize to zero (will be set on first contact)
        XI1_HIST = 0.d0
        XI2_HIST = 0.d0
        TTRIAL1_HIST = 0.d0
        TTRIAL2_HIST = 0.d0
        
        ! Store current size
        MAX_STS_SIZE_HIST = MAX_STS_SIZE_ACTUAL
        IP_MAX_HIST = IP_MAX
      ENDIF

C-----------------------------------------------
C   Allocate/Reallocate Global Gauss Point State Arrays
C-----------------------------------------------
      IF (.NOT. ALLOCATED(GP_XI1_GLOBAL) .OR.
     &    MAX_STS_SIZE_ACTUAL > MAX_STS_SIZE_HIST .OR.
     &    IP_MAX /= IP_MAX_HIST) THEN

        ! Deallocate if already allocated
        IF (ALLOCATED(GP_GLOBAL_ID))        DEALLOCATE(GP_GLOBAL_ID)
        IF (ALLOCATED(GP_SECONDARY_EL_ID))  DEALLOCATE(GP_SECONDARY_EL_ID)
        IF (ALLOCATED(GP_MASTER_EL_ID))     DEALLOCATE(GP_MASTER_EL_ID)
        IF (ALLOCATED(GP_MASTER_EL_ID_PREV)) DEALLOCATE(GP_MASTER_EL_ID_PREV)

        IF (ALLOCATED(GP_XI1_GLOBAL))       DEALLOCATE(GP_XI1_GLOBAL)
        IF (ALLOCATED(GP_XI2_GLOBAL))       DEALLOCATE(GP_XI2_GLOBAL)
        IF (ALLOCATED(GP_XI1_GLOBAL_PREV))  DEALLOCATE(GP_XI1_GLOBAL_PREV)
        IF (ALLOCATED(GP_XI2_GLOBAL_PREV))  DEALLOCATE(GP_XI2_GLOBAL_PREV)

        IF (ALLOCATED(GP_XI1_PERIOD))       DEALLOCATE(GP_XI1_PERIOD)
        IF (ALLOCATED(GP_XI2_PERIOD))       DEALLOCATE(GP_XI2_PERIOD)

        IF (ALLOCATED(GP_TTRIAL1_HIST))     DEALLOCATE(GP_TTRIAL1_HIST)
        IF (ALLOCATED(GP_TTRIAL2_HIST))     DEALLOCATE(GP_TTRIAL2_HIST)

        IF (ALLOCATED(GP_IS_STICKING))      DEALLOCATE(GP_IS_STICKING)
        IF (ALLOCATED(GP_INITIALIZED))      DEALLOCATE(GP_INITIALIZED)

        ! Allocate global Gauss point state arrays
        ALLOCATE(GP_GLOBAL_ID       (MAX_GLOBAL_GP))
        ALLOCATE(GP_SECONDARY_EL_ID (MAX_GLOBAL_GP))
        ALLOCATE(GP_MASTER_EL_ID    (MAX_GLOBAL_GP))
        ALLOCATE(GP_MASTER_EL_ID_PREV(MAX_GLOBAL_GP))

        ALLOCATE(GP_XI1_GLOBAL      (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI2_GLOBAL      (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI1_GLOBAL_PREV (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI2_GLOBAL_PREV (MAX_GLOBAL_GP))

        ALLOCATE(GP_XI1_PERIOD      (MAX_GLOBAL_GP))
        ALLOCATE(GP_XI2_PERIOD      (MAX_GLOBAL_GP))

        ALLOCATE(GP_TTRIAL1_HIST    (MAX_GLOBAL_GP))
        ALLOCATE(GP_TTRIAL2_HIST    (MAX_GLOBAL_GP))

        ALLOCATE(GP_IS_STICKING     (MAX_GLOBAL_GP))
        ALLOCATE(GP_INITIALIZED     (MAX_GLOBAL_GP))

        ! Initialize state
        GP_GLOBAL_ID        = 0
        GP_SECONDARY_EL_ID  = 0
        GP_MASTER_EL_ID     = 0
        GP_MASTER_EL_ID_PREV= 0

        GP_XI1_GLOBAL       = 0.d0
        GP_XI2_GLOBAL       = 0.d0
        GP_XI1_GLOBAL_PREV  = 0.d0
        GP_XI2_GLOBAL_PREV  = 0.d0

        GP_XI1_PERIOD       = 0
        GP_XI2_PERIOD       = 0

        GP_TTRIAL1_HIST     = 0.d0
        GP_TTRIAL2_HIST     = 0.d0

        GP_IS_STICKING      = .FALSE.
        GP_INITIALIZED      = .FALSE.
      ENDIF

C==========================================================
C CONTACT ALGORITHM SELECTION
C   STS_ACTIVE = .TRUE.  -> NEW STS (GAUSS/LOBATTO)
C   STS_ACTIVE = .FALSE. -> OLD NTS
C==========================================================

        IF (.NOT. STS_ACTIVE) THEN
C==========================================================
C OLD NTS ALGORITHM
C==========================================================
        WRITE(*,*) 
        WRITE(*,*) "***********************************"
        WRITE(*,*) "NTS ALGORITHM STARTED"
        WRITE(*,*) "***********************************"

        !----------------------------------------------------------
        ! Distance calculation and projection for contact detection
        !----------------------------------------------------------  
        CALL I7DST3(
     1        JLT          ,CAND_N_N,CAND_E_N,CN_LOC ,CE_LOC       ,
     2        X1           ,X2      ,X3      ,X4     ,Y1           ,
     3        Y2           ,Y3      ,Y4      ,Z1     ,Z2           ,
     4        Z3           ,Z4      ,XI      ,YI     ,ZI           ,
     5        NX1          ,NX2     ,NX3     ,NX4    ,NY1          ,
     6        NY2          ,NY3     ,NY4     ,NZ1    ,NZ2          ,
     7        NZ3          ,NZ4     ,LB1     ,LB2    ,LB3          ,
     8        LB4          ,LC1     ,LC2     ,LC3    ,LC4          ,
     9        P1           ,P2      ,P3      ,P4     ,IX1          ,
     A        IX2          ,IX3     ,IX4     ,NSVG   ,STIF         ,
     B        JLT_NEW      ,GAPV    ,INACTI ,INTBUF_TAB%CAND_P,
     C        INDEX2(NFT+1),VXI     ,VYI     ,
     D        VZI          ,MSI     ,KINI    ,ICURV  ,INTBUF_TAB%IRECTM,
     E        NNX1         ,NNX2    ,NNX3    ,NNX4   ,NNY1         ,
     F        NNY2         ,NNY3    ,NNY4    ,NNZ1   ,NNZ2         ,
     G        NNZ3         ,NNZ4    ,INTBUF_TAB%NODNORM_NORMAL ,IADM,RCURVI       ,
     H        ANGLMI       ,INTTH   ,TEMPI   ,PHI     ,AREASI      ,
     I        IELECI       ,NSMS    ,CMAJ    ,DRAD2   ,
     J        INTFRIC ,IPARTFRICSI,
     K        IPARTFRICMI  ,ITIED   ,JLT_TIED,INTBUF_TAB%CAND_F,IORTHFRIC,
     L        IREP_FRICMI  ,DIR_FRICMI,DGAPLOADP)

            JLT = JLT_NEW
            IF (IMONM > 0 .AND. JTASK == 1) CALL STARTIME(TIMERS,20)
            IF(JLT_NEW /= 0) THEN
              IPARI(29,NIN) = 1
              IF (DEBUG(3)>=1) NB_JLT_NEW = NB_JLT_NEW + JLT_NEW

        !-------------------------------------------------------------------------------
        ! Friction model : computation of friction coefficients based on Material of connected Parts
        !-------------------------------------------------------------------------------
        IF(ITASK==1) CALL STARTIME(TIMERS,MACRO_TIMER_FRIC)
           IF(IORTHFRIC > 0) THEN
              CALL FRICTIONPARTS_MODEL_ORTHO(
     1           INTFRIC        ,JLT          ,IPARTFRICSI  ,IPARTFRICMI  ,ADPARTS_FRIC   ,
     2           NSETPRTS       ,TABCOUPLEPARTS_FRIC,NPARTFRIC,TABPARTS_FRIC,TABCOEF_FRIC ,
     3           FRIC           ,VISCF        ,INTBUF_TAB%FRIC_P,FRIC_COEFS   , FRICC     ,
     4           VISCFFRIC      ,NTY          ,MFROT            ,IORTHFRIC    , FRIC_COEFS2,
     5           FRICC2         ,VISCFFRIC2   ,IFRICORTH        ,NFORTH       , NFISOT    ,
     6           INDEXORTH      ,INDEXISOT    ,JLT_TIED         ,IREP_FRICMI  , DIR_FRICMI,         
     7           IX3            ,IX4           ,X1              ,Y1            , Z1        , 
     8           X2             ,Y2            ,Z2              ,X3            , Y3        ,
     9           Z3             ,X4            ,Y4              ,Z4            ,CE_LOC     ,
     A           DIR1           ,DIR2          )
            ELSE
              NFORTH = 0
              NFISOT = 0
              CALL FRICTIONPARTS_MODEL_ISOT(
     1           INTFRIC        ,JLT          ,IPARTFRICSI  ,IPARTFRICMI  ,ADPARTS_FRIC   ,
     2           NSETPRTS       ,TABCOUPLEPARTS_FRIC,NPARTFRIC,TABPARTS_FRIC,TABCOEF_FRIC ,
     3           FRIC           ,VISCF        ,INTBUF_TAB%FRIC_P,FRIC_COEFS   , FRICC     ,
     4           VISCFFRIC      ,NTY          ,MFROT            ,IORTHFRIC    , IFRIC     ,
     5           JLT_TIED       , TINT        ,TEMPI            ,NPC          ,TF         ,
     6           TEMP           , H1          ,H2               ,H3           ,H4         ,
     7           IX1            , IX2         ,IX3              ,IX4          ,IFORM      )
            ENDIF
            IF(ITASK==1) CALL STOPTIME(TIMERS,MACRO_TIMER_FRIC)

              CALL I7FOR3(OUTPUT,
     1  JLT          ,A         ,V            ,IBC         ,ICODT    ,
     2  FSAV         ,GAP       ,FRIC         ,MS          ,VISC     ,
     3  VISCF        ,NOINT     ,INTBUF_TAB%STFNS,ITAB     ,CN_LOC   ,
     4  STIGLO       ,STIFN     ,STIF         ,FSKYI       ,ISKY     ,
     5  NX1          ,NX2       ,NX3          ,NX4         ,NY1      ,
     6  NY2          ,NY3       ,NY4          ,NZ1         ,NZ2      ,
     7  NZ3          ,NZ4       ,LB1          ,LB2         ,LB3      ,
     8  LB4          ,LC1       ,LC2          ,LC3         ,LC4      ,
     9  P1           ,P2        ,P3           ,P4          ,FCONT    ,
     A  IX1          ,IX2       ,IX3          ,IX4         ,NSVG     ,
     B  IVIS2        ,NELTST    ,ITYPTST      ,DT2T        ,GAPV     ,
     C  INACTI,INTBUF_TAB%CAND_P,INDEX2(NFT+1),KINET       ,NEWFRONT  ,
     D  ISECIN       ,NSTRF     ,X       ,INTBUF_TAB%IRECTM,CE_LOC    ,
     E  MFROT        ,IFQ       ,INTBUF_TAB%FTSAVX  ,INTBUF_TAB%FTSAVY,
     F  INTBUF_TAB%FTSAVZ,XFILTR_FRIC,INTBUF_TAB%IFPEN     ,IBAG     ,
     H  ICONTACT     ,VISCN     ,VXI          ,VYI         ,VZI      ,    
     I  MSI         ,KINI       ,NIN          ,NISUB   ,INTBUF_TAB%LISUB,
     J  INTBUF_TAB%ADDSUBS,INTBUF_TAB%ADDSUBM,INTBUF_TAB%LISUBS     ,
     .                                     INTBUF_TAB%LISUBM,FSAVSUB,
     K  INTBUF_TAB%CAND_N,IPARI(33,NIN),IPARI(39,NIN),INTBUF_TAB%NODNORM_NORMAL ,FNCONT,     
     L  FTCONT       ,X1         ,X2          ,X3         ,X4         ,
     M  Y1           , Y2        ,Y3          ,Y4         ,Z1         ,
     N  Z2           ,Z3         ,Z4          ,XI         ,YI         ,
     O  ZI           ,IADM       ,RCURVI      ,RCONTACT   ,ACONTACT   ,
     P  PCONTACT     ,ANGLMI     ,PADM        ,INTTH      ,TEMP       , 
     Q  TEMPI        ,IFORM      ,NPC         ,TF         ,CMAJ       ,
     R  DTMINI       ,DRAD       ,FHEATS      ,FHEATM     ,EFRICT     ,
     S  QFRICINT(NIN),FNI        ,IFRIC       ,JTASK      ,H1         ,
     U  H2           ,H3         ,H4          ,KS         ,KT         ,
     V  K1           ,K2         ,K3          ,K4         ,C1         ,
     W  C2           ,C3         ,C4          ,CS         ,C          ,
     X  CF           ,TINT       ,XFRIC       ,FXI        ,FYI        ,
     Y  FZI          ,FX1        ,FY1         ,FZ1        ,FX2        ,
     Z  FY2          ,FZ2        ,FX3         ,FY3        ,FZ3        ,
     1  FX4          ,FY4        ,FZ4         ,ISENSINT    ,FSAVPARIT  ,
     5  NFT          ,SYM_FLAG_TYPE19,H3D_DATA,FRICC      ,VISCFFRIC  , 
     6  FRIC_COEFS   ,ITIED      ,JLT_TIED   ,INTBUF_TAB%CAND_F,IORTHFRIC,
     7  FRIC_COEFS2  ,FRICC2     ,VISCFFRIC2  ,NFORTH      ,NFISOT     ,
     8  INDEXORTH    ,INDEXISOT  ,DIR1        ,DIR2        ,TAGNCONT   ,
     9  KLOADPINTER  ,LOADPINTER ,LOADP_HYD_INTER,INTBUF_TAB%TYPSUB    ,  
     .                                            INTBUF_TAB%INFLG_SUBS,
     A  INTBUF_TAB%INFLG_SUBM,NINLOADP,DGAPLOADINT,S_LOADPINTER,DGAPLOADP,
     B  INTEREFRIC   )

              IF( INTTH > 0 ) THEN 
                 CALL I7THERM(JLT   ,IPARG ,PM  ,IXS   ,IFORM ,X      ,   
     1                        XI    ,YI    ,ZI  ,X1   ,Y1    ,Z1      ,
     2                        X2    ,Y2    ,Z2  ,X3   ,Y3    ,Z3      ,
     3                        X4    ,Y4    ,Z4  ,IX1  ,IX2   ,IX3     ,
     4                        IX4 ,RSTIF  ,TEMPI, INTBUF_TAB%IELEC    ,
     5                        PHI   ,TINT , AREASI,IELECI  ,FRAD,DRAD ,
     6                        GAPV  ,FNI  ,IRSTH,XTHE  ,NPC    ,TF    ,
     6                        CONDINT,PHI1,PHI2  ,PHI3   ,PHI4  ,FHEATS,
     7                        FHEATM,EFRICT,TEMP ,H1     ,H2     ,H3   ,
     8                        H4    ,theaccfact)      
              ENDIF

      CALL I7ASS33(
     1  JLT       ,A         ,NOINT     ,ITAB      ,STIFN     ,
     2 STIF      ,FSKYI      ,ISKY      ,FCONT     ,IX1       ,
     3  IX2       ,IX3       ,IX4       ,NSVG      ,NELTST    ,
     4  ITYPTST   ,DT2T      ,NISKYFI   ,ISECIN    ,NSTRF     ,
     5  SECFCUM   ,VISCN     ,NIN       ,FXI       ,FYI       , 
     6  FZI       ,FX1       ,FY1       ,FZ1       ,FX2       ,
     7  FY2       ,FZ2       ,FX3       ,FY3       ,FZ3       ,
     8  FX4       ,FY4       ,FZ4       ,H1        ,H2        ,
     9  H3        ,H4        ,KS        ,KT        ,K1        ,
     A  K2        ,K3        ,K4        ,CS        ,CF        ,
     B  C1        ,C2        ,C3        ,C4        ,C         ,
     C  INTTH     ,PHI       ,PHI1      ,PHI2      ,PHI3      ,
     D  PHI4      ,FTHE      ,FTHESKYI  ,MSKYI_SMS ,ISKYI_SMS ,
     E  NSMS      ,DTMINI    ,JTASK     ,
     F  CONDN     ,CONDINT   ,CONDNSKYI ,IXIG3D    ,KXIG3D    ,
     J  WIGE      ,KNOT      ,IGEO      ,INTBUF_TAB%NIGE      ,
     K  INTBUF_TAB%RIGE      ,X         ,H3D_DATA  ,KNOTLOCPC ,
     L  KNOTLOCEL ,IFORM     ,NODADT_THERM) 

           ENDIF !(JLT_NEW/=0) THEN
      ELSE !STS_ACTIVE = .TRUE.
C==========================================================
C NEW STS ALGORITHM - UNIFIED FOR GAUSS AND LOBATTO
C==========================================================
C


        ! Find contact partners/candidates only once
        CALL REMAP_SEGMENTS_STS(INTBUF_TAB, ITAB, X, CAND_SEC_SEG,
     &   INTBUF_TAB%IRECTM, CONT_ELEMENT, COUNT, IGRSURF,
     &   CAND_SEC_SEG_ID, CAND_MST_SEG_ID, MAX_STS_SIZE_ACTUAL)

        WRITE(*,*) "---------------------------------------------------"
        WRITE(*,*) "STS: Found ", COUNT, " potentional contact element pairs"
        WRITE(*,*)

        !-------------------------------------------------------------------------------
        ! Friction model : computation of friction coefficients based on Material of connected Parts
        !-------------------------------------------------------------------------------
        IF(ITASK==1) CALL STARTIME(TIMERS,MACRO_TIMER_FRIC)

        IF (COUNT > 0) THEN
          NFORTH = 0
          NFISOT = 0
          ! Map STS contact pairs to friction arrays
          ! For now, use simplified approach: map each contact pair
          ! TODO: Properly map part IDs from CAND_SEC_SEG_ID and CAND_MST_SEG_ID
          DO I = 1, MIN(COUNT, MVSIZ)
            ! Initialize friction arrays for this contact pair
            ! This is a placeholder - proper part ID mapping needed
            IPARTFRICSI(I) = 0
            IPARTFRICMI(I) = 0
          ENDDO
          
          ! Call friction model (isotropic)
          ! It will use default values from FRIC and INTBUF_TAB%FRIC_P

          IF (IORTHFRIC == 0) THEN
            CALL FRICTIONPARTS_MODEL_ISOT(
     1           INTFRIC        ,MIN(COUNT,MVSIZ),IPARTFRICSI  ,IPARTFRICMI  ,ADPARTS_FRIC   ,
     2           NSETPRTS       ,TABCOUPLEPARTS_FRIC,NPARTFRIC,TABPARTS_FRIC,TABCOEF_FRIC ,
     3           FRIC           ,VISCF        ,INTBUF_TAB%FRIC_P,FRIC_COEFS   , FRICC     ,
     4           VISCFFRIC      ,NTY          ,MFROT            ,IORTHFRIC    , IFRIC     ,
     5           JLT_TIED       , TINT        ,TEMPI            ,NPC          ,TF         ,
     6           TEMP           , H1          ,H2               ,H3           ,H4         ,
     7           IX1            , IX2         ,IX3              ,IX4          ,IFORM      )

          ENDIF
        ENDIF

        IF(ITASK==1) CALL STOPTIME(TIMERS,MACRO_TIMER_FRIC)

        ! -----------------------
        ! Try GAUSS SEARCH ALGORITHM
        ! -----------------------
        WRITE(*,*) "***********************************"
        WRITE(*,*) "STS - GAUSS SEARCH ALGORITHM STARTED"
        WRITE(*,*) "***********************************"

        OPTION = 0 !FOR GAUSS STS ALGORITHM
        QFRICT = 0.d0 !Friction energy
        ALPHA0 = XFILTR_FRIC !Parameter for friction force calculation
               
        CALL PROCESS_STS_CONTACTS(CONT_ELEMENT, COUNT, OPTION,
     &   CAND_MST_SEG_ID, CAND_SEC_SEG_ID,
     &   load_arr, node_id_load, L, IMPACT_glob, STIF, MAX_STS_SIZE_ACTUAL,
     &   FRICC, FRIC_COEFS, VISCFFRIC, XMU, MFROT, IFQ, DT1, DT12,
     &   V, MS, INTBUF_TAB%CAND_F, ALPHA0, INTBUF_TAB%IFPEN, INTTH, QFRICT,
     &   INTBUF_TAB, GAP, XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST,
     &   GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     &   GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     &   GP_XI1_PERIOD, GP_XI2_PERIOD,
     &   GP_INITIALIZED,
     &   GP_TTRIAL1_HIST, GP_TTRIAL2_HIST,
     &   GP_IS_STICKING)
        !IMPACT_glob indicates if any contact was found

        ! If no contact found with GAUSS algorithm, try LOBATTO algorithm
        IF (IMPACT_glob == 0) THEN
            WRITE(*,*) "NO IMPACT FOUND WITH GAUSS SEARCH ALGORITHM"

            ! -----------------------
            ! Try LOBATTO SEARCH ALGORITHM
            ! -----------------------
            WRITE(*,*) "***********************************"
            WRITE(*,*) "STS - LOBATTO SEARCH ALGORITHM STARTED"
            WRITE(*,*) "***********************************"

            OPTION = 1 !FOR LOBATTO STS ALGORITHM
            QFRICT = 0.d0
            ALPHA0 = XFILTR_FRIC
            CALL PROCESS_STS_CONTACTS(CONT_ELEMENT, COUNT, OPTION,
     &      CAND_MST_SEG_ID, CAND_SEC_SEG_ID,
     &      load_arr, node_id_load, L, IMPACT_glob, STIF, MAX_STS_SIZE_ACTUAL,
     &      FRICC, FRIC_COEFS, VISCFFRIC, XMU, MFROT, IFQ, DT1, DT12,
     &      V, MS, INTBUF_TAB%CAND_F, ALPHA0, INTBUF_TAB%IFPEN, INTTH, QFRICT,
     &      INTBUF_TAB, GAP, XI1_HIST, XI2_HIST, TTRIAL1_HIST, TTRIAL2_HIST,
     &      GP_XI1_GLOBAL, GP_XI2_GLOBAL,
     &      GP_XI1_GLOBAL_PREV, GP_XI2_GLOBAL_PREV,
     &      GP_XI1_PERIOD, GP_XI2_PERIOD,
     &      GP_INITIALIZED,
     &      GP_TTRIAL1_HIST, GP_TTRIAL2_HIST,
     &      GP_IS_STICKING)
          
            IF (IMPACT_glob == 0) THEN
                WRITE(*,*) "NO IMPACT FOUND WITH LOBATTO SEARCH ALGORITHM"
            ENDIF  
        ENDIF

        ! Assemble forces if contact was found with either algorithm  
        IF (IMPACT_glob /= 0) THEN
          CALL I7ASS33(
     1  JLT       ,A         ,NOINT     ,ITAB      ,STIFN     ,
     2  STIF      ,FSKYI      ,ISKY      ,FCONT     ,IX1       ,
     3  IX2       ,IX3       ,IX4       ,NSVG      ,NELTST    ,
     4  ITYPTST   ,DT2T      ,NISKYFI   ,ISECIN    ,NSTRF     ,
     5  SECFCUM   ,VISCN     ,NIN       ,FXI       ,FYI       , 
     6  FZI       ,FX1       ,FY1       ,FZ1       ,FX2       ,
     7  FY2       ,FZ2       ,FX3       ,FY3       ,FZ3       ,
     8  FX4       ,FY4       ,FZ4       ,H1        ,H2        ,
     9  H3        ,H4        ,KS        ,KT        ,K1        ,
     A  K2        ,K3        ,K4        ,CS        ,CF        ,
     B  C1        ,C2        ,C3        ,C4        ,C         ,
     C  INTTH     ,PHI       ,PHI1      ,PHI2      ,PHI3      ,
     D  PHI4      ,FTHE      ,FTHESKYI  ,MSKYI_SMS ,ISKYI_SMS ,
     E  NSMS      ,DTMINI    ,JTASK     ,
     F  CONDN     ,CONDINT   ,CONDNSKYI ,IXIG3D    ,KXIG3D    ,
     J  WIGE      ,KNOT      ,IGEO      ,INTBUF_TAB%NIGE      ,
     K  INTBUF_TAB%RIGE      ,X         ,H3D_DATA  ,KNOTLOCPC ,
     L  KNOTLOCEL ,IFORM     ,NODADT_THERM, ASS_FLAG , load_arr, node_id_load, L, MAX_STS_SIZE_ACTUAL ) 
        ENDIF
        
        ENDIF ! STS_ACTIVE = .TRUE.
        ENDIF ! TYPE18
C-----------------------------------------------
C
          IF (IMONM > 0 .AND. JTASK == 1) CALL STOPTIME(TIMERS,20)
          IF(IMPL_S = =1) THEN
            DO I = 1 ,JLT_NEW
             NS_IMP(I+NUM_IMP)=CN_LOC(I)
             NE_IMP(I+NUM_IMP)=CE_LOC(I)
             IND_IMP(I+NUM_IMP)=INDEX2(I+NFT)
            ENDDO
            NUM_IMP=NUM_IMP+JLT_NEW
          ENDIF
        ENDDO
C!$OMP END ORDERED
        IF (SFSAVPARIT /= 0)THEN
          CALL SUM_6_FLOAT_SENS(FSAVPARIT, NISUB+1, 11, I_STOK,1,I_STOK, FBSAV6, 12, 6, DIMFB, ISENSINT )
        ENDIF
        IF (ALLOCATED(FSAVPARIT)) DEALLOCATE (FSAVPARIT)

        IF(ICURV == 3 .OR. IADM /= 0)THEN
          CALL MY_BARRIER
          IF(JTASK==1)DEALLOCATE(INTBUF_TAB%NODNORM_NORMAL)
        END IF
        IF(IADM/=0)THEN
          CALL MY_BARRIER
          IF(JTASK==1)DEALLOCATE(INTBUF_TAB%MODRCURV,INTBUF_TAB%MODANGLM)
        END IF
      ELSEIF(NTY==17)THEN
        ILAGM = IPARI(33,NIN)
C       I7KGLO = 1 => mis ds inttri
        IF(ILAGM == 0)THEN
          CALL I17MAIN_PENA(output,
     1          NIN    ,IPARI  ,INTBUF_TAB    ,X           ,V  ,
     2          A      ,JTASK-1,IGRBRIC ,EMINX(IAD17),MS ,
     3          IXS    ,IXS16  ,IXS20   ,STIFN  ,FSKYI  ,ISKY    ,
     4          FSAV   ,FCONT  ,NISKYFI ,H3D_DATA)
        ENDIF
      ENDIF
C
      RETURN
      END

