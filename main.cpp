#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include "src/model.hpp"
#include "src/functions.hpp"
#include "src/equilibrium.hpp"
#include "src/fluid_prop.hpp"
#include "src/permeability.hpp"
#include "src/transmissibility.hpp"
#include "src/pertubation.hpp"
#include "src/slope-calc.hpp"
#include "src/accumulation.hpp"
#include "src/transport.hpp"
#include "src/well.hpp"
#include "src/solver.hpp"

int main (void)
{

    gsl_vector * SO = gsl_vector_alloc(GMIJK+2); gsl_vector * SG = gsl_vector_alloc(GMIJK+2); gsl_vector * SW = gsl_vector_alloc(GMIJK+2);
    gsl_vector * SO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * SG_V = gsl_vector_alloc(GMIJK+2); gsl_vector * SW_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * PO = gsl_vector_alloc(GMIJK+2); gsl_vector * PG = gsl_vector_alloc(GMIJK+2); gsl_vector * PW = gsl_vector_alloc(GMIJK+2);
    gsl_vector * PO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * PG_V = gsl_vector_alloc(GMIJK+2); gsl_vector * PW_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * PCGO = gsl_vector_alloc(GMIJK+2); gsl_vector * PCOW = gsl_vector_alloc(GMIJK+2); gsl_vector * ZB = gsl_vector_alloc(GMIJK+2);
    gsl_vector * PCGO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * PCOW_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * PB = gsl_vector_alloc(GMIJK+2); gsl_vector * PB_V = gsl_vector_alloc(GMIJK+2);

    gsl_vector * BO = gsl_vector_alloc(GMIJK+2); gsl_vector * BG = gsl_vector_alloc(GMIJK+2); gsl_vector * BW = gsl_vector_alloc(GMIJK+2); gsl_vector * RS = gsl_vector_alloc(GMIJK+2);
    gsl_vector * RHO = gsl_vector_alloc(GMIJK+2); gsl_vector * RHG = gsl_vector_alloc(GMIJK+2); gsl_vector * RHW = gsl_vector_alloc(GMIJK+2);
    gsl_vector * VISO = gsl_vector_alloc(GMIJK+2); gsl_vector * VISG = gsl_vector_alloc(GMIJK+2); gsl_vector * BPOR = gsl_vector_alloc(GMIJK+2);
    gsl_vector * KI = gsl_vector_alloc(GMIJK+2); gsl_vector * KJ = gsl_vector_alloc(GMIJK+2); gsl_vector * KK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * KRO = gsl_vector_alloc(GMIJK+2); gsl_vector * KRG = gsl_vector_alloc(GMIJK+2); gsl_vector * KRW = gsl_vector_alloc(GMIJK+2);
    gsl_vector * KRO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * KRG_V = gsl_vector_alloc(GMIJK+2); gsl_vector * KRW_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * SLPOR = gsl_vector_alloc(GMIJK+2); gsl_vector * SLBO = gsl_vector_alloc(GMIJK+2); gsl_vector * SLBG = gsl_vector_alloc(GMIJK+2); gsl_vector * SLBW = gsl_vector_alloc(GMIJK+2); gsl_vector * SLRS = gsl_vector_alloc(GMIJK+2);
    gsl_vector * BO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * BG_V = gsl_vector_alloc(GMIJK+2); gsl_vector * BW_V = gsl_vector_alloc(GMIJK+2); gsl_vector * RS_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * RHO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * RHG_V = gsl_vector_alloc(GMIJK+2); gsl_vector * RHW_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * VISO_V = gsl_vector_alloc(GMIJK+2); gsl_vector * VISG_V = gsl_vector_alloc(GMIJK+2); gsl_vector * BPOR_V = gsl_vector_alloc(GMIJK+2);

    gsl_vector * TGI = gsl_vector_alloc(GMIJK+2); gsl_vector * TGJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TGK = gsl_vector_alloc(GMIJK+2);

    gsl_vector * TFP_OI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_OJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_OK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFP_GI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_GJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_GK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFP_WI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_WJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_WK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFP_ORSI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_ORSJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_ORSK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFS_OI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_OJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_OK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFS_GI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_GJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_GK = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFS_WI = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_WJ = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_WK = gsl_vector_alloc(GMIJK+2);

    gsl_vector * TFP_OI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_OJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_OK_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFP_GI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_GJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_GK_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFP_WI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_WJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_WK_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFP_ORSI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_ORSJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFP_ORSK_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFS_OI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_OJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_OK_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFS_GI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_GJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_GK_V = gsl_vector_alloc(GMIJK+2);
    gsl_vector * TFS_WI_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_WJ_V = gsl_vector_alloc(GMIJK+2); gsl_vector * TFS_WK_V = gsl_vector_alloc(GMIJK+2);

    gsl_vector * CN_WGO = gsl_vector_alloc(GMIJK*3);
    gsl_vector * CV_WGO = gsl_vector_alloc(GMIJK*3);
    gsl_matrix * CM_WGO =  gsl_matrix_alloc(GMIJK*3,GMIJK*3);

    gsl_vector * QV = gsl_vector_alloc(GMIJK*3);
    gsl_matrix * QRN =  gsl_matrix_alloc(GMIJK*3,GMIJK*3);

    gsl_vector * FNV = gsl_vector_alloc(GMIJK*3);
    gsl_matrix * JACOBIAN = gsl_matrix_alloc(GMIJK*3,GMIJK*3);

    gsl_vector * X_N = gsl_vector_alloc(GMIJK*3);
    gsl_vector * X_NM = gsl_vector_alloc(GMIJK*3);
    gsl_vector * D_NV = gsl_vector_alloc(GMIJK*3);
    gsl_vector * X_VP = gsl_vector_alloc(GMIJK*3);
    gsl_vector * X_V = gsl_vector_alloc(GMIJK*3);
    gsl_vector * residual = gsl_vector_alloc(1);
    gsl_matrix * RESULT = gsl_matrix_alloc(GMIJK*3,30);

    EQUILIBRIUM(Gi, Gj, Gk, Ks, DWOC, DGOC, Dk,
                REFDP, REFPRES, PRBB, TSW, TSL, TPCGO, TPCOW, SO, SG, SW,
                PO, PG, PW, PB, PCGO, PCOW, ZB);

    PROP_2_VEC_WGP(Gi, Gj, Gk, SW, SG, PO, X_N);
    gsl_matrix_set_col(RESULT,0,X_N);

    PERM_DIST(Gi, Gj, Gk, PERMI, PERMJ, PERMK, KI, KJ, KK);

    GEO_REC_CTER_BLCK(Gi, Gj, Gk, Li, Lj, Lk, KI, KJ, KK,
                      TGI, TGJ, TGK);

    /** inciar/chute */
    int Tsz = 1;
    double sss;

    double DIF_G_SW = 0.03;
    double DIF_G_SG = 0.005;
    double DIF_G_PO = -100;

    /** Generate Guess for t = n+1 Props*/
    PERTUB_SAT(Gi, Gj, Gk,
               DIF_G_SW, DIF_G_SW, DIF_G_SW, DIF_G_SW, DIF_G_SW,
               DIF_G_SG, DIF_G_SG, DIF_G_SG, DIF_G_SG, DIF_G_SG,
               SO, SG, SW,
               SO_V, SG_V, SW_V);

    PERTUB_PRES(Gi, Gj, Gk,
                PRBB, DIF_G_PO,
                SG, PO, PG, PW, PB, PCGO, PCOW,
                PO_V, PG_V, PW_V, PB_V, PCGO_V, PCOW_V);

    do { /** Time step */
        printf("\n\nInicio Time: %d\n\n",Tsz);
        /** Hold t = n, Props */
        gsl_vector_memcpy(X_NM, X_N);

        PROP_2_VEC_WGP(Gi, Gj, Gk, SW_V, SG_V, PO_V, X_V);
        gsl_vector_memcpy(D_NV, X_N);
        gsl_vector_sub(D_NV, X_V);
        /** Prop t = n */
        FLUID_PROP(Gi, Gj, Gk, Dk, PRBB,
                   SO, SG, SW, PO, PG, PW, PB,
                   BO, BG, BW, RS, VISO, VISG,
                   RHO, RHG, RHW, BPOR);

        REL_PERM_DIST(Gi, Gj, Gk,
                      SO, SG, SW,
                      KRO, KRG, KRW);

        TRANS_1PW_PRESS(Gi, Gj, Gk, PO, PG, PW, BO, BG, BW, RS, VISO, VISG, VISW,
                    TFP_OI, TFP_OJ, TFP_OK,
                    TFP_GI, TFP_GJ, TFP_GK,
                    TFP_WI, TFP_WJ, TFP_WK,
                    TFP_ORSI, TFP_ORSJ, TFP_ORSK);

        TRANS_1PW_SAT(Gi, Gj, Gk, KRO, KRG, KRW,
                  TFS_OI, TFS_OJ, TFS_OK,
                  TFS_GI, TFS_GJ, TFS_GK,
                  TFS_WI, TFS_WJ, TFS_WK);
        /** end Prop t = n */

        /** Prop t = n+1 */
        /** Guess Press for t = n+1 */
        if ( Tsz > 1 ) {
            DIF_G_SW = 0.03;
            DIF_G_SG = 0.005;
            DIF_G_PO = -55;

            PERTUB_SAT(Gi, Gj, Gk,
                       DIF_G_SW, DIF_G_SW, DIF_G_SW, DIF_G_SW, DIF_G_SW,
                       DIF_G_SG, DIF_G_SG, DIF_G_SG, DIF_G_SG, DIF_G_SG,
                       SO, SG, SW,
                       SO_V, SG_V, SW_V);

            PERTUB_PRES(Gi, Gj, Gk,
                        PRBB, DIF_G_PO,
                        SG, PO, PG, PW, PB, PCGO, PCOW,
                        PO_V, PG_V, PW_V, PB_V, PCGO_V, PCOW_V);

        }

        FLUID_PROP(Gi, Gj, Gk, Dk, PRBB,
                   SO_V, SG_V, SW_V, PO_V, PG_V, PW_V, PB_V,
                   BO_V, BG_V, BW_V, RS_V, VISO_V, VISG_V,
                   RHO_V, RHG_V, RHW_V, BPOR_V);

        REL_PERM_DIST(Gi, Gj, Gk,
                      SO_V, SG_V, SW_V,
                      KRO_V, KRG_V, KRW_V);

        TRANS_1PW_PRESS(Gi, Gj, Gk, PO_V, PG_V, PW_V, BO_V, BG_V, BW_V, RS_V, VISO_V, VISG_V, VISW,
                    TFP_OI_V, TFP_OJ_V, TFP_OK_V,
                    TFP_GI_V, TFP_GJ_V, TFP_GK_V,
                    TFP_WI_V, TFP_WJ_V, TFP_WK_V,
                    TFP_ORSI_V, TFP_ORSJ_V, TFP_ORSK_V);

        TRANS_1PW_SAT(Gi, Gj, Gk, KRO_V, KRG_V, KRW_V,
                  TFS_OI_V, TFS_OJ_V, TFS_OK_V,
                  TFS_GI_V, TFS_GJ_V, TFS_GK_V,
                  TFS_WI_V, TFS_WJ_V, TFS_WK_V);
        /** end Prop t = n+1 */
        /** Calc Slope Accumulation */
        SLOPE_CALC(Gi, Gj, Gk, ZB,
                   PO, PO_V, PG, PG_V, PW, PW_V,
                   BPOR, BPOR_V, SLPOR,
                   BO, BO_V, SLBO,
                   BG, BG_V, SLBG,
                   BW, BW_V, SLBW,
                   RS, RS_V, SLRS);
        /** Calc Jacobian Transport */
        TRANSPORT_1PW_2(Gi, Gj, Gk, D_NV, ZB, SO_V, SW_V, SG_V, PO_V, PW_V, PG_V, PCGO_V, PCOW_V,
                        RHO, RHG, RHW, RHO_V, RHG_V, RHW_V, RS, RS_V,
                        TGI,  TGJ,  TGK,
                     TFS_OI,  TFS_OJ,  TFS_OK, TFS_GI,  TFS_GJ,  TFS_GK, TFS_WI,  TFS_WJ,  TFS_WK,
                     TFP_OI,  TFP_OJ,  TFP_OK, TFP_GI,  TFP_GJ,  TFP_GK, TFP_WI,  TFP_WJ,  TFP_WK,
                     TFP_ORSI,  TFP_ORSJ,  TFP_ORSK,
                     TFS_OI_V,  TFS_OJ_V,  TFS_OK_V, TFS_GI_V,  TFS_GJ_V,  TFS_GK_V, TFS_WI_V,  TFS_WJ_V,  TFS_WK_V,
                     TFP_OI_V,  TFP_OJ_V,  TFP_OK_V, TFP_GI_V,  TFP_GJ_V,  TFP_GK_V, TFP_WI_V,  TFP_WJ_V,  TFP_WK_V,
                     TFP_ORSI_V,  TFP_ORSJ_V,  TFP_ORSK_V,
                     FNV, JACOBIAN);
        /** Accumulation Vector t = n*/
        ACCUMULATION_RESID(1, Gi, Gj, Gk, Is, Js, Ks,
                           SO, SG, SW, BPOR, BO, BG, BW, RS, CN_WGO);
        /** Accumulation Vector t = n+1*/
        ACCUMULATION_RESID(1, Gi, Gj, Gk, Is, Js, Ks,
                           SO_V, SG_V, SW_V, BPOR_V, BO_V, BG_V, BW_V, RS_V, CV_WGO);
        /** Calc Accumulation Matrix vt = (n+1) - n */
        ACCUMULATION_DIF(1, Gi, Gj, Gk, Is, Js, Ks,
                          SO, SO_V, SG, SG_V, SW, SW_V, BPOR, BPOR_V, SLPOR, BO, BO_V, SLBO, BG, BG_V, SLBG,
                          BW, BW_V, SLBW, RS, RS_V, SLRS, CM_WGO);
        /** Calc Well Stats */
        WELL_PROD(Gi, Gj, Gk, 2, 2, 2, 1.1300643, 5250, 10000, D_NV,
                   PO, SW, SG, KRO, KRG, KRW, BO, BG, BW, RS, VISO, VISG,
                   PO_V, SW_V, SG_V, KRO_V, KRG_V, KRW_V, BO_V, BG_V, BW_V, RS_V, VISO_V, VISG_V,
                   VISW, QV, QRN);
        /** Solve [T(n,v)]X(v)=[A(n,v)] */
        SOLVER(Gi, Gj, Gk, JACOBIAN, CM_WGO, QRN, FNV, CN_WGO, CV_WGO, QV, X_V, residual);
        /** t = n+2 => t = n+1 = t = n, Props to X_N */
        PROP_2_VEC_WGP(Gi, Gj, Gk, SW_V, SG_V, PO_V, X_N);

        for (size_t i=3; i<4; i++) {
            printf("\nTime = n, Props - Time: %d\n",Tsz);
            printf("Zg  | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i),gsl_vector_get(SG,i),gsl_vector_get(SO,i),gsl_vector_get(PO,i));
            printf("Zgo | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+100),gsl_vector_get(SG,i+100),gsl_vector_get(SO,i+100),gsl_vector_get(PO,i+100));
            printf("Zo  | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+200),gsl_vector_get(SG,i+200),gsl_vector_get(SO,i+200),gsl_vector_get(PO,i+200));
            printf("Zow | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+300),gsl_vector_get(SG,i+300),gsl_vector_get(SO,i+300),gsl_vector_get(PO,i+300));
            printf("Zw  | SW: %f | SG: %f | SO: %f | PO: %f\n\n",gsl_vector_get(SW,i+400),gsl_vector_get(SG,i+400),gsl_vector_get(SO,i+400),gsl_vector_get(PO,i+400));
        }

        for (size_t i=3; i<4; i++) {
            printf("Time: n+1, Props - Time: %d\n",Tsz);
            printf("Zg  | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i),gsl_vector_get(SG_V,i),gsl_vector_get(SO_V,i),gsl_vector_get(PO_V,i));
            printf("Zgo | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i+100),gsl_vector_get(SG_V,i+100),gsl_vector_get(SO_V,i+100),gsl_vector_get(PO_V,i+100));
            printf("Zo  | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i+200),gsl_vector_get(SG_V,i+200),gsl_vector_get(SO_V,i+200),gsl_vector_get(PO_V,i+200));
            printf("Zow | SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i+300),gsl_vector_get(SG_V,i+300),gsl_vector_get(SO_V,i+300),gsl_vector_get(PO_V,i+300));
            printf("Zw  | SW: %f | SG: %f | SO: %f | PO: %f\n\n",gsl_vector_get(SW_V,i+400),gsl_vector_get(SG_V,i+400),gsl_vector_get(SO_V,i+400),gsl_vector_get(PO_V,i+400));
        }

        size_t dn = 0;
        do { /** Newton Iter */
            std::cout << std::endl;
            std::cout << "Inicio NR Iters #" << dn << std::endl;
            std::cout << std::endl;
            double rck = gsl_vector_get(residual,0);
            /** Get New Xv to new iter */
            gsl_vector_add(X_V, X_N);
            gsl_vector_memcpy(X_VP, X_V); /// Save for nex vT

            /** Get change in Props */
            gsl_vector_memcpy(D_NV, X_NM);
            gsl_vector_sub(D_NV, X_V);

            SRT_VEC_2_OGW(PRBB, SW_V, SG_V, PO_V, X_V, SO_V, SG_V, SW_V, PO_V, PG_V, PW_V, PB_V);
            //SRT_VEC_2_OGW(PRBB, SW, SG, PO, X_N, SO, SG, SW, PO, PG, PW, PB);

            for (size_t i=2; i<3; i++) {
                printf("Vec NR - Time: %d\n",Tsz);
                printf("SW: %f | SG: %f | PO: %f -|-\n",gsl_vector_get(X_V,i*3),gsl_vector_get(X_V,((i*3)+1)),gsl_vector_get(X_V,((i*3)+2)));
                printf("SW: %f | SG: %f | PO: %f -|-\n",gsl_vector_get(X_V,((i*3)+300)),gsl_vector_get(X_V,((i*3)+301)),gsl_vector_get(X_V,((i*3)+302)));
                printf("SW: %f | SG: %f | PO: %f -|-\n",gsl_vector_get(X_V,((i*3)+600)),gsl_vector_get(X_V,((i*3)+601)),gsl_vector_get(X_V,((i*3)+602)));
                printf("SW: %f | SG: %f | PO: %f -|-\n",gsl_vector_get(X_V,((i*3)+900)),gsl_vector_get(X_V,((i*3)+901)),gsl_vector_get(X_V,((i*3)+902)));
                printf("SW: %f | SG: %f | PO: %f\n\n",gsl_vector_get(X_V,((i*3)+1200)),gsl_vector_get(X_V,((i*3)+1201)),gsl_vector_get(X_V,((i*3)+1202)));
            }

            for (size_t i=3; i<4; i++) {
                printf("Vec Run - Time: %d\n",Tsz);
                printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i),gsl_vector_get(SG_V,i),gsl_vector_get(SO_V,i),gsl_vector_get(PO_V,i));
                printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i+100),gsl_vector_get(SG_V,i+100),gsl_vector_get(SO_V,i+100),gsl_vector_get(PO_V,i+100));
                printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i+200),gsl_vector_get(SG_V,i+200),gsl_vector_get(SO_V,i+200),gsl_vector_get(PO_V,i+200));
                printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW_V,i+300),gsl_vector_get(SG_V,i+300),gsl_vector_get(SO_V,i+300),gsl_vector_get(PO_V,i+300));
                printf("SW: %f | SG: %f | SO: %f | PO: %f\n\n",gsl_vector_get(SW_V,i+400),gsl_vector_get(SG_V,i+400),gsl_vector_get(SO_V,i+400),gsl_vector_get(PO_V,i+400));
            }

            /*FLUID_PROP(Gi, Gj, Gk, Dk, PRBB,
                       SO, SG, SW, PO, PG, PW, PB,
                       BO, BG, BW, RS, VISO, VISG,
                       RHO, RHG, RHW, BPOR);*/

            FLUID_PROP(Gi, Gj, Gk, Dk, PRBB,
                       SO_V, SG_V, SW_V, PO_V, PG_V, PW_V, PB_V,
                       BO_V, BG_V, BW_V, RS_V, VISO_V, VISG_V,
                       RHO_V, RHG_V, RHW_V, BPOR_V);

            /*REL_PERM_DIST(Gi, Gj, Gk,
                          SO, SG, SW,
                          KRO, KRG, KRW);*/

            REL_PERM_DIST(Gi, Gj, Gk,
                          SO_V, SG_V, SW_V,
                          KRO_V, KRG_V, KRW_V);

            TRANS_1PW_PRESS(Gi, Gj, Gk, PO_V, PG_V, PW_V, BO_V, BG_V, BW_V, RS_V, VISO_V, VISG_V, VISW,
                        TFP_OI_V, TFP_OJ_V, TFP_OK_V,
                        TFP_GI_V, TFP_GJ_V, TFP_GK_V,
                        TFP_WI_V, TFP_WJ_V, TFP_WK_V,
                        TFP_ORSI_V, TFP_ORSJ_V, TFP_ORSK_V);

            /*TRANS_1PW_PRESS(Gi, Gj, Gk, PO, PG, PW, BO, BG, BW, RS, VISO, VISG, VISW,
                        TFP_OI, TFP_OJ, TFP_OK,
                        TFP_GI, TFP_GJ, TFP_GK,
                        TFP_WI, TFP_WJ, TFP_WK,
                        TFP_ORSI, TFP_ORSJ, TFP_ORSK);*/

            TRANS_1PW_SAT(Gi, Gj, Gk, KRO_V, KRG_V, KRW_V,
                      TFS_OI_V, TFS_OJ_V, TFS_OK_V,
                      TFS_GI_V, TFS_GJ_V, TFS_GK_V,
                      TFS_WI_V, TFS_WJ_V, TFS_WK_V);

            /*TRANS_1PW_SAT(Gi, Gj, Gk, KRO, KRG, KRW,
                      TFS_OI, TFS_OJ, TFS_OK,
                      TFS_GI, TFS_GJ, TFS_GK,
                      TFS_WI, TFS_WJ, TFS_WK);*/

            SLOPE_CALC(Gi, Gj, Gk, ZB,
                       PO, PO_V, PG, PG_V, PW, PW_V,
                       BPOR, BPOR_V, SLPOR,
                       BO, BO_V, SLBO,
                       BG, BG_V, SLBG,
                       BW, BW_V, SLBW,
                       RS, RS_V, SLRS);

            TRANSPORT_1PW_2(Gi, Gj, Gk, D_NV, ZB, SO_V, SW_V, SG_V, PO_V, PW_V, PG_V, PCGO_V, PCOW_V,
                            RHO, RHG, RHW, RHO_V, RHG_V, RHW_V, RS, RS_V,
                            TGI, TGJ,  TGK,
                            TFS_OI, TFS_OJ,  TFS_OK, TFS_GI,  TFS_GJ,  TFS_GK, TFS_WI,  TFS_WJ,  TFS_WK,
                            TFP_OI, TFP_OJ,  TFP_OK, TFP_GI,  TFP_GJ,  TFP_GK, TFP_WI,  TFP_WJ,  TFP_WK,
                            TFP_ORSI, TFP_ORSJ,  TFP_ORSK,
                            TFS_OI_V, TFS_OJ_V,  TFS_OK_V, TFS_GI_V,  TFS_GJ_V,  TFS_GK_V, TFS_WI_V,  TFS_WJ_V,  TFS_WK_V,
                            TFP_OI_V, TFP_OJ_V,  TFP_OK_V, TFP_GI_V,  TFP_GJ_V,  TFP_GK_V, TFP_WI_V,  TFP_WJ_V,  TFP_WK_V,
                            TFP_ORSI_V, TFP_ORSJ_V,  TFP_ORSK_V,
                            FNV, JACOBIAN);

            ACCUMULATION_RESID(1, Gi, Gj, Gk, Is, Js, Ks,
                               SO_V, SG_V, SW_V, BPOR_V, BO_V, BG_V, BW_V, RS_V, CV_WGO);

            ACCUMULATION_DIF(1, Gi, Gj, Gk, Is, Js, Ks,
                              SO, SO_V, SG, SG_V, SW, SW_V, BPOR, BPOR_V, SLPOR, BO, BO_V, SLBO, BG, BG_V, SLBG,
                              BW, BW_V, SLBW, RS, RS_V, SLRS, CM_WGO);

            WELL_PROD(Gi, Gj, Gk, 2, 2, 2, 1.1300643, 5250, 10000, D_NV,
                       PO, SW, SG, KRO, KRG, KRW, BO, BG, BW, RS, VISO, VISG,
                       PO_V, SW_V, SG_V, KRO_V, KRG_V, KRW_V, BO_V, BG_V, BW_V, RS_V, VISO_V, VISG_V,
                       VISW, QV, QRN);

            SOLVER(Gi, Gj, Gk, JACOBIAN, CM_WGO, QRN, FNV, CN_WGO, CV_WGO, QV, X_V, residual);

            PROP_2_VEC_WGP(Gi, Gj, Gk, SW_V, SG_V, PO_V, X_N);
            sss = rck - gsl_vector_get(residual,0);

            ++dn;
        } while ((sss > 0)&&(dn < 1));

        for (size_t i=3; i<4; i++) {
            printf("\nStats - Time: %d\n",Tsz-1);
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i),gsl_vector_get(SG,i),gsl_vector_get(SO,i),gsl_vector_get(PO,i));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+100),gsl_vector_get(SG,i+100),gsl_vector_get(SO,i+100),gsl_vector_get(PO,i+100));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+200),gsl_vector_get(SG,i+200),gsl_vector_get(SO,i+200),gsl_vector_get(PO,i+200));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+300),gsl_vector_get(SG,i+300),gsl_vector_get(SO,i+300),gsl_vector_get(PO,i+300));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n\n",gsl_vector_get(SW,i+400),gsl_vector_get(SG,i+400),gsl_vector_get(SO,i+400),gsl_vector_get(PO,i+400));
        }

        NEW_PROP_NP(X_NM, X_N, X_VP); /** X_VP -> X_N */
        gsl_matrix_set_col(RESULT,Tsz,X_N);
        UPD_SAT_PS(X_N, PRBB, SO, SG, SW, PO, PG, PW, PB);
        for (size_t i=3; i<4; i++) {
            printf("\nStats - Time: %d\n",Tsz);
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i),gsl_vector_get(SG,i),gsl_vector_get(SO,i),gsl_vector_get(PO,i));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+100),gsl_vector_get(SG,i+100),gsl_vector_get(SO,i+100),gsl_vector_get(PO,i+100));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+200),gsl_vector_get(SG,i+200),gsl_vector_get(SO,i+200),gsl_vector_get(PO,i+200));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n",gsl_vector_get(SW,i+300),gsl_vector_get(SG,i+300),gsl_vector_get(SO,i+300),gsl_vector_get(PO,i+300));
            printf("SW: %f | SG: %f | SO: %f | PO: %f\n\n",gsl_vector_get(SW,i+400),gsl_vector_get(SG,i+400),gsl_vector_get(SO,i+400),gsl_vector_get(PO,i+400));
        }
        printf("\nTermino Time: %d\n",Tsz);
        ++Tsz;

    } while (Tsz < 31);

    /**printf("Ini Per\n");
    for (size_t i=1; i<501; i++) {
        printf("B: %d  | SW: %f | SG: %f | SO: %f | PO: %f\n",i,gsl_vector_get(SW,i),gsl_vector_get(SG,i),gsl_vector_get(SO,i),gsl_vector_get(PO,i));
    }*/

    std::cout << 'A' << std::setw(6) << 'B' << std::endl;
    std::cout << std::setprecision(13);
    std::cout << "IMPSAT - Simulacao Terminada" << std::endl;

}
