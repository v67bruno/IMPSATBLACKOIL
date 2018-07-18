#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <list>
#include <chrono>
#include <utility>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include "model.hpp"
#include "functions.hpp"

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

void TRANSPORT_1PW_2(int Gi, int Gj, int Gk, gsl_vector * D_NV,
                   gsl_vector * ZB,
                   gsl_vector * SO_V, gsl_vector * SW_V, gsl_vector * SG_V,
                   gsl_vector * PO_V, gsl_vector * PW_V, gsl_vector * PG_V, gsl_vector * PCGO_V, gsl_vector * PCOW_V,
                   gsl_vector * RHO_N, gsl_vector * RHG_N, gsl_vector * RHW_N,
                   gsl_vector * RHO_V, gsl_vector * RHG_V, gsl_vector * RHW_V,
                   gsl_vector * RS_N, gsl_vector * RS_V,
                   gsl_vector * TGI, gsl_vector * TGJ, gsl_vector * TGK,
                gsl_vector * TFS_OI_N, gsl_vector * TFS_OJ_N, gsl_vector * TFS_OK_N,
                gsl_vector * TFS_GI_N, gsl_vector * TFS_GJ_N, gsl_vector * TFS_GK_N,
                gsl_vector * TFS_WI_N, gsl_vector * TFS_WJ_N, gsl_vector * TFS_WK_N,
                gsl_vector * TFP_OI_N, gsl_vector * TFP_OJ_N, gsl_vector * TFP_OK_N,
                gsl_vector * TFP_GI_N, gsl_vector * TFP_GJ_N, gsl_vector * TFP_GK_N,
                gsl_vector * TFP_WI_N, gsl_vector * TFP_WJ_N, gsl_vector * TFP_WK_N,
                gsl_vector * TFP_ORSI_N, gsl_vector * TFP_ORSJ_N, gsl_vector * TFP_ORSK_N,
                gsl_vector * TFS_OI_V, gsl_vector * TFS_OJ_V, gsl_vector * TFS_OK_V,
                gsl_vector * TFS_GI_V, gsl_vector * TFS_GJ_V, gsl_vector * TFS_GK_V,
                gsl_vector * TFS_WI_V, gsl_vector * TFS_WJ_V, gsl_vector * TFS_WK_V,
                gsl_vector * TFP_OI_V, gsl_vector * TFP_OJ_V, gsl_vector * TFP_OK_V,
                gsl_vector * TFP_GI_V, gsl_vector * TFP_GJ_V, gsl_vector * TFP_GK_V,
                gsl_vector * TFP_WI_V, gsl_vector * TFP_WJ_V, gsl_vector * TFP_WK_V,
                gsl_vector * TFP_ORSI_V, gsl_vector * TFP_ORSJ_V, gsl_vector * TFP_ORSK_V,
                gsl_vector * FNV, gsl_matrix * JACOBIAN)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;

    int N = 1; int tbnIP, tbnJP, tbnKP, tbnIM, tbnJM, tbnKM;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {

                if (i == (Gi-1)){ tbnIP = Gijk+1; } else { tbnIP = N + 1; }
                if (i == 0) { tbnIM = 0; } else { tbnIM = N - 1; }
                if (j == (Gj-1)){ tbnJP =  Gijk+1; } else { tbnJP = N + Gi; }
                if (j == 0) { tbnJM = 0; } else { tbnJM = N - Gi; }
                if (k == (Gk-1)){ tbnKP = Gijk+1; } else { tbnKP = N + Gij; }
                if (k == 0) { tbnKM = 0; } else { tbnKM = N - Gij; }
                //printf("\nBloco IJK: %d\n",N);
                double FWN = 0, FGN = 0, FON = 0;
                double RWSW_N = 0, RWPO_N = 0,
                       RGSW_N = 0, RGSG_N = 0, RGPO_N = 0,
                       ROSW_N = 0, ROSG_N = 0, ROPO_N = 0,
                       RRSSW_N = 0, RRSSG_N = 0, RRSPO_N = 0;

                std::list<int> dfb = { tbnIP, tbnJP, tbnKP, tbnIM, tbnJM, tbnKM };
                int MC = 0;
                for (int M : dfb) {
                    if ((M != 0) && (M != Gijk+1)) {
                        double TO_UPS_V = 0, TG_UPS_V = 0, TW_UPS_V = 0, TORS_UPS_V = 0,
                        TO_UPS_N = 0, TG_UPS_N = 0, TW_UPS_N = 0, TORS_UPS_N = 0,
                        PCGO_UPS_V = 0, PCGO_UPS_N = 0, PCOW_UPS_V = 0, PCOW_UPS_N = 0,
                        DFTG = 0, DFTW = 0, DFTO = 0, DFTRS = 0,
                        DFTO_SW = 0, DFTO_SG = 0, DFTO_PO = 0,
                        DFTO_SW_N = 0, DFTO_SG_N = 0, DFTO_PO_N = 0,
                        DFTG_SW = 0, DFTG_SG = 0, DFTG_PO = 0,
                        DFTG_SW_N = 0, DFTG_SG_N = 0, DFTG_PO_N = 0,
                        DFTW_SW = 0, DFTW_SG = 0, DFTW_PO = 0,
                        DFTW_SW_N = 0, DFTW_SG_N = 0, DFTW_PO_N = 0,
                        DFTORS_SW = 0, DFTORS_SG = 0, DFTORS_PO = 0,
                        DFTORS_SW_N = 0, DFTORS_SG_N = 0, DFTORS_PO_N = 0,
                        DFPCGO_N = 0, DFPCOW_N = 0,
                        TGFW = 0, TGFO = 0, TGFG = 0,
                        W_TFS = 0, W_TFP = 0, O_TFS = 0, O_TFP = 0, G_TFS = 0, G_TFP = 0, RS_TFP = 0,
                        W_TFS_DF_SW = 0, W_TFP_DF_PO = 0, G_TFS_DF_SG = 0, G_TFP_DF_PO = 0,
                        O_TFS_DF_SW = 0, O_TFS_DF_SG = 0, O_TFP_DF_PO = 0,
                        RS_TFS_DF_SW = 0, RS_TFS_DF_SG = 0, RS_TFP_DF_PO = 0;

                        int UPS_O, DWS_O, UPS_W, DWS_W, UPS_G, DWS_G;

                        double DFPO = gsl_vector_get(PO_V, M) - gsl_vector_get(PO_V, N);
                        double DFPCGO = gsl_vector_get(PCGO_V, M) - gsl_vector_get(PCGO_V, N);
                        double DFPCOW = gsl_vector_get(PCOW_V, M) - gsl_vector_get(PCOW_V, N);
                        double GMO_M = 0.5 * 9.8066352 * (gsl_vector_get(RHO_N, M) + gsl_vector_get(RHO_N, N));
                        double GMG_M = 0.5 * 9.8066352 * (gsl_vector_get(RHG_N, M) + gsl_vector_get(RHG_N, N));
                        double GMW_M = 0.5 * 9.8066352 * (gsl_vector_get(RHW_N, M) + gsl_vector_get(RHW_N, N));
                        double DFZB = gsl_vector_get(ZB, M) - gsl_vector_get(ZB, N);

                        double RPW = DFPO - DFPCOW - (GMW_M * DFZB);
                        double RPG = DFPO + DFPCGO - (GMG_M * DFZB);
                        double RPO = DFPO - (GMO_M * DFZB);

                        std::tie(UPS_O, DWS_O) = FASE_POT(N, M, PO_V, RHO_V, ZB);
                        std::tie(UPS_W, DWS_W) = FASE_POT(N, M, PW_V, RHW_V, ZB);
                        std::tie(UPS_G, DWS_G) = FASE_POT(N, M, PG_V, RHG_V, ZB);

                        if ((MC == 0) || (MC == 3)) {
                            TW_UPS_V = UPS_TRANS_NM(UPS_W, DWS_W, TGI, TFP_WI_V, TFS_WI_V);
                            TG_UPS_V = UPS_TRANS_NM(UPS_G, DWS_G, TGI, TFP_GI_V, TFS_GI_V);
                            TO_UPS_V = UPS_TRANS_NM(UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V);
                            TORS_UPS_V = UPS_TRANS_NM(UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V);
                        }
                        if ((MC == 1) || (MC == 4)) {
                            TW_UPS_V = UPS_TRANS_NM(UPS_W, DWS_W, TGJ, TFP_WJ_V, TFS_WJ_V);
                            TG_UPS_V = UPS_TRANS_NM(UPS_G, DWS_G, TGJ, TFP_GJ_V, TFS_GJ_V);
                            TO_UPS_V = UPS_TRANS_NM(UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V);
                            TORS_UPS_V = UPS_TRANS_NM(UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V);
                        }
                        if ((MC == 2) || (MC == 5)) {
                            TW_UPS_V = UPS_TRANS_NM(UPS_W, DWS_W, TGK, TFP_WK_V, TFS_WK_V);
                            TG_UPS_V = UPS_TRANS_NM(UPS_G, DWS_G, TGK, TFP_GK_V, TFS_GK_V);
                            TO_UPS_V = UPS_TRANS_NM(UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V);
                            TORS_UPS_V = UPS_TRANS_NM(UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V);
                        }

                        FWN += (TW_UPS_V * RPW);
                        FGN += ((TG_UPS_V * RPG) + (TORS_UPS_V * RPO));
                        FON += (TO_UPS_V * RPO);

                        double DIF_G_SW = 0, DIF_G_SG = 0, DIF_G_PO = 0;

                        std::list<int> flow = { UPS_O, UPS_G, UPS_W };
                        int flc = 0;
                        for (int dirc : flow) {
                            DIF_G_SW = abs(gsl_vector_get(D_NV, dirc*3));
                            DIF_G_SG = abs(gsl_vector_get(D_NV, ((dirc*3)+1)));
                            DIF_G_PO = abs(gsl_vector_get(D_NV, ((dirc*3)+2)));
                            /** N/M */
                            if (dirc != N) {
                                //printf("dirc != N - %d,%d\n",dirc,N);
                                if (flc == 0) {
                                    if ((MC == 0) || (MC == 3)) {
                                        DFTO_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V, TFP_OI_N, TFS_OI_N);
                                        DFTO_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V, TFP_OI_N, TFS_OI_N);
                                        DFTO_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V, TFP_OI_N, TFS_OI_N);
                                    }
                                    if ((MC == 1) || (MC == 4)) {
                                        DFTO_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V, TFP_OJ_N, TFS_OJ_N);
                                        DFTO_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V, TFP_OJ_N, TFS_OJ_N);
                                        DFTO_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V, TFP_OJ_N, TFS_OJ_N);
                                    }
                                    if ((MC == 2) || (MC == 5)) {
                                        DFTO_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V, TFP_OK_N, TFS_OK_N);
                                        DFTO_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V, TFP_OK_N, TFS_OK_N);
                                        DFTO_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V, TFP_OK_N, TFS_OK_N);
                                    }
                                }
                                if (flc == 1) {
                                    if ((MC == 0) || (MC == 3)) {
                                        DFTG_SG = DS_TRANS_NM(DIF_G_SG, UPS_G, DWS_G, TGI, TFP_GI_V, TFS_GI_V, TFP_GI_N, TFS_GI_N);
                                        DFTG_PO = DP_TRANS_NM(DIF_G_PO, UPS_G, DWS_G, TGI, TFP_GI_V, TFS_GI_V, TFP_GI_N, TFS_GI_N);
                                        DFTORS_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V, TFP_ORSI_N, TFS_OI_N);
                                        DFTORS_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V, TFP_ORSI_N, TFS_OI_N);
                                        DFTORS_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V, TFP_ORSI_N, TFS_OI_N);
                                    }
                                    if ((MC == 1) || (MC == 4)) {
                                        DFTG_SG = DS_TRANS_NM(DIF_G_SG, UPS_G, DWS_G, TGJ, TFP_GJ_V, TFS_GJ_V, TFP_GJ_N, TFS_GJ_N);
                                        DFTG_PO = DP_TRANS_NM(DIF_G_PO, UPS_G, DWS_G, TGJ, TFP_GJ_V, TFS_GJ_V, TFP_GJ_N, TFS_GJ_N);
                                        DFTORS_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V, TFP_ORSJ_N, TFS_OJ_N);
                                        DFTORS_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V, TFP_ORSJ_N, TFS_OJ_N);
                                        DFTORS_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V, TFP_ORSJ_N, TFS_OJ_N);
                                    }
                                    if ((MC == 2) || (MC == 5)) {
                                        DFTG_SG = DS_TRANS_NM(DIF_G_SG, UPS_G, DWS_G, TGK, TFP_GK_V, TFS_GK_V, TFP_GK_N, TFS_GK_N);
                                        DFTG_PO = DP_TRANS_NM(DIF_G_PO, UPS_G, DWS_G, TGK, TFP_GK_V, TFS_GK_V, TFP_GK_N, TFS_GK_N);
                                        DFTORS_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V, TFP_ORSK_N, TFS_OK_N);
                                        DFTORS_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V, TFP_ORSK_N, TFS_OK_N);
                                        DFTORS_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V, TFP_ORSK_N, TFS_OK_N);
                                    }
                                }
                                if (flc == 2) {
                                    if ((MC == 0) || (MC == 3)) {
                                        DFTW_SW = DS_TRANS_NM(DIF_G_SW, UPS_W, DWS_W, TGI, TFP_WI_V, TFS_WI_V, TFP_WI_N, TFS_WI_N);
                                        DFTW_PO = DP_TRANS_NM(DIF_G_PO, UPS_W, DWS_W, TGI, TFP_WI_V, TFS_WI_V, TFP_WI_N, TFS_WI_N);
                                    }
                                    if ((MC == 1) || (MC == 4)) {
                                        DFTW_SW = DS_TRANS_NM(DIF_G_SW, UPS_W, DWS_W, TGJ, TFP_WJ_V, TFS_WJ_V, TFP_WJ_N, TFS_WJ_N);
                                        DFTW_PO = DP_TRANS_NM(DIF_G_PO, UPS_W, DWS_W, TGJ, TFP_WJ_V, TFS_WJ_V, TFP_WJ_N, TFS_WJ_N);
                                    }
                                    if ((MC == 2) || (MC == 5)) {
                                        DFTW_SW = DS_TRANS_NM(DIF_G_SW, UPS_W, DWS_W, TGK, TFP_WK_V, TFS_WK_V, TFP_WK_N, TFS_WK_N);
                                        DFTW_PO = DP_TRANS_NM(DIF_G_PO, UPS_W, DWS_W, TGK, TFP_WK_V, TFS_WK_V, TFP_WK_N, TFS_WK_N);
                                    }
                                }
                                //printf("N/M | Fluid: [%d] | Flow M[%d] -> N[%d] = DFT/DFX\n",flc,nb,tbn);
                            }
                            else {  if (flc == 0) {DFTO_SW = 0; DFTO_SG = 0; DFTO_PO = 0;}
                                    if (flc == 1) {DFTG_SW = 0; DFTG_SG = 0; DFTG_PO = 0; DFPCGO = 0; DFTORS_PO = 0; DFTORS_SG = 0; DFTORS_SW = 0;}
                                    if (flc == 2) {DFTW_SW = 0; DFTW_SG = 0; DFTW_PO = 0; DFPCOW = 0;}
                                    //printf("N/M | Fluid: [%d] | Flow N[%d] -> M[%d] = 0\n",flc,tbn,nb);
                            }
                            flc++;
                        } /** end check DF N/M - N/N */

                        std::list<int> flown = { UPS_O, UPS_G, UPS_W };
                        int flcn = 0;
                        for (int dircn : flown) {
                            DIF_G_SW = abs(gsl_vector_get(D_NV, dircn*3));
                            DIF_G_SG = abs(gsl_vector_get(D_NV, ((dircn*3)+1)));
                            DIF_G_PO = abs(gsl_vector_get(D_NV, ((dircn*3)+2)));
                            /** N/N */
                            if (dircn == N) {
                                //printf("dircn == N - %d,%d\n",dircn,N);
                                if (flcn == 0) {
                                    if ((MC == 0) || (MC == 3)) {
                                        DFTO_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V, TFP_OI_N, TFS_OI_N);
                                        DFTO_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V, TFP_OI_N, TFS_OI_N);
                                        DFTO_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGI, TFP_OI_V, TFS_OI_V, TFP_OI_N, TFS_OI_N);
                                    }
                                    if ((MC == 1) || (MC == 4)) {
                                        DFTO_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V, TFP_OJ_N, TFS_OJ_N);
                                        DFTO_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V, TFP_OJ_N, TFS_OJ_N);
                                        DFTO_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGJ, TFP_OJ_V, TFS_OJ_V, TFP_OJ_N, TFS_OJ_N);
                                    }
                                    if ((MC == 2) || (MC == 5)) {
                                        DFTO_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V, TFP_OK_N, TFS_OK_N);
                                        DFTO_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V, TFP_OK_N, TFS_OK_N);
                                        DFTO_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGK, TFP_OK_V, TFS_OK_V, TFP_OK_N, TFS_OK_N);
                                    }
                                }
                                if (flcn == 1) {
                                    if ((MC == 0) || (MC == 3)) {
                                        DFTG_SG = DS_TRANS_NM(DIF_G_SG, UPS_G, DWS_G, TGI, TFP_GI_V, TFS_GI_V, TFP_GI_N, TFS_GI_N);
                                        DFTG_PO = DP_TRANS_NM(DIF_G_PO, UPS_G, DWS_G, TGI, TFP_GI_V, TFS_GI_V, TFP_GI_N, TFS_GI_N);
                                        DFTORS_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V, TFP_ORSI_N, TFS_OI_N);
                                        DFTORS_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V, TFP_ORSI_N, TFS_OI_N);
                                        DFTORS_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGI, TFP_ORSI_V, TFS_OI_V, TFP_ORSI_N, TFS_OI_N);
                                    }
                                    if ((MC == 1) || (MC == 4)) {
                                        DFTG_SG = DS_TRANS_NM(DIF_G_SG, UPS_G, DWS_G, TGJ, TFP_GJ_V, TFS_GJ_V, TFP_GJ_N, TFS_GJ_N);
                                        DFTG_PO = DP_TRANS_NM(DIF_G_PO, UPS_G, DWS_G, TGJ, TFP_GJ_V, TFS_GJ_V, TFP_GJ_N, TFS_GJ_N);
                                        DFTORS_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V, TFP_ORSJ_N, TFS_OJ_N);
                                        DFTORS_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V, TFP_ORSJ_N, TFS_OJ_N);
                                        DFTORS_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGJ, TFP_ORSJ_V, TFS_OJ_V, TFP_ORSJ_N, TFS_OJ_N);
                                    }
                                    if ((MC == 2) || (MC == 5)) {
                                        DFTG_SG = DS_TRANS_NM(DIF_G_SG, UPS_G, DWS_G, TGK, TFP_GK_V, TFS_GK_V, TFP_GK_N, TFS_GK_N);
                                        DFTG_PO = DP_TRANS_NM(DIF_G_PO, UPS_G, DWS_G, TGK, TFP_GK_V, TFS_GK_V, TFP_GK_N, TFS_GK_N);
                                        DFTORS_SW = DS_TRANS_NM(DIF_G_SW, UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V, TFP_ORSK_N, TFS_OK_N);
                                        DFTORS_SG = DS_TRANS_NM(DIF_G_SG, UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V, TFP_ORSK_N, TFS_OK_N);
                                        DFTORS_PO = DP_TRANS_NM(DIF_G_PO, UPS_O, DWS_O, TGK, TFP_ORSK_V, TFS_OK_V, TFP_ORSK_N, TFS_OK_N);
                                    }
                                }
                                if (flcn == 2) {
                                    if ((MC == 0) || (MC == 3)) {
                                        DFTW_SW = DS_TRANS_NM(DIF_G_SW, UPS_W, DWS_W, TGI, TFP_WI_V, TFS_WI_V, TFP_WI_N, TFS_WI_N);
                                        DFTW_PO = DP_TRANS_NM(DIF_G_PO, UPS_W, DWS_W, TGI, TFP_WI_V, TFS_WI_V, TFP_WI_N, TFS_WI_N);
                                    }
                                    if ((MC == 1) || (MC == 4)) {
                                        DFTW_SW = DS_TRANS_NM(DIF_G_SW, UPS_W, DWS_W, TGJ, TFP_WJ_V, TFS_WJ_V, TFP_WJ_N, TFS_WJ_N);
                                        DFTW_PO = DP_TRANS_NM(DIF_G_PO, UPS_W, DWS_W, TGJ, TFP_WJ_V, TFS_WJ_V, TFP_WJ_N, TFS_WJ_N);
                                    }
                                    if ((MC == 2) || (MC == 5)) {
                                        DFTW_SW = DS_TRANS_NM(DIF_G_SW, UPS_W, DWS_W, TGK, TFP_WK_V, TFS_WK_V, TFP_WK_N, TFS_WK_N);
                                        DFTW_PO = DP_TRANS_NM(DIF_G_PO, UPS_W, DWS_W, TGK, TFP_WK_V, TFS_WK_V, TFP_WK_N, TFS_WK_N);
                                    }
                                }
                                //printf("N/N | Fluid: [%d] | Flow N[%d] -> M[%d] = DFT/DFX\n",flc,tbn,nb);
                            }
                            else {  if (flcn == 0) {DFTO_SW_N = 0; DFTO_SG_N = 0; DFTO_PO_N = 0;}
                                    if (flcn == 1) {DFTG_SW_N = 0; DFTG_SG_N = 0; DFTG_PO_N = 0; DFPCGO_N = 0; DFTORS_PO_N = 0; DFTORS_SG_N = 0; DFTORS_SW_N = 0;}
                                    if (flcn == 2) {DFTW_SW_N = 0; DFTW_SG_N = 0; DFTW_PO_N = 0; DFPCOW_N = 0;}
                                    //printf("N/N | Fluid: [%d] | Flow M[%d] -> N[%d] = 0\n",flc,nb,tbn);
                            }
                            flcn++;
                        }

                        /**printf("\nN/M - Bloco M: %d\n",M);
                        printf("UPS_O: %d, UPS_G: %d, UPS_W: %d,\n",UPS_O, UPS_G, UPS_W);
                        printf("TO - SW: %f, SG: %f, PO: %f\n",DFTO_SW,DFTO_SG,DFTO_PO);
                        printf("TG - SW: %f, SG: %f, PO: %f\n",DFTG_SW,DFTG_SG,DFTG_PO);
                        printf("TW - SW: %f, SG: %f, PO: %f\n",DFTW_SW,DFTW_SG,DFTW_PO);
                        printf("TORS - SW: %f, SG: %f, PO: %f\n",DFTORS_SW,DFTORS_SG,DFTORS_PO);
                        printf("N/N - Bloco N: %d\n",N);
                        printf("UPS_O: %d, UPS_G: %d, UPS_W: %d,\n",UPS_O, UPS_G, UPS_W);
                        printf("TO - SW: %f, SG: %f, PO: %f\n",DFTO_SW_N,DFTO_SG_N,DFTO_PO_N);
                        printf("TG - SW: %f, SG: %f, PO: %f\n",DFTG_SW_N,DFTG_SG_N,DFTG_PO_N);
                        printf("TW - SW: %f, SG: %f, PO: %f\n",DFTW_SW_N,DFTW_SG_N,DFTW_PO_N);
                        printf("TORS - SW: %f, SG: %f, PO: %f\n",DFTORS_SW_N,DFTORS_SG_N,DFTORS_PO_N);*/

                        double RWSW_M = 0, RWPO_M = 0,
                               RGSW_M = 0, RGSG_M = 0, RGPO_M = 0,
                               ROSW_M = 0, ROSG_M = 0, ROPO_M = 0,
                               RRSSW_M = 0, RRSSG_M = 0, RRSPO_M = 0;

                        if (gsl_vector_get(SW_V, M) > 0){
                            RWSW_M = (RPW * DFTW_SW) - (TW_UPS_V * DFPCOW);
                            RWPO_M = TW_UPS_V + (RPW * DFTW_PO);
                        }
                        if (gsl_vector_get(SG_V, M) > 0){
                            RGSW_M = (RPO * DFTORS_SW);
                            RGSG_M = (RPG * DFTG_SG) + (TG_UPS_V * DFPCGO) + (RPO * DFTORS_SG);
                            RGPO_M = TG_UPS_V + (RPG * DFTG_PO) + TORS_UPS_V + (RPO * DFTORS_PO);
                        } else if ((gsl_vector_get(SG_V, M) < 0.0001) && (gsl_vector_get(SO_V, M) > 0)){
                            RGSW_M = (RPO * DFTORS_SW);
                            RGSG_M = (RPO * DFTORS_SG);
                            RGPO_M = TORS_UPS_V + (RPO * DFTORS_PO);
                        }
                        if (gsl_vector_get(SO_V, M) > 0){
                            ROSW_M = (RPO * DFTO_SW);
                            ROSG_M = (RPO * DFTO_SG);
                            ROPO_M = TO_UPS_V + (RPO * DFTO_PO);
                        }

                        /*printf("\nRES M: %d, N: %d\n",M,N);
                        printf("UPS_O: %d, UPS_G: %d, UPS_W: %d,\n",UPS_O, UPS_G, UPS_W);
                        printf("RWSW_M: %f, RWPO_M: %f\n",RWSW_M,RWPO_M);
                        printf("RGSW_M: %f, RGSG_M: %f, RGPO_M: %f\n",RGSW_M,RGSG_M,RGPO_M);
                        printf("ROSW_M: %f, ROSG_M: %f, ROPO_M: %f\n",ROSW_M,ROSG_M,ROPO_M);*/

                        if (gsl_finite(RWSW_M) !=1){RWSW_M = 0;}
                        if (gsl_finite(RWPO_M) !=1){RWPO_M = 0;}
                        if (gsl_finite(RGSW_M) !=1){RGSW_M = 0;}
                        if (gsl_finite(RGSG_M) !=1){RGSG_M = 0;}
                        if (gsl_finite(RGPO_M) !=1){RGPO_M = 0;}
                        if (gsl_finite(ROSW_M) !=1){ROSW_M = 0;}
                        if (gsl_finite(ROSG_M) !=1){ROSG_M = 0;}
                        if (gsl_finite(ROPO_M) !=1){ROPO_M = 0;}

                        RESID_2_SPM_TRIP(N, M, RWSW_M, 0, RWPO_M, RGSW_M, RGSG_M, RGPO_M, ROSW_M, ROSG_M, ROPO_M, JACOBIAN);


                        if (gsl_vector_get(SW_V, N) > 0){
                            RWSW_N += (RPW * DFTW_SW_N) - (TW_UPS_V * DFPCOW_N);
                            RWPO_N += - TW_UPS_V + (RPW * DFTW_PO_N);
                        }
                        if (gsl_vector_get(SG_V, N) > 0){
                            RGSW_N += (RPO * DFTORS_SW_N);
                            RGSG_N += (RPG * DFTG_SG_N) - (TG_UPS_V * DFPCGO_N) + (RPO * DFTORS_SG_N);
                            RGPO_N += - TG_UPS_V + (RPG * DFTG_PO_N) - TORS_UPS_V + (RPO * DFTORS_PO_N);
                        } else if ((gsl_vector_get(SG_V, N) < 0.01) && (gsl_vector_get(SO_V, N) > 0)){
                            RGSW_N += (RPO * DFTORS_SW_N);
                            RGSG_N += (RPO * DFTORS_SG_N);
                            RGPO_N += TORS_UPS_V + (RPO * DFTORS_PO_N);
                        }
                        if (gsl_vector_get(SO_V, N) > 0){
                            ROSW_N += (RPO * DFTO_SW_N);
                            ROSG_N += (RPO * DFTO_SG_N);
                            ROPO_N += - TO_UPS_V + (RPO * DFTO_PO_N);
                        }

                        /*printf("\nN/M - Bloco N: %d, Bloco M: %d\n",M,N);
                        printf("UPS_O: %d, UPS_G: %d, UPS_W: %d,\n",UPS_O, UPS_G, UPS_W);
                        printf("RWPO_M: %f, RWPO_M: %f\n",RWSW_N,RWPO_N);
                        printf("RGSW_M: %f, RGSG_M: %f, RGPO_M: %f\n",RGSW_N,RGSG_N,RGPO_M);
                        printf("ROSW_M: %f, ROSG_M: %f, ROPO_M: %f\n",ROSW_N,ROSG_N,ROPO_N);*/

                    }/** corner checker*/
                    MC++;
                }/** end tbnIP, tbnJP, tbnKP, tbnIM, tbnJM, tbnKM */
                /** {F} */
                if (gsl_finite(FWN) !=1){FWN = 0;}
                if (gsl_finite(FGN) !=1){FGN = 0;}
                if (gsl_finite(FON) !=1){FON = 0;}
                RESID_2_VEC_WGO(N, FWN, FGN, FON, FNV);
                /** [F] */
                if (gsl_finite(RWSW_N) !=1){RWSW_N = 0;}
                if (gsl_finite(RWPO_N) !=1){RWPO_N = 0;}
                if (gsl_finite(RGSW_N) !=1){RGSW_N = 0;}
                if (gsl_finite(RGSG_N) !=1){RGSG_N = 0;}
                if (gsl_finite(RGPO_N) !=1){RGPO_N = 0;}
                if (gsl_finite(ROSW_N) !=1){ROSW_N = 0;}
                if (gsl_finite(ROSG_N) !=1){ROSG_N = 0;}
                if (gsl_finite(ROPO_N) !=1){ROPO_N = 0;}

                /*printf("\nRES N: %d\n",N);
                /*printf("UPS_O: %d, UPS_G: %d, UPS_W: %d,\n",UPS_O, UPS_G, UPS_W);
                printf("RWSW_N: %f, RWPO_N: %f\n",RWSW_N,RWPO_N);
                printf("RGSW_N: %f, RGSG_N: %f, RGPO_N: %f\n",RGSW_N,RGSG_N,RGPO_N);
                printf("ROSW_N: %f, ROSG_N: %f, ROPO_N: %f\n",ROSW_N,ROSG_N,ROPO_N);*/

                RESID_2_SPM_TRIP(N, N, RWSW_N, 0, RWPO_N, RGSW_N, RGSG_N, RGPO_N, ROSW_N, ROSG_N, ROPO_N, JACOBIAN);

                /*printf("RW - SW: %f, PO: %f\n",RWSW_N,RWPO_N);
                printf("RG - SW: %f, SG: %f, PO: %f\n",RGSW_N, RGSG_N, RGPO_N);
                printf("RO - SW: %f, SG: %f, PO: %f\n",ROSW_N, ROSG_N, ROPO_N);*/

                N++;
            }
        }
    }/** end for transporte*/
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Transport Calc ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}/**Full Implicit*/
