#include <iostream>
#include <iomanip>
#include <cmath>
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

void GEO_REC_CTER_BLCK(int Gi, int Gj, int Gk, int Li, int Lj, int Lk,
                         gsl_vector * KI, gsl_vector * KJ, gsl_vector * KK,
                         gsl_vector * GTI, gsl_vector * GTJ, gsl_vector * GTK)
{
    TimeVar t1=timeNow();
    for (int bn=1; bn<(GMIJK+1); bn++){
        /** A, K e Delta X = cte e retangular */
        double Isz = Li/Gi;
        double Jsz = Lj/Gj;
        double Ksz = Lk/Gk;

        double xAI = Jsz * Ksz;
        double xAJ = Isz * Ksz;
        double xAK = Isz * Jsz;

        double xKI = gsl_vector_get(KI, bn);
        double xKJ = gsl_vector_get(KJ, bn);
        double xKK = gsl_vector_get(KK, bn);

        double xMAKI = xAI * xKI;
        double xMAKJ = xAJ * xKJ;
        double xMAKK = xAK * xKK;

        double xDAKI = xMAKI / Isz;
        double xDAKJ = xMAKJ / Jsz;
        double xDAKK = xMAKK / Ksz;

        gsl_vector_set (GTI, bn, xDAKI);
        gsl_vector_set (GTJ, bn, xDAKJ);
        gsl_vector_set (GTK, bn, xDAKK);

    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Geometria dos blocos Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}  /** Geo Blockss end */

void TRANS_1PW_PRESS(int Gi, int Gj, int Gk,
                gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                gsl_vector * BO, gsl_vector * BG, gsl_vector * BW,
                gsl_vector * RS,
                gsl_vector * VISO, gsl_vector * VISG, double VISW,
                gsl_vector * TFP_OI, gsl_vector * TFP_OJ, gsl_vector * TFP_OK,
                gsl_vector * TFP_GI, gsl_vector * TFP_GJ, gsl_vector * TFP_GK,
                gsl_vector * TFP_WI, gsl_vector * TFP_WJ, gsl_vector * TFP_WK,
                gsl_vector * TFP_ORSI, gsl_vector * TFP_ORSJ, gsl_vector * TFP_ORSK)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    int bn = 1;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {
                double GF, VBO_M, VBG_M, VBW_M, RS_M, KROX, KRGX, KRWX;

                VBO_M = (1/gsl_vector_get(VISO, bn)/gsl_vector_get(BO, bn));
                VBG_M = (1/gsl_vector_get(VISG, bn)/gsl_vector_get(BG, bn));
                VBW_M = (1/VISW/gsl_vector_get(BW, bn));
                RS_M = gsl_vector_get(RS, bn);

                std::list<int> IJK {0,1,2};
                for (int FC : IJK) {
                    double TBO = 0, TBG = 0, TBW = 0, TBORS = 0;

                    TBO = VBO_M;
                    TBG = VBG_M;
                    TBW = VBW_M;
                    TBORS = VBO_M * RS_M;

                    if(FC == 0){
                        gsl_vector_set(TFP_OI, bn, TBO);
                        gsl_vector_set(TFP_GI, bn, TBG);
                        gsl_vector_set(TFP_WI, bn, TBW);
                        gsl_vector_set(TFP_ORSI, bn, TBORS);
                    }
                    if(FC == 1){
                        gsl_vector_set(TFP_OJ, bn, TBO);
                        gsl_vector_set(TFP_GJ, bn, TBG);
                        gsl_vector_set(TFP_WJ, bn, TBW);
                        gsl_vector_set(TFP_ORSJ, bn, TBORS);
                    }
                    if(FC == 2){
                        gsl_vector_set(TFP_OK, bn, TBO);
                        gsl_vector_set(TFP_GK, bn, TBG);
                        gsl_vector_set(TFP_WK, bn, TBW);
                        gsl_vector_set(TFP_ORSK, bn, TBORS);
                    }
                    FC++;
                }
                bn++;
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Trans-1PW-Press Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}

void TRANS_1PW_SAT(int Gi, int Gj, int Gk,
                gsl_vector * KRO, gsl_vector * KRG, gsl_vector * KRW,
                gsl_vector * TFS_OI, gsl_vector * TFS_OJ, gsl_vector * TFS_OK,
                gsl_vector * TFS_GI, gsl_vector * TFS_GJ, gsl_vector * TFS_GK,
                gsl_vector * TFS_WI, gsl_vector * TFS_WJ, gsl_vector * TFS_WK)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    int bn = 1;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {
                double GF, VBO_M, VBG_M, VBW_M, RS_M, KROX, KRGX, KRWX;

                KROX = gsl_vector_get(KRO, bn);
                KRGX = gsl_vector_get(KRG, bn);
                KRWX = gsl_vector_get(KRW, bn);

                std::list<int> IJK {0,1,2};
                for (int FC : IJK) {
                    double TBO = 0, TBG = 0, TBW = 0, TBORS = 0;

                    TBO = KROX;
                    TBG = KRGX;
                    TBW = KRWX;

                    if(FC == 0){
                        gsl_vector_set(TFS_OI, bn, TBO);
                        gsl_vector_set(TFS_GI, bn, TBG);
                        gsl_vector_set(TFS_WI, bn, TBW);
                    }
                    if(FC == 1){
                        gsl_vector_set(TFS_OJ, bn, TBO);
                        gsl_vector_set(TFS_GJ, bn, TBG);
                        gsl_vector_set(TFS_WJ, bn, TBW);
                    }
                    if(FC == 2){
                        gsl_vector_set(TFS_OK, bn, TBO);
                        gsl_vector_set(TFS_GK, bn, TBG);
                        gsl_vector_set(TFS_WK, bn, TBW);
                    }
                    FC++;
                }
                bn++;
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Trans-1PW-Sat Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}
