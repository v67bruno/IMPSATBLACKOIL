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

void WELL_PROD(int Gi, int Gj, int Gk,
               int CI, int CJ, int CK, double WI, double BHP, double RATE,
               gsl_vector * D_NV,
               gsl_vector * PO_N, gsl_vector * SW_N, gsl_vector * SG_N,
               gsl_vector * KRO_N, gsl_vector * KRG_N, gsl_vector * KRW_N,
               gsl_vector * BO_N, gsl_vector * BG_N, gsl_vector * BW_N,
               gsl_vector * RS_N, gsl_vector * VISO_N, gsl_vector * VISG_N,
               gsl_vector * PO_V, gsl_vector * SW_V, gsl_vector * SG_V,
               gsl_vector * KRO_V, gsl_vector * KRG_V, gsl_vector * KRW_V,
               gsl_vector * BO_V, gsl_vector * BG_V, gsl_vector * BW_V,
               gsl_vector * RS_V, gsl_vector * VISO_V, gsl_vector * VISG_V,
               double VISW, gsl_vector * QV, gsl_matrix * QRN)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;

    /** for transporte */
    int tbn = 1;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {

                double QWSW, QWSG, QWPO,
                       QGSW, QGSG, QGPO,
                       QOSW, QOSG, QOPO,
                       P_PROD, P_PRODN,
                       UBW, UBG, UBO,
                       UBWN, UBGN, UBON,
                       QWV, QGV, QOV,
                       DFQW, DFQG, DFQO, DFQRS;
                double DIF_G_SW = 0, DIF_G_SG = 0, DIF_G_PO = 0;

                if ((CI == i) && (CJ == j) && (CK == k)) {

                    double UWI, UGI, UOI, UBWI, UBWIN, BGI, BOI, RSI, KRWI, KRGI, KROI;

                    DIF_G_SW = abs(gsl_vector_get(D_NV, tbn*3));
                    DIF_G_SG = abs(gsl_vector_get(D_NV, ((tbn*3)+1)));
                    DIF_G_PO = abs(gsl_vector_get(D_NV, ((tbn*3)+2)));

                    P_PROD = gsl_vector_get(PO_V, tbn) - BHP;
                    P_PRODN = gsl_vector_get(PO_N, tbn) - BHP;

                    UBW = gsl_vector_get(KRW_V, tbn)/(VISW * gsl_vector_get(BW_V, tbn));
                    UBWN = gsl_vector_get(KRW_N, tbn)/(VISW*gsl_vector_get(BW_N, tbn));
                    UBO = gsl_vector_get(KRO_V, tbn)/(gsl_vector_get(VISO_V, tbn)*gsl_vector_get(BO_V, tbn));
                    UBON = gsl_vector_get(KRO_N, tbn)/(gsl_vector_get(VISO_N, tbn)*gsl_vector_get(BO_N, tbn));
                    UBG = gsl_vector_get(KRG_V, tbn)/(gsl_vector_get(VISG_V, tbn)*gsl_vector_get(BG_V, tbn));
                    UBGN = gsl_vector_get(KRG_N, tbn)/(gsl_vector_get(VISG_N, tbn)*gsl_vector_get(BG_N, tbn));

                    printf("\nWell Prod Cal\nWb P: %f - dBHP: %f\n",gsl_vector_get(PO_N, tbn),P_PRODN);
                    printf("Wb UBO: %f\n",UBO);

                    DFQW = ((WI * UBW) - (WI * UBWN));
                    DFQO = ((WI * UBO) - (WI * UBON));
                    DFQG = ((WI * UBG) - (WI * UBGN));
                    DFQRS = ((WI * UBO * gsl_vector_get(RS_V, tbn)) - (WI * UBON * gsl_vector_get(RS_N, tbn)));


                    QWSW = abs(DFQW)/DIF_G_SW; QWSW *= P_PROD; QWSW *= -1;
                    QWSG = 0;
                    QWPO = abs(DFQW)/DIF_G_PO; QWPO *= P_PROD; QWPO *= -1;

                    QOSW = abs(DFQO)/DIF_G_SW; QOSW *= P_PROD; QOSW *= -1;
                    QOSG = abs(DFQO)/DIF_G_SG; QOSG *= P_PROD; QOSG *= -1;
                    QOPO = abs(DFQO)/DIF_G_PO; QOPO *= P_PROD; QOPO *= -1;

                    QGSW = abs(DFQRS)/DIF_G_SW; QGSW *= P_PROD; QGSW *= -1;
                    QGSG = ((abs(DFQG)/DIF_G_SG) * P_PROD);
                    QGSG += ((abs(DFQRS)/DIF_G_SG) * P_PROD);
                    QGSG *= -1;

                    QGPO = ((abs(DFQG)/DIF_G_PO) * P_PROD);
                    QGPO += ((abs(DFQRS)/DIF_G_PO) * P_PROD);
                    QGPO *= -1;

                    if (gsl_finite(QWSW) !=1){QWSW = 0;}
                    if (gsl_finite(QWPO) !=1){QWPO = 0;}
                    if (gsl_finite(QOSW) !=1){QOSW = 0;}
                    if (gsl_finite(QOSG) !=1){QOSG = 0;}
                    if (gsl_finite(QOPO) !=1){QOPO = 0;}
                    if (gsl_finite(QGSW) !=1){QGSW = 0;}
                    if (gsl_finite(QGSG) !=1){QGSG = 0;}
                    if (gsl_finite(QGPO) !=1){QGPO = 0;}

                    RESID_2_SPM_TRIP(tbn, tbn, QWSW, QWSG, QWPO, QGSW, QGSG, QGPO, QOSW, QOSG, QOPO, QRN);

                    printf("\nProd. Poco Bloco: %d\n",tbn);
                    QOV = WI * UBON * P_PRODN; QOV *= -1;
                    printf("Prod. O: %f\n",QOV);
                    QGV = (WI * UBGN * P_PRODN);
                    printf("Prod. G: %f\n",QGV);
                    QGV += (WI * UBON * P_PRODN * gsl_vector_get(RS_N, tbn));
                    QGV *= -1;
                    printf("Prod. G+RS: %f\n",QGV);
                    QWV = WI * UBWN * P_PRODN; QWV *= -1;
                    printf("Prod. W: %f\n",QWV);

                    printf("Derivative Q: %d\n",tbn);
                    printf("O: %f | %f | %f\n",QOSW, QOSG, QOPO);
                    printf("G: %f | %f | %f\n",QGSW, QGSG, QGPO);
                    printf("W: %f | %f | %f\n\n",QWSW, QWSG, QWPO);



                    RESID_2_VEC_WGO(tbn, QWV, QGV, QOV, QV);

                }

                tbn++;
            }
        }
    }/** end for transporte */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Well Prod ALL Rate Calc ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}/**Full Implicit */


void WELL_INJ_W(int Gi, int Gj, int Gk,
               int CI, int CJ, int CK, double WI, double BHP, double RATE,
               gsl_vector * PO_N, gsl_vector * PO_V,
               gsl_vector * SW_N, gsl_vector * SW_V,
               gsl_vector * SG_N, gsl_vector * SG_V,
               gsl_vector * KRO_N, gsl_vector * KRO_V,
               gsl_vector * KRG_N, gsl_vector * KRG_V,
               gsl_vector * KRW_N, gsl_vector * KRW_V,
               gsl_vector * BO_N, gsl_vector * BO_V,
               gsl_vector * BG_N, gsl_vector * BG_V,
               gsl_vector * BW_N, gsl_vector * BW_V,
               gsl_vector * RS_N, gsl_vector * RS_V,
               gsl_vector * VISO_N, gsl_vector * VISO_V,
               gsl_vector * VISG_N, gsl_vector * VISG_V,
               double VISW,
               gsl_vector * QV,
               gsl_matrix * QRN)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;

    /** for transporte */
    int tbn = 1;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {
                double QWSW, QWSG, QWPO,
                       QGSW, QGSG, QGPO,
                       QOSW, QOSG, QOPO,
                       PPROD_V, PPROD_N, DFPO, DFSW, DFSG,
                       UBW_V, UBW_N, DFWSW, DFWSG, DFWPO,
                       UBG_V, UBG_N, DFGSW, DFGSG, DFGPO,
                       UBO_V, UBO_N, DFOSW, DFOSG, DFOPO,
                       QW_V, QW_N, QOV, DFQW, DFQG, DFQO, DFQRS, WIB_V, WIB_N;

                if ((CI == i) && (CJ == j) && (CK == k)) {

                    PPROD_N = gsl_vector_get(PO_N, tbn) - BHP;
                    UBW_N = gsl_vector_get(KRW_N, tbn)/VISW;
                    UBG_N = gsl_vector_get(KRG_N, tbn)/gsl_vector_get(VISG_N, tbn);
                    UBO_N = gsl_vector_get(KRO_N, tbn)/gsl_vector_get(VISO_N, tbn);
                    WIB_N = WI / gsl_vector_get(BW_N, tbn);

                    PPROD_V = gsl_vector_get(PO_V, tbn) - BHP;
                    UBW_V = gsl_vector_get(KRW_V, tbn)/VISW;
                    UBG_V = gsl_vector_get(KRG_V, tbn)/gsl_vector_get(VISG_V, tbn);
                    UBO_V = gsl_vector_get(KRO_V, tbn)/gsl_vector_get(VISO_V, tbn);
                    WIB_V = WI / gsl_vector_get(BW_V, tbn);

                    QW_N = - WIB_N * (UBW_N+UBG_N+UBO_N) * PPROD_N;
                    QW_V = - WIB_V * (UBW_V+UBG_V+UBO_V) * PPROD_V;

                    QWSW = QW_N - QW_V;
                    QWSW /= (gsl_vector_get(SW_V, tbn) - gsl_vector_get(SW_N, tbn));

                    QWPO = QW_N - QW_V;
                    QWPO /= (gsl_vector_get(PO_V, tbn) - gsl_vector_get(PO_N, tbn));

                    RESID_2_VEC_WGO(tbn, QW_V, 0, 0, QV);

                    RESID_2_SPM_TRIP(tbn, tbn, QWSW, 0, QWPO, 0, 0, 0, 0, 0, 0, QRN);

                }
                tbn++;
            }
        }
    }/** end for transporte */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Well Inj W Rate Calc ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}/**Full Implicit */
