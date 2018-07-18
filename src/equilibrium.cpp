#include <iostream>
#include <iomanip>
#include <cmath>
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

void EQUILIBRIUM (int Gi, int Gj, int Gk, int Ks,
                  int DWOC, int DGOC, int * Dk,
                  double REFDP, double REFPRES, double PRBB,
                  double * TSW, double * TSL, double * TPCGO, double * TPCOW,
                  gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                  gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                  gsl_vector * PB, gsl_vector * PCGO, gsl_vector * PCOW,
                  gsl_vector * ZB)
{

    TimeVar t1=timeNow();
    /** Inicia Modelo Equilibro Vertical*/
    /**szT -> Tempo total do modelo
    //Calcula das pressoes de contato
    //undersaturated state
    //contato gas/oleo
    //PG=PO -> PGO=PG-PO=0 -> PGO = PO(conhece a pressao)
    //Ph=Pi+PGO
    //contato oleo/agua
    //PO=PW -> POW=PO-PW=0 -> POW = PO(conhece a pressao)
    //Ph=Pi+POW*/

    double BOPB = INT_POL(TPVT, TBO, PRBB, sizeof(TPVT)/sizeof(double));
    double RSPB = INT_POL(TPVT, TRS, PRBB, sizeof(TPVT)/sizeof(double));
    double RhoOx = RHO_OOUS(RHOO, BOPB, CO, REFPRES, PRBB);
    double RhoGx = RHO_GOUS(RSPB, RHOG, BOPB, CO, REFPRES, PRBB);
    double PhGO = RF_PRES(REFPRES, REFDP, RhoOx+RhoGx, Dk[DGOC-1]);
    double PhOW = RF_PRES(REFPRES, REFDP, RhoOx+RhoGx, Dk[DWOC-1]);

    int bn = 1;
    for (int k=0; k<Gk; k++){/** Equilibirum Vertical */
        for (int j=0; j<Gj; j++){
            for (int i=0; i<Gi; i++){

                /**bloco centralizado*/
                double Zk = Ks * (k+1); Zk = Zk - (Ks/2);

                gsl_vector_set(ZB, bn, Zk);

                if (k == (DGOC-1)) {/** contato gas/oleo */
                     gsl_vector_set (SW, bn, 0.5);
                     gsl_vector_set (SG, bn, 0.1);
                     gsl_vector_set (SO, bn, 0.4);
                     gsl_vector_set (PO, bn, 10500);
                     gsl_vector_set (PW, bn, 10500);
                     gsl_vector_set (PG, bn, 10500);
                     //gsl_vector_set (PO, bn, PhGO);

                     //double xSL = gsl_vector_get(SW, bn) + gsl_vector_get(SO, bn);
                     //double xSW = gsl_vector_get(SW, bn);
                     //double xPCGO = INT_POL(TSL, TPCGO, xSL, sizeof(TSL)/sizeof(double));
                     //gsl_vector_set(PCGO, bn, 0);
                     //double xPCOW = INT_POL(TSW, TPCOW, xSW, sizeof(TSW)/sizeof(double));
                     //gsl_vector_set(PCOW, bn, 0);

                     //gsl_vector_set (PG, bn, (PhGO + gsl_vector_get(PCGO, bn)));
                     //gsl_vector_set (PW, bn, (PhGO - gsl_vector_get(PCOW, bn)));

                     if (gsl_vector_get(SG, bn) > 0) {
                        gsl_vector_set(PB, bn, gsl_vector_get(PG, bn));
                     } else {
                        gsl_vector_set(PB, bn, PRBB);
                     }

                } else if (k == (DWOC-1)) {/** contato oleo/agua */

                     gsl_vector_set (SW, bn, 0.5);
                     gsl_vector_set (SG, bn, 0.1);
                     gsl_vector_set (SO, bn, 0.4);
                     gsl_vector_set (PO, bn, 10500);
                     gsl_vector_set (PW, bn, 10500);
                     gsl_vector_set (PG, bn, 10500);
                     //gsl_vector_set (PO, bn, PhOW);
                     //gsl_vector_set (PG, bn, (PhOW + gsl_vector_get(PCGO, bn)));

                     //double xSL = gsl_vector_get(SW, bn) + gsl_vector_get(SO, bn);
                     //double xSW = gsl_vector_get(SW, bn);
                     //double xPCGO = INT_POL(TSL, TPCGO, xSL, sizeof(TSL)/sizeof(double));
                     //gsl_vector_set(PCGO, bn, 0);
                     //double xPCOW = INT_POL(TSW, TPCOW, xSW, sizeof(TSW)/sizeof(double));
                     //gsl_vector_set(PCOW, bn, 0);

                     //gsl_vector_set (PW, bn, (PhOW - gsl_vector_get(PCOW, bn)));

                     if (gsl_vector_get(SG, bn) > 0) {
                        gsl_vector_set(PB, bn, gsl_vector_get(PG, bn));
                     } else {
                        gsl_vector_set(PB, bn, PRBB);
                     }

                } else if (k < (DGOC-1)) {/** capa de gas */

                     gsl_vector_set (SW, bn, 0.5);
                     gsl_vector_set (SG, bn, 0.1);
                     gsl_vector_set (SO, bn, 0.4);

                     double Bgx = INT_POL(TPVT, TBG, PhGO, sizeof(TPVT)/sizeof(double));
                     double RhoGx = RHO_FVF(RHOG, Bgx);
                     double PhGp = IT_RFPRES(Dk[DGOC-1], RhoGx, Dk[k]);

                     //double xSL = gsl_vector_get(SW, bn) + gsl_vector_get(SO, bn);
                     //double xSW = gsl_vector_get(SW, bn);
                     //double xPCGO = INT_POL(TSL, TPCGO, xSL, sizeof(TSL)/sizeof(double));
                     //gsl_vector_set(PCGO, bn, 0);
                     //double xPCOW = INT_POL(TSW, TPCOW, xSW, sizeof(TSW)/sizeof(double));
                     //gsl_vector_set(PCOW, bn, 0);
                     //gsl_vector_set (PG, bn, (PhGO + gsl_vector_get(PCGO, bn) + PhGp));
                     //gsl_vector_set (PW, bn, (PhGO - gsl_vector_get(PCOW, bn) + PhGp));
                     //gsl_vector_set (PO, bn, gsl_vector_get(PG, bn));
                     gsl_vector_set (PO, bn, 10500);
                     gsl_vector_set (PW, bn, 10500);
                     gsl_vector_set (PG, bn, 10500);

                     if (gsl_vector_get(SG, bn) > 0) {
                        gsl_vector_set(PB, bn, gsl_vector_get(PG, bn));
                     } else {
                        gsl_vector_set(PB, bn, PRBB);
                     }

                } else if (k > (DWOC-1)) {/** zona de agua */

                     gsl_vector_set (SW, bn, 0.5);
                     gsl_vector_set (SG, bn, 0.1);
                     gsl_vector_set (SO, bn, 0.4);

                     double RhoWx = RHO_WX(BWI, CW, RHOW, REFPVTSTOCK, PhOW);
                     double PhWp = IT_RFPRES(Dk[DWOC-1], RhoWx, Dk[k]);
                     gsl_vector_set (PO, bn, 10500);
                     gsl_vector_set (PW, bn, 10500);
                     gsl_vector_set (PG, bn, 10500);
                     //gsl_vector_set (PW, bn, PhOW + PhWp);
                     //gsl_vector_set (PO, bn, gsl_vector_get(PW, bn));
                     //gsl_vector_set (PG, bn, gsl_vector_get(PW, bn));

                     if (gsl_vector_get(SG, bn) > 0) {
                        gsl_vector_set(PB, bn, gsl_vector_get(PG, bn));
                     } else {
                        gsl_vector_set(PB, bn, PRBB);
                     }

                } else if (k > (DGOC-1) && k < (DWOC-1)) {/** zona de oleo */

                     gsl_vector_set (SW, bn, 0.5);
                     gsl_vector_set (SG, bn, 0.1);
                     gsl_vector_set (SO, bn, 0.4);

                     double xBox = INT_POL(TPVT, TBO, PRBB, sizeof(TPVT)/sizeof(double));
                     double xRsx = INT_POL(TPVT, TRS, PRBB, sizeof(TPVT)/sizeof(double));
                     double RhoOx = RHO_OOUS(RHOO, xBox, CO, REFPRES, PRBB);
                     double RhoGx = RHO_GOUS(xRsx, RHOG, xBox, CO, REFPRES, PRBB);

                     double PhOp = IT_RFPRES(REFDP, (RhoOx+RhoGx), Dk[k]);

                     //double xSL = gsl_vector_get(SW, bn) + gsl_vector_get(SO, bn);
                     //double xSW = gsl_vector_get(SW, bn);
                     //double xPCGO = INT_POL(TSL, TPCGO, xSL, sizeof(TSL)/sizeof(double));
                     //gsl_vector_set(PCGO, bn, 0);
                     //double xPCOW = INT_POL(TSW, TPCOW, xSW, sizeof(TSW)/sizeof(double));
                     //gsl_vector_set(PCOW, bn, 0);
                     //gsl_vector_set (PO, bn, (REFPRES + PhOp));
                     //gsl_vector_set (PG, bn, (REFPRES + PhOp - gsl_vector_get(PCGO, bn)));
                     //gsl_vector_set (PW, bn, (REFPRES + PhOp - gsl_vector_get(PCOW, bn)));
                     gsl_vector_set (PO, bn, 10500);
                     gsl_vector_set (PW, bn, 10500);
                     gsl_vector_set (PG, bn, 10500);

                     if (gsl_vector_get(SG, bn) > 0) {
                        gsl_vector_set(PB, bn, gsl_vector_get(PG, bn));
                     } else {
                        gsl_vector_set(PB, bn, PRBB);
                     }

                }

                bn++; /** numero do bloco */
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Equilibrium Vertical Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}
