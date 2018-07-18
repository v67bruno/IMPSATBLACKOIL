#include <iostream>
#include <iomanip>
#include<cmath>
#include <chrono>
#include <utility>
#include<gsl/gsl_math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_spmatrix.h>
#include "model.hpp"
#include "functions.hpp"

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

void FLUID_PROP (int Gi, int Gj, int Gk, int * Dk,
                 double PRBB,
                 gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                 gsl_vector * PO, gsl_vector * PG, gsl_vector * PW, gsl_vector * PB,
                 gsl_vector * BO, gsl_vector * BG, gsl_vector * BW, gsl_vector * RS,
                 gsl_vector * VISO, gsl_vector * VISG,
                 gsl_vector * RHO, gsl_vector * RHG, gsl_vector * RHW,
                 gsl_vector * BPOR)
{
    TimeVar t1=timeNow();

    for (int bn=1; bn<(GMIJK+1); bn++){ /** Propriedades dos Fluidos */

        double Pbx = gsl_vector_get(PB, bn);
        double Pgx = gsl_vector_get(PG, bn);
        double Pox = gsl_vector_get(PO, bn);
        double Pwx = gsl_vector_get(PW, bn);
        double Sox = gsl_vector_get(SO, bn);
        double Swx = gsl_vector_get(SW, bn);
        double Sgx = gsl_vector_get(SG, bn);

        if (Pbx >= PRBB) {/** PO >= Pb Undersaturated */
            /** SG < 1e-2, Undersaturated*/
            /** sg < 1% -> Pb = Pb */
            if (Sgx < 0.0001) {

                if (Pgx > 0){
                    double fBgx = INT_POL(TPVT, TBG, Pgx, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (BG, bn, fBgx);
                    double xVsG = INT_POL(TPVT, TVISG, Pgx, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (VISG, bn, xVsG);
                    gsl_vector_set (RHG, bn, (RHO_FVF(RHOG, gsl_vector_get(BG, bn))));
                }

                if (Pox > 0){
                    double fBoPb = INT_POL(TPVT, TBO, Pbx, sizeof(TPVT)/sizeof(double));
                    double fBoX = XBO_PB(fBoPb, CO, Pbx, Pox);
                    gsl_vector_set (BO, bn, fBoX);

                    double fRsX = INT_POL(TPVT, TRS, Pbx, sizeof(TPVT)/sizeof(double));
                    fRsX *= (Sgx/(Sgx+0.0001));
                    gsl_vector_set (RS, bn, fRsX);

                    double doo = RHO_OOUS(RHOO, gsl_vector_get(BO, bn), CO, REFPRES, Pbx);
                    double dgo = RHO_GOUS(gsl_vector_get(RS, bn), RHOG, gsl_vector_get(BO, bn), CO, REFPRES, Pbx);
                    gsl_vector_set (RHO, bn, (doo+dgo));
                    double xVsO = INT_POL(TPVT, TVISO, Pox, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (VISO, bn, xVsO);
                }

                if (Pwx > 0){
                    double fBwX = XFVF_PR(BWI, CW, REFPVTSTOCK, Pwx);
                    gsl_vector_set (BW, bn, fBwX);
                    double fRhwX = RHO_WX(BWI, CW, RHOW, REFPVTSTOCK, gsl_vector_get(PW, bn));
                    gsl_vector_set (RHW, bn, fRhwX);
                }

                if (Pgx > Pox){
                    double xPor = POR_P(POR, CR, REFPVTSTOCK, Pgx);
                    gsl_vector_set (BPOR, bn, xPor);

                } else if (Pwx > Pox) {
                    double xPor = POR_P(POR, CR, REFPVTSTOCK, Pwx);
                    gsl_vector_set (BPOR, bn, xPor);

                } else {
                    double xPor = POR_P(POR, CR, REFPVTSTOCK, Pox);
                    gsl_vector_set (BPOR, bn, xPor);
                }
            /** sg > 1% -> Pb = Pf */
            /** Saturated Sg > 0 */
            } else {

                if (Pgx > 0){
                    double fBgx = INT_POL(TPVT, TBG, Pgx, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (BG, bn, fBgx);
                    double xVsG = INT_POL(TPVT, TVISG, Pgx, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (VISG, bn, xVsG);
                    gsl_vector_set (RHG, bn, (RHO_FVF(RHOG, gsl_vector_get(BG, bn))));
                }

                if (Pox > 0){
                    double fBoPb = INT_POL(TPVT, TBO, Pox, sizeof(TPVT)/sizeof(double));
                    double fBoX = XBO_PB(fBoPb, CO, Pbx, Pox);
                    gsl_vector_set (BO, bn, fBoX);
                    double fRsX = INT_POL(TPVT, TRS, Pox, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (RS, bn, fRsX);
                    double doo = RHO_FVF(RHOO, gsl_vector_get(BO, bn));
                    double dgo = RHO_FVF((RHOG * gsl_vector_get(RS, bn)), gsl_vector_get(BO, bn));
                    gsl_vector_set (RHO, bn, (doo+dgo));
                    double xVsO = INT_POL(TPVT, TVISO, Pox, sizeof(TPVT)/sizeof(double));
                    gsl_vector_set (VISO, bn, xVsO);
                }

                if (Pwx > 0){
                    double fBwX = XFVF_PR(BWI, CW, REFPVTSTOCK, Pwx);
                    gsl_vector_set (BW, bn, fBwX);
                    double fRhwX = RHO_WX(BWI, CW, RHOW, REFPVTSTOCK, Pwx);
                    gsl_vector_set (RHW, bn, fRhwX);
                }

                if (Pgx > Pox){
                    double xPor = POR_P(POR, CR, REFPVTSTOCK, Pgx);
                    gsl_vector_set (BPOR, bn, xPor);

                } else if (Pwx > Pox) {
                    double xPor = POR_P(POR, CR, REFPVTSTOCK, Pwx);
                    gsl_vector_set (BPOR, bn, xPor);

                } else {
                    double xPor = POR_P(POR, CR, REFPVTSTOCK, Pox);
                    gsl_vector_set (BPOR, bn, xPor);
                }

            }/** sg check */

        } else {/** PO < Pb Saturated */

            if (Pgx > 0){
                double fBgx = INT_POL(TPVT, TBG, Pgx, sizeof(TPVT)/sizeof(double));
                gsl_vector_set (BG, bn, fBgx);
                double xVsG = INT_POL(TPVT, TVISG, Pgx, sizeof(TPVT)/sizeof(double));
                gsl_vector_set (VISG, bn, xVsG);
                gsl_vector_set (RHG, bn, (RHO_FVF(RHOG, gsl_vector_get(BG, bn))));
            }

            if (Pox > 0){
                double fBoPb = INT_POL(TPVT, TBO, Pox, sizeof(TPVT)/sizeof(double));
                double fBoX = XBO_PB(fBoPb, CO, Pbx, Pox);
                gsl_vector_set (BO, bn, fBoX);
                double fRsX = INT_POL(TPVT, TRS, Pox, sizeof(TPVT)/sizeof(double));
                gsl_vector_set (RS, bn, fRsX);
                double doo = RHO_FVF(RHOO, gsl_vector_get(BO, bn));
                double dgo = RHO_FVF((RHOG * gsl_vector_get(RS, bn)), gsl_vector_get(BO, bn));
                gsl_vector_set (RHO, bn, (doo+dgo));
                double xVsO = INT_POL(TPVT, TVISO, Pox, sizeof(TPVT)/sizeof(double));
                gsl_vector_set (VISO, bn, xVsO);
            }

            if (Pwx > 0){
                double fBwX = XFVF_PR(BWI, CW, REFPVTSTOCK, Pwx);
                gsl_vector_set (BW, bn, fBwX);
                double fRhwX = RHO_WX(BWI, CW, RHOW, REFPVTSTOCK, Pwx);
                gsl_vector_set (RHW, bn, fRhwX);
            }

            if (Pgx > Pox){
                double xPor = POR_P(POR, CR, REFPVTSTOCK, Pgx);
                gsl_vector_set (BPOR, bn, xPor);

            } else if (Pwx > Pox) {
                double xPor = POR_P(POR, CR, REFPVTSTOCK, Pwx);
                gsl_vector_set (BPOR, bn, xPor);

            } else {
                double xPor = POR_P(POR, CR, REFPVTSTOCK, Pox);
                gsl_vector_set (BPOR, bn, xPor);
            }
        } /** Termina for */

    } /** Termina for */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Propriedades das Fases Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;

} /** Termina funcao */
