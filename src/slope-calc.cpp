#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <utility>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include "functions.hpp"

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

void SLOPE_CALC(int Gi, int Gj, int Gk,
                gsl_vector * ZB,
                gsl_vector * PO, gsl_vector * PO_V,
                gsl_vector * PG, gsl_vector * PG_V,
                gsl_vector * PW, gsl_vector * PW_V,
                gsl_vector * BPOR, gsl_vector * BPOR_V, gsl_vector * SLPOR,
                gsl_vector * BO, gsl_vector * BO_V, gsl_vector * SLBO,
                gsl_vector * BG, gsl_vector * BG_V, gsl_vector * SLBG,
                gsl_vector * BW, gsl_vector * BW_V, gsl_vector * SLBW,
                gsl_vector * RS, gsl_vector * RS_V,  gsl_vector * SLRS)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    /** for slope*/
    for (int sbn = 1; sbn<(Gijk+1); sbn++){/** for slope*/
        /** pg 211 - 236 */
        /** Slope Calc*/

        double DFPO;
        if ((gsl_vector_get(PO_V, sbn) - gsl_vector_get(PO, sbn)) != 0){
            DFPO = 1 / (gsl_vector_get(PO_V, sbn) - gsl_vector_get(PO, sbn));
        } else { DFPO = 0; }

        double DFPW;
        if ((gsl_vector_get(PW_V, sbn) - gsl_vector_get(PW, sbn)) != 0) {
            DFPW = 1 / (gsl_vector_get(PW_V, sbn) - gsl_vector_get(PW, sbn));
        } else { DFPW = 0; }

        double DFPG;
        if ((gsl_vector_get(PG_V, sbn) - gsl_vector_get(PG, sbn)) != 0) {
            DFPG = 1 / (gsl_vector_get(PG_V, sbn) - gsl_vector_get(PG, sbn));
        } else { DFPG = 0; }

        double DFPOR;
        if (gsl_vector_get(PO_V, sbn) >= gsl_vector_get(PG_V, sbn) && gsl_vector_get(PO_V, sbn) >= gsl_vector_get(PW_V, sbn)) {
            DFPOR = DFPO;
        } else if (gsl_vector_get(PW_V, sbn) >= gsl_vector_get(PO_V, sbn) && gsl_vector_get(PW_V, sbn) >= gsl_vector_get(PG_V, sbn)) {
            DFPOR = DFPW;
        } else { DFPOR = DFPG; }

        double SLPOR_ITER = (gsl_vector_get(BPOR_V, sbn) - gsl_vector_get(BPOR, sbn)) * DFPOR;

        double SLBO_ITER;
        if (gsl_vector_get(BO_V, sbn) != 0 || gsl_vector_get(BO, sbn) != 0 ) {
            SLBO_ITER = ((1/gsl_vector_get(BO_V, sbn)) - (1/gsl_vector_get(BO, sbn))) * DFPO;
        } else { SLBO_ITER = 0; }

        double SLBG_ITER;
        if (gsl_vector_get(BG_V, sbn) != 0 || gsl_vector_get(BG, sbn) != 0 ) {
            SLBG_ITER = ((1/gsl_vector_get(BG_V, sbn)) - (1/gsl_vector_get(BG, sbn))) * DFPG;
        } else { SLBG_ITER = 0; }

        double SLBW_ITER;
        if (gsl_vector_get(BW_V, sbn) != 0 || gsl_vector_get(BW, sbn) != 0 ) {
            SLBW_ITER = ((1/gsl_vector_get(BW_V, sbn)) - (1/gsl_vector_get(BW, sbn))) * DFPW;
        } else { SLBW_ITER = 0; }

        double SLRS_ITER;
        if (gsl_vector_get(RS_V, sbn) != 0 || gsl_vector_get(RS, sbn) != 0 ) {
            SLRS_ITER = (gsl_vector_get(RS_V, sbn) - gsl_vector_get(RS, sbn)) * DFPO;
        } else { SLRS_ITER = 0; }

        gsl_vector_set(SLPOR, sbn, SLPOR_ITER);
        gsl_vector_set(SLBO, sbn, SLBO_ITER);
        gsl_vector_set(SLBG, sbn, SLBG_ITER);
        gsl_vector_set(SLBW, sbn, SLBW_ITER);
        gsl_vector_set(SLRS, sbn, SLRS_ITER);

    }/** end for slope */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Slope Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}/**Full Implicit */
