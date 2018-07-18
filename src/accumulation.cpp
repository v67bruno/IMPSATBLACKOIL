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

void ACCUMULATION_RESID(int tt, int Gi, int Gj, int Gk, int Is, int Js, int Ks,
                   gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                   gsl_vector * BPOR, gsl_vector * BO, gsl_vector * BG,
                   gsl_vector * BW, gsl_vector * RS, gsl_vector * C_WGO)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    /** for acumulation */
    for (int bn = 1; bn<(Gijk+1); bn++){ /** for acumulation */
        double VbTt = Is * Js * Ks / tt;
        double CON_ITER, CGN_ITER, CWN_ITER;

        CWN_ITER = gsl_vector_get(BPOR, bn) * gsl_vector_get(SW, bn) / gsl_vector_get(BW, bn);
        CWN_ITER *= VbTt;

        CGN_ITER = gsl_vector_get(BPOR, bn) * gsl_vector_get(RS, bn) * gsl_vector_get(SO, bn) / gsl_vector_get(BO, bn);
        CGN_ITER += (gsl_vector_get(BPOR, bn) * gsl_vector_get(SG, bn) / gsl_vector_get(BG, bn));
        CGN_ITER *= VbTt;

        CON_ITER = gsl_vector_get(BPOR, bn) * gsl_vector_get(SO, bn) / gsl_vector_get(BO, bn);
        CON_ITER *= VbTt;

        if (gsl_finite(CWN_ITER) !=1){CWN_ITER = 0;}
        if (gsl_finite(CGN_ITER) !=1){CGN_ITER = 0;}
        if (gsl_finite(CON_ITER) !=1){CON_ITER = 0;}

        RESID_2_VEC_WGO(bn, CWN_ITER, CGN_ITER, CON_ITER, C_WGO);

    }/** end for acumulation */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Accumulation Resid Terms Calc. ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}/**Full Implicit */

void ACCUMULATION_DIF(int tt, int Gi, int Gj, int Gk, int Is, int Js, int Ks,
                   gsl_vector * SO_N, gsl_vector * SO_V, gsl_vector * SG_N,
                   gsl_vector * SG_V, gsl_vector * SW_N, gsl_vector * SW_V,
                   gsl_vector * BPOR_N, gsl_vector * BPOR_V, gsl_vector * SLPOR,
                   gsl_vector * BO_N, gsl_vector * BO_V, gsl_vector * SLBO,
                   gsl_vector * BG_N, gsl_vector * BG_V, gsl_vector * SLBG,
                   gsl_vector * BW_N, gsl_vector * BW_V, gsl_vector * SLBW,
                   gsl_vector * RS_N, gsl_vector * RS_V,  gsl_vector * SLRS,
                   gsl_matrix * CM_WGO)
{
    TimeVar t1=timeNow();
    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    /** for acumulation */
    for (int bn = 1; bn<(Gijk+1); bn++){ /** for acumulation */

        double COP_ITER, COW_ITER, COG_ITER, CWP_ITER, CWW_ITER, CWG_ITER, CGP_ITER, CGW_ITER, CGG_ITER, CGP_RS_ITER, CGP_G_ITER;
        double VbTt = Is * Js * Ks / tt;

        CWP_ITER = gsl_vector_get(SLPOR, bn) / gsl_vector_get(BW_N, bn);
        CWP_ITER += (gsl_vector_get(BPOR_V, bn) * gsl_vector_get(SLBW, bn));
        CWP_ITER *= gsl_vector_get(SW_N, bn) * VbTt;

        CWW_ITER = gsl_vector_get(BPOR_V, bn) / gsl_vector_get(BW_V, bn);
        CWW_ITER *= VbTt;

        CWG_ITER = 0;

        if (gsl_finite(CWP_ITER) !=1){CWP_ITER = 0;}
        if (gsl_finite(CWW_ITER) !=1){CWW_ITER = 0;}
        if (gsl_finite(CWG_ITER) !=1){CWG_ITER = 0;}

        CGP_RS_ITER = gsl_vector_get(SLPOR, bn) / gsl_vector_get(BO_N, bn);
        CGP_RS_ITER += (gsl_vector_get(BPOR_V, bn) * gsl_vector_get(SLBO, bn));
        CGP_RS_ITER *= gsl_vector_get(RS_N, bn);
        CGP_RS_ITER += (gsl_vector_get(SLRS, bn) * gsl_vector_get(BPOR_V, bn) / gsl_vector_get(BO_V, bn));
        CGP_RS_ITER *= gsl_vector_get(SO_N, bn);
        CGP_G_ITER = gsl_vector_get(SLPOR, bn) / gsl_vector_get(BG_N, bn);
        CGP_G_ITER += (gsl_vector_get(BPOR_V, bn) * gsl_vector_get(SLBG, bn));
        CGP_G_ITER *= gsl_vector_get(SG_N, bn);
        CGP_ITER = CGP_RS_ITER + CGP_G_ITER;

        CGW_ITER = gsl_vector_get(BPOR_V, bn) / gsl_vector_get(BO_V, bn);
        CGW_ITER *= VbTt * gsl_vector_get(RS_V, bn);
        CGW_ITER *= - 1;

        CGG_ITER = gsl_vector_get(BPOR_V, bn) / gsl_vector_get(BG_V, bn);
        CGG_ITER -= (gsl_vector_get(RS_V, bn) * gsl_vector_get(BPOR_V, bn) / gsl_vector_get(BO_V, bn));
        CGG_ITER *= VbTt;

        if (gsl_finite(CGP_ITER) !=1){CGP_ITER = 0;}
        if (gsl_finite(CGW_ITER) !=1){CGW_ITER = 0;}
        if (gsl_finite(CGG_ITER) !=1){CGG_ITER = 0;}

        COP_ITER = gsl_vector_get(SLPOR, bn) / gsl_vector_get(BO_N, bn);
        COP_ITER += (gsl_vector_get(BPOR_V, bn) * gsl_vector_get(SLBO, bn));
        COP_ITER *= gsl_vector_get(SO_N, bn) * VbTt;

        COW_ITER = gsl_vector_get(BPOR_V, bn) / gsl_vector_get(BO_V, bn);
        COW_ITER *= VbTt;
        COW_ITER *= -1;

        COG_ITER = COW_ITER;

        if (gsl_finite(COP_ITER) !=1){COP_ITER = 0;}
        if (gsl_finite(COW_ITER) !=1){COW_ITER = 0;}
        if (gsl_finite(COG_ITER) !=1){COG_ITER = 0;}

        RESID_2_SPM_TRIP(bn, bn,
                         CWW_ITER, CWG_ITER, CWP_ITER,
                         CGW_ITER, CGG_ITER, CGP_ITER,
                         COW_ITER, COG_ITER, COP_ITER,
                         CM_WGO);

    }/** end for acumulation */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Accumulation Dif Terms Calc. ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}/**Full Implicit */
