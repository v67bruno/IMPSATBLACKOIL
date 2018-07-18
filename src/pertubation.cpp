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

void PERTUB_SAT(int Gi, int Gj, int Gk,
                double GC_W, double GO_W, double OZ_W, double OW_W, double WZ_W,
                double GC_G, double GO_G, double OZ_G, double OW_G, double WZ_G,
                gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                gsl_vector * SO_V, gsl_vector * SG_V, gsl_vector * SW_V)
{
    TimeVar t1=timeNow();

    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    int bn = 1;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {
                /** cantos*/
                double SO_b =  gsl_vector_get(SO, bn);
                double SG_b =  gsl_vector_get(SG, bn);
                double SW_b =  gsl_vector_get(SW, bn);
                double SO_p, SG_p, SW_p;
                /** SO = 0, SG = 1, SW = 2 */
                auto [smax, smed, smin] = SAT_SORT(SO_b, SG_b, SW_b);
                //printf("B[%d]\nSat Sort: %d, %d, %d\n",bn,smax, smed, smin);

                if (k == 0) {
                    SW_p = SW_b - GC_W;
                    SG_p = SG_b + GC_G;
                }
                if (k == 1) {
                    SW_p = SW_b - GO_W;
                    SG_p = SG_b + GO_G;
                }
                if (k == 2) {
                    SW_p = SW_b - OZ_W;
                    SG_p = SG_b - OZ_G;
                }
                if (k == 3) {
                    SW_p = SW_b + OW_W;
                    SG_p = SG_b - OW_G;
                }
                if (k == 4) {
                    SW_p = SW_b + WZ_W;
                    SG_p = SG_b - WZ_G;
                }

                SO_p = 1 - SW_p - SG_p;

                //printf("P: SO: %f, SG: %f, SW: %f\n",SO_p,SG_p,SW_p);
                //printf("B: SO: %f, SG: %f, SW: %f\n\n",SO_b,SG_b,SW_b);

                gsl_vector_set(SO_V, bn, SO_p);
                gsl_vector_set(SG_V, bn, SG_p);
                gsl_vector_set(SW_V, bn, SW_p);

                bn++;

                //SW_p = 0; SG_p = 0; SO_p = 0;
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Pertubacao Saturacao Calc. ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}

void PERTUB_PRES(int Gi, int Gj, int Gk,
                 double PRBB, double P_Per,
                 gsl_vector * SG,
                 gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                 gsl_vector * PB, gsl_vector * PCGO, gsl_vector * PCOW,
                 gsl_vector * PO_V, gsl_vector * PG_V, gsl_vector * PW_V,
                 gsl_vector * PB_V, gsl_vector * PCGO_V, gsl_vector * PCOW_V)
{
    TimeVar t1=timeNow();

    int Gijk = Gi*Gj*Gk, Gij = Gi*Gj;
    int bn = 1;
    for (int k = 0; k<Gk; k++) {
        for (int j = 0; j<Gj; j++) {
            for (int i = 0; i<Gi; i++) {

                double PO_b = gsl_vector_get(PO, bn);
                //double PG_b = gsl_vector_get(PG, bn);
                //double PW_b = gsl_vector_get(PW, bn);
                double PB_b = gsl_vector_get(PB, bn);
                double SG_b = gsl_vector_get(SG, bn);
                double PO_p, PG_p, PW_p, PB_p;

                PO_p = PO_b + P_Per;

                if (PB_b > PRBB) {
                    if (SG_b > 0) {
                        PB_p = PO_p;
                    } else {
                        PB_p = PRBB;
                    }
                } else {
                    PB_p = PO_p;
                }

                PG_p = PO_p; PW_p = PO_p;

                gsl_vector_set(PO_V, bn, PO_p);
                gsl_vector_set(PG_V, bn, PG_p);
                gsl_vector_set(PW_V, bn, PW_p);
                gsl_vector_set(PB_V, bn, PB_p);

                bn++;
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Pertubacao Pressao Calc. ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
}
