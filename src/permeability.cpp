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

void REL_PERM_DIST(int Gi, int Gj, int Gk,
                   gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                   gsl_vector * KRO, gsl_vector * KRG, gsl_vector * KRW)
{
    TimeVar t1=timeNow();
    double Gijk = Gi*Gj*Gk;
    for (int bn=1; bn<(Gijk+1); bn++){ /** Perm do bloco e fluidos */
        //std::cout << "krg" << bn << std::endl;
        double xSO = gsl_vector_get(SO, bn);
        double xSG = gsl_vector_get(SG, bn);
        double xSW = gsl_vector_get(SW, bn);
        double xKRW = INT_POL(TSW, TKRW, xSW, sizeof(TSW)/sizeof(double));
        //std::cout << "krg" << bn << std::endl;
        double xKROW = INT_POL(TSW, TKROW, xSW, sizeof(TSW)/sizeof(double));
        double xKRG = INT_POL(TSL, TKRG, (xSO+xSW), sizeof(TSL)/sizeof(double));

        double xKROG = INT_POL(TSL, TKROG, (xSO+xSW), sizeof(TSL)/sizeof(double));
        double xKRC = INT_POL(TSW, TKROW, TSW[0], sizeof(TSW)/sizeof(double));
        //std::cout << "kr" << bn << std::endl;
        double xA = xKROW/xKRC;
        double xB = xKROG/xKRC;
        double xC = xA + xKRW;
        double xD = xB + xKRG;
        double xE = xC * xD;
        double xF = xE - xKRW - xKRG;
        double xKRO = xKRC * xF;

        gsl_vector_set (KRO, bn, xKRO);
        gsl_vector_set (KRG, bn, xKRG);
        gsl_vector_set (KRW, bn, xKRW);

    } /** Termino Calc Perm */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Distribuicao de Kr Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
} /** Perm Distrib end */

void PERM_DIST(int Gi, int Gj, int Gk,
               double PERMI, double PERMJ, double PERMK,
               gsl_vector * KI, gsl_vector * KJ, gsl_vector * KK)
{
    TimeVar t1=timeNow();
    double Gijk = Gi*Gj*Gk;
    for (int bn=1; bn<(Gijk+1); bn++){ /** Perm do bloco e fluidos */

        gsl_vector_set (KI, bn, PERMI);
        gsl_vector_set (KJ, bn, PERMJ);
        gsl_vector_set (KK, bn, PERMK);

    } /** Termino Calc Perm */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Distribuicao de K Calculado ";
    std::cout << std::setprecision(20);
    std::cout << "| Tempo Exec.: " << duration(timeNow()-t1)/1e+9 << " s" << std::endl;
} /** Perm Distrib end */
