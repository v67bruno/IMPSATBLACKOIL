#ifndef PERMEABILITY_HPP_INCLUDED
#define PERMEABILITY_HPP_INCLUDED

void REL_PERM_DIST(int Gi, int Gj, int Gk,
                   gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                   gsl_vector * KRO, gsl_vector * KRG, gsl_vector * KRW);

void PERM_DIST(int Gi, int Gj, int Gk,
               double PERMI, double PERMJ, double PERMK,
               gsl_vector * KI, gsl_vector * KJ, gsl_vector * KK);

#endif // PERMEABILITY_HPP_INCLUDED
