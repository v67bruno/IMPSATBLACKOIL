#ifndef EQUILIBRIUM_HPP_INCLUDED
#define EQUILIBRIUM_HPP_INCLUDED


void EQUILIBRIUM (int Gi, int Gj, int Gk, int Ks,
                  int DWOC, int DGOC, int * Dk,
                  double REFDP, double REFPRES, double PRBB,
                  double * TSW, double * TSL, double * TPCGO, double * TPCOW,
                  gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                  gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                  gsl_vector * PB, gsl_vector * PCGO, gsl_vector * PCOW,
                  gsl_vector * ZB);

#endif // EQUILIBRIUM_HPP_INCLUDED
