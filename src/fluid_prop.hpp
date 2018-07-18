#ifndef FLUID_PROP_HPP_INCLUDED
#define FLUID_PROP_HPP_INCLUDED

void FLUID_PROP (int Gi, int Gj, int Gk, int * Dk,
                 double PRBB,
                 gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                 gsl_vector * PO, gsl_vector * PG, gsl_vector * PW, gsl_vector * PB,
                 gsl_vector * BO, gsl_vector * BG, gsl_vector * BW, gsl_vector * RS,
                 gsl_vector * VISO, gsl_vector * VISG,
                 gsl_vector * RHO, gsl_vector * RHG, gsl_vector * RHW,
                 gsl_vector * BPOR);

#endif // FLUID_PROP_HPP_INCLUDED
