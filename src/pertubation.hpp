#ifndef PERTUBATION_HPP_INCLUDED
#define PERTUBATION_HPP_INCLUDED

void PERTUB_PRES(int Gi, int Gj, int Gk,
                 double PRBB, double P_Per,
                 gsl_vector * SG,
                 gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                 gsl_vector * PB, gsl_vector * PCGO, gsl_vector * PCOW,
                 gsl_vector * PO_V, gsl_vector * PG_V, gsl_vector * PW_V,
                 gsl_vector * PB_V, gsl_vector * PCGO_V, gsl_vector * PCOW_V);

void PERTUB_SAT(int Gi, int Gj, int Gk,
                double GC_W, double GO_W, double OZ_W, double OW_W, double WZ_W,
                double GC_G, double GO_G, double OZ_G, double OW_G, double WZ_G,
                gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                gsl_vector * SO_V, gsl_vector * SG_V, gsl_vector * SW_V);

#endif // PERTUBATION_HPP_INCLUDED
