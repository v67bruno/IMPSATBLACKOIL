#ifndef WELL_HPP_INCLUDED
#define WELL_HPP_INCLUDED

void WELL_PROD(int Gi, int Gj, int Gk,
               int CI, int CJ, int CK, double WI, double BHP, double RATE,
               gsl_vector * D_NV,
               gsl_vector * PO_N, gsl_vector * SW_N, gsl_vector * SG_N,
               gsl_vector * KRO_N, gsl_vector * KRG_N, gsl_vector * KRW_N,
               gsl_vector * BO_N, gsl_vector * BG_N, gsl_vector * BW_N,
               gsl_vector * RS_N, gsl_vector * VISO_N, gsl_vector * VISG_N,
               gsl_vector * PO_V, gsl_vector * SW_V, gsl_vector * SG_V,
               gsl_vector * KRO_V, gsl_vector * KRG_V, gsl_vector * KRW_V,
               gsl_vector * BO_V, gsl_vector * BG_V, gsl_vector * BW_V,
               gsl_vector * RS_V, gsl_vector * VISO_V, gsl_vector * VISG_V,
               double VISW, gsl_vector * QV, gsl_matrix * QRN);

void WELL_INJ_W(int Gi, int Gj, int Gk,
               int CI, int CJ, int CK, double WI, double BHP, double RATE,
               gsl_vector * PO_N, gsl_vector * PO_V,
               gsl_vector * SW_N, gsl_vector * SW_V,
               gsl_vector * SG_N, gsl_vector * SG_V,
               gsl_vector * KRO_N, gsl_vector * KRO_V,
               gsl_vector * KRG_N, gsl_vector * KRG_V,
               gsl_vector * KRW_N, gsl_vector * KRW_V,
               gsl_vector * BO_N, gsl_vector * BO_V,
               gsl_vector * BG_N, gsl_vector * BG_V,
               gsl_vector * BW_N, gsl_vector * BW_V,
               gsl_vector * RS_N, gsl_vector * RS_V,
               gsl_vector * VISO_N, gsl_vector * VISO_V,
               gsl_vector * VISG_N, gsl_vector * VISG_V,
               double VISW,
               gsl_vector * QV,
               gsl_matrix * QRN);
#endif // WELL_HPP_INCLUDED
