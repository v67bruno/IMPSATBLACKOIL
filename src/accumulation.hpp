#ifndef ACCUMULATION_HPP_INCLUDED
#define ACCUMULATION_HPP_INCLUDED

void ACCUMULATION_RESID(int tt, int Gi, int Gj, int Gk, int Is, int Js, int Ks,
                   gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                   gsl_vector * BPOR, gsl_vector * BO, gsl_vector * BG,
                   gsl_vector * BW, gsl_vector * RS, gsl_vector * C_WGO);

void ACCUMULATION_DIF(int tt, int Gi, int Gj, int Gk, int Is, int Js, int Ks,
                   gsl_vector * SO_N, gsl_vector * SO_V, gsl_vector * SG_N,
                   gsl_vector * SG_V, gsl_vector * SW_N, gsl_vector * SW_V,
                   gsl_vector * BPOR_N, gsl_vector * BPOR_V, gsl_vector * SLPOR,
                   gsl_vector * BO_N, gsl_vector * BO_V, gsl_vector * SLBO,
                   gsl_vector * BG_N, gsl_vector * BG_V, gsl_vector * SLBG,
                   gsl_vector * BW_N, gsl_vector * BW_V, gsl_vector * SLBW,
                   gsl_vector * RS_N, gsl_vector * RS_V,  gsl_vector * SLRS,
                   gsl_matrix * CM_WGO);

#endif // ACCUMULATION_HPP_INCLUDED
