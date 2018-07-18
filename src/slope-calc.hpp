#ifndef SLOPE_CALC_HPP_INCLUDED
#define SLOPE_CALC_HPP_INCLUDED

void SLOPE_CALC(int Gi, int Gj, int Gk,
                gsl_vector * ZB,
                gsl_vector * PO, gsl_vector * PO_V,
                gsl_vector * PG, gsl_vector * PG_V,
                gsl_vector * PW, gsl_vector * PW_V,
                gsl_vector * BPOR, gsl_vector * BPOR_V, gsl_vector * SLPOR,
                gsl_vector * BO, gsl_vector * BO_V, gsl_vector * SLBO,
                gsl_vector * BG, gsl_vector * BG_V, gsl_vector * SLBG,
                gsl_vector * BW, gsl_vector * BW_V, gsl_vector * SLBW,
                gsl_vector * RS, gsl_vector * RS_V,  gsl_vector * SLRS);

#endif // SLOPE-CALC_HPP_INCLUDED
