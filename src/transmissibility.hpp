#ifndef TRANSMISSIBILITY_HPP_INCLUDED
#define TRANSMISSIBILITY_HPP_INCLUDED

void GEO_REC_CTER_BLCK(int Gi, int Gj, int Gk, int Li, int Lj, int Lk,
                         gsl_vector * KI, gsl_vector * KJ, gsl_vector * KK,
                         gsl_vector * GTI, gsl_vector * GTJ, gsl_vector * GTK);

void TRANS_1PW_PRESS(int Gi, int Gj, int Gk,
                gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                gsl_vector * BO, gsl_vector * BG, gsl_vector * BW,
                gsl_vector * RS,
                gsl_vector * VISO, gsl_vector * VISG, double VISW,
                gsl_vector * TFP_OI, gsl_vector * TFP_OJ, gsl_vector * TFP_OK,
                gsl_vector * TFP_GI, gsl_vector * TFP_GJ, gsl_vector * TFP_GK,
                gsl_vector * TFP_WI, gsl_vector * TFP_WJ, gsl_vector * TFP_WK,
                gsl_vector * TFP_ORSI, gsl_vector * TFP_ORSJ, gsl_vector * TFP_ORSK);

void TRANS_1PW_SAT(int Gi, int Gj, int Gk,
                gsl_vector * KRO, gsl_vector * KRG, gsl_vector * KRW,
                gsl_vector * TFS_OI, gsl_vector * TFS_OJ, gsl_vector * TFS_OK,
                gsl_vector * TFS_GI, gsl_vector * TFS_GJ, gsl_vector * TFS_GK,
                gsl_vector * TFS_WI, gsl_vector * TFS_WJ, gsl_vector * TFS_WK);

#endif // TRANSMISSIBILITY_HPP_INCLUDED
