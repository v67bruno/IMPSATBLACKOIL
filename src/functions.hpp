#ifndef FUNCTIONS_HPP_INCLUDED
#define FUNCTIONS_HPP_INCLUDED

double INT_POL(double *xx, double *yy, double point, int n);
double RHO_OOUS(double RHOSC, double BOPB, double CO, double P, double PB);
double RHO_GOUS(double RSPB, double RHOGS, double BOPB, double CO, double P, double PB);
double RF_PRES(double PINI, double ZI, double Rho, double ZB);
double RHO_FVF(double Rho, double FVF);
double IT_RFPRES(double ZI, double Rho, double ZB);
double RHO_WX(double Bwi, double Cw, double Rhowsc, double Prx, double P);
double XBO_PB(double BoPb, double CO, double P_REF, double P);
double XFVF_PR(double Bx, double Cx, double Prx, double P);
double POR_P(double POR, double CR, double P_REF, double P);
std::tuple<int, int> FASE_POT(int I, int IP, gsl_vector * PF, gsl_vector * RHF, gsl_vector * ZB);
double FASE_POT_GET(int I, int IP, gsl_vector * PF, gsl_vector * RHF, gsl_vector * ZB);
std::tuple<int, int ,int> SAT_SORT(double SO, double SG, double SW);
void RESID_2_SPM_TRIP(int IJK, int PM,
                      double RWSW, double RWSG, double RWPO,
                      double RGSW, double RGSG, double RGPO,
                      double ROSW, double ROSG, double ROPO,
                      gsl_matrix * JM);
void RESID_2_VEC_WGO(int IJK, double WN, double GN, double ON, gsl_vector * VEC);
void PROP_2_VEC_WGP(int Gi, int Gj, int Gk,
                    gsl_vector * SWN, gsl_vector * SGN, gsl_vector * PON,
                    gsl_vector * VEC);
std::tuple<double, double> TRANS_NM(int UPS, int DWS, gsl_vector * TG, gsl_vector * TFP_N, gsl_vector * TFP_V, gsl_vector * TFS_N, gsl_vector * TFS_V);
std::tuple<double, double> TRANS_NM_2(int UPS, int DWS, gsl_vector * TG, gsl_vector * TFP_V, gsl_vector * TFS_V);
double UPS_TRANS_NM(int UPS, int DWS, gsl_vector * TG, gsl_vector * TFP, gsl_vector * TFS);
double DS_TRANS_NM(double DV, int UPS, int DWS, gsl_vector * TG,
                   gsl_vector * TFP_V, gsl_vector * TFS_V,
                   gsl_vector * TFP_N, gsl_vector * TFS_N);
double DP_TRANS_NM(double DV, int UPS, int DWS, gsl_vector * TG,
                   gsl_vector * TFP_V, gsl_vector * TFS_V,
                   gsl_vector * TFP_N, gsl_vector * TFS_N);
void SRT_VEC_2_OGW(double PRBB, gsl_vector * SW_N, gsl_vector * SG_N, gsl_vector * PO_N, gsl_vector * XV,
                   gsl_vector * SO_V, gsl_vector * SG_V, gsl_vector * SW_V,
                   gsl_vector * PO_V, gsl_vector * PG_V, gsl_vector * PW_V, gsl_vector * PB_V);
void UPD_SAT_PS(gsl_vector * X_N, double PRBB,
                gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                gsl_vector * PB);

void NEW_PROP_NP(gsl_vector * X_NM, gsl_vector * X_N, gsl_vector * X_VP);

#endif // FUNCTIONS_HPP_INCLUDED
