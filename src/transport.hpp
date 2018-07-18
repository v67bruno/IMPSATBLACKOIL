#ifndef TRANSPORT_HPP_INCLUDED
#define TRANSPORT_HPP_INCLUDED

void TRANSPORT_1PW_2(int Gi, int Gj, int Gk, gsl_vector * D_NV,
                   gsl_vector * ZB,
                   gsl_vector * SO_V, gsl_vector * SW_V, gsl_vector * SG_V,
                   gsl_vector * PO_V, gsl_vector * PW_V, gsl_vector * PG_V, gsl_vector * PCGO_V, gsl_vector * PCOW_V,
                   gsl_vector * RHO_N, gsl_vector * RHG_N, gsl_vector * RHW_N,
                   gsl_vector * RHO_V, gsl_vector * RHG_V, gsl_vector * RHW_V,
                   gsl_vector * RS_N, gsl_vector * RS_V,
                   gsl_vector * TGI, gsl_vector * TGJ, gsl_vector * TGK,
                gsl_vector * TFS_OI_N, gsl_vector * TFS_OJ_N, gsl_vector * TFS_OK_N,
                gsl_vector * TFS_GI_N, gsl_vector * TFS_GJ_N, gsl_vector * TFS_GK_N,
                gsl_vector * TFS_WI_N, gsl_vector * TFS_WJ_N, gsl_vector * TFS_WK_N,
                gsl_vector * TFP_OI_N, gsl_vector * TFP_OJ_N, gsl_vector * TFP_OK_N,
                gsl_vector * TFP_GI_N, gsl_vector * TFP_GJ_N, gsl_vector * TFP_GK_N,
                gsl_vector * TFP_WI_N, gsl_vector * TFP_WJ_N, gsl_vector * TFP_WK_N,
                gsl_vector * TFP_ORSI_N, gsl_vector * TFP_ORSJ_N, gsl_vector * TFP_ORSK_N,
                gsl_vector * TFS_OI_V, gsl_vector * TFS_OJ_V, gsl_vector * TFS_OK_V,
                gsl_vector * TFS_GI_V, gsl_vector * TFS_GJ_V, gsl_vector * TFS_GK_V,
                gsl_vector * TFS_WI_V, gsl_vector * TFS_WJ_V, gsl_vector * TFS_WK_V,
                gsl_vector * TFP_OI_V, gsl_vector * TFP_OJ_V, gsl_vector * TFP_OK_V,
                gsl_vector * TFP_GI_V, gsl_vector * TFP_GJ_V, gsl_vector * TFP_GK_V,
                gsl_vector * TFP_WI_V, gsl_vector * TFP_WJ_V, gsl_vector * TFP_WK_V,
                gsl_vector * TFP_ORSI_V, gsl_vector * TFP_ORSJ_V, gsl_vector * TFP_ORSK_V,
                gsl_vector * FNV, gsl_matrix * JACOBIAN);
#endif // TRANSPORT_HPP_INCLUDED
