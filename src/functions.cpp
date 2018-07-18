#include <cmath>
#include <tuple>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include "functions.hpp"

using namespace std;

double INT_POL(double *xx, double *yy, double point, int n)
{
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
    gsl_interp * interpolation = gsl_interp_alloc (gsl_interp_linear,n);
    gsl_interp_init(interpolation, xx, yy, n);
    double VINTPOL = gsl_interp_eval(interpolation, xx, yy, point, accelerator);
    gsl_interp_free(interpolation);
    gsl_interp_accel_free(accelerator);
    return VINTPOL;
}

double RHO_OOUS(double RHOSC, double BOPB, double CO, double P, double PB)
{
    double RSCPB = RHOSC/BOPB;
    double COPPB = 1 + (CO * (P - PB));
    double VSRSPB = RSCPB * COPPB;
    return VSRSPB;
}

double RHO_GOUS(double RSPB, double RHOGS, double BOPB, double CO, double P, double PB)
{
    double RGPB = (RSPB*RHOGS)/BOPB;
    double COPPB = 1 + (CO * (P - PB));
    double VSGSPB = RGPB * COPPB;
    return VSGSPB;
}

double RF_PRES(double PINI, double ZI, double Rho, double ZB)
{
    /**kPa*/
    double NPRES = PINI + (Rho * 9.8066352) * (ZB - ZI)/1000;
    return NPRES;
}

double RHO_FVF(double Rho, double FVF)
{
    double RB = Rho/FVF;
    return RB;
}

double IT_RFPRES(double ZI, double Rho, double ZB)
{
    double ITPRES = Rho * 9.8066352 * (ZB - ZI)/1000;
    return ITPRES;
}

double RHO_WX(double Bwi, double Cw, double Rhowsc, double Prx, double P)
{
    double Bpxi = (Rhowsc/Bwi) * (1 + Cw * (P - Prx));
    return Bpxi;
}

double XBO_PB(double BoPb, double CO, double P_REF, double P)
{
    double XBOPB = BoPb * (1 - CO * (P - P_REF));
    return XBOPB;
}

double XFVF_PR(double Bx, double Cx, double Prx, double P)
{
    double Bpxi = Bx * (1 - Cx * (P - Prx));
    return Bpxi;
}

double POR_P(double POR, double CR, double P_REF, double P)
{
    double XPOR = POR * (1 + CR * (P - P_REF));
    return XPOR;
}

std::tuple<int, int> FASE_POT(int I, int IP,
                gsl_vector * PF, gsl_vector * RHF, gsl_vector * ZB)
{
    int UPS, DWS;
    double DF_PF = gsl_vector_get(PF, IP) - gsl_vector_get(PF, I);
    double GMF_M = 9.8066352 * 0.5 * (gsl_vector_get(RHF, IP) + gsl_vector_get(RHF, I));
    double PH_M = GMF_M * (gsl_vector_get(ZB, IP) - gsl_vector_get(ZB, I));

    double FSPT = DF_PF - PH_M;

    if (FSPT < 0) {
        UPS = I;
        DWS = IP;
    } else if (FSPT > 0) {
        UPS = IP;
        DWS = I;
    } else if (FSPT == 0) {
        if (I > IP) {UPS = IP; DWS = I;}
        else {UPS = I; DWS = IP;}
    }

    return std::make_tuple(UPS,DWS);
}

double FASE_POT_GET(int I, int IP,
                    gsl_vector * PF, gsl_vector * RHF, gsl_vector * ZB)
{
    double DF_PF = gsl_vector_get(PF, IP) - gsl_vector_get(PF, I);
    double GMF_M = 9.8066352 * 0.5 * (gsl_vector_get(RHF, IP) + gsl_vector_get(RHF, I));
    double PH_M = GMF_M * (gsl_vector_get(ZB, IP) - gsl_vector_get(ZB, I));

    double FSPT = DF_PF - PH_M;
    return FSPT;
}

std::tuple<int, int ,int> SAT_SORT(double SO, double SG, double SW)
{

    double maxx, medd, minn;
    int ica, ice, ici;

    if (SG > SO && SG > SW) {maxx = SG;  ica = 1;}
    else if (SG > SO && SG < SW) {medd = SG;  ice = 1;}
    else if (SG < SO && SG > SW) {medd = SG;  ice = 1;}
    else if (SG < SO && SG < SW) {minn = SG; ici = 1;}

    if (SO > SG && SO > SW) {maxx = SO;  ica = 0;}
    else if (SO > SG && SO < SW) {medd = SO;  ice = 0;}
    else if (SO < SG && SO > SW) {medd = SO;  ice = 0;}
    else if (SO < SG && SO < SW) {minn = SO;  ici = 0;}

    if (SW > SG && SW > SO) {maxx = SW;  ica = 2;}
    else if (SW > SG && SW < SO) {medd = SW;  ice = 2;}
    else if (SW < SG && SW > SO) {medd = SW;  ice = 2;}
    else if (SW < SG && SW < SO) {minn = SW;  ici = 2;}

    return std::make_tuple(ica, ice, ici);

}

void RESID_2_SPM_TRIP(int IJK, int PM,
                      double RWSW, double RWSG, double RWPO,
                      double RGSW, double RGSG, double RGPO,
                      double ROSW, double ROSG, double ROPO,
                      gsl_matrix * JM)
{

    /** ijk to matrix */
    int fl = 0; //fluid control
    for (int bk = ((IJK*3)-3); bk<((IJK*3)); bk++){
        /** W */
        if (fl == 0 ){
            int unk = 0; //unknow control
            for (int fk = ((PM*3)-3); fk<((PM*3)); fk++){
                if (unk == 0){gsl_matrix_set(JM,bk,fk,RWSW);}
                if (unk == 1){gsl_matrix_set(JM,bk,fk,RWSG);}
                if (unk == 2){gsl_matrix_set(JM,bk,fk,RWPO);}
                //unknow control
                unk++;
            }
        }
        /** G */
        if (fl == 1 ){
            int unk = 0; //unknow control
            for (int fk = ((PM*3)-3); fk<((PM*3)); fk++){
                if (unk == 0){gsl_matrix_set(JM,bk,fk,RGSW);}
                if (unk == 1){gsl_matrix_set(JM,bk,fk,RGSG);}
                if (unk == 2){gsl_matrix_set(JM,bk,fk,RGPO);}
                //unknow control
                unk++;
            }
        }
        /** O */
        if (fl == 2 ){
            int unk = 0; //unknow control
            for (int fk = ((PM*3)-3); fk<((PM*3)); fk++){
                if (unk == 0){gsl_matrix_set(JM,bk,fk,ROSW);}
                if (unk == 1){gsl_matrix_set(JM,bk,fk,ROSG);}
                if (unk == 2){gsl_matrix_set(JM,bk,fk,ROPO);}
                //unknow control
                unk++;
            }
        }
        //fluid control
        fl++;
        // main for
    }/** ijk to matrix */

    //return gsl_spmatrix * JM;

}

void RESID_2_VEC_WGO(int IJK, double WN, double GN, double ON,
                     gsl_vector * VEC)
{

    int fnl = 0;
    for (int ik = ((IJK*3)-3); ik<((IJK*3)); ik++){
        if (fnl == 0){gsl_vector_set(VEC, ik, WN);}
        if (fnl == 1){gsl_vector_set(VEC, ik, GN);}
        if (fnl == 2){gsl_vector_set(VEC, ik, ON);}
        fnl++;
    }

    //return gsl_vector * VEC;
}

void PROP_2_VEC_WGP(int Gi, int Gj, int Gk,
                    gsl_vector * SWN, gsl_vector * SGN, gsl_vector * PON,
                    gsl_vector * VEC)
{
    int bn = 1;
    for (int k = 0; k<Gk; k++){
        for (int j = 0; j<Gj; j++){
            for (int i = 0; i<Gi; i++){

                double WN, GN, ON;
                WN = gsl_vector_get(SWN,bn);
                GN = gsl_vector_get(SGN,bn);
                ON = gsl_vector_get(PON,bn);

                int fnl = 0;
                for (int ik = ((bn*3)-3); ik<((bn*3)); ik++){
                    if (fnl == 0){gsl_vector_set(VEC, ik, WN);}
                    if (fnl == 1){gsl_vector_set(VEC, ik, GN);}
                    if (fnl == 2){gsl_vector_set(VEC, ik, ON);}
                    fnl++;
                }
                bn++;
            }
        }
    }
}

std::tuple<double, double> TRANS_NM(int UPS, int DWS, gsl_vector * TG, gsl_vector * TFP_N, gsl_vector * TFP_V, gsl_vector * TFS_N, gsl_vector * TFS_V)
{
    double T_N, T_V, D_T;
    T_N = gsl_vector_get(TG, UPS);
    T_N *= (0.5 * (gsl_vector_get(TFP_N, UPS) + gsl_vector_get(TFP_N, DWS)));
    T_N *= gsl_vector_get(TFS_N, UPS);

    T_V = gsl_vector_get(TG, UPS);
    T_V *= (0.5 * (gsl_vector_get(TFP_V, UPS) + gsl_vector_get(TFP_V, DWS)));
    T_V *= gsl_vector_get(TFS_V, UPS);

    D_T = T_V - T_N;

    return std::make_tuple(T_V, D_T);
}

std::tuple<double, double> TRANS_NM_2(int UPS, int DWS, gsl_vector * TG, gsl_vector * TFP_V, gsl_vector * TFS_V)
{
    double T_U, T_D, D_T;
    T_U = gsl_vector_get(TG, UPS);
    T_U *= (0.5 * (gsl_vector_get(TFP_V, UPS) + gsl_vector_get(TFP_V, DWS)));
    T_U *= gsl_vector_get(TFS_V, UPS);

    T_D = gsl_vector_get(TG, DWS);
    T_D *= (0.5 * (gsl_vector_get(TFP_V, UPS) + gsl_vector_get(TFP_V, DWS)));
    T_D *= gsl_vector_get(TFS_V, DWS);

    D_T = T_U - T_D;

    return std::make_tuple(T_U, D_T);
}

double UPS_TRANS_NM(int UPS, int DWS, gsl_vector * TG, gsl_vector * TFP, gsl_vector * TFS)
{
    double T_U = gsl_vector_get(TG, UPS);
    T_U *= (0.5 * (gsl_vector_get(TFP, UPS) + gsl_vector_get(TFP, DWS)));
    T_U *= gsl_vector_get(TFS, UPS);

    return T_U;
}

double DS_TRANS_NM(double DV, int UPS, int DWS, gsl_vector * TG,
                   gsl_vector * TFP_V, gsl_vector * TFS_V,
                   gsl_vector * TFP_N, gsl_vector * TFS_N)
{
    double T_V, T_N, D_V;

    T_V = gsl_vector_get(TG, UPS);
    T_V *= (0.5 * (gsl_vector_get(TFP_V, UPS) + gsl_vector_get(TFP_V, DWS)));
    T_V *= gsl_vector_get(TFP_V, UPS);

    T_N = gsl_vector_get(TG, UPS);
    T_N *= (0.5 * (gsl_vector_get(TFP_N, UPS) + gsl_vector_get(TFP_N, UPS)));
    T_N *= gsl_vector_get(TFP_N, UPS);

    D_V = abs(T_V - T_N)/DV;

    return D_V;

}

double DP_TRANS_NM(double DV, int UPS, int DWS, gsl_vector * TG,
                   gsl_vector * TFP_V, gsl_vector * TFS_V,
                   gsl_vector * TFP_N, gsl_vector * TFS_N)
{
    double T_V, T_N, D_V;

    T_V = gsl_vector_get(TG, UPS);
    T_V *= (0.5 * (gsl_vector_get(TFP_V, UPS) + gsl_vector_get(TFP_V, DWS)));
    T_V *= gsl_vector_get(TFP_V, UPS);

    T_N = gsl_vector_get(TG, UPS);
    T_N *= (0.5 * (gsl_vector_get(TFP_N, UPS) + gsl_vector_get(TFP_N, UPS)));
    T_N *= gsl_vector_get(TFP_N, UPS);

    D_V = abs(T_V - T_N)/DV;

    return D_V;
}

void SRT_VEC_2_OGW(double PRBB, gsl_vector * SW_N, gsl_vector * SG_N, gsl_vector * PO_N, gsl_vector * XV,
                   gsl_vector * SO_V, gsl_vector * SG_V, gsl_vector * SW_V,
                   gsl_vector * PO_V, gsl_vector * PG_V, gsl_vector * PW_V, gsl_vector * PB_V)
{
    size_t ij=1;
    for(size_t i = 0; i<(((XV->size)/3)); i++){

        //printf("%d,%d\n",i,ij);
        double nPO = gsl_vector_get(PO_N, ij);
        double nSW = gsl_vector_get(SW_N, ij);
        double nSG = gsl_vector_get(SG_N, ij);
        double rPB, rSO;
        double rSW = gsl_vector_get(XV, i*3);
        double rSG = gsl_vector_get(XV, ((i*3)+1));
        double rPO = gsl_vector_get(XV, ((i*3)+2));

        if (rSW < 0.25) {rSW = nSW;}
        else if (rSW > 1) {rSW = nSW;}
        else {rSW = rSW;}

        if (rPO < 0) {rPO = nPO;}
        else if (rPO > 11000) {rPO = nPO;}
        else {rPO = rPO;}

        if ((nPO == PRBB)&&(rPO > PRBB)) /// Under -> Under
        {
            if (rSG < 0.0001) { rSG = pow((nSG*0.0001),0.5); }
            if (rSG > 0.7499) { rSG = nSG; }
            else { rSG = rSG; }
            rSO = 1 - rSW - rSG;
            rPB = PRBB;
        }
        if ((nPO == PRBB)&&(rPO < PRBB)) /// Under -> Sat
        {
            if (rSG < 0.0001) { rSG = pow((nSG*0.0001),0.5); }
            if (rSG > 0.7499) { rSG = nSG; }
            else { rSG = rSG; }
            rSO = 1 - rSW - rSG;
            rPB = rPO;
        }
        if ((nPO < PRBB)&&(rPO > PRBB)) /// Sat -> Under
        {
            if (rSG < 0.0001) { rSG = 0.00001; }
            if (rSG > 0.7499) { rSG = nSG; }
            else { rSG = 0.00001; }
            rSO = 1 - rSW - rSG;
            rPB = PRBB;
        }
        if ((nPO < PRBB)&&(rPO < PRBB)) /// Sat - Sat
        {
            if (rSG < 0.0001) { rSG = pow((nSG*0.0001),0.5); }
            if (rSG > 0.7499) { rSG = nSG; }
            else { rSG = rSG; }
            rSO = 1 - rSW - rSG;
            rPB = rPO;
        }
        if ((nPO > PRBB)&&(rPO < PRBB)) /// Under -> Sat
        {
            if (rSG < 0.0001) { rSG = pow((nSG*0.0001),0.5); }
            if (rSG > 0.7499) { rSG = nSG; }
            else { rSG = rSG; }
            rSO = 1 - rSW - rSG;
            rPB = rPO;
        }
        if ((nPO > PRBB)&&(rPO > PRBB)) /// Under - Under
        {
            if (rSG < 0.0001) { rSG = 0.00001; }
            if (rSG > 0.7499) { rSG = pow((nSG*0.0001),0.5); }
            else { rSG = rSG; }
            rSO = 1 - rSW - rSG;
            rPB = PRBB;
        }

        gsl_vector_set(PB_V, ij, rPB);
        gsl_vector_set(SW_V, ij, rSW);
        gsl_vector_set(SG_V, ij, rSG);
        gsl_vector_set(SO_V, ij, rSO);
        gsl_vector_set(PO_V, ij, rPO);
        gsl_vector_set(PG_V, ij, rPO);
        gsl_vector_set(PW_V, ij, rPO);

        ++ij;

    }
}

void NEW_PROP_NP(gsl_vector * X_NM, gsl_vector * X_N, gsl_vector * X_VP)
{
    for (size_t i = 0; i<(((X_NM->size)/3)); i++) {
        double nSW, nSG, nPO,
               vSW, vSG, vPO,
               SW, SG, PO;

        nSW = gsl_vector_get(X_NM, (i*3));
        nSG = gsl_vector_get(X_NM, ((i*3)+1));
        nPO = gsl_vector_get(X_NM, ((i*3)+2));
        vSW = gsl_vector_get(X_VP, (i*3));
        vSG = gsl_vector_get(X_VP, ((i*3)+1));
        vPO = gsl_vector_get(X_VP, ((i*3)+2));

        if ((vSW < 0.25)||(vSW > 1)) {SW = nSW;}
        else {SW = vSW;}

        if ((vSG < 0)||(vSG > 1)) {SG = nSG;}
        else {SG = vSG;}

        if ((vPO < 0)||(vPO > 11000)) {PO = nPO;}
        else {PO = vPO;}

        gsl_vector_set(X_N,(i*3),SW);
        gsl_vector_set(X_N,((i*3)+1),SG);
        gsl_vector_set(X_N,((i*3)+2),PO);
    }
}

void UPD_SAT_PS(gsl_vector * X_N, double PRBB,
                gsl_vector * SO, gsl_vector * SG, gsl_vector * SW,
                gsl_vector * PO, gsl_vector * PG, gsl_vector * PW,
                gsl_vector * PB)
{
    size_t ij = 1;
    for(size_t i = 0; i<(((X_N->size)/3)); i++){

        double nSO, nSG, nSW, nPO, nPB;
        nSW = gsl_vector_get(X_N, (i*3));
        nSG = gsl_vector_get(X_N, ((i*3)+1));
        nPO = gsl_vector_get(X_N, ((i*3)+2));
        nSO = 1 - nSW - nSG;

        if (nPO > PRBB) {
            if (nSG > 0.0001) { nPB = nPO; }
            else { nPB = PRBB; }
        }
        else { nPB = nPO; }

        gsl_vector_set(SO, ij, nSO);
        gsl_vector_set(SG, ij, nSG);
        gsl_vector_set(SW, ij, nSW);
        gsl_vector_set(PO, ij, nPO);
        gsl_vector_set(PG, ij, nPO);
        gsl_vector_set(PW, ij, nPO);
        gsl_vector_set(PB, ij, nPB);

        ++ij;
    }
}
/**void VEC_2_PROP(int Gi, int Gj, int Gk,
                gsl_vector * VEC,
                gsl_vector * SO, gsl_vector * SW, gsl_vector * SG,
                gsl_vector * PO, gsl_vector * PW, gsl_vector * PG)
{
    size_t ij = 1;
    for(size_t i = 0; i<(((VEC->size)/3)); i++){

        double SW = gsl_vector_get(VEC, i*3);
        double SG = gsl_vector_get(VEC, ((i*3)+1));
        double PO = gsl_vector_get(VEC, ((i*3)+2));


        gsl_vector_set(PW_V, ij, rPO);

        ++ij;

    }
}*/
