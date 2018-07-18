#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

extern int Gi, Gj, Gk,
           Li, Lj, Lk, Tt,
           Is, Js, Ks,
           Di[10],
           Dj[10],
           Dk[5],
           GMIJK,
           DWOC, DGOC; //Coord em Bloco
extern double POR, PERMI, PERMJ, PERMK, CR,
           RHOW, BWI, CW, VISW,
           RHOO, CO,
           RHOG,
           REFPVTSTOCK, REFDP, REFPRES, atm,
           PRBB, ATM,
           TSW[4], TKRW[4], TKROW[4], TPCOW[4],
           TSL[11], TKRG[11], TKROG[11], TPCGO[11],
           TPVT[8], TRS[8], TBO[8], TBG[8],
           TVISO[8], TVISG[8];
#endif // MODEL_H_INCLUDED
