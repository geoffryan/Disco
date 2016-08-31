#include <stdio.h>
#include <math.h>
#include "paul.h"
#include "Hydro/metric.h"
#include "Hydro/frame.h"

#define ACC 1e-10
#define MAX_ITER 10

enum{DD,EN,MX,MY,MZ,BX,BY,BZ};

static double GAMMA_LAW = 0.0;

//Global Functions
void get_Ustar_HLLD(double w, double *pL, double *pR, double *F, double *U, 
                    double *x, double *n);

//Local Functions
void solve_HLLD_SR(double sL, double sR, double Bx, double *UL, double *FR,
                    double *UR, double *FR, double w, double *U, double *F);
void calc_va(double p, double s, double Bx, double *R, double *va, 
                double *dva);
void calc_Ba(double s, double Bx, double *R, double *va, double *Ba, 
                double *dva, double *dBa);
void calc_entha(double p, double s, double *R, double *va, double *enth,
                    double *dva, double *denth);
void calc_Ka(double p, double s, double Bx, double enth, double *R, double *K,
                double denth, double *dKa, int LR);
void calc_Bc(double Bx, double sL, double sR, double *BaL, double *BaR, 
                double *vaL, double *vaR, double saL, double saR, double *Bc,
                double *dBaL, double *dBaR, double *dvaL, double *dvaR,
                double dsaL, double dsaR, double *dBc);
void cont_func(double Bx, double *BaL, double *BaR,
                double *vaL, double *vaR, double *KaL, double *KaR,
                double enthL, double enthR, double *f,
                double *dBaL, double *dBaR, double *dvaL, double *dvaR,
                double *dKaL, double *dKaR, double denthL, double denthR,
                double *df);
void calc_vc(double Bx, double *Ka, double *Bc, double enth, double *vc,
                int LR);
int calc_Ps(double *RL, double *RR, double Bx, double sL, double sR, 
                double *Ps);
void calc_state(double p, double *RL, double *RR, double Bx, double sL, 
                    double sR, double *vaL, double *vaR, double *BaL, 
                    double *BaR, double *enthL, double *enthR, double *KaL,
                    double *KaR, double *Bc, double *vcL, double *vcR);
void calc_Ua(double p, double s, double Bx, double *R, double *va, double *Ba,
                double *Ua);
void calc_Uc(double p, double sa, double Bx, double *va, double *vc, 
                double *Bc, double *Ua, double *Uc);

void setHlldParams( struct domain * theDomain )
{
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
}

void get_Ustar_HLLD(double w, double *pL, double *pR, double *F, double *U, 
                    double *x, double *n)
{
    double r = x[0];
    double al, be[3], gam[9], igam[9], jac;
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    double sqrtgam = jac / al;

    double rhoL = pL[RHO];
    double PL = pL[PPP];
    double lL[3] = {pL[URR], pL[UPP], pL[UZZ]};
    double BL[3] = {pL[BRR]/sqrtgam,pL[BPP]/(r*sqrtgam),pL[BZZ]/sqrtgam};

    double rhoR = pR[RHO];
    double PR = pR[PPP];
    double lR[3] = {pR[URR], pR[UPP], pR[UZZ]};
    double BR[3] = {pR[BRR]/sqrtgam,pR[BPP]/(r*sqrtgam),pR[BZZ]/sqrtgam};

    double U[4];
    frame_U(x, U);

    //Pseudo code
    L,R = transform(pL, pR);

    Ps = getP(L, R);

    saL, sC, saR = calc_vel(Ps, L, R);

    Us, Fs = flux_for_region(saL, sC, saR, w);

    U, F = transform_back(Us, Fs);

}

void solve_HLLD_SR(double sL, double sR, double Bx, double *UL, double *FR,
                    double *UR, double *FR, double w, double *U, double *F)
{
    double Ps;

    double RL[NUM_Q], RR[NUM_Q], Uhll[NUM_Q], Fhll[NUM_Q];
    int q;
    for(q = 0; q<NUM_Q; q++)
    {
        RL[q] = sL*UL[q] - FL[q];
        RR[q] = sR*UR[q] - FR[q];
        Uhll[NUM_Q] = (sR*UR[q] - sL*UL[q] + FL[q]-FR[q]) / (sR - sL);
        Fhll[NUM_Q] = (sR*FL[q] - sL*FR[q] + sL*sR*(UR[q]-UL[q])) / (sR - sL);
    }

    int err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);

    double BaL[3], Bc[3], BaR[3];
    double vaL[3], vcL[3], vcR[3], vaR[3];
    double KL[3], KR[3];
    double enthL, enthR;
    double saL, sc, saR;

    calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, &enthL, &enthR, 
                KL, KR, Bc, vcL, vcR);

    saL = KL[0];
    saR = KR[0];
    sc = 0.5*(vcL[0] + vcR[0]);

    if(w < sL)
    {
        for(q=0; q<NUM_Q; q++)
        {
            U[q] = UL[q];
            F[q] = FL[q];
        }
    }
    else if (w < saL)
    {
        calc_Ua(p, sL, Bx, RL, vaL, BaL, U);
        for(q=0; q<NUM_Q; q++)
            F[q] = FL[q] + sL*(U[q]-UL[q]);
    }
    else if(w < sc)
    {
        double UaL[NUM_Q];
        calc_Ua(p, sL, Bx, RL, vaL, BaL, UaL);
        calc_Uc(p, saL, Bx, vaL, vcL, Bc, UaL, U);
        for(q=0; q<NUM_Q; q++)
            F[q] = FL[q] + sL*(UaL[q]-UL[q]) + saL*(U[q]-UaL[q]);
    }
    else if(w < saR)
    {
        double UaR[NUM_Q];
        calc_Ua(p, sR, Bx, RR, vaR, BaR, UaR);
        calc_Uc(p, saR, Bx, vaR, vcR, Bc, UaR, U);
        for(q=0; q<NUM_Q; q++)
            F[q] = FR[q] + sR*(UaR[q]-UR[q]) + saR*(U[q]-UaR[q]);
    }
    else if(w < sR)
    {
        calc_Ua(p, sR, Bx, RR, vaR, BaR, U);
        for(q=0; q<NUM_Q; q++)
            F[q] = FR[q] + sR*(UR[q]-UR[q]);
    }
    else
    {
        for(q=0; q<NUM_Q; q++)
        {
            U[q] = UR[q];
            F[q] = FR[q];
        }
    }
}

void calc_va(double p, double s, double Bx, double *R, double *va, double *dva)
{
    //Eq. (23-25) of MUB
    double A = R[MX] - s*R[EN] + p*(1.0-s*s);
    double G = R[BY]*R[BY] + R[BZ]*R[BZ];
    double C = R[MY]*R[BY] + R[MZ]*R[BZ];
    double Q = -A - G + Bx*Bx*(1.0-s*s);
    double iX = 1.0/(Bx * (s*A*Bx+C) - (A+G) * (s*p+R[EN]));

    va[0] = iX * (Bx * (A*Bx+s*C) - (A+G) * (p+R[MX]));
    va[1] = iX * (Q*R[MY] + R[BY] * (C+Bx*(s*R[MX]-R[EN])));
    va[2] = iX * (Q*R[MZ] + R[BZ] * (C+Bx*(s*R[MX]-R[EN])));

    if(dva != NULL)
    {
        double dA = 1.0-s*s;
        double dQ = -1.0;
        double diX = iX*iX*s*(A+G);
        dva[0] = diX*(Bx * (A*Bx+s*C) - (A+G) * (p+R[MX]))
                    + iX*(Bx*Bx*dA - dA*(p+R[MX]) - (A+G));
        dva[1] = diX * (Q*R[MY] + R[BY] * (C+Bx*(s*R[MX]-R[EN])))
                    + iX*(dQ*R[MY]);
        dva[2] = diX * (Q*R[MZ] + R[BZ] * (C+Bx*(s*R[MX]-R[EN])))
                    + iX*(dQ*R[MZ]);
    }
}

void calc_Ba(double s, double Bx, double *R, double *va, double *Ba, 
                double *dva, double *dBa)
{
    //Eq. (21) of MUB
    Ba[0] = Bx;
    Ba[1] = (R[BY] - Bx*va[1]) / (s - va[0]);
    Ba[2] = (R[BZ] - Bx*va[2]) / (s - va[0]);

    if(dBA != NULL)
    {
        double ids2 = 1.0/((s-va[0])*(s-va[0]));
        dBa[0] = 0.0;
        dBa[1] = (-Bx*dva[1]*(s-va[0]) + (R[BY]-Bx*va[1])*dva[0]) * ids2;
        dBa[2] = (-Bx*dva[2]*(s-va[0]) + (R[BZ]-Bx*va[2])*dva[0]) * ids2;
    }
}

void calc_entha(double p, double s, double *R, double *va, double *enth,
                    double *dva, double *denth)
{
    //Eq. (31) of MUB
    *enth =  p + (R[EN] - va[0]*R[MX]-va[1]*R[MY]-va[2]*R[MZ]) / (s - va[0]);

    if(denth != NULL)
        *denth = 1.0 + ((-dva[0]*R[MX]-dva[1]*R[MY]-dva[2]*R[MZ])*(s-va[0])
                        + (R[EN]-va[0]*R[MX]-va[1]*R[MY]-va[2]*R[MZ])*dva[0])
                        / ((s-va[0])*(s-va[0]));
}

void calc_Ka(double p, double s, double Bx, double enth, double *R, double *K,
                double denth, double *dKa, int LR)
{
    //Eq. (43) of MUB
    //LR == 1 for Right, LR == -1 for Left
    
    int sgn = Bx > 0 ? LR : -LR;
    double eta = sgn*sqrt(enth);
    double denom = 1.0 / (s*p + R[EN] + Bx*eta)

    Ka[0] = (R[MX] + p + R[BX]*eta) * denom;
    Ka[1] = (R[MY]     + R[BY]*eta) * denom;
    Ka[2] = (R[MZ]     + R[BZ]*eta) * denom;

    if(*dKa != NULL)
    {
        double deta = 0.5/eta * denth;
        double ddenom = -denom*denom*(s + Bx*deta);

        dKa[0] = (1.0 + R[BX]*deta)*denom + (R[MX]+p+R[BX]*eta)*ddenom;
        dKa[1] = (R[BY]*deta)*denom + (R[MY]+R[BY]*eta)*ddenom;
        dKa[2] = (R[BZ]*deta)*denom + (R[MZ]+R[BZ]*eta)*ddenom;
    } 
}

void calc_Bc(double Bx, double sL, double sR, double *BaL, double *BaR, 
                double *vaL, double *vaR, double saL, double saR, double *Bc,
                double *dBaL, double *dBaR, double *dvaL, double *dvaR,
                double dsaL, double dsaR, double *dBc)
{
    //Eq. (45) of MUB

    Bc[0] = Bx;
    Bc[1] = ((BaR[1]*(saR-vaR[0])+Bx*vaR[1])
                - (BaL[1]*(saL-vaL[0]) + Bx*vaL[1]))
            / (saR - saL);
    Bc[2] = ((BaR[2]*(saR-vaR[0])+Bx*vaR[2])
                - (BaL[2]*(saL-vaL[0]) + Bx*vaL[2]))
            / (saR - saL);

    if(dBc != NULL)
    {
        denom = 1.0 / (saR - saL);
        ddenom = -denom*denom*(dsaR - dsaL);

        dBc[0] = 0.0;
        dBc[1] = ((dBaR[1]*(saR-vaR[0]) - BaR[1]*(dsaR-dvaR[0]) + Bx*dvaR[1])
                - (dBaL[1]*(saL-vaL[0]) - BaL[1]*(dsaL-dvaL[0]) + Bx*dvaL[1]))
                 * denom + ((BaR[1]*(saR-vaR[0])+Bx*vaR[1])
                            - (BaL[1]*(saL-vaL[0]) + Bx*vaL[1])) * ddenom;
        dBc[2] = ((dBaR[2]*(saR-vaR[0]) - BaR[2]*(dsaR-dvaR[0]) + Bx*dvaR[2])
                - (dBaL[2]*(saL-vaL[0]) - BaL[2]*(dsaL-dvaL[0]) + Bx*dvaL[2]))
                 * denom + ((BaR[2]*(saR-vaR[0])+Bx*vaR[2])
                            - (BaL[2]*(saL-vaL[0]) + Bx*vaL[2])) * ddenom;
    }
}

void cont_func(double Bx, double *BaL, double *BaR,
                double *vaL, double *vaR, double *KaL, double *KaR,
                double enthL, double enthR, double *f,
                double *dBaL, double *dBaR, double *dvaL, double *dvaR,
                double *dKaL, double *dKaR, double denthL, double denthR,
                double *df)
{
    //Eq. (48) of MUB

    double saL = KaL[0];
    double saR = KaR[0];
    double dsaL = dKaL[0];
    double dsaR = dKaR[0];
    
    // Bch = \hat{B_c}, numerator of Eq. 45 for Bc. Recalculating here is more
    //  stable if B~0.
    double Bch[3];
    Bch[0] = (saR-saL)*Bx;
    Bch[1] =  (BaR[1]*(saR-vaR[0]) + Bx*vaR[1]) 
            - (BaL[1]*(saL-vaL[0]) + Bx*vaL[1]);
    Bch[2] =  (BaR[2]*(saR-vaR[0]) + Bx*vaR[2]) 
            - (BaL[2]*(saL-vaL[0]) + Bx*vaL[2]);

    double KaL2 = KaL[0]*KaL[0] + KaL[1]*KaL[1] + KaL[2]*KaL[2];
    double KaR2 = KaR[0]*KaR[0] + KaR[1]*KaR[1] + KaR[2]*KaR[2];
    double KBL = KaL[0]*Bch[0] + KaL[1]*Bch[1] + KaL[2]*Bch[2];
    double KBR = KaR[0]*Bch[0] + KaR[1]*Bch[1] + KaR[2]*Bch[2];

    int sgn = Bx>0 ? 1 : -1;
    double etaL = -sgn*sqrt(enthL);
    double etaR =  sgn*sqrt(enthR);

    double delKx = (KaR[0]-KaL[0]);

    // Eq. 49
    double YL = (1.0 - KaL2) / (etaL*delKx - KBL);
    double YR = (1.0 - KaR2) / (etaR*delKx - KBR);

    //Eq. 48
    *f = (KaR[0]-KaL[0]) * (1.0 - Bx*(YR-YL));

    if(*df != NULL)
    {
        double dBch[3];
        dBch[0] = (dsaR - dsaL) * Bx;
        dBch[1] = (dBaR[1]*(saR-vaR[1]) + BaR[1]*(dsaR-dvaR[1]) + Bx*dvaR[1])
                - (dBaL[1]*(saL-vaL[1]) + BaL[1]*(dsaL-dvaL[1]) + Bx*dvaL[1]);
        dBch[2] = (dBaR[2]*(saR-vaR[2]) + BaR[2]*(dsaR-dvaR[2]) + Bx*dvaR[2])
                - (dBaL[2]*(saL-vaL[2]) + BaL[2]*(dsaL-dvaL[2]) + Bx*dvaL[2]);
        
        double dKaL2 = 2*(KaL[0]*dKaL[0] + KaL[1]*dKaL[1] + KaL[2]*dKaL[2]);
        double dKaR2 = 2*(KaR[0]*dKaR[0] + KaR[1]*dKaR[1] + KaR[2]*dKaR[2]);
        double dKBL = dKaL[0]*Bch[0] + KaL[0]*dBch[0]
                    + dKaL[1]*Bch[1] + KaL[1]*dBch[1]
                    + dKaL[2]*Bch[2] + KaL[2]*dBch[2];
        double dKBR = dKaR[0]*Bch[0] + KaR[0]*dBch[0]
                    + dKaR[1]*Bch[1] + KaR[1]*dBch[1]
                    + dKaR[2]*Bch[2] + KaR[2]*dBch[2];

        double detaL = 0.5/etaL * denthL;
        double detaR = 0.5/etaR * denthR;

        double ddelKx = dKaR[0] - dKaL[0];
        
        double dYL = (-dKaL2*(etaL*delKx-KBL)
                        - (1-KaL2)*(detaL*delKx + etaL*ddelKx - dKBL))
                    / ((etaL*delKx-KBL)*(etaL*delKx-KBL));
        double dYR = (-dKaR2*(etaR*delKx-KBR)
                        - (1-KaR2)*(detaR*delKx + etaR*ddelKx - dKBR))
                    / ((etaR*delKx-KBR)*(etaR*delKx-KBR));

        *df = ddelKx * (1.0-Bx*(YR-YL)) + delKx * (-Bx*(dYR-dYL));
    }
}

void calc_vc(double Bx, double *Ka, double *Bc, double enth, double *vc,
                int LR)
{
    int sgn = Bx > 0 ? LR : -LR;
    double eta = sgn*sqrt(enth);

    double K2 = Ka[0]*Ka[0] + Ka[1]*Ka[1] + Ka[2]*Ka[2];
    double KB = Ka[0]*Bx + Ka[1]*Bc[1] + Ka[2]*Bc[2];

    vc[0] = Ka[0] - Bx * (1.0-K2) / (eta - KB);
    vc[1] = Ka[1] - Bc[1] * (1.0-K2) / (eta - KB);
    vc[2] = Ka[2] - Bc[2] * (1.0-K2) / (eta - KB);
}

int calc_Ps(double *RL, double *RR, double Bx, double sL, double sR, 
                double *Ps)
{
    

    int i;
    double p0, p, p1, dp;
    double f, df;

    p0 = 1.0;

    double p1 = p0;

    i = 0;
    do{
        p = p1;
        double vaL[3], vaR[3], dvaL[3], dvaR[3];
        double BaL[3], BaR[3], dBaL[3], dBaR[3];
        double enthL, enthR, denthL, denthR;
        double KaL[3], KaR[3], dKaL[3], dKaR[3];

        calc_va(p, sL, Bx, RL, vaL, dvaL);
        calc_va(p, sR, Bx, RR, vaR, dvaR);
        calc_Ba(sL, Bx, RL, vaL, BaL, dvaL, dBaL);
        calc_Ba(sR, Bx, RR, vaR, BaR, dvaR, dBaR);
        calc_entha(p, sL, RL, vaL, &enthL, dvaL, &denthL);
        calc_entha(p, sR, RR, vaR, &enthR, dvaR, &denthR);
        calc_Ka(p, sL, Bx, enthL, RL, KaL, denthL, dKaL, -1);
        calc_Ka(p, sR, Bx, enthR, RR, KaR, denthR, dKaR, +1);

        cont_func(Bx, BaL, BaR, vaL, vaR, KaL, KaR, enthL, enthR, &f
                    dBaL, dBaR, dvaL, dvaR, dKaL, dKaR, denthL, denthR, &df);

        dp = -f/df;

        p1 = p + dp;
        
    }while(fabs(dp)/p > ACC && i < MAX_ITER);

    *Ps = p1;

    return 0;
}

void calc_state(double p, double *RL, double *RR, double Bx, double sL, 
                    double sR, double *vaL, double *vaR, double *BaL, 
                    double *BaR, double *enthL, double *enthR, double *KaL,
                    double *KaR, double *Bc, double *vcL, double *vcR)
{
        calc_va(p, sL, Bx, RL, vaL, NULL);
        calc_va(p, sR, Bx, RR, vaR, NULL);
        calc_Ba(sL, Bx, RL, vaL, BaL, NULL, NULL);
        calc_Ba(sR, Bx, RR, vaR, BaR, NULL, NULL);
        calc_entha(p, sL, RL, vaL, &enthL, NULL, NULL);
        calc_entha(p, sR, RR, vaR, &enthR, NULL, NULL);
        calc_Ka(p, sL, Bx, enthL, RL, KaL, 0.0, NULL, -1);
        calc_Ka(p, sR, Bx, enthR, RR, KaR, 0.0, NULL, +1);
        calc_Bc(Bx, BaL, BaR, vaL, vaR, KaL[0], KaR[0], Bc, 
                NULL,NULL,NULL,NULL,0.0,0.0,NULL);
        calc_vc(Bx, KaL, Bc, enthL, vcL, -1);
        calc_vc(Bx, KaR, Bc, enthR, vcR, +1);
}

void calc_Ua(double p, double s, double Bx, double *R, double *va, double *Ba,
                double *Ua)
{
    //Eq. (32-34) of MUB
    double vB = va[0]*Bx+va[1]*Ba[1]+va[2]*Ba[2];

    Ua[DD] = R[DD] / (s - va[0]);
    Ua[EN] = (R[EN] + p*va[0] - vB*Bx) / (s-va[0]);
    Ua[MX] = (Ua[EN] + p) * va[0] - vB*Bx;
    Ua[MY] = (Ua[EN] + p) * va[1] - vB*Ba[1];
    Ua[MZ] = (Ua[EN] + p) * va[2] - vB*Ba[2];
    Ua[BX] = Bx;
    Ua[BY] = Ba[1];
    Ua[BZ] = Ba[2];
}

void calc_Uc(double p, double sa, double Bx, double *va, double *vc, 
                double *Bc, double *Ua, double *Uc)
{
    //Eq. (50-52) of MUB

    double vB = vc[0]*Bx + vc[1]*Bc[1] + vc[2]*Bc[2];

    Uc[DD] = Ua[DD] * (sa - va[0]) / (sa - vc[0]);
    Uc[EN] = (sa*Ua[EN] - Ua[MX] + p*vc[0] - vB*Bx) / (sa - vc[0]);
    Uc[MX] = (Uc[EN] + p)*vc[0] - vB*Bx;
    Uc[MY] = (Uc[EN] + p)*vc[1] - vB*Bc[1];
    Uc[MZ] = (Uc[EN] + p)*vc[2] - vB*Bc[2];
    Uc[BX] = Bx;
    Uc[BY] = Bc[1];
    Uc[BA] = Bc[2];
}

