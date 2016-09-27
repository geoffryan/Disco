
#include "../paul.h"
#include "../Calc/calc.h"

static double M = 1.0;
static double gam = 0.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);
    double R0 = 10.0;
    double rs = 10.0;
    double Mdot = 1.0;
    double b0 = 1.0e-4;
    
    double us2 = M / (2*rs);
    double as2 = us2 / (1 - 3*us2);
    double us = -sqrt(us2);

    double rhos = -Mdot / (4*M_PI*rs*rs*us);
    double K = as2 / (gam * pow(rhos,gam-1) * (1-as2/(gam-1)));

    double a02 = (gam-1) * (1 - sqrt(1+3*as2) * fabs(1-as2/(gam-1)));
    double rho0 = pow(as2 / (gam * K * (1-as2/(gam-1))), 1.0/(gam-1));
    double h0 = 1.0 + gam/(gam-1) * K * pow(rho0, gam-1);

    double uRsc = bondi_rel_solve(Mdot, M, R, gam, h0, K, sqrt(as2));
    double rhoB = -Mdot / (4*M_PI*R*R*uRsc);
    double PB = K * pow(rhoB, gam);


    double u0sc = sqrt((1+uRsc*uRsc/(1-2/R))/(1-2/R));

    double uRks = uRsc;
    double u0ks = u0sc + uRsc/(R/2-1);
    double lRks = 2/R*u0ks + (1+2/R)*uRks;

    double chi = R>R0 ? 1.0 : 0.5*(1+tanh(tan(M_PI*(R-0.5*R0)/R0)));

    double rho = rhoB;
    double P = PB;
    double lR = lRks;

    if(R < 2.5)
    {
        rho = rhos;
        P = K * pow(rhos, gam);
        lR = -sqrt(2/R) / (1+sqrt(2/R));
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = r/R * lR * chi;
    prim[UPP] = 0.0;
    prim[UZZ] = z/R * lR * chi;

    if(NUM_C > 5)
    {
        prim[BRR] = r/(R*R*R)*b0 * chi;
        prim[BPP] = 0.0;
        prim[BZZ] = z/(R*R*R)*b0 * chi;
    }

    int q;
    for(q=NUM_C; q < NUM_Q; q++)
    {
        prim[q] = 0.0;
    }
}
