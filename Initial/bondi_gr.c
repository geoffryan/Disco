
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

    double P0 = 1.0e-2;
    double b0 = 1.0e-4;

    double r = x[0];
    double z = x[1];
    double R = sqrt(r*r+z*z);
    double R0 = 10.0;
    double chi = R>R0 ? 1.0 : 0.5*(1+tanh(tan(M_PI*(R-0.5*R0)/R0)));
    double uRsc = -0.11;
    double u0sc = sqrt((1+uRsc*uRsc/(1-2/R))/(1-2/R));

    double uRks = uRsc;
    double u0ks = u0sc + uRsc/(R/2-1);
    double lRks = 2/R*u0ks + (1+2/R)*uRks;

    if(R < 2.5)
        lRks = -sqrt(2/R) / (1+sqrt(2/R));

    prim[RHO] = 1.5e-3;
    prim[PPP] = 5.5e-5;
    prim[URR] = r/R * lRks * chi;
    prim[UPP] = 0.0;
    prim[UZZ] = z/R * lRks * chi;

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
