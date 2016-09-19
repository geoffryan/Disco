
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

    double P0 = 1.0e-2;
    double b0 = 1.0e-2;

    prim[RHO] = 1.0;
    prim[PPP] = 1.0e-2;
    prim[URR] = 0.0;
    prim[UPP] = 0.0;
    prim[UZZ] = 0.0;

    if(NUM_C > 5)
    {
        prim[BRR] = b0;
        prim[BPP] = 0.0;
        prim[BZZ] = 0.0;
    }

    int q;
    for(q=NUM_C; q < NUM_Q; q++)
    {
        prim[q] = 0.0;
    }
}
