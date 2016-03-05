
#include "../paul.h"

// Cylindrical relativistic blast wave in a uniform magnetic field,
// cf. Del Zanna et al 2007

static double gamma_law = 0.0;

void setICparams( struct domain * theDomain ){
    gamma_law = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double r0 = 1.0;
   double B0 = 0.1;

   double relfac = 1.0;
   
   double rhoL = 1.0e-2;
   double PL = relfac*relfac*1.0;
   double rhoR = 1.0e-4;
   double PR = relfac*relfac*3.0e-5;

   double rho = r > r0 ? rhoR : rhoL;
   double Pp  = r > r0 ? PR : PL;


   prim[RHO] = rho; 
   prim[PPP] = Pp;
   prim[URR] = 0.0; 
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0; 

   if(NUM_C == 8)
   {
       prim[BRR] =  relfac * B0 * cos(phi);
       prim[BPP] = -relfac * B0 * sin(phi);
       prim[BZZ] = 0.0; 
   }

   double X = cos(phi)>0 ? 1.0 : 0.0;
   if( NUM_N>0) prim[NUM_C] = X;

}
