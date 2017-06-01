
#include "../paul.h"

static double nu = 0.0;
static double r0 = 0.0;
static double dr = 0.0;
static double sig0 = 0.0;
static double sig1 = 0.0;
static double mach = 0.0;
static double gam = 0.0;

void setICparams( struct domain * theDomain ){
   nu = theDomain->theParList.viscosity;
   r0 = theDomain->theParList.initPar1;
   dr = theDomain->theParList.initPar2;
   sig0 = theDomain->theParList.initPar3;
   sig1 = theDomain->theParList.initPar4;
   mach = theDomain->theParList.Disk_Mach;
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double rho = sig0 + (sig1-sig0)*exp(0.5*(r-r0)*(r-r0)/(dr*dr));
   double cs20 = 1.0/r0 / (mach*mach);
   double Pp = sig1 * cs20 / gam;
   double omega = 1.0/pow(r,1.5);
   double vr = -1.5*nu/rho/r;

   double X = 0.0; 
   if( cos(phi) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
