//http://www.aanda.org/articles/aa/pdf/2010/08/aa12443-09.pdf
#include "../paul.h"

static double Mach = 0.0;
static double gam = 0.0;
static double M = 0.0;
static double r0 = 1.0;
static double r1 = 2.0;
static double r2 = 4.0;
static double H0 = 0.2;

void setICparams( struct domain * theDomain ){
   Mach = theDomain->theParList.Disk_Mach;
   gam = theDomain->theParList.Adiabatic_Index;
   M = theDomain->theParList.metricPar2;
}

double I(double r){
   double v=0.0;
   if( r>1.2*r0 && r<3.8/4.0 * r2 ) v = 1.0;
   return(v);
}

double floor( double x ){
   return( (double)(int)x );
}

double Rs( double x ){
   return( 1.2*r0 + (floor((x-r0)/H0)-.5)*H0 );
}

void initial( double * prim , double * x ){

   double rnd1 = (double)rand()/(double)(RAND_MAX);
   double rnd2 = (double)rand()/(double)(RAND_MAX);

   double r   = x[0];
   double phi = x[1];
   double z   = x[2];

   double rho = 100.0;
   double omega = sqrt(M/(r*r*r));
   double cs2 = M/(2*r1*r1*r1) * H0*H0; // actually P / rho*h
   double Pp = (gam-1.0) * cs2 / (gam*(1-cs2)-1) * rho

   double dvp = (2e-2*rnd1 - 1e-2) * omega;
   double vz = (2e-2*rnd2 - 1e-2) * sqrt(M/r0);
   double vr = 0.0;
   double vp = omega +dvp;

   double A0 = 0.37;
   double Bz = A0*sin(2.*M_PI*(r-r0)/H0)/r*I(r)/sqrt(Rs(r));

   double R2 = r*r+z*z;
   double R = sqrt(R2);
   double R3 = R*R2;

   double g00 = -1+2*M/R;
   double g0r = 2*M*r/R2;
   double g0z = 2*M*z/R2;
   double grr = 1+2*M*r*r/R3;
   double grz = 2*M*r*z/R3;
   double gzz = 1+2*M*z*z/R3;
   double gpp = r*r;
   double u0 = 1.0 / sqrt(-g00 - 2*vr*g0r - 2*vz*g0z - grr*vr*vr - 2*grz*vr*vz
                            - gzz*vz*vz - gpp*vp*vp);
   double lr = u0*(g0r+grr*vr+grz*vz); 
   double lp = u0*gpp*vp; 
   double lz = u0*(g0z+grz*vr+gzz*vz); 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = lr;
   prim[UPP] = lp;
   prim[UZZ] = lz;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

