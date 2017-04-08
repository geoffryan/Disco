
#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics_r(void){
   return(5);
}
int num_diagnostics_z(void){
   return(0);
}
int num_diagnostics_rz(void){
   return(0);
}

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

/* 2D Disk */

void get_diagnostics( double * x , double * prim , double * Qr , double * Qz,
                        double * Qrz, struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];
   double z = x[2];

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];

   Qr[0] = rho;
   Qr[1] = 2.*M_PI*r*rho*vr;
   double Fr,Fp,Fz;
   Fp = 0.0;
   if( theDomain->Npl > 1 ){
      struct planet * pl = theDomain->thePlanets+1;
      planetaryForce( pl , r , phi , 0.0 , &Fr , &Fp , &Fz , 0 );
   }
   Qr[2] = 2.*M_PI*r*rho*(r*Fp);
   Qr[3] = 2.*M_PI*r*rho*( r*r*omega*vr );
   Qr[4] = omega;
}

