
#include "paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics_r(void){
   //return(5);
   return(0);
}
int num_diagnostics_z(void){
   return(0);
}
int num_diagnostics_rz(void){
   //return(0);
   return(18);
}

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

/* 2D Disk */
/*
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
*/

/* MHD */

void get_diagnostics( double * x , double * prim , double * Qr , double * Qz,
                        double * Qrz, struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];
   double z = x[2];

   double rho = prim[RHO];
   double Pp = prim[PPP];
   double vr = prim[URR];
   double vp = prim[UPP];
   double vz = prim[UZZ];
   double Br = prim[BRR];
   double Bp = prim[BPP];
   double Bz = prim[BZZ];

   double cs2 = gamma_law * Pp / rho;
   double B2 = Br*Br + Bp*Bp + Bz*Bz;

   double sr = rho*vr;
   double sp = r*r*rho*vp;
   double sz = rho*vz;
   double e = 0.5*rho*(vr*vr + r*r*vp*vp + vz*vz) + Pp/(gamma_law-1.0) + 0.5*B2;
   Qrz[0] = rho;
   Qrz[1] = Pp;
   Qrz[2] = vr;
   Qrz[3] = vp;
   Qrz[4] = vz;
   Qrz[5] = Br;
   Qrz[6] = Bp;
   Qrz[7] = Bz;
   Qrz[8] = rho;
   Qrz[9] = e;
   Qrz[10] = sr;
   Qrz[11] = sp;
   Qrz[12] = sz;
   Qrz[13] = rho*vr*vp;
   Qrz[14] = Br*Bp;
   Qrz[15] = B2;
   Qrz[16] = sqrt(cs2);
   Qrz[17] = sqrt(B2/rho);

}

void zero_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int Nqr = theDomain->num_tools_r;
   int Nqz = theDomain->num_tools_z;
   int Nqrz = theDomain->num_tools_rz;
   struct diagnostic_avg * theTools = &(theDomain->theTools);

   int i,k,q;
   for( i=0 ; i<Nr ; ++i ){
      for( q=0 ; q<Nqr ; ++q ){
         int iq = i*Nqr + q;
         theTools->Qr[iq] = 0.0;
      }
   }
   for( k=0 ; k<Nz ; ++k ){
      for( q=0 ; q<Nqz ; ++q ){
         int iq = k*Nqz + q;
         theTools->Qz[iq] = 0.0;
      }
   }
   for( k=0 ; k<Nz ; ++k ){
      for( i=0 ; i<Nr ; ++i ){
         for( q=0 ; q<Nqrz ; ++q ){
            int iq = k*Nr*Nqrz + i*Nqrz + q;
               theTools->Qrz[iq] = 0.0;
         }
      }
   }
   theTools->t_avg = 0.0;

}

void avg_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int Nqr = theDomain->num_tools_r;
   int Nqz = theDomain->num_tools_z;
   int Nqrz = theDomain->num_tools_rz;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   double dt = theTools->t_avg;

   int i,k,q; 
   for( i=0 ; i<Nr ; ++i ){
      for( q=0 ; q<Nqr ; ++q ){
         int iq = i*Nqr + q; 
         theTools->Qr[iq] /= dt; 
      }    
   }
   for(k=0; k<Nz; k++){
      for( q=0 ; q<Nqz ; ++q ){
         int iq = k*Nqz + q; 
         theTools->Qz[iq] /= dt; 
      }
   }
   for(k=0; k<Nz; k++){
      for( i=0 ; i<Nr ; ++i ){
         for( q=0 ; q<Nqrz ; ++q ){
            int iq = k*Nqrz*Nr + i*Nqrz + q; 
            theTools->Qrz[iq] /= dt; 
         }
      }
   }
       
   theTools->t_avg = 0.0; 

}

double get_dV( double * , double * );

void add_diagnostics( struct domain * theDomain , double dt ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Nqr = theDomain->num_tools_r;
   int Nqz = theDomain->num_tools_z;
   int Nqrz = theDomain->num_tools_rz;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   struct cell ** theCells = theDomain->theCells;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double temp_sum_r[Nr*Nqr];
   double temp_sum_z[Nz*Nqz];
   double temp_sum_rz[Nr*Nz*Nqrz];
   double temp_vol_r[Nr];
   double temp_vol_z[Nz];
   double temp_vol_rz[Nr*Nz];
   memset( temp_sum_r , 0 , Nr*Nqr*sizeof(double) );
   memset( temp_sum_z , 0 , Nz*Nqz*sizeof(double) );
   memset( temp_sum_rz , 0 , Nr*Nz*Nqrz*sizeof(double) );
   memset( temp_vol_r , 0 , Nr*sizeof(double) );
   memset( temp_vol_z , 0 , Nz*sizeof(double) );
   memset( temp_vol_rz , 0 , Nr*Nz*sizeof(double) );
   int i,j,k,q;
   int kmin = 0;
   int kmax = Nz;
   int jmin = 0;
   int jmax = Nr;

    for(k=kmin; k<kmax; k++)
    {
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = theCells[jk]+i;
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
                double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
                double xc[3];
                for( q=0 ; q<3 ; ++q ) xc[q] = .5*(xp[q]+xm[q]);
                double dV = get_dV(xp,xm);
                double Qr[Nqr];
                double Qz[Nqz];
                double Qrz[Nqrz];
                get_diagnostics( xc , c->prim , Qr, Qz, Qrz , theDomain );
                for( q=0 ; q<Nqr ; ++q ) 
                    temp_sum_r[ j*Nqr + q ] += Qr[q]*dV;
                for( q=0 ; q<Nqz ; ++q ) 
                    temp_sum_z[ k*Nqz + q ] += Qz[q]*dV;
                for( q=0 ; q<Nqrz ; ++q ) 
                    temp_sum_rz[ jk*Nqrz + q ] += Qrz[q]*dV;
                temp_vol_r[j] += dV;
                temp_vol_z[k] += dV;
                temp_vol_rz[jk] += dV;
            }
        }
    }

    for(j=jmin; j<jmax; j++)
        for(q=0; q<Nqr; q++)
            theTools->Qr[j*Nqr+q] += temp_sum_r[j*Nqr+q]*dt/temp_vol_r[j];
    for(k=kmin; k<kmax; k++)
        for(q=0; q<Nqz; q++)
            theTools->Qz[k*Nqz+q] += temp_sum_z[k*Nqz+q]*dt/temp_vol_z[k];
    for(k=kmin; k<kmax; k++)
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(q=0; q<Nqrz; q++)
            {
                theTools->Qrz[jk*Nqrz+q] += temp_sum_rz[jk*Nqrz+q]*dt
                                                /temp_vol_rz[jk];
            }
        }

   theTools->t_avg += dt;
}

