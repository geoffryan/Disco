
#include "paul.h"

int num_diagnostics_r(void);
int num_diagnostics_z(void);
int num_diagnostics_rz(void);
void get_diagnostics( double * x , double * prim , double * Qr , double * Qz,
                        double * Qrz, struct domain * theDomain );

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

