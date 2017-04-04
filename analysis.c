
#include "paul.h"
#include "Hydro/metric.h"
#include "Hydro/frame.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics_r(void){
   return(5);
   //return(0);
}
int num_diagnostics_z(void){
   return(0);
}
int num_diagnostics_rz(void){
   return(0);
   //return(18);
   //return(48);
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


/* MHD */
/*
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
*/


/* GRMHD */
/*
void get_diagnostics( double * x , double * prim , double * Qr , double * Qz,
                        double * Qrz, struct domain * theDomain )
{
    double r = x[0];
    double al, be[3], gam[9], igam[9], jac;
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    double sqrtgam = jac / al;

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};
    double B[3] = {prim[BRR]/sqrtgam,prim[BPP]/(r*sqrtgam),prim[BZZ]/sqrtgam};

    double U[4];
    frame_U(x, U);

    double w, u0, u2;
    double u[3];
    double igaml[3];
    int i,j;
    for(i=0; i<3; i++)
    {
        igaml[i] = 0.0;
        for(j=0; j<3; j++)
            igaml[i] += igam[3*i+j]*l[j];
    }
    u2 = l[0]*igaml[0] + l[1]*igaml[1] + l[2]*igaml[2];
    w = sqrt(1.0 + u2);
    u0 = w/al;
    for(i=0; i<3; i++)
        u[i] = igaml[i] - be[i]*u0;

    double l0 = -al*w + be[0]*l[0] + be[1]*l[1] + be[2]*l[2];
    double uU = U[0]*l0 + U[1]*l[0] + U[2]*l[1] + U[3]*l[2];

    double uB, b0, bd[3], B2, b2, UB, Ub;
    uB = l[0]*B[0] + l[1]*B[1] + l[2]*B[2];
    b0 = uB / al;
    bd[0] = (gam[0]*B[0] + gam[1]*B[1] + gam[2]*B[2] + l[0]*uB) / w;
    bd[1] = (gam[3]*B[0] + gam[4]*B[1] + gam[5]*B[2] + l[1]*uB) / w;
    bd[2] = (gam[6]*B[0] + gam[7]*B[1] + gam[8]*B[2] + l[2]*uB) / w;
    double b[3] = {(B[0]+uB*u[0])/w, (B[1]+uB*u[1])/w, (B[2]+uB*u[2])/w};
    B2 = 0.0;
    UB = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {    
            B2 += gam[3*i+j]*B[i]*B[j];
            UB += gam[3*i+j] * B[i] * (U[0]*be[j] + U[j+1]);
        }
    b2 = (B2 + uB*uB) / (w*w);
    Ub = (UB + uU*uB) / w;
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp + b2;
    double rhoe = Pp / (gamma_law-1.0);
    double rhohs = rhoh + b2;
    double Ps = Pp + 0.5*b2;

    double D = jac*rho*u0;
    double tauU = jac * (-(rhoe+0.5*b2)*uU*u0 - (Pp+0.5*b2)*(uU*u0+U[0])
                            - rho*(uU+1.0)*u0 + Ub*b0);
    double tau = sqrtgam * ((rhoe+0.5*b2)*w*w + (Pp+0.5*b2)*u2
                            + rho*w*u2/(w+1) - b0*b0);
    double Sr = jac * (rhohs*u0*l[0] - b0*bd[0]);
    double Sp = jac * (rhohs*u0*l[1] - b0*bd[1]);
    double Sz = jac * (rhohs*u0*l[2] - b0*bd[2]);

    double bd0 = -al*al*b0 + be[0]*bd[0] + be[1]*bd[1] + be[2]*bd[2];
    double uu[4] = {u0, u[0], u[1], u[2]};
    double uud[4] = {l0, l[0], l[1], l[2]};
    double bb[4] = {b0, b[0], b[1], b[2]};
    double bbd[4] = {bd0, bd[0], bd[1], bd[2]};

    double Tgasuu[4][4];
    double Tmaguu[4][4];
    double Tmagbb[4][4];
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
        {
            Tgasuu[i][j] = rhoh*uud[i]*uu[j];
            Tmaguu[i][j] = b2*uud[i]*uu[j];
            Tmagbb[i][j] = -bbd[i]*bb[j];
        }

    // Prims
    Qrz[0] = rho;
    Qrz[1] = Pp;
    Qrz[2] = l[0];
    Qrz[3] = l[1];
    Qrz[4] = l[2];
    Qrz[5] = B[0];
    Qrz[6] = B[1];
    Qrz[7] = B[2];

    // Cons
    Qrz[8] = D;
    Qrz[9] = tauU;
    Qrz[10] = Sr;
    Qrz[11] = Sp;
    Qrz[12] = Sz;

    //Fluxes
    Qrz[13] = jac*rho*uu[1];
    Qrz[14] = jac*rho*uu[2];
    Qrz[15] = jac*rho*uu[3];

    Qrz[16] = jac*Tgasuu[1][1];
    Qrz[17] = jac*Tgasuu[1][2];
    Qrz[18] = jac*Tgasuu[1][3];
    Qrz[19] = jac*Tgasuu[2][1];
    Qrz[20] = jac*Tgasuu[2][2];
    Qrz[21] = jac*Tgasuu[2][3];
    Qrz[22] = jac*Tgasuu[3][1];
    Qrz[23] = jac*Tgasuu[3][2];
    Qrz[24] = jac*Tgasuu[3][3];
    
    Qrz[25] = jac*Tmaguu[1][1];
    Qrz[26] = jac*Tmaguu[1][2];
    Qrz[27] = jac*Tmaguu[1][3];
    Qrz[28] = jac*Tmaguu[2][1];
    Qrz[29] = jac*Tmaguu[2][2];
    Qrz[30] = jac*Tmaguu[2][3];
    Qrz[31] = jac*Tmaguu[3][1];
    Qrz[32] = jac*Tmaguu[3][2];
    Qrz[33] = jac*Tmaguu[3][3];
    
    Qrz[34] = jac*Tmagbb[1][1];
    Qrz[35] = jac*Tmagbb[1][2];
    Qrz[36] = jac*Tmagbb[1][3];
    Qrz[37] = jac*Tmagbb[2][1];
    Qrz[38] = jac*Tmagbb[2][2];
    Qrz[39] = jac*Tmagbb[2][3];
    Qrz[40] = jac*Tmagbb[3][1];
    Qrz[41] = jac*Tmagbb[3][2];
    Qrz[42] = jac*Tmagbb[3][3];

    // Miscellaneous
    Qrz[43] = B2;
    Qrz[44] = b2;
    Qrz[45] = sqrt(gamma_law*Pp / rhoh);
    Qrz[46] = sqrt(b2 / rhohs);
    Qrz[47] = tau;
}
*/

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

