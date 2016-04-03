#include "../../paul.h"
#include "../metric.h"

/*
 * This is an approximate metric for binary black holes, constructed as the
 * superposition of two kerr-schild metrics.  It is the initial condition
 * for some iterative schemes that fully solve the constraint equations to
 * generate initial data for numerical relativity codes.
 * See: Marronnetti et al. 2000 (gr-qc/0001077v2)
 *      Yo et al. 2004 (PRD.10.084033)
 *
 * An improvement would be to incorporate:
 * Tichy et al. 2002 (gr-qc/0207011) and/or
 * Nissanke 2006 (gr-qc/0509128)
 *
 * And finally to use:
 * Ireland et al. 2016 (gr-qc/1512.05650)
 */

static double M = 1.0; 
static double q = 1.0;
static double a = 100.0;

void setMetricParams(struct domain *theDomain)
{
   //om = theDomain->theParList.MetricPar1;
   //M = theDomain->theParList.MetricPar2;
}

double metric_lapse(double x[3])
{
    double x = x[0]*cos(x[1]);
    double y = x[0]*sin(x[1]);
    double z = x[2];

    double M1 = M/(1+q);
    double M2 = q*M/(1+q);
    double a1 = q*a/(1+q);
    double a2 = -a/(1+q);

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double al1 = 1.0/sqrt(1.0 + 2*M1/R1);
    double al2 = 1.0/sqrt(1.0 + 2*M2/R2);

    return al1*al2;
}

void metric_shift(double x[3], double b[3])
{
    double r = x[0];
    double cosp = cos(x[1]);
    double sinp = sin(x[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = x[2];

    double M1 = M/(1+q);
    double M2 = q*M/(1+q);
    double a1 = q*a/(1+q);
    double a2 = -a/(1+q);

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double k1x = (a1-x)/R1;
    double k1y = -y/R1;
    double k1z = -z/R1;
    double k2x = (a2-x)/R2;
    double k2y = -y/R2;
    double k2z = -z/R2;

    double be1x = -2*H1*k1x/(1+2*H1);
    double be1y = -2*H1*k1y/(1+2*H1);
    double be1z = -2*H1*k1z/(1+2*H1);
    double be2x = -2*H2*k2x/(1+2*H2);
    double be2y = -2*H2*k2y/(1+2*H2);
    double be2z = -2*H2*k2z/(1+2*H2);
    
    double om = sqrt(M/(a*a*a));

    b[0] = cosp*(be1x+be2x) + sinp*(be1y+be2y);
    b[1] = (-sinp*(be1x+be2x) + cosp*(be1y+be2y))/r + om;
    b[2] = be1z + be2z;
}

void metric_gam(double x[3], double gam[9])
{
    double r = x[0];
    double cosp = cos(x[1]);
    double sinp = sin(x[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = x[2];

    double M1 = M/(1+q);
    double M2 = q*M/(1+q);
    double a1 = q*a/(1+q);
    double a2 = -a/(1+q);

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    //Covariant (w.r.t. Minkowski) k components in r/phi
    double k1r = (cp*a1 - r) / R1;
    double k1p = r*(-sp*a1) / R1;
    double k1z = -z/R1;
    double k2r = (cp*a2 - r) / R2;
    double k2p = r*(-sp*a2) / R2;
    double k2z = -z/R2;

    gam[0] = 1.0 + 2*H1*k1r*k1r + 2*H2*k2r*k2r;
    gam[1] = 2*H1*k1r*k1p + 2*H2*k2r*k2p;
    gam[2] = 2*H1*k1r*k1z + 2*H2*k2r*k2z;
    gam[3] = gam[1];
    gam[4] = r*r + 2*H1*k1p*k1p + 2*H2*k2p*k2p;
    gam[5] = 2*H1*k1p*k1z + 2*H2*k2p*k2z;
    gam[6] = gam[2];
    gam[7] = gam[5];
    gam[8] = 1.0 + 2*H1*k1z*k1z + 2*H2*k2z*k2z;
}

void metric_igam(double x[3], double igam[9])
{
    double r = x[0];
    double cosp = cos(x[1]);
    double sinp = sin(x[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = x[2];

    double M1 = M/(1+q);
    double M2 = q*M/(1+q);
    double a1 = q*a/(1+q);
    double a2 = -a/(1+q);

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR12 = 1.0/(R1*R1);
    double iR22 = 1.0/(R2*R2);

    double Px2 = y*y+z*z;

    double idetgam = 1.0/(1.0 + 2.0*H1 + 2.0*H2
                            + 4.0*H1*H2*a*a*Px2*iR12*iR22);

    double igxx = idetgam*(1.0 + 2*H1*Px2*iR12 + 2*H2*Px2*iR22);
    double igxy = idetgam*(2*H1*(a1-x)*y*iR12 + 2*H2*(a2-x)*y*iR22);
    double igxy = idetgam*(2*H1*(a1-x)*z*iR12 + 2*H2*(a2-x)*z*iR22);
    double igyy = idetgam*(1 + 2*H1*(R1*R1-y*y)*iR12 + 2*H2*(R2*R2-y*y)*iR22
                            + 4*H1*H2*a*a*z*z*iR12*iR22);
    double igyz = idetgam*(-2*H1*y*z*iR12 - 2*H2*y*z*iR22
                            - 4*H1*H2*a*a*y*z*iR12*iR22);
    double igzz = idetgam*(1 + 2*H1*(R1*R1-z*z)*iR12 + 2*H2*(R2*R2-z*z)*iR22
                            + 4*H1*H2*a*a*y*y*iR12*iR22);

    
    igam[0] = cosp*cosp*igxx + 2*cosp*sinp*igxy + sinp*sinp*igyy;
    igam[1] = (cosp*sinp*(-igxx+igyy) + (cosp*cosp-sinp*sinp)*igxy) / r;
    igam[2] = cosp*igxz + sinp*igyz;
    igam[3] = igam[1];
    igam[4] = (sinp*sinp*igxx - 2*cosp*sinp*igxy + sinp*sinp*igyy) / (r*r);
    igam[5] = -sinp*igxz + cosp*igyz;
    igam[6] = igam[2];
    igam[7] = igam[5];
    igam[8] = igzz;
}

double metric_jacobian(double x[3])
{
    double r = x[0];
    double cosp = cos(x[1]);
    double sinp = sin(x[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = x[2];

    double M1 = M/(1+q);
    double M2 = q*M/(1+q);
    double a1 = q*a/(1+q);
    double a2 = -a/(1+q);

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR12 = 1.0/(R1*R1);
    double iR22 = 1.0/(R2*R2);

    double Px2 = y*y+z*z;

    double detgam = 1.0 + 2.0*H1 + 2.0*H2 + 4.0*H1*H2*a*a*Px2*iR12*iR22;
    double al1 = 1.0/sqrt(1.0 + 2*M1/R1);
    double al2 = 1.0/sqrt(1.0 + 2*M2/R2);
    double al = al1*al2;

    return al * detgam * r;
}

void metric_der_g(double x[3], int i, double dg[16])
{
    int mu;
    for(mu=0; mu<16; mu++)
        dg[mu] = 0.0;

    if(i == 0)
        return;

    double r = x[0];
    double cp = cos(x[1]);
    double sp = sin(x[1]);

    double x = r * cp;
    double y = r * sp;
    double z = x[2];

    double M1 = M/(1+q);
    double M2 = q*M/(1+q);
    double a1 = q*a/(1+q);
    double a2 = -a/(1+q);

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR12 = 1.0/(R1*R1);
    double iR22 = 1.0/(R2*R2);

    //Covariant (w.r.t. Minkowski) k components in r/phi
    double k1r = (cp*a1 - r) / R1;
    double k1p = r*(-sp*a1) / R1;
    double k1z = -z/R1;
    double k2r = (cp*a2 - r) / R2;
    double k2p = r*(-sp*a2) / R2;
    double k2z = -z/R2;

    double al1 = 1.0/sqrt(1.0 + 2*M1/R1);
    double al2 = 1.0/sqrt(1.0 + 2*M2/R2);

    double dal1, dal2, dr2, dH1, dH2, dk1r, dk1p, dk1z, dk2r, dk2p, dk2z;

    if(i == 1) // d/dr
    {
        dr2 = 2*r;
        double dR1 = (r - cp*a1) / R1;
        double dR2 = (r - cp*a2) / R2;

        dal1 = -0.5 * al1*al1*al1 * (-2*M1*iR12*dR1);
        dal2 = -0.5 * al2*al2*al2 * (-2*M2*iR22*dR2);

        dH1 = -M1 * iR12 * dR1;
        dH2 = -M2 * iR22 * dR2;

        dk1r = (-R1-(cp*a1-r)*dR1) * iR12;
        dk2r = (-R2-(cp*a2-r)*dR2) * iR22;
        dk1p = -sp*a1 * (R1 - r*dR12) * iR12;
        dk2p = -sp*a2 * (R2 - r*dR22) * iR22;
        dk1z = z * iR12 * dR1;
        dk2z = z * iR22 * dR2;
    }
    else if(i == 2) // d/dphi
    {
        dr2 = 0.0;
        double dR1 = r*a1*sp / R1;
        double dR2 = r*a2*sp / R2;

        dal1 = -0.5 * al1*al1*al1 * (-2*M1*iR12*dR1);
        dal2 = -0.5 * al2*al2*al2 * (-2*M2*iR22*dR2);

        dH1 = -M1 * iR12 * dR1;
        dH2 = -M2 * iR22 * dR2;

        dk1r = (-sp*a1*R1-(cp*a1-r)*dR1) * iR12;
        dk2r = (-sp*a2*R2-(cp*a2-r)*dR2) * iR22;
        dk1p = -r*a1 * (cp*R1 - sp*dR12) * iR12;
        dk2p = -r*a2 * (cp*R2 - sp*dR22) * iR22;
        dk1z = z * iR12 * dR1;
        dk2z = z * iR22 * dR2;
    }
    else if(i == 3) // d/dz
    {
        dr2 = 0.0;
        double dR1 = z / R1;
        double dR2 = z / R2;

        dal1 = -0.5 * al1*al1*al1 * (-2*M1*iR12*dR1);
        dal2 = -0.5 * al2*al2*al2 * (-2*M2*iR22*dR2);

        dH1 = -M1 * iR12 * dR1;
        dH2 = -M2 * iR22 * dR2;

        dk1r = -(cp*a1-r) * iR12 * dR1;
        dk2r = -(cp*a2-r) * iR22 * dR2;
        dk1p = r*sp*a1 * iR12 * dR1;
        dk2p = r*sp*a2 * iR22 * dR2;
        dk1z = -(R1 - z*dR1) * iR12;
        dk2z = -(R2 - z*dR2) * iR22;
    }

    /*
    gam[0] = 1.0 + 2*H1*k1r*k1r + 2*H2*k2r*k2r;
    gam[1] = 2*H1*k1r*k1p + 2*H2*k2r*k2p;
    gam[2] = 2*H1*k1r*k1z + 2*H2*k2r*k2z;
    gam[3] = gam[1];
    gam[4] = r*r + 2*H1*k1p*k1p + 2*H2*k2p*k2p;
    gam[5] = 2*H1*k1p*k1z + 2*H2*k2p*k2z;
    gam[6] = gam[2];
    gam[7] = gam[5];
    gam[8] = 1.0 + 2*H1*k1z*k1z + 2*H2*k2z*k2z;
    */

    double dgrr, dgrp, dgrz, dgpp, dgpz, dgzz;

    dgrr = 2*dH1*k1r*k1r + 4*H1*k1r*dk1r + 2*dH2*k2r*k2r + 4*H2*k2r*dk2r;
    dgrp = 2*dH1*k1r*k1p + 2*H1*dk1r*k1p + 2*H1*k1r*dk1p
             + 2*dH2*k2r*k2p + 2*H2*dk2r*k2p + 2*H2*k2r*dk2p;
    dgrz = 2*dH1*k1r*k1z + 2*H1*dk1r*k1z + 2*H1*k1r*dk1z
             + 2*dH2*k2r*k2z + 2*H2*dk2r*k2z + 2*H2*k2r*dk2z;
    dgpp = d2r + 2*dH1*k1p*k1p + 4*H1*k1p*dk1p + 2*dH2*k2p*k2p
                + 4*H2*k2p*dk2p;
    dgpz = 2*dH1*k1p*k1z + 2*H1*dk1p*k1z + 2*H1*k1p*dk1z
                + 2*dH2*k2p*k2z + 2*H2*dk2p*k2z + 2*H2*k2p*dk2z;
    dgzz = 2*dH1*k1z*k1z + 4*H1*k1z*dk1z + 2*dH2*k2z*k2z + 4*H2*k2z*dk2z;
    dg[0] = ;
    dg[1] = ;
    dg[2] = ;
    dg[3] = ;
    dg[4] = dg[1];
    dg[5] = dgrr;
    dg[6] = dgrp;
    dg[7] = dgrz;
    dg[8] = dg[2];
    dg[9] = dgrp;
    dg[10] = dgpp;
    dg[11] = dgpz;
    dg[12] = dg[3];
    dg[13] = dgrz;
    dg[14] = dgpz;
    dg[15] = dgzz;
}

void metric_der_lapse(double x[3], double da[4])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double st = r/R;
    double ct = z/R;

    double MoR2 = M/(R*R);
    double a = 1.0/sqrt(1.0+2*M/R);

    da[0] = 0.0;
    da[1] = a*a*a * MoR2 * st;
    da[2] = 0.0;
    da[3] = a*a*a * MoR2 * ct;
}

void metric_der_shift(double x[3], double db[12])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double st = r/R;
    double ct = z/R;

    double denom = 1.0/(R*R + 2*M*R);
    double d_denom = -denom*denom * 2*(R+M);

    int i;

    //dt
    for(i=0; i<3; i++)
        db[i] = 0;

    //dr
    db[3] = 2*M * denom + 2*M*r * d_denom*st;
    db[4] = 0.0;
    db[5] = 2*M*z * d_denom*st;

    //dphi
    for(i=0; i<3; i++)
        db[6+i] = 0;

    //dz
    db[9] = 2*M*r * d_denom*ct;
    db[10] = 0.0;
    db[11] = 2*M * denom + 2*M*z * d_denom*ct;
}

int metric_killing(int mu)
{
    if(mu == 0)
        return 1;
    return 0;
}
