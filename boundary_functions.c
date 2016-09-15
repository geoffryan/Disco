
#include "paul.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
void cons2prim( double * , double * , double * , double );
void prim2cons( double * , double * , double * , double );
void subtract_omega( double * );
void reflect_prims(double *, double *, int);

void set_cell_init(struct cell *c, double *r_jph, double *z_kph, int j, int k)
{
    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
    double r = get_moment_arm(xp, xm);
    double phi = c->piph - 0.5*c->dphi;
    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
    initial(c->prim, x);
    subtract_omega(c->prim);
}

void set_cells_copy(struct cell *c, int Np, struct face *theFaces, 
                    int n0, int n1, int LR)
{
    int i,n,q;

    // Clear annulus.
    for(i=0; i<Np; i++)
    {
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] = 0.0;
        c[i].tempDoub = 0.0;
    }

    // Add outer strip.
    for(n=n0; n<n1; n++)
    {
        struct face *f = theFaces+n;
        struct cell *cDst, *cSrc;
        if(LR > 0)
        {
            cSrc = f->L;
            cDst = f->R;
        }
        else
        {
            cSrc = f->R;
            cDst = f->L;
        }
        for(q=0; q<NUM_Q; q++)
            cDst->prim[q] += cSrc->prim[q] * f->dA;
        cDst->tempDoub += f->dA;
    }
    
    // Divide by total area.
    for(i=0; i<Np; i++)
    {
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] /= c[i].tempDoub;
    }
}

void set_cells_copy_distant(struct cell *c, int Np, struct cell *c1, int Np1)
{
    // Copies annulus c1 to c. Then calls reflect_prim() on each cell in c.

    int i, q;

    // Clear annulus.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] = 0.0;

    // Find first cell in c1 that intersects with c[0].
    int i1 = 0;
    double phi0 = c[0].piph - c[0].dphi;
    for(i1=0; i1 < Np1; i1++)
    {
        double dphip = c1[i1].piph - phi0;
        while(dphip > M_PI)
            dphip -= 2*M_PI;
        while(dphip < -M_PI)
            dphip += 2*M_PI;
        double dphim = dphip - c1[i1].dphi;

        if(dphip > 0 && dphim <= 0.0)
            break;
    }

    //Loop through all cells c, add contributions from all neighbours in c1.
    for(i=0; i<Np; i++)
    {
        double phip = c[i].piph;
        double phim = c[i].piph - c[i].dphi;
        double phip1, phim1;
        do
        {
            phip1 = c1[i1].piph;
            phim1 = c1[i1].piph - c1[i1].dphi;
            while(phip1 < phim)
            {
                phip1 += 2*M_PI;
                phim1 += 2*M_PI;
            }
            while(phim1 > phip)
            {
                phip1 -= 2*M_PI;
                phim1 -= 2*M_PI;
            }
            double phi1 = phip < phip1? phip : phip1;
            double phi2 = phim > phim1? phim : phim1;
            double dphi = phi1-phi2;
            for(q=0; q<NUM_Q; q++)
                c[i].prim[q] += c1[i1].prim[q] * dphi;
            i1++;
        }
        while(phip1 < phip);
    }

    // Divide by total area.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] /= c[i].dphi;
}

void boundary_fixed_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(j=0; j<Ng; j++)
            for(k=0; k<Nz; k++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(j=Nr-Ng; j<Nr; j++)
            for(k=0; k<Nz; k++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(j=0; j<Nr; j++)
            for(k=0; k<Ng; k++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}
void boundary_fixed_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(j=0; j<Nr; j++)
            for(k=Nz-Ng; k<Nz; k++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_zerograd_rinn( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_1;
    int *fIndex = theDomain->fIndex_r;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k,n,q;

    if(dim_rank[0] == 0 )
    {
        for(j=Ng-1; j>=0; j--)
            for(k=0; k<Nz; k++)
            {
                int jk = j+Nr*k;
                int n0 = fIndex[jk];
                int n1 = fIndex[jk+1];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, -1);
            }
    }
}

void boundary_zerograd_rout( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_1;
    int *fIndex = theDomain->fIndex_r;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k,n,q;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(j=Nr-Ng; j<Nr; j++)
            for(k=0; k<Nz; k++)
            {
                int jk = j+Nr*k;
                int n0 = fIndex[jk-1];
                int n1 = fIndex[jk];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, +1);
            }
    }
}

void boundary_zerograd_zbot( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_2;
    int *fIndex = theDomain->fIndex_z;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k,n,q;

    if(dim_rank[1] == 0)
    {
        for(j=0; j<Nr; j++)
            for(k=Ng-1; k>=0; k--)
            {
                int jk = j+Nr*k;
                int n0 = fIndex[jk];
                int n1 = fIndex[jk+1];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, -1);
            }
    }
}

void boundary_zerograd_ztop( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_2;
    int *fIndex = theDomain->fIndex_z;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k,n,q;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(j=0; j<Nr; j++)
            for(k=Nz-Ng; k<Nz; k++)
            {
                int jk = j+Nr*k;
                int n0 = fIndex[jk-Nr];
                int n1 = fIndex[jk-Nr+1];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, +1);
            }
    }
}

void boundary_reflect_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(j=Ng-1; j>=0; j--)
            for(k=0; k<Nz; k++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*Ng-j-1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                }
            }
    }
}

void boundary_reflect_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(j=Nr-Ng; j<Nr; j++)
            for(k=0; k<Nz; k++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*(Nr-Ng) - j - 1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                }
            }
    }
}

void boundary_reflect_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(j=0; j<Nr; j--)
            for(k=Ng-1; k>=0; k--)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*Ng-k-1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
    }
}

void boundary_reflect_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(j=0; j<Nr; j++)
            for(k=Nz-Ng; k<Nz; k++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*(Nz-Ng) - k - 1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
    }
}