#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define NN       70

#define alpha    -3.0

double  AA[NN + 2];
double  CC[2 * NN + 2 + 1];  /* factorial is absorbed in the coeffs */

#define sqr(x)  ((x) * (x))

#define NPOINTS 400
#define x0      0.0
#define x1      4.0
#define npoints ((double)(NPOINTS))

/* the value of \partial_r\phi_u|_0 */
#define DPHIU_0  0.4
/* the value of \partial_r\phi_d|_0 */
#define DPHID_0  0.4

/* OMEGA: the parameter of masses, see [1], eq.(47) */
#define omega  1.0

/* the winding numbers of the u- and d- quarks */
#define _n  1.0
#define _k  1.0

/* the coupling constant */
#define _g 1.0

/*
 * returns -1^k
 */
double
minus (int k)
{
    return (k % 2 == 0 ? 1.0 : -1.0);
}

double 
epsilon()
{
    if ( _n + omega * _k > 0)
        return 1.0;
    else 
        if ( _n + omega * _k < 0)
            return -1.0;
        else
        {
            fprintf(stderr, "epsilon: n + omega * k = 0\n");
            exit(1);
        }
}

int
main(int argc, char *argv[])
{
    int    i, m;
    FILE  *ff_1, *ff_d, *fphi_u, *fphi_d;
    double delta;
    double f_1 = 1.0;
    double f_d = 1.0;
    double phi_u = 0.0;
    double dphi_u = 0.0;
    double phi_d = 0.0;
    double dphi_d = 0.0;
    double df_1 = 0.0;
    double df_d = 0.0;

    AA[0] = 1.0;
    AA[1] = - 0.25;

    AA[2] = alpha;
#if 0
    AA[2] = alpha * sqr(AA[1]);
#endif

    /* now calculate the profile */
    ff_1 = fopen("profile.f_1", "w+");
    if (ff_1 == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    ff_d = fopen("profile.f_d", "w+");
    if (ff_d == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fphi_u = fopen("profile.phi_u", "w+");
    if (fphi_u == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fphi_d = fopen("profile.phi_d", "w+");
    if (fphi_d == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fprintf(ff_1, "#  r   f_1(r) \n");
    fprintf(ff_d, "#  r   f_d(r) \n");
    fprintf(fphi_u, "#  r   phi_u(r) \n");
    fprintf(fphi_d, "#  r   phi_d(r) \n");
    
    delta = (x1 - x0)/(npoints - 1.0);

    /* the value at x = 0 */
    phi_u = 0.0;
    phi_d = 0.0;
    f_1 = epsilon() * ( _n + _k * 0.5 );
    f_d = epsilon() * _k;
    for (i = 0; i < NPOINTS; i++)
    {
        int    k;
        double ii = (double)(i);
        double x = x0 + delta * ii;

        double res_f = 0.0;
        double res_phi = 0.0;
        double tt;

        double df, dphi;


        fprintf(ff_1, " %g %g\n", x, f_1);
        fprintf(ff_d, " %g %g\n", x, f_d);
        fprintf(fphi_u, " %g %g\n", x, phi_u);
        fprintf(fphi_d, " %g %g\n", x, phi_d);

        if (fabs(x) < 1.00e-5)
        {
            dphi_u = DPHIU_0;
            dphi_d = DPHID_0;
            printf("initial value of dphi_u is %g, x = %g\n", dphi_u, x);
            printf("initial value of dphi_d is %g, x = %g\n", dphi_d, x);
        }
        else
        {
            dphi_u = (1.0 / x) * ( f_1 * phi_u - 0.5 * f_d * phi_u );
            dphi_d = (1.0 / x) * f_d * phi_d;
        }

        df_1 = (1./4.) * (_g*_g) * x * (phi_u * phi_u - 1.0);
        df_d = (1./3.) * (_g*_g) * x * (phi_d * phi_d - 0.5 * phi_u * phi_u +
                                        0.5 - omega);

        /* increment the functions */
        f_1 += df_1 * delta;
        f_d += df_d * delta;
        phi_u += dphi_u * delta;
        phi_d += dphi_d * delta;
    }
    printf("f_1 = %g, df_1 = %g, f_d = %g, df_d = %g\n",
           f_1, df_1, f_d, df_d);
    printf("phi_u = %g, dphi_u = %g, phi_d = %g, dphi_d = %g\n",
           phi_u, dphi_u, phi_d, dphi_d);

    fclose(ff_1);
    fclose(ff_d);
    fclose(fphi_u);
    fclose(fphi_d);
}

