#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define NN       70

#define alpha    -3.0

double  AA[NN + 2];
double  CC[2 * NN + 2 + 1];  /* factorial is absorbed in the coeffs */

#define sqr(x)  ((x) * (x))

#define NPOINTS 2000
#define x0      0.0
#define x1      20.0
#define npoints ((double)(NPOINTS))

/* the value of \partial_r\phi|_0 */
#define DPHI_0  0.5983392383086577


double
fact (int k)
{
    int    i;
    double kk, ll;

    if (k < 0)
    {
        fprintf(stderr, "factorial of a negative %d\n", k);
        exit(2);
    }
    if (k == 0 || k == 1)
        return 1.0;

    ll = 1.0;
    for (i = 2; i <= k; i++)
        ll *= (double)(i);

    return ll;
}

/*
 * returns -1^k
 */
double
minus (int k)
{
    return (k % 2 == 0 ? 1.0 : -1.0);
}

int
main(int argc, char *argv[])
{
    int    i, m;
    FILE  *fa, *fphi;
    double delta;
    double f = 1.0;
    double phi = 0.0;
    double dphi = 0.0;
    double df = 0.0;

    AA[0] = 1.0;
    AA[1] = - 0.25;

    AA[2] = alpha;
#if 0
    AA[2] = alpha * sqr(AA[1]);
#endif

    /* now calculate the profile */
    fa = fopen("profile.a", "w+");
    if (fa == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fphi = fopen("profile.phi", "w+");
    if (fphi == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fprintf(fa, "#  r   f(r) \n");
    fprintf(fphi, "#  r   phi(r) \n");
    
    delta = (x1 - x0)/(npoints - 1.0);

    /* the value at x = 0 */
    phi = 0.0;
    f = 1.0;
    for (i = 0; i < NPOINTS; i++)
    {
        int    k;
        double ii = (double)(i);
        double x = x0 + delta * ii;

        double res_f = 0.0;
        double res_phi = 0.0;
        double tt;

        double df, dphi;


        fprintf(fa, " %g %g\n", x, f);

        fprintf(fphi, " %g %g\n", x, phi);

        if (fabs(x) < 1.00e-5)
        {
            dphi = DPHI_0;
            printf("initial value of dphi is %g, x = %g\n", dphi, x);
        }
        else
            dphi = (1.0 / x) * f * phi;

        df = 0.5 * x * (phi * phi - 1.0);


        /* increment the functions */
        f += df * delta;
        phi += dphi * delta;
    }
    printf("f = %g, df = %g, phi = %g, dphi = %g\n",
           f, df, phi, dphi);

    fclose(fa);
    fclose(fphi);
}

