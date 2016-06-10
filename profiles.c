#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define NN       70

#define alpha    0.033

double  AA[NN + 2];
double  CC[2 * NN + 2 + 1];  /* factorial is absorbed in the coeffs */

#define sqr(x)  ((x) * (x))

#define NPOINTS 500
#define x0      0.0
#define x1      2.0
#define npoints ((double)(NPOINTS))


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

    AA[0] = 1.0;
    AA[1] = - 0.25;

    AA[2] = alpha;
#if 0
    AA[2] = alpha * sqr(AA[1]);
#endif

    printf("fact(5) = %g\n", fact(5));
    printf("fact(6) = %g\n", fact(6));
    printf("minus(10) = %g\n", minus(10));
    printf("minus(11) = %g\n", minus(11));
    printf("minus(12) = %g\n", minus(12));

    /* resolve the coefficients */
    for (i = 3; i <= NN + 1; i++)
    {
        int     k;
        double  tt = 0.0;
        double  ss = (double)(i);

        for (k = 2; k <= i - 1; k++)
        {
            double pp = (double)(k);

            tt += pp * AA[k] * AA[i - k];
        }
        tt /= ss * (ss - 2.0);

        AA[i] = tt;
    }

    /* resolve the CC coefficients */
    CC[0] = -1.0;
    CC[1] = -1.0;
    for (m = 2; m <= 2 * NN + 2 + 1; m++)
    {
        int k;

        double tt = 0;

        for (k = 0; k <= m - 1; k++)
        {
            double kk = (double)(k);

            tt += CC[k] * minus(m - k + 1) / fact(m - k);
        }
        if (m % 2 == 0)
        {
            tt += AA[m / 2];
        }
        CC[m] = tt;
        printf("got CC[%d] = %g\n", m, CC[m]);
    }
    printf("CC[2] = %g, CC[3] = %g\n", CC[2], CC[3]);

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
    fprintf(fa, "#  r   f(r) <residual term>\n");
    fprintf(fphi, "#  r   phi(r) <residual term>\n");
    
    delta = (x1 - x0)/(npoints - 1.0);
    for (i = 0; i < NPOINTS; i++)
    {
        int    k;
        double ii = (double)(i);
        double x = x0 + delta * ii;

        double f = 0.0;
        double phi = 0.0;
        double res_f = 0.0;
        double res_phi = 0.0;
        double tt;

        for (k = 0; k <= 2 * NN + 2 + 1; k++)
        {
            double ss = (double)k;
            double power = pow(x, ss);

            if (k == 2 * NN + 2 + 1)    /* the residual term */
            {
                res_f = CC[k] * power;
#if 0
                res_phi = (4.0 * ss + 4.0) * AA[k + 1] * power;
#endif
            }
            else
            {
                f   += CC[k] * power;
#if 0
                if (k > 0)
                    phi += (4.0 * ss + 2.0) * AA[k + 1] * power;
#endif
            }
        }
        tt = exp ( - x );   
/* printf( " exp (-x ) = %g\n", tt); */
        res_f *= tt;
        f *= tt;
        f += 2.0;

#if 0
        phi = sqrt(phi);
#endif

        fprintf(fa, " %g %g %g\n", x, f, res_f);
#if 0
        fprintf(fphi, " %g %g %g\n", x, phi, res_phi);
#endif
    }

    fclose(fa);
    fclose(fphi);
}

