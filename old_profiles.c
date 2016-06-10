#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define NN       70

#define alpha    0.5

double  AA[NN + 2];
#define sqr(x)  ((x) * (x))

#define NPOINTS 100
#define x0      0.0
#define x1      2.7
#define npoints ((double)(NPOINTS))

int
main(int argc, char *argv[])
{
    int    i;
    FILE  *fa, *fphi;
    double delta;

    double res_f, res_phi;    /* residual terms of the series */

    AA[0] = 1.0;
    AA[1] = - 0.25;
    AA[2] = alpha * sqr(AA[1]);

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

        for (k = 0; k <= NN + 1; k++)
        {
            double ss = (double)k;
            double power = pow(x, 2.0 * ss);

            if (k == NN + 1)    /* the residual term */
            {
                res_f = AA[k] * power;
                res_phi = (4.0 * ss + 4.0) * AA[k + 1] * power;
            }
            else
            {
                f   += AA[k] * power;
                if (k > 0)
                    phi += (4.0 * ss + 2.0) * AA[k + 1] * power;
            }
        }
        phi = sqrt(phi);

        fprintf(fa, " %g %g %g\n", x, f, res_f);
        fprintf(fphi, " %g %g %g\n", x, phi, res_phi);
    }

    fclose(fa);
    fclose(fphi);
}

