#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>


/* 
   The method of the variable 1/r.
   We change to the variable t = 1/r and write the equations
   with respect to this variable. And then solve them
 */

#define NN       70

#define alpha    -3.0

double  AA[NN + 2];
double  CC[2 * NN + 2 + 1];  /* factorial is absorbed in the coeffs */

#define sqr(x)  ((x) * (x))

#define NPOINTS 2000
#define t0      0.1
#define npoints ((double)(NPOINTS))

/* the border that approaches real zero */
#define t2      100.0


#define FLUX    3.0
#define t1      (1.0 / 0.08)
#define PHI_INIT 0.00003906205

#if 0
#define FLUX   2.0
#define t1      20.0
#define PHI_INIT 0.00060533
#endif

#if 0
/* These values are good for n = 1 */
#define t1      20.0
#define F_INIT   (1.0 - 0.01)
#define PHI_INIT 0.03137506
#endif 

double f_init = 0.0;
double phi_init = 0.0;


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


    /* 
     * first we calculate the initial values
     * using very basic asymptotics:
     *   f  = n  - 1/4 * x^2
     */
    f_init   = FLUX - 0.25 * (1./t1) * (1/t1);
    phi_init = PHI_INIT;

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
    

    /* the value at t = t1 */
    delta = (t1 - t0)/(npoints - 1.0);
    phi = phi_init;
    f = f_init;
    for (i = NPOINTS - 1; i > 0; i--)
    {
        int    k;
        double ii = (double)(i);
        double t = t0 + delta * ii;

        double res_f = 0.0;
        double res_phi = 0.0;
        double tt;

        double df, dphi;


        fprintf(fa, " %g %g\n", 1.0/t, f);

        fprintf(fphi, " %g %g\n", 1.0/t, phi);

        dphi = - (1.0 / t) * f * phi;

        df = - 0.5 * (1.0 /(t*t*t)) * (phi * phi - 1.0);


        /* increment the functions */
        f -=  df * delta;
        phi -= dphi * delta;
    }
    printf("f = %g, df = %g, phi = %g, dphi = %g\n",
           f, df, phi, dphi);

    /* we try to write in the same file */
    fprintf(fa, "##\n");
    fprintf(fa, "## going another direction\n");
    fprintf(fa, "##\n");
    fprintf(fphi, "##\n");
    fprintf(fphi, "## going another direction\n");
    fprintf(fphi, "##\n");

#if 0
    fclose(fa);
    fclose(fphi);


    /*****
     * Calculate in the opposite direction
     *****/    
    fa = fopen("profile_zero.a", "w+");
    if (fa == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fphi = fopen("profile_zero.phi", "w+");
    if (fphi == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fprintf(fa, "#  r   f(r) \n");
    fprintf(fphi, "#  r   phi(r) \n");
#endif


    printf("============ going another direction ============\n");
    /* the value at t = t1 */
    phi = phi_init;
    f = f_init;
    delta = (t2 - t1)/(npoints - 1.0);
    for (i = 1; i <= npoints; i++)
    {
      double ii = (double)(i);

      double t = t1 + delta * (ii - 1.0);
      double df, dphi;


      dphi = - (1.0 / t) * f * phi;

      df = - 0.5 * (1.0 /(t*t*t)) * (phi * phi - 1.0);


      /* increment the functions */
      f +=  df * delta;
      phi += dphi * delta;

      t += delta;  /* we actually calculated in a further point */
      fprintf(fa, " %g %g\n", 1.0/t, f);
      fprintf(fphi, " %g %g\n", 1.0/t, phi);
    }
    fclose(fa);
    fclose(fphi);
}
