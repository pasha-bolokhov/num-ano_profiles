#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define NN       70

#define alpha    -3.0

double  AA[NN + 2];
double  CC[2 * NN + 2 + 1];  /* factorial is absorbed in the coeffs */

#define sqr(x)  ((x) * (x))


#define NPOINTS      2000
#define NPOINTS_INF  700
#define x0      0.0
#define x1      20.0
#define npoints ((double)(NPOINTS))

/*****************************************/
/* the value of f and \phi at the \infty */
/*****************************************/

/* PLAYING */

/*
 *  Hypothesis: it may approach a semi-integer number
 */
/* a(0) = 3.76 */
#define F_INIT   0.000030977
#define PHI_INIT 1.0 - 0.000002
#define x1       15.0


#if 0
/* a(0) = 4.997 */
#define F_INIT   0.000061955
#define PHI_INIT 1.0 - 0.000004
#define x1       15.0
#endif

#if 0
/* a(0) = 4.46 */
#define F_INIT   0.000046465
#define PHI_INIT 1.0 - 0.000003
#define x1       15.0
#endif

#if 0
/* a(0) = 5.43938 */
#define F_INIT   0.00007744
#define PHI_INIT 1.0 - 0.000005
#define x1       15.0
#endif

#if 0
/* 
 * Good result -- shows approaching to 12
 * a(0) = 12.9259
 */
#define F_INIT   0.00001024453
#define PHI_INIT 1.0 - 0.0000005
#define x1      20.0
#endif

#if 0
/* This value fits for n = 1, x1 = 20.0 */
#define F_INIT   0.0000000427545
#define PHI_INIT 1.0
#define x1      20.0
#endif

#if 0
/* This value fits for n = 2, phi_init = 1.0 */
#define F_INIT   0.000000133
#define PHI_INIT 1.0
#define x1       20.0
#endif

#if 0
/* the precision is not enough evidently to handle x1 = 40.0 */
#define x1       40.0
#define F_INIT   0.000000000000095192575
#define PHI_INIT 1.0
#endif

/* for n = 2 */
#if 0
#define x1       10.0
#define F_INIT   0.0019265
#define PHI_INIT 1.0
#endif


#if 0
/* for n = 3 */
#define x1       10.0
#define F_INIT   0.004
#define PHI_INIT 1.0
#endif


/* 
 * Here we'll store the profiles calculated
 */
double p_x[NPOINTS];    /* the x-coordinate   */
double p_phi[NPOINTS];  /* '\phi' profile     */
double p_f[NPOINTS];    /* 'f' profile        */

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
    double f    = 1.0;
    double phi  = 0.0;
    double dphi = 0.0;
    double df   = 0.0;

    AA[0] = 1.0;
    AA[1] = - 0.25;

    AA[2] = alpha;
#if 0
    AA[2] = alpha * sqr(AA[1]);
#endif

    /*****************************/
    /* now calculate the profile */
    /*****************************/

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

    /* the value at x = x1 */
    phi = PHI_INIT;
    f = F_INIT;
    for (i = NPOINTS - 1; i >= 0; i--)
    {
        int    k;
        double ii = (double)(i);
        double x = x0 + delta * ii;

        double res_f = 0.0;
        double res_phi = 0.0;
        double tt;

        double df, dphi;

	/* store the profile in memory */
	p_x[i] = x;
	p_f[i] = f;
	p_phi[i] = phi;

	/* store the profile to file */
        fprintf(fa, " %g %g\n", x, f);

        fprintf(fphi, " %g %g\n", x, phi);

        dphi = (1.0 / x) * f * phi;

        df = 0.5 * x * (phi * phi - 1.0);


        /* increment the functions */
        f -= df * delta;
        phi -= dphi * delta;
    }
    fclose(fa);
    fclose(fphi);
    printf("f = %g, df = %g, phi = %g, dphi = %g\n",
           f, df, phi, dphi);

    printf("the first 5 estimates of phi's power, x, est1, est2...\n");
    for (i = 1; i <= 5; i++)
    {
      double diff, dlogx, dphi;

      if (p_phi[i] < 0 ||
	  p_phi[i + 1] < 0)
      {
	printf(" BAD VALUE\n");
      }
      else
      {
	/* one method */
          diff = log(p_phi[i + 1]) - log(p_phi[i]);
          dlogx = log(p_x[i + 1]) - log(p_x[i]);

	  /* another method */
	  dphi = (p_phi[i + 1] - p_phi[i]) /
	    (p_x[i + 1] - p_x[i]);

          printf(" %g %g  %g %g %g %g\n", p_x[i], 
		 diff/dlogx,
		 p_x[i] * dphi / p_phi[i],
		 p_x[i] * dphi / p_phi[i + 1], 
		 p_x[i + 1] * dphi / p_phi[i],
		 p_x[i + 1] * dphi / p_phi[i + 1]
		 );
      }
    }

    /**********************/
    /** now calculate in the other (positive) direction **/
    /**********************/
    printf("*******************************\n");
    printf("FINAL RANGE\n");
    fa = fopen("profile_inf.a", "w+");
    if (fa == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fphi = fopen("profile_inf.phi", "w+");
    if (fphi == NULL)
    {
        fprintf(stderr, "can't open for write: %s\n", strerror(errno));
        exit(1);
    }
    fprintf(fa, "#  r   f(r) \n");
    fprintf(fphi, "#  r   phi(r) \n");

    /* the value at x = x1 */
    phi = PHI_INIT;
    f = F_INIT;
    printf("PHI_INIT = %e, F_INIT = %e\n", PHI_INIT, F_INIT);
    for (i = 0; i < NPOINTS_INF; i++)
      {
	double ii = (double)(i);
	double x = x1 + delta * ii;
	
	/* store the profile to file */
        fprintf(fa, " %g %e\n", x, f);

        fprintf(fphi, " %g %e\n", x, phi);

        dphi = (1.0 / x) * f * phi;

        df = 0.5 * x * (phi * phi - 1.0);
#if 0
	printf("x = %g, dphi = %g, df = %g\n", x, dphi, df);
#endif


        /* increment the functions */
        f += df * delta;
        phi += dphi * delta;
	
      }

    fclose(fa);
    fclose(fphi);
    printf("f = %g, df = %g, phi = %g, dphi = %g\n",
           f, df, phi, dphi);
}

