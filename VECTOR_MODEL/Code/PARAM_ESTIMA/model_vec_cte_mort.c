/* file age3classp.c */
#include <R.h>
static double parms[5];
static double forc[3];

#define fecun parms[0]
#define Ka parms[1]
#define del_L parms[2]
#define del_A parms[3]
#define Hum parms[4]

#define i1 forc[0] 	// rho(A,t)
#define i2 forc[1]	// dev_L
#define i3 forc[2]	// gon


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=5;
  odeparms(&N, parms);
}

/* Forcings initializers */
void forcc(void (* odeforcs)(int *, double *))
{
  int N=3;
  odeforcs(&N, forc);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
    	ydot[0] =  i3*fecun*y[2]*(1-(y[0]/Ka))-(i2+del_L)*y[0];
	ydot[1] =  i2*y[0] - i1*Hum + del_A*y[1];
	ydot[2] =  i1*Hum - (i3 + del_A)*y[2];
}
/* END file age3classp.c */
