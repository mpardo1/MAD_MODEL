/* file age3classp.c */
#include <R.h>
static double parms[7];
static double forc[1];

#define fecun parms[0]
#define Ka parms[1]
#define Hu parms[2]
#define del_L parms[3]
#define del_A parms[4]
#define dev_L parms[5]
#define gon parms[6]

#define i1 forc[0] // rho(A,t)


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=7;
  odeparms(&N, parms);
}

/* Forcings initializers */
void forcc(void (* odeforcs)(int *, double *))
{
  int N=1;
  odeforcs(&N, forc);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
    	ydot[0] =  gon*fecun *y[2]*(1-(y[0]/Ka))-(dev_L+del_L)*y[1];
	ydot[1] =  dev_L*y[0] - (i1*Hu + del_A)*y[1];
	ydot[2] =  i1*Hu*y[1] - (gon + del_A)*y[2];
}
/* END file age3classp.c */
