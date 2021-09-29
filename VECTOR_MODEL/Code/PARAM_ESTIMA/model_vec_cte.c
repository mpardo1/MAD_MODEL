/* file age3classp.c */
#include <R.h>
static double parms[7];
static double forc[1];

#define fecun parms[0]
#define Ka parms[1]
#define del_L parms[2]
#define del_A parms[3]
#define dev_L parms[4]
#define gon parms[5]
#define Hum parms[6]

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
    	ydot[0] =  gon*fecun *y[2]*(1-(y[0]/Ka))-(dev_L+del_L)*y[0]; // L
	ydot[1] =  dev_L*y[0] - i1*Hum - del_A*y[1];   // A
	ydot[2] =   i1*Hum- (gon + del_A)*y[2];	// Ah	
	
}
/* END file age3classp.c */
