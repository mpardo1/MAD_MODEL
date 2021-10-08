/* file age3classp.c */
#include <R.h>
static double parms[4];
static double forc[4];

#define fecun parms[0]
#define Ka parms[1]
#define Hum parms[2]
#define rho parms[3]


#define i1 forc[1]  // del_L
#define i2 forc[2]  // del_A
#define i3 forc[3]  // dev_L
#define i4 forc[4]  // gono


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=4;
  odeparms(&N, parms);
}

/* Forcings initializers */
void forcc(void (* odeforcs)(int *, double *))
{
  int N=4;
  odeforcs(&N, forc);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
    	ydot[0] =  i4*fecun*y[2]*(1-(y[0]/Ka))-(i3+i1)*y[0]; //L 
	ydot[1] =  i3*y[0] - rho*Hum*y[1] - i2*y[1];	//A
	ydot[2] =  rho*Hum*y[1] - (i4 + i2)*y[2];	        //Ah
}
/* END file age3classp.c */
