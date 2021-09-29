/* file age3classp.c */
#include <R.h>
static double parms[3];
static double forc[5];

#define fecun parms[0]
#define Ka parms[1]
#define Hum parms[2]

#define i1 forc[0]   // rho(A,t)
#define i2 forc[1]  // del_L
#define i3 forc[2]  // del_A
#define i4 forc[3]  // dev_L
#define i5 forc[4]  // gono


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=3;
  odeparms(&N, parms);
}

/* Forcings initializers */
void forcc(void (* odeforcs)(int *, double *))
{
  int N=5;
  odeforcs(&N, forc);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
    	ydot[0] =  i5*fecun*y[2]*(1-(y[0]/Ka))-(i4+i2)*y[0]; //L 
	ydot[1] =  i4*y[0] - i1*Hum - i3*y[1];	//A
	ydot[2] =  i1*Hum - (i5 + i3)*y[2];	        //Ah
}
/* END file age3classp.c */
