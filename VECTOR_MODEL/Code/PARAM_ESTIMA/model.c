/* file age3classp.c */
#include <R.h>
static double parms[4];
static double forc[4];

#define f parms[0] 
#define K parms[1]
#define H parms[2]
#define omega parms[3]

#define i1 forc[0] // a(t)
#define i2 forc[1] // d_L(t)
#define i3 forc[2] // delta_L(t)
#define i4 forc[3] // delta_A(t)


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=4;
  odeparms(&N, parms);
}
void forcc(void (* odeforcs)(int *, double *))
{
  int N=4;
  odeforcs(&N, forc);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
/*
  ydot[0] = i1*f*y[2]*(1 - (y[0]/K)) - 0.04*y[0] - i3*y[0];     //L
  ydot[1] = 0.04*y[0] - i4*H - i5*y[1];				//A
  ydot[2] = i4*H - i1*y[2] - i5*y[2];				//A^(h)
  */
  ydot[0] = i1*f*y[2]*(1 - (y[0]/K)) - i2*y[0] - i3*y[0];     //L
  ydot[1] = i2*y[0] - omega*y[1]*H - i4*y[1];				//A
  ydot[2] = omega*y[1]*H  - i1*y[2] - i4*y[2];				//A^(h)
}
/* END file age3classp.c */
