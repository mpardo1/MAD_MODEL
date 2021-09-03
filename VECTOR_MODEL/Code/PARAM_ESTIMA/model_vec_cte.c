/* file age3classp.c */
#include <R.h>
static double parms[4];

#define f parms[0] 
#define K parms[1]
#define H parms[2]
#define omega parms[3]


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=4;
  odeparms(&N, parms);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
  ydot[0] = 0.05*f*y[2];//*(1 - (y[0]/K)) - 0.05*y[0] - 0.05*y[0];     //L
  ydot[1] = 0.05*y[0] - omega*y[1]*H - 0.05*y[1];		    //A
  ydot[2] = omega*y[1]*H  - 0.05*y[2] - 0.05*y[2];	            //A^(h)
}
/* END file age3classp.c */
