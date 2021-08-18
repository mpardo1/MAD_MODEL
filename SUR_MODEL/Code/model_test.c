/* file age3classp.c */
#include <R.h>
static double parms[1];
static double forc[1];

#define gam1 parms[0]

#define i1 forc[0] // a(t)


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
    
  int N=1;
  odeparms(&N, parms);
}
void forcc(void (* odeforcs)(int *, double *))
{
  int N=1;
  odeforcs(&N, forc);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
    ydot[0] = i1 - (gam1 + 1)*y[0];
    
}
/* END file age3classp.c */
