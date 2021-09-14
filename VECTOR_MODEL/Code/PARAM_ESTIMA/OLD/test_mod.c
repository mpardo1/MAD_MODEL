/* file age3classp.c */
#include <R.h>
static double parms[4];

#define a parms[0]
#define b parms[1]
#define c parms[2]
#define d parms[3]


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=4;
  odeparms(&N, parms);
}


/* Derivatives */
  void derivs(double *t, double *y, double *ydot,
               double *yout)
{

    ydot[0] = a*y[0] - b*y[1] + c;  //L
    ydot[1] = d*y[1];	     //A
}
/* END file age3classp.c */
