/* file age3classp.c */
#include <R.h>
static double parms[8];

#define fec parms[0]
#define Ka parms[1]
#define H parms[2]
#define omeg parms[3]
#define del_L parms[4]
#define del_A parms[5]
#define dev_L parms[6]
#define gon parms[7]

/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=8;
  odeparms(&N, parms);
}


/* Derivatives */
  void derivs(double *t, double *y, double *ydot,
               double *yout)
{

    ydot[0] = gon*fec*y[2]*(1-(y[0]/Ka))-(dev_L+del_L)*y[0];  //L
    ydot[1] = dev_L*y[0] - (omeg*H + del_A)*y[1];	     //A
    ydot[2] = omeg*H*y[1] - (gon + del_A)*y[2];	     //Ah

}
/* END file age3classp.c */
