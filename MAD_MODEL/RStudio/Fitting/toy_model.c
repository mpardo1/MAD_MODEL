/* file toy_model.c */
#include <R.h>


// Variables to change depending on the number of equations and where and how many the breakpoints are.
// Number of "ages"
#define NUM_AGE 4
// Number of breakpoints.
#define NUM_BREAKS 2
// Location of the breakpoints.
#define LOC_BREAK_1 2
#define LOC_BREAK_2 20
#define LOC_BREAK_3 30

// Number of parameters (mus and gammas), warning change the R code aswell the vector parms.
#define NUM_PARAMS 4

static double parms[NUM_PARAMS];
static double forc[1];
#define mu1 parms[0]
#define gamma1 parms[1]
#define mu2 parms[2]
#define gamma2 parms[3]
#define mu3 parms[4]
#define gamma3 parms[5]

#define A forc[0]


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
  int N=NUM_PARAMS;
  odeparms(&N, parms);
}

void forcc(void (* odeforcs)(int *, double *))
{
  int N=1;
  odeforcs(&N, forc);
}


// Vector_field with 3 equations:
/*
void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip){
  ydot[0] = A - mu1*y[0] - gamma1*y[0];
  ydot[1] = mu1*y[0]  -  mu2*y[1] - gamma2*y[1];
  ydot[2] = mu1*y[1]  - gamma2*y[2];
}
*/

// Vector Field with n equations
void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip){
  int i;
  int age = NUM_AGE;
  int num_break = NUM_BREAKS;
  int loc_break1 = LOC_BREAK_1;
  int loc_break2 = LOC_BREAK_2;
  int loc_break3 = LOC_BREAK_3;
  
  ydot[0] = A - mu1*y[0] - gamma1*y[0];
  for (i = 1; i < age; i++){
  	if( num_break >= 1 && i <= loc_break1 ){
  	
  		ydot[i] = mu1*y[i-1]  -  mu1*y[i] - gamma1*y[i];
  		
  	}else if( num_break >= 2 && i > loc_break1 && i <= loc_break2 ){
  	
  		ydot[i] = mu1*y[i-1]  -  mu2*y[i] - gamma2*y[i];
  	}else{
  		
  		ydot[i] = mu2*y[i-1]  -  mu3*y[i] - gamma3*y[i];  	
  	}
  }
  ydot[age] = mu3*y[age-1]  - gamma3*y[age];
}

/* END file toy_model.c */
