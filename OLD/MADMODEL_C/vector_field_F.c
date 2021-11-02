#include "malloc.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double *sub_diag, * sup_diag, *diag, *d ;

// Function which allocates memory for the variables used inside the interpolation method
void inimem_spline(int n)
{
    int j;

    if (n<1) {
        printf("\ninimem_spline: dimension must be at least 1\n");
        exit(15);
       }

    d=(double *) calloc(n,sizeof(double));
    if(d==NULL) {
        printf("\ninimem_spline:It is not possible to allocate memory\n");
        exit(16);
    }
    for(j=0;j<n;j++) {
        d[j]=0.e+0;
        //printf("%lf ", d[j]);
    }

    //printf("\n Dentro de la subrutina INIMEN despues de allocatar d\n");

    sub_diag=(double *) calloc(n,sizeof(double));
    if(sub_diag==NULL) {
        printf("\ninimem_spline:It is not possible to allocate memory\n");
        exit(17);
    }
    for(j=0;j<n;j++) sub_diag[j]=0.e+0;

    sup_diag=(double *) calloc(n,sizeof(double));
    if(sup_diag==NULL) {
        printf("\ninimem_spline:It is not possible to allocate memory\n");
        exit(18);
    }
    for(j=0;j<n;j++) sup_diag[j]=0.e+0;
    
    diag=(double *) calloc(n,sizeof(double));
    if(diag==NULL) {
        printf("\ninimem_spline:It is not possible to allocate memory\n");
        exit(19);
    }
    for(j=0;j<n;j++) diag[j]=0.e+0;

    return;

}

// Free Space of the Variables used inside the interpolation method.

void freemem_spline(int n)
{
    free(sub_diag);
    free(sup_diag);
    free(diag);
    free(d);
    printf("\freemem_spline:Memmory free\n");

}

// Función que resuelve un sistema tridiagonal.
/*
   solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
     x - initially contains the input vector v, and returns the solution x. indexed from 0 to X - 1 inclusive
     X - number of equations (length of vector x)
     a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive
     b - the main diagonal, indexed from 0 to X - 1 inclusive
     c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive
*/

void solve_tridiagonal_in_place_reusable(double *x, int N, double *a, double *b, double *c)
{
    int n;

    printf("\nDim of the system: %d\n", N);
    /* Allocate scratch space. */
    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    /* loop from 1 to N - 1 inclusive */
    for (n = 1; n < N; n++) {
        float m = 1.0f / (b[n] - a[n] * c[n - 1]);
        if(n < (N-1)) c[n] = c[n] * m;
        x[n] = (x[n] - a[n] * x[n - 1]) * m;
    }
    
    //for (n = 0; n < N; n++) printf("\n sol dentro x[%d] = %lf", n, x[n]);
    /* loop from N - 2 to 0 inclusive */
    for (n = N - 2; n >= 0; --n) {
        x[n] = x[n] - c[n] * x[n + 1];
    }
    //for (n=0; n< N; n++) printf("\n Sol x[%d] = %lf \n",n,x[n]);
    
}

// Function which compute the spline matrix with all the coefficients needed to compute the interpolate the curve.

int spline(int dim, double *fx, double **arr_spline)
{
    int i;

    printf("\n dim en el spline = %d",dim);
    for(i=0;i<dim+1;i++) printf("\n i:%d, h[i] :\n %lf ",i,  arr_spline[i][4]);
    for(i=0;i<dim+1;i++){
        if(arr_spline[i][4] == 0){
            printf("spline: There is a coefficient of h = 0, diving between zero\n");
            exit(21);
        }
    }
    
    // Compute the D vector.
    for(i=1; i < dim; i++){
        sub_diag[i] = arr_spline[i][4];
    }
    for(i=0; i < dim-1; i++){
        sup_diag[i] = arr_spline[i][4];
    }

    for(i=0; i < dim; i++){
        d[i] = 6 * ((( fx[i+2] - fx[i+1] ) / arr_spline[i+1][4] ) - (( fx[i+1] - fx[i] ) / arr_spline[i][4] )) ;
        diag[i] = 2*(arr_spline[i][4] + arr_spline[i+1][4]);
    }
    solve_tridiagonal_in_place_reusable(d, dim, sub_diag, diag, sup_diag);

    for(i=0; i < dim; i++) arr_spline[i+1][2] = d[i];

    for(i=0; i < dim; i++) arr_spline[i][3] = d[i];
    


    for(i=0; i < dim+1; i++){
        arr_spline[i][1] = (( fx[i+1] - fx[i] ) / arr_spline[i][4] ) - (arr_spline[i][3] - arr_spline[i][2])*(arr_spline[i][4]/6) - arr_spline[i][2]*(arr_spline[i][4]/2);
    }
    
    return i;
}
double interp_temp(double t, int dim, double *x, double **arr_spline){
    // dim: is the dimension of the input data.
    double temp;
    int i;
    int time1;
    double time;
    
    double part_dec;
    double part_int;
    int integ;
    
    
  //  printf("\n");
  //  printf("\n Dentro del interp_temp");
    /*
    printf("Dentro del interp_temp ***");
    for (int i = 0; i < dim-1; i++) {
        for (int j = 0; j < 5; j++) {
            printf("\n arr_spline[%d,%d] = %lf ",i,j,arr_spline[i][j]);
        }
        printf("\n");
    }
     */
    
    //printf("\n t%15.12e in interpret", t);
    part_dec = modf(t, &part_int);
    integ = (int) part_int;
    //printf("-----> Interp: integ:%d, part_dec:%lf", integ,part_dec);
    
    time1 = integ % 366;
    time = time1 + part_dec;
    //printf("-----> Interp: time:%lf", time);
    time = t;
  //  printf("\n dim-1:%d",dim-1);
    temp = 1.0;
    for(i=0; i < (dim); i++){
        //printf("\n x[%d]:%lf", i, x[i]);
        if(time<=x[i+1] && time>=x[i]){
            
            //printf("\n time:%lf, x[i+1]:%lf,x[i]%lf",time,x[i+1],x[i]);
            temp = arr_spline[i][0]+arr_spline[i][1]*(time-x[i])+arr_spline[i][2]*(pow(time-x[i],2.)/2.)+((arr_spline[i][3]-arr_spline[i][2])/arr_spline[i][4])*((pow(time-x[i],3.))/6.);
            break;
        }
    }
    //if (temp_min == 0){
     //   temp = arr_spline[0][0]+arr_spline[0][1]*(t-x[0])+arr_spline[0][2]*(pow(t-x[0],3)/2)+((arr_spline[0][3]-arr_spline[0][2])/arr_spline[0][4])*((pow(t-x[0],3))/6);
   // }
   // printf("\n temp: %15.12e, t:%15.12e", temp,t);

    if (temp < 0 ){
        temp = ((arr_spline[i+1][0]-arr_spline[i][0])/(x[i+1]-x[i]))*(time-x[i])+arr_spline[i][0];
    //    printf("\n temp: %15.12e, t:%15.12e", temp,t);
    }
    

    return temp;
}
/*

// Function with the vector field.
void mapF(int dim, double t, double *x, int numparam, double *param, double *y, double **arr_spline, int dim_mos, double *Fx)
{
    double A;   			// Download rate, Total human population, turnover rate
    double mu1, mu, mu_prev;			// Transitions rate
    double gamma1, gamma, gamma4, gamma_last;
    double a, b0, b1, b2;			// Probabilities
    double M;					// Mosquito Density (As a first approach as a cte)
    double **array_break;
    int num_breaks; 				// Number of breakpoints 
    int i, j;
      
    A = param[0];	// Download rate
    a = param[1]; 	// Encounter rate
    b0 = param[2]; 	// p_0(A/B)
    b1 = param[3]; 	// p_1(A/B)
    b2 = param[4]; 	// p_2(A/B)
    mu1 = param[5];
    gamma1 = param[6];
    num_breaks = param[7]; // Number of breakpoints.
    gamma_last = param[numparam-1];
    
    
    // Allocate memory for breakpoints: first colum x_i, second: mu, third: gamma.
    array_break = (double **)malloc(num_breaks * sizeof(double *));
        for (i=0; i<dim_mos; i++)
            array_break[i] = (double *)malloc(3 * sizeof(double));

        if (array_break==NULL) {
            printf("\nIt is not possible to allocate memory for array_break.\n");
            exit(8);
        }

        for (i = 0; i < num_breaks; i++)
        for (j = 0; j < 3; j++)
            array_break[i][j] = 0.;
    
    
    // Save the values of x for the breakpoints.
    for(i=0; i < num_breaks; i++){
    	array_break[i][0] = param[8 + 3*i];
    	array_break[i][1] = param[9 + 3*i];
    	array_break[i][2] = param[10 + 3*i];
    }
     
    printf("Number of breakpoints: %d",num_breaks);
    for (i = 0; i < num_breaks; i++)
        printf("\n i: %d, x_i =%lf, mu_i=%lf, gamma_i = %lf ",i, array_break[i][0],array_break[i][1],array_break[i][2]);
    printf("\n");
    // Compute the value of the mosquito density by the spline cubics method for time t.
    M = interp_temp(t, dim_mos, y, arr_spline);
    
    
    Fx[0] = A - mu1*x[0] - gamma1*x[0];  	   // P1 : Participants 1º age group.
    for(i=1; i < dim-2; i++){
    	if(num_breaks >= 1 && i <= array_break[0][0] ){
    		mu = array_break[0][1];
    		gamma = array_break[0][2];
    		Fx[i] = mu1*x[i-1] - mu*x[i] - gamma*x[i];  	          // P2 : Participants 2º age group.
    	}else if(num_breaks >= 2 && i > array_break[0][0] && i <= array_break[1][0] ){
    		mu_prev = array_break[0][1];
    		mu = array_break[1][1];
    		gamma = array_break[1][2];
    		Fx[i] = mu_prev*x[i-1] - mu*x[i] - gamma*x[i];  	   // P2 : Participants 2º age group.
    	}else if(num_breaks >= 3 && i > array_break[1][0] && i <= array_break[2][0] ){
    		mu_prev = array_break[1][1];
    		mu = array_break[2][1];
    		gamma = array_break[2][2];
    		Fx[i] = mu_prev*x[i-1] - mu*x[i] - gamma*x[i];  	   // P2 : Participants 2º age group.
    	}
    }
    Fx[dim-2] = mu*x[2] - gamma_last*x[3];                                 // P4 : Participants 4º age group.
    Fx[dim-1] = b0*a*M;
    //x[] + b1*a*M*x[2] + b2*a*M*x[3] ;    // R : Reports.
}



// Function with the ODE solution.
void solF(double t, int numparam, double *param, double *Gx){
	
    double A;   					// Download rate, Total human population, turnover rate
    double mu1, mu2, mu3;				// Transitions rate
    double gamma1, gamma2, gamma3, gamma4;	// Mortality rate
    double C1, C2, C3, C4;	
    double alpha1, alpha2, alpha3;
    
    A = param[0];	// Download rate
    mu1 = param[1];
    mu2 = param[2];
    mu3 = param[3];
    gamma1 = param[4];
    gamma2 = param[5];
    gamma3 = param[6];
    gamma4 = param[7];
    
    alpha1 = mu1+gamma1;
    alpha2 = mu2+gamma2;
    alpha3 = mu3+gamma3;
    
    C1 = (A*mu1*mu2*mu3)/(alpha1*(alpha1-alpha2)*(alpha1-alpha3)*(alpha1-gamma4));
    C2 = (A*mu1*mu2*mu3)/(alpha2*(alpha1-alpha2)*(alpha2-alpha3)*(-alpha2+gamma4));
    C3 = (A*mu1*mu2*mu3)/(alpha3*(alpha1-alpha3)*(alpha2-alpha3)*(alpha3-gamma4));
    C4 = (A*mu1*mu2*mu3)/(gamma4*(gamma4-alpha3)*(alpha2-gamma4)*(alpha1-gamma4));
    
    //printf("\n t = %lf, C1:%lf, C2:%lf, C3:%lf, C4:%lf",t,C1,C2,C3,C4);
    // Equations for the solution:
    if( fabs(alpha1 - alpha2) < 1.e-10){
        printf("\n Parameters are  equal");
        Gx[0] =  (A- A*exp(-(gamma1 + mu1)*t))/(gamma1 + mu1);// P1 : Participants 1º age group. 
        Gx[1] =  (A*exp(-(gamma1 + mu1)*t)*mu1*(-1 + exp((gamma1 + mu1)*t) - gamma1*t - mu1*t))/pow((gamma1 + mu1),2);  // P2 : Participants 2º age group. 
        Gx[2] =  -(A*exp(-(gamma1 + mu1)*t)*pow(mu1,2)*(2 - 2*exp((gamma1 + mu1)*t) + 2*mu1*t + pow(gamma1*t,2) + pow(mu1*t,2) + 2*gamma1*t*(1 + mu1*t)))/(2*pow((gamma1 + mu1),3));
        Gx[3] =  (A*exp(-(gamma1 + mu1)*t)*(2*exp(mu1*t)*((-1 + exp(gamma1*t))*pow(mu1,3)) + (pow(gamma1,3)*(2 - 2*exp(mu1*t) + 2*mu1*t + pow(mu1*t,2)) + 2*pow(gamma1,2)*mu1*(3 - 3*exp(mu1*t) + 3*mu1*t + pow(mu1*t,2)) + gamma1*pow(mu1,2)*(6 - 6*exp(mu1*t) + 4*mu1*t + pow(mu1*t,2)))))/(2*gamma1*pow((gamma1 + mu1),3));
    }else{

        Gx[0] = C1*exp(-alpha1*t)*((alpha2-alpha1)*(alpha1-alpha3)*(alpha1-gamma4))/(mu1*mu2*mu3) + A/alpha1;									   // P1 : Participants 1º age group. 
        Gx[1] = C1*exp(-alpha1*t)*(((alpha1-alpha3)*(alpha1-gamma4))/(mu2*mu3)) + C2*exp(-alpha2*t)*(((alpha2-alpha3)*(alpha2-gamma4))/(mu2*mu3)) + (mu1/(alpha1*alpha2))*A;     // P2 : Participants 2º age group.
        Gx[2] = C1*exp(-alpha1*t)*((gamma4-alpha1)/mu3) + C2*exp(-alpha2*t)*((gamma4-alpha2)/mu3) + C3*exp(-alpha3*t)*((gamma4-alpha3)/mu3) + (mu1*mu2/(alpha1*alpha2*alpha3))*A;// P3 : Participants 3º age group.
        Gx[3] = C1*exp(-alpha1*t) + C2*exp(-alpha2*t) + C3*exp(-alpha3*t) + C4*exp(-gamma4*t) + (mu1*mu2*mu3/(alpha1*alpha2*alpha3*gamma4))*A;                         	    // P4 : Participants 4º age group.
    }
	
}
*/

void mapF(int dim, double t, double *x, int numparam, double *param, double *y, double **arr_spline, int dim_mos, double *Fx){
	

    double A;   					// Transitions rate
    double gamma;	// Mortality rate
  
    
    A = param[0];	// Download rate
    gamma = param[1];
    

    
    //printf("\n t = %lf, C1:%lf, C2:%lf, C3:%lf, C4:%lf",t,C1,C2,C3,C4);
    // Equations for the solution:
    
       Fx[0] = A - (gamma + 1)*x[0];
        Fx[1] = x[0] - (1 + gamma)*x[1];   
        Fx[2] = x[1] - (1 + gamma)*x[2];
        Fx[3] = x[2] - (1 + gamma)*x[3];  
	Fx[4] = x[3] - (1 + gamma)*x[4];
    	Fx[5] = x[4] - (1 + gamma)*x[5]; 
}

