 // File with the vector field F de x'=F(t,x)
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
    for(i=0;i<dim+1;i++) printf("\n h[i] :\n %lf ", arr_spline[i][4]);
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
    double x_i1,x_i;
    
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
            x_i = x[i];
            x_i1 = x[i+1];
      //      printf("\n time:%lf, x[i+1]:%lf,x[i]%lf",time,x[i+1],x[i]);
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

// function which computes the moving average for rainfall taking into account the last 12 months, if the time is less than 12 months it takes the average of the pasts months.
double moving_average(double *fx, double t){
    int i;
    double mov_avg;
    int for_st, for_end;
    
    mov_avg = 0;
    
    if( t < 1){
        mov_avg = fx[0];
    }
    else if ( t < 365 && t >= 1){
        //printf("\n ****Mov avg: time between 12 and 1");
        mov_avg = 0;
        for_end = round(t);
        for(i=0; i < for_end; i++){
            mov_avg = mov_avg + fx[i];
        }
        mov_avg = mov_avg /for_end;
    }
    else if( t >= 365){
        //printf("\n ****Mov avg: time greater than 12");
        for_st = round(t-365);
        for_end = round(t);
        //printf("\n part_int_1:%lf", part_int_1);
        mov_avg = 0;
        //printf("\n for_st:%lf ,for_end:%lf", for_st, for_end);
        for(i=for_st; i<for_end; i++){
            mov_avg = mov_avg + fx[i];
        }
        mov_avg = mov_avg /for_end;
    }
    //printf("\n mov_avg: %lf", mov_avg);

    return mov_avg;
}
/*
// Function with the vector field.
void mapF(int dim, double t, double *x, double *y, double **arr_spline,double *y2, double **arr_spl_rain, int numparam, double *param, int dim_temp, int dim_rain, double *Fx,double seasons)
{
    double sigma, rho, beta, delta_T;
    double a, beta_e, sigma_0, LAMBDA, delta_H, gamma;
    double n_H, epsilon , eta, ni, F, delta_L;
    double N, M, lambda_T, delta_M, B, r, r_0;
    double d_L, c, y_p, gamma_P, n_P, k_A, K_E;
    double Tem,P, delta_0, delta_P, delta_R,b_min;
    double omega;
    double avg_rain;
    int i;
    double norm;
    // Fixed parameters from a paper: "Temperature-related duration of aquatic stages of the
    // Afrotropical malaria vector mosquito Anopheles gambiae in the laboratory"
    double delta_10 = 1/2.7;
    double delta_12 = 1/3.7;
    double delta_14 = 1/20.5;
    double delta_16 = 1/25.5;
    double delta_18 = 1/24.9;
    double delta_20 = 1/24.9;
    double delta_22 = 1/18.1;
    double delta_24 = 1/16.4;
    double delta_26 = 1/13.5;
    double delta_28 = 1/11;
    double delta_30 = 1/11.2;
    double delta_32 = 1/10.2;
    double delta_34 = 1/8.9;
    double delta_36 = 1/6.9;
    double delta_38 = 1/4.8;
    double delta_40 = 1/2.8;
        
    
    //printf("\n\n Dentro del vector Field al principios");
    //for(i=0; i<10;i++){
    //    printf("\n x[%d] = %lf", i, x[i]);
    //}
    //printf("\n\n");
    // Initializes variables to zero.
    
    if( t != t ){
           printf("\nTime t is Nan, t:%15.12e",t);
           for(i=0; i<10;i++){
               printf("\n x[%d] = %lf", i, x[i]);
           }
           printf("\n Total de la poblacion N:%lf",x[0]+x[1]+x[2]+x[3]+x[4]);
           exit(35);
       }
    
    
    Tem = 0;
    omega = 0;
    delta_T = 0;
    delta_P = 0;
    a = 0;
    
    Tem = interp_temp(t, dim_temp, y, arr_spline);
    //printf("\n t:%15.12e  T: %lf",t, Tem);
    
    d_L = 0;
    if (0.00554*Tem-0.06737 > 0){
        d_L = 0.00554*Tem-0.06737;
    }
    
    //d_L = 0.0988299;
    
    gamma_P = 0;
    if (0.009*Tem - 0.1441 > 0){
        gamma_P = 0.009*Tem - 0.1441;
    }
      
    
    //gamma_P = 0.1258999;
    
    // First cont parameters.
    b_min = param[0]; //Probability M-->H.
    beta_e = param[1]; // External force of infection.
    sigma_0 = param[2]; // Initial Value of loss of inmunity.
    K_E = param[3]; //Loss rate.
    delta_H = param[4]; // Mortality rate humans.
    rho = param[5]; // Rate of recovery.
    gamma = param[6]; // Average time in the exposed phase.
    n_H = param[7]; //Exposed number.
    epsilon = param[8]; // Franctions of humas which develops severe symtoms.
    eta = param[9]; // Case Probability.
    ni = param[10]; // Recovery rate.
    F = param[11]; // Mosquito fecundity factor.
    delta_0 = param[12]; // Basic mortality Rate for Larvae. ---> (Tomás: tenías un 0. Conviene buscar valores reales.)
                         // De https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4233045/ he sacado delta_0 aprox DMR=daily mortality rate
                         // de la larva del anopheles gambiae igual a 0.1)  
    delta_R = param[13];//Death factor due to rain peaks.
    r_0 = param[14]; //Recovery rate.
    c = param[15]; //Probability H-->M.
   
    k_A = param[17]; //Conversion factor.
    n_P = param[18]; // Exposed number P
    
    

    //Second parameters which are parametrized by a function.
    N = x[0] + x[1] + x[2] + x[3] + x[4]; // Total human Population.
    //printf("\n N:%20.17lf", N);
    y_p = (x[2]+x[4])/N; // fractions of infectious humans.
    M = x[6] + x[7] + x[8]; // Total mosquito Population,
    lambda_T = -4.4 + 1.31*Tem - 0.03*pow(Tem,2.);
    
    if( fabs(lambda_T) < 1.e-8 ){
        printf("\n lambda_T close to zero :%15.12e ",lambda_T);
        for(i=0; i<10;i++){
            printf("\n x[%d] = %lf", i, x[i]);
        }
        printf("\n Total de la poblacion N:%lf",N );
        exit(50);
    }
    
    delta_M = 1./lambda_T; // Mortality rate mosquitos.
    //delta_M = 0.126582278;
    
    // Seasons 0: dry, 1:rainy
    if (seasons == 0){
        a = 1.2;
    }else{
        a = 3.7;
    }
    //a = 0.95214;
    
    LAMBDA = a*x[8]/N; // Rate of infectious bites per human.
    
    if( fabs(N) < 1.e-8 ){
        printf("\n Total de la poblacion close to zero N:%lf",N );
        for(i=0; i<10;i++){
            printf("\n x[%d] = %lf", i, x[i]);
        }
        exit(51);
    }
    
    //LAMBDA = 0.190428;
    //printf("\n");
    //printf("\n-----****** LAMBDA:%lf, a:%lf, x[8]:%15.12e, N:%lf\n", LAMBDA,a,x[8],N);
    B = delta_H*N; // Births and immigrations.
    beta = b_min*(1/N)*x[8]*a + beta_e; // Transmision rate.
    //beta = 57.93179;
    //printf("\n exp(LAMBDA/sigma_0) :%lf", exp(LAMBDA/sigma_0));
    sigma = LAMBDA/(exp(LAMBDA/sigma_0)-1); // Loss of immunity.
    if( fabs((exp(LAMBDA/sigma_0)-1)) < 1.e-8 ){
       
        sigma = sigma_0;
    }
    //sigma = 0.019733;
    r = LAMBDA/(exp(LAMBDA/r_0)-1); // Recovery rate.
    //r = 0.42273;
     if((exp(LAMBDA/r_0)-1) < 1.e-8 ){
      
         r=r_0;
     }
    
    if( Tem >= 10 && Tem <= 12 ){
        delta_T = (delta_10+delta_12)/2;
    }else if( Tem > 12 && Tem<=14 ){
        delta_T = (delta_12+delta_14)/2;
    }else if( Tem > 14 && Tem<=16 ){
        delta_T = (delta_14+delta_16)/2;
    }else if(Tem> 16 &&Tem<=18 ){
        delta_T = (delta_16+delta_18)/2;
    }else if(Tem> 18 &&Tem<=20 ){
        delta_T = (delta_18+delta_20)/2;
    }else if(Tem> 20 &&Tem<=22 ){
        delta_T = (delta_20+delta_22)/2;
    }else if(Tem> 22 &&Tem<=24 ){
        delta_T = (delta_22+delta_24)/2;
    }else if(Tem> 24 &&Tem<=26 ){
        delta_T = (delta_24+delta_26)/2;
    }else if(Tem> 26 &&Tem<=28 ){
        delta_T = (delta_26+delta_28)/2;
    }else if(Tem> 28 &&Tem<=30 ){
        delta_T = (delta_28+delta_30)/2;
    }else if(Tem> 30 &&Tem<=32 ){
        delta_T = (delta_30+delta_32)/2;
    }else if(Tem> 32 &&Tem<=34 ){
        delta_T = (delta_32+delta_34)/2;
    }else if(Tem> 34 &&Tem<=36 ){
        delta_T = (delta_34+delta_36)/2;
    }else if(Tem> 36 &&Tem<=38 ){
        delta_T = (delta_36+delta_38)/2;
    }else if(Tem> 38 &&Tem<=40 ){
        delta_T = (delta_38+delta_40)/2;
    }
    
    P = interp_temp(t, dim_rain, y2, arr_spl_rain);
    if (P<0) {  // Tomás
      printf("\n P negativo! \n");
      exit(666); 
      }
    
    avg_rain = moving_average(y2,t);

    
    if((P-avg_rain)>0){
        omega = (P-avg_rain);
    }
 
    
    
    if( fabs(x[8]) < 1.e-8 ){
        printf("\n No infectious mosquitos Tem=%lf,a=%15.12e, W =%15.12e",Tem,a, x[8]);
        for(i=0; i<10;i++){
            printf("\n x[%d] = %lf", i, x[i]);
        }
        printf("\n Total de la poblacion N:%lf",N );
        //exit(48);
    }
    
    delta_P = delta_R*omega;
    delta_L = delta_0 + delta_T + delta_P;
    
   
    Fx[0] = B - beta*x[0] + sigma*x[3]- delta_H*x[0]+rho*x[4];
    //printf("\n *****x[0]:\tB:%lf, beta=%lf,\n\t   sigma=%lf , LAMBDA:%15.12e ,\n\t  delta_H=%lf, rho:%lf",B,beta,   sigma,LAMBDA,delta_H,rho);
     //printf("\n\t  B:%lf,x[0]:%lf, \n\t x[3]=%lf, x[4]:%lf\n", B,x[0],x[3],x[4]);
    // S : Susceptible Patients.
    Fx[1] = beta*x[0] - delta_H*x[1] - gamma*n_H*x[1];
     //printf("\n *****x[1]:\t beta=%lf,   delta_H=%lf,\n\t  gamma=%lf, n_H:%lf",beta,delta_H,gamma,n_H);
     //printf("\n\t  x[0]:%lf, x[1]=%lf\n", x[0],x[1]);
    // E : Population Infected.
    Fx[2] = (1 - epsilon)*gamma*n_H*x[1] - eta * beta*x[2] + ni*x[4] - r*x[2] - delta_H*x[2]; // I : E without symptoms.
     //printf("\n *****x[2]:\t epsilon=%lf,   gamma=%lf  ,\n\t  eta=%lf, beta:%lf, ni:%lf",epsilon,gamma,eta,beta,ni);
     //printf("\n\t  x[1]:%lf, x[2]=%lf,\n\t  x[4]=%lf\n", x[1],x[2],x[4]);
    Fx[3] = -sigma*x[3] + r*x[2] - delta_H*x[3]; // R : Recovered.
     //printf("\n *****x[3]:\t sigma=%lf,   r=%lf  , \n\t delta_H=%lf",sigma,r,delta_H);
     //printf(" x[3]:%lf,\n\t  x[2]=%lf\n", x[3],x[2]);
    Fx[4] = epsilon*gamma*n_H*x[1] + eta*beta*x[2] - ni*x[4] - rho*x[4] - delta_H*x[4]; // C : E with symptoms.
    Fx[5] = F*a*M*((x[9] - x[5])/x[9]) - delta_L*x[5] - d_L*x[5]; // L : Larvae Mosquito.
    //printf("\n\n eq-------> Fx[5] = F*a*M*((x[9] - x[5])/x[9]) - delta_L*x[5] - d_L*x[5]");
    if (a == 0){
        printf("\n");
        printf("\n temp:%lf",Tem);}
        
    Fx[6] = -c*a*y_p*x[6] - delta_M*x[6] + d_L*x[5]; // X : Non-infected mosquitos.
     // printf("\n");
     // printf("\n *****x[6]:\t x[6]=%lf, x[5]=%lf,\n\t  c=%lf,a=%lf,\n\t  y_p=%lf,delta_M:%lf,\n\t  d_L=%lf\n",x[6], x[5],c,a,y_p,delta_M,d_L);
    Fx[7] = c*a*y_p*x[6] - gamma_P*n_P*x[7] - delta_M*x[7]; // V : Infected non infectious mosquitos.
      //printf("\n *****x[7]:\tgamma_P=%lf,   n_P=%lf,\n\t  delta_M=%lf, Temp:%lf\n",gamma_P,n_P,delta_M,Tem);
    Fx[8] = gamma_P*n_P*x[7] - delta_M*x[8]; // W : Infectious mosquitos.
   
    Fx[9] = k_A*P - K_E*x[9]; // K : Carrying capacity.
    
    norm = 0;
    
        if (fabs(Fx[i]) > norm){
            norm = fabs(Fx[i]);
        }
    }
    if (norm < 1.e-06 ){
        printf("\n *****Norma del campo inferior a 0.1");
        exit(123);
    }
    
}*/

// Function with the vector field.
void mapF(int dim, double t, double *x, double *y, double **arr_spline,double *y2, double **arr_spl_rain, int numparam, double *param, int dim_temp, int dim_rain, double *Fx,double seasons)
{
    
    double a, gamma;
    int i;
 	
    a = param[0]; //Probability M-->H.
    gamma = param[1]; 
    //P = interp_temp(t, dim_rain, y2, arr_spl_rain);
   
    
    Fx[0] = a - (1 + gamma)*x[0];
    Fx[1] = x[0] - (1 + gamma)*x[1];
    Fx[2] = x[1] - (1 + gamma)*x[2];
    Fx[3] = x[2] - (1 + gamma)*x[3];
    Fx[4] = x[3] - (1 + gamma)*x[4]; 
    Fx[5] = x[4] - (1 + gamma)*x[5];
    
}
