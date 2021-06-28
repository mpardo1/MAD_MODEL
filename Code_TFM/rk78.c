/*=====================================================================\
|     function solving  o.d.e. using  Runge-Kutta-Fehlberg method      |
|                          of order 7 and 8                            |
|----------------------------------------------------------------------|
|                                                                      |
|   double rk78 ( *at, x[], *ah, tol, hmin, hmax, n, deriv )           |
|                                                                      |
|     at          pointer to the variable TIME                         |
|     x[]         vector POSITION                                      |
|     y[]         vector Month                                         |
|     arr_spline  array with spline coeficients                        |
|     ah          pointer to the variable STEP                         |
|     tol         TOLERANCE                                            |
|     hmin, hmax  minimal and maximal step                             |
|     n           dimension                                            |
|     deriv       pointer to the function containing the VECTOR FIELD  |
|                                                                      |
|                 void  ( * deriv ) ( n, t, x[], dx[] )                |
|                                                                      |
|     it returns an estimation of the ERROR done                       |
|                 or -1 if the dimension is greater than NMAX          |
|                                                                      |
\=====================================================================*/


// inimem_rk78: allocates at the beginning the memory used for the rk78
// freemem_rk78: liberates such a memory
// rk78: integrator

// This method is two RK methods one of 7 order and another of 8 order, which computes de time step at each iteration
// with the difference between the x at RK7 and RK8 and take the better solution
//which will  be always the RK8 since it is high order.


#include "math.h"
#include "malloc.h"
#include "stdlib.h"
#include "stdio.h"
#include "subroutines.h"

#define sgn(a)   ((a)<0 ? -1. : 1. )

 static double alfa[ 13 ] =
     {            0.,      2. / 27.,       1. / 9.,       1. / 6.,
            5. / 12.,           .5 ,       5. / 6.,       1. / 6.,
             2. / 3.,       1. / 3.,            1.,            0.,
                  1.};

  static double beta[ 79 ] =
     {            0.,      2. / 27.,      1. / 36.,      1. / 12.,
            1. / 24.,            0.,       1. / 8.,       5. /12.,
                  0.,   - 25. / 16.,     25. / 16.,        .5e-1 ,
                  0.,            0.,          .25 ,           .2 ,
        - 25. / 108.,            0.,            0.,   125. / 108.,
         - 65. / 27.,    125. / 54.,    31. / 300.,            0.,
                  0.,            0.,    61. / 225.,     - 2. / 9.,
          13. / 900.,            2.,            0.,            0.,
          - 53. / 6.,    704. / 45.,   - 107. / 9.,     67. / 90.,
                  3.,  - 91. / 108.,            0.,            0.,
          23. / 108., - 976. / 135.,    311. / 54.,   - 19. / 60.,
            17. / 6.,    - 1. / 12., 2383. / 4100.,            0.,
                  0., - 341. / 164., 4496. / 1025.,  - 301. / 82.,
       2133. / 4100.,     45. / 82.,    45. / 164.,     18. / 41.,
           3. / 205.,            0.,            0.,            0.,
                  0.,    - 6. / 41.,   - 3. / 205.,    - 3. / 41.,
            3. / 41.,      6. / 41.,            0.,-1777. / 4100.,
                  0.,            0., - 341. / 164., 4496. / 1025.,
        - 289. / 82.,  2193. /4100.,     51. / 82.,    33. / 164.,
           12. / 41.,            0.,            1.};

  static double c7[ 11 ] =
     {    41. / 840.,            0.,            0.,            0.,
                  0.,    34. / 105.,      9. / 35.,      9. / 35.,
           9. / 280.,     9. / 280.,    41. / 840.};

  static double c8[ 13 ] =
     {            0.,            0.,            0.,            0.,
                  0.,    34. / 105.,      9. / 35.,      9. / 35.,
           9. / 280.,     9. / 280.,            0.,    41. / 840.,
          41. / 840.};




static double **krk, *x7, *x8, *xpon, *dx_rk;
static int n_equations=0;  // it controls the number of equations


// It allocates memory for the rk78 routine. n is the dimension

void inimem_rk78(int n)
{
  int i,j;
  if (n<1) {
    printf("\ninimem_rk78: dimension must be at least 1\n");
    exit(11);
  }
  if (n_equations!=0) {
    free(x7); // Deallocates memory, must be done each component if it is an array
    free(x8);
    free(xpon);
    free(dx_rk);
    for(j=0;j<n;j++) free(krk[j]); 
  }
  n_equations=n;

  //x7 [n]
  x7=(double *) calloc(n,sizeof(double));
  if(x7==NULL) {
    printf("\ninimem_rk78:It is not possible to allocate memory\n");
    exit(12);
  }
  for(j=0;j<n;j++) x7[j]=0.e+0;

  //x8 [n]
  x8=(double *) calloc(n,sizeof(double));
  if(x8==NULL) {
    printf("\ninimem_rk78: It is not possible to allocate memory\n");
    exit(13);
  }
  for(j=0;j<n;j++) x8[j]=0.e+0;

  //xpon [n]
  xpon=(double *) calloc(n,sizeof(double));
  if(xpon==NULL) {
    printf("\ninimem_rk78: It is not possible to allocate memory\n");
    exit(14);
  }
  for(j=0;j<n;j++) xpon[j]=0.e+0;

  //dx_rk [n]
  dx_rk=(double *) calloc(n,sizeof(double));
  if(dx_rk==NULL) {
    printf("\ninimem_rk78:It is not possible to allocate memory\n");
    exit(15);
  }
  for(j=0;j<n;j++) dx_rk[j]=0.e+0;


  //krk [n][13] 
  krk=(double **) calloc(n,sizeof(double*));
  if(krk==NULL) {
     printf("\ninimem_rk78:It is not possible to allocate memory\n");
     exit(16);
   }
  for(i=0;i<n;i++) {
    krk[i]=(double *) calloc(13,sizeof(double));
    if(krk[i]==NULL) {
     printf("\ninimem_rk78:It is not possible to allocate memory\n");
     exit(17);
   }
  }
  for(i=0;i<n;i++) {
    for(j=0;j<13;j++) {
      krk[i][j]=0.e+0;
    }
  }

  return;

}






// It liberates the memory allocated for the rk78

// n is the dimension of the system; it must coincide with the one used
// at inimem_rk78.c

void freemem_rk78(int n)
{
  int j;
  if (n!=n_equations) {
    printf("\nfreemem_rk78: dimensions do not coincide ! \n");
    exit(1);
  }
  free(x7);
  free(x8);
  free(xpon);
  free(dx_rk);
  for(j=0;j<n;j++) free(krk[j]); 
  n_equations=0;
}



// %%%%%%%%%%%%%%%%%%%%%%% RK78 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



double rk78( double *at, double x[], double *ah, double tol,
              double hmin, double hmax, int n, int numparam_rk, double *y_rk, double **arr_spline_rk,double *y2_rk, double **arr_spl_rain_rk, double *param_rk, int dim_temp_rk, int dim_rain_rk,int seasons,
              void ( *deriv ) ( int, double, double *,double *, double **,double *, double **, int, double * , int, int, double *,double))
{
     double tpon, tol1, err, nor, kh, beth, h1;
     double aux;
     register int j, l;
     int i, m;
    //int k;

     aux=0.e+0;
     tpon=0.e+0;
     tol1=0.e+0;
     kh=0.e+0;
     beth=0.e+0;
     h1=0.e+0;
    
     if (n>n_equations) {
       printf("\nrk78:error in dimensions (%d and %d)\n",n,n_equations);
       exit(1);
     }

     for(j=0;j<n;j++) {
       x7[j]=0.e+0;
       x8[j]=0.e+0;
       xpon[j]=0.e+0;
       dx_rk[j]=0.e+0;
     }
     for(i=0;i<n;i++) {
       for(j=0;j<13;j++) {
         krk[i][j]=0.e+0;
       }
     }
    
    // Check if pointers works:
    /*
    printf("\n\n -------> Comprobación Pointers antes llamar vector field rk79:");
    printf("\n dim:%d",n );
    printf("\n t:%15.12e",tpon );
    for(k=0; k<10;k++){
        printf("\n x[%d] = %lf", k, x[k]);
    }
    for(k=0; k<10;k++){
        printf("\n Fx[%d] = %lf", k, dx_rk[k]);
    }
     */


     do {
         /*---------------> Computing K(i,j) <-----------------------*/
         m = 0;
         for ( i = 0; i < 13; i++ )
            {
            tpon = *at + alfa[i] * *ah;
                for ( j = 0; j < n; j++){ xpon[j] = x[j];
                    //printf("\n Dentro del for donde xpon[%d]=%lf=x[%d]=%lf",j,xpon[j], j, x[j]);
                }
            for ( l = 0; l < i; l++ )
               {
           m=m+1;
               beth = *ah * beta[m];
               for ( j = 0; j < n; xpon[j] += beth * krk[j][l], j=j+1);
               }
            
                /*
            printf("\n\n -------> Comprobación Pointers antes de llamar vector field rk79:");
            printf("\n dim:%d",n );
            printf("\n i:%d",i);
            printf("\n t:%15.12e",tpon );
            for(k=0; k<10;k++){
                printf("\n xpon[%d] = %lf", k, xpon[k]);
                printf("\n x[%d] = %lf", k, x[k]);
            }
            for(k=0; k<10;k++){
                printf("\n Fx[%d] = %lf", k, dx_rk[k]);
            }
                 */
            ( *deriv )( n, tpon, xpon, y_rk, arr_spline_rk, y2_rk, arr_spl_rain_rk, numparam_rk, param_rk ,dim_temp_rk, dim_rain_rk, dx_rk,seasons);
            // Check if pointers works:
            /*
            printf("\n\n -------> Comprobación Pointers despues de llamar vector field rk79:");
            printf("\n dim:%d",n );
            printf("\n t:%15.12e",tpon );
            for(k=0; k<10;k++){
                printf("\n xpon[%d] = %lf", k, xpon[k]);
                printf("\n x[%d] = %lf", k, x[k]);
            }
            for(k=0; k<10;k++){
                printf("\n Fx[%d] = %lf", k, dx_rk[k]);
            }
                printf("\n Donde evalua krk");
             */
            for ( j = 0; j < n; krk[j][i] = dx_rk[j], j=j+1 );
            //for(k=0; k<10;k++){
            //    printf("\n krk[%d] = %lf", k, krk[k][1]);
            //}
            }


         /*------> Computing the 2 points and associated values <--------------*/
         err = 0.e+0;
         nor = 0.e+0;
         for ( j = 0; j < n; j++ )
            {
            x7[j] = x[j];
        x8[j] = x[j];
            for ( l = 0; l < 11; l++ )
               {
               kh = *ah * krk[j][l];
               x7[j] += kh * c7[l];
               x8[j] += kh * c8[l];
               }
            x8[j] += *ah * ( c8[11] * krk[j][11] + c8[12] * krk[j][12] );
            err += fabs ( x8[j] - x7[j] );
            nor += fabs ( x8[j] );
            }
         err /= n;
         /*---------------> Computing the new step size h <-------------------*/
         tol1 = tol * ( 1. + nor / 100. );
         if ( err < tol1 ) err = fmax ( err, tol1 / 256. );
         h1 = *ah;

         //printf("\rRK78: err=%15.12e,  tol1=%15.12e ,  h1=%15.12e",err,tol1,h1);


         *ah *= 0.9 * pow ( (double)(tol1 / err), 0.125 );
         if ( fabs ( *ah ) < hmin ) *ah = hmin * sgn ( *ah );
         if ( fabs ( *ah ) > hmax ) *ah = hmax * sgn ( *ah );
         /*
         printf("\n\n -------> Comprobación Pointers antes de acabar el bucle while rk79:");
         printf("\n dim:%d",n );
         printf("\n t:%15.12e",tpon );
         for(i=0; i<10;i++){
             printf("\n x[%d] = %lf", i, xpon[i]);
             printf("\n xpon[%d] = %lf", i, x[i]);
         }
         for(i=0; i<10;i++){
             printf("\n Fx[%d] = %lf", i, dx_rk[i]);
         }*/
           } while (( err >= tol1 ) && ( fabs ( *ah ) > hmin ));
      *at += h1;
      for ( j = 0; j < n; j++) x[j] = x8[j];
      return err;

}


