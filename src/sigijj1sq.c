#include <stdlib.h>
#include <R.h>
#include <R_ext/Rdynload.h>

//' Squared covariance and product of variance in C
 //'
 //' Calculation of unbiased estimate of square of sigmaijj1 and product sigma2_x sigma2_y.
 //' @param n an integer equal to length of \code{y} and \code{x}
 //' @param x a numeric vector
 //' @param y a numeric vector
 //' @param sigmaijj12 a numeric value
 //' @param ssijj1 a numeric value
 //'
 //' @export
 void sigijj12(int* n, double* x, double* y, double* sigmaijj12, double* ssijj1){
   *sigmaijj12 = 0;
   *ssijj1 = 0;
   for (int m1=0; m1< *n; m1++){
     for (int m2=0; m2< *n; m2++){
       for (int m3=0; m3< *n; m3++){
         for (int m4=0; m4< *n; m4++ ){
           if ((m1!=m2)&(m1 !=m3)&(m1 !=m4) &(m2 !=m3) &(m2 !=m4)  & (m3 !=m4)){
             *sigmaijj12 +=  ( (x[m1]-x[m2])*(y[m1]-y[m2])* (x[m3]-x[m4])*(y[m3]-y[m4]) );
             *ssijj1 += (x[m1]-x[m2])*(x[m1]-x[m2])*(y[m3]-y[m4])*(y[m3]-y[m4]);
           }
         }
       }
     }
   }
     *sigmaijj12 /= (4* (*n)*( *n-1)*( *n-2)*( *n-3) );
     *ssijj1 /= (4* (*n)*(*n-1)*(*n-2)*( *n -3) );
   }
