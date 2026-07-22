#include <stdlib.h>
#include <R.h>
#include <R_ext/Rdynload.h>

void sigma4c(int* n, double* y, double* res);

//' sigma4 in C
//'
//' Calculation of unbiased estimate of $\sigma^4$ in C to speed up the loop
//'
//' @param n an integer equal to length of y
//' @param y a numeric vector
//' @param res a numeric value
//'
//' @export
void sigma4c(int* n, double* y, double* res){
  *res =0;
  int count=0;
  int i, j, k, l;
  for (i=0; i< *n-3; i++){
    for (j =i+1; j<*n-2; j++){
      for (k = j+1; k<*n-1; k++){
        for (l = k+1; l<*n; l++){
          *res += (y[i]-y[j])*(y[i]-y[j])*(y[k]-y[l])*(y[k]-y[l]);
          *res += (y[i]-y[k])*(y[i]-y[k])*(y[j]-y[l])*(y[j]-y[l]);
          *res += (y[i]-y[l])*(y[i]-y[l])*(y[k]-y[j])*(y[k]-y[j]);
          count +=3;
        }
      }
    }

  }
  *res = *res/count/4;
}

void sigijj12(int* n, double* x, double* y, double* sigmaijj12, double* ssijj1);



/* Array defining routines */

static const R_CMethodDef cMethods[] = {
  {"sigma4c", (DL_FUNC) &sigma4c, 3},  // 3 arguments in the c function
  {"sigijj12", (DL_FUNC) &sigijj12, 5},
  {NULL, NULL, 0}  /* Sentinel */
};


/* Package initialization function */
/*the C initialization function to match the package name */
/* (e.g., R_init_jjvec1 instead of R_init_Jjvec3) */

 void R_init_HeterFunctionalData(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
 }

