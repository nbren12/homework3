#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define STOP_ITER_RAT 10e-6


double calc_resid(int N, double* f, double *u){
  int i, j;
  double resid = 0.0;
  double ai;

  double h2 = 1.0/(N+1)/(N+1);
  
#pragma omp parallel for reduction(+:resid) 
  for (i = 1; i < N+1; i++) {
    ai = (- u[i+1] + 2 * u[i] - u[i-1])/h2 -f[i];
    resid += ai*ai;
  }
  
  return sqrt(resid);
}

void gauss_seidel_laplace_redblack(int N, int start, double *f, double *u){
  int i;
  double h2 = 1.0/(N+1)/(N+1);

  #pragma omp parallel for
  for (i = 1 + start; i < N+1; i+=2) {
    u[i] = (h2 * f[i] + u[i-1] + u[i+1])/2.0;
  }
}

void test_resid(){
  int N = 100;
  
  double *u, *f;

  // allocate arrays
  u = (double *) malloc(N*sizeof (double));
  f = (double *) malloc(N*sizeof (double));
  
  // initialize f and u
  int i;
  for (i = 0; i < N; i++) {
    f[i] = 2.0;
    u[i] = 0.0;
  }
  
  double resid;
  
  resid = calc_resid(N, f, u);
  printf("Resid is %f", resid);
  
}

int main(int argc, char *argv[])
{
  

  int N = atoi(argv[1]);
  int maxIter, iter;
  maxIter = atoi(argv[2]);
  iter = 0;
  
  double *u, *f, *w;

  // allocate arrays
  u = (double *) malloc((N+2)*sizeof (double));
  f = (double *) malloc((N+2)*sizeof (double));
  
  // initialize f and u
  int i;
  for (i = 0; i < N+2; i++) {
    f[i] = 1.0;
    u[i] = 0.0;
  }

  // Begin iterations
  double resid_init, resid_cur;
  resid_init = calc_resid(N, f, u);
  resid_cur = resid_init;
  
  timestamp_type t1, t2;

  get_timestamp(&t1);
  while (resid_cur / resid_init > STOP_ITER_RAT){
    u[0] = u[N+1] = 0.0;
    
    gauss_seidel_laplace_redblack(N, 0, f, u);
    gauss_seidel_laplace_redblack(N, 1, f, u);

    resid_cur = calc_resid(N, f, u);
    printf("Resid is %f\n", resid_cur );
    if (++iter >= maxIter ) break;
  }
  get_timestamp(&t2);


  printf("Total time: %f\n", timestamp_diff_in_seconds(t1,t2));

  // deallocate
  free(f);
  free(u);
  
  return 0;
}
