/******************************************************************************
* FILE: omp_bug4.c
* DESCRIPTION:
*   This very simple program causes a segmentation fault.
* AUTHOR: Blaise Barney  01/09/04
* LAST REVISED: 04/06/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1048

// A function for allocating 2d array
double** allocate2dArray(){
  double ** a;
  int i;
  a = (double **) calloc(N, sizeof(double*));
  for (i = 0; i < N; i++) {
    a[i] = (double*) calloc(N, sizeof(double));
  }

  return a;
}

int main (int argc, char *argv[]) 
{
int nthreads, tid, i, j;

/* Noah: each thread must dynamically allocate a because it is too
	big to fit in the stack */
double **a;




/* Fork a team of threads with explicit variable scoping */
#pragma omp parallel shared(nthreads) private(i,j,tid,a)
  {
  /* Noah: each thread must dynamically allocate a because it is too
  	    big to fit in the stack */
  
  a = allocate2dArray();
  /* Obtain/print thread info */
  tid = omp_get_thread_num();
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);

  /* Each thread works on its own private copy of the array */
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j] = tid + i + j;

  /* For confirmation */
  printf("Thread %d done. Last element= %f\n",tid,a[N-1][N-1]);

  }  /* All threads join master thread and disband */

}


