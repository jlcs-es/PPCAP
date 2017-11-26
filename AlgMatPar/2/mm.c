#include <stdlib.h>
#include <omp.h>

#define NUM_THREADS 2
    
void  sec(int n,double *a,double *b,double *c)
{
  int i, j, k;
  double s;

   omp_set_num_threads(NUM_THREADS);
   #pragma omp parallel for private(i, j, k, s) collapse(2) schedule(dynamic)
// cada thread trabaja con un bloque de filas de la matriz a
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            s = 0.;
            for (k = 0; k < n; k++)
                s += a[i * n + k] * b[k * n + j];
            c[i * n + j] = s;
        }
}

