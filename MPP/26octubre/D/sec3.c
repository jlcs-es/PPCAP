/*
  CPP_CONTEST=PP-FP2014
  CPP_PROBLEM=D
  CPP_LANG=C+OPENMP+MPI
  CPP_PROCESSES_PER_NODE=calisto 4
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

#define NUM_THREADS 4

 
// multiplicación de matrices secuencial
// por cada matriz aparece la zona de datos (a, b y c)
// y el número de filas, de columnas y el leading dimension
void mms(double *a, int fa, int ca, int lda, double *b, int fb, int cb, int ldb, double *c, int fc, int cc, int ldc) {
    int i, j, k;
    double s;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for private(i, j, k, s) schedule(dynamic)
    for (i = 0; i < fa; i++) 
        for (j = 0; j < cb; j++) {
            s = 0.;
            for (k = 0; k < ca; k++)
                s += a[i * lda + k] * b[k * ldb + j];
            c[i * ldc + j] = s;
        }
}

// nodo es un identificador del proceso
// y np el número total de procesos
void sec(int n,double *a, double *b,double *c,int nodo, int np) {
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    int i,fa=n,ca=n,lda=n,fb=n,cb=n,ldb=n,fc=n,cc=n,ldc=n,position, prevposition;
    double s;

    if(nodo!=0)
    {
        a=(double *) malloc(sizeof(double)*n*n/np);
        b=(double *) malloc(sizeof(double)*n*n);
        c=(double *) malloc(sizeof(double)*n*n/np);
    }

    char * buffer = (char *) malloc(sizeof(double)*(n*n + n*n/np));

    if (nodo == 0) {
        position = 0;
        MPI_Pack(b, fb*cb, MPI_DOUBLE, buffer, sizeof(double)*(n*n + n*n/np), &position, MPI_COMM_WORLD);
        prevposition = position;
        for (i = 1; i < np; i++){
            MPI_Pack(&a[i * lda * fa / np], fa / np * ca, MPI_DOUBLE, buffer, sizeof(double)*(n*n + n*n/np), &position, MPI_COMM_WORLD);
            MPI_Send(buffer, sizeof(double)*(n*n + n*n/np), MPI_PACKED, i, 20, MPI_COMM_WORLD);
            position = prevposition;
        }
    } else {
        MPI_Recv(buffer, sizeof(double)*(n*n + n*n/np), MPI_PACKED, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        position = 0;
        MPI_Unpack(buffer, sizeof(double)*(n*n + n*n/np), &position, b, fb*cb, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(buffer, sizeof(double)*(n*n + n*n/np), &position, a, fa / np * ca, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    free(buffer);

    mms(a, fa / np, ca, lda, b, fb, cb, ldb, c, fc / np, cc, ldc);
    MPI_Gather(c, fc / np * cc, MPI_DOUBLE,c, fc / np * cc, MPI_DOUBLE,0, MPI_COMM_WORLD);
   if(nodo!=0)
   {
       free(a);
       free(b);
       free(c);
   }
}

