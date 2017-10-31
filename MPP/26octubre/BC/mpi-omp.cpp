/*
CPP_CONTEST=2017
CPP_PROBLEM=C
CPP_LANG=C+OPENMP+MPI
CPP_PROCESSES_PER_NODE=calisto 2
*/

#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#define NUM_THREADS 4


void obtainzerosrows(int n,double *a,int *zeros) {
    for(int i=0;i<n;i++) {
        int j=n-1;
        while(j>=0 and a[i*n+j]==0.)
            j--;
        zeros[i]=j+1;
    }
}


void trasponer(int n,double *a) {
    double temp;
    for(int i=0;i<n;i++) {
        for(int j=i+1;j<n;j++) {
            temp=a[i*n+j];
            a[i*n+j]=a[j*n+i];
            a[j*n+i]=temp;
        }
    }
}


void sec(int n,double *a,double *b,double *c,int nodo,int np) {
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    int *zerosrows, *zeroscolumns;
    if (nodo == 0) {
        zerosrows = new int[n];
        zeroscolumns = new int[n];
        trasponer(n,b);
        obtainzerosrows(n,b,zeroscolumns);
        obtainzerosrows(n,a,zerosrows);
        // for (int i = 1; i < np; i++){
        //     MPI_Send(&zerosrows[i * n / np], n / np, MPI_INT, i, 20, MPI_COMM_WORLD);
        //     MPI_Send(&a[i * n * n / np], n / np * n, MPI_DOUBLE, i, 20, MPI_COMM_WORLD);
        // }
    }
    if(nodo!=0)
    {
        a = new double[n*n/np];
        b = new double[n*n];
        c = new double[n*n/np]();
        zerosrows = new int[n/np];
        zeroscolumns = new int[n];
        // MPI_Recv(zerosrows, n / np, MPI_INT, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Recv(a, n / np * n, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Scatter(zerosrows, n / np, MPI_INT, zerosrows, n/np, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(a, n*n / np, MPI_DOUBLE, a, n*n/np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(zeroscolumns, n, MPI_INT, 0, MPI_COMM_WORLD);

    int i,j,k; double s;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for private(i, j, k, s) schedule(dynamic)
    for (i = 0; i < n/np; i++)
        for (j = 0; j < n; j++) {
            s = 0.;
            for (k = 0; k < (zerosrows[i]<zeroscolumns[j]?zerosrows[i]:zeroscolumns[j]); k++)
                s += a[i * n + k] * b[j * n + k];
            c[i * n + j] = s;
        }

    MPI_Gather(c, n / np * n, MPI_DOUBLE, c, n / np * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] zerosrows;
    delete[] zeroscolumns;
    if(nodo!=0){
        delete[] a;
        delete[] b;
        delete[] c;
    }
}
