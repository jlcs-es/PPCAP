/*
CPP_CONTEST=2017
CPP_PROBLEM=B
CPP_LANG=C+MPI
CPP_PROCESSES_PER_NODE=calisto 4
*/

#include <stdlib.h>
#include <mpi.h>


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
    if (nodo == 0) {
        int *zerosrows=new int[n],*zeroscolumns=new int[n];

        trasponer(n,b);
        obtainzerosrows(n,a,zerosrows);
        obtainzerosrows(n,b,zeroscolumns);
    }
    if(nodo!=0)
    {
        a=(double *) malloc(sizeof(double)*n*n/np);
        b=(double *) malloc(sizeof(double)*n*n);
        c=(double *) malloc(sizeof(double)*n*n/np);
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.;
            for (int k = 0; k < (zerosrows[i]<zeroscolumns[j]?zerosrows[i]:zeroscolumns[j]); k++)
                s += a[i * n + k] * b[j * n + k];
            c[i * n + j] = s;
        }
    delete[] zerosrows;
    delete[] zeroscolumns;
}

