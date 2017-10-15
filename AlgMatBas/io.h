#include <stdio.h>
#include <stdlib.h>

void generar_vector(double *v,int n,double l,double u);

void escribir_vector(double *v,int n);

void generar_matriz(double **a,int n,int m,double l,double u);

void escribir_matriz(double **a,int n,int m);

void escribir_matriz_ld(double *a,int n,int m,int ld);

void generar_matriz_dispersa_fast(double *m,int filas,int columnas,int *fm,int *cm,int ndm,double l,double u);

void generar_matriz_dispersa(double *m,int filas,int columnas,int *fm,int *cm,int ndm,double l,double u);

void generar_vector_disperso(double *v,int filas,int *fv,int ndv,double l,double u);

void generar_matriz_ld(double *a,int n,int m,int ld,double l,double u);

void escribir_matriz_dispersa(double *m,int *fm,int *cm,int ndm);

void escribir_vector_disperso(double *v,int *fv,int ndv);

