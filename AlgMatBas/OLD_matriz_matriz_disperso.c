#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "io.h"


void trasponer_especial(double *db,int *fb,int *cb,int ndb, double *dr, int *fr, int *cr){
  int i, cact, csig, csigsig;
  cact = cb[0];
  for(int j=1; j<ndb; j++) {
    if(cact > cb[j]){
      cact = cb[j];
    }
  }
  i=0;
  while(i<ndb & cact < INT_MAX) {
    csig = INT_MAX;
    for(int j=0; j<ndb; j++) {
      if(cact == cb[j]){
        dr[i] = db[j];
        fr[i] = cact;
        cr[i] = fb[j];
        i++;
        if(i>=ndb)
          return;
      } else if ( cb[j] > cact && cb[j] < csig) {
        csig = cb[j];  // siguiente columna a buscar
      }
    }
    cact = csig;
  }

}


void matriz_vector_disperso_ld(double *dm,int *fm,int *cm,int ndm,double *dv,int *fv,int ndv,double *r,int ldr)
{
  int i,j,fact;
  double s;

  i=0;
  while(i<ndm)
  {
    fact=fm[i]; // fila actual
    j=0;
    s=0.;
    while(i<ndm && fm[i]==fact)
    {
      while(j<ndv && fv[j]<cm[i])
        j++;
      if(j<ndv && fv[j]==cm[i])
        s+=dm[i]*dv[j];
      r[fact*ldr]=s;
      i++;
    }
  }
}


void matriz_matriz_disperso(double *da,int *fa,int *ca,int nda, 
                            double *db,int *fb,int *cb,int ndb,
                            double *c,int fc,int cc,int ldc) 
{

  // Suponemos fm ordenado y cm ordenado en bloques respecto de fm
  double *dbt = (double *) malloc(sizeof(double)*ndb);
  int *fbt = (int *) malloc(sizeof(int)*ndb);
  int *cbt = (int *) malloc(sizeof(int)*ndb);
  trasponer_especial(db, fb, cb, ndb, dbt, fbt, cbt);

  int col_act_init, col_act, col_act_fin;
  int j=0;
  while(j<ndb){
    col_act_init=j;
    col_act = fbt[j];
    while(j<ndb && fbt[j]==col_act)
      j++;
    col_act_fin=j; //points to the last element + 1

    matriz_vector_disperso_ld(da, fa, ca, nda, &dbt[col_act_init], &cbt[col_act_init], (col_act_fin-col_act_init), &c[col_act], ldc);

  }

}


int main(int argc,char **argv)
{
  int filas_a,col_a,filas_b,col_b;
  int nda, ndb;
  double l,u;

  int ti1,tf1,si,sf;
  struct timeval *tv;
  struct timezone *tz;

  tv=(struct timeval *) malloc(sizeof(struct timeval));
  tz=(struct timezone *) malloc(sizeof(struct timezone));


  printf("De los tamaños de la matriz A: ");
  scanf("%d %d",&filas_a,&col_a);
  filas_b=col_a;
  printf("De las columnas de la matriz B: ");
  scanf("%d",&col_b);
  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);
  printf("De el numero de datos distintos de cero en la matriz A: ");
  scanf("%d",&nda);
  printf("De el numero de datos distintos de cero en la matriz B: ");
  scanf("%d",&ndb);

  double *da = (double *) malloc(sizeof(double)*nda);
  int *fa = (int *) malloc(sizeof(int)*nda);
  int *ca = (int *) malloc(sizeof(int)*nda);
  double *db = (double *) malloc(sizeof(double)*ndb);
  int *fb = (int *) malloc(sizeof(int)*ndb);
  int *cb = (int *) malloc(sizeof(int)*ndb);
  double *c = (double *) calloc(sizeof(double),filas_a*col_b); // init with 0's

  generar_matriz_dispersa(da,filas_a,col_a,fa,ca,nda,l,u);
  generar_matriz_dispersa(db,filas_b,col_b,fb,cb,ndb,l,u);
#ifdef DEBUG
  escribir_matriz_dispersa(da,fa,ca,nda);
  escribir_matriz_dispersa(db,fb,cb,ndb);
#endif


gettimeofday(tv,tz);
si=(tv->tv_sec);
ti1=(tv->tv_usec);
ti1=si*1000000+ti1;

matriz_matriz_disperso(da, fa, ca, nda, db, fb, cb, ndb, c, filas_a, col_b, col_b);


gettimeofday(tv,tz);
sf=(tv->tv_sec);
tf1=(tv->tv_usec);
tf1=sf*1000000+tf1;
printf("Tamaño %d: %.6lf seg \n",filas_a,(tf1-ti1)/1000000.);

#ifdef DEBUG
  escribir_matriz_ld(c,filas_a,col_b,col_b);
#endif

  free(da);
  free(fa);
  free(ca);
  free(db);
  free(fb);
  free(cb);
  free(c);

  return 0;
}
