#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "io.h"


void fast_traspose(double *db,int *fb,int *cb,int ndb, int filas, int col, double *dr, int *fr, int *cr){

  int i,colNum,j;
  int total[col], index[col];
  for(i=0;i<col;i++)
    total[i]=0;
  for(i=0;i<ndb;i++)
  {
    colNum=cb[i];
    total[colNum]++;
  }
  index[0]=0;
  for(i=1;i<col;i++)
    index[i]=index[i-1]+total[i-1];

  for(i=0;i<ndb;i++)
  {
    colNum=cb[i];
    j=index[colNum];
    index[colNum]++;
    fr[j]=cb[i];
    cr[j]=fb[i];
    dr[j]=db[i];
  }

}


int matriz_matriz_disperso(double *da,int *fa,int *ca,int nda, int filas_a, int col_a_filas_b,
  double *db,int *fb,int *cb,int ndb, int col_b,
  double *dc,int *fc,int *cc) 
{
  int i, j, k, ii, fact, cact;
  double sum;

  double *dbt = (double *) malloc(sizeof(double)*ndb);
  int *fbt = (int *) malloc(sizeof(int)*ndb);
  int *cbt = (int *) malloc(sizeof(int)*ndb);
  fast_traspose(db, fb, cb, ndb, col_a_filas_b, col_b, dbt, fbt, cbt);

  k=0; // iterar sparce matrix c
  i=0; // iterar sparse a
  while(i<nda) // Recorro filas de a
  {
    ii=i; // mantener fila de la multiplicación
    j=0;
    while(j<ndb) // Recorro filas de bt == columnas de b
    {
      fact=fa[i];
      cact=fbt[j];
      sum=0.;
      while(i<nda && j<ndb && fact==fa[i] && cact==fbt[j]) // Me mantengo en fila fact de a, columna cact de b (fila de bt)
      {
        if(ca[i]==cbt[j]) {
          sum += da[i] * dbt[j];
          i++; j++;
        } 
        else if(ca[i]<cbt[j])
                i++;
            else 
                j++;
      }
      if(sum>=0.000001) {
        dc[k]=sum;
        fc[k]=fact;
        cc[k]=cact;
        k++;
      }
      if(j<ndb)
        i=ii; // volver al inicio de la fila fact de a
      while(cact==fbt[j] && j<ndb) 
        j++; // ir a la siguiente columna de b / fila de bt
    }
    while(fact==fa[i] && i<nda)
        i++;
  }

  free(dbt);
  free(fbt);
  free(cbt);

  return k; // ndc

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
  double *dc = (double *) malloc(sizeof(double)*filas_a*col_b); // worst case scenario
  int *fc = (int *) malloc(sizeof(int)*filas_a*col_b);
  int *cc = (int *) malloc(sizeof(int)*filas_a*col_b);
  int ndc = filas_a*col_b;

  printf("Generando matrices \n");
  gettimeofday(tv,tz);
  si=(tv->tv_sec);
  ti1=(tv->tv_usec);
  ti1=si*1000000+ti1;
  
  generar_matriz_dispersa_fast(da,filas_a,col_a,fa,ca,nda,l,u);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);
  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tamaño %d: %.6lf seg \n",filas_a,(tf1-ti1)/1000000.);
  printf("Generada matriz A \n");  

  gettimeofday(tv,tz);
  si=(tv->tv_sec);
  ti1=(tv->tv_usec);
  ti1=si*1000000+ti1;

  generar_matriz_dispersa_fast(db,filas_b,col_b,fb,cb,ndb,l,u);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);
  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tamaño %d: %.6lf seg \n",filas_a,(tf1-ti1)/1000000.);
  printf("Generada matriz B \n");
  
#ifdef DEBUG
  escribir_matriz_dispersa(da,fa,ca,nda);
  escribir_matriz_dispersa(db,fb,cb,ndb);
#endif

printf("Entrando a matriz_matriz_disperso \n");

gettimeofday(tv,tz);
si=(tv->tv_sec);
ti1=(tv->tv_usec);
ti1=si*1000000+ti1;

ndc = matriz_matriz_disperso(da, fa, ca, nda, filas_a, col_a, db, fb, cb, ndb, col_b, dc, fc, cc);


gettimeofday(tv,tz);
sf=(tv->tv_sec);
tf1=(tv->tv_usec);
tf1=sf*1000000+tf1;
printf("Tamaño %d: %.6lf seg \n",filas_a,(tf1-ti1)/1000000.);

#ifdef DEBUG
  //escribir_matriz_ld(c,filas_a,col_b,col_b);
  escribir_matriz_dispersa(dc,fc,cc,ndc);
#endif

  free(da);
  free(fa);
  free(ca);
  free(db);
  free(fb);
  free(cb);
  free(dc);
  free(fc);
  free(cc);

  return 0;
}
