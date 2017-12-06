#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "io.h"





void trasponer_matriz_esp(double *a,int fa,int ca,int lda,double *at,int fat, int cat,int ldat)
{
  int i,j;
    double t;
    
      for(i=0;i<fat;i++)
          for(j=0;j<cat;j++)
              {
                   at[i*ldat+j]=a[j*lda+i];
               }
}
                                   












void matriz_matriz_tras(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc)
{
  int i,j,k;
  double s;
  double *bt,*da,*db;

  bt=(double *) malloc(sizeof(double)*cb*fb);

  trasponer_matriz_esp(b,fb,cb,ldb,bt,cb,fb,fb);

  for(i=0;i<fa;i++)
  {
    for(j=0;j<cb;j++)
    {
      s=0.;
      da=&a[i*lda];
      db=&bt[j*fb];
      for(k=0;k<ca;k++,da++,db++)
        s+=da[0]*db[0];
      c[i*ldc+j]=s;
    }
  }
  
  free(bt);
}





main()
{
  double *a,*b,*c;
  int fa,ca,lda,fb,cb,ldb,fc,cc,ldc;
  double l,u;
  int i;

  int ti1,tf1,si,sf;
  struct timeval *tv;
  struct timezone *tz;

  tv=(struct timeval *) malloc(sizeof(struct timeval));
  tz=(struct timezone *) malloc(sizeof(struct timezone));

  printf("De las filas y columnas de la primera matriz y las columnas de la segunda: ");
  scanf("%d",&fa);
  scanf("%d",&ca);
  lda=fb=ca;
  scanf("%d",&cb);
  ldb=cb;
  fc=fa;
  ldc=cc=cb;

  a=(double *) malloc(sizeof(double)*fa*ca);
  b=(double *) malloc(sizeof(double)*fb*cb);
  c=(double *) malloc(sizeof(double)*fc*cc);

  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);

  generar_matriz_ld(a,fa,ca,lda,l,u);
  generar_matriz_ld(b,fb,cb,ldb,l,u);
#ifdef DEBUG
  escribir_matriz_ld(a,fa,ca,lda);
  escribir_matriz_ld(b,fb,cb,ldb);
#endif

  gettimeofday(tv,tz);
  si=(tv->tv_sec);
  ti1=(tv->tv_usec);
  ti1=si*1000000+ti1;

  matriz_matriz_tras(a,fa,ca,lda,b,fb,cb,ldb,c,fc,cc,ldc);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);
  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tamaño %d: %.6lf seg \n",fa,(tf1-ti1)/1000000.);

#ifdef DEBUG
  escribir_matriz_ld(c,fc,cc,ldc);
#endif
 
  free(tv);
  free(tz);
  free(a);
  free(b);
  free(c);
}




























