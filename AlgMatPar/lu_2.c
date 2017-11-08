#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "io.h"
#include <omp.h>

#define NUM_THREADS 4


double fabs();




void matriz_matriz_ld(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc)
{
  int i,j,k;
  double s;

  for(i=0;i<fa;i++)
  {
    for(j=0;j<cb;j++)
    {
      s=0.;
      for(k=0;k<ca;k++)
        s+=a[i*lda+k]*b[k*ldb+j];
      c[i*ldc+j]=s;
    }
  }
}



void lu(double *a,int fa,int ca,int lda)
{
  int i,j,k;
  omp_set_num_threads(NUM_THREADS);

  for(i=0;i<fa;i++)
  {
    #pragma omp parallel for private(j) schedule(static, 64)
    for(j=i+1;j<ca;j++)
      a[i*lda+j]/=a[i*lda+i];
    #pragma omp parallel for private(j,k) schedule(static, 64)
    for(j=i+1;j<fa;j++)
      for(k=i+1;k<ca;k++)
        a[j*lda+k]-=a[j*lda+i]*a[i*lda+k];
  }
}

copiar_matriz(double *mo,int fo,int co,int ldo,double *md,int fd,int cd,int ldd)
{
  int i,j;

  for(i=0;i<fo;i++)
    for(j=0;j<co;j++)
      md[i*ldd+j]=mo[i*ldo+j];
}

comprobar_lu(double *b,int fb,int cb,int ldb,double *a,int fa,int ca,int lda)
{
  int i,j,k;
  double s,*l,*u,*p;

  l=(double *) calloc(fa*fa,sizeof(double));
  u=(double *) calloc(fa*ca,sizeof(double));
  p=(double *) malloc(sizeof(double)*fa*ca);

  for(i=0;i<fa;i++)
    for(j=0;j<=i;j++)
      l[i*lda+j]=a[i*lda+j];
  for(i=0;i<fa;i++)
  {
    u[i*lda+i]=1.;
    for(j=i+1;j<ca;j++)
      u[i*lda+j]=a[i*lda+j];
  }
  
  matriz_matriz_ld(l,fa,fa,fa,u,fa,ca,ca,p,fa,ca,ca);

  for(i=0;i<fa;i++)
    for(j=0;j<ca;j++)
      if(fabs(b[i*lda+j]-p[i*lda+j])>0.00000001)
        printf("Error en %d, %d: %.6lf, %.6lf\n",i,j,b[i*lda+j],p[i*lda+j]);

  free(l);
  free(u);
  free(p);
}

main()
{
  double *a,*copia;
  int fa,ca,lda;
  double l,u;

  int ti1,tf1,si,sf;
  struct timeval *tv;
  struct timezone *tz;

  tv=(struct timeval *) malloc(sizeof(struct timeval));
  tz=(struct timezone *) malloc(sizeof(struct timezone));

  printf("De la dimension de la matriz: ");
  scanf("%d",&fa);
  lda=ca=fa;

  a=(double *) malloc(sizeof(double)*fa*ca);
  copia=(double *) malloc(sizeof(double)*fa*ca);

  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);

  generar_matriz_ld(a,fa,ca,lda,l,u);
  copiar_matriz(a,fa,ca,lda,copia,fa,ca,lda);

#ifdef DEBUG
  //escribir_matriz_ld(a,fa,ca,lda);
#endif

  gettimeofday(tv,tz);
  si=(tv->tv_sec);
  ti1=(tv->tv_usec);
  ti1=si*1000000+ti1;

  lu(a,fa,ca,lda);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);
  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tama�o %d: %.6lf seg \n",fa,(tf1-ti1)/1000000.);

#ifdef DEBUG
  //escribir_matriz_ld(a,fa,ca,lda);
  comprobar_lu(copia,fa,ca,lda,a,fa,ca,lda);
#endif
 
  free(tv);
  free(tz);
  free(a);
}







