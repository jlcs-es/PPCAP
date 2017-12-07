#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "io.h"

double fabs();

/*void generar_matriz_ld(double *a,int n,int m,int ld,double l,double u)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      a[i*ld+j]=((double) rand()/RAND_MAX)*(u-l)+l;
}*/

void lu(double *a,int fa,int ca,int lda)
{
  int i,j,k;

  for(i=0;i<fa;i++)
  {
    for(j=i+1;j<ca;j++)
      a[i*lda+j]/=a[i*lda+i];
    for(j=i+1;j<fa;j++)
      for(k=i+1;k<ca;k++)
        a[j*lda+k]-=a[j*lda+i]*a[i*lda+k];
  }
}

void matriz_cero(double *m,int fm,int cm,int ldm)
{
  int i,j;
  double *dm;

  for(i=0;i<fm;i++)
  {
    dm=&m[i*ldm];
    for(j=0;j<cm;j++)
    {
      dm[j]=0.;
    }
  }
}

void copiar_matriz(double *mo,int fmo,int cmo,int ldmo,double *md,int fmd,int cmd,int ldmd)
{
  int i,j;
  double *dmo,*dmd;

  for(i=0;i<fmo;i++)
  {
    dmo=&mo[i*ldmo];
    dmd=&md[i*ldmd];
    for(j=0;j<cmo;j++)
    {
      dmd[j]=dmo[j];
    }
  }
}

void multiplicar_restar_matrices(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc)
{
  int i,j,k,kb;
  double *da,*db,s;

  for(i=0;i<fa;i++)
  {
    da=&a[i*lda];
    for(j=0;j<cb;j++)
    {
      db=&b[j];
      s=c[i*ldc+j];
      for(k=0,kb=0;k<ca;k++,kb=kb+ldb)
      {
        s=s-da[k]*db[kb];
      }
      c[i*ldc+j]=s;
    }
  }
}

void sistema_triangular_inferior(double *l,int fl,int cl,int ldl,double *x,int fx,int cx,int ldx)
{
  int i,j,k;

  for(i=0;i<fl;i++)
  {
    for(j=0;j<cx;j++)
      x[i*ldx+j]/=l[i*ldl+i];
    for(j=i+1;j<fl;j++)
      for(k=0;k<cx;k++)
        x[j*ldx+k]-=x[i*ldx+k]*l[j*ldl+i];
  }
}

void sistema_triangular_superior(double *u,int fu,int cu,int ldu,double *x,int fx,int cx,int ldx)
{
  int i,j,k;

  for(i=0;i<cu;i++)
  {
    for(j=i+1;j<cu;j++)
      for(k=0;k<fx;k++)
        x[k*ldx+j]-=x[k*ldx+i]*u[i*ldu+j];
  }
}

void lu_bloques(double *a,int fa,int ca,int lda,int tb)
{
  int i,j,k,f,c;

  for(i=0;i<fa;i=i+tb)
  {
    f=(tb<fa-i ? tb : fa-i);
    c=(tb<ca-i ? tb : ca-i);
    lu(&a[i*lda+i],f,c,lda);
    if(i+tb<fa)
    {
      sistema_triangular_inferior(&a[i*lda+i],f,c,lda,&a[i*lda+i+c],f,ca-i-c,lda);
#ifdef DEBUG
      printf("Tras sistema triangular inferior, %d, %d, %d, %d\n",i*lda+i,f,c,ca-i-c);
      escribir_matriz_ld(a,fa,ca,lda);
#endif
      sistema_triangular_superior(&a[i*lda+i],f,c,lda,&a[(i+f)*lda+i],fa-i-f,c,lda);
#ifdef DEBUG
      printf("Tras sistema triangular superior, %d, %d, %d, %d\n",i*lda+i,fa-i-f,f,c);
      escribir_matriz_ld(a,fa,ca,lda);
#endif
      multiplicar_restar_matrices(&a[(i+f)*lda+i],fa-i-f,c,lda,&a[i*lda+i+f],f,ca-i-c,lda,&a[(i+f)*lda+i+c],fa-i-f,ca-i-c,lda);
#ifdef DEBUG
      printf("Tras multiplicar matrices\n");
      escribir_matriz_ld(a,fa,ca,lda);
#endif
    }
  }
}

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

void comparar_matrices(double *b,int fb,int cb,int ldb,double *a,int fa,int ca,int lda)
{
  int i,j;

  for(i=0;i<fa;i++)
    for(j=0;j<ca;j++)
      if(fabs(b[i*lda+j]-a[i*lda+j])>0.00000001)
        printf("Error en %d, %d: %.6lf, %.6lf\n",i,j,b[i*lda+j],a[i*lda+j]);
}

main()
{
  double *a,*copia;
  int fa,ca,lda,tb;
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

  printf("De el tamaño de bloque: ");
  scanf("%d",&tb);

  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);

  generar_matriz_ld(a,fa,ca,lda,l,u);
  copiar_matriz(a,fa,ca,lda,copia,fa,ca,lda);

#ifdef DEBUG
  escribir_matriz_ld(a,fa,ca,lda);
#endif

  gettimeofday(tv,tz);
  si=(tv->tv_sec);
  ti1=(tv->tv_usec);
  ti1=si*1000000+ti1;

  lu_bloques(a,fa,ca,lda,tb);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);

  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tamaño matriz %d, bloque %d: %.6f seg\n",fa,tb,(tf1-ti1)/1000000.);

#ifdef DEBUG
  printf("Factorizacion LU por bloques:\n");
  escribir_matriz_ld(a,fa,ca,lda);
  lu(copia,fa,ca,lda);
  printf("Factorizacion LU sin bloques:\n");
  escribir_matriz_ld(copia,fa,ca,lda);
  comparar_matrices(a,fa,ca,lda,copia,fa,ca,lda);
#endif
 
  free(tv);
  free(tz);
  free(a);
  free(copia);
}





























