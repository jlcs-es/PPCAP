#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "io.h"
#include <omp.h>


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

void multiplicar_acumular_matrices(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc)
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
        s=s+da[k]*db[kb];
      }
      c[i*ldc+j]=s;
    }
  }
}

void matriz_matriz_bloques(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc,int tb)
{
  int i,j,k, fblocka, cblocka, cblockb;
  double *s, *ds;

  s=(double *) malloc(sizeof(double)*tb*tb*omp_get_max_threads());

  omp_set_num_threads(omp_get_max_threads());

  #pragma omp parallel for private(i,j,k,fblocka, cblocka, cblockb, ds)
  for(i=0;i<fa;i=i+tb)
  {
    ds = &s[tb*tb*omp_get_thread_num()];
    fblocka = ( (i+tb) < fa ) ? tb : (fa-i);
    for(j=0;j<cb;j=j+tb)
    {
      cblockb = ( (j+tb) < cb ) ? tb : (cb-j);
      matriz_cero(ds,tb,tb,tb);
      for(k=0;k<ca;k=k+tb)
      {
        cblocka = ( (k+tb) < ca ) ? tb : (ca-k); // cblocka == fblockb
        multiplicar_acumular_matrices(&a[i*lda+k],fblocka,cblocka,lda,&b[k*ldb+j],cblocka,cblockb,ldb,ds,fblocka,cblockb,tb);
      }
      copiar_matriz(ds,fblocka,cblockb,tb,&c[i*ldc+j],fblocka,cblockb,ldc);
    }
  }
  free(s);
}

int main()
{
  double *a,*b,*c,l,u;
  int fa,ca,lda,fb,cb,ldb,fc,cc,ldc,tb;
  
  int ti1,tf1,si,sf;
  struct timeval *tv;
  struct timezone *tz;

  tv=(struct timeval *) malloc(sizeof(struct timeval));
  tz=(struct timezone *) malloc(sizeof(struct timezone));

  printf("De las filas y columnas de la primera matriz y las columnas de la segunda :");
  scanf("%d",&fa);
  scanf("%d",&ca);
  lda=ca;
  fb=ca;
  scanf("%d",&cb);
  ldb=cb;
  fc=fa;
  cc=cb;
  ldc=cc;
  printf("De el tama�o de bloque (divisor de los tama�os de las matrices): ");
  scanf("%d",&tb);
  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);

  a=(double *) malloc(sizeof(double)*fa*ca);
  b=(double *) malloc(sizeof(double)*fb*cb);
  c=(double *) malloc(sizeof(double)*fc*cc);

  generar_vector(a,fa*ca,l,u);
  generar_vector(b,fb*cb,l,u);
  
#ifdef DEBUG
  escribir_matriz_ld(a,fa,ca,lda);
  escribir_matriz_ld(b,fb,cb,ldb);
#endif

  gettimeofday(tv,tz);
  si=(tv->tv_sec);
  ti1=(tv->tv_usec);
  ti1=si*1000000+ti1;

  matriz_matriz_bloques(a,fa,ca,lda,b,fb,cb,ldb,c,fc,cc,ldc,tb);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);

  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tama�o matriz %d, bloque %d: %.6f seg\n",fa,tb,(tf1-ti1)/1000000.);

#ifdef DEBUG
  escribir_matriz_ld(c,fc,cc,ldc);
#endif

 
  free(tv);
  free(tz);
  free(a);
  free(b);
  free(c);

  return 0;
}





























