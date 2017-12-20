#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include <mkl.h>
#include <omp.h>

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
  double alpha = -1.0; 
  double beta = 1.0;
  //fprintf(stderr, "0.-->\n");

  #pragma omp parallel
  { 
    int t = omp_get_thread_num();
    int numt = omp_get_num_threads();
    
    double * da = &a[ t * lda * (fa / numt) ];
    double * dc = &c[ t * ldc * (fc / numt) ];
    int dfa = (t != numt-1) ? (fa / numt) : (fa / numt) + fa%numt;
    
    //fprintf(stderr, "1.--> t = %d - numt = %d - dfa = %d - row %d \n", t, numt, dfa, (t * (fa / numt)));
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dfa, cb, ca, alpha, da, lda, b, ldb, beta, dc, ldc);
    
    //fprintf(stderr, "2.--> t = %d\n", t);
  }
  
}


void sistema_triangular_inferior(double* a, int fa, int ca, int lda, double* b, int fb, int cb, int ldb){
		
	int i,j,k;
	#pragma omp parallel for private(i,j,k)
	for(j = 0; j < cb; ++j){
		for(i = 0; i < fa; ++i){
			double aux = 0.;
			for(k = 0; k < i; k++)
				aux += a[i*lda + k]*b[k*ldb+j];
			b[i*ldb + j] = (b[i*ldb + j] - aux)/a[i*lda + i];
		}
	}
}


void sistema_triangular_superior(double* a, int fa, int ca, int lda, double* b, int fb, int cb, int ldb){

	int i,j,k;
	#pragma omp parallel for private(i,j,k)
	for(j = 0; j < fb; ++j){
		for(i = 0; i < cb; ++i){	
			double aux = 0.;
			for(k = 0; k < i; ++k)
				aux += b[j*ldb + k]*a[k*lda+i];
			b[j*ldb + i] = (b[j*ldb + i] - aux);
		}
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

// Cortesía de Adrian:
void comprobar_lu(double *b,int fb,int cb,int ldb,double *a,int fa,int ca,int lda)
{
  int i,j,k;
  double s,*l,*u,*p;

  l=(double *) mkl_calloc(fa*fa,sizeof(double),64);
  u=(double *) mkl_calloc(fa*ca,sizeof(double),64);
  p=(double *) mkl_calloc(fa*ca,sizeof(double),64);

  for(i=0;i<fa;i++)
    memcpy(l+fa*i, a+lda*i, (i+1)*sizeof(double));
    
  for(i=0;i<fa;i++)
  {
    u[i*lda+i]=1.;
    for(j=i+1;j<ca;j++)
      u[i*lda+j]=a[i*lda+j];
  }

#ifdef DEBUG
  printf("L:\n");
  escribir_matriz_ld(l, fa, fa, fa);
  
  printf("U:\n");
  escribir_matriz_ld(u, fa, ca, ca);
#endif
  
  //multiplicar_acumular_matrices(l, fa, fa, fa, u, fa, ca, ca, p, fa, ca, ca);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        fa, ca, fa,
        1., l, fa, u, ca,
        0., p, ca
      );
  
#ifdef DEBUG
  printf("A:\n");
  escribir_matriz_ld(p, fa, ca, ca);
#endif

  for(i=0;i<fa;i++)
    for(j=0;j<ca;j++)
      if(fabs(b[i*lda+j]-p[i*lda+j])>0.00000001)
        fprintf(stderr, "lu_bloques_mkl_omp: Error en %d, %d: %.6lf, %.6lf\n",i,j,b[i*lda+j],p[i*lda+j]);

  mkl_free(l);
  mkl_free(u);
  mkl_free(p);
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

  a = (double *)mkl_malloc( fa*ca*sizeof( double ), 64 );
  copia = (double *)mkl_malloc( fa*ca*sizeof( double ), 64 );
  if (a == NULL || copia == NULL) {
    printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
    mkl_free(a);
    mkl_free(copia);
    return 1;
  }


  printf("De el tamaño de bloque: ");
  scanf("%d",&tb);

  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);

  generar_matriz_ld(a,fa,ca,lda,l,u);
  copiar_matriz(a,fa,ca,lda,copia,fa,ca,lda);

#ifdef DEBUG
  escribir_matriz_ld(a,fa,ca,lda);
#endif

  omp_set_num_threads( omp_get_max_threads() );
  // mkl_set_num_threads( mkl_get_max_threads() );
  
  //mkl_set_dynamic( 0 );
  //mkl_set_num_threads( 1 );
  //printf("omp_get_max_threads() = %d \n", omp_get_max_threads());


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
  // int nt = mkl_get_max_threads();
  // printf("mkl_get_max_threads() = %d", nt);

#ifdef DEBUG
  //printf("Factorizacion LU por bloques:\n");
  //escribir_matriz_ld(a,fa,ca,lda);
  lu(copia,fa,ca,lda);
  printf("Factorizacion LU sin bloques:\n");
  //escribir_matriz_ld(copia,fa,ca,lda);
  comparar_matrices(a,fa,ca,lda,copia,fa,ca,lda);
#endif

comprobar_lu(copia, fa, ca, lda, a,fa,ca,lda);
 
  free(tv);
  free(tz);
  mkl_free(a);
  mkl_free(copia);
}



