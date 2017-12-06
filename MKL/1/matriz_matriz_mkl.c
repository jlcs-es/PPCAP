#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "io.h"
#include <mkl.h>



main()
{
  double *a,*b,*c;
  int fa,ca,lda,fb,cb,ldb,fc,cc,ldc;
  double l,u;
  int i;
  
  double alpha = 1.0; 
  double beta = 0.0;

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

  a = (double *)mkl_malloc( fa*ca*sizeof( double ), 64 );
  b = (double *)mkl_malloc( fb*cb*sizeof( double ), 64 );
  c = (double *)mkl_malloc( fc*cc*sizeof( double ), 64 );
  if (a == NULL || b == NULL || c == NULL) {
    printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
    mkl_free(a);
    mkl_free(b);
    mkl_free(c);
    return 1;
  }

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

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, fa, cb, ca, alpha, a, lda, b, ldb, beta, c, ldc);

  gettimeofday(tv,tz);
  sf=(tv->tv_sec);
  tf1=(tv->tv_usec);
  tf1=sf*1000000+tf1;
  printf("Tama�o %d: %.6lf seg \n",fa,(tf1-ti1)/1000000.);

#ifdef DEBUG
  escribir_matriz_ld(c,fc,cc,ldc);
#endif
 
  free(tv);
  free(tz);
  mkl_free(a);
  mkl_free(b);
  mkl_free(c);
}




























