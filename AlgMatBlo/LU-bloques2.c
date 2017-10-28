#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "io.h"

double fabs();
void lu(double *a,int fa,int ca,int lda);
void comprobar_lu(double *b,int fb,int cb,int ldb,double *a,int fa,int ca,int lda);
void matriz_matriz_ld(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc);
void copiar_matriz(double *mo,int fo,int co,int ldo,double *md,int fd,int cd,int ldd);
void matriz_cero(double *m,int fm,int cm,int ldm);
void multiplicar_restar_matrices(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc);
void sistema_triangular_inferior(double *a,int fa,int ca,int lda,double *x,int fx,int cx,int ldx);
void sistema_triangular_superior(double *a,int fa,int ca,int lda,double *x,int fx,int cx,int ldx);

int tbm;


void lu_bloques(double *a,int fa,int ca,int lda, int tb)
{
    int i,j,k,f,c;

    for(i=0;i<fa;i=i+tb) {
        f= ( (i+tb) < fa ) ? tb : (fa-i);
        c= ( (i+tb) < ca ) ? tb : (ca-i);
        // Paso 1:
        lu(&a[i*lda+i],f,c,lda);
        if(i+tb<fa) {
            // Paso 2:
            sistema_triangular_inferior(&a[i*lda+i],f,c,lda,&a[i*lda+i+c],f,ca-i-c,lda);
            // Paso 3:
            sistema_triangular_superior(&a[i*lda+i],f,c,lda,&a[(i+f)*lda+i],fa-i-f,c,lda);
            // Paso 4:
            multiplicar_restar_matrices(&a[(i+f)*lda+i],fa-i-f,c,lda,&a[i*lda+i+c],f,ca-i-c,lda,&a[(i+f)*lda+i+c],fa-i-f,ca-i-c,lda);
        }
    }

}

int main()
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

  printf("De el tamaño de bloque mult mat: ");
  scanf("%d",&tbm);

  printf("De los valores inferior y superior: ");
  scanf("%lf %lf",&l,&u);

  generar_matriz_ld(a,fa,ca,lda,l,u);
  copiar_matriz(a,fa,ca,lda,copia,fa,ca,lda);

#ifdef DEBUG
  escribir_matriz_ld(a,fa,ca,lda);
#endif

// LU con bloques
    gettimeofday(tv,tz);
    si=(tv->tv_sec);
    ti1=(tv->tv_usec);
    ti1=si*1000000+ti1;

    lu_bloques(a,fa,ca,lda,tb);

    gettimeofday(tv,tz);
    sf=(tv->tv_sec);
    tf1=(tv->tv_usec);
    tf1=sf*1000000+tf1;
    printf("LU con bloques de %d:%d: Tamaño %d: %.6lf seg \n",tb,tbm,fa,(tf1-ti1)/1000000.);

    #ifdef DEBUG
    escribir_matriz_ld(a,fa,ca,lda);
    comprobar_lu(copia,fa,ca,lda,a,fa,ca,lda);
    #endif

 
  free(tv);
  free(tz);
  free(a);


}









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

void comprobar_lu(double *b,int fb,int cb,int ldb,double *a,int fa,int ca,int lda)
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





void copiar_matriz(double *mo,int fo,int co,int ldo,double *md,int fd,int cd,int ldd)
{
  int i,j;

  for(i=0;i<fo;i++)
    for(j=0;j<co;j++)
      md[i*ldd+j]=mo[i*ldo+j];
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

// md = md - mo
void restar_a_matriz(double *mo,int fo,int co,int ldo,double *md,int fd,int cd,int ldd)
{
  int i,j;

  for(i=0;i<fo;i++)
    for(j=0;j<co;j++)
      md[i*ldd+j]-=mo[i*ldo+j];
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

void matriz_matriz_resta_bloques(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc,int tb)
{
  int i,j,k, fblocka, cblocka, cblockb;
  double *s;

  s=(double *) malloc(sizeof(double)*tb*tb);

  for(i=0;i<fa;i=i+tb)
  {
    fblocka = ( (i+tb) < fa ) ? tb : (fa-i);
    for(j=0;j<cb;j=j+tb)
    {
      cblockb = ( (j+tb) < cb ) ? tb : (cb-j);
      matriz_cero(s,tb,tb,tb);
      for(k=0;k<ca;k=k+tb)
      {
        cblocka = ( (k+tb) < ca ) ? tb : (ca-k); // cblocka == fblockb
        multiplicar_acumular_matrices(&a[i*lda+k],fblocka,cblocka,lda,&b[k*ldb+j],cblocka,cblockb,ldb,s,fblocka,cblockb,tb);
      }
      restar_a_matriz(s,fblocka,cblockb,tb,&c[i*ldc+j],fblocka,cblockb,ldc);
    }
  }
  free(s);
}


void multiplicar_restar_matrices(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc)
{
  matriz_matriz_resta_bloques(a,fa,ca,lda,b,fb,cb,ldb,c,fc,cc,ldc,tbm);
}


// void multiplicar_restar_matrices(double *a,int fa,int ca,int lda,double *b,int fb,int cb,int ldb,double *c,int fc,int cc,int ldc)
// {
//   int i,j,k,kb;
//   double *da,*db,s;
//   for(i=0;i<fa;i++)
//   {
//     da=&a[i*lda];
//     for(j=0;j<cb;j++)
//     {
//       db=&b[j];
//       s=c[i*ldc+j];
//       for(k=0,kb=0;k<ca;k++,kb=kb+ldb)
//       {
//         s=s-da[k]*db[kb];
//       }
//       c[i*ldc+j]=s;
//     }
//   }
// }


void sistema_triangular_inferior(double *a,int fa,int ca,int lda,double *x,int fx,int cx,int ldx)
{
  int i,j,k;
  for(i=0;i<fa;i++)
  {      
    for(j=0;j<cx;j++)
      x[i*ldx+j]/=a[i*lda+i];
    for(j=i+1;j<fa;j++)
      for(k=0;k<cx;k++)
        x[j*ldx+k]-=x[i*ldx+k]*a[j*lda+i];
  }
}

void sistema_triangular_superior(double *a,int fa,int ca,int lda,double *x,int fx,int cx,int ldx)
{
  int i,j,k;
  for(i=0;i<ca;i++)
  {
    for(j=i+1;j<ca;j++)
        for(k=0;k<fx;k++)
            x[k*ldx+j]-=x[k*ldx+i]*a[i*lda+j];
  }
}
