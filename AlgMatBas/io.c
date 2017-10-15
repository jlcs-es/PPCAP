#include <stdio.h>
#include <stdlib.h>

void generar_vector(double *v,int n,double l,double u)
{
  int i;

  for(i=0;i<n;i++)
    v[i]=((double) rand()/RAND_MAX)*(u-l)+l;
}

void escribir_vector(double *v,int n)
{
  int i;

  for(i=0;i<n;i++)
    printf("%.6lf ",v[i]);
  printf("\n");
}


void generar_matriz(double **a,int n,int m,double l,double u)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      a[i][j]=((double) rand()/RAND_MAX)*(u-l)+l;
}

void escribir_matriz(double **a,int n,int m)
{
  int i,j;

  for(i=0;i<n;i++)
  {
    for(j=0;j<m;j++)
      printf("%.6lf ",a[i][j]);
    printf("\n");
  }
  printf("\n");
}

void generar_matriz_ld(double *a,int n,int m,int ld,double l,double u)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      a[i*ld+j]=((double) rand()/RAND_MAX)*(u-l)+l;
}

void escribir_matriz_ld(double *a,int n,int m,int ld)
{
  int i,j;

  for(i=0;i<n;i++)
  {
    for(j=0;j<m;j++)
      printf("%.6lf ",a[i*ld+j]);
    printf("\n");
  }
  printf("\n");
}


void generar_matriz_dispersa_fast(double *m,int filas,int columnas,int *fm,int *cm,int ndm,double l,double u)
{
  int i,rowNum,j,k,c;
  int total[filas], index[filas];
  for(i=0;i<filas;i++) 
    total[i]=0; 
  for(i=0;i<ndm;i++) 
  { 
    rowNum=(int) ((double) rand()/RAND_MAX*filas);;
    total[rowNum]++;
    m[i]=((double) rand()/RAND_MAX)*(u-l)+l; // creamos de paso los datos aleatorios
  }
  index[0]=0;
  for(i=1;i<filas;i++){
    index[i]=index[i-1]+total[i-1];
  }
  for(i=0;i<filas;i++){
    if(total[i]>0)
      cm[index[i]]=(int) ((double) rand()/RAND_MAX*columnas);
    for(j=1;j<total[i];j++){
      fm[index[i]+j] = i;
      k=j;
      c=(int) ((double) rand()/RAND_MAX*columnas);
      // TODO reordenar dentro de la fila actual
      // ir hacia atrás mientras fila actual && get ; comprobar == o <
    }
  }
  

  for(i=0;i<ndb;i++) 
  { 
    rowNum=cb[i];
    j=index[rowNum];
    index[rowNum]++;
    fr[j]=cb[i];
    cr[j]=fb[i];
    dr[j]=db[i];
  }


  int i,j,k,f,c;
 
  for(i=0;i<ndm;i++)
  {
    f=(int) ((double) rand()/RAND_MAX*filas);
    c=(int) ((double) rand()/RAND_MAX*columnas);
    j=0;
    while(j<i && (fm[j]<f || (fm[j]==f && cm[j]<c)))
      j++;
    if(j<i && fm[j]==f && cm[j]==c)
      i--;
    else
    {
      for(k=i;k>j;k--)
      {
        m[k]=m[k-1];
        fm[k]=fm[k-1];
        cm[k]=cm[k-1];
      }
      m[j]=((double) rand()/RAND_MAX)*(u-l)+l;
      fm[j]=f;
      cm[j]=c;
    }
  }
}

void generar_matriz_dispersa(double *m,int filas,int columnas,int *fm,int *cm,int ndm,double l,double u)
{
  int i,j,k,f,c;
 
  for(i=0;i<ndm;i++)
  {
    f=(int) ((double) rand()/RAND_MAX*filas);
    c=(int) ((double) rand()/RAND_MAX*columnas);
    j=0;
    while(j<i && (fm[j]<f || (fm[j]==f && cm[j]<c)))
      j++;
    if(j<i && fm[j]==f && cm[j]==c)
      i--;
    else
    {
      for(k=i;k>j;k--)
      {
        m[k]=m[k-1];
        fm[k]=fm[k-1];
        cm[k]=cm[k-1];
      }
      m[j]=((double) rand()/RAND_MAX)*(u-l)+l;
      fm[j]=f;
      cm[j]=c;
    }
  }
}

void generar_vector_disperso(double *v,int filas,int *fv,int ndv,double l,double u)
{
  int i,j,k,f;
 
  for(i=0;i<ndv;i++)
  {
    f=(int) ((double) rand()/RAND_MAX*filas);
    j=0;
    while(j<i && fv[j]<f)
      j++;
    if(j<i && fv[j]==f)
      i--;
    else
    {
      for(k=i;k>j;k--)
      {
        v[k]=v[k-1];
        fv[k]=fv[k-1];
      }
      v[j]=((double) rand()/RAND_MAX)*(u-l)+l;
      fv[j]=f;
    }
  }
}

void escribir_matriz_dispersa(double *m,int *fm,int *cm,int ndm)
{
  int i;

  for(i=0;i<ndm;i++)
    printf("%d, %d: %.6lf\n",fm[i],cm[i],m[i]);
  printf("\n");
}

void escribir_vector_disperso(double *v,int *fv,int ndv)
{
  int i;

  for(i=0;i<ndv;i++)
    printf("%d: %.6lf\n",fv[i],v[i]);
  printf("\n");
}

